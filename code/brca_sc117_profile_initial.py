"""Server-only, streamed raw-count descriptive profiles; no hypothesis tests.

Run after inspect/configuration. Cell/donor-level outputs remain in run/private.
Public output contains only aggregates with >=3 source donor labels.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse


def column(g, key):
    x = g[key]
    if isinstance(x, h5py.Dataset):
        return x.asstr()[:] if x.dtype.kind in 'OS' else x[:]
    if 'categories' in x:
        cat = column(x, 'categories')
        codes = x['codes'][:]
        result = np.full(len(codes), 'NA', dtype=object)
        result[codes >= 0] = cat[codes[codes >= 0]]
        return result
    if 'values' in x:
        vals = column(x, 'values').astype(object)
        vals[x['mask'][:]] = 'NA'
        return vals
    raise ValueError('Unsupported column: ' + key)


def dataframe(g):
    key = g.attrs.get('_index', '_index')
    if isinstance(key, bytes):
        key = key.decode()
    return pd.DataFrame({k: column(g, k) for k in g if k != key},
                        index=column(g, key))


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', lineterminator='\n')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run(root, config):
    cfg = json.loads(config.read_text())
    cohort = cfg['cohort']
    pub, private = root / 'public', root / 'private'
    genes = pd.read_csv(root / 'source/candidate_comparison_v1.tsv', sep='\t')['gene'].tolist()
    assert len(genes) == len(set(genes)) == 117
    source = root / 'source' / cfg['file']
    with h5py.File(source, 'r') as h:
        obs = dataframe(h['obs'])
        assert obs.index.is_unique, 'Nonunique cell identifiers'
        assert cfg['raw_path'] in h
        mat = h[cfg['raw_path']]
        assert mat.attrs['encoding-type'] == 'csr_matrix', 'Inspect non-CSR before use'
        shape = tuple(mat.attrs['shape'])
        assert shape[0] == len(obs)
        var = dataframe(h[cfg['var_path']])
        names = var[cfg['symbol_column']].astype(str)
        assert var.index.is_unique
        idx = []
        coverage = []
        for gene in genes:
            positions = np.flatnonzero(names.values == gene)
            ok = len(positions) == 1
            idx.append(int(positions[0]) if ok else -1)
            coverage.append({'cohort': cohort, 'gene': gene, 'n_feature_matches': len(positions),
                             'status': 'DONE' if ok else 'NOT_EVALUABLE',
                             'reason': '' if ok else 'absent_or_ambiguous_feature'})
        write_tsv(pd.DataFrame(coverage), pub / (cohort + '_feature_coverage.tsv'))
        valid_genes = np.array(idx) >= 0
        selected = np.ones(len(obs), dtype=bool)
        if cfg.get('filter_column'):
            selected &= obs[cfg['filter_column']].astype(str).eq(cfg['filter_value']).values
        donor = obs[cfg['donor_column']].astype(str)
        types = obs[cfg['celltype_column']].astype(str)
        mapping = cfg['celltype_map']
        assert set(types[selected]).issubset(mapping), sorted(set(types[selected]) - set(mapping))
        broad = types.map(mapping).fillna('EXCLUDED')
        strata = obs[cfg['stratum_column']].astype(str) if cfg.get('stratum_column') else pd.Series('all', index=obs.index)
        treatment = obs[cfg['treatment_column']].astype(str) if cfg.get('treatment_column') else pd.Series('not_available', index=obs.index)
        invalid = donor.isin(['NA', 'nan', 'unknown', '']).values
        assert not (selected & invalid).any(), 'Missing donor labels require explicit handling'
        metadata = pd.DataFrame({'donor': donor, 'celltype': broad, 'stratum': strata, 'treatment': treatment})
        # Pool repeated cells/fractions within the same source donor x lineage x stratum.
        keys = pd.MultiIndex.from_frame(metadata[['donor', 'celltype', 'stratum', 'treatment']])
        codes, uniques = pd.factorize(keys)
        codes[~selected] = -1
        group = uniques.to_frame(index=False)
        group.columns = ['donor', 'celltype', 'stratum', 'treatment']
        ng = len(group)
        n = np.zeros(ng, dtype=np.int64)
        logsum = np.zeros((ng, 117), dtype=np.float64)
        detected = np.zeros((ng, 117), dtype=np.int64)
        rawsum = np.zeros((ng, 117), dtype=np.float64)
        librarysum = np.zeros(ng, dtype=np.float64)
        ip = mat['indptr'][:]
        checked = 0
        zero_libraries = 0
        for start in range(0, len(obs), cfg['chunk_cells']):
            end = min(start + cfg['chunk_cells'], len(obs))
            keep = selected[start:end]
            if not keep.any():
                continue
            lo, hi = ip[start], ip[end]
            values = mat['data'][lo:hi]
            assert np.isfinite(values).all() and (values >= 0).all()
            assert np.allclose(values, np.rint(values), atol=1e-6), 'Raw layer is not integer counts'
            checked += len(values)
            x = sparse.csr_matrix((values, mat['indices'][lo:hi], ip[start:end+1] - lo),
                                  shape=(end-start, shape[1]))[keep]
            lib = np.asarray(x.sum(axis=1)).ravel().astype(float)
            zero_libraries += int((lib <= 0).sum())
            assert (lib > 0).all()
            y = np.zeros((len(lib), 117), dtype=np.float64)
            y[:, valid_genes] = x[:, np.array(idx)[valid_genes]].toarray()
            z = np.log1p(y * (10000.0 / lib)[:, None])
            c = codes[start:end][keep]
            np.add.at(n, c, 1)
            np.add.at(librarysum, c, lib)
            np.add.at(logsum, c, z)
            np.add.at(detected, c, (y > 0).astype(np.int64))
            np.add.at(rawsum, c, y)
            if start % 50000 == 0:
                print(cohort, start, '/', len(obs), flush=True)
        assert n.sum() == selected.sum()
        # Metadata and individual measurements never exported to public/.
        private_rows = []
        for j, g in group.iterrows():
            if n[j] == 0:
                continue
            for k, gene in enumerate(genes):
                private_rows.append(dict(g, gene=gene, n_cells=int(n[j]),
                    mean_log1p10k=logsum[j,k]/n[j] if valid_genes[k] else np.nan,
                    detection_fraction=detected[j,k]/n[j] if valid_genes[k] else np.nan,
                    raw_count=rawsum[j,k] if valid_genes[k] else np.nan,
                    all_gene_library_sum=librarysum[j]))
        d = pd.DataFrame(private_rows)
        write_tsv(d, private / (cohort + '_donor_profiles.tsv'))
        write_tsv(metadata[selected].reset_index(names='barcode'), private / (cohort + '_cell_metadata.tsv'))
        # Subtype and treatment partitions plus an ALL partition. Pool same donor if strata differ.
        summary_rows = []
        partitions = [('ALL', d)]
        partitions += [('subtype:' + str(k), v) for k, v in d.groupby('stratum')]
        if cfg.get('treatment_column'):
            partitions += [('treatment:' + str(k), v) for k, v in d.groupby('treatment')]
        for partition, block in partitions:
            pooled = []
            for (don, celltype, gene), v in block.groupby(['donor', 'celltype', 'gene'], sort=True):
                nn = v.n_cells.sum()
                pooled.append({'donor': don, 'celltype': celltype, 'gene': gene, 'n_cells': nn,
                    'mean': np.average(v.mean_log1p10k, weights=v.n_cells),
                    'detect': np.average(v.detection_fraction, weights=v.n_cells)})
            pooled = pd.DataFrame(pooled)
            for (celltype, gene), v in pooled.groupby(['celltype', 'gene'], sort=True):
                usable = v[(v.n_cells >= cfg['min_cells_per_donor_group']) & v['mean'].notna()]
                ok = len(usable) >= cfg['min_donors_per_public_group']
                why = '' if ok else 'fewer_than_3_eligible_source_donor_labels_or_missing_gene'
                summary_rows.append({'cancer':'BRCA', 'cohort':cohort, 'stage_id':'06_EXTERNAL',
                    'run_id':root.name, 'analysis_version':'sc117_v1', 'analysis_type':'descriptive_cell_source',
                    'metabolite_key':'NA','metabolite_name':'NA','gene':gene,'unit':cfg['unit'],
                    'n':len(usable), 'n_reference':'NA','effect_type':'equal_donor_mean_log1p_counts_per_10k',
                    'effect':usable['mean'].mean() if ok else np.nan,
                    'ci_lower':np.nan,'ci_upper':np.nan,'p_value':np.nan,'q_value':np.nan,
                    'test_family':'NA_descriptive_no_tests','family_n_evaluable':'NA',
                    'status':'DONE' if ok else 'NOT_EVALUABLE','reason':why,'source_id':cfg['source_id'],
                    'partition':partition,'celltype':celltype,'n_source_labels_total':len(v),
                    'n_cells_total':int(v.n_cells.sum()),'n_cells_eligible':int(usable.n_cells.sum()),
                    'median_donor_mean':usable['mean'].median() if ok else np.nan,
                    'donor_mean_q25':usable['mean'].quantile(.25) if ok else np.nan,
                    'donor_mean_q75':usable['mean'].quantile(.75) if ok else np.nan,
                    'mean_detection_fraction':usable.detect.mean() if ok else np.nan,
                    'median_detection_fraction':usable.detect.median() if ok else np.nan,
                    'annotation_origin':cfg['annotation_origin']})
        result = pd.DataFrame(summary_rows)
        write_tsv(result, pub / (cohort + '_celltype_profiles.tsv'))
        audit = {'cohort':cohort,'status':'DONE','source_shape':list(map(int,shape)),
            'cells_selected':int(selected.sum()),'source_donor_labels':int(donor[selected].nunique()),
            'gene_count':117,'genes_unique_in_dictionary':int(valid_genes.sum()),
            'cells_by_celltype':broad[selected].value_counts().to_dict(),
            'source_labels_by_stratum':metadata[selected].groupby('stratum').donor.nunique().to_dict(),
            'source_labels_by_treatment':metadata[selected].groupby('treatment').donor.nunique().to_dict(),
            'unique_barcodes':True,'nonnegative_integer_raw_entries_checked':int(checked),
            'zero_libraries':zero_libraries,'aggregate_cell_count_verified':True,
            'public_rows':len(result),'all117_retained':set(result.gene)==set(genes),
            'new_hypothesis_tests':0,'normalization':'log1p(raw gene count/all retained raw gene counts*10000)',
            'statistical_independence_limit':cfg['independence_limit'],
            'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
            'scipy':scipy.__version__,'h5py':h5py.__version__,
            'script_sha256':sha256(__file__), 'config_sha256':sha256(config)}
        (pub / (cohort + '_validation.json')).write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
        print(json.dumps(audit,ensure_ascii=False),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--config', type=Path, required=True)
    a = p.parse_args()
    run(a.root, a.config)
