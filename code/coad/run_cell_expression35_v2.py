"""Expand fixed cell-source summaries to all35 genes; raw/private data stay server165."""
import argparse
import csv
import gzip
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
from cell_origin_v1 import ROOT, KEYS, PREFIX, sha, vector, load_metadata, levels, base, write_table

VERSION = 'COAD_cell_origin35_v2'
OLD_RUN = '20260920T072000Z_cell_expression_v1'

def symbol_choice(gene, names, aliases):
    hits = [i for i, name in enumerate(names) if name == gene]
    if len(hits) == 1:
        return hits[0], gene, 'EXACT_SOURCE_SYMBOL'
    if len(hits) > 1:
        return None, 'NA', 'AMBIGUOUS_SOURCE_SYMBOL'
    options = [(i, name) for i, name in enumerate(names) if name in aliases.get(gene, [])]
    if len(options) == 1:
        return options[0][0], options[0][1], 'HGNC_CONFIRMED_PREVIOUS_SYMBOL'
    return None, 'NA', 'SOURCE_GENE_NOT_FOUND' if not options else 'AMBIGUOUS_ALIAS'

def extract_lee(path, ids, genes, aliases):
    # Only target rows need parsing; full-gene denominator was validated in v1.
    candidates = set(genes) | {v for g in genes for v in aliases.get(g, [])}
    names, selected = [], {}
    with gzip.open(path, 'rt') as f:
        h = [v.strip('"') for v in f.readline().rstrip('\r\n').split('\t')]
        if len(h) == len(ids) + 1:
            h = h[1:]
        assert len(h) == len(ids) == len(set(h)) and set(h) == set(ids)
        order = pd.Index(h).get_indexer(ids)
        for line in f:
            name, text = line.rstrip('\r\n').split('\t', 1)
            name = name.strip('"')
            names.append(name)
            if name in candidates:
                v = np.fromstring(text, dtype=float, sep='\t')
                assert len(v) == len(ids) and np.isfinite(v).all() and (v >= 0).all() and (v == np.floor(v)).all()
                selected[len(names) - 1] = v.astype(np.int64)[order]
    counts, mapping = {}, []
    for gene in genes:
        index, source, status = symbol_choice(gene, names, aliases)
        mapping.append(dict(gene=gene, source_symbol=source, mapping_status=status, source_feature_index_0based=index))
        if index is not None:
            counts[gene] = selected[index]
    return counts, mapping, len(names)

def extract_pelka(d, ids, genes, aliases):
    import h5py
    names = vector(d / 'colon10x_default_dDvec_geneID.csv.gz').tolist()
    ens = vector(d / 'colon10x_default_dDvec_ensgID.csv.gz')
    assert len(names) == len(ens)
    assert np.array_equal(ids, vector(d / 'colon10x_default_dDvec_sampleID.csv.gz'))
    mapping, found, lookup = [], [], np.full(len(names) + 1, -1, dtype=np.int32)
    for gene in genes:
        index, source, status = symbol_choice(gene, names, aliases)
        mapping.append(dict(gene=gene, source_symbol=source, mapping_status=status,
            source_feature_index_0based=index, source_ensembl=ens[index] if index is not None else 'NA'))
        if index is not None:
            lookup[index + 1] = len(found)
            found.append(gene)
    counts = np.zeros((len(found), len(ids)), dtype=np.int64)
    with h5py.File(d / 'colon10x_default_dSp_rawCount.h5') as f:
        nnz = f['i'].shape[1]
        for start in range(0, nnz, 2000000):
            stop = min(start + 2000000, nnz)
            i = f['i'][0, start:stop].astype(np.int64)
            assert i.min() >= 1 and i.max() <= len(names)
            target = lookup[i]
            mask = target >= 0
            if mask.any():
                j = f['j'][0, start:stop][mask].astype(np.int64) - 1
                v = f['v'][0, start:stop][mask]
                assert j.min() >= 0 and j.max() < len(ids)
                assert np.isfinite(v).all() and (v >= 0).all() and (v == np.floor(v)).all()
                np.add.at(counts, (target[mask], j), v.astype(np.int64))
            if start % 100000000 == 0:
                print(json.dumps({'triplets_read': stop, 'triplets_total': nnz}), flush=True)
    return dict(zip(found, counts)), mapping, len(names)

def summarize(meta, total, counts, genes, out):
    found = [g for g in genes if g in counts]
    rows = {g: [] for g in genes}
    private = []
    for frame in levels(meta):
        values = {'library_UMI': total, 'cells': np.ones(len(meta), dtype=np.int64)}
        for gene in found:
            values[gene + '_UMI'] = counts[gene]
            values[gene + '_detected'] = (counts[gene] > 0).astype(np.int64)
        f = pd.concat([frame.reset_index(drop=True), pd.DataFrame(values)], axis=1)
        p = f.groupby(KEYS + ['patient'], dropna=False)[list(values)].sum().reset_index()
        assert int(p.cells.sum()) == len(meta) and int(p.library_UMI.sum()) == int(total.sum())
        private.append(p)
        for key, group in p.groupby(KEYS, sort=True, dropna=False):
            for subset in ['all_covered_patients', 'patients_with_at_least_20_cells']:
                q = group if subset == 'all_covered_patients' else group[group.cells >= 20]
                n = len(q)
                quantiles = {}
                if n >= 3 and found:
                    raw = q[[g + '_UMI' for g in found]].to_numpy(dtype=float)
                    detected = q[[g + '_detected' for g in found]].to_numpy(dtype=float)
                    arrays = {'detection_fraction': detected / q.cells.to_numpy()[:, None],
                        'pseudobulk_CPM': 1e6 * raw / q.library_UMI.to_numpy()[:, None],
                        'target_UMI': raw, 'library_UMI': np.tile(q.library_UMI.to_numpy()[:, None], (1, len(found)))}
                    quantiles = {m: np.quantile(a, [0, .25, .5, .75, 1], axis=0) for m, a in arrays.items()}
                for gene in genes:
                    status = 'DONE' if n >= 5 else 'PARTIAL'
                    reason = 'Descriptive expression; no association or functional validation' if n >= 5 else 'LIMITED_COVERAGE; descriptive only'
                    if n < 3:
                        status, reason = 'NOT_EVALUABLE', 'Expression summary suppressed: fewer than 3 patients; not gene-negative'
                    if gene not in counts:
                        status, reason = 'NOT_EVALUABLE', 'SOURCE_GENE_NOT_FOUND_OR_AMBIGUOUS; not measured zero; no paralog substitution'
                    row = base(out.name, key, 'donor_equal_expression', n, status, reason)
                    row.update(analysis_version=VERSION, gene=gene, n_reference=len(group),
                        effect_type='median_patient_pseudobulk_CPM', support_set=subset, n_cells=int(q.cells.sum()),
                        n_patients_ge20=int((group.cells >= 20).sum()))
                    for metric in ['detection_fraction', 'pseudobulk_CPM', 'target_UMI', 'library_UMI']:
                        for j, stat in enumerate(['min', 'q25', 'median', 'q75', 'max']):
                            row[metric + '_' + stat] = float(quantiles[metric][j, found.index(gene)]) if n >= 3 and gene in counts else 'NA'
                    row['effect'] = row['pseudobulk_CPM_median']
                    row.update(computation='NEW_GENE_SUMMARY', reuse_run_id='NA')
                    rows[gene].append(row)
    pd.concat(private, ignore_index=True).to_csv(out / ('private_' + meta.study.iloc[0] + '_patient_expression.tsv'), sep='\t', index=False)
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--spec', type=Path, required=True)
    p.add_argument('--lock-commit', required=True)
    p.add_argument('--study', choices=['Lee', 'Pelka'], required=True)
    a = p.parse_args()
    assert a.out.resolve().parent == ROOT / 'results/collaborative/COAD/B' and (a.out / '.running').exists()
    assert a.data_dir.resolve().parent == ROOT / 'data/candidates'
    spec = json.loads(a.spec.read_text()); genes = spec['targets']; new = spec['new_targets']
    assert len(genes) == len(set(genes)) == 35 and len(new) == 33
    assert spec['minimum_cells_per_patient_type'] == 20 and spec['minimum_qualified_patients'] == 5
    old = a.out.parent / OLD_RUN
    oldcheck = json.loads((old / (a.study + '_validation.json')).read_text())
    matrix = a.data_dir / ('lee_raw_UMI.txt.gz' if a.study == 'Lee' else 'colon10x_default_dSp_rawCount.h5')
    assert sha(matrix) == oldcheck['matrix_sha256']
    metadata_hashes = json.loads((a.out.parent / '20260920T050000Z_cell_origin_v1/source_hashes_metadata.json').read_text())
    for record in metadata_hashes:
        if record['name'].startswith(('lee_annotation', 'pelka_author', 'colon10x_default_dDvec')):
            assert sha(a.data_dir / record['name']) == record['sha256']
    data = load_metadata(a.data_dir); meta = data[next(k for k in data if k.startswith(a.study))]
    cache = old / ('private_' + a.study + '_cell_UMI.npz')
    with np.load(cache) as archive:
        total = archive['total']
    assert len(total) == len(meta) and (total > 0).all() and int(total.sum()) == oldcheck['all_gene_total_UMI']
    counts, mapping, features = extract_lee(matrix, meta.cell_id.to_numpy(), new, spec['confirmed_symbol_aliases']) if a.study == 'Lee' else extract_pelka(a.data_dir, meta.cell_id.to_numpy(), new, spec['confirmed_symbol_aliases'])
    assert all((v <= total).all() for v in counts.values())
    np.savez_compressed(a.out / ('private_' + a.study + '_new33_cell_UMI.npz'), **counts)
    public = a.out / ('public_' + a.study); public.mkdir()
    (public / 'by_gene').mkdir()
    rows = summarize(meta, total, counts, new, a.out)
    for gene, records in rows.items():
        write_table(records, public / 'by_gene' / (gene + '.tsv'))
    with (old / (a.study + '_source_expression_summary.tsv')).open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t'); header = reader.fieldnames
        reused = {g: [] for g in spec['reuse_targets']}
        for row in reader:
            gene = row['gene']; row.update(run_id=a.out.name, analysis_version=VERSION,
                computation='EXACT_V1_SUMMARY_REUSE', reuse_run_id=OLD_RUN); reused[gene].append(row)
    for gene, records in reused.items():
        with (public / 'by_gene' / (gene + '.tsv')).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=header + ['computation', 'reuse_run_id'], delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(records)
        mapping.append(dict(gene=gene, source_symbol=gene, mapping_status='EXACT_V1_REUSE', source_feature_index_0based='NA'))
    pd.DataFrame(mapping).sort_values('gene').to_csv(public / 'gene_availability.tsv', sep='\t', index=False, na_rep='NA')
    checks = dict(status='PASS', analysis_version=VERSION, study=a.study, cells=len(meta), genes_requested=35,
        new_genes_available=len(counts), missing_new_genes=sorted(set(new) - set(counts)), reused_genes=spec['reuse_targets'],
        rows=sum(map(len, rows.values())) + sum(map(len, reused.values())), source_features=features,
        matrix_sha256=oldcheck['matrix_sha256'], previous_denominator_npz_sha256=sha(cache),
        denominator_reused_from=OLD_RUN, metadata_hashes_unchanged=True, lock_commit=a.lock_commit,
        spec_sha256=sha(a.spec), code_sha256=sha(__file__), python=platform.python_version(),
        numpy=np.__version__, pandas=pd.__version__, new_hypothesis_tests=False, patient_rows_exported=False)
    (public / 'validation.json').write_text(json.dumps(checks, indent=2))
    print(json.dumps(checks), flush=True)

if __name__ == '__main__':
    main()
