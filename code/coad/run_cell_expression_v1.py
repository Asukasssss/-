"""Read complete raw-count inputs on server165, aggregate donor-equal summaries."""
import argparse
import gzip
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
from cell_origin_v1 import ROOT, KEYS, sha, vector, load_metadata, levels, base, write_table

TARGETS = ['HDC', 'GSTA4']

def read_lee(path, ids):
    targets, totals, n_genes, seen = {}, np.zeros(len(ids), dtype=np.int64), 0, set()
    with gzip.open(path, 'rt') as f:
        header = [s.strip('"') for s in f.readline().rstrip('\r\n').split('\t')]
        if len(header) == len(ids) + 1:
            header = header[1:]
        assert len(header) == len(ids) == len(set(header))
        assert set(header) == set(ids)
        order = pd.Index(header).get_indexer(ids)
        for line in f:
            name, values = line.rstrip('\r\n').split('\t', 1)
            name = name.strip('"')
            v = np.fromstring(values, dtype=np.float64, sep='\t')
            assert len(v) == len(ids) and np.isfinite(v).all()
            assert (v >= 0).all() and (v == np.floor(v)).all()
            assert name not in seen
            seen.add(name)
            v = v.astype(np.int64)[order]
            totals += v
            if name in TARGETS:
                targets[name] = v
            n_genes += 1
    assert set(targets) == set(TARGETS) and n_genes > 10000 and (totals > 0).all()
    return totals, targets, {'genes': n_genes, 'cells': len(ids), 'id_set_and_order_alignment': 'PASS', 'all_counts_nonnegative_integer': True, 'gzip_crc_complete': True}

def read_pelka(d, ids):
    import h5py
    names = vector(d / 'colon10x_default_dDvec_geneID.csv.gz')
    ens = vector(d / 'colon10x_default_dDvec_ensgID.csv.gz')
    assert len(names) == len(ens) and len(names) > 10000
    assert np.array_equal(ids, vector(d / 'colon10x_default_dDvec_sampleID.csv.gz'))
    pos = {}
    for gene in TARGETS:
        idx = np.flatnonzero(names == gene)
        assert len(idx) == 1, (gene, len(idx))
        pos[gene] = int(idx[0]) + 1
    totals = np.zeros(len(ids), dtype=np.int64)
    targets = {g: np.zeros(len(ids), dtype=np.int64) for g in TARGETS}
    with h5py.File(d / 'colon10x_default_dSp_rawCount.h5') as f:
        assert f['i'].shape == f['j'].shape == f['v'].shape
        assert f['i'].shape[0] == 1
        nnz = f['i'].shape[1]
        for start in range(0, nnz, 2000000):
            end = min(start + 2000000, nnz)
            i, j, v = [f[k][0, start:end] for k in ['i', 'j', 'v']]
            for a in [i, j, v]:
                assert np.isfinite(a).all() and np.equal(a, np.floor(a)).all()
            assert i.min() >= 1 and i.max() <= len(names)
            assert j.min() >= 1 and j.max() <= len(ids) and (v >= 0).all()
            j = j.astype(np.int64) - 1
            totals += np.bincount(j, weights=v, minlength=len(ids)).astype(np.int64)
            for gene, index in pos.items():
                mask = i == index
                targets[gene] += np.bincount(j[mask], weights=v[mask], minlength=len(ids)).astype(np.int64)
            if start % 100000000 == 0:
                print(json.dumps({'triplets_read': end, 'triplets_total': nnz}), flush=True)
    assert (totals > 0).all()
    return totals, targets, {'genes': len(names), 'cells': len(ids), 'sparse_triplets': nnz, 'source_index_base': 1, 'all_counts_nonnegative_integer': True, 'all_gene_denominator': True}

def summarize(meta, total, targets, out, min_cells, min_patients):
    rows, private = [], []
    for frame in levels(meta):
        frame['library_UMI'] = total
        frame['cells'] = 1
        for gene, values in targets.items():
            frame[gene + '_UMI'] = values
            frame[gene + '_detected'] = (values > 0).astype(int)
        value_cols = ['library_UMI', 'cells'] + [g + s for g in TARGETS for s in ['_UMI', '_detected']]
        p = frame.groupby(KEYS + ['patient'], dropna=False)[value_cols].sum().reset_index()
        assert int(p.cells.sum()) == len(meta)
        assert int(p.library_UMI.sum()) == int(total.sum())
        private.append(p)
        for key, group in p.groupby(KEYS, sort=True, dropna=False):
            for subset in ['all_covered_patients', 'patients_with_at_least_20_cells']:
                q = group if subset == 'all_covered_patients' else group[group.cells >= min_cells]
                for gene in TARGETS:
                    n = len(q)
                    status = 'DONE' if n >= min_patients else 'PARTIAL'
                    reason = 'Descriptive expression; no association or functional validation' if n >= min_patients else 'LIMITED_COVERAGE; descriptive only'
                    if n < 3:
                        status, reason = 'NOT_EVALUABLE', 'Expression summary suppressed: fewer than 3 patients; not gene-negative'
                    row = base(out.name, key, 'donor_equal_expression', n, status, reason)
                    row.update(gene=gene, n_reference=len(group), effect_type='median_patient_pseudobulk_CPM',
                        support_set=subset, n_cells=int(q.cells.sum()), n_patients_ge20=int((group.cells >= min_cells).sum()))
                    for metric in ['detection_fraction', 'pseudobulk_CPM', 'target_UMI', 'library_UMI']:
                        for stat in ['min', 'q25', 'median', 'q75', 'max']:
                            row[metric + '_' + stat] = 'NA'
                    if n >= 3:
                        values = {'detection_fraction': q[gene + '_detected'] / q.cells,
                            'pseudobulk_CPM': 1e6 * q[gene + '_UMI'] / q.library_UMI,
                            'target_UMI': q[gene + '_UMI'], 'library_UMI': q.library_UMI}
                        assert (values['detection_fraction'].between(0, 1)).all()
                        for metric, v in values.items():
                            for stat, quant in [('min', 0), ('q25', .25), ('median', .5), ('q75', .75), ('max', 1)]:
                                row[metric + '_' + stat] = float(v.quantile(quant))
                        row['effect'] = row['pseudobulk_CPM_median']
                    rows.append(row)
    # Private patient-level quantities never go to GitHub or local handoff.
    pd.concat(private, ignore_index=True).to_csv(out / ('private_' + meta.study.iloc[0] + '_patient_expression.tsv'), sep='\t', index=False)
    return rows

def main():
    import h5py
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--spec', type=Path, required=True)
    p.add_argument('--lock-commit', required=True)
    p.add_argument('--study', choices=['Lee', 'Pelka'], required=True)
    a = p.parse_args()
    assert a.data_dir.resolve().parent == ROOT / 'data/candidates'
    assert a.out.resolve().parent == ROOT / 'results/collaborative/COAD/B'
    assert (a.out / '.running').exists() and len(a.lock_commit) == 40
    spec = json.loads(a.spec.read_text())
    assert spec['targets'] == TARGETS and spec['minimum_cells_per_patient_type'] == 20 and spec['minimum_qualified_patients'] == 5
    file = a.out / (a.study + '_source_expression_summary.tsv')
    assert not file.exists()
    metadata = load_metadata(a.data_dir)
    study = next(k for k in metadata if k.startswith(a.study))
    meta = metadata[study]
    matrix = a.data_dir / ('lee_raw_UMI.txt.gz' if a.study == 'Lee' else 'colon10x_default_dSp_rawCount.h5')
    total, targets, checks = read_lee(matrix, meta.cell_id.to_numpy()) if a.study == 'Lee' else read_pelka(a.data_dir, meta.cell_id.to_numpy())
    assert all((v <= total).all() for v in targets.values())
    np.savez_compressed(a.out / ('private_' + a.study + '_cell_UMI.npz'), total=total, **targets)
    rows = summarize(meta, total, targets, a.out, 20, 5)
    write_table(rows, file)
    checks.update(study=study, status='PASS', rows=len(rows), matrix_sha256=sha(matrix), spec_sha256=sha(a.spec),
        pre_expression_remote_lock=a.lock_commit, code_sha256=sha(__file__), python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, h5py=h5py.__version__,
        target_total_UMI={g:int(v.sum()) for g,v in targets.items()}, all_gene_total_UMI=int(total.sum()),
        candidate_or_covariate_recalculation=False, patient_rows_exported=False)
    (a.out / (a.study + '_validation.json')).write_text(json.dumps(checks, indent=2))
    print(json.dumps(checks), flush=True)

if __name__ == '__main__':
    main()
