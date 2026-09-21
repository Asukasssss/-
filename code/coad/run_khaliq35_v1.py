"""Server-only author compartment counts; exact cell joins, donor-equal summaries."""
import argparse
import csv
import gzip
import io
import json
import platform
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl
from cell_origin_v1 import ROOT, sha, coverage, write_table
from run_cell_expression35_v2 import summarize, symbol_choice
import run_cell_expression35_v2 as summary_module

VERSION = 'COAD_Khaliq35_v1'

def inner_csv(outer, name):
    inner = zipfile.ZipFile(io.BytesIO(outer.read(name)))
    names = [n for n in inner.namelist() if n.endswith('.csv') and not n.startswith('__MACOSX/') and not Path(n).name.startswith('._')]
    assert len(names) == 1
    return inner, io.TextIOWrapper(inner.open(names[0]))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--spec', type=Path, required=True)
    p.add_argument('--lock-commit', required=True)
    a = p.parse_args()
    assert a.out.resolve().parent == ROOT / 'results/collaborative/COAD/B'
    assert (a.out / '.running').is_dir()
    spec = json.loads(a.spec.read_text()); genes = spec['targets']
    assert len(genes) == len(set(genes)) == 35
    assert spec['minimum_cells_per_patient_type'] == 20 and spec['minimum_qualified_patients'] == 5
    out = a.out / 'public'; out.mkdir(exist_ok=True); (out / 'by_gene').mkdir(exist_ok=True)
    geo = pd.read_csv(a.data_dir / 'annotation.csv.gz', index_col=0, keep_default_na=False)
    assert geo.index.is_unique and len(geo) == 49859
    table1 = list(openpyxl.load_workbook(a.data_dir / 'MOESM2.xlsx', read_only=True, data_only=True).active.values)
    table2 = list(openpyxl.load_workbook(a.data_dir / 'MOESM3.xlsx', read_only=True, data_only=True).active.values)
    assert table2[2][0] == 'Patient'
    samples = set(geo.samples)
    t1 = {r[0]: r for r in table1 if r[0] in samples}
    t2 = {r[0]: r for r in table2 if r[0] in samples}
    assert set(t1) == set(t2) == samples
    # Table S2 explicitly calls these labels Patient; analyze each tissue separately.
    # No cross-tissue pairing or suffix-derived linking is performed.
    sample_tissue = geo.groupby('samples').Condition.agg(lambda s: sorted(set(s)))
    assert all(len(v) == 1 for v in sample_tissue)
    study = 'Khaliq2022_GSE200997'
    z = zipfile.ZipFile(a.data_dir / 'celltypes.joined.zip')
    lookup = {
        'Bcells': ('Bcells', 1), 'Tcells': ('Tcells_including_NK_ILC', 2),
        'Myeloid': ('Myeloid', 3), 'Fibroblast': ('Fibroblast', 5),
        'Endothelial': ('Endothelial', 6), 'Epithelial': ('Epithelial', 8),
    }
    frames, totals, target_arrays, checks, maps = [], [], [], [], []
    for folder, (lineage, table_col) in lookup.items():
        names = [n for n in z.namelist() if n.startswith(folder + '/')]
        mn = [n for n in names if n.endswith('Metadata.csv.zip')]
        cn = [n for n in names if n.endswith('counts.csv.zip')]
        assert len(mn) == len(cn) == 1
        q, stream = inner_csv(z, mn[0]); m = pd.read_csv(stream, index_col=0, keep_default_na=False); stream.close(); q.close()
        assert m.index.is_unique and m.index.isin(geo.index).all()
        assert (m['orig.ident'] == geo.loc[m.index, 'samples']).all()
        for field in ['Condition', 'Location', 'MSI_Status']:
            assert (m[field] == geo.loc[m.index, field]).all()
        actual = m.groupby('orig.ident').size()
        expected = {s: int(t2[s][table_col] or 0) for s in samples}
        table_matches = all(int(actual.get(s, 0)) == v for s, v in expected.items())
        # A discrepancy must be visible; do not silently replace author file values.
        qc = dict(lineage=lineage, cells=len(m), samples=m['orig.ident'].nunique(), table_S2_exact_count_match=table_matches)
        sub = m.seurat_clusters.astype(str)
        # Preserve author codes; do not invent biological subtypes for numeric clusters.
        meta = pd.DataFrame(dict(cell_id=m.index, patient=m['orig.ident'].to_numpy(), specimen=m['orig.ident'].to_numpy(),
            library='NA', study=study, tissue=m.Condition.to_numpy(), enrichment='DAPI_VIABLE_NO_CD45_REPORTED',
            technology='10x_5prime', lineage=lineage, subtype=[lineage + '::' + x for x in sub]))
        q, stream = inner_csv(z, cn[0])
        h = next(csv.reader([stream.readline()]))[1:]
        assert len(h) == len(set(h)) == len(m) and set(h) == set(m.index)
        order = pd.Index(h).get_indexer(m.index)
        total = np.zeros(len(h), dtype=np.int64)
        found = {}; gene_names = []
        choices = set(genes) | {x for g in genes for x in spec['confirmed_symbol_aliases'].get(g, [])}
        for line in stream:
            name, text = line.rstrip('\r\n').split(',', 1); name = name.strip('"'); gene_names.append(name)
            v = np.fromstring(text, sep=',', dtype=float)
            assert len(v) == len(h) and np.isfinite(v).all() and (v >= 0).all() and (v == np.floor(v)).all()
            total += v.astype(np.int64)
            if name in choices:
                found.setdefault(name, []).append(v.astype(np.int64)[order])
        stream.close(); q.close(); total = total[order]
        assert (total > 0).all()
        match = total == pd.to_numeric(m.nCount_RNA).to_numpy()
        qc.update(features=len(gene_names), full_gene_UMI=int(total.sum()), nCount_RNA_exact_matches=int(match.sum()), nCount_RNA_exact_all=bool(match.all()))
        # Complete source UMI sum must agree with author's RNA library count.
        assert match.all(), 'Full-gene denominator disagrees with author nCount_RNA: ' + lineage
        count = {}
        for gene in genes:
            ix, symbol, status = symbol_choice(gene, gene_names, spec['confirmed_symbol_aliases'])
            maps.append(dict(gene=gene,lineage=lineage,source_symbol=symbol,mapping_status=status,source_feature_index_0based=ix))
            if ix is not None: count[gene] = found[symbol][0]
        assert all((v <= total).all() for v in count.values())
        frames.append(meta); totals.append(total); target_arrays.append(count); checks.append(qc)
        print(json.dumps(qc), flush=True)
    meta = pd.concat(frames, ignore_index=True); total = np.concatenate(totals)
    assert meta.cell_id.is_unique
    # Target availability must be identical across compartments before concatenation.
    assert all(set(c) == set(target_arrays[0]) for c in target_arrays)
    counts = {g: np.concatenate([c[g] for c in target_arrays]) for g in target_arrays[0]}
    np.savez_compressed(a.out / 'private_cell_counts.npz', total=total, **counts)
    meta.to_csv(a.out / 'private_cell_metadata.tsv', sep='\t', index=False)
    coverage({study:meta}, out)
    summary_module.VERSION = VERSION
    # summarize writes patient-level values to this private directory.
    rows = summarize(meta, total, counts, genes, a.out)
    for g in genes:
        write_table(rows[g], out / 'by_gene' / (g + '.tsv'))
    pd.DataFrame(maps).to_csv(out / 'gene_availability.tsv', sep='\t', index=False)
    pd.DataFrame(checks).to_csv(out / 'compartment_checks.tsv', sep='\t', index=False)
    # Coverage includes the unassigned remainder, not incorrectly imputed expression.
    remaining = geo.loc[~geo.index.isin(meta.cell_id)]
    missing = [dict(tissue=t,cells=len(f),samples=f.samples.nunique(),status='NOT_EVALUABLE',reason='Not present in author compartment exports;cell type not assigned') for t,f in remaining.groupby('Condition')]
    pd.DataFrame(missing).to_csv(out / 'unassigned_coverage.tsv',sep='\t',index=False)
    audit = dict(status='PASS',analysis_version=VERSION,lock_commit=a.lock_commit,cells_geo=len(geo),cells_analyzed=len(meta),cells_unassigned=len(remaining),genes_requested=len(genes),genes_available=len(counts),samples=meta.specimen.nunique(),tumor_samples=meta.loc[meta.tissue=='Tumor','specimen'].nunique(),exact_cell_join=True,unique_cells=True,full_gene_denominator_matches_author=True,table_S1_S2_all_samples_matched=True,table_S2_all_compartment_counts_match=all(c['table_S2_exact_count_match'] for c in checks),patient_unit='Author Table S2 Patient labels within separate tissue strata;no cross-tissue pairing inferred',p_q_computed=False,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,openpyxl=openpyxl.__version__))
    (out / 'validation.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    manifest=[]
    for name in ['annotation.csv.gz','CRC_Metadata.zip','celltypes.zip.001','celltypes.zip.002','celltypes.zip.003','celltypes.joined.zip','MOESM2.xlsx','MOESM3.xlsx','geo_family.soft.gz']:
        path=a.data_dir/name;manifest.append(dict(source_id=name,server_path=str(path),sha256=sha(path),bytes=path.stat().st_size))
    for name in ['run_khaliq35_v1.py','cell_origin_v1.py','run_cell_expression35_v2.py','analysis_spec.json']:
        path=a.out/name;manifest.append(dict(source_id=name,server_path=str(path),sha256=sha(path),bytes=path.stat().st_size))
    pd.DataFrame(manifest).to_csv(out/'source_manifest.tsv',sep='\t',index=False)
    (out/'analysis_spec.json').write_text(a.spec.read_text())
    print(json.dumps(audit),flush=True)

if __name__ == '__main__':
    main()
