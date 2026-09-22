"""Server-only GSE166555 extraction and patient-equal source summaries."""
import argparse
import gzip
import io
import json
import platform
import re
import tarfile
from pathlib import Path
import numpy as np
import pandas as pd
from cell_origin_v1 import ROOT, sha, coverage, write_table
from run_cell_expression35_v2 import summarize, symbol_choice
import run_cell_expression35_v2 as summary_module

VERSION = 'COAD_Uhlitz_union_extract_v2'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--spec', type=Path, required=True)
    p.add_argument('--lock-commit', required=True)
    p.add_argument('--mode', choices=['extract', 'summarize'], required=True)
    a = p.parse_args()
    assert a.out.resolve().parent == ROOT / 'results/collaborative/COAD/B'
    assert (a.out / '.running').is_dir()
    spec = json.loads(a.spec.read_text()); genes = spec['targets']
    assert len(genes) == len(set(genes)) == 808
    public = a.out / 'public'; public.mkdir(exist_ok=True)
    (public / 'by_gene').mkdir(exist_ok=True)
    if a.mode == 'extract':
        columns = ['cell','cell_id','sample_id','case_id','source_id','sample_origin',
                   'main_cell_type','minor_cell_type','cell_type_imm_simple',
                   'cell_type_str_simple','Dissociation Method','nCount_RNA']
        m = pd.read_csv(a.data_dir / 'metadata.tsv.gz', sep='\t', usecols=columns, keep_default_na=False)
        assert m.cell.is_unique and (m.cell == m.cell_id).all()
        assert not m[['case_id','sample_id','sample_origin']].isin(['NA','']).any().any()
        assert m.groupby('sample_id')[['case_id','source_id','sample_origin']].nunique().max().max() == 1
        assert m.groupby('case_id').source_id.nunique().max() == m.groupby('source_id').case_id.nunique().max() == 1
        assert set(m['Dissociation Method']) == {'37C_h_TDK_1'}
        # GEO sample titles give the patient explicitly. Do not infer from barcode suffixes.
        soft = gzip.open(a.data_dir / 'family.soft.gz', 'rt').read()
        sample_info = {}
        for block in soft.split('^SAMPLE = ')[1:]:
            fields = {}
            for line in block.splitlines():
                if line.startswith('!Sample') and ' = ' in line:
                    k, v = line.split(' = ', 1); fields.setdefault(k, []).append(v)
            title = fields['!Sample_title'][0]
            match = re.fullmatch(r'Patient(\d+) (normal|tumor)(?: [12])?', title)
            assert match, title
            name = Path(fields['!Sample_supplementary_file_1'][0]).name
            sample_info[name] = (match[1], match[2].title())
        tar = tarfile.open(a.data_dir / 'counts.tar')
        members = [x for x in tar.getmembers() if x.isfile()]
        assert {x.name for x in members} == set(sample_info)
        counts = {g: np.full(len(m), -1, dtype=np.int64) for g in genes}
        total = np.zeros(len(m), dtype=np.int64); seen = np.zeros(len(m), dtype=bool)
        mapping, matrix_checks, sample_links = [], [], []
        names_by_sample = []
        choices = set(genes) | {v for g in genes for v in spec['confirmed_symbol_aliases'].get(g, [])}
        for member in members:
            stream = io.TextIOWrapper(gzip.GzipFile(fileobj=tar.extractfile(member)))
            header = stream.readline().rstrip('\r\n').split('\t')
            assert header[0] == 'gene'
            ids = header[1:]; ix = pd.Index(m.cell).get_indexer(ids)
            assert len(ids) == len(set(ids)) and (ix >= 0).all() and not seen[ix].any()
            part = m.iloc[ix]
            assert part.sample_id.nunique() == part.case_id.nunique() == part.sample_origin.nunique() == 1
            patient, tissue = sample_info[member.name]
            assert set(part.sample_origin) == {tissue}
            sample_links.append(dict(sample=part.sample_id.iloc[0], case=part.case_id.iloc[0], geo_patient=patient))
            den = np.zeros(len(ids), dtype=np.int64); names = []; selected = {}
            for line in stream:
                name, text = line.rstrip('\r\n').split('\t', 1)
                v = np.fromstring(text, sep='\t', dtype=float)
                assert len(v) == len(ids) and np.isfinite(v).all() and (v >= 0).all() and (v == np.floor(v)).all()
                den += v.astype(np.int64); names.append(name)
                if name in choices: selected[len(names)-1] = v.astype(np.int64)
            stream.close()
            assert len(names) == len(set(names)) and (den > 0).all()
            matches = den == pd.to_numeric(part.nCount_RNA).to_numpy()
            # Exact author denominator disagreement is a failed admission, not a warning.
            assert matches.all(), 'Full-gene UMI differs from source nCount_RNA'
            for gene in genes:
                j, symbol, status = symbol_choice(gene, names, spec['confirmed_symbol_aliases'])
                mapping.append(dict(gene=gene, sample=part.sample_id.iloc[0], source_symbol=symbol, mapping_status=status))
                if j is not None: counts[gene][ix] = selected[j]
            total[ix] = den; seen[ix] = True; names_by_sample.append(set(names))
            matrix_checks.append(dict(matrix_ordinal=len(matrix_checks)+1, cells=len(ids), genes=len(names),
                                      full_gene_UMI=int(den.sum()), nCount_RNA_exact_all=bool(matches.all())))
            print(json.dumps(matrix_checks[-1]), flush=True)
        tar.close(); assert seen.all()
        links = pd.DataFrame(sample_links)
        assert links.groupby('case').geo_patient.nunique().max() == links.groupby('geo_patient').case.nunique().max() == 1
        assert links.case.nunique() == 12
        links.to_csv(a.out / 'private_sample_patient_link.tsv', sep='\t', index=False)
        source_type = np.where(m.main_cell_type == 'Epithelial', 'Epithelial',
                      np.where(m.main_cell_type == 'Immune', m.cell_type_imm_simple, m.cell_type_str_simple))
        assert set(source_type) == set(spec['compartment_mapping'])
        meta = pd.DataFrame(dict(cell_id=m.cell, patient=m.case_id, specimen=m.sample_id,
            library=m.sample_id, study='Uhlitz2021_GSE166555', tissue=m.sample_origin,
            enrichment='NO_CD45_REPORTED_37C_h_TDK_1', technology='10x_3prime_VERSION_DISCREPANCY',
            lineage=[spec['compartment_mapping'][x] for x in source_type],
            midway=m.main_cell_type, subtype=m.main_cell_type+'::'+m.minor_cell_type))
        meta.to_csv(a.out / 'private_cell_metadata.tsv', sep='\t', index=False)
        np.savez_compressed(a.out / 'private_cell_counts.npz', total=total, **counts)
        pd.DataFrame(mapping).to_csv(a.out / 'private_gene_mapping.tsv', sep='\t', index=False)
        pd.DataFrame(matrix_checks).to_csv(public / 'matrix_checks.tsv', sep='\t', index=False)
        rows=[]
        for g,f in pd.DataFrame(mapping).groupby('gene',sort=True):
            ok=~f.mapping_status.str.startswith(('SOURCE_GENE_NOT','AMBIGUOUS'))
            rows.append(dict(gene=g,measured_sample_matrices=int(ok.sum()),total_sample_matrices=len(f),
                source_symbols=';'.join(sorted(set(f.loc[ok,'source_symbol']))) or 'NA',
                status='DONE' if ok.all() else 'NOT_EVALUABLE' if not ok.any() else 'PARTIAL',
                reason='Missing rows remain unmeasured;no zero or paralog substitution'))
        pd.DataFrame(rows).to_csv(public/'gene_availability.tsv',sep='\t',index=False)
        coverage({'Uhlitz2021_GSE166555':meta},public)
        audit=dict(status='PASS',cells=len(meta),patients=meta.patient.nunique(),samples=meta.specimen.nunique(),
            tumor_patients=meta.loc[meta.tissue=='Tumor','patient'].nunique(),
            normal_patients=meta.loc[meta.tissue=='Normal','patient'].nunique(),
            exact_cell_join=True,all_metadata_cells_in_counts=True,explicit_geo_patient_case_bijection=True,
            full_gene_denominator_matches_author=True,all_matrices_same_gene_universe=all(x==names_by_sample[0] for x in names_by_sample),
            source_hashes={n:sha(a.data_dir/n) for n in ['metadata.tsv.gz','counts.tar','family.soft.gz']})
        (public/'extraction_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
        print(json.dumps(audit),flush=True)
        return
    validation=json.loads((public/'extraction_validation.json').read_text())
    assert validation['status']=='PASS'
    assert all(sha(a.data_dir/n)==h for n,h in validation['source_hashes'].items())
    meta=pd.read_csv(a.out/'private_cell_metadata.tsv',sep='\t',keep_default_na=False)
    with np.load(a.out/'private_cell_counts.npz') as z: data={k:z[k] for k in z.files}
    total=data.pop('total')
    # Per-sample filtering could bias donor summaries: require full coverage or explicitly stop.
    assert all((v>=0).all() or (v<0).all() for v in data.values()), 'Partial gene coverage needs a separately specified missing-measurement analysis'
    counts={g:v for g,v in data.items() if (v>=0).all()}
    summary_module.VERSION=VERSION
    rows=summarize(meta,total,counts,genes,a.out)
    for g in genes: write_table(rows[g],public/'by_gene'/(g+'.tsv'))
    cov=pd.read_csv(public/'coverage_by_lineage.tsv',sep='\t',keep_default_na=False)
    cov['run_id']=a.out.name;cov['analysis_version']=VERSION;cov.to_csv(public/'coverage_by_lineage.tsv',sep='\t',index=False)
    validation.update(analysis_version=VERSION,lock_commit=a.lock_commit,genes_requested=35,genes_available=len(counts),
        missing_genes=sorted(set(genes)-set(counts)),rows=sum(map(len,rows.values())),p_q_computed=False,
        software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__))
    (public/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    (public/'analysis_spec.json').write_text(a.spec.read_text())
    sources=[]
    urls=spec['source_urls']
    for folder,names in [(a.data_dir,['metadata.tsv.gz','counts.tar','family.soft.gz']),
                         (a.out,['run_uhlitz35_v1.py','check_uhlitz35_v1.py','cell_origin_v1.py','run_cell_expression35_v2.py','analysis_spec.json'])]:
        for n in names:
            f=folder/n;sources.append(dict(source_id=n,server_path=str(f),url=urls.get(n,'NA:versioned code/config'),sha256=sha(f),bytes=f.stat().st_size))
    pd.DataFrame(sources).to_csv(public/'source_manifest.tsv',sep='\t',index=False)
    print(json.dumps(validation),flush=True)

if __name__=='__main__': main()
