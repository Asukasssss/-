"""Export only server-produced aggregate shards; validate exact v1 reuse."""
import argparse, csv, hashlib, json, shutil
from pathlib import Path
import pandas as pd
from cell_origin_v1 import PREFIX

RUN = '20260920T081200Z_cell_expression35_v2'
OLD = '20260920T072000Z_cell_expression_v1'
ROOT = Path(__file__).resolve().parents[2]
KEY = ['study','tissue','enrichment','technology','annotation_level','cell_type','support_set']

def read(p):
    with p.open(encoding='utf-8',newline='') as f:
        r=csv.DictReader(f,delimiter='\t'); return r.fieldnames,list(r)

def write(p,header,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=header,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,required=True);a=parser.parse_args()
    spec=json.loads((ROOT/'config/coad_cell_origin35_v2.json').read_text())
    out=ROOT/'results/COAD/06_EXTERNAL'/RUN;out.mkdir(exist_ok=True);(out/'by_gene').mkdir(exist_ok=True)
    overview=[];availability=[];checks={};nrows=0
    for study in ['Lee','Pelka']:
        checks[study]=json.loads((a.input/('public_'+study)/'validation.json').read_text())
        assert checks[study]['status']=='PASS' and checks[study]['genes_requested']==35
        _,rows=read(a.input/('public_'+study)/'gene_availability.tsv')
        availability.extend([dict(study=study,**r) for r in rows])
        shutil.copy2(a.input/('public_'+study)/'validation.json',out/(study+'_validation.json'))
    for gene in spec['targets']:
        combined=[];header=None
        for study,expected in [('Lee',148),('Pelka',2304)]:
            h,rows=read(a.input/('public_'+study)/'by_gene'/(gene+'.tsv'))
            assert len(rows)==expected and (header is None or h==header);header=h
            assert all(r['gene']==gene and r['p_value']==r['q_value']=='NA' for r in rows)
            assert len({tuple(r[k] for k in KEY) for r in rows})==len(rows)
            if gene in spec['reuse_targets']:
                oh,old=read(ROOT/'results/COAD/06_EXTERNAL'/OLD/(study+'_source_expression_summary.tsv'))
                old={tuple(r[k] for k in KEY):r for r in old if r['gene']==gene}
                for r in rows:
                    prior=old[tuple(r[k] for k in KEY)]
                    assert all(r[k]==prior[k] for k in oh if k not in ['run_id','analysis_version'])
            for r in rows:
                assert int(r['n'])<=int(r['n_reference'])
                if r['support_set']=='patients_with_at_least_20_cells':
                    assert int(r['n'])==int(r['n_patients_ge20']) and int(r['n_cells'])>=20*int(r['n'])
                if int(r['n'])<3 or gene in checks[study]['missing_new_genes']:
                    assert r['effect']=='NA' and r['detection_fraction_median']=='NA' and r['status']=='NOT_EVALUABLE'
                if r['detection_fraction_median']!='NA':
                    assert int(r['n'])>=3
                    for metric in ['detection_fraction','pseudobulk_CPM','target_UMI','library_UMI']:
                        v=[float(r[metric+'_'+s]) for s in ['min','q25','median','q75','max']]
                        assert v==sorted(v) and min(v)>=0
                        if metric=='detection_fraction':assert max(v)<=1
            assert sum(int(r['n_cells']) for r in rows if r['annotation_level']=='lineage' and r['support_set']=='all_covered_patients')==checks[study]['cells']
            combined.extend(rows)
        combined.sort(key=lambda r:tuple(r[k] for k in KEY));write(out/'by_gene'/(gene+'.tsv'),header,combined);nrows+=len(combined)
        f=pd.DataFrame(combined)
        chosen=f[f.tissue.isin(['Tumor','T']) & (f.annotation_level=='lineage') & (f.support_set=='patients_with_at_least_20_cells') & f.enrichment.isin(['unsorted','NO_CD45_REPORTED_FICOLL_PURIFIED'])].copy()
        for key,g in chosen.groupby(['study','technology'],sort=True):
            valid=g[(g.status=='DONE') & (g.pseudobulk_CPM_median!='NA')].copy()
            valid['value']=pd.to_numeric(valid.pseudobulk_CPM_median)
            valid=valid.sort_values(['value','cell_type'],ascending=[False,True])
            top=valid.iloc[0] if len(valid) else None
            top_names='NA' if top is None else ('NO_POSITIVE_MEDIAN' if top.value==0 else '|'.join(valid.loc[valid.value==top.value,'cell_type']))
            display=dict.fromkeys(PREFIX,'NA')
            display.update(cancer='COAD',cohort=key[0],stage_id='06_EXTERNAL',run_id=RUN,
                analysis_version=spec['analysis_version'],analysis_type='descriptive_lineage_maximum',gene=gene,
                unit='patient',n=top.n if top is not None else 'NA',effect_type='median_patient_pseudobulk_CPM',
                effect=top.pseudobulk_CPM_median if top is not None else 'NA',status='DONE' if top is not None else 'NOT_EVALUABLE',
                reason='Display of existing summary; no independent evidence; not cell-type specificity',source_id=key[0])
            display.update(dict(study=key[0],technology=key[1],eligible_lineages=len(valid),
                highest_median_CPM_lineage=top_names,
                n_qualified_patients=top.n if top is not None and top.value>0 else 'NA',
                median_CPM=top.pseudobulk_CPM_median if top is not None else 'NA',
                median_detection=top.detection_fraction_median if top is not None else 'NA',
                interpretation='Descriptive maximum among sufficiently covered lineages; not specificity, significance, or function' if top is not None else 'NOT_EVALUABLE: gene unavailable or insufficient coverage'))
            overview.append(display)
    assert nrows==85820
    pd.DataFrame(availability).fillna('NA').sort_values(['gene','study']).to_csv(out/'gene_availability.tsv',sep='\t',index=False)
    pd.DataFrame(overview).to_csv(out/'tumor_lineage_overview.tsv',sep='\t',index=False)
    shutil.copy2(ROOT/'config/coad_cell_origin35_v2.json',out/'analysis_spec.json')
    validation=dict(status='PASS',genes=35,rows=nrows,per_gene_rows=2452,exact_v1_reuse=['HDC','GSTA4'],
        checks=['unique stratum keys','35 genes retained','all p/q NA','quantile order and bounds','privacy minimum3','missing gene all NA','coverage conservation','20-cell support rule','HDC/GSTA4 old fields string-identical'],studies=checks)
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    files=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(out.rglob('*')) if p.is_file()]
    assert all(r['bytes']<5_000_000 for r in files)
    (out/'aggregate_manifest.json').write_text(json.dumps(files,indent=2)+'\n')
    print(json.dumps(validation))

if __name__=='__main__':main()
