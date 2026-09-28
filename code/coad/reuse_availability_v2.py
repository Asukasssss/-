"""Apply predeclared identical-input reuse to public aggregate results on server."""
import argparse,csv,hashlib,json
from datetime import datetime,timezone
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def dump(p,x):
    with p.open('x',encoding='utf-8') as f:json.dump(x,f,ensure_ascii=False,indent=2);f.write('\n')
def table(p,fields,rows):
    with p.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def key(r):return tuple(r[k] for k in ['cohort','metabolite_key','human_gene_id'])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base-dir',type=Path,required=True);ap.add_argument('--run-dir',type=Path,required=True);ap.add_argument('--code-commit',required=True);a=ap.parse_args()
    base,out=a.base_dir.resolve(),a.run_dir.resolve()
    assert base.parent==out.parent and out.parent.parts[-4:]==('results','collaborative','COAD','B')
    assert (base/'DONE').exists() and (out/'.running').read_text().strip()==out.name
    spec=json.loads((out/'locked_spec.json').read_text());assert base.name==spec['base_run']
    v=json.loads((base/'validation.json').read_text());assert v['status']=='PASS'
    inputs=[base/'results.tsv',base/'availability33.tsv',base/'validation.json',base/'analysis_spec.json',out/'locked_spec.json',Path(__file__)]
    hashes={str(p):sha(p) for p in inputs}
    for name in ['results.tsv','availability33.tsv','analysis_spec.json']:assert sha(base/name)==v['public_output_sha256'][name]
    srcspec=json.loads((base/'analysis_spec.json').read_text());assert srcspec['primary_population'].endswith('n=33 unique individuals')
    src={key(r):r for r in read(base/'results.tsv')};rows=read(base/'availability33.tsv');assert len(src)==len(rows)==974
    fields=list(rows[0])+['reuse_status','source_analysis_version','source_run_id']
    reuse_fields=['effect','p_value','ci_lower','ci_upper','bootstrap_valid','permutations','random_seed','n_unique_metabolite','n_unique_rna']
    reused=reduced=0
    for r in rows:
        original=src[key(r)]
        r.update(run_id=out.name,analysis_version=spec['version'],source_analysis_version=srcspec['version'],source_run_id=base.name,reuse_status='NOT_EVALUABLE_OR_OUTSIDE_FAMILY')
        if r['test_family']!='NOT_IN_CURRENT_FAMILY':r['test_family']=spec['version']+'_availability33'
        if r['status']=='DONE':
            assert original['status']=='DONE' and original['n']=='33'
            if r['n']=='33':
                assert r['effect']==original['effect']
                for name in reuse_fields:r[name]=original[name]
                r['reuse_status']='IDENTICAL_INPUT_PRIMARY_REUSED';reused+=1
            else:r['reuse_status']='REDUCED_SAMPLE_V1_REUSED';reduced+=1
    assert reused==spec['expected_identical_input_rows']==437 and reduced==spec['expected_reduced_sample_rows']==215
    good=[r for r in rows if r['status']=='DONE'];ranked=sorted(good,key=lambda r:float(r['p_value']));tail=1.0
    for i in range(len(ranked)-1,-1,-1):
        tail=min(tail,len(ranked)*float(ranked[i]['p_value'])/(i+1));ranked[i]['q_value']=str(tail)
    spec.update(run_id=out.name,code_commit=a.code_commit,created_utc=datetime.now(timezone.utc).isoformat(),input_hashes=hashes)
    dump(out/'analysis_spec.json',spec);table(out/'results.tsv',fields,rows)
    sig=[r for r in good if float(r['q_value'])<.05]
    summary=dict(planned=674,evaluable=len(good),catalog_rows=len(rows),reused_primary=reused,reused_reduced_sample=reduced,q_lt_005=len(sig),significant_positive=sum(float(r['effect'])>0 for r in sig),significant_negative=sum(float(r['effect'])<0 for r in sig),significant_unique_genes=len({r['human_gene_id'] for r in sig}))
    dump(out/'summary.json',summary)
    table(out/'source_manifest.tsv',['source_path','sha256'],[dict(source_path=p,sha256=h) for p,h in hashes.items()])
    for p,h in hashes.items():assert sha(Path(p))==h
    dump(out/'validation.json',dict(status='PASS',checks=['base batch completed and hashed','437 identical-input statistics reused exactly','215 reduced-sample statistics reused','BH recomputed across same652 evaluable tests','all974 rows retained','all inputs unchanged'],public_output_sha256={n:sha(out/n) for n in ['results.tsv','summary.json','analysis_spec.json','source_manifest.tsv']}))
    (out/'.running').unlink();(out/'DONE').write_text(datetime.now(timezone.utc).isoformat()+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
