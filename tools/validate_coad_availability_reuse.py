"""Verify aggregate-only reuse and summarize prespecified sensitivity comparisons."""
import argparse,csv,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def key(r):return tuple(r[k] for k in ['cohort','metabolite_key','human_gene_id'])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sig(r):return r['status']=='DONE' and float(r['q_value'])<.05
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-id',required=True);a=ap.parse_args()
    out=(ROOT/'results/COAD/04_ROBUSTNESS'/a.run_id).resolve();assert out.parent==(ROOT/'results/COAD/04_ROBUSTNESS').resolve()
    spec=json.loads((out/'analysis_spec.json').read_text());base=ROOT/'results/COAD/03_PATIENT'/spec['base_run']
    validation=json.loads((out/'validation.json').read_text());assert validation['status']=='PASS'
    for name,h in validation['public_output_sha256'].items():assert sha(out/name)==h
    rows=read(out/'results.tsv');original=read(base/'availability33.tsv');primary=read(base/'results.tsv');all37=read(base/'all37.tsv');rna=read(base/'paired_RNA33.tsv')
    inputs=[out/'results.tsv',base/'availability33.tsv',base/'results.tsv',base/'all37.tsv',base/'paired_RNA33.tsv']
    pr={key(r):r for r in primary};old={key(r):r for r in original};a37={key(r):r for r in all37};av={key(r):r for r in rows}
    assert len(rows)==len(av)==974 and set(av)==set(pr)==set(old)==set(a37)
    prefix=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t');assert list(rows[0])[:len(prefix)]==prefix
    fields=['effect','p_value','ci_lower','ci_upper','bootstrap_valid','permutations','random_seed','n_unique_metabolite','n_unique_rna']
    counts={}
    for r in rows:
        k=key(r);counts[r['reuse_status']]=counts.get(r['reuse_status'],0)+1
        for f in ['status','reason','n','original_effect','original_q','family_n_evaluable','family_n_planned']:assert r[f]==old[k][f]
        if r['status']=='DONE':
            source=pr[k] if r['n']=='33' else old[k]
            for f in fields:assert r[f]==source[f]
        else:assert r['effect']==r['p_value']==r['q_value']=='NA'
    assert counts['IDENTICAL_INPUT_PRIMARY_REUSED']==437 and counts['REDUCED_SAMPLE_V1_REUSED']==215
    good=sorted([r for r in rows if r['status']=='DONE'],key=lambda r:float(r['p_value']));assert len(good)==652
    tail=1.0
    for i in range(651,-1,-1):
        tail=min(tail,float(good[i]['p_value'])*652/(i+1));assert math.isclose(float(good[i]['q_value']),tail,rel_tol=1e-12,abs_tol=1e-14)
    headers=[x[3:] for x in (out/'README_CN.md').read_text(encoding='utf-8').splitlines() if x.startswith('## ')]
    assert headers==['本轮问题','输入与范围','实际结果','新手解释','限制/反证','当前决定','下一步','复现命令']
    ps={k for k,r in pr.items() if sig(r)};ss={k for k,r in av.items() if sig(r)};ts={k for k,r in a37.items() if sig(r)}
    rs={r['human_gene_id'] for r in rna if sig(r)}
    summary={'status':'PASS','checks':['all hashes','all974 unique keys and common prefix','437 exact primary reuses','215 exact reduced-sample reuses','all family BH independently verified','frozen values and statuses unchanged','report section order'],
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'reuse_counts':counts,
      'primary_sig':len(ps),'all37_sig':len(ts),'availability_v2_sig':len(ss),'primary_sig_unique_genes':len({pr[k]['human_gene_id'] for k in ps}),'primary_sig_unique_features':len({pr[k]['metabolite_key'] for k in ps}),
      'primary_and_all37_sig_same_direction':sum(k in ts and float(pr[k]['effect'])*float(a37[k]['effect'])>0 for k in ps),
      'primary_and_availability_sig_same_direction':sum(k in ss and float(pr[k]['effect'])*float(av[k]['effect'])>0 for k in ps),
      'all_three_sig_same_direction':sum(k in ts and k in ss and float(pr[k]['effect'])*float(a37[k]['effect'])>0 and float(pr[k]['effect'])*float(av[k]['effect'])>0 for k in ps),
      'primary_sig_with_paired_RNA_q_lt_005':sum(pr[k]['human_gene_id'] in rs for k in ps),
      'primary_sig_loo_sign_change':sum(pr[k]['loo_any_sign_change']=='True' for k in ps),
      'primary_sig_largest_loo_absolute_change':max([float(pr[k]['loo_max_abs_delta']) for k in ps],default=None),
      'availability_q_classification_changes_from_v1':sum(sig(r)!=sig(old[key(r)]) for r in rows)}
    (out/'delivery_validation.json').write_bytes((json.dumps(summary,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(summary,ensure_ascii=False))
if __name__=='__main__':main()
