"""Verify aggregate delivery and append sensitivity fields without changing source cells."""
import argparse,csv,hashlib,json,io
from pathlib import Path
R=Path(__file__).resolve().parents[2]
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):
 with p.open('x',encoding='utf-8') as f:json.dump(x,f,ensure_ascii=False,indent=2);f.write('\n')
def table(p,fields,rows):
 with p.open('x',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def key(r):return r['metabolite_key'],r['gene']
def main():
 p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);a=p.parse_args();o=R/'results/COAD/04_ROBUSTNESS'/a.run_id
 v=json.loads((o/'validation.json').read_text());assert v['status']=='PASS'
 for n,h in v['public_output_sha256'].items():assert sha(o/n)==h,n
 src=R/'results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1'
 orig=R/'results/COAD/03_PATIENT/20260919T122034Z_patient_v1/results.tsv';_,old=read(orig);old={key(r):r for r in old}
 fh,full=read(src/'catalog974_annotated_exact.tsv');nh,nominal=read(src/'nominal39_annotated_exact.tsv')
 assert len(full)==974 and len(nominal)==39
 prefix=(R/'templates/statistical_result.tsv').read_text().strip().split('\t');models={};maxerror=0.;control_max=0.
 for model in ['M1','M2','M3']:
  h,rs=read(o/(model+'.tsv'));assert h[:len(prefix)]==prefix and len(rs)==len({key(r) for r in rs})==974
  assert {key(r) for r in rs}==set(old);models[model]={key(r):r for r in rs}
  assert sum(r['test_family']!='NOT_IN_CURRENT_FAMILY' for r in rs)==674
  for name in [model,model+'_M0_CC']:
   _,rows=read(o/(name+'.tsv'));good=[r for r in rows if r['status']=='DONE'];ordered=sorted(good,key=lambda r:float(r['p_value']));best=1.
   for j in range(len(ordered)-1,-1,-1):
    rr=ordered[j];best=min(best,float(rr['p_value'])*len(ordered)/(j+1));maxerror=max(maxerror,abs(best-float(rr['q_value'])));assert abs(best-float(rr['q_value']))<1e-12
   for rr in rows:
    if rr['status']!='DONE':assert rr['effect']=='NA' and rr['p_value']=='NA' and rr['q_value']=='NA'
    else:
     assert int(float(rr['n']))==33 and -1<=float(rr['effect'])<=1 and 0<=float(rr['p_value'])<=1 and 0<=float(rr['q_value'])<=1
     if name.endswith('M0_CC'):
      assert rr['m0_cc_origin']=='IDENTICAL_INPUT_PRIMARY_REUSED'
      for k in ['effect','ci_lower','ci_upper','p_value','q_value','n']:
       e=abs(float(rr[k])-float(old[key(rr)][k]));control_max=max(control_max,e);assert e<1e-12
    assert rr['original_q']==old[key(rr)]['original_q'] and rr['original_effect']==old[key(rr)]['original_effect']
 additions=['status','reason','n','effect','ci_lower','ci_upper','p_value','q_value','delta_adjustment','delta_case_selection','loo_max_abs_delta','loo_any_sign_change','bootstrap_valid','ci_status']
 ext=[m+'_'+k for m in models for k in additions]
 def enrich(rs):return [dict(r,**{m+'_'+k:models[m][key(r)][k] for m in models for k in additions}) for r in rs]
 outfull=enrich(full);outnom=enrich(nominal)
 table(o/'catalog974_with_covariates.tsv',fh+ext,outfull);table(o/'nominal39_with_covariates.tsv',nh+ext,outnom)
 priority={'AQP9','BCAT2','GSTA4','HDC','KMT5A','PRMT7','SLC38A3','SLC6A6','UCKL1'}
 pr=[r for r in outnom if r['gene'] in priority];table(o/'priority9_relations.tsv',nh+ext,pr)
 cross=R/'results/COAD/07_INTEGRATION/20260919T142700Z_next_stage_accept_v1/shared_relations_exact.tsv';ch,cr=read(cross)
 table(o/'cross_cancer3_updated.tsv',ch+ext,enrich(cr))
 for original,after in [(full,outfull),(nominal,outnom)]:assert all(all(r[k]==a[k] for k in r) for r,a in zip(original,after))
 interpret={'priority9_relations':len(pr),'nominal39':{},'whole_family_significant':{},'original_fields_preserved':True}
 for m in models:
  good=[r for r in outnom if r[m+'_status']=='DONE'];allgood=[r for r in models[m].values() if r['status']=='DONE']
  interpret['nominal39'][m]={'evaluable':len(good),'p_lt_005':sum(float(r[m+'_p_value'])<.05 for r in good),'q_lt_005':sum(float(r[m+'_q_value'])<.05 for r in good),'same_direction_as_M0':sum(float(r[m+'_effect'])*float(r['effect'])>0 for r in good),'loo_sign_change':sum(r[m+'_loo_any_sign_change']=='True' for r in good),'max_abs_delta_adjustment':max((abs(float(r[m+'_delta_adjustment'])) for r in good),default=None)}
  interpret['whole_family_significant'][m]=[{'gene':r['gene'],'metabolite':r['metabolite_name'],'rho':float(r['effect']),'p':float(r['p_value']),'q':float(r['q_value'])} for r in allgood if float(r['q_value'])<.05]
 dump(o/'interpretation_summary.json',interpret)
 dump(o/'delivery_validation.json',{'status':'PASS','independent_BH_max_error':maxerror,'identical_M0_CC_max_error':control_max,'checks':['all server output hashes','all974 keys per model;674 in fixed family','public field prefix','independent BH on complete evaluable family','M0_CC identical primary statistics','non-estimable rows have no effect/P/q','all original974 and39 cells unchanged in annotated catalogs','priority9 and shared3 only display subsets'],'limitations':v['limitations'],'additional_output_hashes':{n:sha(o/n) for n in ['catalog974_with_covariates.tsv','nominal39_with_covariates.tsv','priority9_relations.tsv','cross_cancer3_updated.tsv','interpretation_summary.json']}})
 print(json.dumps(interpret,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
