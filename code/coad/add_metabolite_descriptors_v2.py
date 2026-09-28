"""Add missing mean/rank descriptors without changing existing Wilcoxon P/q."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from patient_statistics import stable_seed
from run_source_contract_patient_v2 import ROOT,sha
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out
 assert out.parent==ROOT/'results/collaborative/COAD/B' and (out/'.running').exists();pub=out/'public';pub.mkdir(exist_ok=True)
 old=out.parent/'20260921T123200Z_paired_metabolites33_v1';check=json.loads((old/'validation.json').read_text())
 xp=ROOT/'data/candidates/camp_primary_tissue_multicancer/processed_metabolomics/PreprocessedData_COAD.xlsx';cp=out.parent/'20260919T115842Z_identity_units_v1/private_clinical_join.tsv'
 oldspec=json.loads((old/'analysis_spec.json').read_text());sources=[xp,cp]
 for p in sources:assert sha(p)==oldspec['input_sha256'][str(p.relative_to(ROOT))]
 clin=pd.read_csv(cp,sep='\t',dtype=str);t=clin[clin.TN.eq('Tumor')&clin.stage.isin(['stage I','stage II','stage III','stage IV'])];n=clin[clin.TN.eq('Normal')].set_index('individual').loc[t.individual];assert len(t)==33
 sheets=pd.read_excel(xp,sheet_name=['data','data_imputed'],index_col=0);raw=sheets['data'];full=sheets['data_imputed'];rows=[]
 spec=dict(version='COAD_metabolite_descriptor_v2',scope='159 features; primary and available subsets',new_hypothesis_tests=False,bootstrap=4000,effect='mean paired T-minus-N author normalized log scale; log base unverified',historical_P_q='exact strings retained including fixed159 sensitivity BH',historical_median_CI='retained separate fields; not relabeled mean CI')
 (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 for kind in ['primary33','observed_pairs']:
  f=old/(kind+'.tsv');assert sha(f)==check['output_sha256'][f.name];sources.append(f)
  tab=pd.read_csv(f,sep='\t',dtype=str,keep_default_na=False)
  for r in tab.to_dict('records'):
   name=r['metabolite_name'];d=full.loc[name,t.MetabID].to_numpy(float)-full.loc[name,n.MetabID].to_numpy(float)
   if kind=='observed_pairs':d=d[np.isfinite(raw.loc[name,t.MetabID].to_numpy(float))&np.isfinite(raw.loc[name,n.MetabID].to_numpy(float))]
   assert len(d)==int(r['n']) and int((d>0).sum())==int(r['up_pairs']) and int((d<0).sum())==int(r['down_pairs'])
   mean=float(d.mean()) if len(d) else np.nan;nz=d[d!=0];rank=rankdata(abs(nz));rb=float((rank[nz>0].sum()-rank[nz<0].sum())/rank.sum()) if len(nz) else 0.
   lo=hi=np.nan;seed=stable_seed('COAD_metabolite_descriptor_v2|'+kind+'|'+name)
   if r['status']=='DONE':
    rng=np.random.default_rng(seed);boot=d[rng.integers(len(d),size=(4000,len(d)))].mean(1);lo,hi=np.quantile(boot,[.025,.975])
   r.update(historical_run_id=r['run_id'],historical_effect=r['effect'],historical_effect_type=r['effect_type'],historical_ci_lower=r['ci_lower'],historical_ci_upper=r['ci_upper'],historical_test_family=r['test_family'],mean_difference=mean,rank_biserial=rb,mean_ci_lower=lo,mean_ci_upper=hi,mean_bootstrap_seed=seed,mean_majority_conflict=bool(mean*(int(r['up_pairs'])-int(r['down_pairs']))<0),pairs_higher=r['up_pairs'],pairs_lower=r['down_pairs'],pairs_equal=r['equal_pairs'])
   r.update(effect=mean if r['status']=='DONE' else np.nan,effect_type='mean_paired_T_minus_N_author_normalized_log_scale',ci_lower=lo,ci_upper=hi,run_id=out.name,analysis_version='COAD_metabolite_descriptor_v2',analysis_type='metabolite_primary' if kind=='primary33' else 'metabolite_available')
   rows.append(r)
 pd.DataFrame(rows).to_csv(pub/'metabolite_results.tsv',sep='\t',index=False,na_rep='NA')
 pd.DataFrame([dict(source_path=str(p),sha256=sha(p)) for p in sources+[Path(__file__)]]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 (pub/'validation.json').write_text(json.dumps(dict(status='PASS',rows=len(rows),all_existing_P_q_preserved=True,counts_verified=True,new_tests=0,output_sha256=sha(pub/'metabolite_results.tsv')),indent=2));(out/'.running').unlink();(out/'DONE').write_text('DONE\n')
if __name__=='__main__':main()
