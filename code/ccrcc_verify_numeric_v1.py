"""Independent aggregate-matrix recomputation; public output has no patient identifiers."""
import json,sys
from pathlib import Path
import pandas as pd,numpy as np
from scipy import stats
R=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
B=R/'results/collaborative/ccRCC/B';rows=[]
for cohort,run in [('ccRCC3','20260925T151000Z_ccrcc3_v1'),('ccRCC4','20260925T155000Z_ccrcc4_histology_v2')]:
 out=B/run;m=pd.read_csv(out/'private/sample_identity_audit_private.tsv',sep='\t',dtype=str)
 x=pd.read_excel(R/f'data/candidates/camp_primary_tissue_multicancer/processed_metabolomics/PreprocessedData_{cohort}.xlsx',sheet_name='data_imputed',index_col=0);x.index=x.index.astype(str).str.strip()
 joined=m[['MetabID','SUBJECT_ID','TN']].merge(x.T,left_on='MetabID',right_index=True,validate='one_to_one')
 means=joined.drop(columns='MetabID').groupby(['SUBJECT_ID','TN']).mean(numeric_only=True)
 tumor=means.xs('Tumor',level='TN');normal=means.xs('Normal',level='TN');common=sorted(set(tumor.index)&set(normal.index));delta=tumor.loc[common]-normal.loc[common]
 r=pd.read_csv(out/'public/04_ROBUSTNESS/metabolite_paired.tsv',sep='\t').set_index('metabolite_name');check=delta[r.index]
 assert np.allclose(check.mean(0),r.effect,rtol=0,atol=1e-12)
 assert np.array_equal((check>0).sum(0),r.n_up) and np.array_equal((check<0).sum(0),r.n_down)
 exact=[]
 for name in sorted(r.index)[:12]:
  v=delta[name].to_numpy();nz=v[v!=0]
  if len(nz)==len(v) and len(np.unique(abs(nz)))==len(nz):
   p=stats.wilcoxon(v,method='exact',zero_method='wilcox',correction=False).pvalue
   tolerance=max(6*np.sqrt(p*(1-p)/100000),.00006) if r.loc[name,'p_method']=='Monte_Carlo_sign_flip_plus1' else 1e-14
   assert abs(p-r.loc[name,'p_value'])<=tolerance,(cohort,name,p,r.loc[name,'p_value'])
   exact.append(dict(feature=name,reference_exact_p=float(p),reported_p=float(r.loc[name,'p_value']),tolerance=tolerance))
 rows.append(dict(cohort=cohort,status='DONE',all_features_mean_and_direction_checked=len(r),n_pairs=len(common),independent_method='wide matrix groupby patient/tissue mean, contrasted with per-feature per-patient loop',exact_P_checks=exact))
p=Path(sys.argv[1]);p.write_text(json.dumps(dict(status='DONE',cohorts=rows),indent=2));print(json.dumps(dict(status='DONE',cohorts=[{k:v for k,v in z.items() if k!='exact_P_checks'} for z in rows])))
