"""Independent rank-sum dynamic-program validation; server-only, aggregate output."""
from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import rankdata
from collections import Counter
root=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');run=Path(__import__('sys').argv[1]);assert run.parent==root/'results/collaborative/PDAC/B';x=root/'data/candidates/camp_primary_tissue_multicancer/processed_metabolomics/PreprocessedData_PDAC.xlsx'
t=pd.read_excel(x,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);n=pd.read_excel(x,sheet_name='metabo_imputed_filtered_Normal',index_col=0);raw=pd.read_excel(x,sheet_name='data',index_col=0);pairs=pd.read_csv(run/'private_complete_pairs.tsv',sep='\t');r=pd.read_csv(run/'public/results.tsv',sep='\t');checked=0
for _,row in r.iterrows():
 f=row.metabolite_name;a=t.loc[f,pairs.tumor_metab_id].to_numpy(float);b=n.loc[f,pairs.normal_metab_id].to_numpy(float);keep=np.isfinite(a)&np.isfinite(b)
 if row.analysis_type=='both_preimputation_available':keep &= np.isfinite(raw.loc[f,pairs.tumor_metab_id].to_numpy(float))&np.isfinite(raw.loc[f,pairs.normal_metab_id].to_numpy(float))
 assert int(keep.sum())==row.n
 if keep.sum()<6:assert row.status=='NOT_EVALUABLE' and pd.isna(row.p_value);continue
 d=(a-b)[keep];assert np.isclose(d.mean(),row.effect);d=d[d!=0];ranks=(rankdata(abs(d))*2).astype(int);obs=int(ranks[d>0].sum());total=int(ranks.sum());dist=abs(2*obs-total);counts=Counter({0:1})
 for rank in ranks:
  new=Counter(counts)
  for k,v in counts.items():new[k+int(rank)]+=v
  counts=new
 p=sum(v for k,v in counts.items() if abs(2*k-total)>=dist)/2**len(d)
 assert np.isclose(p,row.p_value,rtol=0,atol=1e-15);checked+=1
for mode,rr in r.groupby('analysis_type'):
 ps=rr.p_value.fillna(1).to_numpy();ordered=np.sort(ps);m=len(ps)
 for _,row in rr.iterrows():
  if pd.isna(row.p_value):assert pd.isna(row.q_value);continue
  expected=min([1.]+[float(ordered[j]*m/(j+1)) for j in range(m) if ordered[j]>=row.p_value])
  assert np.isclose(expected,row.q_value,rtol=0,atol=1e-14)
summary={'status':'PASS','rows_checked':len(r),'p_values_checked_independent_rank_sum_dynamic_program':checked,'BH_checked_both_fixed307_families':True,'pairwise_n_and_mean_delta_checked':True,'patient_level_data_exported':False}
(run/'public/independent_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
