"""Independent aggregate/statistic checks; only cohort output is published."""
import json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
p=Path(sys.argv[1]);pub=p/'public'
r=pd.read_csv(pub/'results.tsv',sep='\t');s=pd.read_csv(p/'patient_region_summary.tsv',sep='\t');q=pd.read_csv(p/'sample_qc.tsv',sep='\t')
assert len(q)==6 and q.patientid.nunique()==6
assert s.n_spots.sum()==q.n_spots.sum()
checks=[]
for i,row in r.iterrows():
 a=s[(s.n_spots>=row.min_spots)&(s.region=='Tumor_containing')].set_index('patientid')
 b=s[(s.n_spots>=row.min_spots)&(s.region==row.reference)].set_index('patientid')
 ix=a.index.intersection(b.index);delta=a.loc[ix,'mean_log1p_CP10k']-b.loc[ix,'mean_log1p_CP10k']
 assert len(ix)==row.n and (delta>0).sum()==row.n_up and (delta<0).sum()==row.n_down
 if len(ix): assert np.isclose(delta.mean(),row.effect)
 if len(ix)>=3:
  n=(delta!=0).sum();k=min((delta>0).sum(),(delta<0).sum());exact=min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n)
  assert np.isclose(exact,row.p_value)
  assert row.ci_lower<=row.ci_upper
 # Library-weighted pseudobulk direction is an additional descriptive normalization sensitivity.
 pdiff=a.loc[ix,'pseudobulk_CP10k']-b.loc[ix,'pseudobulk_CP10k']
 checks.append(dict(reference=row.reference,min_spots=int(row.min_spots),n=len(ix),n_pseudobulk_up=int((pdiff>0).sum()),n_pseudobulk_down=int((pdiff<0).sum()),same_direction_as_primary_metric=int(((pdiff>0)==(delta>0)).sum())))
 r.loc[i,'test_family']='primary_Tumor_vs_Stroma' if row.reference=='Stroma' else 'secondary_4_region_contrasts_exploratory'
 r.loc[i,'family_n_evaluable']=int(((r.min_spots==row.min_spots)&(r.reference!='Stroma')&(r.n>=3)).sum()) if row.reference!='Stroma' else int(row.n>=3)
r.to_csv(pub/'results.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame(checks).to_csv(pub/'normalization_sensitivity.tsv',sep='\t',index=False)
val=json.loads((pub/'validation.json').read_text());val.update(independent_recomputed_region_effects=True,independent_combinatorial_sign_P=True,region_spot_totals_reconcile=True,normalization_sensitivity='pseudobulk ratio directions retained; no extra significance tests')
(pub/'validation.json').write_text(json.dumps(val,ensure_ascii=False,indent=2))
print(r[r.min_spots==20][['reference','n','n_up','effect','p_value']].to_string(index=False));print(pd.DataFrame(checks).to_string(index=False));print('AUDIT_PASS')
