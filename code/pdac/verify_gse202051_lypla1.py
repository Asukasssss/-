"""Independent server-side verification; exports aggregates only."""
from pathlib import Path
import itertools,json
import numpy as np,pandas as pd
from scipy.stats import rankdata
R=Path(__file__).resolve().parent;P=R/'public'
c=pd.read_csv(R/'private/cell_values.tsv.gz',sep='\t')
r=pd.read_csv(P/'results.tsv',sep='\t');s=pd.read_csv(P/'expression_summary.tsv',sep='\t');h=pd.read_csv(P/'histogram_source.tsv',sep='\t')
checked=0;exact=[]
for _,a in r.iterrows():
 def subset(ct):
  z=c[(c.treatment==a.treatment)&(c.level2==ct)]
  return z if a.population=='All_nuclei' else z[z.gene_count>0]
 x,y=subset('Malignant'),subset(a.reference)
 if a.unit=='pooled_nucleus':
  assert (len(x),len(y))==(a.n,a.n_reference)
  if a.status=='DONE':
   ranks=rankdata(np.r_[x.expression,y.expression]);u=ranks[:len(x)].sum()-len(x)*(len(x)+1)/2
   assert np.isclose(2*u/(len(x)*len(y))-1,a.effect,atol=1e-12)
 else:
  x=x.groupby('pid').expression.agg(['size','mean']);y=y.groupby('pid').expression.agg(['size','mean'])
  x=x[x['size']>=20];y=y[y['size']>=20];ids=sorted(set(x.index)&set(y.index))
  assert len(ids)==a.n
  if a.status=='DONE':
   d=(x.loc[ids,'mean']-y.loc[ids,'mean']).to_numpy();assert np.isclose(d.mean(),a.effect)
   if len(d)<=18:
    perm=np.array(list(itertools.product([-1,1],repeat=len(d))));p=np.mean(abs(perm@d/len(d))>=abs(d.mean())-1e-12)
    assert abs(p-a.p_value)<max(.001,5*np.sqrt(p*(1-p)/100000))
    exact.append(dict(treatment=a.treatment,reference=a.reference,population=a.population,n=int(a.n),exact_sign_permutation_p=p,recorded_monte_carlo_p=a.p_value))
 checked+=1
for fam,b in r.groupby('test_family'):
 b=b[b.p_value.notna()].sort_values('p_value');q=np.minimum.accumulate((b.p_value.to_numpy()*len(b)/np.arange(1,len(b)+1))[::-1])[::-1].clip(0,1)
 assert np.allclose(q,b.q_value)
for _,a in s.iterrows():
 b=h[(h.treatment==a.treatment)&(h.celltype==a.celltype)&(h.population==a.population)]
 assert b.n.sum()==a.n_cells
for (tr,ct),b in s.groupby(['treatment','celltype']):
 b=b.set_index('population');assert np.isclose(b.loc['All_nuclei','pooled_mean'],b.loc['LYPLA1_positive_nuclei','pooled_mean']*b.loc['All_nuclei','detection_fraction'])
v=dict(status='PASS',comparison_rows_checked=checked,bh_families_checked=r.test_family.nunique(),histogram_totals_checked=len(s),all_mean_equals_positive_mean_times_detection=True,exact_permutation_checks=exact)
(P/'independent_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v,indent=2))
