"""Audit pooled-cell ranks and tie-corrected normal approximation independently."""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,norm
R=Path(sys.argv[1]);d=np.load(R/'private/Wu2021_cell_values.npz');m=d['malignant'];t=pd.read_csv(R/'public/results.tsv',sep='\t');out=[]
for metric,key in [('LYPLA1','expression'),('KEGG92_score','score')]:
 v=d[key];n=int(m.sum());b=int((~m).sum());N=n+b;r=rankdata(v);u=r[m].sum()-n*(n+1)/2
 _,counts=np.unique(v,return_counts=True);ties=np.sum(counts.astype(float)**3-counts)
 sigma=np.sqrt(n*b/12*((N+1)-ties/(N*(N-1))));z=(abs(u-n*b/2)-.5)/sigma;p=2*norm.sf(z);logp=(np.log(2)+norm.logsf(z))/np.log(10)
 row=t[(t.cohort=='Wu2021')&(t.metric==metric)].iloc[0]
 assert np.isclose(2*u/(n*b)-1,row.effect,atol=1e-12)
 assert np.isclose(p,row.p_value,rtol=1e-8,atol=1e-300)
 assert n==row.n and b==row.n_reference
 out.append(dict(metric=metric,n=n,n_reference=b,rank_biserial=2*u/(n*b)-1,p_value=p,log10_p=logp,status='PASS'))
pd.DataFrame(out).to_csv(R/'public/independent_rank_audit.tsv',sep='\t',index=False)
j=json.loads((R/'public/validation.json').read_text());j['independent_recalculation']='PASS ranks, effect, tie correction, continuity correction, P, cell counts';(R/'public/validation.json').write_text(json.dumps(j,indent=2))
print(pd.DataFrame(out).to_string(index=False))
