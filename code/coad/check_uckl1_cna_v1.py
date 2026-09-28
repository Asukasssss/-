"""Independent arithmetic checks of server-private paired inputs and public results."""
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
d=pd.read_csv(a.out/'private_donor_values.tsv',sep='\t');d=d[d.n_cells>=20]
r=pd.read_csv(a.out/'public/results.tsv',sep='\t');checks=[]
assert len(r)==3 and not d.duplicated(['patient','group']).any()
for x in r.to_dict('records'):
    z=d.pivot(index='patient',columns='group',values='mean_log1p_CP10K')[[x['group_a'],x['group_b']]].dropna()
    delta=(z.iloc[:,0]-z.iloc[:,1]).to_numpy();n=len(delta)
    assert n==x['n'] and np.isclose(np.median(delta),x['effect'],rtol=0,atol=1e-12)
    up=int(sum(delta>0));down=int(sum(delta<0));nz=up+down
    exact=min(1.,2*sum(math.comb(nz,k) for k in range(min(up,down)+1))/2**nz) if nz else 1.
    assert up==x['n_up'] and down==x['n_down'] and n-nz==x['n_equal']
    if x['status']=='DONE':
        assert np.isclose(exact,x['p_value'],rtol=0,atol=1e-12)
        assert x['ci_lower']<=x['effect']<=x['ci_upper'] and x['ci_coverage']>=.95
    checks.append(dict(contrast=x['contrast'],n=n,independent_sign_probability=exact,status='PASS'))
q=multipletests(r.p_value.fillna(1).to_numpy(),method='fdr_bh')[1]
assert np.allclose(q,r.q_value.to_numpy(),atol=1e-12,rtol=0)
obj=dict(status='PASS',checks=checks,independent_BH='statsmodels fdr_bh',
         paired_effect_and_sign_counts_recomputed=True,patient_values_stayed_on_server=True)
(a.out/'public/independent_validation.json').write_text(json.dumps(obj,indent=2)+'\n')
print(json.dumps(obj))
