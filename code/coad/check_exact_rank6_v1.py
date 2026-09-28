"""Synthetic validation only; no real patient data."""
import json,random,sys
from pathlib import Path
from exact_rank6_v1 import exact_test,bh
import numpy as np
import scipy
from scipy.stats import rankdata,permutation_test,false_discovery_control

random.seed(20260921)
cases=[(list(range(6)),list(range(6))),([0,0,1,2,3,4],[0,1,1,2,4,3])]
while len(cases)<22:
 x=[random.randrange(4)for _ in range(6)];y=[random.randrange(5)for _ in range(6)]
 if len(set(x))>1 and len(set(y))>1:cases.append((x,y))
for x,y in cases:
 got=exact_test(x,y)
 a,b=rankdata(x),rankdata(y)
 def statistic(v):return abs(float(np.corrcoef(v,b)[0,1]))
 ref=permutation_test((a,),statistic,permutation_type='pairings',alternative='greater',n_resamples=np.inf,vectorized=False)
 assert abs(got['p']-ref.pvalue)<1e-12
 assert abs(got['rho']-float(np.corrcoef(a,b)[0,1]))<1e-12
assert exact_test(*cases[0])['p']==2/720
for _ in range(20):
 ps=[random.randrange(721)/720 for _ in range(30)]
 assert np.allclose(bh(ps),false_discovery_control(ps),rtol=0,atol=1e-15)
try:exact_test([1]*6,list(range(6)))
except ValueError:pass
else:raise AssertionError('Constant must fail')
record=dict(status='PASS',scope='Synthetic only',permutation_cases=len(cases),bh_cases=20,scipy=scipy.__version__,numpy=np.__version__,real_patient_tests=0,tail='Absolute statistic greater-or-equal; NOT scipy default doubled smaller tail',ties='All 720 index permutations retain multiplicity',numerics='Integer centered doubled ranks: exact tail comparison')
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
