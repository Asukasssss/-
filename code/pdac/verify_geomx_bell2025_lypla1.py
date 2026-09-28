"""Independent aggregate/paired-test replay. Patient data remain server-side."""
import sys,json,itertools,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import gmean
run=Path(sys.argv[1]);data=Path(sys.argv[2]);public=run/'public'
r=pd.read_csv(public/'results.tsv',sep='\t');a=pd.read_csv(run/'private/patient_pairs.tsv',sep='\t')
roi=pd.read_csv(run/'private/roi_lypla1.tsv',sep='\t',index_col=0)
x=pd.read_csv(data/'ProbeQC_merged_batches.csv',index_col=0)[roi.index]
threshold=gmean(x.loc['NegProbe-WTX'])*np.exp(np.std(np.log(x.loc['NegProbe-WTX']),ddof=1))**2
eligible=[g for g in x.index if not g.startswith('NegProbe') and sum(x.loc[g]>threshold)>10]
q3=x.loc[eligible].apply(lambda col: np.percentile(col.to_numpy(),75))
recomputed=np.log2(x.loc['LYPLA1']*gmean(q3)/q3+1)
assert np.allclose(recomputed,roi.q3_log2,rtol=0,atol=1e-12)
checks=[]
for row in r.itertuples():
    subset=a[(a['scale']==row.scale)&(a.contrast==row.contrast)]
    d=(subset.cancer-subset.reference).to_numpy();n=len(d)
    simulated=[]
    for bits in range(2**n):
        simulated.append(sum(value if bits&(1<<j) else -value for j,value in enumerate(d))/n)
    pv=sum(abs(z)>=abs(np.mean(d))-1e-12 for z in simulated)/len(simulated)
    assert abs(pv-row.p_value)<1e-12 and abs(d.mean()-row.effect)<1e-12
    assert n==row.n and sum(d>0)==row.positive_patients and sum(d<0)==row.negative_patients
    checks.append(dict(contrast=row.contrast,scale=row.scale,n=n,exact_p=pv,status='PASS'))
for scale,g in r.groupby('scale'):
    p=g.p_value.to_numpy();expected=np.array([min(1,min(len(p)*pj/sum(p<=pj) for pj in p if pj>=pi)) for pi in p])
    assert np.allclose(expected,g.q_value)
g=pd.read_csv(public/'group_summary.tsv',sep='\t');assert g.regions.sum()==95
assert len(roi)==95 and roi['Patient.Alias'].nunique()==8
assert sum(roi.Type=='ND')==14 and sum(roi.Type=='PDAC')==38
sources=pd.read_csv(public/'source_manifest.tsv',sep='\t')
assert all(hashlib.sha256((data/z.file).read_bytes()).hexdigest()==z.sha256 for z in sources.itertuples())
out=dict(status='PASS',q3_recomputed=True,exact_tests_recomputed=checks,bh_recomputed=True,roi_and_patient_counts_verified=True,source_sha256_verified=True,
    limits=['Bootstrap intervals not independently replicated','Pathology annotations trusted from authors, not independently re-read from tissue images','No pure stromal comparator'])
(public/'independent_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
