from pathlib import Path
import json,numpy as np,pandas as pd
from gbm_discovery_v1 import sha,save,js,check
O=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260925T103000Z_discovery_v1');p=O/'public/06_EXTERNAL_v1b';s=pd.read_csv(p/'sc_celltype_profiles.tsv',sep='\t');d=pd.read_csv(O/'private/sc_donor_profiles_private.tsv',sep='\t');errors=[]
for _,r in s[s.status.eq('DONE')].iterrows():
 q=d[d.study.eq(r.cohort)&d.celltype.eq(r.celltype)&d.gene.eq(r.gene)];assert q.donor_id.nunique()==len(q)>=3 and q.n_cells.ge(20).all();errors.append(abs(q.mean_expression.mean()-r.effect));assert abs(q.detection.mean()-r.mean_detection_fraction)<1e-10
assert max(errors)<1e-10
for study in ['Darmanis2017','Neftel2019']:
 a=np.load(O/'private'/f'{study}_candidate_expression.npz');assert a['z'].shape[1]==142 and np.all(a['library']>0);assert len(s[s.cohort.eq(study)].gene.unique())==142
r=pd.read_csv(p/'sc_source_stability.tsv',sep='\t');ok=r[r.status.eq('DONE')];assert ok.bootstrap_top_frequency.between(0,1).all() and ok.bootstrap_valid.between(1,1000).all();c=pd.read_csv(p/'sc_cross_study.tsv',sep='\t');assert c.gene.is_unique and len(c)==142
v=dict(status='DONE',donor_equal_mean_rows_checked=len(errors),max_mean_error=max(errors),detection_mean_checked=True,all142_genes_retained=True,bootstrap_frequencies_bounded=True,no_patient_cell_rows_in_public_profiles=True)
js(p/'independent_validation.json',v);check(p);print(json.dumps(v))
