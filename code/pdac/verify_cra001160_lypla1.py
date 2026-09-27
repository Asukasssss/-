"""Server-only independent checks; exports verification summary, never unit data."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import io
RUN=Path(__file__).resolve().parent
DATA=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/CRA001160_author_v1')
c=pd.read_csv(RUN/'private/cell_values.tsv.gz',sep='\t');u=pd.read_csv(RUN/'private/unit_values.tsv',sep='\t')
r=pd.read_csv(RUN/'public/results.tsv',sep='\t')
assert len(c)==57530 and c.cell.is_unique and c.unit.nunique()==35
assert c[c.tissue=='Tumor'].unit.nunique()==24 and c[c.tissue=='Control'].unit.nunique()==11
assert np.allclose(c.log1p_10k,np.log1p(10000*c['count']/c.library))
assert np.array_equal(c.detected,(c['count']>0).astype(float))
for sample in ['N1','T1']:
 d=DATA/'extracted'/sample
 genes=pd.read_csv(d/'genes.tsv',sep='\t',header=None)
 bars=pd.read_csv(d/'barcodes.tsv',sep='\t',header=None)[0].str.replace('-1','',regex=False)
 m=io.mmread(str(d/'matrix.mtx')).tocsc()
 j=genes.index[genes[1]=='LYPLA1'][0]
 direct=pd.DataFrame({'cell':sample+'_'+bars,'count_check':m[j,:].toarray().ravel(),'lib_check':np.asarray(m.sum(axis=0)).ravel()})
 joined=c[c.unit==sample].merge(direct,on='cell',validate='one_to_one')
 assert len(joined)==sum(c.unit==sample)
 assert np.array_equal(joined['count'],joined.count_check)
 assert np.array_equal(joined.library,joined.lib_check)
for _,row in r.iterrows():
 g=u[u.n_cells>=row.min_cells]
 a=g[(g.tissue=='Tumor')&(g.cell_type=='Ductal cell type 2')].set_index('unit')
 b=g[(g.tissue==('Tumor' if row.paired else 'Control'))&(g.cell_type=='Ductal cell type 1')].set_index('unit')
 if row.paired:
  ids=sorted(set(a.index)&set(b.index));a=a.loc[ids];b=b.loc[ids]
 assert len(a)==row.n and len(b)==row.n_reference
 assert a.n_cells.sum()==row.cells_test and b.n_cells.sum()==row.cells_reference
 if row.status=='DONE':
  assert np.isclose(a[row.metric].mean()-b[row.metric].mean(),row.effect)
  if row.paired:
   assert int((a[row.metric]>b[row.metric]).sum())==row.n_higher
primary=r[r.test_family=='LYPLA1_two_primary_contrasts'].sort_values('p_value')
assert len(primary)==2
expected=np.minimum.accumulate((primary.p_value.values*2/np.arange(1,3))[::-1])[::-1].clip(0,1)
assert np.allclose(primary.q_value,expected)
v={'status':'PASS','annotation_count_and_unique_keys':'PASS','normalization_from_saved_raw_counts':'PASS',
 'independent_csc_reextraction_two_original_matrices':'PASS','all_contrast_eligibility_counts_and_effects':'PASS',
 'paired_direction_counts':'PASS','independent_BH_two_primary_contrasts':'PASS',
 'not_validated':['No DNA-level malignancy verification','No complete original clinical patient audit','Monte Carlo P values are randomization estimates']}
(RUN/'public/independent_validation.json').write_text(json.dumps(v,indent=2))
print(json.dumps(v,indent=2))
