from pathlib import Path
import h5py,fsspec,numpy as np,pandas as pd,json
O=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260925T103000Z_discovery_v1')
def read(q):
 if isinstance(q,h5py.Group):
  c=np.array([x.decode() if isinstance(x,bytes) else x for x in q['categories'][:]]);codes=q['codes'][:];return np.array([c[i] if i>=0 else '' for i in codes])
 return np.array([x.decode() if isinstance(x,bytes) else x for x in q[:]])
u='https://datasets.cellxgene.cziscience.com/861acfd8-25f0-418b-a445-aa96da232827.h5ad'
with fsspec.open(u,'rb',block_size=2**20).open() as f:
 with h5py.File(f,'r') as h:
  author=read(h['obs/author']);keep=author=='Neftel2019';print('selected',int(keep.sum()),'blocks',int((np.diff(np.flatnonzero(keep))>1).sum()+1),flush=True)
  for k in ['donor_id','stage','cell_type','annotation_level_2','annotation_level_3']:
   a=read(h['obs/'+k])[keep];print(k,dict(pd.Series(a).value_counts()) if k!='donor_id' else len(set(a)),flush=True)
  print('raw',list(h['raw']),'var',list(h['raw/var']),dict(h['raw/X'].attrs),flush=True)
  print('raw_examples',h['raw/X/data'][:20],flush=True);print('Xdist',h['uns/X_approximate_distribution'][()],flush=True)
