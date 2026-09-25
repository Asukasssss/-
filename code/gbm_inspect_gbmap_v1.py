from pathlib import Path
import requests,json,h5py,fsspec
u='https://datasets.cellxgene.cziscience.com/861acfd8-25f0-418b-a445-aa96da232827.h5ad'
with fsspec.open(u,'rb',block_size=2**20).open() as remote:
 with h5py.File(remote,'r') as h:
  print('keys',list(h),flush=True);print('obs',list(h['obs']),flush=True);print('uns',list(h['uns']),flush=True)
  for k in h['obs']:
   if any(t in k.lower() for t in ['dataset','study','source','cohort','author','annotation_level_1','stage']):
    q=h['obs'][k];print(k,q['categories'][:].tolist() if isinstance(q,h5py.Group) else 'dataset',flush=True)


