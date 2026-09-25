from pathlib import Path
import requests,json,h5py
O=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260925T103000Z_discovery_v1');p=O/'source/Darmanis2017_cellxgene.h5ad';u='https://datasets.cellxgene.cziscience.com/2e357b12-781b-4c90-9321-30206529cc4a.h5ad'
if not p.exists():
 with requests.get(u,stream=True,timeout=(20,60)) as r:
  r.raise_for_status()
  with p.with_suffix('.part').open('wb') as f:
   for b in r.iter_content(1024*1024):f.write(b)
 p.with_suffix('.part').rename(p)
with h5py.File(p,'r') as h:
 print('keys',list(h),flush=True);print('obs',list(h['obs']),flush=True);print('var',list(h['var']),flush=True);print('X',dict(h['X'].attrs),flush=True)
 for k in h['obs']:
  q=h['obs'][k]
  if isinstance(q,h5py.Group) and 'categories' in q and ('donor' not in k and 'sample' not in k):print(k,q['categories'][:].tolist()[:25],flush=True)
 print('raw',list(h['raw']) if 'raw' in h else None,flush=True)
