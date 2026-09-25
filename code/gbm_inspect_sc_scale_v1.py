from pathlib import Path
import h5py,numpy as np
p=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260925T103000Z_discovery_v1/source/Darmanis2017_cellxgene.h5ad')
with h5py.File(p,'r') as h:
 print('uns',list(h['uns']));print('layers',list(h['layers']))
 for k in h['uns']:
  q=h['uns'][k]
  if isinstance(q,h5py.Dataset):print(k,str(q[()])[:1400])
 v=h['X/data'][:];print('X_range',v.min(),v.max(),'integer_fraction',np.mean(v==np.floor(v)))
 print('donor_count',len(h['obs/donor_id/categories']));print('feature_examples',h['var/feature_name/categories'][:5] if isinstance(h['var/feature_name'],h5py.Group) else h['var/feature_name'][:5])
