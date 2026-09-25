"""Inspect depositor myeloid labels without exporting cell identifiers."""
import json
from pathlib import Path
import anndata as ad
import pandas as pd
src=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/PRAD/B/20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad')
a=ad.read_h5ad(src,backed='r');o=a.obs
print('VERSIONS',ad.__version__,pd.__version__)
print('COLUMNS',list(o.columns));print('UNS_KEYS',list(a.uns.keys()))
m=o[o.celltype_major_v2.astype(str).eq('Myeloid')]
print('MYELOID_CELLS',len(m))
for c in o.columns:
 if any(w in c.lower() for w in ['celltype','cell_type','cluster','annot','leiden','louvain']):
  print('FIELD',c,'MYELOID_LEVELS',json.dumps(m[c].astype(str).value_counts().to_dict()))
print('TOTAL_TISSUES',m['type'].value_counts().to_dict())
for field in ['celltype_minor_v2','celltype_subset_v2']:
 counts=m.groupby(['donor_id','type',field],observed=True).size().unstack('type',fill_value=0)
 for threshold in [20,10]:
  ok=counts[(counts.cancer>=threshold)&(counts.adj_benign>=threshold)]
  print('PAIRED_COVERAGE',field,threshold,ok.groupby(field,observed=True).size().to_dict())
for k in a.uns:
 if any(w in k.lower() for w in ['citation','title','description','annotation','method']):print('UNS',k,str(a.uns[k])[:4000])
a.file.close()
