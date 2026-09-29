from pathlib import Path
import h5py,numpy as np,pandas as pd
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');C=B/'results/collaborative'
files={'BRCA':C/'BRCA/A/20260919T140105Z_scRNA117_v1/source/wu_curated.h5ad','PDAC':B/'data/candidates/SCP1089/GSE202051_totaldata-final-toshare.h5ad','PRAD':C/'PRAD/B/20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad'}
for c,p in files.items():
 with h5py.File(p) as h:
  print(c,'obs',list(h['obs']),flush=True)
  for k in ['X','raw/X','var','raw/var','layers']:
   if k in h:
    g=h[k];print(c,k,dict(g.attrs),list(g)[:12] if isinstance(g,h5py.Group) else g.shape,flush=True)
p=C/'COAD/B/20260922T141755Z_source_contract_sc_v2/private_cell_counts.npz'
z=np.load(p,allow_pickle=True);print('COAD NPZ',[(k,z[k].shape) for k in z.files][:20],flush=True)
for p in [C/'BRCA/A/20260923T040012Z_epithelial_paired20_v1/public/results.tsv',C/'COAD/B/20260925T140257Z_uckl1_cna_v1/public/results.tsv',C/'PRAD/B/20260925T133653Z_slc6a6_sc_paired_v1/public/06_EXTERNAL/results.tsv']:
 if p.exists():
  d=pd.read_csv(p,sep='\t');print(str(p),d[d.gene.isin(['ASNS','UCKL1','SLC6A6'])].to_csv(sep='\t',index=False),flush=True)
for c,p in [('PDAC',C/'PDAC/B/20260928T031629Z_gse202051_lypla1_v1/private/cell_values.tsv.gz'),('PRAD',C/'PRAD/B/20260928T101700Z_lypla1_erickson_v1/private/spot_values.tsv.gz')]:print(c,'cache columns',list(pd.read_csv(p,sep='\t',nrows=0)),flush=True)
