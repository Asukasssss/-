from pathlib import Path
import tarfile,pandas as pd
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');s=r/'data/candidates/camp_primary_tissue_multicancer';o=r/'results/collaborative/GBM/A/20260925T100000Z_discovery_v1/source';p=o/'GBM_original.xlsx'
if not p.exists():
 with tarfile.open(s/'pancancer_metabolomics_v.0.3.4.tar.gz','r|gz') as t:
  for member in t:
   if member.name=='pancancer_metabolomics/data/metabolomics_original/GBM.xlsx':
    p.write_bytes(t.extractfile(member).read());break
x=pd.ExcelFile(p);print(x.sheet_names,flush=True)
for sh in x.sheet_names:
 d=pd.read_excel(x,sheet_name=sh);print(sh,d.shape,'columns',d.columns[:12].tolist(),flush=True)
