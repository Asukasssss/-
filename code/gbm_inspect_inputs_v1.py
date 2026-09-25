from pathlib import Path
import pandas as pd,json
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
s=r/'data/candidates/camp_primary_tissue_multicancer'
m=pd.read_csv(s/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv',dtype=str);g=m[m.Dataset.eq('GBM')]
print('mapping_columns',list(m.columns));print('GBM_n',len(g))
for c in ['TN','RNAFile','MetabFile','MetabFile_sheet','Dataset']:
 print(c,g[c].value_counts(dropna=False).to_dict())
x=pd.ExcelFile(s/'processed_metabolomics/PreprocessedData_GBM.xlsx');print('sheets',x.sheet_names)
for sh in x.sheet_names:
 d=pd.read_excel(x,sheet_name=sh,index_col=0); print('sheet',sh,'shape',d.shape,'first_feature_keys',d.index[:5].tolist(),'missing',int(d.isna().sum().sum()))
print('mapping GBM missing RNA',int(g.RNAID.isna().sum()))
for p in (r/'results/CAMP_Phase1_v1.0').rglob('*'):
 if p.is_file() and ('effect' in p.name or 'manifest' in p.name): print('frozen',p)
