"""Read-only server inspection; print schema and aggregate coverage only."""
import json, tarfile, io
from pathlib import Path
import pandas as pd
R=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
S=R/'data/candidates/camp_primary_tissue_multicancer'
m=pd.read_csv(S/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv',dtype=str)
print('MAPPING_COLUMNS',m.columns.tolist())
for cohort in ['ccRCC3','ccRCC4']:
 d=m[m.Dataset.eq(cohort)]
 print('MAPPING',cohort,len(d),d.TN.value_counts().to_dict(), 'RNA_FILES',d.RNAFile.unique().tolist())
 x=pd.ExcelFile(S/f'processed_metabolomics/PreprocessedData_{cohort}.xlsx')
 print('PROCESSED',cohort,x.sheet_names)
 for sheet in ['sampleanno','data','data_imputed']:
  a=pd.read_excel(x,sheet_name=sheet,index_col=0)
  print(sheet,a.shape,'columns' if sheet=='sampleanno' else 'unique_axes',a.columns.tolist() if sheet=='sampleanno' else [a.index.is_unique,a.columns.is_unique])
with tarfile.open(S/'pancancer_metabolomics_v.0.3.4.tar.gz') as t:
 names=[n for n in t.getnames() if ('ccRCC' in n or 'KIRC' in n) and not n.endswith('/')]
 print('ARCHIVE_FILES',names)
 for name in names:
  if 'metabolomics_original' in name and name.endswith('.xlsx'):
   x=pd.ExcelFile(io.BytesIO(t.extractfile(name).read()))
   print('ORIGINAL',name,x.sheet_names)
   for sheet in x.sheet_names:
    if any(k in sheet.lower() for k in ['sample','clinical','metadata']):
     a=pd.read_excel(x,sheet_name=sheet)
     print('METADATA',sheet,a.shape,a.columns.tolist())
     for c in a:
      if any(k in str(c).lower() for k in ['patient','case','tissue','group','region','tumor','normal']):
       print('FIELD',c,'unique',a[c].nunique(),'missing',int(a[c].isna().sum()),'max_repeat',int(a[c].value_counts().max()) if a[c].notna().any() else 0)
