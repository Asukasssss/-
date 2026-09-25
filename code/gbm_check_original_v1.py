from pathlib import Path
import pandas as pd
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');s=r/'data/candidates/camp_primary_tissue_multicancer';o=r/'results/collaborative/GBM/A/20260925T100000Z_discovery_v1'
x=o/'source/GBM_original.xlsx';print(pd.read_excel(x,sheet_name='README').to_string(index=False))
a=pd.read_excel(x,sheet_name='sampleanno');print('clinical_columns',a.columns.tolist())
for c in a.columns:
 if a[c].nunique()<12:print(c,a[c].value_counts(dropna=False).to_dict())
m=pd.read_csv(s/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv',dtype=str);m=m[m.Dataset.eq('GBM')].copy()
print('RNA_metab_exact_equal',int(m.RNAID.eq(m.MetabID).sum()),'Common_metab_equal',int(m.CommonID.eq(m.MetabID).sum()))
print('clinical exact coverage', {c:int(m[c].isin(a.case_id.astype(str)).sum()) for c in ['RNAID','MetabID','CommonID']})
for c in a.columns:
 if a[c].nunique()<12:
  v=m.MetabID.map(a.set_index('case_id')[c]);print('clinical cross',c,pd.crosstab(m.TN,v).to_dict())
