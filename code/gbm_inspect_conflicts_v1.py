from pathlib import Path
import pandas as pd,json
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');s=r/'data/candidates/camp_primary_tissue_multicancer';o=r/'results/collaborative/GBM/A/20260925T100000Z_discovery_v1';(o/'private').mkdir(exist_ok=True)
x=s/'processed_metabolomics/PreprocessedData_GBM.xlsx';m=pd.read_csv(s/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv',dtype=str);m=m[m.Dataset.eq('GBM')].copy();sa=pd.read_excel(x,sheet_name='sampleanno',index_col=0);m['sampleanno']=m.MetabID.map(sa.GROUP);m['consistent']=m.TN.str.upper().eq(m.sampleanno)
m.to_csv(o/'private/initial_identity_audit.tsv',sep='\t',index=False)
print(pd.crosstab(m.TN,m.sampleanno,dropna=False).to_string());print('GTEx normal',int(m.loc[m.TN.eq('Normal'),'MetabID'].str.startswith('GTEX').sum()));print('GTEx labeled tumor',int(m.loc[m.TN.eq('Tumor'),'MetabID'].str.startswith('GTEX').sum()))
print('sampleanno GTEx bylabel',sa.loc[sa.index.astype(str).str.startswith('GTEX'),'GROUP'].value_counts().to_dict())
print('data duplicate feature count',pd.read_excel(x,sheet_name='data',index_col=0).index.duplicated().sum())
