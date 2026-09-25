from pathlib import Path
import pandas as pd,json,subprocess
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');s=r/'data/candidates/camp_primary_tissue_multicancer'
x=s/'processed_metabolomics/PreprocessedData_GBM.xlsx';m=pd.read_csv(s/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv',dtype=str);m=m[m.Dataset.eq('GBM')]
a=pd.read_excel(x,sheet_name='metanno');print('annotation',a.head(6).to_json(orient='records'));sa=pd.read_excel(x,sheet_name='sampleanno',index_col=0);print('sampleanno_cols',sa.columns.tolist(), 'types',sa.iloc[:,0].value_counts().to_dict())
for name in m.RNAFile.unique():
 p=s/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/name
 d=pd.read_csv(p,index_col=0);print('RNA',d.shape,'index_unique',d.index.is_unique,'gene_examples',d.index[:5].tolist(),'ID_coverage',m.RNAID.isin(d.columns).sum(),'range',d.min().min(),d.max().max(),'column_sums_range',d.sum().min(),d.sum().max())
for sh in ['data','data_imputed','metabo_imputed_filtered_Tumor']:
 d=pd.read_excel(x,sheet_name=sh,index_col=0);print(sh,'range',d.min().min(),d.max().max())
print('fields_unique', {c:int(m[c].nunique()) for c in ['CommonID','MetabID','RNAID','Identifier']})
print(subprocess.run(['Rscript','-e','cat(requireNamespace("edgeR",quietly=TRUE))'],capture_output=True,text=True).stdout)
