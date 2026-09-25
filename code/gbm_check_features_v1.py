from pathlib import Path
import pandas as pd
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');x=r/'data/candidates/camp_primary_tissue_multicancer/processed_metabolomics/PreprocessedData_GBM.xlsx';o=pd.read_csv(r/'results/tables/cohort_effects.tsv',sep='\t');o=o[o.dataset.eq('GBM')];a=pd.read_excel(x,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);b=pd.read_excel(x,sheet_name='metabo_imputed_filtered_Normal',index_col=0);common=set(a.index)&set(b.index)
print('old',len(o),'intersection',len(common));print('no_old',sorted(common-set(o.feature_name)));print('old_only',sorted(set(o.feature_name)-common));print('old_missing_p',o.wilcoxon_p.isna().sum())
