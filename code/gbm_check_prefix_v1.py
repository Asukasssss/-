from pathlib import Path
import pandas as pd
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');s=r/'data/candidates/camp_primary_tissue_multicancer';m=pd.read_csv(s/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv',dtype=str);m=m[m.Dataset.eq('GBM')]
for tn,z in m.groupby('TN'):print(tn,'prefixes',z.MetabID.str.split('-').str[0].value_counts().to_dict())
