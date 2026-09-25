"""Covariate-only audit; never export patient rows."""
import argparse,hashlib,json,platform
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--resolution-dir',type=Path,required=True);p.add_argument('--run-dir',type=Path,required=True);a=p.parse_args()
o=a.run_dir.resolve();assert o.parent==a.resolution_dir.resolve().parent and not (o/'covariate_audit.json').exists()
assert (o/'.running').read_text().strip()==o.name
src=a.resolution_dir/'private_clinical_join.tsv'
c=pd.read_csv(src,sep='\t',dtype=str,keep_default_na=False)
c=c[c.TN.eq('Tumor') & c.stage.isin(['stage I','stage II','stage III','stage IV'])].copy()
assert len(c)==33 and c.individual.is_unique
age=pd.to_numeric(c.age,errors='coerce');cats={k:c[k].value_counts(dropna=False).to_dict() for k in ['gender','stage','tissue']}
models={}
for name,extra in [('M1',[]),('M2',['stage']),('M3',['tissue'])]:
 cols=['gender']+extra
 ok=age.notna() & c[cols].ne('').all(axis=1)
 z=np.column_stack([np.ones(ok.sum()),age[ok].to_numpy(dtype=float)-69,pd.get_dummies(c.loc[ok,cols],drop_first=True).to_numpy(dtype=float)])
 models[name]={'n':int(ok.sum()),'nuisance_columns':z.shape[1],'rank':int(np.linalg.matrix_rank(z)),'full_model_residual_df':int(ok.sum()-np.linalg.matrix_rank(z)-1),'min_category_n':min(int(c.loc[ok,k].value_counts().min()) for k in cols)}
r={'primary_n':33,'age':{'available':int(age.notna().sum()),'missing':int(age.isna().sum()),'min':float(age.min()),'max':float(age.max()),'median':float(age.median())},'categories':cats,'models':models,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source_path':str(src),'patient_statistics_run':False,'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}}
(o/'covariate_audit.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
(o/'.running').unlink();(o/'DONE').write_text('Covariate-only audit completed; no relation statistics.\n')
print(json.dumps(r))
