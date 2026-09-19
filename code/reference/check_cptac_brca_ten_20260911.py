"""Ten CAMP association candidates; server-only source matrix; exploratory protein contrasts."""
import json,sys,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
sys.path.insert(0,str(ROOT/'scripts'))
from run_cptac_multicancer_inosine_v0_1 import make_protein_long,paired_effects,extract_feature_table
SRC=ROOT/'data/candidates/cptac_inosine_pan_cancer_v0.1/BRCA/pdc'
OUT=ROOT/'results/CPTAC_BRCA_ten_20260911'
OUT.mkdir(exist_ok=True)
genes=['GPCPD1','PNP','GPI','ASNS','KYNU','SLC7A11','SORD','ETNK1','PCYT2','NNMT']
files=[SRC/'pdc000120_biospecimen.json',SRC/'pdc000120_log2_ratio_matrix.tsv.gz']
sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
before={str(p):sha(p) for p in files}
m=pd.read_csv(files[1],sep='\t'); b=pd.DataFrame(json.load(open(files[0])))
assert b.aliquot_id.is_unique
suffix=[str(c).split(':',1)[1] for c in m.columns[1:] if ':' in str(c)]
assert len(suffix)==len(set(suffix))
audit,values=make_protein_long(m,b,'Breast Invasive Carcinoma',genes)
assert not values.duplicated(['gene','case_submitter_id','sample_type']).any()
effects=paired_effects(values,pd.DataFrame({'gene':genes}),'protein','PDC supplied log2_ratio; exploratory unadjusted contrast')
coverage=[]
for _,f in extract_feature_table(m,genes).iterrows():
 row={'gene':f.gene,'feature_status':f.feature_status}
 for typ,key in [('Primary Tumor','tumor'),('Solid Tissue Normal','normal')]:
  a=audit[(audit.sample_type==typ)&audit.matrix_column.notna()]
  a=a[~a.duplicated('case_submitter_id',keep=False)]
  v=pd.to_numeric(m.loc[m.iloc[:,0]==f.protein_feature_label,a.matrix_column].iloc[0],errors='coerce') if f.feature_status=='UNIQUE' else pd.Series(dtype=float)
  row[key+'_unique_cases_mapped']=len(a); row[key+'_finite']=int(np.isfinite(v).sum())
 coverage.append(row)
result=effects.merge(pd.DataFrame(coverage),on='gene',validate='one_to_one')
for g in genes:
 if result.loc[result.gene==g,'protein_status'].iloc[0]!='ANALYSED': continue
 p=values[values.gene==g].pivot(index='case_submitter_id',columns='sample_type',values='value').dropna()
 delta=p['Primary Tumor']-p['Solid Tissue Normal']
 r=result[result.gene==g].iloc[0]
 assert np.isclose(delta.mean(),r.protein_paired_mean_delta_tumor_minus_normal)
 assert np.isclose(stats.wilcoxon(delta).pvalue,r.protein_paired_wilcoxon_p_value)
valid=result[result.protein_paired_wilcoxon_p_value.notna()]
assert np.allclose(multipletests(valid.protein_paired_wilcoxon_p_value,method='fdr_bh')[1],valid.protein_paired_wilcoxon_fdr)
assert before=={str(p):sha(p) for p in files}
result.to_csv(OUT/'ten_protein_results.tsv',sep='\t',index=False)
audit.to_csv(OUT/'sample_mapping.tsv',sep='\t',index=False)
spec={'genes':genes,'study':'PDC000120','source_sha256':before,'test_family':'ten CAMP association genes; Wilcoxon primary; paired t secondary; separate BH families','pairing':'explicit PDC case and tissue types; unique exact matrix aliquot only','limitations':['unadjusted PDC processed log2 ratios; plex and tissue confounding not resolved','not CAMP metabolite replication','CAMP cohort overlap not assessed','RNA protein correlation not computed; cached RNA matrix absent','no new transformation or imputation'],'validation':'unique joins, direct mean and Wilcoxon recalculation, independent BH, unchanged source hashes PASS'}
(OUT/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
print(result.to_string(index=False))
