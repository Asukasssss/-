"""Independent aggregate validation; no access to or recomputation of patient data."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
repo=Path(__file__).resolve().parents[1];R=repo/'results/BRCA/07_INTEGRATION/20260921T154000Z_mapping262_integration_v1'
old=pd.read_csv(repo/'results/BRCA/04_ROBUSTNESS/20260921T103911Z_camp_pair_sensitivity_v1/all117_comparison_identity_appended.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('gene')
new=pd.read_csv(R/'genes156_all_history_comparison.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('gene');assert len(new)==156 and new.index.is_unique;assert new.loc[old.index,old.columns].equals(old)
rel=pd.read_csv(R/'relations262_comparison.tsv',sep='\t');assert len(rel)==262 and rel.relation_key.is_unique and rel.gene.nunique()==150
maxerr=0.
for pre in ['CAMP','CAMP_available','FUSCC','FUSCC_warning_excluded','Tang']:
 ok=rel[pre+'_p_value'].notna();ps=rel[pre+'_p_value'].fillna(1).to_numpy() if pre.startswith('FUSCC') else rel.loc[ok,pre+'_p_value'].to_numpy();o=np.argsort(ps);qs=np.empty(len(ps));qs[o]=np.minimum(1,np.minimum.accumulate((ps[o]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]);qs=qs[ok.to_numpy()] if pre.startswith('FUSCC') else qs;err=np.max(abs(qs-rel.loc[ok,pre+'_q_value'].to_numpy()));maxerr=max(maxerr,float(err));assert err<1e-12
 assert rel.loc[~ok,pre+'_q_value'].isna().all()
camp=pd.read_csv(repo/'results/BRCA/03_PATIENT/20260921T150000Z_mapping262_patient_v1/CAMP_relations262.tsv',sep='\t');oldc=pd.read_csv(repo/'results/BRCA/04_ROBUSTNESS/20260921T103911Z_camp_pair_sensitivity_v1/association_174_sensitivities.tsv',sep='\t').set_index(['relation_id','test_family'])
for t in camp[camp.statistics_reused].itertuples():
 o=oldc.loc[(t.relation_id,t.test_family.replace('262',''))]
 for k in ['n','effect','ci_lower','ci_upper','p_value']:assert np.isclose(getattr(t,k),o[k],rtol=1e-13,atol=1e-15)
checks=0
for path in ['03_PATIENT/20260921T150000Z_mapping262_patient_v1','05_FUNCTION/20260921T152000Z_new39_literature_v1','06_EXTERNAL/20260921T151000Z_mapping262_external_v1','07_INTEGRATION/20260921T154000Z_mapping262_integration_v1']:
 p=repo/'results/BRCA'/path;d=pd.read_csv(p/'checksums.tsv',sep='\t')
 for t in d.itertuples():assert hashlib.sha256((p/t.file).read_bytes()).hexdigest()==t.sha256;checks+=1
v=dict(historical117_by249_exact_on_disk=True,current150_union156_checked=True,relations262_exact=True,old_CAMP318_reused_statistic_rows_verified=True,all5_relation_BH_families_independently_verified=True,max_BH_error=maxerr,published_file_checksums_verified=checks,raw_patient_model_refit=False)
(R/'independent_delivery_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
