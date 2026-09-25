"""Refine treatment strata from linked author master metadata; reuse primary statistics."""
import argparse,json,shutil
from pathlib import Path
import numpy as np,pandas as pd
from ccrcc3_patient_v1 import corr,row,R,S
from ccrcc_stats_v1 import save,js,sha,checksums,adjust
def main(out,source,commit):
 out.resolve().relative_to(R/'results/collaborative/ccRCC/B')
 out.mkdir(exist_ok=True);(out/'private').mkdir(exist_ok=False);pub=out/'public/03_PATIENT';pub.mkdir(parents=True,exist_ok=False)
 with (out/'.running').open('x') as f:f.write('ccrcc3_treatment_sensitivity_v2')
 m=pd.read_csv(source/'private/sample_identity_audit_private.tsv',sep='\t',dtype=str)
 ap=out/'source/master_mapping.tsv';a=pd.read_csv(ap,sep='\t',dtype=str);a=a[a.METABOLON_COHORT.eq('2018_M4')]
 m=m.merge(a[['ID_METABOLON_CLIENT','ID_RNASEQ','RCC_SUBTYPE','THERAPY_EXPOSURE']],left_on='CLIENT_IDENTIFIER',right_on='ID_METABOLON_CLIENT',validate='one_to_one')
 assert len(m)==114 and m.RNAID.eq(m.ID_RNASEQ).all()
 assert m.loc[m.TN.eq('Tumor'),'RCC_SUBTYPE'].eq('clear cell').all()
 m['treatment_stratum']=m.TREATMENT+' | '+m.THERAPY_EXPOSURE
 assert m.groupby('SUBJECT_ID').treatment_stratum.nunique().max()==1
 save(m,out/'private/sample_identity_audit_private.tsv')
 coverage=m[['SUBJECT_ID','TREATMENT','THERAPY_EXPOSURE','treatment_stratum']].drop_duplicates().groupby(['TREATMENT','THERAPY_EXPOSURE','treatment_stratum']).size().reset_index(name='n_author_patients');save(coverage,pub/'treatment_source_audit.tsv')
 oldp=source/'public/03_PATIENT/association_results.tsv';old=pd.read_csv(oldp,sep='\t');ep=source/'source/direct_relations.tsv';edges=pd.read_csv(ep,sep='\t')
 rp=S/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0];rna=pd.read_csv(rp,index_col=0)
 metp=S/'processed_metabolomics/PreprocessedData_ccRCC3.xlsx';mat=pd.read_excel(metp,sheet_name='data_imputed',index_col=0);mat.index=mat.index.astype(str).str.strip()
 groups=m[m.TN.eq('Tumor')].groupby('SUBJECT_ID');patients=sorted(groups.groups);metids={k:z.MetabID.tolist() for k,z in groups};rnaids={k:z.RNAID.tolist() for k,z in groups};strata=groups.treatment_stratum.first().reindex(patients).to_numpy();rows=[]
 for _,e in edges.iterrows():
  fam='REL_TUMOR_TREATMENT_STRATIFIED';d=row(out,fam,e.gene,e.metabolite_key,e.metabolite_name);d.update(analysis_version='ccrcc3_treatment_sensitivity_v2',feature_id=e.feature_id,relation_id=e.relation_id,stable_gene_id=e.stable_gene_id,n_tumor_linked=len(patients),same_cohort_postselection=True)
  if e.gene in rna.index:
   x=np.array([mat.loc[e.metabolite_name,metids[p]].mean() for p in patients]);y=np.array([np.log2(1+rna.loc[e.gene,rnaids[p]].to_numpy(float)).mean() for p in patients]);mask=np.isfinite(x)&np.isfinite(y);d.update(corr(x[mask],y[mask],'v2|'+fam+'|'+e.relation_id,strata[mask]))
  rows.append(d)
 revised=adjust(pd.DataFrame(rows));retained=old[~old.test_family.eq('REL_TUMOR_TREATMENT_STRATIFIED')].copy();retained['statistics_reused']=True;retained['reused_from']=str(oldp);revised['statistics_reused']=False
 save(pd.concat([retained,revised],ignore_index=True),pub/'association_results.tsv');shutil.copyfile(source/'public/03_PATIENT/rna_results.tsv',pub/'rna_results.tsv')
 spec=dict(version='ccrcc3_treatment_sensitivity_v2',code_commit=commit,primary_and_available='exact original table rows reused,not recomputed',RNA='exact file copied,not recomputed',refinement='strata combine both author MetaData_M4 TREATMENT and original master_mapping THERAPY_EXPOSURE;NO includes one cabozantinib-exposed patient',unit='17 author patients;all tumor pathology clear cell and all114 RNA IDs independently concordant',method='same Spearman;9999 within refined strata permutations;4000 within-strata whole-patient bootstrap',BH='all221 current relationships,all evaluable',independent_validation=False)
 js(pub/'analysis_spec.json',spec);v=dict(status='DONE',refined_strata=coverage.to_dict('records'),relations=len(revised),p005=int(revised.p_value.lt(.05).sum()),q005=int(revised.q_value.lt(.05).sum()),RNA_reuse_sha256_equal=sha(pub/'rna_results.tsv')==sha(source/'public/03_PATIENT/rna_results.tsv'),all_67_tumor_histologies_verified=True,all114_RNA_links_verified=True)
 js(pub/'validation.json',v);save(pd.DataFrame([dict(path_or_url=str(p),sha256=sha(p)) for p in [ap,ep,oldp,source/'public/03_PATIENT/rna_results.tsv',rp,metp,Path(__file__)]]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text('''# ccRCC3治疗资料细化敏感性\n\n问题：原始MetaData_M4的NO标签是否都代表未治疗？\n\n输入：原研究master_mapping，以CLIENT_IDENTIFIER一对一连接，并用RNA ID交叉核实114标本。\n\n实际结果：67肿瘤均有clear cell病理；原NO组有1位患者在THERAPY_EXPOSURE记为Cabozantinib。采用两个作者字段的组合分层重算完整关系敏感性，计数见validation.json。\n\n新手解释：更细治疗分层仅用于敏感性，不是治疗效果检验。\n\n限制：治疗历史不完全等于采样时用药；保留两个字段，不推断未记录时间点。\n\n当前决定：主关联、可用值敏感性和RNA原数值原样复用；仅以本版替代旧粗分层敏感性。\n\n下一步：整合。\n\n复现：python ccrcc3_treatment_sensitivity_v2.py --out NEW_RUN --source ORIGINAL_CCRCC3_RUN --code-commit COMMIT。\n''',encoding='utf-8')
 checksums(pub);(out/'.running').rename(out/'.done');print(json.dumps(v),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.source,a.code_commit)
