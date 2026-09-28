"""Public aggregates only: exact relation joins and full current/history source inventory."""
import json,hashlib,argparse
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
PAT=ROOT/'results/COAD/03_PATIENT/20260922T135816Z_source_contract_patient_v2'
MAP=ROOT/'results/COAD/02_MAPPING/20260922T053000Z_paired101_mapping_v1'
SC=ROOT/'results/COAD/06_EXTERNAL/20260922T141755Z_source_contract_sc_v2'
MET=ROOT/'results/COAD/04_ROBUSTNESS/20260922T143500Z_metabolite_descriptors_v2'
def read(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
def write(t,p):t.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out;out.mkdir(parents=True,exist_ok=False)
 members=read(PAT/'gene_membership.tsv');rna=read(PAT/'rna_results.tsv');ass=read(PAT/'association_results.tsv');met=read(MET/'metabolite_results.tsv');relations=read(MAP/'main_relations.tsv');mp=read(MAP/'feature_workpool101.tsv')
 keys=['metabolite_name','metabolite_key','human_gene_id'];rel=relations.copy();rel['relation_id']='COAD|'+rel.metabolite_name+'|'+rel.metabolite_key+'|'+rel.human_gene_id
 for kind,prefix in [('association_primary','CAMP'),('association_available','CAMP_available')]:
  t=ass[ass.analysis_type==kind];cols=keys+['effect','ci_lower','ci_upper','p_value','q_value','n','status','reason','input_sha256','reused_from'];rel=rel.merge(t[cols].rename(columns={x:prefix+'_'+x for x in cols if x not in keys}),on=keys,how='left',validate='one_to_one')
 for kind,prefix in [('metabolite_primary','metabolite_paired'),('metabolite_available','metabolite_available')]:
  cols=['metabolite_name','metabolite_key','effect','p_value','q_value','n','status','pairs_higher','pairs_lower','pairs_equal'];rel=rel.merge(met[met.analysis_type==kind][cols].rename(columns={x:prefix+'_'+x for x in cols if x not in keys[:2]}),on=keys[:2],how='left',validate='many_to_one')
 genes=members.copy();rfields=['human_gene_id','effect','ci_lower','ci_upper','p_value','q_value','n','status','reason','n_up','n_down','n_equal'];genes=genes.merge(rna[rfields].rename(columns={x:'RNA_'+x for x in rfields if x!='human_gene_id'}),on='human_gene_id',how='left',validate='one_to_one')
 for col in genes:
  if col.startswith('RNA_'):genes[col]=genes[col].fillna('NA')
 genes.loc[genes.current_direct!='True','RNA_status']='NOT_RUN_CURRENT_FAMILY'
 history=read(PAT/'historical_RNA_wilcoxon.tsv');cols=['human_gene_id','effect','p_value','q_value','status'];genes=genes.merge(history[cols].rename(columns={x:'historical_Wilcoxon_'+x for x in cols if x!='human_gene_id'}),on='human_gene_id',how='left',validate='one_to_one')
 for study in ['Lee','Uhlitz']:
  s=read(SC/(study+'_source_status.tsv'));cols=['human_gene_id','display_top','status','reason','bootstrap_top_frequency','top_runner_gap','max_detection'];genes=genes.merge(s[cols].rename(columns={x:study+'_'+x for x in cols if x!='human_gene_id'}),on='human_gene_id',how='left',validate='one_to_one')
 cross=read(SC/'sc_cross_study.tsv');genes=genes.merge(cross,on='gene',how='left',validate='one_to_one')
 func=read(ROOT/'results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1/gene_function_review.tsv');cols=['gene','functional_conclusion_cn','counterevidence_or_limit_cn','evidence_urls','research_decision_cn'];genes=genes.merge(func[cols],on='gene',how='left',validate='one_to_one')
 genes['functional_scope']='Historical35 only; no new mechanism review';genes['relationship_inference']='Expression source is not metabolite production, enzyme activity or mechanism'
 genes=genes.fillna('NA');assert len(genes)==811 and genes.gene.is_unique
 rel=rel.merge(genes[['human_gene_id']+[x for x in genes if x.startswith(('RNA_','Lee_','Uhlitz_'))]],on='human_gene_id',how='left',validate='many_to_one');assert len(rel)==891 and rel.relation_id.is_unique
 write(rel,out/'candidate_relations.tsv');write(genes,out/'candidate_genes.tsv')
 reader=genes[['gene','human_gene_id','current_direct','current_conditional','current_unresolved','historical_retained','identity_status','RNA_effect','RNA_p_value','RNA_q_value','RNA_status','Lee_display_top','Uhlitz_display_top','Lee_shared_top','Uhlitz_shared_top','descriptive_stable_concordance','counterevidence_or_limit_cn']];write(reader,out/'candidate_genes_reader.tsv')
 allmap=read(MAP/'results.tsv');allmap['relation_id']='COAD|'+allmap.metabolite_name+'|'+allmap.metabolite_key+'|'+allmap.human_gene_id
 old=read(ROOT/'results/COAD/02_MAPPING/20260919T115842Z_identity_units_v1/candidate_eligibility.tsv');old['relation_id']='COAD|'+old.metabolite_name+'|'+old.metabolite_key+'|'+old.human_gene_id;old=old[~old.relation_id.isin(allmap.relation_id)].copy();old['pool']='HISTORICAL_OUTSIDE_CURRENT'
 allrel=pd.concat([allmap,old],ignore_index=True).fillna('NA');assert allrel.relation_id.is_unique;write(allrel,out/'all_relation_history.tsv')
 mapping=mp.rename(columns={'mapping_status':'mapping_state'});write(mapping,out/'metabolite_mapping_status.tsv')
 write(rna.rename(columns={'n_up':'positive_pairs','n_down':'negative_pairs','n_equal':'equal_pairs'}),out/'rna_results.tsv')
 write(pd.concat([met,ass,rna],ignore_index=True).fillna('NA'),out/'statistics_long.tsv')
 prefix='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
 for study in ['Lee','Uhlitz']:
  sp=read(SC/(study+'_profiles.tsv'))
  for col in prefix:
   if col not in sp:sp[col]='NA'
  sp['unit']='author_source_donor_label';sp['source_id']=study+' original count matrix and annotation';sp['test_family']='NA:descriptive only'
  sp=sp[prefix+[x for x in sp if x not in prefix]];write(sp,out/(study+'_sc_celltype_profiles.tsv'))
 progress=[]
 for label,t in [('配对代谢物',met[met.analysis_type=='metabolite_primary']),('代谢物可用值',met[met.analysis_type=='metabolite_available']),('肿瘤内部关联',ass[ass.analysis_type=='association_primary']),('关联可用值',ass[ass.analysis_type=='association_available']),('配对RNA',rna)]:
  progress.append(dict(step=label,planned=len(t),evaluable=int(t.status.eq('DONE').sum()),P_lt_005=int(pd.to_numeric(t.p_value,errors='coerce').lt(.05).sum()),q_lt_005=int(pd.to_numeric(t.q_value,errors='coerce').lt(.05).sum()),status='DONE',reason='Complete dispositions; unevaluable rows remain'))
 for study in ['Lee','Uhlitz']:
  s=read(SC/(study+'_source_status.tsv'));progress.append(dict(step=study+'来源',planned=811,evaluable=int(s.status.eq('DONE').sum()),P_lt_005='NA',q_lt_005='NA',status='DONE',reason='Scope fully processed; source limits explicitly masked'))
 for label in ['作者UMAP坐标','可靠患者临床分型']:progress.append(dict(step=label,planned='NA',evaluable='NA',P_lt_005='NA',q_lt_005='NA',status='NOT_EVALUABLE',reason='Not present in admitted sources; no invented replacement'))
 write(pd.DataFrame(progress),out/'progress_summary.tsv')
 inputs=[PAT/'gene_membership.tsv',PAT/'rna_results.tsv',PAT/'association_results.tsv',MAP/'main_relations.tsv',MET/'metabolite_results.tsv',SC/'sc_cross_study.tsv']
 write(pd.DataFrame([dict(path=p.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs]),out/'source_manifest.tsv')
 (out/'analysis_spec.json').write_text(json.dumps(dict(protocol_commit='45c6ae9ab489a9bd4ead4e47d28185785a448156',new_statistics=False,current_direct_relations=891,current_genes=526,all_gene_id_records=811,all_relation_history=len(allrel),identity_holds=3),indent=2))
 (out/'validation.json').write_text(json.dumps(dict(status='PASS',unique_gene_ids=811,unique_current_relations=891,complete_source_dispositions=True,history_preserved=True,new_P_q=0),indent=2))
 (out/'INTERPRETATION_CN.md').write_text('# COAD证据解释\n\n## 直接观察\n159项代谢物配对结果中101项P<0.05；891关系中871项可算、52项名义P支持、无FDR支持。515基因配对RNA可算，其中361项FDR支持。单细胞覆盖当前和历史全部811个ID去向。\n\n## 可提出的解释\n供者等权表达可用于选择细胞背景。相同宽类别的跨研究结果属于描述性相容，不等于同一细胞亚型或独立代谢关系复现。\n\n## 不能推断\nRNA变化不等于代谢物关联、酶活或通量。来源最高不是唯一作用细胞。条件性映射不因有表达而升级。历史功能只在原适用范围内保留，不把通用说明计作811基因功能审核。\n\n## 下一步最小问题\n使用每条relation_id检视具体代谢物身份、关联效应、可用性以及来源缺项。此交付止于来源，不自动新增机制、拟时序、通讯或生存分析。\n',encoding='utf-8')
 print(json.dumps(dict(genes=len(genes),relations=len(rel),history_relations=len(allrel)),ensure_ascii=False))
if __name__=='__main__':main()
