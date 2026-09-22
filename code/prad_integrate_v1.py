"""Integrate only released PRAD aggregate evidence; preserve every candidate and old results."""
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='20260922T142000Z_discovery_v1';BASE=ROOT/'results/PRAD';OUT=BASE/'07_INTEGRATION'/RUN
def read(stage,name):return pd.read_csv(BASE/stage/RUN/name,sep='\t')
def save(d,name):d.to_csv(OUT/name,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 met=read('04_ROBUSTNESS','metabolite_paired.tsv');edges=read('02_MAPPING','direct_relations.tsv');assoc=read('03_PATIENT','tumor_association_all_families.tsv');rna=read('03_PATIENT','paired_RNA_all_families.tsv');sc=read('06_EXTERNAL','sc_source_stability.tsv');coverage=read('02_MAPPING','metabolite_mapping_status.tsv')
 old=read('01_CAMP','historical_cohort_effects_preserved.tsv');old['metabolite_name']=old.feature_name.str.strip()
 allmet=met.merge(old[['metabolite_name','hedges_g','wilcoxon_p','wilcoxon_fdr','analysis_design']].rename(columns={c:'historical_'+c for c in ['hedges_g','wilcoxon_p','wilcoxon_fdr','analysis_design']}),on='metabolite_name',validate='one_to_one')
 for file,label in [('metabolite_available_sensitivity.tsv','available'),('metabolite_CAPT_concordant_sensitivity.tsv','CAPT36')]:
  z=read('04_ROBUSTNESS',file)[['feature_id','n','effect','p_value','q_value','status']];allmet=allmet.merge(z.rename(columns={c:label+'_'+c for c in z if c!='feature_id'}),on='feature_id',validate='one_to_one')
 save(allmet,'metabolite_all361_integrated.tsv');save(allmet[allmet.q_value.lt(.05)],'metabolite_primary_FDR9.tsv')
 relations=edges.copy()
 z=met[['feature_id','n','effect','p_value','q_value','n_up','n_down','up_fraction','down_fraction']]
 relations=relations.merge(z.rename(columns={c:'metab_'+c for c in z if c!='feature_id'}),on='feature_id',validate='many_to_one')
 for family,label in [('REL_TUMOR_PRIMARY','patient_primary'),('REL_TUMOR_AVAILABLE','patient_available'),('REL_TUMOR_BATCH','patient_batch')]:
  z=assoc[assoc.test_family.eq(family)][['relation_id','n','effect','ci_lower','ci_upper','p_value','q_value','status','reason']]
  relations=relations.merge(z.rename(columns={c:label+'_'+c for c in z if c!='relation_id'}),on='relation_id',validate='one_to_one')
 for family,label in [('RNA_PAIRED_CURRENT','RNA_primary'),('RNA_PAIRED_CAPT','RNA_CAPT36')]:
  z=rna[rna.test_family.eq(family)][['stable_gene_id','n','effect','p_value','q_value','status']]
  relations=relations.merge(z.rename(columns={c:label+'_'+c for c in z if c!='stable_gene_id'}),on='stable_gene_id',validate='many_to_one')
 z=sc[sc.partition.eq('cancer')][['stable_gene_id','top_celltype','bootstrap_top_frequency','status','reason']]
 relations=relations.merge(z.rename(columns={c:'SC_'+c for c in z if c!='stable_gene_id'}),on='stable_gene_id',validate='many_to_one')
 relations['same_cohort_postselection']=True;relations['independent_metabolite_RNA_validation']='NOT_RUN';relations['functional_validation']='NOT_RUN'
 assert len(relations)==len(edges) and relations.relation_id.is_unique
 save(relations,'candidate_relations_integrated.tsv')
 gene_rows=[]
 for (gene,gid),z in relations.groupby(['gene','stable_gene_id']):
  r=z.iloc[0];rnasig=bool(r.RNA_primary_q_value<.05)
  gene_rows.append(dict(cancer='PRAD',gene=gene,stable_gene_id=gid,current_pool=True,history_only=False,n_direct_relations=len(z),
   relation_ids=';'.join(sorted(z.relation_id)),metabolites=';'.join(sorted(set(z.metabolite_name))),
   n_metabolites_primary_FDR=int(z.loc[z.metab_q_value.lt(.05),'feature_id'].nunique()),
   n_patient_nominal_relations=int(z.patient_primary_p_value.lt(.05).sum()),n_patient_FDR_relations=int(z.patient_primary_q_value.lt(.05).sum()),
   RNA_n=r.RNA_primary_n,RNA_effect=r.RNA_primary_effect,RNA_p=r.RNA_primary_p_value,RNA_q=r.RNA_primary_q_value,RNA_status=r.RNA_primary_status,
   cell_source=r.SC_top_celltype,cell_source_status=r.SC_status,cell_source_bootstrap_frequency=r.SC_bootstrap_top_frequency,
   RNA_background='FDR_SUPPORTED' if rnasig else 'NOMINAL_EXPLORATORY' if r.RNA_primary_p_value<.05 else 'NOT_SUPPORTED_THIS_TEST',
   patient_association_evidence='NO_RELATION_FDR_SUPPORTED_CURRENT_FAMILY',cell_background='ONE_STUDY_DESCRIPTION' if r.SC_status=='DONE' else 'NOT_EVALUABLE',
   data_limitations='same_CAMP_postselection;CAPT_discordance;one_sc_study;no_causal_or_functional_validation',
   next_action='保留RNA差异支持线索，优先补独立关系证据' if rnasig else '保留全候选背景，暂不宣称代谢调控作用'))
 genes=pd.DataFrame(gene_rows)
 hist=pd.read_csv(OUT/'historical_locked_relations.tsv',sep='\t');historical_genes=set(hist.gene)
 assert historical_genes.issubset(set(genes.gene)), 'Preserve history-only genes explicitly before finalizing'
 genes['present_in_historical_PRAD_panel']=genes.gene.isin(historical_genes)
 assert genes.stable_gene_id.is_unique;save(genes,'candidate_genes_integrated.tsv');save(coverage,'workpool_all69_mapping_coverage.tsv')
 view=genes[['gene','metabolites','n_direct_relations','n_metabolites_primary_FDR','n_patient_nominal_relations','n_patient_FDR_relations','RNA_effect','RNA_p','RNA_q','cell_source','cell_source_bootstrap_frequency','cell_source_status','next_action']].copy()
 view.columns=['基因','对应代谢物','直接关系数','代谢物FDR支持项数','肿瘤关联名义显著数','肿瘤关联FDR支持数','配对RNA均值差','RNA_P','RNA_q','肿瘤组织表达来源','来源排名保持率','来源可评估状态','当前解释与下一步']
 save(view,'全部101基因_易读比较表.tsv')
 summary=dict(status='PARTIAL',completed='paired discovery;bounded direct mapping;all candidate patient statistics;one tumor-study cell source;complete current candidate tables',
  retained_metabolites=len(met),metabolite_workpool=int(met.p_value.lt(.05).sum()),metabolite_FDR=int(met.q_value.lt(.05).sum()),
  direct_features=int(edges.feature_id.nunique()),direct_relations=len(edges),candidate_genes=len(genes),RNA_FDR=int(genes.RNA_q.lt(.05).sum()),
  patient_relation_FDR=int(relations.patient_primary_q_value.lt(.05).sum()),single_cell_source_evaluable=int(genes.cell_source_status.eq('DONE').sum()),
  historical_genes_preserved=len(historical_genes),mapping_unresolved=int(coverage.n_direct_relations.eq(0).sum()),second_sc_study_complete=False,functional_validation_complete=False)
 (OUT/'validation.json').write_text(json.dumps(dict(summary,all_join_cardinalities_checked=True,all_current_and_historical_genes_preserved=True,numerical_statistics_recomputed=False),indent=2),encoding='utf-8')
 (OUT/'analysis_spec.json').write_text(json.dumps(dict(version='prad_integration_v1',run_id=RUN,key='feature_id x stable_gene_id x version',statistics='existing released values only',ranking='no opaque score,no best-relation substitution',missing='NA and explicit states',source_boundaries='same-cohort association;RNA background;one-study descriptive SC;no independent mechanism claim'),indent=2),encoding='utf-8')
 paths=[]
 for stage in ['01_CAMP','02_MAPPING','03_PATIENT','04_ROBUSTNESS','06_EXTERNAL']:
  for p in (BASE/stage/RUN).glob('*.tsv'):paths.append(p)
 save(pd.DataFrame([dict(source_path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in paths]),'source_manifest.tsv')
 readme=f'''# PRAD参考BRCA主线：第一轮整合交付\n\n问题：PRAD代谢差异能否连接到直接基因、同肿瘤RNA关联与细胞表达背景？\n\n输入与范围：91个作者肿瘤case，43对肿瘤/正常；全361个作者保留特征。原351项冻结结果与旧5条关系/10条统计原样保留。\n\n实际结果：\n\n|环节|结果|\n|---|---|\n|代谢物主配对|361项可评估；69项P<0.05；9项q<0.05|\n|作者可用值敏感性|344项可评估；66项P<0.05；9项q<0.05|\n|CAPT一致36对敏感性|361项可评估；57项P<0.05；2项q<0.05|\n|直接映射|48项代谢物；161条关系；101个基因；21项待身份/证据补充|\n|肿瘤内关联|161条均可评估；18条名义P<0.05，0条q<0.05|\n|配对RNA|101基因可评估；35个P<0.05；21个q<0.05|\n|单细胞来源|一项24供者肿瘤研究；{summary['single_cell_source_evaluable']}个基因可评价主要来源；第二研究未完成|\n\n新手解释：PRAD有代谢和RNA差异，但本轮没有代谢物—RNA关系通过其完整检验族的FDR校正。因此不能把RNA差异或单细胞高表达补成已验证的调控轴。全部候选均保留供比较，而不是只交显著项。\n\n限制与反证：7对CAPT字段不一致；361个特征中10个化学键碰撞已分峰保留；作者可用值掩码不等于确认检出；关联属同队列事后选择；单细胞只有一研究，癌旁是同研究参考；CELLxGENE网页GEO外链冲突在单细胞阶段记录；没有功能验证。21项未充分映射不是阴性。\n\n当前决定：第一轮作为可讨论的探索性证据包。21个RNA FDR基因提供表达背景，不自动升级治疗靶点。未完成项：21特征映射补证、第二独立单细胞研究、独立患者关系验证；功能深入不是本轮来源主线强制项。\n\n下一步：优先补可改变判断的独立关系证据与第二研究可靠注释；不改变当前阈值、不删除未显著候选、不追加机制分析。\n\n文件入口：\n\n- [全部101基因易读比较表](全部101基因_易读比较表.tsv)\n- [完整161关系表](candidate_relations_integrated.tsv)\n- [完整101基因表](candidate_genes_integrated.tsv)\n- [全361代谢物及新旧统计](metabolite_all361_integrated.tsv)\n- [9项代谢物主分析FDR支持](metabolite_primary_FDR9.tsv)\n- [69项完整映射覆盖](workpool_all69_mapping_coverage.tsv)\n- [单细胞来源、热图与点图](../../06_EXTERNAL/{RUN}/README_CN.md)\n\n复现：python code/prad_integrate_v1.py；只读取已发布汇总，不读取患者矩阵或重算P/q。完整阶段代码与参数见coordination/stages/PRAD.tsv。\n'''
 (OUT/'README_CN.md').write_text(readme,encoding='utf-8')
 save(pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='checksums.tsv']),'checksums.tsv')
 print(json.dumps(summary,ensure_ascii=False))
if __name__=='__main__':main()
