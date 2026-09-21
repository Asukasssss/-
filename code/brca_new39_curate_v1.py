"""Explicit human-readable evidence decisions; no automated functional hit calling."""
from pathlib import Path
import sys,json,hashlib,pandas as pd
R=Path(sys.argv[1]).resolve();repo=Path(__file__).resolve().parents[1]
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
# gene, PMID, evidence class, model, intervention, observed outcome, limit, depth
E=[
('ABHD12','32366405','BRCA_GENETIC','MCF7;MDA-MB-231','siRNA','降低后增殖、迁移、侵袭受抑','摘要层；不能证明当前脂质介导','abstract'),
('ACSL1','31581558','BRCA_GENETIC','MDA-MB-231;TNFalpha条件','siRNA;triacsin C','降低ACSL1后TNF诱导GM-CSF减少','直接终点是炎症因子；triacsin C并非单基因特异证明','abstract_and_methods_excerpt'),
('ACSL3','38953696','BRCA_GENETIC_DIRECTION_CAUTION','MCF7;MDA-MB-231','siRNA;过表达','降低后增殖、迁移、EMT增强','不能直接建议抑制；原研究新数据需申请','abstract_and_methods_excerpt'),
('ACSL4','39700137','BRCA_GENETIC','TNBC;MDA-MB-231 RNA实验','遗传降低;PRGL493','迁移、侵袭及转移相关终点降低','脂肪酸氧化与表观调控背景；不是全部BRCA通用结论','abstract_and_data_availability'),
('ACSL4','38078907','BRCA_GENETIC','乳腺癌细胞及肺转移模型','敲低;过表达救援','降低后脂滴、迁移及转移结节减少','与ZEB2共同背景；不是独立患者复现','abstract'),
('ACSL5','42755793','BRCA_GENETIC_DIRECTION_CAUTION','MDA-MB-231;BT-549;COL I','过表达','部分逆转胶原I诱导的恶性表型','体外COL I环境；不推广至肺转移营养背景','abstract'),
('ACSL5','41570334','BRCA_FUNCTION_CONTEXT','肺偏好转移乳腺癌;富棕榈酸环境','靶向ACSL5/COX2/EP4轴','支持肺转移存活与药物增敏问题','摘要未逐项区分基因构建；与COL I模型方向不可简单合并','abstract'),
('ACSL6','41271644','FUNCTION_MODEL_DETAIL_PENDING','结直肠癌与乳腺癌;DHA条件','降低;过表达','ACSL6参与DHA敏感性与铁死亡','需分开乳腺与结直肠各实验；不能把肠癌小鼠疗效转到乳腺癌','abstract_and_selected_figure_text'),
('ALOX5','33224971','BRCA_FUNCTION_CONTEXT','HER2背景;SKBR3;BT-474','ALOX5功能抑制相关实验','生长、迁移与HER2/ALOX5活性相关','本轮不把HER2 siRNA当作ALOX5自身遗传验证；干预细节待核对','abstract_and_selected_methods'),
('CYP4F2','38646426','BRCA_GENETIC','ER阳性;MCF7;雌激素背景','过表达;shRNA','促进增殖与抗凋亡;研究涉及20-HETE','CYP4F11并行；不把其结果全归CYP4F2','abstract_and_methods_excerpt'),
('ENPP2','37296899','BRCA_PHARMACOLOGY_WITH_GENETIC_NULL','E0771;MMTV-PyMT;脂肪细胞特异敲除','脂肪细胞敲除;IOA-289','脂肪细胞敲除未降低生长或肺转移；全身抑制在E0771降低生长','细胞来源重要；负结果必须保留；GSE176078是背景不是自身干预','abstract_and_data_availability'),
('LPCAT1','36909375','BRCA_GENETIC','MDA-MB-231;紫杉醇耐药株','shRNA','降低后增殖、迁移、侵袭和紫杉醇耐受降低','FOXA1救援不等于证明CAMP特定LPC介导','abstract_and_methods_excerpt'),
('LYPLA1','41992085','BRCA_PHARMACOLOGY','MDA-MB-468 TNBC','ML348抑制','增殖下降;以细胞周期停滞为主','不是本项目新DepMap分析；主要药理证据，未证明当前LPC关联介导','abstract'),
('LYPLA2','28065656','EPITHELIAL_MODEL_PENDING','Snail转化的上皮细胞;具体来源待核对','APT2抑制与敲低','Scrib定位与棕榈酰化恢复;MEK降低','不能仅凭摘要称为已验证乳腺癌模型；全文取得失败','abstract'),
('MGLL','35351947','BRCA_PHARMACOLOGY','TNBC-HBMEC共培养;4T1-BrM5','AM9928','减少脑内皮穿越、相关炎症与脑定植','药理为主；不能由CAMP整体RNA下降判定抑制方向','abstract_and_selected_methods'),
('PTGS2','40370195','BRCA_GENETIC','DANCR/脑转移背景乳腺癌细胞','PTGS2沉默','逆转DANCR诱导的部分体外促肿瘤效应','公开GEO多为背景队列，不能当PTGS2直接干预数据','abstract_and_data_availability'),
('SDHA','40821099','BRCA_GENETIC','MCF7;SKBR3','敲低','增殖、克隆形成与迁移下降','全文取得失败；不据此把所有SDH亚基统一抑制','abstract'),
('SDHC','31164982','BRCA_GENETIC_DIRECTION_CAUTION','MCF7','CRISPR杂合降低','出现EMT表型及细胞形态改变','球体生长/稳定性也下降；不能只用一种终点判断促癌抑癌','abstract_and_methods_results_excerpt'),
('SDHD','37974462','BRCA_GENETIC_DIRECTION_CAUTION','MCF7','siRNA','降低后增殖、迁移、侵袭增强并伴EMT标记变化','不能把SDHA的干预方向直接移到SDHD','abstract_and_results_excerpt'),
('SLC13A3','39515327','FUNCTION_MODEL_DETAIL_PENDING','肿瘤-巨噬细胞;模型谱待核对','基因缺失;候选抑制剂','衣康酸摄取与铁死亡抵抗/免疫治疗响应相关','摘要未确认具体乳腺模型；此衣康酸关系不等于CAMP当前映射底物','abstract'),
('SLC25A10','30205167','FUNCTION_MODEL_DETAIL_PENDING','慢性循环低氧;肿瘤模型待核对','遗传或药理抑制','增强放射杀伤与抗氧化背景有关','乳腺患者生存证据不等于乳腺模型直接干预；全文未取得','abstract'),
('SLC47A1','37583829','BRCA_DRUG_TRANSPORT_CONTEXT','MDA-MB-436及其他癌细胞','MATE1抑制;铂-吖啶化合物','转运抑制降低药物活性','涉及药物摄取敏感性，非当前肌酐异常的功能解释','abstract'),
('SLC47A2','27959931','ASSOCIATION_NOT_SELF_PERTURBATION','19癌细胞系;非乳腺专属','二甲双胍处理','MATE2表达与药物反应相关','未证明SLC47A2自身干预；不能借用为直接遗传功能支持','abstract'),
('SUCLA2','37253003','BRCA_GENETIC','MDA-MB231;MDA-MB468;脱附条件','稳定敲低;过表达;构建救援','降低后失巢凋亡增强，贴壁增殖/凋亡未显著改变','乳腺细胞与肺癌小鼠实验分开；琥珀酸未变且补充不救援，提示非经典功能','abstract_and_results_figure_excerpt'),
('SUCLG1','37253003','BRCA_GENETIC_NULL_SPECIFIC_ENDPOINT','A549及MDA-MB231;脱附','敲低','未像SUCLA2降低那样诱导失巢凋亡','只限定该研究终点；不是无任何生物学作用','results_figure3_excerpt'),
('SUCLG2','37253003','BRCA_GENETIC_NULL_SPECIFIC_ENDPOINT','A549及MDA-MB231;脱附','敲低','未像SUCLA2降低那样诱导失巢凋亡','不能从同一复合体成员推断相同功能','results_figure3_excerpt')]
d=pd.concat([pd.read_csv(R/f,sep='\t',dtype=str) for f in ['literature_hits.tsv','literature_hits_focused.tsv']]).drop_duplicates(['gene','pmid','doi']);records=[]
for i,(g,pmid,cls,model,pert,result,limit,depth) in enumerate(E):
 hit=d[d.pmid.eq(pmid)];assert len(hit)>0,pmid
 records.append(dict(evidence_id=f'N39E{i+1:03}',gene=g,pmid=pmid,source_url='https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/',title=hit.iloc[0].title,evidence_class=cls,applicable_model=model,intervention=pert,outcome_cn=result,counterevidence_and_limits_cn=limit,reading_depth=depth,own_patient_statistics=False,resource_status='NOT_YET_ELIGIBLE',full_figure_audit=False))
ev=pd.DataFrame(records);save(ev,R/'functional_evidence.tsv')
genes=pd.read_csv(repo/'results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/genes_unique.tsv',sep='\t');genes=genes.loc[~genes.in_original117,'gene'].tolist()
logs=pd.concat([pd.read_csv(R/f,sep='\t') for f in ['search_log.tsv','search_log_focused.tsv']]);summ=[]
for g in genes:
 es=ev[ev.gene.eq(g)];h=d[d.gene.eq(g)];ll=logs[logs.gene.eq(g)]
 status='DONE' if len(es) else 'NEEDS_REVIEW';classes=';'.join(sorted(set(es.evidence_class))) if len(es) else 'NO_CURATED_BREAST_SELF_PERTURBATION_THIS_PASS'
 reason='本轮已摘录指定来源，未穷尽全部文献' if len(es) else '定向检索已记录；尚未提取可核对的直接乳腺自身干预证据，不是生物学阴性'
 if g=='LPCAT4':reason='Q643R3与旧LPCAT4/MBOAT2命名歧义；按稳定ID核对前不合并文献'
 summ.append(dict(gene=g,searches=len(ll),searches_done=int(ll.status.eq('DONE').sum()),returned_hits=len(h),curated_evidence_records=len(es),evidence_ids=';'.join(es.evidence_id),functional_evidence_class=classes,status=status,reason_cn=reason,model_and_limits_cn=' | '.join(es.applicable_model+'：'+es.counterevidence_and_limits_cn),next_minimum_work_cn='先核对具体模型及可用公开干预数据；不直接开始机制分析' if len(es) else '保留候选；需要深入时扩展别名或全文检索',depmap_status='NOT_RUN',singlecell_extension_status='NOT_RUN',public_intervention_matrix_status='SOURCE_IDENTIFIED_METADATA_CHECKED' if g=='ACSL4' else 'NOT_ESTABLISHED_THIS_PASS'))
su=pd.DataFrame(summ);assert len(su)==39 and su.gene.is_unique;save(su,R/'new39_functional_summary.tsv')
save(d.drop(columns=['abstract']),R/'search_hit_index.tsv')
resources=pd.DataFrame([
 dict(gene='ACSL4',source_id='GSE283282',source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283282',status='NEEDS_REVIEW',actual_scope='MDA-MB-231;control vs shACSL4 vs PRGL493;author reports three batches per group;9 samples',matrix='GSE283282_Whole-Genome_Genes_Read_Count_All_Samples.txt.gz',next_step='processed counts on server only;verify GSM groups,units,independent batches before analysis',why_not_done='metadata read;matrix not downloaded or modeled'),
 dict(gene='ACSL3',source_id='PMID38953696',source_url='https://pubmed.ncbi.nlm.nih.gov/38953696/',status='NOT_EVALUABLE',actual_scope='new data available upon request;GSE22226/GSE123833 are background',matrix='NA',next_step='do not treat background series as ACSL3 perturbation',why_not_done='no public self-perturbation matrix established'),
 dict(gene='ENPP2',source_id='PMID37296899',source_url='https://pubmed.ncbi.nlm.nih.gov/37296899/',status='NOT_EVALUABLE',actual_scope='GSE176078 is existing patient singlecell background',matrix='NA',next_step='use functional/model context;no repeat singlecell run',why_not_done='not ENPP2 perturbation dataset'),
 dict(gene='PTGS2',source_id='PMID40370195',source_url='https://pubmed.ncbi.nlm.nih.gov/40370195/',status='NOT_EVALUABLE',actual_scope='listed GEO series are metastasis/background comparisons',matrix='NA',next_step='verify dedicated PTGS2 perturbation matrix before reuse',why_not_done='no eligible own-gene perturbation matrix established'),
 dict(gene='SUCLA2',source_id='PMID37253003',source_url='https://pubmed.ncbi.nlm.nih.gov/37253003/',status='NEEDS_REVIEW',actual_scope='breast detachment endpoints;metabolomics described in lung cells',matrix='NA',next_step='separate breast endpoints from lung metabolomics;do not transfer model',why_not_done='public matrix eligibility not verified')]);save(resources,R/'intervention_resource_inventory.tsv')
# Only concise evidence and bibliographic metadata are published; cached abstracts/full texts stay ignored.
(R/'.gitignore').write_text('sources/\nfulltext/\nliterature_hits*.tsv\npublic_resource_lookup.tsv\n')
(R/'.gitattributes').write_text('* -text\n')
v=dict(genes=39,queries=len(logs),queries_done=int(logs.status.eq('DONE').sum()),unique_returned_sources=int(d[['pmid','doi']].drop_duplicates().shape[0]),curated_records=len(ev),genes_with_curated_records=ev.gene.nunique(),fulltext_audit_all=False,new_patient_tests=0,resource_matrices_analyzed=0,curated_source_classes=ev.evidence_class.value_counts().to_dict());(R/'validation.json').write_text(json.dumps(v,indent=2))
spec=json.loads((R/'analysis_spec.json').read_text());spec.update(query_families=['symbol-title_abstract plus breast plus intervention','gene-title/selected_aliases plus breast mentions;reviews excluded;off-context hits expected'],manual_review='26 selected evidence records;not all returned papers full-text audited',resource_lookup='20 selected sources;metadata/specified paragraphs only;no intervention matrix calculation');(R/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
paths=list((R/'sources').glob('*.json'))+list((R/'sources').glob('*.txt'))+list((R/'fulltext').glob('*.xml'))+[repo/'code/brca_new39_literature_v1.py',repo/'code/brca_new39_curate_v1.py',repo/'code/brca_new39_resource_lookup_v1.py']
save(pd.DataFrame([dict(file=str(p.relative_to(repo)).replace('\\','/'),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),published='code' in p.parts) for p in paths]),R/'source_manifest.tsv')
text=f'''# 新39基因功能证据第一批

## 本轮问题
新候选是否有适用干预证据、相反结果和可以复用的数据？不以CAMP关联显著作为检索门槛。

## 输入与范围
全部39新基因；两类定向查询，共{len(logs)}条，成功{int(logs.status.eq('DONE').sum())}条；每式最多20条。选取{len(ev)}条证据，覆盖{ev.gene.nunique()}个基因，其他候选仍保留。它是有边界的第一轮取证，不是系统综述或全量全文审计。

## 实际结果
ABHD12、ACSL1、ACSL3、ACSL4、CYP4F2、LPCAT1、SDHA、SDHC、SDHD、SUCLA2等有具体干预信息。LYPLA1、MGLL有相关药理研究；药理不冒充遗传敲除。
ACSL3/SDHC/SDHD降低后某些恶性终点增强，不能统一建议抑制。ENPP2脂肪细胞特异敲除的无效结果保留。SUCLA2在脱附条件下有作用，贴壁生长终点不相同；SUCLG1/SUCLG2在该研究的失巢凋亡终点未呈现SUCLA2式效应。
找到ACSL4直接干预转录组GSE283282，官方说明是MDA-MB-231对照、shACSL4、药物PRGL493三组，各三个作者批次；有处理后计数文件。尚未下载矩阵或核实独立性并计算。

## 新手解释
有功能研究代表前人改变过该基因或药物靶向相关蛋白并观察终点。它既不意味着我们新发现，也不意味着研究中的底物就是CAMP中关联的代谢物。

## 限制与反证
未提取证据的基因标NEEDS_REVIEW，不能写成无作用。部分证据只核到摘要、指定正文或图注，阅读深度逐行记录。SLC13A3/SLC25A10等的具体乳腺模型仍待核对；不能把乳腺患者关联与其他癌模型拼成乳腺直接功能证明。LPCAT4身份仍需核对。没有新增DepMap或单细胞分析。

## 当前决定
先把多层证据回接156基因总表。保留方向矛盾和模型限制；不会因为某个q较小就宣布第一靶点。

## 下一步
需要选择下一项真正的新计算时，ACSL4已有明确直接干预数据入口；先做样本设计与计数单位的局部核对。其余候选按适用模型和资源明确度安排，不为凑满资料重新跑全基因流程。

## 复现
code/brca_new39_literature_v1.py <result_dir> 及 --focused；随后 brca_new39_resource_lookup_v1.py、brca_new39_curate_v1.py。后者固定人工取证结论，非关键词自动认定功能。摘要和全文缓存不随Git发布；公开元数据索引、查询、简要判断与来源哈希可追溯。
'''
(R/'README_CN.md').write_text(text,encoding='utf-8');print(json.dumps(v))
