"""Rebuild scoped mechanism review from curated primary-source readings;no new statistics."""
import csv,json,hashlib,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[2];RUN='20260921T031711Z_mediation4_v1';OUT=R/'results/COAD/05_FUNCTION'/RUN
SOURCE=R/'results/COAD/06_EXTERNAL/20260920T091100Z_cell_paired35_v1/relations39_integrated.tsv'
def read(p):
 with p.open(encoding='utf-8-sig',newline='')as f:return list(csv.DictReader(f,delimiter='\t'))
def save(name,rows):
 with (OUT/name).open('w',encoding='utf-8',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def write(name,text):(OUT/name).write_text(text,encoding='utf-8')
def main():
 OUT.mkdir(parents=True,exist_ok=True);old=read(SOURCE)
 fields=['review_id','prior_id','gene','year','url','new_research_to_project','reading_depth','model','gene_intervention','exact_metabolite_readout','rescue_or_mediation','limit','next_question']
 records=[
 ['M01','E01','SLC6A6',2014,'https://www.nature.com/articles/srep04852',False,'Publisher indexed Results/Fig2/5 and uptake description','DLD1/HT29/HCT15;Colo320DM overexpression','KD/OE','Taurine uptake changes with gene perturbation','Survival/drug phenotypes reported;exact metabolite rescue not established by this reading','Uptake is not whole-tissue taurine or mitochondrial pool;no clinical causality','Distinguish membrane uptake,subcellular pool and survival under treatment'],
 ['M02','E07','NNMT',2016,'https://www.oncotarget.com/article/9962/text/',False,'Publisher Results/Fig6-8 focused full-text sections','HT29 KD;SW480 OE;5-FU context','KD/OE','Intracellular1-MNA HPLC-UV','HT29 NNMT-KD+5-FU:2/4mM1-MNA partly rescues ROS/apoptosis and5-FU IC50','Pharmacologic addition and treatment context;patient concentration relevance unverified;not SAM rescue','Confirm relevant local1-MNA dose and cell compartment;separate SAM hypothesis'],
 ['M03','E08','NNMT',2021,'https://www.jcancer.org/v12p6170.htm',False,'Publisher methods/results/Fig4-6 focused sections','CCD-18Co colon fibroblast line;SW480/HCT116 cancer cells','NNMT downregulation in fibroblast model;CRC KD/OE','LC-MS/MS intracellular and culture-medium1-MNA','1/2mM exogenous1-MNA promotes migration/invasion;not equivalent to complete conditioned-medium mediation/rescue','CCD-18Co is a colon fibroblast line,not proof of patient-derived CAF-specific mechanism;SAM not tested as mediator here','Test cell-source transfer and physiological-dose add-back;do not transfer1-MNA support to SAM'],
 ['M04','E04','UCKL1',2023,'https://pmc.ncbi.nlm.nih.gov/articles/PMC10363437/',False,'Original indexed Results and Fig4 caption;PMC page challenged','RKO metabolomics;HCT116 supplement experiment;DLD1 motif mutant','KD;deletion aa132-136 uridine-binding motif','UMP/CMP not significantly decreased in tested KD condition;uridine-substrate kinase activity reduced for mutant','UMP/CMP did not rescue measured phenotypes;mutant retained suppression of erastin-related lipid peroxidation/ROS','Reduced activity does not mean zero activity;negative UMP/CMP rescue does not negate all uridine biology','Separate catalytic uridine flux from noncanonical NRF2/ferroptosis function'],
 ['M05','NA','UCKL1',2022,'https://pmc.ncbi.nlm.nih.gov/articles/PMC9246348/',True,'Original indexed abstract/results overview;not complete figure audit','Purified UCKL1;YAC1 lymphoma/K562 leukemia models','Purified-enzyme assay;KD in non-CRC models','Direct uridine/cytidine phosphorylation reported','Enzymatic capability supported;CRC uridine-mediated phenotype not established','Biochemical capacity and necessity for a particular CRC phenotype are different claims','Use catalytic versus noncatalytic evidence as separate questions'],
 ['M06','NA','SLC6A6',2026,'https://www.nature.com/articles/s42255-026-01455-6',True,'Publisher abstract/figure titles/data statement;subscription-limited full text','Liver/lung cancer and hepatocyte context per figure titles','SLC6A6 deficiency;localization regulation reported','Mitochondrial taurine and translation reported','Subcellular import hypothesis;CRC mediation not established','Not COAD replication;Fig7a/7j source data corrected2026-07-14;not independently reanalyzed','Assess compartment relevance before applying model to COAD']]
 save('mechanism_evidence.tsv',[dict(zip(fields,x))for x in records])
 correction=dict(parent_review_id='M06',url='https://doi.org/10.1038/s42255-026-01585-x',published='2026-07-14',reading_depth='Publisher-indexed correction text',scope='Fig7a,Fig7j and associated Source Data corrected for data preparation errors',effect_on_use='Use corrected version;no independent source-data or figure reanalysis performed',independent_study=False)
 write('correction_record.json',json.dumps(correction,indent=2)+'\n')
 decisions={
 ('NNMT','KEGG:C02918'):('M02;M03','CRC_PRODUCT_MEASUREMENT_AND_CONDITIONAL_PARTIAL_RESCUE','5-FU-treated NNMT-KD cells have reported1-MNA partial rescue;stromal-line secretion/addition evidence is distinct','Patient-relevant dose,cell-source mediation and independent matched patient relationship'),
 ('NNMT','KEGG:C00019'):('M02;M03','SAM_MEDIATION_NOT_ESTABLISHED_BY_PRODUCT_STUDIES','1-MNA rescue cannot be assigned to SAM;reaction membership is not SAM-pool mediation','Measure SAM/SAH and relevant methylation after perturbation in matched model'),
 ('SLC6A6','KEGG:C00245'):('M01;M06','CRC_UPTAKE_SUPPORTED_COMPARTMENT_MECHANISM_UNRESOLVED','CRC uptake/survival evidence retained;2026 mitochondrial mechanism is other-context and corrected','CRC-specific subcellular pool,localization/transport control and metabolite-mediated phenotype'),
 ('UCKL1','KEGG:C00299'):('M04;M05','CATALYTIC_CAPACITY_SUPPORTED_CRC_FERROPTOSIS_NONCANONICAL','Enzyme can phosphorylate uridine;CRC ferroptosis findings do not establish uridine-consumption mediation','Quantify uridine flux and separate catalytic from NRF2-related effects')}
 delta=[];coverage=[]
 for i,x in enumerate(old,1):
  key=(x['gene'],x['metabolite_key']);reviewed=key in decisions
  coverage.append(dict(source_row_1based=i,gene=x['gene'],metabolite_key=x['metabolite_key'],this_batch_status='DONE'if reviewed else'NOT_RUN',scope='Focused mechanism interpretation'if reviewed else'Prior evidence retained;not newly reread',prior_delta='results/COAD/07_INTEGRATION/20260920T105400Z_relation39_feasibility_v1/relation39_next_evidence.tsv'))
  if reviewed:
   ids,label,meaning,gap=decisions[key];delta.append(dict(source_row_1based=i,gene=x['gene'],metabolite_key=x['metabolite_key'],metabolite_name=x['metabolite_name'],review_ids=ids,reviewed_label=label,interpretation=meaning,remaining_gap=gap,independent_patient_validation='NOT_RUN',arrangement='UNCHANGED_9_21_5',new_p_value='NA',new_q_value='NA'))
 save('relation4_mediation_delta.tsv',delta);save('relation39_review_coverage.tsv',coverage)
 write('README_CN.md','''# COAD：四条关系的代谢介导证据定向核查

## 本轮问题

取数尚未取得新矩阵，本批独立推进NNMT—1-MNA、NNMT—SAM、SLC6A6—牛磺酸、UCKL1—尿苷。要分清基因干预、准确代谢物读数和代谢物介导/救援，不追加表达显著性筛选。

## 输入与范围

原39关系仍全保留。本批3基因4关系，4篇既有研究定向复读（E01/E04/E07/E08）、2篇新纳入项目研究及1份更正通知。不是7篇独立研究，也不是4条新验证关系。其余35关系本批NOT_RUN表示未新增阅读，不取消已有证据；覆盖表逐行记录。不同文章可能共享研究系列，不按论文数推断独立临床复现。

## 实际结果

|关系|本次更明确的证据|解释边界|
|---|---|---|
|NNMT—1-MNA|2016原文Fig6–8有产物测量；HT29敲低+5-FU背景下，2/4mM1-MNA部分救援ROS、凋亡和药物IC50|不是仅有外源添加，但剂量与治疗背景不能直接等同当前患者浓度|
|NNMT—SAM|上述1-MNA证据不提供SAM特异介导或救援|不得把同一酶的产物证据复制给甲基供体关系|
|SLC6A6—牛磺酸|CRC研究有直接摄取测量；新纳入2026研究提示线粒体牛磺酸区室机制|2026主要肝/肺背景，非COAD复现；有更正记录，未独立重审源数据|
|UCKL1—尿苷|2022纯化蛋白研究支持尿苷磷酸化能力；2023 CRC研究提供UMP/CMP未救援与低催化活性突变体仍保留铁死亡保护的证据|“能催化尿苷”与“该CRC表型依赖尿苷代谢”不可混同；突变体是活性降低，不能擅自写成完全失活|

NNMT 2021全文还明确了模型：CCD-18Co结肠成纤维细胞系，包含1-MNA测量和外源添加实验。不能把它直接改写成患者来源CAF特异性因果验证，也不能把分泌/添加当作已完成全部条件培养基介导链。

SLC6A6 2026论文的2026-07-14更正涉及Fig7a/7j及对应Source Data，原因是数据准备错误。这里记录更正与读取限制，不据此认定整篇无效或已经独立验证修正无误。

出处：[NNMT 2016](https://www.oncotarget.com/article/9962/text/)、[NNMT 2021](https://www.jcancer.org/v12p6170.htm)、[SLC6A6 2014](https://www.nature.com/articles/srep04852)、[SLC6A6 2026](https://www.nature.com/articles/s42255-026-01455-6)、[更正](https://doi.org/10.1038/s42255-026-01585-x)、[UCKL1 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10363437/)、[UCKL1 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9246348/)。

## 新手解释

部分救援可以增强特定条件下的机制支持，但不会自动证明组织中的相关系数由同一机制产生。组织总量也不能替代细胞内某区室的代谢物含量。这一批增加的是“该怎样解释、还缺哪一步”，不是新的患者阳性名单。

## 限制/反证

NNMT两篇核查了可读正文重点段落；UCKL1等部分核查依赖原始来源索引正文/图注；SLC6A6 2026完整正文受订阅限制，只读取摘要、图题、数据声明与更正。没有逐图/原始数据重分析、完整系统综述或全面撤稿数据库审计。没有新患者矩阵、P/q或实验。不能从未在所读段落确认的救援推断文献绝无救援。

## 当前决定

四关系解释核查批次DONE，独立患者关联与精确机制验证仍PARTIAL。保留974候选、39关系、35基因及9/21/5。原统计和既有证据不覆盖；本轮无作者联系，取数草稿仍未发送。B4GALT2上一批的更新继续有效。

## 下一步

NNMT优先区分治疗/细胞来源与1-MNA生理剂量，并将SAM独立处理；SLC6A6先核对COAD模型的膜/线粒体区室适用性；UCKL1区分催化通量与非催化保护。其余关系按39条台账逐批推进，HDC宿主/肿瘤来源和BCAT2营养条件可作为下一批。数据到达后按既定取数合同锁定可评估集合，不能按新P值选择。

## 复现命令

```text
python code/coad/build_mediation_review_v1.py
python tools/check_repository.py
```

脚本仅重建人工核实的公开解释表，不重跑网络检索、患者统计或实验。NA表示本批无新统计；输入与代码哈希见source_manifest.tsv。
''')
 manifest=[dict(source_id=x[0],path_or_url=x[4],version='2026-09-21 focused reading',sha256='NA',depth=x[6])for x in records]
 manifest.append(dict(source_id='M06_correction',path_or_url=correction['url'],version=correction['published'],sha256='NA',depth=correction['reading_depth']))
 for p in [SOURCE,Path(__file__)]:manifest.append(dict(source_id='local',path_or_url=str(p.relative_to(R)),version='exact bytes',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),depth='Input/code'))
 save('source_manifest.tsv',manifest)
 spec=dict(run_id=RUN,analysis_version='COAD_mediation4_v1',statistical_unit='NA',test_family='NA',seed='NA',new_tests=False,new_data=False,email_sent=False,genes=3,relations=4,previous_studies_reread=4,new_research_records=2,correction_records=1,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip())
 write('analysis_spec.json',json.dumps(spec,indent=2)+'\n')
 assert len(delta)==4 and len(coverage)==39 and len({x['gene']for x in delta})==3
 assert [(x['gene'],x['metabolite_key'])for x in old]==[(x['gene'],x['metabolite_key'])for x in coverage]
 assert all(x['new_p_value']==x['new_q_value']=='NA'for x in delta)
 write('validation.json',json.dumps(dict(status='PASS',scope='Review coverage/keys/order/no new statistical values;not biological validation',original_relations=39,reviewed_relations=4,not_newly_reviewed_relations=35,reviewed_genes=3,research_records=6,correction_records=1,arrangement_unchanged=True,input_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest()),indent=2)+'\n')
 print(str(OUT))
if __name__=='__main__':main()
