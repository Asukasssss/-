"""Build public acquisition contract and versioned B4GALT2 interpretation;no patient tests."""
import csv,json,hashlib,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[2];RUN='20260921T025931Z_acquisition_review_v1'
OUT=R/'results/COAD/07_INTEGRATION'/RUN
OLD=R/'results/COAD/06_EXTERNAL/20260920T091100Z_cell_paired35_v1/relations39_integrated.tsv'
def read(p):
 with p.open(encoding='utf-8-sig',newline='')as f:return list(csv.DictReader(f,delimiter='\t'))
def save(name,rows):
 with (OUT/name).open('w',encoding='utf-8',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def write(name,text):(OUT/name).write_text(text,encoding='utf-8')
def main():
 OUT.mkdir(parents=True,exist_ok=True);old=read(OLD);assert len(old)==39 and len({x['gene']for x in old})==35
 sources=[
 ('D01','Chen2025','https://link.springer.com/article/10.1186/s12943-025-02359-x','Publisher correspondence/data declaration;prior two-DOCX audit reused','PARTIAL','AUTHOR_DATA_REQUEST_READY_NOT_SENT','13 RNA patients;11 protein/metabolite patients reported;actual tumor overlap unverified','Processed full matrices and exact specimen/technical-replicate map required'),
 ('D10','Li2026','https://link.springer.com/article/10.1186/s12885-026-15978-4','Publisher methods/results/data availability/correspondence','PARTIAL','AUTHOR_CLARIFICATION_READY_NOT_SENT','40 patients overall;proteomics selected24 paired tissue sets,12CC/12RC per article;not40 matched COAD cases','Clarify treatment timing/no-chemotherapy group;exact matching and colon subset;available on reasonable request'),
 ('D11','Wang2025','https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00010','Original publisher indexed methods/data availability;direct page challenge','NOT_EVALUABLE','EXCLUDE_FROM_PRIMARY_TUMOR_VALIDATION','Serum;OMIX008632 and PXD059480 reported raw-data deposits','Different biological material;no download;accessions do not prove processed tumor matrix availability')]
 save('dataset_acquisition_delta.tsv',[dict(source_id=a,study=b,url=c,read_depth=d,status=e,action=f,reported_design=g,reason=h,actual_feature_coverage='NOT_VERIFIED',new_association='NOT_RUN',contact_sent=False)for a,b,c,d,e,f,g,h in sources])
 contacts=[('D01','Ningning Zhao','nnzhao@cmpt.ac.cn'),('D01','Sheng Dai','daimd@zju.edu.cn'),('D01','Fan Zhang','fzhang@cmpt.ac.cn'),('D01','Hongcang Gu','gu_hongcang@cmpt.ac.cn'),('D10','Tianbao Xiao','prof_xiaotianbao@163.com'),('D10','Jiang Chen','chenjiang05v92@163.com'),('D10','Tao Yang','yangtao901212@163.com')]
 save('contact_register.tsv',[dict(source_id=a,name=b,published_email=c,source_url=next(x[2]for x in sources if x[0]==a),verified_date='2026-09-21',status='NOT_RUN',reason='Public paper contact;draft prepared;not sent;no agreement accepted')for a,b,c in contacts])
 fields=[('metabolite_matrix','feature_id,specimen_id,value,missing_flag','Full processed matrix without differential-significance filtering;no raw mass spectra'),('feature_annotation','feature_id,name,structure_id,annotation_confidence,adduct_or_peak_context','Exact feature identity;do not substitute related metabolites/isomers'),('RNA_or_protein_matrix','gene_or_protein_id,specimen_id,value,scale','At least one modality;RNA and protein analyzed and reported separately'),('specimen_map','anonymous_patient_id,tissue_id,omics_sample_id,aliquot_id,technical_repeat_id','Explicit correspondence;never infer pairing from order/names;one independent patient cannot be multiplied by repeats'),('metadata','tumor_or_adjacent,colon_or_rectum,site,pretreatment,collection_time,batch','Establish tissue scope,treatment timing,independence and technical context'),('processing','normalization,missing_value_handling,replicate_aggregation,QC,units','Do not invent transformations or repeat pooling rules'),('provenance','repository_accession,version,license_or_access_terms,file_checksum','Server-only source storage;access agreement if supplied requires review')]
 save('requested_data_contract.tsv',[dict(file_type=a,minimum_fields=b,purpose_or_rule=c,status='NOT_RUN')for a,b,c in fields])
 save('relation39_request_scope.tsv',[dict(source_row_1based=i,metabolite_key=x['metabolite_key'],metabolite_name=x['metabolite_name'],gene=x['gene'],actual_external_coverage='NOT_VERIFIED',include_by_new_p_value=False)for i,x in enumerate(old,1)])
 oldfunction=next(x for x in read(R/'results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1/gene_function_review.tsv')if x['gene']=='B4GALT2')
 # Exact historical row retained separately as audit reference, without changing its source.
 save('B4GALT2_original_record.tsv',[oldfunction])
 save('B4GALT2_review.tsv',[dict(gene='B4GALT2',metabolite_key='KEGG:C00015',original_arrangement='暂挂',current_arrangement='暂挂',arrangement_9_21_5='UNCHANGED',reviewed_label='CRC_PERTURBATION_SUPPORTED_CONTEXT_DEPENDENT_UDP_MEDIATION_UNRESOLVED',old_evidence='E23;2016 SW480 wild-type/A146V expression;measured migration not changed',new_evidence='F05 reused;2026 CRC interference and MUC20 glycosylation/protein findings;not a newly discovered source in this batch',comparison='Different cells,perturbation and endpoints;not a matched-condition contradiction',reviewed_hold_reason='Do not claim absence of CRC intervention evidence;retain hold on general therapeutic direction and exact UDP axis',UDP_measurement='NOT_VERIFIED',UDP_rescue='NOT_VERIFIED',cell_context='CRC epithelial perturbation does not establish stromal mechanism',review_depth='Original indexed results/figure captions;no independent raw figure/source-data audit',next_action='Resolve model-specific direction and on-target/catalytic rescue;measure exact UDP separately from UDP-galactose',old_url='https://www.nature.com/articles/srep23642',new_url='https://www.nature.com/articles/s42003-026-10017-1')])
 common='''我们正在开展结直肠癌组织中代谢物与基因表达关系的研究，希望在独立标本中检验预先确定的关系。请问能否提供存储库编号或可申请获取的处理后配套数据？

希望取得未经显著性筛选的代谢特征定量矩阵及身份注释、配套RNA或蛋白矩阵，以及连接去标识化患者代号、组织标本、各组学样本和技术重复的对应表。RNA和蛋白不必同时提供。另请提供组织类别、结肠/直肠部位、取材前治疗和取材时间、分析批次、数值尺度、归一化、缺失值处理及重复合并说明。

我们不需要姓名等直接身份信息或原始质谱文件。如果完整矩阵暂不便共享，可先提供全特征注释和字段字典，以判断覆盖和申请要求。如需数据使用协议或机构审批，请告知流程。感谢支持。

【真实姓名】
【单位】
【机构邮箱】
'''
 write('CHEN_DATA_REQUEST_CN.md','''# Chen 2025取数邮件草稿

状态：草稿已完成，未发送。发件人身份须由真实发送者补齐；未同意任何协议。优先联系论文通讯作者，其他通讯作者仅记录，不自动群发。

收件参考：Ningning Zhao <nnzhao@cmpt.ac.cn>；Sheng Dai <daimd@zju.edu.cn>

主题：申请独立研究所需的配套组织多组学数据（10.1186/s12943-025-02359-x）

赵老师／戴老师您好：

已阅读贵文及两份补充材料。希望进一步确认13名RNA患者与11名蛋白/代谢患者的对应关系、结肠/直肠构成，以及代谢组三次技术重复的标记与处理规则。

'''+common)
 write('LI_DATA_REQUEST_CN.md','''# Li 2026取数与设计澄清邮件草稿

状态：草稿已完成，未发送。收件参考：Tianbao Xiao <prof_xiaotianbao@163.com>；另两位通讯作者见contact_register.tsv，未自动抄送。

主题：咨询配套组织代谢组/蛋白组数据及治疗分组（10.1186/s12885-026-15978-4）

肖老师您好：

已阅读贵文。为准确复用数据，想请教方法中术前新辅助治疗后的取材描述与结果中未化疗组的对应关系：各组具体治疗、取材时点和患者数如何定义？正文报告蛋白组选择24组配对组织，其中结肠、直肠各12组；这些与40人代谢组池的精确交集、技术重复及原发结肠病例范围也希望获得说明。

'''+common)
 write('SERVER165_NEXT_TASK_CN.md','''# COAD配套数据到达后的执行入口

不再拆读已核查的Chen两份DOCX。源矩阵、临床个体行和完整映射仅存server165。收到合法可用的数据入口后，在新的COAD/B运行目录建独占锁，下载处理后数据到新的data/candidates子目录，记录来源、许可和哈希；不下载原始质谱。

1. 先核实身份、材料、精确标本对应和独立患者。技术重复按明确的作者处理说明与预定规则处理；肿瘤/癌旁不混合扩大肿瘤内样本量；结肠与直肠分开。
2. 面向原39关系/23特征/35基因生成完整覆盖去向。至少需要代谢物加RNA或蛋白之一；两层分别报告，不挑较小P值。仅有注释/字段字典时标记覆盖预审，不标数据可分析。
3. 在查看新相关系数/P值前，依据身份、缺失、质量及有效独立样本规则，锁定可评估集合、技术重复规则、统计方法和各层BH范围。其余关系保留不可评估原因；不要求39条齐全才开跑。
4. Chen按小样本复核规划，不能把11人膨胀为33/66。Li的24组配对及12组结肠仍需原始对应与治疗澄清，不能直接当作已核实有效n。
5. 数值完成后验收，公开汇总/代码/参数/状态，更新统一阶段索引并推送核对远端；原CAMP与临床/表达统计不覆盖。

当前数据接收、覆盖核验、关联计算均NOT_RUN。邮件尚未发送，等待发送者身份和明确发送指令；本文件不是定时或后台自动执行任务。
''')
 write('README_CN.md','''# COAD取数准备与B4GALT2暂挂理由评议

## 本轮问题

Chen补充文件核查已经结束，本批转为具体取数准备，登记两个入口，并对B4GALT2的既有新旧功能记录形成解释决定。不追加同类表达筛选。

## 输入与范围

复用原39关系、35基因、23特征及两份补充文件核查。用户提供的是文字交接，未收到实体取数交接ZIP；本包根据文字要求和本轮核实的原始来源重新整理。来源读取深度见dataset_acquisition_delta.tsv；B4GALT2为E23/F05复核，不计作新增独立文献。未登录服务器重跑、未取得新矩阵。

## 实际结果

05_FUNCTION：B4GALT2原记录与本次评议并列保存。2016的SW480野生型/A146V表达未改变所测迁移；2026 CRC干扰支持表型及MUC20相关功能。两者模型/操纵不同，不构成同条件反向重复。结论是不能再笼统说缺少CRC干预支持；仍暂挂通用治疗方向和精确UDP介导判断，9/21/5不变。新研究本轮仅复核原始来源索引正文/图注，未重新审计图像原始数据。来源：[2016](https://www.nature.com/articles/srep23642)、[2026](https://www.nature.com/articles/s42003-026-10017-1)。

06_EXTERNAL：Chen转作者取数，草稿和全字段请求已准备，未发送。Li2026列为第二取数候选：总体40人，正文蛋白组选择24组配对组织、结肠/直肠各12组；不是40名配套COAD患者。治疗后取材与未化疗组描述需澄清，实际同标本交集未知。Wang2025为血清组学，论文报告OMIX008632/PXD059480原始数据编号；排除出当前原发肿瘤内部验证，无下载。来源：[Chen](https://link.springer.com/article/10.1186/s12943-025-02359-x)、[Li](https://link.springer.com/article/10.1186/s12885-026-15978-4)、[Wang](https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00010)。

07_INTEGRATION：39关系的取数范围保持原键/顺序；7类字段请求与服务器接收后流程已形成。原9入口保持历史，新增D10/D11及D01状态增量，不把它们称新增3个可用队列。目前新增可直接运行的验证集合仍为0。

## 新手解释

取数准备完成不等于已取得数据。RNA与蛋白有一层配套即可判断该层的可评估关系，两层不能互相替代；也不必凑齐39条，但不能看过新P值再挑关系。B4GALT2有基因功能证据与UDP是否解释表型是两个问题。

## 限制/反证

没有发送邮件、同意协议或取得新矩阵。Li实际样本对应/治疗/覆盖未核实，Wang材料不符。B4GALT2的UDP含量、通量和UDP依赖救援仍未核实，CRC上皮模型不能直接解释基质来源。全部临床和表达P/q保持不变；本批无新检验。没有对其他34基因重新系统检索。

## 当前决定

本批取数文稿、入口登记及解释评议DONE；独立代谢验证PARTIAL/NOT_RUN。保持974候选、原39关系和9/21/5。旧暂挂表不覆盖，通过reviewed_label/reviewed_hold_reason给出更新理由。COAD仍由B负责，不改PDAC。

## 下一步

由真实发送者补齐姓名、单位、机构邮箱后发送取数请求；本轮未获明确发送指令，因此只交付草稿。收到资料后按SERVER165_NEXT_TASK_CN.md核验并预先锁定实际检验集合。功能推进可继续处理具体代谢读数、细胞模型及救援缺口，不依赖Chen回复才开展。

## 复现命令

```text
python code/coad/build_acquisition_handoff_v1.py
python tools/check_repository.py
```

脚本重建已人工整理的公开核查与请求材料，不执行网络搜索、邮件发送、下载或患者分析。输入哈希、读取范围和验证见随附manifest/spec/validation。
''')
 spec=dict(run_id=RUN,analysis_version='COAD_acquisition_review_v1',scope='D01 acquisition transition,D10/D11 screening,B4GALT2 interpretation delta',statistical_unit='NA',test_family='NA',seed='NA',new_tests=False,email_sent=False,agreements_accepted=False,matrices_acquired=False,original_arrangement='UNCHANGED_9_21_5',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip())
 write('analysis_spec.json',json.dumps(spec,indent=2)+'\n')
 manifests=[dict(source_id=x[0],path_or_url=x[2],version='2026-09-21 read',sha256='NA',depth=x[3])for x in sources]
 for ident,url in [('E23','https://www.nature.com/articles/srep23642'),('F05','https://www.nature.com/articles/s42003-026-10017-1')]:manifests.append(dict(source_id=ident,path_or_url=url,version='2026-09-21 focused reread',sha256='NA',depth='Original indexed results/figure captions;no new source-data audit'))
 for p in [OLD,Path(__file__),R/'results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1/gene_function_review.tsv']:
  manifests.append(dict(source_id='local',path_or_url=str(p.relative_to(R)),version='local exact bytes',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),depth='Versioned input/code'))
 save('source_manifest.tsv',manifests)
 assert len(read(OUT/'relation39_request_scope.tsv'))==39
 assert [(x['gene'],x['metabolite_key'])for x in old]==[(x['gene'],x['metabolite_key'])for x in read(OUT/'relation39_request_scope.tsv')]
 assert read(OUT/'B4GALT2_original_record.tsv')==[oldfunction]
 write('validation.json',json.dumps(dict(status='PASS',scope='Keys/order/counts and exact original B4GALT2 row;not biological validation',relations=39,genes=35,metabolites=len({x['metabolite_key']for x in old}),contact_records=7,emails_sent=0,new_patient_tests=0,new_ready_datasets=0,old_arrangement_preserved=True,input_sha256=hashlib.sha256(OLD.read_bytes()).hexdigest()),indent=2)+'\n')
 print(str(OUT))
if __name__=='__main__':main()
