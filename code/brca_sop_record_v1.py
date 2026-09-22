"""Index published BRCA aggregates and append interpretations; no patient tests."""
from pathlib import Path
import csv,json,hashlib,subprocess,math,collections
R=Path(__file__).resolve().parents[1]
RUN='20260922T112000Z_sop_record_v1'
O=R/'results/BRCA/07_INTEGRATION'/RUN
O.mkdir(parents=True,exist_ok=True)
BASE=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
B='https://github.com/Asukasssss/-/blob/'+BASE+'/'
sources={
 'metabolites':'results/BRCA/04_ROBUSTNESS/20260921T112611Z_metabolite318_paired_v1/paired_metabolite318.tsv',
 'mapping':'results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/direct_relations.tsv',
 'mapping_status':'results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/metabolite190_mapping_status.tsv',
 'RNA':'results/BRCA/03_PATIENT/20260921T150000Z_mapping262_patient_v1/paired_RNA150.tsv',
 'association':'results/BRCA/03_PATIENT/20260921T150000Z_mapping262_patient_v1/CAMP_relations262.tsv',
 'relations':'results/BRCA/07_INTEGRATION/20260921T154000Z_mapping262_integration_v1/relations262_comparison.tsv',
 'genes':'results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/genes156_subtype_appended.tsv',
 'sc39':'results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_cross_study_comparison.tsv',
 'subtype':'results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/genes156_subtype_comparison.tsv',
}
def read(p):
 with (R/p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows,fields=None):
 fields=fields or list(rows[0])
 with (O/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def txt(name,s):(O/name).write_bytes((s.rstrip()+'\n').encode('utf-8'))
def js(name,obj):txt(name,json.dumps(obj,ensure_ascii=False,indent=2))
def number(x):
 try:return float(x)
 except (ValueError,TypeError):return float('nan')
def tier(p,q,status):
 if status!='DONE':return status or 'NOT_EVALUABLE'
 if number(q)<.05:return 'FDR_SUPPORTED'
 if number(p)<.05:return 'NOMINAL_P_ONLY'
 return 'NOT_NOMINAL_SIGNIFICANT'
D={k:read(v) for k,v in sources.items()}
genes=D['genes'];relations=D['relations'];rn={x['gene']:x for x in D['RNA']}
assert len(genes)==156 and len({x['gene'] for x in genes})==156
assert len(relations)==262 and len({x['relation_key'] for x in relations})==262
assert len(D['mapping_status'])==190 and len(rn)==150
notes={
 'ASNS':('与GLS共同保留为谷氨酰胺专题候选；已有患者关联与Asns自身干预RNA资料。','条件相关不是独立通量；小鼠再培养细胞RNA不能证明人类谷氨酰胺介导或GLS功能补偿。','保留专题及未支持结果；当前主线停在来源，不再追加同队列模型。'),
 'GLS':('FUSCC条件模型中负关联较突出，但不据q大小判定优于ASNS。','组织RNA不是酶活；单细胞内皮较高不证明患者关联来自内皮。','将自身干预与患者关系分开记载，等待另行授权的功能问题。'),
 'GLUL':('保留弱正向条件关联背景。','后选择六项模型的q不替代原全关系q，不表示突然发现强机制。','不单独升级主线。'),
 'GLS2':('保留当前患者关联支持不足的结果。','不能称生物学阴性对照。','按全候选保留来源与功能背景。'),
 'SLC6A8':('肌酸配对背景与FUSCC关联构成值得讨论的运输线索。','丰度和RNA不能确定运输方向、活性或治疗干预方向。','保留为运输类候选。'),
 'GPCPD1':('CAMP内部负关联对已做组成代理调整较稳定，可作为内部患者联系参照。','外部关联支持不足；髓系表达较高不等于唯一作用细胞；调整不是去除全部混杂。','保留内部支持与外部未支持并列。'),
 'GPI':('保留CAMP内部关联，但Tang已可评估且未支持正向关联。','不能再统一标外部缺测；小队列宽区间也不能证明反向作用。','不以换模型追求显著，保留队列不一致。'),
 'NNMT':('1-MNA患者线索需结合成纤维背景讨论。','组成调整后关联减弱；不能说已证明细胞比例导致，也不能把表达来源当代谢物来源。','保留跨细胞背景问题。'),
 'MDH1':('Tang新增较强苹果酸关联，作为小队列新线索保留。','FUSCC点估计接近零；RNA肿瘤较低与肿瘤内部正相关是不同问题；不可混用MDH2证据。','不升为全池第一名。'),
 'PNP':('次黄嘌呤和鸟嘌呤关系分别保留正向线索。','同队列两代谢物不是独立验证；新版q变化不是新实验否定。','保留效应、区间与全部q版本。'),
 'PYCR1':('已有恶性上皮表达背景及CAF功能文献，作为限定专题。','最高表达不是唯一功能细胞；既有论文端点摘录不是新独立验证。','已有专题收束，需新问题再深入。'),
 'LYPLA1':('新增候选中患者线索与恶性上皮表达背景较清楚，值得优先讨论。','1-palmitoyl-GPC不是游离GPC；当前精确关系尚无可用外部检验，来源不证明脂质介导。','保留为讨论对象，不宣布唯一主线。'),
 'ENPP2':('配对RNA降低、脂质负关联和内皮来源背景均应保留。','可用值相关未获FDR支持；功能取决于作用细胞，不能从RNA降低推断应提高它。','标记敏感性限制及细胞背景。'),
 'ABHD12':('配对RNA升高，Wu/Pal均髓系最高且排名稳定。','表达来源不是乳腺癌功能验证。','作为髓系背景候选并列讨论。'),
 'LYPLA2':('配对RNA升高，两研究恶性上皮来源一致。','已有功能记录的模型身份限制仍在；不能自动继承LYPLA1证据。','按自身证据讨论。'),
 'LPCAT1':('配对RNA升高并有既有功能背景。','Wu髓系、Pal内皮最高；HER2+中的B/髓系均值接近，不是稳定细胞特异。','不强制指定一个唯一作用细胞。'),
 'LPCAT4':('单细胞稳定ID匹配提供了来源背景。','单细胞身份解决不自动解除CAMP微阵列探针身份待核。','患者RNA及相关继续沿用待核状态。'),
}
gout=[];reader=[];judgments=[]
for g in genes:
 gene=g['gene'];a=rn.get(gene);edges=[x for x in relations if x['gene']==gene]
 current=gene in rn
 interpretation,limit,nextstep=notes.get(gene,('完整保留候选，按自身配对RNA、关系级患者结果和细胞来源讨论。','未被选作示例不表示无作用；来源与相关不证明功能。','当前完成来源记录，尚无足够依据安排独立机制课题。'))
 if not current:interpretation='历史候选保留，不在当前190项入口直接池；原证据不删除。'
 wu=g.get('sc39_Wu_top') if g.get('sc39_Wu_top') not in ('','NA',None) else g.get('rob_Wu2021_top_lineage') or g.get('sc_Wu2021_top_lineage')
 pal=g.get('sc39_Pal_top') if g.get('sc39_Pal_top') not in ('','NA',None) else g.get('rob_Pal2021_reprocessed_top_lineage') or g.get('sc_Pal2021_reprocessed_top_lineage')
 summary=dict(gene=gene,current_pool=current,history_only=not current,n_current_relations=len(edges),RNA_status=a['status'] if a else 'HISTORICAL_ONLY',RNA_n=a['n'] if a else 'NA',RNA_effect=a['effect'] if a else 'NA',RNA_p=a['p_value'] if a else 'NA',RNA_q=a['q_value'] if a else 'NA',RNA_positive_pairs=a['positive_pairs'] if a else 'NA',RNA_negative_pairs=a['negative_pairs'] if a else 'NA',RNA_tier=tier(a['p_value'],a['q_value'],a['status']) if a else 'HISTORICAL_ONLY',CAMP_evaluable=sum(x['CAMP_status']=='DONE' for x in edges),CAMP_nominal=sum(number(x['CAMP_p_value'])<.05 for x in edges),CAMP_FDR=sum(number(x['CAMP_q_value'])<.05 for x in edges),Wu_top=wu or 'NA',Pal_top=pal or 'NA',Wu_ER_top=g.get('Wu_subtype_ER+_top','NA'),Wu_HER2_top=g.get('Wu_subtype_HER2+_top','NA'),Wu_TNBC_top=g.get('Wu_subtype_TNBC_top','NA'),interpretation_cn=interpretation,limitation_cn=limit,next_action_cn=nextstep,RNA_source=sources['RNA'] if a else sources['genes'],relation_source=sources['relations'],cell_source=sources['genes'],judgment_type='INTERPRETATION_NOT_NEW_TEST',new_patient_statistics=False)
 reader.append(summary)
 z=dict(g);z.update({'record_v1_'+k:v for k,v in summary.items() if k!='gene'});gout.append(z)
 judgments.append({k:summary[k] for k in ['gene','current_pool','history_only','interpretation_cn','limitation_cn','next_action_cn','RNA_source','relation_source','cell_source','judgment_type']})
write('candidate_genes156_history_preserved.tsv',gout)
write('candidate_genes156_reader.tsv',reader)
write('interpretations156.tsv',judgments)
gi={x['gene']:x for x in reader};rout=[]
for r in relations:
 z=dict(r);s=gi[r['gene']]
 z.update(record_v1_mapping_stage_patient_status=r['patient_analysis_this_round'],record_v1_current_patient_status=r['CAMP_status'],record_v1_association_tier=tier(r['CAMP_p_value'],r['CAMP_q_value'],r['CAMP_status']),record_v1_source=sources['relations'])
 for k in ['RNA_status','RNA_n','RNA_effect','RNA_p','RNA_q','RNA_positive_pairs','RNA_negative_pairs','RNA_tier','Wu_top','Pal_top','interpretation_cn','limitation_cn']:z['record_v1_'+k]=s[k]
 rout.append(z)
write('candidate_relations262_history_preserved.tsv',rout)
# Standard prefix already exists in all three original scientific tables. Preserve strings exactly.
prefix=(R/'templates/statistical_result.tsv').read_text().strip().split('\t')
stats=[];mapping=[]
for name in ['metabolites','association','RNA']:
 for x in D[name]:
  assert all(k in x for k in prefix)
  z={k:x[k] for k in prefix};z.update(record_source_path=sources[name],record_source_commit=BASE,record_operation='EXACT_COLUMN_COPY_NO_RECALCULATION');stats.append(z)
 for k in prefix:mapping.append(dict(output='statistics_long.tsv',output_field=k,source_path=sources[name],source_field=k,transformation='identity_string_copy'))
write('statistics_long.tsv',stats)
write('paired_RNA85_P005_view.tsv',[x for x in D['RNA'] if number(x['p_value'])<.05])
workpool=[x for x in D['metabolites'] if x['analysis_type']=='paired_processed' and number(x['p_value'])<.05]
workpool.sort(key=lambda x:-max(number(x['higher_fraction']),number(x['lower_fraction'])))
assert len(workpool)==190
write('metabolite190_P005_direction_view.tsv',workpool)
write('field_mapping.tsv',mapping)
counts=[]
for key,rows in __import__('itertools').groupby(sorted(stats,key=lambda x:(x['record_source_path'],x['analysis_type'])),lambda x:(x['record_source_path'],x['analysis_type'])):
 a=list(rows);counts.append(dict(source_path=key[0],analysis_type=key[1],planned_rows=len(a),evaluable=sum(x['status']=='DONE' for x in a),P_lt005=sum(number(x['p_value'])<.05 for x in a),q_lt005=sum(number(x['q_value'])<.05 for x in a),status_counts=json.dumps(dict(collections.Counter(x['status'] for x in a))),family_n_evaluable=';'.join(sorted(set(x['family_n_evaluable'] for x in a)))))
write('statistical_summary.tsv',counts)
# Full historical inventory references existing tracked artifacts, not raw server files.
tracked=subprocess.check_output(['git','-c','core.quotepath=false','ls-files','results/BRCA'],cwd=R,encoding='utf-8').splitlines()
inventory=[]
for p in tracked:
 f=R/p
 if not f.is_file():raise RuntimeError('Missing tracked artifact '+p)
 inventory.append(dict(path=p,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),source_snapshot=BASE,url=B+p,record_scope='EXISTING_TRACKED_ARTIFACT_NO_NEW_SCIENTIFIC_REVIEW'))
write('historical_artifact_inventory.tsv',inventory)
stage_specs=[
 ('01','身份和配对','03_PATIENT','20260921T102429Z_camp_sample_identity_v1','作者表追溯45对；冲突标本单列；新相关使用60肿瘤作者病例','NEEDS_REVIEW','身份追溯完成，1个组织标签冲突仍未解决'),
 ('02','配对代谢物','04_ROBUSTNESS','20260921T112611Z_metabolite318_paired_v1','318项全量；45对；Wilcoxon＋整对bootstrap；P入口190项','DONE','作者data可用值不是独立核验仪器检测掩码'),
 ('03','直接映射','02_MAPPING','20260921T124136Z_mapping190_v2','190入口；93项直接映射；262关系；150当前＋6历史','DONE','映射不是功能或患者统计'),
 ('04','患者相关和配对RNA','03_PATIENT','20260921T150000Z_mapping262_patient_v1','Spearman置换P、bootstrap；RNA配对t；分别BH','PARTIAL','关系255/262、RNA148/150；缺失或身份待核保留'),
 ('05a','旧117单细胞来源','06_EXTERNAL','20260919T140105Z_scRNA117_v1','供者内均值后等权；全部旧候选','DONE','注释来源与队列覆盖见原报告'),
 ('05b','旧117来源稳定性','04_ROBUSTNESS','20260919T151200Z_all117_robustness_v2','供者标签bootstrap及两研究对照','DONE','描述性稳定性不是患者百分比'),
 ('05c','新增39来源','06_EXTERNAL','20260922T023000Z_sc39_v1','Wu/Pal；20细胞、3来源标签；1000次重采样','DONE','LPCAT4单细胞身份不能解除微阵列待核'),
 ('05d','全156分型背景','06_EXTERNAL','20260922T064625Z_sc156_subtypes_v1','供者等权的分层描述，无新增P/q','DONE','临床分型不是上皮细胞状态；Pal无可靠分型连接'),
 ('05e','全156表达UMAP','06_EXTERNAL','20260922T082327Z_sc156_umap_v1','表达可视化','DONE','UMAP不检验细胞类型富集'),
 ('05f','细胞类型对照UMAP','06_EXTERNAL','20260922T085000Z_sc156_umap_pairs_v1','相同坐标表达与细胞类型并排','DONE','图示不是机制证据'),
 ('05g','上皮专用UMAP','06_EXTERNAL','20260922T092000Z_epithelial_umap_v1','既有作者注释＋上皮探索嵌入','DONE','可选展示；不称新聚类或患者亚型发现'),
 ('06','整合','07_INTEGRATION','20260921T154000Z_mapping262_integration_v1','全量关系与基因历史保留，本轮再接来源及解释','DONE','外部与功能属于已存附加证据，不是来源主线门槛'),
 ('附加','功能文献117','05_FUNCTION','20260919T120437Z_functional_review_v4','全候选定向检索及证据分层','DONE','不等于所有原文全文逐图审计'),
 ('附加','功能文献新增39','05_FUNCTION','20260921T152000Z_new39_literature_v1','新增基因功能背景','DONE','遵循模型适用范围，不据文献命中证明功能'),
 ('附加','外部262关系','06_EXTERNAL','20260921T151000Z_mapping262_external_v1','FUSCC126、Tang235可评估；旧q族保留','PARTIAL','FUSCC64条主关系RNA下载受阻；Oslo连接未解决'),
]
stage=[]
for step,title,sid,run,method,status,limit in stage_specs:
 p=f'results/BRCA/{sid}/{run}/README_CN.md';assert (R/p).exists(),p
 stage.append(dict(workflow_step=step,title_cn=title,stage_id=sid,source_run=run,method_and_scope_cn=method,status=status,limitation_cn=limit,report_path=p,report_url=B+p))
write('workflow_stage_index.tsv',stage)
write('historical_stage_registry.tsv',read('coordination/stages/BRCA.tsv'))
manifest=[dict(source_id=k,path=p,sha256=hashlib.sha256((R/p).read_bytes()).hexdigest(),source_commit=BASE,operation='READ_AGGREGATE_ONLY') for k,p in sources.items()]
manifest.append(dict(source_id='record_script',path='code/brca_sop_record_v1.py',sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_commit='DELIVERY_COMMIT_SEE_GIT_HISTORY',operation='INDEX_AND_INTERPRET_NO_PATIENT_TESTS'))
write('source_manifest.tsv',manifest)
js('analysis_spec.json',dict(version='sop_record_v1',run_id=RUN,source_commit=BASE,operation='aggregate_index_and_interpretation',new_statistical_tests=0,new_P_q=0,method_defaults_applied_retroactively=False,source_sop_branch='docs/camp-source-sop-v1-20260922',source_sop_commit='6dd0a73',source_sop_merge_status='NOT_MERGED',scope='BRCA only; single-cell stops at expression source',historical_q_preserved=True,FUSCC_BH='262 planned with internal missing P=1; published missing NA',Tang_BH='235 evaluable',CAMP_association_BH='255 evaluable separately per primary/availability',RNA_BH='148 evaluable',private_data_read=False))
for old,new in zip(genes,gout):assert all(new[k]==v for k,v in old.items())
for old,new in zip(relations,rout):assert all(new[k]==v for k,v in old.items())
assert sum(x['current_pool'] for x in reader)==150
assert sum(x['history_only'] for x in reader)==6
assert sum(number(x['p_value'])<.05 for x in D['RNA'])==85
js('validation.json',dict(status='DONE',scope='structural_and_value_preservation_not_patient_reanalysis',historical_artifacts_indexed=len(inventory),all_indexed_files_exist=True,gene_rows=156,relation_rows=262,gene_source_columns=len(genes[0]),relation_source_columns=len(relations[0]),all_historical_cells_preserved=True,standard_prefix_exact_copy=True,RNA_P_lt005=85,no_new_patient_tests=True,not_verified=['server patient matrices','all historical literature claims','full historical statistical recomputation','exhaustiveness of uncommitted server-only outputs']))
txt('INTERPRETATION_CN.md','''# BRCA现有结果：证据判断与待解决问题

本记录是对仓库已发布结果的解释，不是新增机制验证，也不是逐篇文献重新审核。逐基因事实见candidate_genes156_reader.tsv和来源表；解释单独见interpretations156.tsv。没有综合打分或强行排第一名。

## 我对整条主线的判断

候选产生过程已从原117基因扩展到190项配对P<0.05代谢物对应的150个直接候选，历史并集156。入口是探索性P值，不代表190项均有FDR支持；映射成功也不升级统计证据。配对变化、肿瘤内部相关和单细胞来源回答不同问题，不能作为三次独立验证相乘。

配对RNA148可评估、85个P<0.05、79个q<0.05。展示可按P保留，但单细胞覆盖完整候选，不能把RNA或相关不显著的对象排除。真正有用的产物是完整关系与细胞背景对照，而不是新的显著基因交集。

## 值得讨论的方向及我为什么这样判断

- LYPLA1：新关系的患者关联、可用值敏感性和RNA升高，加上两研究恶性上皮表达背景，形成较清楚的讨论基础。但1-palmitoyl-GPC不是游离GPC；没有精确关系外部检验，不能称乳腺癌脂质机制成立。
- ABHD12、LYPLA2：分别提供较一致的髓系和恶性上皮表达背景，应与LYPLA1并列讨论。RNA变化与定位不代替各自功能证据。
- ENPP2：RNA降低和内皮来源值得保留；脂质负关联的可用值分析未获FDR支持。不能因为RNA降低就推断需要激活。
- LPCAT1：RNA升高和已有功能背景有价值，细胞来源却不完全一致；HER2+下第一、第二类别差距极小，不能机械指定B细胞特异。
- ASNS/GLS：已有FUSCC患者联系及条件模型信息；Asns干预捕捉到了目标RNA降低和广泛转录反应，但未证明谷氨酰胺介导、两酶独立通量或GLS补偿。二者留作共同问题，不用q大小选胜者。
- SLC6A8：肌酸配对变化和外部关联是运输方向线索，尚不能决定干预方向。
- GPCPD1：CAMP内部关系对有限组成代理调整稳定，外部支持不足并列保留。髓系表达高不代表已解释该关联。
- GPI：CAMP支持仍在；Tang已可评估却未支持正相关，应降低“跨队列稳定关系”的表述，不能继续只写外部缺测。
- NNMT：患者线索与成纤维背景一起读，调整后减弱是背景敏感性信息，不是证明细胞比例造成。
- MDH1、PNP：Tang带来小队列线索；保留不一致与区间，不升为全池第一。不同PNP代谢物关系分开；MDH1不能继承MDH2证据。
- PYCR1：已有功能文献及不同细胞背景；恶性上皮表达最高不否定CAF作用。已有文献附件整理已足够，重复整理不构成新验证。
- 其他基因：全部在156行中保留。未列入上述例子不是阴性或淘汰，没有依据为每个基因编造独立机制问题。

## 不能被汇总掩盖的限制

1. 作者病例表支持45对，另有组织标签冲突；并非实物或基因型核验，新关联已排除冲突标本。配对数与肿瘤相关病例数不同。
2. 同CAMP内发现后再相关属于同队列探索。源data表有值的敏感性不是重新核验仪器缺失。
3. LPCAT4单细胞稳定ID匹配只解决该数据层，微阵列探针仍待核。CHKB覆盖缺项不是阴性。
4. FUSCC已存覆盖受下载限制；Tang仅20例且新版235项校正。FUSCC与Tang的q定义不能混作排名。
5. 单细胞供者标签沿用研究定义；表达最高、重采样稳定、亚型描述和UMAP均不证明功能或代谢物来源。Pal可靠亚型连接未接入。
6. 旧NOT_RUN等字段是旧阶段快照。新的record_v1字段明确本轮视图，不能把旧字段当当前运行状态。
7. 本轮不宣称检查服务器所有未发布文件；历史清单覆盖本次基准提交中所有已跟踪BRCA结果。源矩阵和患者级数据仍只存服务器。

## 当前决定

将“发现至单细胞来源”的结果固定为可讨论的一版。当前150与历史6全部保留；P和q分别呈现；缺测与未支持分开。外部/功能历史作为补充证据索引，不阻塞本主线。当前不新增通讯、拟时序、通路机制、分子对接或无止境协变量模型。
''')
lines=['# BRCA按统一规范归档与解释 v1','','## 本轮问题','将已发布BRCA结果按六步发现至单细胞来源规范记录，保留全部历史并追加可追溯解释。','','## 输入与范围',f'基准提交：`{BASE}`。只读仓库公开汇总；不读服务器矩阵，不新增P/q。公共SOP为6dd0a73分支版本，未合并main。','','## 实际结果',f'历史文件索引{len(inventory)}项；156基因完整历史表、262关系完整历史表、标准统计长表及逐基因解释已生成。','', '|分析|计划行|可评估|P<0.05|q<0.05|','|---|---:|---:|---:|---:|']
for x in counts:lines.append(f"|{x['analysis_type']}|{x['planned_rows']}|{x['evaluable']}|{x['P_lt005']}|{x['q_lt005']}|")
lines+=['','## 新手解释','先看candidate_genes156_reader.tsv，再按关系表核对具体代谢物。统计长表保留来源run与原P/q；不同分析行不代表独立证据。单细胞来源仅描述RNA在哪些细胞背景表达。','','## 限制/反证','见INTERPRETATION_CN.md；LPCAT4/CHKB、组织冲突、FUSCC受阻、Oslo连接及来源标签限制全部保留。仅对已发布文件建立索引，不宣称服务器所有工作均已发布。','','## 当前决定','当前池150＋历史6全部保留；RNA按P<0.05展示但不作为单细胞入场门槛；不排名或宣称新靶点。','','## 下一步','用本版与导师讨论有限问题；不自动启动机制分析。','','## 复现命令','在同一基准结果仓库执行 `python code/brca_sop_record_v1.py`。脚本读取当前HEAD作来源快照，因此新提交上重跑会生成新的来源指针；不能覆盖已冻结交付，正式重跑需修改RUN_ID。','','## 六步及历史结果入口']
for x in stage:lines.append(f"- {x['workflow_step']} {x['title_cn']}：[原报告]({x['report_url']})；{x['status']}。{x['limitation_cn']}")
txt('README_CN.md','\n'.join(lines))
txt('.gitattributes','* -text')
write('checksums.tsv',[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(O.iterdir()) if p.is_file() and p.name!='checksums.tsv'])
registry=R/'coordination/stages/BRCA.tsv'
old=read('coordination/stages/BRCA.tsv')
if not any(x['run_id']==RUN for x in old):
 old.append(dict(cancer='BRCA',stage_id='07_INTEGRATION',run_id=RUN,analysis_version='sop_record_v1',status='DONE',scope='564 historical artifacts indexed;156 genes262 relations preserved;1310 statistics rows;interpretations recorded;no new tests',result_path=O.relative_to(R).as_posix(),code_path='code/brca_sop_record_v1.py',git_branch='analysis/brca-functional-review-20260919',reason='Aggregate recording only;unresolved scientific gaps retained;shared SOP branch not merged',next_action='Use full reader and evidence limits for discussion;single-cell stops at source'))
 old.sort(key=lambda x:(x['stage_id'],x['run_id']))
 with registry.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(old[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(old)
print(json.dumps(dict(output=str(O),historical_files=len(inventory),statistics_rows=len(stats),summary=counts),ensure_ascii=False))
