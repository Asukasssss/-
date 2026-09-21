"""Build a current BRCA aggregate evidence overview; no new hypothesis tests."""
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
RUN='20260921T114912Z_current_overview_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION'/RUN
OUT.mkdir(parents=True,exist_ok=True)
MET=ROOT/'results/BRCA/04_ROBUSTNESS/20260921T112611Z_metabolite318_paired_v1'
PAT=ROOT/'results/BRCA/04_ROBUSTNESS/20260921T103911Z_camp_pair_sensitivity_v1'
EXT=ROOT/'results/BRCA/06_EXTERNAL/20260920T090922Z_all174_external117_resources_v1'
TANG=ROOT/'results/BRCA/06_EXTERNAL/20260921T072503Z_oslo_tang_external_v1'
sources=[MET/'comparison318.tsv',MET/'comparison_summary.json',PAT/'paired_RNA117.tsv',PAT/'family_before_after.tsv',PAT/'all117_comparison_identity_appended.tsv',EXT/'external_validation.json',TANG/'analysis_summary.json']
m=pd.read_csv(sources[0],sep='\t');e=pd.read_csv(sources[2],sep='\t');master=pd.read_csv(sources[4],sep='\t')
assert len(m)==318 and len(e)==117 and len(master)==117
sources += [TANG/'GSE42568_summary.json',ROOT/'results/BRCA/05_FUNCTION/20260921T034354Z_asns117_programs_v1/README_CN.md',ROOT/'results/BRCA/05_FUNCTION/20260921T034354Z_asns117_programs_v1/mapping_validation.json',ROOT/'results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/integration_validation.json',ROOT/'results/BRCA/03_PATIENT/20260921T102429Z_camp_sample_identity_v1/validation.json']
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
save(m[m.new_significant&~m.old_significant],OUT/'new11_paired_metabolites_mapping_pending.tsv')
save(m[~m.new_significant&m.old_significant],OUT/'old21_no_longer_paired_significant.tsv')
rows=[
['ASNS / GLS—谷氨酰胺','谷氨酰胺39/45对降低；45对原数据完整','ASNS配对RNA升高；GLS无明确组间RNA差异','FUSCC负关联及条件关联；Asns干预RNA已分析','未证明谷氨酰胺通量/介导或ASNS—GLS补偿','保留共同专题；不以RNA差异决定两者去留'],
['SLC6A8—肌酸','肌酸36/45对降低；45对原数据完整','SLC6A8配对RNA差异未通过校正','FUSCC正关联','丰度变化不等于运输活性；未证明抑制方向','保留功能/适用背景方向'],
['GPCPD1—GPC','GPC31/45对降低；45对原数据完整','CAMP调整负关联效应稳定；配对RNA差异未通过校正','既有功能文献；外部关联未形成稳定支持','主模型q0.0513，次模型q0.0171；并非两模型均显著','保留CAMP内部线索，明确外部限制'],
['GPI—G6P/F6P','G6P39/45对升高；F6P37升7降1平；可用性子集40/19对','G6P调整效应稳定；GPI配对RNA无明确差异','Tang可评估但未支持正关联','代谢物配对变化不等于GPI调控机制；F6P高缺失','保留内部线索和外部未支持结果'],
['NNMT—1-MNA','1-MNA32/45对升高；原数据可用39对','NNMT配对RNA降低；组成调整后关联支持有限','FUSCC正关联','组间RNA较低与肿瘤内部正相关回答不同问题','保留细胞背景和适用模型问题'],
['MDH1—苹果酸','苹果酸35/45对升高；原数据可用41对','MDH1配对RNA降低','Tang较强正关联；FUSCC接近零；GSE42568中RNA较低','跨队列不一致；不借用MDH2功能证据','限定为值得核查的新外部线索'],
['PNP—鸟嘌呤/次黄嘌呤','鸟嘌呤40/45对升高；次黄嘌呤37/45对升高','配对RNA未通过校正；CAMP鸟嘌呤正关联','Tang次黄嘌呤q0.0489；鸟嘌呤q0.07335','两关系不是两次独立队列验证；Tang仅20例','保留关系级比较，不能合并成一条已验证机制'],
['PYCR1','既有CAMP生化映射保留；不以其他代谢物替代','配对RNA升高37/45；已有单细胞细胞背景','已有CAF功能文献和有限终点整理','重整理同一论文不是独立验证；不能凭表达最高限定唯一作用细胞','保留限定专题，不提升为全项目唯一主线']]
cols=['候选方向','配对代谢物层','CAMP基因层','外部或功能层','限制与反证','当前安排']
save(pd.DataFrame(rows,columns=cols),OUT/'candidate_directions_CN.tsv')
text='''# BRCA当前结果总览（2026-09-21）

## 本轮问题
将发现、直接映射、样本身份、最新配对代谢物/RNA、患者关联、外部队列、单细胞及功能证据放在同一张研究路线图中。此轮仅整理已完成结果，新增患者统计为0，不重新筛基因。

## 输入与范围
原冻结BRCA为318条代谢物效应、186条q<0.05。首轮73个代谢特征建立174条直接关系、117个基因；其他映射状态及未知特征保留。当前117基因池来源于旧186条，不是新176条配对显著代谢物的完整映射。

## 实际结果

|阶段|当前结果|含义与边界|
|---|---|---|
|样本身份|108个跨组学LHC一致；63个作者病例记录；45对标签一致组织；1个标签冲突|依据作者明确对应表，不按编号猜配；冲突标本排除但未判定哪方标签错误|
|配对代谢物|318项均有45对作者处理值；176项q<0.05|旧186条中165保留、21本轮未过线；另11条新通过。原统计仍冻结|
|缺失值敏感性|288项可评估、121项显著；93项两新家族均显著且均值方向一致|57项在可评估子集中均值方向不同；不能把176项都称稳健，也不强制只保留93项|
|配对RNA|117保留、116可评估、60项q<0.05：39较高、21较低|每项45对；CHKB缺测。此处是RNA差异，不是代谢物差异|
|CAMP肿瘤内关联|排除冲突后60病例；174关系中171可评估，未调整主分析12条q<0.05|ER处理值模型7条、ER可用子集8条；组成主模型0条、次模型2条。不同家族不能按q直接排名|
|FUSCC外部关联|258作者患者ID交集；133/174可评估；主模型16条q<0.05|不同关系分别记录；不等于16条都同向复现CAMP；原BH174范围保留|
|Tang外部关联|20例连接；163/174可评估；2条q<0.05|MDH1—苹果酸、PNP—次黄嘌呤；小队列，跨队列一致性有限|
|GSE42568 RNA背景|112/117可评估，57条差异q<0.05|只补RNA层，不能替代代谢关联验证|
|单细胞背景|全117已分析；43个达到既定两研究细胞来源描述稳定规则|描述表达背景，不是43个靶点通过功能验证；旧82与43口径不同|
|功能证据|全117已有文献检索/证据整理，保留反证与模型限制|不宣称全篇全文审计，也不把历史24/57/36工作安排当最新疗效排名|
|Asns干预|四个比较均捕捉到目标RNA降低；其他三个节点未检出明确联合校正补偿上调|小鼠、选择后再培养模型；没有测到本批次谷氨酰胺通量|
|Asns全候选回接|117均保留；111人鼠一对一映射；原发瘤来源23基因两构建同向显著，含Asns|其余22是响应基因，不是各自获得功能验证；肺来源仅Asns达同一规则|

### 最直接的患者配对结果

|代谢物|肿瘤较高/较低/持平|两侧原数据可用对数|
|---|---|---:|
|谷氨酰胺|6 / 39 / 0|45|
|谷氨酸|39 / 6 / 0|45|
|肌酸|9 / 36 / 0|45|
|GPC|14 / 31 / 0|45|
|G6P|39 / 6 / 0|40|
|F6P|37 / 7 / 1|19|
|天冬酰胺|23 / 22 / 0|45|

除天冬酰胺外，上述重点代谢物均通过本轮主分析及可用性敏感性校正；各家族q分列。数量来自作者处理值，F6P等必须带上缺失限制。不能从ASNS升高推导天冬酰胺必然升高，或从谷氨酰胺降低推导ASNS/GLS反应通量增加。

### 重点候选当前判断
详见candidate_directions_CN.tsv，列出8个讨论方向，不是新的前8名。GPCPD1—GPC组成主模型rho=-0.466/q0.0513，次模型rho=-0.529/q0.0171；GPI—G6P分别+0.447/q0.0513和+0.459/q0.0342。与旧效应接近，不能把阈值跨越写成关联消失。

ASNS与PYCR1的RNA分别36/45、37/45对升高；NNMT与MDH1分别30/45、26/45对降低。方向一致性与均值差异显著性分开解释。GLS、GPCPD1、GPI、SLC6A8、PNP的配对RNA未通过校正，不取消其其他证据。

## 新手解释
目前可以回答三层问题：代谢物在同患者肿瘤/正常中怎么变；哪些直接相关基因与其在肿瘤中一起变化；外部患者和干预研究是否提供支持。单细胞帮助判断在哪里研究，RNA变化不等于酶活、代谢通量或因果机制。

## 限制/反证
- 45对并非新增独立队列；不把旧非配对、配对和敏感性相加为多次验证。
- 原数据可用性子集改变病例构成；方向改变不能单独证明作者填补错误。
- 组成评分不是细胞百分比；临床治疗价值、正常组织风险、代谢介导仍未证明。
- Oslo2当前没有可追溯编号对照，真实连接0；CPTAC原肿瘤正常对照存在组织/批次混杂；DepMap数值筛选尚未完成。
- 文献已研究的功能不能作为本项目首次发现；各基因自己的干预与响应其他基因干预不能混写。

## 当前决定
保留全117及历史结果，以新配对和排除冲突版本解释当前CAMP证据。谷氨酰胺ASNS/GLS可继续作为专题；SLC6A8、GPCPD1、NNMT、MDH1、PNP、GPI、PYCR1分别保留适用支持和限制，不机械缩成两个基因。

## 下一步
1. 先回接原174关系对应的最新代谢物配对与缺失敏感性标签，避免继续展示过期发现标签。
2. 对新增11条配对显著代谢物核对现有映射与身份：其中未知X特征不猜基因；可靠新增关系另建映射版本，不混进旧q。
3. 原21条未通过配对校正的代谢物及其基因保留功能、外部支持；不自动删除。
4. 深入限于能回答具体问题的少数方向；不追加无目的全117单细胞、通路或生存流水线。

## 复现命令
`python code/brca_current_overview_v1.py`
本轮没有重算统计。关键数值读取已提交汇总表；来源和哈希见source_manifest.tsv。详细方法、完整318/117/174结果以各阶段路径为准。
'''
(OUT/'README_CN.md').write_bytes(text.encode())
(OUT/'.gitattributes').write_bytes(b'* -text\n')
(OUT/'analysis_spec.json').write_bytes(json.dumps(dict(version='current_overview_v1',new_tests=0,scope='aggregate evidence organization;not a new ranking',historical_results_modified=False),indent=2).encode())
save(pd.DataFrame([dict(path=p.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources+[Path(__file__)]]),OUT/'source_manifest.tsv')
(OUT/'validation.json').write_bytes(json.dumps(dict(metabolites=318,genes=117,master_columns=len(master.columns),new11=int((m.new_significant&~m.old_significant).sum()),old21=int((~m.new_significant&m.old_significant).sum()),RNA_significant=int(e.q_value.lt(.05).sum()),RNA_up=int((e.q_value.lt(.05)&e.effect.gt(0)).sum()),RNA_down=int((e.q_value.lt(.05)&e.effect.lt(0)).sum()),new_patient_statistics=0),indent=2).encode())
print(OUT)
