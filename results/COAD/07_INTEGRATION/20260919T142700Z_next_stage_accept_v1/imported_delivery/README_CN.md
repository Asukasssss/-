# COAD 下一阶段：9基因任务细化、跨癌候选对照与协变量方案

日期：2026-09-19

## 本轮完成了什么

读取固定版本的COAD已接收报告、35基因功能表、原患者分析配置，以及BRCA较新的ER调整/可用性敏感性结果。本包完成既有证据的下一步任务细化、候选身份集合匹配和3条同名代谢物关系的描述性对照。

本轮没有读取server165患者矩阵，没有运行新的患者统计，没有完成新的系统文献复核，也没有向GitHub写入或合并PR。协变量文件均标为DRAFT / NOT RUN，不能将它们当作完成结果或已经注册的方案。

## 1. 实际集合比较

比较范围是 **COAD已核查35基因 vs BRCA完整117基因**，不是COAD全部458个计划基因与BRCA的全候选比较。两边的筛选阶段不同，不做交集富集显著性检验，也不据不在交集中判定癌种特异性。

共同基因：**BCAT2、NNMT、SLC6A6、UPP2**。
COAD原9个优先基因中，进入该交集的是：**BCAT2、SLC6A6**。

|共同基因|可以对应的同名代谢物|本轮解释|
|---|---|---|
|BCAT2|异亮氨酸|BRCA不只有谷氨酸关系，也包含异亮氨酸；两癌种该关系方向均为负|
|NNMT|1-MNA|两癌种该关系方向均为正；仍保留参照候选定位，不自动升为最高优先|
|SLC6A6|牛磺酸|两癌种该关系方向均为正；不能据此认定共同耐药机制|
|UPP2|无；COAD为尿苷，BRCA为尿嘧啶|共同基因，不是相同代谢物关系；不把符号差异称为反向复现|

三个同名关系是已读表中相同代谢物名称的对应，不构成新的化学身份、区室、通量或介导验证。COAD的SAM与BRCA的SAH也不合并为同一种代谢物。

## 2. 统计对照

|关系|COAD未调整rho（展示值）|COAD q（展示值）|BRCA ER调整rho|BRCA ER主分析q|BRCA ER+可用性q|
|---|---:|---:|---:|---:|---:|
|BCAT2—异亮氨酸|-0.4051|0.6878|-0.3355|0.1005|0.1148|
|NNMT—1-MNA|0.4071|0.6878|0.3844|0.0534|0.0534|
|SLC6A6—牛磺酸|0.3570|0.7335|0.1672|0.4351|0.4828|

COAD数字来自已接收功能表的四位小数展示字段，不是重新读取其原高精度患者统计。BRCA提取高精度数值保存在TSV，工作簿显示四位小数。原统计不做修改。

COAD未调整Spearman与BRCA ER调整偏秩相关不是同口径模型。这里只能描述符号和数值背景，不能据此给出跨癌效应差异检验、合并效应或外部复现结论。三条均没有达到“两癌种均通过BH校正”的证据标准。BRCA的ER+可用性仍来自同一批标本；这三条关系在该轮输入相同，复用统计，不增加独立证据次数。

BRCA更新后的两轮ER分析均获支持的关系涉及GPCPD1和GPI，不能用历史11条未调整显著关系直接充当当前结论。它们没有出现在本次COAD名义候选35基因中，并不说明其在COAD全部目录中不存在。

## 3. 9个优先基因怎样推进

完整任务见工作簿“9基因任务”及COAD_priority9_action_matrix.tsv。原9/21/5安排不变。

BCAT2、SLC6A6有本轮跨癌关系对照入口；BCAT2保留饮食/circRNA背景及正常饮食单独敲低阴性条件，SLC6A6优先寻找治疗或细胞状态信息，而不是重复已经报道的通用耐药故事。

HDC、GSTA4先解决细胞来源。公开单细胞或空间表达只能提供定位/细胞组成线索，不能替代本CAMP样本的纯度或功能验证；按患者汇总，不将每个细胞视作独立患者。GSTA4宿主构成性全身缺失不改写成髓系特异缺失，HDC正相关不推出应抑制HDC。

AQP9、KMT5A先按模型、干预、处理条件、底物或位点和终点整理方向差异；PRMT7先确认能否测到具体剪接变体；SLC38A3区分HSP70相关功能与天冬酰胺特异运输/救援；UCKL1保留UMP/CMP未救援和尿苷结合区缺失实验的反证，不能概括为尿苷消耗介导抗肿瘤作用。

这些是对已接收证据的执行问题细化，不是本轮独立重审所有论文后新增的阳性功能结论。

## 4. 协变量分析：已形成方案，尚未运行

主要人群沿用33名明确个体的I-IV期肿瘤样本。原37例扩展敏感性含3例腺瘤和1例stage0，不混入主分析。

原974条关系目录全部保留；新调整检验应覆盖原674条预定可纳入关系，逐模型报告实际可估计数和缺失原因。652是此前未调整模型的可计算数，不承诺调整后仍是652。

建议模型：M1=年龄+性别；M2=M1+分期（分类）；M3=M1+取材部位（原始类别）。M2与M3分别作为敏感性分析，不根据哪个P更小选择结论。M0保持原结果，每个模型另计算M0_CC：与该调整模型完全相同病例子集上的未调整参照。

先只看协变量的频数、缺失与模型设计矩阵，不按候选显著性决定分组。编码、稀疏类别及可估计性规则在真实关系分析前锁定；秩亏或样本支持不足时明确标记，而不是暗中删变量。无新插补，不从其他癌种借用ER，也不生成伪纯度或伪技术批次。

每个模型在所有可估计关系内单独BH校正，保留全674个关系ID，包括不可估计项。报告N、效应、区间、P/q、与M0_CC的变化及影响诊断，不能只报告过阈值的39条或9基因。

偏秩相关的置换必须兼容协变量结构，不能机械复制未调整随机打乱流程。具体重采样实现、次数、诊断和验证尚待锁定与执行，本包不声称提供已验证的患者统计脚本。

## 5. 文件说明

- COAD_next_stage_workbook.xlsx：7个工作表，含公式计算的集合统计。
- COAD_priority9_action_matrix.tsv：9基因任务表。
- BRCA_COAD_shared4_gene_comparison.tsv：4共同基因。
- BRCA_COAD_same_named_metabolite3_comparison.tsv：3同名关系及BRCA高精度提取值。
- COAD_covariate_plan_NOT_RUN.tsv / covariate_spec_DRAFT_NOT_RUN.json：待执行方案。
- COAD_35_gene_identities_extracted.tsv、BRCA_117_gene_identities_extracted.tsv、BRCA_174_relation_identities_extracted.tsv：从连接器读取内容整理的身份字段，不是源TSV的逐字节副本。
- reproduce_overlap.py：仅复现集合比较，不进行患者统计。
- SERVER165_NEXT_TASK_CN.md：下一轮服务器执行交接。
- provenance.json、validation.json、MANIFEST_SHA256.tsv：固定来源、范围与新文件校验。

复现集合比较：在解压目录运行 `python reproduce_overlap.py`。标准库即可。验证计数与交集，不验证原患者数据或原文结论。

## 来源版本

COAD提交：6f1d28a5607933e34886191b820fd16a4a69b832
BRCA提交：227245e3171ac3b8530835624711c1a0a69363cf

COAD功能表：https://github.com/Asukasssss/-/blob/6f1d28a5607933e34886191b820fd16a4a69b832/results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1/gene_function_review.tsv
BRCA完整关系表：https://github.com/Asukasssss/-/blob/227245e3171ac3b8530835624711c1a0a69363cf/results/BRCA/20260919_ER_availability_v3/er_availability_all_174.tsv
COAD原分析配置：https://github.com/Asukasssss/-/blob/6f1d28a5607933e34886191b820fd16a4a69b832/config/coad_patient_v1.json
BRCA更新报告：https://github.com/Asukasssss/-/blob/227245e3171ac3b8530835624711c1a0a69363cf/results/BRCA/20260919_ER_availability_v3/README_CN.md
COAD已接收报告：https://github.com/Asukasssss/-/blob/6f1d28a5607933e34886191b820fd16a4a69b832/results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1/README_CN.md

方法参考：Winkler等，Permutation inference for the general linear model (2014)，PMID 24530839，https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/。用于协变量置换方案的进一步实现与核查，不表示本包完成了该方法的统计验证。
