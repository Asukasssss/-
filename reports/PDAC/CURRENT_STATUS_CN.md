# PDAC 当前进度：按 GitHub SOP v3 执行，服务器计算待恢复

## 本轮问题

用户要求按GitHub的CAMP发现至单细胞来源规范执行PDAC。本轮采用[固定提交规范](https://github.com/Asukasssss/-/blob/6dd0a73fb77a392e428254ad265b8710ce967c73/docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md)，单细胞止于表达来源，不新增机制分析。

## 输入与范围

作者明确配对11对，关联使用21个作者编号互不重复的肿瘤单位；6个缺编号肿瘤不默认独立。全307代谢物主P<0.05工作池51项，配对q<0.05仍为0。现有357直接与358条件关系复用；按SOP当前DIRECT池为250基因，历史并集为687基因。此前538是配对51项直接加条件池，不是新规范当前DIRECT池。

## 实际结果

|部分|状态|本轮实际完成|
|---|---|---|
|全量代谢物与原有值敏感性|PARTIAL|保留原效应/P/CI，最低8对与新可评估BH；307主可评估，193敏感性可评估。缺每侧可用性边际比例|
|映射与池定义|DONE（本批整理）|51项去向、357直接、358条件；250当前基因、687历史并集；稳定UniProtID唯一|
|肿瘤内部关联|DONE（统计复用批）|四个直接/条件×主/可用性集合重新BH，效应/P和种子保持不变|
|新配对RNA|ACCESS_BLOCKED|已冻结配对t、4000整对bootstrap；未实际执行。旧Wilcoxon全表及180项P<0.05展示另存历史|
|新单细胞来源|ACCESS_BLOCKED|三队列×687基因拟按逐细胞log1p10k、供者等权、80%描述标签运行；没有把旧pseudobulk结果改名|
|两张最终比较表|PARTIAL|715关系和687基因保留，新增RNA/来源明确待算；没有A/B/C/D总分类|

GABA—ABAT相同输入关联rho仍为0.84866；新直接主q=0.0354，可用性rho=0.84412、q=0.0342。改变来自357→354、357→342的BH有效检验数，旧q=0.0357保留；这不是新验证。GABA配对q=0.2578仍未通过FDR。

## 新手解释

当前250基因是直接关系集合；其余条件和历史基因继续保留，不按P值删去。新RNA展示表目前没有结果行，意思是未运行，不是“0项P<0.05”。本地五项方法单元测试及六个BH族的独立复核通过，不等于服务器流程回归。

## 限制/反证

2026-09-22访问172.22.148.165:22多次超时，未到认证阶段。原矩阵和真实连接表按规定留服务器，无法在本机补算。临床身份、分装层级、临床协变量和跨研究患者独立性仍有未确认项。旧303/42、旧362关联及v2配对RNA/来源不覆盖。

## 当前决定

已完成部分逐批上传；新规范整流程尚未完成。保持PR #2草稿、未合并。外部代谢关系验证及功能深入不属于本次完成前置。

## 下一步

恢复可访问服务器的内网/VPN后，执行冻结的新目录代码，补齐新配对RNA、单细胞来源、必备供者等权点图/全候选热图，再更新整合表。现成UMAP在SOP为可选，本轮不重新嵌入。

## 复现命令

本地：`python code/pdac/prepare_sop_v3.py`，`python code/pdac/verify_sop_v3_local.py`，`python code/pdac/finalize_sop_v3_preflight.py`。
服务器：新的独占PDAC/B目录执行`run_sop_v3_server.py --mode internal`或`--mode source`，各自目录传入`--code-commit`。仅代码和注释包可传输，患者/细胞级数据不得回传。

交付入口：

- [逐表完成状态](../../results/PDAC/07_INTEGRATION/20260922T112900Z_sop_v3_integration/SOP_DELIVERY_STATUS.tsv)
- [关系级表（部分完成）](../../results/PDAC/07_INTEGRATION/20260922T112900Z_sop_v3_integration/candidate_relations_integrated.tsv)
- [基因级表（部分完成）](../../results/PDAC/07_INTEGRATION/20260922T112900Z_sop_v3_integration/candidate_genes_integrated.tsv)
- [当前基因池](../../results/PDAC/02_MAPPING/20260922T112600Z_sop_v3_mapping/gene_pool_current.tsv)
- [原v2结果解释，保留历史口径](../../results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2/FINDINGS_CN.md)
