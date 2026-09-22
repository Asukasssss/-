# BRCA按统一规范归档与解释 v1

## 本轮问题
将已发布BRCA结果按六步发现至单细胞来源规范记录，保留全部历史并追加可追溯解释。

## 输入与范围
基准提交：`5f78d08a696d5fac3d95c697b0eb05358cfa1868`。只读仓库公开汇总；不读服务器矩阵，不新增P/q。公共SOP为6dd0a73分支版本，未合并main。

## 实际结果
历史文件索引564项；156基因完整历史表、262关系完整历史表、标准统计长表及逐基因解释已生成。

|分析|计划行|可评估|P<0.05|q<0.05|
|---|---:|---:|---:|---:|
|available_Spearman262|262|255|56|9|
|processed_Spearman262|262|255|69|14|
|paired_RNA150|150|148|85|79|
|paired_author_available|318|288|148|121|
|paired_processed|318|318|190|176|

## 新手解释
先看candidate_genes156_reader.tsv，再按关系表核对具体代谢物。统计长表保留来源run与原P/q；不同分析行不代表独立证据。单细胞来源仅描述RNA在哪些细胞背景表达。

## 限制/反证
见INTERPRETATION_CN.md；LPCAT4/CHKB、组织冲突、FUSCC受阻、Oslo连接及来源标签限制全部保留。仅对已发布文件建立索引，不宣称服务器所有工作均已发布。

## 当前决定
当前池150＋历史6全部保留；RNA按P<0.05展示但不作为单细胞入场门槛；不排名或宣称新靶点。

## 下一步
用本版与导师讨论有限问题；不自动启动机制分析。

## 复现命令
在同一基准结果仓库执行 `python code/brca_sop_record_v1.py`。脚本读取当前HEAD作来源快照，因此新提交上重跑会生成新的来源指针；不能覆盖已冻结交付，正式重跑需修改RUN_ID。

## 六步及历史结果入口
- 01 身份和配对：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/03_PATIENT/20260921T102429Z_camp_sample_identity_v1/README_CN.md)；NEEDS_REVIEW。身份追溯完成，1个组织标签冲突仍未解决
- 02 配对代谢物：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/04_ROBUSTNESS/20260921T112611Z_metabolite318_paired_v1/README_CN.md)；DONE。作者data可用值不是独立核验仪器检测掩码
- 03 直接映射：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/README_CN.md)；DONE。映射不是功能或患者统计
- 04 患者相关和配对RNA：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/03_PATIENT/20260921T150000Z_mapping262_patient_v1/README_CN.md)；PARTIAL。关系255/262、RNA148/150；缺失或身份待核保留
- 05a 旧117单细胞来源：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/README_CN.md)；DONE。注释来源与队列覆盖见原报告
- 05b 旧117来源稳定性：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/README_CN.md)；DONE。描述性稳定性不是患者百分比
- 05c 新增39来源：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/README_CN.md)；DONE。LPCAT4单细胞身份不能解除微阵列待核
- 05d 全156分型背景：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/README_CN.md)；DONE。临床分型不是上皮细胞状态；Pal无可靠分型连接
- 05e 全156表达UMAP：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/README_CN.md)；DONE。UMAP不检验细胞类型富集
- 05f 细胞类型对照UMAP：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/README_CN.md)；DONE。图示不是机制证据
- 05g 上皮专用UMAP：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/README_CN.md)；DONE。可选展示；不称新聚类或患者亚型发现
- 06 整合：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/07_INTEGRATION/20260921T154000Z_mapping262_integration_v1/README_CN.md)；DONE。外部与功能属于已存附加证据，不是来源主线门槛
- 附加 功能文献117：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/05_FUNCTION/20260919T120437Z_functional_review_v4/README_CN.md)；DONE。不等于所有原文全文逐图审计
- 附加 功能文献新增39：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/05_FUNCTION/20260921T152000Z_new39_literature_v1/README_CN.md)；DONE。遵循模型适用范围，不据文献命中证明功能
- 附加 外部262关系：[原报告](https://github.com/Asukasssss/-/blob/5f78d08a696d5fac3d95c697b0eb05358cfa1868/results/BRCA/06_EXTERNAL/20260921T151000Z_mapping262_external_v1/README_CN.md)；PARTIAL。FUSCC64条主关系RNA下载受阻；Oslo连接未解决
