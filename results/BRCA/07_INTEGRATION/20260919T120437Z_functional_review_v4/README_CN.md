# BRCA v4完整候选整合索引

## 本轮问题
将本轮功能证据、反证和研究工作安排回接全部117个基因，不删选其他107个。

## 输入与范围
原BRCA快照：227245e3171ac3b8530835624711c1a0a69363cf；174关系、171可评估。本次不新增患者统计。

## 实际结果
完整结果仅保存于`results/BRCA/05_FUNCTION/20260919T120437Z_functional_review_v4/`，本目录是07_INTEGRATION索引，不复制为独立证据。主要文件：`candidate_comparison_117_functional_v4.tsv`、`relations_174_annotated_v4.tsv`、`functional_studies_v4.tsv`、`priority_work_queue_v4.tsv`。

## 新手解释
117全保留；24优先核查、57保留、36暂挂是工作预算安排，不是治疗靶点认证。

## 限制／反证
功能取证未完成所有全文质量审计；患者独立性、纯度/批次及代谢介导未解决。DepMap为空，CPTAC混杂不加分。GitHub写入403，本批未上传。

## 当前决定
当前有界证据下的完整比较已完成；各候选待补事项见主表。本阶段完成不表示整个BRCA完成。

## 下一步
不再围绕原相关性循环筛选；优先补每个候选最关键的功能、模型或公开数据证据。

## 复现命令
`python code/brca_functional_review_v4.py --repo-root . --register-stage`
