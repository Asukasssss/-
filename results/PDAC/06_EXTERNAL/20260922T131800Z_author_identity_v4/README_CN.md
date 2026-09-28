# PDAC 作者恶性身份接回 v4

## 本轮问题

执行GitHub固定规范6dd0a73fb77a392e428254ad265b8710ce967c73。只处理PDAC，单细胞止于表达来源。

## 输入与范围

同一687基因和148777细胞；只改变有依据的身份划分，不按表达结果推断恶性。

## 实际结果

GSE242230旧Malignant大类11133细胞拆为10783作者恶性（5660 Classical、5123 Basal）和350作者正常上皮。该队列全687来源统计重新计算；另两队列只将导管显示为身份未定并复用数值。

## 新手解释

不同分析的P/q分别保留。当前DIRECT池250基因，历史并集687基因；RNA差异和关联显著性均不是进入单细胞的门槛。不可评估是缺项，不是阴性。

## 限制/反证

细分Normal Epithelial与宽类Malignant冲突，优先使用作者更具体身份。另两队列论文恶性分析不能代替逐细胞身份连接。作者注释不等于本轮独立CNV验证。

## 当前决定

本批已完成范围登记DONE；不覆盖历史统计，不将表达来源升级为同一代谢关系的外部验证或因果机制。

## 下一步

全候选图与比较表已更新；另两队列恶性标签接回仍为NEEDS_REVIEW。

[全部来源图](figures/README_CN.md) · [身份取证](author_identity_evidence.tsv) · [作者标签计数](author_identity_counts.tsv)

## 复现命令

`python code/pdac/run_author_identity_v4_server.py`（服务器新目录）；本地`python code/pdac/deliver_author_identity_v4.py --mode record`。


完整来源表为18,549行，按队列无损分为3个TSV以遵守单文件5 MB限制。字段、精度、顺序保持；分片清单及原逻辑表哈希见[sc_celltype_profile_parts.tsv](sc_celltype_profile_parts.tsv)。规范中的sc_celltype_profiles.tsv指这些分片合并的逻辑表。
