# PDAC SOP v3：03_PATIENT

## 本轮问题

按 GitHub 规范 6dd0a73fb77a392e428254ad265b8710ce967c73 整理 PDAC，单细胞终点仅为基因表达来源。

## 输入与范围

307项代谢物，11个作者明确配对，21个作者明确且互不重复编号的肿瘤单位。配对P<0.05工作池51项；357直接关系、358条件关系。当前DIRECT池250基因，历史并集687基因。

## 实际结果

本阶段状态 **PARTIAL**。Server165:22 timeout;local reuse and mapping completed;matrix-dependent outputs not run。完整计数见summary.json。
新BH只在预定集合的可评估P上计算，旧q保存在q_original。主代谢物及关联效应/P/CI沿用同输入历史统计，不称新增验证。低于8对的新P/q留NA，旧P另存。

## 新手解释

当前基因池只包含直接关系基因；条件与历史基因没有删除，进入历史并集。RNA显著性和关联显著性不作为进入单细胞的门槛。

## 限制/反证

服务器连接超时，新配对t检验和逐细胞log1p10k汇总未运行。旧Wilcoxon结果和旧供者合并计数结果不能换名冒充新方法结果。缺项以ACCESS_BLOCKED保留；新RNA P<0.05空表是尚未执行，不是零个显著基因。687身份采用reviewed human UniProt稳定accession，矩阵按规范将另保存实际符号/稳定ID连接。独立患者临床身份和同一分装层级未重新验证。

## 当前决定

保留旧结果，按新规范分开证据维度，不沿用A/B/C/D作为总分或硬筛选。

## 下一步

服务器连接恢复后运行冻结的PDAC SOP v3，补RNA、原有值边际覆盖、三队列来源和必备点图/全基因热图。停止在来源描述，不追加机制任务。

## 复现命令

`python code/pdac/prepare_sop_v3.py`；模板来自固定GitHub提交，公共模板未改。服务器计算使用新的独占运行目录，患者与细胞级数据不进入本机或GitHub。
