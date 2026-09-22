# PRAD 全量配对代谢物发现

问题：作者case定义的肿瘤与正常配对有哪些代谢差异？

输入：审计确认的43对，作者两组织过滤交集361特征；CAPT一致36对另列敏感性。使用作者尺度，无新增填补或变换。

实际结果：

- METAB_PAIRED_AVAILABLE: {'planned': 361, 'evaluable': 344, 'p_lt005': 66, 'q_lt005': 9}
- METAB_PAIRED_CAPT_CONCORDANT: {'planned': 361, 'evaluable': 361, 'p_lt005': 57, 'q_lt005': 2}
- METAB_PAIRED_PRIMARY: {'planned': 361, 'evaluable': 361, 'p_lt005': 69, 'q_lt005': 9}

新手解释：P<0.05是探索工作池，q<0.05另标FDR支持；升降比例是配对内方向人数比例，不是显著患者比例。效应是作者log2尺度均值差，不称浓度倍数。

限制：作者case身份；CAPT不一致7对保留敏感性；可用值掩码不是已验证原始检出；同队列敏感性不是独立验证。

当前决定：用主分析P<0.05全池进入直接生化映射，不按敏感性显著性删候选。

下一步：全部工作特征映射、全关系肿瘤相关、全基因RNA与单细胞来源。

复现：python prad_discovery_v1.py --mode paired --out AUDITED_RUN。
