# PRAD 输入、真实配对与冻结结果审计

问题：能否按BRCA主线进行PRAD癌种内分析？

输入：CAMP v0.3.4原始PRAD.xlsx的sampleinfo、处理矩阵、MasterMapping与作者RNA矩阵；全部源数据留server165。

实际结果：137个标本，91肿瘤病例、46正常病例；作者case明确连接43对。组织标签全部一致，每组织内case唯一，配对均同批次。另一个CAPT字段有7对不一致，不能隐去，预设仅CAPT一致的36对敏感性。

原矩阵476个特征；同时通过作者两组织过滤的361项为新分析全族。旧癌种冻结表351项，原q<0.05为34。10项具有碰撞化学键，保留不同峰名与feature_id，不合并或遗漏。

新手解释：43对指同一作者病例有两种组织，不是把两列按顺序配起来。新配对分析尚未在本审计批执行，旧非配对P/q完整保留。

限制：依据作者case字段，不是基因型身份验证；CAPT差异含义未明确。data中有限值只称作者可用值，不能宣称原始检出掩码。

当前决定：冻结analysis_spec，按作者case进行配对主分析，并执行可用值与CAPT一致配对敏感性。

下一步：全361特征配对分析、P<0.05工作池、直接生化映射；相关/RNA/单细胞保留各自状态。

复现：python prad_discovery_v1.py --mode audit --out NEW_SERVER_PRAD_RUN --code-commit COMMIT；然后 --mode paired --out SAME_RUN。
