# PRAD全关系肿瘤关联与全候选RNA背景

问题：直接关系是否在PRAD肿瘤内共变，候选RNA是否在配对肿瘤中改变？

输入：全部161关系、101基因、91个唯一作者肿瘤case与43对RNA；CAPT一致子集36对。

实际结果：

- REL_TUMOR_AVAILABLE: {'planned': 161, 'evaluable': 161, 'p_lt005': 15, 'q_lt005': 0}
- REL_TUMOR_BATCH: {'planned': 161, 'evaluable': 161, 'p_lt005': 15, 'q_lt005': 0}
- REL_TUMOR_PRIMARY: {'planned': 161, 'evaluable': 161, 'p_lt005': 18, 'q_lt005': 0}
- RNA_PAIRED_CAPT: {'planned': 101, 'evaluable': 101, 'p_lt005': 34, 'q_lt005': 20}
- RNA_PAIRED_CURRENT: {'planned': 101, 'evaluable': 101, 'p_lt005': 35, 'q_lt005': 21}

新手解释：相关系数描述两项测量共同变化；正负相关不是合成或消耗方向。RNA变化不是酶活或通量。每个检验族单独BH，原CAMP q没有移植。

限制：同队列筛选后的探索性分析，不是独立验证；作者case无基因型核实。作者批次校正后仍给批次内置换的partial Spearman敏感性，不能将批次影响隐藏。

当前决定：全部基因继续单细胞来源，不按本表显著性删除；未测项保持不可评估。

下一步：公开单细胞数据来源与供者标签核查；整合全部候选。

复现：python prad_patient_v1.py --out AUDITED_RUN --code-commit COMMIT；读取同一服务器私有连接和固定映射输入。
