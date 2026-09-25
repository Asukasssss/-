# ccRCC3治疗资料细化敏感性

问题：原始MetaData_M4的NO标签是否都代表未治疗？

输入：原研究master_mapping，以CLIENT_IDENTIFIER一对一连接，并用RNA ID交叉核实114标本。

实际结果：67肿瘤均有clear cell病理；原NO组有1位患者在THERAPY_EXPOSURE记为Cabozantinib。采用两个作者字段的组合分层重算完整关系敏感性，计数见validation.json。

新手解释：更细治疗分层仅用于敏感性，不是治疗效果检验。

限制：治疗历史不完全等于采样时用药；保留两个字段，不推断未记录时间点。

当前决定：主关联、可用值敏感性和RNA原数值原样复用；仅以本版替代旧粗分层敏感性。

下一步：整合。

复现：python ccrcc3_treatment_sensitivity_v2.py --out NEW_RUN --source ORIGINAL_CCRCC3_RUN --code-commit COMMIT。
