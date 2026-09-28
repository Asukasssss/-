# GSE272362 作者空间对象：LYPLA1 覆盖核查

## 本轮问题与来源
寻找可判断 PDAC LYPLA1 癌上皮富集、带病理注释及原始组织图的空间队列。使用 [原作者 GSE272362 研究](https://doi.org/10.1038/s41588-024-01914-4) 链接的 [Zenodo v2 对象](https://zenodo.org/records/22695225)，下载完整 PDAC_Updated_V2.rds 并通过指定 MD5 校验。来源尺寸、SHA-256 和版本见 target_coverage.json。

## 实际结果
对象为 Seurat，成功读到 Spatial 表达层：17,893 个特征 × 91,496 个点；其行名未匹配 LYPLA1 或 ENSG00000120992。没有开展 LYPLA1 表达差异分析，当前提供对象的该层无法支持目标比较。

随后检查 SCT 等其他层时，因服务器没有完整 Seurat 包而中止。inspection_log.txt 保留原始错误。**这是部分对象核查，不是全对象验证通过；原始 counts 槽的独立核对未完成，其他分析层和患者映射也未完成核查。**

## 新手解释与限制
有大量空间点和癌区注释的论文，不意味着其可用表达对象一定包含所需基因。这里是已读到的 Spatial 表达层没有目标行，不能写成 LYPLA1 零表达，也不能证明整项研究的原始测量完全没有该基因。

## 决定与下一步
本轮 LYPLA1 癌上皮比较为 NOT_EVALUABLE，对象审查范围为 PARTIAL。保留该候选及缺项，不以补零、替换同家族基因或模型预测表达代替测量。本轮实际计算与配准图由 GSE278687 五张新鲜冷冻切片提供；上皮参考关联仍与恶性身份分开解释。

## 复现
code/pdac/inspect_mets_spatial.R 为实际受测检查代码；code/pdac/record_mets_coverage.py 只记录成功读到的覆盖结果、后续失败和来源校验，不声称恢复了未执行部分。参数见 analysis_spec.json。源对象留 server165。
