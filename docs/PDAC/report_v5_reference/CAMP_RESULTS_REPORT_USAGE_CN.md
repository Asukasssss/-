# 配置驱动的结果展示程序

本程序读取指定Git提交中的公开汇总，生成真实PDF、Markdown、Excel、独立图件、图源TSV与ZIP。不运行患者模型、不计算P/q、不重新聚类或生成UMAP。默认拒绝覆盖非空输出目录。

```bash
python code/camp_results_report.py --config configs/BRCA_report_v1.yaml --validate-only
python code/camp_results_report.py --config configs/BRCA_report_v1.yaml --out results/BRCA/07_INTEGRATION/NEW_RUN_ID
python -m unittest discover -s tests -p test_camp_results_report.py -v
```

BRCA输入完整提交固定为8663004f8d98acaf1368badc7436c84b059e4462；分析规范参照6dd0a73。展示分支的HEAD不是数据输入。主报告的来源纠正直接读取该基线的原汇总，不借后续提交替换数值。

## 可复用的部分

1. **输入**：`files`为逻辑表名到仓库路径；`input_commit`必须为完整不可变SHA。文件以Git对象读取，不采用当前工作目录中同名表。
2. **薄适配**：统一统计列为gene、cohort、metabolite_key、analysis_type、effect、effect_type、p_value、q_value、status、n；异名通过field_aliases映射。关系表保留完整关系ID，阶段特有字段参照BRCA配置与示例TSV；新的结构需显式适配，不猜含义。
3. **配置**：癌种、研究、实际设计、源表、家族名称、类别顺序/中文标签、分型、示例关系、预期计数、正文解释均在配置。源表中未列入首选排序的细胞类别自动追加，不静默丢弃。
4. **分队列**：代谢物、相关和RNA按cohort分开出图，不汇成新检验。示例关系ID必须能唯一定位到具体队列关系；多队列同名关系需要适配为带队列的完整ID，不可跨队列复用一条ID。
5. **可选模块**：没有正常RNA时省略rna源表并关闭modules.rna；没有分型时省略subtypes/subtype_profiles并关闭modules.subtypes；只有一套单细胞时仅列该研究。没有可用值敏感性时不提供相应analysis_type，记录缺项，不编造数据。
6. **展示边界**：新癌种必须同步改narratives与design，不能复制BRCA结论。预期计数仅用于验收，不决定图上结果。示例在display_selection.tsv固定，并非显著性排序。
7. **只读验证**：validate-only在载入/检查后直接退出，不创建指定输出目录。单独合成测试覆盖不同数量、单研究、缺RNA、缺分型、缺敏感性、不可解释来源、真零、错误关系和拒绝覆盖。

## 输出

- 主报告：问题—图表—结果—解释—限制，精选仅发生在展示层。
- 全量附录与Excel：包括历史列、来源状态、完整关系、全代谢物与RNA、单细胞供者汇总和已有分型。
- 新图：PNG、PDF、SVG、绘图源TSV、中文图注。全基因热图自动分页，行缩放标注只用于同基因分布。
- 原图复用：只复制原字节并给页码/尺度索引，不从栅格图恢复矢量，不伪造原始坐标。复用图册的源表不是逐细胞坐标。
- figure_manifest、source_manifest、校验清单及版本：能回到实际源文件和输入提交。
- ZIP：含全部实际图册。大图册已在基线公开，展示分支不重复提交，本地ZIP完整包含；提供原图GitHub链接。

## 限制

已经完成BRCA真实运行和独立合成测试，不声称完成其他癌种真实数据验证。历史统计假设和数据限制不会因排版而消失。程序不能取代新的数据身份核查。字体从本机读取，字体文件不随包分发；实际依赖版本见交付repro/requirements.txt。
