# BRCA单细胞结果：表格、图和逐基因UMAP总入口

这里集中展示全部候选的既有单细胞表达来源结果。当前150个基因＋6个历史基因均保留；新目录只做链接与清单，不新增分析。所有链接固定到已上传提交，便于另一个账号直接使用。

## 先打开这些图

- [全156基因表达图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/all156_expression_UMAP.pdf)
- [全156基因表达—细胞类型对照图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/all156_expression_celltype_pairs.pdf)
- [全156基因上皮专用三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/all156_epithelial_triptychs.pdf)
- [全156基因分型热图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/all156_subtype_expression_patterns.pdf)

PDF有全部基因；重点基因另提供PNG。页码请查对应目录的gene_page_index.tsv或gene_scales_and_pages.tsv。GitHub若不能在线预览大PDF，可下载原文件查看。

## 来源和图形的边界

Wu沿用作者注释，Pal采用既有重注释；Reed为正常背景，不能计为第三套独立肿瘤验证。供者标签不是本次新核验的独立患者。分型结果以Wu为主，Pal未接入可靠分型连接。表达图显示RNA位置，不证明代谢物来源、功能或通量。上皮细胞状态不是临床患者亚型。本次不新增差异检验、聚类、通讯或拟时序。


## 原117基因：两肿瘤研究与正常背景

[说明与方法](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/README_CN.md)

|文件|用途|
|---|---|
|[all117_cell_source.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/all117_cell_source.png)|图片|
|[analysis_spec.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/analysis_spec.json)|参数或验证|
|[candidate_comparison_117_sc_appended.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/candidate_comparison_117_sc_appended.tsv)|数值或来源表|
|[cell_coverage_sensitivity.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/cell_coverage_sensitivity.tsv)|数值或来源表|
|[celltype_labels_CN.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/celltype_labels_CN.tsv)|数值或来源表|
|[checksums.sha256](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/checksums.sha256)|校验/记录|
|[gene117_cell_source_CN.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/gene117_cell_source_CN.tsv)|数值或来源表|
|[gene117_cell_source_comparison.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/gene117_cell_source_comparison.tsv)|数值或来源表|
|[integration_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/integration_validation.json)|参数或验证|
|[key_genes_cell_source.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/key_genes_cell_source.png)|图片|
|[numeric_fixture_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/numeric_fixture_validation.json)|参数或验证|
|[Pal2021_reprocessed_celltype_profiles.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Pal2021_reprocessed_celltype_profiles.tsv)|数值或来源表|
|[Pal2021_reprocessed_feature_coverage.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Pal2021_reprocessed_feature_coverage.tsv)|数值或来源表|
|[Pal2021_reprocessed_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Pal2021_reprocessed_validation.json)|参数或验证|
|[README_CN.md](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/README_CN.md)|中文说明|
|[Reed2024_celltype_profiles.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Reed2024_celltype_profiles.tsv)|数值或来源表|
|[Reed2024_feature_coverage.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Reed2024_feature_coverage.tsv)|数值或来源表|
|[Reed2024_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Reed2024_validation.json)|参数或验证|
|[scope_status.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/scope_status.tsv)|数值或来源表|
|[source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/source_manifest.tsv)|数值或来源表|
|[validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/validation.json)|参数或验证|
|[Wu2021_celltype_profiles.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Wu2021_celltype_profiles.tsv)|数值或来源表|
|[Wu2021_feature_coverage.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Wu2021_feature_coverage.tsv)|数值或来源表|
|[Wu2021_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Wu2021_validation.json)|参数或验证|
|[Wu_treatment_descriptive_sensitivity.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/Wu_treatment_descriptive_sensitivity.tsv)|数值或来源表|

## 原117基因：来源稳定性

[说明与方法](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/README_CN.md)

|文件|用途|
|---|---|
|[all117_bootstrap_source_stability.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/all117_bootstrap_source_stability.png)|图片|
|[all117_comparison_robustness_appended.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/all117_comparison_robustness_appended.tsv)|数值或来源表|
|[all117_next_actions_CN.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/all117_next_actions_CN.tsv)|数值或来源表|
|[all117_paired_source_distributions.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/all117_paired_source_distributions.png)|图片|
|[all117_source_distributions.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/all117_source_distributions.tsv)|数值或来源表|
|[all117_source_stability.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/all117_source_stability.tsv)|数值或来源表|
|[analysis_spec.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/analysis_spec.json)|参数或验证|
|[checksums.sha256](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/checksums.sha256)|校验/记录|
|[integration_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/integration_validation.json)|参数或验证|
|[numeric_implementation_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/numeric_implementation_validation.json)|参数或验证|
|[README_CN.md](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/README_CN.md)|中文说明|
|[sc_input_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/sc_input_manifest.tsv)|数值或来源表|
|[sc_stability_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/sc_stability_validation.json)|参数或验证|
|[source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/source_manifest.tsv)|数值或来源表|
|[validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/validation.json)|参数或验证|

## 新增39基因：表达来源、点图与热图

[说明与方法](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/README_CN.md)

|文件|用途|
|---|---|
|[analysis_spec.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/analysis_spec.json)|参数或验证|
|[BRCA_新增39基因_单细胞来源.xlsx](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/BRCA_新增39基因_单细胞来源.xlsx)|Excel汇总|
|[checksums.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/checksums.tsv)|数值或来源表|
|[delivery_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/delivery_validation.json)|参数或验证|
|[genes156_sc39_appended.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/genes156_sc39_appended.tsv)|数值或来源表|
|[integration_source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/integration_source_manifest.tsv)|数值或来源表|
|[LPCAT4_stable_identity.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/LPCAT4_stable_identity.tsv)|数值或来源表|
|[new39_cell_source_reader_cn.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_cell_source_reader_cn.tsv)|数值或来源表|
|[new39_cross_study_comparison.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_cross_study_comparison.tsv)|数值或来源表|
|[new39_expression_dotplot.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_expression_dotplot.pdf)|完整图册|
|[new39_expression_dotplot.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_expression_dotplot.png)|图片|
|[new39_expression_heatmap.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_expression_heatmap.pdf)|完整图册|
|[new39_expression_heatmap.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_expression_heatmap.png)|图片|
|[new39_relative_pattern_heatmap.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_relative_pattern_heatmap.pdf)|完整图册|
|[new39_relative_pattern_heatmap.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_relative_pattern_heatmap.png)|图片|
|[new39_source_distributions.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_source_distributions.tsv)|数值或来源表|
|[new39_source_stability.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/new39_source_stability.tsv)|数值或来源表|
|[Pal2021_reprocessed_celltype_profiles.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/Pal2021_reprocessed_celltype_profiles.tsv)|数值或来源表|
|[Pal2021_reprocessed_feature_coverage.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/Pal2021_reprocessed_feature_coverage.tsv)|数值或来源表|
|[Pal2021_reprocessed_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/Pal2021_reprocessed_validation.json)|参数或验证|
|[profile39_adapter.py](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/profile39_adapter.py)|复现脚本|
|[README_CN.md](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/README_CN.md)|中文说明|
|[source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/source_manifest.tsv)|数值或来源表|
|[stability39_adapter.py](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/stability39_adapter.py)|复现脚本|
|[validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/validation.json)|参数或验证|
|[Wu2021_celltype_profiles.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/Wu2021_celltype_profiles.tsv)|数值或来源表|
|[Wu2021_feature_coverage.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/Wu2021_feature_coverage.tsv)|数值或来源表|
|[Wu2021_validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/Wu2021_validation.json)|参数或验证|

## 全156基因：Wu分型背景

[说明与方法](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/README_CN.md)

|文件|用途|
|---|---|
|[all156_subtype_expression_patterns.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/all156_subtype_expression_patterns.pdf)|完整图册|
|[analysis_spec.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/analysis_spec.json)|参数或验证|
|[BRCA_全156基因_单细胞分型比较.xlsx](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/BRCA_全156基因_单细胞分型比较.xlsx)|Excel汇总|
|[checksums.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/checksums.tsv)|数值或来源表|
|[genes156_subtype_appended.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/genes156_subtype_appended.tsv)|数值或来源表|
|[genes156_subtype_comparison.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/genes156_subtype_comparison.tsv)|数值或来源表|
|[README_CN.md](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/README_CN.md)|中文说明|
|[source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/source_manifest.tsv)|数值或来源表|
|[subtype_celltype_coverage.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/subtype_celltype_coverage.tsv)|数值或来源表|
|[subtype_celltype_profiles.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/subtype_celltype_profiles.tsv)|数值或来源表|
|[subtype_expression_page1.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/subtype_expression_page1.png)|图片|
|[subtype_expression_page2.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/subtype_expression_page2.png)|图片|
|[subtype_expression_page3.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/subtype_expression_page3.png)|图片|
|[subtype_expression_page4.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/subtype_expression_page4.png)|图片|
|[subtype_top_and_runner.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/subtype_top_and_runner.tsv)|数值或来源表|
|[validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/validation.json)|参数或验证|
|[within_celltype_subtype_comparison.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/within_celltype_subtype_comparison.tsv)|数值或来源表|

## 全156基因：表达UMAP

[说明与方法](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/README_CN.md)

|文件|用途|
|---|---|
|[all156_expression_UMAP.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/all156_expression_UMAP.pdf)|完整图册|
|[analysis_spec.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/analysis_spec.json)|参数或验证|
|[ASNS_LYPLA1_subtypes_UMAP.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/ASNS_LYPLA1_subtypes_UMAP.pdf)|完整图册|
|[ASNS_LYPLA1_subtypes_UMAP.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/ASNS_LYPLA1_subtypes_UMAP.png)|图片|
|[celltype_reference_UMAP.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/celltype_reference_UMAP.pdf)|完整图册|
|[celltype_reference_UMAP.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/celltype_reference_UMAP.png)|图片|
|[checksums.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/checksums.tsv)|数值或来源表|
|[delivery_provenance.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/delivery_provenance.json)|参数或验证|
|[gene_identity.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/gene_identity.tsv)|数值或来源表|
|[gene_scales_and_pages.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/gene_scales_and_pages.tsv)|数值或来源表|
|[key_genes_umap_1.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/key_genes_umap_1.pdf)|完整图册|
|[key_genes_umap_1.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/key_genes_umap_1.png)|图片|
|[key_genes_umap_2.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/key_genes_umap_2.pdf)|完整图册|
|[key_genes_umap_2.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/key_genes_umap_2.png)|图片|
|[README_CN.md](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/README_CN.md)|中文说明|
|[source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/source_manifest.tsv)|数值或来源表|
|[validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/validation.json)|参数或验证|

## 全156基因：表达与细胞类型对照UMAP

[说明与方法](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/README_CN.md)

|文件|用途|
|---|---|
|[ABHD12_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ABHD12_expression_celltype_pair.png)|图片|
|[ACSL4_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ACSL4_expression_celltype_pair.png)|图片|
|[all156_expression_celltype_pairs.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/all156_expression_celltype_pairs.pdf)|完整图册|
|[analysis_spec.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/analysis_spec.json)|参数或验证|
|[ASNS_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ASNS_expression_celltype_pair.png)|图片|
|[celltype_legend.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/celltype_legend.tsv)|数值或来源表|
|[checksums.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/checksums.tsv)|数值或来源表|
|[delivery_provenance.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/delivery_provenance.json)|参数或验证|
|[ENPP2_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ENPP2_expression_celltype_pair.png)|图片|
|[gene_identity.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/gene_identity.tsv)|数值或来源表|
|[gene_page_index.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/gene_page_index.tsv)|数值或来源表|
|[GLS_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/GLS_expression_celltype_pair.png)|图片|
|[GPCPD1_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/GPCPD1_expression_celltype_pair.png)|图片|
|[GPI_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/GPI_expression_celltype_pair.png)|图片|
|[LPCAT1_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/LPCAT1_expression_celltype_pair.png)|图片|
|[LYPLA1_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/LYPLA1_expression_celltype_pair.png)|图片|
|[LYPLA2_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/LYPLA2_expression_celltype_pair.png)|图片|
|[NNMT_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/NNMT_expression_celltype_pair.png)|图片|
|[PYCR1_expression_celltype_pair.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/PYCR1_expression_celltype_pair.png)|图片|
|[README_CN.md](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/README_CN.md)|中文说明|
|[source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/source_manifest.tsv)|数值或来源表|
|[validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/validation.json)|参数或验证|

## 全156基因：上皮专用三联图

[说明与方法](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/README_CN.md)

|文件|用途|
|---|---|
|[ABHD12_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ABHD12_epithelial_triptych.png)|图片|
|[ACSL4_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ACSL4_epithelial_triptych.png)|图片|
|[all156_epithelial_triptychs.pdf](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/all156_epithelial_triptychs.pdf)|完整图册|
|[analysis_spec.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/analysis_spec.json)|参数或验证|
|[ASNS_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ASNS_epithelial_triptych.png)|图片|
|[checksums.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/checksums.tsv)|数值或来源表|
|[delivery_provenance.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/delivery_provenance.json)|参数或验证|
|[embedding_genes.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/embedding_genes.tsv)|数值或来源表|
|[ENPP2_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ENPP2_epithelial_triptych.png)|图片|
|[epithelial_group_counts.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/epithelial_group_counts.tsv)|数值或来源表|
|[gene_identity.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/gene_identity.tsv)|数值或来源表|
|[gene_page_index.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/gene_page_index.tsv)|数值或来源表|
|[GLS_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/GLS_epithelial_triptych.png)|图片|
|[GPCPD1_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/GPCPD1_epithelial_triptych.png)|图片|
|[GPI_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/GPI_epithelial_triptych.png)|图片|
|[LPCAT1_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/LPCAT1_epithelial_triptych.png)|图片|
|[LYPLA1_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/LYPLA1_epithelial_triptych.png)|图片|
|[LYPLA2_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/LYPLA2_epithelial_triptych.png)|图片|
|[NNMT_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/NNMT_epithelial_triptych.png)|图片|
|[PYCR1_epithelial_triptych.png](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/PYCR1_epithelial_triptych.png)|图片|
|[README_CN.md](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/README_CN.md)|中文说明|
|[source_manifest.tsv](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/source_manifest.tsv)|数值或来源表|
|[validation.json](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/validation.json)|参数或验证|

## 重点基因图：直接查看

这些只是方便浏览的示例，不是筛选后的唯一名单。全部156基因见上方PDF。

- **ASNS**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ASNS_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ASNS_epithelial_triptych.png)
- **LYPLA1**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/LYPLA1_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/LYPLA1_epithelial_triptych.png)
- **ABHD12**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ABHD12_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ABHD12_epithelial_triptych.png)
- **LYPLA2**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/LYPLA2_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/LYPLA2_epithelial_triptych.png)
- **ENPP2**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ENPP2_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ENPP2_epithelial_triptych.png)
- **LPCAT1**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/LPCAT1_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/LPCAT1_epithelial_triptych.png)
- **GLS**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/GLS_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/GLS_epithelial_triptych.png)
- **GPCPD1**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/GPCPD1_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/GPCPD1_epithelial_triptych.png)
- **GPI**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/GPI_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/GPI_epithelial_triptych.png)
- **NNMT**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/NNMT_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/NNMT_epithelial_triptych.png)
- **PYCR1**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/PYCR1_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/PYCR1_epithelial_triptych.png)
- **ACSL4**：[全细胞表达与类型对照](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/ACSL4_expression_celltype_pair.png) · [上皮三联图](https://github.com/Asukasssss/-/blob/8663004f8d98acaf1368badc7436c84b059e4462/results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/ACSL4_epithelial_triptych.png)

## 与统一候选表连接

[BRCA全部结果入口](BRCA_CURRENT_RESULTS_CN.md)保留配对代谢物、患者相关、RNA及解释。不能因为只想看图片，就把缺测、相反结果或历史候选删掉。
