# WPS第4页样本数补充

## 本轮问题
用户要求直接在WPS补充各队列样本数与患者配对数。
## 输入与范围
当前打开的25页生化插图精修PPT第4页，使用WPS COM原生文本对象修改。
## 实际结果
BRCA 108份/45对；COAD 66份/33对；PDAC 39份/11对；PRAD 137份/43对；ccRCC3 114份/17位配对患者，ccRCC4 71份/12位；GBM 80份，74肿瘤和6非配对参照，0对。已保存当前PPT并更新预览图和PDF，另存此目录快照。
## 新手解释
标本份数不同于患者人数。普通两组织配对分别为90、66、22、86份；ccRCC为患者内多区域采样，不按两倍配对人数计算标本数。
## 限制/反证
BRCA 108份包含1份组织标签冲突。PDAC全部39份与11对可配对范围不同。队列规模不是每一检验的有效n；ccRCC两队列不合并独立患者数。
## 当前决定
直接保存用户当前WPS文件；不修改分析或统计结果。
## 下一步
组会放映审阅。
## 复现命令
打开指定PPT后运行code/ppt_wps_sample_counts_v5.ps1。
## 来源
- BRCA：本仓库results/BRCA/07_INTEGRATION/20260922T140000Z_report_v1/figure_data/01_design.tsv与BRCA_主报告.md。
- COAD：Git提交6214e20636273c43b14d033f75d3e3444d0cee35，results/COAD/07_INTEGRATION/20260922T150804Z_source_report_v2/figure_data/01_design.tsv。
- PDAC：Git提交adedb6ba47bb7b20bb0780f4a96f9e7d3e1b167e，results/PDAC/07_INTEGRATION/20260922T140000Z_unified_report_v5/design_summary.tsv与sample_audit_summary.tsv。
- PRAD：../camp-prad-analysis/results/PRAD/01_CAMP/20260922T142000Z_discovery_v1/sample_audit_summary.tsv。
- ccRCC：../camp-ccrcc-analysis/results/ccRCC/01_CAMP/20260925T151000Z_ccrcc3_v1/sample_audit_summary.tsv与20260925T155000Z_ccrcc4_histology_v2/sample_audit_summary.tsv。
- GBM：../camp-gbm-analysis/results/GBM/01_CAMP/20260925T103000Z_discovery_v1/sample_audit_summary.tsv（74/6），20260925T144100Z_identity_resolution_v2/sample_audit_summary.tsv（组织身份修正）。
