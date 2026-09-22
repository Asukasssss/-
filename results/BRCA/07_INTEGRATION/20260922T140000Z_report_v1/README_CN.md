# BRCA首版结果展示

## 本轮问题
将现有结果生成可读报告和完整附录，不重算。
## 输入与范围
固定提交 8663004f8d98acaf1368badc7436c84b059e4462。
## 实际结果
主报告15页；新图26张；复用4件。详见validation与清单。
## 新手解释
先读主报告，再按完整关系ID查Excel和附录；精选不是排名。
## 限制/反证
来源状态优先机械排名。外部未支持、身份歧义和覆盖缺项保留；其他癌种未实际运行。
## 当前决定
单细胞停在来源，不新增机制分析。
## 下一步
导师讨论；不自动加统计。
## 复现
python code/camp_results_report.py --config configs/BRCA_report_v1.yaml --validate-only

python code/camp_results_report.py --config configs/BRCA_report_v1.yaml --out <全新目录>

ZIP见同级目录。附录原图册在本地完整包；GitHub按来源提交链接复用，避免重复上传大文件。
