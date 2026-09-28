# 实际展示复现顺序
保留仓库代码目录及配置引用的固定输入提交3c64c24ec171b2753b0226fdd5aafb7766e45b00。此处脚本副本用于审阅，执行时放回原仓库目录。
1. Python code/camp_results_report.py --config configs/COAD_source_report_v2.yaml --out <全新目录>。
2. Node code/coad/build_source_contract_workbook.mjs <同一输出目录>，使用@oai/artifact-tool生成完整Excel。
3. Python code/coad/refine_source_report_text_v2.py --out <同一输出目录>，嵌入本机Microsoft YaHei字体。
4. 完成页面、工作表及数据回读检查，更新验收；最后调用展示脚本package函数重建ZIP及逐文件校验。
Python需要pandas、numpy、matplotlib、PyMuPDF、PyYAML；Excel构建需要Node与@oai/artifact-tool。本次Python额外依赖位于任务runtime/report_dependencies，经sys.path引入；未改写统计。openpyxl仅用于独立读取验证，不用于制作工作簿。字体不随包分发。
患者、细胞和供者层原始数据不在本包。来源计算复现脚本与规范位于各阶段Git提交中，报告只读其公开汇总。
