# PDAC 功能阶段：NNMT、SLC6A8 与肌酸相关基因（DepMap 24Q4）

## 本轮问题
PDAC 癌细胞是否在常规培养下普遍依赖 SLC6A8？CKB 是否可直接当成其功能上游？

## 输入与范围
[Broad 官方 24Q4 归档](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248)的 Model、CRISPR gene effect/依赖概率及 RNA 四文件；官方 MD5 与服务器 SHA256 见 `source_manifest.tsv`。只取标为 Pancreatic Adenocarcinoma 且 ModelType=Cell Line 的 45 个 ModelID。缺失 CKMT1A 敲除特征单列，不填零。

## 实际结果
SLC6A8 gene effect 中位数 −0.178，效应≤−0.5 为 2/45，依赖概率≥0.9 为 0/45；NNMT、CKB、GATM、GAMT 均无普遍强基线敲除效应。阳性参照 KRAS/RPA3 效应≤−0.5 分别 42/45、45/45。CKB RNA 与 SLC6A8 敲除效应在 43 个 RNA 匹配模型的 Spearman ρ=−0.026。完整范围见 `results.tsv`，关键图见 `figures/15_depmap_24q4_pdac_functional.png`。

## 新手解释
数值越负，敲掉该基因后细胞生长受影响越大；但这只测试标准培养下的癌细胞基线生长，不测试 CAF、低氧、转移定植、药物反应或代谢通量。表达高不能直接等同“敲除必死”。

## 限制与反证
这批模型不支持“PDAC 普遍依赖 SLC6A8”的主张，也没有证据将 NNMT 直接接到 SLC6A8。已发表的 PANC-1 RGX-202 肝定植实验与 DepMap 的基线生长读数终点不同；其 CKB 疗效标志物证据主要来自结直肠癌。[原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC8494442/)。[26Q3 发布说明](https://forum.depmap.org/t/announcing-the-26q3-release/4692)指出文库校正已更新；本轮不是 26Q3。原始矩阵、逐模型明细和患者测量未上传 Git。

## 当前决定
保留“PDAC 基质 NNMT 与癌上皮肌酸利用的空间分工”作为待验证假说；CKB–SLC6A8 可作为癌细胞侧情境性分支，不能当作 PDAC 已证实机制。若后续没有条件性摄取和双扰动交互，分开叙述 NNMT 与肌酸通路。

## 下一步
在选定 PDAC 模型中复测 SLC6A8/CKB 蛋白、不同肌酸供给及低氧条件下的同位素摄取；患者来源 CAF 中扰动 NNMT，随后做 2×2 双扰动与化学读数。26Q3 官方文件可取后复算。

## 复现命令
将四份已校验的官方 CSV 放在同一服务器目录，与 `code/extract_depmap_24q4_targets.py` 同级，运行 `python code/extract_depmap_24q4_targets.py`；把输出的 `depmap_24q4_target_models.tsv` 放到本运行目录 `tables/`，再运行 `python code/analyze_depmap24q4_functional.py` 和 `python code/analyze_depmap24q4_creatine_context.py`。本仓库只保留允许公开的汇总；重算完整报告需取得原始文件。
