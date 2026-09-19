# 乳腺癌第二轮：样本关联稳健性

本轮已经实际运行服务器计算；原 CAMP、原 RNA 差异和原关联的 P/q 均未覆盖。

## 结果

- 原 BH 校正经 statsmodels 独立核对，最大误差小于 4e-16。主关联为 171 项可计算、11 项 q<0.05；缺失可用性敏感性为 171 项、9 项；RNA 差异为 116 项、61 项。不是只对选中的10个基因校正。
- 全部171条可计算关系完成逐标本留一法。原11条显著关系的每次留一计算均保持原方向。留一法仅描述效应敏感性，没有反复删除标本寻找显著结果，也没有声称每次删除后仍通过FDR。
- GEO GSE37751 注释通过精确 GSM accession 接上全部61个肿瘤标本：31个ER阳性、30个ER阴性；108个肿瘤/非肿瘤标签均与CAMP一致。GEO的 Non-tumor 与CAMP的 Normal 明确对应。
- ER调整分析覆盖原174条关系中的171条，单独进行BH校正。GPCPD1—GPC：偏秩相关 -0.5152、q=0.0171；GPI—葡萄糖6磷酸：+0.4492、q=0.0399；GPI—果糖6磷酸：+0.4041、q=0.0399。这是3条关系、2个基因。
- ASNS、SORD、ETNK1、NNMT 的ER调整q约0.0534；PNP约0.0608。不能把阈值附近结果解释成确定有/无作用。9999次置换存在蒙特卡洛误差，本轮不为跨过0.05继续试种子。

## 怎样理解

当前结果支持优先核查GPCPD1、GPI的患者关联是否在更多背景下稳定。它们不是新确认的治疗靶点；ER调整仍使用同一批CAMP标本，不是外部重复验证。GPI的两种相关糖磷酸不算两份独立队列证据。

ASNS等其余候选及全117个基因均保留。不能只用ER调整显著性筛靶点；SORD仍保留原缺失可用性敏感性较弱的标记。已有功能文献结论原样复用，未声称本轮补完117个基因的人工审查。

## 尚未完成

- DepMap官方下载清单接口返回HTTP403，官方历史Figshare接口本次同样403。没有可用基因效应/依赖概率矩阵，本轮所有117个基因的新增DepMap分析均为ACCESS_BLOCKED，不能写成阴性。
- 独立患者身份仍未证实。61个唯一LHC样本标签不自动证明61位独立患者；没有按编号推断肿瘤/正常配对。
- 只调整ER状态，没有调整完整分子亚型、纯度、细胞组成和批次；这些分析尚未完成。本轮未新增生存、药物预测或对接。
- CPTAC原肿瘤/正常对比的组织与批次混杂限制仍有效，不纳入候选加分。

## 文件

- `candidate_comparison_117_v2.tsv`：全部117个基因、旧功能证据和新增ER显著关系数。
- `original_eleven_relations_followup.tsv`：原11条关系的留一法及ER结果；original_q为旧值，q_value为ER调整后的新检验族。
- `robustness_all_relations.tsv`、`er_adjusted_all_relations.tsv`：完整174条记录，不可评估项保留。
- `multiplicity_verification.tsv`、`independent_validation.json`：独立校正与线性代数复核。
- `source_manifest.tsv`、`analysis_spec_v2.json`、`metadata_summary.json`：输入哈希、方法与样本汇总。

## 复现

服务器原项目根：`/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716`。
本轮运行目录：`results/collaborative/BRCA/A/20260919_followup_ER_v2/`。

创建新的运行目录及独占`.running`锁，创建`source`和`aggregate`子目录，复制以下3个脚本到该目录。先执行prepare，再执行分析和验证；不可覆盖已有运行。

```bash
python3 prepare_brca_geo_metadata.py --out "$RUN/source"
OPENBLAS_NUM_THREADS=1 python3 brca_followup_v2.py --root "$ROOT" --out "$RUN"
OPENBLAS_NUM_THREADS=1 python3 verify_brca_followup_v2.py --root "$ROOT" --out "$RUN"
```

本轮服务器直连NCBI超时，最终从可联网电脑获取GEO文本并经SSH内存传输到服务器；没有把患者表达/代谢矩阵下载到本地。prepare脚本支持`--relay-host server165`。初始仅完成留一法的运行目录为`20260919_followup_v2`，以本次有完整GEO注释的ER运行作为交付版本。

参考来源：[GEO GSE37751](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37751)；[DepMap官方批量下载说明](https://forum.depmap.org/t/bulk-download-of-25q2-files/4443)。
