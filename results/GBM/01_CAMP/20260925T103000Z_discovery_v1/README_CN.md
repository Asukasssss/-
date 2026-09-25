# GBM发现与身份核查

## 本轮问题
按BRCA主线建立GBM全量发现与候选工作池。

## 输入与范围
CAMP作者处理矩阵、明确MasterMapping、原始GBM.xlsx病例资料；仅服务器读取。

## 实际结果
{
  "status": "PARTIAL",
  "total_input_features": 713,
  "retained_both_tissues": 697,
  "historical_cancer_features": 589,
  "workpool_p005": 480,
  "current_family_q005": 460,
  "tumor_cases": 74,
  "normal_cases": 6,
  "normal_label_conflicts": 6,
  "tumor_pathology_PASS": 68,
  "tumor_pathology_FAIL": 6,
  "RNA_exact_ID_coverage": 80,
  "RNA_gene_symbols": 39309,
  "historical_P_independently_reproduced": true,
  "BH_independently_checked": true,
  "all_source_feature_rows_preserved": true,
  "original_statistics_modified": false,
  "genotype_check": "NOT_EVALUABLE_no_genotype_data",
  "sex_age": "available_from_original_clinical;used_in_next_stage"
}

## 新手解释
原始P/效应精确复用，当前完整特征族的BH另列。不同脂质峰即使同KEGG仍保留，不合并。P<0.05工作池是探索入口。

## 限制/反证
74肿瘤与6正常来自非配对设计。6正常的MasterMapping和processed sampleanno标签冲突，尚未以独立明确组织标签解决。原临床的病例ID均能精确连接，但缺少肿瘤字段不能单独证明正常身份。因此所有依赖正常组的结果为NEEDS_REVIEW，不升级确认性证据。原研究正常来自GTEx，与肿瘤来源有混杂；缺失掩码仅为作者可用值。

## 当前决定
保留冻结结果和条件性全量工作池；不伪造配对、不自动更改组织标签。74个组织一致的肿瘤可做内部关系，另给68个病理PASS子集敏感性。

## 下一步
直接映射、肿瘤内关联、全候选细胞来源；正常依赖分析待身份解决。

## 复现
python gbm_discovery_v1.py --out NEW_SERVER_RUN --code-commit COMMIT。source放GBM_original.xlsx。
