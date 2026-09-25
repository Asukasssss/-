# GBM正常身份修订与非配对RNA背景

## 本轮问题
解决正常标签冲突并按真实非配对设计补充RNA背景。

## 输入与范围
{
  "status": "DONE",
  "mapped_author_cases": 80,
  "master_TN_matches_explicit_CPTAC_loader": 80,
  "processed_sampleanno_conflicts": 6,
  "conflict_resolution": "CAMP explicit master TN retained; independent CPTAC maintained package explicit PT tissue rule corroborates; derived sampleanno GROUP not used",
  "source_loader_url": "https://raw.githubusercontent.com/PayneLab/cptac/v0.9.7/cptac/gbm.py",
  "source_loader_sha256": "10769ec041b5f03add366e7b9ab514b5cfa91145e7a8cdda8f76b5e0315a9dd0",
  "paired_design": false,
  "normal_reference": "unmatched GTEx normal brain per Wang2021",
  "original_files_modified": false,
  "genotype_verification": false
}

## 实际结果
{
  "RNA_UNPAIRED_AGE_SEX": {
    "planned": 142,
    "evaluable": 139,
    "P005": 87,
    "q005": 81
  },
  "RNA_UNPAIRED_CURRENT": {
    "planned": 142,
    "evaluable": 139,
    "P005": 85,
    "q005": 75
  }
}

## 新手解释
CPTAC维护包对PT标本有明确正常定义；本轮沿用已有作者MasterMapping，记录processed sampleanno页错误，不更改源文件。不按编号猜配对。RNA效应是log2(TPM+1)均值差，不称酶活或精确表达倍数。

## 限制/反证
正常仅6例、GTEx与肿瘤来源混杂。Welch/年龄性别HC3分析仅为探索背景，不能消除来源混杂，不是原始counts负二项DE。无配对RNA设计。

## 当前决定
上游条件性身份标记在本修订审计中已解决；此前文件保持历史记录，整合表引用本次结论。

## 下一步
完成全基因两研究细胞来源并整合。

## 复现
python gbm_identity_RNA_v2.py --out AUDITED_RUN --code-commit COMMIT。
