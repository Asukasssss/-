# GBM全直接关系肿瘤内关联

## 本轮问题
直接代谢物—基因关系在肿瘤病例中是否共变？

## 输入与范围
171条关系、142基因、74个组织标签一致的独立作者case；原处理代谢值与TPM不再变换。

## 实际结果
{
  "REL_TUMOR_AGE_SEX": {
    "planned": 171,
    "evaluable": 168,
    "P005": 13,
    "q005": 0
  },
  "REL_TUMOR_AVAILABLE": {
    "planned": 171,
    "evaluable": 168,
    "P005": 14,
    "q005": 0
  },
  "REL_TUMOR_PATH_PASS": {
    "planned": 171,
    "evaluable": 168,
    "P005": 10,
    "q005": 0
  },
  "REL_TUMOR_PRIMARY": {
    "planned": 171,
    "evaluable": 168,
    "P005": 14,
    "q005": 0
  }
}

## 新手解释
效应为秩相关，不是倍数、酶活或合成方向。全部关系各家族独立BH。

## 限制/反证
同队列选择后分析；上游候选入口仍依赖身份待核实的正常组。关联本身不使用正常标本。病理PASS敏感性排除6个低质量病例；年龄性别partial rank为近似t检验。未控制IDH/纯度等未冻结协变量，不称独立验证。RNA配对和正常差异暂不可评估，不填假P。

## 当前决定
所有142基因保留进入单细胞，不按P筛减。

## 下一步
全候选细胞来源描述和整合。

## 复现
python gbm_patient_v1.py --out AUDITED_RUN --code-commit COMMIT；同目录需要gbm_discovery_v1.py。
