# GBM直接生化映射

## 本轮问题
将全部条件性工作特征按直接生化证据映射到人类基因。

## 输入与范围
当前480个P<0.05特征，与已审定BRCA/PRAD/COAD/PDAC数据库关系按精确名称及化学键双重匹配。只复用生化证据，不复用其他癌种统计。

## 实际结果
{
  "status": "PARTIAL",
  "workpool_features": 480,
  "mapped_features": 23,
  "direct_relations": 171,
  "unique_genes": 142,
  "unmapped_features": 457,
  "unique_relation_keys": true,
  "unique_gene_ids": true,
  "statistics_reused_from_other_cancers": false,
  "normal_discovery_conflict_unresolved": true
}

## 新手解释
多篇出处合并为一条特征—基因检验；不同峰不按化学键合并。

## 限制/反证
这不是穷尽式数据库新检索；未映射特征全保留NEEDS_REVIEW。上游正常标签冲突仍在。生化映射不证明GBM功能、通量或干预效果。

## 当前决定
全部直接关系进入肿瘤分析，全部基因进入单细胞，不按后续P二次筛选。

## 下一步
肿瘤内部关系与全基因细胞来源。

## 复现
python code/gbm_mapping_v1.py；使用source_manifest冻结提交。
