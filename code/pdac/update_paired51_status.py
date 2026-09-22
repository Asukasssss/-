from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
def main():
 folder=ROOT/'results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2';df=pd.read_csv(folder/'candidate_scRNA_appended.tsv',sep='\t');summary=json.loads((folder/'summary.json').read_text());sc=json.loads((ROOT/'results/PDAC/06_EXTERNAL/20260922T054000Z_sc_paired51_v2/summary.json').read_text());a=df[(df.gene=='ABAT')&(df.metabolite_name=='gamma-aminobutyrate (GABA)')].iloc[0]
 text=f'''# PDAC 当前进度：配对工作池内部证据与三队列来源已完成

## 本轮问题
采用用户批准路线：真实配对→配对代谢物发现→直接关系映射→CAMP内部关联与配对RNA→全候选单细胞来源→内部候选分层。外部matched-omics与功能干预后置。本轮只推进PDAC，不修改其他癌种或共享阶段编号。

## 输入与范围
- CAMP 27T/12N；作者配对表证明完整CAMP代谢物配对为11对。21个肿瘤有作者明确且互不重复配对编号用于关联，另6个不默认独立。
- 全307项配对结果保留，51项原始P<0.05进入探索池；其中0项配对q<0.05。
- 51项去向：30直接、12仅条件、3身份暂挂、6当前范围未找到可靠精确关系。锁定715关系=357直接+358条件，共538基因；7项没有可检验关系的特征保留占位。

## 实际结果
|批次|结果|边界|
|---|---|---|
|直接关系关联|主354/357可评估，1条q<0.05；可用性342/357可评估，1条q<0.05|21个作者确认单位，仍为内部探索|
|条件关系关联|主350/358可评估，4条q<0.05；可用性342/358可评估，3条q<0.05|独立条件检验族，显著不升级直接生化证据|
|配对RNA|528/538可评估，180条P<0.05，55条q<0.05|同11对；RNA显著性不是入场门槛|
|三队列来源|538基因全部分析；矩阵唯一覆盖{'、'.join(str(x['genes_available']) for x in sc['cohorts'])}个|GSE263733、GSE278688、GSE242230；{sum(x['cells_included'] for x in sc['cohorts'])}肿瘤细胞|
|来源稳定|全类别同top {sc['three_cohort_same_top']}基因；通过bootstrap/检出规则{sc['stable_source_genes']}基因|描述来源，不是疾病差异或关系验证|
|完整整合|722行=715关系+7未解决占位；分层{json.dumps(summary['candidate_classes'],ensure_ascii=False)}|未删弱关联、条件关系或缺测条目|

直接池主q<0.05的关系为GABA—ABAT：rho=0.8487、q=0.0357；可用性n=16、rho=0.8441、q=0.0357。ABAT RNA在11对中10对降低、q=0.0438；GABA配对q=0.2578仍未通过FDR。ABAT当前来源分层{a.candidate_class}，稳定来源规则={a.stable_source}，具体三队列结果及解释见下方发现报告。

交付入口：
- [发现与解释](../../results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2/FINDINGS_CN.md)
- [完整整合说明](../../results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2/README_CN.md)
- [单细胞前候选表](../../results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2/candidate_pre_scRNA.tsv)
- [附加三队列后的候选表](../../results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2/candidate_scRNA_appended.tsv)
- [患者关联及RNA](../../results/PDAC/03_PATIENT/20260922T053500Z_internal_paired51_v2/README_CN.md)
- [三队列来源](../../results/PDAC/06_EXTERNAL/20260922T054000Z_sc_paired51_v2/README_CN.md)

来源级1388条关联n/rho、528条RNA动态规划P及5个BH族核对通过；scRNA独立供者汇总/第一来源/一致性核对详见对应independent_validation.json。预设来源与RNA标签防护已接入本轮实际关联入口。源码、参数、来源哈希、完整NA原因均随各批次保存，患者级矩阵和配对映射仅在服务器。

## 新手解释
统计关系不等于因果。ABAT同向降低与GABA正相关不能直接解释为酶驱动代谢物降低，也可能与组织组成有关。Epithelial/Ductal仅是作者标签；Acinar/Endocrine属于胰腺上皮实质，不等同免疫/间质。A/B/C/D是透明的探索分层，C/D不是删除名单。

## 限制/反证
51项配对代谢物均未通过FDR；21单位样本量小，置换P精度有限。部分化学身份、条件底物关系、临床协变量、恶性CNV身份与跨研究临床身份去重尚未补齐。原冻结CAMP对比编码仍未由本轮重建；新配对分析明确采用肿瘤减同患者正常，不依赖旧效应符号。

三套单细胞仅定位表达背景，不构成同一代谢物—基因关系的外部复现。exact_relation_external_validation、functional_perturbation均NOT_PERFORMED。

## 当前决定
本轮02映射去向、03内部关联/配对RNA、04可用性/留一法、06三队列来源、07候选整合各批次DONE。03/04/05/06/07总阶段仍按PARTIAL理解；整癌种与机制验证尚未完成。PR #2保持草稿、未合并。

历史冻结303/42、旧362关系和所有旧P/q保持原样。旧27标本批次与新21单位/51项工作池的检验范围不同，不把同队列重新分析称独立验证。历史三项单细胞标签整改及历史证据范围台账继续有效，不能混入新关系作为自动验证链。

## 下一步
围绕GABA—ABAT先核查组织组成解释与适合的功能模型；条件关系先补精确底物证据，再针对具体关系选择外部matched-omics与功能干预。其余候选全部保留。

## 复现命令
协议：docs/PDAC/PAIRED51_INTERNAL_V2_LOCK.md。各结果目录README提供独立运行命令；code/pdac中的internal_paired51_v2.py、sc_source_paired51_v2.py仅在新的server165隔离目录执行。禁止覆盖旧批次。
'''
 (ROOT/'reports/PDAC/CURRENT_STATUS_CN.md').write_text(text,encoding='utf-8')
if __name__=='__main__':main()
