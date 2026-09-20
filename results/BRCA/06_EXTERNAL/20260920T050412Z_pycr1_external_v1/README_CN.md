# BRCA 外部关联与PYCR1细胞背景专题（阶段内交付）

## 本轮问题
检查GPCPD1—GPC、GPI—G6P的外部可用性，并用现有单细胞数据检验PYCR1与E2F及胶原形成转录程序的来源层面关联。

## 输入与范围
用户20260920方向决策包已核对10项文件SHA256，六条关系36项数值与65fdb0f一致。FUSCC作者处理代谢组及2022门户TPM；Wu与Pal现有原始计数和已复用注释。原117基因与CAMP统计不修改。

## 实际结果
FUSCC GPCPD1—GPC按作者患者键连接258例：rho=-0.040575，置换P=0.5224，计划两关系BH q=1。排除既有WES身份警示病例后257例：rho=-0.046428，P=0.4602，q=0.9204。精确G6P在594项极性注释中未找到，GPI关系不可评估，不替换为G1P或F6P。现有CPTAC目录没有对应代谢矩阵，不能做这两条代谢物关联。

PYCR1主分析（来源标签为单位；8项统一BH）：

|数据与细胞|程序|n|rho|q|
|---|---|---:|---:|---:|
|Wu恶性上皮|E2F细胞周期|20|0.484|0.0674|
|Pal恶性上皮|E2F细胞周期|33|0.385|0.0674|
|Wu成纤维|胶原形成|23|0.514|0.0674|
|Pal成纤维|胶原形成|30|0.423|0.0674|

全部8项及区间见pycr1_program_associations.tsv，主分析没有q<0.05。预定Wu未治疗敏感性分析另校正4项：恶性上皮E2F n=16、rho=0.700、q=0.0132；恶性上皮胶原形成rho=0.597、q=0.0322。成纤维两项未达到阈值。不能只展示E2F而忽略恶性上皮的胶原关联。

同来源细胞类型关联差值及点区间见pycr1_celltype_contrasts.tsv，不作正式差异检验。两研究并未一致证明细胞类型差异。全117基因历史列原样保留，只追加本批三基因解释；其余为NOT_RUN，不是阴性。

FUSCC主分析95%区间为[-0.161831,0.085840]，排除警示病例后为[-0.173036,0.080619]。

## 新手解释
本次FUSCC没有支持CAMP中明显的GPCPD1—GPC负相关；不是发现正相关，也不是证明GPCPD1无功能。FUSCC局限于TNBC，和CAMP总体背景不同。外部BH按原计划两关系保留范围，缺项仅在校正计算中作P=1占位，公开表该缺项仍为NA。

PYCR1的正相关表示：在同一细胞类别中，某来源的PYCR1表达较高时，其程序评分往往也较高。E2F评分是细胞周期相关RNA特征；胶原形成评分也是RNA特征，不是实测增殖、胶原沉积、脯氨酸通量或PYCR1功能。方向相同值得保留，但主分析尚未达到预定校正阈值；未治疗子集来自同一Wu数据，不能算独立复现。

## 限制/反证
作者患者键及RNA样本映射明确，仍未核验到相同组织分装/区域。代谢值是作者处理后数据，原缺失填补未逆转；没有额外填补。FUSCC原文报告258份样本有转录组，和实测连接数一致（https://www.nature.com/articles/s41422-022-00614-0）。CAMP的Terunuma队列来源与FUSCC不同（https://www.jci.org/articles/view/71180），但本轮未用跨队列基因型证明完全独立。不是因果或酶活验证。CPTAC缺项结论仅适用于服务器现有资料。

外部未调整Spearman与CAMP调整ER/组成后的部分秩相关不是完全相同估计目标，本轮没有正式比较两者相关系数差异。不能据此断言TNBC特异机制。

单细胞沿用现有计数与注释，Pal使用Chen2026对Pal2021的再注释。仅按作者来源标签汇总，不宣称全部为经核实独立患者；同一标签内合并，未按编号猜测患者。每细胞类型至少20细胞、推断至少10来源。原始计数按来源×细胞类型求和，全基因库量归一为log1p CPM；程序为程序内基因跨来源z值均值，剔除PYCR1，常量基因不参与。两程序来自MSigDB HALLMARK_E2F_TARGETS（200基因）和REACTOME_COLLAGEN_FORMATION（90基因），预先冻结。9999次置换和2000次来源bootstrap；bootstrap固定已算程序及缩放，区间为条件性逐项区间。未新增亚型调整；治疗/亚型覆盖仅作描述。未治疗标签按作者Naïve修正，变更记录见implementation_notes.json，未改变主分析。

## 当前决定
记录GPCPD1外部未支持结果，不增加模型挽救显著性；GPI外部缺项保留。PYCR1只做已固定的有限专题，不扩大全117筛选。定时GitHub监控未恢复。

## 下一步
本批已完成。GPCPD1保留CAMP及功能证据，但降低“普遍患者关联可复现”的判断；GPI等待精确代谢物匹配的独立数据。PYCR1保留为情境相关探索线索，下一步若深入应选择直接干预公开数据检验对应程序，避免继续在同一数据尝试模型获取显著性。尚未执行新的干预数据分析。

## 复现命令
python3 code/brca_pycr1_context_v1.py --root "$RUN_ROOT" --previous "$PREVIOUS_SC_ROOT"
python3 code/brca_external_two_relations_v1.py --root "$RUN_ROOT" --project "$PROJECT_ROOT"

脚本依赖同目录brca_sc117_profile_v1.py、两份sc117配置和固定programs文件；患者与来源标签测量留server165。参数先于结果保存在analysis_spec.json。

复现计算后运行：
```
python3 code/brca_verify_pycr1_external_v1.py --root "$RUN_ROOT" --previous "$PREVIOUS_SC_ROOT"
python code/brca_finalize_pycr1_external_v1.py
python tools/check_repository.py
```
运行路径与源文件哈希见输入JSON、source_manifest.tsv及external_source_manifest.tsv；软件版本见analysis_runtime.json。数值独立核查通过，111个来源×细胞类型聚合与旧分析完全一致。原患者数据仅留server165。森林图：pycr1_program_forest.png、external_GPCPD1_forest.png。
