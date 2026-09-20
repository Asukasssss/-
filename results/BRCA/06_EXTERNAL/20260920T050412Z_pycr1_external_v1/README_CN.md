# BRCA 外部关联与PYCR1细胞背景专题（阶段内交付）

## 本轮问题
检查GPCPD1—GPC、GPI—G6P的外部可用性，并用现有单细胞数据检验PYCR1与E2F及胶原形成转录程序的来源层面关联。

## 输入与范围
用户20260920方向决策包已核对10项文件SHA256，六条关系36项数值与65fdb0f一致。FUSCC作者处理代谢组及2022门户TPM；Wu与Pal现有原始计数和已复用注释。原117基因与CAMP统计不修改。

## 实际结果
FUSCC GPCPD1—GPC按作者患者键连接258例：rho=-0.040575，置换P=0.5224，计划两关系BH q=1。排除既有WES身份警示病例后257例：rho=-0.046428，P=0.4602，q=0.9204。精确G6P在594项极性注释中未找到，GPI关系不可评估，不替换为G1P或F6P。现有CPTAC目录没有对应代谢矩阵，不能做这两条代谢物关联。
PYCR1细胞计数已汇总，完整数值与最终说明将在同目录追加；本批不冒充完成该专题。

## 新手解释
本次FUSCC没有支持CAMP中明显的GPCPD1—GPC负相关；不是发现正相关，也不是证明GPCPD1无功能。FUSCC局限于TNBC，和CAMP总体背景不同。外部BH按原计划两关系保留范围，缺项仅在校正计算中作P=1占位，公开表该缺项仍为NA。

## 限制/反证
作者患者键及RNA样本映射明确，仍未核验到相同组织分装/区域。代谢值是作者处理后数据，原缺失填补未逆转；没有额外填补。FUSCC原文报告258份样本有转录组，和实测连接数一致（https://www.nature.com/articles/s41422-022-00614-0）。CAMP的Terunuma队列来源与FUSCC不同（https://www.jci.org/articles/view/71180），但本轮未用跨队列基因型证明完全独立。不是因果或酶活验证。CPTAC缺项结论仅适用于服务器现有资料。

## 当前决定
记录GPCPD1外部未支持结果，不增加模型挽救显著性；GPI外部缺项保留。PYCR1只做已固定的有限专题，不扩大全117筛选。定时GitHub监控未恢复。

## 下一步
完成PYCR1来源层面程序相关和癌细胞—成纤维细胞关联差值描述，合并研究判断。

## 复现命令
python3 code/brca_pycr1_context_v1.py --root "$RUN_ROOT" --previous "$PREVIOUS_SC_ROOT"
python3 code/brca_external_two_relations_v1.py --root "$RUN_ROOT" --project "$PROJECT_ROOT"

脚本依赖同目录brca_sc117_profile_v1.py、两份sc117配置和固定programs文件；患者与来源标签测量留server165。参数先于结果保存在analysis_spec.json。
