"""Assemble cell-state delivery report from public aggregates only."""
from pathlib import Path
import pandas as pd,json,hashlib
P=Path('results/BRCA/06_EXTERNAL/20260928T131206Z_lypla1_cellstates_v1')
def table(df):
 return '|'+ '|'.join(df.columns)+'|\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('|'+ '|'.join(str(x) for x in r)+'|' for r in df.itertuples(index=False,name=None))
def main():
 c=pd.read_csv(P/'cross_cohort_companions.tsv',sep='\t');s=pd.read_csv(P/'state_summary.tsv',sep='\t');cv=pd.read_csv(P/'patient_coverage_summary.tsv',sep='\t');hits=c[c.replicated_descriptive];top=pd.read_csv(P/'displayed_companions.tsv',sep='\t').head(15)[['gene','median_r_Wu','median_r_Pal','state_adjusted_median_r_Wu','state_adjusted_median_r_Pal']].round(3);lead=s[s.status.eq('DONE')].sort_values(['cohort','mean_log'],ascending=[True,False]).groupby('cohort').head(5)[['cohort','state','n_patients_eligible','mean_log','detection','fraction_above_own_malignant_mean']].round(3)
 text=f'''# LYPLA1癌细胞亚群与患者内伴随表达

## 1. 本轮问题
LYPLA1是否集中于特定癌细胞亚群？同一患者内伴随LYPLA1变化的基因，能否在Wu与Pal中重复？不新增富集或机制检验。

## 2. 输入与范围
复用Wu2021与Pal2021再注释原始计数、既有亚群标签。Wu为作者Cancer LumA/B、Basal、Her2、Cycling；Pal为Chen再注释epi_*，不称为原Pal作者命名、不强行对应Wu亚群。
Pal使用上一轮GEO患者对应，两个明确区域合并；一个分型/样本身份冲突标签排除；PR-only因患者身份明确保留（本轮不比较临床分型）。全源矩阵与逐细胞/逐患者值留server165。
{table(cv)}

## 3. 实际结果
完整亚群记录见state_summary.tsv。下列为各队列可评估亚群表达最高的五项，不是新筛选名单。
{table(lead)}

达到预先固定描述性重复规则的基因共{len(hits)}个。规则：两研究深度调整中位相关同向且绝对值均>=0.1；各自>=70%患者同方向；各自>=5位可评估患者。不以P/q作该标签门槛，不称靶点验证。
下表是图中展示的共同方向效应较大项；若固定规则计数为0，这些仍只是弱关联，不能称为通过复现的候选。
{table(top)}
完整跨队列比较见cross_cohort_companions.tsv；所有基因、三种模型按队列分别存档。图中只展示前20项；未显示不等于无结果。

## 4. 新手解释
亚群点图颜色为供者/患者等权平均log1p(CP10k)，面积反映检出比例；n是有>=20该亚群细胞的患者标签数。少于3人的类别不可作公共来源定位。
患者×亚群热图为“该亚群均值减去同一患者全部恶性细胞均值”，红色表示相对这位患者自身较高。灰色是缺失/不足20细胞，不是低表达。匿名P编号在两研究间不对应。没有新增聚类。
患者内关联先在每个人中计算，随后对相关系数等权取中位数；一个人细胞多不会在汇总中自动获得更高权重。Pearson r不是倍数。主分析depth调整总UMI与检出基因数；depth_state再加入现有亚群指标。两协变量排除LYPLA1自身计数/检出贡献。
表中P是患者间正负方向的双侧二项符号检验，不是把细胞当独立患者得到的P。q在每队列全基因×三模型范围计算。非等效性检验，未估计中位数CI。

## 5. 限制与反证
这是相同历史队列的新探索，不是新增独立实验。亚群使用现有注释，不能自动当生物学稳定细胞类型；Pal epi_*未获独立命名。亚群不均衡可影响未调整相关；调整亦不能排除所有状态、平台、批次、性别和归一化组成效应。
相关可以包含LYPLA1高低背后的共同调控，不能解释为LYPLA1调控这些基因。患者标签沿用作者/GEO，无基因型核验。零值保留；筛选为基因每患者至少max(10,1%细胞)检出，LYPLA1至少10阳性细胞、每患者>=100恶性细胞。低覆盖为不可评估。
源矩阵SHA沿用已冻结历史manifest，并记录当前文件大小，本轮没有为归档再读全矩阵重新hash。数值计算实际读取raw计数；审计另用显式最小二乘残差与相关核对选定基因，未独立重做所有基因。
公共表字段映射：cohort为队列、gene为基因、n_patients为实际可评估患者标签、median_r为效应、p_value为方向一致性检验、q_cohort_all3_modes为本轮家族q；所有表属BRCA/06_EXTERNAL/{P.name}。单患者亚群和相关矩阵保留private，不上传。

## 6. 当前决定
以完整跨队列表和亚群图选择可重复的伴随特征，不据此宣布靶点有效。若加入亚群后的效应明显变小，保留这项限制，不再调参追显著。

## 7. 下一步
将本轮冻结的特征用于独立队列；不把现有两队列再次作为新的验证。不自动增加通路富集、细胞通讯或机制推断。

## 8. 复现命令
server165新目录建立.running后运行brca_lypla1_cellstates_v1.py；数值完成后运行brca_lypla1_cellstates_present_v1.py。运行目录常量须改为新的唯一目录，禁止覆盖本轮。报告脚本在仓库运行：python code/brca_lypla1_cellstates_report_v1.py。
'''
 (P/'README_CN.md').write_text(text,encoding='utf-8')
 doc=Path('docs/BRCA_CURRENT_RESULTS_CN.md');doc.write_text('最新补充：[LYPLA1癌细胞亚群与患者内伴随表达](../'+P.as_posix()+'/README_CN.md)。复用Wu/Pal亚群，患者内深度调整及亚群敏感性，完整保留跨队列一致与不一致结果。\n\n'+doc.read_text(encoding='utf-8'),encoding='utf-8')
 f=Path('coordination/stages/BRCA.tsv');d=pd.read_csv(f,sep='\t');d.loc[len(d)]=['BRCA','06_EXTERNAL',P.name,'lypla1_cellstates_v1','DONE','Existing malignant states; within-patient raw/depth/state correlations',P.as_posix(),'code/brca_lypla1_cellstates_v1.py','analysis/brca-lypla1-cellstates-20260928','Pal conflict excluded; exploratory associations not target validation','Independent cohort replication of frozen features'];d.sort_values(['stage_id','run_id']).to_csv(f,sep='\t',index=False)
 for f in P.iterdir():
  if f.suffix in ['.tsv','.json','.md']:f.write_bytes(f.read_bytes().replace(b'\r\n',b'\n'))
  assert f.stat().st_size<5_000_000,(f,f.stat().st_size)
 pd.DataFrame([dict(file=f.name,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted(P.iterdir()) if f.is_file() and f.name!='DELIVERY_SHA256.tsv']).to_csv(P/'DELIVERY_SHA256.tsv',sep='\t',index=False)
 print('REPORT_DONE')
if __name__=='__main__':main()
