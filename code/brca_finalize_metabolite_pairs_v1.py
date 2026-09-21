"""Publish aggregate paired318 results and descriptive direction plots."""
from pathlib import Path
import argparse,json
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--run',required=True);a=p.parse_args();out=a.repo/'results/BRCA/04_ROBUSTNESS'/a.run
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
d=pd.read_csv(out/'paired_metabolite318.tsv',sep='\t');x=d[d.test_family.eq('paired_processed')].copy();y=d[d.test_family.eq('paired_author_available')].copy()
assert len(x)==318 and len(y)==318
cols=['metabolite_name','n','effect','median_paired_difference','paired_rank_biserial','q_value','pairs_higher','pairs_lower','pairs_equal','normal_approximation_warning']
c=x.merge(y[cols],on='metabolite_name',suffixes=('','_available'),validate='one_to_one')
c['available_direction_matches_primary']=np.sign(c.effect)==np.sign(c.effect_available)
c['mean_direction_matches_signed_rank']=np.sign(c.effect)==np.sign(c.paired_rank_biserial)
c['old_significant']=c.old_q.lt(.05);c['new_significant']=c.q_value.lt(.05)
save(c,out/'comparison318.tsv')
names=['glutamine','glutamate','asparagine','creatine','glycerophosphorylcholine (GPC)','glucose 6-phosphate','fructose-6-phosphate','guanine','hypoxanthine','malate','1-methylnicotinamide']
f=c.set_index('metabolite_name').loc[names].reset_index();save(f,out/'focus11_metabolites.tsv')
summary=dict(primary_significant=int(c.new_significant.sum()),old_significant=int(c.old_significant.sum()),retained_significant=int((c.new_significant&c.old_significant).sum()),new_only=int((c.new_significant&~c.old_significant).sum()),old_only=int((~c.new_significant&c.old_significant).sum()),
 primary_mean_direction_vs_old_changed=int((np.sign(c.effect)!=np.sign(c.old_hedges_g)).sum()),available_evaluable=int(y.p_value.notna().sum()),available_significant=int(y.q_value.lt(.05).sum()),
 both_significant=int((c.new_significant&c.q_value_available.lt(.05)).sum()),both_significant_same_mean_direction=int((c.new_significant&c.q_value_available.lt(.05)&c.available_direction_matches_primary).sum()),
 available_vs_primary_mean_direction_changed=int((c.q_value_available.notna()&~c.available_direction_matches_primary).sum()),primary_approximation_warnings=int(c.normal_approximation_warning.sum()),available_approximation_warnings=int(c.normal_approximation_warning_available.sum()),
 primary_significant_mean_rank_direction_disagree=int((c.new_significant&~c.mean_direction_matches_signed_rank).sum()))
(out/'comparison_summary.json').write_bytes((json.dumps(summary,indent=2)+'\n').encode())
fig,ax=plt.subplots(figsize=(10,6));pos=np.arange(len(f));labels=['Glutamine','Glutamate','Asparagine','Creatine','GPC','G6P','F6P','Guanine','Hypoxanthine','Malate','1-MNA']
ax.barh(pos,f.pairs_higher,color='#CF6354',label='Tumor higher');ax.barh(pos,f.pairs_lower,left=f.pairs_higher,color='#397FAD',label='Tumor lower');ax.barh(pos,f.pairs_equal,left=f.pairs_higher+f.pairs_lower,color='#BBBBBB',label='Equal')
for i,r in f.iterrows():
 ax.text(45.5,i,f'{int(r.pairs_higher)} / {int(r.pairs_lower)} / {int(r.pairs_equal)}',va='center',fontsize=9)
ax.set(yticks=pos,yticklabels=labels,xlim=(0,55),xlabel='Number of author-matched pairs (45 total)',title='Paired metabolite directions — author processed values');ax.invert_yaxis();ax.spines[['top','right']].set_visible(False);ax.legend(loc='upper center',bbox_to_anchor=(.5,-.12),ncol=3,frameon=False);fig.tight_layout();fig.savefig(out/'focus_pair_counts.png',dpi=180);plt.close(fig)
lines=['# BRCA 全318代谢物：45对配对分析','','## 本轮问题','原代谢物变化在同患者肿瘤—正常比较中是否保留？有多少对支持升高/降低？','','## 输入与范围','固定原癌种表318个特征；队列底表有322行，其中4行在原流程已因重复代谢物键排除，本轮未新增筛选。沿用作者PDF明确对应的45对，排除标签冲突标本。仅去除名称首尾空白，与原脚本trimws一致，未改变化学身份。',
'肿瘤、正常处理表与同一data_imputed矩阵逐项一致，最大误差0；不新增变换、标准化或填补。主分析保留作者填补值，敏感性仅取两侧在data表均有限值的配对；data可用性不是经过核实的原始检出掩码。',
'配对Wilcoxon符号秩检验，双侧、零差排除、并列秩与连续性校正、正态近似；至少8个完整对。少于10个非零差值的近似风险单列警示。两个318项家族分别对可评估项BH。效应为作者处理尺度的平均配对差，4000次整对bootstrap区间；区间针对均值，不是符号秩检验反演区间。',
'', '## 实际结果',f"主分析318项全部可评估、每项45对，{summary['primary_significant']}项q<0.05；原非配对显著{summary['old_significant']}项。其中两版本均显著{summary['retained_significant']}项，新增{summary['new_only']}项，原显著而本轮未过线{summary['old_only']}项。",f"原数据可用性敏感性{summary['available_evaluable']}项可评估、{summary['available_significant']}项显著。两套均显著且均值方向一致{summary['both_significant_same_mean_direction']}项；57项在可评估子集中的均值方向不同，不能把176项都称为稳健发现。敏感性改变了样本子集，因此不能直接归因为填补错误。",'','|代谢物|升高对数|降低对数|持平|主q|原数据可用对数|敏感性q|','|---|---:|---:|---:|---:|---:|---:|']
for _,r in f.iterrows():lines.append(f'|{r.metabolite_name}|{r.pairs_higher}|{r.pairs_lower}|{r.pairs_equal}|{r.q_value:.4g}|{r.n_available}|{r.q_value_available:.4g}|')
lines += ['', '## 新手解释','例如谷氨酰胺39/45对降低，表示39位患者的肿瘤值低于自己的正常组织，不表示39个患者各自检验显著。配对RNA、肿瘤内代谢物RNA相关、配对代谢物差异是不同问题，不能互相替代。',
'', '## 限制/反证','同队列补充分析，不是外部验证。癌旁组织不是健康人组织；作者病例身份并非基因型核验。高缺失、填补和可用子集选择均可能影响方向；必须同时查看两套结果。符号秩位置解释依赖差值分布对称；均值、数量多数及秩方向不必相同，方向分歧另列，不把q显著自动称为平均升高/降低得到证明。原q和新q的检验范围、设计、样本量不同。',
f"主分析近似警示{summary['primary_approximation_warnings']}项，敏感性警示{summary['available_approximation_warnings']}项；主显著但均值与符号秩方向不同{summary['primary_significant_mean_rank_direction_disagree']}项。",
'', '## 当前决定','保留原318项及原统计，不重排117个基因为新靶点名单。方向计数是作者处理值结果；缺失敏感条目需单独解释。先将配对证据作为新的一层，与已有患者关联和外部结果对照。',
'', '## 下一步','优先解释原数据完整且配对方向较一致的代谢物；依赖填补或子集变化的条目保留限制，不反复换检验追求显著。',
'', '## 复现命令',f'`python3 brca_metabolite318_paired_v1.py --out SERVER_RUN_DIR`',f'`Rscript brca_check_metabolite_pairs_v1.R SERVER_RUN_DIR`',f'`python3 brca_validate_metabolite_pairs_v1.py --out SERVER_RUN_DIR`',f'`python code/brca_finalize_metabolite_pairs_v1.py --repo . --run {a.run}`','新目录需本任务独占锁；来源及SHA256、参数、软件版本和独立R核对见相应文件。完整样本/差值表仅在server165。']
(out/'README_CN.md').write_bytes(('\n'.join(lines)+'\n').encode());print(json.dumps(summary,indent=2));print(f[['metabolite_name','n_available','effect_available','q_value_available']].to_string(index=False))
