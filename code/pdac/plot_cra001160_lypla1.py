"""Plot only public patient-aggregated CRA001160 summaries."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'runtime/plot_dependencies'))
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
P=ROOT/'results/PDAC/06_EXTERNAL/20260927T130900Z_cra001160_lypla1_v1'
r=pd.read_csv(P/'results.tsv',sep='\t');g=pd.read_csv(P/'celltype_summary.tsv',sep='\t')
font=FontProperties(fname='C:/Windows/Fonts/msyh.ttc');plt.rcParams['font.family']=font.get_name();plt.rcParams['axes.unicode_minus']=False
plt.rcParams['pdf.fonttype']=42
fig,axs=plt.subplots(1,3,figsize=(17,6));fig.subplots_adjust(left=.07,right=.97,top=.75,bottom=.26,wspace=.48)
fig.suptitle('CRA001160 · LYPLA1 导管细胞表达',x=.055,y=.96,ha='left',fontsize=23,fontweight='bold')
fig.text(.055,.865,'原作者 type 2：恶性相关导管；type 1：相对非恶性导管。每个患者／标本单位等权。',fontsize=12,color='#475569')
keys=[('Tumor','Ductal cell type 2'),('Tumor','Ductal cell type 1'),('Control','Ductal cell type 1')]
labels=['肿瘤内\ntype 2','肿瘤内\ntype 1','对照胰腺\ntype 1'];colors=['#c65b46','#e5af5a','#44839b']
sel=pd.DataFrame([g[(g.tissue==t)&(g.cell_type==c)].iloc[0] for t,c in keys])
for ax,metric,title,mult in [(axs[0],'log1p_10k','A  患者等权平均表达',1),(axs[1],'detection','B  患者等权检出比例',100)]:
 y=sel[metric+'_mean'].to_numpy()*mult
 err=np.vstack([y-sel[metric+'_ci_lower'].to_numpy()*mult,sel[metric+'_ci_upper'].to_numpy()*mult-y])
 ax.bar(range(3),y,color=colors,width=.62);ax.errorbar(range(3),y,yerr=err,fmt='none',ecolor='#334155',capsize=4)
 ax.set_xticks(range(3),labels,fontsize=11);ax.set_title(title,loc='left',fontsize=13,pad=14)
 for k,(_,row) in enumerate(sel.iterrows()):ax.text(k,0,f"n={int(row.n_units_eligible)}",ha='center',va='bottom',fontsize=11,color='white',fontweight='bold')
 ax.set_ylabel('平均 log1p(每万 UMI)' if mult==1 else '表达大于零的细胞比例（%）');ax.margins(y=.15)
p=r[r.test_family=='LYPLA1_two_primary_contrasts'].reset_index(drop=True);a=axs[2]
for j,row in p.iterrows():
 a.errorbar(row.effect,1-j,xerr=[[row.effect-row.ci_lower],[row.ci_upper-row.effect]],fmt='o',color=colors[j],capsize=5,ms=8)
 a.annotate(f"n={int(row.n)} / {int(row.n_reference)}\nP={row.p_value:.4g}；q={row.q_value:.4g}",(row.effect,1-j),xytext=(0,18),textcoords='offset points',ha='center',fontsize=10)
a.axvline(0,color='#94a3b8',ls='--');a.set_yticks([1,0],['同肿瘤内\ntype 2 − type 1','跨组织来源\n肿瘤 type 2 − 对照 type 1'],fontsize=10)
a.set_ylim(-.65,1.75);a.set_title('C  两项预设主比较',loc='left',fontsize=13,pad=14);a.set_xlabel('患者平均表达差及 95% bootstrap 区间')
for ax in axs:
 ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.text(.055,.105,'每组每单位 ≥20 个细胞；A/B 使用该组所有合格单位，C 的肿瘤内比较仅使用两类细胞均合格的单位。',fontsize=11,color='#475569')
fig.text(.055,.055,'P：患者层面双侧置换；q：两项主比较的 BH 校正。type 1 不等同健康细胞；未新增 CNV 验证。',fontsize=11,color='#475569')
fig.savefig(P/'LYPLA1_CRA001160.png',dpi=180);fig.savefig(P/'LYPLA1_CRA001160.pdf');plt.close(fig)
print(P/'LYPLA1_CRA001160.png')
