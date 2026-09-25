"""Plot public aggregates only after actual inferCNV inference has completed."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'runtime/plot_dependencies'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
P=ROOT/'results/PDAC/06_EXTERNAL/20260925T141744Z_myeloid_infercnv_v1'
g=pd.read_csv(P/'group_summary.tsv',sep='\t');t=pd.read_csv(P/'SLC6A6_comparison.tsv',sep='\t')
states=['CNV_ABNORMAL_CANDIDATE','REFERENCE_LIKE','UNCERTAIN']
labels=['CNV 异常候选','参考样','不确定'];colors=['#c35442','#35759b','#abb1ba']
cohorts=['GSE263733','GSE278688','GSE242230']
font=FontProperties(fname='C:/Windows/Fonts/msyh.ttc');plt.rcParams['font.family']=font.get_name();plt.rcParams['axes.unicode_minus']=False
fig,axes=plt.subplots(2,3,figsize=(15,9));fig.subplots_adjust(left=.07,right=.97,top=.80,bottom=.16,wspace=.35,hspace=.62)
fig.suptitle('PDAC 髓系 inferCNV 分组与 SLC6A6',fontsize=22,x=.06,ha='left',y=.96,fontweight='bold')
fig.text(.06,.90,'表达推断的 CNV 分组，不等于已确认的恶性／正常身份',fontsize=14,color='#475569')
for k,c in enumerate(cohorts):
    d=g[g.cohort==c].set_index('state').reindex(states);a=axes[0,k]
    a.bar(range(3),d.n_cells,color=colors,width=.65)
    for j,n in enumerate(d.n_cells):a.text(j,n,f'{int(n):,}',ha='center',va='bottom',fontsize=11)
    a.set_xticks(range(3),labels,fontsize=10);a.set_title(c,fontsize=14);a.set_ylabel('纳入分析的髓系细胞数');a.margins(y=.18)
    b=axes[1,k]
    for j,(_,r) in enumerate(d.iterrows()):
        if pd.notna(r.donor_equal_mean_expression):
            b.bar(j,r.donor_equal_mean_expression,color=colors[j],width=.65)
            b.text(j,r.donor_equal_mean_expression,f"n={int(r.n_units_ge20cells)}",ha='center',va='bottom',fontsize=10)
        else:b.text(j,.01,'不可评估',ha='center',va='bottom',fontsize=10,rotation=90,color='#6b7280')
    b.set_xticks(range(3),labels,fontsize=10);b.set_ylabel('SLC6A6 供者等权均值');b.margins(y=.25)
    r=t[t.cohort==c].iloc[0]
    msg=(f"配对单位 n={int(r.n_paired_author_units)}；P={r.p_value:.4g}；q={r.q_value:.4g}" if r.status=='DONE' else f"配对单位 n={int(r.n_paired_author_units)}：不足 5，未计算 P/q")
    b.set_title(msg,fontsize=10,pad=12)
    for ax in [a,b]:
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.text(.06,.055,'CNV 分组要求三种参考设置一致；chr3（含 SLC6A6）不参与 CNV 评分。\n表达尺度为 log1p(每万计数)。n 为每组至少 20 个细胞的作者样本单位数，临床独立性未再次认证。\n不足 3 个单位不展示组均值；组间检验要求至少 5 个配对单位。各研究分别解释，不能用细胞数替代供者数。',fontsize=10,color='#475569',linespacing=1.5)
for ext in ['png','svg']:fig.savefig(P/('SLC6A6_infercnv_comparison.'+ext),dpi=180,facecolor='white')
plt.close(fig)
print(P/'SLC6A6_infercnv_comparison.png')
