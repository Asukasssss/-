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
        else:b.text(j,.08,'不可评估',transform=b.get_xaxis_transform(),ha='center',va='bottom',fontsize=10,color='#6b7280')
    b.set_xlim(-.6,2.6)
    b.set_xticks(range(3),labels,fontsize=10);b.set_ylabel('SLC6A6 供者等权均值');b.margins(y=.25)
    r=t[t.cohort==c].iloc[0]
    msg=(f"配对单位 n={int(r.n_paired_author_units)}；P={r.p_value:.4g}；q={r.q_value:.4g}" if r.status=='DONE' else f"配对单位 n={int(r.n_paired_author_units)}：不足 5，未计算 P/q")
    b.set_title(msg,fontsize=10,pad=12)
    for ax in [a,b]:
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.text(.06,.035,'CNV 分组要求三种参考设置一致；chr3（含 SLC6A6）不参与 CNV 评分。\n表达尺度为 log1p(每万计数)。n 为每组至少 20 个细胞的作者样本单位数，临床独立性未再次认证。\n柱形使用各组合格单位；检验仅使用两组均合格的配对单位，因此柱高差不一定等于配对效应。\n不足 3 个单位不展示组均值；检验要求至少 5 个配对单位。各研究分别解释，不能用细胞数替代供者数。',fontsize=10,color='#475569',linespacing=1.5)
for ext in ['png','svg']:fig.savefig(P/('SLC6A6_infercnv_comparison.'+ext),dpi=180,facecolor='white')
plt.close(fig)
print(P/'SLC6A6_infercnv_comparison.png')
if (P/'aggregate_chr_profiles.tsv').exists():
    ch=pd.read_csv(P/'aggregate_chr_profiles.tsv',sep='\t')
    idx=pd.MultiIndex.from_product([cohorts,states],names=['cohort','state'])
    columns=['chr'+str(i) for i in range(1,23)]
    mat=ch.pivot(index=['cohort','state'],columns='chr',values='mean_residual').reindex(index=idx,columns=columns)
    assert mat['chr3'].isna().all()
    a=mat.to_numpy(float);finite=np.abs(a[np.isfinite(a)]);span=max(.05,float(finite.max())) if finite.size else .05
    fig,ax=plt.subplots(figsize=(14,6));fig.subplots_adjust(left=.27,right=.88,top=.82,bottom=.29)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d1d5db')
    im=ax.imshow(np.ma.masked_invalid(a),aspect='auto',cmap=cmap,vmin=-span,vmax=span,interpolation='nearest')
    ax.set_xticks(range(22),[str(i)+('*' if i==3 else '') for i in range(1,23)])
    ax.set_yticks(range(9),[c+'  '+labels[states.index(s)] for c,s in idx])
    ax.set_xlabel('染色体');ax.axhline(2.5,color='white',lw=3);ax.axhline(5.5,color='white',lw=3)
    fig.colorbar(im,ax=ax,fraction=.035,pad=.04,label='平均相对表达残差（参考基线为 0）')
    fig.suptitle('髓系 CNV 分组的染色体表达模式',fontsize=20,x=.05,ha='left',y=.95,fontweight='bold')
    fig.text(.05,.045,'每个样本单位先取组内均值，再在研究内等权汇总；采用 T＋B 参考的 inferCNV 残差。\n灰色表示不可评估或未参与：chr3 全部排除；每组至少 20 个细胞、至少 3 个作者单位才展示。\n不同细胞的相反信号可能被均值抵消；该图不是逐细胞克隆图，也不证明共同克隆。\n颜色不是 DNA 拷贝数。该图不能确认髓系恶性，跨研究不作绝对量比较。',fontsize=10,color='#475569',linespacing=1.5)
    for ext in ['png','svg']:fig.savefig(P/('myeloid_CNV_aggregate_heatmap.'+ext),dpi=180,facecolor='white')
    plt.close(fig)
for name in ['SLC6A6_infercnv_comparison.svg','myeloid_CNV_aggregate_heatmap.svg']:
    svg=P/name
    if svg.exists():svg.write_bytes(('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n').encode('utf-8'))
