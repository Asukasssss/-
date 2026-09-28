"""Scientific plots from aggregate results; UMAP uses server-only admitted cells."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();o=a.out;p=o/'public';f=p/'figures';f.mkdir(exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc');plt.rcParams.update({'font.family':'Noto Sans CJK JP','axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
    fig.savefig(f/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(f/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
d=pd.read_csv(p/'Uhlitz_CNA_subtype_summary.tsv',sep='\t');d=d[d.status=='DONE'].sort_values('donor_equal_mean',ascending=False)
fig,axs=plt.subplots(1,2,figsize=(12,5.5));labels=[r.subtype+' (n='+str(r.donors_ge20)+')' for r in d.itertuples()]
axs[0].barh(labels,d.donor_equal_mean,color='#d95555');axs[0].invert_yaxis();axs[0].set_xlabel('供者等权平均 log1p(CP10K)');axs[0].set_title('全部CNA细胞中的表达')
axs[1].barh(labels,100*d.donor_equal_detection,color='#398c97');axs[1].invert_yaxis();axs[1].set_xlim(0,100);axs[1].set_xlabel('供者等权检出率（%）');axs[1].set_title('LYPLA1被检出的比例')
fig.suptitle('LYPLA1在作者上皮亚群中的分布｜仅限已纳入CNA细胞');fig.text(.05,-.02,'n：该亚群至少20个细胞的供者数；至少3位供者才展示。各亚群供者集合不完全相同。',fontsize=10);fig.tight_layout();save(fig,'01_CNA_subtypes')
d=pd.read_csv(p/'cross_CRC_comparison.tsv',sep='\t');z=d[d.depth_supported_both].copy();z['minimum_depth_abs']=z[['depth_rho_Uhlitz','depth_rho_Lee']].abs().min(axis=1);z=z.sort_values(['minimum_depth_abs','gene'],ascending=[False,True]).head(20)
fig,axs=plt.subplots(1,2,figsize=(13,7));yy=np.arange(len(z))
for ax,co,title in zip(axs,['Uhlitz','Lee'],['Uhlitz：作者CNA上皮','Lee：肿瘤CMS注释上皮']):
    for mode,label,col,off in [('effect','未调整','#b9b9b9',-.17),('depth_rho','调整UMI','#c74949',0),('depth_subtype_rho','调整UMI及亚群','#287b89',.17)]:ax.scatter(z[mode+'_'+co],yy+off,label=label,c=col,s=26)
    ax.set_yticks(yy);ax.set_yticklabels(z.gene);ax.invert_yaxis();ax.axvline(0,color='gray',lw=.6);ax.set_xlim(-.05,.4);ax.set_title(title);ax.set_xlabel('患者内Spearman相关的患者等权中位数');ax.legend(fontsize=9)
fig.suptitle('伴随基因示例：调整后效应往往减弱');fig.text(.02,-.015,'按两队列较小的调整后绝对相关排序展示20项；完整结果保留。不是因果基因或通路验证。',fontsize=10);fig.tight_layout();save(fig,'02_CRC_companions')
focus='CDS1 ETNK2 GPAT4 GPD1L PISD PLA2G10 PLPP5 ASNS UCKL1 NNMT SLC6A6'.split();z=d.set_index('gene').reindex(focus)
cols=['effect_Uhlitz','depth_rho_Uhlitz','depth_subtype_rho_Uhlitz','effect_Lee','depth_rho_Lee','depth_subtype_rho_Lee'];x=z[cols].to_numpy(float)
fig,ax=plt.subplots(figsize=(11,6));im=ax.imshow(np.ma.masked_invalid(x),cmap='RdBu_r',vmin=-.35,vmax=.35,aspect='auto');ax.set_yticks(range(len(focus)));ax.set_yticklabels(focus);ax.set_xticks(range(6));ax.set_xticklabels(['Uhlitz\n未调整','Uhlitz\n调整UMI','Uhlitz\nUMI+亚群','Lee\n未调整','Lee\n调整UMI','Lee\nUMI+亚群']);ax.set_title('关注基因与LYPLA1的患者内相关｜两CRC队列')
for i in range(len(focus)):
    for j in range(6):ax.text(j,i,'NA' if not np.isfinite(x[i,j]) else '%.3f'%x[i,j],ha='center',va='center',fontsize=9)
fig.colorbar(im,ax=ax,label='患者等权中位rho');fig.text(.04,-.02,'NA表示缺少可计算结果；数值可为低于正式供者门槛的描述，实际n与P请读完整表。',fontsize=9);fig.tight_layout();save(fig,'03_focus_genes')
# Frozen UMAP, exact cell join. No generated histology or inferred new labels.
prev=o.parent/'20260928T031628Z_lypla1_umap_v1/private_lypla1_overlay.tsv';allcells=pd.read_csv(prev,sep='\t',index_col=0,keep_default_na=False);m=pd.read_csv(o/'private_Uhlitz_metadata.tsv',sep='\t',keep_default_na=False).set_index('cell_id');assert m.index.is_unique and set(m.index)<=set(allcells.index)
z=allcells.loc[m.index].copy();assert len(z)==4477;z['subtype']=m.subtype;epi=allcells[allcells.main_cell_type=='Epithelial'];pal={'TC1':'#e2a32f','TC2':'#3876b2','TC3':'#6cae7d','TC4':'#be4261','其他作者亚群':'#aaa6ad'};z['display']=z.subtype.where(z.subtype.isin(['TC1','TC2','TC3','TC4']),'其他作者亚群')
fig,axs=plt.subplots(1,2,figsize=(12,6));
for ax in axs:ax.scatter(epi.UMAP1,epi.UMAP2,c='#f0f0f0',s=1,linewidths=0);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
for label,c in pal.items():
    s=z[z.display==label];axs[0].scatter(s.UMAP1,s.UMAP2,c=c,s=4,linewidths=0,label=label+' (n='+str(len(s))+')')
axs[0].legend(fontsize=9);axs[0].set_title('同一批4477个CNA细胞：作者亚群')
zz=z.sort_values('expression');im=axs[1].scatter(zz.UMAP1,zz.UMAP2,c=zz.expression,cmap='Reds',vmin=0,vmax=float(z.expression.max()),s=4,linewidths=0);axs[1].set_title('同一坐标：LYPLA1表达');fig.colorbar(im,ax=axs[1],label='log1p(CP10K)',fraction=.035);fig.suptitle('LYPLA1与CNA细胞亚群的对应');fig.text(.03,.015,'浅灰背景是其他上皮细胞；UMAP不是组织位置，亚群分布不等于患者间普遍升高。',fontsize=10);fig.tight_layout(rect=(0,.04,1,.94));save(fig,'04_CNA_subtypes_UMAP')
print('FOUR_FIGURES_DONE',flush=True)
