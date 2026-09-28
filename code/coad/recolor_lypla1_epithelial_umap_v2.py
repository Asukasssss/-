"""Recolor frozen epithelial LYPLA1 overlay; preserve values and shared range."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/COAD/B')
OLD=ROOT/'20260928T034349Z_lypla1_epithelial_umap_v1'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--commit',required=True);a=ap.parse_args()
out=a.out;assert out.parent==ROOT and (out/'.running').is_dir();pub=out/'public';pub.mkdir()
f=pd.read_csv(OLD/'private_epithelial_overlay.tsv',sep='\t',index_col=0,keep_default_na=False)
spec=json.loads((OLD/'public/analysis_spec.json').read_text());norm=Normalize(spec['vmin'],spec['vmax'])
assert len(f)==31247 and f.index.is_unique and np.isfinite(f[['UMAP1','UMAP2','expression']]).all().all()
palette=['#f2f2f2','#fcbba1','#fb6a4a','#cb181d','#67000d'];cmap=LinearSegmentedColormap.from_list('gray_red',palette,N=256)
groups=['tumor_CNA','tumor_CNN','normal_reference'];names=['Tumor CNA','Tumor CNN','Normal reference'];expected=[4477,7885,16138]
for g,n in zip(groups,expected):assert (f.group==g).sum()==n
assert np.isclose(f.loc[f.group.isin(groups),'expression'].max(),spec['vmax'])
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
xlim=(f.UMAP1.min()-.7,f.UMAP1.max()+.7);ylim=(f.UMAP2.min()-.7,f.UMAP2.max()+.7)
def clean(ax,title):ax.set(xlim=xlim,ylim=ylim,xlabel='UMAP 1',ylabel='UMAP 2',xticks=[],yticks=[],title=title);ax.set_aspect('equal')
def expr(ax,z,title):
    ax.scatter(f.UMAP1,f.UMAP2,s=1,c='#eeeeee',linewidths=0,rasterized=True)
    zero=z.loc[~z.detected];ax.scatter(zero.UMAP1,zero.UMAP2,s=2,c='#bdbdbd',linewidths=0,rasterized=True)
    p=z.loc[z.detected].sort_values('expression',kind='stable');ax.scatter(p.UMAP1,p.UMAP2,c=p.expression,s=3,cmap=cmap,norm=norm,linewidths=0,rasterized=True);clean(ax,title)
def save(fig,name):
    fig.savefig(pub/(name+'.png'),dpi=250,bbox_inches='tight',facecolor='white');p=pub/(name+'.svg');fig.savefig(p,dpi=180,bbox_inches='tight',facecolor='white');plt.close(fig)
    p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
fig,axes=plt.subplots(1,3,figsize=(16,5.7))
for ax,g,name,n in zip(axes,groups,names,expected):expr(ax,f.loc[f.group==g],'%s | n=%s'%(name,format(n,',')))
fig.subplots_adjust(left=.04,right=.9,bottom=.2,top=.86,wspace=.15);cax=fig.add_axes([.925,.29,.013,.45]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='LYPLA1 log1p(CP10K)')
fig.suptitle('LYPLA1 in epithelial cells | shared scale | gray-to-red',fontsize=16)
fig.text(.04,.035,'Darker red: higher expression. Gray: zero counts. Pale background: other epithelial cells.\nSame cells, coordinates, point sizes and expression range as v1; no clipping or group-specific scaling.',fontsize=9);save(fig,'LYPLA1_epithelial_UMAP_groups_red')
fig,axes=plt.subplots(1,2,figsize=(13,6.5));z=f.sample(frac=1,random_state=20260928)
colors={'tumor_CNA':'#d1495b','tumor_CNN':'#4477aa','normal_reference':'#44aa99','excluded_donor_identity_conflict':'#aaaaaa','not_selected':'#dddddd'}
labels=dict(zip(groups,names));labels.update(excluded_donor_identity_conflict='Excluded: donor ID conflict',not_selected='Other / not admitted')
axes[0].scatter(z.UMAP1,z.UMAP2,c=z.group.map(colors),s=2,linewidths=0,rasterized=True)
for g in colors:
    n=int((f.group==g).sum())
    if n:axes[0].scatter([],[],c=colors[g],s=20,label='%s (n=%s)'%(labels[g],format(n,',')))
clean(axes[0],'Author CNA groups | epithelial cells');axes[0].legend(loc='lower center',bbox_to_anchor=(.5,-.27),frameon=False,fontsize=8,ncol=2)
expr(axes[1],f.loc[f.group.isin(groups)],'LYPLA1 | gray-to-red expression')
fig.subplots_adjust(left=.05,right=.88,bottom=.25,top=.85,wspace=.17);cax=fig.add_axes([.91,.34,.015,.4]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='LYPLA1 log1p(CP10K)')
fig.suptitle('Uhlitz / GSE166555 | epithelial LYPLA1 UMAP',fontsize=16)
fig.text(.05,.025,'Color-only revision; same frozen UMAP and expression. CNN is not proven benign.\nDarker red: higher expression. Gray: zero counts; excluded / unassigned cells are pale background.',fontsize=9);save(fig,'LYPLA1_epithelial_UMAP_red')
shutil.copyfile(OLD/'public/coverage.tsv',pub/'coverage.tsv')
paths=[OLD/'private_epithelial_overlay.tsv',OLD/'public/analysis_spec.json',OLD/'public/coverage.tsv',Path(__file__)]
pd.DataFrame([dict(source_id=p.name,server_path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
spec.update(analysis_version='COAD_LYPLA1_epithelial_UMAP_color_v2',code_lock_commit=a.commit,prior_run=OLD.name,palette=palette,palette_interpolation='linear evenly spaced stops',color_only_change=True,point_size_unchanged=True,software_render=dict(matplotlib=matplotlib.__version__,numpy=np.__version__,pandas=pd.__version__))
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
v=dict(status='PASS',epithelial_cells=len(f),admitted_cells=28500,source_values_and_coordinates_reused=True,group_counts_match_v1=True,color_range_matches_v1=True,no_clipping=True,new_statistics=False,private_records_exported=False)
(pub/'validation.json').write_text(json.dumps(v,indent=2)+'\n')
readme='''# COAD：LYPLA1上皮UMAP红色配色版

## 本轮问题

按用户要求换用更易分辨的表达配色。

## 输入与范围

完整复用上一批20260928T034349Z_lypla1_epithelial_umap_v1的私有上皮坐标、表达及CNA分组，不重算降维、表达或统计。CNA/CNN/正常参照分别4477/7885/16138细胞。

## 实际结果

![三组红色表达图](LYPLA1_epithelial_UMAP_groups_red.png)

![注释与表达](LYPLA1_epithelial_UMAP_red.png)

表达由浅灰渐变至红、深红；颜色越深表示表达越高。零计数仍单独用灰色表示，其他上皮作浅灰背景。点大小、绘制顺序、坐标、三组共享色标上下限与前版完全一致，未做分位数截断或按组缩放。

## 新手解释

从左到右为肿瘤CNA、肿瘤CNN、正常参照。重点观察哪些区域检出表达及颜色强度，不能由亮点覆盖面积推算检出率。

## 限制/反证

换色只改善可读性，不增加组间差异或统计证据。CNN并非已证实非恶性；供者差异、测序深度及降维限制沿用前版。未做批次校正。

## 当前决定

新增两张PNG/SVG，保留前版。不新增P/q或机制解释。逐细胞表留服务器。

## 下一步

使用红色版对照既有表达统计。

## 复现命令

在server165的新独占目录运行`python3 recolor_lypla1_epithelial_umap_v2.py --out <目录> --commit <锁定提交>`。脚本为code/coad/recolor_lypla1_epithelial_umap_v2.py，输入及代码哈希见source_manifest.tsv。
'''
(pub/'README_CN.md').write_text(readme,encoding='utf-8');(out/'DONE').write_text('DONE\n');(out/'.running').rmdir();print(json.dumps(v),flush=True)
