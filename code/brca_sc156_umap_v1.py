"""Server-only expression UMAP rendering; export figures and aggregate provenance only."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np,pandas as pd,h5py
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
R=Path(sys.argv[1]);P=R/'public';P.mkdir(exist_ok=True)
sp=importlib.util.spec_from_file_location('helper',R/'brca_sc117_profile_v1.py');helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)
src=R.parent/'20260919T140105Z_scRNA117_v1/source/wu_curated.h5ad'
genes=sorted(pd.read_csv(R/'genes156.tsv',sep='\t').gene);assert len(set(genes))==156
def sha(p):return helper.sha256(p)
digest=sha(src)
with h5py.File(src,'r') as h:
 obs=helper.dataframe(h['obs']);var=helper.dataframe(h['raw/var']);xy=h['obsm/X_umap'][:]
 assert len(obs)==len(xy)==100064 and xy.shape[1]==2 and np.isfinite(xy).all() and obs.index.is_unique
 names=var.feature_name.astype(str).to_numpy();ids=var.index.astype(str).str.split('.').str[0];positions=[];mapping=[]
 for g in genes:
  ix=np.flatnonzero(ids.isin(['ENSG00000176454','ENSG00000291994'])) if g=='LPCAT4' else np.flatnonzero(names==g)
  assert len(ix)==1,(g,len(ix));positions.append(ix[0]);mapping.append(dict(gene=g,feature_id=var.index[ix[0]],source_symbol=names[ix[0]]))
 raw=h['raw/X'];ptr=raw['indptr'][:];ng=len(var);expr=np.zeros((len(obs),156),dtype=np.float32)
 for start in range(0,len(obs),2000):
  stop=min(start+2000,len(obs));lo,hi=ptr[start],ptr[stop]
  a=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[start:stop+1]-lo),shape=(stop-start,ng))
  total=np.asarray(a.sum(axis=1)).ravel();assert (total>0).all()
  expr[start:stop]=np.log1p(a[:,positions].toarray()/total[:,None]*10000)
  if start%20000==0:print('Read',start,flush=True)
assert digest==sha(src)
pd.DataFrame(mapping).to_csv(P/'gene_identity.tsv',sep='\t',index=False)
cm=LinearSegmentedColormap.from_list('expression',['#e4e4e4','#ddd4ef','#a187ca','#65339b','#30045e'])
caps=np.array([np.quantile(expr[:,i][expr[:,i]>0],.99) if (expr[:,i]>0).any() else 1 for i in range(156)])
limits=[xy[:,0].min()-.4,xy[:,0].max()+.4,xy[:,1].min()-.4,xy[:,1].max()+.4]
def style(ax):
 ax.set_xlim(limits[:2]);ax.set_ylim(limits[2:]);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
 for s in ax.spines.values():s.set_color('#dedede');s.set_linewidth(.5)
 ax.set_facecolor('white')
def feature(ax,g,mask=None):
 i=genes.index(g);mask=np.ones(len(obs),dtype=bool) if mask is None else mask;sel=np.flatnonzero(mask);sel=sel[np.argsort(expr[sel,i],kind='stable')]
 im=ax.scatter(xy[sel,0],xy[sel,1],c=expr[sel,i],cmap=cm,vmin=0,vmax=caps[i],s=.45,linewidths=0,rasterized=True)
 style(ax);ax.set_title(g,fontsize=12);return im
focus=['ASNS','GLS','LYPLA1','LYPLA2','ABHD12','ENPP2','LPCAT1','GPCPD1','GPI','NNMT','PYCR1','ACSL4']
for page in range(2):
 fig,axs=plt.subplots(2,3,figsize=(14,9))
 for ax,g in zip(axs.flat,focus[page*6:(page+1)*6]):
  im=feature(ax,g);fig.colorbar(im,ax=ax,fraction=.035,pad=.02)
 fig.suptitle('Wu BRCA | gene expression on existing UMAP',fontsize=16);fig.text(.03,.015,'All 100,064 cells; log1p(counts per10k). Per-gene positive-expression 99th percentile cap; gray = zero/low.\nHigh-expression points drawn last. Cell-level display, not donor-weighted evidence or subtype significance.',fontsize=9)
 fig.tight_layout(rect=[0,.07,1,.94]);fig.savefig(P/('key_genes_umap_'+str(page+1)+'.png'),dpi=180);fig.savefig(P/('key_genes_umap_'+str(page+1)+'.pdf'));plt.close(fig)
with PdfPages(P/'all156_expression_UMAP.pdf') as pdf:
 for page in range(26):
  fig,axs=plt.subplots(2,3,figsize=(14,9))
  for ax,g in zip(axs.flat,genes[page*6:(page+1)*6]):
   im=feature(ax,g);fig.colorbar(im,ax=ax,fraction=.035,pad=.02)
  fig.suptitle('Wu BRCA | all156 expression UMAP | page '+str(page+1));fig.text(.03,.015,'Existing embedding; all cells; high values plotted last. Gray = zero/low; each gene has own color cap.\nNo new clustering or statistical tests. log1p(counts per10k), clipped at positive-value 99th percentile.',fontsize=9)
  fig.tight_layout(rect=[0,.07,1,.94]);pdf.savefig(fig,dpi=130);plt.close(fig);print('Page',page+1,flush=True)
fig,axs=plt.subplots(2,3,figsize=(14,9))
for row,g in enumerate(['ASNS','LYPLA1']):
 for col,st in enumerate(['ER+','HER2+','TNBC']):
  ax=axs[row,col];mask=obs.subtype.eq(st).to_numpy();im=feature(ax,g,mask);ax.set_title(g+' | '+st+' | '+str(mask.sum())+' cells');fig.colorbar(im,ax=ax,fraction=.035,pad=.02)
fig.suptitle('ASNS / LYPLA1 by subtype | same UMAP and color scale within each gene');fig.text(.03,.015,'Cells retain their original coordinates; not separate subtype embeddings. Point density is not donor-weighted.\nSame gene uses same color cap across subtypes; no formal subtype test or treatment adjustment.',fontsize=9);fig.tight_layout(rect=[0,.07,1,.94]);fig.savefig(P/'ASNS_LYPLA1_subtypes_UMAP.png',dpi=180);fig.savefig(P/'ASNS_LYPLA1_subtypes_UMAP.pdf');plt.close(fig)
# Broad author-derived annotation reference; labels at medians, not newly inferred identities.
fig,ax=plt.subplots(figsize=(11,9));colors=plt.get_cmap('tab10');cats=sorted(obs.celltype_major.unique())
for j,c in enumerate(cats):
 sel=obs.celltype_major.eq(c).to_numpy();color=colors(j);ax.scatter(xy[sel,0],xy[sel,1],s=.55,c=[color],linewidths=0,label=c,rasterized=True)
 center=np.median(xy[sel],axis=0);ax.text(*center,c,fontsize=8,ha='center',bbox=dict(facecolor='white',edgecolor=color,alpha=.85,boxstyle='round,pad=.2'))
style(ax);ax.set_title('Wu BRCA | Author-derived celltype_major annotation');ax.legend(loc='center left',bbox_to_anchor=(1,0.5),markerscale=5,frameon=False,fontsize=9);fig.tight_layout();fig.savefig(P/'celltype_reference_UMAP.png',dpi=180);fig.savefig(P/'celltype_reference_UMAP.pdf');plt.close(fig)
scales=pd.DataFrame([dict(gene=g,color_cap=float(caps[i]),positive_cells=int((expr[:,i]>0).sum()),all_cells=len(obs),page=i//6+1) for i,g in enumerate(genes)]);scales.to_csv(P/'gene_scales_and_pages.tsv',sep='\t',index=False)
spec=dict(run_id=R.name,version='sc156_umap_v1',source='Wu2021 curated H5AD; frozen obsm/X_umap',genes=156,cells=100064,normalization='log1p(raw count / all raw-gene library *10000)',color_cap='99th percentile of positive cells for each gene; fixed across subtype panels',all_cells_preserved=True,draw_order='ascending expression',new_embedding=False,new_clustering=False,new_P_q=0,source_sha256=digest,matplotlib=matplotlib.__version__,h5py=h5py.__version__,style='white;small dots;no ticks;light borders;no contour on expression panels to avoid obscuring gradients')
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));(P/'validation.json').write_text(json.dumps(dict(status='DONE',finite_coordinates=True,unique156_gene_ids=True,all100064_cells_preserved=True,source_hash_unchanged=True,pdf_pages=26),indent=2))
pd.DataFrame([dict(path=str(p),sha256=sha(p)) for p in [src,R/'genes156.tsv',R/'brca_sc117_profile_v1.py',Path(__file__)]]).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
(P/'README_CN.md').write_text('''# BRCA全156基因表达UMAP

## 本轮问题
展示候选在Wu单细胞中的表达位置，提供细胞类型参照及ASNS、LYPLA1分型展示。
## 输入与范围
全156基因、全部100064细胞，复用冻结X_umap坐标，不重新降维或聚类。只导出图片与群组汇总，不导出逐细胞数据。
## 实际结果
26页全基因PDF、两张重点基因图、细胞类型参照图及两基因三分型图。gene_scales_and_pages列出每基因页码与色标。
## 新手解释
每个点是细胞，紫色越深表示该基因RNA相对较高，灰色为未检出或低表达。每基因按阳性细胞99分位截断色标，所有细胞保留，超上限仍显示最深色。同基因跨分型共用色标；不同基因色标不同，不直接比较深浅。热图是供者等权，本图是细胞展示，两者不是同一统计权重。
## 限制与反证
UMAP位置和点密度不代表效应、细胞比例或代谢通量；同一供者的很多细胞不是独立患者。高表达点最后绘制可遮盖低值，不用于估算阳性比例。作者标签沿用celltype_major，不重新推断。无新P/q、无亚型差异检验。坐标来自处理对象，未重新核验作者UMAP算法参数。LPCAT4按稳定ID匹配，不解除其他平台身份待核。
## 当前决定
全部候选展示，不按表达值删基因。
## 下一步
用图理解表达位置，统计结论仍查供者级汇总。
## 复现
server165新目录含genes156.tsv、brca_sc117_profile_v1.py及本脚本，执行python3 brca_sc156_umap_v1.py RUN_DIR。源H5AD留服务器。
''',encoding='utf-8')
(P/'.gitattributes').write_text('* -text\n');pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(P.iterdir()) if p.name!='checksums.tsv']).to_csv(P/'checksums.tsv',sep='\t',index=False)
(R/'DONE').write_text('DONE');print('DONE',flush=True)
