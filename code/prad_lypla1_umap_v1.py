"""Plot source UMAP and LYPLA1; source matrices and cell coordinates stay server165."""
import argparse,json,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter,label
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--commit',required=True);args=p.parse_args();r=args.out;pub=r/'public/06_EXTERNAL'
with (r/'.running').open('x') as f:f.write('prad_lypla1_umap_v1')
src=r.parents[0]/'20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
a=ad.read_h5ad(src,backed='r');obs=a.obs.copy();xy=np.asarray(a.obsm['X_umap']).copy();assert obs.index.is_unique and np.isfinite(xy).all()
hits=np.flatnonzero(a.raw.var.feature_name.astype(str).eq('LYPLA1'));assert len(hits)==1 and a.raw.var.index[hits[0]]=='ENSG00000120992'
expr=np.empty(len(obs));det=np.zeros(len(obs),dtype=bool)
for st in range(0,len(obs),2000):
 x=a.raw.X[st:st+2000].tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.round(x.data)).all()
 lib=np.asarray(x.sum(1)).ravel();assert (lib>0).all();v=x[:,hits].toarray().ravel();expr[st:st+len(lib)]=np.log1p(10000*v/lib);det[st:st+len(lib)]=v>0
obs['plot_celltype']=obs.celltype_major_v2.astype(str)
ep=obs.plot_celltype.eq('Epithelial');obs.loc[ep,'plot_celltype']='Epithelial_'+obs.loc[ep,'malignant_anno_merged'].astype(str)
types=sorted(obs.plot_celltype.unique());colors=dict(zip(types,plt.get_cmap('tab20').colors[:len(types)]));colors.update(Epithelial_malignant='#c33248',Epithelial_normal='#3a9b8f',Epithelial_altered_benign='#e4aa37')
names={t:t.replace('Epithelial_','').replace('_',' ')+' epithelium' if t.startswith('Epithelial_') else t.replace('_',' ') for t in types}
vmax=float(np.quantile(expr[expr>0],.99));norm=Normalize(0,vmax);cmap=LinearSegmentedColormap.from_list('lypla1',['#dfdfdf','#f5c3c2','#d64860','#8c1248','#3c003d'])
lo=xy.min(0);hi=xy.max(0);pad=(hi-lo)*.035;lims=[lo-pad,hi+pad]
def base(ax,title):
 ax.set(xlim=(lims[0][0],lims[1][0]),ylim=(lims[0][1],lims[1][1]),xticks=[],yticks=[],title=title);ax.set_aspect('equal');ax.set_facecolor('white')
 for sp in ax.spines.values():sp.set_color('#d8d8d8');sp.set_linewidth(.6)
def feature(ax,mask,title):
 ix=np.flatnonzero(mask);ix=ix[np.argsort(expr[ix],kind='stable')]
 im=ax.scatter(xy[ix,0],xy[ix,1],c=expr[ix],cmap=cmap,norm=norm,s=2.2,linewidths=0,rasterized=True);base(ax,title);return im
def annotation(ax,mask):
 rng=np.random.default_rng(20260928);ix=rng.permutation(np.flatnonzero(mask))
 ax.scatter(xy[ix,0],xy[ix,1],c=[colors[t] for t in obs.plot_celltype.iloc[ix]],s=2.2,linewidths=0,rasterized=True)
 handles=[]
 for n,t in enumerate(types,1):
  z=xy[mask&obs.plot_celltype.eq(t).to_numpy()]
  if not len(z):continue
  handles.append(Line2D([],[],marker='o',ls='',color=colors[t],markersize=4,label=f'{n}  {names[t]}'))
  if len(z)>=150:
   hist,xe,ye=np.histogram2d(z[:,0],z[:,1],bins=100,range=[(lims[0][0],lims[1][0]),(lims[0][1],lims[1][1])]);h=gaussian_filter(hist,1.1);s=np.sort(h.ravel())[::-1];level=s[np.searchsorted(np.cumsum(s),.97*s.sum())]
   if level>0:
    components,ncomp=label(h>=level)
    for cc in range(1,ncomp+1):
     if hist[components==cc].sum()<max(150,.03*len(z)):h[components==cc]=0
    if h.max()>level:ax.contour((xe[1:]+xe[:-1])/2,(ye[1:]+ye[:-1])/2,h.T,levels=[level],colors=[colors[t]],linewidths=.65,linestyles='--',alpha=.65)
   imax=np.unravel_index(np.argmax(h),h.shape);pos=((xe[imax[0]]+xe[imax[0]+1])/2,(ye[imax[1]]+ye[imax[1]+1])/2)
   ax.text(*pos,str(n),ha='center',va='center',fontsize=7,color=colors[t],bbox=dict(boxstyle='round,pad=.17',fc='white',ec=colors[t],lw=.5,alpha=.9))
 base(ax,'Cell types | author-derived display labels');return handles
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(pub/(name+'.'+ext),dpi=320,facecolor='white',bbox_inches='tight')
 plt.close(fig)
cancer=obs['type'].eq('cancer').to_numpy();adj=obs['type'].eq('adj_benign').to_numpy();assert (cancer|adj).all()
for mask,name,title in [(cancer,'LYPLA1_UMAP_cancer','PRAD | tumor tissue'),(np.ones(len(obs),dtype=bool),'LYPLA1_UMAP_all_cells','PRAD | tumor + adjacent tissue')]:
 fig,axes=plt.subplots(1,2,figsize=(13,6));handles=annotation(axes[0],mask);im=feature(axes[1],mask,'LYPLA1 expression');fig.subplots_adjust(left=.02,right=.78,wspace=.1,top=.88,bottom=.11)
 fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.795,.54),frameon=False,fontsize=8,handletextpad=.3)
 cb=fig.colorbar(im,cax=fig.add_axes([.80,.20,.16,.023]),orientation='horizontal',extend='max');cb.set_label('log1p(CP10K)',fontsize=9);cb.ax.tick_params(labelsize=8)
 fig.suptitle(title+f' | {int(mask.sum()):,} cells',fontsize=14)
 fig.text(.025,.035,'Source X_umap reused; no reclustering. Dashed outlines are display guides, not inferred cell states.',fontsize=8,color='#555555')
 fig.text(.80,.075,f'Common scale; upper cap {vmax:.2f}\n99th percentile of positive cells',fontsize=7,color='#555555');save(fig,name)
fig,axes=plt.subplots(1,2,figsize=(11,5.6));feature(axes[0],adj,f'Adjacent tissue | {adj.sum():,} cells');im=feature(axes[1],cancer,f'Tumor tissue | {cancer.sum():,} cells');fig.subplots_adjust(right=.9,top=.88,bottom=.1,wspace=.08)
fig.colorbar(im,cax=fig.add_axes([.93,.23,.015,.5]),extend='max',label='LYPLA1: log1p(CP10K)');fig.suptitle('LYPLA1 | same UMAP coordinates and expression scale',fontsize=13);fig.text(.12,.025,'All cells retained within each tissue; visualization is not a differential-expression test.',fontsize=8);save(fig,'LYPLA1_UMAP_tissue_comparison')
summary=obs[['type','plot_celltype']].copy();summary['expression']=expr;summary['detected']=det
summary.groupby(['type','plot_celltype'],observed=True).agg(n_cells=('expression','size'),pooled_mean_log1p_cp10k=('expression','mean'),detected_fraction=('detected','mean')).reset_index().to_csv(pub/'plot_group_summary.tsv',sep='\t',index=False)
pd.DataFrame({'cell':obs.index,'UMAP1':xy[:,0],'UMAP2':xy[:,1],'LYPLA1_log1p_cp10k':expr}).to_csv(r/'private/plot_values.tsv.gz',sep='\t',index=False)
spec=dict(version='prad_lypla1_umap_v1',parent_commit=args.commit,dataset_version='68b23fda-7191-46a5-8870-819feca3e66e',gene='LYPLA1',ensembl='ENSG00000120992',embedding='depositor X_umap reused unchanged',embedding_parameters='not re-estimated; original hyperparameters not verified',normalization='log1p(10000*raw_gene_count/sum_all_raw_genes_per_cell)',gene_filter='none',cell_filter='none; panel tissue masks only',color_max=vmax,color_cap='99th percentile across all positive cells; common scale; values above cap retained',plot_order='ascending expression; positive cells drawn last',annotation='celltype_major_v2; epithelial split by author malignant_anno_merged; combined display labels author-derived, no new inference',contours='smoothed 100x100 histogram; sigma1.1; outer97pct mass; groups>=150; display only',statistics='NOT_RUN: visualization only; pooled means are not donor-equal estimates or fold changes',software=dict(python=platform.python_version(),anndata=ad.__version__,numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__))
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2));(pub/'validation.json').write_text(json.dumps(dict(status='DONE',n_cells=len(obs),n_cancer=int(cancer.sum()),n_adjacent=int(adj.sum()),gene_unique=True,raw_counts_integer=True,coordinates_finite=True,coordinates_unchanged=True,no_subsampling=True,common_color_scale=True,expression_max=float(expr.max()),clipped_display_cells=int((expr>vmax).sum())),indent=2))
pd.DataFrame([dict(source=str(f),sha256=sha(f)) for f in [src,Path(__file__)]]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
a.file.close();print(json.dumps(spec),flush=True)
