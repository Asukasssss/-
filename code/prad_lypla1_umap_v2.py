"""Readable LYPLA1 views; reuse exact v1 values, no expression recomputation."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);r=p.parse_args().out;pub=r/'public/06_EXTERNAL';(r/'.running').write_text('prad_lypla1_umap_v2')
old=r.parent/'20260928T030756Z_lypla1_umap_v1';source=r.parent/'20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad'
a=ad.read_h5ad(source,backed='r');o=a.obs.copy();v=pd.read_csv(old/'private/plot_values.tsv.gz',sep='\t').set_index('cell');assert v.index.is_unique and set(v.index)==set(o.index);v=v.loc[o.index]
xy=v[['UMAP1','UMAP2']].to_numpy();assert np.allclose(xy,a.obsm['X_umap'],rtol=0,atol=1e-12)
ex=v.LYPLA1_log1p_cp10k.to_numpy();tum=o['type'].eq('cancer').to_numpy();ep=o.celltype_major_v2.eq('Epithelial').to_numpy()
normal=tum&ep&o.malignant_anno_merged.eq('normal').to_numpy();malignant=tum&ep&o.malignant_anno_merged.eq('malignant').to_numpy()
cap=json.loads((old/'public/06_EXTERNAL/analysis_spec.json').read_text())['color_max'];cmap=LinearSegmentedColormap.from_list('expr',['#e4e4e4','#bfd7ed','#619ccc','#245d9b','#082b60']);norm=Normalize(0,cap)
lo=xy.min(0);hi=xy.max(0);pad=(hi-lo)*.03;rg=np.random.default_rng(20260928)
def axes_clean(ax):
 ax.set(xlim=(lo[0]-pad[0],hi[0]+pad[0]),ylim=(lo[1]-pad[1],hi[1]+pad[1]),xticks=[],yticks=[]);ax.set_aspect('equal')
 for s in ax.spines.values():s.set_color('#dddddd');s.set_linewidth(.6)
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(pub/(name+'.'+ext),dpi=320,bbox_inches='tight',facecolor='white')
 plt.close(fig)
fig=plt.figure(figsize=(11,10));gs=fig.add_gridspec(2,2,height_ratios=[1.4,1],hspace=.35,wspace=.32);colors=['#329c90','#c43a52'];masks=[normal,malignant];names=['Normal epithelium','Malignant epithelium']
for k,(mask,name) in enumerate(zip(masks,names)):
 ax=fig.add_subplot(gs[0,k]);ax.scatter(xy[tum,0],xy[tum,1],s=1,c='#f0f0f0',linewidths=0,rasterized=True)
 ix=rg.permutation(np.flatnonzero(mask));im=ax.scatter(xy[ix,0],xy[ix,1],s=5,c=ex[ix],norm=norm,cmap=cmap,linewidths=0,rasterized=True);axes_clean(ax);ax.set_title(name+f'\n{mask.sum():,} cells',fontsize=12)
cax=fig.add_axes([.92,.57,.015,.26]);fig.colorbar(im,cax=cax,label='LYPLA1 log1p(CP10K)',extend='max')
ax=fig.add_subplot(gs[1,0]);rows=[]
for k,(mask,name,col) in enumerate(zip(masks,names,colors)):
 d=pd.DataFrame({'donor':o.donor_id.to_numpy()[mask],'expr':ex[mask]});stats=d.groupby('donor').agg(n=('expr','size'),mean=('expr','mean'));eligible=stats[stats.n>=20]['mean'].to_numpy()
 ax.boxplot([eligible],positions=[k],widths=.4,showfliers=False,medianprops={'color':'black'},boxprops={'color':col},whiskerprops={'color':col},capprops={'color':col})
 ax.scatter(k+rg.uniform(-.1,.1,len(eligible)),eligible,s=28,color=col,alpha=.85,zorder=3)
 ax.text(k,1.04,f'mean = {eligible.mean():.3f}\nn = {len(eligible)} patients',transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
 rows.append(dict(celltype=name,n_cells=int(mask.sum()),positive_fraction=float((ex[mask]>0).mean()),pooled_mean=float(ex[mask].mean()),n_eligible_patients=len(eligible),patient_equal_mean=float(eligible.mean())))
ax.set_xticks([0,1],['Normal','Malignant']);ax.set_ylabel('Patient mean LYPLA1 log1p(CP10K)');ax.set_ylim(bottom=0);ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[1,1]);grid=np.linspace(0,float(ex.max()),250)
for mask,name,col in zip(masks,names,colors):
 frac=np.array([(ex[mask]>t).mean()*100 for t in grid]);ax.plot(grid,frac,color=col,lw=2,label=name.replace(' epithelium',''))
 ax.scatter([0],[frac[0]],color=col,s=25)
ax.set(xlabel='Expression threshold: log1p(CP10K)',ylabel='Cells above threshold (%)',ylim=(0,100),xlim=(0,ex.max()));ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
fig.suptitle('LYPLA1 | normal versus malignant epithelium in PRAD tumors',fontsize=14,y=.98)
fig.text(.12,.025,'Top: same source UMAP and color scale; random draw order. Bottom left: patients with >=20 cells.\nBottom right: all cells, descriptive distribution. No differential-expression P value is implied.',fontsize=9,color='#555555');fig.subplots_adjust(bottom=.15,top=.91,right=.88)
save(fig,'LYPLA1_epithelial_comparison_v2');pd.DataFrame(rows).to_csv(pub/'epithelial_summary.tsv',sep='\t',index=False)
fig,ax=plt.subplots(figsize=(8,7));axes_clean(ax)
hb=ax.hexbin(xy[tum,0],xy[tum,1],C=ex[tum],reduce_C_function=np.mean,gridsize=100,mincnt=5,cmap=cmap,vmin=0,vmax=1.2,linewidths=0,rasterized=True)
fig.colorbar(hb,ax=ax,label='Mean log1p(CP10K) within hexagon',extend='max',shrink=.65)
ax.set_title('LYPLA1 | local mean expression in tumor UMAP',fontsize=13);fig.text(.1,.02,'Includes zero-expression cells; bins with >5 cells shown. Spatial bin means are not patient means.\nLocal averaging removes high-expression-on-top bias; this is not a differential test.',fontsize=8);fig.subplots_adjust(bottom=.13);save(fig,'LYPLA1_UMAP_local_mean_v2')
spec=dict(status='DONE',version='prad_lypla1_umap_v2',parent_commit='1bf00ee37c5f19e4639bba0b8c59f74fcde16b58',source_values=str(old/'private/plot_values.tsv.gz'),normalization='unchanged log1p(CP10K) from v1',embedding='unchanged source X_umap',cell_selection='tumor tissue; author Epithelial with normal or malignant annotation; altered_benign excluded from two-group comparison',draw_order='seeded random, seed20260928; not expression-sorted',expression_cap=cap,local_mean_map='all tumor cells including zeros; hexbin gridsize100; >5cells per bin; display cap1.2; no use as differential test',patient_panel='mean per patient and class; >=20cells; groups not restricted to matched patients',statistics='descriptive only; no P values',validation='exact cell-ID join; source coordinates equal within 1e-12; cached expression unchanged')
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
pd.DataFrame([{'source':str(f),'sha256':sha(f)} for f in [old/'private/plot_values.tsv.gz',Path(__file__)] ]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
(pub/'validation.json').write_text(json.dumps(dict(status='DONE',n_normal=int(normal.sum()),n_malignant=int(malignant.sum()),n_tumor=int(tum.sum()),unchanged_expression=True,coordinates_verified=True),indent=2));a.file.close();print(pd.DataFrame(rows).to_string(index=False))
