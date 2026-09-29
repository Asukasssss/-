"""All-cell reference-layout slide, server-only raw access; frozen donor statistics."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
R=Path(sys.argv[1]);P=R/'public';P.mkdir(exist_ok=True)
B=R.parent;H=B/'20260919T140105Z_scRNA117_v1/source/wu_curated.h5ad';D=B/'20260923T040012Z_epithelial_paired20_v1'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def col(g,k):
 x=g[k]
 def dec(a):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])
 if isinstance(x,h5py.Group):
  codes=x['codes'][:];assert (codes>=0).all();return dec(x['categories'][:])[codes]
 return dec(x[:])
cache=R/'private_display_cache.npz'
if cache.exists():
 z=np.load(cache);xy=z['xy'];expr=z['expr'];types=z['types']
else:
 with h5py.File(H) as h:
  xy=h['obsm/X_umap'][:];types=col(h['obs'],'celltype_major');names=col(h['raw/var'],'feature_name')
  ix=np.flatnonzero(names=='LYPLA1');assert len(ix)==1
  raw=h['raw/X'];ptr=raw['indptr'][:];expr=np.empty(len(xy))
  for st in range(0,len(xy),4000):
   en=min(st+4000,len(xy));lo,hi=ptr[st],ptr[en]
   a=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[st:en+1]-lo),shape=(en-st,len(names)))
   expr[st:en]=np.log1p(a[:,ix[0]].toarray().ravel()/np.asarray(a.sum(axis=1)).ravel()*10000)
 np.savez_compressed(cache,xy=xy,expr=expr,types=types)
assert len(xy)==len(expr)==len(types)==100064 and np.isfinite(expr).all()
colors={'B-cells':'#d73027','CAFs':'#4575b4','Cancer Epithelial':'#78b436','Endothelial':'#7b3294','Myeloid':'#ef8a0c','Normal Epithelial':'#43a2ca','Plasmablasts':'#d95f9f','PVL':'#888888','T-cells':'#66a98c'}
labels={'B-cells':'B cells','CAFs':'CAFs','Cancer Epithelial':'Cancer epithelial','Endothelial':'Endothelial','Myeloid':'Myeloid','Normal Epithelial':'Normal epithelial','Plasmablasts':'Plasmablasts','PVL':'PVL','T-cells':'T cells'}
counts={c:int((types==c).sum()) for c in colors};assert sum(counts.values())==100064 and counts['Cancer Epithelial']==24489
paired=pd.read_csv(D/'private/paired_measurements.tsv',sep='\t');paired=paired[(paired.cohort=='Wu2021')&(paired.partition=='ALL')&(paired.gene=='LYPLA1')]
frozen=pd.read_csv(R/'frozen_results.tsv',sep='\t');frozen=frozen[(frozen.cohort=='Wu2021')&(frozen.gene=='LYPLA1')];r=frozen[frozen.analysis_type=='ALL'].iloc[0]
assert len(paired)==r.n==8 and (paired.difference>0).sum()==r.positive_pairs==7 and np.isclose(paired.difference.mean(),r.effect)
frozen.to_csv(P/'results.tsv',sep='\t',index=False)
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont(str(R/n))
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0';RED='#DD3434'
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=13,c=INK,bold=False,ha='left'):f.text(x,y,s,fontsize=size,color=c,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(R/'sysu_logo.png'));ax.axis('off')
text(.95,.930,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.883);line(.05,.103,.883,G,2)
text(.05,.828,'BRCA：从全细胞表达定位到同患者证据',27,G,True)
text(.05,.777,'Wu 2021 · GSE176078 · 全部 100,064 个细胞，沿用原有 UMAP 坐标与细胞注释',12,M)
text(.05,.717,'LYPLA1 expression in all cells',15,INK,True)
text(.417,.717,'All cells · cell types',15,INK,True)
text(.790,.717,'同患者上皮比较',15,G,True)
ax1=f.add_axes([.05,.263,.270,.413]);ax2=f.add_axes([.414,.263,.270,.413])
cap=float(np.quantile(expr[expr>0],.99));order=np.argsort(expr,kind='stable')
ax1.scatter(xy[:,0],xy[:,1],s=.55,c='#e6e9e9',lw=0,rasterized=True)
pos=order[expr[order]>0]
im=ax1.scatter(xy[pos,0],xy[pos,1],c=expr[pos],cmap='viridis',vmin=0,vmax=cap,s=.60,lw=0,rasterized=True)
for c in colors:
 mask=types==c;ax2.scatter(xy[mask,0],xy[mask,1],s=.60,c=colors[c],lw=0,rasterized=True)
# A red outline is a visual locator for the central cancer epithelial cloud, not a gate.
cxy=xy[types=='Cancer Epithelial'];q1,q2=np.quantile(cxy,[.03,.97],axis=0);center=(q1+q2)/2;width,height=(q2-q1)*1.12
for ax in [ax1,ax2]:
 ax.set_aspect('equal');ax.set_xlim(xy[:,0].min()-.6,xy[:,0].max()+.6);ax.set_ylim(xy[:,1].min()-.6,xy[:,1].max()+.6)
 ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1',fontsize=10,color=M);ax.set_ylabel('UMAP 2',fontsize=10,color=M)
 for s in ax.spines.values():s.set_visible(False)
 for k in ['left','bottom']:ax.spines[k].set_visible(True);ax.spines[k].set_color('#ADB7B1');ax.spines[k].set_linewidth(.7)
 ax.add_patch(Ellipse(center,width,height,angle=0,fill=False,ec=RED,lw=1.25,zorder=7))
 ax.annotate('Cancer epithelial',xy=(center[0]-width*.25,center[1]+height*.42),xytext=(.10,.75),textcoords='axes fraction',fontsize=9,color=RED,arrowprops=dict(arrowstyle='->',color=RED,lw=1),zorder=8)
cb=f.add_axes([.334,.298,.009,.344]);f.colorbar(im,cax=cb,extend='max');cb.tick_params(labelsize=8,length=2);cb.set_ylabel('LYPLA1 · log1p(CP10K)',fontsize=9,color=M,labelpad=7)
# Cell type legend under the two all-cell panels avoids squeezing either UMAP.
handles=[Line2D([],[],marker='o',ls='',color=colors[c],markersize=4,label=f'{labels[c]}  ({counts[c]:,})') for c in colors]
leg=f.legend(handles=handles,loc='center',bbox_to_anchor=(.375,.172),ncol=3,frameon=False,fontsize=8.6,handletextpad=.4,columnspacing=1.5,labelspacing=.5)
for t in leg.get_texts():
 if t.get_text().startswith('Cancer epithelial'):t.set_color(colors['Cancer Epithelial']);t.set_weight('bold')
f.add_artist(Line2D([.753,.753],[.150,.684],transform=f.transFigure,color=LINE,lw=.8))
text(.868,.661,'8 位患者 · 配对 P = 0.0391',11.5,G,True,ha='center')
ax=f.add_axes([.801,.350,.153,.260])
for row in paired.itertuples():ax.plot([0,1],[row.nonmalignant,row.malignant],c=UP if row.difference>0 else DOWN,lw=1.2,alpha=.60)
ax.scatter(np.zeros(8),paired.nonmalignant,c=DOWN,s=23,edgecolor='white',lw=.4,zorder=3)
ax.scatter(np.ones(8),paired.malignant,c=UP,s=23,edgecolor='white',lw=.4,zorder=3)
for x,v in [(0,paired.nonmalignant.mean()),(1,paired.malignant.mean())]:ax.plot([x-.12,x+.12],[v,v],c=G,lw=2.7,zorder=4)
ax.set_xlim(-.3,1.3);ax.set_ylim(0,1.28);ax.set_xticks([0,1],['非恶性上皮','恶性上皮'],fontsize=9,color=INK)
ax.set_ylabel('患者平均表达',fontsize=9,color=M);ax.tick_params(length=2,labelsize=9,color=M,labelcolor=M)
for k in ['top','right']:ax.spines[k].set_visible(False)
for k in ['bottom','left']:ax.spines[k].set_color(LINE)
text(.873,.282,'7 / 8',27,G,True,ha='center');text(.873,.239,'患者恶性上皮表达更高',11.5,G,True,ha='center')
text(.873,.193,'平均差值 +0.215',10.5,M,ha='center');text(.873,.162,'每位患者每类 ≥20 个细胞',9,M,ha='center')
line(.05,.95,.115)
text(.05,.081,'同患者比较支持恶性上皮上调；全细胞图同时保留其他细胞的表达背景。',14,G,True)
text(.05,.045,'参照为肿瘤内非恶性上皮；红圈仅定位，不用于划分细胞。名义 P；Naive 子集 6/7 同向，P=0.078。',8.8,M)
text(.96,.045,'11',14,G,True,ha='right')
name='11_BRCA_LYPLA1_全细胞与患者'
for ext in ['png','pdf','svg']:f.savefig(P/(name+'.'+ext),dpi=220)
plt.close(f)
pd.DataFrame([dict(cell_type=c,cells=counts[c],color=colors[c]) for c in colors]).to_csv(P/'coverage.tsv',sep='\t',index=False)
manifest=[dict(path=str(p),sha256=sha(p)) for p in [H,D/'private/paired_measurements.tsv',R/'frozen_results.tsv',Path(__file__)]]
pd.DataFrame(manifest).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
spec=dict(version='brca_allcells_slide_v2',scope='all 100064 cells; nine cell classes',embedding='source obsm/X_umap; unchanged',annotation='Author-derived celltype_major; unchanged',expression='log1p(raw/library*10000)',color='viridis; zero gray; positive99th percentile cap',cap=cap,normalization_unchanged=True,circle='3rd-97th percentile bounding ellipse of cancer epithelial coordinates, expanded1.12; locator only; not cell selection or enrichment',paired='frozen 8 donors, min20 cells per type',new_tests=0,raw_data_exported=False,software=dict(matplotlib=matplotlib.__version__,numpy=np.__version__))
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
(P/'validation.json').write_text(json.dumps(dict(status='PASS',all_100064_retained=True,nine_classes_counts_match=True,exact_same_coordinates_both_panels=True,frozen_paired_mean_and_counts_match=True,no_new_tests=True,visual_review='PENDING'),indent=2))
print(json.dumps(spec),flush=True)
