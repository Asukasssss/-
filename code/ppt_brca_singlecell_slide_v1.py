"""Server-only BRCA slide: frozen epithelial embedding and donor summary, no new tests."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy import sparse,ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
R=Path(sys.argv[1]); P=R/'public'; P.mkdir(exist_ok=False)
B=R.parent; E=B/'20260922T092000Z_epithelial_umap_v1'; D=B/'20260923T040012Z_epithelial_paired20_v1'
H=B/'20260919T140105Z_scRNA117_v1/source/wu_curated.h5ad'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def col(g,k):
 x=g[k]
 def dec(v):return np.array([a.decode() if isinstance(a,bytes) else str(a) for a in v])
 if isinstance(x,h5py.Group):return dec(x['categories'][:])[x['codes'][:]]
 return dec(x[:])
meta=pd.read_csv(E/'private/epithelial_cell_metadata.tsv',sep='\t',index_col=0)
xy=np.load(E/'private/epithelial_umap.npy'); assert xy.shape==(28844,2)
with h5py.File(H) as h:
 obs=h['obs']; ids=col(obs,obs.attrs.get('_index','_index')); raw=h['raw/X']; names=col(h['raw/var'],'feature_name')
 ix=np.flatnonzero(names=='LYPLA1');assert len(ix)==1
 order=pd.Index(ids).get_indexer(meta.index);assert (order>=0).all()
 vals=np.empty(len(ids));ptr=raw['indptr'][:]
 for st in range(0,len(ids),4000):
  en=min(st+4000,len(ids));lo,hi=ptr[st],ptr[en]
  a=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[st:en+1]-lo),shape=(en-st,len(names)))
  vals[st:en]=np.log1p(a[:,ix[0]].toarray().ravel()/np.asarray(a.sum(axis=1)).ravel()*10000)
 expr=vals[order]
normal=meta.group.eq('Normal epithelial').to_numpy();assert normal.sum()==4355
paired=pd.read_csv(D/'private/paired_measurements.tsv',sep='\t');paired=paired[(paired.cohort=='Wu2021')&(paired.partition=='ALL')&(paired.gene=='LYPLA1')]
frozen=pd.read_csv(R/'frozen_results.tsv',sep='\t'); frozen=frozen[(frozen.cohort=='Wu2021')&(frozen.gene=='LYPLA1')]
r=frozen[frozen.analysis_type=='ALL'].iloc[0];sens=frozen[frozen.analysis_type=='AUTHOR_NAIVE'].iloc[0]
assert len(paired)==r.n==8 and (paired.difference>0).sum()==r.positive_pairs==7
assert np.isclose(paired.difference.mean(),r.effect) and np.isclose(paired.nonmalignant.mean(),r.mean_nonmalignant)
frozen.to_csv(P/'results.tsv',sep='\t',index=False)
for name in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont(str(R/name))
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B'; INK='#183D34'; M='#6A7973'; UP='#C77462'; DOWN='#468B9E'; LINE='#DDE7E0'
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=13,c=INK,bold=False,ha='left'):f.text(x,y,s,fontsize=size,color=c,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(R/'sysu_logo.png'));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'BRCA：LYPLA1 在恶性上皮中是否更高？',29,G,True)
text(.052,.754,'Wu 2021 · 单细胞上皮定位与同患者比较 · GSE176078',13.5,M)
text(.052,.690,'A  上皮细胞身份',16,G,True);text(.364,.690,'B  LYPLA1 表达',16,G,True);text(.711,.690,'C  同患者表达比较',16,G,True)
text(.052,.651,'28,844 个上皮细胞 · 沿用已有注释',10.5,M)
text(.364,.651,'相同细胞、相同坐标',10.5,M)
text(.711,.651,'8 位患者 · 每位每类 ≥20 个细胞',10.5,M)
axes=[f.add_axes([.05,.293,.286,.321]),f.add_axes([.363,.293,.286,.321])]
idx=np.random.default_rng(20260929).permutation(len(xy))
axes[0].scatter(xy[idx,0],xy[idx,1],c=np.where(normal[idx],DOWN,UP),s=1,lw=0,rasterized=True)
def outline(ax,z,color):
 hist,xe,ye=np.histogram2d(z[:,0],z[:,1],bins=95);smooth=ndimage.gaussian_filter(hist.T,1.4);mask=smooth>smooth.max()*.025
 labs,n=ndimage.label(mask);sizes=np.bincount(labs.ravel());sizes[0]=0;mask=np.isin(labs,np.flatnonzero(sizes>max(25,sizes.max()*.17)))
 ax.contour((xe[:-1]+xe[1:])/2,(ye[:-1]+ye[1:])/2,mask.astype(float),levels=[.5],colors=[color],linewidths=.65,linestyles='--')
outline(axes[0],xy[normal],DOWN);outline(axes[0],xy[~normal],UP)
cap=3.;cm=LinearSegmentedColormap.from_list('exp',['#E5E8EC','#F8D3CE','#E9897F','#C94B54','#7B203A'])
sort=np.argsort(expr,kind='stable');im=axes[1].scatter(xy[sort,0],xy[sort,1],c=expr[sort],cmap=cm,vmin=0,vmax=cap,s=1,lw=0,rasterized=True)
for ax in axes:
 ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
 ax.set_xlim(xy[:,0].min()-.5,xy[:,0].max()+.5);ax.set_ylim(xy[:,1].min()-.5,xy[:,1].max()+.5)
 for s in ax.spines.values():s.set_color(LINE);s.set_linewidth(.7)
handles=[Line2D([],[],marker='o',ls='',color=c,markersize=6,label=t) for c,t in [(UP,'恶性上皮  24,489'),(DOWN,'非恶性上皮  4,355')]]
f.legend(handles=handles,loc='center',bbox_to_anchor=(.19,.253),frameon=False,fontsize=10,ncol=1,labelspacing=.5)
cb=f.add_axes([.396,.260,.213,.012]);f.colorbar(im,cax=cb,orientation='horizontal',extend='max',ticks=[0,1,2,3]);cb.tick_params(labelsize=9,color=M)
text(.505,.213,'单细胞表达 · log1p(CP10K)',10,M,ha='center')
ax=f.add_axes([.741,.330,.194,.264])
for row in paired.itertuples():
 ax.plot([0,1],[row.nonmalignant,row.malignant],c=UP if row.difference>0 else DOWN,lw=1.25,alpha=.62,zorder=1)
ax.scatter(np.zeros(8),paired.nonmalignant,c=DOWN,s=26,zorder=3,edgecolor='white',lw=.5)
ax.scatter(np.ones(8),paired.malignant,c=UP,s=26,zorder=3,edgecolor='white',lw=.5)
for x,v in [(0,paired.nonmalignant.mean()),(1,paired.malignant.mean())]:ax.plot([x-.12,x+.12],[v,v],c=G,lw=3,zorder=4)
ax.set_xlim(-.35,1.35);ax.set_ylim(0,max(paired.malignant.max(),paired.nonmalignant.max())*1.17)
ax.set_xticks([0,1],['非恶性上皮','恶性上皮'],fontsize=10,color=INK);ax.tick_params(axis='both',length=3,color=M,labelcolor=M,labelsize=10)
ax.set_ylabel('患者均值 · log1p(CP10K)',fontsize=10,color=M)
ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
for k in ['left','bottom']:ax.spines[k].set_color(LINE)
text(.839,.607,'配对 P = 0.0391',13,G,True,ha='center')
text(.839,.272,'7 / 8 位患者恶性上皮更高',15,G,True,ha='center')
text(.839,.229,'平均差值 +0.215',12,INK,ha='center')
line(.05,.95,.175)
text(.05,.140,'支持：LYPLA1 在该 BRCA 队列的恶性上皮中上调。',17,G,True)
text(.05,.100,'参照为肿瘤标本内非恶性上皮，并非健康乳腺；作者 Naive 子集同向 6/7，P = 0.078。',10.3,M)
text(.05,.059,'UMAP 复用既有上皮嵌入（未做批次校正）；配对线按患者汇总，横杠为组均值。双侧配对检验，展示名义 P 值。',8.8,M)
text(.95,.045,'11',14,G,True,ha='right')
name='11_BRCA_LYPLA1_单细胞'
for ext in ['png','pdf','svg']:f.savefig(P/(name+'.'+ext),dpi=200)
plt.close(f)
sources=[H,E/'private/epithelial_umap.npy',E/'private/epithelial_cell_metadata.tsv',D/'private/paired_measurements.tsv',R/'frozen_results.tsv',Path(__file__)]
pd.DataFrame([dict(path=str(p),sha256=sha(p)) for p in sources]).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
pd.DataFrame([dict(cell_type='Malignant_epithelial',cells=24489),dict(cell_type='Nonmalignant_epithelial',cells=4355)]).to_csv(P/'coverage.tsv',sep='\t',index=False)
spec=dict(version='brca_sc_slide_v1',source='Wu2021 GSE176078',scope='BRCA only; epithelial subset',cells=28844,coordinates='frozen analyst epithelial UMAP from 20260922T092000Z; no new embedding or clustering',annotation='Author-derived; original major epithelial identities; no new CNV calls',expression='log1p(raw/library*10000)',cap=cap,cap_fraction=float((expr>cap).mean()),paired_donors=8,min_cells_each_type=20,statistics='Reused frozen paired donor results, two-sided exact signed-rank enumeration; nominal P',new_tests=0,software=dict(matplotlib=matplotlib.__version__,numpy=np.__version__),raw_data_exported=False)
(P/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2))
(P/'validation.json').write_text(json.dumps(dict(status='PASS',paired_count_matches=True,paired_mean_matches=True,positive_pairs_matches=True,all_epithelial_cells_retained=True,no_new_embedding=True,no_new_tests=True,visual_review='PENDING'),indent=2))
print(json.dumps(spec,ensure_ascii=False),flush=True)
