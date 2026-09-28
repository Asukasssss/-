"""Server-only display redraw. Reuse coordinates, cell admission and stored expression."""
from pathlib import Path
import os,json,hashlib,importlib.util
import numpy as np,pandas as pd,h5py
from scipy import sparse,ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
from matplotlib import font_manager
from PIL import Image
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
C=B/'results/collaborative';O=C/'BRCA/A/20260928T150000Z_ppt_redraw_v1';P=O/'public';P.mkdir(exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
plt.rcParams.update({'font.family':'Noto Sans CJK JP','font.size':11,'axes.titlesize':14,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False,'savefig.facecolor':'white'})
RED='#C94B54';BLUE='#3276A5';GRAY='#D7DCE2';DARK='#172A3A';TEAL='#258F88'
CM=LinearSegmentedColormap.from_list('LYPLA1',['#E5E8EC','#F8D3CE','#E9897F','#C94B54','#7B203A'])
CAP=3.0;SEED=20260928;checks=[];sources=[]
def src(p):
 sources.append({'server_path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<100_000_000 else 'Large frozen source: inherited SHA in original manifest'})
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(P/f'{name}.{ext}',dpi=220,bbox_inches='tight')
 plt.close(fig)
def decode(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
def hcol(g,key):
 x=g[key]
 if isinstance(x,h5py.Group):
  cats=decode(x['categories'][:]);codes=x['codes'][:];return np.where(codes>=0,cats[np.maximum(codes,0)],'')
 vals=x[:]
 if '__categories' in g and key in g['__categories']:
  cats=decode(g['__categories'][key][:]);return np.where(vals>=0,cats[np.maximum(vals,0)],'')
 return decode(vals)
def border(ax):
 ax.set(xticks=[],yticks=[]);ax.set_aspect('equal')
 for s in ax.spines.values():s.set_color('#D7DCE2');s.set_linewidth(.6)
def outline(ax,xy):
 hist,xe,ye=np.histogram2d(xy[:,0],xy[:,1],bins=100);z=ndimage.gaussian_filter(hist.T,1.6);mask=z>z.max()*.035
 labs,n=ndimage.label(mask);sizes=np.bincount(labs.ravel());sizes[0]=0
 mask=np.isin(labs,np.flatnonzero(sizes>max(20,sizes.max()*.15)))
 if mask.any():ax.contour((xe[:-1]+xe[1:])/2,(ye[:-1]+ye[1:])/2,mask.astype(float),levels=[.5],colors=['#A4ADB6'],linewidths=.55,linestyles='--')
def umap(name,xy,exp,groups,labels,note):
 assert len(xy)==len(exp)==len(groups) and np.isfinite(xy).all() and np.isfinite(exp).all()
 fig,axs=plt.subplots(1,2,figsize=(11.8,4.1));fig.subplots_adjust(left=.02,right=.89,bottom=.13,top=.88,wspace=.10)
 idx=np.random.default_rng(SEED).permutation(len(xy));palette={k:v[1] for k,v in labels.items()}
 axs[0].scatter(xy[idx,0],xy[idx,1],c=[palette.get(g,GRAY) for g in groups[idx]],s=1.0,lw=0,rasterized=True)
 axs[1].scatter(xy[idx,0],xy[idx,1],c=GRAY,s=1.0,lw=0,rasterized=True)
 shown=idx[groups[idx]!='other']
 h=axs[1].scatter(xy[shown,0],xy[shown,1],c=exp[shown],cmap=CM,vmin=0,vmax=CAP,s=1.0,lw=0,rasterized=True)
 for ax in axs:border(ax);outline(ax,xy)
 axs[0].set_title('细胞身份 / 原有注释',loc='left',color=DARK);axs[1].set_title('LYPLA1 实测表达',loc='left',color=DARK)
 handles=[Line2D([0],[0],marker='o',ls='',color=v[1],markersize=5,label=v[0]+f'  n={sum(groups==k):,}') for k,v in labels.items() if (groups==k).any()]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.47,.035),ncol=min(3,len(handles)),frameon=False,fontsize=9)
 cb=fig.add_axes([.92,.23,.013,.49]);fig.colorbar(h,cax=cb,extend='max',label='log1p(CP10K)')
 save(fig,'umap_'+name)
 checks.append(dict(plot='umap_'+name,cells=len(xy),positive=int(sum(exp>0)),expression_cap=CAP,clipped_fraction=float(np.mean(exp>CAP)),seed=SEED,note=note,coordinate_reuse=True,new_embedding=False))

# BRCA: exact old epithelial coordinates + original raw normalization, no new embedding.
r=C/'BRCA/A/20260922T092000Z_epithelial_umap_v1';meta=pd.read_csv(r/'private/epithelial_cell_metadata.tsv',sep='\t',index_col=0);xy=np.load(r/'private/epithelial_umap.npy')
hp=C/'BRCA/A/20260919T140105Z_scRNA117_v1/source/wu_curated.h5ad';src(hp);src(r/'private/epithelial_umap.npy');src(r/'private/epithelial_cell_metadata.tsv')
with h5py.File(hp) as h:
 obs=h['obs'];ids=hcol(obs,obs.attrs.get('_index','_index'));raw=h['raw/X'];var=h['raw/var'];names=hcol(var,'feature_name');ix=np.flatnonzero(names=='LYPLA1');assert len(ix)==1
 order=pd.Index(ids).get_indexer(meta.index);assert (order>=0).all();vals=np.empty(len(ids));ptr=raw['indptr'][:]
 for st in range(0,len(ids),4000):
  en=min(st+4000,len(ids));lo,hi=ptr[st],ptr[en];a=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[st:en+1]-lo),shape=(en-st,len(names)))
  vals[st:en]=np.log1p(a[:,ix[0]].toarray().ravel()/np.asarray(a.sum(axis=1)).ravel()*10000)
 exp=vals[order]
groups=np.where(meta.group.eq('Normal epithelial'),'normal','malignant')
umap('BRCA',xy,exp,groups,{'malignant':('恶性上皮',RED),'normal':('非恶性上皮',BLUE)},'Wu; previous analyst-computed epithelial embedding; no batch correction; original identity labels')
print('BRCA UMAP done',flush=True)
# COAD cached admitted cells. Excluded cells remain gray background, never assigned malignancy.
hp=C/'COAD/B/20260928T034349Z_lypla1_epithelial_umap_v1/private_epithelial_overlay.tsv';src(hp);d=pd.read_csv(hp,sep='\t');groups=d.group.fillna('not_selected').to_numpy();groups=np.where(np.isin(groups,['tumor_CNA','tumor_CNN','normal_reference']),groups,'other')
umap('COAD',d[['UMAP1','UMAP2']].to_numpy(),d.expression.to_numpy(),groups,{'tumor_CNA':('CNA上皮',RED),'tumor_CNN':('CNN上皮 / 身份未定',TEAL),'normal_reference':('正常组织参照',BLUE),'other':('未纳入比较',GRAY)},'Uhlitz; cached previous coordinates; CNA/CNN author inferred labels; CNN is not proven nonmalignant')
print('COAD UMAP done',flush=True)
# PRAD: original coordinates and cached expression, exact identifier join.
hp=C/'PRAD/B/20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad';vp=C/'PRAD/B/20260928T030756Z_lypla1_umap_v1/private/plot_values.tsv.gz';src(hp);src(vp);v=pd.read_csv(vp,sep='\t').set_index('cell')
with h5py.File(hp) as h:
 o=h['obs'];ids=hcol(o,o.attrs.get('_index','_index'));v=v.loc[ids];assert np.allclose(v[['UMAP1','UMAP2']],h['obsm/X_umap'][:]);typ=hcol(o,'type');major=hcol(o,'celltype_major_v2');mal=hcol(o,'malignant_anno_merged');keep=(typ=='cancer')&(major=='Epithelial')&np.isin(mal,['normal','malignant'])
umap('PRAD',v[['UMAP1','UMAP2']].to_numpy()[keep],v.LYPLA1_log1p_cp10k.to_numpy()[keep],mal[keep],{'malignant':('恶性上皮',RED),'normal':('非恶性上皮',BLUE)},'PRAD24; author embedding; tumor tissue only; altered_benign not included')
print('PRAD UMAP done',flush=True)
# PDAC: raw author embedding + frozen cell values; untreated epithelial categories.
hp=B/'data/candidates/SCP1089/GSE202051_totaldata-final-toshare.h5ad';vp=C/'PDAC/B/20260928T031629Z_gse202051_lypla1_v1/private/cell_values.tsv.gz';src(hp);src(vp);d=pd.read_csv(vp,sep='\t');d=d[(d.treatment=='Untreated')&d.level2.isin(['Malignant','Ductal','Ductal (atypical)'])].copy()
with h5py.File(hp) as h:
 o=h['obs'];ids=hcol(o,o.attrs.get('_index','_index'));ix=pd.Index(ids).get_indexer(d.cell);assert (ix>=0).all();xy=h['obsm/X_umap'][:][ix]
umap('PDAC',xy,d.expression.to_numpy(),d.level2.to_numpy(),{'Malignant':('恶性上皮',RED),'Ductal':('非恶性导管',BLUE),'Ductal (atypical)':('非典型导管',TEAL)},'GSE202051 untreated nuclei; author X_umap; frozen expression X log1p(CP10K)')
print('PDAC UMAP done',flush=True)

def spatial(name,img,xy,expr,regions,cols,legend,poly=None):
 fig,axs=plt.subplots(1,3,figsize=(12,4));fig.subplots_adjust(left=.01,right=.91,bottom=.16,top=.90,wspace=.025)
 for ax in axs:ax.imshow(img);border(ax)
 axs[0].set_title('H&E',loc='left');axs[1].set_title('作者病理区域',loc='left');axs[2].set_title('LYPLA1 实测表达',loc='left')
 for g,c in cols.items():
  m=regions==g
  if m.any():axs[1].scatter(xy[m,0],xy[m,1],s=3.5,c=c,lw=0,alpha=.82,rasterized=True)
 if poly:
  for q in poly:
   if q['label']=='Tumor':
    z=np.asarray(q['xy']);axs[1].plot(z[:,0],z[:,1],c=RED,lw=.55)
 ix=np.random.default_rng(SEED).permutation(len(xy));h=axs[2].scatter(xy[ix,0],xy[ix,1],c=expr[ix],s=3.5,cmap=CM,vmin=0,vmax=CAP,lw=0,alpha=.90,rasterized=True)
 cb=fig.add_axes([.935,.25,.012,.44]);fig.colorbar(h,cax=cb,extend='max',label='log1p(CP10K)')
 handles=[Line2D([0],[0],marker='o',ls='',color=cols[g],markersize=5,label=lab) for g,lab in legend.items() if (regions==g).any()]
 fig.legend(handles=handles,loc='lower center',ncol=4,frameon=False,fontsize=9)
 save(fig,name);checks.append(dict(plot=name,spots=len(xy),expression_cap=CAP,clipped_fraction=float(np.mean(expr>CAP)),note='Fixed cap for display only; author registration and existing QC preserved; no interpolation'))
# CTA all 14; select examples by patient then lexical sample order, not expression.
r=C/'BRCA/A/20260928T060000Z_cta_lypla1_v1';dat=B/'data/candidates/BRCA_CTA2025';selection=[]
for hp in sorted(r.glob('*_spots_SERVER_ONLY.tsv.gz')):
 src(hp);d=pd.read_csv(hp,sep='\t');d=d[d.qc.astype(str).str.lower().isin(['true','1'])];s=d['sample'].iloc[0];patient=d.patient.iloc[0]
 ip=dat/'spaceranger_output'/s/'outs/spatial/tissue_lowres_image.png';src(ip);img=Image.open(ip)
 regs=list((r/'registration').glob(s+'__*.json'));poly=None
 # Spot regions are the frozen result of reciprocal registration, no registration rerun.
 spatial('spatial_BRCA_'+s,img,d[['x','y']].to_numpy(),d.lognorm.to_numpy(),d.region.to_numpy(),{'Tumor':RED,'Immune':'#B79A34','DCIS':'#836BB1','Unannotated':GRAY,'Mixed':'#EEF0F3','Vessel':TEAL,'Necrosis':'#777777'},{'Tumor':'癌区','Immune':'免疫区','DCIS':'原位癌','Unannotated':'未标注混合区'})
 selection.append(dict(cancer='BRCA',sample=s,patient=patient,plot='spatial_BRCA_'+s))
print('CTA 14 sections done',flush=True)
# PRAD all seven source sections, same panel order and cap.
r=C/'PRAD/B/20260928T101700Z_lypla1_erickson_v1';hp=r/'private/spot_values.tsv.gz';src(hp);d=pd.read_csv(hp,sep='\t')
for s,x in d.groupby('section',sort=True):
 folder=r/'source/Patient1'/s;sf=json.loads((folder/'scalefactors_json.json').read_text());ip=folder/(s+'_tissue_hires_image.png');src(ip);img=Image.open(ip)
 xy=x[['pixel_col','pixel_row']].to_numpy()*sf['tissue_hires_scalef']
 spatial('spatial_PRAD_'+s,img,xy,x.expression.to_numpy(),x.region.to_numpy(),{'Cancer':RED,'Benign':BLUE,'Stroma':GRAY,'PIN':'#B79A34','Transition_State':TEAL},{'Cancer':'癌区','Benign':'良性腺体','Stroma':'间质','PIN':'PIN','Transition_State':'过渡区'})
 selection.append(dict(cancer='PRAD',sample=s,patient='Patient1',plot='spatial_PRAD_'+s,primary_evaluable=bool(sum(x.region=='Cancer')>=20 and sum(x.region=='Benign')>=20)))
pd.DataFrame(selection).to_csv(P/'spatial_display_inventory.tsv',sep='\t',index=False)
(P/'render_validation.json').write_text(json.dumps(dict(status='PASS',checks=checks,new_statistics=False,new_embedding=False,shared_cap=CAP,seed=SEED,zero_color='#E5E8EC',notes='Contours are decorative smoothed outer envelopes; no cells removed. Standardization units agree but platforms differ; do not infer cross-study absolute expression differences.'),indent=2))
pd.DataFrame(sources).to_csv(P/'server_render_source_manifest.tsv',sep='\t',index=False)
print('All server redraws complete',flush=True)
