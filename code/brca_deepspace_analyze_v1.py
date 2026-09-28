"""Expert-defined invasive regions versus normal ducts; descriptive section effects.

No patient P values: patient independence is not verified for GSE242311.
No inferred labels. Polygon y coordinates use the documented image-height flip.
"""
from pathlib import Path
import json,gzip,sys,re,hashlib,platform
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as Polygon
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_DeepSpaceDB')
OUT=Path(sys.argv[1]);OUT.mkdir(exist_ok=True,parents=True);(OUT/'figures').mkdir(exist_ok=True)
SPEC={'gene':'LYPLA1','version':'deepspace_v1','selection':'Expert Normal duct and Invasive labels; GEO explicitly One sample; no selection by LYPLA1',
 'included':['DSID000277','DSID000582','DSID000593','DSID000595'],
 'excluded':{'DSID000592':'two specimens; spatial specimen identities not resolved','DSID000594':'two specimens; spatial specimen identities not resolved','DSID000340':'conflicting invasive versus original DCIS pathology, not resolved'},
 'QC':{'min_counts':500,'min_genes':200,'max_mito_fraction':0.25},'normalization':'log1p(count/total_counts*10000)',
 'region_assignment':'spot center in expert polygon, conflicting labels excluded; y=image_height-polygon_y',
 'sensitivity':'24 perimeter points and center within same polygon; radius=spot_diameter_fullres*tissue_hires_scalef/2',
 'minimum_region_spots':20,'secondary_minimum_spots':10,'unit':'section descriptive; no claim of patient independence',
 'P':'NOT_EVALUABLE for patient inference; do not use individual spots as independent biological replicates',
 'normal_reference':'morphologically normal duct region, not purified epithelial cells or healthy donor',
 'GSE210616':'reanalysis of an existing cohort, not independent validation','seed':20260928,
 'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__}}
(OUT/'analysis_spec.json').write_text(json.dumps(SPEC,indent=2))
records=json.loads((ROOT/'expert_annotations.json').read_text())['records']
decode=lambda a:np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
def locate(d,name):
 p=d/name
 return p if p.exists() else Path(str(p)+'.gz')
def fopen(p):return gzip.open(p,'rb') if str(p).endswith('.gz') else open(p,'rb')
rows=[];checks=[];hashes=[];allfig=[]
for rec in records:
 sid=rec['metadata']['ID']
 if sid not in SPEC['included']:continue
 d=ROOT/sid;gsm=rec['metadata']['sample_id'];meta=(d/'metadata.soft').read_text();title=re.search(r'!Sample_title = (.+)',meta).group(1)
 if rec['metadata']['series_id']=='GSE242311':assert 'One sample' in meta
 hp=d/'filtered_feature_bc_matrix.h5'
 if hp.exists():
  with h5py.File(hp) as f:
   g=f['matrix'];m=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))
   genes=decode(g['features']['name'][:]);bars=decode(g['barcodes'][:])
 else:
  m=io.mmread(d/'matrix.mtx.gz').tocsc();m.eliminate_zeros();genes=pd.read_csv(d/'features.tsv.gz',sep='\t',header=None)[1].to_numpy();bars=pd.read_csv(d/'barcodes.tsv.gz',sep='\t',header=None)[0].to_numpy()
 assert len(set(bars))==len(bars)
 pos=pd.read_csv(fopen(locate(d,'spatial/tissue_positions_list.csv')),header=None,names=['barcode','in_tissue','ar','ac','py','px'])
 if pos.iloc[0,0]=='barcode':pos=pos.iloc[1:]
 pos=pos.set_index('barcode');assert pos.index.is_unique and set(bars)<=set(pos.index);pos=pos.loc[bars].astype(float)
 sc=json.load(fopen(locate(d,'spatial/scalefactors_json.json')));img=plt.imread(fopen(locate(d,'spatial/tissue_hires_image.png')),format='png')
 p=json.loads(rec['response']['manual_annotation_plot']);assert (img.shape[1],img.shape[0])==(p['image_width'],p['image_height'])
 xy=np.column_stack([pos.px, pos.py])*sc['tissue_hires_scalef'];radius=sc['spot_diameter_fullres']*sc['tissue_hires_scalef']/2
 total=np.asarray(m.sum(axis=0)).ravel();ng=np.diff(m.indptr);mt=np.asarray(m[np.char.startswith(genes.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(total,1)
 ok=(total>=500)&(ng>=200)&(mt<=.25)&(pos.in_tissue.to_numpy()==1)
 ii=np.flatnonzero(genes=='LYPLA1');assert len(ii)==1;counts=m[ii[0]].toarray().ravel();y=np.log1p(counts/np.maximum(total,1)*1e4)
 paths=[];names=[];center=[];interior=[]
 ang=np.arange(24)*2*np.pi/24
 for t in p['data']:
  pts=np.column_stack([t['x'],p['image_height']-np.array(t['y'])]);path=Polygon(pts);mask=path.contains_points(xy)
  core=mask.copy()
  for a in ang:core &= path.contains_points(xy+radius*np.array([np.cos(a),np.sin(a)]))
  paths.append(pts);names.append(t['name']);center.append(mask);interior.append(core)
 def assign(masks):
  grouped={name:np.any([masks[i] for i,n in enumerate(names) if n==name],axis=0) for name in set(names)}
  lab=np.full(len(bars),'Unannotated',dtype=object);n=np.zeros(len(bars),int)
  for name,mask in grouped.items():lab[mask]=name;n+=mask
  lab[n>1]='Conflicting';return lab
 lab=assign(center);corelab=assign(interior);corelab[lab=='Conflicting']='Conflicting'
 api=d/'spatial_api.json';coord_error=None
 if api.exists():
  j=json.loads(api.read_text());coords=pd.DataFrame(j['coordinates'],columns=['barcode','x','y']).set_index('barcode');common=coords.index.intersection(bars)
  # API uses a fixed 500x500 canvas; compare both axes in its bottom-origin convention.
  ix=pd.Index(bars).get_indexer(common);expected=np.column_stack([xy[ix,0]/img.shape[1]*500,(img.shape[0]-xy[ix,1])/img.shape[0]*500])
  coord_error=float(np.abs(coords.loc[common,['x','y']].astype(float).to_numpy()-expected).max())
  # Record rather than silently changing coordinate transform; H&E overlays checked separately.
 checks.append({'sample':sid,'gene_unique':True,'barcodes_unique':True,'polygon_image_dimensions_match':True,'api_coordinate_max_error':coord_error,'n_qc':int(ok.sum()),'n_conflicting':int(((lab=='Conflicting')&ok).sum())})
 for mode,labels in [('center',lab),('interior',corelab)]:
  a=(labels=='Invasive')&ok;b=(labels=='Normal duct')&ok
  row={'sample':sid,'gsm':gsm,'cohort':rec['metadata']['series_id'],'title':title,'mode':mode,'n_cancer':int(a.sum()),'n_normal':int(b.sum()),'p_value':None,'patient_independence':'UNVERIFIED' if rec['metadata']['series_id']=='GSE242311' else 'one known patient'}
  row['status']='DONE' if min(a.sum(),b.sum())>=20 else 'NOT_EVALUABLE';row['reason']='descriptive section effect only' if row['status']=='DONE' else 'fewer than 20 QC spots in a region'
  if min(a.sum(),b.sum())>0:
   row.update(mean_log_cancer=float(y[a].mean()),mean_log_normal=float(y[b].mean()),delta_log=float(y[a].mean()-y[b].mean()),detection_cancer=float((counts[a]>0).mean()),detection_normal=float((counts[b]>0).mean()))
   ca=counts[a].sum()/total[a].sum()*1e4;cb=counts[b].sum()/total[b].sum()*1e4
   row.update(pseudobulk_cp10k_cancer=float(ca),pseudobulk_cp10k_normal=float(cb),pseudobulk_ratio=float(ca/cb) if cb>0 else None)
   sel=a|b;z=np.log1p(total[sel]);X=np.column_stack([np.ones(sel.sum()),a[sel].astype(float),z-z.mean()]);row['depth_adjusted_delta']=float(np.linalg.lstsq(X,y[sel],rcond=None)[0][1])
  rows.append(row)
 spot=pd.DataFrame({'barcode':bars,'LYPLA1_count':counts,'total_counts':total,'qc':ok,'pathology':lab,'interior_pathology':corelab,'LYPLA1_log1p':y})
 spot.to_csv(OUT/(sid+'_spots.tsv.gz'),sep='\t',index=False)
 fig,axs=plt.subplots(1,3,figsize=(15,5));colors={'Invasive':'#d62728','Normal duct':'#1976d2','Non-invasive':'#e9a323'}
 for ax in axs:ax.imshow(img);ax.axis('off')
 axs[0].set_title('H&E + expert boundaries')
 for pts,name in zip(paths,names):axs[0].plot(pts[:,0],pts[:,1],color=colors.get(name,'#777777'),lw=.8)
 for name,col in colors.items():
  z=(lab==name)&ok;axs[1].scatter(xy[z,0],xy[z,1],s=7,c=col,label=name)
 axs[1].set_title('QC spot centers by expert region');axs[1].legend(fontsize=7)
 sca=axs[2].scatter(xy[ok,0],xy[ok,1],s=7,c=y[ok],cmap='magma',vmin=0,vmax=max(.1,float(np.quantile(y[ok],.99))))
 for pts,name in zip(paths,names):
  if name in colors:axs[2].plot(pts[:,0],pts[:,1],color=colors[name],lw=.6)
 axs[2].set_title('LYPLA1 log1p(CP10K); per-section scale');fig.colorbar(sca,ax=axs[2],shrink=.6)
 fig.suptitle(sid+' | '+gsm+' | '+title);fig.tight_layout();fig.savefig(OUT/'figures'/(sid+'.png'),dpi=170);plt.close(fig)
 for fp in sorted(d.rglob('*')):
  if fp.is_file():hashes.append({'sample':sid,'path':str(fp),'sha256':hashlib.sha256(fp.read_bytes()).hexdigest()})
pd.DataFrame(rows).to_csv(OUT/'section_effects_SERVER_ONLY.tsv',sep='\t',index=False)
(OUT/'validation.json').write_text(json.dumps(checks,indent=2));pd.DataFrame(hashes).to_csv(OUT/'input_hashes.tsv',sep='\t',index=False)
print(pd.DataFrame(rows)[['sample','mode','n_cancer','n_normal','status','delta_log','pseudobulk_ratio','depth_adjusted_delta']].to_string(index=False))
