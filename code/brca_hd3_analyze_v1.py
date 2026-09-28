"""HD 3-prime Kinnex LYPLA1 spatial feasibility; no inferred malignant/normal claims."""
from pathlib import Path
import json,sys,base64,io,re,hashlib,platform
import numpy as np,pandas as pd,scipy
from scipy import io as sio,sparse,stats
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_HD3_Kinnex')
O=Path(sys.argv[1]);(O/'figures').mkdir(parents=True,exist_ok=True);(O/'public').mkdir(exist_ok=True)
SPEC={'version':'hd3_v1','gene':'LYPLA1','sample':'PacBio public humanBreast HD3-Kinnex','donors':1,'histology':'Infiltrating ductal carcinoma per source README',
 'source_matrix':'provider gene-level Kinnex processed count matrix; not short-read Space Ranger matrix',
 'target':'exact LYPLA1 feature only; LYPLA1+TCEA1 excluded from target numerator and reported separately',
 'bins':{'8':{'min_counts':20,'min_genes':10},'16':{'min_counts':50,'min_genes':25},'32':{'min_counts':100,'min_genes':50}},
 'primary_bin_um':16,'tissue_mask':'provider web_summary 8um mask (0=tissue), coarse bins require >=75% covered area',
 'normalization':'log1p(target count / all-feature count * 10000), no imputation',
 'annotation':'Only source tissue mask available. Marker maps and source unsupervised clusters are not pathological tumor/normal annotations',
 'tests':'descriptive within-section Spearman and partial rank association with depth; no inferential P values',
 'primary_question_status':'NOT_EVALUABLE: no verified cancer versus normal epithelial spatial labels',
 'seed':20260928,'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__},
 'markers':['EPCAM','KRT8','KRT18','KRT19','COL1A1','DCN','PTPRC']}
(O/'analysis_spec.json').write_text(json.dumps(SPEC,indent=2))
genes=pd.read_csv(R/'genes.tsv',sep='\t',header=None)[0].to_numpy();gi=np.flatnonzero(genes=='LYPLA1');assert len(gi)==1
bars=pd.read_csv(R/'barcodes.tsv',sep='\t',header=None)[0]
z=bars.str.extract(r'^s_002um_(\d+)_(\d+)-1$');assert z.notna().all().all();rr=z[0].to_numpy(dtype=np.int32);cc=z[1].to_numpy(dtype=np.int32)
m=sio.mmread(R/'matrix.mtx').tocoo();assert m.shape==(len(genes),len(bars));assert (m.data>=0).all() and np.allclose(m.data,np.round(m.data))
m.data=m.data.astype(np.int32);print('MATRIX',m.shape,m.nnz,flush=True)
j=json.JSONDecoder().raw_decode((R/'web_summary.html').read_text().split('const data = ',1)[1])[0]
def image_data(v):
 if v.startswith('_resources_'):v=j['_resources'][v.split('_')[-1]]
 return np.array(Image.open(io.BytesIO(base64.b64decode(v.split(',',1)[1]))))
a=j['tabs']['tab_data'][0]['top']['right']['end_to_end_alignment']['card']['inner']['alignment'];mask=image_data(a['tissueMaskImage'])==0
assert mask.shape==(838,838) and mask.sum()==427917
assert rr.max()<3352 and cc.max()<3352
he=image_data(a['tissueImage']);gray=image_data(a['grayscaleUmiImage']);assert gray.shape[:2]==mask.shape
props=j['tabs']['tab_data'][2]['clustering']['card']['inner']['per_bin_level']['options'][0]['component']['plot']['spatial_plot_props']
T=np.array(props['spot_css_transform']).reshape(4,4,order='F');Tissue=np.array(props['tissue_css_transform']).reshape(4,4,order='F');Tx=np.linalg.inv(Tissue)@T
def xy_hires(row,col,step):
 v=np.vstack([(col+.5)*step,(row+.5)*step,np.zeros(len(row)),np.ones(len(row))]);out=Tx@v;return out[:2]/out[3]
def partial(a,b,depth):
 X=np.c_[np.ones(len(depth)),stats.rankdata(depth)];a=stats.rankdata(a);b=stats.rankdata(b);a-=X@np.linalg.lstsq(X,a,rcond=None)[0];b-=X@np.linalg.lstsq(X,b,rcond=None)[0]
 return float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else None
summary=[];assoc=[];checks={'source_barcode_rows':len(bars),'source_matrix_shape':list(m.shape),'nonzero_entries':m.nnz,'unique_spatial_barcodes':int(bars.nunique()),'mask_tissue_8um':int(mask.sum()),'feature_count':len(genes),'exact_target_rows':len(gi),'normal_epi_labels':'NOT_AVAILABLE'}
amb=np.flatnonzero(genes=='LYPLA1+TCEA1');checks['ambiguous_target_feature_present']=bool(len(amb));checks['ambiguous_target_counts_excluded']=int(m.data[np.isin(m.row,amb)].sum())
checks['exact_LYPLA1_source_counts']=int(m.data[m.row==gi[0]].sum())
for size,param in SPEC['bins'].items():
 size=int(size);factor=size//2;dim=int(np.ceil(3352/factor));binid=(rr//factor)*dim+cc//factor
 b=sparse.coo_matrix((m.data,(m.row,binid[m.col])),shape=(len(genes),dim*dim)).tocsr();b.sum_duplicates();b.eliminate_zeros()
 total=np.asarray(b.sum(axis=0)).ravel();nfeat=np.bincount(b.indices,minlength=dim*dim);cnt=b[gi[0]].toarray().ravel()
 step=size//8
 padded=np.zeros((dim*step,dim*step),bool);padded[:mask.shape[0],:mask.shape[1]]=mask
 tm=padded.reshape(dim,step,dim,step).mean(axis=(1,3)).ravel()>=.75
 good=tm&(total>=param['min_counts'])&(nfeat>=param['min_genes']);y=np.log1p(cnt/np.maximum(total,1)*1e4)
 row={'bin_um':size,'donors':1,'n_tissue_bins':int(tm.sum()),'n_qc_bins':int(good.sum()),'n_detected_bins':int(((cnt>0)&good).sum()),'detected_fraction':float((cnt[good]>0).mean()),'LYPLA1_counts_qc':int(cnt[good].sum()),'median_total_counts_qc':float(np.median(total[good])),'median_features_qc':float(np.median(nfeat[good]))}
 summary.append(row);print('BIN',row,flush=True)
 if size==8:
  intensity=gray if gray.ndim==2 else gray[:,:,:3].mean(axis=2)
  checks['counts_vs_provider_UMI_image_spearman']=float(stats.spearmanr(total[tm],intensity.ravel()[tm]).statistic)
 maps={'LYPLA1':y};rawmaps={'LYPLA1':cnt}
 for gene in SPEC['markers']:
  ii=np.flatnonzero(genes==gene)
  if len(ii)!=1:continue
  c=b[ii[0]].toarray().ravel();v=np.log1p(c/np.maximum(total,1)*1e4);maps[gene]=v;rawmaps[gene]=c
  assoc.append({'bin_um':size,'marker':gene,'n_bins':int(good.sum()),'spearman':float(stats.spearmanr(y[good],v[good]).statistic),'partial_spearman_depth':partial(y[good],v[good],total[good]),'p_value':None,'interpretation':'within-section descriptive only; not independent patients'})
 if all(g in maps for g in ['EPCAM','KRT8','KRT18','KRT19']):
  ep=np.mean([maps[g] for g in ['EPCAM','KRT8','KRT18','KRT19']],axis=0);maps['Epithelial marker mean']=ep
  assoc.append({'bin_um':size,'marker':'Epithelial marker mean','n_bins':int(good.sum()),'spearman':float(stats.spearmanr(y[good],ep[good]).statistic),'partial_spearman_depth':partial(y[good],ep[good],total[good]),'p_value':None,'interpretation':'not a malignant epithelial label'})
 grid=np.arange(dim*dim);br=grid//dim;bc=grid%dim
 pd.DataFrame({'bin_row':br,'bin_col':bc,'tissue':tm,'qc':good,'total_counts':total,'n_features':nfeat,**{g+'_count':v for g,v in rawmaps.items()}}).to_csv(O/f'bin_{size}um_SERVER_ONLY.tsv.gz',sep='\t',index=False)
 fig,axs=plt.subplots(2,3,figsize=(14,10));keys=['LYPLA1','Epithelial marker mean','COL1A1','DCN','PTPRC','Total counts']
 for ax,gene in zip(axs.flat,keys):
  v=np.log1p(total) if gene=='Total counts' else maps.get(gene,np.zeros_like(y));vmax=max(.1,float(np.quantile(v[good],.995)));data=np.where(good,v,np.nan).reshape(dim,dim)
  im=ax.imshow(data,cmap='magma',vmin=0,vmax=vmax,interpolation='nearest');ax.set_title(gene);ax.axis('off');fig.colorbar(im,ax=ax,shrink=.65)
 fig.suptitle(f'Breast IDC | HD 3-prime / Kinnex | {size} um | one specimen\nMarker maps are not cancer/normal annotations');fig.tight_layout();fig.savefig(O/'figures'/f'maps_{size}um.png',dpi=160);plt.close(fig)
 if size==16:
  xy=xy_hires(br,bc,step);checks['transformed_points_in_image_fraction']=float(((xy[0]>=0)&(xy[0]<he.shape[1])&(xy[1]>=0)&(xy[1]<he.shape[0]))[good].mean())
  fig,axs=plt.subplots(1,3,figsize=(15,5))
  for ax in axs:ax.imshow(he);ax.axis('off')
  axs[0].set_title('H&E (provider image)')
  for ax,gene in zip(axs[1:],['LYPLA1','Epithelial marker mean']):
   v=maps[gene];im=ax.scatter(xy[0,good],xy[1,good],c=v[good],s=2,cmap='magma',vmin=0,vmax=max(.1,np.quantile(v[good],.995)),alpha=.8,rasterized=True);ax.set_title(gene+' | 16 um');fig.colorbar(im,ax=ax,shrink=.65)
  fig.suptitle('Exploratory spatial overlay; source CSS homography, no pathology boundaries');fig.tight_layout();fig.savefig(O/'figures'/'he_overlay_16um.png',dpi=170);plt.close(fig)
  # Independent display-space alignment check: provider UMI raster plus our count raster.
  fig,axs=plt.subplots(1,3,figsize=(13,4))
  axs[0].imshow(gray,cmap='gray');axs[0].set_title('Provider 8um UMI image')
  axs[1].imshow(np.where(good,np.log1p(total),np.nan).reshape(dim,dim),cmap='gray');axs[1].set_title('Kinnex 16um counts')
  axs[2].imshow(mask,cmap='gray');axs[2].set_title('Provider tissue mask')
  for ax in axs:ax.axis('off')
  fig.tight_layout();fig.savefig(O/'figures'/'coordinate_qc.png',dpi=130);plt.close(fig)
pd.DataFrame(summary).to_csv(O/'public'/'coverage.tsv',sep='\t',index=False);pd.DataFrame(assoc).to_csv(O/'public'/'marker_associations.tsv',sep='\t',index=False)
(O/'public'/'validation.json').write_text(json.dumps(checks,indent=2));(O/'public'/'analysis_spec.json').write_text(json.dumps(SPEC,indent=2))
manifest=json.loads((R/'manifest.json').read_text());pd.DataFrame(manifest).to_csv(O/'public'/'source_manifest.tsv',sep='\t',index=False)
for f in manifest:assert hashlib.sha256((R/f['file']).read_bytes()).hexdigest()==f['sha256']
print('CHECKS',checks,flush=True)
