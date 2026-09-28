"""Pathologist-labelled CRC sections: descriptive LYPLA1 maps, no inferred malignancy."""
import argparse, hashlib, io, json, platform, time, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Circle
from matplotlib.collections import PatchCollection
from PIL import Image
from acquire_lypla1_pathology_v1 import DATA, SAMPLES

ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--commit',required=True);a=ap.parse_args()
assert (a.out/'.running').is_dir()
pub=a.out/'public';pub.mkdir(exist_ok=True)
fig=pub/'figures';fig.mkdir(exist_ok=True)
private=a.out/'private';private.mkdir(exist_ok=True)
spec=dict(version='COAD_LYPLA1_pathology_v1',code_base=a.commit,gene='LYPLA1',gene_id='ENSG00000120992',
 source='Zenodo 7760264; Valdeolivas 2024; DOI 10.1038/s41698-023-00488-4',
 selection='Author Methods selected one section per case for histological grading; fixed before LYPLA1 inspection',
 samples=SAMPLES,unit='Multi-cell Visium spot; seven cases, one section each; CRC includes rectum',
 qc='Source in_tissue=1; pathologist annotation present and not exclude; 500<=Gene Expression UMI<=45000; MT- fraction<=0.5',
 excluded='Unannotated/author-excluded/numeric-QC-failed retained in coverage audit and map as gray context',
 normalization='log1p(10000*LYPLA1 UMI/total Gene Expression UMI); descriptive, not source SCTransform',
 region_summary='Per section, original labels plus explicit broad groups; no merging mixed tumor-stroma into tumor',
 comparison='Tumor versus fibroblastic/desmoplastic stroma, and tumor versus non-neoplastic epithelium separately; each region >=20 admitted spots',
 tests='NONE; no spot-level P/q; no pooled-cell patient inference',imputation=False,smoothing=False,
 image='Original low-resolution H&E and matched full-resolution pixel coordinates scaled by author factor',
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__))
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
colors={'Tumor':'#cf243c','Tumor-stroma mixed':'#f39935','Stroma':'#237ec5','Non-neoplastic epithelium':'#2b9c62','Other annotated tissue':'#a481b1','Excluded / unannotated / QC fail':'#bdbdbd'}
def broad(label):
 s=label.lower().strip()
 if s=='tumor':return 'Tumor'
 if s.startswith('tumor&stroma'):return 'Tumor-stroma mixed'
 if s.startswith('stroma'):return 'Stroma'
 if s=='non neo epithelium':return 'Non-neoplastic epithelium'
 if s=='exclude' or not s:return 'Excluded / unannotated / QC fail'
 return 'Other annotated tissue'
def decode(x):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in x])
def member(z,suffix):
 names=[n for n in z.namelist() if n.endswith(suffix) and '__MACOSX' not in n]
 assert len(names)==1,(suffix,names)
 return names[0]
annzip=zipfile.ZipFile(DATA/'Pathology_SpotAnnotations.zip')
assert hashlib.md5((DATA/'Pathology_SpotAnnotations.zip').read_bytes()).hexdigest()=='92882eb0db7c2cacc5077cb5fc9578fd'
records=[];summ=[];audit=[];manifest=[];maps=[];alllabels=[];contrasts=[]
for alias,name,md5 in SAMPLES:
 for attempt in range(180):
  if (DATA/(name+'.selected_manifest.json')).exists():break
  time.sleep(10)
 else:raise TimeoutError('Source not ready: '+name)
 p=DATA/(name+'.selected.zip')
 proof=json.loads((DATA/(name+'.selected_manifest.json')).read_text())
 assert hashlib.sha256(p.read_bytes()).hexdigest()==proof['selected_zip_sha256']
 assert proof['member_crc_verified'] and not proof['source_whole_md5_verified']
 with zipfile.ZipFile(p) as z:
  for item in proof['members']:
   assert hashlib.sha256(z.read(item['member'])).hexdigest()==item['sha256']
  hname=member(z,'filtered_feature_bc_matrix.h5')
  with h5py.File(io.BytesIO(z.read(hname)),'r') as h:
   g=h['matrix'];x=csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:])
   ids=decode(g['features/id'][:]);genes=decode(g['features/name'][:]);bc=decode(g['barcodes'][:]);types=decode(g['features/feature_type'][:])
  assert len(set(bc))==len(bc)
  ix=np.flatnonzero((genes=='LYPLA1')&(ids=='ENSG00000120992')&(types=='Gene Expression'));assert len(ix)==1
  gx=types=='Gene Expression';total=np.asarray(x[gx].sum(axis=0)).ravel()
  assert np.count_nonzero(gx&np.char.startswith(genes,'MT-'))>0
  assert np.all(x.data>=0) and np.all(np.isfinite(x.data))
  mt=np.asarray(x[gx&np.char.startswith(genes,'MT-')].sum(axis=0)).ravel()
  ly=x[ix[0]].toarray().ravel()
  posname=member(z,'tissue_positions_list.csv')
  pos=pd.read_csv(io.BytesIO(z.read(posname)),header=None,names=['Barcode','in_tissue','array_row','array_col','px_row','px_col'])
  assert pos.Barcode.is_unique
  scale=json.loads(z.read(member(z,'scalefactors_json.json')))
  img=np.asarray(Image.open(io.BytesIO(z.read(member(z,'tissue_lowres_image.png')))))
  an=pd.read_csv(annzip.open(member(annzip,'Pathologist_Annotations_'+name+'.csv'))).fillna('')
  assert list(an.columns)[0]=='Barcode' and an.shape[1]==2 and an.Barcode.is_unique
  an.columns=['Barcode','original_label']
  frame=pd.DataFrame({'Barcode':bc,'total_umi':total,'mt_umi':mt,'LYPLA1_umi':ly}).merge(pos,on='Barcode',how='left',validate='one_to_one').merge(an,on='Barcode',how='left',validate='one_to_one')
  assert frame.px_row.notna().all() and frame.px_col.notna().all()
  frame.original_label=frame.original_label.fillna('')
  frame['mt_fraction']=np.divide(mt,total,out=np.zeros_like(total,dtype=float),where=total>0)
  frame['log1pCP10K']=np.log1p(np.divide(10000*ly.astype(float),total,out=np.zeros_like(total,dtype=float),where=total>0))
  frame['original_broad']=frame.original_label.map(broad)
  frame['admitted']=(frame.in_tissue==1)&frame.total_umi.between(500,45000)&(frame.mt_fraction<=.5)&(frame.original_label!='')&(frame.original_label.str.lower()!='exclude')
  assert np.isfinite(frame.log1pCP10K).all()
  for j in [0,len(frame)//2,len(frame)-1]:
   dense=x[:,j].toarray().ravel()
   expected=np.log1p(10000.*dense[ix[0]]/dense[gx].sum()) if dense[gx].sum()>0 else 0.
   assert abs(frame.iloc[j].log1pCP10K-expected)<1e-12
  frame['display_group']=frame.original_broad.where(frame.admitted,'Excluded / unannotated / QC fail')
  frame['x']=frame.px_col*scale['tissue_lowres_scalef'];frame['y']=frame.px_row*scale['tissue_lowres_scalef']
  frame.to_csv(private/(alias+'_spots.tsv'),sep='\t',index=False)
  assert len(set(an.Barcode)-set(pos.Barcode))==0
  radius=scale['spot_diameter_fullres']*scale['tissue_lowres_scalef']/2
  maps.append((alias,frame,img,radius))
  audit.append(dict(section=alias,matrix_spots=len(frame),annotation_spots=len(an),annotation_outside_filtered_matrix=len(set(an.Barcode)-set(bc)),in_tissue=int((frame.in_tissue==1).sum()),unannotated=int((frame.original_label=='').sum()),author_exclude=int((frame.original_label.str.lower()=='exclude').sum()),numeric_qc_failed=int((~frame.total_umi.between(500,45000)|(frame.mt_fraction>.5)).sum()),admitted_spots=int(frame.admitted.sum()),lypla1_detected=int(((frame.LYPLA1_umi>0)&frame.admitted).sum())))
  for lab,bg in an.groupby('original_label'):
   alllabels.append(dict(section=alias,original_label=lab,broad_group=broad(lab),annotation_spots=len(bg)))
  for level,col in [('original','original_label'),('broad','original_broad')]:
   for label,sub in frame.groupby(col):
    good=sub[sub.admitted];n=len(good)
    summ.append(dict(section=alias,level=level,region=label,total_mapped_spots=len(sub),admitted_spots=n,detected_spots=int((good.LYPLA1_umi>0).sum()),detected_fraction=float((good.LYPLA1_umi>0).mean()) if n else np.nan,mean_log1pCP10K=float(good.log1pCP10K.mean()) if n else np.nan,median_log1pCP10K=float(good.log1pCP10K.median()) if n else np.nan,mean_raw_umi=float(good.LYPLA1_umi.mean()) if n else np.nan,median_total_umi=float(good.total_umi.median()) if n else np.nan,status='DONE' if n else 'NOT_EVALUABLE'))
  for ref in ['Stroma','Non-neoplastic epithelium']:
   t=frame[frame.admitted&(frame.original_broad=='Tumor')];r=frame[frame.admitted&(frame.original_broad==ref)]
   ok=len(t)>=20 and len(r)>=20
   contrasts.append(dict(cancer='COAD',cohort='Valdeolivas2024_CRC',stage_id='06_EXTERNAL',run_id=a.out.name,analysis_version=spec['version'],analysis_type='descriptive_region_comparison',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit='spot within section',n=len(t),n_reference=len(r),effect_type='difference_of_mean_log1pCP10K',effect=float(t.log1pCP10K.mean()-r.log1pCP10K.mean()) if ok else np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family='NONE',family_n_evaluable=0,status='DONE' if ok else 'NOT_EVALUABLE',reason='Descriptive only' if ok else 'Fewer than 20 admitted spots in one or both regions',source_id='Zenodo7760264',section=alias,reference_region=ref))
 manifest.append(dict(file=p.name,url=proof['source_url'],version='7760264',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),official_md5=md5,md5_verified=False,verification='Selected member source CRC32 and SHA256 verified; whole source archive not downloaded'))
 for item in proof['members']:
  manifest.append(dict(file=name+'/'+item['member'],url=proof['source_url'],version='7760264',sha256=item['sha256'],verification='Source ZIP member CRC32='+str(item['crc32'])))
 print('PROCESSED',alias,audit[-1],flush=True)
for fn,rows in [('region_summary.tsv',summ),('coverage.tsv',audit),('annotation_dictionary.tsv',alllabels),('region_comparisons.tsv',contrasts)]:
 pd.DataFrame(rows).to_csv(pub/fn,sep='\t',index=False,na_rep='NA')
pa=DATA/'Pathology_SpotAnnotations.zip';manifest.append(dict(file=pa.name,url='https://zenodo.org/api/records/7760264/files/'+pa.name+'/content',version='7760264',sha256=hashlib.sha256(pa.read_bytes()).hexdigest(),official_md5='92882eb0db7c2cacc5077cb5fc9578fd',md5_verified=True))
pd.DataFrame(manifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
cmap=LinearSegmentedColormap.from_list('LYPLA1',['#f0f0f0','#ffd074','#ef703c','#980019'])
vmax=max(f.loc[f.admitted,'log1pCP10K'].max() for _,f,_,_ in maps);norm=Normalize(0,vmax)
def draw(ax,f,img,r,kind,background=True):
 if background:ax.imshow(img,alpha=.38)
 ff=f[f.in_tissue==1]
 if kind=='pathology':
  col=[colors[v] for v in ff.display_group];ec='none'
 else:
  col=[cmap(norm(v)) if ok else '#bdbdbd' for v,ok in zip(ff.log1pCP10K,ff.admitted)];ec='none'
 patches=[Circle((x,y),r*.92) for x,y in zip(ff.x,ff.y)]
 ax.add_collection(PatchCollection(patches,facecolors=col,edgecolors=ec,alpha=.95,rasterized=True))
 if kind=='expression':
  tumor=ff[ff.admitted&(ff.original_broad=='Tumor')]
  ax.add_collection(PatchCollection([Circle((x,y),r*1.04) for x,y in zip(tumor.x,tumor.y)],facecolors='none',edgecolors='#00a4b5',linewidths=.55,rasterized=True))
 ax.set_xlim(0,img.shape[1]);ax.set_ylim(img.shape[0],0);ax.set_aspect('equal');ax.axis('off')
for alias,f,img,r in maps:
 ff,axs=plt.subplots(1,3,figsize=(15,5.3))
 axs[0].imshow(img);axs[0].axis('off');axs[0].set_title('Original H&E')
 draw(axs[1],f,img,r,'pathology');axs[1].set_title('Author pathologist regions')
 draw(axs[2],f,img,r,'expression');axs[2].set_title('LYPLA1 | cyan rings = tumor spots')
 ff.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axs[2],shrink=.6,label='log1p(CP10K)')
 ff.suptitle(alias+' | exact barcode/coordinate match',fontsize=14)
 ff.legend(handles=[Patch(color=c,label=l) for l,c in colors.items()],loc='lower center',ncol=3,fontsize=9)
 ff.subplots_adjust(bottom=.17,top=.88,wspace=.12)
 for ext in ['png','pdf']:ff.savefig(fig/(alias+'_pathology_LYPLA1.'+ext),dpi=170,bbox_inches='tight')
 plt.close(ff)
ff,axs=plt.subplots(7,2,figsize=(10,25))
for i,(alias,f,img,r) in enumerate(maps):
 draw(axs[i,0],f,img,r,'pathology');draw(axs[i,1],f,img,r,'expression')
 axs[i,0].set_title(alias+' | pathology',fontsize=11);axs[i,1].set_title('LYPLA1 | cyan rings = tumor',fontsize=11)
ff.legend(handles=[Patch(color=c,label=l) for l,c in colors.items()],loc='lower center',ncol=2,fontsize=10)
ff.subplots_adjust(bottom=.05,top=.97,wspace=.08,hspace=.15)
cb=ff.add_axes([.92,.4,.018,.2]);ff.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cb,label='LYPLA1 log1p(CP10K)')
ff.suptitle('CRC: pathologist regions and LYPLA1 on identical coordinates',fontsize=15)
for ext in ['png','pdf']:ff.savefig(fig/('all7_pathology_LYPLA1.'+ext),dpi=150,bbox_inches='tight')
plt.close(ff)
checks=dict(status='PASS',sections=len(maps),unique_cases=7,admitted_spots=sum(v['admitted_spots'] for v in audit),annotation_archive_official_md5_verified=True,selected_matrix_members_source_crc32_and_sha256_verified=True,whole_matrix_archives_md5_verified=False,exact_unique_barcode_joins=True,lypla1_unique_ensembl_match=True,unmatched_annotation_barcodes_outside_positions=0,shared_expression_range=[0,float(vmax)],no_new_tests=True,raw_data_location='server165 only',limitations=['CRC not pure COAD','Pathologist tumor spots contain multiple cells','Descriptive region contrasts, no causal or population inference','S7 mainly non-neoplastic; missing tumor comparison retained','Between-study independence from CAMP not independently certified'])
(pub/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print('DONE',json.dumps(checks),flush=True)
