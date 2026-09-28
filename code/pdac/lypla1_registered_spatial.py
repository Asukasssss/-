"""Measured Visium overlays and epithelial-marker association; no malignant calling."""
import argparse,json,gzip,io,tarfile,hashlib,shutil,platform
from pathlib import Path
import numpy as np,pandas as pd,h5py,scipy
from scipy.sparse import csc_matrix
from scipy.stats import spearmanr,rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib import font_manager
from PIL import Image

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def archive_files(path):
 b=path.read_bytes();layers=0
 while b[:2]==b'\x1f\x8b':b=gzip.decompress(b);layers+=1
 with tarfile.open(fileobj=io.BytesIO(b),mode='r:') as t:
  files={x.name:t.extractfile(x).read() for x in t.getmembers() if x.isfile()}
 return files,layers

def one_file(files,suffixes):
 found=[k for k in files if any(k.endswith(s) for s in suffixes) and not Path(k).name.startswith('._')]
 if len(found)!=1:raise ValueError('Expected unique spatial file: '+str(found))
 return files[found[0]]

def corr(x,y):
 return float(spearmanr(x,y).statistic) if np.std(x)>0 and np.std(y)>0 else np.nan

def partial(x,y,cov):
 a=np.column_stack([np.ones(len(x))]+[rankdata(v) for v in cov]);x=rankdata(x);y=rankdata(y)
 rx=x-a@np.linalg.lstsq(a,x,rcond=None)[0];ry=y-a@np.linalg.lstsq(a,y,rcond=None)[0]
 return float(np.corrcoef(rx,ry)[0,1]) if np.std(rx)>1e-10 and np.std(ry)>1e-10 else np.nan

def run(args):
 out=Path(args.run);pub=out/'public';private=out/'private';spec=json.loads((out/'analysis_spec.json').read_text())
 font_manager.fontManager.addfont(str(private/'msyh.ttc'));plt.rcParams['font.family']=font_manager.FontProperties(fname=str(private/'msyh.ttc')).get_name();plt.rcParams['axes.unicode_minus']=False
 all_data=[];summary=[];effects=[];sources=[];checks=[];groups=spec['reference_markers']
 samples=['GSM8552941_PA02','GSM8552942_PA03','GSM8552944_PA04-2','GSM8552948_PA08','GSM8552949_PA11-1']
 for sample in samples:
  hp=Path(args.data)/(sample+'_filtered_feature_bc_matrix.h5');ap=Path(args.data)/(sample+'_spatial.tar.gz')
  for p in [hp,ap]:
   if not p.exists():raise FileNotFoundError(str(p))
   gsm=sample.split('_')[0];sources.append(dict(file=p.name,url='https://ftp.ncbi.nlm.nih.gov/geo/samples/'+gsm[:-3]+'nnn/'+gsm+'/suppl/'+p.name,bytes=p.stat().st_size,sha256=sha(p)))
  with h5py.File(hp) as f:
   m=f['matrix'];names=m['features/name'][:].astype(str);ids=m['features/id'][:].astype(str);bc=m['barcodes'][:].astype(str)
   mat=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).tocsr()
   assert len(set(bc))==len(bc)
   hit=np.where(names=='LYPLA1')[0];ihit=np.where(ids=='ENSG00000120992')[0];assert len(hit)==1 and np.array_equal(hit,ihit)
   if 'feature_type' in m['features']:assert set(m['features/feature_type'][:].astype(str))=={'Gene Expression'}
  files,layers=archive_files(ap);image=Image.open(io.BytesIO(one_file(files,['tissue_hires_image.png']))).convert('RGB')
  sf=json.loads(one_file(files,['scalefactors_json.json']));position_tables=[]
  for key,pbytes in files.items():
   if Path(key).name not in ['tissue_positions_list.csv','tissue_positions.csv']:continue
   first=pbytes.splitlines()[0].decode();pos=pd.read_csv(io.BytesIO(pbytes),header=0 if 'barcode' in first else None)
   if pos.shape[1]!=6:raise ValueError('Unexpected position schema')
   pos.columns=['barcode','in_tissue','array_row','array_col','pxl_row','pxl_col'];pos=pos.set_index('barcode');assert pos.index.is_unique;position_tables.append(pos.sort_index())
  assert position_tables,'No position file'
  for other in position_tables[1:]:pd.testing.assert_frame_equal(position_tables[0],other,check_dtype=False,check_exact=True)
  pos=position_tables[0]
  assert set(bc)<=set(pos.index);pos=pos.loc[bc];total=np.asarray(mat.sum(axis=0)).ravel();ng=np.asarray((mat>0).sum(axis=0)).ravel()
  mt=np.asarray(mat[np.char.startswith(names,'MT-')].sum(axis=0)).ravel()/np.maximum(total,1)
  keep=(pos.in_tissue.values==1)&(total>0);mat=mat[:,keep];pos=pos.iloc[np.where(keep)[0]];total=total[keep];ng=ng[keep];mt=mt[keep]
  def gene(g):
   idx=np.where(names==g)[0];assert len(idx)==1,g
   return mat[idx[0]].toarray().ravel()
  raw=gene('LYPLA1');expr=np.log1p(raw/total*10000);x=pos.pxl_col.values*sf['tissue_hires_scalef'];y=pos.pxl_row.values*sf['tissue_hires_scalef']
  assert np.all((x>=0)&(x<image.width)&(y>=0)&(y<image.height))
  scores={};used={}
  for label,genes in groups.items():
   z=[];used[label]=[]
   for g in genes:
    a=np.log1p(gene(g)/total*10000)
    if np.std(a)>0:z.append((a-a.mean())/a.std());used[label].append(g)
   assert len(z)>=3;scores[label]=np.mean(z,axis=0)
  frame=pd.DataFrame(dict(barcode=pos.index,x=x,y=y,raw=raw,log1p=expr,total=total,genes=ng,mito_fraction=mt,**scores));frame.to_csv(private/(sample+'_spots.tsv.gz'),sep='\t',index=False)
  for subset,mask in [('author_filtered',np.ones(len(raw),bool)),('qc_sensitivity',(ng>=200)&(total>=500)&(mt<=.2))]:
   n=int(mask.sum())
   if n<20:
    summary.append(dict(sample=sample,subset=subset,n=n,status='NOT_EVALUABLE'));continue
   ep=scores['epithelial'][mask];lo,hi=np.quantile(ep,[.25,.75]);low=ep<=lo;high=ep>=hi
   assert not np.any(low&high),'Quartile tie overlap'
   a=expr[mask];r=raw[mask]
   summary.append(dict(sample=sample,subset=subset,status='DONE',n=n,positive=int((r>0).sum()),detection_pct=100*np.mean(r>0),mean_log1p=a.mean(),low_n=int(low.sum()),high_n=int(high.sum()),low_detection_pct=100*np.mean(r[low]>0),high_detection_pct=100*np.mean(r[high]>0),low_mean=a[low].mean(),high_mean=a[high].mean(),high_minus_low=a[high].mean()-a[low].mean()))
   for label,sc in scores.items():effects.append(dict(sample=sample,subset=subset,reference=label,n=n,spearman_rho=corr(a,sc[mask]),partial_rank_rho=partial(a,sc[mask],[np.log1p(total[mask]),ng[mask],mt[mask]]),p='NA',reason='Descriptive section effect; spots spatially dependent; donor independence not verified'))
  checks.append(dict(sample=sample,target_unique=True,barcodes_unique=True,exact_coordinate_match=True,coordinates_in_bounds=True,gzip_layers=layers,author_filtered_spots=len(bc),retained_in_tissue_nonempty=int(keep.sum()),image_width=image.width,image_height=image.height,marker_genes_used=used))
  # Independent CSC column access checks the target vector without CSR slicing.
  with h5py.File(hp) as f:
   m=f['matrix'];ind=m['indices'][:];ptr=m['indptr'][:];vals=m['data'][:];manual=np.array([vals[ptr[k]:ptr[k+1]][ind[ptr[k]:ptr[k+1]]==hit[0]].sum() for k in range(len(bc))])[keep]
   assert np.array_equal(manual,raw)
  all_data.append(dict(sample=sample,image=image,x=x,y=y,radius=sf['spot_diameter_fullres']*sf['tissue_hires_scalef']/2,expr=expr,raw=raw,scores=scores))
  print(sample,'spots',len(raw),'positive',int((raw>0).sum()),flush=True)
 vmax=max(z['expr'].max() for z in all_data)
 def panel(ax,z,values=None,vmax=None,cmap='magma'):
  ax.imshow(z['image']);ax.set_xlim(0,z['image'].width);ax.set_ylim(z['image'].height,0);ax.axis('off')
  if values is not None:
   patches=[Circle((a,b),z['radius']) for a,b in zip(z['x'],z['y'])]
   pc=PatchCollection(patches,cmap=cmap,linewidth=0,alpha=.75);pc.set_array(values);pc.set_clim(0 if cmap=='magma' else -2,vmax);ax.add_collection(pc);return pc
 for z in all_data:
  fig,axs=plt.subplots(1,5,figsize=(22,5));labels=['H&E 原始切片','LYPLA1 实测表达','上皮参考分数','间质参考分数','免疫参考分数']
  panel(axs[0],z)
  for j,values in enumerate([z['expr'],z['scores']['epithelial'],z['scores']['stromal'],z['scores']['immune']],1):
   pc=panel(axs[j],z,values,vmax if j==1 else 2,'magma' if j==1 else 'coolwarm');fig.colorbar(pc,ax=axs[j],fraction=.035,pad=.01)
  for ax,title in zip(axs,labels):ax.set_title(title,fontsize=13)
  fig.suptitle('GSE278687 | '+z['sample']+' | 原始组织像素坐标配准',fontsize=17)
  fig.text(.5,.035,'LYPLA1：log1p(每万计数)，含零值；参考分数不含 LYPLA1，不能代替恶性标签。实测点，无插值。',ha='center',fontsize=12)
  fig.tight_layout(rect=[0,.07,1,.93]);fig.savefig(pub/(z['sample']+'_registered.png'),dpi=170);fig.savefig(pub/(z['sample']+'_registered.pdf'));plt.close(fig)
 fig,axs=plt.subplots(len(all_data),3,figsize=(12,4*len(all_data)))
 for i,z in enumerate(all_data):
  for j,title in enumerate(['H&E','LYPLA1：共同标尺','上皮参考分数']):
   pc=panel(axs[i,j],z,None if j==0 else z['expr'] if j==1 else z['scores']['epithelial'],vmax if j==1 else 2,'magma' if j==1 else 'coolwarm');axs[i,j].set_title(z['sample']+' | '+title,fontsize=10)
   if j:fig.colorbar(pc,ax=axs[i,j],fraction=.03,pad=.01)
 fig.suptitle('PDAC LYPLA1：5 张新鲜冷冻切片的真实配准叠图',fontsize=17);fig.tight_layout(rect=[0,.02,1,.98]);fig.savefig(pub/'all_five_registered_sections.png',dpi=150);fig.savefig(pub/'all_five_registered_sections.pdf');plt.close(fig)
 pd.DataFrame(summary).to_csv(pub/'section_summary.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(pub/'reference_associations.tsv',sep='\t',index=False);pd.DataFrame(sources).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 shutil.copyfile(out/'analysis_spec.json',pub/'analysis_spec.json')
 (pub/'validation.json').write_text(json.dumps(dict(status='PARTIAL',completed='Five registered measured-expression sections; epithelial/stromal/immune proxy associations',malignant_enrichment='NOT_EVALUABLE: author spot-linked malignant labels unavailable in downloaded GEO files',independent_target_count_check=True,samples=checks,versions={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},code_sha256=sha(__file__)),indent=2))
 (out/'.running').rmdir();print('COMPLETED_REGISTERED_MAPS',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--run',required=True);run(p.parse_args())
