"""LYPLA1 in all untreated STT0000036 sections; author labels, descriptive only."""
import gzip,json,hashlib,sys,platform
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
DATA=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/coad_lypla1_spatial_alternative_20260928')
RUN=Path(__file__).resolve().parent
OUT=RUN/'public';OUT.mkdir(exist_ok=True);(OUT/'figures').mkdir(exist_ok=True)
COLORS={'Tumor-associated epithelium':'#d6334a','Normal epithelium':'#2379bd','Tumor-stroma boundary':'#e3af28','Stroma':'#9ba7b4','Immune-rich stroma':'#49a582','Muscle':'#9671ad','Fiber/cavity':'#e2e2e2','Other':'#555555'}
def region(r):
 if r.level2 in ['epi_CEA','epi_MKI67']:return 'Tumor-associated epithelium'
 if r.level2=='epi_normal':return 'Normal epithelium'
 return {'boundary':'Tumor-stroma boundary','stroma':'Stroma','immune-infiltrated stroma':'Immune-rich stroma','smooth_muscle':'Muscle','fiber and cavity':'Fiber/cavity'}.get(r.level1,'Other')
def analyze_section(row):
 from types import SimpleNamespace
 row=SimpleNamespace(**row)
 results=[];checks=[];saved={}
 sid=row.id;m=pd.read_csv(DATA/row.file);assert not m.duplicated(['x','y']).any()
 index=pd.MultiIndex.from_frame(m[['x','y']]);tot=np.zeros(len(m));nnz=np.zeros(len(m),dtype=int);gene=np.zeros(len(m));seen=np.zeros(len(m),dtype=bool);max_round_error=0.;genes=set();rows=0
 p=DATA/('STexpression_'+sid+'.gem.gz');assert p.exists(),str(p)
 for d in pd.read_csv(p,sep='\t',chunksize=500000):
  assert list(d.columns)==['geneID','x','y','MIDCounts'];v=d.MIDCounts.to_numpy(dtype=float);assert np.isfinite(v).all() and (v>0).all()
  ix=index.get_indexer(pd.MultiIndex.from_frame(d[['x','y']]));assert (ix>=0).all()
  linear=np.expm1(v);np.add.at(tot,ix,linear);np.add.at(nnz,ix,1)
  inferred=linear*m.nCount_Spatial.to_numpy()[ix]/10000
  max_round_error=max(max_round_error,float(np.max(abs(inferred-np.rint(inferred)))))
  genes.update(d.geneID.unique());rows+=len(d)
  sel=d.geneID.eq('LYPLA1').to_numpy();gix=ix[sel];assert len(np.unique(gix))==len(gix) and not seen[gix].any();gene[gix]=v[sel];seen[gix]=True
 # Prove sparse omission is zero and identify the transformation numerically before display.
 assert 'LYPLA1' in genes
 assert np.array_equal(nnz,m.nFeature_Spatial.to_numpy())
 assert np.max(abs(tot-10000))<0.05
 assert max_round_error<0.05
 m['LYPLA1']=gene;m['region']=m.apply(region,axis=1);saved[sid]=m
 m[['index','x','y','level1','level2','level3','region','LYPLA1']].to_csv(RUN/(sid+'_private.tsv.gz'),sep='\t',index=False)
 checks.append(dict(section=sid,n_bins=len(m),gene_features=len(genes),gene_nonzero=int(seen.sum()),source_entries=rows,nnz_matches_nFeature=True,max_library_sum_deviation=float(max(abs(tot-10000))),max_implied_count_rounding_error=max_round_error))
 for (reg,lev2),d in m.groupby(['region','level2']):
  results.append(dict(cancer='COAD',stage_id='06_EXTERNAL',run_id=RUN.name,cohort='STT0000036',section=sid,mmr=row.Type,gene='LYPLA1',region=reg,author_level2=lev2,n_bins=len(d),n_detected=int((d.LYPLA1>0).sum()),detection_fraction=float((d.LYPLA1>0).mean()),mean_expression=d.LYPLA1.mean(),median_expression=d.LYPLA1.median(),effect_scale='source log1p(CP10K), numerically verified',p_value=None,q_value=None,status='DONE'))
 print('ANALYZED',sid,len(m),int(seen.sum()),flush=True)
 return results,checks,saved

def main():
 cat=pd.read_csv(DATA/'section_catalog.tsv',sep='\t');cat['included']=cat.Treatment.eq('naive');cat['reason']=np.where(cat.included,'Untreated section','Treatment stratum excluded from this batch')
 cat.to_csv(OUT/'section_scope.tsv',sep='\t',index=False)
 results=[];checks=[];saved={}
 from concurrent.futures import ProcessPoolExecutor
 with ProcessPoolExecutor(3) as pool:
  for rr,cc,ss in pool.map(analyze_section,cat[cat.included].to_dict('records')):
   results.extend(rr);checks.extend(cc);saved.update(ss)
 pd.DataFrame(results).to_csv(OUT/'author_region_summary.tsv',sep='\t',index=False,na_rep='NA')
 coarse=[]
 for sid,m in saved.items():
  for reg,d in m.groupby('region'):
   coarse.append(dict(section=sid,region=reg,n_bins=len(d),n_detected=int((d.LYPLA1>0).sum()),detection_fraction=(d.LYPLA1>0).mean(),mean_expression=d.LYPLA1.mean(),median_expression=d.LYPLA1.median()))
 coarse=pd.DataFrame(coarse);coarse.to_csv(OUT/'region_summary.tsv',sep='\t',index=False)
 vmax=max(m.LYPLA1.max()for m in saved.values());cmap=LinearSegmentedColormap.from_list('lypla1',['#eeeeee','#ffcc80','#ed7038','#a60026'])
 for sid,m in saved.items():
  fig,axs=plt.subplots(1,2,figsize=(12,6.7));fig.suptitle(sid+' | untreated CRC | LYPLA1',fontsize=17)
  for reg,c in COLORS.items():
   d=m[m.region.eq(reg)];axs[0].scatter(d.x,d.y,c=c,s=5,marker='s',linewidths=0,rasterized=True)
  z=axs[1].scatter(m.x,m.y,c=m.LYPLA1,cmap=cmap,vmin=0,vmax=vmax,s=5,marker='s',linewidths=0,rasterized=True)
  # Author tumor-associated epithelial boundary; no expression-derived ROI.
  xs=np.arange(m.x.min()-1,m.x.max()+2);ys=np.arange(m.y.min()-1,m.y.max()+2);mask=np.zeros((len(ys),len(xs)))
  t=m[m.region.eq('Tumor-associated epithelium')];mask[(t.y-ys[0]).astype(int),(t.x-xs[0]).astype(int)]=1
  if len(t):axs[1].contour(xs,ys,mask,levels=[.5],colors=['#111111'],linewidths=.35)
  axs[0].set_title('1. Find the tissue regions');axs[1].set_title('2. Read LYPLA1 in the same locations')
  for ax in axs:ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off')
  fig.colorbar(z,ax=axs[1],fraction=.032,pad=.02,label='Source log1p(CP10K), shared scale')
  present=[r for r in COLORS if r in set(m.region)];fig.legend(handles=[Patch(color=COLORS[r],label=r)for r in present],loc='lower center',ncol=3,fontsize=9,bbox_to_anchor=(.5,.055))
  fig.text(.5,.018,'Black outline = author tumor-associated epithelium. Red expression is NOT a cancer label. Each point is a 50-um bin.',ha='center',fontsize=9)
  fig.subplots_adjust(left=.015,right=.965,bottom=.22,top=.88,wspace=.08)
  for ext in ['png','pdf']:fig.savefig(OUT/'figures'/(sid+'_LYPLA1.'+ext),dpi=190,bbox_inches='tight',pad_inches=.15)
  plt.close(fig)
 pd.DataFrame(checks).to_csv(OUT/'scale_and_join_checks.tsv',sep='\t',index=False)
 manifest=json.load(open(DATA/'stereoseq_download_manifest.json'));pd.DataFrame(manifest).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
 spec=dict(version='COAD_LYPLA1_Stereo_seq_v1',gene='LYPLA1',scope='All source Treatment=naive sections; selected before inspecting expression',statistical_unit='Within-section spatial bins, descriptive; potentially repeated CRCP59 sections retained separately; patient independence not inferred',region_rule='epi_CEA and epi_MKI67 = author tumor-associated epithelium; epi_normal = author normal epithelium; retain original level2',expression='Source MIDCounts values retained unchanged; log1p(CP10K) confirmed by inverse transform library totals, integer consistency and nFeature match',tests='None; no P or FDR',color_range=[0,float(vmax)],seed=None,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__))
 json.dump(spec,open(OUT/'analysis_spec.json','w'),indent=2)
 json.dump(dict(status='PASS',sections=len(saved),bins=sum(len(x)for x in saved.values()),all_gene_present=True,explicit_xy_join=True,scale_checks=checks,limitations=['Author expression/histology-supported regions, not independent pathology re-review','CRC; no claim of pure COAD','CRCP59 and CRCP59_T_2 have unresolved patient linkage; no patient-level inference','No histology image supplied in this processed export','No differential expression or causal inference']),open(OUT/'validation.json','w'),indent=2)
 print('DONE',len(saved),sum(len(x)for x in saved.values()),flush=True)
if __name__=='__main__':main()
