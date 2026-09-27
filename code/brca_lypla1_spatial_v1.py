"""LYPLA1 Wu2021 Visium: region summaries with patient-level inference.

Source matrices and patient/spot outputs remain on server165. Only public/
contains cohort aggregates suitable for the repository.
"""
import argparse, gzip, hashlib, json, platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import h5py
from scipy.io import mmread
from scipy.stats import binomtest, wilcoxon
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VERSION='LYPLA1_Wu_spatial_v1'
GENES=['LYPLA1','EPCAM','KRT8','KRT18','KRT19','PTPRC','COL1A1']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def op(p):
 with open(p,'rb') as f: magic=f.read(2)
 return gzip.open(p,'rb') if magic==b'\x1f\x8b' else open(p,'rb')
def group(x):
 if pd.isna(x):return 'Excluded_missing'
 if 'Invasive cancer' in x or 'Cancer trapped' in x:return 'Tumor_containing'
 if x=='DCIS':return 'DCIS'
 if x=='Normal duct':return 'Normal_duct'
 if x.startswith('Normal'):return 'Normal_mixed'
 if x=='Stroma':return 'Stroma'
 if x in ['Lymphocytes','TLS']:return 'Lymphoid'
 if x=='Necrosis':return 'Necrosis'
 if x in ['Artefact','Uncertain']:return 'Excluded_'+x
 return 'Other_'+x
def dump(x,p):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
def main():
 a=argparse.ArgumentParser();a.add_argument('--root',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--commit',required=True);args=a.parse_args()
 out=args.out; pub=out/'public';pub.mkdir(exist_ok=True);(out/'spatial_figures').mkdir(exist_ok=True)
 specs=dict(version=VERSION,gene='LYPLA1',source='https://zenodo.org/records/4739739',code_commit=args.commit,unit='author patientid; one section per patient verified',normalization='log1p(count/total_UMI*10000); region mean across spots, then paired patient differences',primary='Tumor_containing vs Stroma',secondary=['Tumor_containing vs Normal_mixed','Tumor_containing vs Normal_duct','Tumor_containing vs DCIS','Tumor_containing vs Lymphoid'],min_spots=20,sensitivity_min_spots=[10,30],test='two-sided exact sign test on nonzero patient differences; descriptive only if <3 patients',wilcoxon='supplementary two-sided exact signed-rank when no zeros/ties',p_threshold=0.05,FDR='not used for significance per user; one prespecified primary contrast; secondary exploratory',depth_sensitivity='within patient OLS log1pCP10k ~ tumor_indicator + centered log1p(total_UMI), descriptive; no spot-level P values',composition='EPCAM/KRT8/KRT18/KRT19 and PTPRC/COL1A1 measured descriptively; no validated tumor-fraction adjustment',excluded='Artefact, Uncertain, missing Classification from inferential contrasts; retain QC counts',spatial='coordinates joined by barcode, not row order; no imputation/smoothing',bootstrap='patient resampling CI for mean effect, 20000 draws seed 20260927; unstable with small n',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__))
 dump(specs,pub/'analysis_spec.json')
 audits=[]; manifest=[]; summaries=[]; qc=[]; allspots=[]; depths=[]; patients=[]
 for mf in sorted((args.root/'metadata').glob('*_metadata.csv')):
  sample=mf.name.replace('_metadata.csv','');md=args.root/'filtered_count_matrices'/f'{sample}_filtered_count_matrix'
  paths=[mf]+[md/x for x in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']]
  for p in paths:manifest.append(dict(source_id='Wu2021_Zenodo4739739',file_role=p.name,source_path=str(p),sha256=sha(p)))
  d=pd.read_csv(mf).rename(columns={'Unnamed: 0':'barcode'});d.barcode=d.barcode.astype(str)
  with op(md/'features.tsv.gz') as f:features=pd.read_csv(f,sep='\t',header=None).iloc[:,0].astype(str)
  with op(md/'barcodes.tsv.gz') as f:barcodes=pd.read_csv(f,sep='\t',header=None).iloc[:,0].astype(str)
  with op(md/'matrix.mtx.gz') as f:m=mmread(f).tocsr()
  assert m.shape==(len(features),len(barcodes)) and barcodes.is_unique and d.barcode.is_unique
  assert set(barcodes)==set(d.barcode)
  d=d.set_index('barcode').loc[barcodes].reset_index();d=d.rename(columns={d.columns[0]:'barcode'})
  lib=np.asarray(m.sum(0)).ravel();ng=np.asarray((m>0).sum(0)).ravel()
  assert np.all(lib>0) and np.array_equal(lib,d.nCount_RNA) and np.array_equal(ng,d.nFeature_RNA)
  assert np.all(m.data>=0) and np.all(m.data==np.rint(m.data))
  assert d.patientid.nunique()==1 and d.subtype.nunique()==1
  patients.append(str(d.patientid.iloc[0]));d['sample']=sample;d['region']=d.Classification.map(group)
  for g in GENES:
   ii=np.flatnonzero(features.values==g);assert len(ii)==1,(sample,g,len(ii))
   v=m[ii[0]].toarray().ravel();d[g+'_count']=v;d[g]=np.log1p(v/lib*10000)
  d['total_UMI']=lib
  qc.append(dict(sample=sample,patientid=d.patientid.iloc[0],subtype=d.subtype.iloc[0],n_spots=len(d),n_valid=int((~d.region.str.startswith('Excluded')).sum()),detected=int((d.LYPLA1_count>0).sum()),detection_fraction=float((d.LYPLA1_count>0).mean()),median_UMI=float(np.median(lib))))
  for region,z in d.groupby('region'):
   summaries.append(dict(sample=sample,patientid=d.patientid.iloc[0],subtype=d.subtype.iloc[0],region=region,n_spots=len(z),detection_fraction=float((z.LYPLA1_count>0).mean()),mean_log1p_CP10k=float(z.LYPLA1.mean()),median_log1p_CP10k=float(z.LYPLA1.median()),pseudobulk_CP10k=float(z.LYPLA1_count.sum()/z.total_UMI.sum()*10000),median_UMI=float(z.total_UMI.median()),**{g+'_mean':float(z[g].mean()) for g in GENES[1:]}))
  z=d[d.region.isin(['Tumor_containing','Stroma'])]
  if z.region.value_counts().min()>=20 and z.region.nunique()==2:
   x=np.column_stack([np.ones(len(z)),(z.region=='Tumor_containing').astype(float),np.log1p(z.total_UMI)-np.log1p(z.total_UMI).mean()]);b=np.linalg.lstsq(x,z.LYPLA1,rcond=None)[0]
   depths.append(dict(sample=sample,patientid=d.patientid.iloc[0],adjusted_tumor_effect=float(b[1])))
  # Barcode based joins; positions can contain off-tissue spots absent from filtered matrix.
  spatial_paths=list((args.root/'spatial').rglob('*'))
  pos=[p for p in spatial_paths if sample in str(p) and p.is_file() and 'tissue_positions' in p.name]
  imgs=[p for p in spatial_paths if sample in str(p) and p.is_file() and 'tissue_lowres_image' in p.name]
  sf=[p for p in spatial_paths if sample in str(p) and p.is_file() and 'scalefactors_json' in p.name]
  spatial_ok=False
  hpath=Path('/public3/xuzx/Cancer/breast_cancer_full_workflow_20260710/04_single_cell/cellxgene_wu2021_spatial_h5ad')/f'wu2021_visium_{sample}.h5ad'
  if not pos and hpath.exists():
   # Reuse existing same-study CELLxGENE coordinates/images only; raw counts remain Zenodo.
   with h5py.File(hpath) as f:
    ix=f['obs'][f['obs'].attrs.get('_index','_index')].asstr()[:]
    coords=np.asarray(f['obsm']['spatial']);pp=pd.DataFrame({'barcode':ix,'pixel_col':coords[:,0],'pixel_row':coords[:,1]})
    assert pp.barcode.is_unique
    joined=d.merge(pp,on='barcode',how='left',validate='one_to_one');assert joined.pixel_row.notna().all()
    im=np.asarray(f['uns']['spatial'][sample]['images']['hires']);scale=float(f['uns']['spatial'][sample]['scalefactors']['tissue_hires_scalef'][()])
    cl=f['obs']['Classification']
    if isinstance(cl,h5py.Group):
     cats=cl['categories'].asstr()[:];codes=cl['codes'][:];labs=[cats[c] if c>=0 else None for c in codes]
    else:labs=cl.asstr()[:]
    original=pd.Series(labs,index=ix).reindex(d.barcode).fillna('NA').to_numpy()
    assert np.array_equal(original,d.Classification.fillna('NA').to_numpy()), 'H5AD vs Zenodo pathology mismatch'
   xx=joined.pixel_col*scale;yy=joined.pixel_row*scale
   manifest.append(dict(source_id='CELLxGENE_Wu2021_existing_cache_geometry_only',file_role='spatial_h5ad',source_path=str(hpath),sha256=sha(hpath)))
   spatial_ok=True
  elif len(pos)==1 and len(imgs)==1 and len(sf)==1:
   pp=pd.read_csv(pos[0],header=None)
   if str(pp.iloc[0,0])=='barcode':pp=pp.iloc[1:].copy()
   pp.columns=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col'];pp.barcode=pp.barcode.astype(str)
   assert pp.barcode.is_unique
   joined=d.merge(pp,on='barcode',how='left',validate='one_to_one');assert joined.pixel_row.notna().all()
   scale=json.loads(sf[0].read_text())['tissue_lowres_scalef'];xx=joined.pixel_col.astype(float)*scale;yy=joined.pixel_row.astype(float)*scale;im=plt.imread(imgs[0])
   for p in [pos[0],imgs[0],sf[0]]:manifest.append(dict(source_id='Wu2021_Zenodo4739739',file_role=p.name,source_path=str(p),sha256=sha(p)))
   spatial_ok=True
  if spatial_ok:
   fig,ax=plt.subplots(1,2,figsize=(12,6))
   for aa in ax:aa.imshow(im);aa.axis('off')
   ax[0].set_title('Author pathology (mixed spots)');colors={'Tumor_containing':'#cd4e42','Stroma':'#4275aa','Normal_mixed':'#4ba574','Normal_duct':'#16734a','DCIS':'#d49b37','Lymphoid':'#8964b5','Necrosis':'#777777'}
   for reg,col in colors.items():
    mask=joined.region.eq(reg)
    if mask.any():ax[0].scatter(xx[mask],yy[mask],s=5,c=col,label=reg,alpha=.7)
   ax[0].legend(loc='upper left',bbox_to_anchor=(0,-.01),fontsize=7,ncol=2,frameon=False)
   sc=ax[1].scatter(xx,yy,c=joined.LYPLA1,s=5,cmap='magma',vmin=0,alpha=.8);ax[1].set_title('LYPLA1 observed log1p(CP10k)');fig.colorbar(sc,ax=ax[1],shrink=.6)
   fig.suptitle(sample+' | '+str(d.subtype.iloc[0]));fig.tight_layout();fig.savefig(out/'spatial_figures'/f'{sample}_LYPLA1.png',dpi=180);plt.close(fig)
  audits.append(dict(sample=sample,barcodes_unique_and_matched=True,library_sums_match=True,feature_counts_match=True,integer_counts=True,LYPLA1_unique=True,spatial_join=spatial_ok))
  allspots.append(d)
 assert len(patients)==len(set(patients)), 'Repeated patient: aggregate before inference'
 summary=pd.DataFrame(summaries);summary.to_csv(out/'patient_region_summary.tsv',sep='\t',index=False)
 pd.concat(allspots).to_csv(out/'spot_LYPLA1_markers.tsv.gz',sep='\t',index=False)
 pd.DataFrame(qc).to_csv(out/'sample_qc.tsv',sep='\t',index=False);pd.DataFrame(depths).to_csv(out/'patient_depth_adjusted.tsv',sep='\t',index=False)
 dump(audits,out/'sample_validation.json')
 rows=[];contrasts=[];rng=np.random.default_rng(20260927)
 for minimum in [20,10,30]:
  for ref in ['Stroma','Normal_mixed','Normal_duct','DCIS','Lymphoid']:
   s=summary[summary.n_spots>=minimum];w=s.pivot(index='patientid',columns='region',values='mean_log1p_CP10k')
   if ref in w.columns and 'Tumor_containing' in w.columns:w=w[['Tumor_containing',ref]].dropna()
   else:w=pd.DataFrame(columns=['Tumor_containing',ref])
   delta=(w.Tumor_containing-w[ref]).to_numpy();n=len(delta);nz=delta[delta!=0];p=float(binomtest(int((nz>0).sum()),len(nz)).pvalue) if n>=3 and len(nz) else np.nan
   ci=np.quantile(rng.choice(delta,(20000,n),replace=True).mean(1),[.025,.975]) if n>=3 else [np.nan,np.nan]
   row=dict(cancer='BRCA',cohort='Wu2021_Visium',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=VERSION,analysis_type='spatial_region_comparison',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit='patient',n=n,n_reference=n,effect_type='paired_difference_mean_log1p_CP10k',effect=float(delta.mean()) if n else np.nan,ci_lower=ci[0],ci_upper=ci[1],p_value=p,q_value=np.nan,test_family='one_primary_stroma; other_contrasts_exploratory',family_n_evaluable=1,status='DONE' if n>=3 else 'NOT_EVALUABLE',reason='' if n>=3 else 'Fewer than 3 patients with >= minimum spots in both regions',source_id='Zenodo4739739',reference=ref,min_spots=minimum,n_up=int((delta>0).sum()),n_down=int((delta<0).sum()),median_difference=float(np.median(delta)) if n else np.nan,wilcoxon_p=float(wilcoxon(delta,method='exact').pvalue) if n>=3 and len(nz)==n and len(set(abs(delta)))==n else np.nan)
   rows.append(row)
   for pat,dd in zip(w.index,delta):contrasts.append(dict(patientid=pat,reference=ref,min_spots=minimum,difference=dd))
 results=pd.DataFrame(rows);results.to_csv(pub/'results.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(contrasts).to_csv(out/'patient_contrasts.tsv',sep='\t',index=False)
 aggregates=summary.groupby('region').agg(n_patients=('patientid','nunique'),n_spots=('n_spots','sum'),patient_mean_expression=('mean_log1p_CP10k','mean'),patient_mean_detection=('detection_fraction','mean')).reset_index();aggregates.to_csv(pub/'region_aggregate.tsv',sep='\t',index=False)
 q=pd.DataFrame(qc);info=dict(n_patients=len(patients),n_sections=len(patients),n_spots=int(q.n_spots.sum()),n_valid_spots=int(q.n_valid.sum()),LYPLA1_positive_spots=int(q.detected.sum()),overall_detection=float(q.detected.sum()/q.n_spots.sum()),sample_detection_min=float(q.detection_fraction.min()),sample_detection_max=float(q.detection_fraction.max()),subtypes=q.subtype.value_counts().to_dict(),n_spatial_maps=sum(x['spatial_join'] for x in audits),depth_adjusted_n=len(depths),depth_adjusted_n_positive=sum(x['adjusted_tumor_effect']>0 for x in depths),depth_adjusted_mean=float(np.mean([x['adjusted_tumor_effect'] for x in depths])))
 dump(info,pub/'cohort_summary.json')
 pd.DataFrame(manifest).to_csv(out/'source_manifest_private.tsv',sep='\t',index=False)
 publicmanifest=[]
 for name in ['filtered_count_matrices.tar.gz','metadata.tar.gz']:
  p=args.root/name;publicmanifest.append(dict(source_id='Zenodo4739739',url='https://zenodo.org/records/4739739/files/'+name,sha256=sha(p),bytes=p.stat().st_size))
 for row in manifest:
  if row['file_role']=='spatial_h5ad':publicmanifest.append(dict(source_id=row['source_id'],url='same-study existing CELLxGENE spatial cache; original file path retained in server manifest',sha256=row['sha256'],bytes=Path(row['source_path']).stat().st_size))
 pd.DataFrame(publicmanifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 dump(dict(status='PASS',n_samples_checked=len(audits),all_matrix_checks_pass=True,all_spatial_joins_pass=all(x['spatial_join'] for x in audits),patient_ids_unique=True,private_measurements_server_only=True,limitations=['mixed spots are not purified epithelial cells','normal duct has only three spots; no adequate pure-normal comparison','no independent cohort validation','no tumor-fraction or spatial-autocorrelation model; no spot-level P values','small patient sample; bootstrap CI descriptive']),pub/'validation.json')
 rr=results[results.min_spots.eq(20)&results.n.gt(0)].copy();fig,axs=plt.subplots(1,2,figsize=(11,4.5))
 z=aggregates[~aggregates.region.str.startswith('Excluded')];axs[0].barh(z.region,z.patient_mean_expression,color='#4e79a7');axs[0].set_xlabel('Mean across patient-region means: log1p(CP10k)');axs[0].set_title('LYPLA1 by author pathology group')
 axs[1].barh(rr.reference,rr.effect,color=['#c96851' if x>0 else '#4e79a7' for x in rr.effect]);axs[1].axvline(0,color='grey',lw=.8);axs[1].set_xlabel('Tumor-containing minus reference');axs[1].set_title('Within-patient differences (region >=20 spots)')
 for i,(_,r) in enumerate(rr.iterrows()):axs[1].text(r.effect,i,f' n={r.n}, up={r.n_up}',va='center',fontsize=8)
 fig.tight_layout();fig.savefig(pub/'LYPLA1_region_overview.png',dpi=180,bbox_inches='tight');fig.savefig(pub/'LYPLA1_region_overview.pdf',bbox_inches='tight');plt.close(fig)
 print(json.dumps(info));print(results[results.min_spots.eq(20)][['reference','n','effect','n_up','p_value','status']].to_string(index=False))
if __name__=='__main__':main()
