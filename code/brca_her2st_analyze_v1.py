"""LYPLA1 comparisons using author pathologist labels, one section per patient."""
from pathlib import Path
import sys,json,hashlib,platform
import numpy as np,pandas as pd,scipy
from scipy import stats
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Image.MAX_IMAGE_PIXELS=200000000
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_Andersson_HER2ST')
OUT=Path(sys.argv[1]);(OUT/'figures').mkdir(parents=True,exist_ok=True)
SAMPLES=['A1','B1','C1','D1','E1','F1','G2','H1']
CONTRASTS={'invasive_vs_connective':'connective tissue','invasive_vs_glands':'breast glands','invasive_vs_in_situ':'cancer in situ'}
spec={'version':'v1','annotation':'author pathologist labels; no marker inference','patient_identity':'author naming A-H; one labeled section each, G2 for G',
 'normalization':'log1p(count/total counts*10000)','QC':'author processed count matrix and selected tissue spots; require total>0; no imputation',
 'main_min_spots_per_region':20,'sensitivity_min_spots':10,'contrasts':CONTRASTS,'statistics':'two-sided exact patient sign test; >=3 evaluable patients for inferential reporting; nominal P, no FDR',
 'sensitivities':['regional sum counts normalization','within-patient OLS log normalized expression ~ invasive indicator + log1p(total counts)'],
 'not_claimed':['pure malignant cells','normal epithelial validation from breast glands label alone','activity or mechanism'],
 'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__}}
(OUT/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
manifest=json.loads((ROOT/'manifest.json').read_text())
for f in manifest['files']:assert hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()==f['sha256']
pd.DataFrame(manifest['files']).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
coverage=[];effects=[];audit=[]
colors={'invasive cancer':'#d73027','cancer in situ':'#fc8d59','breast glands':'#1a9850','connective tissue':'#4575b4','immune infiltrate':'#e6b800','adipose tissue':'#66c2a5','undetermined':'#aaaaaa','unmapped':'#dddddd'}
for s in SAMPLES:
 m=pd.read_csv(ROOT/f'data/ST-cnts/{s}.tsv.gz',sep='\t',index_col=0);assert m.index.is_unique and 'LYPLA1' in m.columns
 assert np.isfinite(m.values).all() and (m.values>=0).all() and np.allclose(m.values,np.round(m.values))
 p=pd.read_csv(ROOT/f'data/ST-spotfiles/{s}_selection.tsv',sep='\t');p.index=p.x.astype(int).astype(str)+'x'+p.y.astype(int).astype(str);assert p.index.is_unique
 raw=pd.read_csv(ROOT/f'data/ST-pat/lbl/{s}_labeled_coordinates.tsv',sep='\t');l=raw.dropna(subset=['x','y','pixel_x','pixel_y']).copy()
 l.index=l.x.round().astype(int).astype(str)+'x'+l.y.round().astype(int).astype(str);assert l.index.is_unique
 j=l.join(p,rsuffix='_s',how='inner');assert len(j)==len(l)
 resid=[]
 for axis in ['x','y']:
  X=np.c_[j['pixel_'+axis],np.ones(len(j))];b=j['pixel_'+axis+'_s'];coef=np.linalg.lstsq(X,b,rcond=None)[0]
  resid.append(float(np.max(np.abs(X@coef-b))))
 assert max(resid)<.02
 assert set(m.index).issubset(p.index)
 total=m.sum(axis=1);keep=(total>0)&(p.loc[m.index,'selected'].values==1);m=m.loc[keep];total=total.loc[keep]
 d=p.loc[m.index].copy();d['label']=l.label.reindex(m.index).fillna('unmapped');d['count']=m.LYPLA1;d['total']=total;d['logexpr']=np.log1p(d['count']/total*1e4)
 d.to_csv(OUT/f'{s}_spot_measurements.tsv.gz',sep='\t')
 for label,g in d.groupby('label'):coverage.append({'patient':s[0],'section':s,'label':label,'n_spots':len(g),'detected':int((g['count']>0).sum()),'mean_logexpr':g.logexpr.mean()})
 for contrast,ref in CONTRASTS.items():
  a=d[d.label=='invasive cancer'];b=d[d.label==ref]
  for threshold in [20,10]:
   row={'patient':s[0],'section':s,'contrast':contrast,'threshold':threshold,'n_invasive':len(a),'n_reference':len(b),'delta':np.nan,'pseudobulk_delta':np.nan,'depth_beta':np.nan}
   if min(len(a),len(b))>=threshold:
    row['delta']=a.logexpr.mean()-b.logexpr.mean();row['pseudobulk_delta']=np.log1p(a['count'].sum()/a.total.sum()*1e4)-np.log1p(b['count'].sum()/b.total.sum()*1e4)
    ab=pd.concat([a,b]);X=np.c_[np.ones(len(ab)),(ab.label=='invasive cancer').astype(float),np.log1p(ab.total)]
    row['depth_beta']=np.linalg.lstsq(X,ab.logexpr,rcond=None)[0][1]
   effects.append(row)
 image=Image.open(next((ROOT/f'data/ST-imgs/{s[0]}/{s}').glob('*.jpg')));w,h=image.size
 assert ((d.pixel_x>=0)&(d.pixel_x<w)&(d.pixel_y>=0)&(d.pixel_y<h)).all()
 image.thumbnail((1800,1800));scale=image.width/w;xx=d.pixel_x*scale;yy=d.pixel_y*scale
 fig,axs=plt.subplots(1,3,figsize=(15,5))
 for ax in axs:ax.imshow(image);ax.axis('off')
 axs[0].set_title('H&E');axs[1].set_title('Author pathologist regions')
 for lab,c in colors.items():
  k=d.label==lab
  if k.any():axs[1].scatter(xx[k],yy[k],c=c,s=16,label=lab,alpha=.85,linewidths=0)
 axs[1].legend(loc='upper center',bbox_to_anchor=(.5,-.02),fontsize=7,ncol=2)
 sc=axs[2].scatter(xx,yy,c=d.logexpr,s=16,cmap='magma',vmin=0,vmax=2,linewidths=0);axs[2].set_title('LYPLA1 log1p(CP10K), common scale')
 fig.colorbar(sc,ax=axs[2],shrink=.6,extend='max');fig.suptitle(f'HER2ST | Patient {s[0]} | Section {s} | {len(d)} spots')
 fig.savefig(OUT/'figures'/f'{s}.png',dpi=150,bbox_inches='tight');plt.close(fig)
 audit.append({'section':s,'matrix_spots':len(m),'label_rows':len(raw),'label_rows_missing_coordinates':len(raw)-len(l),'unmapped_matrix_spots':int((d.label=='unmapped').sum()),'max_coordinate_affine_residual_px':max(resid),'detected_spots':int((d['count']>0).sum())})
 print(s,d.label.value_counts().to_dict(),flush=True)
cov=pd.DataFrame(coverage);cov.to_csv(OUT/'patient_region_coverage.tsv',sep='\t',index=False)
eff=pd.DataFrame(effects);eff.to_csv(OUT/'patient_effects.tsv',sep='\t',index=False)
results=[]
for (contrast,threshold),g in eff.groupby(['contrast','threshold']):
 for measure in ['delta','pseudobulk_delta','depth_beta']:
  a=g[measure].dropna();nz=a[a!=0];k=int((nz>0).sum());n=len(a)
  results.append({'cancer':'BRCA','cohort':'Andersson_HER2ST','stage_id':'06_EXTERNAL','run_id':OUT.name,'analysis_version':'v1','analysis_type':contrast,'metabolite_key':'NA','metabolite_name':'NA','gene':'LYPLA1','unit':'patient','n':n,'n_reference':n,'effect_type':measure,'effect':a.mean() if n else np.nan,'ci_lower':'NA','ci_upper':'NA','p_value':stats.binomtest(k,len(nz)).pvalue if n>=3 and len(nz) else np.nan,'q_value':'NA','test_family':'3 pathology contrasts plus 3 measures x 2 thresholds; nominal exploratory P','family_n_evaluable':'NA','status':'DONE' if n>=3 else 'NOT_EVALUABLE','reason':'' if n>=3 else 'fewer than 3 eligible patients; descriptive effect only','source_id':'zenodo.4751624;almaan/her2st','min_spots':threshold,'positive':k,'negative':int((a<0).sum())})
res=pd.DataFrame(results);res['family_n_evaluable']=int((res.status=='DONE').sum());res.to_csv(OUT/'results.tsv',sep='\t',index=False)
cov.groupby('label').agg(patients=('patient','nunique'),spots=('n_spots','sum')).to_csv(OUT/'region_coverage_summary.tsv',sep='\t')
validation={'status':'PASS','patients':8,'sections':8,'spots':int(sum(x['matrix_spots'] for x in audit)),'detected_spots':int(sum(x['detected_spots'] for x in audit)),'unmapped_spots':int(sum(x['unmapped_matrix_spots'] for x in audit)),'missing_coordinate_label_rows':int(sum(x['label_rows_missing_coordinates'] for x in audit)),'input_files_verified':len(manifest['files']),'git_tree':manifest['commit_tree'],'annotation':'original pathologist labels; ambiguous/unmapped excluded from contrasts','coordinates':'rounded array coordinates unique; author pixel coordinate systems related by affine transform residual <0.02 px; all points within image','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'validation.json').write_text(json.dumps(validation,indent=2));pd.DataFrame(audit).to_csv(OUT/'patient_audit.tsv',sep='\t',index=False)
fig,axs=plt.subplots(1,3,figsize=(13,4))
for ax,(contrast,ref) in zip(axs,CONTRASTS.items()):
 g=eff[(eff.contrast==contrast)&(eff.threshold==20)].dropna(subset=['delta']);x=np.arange(len(g));ax.bar(x-.18,g.delta,.36,label='Mean log difference');ax.bar(x+.18,g.depth_beta,.36,label='Depth-adjusted beta');ax.set_xticks(x,g.patient);ax.axhline(0,color='black',lw=.8);ax.set_title('Invasive vs '+ref,fontsize=10);ax.set_xlabel('Patient');ax.legend(fontsize=7)
fig.suptitle('Pathologist-defined regions | Minimum 20 spots per region');fig.tight_layout();fig.savefig(OUT/'figures'/'summary.png',dpi=160);plt.close(fig)
print(json.dumps(validation),flush=True)
