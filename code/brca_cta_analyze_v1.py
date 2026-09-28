"""LYPLA1 in morphology-annotated CTA regions; all per-spot/patient values stay server-side."""
import json,pathlib,hashlib,sys,platform
import numpy as np,pandas as pd,h5py,scipy
from scipy.sparse import csc_matrix
from scipy.stats import binomtest
from matplotlib.path import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from brca_cta_register_v1 import D,O
RUN=O.name
COL={'Tumor':'#cf3232','Immune':'#d8ad22','DCIS':'#3030a0','Vessel':'#3fcbea','Necrosis':'#333333','Unannotated':'#bbbbbb','Mixed':'#dddddd'}
def main():
 (O/'figures').mkdir(exist_ok=True);(O/'public').mkdir(exist_ok=True)
 checks=[]
 for x in json.loads((D/'manifest.json').read_text()):
  p=D/x['file'];h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==x['sha256'];checks.append(x['file'])
 regs=json.loads((O/'registration_summary.json').read_text());regpass=[r for r in regs if r['status']=='PASS']
 # Multiple acceptable SVGs for one expression image would be ambiguous; do not choose by gene expression.
 counts=pd.Series([r['sample'] for r in regpass]).value_counts();regpass=[r for r in regpass if counts[r['sample']]==1]
 coverage=[];effects=[];allspots=[]
 for r in regpass:
  name=r['sample'];d=D/'spaceranger_output'/name/'outs';hp=d/'filtered_feature_bc_matrix.h5'
  if not hp.exists():continue
  with h5py.File(hp) as f:
   g=f['matrix'];genes=np.array([s.decode() for s in g['features/name'][:]]);ids=np.array([s.decode() for s in g['features/id'][:]]);bar=[s.decode() for s in g['barcodes'][:]]
   mat=csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))
  ix=np.where((genes=='LYPLA1')|(ids=='ENSG00000120992'))[0];assert len(ix)==1
  total=np.asarray(mat.sum(axis=0)).ravel();ng=np.diff(mat.indptr);mt=np.asarray(mat[np.char.startswith(genes,'MT-')].sum(axis=0)).ravel()/np.maximum(total,1)
  ly=mat[ix[0]].toarray().ravel();expr=np.log1p(ly/np.maximum(total,1)*1e4)
  co=pd.read_csv(d/'spatial/tissue_positions_list.csv',header=None,names=['barcode','tissue','ar','ac','py','px']).set_index('barcode').loc[bar]
  sf=json.loads((d/'spatial/scalefactors_json.json').read_text());scale=sf['tissue_lowres_scalef'];xy=co[['px','py']].values*scale;radius=sf['spot_diameter_fullres']*scale/2
  poly=json.loads((O/'registration'/f"{name}__{pathlib.Path(r['svg']).stem}.json").read_text())['polygons']
  offs=np.array([(a,b) for a in np.linspace(-1,1,9) for b in np.linspace(-1,1,9) if a*a+b*b<=1])*radius
  pts=(xy[:,None,:]+offs[None,:,:]).reshape(-1,2);masks={}
  for label in ['Tumor','Immune','DCIS','Vessel','Necrosis','Unknown']:
   mask=np.zeros(len(pts),bool)
   for q in poly:
    if q['label']==label:mask |= Path(np.array(q['xy'])).contains_points(pts)
   masks[label]=mask.reshape(len(xy),-1)
  overlap=sum(v.astype(int) for v in masks.values())>1
  unassigned=~np.logical_or.reduce(list(masks.values()));masks['Unannotated']=unassigned
  fractions={k:(v&~overlap).mean(axis=1) for k,v in masks.items()}
  region=np.full(len(xy),'Mixed',dtype=object)
  for label,f in fractions.items():region[f>=.8]=label
  qc=(total>=500)&(ng>=200)&(mt<=.25)&(co.tissue.values==1)
  # Patient source: BC_SA2/3 filenames and embedded BCSA1TumB1.svg on the ST-AR slide.
  patient='BCSA1' if name.startswith('V19T26-012') else 'BCSA2' if 'SA2' in r['svg'] else 'BCSA3' if 'SA3' in r['svg'] else 'UNCONFIRMED'
  frame=pd.DataFrame({'barcode':bar,'sample':name,'patient':patient,'x':xy[:,0],'y':xy[:,1],'total':total,'genes':ng,'mt_fraction':mt,'LYPLA1':ly,'lognorm':expr,'region':region,'qc':qc})
  for label,v in fractions.items():frame[label+'_area_fraction']=v
  frame.to_csv(O/f'{name}_spots_SERVER_ONLY.tsv.gz',sep='\t',index=False);allspots.append(frame[qc])
  coverage.append({'sample':name,'patient':patient,'spots_input':len(bar),'spots_qc':int(qc.sum()),'LYPLA1_positive_qc':int(((ly>0)&qc).sum()),'LYPLA1_counts_qc':int(ly[qc].sum()),'registration_inliers':r['n_inliers'],'registration_error':r['median_error_lowres_px'],**{k:int(((region==k)&qc).sum()) for k in COL}})
  for purity in [.8,1.]:
   for ref in ['Immune','Unannotated']:
    a=qc&(fractions['Tumor']>=purity);b=qc&(fractions[ref]>=purity)
    n=int(a.sum());nr=int(b.sum());ok=n>=20 and nr>=20
    row={'sample':name,'patient':patient,'reference':ref,'purity':purity,'n_tumor':n,'n_reference':nr,'status':'DONE' if ok else 'NOT_EVALUABLE'}
    if ok:
     delta=float(expr[a].mean()-expr[b].mean());ratea=ly[a].sum()/total[a].sum()*1e4;rateb=ly[b].sum()/total[b].sum()*1e4
     sel=a|b;x=np.column_stack([np.ones(sel.sum()),a[sel].astype(float),np.log1p(total[sel]),mt[sel]])
     coef=float(np.linalg.lstsq(x,expr[sel],rcond=None)[0][1])
     row.update(mean_lognorm_delta=delta,pseudobulk_log2ratio=float(np.log2(ratea/rateb)) if ratea>0 and rateb>0 else None,depth_mt_adjusted_delta=coef,tumor_rate_10k=float(ratea),reference_rate_10k=float(rateb))
    effects.append(row)
  image=np.array(Image.open(d/'spatial/tissue_lowres_image.png').convert('RGB'))
  fig,axes=plt.subplots(1,3,figsize=(15,5))
  axes[0].imshow(image)
  for q in poly:
   if q['label'] in COL:
    pts1=np.array(q['xy']);axes[0].plot(*np.r_[pts1,pts1[:1]].T,color=COL[q['label']],lw=.6)
  axes[0].set_title('Author pathology: red tumor, yellow immune')
  axes[1].imshow(image);c=axes[1].scatter(*xy[qc].T,c=expr[qc],s=6,cmap='magma',vmin=0,vmax=max(float(np.quantile(expr[qc],.99)),.1));fig.colorbar(c,ax=axes[1],fraction=.035,label='LYPLA1 log1p(count / total * 10,000)');axes[1].set_title('LYPLA1 expression (QC spots)')
  for k in ['Tumor','Immune','Unannotated','DCIS','Mixed']:
   m=qc&(region==k)
   axes[2].scatter(*xy[m].T,s=7,c=COL[k],label=f'{k} ({m.sum()})')
  axes[2].set_title('Region assignment: >=80% sampled spot area');axes[2].legend(fontsize=7,loc='upper left',bbox_to_anchor=(1,1))
  for ax in axes:ax.set_xlim(0,image.shape[1]);ax.set_ylim(image.shape[0],0);ax.set_aspect('equal');ax.axis('off')
  fig.suptitle(name+' | '+patient);fig.tight_layout();fig.savefig(O/'figures'/f'{name}.png',dpi=150,bbox_inches='tight');plt.close(fig)
  print(name,'QC',int(qc.sum()),'DETECT',int(((ly>0)&qc).sum()),flush=True)
 cov=pd.DataFrame(coverage);eff=pd.DataFrame(effects);cov.to_csv(O/'coverage_SERVER_ONLY.tsv',sep='\t',index=False);eff.to_csv(O/'section_effects_SERVER_ONLY.tsv',sep='\t',index=False)
 agg=[];patients=[]
 for purity in [.8,1.]:
  for ref in ['Immune','Unannotated']:
   for metric in ['mean_lognorm_delta','pseudobulk_log2ratio','depth_mt_adjusted_delta']:
    s=eff[(eff.purity==purity)&(eff.reference==ref)&(eff.status=='DONE')&(eff.patient!='UNCONFIRMED')]
    p=s.groupby('patient')[metric].mean().dropna();nt=int((p!=0).sum());positive=int((p>0).sum());pv=float(binomtest(positive,nt).pvalue) if nt>=3 else None
    agg.append({'reference':ref,'purity':purity,'metric':metric,'patients':len(p),'sections':int(s[metric].notna().sum()),'positive_patients':positive,'mean_patient_effect':float(p.mean()) if len(p) else None,'p_value':pv,'status':'DONE' if len(p) else 'NOT_EVALUABLE'})
    patients.extend({'patient':pid,'reference':ref,'purity':purity,'metric':metric,'effect':v} for pid,v in p.items())
 pd.DataFrame(patients).to_csv(O/'patient_effects_SERVER_ONLY.tsv',sep='\t',index=False)
 pd.DataFrame(agg).to_csv(O/'public/aggregate_effects.tsv',sep='\t',index=False)
 main=[a for a in agg if a['purity']==.8]
 fig,axes=plt.subplots(1,3,figsize=(12,4))
 for ax,metric in zip(axes,['mean_lognorm_delta','pseudobulk_log2ratio','depth_mt_adjusted_delta']):
  for j,ref in enumerate(['Immune','Unannotated']):
   rows=[r for r in patients if r['metric']==metric and r['purity']==.8 and r['reference']==ref]
   vals=[r['effect'] for r in rows];ax.scatter(np.full(len(vals),j)+np.linspace(-.08,.08,len(vals)),vals,c='#bb3333',s=45)
   if vals:ax.plot([j-.15,j+.15],[np.mean(vals)]*2,c='black',lw=2)
  ax.axhline(0,color='gray',ls='--');ax.set_xticks([0,1],['vs immune','vs unannotated']);ax.set_title(metric.replace('_',' '),fontsize=10)
 fig.suptitle('Each dot = one patient; equal weighting of eligible sections');fig.tight_layout();fig.savefig(O/'figures/patient_summary.png',dpi=170);plt.close(fig)
 summary={'matrix_files_available':len(list(D.glob('**/filtered_feature_bc_matrix.h5'))),'manual_svg_available':len(list(D.glob('Images/Manual_annotation/*.svg'))),'sections_analyzed':len(cov),'patients_analyzed':int(cov.patient.nunique()),'spots_input':int(cov.spots_input.sum()),'spots_qc':int(cov.spots_qc.sum()),'LYPLA1_positive_qc':int(cov.LYPLA1_positive_qc.sum()),'normal_epithelium_comparison':'NOT_EVALUABLE','normal_epithelium_reason':'Author manual classes do not include normal epithelium. Unannotated regions are mixed and not normal controls.','aggregate_effects':agg}
 (O/'public/summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
 validation={'source_sha256_verified':len(checks),'unique_gene_hit_each_matrix':True,'registration_used':len(cov),'registration_max_median_error_px':float(cov.registration_error.max()),'patient_id_rule':'Same slide barcode and author SVG embedded BCSA1 identity; SA2/SA3 explicit author filenames; no numeric-order pairing','normal_epithelium':'NOT_EVALUABLE','software':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__},'visual_overlay_check':'PENDING'}
 (O/'public/validation.json').write_text(json.dumps(validation,indent=2))
 print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
