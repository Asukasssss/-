"""Restore Pal clinical labels using exact original GEO titles; source-side only."""
from pathlib import Path
import sys,json,gzip,hashlib,importlib.util,itertools
import requests,h5py,numpy as np,pandas as pd
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
R=B/'20260927T122639Z_pal_subtypes_v1';P=R/'public';D=R/'private';S=R/'source'
def save(x,p):x.to_csv(p,sep='\t',index=False,na_rep='NA')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for z in iter(lambda:f.read(8*1024*1024),b''):h.update(z)
 return h.hexdigest()
def main():
 assert (R/'.running').exists()
 for d in [P,D,S]:d.mkdir(exist_ok=True)
 url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/soft/GSE161529_family.soft.gz'
 res=requests.get(url,timeout=90);res.raise_for_status();sf=S/'GSE161529_family.soft.gz';sf.write_bytes(res.content)
 rows=[]
 for b in gzip.decompress(res.content).decode().split('^SAMPLE = ')[1:]:
  ls=b.splitlines();r={'gsm':ls[0]}
  for l in ls:
   if l.startswith('!Sample_title = '):r['title']=l.split(' = ',1)[1]
   if l.startswith('!Sample_characteristics_ch1 = '):
    k,v=l.split(' = ',1)[1].split(': ',1);r[k]=v
  if 'tumour' in r.get('cancer type','') and r.get('cell population')=='Total':
   r['donor_id']='pal_Patient '+r['title'].split('from Patient ',1)[1];rows.append(r)
 geo=pd.DataFrame(rows);save(geo,S/'GEO_primary_tumor_metadata.tsv')
 hp=B/'20260925T133935Z_lypla1_KEGGscore_v1/brca_sc117_profile_v1.py'
 sp=importlib.util.spec_from_file_location('h',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 src=B/'20260919T140105Z_scRNA117_v1/source/pal_atlas_container.h5ad'
 with h5py.File(src) as f:
  obs=h.dataframe(f['obs']);var=h.dataframe(f['raw/var']);sel=obs.batch.eq('pal_2021');mapping=[]
  for name in obs.loc[sel,'donor_id'].unique():
   g=geo[geo.donor_id.eq(name)];types=g['cancer type'].unique();ok=len(g)==1
   subtype='UNRESOLVED'
   if ok:
    t=types[0];subtype='TNBC' if t.startswith('Triple negative') else ('ER+' if t=='ER+ tumour' else ('HER2+' if t=='HER2+ tumour' else 'PR+_only_label'))
   mapping.append(dict(donor_id=name,subtype=subtype,patient=g.patient.iloc[0] if ok else name,status='DONE' if ok else 'NEEDS_REVIEW',reason='exact_title_join' if ok else 'conflicting_original_GEO_records',gsm=';'.join(g.gsm),original_type=';'.join(types)))
  m=pd.DataFrame(mapping);assert m.donor_id.is_unique;save(m,D/'sample_subtype_mapping.tsv')
  o=obs.loc[sel].copy();o=o.join(m.set_index('donor_id')[['subtype','patient','status']],on='donor_id',validate='many_to_one');assert o.subtype.notna().all()
  save(o[['donor_id','subtype','patient','status']].reset_index(),D/'cell_subtype_sidecar.tsv')
  coverage=o.groupby('subtype',observed=True).agg(n_cells=('donor_id','size'),n_source_labels=('donor_id','nunique'),n_GEO_patient_labels=('patient','nunique')).reset_index();save(coverage,P/'mapping_coverage.tsv')
  mal=(sel & obs.cell_type.eq('malignant cell')).to_numpy();n=len(obs);p=len(var);target=np.flatnonzero(var.feature_name.astype(str).eq('LYPLA1').to_numpy());assert len(target)==1
  raw=f['raw/X'];ptr=raw['indptr'][:];lin=np.full(n,np.nan)
  for st in range(0,n,2000):
   en=min(n,st+2000);mask=mal[st:en]
   if not mask.any():continue
   lo,hi=ptr[st],ptr[en];x=sparse.csr_matrix((raw['data'][lo:hi].astype(float),raw['indices'][lo:hi],ptr[st:en+1]-lo),shape=(en-st,p))[mask];tot=np.asarray(x.sum(1)).ravel();assert (tot>0).all();lin[np.flatnonzero(mask)+st]=x[:,target[0]].toarray().ravel()/tot*10000
  a=obs.loc[mal,['donor_id']].copy().join(m.set_index('donor_id'),on='donor_id',rsuffix='_map');a['linear']=lin[mal];a['log']=np.log1p(lin[mal]);a['detected']=lin[mal]>0;save(a.reset_index(),D/'LYPLA1_malignant_cell_values.tsv')
  d=a[a.status.eq('DONE')].groupby(['subtype','patient'],observed=True).agg(n_cells=('log','size'),mean_log=('log','mean'),mean_CP10k=('linear','mean'),detection=('detected','mean')).reset_index();save(d,D/'LYPLA1_patient_values.tsv');d=d[d.n_cells>=20]
  summary=[]
  for typ,g in a.groupby('subtype',observed=True):
   dd=d[d.subtype.eq(typ)];summary.append(dict(subtype=typ,n_cells=len(g),n_source_labels=g.donor_id.nunique(),n_patients_eligible=len(dd),equal_patient_mean_log=dd.mean_log.mean(),equal_patient_mean_CP10k=dd.mean_CP10k.mean(),equal_patient_detection=dd.detection.mean(),pooled_cell_mean_log=g.log.mean(),pooled_cell_mean_CP10k=g.linear.mean(),pooled_detection=g.detected.mean(),status='NEEDS_REVIEW' if typ=='UNRESOLVED' else 'DONE'))
  su=pd.DataFrame(summary);save(su,P/'LYPLA1_subtype_summary.tsv')
  tests=[]
  for unit in ['cell_exploratory','GEO_patient_label']:
   for x,y in itertools.combinations(['ER+','HER2+','TNBC'],2):
    v=a if unit=='cell_exploratory' else d;col='log' if unit=='cell_exploratory' else 'mean_log';xx=v[v.subtype.eq(x)][col].to_numpy();yy=v[v.subtype.eq(y)][col].to_numpy();u,pv=stats.mannwhitneyu(xx,yy,alternative='two-sided',method='asymptotic');tests.append(dict(cancer='BRCA',cohort='Pal2021_reprocessed',stage_id='06_EXTERNAL',run_id=R.name,analysis_version='pal_subtypes_v1',analysis_type='LYPLA1_subtype_comparison',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit=unit,n=len(xx),n_reference=len(yy),effect_type='rank_biserial_positive_first_group_higher',effect=2*u/(len(xx)*len(yy))-1,ci_lower=np.nan,ci_upper=np.nan,p_value=pv,q_value=np.nan,test_family=unit+'_3_pairs',family_n_evaluable=3,status='DONE',reason='two_sided_MWU_asymptotic_tie_and_continuity_corrected',source_id='GSE161529',group=x,reference=y,mean_difference=xx.mean()-yy.mean()))
  te=pd.DataFrame(tests)
  for unit,ix in te.groupby('unit').groups.items():te.loc[ix,'q_value']=multipletests(te.loc[ix,'p_value'],method='fdr_bh')[1]
  save(te,P/'results.tsv')
  fig,axs=plt.subplots(1,2,figsize=(10,4));order=['ER+','HER2+','TNBC'];colors=['#4676ac','#dd9d42','#b35777']
  for ax,col,title in zip(axs,['mean_log','detection'],['LYPLA1: malignant epithelial expression','LYPLA1: fraction of detected cells']):
   vals=[d[d.subtype.eq(t)][col].to_numpy() for t in order];bp=ax.boxplot(vals,positions=range(3),widths=.45,patch_artist=True,showfliers=False)
   for b,c in zip(bp['boxes'],colors):b.set_facecolor(c);b.set_alpha(.25)
   for i,v in enumerate(vals):ax.scatter(i+np.linspace(-.12,.12,len(v)),v,color=colors[i],s=24)
   ax.set_xticks(range(3),[t+'\nn='+str(len(v)) for t,v in zip(order,vals)]);ax.set_title(title,fontsize=11);ax.spines[['top','right']].set_visible(False)
  axs[0].set_ylabel('Patient mean log1p(CP10k)');axs[1].set_ylabel('Patient detection fraction');fig.suptitle('Pal2021: restored GEO subtype labels');fig.text(.5,.01,'Each point is a GEO patient label; >=20 malignant cells. Conflicting and PR-only labels excluded.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.05,1,.95]);fig.savefig(P/'LYPLA1_Pal_subtypes.png',dpi=180);fig.savefig(P/'LYPLA1_Pal_subtypes.pdf');plt.close(fig)
  assert len(a)==130078
  audit=dict(status='PASS',source_labels=len(m),matched_labels=int(m.status.eq('DONE').sum()),conflicting_labels=int(m.status.ne('DONE').sum()),malignant_cells=len(a),one_to_many_join_rejected=True,patient_identity='GEO patient field; no genotype verification',same_patient_regions_merged=True,old_matrix_unchanged=True,old_statistics_not_overwritten=True,not_checked='Clinical subtype conflict not resolved; tests are exploratory; no subtype-by-cancer-normal comparison')
  (P/'validation.json').write_text(json.dumps(audit,indent=2))
  spec=dict(version='pal_subtypes_v1',source_commit='6f33f047a8a94b9295cd39cced30a4c3507acfad',join='pal_Patient + exact GEO title suffix, primary Total tumor only; unique GSM required',PR_only='separate, not assumed ER+',conflicts='retain unresolved, exclude from subtype tests',patient_unit='GEO patient characteristic; explicit shared patient regions pooled before >=20 cell filter',expression='raw counts/library size*1e4; log1p; patient mean log expression and detection',tests='3 pairwise two-sided asymptotic Mann Whitney per unit, BH3 separate; cells and GEO patient labels',limitations='No donor-adjusted cell test, batch/depth/sex adjustment, or causal interpretation; PR-only n=1 descriptive; sources not new independent cohort',software=dict(numpy=np.__version__,pandas=pd.__version__))
  (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
  save(pd.DataFrame([dict(path=str(q),sha256=sha(q),url=url if q==sf else 'NA') for q in [sf,src,hp,Path(__file__)]]),P/'source_manifest.tsv')
  print(su.to_string(index=False));print(te[['unit','group','reference','effect','p_value']].to_string(index=False));print(audit)
 (R/'NUMERICAL_DONE').write_text('DONE')
if __name__=='__main__':main()
