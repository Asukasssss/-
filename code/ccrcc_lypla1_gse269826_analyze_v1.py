"""Prespecified patient-paired LYPLA1 analysis; all cell/donor data stay on server."""
import argparse,json,itertools,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def bh(x):
 x=np.asarray(x,float);q=np.full(len(x),np.nan);ok=np.flatnonzero(np.isfinite(x));o=ok[np.argsort(x[ok])];q[o]=np.minimum(1,np.minimum.accumulate((x[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1]);return q

def exact(d):
 d=np.asarray(d);s=np.array(list(itertools.product([-1,1],repeat=len(d))),float)
 return float(np.mean(np.abs(s@d/len(d))>=abs(d.mean())-1e-12))

def main(out):
 pub=out/'public';priv=out/'private';rng=np.random.default_rng(20260928)
 spec=json.loads((out/'source/ccrcc_lypla1_gse269826_spec.json').read_text(encoding='utf-8-sig'));selection=json.loads((out/'source/ccrcc_gse269826_selection.json').read_text(encoding='utf-8-sig'))
 d=pd.read_csv(priv/'cells_analysis.tsv',sep='\t');assert d.cell.is_unique
 # Conditions and lineage require exact author fields resolved and locked before inspecting LYPLA1 effects.
 cancer=d.tissue.eq('ccRCC')&d.annotation.isin(selection['tumor_labels']);normal=d.tissue.eq('Normal')&d.annotation.isin(selection['normal_labels'])
 patients=set(d.loc[cancer,'patient']);normal &= d.patient.isin(patients)
 d['condition']=np.where(cancer,'Tumor',np.where(normal,'Normal','Exclude'))
 save(d.groupby(['tissue','annotation','condition'],observed=True).size().rename('n_cells').reset_index(),pub/'annotation_selection_counts.tsv')
 use=d[d.condition.ne('Exclude')].copy();assert use.patient.notna().all();assert np.isfinite(use.expression).all()
 qc=[]
 for condition,z in use.groupby('condition'):
  qc.append(dict(condition=condition,n_cells=len(z),n_patients=z.patient.nunique(),n_regions=z['sample'].nunique(),median_UMI=z.library.median(),median_nFeature=z.nFeature.median(),median_mito_fraction=z.mito.median(),mito_over_025_fraction=(z.mito>.25).mean(),doublet_flagged=int(z.doublet.ne('No').sum())))
 save(pd.DataFrame(qc),pub/'cell_quality_summary.tsv')
 modes=['primary','pseudobulk','region_equal','min100cells','strict_QC','chr3p_loss'];rows=[];pairrows=[];lorows=[];summary=[]
 for mode in modes:
  x=use.copy();minimum=100 if mode=='min100cells' else 20
  if mode=='strict_QC':x=x[x.nFeature.ge(1000)&x.mito.le(.25)]
  if mode=='chr3p_loss':x=x[x.condition.eq('Normal')|x.chr3p.eq('chr3p_loss')]
  p=x.groupby(['patient','condition']).agg(n_cells=('cell','size'),expression=('expression','mean'),detection=('detected','mean'),counts=('counts','sum'),library=('library','sum'),n_regions=('sample','nunique')).reset_index()
  if mode=='pseudobulk':p['expression']=np.log2(1+1e6*p.counts/p.library)
  if mode=='region_equal':
   reg=x.groupby(['patient','condition','sample']).agg(expression=('expression','mean'),detection=('detected','mean')).groupby(['patient','condition']).mean().reset_index()
   p=p.drop(columns=['expression','detection']).merge(reg,on=['patient','condition'],validate='one_to_one')
  p['mode']=mode;save(p,priv/('patient_profiles_'+mode+'.tsv'));e=p[p.n_cells.ge(minimum)]
  t=e[e.condition.eq('Tumor')].set_index('patient');n=e[e.condition.eq('Normal')].set_index('patient');shared=sorted(set(t.index)&set(n.index));t=t.loc[shared];n=n.loc[shared];delta=t.expression.to_numpy()-n.expression.to_numpy();k=len(shared)
  assert k<=20
  base=dict(cancer='ccRCC',cohort='GSE269826',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=spec['analysis_version'],analysis_type='paired_malignant_epithelium_vs_normal_PT',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit='author_patient',n=k,n_reference=k,effect_type='Tumor-minus-Normal '+('log2(1+CPM)_pseudobulk' if mode=='pseudobulk' else 'mean_log1p_CP10K'),effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family='LYPLA1_primary_single_contrast' if mode=='primary' else 'LYPLA1_five_sensitivities',family_n_evaluable=1 if mode=='primary' else 5,status='NOT_EVALUABLE' if k<3 else 'DONE',reason='fewer_than_three_pairs' if k<3 else 'patient_paired_exploratory',source_id='GSE269826',mode=mode,n_tumor_all_eligible=len(e[e.condition.eq('Tumor')]),n_normal_all_eligible=len(e[e.condition.eq('Normal')]),tumor_cells=int(t.n_cells.sum()),normal_cells=int(n.n_cells.sum()),tumor_regions=int(t.n_regions.sum()),normal_regions=int(n.n_regions.sum()),tumor_mean=t.expression.mean(),normal_mean=n.expression.mean(),tumor_detection=t.detection.mean(),normal_detection=n.detection.mean(),n_up=int((delta>0).sum()),n_down=int((delta<0).sum()))
  if k>=3:
   boot=delta[rng.integers(k,size=(10000,k))].mean(1);base.update(effect=delta.mean(),ci_lower=np.quantile(boot,.025),ci_upper=np.quantile(boot,.975),p_value=exact(delta),median_difference=np.median(delta),wilcoxon_p=stats.wilcoxon(delta,method='exact').pvalue)
   for j,patient in enumerate(shared):
    pairrows.append(dict(mode=mode,patient=patient,tumor=t.loc[patient,'expression'],normal=n.loc[patient,'expression'],delta=delta[j]))
    if mode=='primary':
     v=np.delete(delta,j);lorows.append(dict(omitted_patient=patient,n=len(v),effect=v.mean(),p_value=exact(v)))
  rows.append(base)
 df=pd.DataFrame(rows);df.loc[df['mode'].eq('primary'),'q_value']=df.loc[df['mode'].eq('primary'),'p_value'];df.loc[df['mode'].ne('primary'),'q_value']=bh(df.loc[df['mode'].ne('primary'),'p_value'])
 save(df,pub/'results.tsv');save(pd.DataFrame(pairrows),priv/'paired_measurements.tsv');save(pd.DataFrame(lorows),priv/'leave_one_patient_out.tsv')
 lo=pd.DataFrame(lorows);save(pd.DataFrame([dict(n_analyses=len(lo),effect_min=lo.effect.min(),effect_max=lo.effect.max(),p_min=lo.p_value.min(),p_max=lo.p_value.max(),all_effects_negative=bool(lo.effect.lt(0).all()),all_effects_positive=bool(lo.effect.gt(0).all()))]),pub/'leave_one_patient_out_summary.tsv')
 # Export marker summaries, never individual donor expression.
 markers=[c for c in use if c.startswith('marker_')];ms=[]
 for marker in markers:
  p=use.groupby(['patient','condition'])[marker].mean().reset_index()
  for c,z in p.groupby('condition'):ms.append(dict(gene=marker[7:],condition=c,n_patients=z.patient.nunique(),mean_log1p_CP10K=z[marker].mean()))
 save(pd.DataFrame(ms),pub/'marker_context.tsv')
 a=df.iloc[0];fig,axes=plt.subplots(1,2,figsize=(10,4));axes[0].bar(['Normal PT','Cancer epithelium'],[a.normal_mean,a.tumor_mean],color=['#729eaa','#c46c63']);axes[0].set_ylabel('Patient-equal mean log1p(CP10K)');axes[0].set_title('LYPLA1 | %d paired patients\nExact paired P=%.4g'%(a['n'],a.p_value))
 f=df[df['mode'].ne('pseudobulk')].reset_index(drop=True);yy=np.arange(len(f));axes[1].errorbar(f.effect,yy,xerr=[f.effect-f.ci_lower,f.ci_upper-f.effect],fmt='o',color='#34495e',capsize=3);axes[1].axvline(0,color='#aaaaaa',ls='--');axes[1].set_yticks(yy);axes[1].set_yticklabels(f['mode']);axes[1].invert_yaxis();axes[1].set_xlabel('Cancer minus normal PT\n95% paired bootstrap interval');axes[1].set_title('Prespecified sensitivity analyses')
 for ax in axes:ax.spines[['top','right']].set_visible(False)
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(pub/('LYPLA1_new_cohort.'+ext),dpi=200)
 plt.close(fig)
 spec.update(selection=selection,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__));(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 validation=dict(status='DONE',cell_keys_unique=True,independent_unit='author patient',paired_conditions=True,normal_patients_restricted_to_ccRCC=True,all_stats_recomputed_from_full_RNA_counts=True,primary_exact_sign_combinations=2**int(a['n']),primary_family=1,sensitivity_family=5,source_study_independent_of_GSE159115=True,limitations=['normal-adjacent not healthy donors','author inferred CNV annotation, no new genotype verification','bootstrap CI and exact P are different procedures'])
 (pub/'validation.json').write_text(json.dumps(validation,indent=2));print(df[['mode','n','effect','p_value','q_value','tumor_mean','normal_mean','n_up','n_down']].to_string(index=False),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)
