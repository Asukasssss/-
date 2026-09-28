"""Server-only all-candidate expression source, donor-equal log-normalized counts."""
import argparse,json,platform
from pathlib import Path
import numpy as np
import pandas as pd
from cell_origin_v1 import ROOT,sha,load_metadata
from run_cell_expression35_v2 import extract_lee
from patient_statistics import stable_seed

VERSION='COAD_source_contract_sc_v2'
def write(t,p):t.to_csv(p,sep='\t',index=False,na_rep='NA')
def profiles(meta,total,counts,membership,out,study):
 genes=membership.gene.tolist();valid=[g for g in genes if g in counts]
 sel=meta.tissue.eq('Tumor').to_numpy();m=meta.loc[sel].reset_index(drop=True);den=total[sel]
 tumor_counts={g:counts[g][sel] for g in valid}
 assert (den>0).all()
 types=sorted(m.lineage.unique());donors=sorted(m.patient.unique());G=len(genes);D=len(donors);C=len(types)
 expr=np.full((D,C,G),np.nan);det=np.full_like(expr,np.nan);cells=np.zeros((D,C),int)
 private=[];gix={g:genes.index(g) for g in valid}
 for (d,c),ix in m.groupby(['patient','lineage']).groups.items():
  ix=np.array(list(ix));di=donors.index(d);ci=types.index(c);cells[di,ci]=len(ix)
  if len(ix)<20:continue
  for g in valid:
   v=tumor_counts[g][ix];ex=float(np.log1p(v/den[ix]*10000).mean());de=float((v>0).mean());expr[di,ci,gix[g]]=ex;det[di,ci,gix[g]]=de
   private.append(dict(patient=d,celltype=c,gene=g,n_cells=len(ix),mean_expression=ex,mean_detection=de))
 write(pd.DataFrame(private),out/'private'/('donor_'+study+'.tsv'))
 n=(cells>=20).sum(0);eligible=n>=3
 with np.errstate(invalid='ignore',divide='ignore'):
  mean=np.nansum(expr,axis=0)/np.sum(np.isfinite(expr),axis=0);dmean=np.nansum(det,axis=0)/np.sum(np.isfinite(det),axis=0)
 rows=[]
 for ci,c in enumerate(types):
  for gi,g in enumerate(genes):
   good=eligible[ci] and g in valid
   rows.append(dict(cancer='COAD',cohort=study,study=study,stage_id='06_EXTERNAL',run_id=out.name,analysis_version=VERSION,analysis_type='sc_expression_source',partition='ALL',celltype=c,gene=g,human_gene_id=membership.iloc[gi].human_gene_id,n=int(n[ci]),n_donors=int(n[ci]),n_cells=int(cells[cells[:,ci]>=20,ci].sum()),mean_expression=mean[ci,gi] if good else np.nan,mean_detection_fraction=dmean[ci,gi] if good else np.nan,effect=mean[ci,gi] if good else np.nan,effect_type='donor_equal_mean_log1p_CP10K',expression_scale='per_cell_log1p_10000_full_gene_UMI_then_donor_mean',p_value='NA',q_value='NA',status='DONE' if good else 'NOT_EVALUABLE',reason='DESCRIPTIVE_SOURCE_ONLY' if good else 'GENE_MISSING_OR_IDENTITY_HOLD' if g not in valid else 'FEWER_THAN_3_DONORS_WITH_20_CELLS'))
 # Joint donor resampling preserves missing types. A bootstrap needs two eligible observed types.
 rng=np.random.default_rng(stable_seed(VERSION+'|'+study));draws=rng.integers(D,size=(1000,D));freq=np.zeros((C,G));bn=np.zeros(G,int)
 for sample in draws:
  a=expr[sample];denom=np.isfinite(a).sum(0);b=np.divide(np.nansum(a,axis=0),denom,out=np.full((C,G),np.nan),where=denom>0);b[~eligible]=np.nan
  good=np.isfinite(b).sum(0)>=2;bn+=good
  mx=np.max(np.where(np.isfinite(b),b,-np.inf),axis=0);tie=np.isfinite(b)&np.isclose(b,mx[None,:],atol=1e-12,rtol=0);cnt=tie.sum(0)
  freq+=np.divide(tie,cnt[None,:],out=np.zeros_like(b),where=(cnt>0)[None,:])*good[None,:]
 states=[]
 for gi,g in enumerate(genes):
  choices=[ci for ci in range(C) if eligible[ci] and np.isfinite(mean[ci,gi])];choices.sort(key=lambda ci:(-mean[ci,gi],types[ci]))
  top=choices[0] if choices else None;runner=choices[1] if len(choices)>1 else None
  good=runner is not None and dmean[top,gi]>=.01
  both=np.isfinite(expr[:,top,gi])&np.isfinite(expr[:,runner,gi]) if runner is not None else np.zeros(D,bool)
  top_tied=len(choices)>1 and np.isclose(mean[top,gi],mean[runner,gi],atol=1e-12,rtol=0)
  states.append(dict(cancer='COAD',cohort=study,study=study,partition='ALL',gene=g,human_gene_id=membership.iloc[gi].human_gene_id,run_id=out.name,analysis_version=VERSION,top_lineage=types[top] if top is not None else 'NA',runner=types[runner] if runner is not None else 'NA',display_top=types[top] if good else '暂不可定位',top_raw=types[top] if top is not None else 'NA',top_tied=top_tied,max_detection=dmean[top,gi] if top is not None else np.nan,top_expression=mean[top,gi] if top is not None else np.nan,top_runner_gap=mean[top,gi]-mean[runner,gi] if runner is not None else np.nan,bootstrap_top_frequency=freq[top,gi]/bn[gi] if top is not None and bn[gi] else np.nan,valid_bootstrap_n=int(bn[gi]),n_shared_donors=int(both.sum()),top_greater_in_shared_donors=float((expr[both,top,gi]>expr[both,runner,gi]).mean()) if both.sum()>=3 else np.nan,eligible_categories=len(choices),status='DONE' if good else 'NOT_EVALUABLE',reason='DESCRIPTIVE_ONLY;TIE_RETAINED' if good and top_tied else 'DESCRIPTIVE_ONLY' if good else 'GENE_MISSING_OR_IDENTITY_HOLD' if g not in valid else 'INSUFFICIENT_TYPES_OR_TOP_DETECTION_BELOW_1_PERCENT'))
 p=pd.DataFrame(rows);s=pd.DataFrame(states);write(p,out/'public'/(study+'_profiles.tsv'));write(s,out/'public'/(study+'_source_status.tsv'))
 np.savez_compressed(out/'private'/('donor_arrays_'+study+'.npz'),expression=expr,detection=det,cells=cells,genes=np.array(genes),types=np.array(types),donors=np.array(donors))
 return dict(study=study,all_cells=len(meta),tumor_cells=len(m),tumor_donors=D,targets=G,measured=len(valid),source_evaluable=int(s.status.eq('DONE').sum()),stable_top_80=int((s.status.eq('DONE')&(s.bootstrap_top_frequency>=.8)&~s.top_tied).sum()),profile_rows=len(p))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--study',choices=['Lee','Uhlitz'],required=True);a=ap.parse_args();out=a.out
 assert out.parent==ROOT/'results/collaborative/COAD/B' and (out/'.running').exists()
 for d in ['public','private']:(out/d).mkdir(exist_ok=True)
 membership=pd.read_csv(out/'gene_membership.tsv',sep='\t',dtype=str,keep_default_na=False);assert len(membership)==811
 valid=membership[membership.identity_status.eq('DONE')].gene.tolist();assert len(valid)==808
 spec=json.loads((out/'analysis_spec.json').read_text());aliases=spec['confirmed_symbol_aliases']
 if a.study=='Lee':
  d=ROOT/'data/candidates/coad_cell_origin_20260920';meta=load_metadata(d)['Lee2020_GSE132465'];old=out.parent/'20260920T072000Z_cell_expression_v1';check=json.loads((old/'Lee_validation.json').read_text())
  assert sha(d/'lee_raw_UMI.txt.gz')==check['matrix_sha256'];counts,mapping,features=extract_lee(d/'lee_raw_UMI.txt.gz',meta.cell_id.to_numpy(),valid,aliases)
  with np.load(old/'private_Lee_cell_UMI.npz') as z:total=z['total']
  assert len(total)==len(meta) and int(total.sum())==check['all_gene_total_UMI']
  np.savez_compressed(out/'private/Lee_counts.npz',total=total,**counts)
  write(meta,out/'private/Lee_metadata.tsv');write(pd.DataFrame(mapping),out/'public/Lee_feature_coverage.tsv')
  inputs=[d/'lee_raw_UMI.txt.gz',d/'lee_annotation.txt.gz',old/'private_Lee_cell_UMI.npz']
 else:
  meta=pd.read_csv(out/'private_cell_metadata.tsv',sep='\t',keep_default_na=False)
  with np.load(out/'private_cell_counts.npz') as z:
   total=z['total'];counts={k:z[k] for k in z.files if k!='total' and (z[k]>=0).all()}
  inputs=[out/'private_cell_metadata.tsv',out/'private_cell_counts.npz']
 assert all(len(v)==len(total) and (v<=total).all() for v in counts.values())
 summary=profiles(meta,total,counts,membership,out,a.study)
 broad=spec['broad_mapping'][a.study]
 assert set(meta.lineage).issubset(broad)
 meta2=meta.copy();meta2['lineage']=meta2.lineage.map(broad)
 summary['broad_scope']=profiles(meta2,total,counts,membership,out,a.study+'_broad')
 write(pd.DataFrame([dict(study=a.study,original_type=k,broad_type=v,interpretation='Approximate broad compartment; not identical subtypes') for k,v in broad.items()]),out/'public'/(a.study+'_annotation_mapping.tsv'))
 summary.update(status='PASS',p_q_computed=False,source_hashes={str(f):sha(f) for f in inputs},software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),UMAP='NOT_EVALUABLE: no coordinates in admitted metadata; no invented embedding',clinical_subtypes='NOT_EVALUABLE: no verified patient clinical subtype table; cell expression CMS not substituted')
 (out/'public'/(a.study+'_validation.json')).write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
