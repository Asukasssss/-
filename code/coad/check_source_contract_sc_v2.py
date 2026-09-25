"""Independent aggregate reconstruction, shared-category bootstrap, public delivery checks."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from cell_origin_v1 import sha
from patient_statistics import stable_seed
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);out=p.parse_args().out;pub=out/'public'
 arrays={};checked=0
 for study in ['Lee','Uhlitz','Lee_broad','Uhlitz_broad']:
  with np.load(out/'private'/('donor_arrays_'+study+'.npz')) as z:arrays[study]={k:z[k] for k in z.files}
  a=arrays[study];prof=pd.read_csv(pub/(study+'_profiles.tsv'),sep='\t');state=pd.read_csv(pub/(study+'_source_status.tsv'),sep='\t',keep_default_na=False)
  assert len(state)==811 and state.gene.is_unique and not prof.duplicated(['gene','celltype']).any()
  for row in prof.itertuples():
   if row.status!='DONE':assert pd.isna(row.effect);continue
   ci=list(a['types']).index(row.celltype);gi=list(a['genes']).index(row.gene);v=a['expression'][:,ci,gi];dv=a['detection'][:,ci,gi];v=v[np.isfinite(v)];dv=dv[np.isfinite(dv)]
   assert len(v)==row.n>=3 and np.isclose(sum(v)/len(v),row.effect,atol=1e-12) and np.isclose(sum(dv)/len(dv),row.mean_detection_fraction,atol=1e-12);checked+=1
  assert state.loc[state.status!='DONE','display_top'].eq('暂不可定位').all()
  assert pd.to_numeric(state.loc[state.status=='DONE','max_detection']).ge(.01).all()
 # Compare exactly the same broad categories that satisfy coverage in BOTH studies.
 common=set(arrays['Lee_broad']['types'][(arrays['Lee_broad']['cells']>=20).sum(0)>=3])&set(arrays['Uhlitz_broad']['types'][(arrays['Uhlitz_broad']['cells']>=20).sum(0)>=3]);common=sorted(common);results=[]
 for study in ['Lee','Uhlitz']:
  a=arrays[study+'_broad'];ix=[list(a['types']).index(c) for c in common];v=a['expression'][:,ix,:];det=a['detection'][:,ix,:];D,C,G=v.shape
  den=np.isfinite(v).sum(0);mean=np.divide(np.nansum(v,0),den,out=np.full((C,G),np.nan),where=den>0);dm=np.divide(np.nansum(det,0),den,out=np.full((C,G),np.nan),where=den>0)
  rng=np.random.default_rng(stable_seed('COAD_source_contract_sc_v2|'+study+'|COMMON'));freq=np.zeros((C,G));valid=np.zeros(G,int)
  for _ in range(1000):
   ss=v[rng.integers(D,size=D)];n=np.isfinite(ss).sum(0);b=np.divide(np.nansum(ss,0),n,out=np.full((C,G),np.nan),where=n>0);good=np.isfinite(b).sum(0)>=2;valid+=good
   mx=np.where(np.isfinite(b),b,-np.inf).max(0);ties=np.isfinite(b)&np.isclose(b,mx[None,:],atol=1e-12,rtol=0);cnt=ties.sum(0);freq+=np.divide(ties,cnt[None,:],out=np.zeros_like(b),where=(cnt>0)[None,:])*good
  for gi,g in enumerate(a['genes']):
   order=sorted([ci for ci in range(C) if np.isfinite(mean[ci,gi])],key=lambda ci:(-mean[ci,gi],common[ci]));top=order[0] if order else None;ok=len(order)>=2 and dm[top,gi]>=.01;tie=len(order)>=2 and np.isclose(mean[order[0],gi],mean[order[1],gi],atol=1e-12,rtol=0)
   results.append(dict(gene=g,study=study,common_categories=';'.join(common),common_category_n=C,common_top=common[top] if ok else '暂不可定位',common_raw_top=common[top] if top is not None else 'NA',status='DONE' if ok else 'NOT_EVALUABLE',top_tied=tie,bootstrap_top_frequency=freq[top,gi]/valid[gi] if top is not None and valid[gi] else np.nan,valid_bootstrap_n=int(valid[gi]),top_detection=dm[top,gi] if top is not None else np.nan))
 t=pd.DataFrame(results);t.to_csv(pub/'sc_shared_category_status.tsv',sep='\t',index=False,na_rep='NA');a=t[t.study=='Lee'].set_index('gene');b=t[t.study=='Uhlitz'].set_index('gene');rows=[]
 for g in a.index:
  x,y=a.loc[g],b.loc[g];good=x.status==y.status=='DONE';same=good and x.common_top==y.common_top;stable=same and min(x.bootstrap_top_frequency,y.bootstrap_top_frequency)>=.8 and not(x.top_tied or y.top_tied)
  rows.append(dict(gene=g,Lee_shared_top=x.common_top,Uhlitz_shared_top=y.common_top,Lee_status=x.status,Uhlitz_status=y.status,Lee_top_frequency=x.bootstrap_top_frequency,Uhlitz_top_frequency=y.bootstrap_top_frequency,common_category_n=len(common),both_evaluable=good,same_common_top=same,descriptive_stable_concordance=stable,interpretation='BROAD_COMPARTMENT_ONLY;NOT_IDENTICAL_SUBTYPE_OR_MECHANISM'))
 cross=pd.DataFrame(rows);cross.to_csv(pub/'sc_cross_study.tsv',sep='\t',index=False,na_rep='NA')
 spec=json.loads((out/'analysis_spec.json').read_text());(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 source=[]
 for f in list(out.glob('*.py'))+[out/'analysis_spec.json',out/'gene_membership.tsv']:
  source.append(dict(source_path=str(f),sha256=sha(f)))
 for study in ['Lee','Uhlitz']:
  for path,h in json.loads((pub/(study+'_validation.json')).read_text())['source_hashes'].items():source.append(dict(source_path=path,sha256=h))
 extraction=json.loads((pub/'extraction_validation.json').read_text())
 for f,h in extraction['source_hashes'].items():source.append(dict(source_path=str(Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/coad_uhlitz_20260921')/f),sha256=h))
 pd.DataFrame(source).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 summary=dict(status='PASS',reconstructed_profile_rows=checked,common_categories=common,both_evaluable=int(cross.both_evaluable.sum()),same_common_top=int(cross.same_common_top.sum()),descriptive_stable_concordance=int(cross.descriptive_stable_concordance.sum()),new_P_q=False,limitations=['No raw coordinate file admitted;UMAP unavailable','No verified clinical patient subtypes','Broad categories are approximations,not identical subtypes','Public source labels conditional on coverage;no mechanism'],public_hashes={p.name:sha(p) for p in pub.iterdir() if p.is_file()})
 (pub/'validation.json').write_text(json.dumps(summary,indent=2));(out/'.running').rmdir();(out/'DONE').write_text('DONE\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
