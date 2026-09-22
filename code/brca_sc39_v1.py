"""Extend frozen single-cell profiles to new39 only; server165 raw data stay resident."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,importlib.util
import pandas as pd,numpy as np,h5py
R=Path(sys.argv[1]);BASE=R.parent;PREV=BASE/'20260919T140105Z_scRNA117_v1';P=R/'public';V='sc39_v1'
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
assert (R/'.running').read_text().strip()==V
genes=pd.read_csv(R/'source/genes_unique.tsv',sep='\t');genes=genes.loc[~genes.in_original117,['gene']];assert len(genes)==39
save(genes,R/'source/candidate_comparison_v1.tsv')
# Parameter adaptation of frozen helper, preserving all normalization/aggregation logic.
text=(R/'brca_sc117_profile_v1.py').read_text();text=text.replace('117','39').replace("'sc39_v1'","'sc39_v1'")
# Resolve the new LPCAT4 Q643R3 candidate by stable Ensembl ID rather than its ambiguous symbol.
text=text.replace('positions = np.flatnonzero(names.values == gene)',"positions = np.flatnonzero(var.index.astype(str).str.split('.').str[0].isin(['ENSG00000176454', 'ENSG00000291994'])) if gene == 'LPCAT4' else np.flatnonzero(names.values == gene)")
(R/'profile39_adapter.py').write_text(text);prof=module('profile39',R/'profile39_adapter.py')
st=(R/'brca_all117_stability_v2.py').read_text().replace('117','39');(R/'stability39_adapter.py').write_text(st);stab=module('stability39',R/'stability39_adapter.py')
spec=dict(version=V,genes=genes.gene.tolist(),cohorts=['Wu2021','Pal2021_reprocessed'],counts_normalization='mean cellular log1p(count/library*10000) within source donor;equal donor weights across donors',min_cells_per_donor_type=20,min_source_labels_per_group=3,bootstrap=1000,ranking='all evaluable types per study;also compare shared types separately',partition=['ALL','Wu treatment-naive sensitivity'],no_new_hypothesis_tests=True,no_reclustering=True,annotation='frozen Wu author labels and Chen reannotation of Pal;not newly inferred',LPCAT4='Q643R3;stable ENSG00000176454 or ENSG00000291994;require unique hit;not MBOAT2 ENSG00000143797',source_commit='7d22872',raw_patient_export=False)
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));manifest=[];allr=[];allz=[];ident=[]
for cf in ['brca_sc117_Wu2021_config.json','brca_sc117_Pal2021_config.json']:
 cfg=json.loads((R/cf).read_text());src=PREV/'source'/cfg['file'];link=R/'source'/cfg['file'];link.symlink_to(src) if not link.exists() else None
 digest=sha(src);manifest.append(dict(path=str(src),sha256=digest,kind='source_h5ad'))
 with h5py.File(src,'r') as h:
  var=prof.dataframe(h[cfg['var_path']]);name=var[cfg['symbol_column']].astype(str);ids=var.index.astype(str).str.split('.').str[0];mask=ids.isin(['ENSG00000176454','ENSG00000291994','ENSG00000143797'])|name.isin(['LPCAT4','AGPAT7','LPEAT2','MBOAT2'])
  for gid,n in zip(var.index[mask],name[mask]):ident.append(dict(cohort=cfg['cohort'],requested_gene='LPCAT4',feature_id=gid,source_symbol=n,accepted_candidate_id=str(gid).split('.')[0] in ['ENSG00000176454','ENSG00000291994'],source='UniProt Q643R3 cached mapping190 source'))
 prof.run(R,R/cf)
 d=pd.read_csv(R/'private'/(cfg['cohort']+'_donor_profiles.tsv'),sep='\t');parts=[('ALL',d)]
 if cfg['cohort']=='Wu2021':parts.append(('Naive',d[d.treatment.eq('Naïve')]))
 for part,block in parts:
  seed=int(hashlib.sha256((V+cfg['cohort']+part).encode()).hexdigest()[:8],16);r,z=stab.one(block,cfg['cohort'],part,seed);allr.append(r);allz.append(z)
 assert digest==sha(src)
 manifest.append(dict(path=str(R/cf),sha256=sha(R/cf),kind='frozen_config'))
r=pd.concat(allr,ignore_index=True);z=pd.concat(allz,ignore_index=True);save(r,P/'new39_source_stability.tsv');save(z,P/'new39_source_distributions.tsv');save(pd.DataFrame(ident),P/'LPCAT4_stable_identity.tsv')
rows=[]
for g in genes.gene:
 w=r[(r.gene==g)&(r.cohort=='Wu2021')&(r.partition=='ALL')].iloc[0];p=r[(r.gene==g)&(r.cohort=='Pal2021_reprocessed')&(r.partition=='ALL')].iloc[0];b=dict(gene=g,Wu_status=w.status,Pal_status=p.status,Wu_top=w.top_lineage,Pal_top=p.top_lineage,Wu_top_frequency=w.bootstrap_top_frequency,Pal_top_frequency=p.bootstrap_top_frequency,Wu_paired_source_labels=w.paired_source_labels,Pal_paired_source_labels=p.paired_source_labels,Wu_paired_positive_fraction=w.paired_positive_fraction,Pal_paired_positive_fraction=p.paired_positive_fraction)
 ok=w.status=='DONE' and p.status=='DONE';b['same_top_all_categories']=bool(ok and w.top_lineage==p.top_lineage);b['same_top_and_both_bootstrap_ge080']=bool(b['same_top_all_categories'] and w.bootstrap_top_frequency>=.8 and p.bootstrap_top_frequency>=.8)
 zz=z[(z.gene==g)&(z.partition=='ALL')&z.status.eq('DONE')];a=zz[zz.cohort.eq('Wu2021')].set_index('celltype');c=zz[zz.cohort.eq('Pal2021_reprocessed')].set_index('celltype');shared=sorted(set(a.index)&set(c.index));b['n_shared_categories']=len(shared);b['Wu_top_shared']=a.loc[shared].effect.idxmax() if ok and shared else 'NA';b['Pal_top_shared']=c.loc[shared].effect.idxmax() if ok and shared else 'NA';b['same_top_shared_categories']=bool(ok and shared and b['Wu_top_shared']==b['Pal_top_shared']);b['interpretation']='descriptive expression location;not exclusive source or function';rows.append(b)
save(pd.DataFrame(rows),P/'new39_cross_study_comparison.tsv')
for p in [Path(__file__),R/'brca_sc117_profile_v1.py',R/'brca_all117_stability_v2.py',R/'profile39_adapter.py',R/'stability39_adapter.py',R/'source/genes_unique.tsv']:manifest.append(dict(path=str(p),sha256=sha(p),kind='code_or_gene_scope'))
save(pd.DataFrame(manifest),P/'source_manifest.tsv');co=pd.DataFrame(rows);val=dict(status='DONE',genes=39,same_top=int(co.same_top_all_categories.sum()),same_top_both_bootstrap080=int(co.same_top_and_both_bootstrap_ge080.sum()),same_top_shared=int(co.same_top_shared_categories.sum()),source_hashes_unchanged=True,new_P_or_q=0,cohorts={c:json.loads((P/(c+'_validation.json')).read_text()) for c in spec['cohorts']});(P/'validation.json').write_text(json.dumps(val,indent=2));(R/'NUMERICAL_DONE').write_text('DONE');print(json.dumps({k:v for k,v in val.items() if k!='cohorts'}),flush=True)
