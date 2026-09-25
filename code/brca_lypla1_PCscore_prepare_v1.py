"""Server preparation only: freeze coverage-based membership before correlations."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,importlib.util,platform
import numpy as np,pandas as pd,h5py
from scipy import sparse

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
BASE=ROOT/'results/collaborative/BRCA/A'
VERSION='lypla1_PCscore_v1'


def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()


def main():
 R=Path(sys.argv[1]);P=R/'public';D=R/'private';S=R/'source'
 assert (R/'.running').read_text().strip()==VERSION and not (P/'frozen_score_members.tsv').exists()
 spec=json.loads((S/'analysis_spec.json').read_text());(P/'analysis_spec.json').write_bytes((S/'analysis_spec.json').read_bytes())
 members=pd.read_csv(S/'reactome_members.tsv',sep='\t');assert len(members)==27 and members.gene.is_unique
 genes=sorted(set(members.gene)|{'LYPLA1','LYPLA2'});eligible={};coverage=[];summary=[];manifest=[]
 helper=importlib.util.spec_from_file_location('sc_helper',R/'brca_sc117_profile_v1.py');h=importlib.util.module_from_spec(helper);helper.loader.exec_module(h)
 for co,cfgname in [('Wu2021','brca_sc117_Wu2021_config.json'),('Pal2021_reprocessed','brca_sc117_Pal2021_config.json')]:
  cfg=json.loads((S/cfgname).read_text());src=BASE/'20260919T140105Z_scRNA117_v1/source'/cfg['file']
  print('Preparing',co,flush=True)
  manifest.append(dict(path=str(src),sha256=sha(src)))
  with h5py.File(src,'r') as f:
   obs=h.dataframe(f['obs']);var=h.dataframe(f[cfg['var_path']]);names=var[cfg['symbol_column']].astype(str)
   selected=obs[cfg['celltype_column']].map(cfg['celltype_map']).eq('Malignant_epithelial')
   if cfg.get('filter_column'):selected &= obs[cfg['filter_column']].eq(cfg['filter_value'])
   for key,labels in cfg.get('exclude_labels',{}).items():selected &= ~obs[key].isin(labels)
   donor=obs[cfg['donor_column']].astype(str)
   assert not (selected & donor.isin(['NA','nan','unknown',''])).any()
   nper=donor[selected].value_counts();keep=sorted(nper[nper>=spec['min_malignant_cells_per_donor']].index)
   selected &= donor.isin(keep);codes=pd.Categorical(donor,categories=keep).codes
   idx={g:int(np.flatnonzero(names.eq(g))[0]) for g in genes if names.eq(g).sum()==1}
   present=[g for g in genes if g in idx]
   sums=np.zeros((len(keep),len(present)));detect=np.zeros_like(sums);libsum=np.zeros(len(keep));nc=np.zeros(len(keep),int)
   raw=f[cfg['raw_path']];ptr=raw['indptr'][:];shape=tuple(raw.attrs['shape'])
   assert raw.attrs['encoding-type']=='csr_matrix' and shape==(len(obs),len(var))
   for start in range(0,len(obs),2000):
    end=min(start+2000,len(obs));mask=selected.iloc[start:end].to_numpy()
    if not mask.any():continue
    lo,hi=ptr[start],ptr[end];x=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(var)))[mask]
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data),atol=1e-6)
    lib=np.asarray(x.sum(1)).ravel();assert (lib>0).all()
    y=x[:,[idx[g] for g in present]].toarray();c=codes[start:end][mask]
    np.add.at(sums,c,y);np.add.at(detect,c,y>0);np.add.at(libsum,c,lib);np.add.at(nc,c,1)
   assert nc.sum()==selected.sum() and (nc>=20).all()
   expr=pd.DataFrame(np.log1p(sums/libsum[:,None]*10000),index=keep,columns=present)
   expr.index.name='donor';expr.to_csv(D/(co+'_expression.tsv'),sep='\t')
   meta=[]
   for i,g in enumerate(keep):
    sub=obs.loc[selected&donor.eq(g)]
    label=lambda c:';'.join(sorted(sub[c].astype(str).unique())) if c else 'NOT_AVAILABLE'
    meta.append(dict(donor=g,n_cells=int(nc[i]),library_sum=float(libsum[i]),subtype=label(cfg.get('stratum_column')),treatment=label(cfg.get('treatment_column'))))
   save(pd.DataFrame(meta),D/(co+'_metadata.tsv'))
   eligible[co]=[]
   for g in genes:
    exact=g in idx;k=present.index(g) if exact else None
    nd=int((sums[:,k]>0).sum()) if exact else 0
    sd=float(expr[g].std(ddof=1)) if exact else np.nan
    ok=exact and nd>=3 and np.isfinite(sd) and sd>0
    if ok:eligible[co].append(g)
    coverage.append(dict(cohort=co,gene=g,n_symbol_matches=int(names.eq(g).sum()),n_donors=len(keep),n_donors_detected=nd,
                         mean_donor_detection_fraction=float(np.mean(detect[:,k]/nc)) if exact else np.nan,
                         sd_log1p10k=sd,eligible=ok,reason='OK' if ok else 'absent_ambiguous_or_detected_in_less_than3_donors_or_constant'))
   summary.append(dict(cohort=co,unit='publisher_donor_label',source_malignant_donors=len(nper),eligible_donors=len(keep),malignant_cells=int(nc.sum()),min_cells=int(nc.min()),max_cells=int(nc.max()),annotation_origin=cfg['annotation_origin']))

 identity=BASE/'20260921T102429Z_camp_sample_identity_v1/private/audited_mapping.tsv';m=pd.read_csv(identity,sep='\t',dtype=str)
 m=m[m.TN.eq('Tumor')&m.TN.eq(m.pdf_TN)].sort_values('case_row').reset_index(drop=True)
 assert len(m)==60 and m.case_row.is_unique
 rf=ROOT/'data/candidates/camp_primary_tissue_multicancer/gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
 rna=pd.read_csv(rf,index_col=0);rna.columns=rna.columns.astype(str);assert rna.index.is_unique
 save(m,D/'CAMP_audited_tumors.tsv')
 expr=rna.reindex(genes)[m.RNAID].T;expr.index=m.case_row;expr.index.name='case_row';expr.to_csv(D/'CAMP_expression.tsv',sep='\t')
 eligible['CAMP']=[]
 for g in genes:
  sd=float(expr[g].std(ddof=1));ok=g in rna.index and np.isfinite(expr[g]).all() and sd>0 and g!='LPCAT4'
  if ok:eligible['CAMP'].append(g)
  coverage.append(dict(cohort='CAMP',gene=g,n_symbol_matches=int(g in rna.index),n_donors=60,n_donors_detected=np.nan,
                       mean_donor_detection_fraction=np.nan,sd_log1p10k=np.nan,sd_author_scale=sd,eligible=ok,
                       reason='historical_LPCAT4_symbol_identity_unresolved' if g=='LPCAT4' else ('OK' if ok else 'not_unique_finite_nonconstant_exact_symbol')))
 retained=sorted(set(members.gene).intersection(*[set(v) for v in eligible.values()])-{'LYPLA1','LYPLA2'})
 assert len(retained)>=5 and len(retained)/len(members)>=.5,'Coverage insufficient; stop without association tests'
 freeze=members.copy();freeze['used_in_all3_cohorts']=freeze.gene.isin(retained)
 freeze['excluded_target_or_comparator']=freeze.gene.isin(['LYPLA1','LYPLA2'])
 save(freeze,P/'frozen_score_members.tsv');save(pd.DataFrame(coverage),P/'gene_coverage.tsv');save(pd.DataFrame(summary),P/'sc_cohort_coverage.tsv')
 # No outcome correlations have been run before this coverage-only freeze.
 (P/'membership_freeze.json').write_text(json.dumps(dict(status='FROZEN_BEFORE_ASSOCIATION',genes=retained,n_members=len(members),n_used=len(retained),
  membership_sha256=sha(P/'frozen_score_members.tsv'),rule='exact common coverage,nonconstant;>=3 detected source donors;targets excluded;CAMP LPCAT4 blocked',
  scoring='within each cohort, z=(gene donor/case expression-mean)/sampleSD; unweighted mean of the same retained genes; no cross-cohort score magnitude comparisons'),indent=2))
 for p in [identity,rf,Path(__file__),R/'brca_sc117_profile_v1.py']+list(S.iterdir()):
  if p.is_file():manifest.append(dict(path=str(p),sha256=sha(p)))
 save(pd.DataFrame(manifest),P/'preparation_manifest.tsv')
 (R/'PREPARATION_DONE').write_text('DONE')
 print('FROZEN',len(retained),retained,flush=True);print(pd.DataFrame(summary).to_string(index=False),flush=True)

if __name__=='__main__':main()
