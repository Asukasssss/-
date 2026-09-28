"""LYPLA1 existing malignant states and within-patient sparse partial correlation."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
R=B/'20260928T131206Z_lypla1_cellstates_v1';P=R/'public';D=R/'private'
MODES=['raw','depth','depth_state']
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.8g')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for z in iter(lambda:f.read(8*1024*1024),b''):h.update(z)
 return h.hexdigest()
def block(f,n,p,sel):
 r=f['raw/X'];ptr=r['indptr'][:]
 for st in range(0,n,2000):
  en=min(st+2000,n);m=sel[st:en]
  if not m.any():continue
  lo,hi=ptr[st],ptr[en];x=sparse.csr_matrix((r['data'][lo:hi].astype(float),r['indices'][lo:hi],ptr[st:en+1]-lo),shape=(en-st,p))[m];x.eliminate_zeros();yield np.flatnonzero(m)+st,x

def main():
 assert (R/'.running').exists()
 for d in [P,D]:d.mkdir(exist_ok=True)
 spec=dict(version='lypla1_cellstates_v1',source_commit='d9846c63e158adca5aad0b05c98652e858f17b4d',cohorts=['Wu2021','Pal2021_reprocessed'],unit='author/GEO patient label; Pal conflict excluded and explicit shared patient regions pooled',states='existing Wu celltype_minor and Pal author_cell_type; no cross-cohort state equivalence assumed; no reclustering',min_cells_patient=100,min_LYPLA1_positive=10,min_gene_positive_per_patient='max(10,1% cells)',cohort_gene_detection=.01,min_patients_gene=5,state_summary='>=20 cells per patient-state; >=3 patients for public aggregate; all rows retained',normalization='log1p(raw/library*10000)',correlation='within patient Pearson log expression; depth partial Pearson main; depth+existing state indicators sensitivity',covariates='log1p(total UMI minus target counts), log1p(detected genes minus target detection); per-patient standardization; intercept; sensitivity dummy indicators',aggregation='equal-patient median correlation; positive/negative fraction; two-sided binomial sign test ignoring exact zero; descriptive effects; BH all genes and 3 modes within cohort',replication_rule='depth median same direction with abs>=0.1 in both cohorts, >=70% patients same direction each, >=5 patients each; no significance gate; fixed before results',limitations='Exploratory associations; unmodelled technical/state factors; depth controls do not establish independence/causality; detection filter; no enrichment',software={'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__})
 (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));hp=B/'20260925T133935Z_lypla1_KEGGscore_v1/brca_sc117_profile_v1.py';sp=importlib.util.spec_from_file_location('h',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 manifests=[];states=[];coverage=[];qc=[];sums=[];audits=[]
 for co,fn in [('Wu2021','wu_curated.h5ad'),('Pal2021_reprocessed','pal_atlas_container.h5ad')]:
  print('START',co,flush=True);src=B/'20260919T140105Z_scRNA117_v1/source'/fn
  with h5py.File(src) as f:
   o=h.dataframe(f['obs']);v=h.dataframe(f['raw/var']);names=v.feature_name.astype(str).to_numpy();n,p=len(o),len(v);target=np.flatnonzero(names=='LYPLA1');assert len(target)==1;target=int(target[0]);unique=~pd.Series(names).duplicated(keep=False).to_numpy()
   if co=='Wu2021':
    sel=o.celltype_major.eq('Cancer Epithelial').to_numpy();patient=o.donor_id.astype(str).to_numpy();state=o.celltype_minor.astype(str).to_numpy();excluded=0
   else:
    mp=B/'20260927T122639Z_pal_subtypes_v1/private/sample_subtype_mapping.tsv';m=pd.read_csv(mp,sep='\t',dtype=str).set_index('donor_id');status=o.donor_id.map(m.status);sel=(o.batch.eq('pal_2021')&o.cell_type.eq('malignant cell')&status.eq('DONE')).to_numpy();patient=o.donor_id.map(m.patient).fillna('UNRESOLVED').to_numpy();state=o.author_cell_type.astype(str).to_numpy();excluded=int((o.batch.eq('pal_2021')&o.cell_type.eq('malignant cell')&~status.eq('DONE')).sum());manifests.append(dict(path=str(mp),sha256=sha(mp),hash_basis='current_file'))
   dep=np.full(n,np.nan);det=np.full(n,np.nan);ly=np.full(n,np.nan);linear=np.full(n,np.nan)
   for ii,x in block(f,n,p,sel):
    total=np.asarray(x.sum(1)).ravel();assert (total>0).all();xx=x[:,target].toarray().ravel();linear[ii]=xx/total*10000;ly[ii]=np.log1p(linear[ii]);dep[ii]=np.log1p(total-xx);det[ii]=np.log1p(np.asarray((x>0).sum(1)).ravel()-(xx>0))
   cells=pd.DataFrame({'patient':patient[sel],'state':state[sel],'log':ly[sel],'linear':linear[sel],'detected':ly[sel]>0});ps=cells.groupby(['patient','state']).agg(n_cells=('log','size'),mean_log=('log','mean'),mean_CP10k=('linear','mean'),detection=('detected','mean')).reset_index();save(ps,D/(co+'_patient_states.tsv'))
   for st,g in ps.groupby('state'):
    z=g[g.n_cells>=20];states.append(dict(cohort=co,state=st,n_cells=int(g.n_cells.sum()),n_patients_eligible=len(z),mean_log=z.mean_log.mean() if len(z)>=3 else np.nan,mean_CP10k=z.mean_CP10k.mean() if len(z)>=3 else np.nan,detection=z.detection.mean() if len(z)>=3 else np.nan,largest_patient_cell_fraction=g.n_cells.max()/g.n_cells.sum(),status='DONE' if len(z)>=3 else 'NOT_EVALUABLE',reason='>=3 patients each >=20 cells' if len(z)>=3 else 'insufficient_patient_state_coverage'))
   ids=sorted(set(patient[sel]));acc={};bad=[]
   for pa in ids:
    ii=np.flatnonzero(sel&(patient==pa));ok=len(ii)>=100 and (ly[ii]>0).sum()>=10 and np.std(ly[ii])>0;coverage.append(dict(cohort=co,n_cells=len(ii),target_positive=int((ly[ii]>0).sum()),status='DONE' if ok else 'NOT_EVALUABLE',reason='eligible' if ok else 'insufficient_cells_or_target_detection'))
    if not ok:continue
    q=np.c_[dep[ii],det[ii]];q=(q-q.mean(0))/np.maximum(q.std(0),1e-12);q=np.c_[np.ones(len(ii)),q];dummy=pd.get_dummies(pd.Series(state[ii]),drop_first=True).to_numpy(dtype=float);Q=np.c_[q,dummy];ranks=[1,np.linalg.matrix_rank(q),np.linalg.matrix_rank(Q)];ks=[1,3,Q.shape[1]]
    acc[pa]=dict(ii=ii,Q=Q,ks=ks,N=len(ii),sum=np.zeros(p),square=np.zeros(p),detect=np.zeros(p),xt=np.zeros(p),xq=np.zeros((p,Q.shape[1])),qq=Q.T@Q,qt=Q.T@ly[ii],t2=float(ly[ii]@ly[ii]),audit_y=[],audit_ii=[])
    for mode,k,rank in zip(MODES,ks,ranks):
     rt=ly[ii]-Q[:,:k]@np.linalg.lstsq(Q[:,:k],ly[ii],rcond=None)[0];qc.append(dict(cohort=co,mode=mode,patient_index=ids.index(pa)+1,n_cells=len(ii),design_rank=int(rank),target_residual_variance=float(np.var(rt)),target_depth_r=float(stats.pearsonr(ly[ii],dep[ii])[0])))
   auditix=[int(np.flatnonzero(names==g)[0]) for g in ['LYPLA2','GPCPD1','PCYT2','ETNK1','ACTB','GAPDH'] if (names==g).sum()==1]
   for ii,x in block(f,n,p,sel):
    total=np.asarray(x.sum(1)).ravel();y=x.multiply((10000/total)[:,None]).tocsr();y.data=np.log1p(y.data)
    for pa in np.unique(patient[ii]):
     if pa not in acc:continue
     a=acc[pa];mask=patient[ii]==pa;jj=ii[mask];z=y[mask];pos=np.searchsorted(a['ii'],jj);Q=a['Q'][pos]
     a['sum']+=np.asarray(z.sum(0)).ravel();a['square']+=np.asarray(z.multiply(z).sum(0)).ravel();a['detect']+=np.asarray((z>0).sum(0)).ravel();a['xt']+=np.asarray(z.T@ly[jj]).ravel();a['xq']+=z.T@Q;a['audit_y'].append(z[:,auditix].toarray());a['audit_ii'].extend(jj)
   cor=[];global_detect=sum(a['detect'] for a in acc.values());globaln=sum(a['N'] for a in acc.values());base=unique&(names!='LYPLA1')&(names!='')&(names!='nan')&(global_detect/globaln>=.01)
   save(pd.DataFrame({'gene':names,'eligible_cohort':base,'detection':global_detect/globaln}),P/(co+'_gene_coverage.tsv'))
   for pa,a in acc.items():
    N=a['N'];ok=base&(a['detect']>=max(10,.01*N));axy=np.vstack(a['audit_y']);aidx=np.array(a['audit_ii']);pos=np.searchsorted(a['ii'],aidx)
    for mode,k in zip(MODES,a['ks']):
     inv=np.linalg.pinv(a['qq'][:k,:k]);xq=a['xq'][:,:k];qt=a['qt'][:k];cross=a['xt']-xq@inv@qt;vx=a['square']-np.einsum('ij,ij->i',xq@inv,xq);vt=a['t2']-qt@inv@qt;valid=ok&(vx>1e-10)&(vt>1e-10);rr=np.full(p,np.nan);rr[valid]=np.clip(cross[valid]/np.sqrt(vx[valid]*vt),-1,1)
     cor.append(pd.DataFrame({'patient':pa,'gene':names[base],'mode':mode,'r':rr[base]}))
     Q=a['Q'][pos,:k];t=ly[aidx];rt=t-Q@np.linalg.lstsq(Q,t,rcond=None)[0]
     for j,ix in enumerate(auditix):
      if not valid[ix]:continue
      xx=axy[:,j];rx=xx-Q@np.linalg.lstsq(Q,xx,rcond=None)[0];check=np.corrcoef(rx,rt)[0,1];assert np.isclose(check,rr[ix],atol=1e-7);audits.append(dict(cohort=co,mode=mode,gene=names[ix],status='PASS'))
   c=pd.concat(cor,ignore_index=True);save(c,D/(co+'_all_patient_correlations.tsv'))
   for (mode,gene),g in c.groupby(['mode','gene'],sort=True):
    v=g.r.dropna().to_numpy();nn=len(v);pos=int((v>0).sum());neg=int((v<0).sum());valid=nn>=5;pv=stats.binomtest(pos,pos+neg,.5).pvalue if valid and pos+neg else np.nan;sums.append(dict(cohort=co,mode=mode,gene=gene,n_patients=nn,median_r=np.median(v) if nn else np.nan,mean_r=np.mean(v) if nn else np.nan,positive_fraction=pos/nn if nn else np.nan,negative_fraction=neg/nn if nn else np.nan,p_value=pv,status='DONE' if valid else 'NOT_EVALUABLE',reason='>=5 eligible patients' if valid else 'fewer_than5_patients'))
   print('DONE',co,'patients',len(acc),'excluded_cells',excluded,flush=True)
  # Reuse frozen source hash with explicit provenance, do not reread multi-GB source solely for hashing.
  old=B/('20260927T122639Z_pal_subtypes_v1/public/source_manifest.tsv' if co.startswith('Pal') else '20260925T151104Z_lypla1_highlow_v1/public/source_manifest.tsv');mm=pd.read_csv(old,sep='\t');hashval=mm.loc[mm.path.eq(str(src)),'sha256'].iloc[0];manifests.append(dict(path=str(src),sha256=hashval,hash_basis='inherited_frozen_manifest:'+str(old),current_size_bytes=src.stat().st_size))
 out=pd.DataFrame(sums);out['q_cohort_all3_modes']=np.nan
 for co,idx in out.groupby('cohort').groups.items():
  valid=out.loc[idx,'p_value'].notna();ix=np.asarray(list(idx))[valid];out.loc[ix,'q_cohort_all3_modes']=multipletests(out.loc[ix,'p_value'],method='fdr_bh')[1]
 for (co,mode),g in out.groupby(['cohort','mode']):save(g,P/(co+'_'+mode+'_gene_summary.tsv'))
 a=out[(out.cohort=='Wu2021')&(out['mode']=='depth')];b=out[(out.cohort=='Pal2021_reprocessed')&(out['mode']=='depth')];cols=['gene','n_patients','median_r','positive_fraction','negative_fraction','p_value','q_cohort_all3_modes','status'];cross=a[cols].merge(b[cols],on='gene',how='outer',suffixes=('_Wu','_Pal'));cross['same_direction']=np.sign(cross.median_r_Wu)==np.sign(cross.median_r_Pal);cross['minimum_abs_median_r']=np.minimum(abs(cross.median_r_Wu),abs(cross.median_r_Pal));cross['replicated_descriptive']=cross.same_direction&(cross.minimum_abs_median_r>=.1)&(cross.n_patients_Wu>=5)&(cross.n_patients_Pal>=5)&(((cross.positive_fraction_Wu>=.7)&(cross.positive_fraction_Pal>=.7))|((cross.negative_fraction_Wu>=.7)&(cross.negative_fraction_Pal>=.7)))
 for co,short in [('Wu2021','Wu'),('Pal2021_reprocessed','Pal')]:
  extra=out[(out.cohort==co)&(out['mode']=='depth_state')][['gene','median_r','n_patients']].rename(columns={'median_r':'state_adjusted_median_r_'+short,'n_patients':'state_adjusted_n_'+short});cross=cross.merge(extra,on='gene',how='left')
 cross=cross.sort_values(['replicated_descriptive','minimum_abs_median_r'],ascending=False);save(cross,P/'cross_cohort_companions.tsv');save(cross[cross.replicated_descriptive],P/'replicated_descriptive_candidates.tsv');save(pd.DataFrame(states),P/'state_summary.tsv');save(pd.DataFrame(coverage),P/'patient_coverage_anonymous.tsv');save(pd.DataFrame(qc),D/'depth_qc_patient.tsv');save(pd.DataFrame(audits),P/'independent_residual_audit.tsv')
 q=pd.DataFrame(qc).groupby(['cohort','mode']).agg(n_patients=('n_cells','size'),median_target_depth_r=('target_depth_r','median'),median_target_residual_variance=('target_residual_variance','median')).reset_index();save(q,P/'depth_qc_summary.tsv')
 manifests.append(dict(path=str(Path(__file__)),sha256=sha(__file__),hash_basis='current_script'));save(pd.DataFrame(manifests),P/'source_manifest.tsv');(P/'validation.json').write_text(json.dumps(dict(status='PASS',independent_dense_residual_checks=len(audits),no_source_matrix_modified=True,old_results_unchanged=True,no_reclustering=True,source_hashes='inherited verified historical manifests, not rehashed this run',uncertainty='No patient bootstrap CI; P is direction consistency sign test, not cell-level association P'),indent=2));(R/'NUMERICAL_DONE').write_text('DONE');print('ALL_DONE replicated',int(cross.replicated_descriptive.sum()),flush=True)
if __name__=='__main__':main()
