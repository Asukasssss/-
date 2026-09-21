"""Incremental FUSCC/Tang associations; patient-level data never leave server165."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,re,concurrent.futures
import numpy as np,pandas as pd,requests
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');BASE=ROOT/'results/collaborative/BRCA/A'
R=Path(sys.argv[1]);S=R/'source';P=R/'public';V='mapping262_external_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def norm(x):return re.sub(r'[\s-]+','',str(x).lower().strip().rstrip('*'))
def infer(x,y,seed):
 rng=np.random.default_rng(seed);n=len(x);a=stats.rankdata(x);b=stats.rankdata(y);a-=a.mean();b-=b.mean();den=np.linalg.norm(a)*np.linalg.norm(b);rho=float(a@b/den)
 pp=np.array([rng.permutation(b) for _ in range(9999)])@a/den;p=(1+(np.abs(pp)>=abs(rho)-1e-12).sum())/10000
 ix=rng.integers(0,n,(2000,n));ar=stats.rankdata(x[ix],axis=1);br=stats.rankdata(y[ix],axis=1);ar-=ar.mean(1,keepdims=True);br-=br.mean(1,keepdims=True);dd=np.sqrt((ar*ar).sum(1)*(br*br).sum(1));bs=np.divide((ar*br).sum(1),dd,out=np.full(2000,np.nan),where=dd>0);lo,hi=np.nanquantile(bs,[.025,.975]);return rho,p,lo,hi
def api_more(api,profile,samples,genes,prefix):
 info={};coverage=[]
 def one(g):
  p=S/(prefix+'_gene_'+g+'.json')
  try:
   if not p.exists():
    r=requests.get(api+'/genes/'+g,timeout=30);r.raise_for_status();p.write_text(json.dumps(r.json()))
   j=json.loads(p.read_text());assert j['hugoGeneSymbol']==g
   return g,j,''
  except Exception as e:return g,None,type(e).__name__+':'+str(e)[:100]
 for g,j,err in concurrent.futures.ThreadPoolExecutor(max_workers=3).map(one,genes):
  coverage.append(dict(cohort=prefix,gene=g,status='DONE' if j else 'ACCESS_BLOCKED',reason=err,entrez_id=j['entrezGeneId'] if j else np.nan))
  if j:info[g]=j
 rows=[]
 for i in range(0,len(info),12):
  gs=sorted(info)[i:i+12];p=S/(prefix+'_RNA_'+str(i)+'.json')
  try:
   if not p.exists():
    r=requests.post(api+'/molecular-profiles/'+profile+'/molecular-data/fetch',params={'projection':'DETAILED'},json={'entrezGeneIds':[info[g]['entrezGeneId'] for g in gs],'sampleIds':[x['sampleId'] for x in samples]},timeout=90);r.raise_for_status();p.write_text(json.dumps(r.json()))
   data=json.loads(p.read_text());byid={j['entrezGeneId']:g for g,j in info.items()};samplemap={x['sampleId']:x['patientId'] for x in samples}
   for x in data:
    assert x['sampleId'] in samplemap and samplemap[x['sampleId']]==x['patientId'];eid=x['gene']['entrezGeneId'];assert eid in byid
    rows.append(dict(patient=x['patientId'],gene=byid[eid],value=x['value']))
   print(prefix,'RNA',i,len(data),flush=True)
  except Exception as e:
   for c in coverage:
    if c['gene'] in gs:c.update(status='ACCESS_BLOCKED',reason=type(e).__name__+':RNA_batch_failed')
 save(pd.DataFrame(coverage),P/(prefix+'_new_gene_api_coverage.tsv'))
 if not rows:return pd.DataFrame()
 d=pd.DataFrame(rows);assert not d.duplicated(['patient','gene']).any();return d.pivot(index='patient',columns='gene',values='value')
def main():
 assert (R/'.running').read_text().strip()==V
 rel=pd.read_csv(S/'direct_relations.tsv',sep='\t');assert len(rel)==262 and rel.relation_key.is_unique
 spec=dict(version=V,source_commit='9e4387f',planned_relations_per_family=262,FUSCC_families=['FUSCC_PRIMARY262','FUSCC_EXCLUDE_WARNING262'],Tang_family='Tang262',FUSCC_BH='262 planned;missing p=1 internally;public p/q stay NA',Tang_BH='all evaluable of262',method='Spearman9999 permutation plus1;2000case percentile bootstrap95CI',min_n={'FUSCC':20,'Tang':10},reuse='same source matrices+join;exact relation;n/effect checked;reuse old P/CI then new BH',LPCAT4='pending stable gene identity;no symbol-only substitution',Oslo='not reattempted;true crosswalk still required',patient_values_exported=False)
 (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 fprev=BASE/'20260920T090922Z_all174_external117_resources_v1';tprev=BASE/'20260921T072503Z_oslo_tang_external_v1';tp=tprev/'source'
 annp=ROOT/'results/CAMP_first_analysis_20260910/tables/FUSCC_author_polar_annotation_594.tsv';ann=pd.read_csv(annp,sep='\t');assert ann.fuscc_peak.is_unique
 cv=[];td=pd.read_excel(tp/'Tang_S2.xlsx',header=None);names=td.iloc[19:,1].astype(str).str.strip();ids=td.iloc[19:,0].astype(int).values
 for key,name in rel[['metabolite_key','metabolite_name']].drop_duplicates().values:
  typ,val=key.split(':',1);a=ann[ann.fuscc_kegg.eq(val)] if typ=='KEGG' else (ann[ann.fuscc_hmdb.eq('HMDB'+val[4:].zfill(7))] if typ=='HMDB' else ann.iloc[:0])
  state='DONE' if len(a)==1 else 'NOT_EVALUABLE';reason='unique_author_identifier' if len(a)==1 else ('multiple_features_unresolved' if len(a)>1 else 'no_exact_author_identifier')
  if key=='KEGG:C00186' and len(a)==1:state='NEEDS_REVIEW';reason='L_vs_DL_lactate_scope_unresolved'
  cv.append(dict(cohort='FUSCC',metabolite_key=key,metabolite_name=name,status=state,reason=reason,feature=str(a.iloc[0].fuscc_peak) if len(a)==1 else '',external_name=a.iloc[0].fuscc_author_name if len(a)==1 else '',row_index=-1))
  target={'KEGG:C00668':'glucose-6-phosphate (G6P)'}.get(key,name);hits=[i for i,x in enumerate(names) if norm(x)==norm(target)];state='DONE' if len(hits)==1 else 'NOT_EVALUABLE';reason='unique_author_name;not_spectral_reidentification' if len(hits)==1 else 'no_unique_author_name'
  if key=='KEGG:C00186' and len(hits)==1:state='NEEDS_REVIEW';reason='L_lactate_scope_unresolved'
  cv.append(dict(cohort='Tang',metabolite_key=key,metabolite_name=name,status=state,reason=reason,feature=str(ids[hits[0]]) if len(hits)==1 else '',external_name=names.iloc[hits[0]] if len(hits)==1 else '',row_index=19+hits[0] if len(hits)==1 else -1))
 cv=pd.DataFrame(cv);save(cv,P/'identity_coverage.tsv')
 fwpath=fprev/'private/author_patient_joined_values.tsv';fz=pd.read_csv(fwpath,sep='\t').set_index('patient_id');assert fz.index.is_unique and len(fz)==258
 # Keep old numeric columns in memory unchanged; no extra TSV roundtrip for old values.
 fsamples=json.loads((ROOT/'data/candidates/BRCA_1MNA_followup_v0.1/Fudan_samples.json').read_text())
 fc=cv[cv.cohort.eq('FUSCC')].set_index('metabolite_key');needed=sorted({t.gene for t in rel.itertuples() if fc.loc[t.metabolite_key,'status']=='DONE' and t.gene not in fz and t.gene!='LPCAT4'})
 fx=api_more('https://data.3steps.cn/cdataportal/api','FUSCC_BRCA_2022_mrna_seq_tpm',fsamples,needed,'FUSCC')
 for g in fx:fz[g]=fx[g].reindex(fz.index)
 wp=ROOT/'data/candidates/brca_four_axes_20260908/FUSCC_41422_2022_614_MOESM12_ESM.xlsx';mat=pd.read_excel(wp,sheet_name='Table S1',header=1);mat.columns=[str(c).strip() for c in mat];mat=mat.set_index('Peak');assert mat.index.is_unique
 for peak in fc.loc[fc.status.eq('DONE'),'feature'].unique():
  if peak not in fz:
   values=pd.to_numeric(mat.loc[peak,fz.metabolite_sample_id],errors='coerce').to_numpy();tmp=R/'private'/('new_peak_'+peak+'.tsv');save(pd.DataFrame({'value':values}),tmp);fz[peak]=pd.read_csv(tmp,sep='\t').value.to_numpy()
 oldf=pd.read_csv(fprev/'public/external_all174_associations.tsv',sep='\t').set_index(['relation_id','test_family'])
 tj=json.loads((tp/'Tang_tcga_sample_join.json').read_text());data=[x for p in tp.glob('tcga_rna_*.json') for x in json.loads(p.read_text())];ex=pd.DataFrame([dict(patient=x['patientId'],gene=x['gene']['hugoGeneSymbol'],value=x['value']) for x in data]);assert not ex.duplicated(['patient','gene']).any();tz=ex.pivot(index='patient',columns='gene',values='value');assert len(tz)==20
 tx=api_more('https://www.cbioportal.org/api','brca_tcga_pub2015_rna_seq_v2_mrna',tj,sorted(set(rel.gene)-set(tz.columns)-{'LPCAT4'}),'Tang')
 for g in tx:tz[g]=tx[g].reindex(tz.index)
 tc=cv[cv.cohort.eq('Tang')].set_index('metabolite_key');cols={('TCGA-'+str(td.iloc[0,k]).strip()):k for k in range(3,33) if re.fullmatch(r'[A-Z0-9]{2}-[A-Z0-9]{4}',str(td.iloc[0,k]).strip())};assert len(cols)==25
 for k,c in tc.iterrows():
  if c.status=='DONE':tz['met_'+k]=pd.Series({p:pd.to_numeric(td.iloc[int(c.row_index),col],errors='coerce') for p,col in cols.items()}).reindex(tz.index)
 oldt=pd.read_csv(tprev/'public/Tang_all174_associations.tsv',sep='\t').set_index('relation_id');save(fz.reset_index(),R/'private/FUSCC_joined.tsv');save(tz.reset_index(),R/'private/Tang_joined.tsv')
 results=[]
 for fam,cohort,z,cover,oldfam in [('FUSCC_PRIMARY262','FUSCC',fz,fc,'FUSCC_PRIMARY174'),('FUSCC_EXCLUDE_WARNING262','FUSCC',fz[~fz.index.isin(['FUSCCTNBC030','FUSCCTNBC044','FUSCCTNBC140'])],fc,'FUSCC_EXCLUDE_WARNING174'),('Tang262','Tang',tz,tc,'')]:
  for t in rel.itertuples():
   c=cover.loc[t.metabolite_key];b=dict.fromkeys(PREFIX,np.nan);b.update(cancer='BRCA',cohort='FUSCC_TNBC' if cohort=='FUSCC' else 'Tang2014',stage_id='06_EXTERNAL',run_id=R.name,analysis_version=V,analysis_type='author_patient_linked_spearman',metabolite_key=t.metabolite_key,metabolite_name=t.metabolite_name,gene=t.gene,unit='author_patient_identifier',effect_type='Spearman_rho',test_family=fam,status=c.status,reason=c.reason,source_id='10.1038/s41422-022-00614-0' if cohort=='FUSCC' else '10.1186/s13058-014-0415-9',relation_id=t.relation_key,external_feature=c.feature,statistics_reused=False)
   pr=(oldf.loc[(t.relation_key,oldfam)] if (t.relation_key,oldfam) in oldf.index else None) if cohort=='FUSCC' else (oldt.loc[t.relation_key] if t.relation_key in oldt.index else None)
   if pr is not None:b.update(previous_p=pr.p_value,previous_q=pr.q_value,previous_effect=pr.effect)
   if t.gene=='LPCAT4':b.update(status='NEEDS_REVIEW',reason='historical_gene_symbol_ambiguity_requires_stable_ID')
   elif c.status=='DONE':
    peak=c.feature if cohort=='FUSCC' else 'met_'+t.metabolite_key
    if t.gene not in z:b.update(status='NOT_EVALUABLE',reason='RNA_unavailable;see_gene_API_coverage')
    else:
     a=z[[t.gene,peak]].replace([np.inf,-np.inf],np.nan).dropna();b['n']=len(a)
     if len(a)<(20 if cohort=='FUSCC' else 10) or a.nunique().min()<2:b.update(status='NOT_EVALUABLE',reason='insufficient_or_constant')
     else:
      rho=stats.spearmanr(a.iloc[:,0],a.iloc[:,1]).statistic
      if pr is not None and pd.notna(pr.p_value):
       assert len(a)==pr.n and abs(rho-pr.effect)<1e-12,(fam,t.relation_key,len(a),pr.n,rho,pr.effect)
       b.update({k:pr[k] for k in ['effect','p_value','ci_lower','ci_upper']});b['statistics_reused']=True
      else:
       sd=int(hashlib.sha256((V+fam+t.relation_key).encode()).hexdigest()[:8],16);rho,p,lo,hi=infer(a.iloc[:,0].to_numpy(),a.iloc[:,1].to_numpy(),sd);b.update(effect=rho,p_value=p,ci_lower=lo,ci_upper=hi,seed=sd)
      b.update(status='DONE',reason='unadjusted;author_processed_metabolites;same_aliquot_unverified')
   results.append(b)
  print('complete',fam,flush=True)
 a=pd.DataFrame(results)
 for f,ix in a.groupby('test_family').groups.items():
  ok=a.loc[ix,'p_value'].notna();ii=np.array(list(ix))[ok];a.loc[ix,'family_n_evaluable']=len(ii)
  if f.startswith('FUSCC'):
   qs=multipletests(a.loc[ix,'p_value'].fillna(1),method='fdr_bh')[1];a.loc[ii,'q_value']=qs[ok]
  elif len(ii):a.loc[ii,'q_value']=multipletests(a.loc[ii,'p_value'],method='fdr_bh')[1]
 a=a[PREFIX+[c for c in a if c not in PREFIX]];save(a,P/'external_relations262.tsv')
 paths=[annp,wp,fwpath,tp/'Tang_S2.xlsx',tp/'Tang_tcga_sample_join.json',fprev/'public/external_all174_associations.tsv',tprev/'public/Tang_all174_associations.tsv',S/'direct_relations.tsv',Path(__file__)]+list(tp.glob('tcga_rna_*.json'))+list(S.glob('*.json'));save(pd.DataFrame([dict(path=str(p),sha256=sha(p)) for p in paths]),P/'source_manifest.tsv')
 v=dict(status='DONE',old_overlap_n_and_effect_verified=True,old_statistics_unchanged=True,patient_values_exported=False,families={f:dict(planned=len(d),evaluable=int(d.p_value.notna().sum()),reused=int(d.statistics_reused.sum()),newly_calculated=int((d.p_value.notna()&~d.statistics_reused).sum()),q_lt005=int(d.q_value.lt(.05).sum())) for f,d in a.groupby('test_family')});(P/'validation.json').write_text(json.dumps(v,indent=2));(R/'NUMERICAL_DONE').write_text('DONE');print(json.dumps(v),flush=True)
if __name__=='__main__':main()
