"""Fixed original174 associations; processed-source identifiers; no patient values exported."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
from pathlib import Path
import sys,json,re,hashlib,io
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
R=Path(sys.argv[1]);S=R/'source';P=R/'public'
def save(a,p):a.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def norm(x):return re.sub(r'[\s-]+','',str(x).lower().strip().rstrip('*'))
def row(cohort,g,met='NA',key='NA',kind='tumor_spearman_permutation'):
 return dict(cancer='BRCA',cohort=cohort,stage_id='06_EXTERNAL',run_id=R.name,analysis_version='oslo_tang_external_v1',analysis_type=kind,metabolite_key=key,metabolite_name=met,gene=g,unit='author_case_identifier',n=np.nan,n_reference=np.nan,effect_type='Spearman_rho',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family=cohort+'_original174',family_n_evaluable=0,status='NOT_EVALUABLE',reason='not_measured',source_id=cohort)
def infer(x,y,seed):
 rng=np.random.default_rng(seed);n=len(x);a=stats.rankdata(x);b=stats.rankdata(y);a-=a.mean();b-=b.mean();rho=float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))
 pp=np.array([rng.permutation(b) for _ in range(9999)])@a/(np.linalg.norm(a)*np.linalg.norm(b));p=(1+(np.abs(pp)>=abs(rho)-1e-12).sum())/10000
 ix=rng.integers(0,n,(2000,n));ar=stats.rankdata(x[ix],axis=1);br=stats.rankdata(y[ix],axis=1);ar-=ar.mean(1,keepdims=True);br-=br.mean(1,keepdims=True);den=np.sqrt((ar*ar).sum(1)*(br*br).sum(1));bs=np.divide((ar*br).sum(1),den,out=np.full(2000,np.nan),where=den>0);ci=np.nanquantile(bs,[.025,.975]);assert np.isclose(rho,stats.spearmanr(x,y).statistic);return rho,p,*ci,int(np.isfinite(bs).sum())
def family(rows):
 a=pd.DataFrame(rows);ix=a.status.eq('DONE');a.loc[ix,'q_value']=multipletests(a.loc[ix,'p_value'],method='fdr_bh')[1] if ix.any() else [];a['family_n_evaluable']=int(ix.sum());return a
def main():
 rel=pd.read_csv(S/'relations174.tsv',sep='\t');assert len(rel)==174 and rel.relation_id.is_unique and rel.gene.nunique()==117
 d=pd.read_excel(S/'Tang_S2.xlsx',header=None);assert d.shape==(418,33)
 names=d.iloc[19:,1].astype(str).str.strip();ids=d.iloc[19:,0].astype(int).values;assert len(names)==399 and len(set(ids))==399
 aliases={'KEGG:C00668':'glucose-6-phosphate (G6P)'}
 osm={'KEGG:C00025':'Glutamate','KEGG:C00037':'Glycine','KEGG:C00064':'Glutamine','KEGG:C00082':'L-Tyrosine','KEGG:C00114':'Choline','KEGG:C00245':'Taurine','KEGG:C00300':'Creatine','KEGG:C00670':'Glycerophosphocholine','KEGG:C00186':'Lactate','KEGG:C00051':'Glutathione'}
 od=pd.read_excel(S/'Oslo2_tables.xlsx',sheet_name='Add. Table 2',header=1);assert len(od)==228 and od['Sample number'].is_unique
 cov=[]
 for key,name in rel[['metabolite_key','metabolite']].drop_duplicates().values:
  target=aliases.get(key,name);hits=[i for i,x in enumerate(names) if norm(x)==norm(target)]
  state='DONE' if len(hits)==1 else 'NOT_EVALUABLE';reason='unique_author_name_match;not_spectral_reidentification' if len(hits)==1 else 'no_unique_author_name_match'
  if key=='KEGG:C00186' and len(hits)==1:state='NEEDS_REVIEW';reason='CAMP_L_lactate_key;author_unspecified_stereochemistry'
  cov.append(dict(cohort='Tang2014',metabolite_key=key,CAMP_name=name,external_feature_id=str(ids[hits[0]]) if len(hits)==1 else 'NA',external_name=names.iloc[hits[0]] if len(hits)==1 else 'NA',source_row_index=19+hits[0] if len(hits)==1 else -1,status=state,reason=reason))
  found=key in osm;st='DONE' if found else 'NOT_EVALUABLE';rs='explicit_author_name_correspondence;not_spectral_reidentification' if found else 'not_in_18_metabolite_panel'
  if key in ['KEGG:C00186','KEGG:C00051'] and found:st='NEEDS_REVIEW';rs='stereochemical_or_redox_scope_not_explicit'
  cov.append(dict(cohort='Oslo2',metabolite_key=key,CAMP_name=name,external_feature_id=osm.get(key,'NA'),external_name=osm.get(key,'NA'),source_row_index=-1,status=st,reason=rs))
 cv=pd.DataFrame(cov);save(cv,P/'identity_coverage.tsv')
 # Frozen mapping contains only names and identifiers, no statistical outcomes.
 sp=dict(source_commit='66bc398',relations_planned_per_cohort=174,genes=117,minimum_complete_cases=10,method='Spearman two-sided permutation9999 plus1;pointwise percentile bootstrap2000;no new imputation',family='BH separately per cohort over all evaluable original174;original q unchanged',seed_base=20260921,tcga_profile='brca_tcga_pub2015_rna_seq_v2_mrna',tcga_join='exact author TCGA participant barcode;one primary sample per patient;different aliquot cannot be ruled out',Tang_processing='author normalized and imputed values reused;RNA RSEM from cBioPortal2015 not assumed exact2014 release',Oslo_join='MCOS versus OSL2 identifiers;no order-based inference',RNA117_method='GSE42568 author log2 GCRMA;median all uniquely annotated probes per gene;unpaired Welch mean difference with Welch95CI;BH all evaluable117;no further normalization')
 spec=P/'analysis_spec.json'
 if spec.exists():assert json.loads(spec.read_text())==sp
 else:spec.write_text(json.dumps(sp,indent=2)+'\n')
 j=json.loads((S/'Tang_tcga_sample_join.json').read_text());assert len({x['patientId'] for x in j})==len(j)
 data=[x for p in S.glob('tcga_rna_*.json') for x in json.loads(p.read_text())];ex=pd.DataFrame([dict(patient=x['patientId'],sample=x['sampleId'],gene=x['gene']['hugoGeneSymbol'],value=x['value']) for x in data]);assert not ex.duplicated(['patient','gene']).any();assert set(ex['sample'])<=set(x['sampleId'] for x in j)
 wide=ex.pivot(index='patient',columns='gene',values='value');cols={('TCGA-'+str(d.iloc[0,k]).strip()):k for k in range(3,33) if re.fullmatch(r'[A-Z0-9]{2}-[A-Z0-9]{4}',str(d.iloc[0,k]).strip())};assert len(cols)==25
 result=[];oslo=[]
 for i,r in rel.iterrows():
  c=cv[(cv.cohort=='Tang2014')&(cv.metabolite_key==r.metabolite_key)].iloc[0];rr=row('Tang2014',r.gene,r.metabolite,r.metabolite_key);rr.update(relation_id=r.relation_id,external_feature_id=c.external_feature_id,status=c.status,reason=c.reason)
  if c.status=='DONE':
   if r.gene not in wide:rr.update(status='NOT_EVALUABLE',reason='RNA_gene_not_available')
   else:
    values=pd.Series({p:pd.to_numeric(d.iloc[int(c.source_row_index),col],errors='coerce') for p,col in cols.items()});z=pd.concat([values.rename('met'),wide[r.gene].rename('RNA')],axis=1).replace([np.inf,-np.inf],np.nan).dropna();rr['n']=len(z)
    if len(z)<10 or z.nunique().min()<2:rr.update(status='NOT_EVALUABLE',reason='less_than10_complete_cases_or_constant')
    else:
     rho,p,lo,hi,nb=infer(z.met.values,z.RNA.values,20260921+i);rr.update(status='DONE',reason='author_TCGA_case_join;imputed_metabolomics;different_aliquot_possible;unadjusted',effect=rho,p_value=p,ci_lower=lo,ci_upper=hi,n_bootstrap_valid=nb,seed=20260921+i)
  result.append(rr)
  oc=cv[(cv.cohort=='Oslo2')&(cv.metabolite_key==r.metabolite_key)].iloc[0];oo=row('Oslo2',r.gene,r.metabolite,r.metabolite_key);oo.update(relation_id=r.relation_id,external_feature_id=oc.external_feature_id,status='NEEDS_REVIEW' if oc.status=='DONE' else oc.status,reason='metabolite_present_but_MCOS_to_OSL2_crosswalk_missing' if oc.status=='DONE' else oc.reason);oslo.append(oo)
 a=family(result);save(a,P/'Tang_all174_associations.tsv');save(pd.DataFrame(oslo),P/'Oslo_all174_evaluability.tsv')
 save(pd.DataFrame([dict(gene=g,Tang_RNA_available=g in wide,Tang_RNA_n=int(wide[g].notna().sum()) if g in wide else 0,Oslo_RNA_status='NOT_RUN_sample_crosswalk_missing') for g in sorted(rel.gene.unique())]),P/'all117_RNA_coverage.tsv')
 out=dict(Tang_author_tumors=25,Tang_primary_samples_in_portal=len(j),Tang_RNA_cases=int(wide.shape[0]),Tang_evaluable=int(a.status.eq('DONE').sum()),Tang_q_lt005=int(a.q_value.lt(.05).sum()),Oslo_metabolomics_rows=228,Oslo_joined=0,Oslo_174_preserved=True,original_statistics_unchanged=True)
 (P/'analysis_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(out);print(a[a.gene.isin(['GPI','ASNS','GLS','GPCPD1','SLC6A8'])][['gene','metabolite_name','n','effect','ci_lower','ci_upper','q_value','status']].to_string(index=False))
if __name__=='__main__':main()
