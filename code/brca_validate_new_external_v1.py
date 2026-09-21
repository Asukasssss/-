"""Independent effect/BH and source-identifier checks; public summaries only."""
from pathlib import Path
import sys,json,hashlib,re
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
r=Path(sys.argv[1]);s=r/'source';p=r/'public'
a=pd.read_csv(p/'Tang_all174_associations.tsv',sep='\t');o=pd.read_csv(p/'Oslo_all174_evaluability.tsv',sep='\t');g=pd.read_csv(p/'GSE42568_all117_RNA.tsv',sep='\t');assert len(a)==len(o)==174 and len(g)==117 and a.relation_id.is_unique and g.gene.is_unique
err={}
for name,d in [('Tang',a),('RNA',g)]:
 ix=d.status.eq('DONE');err[name+'_BH_max_error']=float(np.max(abs(multipletests(d.loc[ix,'p_value'],method='fdr_bh')[1]-d.loc[ix,'q_value'])));assert err[name+'_BH_max_error']<1e-12
d=pd.read_excel(s/'Tang_S2.xlsx',header=None);cols={('TCGA-'+str(d.iloc[0,k]).strip()):k for k in range(3,33) if re.fullmatch(r'[A-Z0-9]{2}-[A-Z0-9]{4}',str(d.iloc[0,k]).strip())};ids={str(int(d.iloc[i,0])):i for i in range(19,len(d))}
raw=[x for f in s.glob('tcga_rna_*.json') for x in json.loads(f.read_text())];z=pd.DataFrame([dict(patient=x['patientId'],sample=x['sampleId'],gene=x['gene']['hugoGeneSymbol'],value=x['value']) for x in raw]);assert not z.duplicated(['patient','gene']).any();assert z.groupby('patient')['sample'].nunique().max()==1
errors=[]
for x in a[a.status.eq('DONE')].itertuples():
 key=str(int(x.external_feature_id));v={pat:float(d.iloc[ids[key],col]) for pat,col in cols.items()};zz=z[z.gene.eq(x.gene)].copy();zz['met']=zz.patient.map(v);zz=zz.dropna();assert len(zz)==x.n;errors.append(abs(stats.spearmanr(zz.met,zz.value).statistic-x.effect))
err['Tang_rho_max_error']=max(errors);assert max(errors)<1e-12
e=pd.read_csv(s/'GSE42568_candidate_expression_private.tsv',sep='\t').set_index('sample');md=pd.read_csv(s/'GSE42568_sample_metadata_private.tsv',sep='\t');assert md['sample'].is_unique and md.title.is_unique
errors=[]
for x in g[g.status.eq('DONE')].itertuples():
 t=e.loc[md.loc[md.group.eq('tumor'),'sample'],x.gene].dropna();n=e.loc[md.loc[md.group.eq('normal'),'sample'],x.gene].dropna();assert len(t)==x.n and len(n)==x.n_reference
 errors.extend([abs(t.mean()-n.mean()-x.effect),abs(stats.ttest_ind(t,n,equal_var=False).pvalue-x.p_value)])
err['RNA_effect_and_P_max_error']=max(errors);assert max(errors)<1e-10
res=dict(status='DONE',Tang_rows174=True,Oslo_rows174=True,RNA_rows117=True,patient_order_not_used=True,no_duplicate_patient_gene=True,checks=err,not_recomputed='permutation P and bootstrap CI not regenerated;no spectral reidentification;no Oslo crosswalk;patient aliquot equivalence not proved')
(p/'validation.json').write_text(json.dumps(res,indent=2)+'\n')
manifest=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted(s.iterdir()) if f.is_file()]
pd.DataFrame(manifest).to_csv(p/'source_manifest.tsv',sep='\t',index=False,lineterminator='\n')
print(json.dumps(res,indent=2))
