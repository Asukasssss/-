"""Three fixed descriptive case-leave-out checks; scatter and case rows stay server-only."""
from pathlib import Path
import sys,json,re,hashlib
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(sys.argv[1]);old=r.parent/'20260921T072503Z_oslo_tang_external_v1';s=old/'source';p=r/'public';priv=r/'private'
pairs=[('MDH1','KEGG:C00149'),('PNP','KEGG:C00262'),('PNP','KEGG:C00242')]
spec=dict(source_commit='50dfd7b',relations=[dict(gene=g,metabolite_key=k) for g,k in pairs],method='omit each of20cases once;Spearman rho only;no P/q or CI recalculation;no sample exclusion from main result',selection='post_Tang_results;descriptive_only',scatter='raw author-normalized-imputed metabolomics vs RSEM and ranks;private server only',public='min max max_abs_change sign_retention;no individual values or patient identifiers')
(p/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
d=pd.read_excel(s/'Tang_S2.xlsx',header=None);cols={('TCGA-'+str(d.iloc[0,k]).strip()):k for k in range(3,33) if re.fullmatch(r'[A-Z0-9]{2}-[A-Z0-9]{4}',str(d.iloc[0,k]).strip())};ids={str(int(d.iloc[i,0])):i for i in range(19,len(d))}
raw=[x for f in sorted(s.glob('tcga_rna_*.json')) for x in json.loads(f.read_text())];z=pd.DataFrame([dict(patient=x['patientId'],gene=x['gene']['hugoGeneSymbol'],value=x['value']) for x in raw]);assert not z.duplicated(['patient','gene']).any()
a=pd.read_csv(old/'public/Tang_all174_associations.tsv',sep='\t');summ=[];detail=[];fig,ax=plt.subplots(3,2,figsize=(10,11));maxerr=0
for i,(gene,key) in enumerate(pairs):
 t=a[a.gene.eq(gene)&a.metabolite_key.eq(key)];assert len(t)==1;t=t.iloc[0];source_id=str(int(t.external_feature_id));met={pat:float(d.iloc[ids[source_id],col]) for pat,col in cols.items()};zz=z[z.gene.eq(gene)].copy();zz['met']=zz.patient.map(met);zz=zz.dropna().sort_values('patient');assert len(zz)==t.n==20
 x=zz.value.to_numpy();y=zz.met.to_numpy();base=float(stats.spearmanr(x,y).statistic);maxerr=max(maxerr,abs(base-t.effect));assert abs(base-t.effect)<1e-12
 rr=[]
 for j,patient in enumerate(zz.patient):
  mask=np.arange(len(x))!=j;val=float(stats.spearmanr(x[mask],y[mask]).statistic);rr.append(val);detail.append(dict(relation_id=t.relation_id,omitted_patient=patient,n=19,rho=val))
 rr=np.array(rr);assert np.isfinite(rr).all()
 summ.append(dict(cancer='BRCA',cohort='Tang2014',stage_id='04_ROBUSTNESS',run_id=r.name,analysis_version='tang_three_loo_v1',analysis_type='descriptive_case_leave_one_out',metabolite_key=key,metabolite_name=t.metabolite_name,gene=gene,unit='author_TCGA_case',n=20,n_reference=np.nan,effect_type='original_full_sample_Spearman_rho',effect=base,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family='NA_no_new_inference',family_n_evaluable=0,status='DONE',reason='post_selection_descriptive;no_P_q_recomputed;20_leaveouts_are_not20_independent_replications',source_id='Tang2014;TCGA2015',relation_id=t.relation_id,original_p=t.p_value,original_q=t.q_value,original_q_family='Tang163_evaluable',original_ci_lower=t.ci_lower,original_ci_upper=t.ci_upper,loo_n=19,n_leaveouts=20,loo_rho_min=rr.min(),loo_rho_max=rr.max(),loo_rho_median=np.median(rr),max_abs_rho_change=np.max(abs(rr-base)),sign_retained=int((np.sign(rr)==np.sign(base)).sum()),RNA_unique_values=len(np.unique(x)),metabolite_unique_values=len(np.unique(y))))
 for j,(xx,yy,label) in enumerate([(x,y,'Author scales'),(stats.rankdata(x),stats.rankdata(y),'Ranks')]):
  ax[i,j].scatter(xx,yy,color='#24689b',s=30);ax[i,j].set_title(f'{gene} / {t.metabolite_name}: {label}');ax[i,j].set_xlabel('RNA RSEM' if j==0 else 'RNA rank');ax[i,j].set_ylabel('Author metabolite normalized/imputed' if j==0 else 'Metabolite rank')
fig.tight_layout();fig.savefig(priv/'three_relationship_scatter_PRIVATE.png',dpi=160);plt.close(fig)
pd.DataFrame(detail).to_csv(priv/'leaveout_case_detail_PRIVATE.tsv',sep='\t',index=False)
out=pd.DataFrame(summ);out.to_csv(p/'three_relationship_LOO_summary.tsv',sep='\t',index=False,na_rep='NA',lineterminator='\n')
# Aggregate ranges only; no individual patient points exported.
fig,ax=plt.subplots(figsize=(8,3.3))
for i,t in out.iterrows():ax.plot([t.loo_rho_min,t.loo_rho_max],[i,i],lw=5,color='#76a5c6');ax.plot(t.effect,i,'o',color='#143e59');ax.text(.99,i,f'{t.sign_retained}/20 positive',ha='right',va='center',fontsize=9)
ax.set_xlim(0,1);ax.set_yticks(range(3));ax.set_yticklabels([g+' / '+m for g,m in out[['gene','metabolite_name']].values]);ax.invert_yaxis();ax.set_xlabel('Spearman rho');ax.set_title('Case leave-one-out: ranges of 20 recalculated effects\nDots = original n20 effect; bars are NOT confidence intervals');fig.tight_layout();fig.savefig(p/'three_relationship_LOO_ranges.png',dpi=180);plt.close(fig)
files=[s/'Tang_S2.xlsx',old/'public/Tang_all174_associations.tsv']+sorted(s.glob('tcga_rna_*.json'))
pd.DataFrame([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(p/'source_manifest.tsv',sep='\t',index=False,lineterminator='\n')
(p/'validation.json').write_text(json.dumps(dict(status='DONE',fixed_relations=3,leaveouts_per_relation=20,full_rho_max_error=maxerr,original_q_unchanged=True,new_P_q_tests=0,scatter_and_individual_details_server_only=True),indent=2)+'\n')
print(out[['gene','metabolite_name','effect','loo_rho_min','loo_rho_max','max_abs_rho_change','sign_retained','RNA_unique_values','metabolite_unique_values']].to_string(index=False))
