"""Author processed log2 GCRMA; all117 annotation; unpaired Welch RNA context."""
from pathlib import Path
import sys,json,io,hashlib
import numpy as np,pandas as pd
from scipy import stats
from brca_oslo_tang_analysis_v1 import row,save,family
r=Path(sys.argv[1]);s=r/'source';p=r/'public'
assert (p/'analysis_spec.json').exists()
def table(text,kind):
 return pd.read_csv(io.StringIO(text.split('!'+kind+'_table_begin')[1].split('!'+kind+'_table_end')[0].strip()),sep='\t',dtype=str)
platform=table((s/'GPL570_self_full.soft').read_text(),'platform');assert platform.ID.is_unique
entrez=platform.ENTREZ_GENE_ID.fillna('').str.strip();valid=entrez.str.fullmatch(r'\d+');ann=platform[valid].copy();ann['entrez']=entrez[valid]
rel=pd.read_csv(s/'relations174.tsv',sep='\t');genes=sorted(rel.gene.unique());gm={};coverage=[]
for g in genes:
 j=json.loads((s/('tcga_gene_'+g+'.json')).read_text());assert j['hugoGeneSymbol']==g
 probes=ann.loc[ann.entrez.eq(str(j['entrezGeneId'])),'ID'].tolist();gm[g]=probes
 coverage.append(dict(gene=g,entrez_id=j['entrezGeneId'],n_unambiguous_probes=len(probes),probe_ids=';'.join(probes),status='DONE' if probes else 'NOT_EVALUABLE',reason='median_all_probes_unique_Entrez_mapping;no_outcome_selection' if probes else 'no_unique_Entrez_probes'))
save(pd.DataFrame(coverage),p/'GSE42568_probe_coverage.tsv')
allids=[l.split(' = ')[-1] for l in (s/'GSE42568_series.soft').read_text().splitlines() if l.startswith('!Series_sample_id')];assert len(allids)==121 and len(set(allids))==121
records=[];meta=[]
for acc in allids:
 text=(s/(acc+'_self_full.soft')).read_text();lines=text.splitlines();title=[l.split(' = ',1)[1] for l in lines if l.startswith('!Sample_title')][0];tissue=[l.split('tissue: ',1)[1] for l in lines if l.startswith('!Sample_characteristics_ch1 = tissue: ')];assert len(tissue)==1
 group='normal' if tissue[0]=='normal breast' else 'tumor' if tissue[0]=='breast cancer' else 'UNRESOLVED'
 meta.append(dict(sample=acc,title=title,tissue=tissue[0],group=group));a=table(text,'sample').set_index('ID_REF')['VALUE'].apply(pd.to_numeric,errors='coerce');assert a.index.is_unique
 records.append(dict(sample=acc,**{g:float(a.reindex(probes).median()) if probes else np.nan for g,probes in gm.items()}))
md=pd.DataFrame(meta);print(md.groupby(['tissue','group']).size().to_string(),flush=True);assert md.title.is_unique and set(md.group)=={'normal','tumor'}
assert md.group.eq('normal').sum()==17 and md.group.eq('tumor').sum()==104
wide=pd.DataFrame(records).set_index('sample');save(wide.reset_index(),s/'GSE42568_candidate_expression_private.tsv');save(md,s/'GSE42568_sample_metadata_private.tsv')
rows=[]
for g in genes:
 a=wide.loc[md[md.group.eq('tumor')]['sample'],g].dropna().values;b=wide.loc[md[md.group.eq('normal')]['sample'],g].dropna().values
 rr=row('GSE42568',g,kind='unpaired_Welch_author_log2GCRMA');rr.update(unit='author_sample;no_pairing_inferred',effect_type='mean_log2_expression_tumor_minus_normal',test_family='GSE42568_all117_RNA',n=len(a),n_reference=len(b),n_probes=len(gm[g]),reason='no_unambiguous_probe_or_insufficient_values')
 if len(a)>2 and len(b)>2:
  v1=a.var(ddof=1)/len(a);v2=b.var(ddof=1)/len(b);se=np.sqrt(v1+v2);df=(v1+v2)**2/(v1*v1/(len(a)-1)+v2*v2/(len(b)-1));effect=a.mean()-b.mean();ci=stats.t.ppf(.975,df)*se
  if se>0:rr.update(status='DONE',reason='unadjusted_RNA_context_only;not_metabolite_relation_validation',effect=effect,ci_lower=effect-ci,ci_upper=effect+ci,p_value=stats.ttest_ind(a,b,equal_var=False).pvalue,welch_df=df)
 rows.append(rr)
out=family(rows);save(out,p/'GSE42568_all117_RNA.tsv');summary=dict(author_samples=121,tumor=104,normal=17,evaluable=int(out.status.eq('DONE').sum()),q_lt005=int(out.q_value.lt(.05).sum()),no_pairing_inferred=True,scope='RNA context not metabolite association or function');(p/'GSE42568_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
