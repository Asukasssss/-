"""Recompute selected GSEA ES and ORA probabilities independently."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy.stats import fisher_exact
R=Path(sys.argv[1]);P=R/'public';lib={}
for line in (R/'source/Reactome_2022.gmt').read_text().splitlines():
 a=line.split('\t');lib[a[0]]=set(x.split(',')[0] for x in a[2:] if x)-{'LYPLA1'}
checks=[]
for co in ['Wu2021','Pal2021_reprocessed']:
 for ct in ['high_vs_low_detected','high_vs_zero']:
  stem=co+'__'+ct;de=pd.read_csv(P/(stem+'_DE.tsv'),sep='\t');en=pd.read_csv(P/(stem+'_GSEA.tsv'),sep='\t');ora=pd.read_csv(P/(stem+'_ORA.tsv'),sep='\t');de=de.sort_values(['welch_t','gene'],ascending=[False,True]);rank=de.welch_t.to_numpy();names=de.gene.to_numpy();assert 'LYPLA1' not in set(names)
  for _,r in en.sort_values('Term').iloc[::max(1,len(en)//12)].iterrows():
   hit=np.isin(names,list(lib[r.Term]));w=np.where(hit,abs(rank),0);running=np.cumsum(w/w.sum()-(~hit)/(~hit).sum());maximum=running.max();minimum=running.min();es=maximum if abs(maximum)>abs(minimum) else minimum
   assert np.isclose(es,r.ES,rtol=1e-6,atol=1e-7),(stem,r.Term,es,r.ES)
   checks.append(dict(cohort=co,contrast=ct,test='GSEA_ES',term=r.Term,status='PASS'))
  for _,r in ora.sort_values(['library','term','direction']).iloc[::max(1,len(ora)//15)].iterrows():
   N,K,n,k=map(int,[r.universe_n,r.term_n,r.selected_n,r.overlap_n]);pv=fisher_exact([[k,n-k],[K-k,N-K-n+k]],alternative='greater')[1]
   assert np.isclose(pv,r.p_value,rtol=1e-7,atol=1e-200)
   checks.append(dict(cohort=co,contrast=ct,test='ORA_Fisher',term=r.term,status='PASS'))
pd.DataFrame(checks).to_csv(P/'independent_enrichment_audit.tsv',sep='\t',index=False)
j=json.loads((P/'validation.json').read_text());j['independent_enrichment_checks']=len(checks);j['enrichment_audit']='PASS selected ES and ORA exact P; permutation NES/FDR not independently rerun';(P/'validation.json').write_text(json.dumps(j,indent=2));print('PASS',len(checks),'enrichment checks')
