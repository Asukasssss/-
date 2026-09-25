"""KEGG supplement using frozen LYPLA1 high/low DE, without refitting DE."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib
import pandas as pd,numpy as np,requests,gseapy
from scipy import stats
from statsmodels.stats.multitest import multipletests
BASE=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
OLD=BASE/'20260925T151104Z_lypla1_highlow_v1/public'
LIPID=['Glycerophospholipid metabolism','Glycerolipid metabolism','Fatty acid metabolism','Fatty acid biosynthesis','Fatty acid elongation','Fatty acid degradation','Biosynthesis of unsaturated fatty acids','Sphingolipid metabolism','Ether lipid metabolism','Steroid biosynthesis','Primary bile acid biosynthesis','Cholesterol metabolism','Choline metabolism in cancer','PPAR signaling pathway','Adipocytokine signaling pathway']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.10g')
def bh(p):return multipletests(p,method='fdr_bh')[1]
def main():
 R=Path(sys.argv[1]);assert (R/'.running').exists()
 for d in ['public','source','private']:(R/d).mkdir()
 P=R/'public';S=R/'source';D=R/'private'
 spec=dict(version='LYPLA1_KEGG_enrichment_v1',parent_run=str(OLD),parent_commit='090c8aa4b746871830a8e41a285491db19cac5bd',library='KEGG_2021_Human via Enrichr, fixed 2021 snapshot, not current KEGG release',contrasts=['high_vs_low_detected','high_vs_zero'],cohorts=['Wu2021','Pal2021_reprocessed'],DE='Reuse frozen full tables unchanged; no new DE/grouping',GSEA='all tested genes signed Welch t rank; gene-set permutations1000; weight1; seed20260925; sizes15..500',ORA='up/down P<0.05 and descriptive log2ratio>=0.25 or <=-0.25; full tested gene universe; hypergeometric right tail; sizes15..500',LYPLA1_excluded=True,nominalP_focus=True,multiplicity='ORA BH all up/down KEGG terms per contrast and joint4; GSEA original empirical FDR; no mixing with prior GO/Reactome q',lipid_display_terms=LIPID,software=dict(gseapy=gseapy.__version__,numpy=np.__version__,pandas=pd.__version__),limitation='Pooled-cell exploratory; inherited donor/subtype/depth imbalance; no metabolic flux or causal inference')
 (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));url='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human';r=requests.get(url,timeout=90);r.raise_for_status();f=S/'KEGG_2021_Human.gmt';f.write_bytes(r.content)
 lib={}
 for line in r.text.splitlines():
  a=line.split('\t')
  if len(a)>3:lib[a[0]]=sorted(set(x.split(',')[0] for x in a[2:] if x)-{'LYPLA1'})
 assert len(lib)>200
 manifest=[dict(path=str(f),url=url,sha256=sha(f))];enall=[];orall=[];coverage=[];checks=[];summ=[];oldhash=pd.read_csv(OLD/'DELIVERY_SHA256.tsv',sep='\t').set_index('file').sha256.to_dict()
 for co in spec['cohorts']:
  for ct in spec['contrasts']:
   stem=co+'__'+ct;source=OLD/(stem+'_DE.tsv');assert sha(source)==oldhash[source.name];manifest.append(dict(path=str(source),url='',sha256=sha(source)));de=pd.read_csv(source,sep='\t');assert de.gene.is_unique and 'LYPLA1' not in set(de.gene)
   uni=set(de.gene);rank=de[['gene','welch_t']].sort_values(['welch_t','gene'],ascending=[False,True]);filtered={}
   for term,gg in lib.items():
    g=uni.intersection(gg);ok=15<=len(g)<=500
    coverage.append(dict(cohort=co,contrast=ct,pathway=term,library_genes=len(gg),tested_overlap=len(g),status='DONE' if ok else 'NOT_EVALUABLE',reason='within_size15_to500' if ok else 'outside_size15_to500'))
    if ok:filtered[term]=g
   print('START',stem,'terms',len(filtered),flush=True)
   pre=gseapy.prerank(rnk=rank,gene_sets=lib,outdir=str(D/stem),min_size=15,max_size=500,permutation_num=1000,weighted_score_type=1,threads=4,no_plot=True,seed=20260925,verbose=False)
   en=pre.res2d.copy();en.insert(0,'contrast',ct);en.insert(0,'cohort',co);en['library']='KEGG_2021_Human';en['p_resolution_note']='zero means below permutation resolution, not exact zero';save(en,P/(stem+'_GSEA.tsv'));enall.append(en)
   # Independent ES check on evenly sampled terms.
   vals=rank.welch_t.to_numpy();names=rank.gene.to_numpy()
   for _,z in en.sort_values('Term').iloc[::max(1,len(en)//10)].iterrows():
    hit=np.isin(names,lib[z.Term]);w=np.where(hit,abs(vals),0);walk=np.cumsum(w/w.sum()-(~hit)/(~hit).sum());es=walk.max() if abs(walk.max())>abs(walk.min()) else walk.min();assert np.isclose(es,z.ES,rtol=1e-6,atol=1e-7);checks.append(dict(cohort=co,contrast=ct,type='GSEA_ES',pathway=z.Term,status='PASS'))
   rows=[]
   for direction in ['up','down']:
    sel=(de.p_value<.05)&((de.descriptive_log2ratio>=.25) if direction=='up' else (de.descriptive_log2ratio<=-.25));genes=set(de.loc[sel,'gene'])
    for term,g in filtered.items():
     hits=genes&g;pv=stats.hypergeom.sf(len(hits)-1,len(uni),len(g),len(genes));rows.append(dict(cohort=co,contrast=ct,direction=direction,pathway=term,universe_n=len(uni),selected_n=len(genes),term_n=len(g),overlap_n=len(hits),p_value=pv,genes=';'.join(sorted(hits))))
   ora=pd.DataFrame(rows);ora['q_contrast']=bh(ora.p_value);orall.append(ora)
   for _,z in ora.iloc[::max(1,len(ora)//10)].iterrows():
    N,K,n,k=map(int,[z.universe_n,z.term_n,z.selected_n,z.overlap_n]);pv=stats.fisher_exact([[k,n-k],[K-k,N-K-n+k]],alternative='greater')[1];assert np.isclose(pv,z.p_value,rtol=1e-7,atol=1e-200);checks.append(dict(cohort=co,contrast=ct,type='ORA_Fisher',pathway=z.pathway,status='PASS'))
   summ.append(dict(cohort=co,contrast=ct,n_high=int(de.n_high.iloc[0]),n_reference=int(de.n_reference.iloc[0]),tested_genes=len(de),GSEA_terms=len(en),GSEA_nominal_P_lt05=int((en['NOM p-val']<.05).sum()),ORA_terms_both_directions=len(ora),ORA_nominal_P_lt05=int((ora.p_value<.05).sum())))
   print('DONE',stem,flush=True)
 oa=pd.concat(orall,ignore_index=True);oa['q_joint4']=bh(oa.p_value)
 for (co,ct),d in oa.groupby(['cohort','contrast']):save(d,P/(co+'__'+ct+'_ORA.tsv'))
 ea=pd.concat(enall,ignore_index=True);save(ea,P/'GSEA_all4.tsv');save(oa,P/'ORA_all4.tsv');save(pd.DataFrame(coverage),P/'pathway_coverage.tsv');save(pd.DataFrame(summ),P/'analysis_summary.tsv');save(pd.DataFrame(checks),P/'independent_enrichment_audit.tsv')
 save(ea[ea.Term.isin(LIPID)],P/'lipid_GSEA_all4.tsv');save(oa[oa.pathway.isin(LIPID)],P/'lipid_ORA_all4.tsv')
 # Join only matching pathway, direction and primary contrast; preserve non-support.
 for mode,all_,key,cols in [('GSEA',ea,['Term'],['NES','NOM p-val','FDR q-val','Lead_genes']),('ORA',oa,['pathway','direction'],['p_value','q_contrast','overlap_n','genes'])]:
  primary=all_[all_.contrast.eq('high_vs_low_detected')];a=primary[primary.cohort.eq('Wu2021')];b=primary[primary.cohort.eq('Pal2021_reprocessed')];c=a[key+cols].merge(b[key+cols],on=key,how='outer',suffixes=('_Wu','_Pal'));save(c,P/('cross_cohort_primary_'+mode+'.tsv'))
 manifest.append(dict(path=str(Path(__file__)),url='',sha256=sha(__file__)));save(pd.DataFrame(manifest),P/'source_manifest.tsv')
 (P/'validation.json').write_text(json.dumps(dict(status='PASS',old_DE_hashes_verified=4,no_DE_refit=True,target_excluded=True,independent_enrichment_checks=len(checks),audit_limit='Selected ES and ORA P checked; permutation NES/FDR not independently rerun'),indent=2));(R/'NUMERICAL_DONE').write_text('DONE')
 print('ALL_DONE',flush=True)
if __name__=='__main__':main()
