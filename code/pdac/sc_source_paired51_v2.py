"""All 538 candidates; reuse identical donor summaries, add new genes and donor bootstrap."""
import argparse,hashlib,json,traceback,platform
from pathlib import Path
import numpy as np
import pandas as pd
import sc_source_three_cohorts_v1 as old
from internal_paired51_v2 import resolve_genes
ROOT=old.ROOT
PREVIOUS=ROOT/'results/collaborative/PDAC/B/20260921T151600Z_sc_source_three_cohorts_v1'
COHORTS=['GSE263733','GSE278688','GSE242230']
# These historical synonyms also name other reviewed human canonical genes.
# Do not resolve a missing canonical target to another gene's current symbol.
ALIAS_COLLISIONS={'AK4':['AK3'],'AOC1':['DAO'],'AQP7':['AQP9'],'CES1':['CES2'],'KDM5A':['RBP2'],'PEMT':['PNMT'],'SLC25A15':['ORC1'],'SLC25A2':['ORC2'],'SLC38A2':['SAT2'],'SLC6A20':['SIT1'],'SLC7A8':['LAT2']}
def sha(p):return old.sha(Path(p))
def save(p,v):old.save(p,v)
def summarize_pb(pb,cohort,genes,seed):
 pb=pb[pb.scheme=='broad'].copy();units=sorted(pb.unit.unique());rng=np.random.default_rng(seed);weights=rng.multinomial(len(units),np.full(len(units),1/len(units)),size=1000);summaries=[];tops=[]
 for gene in genes:
  d=pb[pb.gene==gene];cats=sorted(d.celltype.unique());means=[];eligible=[];mat=[]
  for ct in cats:
   a=d[(d.celltype==ct)&(d.n_cells>=20)];a=a[a.log1p_cpm.notna()];nd=int((a.sum_counts>0).sum());status='EVALUABLE' if len(a)>=3 and nd>=3 else 'GENE_ABSENT_OR_AMBIGUOUS' if d.log1p_cpm.notna().sum()==0 else 'INSUFFICIENT_UNITS_OR_DETECTION'
   mean=float(a.log1p_cpm.mean()) if len(a) else np.nan;means.append(mean);eligible.append(status=='EVALUABLE');mat.append(a.set_index('unit').log1p_cpm.reindex(units).to_numpy(float))
   summaries.append({'cohort':cohort,'gene':gene,'celltype':ct,'mean_expression':mean,'detection_fraction':float(a.positive_fraction.mean()) if len(a) else np.nan,'donor_count':len(a),'detected_donor_count':nd,'status':status})
  x=np.array(mat).T;mask=np.isfinite(x);den=weights@mask.astype(float);num=weights@np.nan_to_num(x);boot=np.full(den.shape,np.nan);np.divide(num,den,out=boot,where=den>0);boot[:,~np.array(eligible)]=np.nan
  indices=[i for i,e in enumerate(eligible) if e];top=[];second=[];freq=np.full(len(cats),np.nan);valid=np.isfinite(boot).any(axis=1)
  if indices:
   ordered=sorted(indices,key=lambda i:(-means[i],cats[i]));top=[i for i in ordered if np.isclose(means[i],means[ordered[0]],rtol=0,atol=1e-12)];others=[i for i in ordered if i not in top]
   if others:second=[i for i in others if np.isclose(means[i],means[others[0]],rtol=0,atol=1e-12)]
   b=boot[valid];mx=np.nanmax(b,axis=1);wins=np.isclose(b,mx[:,None],rtol=0,atol=1e-12);share=wins/wins.sum(axis=1)[:,None];freq=share.mean(axis=0)
   for j in range(len(cats)):
    if not eligible[j]:freq[j]=np.nan
  topct=';'.join(cats[i] for i in top) or 'NOT_EVALUABLE';found=next((r for r in summaries[::-1] if r['gene']==gene and r['celltype']==topct),None)
  tops.append({'cohort':cohort,'gene':gene,'top_celltype':topct,'second_celltype':';'.join(cats[i] for i in second) or 'NOT_EVALUABLE','mean_expression':found['mean_expression'] if found else np.nan,'detection_fraction':found['detection_fraction'] if found else np.nan,'donor_count':found['donor_count'] if found else 0,'bootstrap_top_frequency':float(freq[top[0]]) if len(top)==1 else np.nan,'bootstrap_valid_draws':int(valid.sum()),'eligible_categories':';'.join(cats[i] for i in indices),'seed':seed})
  for row,f in zip(summaries[-len(cats):],freq):row['bootstrap_first_frequency']=f
 return pd.DataFrame(summaries),pd.DataFrame(tops)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--code-commit',required=True);args=ap.parse_args();run=Path(__file__).resolve().parent;assert run.parent==ROOT/'results/collaborative/PDAC/B'
 with (run/'.running').open('x') as f:f.write('538 candidate sources')
 try:
  public=run/'public';public.mkdir(exist_ok=False);panel=json.loads((run/'panel.json').read_text());genes=panel['genes'];assert len(genes)==538;aliases=json.loads((run/'gene_aliases.json').read_text());outputs=[];tops=[];checks=[];resolutions=[];extra_sources=[]
  for g,bad in ALIAS_COLLISIONS.items():aliases[g]=[a for a in aliases.get(g,[]) if a not in bad]
  prior_manifest=pd.read_csv(PREVIOUS/'public/source_manifest.tsv',sep='\t');expected=dict(zip(prior_manifest.path,prior_manifest.sha256))
  for cohort,reader in zip(COHORTS,[old.dense263,old.sparse278,old.sparse242]):
   print('START',cohort,flush=True);prior_path=PREVIOUS/(cohort+'_private_unit_pseudobulk.tsv');prior=pd.read_csv(prior_path,sep='\t');extra_sources.append(prior_path)
   reused=sorted(set(prior.loc[prior.log1p_cpm.notna(),'gene'])&set(genes));new=sorted(set(genes)-set(reused));sentinels=sorted(set(reused)&{'MGLL','PNP','PEPD','GGT5'});search=sorted(set(new)|{s for g in new for s in aliases.get(g,[])}|set(sentinels))
   meta,target,lib,present=reader(search);resolved,why=resolve_genes(new,present,aliases)
   # Re-read four retained genes as controls for cached donor-summary reuse.
   for gene in sentinels:
    gi=search.index(gene);assert gene in present
    for _,row in prior[prior.gene==gene].iterrows():
     celllabels=meta.celltype.map(old.broad) if row.scheme=='broad' else meta.celltype
     idx=np.flatnonzero((meta.unit.astype(str)==str(row.unit))&(celllabels==row.celltype));values=target[gi,idx]
     assert len(idx)==row.n_cells and int(values.sum())==row.sum_counts and int(lib[idx].sum())==row.library_umi
     assert abs(float((values>0).mean())-row.positive_fraction)<1e-12
   # Prevent an alias for a new canonical target from reusing an already represented gene.
   for g,l in list(resolved.items()):
    if l in reused:resolved[g]=None;why[g]='COLLISION_WITH_REUSED_CANONICAL'
   new_present={g for g in new if resolved[g]};lookup={g:i for i,g in enumerate(search)};a=np.zeros((len(new),len(meta)),dtype=np.int64)
   for i,g in enumerate(new):
    if resolved[g]:a[i]=target[lookup[resolved[g]]]
   del target
   new_out,cells,check=old.summarize(cohort,meta,a,lib,new_present,new,run);newpb=pd.read_csv(run/(cohort+'_private_unit_pseudobulk.tsv'),sep='\t');allpb=pd.concat([prior[prior.gene.isin(reused)],newpb],ignore_index=True)
   assert set(allpb.gene)==set(genes) and not allpb.duplicated(['scheme','unit','celltype','gene']).any()
   # The previous units, cell groups and denominators must agree with the current extraction.
   k=['scheme','unit','celltype'];base=prior[k+['n_cells','library_umi']].drop_duplicates();now=newpb[k+['n_cells','library_umi']].drop_duplicates();pd.testing.assert_frame_equal(base.sort_values(k).reset_index(drop=True),now.sort_values(k).reset_index(drop=True),check_dtype=False)
   allpb.to_csv(run/(cohort+'_private_unit_pseudobulk.tsv'),sep='\t',index=False,na_rep='NA')
   check.update(target_genes=538,genes_reused=len(reused),genes_new_or_alias_rechecked=len(new),genes_available=len(reused)+len(new_present),missing_genes=sorted(set(new)-new_present),reused_gene_sentinels_checked=sentinels);checks.append(check)
   for g in genes:resolutions.append({'cohort':cohort,'gene':g,'matrix_label':g if g in reused else resolved[g],'method':'REUSED_IDENTICAL_DONOR_SUMMARIES' if g in reused else why[g]})
   seed=int(hashlib.sha256(('PDAC_sc_paired51_v2|'+cohort).encode()).hexdigest()[:8],16);agg,top=summarize_pb(allpb,cohort,genes,seed);outputs.append(agg);tops.append(top);cells.to_csv(public/(cohort+'_cell_counts.tsv'),sep='\t',index=False)
   print('DONE',cohort,json.dumps(check),flush=True)
  allrows=pd.concat(outputs,ignore_index=True);toprows=pd.concat(tops,ignore_index=True);consensus=[]
  for gene in genes:
   d=toprows[toprows.gene==gene].set_index('cohort');validtops=[d.loc[c,'top_celltype'] for c in COHORTS];sets=[set(allrows.loc[(allrows.gene==gene)&(allrows.cohort==c)&(allrows.status=='EVALUABLE'),'celltype']) for c in COHORTS];common=set.intersection(*sets);shared=[]
   for c in COHORTS:
    sub=allrows[(allrows.gene==gene)&(allrows.cohort==c)&allrows.celltype.isin(common)]
    st=';'.join(sorted(sub.loc[np.isclose(sub.mean_expression,sub.mean_expression.max(),rtol=0,atol=1e-12),'celltype'])) if len(common)>=2 else 'NOT_EVALUABLE';shared.append(st)
   same=len(set(validtops))==1 and validtops[0]!='NOT_EVALUABLE' and ';' not in validtops[0];sharedsame=len(common)>=2 and len(set(shared))==1 and ';' not in shared[0]
   stable=bool(same and (d.bootstrap_top_frequency>=.7).all() and (d.detection_fraction>=.05).all())
   row={'gene':gene,'cross_dataset_same_top':same,'shared_category_same_top':sharedsame if len(common)>=2 else 'NOT_EVALUABLE','n_shared_categories':len(common),'shared_categories':';'.join(sorted(common)),'stable_source':stable,'consensus_top':validtops[0] if same else 'COHORT_DEPENDENT_OR_NOT_EVALUABLE'}
   for c,st in zip(COHORTS,shared):
    for key in ['top_celltype','second_celltype','mean_expression','detection_fraction','donor_count','bootstrap_top_frequency']:row[c+'_'+key]=d.loc[c,key]
    row[c+'_shared_top']=st
   consensus.append(row)
  allrows.to_csv(public/'cell_source_all.tsv',sep='\t',index=False,na_rep='NA');toprows.to_csv(public/'cohort_top_sources.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(consensus).to_csv(public/'cross_cohort_source.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(resolutions).to_csv(public/'gene_resolution.tsv',sep='\t',index=False,na_rep='NA')
  for p in old.SOURCE_FILES:
   if str(p) in expected and sha(p)!=expected[str(p)]:raise ValueError('Historical source changed: '+str(p))
  extra_sources += [run/f for f in ['panel.json','gene_aliases.json','sc_source_three_cohorts_v1.py','internal_paired51_v2.py','patient_first_round.py','paired_metabolites_v1.py']]+[Path(__file__),PREVIOUS/'public/source_manifest.tsv']
  pd.DataFrame([{'path':str(p),'sha256':sha(p)} for p in sorted(old.SOURCE_FILES|set(extra_sources))]).to_csv(public/'source_manifest.tsv',sep='\t',index=False)
  summary={'candidate_genes':538,'cohorts':checks,'three_cohort_same_top':sum(r['cross_dataset_same_top'] for r in consensus),'stable_source_genes':sum(r['stable_source'] for r in consensus),'shared_category_same_top':sum(r['shared_category_same_top'] is True for r in consensus),'source_rows':len(allrows)}
  save(public/'summary.json',summary);save(public/'analysis_spec.json',{'code_commit':args.code_commit,'protocol':'docs/PDAC/PAIRED51_INTERNAL_V2_LOCK.md','normalization':'donor-celltype log1p(sum_gene_UMI*1e6/sum_all_gene_UMI);equal donor mean','bootstrap':'1000 joint donor draws;ties fractional;no P/q','minimum_cells':20,'minimum_units':3,'minimum_detected_units':3,'shared_categories':'Intersection of gene-eligible categories across all3 cohorts;at least2','source_rule':'All3 same unique top,each bootstrap>=0.70 and detection>=0.05','canonical_alias_policy':'Exact canonical first;unique reviewed alias only;collision excluded','versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}})
  save(public/'validation.json',{'status':'PASS','all538_in_each_cohort':True,'historical_source_hashes_match':True,'reused_unit_cell_groups_denominators_match_new_extraction':True,'no_patient_rows_exported':True,'unverified':['Malignant CNV reannotation','Clinical patient-level cross-study identity de-duplication','Ambient RNA correction','Disease differential expression','Metabolite-gene external validation']})
  save(run/'DONE.json',summary);(run/'.running').unlink()
 except Exception:(run/'FAILED.txt').write_text(traceback.format_exc());raise
if __name__=='__main__':main()
