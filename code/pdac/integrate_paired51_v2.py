"""Aggregate-only integration; no significance gate for candidate retention."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
MAPPING=ROOT/'results/PDAC/02_MAPPING/20260922T055000Z_paired51_mapping_v2'
INTERNAL=ROOT/'results/PDAC/03_PATIENT/20260922T053500Z_internal_paired51_v2'
SC=ROOT/'results/PDAC/06_EXTERNAL/20260922T054000Z_sc_paired51_v2'
PAIRED=ROOT/'results/PDAC/04_ROBUSTNESS/20260921T121509Z_paired_metabolites_v1'
DIRECTION=ROOT/'results/PDAC/04_ROBUSTNESS/20260921T122353Z_pair_direction_counts_v1'
def relation_label(main,avail):
 if main['status']!='DONE':return 'R-not-evaluable'
 sensitive=main['loo_sign_changes']>0
 if avail['status']=='DONE':sensitive|=np.sign(main['effect'])!=np.sign(avail['effect']) or abs(main['effect']-avail['effect'])>.20
 if sensitive:return 'R-sensitive'
 if main['q_value']<.05 and avail['status']=='DONE':return 'R-strong'
 if main['p_value']<.05 or abs(main['effect'])>=.4:return 'R-suggestive'
 return 'R-weak'
def write(path,df):df.to_csv(path,sep='\t',index=False,na_rep='NA')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-id',required=True);ap.add_argument('--append-sc',action='store_true');args=ap.parse_args();out=ROOT/'results/PDAC/07_INTEGRATION'/args.run_id;out.mkdir(parents=True,exist_ok=True)
 planned=pd.read_csv(MAPPING/'planned_relations.tsv',sep='\t');features=pd.read_csv(MAPPING/'feature_coverage_51.tsv',sep='\t');met=pd.read_csv(PAIRED/'results.tsv',sep='\t');di=pd.read_csv(DIRECTION/'direction_counts.tsv',sep='\t');assoc=pd.read_csv(INTERNAL/'associations.tsv',sep='\t');rna=pd.read_csv(INTERNAL/'paired_RNA.tsv',sep='\t').set_index('gene');main=assoc[assoc.analysis_type=='primary'].set_index('relation_id');av=assoc[assoc.analysis_type=='availability'].set_index('relation_id');rows=[]
 for _,f in features.iterrows():
  p=met[(met.metabolite_name==f.feature_name)&(met.analysis_type=='primary')].iloc[0];s=met[(met.metabolite_name==f.feature_name)&(met.analysis_type=='both_preimputation_available')].iloc[0];d=di[(di.feature_name==f.feature_name)&(di.analysis_type=='primary')].iloc[0]
  base={'cancer':'PDAC','metabolite_name':f.feature_name,'metabolite_key':f.metabolite_key,'metabolite_paired_n':p.n,'metabolite_n_up':d.n_up,'metabolite_n_down':d.n_down,'metabolite_n_equal':d.n_equal,'metabolite_direction':d.majority,'metabolite_paired_mean_delta':p.effect,'metabolite_paired_median_delta':p.median_delta,'metabolite_paired_ci_lower':p.ci_lower,'metabolite_paired_ci_upper':p.ci_upper,'metabolite_paired_rank_biserial':p.rank_biserial,'metabolite_paired_p':p.p_value,'metabolite_paired_q':p.q_value,'metabolite_evidence':'NOMINAL_EXPLORATORY','original_camp_q':p.original_camp_q,'metabolite_available_pairs':s.n,'metabolite_both_available_pair_fraction':s.n/11,'metabolite_available_effect':s.effect,'metabolite_available_p':s.p_value,'metabolite_available_q':s.q_value,'metabolite_available_status':s.status,'metabolite_available_direction_agrees':bool(np.sign(s.effect)==np.sign(p.effect)) if s.status=='DONE' else 'NOT_EVALUABLE','identity_note':f.chemical_identity_note,'mapping_disposition':f.mapping_status}
  rels=planned[planned.feature_name==f.feature_name]
  if not len(rels):rows.append({**base,'gene':'UNRESOLVED','relation_id':'UNRESOLVED:'+str(f.metabolite_key),'mapping_tier':'UNRESOLVED','relation_type':'UNRESOLVED','mapping_reason':f.mapping_status,'relation_label':'R-not-evaluable'});continue
  for _,r in rels.iterrows():
   a=main.loc[r.relation_id];b=av.loc[r.relation_id];e=rna.loc[r.gene];row={**base,'gene':r.gene,'relation_id':r.relation_id,'mapping_tier':r.mapping_tier,'relation_type':r.relation_type,'mapping_reason':r.mapping_reason,'mapping_source_urls':r.source_urls,'reaction_roles':r.reaction_roles,'compartment':r.compartment,'relation_label':relation_label(a,b),'CAMP_n':a.n,'CAMP_rho':a.effect,'CAMP_ci_lower':a.ci_lower,'CAMP_ci_upper':a.ci_upper,'CAMP_p':a.p_value,'CAMP_q':a.q_value,'CAMP_status':a.status,'CAMP_family':a.test_family,'CAMP_family_n':a.family_n_planned,'CAMP_available_n':b.n,'CAMP_available_rho':b.effect,'CAMP_available_ci_lower':b.ci_lower,'CAMP_available_ci_upper':b.ci_upper,'CAMP_available_p':b.p_value,'CAMP_available_q':b.q_value,'CAMP_available_status':b.status,'CAMP_available_fraction':b.available_fraction,'CAMP_rho_available_abs_change':abs(a.effect-b.effect),'CAMP_available_same_input':a.input_hash==b.input_hash,'CAMP_loo_rho_min':a.loo_rho_min,'CAMP_loo_rho_max':a.loo_rho_max,'CAMP_loo_max_abs_change':a.loo_max_abs_change,'CAMP_loo_sign_changes':a.loo_sign_changes,'RNA_n':e.n,'RNA_n_up':e.n_up,'RNA_n_down':e.n_down,'RNA_n_equal':e.n_equal,'RNA_mean_delta':e.effect,'RNA_median_delta':e.median_delta,'RNA_ci_lower':e.ci_lower,'RNA_ci_upper':e.ci_upper,'RNA_p':e.p_value,'RNA_q':e.q_value,'RNA_status':e.status,'RNA_direction':'UP' if e.effect>0 else 'DOWN' if e.effect<0 else 'TIE' if e.effect==0 else 'NOT_EVALUABLE'};rows.append(row)
 pre=pd.DataFrame(rows);assert len(pre)==722 and pre.metabolite_name.nunique()==51 and pre.relation_id.is_unique
 # Both invocations construct exactly the same pre-scRNA table; reject unintended drift.
 path=out/'candidate_pre_scRNA.tsv'
 if path.exists():assert path.read_text(encoding='utf-8')==pre.to_csv(sep='\t',index=False,na_rep='NA'), 'Pre-scRNA content changed'
 else:write(path,pre)
 if args.append_sc:
  sc=pd.read_csv(SC/'cross_cohort_source.tsv',sep='\t');assert sc.gene.is_unique and set(sc.gene)==set(planned.gene);merged=pre.merge(sc,on='gene',how='left',validate='many_to_one');grades=[];reasons=[]
  for _,r in merged.iterrows():
   reason=[]
   if r.mapping_tier!='DIRECT':reason.append('MAPPING_CONDITIONAL_OR_UNRESOLVED')
   if r.get('CAMP_status')!='DONE':reason.append('ASSOCIATION_NOT_EVALUABLE')
   if r.get('CAMP_available_status')!='DONE':reason.append('AVAILABLE_SUBSET_NOT_EVALUABLE')
   if r.relation_label=='R-sensitive':reason.append('ASSOCIATION_SENSITIVE')
   if r.get('RNA_status')!='DONE':reason.append('RNA_BACKGROUND_NOT_EVALUABLE')
   if r.metabolite_available_status!='DONE':reason.append('PAIRED_METABOLITE_AVAILABLE_SUBSET_NOT_EVALUABLE')
   if str(r.metabolite_available_direction_agrees)=='False':reason.append('PAIRED_METABOLITE_AVAILABLE_DIRECTION_CHANGED')
   for c in ['GSE263733','GSE278688','GSE242230']:
    if pd.isna(r.get(c+'_detection_fraction')) or r[c+'_detection_fraction']<.05:reason.append(c+'_LOW_OR_UNAVAILABLE_DETECTION')
   if reason:grade='D'
   elif r.relation_label=='R-strong' and r.stable_source==True:grade='A' if r.consensus_top=='Epithelial/Ductal' else 'B';reason=['EXPLORATORY_INTERNAL_PRIORITY;NOT_EXTERNAL_VALIDATION']
   else:grade='C';reason=['BIOCHEMICAL_CANDIDATE_RETAINED',r.relation_label,'STABLE_SOURCE' if r.stable_source==True else 'SOURCE_UNSTABLE_OR_COHORT_DEPENDENT']
   grades.append(grade);reasons.append(';'.join(reason))
  merged['candidate_class']=grades;merged['class_reason']=reasons;merged['exact_relation_external_validation']='NOT_PERFORMED';merged['functional_perturbation']='NOT_PERFORMED';write(out/'candidate_scRNA_appended.tsv',merged)
  summary={'metabolites':51,'relations':715,'unresolved_feature_rows':7,'genes':538,'rows':722,'relation_labels':pre.relation_label.value_counts().to_dict(),'candidate_classes':merged.candidate_class.value_counts().to_dict(),'candidate_class_scope':'Relation rows plus7 unresolved placeholders;not gene counts','all_metabolites_nominal_exploratory':True}
  (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
  (out/'validation.json').write_text(json.dumps({'status':'PASS','all51_metabolites_retained':True,'all715_relations_retained':True,'all538_genes_appended':True,'pre_scRNA_unchanged_after_append':True,'no_significance_gate':True,'external_validation':'NOT_PERFORMED'},indent=2)+'\n')
 else:print('pre_scRNA ready',len(pre),'rows;all538 genes',flush=True)
 sourcepaths=[MAPPING/'planned_relations.tsv',MAPPING/'feature_coverage_51.tsv',PAIRED/'results.tsv',DIRECTION/'direction_counts.tsv',INTERNAL/'associations.tsv',INTERNAL/'paired_RNA.tsv',Path(__file__)]
 if args.append_sc:sourcepaths.append(SC/'cross_cohort_source.tsv')
 write(out/'source_manifest.tsv',pd.DataFrame([{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sourcepaths]));(out/'analysis_spec.json').write_text(json.dumps({'protocol':'docs/PDAC/PAIRED51_INTERNAL_V2_LOCK.md','version':'PDAC_paired51_internal_sc_v2','scope':'Full relation pool;descriptive tiers;not causal score','paired_metabolite_D_rule':'Available paired subset not evaluable or mean direction changed => D','RNA_significance_filter':False,'rows':722},indent=2)+'\n')
if __name__=='__main__':main()
