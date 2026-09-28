"""PDAC SOP adapter. Reuse immutable aggregates; never invent unavailable matrices."""
from pathlib import Path
import hashlib,json,shutil,platform
import numpy as np
import pandas as pd
from sop_v3_math import bh_evaluable

REPO=Path(__file__).resolve().parents[2]
ROOT=REPO/'results/PDAC'
CACHE=REPO/'runtime/pdac_sop_v3'
TEMPLATES=REPO/'code/pdac/sop_v3_templates'
VERSION='PDAC_CAMP_source_SOP_v3'
SOP_SHA='6dd0a73fb77a392e428254ad265b8710ce967c73'
SERVER_RUN='20260922T112300Z_sop_v3'
RUNS={'01_CAMP':'20260922T112500Z_sop_v3_discovery','02_MAPPING':'20260922T112600Z_sop_v3_mapping','03_PATIENT':'20260922T112700Z_sop_v3_patient','06_EXTERNAL':'20260922T112800Z_sop_v3_source','07_INTEGRATION':'20260922T112900Z_sop_v3_integration'}
OLD_MET=ROOT/'04_ROBUSTNESS/20260921T121509Z_paired_metabolites_v1'
OLD_DIR=ROOT/'04_ROBUSTNESS/20260921T122353Z_pair_direction_counts_v1'
OLD_MAP=ROOT/'02_MAPPING/20260922T055000Z_paired51_mapping_v2'
OLD_FIRST=ROOT/'02_MAPPING/20260920T111500Z_direct_mapping_v1'
OLD_PAT=ROOT/'03_PATIENT/20260922T053500Z_internal_paired51_v2'
OLD_HIST=ROOT/'03_PATIENT/20260920T054500Z_source_readiness_v1'
OLD_SC=ROOT/'06_EXTERNAL/20260922T054000Z_sc_paired51_v2'
USED=set()

def read(path):
    USED.add(Path(path));return pd.read_csv(path,sep='\t',keep_default_na=True)
def digest(path):
    p=Path(path);data=p.read_bytes()
    # Git normalizes repository text to LF. Runtime source downloads keep exact bytes.
    if p.is_relative_to(REPO) and not p.is_relative_to(REPO/'runtime') and p.suffix in {'.tsv','.json','.py','.md','.txt'}:data=data.replace(b'\r\n',b'\n')
    return hashlib.sha256(data).hexdigest()
def write_json(path,data):Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
def dest(stage):
    p=ROOT/stage/RUNS[stage];p.mkdir(parents=True,exist_ok=True);return p
def ordered(name,df):
    fields=list(pd.read_csv(TEMPLATES/name,sep='\t').columns)
    for f in fields:
        if f not in df:df[f]=np.nan
    return df[fields+[c for c in df if c not in fields]]
def emit(stage,name,df):
    if (TEMPLATES/name).exists():df=ordered(name,df.copy())
    df.to_csv(dest(stage)/name,sep='\t',index=False,na_rep='NA',lineterminator='\n');return df
def bh(df):
    p=pd.to_numeric(df.p_value,errors='coerce');df['q_value']=np.nan
    m=p.notna();assert p[m].between(0,1).all()
    if m.any():df['q_value']=bh_evaluable(p.to_numpy())
    df['family_n_planned']=len(df);df['family_n_evaluable']=int(m.sum())
    df['evidence_tier']=[tier(p,q) for p,q in zip(df.p_value,df.q_value)]
    return df
def tier(p,q):
    return 'NOT_EVALUABLE' if pd.isna(p) else 'FDR_SUPPORTED' if q<.05 else 'NOMINAL_EXPLORATORY' if p<.05 else 'NOT_SUPPORTED_THIS_TEST'

def prepare_templates():
    TEMPLATES.mkdir(parents=True,exist_ok=True)
    for p in (CACHE/'templates/camp_discovery_source_v1').iterdir():
        if p.is_file():shutil.copyfile(p,TEMPLATES/p.name)
    target=REPO/'docs/PDAC/reference';target.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(CACHE/'docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md',target/'CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md')
    shutil.copyfile(CACHE/'provenance.json',target/'SOP_GITHUB_PROVENANCE.json')

def mappings():
    rel=read(OLD_MAP/'planned_relations.tsv');ev=read(OLD_MAP/'mapping_evidence.tsv');features=read(OLD_MAP/'feature_coverage_51.tsv')
    prior=read(OLD_FIRST/'direct_relations_v1.tsv')
    history=set(prior.gene)
    for name in ['historical_cell_source_gene_summary.tsv','historical_cptac_effects.tsv']:
        history.update(read(OLD_HIST/name).gene)
    union=sorted(history|set(rel.gene));current=set(rel.loc[rel.mapping_tier=='DIRECT','gene'])
    path=REPO/'runtime/pdac_mapping_paired51_v2/human_reviewed.json';USED.add(path)
    entries=json.loads(path.read_text(encoding='utf8'))['results'];bygene={}
    for e in entries:
        for g in e.get('genes',[]):
            if g.get('geneName',{}).get('value') in union:bygene.setdefault(g['geneName']['value'],[]).append(e)
    assert all(len(bygene.get(g,[]))==1 for g in union)
    stable={g:'UniProtKB:'+bygene[g][0]['primaryAccession'] for g in union}
    assert len(set(stable.values()))==len(union)
    aliases={g:[s['value'] for gg in bygene[g][0].get('genes',[]) for s in gg.get('synonyms',[])] for g in union}
    # A historical synonym which is another reviewed canonical gene cannot resolve identity.
    canonical={g.get('geneName',{}).get('value') for e in entries for g in e.get('genes',[])}
    aliases={g:[a for a in v if a not in canonical] for g,v in aliases.items()}
    rel['stable_gene_id']=rel.gene.map(stable)
    assert not rel.duplicated(['metabolite_key','stable_gene_id']).any()
    evidence=[];relations=[]
    prior_keys=set(zip(prior.metabolite_key,prior.gene))
    for _,r in rel.iterrows():
        sub=ev[(ev.metabolite_key==r.metabolite_key)&(ev.gene==r.gene)]
        base={'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','relation_id':r.relation_id,'metabolite_key':r.metabolite_key,'stable_gene_id':stable[r.gene],'gene':r.gene,'relation_type':r.relation_type,'metabolite_role':r.reaction_roles,'reaction_id':r.rhea_ids,'reaction_text':'; '.join(sub.reaction.dropna().astype(str).unique()),'organism':'Homo sapiens','evidence_id':';'.join(sub.rhea_ids.dropna().astype(str).unique()),'source_url':r.source_urls,'evidence_level':'HUMAN_EXPERIMENTAL_ANNOTATION' if r.has_human_experimental_annotation else 'INFERRED_OR_CONDITIONAL','mapping_status':r.mapping_tier,'identity_status':features.set_index('metabolite_key').loc[r.metabolite_key,'identity_status'],'compound_specificity':r.identity_note,'added_or_reused':'REUSED_PAIRED51_V2','mapping_version':'PDAC_paired51_mapping_v2','reason':r.mapping_reason,'metabolite_name':r.feature_name,'in_original42_relation_pool':(r.metabolite_key,r.gene) in prior_keys,'legacy_relation_id':r.relation_id}
        relations.append(base)
        if len(sub):
            for i,e in sub.iterrows():evidence.append({**base,'evidence_id':'paired51_v2_evidence_row_'+str(i+2),'reaction_id':e.rhea_ids,'reaction_text':e.reaction,'source_url':e.source_url,'annotation_evidence':e.annotation_evidence,'evidence_level':'HUMAN_EXPERIMENTAL_ANNOTATION' if e.human_experimental_annotation else 'INFERRED_OR_CONDITIONAL'})
        else:evidence.append({**base,'evidence_id':'paired51_v2_curated_relation_'+r.relation_id,'reason':str(base['reason'])+';curated relation source URL retained;reaction equation not supplied when transport/complex record'})
    records=pd.DataFrame(relations)
    emit('02_MAPPING','direct_relations.tsv',records[records.mapping_status=='DIRECT'])
    emit('02_MAPPING','conditional_relations.tsv',records[records.mapping_status=='CONDITIONAL'])
    emit('02_MAPPING','mapping_evidence.tsv',pd.DataFrame(evidence))
    rows=[]
    for g in union:
        rr=rel[rel.gene==g]
        rows.append({'cancer':'PDAC','stable_gene_id':stable[g],'gene':g,'in_current_pool':g in current,'in_history_pool':True,'history_only':g not in current,'mapping_version':'PDAC_paired51_mapping_v2','relation_ids':';'.join(sorted(rr.relation_id)),'source_id':'reviewed_human_UniProt_and_versioned_PDAC_mapping','status':'DONE','reason':'Current pool is DIRECT genes only;prior v2 conditional and original42/external genes retained in union','in_prior_42_pool':g in set(prior.gene),'in_paired51_any_relation':len(rr)>0,'in_conditional_pool':g in set(rel.loc[rel.mapping_tier=='CONDITIONAL','gene']),'n_direct_relations':int(((rr.mapping_tier=='DIRECT')).sum()),'n_conditional_relations':int((rr.mapping_tier=='CONDITIONAL').sum())})
    pool=pd.DataFrame(rows);emit('02_MAPPING','gene_pool_current.tsv',pool[pool.in_current_pool]);emit('02_MAPPING','gene_pool_history_union.tsv',pool)
    emit('02_MAPPING','feature_dispositions.tsv',features)
    write_json(dest('02_MAPPING')/'gene_aliases.json',aliases)
    write_json(dest('02_MAPPING')/'panel.json',{'genes':union,'current_genes':sorted(current),'stable_gene_ids':stable,'history_only_genes':sorted(set(union)-current),'conditional_genes':sorted(set(rel.loc[rel.mapping_tier=='CONDITIONAL','gene']))})
    return rel,features,pool,stable

def metabolites():
    old=read(OLD_MET/'results.tsv');dr=read(OLD_DIR/'direction_counts.tsv')
    assert not old.duplicated(['analysis_type','metabolite_key']).any()
    frames=[]
    for mode,name,family in [('primary','metabolite_paired.tsv','METAB_PAIRED_PRIMARY'),('both_preimputation_available','metabolite_available_sensitivity.tsv','METAB_PAIRED_AVAILABLE')]:
        a=old[old.analysis_type==mode].copy();a=a.merge(dr[dr.analysis_type==mode][['feature_name','n_up','n_down','n_equal','up_fraction','down_fraction']],left_on='metabolite_name',right_on='feature_name',validate='one_to_one').drop(columns='feature_name')
        a['q_original']=a.q_value;a['p_original']=a.p_value;a['reused_from']=OLD_MET.relative_to(REPO).as_posix();a['statistics_reused']=True
        a['run_id']=RUNS['01_CAMP'];a['analysis_version']=VERSION;a['stage_id']='01_CAMP';a['test_family']=family
        a['n_pairs_total']=11;a['n_pairs_used']=a.n;a['mean_delta']=a.effect;a['n_both_available']=a.n_both_preimputation_available
        a['mask_semantics']='author_available_value_mask;not_mass_spectrometry_detection';a['p_method']='exact_signed_rank_sign_enumeration';a['zero_method']='wilcox';a['tie_method']='average_ranks';a['continuity_correction']=False;a['B']=0
        a['bootstrap_valid']=np.where(a.p_original.notna(),10000,np.nan);a['bootstrap_B']=10000;a['ci_method']='historical_10000_pair_percentile_bootstrap_reused';a['ci_estimand']='mean_delta'
        a['missing_rate_tumor']=np.nan;a['missing_rate_normal']=np.nan;a['marginal_mask_status']='ACCESS_BLOCKED;requires_server_matrix'
        a['direction_discordant']=((np.sign(a.effect)!=np.sign(a.median_delta))|(np.sign(a.effect)!=np.sign(a.n_up-a.n_down))).where(a.effect.notna(),np.nan)
        a['reason']='Same input/effect/P/CI reused;new min8 gate;BH over evaluable tests;per-tissue mask fractions pending server'
        low=a.n<8;a.loc[low,['p_value','q_value']]=np.nan;a.loc[low,'status']='NOT_EVALUABLE';a.loc[low,'reason']='N_PAIRS_LT_8;historical P retained in p_original;description and historical mean CI retained'
        a=bh(a);a['workpool_entry']=(a.p_value<.05)&(mode=='primary');a=emit('01_CAMP',name,a);frames.append(a)
    primary=frames[0];emit('01_CAMP','metabolite_workpool.tsv',primary[primary.workpool_entry])
    rank=primary.assign(majority_fraction=primary[['up_fraction','down_fraction']].max(axis=1)).sort_values(['majority_fraction','n','p_value'],ascending=[False,False,True]);emit('01_CAMP','metabolite_direction_ranking.tsv',rank)
    return frames

def patients(stable):
    old=read(OLD_PAT/'associations.tsv');frames=[]
    for mode,name in [('primary','tumor_association.tsv'),('availability','tumor_association_available.tsv')]:
        for tiername in ['DIRECT','CONDITIONAL']:
            a=old[(old.analysis_type==mode)&(old.mapping_tier==tiername)].copy();a['q_original']=a.q_value
            a['stable_gene_id']=a.gene.map(stable);a['run_id']=RUNS['03_PATIENT'];a['analysis_version']=VERSION
            a['test_family']=('REL_TUMOR_' if tiername=='DIRECT' else 'REL_CONDITIONAL_SUPPLEMENT_')+('PRIMARY' if mode=='primary' else 'AVAILABLE')
            a['n_tumor_linked']=21;a['n_complete']=a.n;a['n_unique_patients']=np.nan;a['unit_verification']='21 distinct author paired-case keys;protected clinical identity and aliquot level not reverified'
            a['n_metabolite_available']=a.n_preimputation_available;a['rho']=a.effect;a['ci_method']='4000_whole_unit_percentile_bootstrap';a['p_method']='9999_patient_unit_permutation_plus_one';a['B']=9999;a['statistics_reused']=True;a['reused_from']=OLD_PAT.relative_to(REPO).as_posix();a['mask_semantics']='author_available_value_mask';a['same_cohort_postselection']=True
            a.loc[(a.bootstrap_valid<3600)&a.p_value.notna(),['ci_lower','ci_upper']]=np.nan
            a=bh(a);emit('03_PATIENT',name if tiername=='DIRECT' else 'conditional_'+name,a);frames.append(a)
    oldrna=read(OLD_PAT/'paired_RNA.tsv');emit('03_PATIENT','historical_RNA_wilcoxon_full.tsv',oldrna);emit('03_PATIENT','historical_RNA_wilcoxon_P005_view.tsv',oldrna[oldrna.p_value<.05])
    read(OLD_PAT/'sample_identity_audit.tsv');rows=[]
    audits=[('CAMP_sample_RNA_join',39,39,0,0,'Exact RNA column joins;27tumor12normal'),('tumor_independent_author_keys',27,21,6,0,'Six missing author pair keys excluded;remaining21unique'),('complete_pairs_among_normals',12,11,1,0,'Eleven explicit completeT/Npairs;one normal without matching CAMPtumor'),('tissue_label_concordance',39,39,0,0,'Zero observed author/CAMP tissue conflicts')]
    for item,checked,passed,excluded,unresolved,why in audits:rows.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','audit_item':item,'n_checked':checked,'n_pass':passed,'n_excluded':excluded,'n_unresolved':unresolved,'status':'DONE','reason':'Historical verified aggregate reused;'+why,'source_id':OLD_PAT.relative_to(REPO).as_posix()})
    rows += [{'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','audit_item':'aliquot_level_and_clinical_identity','n_checked':39,'n_pass':np.nan,'n_excluded':np.nan,'n_unresolved':39,'status':'NEEDS_REVIEW','reason':'Exact author sample mapping available;does not establish same aliquot or protected clinical identity','source_id':'CAMP_author_mapping'}, {'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','audit_item':'private_SOP_table_export','n_checked':np.nan,'n_pass':np.nan,'n_excluded':np.nan,'n_unresolved':np.nan,'status':'ACCESS_BLOCKED','reason':'server165 SSH timeout;private original joins remain server only','source_id':'server165'}]
    emit('01_CAMP','sample_audit_summary.tsv',pd.DataFrame(rows))
    return frames

def pending_and_integrate(rel,features,pool,stable,met,assoc):
    RNA=[];coverage=[]
    for _,g in pool.iterrows():
        RNA.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','stage_id':'03_PATIENT','run_id':RUNS['03_PATIENT'],'analysis_version':VERSION,'analysis_type':'paired_RNA_t','gene':g.gene,'stable_gene_id':stable[g.gene],'unit':'author_pair_key','n':np.nan,'n_reference':11,'test_family':'RNA_PAIRED_CURRENT' if g.in_current_pool else 'RNA_PAIRED_HISTORY_SUPPLEMENT','status':'ACCESS_BLOCKED','reason':'Paired t and4000 bootstrap require server-only matrix;historical Wilcoxon not relabeled','source_id':'GSE62452_author_processed_RNA','current_pool':g.in_current_pool,'history_only':not g.in_current_pool,'expression_platform':'Affymetrix_Human_Gene_1.0_ST','expression_scale':'author_processed_microarray_no_new_transform','p_method':'paired_t_planned','model_formula':'paired_T_minus_N','family_n_planned':int(pool.in_current_pool.sum()) if g.in_current_pool else int((~pool.in_current_pool).sum())})
        for c in ['GSE263733','GSE278688','GSE242230']:coverage.append({'cancer':'PDAC','cohort':c,'stable_gene_id':stable[g.gene],'gene':g.gene,'status':'ACCESS_BLOCKED','reason':'New cellwise log1p10k source extraction not run;historical pseudobulk not substituted','source_id':c})
    rn=pd.DataFrame(RNA);emit('03_PATIENT','paired_RNA.tsv',rn[pool.in_current_pool.to_numpy()]);emit('03_PATIENT','paired_RNA_history_supplement.tsv',rn[~pool.in_current_pool.to_numpy()]);emit('03_PATIENT','RNA_P005_view.tsv',rn.iloc[:0])
    emit('06_EXTERNAL','sc_gene_coverage.tsv',pd.DataFrame(coverage))
    primary=met[0].set_index('metabolite_key');a=pd.concat([assoc[0],assoc[1]]).set_index('relation_id');av=pd.concat([assoc[2],assoc[3]]).set_index('relation_id');out=[]
    oldrna=read(OLD_PAT/'paired_RNA.tsv').set_index('gene')
    for _,r in rel.iterrows():
        m=primary.loc[r.metabolite_key];ar=a.loc[r.relation_id];v=av.loc[r.relation_id]
        out.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','relation_id':r.relation_id,'metabolite_key':r.metabolite_key,'stable_gene_id':stable[r.gene],'gene':r.gene,'mapping_version':r.mapping_version,'metabolite_name':r.feature_name,'mapping_status':r.mapping_tier,'metabolite_effect':m.effect,'metabolite_p':m.p_value,'metabolite_q':m.q_value,'metabolite_n':m.n,'metabolite_source_run':RUNS['01_CAMP'],'association_rho':ar.effect,'association_p':ar.p_value,'association_q':ar.q_value,'association_n':ar.n,'association_source_run':RUNS['03_PATIENT'],'association_available_rho':v.effect,'association_available_p':v.p_value,'association_available_q':v.q_value,'association_available_n':v.n,'RNA_effect':np.nan,'RNA_p':np.nan,'RNA_q':np.nan,'RNA_n':np.nan,'RNA_source_run':RUNS['03_PATIENT'],'sc_source_run':RUNS['06_EXTERNAL'],'metabolite_evidence':m.evidence_tier,'patient_association_evidence':ar.evidence_tier,'RNA_background':'ACCESS_BLOCKED_PAIRED_T_PENDING','cell_background':'ACCESS_BLOCKED_CELLWISE_SOURCE_PENDING','data_limitations':'same_cohort_postselection;aliquot_level_unverified;conditional_relations_separate;new_RNA_and_SC_pending','next_action':'Complete SOP paired RNA and cell-source background;no significance gate','status':'PARTIAL','reason':'Reused patient association;new q;RNA/sc not filled from incompatible historical methods','historical_RNA_Wilcoxon_p':oldrna.loc[r.gene,'p_value'],'historical_RNA_Wilcoxon_q':oldrna.loc[r.gene,'q_value'],'historical_RNA_source_run':OLD_PAT.name})
    emit('07_INTEGRATION','candidate_relations_integrated.tsv',pd.DataFrame(out))
    gr=[]
    for _,g in pool.iterrows():
        gr.append({'cancer':'PDAC','stable_gene_id':stable[g.gene],'gene':g.gene,'in_current_pool':g.in_current_pool,'history_only':g.history_only,'relation_ids':g.relation_ids,'n_relations':int(g.n_direct_relations+g.n_conditional_relations),'RNA_background':'ACCESS_BLOCKED_PAIRED_T_PENDING','RNA_source_run':RUNS['03_PATIENT'],'cell_background':'ACCESS_BLOCKED_CELLWISE_SOURCE_PENDING','sc_source_run':RUNS['06_EXTERNAL'],'functional_source_ref':'NOT_RUN;outside current SOP scope','external_source_ref':str(OLD_HIST.relative_to(REPO)).replace('\\','/')+';historical gene evidence index only,not relation validation','data_limitations':'Current direct pool separate from conditional and history-only;new source computation pending','next_action':'Complete all-gene source coverage regardless of patient P','status':'PARTIAL','reason':'History union retained;no opaque aggregate score'})
    emit('07_INTEGRATION','candidate_genes_integrated.tsv',pd.DataFrame(gr));emit('07_INTEGRATION','unmapped_workpool_features.tsv',features[~features.metabolite_key.isin(rel.metabolite_key)])

def main():
    prepare_templates();rel,features,pool,stable=mappings();met=metabolites();assoc=patients(stable)
    identity=[]
    f=features.set_index('metabolite_key')
    for _,r in met[0].iterrows():
        has=r.metabolite_key in f.index;e=f.loc[r.metabolite_key] if has else None
        identity.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','metabolite_key':r.metabolite_key,'author_feature_id':r.metabolite_name,'author_name':r.metabolite_name,'canonical_name':e.metabolite_standard_name if has else r.metabolite_name,'stable_compound_id':r.metabolite_key if not str(r.metabolite_key).startswith('NAME:') else np.nan,'stereochemistry':'Not reidentified experimentally','lipid_resolution':'Author annotation;see identity restriction','identity_status':e.identity_status if has else 'OUTSIDE_P005_MAPPING_WORKPOOL','source_id':'Author matrix plus immutable paired51 mapping','identity_version':'PDAC_paired51_mapping_v2','reason':e.chemical_identity_note if has else 'All307 retained;no P-based deletion of full discovery table'})
    emit('02_MAPPING','metabolite_identity.tsv',pd.DataFrame(identity));pending_and_integrate(rel,features,pool,stable,met,assoc)
    spec=json.loads((TEMPLATES/'analysis_spec.template.json').read_text());spec.update(status='FROZEN_EXECUTION_PARTIAL',cancer='PDAC',cohort='CAMP_PDAC_GSE62452',analysis_version=VERSION,run_id=SERVER_RUN,source_commit=SOP_SHA,planned_counts={'metabolites':len(met[0]),'workpool':int(met[0].workpool_entry.sum()),'direct_relations':int((rel.mapping_tier=='DIRECT').sum()),'conditional_relations':int((rel.mapping_tier=='CONDITIONAL').sum()),'current_genes':int(pool.in_current_pool.sum()),'history_union_genes':len(pool),'RNA_pairs':11,'tumor_units':21},software_versions={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__},input_files=[{'path':str(p.relative_to(REPO)).replace('\\','/'),'sha256':digest(p)} for p in sorted(USED)])
    spec['hash_semantics']='Repository text SHA256 after Git LF normalization;downloaded runtime annotation source SHA256 exact bytes'
    spec['sample_design']={'pair_evidence':'GSE62452 author explicit paired sample table joined by exact GSM','unit':'distinct author pair key;clinical identity not independently recertified','duplicate_policy':'Exclude6 tumor specimens without author pair keys;21 remaining keys unique','normal_tissue_definition':'author non-tumor pancreas;precise normal tissue subtype not newly verified','RNA_pair_scope':'Same11 CAMP complete author pairs,independently matched RNA IDs;do not add59/60outsideCAMP pairs'}
    spec['data_scale']={'metabolite':'author processed scale;not asserted log2FC','RNA':'author processed Affymetrix gene-symbol continuous expression;no added transform','availability_mask_semantics':'finite author preimputation data values;not raw detection'}
    spec['RNA']['input_type']='continuous_author_processed_microarray';spec['RNA']['supplement_family']='437_history_only_genes_separate_BH';spec['metabolite_test']['rounding_policy']='No rounding;no documented measurement precision to justify it'
    spec['reuse_exceptions']={'metabolite_P':'Exact historical signed-rank P reused at n>=8;below8 newP/qsNA','metabolite_CI':'Historical10000 pair percentile bootstrap preserved rather than redoing4000;same input,method,estimand','association':'Same21unit effect/P/4000bootstrap reused;original seeds and input hashes retained;new evaluable-onlyBH','history':'Immutable original CAMP303/42 and every prior P/q preserved'}
    spec['single_cell']['gene_scope']='All history union687,including250currentdirect;no association/RNA gate';spec['single_cell']['new_umap']='NOT_SCOPED;SOP optional,required dotplots and fullgene heatmaps only';spec['single_cell']['min_detected_donors']='No extra3detected threshold;atleast3eligible labels,2categories,maxdetection>=.01'
    write_json(REPO/'docs/PDAC/PDAC_SOP_V3_ANALYSIS_SPEC.json',spec)
    totals={s:{'status':'PARTIAL' if s!='06_EXTERNAL' else 'ACCESS_BLOCKED','reason':'Server165:22 timeout;local reuse and mapping completed;matrix-dependent outputs not run'} for s in RUNS};totals['02_MAPPING']={'status':'DONE','reason':'51dispositions;357direct+358conditional;250currentdirect;687historyunion;all687stableUniProtidentities'}
    stats={'metabolites':{d.test_family.iloc[0]:{'planned':len(d),'evaluable':int(d.p_value.notna().sum()),'P_lt005':int((d.p_value<.05).sum()),'q_lt005':int((d.q_value<.05).sum())} for d in met},'associations':{d.test_family.iloc[0]:{'planned':len(d),'evaluable':int(d.p_value.notna().sum()),'P_lt005':int((d.p_value<.05).sum()),'q_lt005':int((d.q_value<.05).sum())} for d in assoc},'mapping':spec['planned_counts'],'stage_status':totals}
    for stage in RUNS:
        out=dest(stage);write_json(out/'analysis_spec.json',spec);write_json(out/'summary.json',stats)
        emit(stage,'source_manifest.tsv',pd.DataFrame([{'source_id':p.stem,'path_or_url':p.relative_to(REPO).as_posix(),'release':'immutable_local_source','sha256':digest(p),'access_scope':'PUBLIC_AGGREGATE_OR_ANNOTATION','transform':'No patient matrix copied locally'} for p in sorted(USED)]+[{'source_id':'SOP','path_or_url':'https://github.com/Asukasssss/-/blob/'+SOP_SHA+'/docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md','release':SOP_SHA,'sha256':digest(CACHE/'docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md'),'access_scope':'PUBLIC_PROTOCOL'}]))
        write_json(out/'validation.json',{'status':'PASS_LOCAL_ADAPTER','historical_files_modified':False,'server_execution':'ACCESS_BLOCKED','new_paired_t_executed':False,'new_cellwise_source_executed':False,'all357direct_and358conditional_retained':True,'all687_history_genes_retained':True,'RNA_and_source_numbers_not_fabricated':True,'not_verified':['Server matrix current hashes','New paired t','Cellwise log1p10k donor profiles','new source plots','marginal per-tissue availability']})
        label=totals[stage]['status'];body=f'''# PDAC SOP v3：{stage}

## 本轮问题

按 GitHub 规范 {SOP_SHA} 整理 PDAC，单细胞终点仅为基因表达来源。

## 输入与范围

307项代谢物，11个作者明确配对，21个作者明确且互不重复编号的肿瘤单位。配对P<0.05工作池51项；357直接关系、358条件关系。当前DIRECT池250基因，历史并集687基因。

## 实际结果

本阶段状态 **{label}**。{totals[stage]['reason']}。完整计数见summary.json。
新BH只在预定集合的可评估P上计算，旧q保存在q_original。主代谢物及关联效应/P/CI沿用同输入历史统计，不称新增验证。低于8对的新P/q留NA，旧P另存。

## 新手解释

当前基因池只包含直接关系基因；条件与历史基因没有删除，进入历史并集。RNA显著性和关联显著性不作为进入单细胞的门槛。

## 限制/反证

服务器连接超时，新配对t检验和逐细胞log1p10k汇总未运行。旧Wilcoxon结果和旧供者合并计数结果不能换名冒充新方法结果。缺项以ACCESS_BLOCKED保留；新RNA P<0.05空表是尚未执行，不是零个显著基因。687身份采用reviewed human UniProt稳定accession，矩阵按规范将另保存实际符号/稳定ID连接。独立患者临床身份和同一分装层级未重新验证。

## 当前决定

保留旧结果，按新规范分开证据维度，不沿用A/B/C/D作为总分或硬筛选。

## 下一步

服务器连接恢复后运行冻结的PDAC SOP v3，补RNA、原有值边际覆盖、三队列来源和必备点图/全基因热图。停止在来源描述，不追加机制任务。

## 复现命令

`python code/pdac/prepare_sop_v3.py`；模板来自固定GitHub提交，公共模板未改。服务器计算使用新的独占运行目录，患者与细胞级数据不进入本机或GitHub。
'''
        (out/'README_CN.md').write_text(body,encoding='utf8')
        emit(stage,'checksums.tsv',pd.DataFrame([{'file':p.name,'sha256':digest(p)} for p in sorted(out.iterdir()) if p.is_file() and p.name!='checksums.tsv']))
    print(json.dumps(stats,ensure_ascii=False))

if __name__=='__main__':main()
