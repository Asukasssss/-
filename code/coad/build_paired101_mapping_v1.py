"""Build scoped direct/conditional relation pools from paired101 source evidence."""
import argparse,csv,hashlib,json,re,subprocess
from collections import Counter,defaultdict
from pathlib import Path
from datetime import datetime,timezone
from review_paired101_mapping_v1 import ROOT,WORK,read

COMPLEX={'SDHA','SDHB','SDHC','SDHD','SUCLG1','SUCLG2','SUCLA2','FARSA','FARSB','SPTLC1','SPTLC2','SPTLC3','SPTSSA','SPTSSB','GCLM','SLC3A1','SLC3A2'}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def table(p,rows,fields=None):
    fields=fields or list(rows[0])
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:w.writerow({k:json.dumps(r[k],ensure_ascii=False) if isinstance(r.get(k),(list,dict)) else r.get(k,'NA') for k in fields})


def modification_context(e):
    t=e['reaction'].lower()
    if re.search(r'protein|histone|collagen|calmodulin|peptide chain release factor|translation elongation factor|hypoxia-inducible factor|peptidyl|wyosine|wybutosine',t):return True
    # Retain ordinary aminoacyl-tRNA charging; separate nucleotide-residue modification.
    if 'rna' in t and not ('-trna(' in t and 'diphosphate' in t):return True
    return ' in dna' in t


def gtp_cycle_context(e):
    symbol=e['gene_symbol'].split(';')[0]
    return e['feature_name']=="guanosine 5'- diphosphate (GDP)" and bool(re.match(r'^(RAB\d|ARF\d|ARL\d|GNA|RHO|RAS|RAN$|SAR\d|CRACR|LSG1$)',symbol))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out
    out.mkdir(parents=True,exist_ok=False)
    ids=json.loads((WORK/'identity_review.json').read_text());by_name={r['feature_name']:r for r in ids}
    ev=json.loads((WORK/'human_evidence.json').read_text());chains=json.loads((WORK/'reaction_chains.json').read_text())
    db=json.loads((WORK/'evidence_extended.json').read_text())
    pairs=read(ROOT/'results/COAD/03_PATIENT/20260921T123200Z_paired_metabolites33_v1/primary33.tsv')
    pmap={r['metabolite_name']:r for r in pairs};assert len(by_name)==101
    old=json.loads((ROOT/'reports/COAD/current_catalog_v0_1/candidate_catalog.json').read_text())
    key=lambda r:(r['feature_name'],r['metabolite_key'],r['human_gene_id'])
    oldmap={key(r):r for r in old}
    legacy=read(ROOT/'results/COAD/03_PATIENT/20260919T122034Z_patient_v1/results.tsv')
    legacy={(r['metabolite_name'],r['metabolite_key'],r['human_gene_id']):r for r in legacy}
    rna={r['human_gene_id']:r for r in read(ROOT/'results/COAD/03_PATIENT/20260919T122034Z_patient_v1/paired_RNA33.tsv')}
    old35={r['gene'] for r in read(ROOT/'results/COAD/06_EXTERNAL/20260920T091100Z_cell_paired35_v1/relations39_integrated.tsv')}
    old39={(r['metabolite_name'],r['gene']) for r in read(ROOT/'results/COAD/06_EXTERNAL/20260920T091100Z_cell_paired35_v1/relations39_integrated.tsv')}
    grouped=defaultdict(list);chainmap=defaultdict(list)
    for r in chains:chainmap[key(r)].append(r)
    for e in ev:grouped[key(e)].append(e)
    allkeys=set(grouped)|set(chainmap)|{k for k in oldmap if k[0] in by_name}
    prefix=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t')
    result=[];support=[]
    for k in sorted(allkeys):
        feature,mkey,gid=k;i=by_name[feature];prs=pmap[feature];previous=oldmap.get(k);events=grouped[k]
        # Deduplicate evidence irrespective of whether it was retrieved by KEGG or direct Rhea.
        unique={}
        for e in events:
            ek=(e['uniprot_accession'],e['compound_hypothesis'],e['reaction'])
            unique.setdefault(ek,e)
        events=list(unique.values())
        exact=[e for e in events if e['target_chebi_matches']]
        eligible=[e for e in exact if not modification_context(e) and not gtp_cycle_context(e)]
        notes=[]
        if i['identity_review_status'].startswith('HOLD_'):pool='CONDITIONAL';notes.append(i['identity_review_status'])
        elif eligible:pool='DIRECT'
        elif exact:
            pool='CONDITIONAL'
            if any(modification_context(e) for e in exact):notes.append('MACROMOLECULE_MODIFICATION_CONTEXT')
            if any(gtp_cycle_context(e) for e in exact):notes.append('SIGNALING_OR_GTPASE_CYCLE_CONTEXT')
        elif events:pool='CONDITIONAL';notes.append('REACTION_MATCH_WITHOUT_EXACT_TARGET_CHEBI')
        else:pool='UNRESOLVED';notes.append('NO_EXACT_REVIEWED_HUMAN_REACTION_IN_SCOPED_SEARCH')
        desc=db['human_gene_names'].get(gid,'');symbol=desc.split(';')[0].split(',')[0].strip() if ';' in desc else (events[0]['gene_symbol'] if events else previous['gene_symbol'] if previous else chainmap[k][0]['gene_symbol'])
        types=[]
        if any(e['kind']=='TRANSPORT_REACTION' for e in eligible):types.append('DIRECT_TRANSPORT')
        if any(e['kind']!='TRANSPORT_REACTION' for e in eligible):types.append('DIRECT_SUBSTRATE_OR_PRODUCT')
        if symbol in COMPLEX:types.append('COMPLEX_MEMBER_NOT_INDEPENDENT_ENZYME_OR_CARRIER')
        if any(modification_context(e) for e in exact):types.append('HAS_SEPARATE_MACROMOLECULE_MODIFICATION_CONTEXT')
        oldstat=legacy.get(k,{});oldrna=rna.get(gid,{})
        x=dict.fromkeys(prefix,'NA');x.update(cancer='COAD',cohort='CAMP_COAD',stage_id='02_MAPPING',run_id=out.name,
            analysis_version='paired101_direct_mapping_v1',analysis_type='scoped_direct_annotation_mapping',metabolite_key=mkey,metabolite_name=feature,gene=symbol,
            unit='metabolite_feature_x_human_gene',status='DONE' if pool=='DIRECT' else 'NEEDS_REVIEW',reason=';'.join(notes) or 'EXACT_HUMAN_REACTION_TARGET_ANNOTATION',source_id='KEGG;reviewed_human_UniProt;Rhea')
        x.update(human_gene_id=gid,pool=pool,relation_types=';'.join(types) or 'NA',identity_review_status=i['identity_review_status'],identity_note=i['identity_note'],
            compound_hypotheses=i['compound_hypotheses'],source_urls=sorted({'https://www.uniprot.org/uniprotkb/'+e['uniprot_accession']+'/entry' for e in events}|{'https://www.kegg.jp/entry/'+r['reaction_id'] for r in chainmap[k]}),
            exact_human_support_count=len(exact),eligible_direct_support_count=len(eligible),experimental_annotation=any(z.get('evidenceCode')=='ECO:0000269' for e in eligible for z in e['annotation_evidence']),
            original_camp_effect=prs['original_effect'],original_camp_q=prs['original_q'],paired_metabolite_P=prs['p_value'],paired_metabolite_q=prs['q_value'],
            paired_up=prs['up_pairs'],paired_down=prs['down_pairs'],paired_equal=prs['equal_pairs'],both_observed_pairs=prs['both_observed_pairs'],
            previous_catalog_relation=k in oldmap,previous_catalog_disposition=previous['disposition'] if previous else 'NA',
            legacy_patient_status=oldstat.get('status','NOT_RUN'),legacy_patient_n=oldstat.get('n','NA'),legacy_patient_rho=oldstat.get('effect','NA'),legacy_patient_P=oldstat.get('p_value','NA'),legacy_patient_q=oldstat.get('q_value','NA'),
            current_patient_status='NOT_RUN',current_patient_P='NA',current_patient_q='NA',reuse_requirement='CHECK_INPUT_HASHES_BEFORE_REUSE;NEW_FAMILY_Q_REQUIRED',
            legacy_RNA_status=oldrna.get('status','NOT_RUN'),legacy_RNA_effect=oldrna.get('effect','NA'),legacy_RNA_P=oldrna.get('p_value','NA'),legacy_RNA_q=oldrna.get('q_value','NA'),
            prior_35gene_scRNA_panel=symbol in old35,prior_39relation=(feature,symbol) in old39,coad_functional_proof='NOT_ESTABLISHED_BY_MAPPING')
        result.append(x)
        for e in events:
            support.append(dict(metabolite_name=feature,metabolite_key=mkey,human_gene_id=gid,gene=symbol,compound_hypothesis=e['compound_hypothesis'],
                uniprot_accession=e['uniprot_accession'],retrieval_kind=e['kind'],reaction=e['reaction'],rhea_ids=e['rhea_ids'],target_chebi_matches=e['target_chebi_matches'],
                modification_context=modification_context(e),gtp_cycle_context=gtp_cycle_context(e),annotation_evidence=e['annotation_evidence'],subunit_notes=e['subunit_notes'],source_url='https://www.uniprot.org/uniprotkb/'+e['uniprot_accession']+'/entry'))
    assert len(result)==len({(r['metabolite_name'],r['metabolite_key'],r['human_gene_id']) for r in result})
    direct=[r for r in result if r['pool']=='DIRECT']
    feature_rows=[]
    for name,p in pmap.items():
        rr=[r for r in result if r['metabolite_name']==name];i=by_name.get(name)
        counts=Counter(r['pool'] for r in rr)
        state='OUTSIDE_CURRENT_P_WORK_POOL' if i is None else 'HAS_DIRECT_RELATIONS' if counts['DIRECT'] else 'CONDITIONAL_IDENTITY' if i['identity_review_status'].startswith('HOLD_') else 'CONDITIONAL_RELATIONS_ONLY' if counts['CONDITIONAL'] else 'NO_DIRECT_HIT_IN_SCOPED_SEARCH'
        feature_rows.append(dict(metabolite_name=name,metabolite_key=p['metabolite_key'],selected_paired_P_lt_005=i is not None,
            paired_P=p['p_value'],paired_q=p['q_value'],up_pairs=p['up_pairs'],down_pairs=p['down_pairs'],equal_pairs=p['equal_pairs'],both_observed_pairs=p['both_observed_pairs'],
            original_effect=p['original_effect'],original_q=p['original_q'],mapping_status=state,direct_relations=counts['DIRECT'],conditional_relations=counts['CONDITIONAL'],unresolved_relations=counts['UNRESOLVED'],
            identity_status=i['identity_review_status'] if i else 'NOT_REVIEWED_THIS_POOL',identity_note=i['identity_note'] if i else 'Outside paired-P work pool; not deleted',
            legacy_feature=i['identity_evidence_origin'].startswith('REUSED') if i else name in {r['feature_name'] for r in old},
            genes=';'.join(sorted({r['gene'] for r in rr if r['pool']=='DIRECT'})) or 'NA'))
    assert len(feature_rows)==159 and sum(r['selected_paired_P_lt_005'] for r in feature_rows)==101
    table(out/'results.tsv',result);table(out/'candidate_pre_scRNA.tsv',result)
    table(out/'feature_ledger159.tsv',feature_rows);table(out/'feature_workpool101.tsv',[r for r in feature_rows if r['selected_paired_P_lt_005']])
    table(out/'main_relations.tsv',direct);table(out/'supporting_evidence.tsv',support)
    genes=[]
    for gid in sorted({r['human_gene_id'] for r in direct}):
        rr=[r for r in direct if r['human_gene_id']==gid]
        genes.append(dict(human_gene_id=gid,gene=rr[0]['gene'],direct_relations=len(rr),metabolites=';'.join(r['metabolite_name'] for r in rr),
            legacy_RNA_status=rr[0]['legacy_RNA_status'],prior_35gene_scRNA_panel=rr[0]['prior_35gene_scRNA_panel'],
            current_RNA_analysis='NOT_RUN',current_scRNA_analysis='REUSE_OUTCOME_CHECK_REQUIRED' if rr[0]['prior_35gene_scRNA_panel'] else 'NOT_RUN'))
    table(out/'gene_coverage526.tsv',genes)
    table(out/'legacy_outside_current_pool.tsv',[r for r in old if r['feature_name'] not in by_name])
    old39status=[]
    for feature,gene in sorted(old39):
        found=[r for r in result if r['metabolite_name']==feature and r['gene']==gene]
        old39status.append(dict(metabolite_name=feature,gene=gene,new_pool_status=found[0]['pool'] if found else 'OUTSIDE_PAIRED_WORK_POOL',reason=found[0]['reason'] if found else 'METABOLITE_NOT_PAIRED_P_LT_005',old_result_preserved=True))
    table(out/'legacy39_transition.tsv',old39status)
    summary=dict(features_all=159,paired_workpool=101,features_reused=72,features_new=29,pool_counts=dict(Counter(r['pool'] for r in result)),
        total_relations=len(result),direct_genes=len({r['human_gene_id'] for r in direct}),direct_features=len({r['metabolite_name'] for r in direct}),
        feature_status_counts=dict(Counter(r['mapping_status'] for r in feature_rows if r['selected_paired_P_lt_005'])),
        direct_relations_reused=sum(r['previous_catalog_relation'] for r in direct),direct_relations_new=sum(not r['previous_catalog_relation'] for r in direct),
        direct_relations_with_legacy_patient_results=sum(r['legacy_patient_status']=='DONE' for r in direct),
        direct_genes_with_prior_scRNA_panel=len({r['human_gene_id'] for r in direct if r['prior_35gene_scRNA_panel']}),
        legacy39_status_counts=dict(Counter(r['new_pool_status'] for r in old39status)),new_patient_tests=False)
    dump(out/'summary.json',summary)
    dump(out/'analysis_spec.json',dict(version='paired101_direct_mapping_v1',created_utc=datetime.now(timezone.utc).isoformat(),base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        selection='Existing paired P<0.05; all159 retained',key='CAMP_COAD x actual feature x human GeneID',scope='Archived human reviewed UniProt enzyme/transport records plus81 additional KEGG-chain gene queries and targeted missing-path queries',
        direct_rule='No known identity hold; reviewed human GeneID; explicit target ChEBI or official pH form in catalytic/transport reaction; no macromolecule-modification-only or signaling/GTPase-cycle-only support',
        condition_rule='Known identity conflict/combined signal; modification-only context; reaction crosslink without exact target chemical',
        complex_rule='Curated necessary enzyme/carrier complex members are labelled, not independent enzymes/carriers; annotation support is not COAD functional proof',
        family_locked=True,planned_patient_relations=len(direct),planned_RNA_genes=len(genes),main_relation_file_sha256=sha(out/'main_relations.tsv'),
        next_BH_rule='Primary and available-value sensitivity each fixed891 planned relationships; pairedRNA fixed526 genes separately. Nonestimable P=1 internally, public P/q=NA.',
        limitations=['Scoped annotation search, not exhaustive biochemical literature','No primary measurement chemical authentication','No patient statistics recomputed','Legacy association/RNA q applies only to legacy families','Same metabolite used in multiple reaction contexts is not a specific causal target'],
        next_analysis='All direct relations; input-hash checked reuse of matching old statistics; recalculate BH for newly locked family; RNA and observed-value sensitivity separate'))
    inputs=[WORK/'identity_review.json',WORK/'human_evidence.json',WORK/'reaction_chains.json',WORK/'evidence_extended.json',ROOT/'reference/camp/cancer_effects.tsv',ROOT/'reports/COAD/current_catalog_v0_1/candidate_catalog.json',ROOT/'results/COAD/03_PATIENT/20260921T123200Z_paired_metabolites33_v1/primary33.tsv',Path(__file__)]
    table(out/'source_manifest.tsv',[dict(source_path=str(p.relative_to(ROOT) if p.is_absolute() else p),version='read_only_input_or_code',sha256=sha(p)) for p in inputs])
    # Preserve provenance of reused and newly retrieved database snapshots.
    registry={}
    for p in [ROOT/'runtime/coad_fresh_kegg_20260919/sources.json',ROOT/'runtime/coad_uniprot_rhea_20260919/sources.json',WORK/'kegg/sources.json',WORK/'uniprot/sources.json']:
        for url,item in json.loads(p.read_text()).items():registry[url]=dict(item,registry_path=str(p.relative_to(ROOT)))
    dump(out/'database_source_registry.json',registry)
    dump(out/'identity_review.json',ids)
    dump(out/'validation.json',dict(status='PASS',all159_retained=True,selected101=True,unique_relation_keys=True,
        no_identity_hold_in_direct=all(not r['identity_review_status'].startswith('HOLD_') for r in direct),
        direct_requires_exact_target_support=all(r['eligible_direct_support_count']>0 for r in direct),legacy_not_overwritten=True,
        new_patient_stats_not_run=True,output_sha256={p.name:sha(p) for p in out.glob('*.tsv')}))
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
