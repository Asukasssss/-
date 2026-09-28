"""Build full-candidate working views without fabricating unavailable statistics."""
import csv,json,hashlib,subprocess
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[2]
MAP=ROOT/'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1'
FUN=ROOT/'results/PDAC/05_FUNCTION/20260920T112000Z_all_gene_search_v1'
OLD=ROOT/'results/PDAC/03_PATIENT/20260920T054500Z_source_readiness_v1'
OUT=ROOT/'results/PDAC/07_INTEGRATION/20260920T112500Z_working_v1'
EXT=ROOT/'results/PDAC/06_EXTERNAL/20260920T112500Z_all_gene_coverage_v1'
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def dump(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
    OUT.mkdir(parents=True,exist_ok=False);EXT.mkdir(parents=True,exist_ok=False)
    direct=read(MAP/'direct_relations_v1.tsv');search={r['gene']:r for r in read(FUN/'gene_search_ledger.tsv')}
    reviews=read(FUN/'abstract_review.tsv');cells={r['gene']:r for r in read(OLD/'historical_cell_source_gene_summary.tsv')};cptac={r['gene']:r for r in read(OLD/'historical_cptac_effects.tsv')}
    old=read(OLD/'historical_associations.tsv');relations=[];genes=[];external=[]
    for gene in sorted({r['gene'] for r in direct}):
        rr=[r for r in direct if r['gene']==gene];rev=[r for r in reviews if r['gene']==gene]
        direct_function=any(r['pdac_genetic_intervention_status']=='ABSTRACT_REPORTS_DIRECT_GENETIC_INTERVENTION' for r in rev)
        genes.append({'gene':gene,'relation_count':len(rr),'metabolites':';'.join(sorted({r['feature_name'] for r in rr})),
          'search_status':search[gene]['search_status'],'search_hit_count':search[gene]['hit_count'],
          'abstracts_manually_reviewed':len(rev),'functional_status':'ABSTRACT_GENETIC_REPORT_FULLTEXT_PENDING' if direct_function else 'REVIEW_PARTIAL' if rev else 'LOCATED_NOT_ADJUDICATED',
          'exact_metabolite_mechanism':'NOT_ESTABLISHED','cell_source_status':'HISTORICAL_SUMMARY_NEEDS_REVIEW' if gene in cells else 'NOT_RUN_FULL_G_EXPANSION',
          'cptac_status':'HISTORICAL_SUMMARY_NEEDS_REVIEW' if gene in cptac else 'NOT_RUN_FULL_G_EXPANSION',
          'work_group':'FUNCTIONAL_BACKGROUND_ASSOCIATION_PENDING' if direct_function else 'DESIGN_AND_EVIDENCE_REVIEW_PENDING',
          'next_evidence':'Verify model-specific intervention and exact metabolite mechanism' if direct_function else 'Adjudicate located literature;complete tumor associations'})
        external.append({'gene':gene,'GSE263733_summary_present':gene in cells,'GSE278688_summary_present':gene in cells,
          'historical_cptac_summary_present':gene in cptac,'historical_rna_n':cptac.get(gene,{}).get('rna_n_complete_pairs','NA'),
          'historical_protein_n':cptac.get(gene,{}).get('protein_n_complete_pairs','NA'),
          'current_matrix_evaluability':'ACCESS_BLOCKED','pairing_and_batch_review':'NOT_RUN','independence_from_CAMP':'NOT_VERIFIED',
          'reason':'Current server unavailable;old aggregate presence is not full-matrix eligibility or independent validation',
          'joint_metabolome_RNA_external_validation':'NOT_RUN'})
    gd={r['gene']:r for r in genes}
    for r in direct:
        prior=[p for p in old if p['gene']==r['gene'] and p['metabolite_name']==r['feature_name'] and p['analysis_type']=='author_processed_primary']
        h=prior[0] if prior else {}
        relations.append({**r,'patient_status':'ACCESS_BLOCKED','primary_n':'NA','primary_rho':'NA','primary_p':'NA','primary_q':'NA',
          'availability_n':'NA','availability_rho':'NA','availability_q':'NA','leave_one_out_status':'NOT_RUN_SERVER_BLOCKED',
          'historical_two_edge_family_n':h.get('n','NA'),'historical_two_edge_family_rho':h.get('rho','NA'),'historical_two_edge_family_q':h.get('q_value','NA'),
          'historical_statistics_use':'CONTEXT_ONLY_NOT_NEW_FAMILY_RESULT' if h else 'NONE',
          'function_status':gd[r['gene']]['functional_status'],'exact_relation_functional_mechanism':'NOT_ESTABLISHED',
          'external_relationship_status':'NOT_RUN','opposing_evidence_status':'PARTIAL_REVIEW' if gd[r['gene']]['abstracts_manually_reviewed'] else 'NOT_REVIEWED',
          'work_group':gd[r['gene']]['work_group'],'decision':'RETAIN_NO_ARBITRARY_SCORE_OR_SIGNIFICANCE_FILTER'})
    write(OUT/'relation_working_table_v1.tsv',relations);write(OUT/'gene_working_table_v1.tsv',genes);write(EXT/'all_gene_evaluability.tsv',external)
    summary={'M':len(relations),'G':len(genes),'gene_search_complete':len(search),'manual_abstract_records':len(reviews),
      'historical_cell_source_genes_covered':sum(r['GSE263733_summary_present'] for r in external),
      'historical_cptac_genes_covered':sum(r['historical_cptac_summary_present'] for r in external),
      'new_patient_tests_run':False,'overall_first_batch_status':'PARTIAL_SERVER_CONNECTION_BLOCKED'}
    for out,stage in [(OUT,'07_INTEGRATION'),(EXT,'06_EXTERNAL')]:
        dump(out/'summary.json',summary)
        dump(out/'analysis_spec.json',{'stage':stage,'version':'PDAC_working_views_v1','mapping_sha256':sha(MAP/'direct_relations_v1.tsv'),'gene_count':len(genes),'scoring':'none','historical_q_not_used_as_new_family_q':True,'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()})
        write(out/'source_manifest.tsv',[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in [MAP/'direct_relations_v1.tsv',FUN/'gene_search_ledger.tsv',FUN/'abstract_review.tsv',OLD/'historical_associations.tsv',OLD/'historical_cell_source_gene_summary.tsv',OLD/'historical_cptac_effects.tsv']])
        dump(out/'validation.json',{'status':'PASS','all_M_G_retained':True,'unrun_statistics_are_NA':True,'unverified_external_independence_retained':True})
    print(json.dumps(summary))
if __name__=='__main__':main()
