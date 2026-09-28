"""Verify and integrate public first-round summaries, without patient-level data."""
import csv,json,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from patient_first_round import fixed_bh
ROOT=Path(__file__).resolve().parents[2]
PAT=ROOT/'results/PDAC/03_PATIENT/20260920T113000Z_first_round_v1'
MAP=ROOT/'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1'
PREV=ROOT/'results/PDAC/07_INTEGRATION/20260920T112500Z_working_v1'
OUT=ROOT/'results/PDAC/07_INTEGRATION/20260920T114000Z_first_round_v1'
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def write(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def dump(p,o):p.write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    rs=read(PAT/'associations.tsv');direct=read(MAP/'direct_relations_v1.tsv');resolution=read(PAT/'gene_resolution.tsv');loo=read(PAT/'leave_one_out.tsv')
    assert len(rs)==len(loo)==724 and len(direct)==362 and len(resolution)==289
    bylabel=defaultdict(list)
    for r in resolution:
        if r['status']=='DONE':bylabel[r['rna_label']].append(r['gene'])
    assert all(len(v)==1 for v in bylabel.values()),'Different genes collapse to same RNA label;review before accepting'
    result={(r['relation_id'],r['analysis_type']):r for r in rs};checks=[]
    assert len(result)==724
    for kind in ['primary','availability']:
        subset=[r for r in rs if r['analysis_type']==kind];expected=fixed_bh([float(r['p_value']) if r['status']=='DONE' else None for r in subset])
        for r,q in zip(subset,expected):
            assert int(r['family_n_planned'])==362
            if q is None:assert r['p_value']==r['q_value']=='NA'
            else:assert np.isclose(float(r['q_value']),q,rtol=0,atol=1e-12)
    for rel in direct:
        a=result[(rel['relation_id'],'primary')];b=result[(rel['relation_id'],'availability')]
        assert a['original_camp_q']==rel['original_q'] and a['original_camp_effect']==rel['original_effect']
        if a['input_hash']==b['input_hash']:
            assert all(a[k]==b[k] for k in ['effect','p_value','ci_lower','ci_upper','n','bootstrap_valid'])
            checks.append(rel['relation_id'])
    assert sha(ROOT/'code/pdac/patient_first_round.py')==json.loads((PAT/'analysis_spec.json').read_text())['script_sha256']
    primary=[r for r in rs if r['analysis_type']=='primary'];sig=[r for r in primary if r['status']=='DONE' and float(r['q_value'])<.05]
    validation={'status':'PASS','all362x2_rows_present':True,'unique_RNA_resolution':True,'fixed_M_BH_recomputed_locally_matches':True,
                'identical_input_records_reused':len(checks),'computable_identical_input_statistics_reused':sum(result[(k,'primary')]['status']=='DONE' for k in checks),'frozen_effect_q_strings_preserved':True,'server_script_hash_matches':True,
                'patient_independence':'NOT_SEPARATELY_VERIFIED'}
    dump(PAT/'acceptance_validation.json',validation)
    OUT.mkdir(parents=True,exist_ok=False)
    work=read(PREV/'relation_working_table_v1.tsv');genes=read(PREV/'gene_working_table_v1.tsv');ld={(r['relation_id'],r['analysis_type']):r for r in loo}
    comp=[]
    for r in work:
        a=result[(r['relation_id'],'primary')];b=result[(r['relation_id'],'availability')];l=ld[(r['relation_id'],'primary')]
        r.update(patient_status=a['status'],primary_n=a['n'],primary_rho=a['effect'],primary_p=a['p_value'],primary_q=a['q_value'],
                 availability_n=b['n'],availability_rho=b['effect'],availability_q=b['q_value'],leave_one_out_status='DONE' if int(l['loo_n_valid'])>0 else 'NOT_EVALUABLE',
                 primary_ci_lower=a['ci_lower'],primary_ci_upper=a['ci_upper'],availability_ci_lower=b['ci_lower'],availability_ci_upper=b['ci_upper'],
                 primary_reason=a['reason'],availability_reason=b['reason'],loo_rho_min=l['loo_rho_min'],loo_rho_max=l['loo_rho_max'],loo_max_abs_change=l['loo_max_abs_change'],loo_sign_changes=l['loo_sign_changes'],
                 n_preimputation_available=a['n_preimputation_available'],current_test_family=a['test_family'])
        support=a['status']=='DONE' and float(a['q_value'])<.05
        if r['function_status']=='ABSTRACT_GENETIC_REPORT_FULLTEXT_PENDING':r['work_group']='FUNCTIONAL_BACKGROUND_ASSOCIATION_TO_INTERPRET'
        else:r['work_group']='RELATION_EVIDENCE_REVIEW_PENDING'
        r['decision']='RETAIN;no significance-based exclusion;functional and source evidence still incomplete'
        comp.append({'relation_id':r['relation_id'],'feature_name':r['feature_name'],'gene':r['gene'],'primary_n':a['n'],'availability_n':b['n'],
          'primary_rho':a['effect'],'availability_rho':b['effect'],'rho_difference':float(b['effect'])-float(a['effect']) if a['status']==b['status']=='DONE' else 'NA',
          'primary_ci_lower':a['ci_lower'],'primary_ci_upper':a['ci_upper'],'availability_ci_lower':b['ci_lower'],'availability_ci_upper':b['ci_upper'],
          'same_inputs':a['input_hash']==b['input_hash'],'same_sign':str(np.sign(float(a['effect']))==np.sign(float(b['effect']))) if a['status']==b['status']=='DONE' else 'NA',
          'primary_status':a['status'],'availability_status':b['status'],'availability_reason':b['reason']})
    for g in genes:
        rr=[r for r in primary if r['gene']==g['gene']]
        g.update(patient_relations_evaluable=sum(r['status']=='DONE' for r in rr),patient_relations_q_lt_005=sum(r['status']=='DONE' and float(r['q_value'])<.05 for r in rr),
                 next_evidence='Review full candidate evidence;continue all-G cell source/CPTAC feasibility;do not filter solely by significance')
    write(OUT/'relation_working_table_v1.tsv',work);write(OUT/'gene_working_table_v1.tsv',genes);write(PAT/'primary_availability_comparison.tsv',comp)
    summary=json.loads((PAT/'summary.json').read_text());summary.update(identical_inputs_reused=len(checks),primary_nominal_p_lt_005=sum(r['status']=='DONE' and float(r['p_value'])<.05 for r in primary),
      primary_q_significant_genes=len({r['gene'] for r in sig}),primary_q_significant_features=len({r['metabolite_name'] for r in sig}),
      RNA_missing_genes=[r['gene'] for r in resolution if r['status']!='DONE'],loo_primary_evaluable=sum(int(r['loo_n_valid'])>0 and r['analysis_type']=='primary' for r in loo))
    dump(OUT/'summary.json',summary);dump(OUT/'validation.json',validation)
    dump(OUT/'analysis_spec.json',{'version':'PDAC_first_round_integrated_v1','arbitrary_score':False,'no_significance_exclusion':True,'M':362,'G':289,'original_q_separate':True})
    write(OUT/'source_manifest.tsv',[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in [PAT/'associations.tsv',PAT/'leave_one_out.tsv',PAT/'gene_expression_descriptive.tsv',PREV/'relation_working_table_v1.tsv',PREV/'gene_working_table_v1.tsv']])
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
