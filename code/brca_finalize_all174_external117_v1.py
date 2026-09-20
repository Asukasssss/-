"""Aggregate-only integration; prior117 columns and CAMP statistics remain unchanged."""
from pathlib import Path
import csv,json,hashlib,math
ROOT=Path(__file__).resolve().parents[1]
RUN='20260920T090922Z_all174_external117_resources_v1'
OUT=ROOT/'results/BRCA/06_EXTERNAL'/RUN
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)
def write(p,cols,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def num(x):
    try:return float(x)
    except:return float('nan')
def main():
    cols,rows=read(OUT/'external_all174_associations.tsv')
    _,old=read(ROOT/'reference/brca/association_summary_v1.tsv')
    old={r['relation_id']:r for r in old if r['analysis_type']=='author_processed_primary'};assert len(old)==174
    append=['CAMP_unadjusted_rho','CAMP_unadjusted_q','direction_vs_CAMP_unadjusted','comparison_interpretation']
    for r in rows:
        c=old[r['relation_id']];a=num(r['effect']);b=num(c['rho'])
        sign='NOT_EVALUABLE' if not(math.isfinite(a) and math.isfinite(b)) else ('SAME_SIGN' if a*b>0 else 'OPPOSITE_OR_ZERO')
        r.update(CAMP_unadjusted_rho=c['rho'],CAMP_unadjusted_q=c['q_value'],direction_vs_CAMP_unadjusted=sign,comparison_interpretation='descriptive_sign_comparison;CAMP_specimens_FUSCC_author_patients;TNBC_and_platform_differences;no_formal_effect_difference_test')
    write(OUT/'external_all174_with_CAMP_context.tsv',cols+append,rows)
    prim=[r for r in rows if r['test_family']=='FUSCC_PRIMARY174']
    assert len(prim)==174
    prior=ROOT/'results/BRCA/06_EXTERNAL/20260920T050412Z_pycr1_external_v1/all117_comparison_deep_appended.tsv'
    gc,gr=read(prior);assert len(gr)==117
    _,inv=read(OUT/'all117_intervention_resource_inventory.tsv');iv={r['gene']:r for r in inv}
    added=['external174_status','external174_n_evaluable','external174_n_q_lt005','external174_significant_relations','external174_significant_directions_vs_CAMP','resource117_status','resource117_accession_leads','resource117_scope','external174_source']
    for g in gr:
        a=[r for r in prim if r['gene']==g['gene']];ok=[r for r in a if r['status']=='DONE'];sig=[r for r in ok if num(r['q_value'])<.05]
        g.update(external174_status='DONE' if len(ok)==len(a) else ('PARTIAL' if ok else 'NOT_EVALUABLE'),external174_n_evaluable=str(len(ok)),external174_n_q_lt005=str(len(sig)),external174_significant_relations=';'.join(r['relation_id'] for r in sig),external174_significant_directions_vs_CAMP=';'.join(r['relation_id']+':'+r['direction_vs_CAMP_unadjusted'] for r in sig),resource117_status=iv[g['gene']]['status'],resource117_accession_leads=iv[g['gene']]['accession_leads'],resource117_scope='existing_literature_machine_discovery;see manual_resource_design_notes.tsv;not 117 verified intervention datasets',external174_source=OUT.relative_to(ROOT).as_posix())
    write(OUT/'all117_external_resources_appended.tsv',gc+added,gr)
    _,original=read(prior);assert [{c:r[c] for c in gc} for r in gr]==original
    summary={'all117_retained':True,'historical_columns_unchanged':len(gc),'primary_evaluable':sum(r['status']=='DONE' for r in prim),'primary_q_lt005':sum(num(r['q_value'])<.05 for r in prim),'primary_significant_genes':len({r['gene'] for r in prim if num(r['q_value'])<.05}),'significant_same_sign_CAMP_unadjusted':sum(num(r['q_value'])<.05 and r['direction_vs_CAMP_unadjusted']=='SAME_SIGN' for r in prim),'significant_opposite_or_zero_CAMP_unadjusted':sum(num(r['q_value'])<.05 and r['direction_vs_CAMP_unadjusted']=='OPPOSITE_OR_ZERO' for r in prim)}
    (OUT/'integration_validation.json').write_bytes((json.dumps(summary,indent=2)+'\n').encode());print(json.dumps(summary))
    sf,sr=read(ROOT/'coordination/stages/BRCA.tsv');sr=[r for r in sr if r['run_id']!=RUN]
    for stage,status,scope in [('05_FUNCTION','PARTIAL','All117 existing-literature resource discovery; design verification incomplete; PYCR1 CAF review closed'),('06_EXTERNAL','DONE','All174 feasibility and 2 fixed-family FUSCC analyses; missing entries retained'),('07_INTEGRATION','DONE','All117 aggregate comparison appended; historical fields unchanged')]:
        vals=['BRCA',stage,RUN,'external174_resources_v1',status,scope,OUT.relative_to(ROOT).as_posix(),'code/brca_finalize_all174_external117_v1.py','analysis/brca-functional-review-20260919','Author-annotation identity; TNBC context; accession is not intervention proof','Use exact relationship and model context; verify resource design before new perturbation tests']
        sr.append(dict(zip(sf,vals)))
    write(ROOT/'coordination/stages/BRCA.tsv',sf,sorted(sr,key=lambda r:(r['stage_id'],r['run_id'])))
    for name in ['identity_coverage.tsv','gene_identity_coverage.tsv']:
        fields, records = read(OUT/name)
        write(OUT/name, fields, [{k: ('NA' if v == '' else v) for k,v in r.items()} for r in records])
    scripts=['brca_all174_external_v1.py','brca117_perturb_resource_inventory_v1.py','brca_verify_all174_external_v1.py','brca_finalize_all174_external117_v1.py','brca_pycr1_context_v1.py','brca_sc117_profile_v1.py']
    write(OUT/'code_manifest.tsv',['path','sha256'],[{'path':'code/'+s,'sha256':hashlib.sha256((ROOT/'code'/s).read_bytes()).hexdigest()} for s in scripts])
    for p in OUT.rglob('*'):
        if p.is_file() and p.suffix in {'.md','.tsv','.json','.txt'}:p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
    ch=OUT/'checksums.sha256';ch.write_bytes(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(OUT).as_posix()}\n' for p in sorted(OUT.rglob('*')) if p.is_file() and p!=ch).encode())
if __name__=='__main__':main()
