"""Finalize aggregate-only BRCA delivery; no patient calculations."""
import csv, hashlib, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = '20260920T050412Z_pycr1_external_v1'
OUT = ROOT / 'results/BRCA/06_EXTERNAL' / RUN
def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        r = csv.DictReader(f, delimiter='\t'); return list(r.fieldnames), list(r)
def write(p, fields, rows):
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    fields, rows = read(OUT/'pycr1_celltype_contrasts.tsv')
    prefix, _ = read(ROOT/'templates/statistical_result.tsv')
    for r in rows:
        r.update(cancer='BRCA', stage_id='06_EXTERNAL', run_id=RUN, analysis_version='pycr1_context_v1', analysis_type='paired_source_correlation_difference_descriptive', gene='PYCR1', unit='source_donor_label', n=r['n_common_sources'], test_family='DESCRIPTIVE_NO_P', source_id=r['cohort'])
    write(OUT/'pycr1_celltype_contrasts.tsv', prefix+[x for x in fields if x not in prefix], [{k:r.get(k,'NA') for k in prefix+[x for x in fields if x not in prefix]} for r in rows])
    old = ROOT/'results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/all117_comparison_robustness_appended.tsv'
    cols, rows = read(old)
    assert len(rows)==117 and len({r['gene'] for r in rows})==117
    updates = {
        'GPCPD1':('DONE','FUSCC n258 rho=-0.040575 q=1; strong negative CAMP association not reproduced in this TNBC context','Retain CAMP and functional evidence; do not claim general external replication'),
        'GPI':('NOT_EVALUABLE','Exact glucose-6-phosphate absent from available FUSCC polar annotation; available CPTAC lacks target metabolite matrix','Retain candidate; missing coverage is not a negative association'),
        'PYCR1':('DONE','Primary8: no q<0.05; cancer E2F and fibroblast collagen directions concordant across cohorts; Wu untreated cancer E2F and collagen q<0.05','Exploratory transcriptional context only; no demonstrated cell-specific metabolic function')}
    newcols=['deep_v1_status','deep_v1_result','deep_v1_decision','deep_v1_source']
    for r in rows:
        values=updates.get(r['gene'],('NOT_RUN','Not tested in this bounded batch','Retain previous evidence and priority'))
        r.update(dict(zip(newcols,(*values,str(OUT.relative_to(ROOT)).replace('\\','/')))))
    target=OUT/'all117_comparison_deep_appended.tsv';write(target,cols+newcols,rows)
    _,check=read(target);_,original=read(old)
    assert [{k:r[k] for k in cols} for r in check]==original
    stages=ROOT/'coordination/stages/BRCA.tsv';sf,sr=read(stages)
    sr=[r for r in sr if r['run_id']!=RUN]
    for stage,scope in [('06_EXTERNAL','Two external relations and limited PYCR1 program analysis completed'),('07_INTEGRATION','All117 retained; all historical columns exactly preserved; three genes annotated')]:
        vals=['BRCA',stage,RUN,'pycr1_external_v1','DONE',scope,str(OUT.relative_to(ROOT)).replace('\\','/'),'code/brca_finalize_pycr1_external_v1.py','analysis/brca-functional-review-20260919','Bounded batch complete; missing metabolite coverage and exploratory inference retained','Read README_CN.md; no new all117 screening']
        sr.append(dict(zip(sf,vals)))
    write(stages,sf,sorted(sr,key=lambda r:(r['stage_id'],r['run_id'])))
    numerical=json.loads((OUT/'numeric_validation.json').read_text())
    validation={'status':'DONE','scope':'bounded external and PYCR1 batch','numerical_validation':numerical,'all117_retained':True,'historical_columns_exactly_preserved':len(cols),'primary_tests':8,'untreated_sensitivity_tests':4,'external_planned_tests':2,'not_verified':['genotype-confirmed patient independence','same tissue aliquot for external modalities','causal function or metabolite flux'],'patient_measurements_exported':False}
    (OUT/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    programs=ROOT/'code/brca_pycr1_programs_v1.json'
    programs.write_bytes(programs.read_bytes().replace(b'\r\n',b'\n'))
    files=[ROOT/'code'/n for n in ['brca_pycr1_context_v1.py','brca_external_two_relations_v1.py','brca_verify_pycr1_external_v1.py','brca_finalize_pycr1_external_v1.py','brca_pycr1_programs_v1.json','brca_sc117_profile_v1.py','brca_sc117_Wu2021_config.json','brca_sc117_Pal2021_config.json']]
    files += [old,OUT/'analysis_spec.json',OUT/'Wu2021_input.json',OUT/'Pal2021_reprocessed_input.json',OUT/'external_source_manifest.tsv']
    write(OUT/'source_manifest.tsv',['kind','path','sha256'],[{'kind':'code_or_input_provenance','path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':sha(p)} for p in files])
    manifest=OUT/'checksums.sha256'
    for p in OUT.rglob('*'):
        if p.is_file() and p.suffix in {'.json','.tsv','.md','.txt'}:
            p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
    # Recompute provenance after text normalization.
    write(OUT/'source_manifest.tsv',['kind','path','sha256'],[{'kind':'code_or_input_provenance','path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':sha(p)} for p in files])
    manifest.write_bytes(''.join(f'{sha(p)}  {p.relative_to(OUT).as_posix()}\n' for p in sorted(OUT.rglob('*')) if p.is_file() and p!=manifest).encode('utf-8'))
    print(json.dumps(validation,ensure_ascii=False))
if __name__=='__main__':main()
