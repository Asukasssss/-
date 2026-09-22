"""Independent relational and boundary audit of public mapping outputs."""
import argparse,csv,json,hashlib
from pathlib import Path
from collections import Counter


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);a=ap.parse_args();p=a.results
    root=Path(__file__).resolve().parents[2]
    original={r['metabolite_name']:r for r in read(root/'results/COAD/03_PATIENT/20260921T123200Z_paired_metabolites33_v1/primary33.tsv')}
    allrows=read(p/'results.tsv');mainrows=read(p/'main_relations.tsv');ledger=read(p/'feature_ledger159.tsv');support=read(p/'supporting_evidence.tsv')
    key=lambda x:(x['metabolite_name'],x['metabolite_key'],x['human_gene_id'])
    assert len(allrows)==len({key(r) for r in allrows})==1466
    assert {key(r) for r in mainrows}=={key(r) for r in allrows if r['pool']=='DIRECT'}
    assert len(mainrows)==891
    assert len(ledger)==159 and sum(r['selected_paired_P_lt_005']=='True' for r in ledger)==101
    for r in ledger:
        old=original[r['metabolite_name']]
        for f,g in [('paired_P','p_value'),('paired_q','q_value'),('original_q','original_q'),('original_effect','original_effect'),('up_pairs','up_pairs'),('down_pairs','down_pairs'),('equal_pairs','equal_pairs')]:assert r[f]==old[g]
    supported={key(s) for s in support if json.loads(s['target_chebi_matches']) and s['modification_context']=='False' and s['gtp_cycle_context']=='False'}
    assert {key(r) for r in mainrows}<=supported
    assert all(not r['identity_review_status'].startswith('HOLD_') for r in mainrows)
    assert all(r['current_patient_status']=='NOT_RUN' and r['current_patient_P']==r['current_patient_q']=='NA' for r in allrows)
    by={(r['metabolite_name'],r['gene']):r for r in allrows}
    positive=[('pantothenate','PANK1'),('pantothenate','SLC5A6'),('5-oxoproline','OPLAH'),('xanthine','XDH'),('succinate','SDHB'),('malonate','ACSF3'),('phenylalanine','FARSB')]
    for k in positive:assert by[k]['pool']=='DIRECT'
    negative=[('S-adenosylmethionine (SAM)','KMT5A'),('S-adenosylmethionine (SAM)','PRMT7'),('S-adenosylmethionine (SAM)','ICMT'),("guanosine 5'- diphosphate (GDP)",'RAB1A'),('succinate','PLOD1'),('serine','SRR')]
    for k in negative:assert by[k]['pool']=='CONDITIONAL'
    for k in [('succinate','SDHB'),('phenylalanine','FARSB')]:assert 'COMPLEX_MEMBER' in by[k]['relation_types']
    for r in mainrows:
        old=original[r['metabolite_name']]
        assert r['paired_metabolite_P']==old['p_value'] and r['paired_metabolite_q']==old['q_value']
    # Recount all feature dispositions from relation rows, independently of summary values.
    for r in ledger:
        c=Counter(x['pool'] for x in allrows if x['metabolite_name']==r['metabolite_name'])
        assert int(r['direct_relations'])==c['DIRECT'] and int(r['conditional_relations'])==c['CONDITIONAL'] and int(r['unresolved_relations'])==c['UNRESOLVED']
    spec=json.loads((p/'analysis_spec.json').read_text());assert spec['planned_patient_relations']==891 and spec['main_relation_file_sha256']==hashlib.sha256((p/'main_relations.tsv').read_bytes()).hexdigest()
    audit=dict(status='PASS',all159_source_statistic_strings_preserved=True,all1466_relation_keys_unique=True,main891_exact_support_join=True,
        all_feature_counts_reconciled=True,no_current_patient_stats_claimed=True,positive_boundary_controls=len(positive),conditional_boundary_controls=len(negative),complex_labels_checked=True)
    (p/'independent_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit))


if __name__=='__main__':main()
