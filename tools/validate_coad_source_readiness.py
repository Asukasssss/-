"""Verify aggregate source-audit delivery without copying patient matrices locally."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);a=p.parse_args()
    batch=(ROOT/'results/COAD/03_PATIENT'/a.run_id).resolve()
    assert batch.parent==(ROOT/'results/COAD/03_PATIENT').resolve()
    validation=json.loads((batch/'validation.json').read_text())
    assert validation['status']=='PASS'
    for name,expected in validation['output_sha256'].items():
        assert sha(batch/name)==expected
    records=read(batch/'results.tsv')
    template=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t')
    assert list(records[0])[:len(template)]==template
    frozen={(r['feature_name'],r['metabolite_key']):r for r in read(ROOT/'reference/camp/cancer_effects.tsv') if r['cancer']=='COAD'}
    assert len(records)==974
    keys=[]
    for row in records:
        source=frozen[(row['metabolite_name'],row['metabolite_key'])]
        assert (row['original_effect'],row['original_q'])==(source['hedges_g'],source['effect_fdr'])
        assert row['effect']==row['p_value']==row['q_value']=='NA'
        assert row['ready_for_patient_testing']=='False'
        assert row['stage_id']=='03_PATIENT' and row['analysis_type']=='source_readiness_no_association'
        if row['input_coverage']=='BOTH_PRESENT':
            assert row['n']=='37' and row['status']=='NOT_RUN'
        else:
            assert row['n']=='NA' and row['status']=='NOT_EVALUABLE'
        if row['gene']=='NA':assert row['input_coverage']=='GENE_SYMBOL_UNRESOLVED'
        keys.append(tuple(row[k] for k in ['cohort','metabolite_name','metabolite_key','gene','human_gene_id']))
    assert keys==sorted(keys) and len(set(keys))==len(keys)
    summary=json.loads((batch/'summary.json').read_text())
    assert dict(Counter(r['input_coverage'] for r in records))==summary['pair_input_coverage_counts']
    assert len({r['human_gene_id'] for r in records})==summary['n_candidate_genes']
    identity=read(batch/'identity_source_review.tsv')
    assert len(identity)==73 and all(r['exact_source_ids_equal_frozen']=='True' for r in identity)
    features=read(batch/'feature_coverage.tsv');assert len(features)==159
    assert all(r['tumor_matrix_found']=='True' for r in features)
    headings=[line[3:] for line in (batch/'README_CN.md').read_text(encoding='utf-8').splitlines() if line.startswith('## ')]
    assert headings==['本轮问题','输入与范围','实际结果','新手解释','限制/反证','当前决定','下一步','复现命令']
    index=read(ROOT/'coordination/stages/COAD.tsv')
    assert list(index[0])==(ROOT/'templates/stage_index.tsv').read_text().strip().split('\t')
    assert [(r['stage_id'],r['run_id']) for r in index]==sorted((r['stage_id'],r['run_id']) for r in index)
    print('PASS: source-audit hashes, common format, counts, coverage categories, frozen values, unique ordered keys and stage index.')


if __name__=='__main__':main()
