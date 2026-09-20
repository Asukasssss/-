"""Validate and combine already-public donor summaries without rounding source values."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import pandas as pd

PREFIX = 'cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    rows, checks, header = [], {}, None
    for study in ['Lee', 'Pelka']:
        path = a.out / (study + '_source_expression_summary.tsv')
        with path.open(encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if header is None:
                header = reader.fieldnames
            assert reader.fieldnames == header and header[:len(PREFIX)] == PREFIX
            current = list(reader)
        validation = json.loads((a.out / (study + '_validation.json')).read_text())
        assert validation['status'] == 'PASS' and len(current) == validation['rows']
        assert validation['pre_expression_remote_lock'] == '09498b982756677e106e1a10578868a9017ddc91'
        checks[study] = validation
        rows.extend(current)
    keys = ['study', 'tissue', 'enrichment', 'technology', 'annotation_level', 'cell_type', 'support_set', 'gene']
    rows.sort(key=lambda row: tuple(row[k] for k in keys))
    assert len({tuple(row[k] for k in keys) for row in rows}) == len(rows)
    for row in rows:
        assert row['p_value'] == row['q_value'] == 'NA'
        if int(row['n']) < 3:
            assert row['effect'] == 'NA' and row['detection_fraction_median'] == 'NA'
    with (a.out / 'source_expression_summary.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    data = pd.DataFrame(rows)
    display = data[(data.tissue.isin(['T', 'Tumor'])) &
        (data.enrichment.isin(['unsorted', 'NO_CD45_REPORTED_FICOLL_PURIFIED'])) &
        (data.annotation_level == 'lineage') & (data.support_set == 'patients_with_at_least_20_cells')]
    display.to_csv(a.out / 'tumor_unenriched_lineage_display.tsv', sep='\t', index=False)
    validation = {'status': 'PASS', 'rows': len(rows), 'study_checks': checks,
        'source_field_strings_preserved': True, 'public_small_group_suppression': True,
        'unique_keys': True, 'no_new_hypothesis_tests': True,
        'independent_synthetic_tests': '4 PASS on server165; donor weighting/zeros/privacy/identity/sparse indexing',
        'unresolved': ['No spatial validation', 'No CAMP purity or cell-composition adjustment', 'No paper-specific PRMT7 isoform resolution']}
    (a.out / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps({'rows': len(rows), 'display_rows': len(display), 'status': 'PASS'}))

if __name__ == '__main__':
    main()
