"""Reproduce bounded receipt checks and append receiver notes; no patient tests.

Default outputs must not exist. Use --output-root NEW_DIRECTORY to reproduce
without overwriting this receipt. Bibliographic metadata and manual notes are
versioned inputs; this script does not certify experimental claims.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import re
from pathlib import Path

IMPORTED = '20260919T120437Z_functional_review_v4'
RECEIPT = '20260919T130605Z_v4_reception'
ORIGINAL_TEN = set('GPCPD1 PNP GPI ASNS KYNU PCYT2 ETNK1 NNMT SORD SLC7A11'.split())


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader), reader.fieldnames


def write(path, rows, fields):
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(value):
    return None if value in ('', 'NA', 'nan', None) else float(value)


def same(a, b):
    a, b = number(a), number(b)
    return a is None and b is None or a is not None and b is not None and math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    incoming = root / 'results/BRCA/05_FUNCTION' / IMPORTED
    receipt = root / 'results/BRCA/05_FUNCTION' / RECEIPT
    merged = root / 'results/BRCA/07_INTEGRATION' / RECEIPT
    output = args.output_root or merged
    legacy_path = root / 'results/BRCA/20260919_ER_availability_v3/candidate_comparison_117_integrated_v3.tsv'
    edge_path = root / 'results/BRCA/20260919_ER_availability_v3/er_availability_all_174.tsv'
    old, old_fields = read(legacy_path)
    reviews, _ = read(incoming / 'candidate_comparison_117_functional_v4.tsv')
    full, fields = read(merged / 'candidate_comparison_117_legacy_plus_v4.tsv')
    edges, _ = read(edge_path)
    compact, _ = read(incoming / 'relations_174_annotated_v4.tsv')
    studies, _ = read(incoming / 'functional_studies_v4.tsv')
    searches, _ = read(incoming / 'search_log_v4.tsv')
    manifest, _ = read(receipt / 'incoming_PACKAGE_SHA256.tsv')
    metadata = json.loads((receipt / 'bibliographic_metadata_snapshot.json').read_text(encoding='utf8'))
    notes, _ = read(receipt / 'receiver_source_notes.tsv')
    claims, _ = read(receipt / 'key_claim_checks.tsv')
    by_gene = {row['gene']: row for row in full}
    assert len(full) == len(by_gene) == len(old) == len(reviews) == 117
    assert [x['gene'] for x in old] == [x['gene'] for x in full], 'Legacy row order changed'
    assert set(by_gene) == {x['gene'] for x in reviews}
    assert all(row[k] == by_gene[row['gene']][k] for row in old for k in old_fields), 'Legacy cell changed'
    edge_by = {row['relation_id']: row for row in edges}
    compact_by = {row['relation_id']: row for row in compact}
    assert len(edges) == len(compact) == len(edge_by) == len(compact_by) == 174
    assert set(edge_by) == set(compact_by)
    for key, row in compact_by.items():
        original = edge_by[key]
        for field in ['gene', 'metabolite']:
            assert row[field] == original[field], (key, field)
        for a, b in [('n_specimens', 'n'), ('er_available_partial_rank_rho', 'partial_rank_rho'), ('er_available_q', 'q_value'), ('er_primary_q', 'primary_er_q')]:
            assert same(row[a], original[b]), (key, a)
    imported_checked = 0
    for row in manifest:
        if row['path'].startswith(('code/', 'results/')):
            path = root / row['path']
            assert path.stat().st_size == int(row['bytes']) and digest(path) == row['sha256'], row['path']
            imported_checked += 1
    assert len(studies) == len({row['study_id'] for row in studies}) == 122
    assert len({row['source_url'] for row in studies}) == 121
    assert len(searches) == 127 and {row['gene'] for row in searches} == set(by_gene)
    genetic = {g.strip() for row in studies if row['evidence_class'] == 'GENETIC_BRCA' for g in row['genes'].split(';')}
    # Input uses semicolon-separated gene symbols; a source may support >1 gene.
    assert genetic <= set(by_gene) and len(genetic) == 52
    assert len(genetic - ORIGINAL_TEN) == 44
    priorities = dict(collections.Counter(row['work_priority'] for row in reviews))
    assert priorities == {'保留': 57, '暂挂': 36, '优先核查': 24}
    assert all(row['depmap_score'] == row['depmap_n_models'] == '' for row in reviews)
    meta = {row['pmid']: row for row in metadata}
    audit = []
    for row in studies:
        match = re.search(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)', row['source_url'])
        pmid = match[1] if match else ''
        record = meta.get(pmid, {})
        if pmid:
            assert record, pmid
        audit.append(dict(cancer='BRCA', study_id=row['study_id'], genes=row['genes'], source_url=row['source_url'], pmid=pmid or 'NA', doi=record.get('doi', 'NA'), canonical_title=record.get('title', 'NA'), status='IDENTIFIER_RESOLVED_NOT_CLAIM_VERIFIED' if pmid else 'NOT_CHECKED_IN_PMID_BATCH', publication_types=';'.join(record.get('pubtypes', [])) or 'NA', linked_notices=json.dumps(record.get('related_notices', []), ensure_ascii=False), statistical_unit='publication_record', n_patient='NA', new_patient_tests=0))
    note_by = collections.defaultdict(list)
    for row in notes:
        for gene in row['genes'].split(';'):
            assert gene in by_gene
            note_by[gene].append(row)
    claim_by = collections.defaultdict(list)
    for row in claims:
        claim_by[row['gene']].append(row)
    new_fields = ['receiver_v4_status', 'receiver_v4_claim_check', 'receiver_v4_source_note_ids', 'receiver_v4_source_cautions', 'receiver_v4_source_urls', 'receiver_v4_next_action', 'receiver_v4_genetic_report_count_basis', 'receiver_v4_original_ten_member']
    assert not set(fields) & set(new_fields)
    updated = []
    for row in full:
        gene = row['gene']; ns = note_by[gene]; cs = claim_by[gene]
        extra = dict(zip(new_fields, [
            'RECEIVED_WITH_SOURCE_CAUTION' if ns else 'RECEIVED_BOUNDED_REVIEW',
            ' | '.join(c['finding'] for c in cs) or 'NOT_IN_TARGETED_CLAIM_CHECK;imported interpretation not independently certified',
            ';'.join(n['note_id'] for n in ns) or 'NA',
            ' | '.join(n['finding'] for n in ns) or 'No additional issue identified within limited receipt checks;not a clean bill for all literature',
            ';'.join(n['notice_url'] for n in ns) or 'NA',
            ' | '.join(n['action'] for n in ns) or row['review_v4_next_minimal_action'],
            'IMPORTED_GENETIC_BRCA_CLASS' if gene in genetic else 'NOT_IN_IMPORTED_52;not_negative',
            str(gene in ORIGINAL_TEN),
        ]))
        updated.append({**row, **extra})
    write(output / 'candidate_comparison_117_receiver_v4.tsv', updated, fields + new_fields)
    reading = [dict(cancer='BRCA', gene=r['gene'], direct_metabolites=r['review_v4_linked_direct_metabolites'], imported_functional_class=r['review_v4_functional_class'], imported_work_priority=r['review_v4_work_priority'], patient_evidence=r['review_v4_patient_evidence_summary'], original_best_q=r['best_q'], receiver_source_cautions=r['receiver_v4_source_cautions'], next_minimal_action=r['receiver_v4_next_action'], study_ids=r['review_v4_study_ids'], source_urls=r['review_v4_source_urls']) for r in updated]
    write(output / 'candidate_reading_view_117.tsv', reading, list(reading[0]))
    write(output / 'bibliographic_id_check_122.tsv', audit, list(audit[0]))
    inputs = [legacy_path, edge_path, incoming / 'functional_studies_v4.tsv', incoming / 'candidate_comparison_117_functional_v4.tsv', incoming / 'relations_174_annotated_v4.tsv', receipt / 'bibliographic_metadata_snapshot.json', receipt / 'receiver_source_notes.tsv', receipt / 'key_claim_checks.tsv', root / 'code/brca_receive_functional_v4.py']
    validation = dict(passed=True, scope='Structural/numerical receipt and bounded literature annotation;not full biological verification', genes=117, relations=174, legacy_columns_preserved=len(old_fields), legacy_cells_preserved=117 * len(old_fields), imported_review_columns=len(fields)-len(old_fields), receiver_columns=len(new_fields), retained_package_files_hash_checked=imported_checked, studies=122, unique_source_urls=121, directed_search_records=127, searched_genes=117, genetic_effect_genes_by_imported_class=len(genetic), genetic_effect_genes_outside_original_ten=len(genetic-ORIGINAL_TEN), any_genetic_report_genes=sum(x['genetic_intervention_reported']=='True' for x in reviews), priority_counts=priorities, pmid_rows=sum(x['pmid']!='NA' for x in audit), unique_pmids=len(meta), targeted_claim_records=len(claims), source_notice_records=len(notes), new_patient_statistics=0, depmap_calculations=0, source_sha256={p.relative_to(root).as_posix():digest(p) for p in inputs}, unverified=['All 122 experimental claims and full texts', 'Full publication integrity for non-PMID sources or unindexed notices', 'Independent patient identity, subtype/purity/batch and original detection masks', 'CAMP metabolite mediation', 'Consequences of notices explicitly marked NEEDS_REVIEW'])
    vp = output / 'validation.json'
    if vp.exists():
        raise FileExistsError(vp)
    vp.write_text(json.dumps(validation, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k:v for k,v in validation.items() if k not in ('source_sha256','unverified')}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
