"""Server-only COAD source readiness audit; export aggregate counts, never sample values.

Run in a new B-owned directory with an exclusive .running lock. Source files are read-only.
No effect/q recalculation, identifier inference, imputation, correlation or patient pairing.
"""
import argparse
import csv
import hashlib
import json
import platform
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import openpyxl


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, obj):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(obj, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def tsv(path, fields, rows):
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def clean(value):
    return 'NA' if pd.isna(value) or str(value).strip() == '' else str(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project-root', type=Path, required=True)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--code-commit', required=True)
    args = parser.parse_args()
    root, out = args.project_root.resolve(), args.run_dir.resolve()
    expected_parent = root / 'results/collaborative/COAD/B'
    assert str(root).startswith('/public3/') and out.parent == expected_parent
    assert (out / '.running').read_text().strip() == out.name
    assert not (out / 'summary.json').exists(), 'Never overwrite a completed run'
    src = root / 'data/candidates/camp_primary_tissue_multicancer'
    mapping_path = src / 'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    mapping = pd.read_csv(mapping_path, dtype=str)
    cohort = mapping[mapping.Dataset.eq('COAD')].copy()
    tumor = cohort[cohort.TN.eq('Tumor')].copy()
    normal = cohort[cohort.TN.eq('Normal')].copy()
    assert len(tumor) > 0 and set(cohort.TN) == {'Tumor', 'Normal'}
    for col in ['MetabFile', 'RNAFile']:
        assert cohort[col].nunique() == 1
    assert tumor.MetabFile_sheet.nunique() == normal.MetabFile_sheet.nunique() == 1
    for col in ['CommonID', 'MetabID', 'RNAID']:
        assert not cohort[col].isna().any() and cohort[col].is_unique
    xp = src / 'processed_metabolomics' / cohort.MetabFile.iloc[0]
    rp = src / 'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed' / cohort.RNAFile.iloc[0]
    inputs = [mapping_path, xp, rp, out/'effects.tsv', out/'candidate_catalog.json', out/'identity_review.json', Path(__file__), out/'statistical_result_template.tsv']
    input_hashes = {p: sha(p) for p in inputs}
    # Read workbook once. No values or sample identifiers are exported.
    sheets = pd.read_excel(xp, sheet_name=[tumor.MetabFile_sheet.iloc[0], normal.MetabFile_sheet.iloc[0], 'data', 'metanno', 'sampleanno'])
    met = sheets[tumor.MetabFile_sheet.iloc[0]].set_index(sheets[tumor.MetabFile_sheet.iloc[0]].columns[0])
    metn = sheets[normal.MetabFile_sheet.iloc[0]].set_index(sheets[normal.MetabFile_sheet.iloc[0]].columns[0])
    observed = sheets['data'].set_index(sheets['data'].columns[0])
    annotation = sheets['metanno']
    sampleanno = sheets['sampleanno']
    rna = pd.read_csv(rp, index_col=0)
    for matrix in [met, metn, observed, rna]:
        matrix.index = matrix.index.astype(str)
        matrix.columns = matrix.columns.astype(str)
        assert matrix.index.is_unique and matrix.columns.is_unique
    assert tumor.MetabID.isin(met.columns).all() and tumor.MetabID.isin(observed.columns).all()
    assert tumor.RNAID.isin(rna.columns).all()
    assert normal.MetabID.isin(metn.columns).all() and normal.RNAID.isin(rna.columns).all()
    effects = [r for r in csv.DictReader((out/'effects.tsv').open(encoding='utf-8-sig'), delimiter='\t') if r['cancer'] == 'COAD']
    by_feature = {r['feature_name']: r for r in effects}
    assert len(by_feature) == len(effects) == 159
    candidates = json.loads((out/'candidate_catalog.json').read_text(encoding='utf-8'))
    identities = json.loads((out/'identity_review.json').read_text(encoding='utf-8'))
    assert len(candidates) == 974 and len(identities) == 73
    identity_by_feature = {r['feature_name']: r for r in identities}
    features, feature_vectors, identity_rows = [], {}, []
    for name, effect in sorted(by_feature.items()):
        found = name in met.index
        x = pd.to_numeric(met.loc[name, tumor.MetabID], errors='coerce').to_numpy(dtype=float) if found else None
        obs = pd.to_numeric(observed.loc[name, tumor.MetabID], errors='coerce').to_numpy(dtype=float) if name in observed.index else None
        feature_vectors[name] = x
        features.append(dict(feature_name=name, metabolite_key=effect['metabolite_key'],
            original_effect=effect['hedges_g'], original_q=effect['effect_fdr'],
            tumor_matrix_found=found, normal_matrix_found=name in metn.index,
            source_data_sheet_found=obs is not None,
            n_tumor_finite=int(np.isfinite(x).sum()) if x is not None else 'NA',
            n_tumor_source_data_finite=int(np.isfinite(obs).sum()) if obs is not None else 'NA',
            n_tumor_unique=int(len(np.unique(x[np.isfinite(x)]))) if x is not None else 'NA'))
        if name not in identity_by_feature:
            continue
        matches = annotation[annotation.H_name.eq(name)]
        a = matches.iloc[0] if len(matches) == 1 else None
        same = a is not None and clean(a.H_KEGG) == clean(effect['kegg_id']) and clean(a.H_HMDB) == clean(effect['hmdb_id'])
        identity_rows.append(dict(feature_name=name, metabolite_key=effect['metabolite_key'],
            source_annotation_matches=len(matches), source_kegg=clean(a.H_KEGG) if a is not None else 'NA',
            source_hmdb=clean(a.H_HMDB) if a is not None else 'NA',
            frozen_kegg=clean(effect['kegg_id']), frozen_hmdb=clean(effect['hmdb_id']),
            exact_source_ids_equal_frozen=same,
            prior_identity_status=identity_by_feature[name]['identity_review_status'],
            current_decision='SOURCE_IDS_CONFIRMED_HOLDS_UNCHANGED' if same else 'NEEDS_REVIEW_SOURCE_ID_MISMATCH',
            reason='Matching author annotations establish provenance, not stereochemistry, structural confirmation or valid database entity type.'))
    genes = sorted({(r['human_gene_id'],r['gene_symbol']) for r in candidates}, key=lambda g:(g[1],g[0]))
    gene_rows, gene_vectors = [], {}
    for gene_id,gene in genes:
        found = bool(gene) and gene in rna.index
        y = pd.to_numeric(rna.loc[gene, tumor.RNAID], errors='coerce').to_numpy(dtype=float) if found else None
        gene_vectors[gene] = y
        gene_rows.append(dict(human_gene_id=gene_id,gene=gene or 'NA', exact_symbol_found=found,
            n_tumor_finite=int(np.isfinite(y).sum()) if y is not None else 'NA',
            n_tumor_unique=int(len(np.unique(y[np.isfinite(y)]))) if y is not None else 'NA',
            status='DONE' if found else 'NOT_EVALUABLE',
            reason='Exact symbol coverage only; no RNA differential test' if found else 'GENE_SYMBOL_UNRESOLVED' if not gene else 'RNA_SYMBOL_ABSENT_NO_ALIAS_GUESSING'))
    prefix=(out/'statistical_result_template.tsv').read_text().strip().split('\t')
    extra=['human_gene_id','original_effect','original_q','n_metabolite_finite','n_rna_finite','n_metabolite_unique','n_rna_unique','biochemical_disposition','input_coverage','patient_identity','ready_for_patient_testing']
    records=[]
    for edge in candidates:
        name,gene=edge['feature_name'],edge['gene_symbol']
        effect=by_feature[name]
        x,y=feature_vectors[name],gene_vectors[gene]
        both=x is not None and y is not None
        n=int((np.isfinite(x)&np.isfinite(y)).sum()) if both else None
        coverage='GENE_SYMBOL_UNRESOLVED' if not gene else 'BOTH_PRESENT' if both else 'RNA_SYMBOL_ABSENT' if y is None else 'METABOLITE_FEATURE_ABSENT'
        row=dict.fromkeys(prefix+extra,'NA')
        row.update(cancer='COAD',cohort='COAD',stage_id='03_PATIENT',run_id=out.name,
            analysis_version='source_readiness_v2',analysis_type='source_readiness_no_association',
            metabolite_key=edge['metabolite_key'],metabolite_name=name,gene=gene or 'NA',
            unit='author_mapped_tumor_specimen',n=n if n is not None else 'NA',
            test_family='NOT_LOCKED_NO_NEW_TESTS',status='NOT_RUN' if both else 'NOT_EVALUABLE',
            reason='Association not run; patient independence and analysis family not fixed' if both else coverage,
            source_id=rp.name+';'+xp.name,human_gene_id=edge['human_gene_id'],
            original_effect=effect['hedges_g'],original_q=effect['effect_fdr'],
            n_metabolite_finite=int(np.isfinite(x).sum()) if x is not None else 'NA',
            n_rna_finite=int(np.isfinite(y).sum()) if y is not None else 'NA',
            n_metabolite_unique=int(len(np.unique(x[np.isfinite(x)]))) if x is not None else 'NA',
            n_rna_unique=int(len(np.unique(y[np.isfinite(y)]))) if y is not None else 'NA',
            biochemical_disposition=edge['disposition'],input_coverage=coverage,
            patient_identity='NOT_SEPARATELY_VERIFIED',ready_for_patient_testing=False)
        records.append(row)
    records.sort(key=lambda r: (r['cohort'],r['metabolite_name'],r['metabolite_key'],r['gene'],r['human_gene_id']))
    keys=[(r['cohort'],r['metabolite_key'],r['gene'],r['human_gene_id']) for r in records]
    assert len(keys)==len(set(keys))==974
    assert all(r['p_value']==r['q_value']==r['effect']=='NA' for r in records)
    assert sum(float(r['effect_fdr'])<.05 for r in effects)==73
    for p,h in input_hashes.items():
        assert sha(p)==h, 'Input changed during audit'
    tsv(out/'results.tsv',prefix+extra,records)
    for name,rows in [('feature_coverage.tsv',features),('gene_coverage.tsv',gene_rows),('identity_source_review.tsv',identity_rows)]:
        tsv(out/name,list(rows[0]),rows)
    sample_id_fields=[c for c in sampleanno.columns if any(k in c.lower() for k in ['patient','subject','individual','donor'])]
    summary=dict(cancer='COAD',stage_id='03_PATIENT',run_id=out.name,status='PARTIAL',
        connection='SUCCESS',source_readiness_audit='DONE',patient_association='NOT_RUN',
        n_mapping_rows=len(cohort),n_tumor_specimens=len(tumor),n_normal_specimens=len(normal),
        mapping_missing_ids={c:int(cohort[c].isna().sum()) for c in ['CommonID','MetabID','RNAID']},
        mapping_duplicate_excess={c:int(cohort[c].duplicated().sum()) for c in ['CommonID','MetabID','RNAID']},
        n_tumor_rna_ids_found=int(tumor.RNAID.isin(rna.columns).sum()),
        n_tumor_metab_ids_found=int(tumor.MetabID.isin(met.columns).sum()),
        matrix_shapes={'tumor_metabolomics':list(met.shape),'normal_metabolomics':list(metn.shape),'source_data_sheet':list(observed.shape),'rna':list(rna.shape)},
        all_matrix_rows_columns_unique=True,
        sample_annotation={'rows':len(sampleanno),'nonmissing_sample_names':int(sampleanno.SAMPLE_NAME.notna().sum()),'group_counts':{clean(k):int(v) for k,v in sampleanno.GROUP.value_counts(dropna=False).items()},'explicit_patient_id_fields':sample_id_fields},
        patient_independence='NOT_SEPARATELY_VERIFIED',pairing='NOT_INFERRED_FROM_IDS',
        n_frozen_features_found=sum(r['tumor_matrix_found'] for r in features),
        n_significant_source_id_exact_matches=sum(r['exact_source_ids_equal_frozen'] for r in identity_rows),
        n_significant_features=73,n_candidate_pairs=len(records),n_candidate_genes=len(genes),
        n_candidate_genes_found=sum(r['exact_symbol_found'] for r in gene_rows),
        missing_gene_symbols=[r['gene'] for r in gene_rows if not r['exact_symbol_found'] and r['gene']!='NA'],
        unresolved_gene_symbol_ids=[r['human_gene_id'] for r in gene_rows if r['gene']=='NA'],
        pair_input_coverage_counts=dict(Counter(r['input_coverage'] for r in records)),
        pairwise_finite_n_distribution={str(k):v for k,v in Counter(r['n'] for r in records).items()},
        supported_no_identity_hold_input_coverage=dict(Counter(r['input_coverage'] for r in records if r['biochemical_disposition']=='REACTION_ANNOTATION_SUPPORTED_PROVISIONAL')))
    dump(out/'summary.json',summary)
    dump(out/'analysis_spec.json',dict(cancer='COAD',stage_id='03_PATIENT',run_id=out.name,
        analysis_version='source_readiness_v2',code_commit=args.code_commit,generator_sha256=sha(Path(__file__)),
        created_utc=datetime.now(timezone.utc).isoformat(),status='PARTIAL',
        parameters={'mapping_dataset':'COAD','association_subset':'Tumor','new_transform':False,'imputation':False,'alias_matching':False,'identity_holds_changed':False},
        unit='author_mapped_tumor_specimen; independent patients not established',test_family='NOT_LOCKED_NO_NEW_TESTS',seed='NA',
        software={'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'openpyxl':openpyxl.__version__},
        limitations=['NA is unavailable or not applicable, never zero.','Source data sheet is author-provided; raw-detection meaning not established.','RNA filename states log2_transformed; this audit does not establish all preprocessing choices.','No original patient values or complete sample mapping exported.'],
        revision_note='v2 separates one unresolved gene ID (three pairs) from absent RNA symbols. v1 server audit retained; no statistics or frozen source changed.'))
    tsv(out/'source_manifest.tsv',['source_path','sha256'],[dict(source_path=str(p.relative_to(root)),sha256=h) for p,h in input_hashes.items()])
    dump(out/'validation.json',dict(status='PASS',checked_utc=datetime.now(timezone.utc).isoformat(),
        checks=['exclusive B run directory','all source hashes unchanged during audit','all mapping IDs nonmissing and unique','all mapped tumor/normal IDs found','matrix row and column uniqueness','159 frozen features and 974 unique candidate pairs','73 original significant features','no new effects or P/q','original effect/q strings carried through'],
        output_sha256={p.name:sha(p) for p in out.iterdir() if p.name in ['results.tsv','feature_coverage.tsv','gene_coverage.tsv','identity_source_review.tsv','summary.json','analysis_spec.json','source_manifest.tsv']},
        not_validated=['independent patients','tumor-normal pairing','raw chemical identification','raw detection mask','patient association','functional causality']))
    (out/'.running').unlink()
    (out/'DONE').write_text(datetime.now(timezone.utc).isoformat()+'\n')
    print(json.dumps(summary,ensure_ascii=False,allow_nan=False),flush=True)


if __name__=='__main__':
    main()
