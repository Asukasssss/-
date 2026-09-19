"""CAMP per-cancer reproduction. Source matrices are read only on server165.

export: split existing evidence, without recomputing statistics.
preflight: inspect exact author mapping and source matrix coverage.
associate: reuse existing direct-edge statistics; --recompute verifies reproduction.
depmap: describe supplied processed gene effects using explicitly reviewed model IDs.
"""
import argparse
import datetime
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

CANCERS = ['BRCA', 'COAD', 'GBM', 'PDAC', 'PRAD', 'ccRCC']
COHORTS = {'BRCA': ['BRCA1'], 'COAD': ['COAD'], 'GBM': ['GBM'],
           'PDAC': ['PDAC'], 'PRAD': ['PRAD'], 'ccRCC': ['ccRCC3', 'ccRCC4']}
ORIGINAL = 'results/CAMP_first_analysis_20260910'
SOURCE = 'data/candidates/camp_primary_tissue_multicancer'
SERVER_ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
ASSOC_FILE = 'CAMP_direct_metabolite_gene_associations.tsv'
STAT_COLUMNS = ['n', 'rho', 'ci_lower', 'ci_upper', 'p_value',
                'p_asymptotic_descriptive', 'bootstrap_valid', 'status']
REL_REQUIRED = ['cancer', 'metabolite_key', 'canonical_metabolite', 'kegg_id',
                'gene', 'relation_id', 'relation_type', 'human_relation_evidence',
                'source_url', 'compartment']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding='utf-8')


def fresh(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError('Output is not empty; use a NEW output directory: ' + str(path))
    path.mkdir(parents=True, exist_ok=True)
    return path


def table(df, path):
    df.to_csv(path, sep='\t', index=False)


def manifest(paths, out):
    rows = [{'path': str(p), 'bytes': Path(p).stat().st_size, 'sha256': sha(p)}
            for p in dict.fromkeys(paths)]
    table(pd.DataFrame(rows), out / 'input_manifest.tsv')
    write_json(out / 'run_environment.json', {
        'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'python': platform.python_version(), 'numpy': np.__version__,
        'pandas': pd.__version__, 'scipy': scipy.__version__, 'script_sha256': sha(__file__),
        'frozen_CAMP_effects_recomputed': False})


def selected(cancer):
    return CANCERS if cancer == 'ALL' else [cancer]


def server_only(root):
    if Path(root).resolve() != SERVER_ROOT:
        raise ValueError('Matrix operations must run under the verified server165 project root.')


def locate(root, relative):
    p = root / ORIGINAL / relative
    if not p.exists():
        raise FileNotFoundError(str(p))
    return p


def export(a):
    out = fresh(a.out)
    paths = [locate(a.project_root, 'tables/disease_association_comparison.tsv'),
             locate(a.project_root, 'tables/' + ASSOC_FILE),
             locate(a.project_root, 'tables/reused_gene_expression_evidence.tsv'),
             locate(a.project_root, 'parameters/CAMP_locked_direct_relations.tsv')]
    bp = a.project_root / ORIGINAL / 'inputs/first_analysis_batch.tsv'
    if not bp.exists():
        bp = locate(a.project_root, 'parameters/first_analysis_batch.tsv')
    paths.append(bp)
    disease, association, expression, relations, batch = [pd.read_csv(p, sep='\t') for p in paths]
    if disease.duplicated(['cancer', 'metabolite_key']).any():
        raise ValueError('Duplicate frozen cancer-metabolite keys')
    summaries = []
    for cancer in selected(a.cancer):
        dest = fresh(out / cancer)
        d, w = disease[disease.cancer.eq(cancer)], batch[batch.cancer.eq(cancer)]
        r = relations[relations.cancer.eq(cancer)]
        v = association[association.cancer.eq(cancer)]
        e = expression[expression.cancer.eq(cancer)]
        for name, frame in [('all_frozen_effects.tsv', d), ('work_batch.tsv', w),
                            ('locked_relations.tsv', r), ('reused_associations.tsv', v),
                            ('reused_expression.tsv', e)]:
            table(frame, dest / name)
        table(r[['cancer', 'gene']].drop_duplicates(), dest / 'gene_cancer_list.tsv')
        if cancer == 'BRCA':
            ep = locate(a.project_root, 'tables/BRCA_FUSCC_metabolite_summary.tsv')
            xp = locate(a.project_root, 'tables/BRCA_FUSCC_evidence_long.tsv')
            paths.extend([ep, xp])
            external = pd.read_csv(ep, sep='\t')
            table(external, dest / 'external_all318_summary.tsv')
            table(pd.read_csv(xp, sep='\t'), dest / 'external_all318_long.tsv')
            ready = external.external_evidence_status.eq('SAME_DIRECTION_Q_SUPPORTED_CONDITIONAL') & external.camp_original_q.lt(.05) & external.first_analysis_work_batch.astype(str).str.lower().eq('true')
            table(external[ready], dest / 'external_supported_work22.tsv')
        p = v[v.analysis_type.eq('author_processed_primary')]
        summaries.append({'cancer': cancer, 'cohorts': ';'.join(COHORTS[cancer]),
                          'all_effects': len(d), 'work_batch': len(w),
                          'locked_edges': len(r), 'primary_tests': len(p),
                          'primary_q_lt_005': int(p.q_value.lt(.05).sum()),
                          'reused_expression_rows': len(e),
                          'association_status': 'NOT_IN_ORIGINAL_PANEL' if cancer == 'GBM' else
                          'MULTIREGION_PATIENT_UNIT_UNVERIFIED' if cancer == 'ccRCC' else 'EXISTING_RESULTS_REUSED'})
    table(pd.DataFrame(summaries), out / 'cancer_summary.tsv')
    manifest(paths, out)
    print(pd.DataFrame(summaries).to_string(index=False))


def matrix_paths(root, mp):
    src = root / SOURCE
    if not all(mp[c].nunique() == 1 for c in ['RNAFile', 'MetabFile', 'MetabFile_sheet']):
        raise ValueError('Author mapping refers to multiple files or sheets within tumor cohort')
    rp = src / 'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed' / mp.RNAFile.iloc[0]
    xp = src / 'processed_metabolomics' / mp.MetabFile.iloc[0]
    return rp, xp


def mapping_ok(mp):
    return len(mp) > 0 and not mp[['CommonID', 'MetabID', 'RNAID']].isna().any().any() and all(mp[c].is_unique for c in ['CommonID', 'MetabID', 'RNAID'])


def preflight(a):
    server_only(a.project_root)
    out = fresh(a.out)
    mapping_path = a.project_root / SOURCE / 'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    mapping = pd.read_csv(mapping_path, dtype=str)
    paths, rows = [mapping_path], []
    for cancer in selected(a.cancer):
        for cohort in COHORTS[cancer]:
            m = mapping[mapping.Dataset.eq(cohort)]
            t = m[m.TN.eq('Tumor')]
            row = {'cancer': cancer, 'cohort': cohort, 'n_tumor_specimens': len(t),
                   'n_normal_specimens': int(m.TN.eq('Normal').sum()),
                   'exact_mapping_unique': mapping_ok(t), 'patient_identity': 'NOT_SEPARATELY_VERIFIED',
                   'genotype_check': 'UNAVAILABLE_AUTHOR_PROCESSED_DATA_ONLY'}
            try:
                rp, xp = matrix_paths(a.project_root, t)
                rna = pd.read_csv(rp, index_col=0)
                met = pd.read_excel(xp, sheet_name=t.MetabFile_sheet.iloc[0], index_col=0)
                row.update({'rna_path': str(rp), 'metabolomics_path': str(xp),
                            'tumor_sheet': t.MetabFile_sheet.iloc[0],
                            'n_tumor_rna_ids_found': int(t.RNAID.isin(rna.columns.astype(str)).sum()),
                            'n_tumor_metab_ids_found': int(t.MetabID.isin(met.columns.astype(str)).sum()),
                            'rna_columns_unique': rna.columns.is_unique,
                            'metab_columns_unique': met.columns.is_unique})
                row['status'] = ('MULTIREGION_PATIENT_UNIT_UNVERIFIED' if cancer == 'ccRCC' else 'AUTHOR_MAPPED_SPECIMEN_DATA_AVAILABLE') if mapping_ok(t) and row['n_tumor_rna_ids_found'] == len(t) and row['n_tumor_metab_ids_found'] == len(t) and rna.columns.is_unique and met.columns.is_unique else 'JOIN_REQUIRES_REVIEW'
                paths.extend([rp, xp])
            except (OSError, ValueError, KeyError, IndexError) as exc:
                row['status'], row['detail'] = 'SOURCE_CHECK_FAILED', str(exc)
            rows.append(row)
            print(cohort, row['status'], flush=True)
    table(pd.DataFrame(rows), out / 'cohort_inputs.tsv')
    manifest(paths, out)


def calculate(x, y, seed):
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    n = len(x)
    blank = {k: np.nan for k in STAT_COLUMNS}
    if n < 8:
        return dict(blank, n=n, status='SKIPPED_N_LT_8')
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return dict(blank, n=n, status='SKIPPED_CONSTANT_VALUES')
    rho, p_asym = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    den = np.sqrt((rx @ rx) * (ry @ ry))
    perms = rng.random((9999, n)).argsort(axis=1)
    null = (ry[perms] * rx).sum(axis=1) / den
    p = (1 + (np.abs(null) >= abs(rho) - 1e-14).sum()) / 10000
    inds = rng.integers(0, n, size=(4000, n))
    bx, by = stats.rankdata(x[inds], axis=1), stats.rankdata(y[inds], axis=1)
    bx -= bx.mean(axis=1, keepdims=True)
    by -= by.mean(axis=1, keepdims=True)
    denb = np.sqrt((bx * bx).sum(axis=1) * (by * by).sum(axis=1))
    valid = denb > 0
    boot = (bx[valid] * by[valid]).sum(axis=1) / denb[valid]
    lo, hi = np.quantile(boot, [.025, .975]) if len(boot) else (np.nan, np.nan)
    return dict(n=n, rho=rho, ci_lower=lo, ci_upper=hi, p_value=p,
                p_asymptotic_descriptive=p_asym, bootstrap_valid=len(boot), status='COMPUTED_EXPLORATORY')


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    qsort = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    answer = np.empty(len(values))
    answer[order] = np.clip(qsort, 0, 1)
    return answer


def associate(a):
    server_only(a.project_root)
    if a.check_only and a.recompute:
        raise ValueError('--check-only and --recompute are mutually exclusive.')
    if a.relations and a.family_id == 'locked_direct_panel_v1':
        raise ValueError('An expanded panel needs a NEW --family-id.')
    out = fresh(a.out)
    relation_path = a.relations or locate(a.project_root, 'parameters/CAMP_locked_direct_relations.tsv')
    relations = pd.read_csv(relation_path, sep='\t')
    for col in REL_REQUIRED:
        if col not in relations or relations[col].isna().any():
            raise ValueError('Missing direct-relation evidence column/value: ' + col)
    if relations.duplicated(['cancer', 'relation_id']).any():
        raise ValueError('Duplicate cancer-relation keys')
    if not relations.relation_id.eq(relations.metabolite_key + '|' + relations.gene).all():
        raise ValueError('relation_id must be metabolite_key|gene')
    if not relations.metabolite_key.eq('KEGG:' + relations.kegg_id).all():
        raise ValueError('This implementation requires exact KEGG mapping; no inferred HMDB conversion.')
    if not relations.cancer.isin(CANCERS).all():
        raise ValueError('Unknown cancer code')
    relations = relations[relations.cancer.isin(selected(a.cancer))].copy()
    table(relations, out / 'locked_relations.tsv')  # lock before reading outcomes
    mapping_path = a.project_root / SOURCE / 'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    effect_path = locate(a.project_root, 'parameters/cohort_level_effects.tsv')
    old_path = locate(a.project_root, 'tables/' + ASSOC_FILE)
    provenance_path = locate(a.project_root, 'tables/CAMP_direct_gene_source_provenance.tsv')
    mapping = pd.read_csv(mapping_path, dtype=str)
    effects = pd.read_csv(effect_path, sep='\t')
    old = pd.read_csv(old_path, sep='\t')
    old_prov = pd.read_csv(provenance_path, sep='\t').set_index('path').sha256.to_dict()
    paths = [relation_path, mapping_path, effect_path, old_path, provenance_path]
    rows, exclusions = [], []
    for cancer in selected(a.cancer):
        group = relations[relations.cancer.eq(cancer)]
        if group.empty:
            exclusions.append(dict(cancer=cancer, reason='NO_LOCKED_RELATIONS_FOR_THIS_CANCER'))
            continue
        if cancer == 'ccRCC':
            for _, r in group.iterrows():
                exclusions.append(dict(cancer=cancer, relation_id=r.relation_id, reason='MULTIREGION_PATIENT_UNIT_UNVERIFIED'))
            continue
        cohort = COHORTS[cancer][0]
        mp = mapping[mapping.Dataset.eq(cohort) & mapping.TN.eq('Tumor')].copy()
        if not mapping_ok(mp):
            exclusions.append(dict(cancer=cancer, reason='AUTHOR_MAPPING_MISSING_OR_DUPLICATE_KEYS'))
            continue
        rp, xp = matrix_paths(a.project_root, mp)
        paths.extend([rp, xp])
        rna = pd.read_csv(rp, index_col=0)
        met = pd.read_excel(xp, sheet_name=mp.MetabFile_sheet.iloc[0], index_col=0)
        obs = pd.read_excel(xp, sheet_name='data', index_col=0)
        anno = pd.read_excel(xp, sheet_name='metanno')
        for frame in [rna, met, obs]:
            frame.columns = frame.columns.astype(str)
        if not all(x.columns.is_unique for x in [rna, met, obs]) or not mp.RNAID.isin(rna.columns).all() or not mp.MetabID.isin(met.columns).all() or not mp.MetabID.isin(obs.columns).all():
            exclusions.append(dict(cancer=cancer, reason='MATRIX_COLUMNS_DUPLICATE_OR_MAPPED_IDS_ABSENT'))
            continue
        can_reuse = all(old_prov.get(str(p)) == sha(p) for p in [rp, xp, mapping_path])
        if a.recompute and not a.relations and not can_reuse:
            raise ValueError('Historical reproduction input hashes differ from original: ' + cohort)
        for _, r in group.iterrows():
            e = effects[effects.dataset.eq(cohort) & effects.metabolite_key.eq(r.metabolite_key)]
            reason = ''
            if len(e) != 1:
                reason = 'NON_UNIQUE_OR_ABSENT_FROZEN_COHORT_KEY'
            elif str(e.iloc[0].duplicate_metabolite_key_within_dataset) != 'No' or str(e.iloc[0].meta_eligible) != 'Yes':
                reason = 'FROZEN_KEY_DUPLICATE_OR_INELIGIBLE'
            else:
                feature = e.iloc[0].feature_name
                if feature.strip().lower() != r.canonical_metabolite.strip().lower() or '/' in feature or ';' in feature:
                    reason = 'CANONICAL_SOURCE_NAME_CONFLICT_OR_MIXTURE'
                elif list(met.index).count(feature) != 1 or list(obs.index).count(feature) != 1 or list(rna.index).count(r.gene) != 1:
                    reason = 'NON_UNIQUE_OR_ABSENT_MATRIX_FEATURE_OR_GENE'
                elif len(anno[anno.H_name.eq(feature) & anno.H_KEGG.eq(r.kegg_id)]) != 1:
                    reason = 'ANNOTATION_NAME_KEGG_NOT_ONE_TO_ONE'
            if reason:
                exclusions.append(dict(cancer=cancer, relation_id=r.relation_id, reason=reason))
                continue
            x = pd.to_numeric(met.loc[feature, mp.MetabID], errors='coerce').values.astype(float)
            y = pd.to_numeric(rna.loc[r.gene, mp.RNAID], errors='coerce').values.astype(float)
            available = pd.to_numeric(obs.loc[feature, mp.MetabID], errors='coerce').values.astype(float)
            for atype, xx in [('author_processed_primary', x), ('author_data_available_sensitivity', np.where(np.isfinite(available), x, np.nan))]:
                seed = int(hashlib.sha256((cohort + r.relation_id + atype).encode()).hexdigest()[:8], 16)
                match = old[old.cohort.eq(cohort) & old.relation_id.eq(r.relation_id) & old.analysis_type.eq(atype)]
                reuse = not a.recompute and can_reuse and len(match) == 1 and match.iloc[0].metabolite_name == feature
                if a.check_only:
                    result = {k: np.nan for k in STAT_COLUMNS}
                    result.update(n=int((np.isfinite(xx) & np.isfinite(y)).sum()), status='NOT_COMPUTED_PANEL_CHECK')
                else:
                    result = match.iloc[0][STAT_COLUMNS].to_dict() if reuse else calculate(xx, y, seed)
                rows.append(dict(cancer=cancer, cohort=cohort, metabolite_key=r.metabolite_key,
                                 metabolite_name=feature, gene=r.gene, relation_id=r.relation_id,
                                 relation_type=r.relation_type, source_url=r.source_url,
                                 analysis_type=atype, test_family=cohort+'|'+atype+'|'+a.family_id,
                                 n_unit='author_mapped_tumor_specimens', patient_identity='NOT_SEPARATELY_VERIFIED',
                                 evidence_origin='PANEL_CHECK_ONLY' if a.check_only else 'REUSED_IDENTICAL_INPUT_STATISTICS' if reuse else 'RECOMPUTED_OR_NEW_EXPLORATORY_ASSOCIATION',
                                 seed=seed, previous_q=match.iloc[0].q_value if len(match)==1 else np.nan,
                                 camp_hedges_g=e.iloc[0].hedges_g, camp_original_q=e.iloc[0].wilcoxon_fdr,
                                 **result))
        print(cohort, 'finished', flush=True)
    df = pd.DataFrame(rows)
    if not df.empty:
        df['q_value'] = np.nan
        df['test_family_n_evaluable'] = 0
        for family, indices in df.groupby('test_family').groups.items():
            valid = [i for i in indices if np.isfinite(df.loc[i, 'p_value'])]
            if valid:
                df.loc[valid, 'q_value'] = bh(df.loc[valid, 'p_value'])
                df.loc[indices, 'test_family_n_evaluable'] = len(valid)
    else:
        df = pd.DataFrame(columns=['cancer','cohort','relation_id','gene','analysis_type']+STAT_COLUMNS+['q_value'])
    table(df, out / 'associations.tsv')
    table(pd.DataFrame(exclusions, columns=['cancer','relation_id','reason']), out / 'exclusions.tsv')
    write_json(out / 'analysis_spec.json', {
        'family_id': a.family_id, 'relations_sha256': sha(out/'locked_relations.tsv'),
        'unit': 'author mapped tumor specimen; no independent genotype or patient uniqueness verification',
        'ccRCC': 'blocked until patient-region design is explicitly implemented',
        'method': 'Spearman; 9999 two-sided rank permutations plus one; 4000 paired-specimen percentile bootstrap',
        'minimum_n': 8, 'minimum_n_meaning': 'computational floor, not adequate-power certification',
        'seed': 'first 8 hex of SHA256(cohort+relation_id+analysis_type)',
        'BH': 'within full locked cohort and analysis type; all finite P, no direction/significance filtering',
        'missing_values': 'complete pairs; sensitivity uses author data-sheet availability, not verified detection',
        'results': len(df), 'exclusions': len(exclusions), 'recompute_requested': a.recompute,
        'check_only': a.check_only})
    if a.recompute and not a.relations and not df.empty:
        compare = df.merge(old, on=['cohort','relation_id','analysis_type'], suffixes=('_new','_old'), validate='one_to_one')
        checks = {col: bool(np.allclose(compare[col+'_new'], compare[col+'_old'], rtol=1e-10, atol=1e-12, equal_nan=True)) for col in ['n','rho','ci_lower','ci_upper','p_value','q_value']}
        checks['row_coverage'] = len(compare) == len(df) == len(old[old.cancer.isin(selected(a.cancer))])
        write_json(out / 'reproduction_check.json', checks)
        if not all(checks.values()):
            raise ValueError('Historical numeric reproduction mismatch; inspect reproduction_check.json')
    manifest(paths, out)
    print('association_rows', len(df), 'exclusions', len(exclusions))


def depmap(a):
    server_only(a.project_root)
    out = fresh(a.out)
    genes = pd.read_csv(a.genes, sep='\t')[['cancer','gene']].drop_duplicates()
    model_map = pd.read_csv(a.model_map, sep='\t')
    required = ['ModelID','cancer','review_status','source_label','source_url']
    if not set(required).issubset(model_map.columns) or model_map[required].isna().any().any():
        raise ValueError('Model map requires nonmissing '+','.join(required))
    if not model_map.ModelID.is_unique or not model_map.cancer.isin(CANCERS).all() or not model_map.review_status.eq('ACCEPT').all():
        raise ValueError('Models must be uniquely assigned and explicitly ACCEPTed; broad lineage alone is insufficient.')
    if not genes.cancer.isin(CANCERS).all() or genes.gene.isna().any():
        raise ValueError('Invalid cancer-gene list')
    models = pd.read_csv(a.model)
    if not models.ModelID.is_unique or not model_map.ModelID.isin(models.ModelID).all():
        raise ValueError('Model metadata missing or duplicate reviewed ModelIDs')
    header = pd.read_csv(a.effects, nrows=0).columns
    colmap = {}
    for gene in genes.gene.unique():
        hits = [c for c in header[1:] if c.split(' (')[0] == gene]
        if len(hits) > 1:
            raise ValueError('Ambiguous gene symbol in gene-effect file: '+gene)
        if hits:
            colmap[hits[0]] = gene
    effect = pd.read_csv(a.effects, usecols=[header[0]]+list(colmap)).rename(columns={header[0]:'ModelID', **colmap})
    if not effect.ModelID.is_unique:
        raise ValueError('Duplicate ModelID in gene effects')
    frame = model_map.merge(effect, on='ModelID', how='left', validate='one_to_one')
    essentials = set(pd.read_csv(a.essentials, sep='\t').gene) if a.essentials else None
    rows, long = [], []
    for _, g in genes.iterrows():
        sub = frame[frame.cancer.eq(g.cancer)]
        values = pd.to_numeric(sub[g.gene], errors='coerce') if g.gene in sub else pd.Series(np.nan, index=sub.index)
        v = values[np.isfinite(values)]
        rows.append(dict(cancer=g.cancer, gene=g.gene, release=a.release, n_reviewed_models=len(sub), n_gene_effect_available=len(v),
                         median=v.median(), q25=v.quantile(.25), q75=v.quantile(.75),
                         n_below_minus_0_5=int((v<-.5).sum()), fraction_below_minus_0_5=float((v<-.5).mean()) if len(v) else np.nan,
                         n_below_minus_1=int((v<-1).sum()),
                         common_essential='NOT_ASSESSED' if essentials is None else str(g.gene in essentials),
                         status='DESCRIPTIVE_SCREEN_CONTEXT' if len(v) else 'UNEVALUABLE_NO_GENE_EFFECT_OR_MODELS',
                         normal_tissue_safety='NOT_ASSESSED', metabolite_mediation='NOT_ESTABLISHED', automatic_L2_assignment=False))
        for index, model in sub.iterrows():
            long.append(dict(cancer=g.cancer, gene=g.gene, ModelID=model.ModelID,
                             source_label=model.source_label, gene_effect=values.loc[index], release=a.release))
    table(pd.DataFrame(rows), out/'dependency_summary.tsv')
    table(pd.DataFrame(long), out/'dependency_model_values.tsv')
    manifest([a.effects,a.model,a.model_map,a.genes]+([a.essentials] if a.essentials else []), out)
    print('dependency_summary_rows',len(rows))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for name in ['export','preflight','associate','depmap']:
        q = sub.add_parser(name)
        q.add_argument('--project-root', type=Path, default=SERVER_ROOT)
        q.add_argument('--out', type=Path, required=True)
        if name != 'depmap':
            q.add_argument('--cancer', choices=CANCERS+['ALL'], default='ALL')
        if name == 'associate':
            q.add_argument('--relations', type=Path)
            q.add_argument('--family-id', default='locked_direct_panel_v1')
            q.add_argument('--recompute', action='store_true')
            q.add_argument('--check-only', action='store_true', help='Check mapping and finite-pair counts; do not calculate rho/P/q/CI.')
        if name == 'depmap':
            for flag in ['effects','model','model-map','genes']:
                q.add_argument('--'+flag, type=Path, required=True)
            q.add_argument('--release', required=True)
            q.add_argument('--essentials', type=Path)
    args = p.parse_args()
    globals()[args.command](args)


if __name__ == '__main__':
    main()
