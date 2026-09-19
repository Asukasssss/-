"""Versioned BRCA 117-gene screen. Read patient matrices on server165 only.

Does not run or modify the attached package's executable code or frozen CAMP tests.
"""
from pathlib import Path
import json
import hashlib
import platform
import datetime
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from camp_per_cancer_repro import calculate, bh, sha, STAT_COLUMNS

ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
OUT = ROOT / 'results/BRCA_117_screen_20260911'
SRC = ROOT / 'data/candidates/camp_primary_tissue_multicancer'


def seed(text):
    return int(hashlib.sha256(text.encode()).hexdigest()[:8], 16)


def hedges(a, b):
    na, nb = a.shape[-1], b.shape[-1]
    sp = np.sqrt(((na-1)*a.var(axis=-1, ddof=1)+(nb-1)*b.var(axis=-1, ddof=1))/(na+nb-2))
    return (1-3/(4*(na+nb)-9))*(a.mean(axis=-1)-b.mean(axis=-1))/sp


def expression(a, b, random_seed):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    base = dict(n_tumor=len(a), n_normal=len(b))
    if min(len(a), len(b)) < 8:
        return dict(base, status='SKIPPED_N_LT_8')
    if np.ptp(np.r_[a, b]) == 0:
        return dict(base, status='SKIPPED_CONSTANT_VALUES')
    test = stats.mannwhitneyu(a, b, alternative='two-sided', method='asymptotic')
    rng = np.random.default_rng(random_seed)
    aa = a[rng.integers(0, len(a), (4000, len(a)))]
    bb = b[rng.integers(0, len(b), (4000, len(b)))]
    with np.errstate(divide='ignore', invalid='ignore'):
        point = hedges(a, b)
        boot = hedges(aa, bb)
    boot = boot[np.isfinite(boot)]
    lo, hi = np.quantile(boot, [.025, .975]) if len(boot) else (np.nan, np.nan)
    return dict(base, status='COMPUTED_EXPLORATORY_UNPAIRED',
                tumor_median=np.median(a), normal_median=np.median(b),
                median_difference_author_scale=np.median(a)-np.median(b),
                rna_hedges_g=point, ci_lower=lo, ci_upper=hi,
                rank_biserial=2*test.statistic/(len(a)*len(b))-1,
                p_value=test.pvalue, bootstrap_valid=len(boot))


def add_q(df, family):
    df['q_value'] = np.nan
    good = df.p_value.notna() if 'p_value' in df else pd.Series(False, index=df.index)
    if good.any():
        df.loc[good, 'q_value'] = bh(df.loc[good, 'p_value'].to_numpy())
    df['test_family'] = family
    df['family_n_planned'] = len(df)
    df['family_n_evaluable'] = int(good.sum())
    return df


def main():
    if not ROOT.exists():
        raise RuntimeError('Patient matrices must remain on server165')
    tab = OUT / 'tables'
    if tab.exists() and any(tab.iterdir()):
        raise RuntimeError('Refusing to overwrite an existing analysis')
    tab.mkdir(exist_ok=True)
    ip = OUT / 'inputs'
    edges = pd.read_csv(ip/'BRCA_direct_human_gene_edges.tsv', sep='\t')
    sig = pd.read_csv(ip/'BRCA_significant_186_with_mapping.tsv', sep='\t')
    genes = pd.read_csv(ip/'BRCA_deduplicated_direct_gene_pool.tsv', sep='\t')
    assert len(edges)==174 and len(genes)==117 and len(sig)==186
    assert not edges.duplicated(['metabolite_key','gene_symbol']).any()
    assert sig.metabolite_key.is_unique and genes.gene_symbol.is_unique
    assert set(edges.gene_symbol)==set(genes.gene_symbol)
    spec = dict(version='BRCA117_v1_20260911', locked_edges=174, locked_genes=117,
        mapping='input curated biochemical candidates, not independently fully revalidated',
        association='tumor-only Spearman; 9999 rank permutations plus one; 4000 specimen bootstrap',
        expression='unpaired tumor vs normal Mann-Whitney U asymptotic tie correction; Hedges g with 4000 independent group bootstrap',
        expression_scale='author processed microarray, no new transform; g is not fold change',
        multiplicity='BH separately for all evaluable of 174 primary edges, 174 sensitivity edges, and 117 expression genes',
        minimum_n=8, primary_missing_handling='pairwise finite author processed values; no new imputation',
        sensitivity='author data-sheet nonmissing availability only; not a verified detection mask',
        statistical_unit='author-mapped specimens; unique patient ID not separately provided',
        sample_design='unpaired exploratory RNA contrast, no patient pairing inferred; patient independence not independently verified',
        new_gene_additions='separate mapping version; not included in these test families',
        covariates='first-pass unadjusted; subtype, batch and tissue composition not controlled',
        frozen_CAMP_statistics_recomputed=False,
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (OUT/'analysis_spec.json').write_text(json.dumps(spec, indent=2))
    mp_path=SRC/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    mp=pd.read_csv(mp_path,dtype=str)
    mp=mp[mp.Dataset.eq('BRCA1') & mp.TN.isin(['Tumor','Normal'])].copy()
    assert not mp[['CommonID','RNAID','MetabID','TN']].isna().any().any()
    assert all(mp[k].is_unique for k in ['CommonID','RNAID','MetabID'])
    assert mp.RNAFile.nunique()==1 and mp.MetabFile.nunique()==1
    rp=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/mp.RNAFile.iloc[0]
    xp=SRC/'processed_metabolomics'/mp.MetabFile.iloc[0]
    rna=pd.read_csv(rp,index_col=0)
    tm=mp[mp.TN.eq('Tumor')]; nm=mp[mp.TN.eq('Normal')]
    assert tm.MetabFile_sheet.nunique()==1
    met=pd.read_excel(xp,sheet_name=tm.MetabFile_sheet.iloc[0],index_col=0)
    obs=pd.read_excel(xp,sheet_name='data',index_col=0)
    anno=pd.read_excel(xp,sheet_name='metanno')
    for frame in [rna,met,obs]:
        frame.columns=frame.columns.astype(str)
        assert frame.columns.is_unique
    assert mp.RNAID.isin(rna.columns).all()
    assert tm.MetabID.isin(met.columns).all() and tm.MetabID.isin(obs.columns).all()
    source_paths=[mp_path,rp,xp,ip/'BRCA_direct_human_gene_edges.tsv',ip/'BRCA_significant_186_with_mapping.tsv',ip/'BRCA_deduplicated_direct_gene_pool.tsv',Path(__file__),Path(__file__).parent/'camp_per_cancer_repro.py']
    hashes={str(p):sha(p) for p in source_paths}
    pd.DataFrame([dict(path=p,sha256=h) for p,h in hashes.items()]).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
    mp.to_csv(OUT/'author_sample_mapping.tsv',sep='\t',index=False)
    old_path=ROOT/'results/CAMP_first_analysis_20260910/tables/CAMP_direct_metabolite_gene_associations.tsv'
    old_prov=ROOT/'results/CAMP_first_analysis_20260910/tables/CAMP_direct_gene_source_provenance.tsv'
    old=pd.read_csv(old_path,sep='\t')
    old_hash=pd.read_csv(old_prov,sep='\t').set_index('path').sha256.to_dict()
    reuse_ok=all(old_hash.get(str(p))==hashes[str(p)] for p in [mp_path,rp,xp])
    results=[]; idchecks=[]
    sm=sig.set_index('metabolite_key')
    for _,edge in edges.iterrows():
        key,gene=edge.metabolite_key,edge.gene_symbol
        assert key in sm.index
        sr=sm.loc[key]; name=sr.feature_name
        reason=''
        ar=anno[anno.H_name.eq(name)]
        if len(ar)!=1 or list(met.index).count(name)!=1 or list(obs.index).count(name)!=1:
            reason='NON_UNIQUE_OR_ABSENT_EXACT_AUTHOR_FEATURE'
        elif name!=edge.metabolite_name:
            reason='EDGE_SOURCE_FEATURE_NAME_CONFLICT'
        elif key.startswith('KEGG:') and str(ar.iloc[0].H_KEGG)!=key[5:]:
            reason='AUTHOR_KEGG_CONFLICT'
        elif key.startswith('HMDB:') and str(ar.iloc[0].H_HMDB)!=key[5:]:
            reason='AUTHOR_HMDB_CONFLICT'
        elif not key.startswith(('KEGG:','HMDB:','NAME:')):
            reason='UNSUPPORTED_KEY_TYPE'
        if str(sr.duplicate_metabolite_key_within_dataset)!='No':
            reason='FROZEN_DUPLICATE_FEATURE_KEY'
        if not reason and list(rna.index).count(gene)!=1:
            reason='NON_UNIQUE_OR_ABSENT_RNA_GENE'
        rid=key+'|'+gene
        idchecks.append(dict(relation_id=rid,reason=reason or 'EXACT_AUTHOR_FEATURE_AND_KEY',metabolite_key=key,gene=gene))
        base=dict(relation_id=rid,metabolite_key=key,metabolite_name=name,gene=gene,
                  source_relation_id=edge.relation_id,source_url=edge.source_url,
                  camp_g=sr.hedges_g,camp_q=sr.effect_fdr,
                  max_input_raw_missing_rate=sr.max_input_raw_missing_rate,
                  mapped_tumor_specimens=len(tm),n_unit='author_mapped_specimens',
                  independent_patient_id='NOT_SEPARATELY_PROVIDED',
                  identity_status='INPUT_IDENTITY_NOT_REIDENTIFIED',
                  interpretation='within CAMP exploratory association; not external replication or mediation')
        if not reason:
            x=pd.to_numeric(met.loc[name,tm.MetabID],errors='coerce').to_numpy(float)
            y=pd.to_numeric(rna.loc[gene,tm.RNAID],errors='coerce').to_numpy(float)
            avail=np.isfinite(pd.to_numeric(obs.loc[name,tm.MetabID],errors='coerce').to_numpy(float))
            base['n_author_data_available']=int(avail.sum())
        for kind in ['author_processed_primary','author_data_available_sensitivity']:
            sd=seed('BRCA1'+rid+kind)
            if reason:
                res=dict(n=0,status='SKIPPED_'+reason,p_value=np.nan)
                origin='NOT_COMPUTED'
            else:
                prior=old[old.cohort.eq('BRCA1') & old.relation_id.eq(rid) & old.analysis_type.eq(kind)]
                if reuse_ok and len(prior)==1:
                    res={c:prior.iloc[0][c] for c in STAT_COLUMNS}
                    origin='REUSED_IDENTICAL_INPUT_STATISTICS_NEW_Q_FAMILY'
                else:
                    xx=x if kind=='author_processed_primary' else np.where(avail,x,np.nan)
                    res=calculate(xx,y,sd); origin='NEW_CALCULATION'
            results.append(dict(base,analysis_type=kind,seed=sd,origin=origin,**res))
    assoc=pd.DataFrame(results)
    assoc=pd.concat([add_q(v.copy(),'BRCA117_v1|'+k) for k,v in assoc.groupby('analysis_type')],ignore_index=True)
    assoc.to_csv(tab/'BRCA_174_patient_associations.tsv',sep='\t',index=False)
    pd.DataFrame(idchecks).to_csv(tab/'BRCA_local_identity_and_gene_checks.tsv',sep='\t',index=False)
    print('ASSOCIATIONS_SAVED',len(assoc),flush=True)
    erows=[]
    for gene in genes.gene_symbol:
        sd=seed('BRCA117_v1|RNA_TumorMinusNormal|'+gene)
        if list(rna.index).count(gene)!=1:
            res=dict(status='SKIPPED_NON_UNIQUE_OR_ABSENT_RNA_GENE',n_tumor=0,n_normal=0,p_value=np.nan)
        else:
            a=pd.to_numeric(rna.loc[gene,tm.RNAID],errors='coerce').to_numpy(float)
            b=pd.to_numeric(rna.loc[gene,nm.RNAID],errors='coerce').to_numpy(float)
            res=expression(a,b,sd)
        erows.append(dict(gene=gene,seed=sd,design='UNPAIRED_AUTHOR_MAPPED_SPECIMENS',**res))
    expr=add_q(pd.DataFrame(erows),'BRCA117_v1|RNA_TumorMinusNormal')
    expr.to_csv(tab/'BRCA_117_RNA_differences.tsv',sep='\t',index=False)
    primary=assoc[assoc.analysis_type.eq('author_processed_primary')]
    summary=dict(spec,n_tumor_mapped=len(tm),n_normal_mapped=len(nm),
        association_planned=len(primary),association_evaluable=int(primary.p_value.notna().sum()),
        association_q_lt_005=int(primary.q_value.lt(.05).sum()),
        expression_planned=len(expr),expression_evaluable=int(expr.p_value.notna().sum()),
        expression_q_lt_005=int(expr.q_value.lt(.05).sum()),
        reused_association_rows=int(assoc.origin.str.startswith('REUSED').sum()),
        source_hashes_unchanged=all(sha(Path(p))==h for p,h in hashes.items()),
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
    assert summary['source_hashes_unchanged']
    (OUT/'patient_analysis_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
