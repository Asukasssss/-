"""ccRCC4 server-only author-case audit and paired discovery. Never overwrites history."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '2'
import argparse
import hashlib
import json
import platform
import tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
SRC = ROOT / 'data/candidates/camp_primary_tissue_multicancer'
V = 'ccrcc_stats_v1'
PREFIX = 'cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def save(d, p):
    d.to_csv(p, sep='\t', index=False, na_rep='NA', lineterminator='\n')

def js(p, d):
    p.write_text(json.dumps(d, ensure_ascii=False, indent=2, default=str), encoding='utf-8')

def manifest(paths, out):
    save(pd.DataFrame([dict(source_id=p.name, path_or_url=str(p), sha256=sha(p), access_scope='SERVER_PRIVATE_INPUT' if p.suffix in ['.xlsx', '.csv'] else 'PROVENANCE', version='existing_input_no_modification') for p in paths]), out/'source_manifest.tsv')

def checksums(out):
    save(pd.DataFrame([dict(file=p.name, sha256=sha(p)) for p in sorted(out.iterdir()) if p.is_file() and p.name != 'checksums.tsv']), out/'checksums.tsv')

def base(out, stage, family, key='NA', name='NA'):
    d = dict.fromkeys(PREFIX, np.nan)
    d.update(cancer='ccRCC4', cohort='CAMP_ccRCC4', stage_id=stage, run_id=out.name,
             analysis_version=V, analysis_type=family, metabolite_key=key,
             metabolite_name=name, gene='NA', unit='author_case', test_family=family,
             status='NOT_EVALUABLE', reason='NA', source_id='CAMP_v0.3.4_ccRCC4')
    return d

def adjust(d):
    for f, ix in d.groupby('test_family').groups.items():
        ok = d.index.isin(ix) & d.p_value.notna()
        p = d.loc[ok, 'p_value'].to_numpy(float)
        d.loc[ix, 'family_n_planned'] = len(ix)
        d.loc[ix, 'family_n_evaluable'] = len(p)
        if len(p):
            q = multipletests(p, method='fdr_bh')[1]
            order = np.argsort(p)
            independent = np.minimum(1, np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
            assert np.allclose(q[order], independent, atol=1e-14)
            d.loc[ok, 'q_value'] = q
    d['evidence_label'] = np.select([d.p_value.isna(), d.q_value.lt(.05), d.p_value.lt(.05)], ['NOT_EVALUABLE', 'FDR_SUPPORTED', 'NOMINAL_EXPLORATORY'], default='NOT_SUPPORTED_THIS_TEST')
    return d[PREFIX + [c for c in d if c not in PREFIX]]

def seed(family, key):
    return int(hashlib.sha256((V+'|20260925|ccRCC4|'+family+'|'+key).encode()).hexdigest()[:8], 16)

def paired_stats(delta, family, key):
    n = len(delta)
    s = seed(family, key)
    rng = np.random.default_rng(s)
    nz = delta[delta != 0]
    ranks = stats.rankdata(abs(nz))
    wp, wm = ranks[nz > 0].sum(), ranks[nz < 0].sum()
    d = dict(n=n, n_reference=n, n_pairs_used=n, n_nonzero=len(nz), n_up=int((delta>0).sum()), n_down=int((delta<0).sum()), n_equal=int((delta==0).sum()),
             up_fraction=float((delta>0).mean()) if n else np.nan, down_fraction=float((delta<0).mean()) if n else np.nan,
             mean_delta=float(delta.mean()) if n else np.nan, median_delta=float(np.median(delta)) if n else np.nan,
             rank_biserial=float((wp-wm)/(wp+wm)) if len(nz) else 0.0 if n else np.nan,
             seed=s, zero_method='wilcox', tie_method='average_ranks', continuity_correction=False, B=0, bootstrap_valid=0,
             ci_estimand='mean_paired_tumor_minus_normal', ci_method='4000_whole_case_percentile_bootstrap_pointwise')
    if n < 8:
        d.update(status='NOT_EVALUABLE', reason='fewer_than_8_complete_pairs')
        return d
    obs = abs(wp-wm)
    if not len(nz):
        p, method = 1., 'all_zero_project_convention'
    elif len(nz) <= 16:
        bits = (np.arange(2**len(nz))[:,None] >> np.arange(len(nz))) & 1
        p = float(np.mean(abs((bits*2-1) @ ranks) >= obs-1e-12)); method='exhaustive_sign_flip'
        d['B'] = 2**len(nz)
    elif len(nz) < 30:
        b = 99999; extreme = 0
        for k in range(0,b,5000):
            signs = rng.integers(0,2,size=(min(5000,b-k),len(nz)))*2-1
            extreme += int((abs(signs @ ranks) >= obs-1e-12).sum())
        p = (extreme+1)/(b+1); method='Monte_Carlo_sign_flip_plus1'; d['B']=b
    else:
        p = float(stats.wilcoxon(delta,zero_method='wilcox',correction=True,alternative='two-sided',method='approx').pvalue)
        method='normal_tie_corrected_continuity_corrected'; d['continuity_correction']=True
        # Independent signed-rank normal formula, distinct from scipy implementation.
        _, counts = np.unique(abs(nz),return_counts=True)
        nn=len(nz); var=(nn*(nn+1)*(2*nn+1)-.5*np.sum(counts**3-counts))/24
        z=(abs(wp-nn*(nn+1)/4)-.5)/np.sqrt(var)
        independent_p=2*stats.norm.sf(max(0,z))
        assert abs(p-independent_p)<1e-10, (key,p,independent_p)
    bs = delta[rng.integers(n,size=(4000,n))].mean(axis=1)
    d.update(effect=float(delta.mean()),effect_type='mean_paired_difference_author_log2_scale',
             ci_lower=float(np.quantile(bs,.025)),ci_upper=float(np.quantile(bs,.975)),p_value=p,p_method=method,
             status='DONE',reason='NA',bootstrap_valid=4000)
    return d


