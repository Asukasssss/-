"""PRAD server-only author-case audit and paired discovery. Never overwrites history."""
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
V = 'prad_discovery_v1'
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
    d.update(cancer='PRAD', cohort='CAMP_PRAD', stage_id=stage, run_id=out.name,
             analysis_version=V, analysis_type=family, metabolite_key=key,
             metabolite_name=name, gene='NA', unit='author_case', test_family=family,
             status='NOT_EVALUABLE', reason='NA', source_id='CAMP_v0.3.4_PRAD')
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
    return int(hashlib.sha256((V+'|20260922|PRAD|'+family+'|'+key).encode()).hexdigest()[:8], 16)

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

def audit(out, commit):
    out.resolve().relative_to(ROOT/'results/collaborative/PRAD/B')
    out.mkdir(parents=True,exist_ok=True)
    assert not (out/'public').exists() and not (out/'private').exists()
    with (out/'.running').open('x') as lock: lock.write(V)
    for p in ['source','private','public/01_CAMP']: (out/p).mkdir(parents=True,exist_ok=True)
    pub=out/'public/01_CAMP'
    source=out/'source/PRAD_original.xlsx'
    if not source.exists():
        with tarfile.open(SRC/'pancancer_metabolomics_v.0.3.4.tar.gz') as t:
            with t.extractfile('pancancer_metabolomics/data/metabolomics_original/PRAD.xlsx') as f:
                source.write_bytes(f.read())
    a=pd.read_excel(source,sheet_name='sampleinfo',dtype=str)
    mp=SRC/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    m=pd.read_csv(mp,dtype=str);m=m[m.Dataset.eq('PRAD')].copy()
    assert all(m[c].is_unique for c in ['MetabID','RNAID','CommonID']) and a.SAMPLE_NAME.is_unique
    m=m.merge(a,left_on='MetabID',right_on='SAMPLE_NAME',validate='one_to_one')
    assert len(m)==137 and m['case'].notna().all()
    assert m.TN.eq(m['T.or.N.'].map({'T':'Tumor','N':'Normal'})).all()
    assert (m.groupby('case').Identifier.nunique()==1).all()
    t=m[m.TN.eq('Tumor')];n=m[m.TN.eq('Normal')]
    assert t['case'].is_unique and n['case'].is_unique
    pairs=t[['case','MetabID','RNAID','CAPT','Identifier']].merge(n[['case','MetabID','RNAID','CAPT','Identifier']],on='case',suffixes=('_tumor','_normal'),validate='one_to_one')
    pairs['CAPT_concordant']=pairs.CAPT_tumor.eq(pairs.CAPT_normal)
    assert pairs.Identifier_tumor.eq(pairs.Identifier_normal).all()
    metp=SRC/'processed_metabolomics/PreprocessedData_PRAD.xlsx'
    rp=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
    rna=pd.read_csv(rp,index_col=0);x=pd.ExcelFile(metp)
    matrix=pd.read_excel(x,sheet_name='data_imputed',index_col=0)
    raw=pd.read_excel(x,sheet_name='data',index_col=0)
    author_names=matrix.index.astype(str).tolist()
    matrix.index=matrix.index.astype(str).str.strip();raw.index=raw.index.astype(str).str.strip()
    anno=pd.read_excel(x,sheet_name='sampleanno',index_col=0)
    assert matrix.index.is_unique and matrix.columns.is_unique and rna.columns.is_unique
    assert m.MetabID.isin(matrix.columns).all() and m.RNAID.isin(rna.columns).all()
    assert m.MetabID.map(anno.GROUP).eq(m.TN.str.upper()).all()
    assert np.allclose(matrix.values[np.isfinite(raw.values)],raw.values[np.isfinite(raw.values)])
    indices=[]
    for tissue in ['Tumor','Normal']:
        d=pd.read_excel(x,sheet_name='metabo_imputed_filtered_'+tissue,index_col=0)
        d.index=d.index.astype(str).str.strip();assert d.index.is_unique
        assert np.array_equal(d.values,matrix.loc[d.index,d.columns].values)
        assert set(d.columns)==set(m.loc[m.TN.eq(tissue),'MetabID'])
        indices.append(set(d.index))
    retained=indices[0]&indices[1]
    oldp=ROOT/'results/tables/cohort_effects.tsv';old=pd.read_csv(oldp,sep='\t');old=old[old.dataset.eq('PRAD')].copy()
    assert set(old.feature_name)==retained and old.feature_name.is_unique
    frozenp=ROOT/'results/CAMP_Phase1_v1.0'
    # The released repository snapshot is compared with the actual unchanged server effects.
    effects=ROOT/'results/tables/cancer_effects.tsv';frozen=pd.read_csv(effects,sep='\t');frozen=frozen[frozen.cancer.eq('PRAD')]
    save(m,out/'private/sample_identity_audit_private.tsv');save(pairs,out/'private/metabolite_pairs_private.tsv')
    save(pairs,out/'private/RNA_pairs_private.tsv');save(t,out/'private/tumor_multiomics_map_private.tsv')
    save(old,pub/'historical_cohort_effects_preserved.tsv');save(frozen,pub/'historical_cancer_effects_preserved.tsv')
    feature=pd.DataFrame({'feature_name':matrix.index,'author_feature_name':author_names})
    feature['retained_both_author_tissue_filters']=feature.feature_name.isin(retained)
    feature['feature_id']=feature.feature_name.map(lambda z:'PRAD_FEATURE:'+hashlib.sha256(z.encode()).hexdigest()[:16])
    feature=feature.merge(old[['feature_name','metabolite_key','identity_confidence']],on='feature_name',how='left',validate='one_to_one')
    feature['chemical_key_collision']=feature.metabolite_key.notna() & feature.metabolite_key.duplicated(False)
    feature['status']=np.where(feature.retained_both_author_tissue_filters,'DONE','NOT_EVALUABLE')
    feature['reason']=np.where(feature.retained_both_author_tissue_filters,'NA','not_in_both_author_filtered_tissue_matrices')
    save(feature,pub/'feature_inventory.tsv')
    summary=dict(status='DONE',mapped_specimens=len(m),tumor_author_cases=len(t),normal_author_cases=len(n),author_cases=m['case'].nunique(),
        author_case_pairs=len(pairs),CAPT_concordant_pairs=int(pairs.CAPT_concordant.sum()),case_pairs_CAPT_discordant=int((~pairs.CAPT_concordant).sum()),
        tissue_conflicts=0,duplicate_tumor_cases=0,duplicate_normal_cases=0,case_cross_batch=0,
        total_input_features=len(matrix),retained_both_tissues=len(retained),frozen_cancer_features=len(frozen),frozen_q_lt005=int(frozen.effect_fdr.lt(.05).sum()),
        chemical_key_collision_features=int(feature.chemical_key_collision.sum()),RNA_rows=len(rna),RNA_unique_symbols=bool(rna.index.is_unique),
        numerical_results_this_audit=0,original_statistics_modified=False,genotype_identity='NOT_EVALUABLE_no_genotypes',
        patient_identity_basis='Explicit author case field; not inferred from suffix/order; CAPT discrepancy retained for conservative subset sensitivity',
        available_mask_semantics='finite_value_in_author_data_sheet; not an independently verified detection mask')
    js(pub/'validation.json',summary)
    save(pd.DataFrame([dict(audit_item=k,value=v,status='DONE' if k!='genotype_identity' else 'NOT_EVALUABLE') for k,v in summary.items()]),pub/'sample_audit_summary.tsv')
    save(m.groupby(['Identifier','TN']).size().rename('n_specimens').reset_index(),pub/'batch_tissue_counts.tsv')
    spec=dict(version=V,code_commit=commit,author_pair_field='PRAD.xlsx/sampleinfo/case',same_batch_pairs=True,
        unit='author case, unique within tissue; no genotype verification',primary_pairs=len(pairs),CAPT_concordant_sensitivity=int(pairs.CAPT_concordant.sum()),
        full_input_inventory=len(matrix),planned_retained_features=len(retained),feature_key='feature_id hash of outer-whitespace-trimmed author harmonized feature name; original spelling retained; uniqueness checked; never collapse chemical-key collisions',
        main_test='two-sided signed rank; exact sign flips n_nonzero<=16;99999 MC n_nonzero17-29 plus1;normal ties+continuity>=30;all zero p1',
        minimum_n=8,bootstrap=4000,master_seed=20260922,seed_derivation='SHA256(version|20260922|PRAD|family|feature_id) first8hex',
        families=['METAB_PAIRED_PRIMARY','METAB_PAIRED_AVAILABLE','METAB_PAIRED_CAPT_CONCORDANT'],BH='separate all evaluable retained features in each family',
        workpool='primary nominal P<0.05; FDR label separately retained',data_scale='CAMP author batch median scaling, PQN, log2, minimum imputation; no new transformations',
        observed_mask='finite author data values; unchanged observed values verified against data_imputed',
        assumptions='sign-flip exchangeability; signed-rank location interpretation assumes symmetric differences; pointwise bootstrap intervals',
        software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
        source_paper='https://doi.org/10.1038/s42255-023-00817-8',BRCA_reference_commit='f18a215dc8a269e4f16c55011633d557f069f631',
        scope='Discovery to all-candidate cell expression source; no automatic mechanism extensions')
    js(pub/'analysis_spec.json',spec)
    manifest([source,mp,metp,rp,oldp,effects,Path(__file__)],pub)
    (pub/'README_CN.md').write_text(f'''# PRAD 输入、真实配对与冻结结果审计\n\n问题：能否按BRCA主线进行PRAD癌种内分析？\n\n输入：CAMP v0.3.4原始PRAD.xlsx的sampleinfo、处理矩阵、MasterMapping与作者RNA矩阵；全部源数据留server165。\n\n实际结果：{len(m)}个标本，{len(t)}肿瘤病例、{len(n)}正常病例；作者case明确连接{len(pairs)}对。组织标签全部一致，每组织内case唯一，配对均同批次。另一个CAPT字段有{int((~pairs.CAPT_concordant).sum())}对不一致，不能隐去，预设仅CAPT一致的{int(pairs.CAPT_concordant.sum())}对敏感性。\n\n原矩阵{len(matrix)}个特征；同时通过作者两组织过滤的{len(retained)}项为新分析全族。旧癌种冻结表{len(frozen)}项，原q<0.05为{int(frozen.effect_fdr.lt(.05).sum())}。{int(feature.chemical_key_collision.sum())}项具有碰撞化学键，保留不同峰名与feature_id，不合并或遗漏。\n\n新手解释：43对指同一作者病例有两种组织，不是把两列按顺序配起来。新配对分析尚未在本审计批执行，旧非配对P/q完整保留。\n\n限制：依据作者case字段，不是基因型身份验证；CAPT差异含义未明确。data中有限值只称作者可用值，不能宣称原始检出掩码。\n\n当前决定：冻结analysis_spec，按作者case进行配对主分析，并执行可用值与CAPT一致配对敏感性。\n\n下一步：全{len(retained)}特征配对分析、P<0.05工作池、直接生化映射；相关/RNA/单细胞保留各自状态。\n\n复现：python prad_discovery_v1.py --mode audit --out NEW_SERVER_PRAD_RUN --code-commit COMMIT；然后 --mode paired --out SAME_RUN。\n''',encoding='utf-8')
    checksums(pub);print(json.dumps(summary),flush=True)

def paired(out):
    assert (out/'.running').read_text()==V
    pub=out/'public/04_ROBUSTNESS';pub.mkdir(exist_ok=False)
    spec=json.loads((out/'public/01_CAMP/analysis_spec.json').read_text())
    js(pub/'analysis_spec.json',spec)
    pp=out/'private/metabolite_pairs_private.tsv';pairs=pd.read_csv(pp,sep='\t',dtype=str)
    ip=out/'public/01_CAMP/feature_inventory.tsv';inv=pd.read_csv(ip,sep='\t');inv=inv[inv.retained_both_author_tissue_filters]
    metp=SRC/'processed_metabolomics/PreprocessedData_PRAD.xlsx';matrix=pd.read_excel(metp,sheet_name='data_imputed',index_col=0);raw=pd.read_excel(metp,sheet_name='data',index_col=0)
    matrix.index=matrix.index.astype(str).str.strip();raw.index=raw.index.astype(str).str.strip()
    assert matrix.index.is_unique and raw.index.is_unique
    rows=[]
    for i,r in inv.iterrows():
        name=r.feature_name; x=matrix.loc[name,pairs.MetabID_tumor].to_numpy(float); y=matrix.loc[name,pairs.MetabID_normal].to_numpy(float)
        finite=np.isfinite(x)&np.isfinite(y)
        av=finite & np.isfinite(raw.loc[name,pairs.MetabID_tumor].to_numpy(float)) & np.isfinite(raw.loc[name,pairs.MetabID_normal].to_numpy(float))
        capt=finite & pairs.CAPT_concordant.eq('True').to_numpy()
        for fam,mask in [('METAB_PAIRED_PRIMARY',finite),('METAB_PAIRED_AVAILABLE',av),('METAB_PAIRED_CAPT_CONCORDANT',capt)]:
            d=base(out,'04_ROBUSTNESS',fam,r.metabolite_key,name)
            d.update(feature_id=r.feature_id,n_pairs_total=len(pairs),chemical_key_collision=r.chemical_key_collision,n_both_available=int(av.sum()),mask_semantics=spec['observed_mask'])
            d.update(paired_stats(x[mask]-y[mask],fam,r.feature_id));rows.append(d)
        if len(rows)%90==0:print('processed_features',len(rows)//3,flush=True)
    d=adjust(pd.DataFrame(rows));assert len(d)==len(inv)*3
    assert (d.n_up+d.n_down+d.n_equal).eq(d.n_pairs_used).all()
    assert not d.duplicated(['feature_id','test_family']).any()
    for fam,fn in [('METAB_PAIRED_PRIMARY','metabolite_paired.tsv'),('METAB_PAIRED_AVAILABLE','metabolite_available_sensitivity.tsv'),('METAB_PAIRED_CAPT_CONCORDANT','metabolite_CAPT_concordant_sensitivity.tsv')]:
        save(d[d.test_family.eq(fam)],pub/fn)
    p=d[d.test_family.eq('METAB_PAIRED_PRIMARY')].copy()
    p['majority_direction_fraction']=p[['up_fraction','down_fraction']].max(axis=1)
    p['direction_discordant']=np.sign(p.mean_delta).ne(np.sign(p.median_delta)) | np.sign(p.mean_delta).ne(np.sign(p.n_up-p.n_down))
    ranked=p.sort_values(['majority_direction_fraction','n','p_value'],ascending=[False,False,True])
    save(ranked,pub/'metabolite_direction_ranking.tsv');save(ranked[ranked.p_value.lt(.05)],pub/'metabolite_workpool.tsv')
    comp=p[['feature_id','metabolite_name','effect','p_value','q_value']].copy()
    for fam,label in [('METAB_PAIRED_AVAILABLE','available'),('METAB_PAIRED_CAPT_CONCORDANT','CAPT')]:
        sub=d[d.test_family.eq(fam)][['feature_id','n','effect','p_value','q_value','status']].rename(columns={c:label+'_'+c for c in ['n','effect','p_value','q_value','status']})
        comp=comp.merge(sub,on='feature_id',validate='one_to_one')
        comp[label+'_same_mean_direction']=np.where(comp[label+'_effect'].notna(),np.sign(comp.effect)==np.sign(comp[label+'_effect']),np.nan)
    save(comp,pub/'sensitivity_comparison.tsv')
    summary={f:dict(planned=len(z),evaluable=int(z.p_value.notna().sum()),p_lt005=int(z.p_value.lt(.05).sum()),q_lt005=int(z.q_value.lt(.05).sum())) for f,z in d.groupby('test_family')}
    js(pub/'validation.json',dict(status='DONE',families=summary,counts_checked=True,feature_identity_unique=True,BH_independent_formula_checked=True,normal_P_independent_formula_checked=True,original_statistics_modified=False))
    manifest([pp,ip,metp,out/'public/01_CAMP/analysis_spec.json',Path(__file__)],pub)
    report(out)

def report(out):
    pub=out/'public/04_ROBUSTNESS'
    summary=json.loads((pub/'validation.json').read_text())['families']
    (pub/'README_CN.md').write_text('# PRAD 全量配对代谢物发现\n\n问题：作者case定义的肿瘤与正常配对有哪些代谢差异？\n\n输入：审计确认的43对，作者两组织过滤交集361特征；CAPT一致36对另列敏感性。使用作者尺度，无新增填补或变换。\n\n实际结果：\n\n'+'\n'.join(f'- {k}: {v}' for k,v in summary.items())+'\n\n新手解释：P<0.05是探索工作池，q<0.05另标FDR支持；升降比例是配对内方向人数比例，不是显著患者比例。效应是作者log2尺度均值差，不称浓度倍数。\n\n限制：作者case身份；CAPT不一致7对保留敏感性；可用值掩码不是已验证原始检出；同队列敏感性不是独立验证。\n\n当前决定：用主分析P<0.05全池进入直接生化映射，不按敏感性显著性删候选。\n\n下一步：全部工作特征映射、全关系肿瘤相关、全基因RNA与单细胞来源。\n\n复现：python prad_discovery_v1.py --mode paired --out AUDITED_RUN。\n',encoding='utf-8')
    checksums(pub);print(json.dumps(summary),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['audit','paired','report'],required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--code-commit',default='UNRECORDED');a=ap.parse_args()
    if a.mode=='audit':audit(a.out,a.code_commit)
    elif a.mode=='paired':paired(a.out)
    else:report(a.out)
