"""Server-only paired change/change analysis; public aggregate outputs only.

Patient measurements and the audited pairing never leave server165. Figures
contain unlabeled rasterized points, not an export of sample identifiers/values.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '2'
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
BASE = ROOT / 'results/collaborative/BRCA/A'
VERSION = 'lypla1_paired_delta_v1'
FEATURES = ['1-palmitoyl-GPC (16:0)', 'glycerophosphorylcholine (GPC)']
PREFIX = 'cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(frame, p):
    frame.to_csv(p, sep='\t', index=False, na_rep='NA', lineterminator='\n')


def rank_rho(x, y):
    a = stats.rankdata(x) - (len(x) + 1) / 2
    b = stats.rankdata(y) - (len(y) + 1) / 2
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def calculate(x, y, seed, nperm, nboot):
    rho = rank_rho(x, y)
    assert abs(rho - stats.spearmanr(x, y).statistic) < 1e-12
    rng = np.random.default_rng(seed)
    a = stats.rankdata(x) - (len(x) + 1) / 2
    b = stats.rankdata(y) - (len(y) + 1) / 2
    den = np.sqrt(np.dot(a, a) * np.dot(b, b))
    extreme = 0
    for start in range(0, nperm, 1000):
        perm = np.stack([rng.permutation(b) for _ in range(min(1000, nperm - start))])
        extreme += int((np.abs(perm @ a / den) >= abs(rho) - 1e-12).sum())
    p = (extreme + 1) / (nperm + 1)
    # Resample complete paired-case rows, keeping both delta values together.
    rng = np.random.default_rng(seed + 1)
    ix = rng.integers(len(x), size=(nboot, len(x)))
    bx, by = stats.rankdata(x[ix], axis=1), stats.rankdata(y[ix], axis=1)
    bx -= bx.mean(1, keepdims=True)
    by -= by.mean(1, keepdims=True)
    den_b = np.sqrt((bx * bx).sum(1) * (by * by).sum(1))
    valid = den_b > 0
    bs = (bx[valid] * by[valid]).sum(1) / den_b[valid]
    lo, hi = np.quantile(bs, [.025, .975])
    return dict(effect=rho, p_value=p, ci_lower=lo, ci_upper=hi,
                permutation_extreme=extreme, permutation_n=nperm,
                permutation_mc_se=np.sqrt(p * (1-p) / (nperm+1)),
                bootstrap_valid=int(valid.sum()), bootstrap_n=nboot, seed=seed)


def plots(result, data, counts, out):
    plt.rcParams.update({'pdf.fonttype':42, 'ps.fonttype':42, 'font.size':10,
                         'axes.spines.top':False, 'axes.spines.right':False,
                         'savefig.dpi':180, 'font.family':'DejaVu Sans'})
    figures = []
    fig, axs = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    for j, name in enumerate(FEATURES):
        d = data[name]
        for i, mode in enumerate(['processed', 'author_available']):
            ax = axs[i,j]
            mask = d['finite'] if mode == 'processed' else d['available']
            row = result[(result.metabolite_name == name) & (result.analysis_type == mode)].iloc[0]
            ax.axhline(0, color='#999999', lw=.8)
            ax.axvline(0, color='#999999', lw=.8)
            if mode == 'processed':
                missing = mask & ~d['available']
                ax.scatter(d['dx'][missing], d['dy'][missing], s=32, color='#999999',
                           alpha=.85, edgecolors='white', linewidths=.3, rasterized=True,
                           label='At least one author value unavailable')
            observed = mask & d['available']
            ax.scatter(d['dx'][observed], d['dy'][observed], s=32, color='#0072B2',
                       alpha=.85, edgecolors='white', linewidths=.3, rasterized=True,
                       label='Both author values available')
            for axis, values in [('x',d['dx'][d['finite']]), ('y',d['dy'][d['finite']])]:
                low, high = min(0, values.min()), max(0, values.max())
                pad = max((high-low)*.08, .05)
                getattr(ax, 'set_'+axis+'lim')(low-pad,high+pad)
            title = 'LPC 16:0' if j == 0 else 'Free GPC'
            subset = 'Processed main' if i == 0 else 'Author-available subset'
            ax.set_title(f'{title} | {subset} | n={int(row.n)}', weight='bold')
            ax.text(.03,.97, f"rho = {row.effect:+.3f}\n95% CI [{row.ci_lower:+.3f}, {row.ci_upper:+.3f}]\nP = {row.p_value:.4g}; BH4 q = {row.q_value:.4g}",
                    transform=ax.transAxes, va='top', fontsize=9,
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85))
            ax.set_xlabel('LYPLA1 RNA change (tumor - own normal; author scale)')
            ax.set_ylabel('Metabolite change (tumor - own normal; author scale)')
            if i == 0 and j == 0: ax.legend(loc='lower right',fontsize=7,frameon=False)
    fig.suptitle('LYPLA1: paired RNA change versus metabolite change\nSame cohort exploration; differences are not enzyme activity or metabolic flux',fontsize=15,weight='bold')
    fig.savefig(out/'paired_delta_scatter.png')
    figures.append(fig)

    fig, axs = plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for j,name in enumerate(FEATURES):
        for i,mode in enumerate(['processed','author_available']):
            ax = axs[i,j]
            sub=counts[(counts.metabolite_name==name)&(counts.analysis_type==mode)]
            z=sub.pivot(index='metabolite_direction',columns='RNA_direction',values='count').reindex(index=['UP','EQUAL','DOWN'],columns=['DOWN','EQUAL','UP']).to_numpy()
            ax.imshow(z,cmap='Blues',vmin=0,vmax=45,aspect='auto')
            total=int(z.sum())
            for y in range(3):
                for x in range(3):
                    value=int(z[y,x]);ax.text(x,y,f'{value}/{total}\n{value/total:.1%}',ha='center',va='center',color='white' if value>25 else '#172B4D',fontsize=11)
            ax.set_xticks(range(3),['Lower','Equal','Higher']);ax.set_yticks(range(3),['Higher','Equal','Lower'])
            ax.set_xlabel('LYPLA1 RNA in tumor versus own normal')
            ax.set_ylabel('Metabolite in tumor versus own normal')
            ax.set_title(('LPC 16:0' if j==0 else 'Free GPC')+' | '+('Processed main' if i==0 else 'Author-available'),weight='bold')
    fig.suptitle('Joint directions in the same paired cases\nCounts are descriptive: joint increase is not a correlation test or individual significance',fontsize=14,weight='bold')
    fig.savefig(out/'joint_directions.png')
    figures.append(fig)
    with PdfPages(out/'LYPLA1_paired_changes.pdf') as pdf:
        for f in figures: pdf.savefig(f)
    for f in figures: plt.close(f)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--source-commit',required=True)
    a=ap.parse_args(); out=a.out; pub=out/'public'; private=out/'private'
    assert (out/'.running').read_text().strip()==VERSION
    assert not (pub/'results.tsv').exists(), 'No overwrite of existing results'
    spec=dict(analysis_version=VERSION, source_commit=a.source_commit,
              gene='LYPLA1', metabolites=FEATURES,
              question='Across verified paired cases, do RNA tumor-minus-normal differences correlate with metabolite tumor-minus-normal differences?',
              unit='author_case_row with audited tissue-consistent tumor-normal pairing',
              primary='unchanged author processed RNA and metabolite scales; no new imputation or transform',
              sensitivity='same processed differences; restrict metabolite to finite values on both sides in author data sheet; not verified original detection mask',
              planned_tests=4, multiplicity='one BH4 family: 2 metabolites x 2 modes, including identical GPC sensitivity; retain all four; do not overwrite historical q',
              test='two-sided absolute Spearman permutation, 99999 shuffles of metabolite differences across cases, plus-one P',
              interval='10000 whole-case percentile bootstrap; pointwise95%; not multiplicity-adjusted',
              seed='sha256(version|exact feature|selected case IDs); identical data subsets use identical seed',
              direction_counts='all 9 up/equal/down combinations; equality means exact zero in unchanged scale; no direction-count hypothesis tests',
              min_pairs=8, no_new_outlier_exclusion=True,
              postselection_exploratory=True, independent_validation=False,
              patient_tables_server_only=True, public_scatter='unlabeled rasterized points only; no public patient identifiers or numeric arrays',
              python=platform.python_version(), numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__)
    (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
    identity=BASE/'20260921T102429Z_camp_sample_identity_v1/private/audited_mapping.tsv'
    olddir=BASE/'20260921T150000Z_mapping262_patient_v1/public'
    oldmet=BASE/'20260921T112611Z_metabolite318_paired_v1/public/paired_metabolite318.tsv'
    m=pd.read_csv(identity,sep='\t',dtype=str)
    clean=m[m.TN.eq(m.pdf_TN)]
    t=clean[clean.TN.eq('Tumor')]; n=clean[clean.TN.eq('Normal')]
    cols=['case_row','RNAID','MetabID']
    pairs=t[cols].merge(n[cols],on='case_row',suffixes=('_t','_n'),validate='one_to_one').sort_values('case_row').reset_index(drop=True)
    assert len(pairs)==45 and pairs.case_row.is_unique
    assert pairs[['RNAID_t','RNAID_n']].stack().nunique()==90
    assert pairs[['MetabID_t','MetabID_n']].stack().nunique()==90
    src=ROOT/'data/candidates/camp_primary_tissue_multicancer'
    assert m.RNAFile.nunique()==1 and m.MetabFile.nunique()==1
    rf=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
    mf=src/'processed_metabolomics'/m.MetabFile.iloc[0]
    paths=[identity,rf,mf,olddir/'paired_RNA150.tsv',olddir/'CAMP_relations262.tsv',oldmet,Path(__file__)]
    hashes={str(p):sha(p) for p in paths}
    # Check the three actual inputs against the previously delivered manifest.
    previous=pd.read_csv(olddir/'source_manifest.tsv',sep='\t').set_index('path').sha256.to_dict()
    for p in [identity,rf,mf]: assert hashes[str(p)]==previous[str(p)],'Historical input hash changed'
    rna=pd.read_csv(rf,index_col=0)
    mat=pd.read_excel(mf,sheet_name='data_imputed',index_col=0)
    avail=pd.read_excel(mf,sheet_name='data',index_col=0)
    for frame in [rna,mat,avail]:
        frame.columns=frame.columns.astype(str)
        frame.index=frame.index.astype(str).str.strip()
        assert frame.index.is_unique and frame.columns.is_unique
    dx=rna.loc['LYPLA1',pairs.RNAID_t].to_numpy(float)-rna.loc['LYPLA1',pairs.RNAID_n].to_numpy(float)
    oldrna=pd.read_csv(olddir/'paired_RNA150.tsv',sep='\t').set_index('gene').loc['LYPLA1']
    assert np.isfinite(dx).all() and (dx>0).sum()==oldrna.positive_pairs==34
    assert abs(dx.mean()-oldrna.effect)<1e-12
    baseline=pd.read_csv(oldmet,sep='\t')
    oldrel=pd.read_csv(olddir/'CAMP_relations262.tsv',sep='\t')
    oldrel=oldrel[oldrel.gene.eq('LYPLA1') & oldrel.metabolite_name.isin(FEATURES)].copy()
    save(oldrel,pub/'historical_tumor_associations.tsv')
    data={}; rows=[]; counts=[]; private_rows=[]; checks=[]
    for name in FEATURES:
        dy=mat.loc[name,pairs.MetabID_t].to_numpy(float)-mat.loc[name,pairs.MetabID_n].to_numpy(float)
        finite=np.isfinite(dx)&np.isfinite(dy)
        observed=finite&np.isfinite(avail.loc[name,pairs.MetabID_t].to_numpy(float))&np.isfinite(avail.loc[name,pairs.MetabID_n].to_numpy(float))
        assert finite.sum()==45
        data[name]=dict(dx=dx,dy=dy,finite=finite,available=observed)
        for i,pair in pairs.iterrows():
            private_rows.append(dict(**pair.to_dict(),metabolite_name=name,RNA_delta=dx[i],metabolite_delta=dy[i],both_author_available=bool(observed[i])))
        for mode, mask, oldfamily in [('processed',finite,'paired_processed'),('author_available',observed,'paired_author_available')]:
            b=baseline[baseline.metabolite_name.eq(name)&baseline.analysis_type.eq(oldfamily)].iloc[0]
            x,y=dx[mask],dy[mask]; nn=len(x)
            assert nn==b.n and (y>0).sum()==b.pairs_higher and (y<0).sum()==b.pairs_lower
            assert abs(y.mean()-b.effect)<1e-12
            checks.append(dict(metabolite_name=name,analysis_type=mode,n=nn,historical_marginals_match=True))
            row=dict.fromkeys(PREFIX,np.nan)
            row.update(cancer='BRCA',cohort='CAMP_BRCA1_Terunuma',stage_id='03_PATIENT',run_id=out.name,analysis_version=VERSION,analysis_type=mode,
                       metabolite_key=b.metabolite_key,metabolite_name=name,gene='LYPLA1',unit='audited_author_paired_case',n=nn,n_reference=45,
                       effect_type='Spearman_of_paired_RNA_and_metabolite_differences',test_family='LYPLA1_delta_BH4',family_n_evaluable=4,
                       status='NOT_EVALUABLE',reason='less_than8_or_constant',source_id='CAMP_GSE37751_author_case_pairing',
                       rna_higher=int((x>0).sum()),rna_lower=int((x<0).sum()),rna_equal=int((x==0).sum()),
                       metabolite_higher=int((y>0).sum()),metabolite_lower=int((y<0).sum()),metabolite_equal=int((y==0).sum()),
                       historical_metabolite_p=b.p_value,historical_metabolite_q=b.q_value,
                       mean_RNA_delta=float(x.mean()),mean_metabolite_delta=float(y.mean()),
                       historical_RNA_main_p=oldrna.p_value,historical_RNA_main_q=oldrna.q_value)
            for xs,xd in [(1,'UP'),(0,'EQUAL'),(-1,'DOWN')]:
                for ys,yd in [(1,'UP'),(0,'EQUAL'),(-1,'DOWN')]:
                    count=int(((np.sign(x)==xs)&(np.sign(y)==ys)).sum())
                    counts.append(dict(metabolite_name=name,analysis_type=mode,RNA_direction=xd,metabolite_direction=yd,count=count,n=nn,fraction=count/nn))
            if nn>=8 and np.ptp(x)>0 and np.ptp(y)>0:
                sd=int(hashlib.sha256((VERSION+'|'+name+'|'+','.join(pairs.loc[mask,'case_row'])).encode()).hexdigest()[:8],16)
                row.update(calculate(x,y,sd,99999,10000))
                row.update(status='DONE',reason='same_cohort_postselection;not_cell_intrinsic_or_causal')
            rows.append(row)
    result=pd.DataFrame(rows)
    ok=result.p_value.notna();result.loc[ok,'q_value']=multipletests(result.loc[ok,'p_value'],method='fdr_bh')[1]
    result['family_n_evaluable']=int(ok.sum())
    assert ok.all(), 'Stop for explicit review if planned test is not evaluable'
    ordered=np.argsort(result.p_value.to_numpy()); ps=result.p_value.to_numpy()[ordered]
    independent=np.minimum(1,np.minimum.accumulate((ps*4/np.arange(1,5))[::-1])[::-1])
    assert np.allclose(independent,result.q_value.to_numpy()[ordered],atol=1e-12,rtol=0)
    c=pd.DataFrame(counts)
    assert c.groupby(['metabolite_name','analysis_type'])['count'].sum().eq(result.set_index(['metabolite_name','analysis_type']).n).all()
    save(result[PREFIX+[col for col in result if col not in PREFIX]],pub/'results.tsv')
    save(c,pub/'joint_directions.tsv')
    save(pd.DataFrame(private_rows),private/'paired_differences.tsv')
    save(pairs,private/'audited_pairs.tsv')
    plots(result,data,c,pub)
    with pd.ExcelWriter(pub/'LYPLA1_paired_changes.xlsx',engine='openpyxl') as writer:
        result.to_excel(writer,sheet_name='4项变化量关联',index=False)
        c.to_excel(writer,sheet_name='共同变化人数',index=False)
        oldrel.to_excel(writer,sheet_name='既有肿瘤关联_非新检验',index=False)
        pd.DataFrame({'项目':list(spec),'说明':[json.dumps(v,ensure_ascii=False) for v in spec.values()]}).to_excel(writer,sheet_name='方法与边界',index=False)
        for ws in writer.book.worksheets:
            ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
            for col in ws.columns: ws.column_dimensions[col[0].column_letter].width=24
    assert all(sha(p)==h for p,h in hashes.items())
    save(pd.DataFrame([dict(path=p,sha256=h) for p,h in hashes.items()]),pub/'source_manifest.tsv')
    validation=dict(status='DONE',pair_count=45,input_hashes_match_previous=True,input_hashes_unchanged=True,
                    RNA_mean_and_direction_match_history=True,metabolite_checks=checks,rank_rho_independently_checked=True,
                    BH4_independently_checked=True,direction_counts_sum_to_n=True,patient_numeric_exports_server_only=True,
                    historical_statistics_modified=False,scope='Numerical checks on these four tests; no re-audit of full historical project',
                    independent_R_audit='PENDING')
    (pub/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    (out/'NUMERICAL_DONE').write_text('DONE')
    print(result[['metabolite_name','analysis_type','n','effect','ci_lower','ci_upper','p_value','q_value']].to_string(index=False),flush=True)
    print(c[c['count']>0].to_string(index=False),flush=True)


if __name__=='__main__':
    main()
