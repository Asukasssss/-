"""Server-only donor-paired expression contrasts; audit before inferential run."""
import argparse, hashlib, json, platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import binomtest, binom
from cell_origin_v1 import ROOT, PREFIX, write_table, sha

KEY=['study','enrichment','technology','annotation_level','cell_type']
JOIN=['study','tissue','enrichment','technology','annotation_level','cell_type','patient']
OLD='20260920T072000Z_cell_expression_v1'
NEW='20260920T081200Z_cell_expression35_v2'

def paired_stats(t,n):
    d=np.asarray(t)-np.asarray(n);size=len(d)
    positive=int((d>0).sum());negative=int((d<0).sum());nonzero=positive+negative
    p=float(binomtest(positive,nonzero,.5).pvalue) if nonzero else 1.
    # Distribution-free median CI, conservative with atoms/ties. Requires >=6 pairs.
    valid=[k for k in range(1,size//2+1) if 1-2*binom.cdf(k-1,size,.5)>=.95]
    k=max(valid) if valid else None;s=np.sort(d)
    return dict(effect=float(np.median(d)),p_value=p,ci_lower=float(s[k-1]) if k else 'NA',
        ci_upper=float(s[size-k]) if k else 'NA',ci_coverage=float(1-2*binom.cdf(k-1,size,.5)) if k else 'NA',
        n_positive=positive,n_negative=negative,n_zero=size-nonzero,n_informative=nonzero)

def bh(p):
    p=np.asarray(p);order=np.argsort(p);v=p[order]*len(p)/np.arange(1,len(p)+1)
    q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return q

def load(base,study):
    paths=[base/run/('private_'+study+'_patient_expression.tsv') for run in [NEW,OLD]]
    a,b=[pd.read_csv(p,sep='\t',keep_default_na=False) for p in paths]
    assert not a.duplicated(JOIN).any() and not b.duplicated(JOIN).any()
    match=a[JOIN+['cells','library_UMI']].merge(b[JOIN+['cells','library_UMI']],on=JOIN,validate='one_to_one',suffixes=('_a','_b'))
    assert len(match)==len(a)==len(b) and (match.cells_a==match.cells_b).all() and (match.library_UMI_a==match.library_UMI_b).all()
    a=a.merge(b[JOIN+['HDC_UMI','HDC_detected','GSTA4_UMI','GSTA4_detected']],on=JOIN,validate='one_to_one')
    assert (a.library_UMI>0).all()
    a['condition']=a.tissue.map({'Tumor':'T','Normal':'N','T':'T','N':'N'});assert a.condition.notna().all()
    return a,[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in paths]

def strata(a):
    for key,g in a.groupby(KEY,sort=True,dropna=False):
        t=g[g.condition=='T'].set_index('patient');n=g[g.condition=='N'].set_index('patient')
        common=t.index.intersection(n.index)
        good=common[(t.loc[common,'cells'].to_numpy()>=20)&(n.loc[common,'cells'].to_numpy()>=20)]
        yield key,t,n,common,good

def coverage(key,t,n,common,good):
    return dict(zip(KEY,key),n_tumor=len(t),n_normal=len(n),n_paired_covered=len(common),n_paired_ge20=len(good),
        n_tumor_only=len(t.index.difference(n.index)),n_normal_only=len(n.index.difference(t.index)),
        n_pairs_excluded_cell_coverage=len(common)-len(good))

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--spec',type=Path,required=True);p.add_argument('--audit-only',action='store_true');p.add_argument('--lock-commit');args=p.parse_args()
    assert args.out.parent==ROOT/'results/collaborative/COAD/B' and (args.out/'.running').exists()
    spec=json.loads(args.spec.read_text());genes=spec['targets'];assert len(genes)==35
    data={};hashes=[];audit=[]
    for study in ['Lee2020_GSE132465','Pelka2021_GSE178341']:
        a,h=load(args.out.parent,study);data[study]=a;hashes.extend(h)
        for values in strata(a):audit.append(coverage(*values))
    if args.audit_only:
        assert not (args.out/'coverage.tsv').exists()
        pd.DataFrame(audit).to_csv(args.out/'coverage.tsv',sep='\t',index=False)
        (args.out/'input_hashes.json').write_text(json.dumps(hashes,indent=2))
        print(json.dumps(dict(audit_strata=len(audit),primary=[r for r in audit if r['annotation_level']=='lineage' and r['enrichment'] in spec['primary_enrichments']])));return
    assert hashes==json.loads((args.out/'input_hashes.json').read_text()) and args.lock_commit
    assert spec['min_cells_each_side']==20 and spec['min_pairs_test']==6
    out=args.out/'public';out.mkdir();(out/'by_gene').mkdir()
    rows=[]
    for study,a in data.items():
        for key,t,n,common,good in strata(a):
            cov=coverage(key,t,n,common,good)
            primary=key[3]=='lineage' and key[1] in spec['primary_enrichments']
            family=study+('_primary_lineage' if primary else '_secondary_other_strata')
            for gene in genes:
                row=dict.fromkeys(PREFIX,'NA');size=len(good);available=gene+'_UMI' in a
                reason='PAIRED_DIRECTION_TEST; not metabolite validation' if size>=6 and available else ('SOURCE_GENE_NOT_FOUND' if not available else 'FEWER_THAN_6_PAIRED_PATIENTS_GE20; no test')
                row.update(cancer='COAD',cohort=study,stage_id='06_EXTERNAL',run_id=args.out.name,analysis_version=spec['analysis_version'],
                    analysis_type='paired_T_minus_N_within_cell_type',gene=gene,unit='paired_patient',n=size,n_reference=len(common),
                    effect_type='median_paired_CPM_difference_T_minus_N',test_family=family,status='DONE' if size>=6 and available else 'NOT_EVALUABLE',reason=reason,source_id=study)
                row.update(cov,tier='primary' if primary else 'secondary',gene_available=available,ci_method='order_statistic_median_95pct_or_higher',
                    ci_coverage='NA',n_positive='NA',n_negative='NA',n_zero='NA',n_informative='NA',
                    paired_tumor_CPM_median='NA',paired_normal_CPM_median='NA',median_detection_difference='NA')
                if available and size>=3:
                    tc=1e6*t.loc[good,gene+'_UMI'].to_numpy()/t.loc[good,'library_UMI'].to_numpy()
                    nc=1e6*n.loc[good,gene+'_UMI'].to_numpy()/n.loc[good,'library_UMI'].to_numpy()
                    row.update(paired_tumor_CPM_median=float(np.median(tc)),paired_normal_CPM_median=float(np.median(nc)),effect=float(np.median(tc-nc)),
                        median_detection_difference=float(np.median(t.loc[good,gene+'_detected'].to_numpy()/t.loc[good,'cells'].to_numpy()-n.loc[good,gene+'_detected'].to_numpy()/n.loc[good,'cells'].to_numpy())))
                    if size>=6:row.update(paired_stats(tc,nc))
                rows.append(row)
    frame=pd.DataFrame(rows)
    for family,g in frame.groupby('test_family'):
        valid=g.index[g.status=='DONE'];frame.loc[g.index,'family_n_evaluable']=len(valid)
        frame.loc[valid,'q_value']=bh(frame.loc[valid,'p_value'].to_numpy(dtype=float))
    frame=frame.sort_values(['gene']+KEY).reset_index(drop=True)
    for gene,g in frame.groupby('gene',sort=True):write_table(g.to_dict('records'),out/'by_gene'/(gene+'.tsv'))
    write_table(frame[frame.tier=='primary'].to_dict('records'),out/'primary_results.tsv')
    pd.DataFrame(audit).to_csv(out/'coverage.tsv',sep='\t',index=False)
    summaries=[]
    for family,g in frame.groupby('test_family'):
        test=g[g.status=='DONE'];summaries.append(dict(family=family,planned=len(g),evaluable=len(test),nominal_p_lt05=int((pd.to_numeric(test.p_value)<.05).sum()),q_lt05=int((pd.to_numeric(test.q_value)<.05).sum()),n_min=int(test.n.min()) if len(test) else None,n_max=int(test.n.max()) if len(test) else None))
    assert not frame.duplicated(['gene']+KEY).any() and set(frame.gene)==set(genes)
    assert frame.loc[frame.gene=='GSTT2','p_value'].eq('NA').all()
    assert all((float(r['ci_lower'])<=float(r['effect'])<=float(r['ci_upper'])) for r in rows if r['ci_lower']!='NA')
    for name,obj in [('summary.json',summaries),('input_hashes.json',hashes),('analysis_spec.json',spec),('validation.json',dict(status='PASS',rows=len(frame),genes=35,primary_rows=int((frame.tier=='primary').sum()),lock_commit=args.lock_commit,
        code_sha256=sha(__file__),spec_sha256=sha(args.spec),python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,input_hashes_unchanged=True,
        old_new_aggregate_keys_counts_identical=True,individual_rows_exported=False,tests='synthetic paired kernels; unique strata;missing gene;interval order',scope='paired CPM direction consistency; not count GLM or absolute expression'))]:
        (out/name).write_text(json.dumps(obj,indent=2)+'\n')
    print(json.dumps(summaries))

if __name__=='__main__':main()
