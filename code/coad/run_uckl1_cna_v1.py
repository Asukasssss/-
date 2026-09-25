"""Server-only UCKL1 comparison using published Uhlitz CNA calls.

No de novo malignancy inference. Individual cells and donor values stay private.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binom, binomtest
import scipy

ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
OLD = ROOT / 'results/collaborative/COAD/B/20260922T141755Z_source_contract_sc_v2'
PREFIX = 'cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
GROUPS = ['tumor_CNA', 'tumor_CNN', 'normal_reference']

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def js(p, x):
    p.write_text(json.dumps(x, ensure_ascii=False, indent=2)+'\n')

def write(p, rows):
    pd.DataFrame(rows).to_csv(p, sep='\t', index=False, na_rep='NA')

def build_labels():
    paths = dict(metadata=ROOT/'data/candidates/coad_uhlitz_20260921/metadata.tsv.gz',
        calls=ROOT/'data/candidates/coad_uckl1_malignancy_20260925/infercnv_clone_scores.tsv',
        cell_metadata=OLD/'private_cell_metadata.tsv', counts=OLD/'private_cell_counts.npz')
    cols=['cell','cell_id','source_id','case_id','sample_origin','main_cell_type','cell_type_epi_custom']
    m=pd.read_csv(paths['metadata'],sep='\t',usecols=cols,dtype=str,keep_default_na=False)
    a=pd.read_csv(paths['calls'],sep='\t',dtype=str)
    # Original infercnv.R replaces ':' with '_' before assigning cell identifiers.
    m['join_id']=m.cell_id.str.replace(':','_',regex=False)
    assert m.cell_id.is_unique and m.join_id.is_unique and (m.cell==m.cell_id).all()
    valid=a[a.cell_id.notna()].copy()
    assert valid.cell_id.is_unique and set(valid.cna_clone)<= {'CNA','CNN'}
    j=m.merge(valid,left_on='join_id',right_on='cell_id',how='left',suffixes=('','_call'),validate='one_to_one')
    assert j.cell_id_call.notna().sum()==len(valid)
    labelled=j.cell_id_call.notna()
    assert (j.loc[labelled,'main_cell_type']=='Epithelial').all()
    assert (j.loc[labelled,'sample_origin']=='Tumor').all()
    conflict=labelled & (j.source_id!=j.source_id_call)
    excluded=set(j.loc[conflict,'case_id'])
    j['group']='not_selected'
    j.loc[labelled & (j.cna_clone=='CNA'),'group']='tumor_CNA'
    j.loc[labelled & (j.cna_clone=='CNN'),'group']='tumor_CNN'
    normal=(j.main_cell_type=='Epithelial') & (j.sample_origin=='Normal') & ~j.cell_type_epi_custom.str.startswith('TC')
    j.loc[normal,'group']='normal_reference'
    j.loc[j.case_id.isin(excluded),'group']='excluded_donor_identity_conflict'
    coverage=[]
    for group in GROUPS:
        z=j[j.group==group].groupby('case_id').size()
        coverage.append(dict(group=group,cells=int(z.sum()),donors_any=int(len(z)),donors_ge20=int((z>=20).sum()),
            cells_in_qualified_donors=int(z[z>=20].sum()),minimum_cells=int(z.min()) if len(z) else 'NA'))
    audit=dict(status='PASS',metadata_cells=len(m),author_call_rows=len(a),call_rows_with_cell_id=len(valid),
        call_rows_without_cell_id=int(a.cell_id.isna().sum()),matched_call_cells=int(labelled.sum()),
        donor_id_conflict_cells=int(conflict.sum()),excluded_donors=len(excluded),
        join_rule="Exact author documented colon-to-underscore cell ID normalization; exclude entire donor on source_id mismatch",
        all_labelled_cells_tumor_epithelial=True,normal_rule="Normal tissue epithelial; exclude author cell_type_epi_custom starting TC, following infercnv.R reference rule",
        CNN_is_not_proven_nonmalignant=True,no_expression_used_in_admission=True)
    return j,paths,coverage,audit

def median_ci(x):
    x=np.sort(x);n=len(x)
    ks=[k for k in range(1,n//2+1) if 2*binom.cdf(k-1,n,.5)<=.05]
    if not ks:return float('nan'),float('nan'),float('nan')
    k=max(ks)
    return float(x[k-1]),float(x[n-k]),float(1-2*binom.cdf(k-1,n,.5))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--spec',type=Path,required=True);ap.add_argument('--mode',choices=['audit','analyze'],required=True)
    ap.add_argument('--lock-commit',default='NOT_RUN');a=ap.parse_args()
    assert a.out.parent.resolve()==ROOT/'results/collaborative/COAD/B'
    spec=json.loads(a.spec.read_text());pub=a.out/'public';pub.mkdir(exist_ok=True)
    assert (a.out/'.running').is_dir()
    j,paths,coverage,audit=build_labels()
    hashes={k:sha(p) for k,p in paths.items()}
    if a.mode=='audit':
        js(pub/'label_audit.json',audit);write(pub/'coverage.tsv',coverage)
        js(a.out/'admission_hashes.json',hashes)
        print(json.dumps(dict(audit=audit,coverage=coverage)));return
    assert json.loads((a.out/'admission_hashes.json').read_text())==hashes
    assert len(a.lock_commit)==40
    meta=pd.read_csv(paths['cell_metadata'],sep='\t',dtype=str,keep_default_na=False)
    assert meta.cell_id.is_unique
    # Align by explicit IDs rather than row order.
    ix=pd.Index(meta.cell_id).get_indexer(j.cell_id)
    assert (ix>=0).all() and len(set(ix))==len(ix)==len(meta)
    assert (meta.iloc[ix].patient.to_numpy()==j.case_id.to_numpy()).all()
    assert (meta.iloc[ix].tissue.to_numpy()==j.sample_origin.to_numpy()).all()
    with np.load(paths['counts']) as z:
        count=z['UCKL1'][ix].astype(float);total=z['total'][ix].astype(float)
    assert np.isfinite(count).all() and (count>=0).all() and (total>0).all() and (count<=total).all()
    j['value']=np.log1p(count/total*10000);j['detected']=count>0;j['counts']=count;j['total']=total
    donor=[]
    for (patient,group),f in j[j.group.isin(GROUPS)].groupby(['case_id','group']):
        donor.append(dict(patient=patient,group=group,n_cells=len(f),mean_log1p_CP10K=f.value.mean(),
            detection_fraction=f.detected.mean(),pseudobulk_CPM=f.counts.sum()/f.total.sum()*1e6))
    d=pd.DataFrame(donor);write(a.out/'private_donor_values.tsv',donor)
    eligible=d[d.n_cells>=spec['min_cells']]
    summaries=[]
    for group,f in eligible.groupby('group'):
        summaries.append(dict(group=group,n_donors=len(f),n_cells=int(f.n_cells.sum()),
            mean_expression=float(f.mean_log1p_CP10K.mean()),median_expression=float(f.mean_log1p_CP10K.median()),
            expression_q25=float(f.mean_log1p_CP10K.quantile(.25)),expression_q75=float(f.mean_log1p_CP10K.quantile(.75)),
            mean_detection_fraction=float(f.detection_fraction.mean()),median_pseudobulk_CPM=float(f.pseudobulk_CPM.median())))
    write(pub/'group_summary.tsv',summaries)
    rows=[]
    for c in spec['contrasts']:
        x=eligible[eligible.group==c['a']].set_index('patient');y=eligible[eligible.group==c['b']].set_index('patient')
        ids=x.index.intersection(y.index);n=len(ids)
        dx=(x.loc[ids,'mean_log1p_CP10K']-y.loc[ids,'mean_log1p_CP10K']).to_numpy()
        up=int((dx>0).sum());down=int((dx<0).sum());zero=int((dx==0).sum())
        row=dict.fromkeys(PREFIX,'NA')
        row.update(cancer='COAD',cohort='Uhlitz_GSE166555',stage_id='06_EXTERNAL',run_id=a.out.name,
            analysis_version=spec['analysis_version'],analysis_type='UCKL1_author_CNA_patient_paired',gene='UCKL1',unit='patient',
            n=n,n_reference=n,effect_type='median_patient_difference_mean_log1p_CP10K',test_family='UCKL1_Uhlitz_three_prespecified_contrasts',
            status='DONE' if n>=spec['min_pairs'] else 'NOT_EVALUABLE',reason='AUTHOR_CNA_VS_CNN_NOT_DEFINITIVE_MALIGNANCY',
            source_id='GSE166555;author infercnv_clone_scores dbac4154',contrast=c['id'],group_a=c['a'],group_b=c['b'],
            n_up=up,n_down=down,n_equal=zero,n_nonzero=up+down,group_a_unpaired_donors=len(x)-n,group_b_unpaired_donors=len(y)-n)
        if n>=3:
            row.update(effect=float(np.median(dx)),mean_difference=float(np.mean(dx)),
                paired_group_a_mean=float(x.loc[ids,'mean_log1p_CP10K'].mean()),paired_group_b_mean=float(y.loc[ids,'mean_log1p_CP10K'].mean()),
                median_detection_difference=float(np.median(x.loc[ids,'detection_fraction']-y.loc[ids,'detection_fraction'])),
                median_pseudobulk_CPM_difference=float(np.median(x.loc[ids,'pseudobulk_CPM']-y.loc[ids,'pseudobulk_CPM'])))
        if n>=spec['min_pairs']:
            row['p_value']=float(binomtest(up,up+down,p=.5).pvalue) if up+down else 1.
            row['ci_lower'],row['ci_upper'],row['ci_coverage']=median_ci(dx)
            row['leave_one_patient_out_median_min']=float(min(np.median(np.delete(dx,k)) for k in range(n)))
            row['leave_one_patient_out_median_max']=float(max(np.median(np.delete(dx,k)) for k in range(n)))
        else:row['reason']='FEWER_THAN_6_QUALIFIED_PAIRS;CNN_NOT_PROVEN_NONMALIGNANT'
        rows.append(row)
    p=np.array([float(r['p_value']) if r['status']=='DONE' else 1. for r in rows]);order=np.argsort(p)
    q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
    for k,r in enumerate(rows):
        r['family_n_evaluable']=sum(x['status']=='DONE' for x in rows)
        if r['status']=='DONE':r['q_value']=float(q[k])
    write(pub/'results.tsv',rows);write(pub/'coverage.tsv',coverage);js(pub/'label_audit.json',audit)
    js(pub/'analysis_spec.json',dict(spec,code_lock_commit=a.lock_commit))
    write(pub/'source_manifest.tsv',[dict(source_id=k,source_path=str(p),sha256=hashes[k]) for k,p in paths.items()])
    js(pub/'validation.json',dict(status='PASS',input_hashes_match_admission=True,explicit_cell_and_patient_joins=True,
        public_tables_contain_no_patient_or_cell_rows=True,planned_tests=3,evaluable_tests=sum(r['status']=='DONE' for r in rows),
        donor_identity_conflicts_excluded=audit['excluded_donors'],normal_label_missing_ID_not_inferred=True,
        software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
        caveats=['CNA evidence from RNA not per-cell DNA','CNN not proven benign','No independent cohort replication']))
    print(json.dumps(rows,ensure_ascii=False));(a.out/'DONE').write_text('DONE\n');(a.out/'.running').rmdir()

if __name__=='__main__':main()
