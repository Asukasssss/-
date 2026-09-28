"""Read-only consistency checks against private unit outputs; export checks only."""
from pathlib import Path
import itertools,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent
idx=pd.read_csv(R/'private/units_index.tsv',sep='\t')
variants=['mixed','T_only','B_only'];checked=0;cells=0;myeloid=0;paired={}
for _,row in idx.iterrows():
    p=Path(row.path)
    if not (p/'INFERENCE_DONE').exists():continue
    meta=pd.read_csv(p/'metadata.tsv',sep='\t').set_index('cell');con=pd.read_csv(p/'consensus_private.tsv',sep='\t').set_index('cell')
    assert meta.index.is_unique and set(meta.index)==set(con.index)
    assert all(meta.slc6a6_log1p10k>=0)
    implied=np.expm1(meta.slc6a6_log1p10k)*meta.total_umi/10000
    np.testing.assert_allclose(implied,np.rint(implied),atol=1e-7,rtol=1e-10)
    assert ((implied>0)==meta.slc6a6_detected.astype(bool)).all()
    variant_states=[];variant_signs=[]
    for v in variants:
        s=pd.read_csv(p/f'scores_{v}.tsv',sep='\t').set_index('cell')
        ch=pd.read_csv(p/f'chromosome_means_{v}.tsv',sep='\t').set_index('cell')
        th=pd.read_csv(p/f'chr_thresholds_{v}.tsv',sep='\t').set_index('chr').threshold99
        assert s.index.is_unique and set(s.index)==set(meta.index)==set(ch.index)
        assert 'chr3' not in ch and 'SLC6A6' not in (p/'genes.txt').read_text().splitlines()
        n=(ch.abs()>th.clip(lower=1e-8)).sum(axis=1)
        np.testing.assert_array_equal(n.reindex(s.index),s.n_abnormal_chr)
        state=np.where((s.score>s.threshold99)&(n.reindex(s.index)>=2),'CNV_ABNORMAL_CANDIDATE',np.where((s.score<=s.threshold95)&(n.reindex(s.index)==0),'REFERENCE_LIKE','UNCERTAIN'))
        assert (state==s.state).all()
        variant_states.append(s.state.reindex(meta.index))
        ordered=ch.reindex(index=meta.index,columns=sorted(ch.columns))
        variant_signs.append(np.where(ordered.abs()>th.reindex(ordered.columns).clip(lower=1e-8),np.sign(ordered),0))
    st=pd.concat(variant_states,axis=1);a=np.stack(variant_signs)
    shared=((a[0]!=0)&(a[0]==a[1])&(a[1]==a[2])).sum(axis=1)
    expected=np.where(st.eq('CNV_ABNORMAL_CANDIDATE').all(axis=1)&(shared>=2),'CNV_ABNORMAL_CANDIDATE',np.where(st.eq('REFERENCE_LIKE').all(axis=1),'REFERENCE_LIKE','UNCERTAIN'))
    np.testing.assert_array_equal(expected,con.reindex(meta.index).consensus_state)
    m=con[con.cnv_label=='Myeloid'];myeloid+=len(m);cells+=len(meta);checked+=1
    summary=m.groupby('consensus_state').agg(n=('slc6a6_log1p10k','size'),mean=('slc6a6_log1p10k','mean'))
    keys=['CNV_ABNORMAL_CANDIDATE','REFERENCE_LIKE']
    if all(k in summary.index and summary.loc[k,'n']>=20 for k in keys):paired.setdefault(row.cohort,[]).append(summary.loc[keys[0],'mean']-summary.loc[keys[1],'mean'])
pub=pd.read_csv(R/'public/group_summary.tsv',sep='\t');assert int(pub.n_cells.sum())==myeloid
tests=pd.read_csv(R/'public/SLC6A6_comparison.tsv',sep='\t')
for _,t in tests.iterrows():
    d=np.array(paired.get(t.cohort,[]),float);assert len(d)==t.n_paired_author_units
    if len(d)>=5:
        signs=np.array(list(itertools.product([-1,1],repeat=len(d))))
        p=np.mean(np.abs(signs.dot(d)/len(d))>=abs(d.mean())-1e-12)
        np.testing.assert_allclose([p,d.mean()],[t.p_value,t.mean_delta_log1p10k],rtol=1e-9,atol=1e-12)
    else:assert t.status=='NOT_EVALUABLE' and pd.isna(t.p_value) and pd.isna(t.q_value)
e=tests[tests.p_value.notna()].sort_values('p_value');vals=e.p_value.to_numpy();m=len(vals)
q=[min(1,min(m*vals[j]/(j+1) for j in range(i,m))) for i in range(m)]
np.testing.assert_allclose(e.q_value,q,rtol=1e-9,atol=1e-12)
v={'status':'PASS','units_checked':checked,'cells_checked':cells,'myeloid_cells_checked':myeloid,'per_unit_reference_runs':3,'SLC6A6_counts_integrality':'PASS','chr3_excluded':'PASS','cell_keys_and_per_variant_labels':'PASS','public_cell_count_reconciliation':'PASS','paired_donor_statistics_and_BH':'PASS','not_validated':['true malignancy','DNA CNV','doublets','clinical donor independence','myeloid subtype confounding']}
(R/'public/independent_validation.json').write_text(json.dumps(v,indent=2)+'\n')
print(json.dumps(v),flush=True)
