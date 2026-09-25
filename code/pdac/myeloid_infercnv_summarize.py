"""Aggregate private CNV-inference labels and SLC6A6 measurements by author unit."""
from pathlib import Path
import itertools,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent
VARIANTS=['mixed','T_only','B_only']
STATES=['CNV_ABNORMAL_CANDIDATE','REFERENCE_LIKE','UNCERTAIN']
def signflip_p(d):
    d=np.asarray(d,float)
    if len(d)<5:return np.nan
    if len(d)>20:raise ValueError('Exact sign flip limit exceeded')
    observed=abs(d.mean());hits=0;total=0
    for s in itertools.product([-1,1],repeat=len(d)):
        hits+=abs(np.mean(d*np.array(s)))>=observed-1e-12;total+=1
    return hits/total
def bootstrap_ci(d):
    rng=np.random.default_rng(20260925)
    d=np.asarray(d,float)
    if len(d)<5:return [np.nan,np.nan]
    v=np.mean(d[rng.integers(0,len(d),size=(10000,len(d)))],axis=1)
    return np.quantile(v,[.025,.975]).tolist()
def main():
    idx=pd.read_csv(R/'private/units_index.tsv',sep='\t');donor=[];cells=[];calibration=[];chr_records=[];statuses=[]
    for _,entry in idx.iterrows():
        p=Path(entry.path)
        if not (p/'INFERENCE_DONE').exists():
            statuses.append({'cohort':entry.cohort,'status':'NOT_EVALUABLE' if not entry.ready else 'FAILED'});continue
        m=pd.read_csv(p/'metadata.tsv',sep='\t').set_index('cell');states=[];signs=[];chrom=[]
        for v in VARIANTS:
            s=pd.read_csv(p/f'scores_{v}.tsv',sep='\t').set_index('cell').reindex(m.index)
            assert s.state.notna().all();states.append(s.state)
            ch=pd.read_csv(p/f'chromosome_means_{v}.tsv',sep='\t').set_index('cell').reindex(m.index)
            th=pd.read_csv(p/f'chr_thresholds_{v}.tsv',sep='\t').set_index('chr').threshold99
            assert set(ch.columns)==set(th.index) and 'chr3' not in ch.columns
            ch=ch.reindex(columns=sorted(ch.columns));t=th.reindex(ch.columns).clip(lower=1e-8)
            signs.append(np.where(ch.abs().gt(t),np.sign(ch),0));chrom.append(ch)
        a=np.stack(signs);consistent=((a[0]!=0)&(a[0]==a[1])&(a[1]==a[2])).sum(axis=1)
        ss=pd.concat(states,axis=1)
        m['consensus_state']=np.where(ss.eq(STATES[0]).all(axis=1)&(consistent>=2),STATES[0],np.where(ss.eq(STATES[1]).all(axis=1),STATES[1],STATES[2]))
        m['cohort']=entry.cohort;m['unit']=entry.unit
        m.to_csv(p/'consensus_private.tsv',sep='\t')
        statuses.append({'cohort':entry.cohort,'status':'DONE'})
        for role,g in m.groupby('cnv_label'):
            for state,gg in g.groupby('consensus_state'):
                calibration.append({'cohort':entry.cohort,'unit':entry.unit,'role':role,'state':state,'n_cells':len(gg)})
        my=m[m.cnv_label=='Myeloid'];cells.append(my)
        for state,g in my.groupby('consensus_state'):
            donor.append({'cohort':entry.cohort,'unit':entry.unit,'state':state,'n_cells':len(g),'mean_expression':g.slc6a6_log1p10k.mean(),'detection':g.slc6a6_detected.mean()})
            if len(g)>=20:
                for c in chrom[0].columns:chr_records.append({'cohort':entry.cohort,'unit':entry.unit,'state':state,'chr':c,'mean_residual':chrom[0].loc[g.index,c].mean()})
    dd=pd.DataFrame(donor);dd.to_csv(R/'private/donor_group_expression.tsv',sep='\t',index=False)
    if not len(dd):raise RuntimeError('No inferCNV unit completed;no expression comparison created')
    cal=pd.DataFrame(calibration);cal.groupby(['cohort','role','state']).agg(n_cells=('n_cells','sum'),n_units=('unit','nunique')).reset_index().to_csv(R/'public/control_calibration.tsv',sep='\t',index=False)
    summaries=[];tests=[]
    for c in ['GSE263733','GSE278688','GSE242230']:
        e=dd[(dd.cohort==c)&(dd.n_cells>=20)]
        for state in STATES:
            allg=dd[(dd.cohort==c)&(dd.state==state)];g=e[e.state==state]
            summaries.append({'cohort':c,'state':state,'n_cells':int(allg.n_cells.sum()),'n_units_any':len(allg),'n_units_ge20cells':len(g),'donor_equal_mean_expression':g.mean_expression.mean() if len(g)>=3 else np.nan,'donor_equal_detection':g.detection.mean() if len(g)>=3 else np.nan,'status':'DONE' if len(g)>=3 else 'NOT_EVALUABLE','reason':'At least3 author units for public aggregate expression;not malignancy labels'})
        p=e.pivot(index='unit',columns='state',values='mean_expression').reindex(columns=STATES[:2]).dropna()
        delta=p[STATES[0]]-p[STATES[1]];n=len(delta);ci=bootstrap_ci(delta)
        tests.append({'cancer':'PDAC','cohort':c,'stage_id':'06_EXTERNAL','run_id':R.name,'analysis_version':'myeloid_infercnv_v1','gene':'SLC6A6','comparison':'CNV_ABNORMAL_CANDIDATE_minus_REFERENCE_LIKE_myeloid','n_paired_author_units':n,'min_cells_per_group_unit':20,'mean_delta_log1p10k':float(delta.mean()) if n>=5 else np.nan,'ci_lower':ci[0],'ci_upper':ci[1],'p_value':signflip_p(delta),'q_value':np.nan,'test_family':'THREE_COHORT_PRIMARY_SLC6A6','n_planned_tests':3,'status':'DONE' if n>=5 else 'NOT_EVALUABLE','reason':'Paired author-unit exact sign-flip;min5 pairs;CNV grouping excludes chr3;reference sensitivity consensus' if n>=5 else 'Fewer than5 author units with at least20 cells in both groups;no cell-level pseudoreplicated P'})
    tests=pd.DataFrame(tests);ok=tests.p_value.notna();pv=tests.loc[ok,'p_value'];order=np.argsort(pv.to_numpy());q=np.empty(len(pv))
    if len(pv):q[order]=np.minimum.accumulate((pv.to_numpy()[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1].clip(0,1);tests.loc[ok,'q_value']=q
    tests['n_evaluable_tests']=int(ok.sum());tests.to_csv(R/'public/SLC6A6_comparison.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(summaries).to_csv(R/'public/group_summary.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(statuses).groupby(['cohort','status']).size().reset_index(name='n_units').to_csv(R/'public/run_coverage.tsv',sep='\t',index=False)
    ch=pd.DataFrame(chr_records)
    if len(ch):
        ch=ch.groupby(['cohort','state','chr']).agg(n_units=('unit','nunique'),mean_residual=('mean_residual','mean')).reset_index();ch.loc[ch.n_units<3,'mean_residual']=np.nan
        ch.to_csv(R/'public/aggregate_chr_profiles.tsv',sep='\t',index=False,na_rep='NA')
    (R/'SUMMARY_DONE.json').write_text(json.dumps({'tests_evaluable':int(ok.sum()),'new_malignancy_diagnoses':0,'individual_data_exported':False}))
if __name__=='__main__':main()
