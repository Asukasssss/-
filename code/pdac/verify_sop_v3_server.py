"""Independent checks of v3 patient or donor summary results, entirely on server."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

def bh_reference(p):
    p=np.asarray(p,float);idx=np.flatnonzero(np.isfinite(p));out=np.full(len(p),np.nan)
    # Direct definition, independent of producer's cumulative-min implementation.
    order=idx[np.argsort(p[idx])];m=len(order)
    for rank,i in enumerate(order):out[i]=min(1,min(m*p[j]/(k+1) for k,j in enumerate(order) if k>=rank))
    return out
def close(a,b,tol=1e-10):np.testing.assert_allclose(a,b,rtol=0,atol=tol,equal_nan=True)
def internal(run):
    root=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
    rp=root/'data/candidates/camp_primary_tissue_multicancer/gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed/GSE62452.hugene10st.gene_symbol.csv'
    rna=pd.read_csv(rp,index_col=0);pairs=pd.read_csv(run/'private/RNA_pairs_private.tsv',sep='\t')
    assert len(pairs)==11 and pairs.patient_id.is_unique
    files=['paired_RNA.tsv','paired_RNA_history_supplement.tsv'];allrows=[];checked=0
    for name in files:
        d=pd.read_csv(run/'public'/name,sep='\t');assert d.stable_gene_id.is_unique
        for _,r in d.iterrows():
            if pd.isna(r.rna_label):assert pd.isna(r.p_value);continue
            x=rna.loc[r.rna_label,pairs.tumor_RNA_id].to_numpy(float);y=rna.loc[r.rna_label,pairs.normal_RNA_id].to_numpy(float);mask=np.isfinite(x)&np.isfinite(y);x=x[mask];y=y[mask];delta=x-y
            assert len(delta)==r.n and (delta>0).sum()==r.n_up and (delta<0).sum()==r.n_down and (delta==0).sum()==r.n_equal
            close(r.effect,delta.mean());close(r.median_delta,np.median(delta))
            if len(delta)>=8 and np.std(delta,ddof=1)>0:
                test=stats.ttest_rel(x,y);close(r.p_value,test.pvalue,1e-12);se=stats.sem(delta);ci=stats.t.interval(.95,len(delta)-1,loc=np.mean(delta),scale=se);close([r.ci_lower,r.ci_upper],ci)
                # Independent RNG replay of the frozen whole-pair bootstrap.
                boot=delta[np.random.default_rng(int(r.seed)).integers(len(delta),size=(4000,len(delta)))].mean(axis=1);close([r.bootstrap_mean_lower,r.bootstrap_mean_upper],np.quantile(boot,[.025,.975]));checked+=1
            else:assert pd.isna(r.p_value)
        close(d.q_value,bh_reference(d.p_value),1e-12);allrows.append(d)
    df=pd.concat(allrows);assert len(df)==687 and df.stable_gene_id.is_unique
    view=pd.read_csv(run/'public/RNA_P005_view.tsv',sep='\t');assert set(view.stable_gene_id)==set(allrows[0].loc[allrows[0].p_value<.05,'stable_gene_id'])
    for name in ['sample_identity_audit_private.tsv','metabolite_pairs_private.tsv','RNA_pairs_private.tsv','tumor_multiomics_map_private.tsv']:assert not (run/'public'/name).exists()
    return {'status':'PASS','independent_paired_t_and_CI_rows':checked,'independent_bootstrap_replays':checked,'independent_BH_families':2,'current_gene_P005_view_count':len(view),'all687_and_private_separation':True,'not_verified':['Clinical identity recertification','same aliquot','external replication']}
def source(run):
    prof=pd.read_csv(run/'public/sc_celltype_profiles.tsv',sep='\t');top=pd.read_csv(run/'public/sc_source_stability.tsv',sep='\t');cross=pd.read_csv(run/'public/sc_cross_study.tsv',sep='\t');rows=0;topn=0
    for cohort in sorted(prof.cohort.unique()):
        pb=pd.read_csv(run/'private'/(cohort+'_sc_donor_profiles_private.tsv'),sep='\t');assert not pb.duplicated(['donor_id','celltype','stable_gene_id']).any()
        assert pb.mean_log1p10k.dropna().ge(0).all() and pb.detection_fraction.dropna().between(0,1).all()
        for (gene,ct),d in pb.groupby(['gene','celltype']):
            r=prof[(prof.cohort==cohort)&(prof.gene==gene)&(prof.celltype==ct)];assert len(r)==1;r=r.iloc[0]
            e=d[(d.n_cells>=20)&d.mean_log1p10k.notna()]
            assert r.n==len(e) and r.n_cells_total==d.n_cells.sum() and r.n_cells_eligible==e.n_cells.sum()
            if len(e)>=3:
                close(r.effect,e.mean_log1p10k.mean());close(r.mean_detection_fraction,e.detection_fraction.mean());close([r.q25,r.median_donor_mean,r.q75],e.mean_log1p10k.quantile([.25,.5,.75]).to_numpy());rows+=1
            else:assert pd.isna(r.effect)
        for _,r in top[top.cohort==cohort].iterrows():
            p=prof[(prof.cohort==cohort)&(prof.gene==r.gene)&(prof.status=='DONE')]
            if len(p)>=2 and p.mean_detection_fraction.max()>=.01:
                expected=sorted(p.loc[np.isclose(p.effect,p.effect.max(),rtol=1e-10,atol=1e-12),'celltype'])
                assert sorted(r.top_celltype.split(';'))==expected and r.status=='DONE'
                if r.bootstrap_valid:close(p.bootstrap_first_frequency.sum(),1)
            else:assert r.status=='NOT_EVALUABLE'
            topn+=1
    for _,r in cross.iterrows():
        a=top[(top.gene==r.gene)&(top.cohort==r.study_A)].iloc[0];b=top[(top.gene==r.gene)&(top.cohort==r.study_B)].iloc[0]
        valid=a.status==b.status=='DONE' and a.tie_status==b.tie_status=='UNIQUE'
        if valid:
            same=a.top_celltype==b.top_celltype;assert bool(r.same_top_all_categories)==same
            assert bool(r.same_top_both_bootstrap_ge080)==bool(same and a.bootstrap_top_frequency>=.8 and b.bootstrap_top_frequency>=.8)
    assert len(top)==len(cross)==2061 and prof.p_value.isna().all() and prof.q_value.isna().all()
    return {'status':'PASS','independent_equal_donor_numeric_rows':rows,'independent_top_rows':topn,'cross_study_rows':len(cross),'all687_per_cohort':True,'SC_no_P_q':True,'normalization_function_tests':'Synthetic unequal cell libraries verified separately in server5test suite','not_replayed':['All raw single-cell rows independently re-extracted','All donor bootstrap draws independently replayed'],'not_verified':['Clinical cross-study patient overlap','Malignant CNV','Mechanism']}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['internal','source'],required=True);a=ap.parse_args();run=Path(__file__).resolve().parent;assert (run/'NUMERICS_DONE.json').exists()
    v=internal(run) if a.mode=='internal' else source(run);(run/'public/independent_validation.json').write_text(json.dumps(v,indent=2)+'\n');print(json.dumps(v))
if __name__=='__main__':main()
