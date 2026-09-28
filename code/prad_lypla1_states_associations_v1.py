"""Author-defined malignant states and donor-wise LYPLA1 coexpression. Server only inputs."""
import argparse, hashlib, json, platform
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def save(df,path): df.to_csv(path,sep='\t',index=False,na_rep='NA')

def residual_rank(x,cov):
    z=stats.rankdata(x,axis=0).astype(float)
    return z-cov@np.linalg.lstsq(cov,z,rcond=None)[0]

def correlate(x,y,cov):
    a=residual_rank(x,cov); b=residual_rank(y,cov)
    den=np.sqrt((a*a).sum(0)*(b*b).sum())
    return np.divide(b@a,den,out=np.full(x.shape[1],np.nan),where=den>1e-10)

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--commit',required=True);args=p.parse_args()
    r=args.out;r.mkdir(parents=True,exist_ok=True)
    with (r/'.running').open('x') as f:f.write(r.name)
    pub=r/'public';priv=r/'private';pub.mkdir();priv.mkdir()
    spec=dict(version='prad_lypla1_states_associations_v1',parent_commit=args.commit,source=str(args.source),gene='LYPLA1',target_id='ENSG00000120992',selection='type=cancer AND celltype_major_v2=Epithelial AND malignant_anno_merged=malignant',annotation='author celltype_subset_v2 reused; no new clustering',normalization='log1p(10000*raw_count/raw_library_total)',donor_min_cells=30,gene_min_detected_cells_per_donor=10,gene_min_fraction_per_donor=.10,min_donors_for_inference=8,state_min_cells_per_donor_each_arm=10,state_min_donors_for_descriptive=3,association_primary='within-donor Spearman rank residual correlation controlling rank(log1p totalUMI) and rank(percent.mito)',association_sensitivity=['ordinary within-donor Spearman','add author subset fixed effects within donor','LYPLA1-positive cells only with same technical adjustment'],association_test='two-sided exact binomial sign test across nonzero donor correlations; no cell-level pooled P',test_family='all measured non-LYPLA1 genes, separately per mode; BH retained but nominal P<0.05 prioritization',state_test='exact binomial sign test of donor subgroup mean minus remaining malignant cells; >=8 donors',state_contrast='within same donor, subset versus all other malignant subsets; not tumor vs normal',software=dict(python=platform.python_version(),anndata=ad.__version__,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),limitations=['donor identity inherited, no genotype files available for swap validation','coexpression is not regulation or pathway activity','state covariate adjustment can remove genuine between-state biology','sign-test P tests direction consistency, not correlation magnitude','all nominal discoveries exploratory; no independent validation in this script'])
    (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
    a=ad.read_h5ad(args.source,backed='r'); obs=a.obs
    mask=(obs['type'].astype(str)=='cancer')&(obs.celltype_major_v2.astype(str)=='Epithelial')&(obs.malignant_anno_merged.astype(str)=='malignant')
    o=obs.loc[mask].copy();assert o.index.is_unique and o.donor_id.notna().all()
    x=a.raw.X[np.flatnonzero(mask)].tocsr();v=a.raw.var.copy();xy=np.asarray(a.obsm['X_umap'])[mask]
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.round(x.data)).all()
    names=v.feature_name.astype(str).to_numpy();ids=v.index.astype(str).to_numpy();j=np.flatnonzero(ids=='ENSG00000120992');assert len(j)==1 and names[j[0]]=='LYPLA1';j=j[0]
    lib=np.asarray(x.sum(1)).ravel();assert (lib>0).all()
    y=np.log1p(10000*x[:,j].toarray().ravel()/lib);o['LYPLA1']=y;o['detected']=y>0;o['library']=lib
    o['state']=o.celltype_subset_v2.astype(str);o['donor']=o.donor_id.astype(str)
    assert o['state'].notna().all() and np.isfinite(o['percent.mito'].astype(float)).all()
    save(o.reset_index(),priv/'malignant_cells.tsv.gz')
    # Read-only genotype QC unavailable: retain depositor sex labels, do not infer donor identity.
    qc=o.groupby('donor',observed=True).agg(n_cells=('LYPLA1','size'),n_positive=('detected','sum'),mean=('LYPLA1','mean'),n_states=('state','nunique'))
    save(qc.reset_index(),priv/'donor_qc.tsv')
    state_rows=[];paired=[]
    for s,g in o.groupby('state',observed=True):
        dd=[];de=[];pp=[];sg=[]
        for d,t in o.groupby('donor',observed=True):
            u=t[t.state==s];w=t[t.state!=s]
            if len(u)>=10: sg.append((u.LYPLA1.mean(),u.detected.mean(),u.loc[u.detected,'LYPLA1'].mean()))
            if len(u)>=10 and len(w)>=10:
                delta=u.LYPLA1.mean()-w.LYPLA1.mean();dd.append(delta);de.append(u.detected.mean()-w.detected.mean())
                paired.append(dict(state=s,donor=d,n_subset=len(u),n_other=len(w),delta_mean=delta))
        nz=np.array(dd);nz=nz[nz!=0];n=len(dd)
        state_rows.append(dict(state=s,n_cells=len(g),n_donors_any=g.donor.nunique(),largest_donor_fraction=g.donor.value_counts().max()/len(g),pooled_detection=g.detected.mean(),pooled_mean=g.LYPLA1.mean(),positive_mean=g.loc[g.detected,'LYPLA1'].mean(),n_donors_ge10=len(sg),donor_equal_mean=np.mean([q[0] for q in sg]) if len(sg)>=3 else np.nan,donor_equal_detection=np.mean([q[1] for q in sg]) if len(sg)>=3 else np.nan,n_paired_donors=n,n_higher=int((np.array(dd)>0).sum()),mean_paired_delta=np.mean(dd) if n else np.nan,p_value=stats.binomtest(int((nz>0).sum()),len(nz)).pvalue if n>=8 and len(nz) else np.nan,status='DONE' if n>=8 else 'NOT_EVALUABLE',reason='' if n>=8 else 'fewer_than8_donors_with10cells_in_both_arms'))
    st=pd.DataFrame(state_rows);st['q_value']=np.nan;ok=st.p_value.notna()
    if ok.any():st.loc[ok,'q_value']=multipletests(st.loc[ok,'p_value'],method='fdr_bh')[1]
    save(st,pub/'state_summary.tsv');save(pd.DataFrame(paired),priv/'state_donor_comparisons.tsv')
    # Matrix rows remain exact obs order. Rank genes in blocks to bound memory.
    modes=['ordinary','technical_adjusted','state_adjusted','positive_only'];arrays={m:[] for m in modes};dn=[];checks=[]
    for d in sorted(o.donor.unique()):
        ix=np.flatnonzero(o.donor.to_numpy()==d)
        if len(ix)<30 or (y[ix]>0).sum()<10 or np.ptp(y[ix])==0:continue
        dn.append(d);xx=x[ix].toarray();z=np.log1p(10000*xx/lib[ix,None]);yy=y[ix]
        cov=np.column_stack([np.ones(len(ix)),stats.rankdata(np.log1p(lib[ix])),stats.rankdata(o['percent.mito'].to_numpy(dtype=float)[ix])]);cov[:,1:]=(cov[:,1:]-cov[:,1:].mean(0))/(cov[:,1:].std(0)+1e-12)
        dummy=pd.get_dummies(o.state.iloc[ix],drop_first=True,dtype=float).to_numpy();cvstate=np.column_stack([cov,dummy]);pos=yy>0
        eligibility=(xx>0).sum(0)>=max(10,int(np.ceil(.1*len(ix))))
        per={m:np.full(len(v),np.nan) for m in modes}
        for k in range(0,len(v),512):
            sl=slice(k,min(k+512,len(v))); zz=z[:,sl]
            per['ordinary'][sl]=correlate(zz,yy,np.ones((len(ix),1)))
            per['technical_adjusted'][sl]=correlate(zz,yy,cov)
            if len(ix)-np.linalg.matrix_rank(cvstate)>=15: per['state_adjusted'][sl]=correlate(zz,yy,cvstate)
            if pos.sum()>=30:
                pe=(xx[pos,sl]>0).sum(0)>=max(10,int(np.ceil(.1*pos.sum())))
                q=correlate(zz[pos],yy[pos],cov[pos]);q[~pe]=np.nan;per['positive_only'][sl]=q
        for m in modes:
            per[m][~eligibility]=np.nan;per[m][j]=np.nan;arrays[m].append(per[m])
        if len(checks)<3:
            kk=np.flatnonzero(eligibility & (np.arange(len(v))!=j))[:3]
            checks.extend(abs(stats.spearmanr(z[:,k],yy).statistic-per['ordinary'][k]) for k in kk)
        print('processed donor',len(dn),'cells',len(ix),flush=True)
    rows=[]
    for mode,ls in arrays.items():
        mat=np.asarray(ls);save(pd.DataFrame(mat,index=dn,columns=ids).rename_axis('donor').reset_index(),priv/(mode+'_donor_rho.tsv.gz'))
        for k in range(len(v)):
            if k==j:continue
            rr=mat[:,k];rr=rr[np.isfinite(rr)];n=len(rr);nz=rr[rr!=0];status='DONE' if n>=8 and len(nz) else 'NOT_EVALUABLE'
            rows.append(dict(mode=mode,gene=names[k],stable_gene_id=ids[k],n_donors=n,n_positive=int((rr>0).sum()),n_negative=int((rr<0).sum()),mean_rho=rr.mean() if n else np.nan,median_rho=np.median(rr) if n else np.nan,positive_fraction=(rr>0).mean() if n else np.nan,p_value=stats.binomtest(int((nz>0).sum()),len(nz)).pvalue if status=='DONE' else np.nan,status=status,reason='' if status=='DONE' else 'fewer_than8_evaluable_donors',technical_gene=bool(names[k].startswith(('MT-','RPL','RPS')))))
    res=pd.DataFrame(rows);res['q_value']=np.nan
    for mode in modes:
        ok=(res['mode']==mode)&res.p_value.notna();res.loc[ok,'q_value']=multipletests(res.loc[ok,'p_value'],method='fdr_bh')[1]
    save(res,pub/'gene_associations_all.tsv.gz')
    primary=res[res['mode']=='technical_adjusted'].copy()
    for mode in ['ordinary','state_adjusted','positive_only']:
        z=res[res['mode']==mode][['stable_gene_id','n_donors','mean_rho','p_value','positive_fraction']].rename(columns={c:mode+'_'+c for c in ['n_donors','mean_rho','p_value','positive_fraction']})
        primary=primary.merge(z,on='stable_gene_id',validate='one_to_one')
    primary=primary.sort_values(['p_value','mean_rho'],ascending=[True,False])
    pre=dict(cancer='PRAD',cohort='PRAD24_CELLxGENE',stage_id='06_EXTERNAL',run_id=r.name,analysis_version=spec['version'],analysis_type='within_donor_rank_association',metabolite_key='NA',metabolite_name='NA',gene=primary.gene,unit='author_donor',n=primary.n_donors,n_reference=24,effect_type='mean_partial_spearman',effect=primary.mean_rho,ci_lower=np.nan,ci_upper=np.nan,p_value=primary.p_value,q_value=primary.q_value,test_family='all_evaluable_non_LYPLA1_genes_primary',family_n_evaluable=int(primary.p_value.notna().sum()),status=primary.status,reason=primary.reason,source_id='68b23fda-7191-46a5-8870-819feca3e66e')
    formatted=pd.DataFrame(pre,index=primary.index);formatted=pd.concat([formatted,primary.drop(columns=[c for c in formatted if c in primary])],axis=1);save(formatted,pub/'results.tsv.gz')
    top=primary[(primary.p_value<.05)&(primary.mean_rho>0)&(~primary.technical_gene)].head(30);save(top,pub/'positive_associations_priority.tsv')
    eligible=st[st.n_donors_ge10>=3].sort_values('donor_equal_mean')
    fig,ax=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[1.2,1]})
    ax[0].barh(eligible.state,eligible.donor_equal_mean,color='#ac3154');ax[0].set_xlabel('Donor-equal mean log1p(CP10K)');ax[0].set_title('Author malignant subsets\n>=3 donors with >=10 cells each')
    tt=top.head(15).iloc[::-1];ax[1].barh(tt.gene,tt.mean_rho,color='#318d88');ax[1].set_xlabel('Mean within-donor adjusted rank correlation');ax[1].set_title('Positive companions | nominal sign-test P<0.05\nTechnical genes excluded from display only')
    fig.suptitle('PRAD LYPLA1 | malignant-cell states and within-donor associations');fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(pub/('LYPLA1_states_associations.'+ext),dpi=240,bbox_inches='tight')
    validation=dict(status='DONE',n_malignant_cells=len(o),n_donors=o.donor.nunique(),n_association_donors=len(dn),n_states=len(st),n_states_inference=int(st.p_value.notna().sum()),n_primary_genes_evaluable=int(primary.p_value.notna().sum()),n_nominal_positive=int(((primary.p_value<.05)&(primary.mean_rho>0)).sum()),n_nominal_negative=int(((primary.p_value<.05)&(primary.mean_rho<0)).sum()),n_primary_q_lt005=int((primary.q_value<.05).sum()),scalar_spearman_max_error=float(max(checks)),exact_metadata_alignment=True,integer_raw_counts=True,no_reclustering=True,genotype_swap_check='NOT_EVALUABLE: genotype/BAM inputs absent; author donor IDs inherited',sex_labels=obs.sex.astype(str).value_counts().to_dict(),cross_cohort_status='NOT_RUN by this script')
    assert validation['scalar_spearman_max_error']<1e-10
    (pub/'validation.json').write_text(json.dumps(validation,indent=2));save(pd.DataFrame([dict(path=str(p),sha256=sha(p)) for p in [args.source,Path(__file__)]]),pub/'source_manifest.tsv')
    a.file.close();print(json.dumps(validation),flush=True)

if __name__=='__main__': main()
