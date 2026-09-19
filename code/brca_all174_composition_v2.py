"""Exploratory all-relation sensitivity to MCPcounter-derived composition proxies.

Modified panels exclude all117 candidates. Scores are NOT cellular proportions.
Official MCPcounter calculation independently checked in base R before inference.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,csv,json,hashlib,subprocess,platform
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def design(er,cov):
    if cov.shape[1]:
        cv=np.column_stack([stats.rankdata(cov[:,i]) for i in range(cov.shape[1])])
        v=np.column_stack([er,cv])
    else:v=er[:,None]
    sd=v.std(0)
    if (sd==0).any():raise ValueError('constant_covariate')
    v=(v-v.mean(0))/sd
    return np.column_stack([np.ones(len(er)),v])


def estimate(x,y,er,cov):
    z=design(er,cov)
    if np.linalg.matrix_rank(z)<z.shape[1] or np.linalg.cond(z)>30:raise ValueError('ill_conditioned_design')
    q=np.linalg.qr(z,mode='reduced')[0]
    xx=stats.rankdata(x);yy=stats.rankdata(y)
    ex=xx-q@(q.T@xx);ey=yy-q@(q.T@yy)
    norm=np.linalg.norm(ex)*np.linalg.norm(ey)
    if norm<1e-10:raise ValueError('constant_residual')
    return float(ex@ey/norm),q,ex,ey,float(np.linalg.cond(z))


def permutation(rho,q,ex,ey,er,rng):
    groups=[np.flatnonzero(er==g) for g in np.unique(er)];exceed=0;nperm=9999
    for start in range(0,nperm,500):
        count=min(500,nperm-start);per=np.tile(ey,(count,1))
        for ix in groups:per[:,ix]=ey[ix][np.argsort(rng.random((count,len(ix))),axis=1)]
        residual_norm2=np.sum(per*per,axis=1)-np.sum((per@q)**2,axis=1)
        values=(per@ex)/(np.linalg.norm(ex)*np.sqrt(np.maximum(residual_norm2,1e-20)))
        exceed+=int(np.sum(abs(values)>=abs(rho)-1e-12))
    return (exceed+1)/(nperm+1)


def bootstrap(x,y,er,cov,rng):
    groups=[np.flatnonzero(er==g) for g in np.unique(er)];res=[]
    for k in range(500):
        ix=np.concatenate([rng.choice(g,size=len(g),replace=True) for g in groups])
        try:res.append(estimate(x[ix],y[ix],er[ix],cov[ix])[0])
        except ValueError:pass
    return (np.quantile(res,[.025,.975]).tolist() if len(res)>=400 else [np.nan,np.nan]),len(res)


def main(root,project):
    out=root/'public';private=root/'private'
    old=project/'results/collaborative/BRCA/A/20260919_ER_availability_v3/aggregate/er_availability_all_174.tsv'
    oldgeo=project/'results/collaborative/BRCA/A/20260919_followup_ER_v2/source/geo_metadata_lines.txt'
    src=project/'data/candidates/camp_primary_tissue_multicancer';mp=src/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    m=pd.read_csv(mp,dtype=str);m=m[m.Dataset.eq('BRCA1')&m.TN.eq('Tumor')].copy()
    assert len(m)==61 and all(m[c].is_unique for c in ['RNAID','MetabID','CommonID'])
    rf=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
    mf=src/'processed_metabolomics'/m.MetabFile.iloc[0]
    rna=pd.read_csv(rf,index_col=0);met=pd.read_excel(mf,sheet_name=m.MetabFile_sheet.iloc[0],index_col=0)
    obs=pd.read_excel(mf,sheet_name='data',index_col=0)
    for d in [rna,met,obs]:d.columns=d.columns.astype(str);assert d.index.is_unique and d.columns.is_unique
    assert m.RNAID.isin(rna.columns).all() and m.MetabID.isin(met.columns).all()
    prior=pd.read_csv(old,sep='\t');assert len(prior)==174 and prior.relation_id.is_unique
    poolpath=project/'results/collaborative/BRCA/A/20260919T140105Z_scRNA117_v1/source/candidate_comparison_v1.tsv'
    pool=set(pd.read_csv(poolpath,sep='\t').gene);assert len(pool)==117
    lines=list(csv.reader(oldgeo.read_text().splitlines(),delimiter='\t'));md=pd.DataFrame(index=lines[0][1:])
    assert md.index.is_unique
    for row in lines[1:]:md[row[1].split(': ',1)[0]]=[v.split(': ',1)[1] for v in row[1:]]
    gsm=m.RNAID.str.replace('.CEL.gz','',regex=False);assert gsm.isin(md.index).all()
    assert md.loc[gsm,'tissue type'].eq('Tumor').all()
    er=md.loc[gsm,'estrogen receptor status'].map({'Positive':1.,'Negative':0.}).to_numpy()
    assert np.isfinite(er).all()
    markerpath=root/'source/MCPcounter_genes.txt';markers=pd.read_csv(markerpath,sep='\t')
    original={};modified={};coverage=[];markerrows=[]
    for population,g in markers.groupby('Cell population',sort=True):
        full=sorted(set(g['HUGO symbols'].dropna()));present=[x for x in full if x in rna.index]
        kept=[x for x in present if x not in pool]
        original[population]=rna.loc[present,m.RNAID].mean(axis=0).to_numpy()
        modified[population]=rna.loc[kept,m.RNAID].mean(axis=0).to_numpy() if kept else np.full(len(m),np.nan)
        ok=len(kept)>=2 and len(kept)/len(full)>=.5
        coverage.append({'population':population,'original_unique_markers':len(full),'measured':len(present),
            'excluded_candidates':';'.join(sorted(set(present)&pool)),'retained':len(kept),
            'coverage_fraction':len(kept)/len(full),'eligible_PCA':ok,
            'original_modified_score_rho':stats.spearmanr(original[population],modified[population]).statistic,
            'reason':'' if ok else 'single_marker_or_insufficient_coverage'})
        for gene in full:markerrows.append({'population':population,'gene':gene,'measured':gene in present,
            'excluded_candidate':gene in pool,'used_modified_score':gene in kept})
    original=pd.DataFrame(original,index=m.RNAID);modified=pd.DataFrame(modified,index=m.RNAID)
    assert np.isfinite(original.values).all() and np.isfinite(modified.values).all()
    save(original.reset_index(),private/'MCPcounter_original_scores.tsv');save(modified.reset_index(),private/'MCPcounter_modified_scores.tsv')
    cover=pd.DataFrame(coverage);save(cover,out/'composition_marker_coverage.tsv');save(pd.DataFrame(markerrows),out/'composition_marker_decisions.tsv')
    # Compare every original panel/sample with official R implementation, no package install.
    rscript=root/'scripts/check_official_MCPcounter.R'
    rscript.write_text('''args <- commandArgs(trailingOnly=TRUE)
source(args[1])
x <- as.matrix(read.csv(args[2],row.names=1,check.names=FALSE))
g <- read.table(args[3],sep="\\t",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
s <- MCPcounter.estimate(x,featuresType="HUGO_symbols",genes=g,probesets=data.frame())
write.table(t(s),args[4],sep="\\t",quote=FALSE,col.names=NA)
''')
    subprocess.run(['Rscript',str(rscript),str(root/'source/MCPcounter_original.R'),str(rf),str(markerpath),str(private/'MCPcounter_R_check.tsv')],check=True)
    rc=pd.read_csv(private/'MCPcounter_R_check.tsv',sep='\t',index_col=0).loc[original.index,original.columns]
    rerr=float(np.max(abs(rc.values-original.values)));assert rerr<1e-10
    eligible=cover.loc[cover.eligible_PCA,'population'].tolist();assert len(eligible)>=8
    zz=modified[eligible].to_numpy();zz=(zz-zz.mean(0))/zz.std(0,ddof=1)
    u,s,v=np.linalg.svd(zz,full_matrices=False);pcs=u[:,:2]*s[:2]
    ev=s*s/(s*s).sum()
    loadings=pd.DataFrame({'population':eligible,'PC1_loading':v[0],'PC2_loading':v[1]});save(loadings,out/'composition_PCA_loadings.tsv')
    save(pd.DataFrame(zz,columns=eligible).corr().rename_axis('population').reset_index(),out/'composition_score_correlations.tsv')
    models={'ER_PC12':pcs,'ER_MonoEndoFib':modified[['Monocytic lineage','Endothelial cells','Fibroblasts']].values}
    result=[];max_baseline_error=0.;checks=[]
    for i,r in prior.iterrows():
        base={'cancer':'BRCA','cohort':'CAMP_BRCA1','stage_id':'04_ROBUSTNESS','run_id':root.name,
            'analysis_version':'all117_robustness_v2','analysis_type':'partial_rank_composition_proxy_sensitivity',
            'metabolite_key':r.relation_id.split('|')[0],'metabolite_name':r.metabolite,'gene':r.gene,
            'unit':'author_mapped_tumor_specimen_NOT_verified_independent_patient',
            'relation_id':r.relation_id,'n':np.nan,'n_reference':np.nan,'effect_type':'partial_rank_correlation',
            'effect':np.nan,'ci_lower':np.nan,'ci_upper':np.nan,'p_value':np.nan,'q_value':np.nan,
            'test_family':'','family_n_evaluable':np.nan,'status':'NOT_EVALUABLE','reason':'',
            'source_id':'GSE37751_CAMP_processed','ER_v3_rho':r.get('partial_rank_rho',np.nan),
            'ER_v3_p':r.get('p_value',np.nan),'ER_v3_q':r.get('q_value',np.nan),
            'original_q':r.get('original_q',np.nan),'signed_delta_from_ER':np.nan,'abs_effect_change':np.nan}
        if pd.isna(r.get('p_value')) or r.gene not in rna.index or r.metabolite not in met.index:
            for model in models:result.append(dict(base,test_family=model,reason='prior_unevaluable_or_missing_feature'))
            continue
        x=pd.to_numeric(met.loc[r.metabolite,m.MetabID],errors='coerce').to_numpy(float)
        y=pd.to_numeric(rna.loc[r.gene,m.RNAID],errors='coerce').to_numpy(float)
        available=pd.to_numeric(obs.loc[r.metabolite,m.MetabID],errors='coerce').to_numpy(float)
        mask=np.isfinite(x)&np.isfinite(y)&np.isfinite(available)&np.isfinite(er)
        xx=x[mask];yy=y[mask];ee=er[mask];base['n']=len(xx)
        if len(xx)<20 or min((ee==0).sum(),(ee==1).sum())<4:
            for model in models:result.append(dict(base,test_family=model,reason='coverage_n20_ERgroups4'))
            continue
        baseline=estimate(xx,yy,ee,np.empty((len(xx),0)))[0]
        err=abs(baseline-r.partial_rank_rho);max_baseline_error=max(max_baseline_error,err)
        assert len(xx)==r.n and err<1e-10,'Same-subset ER baseline differs from v3'
        for model,c in models.items():
            cc=c[mask];b=dict(base,test_family=model)
            if len(xx)-cc.shape[1]-2<10:b['reason']='residual_df';result.append(b);continue
            try:rho,q,ex,ey,cond=estimate(xx,yy,ee,cc)
            except ValueError as e:b['reason']=str(e);result.append(b);continue
            seed=int(hashlib.sha256((model+'|'+r.relation_id+'|20260919').encode()).hexdigest()[:8],16)
            pp=permutation(rho,q,ex,ey,ee,np.random.default_rng(seed))
            ci,bn=bootstrap(xx,yy,ee,cc,np.random.default_rng(seed+1))
            b.update(effect=rho,ci_lower=ci[0],ci_upper=ci[1],p_value=pp,status='DONE',reason='',
                     signed_delta_from_ER=rho-baseline,abs_effect_change=abs(rho)-abs(baseline),
                     design_condition_number=cond,seed=seed,bootstrap_valid=bn,
                     interval_scope='pointwise_stratified_specimen_bootstrap;covariate_scores_fixed;independence_unverified')
            # Independent residual calculation via least squares.
            z=design(ee,cc);rx=stats.rankdata(xx);ry=stats.rankdata(yy)
            ax=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];ay=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0]
            checks.append(abs(rho-np.corrcoef(ax,ay)[0,1]));result.append(b)
        if i%20==0:print('relations',i+1,'/174',flush=True)
    d=pd.DataFrame(result)
    for model,ix in d.groupby('test_family').groups.items():
        ok=d.index.isin(ix)&d.p_value.notna();d.loc[ix,'family_n_evaluable']=int(ok.sum())
        if ok.any():
            ps=d.loc[ok,'p_value'].to_numpy();qs=multipletests(ps,method='fdr_bh')[1];d.loc[ok,'q_value']=qs
            order=np.argsort(ps);other=np.minimum(1,np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]);assert np.max(abs(qs[order]-other))<1e-12
    assert len(d)==348 and d.groupby('test_family').relation_id.nunique().eq(174).all()
    assert max(checks)<1e-10
    save(d,out/'all174_composition_sensitivity.tsv')
    paths=[mp,rf,mf,old,oldgeo,poolpath,markerpath,root/'source/MCPcounter_original.R',Path(__file__)]
    save(pd.DataFrame([{'path':str(f),'sha256':sha(f)} for f in paths]),out/'composition_input_manifest.tsv')
    result={'status':'DONE_EXPLORATORY','n_specimens':61,'independent_patient_identity':'NOT_VERIFIED',
        'rna_platform':'processed_gene_symbol_expression_microarray;no new log transform',
        'official_R_max_error':rerr,'n_modified_PCA_panels':len(eligible),'candidate_marker_exclusions':['KYNU','HAL'],
        'PCA_explained_variance_PC1':float(ev[0]),'PCA_explained_variance_PC2':float(ev[1]),
        'baseline_v3_max_rho_error':max_baseline_error,'partial_rank_independent_formula_max_error':max(checks),
        'scope':'composition-related marker-score sensitivity,NOT corrected true fractions or causal mediation',
        'families':{k:{'planned':len(g),'evaluable':int(g.p_value.notna().sum()),'q_lt_005':int(g.q_value.lt(.05).sum())} for k,g in d.groupby('test_family')},
        'score_method_modified_not_independently_validated':True,'source_sha256':sha(Path(__file__)),
        'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}
    (out/'composition_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--project',type=Path,required=True)
    a=p.parse_args();main(a.root,a.project)
