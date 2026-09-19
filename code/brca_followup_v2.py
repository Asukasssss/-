"""Server-only BRCA followup. Original values remain immutable; export aggregates only.

Run: python3 brca_followup_v2.py --root ROOT --out NEW_RUN_DIR
Requires source/geo_metadata_lines.txt if GEO metadata retrieval succeeded.
"""
import argparse, csv, hashlib, json, platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args(); root,out=a.root,a.out; agg=out/'aggregate'
    assert (out/'.running').exists()
    assert not (agg/'robustness_all_relations.tsv').exists(), 'No overwrite'
    old=root/'results/BRCA_117_screen_20260911'; src=root/'data/candidates/camp_primary_tissue_multicancer'
    mpfile=src/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    mp=pd.read_csv(mpfile,dtype=str);mp=mp[mp.Dataset.eq('BRCA1')].copy()
    assert not mp[['CommonID','RNAID','MetabID','TN']].isna().any().any()
    assert all(mp[c].is_unique for c in ['CommonID','RNAID','MetabID'])
    tm=mp[mp.TN.eq('Tumor')].copy()
    rp=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/mp.RNAFile.iloc[0]
    xp=src/'processed_metabolomics'/mp.MetabFile.iloc[0]
    rna=pd.read_csv(rp,index_col=0); met=pd.read_excel(xp,sheet_name=tm.MetabFile_sheet.iloc[0],index_col=0)
    obs=pd.read_excel(xp,sheet_name='data',index_col=0)
    for df in [rna,met,obs]:
        df.columns=df.columns.astype(str);assert df.columns.is_unique
    paths=[mpfile,rp,xp,old/'tables/BRCA_174_patient_associations.tsv',old/'tables/BRCA_117_RNA_differences.tsv',old/'tables/BRCA_117_candidate_comparison.tsv',Path(__file__)]
    initial={str(p):sha(p) for p in paths}
    prior=pd.read_csv(paths[3],sep='\t'); expr=pd.read_csv(paths[4],sep='\t')
    checks=[]
    for family,df in list(prior.groupby('analysis_type'))+[('RNA_expression',expr)]:
        ok=df.p_value.notna();q=multipletests(df.loc[ok,'p_value'],method='fdr_bh')[1]
        err=float(np.max(abs(q-df.loc[ok,'q_value'].to_numpy())))
        assert err<1e-12
        checks.append(dict(family=family,planned=len(df),evaluable=int(ok.sum()),significant=int((q<.05).sum()),max_q_error=err))
    pd.DataFrame(checks).to_csv(agg/'multiplicity_verification.tsv',sep='\t',index=False)
    primary=prior[prior.analysis_type.eq('author_processed_primary')].copy()
    assert len(primary)==174 and not primary.relation_id.duplicated().any()
    # Exact GSM accession extraction is a documented RNA filename suffix removal, not fuzzy matching.
    tm['GSM']=tm.RNAID.str.replace('.CEL.gz','',regex=False)
    metadata_path=out/'source/geo_metadata_lines.txt'; er=np.full(len(tm),np.nan)
    meta_info={'status':'NOT_AVAILABLE','independent_patient_id':'NOT_VERIFIED','batch':'NOT_AVAILABLE','purity':'NOT_AVAILABLE'}
    if metadata_path.exists():
        rows=list(csv.reader(metadata_path.read_text().splitlines(),delimiter='\t'))
        accession=[r[1:] for r in rows if r[0]=='!Sample_geo_accession'][0]
        assert len(accession)==len(set(accession))
        d={k:{} for k in accession}
        for row in rows:
            if row[0]=='!Sample_characteristics_ch1':
                assert len(row)-1==len(accession)
                for k,v in zip(accession,row[1:]):
                    if ': ' in v:
                        field,value=v.split(': ',1);d[k][field.lower()]=value
        md=pd.DataFrame.from_dict(d,orient='index');assert tm.GSM.isin(md.index).all()
        joined=md.loc[tm.GSM].reset_index(drop=True)
        er=joined['estrogen receptor status'].str.lower().map({'positive':1.,'negative':0.}).to_numpy()
        meta_info.update(status='EXACT_GSM_JOIN',n_matched=len(joined),er_positive=int((er==1).sum()),er_negative=int((er==0).sum()),er_missing=int(np.isnan(er).sum()),available_fields=list(md.columns))
        # LHC is a sample label; do not silently assert patient identity from uniqueness.
        if 'tumor or normal lhc' in joined:
            meta_info['lhc_unique_tumor_labels']=int(joined['tumor or normal lhc'].nunique())
        paths.append(metadata_path);initial[str(metadata_path)]=sha(metadata_path)
    (agg/'metadata_summary.json').write_text(json.dumps(meta_info,indent=2))
    spec={'scope':'all 174 locked direct relations; no new mappings added',
          'primary':'reuse original p/q; independently verify BH over all evaluable relations',
          'LOO':'descriptive leave-one-specimen-out Spearman for every evaluable primary relation; no LOO hypothesis tests',
          'ER_adjusted':'exploratory partial correlation of ranks conditional on binary author ER status; residualized ranks; 9999 within-ER permutations plus one',
          'ER_family':'BH over all evaluable of 174 relations, separate from original and RNA families',
          'ER_limit':'controls ER alone, not full molecular subtype, purity, age, batch or patient dependence; not independent validation',
          'n_unit':'author mapped specimens; patient independence not verified',
          'permutations':9999,'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}}
    (agg/'analysis_spec_v2.json').write_text(json.dumps(spec,indent=2))
    res=[]; adjusted=[]
    for _,r in primary.iterrows():
        base={'relation_id':r.relation_id,'gene':r.gene,'metabolite':r.metabolite_name,'original_n':r.n,'original_rho':r.rho,'original_q':r.q_value}
        if pd.isna(r.p_value):
            res.append(dict(base,status='NOT_EVALUABLE_ORIGINAL')); adjusted.append(dict(base,status='NOT_EVALUABLE_ORIGINAL'));continue
        assert list(rna.index).count(r.gene)==1 and list(met.index).count(r.metabolite_name)==1
        x=pd.to_numeric(met.loc[r.metabolite_name,tm.MetabID],errors='coerce').to_numpy(float)
        y=pd.to_numeric(rna.loc[r.gene,tm.RNAID],errors='coerce').to_numpy(float)
        valid=np.isfinite(x)&np.isfinite(y);xx,yy=x[valid],y[valid]
        rho=stats.spearmanr(xx,yy).statistic
        assert len(xx)==int(r.n) and abs(rho-r.rho)<1e-10
        loo=np.array([stats.spearmanr(np.delete(xx,i),np.delete(yy,i)).statistic for i in range(len(xx))])
        res.append(dict(base,status='DONE',loo_min=float(np.nanmin(loo)),loo_max=float(np.nanmax(loo)),loo_max_abs_change=float(np.nanmax(abs(loo-rho))),loo_same_sign=int((np.sign(loo)==np.sign(rho)).sum()),loo_n=len(loo)))
        # Describe observed-availability subset without re-running original random P-values.
        avail=np.isfinite(pd.to_numeric(obs.loc[r.metabolite_name,tm.MetabID],errors='coerce').to_numpy(float))
        res[-1]['availability_subset_identical']=bool(np.array_equal(valid,valid&avail))
        v=valid&np.isfinite(er); ev=er[v]; xa,ya=x[v],y[v]
        if len(xa)<8 or min(int((ev==0).sum()),int((ev==1).sum()))<4:
            adjusted.append(dict(base,status='NOT_EVALUABLE_ER_COVERAGE',n=len(xa)));continue
        rx,ry=stats.rankdata(xa),stats.rankdata(ya)
        for group in [0,1]:
            ix=ev==group;rx[ix]-=rx[ix].mean();ry[ix]-=ry[ix].mean()
        norm=np.linalg.norm(rx)*np.linalg.norm(ry)
        if norm==0:
            adjusted.append(dict(base,status='NOT_EVALUABLE_CONSTANT_RESIDUAL'));continue
        observed=float(rx@ry/norm)
        sd=int(hashlib.sha256(('BRCA_ER_v2|'+r.relation_id).encode()).hexdigest()[:8],16)
        rng=np.random.default_rng(sd);exceed=0
        idxs=[np.flatnonzero(ev==g) for g in [0,1]]
        for start in range(0,9999,500):
            count=min(500,9999-start); perms=np.tile(ry,(count,1))
            for idx in idxs:
                order=np.argsort(rng.random((count,len(idx))),axis=1)
                perms[:,idx]=ry[idx][order]
            null=perms@rx/norm;exceed+=int((abs(null)>=abs(observed)-1e-12).sum())
        row=dict(base,status='DONE_EXPLORATORY_ER_ADJUSTED',n=len(xa),n_er_positive=int((ev==1).sum()),n_er_negative=int((ev==0).sum()),partial_rank_rho=observed,p_value=(exceed+1)/10000,seed=sd)
        for group,label in [(0,'negative'),(1,'positive')]:
            ix=ev==group;row['within_er_'+label+'_rho']=float(stats.spearmanr(xa[ix],ya[ix]).statistic)
        adjusted.append(row)
    result=pd.DataFrame(res);adj=pd.DataFrame(adjusted)
    adj['q_value']=np.nan
    if 'p_value' in adj:
        ok=adj.p_value.notna();adj.loc[ok,'q_value']=multipletests(adj.loc[ok,'p_value'],method='fdr_bh')[1]
        adj['family_n_evaluable']=int(ok.sum())
    result.to_csv(agg/'robustness_all_relations.tsv',sep='\t',index=False)
    adj.to_csv(agg/'er_adjusted_all_relations.tsv',sep='\t',index=False)
    focus=result[result.original_q.lt(.05)].merge(adj[['relation_id','partial_rank_rho','p_value','q_value']] if 'partial_rank_rho' in adj else adj[['relation_id','q_value']],on='relation_id',validate='one_to_one')
    focus.to_csv(agg/'original_eleven_relations_followup.tsv',sep='\t',index=False)
    comp=pd.read_csv(paths[5],sep='\t')
    comp['er_adjusted_significant_relations']=comp.gene.map(adj[adj.q_value.lt(.05)].groupby('gene').size()).fillna(0).astype(int)
    comp['followup_status']='Patient robustness added; functional evidence remains previous version until DepMap accessible'
    comp.to_csv(agg/'candidate_comparison_117_v2.tsv',sep='\t',index=False)
    assert all(sha(Path(p))==h for p,h in initial.items())
    pd.DataFrame([{'path':p,'sha256':h} for p,h in initial.items()]).to_csv(agg/'source_manifest.tsv',sep='\t',index=False)
    summary={'relations':len(result),'loo_evaluable':int(result.status.eq('DONE').sum()),'original_significant_relations':len(focus),'ER_significant_relations':int(adj.q_value.lt(.05).sum()),'original_significant_retained_ER':int(focus.q_value.lt(.05).sum()),'all_original_sources_unchanged':True}
    (agg/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    (out/'PATIENT_ANALYSIS_DONE').write_text('Completed; other modules tracked separately')

if __name__=='__main__': main()
