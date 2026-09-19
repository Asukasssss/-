"""ER-adjusted author-data-availability sensitivity; server only, aggregates only.

The availability mask is not a verified raw detection mask. Processed values are
reused; no new imputation. Identical observations reuse the v2 P-value and seed.
"""
import argparse,csv,hashlib,json,platform
from pathlib import Path
import numpy as np,pandas as pd,scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--previous',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);(a.out/'.running').write_text('BRCA_A_ER_availability_v3');agg=a.out/'aggregate';agg.mkdir()
    spec={'analysis':'partial rank correlation after ER adjustment in author-data-available subset','planned_relations':174,'minimum_n':8,'minimum_per_ER_group':4,'permutations':9999,'multiplicity':'BH over all evaluable of 174 sensitivity relations; separate from v2 family','identical_inputs':'reuse v2 rho/P/seed; recompute q across full new family','values':'same processed metabolite and RNA values; mask from finite author data worksheet; NOT verified raw detection','unit':'author mapped tumor specimens; independent patients NOT VERIFIED','limitations':'ER alone, no purity/batch/full subtype adjustment; not external validation','python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}
    (agg/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
    src=a.root/'data/candidates/camp_primary_tissue_multicancer';mf=src/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    m=pd.read_csv(mf,dtype=str);m=m[m.Dataset.eq('BRCA1')&m.TN.eq('Tumor')].copy()
    assert len(m)==61 and all(m[c].is_unique for c in ['CommonID','RNAID','MetabID'])
    rf=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0];xf=src/'processed_metabolomics'/m.MetabFile.iloc[0]
    rna=pd.read_csv(rf,index_col=0);met=pd.read_excel(xf,sheet_name=m.MetabFile_sheet.iloc[0],index_col=0);obs=pd.read_excel(xf,sheet_name='data',index_col=0)
    for d in [rna,met,obs]:d.columns=d.columns.astype(str);assert d.columns.is_unique
    gf=a.previous/'source/geo_metadata_lines.txt';rows=list(csv.reader(gf.read_text().splitlines(),delimiter='\t'));ids=rows[0][1:];assert len(ids)==len(set(ids));md=pd.DataFrame(index=ids)
    for row in rows[1:]:md[row[1].split(': ',1)[0]]=[v.split(': ',1)[1] for v in row[1:]]
    gsm=m.RNAID.str.replace('.CEL.gz','',regex=False);assert gsm.isin(md.index).all();assert md.loc[gsm,'tissue type'].eq('Tumor').all()
    er=md.loc[gsm,'estrogen receptor status'].map({'Positive':1.,'Negative':0.}).to_numpy();assert np.isfinite(er).all()
    pf=a.previous/'aggregate/er_adjusted_all_relations.tsv';prior=pd.read_csv(pf,sep='\t');assert len(prior)==174 and prior.relation_id.is_unique
    cf=a.previous/'aggregate/candidate_comparison_117_v2.tsv'
    paths=[mf,rf,xf,gf,pf,cf,Path(__file__)];hashes={str(p):sha(p) for p in paths}
    previous_hashes=pd.read_csv(a.previous/'aggregate/source_manifest.tsv',sep='\t').set_index('path').sha256.to_dict()
    assert all(hashes[str(p)]==previous_hashes[str(p)] for p in [mf,rf,xf,gf])
    results=[];maxerr=0.
    for _,r in prior.iterrows():
        base=dict(relation_id=r.relation_id,gene=r.gene,metabolite=r.metabolite,primary_er_n=r.get('n'),primary_er_rho=r.get('partial_rank_rho'),primary_er_p=r.get('p_value'),primary_er_q=r.get('q_value'),original_q=r.original_q)
        if pd.isna(r.p_value):results.append(dict(base,status='NOT_EVALUABLE_PREVIOUS',origin='NOT_COMPUTED'));continue
        assert list(rna.index).count(r.gene)==1 and list(met.index).count(r.metabolite)==1 and list(obs.index).count(r.metabolite)==1
        x=pd.to_numeric(met.loc[r.metabolite,m.MetabID],errors='coerce').to_numpy(float);y=pd.to_numeric(rna.loc[r.gene,m.RNAID],errors='coerce').to_numpy(float)
        v=np.isfinite(x)&np.isfinite(y)&np.isfinite(er);mask=v&np.isfinite(pd.to_numeric(obs.loc[r.metabolite,m.MetabID],errors='coerce').to_numpy(float));ev=er[mask];xx=x[mask];yy=y[mask]
        base.update(n=int(mask.sum()),n_er_positive=int((ev==1).sum()),n_er_negative=int((ev==0).sum()),n_excluded=int(v.sum()-mask.sum()))
        if len(xx)<8 or min(sum(ev==0),sum(ev==1))<4:results.append(dict(base,status='NOT_EVALUABLE_COVERAGE',origin='NOT_COMPUTED'));continue
        rx,ry=stats.rankdata(xx),stats.rankdata(yy);design=np.column_stack([np.ones(len(ev)),ev]);ex=rx-design@np.linalg.lstsq(design,rx,rcond=None)[0];ey=ry-design@np.linalg.lstsq(design,ry,rcond=None)[0]
        for group in [0,1]:ix=ev==group;rx[ix]-=rx[ix].mean();ry[ix]-=ry[ix].mean()
        norm=np.linalg.norm(rx)*np.linalg.norm(ry)
        if norm==0:results.append(dict(base,status='NOT_EVALUABLE_CONSTANT',origin='NOT_COMPUTED'));continue
        rho=float(rx@ry/norm);maxerr=max(maxerr,abs(rho-np.corrcoef(ex,ey)[0,1]))
        if np.array_equal(mask,v):
            assert len(xx)==r.n and abs(rho-r.partial_rank_rho)<1e-12
            p=r.p_value;seed=int(r.seed);origin='REUSED_IDENTICAL_INPUT_V2_P'
        else:
            seed=int(hashlib.sha256(('BRCA_ER_AVAILABLE_v3|'+r.relation_id).encode()).hexdigest()[:8],16);rng=np.random.default_rng(seed);exceed=0;idxs=[np.flatnonzero(ev==g) for g in [0,1]]
            for start in range(0,9999,500):
                count=min(500,9999-start);perms=np.tile(ry,(count,1))
                for idx in idxs:perms[:,idx]=ry[idx][np.argsort(rng.random((count,len(idx))),axis=1)]
                exceed+=int((abs(perms@rx/norm)>=abs(rho)-1e-12).sum())
            p=(exceed+1)/10000;origin='NEW_SUBSET_CALCULATION'
        results.append(dict(base,status='DONE_EXPLORATORY',origin=origin,partial_rank_rho=rho,p_value=p,seed=seed))
    d=pd.DataFrame(results);ok=d.p_value.notna();d['q_value']=np.nan;d.loc[ok,'q_value']=multipletests(d.loc[ok,'p_value'],method='fdr_bh')[1];d['family_n_evaluable']=int(ok.sum())
    # Independent sorted cumulative minimum implementation checks the new BH family.
    ps=d.loc[ok,'p_value'].to_numpy();order=np.argsort(ps);q=np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1];q=np.minimum(q,1);qerr=float(max(abs(q-d.loc[ok,'q_value'].to_numpy()[order])));assert qerr<1e-12 and maxerr<1e-12
    d.to_csv(agg/'er_availability_all_174.tsv',sep='\t',index=False);d[d.original_q.lt(.05)].to_csv(agg/'original_eleven_sensitivity.tsv',sep='\t',index=False)
    comp=pd.read_csv(cf,sep='\t');assert len(comp)==117 and comp.gene.is_unique
    comp['er_available_n_evaluable']=comp.gene.map(d[ok].groupby('gene').size()).fillna(0).astype(int)
    comp['er_available_n_q_lt_005']=comp.gene.map(d[d.q_value.lt(.05)].groupby('gene').size()).fillna(0).astype(int)
    best=d[ok].sort_values(['q_value','relation_id']).drop_duplicates('gene').set_index('gene')
    for c in ['relation_id','n','partial_rank_rho','q_value']:comp['er_available_best_'+c]=comp.gene.map(best[c])
    comp['patient_followup_category']=np.where(comp.er_available_n_evaluable.eq(0),'NOT_EVALUABLE',np.where(comp.er_available_n_q_lt_005.gt(0)&comp.er_adjusted_significant_relations.gt(0),'SIGNALS_IN_BOTH_ER_FAMILIES','RETAIN_FULL_POOL_NOT_FUNCTIONALLY_EXCLUDED'))
    comp.to_csv(agg/'candidate_comparison_117_v3.tsv',sep='\t',index=False)
    assert all(sha(Path(p))==h for p,h in hashes.items())
    pd.DataFrame([dict(path=p,sha256=h) for p,h in hashes.items()]).to_csv(agg/'source_manifest.tsv',sep='\t',index=False)
    summary=dict(planned=len(d),evaluable=int(ok.sum()),reused_identical=int(d.origin.eq('REUSED_IDENTICAL_INPUT_V2_P').sum()),new_subset_calculations=int(d.origin.eq('NEW_SUBSET_CALCULATION').sum()),q_lt_005=int(d.q_value.lt(.05).sum()),significant_genes=d.loc[d.q_value.lt(.05),'gene'].unique().tolist(),independent_OLS_max_error=maxerr,independent_BH_max_error=qerr,original_sources_unchanged=True)
    (agg/'summary.json').write_text(json.dumps(summary,indent=2));(a.out/'PATIENT_DONE').write_text('Done; functional evidence status separate');(a.out/'.running').unlink();print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
