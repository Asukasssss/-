"""Server-only source-label PYCR1 program associations; public aggregate output only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import argparse,json,hashlib
from pathlib import Path
import h5py,numpy as np,pandas as pd
from scipy import sparse,stats
from brca_sc117_profile_v1 import dataframe,sha256

def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def bh(p):
    p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def rcorr(x,y):
    a=stats.rankdata(x,axis=-1);b=stats.rankdata(y,axis=-1);a-=a.mean(axis=-1,keepdims=True);b-=b.mean(axis=-1,keepdims=True)
    den=np.sqrt((a*a).sum(-1)*(b*b).sum(-1))
    return np.divide((a*b).sum(-1),den,out=np.full(np.shape(den),np.nan),where=den>0)
def infer(x,y,seed):
    rng=np.random.default_rng(seed);n=len(x);rho=float(rcorr(x,y));a=stats.rankdata(x);b=stats.rankdata(y);a-=a.mean();b-=b.mean()
    perm=np.array([rng.permutation(b) for _ in range(9999)]);pr=perm@a/(np.linalg.norm(a)*np.linalg.norm(b));p=(1+np.sum(np.abs(pr)>=abs(rho)-1e-12))/10000
    ix=rng.integers(0,n,(2000,n));boot=rcorr(x[ix],y[ix]);ci=np.nanquantile(boot,[.025,.975]);return rho,p,ci[0],ci[1],int(np.isfinite(boot).sum())

def aggregate(root,previous,cfg,genes):
    cohort=cfg['cohort'];path=previous/'source'/cfg['file']
    cache=root/'private'/(cohort+'_program_counts.tsv')
    if cache.exists():return pd.read_csv(cache,sep='\t')
    with h5py.File(path,'r') as h:
        obs=dataframe(h['obs']);var=dataframe(h[cfg['var_path']]);names=var[cfg['symbol_column']].astype(str).values
        assert obs.index.is_unique
        indices=[np.flatnonzero(names==g) for g in genes];valid=np.array([len(x)==1 for x in indices]);idx=np.array([int(x[0]) if len(x)==1 else -1 for x in indices])
        save(pd.DataFrame({'cohort':cohort,'gene':genes,'matches':[len(x) for x in indices]}),root/'public'/(cohort+'_program_coverage.tsv'))
        typ=obs[cfg['celltype_column']].map(cfg['celltype_map']);selected=typ.isin(['Malignant_epithelial','Fibroblasts']).values
        if cfg.get('filter_column'):selected &= obs[cfg['filter_column']].astype(str).eq(cfg['filter_value']).values
        donor=obs[cfg['donor_column']].astype(str);assert not donor[selected].isin(['NA','nan','unknown','']).any()
        md=pd.DataFrame({'donor':donor,'celltype':typ,'treatment':obs[cfg['treatment_column']].astype(str) if cfg.get('treatment_column') else 'not_available','subtype':obs[cfg['stratum_column']].astype(str) if cfg.get('stratum_column') else 'not_available'})
        codes,uni=pd.factorize(pd.MultiIndex.from_frame(md));meta=uni.to_frame(index=False);meta.columns=md.columns
        mat=h[cfg['raw_path']];shape=tuple(mat.attrs['shape']);assert mat.attrs['encoding-type']=='csr_matrix';ip=mat['indptr'][:]
        counts=np.zeros((len(meta),len(genes)));detect=np.zeros_like(counts);n=np.zeros(len(meta),int);libs=np.zeros(len(meta))
        for start in range(0,len(obs),2000):
            end=min(start+2000,len(obs));keep=selected[start:end]
            if not keep.any():continue
            lo,hi=ip[start],ip[end];v=mat['data'][lo:hi];assert np.isfinite(v).all() and (v>=0).all() and np.allclose(v,np.rint(v),atol=1e-6)
            x=sparse.csr_matrix((v,mat['indices'][lo:hi],ip[start:end+1]-lo),shape=(end-start,shape[1]))[keep]
            lib=np.asarray(x.sum(1)).ravel();assert (lib>0).all();y=np.zeros((len(lib),len(genes)));y[:,valid]=x[:,idx[valid]].toarray();c=codes[start:end][keep]
            np.add.at(counts,c,y);np.add.at(detect,c,y>0);np.add.at(n,c,1);np.add.at(libs,c,lib)
        rows=[]
        for j,r in meta.iterrows():
            if n[j]:
                for k,g in enumerate(genes):rows.append(dict(r,gene=g,raw_count=counts[j,k] if valid[k] else np.nan,n_detected=detect[j,k] if valid[k] else np.nan,n_cells=n[j],library_sum=libs[j]))
        d=pd.DataFrame(rows);save(d,cache)
        print(cohort,'aggregation complete',int(n.sum()),flush=True)
        (root/'public'/(cohort+'_input.json')).write_text(json.dumps({'source':str(path),'sha256':sha256(path),'cells_selected':int(n.sum()),'source_labels':int(d.donor.nunique()),'gene_count':len(genes),'config':cfg},indent=2)+'\n')
        return d

def main(root,previous):
    spec=json.loads((root/'public/analysis_spec.json').read_text());programs=spec['programs'];genes=sorted({'PYCR1'}|{g for v in programs.values() for g in v});results=[];contrasts=[];covers=[];descriptive=[]
    for cf in ['brca_sc117_Wu2021_config.json','brca_sc117_Pal2021_config.json']:
        cfg=json.loads((root/'scripts'/cf).read_text());cohort=cfg['cohort'];d=aggregate(root,previous,cfg,genes)
        parts=[('ALL',d)]+([('UNTREATED',d[d.treatment.eq('Naive')])] if cohort=='Wu2021' else [])
        for partition,part in parts:
            agg=part.groupby(['donor','celltype','gene'],as_index=False)[['raw_count','n_detected','n_cells','library_sum']].sum(min_count=1)
            agg=agg[agg.n_cells>=20];agg['value']=np.log1p(1e6*agg.raw_count/agg.library_sum)
            pertype={}
            for typ in ['Malignant_epithelial','Fibroblasts']:
                q=agg[agg.celltype.eq(typ)];wide=q.pivot(index='donor',columns='gene',values='value');scores=pd.DataFrame(index=wide.index)
                scores['PYCR1']=wide.get('PYCR1',pd.Series(np.nan,index=wide.index))
                for prog,gs in programs.items():
                    measured=[g for g in gs if g in wide and wide[g].notna().all()];variable=[g for g in measured if wide[g].std(ddof=1)>0]
                    coverage=len(measured)/len(gs);covers.append(dict(cohort=cohort,partition=partition,celltype=typ,program=prog,n_sources=len(wide),n_expected=len(gs),n_measured=len(measured),n_variable=len(variable),coverage=coverage))
                    if coverage>=.8 and len(variable)>=3:
                        a=wide[variable];scores[prog]=((a-a.mean())/a.std(ddof=1)).mean(1)
                    else:scores[prog]=np.nan
                    z=scores[['PYCR1',prog]].dropna();ok=len(z)>=10 and z.nunique().min()>1
                    row=dict(cancer='BRCA',cohort=cohort,stage_id='06_EXTERNAL',run_id=root.name,analysis_version='pycr1_context_v1',analysis_type='source_level_program_spearman',metabolite_key='NA',metabolite_name='NA',gene='PYCR1',unit='source_donor_label',n=len(z),n_reference=np.nan,effect_type='Spearman_rho',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family='PRIMARY8' if partition=='ALL' else 'UNTREATED4',family_n_evaluable=np.nan,status='DONE' if ok else 'NOT_EVALUABLE',reason='' if ok else 'coverage_or_less_than10_sources_or_constant',source_id=cfg['source_id'],celltype=typ,program=prog,partition=partition)
                    if ok:
                        rho,p,lo,hi,nb=infer(z.PYCR1.values,z[prog].values,20260920+len(results));row.update(effect=rho,p_value=p,ci_lower=lo,ci_upper=hi,n_bootstrap_valid=nb)
                        assert np.isclose(rho,stats.spearmanr(z.PYCR1,z[prog]).statistic)
                    results.append(row)
                meta=q[q.gene.eq('PYCR1')].set_index('donor');scores=scores.join(meta[['n_cells','n_detected','subtype']] if 'subtype' in meta else meta[['n_cells','n_detected']])
                save(scores.reset_index(),root/'private'/f'{cohort}_{partition}_{typ}_scores.tsv');pertype[typ]=scores
                for (sub,treat),ss in part[part.celltype.eq(typ)].groupby(['subtype','treatment']):
                    donors=set(ss[ss.gene.eq('PYCR1')].donor)&set(scores.index)
                    descriptive.append(dict(cohort=cohort,partition=partition,celltype=typ,subtype=sub,treatment=treat,n_eligible_source_labels=len(donors)))
            for prog in programs:
                a=pertype['Malignant_epithelial'][['PYCR1',prog]];b=pertype['Fibroblasts'][['PYCR1',prog]];ab=a.join(b,lsuffix='_M',rsuffix='_F',how='inner').dropna();n=len(ab)
                row=dict(cohort=cohort,partition=partition,program=prog,n_common_sources=n,effect_type='rho_malignant_minus_fibroblast_same_sources',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,status='NOT_EVALUABLE',reason='less_than10_common_sources_or_constant',p_value=np.nan,q_value=np.nan)
                if n>=10 and ab.nunique().min()>1:
                    ar=ab.to_numpy();delta=float(rcorr(ar[:,0],ar[:,1])-rcorr(ar[:,2],ar[:,3]));rng=np.random.default_rng(20260920);ix=rng.integers(0,n,(2000,n));bs=rcorr(ar[ix,0],ar[ix,1])-rcorr(ar[ix,2],ar[ix,3]);lo,hi=np.nanquantile(bs,[.025,.975]);row.update(effect=delta,ci_lower=lo,ci_upper=hi,status='DONE',reason='descriptive_paired_source_bootstrap_pointwise_no_formal_difference_test')
                contrasts.append(row)
    r=pd.DataFrame(results)
    for family,ix in r.groupby('test_family').groups.items():
        valid=r.loc[ix,'p_value'].notna();r.loc[ix,'q_value']=bh(r.loc[ix,'p_value'].fillna(1));r.loc[np.array(list(ix))[~valid],'q_value']=np.nan;r.loc[ix,'family_n_evaluable']=int(valid.sum())
    save(r,root/'public/pycr1_program_associations.tsv');save(pd.DataFrame(contrasts),root/'public/pycr1_celltype_contrasts.tsv');save(pd.DataFrame(covers),root/'public/program_coverage.tsv');save(pd.DataFrame(descriptive),root/'public/subtype_treatment_coverage.tsv')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,7));colors={'Malignant_epithelial':'#276b9c','Fibroblasts':'#bd7132'}
    for j,row in r.iterrows():
        if row.status=='DONE':ax.plot([row.ci_lower,row.ci_upper],[j,j],c=colors[row.celltype]);ax.plot(row.effect,j,'o',c=colors[row.celltype])
    ax.set_yticks(range(len(r)));ax.set_yticklabels([f'{x.cohort} {x.partition} | {x.celltype} | '+('E2F' if 'E2F' in x.program else 'Collagen') for x in r.itertuples()],fontsize=8);ax.axvline(0,c='gray',ls='--');ax.set_xlim(-1,1);ax.invert_yaxis();ax.set_xlabel('Source-level Spearman rho; pointwise bootstrap interval');ax.set_title('PYCR1 transcriptional programs; not metabolic flux');fig.tight_layout();fig.savefig(root/'public/pycr1_program_forest.png',dpi=170)
    (root/'public/sc_validation.json').write_text(json.dumps({'status':'DONE','rows':len(r),'primary_planned':8,'untreated_planned':4,'no_patient_identity_inference':True,'no_cell_level_pseudoreplication':True,'programs_exclude_PYCR1':all('PYCR1' not in x for x in programs.values()),'script_sha256':sha256(Path(__file__))},indent=2)+'\n')
    print(r[['cohort','partition','celltype','program','n','effect','q_value','status']].to_string(index=False),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);a=p.parse_args();main(a.root,a.previous)
