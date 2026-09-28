"""Server-only paired PDAC metabolite analysis using explicit GEO author pair IDs."""
import argparse,csv,hashlib,itertools,json,platform,traceback
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PAIR_SOURCE=ROOT/'results/collaborative/PDAC/B/20260921T120900Z_paired_source_v1/GSE62452_paired_sample_information.xlsx'
SOURCE=ROOT/'data/candidates/camp_primary_tissue_multicancer'
MAPPING=SOURCE/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
MATRIX=SOURCE/'processed_metabolomics/PreprocessedData_PDAC.xlsx'
FROZEN=ROOT/'results/collaborative/PDAC/B/20260920T054500Z_source_readiness_v1/cancer_effects.tsv'
EXPECTED={MAPPING:'7a5332cfdeba3cc055d79ef223c95a3c8f18ecf6a7058fafd72b390aec62dce7',MATRIX:'e5f408848824394d6be875d7284b77d0a0aabd433d0f2db8ebe01886c29a7ecd',FROZEN:'1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def write(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def pair_join(mapping,author):
    m=mapping.loc[mapping.Dataset.eq('PDAC')].copy()
    for c in ['CommonID','MetabID','RNAID']:
        if m[c].isna().any() or not m[c].is_unique:raise ValueError('Ambiguous CAMP '+c)
    m['GSE ID']=m.RNAID.str.extract(r'^(GSM\d+)_',expand=False)
    if m['GSE ID'].isna().any() or not m['GSE ID'].is_unique or not author['GSE ID'].is_unique:raise ValueError('Ambiguous accession join')
    for pid,g in author.dropna(subset=['Paired Sample ID']).groupby('Paired Sample ID'):
        if len(g)!=2 or sorted(g['T/N'])!=['N','T']:raise ValueError('Invalid author pair group')
    j=m.merge(author,on='GSE ID',how='left',validate='one_to_one',indicator=True)
    if not j['_merge'].eq('both').all() or not j.TN.map({'Tumor':'T','Normal':'N'}).eq(j['T/N']).all():raise ValueError('Unmatched accession or tissue disagreement')
    pairs=[]
    for pid,g in j.dropna(subset=['Paired Sample ID']).groupby('Paired Sample ID',sort=True):
        if len(g)==2:
            t=g.loc[g.TN.eq('Tumor')].iloc[0];n=g.loc[g.TN.eq('Normal')].iloc[0]
            pairs.append({'author_pair_id':pid,'tumor_metab_id':t.MetabID,'normal_metab_id':n.MetabID})
    if not pairs:raise ValueError('No confirmed complete pairs')
    return j,pd.DataFrame(pairs)
def paired_stat(d,seed):
    d=np.asarray(d,float)
    if len(d)<6:return {'status':'NOT_EVALUABLE','reason':'N_COMPLETE_PAIRS_LT_6'}
    nz=d[d!=0];r=rankdata(abs(nz));observed=float(r[nz>0].sum())
    if len(nz):
        sums=np.array([sum(r*np.asarray(s)) for s in itertools.product([0,1],repeat=len(nz))])
        dist=abs(observed-r.sum()/2);p=float(np.mean(abs(sums-r.sum()/2)>=dist-1e-12));rb=float(2*observed/r.sum()-1)
    else:p=1.;rb=0.
    rng=np.random.default_rng(seed);boot=d[rng.integers(0,len(d),size=(10000,len(d)))].mean(axis=1);lo,hi=np.quantile(boot,[.025,.975])
    return {'status':'DONE','reason':'Author-confirmed pairs;exact signed-rank sign enumeration',
        'mean_delta':float(d.mean()),'median_delta':float(np.median(d)),'ci_lower':float(lo),'ci_upper':float(hi),
        'rank_biserial':rb,'p_value':p,'n_nonzero':len(nz)}
def bh_fixed(values):
    p=np.array([1 if x is None else x for x in values]);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return ['NA' if x is None else float(y) for x,y in zip(values,q)]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--code-commit',required=True);a=ap.parse_args();run=Path(__file__).resolve().parent
    if run.parent!=ROOT/'results/collaborative/PDAC/B':raise ValueError('Server isolated run required')
    with (run/'.running').open('x') as f:f.write('paired_metabolites_v1')
    try:
        out=run/'public';out.mkdir(exist_ok=False)
        for p,h in EXPECTED.items():
            if sha(p)!=h:raise ValueError('Source hash mismatch: '+str(p))
        if sha(PAIR_SOURCE)!='0f16963287abb1bf7bf73dd5a2738f3bdccfb0344b2db3db72598199432ac866':raise ValueError('Author pairing source changed')
        mapping=pd.read_csv(MAPPING,dtype=str);author=pd.read_excel(PAIR_SOURCE,dtype=str)
        joined,pairs=pair_join(mapping,author)
        if len(joined)!=39 or len(pairs)!=11:raise ValueError('Preflight sample counts changed')
        joined.to_csv(run/'private_author_join.tsv',sep='\t',index=False);pairs.to_csv(run/'private_complete_pairs.tsv',sep='\t',index=False)
        t=pd.read_excel(MATRIX,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);n=pd.read_excel(MATRIX,sheet_name='metabo_imputed_filtered_Normal',index_col=0);raw=pd.read_excel(MATRIX,sheet_name='data',index_col=0);im=pd.read_excel(MATRIX,sheet_name='data_imputed',index_col=0)
        for d in [t,n,raw,im]:
            if not d.index.is_unique or not d.columns.is_unique:raise ValueError('Duplicate matrix labels')
        for d in [t,n]:
            if not np.allclose(d.to_numpy(float),im.loc[d.index,d.columns].to_numpy(float),equal_nan=True):raise ValueError('Processed sheets differ from common imputed matrix')
        features=sorted(set(t.index)&set(n.index))
        if len(features)!=307:raise ValueError('Prespecified feature family changed')
        frozen=pd.read_csv(FROZEN,sep='\t');frozen=frozen[frozen.cancer.eq('PDAC')].set_index('feature_name')
        actual_by_trim={f.strip():f for f in features}
        if len(actual_by_trim)!=len(features) or len(set(frozen.index.str.strip()))!=len(frozen):raise ValueError('Ambiguous whitespace normalization')
        frozen['frozen_feature_name']=frozen.index
        frozen.index=[actual_by_trim.get(f.strip(),f) for f in frozen.index]
        write(out/'feature_family.tsv',[{'cancer':'PDAC','feature_name':f,'in_frozen303':f in frozen.index,'frozen_feature_name':frozen.loc[f,'frozen_feature_name'] if f in frozen.index else 'NA','name_bridge':'EXACT_OR_UNIQUE_OUTER_WHITESPACE_ONLY','original_camp_q':float(frozen.loc[f,'effect_fdr']) if f in frozen.index else 'NA'} for f in features])
        rows=[]
        for f in features:
            x=t.loc[f,pairs.tumor_metab_id].to_numpy(float);y=n.loc[f,pairs.normal_metab_id].to_numpy(float)
            finite=np.isfinite(x)&np.isfinite(y);available=finite&np.isfinite(raw.loc[f,pairs.tumor_metab_id].to_numpy(float))&np.isfinite(raw.loc[f,pairs.normal_metab_id].to_numpy(float))
            seed=int(hashlib.sha256(('PDAC_paired_v1|'+f).encode()).hexdigest()[:8],16);cache={}
            for mode,mask in [('primary',finite),('both_preimputation_available',available)]:
                key=mask.tobytes();reused=key in cache
                if not reused:cache[key]=paired_stat((x-y)[mask],seed)
                s=cache[key]
                rows.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','stage_id':'04_ROBUSTNESS','run_id':run.name,'analysis_version':'PDAC_paired_metabolites_v1','analysis_type':mode,'metabolite_key':frozen.loc[f,'metabolite_key'] if f in frozen.index else 'NAME:'+f,'metabolite_name':f,'gene':'NA','unit':'author_confirmed_patient_pair','n':int(mask.sum()),'n_reference':11,'effect_type':'paired_mean_tumor_minus_normal_author_processed_scale','effect':s.get('mean_delta','NA'),'ci_lower':s.get('ci_lower','NA'),'ci_upper':s.get('ci_upper','NA'),'p_value':s.get('p_value','NA'),'q_value':'NA','test_family':'PDAC_paired_v1_'+mode,'family_n_evaluable':'NA','status':s['status'],'reason':s['reason'],'source_id':'CAMP_mapping_plus_GSE62452_author_pair_table','median_delta':s.get('median_delta','NA'),'rank_biserial':s.get('rank_biserial','NA'),'n_nonzero':s.get('n_nonzero','NA'),'n_both_preimputation_available':int(available.sum()),'family_n_planned':307,'original_camp_g':float(frozen.loc[f,'hedges_g']) if f in frozen.index else 'NA','original_camp_q':float(frozen.loc[f,'effect_fdr']) if f in frozen.index else 'NA','same_input_reused':reused,'seed':seed})
        summary={'complete_pairs':11,'unpaired_tumor_excluded':16,'unpaired_normal_excluded':1,'family_n_planned':307,'frozen_overlap':len(set(features)&set(frozen.index)),'frozen_not_in_family':sorted(set(frozen.index)-set(features)),'extra_to_frozen':sorted(set(features)-set(frozen.index))}
        for mode in ['primary','both_preimputation_available']:
            rr=[r for r in rows if r['analysis_type']==mode];q=bh_fixed([None if r['p_value']=='NA' else r['p_value'] for r in rr]);ne=sum(r['status']=='DONE' for r in rr)
            for r,v in zip(rr,q):r['q_value']=v;r['family_n_evaluable']=ne
            summary[mode]={'evaluable':ne,'q_lt_005':sum(v!='NA' and v<.05 for v in q)}
        write(out/'results.tsv',rows);dump(out/'summary.json',summary)
        dump(out/'analysis_spec.json',{'code_commit':a.code_commit,'method':'Two-sided exact paired signed-rank sign enumeration;zero differences dropped;average ranks for ties','assumption':'Symmetric paired-difference distribution for signed-rank location interpretation','family':'All307 shared processed features;fixed307 BH separately per analysis;uncomputable internal p=1 public NA','minimum_pairs':6,'bootstrap':'10000 paired resamples;percentile95% CI for mean difference;same-input reuse','scale':'Author processed scale;no added transform;not asserted log2FC','contrast':'Tumor minus same author-pair Normal','selection':'No use of old42 significance','new_imputation':False,'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}})
        write(out/'source_manifest.tsv',[{'path':str(p),'sha256':sha(p)} for p in [MAPPING,MATRIX,PAIR_SOURCE,FROZEN,Path(__file__)]])
        dump(out/'validation.json',{'status':'PASS','exact_GSM_join':39,'tissue_concordance':True,'explicit_pairs':11,'no_identifier_suffix_pair_inference':True,'shared_processed_matrix_projection_verified':True,'historical_files_modified':False,'private_pair_map_exported':False,'not_verified':['clinical covariates','biochemical identities','original frozen contrast code','all27 tumor patient independence']})
        dump(run/'DONE.json',summary);(run/'.running').unlink();print(json.dumps(summary))
    except Exception:
        (run/'FAILED.txt').write_text(traceback.format_exc());raise
if __name__=='__main__':main()
