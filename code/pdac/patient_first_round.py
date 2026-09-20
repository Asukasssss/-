"""PDAC locked-family first round. Run only in a new server165 run directory."""
import argparse,csv,hashlib,json,platform,traceback,tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
METHOD={'method':'Spearman','permutations':9999,'bootstrap':4000,'minimum_n':8,'permutation_comparison_tolerance':1e-14}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,o):p.write_text(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def write(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def fixed_bh(pvalues):
    values=np.asarray([1.0 if p is None else p for p in pvalues],float)
    order=np.argsort(values);q=np.minimum.accumulate((values[order]*len(values)/np.arange(1,len(values)+1))[::-1])[::-1]
    ans=np.empty(len(values));ans[order]=np.minimum(1,q)
    return [None if p is None else float(v) for p,v in zip(pvalues,ans)]
def compute(x,y,seed):
    if len(x)<8:return {'status':'NOT_EVALUABLE','reason':'N_LT_8'}
    if np.ptp(x)==0 or np.ptp(y)==0:return {'status':'NOT_EVALUABLE','reason':'CONSTANT_VALUES'}
    rho=float(stats.spearmanr(x,y).statistic);rx=stats.rankdata(x);ry=stats.rankdata(y);rx-=rx.mean();ry-=ry.mean()
    rng=np.random.default_rng(seed);inds=rng.random((9999,len(x))).argsort(axis=1)
    null=(ry[inds]*rx).sum(axis=1)/np.sqrt((rx@rx)*(ry@ry))
    p=float((1+(np.abs(null)>=abs(rho)-1e-14).sum())/10000)
    bs=rng.integers(0,len(x),size=(4000,len(x)));bx=stats.rankdata(x[bs],axis=1);by=stats.rankdata(y[bs],axis=1)
    bx-=bx.mean(axis=1,keepdims=True);by-=by.mean(axis=1,keepdims=True);den=np.sqrt((bx*bx).sum(axis=1)*(by*by).sum(axis=1));valid=den>0
    boot=(bx[valid]*by[valid]).sum(axis=1)/den[valid]
    lo,hi=np.quantile(boot,[.025,.975])
    return {'status':'DONE','reason':'Exploratory specimen-level independence assumption','rho':rho,'p':p,'lo':float(lo),'hi':float(hi),'bootstrap_valid':int(valid.sum())}
def leave_one_out(x,y):
    if len(x)<9 or np.ptp(x)==0 or np.ptp(y)==0:return {'loo_n_valid':0,'loo_rho_min':'NA','loo_rho_max':'NA','loo_max_abs_change':'NA','loo_sign_changes':'NA'}
    base=float(stats.spearmanr(x,y).statistic);vals=[]
    for j in range(len(x)):
        a=np.delete(x,j);b=np.delete(y,j)
        if np.ptp(a)>0 and np.ptp(b)>0:vals.append(float(stats.spearmanr(a,b).statistic))
    return {'loo_n_valid':len(vals),'loo_rho_min':min(vals),'loo_rho_max':max(vals),'loo_max_abs_change':max(abs(v-base) for v in vals),'loo_sign_changes':sum(np.sign(v)!=np.sign(base) for v in vals)}
def validate_historical_hashes(expected):
    for p,item in expected.items():
        if not item.get('historical') or sha(p)!=item['historical'] or item['current']!=item['historical']:
            raise ValueError('Historical input hash mismatch; reuse prohibited: '+p)
def run(args):
    out=Path(__file__).resolve().parent
    assert out.parent==ROOT/'results/collaborative/PDAC/B','Server-only isolated directory required'
    lock=out/'.running'
    with lock.open('x') as f:f.write('patient_first_round.py\n')
    try:
        public=out/'public';public.mkdir(exist_ok=False)
        expected=json.loads((out/'input_hash_comparison.json').read_text());validate_historical_hashes(expected)
        direct=read(out/'direct_relations_v1.tsv');spec=json.loads((out/'mapping_spec.json').read_text());aliases=json.loads((out/'gene_aliases.json').read_text())
        assert len(direct)==spec['planned_M'] and len({r['gene'] for r in direct})==spec['planned_G']
        assert len(direct)==len({(r['feature_name'],r['metabolite_key'],r['gene']) for r in direct})
        src=ROOT/'data/candidates/camp_primary_tissue_multicancer'
        mp=src/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
        xp=src/'processed_metabolomics/PreprocessedData_PDAC.xlsx'
        rp=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed/GSE62452.hugene10st.gene_symbol.csv'
        mapping=pd.read_csv(mp,dtype=str);m=mapping[mapping.Dataset=='PDAC'];t=m[m.TN=='Tumor'];n=m[m.TN=='Normal']
        assert all(m[k].notna().all() and m[k].is_unique for k in ['CommonID','RNAID','MetabID'])
        rna=pd.read_csv(rp,index_col=0);met=pd.read_excel(xp,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);raw=pd.read_excel(xp,sheet_name='data',index_col=0)
        assert rna.index.is_unique and rna.columns.is_unique and met.index.is_unique and met.columns.is_unique
        assert t.RNAID.isin(rna.columns).all() and n.RNAID.isin(rna.columns).all() and t.MetabID.isin(met.columns).all()
        old=read(out/'historical_associations.tsv');cache={};allrows=[];loo=[];expression=[];resolution=[]
        for gene in sorted({r['gene'] for r in direct}):
            found=([gene] if gene in rna.index else sorted(set(aliases[gene])&set(rna.index)))
            resolved=found[0] if len(found)==1 else None
            resolution.append({'gene':gene,'rna_label':resolved or 'NA','status':'DONE' if resolved else 'NOT_EVALUABLE','reason':'Exact canonical/unique reviewed alias' if resolved else 'Missing or ambiguous gene label','candidate_aliases':';'.join(found)})
            if resolved:
                yt=pd.to_numeric(rna.loc[resolved,t.RNAID],errors='coerce').to_numpy(float);yn=pd.to_numeric(rna.loc[resolved,n.RNAID],errors='coerce').to_numpy(float)
                yt=yt[np.isfinite(yt)];yn=yn[np.isfinite(yn)]
            else:yt=np.array([]);yn=np.array([])
            expression.append({'gene':gene,'rna_label':resolved or 'NA','tumor_n':len(yt),'normal_n':len(yn),
              'tumor_median':float(np.median(yt)) if len(yt) else 'NA','normal_median':float(np.median(yn)) if len(yn) else 'NA',
              'median_difference':float(np.median(yt)-np.median(yn)) if len(yt) and len(yn) else 'NA',
              'scale':'author_processed_microarray_no_new_transform','p_value':'NA','q_value':'NA','status':'DONE' if len(yt) and len(yn) else 'NOT_EVALUABLE',
              'reason':'Descriptive only; independence and pairing unconfirmed'})
        resolved={r['gene']:r['rna_label'] for r in resolution if r['status']=='DONE'}
        for idx,rel in enumerate(direct):
            f=rel['feature_name'];g=rel['gene'];base_reason=''
            if f not in met.index or f not in raw.index:base_reason='METABOLITE_NOT_UNIQUE_OR_ABSENT'
            if g not in resolved:base_reason='RNA_GENE_NOT_UNIQUE_OR_ABSENT'
            x=pd.to_numeric(met.loc[f,t.MetabID],errors='coerce').to_numpy(float) if f in met.index else np.full(len(t),np.nan)
            y=pd.to_numeric(rna.loc[resolved[g],t.RNAID],errors='coerce').to_numpy(float) if g in resolved else np.full(len(t),np.nan)
            avail=np.isfinite(pd.to_numeric(raw.loc[f,t.MetabID],errors='coerce').to_numpy(float)) if f in raw.index else np.zeros(len(t),bool)
            for kind in ['primary','availability']:
                keep=np.isfinite(x)&np.isfinite(y)
                if kind=='availability':keep &= avail
                xx=x[keep];yy=y[keep]
                key=hashlib.sha256(rel['relation_id'].encode()+json.dumps(METHOD,sort_keys=True).encode()+keep.tobytes()+xx.tobytes()+yy.tobytes()).hexdigest()
                seed=int(key[:8],16);origin='NEW_STATISTICS'
                if key in cache:result=cache[key];origin='REUSED_IDENTICAL_INPUT_CACHE'
                else:
                    previous=[r for r in old if r['gene']==g and r['metabolite_name']==f and r['analysis_type']==('author_processed_primary' if kind=='primary' else 'author_data_available_sensitivity')]
                    if previous and not base_reason and len(xx)>=8 and np.ptp(xx)>0 and np.ptp(yy)>0:
                        p=previous[0];assert len(xx)==int(p['n']) and abs(float(p['rho'])-float(stats.spearmanr(xx,yy).statistic))<1e-12
                        result={'status':'DONE','reason':'Historical rho/P/CI reused after strict source hash check','rho':float(p['rho']),'p':float(p['p_value']),'lo':float(p['ci_lower']),'hi':float(p['ci_upper']),'bootstrap_valid':int(p['bootstrap_valid'])};origin='REUSED_HISTORICAL_STATISTICS'
                    else:result={'status':'NOT_EVALUABLE','reason':base_reason} if base_reason else compute(xx,yy,seed)
                    cache[key]=result
                row={'cancer':'PDAC','cohort':'PDAC','stage_id':'03_PATIENT','run_id':out.name,'analysis_version':'PDAC_first_round_v1',
                  'analysis_type':kind,'metabolite_key':rel['metabolite_key'],'metabolite_name':f,'gene':g,'unit':'author_mapped_tumor_specimens_independence_assumed','n':len(xx),'n_reference':'NA',
                  'effect_type':'Spearman_rho','effect':result.get('rho','NA'),'ci_lower':result.get('lo','NA'),'ci_upper':result.get('hi','NA'),
                  'p_value':result.get('p','NA'),'q_value':'NA','test_family':'PDAC_v1_'+kind,'family_n_evaluable':0,'status':result['status'],'reason':result['reason'],
                  'source_id':'CAMP_PDAC_GSE62452','relation_id':rel['relation_id'],'family_n_planned':len(direct),'n_processed_metabolite':int(np.isfinite(x).sum()),
                  'n_preimputation_available':int(avail.sum()),'metabolite_unique_values':len(np.unique(xx)),'rna_unique_values':len(np.unique(yy)),
                  'input_hash':key,'seed':seed if origin=='NEW_STATISTICS' else 'NA_REUSED','statistic_origin':origin,'original_camp_q':rel['original_q'],
                  'original_camp_effect':rel['original_effect'],'bootstrap_valid':result.get('bootstrap_valid','NA'),'rna_label':resolved.get(g,'NA')}
                allrows.append(row)
                loo.append({'relation_id':rel['relation_id'],'feature_name':f,'gene':g,'analysis_type':kind,**leave_one_out(xx,yy)})
            if idx%30==0:print('relationships',idx+1,'/',len(direct),flush=True)
        for kind in ['primary','availability']:
            selected=[r for r in allrows if r['analysis_type']==kind];ps=[r['p_value'] if r['status']=='DONE' else None for r in selected]
            qs=fixed_bh(ps)
            for r,q in zip(selected,qs):r['q_value']=q if q is not None else 'NA';r['family_n_evaluable']=sum(p is not None for p in ps)
        assert len(allrows)==2*len(direct)
        write(public/'associations.tsv',allrows);write(public/'leave_one_out.tsv',loo);write(public/'gene_expression_descriptive.tsv',expression);write(public/'gene_resolution.tsv',resolution)
        summary={'planned_M':len(direct),'planned_G':len(expression),'rows':len(allrows),'patient_identity':'NOT_SEPARATELY_VERIFIED','historical_hash_gate':'PASS'}
        for kind in ['primary','availability']:
            selected=[r for r in allrows if r['analysis_type']==kind];summary[kind]={'evaluable':sum(r['status']=='DONE' for r in selected),'q_lt_0_05':sum(r['status']=='DONE' and r['q_value']<.05 for r in selected)}
        dump(public/'summary.json',summary);dump(public/'analysis_spec.json',{'method':METHOD,'mapping':spec,'code_base_commit':args.code_commit,'script_sha256':sha(__file__),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__})
        write(public/'source_manifest.tsv',[{'path':str(p),'sha256':sha(p)} for p in [mp,xp,rp,out/'direct_relations_v1.tsv',out/'mapping_spec.json',out/'gene_aliases.json',out/'historical_associations.tsv',Path(__file__)]])
        dump(public/'validation.json',{'status':'PASS','all_planned_rows_retained':True,'fixed_M_BH':True,'historical_hash_gate':True,'not_verified':['patient independence','causal function','direction of frozen g']})
        with tarfile.open(out/'public_delivery.tar','w') as tar:
            for p in sorted(public.iterdir()):tar.add(p,arcname=p.name)
        dump(out/'DONE.json',summary);lock.unlink();print(json.dumps(summary))
    except Exception:(out/'FAILED.txt').write_text(traceback.format_exc());raise
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--code-commit',required=True);run(parser.parse_args())
