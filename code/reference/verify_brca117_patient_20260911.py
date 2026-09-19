"""Independent numerical and coverage checks for the completed BRCA screen."""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from run_brca117_patient_20260911 import expression,hedges,add_q,ROOT,OUT,SRC

def main():
    tab=OUT/'tables'
    assoc=pd.read_csv(tab/'BRCA_174_patient_associations.tsv',sep='\t')
    expr=pd.read_csv(tab/'BRCA_117_RNA_differences.tsv',sep='\t')
    edges=pd.read_csv(OUT/'inputs/BRCA_direct_human_gene_edges.tsv',sep='\t')
    sm=pd.read_csv(OUT/'inputs/BRCA_significant_186_with_mapping.tsv',sep='\t').set_index('metabolite_key')
    mp=pd.read_csv(OUT/'author_sample_mapping.tsv',sep='\t',dtype=str)
    tm=mp[mp.TN.eq('Tumor')];nm=mp[mp.TN.eq('Normal')]
    rp=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/mp.RNAFile.iloc[0]
    xp=SRC/'processed_metabolomics'/mp.MetabFile.iloc[0]
    rna=pd.read_csv(rp,index_col=0);rna.columns=rna.columns.astype(str)
    met=pd.read_excel(xp,sheet_name=tm.MetabFile_sheet.iloc[0],index_col=0);met.columns=met.columns.astype(str)
    obs=pd.read_excel(xp,sheet_name='data',index_col=0);obs.columns=obs.columns.astype(str)
    checks={}
    checks['all_174_primary_and_sensitivity_rows']=len(assoc)==348 and all(len(v)==174 for _,v in assoc.groupby('analysis_type'))
    checks['all_117_gene_rows']=len(expr)==117 and expr.gene.is_unique
    keys=set(zip(edges.metabolite_key,edges.gene_symbol))
    checks['exact_locked_edge_coverage']=all(set(zip(v.metabolite_key,v.gene))==keys for _,v in assoc.groupby('analysis_type'))
    nchecked=0
    for _,r in assoc[assoc.p_value.notna()].iterrows():
        x=pd.to_numeric(met.loc[r.metabolite_name,tm.MetabID],errors='coerce').to_numpy(float)
        y=pd.to_numeric(rna.loc[r.gene,tm.RNAID],errors='coerce').to_numpy(float)
        if r.analysis_type=='author_data_available_sensitivity':
            avail=pd.to_numeric(obs.loc[r.metabolite_name,tm.MetabID],errors='coerce').to_numpy(float)
            x=np.where(np.isfinite(avail),x,np.nan)
        k=np.isfinite(x)&np.isfinite(y)
        ref=np.corrcoef(stats.rankdata(x[k]),stats.rankdata(y[k]))[0,1]
        assert k.sum()==r.n and np.isclose(ref,r.rho,atol=1e-12)
        assert r.camp_g==sm.loc[r.metabolite_key].hedges_g
        assert r.camp_q==sm.loc[r.metabolite_key].effect_fdr
        nchecked+=1
    checks['independent_rank_correlation_rows']=nchecked
    checks['bh_independent_statsmodels']=True
    for _,v in assoc.groupby('test_family'):
        v=v[v.p_value.notna()]
        assert np.allclose(v.q_value,multipletests(v.p_value,method='fdr_bh')[1])
    valid=expr[expr.p_value.notna()]
    assert np.allclose(valid.q_value,multipletests(valid.p_value,method='fdr_bh')[1])
    for _,r in valid.iterrows():
        a=pd.to_numeric(rna.loc[r.gene,tm.RNAID],errors='coerce').to_numpy(float)
        b=pd.to_numeric(rna.loc[r.gene,nm.RNAID],errors='coerce').to_numpy(float)
        a=a[np.isfinite(a)];b=b[np.isfinite(b)]
        ref=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic').pvalue
        assert len(a)==r.n_tumor and len(b)==r.n_normal and np.isclose(ref,r.p_value)
        sd=np.sqrt(((len(a)-1)*np.var(a,ddof=1)+(len(b)-1)*np.var(b,ddof=1))/(len(a)+len(b)-2))
        g=(1-3/(4*(len(a)+len(b))-9))*(np.mean(a)-np.mean(b))/sd
        assert np.isclose(g,r.rna_hedges_g)
    checks['independent_expression_rows']=len(valid)
    x=np.arange(1.,18.);y=np.arange(3.,20.)
    a=expression(x,y,125);b=expression(y,x,125)
    assert a['rna_hedges_g']<0 and np.isclose(a['rna_hedges_g'],-b['rna_hedges_g'])
    assert np.isclose(a['p_value'],b['p_value'])
    assert expression(np.ones(9),np.ones(9),1)['status']=='SKIPPED_CONSTANT_VALUES'
    checks['contrast_sign_symmetry_and_constant_skip']=True
    checks['passed']=all(v is not False for v in checks.values())
    (OUT/'numerical_validation.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
