"""Exclude documented tissue conflict; full174 sensitivities and paired117 RNA.

Server only. Original processed values, annotations and historical statistics stay immutable.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import argparse, csv, hashlib, json, platform
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
from camp_per_cancer_repro import calculate
from brca_all174_composition_v2 import estimate, permutation, bootstrap, design

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
BASE=ROOT/'results/collaborative/BRCA/A'
VERSION='camp_pair_sensitivity_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def seed(s):return int(hashlib.sha256((VERSION+'|'+s).encode()).hexdigest()[:8],16)
def correction(d):
    for family,ix in d.groupby('test_family').groups.items():
        ok=d.index.isin(ix)&d.p_value.notna();ps=d.loc[ok,'p_value'].to_numpy();d.loc[ix,'family_n_evaluable']=len(ps)
        if len(ps):
            qs=multipletests(ps,method='fdr_bh')[1];d.loc[ok,'q_value']=qs
            order=np.argsort(ps);check=np.minimum(1,np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1])
            assert np.max(abs(check-qs[order]))<1e-12
    return d
def base(out,gene,family,rid='NA',met='NA'):
    d=dict.fromkeys(PREFIX,np.nan)
    d.update(cancer='BRCA',cohort='CAMP_BRCA1_Terunuma',stage_id='04_ROBUSTNESS',run_id=out.name,analysis_version=VERSION,
             analysis_type=family,metabolite_key=rid.split('|')[0],metabolite_name=met,gene=gene,
             unit='author_case_row',test_family=family,status='NOT_EVALUABLE',reason='NA',source_id='CAMP_GSE37751_original_author_case_table',relation_id=rid)
    return d

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out;pub=out/'public'
    assert (out/'.running').read_text()==VERSION and not (pub/'association_174_sensitivities.tsv').exists()
    spec=dict(version=VERSION,planned_relations=174,planned_genes=117,excluded_accession='GSM927051',exclusion_reason='GEO/CAMP tumor vs original Excel/PDF normal; exclude, do not relabel',
      families=['processed_Spearman','available_Spearman','processed_ER','available_ER','available_ER_PC12','available_ER_MonoEndoFib','paired_RNA'],
      multiplicity='BH separately over evaluable items of each full174 family; paired RNA separately over evaluable of117; missing stay NA',
      permutation_n=9999,spearman_bootstrap=4000,adjusted_bootstrap=500,paired_mean_bootstrap=4000,
      paired_RNA='two-sided paired t test of tumor-normal mean difference on unchanged author microarray scale; t95CI and paired-case bootstrap interval; no fold-change assumption',
      adjusted='same partial rank and within-ER residual permutation method as prior composition analysis; 500 stratified case bootstraps; scores fixed within bootstrap',
      PCA='same eligible marker panels and excluded117 candidates; refit centering/scaling/PCA on60 retained tumors, no parameter search',
      minimum_n_spearman=8,minimum_n_ER=8,minimum_per_ER_group=4,minimum_n_composition=20,
      original_statistics_changed=False,postselection_same_cohort=True,source_patient_identity='author case table, not genetic fingerprint verification',
      software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
    (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
    identity=BASE/'20260921T102429Z_camp_sample_identity_v1/private/audited_mapping.tsv'
    m=pd.read_csv(identity,sep='\t',dtype=str);assert len(m)==108
    conflict=m.TN.ne(m.pdf_TN);assert conflict.sum()==1 and m.loc[conflict,'GSM'].iloc[0]=='GSM927051'
    tm=m[m.TN.eq('Tumor')&~conflict].copy();nm=m[m.TN.eq('Normal')&~conflict].copy()
    assert len(tm)==60 and tm.case_row.is_unique and nm.case_row.is_unique
    pairs=tm[['case_row','RNAID']].merge(nm[['case_row','RNAID']],on='case_row',suffixes=('_tumor','_normal'),validate='one_to_one')
    assert len(pairs)==45;save(pairs,out/'private/paired_RNA_map.tsv')
    src=ROOT/'data/candidates/camp_primary_tissue_multicancer'
    rf=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
    mf=src/'processed_metabolomics'/m.MetabFile.iloc[0]
    rna=pd.read_csv(rf,index_col=0);met=pd.read_excel(mf,sheet_name=tm.MetabFile_sheet.iloc[0],index_col=0);obs=pd.read_excel(mf,sheet_name='data',index_col=0)
    for d in [rna,met,obs]:d.columns=d.columns.astype(str);assert d.index.is_unique and d.columns.is_unique
    priorfile=ROOT/'results/BRCA_117_screen_20260911/tables/BRCA_174_patient_associations.tsv'
    old=pd.read_csv(priorfile,sep='\t');primary=old[old.analysis_type.eq('author_processed_primary')];assert len(primary)==174
    exprfile=ROOT/'results/BRCA_117_screen_20260911/tables/BRCA_117_RNA_differences.tsv';oldexpr=pd.read_csv(exprfile,sep='\t');assert len(oldexpr)==117
    geofile=BASE/'20260919_followup_ER_v2/source/geo_metadata_lines.txt'
    lines=list(csv.reader(geofile.read_text().splitlines(),delimiter='\t'));ids=next(z[1:] for z in lines if z[0]=='!Sample_geo_accession');md=pd.DataFrame(index=ids)
    for z in lines:
        if z[0]=='!Sample_characteristics_ch1':md[z[1].split(': ',1)[0]]=[v.split(': ',1)[1] for v in z[1:]]
    er=md.loc[tm.GSM,'estrogen receptor status'].map({'Positive':1.,'Negative':0.}).to_numpy();assert np.isfinite(er).all()
    compdir=BASE/'20260919T151200Z_all117_robustness_v2';scorefile=compdir/'private/MCPcounter_modified_scores.tsv';coverfile=compdir/'public/composition_marker_coverage.tsv'
    scores=pd.read_csv(scorefile,sep='\t',index_col=0).loc[tm.RNAID];cover=pd.read_csv(coverfile,sep='\t')
    eligible=cover.loc[cover.eligible_PCA,'population'].tolist();zz=scores[eligible].values;zz=(zz-zz.mean(0))/zz.std(0,ddof=1);u,s,v=np.linalg.svd(zz,full_matrices=False);pcs=u[:,:2]*s[:2]
    save(pd.DataFrame(dict(population=eligible,PC1=v[0],PC2=v[1])),pub/'PCA60_loadings.tsv')
    covs={'processed_ER':np.empty((60,0)),'available_ER':np.empty((60,0)),'available_ER_PC12':pcs,'available_ER_MonoEndoFib':scores[['Monocytic lineage','Endothelial cells','Fibroblasts']].values}
    prior_er_file=BASE/'20260919_followup_ER_v2/aggregate/er_adjusted_all_relations.tsv'
    prior_av_file=BASE/'20260919_ER_availability_v3/aggregate/er_availability_all_174.tsv'
    prior_comp_file=compdir/'public/all174_composition_sensitivity.tsv'
    prior_er=pd.read_csv(prior_er_file,sep='\t').set_index('relation_id');prior_av=pd.read_csv(prior_av_file,sep='\t').set_index('relation_id');prior_comp=pd.read_csv(prior_comp_file,sep='\t').set_index(['relation_id','test_family'])
    paths=[identity,rf,mf,priorfile,exprfile,geofile,scorefile,coverfile,prior_er_file,prior_av_file,prior_comp_file,Path(__file__),Path(__file__).parent/'camp_per_cancer_repro.py',Path(__file__).parent/'brca_all174_composition_v2.py']
    hashes={str(p):sha(p) for p in paths};result=[];errors=[]
    for idx,(_,r) in enumerate(primary.iterrows()):
        rid=r.relation_id;features=r.gene in rna.index and r.metabolite_name in met.index and r.metabolite_name in obs.index
        if features:
            x=met.loc[r.metabolite_name,tm.MetabID].to_numpy(float);y=rna.loc[r.gene,tm.RNAID].to_numpy(float);available=np.isfinite(obs.loc[r.metabolite_name,tm.MetabID].to_numpy(float));finite=np.isfinite(x)&np.isfinite(y)
        for family in spec['families'][:-1]:
            b=base(out,r.gene,family,rid,r.metabolite_name);b['effect_type']='Spearman_rho' if 'Spearman' in family else 'partial_rank_correlation'
            if 'Spearman' in family:
                typ='author_processed_primary' if family.startswith('processed') else 'author_data_available_sensitivity';pr=old[(old.relation_id.eq(rid))&old.analysis_type.eq(typ)].iloc[0];ef='rho'
            elif family=='processed_ER':pr=prior_er.loc[rid];ef='partial_rank_rho'
            elif family=='available_ER':pr=prior_av.loc[rid];ef='partial_rank_rho'
            else:pr=prior_comp.loc[(rid,family.replace('available_',''))];ef='effect'
            b.update(previous_n=pr.get('n',np.nan),previous_effect=pr.get(ef,np.nan),previous_p=pr.get('p_value',np.nan),previous_q=pr.get('q_value',np.nan),original_CAMP_g=r.camp_g,original_CAMP_q=r.camp_q)
            if not features or pd.isna(r.p_value):b['reason']='original_unevaluable_or_feature_missing';result.append(b);continue
            mask=finite & (available if family.startswith('available') else True);xx=x[mask];yy=y[mask];ee=er[mask];b['n']=len(xx);sd=seed(family+'|'+rid);b['seed']=sd
            if 'Spearman' in family:
                z=calculate(xx,yy,sd)
                if pd.notna(z['p_value']):b.update(effect=z['rho'],ci_lower=z['ci_lower'],ci_upper=z['ci_upper'],p_value=z['p_value'],status='DONE',bootstrap_valid=z['bootstrap_valid'])
                else:b['reason']=z['status']
            else:
                c=covs[family][mask];minimum=20 if c.shape[1] else 8
                if len(xx)<minimum or min(sum(ee==0),sum(ee==1))<4:b['reason']='insufficient_coverage';result.append(b);continue
                try:rho,q,ex,ey,cond=estimate(xx,yy,ee,c)
                except ValueError as e:b['reason']=str(e);result.append(b);continue
                pp=permutation(rho,q,ex,ey,ee,np.random.default_rng(sd));ci,bn=bootstrap(xx,yy,ee,c,np.random.default_rng(sd+1))
                b.update(effect=rho,p_value=pp,ci_lower=ci[0],ci_upper=ci[1],bootstrap_valid=bn,status='DONE',design_condition=cond)
                z=design(ee,c);rx=stats.rankdata(xx);ry=stats.rankdata(yy);ax=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];ay=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0];errors.append(abs(rho-np.corrcoef(ax,ay)[0,1]))
            b['p_origin']='NEW_EXCLUDED_CASE_ANALYSIS'
            if b['status']=='DONE' and b['n']==b['previous_n'] and family!='available_ER_PC12' and pd.notna(b['previous_p']):
                assert abs(b['effect']-b['previous_effect'])<1e-10
                b['p_value']=b['previous_p'];b['p_origin']='REUSED_IDENTICAL_RETAINED_OBSERVATIONS'
                for c0 in ['ci_lower','ci_upper']:
                    if pd.notna(pr.get(c0,np.nan)):b[c0]=pr[c0]
            result.append(b)
        if idx%20==0:print('relations',idx+1,'/174',flush=True)
    d=correction(pd.DataFrame(result));d['effect_change']=d.effect-d.previous_effect;d['q_threshold_changed']=d.q_value.lt(.05).ne(d.previous_q.lt(.05)) & d.q_value.notna() & d.previous_q.notna();d=d[PREFIX+[c for c in d if c not in PREFIX]]
    assert len(d)==1044 and d.groupby('test_family').relation_id.nunique().eq(174).all();save(d,pub/'association_174_sensitivities.tsv')
    paired=[]
    for _,r in oldexpr.iterrows():
        b=base(out,r.gene,'paired_RNA');b.update(stage_id='03_PATIENT',effect_type='mean_tumor_minus_normal_author_expression_scale',previous_unpaired_p=r.p_value,previous_unpaired_q=r.q_value,previous_unpaired_hedges_g=r.rna_hedges_g)
        if r.gene not in rna.index:b['reason']='RNA_gene_missing';paired.append(b);continue
        x=rna.loc[r.gene,pairs.RNAID_tumor].to_numpy(float);y=rna.loc[r.gene,pairs.RNAID_normal].to_numpy(float);mask=np.isfinite(x)&np.isfinite(y);delta=x[mask]-y[mask];n=len(delta);b.update(n=n,n_reference=n)
        if n<8 or delta.std(ddof=1)==0:b['reason']='insufficient_or_constant_pair_differences';paired.append(b);continue
        test=stats.ttest_rel(x[mask],y[mask]);mean=float(delta.mean());se=delta.std(ddof=1)/np.sqrt(n);ci=stats.t.interval(.95,n-1,loc=mean,scale=se);sd=seed('paired_RNA|'+r.gene);rng=np.random.default_rng(sd);bs=delta[rng.integers(n,size=(4000,n))].mean(1)
        assert abs(test.statistic-mean/se)<1e-10
        b.update(effect=mean,ci_lower=ci[0],ci_upper=ci[1],p_value=test.pvalue,status='DONE',t_statistic=test.statistic,df=n-1,median_paired_difference=float(np.median(delta)),bootstrap_mean_lower=float(np.quantile(bs,.025)),bootstrap_mean_upper=float(np.quantile(bs,.975)),positive_pair_fraction=float((delta>0).mean()),seed=sd)
        paired.append(b)
    e=correction(pd.DataFrame(paired));e=e[PREFIX+[c for c in e if c not in PREFIX]];assert len(e)==117;save(e,pub/'paired_RNA117.tsv')
    assert max(errors)<1e-10 and all(sha(p)==h for p,h in hashes.items())
    save(pd.DataFrame([dict(path=p,sha256=h) for p,h in hashes.items()]),pub/'source_manifest.tsv')
    val=dict(tumor_cases=60,paired_cases=45,source_hashes_unchanged=True,independent_partial_OLS_max_error=max(errors),BH_independently_checked=True,paired_t_formula_checked=True,original_statistics_modified=False,PCA_variance=(s*s/(s*s).sum())[:2].tolist(),families={})
    for k,g in pd.concat([d,e],ignore_index=True).groupby('test_family'):
        val['families'][k]=dict(planned=len(g),evaluable=int(g.p_value.notna().sum()),q_lt005=int(g.q_value.lt(.05).sum()),n_min=float(g.loc[g.p_value.notna(),'n'].min()),n_max=float(g.loc[g.p_value.notna(),'n'].max()))
    (pub/'validation.json').write_text(json.dumps(val,indent=2));print(json.dumps(val,indent=2),flush=True)
    (out/'NUMERICAL_DONE').write_text('Complete; lock retained until publication')

if __name__=='__main__':main()
