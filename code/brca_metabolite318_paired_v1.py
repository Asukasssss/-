"""Server-only fixed 318-metabolite paired analysis, no historical changes."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,hashlib,json,platform
from pathlib import Path
import numpy as np,pandas as pd,scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
BASE=ROOT/'results/collaborative/BRCA/A'
V='metabolite318_paired_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);out=ap.parse_args().out;pub=out/'public'
 assert (out/'.running').read_text()==V and not (pub/'paired_metabolite318.tsv').exists()
 spec=dict(version=V,planned_features=318,pairs=45,primary='author_processed_all45',sensitivity='both_sides_finite_in_author_data_sheet;not_verified_detection_mask',
  test='two-sided paired Wilcoxon signed-rank on unchanged author processed scale;drop exact-zero differences;asymptotic normal tie correction;continuity correction;minimum8 pairs',
  effect='mean tumor-minus-own-normal on author processed scale, NOT Hedges g or concentration fold change',interval='4000 whole-pair percentile bootstrap of mean difference;pointwise;not inversion of signed-rank test',
  multiplicity='BH separately over all evaluable of318 primary and318 availability sensitivity;missing NA;all-zero difference p1',
  assumptions='signed-rank location interpretation assumes symmetric paired-difference distribution;counts are descriptive and not per-patient significance',
  original_statistics_modified=False,imputation='reuse author processed values only;no new imputation/transform/normalization',
  python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
 (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 ip=BASE/'20260921T102429Z_camp_sample_identity_v1/private/audited_mapping.tsv';m=pd.read_csv(ip,sep='\t',dtype=str);m=m[m.TN.eq(m.pdf_TN)]
 t=m[m.TN.eq('Tumor')];n=m[m.TN.eq('Normal')];pairs=t[['case_row','MetabID']].merge(n[['case_row','MetabID']],on='case_row',suffixes=('_tumor','_normal'),validate='one_to_one');assert len(pairs)==45 and pairs.case_row.is_unique
 save(pairs,out/'private/metabolite_pairs.tsv')
 mf=ROOT/'data/candidates/camp_primary_tissue_multicancer/processed_metabolomics/PreprocessedData_BRCA1.xlsx';xp=pd.ExcelFile(mf)
 matrix=pd.read_excel(xp,sheet_name='data_imputed',index_col=0);raw=pd.read_excel(xp,sheet_name='data',index_col=0)
 for d in [matrix,raw]:d.columns=d.columns.astype(str);d.index=d.index.astype(str).str.strip();assert d.index.is_unique and d.columns.is_unique
 scale={}
 for sh in ['metabo_imputed_filtered_Tumor','metabo_imputed_filtered_Normal']:
  d=pd.read_excel(xp,sheet_name=sh,index_col=0);d.columns=d.columns.astype(str);d.index=d.index.astype(str).str.strip();assert d.index.is_unique;err=float(np.nanmax(abs(d.values-matrix.loc[d.index,d.columns].values)));assert err==0;scale[sh]=err
 oldfile=ROOT/'results/tables/cohort_effects.tsv';lockedfile=ROOT/'results/tables/cancer_effects.tsv'
 locked=pd.read_csv(lockedfile,sep='\t');locked=locked[locked.cancer.eq('BRCA')];assert len(locked)==318 and locked.feature_name.is_unique
 old=pd.read_csv(oldfile,sep='\t');old=old[old.dataset.eq('BRCA1')&old.feature_name.isin(locked.feature_name)];assert len(old)==318 and old.feature_name.is_unique
 assert np.allclose(old.set_index('feature_name').loc[locked.feature_name,'wilcoxon_fdr'],locked.effect_fdr)
 paths=[ip,mf,oldfile,lockedfile,Path(__file__)];hashes={str(p):sha(p) for p in paths};rows=[];deltas=[]
 for _,r in old.iterrows():
  name=r.feature_name;assert name in matrix.index and name in raw.index
  x=matrix.loc[name,pairs.MetabID_tumor].to_numpy(float);y=matrix.loc[name,pairs.MetabID_normal].to_numpy(float);finite=np.isfinite(x)&np.isfinite(y);available=finite&np.isfinite(raw.loc[name,pairs.MetabID_tumor].to_numpy(float))&np.isfinite(raw.loc[name,pairs.MetabID_normal].to_numpy(float))
  assert finite.sum()==45
  for family,mask in [('paired_processed',finite),('paired_author_available',available)]:
   delta=x[mask]-y[mask];nn=len(delta);b=dict.fromkeys(PREFIX,np.nan);b.update(cancer='BRCA',cohort='CAMP_BRCA1_Terunuma',stage_id='04_ROBUSTNESS',run_id=out.name,analysis_version=V,analysis_type=family,metabolite_key=r.metabolite_key,metabolite_name=name,gene='NA',unit='author_paired_case',n=nn,n_reference=nn,effect_type='mean_paired_difference_author_scale',test_family=family,status='NOT_EVALUABLE',reason='NA',source_id='CAMP_processed_original_author_pair_table',
    pairs_higher=int((delta>0).sum()),pairs_lower=int((delta<0).sum()),pairs_equal=int((delta==0).sum()),higher_fraction=float((delta>0).mean()) if nn else np.nan,lower_fraction=float((delta<0).mean()) if nn else np.nan,
    n_both_author_available=int(available.sum()),n_pairs_with_at_least_one_author_missing=int(finite.sum()-available.sum()),old_hedges_g=r.hedges_g,old_p=r.wilcoxon_p,old_q=r.wilcoxon_fdr,old_n_tumor=r.n_tumor,old_n_normal=r.n_normal,old_design=r.analysis_design)
   if nn>=8:
    nz=delta[delta!=0];ranks=stats.rankdata(abs(nz));wplus=float(ranks[nz>0].sum());wminus=float(ranks[nz<0].sum())
    p=stats.wilcoxon(delta,zero_method='wilcox',correction=True,alternative='two-sided',method='approx').pvalue if len(nz) else 1.
    sd=int(hashlib.sha256((V+'|'+family+'|'+name).encode()).hexdigest()[:8],16);rng=np.random.default_rng(sd);bs=delta[rng.integers(nn,size=(4000,nn))].mean(1)
    b.update(effect=float(delta.mean()),median_paired_difference=float(np.median(delta)),ci_lower=float(np.quantile(bs,.025)),ci_upper=float(np.quantile(bs,.975)),p_value=float(p),status='DONE',w_positive=wplus,w_negative=wminus,paired_rank_biserial=(wplus-wminus)/(wplus+wminus) if len(nz) else 0.,seed=sd)
   else:b['reason']='fewer_than8_complete_author_available_pairs'
   rows.append(b)
   deltas.append(dict(metabolite_name=name,family=family,values=delta.tolist()))
 d=pd.DataFrame(rows)
 for family,ix in d.groupby('test_family').groups.items():
  ok=d.index.isin(ix)&d.p_value.notna();p=d.loc[ok,'p_value'].values;qs=multipletests(p,method='fdr_bh')[1];d.loc[ok,'q_value']=qs;d.loc[ix,'family_n_evaluable']=len(p)
  order=np.argsort(p);check=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);assert max(abs(check-qs[order]))<1e-12
 d['nonzero_pairs']=d.pairs_higher+d.pairs_lower
 d['normal_approximation_warning']=d.nonzero_pairs.lt(10)&d.p_value.notna()&d.nonzero_pairs.gt(0)
 d['same_mean_direction_as_old_g']=np.sign(d.effect)==np.sign(d.old_hedges_g)
 d['pairs_in_old_direction']=np.where(d.old_hedges_g>0,d.pairs_higher,d.pairs_lower)
 d['fraction_in_old_direction']=d.pairs_in_old_direction/d.n
 d=d[PREFIX+[c for c in d if c not in PREFIX]];assert len(d)==636 and (d.pairs_higher+d.pairs_lower+d.pairs_equal).eq(d.n).all()
 save(d,pub/'paired_metabolite318.tsv');(out/'private/differences_for_R.json').write_text(json.dumps(deltas))
 # Use a plain long table for independent R signed-rank verification, confined to server.
 save(pd.DataFrame([dict(metabolite_name=z['metabolite_name'],family=z['family'],difference=v) for z in deltas for v in z['values']]),out/'private/differences_for_R.tsv')
 assert all(sha(p)==h for p,h in hashes.items());save(pd.DataFrame([dict(path=p,sha256=h) for p,h in hashes.items()]),pub/'source_manifest.tsv')
 val=dict(status='DONE',pairs=45,planned_per_family=318,scale_max_errors=scale,input_hashes_unchanged=True,BH_independently_checked=True,counts_sum_checked=True,old_statistics_modified=False,families={})
 for k,g in d.groupby('test_family'):val['families'][k]=dict(evaluable=int(g.p_value.notna().sum()),q_lt005=int(g.q_value.lt(.05).sum()),n_min=int(g.n.min()),n_max=int(g.n.max()),higher_significant=int((g.q_value.lt(.05)&g.effect.gt(0)).sum()),lower_significant=int((g.q_value.lt(.05)&g.effect.lt(0)).sum()))
 (pub/'validation.json').write_text(json.dumps(val,indent=2));print(json.dumps(val,indent=2),flush=True)
if __name__=='__main__':main()
