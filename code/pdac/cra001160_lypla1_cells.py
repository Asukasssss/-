"""Exploratory pooled-cell comparisons; no claim of independent cell replicates."""
from pathlib import Path
import hashlib,json,platform
import numpy as np,pandas as pd,scipy
from scipy import stats
RUN=Path(__file__).resolve().parent
BASE=RUN.parent/'20260927T130900Z_cra001160_lypla1_v1'
INPUT=BASE/'private/cell_values.tsv.gz'
EXPECTED='39c164112b49c9807e34b378b912f7ad510e3faacc63ed415f28f624960c8d8c'
D1='Ductal cell type 1';D2='Ductal cell type 2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=RUN/'public';out.mkdir(exist_ok=True)
 spec={'version':'cra001160_lypla1_cells_v1','gene':'LYPLA1','unit':'pooled_cell',
 'source_run':BASE.name,'input_sha256':EXPECTED,'expression':'Inherited log1p(10000 * gene UMI / all-gene UMI)',
 'contexts':['Tumor_all_units','Tumor_vs_control','Tumor_14_shared_units_sensitivity'],
 'populations':['All_cells_including_zero','Detected_cells_raw_UMI_gt0'],
 'test':'Two-sided asymptotic Mann-Whitney U, tie correction and continuity correction',
 'multiplicity':'BH across all six prespecified exploratory comparisons',
 'effect':'Rank-biserial correlation = 2U/(n_test*n_reference)-1; positive means higher type2 distribution',
 'all_units_rule':'No minimum cell threshold; retain every author-annotated ductal cell',
 'shared_units_rule':'Both ductal groups >=20 cells per T unit, before LYPLA1-positive filtering',
 'no_pseudoreplication_correction':True,'CI':'Not supplied: naive cell CIs would imply independent cells',
 'data_residency':'Per-cell/per-unit data remain on server; export aggregate histograms and statistics only'}
 (out/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 assert sha(INPUT)==EXPECTED,'Input changed'
 c=pd.read_csv(INPUT,sep='\t');assert len(c)==57530 and c.cell.is_unique
 assert np.isfinite(c.log1p_10k).all() and (c.library>0).all()
 assert np.allclose(c.log1p_10k,np.log1p(c['count']*10000/c.library))
 assert ((c['count']>0)==(c.log1p_10k>0)).all()
 tumor=c[c.tissue=='Tumor'];counts=tumor.groupby(['unit','cell_type']).size().unstack(fill_value=0)
 shared=counts.index[(counts[D1]>=20)&(counts[D2]>=20)]
 assert len(shared)==14
 rows=[];profiles=[];hist=[];checks=[]
 for context in spec['contexts']:
  a=c[(c.tissue=='Tumor')&(c.cell_type==D2)].copy()
  b=c[(c.tissue==('Control' if context=='Tumor_vs_control' else 'Tumor'))&(c.cell_type==D1)].copy()
  if context=='Tumor_14_shared_units_sensitivity':a=a[a.unit.isin(shared)];b=b[b.unit.isin(shared)]
  for population in spec['populations']:
   aa=a[a['count']>0] if population.startswith('Detected') else a
   bb=b[b['count']>0] if population.startswith('Detected') else b
   x=aa.log1p_10k.to_numpy();y=bb.log1p_10k.to_numpy()
   test=stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic',use_continuity=True)
   # Independent U and tie-corrected normal tail verification from pooled ranks.
   z=np.r_[x,y];u=stats.rankdata(z)[:len(x)].sum()-len(x)*(len(x)+1)/2
   _,ties=np.unique(z,return_counts=True);n=len(z)
   sd=np.sqrt(len(x)*len(y)/12*((n+1)-sum(ties.astype(float)**3-ties)/(n*(n-1))))
   pcheck=min(1.,2*stats.norm.sf((abs(u-len(x)*len(y)/2)-.5)/sd)) if sd else 1.
   assert np.isclose(u,test.statistic) and np.isclose(pcheck,test.pvalue,rtol=1e-8,atol=1e-300)
   checks.append({'context':context,'population':population,'U_rank_reconstruction':True,'P_tie_corrected_reconstruction':True})
   row=dict(cancer='PDAC',cohort='CRA001160',stage_id='06_EXTERNAL',run_id=RUN.name,analysis_version=spec['version'],
    analysis_type=context+'__'+population,metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit='pooled_cell',
    n=len(x),n_reference=len(y),effect_type='rank_biserial_type2_minus_type1',effect=2*u/(len(x)*len(y))-1,
    ci_lower=np.nan,ci_upper=np.nan,p_value=test.pvalue,q_value=np.nan,test_family='six_exploratory_cell_comparisons',
    family_n_evaluable=6,status='DONE',reason='Naive pooled-cell exploratory P; donor dependence uncorrected',source_id=EXPECTED,
    context=context,population=population,units_test=aa.unit.nunique(),units_reference=bb.unit.nunique(),
    mean_test=x.mean(),mean_reference=y.mean(),median_test=np.median(x),median_reference=np.median(y),
    mean_difference=x.mean()-y.mean(),detected_fraction_test=(a['count']>0).mean(),detected_fraction_reference=(b['count']>0).mean(),U=u)
   rows.append(row)
   edges=np.linspace(0,max(z.max(),.1),161)
   for role,df,v in [('type2',aa,x),('type1',bb,y)]:
    shares=df.unit.value_counts()/len(df)
    qs=np.quantile(v,[.0,.25,.5,.75,1.])
    profiles.append(dict(context=context,population=population,role=role,n_cells=len(v),n_units=df.unit.nunique(),mean=v.mean(),
      minimum=qs[0],q25=qs[1],median=qs[2],q75=qs[3],maximum=qs[4],largest_unit_fraction=shares.iloc[0],top3_units_fraction=shares.iloc[:3].sum(),
      zero_fraction=np.mean(v==0)))
    values,_=np.histogram(v,bins=edges);assert values.sum()==len(v)
    for left,right,number in zip(edges[:-1],edges[1:],values):hist.append(dict(context=context,population=population,role=role,left=left,right=right,n_cells=int(number),fraction=number/len(v)))
 r=pd.DataFrame(rows);order=np.argsort(r.p_value.values);adjusted=np.minimum.accumulate((r.p_value.values[order]*6/np.arange(1,7))[::-1])[::-1].clip(0,1)
 r.loc[order,'q_value']=adjusted
 for frame,name in [(r,'results.tsv'),(pd.DataFrame(profiles),'group_profiles.tsv'),(pd.DataFrame(hist),'histogram_source.tsv')]:frame.to_csv(out/name,sep='\t',index=False,na_rep='NA')
 pd.DataFrame([{'source_id':INPUT.name,'server_path':str(INPUT),'sha256':EXPECTED,'origin':'Original GSA counts, extracted in prior run'},
 {'source_id':Path(__file__).name,'server_path':str(Path(__file__)),'sha256':sha(Path(__file__)),'origin':'code/pdac/cra001160_lypla1_cells.py'}]).to_csv(out/'source_manifest.tsv',sep='\t',index=False)
 (out/'validation.json').write_text(json.dumps({'status':'PASS','input_sha256_verified':True,'all_cells_unique':True,
  'raw_count_normalization_checked':True,'positive_filter_checked':True,'histogram_counts_conserved':True,'independent_statistic_checks':checks,
  'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,
  'not_validated':['No donor dependence adjustment for cell tests','No new malignant identity validation','No LYPLA1 depth-adjusted detection analysis']},indent=2))
 print(r[['context','population','n','n_reference','mean_test','mean_reference','effect','p_value','q_value']].to_string(index=False))
if __name__=='__main__':main()
