"""Frozen pooled-cell LYPLA1 grouping, whole-transcriptome DE and enrichment."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,importlib.util,time
import numpy as np,pandas as pd,h5py,scipy,requests,gseapy
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
BASE=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
PREV=BASE/'20260925T133935Z_lypla1_KEGGscore_v1'
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.10g')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def bh(p):return multipletests(p,method='fdr_bh')[1]
def blocks(f,n,p,selection):
 raw=f['raw/X'];ptr=raw['indptr'][:]
 for start in range(0,n,2000):
  end=min(n,start+2000);mask=selection[start:end]
  if not mask.any():continue
  lo,hi=ptr[start],ptr[end]
  x=sparse.csr_matrix((raw['data'][lo:hi].astype(float),raw['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,p))[mask]
  yield np.flatnonzero(mask)+start,x
def main():
 R=Path(sys.argv[1]);assert (R/'.running').exists();P=R/'public';D=R/'private';S=R/'source'
 for d in [P,D,S]:d.mkdir(exist_ok=True)
 spec=dict(version='LYPLA1_highlow_v1',unit='cell_exploratory_no_donor_adjustment',cohorts=['Wu2021','Pal2021_reprocessed'],cell_type='author/reprocessed malignant epithelial only',grouping='within each cohort, normalized LYPLA1>0 detected; high>median of detected, low=(0,median], zero=0; ties not split',contrasts=['high_vs_low_detected','high_vs_zero'],DE='Welch unequal-variance two-sided t on per-cell log1p(CP10k); sparse sufficient statistics',filter='unique nonempty symbol; detected in >=1% of all selected malignant cells; nonzero variance in contrast; exclude LYPLA1 from testing and enrichment',effect='mean log1p difference AND log2((mean CP10k high+0.1)/(mean CP10k ref+0.1)); latter descriptive regularized ratio',GSEA='Reactome_2022; complete tested gene ranking by signed Welch t; sizes15..500; weight1;1000 gene-set permutations;seed20260925',ORA='GO_Biological_Process_2023 and Reactome_2022; up/down P<0.05 and absolute descriptive log2ratio>=0.25; tested genes universe; sizes15..500 after intersecting universe; hypergeometric right tail',multiple_testing='DE BH per contrast and joint4; GSEA package empirical FDR retained; ORA BH across both directions and libraries per contrast plus joint4; nominal P exploratory',seed=20260925,LYPLA1_excluded=True,limitations='Cell dependence not modelled; no patient-level or causal inference; pathway enrichment is transcriptional association, not metabolic flux',software=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,gseapy=gseapy.__version__))
 (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));manifest=[];libs={}
 for name in ['Reactome_2022','GO_Biological_Process_2023']:
  url='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName='+name
  z=requests.get(url,timeout=90);z.raise_for_status();f=S/(name+'.gmt');f.write_bytes(z.content)
  lib={}
  for line in z.text.splitlines():
   a=line.split('\t')
   if len(a)>3:lib[a[0]]=sorted(set(x.split(',')[0] for x in a[2:] if x)-{'LYPLA1'})
  assert len(lib)>500;libs[name]=lib;manifest.append(dict(path=str(f),url=url,sha256=sha(f)))
 (P/'libraries.json').write_text(json.dumps({k:len(v) for k,v in libs.items()},indent=2))
 hp=PREV/'brca_sc117_profile_v1.py';hs=importlib.util.spec_from_file_location('helper',hp);h=importlib.util.module_from_spec(hs);hs.loader.exec_module(h)
 group_rows=[];qc=[];composition=[];audits=[];des=[];es=[];oras=[]
 for co,cn in [('Wu2021','brca_sc117_Wu2021_config.json'),('Pal2021_reprocessed','brca_sc117_Pal2021_config.json')]:
  print('START',co,flush=True);cf=PREV/'source'/cn;cfg=json.loads(cf.read_text());src=BASE/'20260919T140105Z_scRNA117_v1/source'/cfg['file']
  with h5py.File(src,'r') as f:
   obs=h.dataframe(f['obs']);var=h.dataframe(f['raw/var']);names=var[cfg['symbol_column']].astype(str).to_numpy();n,p=len(obs),len(var);sel=obs[cfg['celltype_column']].map(cfg['celltype_map']).eq('Malignant_epithelial').to_numpy()
   if cfg.get('filter_column'):sel &= obs[cfg['filter_column']].eq(cfg['filter_value']).to_numpy()
   assert (names=='LYPLA1').sum()==1;target=int(np.flatnonzero(names=='LYPLA1')[0]);ly=np.full(n,np.nan);depth=np.full(n,np.nan);ng=np.full(n,np.nan)
   for ii,x in blocks(f,n,p,sel):
    lib=np.asarray(x.sum(1)).ravel();assert (lib>0).all();ly[ii]=np.log1p(x[:,target].toarray().ravel()/lib*10000);depth[ii]=lib;ng[ii]=np.asarray((x>0).sum(1)).ravel()
   cut=float(np.median(ly[sel & (ly>0)]));group=np.full(n,-1);group[sel & (ly>cut)]=0;group[sel & (ly>0)&(ly<=cut)]=1;group[sel & (ly==0)]=2;nc=np.array([(group==i).sum() for i in range(3)]);assert (nc>=50).all() and nc.sum()==sel.sum()
   donor=obs[cfg['donor_column']].astype(str);private=obs.loc[sel,[cfg['donor_column']]].copy();private['group']=group[sel];private['LYPLA1_log1p10k']=ly[sel];private['counts']=depth[sel];private['n_genes']=ng[sel];private.to_csv(D/(co+'_cell_groups.tsv'),sep='\t')
   ct=pd.crosstab(donor[sel],group[sel]);chi=stats.chi2_contingency(ct,correction=False)[0];cramer=np.sqrt(chi/(nc.sum()*min(ct.shape[0]-1,ct.shape[1]-1)))
   for i,label in enumerate(['high','low_detected','zero']):
    mask=group==i;dc=donor[mask].value_counts();pr=dc/dc.sum();group_rows.append(dict(cohort=co,group=label,n_cells=int(nc[i]),positive_median_cut=cut,median_LYPLA1=float(np.median(ly[mask])),mean_LYPLA1=float(np.mean(ly[mask])),n_donors=len(dc),largest_donor_fraction=float(pr.max()),effective_donor_number=float(1/(pr**2).sum()),donor_group_CramersV=cramer,median_counts=np.median(depth[mask]),median_detected_genes=np.median(ng[mask])))
    for field in [cfg.get('stratum_column'),cfg.get('treatment_column')]:
     if field:
      for label2,num in obs.loc[mask,field].astype(str).value_counts().items():composition.append(dict(cohort=co,group=label,field=field,label=label2,n_cells=int(num),fraction=num/nc[i]))
   for j,contrast in [(1,'high_vs_low_detected'),(2,'high_vs_zero')]:
    for metric,values in [('log1p_total_counts',np.log1p(depth)),('detected_genes',ng)]:
     a,b=values[group==0],values[group==j];qc.append(dict(cohort=co,contrast=contrast,metric=metric,mean_high=a.mean(),mean_reference=b.mean(),standardized_mean_difference=(a.mean()-b.mean())/np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)))
   print('GROUPS',co,nc.tolist(),'cut',cut,flush=True)
   sums=np.zeros((3,p));squares=np.zeros_like(sums);linear=np.zeros_like(sums);detect=np.zeros_like(sums)
   audit_names=['LYPLA1','LYPLA2','SLC6A6','GPCPD1','PCYT2','ETNK1','ACTB','GAPDH','MKI67','EPCAM'];auditix=[int(np.flatnonzero(names==g)[0]) for g in audit_names if (names==g).sum()==1];auditx=[];auditg=[]
   for ii,x in blocks(f,n,p,sel):
    norm=x.multiply((10000/depth[ii])[:,None]).tocsr();y=norm.copy();y.data=np.log1p(y.data);gr=group[ii]
    auditx.append(y[:,auditix].toarray());auditg.extend(gr)
    for i in range(3):
     q=y[gr==i];sums[i]+=np.asarray(q.sum(0)).ravel();squares[i]+=np.asarray(q.multiply(q).sum(0)).ravel();linear[i]+=np.asarray(norm[gr==i].sum(0)).ravel();detect[i]+=np.asarray((q>0).sum(0)).ravel()
   ax=np.vstack(auditx);ag=np.array(auditg);mu=sums/nc[:,None];va=np.maximum(0,(squares-sums*sums/nc[:,None])/(nc[:,None]-1));lin=linear/nc[:,None];pct=detect/nc[:,None]
   unique=~pd.Series(names).duplicated(keep=False).to_numpy();baseok=unique & (names!='') & (names!='nan') & (names!='LYPLA1') & (detect.sum(0)/nc.sum()>=.01)
   save(pd.DataFrame(dict(gene=names,unique_symbol=unique,all_malignant_detection=detect.sum(0)/nc.sum(),eligible_base=baseok)),P/(co+'_gene_coverage.tsv'))
   for j,contrast in [(1,'high_vs_low_detected'),(2,'high_vs_zero')]:
    se2=va[0]/nc[0]+va[j]/nc[j];ok=baseok & (se2>0);se=np.sqrt(se2[ok]);t=(mu[0,ok]-mu[j,ok])/se;df=se2[ok]**2/((va[0,ok]/nc[0])**2/(nc[0]-1)+(va[j,ok]/nc[j])**2/(nc[j]-1));pv=2*stats.t.sf(abs(t),df)
    de=pd.DataFrame(dict(cohort=co,contrast=contrast,gene=names[ok],n_high=nc[0],n_reference=nc[j],mean_log_high=mu[0,ok],mean_log_reference=mu[j,ok],mean_log_difference=mu[0,ok]-mu[j,ok],mean_CP10k_high=lin[0,ok],mean_CP10k_reference=lin[j,ok],descriptive_log2ratio=np.log2((lin[0,ok]+.1)/(lin[j,ok]+.1)),detection_high=pct[0,ok],detection_reference=pct[j,ok],welch_t=t,welch_df=df,p_value=pv,q_contrast=bh(pv)))
    assert de.gene.is_unique and 'LYPLA1' not in set(de.gene);des.append(de);stem=co+'__'+contrast
    save(de,P/(stem+'_DE.tsv'))
    for k,ix in enumerate(auditix):
     if not ok[ix]:continue
     tt,pp=stats.ttest_ind(ax[ag==0,k],ax[ag==j,k],equal_var=False);r=de[de.gene.eq(names[ix])].iloc[0];assert np.isclose(tt,r.welch_t,rtol=1e-7,atol=1e-9);assert np.isclose(pp,r.p_value,rtol=1e-6,atol=1e-200);audits.append(dict(cohort=co,contrast=contrast,gene=names[ix],status='PASS',welch_t=tt))
    rank=de[['gene','welch_t']].sort_values(['welch_t','gene'],ascending=[False,True]);pre=gseapy.prerank(rnk=rank,gene_sets=libs['Reactome_2022'],outdir=str(D/(stem+'_gsea')),min_size=15,max_size=500,permutation_num=1000,weighted_score_type=1,threads=4,no_plot=True,seed=20260925,verbose=False)
    en=pre.res2d.copy();en.insert(0,'contrast',contrast);en.insert(0,'cohort',co);en['library']='Reactome_2022';en['p_resolution_note']='1000 gene-set permutations; reported zero is below permutation resolution, not exact zero';save(en,P/(stem+'_GSEA.tsv'));es.append(en)
    uni=set(de.gene);ora=[]
    for direction in ['up','down']:
     effect=de.descriptive_log2ratio;hit=set(de.loc[(de.p_value<.05)&((effect>=.25) if direction=='up' else (effect<=-.25)),'gene'])
     for libname,lib in libs.items():
      for term,gg in lib.items():
       g=uni.intersection(gg)
       if not 15<=len(g)<=500:continue
       overlap=hit&g;pv=stats.hypergeom.sf(len(overlap)-1,len(uni),len(g),len(hit));ora.append(dict(cohort=co,contrast=contrast,direction=direction,library=libname,term=term,universe_n=len(uni),selected_n=len(hit),term_n=len(g),overlap_n=len(overlap),p_value=pv,genes=';'.join(sorted(overlap))))
    ora=pd.DataFrame(ora);ora['q_contrast']=bh(ora.p_value);save(ora,P/(stem+'_ORA.tsv'));oras.append(ora)
    print('DONE',stem,'genes',len(de),'GSEA',len(en),'ORA',len(ora),flush=True)
  manifest.extend([dict(path=str(cf),url='',sha256=sha(cf)),dict(path=str(src),url='',sha256=sha(src))])
 save(pd.DataFrame(group_rows),P/'group_summary.tsv');save(pd.DataFrame(qc),P/'depth_confounding.tsv');save(pd.DataFrame(composition),P/'subtype_treatment_composition.tsv');save(pd.DataFrame(audits),P/'independent_DE_audit.tsv')
 allDE=pd.concat(des,ignore_index=True);allDE['q_joint4']=bh(allDE.p_value)
 for (co,ct),d in allDE.groupby(['cohort','contrast']):save(d,P/(co+'__'+ct+'_DE.tsv'))
 allORA=pd.concat(oras,ignore_index=True);allORA['q_joint4']=bh(allORA.p_value)
 for (co,ct),d in allORA.groupby(['cohort','contrast']):save(d,P/(co+'__'+ct+'_ORA.tsv'))
 gsea=pd.concat(es,ignore_index=True);save(gsea,P/'GSEA_all4.tsv')
 a=allDE[allDE.cohort.eq('Wu2021')&allDE.contrast.eq('high_vs_low_detected')];b=allDE[allDE.cohort.eq('Pal2021_reprocessed')&allDE.contrast.eq('high_vs_low_detected')]
 cross=a[['gene','mean_log_difference','descriptive_log2ratio','p_value','q_contrast','welch_t']].merge(b[['gene','mean_log_difference','descriptive_log2ratio','p_value','q_contrast','welch_t']],on='gene',suffixes=('_Wu','_Pal'));cross['same_direction']=np.sign(cross.mean_log_difference_Wu)==np.sign(cross.mean_log_difference_Pal);save(cross,P/'cross_cohort_primary_genes.tsv')
 manifest.extend([dict(path=str(hp),url='',sha256=sha(hp)),dict(path=str(Path(__file__)),url='',sha256=sha(__file__))]);save(pd.DataFrame(manifest),P/'source_manifest.tsv')
 (P/'validation.json').write_text(json.dumps(dict(status='PASS',four_contrasts=True,independent_DE_checks=len(audits),no_target_in_enrichment=True,source_cells_only_server=True,limitations='Cell-level exploratory, no donor adjustment, no causal interpretation'),indent=2))
 (R/'NUMERICAL_DONE').write_text('DONE');print('ALL_DONE',flush=True)
if __name__=='__main__':main()
