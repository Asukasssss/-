"""Locked PRAD24 companions evaluated in author's Hirz2023 Tumor cells."""
import argparse,json,hashlib,platform,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy import io,stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prad_lypla1_states_associations_v1 import correlate,sha,save
p=argparse.ArgumentParser();p.add_argument('--source-run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--commit',required=True);a=p.parse_args()
r=a.out;r.mkdir(exist_ok=False);(r/'.running').write_text(r.name);pub=r/'public';pub.mkdir();priv=r/'private';priv.mkdir()
s=a.source_run/'private/hirz';o=pd.read_csv(s/'cells.tsv',sep='\t');genes=pd.read_csv(s/'genes.txt',header=None)[0].astype(str).to_numpy();x=io.mmread(s/'counts.mtx').T.tocsr()
assert x.shape==(len(o),len(genes)) and o.cell.is_unique and len(set(genes))==len(genes) and (o.author_label=='Tumor').all()
clinical=pd.read_excel(a.source_run/'source/supplementary_data1.xlsx',sheet_name='Prostatectomy cases',header=1);clinical.columns=clinical.columns.str.strip()
# Match full expected source sample names to explicit clinical fname and grade, not numeric ordering.
lookup={}
for _,q in clinical.iterrows():
    if pd.isna(q['fname']) or q['Tumor']!='Y':continue
    grade={'low_grade':'LG','high_grade':'HG'}.get(q['grade'])
    if grade:lookup['SCG-'+str(q['fname'])+'-T-'+grade]=str(q['fname'])
unmatched=sorted(set(o['sample'])-set(lookup))
excluded=o[o['sample'].isin(unmatched)].copy();save(excluded,priv/'excluded_metadata_mismatch.tsv')
if len(excluded):
    keep=~o['sample'].isin(unmatched);x=x[keep.to_numpy()];o=o.loc[keep].reset_index(drop=True)
assert set(o['sample']).issubset(lookup) and len(o)>0
o['donor']=o['sample'].map(lookup);assert o.donor.notna().all()
ly=np.flatnonzero(genes=='LYPLA1');assert len(ly)==1;ly=ly[0]
lib=o.library.to_numpy();assert np.isfinite(lib).all() and (lib>0).all()
# Use pre-export FULL source-gene library denominator, not intersected-gene library sum.
z=np.log1p(x.toarray()*10000/lib[:,None]);y=z[:,ly]
o['LYPLA1']=y;save(o,priv/'malignant_cells.tsv.gz')
base=pd.read_csv(a.source_run/'public/results.tsv.gz',sep='\t');base=base[base.p_value.notna()].copy()
assert base.stable_gene_id.is_unique
base['locked_positive']=(base.p_value<.05)&(base.mean_rho>0);base['locked_negative']=(base.p_value<.05)&(base.mean_rho<0)
symbols=base.gene.value_counts();mapping={g:i for i,g in enumerate(genes)}
selected={g:mapping[g] for g in base.gene if g in mapping and symbols[g]==1}
modes=['ordinary','technical_adjusted','positive_only'];arrays={m:[] for m in modes};donors=[];qc=[]
for d,g in o.groupby('donor',sort=True):
    ix=g.index.to_numpy();yy=y[ix];pos=yy>0
    eligible=len(ix)>=30 and pos.sum()>=10 and np.ptp(yy)>0
    qc.append(dict(donor=d,n_cells=len(ix),n_positive=int(pos.sum()),evaluable=eligible))
    if not eligible:continue
    donors.append(d);cov=np.column_stack([np.ones(len(ix)),stats.rankdata(np.log1p(lib[ix])),stats.rankdata(o.percent_mito.to_numpy()[ix])]);cov[:,1:]=(cov[:,1:]-cov[:,1:].mean(0))/(cov[:,1:].std(0)+1e-12)
    per={m:np.full(len(genes),np.nan) for m in modes};detect=np.asarray((x[ix]>0).sum(0)).ravel()>=max(10,int(np.ceil(.1*len(ix))))
    for k in range(0,len(genes),512):
        sl=slice(k,min(k+512,len(genes)));zz=z[ix,sl]
        per['ordinary'][sl]=correlate(zz,yy,np.ones((len(ix),1)))
        per['technical_adjusted'][sl]=correlate(zz,yy,cov)
        if pos.sum()>=30:
            pe=(zz[pos]>0).sum(0)>=max(10,int(np.ceil(.1*pos.sum())))
            q=correlate(zz[pos],yy[pos],cov[pos]);q[~pe]=np.nan;per['positive_only'][sl]=q
    for m in modes:per[m][~detect]=np.nan;arrays[m].append(per[m])
save(pd.DataFrame(qc),priv/'donor_qc.tsv');rows=[]
for mode,ls in arrays.items():
    mat=np.asarray(ls);save(pd.DataFrame(mat,index=donors,columns=genes).rename_axis('donor').reset_index(),priv/(mode+'_rho.tsv.gz'))
    for _,b in base.iterrows():
        mapped=b.gene in selected;rr=mat[:,selected[b.gene]] if mapped else np.array([]);rr=rr[np.isfinite(rr)];n=len(rr);nz=rr[rr!=0];status='DONE' if n>=8 and len(nz) else 'NOT_EVALUABLE'
        rows.append(dict(cancer='PRAD',cohort='Hirz2023_GSE181294',stage_id='06_EXTERNAL',run_id=r.name,analysis_version='prad_lypla1_hirz_replication_v1',analysis_type='locked_companion_cross_cohort',metabolite_key='NA',metabolite_name='NA',gene=b.gene,unit='clinically_mapped_patient',n=n,n_reference=len(o.donor.unique()),effect_type='mean_within_patient_rank_correlation',effect=rr.mean() if n else np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=stats.binomtest(int((nz>0).sum()),len(nz)).pvalue if status=='DONE' else np.nan,q_value=np.nan,test_family='locked_all_PRAD24_evaluable_genes_'+mode,family_n_evaluable=0,status=status,reason='' if status=='DONE' else 'missing_or_nonunique_symbol' if not mapped else 'fewer_than8_evaluable_patients',source_id='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181294',mode=mode,stable_gene_id=b.stable_gene_id,n_positive=int((rr>0).sum()),n_negative=int((rr<0).sum()),discovery_n=int(b.n_donors),discovery_mean_rho=b.mean_rho,discovery_p=b.p_value,discovery_q=b.q_value,locked_positive=bool(b.locked_positive),locked_negative=bool(b.locked_negative),technical_gene=bool(b.technical_gene),direction_agrees=bool(np.sign(rr.mean())==np.sign(b.mean_rho)) if n else None))
res=pd.DataFrame(rows)
for mode in modes:
    ok=(res['mode']==mode)&res.p_value.notna();res.loc[res['mode']==mode,'family_n_evaluable']=int(ok.sum())
    if ok.any():res.loc[ok,'q_value']=multipletests(res.loc[ok,'p_value'],method='fdr_bh')[1]
save(res,pub/'results.tsv.gz')
focus=['STEAP2','CMPK1','RDH11','FOXA1','KLK3','APLP2','TMEM230','EDEM3','IPO7','IDH1','PDLIM5','ENDOD1']
summary=res[(res['mode']=='technical_adjusted') & res.gene.isin(focus)].copy();summary['order']=summary.gene.map({g:i for i,g in enumerate(focus)});summary=summary.sort_values('order').drop(columns='order');save(summary,pub/'priority_gene_comparison.tsv')
desc=res[(res['mode']=='technical_adjusted')&res.locked_positive&(res.n>=3)]
val=dict(status='DONE',formal_replication_status='NOT_EVALUABLE' if res.p_value.notna().sum()==0 else 'DONE',n_author_tumor_cells=len(o),n_excluded_metadata_mismatch=len(excluded),n_excluded_samples_metadata_mismatch=len(unmatched),n_patients_with_any_tumor_cells=o.donor.nunique(),n_patients_ge30=int(sum(q['n_cells']>=30 for q in qc)),n_patients_association_eligible=len(donors),n_discovery_genes=len(base),n_locked_positive=int(base.locked_positive.sum()),n_locked_positive_descriptive_ge3=len(desc),n_locked_positive_same_mean_direction_ge3=int(desc.direction_agrees.sum()),no_new_subtypes=True,subtype_replication='NOT_EVALUABLE: author Tumor label has no epithelial subtype annotation',clinical_sample_mapping_verified=True,barcode_join_verified=True,independence='Different named study; no cross-study donor genotype/identity crosswalk, so no verified non-overlap claim')
(pub/'validation.json').write_text(json.dumps(val,indent=2))
spec=dict(version='prad_lypla1_hirz_replication_v1',parent_commit=a.commit,discovery_results=str(a.source_run/'public/results.tsv.gz'),gene_selection='All 9057 discovery-evaluable genes locked; exact unique symbol join; positive candidates flag inherited',min_cells=30,min_LYPLA1_positive=10,min_gene_positive=10,min_gene_detection_fraction=.1,min_patients_inference=8,min_patients_descriptive=3,normalization='log1p(CP10K) with totalUMI from complete original sample matrix',methods=modes,primary='technical_adjusted',technical_covariates=['rank log1p totalUMI','rank percent_mito'],annotation='author cell1=Tumor; only clinical tumor samples',patient_mapping='Exact SCG-{clinical fname}-T-{clinical grade code} match to Supplementary Data1, each fname unique',no_pairing_inferred=True,thresholds_unchanged=True,source_article='https://www.nature.com/articles/s41467-023-36325-2',author_repository='https://github.com/shenglinmei/ProstateCancerAnalysis',software=dict(python=platform.python_version(),pandas=pd.__version__,numpy=np.__version__))
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
files=[a.source_run/'source/hirz_scrna_counts.rds',a.source_run/'source/hirz_scrna_annotations.rds',a.source_run/'source/supplementary_data1.xlsx',a.source_run/'public/results.tsv.gz',Path(__file__),Path(__file__).with_name('prad_lypla1_states_associations_v1.py'),a.source_run/'prad_hirz_scrna_extract_v1.R']
save(pd.DataFrame([dict(path=str(f),sha256=sha(f)) for f in files]),pub/'source_manifest.tsv')
fig,ax=plt.subplots(figsize=(8,6));ypos=np.arange(len(summary));ax.barh(ypos+.17,summary.discovery_mean_rho,height=.32,color='#ac3154',label='PRAD24 discovery');ax.barh(ypos-.17,summary.effect,height=.32,color='#318d88',label='Hirz2023 descriptive');ax.set_yticks(ypos,summary.gene);ax.invert_yaxis();ax.axvline(0,color='gray',lw=.7);ax.set_xlabel('Mean within-patient adjusted rank correlation');ax.legend();ax.set_title('LYPLA1 companions across PRAD cohorts\nHirz patient coverage below preset inference threshold')
for k,(_,b) in enumerate(summary.iterrows()):ax.text(max(b.effect,b.discovery_mean_rho)+.005,k,f'{b.n_positive}/{b.n} Hirz positive',va='center',fontsize=8)
ax.set_xlim(min(-.04,float(summary.effect.min())-.02),max(.36,float(summary.effect.max())+.15));fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(pub/('LYPLA1_cross_cohort.'+ext),dpi=240,bbox_inches='tight')
print(json.dumps(val),flush=True);print(summary[['gene','n','n_positive','effect','discovery_mean_rho']].to_string(index=False),flush=True)
