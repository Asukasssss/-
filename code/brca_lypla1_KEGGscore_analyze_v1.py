"""Evaluate eight prespecified transcript-program links after membership freeze."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,csv,hashlib,importlib.util
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
BASE=ROOT/'results/collaborative/BRCA/A'
VERSION='lypla1_KEGGscore_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()

def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def seed(s):return int(hashlib.sha256((VERSION+'|'+s).encode()).hexdigest()[:8],16)

def partial(x,y,cov):
 z=np.column_stack([np.ones(len(x)),cov[:,0],stats.rankdata(cov[:,1]),stats.rankdata(cov[:,2])])
 z[:,1:]=(z[:,1:]-z[:,1:].mean(0))/z[:,1:].std(0)
 if not np.isfinite(z).all() or np.linalg.matrix_rank(z)!=4:raise ValueError('degenerate_covariate_design')
 q=np.linalg.qr(z,mode='reduced')[0]
 xx=stats.rankdata(x);yy=stats.rankdata(y);ex=xx-q@(q.T@xx);ey=yy-q@(q.T@yy)
 den=np.linalg.norm(ex)*np.linalg.norm(ey)
 if den<1e-12:raise ValueError('constant_residual')
 return float(ex@ey/den),q,ex,ey

def adjusted(x,y,cov,sd):
 rho,q,ex,ey=partial(x,y,cov);rng=np.random.default_rng(sd)
 groups=[np.flatnonzero(cov[:,0]==i) for i in np.unique(cov[:,0])];extreme=0;N=99999
 for start in range(0,N,1000):
  count=min(1000,N-start);p=np.tile(ey,(count,1))
  for ix in groups:p[:,ix]=ey[ix][np.argsort(rng.random((count,len(ix))),axis=1)]
  den=np.linalg.norm(ex)*np.sqrt(np.maximum((p*p).sum(1)-((p@q)**2).sum(1),1e-20))
  r=p@ex/den;extreme+=int((abs(r)>=abs(rho)-1e-12).sum())
 pvalue=(extreme+1)/(N+1);boot=[];rng=np.random.default_rng(sd+1)
 for _ in range(4000):
  ix=np.concatenate([rng.choice(g,len(g),replace=True) for g in groups])
  try:boot.append(partial(x[ix],y[ix],cov[ix])[0])
  except ValueError:pass
 ci=np.quantile(boot,[.025,.975]) if len(boot)>=3200 else [np.nan,np.nan]
 return dict(effect=rho,p_value=pvalue,ci_lower=ci[0],ci_upper=ci[1],seed=sd,permutation_n=N,permutation_extreme=extreme,
             permutation_mc_se=np.sqrt(pvalue*(1-pvalue)/(N+1)),bootstrap_n=4000,bootstrap_valid=len(boot))

def plot_results(r,points,out):
 plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':10,'savefig.dpi':180,'axes.spines.top':False,'axes.spines.right':False})
 figures=[]
 fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
 for i,co in enumerate(['Wu2021','Pal2021_reprocessed']):
  for j,gene in enumerate(['LYPLA1','LYPLA2']):
   ax=axs[i,j];key=co+'|'+gene;d=points[key];row=r[r.test_id.eq(key)].iloc[0]
   ax.scatter(d['x'],d['y'],s=42,color=['#0072B2','#009E73'][j],alpha=.8,rasterized=True)
   ax.axhline(0,lw=.7,color='#BBBBBB');ax.set_xlabel(gene+' malignant epithelial pseudobulk log1p(CP10k)')
   ax.set_ylabel('Fixed Glycerophospholipid transcript score\n(mean within-study gene z-scores)')
   ax.set_title(co+' | '+gene+f' | {len(d["x"])} donors',weight='bold')
   ax.text(.03,.97,f"rho={row.effect:+.3f}; P={row.p_value:.4g}\n95% CI [{row.ci_lower:+.3f}, {row.ci_upper:+.3f}]",transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
 fig.suptitle('Malignant epithelium: LYPLA1 / LYPLA2 and fixed Glycerophospholipid score\nOne dot per source donor; score excludes both comparison genes',fontsize=14,weight='bold')
 fig.savefig(out/'sc_gene_score.png');figures.append(fig)
 fig,axs=plt.subplots(1,2,figsize=(12,4.7),constrained_layout=True)
 for j,mode in enumerate(['processed','author_available']):
  key='CAMP|'+mode;d=points[key];row=r[r.test_id.eq(key)].iloc[0];ax=axs[j]
  ax.scatter(d['x'],d['y'],s=34,color='#0072B2',alpha=.8,rasterized=True)
  ax.set_xlabel('Bulk tissue Glycerophospholipid transcript score');ax.set_ylabel('Measured LPC 16:0 (unchanged author scale)')
  ax.set_title(('Processed main' if j==0 else 'Author-available subset')+f' | n={len(d["x"])}',weight='bold')
  ax.text(.03,.97,f"rho={row.effect:+.3f}; P={row.p_value:.4g}\n95% CI [{row.ci_lower:+.3f}, {row.ci_upper:+.3f}]",transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
 fig.suptitle('CAMP: the same fixed gene set versus measured LPC 16:0\nBulk tissue score is not malignant-cell lipid abundance or metabolic flux',fontsize=14,weight='bold')
 fig.savefig(out/'CAMP_score_LPC.png');figures.append(fig)
 fig,ax=plt.subplots(figsize=(10,6),constrained_layout=True)
 for i,row in r.iterrows():
  color='#D55E00' if row.p_value<.05 else '#0072B2'
  ax.plot([row.ci_lower,row.ci_upper],[i,i],color=color,lw=2);ax.scatter(row.effect,i,c=color,s=35)
  ax.text(1.03,i,f'P={row.p_value:.4g}',va='center',transform=ax.get_yaxis_transform(),fontsize=9)
 ax.set_yticks(range(len(r)),r.test_id);ax.invert_yaxis();ax.axvline(0,color='#999999',lw=.8);ax.set_xlim(-1,1)
 ax.set_xlabel('Spearman / partial rank correlation, pointwise 95% bootstrap interval')
 ax.set_title('All eight prespecified tests | nominal P for exploratory retention\nCAMP ER_PC12 adjusts ER and composition proxy PCs; q archived in complete table',weight='bold')
 fig.savefig(out/'all8_effects.png');figures.append(fig)
 with PdfPages(out/'LYPLA_KEGGscore_results.pdf') as pdf:
  for fig in figures:pdf.savefig(fig)
 for fig in figures:plt.close(fig)

def main():
 R=Path(sys.argv[1]);P=R/'public';D=R/'private';S=R/'source'
 assert (R/'PREPARATION_DONE').exists() and not (P/'results.tsv').exists()
 frozen=json.loads((P/'membership_freeze.json').read_text());assert sha(P/'frozen_score_members.tsv')==frozen['membership_sha256']
 genes=frozen['genes'];assert not {'LYPLA1','LYPLA2'}&set(genes)
 helper=importlib.util.spec_from_file_location('statistics_helper',R/'brca_lypla1_paired_delta_v1.py');h=importlib.util.module_from_spec(helper);helper.loader.exec_module(h)
 cohort_scores={};scoremeta=[];manifest=[]
 for co in ['Wu2021','Pal2021_reprocessed','CAMP']:
  f=D/(co+'_expression.tsv');e=pd.read_csv(f,sep='\t',index_col=0);mean=e[genes].mean(0);sd=e[genes].std(0,ddof=1)
  assert np.isfinite(e[genes]).all().all() and (sd>0).all()
  z=(e[genes]-mean)/sd;score=z.mean(1)
  assert np.max(abs(z.mean(0)))<1e-10 and np.max(abs(z.std(0,ddof=1)-1))<1e-10
  save(z.reset_index(),D/(co+'_zscore_matrix.tsv'))
  e['PC_score']=score;e.to_csv(D/(co+'_scored_expression.tsv'),sep='\t');cohort_scores[co]=e
  for g in genes:scoremeta.append(dict(cohort=co,gene=g,mean=mean[g],sd=sd[g],weight=1/len(genes)))
  manifest.append(dict(path=str(f),sha256=sha(f)))
 save(pd.DataFrame(scoremeta),P/'score_standardization.tsv')
 rows=[];points={};private_rows=[]
 def test(co,typ,x,y,labels,gene='NA',cov=None):
  key=co+'|'+typ;nn=len(x);row=dict.fromkeys(PREFIX,np.nan)
  row.update(cancer='BRCA',cohort=co,stage_id='06_EXTERNAL' if co!='CAMP' else ('04_ROBUSTNESS' if cov is not None else '03_PATIENT'),run_id=R.name,
             analysis_version=VERSION,analysis_type=typ,metabolite_key='KEGG:C04102' if co=='CAMP' else 'NA',metabolite_name='1-palmitoyl-GPC (16:0)' if co=='CAMP' else 'NA',
             gene=gene,unit='publisher_donor_label' if co!='CAMP' else 'audited_author_tumor_case',n=nn,n_reference=nn,
             effect_type='partial_rank_correlation' if cov is not None else 'Spearman_rho',test_family='KEGGscore_BH8',family_n_evaluable=8,
             status='NOT_EVALUABLE',reason='less_than8_or_constant',source_id='frozen_KEGG_hsa00564',test_id=key,score_gene_n=len(genes))
  x=np.asarray(x,float);y=np.asarray(y,float)
  if nn>=8 and np.ptp(x)>0 and np.ptp(y)>0:
   if cov is not None:assert nn>=20 and min(np.bincount(cov[:,0].astype(int)))>=4
   sd=seed(key)
   row.update(adjusted(x,y,cov,sd) if cov is not None else h.calculate(x,y,sd,99999,4000))
   row.update(status='DONE',reason='same_fixed_transcript_program;not_metabolic_activity_or_causality')
  rows.append(row);points[key]=dict(x=x,y=y)
  for i,label in enumerate(labels):
   q=dict(test_id=key,source_unit=str(label),x=x[i],y=y[i])
   if cov is not None:q.update(ER=cov[i,0],PC1=cov[i,1],PC2=cov[i,2])
   private_rows.append(q)
  print(key,nn,row['effect'],row['p_value'],flush=True)
 for co in ['Wu2021','Pal2021_reprocessed']:
  e=cohort_scores[co]
  for gene in ['LYPLA1','LYPLA2']:test(co,gene,e[gene],e.PC_score,e.index,gene)
 m=pd.read_csv(D/'CAMP_audited_tumors.tsv',sep='\t',dtype=str);e=cohort_scores['CAMP'];e.index=e.index.astype(str)
 assert list(e.index)==m.case_row.tolist()
 src=ROOT/'data/candidates/camp_primary_tissue_multicancer';mf=src/'processed_metabolomics'/m.MetabFile.iloc[0]
 mfraw=pd.read_excel(mf,sheet_name='data',index_col=0);mat=pd.read_excel(mf,sheet_name='data_imputed',index_col=0)
 for frame in [mat,mfraw]:frame.columns=frame.columns.astype(str);frame.index=frame.index.astype(str).str.strip();assert frame.index.is_unique
 y=mat.loc['1-palmitoyl-GPC (16:0)',m.MetabID].to_numpy(float)
 available=np.isfinite(mfraw.loc['1-palmitoyl-GPC (16:0)',m.MetabID].to_numpy(float))
 assert np.isfinite(y).all() and len(y)==60 and available.sum()==55
 x=e.PC_score.to_numpy()
 for typ,mask in [('processed',np.ones(60,bool)),('author_available',available)]:test('CAMP',typ,x[mask],y[mask],m.case_row[mask])
 # A single prespecified background model, excluding program and candidate genes from proxies.
 rf=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
 rna=pd.read_csv(rf,index_col=0);rna.columns=rna.columns.astype(str)
 markerpath=BASE/'20260919T151200Z_all117_robustness_v2/source/MCPcounter_genes.txt';markers=pd.read_csv(markerpath,sep='\t')
 history=set(pd.read_csv(S/'history_genes156.tsv',sep='\t').gene);assert len(history)==156
 pathway=set(pd.read_csv(S/'kegg_members.tsv',sep='\t').gene)
 excluded=history|pathway|{'LYPLA1','LYPLA2'};panels={};mc=[];md=[]
 for pop,g in markers.groupby('Cell population',sort=True):
  full=sorted(set(g['HUGO symbols'].dropna()));present=[a for a in full if a in rna.index];keep=[a for a in present if a not in excluded]
  ok=len(keep)>=2 and len(keep)/len(full)>=.5
  mc.append(dict(population=pop,original_markers=len(full),measured=len(present),removed=';'.join(sorted(set(present)&excluded)),retained=len(keep),eligible=ok))
  for a in full:md.append(dict(population=pop,gene=a,measured=a in present,excluded=a in excluded,used=ok and a in keep))
  if ok:panels[pop]=rna.loc[keep,m.RNAID].mean(0).to_numpy()
 save(pd.DataFrame(mc),P/'composition_marker_coverage.tsv');save(pd.DataFrame(md),P/'composition_marker_decisions.tsv')
 assert len(panels)>=8
 z=pd.DataFrame(panels);zz=(z-z.mean(0))/z.std(0,ddof=1);assert np.isfinite(zz).all().all()
 u,s,v=np.linalg.svd(zz,full_matrices=False);pcs=u[:,:2]*s[:2]
 save(pd.DataFrame(dict(population=z.columns,PC1_loading=v[0],PC2_loading=v[1])),P/'composition_PC_loadings.tsv')
 (P/'composition_summary.json').write_text(json.dumps(dict(panels=len(panels),PC1_variance_fraction=float(s[0]**2/sum(s*s)),PC2_variance_fraction=float(s[1]**2/sum(s*s)),interpretation='modified RNA proxies;not cell fractions;no epithelial fraction estimated'),indent=2))
 geofile=BASE/'20260919_followup_ER_v2/source/geo_metadata_lines.txt'
 lines=list(csv.reader(geofile.read_text().splitlines(),delimiter='\t'));geo=pd.DataFrame(index=lines[0][1:])
 for line in lines[1:]:geo[line[1].split(': ',1)[0]]=[v.split(': ',1)[1] for v in line[1:]]
 gsm=m.RNAID.str.replace('.CEL.gz','',regex=False)
 assert geo.loc[gsm,'tissue type'].eq('Tumor').all()
 er=geo.loc[gsm,'estrogen receptor status'].map({'Positive':1.,'Negative':0.}).to_numpy();assert np.isfinite(er).all()
 cov=np.column_stack([er,pcs])
 for typ,mask in [('processed_ER_PC12',np.ones(60,bool)),('author_available_ER_PC12',available)]:test('CAMP',typ,x[mask],y[mask],m.case_row[mask],cov=cov[mask])
 result=pd.DataFrame(rows);assert len(result)==8 and result.p_value.notna().all()
 result['q_value']=multipletests(result.p_value,method='fdr_bh')[1];result['nominal_P_lt005']=result.p_value.lt(.05)
 result=result[PREFIX+[c for c in result if c not in PREFIX]]
 save(result,P/'results.tsv');save(pd.DataFrame(private_rows),D/'audit_test_inputs.tsv')
 plot_results(result,points,P)
 with pd.ExcelWriter(P/'LYPLA_KEGGscore_results.xlsx',engine='openpyxl') as w:
  for file,title in [('results.tsv','8项完整结果'),('frozen_score_members.tsv','固定通路成员'),('gene_coverage.tsv','三平台基因覆盖'),('sc_cohort_coverage.tsv','单细胞供者范围'),('composition_marker_coverage.tsv','组成代理覆盖')]:
   pd.read_csv(P/file,sep='\t').to_excel(w,sheet_name=title,index=False)
  for ws in w.book.worksheets:
   ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
   for col in ws.columns:ws.column_dimensions[col[0].column_letter].width=23
 for f in [mf,rf,markerpath,geofile,Path(__file__),R/'brca_lypla1_paired_delta_v1.py',P/'frozen_score_members.tsv']:
  manifest.append(dict(path=str(f),sha256=sha(f)))
 save(pd.DataFrame(manifest),P/'analysis_manifest.tsv')
 val=dict(status='DONE',membership_frozen_before_correlations=True,n_common_genes=len(genes),target_comparator_excluded=True,tests=8,
          no_historical_P_q_overwritten=True,composition_panels=len(panels),patient_or_donor_numeric_tables_server_only=True,
          no_cross_cohort_patient_links=True,score_standardization_checked=True,source_row_n={'Wu':len(cohort_scores['Wu2021']),'Pal':len(cohort_scores['Pal2021_reprocessed']),'CAMP':60,'CAMP_available':55},independent_R_audit='PENDING')
 (P/'validation.json').write_text(json.dumps(val,indent=2));(R/'NUMERICAL_DONE').write_text('DONE')
 print(result[['test_id','n','effect','ci_lower','ci_upper','p_value','q_value']].to_string(index=False),flush=True)

if __name__=='__main__':main()
