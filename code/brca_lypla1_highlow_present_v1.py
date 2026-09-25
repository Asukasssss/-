"""Present all four frozen comparisons, cross-cohort signals and confounding."""
from pathlib import Path
import sys,json,textwrap
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(sys.argv[1]);P=R/'public';assert (R/'NUMERICAL_DONE').exists()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.10g')
plt.rcParams.update({'pdf.fonttype':42,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
cos=['Wu2021','Pal2021_reprocessed'];con='high_vs_low_detected';summary=[];coverage=[]
reactome={}
for line in (R/'source/Reactome_2022.gmt').read_text().splitlines():
 a=line.split('\t');reactome[a[0]]=set(x.split(',')[0] for x in a[2:] if x)-{'LYPLA1'}
groups=pd.read_csv(P/'group_summary.tsv',sep='\t')
for co in cos:
 cells=pd.read_csv(R/'private'/(co+'_cell_groups.tsv'),sep='\t')
 counts=np.expm1(cells.LYPLA1_log1p10k)*cells['counts']/10000
 for i,label in enumerate(['high','low_detected','zero']):
  mask=(groups.cohort==co)&(groups.group==label)
  groups.loc[mask,'median_raw_LYPLA1_counts']=float(np.median(counts[cells.group==i]))
save(groups,P/'group_summary.tsv')
for co in cos:
 for ct in [con,'high_vs_zero']:
  d=pd.read_csv(P/(co+'__'+ct+'_DE.tsv'),sep='\t');en=pd.read_csv(P/(co+'__'+ct+'_GSEA.tsv'),sep='\t')
  universe=set(d.gene)
  for term,genes in reactome.items():
   n=len(genes&universe);coverage.append(dict(cohort=co,contrast=ct,term=term,library_genes=len(genes),tested_overlap=n,status='DONE' if 15<=n<=500 else 'NOT_EVALUABLE',reason='size15_to500' if 15<=n<=500 else 'outside_fixed_size15_to500'))
  summary.append(dict(cohort=co,contrast=ct,n_high=int(d.n_high.iloc[0]),n_reference=int(d.n_reference.iloc[0]),tested_genes=len(d),P_lt05=int((d.p_value<.05).sum()),q_contrast_lt05=int((d.q_contrast<.05).sum()),up_P_effect=int(((d.p_value<.05)&(d.descriptive_log2ratio>=.25)).sum()),down_P_effect=int(((d.p_value<.05)&(d.descriptive_log2ratio<=-.25)).sum()),GSEA_terms=len(en),GSEA_nominalP_lt05=int((en['NOM p-val']<.05).sum()),GSEA_empiricalFDR_lt025=int((en['FDR q-val']<.25).sum())))
save(pd.DataFrame(summary),P/'analysis_summary.tsv')
save(pd.DataFrame(coverage),P/'GSEA_pathway_coverage.tsv')
es=pd.read_csv(P/'GSEA_all4.tsv',sep='\t');pri=es[es.contrast.eq(con)];a=pri[pri.cohort.eq(cos[0])];b=pri[pri.cohort.eq(cos[1])]
cols=['Term','NES','NOM p-val','FDR q-val','Lead_genes']
cross=a[cols].merge(b[cols],on='Term',suffixes=('_Wu','_Pal'));cross['same_direction']=np.sign(cross.NES_Wu)==np.sign(cross.NES_Pal);cross['both_nominal_P_lt05']=(cross['NOM p-val_Wu']<.05)&(cross['NOM p-val_Pal']<.05)
cross['shared_leading_genes']=[';'.join(sorted(set(str(x).split(';'))&set(str(y).split(';')))) for x,y in zip(cross.Lead_genes_Wu,cross.Lead_genes_Pal)]
cross['min_abs_NES']=np.minimum(abs(cross.NES_Wu),abs(cross.NES_Pal));cross=cross.sort_values(['same_direction','both_nominal_P_lt05','min_abs_NES'],ascending=False);save(cross,P/'cross_cohort_primary_pathways.tsv')
lipid=es[es.Term.str.contains('lipid|phosphatid|phospholip|acyl chain|cholesterol|glycerophosph|sphingo',case=False,regex=True)].copy();save(lipid,P/'lipid_named_pathways_all4.tsv')
fig,axes=plt.subplots(1,2,figsize=(11,4.7),constrained_layout=True)
for ax,co in zip(axes,cos):
 d=pd.read_csv(P/(co+'__'+con+'_DE.tsv'),sep='\t');y=-np.log10(d.p_value.clip(lower=1e-300));c=np.where((d.p_value<.05)&(d.descriptive_log2ratio>=.25),'#CC6677',np.where((d.p_value<.05)&(d.descriptive_log2ratio<=-.25),'#4477AA','#BBBBBB'))
 ax.scatter(d.descriptive_log2ratio,y,c=c,s=4,alpha=.5,rasterized=True);ax.axhline(-np.log10(.05),color='grey',ls='--',lw=.6);ax.set_title(co);ax.set_xlabel('Descriptive log2 ratio (mean CP10k + 0.1)');ax.set_ylabel('-log10(cell-level P), capped at 300')
 topgenes=d.assign(abst=d.welch_t.abs()).nlargest(5,'abst');ax.text(.98,.97,'Top |t| genes:\n'+'\n'.join(topgenes.gene),transform=ax.transAxes,ha='right',va='top',fontsize=8,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
fig.suptitle('LYPLA1 high vs low among detected malignant cells | exploratory')
fig.savefig(P/'01_DE_primary.png',dpi=170);fig.savefig(P/'01_DE_primary.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(17,12),constrained_layout=True)
for ax,co in zip(axes,cos):
 d=pri[pri.cohort.eq(co)].copy();d['absNES']=d.NES.abs();top=pd.concat([d[d.NES>0].sort_values(['FDR q-val','absNES'],ascending=[True,False]).head(8),d[d.NES<0].sort_values(['FDR q-val','absNES'],ascending=[True,False]).head(8)]).sort_values('NES')
 labels=['\n'.join(textwrap.wrap(x,46))+'\nP='+('below permutation resolution' if pv==0 else f'{pv:.3g}') for x,pv in zip(top.Term,top['NOM p-val'])];ax.barh(np.arange(len(top)),top.NES,color=np.where(top['NOM p-val']>=.05,'#BBBBBB',np.where(top.NES>0,'#CC6677','#4477AA')));ax.set_yticks(np.arange(len(top)),labels,fontsize=7);ax.axvline(0,color='black',lw=.6);ax.set_xlabel('GSEA normalized enrichment score');ax.set_title(co+' | up to 8 terms per direction')
fig.suptitle('Reactome whole-ranked-list enrichment | high vs low detected\nDisplayed by empirical FDR then |NES|; grey: nominal P>=0.05; full terms and leading genes archived')
fig.savefig(P/'02_GSEA_primary.png',dpi=150);fig.savefig(P/'02_GSEA_primary.pdf');plt.close(fig)
gc=pd.read_csv(P/'cross_cohort_primary_genes.tsv',sep='\t');s=gc[gc.same_direction&(gc.p_value_Wu<.05)&(gc.p_value_Pal<.05)].copy();s['min_abs_t']=np.minimum(abs(s.welch_t_Wu),abs(s.welch_t_Pal));s=s.nlargest(24,'min_abs_t');save(s,P/'shared_primary_genes_display.tsv')
fig,ax=plt.subplots(figsize=(6,9),constrained_layout=True);v=s[['descriptive_log2ratio_Wu','descriptive_log2ratio_Pal']].to_numpy();limit=max(.1,float(np.max(abs(v)))) if len(v) else 1
if len(v):
 im=ax.imshow(v,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_yticks(range(len(s)),s.gene);ax.set_xticks([0,1],['Wu','Pal']);fig.colorbar(im,ax=ax,label='Descriptive log2 ratio (+0.1 CP10k)')
ax.set_title('Up to 24 shared-direction genes\nP<0.05 in both; largest minimum |t|');fig.savefig(P/'03_shared_genes.png',dpi=180);fig.savefig(P/'03_shared_genes.pdf');plt.close(fig)
q=pd.read_csv(P/'depth_confounding.tsv',sep='\t');fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True);y=np.arange(len(q));ax.barh(y,q.standardized_mean_difference,color='#4477AA');ax.set_yticks(y,[f'{a} | {b}\n{c}' for a,b,c in zip(q.cohort,q.contrast,q.metric)],fontsize=7);ax.axvline(0,color='black',lw=.5);ax.set_xlabel('Standardized mean difference: high minus reference');ax.set_title('Technical-depth imbalance (descriptive, not adjusted)');fig.savefig(P/'04_depth_imbalance.png',dpi=160);fig.savefig(P/'04_depth_imbalance.pdf');plt.close(fig)
print(pd.DataFrame(summary).to_string(index=False));print('SHARED PATHWAYS',int((cross.same_direction&cross.both_nominal_P_lt05).sum()));print(cross[['Term','NES_Wu','NES_Pal','NOM p-val_Wu','NOM p-val_Pal']].head(15).to_string(index=False));print(lipid[lipid.contrast.eq(con)][['cohort','Term','NES','NOM p-val','FDR q-val']].sort_values('NOM p-val').head(20).to_string(index=False));print(q.to_string(index=False))
