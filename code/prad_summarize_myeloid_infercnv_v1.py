"""Aggregate inferred profiles without exporting cell/patient measurements."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
heatmap_cmap=plt.get_cmap('RdBu_r').copy();heatmap_cmap.set_bad('#d9d9d9')
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root;pub=r/'public/06_EXTERNAL'
meta=pd.read_csv(r/'private/selected_cell_metadata.tsv',sep='\t').set_index('cell')
plan=pd.read_csv(r/'private/donor_plan.tsv',sep='\t')
scores=[];profiles=[];checks=[];bin_info=[];cell_sens=[]
for donor in plan.donor_alias:
 primary_mat=None
 for mode in ['primary','reference_split']:
  d=r/'private'/donor/mode
  if not (d/'DONE').exists():
   checks.append(dict(donor_alias=donor,mode=mode,status='FAILED_OR_NOT_DONE'));continue
  sc=pd.read_csv(d/'cell_scores.tsv',sep='\t').set_index('cell')
  z=meta.loc[sc.index].copy();assert z.donor_alias.eq(donor).all()
  z['group']=np.where(z.role.eq('Myeloid'),'Myeloid | '+z.subtype+' | '+z.tissue,z.role)
  ref='T_A' if mode=='primary' else 'T_B';hold='T_B' if mode=='primary' else 'T_A'
  z.loc[z.role.eq(ref),'group']='T_reference';z.loc[z.role.eq(hold),'group']='T_holdout'
  sc=sc.join(z[['group']]);sc['mode']=mode;sc['donor_alias']=donor;scores.append(sc.reset_index())
  mat=pd.read_csv(d/'bin_log2_ratios.tsv.gz',sep='\t').set_index('bin_id')
  assert set(mat.columns)==set(sc.index)
  if mode=='primary':primary_mat=mat
  elif primary_mat is not None:
   for cell in z.index[z.role.eq('Myeloid')]:
    w=pd.concat([primary_mat[cell].rename('primary'),mat[cell].rename('reference_split')],axis=1).dropna()
    x=w.primary;y=w.reference_split
    rho=float(spearmanr(x,y).statistic) if len(w)>=3 and x.nunique()>1 and y.nunique()>1 else np.nan
    cell_sens.append(dict(cell=cell,donor_alias=donor,subtype=z.loc[cell,'subtype'],tissue=z.loc[cell,'tissue'],rho=rho,mean_absolute_log2_difference=float(np.mean(abs(x-y)))))
  bi=pd.read_csv(d/'bins.tsv',sep='\t');bin_info.append(bi)
  for group,cells in z.groupby('group').groups.items():
   x=mat[list(cells)].mean(axis=1)
   profiles.append(pd.DataFrame(dict(bin_id=x.index,log2_ratio=x.values,group=group,donor_alias=donor,mode=mode)))
  for tissue in ['cancer','adj_benign']:
   cells=z.index[z.role.eq('Myeloid')&z.tissue.eq(tissue)]
   if len(cells):
    x=mat[list(cells)].mean(axis=1)
    profiles.append(pd.DataFrame(dict(bin_id=x.index,log2_ratio=x.values,group='Myeloid | ALL_MYELOID | '+tissue,donor_alias=donor,mode=mode)))
  checks.append(dict(donor_alias=donor,mode=mode,status='DONE',n_cells=len(sc)))
status=pd.DataFrame(checks);status.to_csv(r/'private/completion_validation.tsv',sep='\t',index=False)
if not status.status.eq('DONE').all():raise RuntimeError('Incomplete jobs; no full-cohort conclusion')
cs=pd.DataFrame(cell_sens);assert len(cs)==int(meta.role.eq('Myeloid').sum())
cs.to_csv(r/'private/myeloid_cell_reference_sensitivity.tsv.gz',sep='\t',index=False,na_rep='NA')
cs.groupby(['subtype','tissue']).agg(n_cells=('cell','size'),n_evaluable_rho=('rho','count'),median_cell_profile_rho=('rho','median'),q10_cell_profile_rho=('rho',lambda x:x.quantile(.1)),median_absolute_log2_difference=('mean_absolute_log2_difference','median')).reset_index().to_csv(pub/'myeloid_cell_reference_sensitivity_summary.tsv',sep='\t',index=False,na_rep='NA')
sc=pd.concat(scores,ignore_index=True);pr=pd.concat(profiles,ignore_index=True)
assert set(sc.cell)==set(meta.index) and sc.groupby('cell').size().eq(2).all()
assert len(status)==2*len(plan)
sc.to_csv(r/'private/all_cell_scores.tsv.gz',sep='\t',index=False);pr.to_csv(r/'private/donor_group_bin_profiles.tsv.gz',sep='\t',index=False)
ds=sc.groupby(['donor_alias','mode','group']).agg(n_cells=('cell','size'),mean_abs=('mean_absolute_ratio_deviation','mean'),rms=('rms_ratio_deviation','mean')).reset_index()
all_sc=sc[sc.group.str.startswith('Myeloid | ')].copy();all_sc['group']='Myeloid | ALL_MYELOID | '+all_sc.group.str.split(' | ',regex=False).str[-1]
all_ds=all_sc.groupby(['donor_alias','mode','group']).agg(n_cells=('cell','size'),mean_abs=('mean_absolute_ratio_deviation','mean'),rms=('rms_ratio_deviation','mean')).reset_index();ds=pd.concat([ds,all_ds],ignore_index=True)
summary=ds.groupby(['mode','group']).agg(n_donors=('donor_alias','nunique'),n_cells=('n_cells','sum'),donor_mean_abs_median=('mean_abs','median'),donor_mean_abs_q25=('mean_abs',lambda x:x.quantile(.25)),donor_mean_abs_q75=('mean_abs',lambda x:x.quantile(.75)),donor_rms_median=('rms','median')).reset_index()
summary.to_csv(pub/'group_score_summary.tsv',sep='\t',index=False)
burden=pr.groupby(['donor_alias','mode','group']).agg(mean_abs_bin_log2=('log2_ratio',lambda x:x.abs().mean()),max_abs_bin_log2=('log2_ratio',lambda x:x.abs().max())).reset_index()
burden.groupby(['mode','group']).agg(n_donors=('donor_alias','nunique'),median_mean_abs_bin_log2=('mean_abs_bin_log2','median'),q25_mean_abs_bin_log2=('mean_abs_bin_log2',lambda x:x.quantile(.25)),q75_mean_abs_bin_log2=('mean_abs_bin_log2',lambda x:x.quantile(.75)),median_max_abs_bin_log2=('max_abs_bin_log2','median')).reset_index().to_csv(pub/'group_profile_burden_summary.tsv',sep='\t',index=False,na_rep='NA')
agg=pr.groupby(['mode','group','bin_id']).agg(mean_log2_ratio=('log2_ratio','mean'),n_donors=('log2_ratio','count')).reset_index()
agg.to_csv(pub/'cohort_bin_profiles.tsv',sep='\t',index=False,na_rep='NA')
paired_scores=[];paired_profiles=[]
for subtype in ['ALL_MYELOID']+sorted(meta.loc[meta.role.eq('Myeloid'),'subtype'].unique()):
 tumor='Myeloid | '+subtype+' | cancer';adj='Myeloid | '+subtype+' | adj_benign'
 for mode in ['primary','reference_split']:
  w=ds[ds['mode'].eq(mode)&ds.group.isin([tumor,adj])].pivot(index='donor_alias',columns='group',values='mean_abs')
  if tumor not in w or adj not in w:continue
  w=w[[tumor,adj]].dropna();delta=w[tumor]-w[adj]
  paired_scores.append(dict(subtype=subtype,mode=mode,n_pairs=len(w),n_tumor_higher=int((delta>0).sum()),n_tumor_lower=int((delta<0).sum()),median_tumor_minus_adjacent=float(delta.median())))
  z=pr[pr['mode'].eq(mode)&pr.group.isin([tumor,adj])&pr.donor_alias.isin(w.index)].pivot(index=['donor_alias','bin_id'],columns='group',values='log2_ratio')
  if tumor in z and adj in z:
   z=z[[tumor,adj]].dropna();x=(z[tumor]-z[adj]).rename('difference').reset_index();x['subtype']=subtype;x['mode']=mode;paired_profiles.append(x)
pd.DataFrame(paired_scores).to_csv(pub/'paired_tumor_adjacent_score_summary.tsv',sep='\t',index=False,na_rep='NA')
paired=pd.concat(paired_profiles,ignore_index=True)
paired.to_csv(r/'private/paired_donor_profile_differences.tsv.gz',sep='\t',index=False,na_rep='NA')
paired.groupby(['mode','subtype','bin_id']).agg(n_pairs=('difference','count'),mean_tumor_minus_adjacent_log2=('difference','mean')).reset_index().to_csv(pub/'paired_tumor_adjacent_bin_summary.tsv',sep='\t',index=False,na_rep='NA')
sens=[]
for (donor,group),z in pr.groupby(['donor_alias','group']):
 # T reference/holdout swap membership between runs; they are not matched groups.
 if group in ['T_reference','T_holdout']:continue
 wide=z.pivot(index='bin_id',columns='mode',values='log2_ratio').dropna()
 if len(wide)<3:continue
 x=wide.primary;y=wide.reference_split
 rho=float(spearmanr(x,y).statistic) if x.nunique()>1 and y.nunique()>1 else np.nan
 sens.append(dict(donor_alias=donor,group=group,n_bins=len(wide),rho=rho,mean_absolute_difference=float(np.mean(abs(x-y)))))
sens=pd.DataFrame(sens);sens.to_csv(r/'private/reference_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
ss=sens.groupby('group').agg(n_donors=('donor_alias','nunique'),n_evaluable_rho=('rho','count'),median_profile_rho=('rho','median'),min_profile_rho=('rho','min'),median_absolute_log2_difference=('mean_absolute_difference','median')).reset_index()
ss.to_csv(pub/'reference_sensitivity_summary.tsv',sep='\t',index=False,na_rep='NA')
bi=pd.concat(bin_info).drop_duplicates('bin_id');bi['chr_number']=bi.chr.str.replace('chr','',regex=False).astype(int);bi=bi.sort_values(['chr_number','start']);bin_order=bi.bin_id.tolist()
groups=sorted(agg.group.unique());groups=[g for g in ['T_reference','T_holdout','B_control','Epithelial_control'] if g in groups]+[g for g in groups if g.startswith('Myeloid')]
fig,axes=plt.subplots(2,1,figsize=(18, max(9,len(groups)*.55)),sharex=True)
vmax=.25
for ax,mode in zip(axes,['primary','reference_split']):
 m=agg[agg['mode'].eq(mode)].pivot(index='group',columns='bin_id',values='mean_log2_ratio').reindex(index=groups,columns=bin_order)
 im=ax.imshow(np.ma.masked_invalid(m.to_numpy()),aspect='auto',cmap=heatmap_cmap,vmin=-vmax,vmax=vmax,interpolation='nearest')
 nd=summary[summary['mode'].eq(mode)].set_index('group').n_donors
 ax.set_yticks(range(len(groups)),[g.replace('Myeloid | ','').replace('_control',' control')+' (n='+str(int(nd[g]))+')' for g in groups],fontsize=8)
 ax.set_title(mode+' reference: donor-equal mean inferred log2 ratio',fontsize=11)
 for c in range(2,23):
  ids=np.where(bi.chr_number.to_numpy()==c)[0]
  if len(ids):ax.axvline(ids[0]-.5,color='black',lw=.35)
 ticks=[np.where(bi.chr_number.to_numpy()==c)[0].mean() for c in range(1,23)]
 axes[-1].set_xticks(ticks,[str(c) for c in range(1,23)])
axes[-1].set_xlabel('Autosome; 10 Mb bins (>=10 retained genes/donor); per-bin coverage varies; gray = unavailable')
fig.suptitle('PRAD myeloid inferCNV | RNA-inferred signal, not DNA CNV or malignancy labels',fontsize=13)
fig.colorbar(im,cax=fig.add_axes([.92,.15,.012,.7]),label='Mean log2 inferred ratio; display clipped at +/-0.25')
fig.subplots_adjust(left=.27,right=.89,top=.94,bottom=.05,hspace=.17)
fig.savefig(pub/'cohort_infercnv_heatmap.png',dpi=180);fig.savefig(pub/'cohort_infercnv_heatmap.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(10,max(5,len(groups)*.32)))
for offset,mode,color in [(-.12,'primary','#2166ac'),(.12,'reference_split','#b2182b')]:
 z=summary[summary['mode'].eq(mode)].set_index('group').reindex(groups)
 y=np.arange(len(groups))+offset;x=z.donor_mean_abs_median
 ax.errorbar(x,y,xerr=np.array([x-z.donor_mean_abs_q25,z.donor_mean_abs_q75-x]),fmt='o',markersize=4,capsize=2,label=mode,color=color)
ax.set_yticks(range(len(groups)),[g.replace('Myeloid | ','') for g in groups],fontsize=8);ax.invert_yaxis();ax.set_xlabel('Mean absolute deviation from inferred ratio 1\nMedian and IQR across donors; no malignancy threshold');ax.legend();fig.tight_layout();fig.savefig(pub/'inferred_deviation_scores.png',dpi=180);plt.close(fig)
pa=pd.read_csv(pub/'paired_tumor_adjacent_bin_summary.tsv',sep='\t')
types=sorted(pa.subtype.unique());fig,axes=plt.subplots(2,1,figsize=(15,6.5),sharex=True)
for ax,mode in zip(axes,['primary','reference_split']):
 m=pa[pa['mode'].eq(mode)].pivot(index='subtype',columns='bin_id',values='mean_tumor_minus_adjacent_log2').reindex(index=types,columns=bin_order)
 im=ax.imshow(np.ma.masked_invalid(m.to_numpy()),aspect='auto',cmap=heatmap_cmap,vmin=-.1,vmax=.1,interpolation='nearest')
 counts={x['subtype']:x['n_pairs'] for x in paired_scores if x['mode']==mode}
 ax.set_yticks(range(len(types)),[t+' (pairs='+str(counts.get(t,0))+')' for t in types],fontsize=8);ax.set_title(mode)
 for c in range(2,23):
  ids=np.where(bi.chr_number.to_numpy()==c)[0]
  if len(ids):ax.axvline(ids[0]-.5,color='black',lw=.35)
axes[-1].set_xticks(ticks,[str(c) for c in range(1,23)]);axes[-1].set_xlabel('Autosome; 10 Mb bins; per-bin coverage varies; gray = unavailable')
fig.suptitle('Within-donor, same-subtype tumor minus adjacent inferred log2 ratio\nDescriptive differences; not DNA CNV calls',fontsize=12)
fig.colorbar(im,cax=fig.add_axes([.92,.15,.012,.65]),label='Difference; display clipped at +/-0.1')
fig.subplots_adjust(left=.18,right=.89,top=.86,bottom=.09,hspace=.35);fig.savefig(pub/'paired_tumor_adjacent_profiles.png',dpi=180);plt.close(fig)
report=dict(status='DONE',donors=int(len(plan)),runs=int(len(status)),myeloid_cells=int(meta.role.eq('Myeloid').sum()),total_cells=int(len(meta)),genotype_matching='NOT_RUN: no genotype data',malignancy_classification='NOT_RUN: RNA inferred profiles are not sufficient',HMM_calls='NOT_RUN: HMM=False',bin_min_genes=10,bin_width_bp=10000000,score='mean absolute deviation from denoised inferCNV ratio 1; donor summaries',reference_sensitivity='Two disjoint same-donor T pools; same study, not independent validation',public_data_scope='cohort/group summaries only; individual cell and donor measurements remain server165')
runqc=[dict(x.split(' ',1) for x in f.read_text().splitlines()) for f in (r/'private').glob('D*/*/DONE')]
report.update(cells_present_in_both_runs=True,retained_genes_min=min(int(x['genes']) for x in runqc),retained_genes_max=max(int(x['genes']) for x in runqc),inferred_ratio_min=min(float(x['min']) for x in runqc),inferred_ratio_max=max(float(x['max']) for x in runqc))
(pub/'completion_validation.json').write_text(json.dumps(report,indent=2))
print(summary.to_string(index=False));print(ss.to_string(index=False))
