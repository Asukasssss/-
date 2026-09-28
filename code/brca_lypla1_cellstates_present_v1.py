"""Create public aggregate figures from private patient summaries on server165."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'pdf.fonttype':42,'font.size':10})
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A');R=B/'20260928T131206Z_lypla1_cellstates_v1';P=R/'public';D=R/'private'
def save(x,p):x.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.8g')
def main():
 s=pd.read_csv(P/'state_summary.tsv',sep='\t');extra=[]
 for co in ['Wu2021','Pal2021_reprocessed']:
  ps=pd.read_csv(D/(co+'_patient_states.tsv'),sep='\t');base=ps.assign(total=lambda x:x.mean_log*x.n_cells).groupby('patient').agg(total=('total','sum'),n=('n_cells','sum'));base['mean']=base.total/base.n;ps['centered']=ps.mean_log-ps.patient.map(base['mean']);ps.loc[ps.n_cells<20,['mean_log','detection','centered']]=np.nan
  for st,g in ps.groupby('state'):
   vv=g.centered.dropna();extra.append(dict(cohort=co,state=st,n_patients_centered=len(vv),fraction_above_own_malignant_mean=(vv>0).mean() if len(vv)>=3 else np.nan))
  order=s[s.cohort.eq(co)].sort_values('mean_log',ascending=False,na_position='last').state.tolist();ss=s[s.cohort.eq(co)].set_index('state').loc[order];fig,ax=plt.subplots(figsize=(max(8,len(order)*.6),3.6));ok=ss.status.eq('DONE').to_numpy();xx=np.arange(len(order));sc=ax.scatter(xx[ok],np.zeros(ok.sum()),s=50+450*ss.detection.to_numpy()[ok],c=ss.mean_log.to_numpy()[ok],cmap='viridis',vmin=0,edgecolors='#666666');ax.scatter(xx[~ok],np.zeros((~ok).sum()),marker='x',color='gray');ax.set_xticks(xx,[z.replace('Cancer ','')+'\nn='+str(int(ss.loc[z,'n_patients_eligible'])) for z in order],rotation=40,ha='right');ax.set_yticks([]);ax.set_ylim(-.5,.5);ax.set_title(co+' | LYPLA1 in existing malignant states');fig.colorbar(sc,ax=ax,label='Equal-patient mean log1p(CP10k)',fraction=.025);ax.spines[['top','right','left']].set_visible(False)
  for f in [.1,.5,.9]:ax.scatter([],[],s=50+450*f,color='#999999',label=str(int(f*100))+'%')
  ax.legend(title='Detection',loc='upper left',bbox_to_anchor=(1.15,1),frameon=False);fig.tight_layout();fig.savefig(P/(co+'_state_dotplot.png'),dpi=170,bbox_inches='tight');fig.savefig(P/(co+'_state_dotplot.pdf'),bbox_inches='tight');plt.close(fig)
  pat=base.sort_values('mean',ascending=False).index;mat=ps.pivot(index='patient',columns='state',values='centered').reindex(index=pat,columns=order);fig,ax=plt.subplots(figsize=(max(8,len(order)*.55),max(5,len(pat)*.22)));cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#eeeeee');lim=max(.1,np.nanmax(abs(mat.to_numpy())));im=ax.imshow(mat,aspect='auto',cmap=cmap,vmin=-lim,vmax=lim);ax.set_xticks(range(len(order)),[z.replace('Cancer ','') for z in order],rotation=45,ha='right');ax.set_yticks(range(len(pat)),['P%02d'%(i+1) for i in range(len(pat))]);ax.set_title(co+' | within-patient LYPLA1 state differences');ax.set_ylabel('Anonymous patient labels; sorted by overall expression');fig.colorbar(im,ax=ax,label='State mean minus own all-malignant mean');fig.text(.5,.01,'Gray: <20 cells / absent. Independent anonymous labels across cohorts. Existing annotations, no reclustering.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.03,1,1]);fig.savefig(P/(co+'_patient_state_heatmap.png'),dpi=170);fig.savefig(P/(co+'_patient_state_heatmap.pdf'));plt.close(fig)
 ss=s.merge(pd.DataFrame(extra),on=['cohort','state'],validate='one_to_one');save(ss,P/'state_summary.tsv')
 c=pd.read_csv(P/'cross_cohort_companions.tsv',sep='\t');e=c[c.replicated_descriptive].head(20).copy();label='Fixed descriptive replication rule'
 if e.empty:e=c[c.same_direction&(c.n_patients_Wu>=5)&(c.n_patients_Pal>=5)].head(20).copy();label='Largest shared effects; none met fixed replication rule'
 cols=['median_r_Wu','median_r_Pal','state_adjusted_median_r_Wu','state_adjusted_median_r_Pal'];fig,ax=plt.subplots(figsize=(7,max(4,len(e)*.28)));z=e[cols].to_numpy();lim=max(.1,np.nanmax(abs(z)));im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim);ax.set_yticks(range(len(e)),e.gene);ax.set_xticks(range(4),['Wu: depth','Pal: depth','Wu: depth+state','Pal: depth+state'],rotation=25,ha='right');ax.set_title('LYPLA1 within-patient companions\n'+label,fontsize=11)
 for i in range(len(e)):
  for j in range(4):ax.text(j,i,'%.2f'%z[i,j],ha='center',va='center',fontsize=8,color='white' if abs(z[i,j])>.65*lim else 'black')
 fig.colorbar(im,ax=ax,label='Equal-patient median partial Pearson r');fig.tight_layout();fig.savefig(P/'cross_cohort_companion_heatmap.png',dpi=180);fig.savefig(P/'cross_cohort_companion_heatmap.pdf');plt.close(fig);save(e,P/'displayed_companions.tsv')
 cv=pd.read_csv(P/'patient_coverage_anonymous.tsv',sep='\t');save(cv,D/'patient_coverage_anonymous.tsv');cv.groupby(['cohort','status','reason']).agg(n_patient_labels=('n_cells','size'),total_cells=('n_cells','sum')).reset_index().to_csv(P/'patient_coverage_summary.tsv',sep='\t',index=False);(P/'patient_coverage_anonymous.tsv').unlink()
 print(ss[['cohort','state','n_patients_eligible','mean_log','fraction_above_own_malignant_mean']].to_string(index=False));print('CANDIDATES',len(c[c.replicated_descriptive]));print(e[['gene']+cols].head(12).to_string(index=False))
if __name__=='__main__':main()
