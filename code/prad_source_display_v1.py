"""Display-only labels for non-rankable source rows; never alter expression or P/q."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/PRAD/06_EXTERNAL/20260922T142000Z_discovery_v1'
def main():
 p=OUT/'sc_celltype_profiles.tsv';d=pd.read_csv(p,sep='\t');r=pd.read_csv(OUT/'sc_source_stability.tsv',sep='\t').query('partition == "cancer"').set_index('gene')
 before=d[['gene','partition','celltype','effect','p_value','q_value']].copy()
 # For an unmeasured gene, n cannot be presented as a measured zero. Cell totals are known from other genes.
 absent=d.reason.eq('gene_not_measured');totals=d.loc[~absent].groupby(['partition','celltype']).n_cells_total.max()
 for i,z in d[absent].iterrows():
  d.loc[i,'n_cells_total']=totals.loc[(z.partition,z.celltype)]
  d.loc[i,['n','n_source_labels_eligible','n_cells_eligible']]=np.nan
 d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
 assert before.equals(d[['gene','partition','celltype','effect','p_value','q_value']])
 c=d.query('partition == "cancer"').pivot(index='gene',columns='celltype',values='effect').sort_index().dropna(axis=1,how='all');values=c.to_numpy(float)
 labels=[g+(' *' if r.loc[g,'status']!='DONE' else '') for g in c.index]
 det=d.query('partition == "cancer"').pivot(index='gene',columns='celltype',values='mean_detection_fraction').reindex(index=c.index,columns=c.columns)
 for kind in ['heatmap','dotplot']:
  fig,ax=plt.subplots(figsize=(12,24))
  if kind=='heatmap':
   z=np.full_like(values,np.nan);ok=np.isfinite(values).any(1);mu=np.nanmean(values[ok],axis=1,keepdims=True);sd=np.nanstd(values[ok],axis=1,keepdims=True);z[ok]=(values[ok]-mu)/np.where(sd>0,sd,1)
   cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d9d9d9');artist=ax.imshow(z,aspect='auto',cmap=cmap,vmin=-2,vmax=2,interpolation='nearest');bar='Row z-score of donor-equal mean'
  else:
   xx,yy=np.meshgrid(np.arange(len(c.columns)),np.arange(len(c)));ok=np.isfinite(values);artist=ax.scatter(xx[ok],yy[ok],s=det.to_numpy()[ok]*90,c=values[ok],cmap='viridis',vmin=0);ax.invert_yaxis();bar='Donor-equal mean log1p(CP10K)'
  ax.set_yticks(np.arange(len(c)));ax.set_yticklabels(labels,fontsize=7);ax.set_xticks(np.arange(len(c.columns)));ax.set_xticklabels(c.columns,rotation=45,ha='right',fontsize=9)
  ax.set_title('PRAD cancer tissue: all 101 candidate genes\n* Source ranking not evaluable (low signal or unmeasured); grey = NA',fontsize=11);fig.colorbar(artist,ax=ax,label=bar)
  if kind=='dotplot':ax.set_xlabel('Dot area = donor-equal raw-count detection fraction')
  fig.tight_layout();fig.savefig(OUT/('all101_source_'+kind+'.png'),dpi=180);fig.savefig(OUT/('all101_source_'+kind+'.pdf'));plt.close(fig)
 (OUT/'display_validation.json').write_text(json.dumps(dict(status='DONE',numeric_expression_and_P_q_unchanged=True,unmeasured_gene_n_now_NA=True,source_rank_none_rows_marked=True,new_statistics=False),indent=2),encoding='utf-8')
 rows=[dict(file=p.relative_to(OUT).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='checksums.tsv'];pd.DataFrame(rows).to_csv(OUT/'checksums.tsv',sep='\t',index=False,lineterminator='\n')
if __name__=='__main__':main()
