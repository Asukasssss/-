"""KEGG lipid panel and full-pathway GSEA presentation."""
from pathlib import Path
import sys,textwrap,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(sys.argv[1]);P=R/'public';assert (R/'NUMERICAL_DONE').exists();spec=json.loads((P/'analysis_spec.json').read_text());terms=spec['lipid_display_terms'];cos=['Wu2021','Pal2021_reprocessed'];ct='high_vs_low_detected'
plt.rcParams.update({'pdf.fonttype':42,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
g=pd.read_csv(P/'GSEA_all4.tsv',sep='\t');o=pd.read_csv(P/'ORA_all4.tsv',sep='\t');cover=pd.read_csv(P/'pathway_coverage.tsv',sep='\t')
pan=[]
for term in terms:
 for co in cos:
  a=g[g.Term.eq(term)&g.cohort.eq(co)&g.contrast.eq(ct)];b=o[o.pathway.eq(term)&o.cohort.eq(co)&o.contrast.eq(ct)]
  cv=cover[cover.pathway.eq(term)&cover.cohort.eq(co)&cover.contrast.eq(ct)];row=dict(pathway=term,cohort=co,status='DONE' if len(a) else 'NOT_EVALUABLE',reason='OK' if len(a) else ('not_in_fixed_library' if len(cv)==0 else 'outside_size15_to500'))
  if len(a):row.update(NES=a.NES.iloc[0],GSEA_P=a['NOM p-val'].iloc[0],GSEA_FDR=a['FDR q-val'].iloc[0],GSEA_leading_genes=a.Lead_genes.iloc[0])
  for direction in ['up','down']:
   z=b[b.direction.eq(direction)]
   if len(z):row.update({direction+'_ORA_P':z.p_value.iloc[0],direction+'_ORA_q':z.q_contrast.iloc[0],direction+'_genes':z.genes.iloc[0]})
  pan.append(row)
panel=pd.DataFrame(pan);panel.to_csv(P/'lipid_primary_comparison.tsv',sep='\t',index=False,na_rep='NA')
fig,axes=plt.subplots(1,2,figsize=(14,8),constrained_layout=True)
for ax,co in zip(axes,cos):
 d=panel[panel.cohort.eq(co)].set_index('pathway').reindex(terms);vals=d.NES.fillna(0).to_numpy();bars=ax.barh(np.arange(len(terms)),vals,color=np.where(d.GSEA_P<.05,np.where(vals>0,'#CC6677','#4477AA'),'#BBBBBB'));ax.set_yticks(np.arange(len(terms)),terms,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='black',lw=.5);ax.set_xlabel('GSEA NES; grey: nominal P>=0.05');ax.set_title(co)
 for i,(_,r) in enumerate(d.iterrows()):
  label='not evaluable' if r.status!='DONE' else ('P below permutation resolution' if r.GSEA_P==0 else f'P={r.GSEA_P:.3g}')

 # P values are in the adjacent downloadable table; unavailable pathways explicitly labelled.
 for i,st in enumerate(d.status):
  if st!='DONE':ax.text(0,i,' not evaluable',va='center',fontsize=7)
fig.suptitle('KEGG lipid-related pathways | LYPLA1 high vs low detected\nExploratory pooled-cell ranking; complete P/FDR and leading genes in table')
fig.savefig(P/'01_KEGG_lipid_GSEA.png',dpi=170);fig.savefig(P/'01_KEGG_lipid_GSEA.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(10,8),constrained_layout=True);cols=[];mat=[]
for co in cos:
 for direction in ['up','down']:
  cols.append(co+'\n'+direction);d=panel[panel.cohort.eq(co)].set_index('pathway').reindex(terms);mat.append(-np.log10(d[direction+'_ORA_P'].clip(lower=1e-10)))
M=np.array(mat).T;im=ax.imshow(np.ma.masked_invalid(M),aspect='auto',cmap='YlOrRd',vmin=0,vmax=max(3,float(np.nanmax(M))));ax.set_xticks(range(4),cols);ax.set_yticks(range(len(terms)),terms);fig.colorbar(im,ax=ax,label='-log10(ORA nominal P)')
for i in range(len(terms)):
 for j in range(4):
  if not np.isfinite(M[i,j]):label='NA'
  else:label=f'{10**(-M[i,j]):.2g}'+('*' if M[i,j]>-np.log10(.05) else '')
  ax.text(j,i,label,ha='center',va='center',fontsize=7,color='white' if np.isfinite(M[i,j]) and M[i,j]>max(3,float(np.nanmax(M)))*.6 else 'black')
ax.set_title('KEGG ORA | high vs low detected\nUp/down selected genes separately; * nominal P<0.05')
fig.savefig(P/'02_KEGG_lipid_ORA.png',dpi=170);fig.savefig(P/'02_KEGG_lipid_ORA.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(16,11),constrained_layout=True)
for ax,co in zip(axes,cos):
 d=g[g.cohort.eq(co)&g.contrast.eq(ct)].copy();d['absNES']=abs(d.NES);top=pd.concat([d[d.NES>0].sort_values(['FDR q-val','absNES'],ascending=[True,False]).head(8),d[d.NES<0].sort_values(['FDR q-val','absNES'],ascending=[True,False]).head(8)]).sort_values('NES');ax.barh(np.arange(len(top)),top.NES,color=np.where(top['NOM p-val']>=.05,'#BBBBBB',np.where(top.NES>0,'#CC6677','#4477AA')));ax.set_yticks(np.arange(len(top)),['\n'.join(textwrap.wrap(x,36)) for x in top.Term],fontsize=8);ax.axvline(0,color='black',lw=.6);ax.set_title(co);ax.set_xlabel('NES; grey: nominal P>=0.05')
fig.suptitle('KEGG whole-ranked-list GSEA | up to 8 terms per direction\nDisplay by empirical FDR then |NES|; enriched disease names are not disease diagnoses')
fig.savefig(P/'03_KEGG_allpathways_GSEA.png',dpi=160);fig.savefig(P/'03_KEGG_allpathways_GSEA.pdf');plt.close(fig)
print(panel[['pathway','cohort','NES','GSEA_P','up_ORA_P','down_ORA_P']].to_string(index=False))
