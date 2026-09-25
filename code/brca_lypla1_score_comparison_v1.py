"""Presentation-only comparison of two completed fixed gene-set analyses."""
from pathlib import Path
import sys,json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

old,new=map(Path,sys.argv[1:3])
a=pd.read_csv(old/'results.tsv',sep='\t');b=pd.read_csv(new/'results.tsv',sep='\t')
assert a.test_id.tolist()==b.test_id.tolist()
a.insert(0,'gene_set','Reactome_PC_remodelling_24');b.insert(0,'gene_set','KEGG_glycerophospholipid_92')
pd.concat([a,b]).to_csv(new/'both_gene_sets_results16.tsv',sep='\t',index=False,na_rep='NA',lineterminator='\n')
fa=pd.read_csv(old/'frozen_score_members.tsv',sep='\t');fb=pd.read_csv(new/'frozen_score_members.tsv',sep='\t')
A=set(fa.loc[fa.used_in_all3_cohorts,'gene']);B=set(fb.loc[fb.used_in_all3_cohorts,'gene'])
overlap=pd.DataFrame([dict(gene=g,in_old24=g in A,in_new92=g in B) for g in sorted(A|B)])
overlap.to_csv(new/'gene_set_overlap.tsv',sep='\t',index=False,lineterminator='\n')
(new/'gene_set_overlap.json').write_bytes(json.dumps(dict(old_n=len(A),new_n=len(B),intersection=len(A&B),union=len(A|B),Jaccard=len(A&B)/len(A|B),new_only=len(B-A),old_only=len(A-B)),indent=2).encode())
plt.rcParams.update({'pdf.fonttype':42,'font.size':9,'savefig.dpi':180,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(12,6.5),constrained_layout=True)
for d,offset,color,label in [(a,-.14,'#0072B2','Previous Reactome 24'),(b,.14,'#D55E00','New KEGG 92')]:
 for i,r in d.iterrows():
  ax.plot([r.ci_lower,r.ci_upper],[i+offset,i+offset],color=color,lw=1.4)
  ax.scatter(r.effect,i+offset,c=color,s=27,label=label if i==0 else None)
  ax.text(1.03,i+offset,f'P={r.p_value:.4g}',transform=ax.get_yaxis_transform(),va='center',color=color,fontsize=8)
ax.set_yticks(range(8),a.test_id);ax.invert_yaxis();ax.set_xlim(-1,1);ax.axvline(0,color='#999999',lw=.7)
ax.set_xlabel('Correlation and pointwise 95% interval; changes between gene sets are not a formal difference test')
ax.set_title('Both fixed gene sets retained | identical scoring and testing rules\nSame cohorts reused: this is a supplementary exploration, not independent replication',weight='bold')
ax.legend(loc='lower right',frameon=False)
fig.savefig(new/'two_gene_sets_comparison.png');fig.savefig(new/'two_gene_sets_comparison.pdf');plt.close(fig)
with pd.ExcelWriter(new/'LYPLA_KEGGscore_results.xlsx',engine='openpyxl',mode='a',if_sheet_exists='replace') as w:
 pd.concat([a,b]).to_excel(w,sheet_name='两基因集16项并排',index=False)
 overlap.to_excel(w,sheet_name='基因集交集',index=False)
print('Comparison complete:',len(A&B),'shared genes;',len(B-A),'new-only')
