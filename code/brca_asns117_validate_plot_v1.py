"""Validate frozen-FDR joins and BH12; plot aggregate-only response summaries."""
from pathlib import Path
import argparse,json,hashlib
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save(x,p):x.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main(root):
 pub=root/'public';a=pd.read_csv(pub/'all117_four_contrast_response.tsv',sep='\t');o=pd.read_csv(pub/'orthology117.tsv',sep='\t');s=pd.read_csv(pub/'all117_two_construct_summary.tsv',sep='\t');assert len(a)==468 and a.human_gene.nunique()==117 and len(s)==234
 old=Path(json.loads((root/'analysis_spec.json').read_text())['previous_public'])
 errors=[]
 for site in ['Tumor','Lung']:
  for c in ['shAsns1','shAsns2']:
   z=pd.read_csv(old/(site+'_'+c+'_all_DE.tsv'),sep='\t').set_index('gene')
   for x in a[(a.site==site)&(a.construct==c)&(a.status=='DONE')].itertuples():
    ref=z.loc[x.mouse_gene];errors.extend([abs(x.logFC-ref.logFC),abs(x.original_q_global4-ref.q_global4),abs(x.p_value-ref.PValue)])
 assert max(errors)<1e-12
 b=pd.read_csv(pub/'program_results.tsv',sep='\t');assert len(b)==12
 pe=float(np.max(abs(multipletests(b.p_value,method='fdr_bh')[1]-b.q_value)));assert pe<1e-12
 # Shared set genes are documented, not assumed independent.
 sets=json.loads((root/'source/frozen_program_sets.json').read_text());names=list(sets)
 overlaps=[dict(program1=g,program2=h,overlap_without_Asns=len((set(sets[g])&set(sets[h]))-{'Asns'})) for i,g in enumerate(names) for h in names[i+1:]];save(pd.DataFrame(overlaps),pub/'program_overlap.tsv')
 cols=[('Tumor','shAsns1'),('Tumor','shAsns2'),('Lung','shAsns1'),('Lung','shAsns2')]
 programs=['GOBP_RESPONSE_TO_AMINO_ACID_STARVATION','GOBP_AMINO_ACID_TRANSMEMBRANE_TRANSPORT','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']
 v=np.array([[b[(b.program==g)&(b.site==t)&(b.construct==c)].effect.iloc[0] for t,c in cols] for g in programs]);q=np.array([[b[(b.program==g)&(b.site==t)&(b.construct==c)].q_value.iloc[0] for t,c in cols] for g in programs])
 fig,ax=plt.subplots(figsize=(8,3.8));im=ax.imshow(v,cmap='RdBu_r',vmin=-.35,vmax=.35,aspect='auto')
 for i in range(3):
  for j in range(4):ax.text(j,i,f'{v[i,j]:+.3f}'+('*' if q[i,j]<.05 else ''),ha='center',va='center')
 ax.set_yticks(range(3));ax.set_yticklabels(['Amino acid starvation response','Amino acid transmembrane transport','EMT-related expression']);ax.set_xticks(range(4));ax.set_xticklabels([t+'\n'+c for t,c in cols]);ax.set_title('Asns perturbation: three fixed transcriptional programs\n* cameraPR BH12 < 0.05; colors are median gene log2FC');fig.colorbar(im,ax=ax,label='Median log2FC (descriptive)');fig.tight_layout();fig.savefig(pub/'program_response.png',dpi=180);plt.close(fig)
 genes=sorted(a.human_gene.unique());mat=np.array([[a[(a.human_gene==g)&(a.site==t)&(a.construct==c)].logFC.iloc[0] for t,c in cols] for g in genes]);qs=np.array([[a[(a.human_gene==g)&(a.site==t)&(a.construct==c)].original_q_global4.iloc[0] for t,c in cols] for g in genes])
 fig,ax=plt.subplots(figsize=(7,23));cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#dddddd');im=ax.imshow(mat,cmap=cm,vmin=-2,vmax=2,aspect='auto')
 for i in range(len(genes)):
  for j in range(4):
   if qs[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=7)
 ax.set_yticks(range(len(genes)));ax.set_yticklabels(genes,fontsize=7);ax.set_xticks(range(4));ax.set_xticklabels([t+'\n'+c for t,c in cols]);ax.set_title('All117 human candidates: mapped mouse RNA response\nGray=unevaluable; *=original global49246 q<0.05\nNot self-gene functional validation',fontsize=10);fig.colorbar(im,ax=ax,fraction=.02,pad=.03,label='Mouse log2FC; color clipped at +/-2');fig.tight_layout();fig.savefig(pub/'all117_response_heatmap.png',dpi=160);plt.close(fig)
 counts=s.groupby(['site','response_state']).size().rename('genes').reset_index();save(counts,pub/'response_counts.tsv')
 selected=s[s.response_state.isin(['BOTH_SUPPORTED_UP','BOTH_SUPPORTED_DOWN'])];save(selected,pub/'both_construct_supported_candidates.tsv')
 val=dict(status='DONE',all117_retained=True,response_rows=468,no_new_candidate_DE_test=True,old_stat_max_abs_error=max(errors),program_BH12_max_abs_error=pe,interpretation='Asns downstream RNA response;not own gene target validation;programs not flux or phenotype')
 (pub/'validation.json').write_text(json.dumps(val,indent=2)+'\n')
 manifest=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in list((root/'source').glob('*.gmt'))+[root/'source/frozen_program_sets.json',root/'analysis_spec.json']+list(root.glob('*.py'))+list(root.glob('*.R'))];save(pd.DataFrame(manifest),pub/'program_source_manifest.tsv')
 print(counts.to_string(index=False));print(selected.to_string(index=False));print(o[o.status!='DONE'].to_string(index=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
