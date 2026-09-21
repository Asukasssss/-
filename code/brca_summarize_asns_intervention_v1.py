"""Validate aggregate DE, preserve resource decisions, and render node heatmap."""
from pathlib import Path
import argparse,json,hashlib
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests

def save(x,p):x.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main(root):
 pub=root/'public';files=sorted(pub.glob('*_all_DE.tsv'));assert len(files)==4
 a=pd.concat([pd.read_csv(p,sep='\t') for p in files],ignore_index=True)
 assert not a.duplicated(['site','construct','gene']).any()
 q=multipletests(a.PValue,method='fdr_bh')[1];err=float(np.max(np.abs(q-a.q_global4)));assert err<1e-10
 target=a[a.gene.eq('Asns')];assert len(target)==4 and (target.logFC<0).all()
 node=a[a.gene.isin(['Asns','Gls','Gls2','Glul'])].copy();assert len(node)==16
 node['fold_change']=2**node.logFC
 rows=[]
 for r in node.itertuples():
  rows.append(dict(cancer='BRCA',cohort='GSE104966_mouse4T1T_'+r.site,stage_id='05_FUNCTION',run_id=root.name,analysis_version='asns_intervention_v1',analysis_type='Asns_knockdown_RNA_edgeR_QL',metabolite_key='NA',metabolite_name='NA',gene=r.gene,unit='author_biological_replicate',n=4,n_reference=4,effect_type='log2_RNA_fold_change',effect=r.logFC,ci_lower=np.nan,ci_upper=np.nan,p_value=r.PValue,q_value=r.q_global4,test_family='all_tested_genes_four_contrasts',family_n_evaluable=len(a),status='DONE',reason='RNA_only;mouse_selected_recultured_cells;CI_not_estimated',source_id='GSE104966',construct=r.construct,site=r.site,fold_change=r.fold_change,q_within_contrast=r.q_within_contrast,perturbed_gene='Asns'))
 save(pd.DataFrame(rows),pub/'results.tsv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 cols=[('Tumor','shAsns1'),('Tumor','shAsns2'),('Lung','shAsns1'),('Lung','shAsns2')];genes=['Asns','Gls','Gls2','Glul']
 vals=np.array([[node[(node.gene==g)&(node.site==s)&(node.construct==c)].logFC.iloc[0] for s,c in cols] for g in genes])
 qs=np.array([[node[(node.gene==g)&(node.site==s)&(node.construct==c)].q_global4.iloc[0] for s,c in cols] for g in genes])
 fig,ax=plt.subplots(figsize=(7,4));im=ax.imshow(vals,cmap='RdBu_r',vmin=-2.5,vmax=2.5,aspect='auto')
 for i in range(4):
  for j in range(4):ax.text(j,i,f'{vals[i,j]:+.2f}'+('*' if qs[i,j]<.05 else ''),ha='center',va='center',color='white' if abs(vals[i,j])>1.5 else 'black')
 ax.set_xticks(range(4));ax.set_xticklabels([s+'\n'+c for s,c in cols]);ax.set_yticks(range(4));ax.set_yticklabels(genes);ax.set_title('Asns knockdown RNA response | GSE104966\n* Global BH across four transcriptome contrasts < 0.05');fig.colorbar(im,ax=ax,label='log2 fold change versus shRenilla');fig.tight_layout();fig.savefig(pub/'asns_node_RNA_heatmap.png',dpi=180);plt.close(fig)
 resources=[('ASNS','GSE104967','direct_Asns_KD_invitro','NEEDS_REVIEW','Second control column all zero;one effective control;no DE'),('ASNS','GSE104966','direct_Asns_KD_recovered_recultured','DONE','24 libraries;four author biological replicates per group;four contrasts;RNA only'),('ASNS','GSE104964','asparagine_supply_in_Asns_silenced_cells','NOT_RUN','2 vs2 author replicates;nutrient contrast;not a separate ASNS KD experiment'),('ASNS','GSE107109','patient_primary_vs_metastatic','NOT_EVALUABLE','not gene intervention'),('GLS','PMID28950000','genetic_GLS_KD_paper','NEEDS_REVIEW','Cached fulltext supplementary list only ARRIVE checklist;no processed perturbation matrix located in inspected attachments'),('GLS','GSE263696','chronic_glutamine_deprivation_adaptation','NOT_RUN','MDA-MB-231 adapted derivative versus parental;not GLS gene knockdown'),('GLS','GSE26370','glutamine_deprivation','NOT_RUN','nutrient context not GLS knockdown'),('GLS','GSE173991','mixed_metabolic_profiling_study','NEEDS_REVIEW','summary mentions drug efficacy;expression sample design not verified as GLS-targeted intervention')]
 save(pd.DataFrame(resources,columns=['gene','source_id','resource_type','status','reason']),pub/'resource_decisions.tsv')
 searches=[]
 for i in [0,1]:
  p=root/'source'/f'search{i}.json';j=json.loads(p.read_text())['esearchresult'];searches.append(dict(query_translation=j['querytranslation'],hits=j['count'],retrieved=len(j['idlist']),scope='bounded_gene_keyword_GEO_series_search;not_exhaustive'))
 save(pd.DataFrame(searches),pub/'search_log.tsv')
 # Whole-transcriptome outputs have native edgeR fields; results.tsv provides the common statistical prefix for the prespecified node.
 (pub/'independent_validation.json').write_text(json.dumps(dict(status='PASS',BH_global_max_abs_error=err,unique_gene_contrast_keys=True,node_rows=16,Asns_down_all_four=True,new_metabolite_measurements=0),indent=2)+'\n')
 old=pd.read_csv(pub/'source_manifest.tsv',sep='\t').to_dict('records')
 for f in sorted(root.glob('*.R'))+sorted(root.glob('*.py'))+[root/'analysis_spec.json',root/'source/author_column_assignment.tsv']+sorted((root/'source').glob('search*.json'))+sorted((root/'source').glob('summary*.json')):
  old.append(dict(path=str(f),url='NA',sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 save(pd.DataFrame(old).drop_duplicates('path'),pub/'source_manifest.tsv')
 print(node[['gene','site','construct','logFC','q_global4','fold_change']].to_string(index=False));print((pub/'construct_concordance.tsv').read_text());print((pub/'validation.json').read_text())
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
