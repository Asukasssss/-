"""Keep215 historical columns; append new external layers and generate aggregate figures."""
from pathlib import Path
import csv,json,hashlib
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];RUN='20260921T072503Z_oslo_tang_external_v1';OUT=ROOT/'results/BRCA/06_EXTERNAL'/RUN
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:
  r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)
def write(p,c,r):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,c,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(r)
def main():
 a=pd.read_csv(OUT/'Tang_all174_associations.tsv',sep='\t');o=pd.read_csv(OUT/'Oslo_all174_evaluability.tsv',sep='\t');rna=pd.read_csv(OUT/'GSE42568_all117_RNA.tsv',sep='\t').set_index('gene')
 camp=pd.read_csv(ROOT/'reference/brca/association_summary_v1.tsv',sep='\t');camp=camp[camp.analysis_type.eq('author_processed_primary')].set_index('relation_id')
 fus=pd.read_csv(ROOT/'results/BRCA/06_EXTERNAL/20260920T090922Z_all174_external117_resources_v1/external_all174_with_CAMP_context.tsv',sep='\t');fus=fus[fus.test_family.eq('FUSCC_PRIMARY174')].set_index('relation_id');assert len(camp)==len(fus)==174
 compare=a.copy()
 for prefix,t,ef,q in [('CAMP_unadjusted',camp,'rho','q_value'),('FUSCC_primary174',fus,'effect','q_value')]:
  compare[prefix+'_rho']=compare.relation_id.map(t[ef]);compare[prefix+'_q']=compare.relation_id.map(t[q])
 compare['interpretation']='different_cohorts_and_processed_versions;no_formal_effect_difference_test;Tang_only20_cases'
 compare.to_csv(OUT/'all174_cross_cohort_context.tsv',sep='\t',index=False,na_rep='NA',lineterminator='\n')
 prior=ROOT/'results/BRCA/05_FUNCTION/20260921T034354Z_asns117_programs_v1/all117_comparison_asns_response_appended.tsv';cols,rows=read(prior);assert len(rows)==117
 added=['Tang_external_status','Tang_evaluable_relations','Tang_significant_relations','Tang_source','Oslo_status','GSE42568_RNA_status','GSE42568_log2_difference','GSE42568_RNA_q','new_external_scope']
 assert not set(cols)&set(added)
 for r in rows:
  g=r['gene'];t=a[a.gene.eq(g)];ok=t[t.status.eq('DONE')];sig=ok[ok.q_value.lt(.05)];rn=rna.loc[g]
  r.update(dict(zip(added,['DONE' if len(ok)==len(t) else 'PARTIAL' if len(ok) else 'NOT_EVALUABLE',str(len(ok)),';'.join(sig.relation_id),OUT.relative_to(ROOT).as_posix(),('NEEDS_REVIEW' if o[o.gene.eq(g)].status.eq('NEEDS_REVIEW').any() else 'NOT_EVALUABLE'),rn.status,str(rn.effect) if pd.notna(rn.effect) else 'NA',str(rn.q_value) if pd.notna(rn.q_value) else 'NA','Tang_patient_association;GSE42568_RNA_only;no_target_reranking'])))
 write(OUT/'all117_comparison_external_appended.tsv',cols+added,rows);_,original=read(prior);assert [{c:r[c] for c in cols} for r in rows]==original
 (OUT/'integration_validation.json').write_bytes((json.dumps(dict(all117_retained=True,historical_columns_preserved=len(cols),new_columns=len(added),old_values_exactly_preserved=True),indent=2)+'\n').encode())
 # Fixed displayed relationships include prior questions plus the two new signals.
 pairs=[('GPI','KEGG:C00668'),('GPI','KEGG:C05345'),('GPCPD1','KEGG:C00670'),('ASNS','KEGG:C00064'),('GLS','KEGG:C00064'),('SLC6A8','KEGG:C00300'),('MDH1','KEGG:C00149'),('PNP','KEGG:C00262')]
 labs=['GPI / G6P','GPI / F6P','GPCPD1 / GPC','ASNS / glutamine','GLS / glutamine','SLC6A8 / creatine','MDH1 / malate','PNP / hypoxanthine'];fig,ax=plt.subplots(figsize=(8,5))
 for i,(g,key) in enumerate(pairs):
  z=a[a.gene.eq(g)&a.metabolite_key.eq(key)].iloc[0];color='#a63603' if z.q_value<.05 else '#2166ac';ax.plot([z.ci_lower,z.ci_upper],[i,i],color=color);ax.plot(z.effect,i,'o',color=color);ax.text(1.02,i,f'q={z.q_value:.3g}',va='center',fontsize=9)
 ax.axvline(0,color='gray',lw=.8);ax.set_yticks(range(len(labs)));ax.set_yticklabels(labs);ax.invert_yaxis();ax.set_xlim(-1,1.35);ax.set_xlabel('Spearman rho; pointwise bootstrap 95% interval');ax.set_title('Tang / TCGA: 20 matched cases\nSelected display; FDR covers all163 evaluable original relationships');fig.tight_layout();fig.savefig(OUT/'Tang_selected_relationships.png',dpi=180);plt.close(fig)
 sf,sr=read(ROOT/'coordination/stages/BRCA.tsv');sr=[r for r in sr if r['run_id']!=RUN]
 for stage,status,scope in [('06_EXTERNAL','PARTIAL','Tang163 evaluable/174;GSE42568 all117 retained;Oslo sample crosswalk missing'),('07_INTEGRATION','DONE','All117 history retained;174 exact-relation cross-cohort context;no new significance gate')]:
  sr.append(dict(zip(sf,['BRCA',stage,RUN,'oslo_tang_external_v1',status,scope,OUT.relative_to(ROOT).as_posix(),'code/brca_finalize_new_external_v1.py','analysis/brca-functional-review-20260919','Tang20 cases;imputed metabolomics;RNA-only background;Oslo no inferred join','Resolve Oslo crosswalk;retain GPI nonconfirmation;do not replace full candidate pool by two Tang hits'])))
 write(ROOT/'coordination/stages/BRCA.tsv',sf,sorted(sr,key=lambda r:(r['stage_id'],r['run_id'])))
 scripts=['brca_oslo_tang_download_v1.py','brca_external_fetch_more_v1.py','brca_oslo_tang_analysis_v1.py','brca_gse42568_background_v1.py','brca_validate_new_external_v1.py','brca_finalize_new_external_v1.py']
 write(OUT/'code_manifest.tsv',['path','sha256'],[dict(path='code/'+n,sha256=hashlib.sha256((ROOT/'code'/n).read_bytes()).hexdigest()) for n in scripts])
 for p in OUT.iterdir():
  if p.suffix in {'.tsv','.md','.txt','.json'}:p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
 ch=OUT/'checksums.sha256';ch.write_bytes(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p!=ch).encode())
 print((OUT/'integration_validation.json').read_text())
if __name__=='__main__':main()
