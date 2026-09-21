"""Finalize only public aggregates; normalize inherited metadata labels."""
import json
import sys
import zipfile
from pathlib import Path
import pandas as pd

run=Path(sys.argv[1]); p=run/'public'
assert json.loads((p/'validation.json').read_text())['status']=='PASS'
assert json.loads((p/'independent_check.json').read_text())['status']=='PASS'
f=p/'coverage_by_lineage.tsv';t=pd.read_csv(f,sep='\t',keep_default_na=False)
t['run_id']=run.name;t['analysis_version']='COAD_Khaliq35_v1';t.to_csv(f,sep='\t',index=False)
allrows=pd.concat([pd.read_csv(f,sep='\t',keep_default_na=False) for f in sorted((p/'by_gene').glob('*.tsv'))],ignore_index=True)
primary=allrows[(allrows.tissue=='Tumor')&(allrows.annotation_level=='lineage')&(allrows.support_set=='patients_with_at_least_20_cells')].copy()
primary.to_csv(p/'tumor_lineage_all35.tsv',sep='\t',index=False)
overview=[]
for g,t in primary.groupby('gene',sort=True):
 eligible=t[(t.status=='DONE')&(t.pseudobulk_CPM_median!='NA')].copy()
 row=dict(gene=g,eligible_lineages=len(eligible),highest_median_CPM_lineage='NA',median_CPM='NA',median_detection='NA',n_qualified_patients='NA',status='NOT_EVALUABLE',reason='SOURCE_GENE_NOT_FOUND_OR_AMBIGUOUS;not zero expression')
 if len(eligible):
  eligible['value']=pd.to_numeric(eligible.pseudobulk_CPM_median)
  best=eligible.loc[eligible.value.idxmax()]; ties=eligible[eligible.value==best.value]
  row.update(highest_median_CPM_lineage=';'.join(ties.cell_type),median_CPM=best.value,median_detection=best.detection_fraction_median,n_qualified_patients=best.n,status='DONE',reason='Descriptive maximum among eligible author compartments;not specificity or function;ties retained')
  if best.value == 0:
   row.update(highest_median_CPM_lineage='NO_POSITIVE_MEDIAN',n_qualified_patients='NA',status='NOT_EVALUABLE',reason='ALL_ELIGIBLE_MEDIANS_ZERO;no localization by median;does not imply all cells zero')
  elif len(ties)>1:
   row.update(n_qualified_patients='NA',median_detection='NA',reason='Positive median tie;see compartment-specific coverage and detection')
 overview.append(row)
pd.DataFrame(overview).to_csv(p/'tumor_lineage_overview.tsv',sep='\t',index=False)
assert len(overview)==35 and len(primary)==210
allowed={'coverage_by_lineage.tsv','metadata_audit.json','gene_availability.tsv','compartment_checks.tsv','unassigned_coverage.tsv','validation.json','source_manifest.tsv','analysis_spec.json','independent_check.json','tumor_lineage_all35.tsv','tumor_lineage_overview.tsv'}
with zipfile.ZipFile(run/'public_delivery.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(p.rglob('*')):
  if f.is_file():
   assert f.name in allowed or (f.parent.name=='by_gene' and f.suffix=='.tsv')
   assert 'private' not in str(f.relative_to(p))
   z.write(f,str(f.relative_to(p)))
print(json.dumps(dict(genes=35,primary_rows=len(primary),all_summary_rows=len(allrows),public_files=len(allowed)+35)))
