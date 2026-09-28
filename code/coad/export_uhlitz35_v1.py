"""Export public aggregates only; all three annotation resolutions stay visible."""
import json
import sys
import zipfile
from pathlib import Path
import pandas as pd

run=Path(sys.argv[1]); p=run/'public'
assert json.loads((p/'validation.json').read_text())['status']=='PASS'
assert json.loads((p/'independent_check.json').read_text())['status']=='PASS'
allrows=pd.concat([pd.read_csv(f,sep='\t',keep_default_na=False) for f in sorted((p/'by_gene').glob('*.tsv'))],ignore_index=True)
primary=allrows[(allrows.tissue=='Tumor')&(allrows.support_set=='patients_with_at_least_20_cells')].copy()
primary.to_csv(p/'tumor_all35.tsv',sep='\t',index=False)
overview=[]
for (g,level),t in primary.groupby(['gene','annotation_level'],sort=True):
    eligible=t[(t.status=='DONE')&(t.pseudobulk_CPM_median!='NA')].copy()
    row=dict(gene=g,annotation_level=level,eligible_cell_types=len(eligible),highest_median_CPM_cell_type='NA',
             median_CPM='NA',median_detection='NA',n_qualified_patients='NA',status='NOT_EVALUABLE',
             reason='Missing source row or insufficient coverage;not measured zero')
    if len(eligible):
        eligible['value']=pd.to_numeric(eligible.pseudobulk_CPM_median)
        best=eligible.loc[eligible.value.idxmax()];ties=eligible[eligible.value==best.value]
        row.update(highest_median_CPM_cell_type=';'.join(ties.cell_type),median_CPM=best.value,
                   median_detection=best.detection_fraction_median,n_qualified_patients=best.n,status='DONE',
                   reason='Descriptive highest median among covered types;not specificity,contribution or function')
        if best.value==0:
            row.update(highest_median_CPM_cell_type='NO_POSITIVE_MEDIAN',n_qualified_patients='NA',status='NOT_EVALUABLE',
                       reason='All eligible medians zero;does not imply no expressing cells')
        elif len(ties)>1:
            row.update(n_qualified_patients='NA',median_detection='NA',reason='Positive ties retained;see type-specific distributions')
    overview.append(row)
pd.DataFrame(overview).to_csv(p/'tumor_overview.tsv',sep='\t',index=False)
assert len(overview)==35*3
allowed={'coverage_by_lineage.tsv','metadata_audit.json','gene_availability.tsv','matrix_checks.tsv',
         'extraction_validation.json','validation.json','source_manifest.tsv','analysis_spec.json',
         'independent_check.json','tumor_all35.tsv','tumor_overview.tsv'}
with zipfile.ZipFile(run/'public_delivery.zip','w',zipfile.ZIP_DEFLATED) as z:
    for f in sorted(p.rglob('*')):
        if f.is_file():
            assert f.name in allowed or (f.parent.name=='by_gene' and f.suffix=='.tsv')
            assert 'private' not in str(f.relative_to(p))
            z.write(f,str(f.relative_to(p)))
print(json.dumps(dict(genes=35,primary_rows=len(primary),all_summary_rows=len(allrows),public_only=True)))
