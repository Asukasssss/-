"""Server verification from private celltype pseudobulk; export aggregate checks only."""
from pathlib import Path
import json,sys
import numpy as np,pandas as pd
run=Path(sys.argv[1]);assert run.parent.name=='B';allrows=pd.read_csv(run/'public/cell_source_all.tsv',sep='\t')
checked=0
for cohort in allrows.cohort.unique():
 p=pd.read_csv(run/(cohort+'_private_unit_pseudobulk.tsv'),sep='\t')
 assert not p.duplicated(['scheme','unit','celltype','gene']).any()
 good=p.sum_counts.notna();expected=np.log1p(1e6*p.loc[good,'sum_counts']/p.loc[good,'library_umi'])
 assert np.allclose(expected,p.loc[good,'log1p_cpm'])
 for r in allrows[allrows.cohort==cohort].itertuples():
  x=p[(p.scheme==r.annotation_level)&(p.gene==r.gene)&(p.celltype==r.celltype)]
  z=x[x.n_cells>=20];assert len(z)==r.n_units_ge20cells and x.n_cells.sum()==r.n_cells_total
  if np.isfinite(r.mean_unit_log1p_cpm):assert np.isclose(z.log1p_cpm.mean(),r.mean_unit_log1p_cpm)
  if np.isfinite(r.mean_unit_positive_fraction):assert np.isclose(z.positive_fraction.mean(),r.mean_unit_positive_fraction)
  assert (z.sum_counts>0).sum()==r.n_units_detected
  checked+=1
cons=pd.read_csv(run/'public/cross_cohort_source.tsv',sep='\t');assert len(cons)==173 and cons.gene.is_unique
for r in cons.itertuples():
 cols=[getattr(r,c) for c in ['GSE263733','GSE278688','GSE242230']]
 valid=[c for c in cols if c!='NOT_EVALUABLE' and ';' not in c];assert len(valid)==r.n_evaluable_cohorts
 assert max([valid.count(c) for c in set(valid)] or [0])==r.max_top_agreement
result={'status':'PASS','aggregate_rows_rechecked':checked,'genes_cross_cohort_checked':173,'scope':'Recompute per-unit normalization and summary arithmetic from private pseudobulk;not independent reannotation','unit_records_exported':False}
(run/'public/independent_summary_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
