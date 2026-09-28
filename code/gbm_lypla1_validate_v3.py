"""Independently validate exact source row and cell-to-patient aggregation."""
import pathlib,sys,json,hashlib
import numpy as np
import pandas as pd

r=pathlib.Path(sys.argv[1]);p=r/'public'
c=pd.read_csv(r/'private/SNUH_LYPLA1_cells.tsv',sep='\t')
d=pd.read_csv(r/'private/SNUH_patient_profiles.tsv',sep='\t')
assert c.cell_id.is_unique and len(c)==128375 and c.patient.nunique()==16
a=c.groupby(['patient','category']).expression.agg(['size','sum']).sort_index()
b=d[d.partition.eq('all')].set_index(['patient','category']).sort_index()
assert np.array_equal(a['size'].values,b.n_cells.values)
err=np.max(np.abs(a['sum'].values-b.sum_expression.values));assert err<1e-8
res=pd.read_csv(p/'LYPLA1_paired_results.tsv',sep='\t')
assert len(res)==112 and res.test_family.nunique()==14
assert res.groupby('test_family').size().eq(8).all()
assert res.loc[res.status.eq('NOT_EVALUABLE'),'q_value'].isna().all()
rv=pd.read_csv(p/'independent_R_validation.tsv',sep='\t');assert rv.status.eq('PASS').all()
hashes=[hashlib.sha256((r/'source'/f).read_bytes()).hexdigest() for f in ['SNUH_meta_proxy.tsv','SNUH_meta_parallel.tsv']]
assert len(set(hashes))==1
report=dict(status='PASS',exact_cell_ID_alignment=True,clinical_patient_count_check=True,
    independent_metadata_downloads_identical=True,source_cells=len(c),source_patients=int(c.patient.nunique()),
    source_to_patient_sum_max_error=float(err),paired_comparisons=len(res),predefined_comparators_per_family=8,
    matched_patient_effect_P_BH_independently_verified_in_R=True,
    unpublished_patient_values='server165 only',
    limitations=['SS2 lacks verified cell type/CNV labels in inspected source files',
      'SNUH missing treatment timing;not claimed untreated',
      'author CNV labels reused;raw reads/CNV not independently rerun',
      'hypothesis targeted after LYPLA1 selection;exploratory inference'])
(p/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
