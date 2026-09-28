"""Independent reaggregation from cached cell vectors, never exports patient values."""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

out = Path(sys.argv[1]); public = out / 'public'
m = pd.read_csv(out / 'private_cell_metadata.tsv', sep='\t', keep_default_na=False)
with np.load(out / 'private_cell_counts.npz') as archive:
    a = {k: archive[k] for k in archive.files if (archive[k] >= 0).all()}
checks = 0
for path in sorted((public / 'by_gene').glob('*.tsv')):
    g = path.stem
    table = pd.read_csv(path, sep='\t', keep_default_na=False)
    assert set(table.p_value) == {'NA'} and set(table.q_value) == {'NA'}
    for row in table.itertuples():
        col = row.annotation_level
        assert col in ['lineage','midway','subtype']
        selected = np.flatnonzero((m.tissue == row.tissue) & (m[col] == row.cell_type))
        groups = [selected[m.patient.to_numpy()[selected] == p] for p in sorted(set(m.patient.to_numpy()[selected]))]
        if row.support_set == 'patients_with_at_least_20_cells':
            groups = [v for v in groups if len(v) >= 20]
        assert len(groups) == int(row.n)
        assert sum(map(len, groups)) == int(row.n_cells)
        if len(groups) < 3 or g not in a:
            assert row.pseudobulk_CPM_median == 'NA' and row.detection_fraction_median == 'NA'
        else:
            cpm = [1e6 * sum(a[g][v].tolist()) / sum(a['total'][v].tolist()) for v in groups]
            detect = [np.count_nonzero(a[g][v]) / len(v) for v in groups]
            for metric,values in [('pseudobulk_CPM',cpm),('detection_fraction',detect)]:
                for stat,q in [('min',0),('q25',.25),('median',.5),('q75',.75),('max',1)]:
                    assert np.isclose(float(getattr(row,metric+'_'+stat)),np.quantile(values,q),rtol=1e-12,atol=1e-12)
        checks += 1
assert len(list((public / 'by_gene').glob('*.tsv'))) == 35
result = dict(status='PASS',rows_independently_reaggregated=checks,cell_vectors_and_group_masks_used=True,patient_values_exported=False)
(public / 'independent_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
