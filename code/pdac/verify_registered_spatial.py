"""Independently replay aggregate calculations from private spot extracts on server."""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd
r=Path(sys.argv[1]);t=pd.read_csv(r/'public/section_summary.tsv',sep='\t');checked=0
for _,z in t.iterrows():
 d=pd.read_csv(r/'private'/(z['sample']+'_spots.tsv.gz'),sep='\t')
 assert d.barcode.is_unique
 assert np.allclose(d.log1p,np.log1p(d.raw*10000/d.total),rtol=1e-12,atol=1e-12)
 if z['subset']=='qc_sensitivity':d=d[(d.genes>=200)&(d.total>=500)&(d.mito_fraction<=.2)]
 assert len(d)==z['n']
 low=d[d.epithelial<=d.epithelial.quantile(.25)];high=d[d.epithelial>=d.epithelial.quantile(.75)]
 assert len(low)==z['low_n'] and len(high)==z['high_n']
 assert int((d.raw>0).sum())==z['positive']
 for field,actual in [('low_mean',low.log1p.mean()),('high_mean',high.log1p.mean()),('high_minus_low',high.log1p.mean()-low.log1p.mean()),('high_detection_pct',100*(high.raw>0).mean()),('low_detection_pct',100*(low.raw>0).mean())]:
  assert np.isclose(float(z[field]),actual,rtol=1e-10,atol=1e-10),(z['sample'],field)
 checked+=1
(r/'public/aggregate_replay_validation.json').write_text(json.dumps({'status':'PASS','section_subset_rows_checked':checked,'normalization_from_raw_rechecked':True,'quartile_groups_rechecked':True,'detection_means_and_differences_rechecked':True,'not_independent_biological_validation':True},indent=2))
print('PASS',checked)
