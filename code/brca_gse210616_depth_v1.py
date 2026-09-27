"""Post-primary depth sensitivity: region indicator + log(total UMI), within section."""
import sys,json
from pathlib import Path
import pandas as pd,numpy as np
from scipy import stats
out=Path(sys.argv[1]);sec=pd.read_csv(out/'section_results.tsv',sep='\t');rows=[]
for r in sec.itertuples():
    d=pd.read_csv(out/f'{r.gsm}_spot_measurements.tsv.gz',sep='\t')
    d=d[d.marker_inferred_region.isin(['Epithelial','Stromal'])]
    counts=d.marker_inferred_region.value_counts()
    if min(counts.get('Epithelial',0),counts.get('Stromal',0))<20:continue
    a=(d.marker_inferred_region=='Epithelial').astype(float).to_numpy()
    X=np.column_stack([np.ones(len(d)),a,np.log1p(d.total_count.to_numpy())])
    fit=np.linalg.lstsq(X,d.LYPLA1_log1p.to_numpy(),rcond=None)[0]
    rows.append({'gsm':r.gsm,'patient':r.patient,'depth_adjusted_region_beta':fit[1]})
s=pd.DataFrame(rows);s.to_csv(out/'depth_section_results.tsv',sep='\t',index=False)
p=s.groupby('patient').depth_adjusted_region_beta.mean();p.to_csv(out/'depth_patient_results.tsv',sep='\t')
n=int((p!=0).sum());k=int((p>0).sum())
summary={'measure':'depth_adjusted_region_beta','n':len(p),'positive':k,'negative':int((p<0).sum()),'mean':float(p.mean()),'median':float(p.median()),'p_value':stats.binomtest(k,n).pvalue,
         'model':'within-section LYPLA1 log1p normalized ~ epithelial-marker indicator + log1p(total UMI), restrict epithelial/stromal inferred spots; equal section mean per patient',
         'status':'DONE','timing':'follow-up sensitivity after primary results; does not establish cell-intrinsic upregulation'}
(out/'depth_sensitivity.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
