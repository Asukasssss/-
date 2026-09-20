import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
def main(root):
    d=pd.read_csv(root/'public/external_all174_associations.tsv',sep='\t');z=pd.read_csv(root/'private/author_patient_joined_values.tsv',sep='\t');assert z.patient_id.is_unique
    assert len(d)==348 and not d.duplicated(['test_family','relation_id']).any()
    errs=[]
    for fam,a in d.groupby('test_family'):
        assert len(a)==174
        q=multipletests(a.p_value.fillna(1),method='fdr_bh')[1];ok=a.p_value.notna();assert np.allclose(q[ok],a.loc[ok,'q_value'])
        zz=z if fam=='FUSCC_PRIMARY174' else z[~z.patient_id.isin(['FUSCCTNBC030','FUSCCTNBC044','FUSCCTNBC140'])]
        for row in a.itertuples():
            if row.status!='DONE':assert pd.isna(row.p_value) and pd.isna(row.q_value);continue
            v=zz[[row.gene,row.external_peak]].dropna();assert len(v)==row.n
            errs.append(abs(stats.spearmanr(v.iloc[:,0],v.iloc[:,1]).statistic-row.effect))
    assert max(errs)<1e-12
    result={'status':'DONE','rows':len(d),'independent_rho_checks':len(errs),'max_rho_error':max(errs),'BH_statsmodels_agrees':True,'non_evaluable_public_P_q_missing':True,'unique_author_patient_ids':True,'no_patient_measurements_exported':True}
    (root/'public/independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
