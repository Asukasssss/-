"""Count per-feature paired directions on server; export aggregates only."""
import importlib.util,json,csv,hashlib,traceback
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PREVIOUS=ROOT/'results/collaborative/PDAC/B/20260921T121509Z_paired_metabolites_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    run=Path(__file__).resolve().parent
    if run.parent!=ROOT/'results/collaborative/PDAC/B':raise ValueError('Server-only isolated run')
    with (run/'.running').open('x') as f:f.write('PDAC direction counts')
    try:
        out=run/'public';out.mkdir(exist_ok=False)
        spec=importlib.util.spec_from_file_location('paired',PREVIOUS/'paired_metabolites_v1.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        for p,h in mod.EXPECTED.items():
            if sha(p)!=h:raise ValueError('Input hash mismatch')
        if sha(mod.PAIR_SOURCE)!='0f16963287abb1bf7bf73dd5a2738f3bdccfb0344b2db3db72598199432ac866':raise ValueError('Pair source changed')
        _,pairs=mod.pair_join(pd.read_csv(mod.MAPPING,dtype=str),pd.read_excel(mod.PAIR_SOURCE,dtype=str))
        assert len(pairs)==11
        t=pd.read_excel(mod.MATRIX,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);n=pd.read_excel(mod.MATRIX,sheet_name='metabo_imputed_filtered_Normal',index_col=0);raw=pd.read_excel(mod.MATRIX,sheet_name='data',index_col=0)
        with (PREVIOUS/'public/results.tsv').open() as f:stats=list(csv.DictReader(f,delimiter='\t'))
        rows=[]
        for r in stats:
            name=r['metabolite_name'];x=t.loc[name,pairs.tumor_metab_id].to_numpy(float);y=n.loc[name,pairs.normal_metab_id].to_numpy(float);mask=np.isfinite(x)&np.isfinite(y)
            if r['analysis_type']!='primary':mask &= np.isfinite(raw.loc[name,pairs.tumor_metab_id].to_numpy(float))&np.isfinite(raw.loc[name,pairs.normal_metab_id].to_numpy(float))
            d=(x-y)[mask];up=int((d>0).sum());down=int((d<0).sum());equal=int((d==0).sum());den=len(d)
            assert den==int(r['n']) and up+down+equal==den
            if r['status']=='DONE':assert np.isclose(d.mean(),float(r['effect']))
            rows.append({'cancer':'PDAC','analysis_type':r['analysis_type'],'feature_name':name,'n_pairs':den,'n_up':up,'n_down':down,'n_equal':equal,'up_fraction':up/den if den else None,'down_fraction':down/den if den else None,'majority':'UP' if up>den/2 else 'DOWN' if down>den/2 else 'NO_STRICT_MAJORITY','old42':r['original_camp_q']!='NA' and float(r['original_camp_q'])<.05,'original_camp_q':r['original_camp_q'],'paired_p':r['p_value'],'paired_q':r['q_value'],'paired_status':r['status'],'direction_basis':'Exact sign of Tumor minus Normal on author-processed values;ties retained in denominator'})
        mod.write(out/'direction_counts.tsv',rows);mod.dump(out/'direction_counts.json',rows)
        mod.dump(out/'validation.json',{'status':'PASS','rows':len(rows),'primary_features':sum(r['analysis_type']=='primary' for r in rows),'author_pairs':11,'counts_sum_to_n':True,'n_and_mean_match_previous':True,'new_p_q_calculated':False,'patient_rows_exported':False})
        mod.write(out/'source_manifest.tsv',[{'path':str(p),'sha256':sha(p)} for p in [mod.MAPPING,mod.MATRIX,mod.PAIR_SOURCE,PREVIOUS/'public/results.tsv',PREVIOUS/'paired_metabolites_v1.py',Path(__file__)]])
        mod.dump(run/'DONE.json',{'scope':'Paired direction counts only'});(run/'.running').unlink();print(json.dumps({'run':str(run),'rows':len(rows)}))
    except Exception:
        (run/'FAILED.txt').write_text(traceback.format_exc());raise
if __name__=='__main__':main()
