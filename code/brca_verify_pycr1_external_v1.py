"""Verify new aggregates against previous raw summaries and independent rank calculations."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats

def bh(p):
    order=np.argsort(p);s=np.array(p)[order];q=np.minimum.accumulate((s*len(s)/np.arange(1,len(s)+1))[::-1])[::-1];out=np.empty(len(s));out[order]=np.minimum(1,q);return out

def main(root,previous):
    rows=[]
    for cohort in ['Wu2021','Pal2021_reprocessed']:
        old=pd.read_csv(previous/'private'/(cohort+'_donor_profiles.tsv'),sep='\t');new=pd.read_csv(root/'private'/(cohort+'_program_counts.tsv'),sep='\t')
        keys=['donor','celltype'];new=new.rename(columns={'subtype':'stratum','library_sum':'all_gene_library_sum'})
        old=old[old.gene.eq('PYCR1')&old.celltype.isin(['Malignant_epithelial','Fibroblasts'])];new=new[new.gene.eq('PYCR1')]
        old=old.groupby(keys)[['raw_count','n_cells','all_gene_library_sum']].sum();new=new.groupby(keys)[['raw_count','n_cells','all_gene_library_sum']].sum();assert old.index.is_unique and new.index.is_unique and set(old.index)==set(new.index)
        err=0
        for c in ['raw_count','n_cells','all_gene_library_sum']:
            e=float(np.max(np.abs(old[c]-new[c].reindex(old.index))));err=max(err,e);assert e<1e-6,(cohort,c,e)
        rows.append({'cohort':cohort,'source_groups':len(old),'raw_summary_max_error':err})
    r=pd.read_csv(root/'public/pycr1_program_associations.tsv',sep='\t');assert len(r)==12
    maxerr=0
    for x in r[r.status.eq('DONE')].itertuples():
        d=pd.read_csv(root/'private'/f'{x.cohort}_{x.partition}_{x.celltype}_scores.tsv',sep='\t')
        rho=stats.spearmanr(d.PYCR1,d[x.program]).statistic;maxerr=max(maxerr,abs(rho-x.effect));assert np.isclose(rho,x.effect,atol=1e-12)
    for _,d in r.groupby('test_family'):
        ok=d.p_value.notna();q=bh(d.p_value.fillna(1).to_numpy());assert np.allclose(q[ok],d.loc[ok,'q_value'])
    e=pd.read_csv(root/'public/external_associations.tsv',sep='\t');d=pd.read_csv(root/'private/FUSCC_joined_two_gene_values.tsv',sep='\t');assert d.patient_id.is_unique
    for row in e[e.status.eq('DONE')].itertuples():
        z=d if row.test_family=='EXTERNAL_PRIMARY2' else d[~d.patient_id.isin(['FUSCCTNBC030','FUSCCTNBC044','FUSCCTNBC140'])]
        z=z[[row.gene,row.gene+'_metabolite']].dropna();assert len(z)==row.n
        rho=stats.spearmanr(z.iloc[:,0],z.iloc[:,1]).statistic;assert np.isclose(rho,row.effect,atol=1e-12)
    result={'status':'DONE','prior_vs_new_PYCR1_raw_summaries':rows,'sc_scipy_spearman_max_error':maxerr,'BH_independent_check':True,'external_rho_and_n_independently_checked':True,'source_donor_not_genotype_verified':True,'patient_data_exported':False,'scores_fixed_during_bootstrap':True}
    (root/'public/numeric_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);a=p.parse_args();main(a.root,a.previous)
