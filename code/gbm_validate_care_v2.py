"""Independent reconstruction of patient/category coverage and published means."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd

def main(out):
    pub=out/'public/06_EXTERNAL';m=pd.read_csv(out/'private/CARE_author_cell_metadata.tsv',sep='\t');s=pd.read_csv(out/'private/CARE_sample_metadata.tsv',sep='\t');s.columns=s.columns.str.strip()
    for c in s.select_dtypes('object'):s[c]=s[c].str.strip()
    m=m.merge(s,left_on='ID',right_on='Sample ID',validate='many_to_one');assert len(m)==429305
    m['partition']=np.where(m['Primary or recurrent'].eq('Primary'),'CARE2025_primary','CARE2025_recurrent')
    exclude=m['Primary or recurrent'].eq('Primary')&~(m['Radiation before the surgery'].eq('No')&m['Alkylate agents before the surgery'].eq('No'))
    assert exclude.sum()==795;m=m[~exclude]
    mp=pd.read_csv(pub/'sc_annotation_map.tsv',sep='\t').set_index('original').coarse.to_dict();m['celltype']=m.CellType.map(mp);assert m.celltype.notna().all()
    actual=m.groupby(['partition','Patient ID','celltype']).size()
    d=pd.read_csv(out/'private/sc_donor_profiles_private.tsv',sep='\t')
    for key,z in d.groupby(['partition','Patient ID','celltype']):
        assert len(z)==142 and z.gene.nunique()==142
        assert z.n_cells.eq(actual.loc[key]).all()
    prof=pd.read_csv(pub/'sc_celltype_profiles.tsv',sep='\t');errs=[];evaluated=0
    for row in prof.itertuples():
        a=d[d.partition.eq(row.cohort)&d.celltype.eq(row.celltype)&d.gene.eq(row.gene)];a=a[a.n_cells.ge(20)&a.mean_expression.notna()]
        assert len(a)==row.n
        if row.status=='DONE':
            assert len(a)>=3
            errs.extend([abs(row.effect-a.mean_expression.mean()),abs(row.mean_detection_fraction-a.detection.mean())]);evaluated+=1
        else:assert pd.isna(row.effect) and pd.isna(row.mean_detection_fraction)
    assert max(errs)<1e-10 and prof.p_value.isna().all() and prof.q_value.isna().all()
    rank=pd.read_csv(pub/'sc_source_stability.tsv',sep='\t');assert len(rank)==284 and not rank.top_celltype.eq('Unresolved').any()
    assert rank.bootstrap_top_frequency.dropna().between(0,1).all() and rank.bootstrap_valid.between(0,1000).all()
    primary=m[m.partition.eq('CARE2025_primary')];assert len(primary)==205880 and primary['Patient ID'].nunique()==55
    v=dict(status='PASS',all_private_patient_category_counts_rebuilt_from_author_barcodes=True,validated_public_profiles=evaluated,maximum_numeric_error=max(errs),primary_patients=55,primary_nuclei=205880,excluded_pretreated_primary_nuclei=795,bootstrap_bounds_checked=True,unresolved_never_ranked=True,primary_recurrent_not_independent=True)
    (pub/'independent_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)
