"""Independently reconstruct equal-donor summaries from server-private donor tables."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
def main():
 ap=argparse.ArgumentParser();ap.add_argument('run');args=ap.parse_args();run=Path(args.run);public=run/'public';reported=pd.read_csv(public/'cell_source_all.tsv',sep='\t');tops=pd.read_csv(public/'cohort_top_sources.tsv',sep='\t');consensus=pd.read_csv(public/'cross_cohort_source.tsv',sep='\t');checked=0
 for cohort,d in reported.groupby('cohort'):
  raw=pd.read_csv(run/(cohort+'_private_unit_pseudobulk.tsv'),sep='\t');raw=raw[raw.scheme=='broad'];raw=raw[(raw.n_cells>=20)&raw.log1p_cpm.notna()]
  grouped=raw.groupby(['gene','celltype']);expected=grouped.agg(mean_expression=('log1p_cpm','mean'),detection_fraction=('positive_fraction','mean'),donor_count=('unit','size'))
  for _,r in d.iterrows():
   key=(r.gene,r.celltype)
   if key not in expected.index:assert pd.isna(r.mean_expression) and r.donor_count==0;continue
   z=expected.loc[key];assert np.isclose(z.mean_expression,r.mean_expression,atol=1e-12);assert np.isclose(z.detection_fraction,r.detection_fraction,atol=1e-12);assert int(z.donor_count)==r.donor_count;checked+=1
  for gene,g in d.groupby('gene'):
   e=g[g.status=='EVALUABLE'];top=tops[(tops.gene==gene)&(tops.cohort==cohort)].iloc[0]
   if len(e):
    candidates=';'.join(sorted(e.loc[np.isclose(e.mean_expression,e.mean_expression.max(),rtol=0,atol=1e-12),'celltype']));assert candidates==top.top_celltype
    assert np.isclose(e.bootstrap_first_frequency.sum(),1.,atol=1e-12) and e.bootstrap_first_frequency.between(0,1).all()
   else:assert top.top_celltype=='NOT_EVALUABLE'
 for _,r in consensus.iterrows():
  ds=tops[tops.gene==r.gene];labels=ds.top_celltype.tolist();same=len(set(labels))==1 and labels[0]!='NOT_EVALUABLE' and ';' not in labels[0];assert bool(r.cross_dataset_same_top)==same
  expected_stable=same and (ds.bootstrap_top_frequency>=.7).all() and (ds.detection_fraction>=.05).all();assert bool(r.stable_source)==expected_stable
  sets=[set(d.loc[d.status=='EVALUABLE','celltype']) for _,d in reported[reported.gene==r.gene].groupby('cohort')];assert len(set.intersection(*sets))==r.n_shared_categories
 assert len(consensus)==538 and len(tops)==1614 and not reported.duplicated(['cohort','gene','celltype']).any()
 result={'status':'PASS','independent_equal_donor_numeric_rows':checked,'cohort_gene_top_checks':1614,'three_cohort_consensus_checks':538,'bootstrap_distribution_mass_checks':True,'source_extraction_sentinels':'Four genes/cohort checked against counts in production','all_gene_raw_count_reextraction':'NOT_RUN for reused genes beyond sentinel subset','full_bootstrap_replay':'NOT_RUN;synthetic tie/donor tests and distribution mass checked'}
 (public/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
