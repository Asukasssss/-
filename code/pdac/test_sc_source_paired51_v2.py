import unittest
import pandas as pd
from sc_source_paired51_v2 import summarize_pb
class BootstrapTests(unittest.TestCase):
 def test_donors_equal_and_top_stable_not_cell_weighted(self):
  rows=[]
  for i in range(4):
   for ct,value in [('A',5.),('B',1.)]:rows.append(dict(scheme='broad',unit=str(i),gene='G',celltype=ct,n_cells=20 if i else 2000,log1p_cpm=value,sum_counts=10,positive_fraction=.1))
  agg,top=summarize_pb(pd.DataFrame(rows),'TEST',['G'],123)
  self.assertEqual(top.iloc[0].top_celltype,'A');self.assertEqual(top.iloc[0].second_celltype,'B');self.assertEqual(top.iloc[0].bootstrap_top_frequency,1.)
  self.assertEqual(agg.loc[agg.celltype=='A','mean_expression'].iloc[0],5.)
 def test_tied_top_splits_bootstrap_mass(self):
  rows=[dict(scheme='broad',unit=str(i),gene='G',celltype=c,n_cells=20,log1p_cpm=2.,sum_counts=1,positive_fraction=.1) for i in range(3) for c in ['A','B']]
  agg,top=summarize_pb(pd.DataFrame(rows),'TEST',['G'],123)
  self.assertEqual(top.iloc[0].top_celltype,'A;B');self.assertTrue((agg.bootstrap_first_frequency==.5).all())
if __name__=='__main__':unittest.main()
