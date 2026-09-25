import gzip,tempfile,unittest
from pathlib import Path
import numpy as np,pandas as pd
from sc_source_three_cohorts_v1 import sparse_extract,summarize
class SourceTests(unittest.TestCase):
 def test_sparse_barcode_reorder_ambiguous_symbol_and_full_library(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'matrix.gz'
   with gzip.open(p,'wt') as f:f.write('%%MatrixMarket matrix coordinate integer general\n4 3 5\n1 1 2\n2 1 3\n3 1 7\n3 3 11\n4 3 13\n')
   meta=pd.DataFrame(index=['c','a']);a,lib,present=sparse_extract(p,['AMB','AMB','TARGET','OTHER'],['a','b','c'],meta,['AMB','TARGET'])
   self.assertEqual(present,{'TARGET'});np.testing.assert_array_equal(lib,[24,12]);np.testing.assert_array_equal(a,[[0,0],[11,7]])
 def test_missing_metadata_barcode_rejected(self):
  with self.assertRaises(AssertionError):sparse_extract(Path('unused'),['g'],['a'],pd.DataFrame(index=['absent']),['g'])
 def test_equal_unit_mean_not_cell_weighted_and_low_detection(self):
  with tempfile.TemporaryDirectory() as d:
   unit=['small1']*20+['small2']*20+['large']*100
   meta=pd.DataFrame({'unit':unit,'celltype':['Ductal']*140},index=[str(i) for i in range(140)])
   a=np.array([[10]*40+[1]*100,[0]*40+[1]*100]);lib=np.full(140,100)
   out,_,_=summarize('SYNTHETIC',meta,a,lib,{'HIGH','RARE'},['HIGH','RARE'],Path(d))
   high=out[(out.annotation_level=='broad')&(out.gene=='HIGH')].iloc[0]
   self.assertAlmostEqual(high.mean_unit_log1p_cpm,(2*np.log1p(100000)+np.log1p(10000))/3)
   self.assertEqual(high.n_units_ge20cells,3);self.assertEqual(high.status,'EVALUABLE')
   self.assertEqual(out[(out.annotation_level=='broad')&(out.gene=='RARE')].iloc[0].status,'LOW_DETECTION')
if __name__=='__main__':unittest.main()
