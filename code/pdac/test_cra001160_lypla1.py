"""Boundary and randomization tests independent of the observed gene result."""
import unittest
import numpy as np
from cra001160_lypla1 import bh,inference
class TestStatistics(unittest.TestCase):
 def test_bh_known_family(self):
  np.testing.assert_allclose(bh([.04,.001,.02]),[.04,.003,.03])
  self.assertTrue(np.isnan(bh([np.nan])[0]))
 def test_insufficient_units(self):
  self.assertEqual(inference([1]*4,[0]*4,True)['status'],'NOT_EVALUABLE')
 def test_zero_paired_difference(self):
  r=inference(np.arange(6),np.arange(6),True)
  self.assertEqual(r['effect'],0);self.assertEqual(r['p_value'],1)
  self.assertEqual(r['ci_lower'],0);self.assertEqual(r['ci_upper'],0)
 def test_exact_sign_probability(self):
  r=inference(np.ones(10),np.zeros(10),True)
  self.assertAlmostEqual(r['p_value'],2/1024,delta=.0006)
  self.assertEqual(r['n_higher'],10)
 def test_orientation(self):
  x=np.array([0,1,3,4,8,10.]);y=np.array([1,1,2,1,4,5.])
  for paired in [True,False]:
   a=inference(x,y,paired);b=inference(-x,-y,paired)
   self.assertAlmostEqual(a['effect'],-b['effect'])
   self.assertEqual(a['p_value'],b['p_value'])
if __name__=='__main__':unittest.main()
