import unittest
import numpy as np
from run_cell_paired35_v1 import paired_stats,bh

class Tests(unittest.TestCase):
    def test_exact_sign_and_reversal(self):
        a=paired_stats(np.arange(1,7),np.zeros(6));b=paired_stats(np.zeros(6),np.arange(1,7))
        self.assertEqual(a['p_value'],.03125);self.assertEqual(a['effect'],3.5)
        self.assertEqual(a['p_value'],b['p_value']);self.assertEqual(a['effect'],-b['effect'])
        self.assertEqual((a['ci_lower'],a['ci_upper']),(1,6));self.assertGreaterEqual(a['ci_coverage'],.95)
    def test_ties_and_bh(self):
        a=paired_stats(np.zeros(10),np.zeros(10));self.assertEqual(a['p_value'],1);self.assertEqual(a['n_informative'],0)
        a=paired_stats([1,2,3,0,0,0],[0]*6);self.assertEqual(a['p_value'],.25);self.assertEqual(a['n_zero'],3)
        np.testing.assert_allclose(bh([.01,.04,.03,1]),[.04,.05333333333333334,.05333333333333334,1])
if __name__=='__main__':unittest.main()
