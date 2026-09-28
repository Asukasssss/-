"""Small independent checks of the donor-unit comparison, no matrices required."""
import unittest
import numpy as np
from myeloid_infercnv_summarize import signflip_p,bootstrap_ci
class DonorComparison(unittest.TestCase):
    def test_exact_extreme_probability(self):
        self.assertEqual(signflip_p([1]*6),2/64)
        self.assertEqual(signflip_p([-1]*6),2/64)
    def test_null_and_minimum_units(self):
        self.assertEqual(signflip_p([0]*6),1)
        self.assertTrue(np.isnan(signflip_p([1]*4)))
        self.assertTrue(all(np.isnan(bootstrap_ci([1]*4))))
    def test_nonfinite_rejected(self):
        for f in [signflip_p,bootstrap_ci]:
            with self.assertRaises(ValueError):f([1,2,3,4,np.nan])
    def test_constant_bootstrap(self):
        self.assertEqual(bootstrap_ci([2]*6),[2,2])
if __name__=='__main__':unittest.main()
