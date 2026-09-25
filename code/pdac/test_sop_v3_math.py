import unittest
import numpy as np
from sop_v3_math import bh_evaluable,paired_t,cellwise_profiles,source_rank_bootstrap

class SOPMath(unittest.TestCase):
    def test_bh_missing_excluded(self):
        np.testing.assert_allclose(bh_evaluable([.01,np.nan,.04,.03]),[.03,np.nan,.04,.04],equal_nan=True)
        with self.assertRaises(ValueError):bh_evaluable([1.1])
    def test_t_matches_scipy_independent(self):
        from scipy.stats import ttest_1samp
        d=np.array([1,2,-1,0,4,6,2,3,4,2,1.])
        z=paired_t(d,13);self.assertAlmostEqual(z['p_value'],ttest_1samp(d,0).pvalue,14)
        self.assertEqual(z['n_up']+z['n_down']+z['n_equal'],len(d));self.assertEqual(z['bootstrap_valid'],4000)
        self.assertEqual(paired_t(d[:7],13)['status'],'NOT_EVALUABLE')
        self.assertEqual(paired_t(np.ones(8),13)['reason'],'ZERO_VARIANCE_PAIRED_DIFFERENCES')
    def test_cellwise_not_pseudobulk(self):
        x=np.array([[1,9,2],[0,0,0]]);lib=np.array([10,100,5]);p=cellwise_profiles(x,lib,[[0,1],[2]],np.array([True,False]))
        self.assertAlmostEqual(p[0][1][0],np.mean(np.log1p([1000,900])))
        self.assertNotAlmostEqual(p[0][1][0],np.log1p(10000*10/110),places=5)
        self.assertTrue(np.isnan(p[0][1][1]));self.assertEqual(p[0][2][0],1)
    def test_source_ties_low_signal_and_coverage(self):
        x=np.ones((5,2));z=source_rank_bootstrap(x,x*.2,11)
        self.assertEqual(z['top'],[0,1]);np.testing.assert_allclose(z['frequency'],[.5,.5])
        self.assertEqual(z['valid'],1000)
        self.assertEqual(source_rank_bootstrap(x,x*.009,11)['status'],'NOT_EVALUABLE')
        y=x.copy();y[:3,1]=np.nan
        self.assertEqual(source_rank_bootstrap(y,y*.2,11)['status'],'NOT_EVALUABLE')
    def test_unequal_cell_counts_still_equal_donor(self):
        x=np.array([[1.,.5],[4.,.5],[4.,.5]])
        z=source_rank_bootstrap(x,np.full_like(x,.1),4)
        self.assertAlmostEqual(z['mean'][0],3);self.assertEqual(z['top'],[0])

if __name__=='__main__':unittest.main()
