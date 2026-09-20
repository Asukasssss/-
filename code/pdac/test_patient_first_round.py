"""Checks for fixed-family multiplicity, input gates and tied-data behavior."""
import tempfile,unittest
from pathlib import Path
import numpy as np
from patient_first_round import fixed_bh,compute,leave_one_out,validate_historical_hashes,sha
class StatisticalChecks(unittest.TestCase):
    def test_planned_family_missing_remains_public_na(self):
        q=fixed_bh([.01,None,.04])
        self.assertAlmostEqual(q[0],.03);self.assertIsNone(q[1]);self.assertAlmostEqual(q[2],.06)
    def test_tied_rank_monotone_input_invariance(self):
        x=np.array([0,0,0,1,1,2,3,4.],float);y=np.array([3,0,1,2,2,5,4,7.],float)
        a=compute(x,y,120);b=compute(2*x+10,4*y-3,120)
        for k in ['rho','p','lo','hi']:self.assertAlmostEqual(a[k],b[k])
    def test_uncomputable_and_loo(self):
        self.assertEqual(compute(np.arange(7.),np.arange(7.),1)['reason'],'N_LT_8')
        self.assertEqual(compute(np.ones(9),np.arange(9.),1)['reason'],'CONSTANT_VALUES')
        r=leave_one_out(np.arange(9.),np.arange(9.))
        self.assertEqual(r['loo_n_valid'],9);self.assertEqual(r['loo_sign_changes'],0);self.assertAlmostEqual(r['loo_rho_min'],1)
    def test_hash_mismatch_blocks_reuse(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source';p.write_text('initial');h=sha(p)
            expected={str(p):{'current':h,'historical':h}};validate_historical_hashes(expected)
            p.write_text('changed')
            with self.assertRaises(ValueError):validate_historical_hashes(expected)
if __name__=='__main__':unittest.main()
