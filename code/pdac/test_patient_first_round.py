"""Checks for fixed-family multiplicity, input gates and tied-data behavior."""
import tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from unittest.mock import patch
import patient_first_round as patient
from patient_first_round import fixed_bh,compute,leave_one_out,validate_historical_hashes,validate_rna_labels,sha
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
            paths=[Path(d)/name for name in ['mapping','metabolomics','rna']]
            for p in paths:p.write_text('initial')
            expected={str(p):{'current':sha(p),'historical':sha(p)} for p in paths}
            validate_historical_hashes(expected,paths)
            paths[0].write_text('changed')
            with self.assertRaises(ValueError):validate_historical_hashes(expected,paths)
    def test_missing_or_wrong_source_keys_block_before_file_read(self):
        paths=[Path(name) for name in ['mapping','metabolomics','rna']]
        valid={str(p):{'current':'abc','historical':'abc'} for p in paths}
        cases=[{}, {k:v for k,v in valid.items() if k!='rna'},
               {**{k:v for k,v in valid.items() if k!='rna'},'wrong':valid['rna']},
               {**valid,'extra':valid['rna']}]
        for expected in cases:
            with self.subTest(keys=list(expected)), patch.object(patient,'sha') as file_hash:
                with self.assertRaises(ValueError):validate_historical_hashes(expected,paths)
                file_hash.assert_not_called()
    def test_missing_hash_fields_block(self):
        paths=[Path(name) for name in ['mapping','metabolomics','rna']]
        for entry in [{},{'current':'abc'},{'historical':'abc'},None]:
            expected={str(p):entry for p in paths}
            with self.subTest(entry=entry),self.assertRaises(ValueError):
                validate_historical_hashes(expected,paths)
    def test_rna_missing_blank_duplicate_labels_block(self):
        for labels in [['PNP',None],['PNP',np.nan],['PNP',''],['PNP','  '],['PNP','PNP']]:
            for axis in ['index','columns']:
                frame=pd.DataFrame(np.ones((2,2)),index=['PNP','OTHER'],columns=['T1','T2'])
                setattr(frame,axis,pd.Index(labels))
                with self.subTest(labels=labels,axis=axis),self.assertRaises(ValueError):
                    validate_rna_labels(frame)
        validate_rna_labels(pd.DataFrame([[1]],index=['PNP'],columns=['T1']))
    def test_run_empty_manifest_stops_before_matrix_read(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'results/collaborative/PDAC/B/test';out.mkdir(parents=True)
            (out/'input_hash_comparison.json').write_text('{}')
            with patch.object(patient,'ROOT',root),patch.object(patient,'__file__',str(out/'patient_first_round.py')),patch.object(patient.pd,'read_csv') as matrix_read:
                with self.assertRaisesRegex(ValueError,'three actual input paths'):patient.run(None)
                matrix_read.assert_not_called()
            self.assertTrue((out/'FAILED.txt').exists())
            self.assertFalse((out/'DONE.json').exists())
if __name__=='__main__':unittest.main()
