import unittest
import numpy as np
import pandas as pd
from paired_metabolites_v1 import paired_stat,bh_fixed,pair_join
class PairedChecks(unittest.TestCase):
    def test_exact_all_positive_and_reverse(self):
        a=paired_stat(np.arange(1.,7),1);b=paired_stat(-np.arange(1.,7),1)
        self.assertEqual(a['p_value'],2/64);self.assertEqual(a['p_value'],b['p_value'])
        self.assertEqual(a['rank_biserial'],1);self.assertEqual(b['rank_biserial'],-1)
        self.assertEqual(a['mean_delta'],-b['mean_delta'])
    def test_ties_zeros_and_low_n(self):
        self.assertEqual(paired_stat(np.ones(6),1)['p_value'],2/64)
        self.assertEqual(paired_stat([0,0,0,1,1,1],1)['p_value'],2/8)
        self.assertEqual(paired_stat(np.zeros(6),1)['p_value'],1)
        self.assertEqual(paired_stat(np.arange(5),1)['status'],'NOT_EVALUABLE')
    def test_fixed_family(self):
        q=bh_fixed([.01,None,.04]);self.assertAlmostEqual(q[0],.03);self.assertEqual(q[1],'NA');self.assertAlmostEqual(q[2],.06)
    def fixture(self):
        m=pd.DataFrame({'Dataset':['PDAC']*3,'CommonID':['a','b','c'],'MetabID':['m1','m2','m3'],'RNAID':['GSM1_file','GSM2_file','GSM3_file'],'TN':['Tumor','Normal','Tumor']})
        a=pd.DataFrame({'GSE ID':['GSM1','GSM2','GSM3','GSM4'],'Sample':['s1','s2','s3','s4'],'T/N':['T','N','T','N'],'Paired Sample ID':['pA','pA','pB','pB']})
        return m,a
    def test_explicit_pair_and_unmatched_exclusion(self):
        m,a=self.fixture();_,p=pair_join(m.sample(frac=1,random_state=1),a.sample(frac=1,random_state=2))
        self.assertEqual(len(p),1);self.assertEqual(p.iloc[0].tumor_metab_id,'m1');self.assertEqual(p.iloc[0].normal_metab_id,'m2')
    def test_missing_pair_id_does_not_infer_from_order(self):
        m,a=self.fixture();a['Paired Sample ID']=None
        with self.assertRaises(ValueError):pair_join(m,a)
    def test_tissue_disagreement_blocks(self):
        m,a=self.fixture();m.loc[0,'TN']='Normal'
        with self.assertRaises(ValueError):pair_join(m,a)
    def test_duplicate_author_accession_blocks(self):
        m,a=self.fixture();a.loc[1,'GSE ID']='GSM1'
        with self.assertRaises(ValueError):pair_join(m,a)
if __name__=='__main__':unittest.main()
