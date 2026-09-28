import hashlib,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from internal_paired51_v2 import validate_sources,resolve_genes
from patient_first_round import fixed_bh,validate_rna_labels
class GuardTests(unittest.TestCase):
 def test_exact_manifest_rejects_empty_missing_wrong_and_content_change(self):
  with tempfile.TemporaryDirectory() as d:
   paths=[Path(d)/str(i) for i in range(4)]
   for p in paths:p.write_text('source'+p.name)
   valid={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
   validate_sources(valid,paths)
   wrong=dict(valid);wrong[str(Path(d)/'wrong')]=wrong.pop(str(paths[0]))
   for bad in [{},{k:v for k,v in valid.items() if k!=str(paths[0])},wrong]:
    with self.assertRaises(ValueError):validate_sources(bad,paths)
   paths[0].write_text('changed')
   with self.assertRaises(ValueError):validate_sources(valid,paths)
 def test_blank_rna_is_not_valid_even_if_unique(self):
  for labels in [['A',''],['A',' '],['A',None],['A','A']]:
   with self.assertRaises(ValueError):validate_rna_labels(pd.DataFrame([[1],[2]],index=labels))
 def test_aliases_do_not_duplicate_one_measured_gene(self):
  found,reasons=resolve_genes(['NEW1','NEW2'],['OLD'],{'NEW1':['OLD'],'NEW2':['OLD']})
  self.assertEqual(found,{'NEW1':None,'NEW2':None})
  found,_=resolve_genes(['LARS1'],['LARS'],{'LARS1':['LARS']});self.assertEqual(found['LARS1'],'LARS')
 def test_fixed_family_keeps_missing_hypothesis(self):
  self.assertEqual(fixed_bh([.01,None,.1]),[.03,None,.15000000000000002])
if __name__=='__main__':unittest.main()
