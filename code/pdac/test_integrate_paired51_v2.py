import unittest
from integrate_paired51_v2 import relation_label
class LabelTests(unittest.TestCase):
 def test_missing_sensitivity_does_not_become_strong(self):
  a=dict(status='DONE',effect=.7,p_value=.001,q_value=.01,loo_sign_changes=0);b=dict(status='NOT_EVALUABLE')
  self.assertEqual(relation_label(a,b),'R-suggestive')
 def test_sensitivity_precedes_significance(self):
  a=dict(status='DONE',effect=.7,p_value=.001,q_value=.01,loo_sign_changes=0)
  self.assertEqual(relation_label(a,dict(status='DONE',effect=.2)),'R-sensitive')
  self.assertEqual(relation_label(a,dict(status='DONE',effect=.65)),'R-strong')
 def test_weak_not_removed_and_not_evaluable_distinct(self):
  a=dict(status='DONE',effect=.1,p_value=.8,q_value=.9,loo_sign_changes=0)
  self.assertEqual(relation_label(a,dict(status='DONE',effect=.15)),'R-weak')
  self.assertEqual(relation_label(dict(status='NOT_EVALUABLE'),{}),'R-not-evaluable')
if __name__=='__main__':unittest.main()
