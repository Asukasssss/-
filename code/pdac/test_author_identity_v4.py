import unittest
from author_identity_v4 import resolve_identity,MALIGNANT,NORMAL,UNRESOLVED
class IdentityTests(unittest.TestCase):
    def test_specific_normal_overrides_broad_malignant(self):
        self.assertEqual(resolve_identity('GSE242230','Malignant','Normal Epithelial')[0],NORMAL)
    def test_explicit_malignant_subtypes(self):
        for v in ['Malignant - Basal','Malignant - Classical']:self.assertEqual(resolve_identity('GSE242230','Malignant',v)[0],MALIGNANT)
    def test_unknown_is_not_malignant(self):
        for c,b,f in [('GSE263733','Ductal cell',None),('GSE278688','Ductal',None),('GSE242230','Malignant','novel subtype')]:self.assertEqual(resolve_identity(c,b,f)[0],UNRESOLVED)
    def test_other_lineages_not_reassigned(self):
        self.assertEqual(resolve_identity('GSE242230','Fibroblast','myCAF')[0],'Fibroblast')
if __name__=='__main__':unittest.main()
