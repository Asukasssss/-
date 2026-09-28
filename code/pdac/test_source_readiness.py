"""Fail closed on source changes and ambiguous RNA rows, including unused genes."""
import unittest
import pandas as pd
from source_readiness import validate_rna_labels, validate_source_hashes

class SourceGuards(unittest.TestCase):
    def test_duplicate_unselected_gene_blocks(self):
        with self.assertRaises(ValueError):
            validate_rna_labels(pd.Index(['PNP', 'OTHER', 'OTHER']))

    def test_incomplete_gene_labels_block(self):
        for labels in [['PNP', None], ['PNP', ' ']]:
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                validate_rna_labels(pd.Index(labels))

    def test_changed_missing_or_inconsistent_hash_blocks(self):
        valid={str(i): dict(current='abc', historical='abc', equal=True) for i in range(3)}
        validate_source_hashes(valid)
        for change in [dict(current='changed'), dict(historical=None), dict(equal=False)]:
            bad={k: dict(v) for k,v in valid.items()};bad['0'].update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_source_hashes(bad)
        with self.assertRaises(ValueError):
            validate_source_hashes({})
        with self.assertRaises(ValueError):
            validate_source_hashes({'0': valid['0']})

    def test_unique_labels_pass(self):
        validate_rna_labels(pd.Index(['PNP', 'OTHER']))

if __name__=='__main__': unittest.main()
