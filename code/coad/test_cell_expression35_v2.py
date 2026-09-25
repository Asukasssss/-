"""Regression against v1 patient-weighting plus explicit unavailable-gene handling."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from run_cell_expression35_v2 import summarize, symbol_choice
from run_cell_expression_v1 import summarize as old_summarize

class ExpansionTests(unittest.TestCase):
    def test_alias_and_paralog_boundaries(self):
        names = ['ADSS', 'GSTT2B', 'X', 'X']
        self.assertEqual(symbol_choice('ADSS2', names, {'ADSS2': ['ADSS']})[:2], (0, 'ADSS'))
        self.assertIsNone(symbol_choice('GSTT2', names, {})[0])
        self.assertEqual(symbol_choice('X', names, {})[2], 'AMBIGUOUS_SOURCE_SYMBOL')
        self.assertEqual(symbol_choice('ADSS2', ['ADSS', 'ADSS2'], {'ADSS2': ['ADSS']})[:2], (1, 'ADSS2'))

    def test_vector_quantiles_match_v1_and_missing_is_not_zero(self):
        meta = pd.DataFrame({'study': ['test'] * 63, 'tissue': ['T'] * 63,
            'enrichment': ['unsorted'] * 63, 'technology': ['v2'] * 63,
            'lineage': ['A'] * 63, 'patient': ['p1'] * 20 + ['p2'] * 21 + ['p3'] * 22})
        counts = {'HDC': np.array([1] * 20 + [0] * 21 + [4] * 22), 'GSTA4': np.arange(63) % 4}
        total = np.ones(63, dtype=np.int64) * 10
        with tempfile.TemporaryDirectory() as tmp:
            newout = Path(tmp) / 'new'; oldout = Path(tmp) / 'old'; newout.mkdir(); oldout.mkdir()
            old = old_summarize(meta, total, counts, oldout, 20, 5)
            new = summarize(meta, total, counts, ['HDC', 'GSTA4', 'GSTT2'], newout)
            for row in old:
                candidate = next(x for x in new[row['gene']] if x['support_set'] == row['support_set'])
                for col, val in row.items():
                    if col in ['run_id', 'analysis_version']:
                        continue
                    if isinstance(val, float):
                        self.assertAlmostEqual(val, candidate[col], places=10)
                    else:
                        self.assertEqual(val, candidate[col])
            self.assertTrue(all(x['effect'] == 'NA' and x['status'] == 'NOT_EVALUABLE' for x in new['GSTT2']))

if __name__ == '__main__':
    unittest.main()
