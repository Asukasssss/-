"""Scientific invariants for donor aggregation, zero handling, and cell alignment."""
import gzip
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from run_cell_expression_v1 import summarize, read_lee

class CellOriginTests(unittest.TestCase):
    def test_donor_equal_and_measured_zero(self):
        meta = pd.DataFrame({'study': ['synthetic'] * 63, 'tissue': ['T'] * 63,
            'enrichment': ['unsorted'] * 63, 'technology': ['v2'] * 63,
            'lineage': ['A'] * 63, 'patient': ['p1'] * 20 + ['p2'] * 21 + ['p3'] * 22})
        values = np.array([1] * 20 + [0] * 21 + [4] * 22)
        with tempfile.TemporaryDirectory() as tmp:
            rows = summarize(meta, np.ones(63, dtype=int) * 10,
                {'HDC': values, 'GSTA4': values}, Path(tmp), 20, 5)
            self.assertEqual(len(rows), 4)
            for row in rows:
                self.assertEqual(row['pseudobulk_CPM_median'], 100000)
                self.assertEqual(row['detection_fraction_min'], 0)
                self.assertEqual(row['status'], 'PARTIAL')

    def test_small_group_not_exported_as_expression(self):
        meta = pd.DataFrame({'study': ['synthetic'] * 2, 'tissue': ['T'] * 2,
            'enrichment': ['unsorted'] * 2, 'technology': ['v2'] * 2,
            'lineage': ['Rare'] * 2, 'patient': ['p1', 'p2']})
        with tempfile.TemporaryDirectory() as tmp:
            rows = summarize(meta, np.array([10, 10]), {'HDC': np.array([0, 1]),
                'GSTA4': np.array([1, 0])}, Path(tmp), 20, 5)
            self.assertTrue(all(row['effect'] == 'NA' for row in rows))
            self.assertTrue(all(row['status'] == 'NOT_EVALUABLE' for row in rows))
            self.assertEqual([row['n'] for row in rows], [2, 2, 0, 0])

    def test_matrix_identity_reordering_and_all_gene_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'synthetic.txt.gz'
            with gzip.open(path, 'wt') as f:
                f.write('Index\tcellB\tcellA\nHDC\t1\t2\nGSTA4\t3\t4\n')
                for i in range(10000):
                    f.write(f'GENE{i}\t1\t0\n')
            total, targets, checks = read_lee(path, np.array(['cellA', 'cellB']))
            np.testing.assert_array_equal(total, [6, 10004])
            np.testing.assert_array_equal(targets['HDC'], [2, 1])
            self.assertEqual(checks['genes'], 10002)

if __name__ == '__main__':
    unittest.main()
