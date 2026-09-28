"""Inspect GSE200997 GEO annotation on server165; export aggregates only."""
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

source, output = map(Path, sys.argv[1:3])
with gzip.open(source, 'rt') as handle:
    reader = csv.DictReader(handle)
    fields = reader.fieldnames
    rows = list(reader)
result = {
    'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
    'fields': fields,
    'cells': len(rows),
    'unique_cell_ids': len({r[''] for r in rows}),
    'samples': len({r['samples'] for r in rows}),
    'condition_cells': dict(Counter(r['Condition'] for r in rows)),
    'condition_samples': {c: len({r['samples'] for r in rows if r['Condition'] == c})
                          for c in {r['Condition'] for r in rows}},
    'cell_type_column_present': False,
    'cell_type_review': 'Listed columns contain CMS predictions, not cell-type annotation.',
    'gene_expression_analysis': 'NOT_RUN',
}
assert result['cells'] == result['unique_cell_ids']
output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(result))
