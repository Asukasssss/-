"""Reproduce metadata version comparison and public platform-probe extraction."""
import argparse
import gzip
import io
import json
from pathlib import Path
import pandas as pd
from cell_origin_v1 import load_metadata, sha, ROOT

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    assert a.data_dir.resolve().parent == ROOT / 'data/candidates'
    assert a.out.resolve().parent == ROOT / 'results/collaborative/COAD/B'
    d, r = a.data_dir, a.out
    z = load_metadata(d)['Pelka2021_GSE178341']
    g = pd.read_csv(d / 'pelka_geo_metadata_retry.csv.gz', dtype=str, keep_default_na=False).set_index('cellID')
    m = pd.read_csv(d / 'pelka_author_metadata.csv.gz', dtype=str, keep_default_na=False)
    m.index = z.cell_id
    assert g.index.is_unique and m.index.is_unique
    common = g.index.intersection(m.index)
    eq = {col: bool((g.loc[common, col] == m.loc[common, col]).all()) for col in g.columns}
    assert all(eq.values()) and len(common) == len(g)
    extra = m.loc[~m.index.isin(g.index)]
    audit = {'author_cells': len(m), 'geo_cells': len(g), 'common': len(common),
        'all_metadata_fields_equal_on_common': eq,
        'author_only_processing': extra.PROCESSING_TYPE.value_counts().to_dict(), 'geo_only': 0}
    (r / 'pelka_version_audit.json').write_text(json.dumps(audit, indent=2))
    lines, on = [], False
    with gzip.open(d / 'GPL16699.soft.gz', 'rt') as f:
        for line in f:
            if line.startswith('!platform_table_begin'):
                on = True
                continue
            if line.startswith('!platform_table_end'):
                break
            if on:
                lines.append(line)
    table = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', dtype=str, keep_default_na=False)
    probes = table[table.GENE_SYMBOL.str.contains(r'(?:^|[ ,;/])PRMT7(?:$|[ ,;/])', regex=True)]
    assert len(probes) > 0 and probes.SEQUENCE.str.fullmatch('[ACGT]+').all()
    probes.to_csv(r / 'PRMT7_platform_probes.tsv', sep='\t', index=False)
    print(json.dumps({'metadata_audit': audit, 'PRMT7_probes': len(probes), 'platform_sha256': sha(d / 'GPL16699.soft.gz')}))

if __name__ == '__main__':
    main()
