"""Public-sequence exact matching; does not infer paper-specific V1/V2 identity."""
import argparse
import hashlib
import json
import urllib.request
from pathlib import Path
import pandas as pd

URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_019023.5,NM_001184824.4&rettype=fasta&retmode=text'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--probes', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    text = urllib.request.urlopen(URL, timeout=30).read().decode()
    assert text.startswith('>')
    sequences = {}
    for block in text.split('>')[1:]:
        lines = block.splitlines()
        sequences[lines[0].split()[0]] = ''.join(lines[1:])
    assert set(sequences) == {'NM_019023.5', 'NM_001184824.4'}
    rows = []
    for _, probe in pd.read_csv(a.probes, sep='\t').iterrows():
        for accession, sequence in sequences.items():
            forward = probe.SEQUENCE
            reverse = forward.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            row = {'probe': probe.NAME, 'refseq_version': accession, 'sequence_length': len(sequence),
                'paper_V1_V2_mapping': 'NEEDS_REVIEW'}
            for label, query in [('forward', forward), ('reverse_complement', reverse)]:
                hits = [str(i + 1) for i in range(len(sequence) - len(query) + 1) if sequence[i:i + len(query)] == query]
                row[label + '_exact_60nt_positions_1based'] = ';'.join(hits) or 'NA'
            rows.append(row)
    pd.DataFrame(rows).to_csv(a.out / 'PRMT7_refseq_exact_matches.tsv', sep='\t', index=False)
    (a.out / 'PRMT7_refseq_source.json').write_text(json.dumps({'url': URL, 'fasta_sha256': hashlib.sha256(text.encode()).hexdigest(),
        'accessed_utc_date': '2026-09-20', 'versions': list(sequences), 'method': 'Exact 60nt forward and reverse-complement search in two explicitly versioned reference transcripts; not genome-wide cross-hybridization analysis',
        'interpretation': 'Both probes match both references; these two probes cannot distinguish these two reference transcripts. Paper V1/V2 mapping remains unresolved.'}, indent=2))

if __name__ == '__main__':
    main()
