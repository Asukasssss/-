"""Check tracked files before sharing. Does not replace manual data review."""
import csv, hashlib, subprocess
from pathlib import Path

root=Path(__file__).resolve().parents[1]
paths=subprocess.check_output(['git','ls-files','-z'],cwd=root).decode().split('\0')
bad=[]
for name in filter(None,paths):
    p=root/name
    if p.stat().st_size>5_000_000: bad.append(name+': exceeds 5 MB')
    if p.suffix.lower() in {'.xlsx','.gz','.zip','.h5ad','.bam','.pem','.key','.bundle'}: bad.append(name+': disallowed format')
    if any(x in name.lower() for x in ['sample_mapping','patient_values','candidate_protein_values']): bad.append(name+': patient-level filename')
    if p.name.startswith('.env'): bad.append(name+': environment credentials risk')
manifest=root/'reference/MANIFEST_SHA256.tsv'
if manifest.exists():
    for row in csv.DictReader(manifest.open(encoding='utf-8'),delimiter='\t'):
        p=root/row['repository_path']
        if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
            bad.append(row['repository_path']+': snapshot hash mismatch')
if bad: raise SystemExit('\n'.join(bad))
print('PASS: tracked file size/type checks and baseline snapshot hashes. Manual review still required for new content.')
