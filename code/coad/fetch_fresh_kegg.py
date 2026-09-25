"""Fetch new KEGG evidence starting solely from frozen COAD significant effects.

Public database responses are cached outside tracked reference snapshots.
This prepares evidence, not a validated gene shortlist or new statistical tests.
"""
import csv
import hashlib
import json
import re
import time
import urllib.request
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / 'runtime/coad_fresh_kegg_20260919'
CACHE.mkdir(parents=True, exist_ok=True)
REGISTRY = CACHE / 'sources.json'
sources = json.loads(REGISTRY.read_text(encoding='utf-8')) if REGISTRY.exists() else {}


def fetch(endpoint):
    url = 'https://rest.kegg.jp/' + endpoint
    path = CACHE / (hashlib.sha256(url.encode()).hexdigest() + '.txt')
    if path.exists() and url in sources:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == sources[url]['sha256']
        return path.read_text(encoding='utf-8')
    for attempt in range(3):
        try:
            time.sleep(.4)
            with urllib.request.urlopen(url, timeout=35) as response:
                data = response.read()
            path.write_bytes(data)
            sources[url] = {'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                            'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data),
                            'cache_file': path.name, 'status': 'DONE'}
            REGISTRY.write_text(json.dumps(sources, indent=2), encoding='utf-8')
            return data.decode('utf-8')
        except Exception as exc:
            if attempt == 2:
                sources[url] = {'url': url, 'status': 'ACCESS_BLOCKED', 'error': str(exc)}
                REGISTRY.write_text(json.dumps(sources, indent=2), encoding='utf-8')
                raise
            time.sleep(1 + attempt)


def parse(text):
    result = {}
    for block in text.split('///'):
        fields = {}
        field = None
        for line in block.strip('\r\n').splitlines():
            label, value = line[:12].strip(), line[12:].strip()
            if label:
                field = label
                fields.setdefault(field, [])
            if field:
                fields[field].append(value)
        if 'ENTRY' in fields:
            ident = fields['ENTRY'][0].split()[0]
            result[ident] = fields
    return result


def get_entries(ids):
    records = {}
    ids = sorted(set(ids))
    for start in range(0, len(ids), 10):
        batch = ids[start:start+10]
        records.update(parse(fetch('get/' + '+'.join(batch))))
        print(f'Fetched {min(start+10, len(ids))}/{len(ids)} entries', flush=True)
    if set(ids) - set(records):
        print('ENTRIES_NOT_RETURNED: ' + ','.join(sorted(set(ids) - set(records))), flush=True)
    return records


def main():
    source = ROOT / 'reference/camp/cancer_effects.tsv'
    rows = [r for r in csv.DictReader(source.open(encoding='utf-8-sig'), delimiter='\t')
            if r['cancer'] == 'COAD' and float(r['effect_fdr']) < .05]
    assert len(rows) == 73
    compounds = get_entries([r['kegg_id'] for r in rows if re.fullmatch(r'C\d{5}', r['kegg_id'])])
    # Independently query exceptional names and the suspicious reaction ID.
    exceptions = {name: fetch('find/compound/' + urllib.parse.quote(name))
                  for name in ['ophthalmate', 'ophthalmic acid', 'glutathione cysteine', 'xylulose 5-phosphate', 'mucate', 'galactarate']}
    exceptions['original_R00900'] = fetch('get/R00900')
    # Follow-up hypotheses are kept separate from original identifiers.
    # C05526 from the written R00900 equation; C21016/C00879/C00231 from name searches.
    # C00217 is the D-glutamate hypothesis revealed by the original HMDB ID.
    supplemental = get_entries(['C05526', 'C21016', 'C00879', 'C00231', 'C00217'])
    compounds.update(supplemental)
    (CACHE / 'compound_evidence.json').write_text(json.dumps(compounds, ensure_ascii=False, indent=2), encoding='utf-8')
    (CACHE / 'exception_evidence.json').write_text(json.dumps(exceptions, ensure_ascii=False, indent=2), encoding='utf-8')
    human_kos = {}
    for line in fetch('link/ko/hsa').splitlines():
        gene, ko = line.split('\t')
        human_kos.setdefault(ko, []).append(gene)
    reaction_kos = {}
    for line in fetch('link/ko/reaction').splitlines():
        reaction, ko = line.split('\t')
        if ko in human_kos:
            reaction_kos.setdefault(reaction.removeprefix('rn:'), []).append(ko)
    all_reactions = set()
    for record in compounds.values():
        all_reactions.update(re.findall(r'R\d{5}', ' '.join(record.get('REACTION', []))))
    selected = sorted(all_reactions & reaction_kos.keys())
    print(f'{len(compounds)} compounds, {len(all_reactions)} reactions, {len(selected)} with human KO links', flush=True)
    reactions = get_entries(selected)
    human_names = {}
    for line in fetch('list/hsa').splitlines():
        fields = line.split('\t')
        human_names[fields[0]] = fields[-1]
    output = {'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'compound_records': compounds, 'reaction_records': reactions,
              'reaction_human_kos': {r: reaction_kos[r] for r in selected},
              'human_ko_genes': human_kos, 'human_gene_names': human_names,
              'exceptions': exceptions}
    (CACHE / 'evidence.json').write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding='utf-8')
    print('DONE: fresh evidence cache', flush=True)


if __name__ == '__main__':
    main()
