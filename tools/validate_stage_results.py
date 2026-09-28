"""Validate common cancer stage layout and compare copied effects to frozen source."""
import argparse,csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--cancer',required=True);a=p.parse_args()
    schema=read(ROOT/'config/stage_result_schema_v1.json')
    root=ROOT/'results'/a.cancer
    assert root.resolve().parent==(ROOT/'results').resolve()
    index=read(root/'stage_index.json')
    assert index['cancer']==a.cancer
    assert [s['stage_id'] for s in index['stages']]==[s['stage_id'] for s in schema['stages']]
    frozen={(r['dataset'],r['cancer'],r['feature_name'],r['metabolite_key']):r for r in csv.DictReader((ROOT/'reference/camp/cancer_effects.tsv').open(encoding='utf-8-sig'),delimiter='\t')}
    n=0
    for stage in index['stages']:
        assert stage['stage_status'] in schema['statuses']
        if stage['path'] is None:continue
        batch=(ROOT/stage['path']).resolve()
        assert batch.is_relative_to(root.resolve())
        for name in schema['required_batch_files']:assert (batch/name).is_file()
        summary=read(batch/'summary.json');records=read(batch/'records.json')
        assert list(summary)==schema['summary_field_order']
        assert summary['stage_id']==stage['stage_id'] and summary['cancer']==a.cancer
        assert summary['status']==stage['batch_status'] and summary['status'] in schema['statuses']
        assert summary['publication']['status'] in schema['publication_statuses']
        assert len(records)==summary['record_count']==summary['planned_count']
        keys=[]
        for r in records:
            assert list(r)==schema['record_prefix']
            assert r['schema_version']==schema['schema_version'] and r['cancer']==a.cancer and r['stage_id']==stage['stage_id']
            assert r['record_status'] in schema['statuses']
            source=frozen[(r['dataset'],a.cancer,r['feature_name'],r['metabolite_key'])]
            assert (r['original_effect'],r['original_q'])==(source['hedges_g'],source['effect_fdr'])
            keys.append(tuple(r[k] or '' for k in schema['record_sort']))
        assert keys==sorted(keys) and len(set(keys))==len(keys)
        for entry in summary['inputs']:assert sha(ROOT/entry['path'])==entry['sha256']
        manifest=read(batch/'manifest.json')['files']
        assert [e['path'] for e in manifest]==['README_CN.md','summary.json','records.json']
        for entry in manifest:
            file=batch/entry['path'];assert file.stat().st_size==entry['bytes'] and sha(file)==entry['sha256']
        headings=[line[3:] for line in (batch/'README_CN.md').read_text(encoding='utf-8').splitlines() if line.startswith('## ')]
        assert headings==schema['report_sections']
        n+=1
    print(f'PASS: {a.cancer}, {n} batches; stage/field/row order, unique keys, counts, frozen values and manifests.')
if __name__=='__main__':main()
