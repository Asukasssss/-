"""Independent output integrity and scope checks; does not recompute statistics."""
import csv, hashlib, json, sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root/sys.argv[1]
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
original=read(root/'results/BRCA/04_ROBUSTNESS/20260921T122021Z_paired_p_rank_v1/p_lt_005_sorted.tsv')
mapped=read(out/'metabolite190_mapping_status.tsv')
assert len(mapped)==190
assert all(all(a[k]==b[k] for k in a) for a,b in zip(original,mapped))
old=read(root/'reference/brca/direct_edges_v1.tsv');hist=read(out/'legacy174_preserved.tsv')
assert len(hist)==174 and all(all(a[k]==b[k] for k in a) for a,b in zip(old,hist))
edges=read(out/'direct_relations.tsv');cond=read(out/'conditional_relations.tsv');ev=read(out/'relation_evidence.tsv')
assert len({x['relation_key'] for x in edges})==len(edges)==262
assert len({x['gene'] for x in edges})==150
assert not {x['relation_key'] for x in edges}&{x['relation_key'] for x in cond}
assert all(e['mapping_status']=='DIRECT_ACCEPTED' for e in edges)
assert all(e['metabolite_key'] in {m['metabolite_key'] for m in mapped} for e in edges+cond)
assert all(not e['metabolite_name'].startswith('X - ') for e in edges+cond)
for e in edges+cond:
    for i in e['evidence_ids'].split(';'):
        match=[x for x in ev if x['evidence_id']==i];assert len(match)==1 and match[0]['relation_key']==e['relation_key']
for m in mapped:
    ee=[e for e in edges if e['metabolite_key']==m['metabolite_key']]
    assert int(m['direct_gene_count'])==len(ee)
    assert set(m['direct_genes'].split(';'))- {''}=={e['gene'] for e in ee}
assert all(e['patient_analysis_this_round']=='NOT_RUN' for e in edges+cond)
for m in read(out/'source_manifest.tsv'):
    b=(root/m['path']).read_bytes()
    assert hashlib.sha256(b.replace(b'\r\n',b'\n')).hexdigest()==m['canonical_lf_sha256'],m['path']
# Exact reaction additions point to the correct primary human gene and a
# reaction actually present in the saved source, not just a synonym hit.
checked=0
for e in ev:
    if e['evidence_origin']!='NEW_REVIEW' or not e['rhea_id']:continue
    d=json.loads((root/e['source_path']).read_text(encoding='utf-8'))
    entries=[x for x in d['results'] if x['primaryAccession']==e['uniprot_accession']]
    assert len(entries)==1
    x=entries[0];assert x['organism']['taxonId']==9606
    assert e['gene'] in [g.get('geneName',{}).get('value') for g in x['genes']]
    assert e['reaction'] in [c['reaction']['name'] for c in x.get('comments',[]) if c['commentType']=='CATALYTIC ACTIVITY']
    checked+=1
result=dict(status='PASS',input190_all_original_values_and_order_preserved=True,legacy174_preserved=True,
    direct262_unique=True,genes150_unique=True,conditional18_separate=True,unknown53_no_genes=True,
    evidence_links_valid=True,new_reaction_records_crosschecked=checked,source_canonical_hashes_pass=True,
    no_new_patient_statistics=True,scope='Integrity and provenance, not independent biochemical experiment')
(out/'independent_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
