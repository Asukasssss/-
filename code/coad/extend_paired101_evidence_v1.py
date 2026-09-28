"""Extend immutable COAD biochemical evidence to the paired-P work pool.

Reuses archived public database evidence; new requests cached separately. No patient data.
"""
import csv,hashlib,json,re,sys,urllib.parse
from pathlib import Path
import fetch_fresh_kegg as kegg
import check_human_reactions as human
import check_transport_reactions as transport
from build_fresh_mapping import SPECIAL

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'runtime/coad_paired101_mapping_20260922'
WORK.mkdir(exist_ok=True)
OLD=ROOT/'runtime/coad_fresh_kegg_20260919'


def immutable_fallback(original,oldcache,registry):
    def get(url):
        item=registry.get(url)
        if item and item.get('status')=='DONE':
            p=oldcache/item['cache_file']
            assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
            return p.read_text(encoding='utf-8')
        return original(url)
    return get


def setup():
    kegg.CACHE=WORK/'kegg';kegg.CACHE.mkdir(exist_ok=True);kegg.REGISTRY=kegg.CACHE/'sources.json'
    kegg.sources=json.loads(kegg.REGISTRY.read_text()) if kegg.REGISTRY.exists() else {}
    human.CACHE=WORK/'uniprot';human.CACHE.mkdir(exist_ok=True);human.REG=human.CACHE/'sources.json'
    human.sources=json.loads(human.REG.read_text()) if human.REG.exists() else {}
    oldcache=ROOT/'runtime/coad_uniprot_rhea_20260919'
    human.fetch=immutable_fallback(human.fetch,oldcache,json.loads((oldcache/'sources.json').read_text()))
    return oldcache


def acquire():
    setup()
    with (ROOT/'results/COAD/03_PATIENT/20260921T123200Z_paired_metabolites33_v1/primary33.tsv').open() as f:
        selected={r['metabolite_name'] for r in csv.DictReader(f,delimiter='\t') if float(r['p_value'])<.05}
    with (ROOT/'reference/camp/cancer_effects.tsv').open() as f:
        effects=[r for r in csv.DictReader(f,delimiter='\t') if r['cancer']=='COAD' and r['feature_name'] in selected]
    assert len(effects)==101
    ev=json.loads((OLD/'evidence.json').read_text())
    oldident=json.loads((ROOT/'reports/COAD/fresh_mapping_v0_1/identity_review.json').read_text())
    oldnames={i['feature_name'] for i in oldident}
    new=[r for r in effects if r['feature_name'] not in oldnames];assert len(new)==29
    compounds=kegg.get_entries([r['kegg_id'] for r in new if re.fullmatch(r'C\d{5}',r['kegg_id']) and r['kegg_id'] not in ev['compound_records']])
    ev['compound_records'].update(compounds)
    # Names lacking a usable identifier are recorded, not silently assigned a near match.
    for r in new:
        if not r['kegg_id']:
            ev.setdefault('new_name_queries',{})[r['feature_name']]=kegg.fetch('find/compound/'+urllib.parse.quote(r['feature_name']))
    allrk={}
    # Reuse the complete archived reaction/KO link response, whose hash is checked.
    reg=json.loads((OLD/'sources.json').read_text());item=reg['https://rest.kegg.jp/link/ko/reaction'];p=OLD/item['cache_file']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
    for line in p.read_text().splitlines():
        reaction,ko=line.split('\t')
        if ko in ev['human_ko_genes']:allrk.setdefault(reaction.removeprefix('rn:'),[]).append(ko)
    rids=set()
    for c in compounds.values():rids.update(re.findall(r'R\d{5}',' '.join(c.get('REACTION',[]))))
    needed=sorted(rids & set(allrk)-set(ev['reaction_records']))
    ev['reaction_records'].update(kegg.get_entries(needed))
    for rid in rids & set(allrk):ev['reaction_human_kos'][rid]=allrk[rid]
    (WORK/'evidence.json').write_text(json.dumps(ev,ensure_ascii=False,indent=2),encoding='utf-8')
    (WORK/'new_features.json').write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8')
    print('NEW COMPOUND NAMES',flush=True)
    for r in new: print(r['feature_name'],r['kegg_id'],ev['compound_records'].get(r['kegg_id'],{}).get('NAME'),flush=True)
    print('NAME QUERIES',ev.get('new_name_queries',{}),flush=True)


if __name__=='__main__':acquire()
