"""Scan reviewed human transport proteins for explicit transport of target compounds.

Require exact ChEBI (or official pH mapping) and same named molecule with (in)/(out).
No free-text mention or regulator alone can produce an edge.
"""
import csv
import hashlib
import io
import json
import re
import urllib.parse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import check_human_reactions as api

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'reports/COAD/transport_check_v0_1'


def main():
    identities = json.loads((ROOT/'reports/COAD/fresh_mapping_v0_1/identity_review.json').read_text(encoding='utf-8'))
    ev = json.loads((ROOT/'runtime/coad_fresh_kegg_20260919/evidence.json').read_text(encoding='utf-8'))
    ph = {r['CHEBI']: r for r in csv.DictReader(io.StringIO(api.fetch('https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv')), delimiter='\t')}
    names = {}
    for line in api.fetch('https://ftp.expasy.org/databases/rhea/tsv/chebiId_name.tsv').splitlines():
        ident, name = line.split('\t', 1)
        names[ident.removeprefix('CHEBI:')] = name.strip()
    compound_chebi = {}
    for compound, record in ev['compound_records'].items():
        ids = set()
        for line in record.get('DBLINKS', []):
            if line.startswith('ChEBI:'):
                ids.update(re.findall(r'\d+', line))
        compound_chebi[compound] = {i: ph.get(i, {}).get('CHEBI_PH7_3', i) for i in ids}
    keyword = json.loads(api.fetch('https://rest.uniprot.org/keywords/KW-0813'))
    assert keyword['keyword']['name'] == 'Transport'
    params = {'query': 'organism_id:9606 AND reviewed:true AND keyword:KW-0813', 'format': 'json',
              'fields': 'accession,gene_names,organism_id,xref_geneid,cc_catalytic_activity,cc_subunit'}
    url = 'https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode(params)
    proteins = json.loads(api.fetch(url))['results']
    assert any(p['primaryAccession']=='P31641' for p in proteins), 'Known taurine transporter missing from query'
    print(f'Read {len(proteins)} reviewed human transport-keyword records', flush=True)
    evidence = []
    for ident in identities:
        for compound in ident['compound_hypotheses']:
            cid_pairs = compound_chebi.get(compound, {})
            allowed = set(cid_pairs) | set(cid_pairs.values())
            if not allowed:
                continue
            for protein in proteins:
                assert protein['organism']['taxonId'] == 9606 and protein['entryType']=='UniProtKB reviewed (Swiss-Prot)'
                gene_ids = [x['id'] for x in protein.get('uniProtKBCrossReferences', []) if x['database']=='GeneID']
                symbols = [g['geneName']['value'] for g in protein.get('genes', []) if 'geneName' in g]
                if not gene_ids:
                    continue
                for comment in protein.get('comments', []):
                    if comment['commentType'] != 'CATALYTIC ACTIVITY':
                        continue
                    reaction = comment['reaction']
                    text = reaction['name']
                    if '(in)' not in text or '(out)' not in text:
                        continue
                    chebis = {x['id'].removeprefix('CHEBI:') for x in reaction.get('reactionCrossReferences',[]) if x['database']=='ChEBI'}
                    matches = allowed & chebis
                    # Only substrate names attached to compartment labels count, not cofactors mentioned elsewhere.
                    matched_names = {c: names[c] for c in matches if c in names and names[c]+'(in)' in text and names[c]+'(out)' in text}
                    if not matched_names:
                        continue
                    for gene_id in gene_ids:
                        evidence.append({'dataset': ident['dataset'], 'feature_name': ident['feature_name'],
                                         'metabolite_key': ident['metabolite_key'], 'compound_hypothesis': compound,
                                         'identity_review_status': ident['identity_review_status'],
                                         'human_gene_id': 'hsa:'+gene_id, 'gene_symbol': ';'.join(symbols),
                                         'uniprot_accession': protein['primaryAccession'], 'reaction': text,
                                         'rhea_ids': [x['id'] for x in reaction.get('reactionCrossReferences',[]) if x['database']=='Rhea'],
                                         'matched_chebi_names': matched_names,
                                         'kegg_chebi_to_ph_mapping': cid_pairs,
                                         'evidence': reaction.get('evidences', []),
                                         'has_experimental_annotation': any(x['evidenceCode']=='ECO:0000269' for x in reaction.get('evidences', [])),
                                         'subunit_notes': [c for c in protein.get('comments',[]) if c['commentType']=='SUBUNIT'],
                                         'source_url': 'https://www.uniprot.org/uniprotkb/'+protein['primaryAccession']+'/entry',
                                         'no_current_identity_hold': not ident['identity_review_status'].startswith('HOLD_'),
                                         'relationship': 'DIRECT_TRANSPORT_REACTION_ANNOTATION',
                                         'coad_transport_or_function': 'NOT_ASSESSED'})
    groups = defaultdict(list)
    for row in evidence:
        groups[(row['dataset'], row['feature_name'], row['metabolite_key'], row['human_gene_id'])].append(row)
    pairs = []
    for key, support in sorted(groups.items()):
        first = support[0]
        pairs.append({k: first[k] for k in ['dataset','feature_name','metabolite_key','human_gene_id','gene_symbol','identity_review_status','no_current_identity_hold']} | {
            'has_experimental_annotation': any(r['has_experimental_annotation'] for r in support),
            'compound_hypotheses': sorted({r['compound_hypothesis'] for r in support}),
            'uniprot_accessions': sorted({r['uniprot_accession'] for r in support}),
            'rhea_ids': sorted({r for s in support for r in s['rhea_ids']}),
            'ready_for_patient_testing': False, 'relationship': 'DIRECT_TRANSPORT_REACTION_ANNOTATION'})
    summary = {'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'DONE_SCOPED_TRANSPORT_REACTION_SCAN',
               'reviewed_human_transport_keyword_entries': len(proteins), 'features_searched':len(identities),
               'evidence_records': len(evidence), 'unique_feature_gene_pairs':len(pairs),
               'features_with_pairs':len({r['metabolite_key'] for r in pairs}),
               'unique_gene_ids':len({r['human_gene_id'] for r in pairs}),
               'pairs_without_current_identity_hold':sum(r['no_current_identity_hold'] for r in pairs),
               'pairs_with_experimental_annotation':sum(r['has_experimental_annotation'] for r in pairs),
               'features_without_pairs':len(identities)-len({r['metabolite_key'] for r in pairs}),
               'scope':'Reviewed human UniProt transport keyword; catalytic reaction explicitly names same ChEBI substrate on both in/out sides; direct or official Rhea pH microspecies mapping. No keyword-only association.',
               'limitations':'No hit is not negative. Transporters lacking explicit reaction annotations, alternative compartment wording, unreviewed proteins and more general substrate classes are outside this scoped scan. Original identity remains provisional. Experimental annotation is not COAD functional evidence.',
               'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    assert any(p['feature_name']=='taurine' and p['gene_symbol']=='SLC6A6' for p in pairs), 'Positive transport control failed'
    OUT.mkdir(parents=True,exist_ok=True)
    for name, data in [('transport_evidence.json',evidence),('transport_pairs.json',pairs),('summary.json',summary),('source_registry.json',api.sources)]:
        (OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    lines=['# COAD 直接转运反应注释检索', '',
           f"检索 {len(identities)} 个显著特征，读取 {len(proteins)} 个带 transport 关键词的 reviewed 人类蛋白条目。得到 {len(pairs)} 条唯一特征—基因组合，覆盖 {summary['features_with_pairs']} 个特征、{summary['unique_gene_ids']} 个基因。", '',
           '只纳入 ChEBI 精确匹配（或 Rhea 官方 pH7.3 映射）、且同一底物名称出现在反应 in/out 两侧的转运注释。关键词命中、抑制运输、调节蛋白或同一通路不足以纳入。', '',
           '| 原特征 | 转运基因 | 身份暂挂 |', '|---|---|---|']
    for name in sorted({r['feature_name'] for r in pairs}):
        subset=[r for r in pairs if r['feature_name']==name]
        lines.append(f"| {name} | {', '.join(sorted({r['gene_symbol'] for r in subset}))} | {'是' if any(not r['no_current_identity_hold'] for r in subset) else '未发现特定暂挂；身份仍待源注释核查'} |")
    lines+=['', '范围限制：无结果不说明没有转运体；只覆盖上述结构化反应、reviewed 条目和 in/out 命名格式，尚未穷尽自由文本、其他区室表述和底物类别注释。保留每条实验/推断证据及复合体说明，不把辅助亚基自动理解成独立转运载体。',
            '', '现阶段仅得到生化注释支持，未计算 COAD 患者相关，未评价功能或确定治疗靶点。', '',
            '来源：[UniProt API](https://www.uniprot.org/help/api)、[Rhea 下载与化学实体映射](https://www.rhea-db.org/help/download)。逐条出处在 transport_evidence.json。']
    (OUT/'RESULTS_CN.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
