"""Check candidate enzyme links against reviewed human UniProt catalytic Rhea annotations.

An exact reaction annotation is biochemical database support, not a COAD experiment.
"""
import csv
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'reports/COAD/fresh_mapping_v0_1'
OUT = ROOT/'reports/COAD/human_reaction_check_v0_1'
CACHE = ROOT/'runtime/coad_uniprot_rhea_20260919'
CACHE.mkdir(parents=True, exist_ok=True)
REG = CACHE/'sources.json'
sources = json.loads(REG.read_text(encoding='utf-8')) if REG.exists() else {}


def fetch(url):
    p = CACHE/(hashlib.sha256(url.encode()).hexdigest()+'.txt')
    if p.exists() and url in sources:
        assert hashlib.sha256(p.read_bytes()).hexdigest() == sources[url]['sha256']
        return p.read_text(encoding='utf-8')
    for attempt in range(3):
        try:
            time.sleep(.4)
            with urllib.request.urlopen(url, timeout=40) as response:
                data = response.read()
                next_link = response.headers.get('Link', '')
                if 'rel="next"' in next_link:
                    raise RuntimeError('Unexpected pagination; reduce batch size')
            p.write_bytes(data)
            sources[url] = {'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                            'sha256': hashlib.sha256(data).hexdigest(), 'cache_file': p.name, 'status': 'DONE'}
            REG.write_text(json.dumps(sources, indent=2), encoding='utf-8')
            return data.decode()
        except Exception as exc:
            if attempt == 2:
                sources[url] = {'url': url, 'status': 'ACCESS_BLOCKED', 'error': str(exc)}
                REG.write_text(json.dumps(sources, indent=2), encoding='utf-8')
                raise
            time.sleep(1+attempt)


def main():
    pairs = json.loads((BASE/'provisional_feature_gene_pairs.json').read_text(encoding='utf-8'))
    support = json.loads((BASE/'reaction_evidence.json').read_text(encoding='utf-8'))
    ev = json.loads((ROOT/'runtime/coad_fresh_kegg_20260919/evidence.json').read_text(encoding='utf-8'))
    rhea_master = {}
    for r in csv.DictReader(io.StringIO(fetch('https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv')), delimiter='\t'):
        for ident in r.values():
            rhea_master[ident] = r['RHEA_ID_MASTER']
    # Current Rhea to KEGG links supplement KEGG's own cross-references.
    kegg_rhea = defaultdict(set)
    crossref = fetch('https://ftp.expasy.org/databases/rhea/tsv/rhea2kegg_reaction.tsv')
    for r in csv.DictReader(io.StringIO(crossref), delimiter='\t'):
        kegg_rhea[r['ID']].add(r['MASTER_ID'])
    symbols = sorted({p['gene_symbol'] for p in pairs if p['gene_symbol']})
    proteins = {}
    for start in range(0, len(symbols), 25):
        batch = symbols[start:start+25]
        query = 'organism_id:9606 AND reviewed:true AND ('+' OR '.join('gene_exact:'+s for s in batch)+')'
        params = {'query': query, 'format': 'json', 'size': 500,
                  'fields': 'accession,gene_names,organism_id,xref_geneid,cc_catalytic_activity,cc_subunit'}
        data = json.loads(fetch('https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode(params)))
        for protein in data['results']:
            proteins[protein['primaryAccession']] = protein
        print(f'Human UniProt search {min(start+25,len(symbols))}/{len(symbols)} gene symbols', flush=True)
    (CACHE/'proteins.json').write_text(json.dumps(proteins, indent=2), encoding='utf-8')
    by_gene = defaultdict(list)
    for protein in proteins.values():
        assert protein['organism']['taxonId'] == 9606 and protein['entryType'] == 'UniProtKB reviewed (Swiss-Prot)'
        for ref in protein.get('uniProtKBCrossReferences', []):
            if ref['database'] == 'GeneID':
                by_gene['hsa:'+ref['id']].append(protein)
    kegg_reactions = ev['reaction_records']
    for rid, record in kegg_reactions.items():
        for line in record.get('DBLINKS', []):
            if line.startswith('RHEA:'):
                for ident in re.findall(r'\d+', line):
                    if ident in rhea_master:
                        kegg_rhea[rid].add(rhea_master[ident])
    records_by_pair = defaultdict(list)
    for r in support:
        records_by_pair[(r['metabolite_key'], r['feature_name'], r['human_gene_id'])].append(r)
    checked = []
    for pair in pairs:
        matches = {}
        for protein in by_gene[pair['human_gene_id']]:
            for comment in protein.get('comments', []):
                if comment['commentType'] != 'CATALYTIC ACTIVITY':
                    continue
                reaction = comment['reaction']
                rhea_ids = [x['id'].removeprefix('RHEA:') for x in reaction.get('reactionCrossReferences', []) if x['database'] == 'Rhea']
                masters = {rhea_master.get(x, x) for x in rhea_ids}
                for record in records_by_pair[(pair['metabolite_key'], pair['feature_name'], pair['human_gene_id'])]:
                    common = masters & kegg_rhea[record['reaction_id']]
                    if not common:
                        continue
                    key = (protein['primaryAccession'], record['compound_hypothesis'], record['reaction_id'], tuple(sorted(common)))
                    evidences = reaction.get('evidences', [])
                    matches[key] = {
                        'uniprot_accession': protein['primaryAccession'], 'compound_hypothesis': record['compound_hypothesis'],
                        'kegg_reaction_id': record['reaction_id'], 'rhea_master_ids': sorted(common),
                        'uniprot_reaction': reaction['name'], 'annotation_evidence': evidences,
                        'has_experimental_annotation': any(e['evidenceCode']=='ECO:0000269' for e in evidences),
                        'subunit_notes': [x for x in protein.get('comments', []) if x['commentType']=='SUBUNIT'],
                        'source_url': 'https://www.uniprot.org/uniprotkb/'+protein['primaryAccession']+'/entry'}
        matches = list(matches.values())
        status = ('NO_REVIEWED_HUMAN_RECORD_LINKED' if not by_gene[pair['human_gene_id']]
                  else 'MATCHED_REACTION_WITH_EXPERIMENTAL_ANNOTATION' if any(x['has_experimental_annotation'] for x in matches)
                  else 'MATCHED_REACTION_OTHER_ANNOTATION' if matches else 'NO_EXACT_REACTION_MATCH_IN_CURRENT_SEARCH')
        checked.append({**pair, 'human_reaction_check': status, 'matched_support': matches,
                        'checked_uniprot_accessions': sorted({p['primaryAccession'] for p in by_gene[pair['human_gene_id']]}),
                        'coad_functional_evidence': 'NOT_ASSESSED',
                        'no_current_identity_hold': not pair['identity_review_status'].startswith('HOLD_'),
                        'original_measurement_identity_confirmed': False,
                        'ready_for_patient_testing': False})
    assert len(checked) == len(pairs)
    matched = [r for r in checked if r['matched_support']]
    clear_matched = [r for r in matched if r['no_current_identity_hold']]
    summary = {'created_utc': datetime.now(timezone.utc).isoformat(),
               'status': 'DONE_BATCH_HUMAN_REACTION_ANNOTATION_CHECK',
               'pairs_checked': len(checked), 'human_genes_searched': len(symbols),
               'reviewed_human_proteins_retrieved': len(proteins),
               'status_counts': dict(Counter(r['human_reaction_check'] for r in checked)),
               'pairs_with_matching_human_reaction_annotation': len(matched),
               'matching_pairs_without_current_identity_hold': len(clear_matched),
               'matching_pairs_with_identity_hold': len(matched)-len(clear_matched),
               'features_with_matching_annotation_without_current_identity_hold': len({r['metabolite_key'] for r in clear_matched}),
               'genes_with_matching_annotation_without_current_identity_hold': len({r['human_gene_id'] for r in clear_matched}),
               'patient_analysis': 'NOT_RUN', 'transporter_review': 'NOT_RUN',
               'final_patient_test_family_locked': False,
               'input_pairs_sha256': hashlib.sha256((BASE/'provisional_feature_gene_pairs.json').read_bytes()).hexdigest(),
               'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'method': 'Exact NCBI GeneID link to reviewed human UniProt; exact Rhea reaction match normalized by official direction table; not EC-only or symbol-only matching.',
               'limitation': 'Annotation support is not a new experimental finding, complete substrate-specificity review, independent validation or proof in COAD. Identity hypotheses remain conditional. No exact match is not a negative result.'}
    OUT.mkdir(parents=True, exist_ok=True)
    for filename, content in [('human_reaction_pairs.json',checked), ('summary.json',summary), ('source_registry.json',sources)]:
        (OUT/filename).write_text(json.dumps(content, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    report = ['# COAD 人类蛋白直接反应注释核对', '',
              f"对新检索的 {len(checked)} 条特征—基因候选逐条核对。{len(matched)} 条能与 reviewed 人类 UniProt 的催化反应通过同一 Rhea 反应对应；其中 {len(clear_matched)} 条当前没有身份暂挂，涉及 {summary['features_with_matching_annotation_without_current_identity_hold']} 个特征、{summary['genes_with_matching_annotation_without_current_identity_hold']} 个人类基因。", '',
              '这里的“没有身份暂挂”只指本次已发现的特定问题，不代表原始鉴定、全部 HMDB 交叉核验或立体化学已经解决。反应对应是已有生化知识支持，不是 COAD 功能证据；当前没有锁定患者检验家族。', '',
              '| 核对结果 | 关系数 |', '|---|---:|']
    report += [f'| {k} | {v} |' for k,v in summary['status_counts'].items()]
    report += ['', 'ECO:0000269 是 UniProt 所记录的实验注释依据，本文未逐篇复核原文或实验条件。其他注释可能包含推断/相似性依据。未找到完全相同 Rhea 反应的关系保留待审，不当作反证。', '',
               '# 按原代谢特征列出没有当前身份暂挂的匹配基因', '', '| 原特征 | 匹配基因数 | 基因 |', '|---|---:|---|']
    by_feature = defaultdict(list)
    for r in clear_matched:
        by_feature[r['feature_name']].append(r['gene_symbol'] or r['human_gene_id'])
    for name, genes in sorted(by_feature.items()):
        report.append(f"| {name} | {len(set(genes))} | {', '.join(sorted(set(genes)))} |")
    report += ['', '# 方法与交付边界', '',
               '从每条 KEGG 反应的 Rhea 交叉引用出发，利用 Rhea 官方方向表统一同一反应的方向编号；只与 NCBI GeneID 精确对应的 reviewed 人类 UniProt 催化反应相连。原反应中必须明确包含待查化合物，不能靠通路、相似基因名或宽泛 EC 编号连接。保留复合体注释，不能把所有亚基都解释为独立催化酶。',
               '此版本仅覆盖第一轮 KEGG 反应候选内的关系；系统转运体检索、KO链未捕获的其他酶关系、全文功能审查、患者分析尚未完成。',
               '逐条反应、UniProt 链接、证据代码和 PubMed 编号均保存在 human_reaction_pairs.json。可公开同步的是这些数据库来源与汇总，不含患者矩阵。',
               '', '来源：[Rhea 官方下载说明](https://www.rhea-db.org/help/download)、[UniProt](https://www.uniprot.org/)。']
    (OUT/'RESULTS_CN.md').write_text('\n'.join(report)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
