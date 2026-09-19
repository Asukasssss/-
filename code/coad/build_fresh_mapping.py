"""Create a COAD evidence review queue from freshly retrieved KEGG records.

All edges remain provisional. No existing mapping, RNA data or candidate list is read.
"""
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / 'runtime/coad_fresh_kegg_20260919'
OUT = ROOT / 'reports/COAD/fresh_mapping_v0_1'

SPECIAL = {
    'glutamate': ('HOLD_ID_CONFLICT', ['C00025', 'C00217'], 'KEGG C00025 为 L-Glutamate；原 HMDB03339 为 D-Glutamic acid。回查作者注释前保留两种假设，不确定为 L 型。'),
    'ribulose/xylulose 5-phosphate': ('HOLD_COMBINED_FEATURE', ['C00199', 'C00231'], '原名合并 ribulose/xylulose，原 KEGG 仅 ribulose。两个假设不拆成两条独立测量或检验。'),
    'cysteine-glutathione disulfide': ('HOLD_REACTION_ID_IN_COMPOUND_FIELD', ['C05526'], '原 R00900 是反应编号，其方程中的 C05526 为拟议化合物标识；不修改原键，需回查源注释。'),
    'Mucate': ('HOLD_ORIGINAL_ENTRY_NOT_RETURNED', ['C00879'], '原 C01807 在本次 KEGG get 批量响应中未返回；名称检索得到 galactarate/mucic acid C00879，仅作替代假设。'),
    'ophthalmate': ('PROVISIONAL_NEW_KEGG_ID', ['C21016'], '原表仅 HMDB；本次按 ophthalmate 独立名称检索得到 C21016。尚未完成原 HMDB 交叉核验。'),
    'Disulfiram': ('HOLD_EXPOSURE_OR_ANNOTATION_REVIEW', ['C01692'], '名称与 KEGG 一致，但应核对外源暴露或原始注释。药物抑制靶点不等于直接生成/消耗该测量物的酶。'),
    'Acetohydroxamate': ('HOLD_EXPOSURE_OR_ANNOTATION_REVIEW', ['C06808'], '名称可对应 acetohydroxamic acid；先核对外源背景/原注释，不将药理靶点当作代谢关系。'),
    'triethanolamine': ('HOLD_EXPOSURE_OR_ANNOTATION_REVIEW', ['C06771'], '名称与 KEGG 一致，原注释/暴露来源需核查；不自动当作内源代谢物映射。'),
    'Phthalate': ('HOLD_EXPOSURE_OR_ANNOTATION_REVIEW', ['C01606'], '名称与 KEGG 一致，原注释/暴露来源需核查；不自动当作内源代谢物映射。'),
    '2-Deoxyglucose 6-phosphate': ('HOLD_EXPOSURE_OR_ANNOTATION_REVIEW', ['C06369'], '名称与 KEGG 的 2-deoxy-D-glucose 6-phosphate 相容，需核查原注释、测量与暴露背景。'),
}


def main():
    source = ROOT / 'reference/camp/cancer_effects.tsv'
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    evidence_file = CACHE / 'evidence.json'
    ev = json.loads(evidence_file.read_text(encoding='utf-8'))
    assert digest == ev['source_sha256']
    effects = [r for r in csv.DictReader(source.open(encoding='utf-8-sig'), delimiter='\t')
               if r['cancer'] == 'COAD' and float(r['effect_fdr']) < .05]
    identities, rows = [], []
    seen = set()
    for effect in effects:
        feature, original = effect['feature_name'], effect['kegg_id']
        status, hypotheses, note = SPECIAL.get(feature, (
            'ORIGINAL_KEGG_NAME_COMPATIBLE', [original],
            '已对照本次 KEGG 名称/同义名；只支持注释相容，不是重新鉴定。原 HMDB 未完成全面交叉核验。'))
        if feature in ['3-hydroxyproline', '5-hydroxylysine']:
            note += ' 只能匹配游离化合物参与的反应，不由名称跳接胶原蛋白残基羟化酶。'
        if feature == 'allantoin':
            note += ' 原名未标手性，KEGG 为 (S)-allantoin；无对应人类反应也不表示不存在非酶生成。'
        identity = {**effect, 'identity_review_status': status, 'compound_hypotheses': hypotheses,
                    'database_names': {c: ev['compound_records'].get(c, {}).get('NAME', []) for c in hypotheses},
                    'identity_note': note, 'hmdb_crosscheck': 'NOT_COMPLETED',
                    'source_urls': ['https://www.kegg.jp/entry/' + c for c in hypotheses]}
        if feature == 'glutamate':
            identity['hmdb_crosscheck'] = 'CONFLICT_CONFIRMED'
            identity['source_urls'].append('https://hmdb.ca/metabolites/HMDB0003339')
        elif feature == '3-hydroxyproline':
            identity['hmdb_crosscheck'] = 'NAME_READ_NOT_STEREOCHEMISTRY_RESOLVED'
            identity['source_urls'].append('https://hmdb.ca/metabolites/HMDB0002113')
        elif feature == 'Acetohydroxamate':
            identity['hmdb_crosscheck'] = 'NAME_COMPATIBLE'
            identity['source_urls'].append('https://hmdb.ca/metabolites/HMDB0014691')
        for compound in hypotheses:
            cr = ev['compound_records'].get(compound, {})
            rids = re.findall(r'R\d{5}', ' '.join(cr.get('REACTION', [])))
            for rid in rids:
                reaction = ev['reaction_records'].get(rid)
                if not reaction:
                    continue
                equation = ' '.join(reaction.get('EQUATION', []))
                sides = re.split(r'\s*(?:<=>|=>|<=)\s*', equation)
                if len(sides) != 2:
                    raise ValueError((rid, equation))
                tokens = [set(re.findall(r'\bC\d{5}\b', side)) for side in sides]
                if not any(compound in t for t in tokens):
                    continue
                role = 'BOTH' if all(compound in t for t in tokens) else 'LEFT' if compound in tokens[0] else 'RIGHT'
                annotated_kos = set(re.findall(r'K\d{5}', ' '.join(reaction.get('ORTHOLOGY', []))))
                for ko in ev['reaction_human_kos'].get(rid, []):
                    assert ko.split(':')[1] in annotated_kos, (rid, ko)
                    for gene in ev['human_ko_genes'][ko]:
                        key = (effect['dataset'], feature, effect['metabolite_key'], compound, rid, ko, gene)
                        if key in seen:
                            continue
                        seen.add(key)
                        name = ev['human_gene_names'].get(gene, '')
                        # KEGG list/hsa uses symbols before ';'; missing symbols stay explicit.
                        symbol = name.split(';')[0].split(',')[0].strip() if ';' in name else ''
                        rows.append({'dataset': effect['dataset'], 'feature_name': feature,
                                     'metabolite_key': effect['metabolite_key'], 'compound_hypothesis': compound,
                                     'identity_review_status': status, 'reaction_id': rid,
                                     'equation': equation, 'written_equation_side': role,
                                     'reaction_definition': ' '.join(reaction.get('DEFINITION', [])),
                                     'ko_id': ko, 'human_gene_id': gene, 'gene_symbol': symbol,
                                     'gene_description': name,
                                     'edge_status': 'PROVISIONAL_REACTION_KO_HUMAN_CHAIN',
                                     'biochemical_specificity_review': 'NOT_DONE_PER_GENE',
                                     'source_urls': ['https://www.kegg.jp/entry/' + x for x in [compound, rid, ko, gene]]})
        identities.append(identity)
    grouped = defaultdict(list)
    for r in rows:
        grouped[(r['dataset'], r['feature_name'], r['metabolite_key'], r['human_gene_id'])].append(r)
    pairs = []
    for key, support in sorted(grouped.items()):
        first = support[0]
        pairs.append({k: first[k] for k in ['dataset', 'feature_name', 'metabolite_key', 'human_gene_id', 'gene_symbol', 'identity_review_status']} | {
            'compound_hypotheses': sorted({r['compound_hypothesis'] for r in support}),
            'reaction_ids': sorted({r['reaction_id'] for r in support}),
            'ko_ids': sorted({r['ko_id'] for r in support}), 'supporting_records': len(support),
            'mapping_version': 'COAD_FRESH_KEGG_v0.1',
            'pair_status': 'HOLD_IDENTITY' if first['identity_review_status'].startswith('HOLD_') else 'PROVISIONAL_BIOCHEMISTRY_REVIEW_REQUIRED',
            'ready_for_patient_testing': False})
    for ident in identities:
        matched = [p for p in pairs if p['metabolite_key'] == ident['metabolite_key'] and p['feature_name'] == ident['feature_name']]
        ident['provisional_pair_count'] = len(matched)
        ident['provisional_genes'] = sorted({p['gene_symbol'] or p['human_gene_id'] for p in matched})
        ident['reaction_search_status'] = 'CANDIDATE_CHAIN_FOUND' if matched else 'NO_EXACT_COMPOUND_REACTION_HUMAN_KO_CHAIN_IN_THIS_SEARCH'
        ident['no_chain_is_negative_evidence'] = False
        ident['transporter_search_status'] = 'NOT_RUN'
    summary = {
        'version': 'COAD_FRESH_KEGG_v0.1', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'DONE_FRESH_REACTION_CANDIDATE_DISCOVERY_NOT_VALIDATED_MAPPING',
        'source_sha256': digest, 'evidence_sha256': hashlib.sha256(evidence_file.read_bytes()).hexdigest(),
        'base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'effects': len(effects), 'identity_status_counts': dict(Counter(i['identity_review_status'] for i in identities)),
        'features_with_candidate_chains': sum(i['provisional_pair_count'] > 0 for i in identities),
        'features_without_candidate_chains': sum(i['provisional_pair_count'] == 0 for i in identities),
        'reaction_support_records': len(rows), 'unique_feature_gene_pairs': len(pairs),
        'unique_human_gene_ids': len({p['human_gene_id'] for p in pairs}),
        'pair_status_counts': dict(Counter(p['pair_status'] for p in pairs)),
        'validated_final_pairs': 0, 'patient_analysis': 'NOT_RUN',
        'server_connection': 'ACCESS_BLOCKED: local SSH alias server165 does not resolve',
        'original_stats_changed': False, 'old_mappings_used': False,
        'transporters': 'NOT_RUN', 'functional_validation': 'NOT_RUN',
        'hmdb_full_crosscheck': 'NOT_COMPLETED; targeted pages partly returned HTTP403',
    }
    assert len(identities) == 73 and len({i['metabolite_key'] for i in identities}) == 73
    assert not any(p['ready_for_patient_testing'] for p in pairs)
    assert len(pairs) == len({(p['dataset'], p['feature_name'], p['metabolite_key'], p['human_gene_id']) for p in pairs})
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    OUT.mkdir(parents=True, exist_ok=True)
    for filename, data in [('identity_review.json', identities), ('reaction_evidence.json', rows),
                           ('provisional_feature_gene_pairs.json', pairs), ('summary.json', summary),
                           ('source_registry.json', json.loads((CACHE/'sources.json').read_text(encoding='utf-8')))]:
        (OUT/filename).write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    report = ['# COAD 从效应记录独立建立的身份与反应候选初稿', '',
              '只以冻结表中 73 条 COAD 原 q<0.05 记录为入口，重新读取 KEGG 化合物、反应、KO 和人类基因注释。未读取旧基因映射或 BRCA 候选。', '',
              f"本次 {summary['features_with_candidate_chains']} 个特征找到候选反应链，{summary['features_without_candidate_chains']} 个未找到；共 {len(pairs)} 个唯一特征—人类基因候选组合、{summary['unique_human_gene_ids']} 个基因，底层有 {len(rows)} 条反应支持记录。", '',
              '**这些数值是数据库候选检索的覆盖，不是经逐基因底物核实的正式映射，更不是显著相关或功能验证。当前进入患者统计的正式锁定关系数为 0。**', '',
              '采用化合物在反应方程中明确出现 → 反应 KO 注释 → KEGG 人类基因的链路。没有把同一通路的全部基因连进来，也没有仅凭一个宽泛 EC 编号批量展开。方程左右侧仅按书写记录，不推断组织内反应方向。多来源/多反应归并到同一个特征—基因组合；合并峰不拆成独立测量。', '',
              '# 身份问题与边界', '',
              '- glutamate：KEGG C00025（L 型）与 HMDB03339（D 型）冲突，两种身份假设保留并暂挂。',
              '- ribulose/xylulose 5-phosphate：合并名称，C00199/C00231 只作为两个假设，不能把一个测量拆成两个。',
              '- cysteine-glutathione disulfide：原 KEGG 字段为 R00900 反应编号；从方程提出 C05526，尚未替换原注释。',
              '- Mucate：本次 KEGG 批量响应未返回原 C01807；名称检索提出 C00879，尚待作者注释核对。',
              '- ophthalmate：本次名称检索新提出 C21016，保留原 HMDB 键，交叉核验未完成。',
              '- Disulfiram、Acetohydroxamate、triethanolamine、Phthalate、2-Deoxyglucose 6-phosphate 暂列原注释/暴露背景核查，不由药物靶点倒推直接代谢关系。',
              '- HIGH_ID 不能消除上述问题。其余名称与 KEGG 注释相容也不等于重新鉴定；HMDB 尚未全面核查，立体化学需结合作者方法。', '',
              '# 全部 73 条的去向', '', '| 原特征 | 身份状态 | 候选基因数 | 候选基因（仅展示前 12 个） |', '|---|---|---:|---|']
    for ident in identities:
        genes = ident['provisional_genes']
        report.append(f"| {ident['feature_name']} | {ident['identity_review_status']} | {len(genes)} | {', '.join(genes[:12]) + (' …' if len(genes)>12 else '')} |")
    report += ['', '# 完成项、未完成项与下一步', '',
               'DONE：73 条的数据库名称/标识初查、KEGG 反应链检索、按实际特征与人类基因去重、出处和访问时间/哈希记录。',
               'NOT_RUN：系统转运体检索、逐人类蛋白底物与复合体角色核对、RNA/患者关联、依赖/功能分析、外部复现。无候选链不是阴性，不能淘汰代谢物。',
               'ACCESS_BLOCKED：本机未配置可解析的 server165 SSH 别名；部分 HMDB 页面/XML 返回 HTTP403。未下载任何患者级矩阵。',
               '下一步：回查有冲突的作者注释；从反应候选中核对人类蛋白的直接催化或复合体角色，并补转运体。身份与底物证据足够后才锁定患者分析关系和 BH 家族。', '',
               '来源：[KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html)、[L-glutamate](https://www.kegg.jp/entry/C00025)、[D-glutamate HMDB](https://hmdb.ca/metabolites/HMDB0003339)、[合并名称的 ribulose 对应条目](https://www.kegg.jp/entry/C00199)。逐条链接见 identity_review.json 和 reaction_evidence.json；完整候选组合见 provisional_feature_gene_pairs.json。',
               '', '机器校验：源哈希不变；73 条全部有去向；每条链的目标化合物确实在方程中、KO 同时存在于反应记录和人类映射；组合无重复。此校验不代替生物学底物审查。']
    (OUT/'RESULTS_CN.md').write_text('\n'.join(report)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
