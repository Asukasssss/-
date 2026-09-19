"""Describe frozen COAD effects only; no mappings or new statistical tests."""
import csv
import hashlib
import json
import math
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'reference/camp/cancer_effects.tsv'
OUT = ROOT / 'reports/COAD/effect_records_v1'


def main():
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    rows = list(csv.DictReader(SOURCE.open(encoding='utf-8-sig', newline=''), delimiter='\t'))
    coad = [r for r in rows if r['cancer'] == 'COAD']
    keys = [(r['dataset'], r['feature_name'], r['metabolite_key']) for r in coad]
    assert coad and len(keys) == len(set(keys)), 'Missing COAD or duplicate actual feature keys'
    for r in coad:
        for col in ['hedges_g', 'effect_fdr', 'hedges_g_se', 'max_input_raw_missing_rate']:
            assert math.isfinite(float(r[col])), (r['feature_name'], col)
        assert 0 <= float(r['effect_fdr']) <= 1
        assert 0 <= float(r['max_input_raw_missing_rate']) <= 1
    sig = [r for r in coad if float(r['effect_fdr']) < .05]
    def counts(rs, col):
        return dict(Counter(r[col] or 'UNANNOTATED' for r in rs))
    def direction(rs):
        return dict(Counter('positive' if float(r['hedges_g']) > 0 else 'negative' if float(r['hedges_g']) < 0 else 'zero' for r in rs))
    summary = {
        'status': 'DONE_DESCRIPTIVE_EFFECT_REVIEW_ONLY',
        'source': SOURCE.relative_to(ROOT).as_posix(), 'source_sha256': digest,
        'base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'analysis_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'all_records': len(coad), 'unique_metabolite_keys': len({r['metabolite_key'] for r in coad}),
        'original_q_lt_005': len(sig), 'all_direction': direction(coad), 'significant_direction': direction(sig),
        'datasets': counts(coad, 'dataset'), 'source_types': counts(coad, 'source_type'),
        'all_identity_labels': counts(coad, 'identity_confidence'),
        'significant_identity_labels': counts(sig, 'identity_confidence'),
        'significant_super_pathways': counts(sig, 'super_pathway'),
        'significant_amino_acid_case_normalized_count': sum(r['super_pathway'].casefold() == 'amino acid' for r in sig),
        'missing_rate_range': [min(float(r['max_input_raw_missing_rate']) for r in coad), max(float(r['max_input_raw_missing_rate']) for r in coad)],
        'missing_rate_gt_050_descriptive_only': sum(float(r['max_input_raw_missing_rate']) > .5 for r in coad),
        'significant_missing_rate_gt_050_descriptive_only': sum(float(r['max_input_raw_missing_rate']) > .5 for r in sig),
        'significant_without_kegg_or_hmdb': sum(not r['kegg_id'] and not r['hmdb_id'] for r in sig),
        'tests_recomputed': False, 'existing_gene_mappings_used': False,
        'patient_sample_size': 'NOT_AVAILABLE_IN_EFFECT_TABLE',
        'patient_independence': 'NOT_ASSESSED', 'direct_gene_mapping': 'NOT_RUN',
        'method': 'Filter cancer=COAD; count original effect_fdr<0.05; describe original g, annotations and missingness. No new inferential test or exclusion threshold.',
    }
    def table(rs):
        lines = ['| 特征（原名） | Hedges g | 原 q | 原最大缺失率 | 原身份标签 |', '|---|---:|---:|---:|---|']
        for r in rs:
            lines.append(f"| {r['feature_name'].replace('|', '/')} | {float(r['hedges_g']):.3f} | {float(r['effect_fdr']):.4g} | {float(r['max_input_raw_missing_rate']):.1%} | {r['identity_confidence']} |")
        return '\n'.join(lines)
    report = [
        '# COAD：从冻结效应记录开始的描述性分析',
        '',
        '账号 B。只读取癌种效应表，不读取既有候选名单或代谢物—基因映射，不重算或更改原效应与 q。',
        '',
        f"共 {len(coad)} 条记录，{summary['unique_metabolite_keys']} 个代谢物键；原 q<0.05 的有 {len(sig)} 条，其余 {len(coad)-len(sig)} 条保留为背景。",
        f"显著记录效应方向：{json.dumps(summary['significant_direction'])}。按项目既有肿瘤减正常解释，正值为肿瘤较高、负值为较低；本轮未回查源矩阵。g 是标准化差异，不是倍数。",
        f"dataset 字段：{json.dumps(summary['datasets'])}；source_type 字段：{json.dumps(summary['source_types'])}。单队列记录不构成独立队列复现。",
        '', '# 原注释与数据限制', '',
        f"原身份标签（全部）：{json.dumps(summary['all_identity_labels'])}；显著记录：{json.dumps(summary['significant_identity_labels'])}。HIGH_ID 仅为沿用标签，不是重新鉴定。",
        f"显著记录中无 KEGG/HMDB 标识：{summary['significant_without_kegg_or_hmdb']} 条。存在标识也不保证异构体、反应底物或峰身份已确定。",
        f"原 max_input_raw_missing_rate 范围为 {summary['missing_rate_range']}。其中 >50% 的共有 {summary['missing_rate_gt_050_descriptive_only']} 条，显著记录中 {summary['significant_missing_rate_gt_050_descriptive_only']} 条。50% 仅作描述性计数，未用于剔除或改阈值；该字段不能代替实际样本量。",
        '表中没有样本数、患者身份、配对关系或原始 P 值，不能由本表核实设计、重新校正 q 或建立患者关联。',
        '', '# 通路注释分布（计数，不是富集分析）', '',
        '| 原大类 | 全部记录 | 原 q<0.05 | 正效应显著 | 负效应显著 |', '|---|---:|---:|---:|---:|',
    ]
    for pathway in sorted({r['super_pathway'] for r in coad}):
        group = [r for r in coad if r['super_pathway'] == pathway]
        hits = [r for r in sig if r['super_pathway'] == pathway]
        report.append(f"| {pathway or '未注释'} | {len(group)} | {len(hits)} | {sum(float(r['hedges_g'])>0 for r in hits)} | {sum(float(r['hedges_g'])<0 for r in hits)} |")
    report += ['', '原表 Amino Acid 与 Amino acid 分别保留显示；仅按大小写统一后的描述性合计为 ' + str(summary['significant_amino_acid_case_normalized_count']) + ' 条显著记录。此合计不是通路富集或通量结论。',
               '', '# 显著但原最大缺失率超过 50% 的记录', '', table([r for r in sig if float(r['max_input_raw_missing_rate']) > .5]),
               '', '名称 Disulfiram 的记录原 q 很小，但大类/子类注释均为空，下一步需回查原注释与身份来源。本轮不更改其名称或标识，也不据此给出基因或干预结论。',
               '', '# 显著记录：按原 q 排序', '', f'此排序用于阅读，不是靶点优先级。所有 {len(sig)} 条均保留，不用 g 或缺失率另设筛选门槛。', '', table(sorted(sig, key=lambda r: float(r['effect_fdr']))),
               '', '# 未达到原 q<0.05 的背景记录', '', table([r for r in coad if float(r['effect_fdr']) >= .05]),
               '', '# 下一步与完成边界', '',
               '效应记录描述已完成。下一步从这批显著特征逐一核对化合物身份和直接生化关系，重新建立有出处的 COAD 映射。现有基因映射未使用，RNA/样本分析未运行，功能与外部验证未运行。',
               '', f"来源：`{summary['source']}`，SHA-256：`{digest}`。可复现脚本：`code/coad/review_effect_records.py`。完整参数与状态见同目录 `summary.json`。"]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    (OUT / 'RESULTS_CN.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == digest
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print('\nLOWEST_Q_RECORDS\n' + table(sorted(sig, key=lambda r: float(r['effect_fdr']))[:10]))


if __name__ == '__main__':
    main()
