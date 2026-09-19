"""Export existing COAD evidence to the main-branch seven-stage delivery standard.

Format conversion only; no mapping discovery, statistical analysis or q recalculation.
"""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGES = ['01_CAMP', '02_MAPPING', '03_PATIENT', '04_ROBUSTNESS',
          '05_FUNCTION', '06_EXTERNAL', '07_INTEGRATION']
HEADINGS = ['本轮问题', '输入与范围', '实际结果', '新手解释',
            '限制/反证', '当前决定', '下一步', '复现命令']
NA = 'NA'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))


def read_tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def header(name):
    return (ROOT / 'templates' / name).read_text(encoding='utf-8').strip().split('\t')


def write_tsv(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8', newline='\n')


def validate_stage(path, stage, frozen):
    records = read_tsv(path / 'results.tsv')
    fields = list(records[0])
    assert fields[:len(header('statistical_result.tsv'))] == header('statistical_result.tsv')
    keys = [(r['cohort'], r['metabolite_name'], r['metabolite_key'], r['gene'], r['human_gene_id']) for r in records]
    assert keys == sorted(keys) and len(keys) == len(set(keys))
    for row in records:
        source = frozen[(row['cohort'], row['metabolite_name'], row['metabolite_key'])]
        assert (row['original_effect'], row['original_q']) == (source['hedges_g'], source['effect_fdr'])
        assert row['stage_id'] == stage and row['cancer'] == 'COAD'
        assert row['status'] in {'DONE', 'NEEDS_REVIEW', 'NOT_EVALUABLE'}
        if stage == '01_CAMP':
            assert (row['effect'], row['q_value']) == (source['hedges_g'], source['effect_fdr'])
        else:
            assert row['effect'] == row['q_value'] == row['p_value'] == NA
    for entry in read_tsv(path / 'source_manifest.tsv'):
        assert sha(ROOT / entry['source_path']) == entry['sha256']
    headings = [line[3:] for line in (path / 'README_CN.md').read_text(encoding='utf-8').splitlines() if line.startswith('## ')]
    assert headings == HEADINGS
    assert len(records) == (159 if stage == '01_CAMP' else 974)
    if (path / 'validation.json').exists():
        validation = json.loads((path / 'validation.json').read_text(encoding='utf-8'))
        for name, expected in validation['artifact_sha256'].items():
            assert sha(path / name) == expected
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--stage', choices=STAGES[:2])
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    if not args.run_id or any(c not in '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-' for c in args.run_id):
        parser.error('Use a UTC timestamp and task name; no path separators.')
    source = 'reference/camp/cancer_effects.tsv'
    effects = [r for r in read_tsv(ROOT / source) if r['cancer'] == 'COAD']
    frozen = {(r['dataset'], r['feature_name'], r['metabolite_key']): r for r in effects}
    catalog_path = 'reports/COAD/current_catalog_v0_1/candidate_catalog.json'
    identity_path = 'reports/COAD/fresh_mapping_v0_1/identity_review.json'
    candidates = read_json(catalog_path)
    stages = [args.stage] if args.stage else STAGES[:2]
    for stage in stages:
        output = ROOT / 'results/COAD' / stage / args.run_id
        if args.validate_only:
            validate_stage(output, stage, frozen)
            print('PASS:', stage, 'counts, shared fields, order, uniqueness, frozen values, source hashes, report order')
            continue
        output.mkdir(parents=True, exist_ok=False)
        rows = []
        extras = ['human_gene_id', 'original_effect', 'original_q', 'identity_review_status', 'detail_source']
        for detail in effects if stage == '01_CAMP' else candidates:
            original = frozen[(detail['dataset'], detail['feature_name'], detail['metabolite_key'])]
            row = dict.fromkeys(header('statistical_result.tsv') + extras, NA)
            reason = detail.get('disposition', 'FROZEN_EFFECT_DESCRIBED')
            status = 'DONE'
            if stage == '02_MAPPING':
                status = 'NOT_EVALUABLE' if reason == 'UNIPROT_RECORD_NOT_EVALUABLE' else 'NEEDS_REVIEW'
            row.update(cancer='COAD', cohort=detail['dataset'], stage_id=stage, run_id=args.run_id,
                       analysis_version='coad_biochemical_v0.1_shared_view_v1',
                       analysis_type='frozen_effect_view' if stage == '01_CAMP' else 'biochemical_annotation_view',
                       metabolite_key=detail['metabolite_key'], metabolite_name=detail['feature_name'],
                       gene=detail.get('gene_symbol') or NA, unit='metabolite_feature' if stage == '01_CAMP' else 'feature_gene_pair',
                       status=status, reason=reason, source_id=source if stage == '01_CAMP' else catalog_path,
                       human_gene_id=detail.get('human_gene_id') or NA, original_effect=original['hedges_g'], original_q=original['effect_fdr'],
                       identity_review_status=detail.get('identity_review_status', original['identity_confidence']),
                       detail_source=source if stage == '01_CAMP' else catalog_path)
            if stage == '01_CAMP':
                row.update(effect_type='Hedges_g', effect=original['hedges_g'], q_value=original['effect_fdr'],
                           test_family='Frozen CAMP effect_fdr; original family details not established here')
            rows.append(row)
        rows.sort(key=lambda r: tuple(r[k] for k in ['cohort', 'metabolite_name', 'metabolite_key', 'gene', 'human_gene_id']))
        write_tsv(output / 'results.tsv', header('statistical_result.tsv') + extras, rows)
        inputs = [source] + ([catalog_path, identity_path] if stage == '02_MAPPING' else [])
        inputs += ['docs/STAGE_STANDARD_CN.md', 'templates/statistical_result.tsv', 'code/coad/export_shared_delivery.py']
        write_tsv(output / 'source_manifest.tsv', ['source_path', 'version', 'sha256'],
                  [dict(source_path=p, version='repository snapshot', sha256=sha(ROOT / p)) for p in inputs])
        spec = dict(cancer='COAD', owner='B', stage_id=stage, run_id=args.run_id,
                    analysis_version='coad_biochemical_v0.1_shared_view_v1',
                    created_utc=datetime.now(timezone.utc).isoformat(),
                    code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    generator_sha256=sha(Path(__file__)), software={'python': platform.python_version()},
                    status='DONE' if stage == '01_CAMP' else 'PARTIAL',
                    parameters={'new_calculation': False, 'original_q_threshold_for_mapping': 0.05},
                    statistical_unit='frozen aggregate feature' if stage == '01_CAMP' else 'feature-human gene annotation pair',
                    test_family={'new_tests': 'NOT_RUN', 'original_q': 'Preserved effect_fdr; original full family size not asserted'},
                    seed=NA, record_count=len(rows),
                    missing_values='NA: no new P/q; n, intervals and original P are unavailable in the input effect snapshot; NA is not zero.',
                    interpretation='A format view of existing evidence, not new or independent evidence. All mapping rows require source/complex-role review before finalization.',
                    patient_analysis='NOT_RUN', mapping_identity_detail=identity_path if stage == '02_MAPPING' else NA)
        dump(output / 'analysis_spec.json', spec)
        result = ('159 条冻结效应记录；原 q<0.05 为 73 条（正向 33、负向 40），86 条背景保留。'
                  if stage == '01_CAMP' else
                  '73 条显著特征身份初查；974 条候选组合，763 条有反应注释支持，其中 674 条当前无特定身份暂挂，涉及 59 个特征和 458 个基因。98 条身份暂挂、190 条反应特异性待查、12 条 reviewed 人类条目不可评估。身份明细见 ' + identity_path + '。')
        command = f'python code/coad/export_shared_delivery.py --run-id {args.run_id} --stage {stage} --validate-only'
        bodies = [
            '将现有 COAD 结果对齐 main 已发布的七阶段规范；不重算或另做候选筛选。',
            '原冻结效应表及本账号独立建立的数据库证据；具体来源和哈希见 source_manifest.tsv。结果表只是同一来源的格式视图。',
            result,
            '原 g 是标准化效应，不是倍数；原 q 不是新的关联 q。数据库反应注释不等于 COAD 患者支持或功能成立。',
            '原始身份、部分底物特异性与复合体角色未全部核实；患者身份和矩阵未读取。server165 别名无法解析。results.tsv 的 NA 及 status/reason 保留缺项，不表示阴性。',
            '效应描述批次完成；映射整体 PARTIAL，保留全部待审和不可评估组合。统一展示不改变历史数值或 A/BRCA 结果。',
            '获得 server165 正确连接信息后核对作者注释、样本身份和 RNA 覆盖；功能取证可并行。',
            f'校验已生成包：`{command}`。重新导出时去掉 --validate-only 并使用新的 UTC run-id；程序拒绝覆盖已有目录。发布状态查 coordination/publications/COAD.json。'
        ]
        report = '# COAD ' + stage + '\n\n' + '\n\n'.join('## ' + h + '\n\n' + b for h, b in zip(HEADINGS, bodies)) + '\n'
        (output / 'README_CN.md').write_text(report, encoding='utf-8', newline='\n')
        validate_stage(output, stage, frozen)
        dump(output / 'validation.json', dict(status='PASS', checked_utc=datetime.now(timezone.utc).isoformat(),
             checks=['shared field prefix and row order', 'unique keys', 'source byte hashes', 'all frozen effects/q unchanged', '159 or 974 rows', 'report section order', 'no new mapping P/q'],
             artifact_sha256={p.name: sha(p) for p in output.iterdir() if p.is_file()},
             not_validated=['original chemical identification', 'patient independence and pairing', 'patient-level RNA or correlation', 'functional causality', 'GitHub publication; see publication ledger']))
        print('EXPORTED and VALIDATED:', output.relative_to(ROOT))


if __name__ == '__main__':
    main()
