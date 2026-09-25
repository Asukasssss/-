"""Export PDAC frozen aggregate effects without recalculation or patient data."""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPECTED = '1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597'

def write_table(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-id', required=True)
    args = parser.parse_args()
    if Path(args.run_id).name != args.run_id or '/' in args.run_id or '\\' in args.run_id:
        raise ValueError('run-id must be a directory name')
    source = ROOT / 'reference/camp/cancer_effects.tsv'
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    assert digest == EXPECTED, 'Frozen source hash mismatch'
    with source.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        rows = [r for r in reader if r['cancer'] == 'PDAC']
    sig = [r for r in rows if float(r['effect_fdr']) < 0.05]
    assert len(rows) == 303 and len(sig) == 42
    keys = [(r['dataset'], r['feature_name'], r['metabolite_key']) for r in rows]
    assert len(keys) == len(set(keys))
    out = ROOT / 'results/PDAC/01_CAMP' / args.run_id
    out.mkdir(parents=True, exist_ok=False)
    write_table(out / 'results.tsv', rows, fields)
    write_table(out / 'significant_features.tsv', sig, fields)
    for name, expected in [('results.tsv', rows), ('significant_features.tsv', sig)]:
        with (out / name).open(encoding='utf-8', newline='') as f:
            assert list(csv.DictReader(f, delimiter='\t')) == expected
    summary = {'all_effects': len(rows), 'original_q_lt_0_05': len(sig),
               'positive_g_significant': sum(float(r['hedges_g']) > 0 for r in sig),
               'negative_g_significant': sum(float(r['hedges_g']) < 0 for r in sig),
               'significant_unknown_x': sum(r['feature_name'].startswith('X-') for r in sig),
               'cohorts': dict(Counter(r['dataset'] for r in rows)),
               'significant_pathways': dict(Counter(r['super_pathway'] for r in sig))}
    spec = {'analysis_version': 'pdac_frozen_export_v1', 'run_id': args.run_id,
            'operation': 'Exact row selection; no statistical recomputation',
            'filter': 'cancer == PDAC; significant subset: original effect_fdr < 0.05',
            'unit': 'cohort x actual metabolomic feature aggregate effect',
            'n': None, 'n_reason': 'Not provided in frozen aggregate snapshot',
            'test_family': 'Original CAMP family retained, no new tests',
            'original_p': None, 'p_reason': 'Raw P not present in snapshot; not inferred from q',
            'seed': None, 'python': platform.python_version(),
            'code_base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'source_sha256': digest}
    for name, obj in [('summary.json', summary), ('analysis_spec.json', spec),
                      ('validation.json', {'status': 'PASS', 'frozen_hash': True,
                        'all_fields_roundtrip_equal': True, 'unique_feature_keys': True,
                        'expected_counts': True, 'not_verified': ['server165 source matrices',
                        'patient independence', 'PDAC historical analyses', 'effect contrast orientation']})]:
        (out / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    write_table(out / 'source_manifest.tsv', [{'path': 'reference/camp/cancer_effects.tsv',
                'version': 'CAMP_Phase1_v1.0 frozen repository snapshot', 'sha256': digest}],
                ['path', 'version', 'sha256'])
    report = f'''# PDAC 冻结发现结果首批交付

## 本轮问题
接手 PDAC，确认冻结发现结果范围，为独立的癌种内分析建立入口。

## 输入与范围
使用仓库冻结汇总表，保留全部原始字段及字符串数值。仅选 PDAC，不重算效应或 q。
本表沿用历史列名：dataset 对应 cohort；hedges_g 对应标准化效应；effect_fdr 对应原始 q；feature_name 是实际特征名。没有原始 P、样本量或区间字段，不反推、不填零。

## 实际结果
全量 {len(rows)} 条，原 q<0.05 共 {len(sig)} 条；显著特征中 g 为正 {summary['positive_g_significant']} 条、负 {summary['negative_g_significant']} 条，未知 X 特征 {summary['significant_unknown_x']} 条。
results.tsv 保留全部背景；significant_features.tsv 为原显著子集。字段逐项回读相等，来源哈希与冻结值一致。

## 新手解释
42 条是代谢特征，不是 42 个基因或治疗靶点。Hedges g 是标准化效应，不是倍数。正负号沿用冻结值；在服务器确认原始对比编码前，不进一步解释为肿瘤升降。

## 限制/反证
本轮没有新增患者关联、基因映射、功能实验或独立验证。HIGH_ID 仅沿用旧身份标签，未知 X 不猜基因。
本机 SSH 的 server165 主机别名无法解析，作者映射、RNA 覆盖、患者独立性和既有 PDAC 结果尚未核对。

## 当前决定
01_CAMP 的冻结快照整理 DONE；PDAC 整体仍处起步阶段。用户已在本任务指定本账号负责 PDAC。

## 下一步
获取有效连接方式后，在服务器读取既有 PDAC 汇总及作者映射；复用适用统计，并核实直接人类酶/转运体关系及来源。功能和外部证据与患者关联并行推进；各批单独保存原 q 与新检验族。

## 复现命令
在仓库根运行 `python code/pdac/export_frozen.py --run-id NEW_UNIQUE_RUN_ID`。已有目录拒绝覆盖；计算只使用允许公开的汇总快照。
'''
    (out / 'README_CN.md').write_text(report, encoding='utf-8')
    print(json.dumps({'result_path': str(out), **summary}, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
