import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/BRCA/04_ROBUSTNESS/20260921T112611Z_metabolite318_paired_v1/comparison318.tsv'

def main():
    run = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_paired_p_rank_v1'
    out = ROOT / 'results/BRCA/04_ROBUSTNESS' / run
    out.mkdir(parents=True, exist_ok=False)
    with SOURCE.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert len(rows) == 318 and len({r['metabolite_key'] for r in rows}) == 318
    for r in rows:
        h, l, e = (int(r[x]) for x in ('pairs_higher', 'pairs_lower', 'pairs_equal'))
        assert h + l + e == 45 and float(r['n']) == 45
        r.update(majority_direction='UP' if h > l else 'DOWN' if l > h else 'TIED',
                 same_direction_pairs=max(h, l), same_direction_fraction=max(h, l)/45,
                 raw_p_lt_005=float(r['p_value']) < .05,
                 majority_matches_mean=(h-l)*float(r['effect']) > 0)
    rows.sort(key=lambda r: (-r['same_direction_pairs'], r['metabolite_key']))
    selected = [r.copy() for r in rows if r['raw_p_lt_005']]
    for i, r in enumerate(selected, 1):
        r['display_order'] = i
        r['proportion_rank'] = 1 + sum(x['same_direction_pairs'] > r['same_direction_pairs'] for x in selected)
    def write(path, data, sep='\t', encoding='utf-8'):
        with path.open('w', encoding=encoding, newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter=sep)
            w.writeheader(); w.writerows(data)
    write(out/'all318_direction_counts.tsv', rows)
    write(out/'p_lt_005_sorted.tsv', selected)
    write(out/'p_lt_005_sorted_excel.csv', selected, ',', 'utf-8-sig')
    counts = {d:sum(r['majority_direction']==d for r in selected) for d in ('UP','DOWN','TIED')}
    summary = dict(n_total=318,n_p_lt005=len(selected),n_q_lt005=sum(float(r['q_value'])<.05 for r in selected),directions=counts,
                   majority_mean_disagreements=sum(not r['majority_matches_mean'] for r in selected))
    table = ['|比例名次|代谢物原名|多数方向|升高/降低/持平|同向比例|原始P|原q|', '|---:|---|---|---|---:|---:|---:|']
    for r in selected:
        d={'UP':'升高','DOWN':'降低','TIED':'方向持平'}[r['majority_direction']]
        table.append(f"|{r['proportion_rank']}|{r['metabolite_name']}|{d}|{r['pairs_higher']}/{r['pairs_lower']}/{r['pairs_equal']}|{r['same_direction_fraction']:.1%}|{float(r['p_value']):.4g}|{float(r['q_value']):.4g}|")
    readme=f'''# BRCA 45对代谢物：原始P筛选及配对方向比例排序

## 本轮问题
按用户要求，45对主分析原始P<0.05，按多数变化方向的配对比例从高到低列出。

## 输入与范围
复用318项已有配对主分析；新统计检验为0，原效应/P/q/检验范围均不改。

## 实际结果
{len(selected)}项原始P<0.05；其中{summary['n_q_lt005']}项原q<0.05。多数升高{counts['UP']}项，多数降低{counts['DOWN']}项，方向持平{counts['TIED']}项。

## 新手解释
同向比例=max(升高对数,降低对数)/45；分母包括持平配对。方向按对数多数判定，不按均值判定。比例相同为并列名次，展示次序按稳定代谢物键排列，不再按P或q排名。比例不是效应幅度，也不是逐患者显著率。

## 限制/反证
使用作者处理值，包含原有缺失处理。原始P筛选属于名义统计支持，不冒充FDR控制。多数方向和均值方向不一致项目数为{summary['majority_mean_disagreements']}，保留原效应及标记。完整表保留可用性敏感性字段与全部318项，X特征不猜身份。

## 当前决定
本表仅为展示排序，不是新的靶点或稳健性排名。

## 下一步
解释具体条目时同时查缺失敏感性、效应及化学身份。

## 复现命令
`python code/brca_paired_p_rank_v1.py`

## 完整排序
''' + '\n'.join(table) + '\n'
    (out/'README_CN.md').write_text(readme,encoding='utf-8')
    (out/'analysis_spec.json').write_text(json.dumps(dict(version='paired_p_rank_v1',filter='existing paired_processed p_value < 0.05',sort='max(higher,lower)/45 descending; key ascending for ties',new_tests=0),indent=2),encoding='utf-8')
    (out/'validation.json').write_text(json.dumps(dict(summary=summary,unique_keys=318,all_pair_counts_sum_45=True,original_fields_preserved=True),indent=2),encoding='utf-8')
    (out/'source_manifest.tsv').write_text('path\tsha256\n'+str(SOURCE.relative_to(ROOT)).replace('\\','/')+'\t'+hashlib.sha256(SOURCE.read_bytes()).hexdigest()+'\n',encoding='utf-8')
    index=ROOT/'coordination/stages/BRCA.tsv'
    with index.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; stages=list(reader)
    vals=['BRCA','04_ROBUSTNESS',run,'paired_p_rank_v1','DONE','318 paired metabolites; nominal P filter; direction count ranking',str(out.relative_to(ROOT)).replace('\\','/'),'code/brca_paired_p_rank_v1.py','analysis/brca-functional-review-20260919','Aggregate display only; no new tests; original q retained','Read missingness alongside directional agreement']
    stages.append(dict(zip(fields,vals))); stages.sort(key=lambda r:(r['stage_id'],r['run_id']))
    with index.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(stages)
    print(json.dumps(dict(output=str(out),**summary)))

if __name__ == '__main__':
    main()
