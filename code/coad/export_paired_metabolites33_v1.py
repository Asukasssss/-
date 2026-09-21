"""Render only public aggregate paired results; no patient-level input."""
import argparse, csv, hashlib, html, json
from collections import defaultdict
from pathlib import Path


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True)
    a=ap.parse_args();out=a.results
    names=json.loads(Path(__file__).with_name('paired_metabolite_names_cn.json').read_text(encoding='utf-8'))
    rows=read(out/'primary33.tsv');obs={r['metabolite_name']:r for r in read(out/'observed_pairs.tsv')}
    assert len(rows)==159 and set(r['metabolite_name'] for r in rows)==set(names)
    original_validation=json.loads((out/'validation.json').read_text())
    for name,h in original_validation['output_sha256'].items():
        assert hashlib.sha256((out/name).read_bytes()).hexdigest()==h
    assert set(obs)==set(names)
    sig=[r for r in rows if float(r['q_value'])<.05]
    lines=['# COAD：33名患者肿瘤与本人正常组织的代谢物配对结果','',
        '共159项特征，94项通过本轮配对Wilcoxon检验的BH校正（新q<0.05）：41项多数升高、52项多数降低，另1项没有超过半数同向变化。',
        '', '“32/33升高”表示33位患者中，32位的肿瘤数值高于自己的正常组织。比例分母保留相等配对，不仅计算非零差值。',
        '', '**† 表示至少1对涉及作者填补值；不是33对全部实测。** 无†的项目33对两侧在作者填补前均有有限值。中文名仅作显示，来源英文名和化学身份限制在完整明细中保留。','']
    for field,title in [('up_pairs','多数配对一致升高：41项'),('down_pairs','多数配对一致降低：52项')]:
        groups=defaultdict(list)
        for r in sig:
            if int(r[field])>16:
                label=names[r['metabolite_name']]+('†' if int(r['both_observed_pairs'])<33 else '')
                groups[int(r[field])].append(label)
        lines += ['## '+title,'','|同向对数|比例|代谢物——该行每一项均达到左侧对数|','|---:|---:|---|']
        for count,labels in sorted(groups.items(),reverse=True):
            lines.append('|{}/33|{:.1f}%|{}|'.format(count,100*count/33,'；'.join(labels)))
        lines.append('')
    lines += ['## 显著但不属于“多数同向”：1项','',
        '葡萄糖醛酸：7对升高、16对降低、10对相等；配对差值中位数为0，q=0.04465。不能写成多数患者降低，也不能把相等配对删除后改写分母。','',
        '## 填补前可用性与统计解释','',
        '77/159项具备33对填补前完整数值，82项至少有一对涉及作者填补。观察值子集敏感性单独锁定159项检验，157项至少8对可评估，97项q<0.05；不能与主分析相加为独立证据。',
        '', '主分析94项中，92项在可用性敏感性中仍通过BH；2-脱氧葡萄糖-6-磷酸（15对）和5-羟赖氨酸（14对）未通过。异柠檬酸敏感性仅9对，使用预定正态近似检验，需保留小样本近似限制。',
        '', '本轮统计是33对组织的代谢物差异，不能替代肿瘤内部代谢物—RNA关联，也不是外部验证。原CAMP q与本轮新配对q分列保留。Wilcoxon比较差值的符号与秩，q并不是单纯对“升高人数比例”做检验。',
        '', '## 完整159项逐项明细','',
        '|代谢物|来源英文名|升高|降低|相等|新配对q|填补前双侧完整对数|子集升高/降低/相等|子集q|',
        '|---|---|---:|---:|---:|---:|---:|---|---:|']
    def qfmt(q):return 'NA' if q=='NA' else '{:.5g}'.format(float(q))
    translated=[]
    for r in rows:
        s=obs[r['metabolite_name']]
        lines.append('|{}|{}|{}|{}|{}|{}|{}|{}/{}/{}|{}|'.format(names[r['metabolite_name']],r['metabolite_name'],r['up_pairs'],r['down_pairs'],r['equal_pairs'],qfmt(r['q_value']),r['both_observed_pairs'],s['up_pairs'],s['down_pairs'],s['equal_pairs'],qfmt(s['q_value'])))
        translated.append(dict(代谢物=names[r['metabolite_name']],来源英文名=r['metabolite_name'],配对数=r['n'],升高对数=r['up_pairs'],降低对数=r['down_pairs'],相等对数=r['equal_pairs'],升高百分比=r['up_percent'],降低百分比=r['down_percent'],本轮配对P=r['p_value'],本轮配对q=r['q_value'],原CAMPq=r['original_q'],中位配对差=r['effect'],差值尺度=r['effect_type'],填补前双侧完整对数=r['both_observed_pairs'],子集升高=s['up_pairs'],子集降低=s['down_pairs'],子集相等=s['equal_pairs'],子集P=s['p_value'],子集q=s['q_value'],子集状态=s['status']))
    md='\n'.join(lines)+'\n';(out/'PAIRED_COUNTS_CN.md').write_text(md,encoding='utf-8')
    with (out/'paired_counts_CN.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(translated[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(translated)
    body=[];table=False
    for line in lines:
        if line.startswith('|'):
            if line.startswith('|---'): continue
            tag='td' if table else 'th'
            if not table:body.append('<table>');table=True
            body.append('<tr>'+''.join('<'+tag+'>'+html.escape(c.strip())+'</'+tag+'>' for c in line.strip('|').split('|'))+'</tr>')
        else:
            if table:body.append('</table>');table=False
            if line.startswith('## '):body.append('<h2>'+html.escape(line[3:])+'</h2>')
            elif line.startswith('# '):body.append('<h1>'+html.escape(line[2:])+'</h1>')
            elif line:body.append('<p>'+html.escape(line.replace('**',''))+'</p>')
    if table:body.append('</table>')
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>COAD 33对代谢物结果</title><style>body{font:16px/1.8 "Microsoft YaHei",sans-serif;max-width:1380px;margin:40px auto;padding:0 24px;color:#20252b}h1{font-size:27px}h2{margin-top:36px;font-size:22px}table{border-collapse:collapse;width:100%;margin:20px 0}th,td{padding:12px 16px;text-align:left;border-bottom:1px solid #ddd}th{background:#f2f6f9}td:nth-child(1){min-width:75px}tr:nth-child(even){background:#fafbfc}p{max-width:1150px}@media print{body{font-size:10px;margin:0}th,td{padding:5px}tr{break-inside:avoid}}</style><body>'+''.join(body)+'</body></html>'
    (out/'PAIRED_COUNTS_CN.html').write_text(doc,encoding='utf-8')
    assert sum(float(obs[r['metabolite_name']]['q_value'])<.05 for r in sig)==92
    print(json.dumps(dict(rows=len(rows),significant=len(sig),still_significant_observed=92)))


if __name__=='__main__': main()
