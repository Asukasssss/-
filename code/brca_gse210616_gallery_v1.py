"""Render a local gallery from exported figures and aggregate results only."""
from pathlib import Path
import sys,json,csv,html
d=Path(sys.argv[1]);v=json.loads((d/'validation.json').read_text())
rows=list(csv.DictReader((d/'results.tsv').open(encoding='utf-8'),delimiter='\t'))
tr=''.join(f"<tr><td>{html.escape(r['measure'])}</td><td>{r['positive']}/{r['n']}</td><td>{float(r['mean']):.3f}</td><td>{float(r['p_value']):.4g}</td></tr>" for r in rows)
cards=''.join(f'<article id="{p.stem}"><h2>{p.stem}</h2><img loading="lazy" src="{p.name}" alt="{p.stem} H&E, marker-inferred regions, LYPLA1"><a href="{p.name}" target="_blank">打开大图</a></article>' for p in sorted(d.glob('GSM*.png')))
page=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>LYPLA1 · GSE210616 空间分析</title>
<style>body{{max-width:1400px;margin:36px auto;padding:0 24px;font:16px/1.7 system-ui;color:#203047;background:#f6f8fb}}h1{{font-size:30px}}article,.box{{background:white;border:1px solid #dde4ec;border-radius:14px;padding:22px;margin:24px 0}}img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;text-align:left;border-bottom:1px solid #ddd}}input{{padding:12px;width:300px;max-width:90%;font-size:16px}}.note{{border-left:5px solid #c68929;padding:12px;background:#fff7e8}}</style>
<h1>LYPLA1 在独立 TNBC 队列中的空间表达</h1><p>GSE210616 · {v['patients']} 位患者 · {v['sections_analyzed']} 张切片 · {v['spots_after_QC']:,} 个质控后空间点</p>
<p class="note">本轮是标记基因辅助的探索性分析。中图区域由表达计算推断，不是病理医生标注；不能据此认定恶性细胞或正常上皮。LYPLA1 未参与区域评分。治疗分层尚不可评估。</p>
<div class="box"><h2>先看患者间是否一致</h2><img src="patient_summary.png"><p>每条代表一位患者；同一患者切片等权汇总。左：控制总 UMI 的秩相关；右：上皮标记富集点减去间质标记富集点的平均标准化表达。正值不是倍数变化。</p><table><thead><tr><th>指标</th><th>正方向患者</th><th>均值</th><th>双侧符号检验 P</th></tr></thead><tbody>{tr}</tbody></table><p>P 为名义 P；敏感性结果相关，不是多个独立验证。rho 是相关系数；delta_10/20/30 表示每区最低点数阈值。delta_pseudobulk 是区域合并计数后标准化的敏感性结果。</p></div>
<div class="box"><h2>逐张看空间位置</h2><p>左：H&E 原图。中：粉色上皮标记富集、蓝色间质、绿色免疫、橙色内皮、灰色不确定。右：LYPLA1，黑色低、紫红中、黄色高。色标按各切片设置，不可凭颜色直接比较患者表达高低。组织外规则圆点是芯片定位标记。</p><input id="q" placeholder="按 GSM 编号筛选" oninput="document.querySelectorAll('article').forEach(x=>x.hidden=!x.id.toLowerCase().includes(this.value.toLowerCase()))"></div>{cards}
<p>来源：<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616">GSE210616</a>。全部切片均展示；没有按 LYPLA1 的结果挑选样本。</p></html>'''
(d/'index.html').write_text(page,encoding='utf-8')
assert len(list(d.glob('GSM*.png')))==43
print(d/'index.html')
