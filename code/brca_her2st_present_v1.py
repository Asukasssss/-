"""Public aggregate report, gallery and stage registration; no patient numeric exports."""
from pathlib import Path
import json,csv,sys,shutil
base=Path('results/BRCA/06_EXTERNAL/20260928T024801Z_her2st_v1');gallery=Path(sys.argv[1])
v=json.loads((base/'validation.json').read_text());v['git_commit']=v.get('git_commit',v['git_tree']);v['git_tree']='a7d6df2e392afde7357f2f9697b2748e95a3fb33';v['analysis_code_commit']='7a3f453'
(base/'validation.json').write_text(json.dumps(v,indent=2))
report='''# Andersson HER2ST：LYPLA1 病理区域分析

## 本轮问题
使用作者病理医生标注，比较浸润癌区与结缔组织、乳腺腺体、原位癌区的LYPLA1表达。

## 输入与范围
8位HER2+患者，各一张有病理标注的切片，共3481个作者处理后组织点。585点检出LYPLA1（16.8%）。患者A–H的映射来自作者命名和论文；G使用G2，其他使用各自第1张。仅分析这8张，不将未标注的其余切片自动归类。作者仓库有额外J1标签，但没有对应本研究A–H矩阵，未纳入。
下载34个文件，使用Git blob SHA1与SHA256检查完整性，来源固定至作者提交85df7411b987a7dec3feff9718b869a7a82091ec。源数据经RAM转发仅存server165。

## 实际结果
每区至少20点为主分析；10点为预设敏感性。所有P为患者方向的双侧精确符号检验；少于3位患者不报告推断P。

|比较|可用患者|较高患者|平均log1p表达差|P|
|---|---:|---:|---:|---:|
|浸润癌 vs 结缔组织，主分析|7|7|+0.1769|0.015625|
|同上，区域合并计数后标准化|7|5|+0.1464|0.453125|
|同上，控制总RNA计数|7|5|回归系数+0.0627|0.453125|
|浸润癌 vs 乳腺腺体，主分析|2|2|+0.3164|NA，患者不足|
|浸润癌 vs 原位癌，主分析|3|1|+0.0445|1.0|

每区最低10点时，浸润癌—结缔组织可纳入8人，8/8同向，P=0.0078125；区域合并计数和深度调整均为6/8，P=0.28906。腺体比较可纳入3人，三种方法均3/3同向，但P=0.25。原位癌比较仍为3人、1/3同向，没有一致升高证据。主分析E患者结缔组织仅17点而不纳入，不是因为LYPLA1方向排除。

## 新手解释
这次红色区域确实是病理医生标出的浸润癌，区域身份比标记推断明确。基础结果支持癌区相对结缔组织较高，但对汇总方法和总RNA计数敏感；不能称稳健的癌细胞内在上调。乳腺腺体比较虽同向，患者少，暂不能确认癌上皮高于正常上皮。

## 限制/反证
- 早期ST为多细胞空间点；病理癌区仍混有免疫/间质细胞。乳腺腺体标签不自动等同纯正常上皮。
- LYPLA1检出点比例仅16.8%；未填补或平滑零值。零计数不证明细胞完全不表达。
- 作者F/G标签有9行缺失坐标，未推断坐标；所有3481个表达矩阵点仍能唯一匹配有效标签，无未匹配矩阵点。额外标签点不加入表达分析。
- 标签与坐标通过唯一阵列键连接，两套像素坐标具有缩放及Y轴翻转；变换残差均小于0.02像素，绘图使用原图对应坐标。已核对示例图组织对齐。
- 调整后的不显著结果、原位癌比较不一致均保留。名义P，没有按结果改阈值、筛患者或挑图。
- 不将8张切片扩称36张病理验证，不将区域差异称酶活、功能机制或泛癌结论。

## 当前决定
可以作为有明确病理依据的有限空间支持；不能取代恶性—正常上皮的患者级单细胞证据。原位癌—浸润癌进展叙事目前没有支持。

## 下一步
若继续外部验证，优先扩大明确非恶性腺体对照和患者数；完整保留本轮反证，不继续更换阈值寻找显著。

## 复现命令
将code/brca_her2st_analyze_v1.py复制到新服务器运行目录，以该目录为唯一参数运行；环境OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1，Python /home/xuzx/miniconda3/envs/ov/bin/python。
源目录：/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_Andersson_HER2ST。
运行目录：/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260928T024801Z_her2st_v1。
患者和空间点数值只存服务器；本地及GitHub仅汇总、代码、图和参数。

来源：[论文](https://www.nature.com/articles/s41467-021-26271-2)、[作者仓库](https://github.com/almaan/her2st)、[Zenodo](https://zenodo.org/records/4751624)。
'''
(base/'README_CN.md').write_text(report,encoding='utf-8')
cards=''.join(f'<section><h2>患者 {s[0]} · {s}</h2><img src="{s}.png" loading="lazy"><a href="{s}.png">打开大图</a></section>' for s in ['A1','B1','C1','D1','E1','F1','G2','H1'])
page=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>LYPLA1 HER2ST 病理区域分析</title><style>body{{max-width:1400px;margin:30px auto;padding:0 24px;background:#f6f8fb;color:#24364b;font:16px/1.8 system-ui}}section{{background:white;padding:22px;border:1px solid #dbe1ea;border-radius:12px;margin:24px 0}}img{{width:100%}}.note{{background:#fff3d9;padding:18px}}table{{border-collapse:collapse;width:100%}}td,th{{text-align:left;padding:10px;border-bottom:1px solid #ddd}}</style><h1>LYPLA1 · HER2+乳腺癌病理区域分析</h1><p>8位患者 · 8张明确病理标注切片 · 3481个空间点 · LYPLA1检出16.8%</p><p class="note">本轮区域来自作者病理医生标注。癌区相对结缔组织的基础差异一致，但区域合并计数、深度调整后证据减弱；乳腺腺体对照患者不足。</p><section><h2>先看完整结论</h2><table><tr><th>比较</th><th>较高患者/可用患者</th><th>患者方向P</th></tr><tr><td>浸润癌 vs 结缔组织（≥20点/区）</td><td>7/7</td><td>0.0156</td></tr><tr><td>同上，区域合并计数</td><td>5/7</td><td>0.4531</td></tr><tr><td>同上，控制总RNA计数</td><td>5/7</td><td>0.4531</td></tr><tr><td>浸润癌 vs 乳腺腺体</td><td>2/2</td><td>患者不足，不作推断</td></tr><tr><td>浸润癌 vs 原位癌</td><td>1/3</td><td>1.0</td></tr></table><p>降低至≥10点/区时，癌—结缔组织8/8同向，P=0.0078；深度调整6/8，P=0.2891。癌—腺体3/3同向，P=0.25。不能把不同阈值结果当独立验证。</p><img src="summary.png"></section><section><h2>读图方法</h2><p>左：H&E；中：作者病理区域；右：LYPLA1标准化表达。红色为浸润癌，橙色原位癌，绿色乳腺腺体，蓝色结缔组织。所有右图使用相同0–2色标，超过2以最高色显示。黑点为零计数，不等于绝对不表达。病理癌区不是纯癌细胞。</p></section>{cards}</html>'''
(gallery/'index.html').write_text(page,encoding='utf-8');assert len(list(gallery.glob('*.png')))==9
(base/'figures').mkdir(exist_ok=True);shutil.copy2(gallery/'summary.png',base/'figures'/'summary.png')
index=Path('coordination/stages/BRCA.tsv')
with index.open(encoding='utf-8-sig') as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;rows=list(r)
rows=[r for r in rows if r['run_id']!=base.name]
rows.append(dict(zip(fields,['BRCA','06_EXTERNAL',base.name,'her2st_v1','DONE','8 pathology annotated sections; 3 contrasts with depth sensitivities',base.as_posix(),'code/brca_her2st_analyze_v1.py','analysis/brca-lypla1-spatial-20260927','Limited gland comparison; depth-adjusted cancer-stroma non-significant','Retain full evidence; expand verified gland controls'])))
with index.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(sorted(rows,key=lambda x:(x['stage_id'],x['run_id'])))
print(gallery/'index.html')
