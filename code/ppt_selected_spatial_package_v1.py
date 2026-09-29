from pathlib import Path
import json,shutil,hashlib
from PIL import Image
from pypdf import PdfWriter
ROOT=Path(__file__).resolve().parents[1]
RUN='20260929T_selected_spatial_v1';O=ROOT/'results/BRCA/07_INTEGRATION'/RUN
V=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_selected_spatial_v1');V.mkdir(exist_ok=True)
pages=['15_BRCA_LYPLA1_空间切片统一版','16_PRAD_LYPLA1_空间切片统一版']
w=PdfWriter()
for name in pages:w.append(O/(name+'.pdf'))
with (O/'LYPLA1_选定空间切片_两页.pdf').open('wb') as f:w.write(f)
for p in O.iterdir():
 if p.suffix in ['.png','.pdf','.svg']:shutil.copy2(p,V/p.name)
cards=''.join(f'<section><h2>{name}</h2><a href="{name}.png" target="_blank"><img src="{name}.png"></a><p><a href="{name}.pdf">PDF</a> · <a href="{name}.svg">SVG</a> · <a href="{name}.png">PNG</a></p></section>' for name in pages)
single=''.join(f'<section><h2>{p.stem}</h2><a href="{p.name}" target="_blank"><img loading="lazy" src="{p.name}"></a></section>' for p in sorted(V.glob('*.png')) if p.stem not in pages)
(V/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>LYPLA1 · 选定四张空间切片</title><style>body{margin:0;background:#e9eeeb;color:#164031;font-family:"Microsoft YaHei",sans-serif}nav{position:sticky;top:0;padding:16px 4%;background:white;z-index:5}main{max-width:1680px;margin:20px auto}img{width:100%;display:block}section{margin:25px 0;background:white;padding:15px;border-radius:10px}a{color:#00553b}h2{font-size:18px}p{padding:8px}</style><nav>选定四张 · 统一重绘　<a href="LYPLA1_选定空间切片_两页.pdf">下载两页 PDF</a>　<a href="#single">四张独立大图</a></nav><main>'+cards+'<h2 id="single">四张独立大图</h2>'+single+'</main>',encoding='utf-8')
for p in O.glob('*.png'):
 with Image.open(p) as im:im.verify()
checks=json.loads((O/'display_checks.json').read_text(encoding='utf-8'))
assert [x['sample'] for x in checks]==['V19T26-012_B1','V19T26-032_B1','H2_5','H2_1']
(O/'validation.json').write_text(json.dumps(dict(status='PASS',selected_sections_verified=True,image_decode='PASS',new_tests=False,visual_review='PASS: both complete slide PNGs inspected'),indent=2))
(O/'README_CN.md').write_text('''# 用户选定四张空间切片统一重绘
## 本轮问题
将用户选定的 BRCA V19T26-012_B1、V19T26-032_B1 与 PRAD H2_5、H2_1 统一为 H&E、病理区域、LYPLA1 实测表达三联图，按癌种制作两页组会图。
## 输入与范围
使用 server165 既有 QC 后坐标、标准化表达和作者区域注释。PRAD 来源目录只读。逐点数据不离开服务器。
## 实际结果
四张独立三联图与两页 16:9 组会页面，提供 PNG/PDF/SVG。两页另合并 PDF。仅重绘，不改变检验、不重跑配准，不插补或平滑表达。
## 新手解释
三列依次为完整 H&E、病理区域、LYPLA1 实测点位。每一行三图使用相同坐标范围和组织方向。癌区统一红色；良性腺体蓝色；免疫区黄色。连续表达使用 viridis，浅灰单独表示零表达。四图共用 0–3 的 log1p(CP10K) 色标，上限截断数记录在 display_checks.json。
## 限制与反证
BRCA 代表不同两位患者；未标注区不能当正常上皮。PRAD 两张均来自一位患者。不同平台表达不可直接据色彩强弱比较；统一色标仅为展示标准。用户选定代表图不替代完整队列证据。图中无新增统计显著性结论。
## 当前决定
完成指定切片的组会版式与独立大图，不生成 PPT。
## 下一步
用户审阅版式后再决定纳入其他切片。
## 复现
server165 新建 BRCA/A 展示运行目录，python ppt_selected_spatial_v1.py RUN_DIR；仅同步 public 图形及汇总到本地；python code/ppt_selected_spatial_package_v1.py。
''',encoding='utf-8')
(O/'checksums.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='checksums.json'},indent=2))
print(V)
