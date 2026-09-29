"""Public delivery of the user-requested full-cell reference layout."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='20260929T065300Z_brca_allcells_slide_v2'
OUT=ROOT/'results/BRCA/07_INTEGRATION'/RUN
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_BRCA_allcells_v2');VIEW.mkdir(parents=True,exist_ok=True)
name='11_BRCA_LYPLA1_全细胞与患者'
svg=OUT/(name+'.svg');svg.write_bytes(('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf8').splitlines())+'\n').encode('utf8'))
(OUT/'.gitattributes').write_bytes(b'* -text\n')
old=ROOT/'results/BRCA/06_EXTERNAL/20260923T040012Z_epithelial_paired20_v1/results.tsv'
d=pd.read_csv(OUT/'results.tsv',sep='\t');o=pd.read_csv(old,sep='\t');o=o[(o.cohort=='Wu2021')&(o.gene=='LYPLA1')]
pd.testing.assert_frame_equal(d.reset_index(drop=True),o.reset_index(drop=True),check_dtype=False)
coverage=pd.read_csv(OUT/'coverage.tsv',sep='\t');assert coverage.cells.sum()==100064 and len(coverage)==9
(OUT/'README_CN.md').write_text('''# 第11页 BRCA：全细胞定位与同患者证据

## 本轮问题
用户要求按所给参考图展示全部细胞的基因表达与细胞类型对照，再加同患者证据；替代上一版仅上皮子集UMAP。

## 输入与范围
Wu2021 / GSE176078，全部100064细胞、9类细胞；同一源H5AD的原有obsm/X_umap与celltype_major。LYPLA1用原始计数计算既有log1p(CP10K)表达尺度，未改变归一化。原始矩阵、逐细胞缓存及患者测量仅存server165。

## 实际结果
左：全细胞LYPLA1表达，零值灰色，阳性表达用viridis；阳性值99分位作为颜色上限，超出显示顶端颜色但不删细胞。中：相同细胞、坐标，多色原有细胞注释。红圈是Cancer Epithelial中心云团的视觉定位，非划分或检验边界，不用于筛选细胞。右：既有8位患者配对结果，7/8恶性上皮较高，平均差0.2153329622，双侧配对精确符号秩P=0.0390625。

## 新手解释
先看中图中癌上皮的位置，再对照左图实测表达；全图仍保留髓系、免疫和基质等表达背景。不能由圈选推断恶性上皮特异表达、在所有细胞中最高或正式富集。支持恶性上皮高于非恶性上皮的定量依据来自右侧患者配对。每条线是一位患者，端点是患者内该类细胞的平均log1p(CP10K)，横杠为患者等权组均值。UMAP全细胞范围与8人可配对范围不同。

## 限制与反证
参照是肿瘤内作者注释的非恶性上皮，不是健康乳腺。每患者每类至少20细胞，身份沿用作者标签；未新做CNV判定。作者Naive子集6/7同向，P=0.078125。按用户要求展示名义P，原q完整保留。单队列探索结果不代表跨队列复现。细胞比例不代表组织真实比例；图中色标不是用户示例中的ASNS log2(CPM+0.5)，不能照抄示例范围。

## 当前决定
替代旧第11页，保留SYSU页眉与配对证据。图由真实数据代码绘制，无IMAGE生成，无新聚类、降维或检验。

## 下一步
先确认全细胞展示方式，再据用户选择继续其他癌种或空间证据。

## 复现
server165新运行目录内放入字体、真实校标、冻结results.tsv，运行python3 ppt_brca_allcells_slide_v2.py RUN_DIR。仅复制public目录回本地，然后python code/ppt_brca_allcells_package_v2.py。原始数据及私有缓存不得上传GitHub。
''',encoding='utf8')
v=json.loads((OUT/'validation.json').read_text());v.update(visual_review='PASS: full-slide image inspected; labels and plots readable',local_frozen_statistics_match=True,public_only=True)
(OUT/'validation.json').write_text(json.dumps(v,indent=2),encoding='utf8')
for ext in ['png','pdf','svg']:shutil.copy2(OUT/(name+'.'+ext),VIEW/(name+'.'+ext))
for n in ['README_CN.md','results.tsv','coverage.tsv']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>BRCA · 全细胞与同患者证据</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1560px;margin:20px auto}}img{{width:100%;display:block;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>BRCA · 全细胞定位与同患者证据</h2><img src="{name}.png"><p><a href="{name}.pdf">高清PDF</a> · <a href="{name}.svg">SVG</a> · <a href="results.tsv">配对统计</a> · <a href="README_CN.md">绘图说明</a></p></main>',encoding='utf8')
pd.DataFrame([dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [ROOT/'code/ppt_brca_allcells_slide_v2.py',Path(__file__),old]]).to_csv(OUT/'local_source_manifest.tsv',sep='\t',index=False)
pd.DataFrame([dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='checksums.tsv']).to_csv(OUT/'checksums.tsv',sep='\t',index=False)
stage=ROOT/'coordination/stages/BRCA.tsv'
if RUN not in stage.read_text(encoding='utf8'):
 with stage.open('a',encoding='utf8') as f:f.write('\t'.join(['BRCA','07_INTEGRATION',RUN,'brca_allcells_slide_v2','DONE','All-cell reference-style UMAP and frozen paired donor evidence',str(OUT.relative_to(ROOT)).replace('\\','/'),'code/ppt_brca_allcells_slide_v2.py','presentation/camp-fourcancer-lypla1-20260928','100064 cells, original coordinates; locator not a gate; nominal P; no new tests','Review revised BRCA all-cell slide'])+'\n')
print(VIEW)
