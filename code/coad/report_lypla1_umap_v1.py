"""Validate public LYPLA1 UMAP outputs and register this COAD batch."""
import csv,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
run='20260928T031628Z_lypla1_umap_v1';rel='results/COAD/06_EXTERNAL/'+run;out=ROOT/rel
v=json.loads((out/'validation.json').read_text());assert v['status']=='PASS' and v['cells']==68702
with (out/'expression_coverage.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert sum(int(r['cells']) for r in rows)==v['cells']
assert sum(int(r['detected_cells']) for r in rows)==v['detected_cells']
for name in ['LYPLA1_UMAP','LYPLA1_UMAP_tissue']:
    for ext in ['png','svg']:assert (out/(name+'.'+ext)).stat().st_size>10000
text=f'''# COAD：LYPLA1单细胞表达UMAP

## 本轮问题

在已完成的Uhlitz/GSE166555细胞类型UMAP上展示LYPLA1表达，方便与注释图直接对照。

## 输入与范围

复用`20260928T025513Z_uhlitz_umap_v1`的全部68,702个细胞、12位供者的UMAP坐标，不重算降维。包含38,233个肿瘤组织来源细胞及30,469个正常组织来源细胞，不是仅恶性上皮子集。新读取同一原始计数文件中的LYPLA1行，以显式细胞ID连接，全部匹配。输入SHA256及坐标来源见source_manifest.tsv。

表达尺度为`log1p(LYPLA1 UMI / 该细胞全基因UMI总数 × 10,000)`。分母沿用已在前批逐细胞核验的作者nCount_RNA；不使用HVG缩放值、插补值或批次校正后的表达。

## 实际结果

![LYPLA1表达UMAP](LYPLA1_UMAP.png)

共{v['detected_cells']:,}个细胞检出LYPLA1计数，{v['zero_cells']:,}个为零计数。这里是被采样细胞的描述，不能当作供者等权表达结果或组织真实细胞比例。

![组织分面](LYPLA1_UMAP_tissue.png)

两组织面板使用同一坐标范围及完整数据共同色标，没有按组单独调色。更浅背景显示另一组织的位置；灰点为本组织零计数，非零表达从低到高绘制，较高值在上层。高值后绘可使其更明显，不应由图上亮点覆盖面积估计检出比例。细胞类型及组织的计数和描述性表达见expression_coverage.tsv，均为合并细胞统计，非患者等权推断。

## 新手解释

每个点是一个细胞；颜色越接近色标的亮黄端，LYPLA1归一化RNA表达越高。灰色表示此次未检出计数，不证明生物学上完全不表达。横纵轴表示表达相似性的UMAP坐标，不是组织真实空间。

## 限制/反证

坐标为本项目前批重算，不是原作者提供的坐标；未做批次校正，供者和技术差异可能影响位置。肿瘤组织来源不等于恶性细胞，不能仅由两面板亮度判断恶性上皮高于非恶性上皮。此图不新增差异检验、P/q、机制或脂代谢通路结论。

## 当前决定

交付全细胞表达图及组织分面图（PNG、SVG），补充既有注释图。逐细胞坐标与表达留server165，不导出本机或GitHub。

## 下一步

结合此前细胞类型注释图阅读来源；既有LYPLA1上皮比较及高低组富集统计保持原样。

## 复现命令

在server165新独占目录中运行`python3 plot_lypla1_umap_v1.py --out <运行目录> --commit <代码锁定提交>`。脚本路径为`code/coad/plot_lypla1_umap_v1.py`，使用前批已保存的坐标；软件、色标和实际输入哈希见随附JSON/TSV。本地运行`python code/coad/report_lypla1_umap_v1.py`核对公开计数并生成说明。
'''
(out/'README_CN.md').write_text(text,encoding='utf-8',newline='\n')
(out/'code_manifest.json').write_text(json.dumps([dict(path='code/coad/'+name,sha256=hashlib.sha256((ROOT/'code/coad'/name).read_bytes()).hexdigest()) for name in ['plot_lypla1_umap_v1.py','report_lypla1_umap_v1.py']],indent=2)+'\n',encoding='utf-8',newline='\n')
p=ROOT/'coordination/stages/COAD.tsv'
with p.open(encoding='utf-8',newline='') as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;ledger=list(r)
if not any(r['run_id']==run for r in ledger):
    for stage,scope in {'06_EXTERNAL':'LYPLA1 expression overlay;68702 cells,12 donors;frozen Uhlitz project UMAP','07_INTEGRATION':'Whole-cohort and tissue-split expression figures;same color scale;no new tests'}.items():
        ledger.append(dict(cancer='COAD',stage_id=stage,run_id=run,analysis_version='COAD_LYPLA1_UMAP_v1',status='DONE',scope=scope,result_path=rel,code_path='code/coad/plot_lypla1_umap_v1.py',git_branch='analysis/coad-initial',reason='Exact cell-ID linkage;log1pCP10K;no batch correction;pooled descriptive plots',next_action='Read alongside author annotations;not malignancy or mechanism inference'))
    with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sorted(ledger,key=lambda r:(r['stage_id'],r['run_id'])))
p=ROOT/'reports/COAD/CURRENT_STATUS_CN.md';old=p.read_text(encoding='utf-8')
if not old.startswith('# 最新：LYPLA1表达UMAP已完成'):
    old=old.replace('# 最新：Uhlitz单细胞UMAP注释图已完成','# Uhlitz单细胞UMAP注释图已完成',1)
    p.write_text('# 最新：LYPLA1表达UMAP已完成\n\n[全细胞及组织分面图](../../'+rel+'/README_CN.md)：复用Uhlitz全68,702细胞、12供者的既有UMAP坐标，叠加LYPLA1 log1pCP10K表达。细胞ID全部匹配，肿瘤/正常共享色标，不重算坐标、不做新P/q；不将组织来源当恶性身份。\n\n'+old,encoding='utf-8',newline='\n')
for p in out.iterdir():
    if p.suffix in ['.tsv','.json','.md','.svg']:p.write_text('\n'.join(x.rstrip() for x in p.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
print(json.dumps(v))
