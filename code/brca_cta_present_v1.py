"""Build public aggregate report and Chinese gallery from completed server results."""
from pathlib import Path
import json,csv,html
ROOT=Path(__file__).resolve().parents[1];RUN='20260928T060000Z_cta_lypla1_v1'
R=ROOT/'results/BRCA/06_EXTERNAL'/RUN
G=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/LYPLA1_CTA_spatial')
s=json.loads((R/'summary.json').read_text());v=json.loads((R/'validation.json').read_text());a=s['aggregate_effects']
def fmt(x):return 'NA' if x is None else f'{x:.4g}'
refnames={'Immune':'免疫区','Unannotated':'未标注混合区域'}
metrics={'mean_lognorm_delta':'平均标准化表达差','pseudobulk_log2ratio':'区域汇总表达 log2 比值','depth_mt_adjusted_delta':'校正总计数和线粒体比例后的表达差'}
primary=[x for x in a if x['purity']==.8]
table='|比较|指标|患者数|方向为正|患者效应均值|P值|\n|---|---|---:|---:|---:|---:|\n'+''.join(f"|癌区 vs {refnames[x['reference']]}|{metrics[x['metric']]}|{x['patients']}|{x['positive_patients']}/{x['patients']}|{fmt(x['mean_patient_effect'])}|{fmt(x['p_value'])}|\n" for x in primary)
cols=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t');rows=[]
base=dict.fromkeys(cols,'NA');base.update(cancer='BRCA',cohort='BreastCancer_CTA_2025',stage_id='06_EXTERNAL',run_id=RUN,analysis_version='cta_v1',gene='LYPLA1',unit='patient',source_id='Zenodo15211538',test_family='LYPLA1 two region contrasts x three metrics x two area thresholds; nominal P',family_n_evaluable=sum(x['p_value'] is not None for x in a))
for x in a:
 row=base.copy();row.update(analysis_type='Tumor_vs_'+x['reference'],n=x['patients'],effect_type=x['metric'],effect=fmt(x['mean_patient_effect']),p_value=fmt(x['p_value']),status=x['status'],reason='Equal-weight mean across eligible sections within patient; two-sided sign test if >=3 nonzero patient effects; reference is not normal epithelium',area_fraction=x['purity'],positive_patients=x['positive_patients'],n_sections=x['sections']);rows.append(row)
row=base.copy();row.update(analysis_type='Tumor_vs_normal_epithelium',status='NOT_EVALUABLE',reason=s['normal_epithelium_reason'],area_fraction='NA',positive_patients='NA',n_sections='NA');rows.append(row)
with (R/'results.tsv').open('w',encoding='utf8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=cols+['area_fraction','positive_patients','n_sections'],delimiter='\t');w.writeheader();w.writerows(rows)
spec={'version':'cta_v1','primary_area_fraction':.8,'sensitivity_area_fraction':1.,'area_sampling':'49 evenly spaced disk grid points within spot; excludes overlapping labels','QC':{'min_UMI':500,'min_genes':200,'max_mitochondrial_fraction':.25},'minimum_spots_per_region_per_section':20,'normalization':'log1p(count / all_gene_UMI * 10000), no imputation','patient_aggregation':'equal mean of eligible section effects; author patient grouping','tests':'two-sided exact sign test across >=3 nonzero patient effects; no spot-level P; nominal P, no FDR per user preference','registration':'H&E SIFT + affine RANSAC; >=20 inliers, >=35% inlier fraction, median error <=1.5 lowres pixels, second-target <65%; reciprocal source best >2x runner-up; explicit capture suffix exact','registration_seed':42,'analysis_code_commit':'6e29956','initial_registration_commit':'918874f','source_residency':'server165 only'}
(R/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
(R/'code_commit.txt').write_text('Analysis: 6e29956\nInitial registration diagnostics: 918874f\n')
readme=f'''# LYPLA1：CTA乳腺癌病理区域空间分析

## 本轮问题
检验 LYPLA1 是否在作者病理标注的癌区富集，并审计能否回答癌上皮高于正常上皮。

## 输入与范围
[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12420830/)研究包括4位患者23张切片。本次[Zenodo公开包](https://zenodo.org/records/15211538)实取14份表达矩阵、23份人工病理SVG；实际可靠匹配并分析{s['sections_analyzed']}张切片、{s['patients_analyzed']}位患者。未把研究总样本数冒充本次分析样本数。未读取患者级数据到本机。

共{s['spots_input']:,}个输入spot，QC后{s['spots_qc']:,}个，其中{s['LYPLA1_positive_qc']:,}个检出LYPLA1（{s['LYPLA1_positive_qc']/s['spots_qc']:.1%}）。每份矩阵均唯一匹配LYPLA1。数据是常规Visium，不是HD；spot不是纯单细胞。

## 实际结果
{table}
癌区相对免疫区，三个患者在三种指标中均为正方向；相对未标注混合区域，校正总计数及线粒体比例后只有2/3患者为正。校正后效应明显缩小，提示检测深度/组织组成影响不可忽略。没有患者层面P<0.05的结果。

正值表示癌区更高。P值来自患者层面双侧精确符号检验；每位患者内对合格切片效应等权平均。单张切片每组至少20个QC spot。q未计算；主阈值80%采样面积，100%阈值敏感性结果完整保留在results.tsv。区域汇总比值以log2表示，不是原始倍数。

## 新手解释
红线是真正的作者人工癌区标注；黄线为免疫区；蓝线为原位癌。中图为LYPLA1表达，右图表示经过spot面积筛选的区域。未标注区域包含混合组织，不能称正常上皮或纯间质。左右图用于定位，患者汇总用于评价重复方向。

## 限制/反证
正常上皮比较为NOT_EVALUABLE：作者人工标签没有正常上皮类别。三个患者即使方向全部一致，双侧符号检验最小P值仍为0.25；不存在因spot多而提高患者重复数的做法。癌区spot仍含其他细胞。患者数有限，且不同区域的细胞密度与检测深度会影响结果，因此同时报告区域汇总表达与深度/线粒体比例校正。相邻切片可产生少量图像匹配，已通过采集区后缀和双向最佳匹配排除。

## 当前决定
作为带病理标签的癌区定位及区域表达补充；不将其升级为癌上皮高于正常上皮的验证。

## 下一步
如需直接正常上皮比较，应取得明确标注正常导管/小叶且有足够覆盖的独立数据。现有结果与之前队列逐项对照，不按结果正负筛除切片。

## 复现命令
本地运行 `python code/brca_cta_fetch_v1.py`，仅通过内存传输到server165；同一自身运行恢复可用`--resume-own-run`。复制注册和分析脚本到该运行目录后，服务器运行 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/xuzx/miniconda3/envs/ov/bin/python brca_cta_register_v1.py`，再运行 `brca_cta_analyze_v1.py`。初始图像注册诊断采用918874f；正式区域分析及双向匹配选择采用6e29956。所有{v['source_sha256_verified']}个下载文件SHA256已核对。配准质控与视觉叠加记录见validation.json。展示脚本 `code/brca_cta_present_v1.py`。
'''
(R/'README_CN.md').write_text(readme,encoding='utf8')
(R/'delivery_status.json').write_text(json.dumps({'pathology_region_analysis':'DONE','normal_epithelium_comparison':'NOT_EVALUABLE','github_upload':'NOT_RUN'},indent=2))
trs=''.join('<tr>'+''.join('<td>'+str(z)+'</td>' for z in [refnames[x['reference']],metrics[x['metric']],f"{x['positive_patients']}/{x['patients']}",fmt(x['mean_patient_effect']),fmt(x['p_value'])])+'</tr>' for x in primary)
figs=sorted(p for p in G.glob('V*.png'))
sections=''.join(f'<details><summary>{html.escape(p.stem)} — 点击展开</summary><a href="{p.name}" target="_blank"><img loading="lazy" src="{p.name}"></a></details>' for p in figs)
page=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>LYPLA1 · CTA病理区域</title><style>body{{font:17px/1.8 system-ui,sans-serif;color:#19304b;background:#f4f7fb;max-width:1350px;margin:auto;padding:30px}}section,details{{background:white;border:1px solid #d9e1ec;border-radius:12px;padding:22px;margin:22px 0}}h1{{font-size:30px}}img{{width:100%}}table{{border-collapse:collapse;width:100%}}th,td{{padding:10px;border-bottom:1px solid #ddd;text-align:left}}.note{{background:#fff2d7;padding:18px}}summary{{cursor:pointer;font-weight:600}}a{{color:#1763a5}}</style><h1>LYPLA1：带病理标注的乳腺癌空间分析</h1><p>CTA 2025 · 常规Visium · {s['sections_analyzed']}张切片 / {s['patients_analyzed']}位患者 · {s['spots_qc']:,}个合格spot</p><p class="note"><b>本次有真实人工癌区标注，但没有正常上皮对照。</b>未标注的混合区域不能等同正常上皮或纯间质。只有{s['patients_analyzed']}位患者，不能将spot数量当患者重复数。</p><section><h2>患者层面的结果</h2><table><tr><th>癌区对照</th><th>指标</th><th>方向为正的患者</th><th>平均效应</th><th>双侧P</th></tr>{trs}</table><p>每位患者内切片等权汇总；正值表示癌区较高。log2比值不是原始倍数。P为名义P值，未做FDR；至少3个非零患者效应才计算。</p><img src="patient_summary.png"><p>每个点代表一位患者，黑线为患者均值。图中三种效应尺度不同，不能直接比较数值大小。</p></section><section><h2>怎样看切片图</h2><p><b>左图：</b>原始组织图叠加作者病理区域，红色癌区、黄色免疫区、蓝色原位癌。<b>中图：</b>LYPLA1标准化表达，颜色越亮表达越高，零值是未检出。<b>右图：</b>至少80%采样spot面积属于同一区域才纳入该区域；灰色未标注区域是混合组织。各切片色标独立，不能凭颜色跨切片比较强弱。</p><p>可点击图片放大。病理图像配准使用图像特征，未用LYPLA1表达来调整区域。</p></section>{sections}<section><h2>来源与范围</h2><p>论文4位患者23张切片，当前公开包有14份矩阵、23份人工标注。本次只使用可靠对应的{s['sections_analyzed']}张切片。LYPLA1检出率为{s['LYPLA1_positive_qc']/s['spots_qc']:.1%}，这是spot检出率，不是阳性癌细胞比例。</p><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC12420830/">论文</a> · <a href="https://zenodo.org/records/15211538">公开数据</a></section></html>'''
G.mkdir(parents=True,exist_ok=True);(G/'index.html').write_text(page,encoding='utf8')
index=ROOT/'coordination/stages/BRCA.tsv'
with index.open(encoding='utf8',newline='') as f:w=csv.DictReader(f,delimiter='\t');fields=w.fieldnames;rr=list(w)
rr=[x for x in rr if not(x['run_id']==RUN and x['stage_id']=='06_EXTERNAL')]
rr.append(dict(zip(fields,['BRCA','06_EXTERNAL',RUN,'cta_v1','PARTIAL',f"CTA manual pathology region analysis: {s['sections_analyzed']} sections / {s['patients_analyzed']} patients",R.relative_to(ROOT).as_posix(),'code/brca_cta_analyze_v1.py','analysis/brca-lypla1-spatial-20260927','Region analysis DONE; normal epithelium NOT_EVALUABLE','Independent annotated normal epithelial comparison'])))
rr.sort(key=lambda x:(x['stage_id'],x['run_id']))
with index.open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rr)
print('CTA public report and gallery built')
