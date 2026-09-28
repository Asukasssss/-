"""Public aggregate report and local gallery; never copies bin measurements."""
from pathlib import Path
import csv,json,html
ROOT=Path(__file__).resolve().parents[1];RUN='20260928T040500Z_hd3_lypla1_v1'
R=ROOT/'results/BRCA/06_EXTERNAL'/RUN
G=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/LYPLA1_HD3_spatial')
coverage=list(csv.DictReader((R/'coverage.tsv').open(),delimiter='\t'))
assoc=list(csv.DictReader((R/'marker_associations.tsv').open(),delimiter='\t'))
rows=[];cols=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t')
base=dict.fromkeys(cols,'NA');base.update(cancer='BRCA',cohort='PacBio_HD3_Kinnex_humanBreast',stage_id='06_EXTERNAL',run_id=RUN,analysis_version='hd3_v1',gene='LYPLA1',source_id='PacBio public HD3-Kinnex humanBreast',unit='specimen',n=1,test_family='descriptive single-specimen HD exploration',family_n_evaluable=0)
for x in coverage:
 row=base.copy();row.update(analysis_type='LYPLA1_detection_coverage',effect_type='fraction_of_QC_bins_detected',effect=x['detected_fraction'],status='DONE',reason='Descriptive coverage only; bin size affects detection probability',bin_um=x['bin_um'],n_qc_bins=x['n_qc_bins']);rows.append(row)
row=base.copy();row.update(analysis_type='cancer_vs_normal_epithelium',status='NOT_EVALUABLE',reason='No verified pathological cancer/normal epithelial region labels; one donor',bin_um=16,n_qc_bins='NA');rows.append(row)
with (R/'results.tsv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=cols+['bin_um','n_qc_bins'],delimiter='\t');w.writeheader();w.writerows(rows)
table='|网格|合格网格|LYPLA1 检出网格|检出率|\n|---|---:|---:|---:|\n'+''.join(f"|{x['bin_um']} μm|{int(x['n_qc_bins']):,}|{int(x['n_detected_bins']):,}|{float(x['detected_fraction']):.1%}|\n" for x in coverage)
readme=f'''# LYPLA1：乳腺癌 Visium HD 3′ 空间探索

## 本轮问题
寻找 HD 空间转录组并实际分析 LYPLA1 的空间表达；优先保留癌区/正常上皮比较所需的证据边界。

## 输入与范围
采用 [PacBio 公开数据](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/DATA-RevioSPRQ-Kinnex-VisiumHD-humanBreast/)，Visium HD 3′ 建库联合 Kinnex 长读长测序，供者为一位70岁女性，来源说明为浸润性导管癌、新鲜冰冻组织。技术来源：[厂商应用说明](https://www.pacb.com/wp-content/uploads/Application-note-Kinnex-single-cell-RNA-kit-for-10x-Visium-HD-3-Spatial-Gene-Expression.pdf)。分析使用基因级 Kinnex 处理矩阵，不将报告中的短读长流程统计冒充本矩阵统计。

输入为22,854个特征、3,672,073个2 μm空间条码、18,721,468个非零项。原始 LYPLA1 特征计数共11,408；另外 LYPLA1+TCEA1 联合特征的1个计数未并入目标基因。

另核查10x常用探针型 HD 乳腺 DCIS 数据：LYPLA1 在参考表被标为 included=FALSE，官方 feature_slice 的 target_sets 也未保留它。即使原始层有读数，也未用其作可靠表达证据，见 probe_coverage_audit.json。此结论针对核查的v2探针及数据，不推及所有HD技术。

## 实际结果
{table}
16 μm为预设主展示尺度，8与32 μm为敏感性展示。网格增大导致检出率上升，不表示基因被进一步上调。

LYPLA1 与上皮标记均值的描述性 Spearman 相关：8 μm为0.126，16 μm为0.257，32 μm为0.363；控制总计数的偏秩相关分别降为0.019、0.042、0.057。完整结果见 marker_associations.tsv。

空间覆盖及可视化完成；正式“癌上皮高于正常上皮”仍为 NOT_EVALUABLE，没有基于网格数量计算显著性 P 值。

## 新手解释
这次确实在HD数据中测到了LYPLA1。较粗尺度能让分布更容易看清，但其与上皮标记的同步变化很大程度随测序深度变化，校正后关联较弱。不能把这些图解释为癌细胞特异富集或癌上皮高于正常上皮。

## 限制/反证
仅一份标本；公开下载包和web_summary未提供明确病理癌区/正常导管标签。EPCAM/KRT等标记不能区分恶性与正常上皮。HD的2 μm条码及8/16/32 μm网格都不是已分割的单细胞；稀疏零值不能当作生物学不表达。上皮关联经深度校正后减弱是必须保留的限制。

## 当前决定
保留为单标本HD空间探索，不升级为正常上皮对照验证。没有继续强制用被过滤探针作阳性图。

## 下一步
若要回答癌上皮高于正常上皮，需获得这一标本可追溯的病理区域标注及足够正常导管，并在更多独立供者重复。当前图可用于查看位置和规划后续标注。

## 复现命令
本地 `python code/brca_hd_fetch_v1.py` 将源文件经RAM传到server165，不保存本地源矩阵。服务器执行 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/xuzx/miniconda3/envs/ov/bin/python brca_hd3_analyze_v1.py RUN_DIR`。随后用 `brca_hd3_zoom_v1.py RUN_DIR` 从保存的16 μm结果绘图，不重跑统计。分析代码提交 bde0827。源文件 SHA256全部核对一致。

图像对齐采用来源web_summary的变换矩阵；本矩阵8 μm总计数与提供方UMI图的秩相关为0.961，组织掩膜427,917个点与原报告完全一致。坐标变换后的QC点均在图像范围内，另完成图像叠加目视检查。
'''
(R/'README_CN.md').write_text(readme,encoding='utf8')
(R/'code_commit.txt').write_text('bde0827\n')
(R/'delivery_status.json').write_text(json.dumps({'spatial_exploration':'DONE','cancer_vs_normal_epithelium':'NOT_EVALUABLE','github_upload':'NOT_RUN'},indent=2))
trs=''.join('<tr>'+''.join(f'<td>{v}</td>' for v in [x['bin_um']+' μm',f"{int(x['n_qc_bins']):,}",f"{int(x['n_detected_bins']):,}",f"{float(x['detected_fraction']):.1%}"])+'</tr>' for x in coverage)
sections=''.join(f'<section><h2>{title}</h2><p>{text}</p><a href="{file}" target="_blank"><img src="{file}" alt="{title}"></a></section>' for file,title,text in [
 ('he_zoom_16um.png','1. 组织图像与 LYPLA1 放大图','左：H&E；中：LYPLA1；右：EPCAM/KRT8/KRT18/KRT19上皮标记均值。标记图不是癌区注释。颜色表示log1p(CP10K)，每幅图各有色标；颜色不能跨基因直接比高低。'),
 ('maps_16um.png','2. 主结果：16 μm空间图','看LYPLA1的位置，再与上皮、间质和免疫标记对照。白色为组织外或未通过QC；黑色为该网格没有检测到。稀疏零值不等同于不表达。'),
 ('maps_8um.png','3. 更细：8 μm','细网格检出率只有6.3%，不适合凭零值判断细胞阴性。'),
 ('maps_32um.png','4. 更粗：32 μm','检出率升至36.9%，属于聚合后的检测变化，不意味着LYPLA1表达上调。'),
 ('coordinate_qc.png','5. 坐标与组织掩膜核查','提供方UMI图、本次聚合计数和组织掩膜。图形布局对应，8 μm总计数与提供方UMI图的秩相关为0.961。')])
page=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>LYPLA1 乳腺癌HD3空间结果</title><style>body{{font:17px/1.8 system-ui,sans-serif;color:#19304b;background:#f4f7fb;max-width:1250px;margin:auto;padding:30px}}section{{background:white;border:1px solid #d9e1ec;border-radius:12px;padding:24px;margin:24px 0}}h1{{font-size:29px}}img{{width:100%}}table{{border-collapse:collapse;width:100%}}th,td{{padding:12px;border-bottom:1px solid #ddd;text-align:left}}.note{{background:#fff2d7;padding:18px}}a{{color:#1763a5}}</style><h1>LYPLA1：乳腺癌 Visium HD 3′ 实际结果</h1><p>PacBio公开HD 3′ + Kinnex长读长数据；一份浸润性导管癌标本。已完成空间表达探索。</p><p class="note"><b>测到了LYPLA1，但不能据此证明癌上皮高于正常上皮。</b>缺少病理区域标签；与上皮标记的关联经测序深度校正后较弱。</p><section><h2>覆盖与信号</h2><table><tr><th>网格</th><th>合格网格数</th><th>LYPLA1检出网格</th><th>检出率</th></tr>{trs}</table><p>16 μm主尺度：LYPLA1与上皮标记均值相关为0.257，校正总计数后为0.042。32 μm分别为0.363和0.057。没有把网格当患者计算P值。</p></section>{sections}<section><h2>为何没有使用常见的探针型HD数据</h2><p>核查的10x探针型乳腺DCIS数据中，LYPLA1被v2参考探针表默认排除，官方保留目标集也没有它。本次改用实际覆盖LYPLA1的HD 3′/Kinnex公开矩阵。没有将LYPLA1+TCEA1联合特征混入LYPLA1。</p><p><a href="https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/DATA-RevioSPRQ-Kinnex-VisiumHD-humanBreast/">数据来源</a> · <a href="https://www.pacb.com/wp-content/uploads/Application-note-Kinnex-single-cell-RNA-kit-for-10x-Visium-HD-3-Spatial-Gene-Expression.pdf">技术说明</a></p></section></html>'''
G.mkdir(parents=True,exist_ok=True);(G/'index.html').write_text(page,encoding='utf8')
index=ROOT/'coordination/stages/BRCA.tsv'
with index.open(encoding='utf8',newline='') as f:w=csv.DictReader(f,delimiter='\t');fields=w.fieldnames;rr=list(w)
rr=[x for x in rr if not(x['run_id']==RUN and x['stage_id']=='06_EXTERNAL')]
rr.append(dict(zip(fields,['BRCA','06_EXTERNAL',RUN,'hd3_v1','PARTIAL','One-sample HD3 Kinnex spatial exploration; probe-based HD LYPLA1 excluded',R.relative_to(ROOT).as_posix(),'code/brca_hd3_analyze_v1.py','analysis/brca-lypla1-spatial-20260927','Spatial plots DONE; pathology cancer-normal comparison NOT_EVALUABLE','Obtain verified pathological regions and independent donors'])))
rr.sort(key=lambda x:(x['stage_id'],x['run_id']))
with index.open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rr)
print('Public report and gallery complete')
