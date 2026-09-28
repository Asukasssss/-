"""Build public Chinese report from server-produced aggregate results."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RUN='20260928T094500Z_lypla1_pathology_v1'
REL=f'results/COAD/06_EXTERNAL/{RUN}'
P=ROOT/REL
def read(name):
 with (P/name).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
v=json.loads((P/'validation.json').read_text())
coverage=read('coverage.tsv');summary=read('region_summary.tsv');comp=read('region_comparisons.tsv')
rows=[]
for c in coverage:
 s=c['section'];d={r['region']:r for r in summary if r['section']==s and r['level']=='broad'}
 def cell(k):
  r=d.get(k)
  return '不可评估' if not r or not int(r['admitted_spots']) else f"{float(r['mean_log1pCP10K']):.3f}（{r['admitted_spots']}点）"
 rows.append('|'+s+'|'+cell('Tumor')+'|'+cell('Stroma')+'|'+cell('Non-neoplastic epithelium')+'|')
eligible=[r for r in comp if r['status']=='DONE' and r['reference_region']=='Stroma']
up=sum(float(r['effect'])>0 for r in eligible)
normal=[r for r in comp if r['status']=='DONE' and r['reference_region']=='Non-neoplastic epithelium']
upn=sum(float(r['effect'])>0 for r in normal)
text=f'''# LYPLA1空间表达：加入作者病理区域标注

## 本轮问题
用户指出旧表达图无法看出癌区。这一意见成立：GSE283052既有图没有取得逐点病理标签，只能说明表达位置，不能据此判断癌区富集。本批使用另一项具有公开病理标签的CRC研究补上区域信息，不覆盖旧结果。

## 输入与范围
[Valdeolivas等原研究](https://www.nature.com/articles/s41698-023-00488-4)；[作者公开数据Zenodo 7760264](https://zenodo.org/records/7760264)。原研究7病例、14张切片；本批在查看LYPLA1数值前沿用作者Methods中组织分级选择的每病例一张切片，共7张。结肠与直肠均有，不能称为纯COAD；S7主要非肿瘤组织。

病理标签直接来自作者Pathology_SpotAnnotations.zip，按唯一Barcode连接表达矩阵和原始坐标。tumor为肿瘤区域；tumor&stroma单列混合区域；非肿瘤上皮、基质及其他组织分别保存。完整细分标签见annotation_dictionary.tsv，不从LYPLA1、EPCAM或聚类推断癌区。

原空间计数与图像只保存在server165。受下载速度限制，按HTTP Range选择读取每个源ZIP的filtered H5、坐标、缩放参数和低分辨率H&E；逐成员源CRC32、SHA256核验通过。没有完整下载7个源ZIP，不能声称它们的官方整包MD5已核验；病理标签整包MD5已核验。

## 实际结果
7张切片共{v['admitted_spots']:,}个点进入区域表达描述。保留作者exclude、无标签和数值QC失败点的去向；这些点不进入区域均值。数字是本批QC后的数量，不等于论文全14张切片总数。

下表为各切片区域内平均log1p(CP10K)，括号内是达标采集点数；不是倍数变化，也不是独立患者数。肿瘤与基质混合区没有并入肿瘤列，完整结果另表保存。

|切片|病理肿瘤区|基质区|非肿瘤上皮|
|---|---:|---:|---:|
{chr(10).join(rows)}

以每个区域至少20个达标点作为描述比较门槛：{len(eligible)}张切片可作肿瘤—基质比较，其中{up}张肿瘤区域均值较高；{len(normal)}张可作肿瘤—非肿瘤上皮比较，其中{upn}张肿瘤区域均值较高。这只是切片内描述，不是显著性检验或跨患者总体结论。未达到条件的比较保留NOT_EVALUABLE。

## 新手解释
- 每张三联图从左到右：原始H&E、作者病理区域、LYPLA1表达。
- 病理图红色=肿瘤，橙色=肿瘤与基质混合，蓝色=基质，绿色=非肿瘤上皮，紫色=其他有注释组织，灰色=排除/缺标签/QC未过。
- LYPLA1图由浅灰至深红表示表达增加，所有切片共享0至{v['shared_expression_range'][1]:.4f}色标。青色圆圈是作者标注的肿瘤采集点；圆圈不是推断出的连续肿瘤边界。只圈入本批达标点。
- 相同坐标使读者能判断红色高表达是否出现在癌区；病理分区颜色与表达颜色含义不同。无平滑、无插值、无表达值截断。
- 总览见[7张切片总览](figures/all7_pathology_LYPLA1.png)，逐切片PNG/PDF均在figures。

## 限制与反证
Visium每点包含多个细胞，病理肿瘤点也可能含免疫和基质细胞，不能将信号全部归于癌细胞。肿瘤—基质差异也不等于恶性—正常上皮差异。S7没有可用病理tumor组时保留缺项，不能把推断少量肿瘤细胞补作新标签。此队列与CAMP的患者独立性未另行认证。

沿用来源数值QC（UMI 500–45000；线粒体比例≤0.5）和作者排除标签，另要求组织内、标签存在。本批表达尺度为log1p(10000×LYPLA1 UMI/总Gene Expression UMI)，不是作者SCTransform输出；点深度和细胞组成仍可能影响比较。无新P/q，不把成千上万个点当独立患者。

## 当前决定
旧图保留为无病理标签的空间表达背景；本批新增明确癌区位置的表达图和逐切片区域描述。是否癌区更高以实际区域数值为准，不预设统一富集，不改变CAMP或既有单细胞结论。

## 下一步
优先阅读有肿瘤、基质、非肿瘤上皮三类覆盖的切片。需要严格比较时再另定以病例为单位、考虑细胞组成和深度的设计，本批不自动扩展。

## 复现命令
从作者Zenodo按code/coad/relay_lypla1_pathology_selected.py选取必要成员并写入165；脚本不保存密码或本地源矩阵。在新的运行目录创建.running锁，复制acquire_lypla1_pathology_v1.py及analyze_lypla1_pathology_v1.py后运行：

```bash
python3 analyze_lypla1_pathology_v1.py --out <new_COAD_B_run_dir> --commit <code_base_commit>
python3 render_lypla1_pathology_overview_v1.py --out <new_COAD_B_run_dir>
```

图源逐点测量留165运行目录private；公开region_summary.tsv、region_comparisons.tsv、coverage.tsv、annotation_dictionary.tsv、analysis_spec.json、source_manifest.tsv与validation.json。公开图件为生成的科学图，不含可下载的源矩阵。
'''
(P/'README_CN.md').write_bytes(text.encode('utf-8'))
status=ROOT/'reports/COAD/CURRENT_STATUS_CN.md'
status.write_bytes((f'# 最新：LYPLA1已加入作者病理癌区标注\n\n[病理与表达对照报告](../../{REL}/README_CN.md)：Valdeolivas CRC研究每病例一张切片，共7张、{v["admitted_spots"]:,}个达标空间点。病理分区和LYPLA1同坐标，表达图青圈标识tumor点；混合区、基质、非肿瘤上皮分别保留。肿瘤—基质{len(eligible)}张可描述、{up}张肿瘤均值较高；肿瘤—非肿瘤上皮{len(normal)}张可描述、{upn}张肿瘤均值较高。无新P/q，不将多细胞点当癌细胞或独立患者。前批无标签图不作癌区富集依据。\n\n'+status.read_text(encoding='utf-8')).encode('utf-8'))
index=ROOT/'coordination/stages/COAD.tsv'
with index.open(encoding='utf-8',newline='') as f:
 reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
assert not any(r['run_id']==RUN for r in rows)
rows.append(dict(zip(fields,['COAD','06_EXTERNAL',RUN,'COAD_LYPLA1_pathology_v1','DONE',f'7 author-selected CRC sections;{v["admitted_spots"]} admitted spots;pathologist regions plus LYPLA1 maps',REL,'code/coad/analyze_lypla1_pathology_v1.py','analysis/coad-initial','Exact barcodes;source pathology;CRC including rectum;no P/q;no inferred tumor labels','Interpret region descriptions with multicellular spot and sampling limits'])))
with index.open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sorted(rows,key=lambda r:(r['stage_id'],r['run_id'])))
print('REPORT_DONE')
