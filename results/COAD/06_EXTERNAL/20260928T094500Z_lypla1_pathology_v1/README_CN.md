# LYPLA1空间表达：加入作者病理区域标注

## 本轮问题
用户指出旧表达图无法看出癌区。这一意见成立：GSE283052既有图没有取得逐点病理标签，只能说明表达位置，不能据此判断癌区富集。本批使用另一项具有公开病理标签的CRC研究补上区域信息，不覆盖旧结果。

## 输入与范围
[Valdeolivas等原研究](https://www.nature.com/articles/s41698-023-00488-4)；[作者公开数据Zenodo 7760264](https://zenodo.org/records/7760264)。原研究7病例、14张切片；本批在查看LYPLA1数值前沿用作者Methods中组织分级选择的每病例一张切片，共7张。结肠与直肠均有，不能称为纯COAD；S7主要非肿瘤组织。

病理标签直接来自作者Pathology_SpotAnnotations.zip，按唯一Barcode连接表达矩阵和原始坐标。tumor为肿瘤区域；tumor&stroma单列混合区域；非肿瘤上皮、基质及其他组织分别保存。完整细分标签见annotation_dictionary.tsv，不从LYPLA1、EPCAM或聚类推断癌区。

原空间计数与图像只保存在server165。受下载速度限制，按HTTP Range选择读取每个源ZIP的filtered H5、坐标、缩放参数和低分辨率H&E；逐成员源CRC32、SHA256核验通过。没有完整下载7个源ZIP，不能声称它们的官方整包MD5已核验；病理标签整包MD5已核验。

## 实际结果
7张切片共9,465个点进入区域表达描述。保留作者exclude、无标签和数值QC失败点的去向；这些点不进入区域均值。数字是本批QC后的数量，不等于论文全14张切片总数。

下表为各切片区域内平均log1p(CP10K)，括号内是达标采集点数；不是倍数变化，也不是独立患者数。肿瘤与基质混合区没有并入肿瘤列，完整结果另表保存。

|切片|病理肿瘤区|基质区|非肿瘤上皮|
|---|---:|---:|---:|
|S1_Cec_Rep1|0.743（375点）|0.758（70点）|不可评估|
|S2_Col_R_Rep1|0.582（233点）|0.520（209点）|0.663（71点）|
|S3_Col_R_Rep1|0.569（669点）|0.533（721点）|不可评估|
|S4_Col_Sig_Rep1|0.866（43点）|0.565（156点）|不可评估|
|S5_Rec_Rep1|0.729（337点）|0.654（249点）|0.496（293点）|
|S6_Rec_Rep2|0.616（163点）|0.422（959点）|0.425（59点）|
|S7_RecSig_Rep1|不可评估|不可评估|0.607（82点）|

以每个区域至少20个达标点作为描述比较门槛：6张切片可作肿瘤—基质比较，其中5张肿瘤区域均值较高；3张可作肿瘤—非肿瘤上皮比较，其中2张肿瘤区域均值较高。这只是切片内描述，不是显著性检验或跨患者总体结论。未达到条件的比较保留NOT_EVALUABLE。

## 新手解释
- 每张三联图从左到右：原始H&E、作者病理区域、LYPLA1表达。
- 病理图红色=肿瘤，橙色=肿瘤与基质混合，蓝色=基质，绿色=非肿瘤上皮，紫色=其他有注释组织，灰色=排除/缺标签/QC未过。
- LYPLA1图由浅灰至深红表示表达增加，所有切片共享0至2.6273色标。青色圆圈是作者标注的肿瘤采集点；圆圈不是推断出的连续肿瘤边界。只圈入本批达标点。
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
