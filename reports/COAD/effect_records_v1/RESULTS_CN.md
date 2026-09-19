# COAD：从冻结效应记录开始的描述性分析

账号 B。只读取癌种效应表，不读取既有候选名单或代谢物—基因映射，不重算或更改原效应与 q。

共 159 条记录，159 个代谢物键；原 q<0.05 的有 73 条，其余 86 条保留为背景。
显著记录效应方向：{"positive": 33, "negative": 40}。按项目既有肿瘤减正常解释，正值为肿瘤较高、负值为较低；本轮未回查源矩阵。g 是标准化差异，不是倍数。
dataset 字段：{"COAD": 159}；source_type 字段：{"single_dataset": 159}。单队列记录不构成独立队列复现。

# 原注释与数据限制

原身份标签（全部）：{"HIGH_ID": 155, "MODERATE_NAME_ONLY": 4}；显著记录：{"HIGH_ID": 73}。HIGH_ID 仅为沿用标签，不是重新鉴定。
显著记录中无 KEGG/HMDB 标识：0 条。存在标识也不保证异构体、反应底物或峰身份已确定。
原 max_input_raw_missing_rate 范围为 [0.0, 0.794871794871795]。其中 >50% 的共有 9 条，显著记录中 6 条。50% 仅作描述性计数，未用于剔除或改阈值；该字段不能代替实际样本量。
表中没有样本数、患者身份、配对关系或原始 P 值，不能由本表核实设计、重新校正 q 或建立患者关联。

# 通路注释分布（计数，不是富集分析）

| 原大类 | 全部记录 | 原 q<0.05 | 正效应显著 | 负效应显著 |
|---|---:|---:|---:|---:|
| 未注释 | 26 | 8 | 3 | 5 |
| Amino Acid | 17 | 10 | 5 | 5 |
| Amino acid | 39 | 18 | 12 | 6 |
| Carbohydrate | 19 | 10 | 3 | 7 |
| Cofactors and Vitamins | 5 | 1 | 1 | 0 |
| Energy | 6 | 3 | 0 | 3 |
| Lipid | 19 | 5 | 2 | 3 |
| Nucleotide | 24 | 15 | 6 | 9 |
| Peptide | 3 | 2 | 1 | 1 |
| Xenobiotics | 1 | 1 | 0 | 1 |

原表 Amino Acid 与 Amino acid 分别保留显示；仅按大小写统一后的描述性合计为 28 条显著记录。此合计不是通路富集或通量结论。

# 显著但原最大缺失率超过 50% 的记录

| 特征（原名） | Hedges g | 原 q | 原最大缺失率 | 原身份标签 |
|---|---:|---:|---:|---|
| N-acetylmethionine | 0.760 | 0.00352 | 59.0% | HIGH_ID |
| 5-hydroxylysine | 1.102 | 0.0001065 | 53.8% | HIGH_ID |
| isocitrate | -0.940 | 0.0007013 | 67.6% | HIGH_ID |
| 4-methyl-2-oxopentanoate | -1.475 | 2.228e-06 | 51.4% | HIGH_ID |
| 2-phosphoglycerate | -1.445 | 2.185e-06 | 64.9% | HIGH_ID |
| serotonin (5HT) | -2.577 | 5.92e-10 | 54.1% | HIGH_ID |

名称 Disulfiram 的记录原 q 很小，但大类/子类注释均为空，下一步需回查原注释与身份来源。本轮不更改其名称或标识，也不据此给出基因或干预结论。

# 显著记录：按原 q 排序

此排序用于阅读，不是靶点优先级。所有 73 条均保留，不用 g 或缺失率另设筛选门槛。

| 特征（原名） | Hedges g | 原 q | 原最大缺失率 | 原身份标签 |
|---|---:|---:|---:|---|
| beta-alanine | 2.226 | 3.399e-10 | 0.0% | HIGH_ID |
| serotonin (5HT) | -2.577 | 5.92e-10 | 54.1% | HIGH_ID |
| glucose | -2.468 | 7.024e-10 | 10.8% | HIGH_ID |
| Disulfiram | -2.024 | 1.467e-08 | 0.0% | HIGH_ID |
| hypotaurine | 1.358 | 3.777e-08 | 2.7% | HIGH_ID |
| proline | 1.832 | 3.777e-08 | 0.0% | HIGH_ID |
| taurine | 1.768 | 3.777e-08 | 0.0% | HIGH_ID |
| N-acetylaspartate (NAA) | -1.609 | 4.01e-08 | 0.0% | HIGH_ID |
| S-adenosylmethionine (SAM) | 1.587 | 9.785e-08 | 5.1% | HIGH_ID |
| 3PG | -1.622 | 2.5e-07 | 0.0% | HIGH_ID |
| cytidine 5'-diphosphocholine | 1.203 | 7.433e-07 | 10.3% | HIGH_ID |
| creatinine | -1.216 | 2.154e-06 | 0.0% | HIGH_ID |
| ribulose/xylulose 5-phosphate | -1.389 | 2.185e-06 | 0.0% | HIGH_ID |
| 2-phosphoglycerate | -1.445 | 2.185e-06 | 64.9% | HIGH_ID |
| 4-methyl-2-oxopentanoate | -1.475 | 2.228e-06 | 51.4% | HIGH_ID |
| allantoin | -1.133 | 2.759e-06 | 16.2% | HIGH_ID |
| urate | -1.382 | 3.037e-06 | 0.0% | HIGH_ID |
| acetylcholine | -1.264 | 3.037e-06 | 27.0% | HIGH_ID |
| adenosine | -1.434 | 3.684e-06 | 0.0% | HIGH_ID |
| histamine | -1.217 | 4.589e-06 | 0.0% | HIGH_ID |
| 3-hydroxy-3-methylglutarate | 1.283 | 4.645e-06 | 15.4% | HIGH_ID |
| citrate | -1.082 | 4.645e-06 | 0.0% | HIGH_ID |
| glycylglycine | 1.094 | 5.325e-06 | 10.3% | HIGH_ID |
| N-acetylglucosamine 6-phosphate | 1.206 | 1.513e-05 | 5.1% | HIGH_ID |
| Acetohydroxamate | -0.521 | 1.513e-05 | 10.3% | HIGH_ID |
| glutamine | -1.185 | 1.532e-05 | 0.0% | HIGH_ID |
| hypoxanthine | 1.134 | 7.858e-05 | 0.0% | HIGH_ID |
| guanosine | -1.106 | 7.945e-05 | 0.0% | HIGH_ID |
| 5-hydroxylysine | 1.102 | 0.0001065 | 53.8% | HIGH_ID |
| adenine | 1.006 | 0.0001719 | 15.4% | HIGH_ID |
| Mucate | -1.073 | 0.0001719 | 8.1% | HIGH_ID |
| kynurenine | 1.036 | 0.0002183 | 12.8% | HIGH_ID |
| UDP-N-acetylglucosamine | 0.939 | 0.0002908 | 0.0% | HIGH_ID |
| cis-aconitate | -0.835 | 0.0003078 | 5.4% | HIGH_ID |
| ophthalmate | -0.854 | 0.0003127 | 2.6% | HIGH_ID |
| isocitrate | -0.940 | 0.0007013 | 67.6% | HIGH_ID |
| N-acetylneuraminate | -0.914 | 0.0009277 | 0.0% | HIGH_ID |
| adenylosuccinate | 0.723 | 0.001036 | 10.3% | HIGH_ID |
| inosine 5'-monophosphate (IMP) | -0.889 | 0.001036 | 0.0% | HIGH_ID |
| 5-methylthioadenosine (MTA) | 0.794 | 0.001107 | 35.9% | HIGH_ID |
| uridine | -0.929 | 0.001112 | 0.0% | HIGH_ID |
| glutathione, reduced (GSH) | 0.724 | 0.00138 | 5.1% | HIGH_ID |
| urea | -0.857 | 0.001917 | 0.0% | HIGH_ID |
| glutamate | 0.753 | 0.002053 | 0.0% | HIGH_ID |
| triethanolamine | -0.763 | 0.002053 | 24.3% | HIGH_ID |
| dihydroxyacetone phosphate (DHAP) | -0.921 | 0.00209 | 5.4% | HIGH_ID |
| uracil | 0.735 | 0.003175 | 30.8% | HIGH_ID |
| inosine | -0.467 | 0.003175 | 0.0% | HIGH_ID |
| guanine | -0.590 | 0.003401 | 16.2% | HIGH_ID |
| N-acetylmethionine | 0.760 | 0.00352 | 59.0% | HIGH_ID |
| creatine | -0.803 | 0.00511 | 0.0% | HIGH_ID |
| phosphoenolpyruvate (PEP) | -0.682 | 0.005825 | 24.3% | HIGH_ID |
| guanidinoacetate | -0.523 | 0.007479 | 2.7% | HIGH_ID |
| 1-methylnicotinamide | 0.737 | 0.01043 | 0.0% | HIGH_ID |
| 1-methyladenosine | 0.745 | 0.01043 | 35.9% | HIGH_ID |
| UDP-glucuronate | 0.514 | 0.01043 | 2.7% | HIGH_ID |
| Phthalate | -0.608 | 0.01043 | 29.7% | HIGH_ID |
| glycerophosphorylcholine (GPC) | -0.672 | 0.01131 | 0.0% | HIGH_ID |
| 2-Deoxyglucose 6-phosphate | 0.660 | 0.01268 | 48.7% | HIGH_ID |
| 3-hydroxyproline | 0.612 | 0.01287 | 0.0% | HIGH_ID |
| asparagine | 0.691 | 0.01635 | 0.0% | HIGH_ID |
| N-Acetylglucosamine 1-phosphate | 0.652 | 0.02332 | 0.0% | HIGH_ID |
| 1-methylhistamine | -0.597 | 0.02587 | 8.1% | HIGH_ID |
| azelate (nonanedioate) | -0.383 | 0.02771 | 12.8% | HIGH_ID |
| 6-phosphogluconate | -0.642 | 0.02912 | 5.4% | HIGH_ID |
| cysteine-glutathione disulfide | 0.592 | 0.02954 | 2.6% | HIGH_ID |
| carnosine | -0.552 | 0.03032 | 21.6% | HIGH_ID |
| spermidine | 0.306 | 0.03827 | 0.0% | HIGH_ID |
| leucine | 0.508 | 0.0399 | 0.0% | HIGH_ID |
| tyrosine | 0.544 | 0.0399 | 0.0% | HIGH_ID |
| isoleucine | 0.568 | 0.04461 | 0.0% | HIGH_ID |
| uridine 5'-diphosphate (UDP) | 0.401 | 0.04585 | 2.7% | HIGH_ID |
| cytidine | -0.485 | 0.04585 | 2.7% | HIGH_ID |

# 未达到原 q<0.05 的背景记录

| 特征（原名） | Hedges g | 原 q | 原最大缺失率 | 原身份标签 |
|---|---:|---:|---:|---|
| xanthine | 0.384 | 0.1219 | 12.8% | HIGH_ID |
| guanosine 5'- diphosphate (GDP) | 0.334 | 0.09032 | 28.2% | HIGH_ID |
| ADP-glucose | 0.389 | 0.06528 | 15.4% | HIGH_ID |
| pantothenate | 0.327 | 0.1666 | 30.8% | HIGH_ID |
| 2,3-diphosphoglycerate | 0.245 | 0.3411 | 64.1% | HIGH_ID |
| guanosine 5'- monophosphate (5'-GMP) | 0.177 | 0.09257 | 2.7% | HIGH_ID |
| threonine | 0.480 | 0.06843 | 0.0% | HIGH_ID |
| S-adenosylhomocysteine (SAH) | 0.176 | 0.08075 | 5.4% | HIGH_ID |
| phenylalanine | 0.463 | 0.05398 | 0.0% | HIGH_ID |
| 5-oxoproline | 0.555 | 0.1276 | 0.0% | HIGH_ID |
| deoxycarnitine | 0.402 | 0.1278 | 0.0% | HIGH_ID |
| succinate | 0.338 | 0.1203 | 0.0% | HIGH_ID |
| aspartate | 0.460 | 0.2679 | 0.0% | HIGH_ID |
| methionine | 0.217 | 0.4647 | 0.0% | HIGH_ID |
| carnitine | 0.351 | 0.1176 | 0.0% | HIGH_ID |
| adenosine 2'-monophosphate (2'-AMP) | 0.205 | 0.1729 | 0.0% | HIGH_ID |
| valine | 0.368 | 0.2484 | 0.0% | HIGH_ID |
| N-acetylglucosamine | 0.409 | 0.1339 | 15.4% | HIGH_ID |
| glycine | 0.250 | 0.338 | 0.0% | HIGH_ID |
| 2-hydroxyglutarate | 0.302 | 0.171 | 0.0% | HIGH_ID |
| cadaverine | 0.113 | 0.6041 | 48.7% | HIGH_ID |
| phosphoethanolamine | 0.328 | 0.5602 | 0.0% | HIGH_ID |
| putrescine | 0.295 | 0.3411 | 0.0% | HIGH_ID |
| glycerol 3-phosphate (G3P) | -0.159 | 0.8994 | 0.0% | HIGH_ID |
| tryptophan | 0.131 | 0.4647 | 0.0% | HIGH_ID |
| Butanoate | 0.245 | 0.34 | 30.8% | MODERATE_NAME_ONLY |
| lactate | 0.381 | 0.2252 | 0.0% | HIGH_ID |
| N-acetylglutamate | 0.396 | 0.1276 | 8.1% | HIGH_ID |
| arginine | 0.158 | 0.9215 | 0.0% | HIGH_ID |
| nicotinamide | -0.234 | 0.5209 | 0.0% | HIGH_ID |
| phosphocholine | 0.344 | 0.1355 | 0.0% | HIGH_ID |
| histidine | 0.098 | 0.7508 | 0.0% | HIGH_ID |
| S7P | -0.191 | 0.5419 | 0.0% | MODERATE_NAME_ONLY |
| lysine | 0.046 | 0.919 | 0.0% | HIGH_ID |
| alanine | -0.014 | 0.9324 | 0.0% | HIGH_ID |
| alanylalanine | 0.320 | 0.2327 | 51.3% | HIGH_ID |
| acetylcarnitine | 0.256 | 0.2484 | 0.0% | HIGH_ID |
| dimethylglycine | -0.030 | 0.7065 | 43.6% | HIGH_ID |
| gamma-aminobutyrate (GABA) | 0.545 | 0.08075 | 0.0% | HIGH_ID |
| GDP-mannose | -0.145 | 0.5667 | 27.0% | HIGH_ID |
| choline | -0.031 | 0.9646 | 0.0% | HIGH_ID |
| serine | -0.357 | 0.1599 | 0.0% | HIGH_ID |
| malate | -0.150 | 0.5069 | 0.0% | HIGH_ID |
| fructose 1,6-diphosphate/glucose 1,6-diphosphate/myo-inositol diphosphates | 0.336 | 0.1339 | 2.7% | HIGH_ID |
| glutathione, oxidized (GSSG) | -0.060 | 0.9917 | 0.0% | HIGH_ID |
| 2AB | -0.039 | 0.9646 | 0.0% | HIGH_ID |
| adenosine 5'-diphosphate (ADP) | 0.079 | 0.5419 | 0.0% | HIGH_ID |
| fumarate | -0.181 | 0.338 | 0.0% | HIGH_ID |
| 5-aminovalerate | 0.119 | 0.7198 | 2.6% | HIGH_ID |
| sarcosine | -0.124 | 0.6041 | 10.8% | HIGH_ID |
| nicotinamide adenine dinucleotide phosphate (NADP+) | -0.110 | 0.6818 | 15.4% | HIGH_ID |
| NADPH | -0.065 | 0.7163 | 12.8% | HIGH_ID |
| betaine | -0.296 | 0.2963 | 0.0% | HIGH_ID |
| nicotinamide adenine dinucleotide (NAD+) | 0.218 | 0.383 | 0.0% | HIGH_ID |
| pipecolate | -0.224 | 0.2328 | 18.9% | HIGH_ID |
| citrulline | -0.208 | 0.2755 | 0.0% | HIGH_ID |
| ornithine | -0.219 | 0.4278 | 0.0% | HIGH_ID |
| 4-Acetylbutyrate | -0.063 | 0.7065 | 45.9% | HIGH_ID |
| trimethylamine N-oxide | -0.046 | 0.803 | 25.6% | HIGH_ID |
| UDP-glucose | -0.480 | 0.119 | 2.7% | HIGH_ID |
| 3-methylhistidine | -0.284 | 0.3557 | 0.0% | HIGH_ID |
| 2-aminoadipate | -0.234 | 0.5419 | 2.7% | HIGH_ID |
| cytidine diphosphate | 0.133 | 0.4647 | 28.2% | HIGH_ID |
| Isethionate | 0.006 | 0.9427 | 25.6% | MODERATE_NAME_ONLY |
| guanosine 5'-triphosphate (GTP) | 0.149 | 0.3954 | 7.7% | HIGH_ID |
| 2-hydroxybutyrate (AHB) | -0.336 | 0.1219 | 17.9% | HIGH_ID |
| 2-Hydroxypentanoate | -0.344 | 0.2963 | 27.0% | MODERATE_NAME_ONLY |
| Cyclohexylamine | 0.186 | 0.4669 | 79.5% | HIGH_ID |
| glucose 1-phosphate | -0.172 | 0.5983 | 28.2% | HIGH_ID |
| glucuronate | -0.445 | 0.1027 | 38.5% | HIGH_ID |
| sebacate (decanedioate) | -0.290 | 0.1219 | 18.9% | HIGH_ID |
| uridine 5'-monophosphate (UMP) | -0.495 | 0.08596 | 0.0% | HIGH_ID |
| Propionate | -0.258 | 0.336 | 21.6% | HIGH_ID |
| 3-hydroxybutyrate (BHBA) | -0.337 | 0.2013 | 13.5% | HIGH_ID |
| cytidine 5'-monophosphate (5'-CMP) | -0.569 | 0.05398 | 0.0% | HIGH_ID |
| Hexanoate | -0.095 | 0.7394 | 5.4% | HIGH_ID |
| malonate | -0.300 | 0.1393 | 27.0% | HIGH_ID |
| glucose 6-phosphate | -0.358 | 0.1618 | 0.0% | HIGH_ID |
| glutarate (C5-DC) | -0.378 | 0.1219 | 7.7% | HIGH_ID |
| uridine 5'-triphosphate (UTP) | 0.179 | 0.5419 | 33.3% | HIGH_ID |
| adenosine 5'-triphosphate (ATP) | -0.205 | 0.5209 | 0.0% | HIGH_ID |
| piperidine | -0.512 | 0.153 | 40.5% | HIGH_ID |
| pelargonate (9:0) | -0.357 | 0.07963 | 35.1% | HIGH_ID |
| nicotinamide adenine dinucleotide reduced (NADH) | -0.346 | 0.2513 | 24.3% | HIGH_ID |
| Octanoate | -0.189 | 0.6292 | 38.5% | HIGH_ID |
| gluconate | -0.429 | 0.169 | 10.3% | HIGH_ID |

# 下一步与完成边界

效应记录描述已完成。下一步从这批显著特征逐一核对化合物身份和直接生化关系，重新建立有出处的 COAD 映射。现有基因映射未使用，RNA/样本分析未运行，功能与外部验证未运行。

来源：`reference/camp/cancer_effects.tsv`，SHA-256：`1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597`。可复现脚本：`code/coad/review_effect_records.py`。完整参数与状态见同目录 `summary.json`。
