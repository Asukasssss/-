# GSE274103：目标基因覆盖筛查

## 本轮问题与输入
寻找可判断 PDAC LYPLA1 癌上皮富集的空间队列。[GSE274103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274103) 报告五位未治疗患者。本次实际下载检查首张 GSM8443449 / PDAC-p1 的作者过滤 H5。

## 结果与解释
矩阵为 17,943 特征 × 4,987 点，LYPLA1 符号和 ENSG00000120992 均无匹配。状态 NOT_EVALUABLE，未做表达差异分析。缺少目标行不等于零表达，不能用 LYPLA2 替代。

## 限制、决定与下一步
只确认这一张切片，未把其他四张推定为缺失。该样本不适合本问题；转用明确包含目标、可配准的 GSE278687。

## 复现
`python3 screen_spatial_h5.py --input <source.h5> --output target_coverage.json`。输入哈希与尺寸见 target_coverage.json；源 H5 仅存 server165。
