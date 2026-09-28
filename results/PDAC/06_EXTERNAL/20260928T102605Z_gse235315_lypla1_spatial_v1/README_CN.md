# GSE235315：LYPLA1 测量覆盖筛查

为补齐 PDAC 的真实空间切片，下载并检查 GSM7498812 / SS1923404 的处理包。该样本包含 17,943 个特征、2,248 个过滤后点位，但基因符号 LYPLA1 和 Ensembl ENSG00000120992 均无匹配。

因此本样本不能用于 LYPLA1 表达制图，状态为 NOT_EVALUABLE。缺失基因行不等于零表达，不使用 LYPLA2 或 LYPLAL1 替代。本结论只适用于实际检查的这一张切片，其他未完成下载的样本不作判断。未运行差异检验或富集分析。

来源：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7498812

源文件与 SHA-256 见 target_coverage.json；代码为 code/pdac/screen_visium_target_gene.py。原始数据仅存于 server165。下一步转用确实包含目标计数的 GSE111672，其独立结果目录为 20260928T104922Z_moncada_lypla1_slice_v1。
