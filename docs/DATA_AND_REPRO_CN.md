# 数据位置与复用

两台电脑均由用户配置访问 server165；仓库不分发SSH密钥或GitHub令牌。账号B本机路径可不同，使用自己的clone目录。

服务器原项目根：`/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716`

相对该根：
- `data/candidates/camp_primary_tissue_multicancer/metadata/MasterMapping_MetImmune_03_16_2022_release.csv`
- `data/candidates/camp_primary_tissue_multicancer/processed_metabolomics/`
- `data/candidates/camp_primary_tissue_multicancer/gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed/`
- `results/CAMP_Phase1_v1.0/`：冻结历史
- `results/BRCA_117_screen_20260911/`：BRCA当前核心结果
- `results/CAMP_first_analysis_20260910/`：FUSCC匹配与旧小面板
- `results/CPTAC_BRCA_ten_20260911/`：蛋白初查

BRCA已核实来源：`PreprocessedData_BRCA1.xlsx`、`GSE37751.hugene10st.gene_symbol.csv`。其他队列先读取作者映射和配置，不按文件名猜患者或RNA尺度。

## 历史代码

`code/reference/camp_per_cancer_repro.py` 为可参考的跨癌统计模块；`run_brca117_patient_20260911.py` 是BRCA特定脚本，有硬编码路径和输出保护，不直接当作其他癌种命令。核对入口和配置后，在自己癌种的新脚本中适配；不要对所有癌种自动重跑。

任何重用代码先确认依赖与版本。BRCA历史为Python3.8.10/numpy1.24.4/pandas1.5.3/scipy1.10.1，服务器还有statsmodels。不是要求每台电脑重建患者数据环境；计算留服务器。

新服务器运行路径：`<原项目根>/results/collaborative/<CANCER>/<A或B>/<UTC日期时间_任务名>/`。每次新建不可覆盖目录并以独占创建方式建`.running`锁；运行正常结束再移除自己锁并写DONE。脚本、日志、参数、版本、哈希和结果均留此目录。账号A/B不要同时改原共享scripts。

只同步到GitHub：代码、参数、摘要、允许共享的小型汇总表与来源清单。不同步：原矩阵、逐患者数值、完整sample mapping、原始质谱或大型下载。上传前检查，不把.gitignore当作唯一防线。

reference内是初始化白名单副本；`reference/MANIFEST_SHA256.tsv` 可核对复制内容。CAMP冻结cancer_effects预期SHA256为1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597。
