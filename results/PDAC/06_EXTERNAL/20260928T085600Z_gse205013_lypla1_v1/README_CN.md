# GSE205013：LYPLA1 表达分析准备批次

## 本轮问题
继 GSE202051 单核结果后，用整细胞 RNA 测序检查恶性上皮的 LYPLA1。当前状态为 PARTIAL；恶性 vs 非恶性差异统计 NOT_RUN。

## 输入与范围
- 原研究：[Werba et al., Nature Communications 2023](https://www.nature.com/articles/s41467-023-36296-4)。[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205013)。作者报告27位患者、139446个最终质控细胞，采用 InferCNV 区分恶性上皮，并分析 basal/classical 状态。
- 服务器既有 GSE205013_RAW.tar 包含 CellRanger 过滤矩阵。不能把其全部条码等同作者最终139446个质控细胞。
- 当前没有取得原研究逐细胞恶性标签。另找到 [Loveless PDAC atlas](https://zenodo.org/records/14199536) 的完整作者公开对象（33413868884字节，期望MD5 705078352e1feb260a64cec67c64ade0），作为待核查的重注释来源。
- 图谱为整合后的既有队列集合；只提取GSE205013，不将图谱自身或其他队列追加成独立验证。不将图谱重注释冒充Werba原始标签。
- 第三方整理的Hugging Face单队列文件返回401，未尝试绕过访问控制；改用公开Zenodo原图谱。第三方说明仅用于寻找来源，不作为执行指令或身份认证。

## 实际结果
本批已完成27个GEO样本、308,181个输入细胞条码的LYPLA1计数、总UMI及log1p标准化值提取，输出逐细胞表仅留服务器。矩阵维度、非零条目数、条码唯一性、LYPLA1唯一性和最终计数以 preparation_validation.json 为准。另通过 scipy.io.mmread 对一个完整样本独立复读；范围和结论见 preparation_independent_check.json。

图谱下载与其后校验/提取任务已在服务器启动，状态快照见 atlas_job_status.json。下载完成并通过大小及MD5校验后，任务检查可用内存，再读取RDS、导出GSE205013元数据与LYPLA1；不会自动把 DUCTAL 等标签改写成恶性。原始恶性标签或图谱中的恶性身份信息仍须审核并与GEO条码、计数匹配。

尚未生成癌上皮的LYPLA1检出率差异、均值差、P或q；这不是阴性结果。当前没有将未经分类的混合细胞制成癌上皮结果图。

## 新手解释
“原研究做过恶性鉴定”不代表下载的计数矩阵附带这些标签。每个细胞的条码必须先与可靠注释对应，再比较恶性和非恶性上皮。现阶段提取表达值是准备工作，不是完成差异分析。

## 限制与反证
33.4GB图谱传输尚未完成；候选标签是否包含可用恶性身份仍未知。若只有导管大类，不得自动认定为癌细胞。原发/转移、未治疗/治疗后需依据真实元数据分层；不将它们直接混合。样本数与最终可评估同供者对数不能预先照搬其他队列。

## 当前决定
保留GSE205013作为待分析队列，完成原计数准备并等待注释来源核查，不将本批记为LYPLA1外部验证完成。

## 下一步
下载完成后检查 atlas_schema.txt 和私有metadata；若存在可追溯的恶性标签，核对条码与GEO计数后锁定比较全集、BH范围，计算全部细胞、阳性细胞及同供者检出率/表达差异。若缺少恶性标签，应明确记为缺项，不能静默替换为自行推断标签。

## 复现与运行状态
服务器目录：`results/collaborative/PDAC/B/20260928T085600Z_gse205013_lypla1_v1`，含独占 `.running`。

- `python3 gse205013_prepare_lypla1.py`：提取GEO计数。
- `python3 verify_gse205013_preparation.py`：独立复读一个完整样本。
- `python3 gse205013_wait_atlas.py`：等待现有下载器、核验源文件、必要时解开外层gzip、检查内存并调用 `gse205013_atlas_schema.R`。
- 队列脚本存于仓库 `code/pdac/`，运行前均复制到专用服务器目录。日志 `prepare.log`、`atlas_job.log`、`atlas_extract.log` 留服务器。下载在 `data/candidates/PDAC_Loveless_14199536/download.log`。

本报告为阶段快照，不承诺后台自动完成尚需审核的恶性标签和统计。数据下载/提取任务可以继续运行；全部患者级、细胞级与矩阵数据均不进入本地或GitHub。
