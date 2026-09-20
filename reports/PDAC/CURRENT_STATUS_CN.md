# PDAC 当前状态

2026-09-20 用户在当前任务明确指定本账号接手 PDAC；使用独立工作树与 analysis/pdac-initial 分支。

01_CAMP 冻结汇总整理已完成：303 条效应、42 条原 q<0.05，其中原 g 正 12 条、负 30 条。全部原列逐项保留，未重新计算 P/q；原始 P 与样本量不在快照中，未反推。42 条显著特征均不是 X- 未知名称，但不代表化学身份已经重新确认。

结果入口：`results/PDAC/01_CAMP/20260920T051000Z_frozen_v1/README_CN.md`；全阶段索引：`coordination/stages/PDAC.tsv`。

03_PATIENT 为 ACCESS_BLOCKED：本机执行 `ssh -o BatchMode=yes -o ConnectTimeout=12 server165 pwd`，返回无法解析主机名，尚未到认证步骤。已向用户请求可用连接方式；不索取密码或私钥。源矩阵和患者资料继续留在服务器。

02/04/05/06/07 尚未开展。接通服务器后先复用既有 PDAC 统计和映射、确认作者样本身份及 RNA 覆盖，再固定候选检验族。不能将此首批整理视为完成 PDAC 研究或独立验证。

分支发布与 main 合并分开记录；实际提交与远程状态见 Git 记录。
