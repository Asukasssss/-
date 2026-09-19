# CAMP 癌种协作项目

本仓库让两台电脑、两个 Codex 账号按癌种协作。用户无湿实验条件，目标是用 CAMP 代谢差异、直接生化关系、患者关联和适用功能资料确定值得深入的候选，不宣称已经证明治疗靶点。

## 新账号阅读顺序

1. `AGENTS.md`：分析与协作边界。
2. `docs/PROJECT_CONTEXT_CN.md`：无需历史对话的项目背景与现状。
3. `docs/DATA_AND_REPRO_CN.md`：server165 位置、已有代码和复现限制。
4. `coordination/assignments.tsv`：癌种归属；未分配不自行认领。
5. `START_ACCOUNT_B_CN.md`：另一个账号的启动指令。

当前账号 A 负责 BRCA，账号 B 负责 COAD。用户指定公共仓库 https://github.com/Asukasssss/- 。仅发布脚本、项目说明与汇总结果；患者矩阵和凭据不上传。实际发布状态以远程提交为准。

各癌种每阶段按 [统一结果格式](docs/STAGE_STANDARD_CN.md) 交付并上传；标准入口为 `coordination/stages/<CANCER>.tsv`。已有历史快照保持原样；COAD 早期本地草案与统一阶段的对应关系见 [兼容说明](docs/STAGE_RESULT_STANDARD_CN.md)。

代码、说明、小型汇总在 GitHub；源矩阵与患者级数据只留 server165。仓库里的 `reference/` 是白名单选取的历史汇总/映射快照，`code/reference/` 是历史脚本，不是自动适配所有癌种的一键流水线。

## 每次工作

先更新 main，读取分工和对应癌种状态；建立 `analysis/<cancer>-<task>` 分支；只写自己癌种的新运行目录。提交前执行 `python tools/check_repository.py`。完成一批后更新该癌种状态和交接说明，通过 Pull Request 合并。不强制推送、不自动合并共享规则改动。

模板见 `docs/HANDOFF_TEMPLATE.md`；GitHub 发布步骤见 `docs/GITHUB_SETUP_CN.md`；实际发布状态见 `docs/DEPLOYMENT_STATUS_CN.md`（当前仅本地完成）。
