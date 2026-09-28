# 发布状态

- 用户指定公共仓库：https://github.com/Asukasssss/- 。
- 账号 B 本机 Git Credential Manager 已完成 Asukasssss 登录，并通过远程分支推送预检查；origin 已配置。发布成功以目标远程分支提交核对为准；PR_OPEN 和 MERGED 分开，不能把分支上传写成 main 已合并。
- Codex GitHub 插件仍是独立连接，本次未修改插件账号。协作可直接使用本机 Git。
- 仓库只包含项目说明、脚本和汇总结果，不含患者矩阵或凭据。
- 账号 B 已通过用户提供的连接信息成功访问 server165，并在独占 COAD/B 新目录完成作者数据覆盖核对。连接信息与凭据不入仓库；源矩阵及患者级数据继续驻留服务器。服务器运行目录见本批 source_manifest.tsv。
- 账号 A 负责 BRCA；账号 B 负责 COAD，按用户要求从效应记录独立分析。
- COAD 阶段索引：`coordination/stages/COAD.tsv`；具体远程提交与 PR：`coordination/publications/COAD.json`。
- COAD 首轮患者关联、配对 RNA 和基础敏感性结果已推送至 `analysis/coad-initial`，PR #1 已核对，尚未合并 main；协变量、功能与外部阶段仍未完成。
- 本机旧 zip/bundle 是初始化快照，后续请以 GitHub 最新提交为准。
