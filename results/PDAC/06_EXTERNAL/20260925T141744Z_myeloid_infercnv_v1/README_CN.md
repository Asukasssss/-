# PDAC 髓系 inferCNV 与 SLC6A6：输入就绪记录

## 本轮问题
按用户请求，对髓系运行 inferCNV，再比较 SLC6A6。CNV异常不直接命名为恶性。
## 输入与范围
三套既有队列，各作者样本单位分别分析。排除chr3参与CNV推断；SLC6A6保留在后续表达比较中。所有逐细胞/逐样本值留在server165。
## 实际结果
输入已整理：GSE263733为11个单位/5078个髓系细胞；GSE278688为11个单位/5123个髓系细胞；GSE242230为9个单位/3932个髓系细胞。535个正常组织髓系对照、603个作者恶性上皮对照已接入。这个检查点尚未产生CNV分组及组间P/q。
## 新手解释
运行T+B、T-only、B-only三种参考。对inferCNV残差信号采用预先记录的探索性阈值；一致异常、参考样和不确定分别保留。自定义阈值不是经过验证的恶性分类器，也不是inferCNV HMM标签。
## 限制/反证
跨谱系表达、髓系亚群构成、双细胞及环境RNA均可能影响结果；本轮没有新双细胞检验或DNA/CNV克隆确认。排除chr3后无法评估chr3本身异常。正常组织髓系只覆盖部分单位。供者临床独立性未再次认证。
## 当前决定
PARTIAL：输入完成，依赖环境修复/安装及推断执行中；无结果不得报告恶性比例或SLC6A6差异。
## 下一步
先验证单样本inferCNV运行，再执行其余样本；数值完成后按原规则汇总。不用单细胞数替代独立供者数。
## 复现命令
在新的PDAC/B服务器目录复制code/pdac/myeloid_infercnv_*脚本及固定参数，依次执行preflight、prepare、schedule、summarize、manifest。schedule需要该运行目录env内可加载infercnv 1.26.0。详细参数见analysis_spec.json。

方法依据：https://github.com/broadinstitute/infercnv/wiki/Running-InferCNV
