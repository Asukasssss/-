# COAD 可用性敏感性 v2：相同输入复用

## 本轮问题

修正可用性分析中相同样本、相同方法重复随机置换的问题，避免将随机波动误作生物学敏感性。此修正由读取关联结果之前的样本可用性核对触发，无显著性驱动的选择。

## 输入与范围

输入为已完成 `20260919T122034Z_patient_v1` 的公开汇总表。完整974行、预设674关系、实际652检验不变。主分析已确认33名I–IV期患者；availability33沿用同一processed代谢物和RNA，仅按作者data表有限值排除样本。因此n=33时输入与主分析完全相同。预先提交的修正规则为 config/coad_availability_reuse_v2.json（06297f664e84a0b65e0c51c453f566d7c078452c）。

## 实际结果

437 条相同输入项精确复用主分析的rho、P、bootstrap区间、种子与重采样次数；215 条样本减少项精确复用首版统计。652 项一起重新计算该可用性家族 BH，仍为 0 条 q<0.05，显著/不显著分类无变化。原主分析、全部37敏感性、配对RNA及冻结CAMP统计均不变。

实际n分布：33人437条、32人95条、31人7条、30人20条、29人33条、28人1条、26人19条、25人7条、16人19条、11人14条。其余300条家族外待审、22条家族内RNA缺项继续保留。

## 新手解释

同一批人、同一组数值、同一种检验不应因重复抽随机数而被当作新的证据。这里直接复用已算结果；由于其余关系用了较少患者，整个家族的P分布仍可能不同，重新BH后，相同输入项的q也可能与主分析不同。这是检验家族分布差异，不是新增独立验证。

## 限制/反证

作者可用性不等于经证实的原始检测状态；样本减少可能引入选择偏差，小n检验力有限。不显著不能证明不存在关系。仍未进行协变量、纯度或技术批次调整。此批只修正统计复用，不解除化合物身份或功能证据限制。

## 当前决定

本版 results.tsv 替代首版 availability33.tsv 作为当前可用性结果，首版永久保留追溯。复用修正 DONE；04_ROBUSTNESS 整阶段仍 PARTIAL。完整候选保留，不按名义P另设事后小家族。

## 下一步

开展预先声明的全家族协变量敏感性与适用功能证据核对；不要将重复汇总计作第二份支持证据。

## 复现命令

在服务器 COAD/B 下新建专属目录和 `.running`，复制 `code/coad/reuse_availability_v2.py` 及 config/coad_availability_reuse_v2.json 为 locked_spec.json 后运行：

```sh
python3 <新目录>/reuse_availability_v2.py --base-dir <项目根>/results/collaborative/COAD/B/20260919T122034Z_patient_v1 --run-dir <新目录> --code-commit 06297f664e84a0b65e0c51c453f566d7c078452c
```

原目录不可覆盖。脚本只读取聚合结果，检查输入/输出哈希并记录每行 reuse_status；本轮无新患者数据传输。公开结果可本地验证：`python tools/validate_coad_availability_reuse.py --run-id 20260919T122940Z_availability_reuse_v2`；delivery_validation.json 包含逐项复用、独立BH复算、完整键和冻结值验证。
