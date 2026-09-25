# PRAD 牛磺酸—SLC6A6 配对变化

本轮问题：同一患者 SLC6A6 的肿瘤-正常变化越低，牛磺酸变化是否也越低？

输入与范围：43对作者case配对；沿用作者处理尺度，不重新填补或变换；患者数据留服务器。单一关系为用户指定、同队列发现后的探索性补充。

实际结果：
- DELTA_PRIMARY: n=43, rho=-0.040320, P=0.802600, bootstrap95%CI=[-0.339152,0.253037]
- DELTA_AVAILABLE: n=43, rho=-0.040320, P=0.802600, bootstrap95%CI=[-0.339152,0.253037]
- DELTA_CAPT: n=36, rho=-0.021879, P=0.899400, bootstrap95%CI=[-0.344614,0.315877]
- DELTA_BATCH: n=43, rho=0.050701, P=0.760100, bootstrap95%CI=[-0.261104,0.357111]

方向组合：
- RNA down, taurine down: 26/43
- RNA down, taurine equal: 0/43
- RNA down, taurine up: 6/43
- RNA equal, taurine down: 0/43
- RNA equal, taurine equal: 0/43
- RNA equal, taurine up: 0/43
- RNA up, taurine down: 10/43
- RNA up, taurine equal: 0/43
- RNA up, taurine up: 1/43

新手解释：每点为一个配对病例；横纵轴均为肿瘤减正常。正相关意味着两种变化量倾向同向排列。共同降低人数多，不自动证明下降幅度相关，尤其两者各自已有较高降低比例。零值单列，不剔除以提高一致性。

限制/反证：P<0.05为本轮展示标准，单关系各检验族q=P，不消除历史筛选偏倚。配对相关不证明因果、酶活或通量，也不能排除细胞组成、批次等影响。CAPT字段有7对不一致，36对一致样本单列敏感性。病例身份沿用作者case，非基因型验证。多个敏感性分析不是独立复现。原作者处理尺度下差值不是浓度倍数。

当前决定：全部结果保留，不按显著性选择图或样本。

下一步：独立样本与细胞组成控制尚未完成。

复现命令：python prad_slc6a6_paired_delta_v1.py --out NEW_RUN --source-run ORIGINAL_RUN --code-commit PARENT_COMMIT；先建立NEW_RUN并复制本脚本及prad_patient_v1.py。
