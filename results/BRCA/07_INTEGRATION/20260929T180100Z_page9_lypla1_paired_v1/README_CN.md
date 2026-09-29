# 第9页：为什么聚焦LYPLA1

本轮问题：将六癌候选概览过渡到LYPLA1四癌专题，展示患者配对RNA变化，而不是重复排名或P值柱图。

输入与范围：CAMP BRCA45对、COAD33对、PDAC11对、PRAD43对。基因来自当前代谢物直接关系候选池。RNA矩阵和配对映射全程留server165，沿用既有身份审计；不按样本编号猜配对。BRCA排除已记录的组织标签冲突，COAD沿用原33患者资格，PDAC使用作者配对表，PRAD沿用43作者case配对。

实际结果：四癌平均差值均正；BRCA、COAD、PRAD名义P<0.05，PDAC P=0.07772254。患者上调/下调数见paired_summary.tsv。原P直接复用；重新提取仅用于核对原样本数、均差及绘图，没有新推断检验。

新手解释：每条细线连接同一患者的正常与肿瘤RNA；暖色升高、蓝色降低；短粗横线为组均值。下面标明多少患者上调。单名患者的一条线不代表其变化达到统计显著。

限制/反证：各面板原始表达尺度及纵轴范围分别保留，不能横向比较绝对表达量或把原差值当统一倍数。四癌同方向指平均变化，不是所有患者都上调。PDAC不显著。六癌概览中的ccRCC反向结果已在第8页展示，这页限定四癌，不据此宣称泛癌一致。CAMP内表达支持不是独立验证，更不是功能或代谢通量证据。

当前决定：关注理由为代谢关系候选加四癌同向表达背景；并不宣称LYPLA1是排名最优靶点。仅交付图片和矢量PDF，不生成PPT。

下一步：单细胞区分恶性上皮与正常上皮，判断组织RNA变化的细胞来源。

复现：把code/ppt_page9_lypla1_paired_server_v1.py、冻结04_lypla1_rna.tsv（命名expected.tsv）、sysu_logo.png及微软雅黑字体置于新的server165 BRCA/A运行目录；运行服务器脚本；使用code/ppt_page9_package_v1.py下载public成品。历史目录不覆盖。绘图输入路径与哈希见source_manifest_server.tsv；原统计Git来源见source_manifest_frozen.tsv。
