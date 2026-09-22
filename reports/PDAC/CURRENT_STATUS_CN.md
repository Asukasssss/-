# PDAC 当前进度：统一报告 v5

## 本轮问题
按45c6ae9交接规范接续已有CAMP至单细胞来源结果，补真实报告、完整表图、验收和ZIP。所有已有统计沿用，不从头重跑。
## 输入与范围
固定展示输入2971100b37e40d66a88cf89102c0fd3a56b7bcee，统计基线bcccc21ff60d6c57111efe72804d50db8020245e。27肿瘤/12正常标本；11作者配对；肿瘤关联21个不同作者病例键。原CAMP303/42冻结版本另存。
## 实际结果
- 配对代谢物307项：51项P<0.05、0项q<0.05；原有值敏感性193项可评估，27项P<0.05、0项q<0.05。
- 357直接关系中354可评估，65项P<0.05、1项q<0.05；358条件关系另族350可评估、4项q<0.05。
- 当前直接池RNA247/250可评估，91项P<0.05、40项q<0.05。437条件/历史补充基因另族429可评估，151项P<0.05、0项q<0.05。
- 687基因全部进入三队列来源，可解释来源数为621/593/610；三队列严格同首位120，且各保持率≥80%为72，当前直接池30。
- GSE242230：10783作者恶性和350作者正常上皮分别计算。另两队列导管身份未定；旧11,133细胞不可全称恶性。
- 本批交付16页主PDF、Markdown、16页签完整Excel、81页全基因热图册、全基因点图、独立图/图源/图注及ZIP。具体数值和验收见下方入口。
## 新手解释
显示“暂不可定位”不等于基因不存在。图册保留全部基因，示例不是排名。直接/条件/主/可用值/RNA的q按原检验族分别解释。MGLL与SLC6A19本轮直接主q分别0.059和0.0531，不挪用旧池q。AASS患者支持弱；ABAT—GABA有内部相关；TDO2有稳定CAF来源线索。
## 限制/反证
作者病例键不等于新认证临床身份；跨研究患者重叠待核。两队列恶性身份PARTIAL，临床分型及UMAP未接入并标NOT_RUN。来源不是机制或同关系外部验证；本轮不追加CNV、拟时序、通讯、生存或功能干预。
## 当前决定
本规范的来源与展示范围实际完成，研究缺项不改成DONE。原历史表不覆盖。PR #2仍为草稿、未合并main。
## 下一步
围绕完整证据表讨论候选；如要求恶性来源结论，先补可逐细胞连接的作者标签。
## 复现
python code/pdac/render_report_v5.py --config configs/PDAC_report_v5.yaml --validate-only
python code/pdac/render_report_v5.py --config configs/PDAC_report_v5.yaml --out <全新目录>

- [统一报告和全部入口](../../results/PDAC/07_INTEGRATION/20260922T140000Z_unified_report_v5/README_CN.md)
- [Markdown主报告](../../results/PDAC/07_INTEGRATION/20260922T140000Z_unified_report_v5/PDAC_主报告.md)
- [逐关系事实、解释与反证](../../results/PDAC/07_INTEGRATION/20260922T140000Z_unified_report_v5/interpretation_by_relation.tsv)
- [验收与缺项](../../results/PDAC/07_INTEGRATION/20260922T140000Z_unified_report_v5/acceptance.tsv)
- [v4作者身份结果](../../results/PDAC/06_EXTERNAL/20260922T131800Z_author_identity_v4/README_CN.md)
