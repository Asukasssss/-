# COAD 从效应记录独立建立的身份与反应候选初稿

只以冻结表中 73 条 COAD 原 q<0.05 记录为入口，重新读取 KEGG 化合物、反应、KO 和人类基因注释。未读取旧基因映射或 BRCA 候选。

本次 60 个特征找到候选反应链，13 个未找到；共 693 个唯一特征—人类基因候选组合、489 个基因，底层有 2199 条反应支持记录。

**这些数值是数据库候选检索的覆盖，不是经逐基因底物核实的正式映射，更不是显著相关或功能验证。当前进入患者统计的正式锁定关系数为 0。**

采用化合物在反应方程中明确出现 → 反应 KO 注释 → KEGG 人类基因的链路。没有把同一通路的全部基因连进来，也没有仅凭一个宽泛 EC 编号批量展开。方程左右侧仅按书写记录，不推断组织内反应方向。多来源/多反应归并到同一个特征—基因组合；合并峰不拆成独立测量。

# 身份问题与边界

- glutamate：KEGG C00025（L 型）与 HMDB03339（D 型）冲突，两种身份假设保留并暂挂。
- ribulose/xylulose 5-phosphate：合并名称，C00199/C00231 只作为两个假设，不能把一个测量拆成两个。
- cysteine-glutathione disulfide：原 KEGG 字段为 R00900 反应编号；从方程提出 C05526，尚未替换原注释。
- Mucate：本次 KEGG 批量响应未返回原 C01807；名称检索提出 C00879，尚待作者注释核对。
- ophthalmate：本次名称检索新提出 C21016，保留原 HMDB 键，交叉核验未完成。
- Disulfiram、Acetohydroxamate、triethanolamine、Phthalate、2-Deoxyglucose 6-phosphate 暂列原注释/暴露背景核查，不由药物靶点倒推直接代谢关系。
- HIGH_ID 不能消除上述问题。其余名称与 KEGG 注释相容也不等于重新鉴定；HMDB 尚未全面核查，立体化学需结合作者方法。

# 全部 73 条的去向

| 原特征 | 身份状态 | 候选基因数 | 候选基因（仅展示前 12 个） |
|---|---|---:|---|
| S-adenosylmethionine (SAM) | ORIGINAL_KEGG_NAME_COMPATIBLE | 85 | AMD1, ASH1L, ASMT, BHMT2, CAMKMT, CARM1, CARNMT1, CDKAL1, CMTR1, CMTR2, CNPY3-GNMT, COMT … |
| glutathione, reduced (GSH) | ORIGINAL_KEGG_NAME_COMPATIBLE | 48 | CHAC1, CHAC2, ESD, ETHE1, GGT1, GGT5, GGT6, GGT7, GLO1, GPX1, GPX2, GPX3 … |
| beta-alanine | ORIGINAL_KEGG_NAME_COMPATIBLE | 17 | ABAT, ALDH1B1, ALDH2, ALDH3A1, ALDH3A2, ALDH3B1, ALDH3B2, ALDH7A1, ALDH9A1, CARNS1, CNDP1, CNDP2 … |
| 2-Deoxyglucose 6-phosphate | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW | 0 |  |
| hypotaurine | ORIGINAL_KEGG_NAME_COMPATIBLE | 10 | ADO, CSAD, FMO1, FMO2, FMO3, FMO4, FMO5, GAD1, GAD2, GADL1 |
| N-acetylglucosamine 6-phosphate | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | AMDHD2, GNPNAT1, NAGK, PGM3 |
| uracil | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | DPYD, TYMP, UPP1, UPP2 |
| adenylosuccinate | ORIGINAL_KEGG_NAME_COMPATIBLE | 3 | ADSL, ADSS1, ADSS2 |
| cysteine-glutathione disulfide | HOLD_REACTION_ID_IN_COMPOUND_FIELD | 0 |  |
| cytidine 5'-diphosphocholine | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | CEPT1, CHPT1, PCYT1A, PCYT1B |
| N-acetylmethionine | ORIGINAL_KEGG_NAME_COMPATIBLE | 0 |  |
| hypoxanthine | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | HPRT1, LACC1, PNP, XDH |
| UDP-N-acetylglucosamine | ORIGINAL_KEGG_NAME_COMPATIBLE | 10 | B3GNT2, B3GNT3, B3GNT4, DPAGT1, EOGT, GNE, OGT, POMGNT2, UAP1, UAP1L1 |
| 1-methylnicotinamide | ORIGINAL_KEGG_NAME_COMPATIBLE | 2 | AOX1, NNMT |
| proline | ORIGINAL_KEGG_NAME_COMPATIBLE | 10 | EPRS1, LAP3, P4HA1, P4HA2, P4HA3, PARS2, PRODH, PYCR1, PYCR2, PYCR3 |
| kynurenine | ORIGINAL_KEGG_NAME_COMPATIBLE | 6 | AADAT, AFMID, KMO, KYAT1, KYAT3, KYNU |
| 5-methylthioadenosine (MTA) | ORIGINAL_KEGG_NAME_COMPATIBLE | 5 | DPH1, LACC1, MTAP, SMS, SRM |
| 3-hydroxy-3-methylglutarate | ORIGINAL_KEGG_NAME_COMPATIBLE | 0 |  |
| glycylglycine | ORIGINAL_KEGG_NAME_COMPATIBLE | 0 |  |
| taurine | ORIGINAL_KEGG_NAME_COMPATIBLE | 15 | BAAT, CSAD, FMO1, FMO2, FMO3, FMO4, FMO5, GAD1, GAD2, GADL1, GGT1, GGT5 … |
| 1-methyladenosine | ORIGINAL_KEGG_NAME_COMPATIBLE | 0 |  |
| N-Acetylglucosamine 1-phosphate | ORIGINAL_KEGG_NAME_COMPATIBLE | 3 | PGM3, UAP1, UAP1L1 |
| UDP-glucuronate | ORIGINAL_KEGG_NAME_COMPATIBLE | 25 | CHPF, CHPF2, CHSY1, CHSY3, UGDH, UGT1A1, UGT1A10, UGT1A3, UGT1A4, UGT1A5, UGT1A6, UGT1A7 … |
| 5-hydroxylysine | ORIGINAL_KEGG_NAME_COMPATIBLE | 1 | HYKK |
| adenine | ORIGINAL_KEGG_NAME_COMPATIBLE | 5 | APRT, LACC1, MTAP, PNP, XDH |
| isoleucine | ORIGINAL_KEGG_NAME_COMPATIBLE | 5 | BCAT1, BCAT2, IARS1, IARS2, IL4I1 |
| leucine | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | BCAT1, BCAT2, LARS1, LARS2 |
| tyrosine | ORIGINAL_KEGG_NAME_COMPATIBLE | 12 | DDC, GOT1, GOT1L1, GOT2, IL4I1, PAH, TAT, TH, TPO, TYR, YARS1, YARS2 |
| asparagine | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | ASNS, ASRGL1, NARS1, NARS2 |
| 3-hydroxyproline | ORIGINAL_KEGG_NAME_COMPATIBLE | 1 | L3HYPDH |
| spermidine | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | PAOX, SMOX, SMS, SRM |
| glutamate | HOLD_ID_CONFLICT | 58 | AADAT, AASS, ABAT, AGXT, AGXT2, ALDH18A1, ALDH4A1, ASNS, BCAT1, BCAT2, CAD, CISD1 … |
| uridine 5'-diphosphate (UDP) | ORIGINAL_KEGG_NAME_COMPATIBLE | 72 | AK3, AK9, ALG5, B3GALNT2, B3GNT2, B3GNT3, B3GNT4, B4GALNT1, B4GALT1, B4GALT2, B4GALT4, B4GALT5 … |
| cytidine | ORIGINAL_KEGG_NAME_COMPATIBLE | 15 | CDA, CDADC1, NT5C, NT5C1A, NT5C1B, NT5C1B-RDH14, NT5C2, NT5C3A, NT5C3B, NT5DC4, NT5E, NT5M … |
| uridine | ORIGINAL_KEGG_NAME_COMPATIBLE | 15 | CDA, CDADC1, NT5C, NT5C1A, NT5C1B, NT5C1B-RDH14, NT5C2, NT5DC4, NT5E, NT5M, UCK1, UCK2 … |
| 1-methylhistamine | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | AOC1, HNMT, MAOA, MAOB |
| guanidinoacetate | ORIGINAL_KEGG_NAME_COMPATIBLE | 2 | GAMT, GATM |
| inosine | ORIGINAL_KEGG_NAME_COMPATIBLE | 12 | ADA, ADA2, LACC1, NT5C, NT5C1A, NT5C1B, NT5C1B-RDH14, NT5C2, NT5DC4, NT5E, NT5M, PNP |
| creatine | ORIGINAL_KEGG_NAME_COMPATIBLE | 6 | CKB, CKM, CKMT1A, CKMT1B, CKMT2, GAMT |
| triethanolamine | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW | 0 |  |
| Phthalate | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW | 0 |  |
| N-acetylneuraminate | ORIGINAL_KEGG_NAME_COMPATIBLE | 7 | CMAS, NANP, NEU1, NEU2, NEU3, NEU4, NPL |
| ribulose/xylulose 5-phosphate | HOLD_COMBINED_FEATURE | 9 | FGGY, PGD, RPE, RPEL1, RPIA, TKT, TKTL1, TKTL2, XYLB |
| azelate (nonanedioate) | ORIGINAL_KEGG_NAME_COMPATIBLE | 0 |  |
| urea | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | AGMAT, ALLC, ARG1, ARG2 |
| creatinine | ORIGINAL_KEGG_NAME_COMPATIBLE | 0 |  |
| glutamine | ORIGINAL_KEGG_NAME_COMPATIBLE | 17 | ASNS, CAD, CTPS1, CTPS2, GATB, GATC, GFPT1, GFPT2, GLS, GLS2, GLUL, GMPS … |
| Acetohydroxamate | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW | 0 |  |
| guanosine | ORIGINAL_KEGG_NAME_COMPATIBLE | 10 | LACC1, NT5C, NT5C1A, NT5C1B, NT5C1B-RDH14, NT5C2, NT5DC4, NT5E, NT5M, PNP |
| phosphoenolpyruvate (PEP) | ORIGINAL_KEGG_NAME_COMPATIBLE | 9 | ENO1, ENO2, ENO3, ENO4, NANS, PCK1, PCK2, PKLR, PKM |
| cis-aconitate | ORIGINAL_KEGG_NAME_COMPATIBLE | 3 | ACO1, ACO2, ACOD1 |
| glycerophosphorylcholine (GPC) | ORIGINAL_KEGG_NAME_COMPATIBLE | 10 | GDPD5, GPCPD1, LYPLA1, LYPLA2, PLA2G15, PLB1, PNPLA6, PNPLA7, TMEM86A, TMEM86B |
| 6-phosphogluconate | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | H6PD, IDNK, PGD, PGLS |
| carnosine | ORIGINAL_KEGG_NAME_COMPATIBLE | 4 | CARNMT1, CARNS1, CNDP1, CNDP2 |
| allantoin | ORIGINAL_KEGG_NAME_COMPATIBLE | 1 | URAD |
| inosine 5'-monophosphate (IMP) | ORIGINAL_KEGG_NAME_COMPATIBLE | 28 | ADSS1, ADSS2, AMPD1, AMPD2, AMPD3, ATIC, CANT1, ENTPD1, ENTPD3, ENTPD4, ENTPD5, ENTPD6 … |
| 3PG | ORIGINAL_KEGG_NAME_COMPATIBLE | 10 | ACYP1, ACYP2, BPGM, PGAM1, PGAM2, PGAM4, PGK1, PGK2, PGM2L1, PHGDH |
| citrate | ORIGINAL_KEGG_NAME_COMPATIBLE | 5 | ACLY, ACO1, ACO2, CS, RIMKLB |
| dihydroxyacetone phosphate (DHAP) | ORIGINAL_KEGG_NAME_COMPATIBLE | 9 | ALDOA, ALDOB, ALDOC, GNPAT, GPD1, GPD1L, GPD2, TKFC, TPI1 |
| N-acetylaspartate (NAA) | ORIGINAL_KEGG_NAME_COMPATIBLE | 5 | ASPA, FOLH1, NAT8L, RIMKLA, RIMKLB |
| guanine | ORIGINAL_KEGG_NAME_COMPATIBLE | 5 | GDA, HPRT1, LACC1, PNP, QTRT1 |
| ophthalmate | PROVISIONAL_NEW_KEGG_ID | 1 | GSS |
| urate | ORIGINAL_KEGG_NAME_COMPATIBLE | 1 | XDH |
| isocitrate | ORIGINAL_KEGG_NAME_COMPATIBLE | 7 | ACO1, ACO2, IDH1, IDH2, IDH3A, IDH3B, IDH3G |
| 4-methyl-2-oxopentanoate | ORIGINAL_KEGG_NAME_COMPATIBLE | 6 | BCAT1, BCAT2, BCKDHA, BCKDHB, DBT, DLD |
| histamine | ORIGINAL_KEGG_NAME_COMPATIBLE | 3 | AOC1, HDC, HNMT |
| adenosine | ORIGINAL_KEGG_NAME_COMPATIBLE | 16 | ADA, ADA2, ADK, AHCY, AHCYL1, AHCYL2, LACC1, NT5C, NT5C1A, NT5C1B, NT5C1B-RDH14, NT5C2 … |
| Mucate | HOLD_ORIGINAL_ENTRY_NOT_RETURNED | 0 |  |
| Disulfiram | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW | 0 |  |
| acetylcholine | ORIGINAL_KEGG_NAME_COMPATIBLE | 2 | ACHE, CHAT |
| 2-phosphoglycerate | ORIGINAL_KEGG_NAME_COMPATIBLE | 10 | BPGM, ENO1, ENO2, ENO3, ENO4, GLYCTK, MINPP1, PGAM1, PGAM2, PGAM4 |
| glucose | ORIGINAL_KEGG_NAME_COMPATIBLE | 24 | ADPGK, B4GALT1, B4GALT2, G6PC1, G6PC2, G6PC3, GAA, GANAB, GANC, GBA1, GBA2, GBA3 … |
| serotonin (5HT) | ORIGINAL_KEGG_NAME_COMPATIBLE | 8 | AANAT, DDC, IDO1, IDO2, INMT, MAOA, MAOB, METTL6 |

# 完成项、未完成项与下一步

DONE：73 条的数据库名称/标识初查、KEGG 反应链检索、按实际特征与人类基因去重、出处和访问时间/哈希记录。
NOT_RUN：系统转运体检索、逐人类蛋白底物与复合体角色核对、RNA/患者关联、依赖/功能分析、外部复现。无候选链不是阴性，不能淘汰代谢物。
ACCESS_BLOCKED：本机未配置可解析的 server165 SSH 别名；部分 HMDB 页面/XML 返回 HTTP403。未下载任何患者级矩阵。
下一步：回查有冲突的作者注释；从反应候选中核对人类蛋白的直接催化或复合体角色，并补转运体。身份与底物证据足够后才锁定患者分析关系和 BH 家族。

来源：[KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html)、[L-glutamate](https://www.kegg.jp/entry/C00025)、[D-glutamate HMDB](https://hmdb.ca/metabolites/HMDB0003339)、[合并名称的 ribulose 对应条目](https://www.kegg.jp/entry/C00199)。逐条链接见 identity_review.json 和 reaction_evidence.json；完整候选组合见 provisional_feature_gene_pairs.json。

机器校验：源哈希不变；73 条全部有去向；每条链的目标化合物确实在方程中、KO 同时存在于反应记录和人类映射；组合无重复。此校验不代替生物学底物审查。
