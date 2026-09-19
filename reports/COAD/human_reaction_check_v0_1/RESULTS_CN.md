# COAD 人类蛋白直接反应注释核对

对新检索的 693 条特征—基因候选逐条核对。482 条能与 reviewed 人类 UniProt 的催化反应通过同一 Rhea 反应对应；其中 424 条当前没有身份暂挂，涉及 58 个特征、322 个人类基因。

这里的“没有身份暂挂”只指本次已发现的特定问题，不代表原始鉴定、全部 HMDB 交叉核验或立体化学已经解决。反应对应是已有生化知识支持，不是 COAD 功能证据；当前没有锁定患者检验家族。

| 核对结果 | 关系数 |
|---|---:|
| MATCHED_REACTION_WITH_EXPERIMENTAL_ANNOTATION | 390 |
| NO_EXACT_REACTION_MATCH_IN_CURRENT_SEARCH | 198 |
| MATCHED_REACTION_OTHER_ANNOTATION | 92 |
| NO_REVIEWED_HUMAN_RECORD_LINKED | 13 |

ECO:0000269 是 UniProt 所记录的实验注释依据，本文未逐篇复核原文或实验条件。其他注释可能包含推断/相似性依据。未找到完全相同 Rhea 反应的关系保留待审，不当作反证。

# 按原代谢特征列出没有当前身份暂挂的匹配基因

| 原特征 | 匹配基因数 | 基因 |
|---|---:|---|
| 1-methylhistamine | 2 | AOC1, HNMT |
| 1-methylnicotinamide | 1 | NNMT |
| 2-phosphoglycerate | 9 | BPGM, ENO1, ENO2, ENO3, ENO4, MINPP1, PGAM1, PGAM2, PGAM4 |
| 3-hydroxyproline | 1 | L3HYPDH |
| 3PG | 8 | BPGM, PGAM1, PGAM2, PGAM4, PGK1, PGK2, PGM2L1, PHGDH |
| 4-methyl-2-oxopentanoate | 4 | BCAT1, BCAT2, BCKDHA, BCKDHB |
| 5-hydroxylysine | 1 | HYKK |
| 5-methylthioadenosine (MTA) | 5 | DPH1, LACC1, MTAP, SMS, SRM |
| 6-phosphogluconate | 4 | H6PD, IDNK, PGD, PGLS |
| N-Acetylglucosamine 1-phosphate | 2 | PGM3, UAP1 |
| N-acetylaspartate (NAA) | 4 | ASPA, NAT8L, RIMKLA, RIMKLB |
| N-acetylglucosamine 6-phosphate | 4 | AMDHD2, GNPNAT1, NAGK, PGM3 |
| N-acetylneuraminate | 2 | NANP, NPL |
| S-adenosylmethionine (SAM) | 33 | AMD1, ASMT, CARM1, CARNMT1, DPH1, GAMT, GNMT, HNMT, ICMT, KMT5A, LCMT2, LIAS, MAT1A, MAT2A, METTL3, MOCS1, NNMT, PEMT, PNMT, PRDM9, PRMT1, PRMT2, PRMT3, PRMT5, PRMT6, PRMT7, PRMT8, PRMT9, SETD2, SETD7, SMYD2, TPMT, TRMT5 |
| UDP-N-acetylglucosamine | 6 | B3GNT3, EOGT, GNE, OGT, POMGNT2, UAP1 |
| UDP-glucuronate | 21 | UGDH, UGT1A1, UGT1A10, UGT1A3, UGT1A4, UGT1A5, UGT1A6, UGT1A7, UGT1A8, UGT1A9, UGT2A1, UGT2A2, UGT2A3, UGT2B10, UGT2B11, UGT2B15, UGT2B17, UGT2B28, UGT2B4, UGT2B7, UXS1 |
| acetylcholine | 2 | ACHE, CHAT |
| adenine | 3 | APRT, LACC1, MTAP |
| adenosine | 8 | ADA, ADA2, ADK, AHCY, LACC1, NT5C1A, NT5C1B, NT5E |
| adenylosuccinate | 3 | ADSL, ADSS1, ADSS2 |
| allantoin | 1 | URAD |
| asparagine | 4 | ASNS, ASRGL1, NARS1, NARS2 |
| beta-alanine | 5 | CARNS1, CNDP1, CSAD, GADL1, UPB1 |
| carnosine | 3 | CARNMT1, CARNS1, CNDP1 |
| cis-aconitate | 1 | ACOD1 |
| citrate | 5 | ACLY, ACO1, ACO2, CS, RIMKLB |
| creatine | 6 | CKB, CKM, CKMT1A, CKMT1B, CKMT2, GAMT |
| cytidine | 7 | CDA, NT5C3A, NT5C3B, NT5E, UCK1, UCK2, UCKL1 |
| cytidine 5'-diphosphocholine | 4 | CEPT1, CHPT1, PCYT1A, PCYT1B |
| dihydroxyacetone phosphate (DHAP) | 8 | ALDOA, ALDOB, ALDOC, GNPAT, GPD1, GPD1L, TKFC, TPI1 |
| glucose | 17 | ADPGK, B4GALT1, B4GALT2, G6PC1, G6PC2, G6PC3, GANAB, GBA1, GBA2, GCK, HK1, HK2, HK3, HKDC1, MGAM, SI, TREH |
| glutamine | 15 | ASNS, CAD, CTPS1, CTPS2, GATB, GFPT1, GFPT2, GLS, GLS2, GLUL, GMPS, NADSYN1, PFAS, PPAT, QARS1 |
| glutathione, reduced (GSH) | 43 | CHAC1, CHAC2, ESD, ETHE1, GGT1, GGT5, GGT6, GGT7, GLO1, GPX1, GPX2, GPX3, GPX4, GPX5, GPX6, GSR, GSS, GSTA1, GSTA2, GSTA3, GSTA4, GSTA5, GSTK1, GSTM1, GSTM2, GSTM3, GSTM4, GSTM5, GSTO1, GSTO2, GSTP1, GSTT1, GSTT2, GSTT2B, GSTT4, HAGH, HPGDS, LANCL1, LTC4S, MGST1, MGST2, SQOR, TXNDC12 |
| glycerophosphorylcholine (GPC) | 8 | GDPD5, GPCPD1, PLA2G15, PLB1, PNPLA6, PNPLA7, TMEM86A, TMEM86B |
| guanidinoacetate | 2 | GAMT, GATM |
| guanine | 5 | GDA, HPRT1, LACC1, PNP, QTRT1 |
| guanosine | 4 | LACC1, NT5C2, NT5E, PNP |
| histamine | 3 | AOC1, HDC, HNMT |
| hypotaurine | 5 | ADO, CSAD, FMO1, FMO3, GADL1 |
| hypoxanthine | 4 | HPRT1, LACC1, PNP, XDH |
| inosine | 8 | ADA, ADA2, LACC1, NT5C, NT5C1A, NT5C2, NT5E, PNP |
| inosine 5'-monophosphate (IMP) | 22 | ADSS1, ADSS2, AMPD1, AMPD2, AMPD3, ATIC, CANT1, ENTPD1, ENTPD3, ENTPD5, ENTPD6, GMPR, GMPR2, HPRT1, IMPDH1, IMPDH2, ITPA, NT5C, NT5C1A, NT5C2, NT5E, NUDT16 |
| isocitrate | 5 | ACO1, ACO2, IDH1, IDH2, IDH3A |
| isoleucine | 4 | BCAT1, BCAT2, IARS1, IARS2 |
| kynurenine | 3 | AFMID, KMO, KYNU |
| leucine | 4 | BCAT1, BCAT2, LARS1, LARS2 |
| ophthalmate | 1 | GSS |
| phosphoenolpyruvate (PEP) | 9 | ENO1, ENO2, ENO3, ENO4, NANS, PCK1, PCK2, PKLR, PKM |
| proline | 6 | EPRS1, PARS2, PRODH, PYCR1, PYCR2, PYCR3 |
| serotonin (5HT) | 2 | DDC, MAOA |
| spermidine | 4 | PAOX, SMOX, SMS, SRM |
| taurine | 5 | BAAT, CSAD, FMO1, FMO3, GADL1 |
| tyrosine | 6 | IL4I1, PAH, TAT, TH, YARS1, YARS2 |
| uracil | 3 | DPYD, UPP1, UPP2 |
| urate | 1 | XDH |
| urea | 3 | AGMAT, ARG1, ARG2 |
| uridine | 9 | CDA, NT5C, NT5E, NT5M, UCK1, UCK2, UCKL1, UPP1, UPP2 |
| uridine 5'-diphosphate (UDP) | 56 | AK9, ALG5, B3GALNT2, B3GNT3, B4GALT1, B4GALT2, B4GALT4, B4GALT5, B4GALT6, CANT1, CHSY1, CHSY3, CMPK1, COLGALT1, COLGALT2, ENTPD1, ENTPD3, ENTPD4, ENTPD5, ENTPD6, ENTPD8, EOGT, GNE, GYG1, GYG2, GYS1, GYS2, NME1, NME2, NME3, NME4, NME6, OGT, PLOD3, POMGNT2, UGCG, UGT1A1, UGT1A10, UGT1A3, UGT1A4, UGT1A5, UGT1A6, UGT1A7, UGT1A8, UGT1A9, UGT2A1, UGT2A2, UGT2A3, UGT2B10, UGT2B11, UGT2B15, UGT2B17, UGT2B28, UGT2B4, UGT2B7, UGT8 |

# 方法与交付边界

从每条 KEGG 反应的 Rhea 交叉引用出发，利用 Rhea 官方方向表统一同一反应的方向编号；只与 NCBI GeneID 精确对应的 reviewed 人类 UniProt 催化反应相连。原反应中必须明确包含待查化合物，不能靠通路、相似基因名或宽泛 EC 编号连接。保留复合体注释，不能把所有亚基都解释为独立催化酶。
此版本仅覆盖第一轮 KEGG 反应候选内的关系；系统转运体检索、KO链未捕获的其他酶关系、全文功能审查、患者分析尚未完成。
逐条反应、UniProt 链接、证据代码和 PubMed 编号均保存在 human_reaction_pairs.json。可公开同步的是这些数据库来源与汇总，不含患者矩阵。

来源：[Rhea 官方下载说明](https://www.rhea-db.org/help/download)、[UniProt](https://www.uniprot.org/)。
