"""Prepare aggregate-only public audit; source and spot measurements stay on server165."""
from pathlib import Path
import json,sys,shutil,hashlib
import pandas as pd
R=Path(sys.argv[1]);P=R/'public';P.mkdir(exist_ok=True)
root=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_DeepSpaceDB')
recs=json.loads((root/'expert_annotations.json').read_text())['records']
coverage=[]
for cohort in sorted(set(r['metadata']['series_id'] for r in recs)):
 rr=[r for r in recs if r['metadata']['series_id']==cohort]
 coverage.append({'cohort':cohort,'expert_annotated_sections':len(rr),'sections_with_normal_duct':sum('Normal duct' in r['labels'] for r in rr),'sections_with_invasive_and_normal_duct':sum('Normal duct' in r['labels'] and 'Invasive' in r['labels'] for r in rr)})
pd.DataFrame(coverage).to_csv(P/'annotation_coverage.tsv',sep='\t',index=False)
df=pd.read_csv(R/'section_effects_SERVER_ONLY.tsv',sep='\t')
records=[]
for cohort in ['GSE210616','GSE242311']:
 x=df[df.cohort==cohort]
 records.append(dict(cancer='BRCA',cohort=cohort,stage_id='06_EXTERNAL',run_id=R.name,analysis_version='deepspace_v1',analysis_type='invasive_vs_normal_duct_eligibility',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit='patient',n=0,n_reference=0,effect_type='NOT_EVALUABLE',effect='NA',ci_lower='NA',ci_upper='NA',p_value='NA',q_value='NA',test_family='patient-level invasive versus normal duct',family_n_evaluable=0,status='NOT_EVALUABLE',reason='One section with 6 normal spots; 2 interior spots; below prespecified 20' if cohort=='GSE210616' else 'LYPLA1 symbol and ENSG00000120992 absent in all three single-specimen matrices',source_id='DeepSpaceDB expert pathology + GEO',n_sections_inspected=x['sample'].nunique()))
pd.DataFrame(records).to_csv(P/'results.tsv',sep='\t',index=False)
for f in ['analysis_spec.json','validation.json']:shutil.copy2(R/f,P/f)
shutil.copy2(R/'input_hashes.tsv',P/'source_manifest.tsv')
audit=[
 ['GSE210616 + DeepSpaceDB','DONE','One previously analyzed section has expert invasive and normal duct annotations; direct descriptive comparison performed. Normal reference only 6 spots; not formal patient evidence.','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616'],
 ['GSE242311 + DeepSpaceDB','NOT_EVALUABLE','5 sections have relevant expert labels; 2 contain multiple specimens. Downloaded 3 single-specimen matrices; all 17943 features lack LYPLA1 and its Ensembl ID.','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242311'],
 ['GSE213688','NOT_EVALUABLE','15 author pathology files audited. Normal-containing sections generally lack definite invasive tumor; M11 has DCIS and 15 normal epithelial spots. DeepSpaceDB invasive label conflicts with original pathology; unresolved.','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213688'],
 ['SMART GSE331245','ACCESS_BLOCKED','Raw data public; required linked processed morphology/clinical data described as controlled via PharosAI. No usable public AOI group map verified.','https://www.biorxiv.org/content/10.1101/2025.08.27.672673v1.full'],
 ['CAF breast spatial Figshare 21591429','ACCESS_BLOCKED','Primary paper reports processed spatial data and normal ducts/lobules. Figshare API returned 403 in both local and server attempts; annotation-to-matrix contents unverified.','https://doi.org/10.6084/m9.figshare.21591429'],
 ['LYPLA1 GeoMx study 2023','ACCESS_BLOCKED','65 tumor-enriched and 6 normal epithelial ROIs reported. Article has no expression matrix accession; supplement is figures and clinical table. Author inquiry required. Prior LYPLA1 spatial/prognostic/IHC research exists.','https://pmc.ncbi.nlm.nih.gov/articles/PMC10463390/'],
 ['TLS-HEV Zenodo 18917100','ACCESS_BLOCKED','Record says public availability expected December 2026. Restricted files not accessed.','https://zenodo.org/records/18917100'],
 ['CAMEO-Breast','ACCESS_BLOCKED','Local Hugging Face connection timed out; LYPLA1 panel coverage and usable annotations not verified.','https://huggingface.co/datasets/theislab/CAMEO-Breast'],
]
pd.DataFrame(audit,columns=['dataset','status','reason','url']).to_csv(P/'candidate_audit.tsv',sep='\t',index=False)
readme='''# LYPLA1：乳腺癌空间正常上皮对照核查

## 本轮问题
公开空间转录组能否检验 LYPLA1 在癌区富集，并高于非恶性上皮？

## 输入与范围
检索病理注释明确且有正常导管的 BRCA 数据。下载 DeepSpaceDB 69 张乳腺癌切片的专家注释，拒绝采用 LLM 标签。先按组织标签筛选，再检查基因覆盖，未按 LYPLA1 结果选择样本。全部源矩阵、测量与区域坐标仅留 server165。

## 实际结果
GSE210616 新找到 1 张有浸润癌及正常导管专家标签的旧队列切片。完成图像坐标核对和 LYPLA1 区域比较，方向为癌区较高。但癌区有 816 个 QC 点，正常导管只有 6 个；要求整个 spot 落在区域内时正常对照只剩 2 个。低于预设每区 20 点，患者统计不可评估，不计算把 spot 当独立患者的 P 值。

GSE242311 有 5 张相关标签切片，其中 2 张含多个组织样本。已下载另 3 张单样本切片的完整矩阵，其 17,943 个特征均不包含 LYPLA1 或 ENSG00000120992，不能用于本基因表达验证。不是零表达，也不是阴性。

其他队列的可获取性及排除原因见 candidate_audit.tsv。此次没有得到新的、满足正常上皮对照要求的多患者验证队列。

## 新手解释
病理标签能说明癌组织和正常导管在哪里，但还必须实际测到 LYPLA1，并有足够正常对照。这次找到的是一张支持方向的示例图，尚不能把“癌上皮高于正常上皮”写成得到充分验证的结论。

## 限制/反证
正常导管区域不是纯化的正常上皮细胞；传统 Visium spot 含混合细胞。GSE210616 为既有队列重新分析，不算独立验证。GSE242311 的患者重复关系尚未核实，且缺目标基因。不能为了显著降低正常点数门槛。

另发现直接先例：[2023 年 GeoMx 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC10463390/)已研究乳腺癌空间区域中的 LYPLA1，并进行蛋白染色及预后分析。因此“首次用乳腺癌空转研究 LYPLA1”不成立；该文不等于已经完成你要的患者内癌上皮—正常上皮比较。

## 当前决定
保留图像作为探索性示例。正式目标状态 NOT_EVALUABLE，而非“显著阳性”或“无差异”。

## 下一步
较贴近目标且值得解决访问问题的是 2024 CAF 研究 Figshare 21591429（正常导管/小叶与癌区）及 SMART 的作者 AOI 元数据。获取后仍须先检查 LYPLA1 覆盖及患者内对照数。没有发送作者邮件或代为申请。

## 复现命令
先在本地运行 code/brca_deepspace_fetch_v1.py，源注释经内存传到服务器。下载器保留服务限流错误；不重试绕过配额，剩余矩阵使用原始 GEO 公开档案 code/brca_deepspace_geo_v1.py。
服务器：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/xuzx/miniconda3/envs/ov/bin/python brca_deepspace_analyze_v1.py RUN_DIR`，随后运行 `brca_deepspace_report_v1.py RUN_DIR`。工作代码提交号见 code_commit.txt。
'''
(P/'README_CN.md').write_text(readme,encoding='utf-8')
(P/'delivery_status.json').write_text(json.dumps({'analysis':'DONE_WITH_UNEVALUABLE_PRIMARY_TARGET','github':'ACCESS_BLOCKED','reason':'HTTPS push failed: TLS close and connection timeout; local commits and server files retained'},indent=2))
(P/'code_commit.txt').write_text('51efae3\n')
print(pd.DataFrame(coverage).to_string(index=False))
