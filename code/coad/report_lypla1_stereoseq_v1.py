"""Publish aggregate interpretation and fixed stage record; never read local source matrices."""
import json,hashlib,subprocess,zipfile
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];RUN='20260928T103900Z_lypla1_stereoseq_v1';REL=Path('results/COAD/06_EXTERNAL')/RUN;OUT=ROOT/REL
v=json.load(open(OUT/'validation.json'));assert v['status']=='PASS'
v['patient_linkage_verified']=False
v['limitations']=[x for x in v['limitations'] if 'Two CRCP59' not in x]+['Similarly named CRCP59 sections have unresolved patient linkage; section count is not patient count']
(OUT/'validation.json').write_bytes((json.dumps(v,indent=2)+'\n').encode())
(OUT/'catalog.json').write_bytes((ROOT/'runtime/stomics_catalog.json').read_bytes())
d=pd.read_csv(OUT/'region_summary.tsv',sep='\t');w=d.pivot(index='section',columns='region',values='mean_expression')
t='Tumor-associated epithelium';n='Normal epithelium';s='Stroma'
a=w[[t,n]].dropna();b=w[[t,s]].dropna()
lines=['|切片|癌相关上皮|正常上皮|基质|','|---|---:|---:|---:|']
for sid,row in w.iterrows():lines.append('|'+sid+'|'+'|'.join('NA'if pd.isna(row.get(k))else '%.3f'%row[k]for k in[t,n,s])+'|')
report=f"""# LYPLA1：更换为STT0000036 Stereo-seq区域表达复核

## 本轮问题
用户要求换其他空间数据，重点解决原图中癌区难以辨识的问题。新批次分析LYPLA1在作者已注释癌相关上皮、正常上皮及基质区域的位置和表达。旧队列及原CAMP统计保持不变。

## 输入与范围
- 原研究：Feng等，Nature Communications 2024，DOI [10.1038/s41467-024-54710-3](https://www.nature.com/articles/s41467-024-54710-3)。[公开数据STT0000036](https://db.cngb.org/stomics/project/STT0000036)。
- 原研究Stereo-seq共16张切片、15名患者。本批先依据来源Treatment=naive纳入全部10张未经治疗切片，6张pMMR、4张dMMR；CRCP59_T和CRCP59_T_2的来源样本isolate不同，但命名相近且论文总样本数多于患者数，具体重复患者对应尚未取得；本批不据编号推断独立性。6张PD-1治疗切片单列范围台账，不进入本批表达比较。
- 本批实际{{v['bins']:,}}个空间单元，分辨率约50×50微米，可含多种细胞。CRC材料，不扩写为纯结肠癌。
- 原始区域字段level2完整保留：epi_CEA、epi_MKI67合并展示为癌相关上皮，epi_normal为正常上皮。论文称相关簇为tumor_CEA、tumor_MKI67，依据表达、组织图像并辅以CNV检查；这里复用作者结果，没有重新做病理诊断。
- 直接读取作者GEM导出的MIDCounts。该字段虽然叫counts，数值实际是变换后的连续值，不能当原始UMI。全部空间点逐项检查反变换后的总和接近10000、按nCount_Spatial还原后接近整数、非零条目数等于nFeature_Spatial，确认与log1p(CP10K)相容。展示保持源数值，不再二次归一化、填补或批次校正。

## 实际结果
全部10张切片均有LYPLA1条目。以下为每个区域的平均log1p(CP10K)，不是倍数变化。完整点数、检出比例、中位数及原细分标签见TSV。

{{chr(10).join(lines)}}

- 同时含癌相关上皮和正常上皮的{{len(a)}}张切片中，{{int((a[t]>a[n]).sum())}}张癌相关上皮平均表达较高。
- 同时含癌相关上皮和基质的{{len(b)}}张切片中，{{int((b[t]>b[s]).sum())}}张癌相关上皮平均表达较高。
- 以上是切片内描述性方向计数，不是独立患者检验，不新增P/q，不据此宣称普遍上调或机制成立。

## 新手解释
优先打开`figures/CRCP59_T_LYPLA1_CN.png`：这一示例在查看表达前已按论文代表性未经治疗pMMR切片选定，未按LYPLA1高低挑选。
1. 左图红色是作者注释的癌相关上皮；蓝色是正常上皮；黄色为肿瘤与基质交界区。
2. 右图与左图位置完全对应。黑线标出左图癌相关上皮的位置；填充色从浅灰到深红表示LYPLA1从低到高。
3. 先找到黑线内，再比较线内、正常上皮和基质的颜色。右图出现红色不自动说明该处是癌区。
4. 全10张图使用共同色标，不裁掉高值、不平滑、不按每张切片单独拉伸颜色；按源坐标展示，没有H&E底图。

## 限制与反证
区域由作者结合多种信息命名，不等于每个空间点都是纯恶性细胞。LYPLA1表达不能反向定义癌区；使用的区域标签也不是本项目独立病理复核。来源平台本次导出没有提供可直接对齐的组织染色图，本批不制作假的H&E叠加图。切片内点数不当作患者样本量；患者对应未完全核实，不把10张切片写成10名独立患者。本批未检验跨队列患者重叠，也不是CAMP代谢关系验证。

另核查GSE294385原发肿瘤M-ST-13：18085项特征中，LYPLA1符号和ENSG00000120992均未匹配，因此未下载该样本完整表达矩阵。这个结论只适用于实际检查的这一份特征表，不扩大为该研究全部样本缺测。

## 当前决定
本批收口为新的空间表达背景补充。保留区域内完整分布及所有切片，不以表达方向挑选结果，不替代前一批空间数据或全候选结果。

## 下一步
先据中文双面板图阅读区域和表达；如需病理图叠加，应取得这批切片实际配套图像和变换参数后另建批次。

## 复现命令
在server165本次独占运行目录执行：
```sh
python3 acquire_lypla1_stereoseq_v1.py catalog.json --expression
python3 analyze_lypla1_stereoseq_v1.py
python3 render_lypla1_stereoseq_cn_v1.py
```
原矩阵和逐点数据仅在服务器；公开交付包含汇总表、图、参数、来源哈希和校验记录。串行首轮因执行效率改为三进程处理，终止日志保留；发布结果以最终完整运行校验为准。
""".replace("{{","{").replace("}}","}")
# Expand only the intended numeric placeholders after building the template.
report=report.replace("{v['bins']:,}",format(v['bins'],',')) .replace('{chr(10).join(lines)}','\n'.join(lines)).replace('{len(a)}',str(len(a))).replace('{len(b)}',str(len(b))).replace('{int((a[t]>a[n]).sum())}',str(int((a[t]>a[n]).sum()))).replace('{int((b[t]>b[s]).sum())}',str(int((b[t]>b[s]).sum())))
(OUT/'README_CN.md').write_bytes(report.encode('utf8'))
cols=(ROOT/'templates/statistical_result.tsv').read_text(encoding='utf8').splitlines()[0].split('\t');out=[]
for row in d.itertuples():
 r={k:None for k in cols};r.update(cancer='COAD',cohort='STT0000036',stage_id='06_EXTERNAL',run_id=RUN,analysis_version='COAD_LYPLA1_Stereo_seq_v1',analysis_type='spatial_region_descriptive',gene='LYPLA1',unit='50um spatial bin',n=row.n_bins,effect_type='mean_source_log1pCP10K',effect=row.mean_expression,test_family='NOT_RUN',family_n_evaluable=0,status='DONE',reason='Descriptive only; no independent-patient inference',source_id=row.section,region=row.region,detection_fraction=row.detection_fraction,n_detected=row.n_detected);out.append(r)
pd.DataFrame(out,columns=cols+['region','detection_fraction','n_detected']).to_csv(OUT/'results.tsv',sep='\t',index=False,na_rep='NA')
spec=json.load(open(OUT/'analysis_spec.json'));spec['statistical_unit']='Within-section spatial bins, descriptive; patient linkage unresolved for similarly named sections';spec['base_git_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();spec['code_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in (ROOT/'code/coad').glob('*lypla1_stereoseq*v1.py')};(OUT/'analysis_spec.json').write_bytes((json.dumps(spec,indent=2)+'\n').encode())
idx=ROOT/'coordination/stages/COAD.tsv';x=pd.read_csv(idx,sep='\t',dtype=str,keep_default_na=False);assert RUN not in set(x.run_id)
r=dict(cancer='COAD',stage_id='06_EXTERNAL',run_id=RUN,analysis_version='COAD_LYPLA1_Stereo_seq_v1',status='DONE',scope=f"10 untreated sections; {v['bins']} bins; author region labels and LYPLA1; descriptive only",result_path=REL.as_posix(),code_path='code/coad/analyze_lypla1_stereoseq_v1.py',git_branch='analysis/coad-initial',reason='Explicit coordinate join; source log1pCP10K verified; repeated sections kept separate; no P/q',next_action='Read Chinese region/expression panels; no mechanism inference')
x=pd.concat([x,pd.DataFrame([r])],ignore_index=True).sort_values(['stage_id','run_id']);idx.write_bytes(x.to_csv(sep='\t',index=False).encode())
p=ROOT/'reports/COAD/CURRENT_STATUS_CN.md';old=p.read_text(encoding='utf8');p.write_bytes((f"## 2026-09-28｜LYPLA1更换空间数据：STT0000036\n\n已完成全部10张未经治疗切片、{v['bins']:,}个空间单元的作者区域与LYPLA1表达描述；原标签、重复切片和治疗范围均保留，无新增P/q。逐点数据仅存165。[本批报告](../../{REL.as_posix()}/README_CN.md)。\n\n"+old).encode())
print('REPORT_READY',v['bins'],len(a),int((a[t]>a[n]).sum()),len(b),int((b[t]>b[s]).sum()))
