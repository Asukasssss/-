"""Preserve actual printed aggregates when server export is temporarily unreachable."""
from pathlib import Path
import csv,json,hashlib,subprocess
ROOT=Path(__file__).resolve().parents[2];RUN='20260927T135000Z_cra001160_lypla1_cells_v1';REL='results/PDAC/06_EXTERNAL/'+RUN;P=ROOT/REL;P.mkdir(parents=True,exist_ok=True)
if (P/'results.tsv').exists():raise SystemExit('Full export exists; use report_cra001160_lypla1_cells.py')
names=['context','population','n','n_reference','mean_test','mean_reference','rank_biserial','p_value','q_value']
data=[['Tumor_all_units','all',11315,2646,.361448,.359252,.018591,1.040993e-1,1.040993e-1],
['Tumor_all_units','positive',5312,1097,.769915,.866527,-.186432,2.116098e-22,3.174148e-22],
['Tumor_vs_control','all',11315,7671,.361448,.327674,.082851,1.827467e-27,5.482401e-27],
['Tumor_vs_control','positive',5312,2543,.769915,.988433,-.354080,1.144984e-142,6.869903e-142],
['Tumor_14_shared_units_sensitivity','all',4019,2616,.371686,.358034,.042525,1.472145e-3,1.766574e-3],
['Tumor_14_shared_units_sensitivity','positive',2028,1081,.736592,.866437,-.234282,4.580406e-27,9.160811e-27]]
with (P/'console_summary.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(names);w.writerows(data)
table='\n'.join('|'+ '|'.join(str(x) for x in a)+'|' for a in data)
text='''# CRA001160 LYPLA1 细胞比较：已计算数值，完整回传待完成

## 本轮问题
直接按细胞比较原作者恶性相关导管type 2与相对非恶性导管type 1；分别纳入全部细胞与仅LYPLA1原始UMI>0的检出细胞。

## 输入与范围
服务器复用上一轮原作者矩阵提取的细胞数值，表达为log1p(10000×LYPLA1 UMI/全基因UMI)。每个细胞等权，不先按患者求均值。包括肿瘤来源全部单位、肿瘤vs对照，以及原先14个共同合格单位的细胞子集。6项双侧Mann–Whitney U检验统一BH，均未校正患者内依赖。

## 实际结果
以下直接记录本轮服务器成功运行后返回的汇总输出，保留终端打印精度；不是完整精度TSV回传。test为type 2，reference为type 1；all含零值，positive仅原始UMI>0。

|比较|口径|细胞数type2|细胞数type1|均值type2|均值type1|秩二列效应|探索性P|探索性q|
|---|---|---:|---:|---:|---:|---:|---:|---:|
'''+table+'''

## 新手解释
仅检出细胞中type 2标准化表达较低，而在所有肿瘤来源细胞（含零）中差异不显著。检出率和阳性细胞表达强度是不同问题，不能用一个指标概括所有结果。

## 限制与反证
这些P/q将细胞作为独立观测，没有患者聚类校正，不能视为患者层面的独立验证。共同14个单位中的全细胞P≈0.00147，与先前患者层面P=0.60445不同，体现单位及权重改变；历史患者结果保留。type 1不等于健康细胞，未新做CNV。平均表达不是倍数；秩检验不是专门的均值检验。

## 当前决定
服务器实际计算已完成，6行汇总已返回。随后SFTP连接超时，完整结果、分布分箱、运行验证JSON尚未回传；本地交付状态PARTIAL，不能声称完整图表已交付。公开此临时汇总与代码，以便接续。

## 下一步
连接恢复后从同一运行目录public回传原结果及图源数据，运行正式报告程序；无需重跑矩阵或改变统计规则。

## 复现命令
服务器新运行目录执行 `python3 cra001160_lypla1_cells.py`。现有计算目录为 `results/collaborative/PDAC/B/20260927T135000Z_cra001160_lypla1_cells_v1`；正式本地图表脚本为 `code/pdac/report_cra001160_lypla1_cells.py`，需先回传public目录。
'''
(P/'README_CN.md').write_text(text,encoding='utf-8')
(P/'validation.json').write_text(json.dumps({'status':'PARTIAL','server_computation':'Completed; successful console table received','input_and_statistics_assertions':'All assertions preceding printed summary passed; source validation JSON not yet transferred','export':'SSH timeout on subsequent SFTP connections','precision':'Terminal display precision only','full_results_and_plot':'Pending transfer'},indent=2))
(P/'analysis_spec.json').write_text(json.dumps({'analysis_version':'cra001160_lypla1_cells_v1','unit':'pooled_cell','populations':['all including zero','raw LYPLA1 UMI>0'],'test':'two-sided asymptotic Mann-Whitney U, tie and continuity correction','family':'BH across six predefined comparisons','clustering_adjustment':False,'input_sha256':'39c164112b49c9807e34b378b912f7ad510e3faacc63ed415f28f624960c8d8c','base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()},indent=2))
with (P/'source_manifest.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['source','sha256','note']);w.writerow(['Prior server private/cell_values.tsv.gz','39c164112b49c9807e34b378b912f7ad510e3faacc63ed415f28f624960c8d8c','Verified before successful cell analysis']);w.writerow(['code/pdac/cra001160_lypla1_cells.py',hashlib.sha256((ROOT/'code/pdac/cra001160_lypla1_cells.py').read_bytes()).hexdigest(),'Executed server code']);w.writerow(['console_summary.tsv',hashlib.sha256((P/'console_summary.tsv').read_bytes()).hexdigest(),'Copied actual printed aggregates; rounded display precision'])
idx=ROOT/'coordination/stages/PDAC.tsv'
with idx.open(encoding='utf-8') as f:d=csv.DictReader(f,delimiter='\t');fields=d.fieldnames;rows=list(d)
rows=[a for a in rows if a['run_id']!=RUN];rows.append(dict(cancer='PDAC',stage_id='06_EXTERNAL',run_id=RUN,analysis_version='cra001160_lypla1_cells_v1',status='PARTIAL',scope='Six cell-level LYPLA1 contrasts computed;console summary delivered',result_path=REL,code_path='code/pdac/cra001160_lypla1_cells.py',git_branch='analysis/pdac-initial',reason='Full public export blocked by subsequent SSH timeouts;no patient-level independence claim',next_action='Retrieve existing server public files and render plots;do not rerun matrices'))
with idx.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sorted(rows,key=lambda a:(a['stage_id'],a['run_id'])))
print(P)
