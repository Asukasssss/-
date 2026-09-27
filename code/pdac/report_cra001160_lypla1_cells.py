"""Public plots/report for exploratory cell-level LYPLA1 analysis."""
from pathlib import Path
import sys,json,csv,hashlib,subprocess
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'runtime/plot_dependencies'))
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
RUN='20260927T135000Z_cra001160_lypla1_cells_v1';REL='results/PDAC/06_EXTERNAL/'+RUN;P=ROOT/REL
r=pd.read_csv(P/'results.tsv',sep='\t');h=pd.read_csv(P/'histogram_source.tsv',sep='\t');g=pd.read_csv(P/'group_profiles.tsv',sep='\t')
# Independently verify BH, histogram counts and weighted group means.
order=np.argsort(r.p_value.to_numpy());expected=np.minimum.accumulate((r.p_value.to_numpy()[order]*6/np.arange(1,7))[::-1])[::-1].clip(0,1)
assert np.allclose(r.q_value.to_numpy()[order],expected,rtol=1e-10,atol=0)
for _,row in r.iterrows():
 for role,n in [('type2',row.n),('type1',row.n_reference)]:assert h[(h.context==row.context)&(h.population==row.population)&(h.role==role)].n_cells.sum()==n
source=pd.read_csv(P/'source_manifest.tsv',sep='\t').iloc[1]
assert hashlib.sha256((ROOT/'code/pdac/cra001160_lypla1_cells.py').read_bytes()).hexdigest()==source.sha256
font=FontProperties(fname='C:/Windows/Fonts/msyh.ttc');plt.rcParams['font.family']=font.get_name();plt.rcParams['axes.unicode_minus']=False;plt.rcParams['pdf.fonttype']=42
fig,axes=plt.subplots(2,2,figsize=(13,9));fig.subplots_adjust(left=.08,right=.97,top=.82,bottom=.13,hspace=.55,wspace=.25)
fig.suptitle('CRA001160 · LYPLA1 细胞层面比较',x=.07,y=.965,ha='left',fontsize=23,fontweight='bold')
fig.text(.07,.905,'type 2：作者恶性相关导管上皮；type 1：相对非恶性导管上皮。每个细胞等权。',fontsize=12,color='#475569')
for i,context in enumerate(['Tumor_all_units','Tumor_vs_control']):
 for j,pop in enumerate(['All_cells_including_zero','Detected_cells_raw_UMI_gt0']):
  ax=axes[i,j];row=r[(r.context==context)&(r.population==pop)].iloc[0]
  for role,color,label in [('type2','#c65b46','肿瘤 type 2'),('type1','#44839b','肿瘤 type 1' if i==0 else '对照 type 1')]:
   d=h[(h.context==context)&(h.population==pop)&(h.role==role)].sort_values('left')
   x=np.r_[d.left.iloc[0],d.right.to_numpy()];y=np.r_[0,np.cumsum(d.fraction)]
   ax.step(x,y,where='post',color=color,lw=2,label=label)
  title=('肿瘤来源细胞' if i==0 else '肿瘤 vs 对照胰腺')+' · '+('全部细胞（含零）' if j==0 else '仅 LYPLA1 检出细胞')
  ax.set_title(title,loc='left',fontsize=13,pad=12);ax.set_xlim(0,2.5);ax.set_ylim(0,1.03)
  ax.set_xlabel('LYPLA1 标准化表达：log1p(每万 UMI)');ax.set_ylabel('累计细胞比例')
  ax.legend(loc='lower right',fontsize=10,frameon=False)
  ax.text(.97,.46,f'n = {int(row.n):,} / {int(row.n_reference):,}\n平均 = {row.mean_test:.3f} / {row.mean_reference:.3f}\nP = {row.p_value:.3g}；q = {row.q_value:.3g}\n秩效应 = {row.effect:.3f}',transform=ax.transAxes,ha='right',va='top',fontsize=10)
  ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
fig.text(.07,.06,'曲线由完整细胞分箱汇总；显示表达 ≤2.5 的区间。P/q 为未校正患者内相关性的探索性细胞检验。',fontsize=11,color='#475569')
fig.text(.07,.025,'q 覆盖 4 项展示比较及 2 项共同患者子集比较；不能替代患者层面证据。',fontsize=11,color='#475569')
fig.savefig(P/'LYPLA1_cell_level.png',dpi=180);fig.savefig(P/'LYPLA1_cell_level.pdf');plt.close(fig)
labels={'Tumor_all_units':'全部肿瘤单位','Tumor_vs_control':'肿瘤 vs 对照','Tumor_14_shared_units_sensitivity':'共同14个肿瘤单位（细胞检验）'}
table=[]
for _,a in r.iterrows():table.append(f"|{labels[a.context]}|{'全部细胞' if a.population.startswith('All') else '仅检出细胞'}|{int(a.n)} / {int(a.n_reference)}|{int(a.units_test)} / {int(a.units_reference)}|{a.mean_test:.4f} / {a.mean_reference:.4f}|{a.median_test:.4f} / {a.median_reference:.4f}|{a.effect:.4f}|{a.p_value:.5g}|{a.q_value:.5g}|")
text='''# CRA001160：LYPLA1 细胞层面探索性比较

## 本轮问题

按用户要求，直接比较细胞分布，同时展示包含零值的全部细胞及仅 LYPLA1 检出细胞。**肿瘤来源全部细胞的 type 2 与 type 1 差异不显著；仅检出细胞中，type 2 的标准化表达反而较低。**

## 输入与范围

复用上一批原作者矩阵提取的57,530细胞结果，输入 SHA-256 固定且核对通过。原作者 type 2 为恶性相关导管上皮，type 1 为相对非恶性导管上皮；未新增 CNV。表达沿用 log1p(10000×LYPLA1 UMI/全基因 UMI)，阳性定义为原始 LYPLA1 UMI>0。

本轮每个细胞等权，不先按患者求平均。主展示保留全部作者注释导管细胞，不使用每患者20细胞门槛；另保留上一轮14个共同合格患者/标本单位的细胞子集作为敏感性比较。共同患者子集仍用细胞检验，不是患者配对检验。

双侧 Mann–Whitney U（Wilcoxon秩和），渐近法、并列秩及连续性校正；6项预设探索性比较共同BH。秩二列效应=2U/(n1×n2)−1，正值为type 2整体分布较高，负值较低；不是倍数。P检验的是分布/秩，而不是专门检验均值。均值和中位数均另行展示。不提供将细胞当独立重复的置信区间。

## 实际结果

下表均以前组type 2、后组type 1的顺序展示。

|比较|细胞口径|细胞数|来源单位数|平均表达|中位数|秩二列效应|探索性P|探索性q|
|---|---|---:|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(table)+'''

![细胞分布](LYPLA1_cell_level.png)

## 新手解释

肿瘤type 2中5312/11315（46.95%）细胞检出LYPLA1，肿瘤type 1中1097/2646（41.46%）检出，对照type 1中2543/7671（33.15%）检出。这是按细胞计数的比例，不是上一轮患者等权检出比例。

type 2能检出的细胞比例较高，但在已经检出的细胞内，相对表达强度较低。不能只看其中一种口径就写成“LYPLA1在恶性上皮整体升高”或“整体降低”。阳性筛选改变了所比较的人群，且标准化表达仍受到捕获深度、总RNA含量和细胞状态的影响。

## 限制与反证

- 同一患者内细胞不独立，本轮P/q没有做患者聚类校正；非常小的P不等于可跨患者复现。BH不能修复这种依赖。
- 每个细胞等权使细胞较多的患者占更大权重；`group_profiles.tsv`保存最大单一单位及前三单位的细胞占比，未导出逐患者测量。
- 在共同14个单位中，全细胞检验P约0.00147，但上一轮患者层面P=0.60445，14个单位7升7降。这正说明两种统计口径的证据含义不同，不是独立验证成功。
- 肿瘤全部单位比较纳入24个type 2来源单位、22个type 1来源单位；单位构成不完全相同。共同14单位的敏感性保留在完整表中。
- type 1不等于健康细胞；跨组织对照不是同患者配对。未重新推断恶性标签，也未校正检测深度。
- 本次q只属于这6个细胞层面检验，不混用上一轮患者q、CAMP q或全部候选q。

## 当前决定

本批细胞层面展示完成，解释为细胞分布的探索性描述。保留已有患者层面结果，不提升为“恶性特异”或机制证据。

## 下一步

若需要稳健的跨患者结论，应保留患者层面的结果，或另行预设包含患者随机效应的模型；不能用本轮更小的细胞P替代。

## 验证与复现命令

输入哈希、细胞唯一性、原始计数标准化、阳性筛选和所有分箱总数检查通过；用独立 pooled ranks 重建6项U与并列秩校正P通过；公开文件另行独立复核6项BH及脚本哈希。逐细胞和逐患者数据留服务器，仅导出汇总。

在新PDAC/B运行目录复制 `code/pdac/cra001160_lypla1_cells.py`，运行 `python3 cra001160_lypla1_cells.py`；仅读既有私有结果。下载public目录后，在本地运行 `python code/pdac/report_cra001160_lypla1_cells.py` 生成报告与图件。
'''
(P/'README_CN.md').write_text(text,encoding='utf-8')
(P/'delivery_validation.json').write_text(json.dumps({'status':'PASS','independent_BH_six_tests':True,'histogram_totals':True,'server_analysis_code_hash_matches':True,'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()},indent=2))
idx=ROOT/'coordination/stages/PDAC.tsv'
with idx.open(encoding='utf-8') as f:reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
rows=[a for a in rows if a['run_id']!=RUN]
rows.append(dict(cancer='PDAC',stage_id='06_EXTERNAL',run_id=RUN,analysis_version='cra001160_lypla1_cells_v1',status='DONE',scope='LYPLA1 pooled-cell all/positive analyses;6 exploratory contrasts',result_path=REL,code_path='code/pdac/cra001160_lypla1_cells.py',git_branch='analysis/pdac-initial',reason='Donor dependence uncorrected;not patient replication;positive type2 expression lower',next_action='Retain donor-level counterevidence;no malignant-specific claim'))
with idx.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sorted(rows,key=lambda a:(a['stage_id'],a['run_id'])))
print(P/'LYPLA1_cell_level.png')
