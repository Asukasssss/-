"""Public aggregate-only overall UCKL1 view."""
import csv,json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'runtime/report_dependencies'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'results/COAD/06_EXTERNAL/20260925T141710Z_uckl1_cell_pooled_v1'
with (OUT/'results.tsv').open(encoding='utf-8') as f:r=list(csv.DictReader(f,delimiter='\t'))
assert [x['group'] for x in r]==['tumor_CNA','tumor_CNN','normal_reference']
labels=['肿瘤CNA上皮','肿瘤CNN上皮','正常组织参照上皮'];colors=['#df7e4d','#649ba9','#8eaa68']
plt.rcParams['font.family']='Microsoft YaHei'
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8))
for ax,key,title,scale in [(axes[0],'effect','整体平均表达\n每细胞 log1p(CP10K)',1),(axes[1],'detection_fraction','检测到UCKL1的细胞比例',100)]:
    vals=[float(x[key])*scale for x in r]
    ax.bar(range(3),vals,color=colors,width=.6)
    for i,v in enumerate(vals):ax.text(i,v+max(vals)*.035,f'{v:.1f}%' if scale==100 else f'{v:.3f}',ha='center',fontsize=12)
    ax.set_xticks(range(3),labels,fontsize=10);ax.set_ylim(0,max(vals)*1.28);ax.set_title(title,fontsize=12)
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('UCKL1整体细胞对比｜Uhlitz队列',fontsize=17)
fig.text(.03,.045,'每个细胞等权；三组细胞数分别为4,477、7,885、16,138。CNN不等于已证实非恶性。',fontsize=11)
fig.subplots_adjust(left=.065,right=.98,top=.76,bottom=.20,wspace=.25)
fig.savefig(OUT/'UCKL1_overall_cells.png',dpi=180);fig.savefig(OUT/'UCKL1_overall_cells.pdf');plt.close(fig)
tab='\n'.join(f"|{lab}|{x['n']}|{float(x['effect']):.4f}|{float(x['mean_CP10K']):.4f}|{x['detected_cells']}（{float(x['detection_fraction'])*100:.1f}%）|" for lab,x in zip(labels,r))
ratio1=float(r[0]['mean_CP10K'])/float(r[1]['mean_CP10K']);ratio2=float(r[0]['mean_CP10K'])/float(r[2]['mean_CP10K'])
txt=f'''# COAD：UCKL1全部细胞等权整体对比

## 本轮问题

按用户要求直接看整体细胞分布，不按患者等权进行主展示。使用同一批已核实CNA/CNN标签及同一供者冲突排除规则。

## 输入与范围

沿用上一批CNA分层的三组细胞：CNA 4,477、CNN 7,885、正常组织参照16,138。源数据与标签哈希逐项相同。每个细胞等权，来自细胞较多的供者会占更大权重；不改变原始UMI、分组或归一化尺度。

## 实际结果

|细胞组|细胞数|平均log1p(CP10K)|平均CP10K|检出细胞数及比例|
|---|---:|---:|---:|---|
{tab}

![全部细胞等权对比](UCKL1_overall_cells.png)

整体描述上，UCKL1在CNA组的平均表达和检出比例较高。在线性CP10K尺度，CNA组平均值约为CNN组的{ratio1:.2f}倍、正常参照组的{ratio2:.2f}倍；这是所采样细胞的归一化表达均值比，不是患者效应倍数或酶活倍数。三组单细胞表达中位数都为0，因此不是所有细胞均普遍高表达。

进一步限定到检测到UCKL1的细胞，平均log1p(CP10K)分别为{float(r[0]['positive_cell_mean_log1p_CP10K']):.3f}、{float(r[1]['positive_cell_mean_log1p_CP10K']):.3f}、{float(r[2]['positive_cell_mean_log1p_CP10K']):.3f}。CNA组在这项条件性描述中并不更高，因此整体均值较高主要伴随更高的检出比例，不能改写成每个表达阳性的细胞都更高。检出也受测序深度等技术因素影响。

## 新手解释

这张图把各组细胞放在一起，回答这批细胞总体看起来谁更高。上一批每位患者等权，回答患者之间是否一致。两种加权方式回答的问题不同，整体均值较高可以与患者配对检验不显著并存。

## 限制/反证

CNA为作者表达推断的拷贝数异常；CNN不是已证实非恶性。整体描述不作把细胞视为独立患者的P/q检验。原患者比较仍为CNA对CNN及正常参照各7/9位较高、P=0.1796875、q=0.26953125；此次不替换它们。

## 当前决定

增加整体细胞视图，保留“CNA组表达及检出较高”的描述性观察，不升级为恶性特异或显著患者效应。

## 下一步

整体展示与[患者配对报告](../20260925T140257Z_uckl1_cna_v1/README_CN.md)一并使用，按研究问题选择解释口径。

## 复现命令

server165新运行目录建立独占`.running`，复制`run_uckl1_cell_pooled_v1.py`，运行`python3 run_uckl1_cell_pooled_v1.py --out <新目录>`。本地从公开汇总运行`python code/coad/report_uckl1_cell_pooled_v1.py`绘图。直接逐细胞均值与按细胞数加权的既有患者均值交叉核对，检出数、细胞数、直方图守恒均PASS；细胞及患者级数据留服务器。
'''
(OUT/'README_CN.md').write_text(txt,encoding='utf-8',newline='\n')
(OUT/'code_manifest.json').write_text(json.dumps([dict(path='code/coad/'+f,sha256=hashlib.sha256((ROOT/'code/coad'/f).read_bytes()).hexdigest()) for f in ['run_uckl1_cell_pooled_v1.py','report_uckl1_cell_pooled_v1.py']],indent=2)+'\n',encoding='utf-8',newline='\n')
print(tab);print('linear mean ratios',ratio1,ratio2)
