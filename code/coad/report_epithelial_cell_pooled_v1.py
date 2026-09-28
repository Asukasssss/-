"""Render a specified gene's public, cell-weighted CNA summary."""
import argparse,csv,json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'runtime/report_dependencies'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out
with (out/'results.tsv').open(encoding='utf-8') as f:r=list(csv.DictReader(f,delimiter='\t'))
gene=r[0]['gene'];assert len(r)==3 and all(x['gene']==gene for x in r)
assert [x['group'] for x in r]==['tumor_CNA','tumor_CNN','normal_reference']
labels=['肿瘤CNA上皮','肿瘤CNN上皮','正常组织参照上皮'];colors=['#df7e4d','#649ba9','#8eaa68']
plt.rcParams['font.family']='Microsoft YaHei'
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8))
for ax,key,title,scale in [(axes[0],'effect','整体平均表达\n每细胞 log1p(CP10K)',1),(axes[1],'detection_fraction','检测到'+gene+'的细胞比例',100)]:
    vals=[float(x[key])*scale for x in r];ax.bar(range(3),vals,color=colors,width=.6)
    for i,v in enumerate(vals):ax.text(i,v+max(vals)*.035,f'{v:.1f}%' if scale==100 else f'{v:.3f}',ha='center',fontsize=12)
    ax.set_xticks(range(3),labels,fontsize=10);ax.set_ylim(0,max(vals)*1.28);ax.set_title(title,fontsize=12)
    ax.spines[['top','right']].set_visible(False)
fig.suptitle(gene+'整体细胞对比｜Uhlitz队列',fontsize=17)
fig.text(.03,.045,'每个细胞等权；三组细胞数分别为4,477、7,885、16,138。CNN不等于已证实非恶性。',fontsize=11)
fig.subplots_adjust(left=.065,right=.98,top=.76,bottom=.20,wspace=.25)
fig.savefig(out/(gene+'_overall_cells.png'),dpi=180);fig.savefig(out/(gene+'_overall_cells.pdf'));plt.close(fig)
tab='\n'.join(f"|{lab}|{x['n']}|{float(x['effect']):.4f}|{float(x['mean_CP10K']):.4f}|{x['detected_cells']}（{float(x['detection_fraction'])*100:.1f}%）|{float(x['positive_cell_mean_log1p_CP10K']):.4f}|" for lab,x in zip(labels,r))
order=' ＞ '.join(labels[k] for k in sorted(range(3),key=lambda k:float(r[k]['effect']),reverse=True))
ratio1=float(r[0]['mean_CP10K'])/float(r[1]['mean_CP10K']);ratio2=float(r[0]['mean_CP10K'])/float(r[2]['mean_CP10K'])
txt=f'''# COAD：{gene}上皮整体细胞对比

## 本轮问题

沿用UCKL1刚完成的整体细胞口径，查看{gene}在作者CNA/CNN及正常组织参照上皮中的表达。整体平均表达排序为：{order}。

## 输入与范围

GSE166555既有全基因UMI分母及{gene}独立条目，与原作者CNA调用连接。完全沿用[标签核查](../20260925T140257Z_uckl1_cna_v1/README_CN.md)：作者明确规则转换细胞ID，排除一位来源编号冲突供者，正常参照排除作者TC标签。三组仍为4,477、7,885、16,138个细胞，来自9、11、10位供者。输入哈希与上一批一致，不按本基因表达选择细胞。

## 实际结果

|组别|细胞数|平均log1p(CP10K)|平均CP10K|检出细胞数及比例|仅检出细胞的平均log1p(CP10K)|
|---|---:|---:|---:|---|---:|
{tab}

![整体细胞对比]({gene}_overall_cells.png)

在线性CP10K尺度，CNA整体均值约为CNN的{ratio1:.2f}倍、正常参照的{ratio2:.2f}倍。该比值只表示采样细胞的归一化表达均值，不是酶活或患者效应。全细胞表达中位数依次为{float(r[0]['median_log1p_CP10K']):.4f}、{float(r[1]['median_log1p_CP10K']):.4f}、{float(r[2]['median_log1p_CP10K']):.4f}。

## 新手解释

平均表达将未检出的零值也算在内；检出比例回答多少细胞有非零计数。应同时查看检出比例和仅检出细胞的表达，区分更多细胞检出与检出细胞中表达更强这两种表现。

## 限制/反证

每个细胞等权，细胞较多的供者贡献较大。CNA是作者推断拷贝数异常，CNN不能改称已证实非恶性，正常参照是患者的正常/邻近组织。没有把数千细胞当独立患者进行P/q检验。本批未计算{gene}患者配对推断，不能把UCKL1的患者P/q借给{gene}。检出比例受测序深度及细胞状态影响，不能直接归因为癌变机制。

## 当前决定

保留实际整体表达及检出比例对照，仅用于描述；原候选安排与历史统计保持不变。

## 下一步

将此整体视图用于阅读候选表达背景；若要判断患者间是否稳定，再以该基因自己的患者值分析。

## 复现命令

server165新独占运行目录运行`python3 run_epithelial_cell_pooled_v1.py --out <新目录> --gene {gene}`；本地从公开表运行`python code/coad/report_epithelial_cell_pooled_v1.py --out <结果目录>`。核对直接逐细胞均值与按细胞数加权的患者均值一致、检出数一致、直方图细胞数守恒。原始和逐患者/细胞值留服务器。
'''
(out/'README_CN.md').write_text(txt,encoding='utf-8',newline='\n')
(out/'code_manifest.json').write_text(json.dumps([dict(path='code/coad/'+f,sha256=hashlib.sha256((ROOT/'code/coad'/f).read_bytes()).hexdigest()) for f in ['run_epithelial_cell_pooled_v1.py','report_epithelial_cell_pooled_v1.py']],indent=2)+'\n',encoding='utf-8',newline='\n')
print(tab);print('mean CP10K ratios',ratio1,ratio2)
