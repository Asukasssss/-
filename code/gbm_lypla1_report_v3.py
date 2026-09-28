"""Render aggregate-only targeted LYPLA1 results; historical all142 files unchanged."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parents[1];RUN='20260928T034033Z_LYPLA1_malignant_v3'
S=R/'results/GBM/06_EXTERNAL'/RUN;O=R/'results/GBM/07_INTEGRATION'/RUN
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')

def main():
    O.mkdir(parents=True,exist_ok=True);figdir=O/'figures';figdir.mkdir(exist_ok=True)
    d=pd.read_csv(S/'LYPLA1_paired_results.tsv',sep='\t')
    studies=[('SNUH','all','SNUH whole-cell RNA\n16 patients; 128,375 cells'),
      ('CARE10x','CARE2025_primary','CARE primary nuclei\n55 patients; 205,880 nuclei'),
      ('CARE10x','CARE2025_recurrent','CARE recurrent nuclei\n59 patients; 222,630 nuclei')]
    comps=['Myeloid','Lymphoid','Oligodendrocyte','Astrocyte','Neuron','OPC','Vascular','All_nonmalignant']
    fig,axes=plt.subplots(1,3,figsize=(16,6),sharey=True)
    colors={'higher':'#bc5b43','lower':'#3273a8','ns':'#626873'}
    for ax,(cohort,partition,title) in zip(axes,studies):
        z=d[d.cohort.eq(cohort)&d.partition.eq(partition)&d.min_cells.eq(20)].set_index('comparator')
        ax.axvline(0,color='#aab1b9',lw=1,zorder=0)
        for y,comp in enumerate(comps):
            a=z.loc[comp]
            if a.status!='DONE':
                ax.text(.04,y,'Not evaluable',transform=ax.get_yaxis_transform(),fontsize=9,color='#90959e');continue
            color=colors['higher' if a.q_value<.05 and a.effect>0 else 'lower' if a.q_value<.05 else 'ns']
            ax.errorbar(a.effect,y,xerr=[[a.effect-a.ci_lower],[a.ci_upper-a.effect]],fmt='o',color=color,capsize=3,ms=6)
            ax.text(.99,y+.27,f'n={int(a.n)}; q={a.q_value:.2g}',transform=ax.get_yaxis_transform(),ha='right',fontsize=8,color=color)
        ax.set_title(title,fontsize=12,pad=13);ax.spines[['top','right']].set_visible(False)
        ax.set_xlabel('Malignant minus comparator\n'+('log2(CPM+1)' if cohort=='SNUH' else 'log1p(CP10k)'),fontsize=10)
        ax.set_yticks(range(len(comps)),[x.replace('_',' ') for x in comps]);ax.set_ylim(7.65,-.65)
        ax.grid(axis='y',alpha=.12);ax.margins(x=.2)
    fig.suptitle('LYPLA1: enrichment depends on the comparison cell type',fontsize=16,y=.99)
    fig.text(.5,.01,'Patient-paired mean differences and 95% patient bootstrap intervals. q: within-family BH Wilcoxon.\nSeparate expression scales; CARE primary/recurrent are the same cohort. Red/blue: q < 0.05; grey: not significant.',ha='center',fontsize=10)
    fig.tight_layout(rect=[0,.09,1,.95]);fig.savefig(figdir/'LYPLA1_paired_enrichment.png',dpi=180,facecolor='white');fig.savefig(figdir/'LYPLA1_paired_enrichment.pdf');plt.close(fig)
    main=d[d.min_cells.eq(20)&((d.cohort.eq('SNUH')&d.partition.eq('all'))|d.cohort.eq('CARE10x'))]
    save(main,S/'LYPLA1_primary_comparisons.tsv')
    rows=[dict(dataset='SNUH',status='DONE',publication_year=2026,technology='whole-cell 10x RNA,including RNA from CITE',
        source_cells=223113,selected_cells=128375,patients=16,malignant_cells=17946,
        main_matched_patients=15,independence='independent author cohort relative to CARE;no individual-level cross-study linkage available',
        quality='author singlet and CNV-backed malignancy;strict GBM IDHwt;patient clinical totals match exactly',
        limits='portal matrices date to 2022,metadata 2024;publication date is not sample collection date;treatment timing missing;35 vascular cells',
        source_url='https://cells.ucsc.edu/multiomic-gbm/scrna/',paper_url='https://www.nature.com/articles/s41467-026-69716-2'),
      dict(dataset='CARE10x_primary',status='DONE',publication_year=2025,technology='10x single-nucleus',source_cells=429305,
        selected_cells=205880,patients=55,malignant_cells=np.nan,main_matched_patients=51,
        independence='reuse previously analyzed CARE cohort;new paired statistical method,not a new dataset',
        quality='author multimodal/CNV annotation;full-library count normalization;untreated-with-RT/TMZ primary subset',
        limits='selected longitudinal reoperation cohort;nuclei versus whole-cell platform differences',
        source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546',paper_url='https://www.nature.com/articles/s41588-025-02168-4'),
      dict(dataset='CARE_Smartseq2',status='NOT_EVALUABLE',publication_year=2025,technology='Smart-seq2 single-nucleus',
        source_cells=14401,selected_cells=np.nan,patients=np.nan,malignant_cells=np.nan,main_matched_patients=np.nan,
        independence='technical complement from CARE;not independent replication',quality='expression matrix successfully downloaded and inspected',
        limits='14424 metadata rows have empty cell-type field;RDS has no cell annotations;verified malignant/CNV labels not found in inspected files',
        source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274548',paper_url='https://www.nature.com/articles/s41588-025-02168-4')]
    save(pd.DataFrame(rows),S/'dataset_quality_registry.tsv')
    report='''# LYPLA1：新 GBM 单细胞队列的恶性细胞富集检验

## 本轮问题
寻找质量更可靠、细胞量足够的独立单细胞队列，判断 LYPLA1 是否在 GBM 恶性细胞富集。GBM 为胶质来源肿瘤，本报告使用“恶性胶质/肿瘤细胞”，不称“癌上皮”。

## 输入与范围
新增 SNUH 图谱（Nature Communications，2026-02-20）：完整细胞 10x RNA，含 CITE-seq 的 RNA。作者恶性标签由 inferCNV 支持；复用作者注释，未重跑 CNV。门户总量223,113含其他疾病与正常脑，不能直接视作 GBM。本次严格选择 GBM、IDH1/2 wildtype、Singlet、肿瘤区域后，保留 **16患者、128,375细胞、17,946恶性细胞**。15患者的恶性细胞达到至少20细胞；进入主配对分析的恶性细胞17,930个。恶性细胞中位检测2,125基因、4,356 UMI。该队列增加完整细胞和独立患者来源，但并非每项质量指标均优于CARE。

患者身份用作者临床表核对；重复解离/测序制备先按临床患者合并，不能算独立患者。细胞条码与表达顺序223,113条完全一致。标准化为作者log2(CPM+1)，用非目标基因RPLP0及原始库大小检查，99.78%阳性细胞的反变换接近整数计数；未对整合矩阵做差异检验，未再次归一化。源矩阵2022年、元数据2024年，2026是发表年份，不是采样年份。

同时复用CARE 2025初发55患者205,880核及复发59患者222,630核，新增患者内配对比较；初发与复发属于同一队列，不能算两次独立重复。其尺度log1p(CP10k)与SNUH不同，效应数值不跨平台直接比较。

## 实际结果
主分析每患者/类别至少20细胞，至少10配对患者才检验；患者等权，双侧配对Wilcoxon，每队列/分区/门槛预设8比较，BH仅对可评估比较校正；10,000次患者配对bootstrap给出均值差区间。区间未作多重比较校正，且均值bootstrap与Wilcoxon检验统计量不同，不能用区间代替q判定。下表效应=恶性细胞减比较细胞；不是倍数。

|队列与比较|配对患者|均值差（95%区间）|q|解释|
|---|---:|---|---:|---|
|SNUH：恶性 vs 髓系|15|+0.263（−0.041～+0.613）|0.303|未检出相对富集；仅8/15患者更高|
|SNUH：恶性 vs 淋巴|15|+0.905（+0.606～+1.287）|0.000122|相对富集，15/15方向一致|
|SNUH：恶性 vs 少突胶质|15|+1.684（+1.457～+1.999）|0.000122|相对富集，15/15方向一致|
|SNUH：恶性 vs 合并非恶性|15|+0.456（+0.141～+0.838）|0.0639|未通过FDR；P=0.0479不能代替q|
|CARE初发：恶性 vs 髓系|51|−0.198（−0.236～−0.158）|1.87×10⁻⁸|恶性细胞更低，46/51患者方向如此|
|CARE复发：恶性 vs 髓系|54|−0.198（−0.230～−0.167）|1.27×10⁻⁹|同队列补充，也更低|

SNUH主分析4个比较可检验，其余4个保留不可评估：无单列正常Astrocyte/Neuron/OPC类别，Vascular仅35细胞、无人达到20细胞。不能将恶性AC-like/OPC-like细胞误作正常胶质对照。

SNUH患者等权描述均值：恶性2.013、髓系1.724、淋巴1.103、少突0.329；检测比例分别27.1%、22.7%、13.2%、4.3%。这些描述使用各类别合格患者集合，不能替代上述同患者配对统计。

稳健性：提高到每类别50细胞，主结论不变；仅核心区髓系比较q=0.599；剔除重复技术制备q=0.252；额外限制UMI≥1000且基因≥500后q=0.762。仅Fresh可配对8患者，未达到预设10患者门槛，不做显著性结论。不能降低门槛追求阳性。

## 新手解释
**LYPLA1不是不在癌细胞表达，而是当前证据不支持它在GBM恶性细胞中稳定、普遍或特异富集。** 新完整细胞队列支持其相对淋巴和少突细胞更高；它与髓系之间的优势不稳定。CARE单核队列反而以髓系更高，因此目前属于明显依赖比较对象和技术/采样背景的表达模式。

## 限制/反证
新SNUH不是初发治疗未暴露专门队列：缺治疗时点，不能推定。该队列以髓系细胞占多数，稀有血管细胞覆盖很弱；不同保存/取样方式与平台可影响检出。RNA表达不等同LYPLA1蛋白、酶活或代谢通量；靶基因系既有结果后选定，P/q按探索性定向复核解释。独立作者队列来自不同研究机构，但无法用公共匿名ID做个体级跨研究排重。

补充筛查CARE Smart-seq2：已获得14,401核表达矩阵；14,424行GEO元数据的细胞类型全部空白，RDS不含恶性/CNV注释，检查的作者仓库对应10x。当前标记NOT_EVALUABLE，不将缺注释当阴性，也不凭LYPLA1自身表达来定义恶性细胞。它原本也是同CARE患者的技术补充，并非独立验证。

## 当前决定
不写“LYPLA1在癌上皮特异高表达”。可写：**“在独立GBM完整细胞队列中，LYPLA1在恶性细胞的表达高于淋巴细胞和少突胶质细胞，但与髓系细胞的差异未达显著；结合CARE单核结果，尚不支持跨队列稳定的恶性细胞富集。”**

历史01_CAMP→02_MAPPING→03_PATIENT→04_ROBUSTNESS→05_FUNCTION结果不变；本轮新增06_EXTERNAL与07_INTEGRATION定向补充。23代谢特征、171关系、142基因的既有全候选交付保留，未因定向LYPLA1分析改写排名或bulk P/q。此前LYPLA1-油酸患者关联仍未通过FDR，不能用单细胞表达补成机制证据。

## 下一步
若拟开展机制实验，应同时纳入恶性细胞与髓系模型，而非先假设癌细胞特异。区分蛋白/活性与RNA；本次没有执行或声称功能验证。

## 复现命令
服务器专用运行目录：`/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260928T034033Z_LYPLA1_malignant_v3`。

在server165准备source数据及预设spec后运行：`gbm_lypla1_source_v3.py <run> qc` → `gbm_lypla1_source_v3.py <run> extract` → `gbm_lypla1_compare_v3.py <run> <CARE_v2_run>` → `gbm_lypla1_validate_v3.R <run>` → `gbm_lypla1_validate_v3.py <run>`。可选下载脚本gbm_lypla1_download_v3.py；预设文件code/gbm_lypla1_v3_spec.json复制为source/analysis_spec_pre_expression.json。验证使用第二次独立元数据下载。

仅复制public汇总目录到本机，再运行`gbm_lypla1_report_v3.py`生成报告与图。细胞/患者级资料仅在server165。R独立重算效应、P、BH q与Python最大误差均<3×10⁻¹⁶；原细胞汇总与患者表一致。输入与代码哈希见source_manifest.tsv。

来源：[SNUH原文](https://www.nature.com/articles/s41467-026-69716-2)；[作者UCSC表达与注释](https://cells.ucsc.edu/multiomic-gbm/scrna/)；[CARE原文](https://www.nature.com/articles/s41588-025-02168-4)；[CARE 10x](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546)；[Smart-seq2筛查](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274548)。
'''
    (O/'README_CN.md').write_text(report,encoding='utf-8')
    (S/'README_CN.md').write_text(report,encoding='utf-8')
    (R/'docs/GBM_LYPLA1_MALIGNANT_V3_CN.md').write_text(report,encoding='utf-8')
    with pd.ExcelWriter(O/'LYPLA1_新单细胞配对结果.xlsx',engine='openpyxl') as w:
        for file,sheet in [('LYPLA1_paired_results.tsv','所有配对比较'),('LYPLA1_celltype_profiles.tsv','患者等权来源描述'),('dataset_quality_registry.tsv','队列筛选与限制'),('SNUH_quality_before_expression.tsv','目标分析前质量核对')]:
            pd.read_csv(S/file,sep='\t').to_excel(w,sheet_name=sheet,index=False)
    idx=R/'coordination/stages/GBM.tsv';t=pd.read_csv(idx,sep='\t',keep_default_na=False);t=t[t.run_id.ne(RUN)]
    add=[]
    for stage in ['06_EXTERNAL','07_INTEGRATION']:
        add.append(dict(cancer='GBM',stage_id=stage,run_id=RUN,analysis_version='gbm_lypla1_malignant_v3',status='DONE',
            scope='targeted LYPLA1 malignant enrichment;SNUH16patients128375cells;CARE paired reanalysis',
            result_path=f'results/GBM/{stage}/{RUN}',code_path='code/gbm_lypla1_compare_v3.py',git_branch='analysis/gbm-discovery-20260925',
            reason='no_general_malignant_enrichment;SS2_NOT_EVALUABLE;historical_all142_bulk_unchanged',next_action='see_GBM_LYPLA1_MALIGNANT_V3_CN'))
    save(pd.concat([t,pd.DataFrame(add)]).sort_values(['stage_id','run_id']),idx)
    # A single canonical numerical table lives in 06_EXTERNAL; integration references it.
    (O/'analysis_spec.json').write_text(json.dumps(dict(analysis_version='gbm_lypla1_malignant_v3',numerical_source=f'results/GBM/06_EXTERNAL/{RUN}',scope='targeted interpretation;historical all142 preserved'),indent=2))
    for filename in ['source_manifest.tsv','validation.json']:
        (O/filename).write_bytes((S/filename).read_bytes())
    print('Report, workbook and forest figure created.')

if __name__=='__main__':main()
