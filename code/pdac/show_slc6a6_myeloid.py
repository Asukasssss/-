"""Display existing public source aggregates; no new inference or patient data."""
import csv, json, hashlib, subprocess, sys
from pathlib import Path
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'runtime/plot_dependencies'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

run=sys.argv[1] if len(sys.argv)>1 else datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_SLC6A6_myeloid_display'
assert '/' not in run and '\\' not in run and run.endswith('_SLC6A6_myeloid_display')
out=ROOT/'results/PDAC/06_EXTERNAL'/run
out.mkdir(parents=True,exist_ok=len(sys.argv)>1)
base=ROOT/'results/PDAC/07_INTEGRATION/20260922T140000Z_unified_report_v5/appendix'
cohorts=['GSE263733','GSE278688','GSE242230']
rows=[]; manifest=[]
for cohort in cohorts:
    p=base/(cohort+'.tsv')
    source=list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
    selected=[r for r in source if r['gene']=='SLC6A6']
    assert len({r['celltype'] for r in selected})==len(selected)
    assert sum(r['celltype']=='Myeloid' for r in selected)==1
    rows.extend(selected)
    manifest.append({'path':p.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
def write_tsv(name,data):
    with (out/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
write_tsv('all_celltype_source.tsv',rows)
write_tsv('myeloid_source.tsv',[r for r in rows if r['celltype']=='Myeloid'])
write_tsv('source_manifest.tsv',manifest)
font=FontProperties(fname='C:/Windows/Fonts/msyh.ttc')
plt.rcParams['font.family']=font.get_name()
plt.rcParams['axes.unicode_minus']=False
labels={'Myeloid':'髓系细胞','Endothelial':'内皮细胞','Fibroblast/CAF':'成纤维细胞 / CAF','Ductal (unresolved)':'导管（身份未定）','Malignant (author)':'恶性（作者标注）','Normal epithelial (author)':'正常上皮（作者标注）','B/Plasma':'B / 浆细胞','T/NK':'T / NK','Acinar':'腺泡细胞','Erythrocyte':'红细胞'}
fig,axes=plt.subplots(1,3,figsize=(16,7.8))
fig.subplots_adjust(left=.125,right=.97,top=.77,bottom=.32,wspace=.85)
fig.suptitle('PDAC · SLC6A6 的髓系表达来源',fontsize=23,x=.04,ha='left',y=.97,fontweight='bold')
fig.text(.04,.89,'两套队列髓系均值排名第一；第三套排名第二，来源排序不稳定',fontsize=14,color='#334155')
summ=[]
for ax,c in zip(axes,cohorts):
    rr=[r for r in rows if r['cohort']==c and r['effect']!='NA']
    rr.sort(key=lambda r:float(r['effect']),reverse=True)
    for i,r in enumerate(rr):
        color='#be4538' if r['celltype']=='Myeloid' else '#94a3b8'
        ax.plot([float(r['q25']),float(r['q75'])],[i,i],color=color,lw=5,alpha=.55,solid_capstyle='round')
        ax.scatter(float(r['effect']),i,s=65,color=color,zorder=3,edgecolor='white',linewidth=.7)
    ax.set_yticks(range(len(rr)),[labels.get(r['celltype'],r['celltype']) for r in rr],fontsize=10)
    ax.invert_yaxis();ax.set_xlim(-.008,.43);ax.set_xticks([0,.1,.2,.3,.4])
    ax.set_title(c,fontsize=14,pad=16,fontweight='bold');ax.grid(axis='x',alpha=.18)
    for s in ['top','right','left']:ax.spines[s].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.set_xlabel('供者等权平均表达',fontsize=10)
    m=next(r for r in rr if r['celltype']=='Myeloid'); rank=rr.index(m)+1
    text=(f"髓系：第 {rank} / {len(rr)} 名\n合格来源标签 {m['n_source_labels_eligible']}；细胞 {int(m['n_cells_eligible']):,}\n平均检出比例 {float(m['mean_detection_fraction']):.1%}\n重抽样首位频率 {float(m['bootstrap_first_frequency']):.1%}")
    ax.text(0,-.18,text,transform=ax.transAxes,va='top',fontsize=10,linespacing=1.5,color='#8d3329')
    summ.append({'cohort':c,'rank':rank,'n_celltypes':len(rr),'myeloid_bootstrap_first':m['bootstrap_first_frequency']})
fig.text(.04,.045,'点：供者等权均值；横线：供者均值的四分位范围（不是置信区间）。表达尺度：log1p(每万计数)。\n灰色为其他可评估细胞类型；不可评估类别保留在图源表。各研究分别解释，不作跨研究绝对表达差异检验。',fontsize=10,color='#475569',linespacing=1.6)
for ext in ['png','svg']:fig.savefig(out/('SLC6A6_myeloid.'+ext),dpi=190,facecolor='white')
plt.close(fig)
spec={'gene':'SLC6A6','scope':'existing donor-equal source aggregates only','new_statistical_tests':0,'source_run':'20260922T131800Z_author_identity_v4','input_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'rank_by':'effect descending within each cohort','NA_policy':'retain in TSV; omit from quantitative plot','units':'author source labels; independent clinical donor identity not recertified','matplotlib_version':matplotlib.__version__}
(out/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(out/'validation.json').write_text(json.dumps({'status':'PASS','source_rows':len(rows),'myeloid_rows':3,'unique_cohort_celltype':True,'new_tests':0,'summary':summ,'not_verified':['clinical donor identity','myeloid subtype','disease differential expression']},indent=2)+'\n',encoding='utf-8')
(out/'README_CN.md').write_text('''# PDAC SLC6A6 髓系表达来源展示

## 本轮问题
查看三队列髓系表达，并与其他细胞类型对照。
## 输入与范围
复用 v4 作者身份修正后的公开供者等权汇总；不读取患者或逐细胞矩阵，不新增检验。
## 实际结果
GSE263733：髓系第1/7，均值0.32595，检出比例41.7%，20个合格作者来源标签、9363细胞，重抽样首位81.7%。
GSE278688：髓系第1/6，均值0.30859，检出比例28.5%，13个合格作者来源标签、5654细胞，重抽样首位79.1%。
GSE242230：髓系第2/8，均值0.20646，检出比例21.1%，20个合格作者来源标签、6224细胞，重抽样首位36.0%。内皮均值第1，但其首位频率只有28.8%；CAF首位频率31.9%。
## 新手解释
点为供者等权表达均值，横线为供者均值的四分位范围，非置信区间。检出比例是供者内表达非零细胞比例的等权均值，不是直接合并全部细胞计算。首位频率来自已有1000次重抽样，不是表达细胞百分比。均值排名与重抽样最常首位类别可不同。
## 限制/反证
来源标签的临床供者身份未重新认证。两队列以髓系为首位，但第三套排序不稳定；没有髓系内部亚群结果，不能指定巨噬细胞或M1/M2。未做髓系肿瘤对正常差异检验，无疾病差异P/q，不证明髓系特异或转运机制。
## 当前决定
DONE仅指本次展示；保留两队列髓系偏高的描述性线索，不提升为三队列一致来源。
## 下一步
若进一步定位髓系亚群，先核实作者原有亚群注释与可连接性。
## 复现命令
python code/pdac/show_slc6a6_myeloid.py
每次写入新运行目录，原统计不覆盖。
''',encoding='utf-8')
index=ROOT/'coordination/stages/PDAC.tsv'
with index.open(encoding='utf-8-sig') as f:
    reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames; existing=[r for r in reader if r['run_id']!=run]
existing.append(dict(zip(fields,['PDAC','06_EXTERNAL',run,'SLC6A6_source_display_v1','DONE','Existing three-cohort myeloid source display;no new tests',out.relative_to(ROOT).as_posix(),'code/pdac/show_slc6a6_myeloid.py','analysis/pdac-initial','Two cohorts myeloid top;third unstable;no subtype inference','Inspect author myeloid annotations only if requested'])))
with index.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sorted(existing,key=lambda r:(r['stage_id'],r['run_id'])))
print(out)
