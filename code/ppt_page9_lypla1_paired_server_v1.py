"""Server165 only: render LYPLA1 paired RNA; export plots and aggregates only."""
from pathlib import Path
import json,hashlib,platform
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--page',type=int,default=9);args=parser.parse_args()
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
C=B/'results/collaborative'
O=Path(__file__).resolve().parent;P=O/'public';P.mkdir(exist_ok=False)
lock=O/'.running';lock.open('x').write('page9 paired display only')
R=B/'data/candidates/camp_primary_tissue_multicancer/gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'
sources=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def track(p):sources.append(dict(path=str(p),sha256=sha(p),bytes=Path(p).stat().st_size))
def read(p):track(p);return pd.read_csv(p,sep='\t',dtype=str)
def values(fn,t,n):
 p=R/fn;track(p);d=pd.read_csv(p,index_col=0);d.columns=d.columns.astype(str)
 assert d.index.is_unique and d.columns.is_unique and len(t)==len(n)
 assert len(set(t))==len(t) and len(set(n))==len(n) and not set(t)&set(n)
 a=d.loc['LYPLA1',list(t)].to_numpy(float);b=d.loc['LYPLA1',list(n)].to_numpy(float)
 assert np.isfinite(a).all() and np.isfinite(b).all()
 return b,a
m=read(C/'BRCA/A/20260921T102429Z_camp_sample_identity_v1/private/audited_mapping.tsv')
bad=m.TN.ne(m.pdf_TN);assert bad.sum()==1 and m.loc[bad,'GSM'].iloc[0]=='GSM927051'
t=m[m.TN.eq('Tumor')&~bad];n=m[m.TN.eq('Normal')&~bad]
q=t[['case_row','RNAID']].merge(n[['case_row','RNAID']],on='case_row',suffixes=('_t','_n'),validate='one_to_one')
data={'BRCA':values(m.RNAFile.iloc[0],q.RNAID_t,q.RNAID_n)}
m=read(C/'COAD/B/20260919T115842Z_identity_units_v1/private_clinical_join.tsv')
t=m[m.TN.eq('Tumor')&m.stage.isin(['stage I','stage II','stage III','stage IV'])];n=m[m.TN.eq('Normal')].set_index('individual').loc[t.individual]
assert t.individual.is_unique and list(t.individual)==list(n.index)
data['COAD']=values('GSE89076.Agilent8x60K.log2_transformed.gene_symbol.csv',t.RNAID,n.RNAID)
q=read(C/'PDAC/B/20260922T120300Z_sop_v3_internal/private/RNA_pairs_private.tsv');assert q.pair_id.is_unique
data['PDAC']=values('GSE62452.hugene10st.gene_symbol.csv',q.tumor_RNA_id,q.normal_RNA_id)
q=read(C/'PRAD/B/20260922T142000Z_discovery_v1/private/RNA_pairs_private.tsv');assert q['case'].is_unique
m=read(C/'PRAD/B/20260922T142000Z_discovery_v1/private/tumor_multiomics_map_private.tsv')
data['PRAD']=values(m.RNAFile.iloc[0],q.RNAID_tumor,q.RNAID_normal)
expected=pd.read_csv(O/'expected.tsv',sep='\t').set_index('cancer');track(O/'expected.tsv')
summary=[]
for c,(n,t) in data.items():
 r=expected.loc[c];delta=t-n
 assert len(n)==int(r.n) and np.isclose(delta.mean(),r.effect,rtol=0,atol=1e-12),(c,len(n),delta.mean())
 summary.append(dict(cancer=c,gene='LYPLA1',n_pairs=len(n),mean_tumor_minus_normal=float(delta.mean()),p_value=float(r.p),up_pairs=int((delta>0).sum()),down_pairs=int((delta<0).sum()),equal_pairs=int((delta==0).sum()),up_fraction=float((delta>0).mean()),expression_scale='Original author processed expression; cohort-specific axis',statistics='Frozen paired t P reused; no new inferential tests'))
summary=pd.DataFrame(summary);summary.to_csv(P/'paired_summary.tsv',sep='\t',index=False)
for fn in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont(str(O/fn))
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0'
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=15,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(O/'sysu_logo.png'));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'LYPLA1 在四癌多数配对患者中上调' if args.page==10 else '为什么聚焦 LYPLA1？',29 if args.page==10 else 31,G,True)
text(.052,.755,'来自 CAMP 代谢物直接关系候选池  ·  四癌 RNA 平均变化均为上调，其中三癌显著',13.5,M)
for i,(c,(n,t)) in enumerate(data.items()):
 x=.073+i*.231;r=summary[summary.cancer==c].iloc[0]
 text(x+.078,.678,c,19,G,True,ha='center')
 text(x+.078,.638,f'{len(n)} 组配对  |  P = {r.p_value:.3g}',11.5,INK,ha='center')
 ax=f.add_axes([x,.322,.168,.266]);ax.set_xlim(-.35,1.35)
 for j in range(len(n)):
  color=UP if t[j]>n[j] else DOWN if t[j]<n[j] else M
  ax.plot([0,1],[n[j],t[j]],color=color,alpha=.33,lw=.85,zorder=1)
 ax.scatter(np.zeros(len(n)),n,s=13,color='#879B90',alpha=.7,linewidths=0,zorder=2)
 ax.scatter(np.ones(len(n)),t,s=13,color=G,alpha=.65,linewidths=0,zorder=2)
 for j,a in enumerate([n,t]):ax.plot([j-.14,j+.14],[a.mean(),a.mean()],color=INK,lw=2.6,zorder=4)
 lo=min(n.min(),t.min());hi=max(n.max(),t.max());pad=(hi-lo)*.13;ax.set_ylim(lo-pad,hi+pad)
 ax.set_xticks([0,1],['正常','肿瘤'],fontsize=12,color=INK);ax.tick_params(axis='y',labelsize=10,color=M,labelcolor=M,length=3);ax.tick_params(axis='x',length=0,pad=9)
 ax.spines[['top','right']].set_visible(False);ax.spines['left'].set_color(LINE);ax.spines['bottom'].set_color(LINE)
 ax.set_ylabel('LYPLA1 表达（作者处理尺度）',fontsize=10.5,color=M,labelpad=8)
 text(x+.078,.255,f'{int(r.up_pairs)} / {int(r.n_pairs)} 位患者上调',13,G,True,ha='center')
 text(x+.078,.215,f'平均差值  +{r.mean_tumor_minus_normal:.3f}'+('  ·  未显著' if r.p_value>=.05 else ''),11,M,ha='center')
 # Keep the four panels separated without boxing their data.
 if i<3:f.add_artist(Line2D([x+.207,x+.207],[.218,.703],transform=f.transFigure,color='#EDF1ED',lw=.8))
text(.052,.163,'连线：同一患者的正常与肿瘤组织；暖色为升高、蓝色为降低；深色短横线为组均值。',10.2,M)
line(.05,.95,.137)
text(.05,.105,'关注理由',13,G,True)
text(.15,.105,'代谢关系候选 + 四癌同向表达背景',13,INK)
text(.61,.105,'下一问：上调来自恶性上皮吗？',13,G,True)
text(.05,.060,'PDAC P = 0.0777；不称“四癌均显著”。这是同一 CAMP 队列的表达支持，不是独立验证或机制证明。',10,M)
text(.05,.031,'各面板保留队列原始处理尺度，不能跨癌比较绝对表达量；配对依据来自既有作者身份审计。',9.4,M)
text(.95,.041,f'{args.page:02}',14,G,True,ha='right')
f.canvas.draw();rend=f.canvas.get_renderer()
for a in f.texts:
 bb=a.get_window_extent(rend);assert bb.x0>=0 and bb.y0>=0 and bb.x1<=f.bbox.width and bb.y1<=f.bbox.height,a.get_text()
for ext in ['png','pdf']:f.savefig(P/f'{args.page:02}_LYPLA1_四癌患者配对.{ext}',dpi=200)
plt.close(f)
pd.DataFrame(sources).to_csv(P/'source_manifest_server.tsv',sep='\t',index=False)
spec=dict(scope='Four cancer LYPLA1 RNA presentation only',source_data_residency='All matrices and sample maps remain server165; export figure and aggregate summaries only',pairing='Inherited audited exact patient/case joins; no inferred pairing',transformation='None',P='Frozen paired t nominal P reused; no new tests; no FDR threshold',display='Separate original-scale axes; thin paired lines; group means',validation='Exact n and mean difference matched to frozen results at atol 1e-12',software=dict(python=platform.python_version(),pandas=pd.__version__,matplotlib=matplotlib.__version__),script_sha256=sha(__file__))
(P/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2))
(P/'validation.json').write_text(json.dumps(dict(status='DONE',source_hashes_unchanged=all(sha(x['path'])==x['sha256'] for x in sources),pair_counts=summary.n_pairs.tolist(),exact_frozen_means_matched=True,patient_values_exported=False,new_inferential_tests=0,visual_review='PENDING'),indent=2))
lock.unlink();print(summary.to_string(index=False))
