"""Two image-only result slides from pinned public aggregates; no new patient tests."""
from pathlib import Path
import hashlib, json, subprocess, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle, FancyBboxPatch

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/BRCA/07_INTEGRATION/20260928T150000Z_ppt_redraw_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T120000Z_pages67_v1'
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'sources').mkdir(exist_ok=True)
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_pages_06_07_v1')
VIEW.mkdir(parents=True,exist_ok=True)
order=['BRCA','COAD','GBM','PDAC','PRAD','ccRCC3','ccRCC4']
manifest=pd.read_csv(OLD/'source_manifest.tsv',sep='\t')
used=[]
def get(name):
 p=OLD/'sources'/f'{name}.tsv';r=manifest[manifest.id==name].iloc[0].to_dict()
 assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],name
 used.append(r);return pd.read_csv(p,sep='\t')
tables=[]
for c in ['brca','coad','gbm','pdac','prad','ccrcc']:
 d=get(c+'_met');d['display']=d.cohort if c=='ccrcc' else c.upper();tables.append(d)
allmet=pd.concat(tables,ignore_index=True)
maps={}
for c in ['brca','coad','pdac','prad']:
 d=get(c+'_rel')
 if c=='pdac':d=d[d.mapping_status=='DIRECT']
 maps[c.upper()]=d
extra={
 'GBM':('2174aa965e29cc0fbb10925bc2c217198960829e','results/GBM/02_MAPPING/20260925T103000Z_discovery_v1/direct_relations.tsv'),
 'ccRCC3':('02593f745fc966c975ab1a25e3e7ce74f770efc6','results/ccRCC/02_MAPPING/20260925T151000Z_ccrcc3_v1/direct_relations.tsv'),
 'ccRCC4':('02593f745fc966c975ab1a25e3e7ce74f770efc6','results/ccRCC/02_MAPPING/20260925T155000Z_ccrcc4_histology_v2/direct_relations.tsv')}
for c,(commit,path) in extra.items():
 b=subprocess.check_output(['git','show',commit+':'+path],cwd=ROOT);p=OUT/'sources'/f'{c}_mapping.tsv';p.write_bytes(b);maps[c]=pd.read_csv(p,sep='\t')
 used.append(dict(id=c+'_mapping',git_commit=commit,git_path=path,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),url=f'https://github.com/Asukasssss/-/blob/{commit}/{path}'))
coverage=[]
for c in order:
 z=allmet[allmet.display==c];ok=z.effect.notna()&z.p_value.notna();sig=ok&(z.p_value<.05);m=maps[c]
 feature='feature_id' if 'feature_id' in m else 'metabolite_key'
 assert not m.duplicated([feature,'gene']).any(),c
 assert set(m[feature]).issubset(set(z.loc[sig,feature])),(c,'mapping_outside_nominal_pool')
 coverage.append(dict(cohort=c,total_features=len(z),evaluable=int(ok.sum()),up=int((sig&(z.effect>0)).sum()),down=int((sig&(z.effect<0)).sum()),ns=int((ok&~sig).sum()),not_evaluable=int((~ok).sum()),mapped_features=m[feature].nunique(),direct_relations=len(m),candidate_genes=m.gene.nunique()))
cov=pd.DataFrame(coverage);cov['significant']=cov.up+cov.down
assert (cov.up+cov.down+cov.ns==cov.evaluable).all()
assert cov.evaluable.tolist()==[318,159,697,307,361,711,904]
cov.to_csv(OUT/'06_overview_data.tsv',sep='\t',index=False)

# Identity-exact display, no alias guesses, no pooling duplicate peaks.
keys=sorted(k for k in allmet.metabolite_key.dropna().unique() if str(k).startswith(('KEGG:','HMDB:')))
cells=[];ranks=[]
for key in keys:
 z=allmet[allmet.metabolite_key==key];observed=[]
 for c in order:
  q=z[z.display==c];status='ABSENT_IN_SOURCE' if len(q)==0 else 'AMBIGUOUS_MULTIPLE_FEATURES' if len(q)>1 else 'NOT_EVALUABLE' if pd.isna(q.iloc[0].p_value) or pd.isna(q.iloc[0].effect) else 'EVALUABLE'
  r=q.iloc[0] if len(q)==1 else None
  row=dict(key=key,cohort=c,status=status,n_source_rows=len(q),source_names=' | '.join(sorted(q.metabolite_name.astype(str).unique())),effect=r.effect if r is not None else np.nan,p=r.p_value if r is not None else np.nan,effect_type=r.effect_type if r is not None else '',n=r['n'] if r is not None else np.nan)
  cells.append(row)
  if status=='EVALUABLE':observed.append(row)
 nc=len(set('ccRCC' if x['cohort'].startswith('ccRCC') else x['cohort'] for x in observed))
 sig=[x for x in observed if x['p']<.05];nsc=len(set('ccRCC' if x['cohort'].startswith('ccRCC') else x['cohort'] for x in sig))
 same=len(set(np.sign(x['effect']) for x in observed))==1
 opposite=any(np.sign(a['effect'])!=np.sign(b['effect']) and a['cohort'].replace('ccRCC3','ccRCC').replace('ccRCC4','ccRCC')!=b['cohort'].replace('ccRCC3','ccRCC').replace('ccRCC4','ccRCC') for a in sig for b in sig)
 ranks.append(dict(key=key,n_cancers=nc,n_cohorts=len(observed),n_significant_cancers=nsc,n_significant_cohorts=len(sig),pattern='same' if same else 'opposite' if opposite else 'mixed',source_names=' | '.join(sorted(z.metabolite_name.astype(str).unique()))))
rank=pd.DataFrame(ranks).sort_values(['n_cancers','n_cohorts','n_significant_cancers','key'],ascending=[False,False,False,True])
# Restrict representative panel to KEGG concrete molecule IDs, excluding generic lipid classes.
eligible=rank[rank.key.str.startswith('KEGG:')&~rank.key.isin(['KEGG:C00157'])&(rank.n_cancers>=4)&(rank.n_significant_cancers>=2)]
same=eligible[eligible.pattern=='same'].head(4);opp=eligible[eligible.pattern=='opposite'].head(6)
selected=pd.concat([same,opp]);assert len(selected)==10
matrix=pd.DataFrame(cells);selected_cells=matrix[matrix.key.isin(selected.key)]
rank.to_csv(OUT/'07_all_identity_ranking.tsv',sep='\t',index=False);matrix.to_csv(OUT/'07_all_identity_matrix.tsv',sep='\t',index=False)
selected.to_csv(OUT/'07_selected_molecules.tsv',sep='\t',index=False);selected_cells.to_csv(OUT/'07_selected_cells.tsv',sep='\t',index=False)
print(cov.to_string(index=False));print(selected.to_string(index=False))

font_manager.fontManager.addfont(str(ROOT/'runtime/ppt_fonts/NotoSansCJK-Regular.ttc'))
for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
GREEN='#00563C';INK='#173D32';MUTED='#697B73';PALE='#EDF4EE';LINE='#D9E6DE';UP='#C36C62';DOWN='#4A8499';GRAY='#E7EBE9'
def txt(fig,x,y,s,size=16,color=INK,weight='normal',ha='left',va='center',**kw):
 return fig.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va=va,**kw)
def rect(fig,x,y,w,h,color,edge='none',radius=.012):
 fig.add_artist(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={radius}',transform=fig.transFigure,fc=color,ec=edge,lw=.8,zorder=-2))
def base(num,title,sub):
 fig=plt.figure(figsize=(16,9),facecolor='white')
 txt(fig,.047,.947,'中 山 大 学',17,GREEN,'bold');txt(fig,.047,.917,'SUN YAT-SEN UNIVERSITY',7.5,GREEN)
 txt(fig,.949,.943,f'{num:02}',17,GREEN,ha='right');fig.add_artist(plt.Line2D([.926,.95],[.918,.918],transform=fig.transFigure,color=GREEN,lw=1.2))
 txt(fig,.047,.846,title,29,GREEN,'bold');txt(fig,.047,.788,sub,13.5,MUTED)
 return fig
def finish(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=200,facecolor='white')
 plt.close(fig)
 for ext in ['png','pdf','svg']:shutil.copy2(OUT/f'{name}.{ext}',VIEW/f'{name}.{ext}')

# Page 6: aligned bars and a compact candidate coverage table.
fig=base(6,'六癌均发现代谢异常，候选覆盖规模不同','7 个队列分别展示：先看检测范围，再看差异代谢物能够连接到多少候选基因')
rect(fig,.043,.193,.914,.54,'#F7FAF7')
txt(fig,.063,.693,'肿瘤相对参照的代谢变化',16,GREEN,'bold')
for x,col,label in [(.063,UP,'升高'),(.145,DOWN,'降低'),(.227,GRAY,'P ≥ 0.05')]:
 rect(fig,x,.65,.012,.012,col,radius=.002);txt(fig,x+.017,.656,label,10.5,MUTED)
txt(fig,.70,.693,'直接关系映射的覆盖',16,GREEN,'bold')
xs=[.638,.749,.825,.909]
for x,s in zip(xs,['差异 / 可评估','已映射\n代谢特征','候选\n基因','直接\n关系']):txt(fig,x,.648,s,10.5,MUTED,ha='center',linespacing=1.4)
ys=np.linspace(.590,.275,7)
barx=.195;barw=.36
design=['45 对','33 对','74 肿瘤 / 6 参照','11 对','43 对','17 对患者均值','12 对患者均值']
for i,r in enumerate(cov.itertuples()):
 y=ys[i]
 if i%2==0:rect(fig,.055,y-.022,.892,.045,'white',radius=.004)
 txt(fig,.068,y+.007,r.cohort,14,GREEN,'bold');txt(fig,.068,y-.012,design[i],8.5,MUTED)
 left=barx
 for n,col in [(r.up,UP),(r.down,DOWN),(r.ns,GRAY)]:
  width=barw*n/950;rect(fig,left,y-.012,width,.024,col,radius=.002)
  if width>.027:txt(fig,left+width/2,y,str(n),8.8,'white' if col!=GRAY else MUTED,ha='center')
  left+=width
 for x,val in zip(xs,[f'{r.significant} / {r.evaluable}',str(r.mapped_features),str(r.candidate_genes),str(r.direct_relations)]):txt(fig,x,y,val,13,INK,ha='center',weight='bold' if x==xs[0] else 'normal')
for t in [0,300,600,900]:
 x=barx+barw*t/950;txt(fig,x,.232,str(t),9,MUTED,ha='center')
txt(fig,(barx+barx+barw)/2,.209,'代谢特征数',9,MUTED,ha='center')
rect(fig,.047,.104,.906,.058,PALE)
txt(fig,.069,.133,'下一步：比较具体分子，寻找跨癌重复与不同方向的代谢变化',16,GREEN,'bold')
txt(fig,.048,.067,'显著性：名义 P < 0.05。映射库覆盖不同、非穷尽；候选数量不用于比较癌种异常强弱。',9,MUTED)
txt(fig,.048,.04,'GBM 为非配对比较，参照组织背景存在混杂；ccRCC 多区域先按患者汇总。原效应与 P 值均复用。',9,MUTED)
finish(fig,'06_六癌代谢异常与候选覆盖')

# Page 7: direction categories rather than a falsely shared effect-size scale.
labels={'KEGG:C00328':'犬尿氨酸','KEGG:C00137':'肌醇','KEGG:C00319':'鞘氨醇','KEGG:C02990':'棕榈酰肉碱','KEGG:C00031':'葡萄糖','KEGG:C00864':'泛酸（维生素 B5）','KEGG:C00025':'谷氨酸','KEGG:C00065':'丝氨酸','KEGG:C00123':'亮氨酸','KEGG:C00148':'脯氨酸','KEGG:C00158':'柠檬酸','KEGG:C00519':'亚牛磺酸'}
labels['KEGG:C00047']='赖氨酸'
assert all(k in labels for k in selected.key),selected.key.tolist()
fig=base(7,'跨癌代谢变化：同向重复与相反方向并存','同一化学标识逐项对照；ccRCC 两队列分列，缺测与身份歧义单独标记')
rect(fig,.043,.185,.705,.544,'#F7FAF7');rect(fig,.775,.185,.18,.544,PALE)
x0=.296;dx=.059;ys=np.r_[np.linspace(.643,.517,4),np.linspace(.439,.229,6)]
for j,c in enumerate(order):txt(fig,x0+j*dx,.695,c,10.8,GREEN,'bold',ha='center')
txt(fig,.06,.648,'同\n向\n趋\n势',10,GREEN,'bold',va='top',linespacing=1.5)
txt(fig,.06,.444,'相\n反\n方\n向',10,GREEN,'bold',va='top',linespacing=1.5)
fig.add_artist(plt.Line2D([.104,.725],[.478,.478],transform=fig.transFigure,color=LINE,lw=1.2))
for i,r in enumerate(selected.itertuples()):
 y=ys[i];txt(fig,.104,y,labels[r.key],12,INK)
 for j,c in enumerate(order):
  a=selected_cells[(selected_cells.key==r.key)&(selected_cells.cohort==c)].iloc[0];x=x0+j*dx
  status=a.status
  if status=='EVALUABLE':
   pos=a.effect>0;col=UP if pos else DOWN;significant=a.p<.05;fill=col if significant else ('#F4E5E1' if pos else '#E4EEF1');label='↑' if pos else '↓';tc='white' if significant else col
  else:fill='#E6EAE7';label={'ABSENT_IN_SOURCE':'—','AMBIGUOUS_MULTIPLE_FEATURES':'?','NOT_EVALUABLE':'×'}[status];tc=MUTED
  rect(fig,x-.024,y-.016,.048,.032,fill,radius=.004);txt(fig,x,y,label,13,tc,ha='center',weight='bold')
  if status=='EVALUABLE' and a.p<.05:txt(fig,x+.016,y+.006,'*',9,'white',ha='center')
txt(fig,.796,.683,'如何读这张图',15,GREEN,'bold')
txt(fig,.796,.624,'犬尿氨酸',14,GREEN,'bold');txt(fig,.796,.565,'可明确对应的队列\n均为升高方向；\n部分队列未显著。',11.5,INK,linespacing=1.7)
txt(fig,.796,.462,'葡萄糖 / 谷氨酸',13,GREEN,'bold');txt(fig,.796,.412,'不同癌种中出现\n显著的相反方向。',11.5,INK,linespacing=1.7)
txt(fig,.796,.315,'解释边界',13,GREEN,'bold');txt(fig,.796,.253,'方向对照提供线索，\n不等于已证实\n癌种特异机制。',10.5,MUTED,linespacing=1.6)
# Compact legend below matrix.
txt(fig,.051,.155,'红 ↑ 升高    蓝 ↓ 降低    * P < 0.05    浅色无 *：P ≥ 0.05    — 未覆盖    ? 多特征同标识    × 不可评估',10,MUTED)
rect(fig,.047,.079,.906,.05,PALE);txt(fig,.069,.104,'跨癌共性需要逐个分子核对；下一步再连接到有生化依据的候选基因',15,GREEN,'bold')
txt(fig,.048,.047,'选图：至少覆盖 4 癌种；按癌种覆盖、队列覆盖、显著癌种数排序，分别取 4 个同向、6 个反向分子。',8.8,MUTED)
txt(fig,.048,.023,'限具体 KEGG 分子；同向组均至少 2 癌种显著。GBM 非配对，其他按配对设计；色深不代表跨队列效应幅度。',8.8,MUTED)
finish(fig,'07_跨癌代谢物方向对照')

pd.DataFrame(used).drop_duplicates('id').to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
spec=dict(version='pages67_v1',nominal_p=.05,new_patient_tests=False,original_p_q_effects_unchanged=True,source_snapshot='20260928 presentation aggregate manifest; latest fetched COAD/PDAC commits have no changes in the relevant stages',effect_display='direction and nominal significance; not comparable effect magnitudes',selection='Concrete KEGG IDs, unique feature per cohort, >=4 cancers and >=2 significant cancers; rank by cancer coverage, cohort coverage, significant cancer count, key; top4 all-observed same direction, top6 opposite significant directions in different cancer types',missing_statuses=['ABSENT_IN_SOURCE','AMBIGUOUS_MULTIPLE_FEATURES','NOT_EVALUABLE'],output='two PNG/PDF/SVG slides; no PPT')
(OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'validation.json').write_text(json.dumps(dict(source_hashes_verified=True,coverage_counts_checked=True,relation_feature_gene_unique=True,selected_count=len(selected),selected_cells=len(selected_cells),selection_keys=selected.key.tolist(),no_new_inference=True),ensure_ascii=False,indent=2),encoding='utf8')
notes='''# 第6–7页：真实结果重绘

## 输入与范围
复用2026-09-28展示包的六癌七队列公开汇总，补充固定提交中的GBM及ccRCC映射关系。源表哈希全部核对；不读取患者矩阵、不重算P/q，不生成PPT。

## 第6页讲法
七个队列均得到名义显著代谢变化。先展示可评估范围、升降方向，再展示直接映射的特征、候选基因与关系数。候选是生化关系覆盖，不是已验证靶点。不同数据库补充范围、检测范围及设计决定数量不可当作癌种强弱排名。BRCA190差异→93已映射特征→150基因；其余队列见06_overview_data.tsv。GBM条件性旧标签说明以当前整合视图为准，仍保留非配对参照混杂。

## 第7页讲法
预先固定排序规则，从完整化学标识矩阵选择同向与反向代表，不按本页故事挑选最小P。犬尿氨酸在可唯一对应的6队列均为升高方向，其中4队列、3癌种P<0.05；GBM当前源表没有对应标识，单独标为未覆盖。葡萄糖、谷氨酸展示不同癌种中显著相反方向。由于效应尺度不同，采用方向/显著性分类而非连续幅度色标。

## 边界
同向指均值效应方向，不表示每位患者同向；显著性为名义P。单癌显著不等于癌特异；未做交互检验。ccRCC两队列在癌种计数中只算一个癌种。未覆盖指当前源表没有对应标识，不能断言平台从未测量。问号代表同队列同标识存在多特征；不按名称猜测合并。其他表中已有HMDB标识仍保留在完整矩阵，本次代表面板限具体KEGG条目。

## 复现
python code/ppt_pages67_v1.py
完整选图排序、长表、来源提交与哈希均随附。所有原始统计值保留。
'''
(OUT/'README_CN.md').write_text(notes,encoding='utf8')
for name in ['README_CN.md','06_overview_data.tsv','07_selected_molecules.tsv','07_selected_cells.tsv','analysis_spec.json','source_manifest.tsv','validation.json']:
 shutil.copy2(OUT/name,VIEW/name)
html='<!doctype html><html lang="zh"><meta charset="utf-8"><title>CAMP · 第6–7页</title><style>body{margin:0;background:#eaf0ec;font-family:system-ui}main{max-width:1440px;margin:24px auto}img{display:block;width:100%;box-shadow:0 8px 30px #183e3320;margin:24px 0}a{color:#00563c}</style><main><h2>CAMP · 组会汇报第6–7页</h2><p>真实汇总数据代码重绘 · 名义 P&lt;0.05 · 原统计不变 · 未生成 PPT</p>'
for name in ['06_六癌代谢异常与候选覆盖','07_跨癌代谢物方向对照']:html+=f'<img src="{name}.png"><a href="{name}.pdf">矢量 PDF</a>'
html+='<p><a href="README_CN.md">讲解与方法</a></p></main>'
(VIEW/'index.html').write_text(html,encoding='utf8')
print('SAVED',VIEW)
