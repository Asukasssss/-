"""Concrete metabolite discovery and biochemical mapping explain the LYPLA1 focus."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/BRCA/07_INTEGRATION/20260928T150000Z_ppt_redraw_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T220000Z_metabolite_bridge_v2'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_LYPLA1_metabolite_bridge_v2')
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
manifest=pd.read_csv(OLD/'source_manifest.tsv',sep='\t');used=[];rows=[]
for c in ['brca','coad','pdac','prad']:
 p=OLD/'sources'/f'{c}_rel.tsv';src=manifest[manifest.id.eq(c+'_rel')].iloc[0].to_dict();assert hashlib.sha256(p.read_bytes()).hexdigest()==src['sha256'];used.append(src)
 d=pd.read_csv(p,sep='\t');d=d[d.gene.eq('LYPLA1')]
 if c=='pdac':d=d[d.mapping_status.eq('DIRECT')]
 me,mp={'brca':('metabolite_paired_effect','metabolite_paired_p_value'),'coad':('metabolite_paired_effect','metabolite_paired_p_value'),'pdac':('metabolite_effect','metabolite_p'),'prad':('metab_effect','metab_p_value')}[c]
 ae,ap,an={'brca':('CAMP_effect','CAMP_p_value','CAMP_n'),'coad':('CAMP_effect','CAMP_p_value','CAMP_n'),'pdac':('association_effect','association_p','association_n'),'prad':('patient_primary_effect','patient_primary_p_value','patient_primary_n')}[c]
 for _,r in d.iterrows():
  if c=='pdac' and ae not in r:ae='association_rho'
  rows.append(dict(cancer=c.upper(),key=r.metabolite_key,metabolite=r.metabolite_name,metabolite_effect=r[me],metabolite_p=r[mp],rho=r[ae],association_p=r[ap],association_n=r[an]))
e=pd.DataFrame(rows);assert len(e)==6 and e.metabolite_p.lt(.05).all() and e.association_p.lt(.05).sum()==3
labels={'KEGG:C04102':'LPC(16:0)','KEGG:C00670':'GPC','KEGG:C00712':'油酸'}
desired=[('BRCA','KEGG:C04102'),('BRCA','KEGG:C00670'),('BRCA','KEGG:C00712'),('COAD','KEGG:C00670'),('PDAC','KEGG:C00670'),('PRAD','KEGG:C00712')]
e=e.set_index(['cancer','key']).loc[desired].reset_index();e.to_csv(OUT/'metabolite_to_LYPLA1_evidence.tsv',sep='\t',index=False)
ref=ROOT/'results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/relation_evidence.tsv';r=pd.read_csv(ref,sep='\t');r=r[r.gene.eq('LYPLA1')];assert len(r)==3
r.to_csv(OUT/'biochemical_mapping_evidence.tsv',sep='\t',index=False);used.append(dict(id='biochemical_evidence',git_path=str(ref.relative_to(ROOT)),sha256=hashlib.sha256(ref.read_bytes()).hexdigest()))
pd.DataFrame(used).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0';PALE='#F2F7F3'
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=15,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(plt.Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
def box(x,y,w,h,fc=PALE):f.add_artist(FancyBboxPatch((x,y),w,h,transform=f.transFigure,boxstyle='round,pad=0,rounding_size=.009',fc=fc,ec='none',zorder=-3))
logo=ROOT/'results/BRCA/07_INTEGRATION/20260929T200000Z_page8_patient_consistency_v3/sysu_logo.png'
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(logo));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'从代谢物异常到 LYPLA1：为什么关注它？',29,G,True)
text(.052,.755,'四癌的差异代谢物关系反复涉及 LYPLA1，再检查患者关联与基因表达',13.5,M)
text(.052,.687,'01  先发现具体代谢物异常',16,G,True)
text(.052,.653,'肿瘤—正常配对比较；以下 6 项均 P < 0.05',10.5,M)
for x,s in [(.088,'癌种'),(.214,'代谢物'),(.335,'方向'),(.418,'P 值')]:text(x,.607,s,10.5,M,ha='center')
line(.05,.463,.583)
for i,a in enumerate(e.itertuples()):
 y=.550-i*.047
 if i in [0,1,2]:box(.05,y-.022,.413,.044,'#F3F7F3')
 elif i%2==0:box(.05,y-.022,.413,.044,'#F7F9F7')
 if i==1:text(.088,y,'BRCA',13,G,True,ha='center')
 elif i>=3:text(.088,y,a.cancer,13,G,True,ha='center')
 text(.214,y,labels[a.key],13,INK,True,ha='center')
 text(.335,y,'↑' if a.metabolite_effect>0 else '↓',17,UP if a.metabolite_effect>0 else DOWN,True,ha='center')
 text(.418,y,f'{a.metabolite_p:.3g}',11.5,INK,ha='center')
text(.05,.249,'四癌均有代谢物线索，映射到同一个候选。',12.4,G,True)
text(.05,.225,'GPC：甘油磷酸胆碱；油酸：18:1n9。',9.4,M)
# Mapping connector is annotation, not a patient-level causal claim.
f.add_artist(FancyArrowPatch((.476,.510),(.510,.510),transform=f.transFigure,arrowstyle='-|>',mutation_scale=15,color='#9AB2A3',lw=1.6))
text(.528,.687,'02  依据生化反应映射到 LYPLA1',16,G,True)
box(.522,.477,.429,.174)
text(.546,.582,'LPC(16:0)',14,G,True)
text(.691,.620,'LYPLA1',14,G,True,ha='center')
text(.691,.581,'→',24,G,ha='center')
text(.785,.582,'GPC + 棕榈酸',14,G,True)
text(.546,.553,'底物',10,M)
text(.546,.530,'油酸：另一条 LYPLA1 注释反应中的产物',12,INK)
text(.546,.494,'反应参与关系 ≠ 患者内净通量或致因关系',9.7,M)
text(.528,.437,'03  肿瘤样本中的代谢物—RNA 关联',14.5,G,True)
sig=e[e.association_p.lt(.05)]
for i,a in enumerate(sig.itertuples()):
 y=.396-i*.040
 text(.535,y,f'{a.cancer} · {labels[a.key]}  (n={int(a.association_n)})',10.8,INK,True)
 text(.752,y,f'ρ = {a.rho:+.3f}',11.5,INK)
 text(.873,y,f'P = {a.association_p:.4g}',11,INK)
text(.528,.266,'其余 3 条未显著：BRCA / COAD 的 GPC、PRAD 的油酸',9.8,M)
line(.05,.95,.214)
box(.05,.127,.90,.068)
text(.067,.161,'因此继续关注',13,G,True)
text(.196,.161,'代谢关系候选反复出现  +  RNA 四癌平均上调、三癌显著',16,G,True)
text(.05,.094,'下一步：看患者配对表达，再定位到恶性上皮；目前尚未证明同一代谢物—LYPLA1 关联跨癌重复。',10.5,INK)
text(.05,.060,'RNA 不能代替酶活或通量；PDAC RNA 未显著，ccRCC RNA 存在反向。不同脂质的显著关联不合并成同一代谢轴。',9.3,M)
text(.05,.030,'映射来源：UniProt O75608；Rhea 40435 / 41720。反应为简写；油酸不是 LPC(16:0) 水解产生的脂肪酸。',9.0,M)
text(.95,.041,'09',14,G,True,ha='right')
f.canvas.draw();rend=f.canvas.get_renderer()
for t in f.texts:
 b=t.get_window_extent(rend);assert b.x0>=0 and b.y0>=0 and b.x1<=f.bbox.width and b.y1<=f.bbox.height,t.get_text()
name='09_代谢物如何指向LYPLA1'
for ext in ['png','pdf','svg']:
 p=OUT/f'{name}.{ext}';f.savefig(p,dpi=200)
 if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
 shutil.copy2(p,VIEW/p.name)
plt.close(f)
spec=dict(version='metabolite_bridge_v2',scope='Four-cancer current LYPLA1 biochemical candidate relations',relation_count=6,metabolites=labels,metabolite_discovery='Tumor-normal within-patient differences',patient_association='Matched metabolite and RNA within tumor samples; not paired tumor-normal deltas',significance='nominal P<0.05',source='Pinned existing aggregates and biochemical evidence, no new tests',reaction_note='LPC16:0 hydrolysis yields GPC and palmitate; oleate belongs to distinct Rhea41720 reaction; water/proton omitted in diagram',no_claims=['one shared metabolic axis across four cancers','expression proves enzyme activity','metabolite-RNA association is causality','all four RNA significant'],web_checked=['https://www.rhea-db.org/rhea/40435','https://www.rhea-db.org/rhea/41720'])
(OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'README_CN.md').write_text('''# 第9页修订：从代谢物异常到LYPLA1

本轮问题：前版仅写通用映射流程，没有具体代谢物，不能说明LYPLA1如何出现。新版直接展示六条当前候选关系对应的代谢差异、注释依据与患者关联。

输入与范围：固定来源版本的BRCA3条、COAD1条、PDAC1条、PRAD1条LYPLA1直接关系。全部已有统计复用，不改阈值、不新增检验。6条对应三种分子：LPC(16:0)、GPC、油酸。

实际结果：BRCA LPC(16:0)升高、GPC降低、油酸升高；COAD与PDAC GPC降低；PRAD油酸升高，均名义P<0.05。这些异常分子经已有直接反应注释映射到同一LYPLA1候选。肿瘤样本内关联显著为BRCA LPC(16:0) rho=0.481712 P=0.0001、BRCA油酸rho=0.384218 P=0.0025、PDAC GPC rho=-0.441558 P=0.0455。BRCA GPC、COAD GPC、PRAD油酸关联均不显著，全部在附表与图中文字保留。

新手解释：第一层是代谢物在肿瘤中改变，第二层是已知生化反应提出候选基因，第三层是同标本配对组学检验相关。最后RNA组织表达支持继续关注，但不能替代代谢证据。四癌重复的是候选基因，不是同一个代谢物—基因相关四癌均显著。

映射核对：冻结UniProt O75608证据中，LPC(16:0)是Rhea40435底物，GPC是产物；同一反应脂肪酸产物为棕榈酸，绝不是油酸。油酸来自另一条Rhea41720所述反应。2026-09-29在线核对Rhea反应式与冻存记录一致，UniProt网页本次仅返回JS提示，基因归属使用既有UniProt注释快照及来源证据。简化图省略水与质子，不表征净通量。

限制/反证：RNA上调并不能推出酶活增强或解释底物/产物升降。尚无同一代谢物—LYPLA1关系在这些癌种中重复显著的证据。不同脂质的显著结果不能合并成同一机制。全部为探索性名义P；RNA仍是四癌平均上调三癌显著，PDAC未显著；ccRCC反向在前页保留。

当前决定：以本页取代旧的泛化动机页，后接第10页四癌RNA配对细节，再接单细胞。暂不生成PPT。

复现：python code/ppt_page9_metabolite_bridge_v2.py。完整六条关系见metabolite_to_LYPLA1_evidence.tsv；反应证据见biochemical_mapping_evidence.tsv。
''',encoding='utf8')
(OUT/'validation.json').write_text(json.dumps(dict(status='DONE',source_hashes_checked=True,six_relations=True,all_six_metabolite_p_lt005=True,three_of_six_associations_nominal=True,no_new_tests=True,reaction_oleate_separated_from_palmitate=True,visual_review='PENDING'),indent=2),encoding='utf8')
for n in ['README_CN.md','metabolite_to_LYPLA1_evidence.tsv','biochemical_mapping_evidence.tsv']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>CAMP · 代谢物如何指向LYPLA1</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1440px;margin:24px auto}}img{{width:100%;display:block;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>第9页：从具体代谢物连接到 LYPLA1</h2><img src="{name}.png"><p><a href="{name}.pdf">矢量 PDF</a> · <a href="metabolite_to_LYPLA1_evidence.tsv">六条关系与原统计</a> · <a href="README_CN.md">逻辑说明</a></p></main>',encoding='utf8')
print(e.to_string(index=False));print(VIEW)
