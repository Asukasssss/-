"""Four-cancer evidence summary: frozen aggregate data only."""
from pathlib import Path
import csv,json,shutil,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/BRCA/07_INTEGRATION'
RUN='20260929T_lypla1_evidence_summary_v1';OUT=BASE/RUN;OUT.mkdir(exist_ok=True)
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_LYPLA1_evidence_summary_v1');VIEW.mkdir(exist_ok=True)
paths={'metabolite':BASE/'20260929T220000Z_metabolite_bridge_v2/metabolite_to_LYPLA1_evidence.tsv','rna':BASE/'20260929T210100Z_page10_paired_layout_v2/paired_summary.tsv','sc':BASE/'20260929T071500Z_four_sc_v3/patient_summary.tsv','spatial':BASE/'20260929T_selected_spatial_v1/display_checks.json'}
def read(p):
 with p.open(encoding='utf-8-sig') as h:return list(csv.DictReader(h,delimiter='\t'))
rna={r['cancer']:r for r in read(paths['rna'])};sc={r['cancer']:r for r in read(paths['sc'])};met=read(paths['metabolite']);sp=json.loads(paths['spatial'].read_text(encoding='utf-8'))
cs=['BRCA','COAD','PDAC','PRAD']
assert [int(rna[c]['up_pairs']) for c in cs]==[34,26,8,33]
assert [int(sc[c]['higher']) for c in cs]==[7,9,9,10]
assert len(met)==6 and all(float(r['metabolite_p'])<.05 for r in met)
assert len({r['patient'] for r in sp if r['cancer']=='BRCA'})==2
assert len({r['patient'] for r in sp if r['cancer']=='PRAD'})==1
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont('C:/Windows/Fonts/'+n)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';LINE='#DDE7E0';AMBER='#997238'
f=plt.figure(figsize=(16,9),facecolor='white')
def t(x,y,s,size=12,col=INK,bold=False,ha='left'):
 return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',va='center',ha=ha,linespacing=1.6)
def line(x1,x2,y,col=LINE,w=.8):f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,color=col,lw=w))
def box(x,y,w,h,c):f.add_artist(Rectangle((x,y),w,h,transform=f.transFigure,facecolor=c,edgecolor='none',zorder=-1))
logo=BASE/'20260929T170000Z_page8_genes_six_v2/sysu_logo.png'
ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(logo));ax.axis('off')
t(.95,.93,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.883);line(.05,.103,.883,G,2)
t(.05,.829,'LYPLA1 多层证据汇总',28,G,True)
t(.05,.778,'从代谢异常指向候选，再以组织、细胞与空间证据逐层核对',13,M)
xs=[.32,.50,.68,.86]
t(.06,.710,'证据层次',12,G,True)
for x,c in zip(xs,cs):t(x,.710,c,16,G,True,ha='center')
line(.05,.95,.683,G,1)
for lo,hi in [(.517,.681),(.397,.516),(.259,.396),(.155,.258)]:
 if lo in [.517,.259]:box(.05,lo,.90,hi-lo,'#F3F7F4')
 line(.05,.95,lo)
labels=[(.621,'01  代谢线索','差异分子 → 关系映射'),(.474,'02  组织 RNA','同患者肿瘤—参照'),(.354,'03  单细胞 / 单核','同患者上皮比较'),(.220,'04  空间定位','本次选定代表切片')]
for y,title,sub in labels:t(.06,y,title,13,G,True);t(.06,y-.038,sub,9.3,M)
met_titles=['LPC(16:0) ↑ · 油酸 ↑\nGPC ↓','GPC ↓','GPC ↓','油酸 ↑']
met_notes=['LPC：ρ = 0.48，P = 0.0001\n油酸：ρ = 0.38，P = 0.0025\nGPC 相关不显著','与 LYPLA1 相关不显著\nP = 0.670','与 LYPLA1 负相关\nρ = −0.44，P = 0.0455','与 LYPLA1 相关不显著\nP = 0.668']
for x,title,note in zip(xs,met_titles,met_notes):t(x,.642,title,11.6,G,True,ha='center');t(x,.562,note,8.6,M,ha='center')
for x,c in zip(xs,cs):
 r=rna[c];p=float(r['p_value']);sig=p<.05
 t(x,.480,f"{r['up_pairs']} / {r['n_pairs']} 患者上调",13,G,True,ha='center')
 t(x,.438,('平均上调 · P = '+format(p,'.3g')) if sig else '平均上调 · P = 0.078（未显著）',9.4,G if sig else AMBER,ha='center')
 r=sc[c];t(x,.360,f"{r['higher']} / {r['n_pairs']} 患者更高",13,G,True,ha='center')
 t(x,.319,'P = '+format(float(r['p']),'.3g') if r['p'] else '描述性配对 · 未做显著性检验',9.3,G if r['p'] else AMBER,ha='center')
 t(x,.281,{'BRCA':'恶性 vs 肿瘤内非恶性上皮','COAD':'肿瘤 CNA vs 正常组织上皮','PDAC':'恶性 vs 非恶性导管（全核）','PRAD':'肿瘤内恶性 vs 非恶性上皮'}[c],8.4,M,ha='center')
for x,c in zip(xs,cs):
 if c=='BRCA':t(x,.222,'2 张切片 · 2 位患者',11.6,G,True,ha='center');t(x,.183,'癌区定位展示；无正常上皮参照',8.7,M,ha='center')
 elif c=='PRAD':t(x,.222,'2 张切片 · 1 位患者',11.6,G,True,ha='center');t(x,.183,'癌区 / 良性腺体定位；非患者重复',8.7,M,ha='center')
 else:t(x,.222,'本次未选入展示',11,M,ha='center');t(x,.183,'不据此判断阴性',8.8,M,ha='center')
box(.05,.075,.90,.060,'#E7F0E9')
t(.065,.105,'当前结论',13,G,True)
t(.175,.105,'LYPLA1 是跨癌反复出现、具有上皮表达支持的候选；尚未建立共同代谢机制。',12.2,G,True)
t(.05,.052,'代谢差异均为名义 P < 0.05；四癌重复的是候选基因，不是同一代谢物—LYPLA1 相关均显著。空间图仅作定位补充。',8.5,M)
t(.05,.027,'人数为方向一致性，不代表每位患者均显著。PDAC 仅阳性核的敏感性分析不支持同向升高；不同参照不合并效应。',8.5,M)
t(.95,.027,'证据汇总',9,G,True,ha='right')
stem='LYPLA1_多层证据汇总'
for ext in ['png','pdf','svg']:f.savefig(OUT/(stem+'.'+ext),dpi=180,facecolor='white');shutil.copy2(OUT/(stem+'.'+ext),VIEW/(stem+'.'+ext))
plt.close(f)
for name,p in paths.items():shutil.copy2(p,OUT/('source_'+name+p.suffix))
with (OUT/'source_manifest.tsv').open('w',encoding='utf-8',newline='') as h:
 wr=csv.writer(h,delimiter='\t');wr.writerow(['kind','source','sha256'])
 for name,p in paths.items():wr.writerow([name,str(p.relative_to(ROOT)),hashlib.sha256(p.read_bytes()).hexdigest()])
(OUT/'analysis_spec.json').write_text(json.dumps(dict(version='v1',scope=cs,new_tests=False,statistics='frozen nominal P; COAD/PRAD single-cell descriptive',spatial='selected sections only, not complete cohort comparison',software=matplotlib.__version__),ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'validation.json').write_text(json.dumps(dict(status='PARTIAL',source_assertions='PASS',visual_review='pending',new_tests=False),indent=2),encoding='utf-8')
(OUT/'README_CN.md').write_text('''# LYPLA1 多层证据汇总
## 本轮问题
将已展示四癌的代谢、组织 RNA、单细胞和选定空间切片组织为一页总结。
## 输入与范围
固定版本汇总表，见 source_manifest.tsv；不含 GBM/ccRCC。只读取已有统计，未重算。空间行仅概括本次选定 BRCA/PRAD 切片，不代表完整空间数据盘点。
## 实际结果
四癌均有显著差异代谢物经既有关系映射指向 LYPLA1，但具体分子不同。六条关联中 BRCA LPC/油酸和 PDAC GPC 显著；COAD GPC、PRAD 油酸及 BRCA GPC 不显著。组织 RNA 四癌平均上调，PDAC 未显著。单细胞方向人数为 BRCA 7/8、COAD 9/9、PDAC 9/9、PRAD 10/10；后两种描述性队列具体为 COAD/PRAD，未展示推断性 P。
## 新手解释
列为癌种，行为证据层次。人数表示同患者比较的方向，不表示每一位患者分别显著。空间行展示实际选入的切片与患者数，没有新增区域差异检验。
## 限制/反证
不存在同一代谢物关系四癌均显著证据；bulk 丰度不等于活性或通量。BRCA/PRAD 单细胞参照为肿瘤内非恶性上皮；COAD 为正常组织上皮；PDAC 为非恶性导管。PDAC 仅阳性核敏感性分析不支持同向增加。空间 BRCA 没有正常上皮参照；PRAD 两张来自同一患者，不能作两患者重复。COAD/PDAC 本次未选图不等于未测或阴性。
## 当前决定
定位为有表达支持的跨癌候选，不称共同代谢机制已证实。不生成 PPTX。
## 下一步
组会汇报用本页归纳支持和边界，再转入癌种特异结果。
## 复现命令
python code/ppt_lypla1_evidence_summary_v1.py
''',encoding='utf-8')
html=f'''<!doctype html><meta charset="utf-8"><title>LYPLA1 多层证据汇总</title><style>body{{margin:0;background:#e9eeeb;color:#183d34;font-family:"Microsoft YaHei",sans-serif}}nav{{padding:18px 4%;background:white}}main{{max-width:1680px;margin:24px auto}}img{{width:100%;display:block}}a{{color:#00553b}}</style><nav>LYPLA1 · 多层证据汇总　<a href="{stem}.pdf">PDF</a>　<a href="{stem}.png">高清 PNG</a>　<a href="{stem}.svg">SVG</a></nav><main><img src="{stem}.png"></main>'''
for d in [OUT,VIEW]:(d/'index.html').write_text(html,encoding='utf-8')
print(VIEW)
