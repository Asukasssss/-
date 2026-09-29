"""Render a source-linked literature slide, without generating scientific evidence."""
from pathlib import Path
import csv, json, hashlib, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch

ROOT=Path(__file__).resolve().parents[1]
RUN='20260929T_literature_recent_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION'/RUN
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_LYPLA1_recent_literature_v1')
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont('C:/Windows/Fonts/'+n)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';LINE='#DDE7E0';BG='#F2F7F3'
papers=[
 dict(journal='Journal of Lipid Research',year='2026',date='2026-01; online 2025-11-28',model='乳腺癌 · MDA-MB-231 / MCF7',title='S-palmitoylation of MTDH regulates\nferroptosis resistance in breast cancer cell',doi='10.1016/j.jlr.2025.100953',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC12794514/',level='全文核对',key='连接脂质代谢与铁死亡',nodes=['APT1 / LYPLA1','MTDH 去棕榈酰化','ACSL4 互作 / 铁死亡敏感性'],body='鉴定 APT1 为 MTDH 去棕榈酰化酶；\nMTDH 修饰状态关联脂质组成与铁死亡。',limit='研究中心为 MTDH；不能将全部表型\n归因于 LYPLA1。'),
 dict(journal='Biochemistry and Cell Biology',year='2026',date='2026-06-30',model='乳腺癌 · MDA-MB-468',title='Acyl-protein thioesterase 1 (LYPLA1)\nactivity promotes the growth of MDA-MB-468\ntriple-negative breast cancer cells',doi='10.1139/bcb-2025-0446',url='https://www.sciencedirect.com/org/science/article/pii/S0829821126000158',level='摘要核对',key='直接研究 LYPLA1 的功能',nodes=['抑制 APT1（ML348）','蛋白 S-酰化状态改变','增殖降低 / G1 期积累'],body='整合患者数据、CRISPR 依赖性与体外实验；\n摘要报告 APT1 抑制使细胞生长受限。',limit='目前依据摘要；效应具有模型依赖性，\n表达丰度不等于酶活性。'),
 dict(journal='npj Precision Oncology',year='2025',date='2025-08-28',model='肺癌 · 奥希替尼耐药模型',title='Restoration of Osimertinib sensitivity in lung\ncancer through BRD4 inhibitor-mediated\ndepalmitoylation of mutant EGFR via APT1',doi='10.1038/s41698-025-01048-8',url='https://www.nature.com/articles/s41698-025-01048-8',level='全文核对',key='参与治疗耐药机制',nodes=['APT1 / LYPLA1','EGFR 去棕榈酰化 / 核转位','CMPK2 / 奥希替尼耐药'],body='敲低 APT1 降低耐药细胞 EGFR 核转位，\n并部分恢复对奥希替尼的敏感性。',limit='来自肺癌耐药情境，不能直接外推至\n本项目四癌或未经治疗的肿瘤。')
]
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=13,col=INK,bold=False,ha='left',**kw):
 return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center',linespacing=1.6,**kw)
def line(x1,x2,y,col=LINE,w=.8):f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,color=col,lw=w))
def rect(x,y,w,h,c):f.add_artist(Rectangle((x,y),w,h,transform=f.transFigure,facecolor=c,edgecolor='none',zorder=-1))
logo=ROOT/'results/BRCA/07_INTEGRATION/20260929T170000Z_page8_genes_six_v2/sysu_logo.png'
ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(logo));ax.axis('off')
text(.95,.93,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.883);line(.05,.103,.883,G,2)
text(.05,.826,'LYPLA1：近期研究已提供肿瘤功能与机制线索',26,G,True)
text(.05,.777,'2025–2026 原始研究  ·  APT1 即 LYPLA1  ·  文献支持研究价值，跨癌共性仍需用数据检验',12,M)
xs=[.05,.36,.67];w=.28
for i,(x,p) in enumerate(zip(xs,papers)):
 text(x,.710,p['journal'],12.1,G,True);text(x+w,.711,p['year'],16,G,True,ha='right')
 text(x,.674,p['model'],11,M)
 line(x,x+w,.650,G,1.2)
 text(x,.607,p['title'],9.8,INK)
 text(x,.535,p['key'],17,G,True)
 rect(x,.323,w,.174,BG)
 for j,node in enumerate(p['nodes']):
  yy=.471-j*.059
  text(x+w/2,yy,node,12,G,j==0,ha='center')
  if j<2:f.add_artist(FancyArrowPatch((x+w/2,yy-.018),(x+w/2,yy-.039),transform=f.transFigure,arrowstyle='-|>',mutation_scale=11,color='#83A396',lw=1.2))
 text(x,.277,p['body'],11.1,INK)
 text(x,.209,p['limit'],9.5,M)
 if i==0:text(x,.158,'2026 年卷期 · 2025 年底在线',8.5,M)
 else:text(x,.158,p['level'],8.5,M)
 text(x,.130,f"[{i+1}] doi: {p['doi']}",8.1,M,url=p['url'])
rect(.05,.047,.90,.062,'#EAF2EC')
text(.064,.078,'本项目切入点',13,G,True)
text(.194,.078,'从 CAMP 代谢关联出发，检验四癌中患者方向、恶性上皮表达与空间定位是否一致。',12,G)
text(.05,.025,'机制链为原文结果概括；本项目尚未验证上述底物、酶活或因果机制。',8.4,M)
text(.95,.025,'研究进展',9,G,True,ha='right')
stem='LYPLA1_近期研究进展'
for ext in ['png','pdf','svg']:f.savefig(OUT/(stem+'.'+ext),dpi=180,facecolor='white');shutil.copy2(OUT/(stem+'.'+ext),VIEW/(stem+'.'+ext))
plt.close(f)
with (OUT/'literature.tsv').open('w',encoding='utf-8',newline='') as h:
 writer=csv.DictWriter(h,fieldnames=['journal','year','date','model','title','doi','url','level','key','nodes','body','limit'],delimiter='\t');writer.writeheader();writer.writerows(papers)
(OUT/'analysis_spec.json').write_text(json.dumps(dict(version='v1',scope='literature presentation only',retrieved='2026-09-29',statistics='NOT_RUN: no new tests',data='published findings; no patient measurements',software=matplotlib.__version__,layout='16:9 SYSU green; three source-linked studies'),ensure_ascii=False,indent=2),encoding='utf-8')
with (OUT/'source_manifest.tsv').open('w',encoding='utf-8',newline='') as h:
 writer=csv.writer(h,delimiter='\t');writer.writerow(['source','version','sha256','reason'])
 for p in papers:writer.writerow([p['url'],p['date'],'NA','Live web text checked; no complete source file downloaded'])
 writer.writerow([str(logo.relative_to(ROOT)),'existing asset',hashlib.sha256(logo.read_bytes()).hexdigest(),'Existing project logo'])
(OUT/'README_CN.md').write_text('''# LYPLA1 近期研究进展页
## 本轮问题
用较新的原始研究替换老文献，制作一页与组会版式一致的文献背景。
## 输入与范围
JLR 2026（2025-11 在线）、BCB 2026、npj Precision Oncology 2025。核对来源见 literature.tsv。BCB 仅获得摘要，另两篇核对全文；未进行穷尽性系统综述。
## 实际结果
一页 PNG/PDF/SVG 与网页预览；文章标题和 DOI 均为真实引用，无生成式论文截图。当前为文献展示，不生成 PPTX。
## 新手解释
左侧连接乳腺癌脂质代谢与铁死亡；中间直接涉及 LYPLA1 功能；右侧提供肺癌耐药机制实例。箭头概括各论文报告的关系，不代表本项目验证。
## 限制/反证
JLR 中 MTDH 为研究中心；BCB 摘要显示模型依赖性；肺癌结果不可直接外推本项目四癌。丰度不等于酶活。所选三篇不称最高分或系统综述排名。没有制造统计值或新增实验结果。
## 当前决定
作为文献背景页，衔接跨癌、患者、细胞与空间证据，不声称首次发现 LYPLA1 与癌症相关。
## 下一步
组会解释研究增量时聚焦已获得的跨癌证据和适用边界。
## 复现命令
python code/ppt_lypla1_recent_literature_v1.py
''',encoding='utf-8')
links=''.join(f'<li><a href="{p["url"]}" target="_blank">{p["journal"]} · {p["year"]} — {p["title"].replace(chr(10)," ")}</a>（{p["level"]}）</li>' for p in papers)
html=f'''<!doctype html><meta charset="utf-8"><title>LYPLA1 近期研究进展</title><style>body{{margin:0;background:#e9eeeb;color:#183d34;font-family:"Microsoft YaHei",sans-serif}}nav{{padding:18px 4%;background:white}}main{{max-width:1680px;margin:24px auto}}img{{width:100%;display:block}}a{{color:#00553b}}section{{padding:20px;background:white;line-height:1.8}}li{{margin:10px 0}}</style><nav>LYPLA1 · 近期研究进展　<a href="{stem}.pdf">PDF</a>　<a href="{stem}.png">高清 PNG</a>　<a href="{stem}.svg">SVG</a></nav><main><img src="{stem}.png"><section><b>原始论文与证据范围</b><ol>{links}</ol><p>机制示意为原文结果概括，不是本项目已验证的机制。BCB 依据摘要；JLR 正式卷期为 2026 年，2025 年底在线。</p></section></main>'''
(VIEW/'index.html').write_text(html,encoding='utf-8');(OUT/'index.html').write_text(html,encoding='utf-8')
(OUT/'validation.json').write_text(json.dumps(dict(status='PARTIAL',render='PASS',visual_review='pending',sources='2 full texts; 1 abstract',new_tests=False),indent=2),encoding='utf-8')
print(VIEW)
