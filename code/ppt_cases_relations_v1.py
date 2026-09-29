"""CAMP-only case pages from fixed public aggregates; no external cohort display."""
from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/BRCA/07_INTEGRATION';S=BASE/'20260928T150000Z_ppt_redraw_v1/sources';O=BASE/'20260929T_representative_cases_v1';O.mkdir(exist_ok=True)
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont('C:/Windows/Fonts/'+n)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';LINE='#DDE7E0';UP='#C77462';DOWN='#468B9E'
rows=[]
for c,g,key,molecule in [('BRCA','ASNS','KEGG:C00064','谷氨酰胺'),('COAD','UCKL1','KEGG:C00299','尿苷'),('PDAC','SLC6A6','KEGG:C00245','牛磺酸'),('PRAD','SLC6A6','KEGG:C00245','牛磺酸')]:
 d=pd.read_csv(S/(c.lower()+'_rel.tsv'),sep='\t');r=d[(d.gene==g)&(d.metabolite_key==key)].iloc[0];genes=pd.read_csv(S/(c.lower()+'_gene.tsv'),sep='\t');gr=genes[genes.gene==g].iloc[0]
 if c in ['BRCA','COAD']:
  m_n=int(r.metabolite_paired_n);m_up=int(r.metabolite_paired_pairs_higher);m_p=r.metabolite_paired_p_value;me=r.metabolite_paired_effect;rho=r.CAMP_effect;ap=r.CAMP_p_value;an=int(r.CAMP_n);lo=r.CAMP_ci_lower;hi=r.CAMP_ci_upper
  if c=='BRCA':rn=int(gr.v2_paired_RNA_n);ru=int(gr.v2_paired_RNA_positive_pairs);rp=gr.v2_paired_RNA_p_value;re=gr.v2_paired_RNA_effect
  else:rn=int(gr.RNA_n);ru=int(gr.RNA_n_up);rp=gr.RNA_p_value;re=gr.RNA_effect
 elif c=='PDAC':m_n=int(r.metabolite_n);m_up=int(r.metabolite_up_pairs);m_p=r.metabolite_p;me=r.metabolite_effect;rho=r.association_rho;ap=r.association_p;an=int(r.association_n);lo=r.association_ci_lower;hi=r.association_ci_upper;rn=int(r.RNA_n);ru=int(r.RNA_up_pairs);rp=r.RNA_p;re=r.RNA_effect
 else:m_n=int(r.metab_n);m_up=int(r.metab_n_up);m_p=r.metab_p_value;me=r.metab_effect;rho=r.patient_primary_effect;ap=r.patient_primary_p_value;an=int(r.patient_primary_n);lo=r.patient_primary_ci_lower;hi=r.patient_primary_ci_upper;rn=int(r.RNA_primary_n);ru=None;rp=r.RNA_primary_p_value;re=r.RNA_primary_effect
 rows.append(dict(cancer=c,gene=g,molecule=molecule,key=key,m_n=m_n,m_up=m_up,m_p=m_p,m_effect=me,rna_n=rn,rna_up=ru,rna_p=rp,rna_effect=re,rho=rho,association_p=ap,association_n=an,lo=lo,hi=hi))
pd.DataFrame(rows).to_csv(O/'CAMP_case_results.tsv',sep='\t',index=False)
def t(f,x,y,s,size=13,col=INK,bold=False,ha='left'):f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center',linespacing=1.6)
def line(f,y):f.add_artist(Line2D([.05,.95],[y,y],transform=f.transFigure,color=LINE,lw=.8))
for slug,title,sub,cs,conclusion in [
 ('01_BRCA_ASNS_CAMP','BRCA：谷氨酰胺降低，ASNS 升高并呈负相关','CAMP 患者匹配数据 · 谷氨酰胺—ASNS 酶反应注释',['BRCA'],'代谢差异与基因表达形成候选线索，下一步检查其细胞来源。'),
 ('04_PDAC_PRAD_SLC6A6_CAMP','牛磺酸—SLC6A6：PDAC 与 PRAD 呈现不同背景','同一代谢物与转运体 · 分别展示组织变化及肿瘤样本内关联',['PDAC','PRAD'],'PDAC 两者升高且正相关；PRAD 两者降低，但主分析相关未显著。'),
 ('11_COAD_UCKL1_CAMP','COAD：尿苷降低，UCKL1 升高并呈负相关','CAMP 患者匹配数据 · 尿苷—UCKL1 直接反应注释',['COAD'],'核苷代谢与 UCKL1 构成研究线索，单细胞与空间用于核对表达背景。')]:
 f=plt.figure(figsize=(16,9),facecolor='white');ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(BASE/'20260929T170000Z_page8_genes_six_v2/sysu_logo.png'));ax.axis('off');t(f,.95,.93,'CAMP  /  组会汇报',10,M,ha='right');line(f,.883);t(f,.05,.828,title,26,G,True);t(f,.05,.776,sub,12,M)
 for x,h in [(.17,'代谢物变化'),(.43,'基因 RNA 变化'),(.72,'肿瘤样本内关联')]:t(f,x,.685,h,17,G,True,ha='center')
 for i,c in enumerate(cs):
  r=next(v for v in rows if v['cancer']==c);y=.495 if len(cs)==1 else .548-i*.252
  if len(cs)>1:t(f,.055,y,c,15,G,True)
  for x,name,eff,p,n,up in [(.17,r['molecule'],r['m_effect'],r['m_p'],r['m_n'],r['m_up']),(.43,r['gene'],r['rna_effect'],r['rna_p'],r['rna_n'],r['rna_up'])]:
   direction='升高' if eff>0 else '降低';color=UP if eff>0 else DOWN;t(f,x,y+.040,name+' '+('↑' if eff>0 else '↓'),23,color,True,ha='center');t(f,x,y-.016,f'配对 P = {p:.3g}',13,INK,ha='center')
   same=(up if eff>0 else n-up) if up is not None else None
   t(f,x,y-.065,f'{same} / {n} 位患者{direction}' if same is not None else f'{n} 组患者配对 · 平均{direction}',11.5,M,ha='center')
  ax=f.add_axes([.59,y-.045,.27,.065]);ax.axvline(0,c='#BDC8C1',lw=.8);ax.plot([r['lo'],r['hi']],[0,0],lw=2.2,c=DOWN if r['rho']<0 else UP);ax.scatter([r['rho']],[0],s=95,edgecolor=DOWN if r['rho']<0 else UP,facecolor=(DOWN if r['rho']<0 else UP) if r['association_p']<.05 else 'white',zorder=3)
  ax.set_xlim(-.85,.85);ax.set_ylim(-1,1);ax.set_yticks([]);ax.set_xticks([-.5,0,.5]);ax.tick_params(labelsize=9);ax.spines[['top','left','right']].set_visible(False);ax.spines['bottom'].set_color(LINE)
  t(f,.725,y+.075,f"ρ = {r['rho']:+.3f}   P = {r['association_p']:.4g}",15,G,True,ha='center');t(f,.725,y-.094,f"{r['association_n']} 个肿瘤样本 · Spearman 相关",11.3,M,ha='center')
  if len(cs)==2 and i==0:line(f,.410)
 if len(cs)==1:
  f.add_artist(Rectangle((.09,.220),.82,.085,transform=f.transFigure,facecolor='#F0F5F1',edgecolor='none',zorder=-1));t(f,.50,.262,'代谢物差异  →  已知关系映射  →  匹配肿瘤样本内的实测关联',16,G,True,ha='center')
 line(f,.150);t(f,.05,.109,conclusion,14,G,True);t(f,.05,.067,'原统计复用，展示名义 P；相关横线为原 95% 区间。相关不等于因果或代谢通量，跨癌方向不同未作交互检验。',9,M)
 for ext in ['png','pdf','svg']:f.savefig(O/(slug+'.'+ext),dpi=180,facecolor='white')
 plt.close(f)
(O/'relation_source_manifest.tsv').write_text('source\tsha256\n'+''.join(f'{p.relative_to(ROOT)}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n' for c in ['brca','coad','pdac','prad'] for p in [S/(c+'_rel.tsv'),S/(c+'_gene.tsv')]),encoding='utf-8')
print('3 CAMP-only pages rendered')
