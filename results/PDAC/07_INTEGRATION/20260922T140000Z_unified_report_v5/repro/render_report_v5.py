"""PDAC presentation adapter over the unchanged pinned cross-cancer renderer.
Layout-only substitutions: missing-module pages, descriptive source page, and
Artifact Tool workbook export instead of the renderer's openpyxl authoring.
"""
from pathlib import Path
import sys,inspect,json,argparse,hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code'))
import camp_results_report as report

def emit_workbook_json(c,d,z,manifest,fm,sel,out):
 profiles=z['profiles'][['gene','cohort','celltype','effect','mean_detection_fraction','n','n_cells_total','n_cells_eligible','status','reason','effect_type','normalization']]
 book={'基因易读表':z['reader'],'完整关系715':d['all_relations'],'基因完整历史687':d['genes'],'代谢物全量':d['metabolites'],'直接关联全量':d['association'],'条件关联全量':d['conditional_association'],'当前RNA':d['rna'],'补充RNA':d['rna_history'],'单细胞来源全量':profiles,'来源状态':d['stability'],'跨研究共同类别':d['cross_study'],'进度与缺项':d['progress'],'缺项说明':d['gaps'],'展示选择':sel,'输入来源':manifest}
 data=[]
 for name,df in book.items():
  df=df.copy();numeric=[]
  for col in df:
   if any(x in col.lower() for x in ['id','key','sha','commit','path','version','gene','family','source','run','relation']):continue
   v=df[col].replace({'NA':np.nan,'':np.nan});n=pd.to_numeric(v,errors='coerce')
   if v.notna().any() and n.notna().sum()==v.notna().sum():
    df[col]=n;numeric.append(col)
  rows=df.astype(object).where(df.notna(),'NA').values.tolist()
  data.append(dict(name=name,headers=df.columns.tolist(),rows=rows,numeric=numeric))
 (out/'workbook_data.json').write_text(json.dumps(dict(input_commit=c['input_commit'],sheets=data),ensure_ascii=False,allow_nan=False),encoding='utf8')

EXTRA='''
 # PDAC missing modules remain visible; no made-up subtype or UMAP image.
 def gapplot(fig):
  ax=fig.add_subplot(111);ax.axis('off')
  rows=[[r.module,r.status,r.reason] for r in d['gaps'].itertuples()]
  tb=ax.table(cellText=rows,colLabels=['模块','状态','实际缺项'],loc='center',cellLoc='left',colWidths=[.21,.16,.63]);tb.auto_set_font_size(False);tb.set_fontsize(8);tb.scale(1,3.2)
 save('13_gaps','分型、UMAP及身份缺项',d['gaps'],gapplot,'gaps',[c['files']['gaps']])
 def sourcecover(fig):
  ax=fig.add_subplot(121);ss=z['source_status'];ct=ss.groupby('cohort').display_top.apply(lambda x:x.ne('暂不可定位').sum()).reindex(c['sc']['cohorts'])
  ax.barh(ct.index,ct.values,color='#087F8C');ax.barh(ct.index,687-ct.values,left=ct.values,color='#D6DADF');ax.set_xlabel('基因数；灰色=暂不可定位')
  for i,n in enumerate(ct):ax.text(n-12,i,str(n),ha='right',va='center',color='white')
  ax=fig.add_subplot(122);ax.bar(['三队列同首位','同首位且均≥80%','其中当前直接池'],[120,72,30],color=['#9AA6B2','#337EAA','#087F8C']);ax.set_ylabel('基因数');ax.tick_params(axis='x',labelsize=8);fig.tight_layout()
 save('13_source_status','来源覆盖与严格身份一致性',d['stability'],sourcecover,'stability',[c['files']['stability'],c['files']['cross_study']])
'''

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',default='configs/PDAC_report_v5.yaml');ap.add_argument('--out');ap.add_argument('--validate-only',action='store_true');args=ap.parse_args()
 c,d,b,m=report.load(args.config);z=report.prepare(c,d);v=report.validate(c,d,z)
 assert len(z['sensitivity'])==307 and len(z['association_sensitivity'])==357
 assert set(d['all_relations'].gene)<=set(z['genes'].gene)
 assert len(z['source_status'])==2061
 if args.validate_only:print(json.dumps(v,ensure_ascii=False));return
 if not args.out:ap.error('--out required')
 # Keep the imported implementation immutable and apply narrowly bounded presentation changes.
 src=inspect.getsource(report.build)
 begin=src.index(' wb=Workbook();wb.remove(wb.active)')
 end=src.index(' # Self-contained code/config package',begin)
 src=src[:begin]+' emit_workbook_json(c,d,z,manifest,fm,sel,out)\n'+src[end:]
 src=src.replace(' # The evidence matrix encodes',EXTRA+' # The evidence matrix encodes',1)
 src=src.replace("'历史保留'","'条件或历史补充'")
 report.emit_workbook_json=emit_workbook_json
 exec(compile(src,str(Path(__file__)), 'exec'),report.__dict__)
 out=Path(args.out);result=report.build(c,d,z,b,m,out,args.config)
 report.dump(out/'renderer_adapter_provenance.json',dict(base_renderer_sha256=hashlib.sha256((ROOT/'code/camp_results_report.py').read_bytes()).hexdigest(),handoff_commit=c['protocol_commit'],changes=['two missing/source coverage pages','non-direct pool label clarified','Artifact Tool workbook JSON instead of openpyxl export'],new_statistics=0))
 print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
