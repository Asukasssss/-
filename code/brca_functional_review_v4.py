#!/usr/bin/env python3
"""Append manually curated evidence to the pinned BRCA aggregate tables.
No patient matrix access, statistical test, new P/q value, or DepMap request.
Run inside a checked-out repository after copying this delivery's files in.
"""
from __future__ import annotations
import argparse,csv,json,math,hashlib,subprocess
from pathlib import Path
from collections import Counter
RUN='20260919T120437Z_functional_review_v4'
BASE_COMMIT='227245e3171ac3b8530835624711c1a0a69363cf'
OLD='results/BRCA/20260919_ER_availability_v3/candidate_comparison_117_integrated_v3.tsv'
EDGES='results/BRCA/20260919_ER_availability_v3/er_availability_all_174.tsv'
DELIVERY=f'results/BRCA/05_FUNCTION/{RUN}'

def read_tsv(p:Path):
 with p.open(encoding='utf-8-sig',newline='') as f:
  rd=csv.DictReader(f,delimiter='\t');return list(rd),rd.fieldnames or []

def write_tsv(p:Path,rows,headers):
 if p.exists():raise FileExistsError(f'不会覆盖既有输出：{p}')
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,headers,delimiter='\t');w.writeheader();w.writerows(rows)

def number(s):
 return None if s is None or s=='' or s=='NA' else float(s)

def same(a,b):
 a,b=number(a),number(b)
 return a is None and b is None or a is not None and b is not None and math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--repo-root',type=Path,required=True)
 ap.add_argument('--out',type=Path,help='必须为不存在的新TSV；默认写当前交付目录新文件')
 ap.add_argument('--register-stage',action='store_true',help='校验通过后只追加本次BRCA阶段记录，不改变其他癌种')
 args=ap.parse_args();root=args.repo_root.resolve();d=root/DELIVERY
 src=root/OLD;edges=root/EDGES
 if not src.is_file() or not edges.is_file():
  raise FileNotFoundError('请先使用含BRCA v3结果的分支；此脚本不下载患者数据或猜测文件。')
 baseline,old_fields=read_tsv(src)
 review,review_fields=read_tsv(d/'candidate_comparison_117_functional_v4.tsv')
 er,_=read_tsv(edges); compact,_=read_tsv(d/'relations_174_annotated_v4.tsv')
 if len(baseline)!=117 or len(review)!=117 or len(er)!=174 or len(compact)!=174:
  raise ValueError('输入不是锁定的117基因/174关系版本，请先核查，禁止静默混版本。')
 oldby={r['gene']:r for r in baseline};newby={r['gene']:r for r in review}
 if len(oldby)!=117 or set(oldby)!=set(newby):raise ValueError('基因范围不一致。')
 mismatch=[]
 for g,r in newby.items():
  for nf,of in [('original_rho','best_rho'),('original_q','best_q'),('original_n','best_n')]:
   if not same(r[nf],oldby[g][of]):mismatch.append(f'{g}:{nf}')
  if g!='CHKB' and r['original_best_metabolite']!=oldby[g]['best_metabolite']:mismatch.append(f'{g}:metabolite')
 eby={r['relation_id']:r for r in er}
 if len(eby)!=174:raise ValueError('原ER表含重复关系。')
 for r in compact:
  oldr=eby.get(r['relation_id'])
  if oldr is None: mismatch.append(r['relation_id']);continue
  for nf,of in [('n_specimens','n'),('er_available_partial_rank_rho','partial_rank_rho'),('er_available_q','q_value'),('er_primary_q','primary_er_q')]:
   if not same(r[nf],oldr[of]):mismatch.append(f"{r['relation_id']}:{nf}")
  if r['gene']!=oldr['gene'] or r['metabolite']!=oldr['metabolite']:mismatch.append(f"{r['relation_id']}:identity")
 if mismatch:raise ValueError('紧凑摘录与仓库原表不一致；不写输出：'+', '.join(mismatch))
 hashes_before={str(p.relative_to(root)):digest(p) for p in [src,edges]}
 # ALL historical columns are preserved verbatim. New review columns carry a namespace.
 append_fields=[k for k in review_fields if k not in ('cancer','gene','snapshot_commit')]
 names=['review_v4_'+k for k in append_fields]
 if set(names)&set(old_fields):raise ValueError('该表已经包含v4列，不重复合并。')
 rows=[]
 for oldr in baseline:
  nr=dict(oldr);r=newby[oldr['gene']]
  nr.update({'review_v4_'+k:r[k] for k in append_fields});rows.append(nr)
 out=args.out or d/'candidate_comparison_117_ALL_LEGACY_COLUMNS_v4.tsv'
 write_tsv(out,rows,old_fields+names)
 for p in [src,edges]:
  if digest(p)!=hashes_before[str(p.relative_to(root))]:raise RuntimeError('源文件发生变化。')
 val={'passed':True,'legacy_rows':117,'legacy_columns_preserved':len(old_fields),'new_evidence_columns':len(names),'relations_checked':174,'numeric_compact_matches_sources':True,'new_statistical_tests':0,'source_sha256':hashes_before,'baseline_commit':BASE_COMMIT,'output':str(out),'output_sha256':digest(out)}
 vp=out.with_suffix('.validation.json');vp.write_text(json.dumps(val,ensure_ascii=False,indent=2),encoding='utf8')
 if args.register_stage:
  index=root/'coordination/stages/BRCA.tsv';prev,head=read_tsv(index)
  updates,_=read_tsv(d/'stage_index_append.tsv')
  existing={(r['stage_id'],r['run_id']) for r in prev}
  additions=[r for r in updates if (r['stage_id'],r['run_id']) not in existing]
  if any(r['cancer']!='BRCA' for r in updates):raise ValueError('拒绝修改其他癌种阶段。')
  combined=sorted(prev+additions,key=lambda r:(r['stage_id'],r['run_id']))
  tmp=index.with_suffix('.v4.tmp')
  with tmp.open('w',encoding='utf-8',newline='') as f:
   w=csv.DictWriter(f,head,delimiter='\t');w.writeheader();w.writerows(combined)
  tmp.replace(index)
 print(json.dumps(val,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
