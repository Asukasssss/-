"""Versioned display correction of published aggregates; no statistical fitting."""
from pathlib import Path
import csv,json,hashlib,ast,subprocess,collections
R=Path(__file__).resolve().parents[1]
RUN='20260922T120000Z_source_display_v2'
O=R/'results/BRCA/07_INTEGRATION'/RUN;O.mkdir(parents=True,exist_ok=True)
P='results/BRCA/07_INTEGRATION/20260922T112000Z_sop_record_v1/'
S='results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/'
N='results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/'
L='results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/'
A='results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/'
inputs=[]
def read(p):
 inputs.append(p)
 with (R/p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
 with (O/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def txt(name,s):(O/name).write_bytes((s.rstrip()+'\n').encode('utf-8'))
def js(name,o):txt(name,json.dumps(o,ensure_ascii=False,indent=2))
old=read(P+'candidate_genes156_reader.tsv')
full=read(P+'candidate_genes156_history_preserved.tsv')
relations=read(P+'candidate_relations262_history_preserved.tsv')
judgments=read(P+'interpretations156.tsv')
subs={x['gene']:x for x in read(S+'genes156_subtype_comparison.tsv')}
stability={};profiles={}
for p in [A+'all117_source_stability.tsv',N+'new39_source_stability.tsv']:
 for x in read(p):
  if x['partition']=='ALL' and x['cohort'] in ('Wu2021','Pal2021_reprocessed'):
   k=(x['cohort'],x['gene']);assert k not in stability,k;stability[k]=(x,p)
for folder in [L,N]:
 for cohort in ['Wu2021','Pal2021_reprocessed']:
  for x in read(folder+cohort+'_celltype_profiles.tsv'):
   if x['partition']=='ALL':
    k=(cohort,x['gene'],x['celltype']);assert k not in profiles,k;profiles[k]=x
node=next(x for x in ast.parse((R/'code/brca_sop_record_v1.py').read_text(encoding='utf-8')).body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='notes' for t in x.targets))
special=set(ast.literal_eval(node.value));inputs.append('code/brca_sop_record_v1.py')
reader=[];source_rows=[];changes=[]
for g in old:
 z=dict(g);gene=g['gene'];z['interpretation_scope']='专门讨论' if gene in special else '通用保留说明'
 z['interpretation_scope_note']='既有汇总中的定向讨论；不等于全文审计或最终排序' if gene in special else '未开展逐基因深入优先级裁决；未列示不表示低价值'
 for label,cohort in [('Wu','Wu2021'),('Pal','Pal2021_reprocessed')]:
  x,p=stability[(cohort,gene)]
  mechanical=x.get('top_lineage','NA');pr=profiles.get((cohort,gene,mechanical),{})
  valid=x['status']=='DONE' and pr.get('status')=='DONE'
  shown=mechanical if valid else '暂不可定位'
  z[label+'_top']=shown
  extra=dict(source_status=x['status'],source_reason=x['reason'],mechanical_top=mechanical,mechanical_top_profile_status=pr.get('status','NOT_EVALUABLE'),mechanical_top_eligible_source_labels=pr.get('n','NA'),mechanical_top_detection_fraction=pr.get('mean_detection_fraction','NA'),max_detection_fraction=x['max_detection'],n_evaluable_celltypes=x['n_evaluable_lineages'],bootstrap_top_frequency=x['bootstrap_top_frequency'],location_interpretable=valid,source_status_path=p)
  for k,v in extra.items():z[label+'_'+k]=v
  source_rows.append(dict(gene=gene,cohort=cohort,display_top=shown,**extra))
  if g[label+'_top']!=shown:changes.append(dict(gene=gene,field=label+'_top',old_display=g[label+'_top'],new_display=shown,reason=x['reason'],source=p))
 for short,subtype in [('ER','ER+'),('HER2','HER2+'),('TNBC','TNBC')]:
  x=subs[gene];field='Wu_'+short+'_top';status=x[subtype+'_status']
  z[field]=x[subtype+'_top'] if status=='DONE' else '暂不可定位'
  z['Wu_'+short+'_source_status']=status
  z['Wu_'+short+'_top_eligible_source_labels']=x[subtype+'_n']
  z['Wu_'+short+'_top_detection_fraction']=x[subtype+'_detection']
  z['Wu_'+short+'_shared_top']=x[subtype+'_shared_top'] if x[subtype+'_shared_status']=='DONE' else '暂不可定位'
  if g[field]!=z[field]:changes.append(dict(gene=gene,field=field,old_display=g[field],new_display=z[field],reason='subtype source status '+status,source=S+'genes156_subtype_comparison.tsv'))
 z['Wu_subtype_pattern']=subs[gene]['pattern']
 z['source_display_version']='source_display_v2'
 reader.append(z)
idx={x['gene']:x for x in reader}
write('candidate_genes156_reader.tsv',reader)
write('source_status312.tsv',source_rows)
write('display_changes.tsv',changes)
full_new=[]
for g in full:
 z=dict(g)
 for k,v in idx[g['gene']].items():
  if k not in old[0] or k.endswith('_top'):z['display_v2_'+k]=v
 assert all(z[k]==v for k,v in g.items())
 full_new.append(z)
write('candidate_genes156_history_preserved.tsv',full_new)
rr=[]
for r in relations:
 z=dict(r);g=idx[r['gene']]
 for k in ['Wu_top','Pal_top','Wu_source_status','Pal_source_status','interpretation_scope']:
  z['display_v2_'+k]=g[k]
 assert all(z[k]==v for k,v in r.items())
 rr.append(z)
write('candidate_relations262_history_preserved.tsv',rr)
for j in judgments:j['interpretation_scope']=idx[j['gene']]['interpretation_scope'];j['interpretation_scope_note']=idx[j['gene']]['interpretation_scope_note']
write('interpretations156.tsv',judgments)
sh=subs['SHMT2'];assert [sh[x+'_shared_top'] for x in ['ER+','HER2+','TNBC']]==['Myeloid','Malignant_epithelial','Malignant_epithelial']
write('SHMT2_source_excerpt.tsv',[sh])
assert idx['CYP4F2']['Wu_top']==idx['CYP4F2']['Pal_top']=='暂不可定位'
allowed={'Wu_top','Pal_top','Wu_ER_top','Wu_HER2_top','Wu_TNBC_top'}
for before,after in zip(old,reader):
 assert all(before[k]==after[k] for k in before if k not in allowed)
for r in source_rows:assert r['location_interpretable'] or r['display_top']=='暂不可定位'
assert len(reader)==156 and len(source_rows)==312 and len(rr)==262
base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
write('source_manifest.tsv',[dict(path=p,sha256=hashlib.sha256((R/p).read_bytes()).hexdigest(),source_commit=base) for p in sorted(set(inputs))])
js('analysis_spec.json',dict(version='source_display_v2',baseline='8663004',source_commit=base,operation='display_only',new_statistics=0,display_rule='source stability DONE and matching top-cell profile DONE; otherwise 暂不可定位',mechanical_top='retained separately; not biological localization',interpretation_scope='17 pre-existing named discussions versus generic retention; determined from original script notes',subtype='source status shown; shared-category top distinct from equal-subtype pooled top'))
js('validation.json',dict(status='DONE',gene_rows=156,relation_rows=262,source_rows=312,display_changes=len(changes),uninterpretable_source_rows=sum(not x['location_interpretable'] for x in source_rows),all_non_display_reader_fields_preserved=True,all_historical_master_columns_preserved=True,all_relation_historical_columns_preserved=True,all_original_P_q_unchanged=True,CYP4F2_masked_both_studies=True,SHMT2_shared_top_verified=True,explanation_counts=dict(collections.Counter(x['interpretation_scope'] for x in reader)),scope='aggregate display and value preservation;not full scientific review'))
txt('README_CN.md',f'''# BRCA来源展示修正 v2

## 本轮问题
避免机械最高类别被误读为可解释定位，并区分专门讨论与通用保留说明。
## 输入与范围
保留8663004基线。读取原156行、262关系及既有单细胞来源/分型汇总，不重算患者数据。
## 实际结果
全156行、262关系保留；{len(changes)}个展示字段更新。两研究312个基因来源记录中，{sum(not x['location_interpretable'] for x in source_rows)}项暂不可定位。实际供者覆盖和检测比例来自机械最高类的原汇总，列名明确，不把它称为已定位来源。
## 新手解释
Wu_top/Pal_top现为可解释展示；mechanical_top只保留原数值排序。NOT_EVALUABLE时即使某类均值最大也显示暂不可定位。状态、原因、有效来源标签数、检测比例和最高类别重采样频率并列。标签数量沿用作者定义，不声称重新确认患者身份。
## 限制/反证
CYP4F2两研究均NOT_EVALUABLE，不能称两研究恶性上皮定位一致。原最大检出比例约0.000514和0.000178（比例单位，不是百分数）。
SHMT2在三亚型自身类别及共同类别下均为ER+髓系、HER2+恶性上皮、TNBC恶性上皮；原表DIFFERENT_TOP_DESCRIPTIVE正确。三亚型等权汇总最高为恶性上皮，不等于每个亚型都如此。本项为解释更正，不修改原分型计算。
## 当前决定
17个专门讨论、139个通用保留说明。该标签仅表示文字来源/深度，不是证据强弱、候选排名或功能认证。AASS、ETNK1的数值与通用说明全部保留。无需给所有基因编机制故事。
## 下一步
使用新版易读表；旧底表的mechanical及record_v1字段仅作历史追溯，当前显示看display_v2。单细胞仍止于来源，不增加分析。
## 复现
python code/brca_source_display_fix_v2.py。重发用新RUN_ID。源哈希、字段差异、逐来源状态与验证一并保存。
''')
txt('.gitattributes','* -text')
write('checksums.tsv',[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(O.iterdir()) if p.is_file() and p.name!='checksums.tsv'])
reg=R/'coordination/stages/BRCA.tsv'
with reg.open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
if not any(x['run_id']==RUN for x in rows):
 rows.append(dict(cancer='BRCA',stage_id='07_INTEGRATION',run_id=RUN,analysis_version='source_display_v2',status='DONE',scope='156 genes262 relations;mask uninterpretable source ranks;scope interpretation labels;SHMT2 correction',result_path=O.relative_to(R).as_posix(),code_path='code/brca_source_display_fix_v2.py',git_branch='analysis/brca-functional-review-20260919',reason='Display only;all historical numeric fields preserved;no new tests',next_action='Read corrected source states and coverage;no mechanism expansion'))
 rows.sort(key=lambda x:(x['stage_id'],x['run_id']))
 with reg.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print((O/'validation.json').read_text(encoding='utf-8'))
