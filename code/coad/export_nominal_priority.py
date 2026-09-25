"""Create a user-requested nominal-P worklist from immutable public aggregates."""
import argparse,csv,hashlib,json,platform,subprocess
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,s):
    with p.open('x',encoding='utf-8',newline='\n') as f:f.write(s)
def dump(p,x):write(p,json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def table(p,fields,rows):
    with p.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def key(r):return tuple(r[k] for k in ['cohort','metabolite_key','human_gene_id'])
def nominal(r):return r['status']=='DONE' and float(r['p_value'])<.05
def same_direction(a,b):return b['status']=='DONE' and float(a['effect'])*float(b['effect'])>0
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-id',required=True);a=ap.parse_args()
    out=(ROOT/'results/COAD/07_INTEGRATION'/a.run_id).resolve();assert out.parent==(ROOT/'results/COAD/07_INTEGRATION').resolve()
    out.mkdir(parents=True,exist_ok=False);write(out/'.running',a.run_id+'\n')
    config=ROOT/'config/coad_nominal_priority_v1.json';spec=json.loads(config.read_text(encoding='utf-8'))
    files={k:ROOT/spec[k+'_source'] for k in ['primary','all37','rna','availability']}
    inputs=list(files.values())+[config,Path(__file__),ROOT/'templates/statistical_result.tsv']
    hashes={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in inputs}
    for p in files.values():
        validation=json.loads((p.parent/'validation.json').read_text(encoding='utf-8'))
        assert validation['status']=='PASS' and sha(p)==validation['public_output_sha256'][p.name]
    primary=read(files['primary']);all37={key(r):r for r in read(files['all37'])};availability={key(r):r for r in read(files['availability'])};rna={r['human_gene_id']:r for r in read(files['rna'])}
    assert len(primary)==len(all37)==len(availability)==974
    assert sum(r['status']=='DONE' for r in primary)==652
    original_fields=list(primary[0]);prefix=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t');assert original_fields[:len(prefix)]==prefix
    extra=['source_run_id','source_analysis_version','source_analysis_type','selection_status','evidence_label','all37_n','all37_rho','all37_p','all37_q','availability_n','availability_rho','availability_p','availability_q','availability_reuse_status','all37_nominal_same_direction','availability_nominal_same_direction','paired_RNA_effect','paired_RNA_p','paired_RNA_q','functional_review_status','external_review_status']
    fields=original_fields+extra;rows=[]
    for source in primary:
        row=dict(source);k=key(source);s37=all37[k];av=availability[k];rr=rna.get(source['human_gene_id'])
        row.update(stage_id=spec['stage_id'],run_id=a.run_id,analysis_version=spec['version'],analysis_type='primary33_nominal_priority_view',source_run_id=source['run_id'],source_analysis_version=source['analysis_version'],source_analysis_type=source['analysis_type'],functional_review_status='NOT_RUN',external_review_status='NOT_RUN')
        if nominal(source):row.update(selection_status='SELECTED_NOMINAL_P',evidence_label='NOMINAL_P_ONLY_NOT_FDR' if float(source['q_value'])>=.05 else 'FDR_SUPPORTED')
        elif source['status']=='DONE':row.update(selection_status='NOT_SELECTED_P_GE_005',evidence_label='NO_NOMINAL_ASSOCIATION_SUPPORT')
        else:row.update(selection_status='NOT_EVALUABLE' if source['status']=='NOT_EVALUABLE' else 'OUTSIDE_LOCKED_FAMILY',evidence_label='NOT_TESTED')
        for name,s in [('all37',s37),('availability',av)]:
            for target,field in [('n','n'),('rho','effect'),('p','p_value'),('q','q_value')]:row[name+'_'+target]=s[field]
            row[name+'_nominal_same_direction']=str(nominal(s) and same_direction(source,s)) if source['status']=='DONE' and s['status']=='DONE' else 'NA'
        row['availability_reuse_status']=av['reuse_status']
        for name,f in [('effect','effect'),('p','p_value'),('q','q_value')]:row['paired_RNA_'+name]=rr[f] if rr else 'NA'
        rows.append(row)
    chosen=[r for r in rows if r['selection_status']=='SELECTED_NOMINAL_P']
    assert len(chosen)==39 and all(r['evidence_label']=='NOMINAL_P_ONLY_NOT_FDR' for r in chosen)
    summary=dict(cancer='COAD',version=spec['version'],run_id=a.run_id,catalog_rows=974,planned_family=674,evaluable_family=652,selected_nominal_p=len(chosen),selected_unique_genes=len({r['human_gene_id'] for r in chosen}),selected_unique_features=len({r['metabolite_key'] for r in chosen}),positive=sum(float(r['effect'])>0 for r in chosen),negative=sum(float(r['effect'])<0 for r in chosen),selected_fdr_supported=0,
      all37_nominal_same_direction=sum(r['all37_nominal_same_direction']=='True' for r in chosen),availability_nominal_same_direction=sum(r['availability_nominal_same_direction']=='True' for r in chosen),both_sensitivities_nominal_same_direction=sum(r['all37_nominal_same_direction']==r['availability_nominal_same_direction']=='True' for r in chosen),selected_loo_sign_change=sum(r['loo_any_sign_change']=='True' for r in chosen),selected_paired_RNA_q_lt_005=sum(r['paired_RNA_q']!='NA' and float(r['paired_RNA_q'])<.05 for r in chosen),stage_status='PARTIAL',functional_review='NOT_RUN',external_review='NOT_RUN')
    spec.update(run_id=a.run_id,created_utc=datetime.now(timezone.utc).isoformat(),code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),software={'python':platform.python_version()},execution_location='local aggregate-only export; no server patient calculation',input_hashes=hashes)
    dump(out/'analysis_spec.json',spec);table(out/'results.tsv',fields,rows);table(out/'nominal_candidates.tsv',fields,chosen);dump(out/'summary.json',summary)
    table(out/'source_manifest.tsv',['source_path','sha256'],[dict(source_path=p,sha256=h) for p,h in hashes.items()])
    candidate_lines=['| 代谢物 | 基因 | rho | 主P | 主q | 全37 P | 可用性P |','|---|---|---:|---:|---:|---:|---:|']
    for r in chosen:candidate_lines.append('| '+r['metabolite_name']+' | '+r['gene']+' | '+' | '.join(format(float(r[f]),'.4g') for f in ['effect','p_value','q_value','all37_p','availability_p'])+' |')
    report=f'''# COAD 名义 P 候选工作批次

## 本轮问题

按用户在查看首轮结果后提出的“Q不显著就用P显著的”，将主分析 P<0.05 的关系用于探索性候选优先核查。记录筛选口径改变的时间和范围，不改变已有检验或数值。

## 输入与范围

复用33名I–IV期患者的主关联聚合表、全部37敏感性、当前可用性v2和配对RNA表。主家族计划674条、652条可计算；完整目录974行继续保留。本批只在本机读取已允许公开的聚合结果，无患者数据读取、无新检验和随机种子。

## 实际结果

选择 **39条主P<0.05关系，涉及{summary['selected_unique_genes']}个基因、{summary['selected_unique_features']}个代谢特征**；{summary['positive']}条正相关，{summary['negative']}条负相关。39条均未通过原652项家族的BH校正，主q最小为0.326。

其中{summary['all37_nominal_same_direction']}条在全部37分析中也为P<0.05且同向；{summary['availability_nominal_same_direction']}条在可用性分析中为P<0.05且同向；两项均满足的是{summary['both_sensitivities_nominal_same_direction']}条。逐一剔除后符号翻转{summary['selected_loo_sign_change']}条。敏感性列仅作描述，未作为额外入选条件；与主分析相同输入的结果直接复用，不能重复计作支持。

有{summary['selected_paired_RNA_q_lt_005']}条候选对应基因的配对RNA q<0.05；这是关系条数，不是独立基因数，也不证明该基因导致代谢物变化。

以下39条按共同的特征/基因顺序展示，不按P值另造最终靶点排名。完整974行见results.tsv；优先核查名单见nominal_candidates.tsv。

'''+ '\n'.join(candidate_lines)+f'''

## 新手解释

“名义P显著”表示单条检验达到P<0.05，未控制同时检验652条带来的多重比较风险；这份清单适合提出假设、安排下一步核查，不能称为已通过FDR的发现。原q保留，不能在39条子集中重新BH制造更小q。

## 限制/反证

这是看过结果后按用户要求采用的探索性优先级，不是预先固定的验证性成功标准。39条可能包含偶然关联，当前无法从这些P值识别哪些为真。样本量33，协变量、纯度/批次仍未完成；bootstrap区间为逐项区间。相同队列的敏感性及RNA差异不是独立验证。数据库反应注释不是COAD功能证据，RNA丰度不是酶活。其他613条可计算但P≥0.05的关系也不作功能阴性判定；22条RNA缺项和300条家族外待审记录完整保留。

## 当前决定

本轮39条列入探索性优先核查名单，标记NOMINAL_P_ONLY_NOT_FDR。首轮统计原值、q、家族和缺项不变。功能证据可以让名单外候选进入后续重点；本批清单导出DONE，07_INTEGRATION阶段PARTIAL，最终候选分层尚未完成。

## 下一步

优先核查这{summary['selected_unique_genes']}个基因的COAD/结直肠适用功能干预、直接代谢关系、反证和既有研究增量，同时保持完整候选池；随后结合预声明的协变量敏感性和适用外部证据判断是否保留为重点候选。

## 复现命令

`python code/coad/export_nominal_priority.py --run-id <新的唯一RUN_ID>`，无需服务器连接。输入与代码哈希见source_manifest.tsv；代码提交、规则及软件版本见analysis_spec.json。输出禁止覆盖已有目录。脚本逐列验证来源数值与状态未改、完整974唯一键、39条筛选条件、字段及排序，结果见validation.json。
'''
    write(out/'README_CN.md',report)
    # Validate exported bytes and all unchanged primary statistics, not only counts.
    check=read(out/'results.tsv');selected=read(out/'nominal_candidates.tsv')
    assert len(check)==974 and len({key(r) for r in check})==974 and len(selected)==39
    assert selected==[r for r in check if nominal(r)]
    assert list(check[0])[:len(prefix)]==prefix
    sortkeys=[tuple(r[f] for f in ['cohort','metabolite_name','metabolite_key','gene','human_gene_id']) for r in check];assert sortkeys==sorted(sortkeys)
    for source,row in zip(primary,check):
        for field in original_fields:
            if field not in ['stage_id','run_id','analysis_version','analysis_type']:assert row[field]==source[field],field
    assert [s[3:] for s in report.splitlines() if s.startswith('## ')]==['本轮问题','输入与范围','实际结果','新手解释','限制/反证','当前决定','下一步','复现命令']
    for p,h in hashes.items():assert sha(ROOT/p)==h
    outputs=['results.tsv','nominal_candidates.tsv','summary.json','analysis_spec.json','source_manifest.tsv','README_CN.md']
    dump(out/'validation.json',dict(status='PASS',checks=['source outputs match previously validated hashes','all974 unique rows retained','39 exact primary P<0.05 selections','all original numerical fields and missing statuses byte-for-byte preserved','common prefix and stable ordering','no BH recalculation or new tests','eight report headings','inputs unchanged'],public_output_sha256={f:sha(out/f) for f in outputs}))
    (out/'.running').unlink();write(out/'DONE',datetime.now(timezone.utc).isoformat()+'\n');print(json.dumps(summary,ensure_ascii=False))
if __name__=='__main__':main()
