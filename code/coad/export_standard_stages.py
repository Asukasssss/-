"""Adapt COAD summaries to the shared ordered schema; no statistical recomputation."""
import csv, hashlib, json, subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path): return json.loads((ROOT/path).read_text(encoding='utf-8'))
def dump(path,obj): path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def main():
    schema=read('config/stage_result_schema_v1.json')
    source=ROOT/'reference/camp/cancer_effects.tsv'
    effects=[r for r in csv.DictReader(source.open(encoding='utf-8-sig'),delimiter='\t') if r['cancer']=='COAD']
    by_key={r['metabolite_key']:r for r in effects}
    batches=[
      ('01_effects',effects,'reference/camp/cancer_effects.tsv','冻结效应记录描述；不重算效应和 q','效应记录',
       {'all_effects':159,'original_q_lt_005':73,'significant_positive':33,'significant_negative':40}),
      ('02_identity',read('reports/COAD/fresh_mapping_v0_1/identity_review.json'),'reports/COAD/fresh_mapping_v0_1/identity_review.json','数据库名称与标识初查；不代表源鉴定全部完成','显著代谢特征',{}),
      ('03_mapping',read('reports/COAD/current_catalog_v0_1/candidate_catalog.json'),'reports/COAD/current_catalog_v0_1/candidate_catalog.json','独立生化候选初稿；保留支持、暂挂、待审和不可评估','特征—人类基因组合',
       {'all_candidate_pairs':974,'reaction_supported_pairs':763,'supported_without_current_identity_hold':674,'features_without_current_identity_hold':59,'genes_without_current_identity_hold':458})]
    created=datetime.now(timezone.utc).isoformat()
    batch_id='20260919_B_v1'
    index=[]
    for sid,detail,input_path,scope,unit,counts in batches:
        output=ROOT/'results/COAD'/sid/batch_id
        if output.exists(): raise RuntimeError('Batch exists; use a new batch version: '+str(output))
        output.mkdir(parents=True)
        records=[]
        for r in detail:
            original=by_key[r['metabolite_key']]
            reason=r.get('disposition',r.get('identity_review_status','FROZEN_EFFECT_DESCRIBED'))
            status='DONE'
            if reason.startswith('HOLD_') or reason in ['IDENTITY_HOLD','REACTION_SPECIFICITY_PENDING','PROVISIONAL_NEW_KEGG_ID']:
                status='IN_PROGRESS'
            if reason=='UNIPROT_RECORD_NOT_EVALUABLE': status='NOT_EVALUABLE'
            record=dict(zip(schema['record_prefix'],[
                '1.0','COAD',sid,batch_id,r.get('dataset'),r['feature_name'],r['metabolite_key'],r.get('gene_symbol'),r.get('human_gene_id'),
                status,reason,original['hedges_g'],original['effect_fdr'],r]))
            records.append(record)
        records.sort(key=lambda r:tuple(r[k] or '' for k in schema['record_sort']))
        counts=counts|{'record_status_counts':dict(Counter(r['record_status'] for r in records))}
        inputs=[{'path':input_path,'sha256':digest(ROOT/input_path)},{'path':'reference/camp/cancer_effects.tsv','sha256':digest(source)}]
        publication={'status':'ACCESS_BLOCKED','branch':'analysis/coad-initial','remote_commit':None,'pr_url':None,
                     'reason':'Local Git has no usable HTTPS credentials; current GitHub connector reports push=false.'}
        interpretation='DONE 表示本批指定的描述/数据库初查完成；不等于原始鉴定、样本关联或功能验证。'
        limitations=['原冻结效应/原 q 不改；身份与底物范围未完全核实。','患者矩阵与凭据不在此包。','server165 本机别名无法解析，RNA/样本分析未做。']
        next_action='核对作者身份与样本覆盖，补查待审关系后确定下一版分析范围。'
        summary=dict(zip(schema['summary_field_order'],[
            '1.0','COAD','B',sid,batch_id,batch_id,'DONE',scope,created,
            subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),digest(Path(__file__)),inputs,unit,
            {'new_tests':'NOT_RUN','original_q':'Frozen effect_fdr; no recalculation'},len(detail),len(records),counts,interpretation,limitations,next_action,publication]))
        dump(output/'records.json',records); dump(output/'summary.json',summary)
        sections=[f"COAD / B / {sid} / {batch_id}\n\n批次状态：DONE（{scope}）。",scope+'。\n\n输入：'+input_path,
                  f'记录数：{len(records)}；单位：{unit}。\n\n'+json.dumps(counts,ensure_ascii=False),interpretation,
                  '公共字段、顺序、记录数、唯一键、输入及输出哈希由 tools/validate_stage_results.py 校验。',
                  '\n'.join('- '+x for x in limitations),next_action,
                  'ACCESS_BLOCKED：本地结果已准备，尚无可用写入凭据，不能称为已上传。后续发布状态查 coordination/publications/COAD.json。']
        report='# '+sid+' '+next(s['title_cn'] for s in schema['stages'] if s['stage_id']==sid)+'\n\n'
        report+='\n\n'.join('## '+title+'\n\n'+body for title,body in zip(schema['report_sections'],sections))+'\n'
        (output/'README_CN.md').write_text(report,encoding='utf-8')
        dump(output/'manifest.json',{'files':[{'path':name,'bytes':(output/name).stat().st_size,'sha256':digest(output/name)} for name in ['README_CN.md','summary.json','records.json']]})
        index.append({'stage_id':sid,'title_cn':next(s['title_cn'] for s in schema['stages'] if s['stage_id']==sid),
                      'stage_status':'DONE' if sid=='01_effects' else 'IN_PROGRESS','latest_batch':batch_id,
                      'batch_status':'DONE','path':str(output.relative_to(ROOT)).replace('\\','/')})
    for stage in schema['stages'][3:]:
        index.append(stage|{'stage_status':'ACCESS_BLOCKED' if stage['stage_id']=='04_patient' else 'NOT_RUN',
                           'latest_batch':None,'batch_status':None,'path':None})
    dump(ROOT/'results/COAD/stage_index.json',{'schema_version':'1.0','cancer':'COAD','owner':'B','stages':index})
    print('Exported three batches: 159 effect rows, 73 identity rows, 974 candidate pairs. No statistical values changed.')
if __name__=='__main__': main()
