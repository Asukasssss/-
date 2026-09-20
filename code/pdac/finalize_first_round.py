"""Publish the accepted first-round scope and limits in standard stage reports."""
import csv,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
PAT=ROOT/'results/PDAC/03_PATIENT/20260920T113000Z_first_round_v1'
INT=ROOT/'results/PDAC/07_INTEGRATION/20260920T114000Z_first_round_v1'
FUN=ROOT/'results/PDAC/05_FUNCTION/20260920T114000Z_targeted_followup_v1'
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def write(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def report(p,title,sections):
    heads=['本轮问题','输入与范围','实际结果','新手解释','限制/反证','当前决定','下一步','复现命令']
    p.write_text('# '+title+'\n\n'+'\n\n'.join('## '+h+'\n'+s for h,s in zip(heads,sections))+'\n',encoding='utf-8')
def main():
    summary=json.loads((INT/'summary.json').read_text());rs=read(PAT/'associations.tsv');by={(r['relation_id'],r['analysis_type']):r for r in rs}
    sig=[r for r in rs if r['analysis_type']=='primary' and r['status']=='DONE' and float(r['q_value'])<.05]
    rows=[]
    for a in sig:
        b=by[(a['relation_id'],'availability')]
        rows.append({'feature_name':a['metabolite_name'],'gene':a['gene'],'primary_n':a['n'],'primary_rho':a['effect'],'primary_q':a['q_value'],'availability_n':b['n'],'availability_rho':b['effect'],'availability_q':b['q_value'],'interpretation':'Same input/rho/P/CI;family q differs' if a['input_hash']==b['input_hash'] else 'Fewer available specimens and weaker rho;do not call independently validated'})
    write(PAT/'primary_q_lt_005_summary.tsv',rows)
    table='| 特征—基因 | 主ρ（n=27） | 主q | 可用子集n/ρ/q |\n|---|---:|---:|---|\n'+'\n'.join(f"| {r['feature_name']}—{r['gene']} | {float(r['primary_rho']):.3f} | {float(r['primary_q']):.5f} | {r['availability_n']} / {float(r['availability_rho']):.3f} / {float(r['availability_q']):.5f} |" for r in rows)
    report(PAT/'README_CN.md','PDAC 全362关系首轮关联与基础稳健性',[
      '将已锁定的完整362条直接关系进行首轮探索，不以基因表达差异或关联P值决定功能检索范围。',
      '27个作者映射肿瘤标本内，作者处理微阵列与代谢矩阵无新增变换/填补。主与填补前可用性分族，各固定M=362；不可计算项公开P/q为NA，仅BH内部P=1占位。9999次双侧置换、4000次标本bootstrap、n≥8技术门槛。正常12标本仅用于独立描述性表达背景，不推断配对或独立两组。',
      '完整724行均保留：主359可计算、4条q<0.05；可用性235可计算、1条q<0.05。3个RNA基因缺失/不能唯一定位：DHRSX、IL4I1、PLA2G4B；另124条可用性联合n<8。\n'+table+'\n\n359条主分析完成逐一剔除；上述4条主关系剔除后均未变号，最大|Δρ|均<0.075。完整结果、实际n、缺失原因、并列值与区间见各表。2条旧主统计精确复用；主/可用性有163条输入相同，其中160条可算统计精确复用、3条为相同不可评估状态。',
      '相关系数不是倍数。主分析显著不能自动称“稳定机制”。MGLL和EPHX2在可用子集样本变少、效应减弱；GPC—PNPLA8的输入/ρ/P/区间完全相同，q从0.02715到0.05430只是其他关系的P分布改变导致的BH变化，不是该关系失稳。SLC6A19两分析输入相同，也不是独立复现。',
      '独立患者身份仍未分别证实；置换和bootstrap采用明确的标本独立性假设，仅作探索。高缺失和选择性可用性不能自动纠偏。纯度、组成、临床协变量未调整；源g方向未作新的组织方向推断。NAD+相关125条具有广泛共因子背景。四条主q显著关系不是四个靶点，其余候选继续保留。',
      '本批全计划关联、可用性及基础影响检查DONE；后续患者身份、协变量、功能/外部证据仍有缺项。此前网络阻塞已解除，准备批次记录作为历史保留。',
      '以完整362关系/289基因进入后续来源、功能和外部解释；依据实际元数据预设低维代理模型，不机械复制BRCA ER模型。不按这4条缩减全候选。',
      '服务器独占目录内python3 patient_first_round.py --code-commit 135c79be925c87ac69e1e5796896114e2e4a5cc6；接回公开汇总后python code/pdac/accept_patient_v1.py。本机合成数据测试：python code/pdac/test_patient_first_round.py。'])
    # Clarify that identical uncomputable states are not new statistical tests.
    for file in [PAT/'acceptance_validation.json',INT/'validation.json']:
        v=json.loads(file.read_text());v.pop('identical_input_primary_availability_statistics_reused',None);v.update(identical_input_records_reused=163,computable_identical_input_statistics_reused=160)
        dump(file,v)
    alias=read(FUN/'alias_search_ledger.tsv');alias_by={r['gene']:r for r in alias}
    genes=read(INT/'gene_working_table_v1.tsv');relations=read(INT/'relation_working_table_v1.tsv')
    for r in genes:
        a=alias_by.get(r['gene'],{})
        r['alias_followup_status']=a.get('status','NOT_REQUIRED_THIS_PASS');r['alias_followup_hits']=a.get('hits','NA')
        if r['gene']=='SIRT6':r['functional_status']='TARGETED_FULLTEXT_GENETIC_FUNCTION_SPECIFIC_MODELS';r['exact_metabolite_mechanism']='NAD_SIRT6_RELATION_NOT_ESTABLISHED'
    for r in relations:
        if r['gene']=='SIRT6':r['function_status']='TARGETED_FULLTEXT_GENETIC_FUNCTION_SPECIFIC_MODELS'
    write(INT/'gene_working_table_v1.tsv',genes);write(INT/'relation_working_table_v1.tsv',relations)
    report(INT/'README_CN.md','PDAC 完整候选首轮工作整合',[
      '把主、可用性、逐一剔除、功能与外部覆盖分别接回完整关系池，不压成任意总分。',
      '362条关系、289个基因；25个原显著特征有主池直接关系，其余17特征按条件/身份/本次未找到精确关系保留。另238个条件特征—基因组合不混入本版检验族。',
      '362行关系工作表与289行基因表完整。主359可算/4条q<0.05，可用性235可算/1条q<0.05；原q、历史2项家族q、新362项家族q分列。全289完成符号检索，另对87个初次零命中基因及PNP同名问题做有边界别名/全名补查：72个查询完成、28个补查对象有命中、定位51条记录；16个对象没有合适别名可查。已读8条摘要，并定向复核SIRT6正文。',
      '工作分组表示下一步去向，未按显著性剔除。MGLL有PANC-1敲低表型摘要，但不能据此认证每个单酰甘油机制；SIRT6存在指定模型中的抑癌遗传背景，即便关联不显著也保留。',
      '全文功能未覆盖全G，标本独立性未证实，组成/临床敏感性未做。全G外部矩阵扩展尚未执行：目前仅6基因有历史单细胞、11基因有历史CPTAC汇总；它们尚非独立联合组学验证。',
      '完整候选首轮数值和工作表已交付，07总体PARTIAL；不提前定稿机制、方向性干预或2—4条重点轴。',
      '继续全G单细胞/CPTAC可评估性与适用分析；有限核查PDAC外部同标本联合组学。全候选功能与统计解释同步完善，再选择少量具体问题深入。DepMap继续暂缓。',
      'python code/pdac/accept_patient_v1.py；python code/pdac/finalize_first_round.py。原数值在03_PATIENT只存一份，其余阶段引用。'])
    ledger=read(ROOT/'coordination/stages/PDAC.tsv');fields=list(ledger[0])
    for stage,run,path,status,scope,code in [
      ('03_PATIENT',PAT.name,PAT,'DONE','362 planned;359 primary evaluable;4 q<0.05;289 descriptive expression rows','patient_first_round.py'),
      ('04_ROBUSTNESS',PAT.name,PAT,'DONE','235 availability evaluable;1 q<0.05;359 primary leave-one-out;covariates not run','patient_first_round.py'),
      ('05_FUNCTION',FUN.name,FUN,'PARTIAL','8 initial abstracts;SIRT6 targeted fulltext;72 alias/full-name queries','search_alias_followup.py'),
      ('07_INTEGRATION',INT.name,INT,'PARTIAL','362 relation and289 gene views updated;full evidence follow-up pending','accept_patient_v1.py')]:
        ledger.append(dict(zip(fields,['PDAC',stage,run,'PDAC_first_round_v1',status,scope,str(path.relative_to(ROOT)).replace('\\','/'),'code/pdac/'+code,'analysis/pdac-initial','Exploratory author specimens;identity and later evidence limits retained','Continue all-candidate evidence and external coverage'])))
    ledger.sort(key=lambda r:(r['stage_id'],r['run_id']));write(ROOT/'coordination/stages/PDAC.tsv',ledger)
    assignments=read(ROOT/'coordination/assignments.tsv')
    for r in assignments:
        if r['cancer']=='PDAC':r['status']='FULL_POOL_FIRST_NUMERICS_DONE_FUNCTION_EXTERNAL_PARTIAL'
    write(ROOT/'coordination/assignments.tsv',assignments)
    report(ROOT/'reports/PDAC/CURRENT_STATUS_CN.md','PDAC 当前进度：完整候选首轮已交付',[
      '执行用户提供的BRCA借鉴方案；学习完整覆盖与证据分层，不复制BRCA候选/ER模型/显著名单。',
      '303条冻结效应、42条原显著特征不变。25条主池直接映射形成362关系/289基因，125关系涉及NAD+共因子背景；9条件、5身份暂挂、3本次未找到精确关系的特征全部保留。',
      '主359/362可算、4条q<0.05；可用性235/362可算、1条q<0.05。四条主结果为2-palmitoylglycerol—MGLL、9,10-DiHOME—EPHX2、GPC—PNPLA8、tryptophan—SLC6A19。\n当前入口：results/PDAC/03_PATIENT/20260920T113000Z_first_round_v1/README_CN.md；整合：results/PDAC/07_INTEGRATION/20260920T114000Z_first_round_v1/README_CN.md。',
      '不是四个治疗靶点。GPC—PNPLA8两种输入和原统计相同，q跨0.05仅由家族其他结果变化导致；MGLL/EPHX2可用子集效应减弱。所有其他候选仍在。',
      '27为作者映射肿瘤标本，不宣称独立患者；正常12仅描述性表达背景。患者身份、少量预设协变量、全G单细胞/CPTAC扩展与外部联合组学仍未完成。',
      '数值首轮及基础稳健性已完成；289基因检索/定位完成，8摘要+1正文已读，功能全面审阅仍PARTIAL。已维护全关系与全基因工作表。DepMap暂缓。网络中断已恢复，历史ACCESS_BLOCKED记录保留不作为当前状态。',
      '进入第二批解释背景：全候选来源/外部评估与适用模型，继续功能证据、反证及精确关系机制边界。完成后再选少数关系深入。',
      '代码、规范、来源哈希、测试及阶段索引均在本癌种目录；患者矩阵与完整映射仅留server165。上传/PR状态见coordination/publications/PDAC.json。'])
    print('Final first-round reports written')
if __name__=='__main__':main()
