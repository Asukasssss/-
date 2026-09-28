"""Document completed public-data work and explicitly blocked patient work."""
import csv,json,hashlib,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def dump(p,o):p.write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def report(p,title,sections):
    heads=['本轮问题','输入与范围','实际结果','新手解释','限制/反证','当前决定','下一步','复现命令']
    p.write_text('# '+title+'\n\n'+'\n\n'.join('## '+h+'\n'+s for h,s in zip(heads,sections))+'\n',encoding='utf-8')
def main():
    m=ROOT/'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1'
    f=ROOT/'results/PDAC/05_FUNCTION/20260920T112000Z_all_gene_search_v1'
    i=ROOT/'results/PDAC/07_INTEGRATION/20260920T112500Z_working_v1'
    e=ROOT/'results/PDAC/06_EXTERNAL/20260920T112500Z_all_gene_coverage_v1'
    p=ROOT/'results/PDAC/03_PATIENT/20260920T112500Z_first_round_prepared_v1';p.mkdir(parents=True,exist_ok=True)
    report(m/'README_CN.md','PDAC 全42特征直接关系池 v1',[
      '从完整42条原显著特征建立PDAC自己的关系池；不以BRCA候选或旧鸟苷两条关系为范围上限。',
      '冻结303条背景不变。读取当前KEGG化合物条目、Rhea小分子名称和官方pH对应、全reviewed人类UniProt注释；无患者数据。原名/键/q逐行保留。来源URL、时间与SHA256见mapping_sources.tsv。',
      '42条均有去向：25条找到进入主池的直接关系，9条仅条件候选，5条身份暂挂，3条本次范围未找到可靠精确关系。主池锁定362条去重关系、289个基因；另保留238个条件特征—基因组合。\n主池要求精确小分子参与reviewed人类反应，并有反应级ECO:0000269实验注释（或人类蛋白功能说明明确实验底物，如Gly-Pro—PEPD）。只有推断/泛类别注释不进入本版主池，保留后续复核。NAD+占125条关系，保留所有合规则记录并单列共因子背景，不作任意优先级截断。',
      'M/G由统一证据规则得到，不是显著性筛选。数据库实验注释不是PDAC干预实验，也不是逐篇原实验已复核。方程左右侧只是书写角色，不能推出肿瘤中实际通量。',
      "2'-O-methylguanosine原C04545实际指修饰tRNA，禁配游离分子到tRNA酶；合并ribulose/xylulose不拆峰；3'-AMP保留历史位置异构体暂挂；methionine sulfoxide未定S/R；1-stearoylglycerol原D01947为药物/混合条目。低特异性/推断关系和复合体角色另表；未找到关系不等于无代谢机制。原始化学身份未重新鉴定。",
      '本版映射和M/G锁定完成。全候选首轮尚未完成：服务器22端口连接超时，本轮未开展新的患者级统计。',
      '按固定362项分别开展主/可用性BH；不能按数据结果再改映射。条件关系新增需新版本。优先解决有明确来源的身份/底物疑点。',
      'python code/pdac/fetch_mapping_sources.py；python code/pdac/build_mapping_v1.py。使用已记录哈希的公共缓存；输出目录必须不存在。'])
    report(f/'README_CN.md','PDAC 全289基因功能文献初查',[
      '全候选都获得一次有记录的检索，不等关联显著后才查功能。',
      '对289个主池基因逐个检索Europe PMC标题/摘要：基因符号与pancreatic或PDAC；每个最多定位5条结果。查询、响应哈希、检索时间全保留；原摘要在本地公共缓存，不批量转载。',
      '289/289查询成功，202个有命中，定位686条基因—文章记录（不是686篇独立论文）。8条摘要经过人工阅读分层：MGLL、SIRT6有直接基因干预表型报告；其他分别保留广谱药物、上游基因干预、临床方案、非肿瘤模型或同名误命中。\n这些数量不代表所有基因已完成全文功能核查。',
      '检索命中、基因干预、当前精确代谢物机制是三个独立层。MGLL表型不能自动证明某个单酰甘油驱动肿瘤；SIRT6资料显示抑癌背景，也不能把所有候选统一解释为应抑制。',
      '短基因符号/旧别名可能漏检或同名命中。PNP-pincer化学配合物不是PNP基因；已显式排除该条支持。87个本次未命中不表示没有功能文献；下一轮需要扩展相关别名。关键词标记仅供定位，不能自动认证遗传干预。重点摘要阅读未替代全文。',
      '全候选统一检索定位DONE，疾病/关系功能证据审阅PARTIAL；DepMap保持暂缓。',
      '继续审阅有清晰遗传干预与代谢测量的重点来源，扩展短符号/旧别名；逐条保留反向、无效与不适用模型。',
      'python code/pdac/search_function_v1.py；摘要人工判断见abstract_review.tsv，不能由关键词自动重建。'])
    dump(f/'analysis_spec.json',{'version':'PDAC_all289_function_search_v1','gene_count':289,'page_size':5,'query_scope':'TITLE_ABS symbol AND pancreatic/PDAC','retrieval_source':'Europe PMC REST','manual_review_records':8,'fulltext_review':'NOT_RUN','DepMap':'DEFERRED','limitation':'Symbol query is not exhaustive;hit count not functional validation'})
    write(f/'source_manifest.tsv',[{'path':str(x.relative_to(ROOT)),'sha256':sha(x)} for x in [m/'direct_relations_v1.tsv',ROOT/'code/pdac/search_function_v1.py',f/'gene_search_ledger.tsv',f/'abstract_review.tsv']])
    dump(f/'validation.json',{'status':'PASS','all289_genes_searched':True,'unreviewed_not_called_functional':True,'manual_abstract_count':8})
    report(i/'README_CN.md','PDAC 关系级与基因级工作主表 v1',[
      '从映射阶段建立工作主表，不等所有数据完成，也不提前定稿机制或靶点。',
      '362条锁定关系、289个基因；关联、缺失、逐一剔除、功能、单细胞与CPTAC分列。历史两条鸟苷关联只放在明确命名的历史列。',
      '完整关系表362行、基因表289行。新患者分析字段均NA，状态ACCESS_BLOCKED；旧2项家族q没有冒充新362项家族q。无任意加权总分，无按关联显著性淘汰。',
      '这是持续更新的研究工作表；有功能背景不等于当前代谢物关系已证明，暂无统计不代表阴性。',
      '第一批仍未满足完整验收：服务器连接阻塞，362项主/可用性/逐一剔除尚未运行；全功能全文与外部扩展尚未完成。',
      '整合工作表建立DONE，整合结论PARTIAL；保留所有主池与独立条件表。',
      '恢复服务器后补齐全计划关系，随后更新工作分组；不提前挑2—4条“成功靶点”。',
      'python code/pdac/build_integration_v1.py；输出目录拒绝覆盖。'])
    report(e/'README_CN.md','PDAC 全候选外部覆盖台账 v1',[
      '将历史单细胞和CPTAC的范围与全289候选对齐，先明确已有与缺失。',
      '只读取已接收的14基因单细胞和13基因CPTAC汇总；本轮无新原矩阵读取，不推断新统计。',
      '289基因均列出：6个有历史两队列单细胞汇总，11个有历史CPTAC汇总；其余保留明确缺项。历史汇总存在不表示当前矩阵可评估。',
      '单细胞来源、肿瘤正常方向、显著性及代谢物—基因独立复现不同；导管细胞不能自动写成恶性细胞。',
      '原矩阵可评估性、配对/批次及CAMP独立性仍待核实。服务器当前不可达，不补造缺失基因的结果。外部联合组学可行性检索NOT_RUN。',
      '全G覆盖台账DONE；06_EXTERNAL总体PARTIAL。',
      '恢复服务器后按原处理和供体设计扩展全G；有限核查PDAC同标本联合组学，不使用BRCA FUSCC冒充验证。',
      'python code/pdac/build_integration_v1.py。'])
    report(p/'README_CN.md','PDAC 首轮统计准备与访问阻塞记录',[
      '按已锁定的362条关系准备完整主分析、可用性和逐一剔除。',
      '只准备代码和参数，本轮未读取患者矩阵。此前27肿瘤/12正常来源核对直接复用，不重复全包审计。',
      '本机连接server165超时，TCP22复测失败；已提示用户检查网络/VPN。四项统计单元检查通过：固定M且NA占位的BH、并列值/单调变换一致性、不可计算/逐一剔除、源哈希不一致阻止复用。没有新患者数值。',
      '代码通过合成数据检查不等于患者数据分析已经完成。独立患者未核实前，只能在显式标本独立性假设下探索。',
      '实际RNA别名覆盖、缺失、置换/bootstrap结果及协变量尚未计算。表达背景先只做描述，不默认肿瘤正常配对或两组完全独立。',
      'ACCESS_BLOCKED，禁止标DONE或生成伪造空数值。',
      '服务器恢复后复制代码和锁定输入到新PDAC/B运行目录，以独占锁执行。哈希不一致立即失败，不认证复用；同输入统计用缓存复用，新检验族单独算q。',
      'python code/pdac/test_patient_first_round.py；在新服务器运行目录执行python3 patient_first_round.py --code-commit <commit>。锁定输入为direct_relations_v1.tsv、mapping_spec.json、gene_aliases.json、历史输入哈希比较和历史关联表。'])
    dump(p/'analysis_spec.json',{'status':'ACCESS_BLOCKED','M':362,'G':289,'minimum_n':8,'permutations':9999,'bootstrap':4000,'BH':'fixed M362 separate primary/availability;public NA retained','independence':'unverified author specimens; explicit exploratory assumption','gene_expression':'descriptive until design confirmed','hash_mismatch':'fail closed','seed_reuse':'relation+actual arrays+mask+method cache'})
    dump(p/'validation.json',{'status':'PASS_CODE_TESTS_ONLY','unit_tests_passed':4,'patient_execution':'NOT_RUN_ACCESS_BLOCKED','python_local':platform.python_version()})
    write(p/'source_manifest.tsv',[{'path':str(x.relative_to(ROOT)),'sha256':sha(x)} for x in [m/'direct_relations_v1.tsv',m/'gene_aliases.json',ROOT/'code/pdac/patient_first_round.py',ROOT/'code/pdac/test_patient_first_round.py']])
    stages=ROOT/'coordination/stages/PDAC.tsv'
    with stages.open(encoding='utf-8') as h:r=csv.DictReader(h,delimiter='\t');fields=r.fieldnames;rows=list(r)
    for stage,path,status,scope,code in [('02_MAPPING',m,'DONE','42 features accounted;362 locked relations;289 genes;238 conditional pairs','build_mapping_v1.py'),('03_PATIENT',p,'ACCESS_BLOCKED','Locked-family statistics prepared;new patient tests not run','patient_first_round.py'),('04_ROBUSTNESS',p,'ACCESS_BLOCKED','Leave-one-out prepared;no patient results','patient_first_round.py'),('05_FUNCTION',f,'PARTIAL','289 searches;686 located gene-article records;8 abstracts reviewed','search_function_v1.py'),('06_EXTERNAL',e,'PARTIAL','All289 coverage ledger;6 historical cell-source and11 CPTAC overlaps','build_integration_v1.py'),('07_INTEGRATION',i,'PARTIAL','362 relationship and289 gene working rows;no arbitrary score','build_integration_v1.py')]:
        rows.append(dict(zip(fields,['PDAC',stage,path.name,'PDAC_complete_pool_v1',status,scope,str(path.relative_to(ROOT)).replace('\\','/'),'code/pdac/'+code,'analysis/pdac-initial','See batch report;first full numerical round remains incomplete','Restore server connection;complete all locked relationships'])))
    rows.sort(key=lambda r:(r['stage_id'],r['run_id']));write(stages,rows)
    print('Reports and stage indices written')
if __name__=='__main__':main()
