"""Create transparent PDAC-only delivery notes and append stage-index batches."""
import argparse,csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RUNS={'internal':('03_PATIENT','20260922T053500Z_internal_paired51_v2','internal_paired51_v2.py'),'sc':('06_EXTERNAL','20260922T054000Z_sc_paired51_v2','sc_source_paired51_v2.py'),'integration':('07_INTEGRATION','20260922T055500Z_internal_sc_paired51_v2','integrate_paired51_v2.py')}
def index(stage,run,path,code,scope,reason,next_action):
 p=ROOT/'coordination/stages/PDAC.tsv';rows=list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'));fields=list(rows[0]);rows=[x for x in rows if not (x['stage_id']==stage and x['run_id']==run)];rows.append(dict(zip(fields,['PDAC',stage,run,'PDAC_paired51_internal_sc_v2','DONE',scope,path,'code/pdac/'+code,'analysis/pdac-initial',reason,next_action])))
 with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(sorted(rows,key=lambda x:(x['stage_id'],x['run_id'])))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('kind',choices=list(RUNS));a=ap.parse_args();stage,run,code=RUNS[a.kind];out=ROOT/'results/PDAC'/stage/run;s=json.loads((out/'summary.json').read_text());relative=str(out.relative_to(ROOT)).replace('\\','/')
 if a.kind=='internal':
  lines=['|关系层与分析|计划关系|可评估|P<0.05|q<0.05|','|---|---:|---:|---:|---:|']+[f"|{k}|{v['planned']}|{v['evaluable']}|{v['P_lt005']}|{v['q_lt005']}|" for k,v in s['families'].items()]
  text=f'''# PDAC 21 个作者确认单位关联及 11 对 RNA 背景

## 本轮问题
按配对51项探索工作池，先完成CAMP内部关系和配对RNA，再进入三队列来源；不用旧27标本结果冒充新样本分析。

## 输入与范围
21个有作者明确且互不重复配对编号的CAMP肿瘤，另6个肿瘤未默认独立。11组同时有CAMP肿瘤/正常代谢物及RNA的作者配对用于RNA背景。固定357直接+358条件关系，538基因。四个实际来源路径及哈希、完整RNA标签均通过入口检查。

## 实际结果
{chr(10).join(lines)}

配对RNA：{s['paired_RNA_evaluable']}/538可评估，P<0.05为{s['paired_RNA_P_lt005']}，q<0.05为{s['paired_RNA_q_lt005']}。全部1430行关联和538行RNA保留，含NA。sample_identity_audit.tsv为汇总核查，完整患者映射只留服务器。

## 新手解释
关联问的是肿瘤患者之间RNA和代谢物是否共同变化；配对RNA问的是同患者肿瘤相对正常怎样变化，两者不能互相替代。RNA差异不显著不删除基因。直接和条件层各自BH，不能把两个层的q拼成一个统一校正家族。同输入的可用性结果是复用，不是独立验证。

## 限制/反证
51个代谢物均为配对原始P<0.05的探索性工作池，没有配对q<0.05。作者处理值可能依赖填补；可用性子集、留一法范围和符号变化逐行保留。作者配对表证明本次单位定义，不等于额外核实全部临床患者身份。协变量、因果机制、原冻结CAMP对比编码、外部联合组学未补齐。

首次预检目录20260922T053000Z_internal_paired51_v2因Windows生成的路径分隔符与实际Linux来源路径不一致被入口阻断，未计算统计。本目录规范化清单后运行。服务器4项入口/别名/BH负例测试通过；independent_validation.json记录重新读取来源的n/rho核对、RNA动态规划P复算与5个BH家族复算范围，不声称独立重放全部置换P。

## 当前决定
该数值批次DONE；03/04总体仍PARTIAL，不能认定机制或整癌种完成。全部候选进入三队列来源。

## 下一步
生成candidate_pre_scRNA.tsv，随后追加538基因三队列来源、第二来源、重采样与共同类别一致性。

## 复现命令
在新隔离server165目录放置源码、三个锁定映射文件及四来源清单：`python internal_paired51_v2.py --code-commit <base_commit>`；`python verify_internal_paired51_v2.py <run_directory>`。source_manifest.tsv记录执行代码哈希，analysis_spec.json记录参数与版本。矩阵不导出。
'''
  scope='715 relations;21 author-key tumors;11 paired RNA;538 genes;all rows retained';reason='Numerical batch complete;covariates,causal/external verification remain PARTIAL'
  index('04_ROBUSTNESS',run,relative,code,'Availability and leave-one-out for all715 relationships',reason,'Keep sensitivity flags;no rule changes from significance')
 elif a.kind=='sc':
  lines=['|队列|肿瘤细胞|作者供者/样本单位|可唯一对应基因|复用基因|重新提取/别名核查基因|','|---|---:|---:|---:|---:|---:|']+[f"|{v['cohort']}|{v['cells_included']}|{v['units_included']}|{v['genes_available']}|{v['genes_reused']}|{v['genes_new_or_alias_rechecked']}|" for v in s['cohorts']]
  text=f'''# PDAC 全538基因三队列来源与稳定性

## 本轮问题
完成CAMP内部关联/配对RNA后，全部538个候选做细胞来源，不以关联或RNA显著性筛选输入。

## 输入与范围
GSE263733原发肿瘤Pm0/Pm1，GSE278688肿瘤组织（不含PBMC/邻近正常），GSE242230作者eusfnb原始肿瘤样本、filtered=False且有细胞注释者。复用前批未改变基因的供者汇总，其余基因/别名重新读原矩阵；每队列MGLL/PNP/PEPD/GGT5哨兵基因重新核对逐单位计数、细胞数、UMI分母和检出率。

## 实际结果
{chr(10).join(lines)}

538个基因三队列均有状态行。三队列全类别唯一第一来源相同：{s['three_cohort_same_top']}；同时满足各队列bootstrap≥0.70及检出≥0.05：{s['stable_source_genes']}；共同合格类别中第一来源一致：{s['shared_category_same_top']}。不能把描述性一致称为疾病差异显著复现。summary.json中的genes_available为复用加新提取的总覆盖，genes_in_matrix仅指新提取子集；二者不能混用。

## 新手解释
先对每个供者×细胞类型汇总，再让供者等权，避免细胞更多的样本支配平均值。重采样以整个供者为单位，保留不同细胞类型来自同一供者的关系。第二来源、检出率和支持供者数同时展示；共同类别交集至少两类才比较。

## 限制/反证
仅沿用作者细胞注释，未重新聚类、做恶性CNV判定或环境RNA校正。Epithelial/Ductal不能直接叫恶性细胞；红细胞等来源可能受注释/环境RNA影响。基因缺失、重复标签或低检出是不可评估，不作零表达。独立研究来源已区分，但没有受保护患者身份级的跨研究去重。没有新的疾病差异P/q，也没有检验代谢物—基因关系外部复现。

## 当前决定
全候选来源批次DONE；06总体仍PARTIAL。数据提取、单位分母与源码哈希检查通过；独立供者汇总/最高来源/一致性核对范围详见independent_validation.json。未重新读取所有复用基因的原始计数。

## 下一步
把来源附加到完整关系表；以A/B/C/D标记探索优先级和数据缺口，随后才有针对性寻找外部与干预证据。

## 复现命令
新server165隔离目录放置固定panel、别名与依赖脚本：`python sc_source_paired51_v2.py --code-commit <commit>`；`python verify_sc_paired51_v2.py <run_directory>`。每批来源哈希、实际单位数和参数均随表提供。
'''
  scope='538 genes across3 cohorts;donor-balanced top/second/bootstrap/shared-category sources';reason='Expression-source batch complete;not same-relation external validation'
 else:
  classes='，'.join(f"{k}类{v}行" for k,v in sorted(s['candidate_classes'].items()));labels='，'.join(f"{k}: {v}" for k,v in sorted(s['relation_labels'].items()))
  text=f'''# PDAC 配对工作池的内部证据与单细胞整合

## 本轮问题
按用户批准顺序完成配对发现、关系映射、CAMP内部关联及配对RNA、三队列来源，形成下一步候选；外部matched-omics与功能深化后置。

## 输入与范围
全部307项配对结果保留；51项P<0.05工作池均属探索性。固定715条关系、538基因，另7项未映射特征占位，共722行。旧303/42、旧362关系和旧P/q不覆盖。

## 实际结果
candidate_pre_scRNA.tsv先保存代谢物、关联、RNA及缺失背景；candidate_scRNA_appended.tsv附加三队列来源后保持原字段不变。

关系标签（含7项未映射占位）：{labels}。

候选分层：{classes}。计数单位是关系/占位行，不是独立基因或发现数。同一基因对应不同代谢物可以属于不同类别。

## 新手解释
A为内部条件满足且稳定上皮/导管来源，B为同样内部条件但稳定非导管来源的组织组成/微环境背景；Acinar/Endocrine属于胰腺上皮实质，不能称为免疫/间质，source_context另作说明。C保留生化候选及提示/弱关联或来源未稳定者，不能统称阴性；D标明条件映射、缺测、敏感性、低检出等具体缺口。RNA差异不显著不会被删除。A/B仍是探索性，不等于治疗靶点、恶性细胞特异或机制验证。

## 限制/反证
所有工作池代谢物的配对q均≥0.05。分层不改变任何原P/q、检验族或入场规则。细胞来源不是同一代谢物—基因关系的外部验证；exact_relation_external_validation和functional_perturbation均为NOT_PERFORMED。CAMP历史效应编码、临床协变量、部分化学身份、研究间临床身份去重仍有缺口。低rho的留一法符号变化也按预设标记敏感，不能把D等同不存在生化关系。

## 当前决定
本轮全候选整合批次DONE；07总体PARTIAL。没有为得到A/B候选调整阈值。完整保留NA、条件关系、反证和同输入复用标记。

## 下一步
查看候选及明确缺口，优先核实具体底物特异性与患者支持；对值得继续的关系寻找外部matched-omics及功能干预，不重跑冻结CAMP。

## 复现命令
`python code/pdac/integrate_paired51_v2.py --run-id {run}` 后执行同命令加 `--append-sc`。首次文件先于单细胞附加保存，重复构建必须逐字一致；全部来源及检验规则有版本记录。
'''
  scope='722 rows;715 relationships+7 unresolved features;51 metabolites and538 genes';reason='Exploratory internal integration complete;causal and external validation not performed'
 (out/'README_CN.md').write_text(text,encoding='utf-8');index(stage,run,relative,code,scope,reason,'Review explicit candidate gaps;targeted external/function follow-up')
 print(relative)
if __name__=='__main__':main()
