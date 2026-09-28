"""Read-only historical interpretation overlay; no statistical recomputation."""
import argparse, csv, hashlib, json, platform, subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
H=ROOT/'results/PDAC/03_PATIENT/20260920T054500Z_source_readiness_v1'
M=ROOT/'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1'
P=ROOT/'results/PDAC/03_PATIENT/20260920T113000Z_first_round_v1'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def table(p, rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader()
        w.writerows({k: ('NA' if v=='' or v is None else v) for k,v in r.items()} for r in rows)
def dump(p,obj): p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def review_label(r):
    original=r['final_gene_evidence_status']
    if original!='REPLICATED_TUMOR_ASSOCIATED_WITHIN_CELLTYPE':
        return original, 'Historical label retained, not newly validated; cell-source preference and disease difference remain separate.'
    qs=[float(r[f'{g}_top_celltype_fdr']) for g in ['GSE263733','GSE278688']]
    ds=[float(r[f'{g}_top_celltype_delta']) for g in ['GSE263733','GSE278688']]
    if ds[0]*ds[1]>0 and all(q>=.05 for q in qs):
        return 'DIRECTION_CONCORDANT_NOT_FDR_SIGNIFICANT', 'Mean effects have the same sign; neither listed within-celltype disease comparison passes FDR<0.05. Not significant replication; not proof of no effect. Median direction may differ.'
    raise ValueError('Unexpected historical replication label; requires explicit review')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-id',required=True);args=parser.parse_args()
    if Path(args.run_id).name!=args.run_id: raise ValueError('run-id must be one directory name')
    out=ROOT/'results/PDAC/06_EXTERNAL'/args.run_id
    out.mkdir(exist_ok=False);lock=out/'.running';lock.write_text('historical interpretation review\n')
    sources=[H/x for x in ['historical_cell_source_gene_summary.tsv','historical_cell_source_relations.tsv',
        'historical_cptac_effects.tsv','historical_hypothesis_lock.yaml','feature_coverage.tsv']]
    sources += [M/'direct_relations_v1.tsv',M/'feature_coverage_42.tsv',M/'analysis_spec.json',P/'associations.tsv',Path(__file__)]
    before={p:sha(p) for p in sources}
    genes=read(sources[0]);sc=read(sources[1]);cp=read(sources[2]);relations=read(M/'direct_relations_v1.tsv')
    reviewed=[]
    for r in genes:
        label,reason=review_label(r)
        row={'cancer':'PDAC','gene':r['gene'],'original_label':r['final_gene_evidence_status'],
             'reviewed_label':label,'review_reason':reason,'celltype_bias_original':r['celltype_bias_replicated']}
        for cohort in ['GSE263733','GSE278688']:
            for field in ['top_celltype','top_celltype_delta','top_celltype_median_delta','top_celltype_fdr','effect_design']:
                row[f'{cohort}_{field}']=r[f'{cohort}_{field}']
        reviewed.append(row)
    table(out/'historical_cell_source_label_review.tsv',reviewed)
    by_gene={r['gene']:r for r in reviewed}
    table(out/'historical_cell_source_relation_review.tsv',[dict(r,original_label=r['final_gene_evidence_status'],
        reviewed_label=by_gene[r['gene']]['reviewed_label'],review_reason=by_gene[r['gene']]['review_reason']) for r in sc])
    sg={r['gene'] for r in genes};cg={r['gene'] for r in cp}
    dump(out/'gene_set_overlap.json',{'singlecell_genes':sorted(sg),'cptac_genes':sorted(cg),
         'n_singlecell':len(sg),'n_cptac':len(cg),'intersection':sorted(sg&cg),
         'n_intersection':len(sg&cg),'meaning':'Historical gene sets only; does not establish cohort independence.'})
    table(out/'historical_cptac_scope_review.tsv',[{'cancer':'PDAC','gene':r['gene'],
        'historical_axis':r['axis'],'historical_relation':r['relation'],'historical_tier':r['tier'],
        'rna_status':r['rna_status'],'protein_status':r['protein_status'],
        'reviewed_scope':'GENE_ABUNDANCE_CONTEXT_ONLY',
        'exact_metabolite_gene_external_replication':'NOT_PERFORMED',
        'review_reason':'No joint metabolite-gene association measured here. Historical axis is not a confirmed feature-level relation. PNP deoxynucleoside hypothesis must not be renamed guanosine replication.'} for r in cp])
    cp_by={r['gene']:r for r in cp};exact={(r['metabolite'],r['gene']) for r in sc}
    evidence=[]
    for r in relations:
        for typ in ['GSE263733_CELL_SOURCE','GSE278688_CELL_SOURCE','CPTAC_RNA','CPTAC_PROTEIN']:
            is_sc=typ.startswith('GSE');covered=r['gene'] in (sg if is_sc else cg)
            c=cp_by.get(r['gene'],{})
            evidence.append({'cancer':'PDAC','cohort':r['cohort'],'feature_name':r['feature_name'],
                'metabolite_key':r['metabolite_key'],'gene':r['gene'],'relation_id':r['relation_id'],'evidence_type':typ,
                'historical_gene_covered':covered,
                'historical_feature_name_link':'EXACT_STRING_LINK_ONLY' if is_sc and (r['feature_name'],r['gene']) in exact else 'NOT_ESTABLISHED',
                'reviewed_cell_source_label':by_gene[r['gene']]['reviewed_label'] if is_sc and covered else 'NA',
                'historical_cptac_axis':c.get('axis','NA') if not is_sc else 'NA',
                'historical_measurement_status':('HISTORICAL_SUMMARY' if is_sc else c.get('rna_status' if typ=='CPTAC_RNA' else 'protein_status','NA')) if covered else 'NOT_COVERED',
                'exact_relation_external_validation':'NOT_PERFORMED',
                'evidence_scope':'GENE_EXPRESSION_CONTEXT_ONLY' if covered else 'NO_HISTORICAL_GENE_COVERAGE',
                'independence':'NOT_VERIFIED','reason':'Gene expression context cannot validate the actual metabolite-gene relation; no source-count score.'})
    table(out/'relation_evidence_review.tsv',evidence)
    coverage={r['feature_name']:r for r in read(H/'feature_coverage.tsv')};elig=[]
    stats=read(P/'associations.tsv')
    for r in read(M/'feature_coverage_42.tsv'):
        c=coverage[r['feature_name']];n=int(c['n_tumor_preimputation_available'])
        elig.append({'cancer':'PDAC','feature_name':r['feature_name'],'metabolite_key':r['metabolite_key'],
            'n_preimputation_available':n,'n_processed_finite':c['n_tumor_processed_finite'],
            'identity_status':r['identity_status'],'mapping_status':r['mapping_status'],
            'identity_note':r['chemical_identity_note'],'locked_primary_relations':r['n_direct_relations'],
            'primary_evaluable':sum(s['metabolite_name']==r['feature_name'] and s['analysis_type']=='primary' and s['status']=='DONE' for s in stats),
            'availability_evaluable':sum(s['metabolite_name']==r['feature_name'] and s['analysis_type']=='availability' and s['status']=='DONE' for s in stats),
            'availability_feature_n_gate':'N_LT_8' if n<8 else 'CHECK_PAIRWISE_FINITE_N_AND_VARIATION',
            'review_reason':'Locked v1 primary uses author-processed finite pairs, n>=8 and variable values; no preimputation minimum for primary. Availability uses preimputation-available subset. Not equal observation support; not MS detection certification.'})
    table(out/'feature_eligibility_review_42.tsv',elig)
    changed=[r['gene'] for r in reviewed if r['original_label']!=r['reviewed_label']]
    assert set(changed)=={'GGT5','MMP14','MGLL'}
    assert len(sg)==14 and len(cg)==13 and not sg&cg
    assert len(relations)==362 and len(elig)==42 and sum(int(r['locked_primary_relations']) for r in elig)==362
    assert len({(r['cohort'],r['feature_name'],r['gene'],r['evidence_type']) for r in evidence})==1448
    assert all(sha(p)==h for p,h in before.items())
    dump(out/'analysis_spec.json',{'version':'PDAC_historical_interpretation_review_v1','run_id':args.run_id,
        'code_base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'reviewed_user_baseline':'c508f94c80f2d68ba664a14d36bdefcd0e9f6cb9','unit':'published aggregate rows',
        'parameters':{'fdr_threshold_for_label':.05,'gene_overlap':'exact symbols','feature_link':'exact historical string only; no inferred synonym'},
        'test_family':'No new statistical tests; original CAMP, historical association and locked362 q unchanged',
        'seed':'NA','python':platform.python_version(),'script_sha256':sha(Path(__file__))})
    table(out/'source_manifest.tsv',[{'path':p.relative_to(ROOT).as_posix(),'version':'sha256',
        'sha256':h,'sha256_lf':hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        'hash_scope':'sha256=local exact bytes; sha256_lf=Git text LF form; no source values changed'} for p,h in before.items()])
    dump(out/'validation.json',{'status':'PASS','changed_labels':changed,'historical_sources_unchanged':True,
        'historical_gene_intersection':0,'features':42,'relations':362,'relation_evidence_rows':1448,
        'not_verified':['patient independence','CAMP effect contrast orientation','external cohort independence',
                        'source pairing or preprocessing','new external matrix validation','chemical reidentification']})
    (out/'README_CN.md').write_text('''# PDAC 历史证据解释复核

## 本轮问题
落实对 c508f94 的审核，新增解释复核层，不改历史统计。当前基线 453025c 已在后续提交中交付完整映射及首轮数值；本批不重跑患者矩阵、不重新开展全包审计。

## 输入与范围
已公开历史单细胞14基因/20关系、CPTAC13基因、42条映射账本、锁定362关系及首轮汇总。全部来源见 source_manifest.tsv。逐患者数据仍留服务器。

## 实际结果
- GGT5、MMP14、MGLL 新 reviewed_label 为 DIRECTION_CONCORDANT_NOT_FDR_SIGNIFICANT；原标签、均值/中位数效应和 q 原样留档。两队列所列比较均不显著，不是“两队列显著复现”，也不是证明无效应。
- 单细胞14与CPTAC13基因交集为0。1448行关系证据表覆盖362关系×4证据类型；表达覆盖与实际关系验证分列，不累计所谓双重支持。
- CPTAC的PNP历史假设为脱氧核苷轴；旧鸟苷—PNP关联不能据基因同名升级为外部复现。CPTAC范围表保留历史axis/relation，避免把模块内每个代谢物自动配给每个基因。
- 42条可用性/身份复核全部保留。dCMP9/27、2′-O-methylguanosine7/27、NAD+5/27的作者填补前可用数，与处理后27/27分开呈现。修饰鸟苷身份暂挂，不强制配特异基因。
- source_readiness 在读取RNA后即检查全体基因标签完整唯一；来源哈希不全、不一致均显式失败。当前版本此前已有哈希失败门，本批补齐可测试防护和负例。

## 新手解释
细胞来源偏好、细胞类型内疾病差异、基因表达、同一代谢物—基因关联是不同证据。表中 exact string link 只是历史名称相同，不是化学身份复核或代谢物实测验证。保留的其他历史标签也不代表本批重新认证。

## 限制/反证
首轮统计在135c79b先锁池、fbf21a1再交付数值：主分析使用作者处理值，最低有效n=8，没有为主分析设置填补前可用率门槛；可用性分析另用填补前可用子集，仍要求成对有效n>=8及有变异。不能把这种探索性规则描述为严格观测主分析。NAD+仅5个填补前可用标本，其子集不可评估；主分析不能替代稳健性支持。未来如采用更严格缺失率门槛，应另建预设版本、单独检验族，不根据已经看到的显著性修改v1。
原CAMP q、旧关联q、v1扩展q独立保留。患者独立性、CAMP原始对比方向、队列独立性及来源配对/预处理均未在本批验证；历史tumor_lower不能代替方向核查。无新的GEO原站核验，也无新增外部矩阵验证。

## 当前决定
01冻结整理和02映射批次完成；03探索性数值首轮、04可用性/留一法批次完成，但03/04总阶段仍PARTIAL（单位/协变量等未闭环）。05功能、06外部、07整合仍PARTIAL。DONE仅指索引中限定范围，本批标签与范围复核DONE不等于外部验证DONE。
不得将对c508f94的旧阶段判断覆盖后来真实完成的批次，也不得用后续数值完成掩盖上述限制。此次接受审核中的解释修正和防护建议，历史产物全部保持原样。

## 下一步
先核实患者单位和原对比编码，再决定确认性关联设计；继续全候选功能取证及逐关系外部覆盖，检索命中和基因同名不计作关系复现。患者显著性不是功能筛选唯一门槛。

## 复现命令
`python code/pdac/review_historical_v1.py --run-id <NEW_UTC_RUN_ID>`；新目录必须不存在。测试：`python -m unittest discover -s code/pdac -p "test_*.py"`。本脚本只读取公开汇总，不计算新P/q。
''',encoding='utf-8')
    lock.unlink();print(out.relative_to(ROOT).as_posix())

if __name__=='__main__':main()
