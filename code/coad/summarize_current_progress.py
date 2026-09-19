"""Join the fresh enzyme/transport evidence without treating repeated sources as tests."""
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'reports/COAD'


def read(name):
    return json.loads((BASE/name).read_text(encoding='utf-8'))


def main():
    identity=read('fresh_mapping_v0_1/identity_review.json')
    enzymes=read('human_reaction_check_v0_1/human_reaction_pairs.json')
    transport=read('transport_check_v0_1/transport_pairs.json')
    source=ROOT/'reference/camp/cancer_effects.tsv'
    source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    original=[r for r in csv.DictReader(source.open(encoding='utf-8-sig'),delimiter='\t') if r['cancer']=='COAD']
    expected={r['metabolite_key'] for r in original if float(r['effect_fdr'])<.05}
    assert len(original)==159 and len(expected)==73 and expected=={r['metabolite_key'] for r in identity}
    catalog={}
    for r in enzymes+transport:
        key=(r['dataset'],r['feature_name'],r['metabolite_key'],r['human_gene_id'])
        if key not in catalog:
            catalog[key]={k:r[k] for k in ['dataset','feature_name','metabolite_key','human_gene_id','gene_symbol','identity_review_status']}
            catalog[key].update({'enzyme_annotation_support':False,'transport_annotation_support':False,
                                 'enzyme_experimental_annotation':False,'transport_experimental_annotation':False,
                                 'human_enzyme_check':'NOT_IN_KEGG_CANDIDATE_SEARCH',
                                 'source_record_files':[],'compound_hypotheses':[],
                                 'patient_association':'NOT_RUN','coad_functional_evidence':'NOT_RUN',
                                 'complex_role_review':'NOT_COMPLETED_PER_GENE','ready_for_patient_testing':False})
        item=catalog[key]
        item['compound_hypotheses']=sorted(set(item['compound_hypotheses'])|set(r['compound_hypotheses']))
        if 'matched_support' in r:
            item['enzyme_annotation_support']=bool(r['matched_support'])
            item['enzyme_experimental_annotation']=any(m['has_experimental_annotation'] for m in r['matched_support'])
            item['human_enzyme_check']=r['human_reaction_check']
            item['source_record_files'].append('human_reaction_check_v0_1/human_reaction_pairs.json')
        else:
            item['transport_annotation_support']=True
            item['transport_experimental_annotation']=r['has_experimental_annotation']
            item['source_record_files'].append('transport_check_v0_1/transport_pairs.json')
    catalog=list(catalog.values())
    for item in catalog:
        item['database_reaction_supported']=item['enzyme_annotation_support'] or item['transport_annotation_support']
        item['no_current_identity_hold']=not item['identity_review_status'].startswith('HOLD_')
        item['disposition']=('IDENTITY_HOLD' if not item['no_current_identity_hold']
                             else 'REACTION_ANNOTATION_SUPPORTED_PROVISIONAL' if item['database_reaction_supported']
                             else 'UNIPROT_RECORD_NOT_EVALUABLE' if item['human_enzyme_check']=='NO_REVIEWED_HUMAN_RECORD_LINKED'
                             else 'REACTION_SPECIFICITY_PENDING')
    supported=[r for r in catalog if r['database_reaction_supported']]
    pending_free=[r for r in supported if r['no_current_identity_hold']]
    features=[]
    for ident in identity:
        related=[r for r in catalog if r['metabolite_key']==ident['metabolite_key']]
        supported_rows=[r for r in related if r['database_reaction_supported']]
        features.append({k:ident[k] for k in ['dataset','feature_name','metabolite_key','kegg_id','hmdb_id','hedges_g','effect_fdr','max_input_raw_missing_rate','identity_review_status','identity_note']} | {
            'candidate_gene_count':len(related),'reaction_supported_gene_count':len(supported_rows),
            'reaction_supported_genes':sorted({r['gene_symbol'] or r['human_gene_id'] for r in supported_rows}),
            'enzyme_supported_count':sum(r['enzyme_annotation_support'] for r in related),
            'transport_supported_count':sum(r['transport_annotation_support'] for r in related),
            'patient_association':'NOT_RUN'})
    summary={'version':'COAD_B_BIOCHEMICAL_DRAFT_20260919','created_utc':datetime.now(timezone.utc).isoformat(),
             'source_sha256':source_hash,'all_coad_effects':159,'significant_features_reviewed':73,
             'candidate_pairs_total':len(catalog),'reaction_supported_pairs_including_identity_holds':len(supported),
             'reaction_supported_pairs_without_current_identity_hold':len(pending_free),
             'supported_features_without_current_identity_hold':len({r['metabolite_key'] for r in pending_free}),
             'supported_genes_without_current_identity_hold':len({r['human_gene_id'] for r in pending_free}),
             'features_without_reaction_support':sum(r['reaction_supported_gene_count']==0 for r in features),
             'disposition_counts':dict(Counter(r['disposition'] for r in catalog)),
             'ready_for_patient_testing':0,'old_mappings_used':False,'original_stats_changed':False,
             'ssh_connection':'ACCESS_BLOCKED_HOST_ALIAS_UNRESOLVED',
             'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    assert len(catalog)==len({(r['dataset'],r['feature_name'],r['metabolite_key'],r['human_gene_id']) for r in catalog})
    assert any(r['feature_name']=='taurine' and r['gene_symbol']=='SLC6A6' and r['transport_annotation_support'] for r in catalog)
    assert any(r['feature_name']=='beta-alanine' and r['gene_symbol']=='UPB1' and r['enzyme_annotation_support'] for r in catalog)
    assert all(not r['enzyme_annotation_support'] for r in catalog if r['feature_name']=='proline' and r['gene_symbol']=='P4HA1')
    assert all(r['disposition']=='IDENTITY_HOLD' for r in catalog if r['feature_name']=='glutamate')
    for f in features:
        raw=next(r for r in original if r['metabolite_key']==f['metabolite_key'])
        assert all(f[k]==raw[k] for k in ['feature_name','kegg_id','hmdb_id','hedges_g','effect_fdr','max_input_raw_missing_rate'])
    out=BASE/'current_catalog_v0_1'
    out.mkdir(parents=True,exist_ok=True)
    for name,obj in [('candidate_catalog.json',catalog),('feature_catalog.json',features),('summary.json',summary)]:
        (out/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    lines=['# COAD 账号 B 当前进展', '',
           '本轮从 159 条冻结效应记录出发，对其中原 q<0.05 的 73 条重新建立数据库证据。未复用原基因映射、BRCA 候选或其他癌种分析。冻结统计不变，86 条非显著背景保留。', '',
           '| 已完成环节 | 本轮结果 |','|---|---|',
           '| 名称与标识初查 | 73 条全部记录去向；63 条原 KEGG 名称相容，另 10 条有身份、暴露背景或新编号核查事项 |',
           '| KEGG 反应候选发现 | 60 个特征，693 条特征—人类基因候选关系；尚不等于直接底物证据 |',
           '| 人类催化反应注释核对 | 482 条候选匹配到 reviewed 人类蛋白的具体 Rhea 反应；其中 58 条存在身份暂挂 |',
           '| 结构化转运反应检索 | 43 个特征，281 条关系；其中 31 条存在身份暂挂 |',
           f"| 两类证据去重合并 | {len(supported)} 条有反应注释支持；其中 {len(pending_free)} 条无当前特定身份暂挂，涉及 {summary['supported_features_without_current_identity_hold']} 个特征、{summary['supported_genes_without_current_identity_hold']} 个基因 |", '',
           '**这是一版有来源的生化候选目录，不是已完成患者关联的结果或最终靶点排名。没有身份暂挂也不等于完成原始鉴定。所有关系的患者分析均为 NOT_RUN。**', '',
           '# 会改变下一步分析的发现', '',
           '1. glutamate 的 KEGG 是 L 型，原 HMDB 是 D 型：保留原记录与两个假设，暂不把相关基因放入正式检验家族。',
           '2. ribulose/xylulose 5-phosphate 是合并特征：不能拆成两个独立测量。',
           '3. cysteine-glutathione disulfide 原 KEGG 栏是反应编号 R00900；Mucate 的原 C01807 本次未返回；另为 ophthalmate 查到拟议 C21016。这些提议均未覆盖原注释。',
           '4. Disulfiram 等五条需先核对原始注释或暴露背景，不能因药物有靶点就直接连入代谢关系。',
           '5. 人类反应核对保留了 beta-alanine—UPB1 等有直接反应注释的组合；游离 proline—P4HA1 没有匹配到相同反应，留在待审区，不能把蛋白残基修饰当作游离代谢物转换。',
           '6. 转运检索补充了 taurine—SLC6A6、creatinine 对应转运关系等，避免仅凭酶反应缺失就淘汰特征。', '',
           '# 全部 73 条的当前覆盖', '',
           '| 原特征 | 原 g | 原 q | 催化注释基因数 | 转运注释基因数 | 身份状态 |', '|---|---:|---:|---:|---:|---|']
    for f in sorted(features,key=lambda r:float(r['effect_fdr'])):
        lines.append(f"| {f['feature_name']} | {float(f['hedges_g']):.3f} | {float(f['effect_fdr']):.4g} | {f['enzyme_supported_count']} | {f['transport_supported_count']} | {f['identity_review_status']} |")
    lines+=['', '# 未完成与接续条件', '',
            '- 作者源注释、样本对应、实际 RNA 覆盖与患者独立单位尚未核查；本机 server165 SSH 别名无法解析，等待连接配置。',
            '- 部分 HMDB 页面访问被拒，未完成全面交叉核验；不以网页不返回作为化合物不存在的证据。',
            '- 此映射依赖当前结构化数据库覆盖。无命中不是阴性；KO 未覆盖的人类酶、其他区室表述的转运、复杂底物范围和复合体角色仍需补查。',
            '- 注释中的实验来源不等于我们已阅读论文，也不等于 COAD 功能证明。尚未执行 RNA、患者关联、依赖性或外部验证。',
            '- 接下来先核对源注释和数据覆盖，再锁定版本化关系与检验范围。候选功能取证与患者分析可并行，不把相关显著当唯一门槛。', '',
            '# 文件与来源', '',
            '- 当前特征目录：current_catalog_v0_1/feature_catalog.json；完整候选（含暂挂和待审）：current_catalog_v0_1/candidate_catalog.json。',
            '- 身份与 KEGG 反应明细：fresh_mapping_v0_1/；人类蛋白逐条核对：human_reaction_check_v0_1/；转运反应：transport_check_v0_1/。每条证据含 URL，来源清单含访问时间和哈希。',
            '- 脚本：code/coad/ 下对应脚本。新内容目前仅本地，未提交或推送。',
            '- 原数据 SHA-256：`'+source_hash+'`。',
            '- 核查通过：原效应/原 q/原标识逐值未改；73 条全覆盖；唯一特征—基因键无重复；UPB1 与 SLC6A6 阳性关系、P4HA1 底物区分和 glutamate 暂挂规则通过。',
            '', '公开来源：[KEGG](https://www.kegg.jp/kegg/rest/keggapi.html)、[HMDB 谷氨酸条目](https://hmdb.ca/metabolites/HMDB0003339)、[UniProt](https://www.uniprot.org/help/api)、[Rhea](https://www.rhea-db.org/help/download)。']
    (BASE/'CURRENT_STATUS_CN.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    assert hashlib.sha256(source.read_bytes()).hexdigest()==source_hash
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
