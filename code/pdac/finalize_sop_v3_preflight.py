"""Record completed reuse/mapping and an explicit server access block."""
from pathlib import Path
import json,hashlib,subprocess,tarfile,shutil
import pandas as pd
from prepare_sop_v3 import REPO,ROOT,RUNS,VERSION,SOP_SHA,TEMPLATES,dest,read,emit,write_json,digest,OLD_PAT,OLD_HIST,OLD_SC

def main():
    spec=json.loads((REPO/'docs/PDAC/PDAC_SOP_V3_ANALYSIS_SPEC.json').read_text())
    summary=json.loads((dest('02_MAPPING')/'summary.json').read_text())
    pool=read(dest('02_MAPPING')/'gene_pool_history_union.tsv')
    resolved=read(OLD_PAT/'gene_resolution.tsv').set_index('gene')
    rows=[]
    for _,r in pool.iterrows():
        exists=r.gene in resolved.index;lab=resolved.loc[r.gene,'rna_label'] if exists else None
        rows.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','stable_gene_id':r.stable_gene_id,'gene':r.gene,'measured_feature_id':lab,'mapping_policy':'Historical exact canonical or unique reviewed alias;not freshly re-executed','mapping_version':VERSION,'status':'DONE' if exists and pd.notna(lab) else 'NOT_EVALUABLE' if exists else 'ACCESS_BLOCKED','reason':resolved.loc[r.gene,'reason'] if exists else 'Additional history-union gene needs server matrix labels','source_id':OLD_PAT.relative_to(REPO).as_posix()})
    emit('03_PATIENT','RNA_identity_coverage.tsv',pd.DataFrame(rows))
    rel=read(dest('02_MAPPING')/'direct_relations.tsv');asoc=read(dest('03_PATIENT')/'tumor_association.tsv').set_index('relation_id');rows=[]
    for _,r in rel.iterrows():
        a=asoc.loc[r.relation_id];rows.append({'cancer':'PDAC','cohort':r.cohort,'relation_id':r.relation_id,'metabolite_key':r.metabolite_key,'stable_gene_id':r.stable_gene_id,'gene':r.gene,'metabolite_identity_status':r.identity_status,'RNA_identity_status':'MATCHED' if pd.notna(a.rna_label) else 'UNRESOLVED','n_linked':21,'n_complete':a.n,'status':a.status,'reason':a.reason,'source_id':OLD_PAT.relative_to(REPO).as_posix()})
    emit('03_PATIENT','relation_coverage.tsv',pd.DataFrame(rows))
    registry=[]
    for c in ['GSE263733','GSE278688','GSE242230']:registry.append({'cancer':'PDAC','cohort':c,'study_accession':c,'release':'Inherited hash-locked files from PDAC_sc_paired51_v2','matrix_layer':'raw_counts','expression_scale':'PLANNED_cellwise_log1p10k_then_equal_source_labels','annotation_origin':'Author annotations with inherited frozen broad mapping','donor_identity_level':'Author sample labels' if c=='GSE242230' else 'Author patient labels','disease_scope':'retained primary tumor subset','treatment_available':'NOT_REVIEWED_THIS_BATCH','subtype_available':'NOT_REVIEWED_THIS_BATCH','independence_group':c,'source_id':OLD_SC.relative_to(REPO).as_posix(),'status':'ACCESS_BLOCKED','reason':'Historical matrix available on server;new cellwise method not run;cross-study protected clinical identity not confirmed'})
    emit('06_EXTERNAL','sc_dataset_registry.tsv',pd.DataFrame(registry))
    # Evidence indices are per-gene; gene overlap never becomes relationship validation.
    hist={}
    for file in ['historical_cell_source_gene_summary.tsv','historical_cptac_effects.tsv']:
        df=read(OLD_HIST/file)
        for g in df.gene.unique():hist.setdefault(g,[]).append((OLD_HIST/file).relative_to(REPO).as_posix())
    genes=read(dest('07_INTEGRATION')/'candidate_genes_integrated.tsv')
    genes['external_source_ref']=genes.gene.map(lambda g:';'.join(hist.get(g,[])) or 'NO_LINKED_LEGACY_EXTERNAL_GENE_RESULT')
    genes['historical_RNA_source_ref']=genes.gene.map(lambda g:(OLD_PAT/'paired_RNA.tsv').relative_to(REPO).as_posix() if g in resolved.index else 'NOT_AVAILABLE_IN_PRIOR538_RNA')
    oldsource=set(read(OLD_SC/'cross_cohort_source.tsv').gene)
    genes['historical_sc_source_ref']=genes.gene.map(lambda g:(OLD_SC/'cross_cohort_source.tsv').relative_to(REPO).as_posix() if g in oldsource else 'NOT_AVAILABLE_IN_PRIOR538_SOURCE')
    emit('07_INTEGRATION','candidate_genes_integrated.tsv',genes)
    catalog=read(TEMPLATES/'table_catalog.tsv');delivery=[]
    for _,r in catalog.iterrows():
        files=[dest(s)/r.file for s in RUNS if (dest(s)/r.file).exists()]
        status='ACCESS_BLOCKED' if r.data_scope=='SERVER_PRIVATE' or not files else 'PARTIAL' if r.file.startswith(('paired_RNA','RNA_P005','sc_','candidate_')) else 'DONE'
        delivery.append({**r.to_dict(),'delivery_status':status,'result_paths':';'.join(f.relative_to(REPO).as_posix() for f in files),'reason':'New matrix-dependent calculation not performed' if status=='ACCESS_BLOCKED' else 'Empty RNA P005 table means not run,not zero significant' if r.file=='RNA_P005_view.tsv' else 'Inherited outputs and explicit pending status;see individual table'})
    emit('07_INTEGRATION','SOP_DELIVERY_STATUS.tsv',pd.DataFrame(delivery))
    # Prepare a code/config-only transfer bundle. Contains no patient values or IDs.
    bundle=REPO/'runtime/pdac_sop_v3/server_bundle';bundle.mkdir(parents=True,exist_ok=True);src=bundle/'source';src.mkdir(exist_ok=True)
    for file in ['panel.json','gene_aliases.json']:shutil.copyfile(dest('02_MAPPING')/file,src/file)
    shutil.copyfile(REPO/'docs/PDAC/PDAC_SOP_V3_ANALYSIS_SPEC.json',src/'analysis_spec.json');shutil.copytree(TEMPLATES,src/'templates',dirs_exist_ok=True);shutil.copyfile(OLD_SC/'source_manifest.tsv',src/'historical_manifest.tsv')
    code=['run_sop_v3_server.py','sop_v3_math.py','sc_source_three_cohorts_v1.py','internal_paired51_v2.py','patient_first_round.py','paired_metabolites_v1.py']
    for name in code:shutil.copyfile(REPO/'code/pdac'/name,bundle/name)
    with tarfile.open(bundle.parent/'server_code_and_config.tar','w') as t:
        for p in bundle.rglob('*'):
            if p.is_file():t.add(p,arcname=p.relative_to(bundle))
    precommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    for stage in RUNS:
        out=dest(stage);write_json(out/'execution_provenance.json',{'base_code_commit':precommit,'SOP_commit':SOP_SHA,'new_code_sha256':{n:digest(REPO/'code/pdac'/n) for n in code+['prepare_sop_v3.py','verify_sop_v3_local.py','test_sop_v3_math.py']},'server_connectivity':'172.22.148.165 TCP22 timed out on2026-09-22','server_run_started':False,'local_math_unit_tests':'5/5 PASS','independent_BH_and_aggregate_verification':'PASS;see independent_validation.json','scope':'Public aggregate reuse and versioned mapping;not server pipeline regression'})
        emit(stage,'checksums.tsv',pd.DataFrame([{'file':p.name,'sha256':digest(p)} for p in sorted(out.iterdir()) if p.is_file() and p.name!='checksums.tsv']))
    index=read(REPO/'coordination/stages/PDAC.tsv');index=index[~index.run_id.isin(RUNS.values())]
    for stage,run in RUNS.items():
        state=summary['stage_status'][stage]
        index=pd.concat([index,pd.DataFrame([{'cancer':'PDAC','stage_id':stage,'run_id':run,'analysis_version':VERSION,'status':state['status'],'scope':{'01_CAMP':'307 full metabolites;51 P005;193 availability evaluable at n>=8','02_MAPPING':'357direct;358conditional;250current;687historyunion','03_PATIENT':'Association reused,new evaluableBH;paired t ACCESS_BLOCKED','06_EXTERNAL':'687genes3cohorts planned;new cellwise source ACCESS_BLOCKED','07_INTEGRATION':'715relation687gene scaffold;RNA/source pending'}[stage],'result_path':dest(stage).relative_to(REPO).as_posix(),'code_path':'code/pdac/prepare_sop_v3.py','git_branch':'analysis/pdac-initial','reason':state['reason'],'next_action':'Restore server network;execute frozen v3 internal/source;validate and append'}])],ignore_index=True)
    index.sort_values(['stage_id','run_id']).to_csv(REPO/'coordination/stages/PDAC.tsv',sep='\t',index=False)
    report=f'''# PDAC 当前进度：按 GitHub SOP v3 执行，服务器计算待恢复

## 本轮问题

用户要求按GitHub的CAMP发现至单细胞来源规范执行PDAC。本轮采用[固定提交规范](https://github.com/Asukasssss/-/blob/{SOP_SHA}/docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md)，单细胞止于表达来源，不新增机制分析。

## 输入与范围

作者明确配对11对，关联使用21个作者编号互不重复的肿瘤单位；6个缺编号肿瘤不默认独立。全307代谢物主P<0.05工作池51项，配对q<0.05仍为0。现有357直接与358条件关系复用；按SOP当前DIRECT池为250基因，历史并集为687基因。此前538是配对51项直接加条件池，不是新规范当前DIRECT池。

## 实际结果

|部分|状态|本轮实际完成|
|---|---|---|
|全量代谢物与原有值敏感性|PARTIAL|保留原效应/P/CI，最低8对与新可评估BH；307主可评估，193敏感性可评估。缺每侧可用性边际比例|
|映射与池定义|DONE（本批整理）|51项去向、357直接、358条件；250当前基因、687历史并集；稳定UniProtID唯一|
|肿瘤内部关联|DONE（统计复用批）|四个直接/条件×主/可用性集合重新BH，效应/P和种子保持不变|
|新配对RNA|ACCESS_BLOCKED|已冻结配对t、4000整对bootstrap；未实际执行。旧Wilcoxon全表及180项P<0.05展示另存历史|
|新单细胞来源|ACCESS_BLOCKED|三队列×687基因拟按逐细胞log1p10k、供者等权、80%描述标签运行；没有把旧pseudobulk结果改名|
|两张最终比较表|PARTIAL|715关系和687基因保留，新增RNA/来源明确待算；没有A/B/C/D总分类|

GABA—ABAT相同输入关联rho仍为0.84866；新直接主q=0.0354，可用性rho=0.84412、q=0.0342。改变来自357→354、357→342的BH有效检验数，旧q=0.0357保留；这不是新验证。GABA配对q=0.2578仍未通过FDR。

## 新手解释

当前250基因是直接关系集合；其余条件和历史基因继续保留，不按P值删去。新RNA展示表目前没有结果行，意思是未运行，不是“0项P<0.05”。本地五项方法单元测试及六个BH族的独立复核通过，不等于服务器流程回归。

## 限制/反证

2026-09-22访问172.22.148.165:22多次超时，未到认证阶段。原矩阵和真实连接表按规定留服务器，无法在本机补算。临床身份、分装层级、临床协变量和跨研究患者独立性仍有未确认项。旧303/42、旧362关联及v2配对RNA/来源不覆盖。

## 当前决定

已完成部分逐批上传；新规范整流程尚未完成。保持PR #2草稿、未合并。外部代谢关系验证及功能深入不属于本次完成前置。

## 下一步

恢复可访问服务器的内网/VPN后，执行冻结的新目录代码，补齐新配对RNA、单细胞来源、必备供者等权点图/全候选热图，再更新整合表。现成UMAP在SOP为可选，本轮不重新嵌入。

## 复现命令

本地：`python code/pdac/prepare_sop_v3.py`，`python code/pdac/verify_sop_v3_local.py`，`python code/pdac/finalize_sop_v3_preflight.py`。
服务器：新的独占PDAC/B目录执行`run_sop_v3_server.py --mode internal`或`--mode source`，各自目录传入`--code-commit`。仅代码和注释包可传输，患者/细胞级数据不得回传。

交付入口：

- [逐表完成状态](../../results/PDAC/07_INTEGRATION/{RUNS['07_INTEGRATION']}/SOP_DELIVERY_STATUS.tsv)
- [关系级表（部分完成）](../../results/PDAC/07_INTEGRATION/{RUNS['07_INTEGRATION']}/candidate_relations_integrated.tsv)
- [基因级表（部分完成）](../../results/PDAC/07_INTEGRATION/{RUNS['07_INTEGRATION']}/candidate_genes_integrated.tsv)
- [当前基因池](../../results/PDAC/02_MAPPING/{RUNS['02_MAPPING']}/gene_pool_current.tsv)
- [原v2结果解释，保留历史口径](../../results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2/FINDINGS_CN.md)
'''
    (REPO/'reports/PDAC/CURRENT_STATUS_CN.md').write_text(report,encoding='utf8')
    print('Recorded v3 partial delivery and server block;code/config bundle prepared;no server launch')

if __name__=='__main__':main()
