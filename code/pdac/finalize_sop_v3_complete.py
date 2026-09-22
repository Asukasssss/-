"""Finish the actual source/integration delivery after figure checks."""
from pathlib import Path
import json,hashlib
import pandas as pd
from record_sop_v3_resume import REPO,ROOT,SOURCE,INTERNAL,DISCOVERY,INTEGRATION,MAP,read,write,write_json,checksum,stage

def main():
    s=json.loads((INTEGRATION/'summary.json').read_text());pv=json.loads((SOURCE/'figures/visual_validation.json').read_text());assert pv['status']=='PASS'
    v=json.loads((SOURCE/'validation.json').read_text());v['figure_status']='PASS';v['figure_validation']=pv;write_json(SOURCE/'validation.json',v)
    sp=json.loads((SOURCE/'analysis_spec.json').read_text());sp['status']='DONE_SOURCE_AND_REQUIRED_PLOTS';sp['figures']={'dotplots':15,'heatmaps':15,'scope':'All687 alphabetical;no significance filter','UMAP':'OPTIONAL_NOT_SCOPED'};write_json(SOURCE/'analysis_spec.json',sp)
    rm=(SOURCE/'README_CN.md').read_text();rm=rm.replace('生成全687基因点图/热图并完成视觉核查，再接入关系/基因两表；不自动开展机制或细胞通讯。','已生成并核查15页点图与15页热图，覆盖全部687基因。接入关系/基因两表；不自动开展机制或细胞通讯。\n\n[全候选来源图](figures/README_CN.md) · [逐图输入索引](figures/figure_manifest.tsv)');(SOURCE/'README_CN.md').write_text(rm,encoding='utf8',newline='\n')
    stage(SOURCE,'06_EXTERNAL','687genes3cohorts;15dotplot+15heatmap pages;all687 retained','Cellwise log1p10k source summaries independently checked;required plots visually verified','Source-only scope complete;no automatic mechanism expansion')
    # Explicitly index every prescribed table, including private server tables and copies/views.
    templates=REPO/'code/pdac/sop_v3_templates';catalog=read(templates/'table_catalog.tsv');registry=[]
    locs={}
    for folder in [MAP,DISCOVERY,ROOT/'03_PATIENT/20260922T112700Z_sop_v3_patient',INTERNAL,SOURCE,INTEGRATION]:
        for f in folder.glob('*.tsv'):locs.setdefault(f.name,[]).append(f.relative_to(REPO).as_posix())
    if (SOURCE/'sc_celltype_profile_parts.tsv').exists():
        locs['sc_celltype_profiles.tsv']=[(SOURCE/f).relative_to(REPO).as_posix() for f in read(SOURCE/'sc_celltype_profile_parts.tsv').file]
    server='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/PDAC/B/'
    for _,r in catalog.iterrows():
        if r.data_scope=='SERVER_PRIVATE':
            p=server+(SOURCE.name+'/private/GSE*_sc_donor_profiles_private.tsv' if r.file=='sc_donor_profiles_private.tsv' else INTERNAL.name+'/private/'+r.file)
            registry.append({**r.to_dict(),'status':'DONE_SERVER_PRIVATE','result_paths':p,'reason':'Not exported to workstation or GitHub;actual computation created these files'})
        else:
            pp=locs.get(r.file,[]);assert pp,r.file
            registry.append({**r.to_dict(),'status':'DONE','result_paths':';'.join(pp),'reason':'Complete planned rows and individual not-evaluable statuses retained;not mechanism/external validation'})
    write(INTEGRATION/'SOP_DELIVERY_STATUS.tsv',pd.DataFrame(registry))
    result=read(INTEGRATION/'candidate_relations_integrated.tsv');sig=result[result.association_q<.05].sort_values(['mapping_status','association_q']);genes=read(INTEGRATION/'candidate_genes_integrated.tsv');abat=genes[genes.gene=='ABAT'].iloc[0]
    findings=['# PDAC SOP v3 本轮结果','', '## 本轮问题','', '从CAMP配对探索工作池建立直接关系、患者背景和单细胞表达来源比较表。止于来源描述。','', '## 输入与范围','', '11个明确作者配对，21个作者病例键互不重复的肿瘤；250当前直接基因，687历史并集。','', '## 实际结果','', '|关系|映射层|rho|关联q|可用性n|可用性q|RNA配对q|','|---|---|---:|---:|---:|---:|---:|']
    for _,r in sig.iterrows():findings.append(f'|{r.metabolite_name}—{r.gene}|{r.mapping_status}|{r.association_rho:.4f}|{r.association_q:.6f}|{int(r.association_available_n)}|{r.association_available_q:.6f}|{r.RNA_q:.6f}|')
    findings+=['','四条条件关系不因关联显著升级为直接生化关系。GPC—GDPD5和leucine—SLC7A7主/可用性使用相同输入，q改变来自其余检验组成，不是独立重复。RNA q分别属于当前或历史补充基因集合，完整字段见关系表。','','### ABAT 的表达来源','','|研究|最高类别|次高类别|最高类别检出率|最高类别保持率|','|---|---|---|---:|---:|']
    for c in ['GSE263733','GSE278688','GSE242230']:findings.append(f"|{c}|{abat[c+'_top_celltype']}|{abat[c+'_runner_celltype']}|{abat[c+'_top_detection_fraction']:.1%}|{abat[c+'_bootstrap_top_frequency']:.1%}|")
    findings+=['',f"ABAT配对RNA q={abat.RNA_q:.6f}，11对中10对降低。GABA配对P=0.032227、q=0.257832，仍是探索入口；不能因后续相关较强而称配对发现已通过FDR。",'', '## 新手解释','', f"三队列同一明确最高类别{s['all3_same_top_genes']}个基因，其中各队列保持率均≥80%为{s['all3_top_bootstrap_ge080_genes']}个。687个都在分析中，没有按RNA或相关显著性取交集。",'', '## 限制/反证','', '作者来源标签不自动证明独立患者；GSE242230按25个作者样本标签汇总。跨研究临床身份去重、同一分装、临床协变量、恶性CNV及机制均没有通过这次分析得到确认。表达最多的类别不是唯一作用细胞，也不证明代谢物在该细胞生成或消耗。', '', '## 当前决定','', '规范内的数值、来源图及比较表交付完成；保留测量缺项和限制，不设必须产出的A/B名单，不覆盖任何历史P/q。','', '## 下一步','', '使用完整比较表讨论候选。额外机制或外部代谢关系验证须作为后续明确任务，不由本轮自动扩展。','', '## 复现命令','', '`python code/pdac/integrate_sop_v3_complete.py`；`python code/pdac/finalize_sop_v3_complete.py`。统计与图形的来源哈希见各批次manifest/checksums。']
    (INTEGRATION/'FINDINGS_CN.md').write_text('\n'.join(findings)+'\n',encoding='utf8',newline='\n')
    report=f'''# PDAC 当前进度：发现至表达来源SOP本轮交付完成

## 本轮问题

按GitHub固定规范6dd0a73执行PDAC，主线为身份连接→配对代谢物→直接映射→内部相关与配对RNA→全候选表达来源→关系/基因比较表。单细胞止于来源描述。

## 输入与范围

11个作者明确配对，21个作者病例键互不重复的肿瘤单位。全307项代谢物、51项探索工作池；357直接关系、358条件关系。当前DIRECT池250基因，含条件及历史的并集687基因。

## 实际结果

|环节|实际结果|
|---|---|
|配对代谢物|307项可评估；51项P<0.05、0项q<0.05；双方原有值敏感性193项可评估|
|内部关联|21作者单位；直接主354项可评估、1条q<0.05；条件主350项可评估、4条q<0.05。效应/P复用，新有效集合BH另存|
|配对RNA当前池|247/250可评估；91项P<0.05，40项q<0.05|
|RNA条件/历史补充|429/437可评估；151项P<0.05，0项q<0.05，独立补充检验族|
|三队列表达来源|全部687基因；细胞逐个log1p10k后供者/样本标签等权。{s['all3_same_top_genes']}基因同最高类别；{s['all3_top_bootstrap_ge080_genes']}基因同最高且各队列保持率≥80%|
|来源图|15页点图、15页热图，全部687基因按字母分页，未按显著性筛图|
|完整整合|715关系表、687基因表、7项未映射特征单列；原42池362关系ID另索引；无综合总分|

676行RNA配对t、t区间和bootstrap重放、两个RNA BH族通过独立核对；来源全体供者汇总和排名通过独立核对；完整关系键连接与图表覆盖验证通过。没有声称独立重读全部原始单细胞计数或重放全部SCbootstrap。

## 新手解释

P<0.05代谢物是探索入口，51项均未通过配对FDR。RNA差异、肿瘤内部相关和单细胞来源回答不同问题，不以全部显著才保留。80%只是来源排名的描述性重采样保持率，不是80%患者符合，也不是机制证据。

## 限制/反证

临床身份重新认证、同一分装、跨研究患者去重、临床协变量、恶性CNV、外部matched-omics和功能干预均未因本次交付完成。GSE242230按作者25样本标签分析，不能提升为25已核实独立患者。原303/42及旧关联/RNA/来源结果完全保留。

## 当前决定

本规范范围的计算、来源图和比较表完成，各交付批次DONE；03/04/05/06/07总研究阶段仍有上述范围外缺项，不宣称PDAC所有研究已完成。PR #2保持草稿、未合并。

## 下一步

以完整比较表进行选题讨论；不自动追加机制任务。

## 复现命令

代码位于code/pdac/run_sop_v3_server.py、verify_sop_v3_server.py、record_sop_v3_resume.py、plot_sop_v3_sources.py与integrate_sop_v3_complete.py。真实患者、细胞级数据和完整映射表继续留服务器。

交付入口：

- [结果解释](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/FINDINGS_CN.md)
- [关系比较表](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/candidate_relations_integrated.tsv)
- [基因比较表](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/candidate_genes_integrated.tsv)
- [RNA P<0.05展示](../../results/PDAC/03_PATIENT/{INTERNAL.name}/RNA_P005_view.tsv)
- [全候选来源图](../../results/PDAC/06_EXTERNAL/{SOURCE.name}/figures/README_CN.md)
- [31项规范表交付索引](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/SOP_DELIVERY_STATUS.tsv)
'''
    (REPO/'reports/PDAC/CURRENT_STATUS_CN.md').write_text(report,encoding='utf8',newline='\n');checksum(SOURCE);checksum(INTEGRATION);print(json.dumps(s))
if __name__=='__main__':main()
