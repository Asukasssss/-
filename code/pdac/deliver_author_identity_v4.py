"""Record the new author-identity source layer without modifying historical results."""
from pathlib import Path
import argparse,json,hashlib,subprocess
import numpy as np
import pandas as pd
import record_sop_v3_resume as records
from record_sop_v3_resume import REPO,ROOT,MAP,INTERNAL,DISCOVERY,read,write,write_json,checksum,readme,stage
import integrate_sop_v3_complete as integrate
from export_sop_v3_workbook_data import records as json_records
SOURCE=ROOT/'06_EXTERNAL/20260922T131800Z_author_identity_v4'
INTEGRATION=ROOT/'07_INTEGRATION/20260922T131900Z_author_identity_v4'
PRE_SOURCE=records.SOURCE;PRE_INTEGRATION=records.INTEGRATION
VERSION='PDAC_author_cell_identity_v4'
OUT=REPO/'outputs/pdac-author-identity-v4-20260922'

def record_source():
    assert json.loads((SOURCE/'independent_validation.json').read_text())['status']=='PASS'
    sp=json.loads((SOURCE/'analysis_spec.json').read_text());sp['status']='DONE_SOURCE_NUMERICS';write_json(SOURCE/'analysis_spec.json',sp)
    registry=read(PRE_SOURCE/'sc_dataset_registry.tsv');registry['identity_resolution']=registry.cohort.map({'GSE263733':'Ductal unresolved;no available exact-cell malignant calls','GSE278688':'Ductal unresolved;no available exact-cell malignant calls','GSE242230':'10783 author malignant;350 author normal;fine label overrides broad'});write(SOURCE/'sc_dataset_registry.tsv',registry)
    readme(SOURCE,'PDAC 作者恶性身份接回 v4','同一687基因和148777细胞；只改变有依据的身份划分，不按表达结果推断恶性。','GSE242230旧Malignant大类11133细胞拆为10783作者恶性（5660 Classical、5123 Basal）和350作者正常上皮。该队列全687来源统计重新计算；另两队列只将导管显示为身份未定并复用数值。','细分Normal Epithelial与宽类Malignant冲突，优先使用作者更具体身份。另两队列论文恶性分析不能代替逐细胞身份连接。作者注释不等于本轮独立CNV验证。','更新全候选图与比较表；不扩展机制。','`python code/pdac/run_author_identity_v4_server.py`（服务器新目录）；本地`python code/pdac/deliver_author_identity_v4.py --mode record`。')
    evidence=[
      {'cohort':'GSE242230','source_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC10587349/','finding':'Authors used CopyKAT and normal epithelial markers;local author fine labels include Normal Epithelial within broad Malignant','cell_level_status':'DONE','analysis_action':'Use fine identity on exactly matching selected cell barcodes;350 normal split from10783 malignant'},
      {'cohort':'GSE263733','source_url':'https://github.com/CompbioLabUnist/PDAC-scRNA-seq','finding':'Repository Source_code includes CopyKAT;tree contains README and code only;GEO metadata has4 columns and no malignant/CNV/subtype field','cell_level_status':'NEEDS_REVIEW','analysis_action':'Ductal unresolved;paper-level malignant claims not assigned to individual cells'},
      {'cohort':'GSE278688','source_url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278688','finding':'GEO metadata has barcode,tissue,patients,all_celltype;no malignant/CNV/subtype field. Author study describes malignant subpopulations;publisher fulltext retrieval403 in this check','cell_level_status':'NEEDS_REVIEW','analysis_action':'Ductal unresolved;no claim that authors never identified malignant cells'}]
    write(SOURCE/'author_identity_evidence.tsv',pd.DataFrame(evidence))
    write_json(SOURCE/'validation.json',{'status':'PASS','numeric_verification':json.loads((SOURCE/'independent_validation.json').read_text()),'cell_join':json.loads((SOURCE/'cell_join_validation.json').read_text()),'all687_retained':True,'original_matrices_not_exported':True,'figure_status':'NOT_RUN','remaining_identity_gaps':['GSE263733 exact-cell malignant labels','GSE278688 exact-cell malignant labels'],'new_CNV_not_run':True})
    import partition_sop_v3_profiles as partition
    partition.SOURCE=SOURCE;partition.main()
    stage(SOURCE,'06_EXTERNAL','Author identity correction;687genes;350normal split;2cohorts unresolved','Numbers independently checked;figures pending','Complete plots and comparisons',status='PARTIAL');checksum(SOURCE)

def integrate_results():
    integrate.SOURCE=SOURCE;integrate.INTEGRATION=INTEGRATION;integrate.VERSION=VERSION;records.VERSION=VERSION;integrate.main()
    old=read(PRE_INTEGRATION/'candidate_relations_integrated.tsv').set_index('relation_id');new=read(INTEGRATION/'candidate_relations_integrated.tsv').set_index('relation_id')
    numeric=[c for c in old if c.startswith(('metabolite_','association_','RNA_')) and pd.api.types.is_numeric_dtype(old[c])]
    np.testing.assert_allclose(old.loc[new.index,numeric].to_numpy(float),new[numeric].to_numpy(float),rtol=0,atol=1e-12,equal_nan=True)
    before=read(PRE_SOURCE/'sc_source_stability.tsv');after=read(SOURCE/'sc_source_stability.tsv');change=before.merge(after,on=['cohort','gene','stable_gene_id'],suffixes=('_v3','_v4'),validate='one_to_one')
    cols=['cohort','gene','stable_gene_id','top_celltype_v3','top_celltype_v4','bootstrap_top_frequency_v3','bootstrap_top_frequency_v4','runner_celltype_v3','runner_celltype_v4','status_v3','status_v4']
    change=change[cols];change['top_label_changed']=change.top_celltype_v3!=change.top_celltype_v4
    def lineage(x):return 'Epithelial lineage' if x in ['Epithelial/Ductal','Malignant (author)','Normal epithelial (author)','Ductal (unresolved)'] else x
    change['broad_lineage_top_changed']=[lineage(a)!=lineage(b) for a,b in zip(change.top_celltype_v3,change.top_celltype_v4)]
    write(INTEGRATION/'source_identity_changes.tsv',change)
    genes=read(INTEGRATION/'candidate_genes_integrated.tsv');oldgenes=read(PRE_INTEGRATION/'candidate_genes_integrated.tsv')
    # Add all-gene malignant vs normal expression comparison, never a disease test.
    prof=read(SOURCE/'sc_celltype_profiles.tsv');view=prof[(prof.cohort=='GSE242230')&prof.celltype.isin(['Malignant (author)','Normal epithelial (author)'])].copy();write(INTEGRATION/'author_malignant_normal_gene_profiles.tsv',view)
    summary=json.loads((INTEGRATION/'summary.json').read_text());summary.update(old_all3_top_bootstrap_ge080_genes=int(oldgenes.all3_top_bootstrap_ge080.sum()),old_current_stable=int(oldgenes.loc[oldgenes.in_current_pool,'all3_top_bootstrap_ge080'].sum()),GSE242230_broad_lineage_top_changes=int(change.loc[change.cohort=='GSE242230','broad_lineage_top_changed'].sum()),GSE242230_new_top_counts=genes.GSE242230_top_celltype.value_counts().to_dict(),patient_numeric_fields_preserved=len(numeric))
    write_json(INTEGRATION/'summary.json',summary)
    catalog=read(PRE_INTEGRATION/'SOP_DELIVERY_STATUS.tsv')
    catalog['result_paths']=catalog.result_paths.str.replace(PRE_SOURCE.name,SOURCE.name,regex=False).str.replace(PRE_INTEGRATION.name,INTEGRATION.name,regex=False)
    catalog.loc[catalog.file.str.startswith('sc_'),'reason']='v4 retains author identity;GSE242230 split;other2 ductal identities unresolved;not all3 malignant validation'
    write(INTEGRATION/'SOP_DELIVERY_STATUS.tsv',catalog)
    v=json.loads((INTEGRATION/'validation.json').read_text());v.update(patient_all_numeric_fields_unchanged=True,patient_numeric_field_count=len(numeric),strict_identity_labels=True);write_json(INTEGRATION/'validation.json',v)
    readme(INTEGRATION,'PDAC v4：保留作者恶性身份的比较表','715关系、687基因完整保留，代谢物/RNA/关联统计沿用上一版。',f"三队列严格同首位且各≥80%为{summary['all3_top_bootstrap_ge080_genes']}基因，当前直接池{summary['all3_top_bootstrap_ge080_current_genes']}基因。GSE242230有{summary['GSE242230_broad_lineage_top_changes']}基因的首位发生跨宽谱系变化；纯身份重命名单列。",'原29个上皮/导管稳定基因不能称三队列恶性来源复现。作者恶性、作者正常、身份未定不等价。严格共同类别比较与原宽谱系比较回答不同问题。','完成新图并核查；原v3结果保留为历史，新的来源解释优先用v4。','`python code/pdac/deliver_author_identity_v4.py --mode integrate`');checksum(INTEGRATION)

def workbook_data():
    OUT.mkdir(exist_ok=True,parents=True);old=json.loads((REPO/'outputs/pdac-sop-v3-20260922/workbook_data.json').read_text());rel=read(INTEGRATION/'candidate_relations_integrated.tsv');rel['order']=rel.mapping_status.map({'DIRECT':0,'CONDITIONAL':1});old['relations']=json_records(rel.sort_values(['order','association_q','metabolite_name','gene']).drop(columns='order'));old['genes']=json_records(read(INTEGRATION/'candidate_genes_integrated.tsv').sort_values('gene'));old['summary']=json.loads((INTEGRATION/'summary.json').read_text());old['source_run']=SOURCE.name;old['integration_run']=INTEGRATION.name;old['identity_v4']=True
    old['identity_profiles']=json_records(read(INTEGRATION/'author_malignant_normal_gene_profiles.tsv').sort_values(['gene','celltype']));old['identity_counts']=json_records(read(SOURCE/'author_identity_counts.tsv'))
    (OUT/'workbook_data.json').write_text(json.dumps(old,ensure_ascii=False,allow_nan=False),encoding='utf8');print(OUT)

def finalize():
    v=json.loads((SOURCE/'figures/visual_validation.json').read_text());assert v['status']=='PASS'
    s=json.loads((INTEGRATION/'summary.json').read_text());genes=read(INTEGRATION/'candidate_genes_integrated.tsv');ab=genes[genes.gene=='ABAT'].iloc[0]
    sp=json.loads((SOURCE/'analysis_spec.json').read_text());sp['status']='DONE_RESOLVABLE_IDENTITY_SCOPE';sp['unresolved_cohorts']=['GSE263733','GSE278688'];write_json(SOURCE/'analysis_spec.json',sp)
    vv=json.loads((SOURCE/'validation.json').read_text());vv['figure_status']='PASS';vv['figure_validation']=v;write_json(SOURCE/'validation.json',vv)
    from sc_source_three_cohorts_v1 import broad
    profiles=read(SOURCE/'sc_celltype_profiles.tsv');counts=read(SOURCE/'author_identity_counts.tsv');ann=read(SOURCE/'sc_annotation_map.tsv');ar=[]
    for _,rr in counts.iterrows():
        row={'cancer':'PDAC','cohort':'GSE242230','author_celltype':rr.celltype,'author_celltype_specific':rr.cell_type_specific,'coarse_celltype':broad(rr.resolved_celltype),'identity_status':rr.identity_status,'n_cells':rr.n_cells,'annotation_origin':'Author supplied broad and specific fields','mapping_evidence':'Specific normal or malignant identity takes precedence','version':VERSION,'status':'DONE','reason':rr.reason};ar.append(row)
    # Idempotent rebuild retains exact original label pairs, not invented author strings.
    write(SOURCE/'sc_annotation_map.tsv',pd.concat([ann[ann.cohort!='GSE242230'],pd.DataFrame(ar)],ignore_index=True))
    coverage=[]
    for ct in ['Malignant (author)','Normal epithelial (author)']:
        dd=profiles[(profiles.cohort=='GSE242230')&(profiles.celltype==ct)];coverage.append({'cohort':'GSE242230','identity':ct,'cells_total':int(dd.n_cells_total.max()),'source_labels_total':int(dd.n_source_labels_total.max()),'eligible_source_labels':int(dd.n.max()),'cells_in_eligible_labels':int(dd.n_cells_eligible.max()),'genes_measured':int(dd.effect.notna().sum()),'reason':'Unequal source-label coverage;descriptive means are not matched tumor-normal comparison'})
    write(SOURCE/'epithelial_identity_coverage.tsv',pd.DataFrame(coverage))
    p=SOURCE/'README_CN.md';p.write_text(p.read_text().replace('更新全候选图与比较表；不扩展机制。','全候选图与比较表已更新；另两队列恶性标签接回仍为NEEDS_REVIEW。\n\n[全部来源图](figures/README_CN.md) · [身份取证](author_identity_evidence.tsv) · [作者标签计数](author_identity_counts.tsv)'),encoding='utf8',newline='\n')
    stage(SOURCE,'06_EXTERNAL','GSE242230 author identity DONE;GSE263733/GSE278688 unresolved;all687 figures complete','350normal separated;10783author malignant;stats verified','Obtain exact-cell author malignancy labels for other2cohorts',status='PARTIAL')
    stage(INTEGRATION,'07_INTEGRATION','715relations687genes updated source identity;patient numbers unchanged','Strict identity comparisons;all rows retained','Use v4 interpretation;do not call unresolved ductal malignant')
    report=f'''# PDAC 当前进度：作者细胞身份整改 v4

## 本轮问题

接回作者恶性/正常标签，并更正此前把三队列统一为上皮/导管后丢失的身份解释。

## 输入与范围

原687基因、148777细胞、715关系保留。GSE242230采用cell_type_specific；其他两队列已核对现有GEO元数据只有导管大类。原始计数和完整逐细胞映射留服务器。

## 实际结果

- GSE242230：宽类Malignant的11133细胞中，细分为5660 Classical、5123 Basal和350 Normal Epithelial。现10783个作者恶性与350个作者正常上皮分开计算；350/11133=3.14%。
- GSE263733与GSE278688：保留Ductal (unresolved)，不能将肿瘤组织来源直接等同恶性身份。作者论文做过恶性分析，但当前没有接回完整逐细胞判定。
- 三队列严格同首位类别{s['all3_same_top_genes']}基因，同首位且各自保持率≥80%为{s['all3_top_bootstrap_ge080_genes']}基因（旧158）；当前250直接基因中为{s['all3_top_bootstrap_ge080_current_genes']}（旧65）。身份类别不等价导致的变化不称生物学证据消失。
- GSE242230拆分后{s['GSE242230_broad_lineage_top_changes']}基因的首位跨宽谱系改变；完整标签变化见source_identity_changes.tsv。正常上皮350细胞中291细胞分布于6个合格来源标签，恶性10783细胞中10759细胞分布于23个合格来源标签，两组不是匹配对照。
- 配对代谢物仍307项/51项P<0.05/0项q<0.05；当前RNA仍247项可评估/91项P<0.05/40项q<0.05。357直接、358条件关系及所有关联P/q不变。
- ABAT在GSE242230的首位为{ab.GSE242230_top_celltype}；完整候选来源与比较表已更新，不以ABAT等少数候选决定规则。

## 新手解释

“作者注释恶性”保留作者的证据来源；“身份未定”是尚未得到逐细胞判定，既不称正常也不称恶性。细分正常上皮标签优先于冲突的宽标签。旧29个上皮/导管稳定来源基因不能写成29个三队列恶性来源复现。

## 限制/反证

本轮未自行运行CNV或做机制分析。两队列身份连接仍待解决；临床身份去重与外部代谢关系验证未因本次整改完成。正常上皮350细胞分布及每类≥3来源标签/每标签≥20细胞阈值会限制可评估性；不可评估不当作零。

## 当前决定

GSE242230已完成作者身份接回、重新计算、独立核对和全候选图。06身份总任务仍PARTIAL；另两队列NEEDS_REVIEW。v4来源解释优先，旧数值不覆盖。PR #2仍草稿、未合并。

## 下一步

取得其他两队列可连接到细胞条码的作者恶性/正常判定后，另建版本接入；不自动用表达高低或肿瘤样本标签替代。

## 复现命令

服务器run_author_identity_v4_server.py；本地deliver_author_identity_v4.py。保留完整参数、输入哈希、逐细胞连接核查和旧新比较。

- [身份整改结果](../../results/PDAC/06_EXTERNAL/{SOURCE.name}/README_CN.md)
- [来源图](../../results/PDAC/06_EXTERNAL/{SOURCE.name}/figures/README_CN.md)
- [关系比较](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/candidate_relations_integrated.tsv)
- [基因比较](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/candidate_genes_integrated.tsv)
- [恶性与正常表达汇总](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/author_malignant_normal_gene_profiles.tsv)
- [完整旧新来源比较](../../results/PDAC/07_INTEGRATION/{INTEGRATION.name}/source_identity_changes.tsv)
'''
    (REPO/'reports/PDAC/CURRENT_STATUS_CN.md').write_text(report,encoding='utf8',newline='\n');checksum(SOURCE);checksum(INTEGRATION)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['record','integrate','workbook','finalize'],required=True);a=ap.parse_args();records.VERSION=VERSION
    {'record':record_source,'integrate':integrate_results,'workbook':workbook_data,'finalize':finalize}[a.mode]()
    ix=REPO/'coordination/stages/PDAC.tsv';dd=read(ix)
    dd.loc[dd.run_id==SOURCE.name,'code_path']='code/pdac/run_author_identity_v4_server.py'
    dd.loc[dd.run_id==INTEGRATION.name,'code_path']='code/pdac/deliver_author_identity_v4.py'
    dd.to_csv(ix,sep='\t',index=False,lineterminator='\n')
