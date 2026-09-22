"""Record actual resumed batches without overwriting the preflight or v2 results."""
from pathlib import Path
import argparse,json,subprocess,hashlib
import numpy as np
import pandas as pd
from prepare_sop_v3 import REPO,ROOT,TEMPLATES,VERSION,SOP_SHA,digest,write_json

INTERNAL=ROOT/'03_PATIENT/20260922T120300Z_sop_v3_internal'
SOURCE=ROOT/'06_EXTERNAL/20260922T120400Z_sop_v3_source'
DISCOVERY=ROOT/'01_CAMP/20260922T121000Z_sop_v3_discovery_complete'
INTEGRATION=ROOT/'07_INTEGRATION/20260922T121100Z_sop_v3_integrated'
PRE_ROOT=ROOT/'01_CAMP/20260922T112500Z_sop_v3_discovery'
MAP=ROOT/'02_MAPPING/20260922T112600Z_sop_v3_mapping'

def read(p):return pd.read_csv(p,sep='\t')
def write(path,df):
    template=TEMPLATES/path.name
    if template.exists():
        cols=list(read(template).columns)
        for c in cols:
            if c not in df:df[c]=np.nan
        df=df[cols+[c for c in df if c not in cols]]
    df.to_csv(path,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def checksum(folder):
    # Canonical LF bytes match the Git object and the manifest on every platform.
    for p in folder.rglob('*'):
        if p.is_file() and p.suffix in ['.tsv','.json','.md','.txt','.svg','.html']:p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
    write(folder/'checksums.tsv',pd.DataFrame([{'file':p.relative_to(folder).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='checksums.tsv']))
def stage(folder,stage_id,scope,reason,next_action,status='DONE'):
    path=REPO/'coordination/stages/PDAC.tsv';d=read(path);d=d[~((d.stage_id==stage_id)&(d.run_id==folder.name))]
    row={'cancer':'PDAC','stage_id':stage_id,'run_id':folder.name,'analysis_version':VERSION,'status':status,'scope':scope,'result_path':folder.relative_to(REPO).as_posix(),'code_path':'code/pdac/run_sop_v3_server.py','git_branch':'analysis/pdac-initial','reason':reason,'next_action':next_action}
    pd.concat([d,pd.DataFrame([row])],ignore_index=True).sort_values(['stage_id','run_id']).to_csv(path,sep='\t',index=False,lineterminator='\n')
def readme(folder,title,scope,result,limits,nextstep,command):
    text=f'''# {title}

## 本轮问题

执行GitHub固定规范{SOP_SHA}。只处理PDAC，单细胞止于表达来源。

## 输入与范围

{scope}

## 实际结果

{result}

## 新手解释

不同分析的P/q分别保留。当前DIRECT池250基因，历史并集687基因；RNA差异和关联显著性均不是进入单细胞的门槛。不可评估是缺项，不是阴性。

## 限制/反证

{limits}

## 当前决定

本批已完成范围登记DONE；不覆盖历史统计，不将表达来源升级为同一代谢关系的外部验证或因果机制。

## 下一步

{nextstep}

## 复现命令

{command}
'''
    (folder/'README_CN.md').write_text(text,encoding='utf8',newline='\n')

def record_internal():
    s=json.loads((INTERNAL/'summary.json').read_text());v=json.loads((INTERNAL/'independent_validation.json').read_text());assert v['status']=='PASS'
    spec=json.loads((INTERNAL/'analysis_spec.json').read_text());spec['status']='DONE_INTERNAL_BATCH';spec['run_id']=INTERNAL.name;spec['frozen_protocol_run_id']='20260922T112300Z_sop_v3';write_json(INTERNAL/'analysis_spec.json',spec)
    write_json(INTERNAL/'validation.json',{'status':'PASS','server_exact_four_source_hash_gate':True,'RNA_labels_complete_unique':True,'private_pair_tables_saved_server_only':True,'gene_scope_current250_history437':True,'independent_validation':v,'server_unit_tests':'5/5PASS','private_tables_downloaded':False})
    readme(INTERNAL,'PDAC SOP v3：配对RNA与原有值覆盖补齐','11个作者明确配对，250当前直接基因及437条件/历史补充基因；两组BH分别计算。连续Affymetrix作者表达尺度，不再log。',f"当前250基因：247可评估，91项P<0.05，40项q<0.05。补充437基因：429可评估，151项P<0.05，0项q<0.05。独立ttest_rel、t区间、4000次整对bootstrap重放核对676行通过。四张真实样本连接表仅存服务器private。",'作者病例配对键不等于临床身份重新认证；同一分装和临床协变量仍未闭环。均值差是作者连续表达尺度，不称倍数。方法/基因检验集合均不同于旧538基因Wilcoxon，不能直接比较显著基因数量。','接入全关系表与全基因表，等待三队列来源。','在独占服务器目录执行`python3 run_sop_v3_server.py --mode internal --code-commit 321cd0f8f6ab8908a0d6daa190da7be57a8ed3a5`；然后`python3 verify_sop_v3_server.py --mode internal`。')
    stage(INTERNAL,'03_PATIENT','250current:247evaluable91P00540q005;437supplement:429evaluable151P0050q005','Actual server paired t and bootstrap completed and independently verified','Integrate source when complete')
    DISCOVERY.mkdir(exist_ok=True);supp=read(INTERNAL/'metabolite_description_supplement.tsv');assert not supp.duplicated(['metabolite_name','analysis_type']).any()
    outputs=[]
    for name,mode in [('metabolite_paired.tsv','primary'),('metabolite_available_sensitivity.tsv','both_preimputation_available')]:
        d=read(PRE_ROOT/name);m=supp[supp.analysis_type==mode].set_index('metabolite_name').loc[d.metabolite_name].reset_index();assert (d.n.to_numpy()==m.n.to_numpy()).all()
        for c in ['n_up','n_down','n_equal']:np.testing.assert_array_equal(d[c],m[c])
        for c in ['effect','median_delta']:
            other='mean_delta' if c=='effect' else c;ok=d[c].notna();np.testing.assert_allclose(d.loc[ok,c],m.loc[ok,other],rtol=0,atol=1e-12)
        # Fill previously unavailable descriptive means for n<6; never change P/q.
        d['effect']=m.mean_delta;d['mean_delta']=m.mean_delta;d['median_delta']=m.median_delta;d['n_nonzero']=m.n_nonzero
        d['missing_rate_tumor']=m.missing_rate_tumor;d['missing_rate_normal']=m.missing_rate_normal;d['n_both_available']=m.n_both_available;d['marginal_mask_status']='DONE_SERVER_SOURCE_CHECKED';d['run_id']=DISCOVERY.name
        d['reason']=np.where(d.n<8,'N_PAIRS_LT_8;description retained;oldP in p_original','Exact same input/effect/P/CI reused;BH over evaluable;availability fractions verified on server')
        d['direction_discordant']=((np.sign(d.effect)!=np.sign(d.median_delta))|(np.sign(d.effect)!=np.sign(d.n_up-d.n_down))).where(d.effect.notna(),np.nan)
        write(DISCOVERY/name,d);outputs.append(d)
    write(DISCOVERY/'metabolite_workpool.tsv',outputs[0][outputs[0].p_value<.05].copy());write(DISCOVERY/'metabolite_direction_ranking.tsv',outputs[0].assign(majority_fraction=outputs[0][['up_fraction','down_fraction']].max(axis=1)).sort_values(['majority_fraction','n','p_value'],ascending=[False,False,True]))
    audit=read(PRE_ROOT/'sample_audit_summary.tsv');idx=audit.audit_item=='private_SOP_table_export';audit.loc[idx,'status']='DONE';audit.loc[idx,'reason']='Four private SOP mapping tables created in server internal run;never exported';write(DISCOVERY/'sample_audit_summary.tsv',audit)
    sp=json.loads((PRE_ROOT/'analysis_spec.json').read_text());sp.update(status='DONE_DISCOVERY_BATCH',run_id=DISCOVERY.name,server_supplement_run=INTERNAL.name);write_json(DISCOVERY/'analysis_spec.json',sp)
    manifest=read(PRE_ROOT/'source_manifest.tsv');extra=pd.DataFrame([{'source_id':'server_description_supplement','path_or_url':(INTERNAL/'metabolite_description_supplement.tsv').relative_to(REPO).as_posix(),'sha256':digest(INTERNAL/'metabolite_description_supplement.tsv'),'access_scope':'PUBLIC_AGGREGATE'}]);write(DISCOVERY/'source_manifest.tsv',pd.concat([manifest,extra],ignore_index=True))
    write_json(DISCOVERY/'summary.json',{'pairs':11,'metabolites':307,'primary_evaluable':307,'workpool_P005':51,'primary_q005':0,'available_evaluable_min8':193,'available_P005':27,'available_q005':0,'scale':'author processed;not log2FC'})
    write_json(DISCOVERY/'validation.json',{'status':'PASS','all614_direction_counts_and_n_match_server':True,'all_existing_means_and_medians_match_server':True,'primary_and_sensitivity_P_q_unchanged_from_v3_preflight':True,'per_tissue_missing_rates_now_from_server':True,'private_maps_server_only':True,'signedrank_original_validation':'Inherited verified exact sign-enumeration batch;not independently rerun now'})
    readme(DISCOVERY,'PDAC SOP v3：全量配对代谢物与缺失敏感性','全307项，11个明确作者配对；最低8对可计算。新q基于全部可评估P。','主307项均可评估，51项P<0.05，0项q<0.05；双方作者原有值敏感性193项可评估，27项P<0.05，0项q<0.05。服务器补齐每侧缺失率及低n描述。','作者data有值不是原始质谱检出认证。原10000次整对均值CI按规范复用，不为了新4000默认值重算。患者独立性局限见样本汇总。','与固定映射、RNA和来源表整合。','`python code/pdac/record_sop_v3_resume.py --mode internal`；P/CI原代码及新q版本见source_manifest与analysis_spec。')
    stage(DISCOVERY,'01_CAMP','307metabolites51P0050q005;availability193evaluable27P0050q005','Actual server marginal availability/descriptive supplement verified;originalP preserved','Integrate complete discovery table')
    stage(DISCOVERY,'04_ROBUSTNESS','Both-author-available paired sensitivity at min8;193evaluable','Same data sensitivity,not independent validation','Keep clinical covariate limitations')
    checksum(INTERNAL);checksum(DISCOVERY)
    status=REPO/'reports/PDAC/CURRENT_STATUS_CN.md'
    status.write_text(f'''# PDAC 当前进度：SOP配对RNA已完成，三队列来源运行中

服务器已恢复。本轮只到基因表达来源，不开展机制分析。旧版本全部保留。

|环节|最新状态|
|---|---|
|样本及配对代谢物|11对、307项，51项P<0.05，0项q<0.05；缺失敏感性193项可评估，已补齐每侧缺失率|
|直接映射|357直接关系/250当前基因；358条件关系；历史并集687基因|
|肿瘤内部关联|原21作者单位统计复用，新BH已完成；直接主1条、条件主4条q<0.05|
|新配对RNA|当前247/250可评估，91项P<0.05，40项q<0.05；补充429/437可评估，151项P<0.05，0项q<0.05|
|新单细胞来源|三队列687基因正在运行；未将旧pseudobulk结果当作新cellwise结果|
|整合与来源图|等待新来源数值；不能称整流程完成|

新RNA的676行t检验、t区间和bootstrap，以及两个BH族已独立复核；不是外部验证。原矩阵与四张完整样本表继续仅存服务器。

- [RNA P<0.05展示](../../results/PDAC/03_PATIENT/{INTERNAL.name}/RNA_P005_view.tsv)
- [配对RNA说明](../../results/PDAC/03_PATIENT/{INTERNAL.name}/README_CN.md)
- [配对代谢物说明](../../results/PDAC/01_CAMP/{DISCOVERY.name}/README_CN.md)

临床身份再认证、同一分装、协变量和跨研究患者去重仍未因此确认；PR保持草稿未合并。
''',encoding='utf8',newline='\n')
    print(json.dumps(s))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['internal'],required=True);a=p.parse_args();record_internal()
