"""Create aggregate reports and stage records; never reads patient matrices."""
from pathlib import Path
import sys,json,hashlib,pandas as pd
repo=Path(__file__).resolve().parents[1];mode=sys.argv[1]
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def register(stage,run,version,path,code,scope,status='DONE',reason='Historical results immutable; new family q separately recorded',next_action='Combine patient and external evidence with functional sources; no hard significance intersection'):
 p=repo/'coordination/stages/BRCA.tsv';d=pd.read_csv(p,sep='\t',keep_default_na=False);assert not ((d.stage_id==stage)&(d.run_id==run)).any()
 row=dict(cancer='BRCA',stage_id=stage,run_id=run,analysis_version=version,status=status,scope=scope,result_path=path,code_path=code,git_branch='analysis/brca-functional-review-20260919',reason=reason,next_action=next_action)
 save(pd.concat([d,pd.DataFrame([row])],ignore_index=True).sort_values(['stage_id','run_id']),p)
def checks(p):
 (p/'.gitattributes').write_text('* -text\n')
 files=[x for x in p.rglob('*') if x.is_file() and x.name!='checksums.tsv']
 save(pd.DataFrame([dict(file=str(x.relative_to(p)).replace('\\','/'),sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in sorted(files)]),p/'checksums.tsv')
if mode=='camp':
 run='20260921T150000Z_mapping262_patient_v1';path='results/BRCA/03_PATIENT/'+run;p=repo/path;d=pd.read_csv(p/'CAMP_relations262.tsv',sep='\t');e=pd.read_csv(p/'paired_RNA150.tsv',sep='\t');v=json.loads((p/'validation.json').read_text())
 assert len(d)==524 and len(e)==150 and d.groupby('test_family').relation_id.nunique().eq(262).all()
 for _,z in d.groupby('test_family'):assert z.p_value.notna().sum()==255
 new=e[~e.statistics_reused & e.p_value.notna()];sig=d[(d.test_family=='processed_Spearman262')&d.q_value.lt(.05)&~d.statistics_reused]
 av=d[d.test_family.eq('available_Spearman262')].set_index('relation_id');lines=[]
 for t in sig.itertuples():a=av.loc[t.relation_id];lines.append(f'|{t.gene}|{t.metabolite_name.strip()}|{t.effect:.3f}|{t.q_value:.4g}|{a.n:g}|{a.effect:.3f}|{a.q_value:.4g}|')
 text='''# 新版候选池：CAMP患者数据第一批

## 本轮问题
把190项代谢物入口形成的262条直接关系、150个基因接到CAMP患者数据。只补新增统计，不把旧117名单重新当作全部候选。

## 输入与范围
映射来源提交9e4387f。排除已有组织类型冲突后，肿瘤内部关联使用60个作者病例；RNA差异使用45对明确作者病例连接。两者是不同分析，不将60称为60对。

## 实际结果
|分析|计划|可评估|复用统计|新计算|新版q<0.05|
|---|---:|---:|---:|---:|---:|
'''
 for f,x in v['families'].items():text+=f"|{f}|{x['planned']}|{x['evaluable']}|{x['reused']}|{x['newly_calculated']}|{x['q_lt005']}|\n"
 text+='\n新增关系中主分析达到新版q<0.05的条目：\n\n|基因|代谢物|主分析ρ|主分析q|可用值n|可用值ρ|可用值q|\n|---|---|---:|---:|---:|---:|---:|\n'+'\n'.join(lines)
 text+=f'\n\n新增39基因中38个完成配对RNA分析，其中{int(new.q_value.lt(.05).sum())}个达到新版RNA q<0.05。全部150行包括不可评估项保留。\n'
 text+='''
## 新手解释
相关系数回答：不同肿瘤之间，某基因RNA较高时，代谢物是否倾向较高或较低。配对RNA差异回答：同一个病例的肿瘤相对正常，RNA是否改变。二者都不是酶活、通量或干预效果。RNA差异效应沿用作者表达尺度，不作为浓度倍数。

## 限制与反证
LPCAT4存在历史命名歧义（Q643R3与Q6ZWT7），暂挂对应4条关系及1个RNA结果，等待探针稳定ID核对；CHKB未覆盖的3条关系及RNA仍不可评估。没有猜别名。此次共有7条关系不可评估或待核对。
BH分别在255条可评估主分析、255条可评估敏感性分析和148个可评估RNA基因中计算。旧q保留previous_q，新q不能覆盖历史q。当前入口从同一CAMP选出，因此属于同队列探索，不是独立验证。敏感性分析使用作者data表可用值，不冒充已核验的原始仪器缺失标志。

## 当前决定
保留全部候选；ENPP2/LYPLA1的新增患者线索进入外部与功能对照，不据此宣布新靶点。不以患者显著性作为功能资料纳入门槛。

## 下一步
接入同批FUSCC/Tang全262关系结果、新39基因文献，以及保留历史6基因的156行比较表。

## 复现
源矩阵仅在server165。创建独占运行目录和`.running`（内容mapping262_patient_v1），将映射direct_relations.tsv放入source，复制本脚本及旧CAMP helper后执行：
`python3 brca_mapping262_patient_v1.py <new_run_directory>`
代码：code/brca_mapping262_patient_v1.py。运行路径及输入SHA256见source_manifest.tsv；脚本核对旧统计n/效应、配对t公式、BH与输入不变性。
'''
 (p/'README_CN.md').write_text(text,encoding='utf-8');checks(p);register('03_PATIENT',run,'mapping262_patient_v1',path,'code/brca_mapping262_patient_v1.py','262 relations x2;150 paired RNA;60 tumor cases/45pairs;old statistics reused;new BH')
