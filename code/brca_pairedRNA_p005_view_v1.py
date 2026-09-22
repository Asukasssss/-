"""User-requested raw-P exploratory RNA retention view; no statistical refitting."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,subprocess
repo=Path(__file__).resolve().parents[1];run=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_pairedRNA_p005_view_v1';version='pairedRNA_p005_view_v1';out=repo/'results/BRCA/03_PATIENT'/run;out.mkdir(parents=True)
src=repo/'results/BRCA/03_PATIENT/20260921T150000Z_mapping262_patient_v1/paired_RNA150.tsv';master=repo/'results/BRCA/07_INTEGRATION/20260921T154000Z_mapping262_integration_v1/genes156_all_history_comparison.tsv'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
hashes={p:sha(p) for p in [src,master]};d=pd.read_csv(src,sep='\t',dtype=str,keep_default_na=False);p=pd.to_numeric(d.p_value,errors='coerce');q=pd.to_numeric(d.q_value,errors='coerce');ef=pd.to_numeric(d.effect,errors='coerce');n=pd.to_numeric(d.n,errors='coerce');pos=pd.to_numeric(d.positive_pairs,errors='coerce');neg=pd.to_numeric(d.negative_pairs,errors='coerce')
original_cols=d.columns.tolist();d['RNA_retained_P005']=np.where(p.notna(),np.where(p<.05,'YES','NO'),'NOT_EVALUABLE');d['RNA_evidence_tier']=np.where(q<.05,'P_AND_Q_LT005',np.where(p<.05,'P_ONLY_EXPLORATORY',np.where(p.notna(),'P_GE005','NOT_EVALUABLE')));d['RNA_direction']=np.where(ef>0,'UP',np.where(ef<0,'DOWN','NOT_EVALUABLE'));d['fraction_pairs_in_mean_direction']=np.where(ef>0,pos/n,np.where(ef<0,neg/n,np.nan));d['retention_view_version']=version;d['retention_view_run']=run
s=d[p<.05].sort_values(['fraction_pairs_in_mean_direction','gene'],ascending=[False,True]);extra=d[(p<.05)&~(q<.05)].sort_values('gene');save(d,out/'paired_RNA150_retention_flags.tsv');save(s,out/'paired_RNA85_P005_retained.tsv');save(extra,out/'additional6_P_only.tsv')
m=pd.read_csv(master,sep='\t',dtype=str,keep_default_na=False);mcols=m.columns.tolist();new=m.merge(d[['gene','RNA_retained_P005','RNA_evidence_tier','retention_view_version']],on='gene',how='left',validate='one_to_one');new.loc[new.RNA_retained_P005.isna(),'RNA_retained_P005']='HISTORICAL_ONLY';assert len(new)==156 and new[mcols].equals(m);save(new,out/'genes156_RNA_retention_appended.tsv')
assert len(s)==85 and len(extra)==6 and d[original_cols].equals(pd.read_csv(src,sep='\t',dtype=str,keep_default_na=False));assert all(sha(f)==h for f,h in hashes.items())
summary=dict(planned=150,evaluable=int(p.notna().sum()),P_lt005=int((p<.05).sum()),up=int(((p<.05)&(ef>0)).sum()),down=int(((p<.05)&(ef<0)).sum()),q_lt005=int((q<.05).sum()),P_only=len(extra),added_genes=extra.gene.tolist(),old_statistics_unchanged=True,historical156_columns_preserved=True,new_tests=0,new_BH=0)
(out/'validation.json').write_text(json.dumps(summary,indent=2));(out/'analysis_spec.json').write_text(json.dumps(dict(version=version,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo).decode().strip(),scope='paired RNA only;current150;history6 retained',entry='two-sided original paired RNA P<0.05;exploratory',q='original full148 family q unchanged;no subset BH',sort='fraction of all evaluable pairs in mean effect direction descending;gene tie-break',other_layers='no change to metabolite-RNA or external thresholds;RNA not mandatory gate for functional work'),indent=2))
save(pd.DataFrame([dict(path=str(f.relative_to(repo)).replace('\\','/'),sha256=h) for f,h in hashes.items()]+[dict(path='code/'+Path(__file__).name,sha256=sha(Path(__file__)))]),out/'source_manifest.tsv')
text=f'''# 配对RNA：按P<0.05探索性保留

## 本轮问题
按用户要求，配对RNA这一层先保留P<0.05的结果。

## 输入与范围
沿用150基因配对RNA表，其中148可评估；通常45对。这里只改变工作清单，不重新拟合、改旧统计或缩小范围重做BH。

## 实际结果
85个保留：50个肿瘤RNA平均升高、35个降低。其中79个也满足原新版q<0.05，另6个仅满足P<0.05：BCAT1、CKM、CMPK1、PNP、SLC13A3、SOAT1。

## 新手解释
85个是探索性保留清单；其中6个标P_ONLY_EXPLORATORY，不能称为通过FDR。方向与同向配对比例均保留；列表按均值方向对应的配对比例从高到低排列，不是功能排名。

## 限制与反证
P>=0.05及缺测基因仍在全150/156表中。RNA筛选不作为其他功能证据的强制入场门槛。原代谢物差异、RNA相关和外部层统计均未修改。

## 当前决定
后续展示配对RNA探索名单使用85个；q值作为证据强弱补充并列显示，旧79个结果保留。

## 下一步
使用附加标记后的156基因表继续比较，不删除其余候选。

## 复现
python code/brca_pairedRNA_p005_view_v1.py；仅读取公开汇总。输出保留上游统计版本字段，新视图另列版本。验证原统计列与完整156历史列不变。
''';(out/'README_CN.md').write_text(text,encoding='utf-8');(out/'.gitattributes').write_text('* -text\n')
save(pd.DataFrame([dict(file=f.name,sha256=sha(f)) for f in sorted(out.iterdir()) if f.is_file()]),out/'checksums.tsv')
ix=repo/'coordination/stages/BRCA.tsv';index=pd.read_csv(ix,sep='\t',keep_default_na=False);row=dict(cancer='BRCA',stage_id='03_PATIENT',run_id=run,analysis_version=version,status='DONE',scope='Paired RNA P<0.05 view:85 genes;50up35down;6 P-only;all150/156 preserved',result_path=str(out.relative_to(repo)).replace('\\','/'),code_path='code/'+Path(__file__).name,git_branch='analysis/brca-functional-review-20260919',reason='User requested exploratory P threshold;no recomputation or subset BH',next_action='Use RNA retention flags alongside original q;not a mandatory functional gate');save(pd.concat([index,pd.DataFrame([row])],ignore_index=True).sort_values(['stage_id','run_id']),ix)
print(run);print(json.dumps(summary))
