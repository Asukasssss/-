"""Adapt immutable PDAC aggregates to the cross-cancer report contract. No tests refitted."""
from pathlib import Path
import io,json,hashlib,subprocess
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
BASE='bcccc21ff60d6c57111efe72804d50db8020245e'
HANDOFF='45c6ae9ab489a9bd4ead4e47d28185785a448156'
OUT=ROOT/'results/PDAC/07_INTEGRATION/20260922T135400Z_report_adapter_v5'
M='results/PDAC/01_CAMP/20260922T121000Z_sop_v3_discovery_complete/'
P='results/PDAC/03_PATIENT/20260922T112700Z_sop_v3_patient/'
R='results/PDAC/03_PATIENT/20260922T120300Z_sop_v3_internal/'
I='results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/'
S='results/PDAC/06_EXTERNAL/20260922T131800Z_author_identity_v4/'
MAP='results/PDAC/02_MAPPING/20260922T112600Z_sop_v3_mapping/'
manifest=[]
def read(p):
 b=subprocess.check_output(['git','show',BASE+':'+p],cwd=ROOT)
 manifest.append(dict(path=p,input_commit=BASE,sha256=hashlib.sha256(b).hexdigest()))
 return pd.read_csv(io.BytesIO(b),sep='\t',dtype=str,keep_default_na=False)
def write(name,df):df.to_csv(OUT/name,sep='\t',index=False,lineterminator='\n')
def main():
 OUT.mkdir(parents=True,exist_ok=False)
 m=pd.concat([read(M+'metabolite_paired.tsv'),read(M+'metabolite_available_sensitivity.tsv')],ignore_index=True)
 for a,b in [('n_up','pairs_higher'),('n_down','pairs_lower'),('n_equal','pairs_equal')]:m[b]=m[a]
 write('metabolite_results.tsv',m)
 a=pd.concat([read(P+'tumor_association.tsv'),read(P+'tumor_association_available.tsv')],ignore_index=True)
 write('association_results.tsv',a)
 ca=pd.concat([read(P+'conditional_tumor_association.tsv'),read(P+'conditional_tumor_association_available.tsv')],ignore_index=True)
 write('conditional_association_results.tsv',ca)
 r=read(R+'paired_RNA.tsv')
 for a0,b in [('n_up','positive_pairs'),('n_down','negative_pairs'),('n_equal','equal_pairs')]:r[b]=r[a0]
 write('rna_results.tsv',r)
 h=read(R+'paired_RNA_history_supplement.tsv');write('rna_history_supplement.tsv',h)
 rel=read(I+'candidate_relations_integrated.tsv')
 for source,prefix in [('metabolite','metabolite_paired'),('association','CAMP'),('association_available','CAMP_available')]:
  for x,y in [('p','p_value'),('q','q_value')]:rel[prefix+'_'+y]=rel[source+'_'+x]
  rel[prefix+'_status']=rel.get(source+'_status',pd.Series('DONE',index=rel.index))
 write('candidate_relations.tsv',rel)
 write('direct_candidate_relations.tsv',rel[rel.mapping_status.eq('DIRECT')])
 g=read(I+'candidate_genes_integrated.tsv');write('candidate_genes.tsv',g)
 reader=g[['gene','stable_gene_id','in_current_pool','history_only','RNA_effect','RNA_p','RNA_q','RNA_n','RNA_background']].copy()
 write('candidate_genes_reader.tsv',reader)
 mp=read(MAP+'feature_dispositions.tsv');mp['mapping_state']=mp.mapping_status;write('metabolite_mapping_status.tsv',mp)
 stats=pd.concat([m,a,ca,r,h],ignore_index=True)
 prefix=list(read('templates/statistical_result.tsv').columns)
 for col in prefix:
  if col not in stats:stats[col]='NA'
 write('statistics_long.tsv',stats[prefix+[x for x in stats if x not in prefix]])
 st=read(S+'sc_source_stability.tsv')
 st['coverage_status']=pd.to_numeric(st.n_evaluable_celltypes).ge(2).map({True:'DONE',False:'NOT_EVALUABLE'})
 st['detection_status']=pd.to_numeric(st.max_detection,errors='coerce').ge(.01).map({True:'DONE',False:'NOT_EVALUABLE'})
 st['display_top']=st.top_celltype.where(st.status.eq('DONE') & st.coverage_status.eq('DONE') & st.detection_status.eq('DONE'),'暂不可定位')
 assert st.loc[st.status.eq('DONE'),'display_top'].ne('暂不可定位').all(), 'Baseline source status does not implement coverage/detection contract'
 write('sc_source_status.tsv',st)
 for name in ['sc_cross_study.tsv','sc_gene_coverage.tsv','sc_dataset_registry.tsv','author_identity_counts.tsv','epithelial_identity_coverage.tsv']:
  write(name,read(S+name))
 audit=read('results/PDAC/03_PATIENT/20260922T053500Z_internal_paired51_v2/sample_identity_audit.tsv');write('sample_audit_summary.tsv',audit)
 write('design_summary.tsv',pd.DataFrame([dict(cancer='PDAC',design='author_confirmed_T_N_pair',metabolite_pairs=11,RNA_pairs=11,tumor_association_units=21,all_tumor_specimens=27,all_normal_specimens=12,excluded_tumor_without_pair_key=6,clinical_identity='NOT_RECERTIFIED',aliquot_identity='NOT_RECERTIFIED',cohort_overlap='NEEDS_REVIEW',source='GSM join to explicit author pairing workbook')]))
 progress=[]
 for kind,frame in [('代谢物',m),('直接关联',a),('条件关联',ca),('当前RNA',r),('补充RNA',h)]:
  for family,t in frame.groupby('test_family',sort=False):
   progress.append(dict(step=kind,test_family=family,planned=len(t),evaluable=int(t.status.eq('DONE').sum()),P_lt005=int(pd.to_numeric(t.p_value,errors='coerce').lt(.05).sum()),q_lt005=int(pd.to_numeric(t.q_value,errors='coerce').lt(.05).sum()),status='DONE',reason='Numerical batch complete;NA rows retained',path='statistics_long.tsv'))
 for cohort,t in st.groupby('cohort',sort=False):
  progress.append(dict(step='单细胞来源',test_family=cohort,planned=687,evaluable=int(t.display_top.ne('暂不可定位').sum()),P_lt005='NA',q_lt005='NA',status='PARTIAL',reason='All candidates processed;individual source gaps and clinical overlap retained',path='sc_source_status.tsv'))
 write('progress_summary.tsv',pd.DataFrame(progress))
 gap=pd.DataFrame([
  dict(module='作者恶性身份',status='PARTIAL',reason='GSE242230已拆分；GSE263733/GSE278688缺可逐细胞连接的作者恶性判定'),
  dict(module='临床分型来源',status='NOT_RUN',reason='当前已接入数据未建立可靠临床分型连接；Classical/Basal细胞状态不等于患者分型'),
  dict(module='UMAP及全基因UMAP图册',status='NOT_RUN',reason='当前锁定来源分析未接入作者坐标；本轮交付完整点图与热图，不另建嵌入'),
  dict(module='队列独立性',status='NEEDS_REVIEW',reason='三套研究来源标签已保留，患者跨研究重叠未完成认证'),
  dict(module='机制及外部同关系验证',status='NOT_RUN',reason='本轮止于表达来源；不将单细胞表达当代谢关系验证')])
 write('module_gaps.tsv',gap)
 write('source_manifest.tsv',pd.DataFrame(manifest).drop_duplicates('path'))
 spec=dict(version='PDAC_report_adapter_v5',input_commit=BASE,handoff_commit=HANDOFF,new_statistical_tests=0,new_P_q=0,units='11 author pairs;21 distinct author case keys;source labels not certified clinical donors',reuse='All numerical statistics reused exactly; aliases and display status only',method_deviations='Metabolite exact signed-rank/10000 pair bootstrap preserved from locked prior run, not replaced by BRCA approximation/4000',scope='357 direct main;358 conditional supplement;250 current DIRECT;687 total genes')
 (OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
 (OUT/'validation.json').write_text(json.dumps(dict(status='DONE',genes=len(g),relations=len(rel),source_records=len(st),raw_data_read=False,numerical_values_reused=True),indent=2)+'\n',encoding='utf8')
 (OUT/'README_CN.md').write_text('# PDAC 统一展示适配\n\n## 本轮问题\n接续已有计算，补统一展示。\n## 输入与范围\n固定输入 '+BASE+'，357直接与358条件关系、687基因。\n## 实际结果\n全部统计沿用，列别名和易读来源状态适配。源表哈希见清单。\n## 新手解释\n主报告以直接关系为主，条件及历史池完整保留；q按原家族分开。\n## 限制/反证\n临床身份、分型、作者坐标和两队列恶性身份缺项见module_gaps.tsv。旧精确配对检验及10000次区间不为排版替换。\n## 当前决定\n统一报告不新增统计，不延伸机制。\n## 下一步\n渲染主报告、完整Excel、全来源图册并验收。\n## 复现\npython code/pdac/prepare_report_v5.py（新目录）\n',encoding='utf8')
 print(json.dumps({'status':'DONE','output':str(OUT),'progress':progress},ensure_ascii=False))
if __name__=='__main__':main()
