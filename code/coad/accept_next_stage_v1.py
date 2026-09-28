"""Accept task refinement and reproduce overlap from pinned real aggregate sources."""
import argparse,csv,hashlib,io,json,subprocess,zipfile
from pathlib import Path
R=Path(__file__).resolve().parents[2]
def read(b):return list(csv.DictReader(io.StringIO(b.decode('utf-8-sig')),delimiter='\t'))
def sha(b):return hashlib.sha256(b).hexdigest()
def write(p,b):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('xb') as f:f.write(b if isinstance(b,bytes) else b.encode('utf-8'))
def dump(p,x):write(p,json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def table(p,rs):
 s=io.StringIO(newline='');w=csv.DictWriter(s,list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs);write(p,s.getvalue())
def main():
 p=argparse.ArgumentParser();p.add_argument('--zip',type=Path,required=True);p.add_argument('--run-id',required=True);a=p.parse_args();o=R/'results/COAD/07_INTEGRATION'/a.run_id;assert not o.exists()
 with zipfile.ZipFile(a.zip) as z:files={Path(i.filename).name:z.read(i) for i in z.infolist() if not i.is_dir()}
 for m in read(files['MANIFEST_SHA256.tsv']):assert len(files[m['file']])==int(m['bytes']) and sha(files[m['file']])==m['sha256']
 prov=json.loads(files['provenance.json']);bc=prov['brca_commit'];bp='results/BRCA/20260919_ER_availability_v3/er_availability_all_174.tsv'
 bb=subprocess.check_output(['git','show',bc+':'+bp],cwd=R);assert subprocess.check_output(['git','rev-parse',bc+':'+bp],cwd=R,text=True).strip()==prov['brca_relation_table_blob']
 cp=R/'results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1';cg=read((cp/'gene_function_review.tsv').read_bytes());cr=read((cp/'nominal39_annotated_exact.tsv').read_bytes());br=read(bb)
 coad={r['gene'] for r in cg};brca={r['gene'] for r in br};common=sorted(coad&brca);assert (len(coad),len(brca),len(br))==(35,117,174) and common==['BCAT2','NNMT','SLC6A6','UPP2']
 assert coad=={r['gene'] for r in read(files['COAD_35_gene_identities_extracted.tsv'])}
 assert {(r['gene'],r['relation_id']) for r in br}=={(r['gene'],r['relation_id']) for r in read(files['BRCA_174_relation_identities_extracted.tsv'])}
 bm={(r['gene'],r['metabolite']):r for r in br};pairs=[]
 for c in cr:
  b=bm.get((c['gene'],c['metabolite_name']))
  if b:
   assert b['relation_id'].split('|')[0]==c['metabolite_key']
   pairs.append({'gene':c['gene'],'metabolite_key':c['metabolite_key'],'metabolite_name':c['metabolite_name'],'COAD_n':c['n'],'COAD_rho':c['effect'],'COAD_p':c['p_value'],'COAD_q':c['q_value'],'BRCA_ER_n':b['primary_er_n'],'BRCA_ER_rho':b['primary_er_rho'],'BRCA_ER_p':b['primary_er_p'],'BRCA_ER_q':b['primary_er_q'],'BRCA_ER_availability_rho':b['partial_rank_rho'],'BRCA_ER_availability_q':b['q_value'],'same_direction':str(float(c['effect'])*float(b['primary_er_rho'])>0),'interpretation':'Descriptive same-name/same-source-key match; different models and populations; no independent mechanism validation'})
 assert len(pairs)==3
 for p0 in pairs:
  src=next(x for x in read(files['BRCA_COAD_same_named_metabolite3_comparison.tsv']) if x['基因']==p0['gene'])
  assert abs(float(src['COAD未调整rho_展示值'])-float(p0['COAD_rho']))<=.000051
  assert abs(float(src['COAD未调整q_展示值'])-float(p0['COAD_q']))<=.000051
  assert float(src['BRCA_ER_rho'])==float(p0['BRCA_ER_rho']) and float(src['BRCA_ER_q'])==float(p0['BRCA_ER_q'])
 for name,data in files.items():
  if not name.endswith('.xlsx'):write(o/'imported_delivery'/name,data)
 table(o/'shared_relations_exact.tsv',sorted(pairs,key=lambda r:r['gene']))
 table(o/'shared_genes.tsv',[{'gene':g,'COAD_metabolites':';'.join(sorted({r['metabolite_name'] for r in cr if r['gene']==g})),'BRCA_metabolites':';'.join(sorted({r['metabolite'] for r in br if r['gene']==g})),'scope':'COAD reviewed35 vs BRCA full117'} for g in common])
 write(o/'priority9_action_matrix.tsv',files['COAD_priority9_action_matrix.tsv'])
 both=sorted({r['gene'] for r in br if r['status']=='DONE_EXPLORATORY' and float(r['primary_er_q'])<.05 and float(r['q_value'])<.05})
 dump(o/'analysis_spec.json',{'version':'COAD_next_stage_accept_v1','scope':prov['scope'],'source_provenance':prov,'matching':'Exact gene, source metabolite name and source feature key; no remapping or tests','BRCA_same_relation_supported_in_both_ER_runs':both,'new_patient_statistics':False,'new_literature_review':False})
 table(o/'source_manifest.tsv',[{'source':str(a.zip),'sha256':sha(a.zip.read_bytes()),'version':'user delivery'}, {'source':bc+':'+bp,'sha256':sha(bb),'version':bc},{'source':str(cp.relative_to(R))+'/nominal39_annotated_exact.tsv','sha256':sha((cp/'nominal39_annotated_exact.tsv').read_bytes()),'version':prov['coad_commit']}])
 dump(o/'validation.json',{'status':'PASS','shared_genes':common,'shared_relations':3,'BRCA_both_ER_supported_genes':both,'checks':['source package manifest','pinned actual BRCA git blob','actual gene and relation identity sets','COAD exact stats instead of rounded display','same source feature keys','supplied numeric precision compatible'],'limitations':['No new independent functional review','No patient statistics','Not COAD458 scope','Excel retained locally unchanged, not rendered or edited']})
 report=f'''# COAD 下一阶段任务与跨癌对照接收

## 本轮问题

接续用户的九基因任务细化与跨癌对照，使用真实仓库高精度记录补齐来源，不重新开展已完成的全包功能审核。

## 输入与范围

COAD已核查35基因、39条名义P关系，对照BRCA完整117基因、174关系。BRCA固定提交{bc}；来源和原交付见analysis_spec.json。仅B账号写入COAD目录，不修改BRCA。

## 实际结果

共同基因为BCAT2、NNMT、SLC6A6、UPP2。前3个分别对应异亮氨酸、1-MNA、牛磺酸，名称及原特征键均相同，方向一致。shared_relations_exact.tsv补入COAD原高精度rho/P/q，未用四位小数替换原值。UPP2分别对应尿苷/尿嘧啶，仅共同基因。BRCA两轮ER分析中同一关系均q<0.05的基因为{', '.join(both)}。九基因任务表接收，9/21/5安排不变。

## 新手解释

共同基因不等于共同代谢关系。方向相同也不等于跨癌机制已证实。COAD未调整与BRCA调整ER的模型不同，不能比较是否跨过阈值来证明癌种差异。

## 限制/反证

交集之外不称癌种特异；这里不是COAD458全基因比较。任务细化复用功能证据，不计新增独立取证。压缩包内尚未运行、上传等状态是其交付时快照。

## 当前决定

任务细化与集合比较接收DONE，07_INTEGRATION整体PARTIAL。协变量执行采用另行锁定的config/coad_covariates_v1.json，草案不是已运行结果。

## 下一步

先完成33人全674关系M1/M2敏感性；M3原部位类别稀疏按锁定规则登记。HDC/GSTA4细胞来源等功能任务继续保留，未声称已做单细胞分析。

## 复现命令

`python code/coad/accept_next_stage_v1.py --zip <用户ZIP> --run-id <新唯一RUN_ID>`。需本地git含固定BRCA提交。仅汇总表与代码公开，原Excel留本机。
'''
 write(o/'README_CN.md',report);write(o/'DONE','Accepted scope complete; no new patient statistics.\n');print(json.dumps({'path':str(o),'shared_genes':common,'same_pairs':3,'BRCA_both_ER':both}))
if __name__=='__main__':main()
