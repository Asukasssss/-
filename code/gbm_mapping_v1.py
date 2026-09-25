from pathlib import Path
import subprocess,io,json,hashlib,re
import pandas as pd
R=Path(__file__).resolve().parents[1];RUN='20260925T103000Z_discovery_v1';O=R/'results/GBM/02_MAPPING'/RUN;O.mkdir(parents=True,exist_ok=False);S=O/'sources';S.mkdir()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
specs=[('BRCA','origin/analysis/brca-functional-review-20260919','results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/relation_evidence.tsv'),('PRAD','origin/analysis/prad-discovery-20260922','results/PRAD/02_MAPPING/20260922T142000Z_discovery_v1/direct_relations.tsv'),('COAD','origin/analysis/coad-initial','results/COAD/02_MAPPING/20260922T053000Z_paired101_mapping_v1/main_relations.tsv'),('PDAC','origin/analysis/pdac-initial','results/PDAC/02_MAPPING/20260922T055000Z_paired51_mapping_v2/direct_relations_v2.tsv')]
frames=[];manifest=[];ids={}
for cancer,ref,path in specs:
 commit=subprocess.check_output(['git','rev-parse',ref],cwd=R,text=True).strip();raw=subprocess.check_output(['git','show',commit+':'+path],cwd=R);p=S/(cancer+'_reviewed_mapping.tsv');p.write_bytes(raw);d=pd.read_csv(io.BytesIO(raw),sep='\t').fillna('');manifest.append(dict(source_id=cancer+'_mapping',path_or_url=path,source_commit=commit,sha256=sha(p)))
 for _,z in d.iterrows():
  if cancer=='BRCA' and z.mapping_status!='DIRECT_ACCEPTED':continue
  if cancer=='PRAD' and z.mapping_status!='DIRECT_ACCEPTED':continue
  if cancer=='COAD' and z.pool!='DIRECT':continue
  gene=z.gene
  if cancer=='COAD':ids[gene]=z.human_gene_id
  elif cancer=='PRAD':ids.setdefault(gene,z.stable_gene_id)
  accession=z.get('uniprot_accession','');urls=z.get('source_urls',z.get('source_url',''));matches=re.findall(r'uniprotkb/([A-Z0-9]+)/',str(urls))
  if not accession and len(set(matches))==1:accession=matches[0]
  if accession:ids.setdefault(gene,'UniProt:'+accession)
  frames.append(dict(metabolite_name=z.get('metabolite_name',z.get('feature_name','')),metabolite_key=z.metabolite_key,gene=gene,uniprot_accession=accession,source_urls=urls,relationship_type=z.get('relationship_type',z.get('relation_types',z.get('relation_type',''))),reaction=z.get('reaction',''),role=z.get('role_in_written_reaction',z.get('reaction_roles','')),evidence_origin=cancer,source_commit=commit,identity_note=z.get('identity_note',z.get('reason','')),source_row_id=z.get('evidence_id',z.get('relation_id',''))))
wpath=R/'results/GBM/01_CAMP'/RUN/'metabolite_workpool.tsv';w=pd.read_csv(wpath,sep='\t');ref=pd.DataFrame(frames);ref.metabolite_name=ref.metabolite_name.str.strip();assert w.feature_id.is_unique
# Exact name AND stable chemical key. No class-level lipid expansion or fuzzy aliases.
ev=w[['feature_id','metabolite_name','metabolite_key']].merge(ref,on=['metabolite_name','metabolite_key'],how='inner',validate='one_to_many');ev['stable_gene_id']=ev.gene.map(ids);ev['cancer']='GBM';ev['mapping_version']='gbm_mapping_v1';save(ev,O/'relation_evidence.tsv')
rows=[]
for (fid,gene),d in ev.groupby(['feature_id','gene'],sort=True):
 sid=ids.get(gene);r=d.iloc[0];rows.append(dict(cancer='GBM',cohort='CAMP_GBM',stage_id='02_MAPPING',run_id=RUN,analysis_version='gbm_mapping_v1',feature_id=fid,relation_id=fid+'|'+str(sid),metabolite_name=r.metabolite_name,metabolite_key=r.metabolite_key,gene=gene,stable_gene_id=sid,mapping_status='DIRECT_ACCEPTED',status='DONE' if sid else 'NEEDS_REVIEW',source_urls=' ; '.join(sorted(set(d.source_urls))),evidence_origins=';'.join(sorted(set(d.evidence_origin))),n_evidence_rows=len(d),discovery_gate='CONDITIONAL_NORMAL_LABEL_CONFLICT',reason='reviewed exact-name plus chemical-key human biochemistry reuse;not functional proof'))
e=pd.DataFrame(rows);assert not e.duplicated(['feature_id','stable_gene_id']).any();save(e,O/'direct_relations.tsv');genes=e[['gene','stable_gene_id']].drop_duplicates();assert genes.gene.is_unique;save(genes,O/'genes_unique.tsv')
c=w[['feature_id','metabolite_key','metabolite_name']].copy();c['n_direct_relations']=c.feature_id.map(e.groupby('feature_id').size()).fillna(0).astype(int);c['status']=c.n_direct_relations.map(lambda n:'DONE' if n else 'NEEDS_REVIEW');c['reason']=c.n_direct_relations.map(lambda n:'exact reviewed evidence reused' if n else 'no exact identity match in current reviewed catalog;not proof of no human enzyme');save(c,O/'mapping_coverage.tsv')
summary=dict(status='PARTIAL',workpool_features=len(w),mapped_features=int(c.n_direct_relations.gt(0).sum()),direct_relations=len(e),unique_genes=len(genes),unmapped_features=int(c.n_direct_relations.eq(0).sum()),unique_relation_keys=True,unique_gene_ids=True,statistics_reused_from_other_cancers=False,normal_discovery_conflict_unresolved=True)
(O/'validation.json').write_text(json.dumps(summary,indent=2));(O/'analysis_spec.json').write_text(json.dumps(dict(version='gbm_mapping_v1',join=['exact_outer_trimmed_name','exact_chemical_key'],pool='all480 conditional discovery P005 features',review_scope='reuse published BRCA PRAD COAD PDAC direct biochemical annotations only;not exhaustive new database mapping',stable_identity_priority='COAD explicit hsa gene ID then PRAD HGNC then source UniProt;exact gene symbols;no expression-based resolution'),indent=2));manifest.append(dict(source_id='GBM_workpool',path_or_url=str(wpath.relative_to(R)),sha256=sha(wpath)));save(pd.DataFrame(manifest),O/'source_manifest.tsv')
(O/'README_CN.md').write_text('# GBM直接生化映射\n\n## 本轮问题\n将全部条件性工作特征按直接生化证据映射到人类基因。\n\n## 输入与范围\n当前480个P<0.05特征，与已审定BRCA/PRAD/COAD/PDAC数据库关系按精确名称及化学键双重匹配。只复用生化证据，不复用其他癌种统计。\n\n## 实际结果\n'+json.dumps(summary,ensure_ascii=False,indent=2)+'\n\n## 新手解释\n多篇出处合并为一条特征—基因检验；不同峰不按化学键合并。\n\n## 限制/反证\n这不是穷尽式数据库新检索；未映射特征全保留NEEDS_REVIEW。上游正常标签冲突仍在。生化映射不证明GBM功能、通量或干预效果。\n\n## 当前决定\n全部直接关系进入肿瘤分析，全部基因进入单细胞，不按后续P二次筛选。\n\n## 下一步\n肿瘤内部关系与全基因细胞来源。\n\n## 复现\npython code/gbm_mapping_v1.py；使用source_manifest冻结提交。\n',encoding='utf-8')
save(pd.DataFrame([dict(file=str(p.relative_to(O)),sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file()]),O/'checksums.tsv');print(json.dumps(summary))
