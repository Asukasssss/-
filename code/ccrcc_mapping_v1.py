"""Reuse exact reviewed chemical relationships, never another cancer's statistics."""
import argparse,json,hashlib,shutil
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main(run,cohort):
 out=R/'results/ccRCC/02_MAPPING'/run;out.mkdir(parents=True,exist_ok=False);src=out/'sources';src.mkdir()
 wp=R/'results/ccRCC/04_ROBUSTNESS'/run/'metabolite_workpool.tsv';work=pd.read_csv(wp,sep='\t').fillna('')
 inputs={
 'BRCA_reviewed_evidence.tsv':R.parent/'camp-cancer-collaboration/results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/relation_evidence.tsv',
 'PRAD_reviewed_evidence.tsv':R.parent/'camp-prad-analysis/results/PRAD/02_MAPPING/20260922T142000Z_discovery_v1/relation_evidence.tsv',
 'BRCA_gene_identity.tsv':R.parent/'camp-cancer-collaboration/results/BRCA/06_EXTERNAL/20260922T082327Z_sc156_umap_v1/gene_identity.tsv',
 'PRAD_gene_identity.tsv':R.parent/'camp-prad-analysis/results/PRAD/02_MAPPING/20260922T142000Z_discovery_v1/gene_identity_dictionary.tsv',
 'ccRCC_history_relations.tsv':R.parent/'pancancer_metabolomics_direct_matrix_20260716/results/CAMP_per_cancer_repro_v1_20260910/reused_by_cancer/ccRCC/locked_relations.tsv'}
 for name,p in inputs.items():shutil.copyfile(p,src/name)
 br=pd.read_csv(src/'BRCA_reviewed_evidence.tsv',sep='\t').fillna('');pr=pd.read_csv(src/'PRAD_reviewed_evidence.tsv',sep='\t').fillna('')
 br['reused_library']='BRCA_mapping190_v2';pr['reused_library']='PRAD_mapping_v1'
 hist=pd.read_csv(src/'ccRCC_history_relations.tsv',sep='\t').fillna('')
 h=hist.rename(columns={'canonical_metabolite':'metabolite_name','reaction_summary':'reaction','relation_type':'relationship_type','human_relation_evidence':'source_evidence','identity_requirement':'reason'})
 h['reused_library']='ccRCC_locked_20260910';h['mapping_status']='DIRECT_ACCEPTED'
 ev=pd.concat([br,pr,h],ignore_index=True).fillna('')
 dictionary={}
 for _,g in pd.read_csv(src/'PRAD_gene_identity.tsv',sep='\t').iterrows():dictionary[g.gene]=g.stable_gene_id
 for _,g in pd.read_csv(src/'BRCA_gene_identity.tsv',sep='\t').iterrows():dictionary[g.gene]=g.feature_id
 rows=[]
 for _,w in work.iterrows():
  hits=ev[ev.metabolite_name.eq(w.metabolite_name)&ev.metabolite_key.eq(w.metabolite_key)] if w.metabolite_key else ev.iloc[:0]
  for _,e in hits.iterrows():
   rows.append(dict(cancer='ccRCC',cohort=cohort,feature_id=w.feature_id,metabolite_key=w.metabolite_key,metabolite_name=w.metabolite_name,gene=e.gene,
    stable_gene_id=dictionary.get(e.gene,''),uniprot_accession=e.uniprot_accession,mapping_status=e.mapping_status,relationship_type=e.relationship_type,
    reaction=e.reaction,role=e.get('role','') or e.get('role_in_written_reaction',''),source_url=e.source_url,source_evidence=e.source_evidence,
    reason=e.reason,reused_library=e.reused_library,original_evidence_id=e.get('evidence_id','') or e.get('original_evidence_id','')))
 evidence=pd.DataFrame(rows).drop_duplicates(['feature_id','gene','reaction','mapping_status','source_url']);save(evidence,out/'relation_evidence.tsv')
 direct=[]
 for (feature,gene),z in evidence[evidence.mapping_status.eq('DIRECT_ACCEPTED')].groupby(['feature_id','gene']):
  e=z.iloc[0];stable=e.stable_gene_id
  direct.append(dict(cancer='ccRCC',cohort=cohort,stage_id='02_MAPPING',run_id=run,analysis_version='ccrcc_mapping_v1',feature_id=feature,relation_id=feature+'|'+gene,metabolite_key=e.metabolite_key,metabolite_name=e.metabolite_name,gene=gene,stable_gene_id=stable,uniprot_accession=';'.join(sorted(set(z.uniprot_accession))),source_urls=' ; '.join(sorted(set(z.source_url))),status='DONE' if stable else 'NEEDS_REVIEW',reason='exact_name_and_chemical_key_reused_biochemistry_not_exhaustive' if stable else 'stable_gene_ID_missing'))
 d=pd.DataFrame(direct);assert not d.duplicated(['feature_id','gene']).any();assert d.stable_gene_id.ne('').all()
 save(d,out/'direct_relations.tsv');save(evidence[~evidence.mapping_status.eq('DIRECT_ACCEPTED')],out/'conditional_relations.tsv')
 coverage=work[['feature_id','metabolite_key','metabolite_name','p_value','q_value']].copy()
 coverage['n_direct_relations']=coverage.feature_id.map(d.groupby('feature_id').size()).fillna(0).astype(int)
 coverage['status']=coverage.n_direct_relations.map(lambda n:'DONE' if n else 'NEEDS_REVIEW')
 coverage['reason']=coverage.apply(lambda z:'DIRECT_MAPPED' if z.n_direct_relations else 'UNKNOWN_FEATURE_IDENTITY' if z.metabolite_name.startswith('X-') else 'NO_EXACT_REVIEWED_RELATION_IN_BOUNDED_LIBRARY',axis=1)
 save(coverage,out/'metabolite_mapping_status.tsv')
 genes=d[['gene','stable_gene_id']].drop_duplicates().sort_values('gene');genes['current_pool']=True;genes['history_only']=False;save(genes,out/'genes_unique.tsv')
 history_only=sorted(set(hist.gene)-set(genes.gene))
 members=pd.concat([genes,pd.DataFrame([dict(gene=g,stable_gene_id=dictionary[g],current_pool=False,history_only=True) for g in history_only])],ignore_index=True)
 save(members,out/'gene_membership.tsv');save(hist,out/'historical_relations_preserved.tsv')
 val=dict(status='DONE',mapping_scope='bounded_exact_reuse_not_exhaustive',work_features=len(work),direct_features=d.feature_id.nunique(),relations=len(d),genes=d.gene.nunique(),history_union_genes=len(members),history_only_genes=history_only,unmapped_features=int(coverage.n_direct_relations.eq(0).sum()),unique_keys=True,statistics_transferred=False,stable_gene_ID_complete=True)
 (out/'validation.json').write_text(json.dumps(val,indent=2),encoding='utf-8')
 spec=dict(version='ccrcc_mapping_v1',entry='all primary patient-paired nominal P<0.05',scope='exact author harmonized name AND chemical key match to frozen reviewed human biochemistry; no fuzzy matches or pathway neighbours',statistics='no BRCA/PRAD numerical evidence transferred',chemical_identity='inherited author annotation, not reidentification',unmapped='retained NEEDS_REVIEW; absence of reviewed mapping is not absence of biochemical relationship',gene_identity='existing exact Ensembl dictionary preferred; otherwise reviewed human HGNC dictionary',source_versions=list(inputs))
 (out/'analysis_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
 (out/'README_CN.md').write_text(f'''# {cohort}直接生化映射\n\n问题：探索工作池有哪些可追溯的人类直接生化关系？\n\n输入范围：{len(work)}项配对P<0.05代谢特征；复用项目已审定BRCA/PRAD生化资料，只做名称及化学键双重精确匹配。\n\n实际结果：{val['direct_features']}个特征、{len(d)}条去重关系、{len(genes)}个基因；{val['unmapped_features']}项无本库充分直接关系或身份未知，保留待审。\n\n新手解释：复用的是人类酶/运输反应知识，不是BRCA或PRAD的统计支持。\n\n限制：本轮有边界的映射不穷尽人类反应；未映射不代表没有相关基因。数据库推断和化学身份限制在证据长表保留。\n\n当前决定：全{len(d)}关系进入患者关联，全{len(genes)}基因进入RNA和单细胞，不能先按相关显著性删基因。\n\n下一步：患者级关联和RNA背景、全候选表达来源。\n\n复现：python code/ccrcc_mapping_v1.py --run {run} --cohort {cohort}；sources保存所用公开库快照。\n''',encoding='utf-8')
 save(pd.DataFrame([dict(path_or_url=str(p.relative_to(R)),sha256=sha(p)) for p in [wp,Path(__file__),*src.iterdir()]]),out/'source_manifest.tsv')
 save(pd.DataFrame([dict(file=p.relative_to(out).as_posix(),sha256=sha(p)) for p in out.rglob('*') if p.is_file() and p.name!='checksums.tsv']),out/'checksums.tsv')
 print(json.dumps(val))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--run',required=True);a.add_argument('--cohort',default='ccRCC4');p=a.parse_args();main(p.run,p.cohort)
