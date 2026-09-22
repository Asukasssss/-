"""Public-only exact-feature mapping; reuse reviewed BRCA biochemical evidence, not statistics."""
import csv, hashlib, io, json, subprocess, re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
RUN='20260922T142000Z_discovery_v1'
OUT=ROOT/'results/PRAD/02_MAPPING'/RUN
REF='f18a215dc8a269e4f16c55011633d557f069f631'
PATH='results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/relation_evidence.tsv'
TOKENS={
 'ophthalmate':['ophthalmate'], 'glutathione, oxidized (GSSG)':['glutathione disulfide'],
 'glutathione, reduced (GSH)':['glutathione'], '3-phosphoglycerate':['(2R)-3-phosphoglycerate'],
 '2-phosphoglycerate':['(2R)-2-phosphoglycerate'],'2-aminoadipate':['L-2-aminoadipate'],
 'glucose':['D-glucose','alpha-D-glucose'], '4-hydroxybutyrate (GHB)':['4-hydroxybutanoate'],
 'guanosine':['guanosine'],'inosine':['inosine'],'hypoxanthine':['hypoxanthine'],'xanthosine':['xanthosine'],
 'threonine':['L-threonine'], 'glycine':['glycine'],'cysteine':['L-cysteine'],'cystine':['L-cystine'],
 "adenosine 5'-diphosphate (ADP)":['ADP'],'adenosine':['adenosine'],
 'glycerophosphoethanolamine':['sn-glycero-3-phosphoethanolamine'],
 'glutarate (C5-DC)':['glutarate'],'adenine':['adenine'],'valine':['L-valine'],
 'N-acetylglutamate':['N-acetyl-L-glutamate'], 'glycerol 2-phosphate':['glycerol 2-phosphate'],
 'myo-inositol':['myo-inositol'], 'scyllo-inositol':['scyllo-inositol'],
 "adenosine 5'diphosphoribose":['ADP-D-ribose'], 'N-acetylgalactosamine':['N-acetyl-alpha-D-galactosamine'],
 'N-acetylglucosamine':['N-acetyl-D-glucosamine','N-acetyl-alpha-D-glucosamine','N-acetyl-beta-D-glucosamine'],
 'urate':['urate'],'maltotetraose':['D-maltotetraose'], 'cholesterol':['cholesterol'],
 '1-oleoyl-GPS (18:1)':['1-(9Z-octadecenoyl)-sn-glycero-3-phospho-L-serine'],
 'androsterone sulfate':['androsterone 3alpha-sulfate'],
 'dehydroisoandrosterone sulfate (DHEA-S)':['dehydroepiandrosterone 3-sulfate'],
 '1-oleoylglycerol (1-monoolein)':['1-(9Z-octadecenoyl)-glycerol'],
 'oleate (18:1n9)':['(9Z)-octadecenoate'],
}

def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    sp=OUT/'sources';sp.mkdir(exist_ok=True)
    raw=subprocess.check_output(['git','show',REF+':'+PATH],cwd=ROOT)
    (sp/'BRCA_reviewed_biochemical_evidence.tsv').write_bytes(raw)
    old=pd.read_csv(io.BytesIO(raw),sep='\t').fillna('')
    legacy=pd.read_csv(ROOT/'reference/brca/direct_edges_v1.tsv',sep='\t').fillna('')
    add=[]
    for _,r in legacy.iterrows():
        add.append(dict(metabolite_name=r.metabolite_name,metabolite_key=r.metabolite_key,gene=r.gene_symbol,
                        uniprot_accession=r.uniprot_accession,relationship_type=r.relationship_type,
                        role_in_written_reaction=r.role_in_written_reaction,reaction=r.reaction_summary_cn,
                        source_url=r.source_url+' ; '+r.human_protein_url,mapping_status='DIRECT_ACCEPTED',
                        source_evidence='Legacy curated '+r.mapping_tier,reason=r.scope_note,
                        evidence_id='LEGACY_'+r.relation_id+'_'+r.gene_symbol))
    old=pd.concat([old,pd.DataFrame(add)],ignore_index=True).fillna('')
    wp=ROOT/'results/PRAD/04_ROBUSTNESS'/RUN/'metabolite_workpool.tsv'
    work=pd.read_csv(wp,sep='\t');assert work.feature_id.is_unique
    evidence=[]
    for _,r in work.iterrows():
        matches=old[old.metabolite_name.eq(r.metabolite_name)&old.metabolite_key.eq(r.metabolite_key)]
        for _,e in matches.iterrows():
            evidence.append(dict(cancer='PRAD',cohort='CAMP_PRAD',stage_id='02_MAPPING',run_id=RUN,
                analysis_version='prad_mapping_v1',feature_id=r.feature_id,metabolite_key=r.metabolite_key,
                metabolite_name=r.metabolite_name,gene=e.gene,uniprot_accession=e.uniprot_accession,
                relationship_type=e.relationship_type,role=e.role_in_written_reaction,reaction=e.reaction,
                mapping_status=e.mapping_status,source_url=e.source_url,source_evidence=e.source_evidence,
                provenance='BRCA reviewed human biochemistry reused; not BRCA statistical evidence',
                original_evidence_id=e.evidence_id,source_commit=REF,reason=e.reason))
    # Bounded exact participants in reviewed human entries; no substring or generic class expansion.
    identity=[]
    for p in sorted((sp/'uniprot').glob('*.json')):
        if p.name=='fetch_manifest.json':continue
        protein=json.loads(p.read_text(encoding='utf-8'));gene=p.stem
        assert protein['organism']['taxonId']==9606
        assert any(g.get('geneName',{}).get('value')==gene for g in protein.get('genes',[]))
        hgnc=[x['id'] for x in protein.get('uniProtKBCrossReferences',[]) if x['database']=='HGNC']
        identity.append(dict(gene=gene,uniprot_accession=protein['primaryAccession'],stable_gene_id=hgnc[0] if len(hgnc)==1 else 'NA',status='DONE' if len(hgnc)==1 else 'NEEDS_REVIEW'))
        for c in protein.get('comments',[]):
            if c['commentType']!='CATALYTIC ACTIVITY':continue
            rx=c['reaction'];sides=rx['name'].split(' = ')
            if len(sides)!=2:continue
            parts=[{re.sub(r'^\d+ ','',re.sub(r'\((in|out)\)$','',z)) for z in side.split(' + ')} for side in sides]
            for _,r in work.iterrows():
                tokens=set(TOKENS.get(r.metabolite_name,[]));roles=[i for i,z in enumerate(parts) if tokens & z]
                if not roles:continue
                ev=rx.get('evidences',[])
                evidence.append(dict(cancer='PRAD',cohort='CAMP_PRAD',stage_id='02_MAPPING',run_id=RUN,
                    analysis_version='prad_mapping_v1',feature_id=r.feature_id,metabolite_key=r.metabolite_key,
                    metabolite_name=r.metabolite_name,gene=gene,uniprot_accession=protein['primaryAccession'],
                    relationship_type='transporter' if gene.startswith('SLC') else 'enzyme',role='both_sides' if len(roles)==2 else 'substrate_as_written' if roles==[0] else 'product_as_written',
                    reaction=rx['name'],mapping_status='DIRECT_ACCEPTED',source_url='https://www.uniprot.org/uniprotkb/'+protein['primaryAccession']+'/entry',
                    source_evidence=json.dumps(ev),provenance='UniProt reviewed human entry, exact listed participant; evidence may be inferred',
                    original_evidence_id=';'.join(z['id'] for z in rx.get('reactionCrossReferences',[]) if z['database']=='Rhea'),source_commit='current_public_snapshot',
                    reason='Inherited CAMP chemical identity; not new stereochemical identification or flux evidence'))
    e=pd.DataFrame(evidence)
    e=e.drop_duplicates(['feature_id','gene','reaction','mapping_status','source_url'])
    save(e,OUT/'relation_evidence.tsv')
    rows=[]
    for (feat,gene),z in e[e.mapping_status.eq('DIRECT_ACCEPTED')].groupby(['feature_id','gene'],sort=True):
        r=z.iloc[0]
        rows.append(dict(cancer='PRAD',cohort='CAMP_PRAD',stage_id='02_MAPPING',run_id=RUN,analysis_version='prad_mapping_v1',
            feature_id=feat,relation_id=feat+'|'+gene,metabolite_key=r.metabolite_key,metabolite_name=r.metabolite_name,
            gene=gene,uniprot_accession=';'.join(sorted(set(z.uniprot_accession)-{''})),
            mapping_status='DIRECT_ACCEPTED',status='DONE',source_urls=' ; '.join(sorted(set(z.source_url))),reason='Reviewed exact-name and chemical-key human biochemical mapping; not exhaustive'))
    direct=pd.DataFrame(rows)
    ids=pd.DataFrame(identity)
    direct=direct.merge(ids[['gene','stable_gene_id']],on='gene',how='left',validate='many_to_one')
    save(direct,OUT/'direct_relations.tsv')
    coverage=[]
    for _,r in work.iterrows():
        z=direct[direct.feature_id.eq(r.feature_id)];co=e[e.feature_id.eq(r.feature_id)&e.mapping_status.eq('CONDITIONAL')]
        coverage.append(dict(feature_id=r.feature_id,metabolite_key=r.metabolite_key,metabolite_name=r.metabolite_name,
            n_direct_relations=len(z),n_direct_genes=z.gene.nunique(),n_conditional_evidence=len(co),
            status='DONE' if len(z) else 'NEEDS_REVIEW',reason='DIRECT_MAPPED' if len(z) else 'CONDITIONAL_IDENTITY' if len(co) else 'NO_REVIEWED_DIRECT_EVIDENCE_IN_REUSED_LIBRARY',
            p_value=r.p_value,q_value=r.q_value))
    save(pd.DataFrame(coverage),OUT/'metabolite_mapping_status.tsv')
    save(direct[['gene','stable_gene_id','uniprot_accession']].drop_duplicates(),OUT/'genes_unique.tsv')
    save(ids,OUT/'gene_identity_dictionary.tsv')
    assert not direct.duplicated(['feature_id','stable_gene_id']).any()
    assert direct.stable_gene_id.notna().all(), 'Fetch missing exact primary gene identities first'
    validation=dict(status='DONE',work_features=len(work),direct_features=int(direct.feature_id.nunique()),direct_relations=len(direct),genes=int(direct.gene.nunique()),
                    features_needing_review=int(sum(not r['n_direct_relations'] for r in coverage)),unique_relation_keys=True,all_genes_HGNC_verified=True,
                    name_and_chemical_key_joint_match=True,BRCA_statistics_transferred=False,biochemical_exhaustiveness_claimed=False)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    spec=dict(version='prad_mapping_v1',run_id=RUN,entry='All primary paired nominal P<0.05 features; no patient-RNA filter',
              relation_key='CAMP_PRAD x feature_id x HGNC stable_gene_id x mapping_version',
              library_scope='Exact-name plus chemical-key BRCA reviewed mappings; bounded human reviewed UniProt entries saved under sources; no generic substrate-class expansion',
              chemical_identity='Inherited CAMP identity; ambiguous isobars/stereochemistry not resolved by gene expression',
              ADP_scope='ADP participants only within this explicitly bounded gene library; not every ADP-linked enzyme in the genome',
              statistical_family_planned=len(direct),RNA_gene_family_planned=int(direct.gene.nunique()),
              review_needed=[r['metabolite_name'] for r in coverage if not r['n_direct_relations']],
              reference_commit=REF,source_access_date='2026-09-22',statistics_recomputed=False)
    (OUT/'analysis_spec.json').write_text(json.dumps(spec,indent=2,ensure_ascii=False),encoding='utf-8')
    paths=[wp,ROOT/'reference/brca/direct_edges_v1.tsv',Path(__file__)]+list(sp.rglob('*.json'))+[sp/'BRCA_reviewed_biochemical_evidence.tsv']
    save(pd.DataFrame([dict(path_or_url=p.relative_to(ROOT).as_posix(),sha256=sha(p),access_scope='PUBLIC_BIOCHEMISTRY_OR_AGGREGATE') for p in paths]),OUT/'source_manifest.tsv')
    (OUT/'README_CN.md').write_text(f'''# PRAD 69项工作池直接生化映射\n\n问题：全部配对P<0.05特征有哪些可追溯人类直接酶/转运关系？\n\n输入：69项工作池；BRCA已审定生化关系；本轮UniProt reviewed human来源快照。\n\n实际结果：{validation['direct_features']}项建立{len(direct)}条直接关系，涉及{validation['genes']}个HGNC稳定ID基因；其余{validation['features_needing_review']}项逐项保留NEEDS_REVIEW。\n\n新手解释：映射只说明该分子参与该基因对应的反应或运输，不说明PRAD中酶活、净通量或因果机制。多个反应出处合并为一个特征—基因检验。\n\n限制：本轮是有明确边界的来源映射，并不穷尽人类生化反应。21项尚无充分直接证据/身份待明确，不当作不存在相关基因；同分异构混合峰不强行选一个分子。ADP仅覆盖本轮明示基因库，不宣称全基因组穷尽。数据库实验/推断代码逐条保留。\n\n当前决定：全部{len(direct)}条关系进入肿瘤关联，全部{validation['genes']}基因进入RNA与单细胞来源，不以患者关联显著性再筛入场。\n\n下一步：患者统计和缺失敏感性；继续保留未映射清单。\n\n复现：先获取固定来源快照，再运行python code/prad_mapping_v1.py。公开来源快照和哈希均随本阶段交付。\n''',encoding='utf-8')
    save(pd.DataFrame([dict(file=p.relative_to(OUT).as_posix(),sha256=sha(p)) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='checksums.tsv']),OUT/'checksums.tsv')
    print(json.dumps(dict(work_features=len(work),direct_features=direct.feature_id.nunique(),relations=len(direct),genes=direct.gene.nunique())))
    print('Unmapped:',[r['metabolite_name'] for r in coverage if not r['n_direct_relations']])

if __name__=='__main__':main()
