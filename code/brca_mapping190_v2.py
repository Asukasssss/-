"""Versioned, exact-key BRCA mapping; only reuses aggregate statistics.
Run: python code/brca_mapping190_v2.py RUN_ID
Public UniProt snapshots must first be fetched with brca_mapping190_fetch.py.
"""
import csv, json, hashlib, re, sys
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/BRCA/04_ROBUSTNESS/20260921T122021Z_paired_p_rank_v1/p_lt_005_sorted.tsv'
OLD=ROOT/'reference/brca/direct_edges_v1.tsv'
SUP=ROOT/'reference/mapping/supplement_v2_candidate_edges.tsv'
CACHE=ROOT/'reference/brca_mapping190_sources_v1'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None,sep='\t'):
    with p.open('w',encoding='utf-8-sig' if sep==',' else 'utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter=sep);w.writeheader();w.writerows(rows)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical_sha(p): return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def rel(p): return p.relative_to(ROOT).as_posix()

# Exact reaction participants, not substring or pathway membership. Names/keys
# remain inherited CAMP identities, never a new chemical identification.
EXACT={
'succinate':'succinate', 'creatinine':'creatinine',
'propionylcarnitine':'O-propanoyl-(R)-carnitine',
'isobutyrylcarnitine':'O-isobutanoyl-(R)-carnitine',
'stearoylcarnitine':'O-octadecanoyl-(R)-carnitine',
'oleoylcarnitine':'O-(9Z)-octadecenoyl-(R)-carnitine',
'laurate (12:0)':'dodecanoate','myristate (14:0)':'tetradecanoate','stearate (18:0)':'octadecanoate',
'oleate (18:1n9)':'(9Z)-octadecenoate','adrenate (22:4n6)':'(7Z,10Z,13Z,16Z)-docosatetraenoate',
'arachidonate (20:4n6)':'(5Z,8Z,11Z,14Z)-eicosatetraenoate',
'eicosapentaenoate (EPA; 20:5n3)':'(5Z,8Z,11Z,14Z,17Z)-eicosapentaenoate',
'docosahexaenoate (DHA; 22:6n3)':'(4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoate',
'1-palmitoyl-GPC (16:0)':'1-hexadecanoyl-sn-glycero-3-phosphocholine',
'1-stearoyl-GPC (18:0)':'1-octadecanoyl-sn-glycero-3-phosphocholine',
'1-oleoylglycerophosphocholine':'1-(9Z-octadecenoyl)-sn-glycero-3-phosphocholine',
'2-palmitoyl-GPC (16:0)':'2-hexadecanoyl-sn-glycero-3-phosphocholine',
'2-oleoylglycerophosphoethanolamine':'2-(9Z-octadecenoyl)-sn-glycero-3-phosphoethanolamine',
'1-palmitoylglycerol (1-monopalmitin)':'1-hexadecanoylglycerol',
'1-stearoylglycerol (18:0)':'1-octadecanoylglycerol',
'xanthosine':'xanthosine',
}
# Participant can be identified only conditionally. Do NOT enter accepted pool.
CONDITIONAL={
'methionine sulfoxide':('L-methionine (S)-S-oxide','MSRA','未区分甲硫氨酸亚砜S/R构型；MSRA只覆盖匹配的游离S-亚砜反应。'),
'N-acetylmethionine':('N-acetyl-L-methionine','ACY1','沿用旧条件性标记；原特征未完成L构型身份核实。'),
'2-hydroxyglutarate':('(R)-2-hydroxyglutarate','D2HGDH','D/L未区分，不把混合或未定构型直接指定给D2HGDH。'),
'erythronate':('D-erythronate','PGP','核实D构型及原峰身份；PGP催化磷酸赤藓糖酸脱磷酸，不是ALDH1A1表达关联的直接替代。'),
}

def main(run):
    out=ROOT/'results/BRCA/02_MAPPING'/run;out.mkdir(parents=True,exist_ok=False)
    rows=read(P); old=read(OLD); sup=read(SUP)
    assert len(rows)==190 and len(old)==174
    bykey={r['metabolite_key']:r for r in rows};byname={r['metabolite_name']:r for r in rows}
    assert len(bykey)==190 and len(byname)==190
    oldgenes={r['gene_symbol'] for r in old}; assert len(oldgenes)==117
    events=[]; proteins={}; alias_exclusions=[]
    for p in sorted(CACHE.glob('*.json')):
        data=json.loads(p.read_text(encoding='utf-8'))
        exact=[x for x in data['results'] if any(g.get('geneName',{}).get('value')==p.stem for g in x.get('genes',[])) and x.get('organism',{}).get('taxonId')==9606]
        assert len(exact)==1,(p,len(exact))
        proteins[p.stem]=exact[0]
        alias_exclusions.extend(dict(query_gene=p.stem,rejected_accession=x['primaryAccession'],reason='Primary human gene symbol differs; synonym hit excluded') for x in data['results'] if x not in exact)
    def add(m,g,kind,role,reaction,url,origin='NEW_REVIEW',accepted=True,note='',accession='',source_path='',evidence='',rhea='',historical=''):
        if m not in byname: return
        r=byname[m]
        events.append(dict(cancer='BRCA',stage_id='02_MAPPING',run_id=run,analysis_version='mapping190_v2',
            metabolite_key=r['metabolite_key'],metabolite_name=m,gene=g,relationship_type=kind,
            role_in_written_reaction=role,reaction=reaction,source_url=url,uniprot_accession=accession,
            evidence_origin=origin,mapping_status='DIRECT_ACCEPTED' if accepted else 'CONDITIONAL',
            status='DONE' if accepted else 'NEEDS_REVIEW',reason=note,
            identity_basis='Inherited CAMP feature; no reidentification; exact key join',
            source_path=source_path,source_evidence=evidence,rhea_id=rhea,historical_reaction_id=historical,
            relation_key=r['metabolite_key']+'|'+g,patient_analysis_this_round='NOT_RUN',functional_analysis_this_round='NOT_RUN'))
    for e in old:
        if e['metabolite_key'] in bykey:
            m=bykey[e['metabolite_key']]['metabolite_name'];assert m==e['metabolite_name']
            add(m,e['gene_symbol'],e['relationship_type'],e['role_in_written_reaction'],e['reaction_summary_cn'],e['source_url']+' ; '+e['human_protein_url'],
                'LEGACY_V1_REUSED',note=e['scope_note']+' 旧关系按既有审定复用，本轮未逐篇重审。',accession=e['uniprot_accession'],source_path=rel(OLD),evidence='Legacy curated tier '+e['mapping_tier'],historical=e['relation_id'])
    for e in sup:
        if e['metabolite_key'] in bykey:
            assert bykey[e['metabolite_key']]['metabolite_name']==e['metabolite_name']
            accepted=e['supplement_status']=='DIRECT_BIOCHEMICAL_RELATION_READY_NEXT_VERSION'
            add(e['metabolite_name'],e['gene'],'transporter' if e['role']=='direct_transport' or e['gene'].startswith('SLC') else 'enzyme',e['role'],e['reaction'],e['source_url'],
                'PREVIOUS_SUPPLEMENT_ACTIVATED' if accepted else 'PREVIOUS_CONDITIONAL_RETAINED',accepted,e['note'],e['uniprot_accession'],rel(SUP),'PubMed:'+e['pubmed_ids'],e['rhea_id'])
    def match_rxn(m,token,genes=None,accepted=True,note=''):
        hits=0
        for g,x in proteins.items():
            if genes and g not in genes: continue
            for c in x.get('comments',[]):
                if c['commentType']!='CATALYTIC ACTIVITY':continue
                rx=c['reaction']; reaction=rx['name']
                sides=reaction.split(' = ')
                if len(sides)!=2: continue
                def parts(s):return [re.sub(r'\((in|out)\)$','',v) for v in s.split(' + ')]
                roles=[i for i,s in enumerate(sides) if token in parts(s)]
                if not roles:continue
                kind='transporter' if g.startswith('SLC') else 'enzyme_complex_member' if g in ['SDHA','SDHB','SUCLA2','SUCLG1','SUCLG2'] else 'enzyme'
                role='transport' if len(roles)==2 else 'substrate_as_written' if roles==[0] else 'product_as_written'
                ev=json.dumps(rx.get('evidences',[]),ensure_ascii=False)
                rh=';'.join(t['id'] for t in rx.get('reactionCrossReferences',[]) if t['database']=='Rhea')
                add(m,g,kind,role,reaction,'https://www.uniprot.org/uniprotkb/'+x['primaryAccession']+'/entry',accepted=accepted,
                    note=note or '精确反应参与物；反应书写方向不是患者内净通量方向。数据库推断/实验证据代码逐条保留；不宣称均为人类直接酶学实验。',
                    accession=x['primaryAccession'],source_path=rel(CACHE/(g+'.json')),evidence=ev,rhea=rh)
                hits+=1
        return hits
    for m,t in EXACT.items():
        if m in byname: assert match_rxn(m,t)>0,(m,t)
    for m,(t,g,n) in CONDITIONAL.items():
        if m in byname: assert match_rxn(m,t,[g],False,n)>0
    match_rxn('2-hydroxyglutarate','(S)-2-hydroxyglutarate',['L2HGDH'],False,'D/L未区分；仅条件性L构型降解关系。')
    # Selected additional direct relations in already mapped compounds, not
    # global expansion through ATP, phosphate or every pathway member.
    for m,t,gs in [('1-methylnicotinamide','1-methylnicotinamide',['SLC22A2']),('choline','choline',['SLC22A2']),('palmitoylcarnitine','O-hexadecanoyl-(R)-carnitine',['SLC25A20','CPT1A']),('glycerophosphorylcholine (GPC)','sn-glycerol 3-phosphocholine',['LYPLA1','LYPLA2']),('inosine','inosine',['NT5C2'])]:
        match_rxn(m,t,gs)
    for g in ['SDHC','SDHD']:
        x=proteins[g];c=next(c for c in x['comments'] if c['commentType']=='FUNCTION')
        add('succinate',g,'necessary_complex_member','membrane_anchor_not_substrate_binding','Succinate dehydrogenase complex II: succinate / fumarate conversion','https://www.uniprot.org/uniprotkb/'+x['primaryAccession']+'/entry',note='复合体膜锚定成员，不标为独立直接催化琥珀酸的酶。',accession=x['primaryAccession'],source_path=rel(CACHE/(g+'.json')),evidence=json.dumps(c,ensure_ascii=False))
    # Class-specific evidence remains conditional unless the exact measured
    # molecular species is explicitly established.
    for m in ['gamma-glutamylglutamine','gamma-glutamylleucine','gamma-glutamyltyrosine']:
        match_rxn(m,'an alpha-(gamma-L-glutamyl)-L-amino acid',['GGCT'],False,'数据库为gamma-glutamyl氨基酸类反应；本轮未核实该具体二肽的人类底物实验及构型。')
    for m in ['2-palmitoyl-GPC (16:0)','2-palmitoleoyl-GPC (16:1)','2-oleoylglycerophosphocholine','2-arachidonoylglycerophosphocholine']:
        match_rxn(m,'a 2-acyl-sn-glycero-3-phosphocholine',['ENPP2'],False,'2-acyl LPC底物类关系；具体酰基和位置异构体底物证据未完成核实。')
    for m in ['1-stearoylglycerophosphoethanolamine']:
        match_rxn(m,'a 1-acyl-sn-glycero-3-phosphoethanolamine',['LPCAT3','LPCAT4','ENPP2'],False,'1-acyl LPE类底物；具体18:0分子种尚未以底物特异性实验确认。')
    match_rxn('1-palmitoylplasmenylethanolamine','a 1-O-(1Z-alkenyl)-sn-glycero-3-phosphoethanolamine',['LPCAT3','LPCAT4'],False,'plasmenyl醚键与acyl酯键不同；本轮仅烯基醚类底物注释，不以普通LPE关系替代。')
    papers=[
      ('stachydrine','SLC22A4','https://pmc.ncbi.nlm.nih.gov/articles/PMC555966/','人OCTN1表达细胞的stachydrine直接摄取实验；植物来源不排除人类转运关系。',True),
      ('1,5-anhydroglucitol (1,5-AG)','SLC5A10','https://pmc.ncbi.nlm.nih.gov/articles/PMC10439028/','人SGLT5直接运输1,5-AG；采用较新直接测量，不能仅凭抑制SGLT4底物摄取推断SLC5A9运输1,5-AG。',True),
      ('7-alpha-hydroxy-3-oxo-4-cholestenoate (7-Hoca)','AKR1D1','https://www.endocrine-abstracts.org/ea/0086/ea0086oc5.2','作者会议摘要报告新底物；尚未核对完整酶学数据及精确异构体，保留为条件性。',False),
    ]
    for m,g,u,n,a in papers:
        x=proteins[g]
        add(m,g,'transporter' if g.startswith('SLC') else 'enzyme','transport' if g.startswith('SLC') else 'candidate_consumption',n,u,accepted=a,note=n,accession=x['primaryAccession'],source_path='reference/brca_mapping190_review_notes.json',evidence='Original study online text reviewed; not a BRCA functional claim')
    # Explicit per-feature decisions for the remaining named entries.
    notes={
      'C-glycosyltryptophan':'糖基种类与游离/蛋白结合形式未核实；DPY19作用于蛋白色氨酸，不直接移植为游离峰的酶。',
      '1-methyladenosine':'RNA中的m1A修饰酶不等于游离1-甲基腺苷直接代谢酶；本轮未建立精确游离底物关系。',
      'N2,N2-dimethylguanosine':'tRNA修饰不等于游离二甲基鸟苷直接代谢；不直接添加TRMT1。',
      'dimethylarginine (SDMA + ADMA)':'SDMA与ADMA合并信号，不把ADMA特异降解酶套给混合峰。',
      'dihomo-linolenate (20:3n3 or n6)':'n3/n6异构体未区分；不指定单一异构体特异酶。',
      'N-acetylserine':'ACY1底物研究中并非高效底物；泛N-acyl类注释不足以升级具体关系。',
      'N-acetylthreonine':'ACY1底物研究显示效率限制；未取得足以确立本项直接关系的精确证据。',
      'N6-acetyllysine':'侧链N6乙酰化不能替换为ACY1的N-alpha乙酰氨基酸；蛋白去乙酰化不能替代游离底物。',
      'mannitol':'人SORD注释未明确支持mannitol底物；植物/微生物甘露醇酶不移植到人类。',
      'pyroglutamine':'未核实精确结构，不将名称近似自动等同pyroglutamate/5-oxoproline。',
      'tryptophan betaine':'未确认该化合物的人类直接酶或转运体；不以stachydrine运输外推。',
      'phenol red':'外源指示剂条目；保留差异，需先核实实验来源和身份，不作为内源人类代谢酶底物推断。',
      'penicillin G':'外源药物条目；需确认暴露/峰身份，暂不按内源癌代谢关系建基因。',
      '13(S)-HODE (prev. X - 11560)':'ALOX15直接产物为13(S)-HpODE，不等于还原后的HODE；本轮未建立精确一步人类关系。',
      'glycerol 2-phosphate':'2位与sn-3位不同，不沿用甘油-3-磷酸关系；PGP已核对反应未覆盖2位。',
      'succinylcarnitine':'未取得该特定二羧酸酰基肉碱的人类直接反应；不因含succinyl就套用SUCL。',
      'threonate':'命中的降解酶主要来自细菌；不将其作为人类基因关系。',
      'aspartylleucine':'具体二肽水解的人类底物特异性未核实，不泛配所有肽酶。',
      'histidylleucine':'具体二肽水解的人类底物特异性未核实，不泛配所有肽酶。',
      '3-(4-hydroxyphenyl)lactate':'LDH类可疑旁反应不足以确认该具体芳香族羟基酸及构型；不按名称直接映射。',
      '2-hydroxybutyrate (AHB)':'尚未核实精确构型及人类直接底物证据；普通lactate注释不直接外推。',
      '2-aminobutyrate':'未核实氨基酸构型及底物特异性；泛转氨反应不作为充分直接证据。',
    }
    # Deduplicate relation keys while retaining ALL reaction-level evidence.
    for i,e in enumerate(events,1):e['evidence_id']=f'M190E{i:04d}'
    write(out/'relation_evidence.tsv',events)
    edge=[]
    for key in sorted({e['relation_key'] for e in events}):
        ev=[e for e in events if e['relation_key']==key]; accepted=[e for e in ev if e['mapping_status']=='DIRECT_ACCEPTED']; pool=accepted or ev
        e=pool[0]
        edge.append(dict(cancer='BRCA',stage_id='02_MAPPING',run_id=run,analysis_version='mapping190_v2',relation_key=key,
            metabolite_key=e['metabolite_key'],metabolite_name=e['metabolite_name'],gene=e['gene'],
            mapping_status='DIRECT_ACCEPTED' if accepted else 'CONDITIONAL',status='DONE' if accepted else 'NEEDS_REVIEW',
            origins=';'.join(sorted({x['evidence_origin'] for x in pool})),relationship_types=';'.join(sorted({x['relationship_type'] for x in pool})),
            evidence_ids=';'.join(x['evidence_id'] for x in ev),source_urls=' ; '.join(sorted({x['source_url'] for x in pool})),
            gene_in_original117=e['gene'] in oldgenes,patient_analysis_this_round='NOT_RUN',functional_analysis_this_round='NOT_RUN'))
    accepted=[e for e in edge if e['mapping_status']=='DIRECT_ACCEPTED']; conditional=[e for e in edge if e['mapping_status']=='CONDITIONAL']
    write(out/'direct_relations.tsv',accepted);write(out/'direct_relations_excel.csv',accepted,sep=',');write(out/'conditional_relations.tsv',conditional)
    statusrows=[]
    for r in rows:
        a=[e for e in accepted if e['metabolite_key']==r['metabolite_key']];c=[e for e in conditional if e['metabolite_key']==r['metabolite_key']]
        unknown=r['metabolite_name'].startswith('X - ')
        state='DIRECT_MAPPED' if a else 'CONDITIONAL_ONLY' if c else 'UNKNOWN_FEATURE' if unknown else 'NAMED_UNRESOLVED'
        reason='已建立至少一条直接关系；非穷尽映射。' if a else '仅条件性证据，不进入直接基因主池。' if c else '未知特征不猜身份或基因。' if unknown else notes.get(r['metabolite_name'],'已核查相关人类蛋白反应注释，但未取得该具体脂肪酸/脂质分子种、位置或构型的充分直接关系；不以通路/家族成员替代。')
        nr=dict(r);nr.update(mapping_version='mapping190_v2',mapping_state=state,mapping_status='DONE' if a else 'NEEDS_REVIEW',
            mapping_reason=reason,direct_gene_count=len(a),direct_genes=';'.join(sorted(e['gene'] for e in a)),conditional_gene_count=len(c),conditional_genes=';'.join(sorted(e['gene'] for e in c)),
            original_v1_feature_mapped=any(x['metabolite_key']==r['metabolite_key'] for x in old))
        statusrows.append(nr)
    write(out/'metabolite190_mapping_status.tsv',statusrows)
    write(out/'metabolite190_mapping_status_excel.csv',statusrows,sep=',')
    genes=[]
    for g in sorted({e['gene'] for e in accepted}):
        ee=[e for e in accepted if e['gene']==g]
        genes.append(dict(gene=g,in_original117=g in oldgenes,n_direct_relations=len(ee),metabolite_keys=';'.join(e['metabolite_key'] for e in ee),metabolite_names=' ; '.join(e['metabolite_name'] for e in ee),relation_keys=';'.join(e['relation_key'] for e in ee),scope='Direct biochemical candidate; no new patient/function test'))
    write(out/'genes_unique.tsv',genes)
    historical=[]
    for r in old:
        nr=dict(r);nr.update(in_current190=r['metabolite_key'] in bykey,preservation='Original columns retained unchanged; no deletion or statistical overwrite');historical.append(nr)
    write(out/'legacy174_preserved.tsv',historical)
    gc=[]
    currentgenes={g['gene'] for g in genes}
    for g in sorted(oldgenes|currentgenes):
        gc.append(dict(gene=g,in_original117=g in oldgenes,in_current190_direct_pool=g in currentgenes,change='RETAINED_CURRENT' if g in oldgenes&currentgenes else 'NEW_TO_DIRECT_POOL' if g in currentgenes else 'HISTORY_RETAINED_OUTSIDE_CURRENT_POOL',historical_evidence='Retained in original version; not deleted' if g in oldgenes else 'No original117 patient statistics transferred'))
    write(out/'gene_old_new_comparison.tsv',gc)
    if alias_exclusions:write(out/'rejected_gene_synonym_hits.tsv',alias_exclusions)
    unresolved=[r for r in statusrows if r['mapping_state']!='DIRECT_MAPPED'];write(out/'unresolved_and_conditional.tsv',unresolved)
    stats=dict(n_input=190,state_counts=dict(Counter(r['mapping_state'] for r in statusrows)),n_direct_relations=len(accepted),n_conditional_relations=len(conditional),n_genes=len(genes),
        n_old_edges_reused=sum('LEGACY_V1_REUSED' in e['origins'] for e in accepted),n_old_metabolites_reused=sum(r['original_v1_feature_mapped'] for r in statusrows),
        n_new_genes=len(currentgenes-oldgenes),n_old_genes_in_pool=len(currentgenes&oldgenes),n_old_genes_history_only=len(oldgenes-currentgenes),n_evidence_records=len(events),n_unknown=sum(r['mapping_state']=='UNKNOWN_FEATURE' for r in statusrows))
    spec=dict(version='mapping190_v2',scope='Mapping only',input_filter='Frozen existing paired P<0.05 list, 190 features',new_patient_tests=0,new_q_values=0,
        accepted='Existing curated exact-key edges; activated prior direct supplement; explicit human reaction participant or direct original transport experiment; named necessary complex member',
        conditional='Unresolved stereochemistry, generic lipid/peptide class extrapolation, conference abstract only',
        not_exhaustive=True,source_species=9606,primary_gene_symbol_required=True,sort='Input paired directional rank preserved',run_id=run)
    (out/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
    assert all(not e['metabolite_name'].startswith('X - ') for e in accepted)
    assert len({e['relation_key'] for e in accepted})==len(accepted)
    assert all(all(n[k]==r[k] for k in r) for n,r in zip(statusrows,rows))
    assert all(all(n[k]==r[k] for k in r) for n,r in zip(historical,old))
    (out/'validation.json').write_text(json.dumps(dict(summary=stats,input_values_unchanged=True,legacy174_all_original_fields_unchanged=True,unknowns_have_no_direct_genes=True,unique_direct_relation_keys=True,independent_source_reaction_assertions=True,not_validated='No claim of exhaustive mapping, chemical reidentification, clinical or functional validation'),ensure_ascii=False,indent=2),encoding='utf-8')
    src=[P,OLD,SUP,Path(__file__),ROOT/'code/brca_mapping190_fetch.py',ROOT/'reference/brca_mapping190_review_notes.json']+sorted(CACHE.glob('*.json'))
    write(out/'source_manifest.tsv',[dict(path=rel(p),sha256=sha(p),canonical_lf_sha256=canonical_sha(p),scope='Public annotation or aggregate only; canonical LF hash supports Git newline normalization') for p in src])
    (out/'.gitattributes').write_text('* -text\n',encoding='utf-8')
    README=f'''# BRCA 190项配对代谢物直接基因映射 v2

## 本轮问题
把45对主分析P<0.05的190项逐一接到直接人类酶、转运体或必要复合体成员。到映射为止。

## 输入与范围
190项及排序/P/q完全复用；原117基因/174关系历史全部保留。继承CAMP化合物身份，不声称标准品重新鉴定。仅公共注释和已有汇总，无患者矩阵。

## 实际结果
|项目|数量|
|---|---:|
|输入代谢物|190|
|至少一条直接关系|{stats['state_counts'].get('DIRECT_MAPPED',0)}|
|仅条件性关系|{stats['state_counts'].get('CONDITIONAL_ONLY',0)}|
|已命名但本轮未建充分直接关系|{stats['state_counts'].get('NAMED_UNRESOLVED',0)}|
|未知X特征|{stats['n_unknown']}|
|直接关系/去重基因|{len(accepted)} / {len(genes)}|
|复用旧关系/新增或激活关系|{stats['n_old_edges_reused']} / {len(accepted)-stats['n_old_edges_reused']}|
|当前主池中原117基因/新增基因|{stats['n_old_genes_in_pool']} / {stats['n_new_genes']}|

原117中有{stats['n_old_genes_history_only']}个目前不在190入口直接池中，其历史分析仍全部保留；历史加当前的基因并集为{len(oldgenes|currentgenes)}。原始候选与当前入口集合不同，不能把150与117的差简单当成新增基因数。

```json
{json.dumps(stats,ensure_ascii=False,indent=2)}
```
metabolite190_mapping_status.tsv逐项保留原统计与映射去向；direct_relations.tsv按特征键×基因去重；relation_evidence.tsv保留多条反应证据。条件性条目在独立表中，不计入直接基因池。

## 新手解释
一个代谢物能连多个基因，一个基因也能连多个代谢物；因此证据条数、关系条数、基因数不同。酶的反应书写方向不代表患者体内净通量。SDHC/SDHD是必要复合体成员，不是独立底物催化酶。运输关系不是乳腺癌功能验证。

## 限制/反证
未把未知X猜成基因；未把SDMA+ADMA混合信号、2HG构型、泛脂质类底物强制放入主池。未把蛋白/RNA修饰酶当成游离小分子直接酶。SLC5A10的1,5-AG关系采用直接运输原研究，SLC5A9旧抑制实验未直接纳入。UniProt命中检查主基因名，排除PNP、LPCAT4查询的同义词异基因命中。
旧关系来源按既有审定复用，不冒充本轮全文逐项重审；数据库的实验、相似性和推断证据代码保留，不一概称人类实验验证。已映射不等于穷尽；未建立关系不等于人类无此反应。公开单项摘要或来源不足者单独暂挂。

## 当前决定
新版本用于后续候选输入。原117基因和174关系不删除、不改历史P/q；本轮没有新增患者相关、单细胞、功能筛选或排名。

## 下一步
本轮停止在映射交付；未解决项列明所需化学身份或底物证据，待后续授权再开展分析。

## 复现命令
`python code/brca_mapping190_v2.py NEW_UNIQUE_RUN_ID`
复现使用已保存注释快照；只有要更新来源时才执行fetch脚本，不混用刷新后的注释冒充同版本。
'''
    (out/'README_CN.md').write_text(README,encoding='utf-8')
    index=ROOT/'coordination/stages/BRCA.tsv';stages=read(index);fields=list(stages[0])
    vals=['BRCA','02_MAPPING',run,'mapping190_v2','DONE','190 features disposition; bounded direct mapping; conditional/unknown retained',rel(out),'code/brca_mapping190_v2.py','analysis/brca-functional-review-20260919','Mapping only; no patient tests; unresolved features separately retained','Stop at mapping delivery']
    stages.append(dict(zip(fields,vals)));stages.sort(key=lambda r:(r['stage_id'],r['run_id']));write(index,stages,fields)
    print(json.dumps(dict(output=rel(out),**stats),ensure_ascii=False))
if __name__=='__main__':main(sys.argv[1])
