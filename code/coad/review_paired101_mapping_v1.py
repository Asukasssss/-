"""Scoped paired101 mapping: exact human reaction/chemical annotation and transport.

Does not infer mechanism from expression or silently replace original identifiers.
"""
import csv,io,json,re,urllib.parse
from collections import defaultdict
from pathlib import Path
import extend_paired101_evidence_v1 as src

ROOT,WORK=src.ROOT,src.WORK
POLICIES={
 'serine':('HOLD_D_L_ID_CONFLICT',['C00065','C00740'],'Original KEGG L-serine vs HMDB0003406 D-serine; https://hmdb.ca/metabolites/HMDB0003406'),
 '2-hydroxyglutarate':('HOLD_STEREOCHEMISTRY',['C02630','C01087','C03196'],'Unresolved measured stereoisomer; R/S are hypotheses, not separate measurements.'),
 'fructose 1,6-diphosphate/glucose 1,6-diphosphate/myo-inositol diphosphates':('HOLD_COMBINED_FEATURE',['C00354'],'Only original fructose-1,6-bisphosphate hypothesis screened here; other combined components unresolved.'),
 'glutarate (C5-DC)':('HOLD_NAME_CLASS_CONFLICT',['C00489'],'Glutarate annotation and C5-DC acylcarnitine-style label require original identity confirmation; not interchangeable.'),
 '3-methylhistidine':('HOLD_METHYL_POSITION_NOMENCLATURE',['C01152'],'Source feature and KEGG historical 1/3-methyl synonyms require position-specific annotation check; no histone residue enzyme jump.'),
 '2-Hydroxypentanoate':('HOLD_NO_EXACT_DATABASE_ID',[],'Original KEGG/HMDB empty; scoped KEGG name search returned no entry. Not negative biological evidence.')
}


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    src.setup();ev=json.loads((WORK/'evidence.json').read_text())
    paired=read(ROOT/'results/COAD/03_PATIENT/20260921T123200Z_paired_metabolites33_v1/primary33.tsv')
    selected={r['metabolite_name'] for r in paired if float(r['p_value'])<.05}
    effects=[r for r in read(ROOT/'reference/camp/cancer_effects.tsv') if r['cancer']=='COAD' and r['feature_name'] in selected]
    oldident={r['feature_name']:r for r in json.loads((ROOT/'reports/COAD/fresh_mapping_v0_1/identity_review.json').read_text())}
    identities=[]
    for r in effects:
        if r['feature_name'] in oldident:i=dict(oldident[r['feature_name']],identity_evidence_origin='REUSED_COAD_20260919')
        else:
            status,ids,note=POLICIES.get(r['feature_name'],('ORIGINAL_KEGG_NAME_COMPATIBLE',[r['kegg_id']],'Original identifier/name compatible; no independent chemical authentication or full HMDB stereochemistry crosscheck.'))
            i=dict(r,identity_review_status=status,compound_hypotheses=ids,identity_note=note,identity_evidence_origin='NEW_PAIRED_POOL_20260922')
        identities.append(i)
    extra={c for i in identities for c in i['compound_hypotheses']}-set(ev['compound_records'])
    ev['compound_records'].update(src.kegg.get_entries(extra))
    # Supplemental stereoisomer hypotheses are always conditional.
    reg=json.loads((src.OLD/'sources.json').read_text());item=reg['https://rest.kegg.jp/link/ko/reaction']
    allrk=defaultdict(list)
    for line in (src.OLD/item['cache_file']).read_text().splitlines():
        rid,ko=line.split('\t')
        if ko in ev['human_ko_genes']:allrk[rid.removeprefix('rn:')].append(ko)
    needed=set()
    for c in extra:needed.update(re.findall(r'R\d{5}',' '.join(ev['compound_records'].get(c,{}).get('REACTION',[]))))
    ev['reaction_records'].update(src.kegg.get_entries(needed&set(allrk)-set(ev['reaction_records'])))
    for rid in needed&set(allrk):ev['reaction_human_kos'][rid]=allrk[rid]
    chains=[]
    for i in identities:
        for c in i['compound_hypotheses']:
            for rid in re.findall(r'R\d{5}',' '.join(ev['compound_records'].get(c,{}).get('REACTION',[]))):
                rr=ev['reaction_records'].get(rid)
                if not rr:continue
                eq=' '.join(rr.get('EQUATION',[]));sides=re.split(r'\s*(?:<=>|=>|<=)\s*',eq)
                if len(sides)!=2:continue
                tokens=[set(re.findall(r'\bC\d{5}\b',s)) for s in sides]
                if c not in tokens[0]|tokens[1]:continue
                for ko in ev['reaction_human_kos'].get(rid,[]):
                    assert ko.split(':')[1] in re.findall(r'K\d{5}',' '.join(rr.get('ORTHOLOGY',[])))
                    for gene in ev['human_ko_genes'][ko]:
                        desc=ev['human_gene_names'].get(gene,'');symbol=desc.split(';')[0].split(',')[0].strip() if ';' in desc else ''
                        chains.append(dict(feature_name=i['feature_name'],metabolite_key=i['metabolite_key'],human_gene_id=gene,gene_symbol=symbol,
                            compound_hypothesis=c,reaction_id=rid,ko_id=ko,equation=eq,written_equation_side='BOTH' if c in tokens[0]&tokens[1] else 'LEFT' if c in tokens[0] else 'RIGHT'))
    proteins=json.loads((ROOT/'runtime/coad_uniprot_rhea_20260919/proteins.json').read_text())
    params={'query':'organism_id:9606 AND reviewed:true AND keyword:KW-0813','format':'json','fields':'accession,gene_names,organism_id,xref_geneid,cc_catalytic_activity,cc_subunit'}
    transport_proteins=json.loads(src.human.fetch('https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode(params)))['results']
    for p in transport_proteins:proteins[p['primaryAccession']]=p
    present={ 'hsa:'+x['id'] for p in proteins.values() for x in p.get('uniProtKBCrossReferences',[]) if x['database']=='GeneID'}
    missing=sorted({r['gene_symbol'] for r in chains if r['human_gene_id'] not in present and r['gene_symbol']})
    for start in range(0,len(missing),25):
        batch=missing[start:start+25];query='organism_id:9606 AND reviewed:true AND ('+' OR '.join('gene_exact:'+s for s in batch)+')'
        params={'query':query,'format':'json','size':500,'fields':'accession,gene_names,organism_id,xref_geneid,cc_catalytic_activity,cc_subunit'}
        for p in json.loads(src.human.fetch('https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode(params)))['results']:proteins[p['primaryAccession']]=p
        print('Additional human genes queried',min(start+25,len(missing)),len(missing),flush=True)
    # Supplement incomplete KEGG->KO paths using reviewed human catalytic records.
    # Search words only locate records; an exact chemical reaction match is still required below.
    for term in ['hydroxyglutarate','hydroxybutyrate','sebacic','nonanoate','piperidine','hydroxypentanoate','CNP','malonate']:
        params={'query':'organism_id:9606 AND reviewed:true AND ('+term+')','format':'json','size':500,
                'fields':'accession,gene_names,organism_id,xref_geneid,cc_catalytic_activity,cc_subunit'}
        for p in json.loads(src.human.fetch('https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode(params)))['results']:
            proteins[p['primaryAccession']]=p
    by_gene=defaultdict(list)
    for p in proteins.values():
        assert p['organism']['taxonId']==9606 and p['entryType']=='UniProtKB reviewed (Swiss-Prot)'
        for x in p.get('uniProtKBCrossReferences',[]):
            if x['database']=='GeneID':by_gene['hsa:'+x['id']].append(p)
    master={}
    for row in csv.DictReader(io.StringIO(src.human.fetch('https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv')),delimiter='\t'):
        for x in row.values():master[x]=row['RHEA_ID_MASTER']
    kr=defaultdict(set)
    for row in csv.DictReader(io.StringIO(src.human.fetch('https://ftp.expasy.org/databases/rhea/tsv/rhea2kegg_reaction.tsv')),delimiter='\t'):kr[row['ID']].add(row['MASTER_ID'])
    for rid,rr in ev['reaction_records'].items():
        for l in rr.get('DBLINKS',[]):
            if l.startswith('RHEA:'):kr[rid].update(master.get(x,x) for x in re.findall(r'\d+',l))
    ph={r['CHEBI']:r for r in csv.DictReader(io.StringIO(src.human.fetch('https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv')),delimiter='\t')}
    cn={}
    for l in src.human.fetch('https://ftp.expasy.org/databases/rhea/tsv/chebiId_name.tsv').splitlines():
        k,v=l.split('\t',1);cn[k.removeprefix('CHEBI:')]=v.strip()
    allowed={}
    for c,rr in ev['compound_records'].items():
        ids={x for l in rr.get('DBLINKS',[]) if l.startswith('ChEBI:') for x in re.findall(r'\d+',l)}
        allowed[c]=ids|{ph.get(x,{}).get('CHEBI_PH7_3',x) for x in ids}
    evidence=[]
    for r in chains:
        for p in by_gene[r['human_gene_id']]:
            for comment in p.get('comments',[]):
                if comment['commentType']!='CATALYTIC ACTIVITY':continue
                reaction=comment['reaction'];refs=reaction.get('reactionCrossReferences',[])
                mids={master.get(x['id'].removeprefix('RHEA:'),x['id'].removeprefix('RHEA:')) for x in refs if x['database']=='Rhea'}
                if not (mids&kr[r['reaction_id']]):continue
                chebi={x['id'].removeprefix('CHEBI:') for x in refs if x['database']=='ChEBI'}
                evidence.append(dict(r,kind='ENZYME_REACTION',uniprot_accession=p['primaryAccession'],reaction=reaction['name'],
                    rhea_ids=sorted(mids&kr[r['reaction_id']]),target_chebi_matches=sorted(chebi&allowed[r['compound_hypothesis']]),
                    annotation_evidence=reaction.get('evidences',[]),subunit_notes=[x for x in p.get('comments',[]) if x['commentType']=='SUBUNIT']))
    # Direct ChEBI-anchored reactions can be valid even without a KEGG KO path.
    for i in identities:
        for c in i['compound_hypotheses']:
            for p in proteins.values():
                for comment in p.get('comments',[]):
                    if comment['commentType']!='CATALYTIC ACTIVITY':continue
                    rx=comment['reaction'];refs=rx.get('reactionCrossReferences',[])
                    cs={x['id'].removeprefix('CHEBI:') for x in refs if x['database']=='ChEBI'}&allowed.get(c,set())
                    rheas=[x['id'] for x in refs if x['database']=='Rhea']
                    if not cs or not rheas or '(in)' in rx['name'] or '(out)' in rx['name']:continue
                    for x in p.get('uniProtKBCrossReferences',[]):
                        if x['database']!='GeneID':continue
                        evidence.append(dict(feature_name=i['feature_name'],metabolite_key=i['metabolite_key'],human_gene_id='hsa:'+x['id'],
                            gene_symbol=';'.join(g['geneName']['value'] for g in p.get('genes',[]) if 'geneName' in g),
                            compound_hypothesis=c,kind='DIRECT_RHEA_CHEBI',uniprot_accession=p['primaryAccession'],reaction=rx['name'],
                            rhea_ids=rheas,target_chebi_matches=sorted(cs),annotation_evidence=rx.get('evidences',[]),
                            subunit_notes=[x for x in p.get('comments',[]) if x['commentType']=='SUBUNIT']))
    for i in identities:
        for c in i['compound_hypotheses']:
            for p in transport_proteins:
                for comment in p.get('comments',[]):
                    if comment['commentType']!='CATALYTIC ACTIVITY':continue
                    rx=comment['reaction'];text=rx['name'];refs=rx.get('reactionCrossReferences',[])
                    if '(in)' not in text or '(out)' not in text:continue
                    cs={x['id'].removeprefix('CHEBI:') for x in refs if x['database']=='ChEBI'}&allowed.get(c,set())
                    cs={c for c in cs if c in cn and cn[c]+'(in)' in text and cn[c]+'(out)' in text}
                    if not cs:continue
                    genes=[x['id'] for x in p.get('uniProtKBCrossReferences',[]) if x['database']=='GeneID']
                    symbols=';'.join(g['geneName']['value'] for g in p.get('genes',[]) if 'geneName' in g)
                    for gene in genes:evidence.append(dict(feature_name=i['feature_name'],metabolite_key=i['metabolite_key'],human_gene_id='hsa:'+gene,gene_symbol=symbols,
                        compound_hypothesis=c,kind='TRANSPORT_REACTION',uniprot_accession=p['primaryAccession'],reaction=text,
                        rhea_ids=[x['id'] for x in refs if x['database']=='Rhea'],target_chebi_matches=sorted(cs),annotation_evidence=rx.get('evidences',[]),
                        subunit_notes=[x for x in p.get('comments',[]) if x['commentType']=='SUBUNIT']))
    for name,obj in [('identity_review.json',identities),('reaction_chains.json',chains),('human_evidence.json',evidence),('evidence_extended.json',ev),('proteins_extended.json',proteins)]:
        (WORK/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(features=len(identities),chains=len(chains),human_evidence=len(evidence),new_gene_queries=len(missing))),flush=True)


if __name__=='__main__':main()
