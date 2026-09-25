"""PDAC-specific exact small-molecule/Rhea annotation mapping, with explicit holds."""
import csv,hashlib,json,re,subprocess
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2];CACHE=ROOT/'runtime/pdac_mapping_paired51_v2'
RUN='20260922T055000Z_paired51_mapping_v2'
OUT=ROOT/'results/PDAC/02_MAPPING'/RUN
# Name equivalents preserve position/chain; no class-to-species expansion here.
ALIASES={
'2-linoleoylglycerol (2-monolinolein)':['2-(9Z,12Z-octadecadienoyl)-glycerol'],
'2-oleoylglycerol (2-monoolein)':['2-(9Z-octadecenoyl)-glycerol'],
'2-palmitoylglycerol (16:0)':['2-hexadecanoylglycerol'],
'2-stearoylglycerol (2-monostearin)':['2-octadecanoylglycerol'],
'glycylproline':['glycyl-L-proline'],
"thymidine 3'-monophosphate":["thymidine 3'-phosphate"],
'methionine sulfoxide':['L-methionine (S)-S-oxide','L-methionine (R)-S-oxide'],
'1-stearoylglycerol (18:0)':['1-octadecanoylglycerol']}
HOLDS={
"2'-O-methylguanosine":'Original C04545 denotes modified tRNA, not free nucleoside; do not map tRNA methyltransferases.',
'ribulose/xylulose 5-phosphate':'Combined isomer feature; original C00199 identifies only ribulose; no splitting into independent measurements.',
"adenosine 3'-monophosphate (3'-AMP)":'Historical positional-isomer identity gate unresolved; do not substitute 5-prime AMP.',
'methionine sulfoxide':'Sulfoxide S/R configuration unspecified; keep both hypotheses conditional; free amino acid not protein residue.',
'1-stearoylglycerol (18:0)':'Original D01947 is glyceryl monostearate drug/mixture entry; name-specific 1-isomer hypothesis needs source confirmation.'}
CLASS_CANDIDATES={
'1-pentadecanoylglycerol (1-monopentadecanoin)': ['MGLL','MOGAT1'],
'2-myristoylglycerol (2-monomyristin)': ['MGLL','ABHD6','MOGAT1','MOGAT2'],
'2-stearoylglycerol (2-monostearin)': ['MGLL','ABHD6','MOGAT1','MOGAT2'],
'gamma-glutamylglutamate':['GGT1','GGT5','GGCT'],
'gamma-glutamyltyrosine':['GGT1','GGT5','GGCT'],
'aspartylphenylalanine':['CNDP2','DPEP1','SLC15A1','SLC15A2'],
'glycylleucine':['CNDP2','DPEP1','SLC15A1','SLC15A2'],
'glycylvaline':['CNDP2','DPEP1','SLC15A1','SLC15A2'],
'pro-hydroxy-pro':['PEPD','CNDP2','DPEP1'],
'glycylproline':['SLC15A1','SLC15A2']}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def write(p,rows,fields=None):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
ALIASES.update({
'1-myristoylglycerol (14:0)':['1-tetradecanoylglycerol'],
'1-oleoylglycerol (1-monoolein)':['1-(9Z-octadecenoyl)-glycerol'],
'2-hydroxypalmitate':['2-hydroxyhexadecanoate']})
HOLDS.update({'pipecolate':'D/L unspecified in measured feature; database enantiomer cannot establish measured stereochemistry.',
'C-glycosyltryptophan':'Free C-glycosyltryptophan is not protein-bound C-mannosyltryptophan; no protein modification gene assignment.'})
CLASS_CANDIDATES.update({'glycylisoleucine':['CNDP2','DPEP1','SLC15A1','SLC15A2'],
'1-myristoylglycerol (14:0)':['MGLL','MOGAT1'],
'1,2-dipalmitoylglycerol':['DGAT1','DGAT2','DGKA','DGKB'],
'succinylcarnitine':['CRAT','SLC25A20']})
def main():
    OUT.mkdir(parents=True,exist_ok=False)
    frozen=ROOT/'reference/camp/cancer_effects.tsv'
    assert sha(frozen)=='1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597'
    paired=ROOT/'results/PDAC/04_ROBUSTNESS/20260921T121509Z_paired_metabolites_v1/results.tsv'
    pool=[r for r in read(paired) if r['analysis_type']=='primary' and float(r['p_value'])<.05]
    lookup={r['feature_name'].strip():r for r in read(frozen) if r['cancer']=='PDAC'}
    assert len(lookup)==303 and len(pool)==51
    effects=[]
    for row in pool:
        old=lookup.get(row['metabolite_name'].strip(),{})
        e={k:old.get(k,'NA') for k in next(iter(lookup.values()))}
        e.update(feature_name=row['metabolite_name'],metabolite_key=row['metabolite_key'],paired_p=row['p_value'],paired_q=row['q_value'],paired_n=row['n'],paired_effect=row['effect'],paired_median_delta=row['median_delta'],paired_evidence='NOMINAL_EXPLORATORY')
        effects.append(e)
    names={line.split('\t',1)[0].replace('CHEBI:',''):line.split('\t',1)[1].strip() for line in (CACHE/'chebiId_name.tsv').read_text().splitlines()}
    byname={v:k for k,v in names.items()}
    ph={r['CHEBI']:r['CHEBI_PH7_3'] for r in read(CACHE/'chebi_pH7_3_mapping.tsv')}
    records={}
    for p in CACHE.glob('kegg_*.txt'):
        for block in p.read_text().split('///'):
            d=defaultdict(list);k=None
            for line in block.splitlines():
                if line[:12].strip():k=line[:12].strip()
                if k:d[k].append(line[12:].strip())
            if d.get('ENTRY'):records[d['ENTRY'][0].split()[0]]=d
    proteins=json.loads((CACHE/'human_reviewed.json').read_text())['results'];bygene={};reactions=[]
    for p in proteins:
        assert p['organism']['taxonId']==9606
        gs=[g['geneName']['value'] for g in p.get('genes',[]) if 'geneName' in g]
        for g in gs:bygene[g]=p
        for c in p.get('comments',[]):
            if c['commentType']=='CATALYTIC ACTIVITY':
                rx=c['reaction'];ids={x['id'].replace('CHEBI:','') for x in rx.get('reactionCrossReferences',[]) if x['database']=='ChEBI'}
                if ids:reactions.append((p,gs,rx,ids))
    evidence=[];features=[];conditional=[]
    for e in effects:
        feature=e['feature_name'];rec=records.get(e['kegg_id'],{})
        ids=set()
        if feature not in ALIASES:
            for line in rec.get('DBLINKS',[]):
                if line.startswith('ChEBI:'):ids.update(re.findall(r'\d+',line))
        for name in ALIASES.get(feature,[]):
            assert name in byname,(feature,name)
            ids.add(byname[name])
        ids|={ph.get(i,i) for i in ids}
        # Never turn a modified-tRNA identifier into free-metabolite candidates.
        if feature=="2'-O-methylguanosine":ids=set()
        for p,gs,rx,rids in reactions:
            matched=ids&rids
            if not matched:continue
            rt=rx['name'];is_transport=any(names.get(i,'MISSING')+'(in)' in rt and names.get(i,'MISSING')+'(out)' in rt for i in matched)
            location='; '.join(x['value'] for c in p.get('comments',[]) if c['commentType']=='SUBCELLULAR LOCATION' for sl in c.get('subcellularLocations',[]) for k,x in sl.items() if k=='location') or 'NOT_SPECIFIED'
            for gene in gs:
                function=' '.join(t['value'] for c in p.get('comments',[]) if c['commentType']=='FUNCTION' for t in c.get('texts',[]))
                complex_hold=gene.startswith(('NDUF','MT-ND')) or gene=='GATB' or bool(re.search(r'non-catalytic|accessory subunit|supernumerary subunit',function,re.I))
                experimental=any(x['evidenceCode']=='ECO:0000269' for x in rx.get('evidences',[]))
                modification_context=feature in ['alpha-ketoglutarate','S-adenosylhomocysteine (SAH)'] and bool(re.search(r'RNA|DNA|\[|histone|protein|peptide',rt,re.I))
                hold=HOLDS.get(feature,'') or ('Macromolecule modification reaction; exact cosubstrate/byproduct annotation retained as conditional context, not primary small-molecule-specific axis' if modification_context else '') or ('Complex-level reaction annotation; independent catalytic role requires review' if complex_hold else '') or ('Exact annotation lacks reaction-specific experimental evidence tag; inference-only candidate retained pending specificity review' if not experimental else '')
                evidence.append({'feature_name':feature,'metabolite_key':e['metabolite_key'],'gene':gene,
                 'uniprot_accession':p['primaryAccession'],'matched_chebi':';'.join(sorted(matched)),
                 'matched_participant_names':';'.join(names.get(i,'') for i in sorted(matched)),
                 'reaction':rt,'relation_type':'COMPLEX_MEMBER' if gene in ['MPC1','MPC1L','MPC2','PDHA1','PDHA2','PDHB','LRRC8A','LRRC8B','LRRC8C','LRRC8D','LRRC8E'] else 'DIRECT_TRANSPORT' if is_transport else 'DIRECT_REACTION',
                 'rhea_ids':';'.join(x['id'] for x in rx.get('reactionCrossReferences',[]) if x['database']=='Rhea'),
                 'compartment':location,'annotation_evidence':json.dumps(rx.get('evidences',[])),
                 'human_experimental_annotation':any(x['evidenceCode']=='ECO:0000269' for x in rx.get('evidences',[])),
                 'subunit_notes':json.dumps([c for c in p.get('comments',[]) if c['commentType']=='SUBUNIT']),
                 'source_url':'https://www.uniprot.org/uniprotkb/'+p['primaryAccession']+'/entry',
                 'status':'CONDITIONAL' if hold else 'DIRECT_ANNOTATION','reason':hold or 'Exact small-molecule ID in reviewed human catalytic annotation; not PDAC function or new chemical identification'})
        for gene in CLASS_CANDIDATES.get(feature,[]):
            p=bygene[gene]
            conditional.append({'feature_name':feature,'metabolite_key':e['metabolite_key'],'gene':gene,
                'relation_type':'SUBSTRATE_CLASS_OR_CONTEXT','status':'NEEDS_REVIEW',
                'reason':'Generic substrate-class or adjacent-context annotation; exact measured species specificity not established; not in primary family',
                'source_url':'https://www.uniprot.org/uniprotkb/'+p['primaryAccession']+'/entry'})
        features.append({**e,'database_names':'; '.join(rec.get('NAME',[])),
          'standard_names':'; '.join(names[i] for i in sorted(ids) if i in names),'aliases_used':'; '.join(ALIASES.get(feature,[])),
          'chemical_identity_note':HOLDS.get(feature,'Inherited author chemical identity; name/ID-compatible database mapping, stereochemistry not experimentally revalidated'),
          'identity_status':'HOLD' if feature in HOLDS else 'ANNOTATION_COMPATIBLE' if ids else 'NO_EXACT_ID_BRIDGE',
          'mapping_status':'PENDING','direct_genes':'','n_direct_relations':0,'n_conditional_relations':0})
    # Exact free Gly-Pro support is in the reviewed human PEPD function comment.
    e=next(e for e in effects if e['feature_name']=='glycylproline');p=bygene['PEPD']
    assert 'preferred dipeptide substrate is Gly-Pro' in str(p['comments'])
    evidence.append(dict(feature_name=e['feature_name'],metabolite_key=e['metabolite_key'],gene='PEPD',uniprot_accession=p['primaryAccession'],matched_chebi='73779',matched_participant_names='glycyl-L-proline',reaction='Gly-Pro + H2O -> glycine + proline (exact substrate explicitly named in human function annotation)',relation_type='DIRECT_REACTION',rhea_ids='',compartment='Cytoplasm',annotation_evidence='PubMed:17081196;35165443',human_experimental_annotation=True,subunit_notes='See source',source_url='https://www.uniprot.org/uniprotkb/P12955/entry',status='DIRECT_ANNOTATION',reason='Exact free dipeptide specificity explicitly named; biochemical support only'))
    groups=defaultdict(list)
    for r in evidence:
        if r['status']=='DIRECT_ANNOTATION':groups[(r['feature_name'],r['metabolite_key'],r['gene'])].append(r)
        else:conditional.append({k:r[k] for k in ['feature_name','metabolite_key','gene','relation_type','status','reason','source_url']})
    direct=[]
    for (f,k,g),rr in sorted(groups.items()):
        e=next(e for e in effects if e['feature_name']==f and e['metabolite_key']==k)
        roles=set()
        for r in rr:
            sides=r['reaction'].split(' = ')
            if len(sides)==2:
                for participant in r['matched_participant_names'].split(';'):
                    if not participant:continue
                    pattern=r'(^| \+ )(?:\d+ )?'+re.escape(participant)+r'(?:\(in\)|\(out\))?(?= \+ |$)'
                    if re.search(pattern,sides[0]):roles.add('WRITTEN_LEFT_SUBSTRATE')
                    if re.search(pattern,sides[1]):roles.add('WRITTEN_RIGHT_PRODUCT')
        direct.append({'cancer':'PDAC','cohort':'PDAC','feature_name':f,'metabolite_key':k,'gene':g,
          'relation_id':hashlib.sha256(('PDAC|'+f+'|'+k+'|'+g).encode()).hexdigest()[:16],
          'relation_type':';'.join(sorted({r['relation_type'] for r in rr})),
          'original_effect':e['hedges_g'],'original_q':e['effect_fdr'],
          'compartment':';'.join(sorted({r['compartment'] for r in rr})),
          'reaction_roles':';'.join(sorted(roles)) or 'EXACT_PARTICIPANT_ROLE_NOT_PARSED_SEE_EVIDENCE',
          'source_urls':';'.join(sorted({r['source_url'] for r in rr})),
          'rhea_ids':';'.join(sorted({r['rhea_ids'] for r in rr if r['rhea_ids']})),
          'has_human_experimental_annotation':any(r['human_experimental_annotation'] for r in rr),
          'identity_note':next(x['chemical_identity_note'] for x in features if x['feature_name']==f),
          'cosubstrate_context':f in ['alpha-ketoglutarate','S-adenosylhomocysteine (SAH)'],
          'mapping_version':'PDAC_paired51_mapping_v2','status':'LOCKED_DIRECT_ANNOTATION'})
    for f in features:
        d=[r for r in direct if r['feature_name']==f['feature_name']];c=[r for r in conditional if r['feature_name']==f['feature_name']]
        f.update(direct_genes=';'.join(r['gene'] for r in d),n_direct_relations=len(d),n_conditional_relations=len({r['gene'] for r in c}),
          mapping_status='IDENTITY_REVIEW_REQUIRED' if f['identity_status']=='HOLD' else 'DIRECT_RELATION_FOUND' if d else 'CONDITIONAL_ONLY' if c else 'NO_RELIABLE_EXACT_RELATION_FOUND_IN_SCOPE')
    # Deduplicate conditional support without dropping reaction provenance evidence.
    conditional=list({tuple(r.values()):r for r in conditional}.values())
    write(OUT/'feature_coverage_51.tsv',features);write(OUT/'direct_relations_v2.tsv',direct)
    write(OUT/'conditional_relations_v2.tsv',conditional);write(OUT/'mapping_evidence.tsv',evidence)
    aliases={}
    for g in sorted({r['gene'] for r in direct+conditional}):
        aliases[g]=sorted({s['value'] for entry in bygene[g].get('genes',[]) if entry.get('geneName',{}).get('value')==g for s in entry.get('synonyms',[])})
    dump(OUT/'gene_aliases.json',aliases)
    sources=json.loads((CACHE/'sources.json').read_text());sources_list=list(sources.values())+json.loads((CACHE/'new_sources.json').read_text());write(OUT/'mapping_sources.tsv',[{k:x.get(k,'NA') for k in ['path','url','sha256','retrieved_utc']} for x in sources_list])
    spec={'version':'PDAC_paired51_mapping_v2','planned_M':len(direct),'planned_G':len({r['gene'] for r in direct}),
      'locked_before_patient_outcomes':True,'source_proteins':len(proteins),'mapping_rule':'Exact KEGG/explicit chemical-name-to-ChEBI bridge, pH mapping, reviewed human catalytic participant AND reaction-specific ECO:0000269 annotation (or explicit experimentally documented exact substrate). Inference-only records retained conditional. No class-to-species inference in primary family.',
      'NAD_rule':'Include exact NAD+ participant including redox/cofactor reactions; flag separately; not NAD-specific target claim. Complex/accessory subunits held.',
      'BH':'Fixed M planned relationships, separately primary/availability; internal p=1 placeholders, public uncomputable P/q NA',
      'permutations':9999,'bootstrap':4000,'minimum_n':8,'new_imputation':False,'RNA_expression':'Author-confirmed 11 pairs;21 author-key tumor units for correlation;all candidate genes retained',
      'historical_reuse':'Require all expected source hashes to match; same finite input arrays+method reuse rho/P/CI; recalculate new-family q',
      'code_base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'script_sha256':sha(Path(__file__)),
      'limitations':['Annotation mapping, not full experimental specificity certification','Original author chemistry not reidentified','No causal claims from annotation','Conditional identities excluded from primary family']}
    dump(OUT/'analysis_spec.json',spec)
    dump(OUT/'summary.json',{'features':len(features),'M':len(direct),'G':spec['planned_G'],'feature_status_counts':dict(Counter(f['mapping_status'] for f in features)),'conditional_unique_pairs':len({(r['feature_name'],r['gene']) for r in conditional})})
    assert len(features)==51 and len(direct)==len({(r['feature_name'],r['metabolite_key'],r['gene']) for r in direct})
    assert not any(r['feature_name'] in HOLDS for r in direct)
    dump(OUT/'validation.json',{'status':'PASS','all51_accounted':True,'unique_relation_keys':True,'identity_holds_excluded':True,'original_fields_retained':True,'patient_statistics':'NOT_RUN'})
    write(OUT/'source_manifest.tsv',[{'path':'reference/camp/cancer_effects.tsv','sha256':sha(frozen)},{'path':str(paired.relative_to(ROOT)),'sha256':sha(paired)},{'path':'code/pdac/build_mapping_paired51_v2.py','sha256':sha(Path(__file__))}])
    planned=[]
    for r in direct:
        planned.append({**r,'mapping_tier':'DIRECT','mapping_reason':'Exact human experimental annotation; source chemistry inherited'})
    cg=defaultdict(list);dk={(r['feature_name'],r['metabolite_key'],r['gene']) for r in direct}
    for r in conditional:
        key=(r['feature_name'],r['metabolite_key'],r['gene'])
        if key not in dk:cg[key].append(r)
    for (f,k,g),rr in sorted(cg.items()):
        e=next(x for x in effects if x['feature_name']==f)
        row={key:'NA' for key in planned[0]}
        row.update(cancer='PDAC',cohort='PDAC',feature_name=f,metabolite_key=k,gene=g,relation_id=hashlib.sha256(('PDAC|'+f+'|'+k+'|'+g).encode()).hexdigest()[:16],mapping_tier='CONDITIONAL',mapping_reason='; '.join(sorted({x['reason'] for x in rr})),relation_type=';'.join(sorted({x['relation_type'] for x in rr})),source_urls=';'.join(sorted({x['source_url'] for x in rr})),original_effect=e['hedges_g'],original_q=e['effect_fdr'],mapping_version='PDAC_paired51_mapping_v2',status='CONDITIONAL')
        planned.append(row)
    write(OUT/'planned_relations.tsv',planned)
    assert len(planned)==len({(x['feature_name'],x['metabolite_key'],x['gene']) for x in planned})
    dump(OUT/'panel.json',{'genes':sorted({x['gene'] for x in planned}),'relations':len(planned),'direct':len(direct),'conditional':len(planned)-len(direct),'features':51})
    print((OUT/'summary.json').read_text());print((OUT/'panel.json').read_text())
if __name__=='__main__':main()
