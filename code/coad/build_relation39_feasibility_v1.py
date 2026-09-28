"""Deterministic relationship-gap delta; preserve existing39 via references, not another copy."""
import csv,json,hashlib
from pathlib import Path
import pandas as pd
from cell_origin_v1 import PREFIX
R=Path(__file__).resolve().parents[2]
SPEC=R/'config/coad_relation39_feasibility_v1.json'
S=json.loads(SPEC.read_text(encoding='utf-8'));RUN=S['run_id']
OUT=R/'results/COAD/07_INTEGRATION'/RUN

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(rows,name):pd.DataFrame(rows).fillna('NA').to_csv(OUT/name,sep='\t',index=False)

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    with (R/S['source_table']).open(encoding='utf-8',newline='')as h:r=csv.DictReader(h,delimiter='\t');original=list(r)
    assert len(original)==39 and len({r['gene']for r in original})==35
    dataset_fields=['source_id','title','url','biological_material','reported_unit','co_measurement','independence','current_decision','reason','next_required']
    datasets=[
    ['D01','Chen2025 PLCD4 multiomics','https://link.springer.com/article/10.1186/s12943-025-02359-x','CRC primary tissue; colon/rectal distinction unverified','13 RNA patients;11 protein/metabolite patients;3 technical metabolome replicates per tissue','Reported matched;actual cross-omics sample keys unverified','Separate institution/publication;patient overlap still needs metadata','PARTIAL','Promising but no usable matched tables verified;data-availability statement conflicts with reported generated data','Processed full feature matrix, RNA/protein table, explicit11-patient pairing/replicate map, colon site and provenance'],
    ['D02','Morris2023 colon liver metastasis','https://pmc.ncbi.nlm.nih.gov/articles/PMC10620771/','Colon liver metastasis and normal liver','90 metabolomic patients;RNA overlapping subset28','Overlap reported;exact specimen join not verified','Potential independent source;metadata unverified','PARTIAL','Metastatic liver context and prior chemotherapy;not primary colon validation','Feature list, paired matrix IDs, treatment metadata;separate metastatic validation design'],
    ['D03','ColoCare adipose multiomics','https://pmc.ncbi.nlm.nih.gov/articles/PMC4515859/','Visceral/subcutaneous adipose and serum','59 CRC patients in metabolomic cohort','Adipose transcriptome/metabolome reported','Distinct material;patient overlap not audited','NOT_EVALUABLE','Wrong tissue for within-colon-tumor relation','Only consider a separately specified adipose hypothesis'],
    ['D04','CPTAC-COAD listed collection','https://www.cancerimagingarchive.net/collection/cptac-coad/','Colon tumor;proteogenomics/imaging','106 subjects on listed imaging collection','Listed genomics/proteomics;matching small-molecule metabolome not established','Not audited for this question','NOT_EVALUABLE','No matching metabolite measurements established at inspected entry;not claim no metabolomics could exist elsewhere','Find actual linked metabolomics study and exact specimen map before use'],
    ['D05','GSE89076 original CRC multiomics','https://pmc.ncbi.nlm.nih.gov/articles/PMC5604037/','Original CRC tumor/normal source','275 in full publication;existing RNA subset already used','Original study supports existing analysis','NOT_INDEPENDENT:original patient analysis source','NOT_EVALUABLE','Must not relabel original data as independent external validation','Retain as original CAMP/patient evidence only'],
    ['D06','Early-onset CRC metabolomics/transcriptomics2026','https://link.springer.com/article/10.1007/s12672-026-05185-9','CRC tissue metabolomics plus external GEO transcriptomics','45 metabolomic tissues;GSE39084 transcriptomic cohort','Different cohorts','No same-subject multimodal pairing','NOT_EVALUABLE','Pathway integration across different patients cannot test individual metabolite-RNA correlation','Obtain truly matched cohort;do not concatenate unrelated samples'],
    ['D07','PXD046882 branched-chain amino acid study','https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD046882','Deposited proteome mouse;human metabolite validation described','Mouse AOM/DSS plus human validation;exact matched human n unverified','Human valine/leucine/isoleucine measurements reported;matching full human RNA/protein not established','Mixed species/design;not established','NOT_EVALUABLE','Mouse proteomics cannot be joined to human metabolites;BCAA presence alone insufficient','Human same-specimen gene and metabolite tables;species and patient keys'],
    ['D08','Patient C CRC PDX trametinib multiomics','https://zenodo.org/records/18471666','Treated patient-derived xenograft','One patient-derived model in entry title','Multiomics model study','Model replicates not independent patients','NOT_EVALUABLE','Potential model resource,not independent human tumor cohort','Separate treatment/model analysis;verify actual feature coverage'],
    ['D09','Sex differences in colon cancer metabolism2020','https://pmc.ncbi.nlm.nih.gov/articles/PMC7078199/','Colon tissue metabolomics with external TCGA expression','Independent transcriptomic validation cohort','Same-subject RNA/metabolite matching not established','TCGA expression is separate','NOT_EVALUABLE','Expression validation of metabolomic pathway is not matched metabolite-RNA validation','Actual matched RNA/protein table and specimen mapping required']]
    save([dict(zip(dataset_fields,x))for x in datasets],'dataset_feasibility.tsv')
    efields=['evidence_id','gene','metabolite_key','url','new_to_project','review_depth','specific_finding','metabolite_perturbation','phenotype_mediation','cell_context_limit','action']
    evidence=[
    ['F01','AQP9','KEGG:C00147','https://journals.physiology.org/doi/10.1152/ajprenal.1999.277.5.F685',True,'Original publisher indexed results and Fig3/method extracts','Human AQP9 expression in Xenopus oocytes increases adenine permeability;labelled substrate uptake','Direct adenine transport assay in heterologous system','CRC adenine-dependent phenotype NOT_SHOWN','Oocyte assay,not CRC myeloid cell','Test adenine uptake after exact-gene manipulation in relevant myeloid/CRC model'],
    ['F02','ENO2','KEGG:C00631','https://pmc.ncbi.nlm.nih.gov/articles/PMC9367517/',False,'New focused reread of existing E10 result3.3/Fig3','RKO/SW480/DLD1 experiments support CYTOR-related phenotype without altered glycolytic readouts;substrate-binding mutant retains migration effect','Glucose/lactate/ECAR measured;exact2-PG mediation not established','Counterevidence to assigning reported phenotype to2-PG','Epithelial models,not observed myeloid compartment','Separate catalytic2-PG hypothesis from noncatalytic CYTOR phenotype'],
    ['F03','ENO2','KEGG:C00631','https://pmc.ncbi.nlm.nih.gov/articles/PMC11080076/',True,'Original abstract/results extracts','MSI-H context reports glycolytic and invasion effects after ENO2 manipulation','Glycolytic phenotype,not exact2-PG measurement verified','Exact2-PG rescue NOT_VERIFIED','Different subtype/model from E10;not universal contradiction','Retain genotype/model-specific evidence;verify exact metabolite'],
    ['F04','ASPA','KEGG:C01042','https://www.nature.com/articles/s41467-026-73002-6',True,'Original results/Fig4/discussion/method extracts and PubMed identity','Fibroblast ASPA loss activates TGFβ-related phenotype;E178D retains tested suppressive effects','NAA measurement/supplementation reported','NAA supplementation ineffective for reported phenotype;supports catalytic-independent function','CRC stromal observations;key manipulations include mammary fibroblast and other cancer contexts,not CRC-specific genetic proof','ASPA function record updated;do not equate NAA association with functional mediation'],
    ['F05','B4GALT2','KEGG:C00015','https://www.nature.com/articles/s42003-026-10017-1',True,'Original indexed results Fig6/7 and supplementary Fig1/2 descriptions;not independent source-data audit','CRC B4GALT2 interference affects phenotype and MUC20 glycosylation/protein;THP1 co-culture described','MUC20 glycosylation assay;UDP pool/flux not verified','Exact UDP-dependent rescue NOT_VERIFIED','HCT116/HT29 epithelial manipulation differs from bulk stromal expression label','New CRC functional evidence requires hold-reason review;keep historical9/21/5 until separate decision version'],
    ['F06','ICMT','KEGG:C00019','https://pmc.ncbi.nlm.nih.gov/articles/PMC555472/',False,'New focused reread of existing E19 results','Cysmethynil target-dependent anchorage-independent growth phenotype reversed by ICMT overexpression','Methylation readout;bulk SAM perturbation/flux not established','Target rescue is not SAM rescue','Colon cancer pharmacology,not stromal genetic evidence','Require SAM/SAH plus prenyl-protein methylation measures;separate target from metabolite rescue'],
    ['F07','ICMT','KEGG:C00019','https://pubmed.ncbi.nlm.nih.gov/42127111/',True,'Original PubMed abstract','ICMT-INPP5E methylation/localization rescue in BRAFV600E melanoma','PI(4,5)P2 change reported;not exact SAM relation','Not CRC SAM-mediated validation','Melanoma,not COAD','Keep as excluded-for-CRC contextual mechanism;do not relabel as CRC evidence']]
    save([dict(zip(efields,x))for x in evidence],'focused_function_evidence.tsv')
    types={'NNMT':'EXACT_METABOLITE_FOLLOWUP','SLC6A6':'EXACT_METABOLITE_FOLLOWUP','HDC':'EXACT_METABOLITE_FOLLOWUP','BCAT2':'CONDITIONAL_BCAA_FOLLOWUP','AQP9':'TRANSPORT_TO_CONTEXT_BRIDGE','ASPA':'SEPARATE_METABOLIC_AND_NONCATALYTIC','ENO2':'SEPARATE_METABOLIC_AND_NONCATALYTIC','UCKL1':'SEPARATE_METABOLIC_AND_NONCATALYTIC','B4GALT2':'GENE_FUNCTION_UPDATE_METABOLITE_GAP','GSTT2':'MEASUREMENT_IDENTITY_FIRST','GSTT2B':'MEASUREMENT_IDENTITY_FIRST','PRMT7':'ISOFORM_AND_MEDIATOR_GAP'}
    actions={
    'NNMT':'For1-MNA quantify exact product after NNMT perturbation and physiological-dose mediation;for SAM separately measure SAM/SAH/methyl-donor effects;do not transfer1-MNA rescue to SAM',
    'SLC6A6':'Use established taurine uptake assay as starting evidence;match treatment/cell state and determine whether transport change mediates phenotype',
    'HDC':'Separate host-myeloid histamine intervention from tumor mast-cell source;cell-specific histamine production/response required',
    'BCAT2':'Resolve individual isoleucine versus pooled BCAA readout;retain nutrition/circRNA conditions and normal-diet negative result',
    'AQP9':'Bridge heterologous adenine permeability to relevant CRC/myeloid uptake and metabolite-dependent phenotype;5-FU effects do not substitute',
    'ENO2':'First confirm2-PG coverage/identity and patient subset sensitivity;compare catalytic versus CYTOR-related function in specified model',
    'ASPA':'Measure NAA turnover separately from E178D/TGFβ noncatalytic effects;CRC stromal genetic/metabolite mediation remains missing',
    'B4GALT2':'Reassess functional hold with2026CRC interference evidence;measure UDP versus UDP-galactose and glycosylation endpoints separately',
    'ICMT':'Distinguish gene/drug target rescue from SAM-mediated rescue;measure exact SAM/SAH and relevant protein methylation in applicable cell context',
    'UCKL1':'Require independent evidence for uridine mediation;existing noncanonical NRF2/ferroptosis function cannot fill that gap',
    'GSTT2':'Resolve gene-specific RNA/protein measurement and GSTT2B cross-assignment before interpreting GSH relation',
    'GSTT2B':'Do not borrow GSTT2 genetic/function evidence;resolve gene-specific measurement and exact GSH phenotype mediation',
    'PRMT7':'Gene-level SAM association remains separate from splice-event function;verify actual isoform measurement and substrate methylation',
    'UPP2':'Keep8-positive/20-zero external-expression finding;verify UPP2-specific uridine turnover without borrowing UPP1',
    'SLC38A3':'Low single-cell detection and missing asparagine-specific uptake/rescue require independent measurement before source-specific mechanism',
    'BPGM':'2-PG usable patient subset n11 is fragile;confirm feature identity and relation before functional escalation',
    'PGAM2':'2-PG coverage and PGAM2 versus PGAM1 identity first;no PGAM1 function substitution'}
    rows=[];matrix=[]
    for i,old in enumerate(original,1):
        gene=old['gene'];key=old['metabolite_key'];ids=[x[0]for x in evidence if x[1]==gene and x[2]==key]
        row=dict.fromkeys(PREFIX,'NA');row.update(cancer='COAD',cohort='existing39_review',stage_id='07_INTEGRATION',run_id=RUN,analysis_version=S['analysis_version'],analysis_type='relation_specific_feasibility_review',
            metabolite_key=key,metabolite_name=old['metabolite_name'],gene=gene,status='PARTIAL',reason='No independently verified paired metabolite-gene dataset ready;not a negative association',source_id=S['source_table'])
        route=types.get(gene,'MEASURE_EXACT_METABOLITE_AND_MEDIATION')
        if gene=='NNMT' and key=='KEGG:C00019':route='METHYL_DONOR_SPECIFIC_GAP'
        row.update(source_row_1based=i,review_route=route,focus5=gene in S['focus_genes'],new_review_ids=';'.join(ids)or'NA',
            inherited_evidence_ids=old['func35_evidence_ids'],inherited_review_note='Reuse prior recorded depth;not fresh exhaustive full-text review',
            inherited_metabolic_link=old['func35_metabolic_link'],inherited_counterevidence=old['func35_counterevidence_or_limit'],
            independent_patient_association='NOT_RUN',actual_metabolite_coverage='NOT_VERIFIED',actual_RNA_or_protein_coverage='NOT_VERIFIED',paired_specimen_identity='NOT_VERIFIED',
            exact_metabolite_identity_requirement=key+'; exact annotated feature required, no related-metabolite substitution',
            functional_decision=actions.get(gene,'Require '+old['metabolite_name']+' readout after exact '+gene+' manipulation, metabolite mediation/rescue, and compatible cell/compartment;retain prior limitations'),
            existing_arrangement='UNCHANGED_9_21_5',new_arrangement_review='NEEDS_REVIEW_NEW_FUNCTIONAL_EVIDENCE' if gene=='B4GALT2'else'NA',
            next_dataset_action='D01 processed tables/replicate map;D02 only under separate metastatic scope',proof_status='NOT_VALIDATED_METABOLITE_AXIS')
        rows.append(row)
        for ds in datasets:matrix.append(dict(cancer='COAD',source_row_1based=i,metabolite_key=key,gene=gene,dataset_id=ds[0],dataset_decision=ds[7],relation_ready=False,feature_presence='NOT_VERIFIED',gene_presence='NOT_VERIFIED',reason=ds[8]))
    save(rows,'relation39_next_evidence.tsv');save(matrix,'relation_dataset_checks.tsv')
    features=[]
    for key,grp in pd.DataFrame(original).groupby('metabolite_key',sort=True):
        features.append(dict(metabolite_key=key,metabolite_name=grp.iloc[0].metabolite_name,genes=';'.join(sorted(grp.gene.unique())),n_relations=len(grp),
            requested_fields='Original feature ID,structure identifier,annotation confidence,adduct/retention context,processed same-specimen values and missingness',actual_external_coverage='NOT_VERIFIED'))
    save(features,'metabolite23_measurement_requirements.tsv')
    assert len(features)==23 and len(matrix)==351 and len(rows)==39
    assert [(x['metabolite_key'],x['gene'])for x in rows]==[(x['metabolite_key'],x['gene'])for x in original]
    assert all(x['p_value']==x['q_value']=='NA'for x in rows)
    (OUT/'analysis_spec.json').write_text(json.dumps(S,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    validation=dict(status='PASS',scope='format/coverage/provenance checks,not biological validation',relations=39,genes=35,metabolites=23,dataset_records=9,relation_dataset_checks=351,
        focused_evidence_records=7,new_to_project_records=5,new_publications_in2026=3,original_relation_order_unchanged=True,original_table_sha256=digest(R/S['source_table']),
        original_statistics_files_modified=False,new_tests=False,independent_ready_relationships=0,server_supplement_audit='ACCESS_BLOCKED_NOT_RUN',code_sha256=digest(Path(__file__)),spec_sha256=digest(SPEC))
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation))
if __name__=='__main__':main()
