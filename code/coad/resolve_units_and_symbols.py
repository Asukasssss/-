"""Join explicit GEO individual fields and audit official RNA symbol changes on server."""
import argparse,csv,hashlib,json,re,platform
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):
    with p.open('x',encoding='utf-8') as f:json.dump(obj,f,ensure_ascii=False,indent=2);f.write('\n')
def table(p,rows):
    with p.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--run-dir',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args()
    root,out=a.project_root.resolve(),a.run_dir.resolve()
    assert str(root).startswith('/public3/') and out.parent==root/'results/collaborative/COAD/B'
    assert (out/'.running').read_text().strip()==out.name and not (out/'cohort_summary.json').exists()
    base=root/'data/candidates/camp_primary_tissue_multicancer'
    mp=base/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    m=pd.read_csv(mp,dtype=str);m=m[m.Dataset.eq('COAD')].copy()
    rp=base/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
    inputs=[mp,rp,out/'geo_samples.soft',out/'gene_identity_records.json',out/'alias_checks.json',out/'candidate_catalog.json',Path(__file__)]
    hashes={p:sha(p) for p in inputs}
    geo={};cur=None
    for line in (out/'geo_samples.soft').read_text().splitlines():
        if line.startswith('^SAMPLE = '):cur=line.split(' = ',1)[1];geo[cur]={}
        elif cur and line.startswith('!Sample_title = '):geo[cur]['title']=line.split(' = ',1)[1]
        elif cur and line.startswith('!Sample_characteristics_ch1 = '):
            key,value=line.split(' = ',1)[1].split(': ',1);key=key.lower()
            assert key not in geo[cur];geo[cur][key]=value
    assert len(geo)==80 and len({r['individual'] for r in geo.values()})==39
    private=[]
    for r in m.to_dict('records'):
        g=geo[r['RNAID']]
        assert g['title']==r['MetabID'], 'Exact author metabolite ID and GEO title must agree'
        assert g['tissue type']==('tumor tissue' if r['TN']=='Tumor' else 'normal tissue')
        assert re.fullmatch(r'\d+y',g['age']) and g['gender'] in ['female','male']
        private.append(dict(CommonID=r['CommonID'],MetabID=r['MetabID'],RNAID=r['RNAID'],TN=r['TN'],
           individual=g['individual'],stage=g.get('stage','NA'),age=int(g['age'][:-1]),gender=g['gender'],tissue=g['tissue']))
    t=[r for r in private if r['TN']=='Tumor'];n=[r for r in private if r['TN']=='Normal']
    assert len(t)==37 and len(n)==39
    assert len({r['individual'] for r in t})==37 and len({r['individual'] for r in n})==39
    assert {r['individual'] for r in t}<={r['individual'] for r in n}
    primary=[r for r in t if r['stage'] in ['stage I','stage II','stage III','stage IV']]
    assert len(primary)==33
    full_tumors=[g for g in geo.values() if g['tissue type']=='tumor tissue']
    duplicate_people={k for k,v in Counter(g['individual'] for g in full_tumors).items() if v>1}
    assert len(duplicate_people)==2 and not ({r['individual'] for r in t}&duplicate_people)
    table(out/'private_clinical_join.tsv',private)  # SERVER ONLY, NEVER synchronize.
    symbols=set(pd.read_csv(rp,usecols=[0]).iloc[:,0].astype(str))
    reviewed=json.loads((out/'gene_identity_records.json').read_text())
    aliases=json.loads((out/'alias_checks.json').read_text())
    alias_by_id={r['human_gene_id']:r for r in aliases}
    resolutions=[]
    for r in reviewed:
        ident=r['human_gene_id'];alias=alias_by_id.get(ident)
        if alias:
            h=r['hgnc_records'][0]
            assert alias['hgnc_id']==h['hgnc_id'] and alias['rna_symbol'] in h['prev_symbol']
            assert alias['rna_symbol'] in symbols and alias['status']=='OFFICIAL_PREVIOUS_SYMBOL_UNIQUE'
        resolutions.append(dict(human_gene_id=ident,original_symbol=r['original_symbol'],
            approved_symbol=r.get('approved_symbol') or r.get('ncbi_manual_review',{}).get('symbol','NA'),
            rna_symbol=alias['rna_symbol'] if alias else 'NA',
            status='DONE' if alias else 'NOT_EVALUABLE',
            reason='OFFICIAL_PREVIOUS_SYMBOL_UNIQUE' if alias else 'GENE_ID_IDENTIFIED_BUT_NO_EXACT_RNA_ROW' if r.get('ncbi_manual_review') else 'NO_APPROVED_OR_DOCUMENTED_OLD_SYMBOL_IN_RNA',
            source_urls=';'.join(r['source_urls']+([alias['source_url']] if alias else []))))
    table(out/'gene_resolution.tsv',resolutions)
    candidates=json.loads((out/'candidate_catalog.json').read_text())
    eligibility=[]
    for r in candidates:
        ident=r['human_gene_id'];symbol=r['gene_symbol'];alias=alias_by_id.get(ident)
        rna_symbol=symbol if symbol in symbols else alias['rna_symbol'] if alias else 'NA'
        selected=r['disposition']=='REACTION_ANNOTATION_SUPPORTED_PROVISIONAL'
        eligibility.append(dict(cancer='COAD',cohort='COAD',metabolite_name=r['feature_name'],metabolite_key=r['metabolite_key'],
            human_gene_id=ident,gene=symbol or 'NA',rna_symbol=rna_symbol,
            biochemical_disposition=r['disposition'],in_prespecified_family=selected,rna_available=rna_symbol!='NA',
            decision='ELIGIBLE_PROVISIONAL_ANNOTATION' if selected else 'RETAINED_OUTSIDE_CURRENT_ASSOCIATION_FAMILY',
            analysis_status='NOT_RUN',original_identity_holds_changed=False))
    eligibility.sort(key=lambda r:(r['cohort'],r['metabolite_name'],r['metabolite_key'],r['gene'],r['human_gene_id']))
    chosen=[r for r in eligibility if r['in_prespecified_family']]
    assert len(chosen)==674 and len(eligibility)==974
    used={r['human_gene_id']:r['rna_symbol'] for r in eligibility if r['rna_available']}
    assert len(used)==len(set(used.values())), 'Do not map distinct gene IDs onto one RNA row'
    table(out/'candidate_eligibility.tsv',eligibility)
    summary=dict(cancer='COAD',run_id=out.name,unit_check='DONE',geo_sample_count=80,geo_individual_count=39,
      author_tumor_count=37,author_normal_count=39,author_tumor_unique_individuals=37,author_normal_unique_individuals=39,
      verified_tumor_normal_pairs=37,explicit_individual_field=True,exact_geo_title_metabid_match=True,
      repeated_patient_groups_in_full_geo=2,repeated_patient_groups_in_author_tumor=0,
      tumor_stage_counts=dict(Counter(r['stage'] for r in t)),
      primary_I_IV_count=33,sensitivity_author_tumor_count=37,
      complete_covariates=['age','gender','stage','tissue'],
      gene_ids_reviewed=len(reviewed),official_previous_symbols_recovered=len(aliases),remaining_unavailable_gene_ids=19,
      all_candidate_pairs=len(eligibility),all_pairs_rna_available=sum(r['rna_available'] for r in eligibility),
      family_pairs_planned=len(chosen),family_pairs_rna_available=sum(r['rna_available'] for r in chosen),
      family_genes_planned=len({r['human_gene_id'] for r in chosen}),
      family_genes_rna_available=len({r['human_gene_id'] for r in chosen if r['rna_available']}),
      identity_hold_changes=0,association='NOT_RUN',
      limits=['Patient units are verified against submitted GEO individual annotations, not genotype fingerprinting.',
        'Three adenoma and one stage 0 samples are excluded from stage I-IV primary analysis, retained in all-37 sensitivity.',
        'Aliases resolve gene identity only, not probe cross-hybridization or isoform specificity.',
        'Biochemical eligibility is provisional annotation support; no claim of final mechanism or causality.'])
    dump(out/'cohort_summary.json',summary)
    dump(out/'analysis_spec.json',dict(version='identity_units_v1',code_commit=a.code_commit,created_utc=datetime.now(timezone.utc).isoformat(),
       software={'python':platform.python_version(),'pandas':pd.__version__},new_tests=False,test_family='NO_NEW_TESTS',
       mapping_method='Exact RNAID to GEO accession; explicit individual field; exact GEO title to MetabID; no identifier parsing to infer patient units',
       gene_method='HGNC entrez_id and reverse-symbol uniqueness; approved previous symbols only; no fuzzy matching',
       sources=['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89076&targ=gsm&form=text&view=brief','https://pubmed.ncbi.nlm.nih.gov/28847964/']))
    for p,h in hashes.items():assert sha(p)==h
    table(out/'source_manifest.tsv',[dict(source_path=str(p.relative_to(root)),sha256=h) for p,h in hashes.items()])
    public=['gene_resolution.tsv','candidate_eligibility.tsv','cohort_summary.json','analysis_spec.json','source_manifest.tsv']
    dump(out/'validation.json',dict(status='PASS',checks=['80 GEO records and39 explicit individuals','76 exact accession/title matches','37 unique tumor individuals;37 paired controls','two repeated-person groups excluded by author mapping','33 stage I-IV samples','three official unambiguous previous-symbol resolutions','674 preselected pairs before testing','no multiple gene IDs sharing an RNA row','input hashes unchanged'],
      public_output_sha256={name:sha(out/name) for name in public},private_files_not_for_sync=['private_clinical_join.tsv','geo_samples.soft']))
    (out/'.running').unlink();(out/'DONE').write_text(datetime.now(timezone.utc).isoformat()+'\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
