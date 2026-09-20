"""Reuse existing117 literature; bounded accession discovery, no fabricated data availability."""
import argparse,json,re,concurrent.futures,xml.etree.ElementTree as ET,hashlib
from pathlib import Path
import requests,pandas as pd

def main(root):
    src=root/'source';pub=root/'public';studies=json.loads((src/'curation_records.json').read_text());genes=json.loads((src/'genes117.json').read_text());assert len(genes)==117
    cache=src/'literature_metadata';cache.mkdir(exist_ok=True)
    def get(url,path):
        if path.exists():return path.read_text()
        x=requests.get(url,timeout=25);x.raise_for_status();path.write_text(x.text);return x.text
    def one(s):
        url=s['source_url'];sid=s['study_id'];pmid=re.search(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)',url);pmc=re.search(r'PMC\d+',url);doi=re.search(r'10\.\d{4,9}/[^\s?#]+',url)
        if not doi and 'nature.com/articles/' in url:doi='10.1038/'+url.split('/articles/')[1].split('?')[0]
        elif doi:doi=doi.group(0)
        query=('EXT_ID:'+pmid.group(1)+' AND SRC:MED') if pmid else ('DOI:'+doi if doi else None)
        row=dict(study_id=sid,genes=';'.join(s['genes']),original_source=url,original_model=s['model'],original_intervention=s['intervention'],original_endpoint=s['outcome'],original_evidence_class=s['evidence_class'],retrieval_status='NOT_EVALUABLE',reason='',pmcid='',accessions='',dataset_design_verified=False)
        try:
            if pmc:pid=pmc.group(0)
            elif query:
                from urllib.parse import quote
                u='https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=5&resultType=core&query='+quote(query)
                hits=json.loads(get(u,cache/(sid+'_search.json')))['resultList']['result']
                hit=[h for h in hits if (str(h.get('id'))==pmid.group(1) if pmid else h.get('doi','').lower()==str(doi).lower())]
                if len(hit)!=1:row['reason']='no_unique_exact_PMID_or_DOI_record';return row
                pid=hit[0].get('pmcid')
                if not pid:row['reason']='no_PMC_fulltext_in_exact_record';return row
            else:row['reason']='source_has_no_resolvable_PMID_PMC_or_DOI';return row
            row['pmcid']=pid;u='https://www.ebi.ac.uk/europepmc/webservices/rest/'+pid+'/fullTextXML'
            xml=get(u,cache/(sid+'_fulltext.xml'));tree=ET.fromstring(xml)
            for parent in tree.iter():
                for child in list(parent):
                    if child.tag in ['ref-list','references']:parent.remove(child)
            text=' '.join(tree.itertext())
            # Accession presence is merely a lead; not assigned to the targeted intervention here.
            acc=sorted(set(re.findall(r'\b(?:GSE\d+|PRJNA\d+|SRP\d{5,}|E-MTAB-\d+|PXD\d+)\b',text)))
            row.update(retrieval_status='DONE',reason='fulltext_read_machine_accession_discovery;intervention_and_replication_design_not_yet_confirmed',accessions=';'.join(acc))
        except Exception as e:row.update(retrieval_status='ACCESS_BLOCKED',reason=type(e).__name__+':'+str(e)[:140])
        return row
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for i,r in enumerate(ex.map(one,studies)):
            rows.append(r)
            if (i+1)%20==0:print('literature',i+1,flush=True);pd.DataFrame(rows).to_csv(pub/'study_resource_discovery.tsv',sep='\t',index=False)
    table=pd.DataFrame(rows);table.to_csv(pub/'study_resource_discovery.tsv',sep='\t',index=False)
    accessions=sorted({a for r in rows for a in r['accessions'].split(';') if a.startswith('GSE')})
    def geo(acc):
        row=dict(accession=acc,status='NOT_EVALUABLE',title='',overall_design='',n_samples='',reason='')
        try:
            u=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full';t=get(u,cache/(acc+'.txt'))
            assert '^SERIES = '+acc in t
            title=re.findall(r'!Series_title = (.*)',t);design=re.findall(r'!Series_overall_design = (.*)',t)
            row.update(status='DONE',title=re.sub(r'\s+',' ',' '.join(title)).strip(),overall_design=re.sub(r'\s+',' ',' '.join(design)).strip(),n_samples=len(re.findall(r'!Series_sample_id = ',t)),reason='repository_metadata_only;do_not_infer_target_gene_perturbation_from_article_link')
        except Exception as e:row.update(status='ACCESS_BLOCKED',reason=str(e)[:150])
        return row
    gr=list(concurrent.futures.ThreadPoolExecutor(max_workers=3).map(geo,accessions));pd.DataFrame(gr).to_csv(pub/'GEO_resource_metadata.tsv',sep='\t',index=False)
    gene_rows=[]
    for g in genes:
        ss=[r for r in rows if g in r['genes'].split(';')];full=[r for r in ss if r['retrieval_status']=='DONE'];acc=sorted({a for r in ss for a in r['accessions'].split(';') if a})
        gene_rows.append(dict(gene=g,existing_study_records=len(ss),fulltexts_retrieved=len(full),study_ids=';'.join(r['study_id'] for r in ss),accession_leads=';'.join(acc),status='NEEDS_REVIEW' if acc else ('NOT_EVALUABLE' if not ss or len(full)==len(ss) else 'PARTIAL'),reason='accession_leads_require_gene_model_endpoint_replicate_verification' if acc else ('no_existing_study_record;not_a_negative_function_result' if not ss else 'no_accession_detected_in_accessible_existing_fulltexts;not_proof_of_absence'),ready_for_target_specific_analysis=False))
    pd.DataFrame(gene_rows).to_csv(pub/'all117_intervention_resource_inventory.tsv',sep='\t',index=False)
    manifest=[dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted(cache.iterdir()) if f.is_file()];pd.DataFrame(manifest).to_csv(pub/'resource_source_manifest.tsv',sep='\t',index=False)
    summary={'status':'PARTIAL','scope':'all117 existing-literature resource discovery, not full manual design validation','genes':len(gene_rows),'studies':len(rows),'retrieval_counts':table.retrieval_status.value_counts().to_dict(),'genes_with_accession_leads':sum(bool(r['accession_leads']) for r in gene_rows),'GEO_accessions':len(gr),'gene_targeted_perturbations_verified_in_this_pass':0}
    (pub/'resource_validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
