"""Read public resource metadata and a gene dictionary; never download matrices.

Run with --out pointing to a new directory. Live access results can change.
Dataset biological eligibility is a separate, versioned manual decision.
"""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

GEO = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=self&form=text&view=full'
FTP = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
COLLECTION = '48259aa8-f168-4bf5-b797-af8e88da6637'
ENDPOINTS = [
    ('Wu_series', 'GET', GEO.format('GSE176078')),
    ('Pal_series', 'GET', GEO.format('GSE161529')),
    ('PDX_reuse_series', 'GET', GEO.format('GSE276609')),
    ('Wu_processed_archive', 'HEAD', FTP+'GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'),
    ('Pal_count_archive', 'HEAD', FTP+'GSE161nnn/GSE161529/suppl/GSE161529_RAW.tar'),
    ('Pal_feature_dictionary', 'GET', FTP+'GSE161nnn/GSE161529/suppl/GSE161529_features.tsv.gz'),
    ('Pal_annotated_objects', 'GET', 'https://api.figshare.com/v2/articles/17058077'),
    ('Reed_collection', 'GET', 'https://api.cellxgene.cziscience.com/curation/v1/collections/'+COLLECTION),
    ('UPP1_mouse_code', 'GET', 'https://raw.githubusercontent.com/chris-mcginnis-ucsf/pymt_atlas/master/README.md'),
    ('metastatic_spatial_portal', 'GET', 'https://singlecell.broadinstitute.org/single_cell/study/SCP3851/a-spatial-single-cell-transcriptomic-atlas-of-metastatic-breast-cancer-progression'),
]


def write_tsv(path, rows, fields=None):
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def fetch(endpoint):
    name, method, url = endpoint
    item = dict(cancer='BRCA', resource_id=name, method=method, url=url, http_status='NA', content_type='NA', declared_bytes='NA', retrieved_metadata_bytes=0, metadata_sha256='NA', checked_utc=datetime.now(timezone.utc).isoformat(), status='ACCESS_BLOCKED', error='')
    try:
        response = requests.request(method, url, timeout=(10, 25), allow_redirects=True, stream=True)
        item.update(http_status=response.status_code, content_type=response.headers.get('Content-Type', 'NA'), declared_bytes=response.headers.get('Content-Length', 'NA'))
        body = b''
        if method == 'GET':
            for chunk in response.iter_content(65536):
                body += chunk
                if len(body) > 3_000_000:
                    response.close()
                    raise ValueError('Metadata exceeds 3 MB limit; stopped')
            item.update(retrieved_metadata_bytes=len(body), metadata_sha256=hashlib.sha256(body).hexdigest())
        response.close()
        if 200 <= item['http_status'] < 300:
            item['status'] = 'DONE'
        return item, body
    except Exception as exc:
        item['error'] = str(exc)
        return item, b''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    pool_path = args.repo_root / 'reference/brca/candidate_comparison_v1.tsv'
    with pool_path.open(encoding='utf-8-sig') as handle:
        genes = [r['gene'] for r in csv.DictReader(handle, delimiter='\t')]
    assert len(genes) == len(set(genes)) == 117
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        responses = list(executor.map(fetch, ENDPOINTS))
    log = [r[0] for r in responses]
    bodies = {r[0]['resource_id']: r[1] for r in responses if r[0]['status'] == 'DONE'}
    write_tsv(args.out/'access_log.tsv', log)
    meta = []
    for name in ['Wu_series', 'Pal_series', 'PDX_reuse_series']:
        text = bodies.get(name, b'').decode('utf8', errors='replace')
        valid = text.startswith('^SERIES = ')
        fields = {}
        for line in text.splitlines():
            if ' = ' in line:
                key, value = line.split(' = ', 1)
                fields.setdefault(key, []).append(value)
        meta.append(dict(resource_id=name, status='DONE' if valid else 'NEEDS_REVIEW', title=';'.join(fields.get('!Series_title', [])), geo_accession=';'.join(fields.get('!Series_geo_accession', [])), listed_GSM_records=len(fields.get('!Series_sample_id', [])) if valid else 'NA', unique_patients='NOT_COUNTED', source_pmid=';'.join(fields.get('!Series_pubmed_id', [])), supplementary_urls=';'.join(fields.get('!Series_supplementary_file', [])), note='GSM count is not donor count;full sample maps not exported'))
    write_tsv(args.out/'series_metadata_summary.tsv', meta)
    coverage = []
    raw = bodies.get('Pal_feature_dictionary')
    features = list(csv.reader(io.StringIO(gzip.decompress(raw).decode('utf8')), delimiter='\t')) if raw else []
    symbols = {row[1] for row in features if len(row) >= 2}
    for gene in genes:
        coverage.append(dict(cancer='BRCA', dataset='GSE161529', gene=gene, in_feature_dictionary=str(gene in symbols) if raw else 'NA', expression_detection='NOT_RUN', cells_expressing='NA', patients_evaluable='NA', interpretation='Dictionary inclusion is not observed expression'))
    write_tsv(args.out/'gene_dictionary_coverage_117.tsv', coverage)
    assets = []
    if 'Reed_collection' in bodies:
        collection = json.loads(bodies['Reed_collection'])
        for dataset in collection.get('datasets', []):
            for asset in dataset.get('assets', []):
                assets.append(dict(cancer='BRCA', collection_id=COLLECTION, dataset_id=dataset['dataset_id'], title=dataset['title'], portal_cell_count=dataset['cell_count'], format=asset['filetype'], declared_bytes=asset['filesize'], url=asset['url'], downloaded='NO', note='Portal version count;not patients;subsets/integrated atlases can overlap'))
    if assets:
        write_tsv(args.out/'normal_reference_assets.tsv', assets)
    validation = dict(scope='Resource discovery and metadata/feature-dictionary checks only', gene_pool_size=117, pal_feature_rows=len(features) if raw else 'NA', pal_genes_in_dictionary=sum(g in symbols for g in genes) if raw else 'NA', pal_genes_missing=sorted(set(genes)-symbols) if raw else 'NOT_EVALUABLE', http_successes=sum(r['status']=='DONE' for r in log), endpoints=len(log), source_matrices_downloaded=0, expression_tests=0, annotation_joins_completed=0, normal_portal_assets=len(assets), gene_pool_sha256=hashlib.sha256(pool_path.read_bytes()).hexdigest(), limitations=['No matrix load or cell annotation join', 'No actual donor counts', 'No claim that all genes are expressed', 'No independent-cohort claim for reused datasets', 'HEAD availability does not verify full download integrity'])
    (args.out/'validation.json').write_text(json.dumps(validation, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
