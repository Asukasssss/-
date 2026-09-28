"""Read-only public metadata probes; store responses on server165 only.

Does not download spectra, individual matrices, or infer paired samples.
"""
import concurrent.futures
import hashlib
import json
import pathlib
import sys
import urllib.request

ROOT = pathlib.Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
run = pathlib.Path(sys.argv[1]).resolve()
assert run.parent == ROOT / 'results/collaborative/COAD/B'
assert (run / '.running').is_dir()
dest = ROOT / 'data/candidates/coad_paired_search_20260921'
dest.mkdir(exist_ok=True)
sources = {
    'integral_article': 'https://pubs.acs.org/doi/10.1021/acs.analchem.4c04421',
    'integral_bioproject': 'https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA026610',
    'integral_gsa': 'https://ngdc.cncb.ac.cn/gsa-human/browse/HRA007600',
    'integral_preprint': 'https://www.biorxiv.org/content/10.1101/2024.09.27.614689v1.supplementary-material',
    'prrx2_preprint': 'https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7030030',
    'massive_metadata': 'https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=9a52b4480ed34cb2b5ac9da0333105f6',
}

def probe(item):
    key, url = item
    record = dict(source_id=key, url=url, processed_matrix_acquired=False)
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            blob = response.read(5_000_001)
            if len(blob) > 5_000_000:
                raise ValueError('metadata response exceeded size limit')
            target = dest / (key + '.html')
            target.write_bytes(blob)
            record.update(http_status=response.status, bytes=len(blob),
                          sha256=hashlib.sha256(blob).hexdigest(), server_path=str(target),
                          status='RETRIEVED_NOT_CONTENT_VALIDATED')
    except Exception as exc:
        record.update(status='ACCESS_BLOCKED', error=str(exc), sha256='NA')
    return record

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    records = list(pool.map(probe, sources.items()))
output = dict(scope='Public metadata retrieval, not dataset eligibility validation', records=records,
              new_patient_tests=0, raw_spectra_downloaded=False)
(run / 'retrieval_log.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps(output))
