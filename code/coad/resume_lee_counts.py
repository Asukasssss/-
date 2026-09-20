"""Resume a known partial GEO file using four verified HTTP range requests."""
import argparse
import concurrent.futures
import gzip
import json
import time
import urllib.request
from pathlib import Path
from cell_origin_v1 import ROOT, sha

URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz'
TOTAL = 128941455

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, required=True)
    a = p.parse_args()
    d = a.data_dir.resolve()
    assert d.parent == ROOT / 'data/candidates'
    src = d / 'lee_raw_UMI.txt.gz'
    start = src.stat().st_size
    assert 0 < start < TOTAL
    parts = d / 'lee_resume_parts'
    parts.mkdir(exist_ok=False)
    bounds = [start + (TOTAL - start) * i // 4 for i in range(5)]
    def fetch(i):
        left, right = bounds[i], bounds[i + 1] - 1
        path = parts / str(i)
        for attempt in range(5):
            offset = path.stat().st_size if path.exists() else 0
            if offset == right - left + 1:
                return path
            try:
                req = urllib.request.Request(URL, headers={'Range': f'bytes={left + offset}-{right}'})
                with urllib.request.urlopen(req, timeout=90) as response:
                    assert response.status == 206
                    assert response.headers['Content-Range'] == f'bytes {left + offset}-{right}/{TOTAL}'
                    with path.open('ab') as f:
                        while True:
                            b = response.read(65536)
                            if not b:
                                break
                            f.write(b)
                assert path.stat().st_size == right - left + 1
                print(json.dumps({'part': i, 'status': 'DONE', 'bytes': path.stat().st_size}), flush=True)
                return path
            except Exception as e:
                print(json.dumps({'part': i, 'attempt': attempt, 'error': str(e)}), flush=True)
                time.sleep(2)
        raise RuntimeError(f'Part {i} incomplete')
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        paths = list(pool.map(fetch, range(4)))
    assert src.stat().st_size == start, 'Original partial file must not be concurrently written'
    dest = d / 'lee_raw_UMI.complete.txt.gz'
    with dest.open('xb') as out:
        for path in [src] + paths:
            with path.open('rb') as f:
                for b in iter(lambda: f.read(1024 * 1024), b''):
                    out.write(b)
    assert dest.stat().st_size == TOTAL
    with gzip.open(dest, 'rb') as f:
        for _ in iter(lambda: f.read(1024 * 1024), b''):
            pass
    (parts / 'validation.json').write_text(json.dumps({'url': URL, 'bytes': TOTAL, 'sha256': sha(dest), 'gzip_crc': 'PASS', 'original_prefix_bytes': start, 'ranges': bounds}, indent=2))
    print('DOWNLOAD_VALIDATED', flush=True)

if __name__ == '__main__':
    main()
