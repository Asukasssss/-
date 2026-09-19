"""Download pinned study assets to a server run; verify against source_manifest.tsv."""
import argparse
import csv
import hashlib
from pathlib import Path
import requests


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
    return h.hexdigest()


def main(root, manifest):
    (root/'source').mkdir(parents=True,exist_ok=True)
    with manifest.open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    for r in rows:
        if r['kind']!='expression_asset':continue
        dest=root/'source'/r['file_name']
        if dest.exists():
            assert digest(dest)==r['sha256'],str(dest)
            continue
        part=dest.with_suffix(dest.suffix+'.part')
        with requests.get(r['url'],stream=True,timeout=(20,90)) as response:
            response.raise_for_status()
            with part.open('xb') as out:
                for block in response.iter_content(4*1024*1024):out.write(block)
        assert part.stat().st_size==int(r['bytes'])
        assert digest(part)==r['sha256']
        part.rename(dest)
        print('verified',r['file_name'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--manifest',type=Path,required=True)
    a=p.parse_args();main(a.root,a.manifest)
