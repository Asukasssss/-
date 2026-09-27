import requests,tarfile,hashlib,json
from pathlib import Path
root=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_Wu2021_Visium_Zenodo4739739')
expected={'filtered_count_matrices.tar.gz':'ea1220f0606c4d4e9307468fb2d5427a','metadata.tar.gz':'1d0b34ead70635c094ac74ac58a88d68','spatial.tar.gz':'29b02ae433d459d19dc52b75a529e8e9'}
for name,md5 in expected.items():
 p=root/name
 if not p.exists():
  with requests.get('https://zenodo.org/records/4739739/files/'+name+'?download=1',stream=True,timeout=90) as r:
   r.raise_for_status()
   with p.with_suffix('.part').open('wb') as f:
    for b in r.iter_content(1048576): f.write(b)
  p.with_suffix('.part').rename(p)
 assert hashlib.md5(p.read_bytes()).hexdigest()==md5,name
 with tarfile.open(p) as t:
  for m in t.getmembers():
   assert (root/m.name).resolve().is_relative_to(root.resolve()) and not m.issym() and not m.islnk()
  t.extractall(root)
 print(name,'MD5_PASS',flush=True)
print('SPATIAL_FILES',[str(x.relative_to(root)) for x in (root/'spatial').rglob('*') if x.is_file()][:25])
