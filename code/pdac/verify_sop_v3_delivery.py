"""Validate repository-bound public delivery and exact staged checksum bytes."""
import hashlib,subprocess,json
from record_sop_v3_resume import SOURCE,INTEGRATION,REPO,read,TEMPLATES

def main():
    n=0
    for folder in [SOURCE,INTEGRATION]:
        for _,r in read(folder/'checksums.tsv').iterrows():
            p=(folder/r.file).relative_to(REPO).as_posix()
            raw=subprocess.check_output(['git','show',':'+p],cwd=REPO)
            assert hashlib.sha256(raw).hexdigest()==r.sha256,p
            n+=1
    catalog=read(INTEGRATION/'SOP_DELIVERY_STATUS.tsv');assert len(catalog)==31
    for _,r in catalog.iterrows():
        if r.data_scope=='SERVER_PRIVATE':continue
        frames=[read(REPO/p) for p in r.result_paths.split(';')]
        for d in frames:assert set(read(TEMPLATES/r.file).columns)<=set(d.columns),r.file
    profiles=read(SOURCE/'sc_celltype_profiles.tsv');assert len(profiles)==17862
    parts=read(SOURCE/'sc_celltype_profile_parts.tsv');lines=[]
    for i,r in parts.iterrows():
        data=(SOURCE/r.file).read_bytes();assert hashlib.sha256(data).hexdigest()==r.sha256
        lines.extend(data.splitlines(keepends=True)[0 if i==0 else 1:])
    assert hashlib.sha256(b''.join(lines)).hexdigest()==parts.logical_table_sha256.iloc[0]
    subprocess.run(['git','diff','--cached','--check'],cwd=REPO,check=True,capture_output=True)
    print(json.dumps({'status':'PASS','staged_hashes':n,'SOP_tables':31,'source_partitions_lossless':True}))
if __name__=='__main__':main()
