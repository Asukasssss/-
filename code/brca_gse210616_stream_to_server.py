"""Stream public data through RAM to server165. No local source-data files."""
import concurrent.futures,subprocess,requests,re,time
BASE='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_GSE210616'
def remote(cmd):return subprocess.check_output(['ssh','server165',cmd],text=True)
def transfer(chunks,path):
    assert re.fullmatch(r'[A-Za-z0-9_./-]+',path)
    proc=subprocess.Popen(['ssh','server165',f"cat > '{path}.relaypart' && mv '{path}.relaypart' '{path}'"],stdin=subprocess.PIPE)
    try:
        for b in chunks:proc.stdin.write(b)
        proc.stdin.close();assert proc.wait()==0
    except Exception:
        proc.kill();raise
def one(gsm):
    for attempt in range(4):
        try:
            r=requests.get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=quick',timeout=60);r.raise_for_status()
            assert '^SAMPLE' in r.text
            d=f'{BASE}/{gsm}';remote(f"mkdir -p '{d}'")
            transfer([r.content],d+'/metadata.soft')
            urls=re.findall(r'!Sample_supplementary_file(?:_\d+)? = (\S+)',r.text)
            for url in urls:
                if '.cloupe' in url:continue
                url=url.replace('ftp://','https://');p=d+'/'+url.rsplit('/',1)[1]
                with requests.get(url,stream=True,timeout=(30,90)) as data:
                    data.raise_for_status();transfer(data.iter_content(1024*1024),p)
            print(gsm,'STREAMED_TO_SERVER',flush=True);return
        except Exception as e:
            print(gsm,'RETRY',attempt,type(e).__name__,flush=True)
            if attempt==3:raise
            time.sleep(3)
if __name__=='__main__':
    s=remote(f"cat '{BASE}/series_metadata.soft'")
    gsms=re.findall(r'!Series_sample_id = (GSM\d+)',s);assert len(gsms)==43
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(one,gsms))
