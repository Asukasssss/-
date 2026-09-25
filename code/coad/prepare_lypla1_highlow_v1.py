"""Server-only full-gene pseudobulk extraction; no cell or patient export."""
import argparse,gzip,io,json,sys,tarfile,urllib.request,zipfile,platform
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PRIOR=ROOT/'results/collaborative/COAD/B/20260925T140257Z_uckl1_cna_v1'
sys.path.insert(0,str(PRIOR))
from run_uckl1_cna_v1 import build_labels,sha
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out;pub=out/'public'
assert out.parent==ROOT/'results/collaborative/COAD/B' and (out/'.running').is_dir()
spec=json.loads((out/'analysis_spec.json').read_text());assert spec['min_cells_per_group']==20
j,paths,_,audit=build_labels()
assert {k:sha(p) for k,p in paths.items()}==json.loads((PRIOR/'admission_hashes.json').read_text())
old=pd.read_csv(paths['cell_metadata'],sep='\t',dtype=str,keep_default_na=False)
ix=pd.Index(old.cell_id).get_indexer(j.cell_id);assert (ix>=0).all()
with np.load(paths['counts']) as z:count=z['LYPLA1'][ix];total=z['total'][ix]
j['count']=count;j['total']=total;j['value']=np.log1p(count/total*10000)
cna=j[j.group=='tumor_CNA'].copy();pos=cna[cna['count']>0].copy()
pos['expression_group']='';coverage=[]
for patient,f in pos.groupby('case_id'):
    med=f.value.median();lab=np.where(f.value>med,'high','low');nhi=int((lab=='high').sum());nlo=int((lab=='low').sum())
    valid=min(nhi,nlo)>=20
    if valid:pos.loc[f.index,'expression_group']=lab
    coverage.append(dict(patient=patient,positive_cells=len(f),high=nhi,low=nlo,threshold=med,included=valid))
pd.DataFrame(coverage).to_csv(out/'private_group_admission.tsv',sep='\t',index=False)
sel=pos[pos.expression_group!=''].copy();donors=sorted(sel.case_id.unique());assert len(donors)>=spec['min_donors']
sample_rows=[]
for i,(patient,grp) in enumerate((p,g) for p in donors for g in ['low','high']):
    f=sel[(sel.case_id==patient)&(sel.expression_group==grp)]
    sample_rows.append(dict(sample='PB'+str(i+1),patient=patient,group=grp,cells=len(f),mean_total_umi=f.total.mean(),mean_lypla1=f.value.mean()))
samples=pd.DataFrame(sample_rows);samples.to_csv(out/'private_samples.tsv',sep='\t',index=False)
key={(r.patient,r.group):r.sample for r in samples.itertuples()};sel['pb']=list(map(key.get,zip(sel.case_id,sel.expression_group)))
sel[['cell_id','case_id','expression_group','pb']].to_csv(out/'private_cell_groups.tsv',sep='\t',index=False)
source=ROOT/'data/candidates/coad_uhlitz_20260921/counts.tar'
oldval=ROOT/'results/collaborative/COAD/B/20260922T141755Z_source_contract_sc_v2/public/extraction_validation.json'
assert sha(source)==json.loads(oldval.read_text())['source_hashes']['counts.tar']
genes=None;bulk=None;seen=set();summed=np.zeros(len(samples),dtype=np.int64)
with tarfile.open(source) as tar:
    for member in tar.getmembers():
        if not member.isfile():continue
        stream=io.TextIOWrapper(gzip.GzipFile(fileobj=tar.extractfile(member)))
        header=stream.readline().rstrip('\r\n').split('\t');assert header[0]=='gene'
        ids=header[1:];idx=pd.Index(sel.cell_id).get_indexer(ids);use=np.flatnonzero(idx>=0)
        if not len(use):stream.close();continue
        matched=sel.iloc[idx[use]];assert not seen.intersection(matched.cell_id);seen.update(matched.cell_id)
        groups=[np.flatnonzero(matched.pb.to_numpy()==s) for s in samples['sample']]
        names=[];part=[];den=np.zeros(len(use),dtype=np.int64);target=None
        for line in stream:
            name,vals=line.rstrip('\r\n').split('\t',1);v=np.fromstring(vals,sep='\t',dtype=np.int64)[use]
            assert (v>=0).all();names.append(name);den+=v
            if name=='LYPLA1':target=v.copy()
            part.append([v[z].sum() for z in groups])
        stream.close();assert np.array_equal(den,matched.total.to_numpy()) and np.array_equal(target,matched['count'].to_numpy())
        part=np.array(part,dtype=np.int64)
        if genes is None:genes=names;bulk=part
        else:assert names==genes;bulk+=part
        print(json.dumps(dict(processed_matrix=member.name,cells=len(use),genes=len(names))),flush=True)
assert seen==set(sel.cell_id) and len(set(genes))==len(genes)
pd.DataFrame(bulk,index=genes,columns=samples['sample']).to_csv(out/'private_pseudobulk_counts.tsv',sep='\t',index_label='gene')
url='https://reactome.org/download/current/ReactomePathways.gmt.zip'
urllib.request.urlretrieve(url,out/'reactome_source.zip')
with zipfile.ZipFile(out/'reactome_source.zip') as z:
    members=[n for n in z.namelist() if n.endswith('.gmt')];assert len(members)==1
    (out/'reactome.gmt').write_bytes(z.read(members[0]))
summary=dict(status='PASS',cna_cells=len(cna),cna_donors=cna.case_id.nunique(),lypla1_positive_cells=len(pos),zero_cells_excluded=int((cna['count']==0).sum()),
    included_donors=len(donors),excluded_donors=sum(not x['included'] for x in coverage),included_positive_cells=len(sel),full_genes=len(genes),
    high_cells=int((sel.expression_group=='high').sum()),low_cells=int((sel.expression_group=='low').sum()),
    min_group_cells=int(samples.cells.min()),max_group_cells=int(samples.cells.max()),
    exact_UMI_denominator_and_LYPLA1_count_rechecked=True,source_hash_matches_history=True)
(pub/'admission.json').write_text(json.dumps(summary,indent=2)+'\n')
diags=[]
for grp,f in sel.groupby('expression_group'):
    diags.append(dict(group=grp,cells=len(f),donors=f.case_id.nunique(),mean_log1p_CP10K=f.value.mean(),median_total_UMI=f.total.median(),mean_total_UMI=f.total.mean(),median_LYPLA1_count=f['count'].median()))
pd.DataFrame(diags).to_csv(pub/'group_summary.tsv',sep='\t',index=False)
paths.update(full_counts=source,reactome_zip=out/'reactome_source.zip',reactome_gmt=out/'reactome.gmt')
pd.DataFrame([dict(source_id=k,source_path=str(p),sha256=sha(p),url=url if k.startswith('reactome') else 'historical source') for k,p in paths.items()]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
print(json.dumps(summary),flush=True)
