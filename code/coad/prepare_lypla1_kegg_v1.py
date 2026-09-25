"""Freeze official KEGG annotations and prepare enrichment without rerunning DE."""
import argparse, collections, datetime, hashlib, json, re, subprocess, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PRIOR=ROOT/'results/collaborative/COAD/B/20260925T153907Z_lypla1_fc025_v3/public'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--commit',required=True);ap.add_argument('--cached-sources',action='store_true');a=ap.parse_args()
out=a.out;pub=out/'public'
assert out.parent==ROOT/'results/collaborative/COAD/B' and (out/'.running').is_dir() and len(a.commit)==40
pub.mkdir(exist_ok=a.cached_sources);src=out/'sources';src.mkdir(exist_ok=a.cached_sources)
if a.cached_sources:assert not (pub/'kegg_ORA.tsv').exists(), 'Do not overwrite computed results'
manifest=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for key,url in [('pathways','https://rest.kegg.jp/list/pathway/hsa'),('genes','https://rest.kegg.jp/list/hsa'),('links','https://rest.kegg.jp/link/pathway/hsa'),('hierarchy','https://rest.kegg.jp/get/br:br08901/json'),('release','https://rest.kegg.jp/info/kegg')]:
    p=src/(key+'.txt')
    if not a.cached_sources:
        subprocess.run(['curl','--fail','--location','--retry','2','--max-time','120','-A','COAD-academic-analysis/1.0',url,'-o',str(p)],check=True)
    assert p.stat().st_size>100
    manifest.append(dict(source_id=key,source_path=str(p),url=url,sha256=sha(p),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),transport='local memory to server SFTP' if a.cached_sources else 'server HTTPS'))
    time.sleep(1)
d=pd.read_csv(PRIOR/'results.tsv',sep='\t')
assert len(d)==11058 and d.gene.is_unique and 'LYPLA1' not in set(d.gene)
selected=(d.p_value<.05)&(d.log2FC.abs()>=.25)
assert selected.sum()==434 and ((d.log2FC>0)&selected).sum()==158
titles={}
for line in (src/'pathways.txt').read_text().splitlines():
    k,title=line.split('\t',1);k=k.split(':')[-1];assert re.fullmatch('hsa[0-9]{5}',k)
    titles[k]=re.sub(r' - Homo sapiens \(human\)$','',title)
primary=collections.defaultdict(set);alias=collections.defaultdict(set)
for line in (src/'genes.txt').read_text().splitlines():
    z=line.split('\t');assert len(z)==4
    names=z[3].split(';',1)[0].split(', ')
    for n in names:alias[n].add(z[0])
    primary[names[0]].add(z[0])
mapping=[]
for gene in d.gene:
    exact=primary.get(gene,set());ids=exact or alias.get(gene,set())
    status=('PRIMARY_SYMBOL' if exact else 'UNIQUE_ALIAS') if len(ids)==1 else ('AMBIGUOUS' if ids else 'UNMAPPED')
    mapping.append(dict(gene=gene,kegg_gene=next(iter(ids)) if len(ids)==1 else '',candidate_ids=';'.join(sorted(ids)),status=status))
m=pd.DataFrame(mapping)
# A KEGG gene must not receive two tested symbols and be counted twice.
dup=m.kegg_gene.ne('')&m.kegg_gene.duplicated(keep=False)
m.loc[dup,'status']='MULTIPLE_TESTED_SYMBOLS_ONE_ID';m.loc[dup,'kegg_gene']=''
id_to_symbol=dict(zip(m.loc[m.kegg_gene.ne(''),'kegg_gene'],m.loc[m.kegg_gene.ne(''),'gene']))
pathgenes=collections.defaultdict(set)
for line in (src/'links.txt').read_text().splitlines():
    gid,pid=line.split('\t');pid=pid.split(':')[-1]
    if pid in titles:pathgenes[pid].add(gid)
categories=collections.defaultdict(list)
def walk(node,trail=()):
    name=node['name'];hit=re.match(r'^(\d{5})\s+',name)
    if hit and 'hsa'+hit[1] in titles:categories['hsa'+hit[1]].append(trail)
    for child in node.get('children',[]):walk(child,trail+(name,))
walk(json.loads((src/'hierarchy.txt').read_text()))
assert any(any(any('Lipid metabolism' in label for label in trail) for trail in paths) for paths in categories.values())
sets={k:{id_to_symbol[g] for g in pathgenes[k] if g in id_to_symbol} for k in titles}
annotated=set().union(*sets.values())
tested=set(d.gene);catrows=[]
for k in sorted(titles):
    K=len(sets[k]);trails=categories.get(k,[])
    broad= k.startswith('hsa011') or k.startswith('hsa012')
    eligible=15<=K<=500 and not broad
    catrows.append(dict(pathway_id=k,pathway=titles[k],source_gene_count=len(pathgenes[k]),tested_gene_count=K,
        status='DONE' if eligible else 'NOT_EVALUABLE',reason='ELIGIBLE' if eligible else ('GLOBAL_OVERVIEW_MAP' if broad else 'OUTSIDE_15_500_TESTED_GENES'),
        lipid_metabolism=any(any('Lipid metabolism' in t for t in trail) for trail in trails),classification=' || '.join(' > '.join(t) for t in trails)))
catalog=pd.DataFrame(catrows);eligible=set(catalog.loc[catalog.status=='DONE','pathway_id'])
assert len(eligible)>100
pd.DataFrame(mapping).assign(kegg_gene=m.kegg_gene,status=m.status,in_any_pathway=m.gene.isin(annotated)).to_csv(pub/'gene_mapping.tsv',sep='\t',index=False)
catalog.to_csv(pub/'pathway_coverage.tsv',sep='\t',index=False)
with (out/'kegg_tested.gmt').open('w') as f:
    for k in sorted(eligible):f.write('\t'.join([k,titles[k]]+sorted(sets[k]))+'\n')
pd.DataFrame([dict(pathway_id=k,gene=g) for k in sorted(eligible) for g in sorted(sets[k])]).to_csv(pub/'tested_pathway_members.tsv',sep='\t',index=False)
# Primary follows the previous all-tested-gene background; annotation-restricted
# sensitivity is reported alongside it, never selected by the smaller P.
rows=[]
for bgname,bg in [('all_tested',tested),('kegg_annotated_sensitivity',annotated)]:
    for direction,mask in [('higher_in_LYPLA1_high',d.log2FC>0),('lower_in_LYPLA1_high',d.log2FC<0)]:
        chosen=set(d.loc[selected&mask,'gene'])&bg;M=len(bg);N=len(chosen)
        for k in sorted(eligible):
            pg=sets[k]&bg;hits=sorted(pg&chosen);K=len(pg);n=len(hits)
            rows.append(dict(pathway_id=k,pathway=titles[k],background=bgname,direction=direction,background_genes=M,selected_genes=N,pathway_tested_genes=K,overlap=n,
                fold_enrichment=n*M/(N*K) if N*K else np.nan,p_value=float(hypergeom.sf(n-1,M,K,N)) if N else np.nan,overlap_genes=';'.join(hits),status='DONE' if N else 'NOT_EVALUABLE',reason='NOMINAL_P_ONLY'))
ora=pd.DataFrame(rows).sort_values(['background','p_value','pathway_id','direction']);ora['selected']=ora.p_value<.05
ora.to_csv(pub/'kegg_ORA.tsv',sep='\t',index=False)
cols=['gene','log2FC','QL_F','p_value','depth_log2FC','depth_p_value','donors_high_greater','donors_high_lower']
d[cols].to_csv(pub/'gene_statistics.tsv',sep='\t',index=False)
manifest.append(dict(source_id='prior_DE',source_path=str(PRIOR/'results.tsv'),url='prior frozen batch',sha256=sha(PRIOR/'results.tsv'),retrieved_utc='REUSED'))
for f in ['prepare_lypla1_kegg_v1.py','run_lypla1_kegg_gsea_v1.R']:
    manifest.append(dict(source_id=f,source_path=str(out/f),url='code lock '+a.commit,sha256=sha(out/f),retrieved_utc='CODE_LOCK'))
pd.DataFrame(manifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
spec=dict(analysis_version='COAD_LYPLA1_KEGG_v1',run_id=out.name,code_lock_commit=a.commit,prior_run=PRIOR.parent.name,
    cohort='Uhlitz_GSE166555',donors=8,cells=2528,genes_tested=len(d),selected_high=158,selected_low=276,
    DE_selection='raw P<0.05 and abs(log2FC)>=0.25',mapping='Exact primary KEGG symbol first, otherwise unique exact alias; ambiguous or duplicate target excluded; no fuzzy matching',
    pathway_scope='All human KEGG pathways; 15-500 measured mapped genes; global/overview hsa011/hsa012 excluded; disease pathways retained and interpreted as gene-set labels',
    ORA='Hypergeometric upper tail; high and low lists separate; all 11058 tested genes primary background; KEGG-annotated background sensitivity',
    GSEA='signed sqrt(QL_F); all 11058 tested genes; tiny negative null QL_F set to zero for ranking only; gene alphabetical order breaks equal-rank input order; no jitter',
    GSEA_parameters=dict(seed=20260925,minSize=15,maxSize=500,eps=0,sampleSize=101,nPermSimple=10000,nproc=1),
    raw_P_threshold=.05,FDR_selection=False,source_data_residency='server165 only',no_new_DE=True)
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
summary=dict(status='PREPARED',genes_tested=len(d),mapping_status=m.status.value_counts().to_dict(),genes_annotated=len(annotated),pathways_total=len(titles),pathways_evaluable=len(eligible),lipid_pathways=int(catalog.lipid_metabolism.sum()),selected_high_annotated=len(set(d.loc[selected&(d.log2FC>0),'gene'])&annotated),selected_low_annotated=len(set(d.loc[selected&(d.log2FC<0),'gene'])&annotated),primary_ORA_P05=int(((ora.background=='all_tested')&ora.selected).sum()))
(pub/'preparation_validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
