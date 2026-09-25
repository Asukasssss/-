"""Reuse completed tests; user-directed nominal-P selection and ORA."""
import argparse,json,sys,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PRIOR=ROOT/'results/collaborative/COAD/B/20260925T151300Z_lypla1_highlow_v1'
sys.path.insert(0,str(ROOT/'results/collaborative/COAD/B/20260925T140257Z_uckl1_cna_v1'))
from run_uckl1_cna_v1 import sha,PREFIX
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--min-logfc',type=float,default=.5);ap.add_argument('--analysis-version',default='COAD_LYPLA1_highlow_nominal_v2');a=ap.parse_args();out=a.out;pub=out/'public'
assert a.min_logfc>0
assert out.parent==ROOT/'results/collaborative/COAD/B' and (out/'.running').is_dir();pub.mkdir()
d=pd.read_csv(PRIOR/'public/results.tsv',sep='\t');base=d.copy()
assert len(d)==11058 and d.gene.is_unique and 'LYPLA1' not in set(d.gene)
d=d.drop(columns=[c for c in d if c in ['q_value','depth_q_value','selected']])
d['selected']=(d.p_value<.05)&(d.log2FC.abs()>=a.min_logfc)
d['depth_selected']=(d.depth_p_value<.05)&(d.depth_log2FC.abs()>=a.min_logfc)
common=dict(cancer='COAD',cohort='Uhlitz_GSE166555',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=a.analysis_version,
    analysis_type='donor_paired_LYPLA1_positive_high_vs_low',metabolite_key='NA',metabolite_name='NA',unit='donor',n=8,n_reference=8,
    effect_type='edgeR_log2FC_high_vs_low',ci_lower='NA',ci_upper='NA',q_value='NA',test_family='nominal_P_no_multiple_testing_threshold',family_n_evaluable=len(d),
    status='DONE',reason='EXPLORATORY_NOMINAL_P;PRIOR_Q_ARCHIVED_NOT_USED',source_id='Uhlitz_GSE166555_author_CNA')
for k,v in common.items():d[k]=v
d['effect']=d.log2FC
d=d[PREFIX+[c for c in d if c not in PREFIX]]
d.to_csv(pub/'results.tsv',sep='\t',index=False);d[d.selected].to_csv(pub/'selected_genes.tsv',sep='\t',index=False)
pathways={};titles={}
for line in (PRIOR/'reactome.gmt').read_text().splitlines():
    z=line.split('\t');assert z[1].startswith('R-HSA-')
    members=set(z[2:])&set(d.gene)
    if 15<=len(members)<=500:pathways[z[1]]=members;titles[z[1]]=z[0]
rows=[];M=len(d)
for direction in ['higher_in_LYPLA1_high','lower_in_LYPLA1_high']:
    selected=set(d.loc[d.selected & (d.direction==direction),'gene']);N=len(selected)
    for k,genes in pathways.items():
        hits=sorted(genes&selected);K=len(genes);n=len(hits)
        rows.append(dict(pathway_id=k,pathway=titles[k],direction=direction,background_genes=M,selected_genes=N,pathway_tested_genes=K,overlap=n,
            fold_enrichment=n/N/(K/M) if N else 'NA',p_value=float(hypergeom.sf(n-1,M,K,N)) if N else 1.,overlap_genes=';'.join(hits),
            status='DONE' if N else 'NOT_EVALUABLE',reason='NOMINAL_P_EXPLORATORY' if N else 'NO_SELECTED_GENES'))
ora=pd.DataFrame(rows).sort_values(['p_value','pathway_id']);ora['selected']=ora.p_value<.05;ora.to_csv(pub/'reactome_ORA.tsv',sep='\t',index=False)
g=pd.read_csv(PRIOR/'public/reactome_GSEA.tsv',sep='\t').drop(columns=['padj'])
g['selected']=g.pval<.05;g=g.sort_values(['pval','pathway_id']);g.to_csv(pub/'reactome_GSEA.tsv',sep='\t',index=False)
for f in ['admission.json','group_summary.tsv','gene_testability.tsv','software_versions.txt']:shutil.copyfile(PRIOR/'public'/f,pub/f)
spec=json.loads((PRIOR/'analysis_spec.json').read_text())
spec.update(analysis_version=a.analysis_version,user_change='User explicitly requested no FDR; nominal P and requested effect threshold',
    min_abs_log2FC=a.min_logfc,DE_selection='P < 0.05 and absolute log2FC >= '+str(a.min_logfc),gene_family='All 11058 tested genes retained; no new BH; prior Q values archived only',
    ORA='One-sided hypergeometric; actual tested-gene background; nominal P<0.05, no FDR threshold',
    GSEA='Reuse prior NES, raw P and leadingEdge; nominal P<0.05; no new permutations',
    prior_run=PRIOR.name,prior_numeric_lock='55e97756d24236f05bd3dca23357fa94f5460ef7')
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
sources=pd.read_csv(PRIOR/'public/source_manifest.tsv',sep='\t')
for f in ['public/results.tsv','public/reactome_GSEA.tsv','analysis_spec.json','reactome.gmt']:
    p=PRIOR/f;sources.loc[len(sources)]=dict(source_id='prior_'+f,source_path=str(p),sha256=sha(p),url='prior frozen batch')
sources.to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
summary=dict(status='PASS',min_abs_log2FC=a.min_logfc,gene_P_below_05=int((d.p_value<.05).sum()),selected_high=int((d.selected&(d.log2FC>0)).sum()),selected_low=int((d.selected&(d.log2FC<0)).sum()),
    depth_selected=int(d.depth_selected.sum()),primary_and_depth_selected=int((d.selected&d.depth_selected).sum()),
    primary_selected_same_direction_depth_P05=int((d.selected&(d.depth_p_value<.05)&d.depth_direction_agrees).sum()),
    pathways_tested=len(pathways),ORA_P05=int(ora.selected.sum()),GSEA_P05=int(g.selected.sum()),
    GSEA_high_P05=int((g.selected&(g.NES>0)).sum()),GSEA_low_P05=int((g.selected&(g.NES<0)).sum()),
    primary_effects_and_P_exactly_reused=bool(np.array_equal(d.p_value,base.p_value)&np.array_equal(d.log2FC,base.log2FC)),
    no_new_differential_tests=True,no_new_GSEA=True,no_new_FDR=True,ORA_recomputed_for_nominal_gene_lists=True,prior_Q_unchanged=True)
(pub/'validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
(out/'DONE').write_text('DONE\n');(out/'.running').rmdir()
