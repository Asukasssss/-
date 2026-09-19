"""Compare published aggregate lineage profiles; no inferential P/q or causal ranking."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def save(d, p):
    d.to_csv(p, sep='\t', index=False, na_rep='NA', lineterminator='\n')


def main(root):
    p = root / 'public'
    cohorts = ['Wu2021', 'Pal2021_reprocessed', 'Reed2024']
    frames = {c: pd.read_csv(p / (c + '_celltype_profiles.tsv'), sep='\t') for c in cohorts}
    allgenes = sorted(frames['Wu2021'].gene.unique())
    assert len(allgenes) == 117
    base = {c: d[(d.partition == 'ALL') & (d.status == 'DONE')].copy() for c,d in frames.items()}
    rows = []
    for gene in allgenes:
        r = {'cancer':'BRCA', 'gene':gene, 'analysis_version':'sc117_v1',
             'interpretation':'descriptive cell source only; RNA is not enzyme activity or metabolite flux'}
        for c,d in base.items():
            z = d[d.gene == gene].sort_values(['effect','celltype'],ascending=[False,True])
            r[c+'_n_evaluable_lineages'] = len(z)
            r[c+'_top_lineage'] = z.iloc[0].celltype if len(z) else 'NOT_EVALUABLE'
            r[c+'_top_mean_log1p10k'] = z.iloc[0].effect if len(z) else np.nan
            r[c+'_top_detection_fraction'] = z.iloc[0].mean_detection_fraction if len(z) else np.nan
            r[c+'_top_n_source_donor_labels'] = z.iloc[0]['n'] if len(z) else 0
            r[c+'_top3_lineages'] = ';'.join(z.celltype.head(3))
        w = base['Wu2021'].query('gene == @gene').set_index('celltype')
        q = base['Pal2021_reprocessed'].query('gene == @gene').set_index('celltype')
        shared = sorted(set(w.index) & set(q.index))
        r['wu_pal_shared_lineages'] = len(shared)
        r['wu_pal_lineage_rank_rho'] = spearmanr(w.loc[shared,'effect'],q.loc[shared,'effect']).statistic if len(shared)>=5 else np.nan
        r['wu_pal_top_lineage_agreement'] = r['Wu2021_top_lineage']==r['Pal2021_reprocessed_top_lineage']
        r['wu_pal_top3_overlap'] = len(set(r['Wu2021_top3_lineages'].split(';')) & set(r['Pal2021_reprocessed_top3_lineages'].split(';')))
        r['status'] = 'DONE' if len(shared)>=5 else 'PARTIAL'
        r['reason'] = 'Pal independent biological source; later atlas reprocessing; annotations not original Pal objects'
        rows.append(r)
    result = pd.DataFrame(rows)
    save(result,p/'gene117_cell_source_comparison.tsv')
    # Sensitivity to a larger within-donor cell minimum uses existing private summaries.
    sensitivity = []
    for c in cohorts:
        d = pd.read_csv(root/'private'/(c+'_donor_profiles.tsv'),sep='\t')
        d['weighted_mean'] = d.mean_log1p10k*d.n_cells
        pooled = d.groupby(['donor','celltype','gene'],as_index=False).agg(n_cells=('n_cells','sum'),
                    weighted_mean=('weighted_mean',lambda x:x.sum(min_count=1)))
        pooled['value'] = pooled.weighted_mean/pooled.n_cells
        for mincells in [20,50]:
            v = pooled[(pooled.n_cells>=mincells) & pooled.value.notna()]
            summary = v.groupby(['gene','celltype'],as_index=False).agg(effect=('value','mean'),n=('value','count'))
            summary = summary[summary.n>=3]
            for gene,g in summary.groupby('gene'):
                g = g.sort_values(['effect','celltype'],ascending=[False,True])
                sensitivity.append({'cohort':c,'gene':gene,'min_cells_per_donor_group':mincells,
                                    'top_lineage':g.iloc[0].celltype,'n_evaluable_lineages':len(g)})
    save(pd.DataFrame(sensitivity),p/'cell_coverage_sensitivity.tsv')
    # Full117 row-standardized display: values are descriptive ranks/scales, not test results.
    lineages = ['Malignant_epithelial','Nonmalignant_epithelial','Fibroblasts','Perivascular',
                'Endothelial','Myeloid','T_NK_cells','B_cells','Plasma_cells','Mast_cells','Adipocytes']
    highlight = ['GPCPD1','PNP','GPI','ASNS','KYNU','PCYT2','ETNK1','NNMT','SORD','SLC7A11',
                 'UPP1','PMM2','PRODH2','LDHA','LDHB']
    for label, genes in [('key_genes',highlight),('all117',allgenes)]:
        fig,axes=plt.subplots(1,3,figsize=(18,max(6,len(genes)*.19)),sharey=True)
        for ax,(c,d) in zip(axes,base.items()):
            z=d.pivot(index='gene',columns='celltype',values='effect').reindex(index=genes,columns=lineages)
            sd=z.std(axis=1).replace(0,np.nan)
            z=z.sub(z.mean(axis=1),axis=0).div(sd,axis=0)
            im=ax.imshow(np.ma.masked_invalid(z.values),aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
            ax.set_title(c+'\nwithin-gene lineage z-score')
            ax.set_xticks(range(len(lineages)));ax.set_xticklabels(lineages,rotation=90,fontsize=8)
            ax.set_yticks(range(len(genes)));ax.set_yticklabels(genes,fontsize=8)
        fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.35,label='Within-study descriptive z-score')
        fig.subplots_adjust(bottom=.23 if label=='key_genes' else .10,right=.87,wspace=.12)
        fig.savefig(p/(label+'_cell_source.png'),dpi=160,bbox_inches='tight');plt.close(fig)
    audits=[json.loads((p/(c+'_validation.json')).read_text()) for c in cohorts]
    summary={'cohorts':[{'cohort':x['cohort'],'n_cells':x['cells_selected'],
                        'n_source_donor_labels':x['source_donor_labels'],
                        'n_genes_in_dictionary':x['genes_unique_in_dictionary']} for x in audits],
             'all117_retained':len(result)==117,
             'wu_pal_same_top_lineage':int(result.wu_pal_top_lineage_agreement.sum()),
             'wu_pal_median_lineage_rank_rho':float(result.wu_pal_lineage_rank_rho.median()),
             'new_hypothesis_tests':0,
             'normal_reference_not_tumor_normal_DE':True,
             'pal_original_normal_preneoplastic_scope':'NOT_RUN_original_annotation_access_blocked',
             'public_individual_expression_exported':False}
    (p/'integration_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True)
    main(parser.parse_args().root)
