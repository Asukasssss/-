"""Run only in a fresh server165 PDAC/B directory. Never writes prior results.

Inputs in source/: panel.json, gene_aliases.json, analysis_spec.json, templates/,
and historical_manifest.tsv. Helpers are copied beside this script.
--mode internal or source; each mode must use its own new run directory.
"""
import argparse,json,platform,traceback,itertools
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import sc_source_three_cohorts_v1 as readers
from internal_paired51_v2 import MP,XP,RP,AP,EXPECTED,validate_sources,resolve_genes
from patient_first_round import validate_rna_labels
from paired_metabolites_v1 import pair_join
from sop_v3_math import paired_t,bh_evaluable,cellwise_profiles,source_rank_bootstrap,seed_for

ROOT=readers.ROOT
VERSION='PDAC_CAMP_source_SOP_v3'
COHORTS=['GSE263733','GSE278688','GSE242230']
def write(run,area,name,df):
    template=run/'source/templates'/name
    if template.exists():
        prefix=list(pd.read_csv(template,sep='\t').columns)
        for c in prefix:
            if c not in df:df[c]=np.nan
        df=df[prefix+[c for c in df if c not in prefix]]
    df.to_csv(run/area/name,sep='\t',index=False,na_rep='NA')
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')

def internal(run,panel,aliases):
    validate_sources(EXPECTED,[MP,XP,RP,AP]);j,pairs=pair_join(pd.read_csv(MP,dtype=str),pd.read_excel(AP,dtype=str))
    rna=pd.read_csv(RP,index_col=0);validate_rna_labels(rna)
    assert j.RNAID.isin(rna.columns).all() and len(j)==39 and len(pairs)==11
    tumors=j[(j.TN=='Tumor')&j['Paired Sample ID'].notna()]
    assert len(tumors)==21 and tumors['Paired Sample ID'].is_unique
    ids=set(pairs.author_pair_id);tp=j[(j.TN=='Tumor')&j['Paired Sample ID'].isin(ids)].sort_values('Paired Sample ID');nn=j[(j.TN=='Normal')&j['Paired Sample ID'].isin(ids)].sort_values('Paired Sample ID')
    assert tp['Paired Sample ID'].tolist()==nn['Paired Sample ID'].tolist()
    # Independently emit the RNA links from exact author RNA IDs, not metabolite order.
    rr=[]
    for (_,t),(_,n) in zip(tp.iterrows(),nn.iterrows()):rr.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','pair_id':t['Paired Sample ID'],'patient_id':t['Paired Sample ID'],'tumor_RNA_id':t.RNAID,'normal_RNA_id':n.RNAID,'pair_source':'Explicit author GSE62452 pairing XLSX;exact GSM join','status':'DONE','reason':'Author case key,not protected clinical identity recertification'})
    write(run,'private','RNA_pairs_private.tsv',pd.DataFrame(rr))
    mm=pairs.rename(columns={'author_pair_id':'pair_id','tumor_metab_id':'tumor_metabolomics_id','normal_metab_id':'normal_metabolomics_id'});mm['patient_id']=mm.pair_id;mm['cancer']='PDAC';mm['cohort']='CAMP_PDAC_GSE62452';mm['pair_source']='Author explicit GSE62452 pairing';mm['status']='DONE';mm['reason']='No ID suffix inference';write(run,'private','metabolite_pairs_private.tsv',mm)
    rows=[]
    for _,r in j.iterrows():rows.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','patient_id':r['Paired Sample ID'],'specimen_id':r.CommonID,'tissue_label_author':r.TN,'tissue_label_verified':r['T/N'],'metabolomics_id':r.MetabID,'RNA_id':r.RNAID,'pair_id':r['Paired Sample ID'],'link_level':'Author mapped sample;aliquot level unknown','link_source':'CAMP mapping and explicit author pair XLSX','duplicate_group':np.nan,'decision':'INCLUDE_AUTHOR_KEY' if pd.notna(r['Paired Sample ID']) else 'EXCLUDE_FROM_INDEPENDENT_UNIT_ANALYSIS','reason':'Author mapped IDs unique;clinical identity not re-certified'})
    write(run,'private','sample_identity_audit_private.tsv',pd.DataFrame(rows))
    tr=[]
    for _,r in tumors.iterrows():tr.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','analysis_unit_id':r['Paired Sample ID'],'patient_id':r['Paired Sample ID'],'specimen_id':r.CommonID,'metabolomics_id':r.MetabID,'RNA_id':r.RNAID,'link_level':'Author corresponding sample;aliquot unknown','unit_verification':'Distinct author case keys','source_id':'CAMP+GSE62452_pair_table','status':'DONE','reason':'Not independently confirmed clinical identity'})
    write(run,'private','tumor_multiomics_map_private.tsv',pd.DataFrame(tr))
    resolved,why=resolve_genes(panel['genes'],rna.index,aliases);out=[];coverage=[];current=set(panel['current_genes'])
    for gene in panel['genes']:
        label=resolved[gene];stable=panel['stable_gene_ids'][gene];seed=seed_for('GSE62452','RNA_PAIRED',stable)
        a=rna.loc[label,tp.RNAID].to_numpy(float) if label else np.full(11,np.nan);b=rna.loc[label,nn.RNAID].to_numpy(float) if label else np.full(11,np.nan)
        z=paired_t(a-b,seed);in_current=gene in current
        row={'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','stage_id':'03_PATIENT','run_id':run.name,'analysis_version':VERSION,'analysis_type':'paired_RNA_t','gene':gene,'stable_gene_id':stable,'unit':'author_case_pair','n_reference':11,'effect_type':'paired_mean_T_minus_N_author_microarray_scale','test_family':'RNA_PAIRED_CURRENT' if in_current else 'RNA_PAIRED_HISTORY_SUPPLEMENT','source_id':'GSE62452_author_processed_RNA','expression_platform':'Affymetrix_Human_Gene_1.0_ST','expression_scale':'author_processed_continuous_no_added_transform','n_pairs_total':11,'model_formula':'paired_T_minus_N','p_method':'paired_t','ci_method':'Student_t_df_n_minus1_for_mean_delta','probe_mapping_status':why[gene],'current_pool':in_current,'history_only':not in_current,'seed':seed,'bootstrap_B':4000,'rna_label':label,**z}
        row.update(effect=z['mean_delta'],n_pairs_used=z['n'],up_fraction=z['n_up']/z['n'] if z['n'] else np.nan)
        if not label:row['reason']=why[gene]
        out.append(row);coverage.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','stable_gene_id':stable,'gene':gene,'measured_feature_id':label,'mapping_policy':'Exact canonical or single reviewed noncolliding alias','mapping_version':VERSION,'status':'DONE' if label else 'NOT_EVALUABLE','reason':why[gene],'source_id':'GSE62452_author_processed_RNA'})
    df=pd.DataFrame(out)
    for family,d in df.groupby('test_family'):
        df.loc[d.index,'q_value']=bh_evaluable(d.p_value);df.loc[d.index,'family_n_evaluable']=d.p_value.notna().sum();df.loc[d.index,'family_n_planned']=len(d)
    df['evidence_tier']=np.where(df.p_value.isna(),'NOT_EVALUABLE',np.where(df.q_value<.05,'FDR_SUPPORTED',np.where(df.p_value<.05,'NOMINAL_EXPLORATORY','NOT_SUPPORTED_THIS_TEST')))
    write(run,'public','paired_RNA.tsv',df[df.current_pool]);write(run,'public','paired_RNA_history_supplement.tsv',df[~df.current_pool]);write(run,'public','RNA_P005_view.tsv',df[df.current_pool&(df.p_value<.05)]);write(run,'public','RNA_identity_coverage.tsv',pd.DataFrame(coverage))
    # Marginal availability and descriptions for all307, never export per-pair values.
    raw=pd.read_excel(XP,sheet_name='data',index_col=0);t=pd.read_excel(XP,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);n=pd.read_excel(XP,sheet_name='metabo_imputed_filtered_Normal',index_col=0)
    for v in [raw,t,n]:validate_rna_labels(v)
    desc=[]
    for f in sorted(set(t.index)&set(n.index)):
        tv=t.loc[f,pairs.tumor_metab_id].to_numpy(float);nv=n.loc[f,pairs.normal_metab_id].to_numpy(float);mt=np.isfinite(raw.loc[f,pairs.tumor_metab_id].to_numpy(float));mn=np.isfinite(raw.loc[f,pairs.normal_metab_id].to_numpy(float))
        for mode,mask in [('primary',np.isfinite(tv)&np.isfinite(nv)),('both_preimputation_available',mt&mn&np.isfinite(tv)&np.isfinite(nv))]:
            d=(tv-nv)[mask];desc.append({'metabolite_name':f,'analysis_type':mode,'n':len(d),'n_nonzero':int((d!=0).sum()),'n_up':int((d>0).sum()),'n_down':int((d<0).sum()),'n_equal':int((d==0).sum()),'mean_delta':float(d.mean()) if len(d) else np.nan,'median_delta':float(np.median(d)) if len(d) else np.nan,'missing_rate_tumor':float((~mt).mean()),'missing_rate_normal':float((~mn).mean()),'n_both_available':int((mt&mn).sum()),'mask_semantics':'author_available_value_mask'})
    write(run,'public','metabolite_description_supplement.tsv',pd.DataFrame(desc))
    return {'RNA':{f:{'planned':len(d),'evaluable':int(d.p_value.notna().sum()),'P_lt005':int((d.p_value<.05).sum()),'q_lt005':int((d.q_value<.05).sum())} for f,d in df.groupby('test_family')},'source_paths':[str(x) for x in EXPECTED]}

def source(run,panel,aliases):
    genes=panel['genes'];allprofiles=[];tops=[];coverage=[];annotations=[];checks=[]
    search=sorted(set(genes)|{a for g in genes for a in aliases.get(g,[])})
    for cohort,reader in zip(COHORTS,[readers.dense263,readers.sparse278,readers.sparse242]):
        print('EXTRACT',cohort,flush=True);meta,raw,lib,present=reader(search);resolved,why=resolve_genes(genes,present,aliases)
        look={g:i for i,g in enumerate(search)};x=np.zeros((len(genes),len(meta)),dtype=np.int64);found=np.zeros(len(genes),bool)
        for i,g in enumerate(genes):
            if resolved[g]:x[i]=raw[look[resolved[g]]];found[i]=True
            coverage.append({'cancer':'PDAC','cohort':cohort,'stable_gene_id':panel['stable_gene_ids'][g],'gene':g,'source_feature_id':resolved[g],'source_symbol':resolved[g],'match_method':why[g],'status':'DONE' if resolved[g] else 'NOT_EVALUABLE','reason':'Source symbol linked to stable UniProt identity;no highest-expression match' if resolved[g] else why[g],'source_id':cohort})
        del raw
        meta=meta.copy();meta['coarse']=meta.celltype.map(readers.broad);meta['position']=np.arange(len(meta));private=[]
        assert meta[['unit','celltype']].notna().all().all() and meta.index.is_unique
        for label in sorted(meta.celltype.unique()):annotations.append({'cancer':'PDAC','cohort':cohort,'author_celltype':label,'coarse_celltype':readers.broad(label),'annotation_origin':'Author supplied','mapping_evidence':'Frozen preexisting PDAC coarse mapping;no output-based relabeling','version':VERSION,'status':'DONE','reason':'Ductal/epithelial does not confirm malignant CNV'})
        grouped=list(meta.groupby(['unit','coarse'],sort=True));summaries=cellwise_profiles(x,lib,[g.position.to_numpy() for _,g in grouped],found)
        for ((unit,ct),group),(nc,mean,det) in zip(grouped,summaries):
            idx=group.position.to_numpy()
            for i,g in enumerate(genes):private.append({'cancer':'PDAC','cohort':cohort,'donor_id':str(unit),'partition':'primary_tumor','celltype':ct,'stable_gene_id':panel['stable_gene_ids'][g],'gene':g,'n_cells':nc,'mean_log1p10k':mean[i],'detection_fraction':det[i],'raw_count':int(x[i,idx].sum()) if found[i] else np.nan,'all_gene_library_sum':int(lib[idx].sum()),'source_id':cohort})
        pb=pd.DataFrame(private);write(run,'private',cohort+'_sc_donor_profiles_private.tsv',pb)
        units=sorted(meta.unit.astype(str).unique());cats=sorted(meta.coarse.unique())
        for gene in genes:
            d=pb[pb.gene==gene];eligible=d[(d.n_cells>=20)&d.mean_log1p10k.notna()]
            mx=eligible.pivot(index='donor_id',columns='celltype',values='mean_log1p10k').reindex(index=units,columns=cats).to_numpy(float);det=eligible.pivot(index='donor_id',columns='celltype',values='detection_fraction').reindex(index=units,columns=cats).to_numpy(float)
            seed=seed_for(cohort,'SC_SOURCE_DESCRIPTIVE',panel['stable_gene_ids'][gene]);s=source_rank_bootstrap(mx,det,seed)
            for k,ct in enumerate(cats):
                dd=d[d.celltype==ct];ee=eligible[eligible.celltype==ct];good=bool(s['eligible'][k]);nd=len(ee)
                allprofiles.append({'cancer':'PDAC','cohort':cohort,'stage_id':'06_EXTERNAL','run_id':run.name,'analysis_version':VERSION,'analysis_type':'sc_expression_source','gene':gene,'stable_gene_id':panel['stable_gene_ids'][gene],'unit':'author_source_label','n':nd,'n_reference':len(units),'effect_type':'equal_donor_mean_cellwise_log1p10k','effect':s['mean'][k] if good else np.nan,'p_value':np.nan,'q_value':np.nan,'test_family':'SC_SOURCE_DESCRIPTIVE','family_n_evaluable':np.nan,'status':'DONE' if good else 'NOT_EVALUABLE','reason':'Donor-equal descriptive expression' if good else 'GENE_MISSING_OR_LT3_ELIGIBLE_LABELS','source_id':cohort,'partition':'primary_tumor','celltype':ct,'n_source_labels_total':len(dd),'n_source_labels_eligible':nd,'n_cells_total':int(dd.n_cells.sum()),'n_cells_eligible':int(ee.n_cells.sum()),'mean_detection_fraction':s['mean_detection'][k] if good else np.nan,'median_donor_mean':float(ee.mean_log1p10k.median()) if good else np.nan,'q25':float(ee.mean_log1p10k.quantile(.25)) if good else np.nan,'q75':float(ee.mean_log1p10k.quantile(.75)) if good else np.nan,'annotation_origin':'Author labels,predeclared coarse mapping','normalization':'mean_donor(mean_cell(log1p(10000*counts/all_gene_library)))','bootstrap_first_frequency':s['frequency'][k]})
            top=s['top'];runner=s['runner'];topnames=';'.join(cats[k] for k in top) if top else 'NOT_EVALUABLE';rnames=';'.join(cats[k] for k in runner) if runner else 'NOT_EVALUABLE';gap=float(s['mean'][top[0]]-s['mean'][runner[0]]) if top and runner else np.nan
            paired=mx[:,top[0]]-mx[:,runner[0]] if len(top)==len(runner)==1 else np.array([]);paired=paired[np.isfinite(paired)]
            tops.append({'cancer':'PDAC','cohort':cohort,'gene':gene,'stable_gene_id':panel['stable_gene_ids'][gene],'partition':'primary_tumor','n_evaluable_celltypes':int(s['eligible'].sum()),'top_celltype':topnames,'runner_celltype':rnames,'top_gap':gap,'max_detection':float(np.nanmax(s['mean_detection'][s['eligible']])) if s['eligible'].any() else np.nan,'tie_status':'TIE' if len(top)>1 else 'UNIQUE' if top else 'NOT_EVALUABLE','bootstrap_top_frequency':float(s['frequency'][top[0]]) if len(top)==1 else np.nan,'bootstrap_valid':s['valid'],'bootstrap_total':1000,'paired_source_labels':len(paired),'paired_positive_fraction':float((paired>0).mean()) if len(paired) else np.nan,'paired_mean_delta':float(paired.mean()) if len(paired) else np.nan,'seed':seed,'status':s['status'],'reason':s['reason'],'source_id':cohort})
        checks.append({'cohort':cohort,'cells':len(meta),'source_labels':len(units),'genes':len(genes),'matrix_resolved':int(found.sum())});print('SUMMARIZED',cohort,flush=True)
        del x,pb
    prof=pd.DataFrame(allprofiles);top=pd.DataFrame(tops);cross=[]
    for gene in genes:
        for ca,cb in itertools.combinations(COHORTS,2):
            ta=top[(top.gene==gene)&(top.cohort==ca)].iloc[0];tb=top[(top.gene==gene)&(top.cohort==cb)].iloc[0];pa=prof[(prof.gene==gene)&(prof.cohort==ca)&(prof.status=='DONE')];pb=prof[(prof.gene==gene)&(prof.cohort==cb)&(prof.status=='DONE')];common=set(pa.celltype)&set(pb.celltype)
            shared=[]
            for d in [pa,pb]:
                sub=d[d.celltype.isin(common)];ok=len(sub)>=2 and sub.mean_detection_fraction.max()>=.01
                shared.append(';'.join(sorted(sub.loc[np.isclose(sub.effect,sub.effect.max(),rtol=1e-10,atol=1e-12),'celltype'])) if ok else 'NOT_EVALUABLE')
            valid=ta.status==tb.status=='DONE' and ta.tie_status==tb.tie_status=='UNIQUE';same=bool(valid and ta.top_celltype==tb.top_celltype)
            cross.append({'cancer':'PDAC','gene':gene,'stable_gene_id':panel['stable_gene_ids'][gene],'study_A':ca,'study_B':cb,'top_A':ta.top_celltype,'top_B':tb.top_celltype,'n_common_categories':len(common),'same_top_all_categories':same if valid else np.nan,'top_A_shared':shared[0],'top_B_shared':shared[1],'same_top_shared_categories':shared[0]==shared[1] if all(x!='NOT_EVALUABLE' and ';' not in x for x in shared) else np.nan,'bootstrap_top_A':ta.bootstrap_top_frequency,'bootstrap_top_B':tb.bootstrap_top_frequency,'same_top_both_bootstrap_ge080':bool(same and ta.bootstrap_top_frequency>=.8 and tb.bootstrap_top_frequency>=.8) if valid else np.nan,'status':'DONE' if valid else 'NOT_EVALUABLE','reason':'Descriptive label;clinical cross-study identity de-duplication remains unverified'})
    for name,df in [('sc_celltype_profiles.tsv',prof),('sc_source_stability.tsv',top),('sc_cross_study.tsv',pd.DataFrame(cross)),('sc_gene_coverage.tsv',pd.DataFrame(coverage)),('sc_annotation_map.tsv',pd.DataFrame(annotations))]:write(run,'public',name,df)
    return {'cohorts':checks,'genes':len(genes),'profile_rows':len(prof),'source_paths':[str(x) for x in readers.SOURCE_FILES]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['internal','source'],required=True);ap.add_argument('--code-commit',required=True);args=ap.parse_args();run=Path(__file__).resolve().parent
    assert run.parent==ROOT/'results/collaborative/PDAC/B','Fresh isolated server directory required'
    (run/'private').mkdir(exist_ok=True);(run/'public').mkdir(exist_ok=True)
    with (run/'.running').open('x') as f:f.write(args.mode)
    try:
        panel=json.loads((run/'source/panel.json').read_text());aliases=json.loads((run/'source/gene_aliases.json').read_text());spec=json.loads((run/'source/analysis_spec.json').read_text())
        assert len(panel['genes'])==len(set(panel['genes']))==687 and len(panel['current_genes'])==250
        if args.mode=='source':
            manifest=pd.read_csv(run/'source/historical_manifest.tsv',sep='\t')
            for _,r in manifest.iterrows():
                if '/data/candidates/' in r.path and readers.sha(Path(r.path))!=r.sha256:raise ValueError('Historical SC input hash changed: '+r.path)
        result=internal(run,panel,aliases) if args.mode=='internal' else source(run,panel,aliases)
        paths=[Path(x) for x in result.pop('source_paths')]+list((run/'source').glob('*.json'))+list(run.glob('*.py'))
        write(run,'public','source_manifest.tsv',pd.DataFrame([{'source_id':p.name,'path_or_url':str(p),'sha256':readers.sha(p),'access_scope':'SERVER_INPUT'} for p in paths]))
        spec.update(code_commit=args.code_commit,execution_mode=args.mode,software_versions={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__})
        save(run/'public/analysis_spec.json',spec);save(run/'public/summary.json',result);save(run/'NUMERICS_DONE.json',result);(run/'.running').unlink()
    except Exception:(run/'FAILED.txt').write_text(traceback.format_exc());raise

if __name__=='__main__':main()
