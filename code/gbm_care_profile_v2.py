"""CARE 2025 donor-equal sources; primary main analysis, recurrence context only."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd

MAP={'Malignant':'Neoplastic','TAM':'Myeloid','Oligodendrocyte':'Oligodendrocyte','Excitatory neuron':'Neuron','Inhibitory neuron':'Neuron','Astrocyte':'Astrocyte','OPC':'OPC','Lymphocyte':'Lymphoid','Endothel':'Vascular','Pericyte':'Vascular','Other':'Unresolved'}
V='gbm_care_sources_v2'
def save(x,p):x.to_csv(p,sep='\t',index=False,na_rep='NA')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out,commit):
    pub=out/'public/06_EXTERNAL';pub.mkdir(exist_ok=False)
    genes=pd.read_csv(out/'source/genes_unique.tsv',sep='\t');ng=len(genes)
    s=pd.read_csv(out/'private/CARE_sample_metadata.tsv',sep='\t');s.columns=s.columns.str.strip()
    for c in s.select_dtypes('object'):s[c]=s[c].str.strip()
    assert s['Sample ID'].is_unique and s['IDH mutation status'].eq('IDHwt').all()
    s['partition']=np.where(s['Primary or recurrent'].eq('Primary'),'CARE2025_primary','CARE2025_recurrent')
    pretreatment=s['Primary or recurrent'].eq('Primary')&~(s['Radiation before the surgery'].eq('No')&s['Alkylate agents before the surgery'].eq('No'))
    s.loc[pretreatment,'partition']='EXCLUDED_previously_treated_primary'
    assert pretreatment.sum()==1
    m=pd.read_csv(out/'private/CARE_author_cell_metadata.tsv',sep='\t');assert m.CellID.is_unique
    m=m.merge(s,left_on='ID',right_on='Sample ID',validate='many_to_one');assert len(m)==429305 and m['Patient ID'].notna().all()
    fs=sorted((out/'private/sample_sums').glob('*.tsv'));assert len(fs)==121
    raw=pd.concat([pd.read_csv(f,sep='\t') for f in fs]);assert set(raw.cell_type)<=set(MAP)
    raw['celltype']=raw.cell_type.map(MAP)
    raw=raw.merge(s[['Sample ID','Patient ID','partition']],left_on='sample_id',right_on='Sample ID',validate='many_to_one')
    raw=raw[~raw.partition.str.startswith('EXCLUDED')]
    # Multiple recurrent specimens from the same patient are pooled within that patient.
    d=raw.groupby(['partition','Patient ID','celltype','gene'],sort=True).agg(n_cells=('n_cells','sum'),sum_expression=('sum_expression',lambda x:x.sum(min_count=len(x))),sum_detected=('sum_detected',lambda x:x.sum(min_count=len(x))),measured=('measured','all')).reset_index()
    d['mean_expression']=d.sum_expression/d.n_cells;d['detection']=d.sum_detected/d.n_cells
    save(d,out/'private/sc_donor_profiles_private.tsv')
    spec=dict(version=V,code_commit=commit,cohort='CARE 2025 GSE274546;one biological cohort',primary='author Primary status;IDHwt;explicit No preoperative radiation AND No alkylating agents;exclude1 treated primary795nuclei before inspecting expression',recurrence='same cohort contextual sensitivity;not independent replication;pool repeated specimens within patient',minimum_cells_per_donor_category=20,minimum_distinct_donors_per_category=3,minimum_categories_for_rank=2,minimum_detection_for_rank=.01,bootstrap=1000,seed=20260925,normalization='log1p(10000*count/full33538_gene_library);nonnegative integer UMI',annotation='author CellType and post-QC barcode list;no reclustering;Other excluded from source ranking',gene_match='exact unique gene symbol;no imputation',statistical_unit='author Patient ID;donor equal after within-patient pooling',P_q='NA;descriptive only',selected_before_expression='2025 large IDHwt multi-center cohort with processed counts,author annotations,explicit patient and timepoint mapping')
    (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
    profiles=[];ranks=[];registry=[];coverage=[]
    for pi,part in enumerate(['CARE2025_primary','CARE2025_recurrent']):
        sub=d[d.partition.eq(part)];mm=m[m.partition.eq(part)];donors=sorted(sub['Patient ID'].unique());types=sorted(sub.celltype.unique());nd=len(donors);nc=len(types)
        idx=pd.MultiIndex.from_product([donors,types,genes.gene],names=['Patient ID','celltype','gene'])
        a=sub.set_index(['Patient ID','celltype','gene']).reindex(idx)
        means=a.mean_expression.to_numpy().reshape(nd,nc,ng);det=a.detection.to_numpy().reshape(nd,nc,ng);counts=a.n_cells.fillna(0).to_numpy().reshape(nd,nc,ng)
        eligible=(counts>=20)&np.isfinite(means);means[~eligible]=np.nan;det[~eligible]=np.nan
        den=eligible.sum(0);avg=np.divide(np.nansum(means,axis=0),den,out=np.full((nc,ng),np.nan),where=den>=3);detection=np.divide(np.nansum(det,axis=0),den,out=np.full((nc,ng),np.nan),where=den>=3)
        rankavg=avg.copy();rankavg[types.index('Unresolved')]=np.nan
        rankdet=detection.copy();rankdet[types.index('Unresolved')]=np.nan
        rng=np.random.default_rng(20260925+pi);wins=np.zeros((nc,ng));validboot=np.zeros(ng,int)
        for draw in rng.integers(nd,size=(1000,nd)):
            w=np.bincount(draw,minlength=nd);unique=(eligible&(w[:,None,None]>0)).sum(0);denom=(eligible*w[:,None,None]).sum(0)
            av=np.divide((np.nan_to_num(means)*w[:,None,None]).sum(0),denom,out=np.full((nc,ng),np.nan),where=(unique>=3)&(den>=3))
            av[types.index('Unresolved')]=np.nan
            good=np.isfinite(av).sum(0)>=2
            for gi in np.flatnonzero(good):
                tie=np.isclose(av[:,gi],np.nanmax(av[:,gi]),atol=1e-10,rtol=0);wins[tie,gi]+=1/tie.sum();validboot[gi]+=1
        for gi,g in enumerate(genes.gene):
            measured=bool(sub[sub.gene.eq(g)].measured.all());coverage.append(dict(study=part,gene=g,stable_gene_id=genes.stable_gene_id.iloc[gi],status='DONE' if measured else 'NOT_EVALUABLE',reason='exact_unique_symbol_all_samples' if measured else 'missing_in_one_or_more_samples'))
            for ci,ct in enumerate(types):
                vals=means[:,ci,gi];vals=vals[np.isfinite(vals)];ok=np.isfinite(avg[ci,gi])
                profiles.append(dict(cancer='GBM',cohort=part,stage_id='06_EXTERNAL',run_id=out.name,analysis_version=V,analysis_type='donor_equal_cell_source',metabolite_key='NA',metabolite_name='NA',gene=g,unit='author_patient',n=len(vals),n_reference=np.nan,effect_type='donor_equal_mean_log1pCP10k',effect=avg[ci,gi],ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family='SC_SOURCE_DESCRIPTIVE',family_n_evaluable=np.nan,status='DONE' if ok else 'NOT_EVALUABLE',reason='descriptive_only' if ok else 'missing_gene_or_insufficient_patient_coverage',source_id='GSE274546',stable_gene_id=genes.stable_gene_id.iloc[gi],partition=part,celltype=ct,n_source_labels_total=nd,n_source_labels_eligible=len(vals),n_cells_total=int(counts[:,ci,gi].sum()),n_cells_eligible=int(counts[:,ci,gi][eligible[:,ci,gi]].sum()),mean_detection_fraction=detection[ci,gi],median_donor_mean=float(np.median(vals)) if ok else np.nan,q25=float(np.quantile(vals,.25)) if ok else np.nan,q75=float(np.quantile(vals,.75)) if ok else np.nan,annotation_origin='CARE author CellType 2025-01-08',normalization=spec['normalization']))
            rr=dict(study=part,gene=g,stable_gene_id=genes.stable_gene_id.iloc[gi],status='NOT_EVALUABLE',reason='fewer_than2_categories_or_detection_below1percent',top_celltype='NA',runner_celltype='NA',top_gap=np.nan,bootstrap_top_frequency=np.nan,bootstrap_valid=int(validboot[gi]),rank_tie=False)
            ix=np.flatnonzero(np.isfinite(rankavg[:,gi]))
            if len(ix)>=2 and np.nanmax(rankdet[:,gi])>=.01:
                order=ix[np.argsort(-rankavg[ix,gi])];tie=ix[np.isclose(rankavg[ix,gi],rankavg[order[0],gi],rtol=0,atol=1e-10)]
                rr.update(status='DONE',reason='expression_source_only;not_function',top_celltype=';'.join(types[j] for j in tie),runner_celltype=types[order[1]],top_gap=rankavg[order[0],gi]-rankavg[order[1],gi],bootstrap_top_frequency=wins[tie,gi].sum()/validboot[gi] if validboot[gi] else np.nan,rank_tie=len(tie)>1)
            ranks.append(rr)
        registry.append(dict(study=part,cells=len(mm),patients=nd,samples=mm.ID.nunique(),unresolved_cells=int(mm.CellType.eq('Other').sum()),genes_measured=int(sum(x['status']=='DONE' for x in coverage if x['study']==part)),candidate_genes=ng,assay='10x Genomics snRNA-seq',source='GSE274546',publication_year=2025,annotation_origin='CARE author CellType',independent_validation=False))
    prof=pd.DataFrame(profiles);rank=pd.DataFrame(ranks)
    for name,tab in [('sc_celltype_profiles',prof),('sc_source_stability',rank),('sc_gene_coverage',pd.DataFrame(coverage)),('sc_dataset_registry',pd.DataFrame(registry)),('sc_annotation_map',pd.DataFrame([dict(original=k,coarse=v,ranking_eligible=v!='Unresolved') for k,v in MAP.items()]))]:save(tab,pub/(name+'.tsv'))
    comp=rank[rank.study.eq('CARE2025_primary')][['gene','top_celltype','bootstrap_top_frequency','status']].merge(rank[rank.study.eq('CARE2025_recurrent')][['gene','top_celltype','bootstrap_top_frequency','status']],on='gene',suffixes=('_primary','_recurrent'),validate='one_to_one')
    comp['same_top']=np.where(comp.status_primary.eq('DONE')&comp.status_recurrent.eq('DONE'),comp.top_celltype_primary.eq(comp.top_celltype_recurrent).astype(float),np.nan);comp['comparison_type']='same_cohort_partitions_not_independent_replication';save(comp,pub/'sc_partition_comparison.tsv')
    audits=pd.concat([pd.read_csv(f,sep='\t') for f in (out/'private/sample_sums').glob('*.audit')]);assert audits.author_qc_cells.sum()==429305 and (audits.author_qc_cells==audits.exact_barcode_matches).all()
    val=dict(status='DONE',samples=121,patients=59,author_qc_cells=429305,excluded_previously_treated_primary_samples=1,excluded_previously_treated_primary_nuclei=795,unique_barcodes=True,exact_metadata_join=True,source_counts_cells=int(audits.source_cells.sum()),all142_retained=True,partitions=registry,minimum_full_library=int(audits.min_library.min()),P_q_all_missing=bool(prof.p_value.isna().all() and prof.q_value.isna().all()),new_clustering=False,independent_replication=False)
    (pub/'validation.json').write_text(json.dumps(val,indent=2));print(json.dumps(val),flush=True)
    save(pd.DataFrame([dict(source_id=p.name,server_relative_path=str(p.relative_to(out)),sha256=sha(p)) for p in [out/'source/celltype_meta_data_2025_01_08.RDS',out/'source/spitzer_supptable1.xlsx',out/'source/genes_unique.tsv',out/'private/CARE_source_file_manifest.tsv']]),pub/'source_manifest.tsv')
    save(pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in pub.iterdir() if p.is_file()]),pub/'checksums.tsv')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--commit',required=True);a=p.parse_args();main(a.out,a.commit)
