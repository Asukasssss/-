"""Bounded FUSCC/CPTAC coverage and author-patient-linked exploratory association."""
import argparse,json,hashlib,datetime
from pathlib import Path
import requests,numpy as np,pandas as pd
from brca_pycr1_context_v1 import infer,bh,save
from brca_sc117_profile_v1 import sha256

def main(root,project):
    old=project/'data/candidates/BRCA_1MNA_followup_v0.1';prior=project/'results/BRCA_1MNA_followup_v0.1/tables'
    workbook=project/'data/candidates/brca_four_axes_20260908/FUSCC_41422_2022_614_MOESM12_ESM.xlsx'
    ann=pd.read_csv(project/'results/CAMP_first_analysis_20260910/tables/FUSCC_author_polar_annotation_594.tsv',sep='\t')
    assert len(ann)==594 and ann.fuscc_peak.is_unique
    requested=[('GPCPD1','C00670','HMDB0000086'),('GPI','C00668','HMDB0001401')]
    hits=[]
    for gene,kegg,hmdb in requested:
        a=ann[ann.fuscc_kegg.eq(kegg)]
        hits.append(dict(cohort='FUSCC_TNBC',gene=gene,metabolite_key='KEGG:'+kegg,n_annotation_hits=len(a),status='DONE' if len(a)==1 and a.iloc[0].fuscc_hmdb==hmdb else 'NOT_EVALUABLE',reason='author_annotation_only_not_spectral_reconfirmation' if len(a)==1 else 'no_exact_feature_in594_polar_annotation',peak=a.iloc[0].fuscc_peak if len(a)==1 else '',author_name=a.iloc[0].fuscc_author_name if len(a)==1 else '',identification_level=a.iloc[0].fuscc_identification_level if len(a)==1 else ''))
    save(pd.DataFrame(hits),root/'public/external_metabolite_coverage.tsv')
    sm=pd.DataFrame(json.loads((old/'Fudan_samples.json').read_text()));assert sm.sampleId.is_unique and sm.patientId.is_unique and sm.sampleType.eq('Primary Solid Tumor').all()
    base='https://data.3steps.cn/cdataportal/api';receipts=[]
    def fetch(endpoint,name,method='GET',body=None):
        p=root/'source'/name
        if not p.exists():
            r=requests.request(method,base+endpoint,json=body,timeout=45);r.raise_for_status();p.write_text(json.dumps(r.json()))
        receipts.append(dict(url=base+endpoint,method=method,path=str(p),sha256=sha256(p),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),request_scope='two target genes; cached author primary tumor sample IDs; individual data server only'))
        return json.loads(p.read_text())
    profiles=fetch('/studies/FUSCC_BRCA_2022/molecular-profiles','profiles.json')
    profile='FUSCC_BRCA_2022_mrna_seq_tpm';assert any(x['molecularProfileId']==profile for x in profiles)
    geneinfo=[fetch('/genes/'+g,g+'_gene.json') for g,_,_ in requested]
    body={'entrezGeneIds':[g['entrezGeneId'] for g in geneinfo],'sampleIds':sm.sampleId.tolist()}
    ex=pd.DataFrame(fetch('/molecular-profiles/'+profile+'/molecular-data/fetch?projection=DETAILED','two_gene_TPM.json','POST',body))
    ex['gene']=ex.gene.map(lambda x:x['hugoGeneSymbol']);assert not ex.duplicated(['patientId','gene']).any();assert ex.gene.isin(['GPCPD1','GPI']).all();assert (ex.value>=0).all()
    chk=ex.merge(sm[['sampleId','patientId']],on='sampleId',suffixes=('_expression','_sample'),validate='many_to_one');assert chk.patientId_expression.eq(chk.patientId_sample).all()
    wide=ex.pivot(index='patientId',columns='gene',values='value').reset_index().rename(columns={'patientId':'patient_id'})
    cl=pd.read_csv(prior/'FUSCC_1MNA_tumor_clinical_values.tsv',sep='\t');rm=pd.read_csv(prior/'FUSCC_author_RNA_sample_patient_mapping.tsv',sep='\t');rm=rm[rm.tissue_type.eq('tumor')]
    join=pd.read_csv(prior/'FUSCC_metabolite_RNA_patient_join_audit.tsv',sep='\t')
    assert cl.patient_id.is_unique and cl.metabolite_sample_id.is_unique and rm.patient_id.is_unique and join.patient_id.is_unique
    z=cl[['patient_id','metabolite_sample_id','mRNA_Subtype']].merge(wide,on='patient_id',validate='one_to_one').merge(rm[['patient_id','sample_id']],on='patient_id',validate='one_to_one')
    z=z.merge(join[['patient_id','author_RNA_sample_id','metabolite_sample_id','join_status']],on=['patient_id','metabolite_sample_id'],validate='one_to_one')
    assert z.sample_id.eq(z.author_RNA_sample_id).all() and z.join_status.eq('PASS_EXPLICIT_AUTHOR_PATIENT_MAPPING').all()
    mat=pd.read_excel(workbook,sheet_name='Table S1',header=1);mat.columns=[str(x).strip() for x in mat.columns];assert mat.Peak.is_unique
    mat=mat.set_index('Peak');assert set(z.metabolite_sample_id).issubset(mat.columns)
    rows=[]
    for h in hits:
        if h['status']=='DONE':z[h['gene']+'_metabolite']=pd.to_numeric(mat.loc[h['peak'],z.metabolite_sample_id],errors='coerce').values
    save(z,root/'private/FUSCC_joined_two_gene_values.tsv')
    warning={'FUSCCTNBC030','FUSCCTNBC044','FUSCCTNBC140'}
    for family,part in [('EXTERNAL_PRIMARY2',z),('EXTERNAL_EXCLUDE_WARNING2',z[~z.patient_id.isin(warning)])]:
        for h in hits:
            gene=h['gene'];row=dict(cancer='BRCA',cohort='FUSCC_TNBC',stage_id='06_EXTERNAL',run_id=root.name,analysis_version='external_two_v1',analysis_type='author_patient_linked_external_spearman',metabolite_key=h['metabolite_key'],metabolite_name=h['author_name'],gene=gene,unit='author_patient_id',n=np.nan,n_reference=np.nan,effect_type='Spearman_rho',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family=family,family_n_evaluable=1,status='NOT_EVALUABLE',reason=h['reason'],source_id='10.1038/s41422-022-00614-0;FUSCC_BRCA_2022_mrna_seq_tpm',same_aliquot_verified=False)
            if h['status']=='DONE':
                a=part[[gene,gene+'_metabolite']].replace([np.inf,-np.inf],np.nan).dropna();row['n']=len(a)
                if len(a)>=20 and a.nunique().min()>1:
                    rho,p,lo,hi,nb=infer(a[gene].values,a[gene+'_metabolite'].values,20260920+len(rows));row.update(effect=rho,p_value=p,ci_lower=lo,ci_upper=hi,n_bootstrap_valid=nb,status='DONE',reason='same_author_patient_primary_tumor;aliquot_unverified;processed_imputed_metabolite;TNBC_context')
            rows.append(row)
    r=pd.DataFrame(rows)
    for fam,ix in r.groupby('test_family').groups.items():
        ok=r.loc[ix,'p_value'].notna();r.loc[ix,'q_value']=bh(r.loc[ix,'p_value'].fillna(1));r.loc[np.array(list(ix))[~ok],'q_value']=np.nan;r.loc[ix,'family_n_evaluable']=int(ok.sum())
    save(r,root/'public/external_associations.tsv')
    cptac=project/'data/candidates/cptac_inosine_pan_cancer_v0.1/BRCA'
    files=[str(p.relative_to(cptac)) for p in cptac.rglob('*') if p.is_file()]
    (root/'public/CPTAC_available_inventory.json').write_text(json.dumps({'scope':'existing server BRCA directory only;not universal claim of CPTAC absence','file_count':len(files),'counts_by_top_directory':{k:sum(x.split('/')[0]==k for x in files) for k in sorted({x.split('/')[0] for x in files})},'metabolite_matrix_identified':False,'direct_metabolite_association_status':'NOT_EVALUABLE','reason':'existing_sources_are_protein_and_RNA_without_target_metabolite_matrix'},indent=2)+'\n')
    receipts += [dict(path=str(f),sha256=sha256(f),url='existing_server_source',method='READ') for f in [workbook,prior/'FUSCC_metabolite_RNA_patient_join_audit.tsv',prior/'FUSCC_author_RNA_sample_patient_mapping.tsv',old/'Fudan_samples.json']]
    save(pd.DataFrame(receipts),root/'public/external_source_manifest.tsv')
    (root/'public/external_validation.json').write_text(json.dumps({'status':'DONE_BOUNDED','metabolomics_patients':len(cl),'target_RNA_patients':len(wide),'author_patient_intersection':len(z),'primary_tumor_only':True,'gene_to_patient_duplicates':0,'explicit_author_mapping_concordant':True,'same_aliquot_verified':False,'strict_independent_validation':'NOT_ESTABLISHED_BY_IDENTIFIERS_ALONE','CAMP_patient_statistics_modified':False,'GPI_G6P_not_substituted':True,'metabolite_preprocessing':'author_processed_log2;original imputation not reversed;no new imputation','script_sha256':sha256(Path(__file__))},indent=2)+'\n')
    # Public plot is aggregate effect/CI, not an export of individual measurements.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,3));a=r[r.status.eq('DONE')]
    for j,t in enumerate(a.itertuples()):ax.plot([t.ci_lower,t.ci_upper],[j,j],c='#276b9c');ax.plot(t.effect,j,'o',c='#276b9c')
    ax.set_yticks(range(len(a)));ax.set_yticklabels(a.test_family);ax.axvline(0,c='gray',ls='--');ax.set_xlim(-1,1);ax.set_xlabel('Spearman rho; pointwise bootstrap CI');ax.set_title('FUSCC: GPCPD1 - GPC (author patient-linked TNBC)');fig.tight_layout();fig.savefig(root/'public/external_GPCPD1_forest.png',dpi=170)
    print(r[['gene','test_family','n','effect','p_value','q_value','status']].to_string(index=False),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--project',type=Path,required=True);a=p.parse_args();main(a.root,a.project)
