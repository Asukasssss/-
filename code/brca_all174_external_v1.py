"""Fixed174 FUSCC external associations; individual values remain on server."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import argparse,json,hashlib,concurrent.futures
from pathlib import Path
import requests,numpy as np,pandas as pd
from scipy import stats
from brca_pycr1_context_v1 import infer,bh,save
from brca_sc117_profile_v1 import sha256

def main(root,project):
    rel=pd.read_csv(root/'source/relations_174_annotated_v4.tsv',sep='\t');assert len(rel)==174 and rel.relation_id.is_unique and rel.gene.nunique()==117
    annpath=project/'results/CAMP_first_analysis_20260910/tables/FUSCC_author_polar_annotation_594.tsv'
    ann=pd.read_csv(annpath,sep='\t');assert len(ann)==594 and ann.fuscc_peak.is_unique
    cov=[]
    for key,name in rel[['metabolite_key','metabolite']].drop_duplicates().values:
        typ,val=key.split(':');a=ann[ann.fuscc_kegg.eq(val)] if typ=='KEGG' else ann[ann.fuscc_hmdb.eq('HMDB'+val[4:].zfill(7))]
        reason='unique_author_identifier;identity_not_spectrally_reconfirmed' if len(a)==1 else ('no_exact_author_identifier' if len(a)==0 else 'multiple_features_no_outcome_based_selection')
        status='DONE' if len(a)==1 else 'NOT_EVALUABLE'
        if key=='KEGG:C00186' and len(a)==1:status='NEEDS_REVIEW';reason='L_lactate_identifier_but_author_name_DL_lactate;stereochemical_scope_unresolved'
        h=dict(metabolite_key=key,CAMP_name=name,annotation_hits=len(a),status=status,reason=reason)
        for field in ['fuscc_peak','fuscc_author_name','fuscc_hmdb','fuscc_kegg','fuscc_identification_level']:h[field]=a.iloc[0][field] if len(a)==1 else ''
        cov.append(h)
    coverage=pd.DataFrame(cov);save(coverage,root/'public/identity_coverage.tsv');rel=rel.merge(coverage,on='metabolite_key',validate='many_to_one',suffixes=('_CAMP','_coverage'))
    ready=rel[rel.status_coverage.eq('DONE')];genes=sorted(ready.gene.unique())
    old=project/'data/candidates/BRCA_1MNA_followup_v0.1';prior=project/'results/BRCA_1MNA_followup_v0.1/tables'
    prev=project/'results/collaborative/BRCA/A/20260920T050412Z_pycr1_external_v1'
    sm=pd.DataFrame(json.loads((old/'Fudan_samples.json').read_text()));assert sm.patientId.is_unique and sm.sampleId.is_unique and sm.sampleType.eq('Primary Solid Tumor').all()
    api='https://data.3steps.cn/cdataportal/api';gene_errors={};info={}
    def fetch_gene(g):
        p=root/'source'/f'gene_{g}.json'
        try:
            if not p.exists():
                x=requests.get(api+'/genes/'+g,timeout=35);x.raise_for_status();p.write_text(json.dumps(x.json()))
            j=json.loads(p.read_text());assert j['hugoGeneSymbol']==g;return g,j,None
        except Exception as e:return g,None,str(e)[:160]
    for g,j,e in concurrent.futures.ThreadPoolExecutor(max_workers=3).map(fetch_gene,genes):
        if j is not None:info[g]=j
        else:gene_errors[g]=e
    save(pd.DataFrame([dict(gene=g,status='DONE' if g in info else 'ACCESS_BLOCKED',reason=gene_errors.get(g,'')) for g in genes]),root/'public/gene_identity_coverage.tsv')
    profile='FUSCC_BRCA_2022_mrna_seq_tpm';parts=[]
    for i in range(0,len(info),15):
        gs=sorted(info)[i:i+15];p=root/'source'/f'RNA_batch_{i:03d}.json'
        if not p.exists():
            body={'entrezGeneIds':[info[g]['entrezGeneId'] for g in gs],'sampleIds':sm.sampleId.tolist()}
            x=requests.post(api+'/molecular-profiles/'+profile+'/molecular-data/fetch?projection=DETAILED',json=body,timeout=90);x.raise_for_status();p.write_text(json.dumps(x.json()))
        q=pd.DataFrame(json.loads(p.read_text()));parts.append(q);print('RNA batch',i,len(q),flush=True)
    ex=pd.concat(parts,ignore_index=True);ex['gene']=ex.gene.map(lambda x:x['hugoGeneSymbol']);assert not ex.duplicated(['patientId','gene']).any();assert (ex.value>=0).all()
    ch=ex.merge(sm[['sampleId','patientId']],on='sampleId',validate='many_to_one',suffixes=('_a','_b'));assert ch.patientId_a.eq(ch.patientId_b).all()
    wide=ex.pivot(index='patientId',columns='gene',values='value').reset_index().rename(columns={'patientId':'patient_id'})
    cl=pd.read_csv(prior/'FUSCC_1MNA_tumor_clinical_values.tsv',sep='\t');rm=pd.read_csv(prior/'FUSCC_author_RNA_sample_patient_mapping.tsv',sep='\t');rm=rm[rm.tissue_type.eq('tumor')];jo=pd.read_csv(prior/'FUSCC_metabolite_RNA_patient_join_audit.tsv',sep='\t')
    assert cl.patient_id.is_unique and cl.metabolite_sample_id.is_unique and rm.patient_id.is_unique and jo.patient_id.is_unique
    z=cl[['patient_id','metabolite_sample_id']].merge(wide,on='patient_id',validate='one_to_one').merge(rm[['patient_id','sample_id']],on='patient_id',validate='one_to_one').merge(jo[['patient_id','author_RNA_sample_id','metabolite_sample_id','join_status']],on=['patient_id','metabolite_sample_id'],validate='one_to_one')
    assert z.sample_id.eq(z.author_RNA_sample_id).all() and z.join_status.eq('PASS_EXPLICIT_AUTHOR_PATIENT_MAPPING').all()
    wp=project/'data/candidates/brca_four_axes_20260908/FUSCC_41422_2022_614_MOESM12_ESM.xlsx';mat=pd.read_excel(wp,sheet_name='Table S1',header=1);mat.columns=[str(c).strip() for c in mat.columns];assert mat.Peak.is_unique;mat=mat.set_index('Peak')
    for peak in ready.fuscc_peak.unique():z[peak]=pd.to_numeric(mat.loc[peak,z.metabolite_sample_id],errors='coerce').values
    save(z,root/'private/author_patient_joined_values.tsv')
    # Compare and analyze the persisted TSV representation used by the previous run.
    # Excel binary floats and default CSV parsing can differ below machine rounding precision.
    z=pd.read_csv(root/'private/author_patient_joined_values.tsv',sep='\t')
    oldz=pd.read_csv(prev/'private/FUSCC_joined_two_gene_values.tsv',sep='\t').set_index('patient_id');zz=z.set_index('patient_id')
    assert set(oldz.index)==set(zz.index)
    assert np.array_equal(oldz.GPCPD1.values,zz.loc[oldz.index,'GPCPD1'].values) and np.array_equal(oldz.GPCPD1_metabolite.values,zz.loc[oldz.index,'M258T238_POS'].values)
    oldres=pd.read_csv(prev/'public/external_associations.tsv',sep='\t');rows=[];warnings={'FUSCCTNBC030','FUSCCTNBC044','FUSCCTNBC140'}
    for family,part,oldfam in [('FUSCC_PRIMARY174',z,'EXTERNAL_PRIMARY2'),('FUSCC_EXCLUDE_WARNING174',z[~z.patient_id.isin(warnings)],'EXTERNAL_EXCLUDE_WARNING2')]:
        for i,t in enumerate(rel.itertuples()):
            row=dict(cancer='BRCA',cohort='FUSCC_TNBC',stage_id='06_EXTERNAL',run_id=root.name,analysis_version='external174_v1',analysis_type='author_patient_linked_spearman',metabolite_key=t.metabolite_key,metabolite_name=t.metabolite,gene=t.gene,unit='author_patient_id',n=np.nan,n_reference=np.nan,effect_type='Spearman_rho',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family=family,family_n_evaluable=np.nan,status=t.status_coverage,reason=t.reason,source_id='FUSCC_BRCA_2022;10.1038/s41422-022-00614-0',relation_id=t.relation_id,external_peak=t.fuscc_peak,external_identification_level=t.fuscc_identification_level,CAMP_ER_rho=t.er_available_partial_rank_rho,CAMP_ER_q=t.er_available_q,old_external_q=np.nan,statistics_reused=False)
            if row['status']=='DONE':
                if t.gene not in part:row.update(status='NOT_EVALUABLE',reason='RNA_not_available')
                else:
                    a=part[[t.gene,t.fuscc_peak]].replace([np.inf,-np.inf],np.nan).dropna();row['n']=len(a)
                    if len(a)<20 or a.nunique().min()<2:row.update(status='NOT_EVALUABLE',reason='less_than20_or_constant')
                    elif t.gene=='GPCPD1' and t.metabolite_key=='KEGG:C00670':
                        o=oldres[(oldres.gene=='GPCPD1')&(oldres.test_family==oldfam)].iloc[0]
                        for col in ['effect','ci_lower','ci_upper','p_value']:row[col]=o[col]
                        row.update(old_external_q=o.q_value,statistics_reused=True)
                    else:
                        seed=int(hashlib.sha256((family+t.relation_id).encode()).hexdigest()[:8],16)
                        rho,p,lo,hi,nb=infer(a.iloc[:,0].to_numpy(),a.iloc[:,1].to_numpy(),seed);row.update(effect=rho,p_value=p,ci_lower=lo,ci_upper=hi)
                    if row['status']=='DONE':assert np.isclose(row['effect'],stats.spearmanr(a.iloc[:,0],a.iloc[:,1]).statistic)
            rows.append(row)
        save(pd.DataFrame(rows),root/'private/checkpoint_results.tsv');print('family complete',family,flush=True)
    res=pd.DataFrame(rows)
    for fam,ix in res.groupby('test_family').groups.items():
        ok=res.loc[ix,'p_value'].notna();res.loc[ix,'q_value']=bh(res.loc[ix,'p_value'].fillna(1));res.loc[np.array(list(ix))[~ok],'q_value']=np.nan;res.loc[ix,'family_n_evaluable']=int(ok.sum())
    res['direction_vs_CAMP_ER']=np.where(res.effect.notna()&res.CAMP_ER_rho.notna(),np.where(res.effect*res.CAMP_ER_rho>0,'SAME_SIGN','OPPOSITE_OR_ZERO'),'NOT_EVALUABLE')
    save(res,root/'public/external_all174_associations.tsv')
    val={'status':'DONE','relations':174,'genes':117,'unique_metabolite_keys':rel.metabolite_key.nunique(),'n_author_patient_intersection':len(z),'identity_counts':coverage.status.value_counts().to_dict(),'families':res.groupby('test_family').apply(lambda a:{'evaluable':int(a.p_value.notna().sum()),'q_lt005':int((a.q_value<.05).sum())}).to_dict(),'GPCPD1_inputs_exactly_same_and_statistics_reused':True,'BH_planned_family_n':174,'missing_p_internal_placeholder':1,'missing_public_p_q':'NA','same_aliquot_verified':False,'raw_missing_mask_verified':False,'CAMP_adjusted_vs_external_unadjusted_not_same_estimand':True,'patient_values_exported':False}
    (root/'public/external_validation.json').write_text(json.dumps(val,indent=2)+'\n')
    save(pd.DataFrame([{'path':str(p),'sha256':sha256(p)} for p in [annpath,wp,root/'source/relations_174_annotated_v4.tsv',prior/'FUSCC_author_RNA_sample_patient_mapping.tsv',prior/'FUSCC_metabolite_RNA_patient_join_audit.tsv',old/'Fudan_samples.json']]+[{'path':str(p),'sha256':sha256(p)} for p in sorted((root/'source').glob('*.json'))]),root/'public/external_source_manifest.tsv')
    print(json.dumps(val),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--project',type=Path,required=True);a=p.parse_args();main(a.root,a.project)
