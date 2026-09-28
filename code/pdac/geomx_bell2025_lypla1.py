"""Focused GeoMx LYPLA1 reanalysis; execute in a fresh server-only run directory."""
import argparse, hashlib, itertools, json, platform, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import scipy

MD5 = {
 'GeoMx_raw_counts.csv':'01fc5162462b18b6ebae8eddfab56bf4',
 'ProbeQC_merged_batches.csv':'b710831172a61d809e4778a0b4974793',
 'fully_batch_corrected_vsd.csv':'211a3d7d3023e61b07a99d9c6ab30a27',
 'metadata_with_VI_subtypes_and_NGS.csv':'fe3b1494285185ed2db34475e65f2b5e',
 'initial_metadata.csv':'68b99d3059b40f2dc222428098d38ff7',
 'README.md':'14676fb545c51435309becf46d6ee47d'}

def bh(values):
    p=np.asarray(values); ix=np.argsort(p); out=np.empty(len(p))
    out[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]); return out

def paired(a,b,rng):
    d=np.asarray(a)-np.asarray(b); n=len(d)
    null=np.asarray(list(itertools.product([-1,1],repeat=n)))@d/n
    p=float(np.mean(np.abs(null)>=abs(d.mean())-1e-12))
    boot=rng.choice(d,size=(20000,n),replace=True).mean(axis=1)
    return dict(effect=float(d.mean()),ci_lower=float(np.quantile(boot,.025)),ci_upper=float(np.quantile(boot,.975)),p_value=p,
        positive_patients=int((d>0).sum()),negative_patients=int((d<0).sum()),zero_patients=int((d==0).sum()),
        cancer_patient_mean=float(np.mean(a)),reference_patient_mean=float(np.mean(b)),
        loo_min=float(((d.sum()-d)/(n-1)).min()),loo_max=float(((d.sum()-d)/(n-1)).max()))

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--run',required=True);p.add_argument('--commit',required=True);args=p.parse_args()
    root=Path(args.run);public=root/'public';private=root/'private';public.mkdir(exist_ok=True);private.mkdir(exist_ok=True)
    assert (root/'.running').exists(), 'Fresh exclusive run lock required'
    spec=json.loads((root/'geomx_bell2025_lypla1_spec.json').read_text()); shutil.copy(root/'geomx_bell2025_lypla1_spec.json',public/'analysis_spec.json')
    data=Path(args.data);manifest=[]
    for name,expected in MD5.items():
        content=(data/name).read_bytes();assert hashlib.md5(content).hexdigest()==expected,name
        manifest.append(dict(file=name,source_url='https://zenodo.org/api/records/16732565/files/'+name+'/content',bytes=len(content),md5=expected,sha256=hashlib.sha256(content).hexdigest()))
    pd.DataFrame(manifest).to_csv(public/'source_manifest.tsv',sep='\t',index=False)
    m=pd.read_csv(data/'metadata_with_VI_subtypes_and_NGS.csv',index_col=0).set_index('Name');x=pd.read_csv(data/'ProbeQC_merged_batches.csv',index_col=0)
    raw=pd.read_csv(data/'GeoMx_raw_counts.csv',index_col=0);v=pd.read_csv(data/'fully_batch_corrected_vsd.csv',index_col=0)
    for f in (x,raw,v):
        assert f.index.is_unique and f.columns.is_unique
        assert not f.index.isna().any() and np.isfinite(f.to_numpy()).all()
    assert m.index.is_unique and not m[['Patient.Alias','Type','Batch']].isna().any().any()
    assert len(m)==95 and m['Patient.Alias'].nunique()==8
    assert m.Type.value_counts().to_dict()=={'PDAC':38,'VI':35,'ND':14,'PNI':8}
    assert set(x.columns)-set(m.index)=={'Stroma_063_003'} and set(m.index)<=set(x.columns)
    assert set(raw.columns)==set(x.columns) and set(v.columns)==set(m.index)
    assert m.groupby('Patient.Alias').Batch.nunique().max()==1, 'Patient spans batches: review normalization'
    x=x.loc[:,m.index];v=v.loc[:,m.index]
    neg=x.loc['NegProbe-WTX'].to_numpy();assert (neg>0).all()
    loq=float(np.exp(np.log(neg).mean()+2*np.log(neg).std(ddof=1)))
    keep=(x.gt(loq).sum(axis=1)>10)&(~x.index.str.startswith('NegProbe'))
    g=x.loc[keep];gene=spec['gene'];present=gene in x.index;retained=gene in g.index
    qc=dict(gene=gene,probe_qc_genes=len(x),retained_endogenous_genes=len(g),author_vst_genes=len(v),
        global_loq=loq,gene_present=present,gene_retained=retained,author_vst_gene_present=gene in v.index,
        all_source_md5_match=True,roi_keys_match=True,patients=8,regions=95,code_commit=args.commit,
        python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,
        unverified=['Pure stromal enrichment','Single-cell localization','Cohort independence from prior datasets','Original raw-probe QC replication'])
    if not present or not retained:
        qc['status']='NOT_EVALUABLE';(public/'validation.json').write_text(json.dumps(qc,indent=2));(root/'.running').unlink();return
    q3=g.quantile(.75,axis=0);assert (q3>0).all()
    norm=x.loc[gene]/q3*np.exp(np.log(q3).mean())
    roi=m[['Patient.Alias','Type','Batch']].copy();roi['raw_count']=x.loc[gene];roi['above_global_loq']=x.loc[gene]>loq
    roi['q3_log2']=np.log2(norm+1);roi['cpm_log2']=np.log2(x.loc[gene]/g.sum(axis=0)*1e6+1)
    roi['author_vst']=v.loc[gene] if gene in v.index else np.nan
    roi.to_csv(private/'roi_lypla1.tsv',sep='\t')
    group=[]
    for typ,f in roi.groupby('Type'):
        z=f.groupby('Patient.Alias').q3_log2.mean()
        group.append(dict(group=typ,regions=len(f),patients=len(z),above_global_loq_regions=int(f.above_global_loq.sum()),
            above_global_loq_percent=100*f.above_global_loq.mean(),roi_mean_q3_log2=f.q3_log2.mean(),patient_mean_q3_log2=z.mean(),patient_sd_q3_log2=z.std()))
    pd.DataFrame(group).to_csv(public/'group_summary.tsv',sep='\t',index=False)
    rng=np.random.default_rng(spec['seed']);rows=[];pair_frames=[]
    for scale in ['q3_log2','cpm_log2','author_vst']:
        for contrast,types in spec['contrasts'].items():
            a=roi[roi.Type.isin(types)].groupby('Patient.Alias')[scale].mean()
            b=roi[roi.Type=='ND'].groupby('Patient.Alias')[scale].mean()
            ab=pd.concat([a.rename('cancer'),b.rename('reference')],axis=1).dropna()
            row=dict(cancer='PDAC',cohort=spec['cohort'],stage_id='06_EXTERNAL',run_id=root.name,analysis_version=spec['analysis_version'],
                analysis_type='paired_patient_spatial_expression',metabolite_key='NA',metabolite_name='NA',gene=gene,unit='paired_patients',n=len(ab),n_reference=len(ab),
                effect_type='mean_paired_difference_'+scale,effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,
                test_family='four_contrasts_'+scale,family_n_evaluable=0,status='NOT_EVALUABLE',reason='Fewer than 3 eligible paired patients',source_id='Zenodo16732565',contrast=contrast,scale=scale)
            if len(ab)>=spec['minimum_pairs']:
                row.update(paired(ab.cancer,ab.reference,rng));row.update(status='DONE',reason='Focused gene comparison; not genome-wide FDR')
            ab['contrast']=contrast;ab['scale']=scale;pair_frames.append(ab.reset_index());rows.append(row)
    out=pd.DataFrame(rows)
    for scale in out.scale.unique():
        mask=(out.scale==scale)&out.p_value.notna();out.loc[mask,'q_value']=bh(out.loc[mask,'p_value']);out.loc[out.scale==scale,'family_n_evaluable']=int(mask.sum())
    out.to_csv(public/'results.tsv',sep='\t',index=False,na_rep='NA');pd.concat(pair_frames).to_csv(private/'patient_pairs.tsv',sep='\t',index=False)
    qc['status']='DONE';qc['n_contrasts']=len(out);qc['n_evaluable']=int(out.p_value.notna().sum())
    (public/'validation.json').write_text(json.dumps(qc,indent=2));(root/'.running').unlink();print(out[['contrast','scale','n','effect','p_value','q_value','positive_patients','negative_patients']].to_string(index=False))

if __name__=='__main__':main()
