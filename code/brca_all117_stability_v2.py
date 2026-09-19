"""Source-label resampling of ALL117 cell-source summaries; no new tests."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')


def one(d,cohort,partition,seed):
    d=d.copy();d['wm']=d.mean_log1p10k*d.n_cells;d['wd']=d.detection_fraction*d.n_cells
    d=d.groupby(['donor','celltype','gene'],as_index=False)[['wm','wd','n_cells']].sum(min_count=1)
    d['mean']=d.wm/d.n_cells;d['detect']=d.wd/d.n_cells
    genes=sorted(d.gene.unique());types=sorted(d.celltype.unique());donors=sorted(d.donor.unique())
    assert len(genes)==117
    index=pd.MultiIndex.from_product([types,genes],names=['celltype','gene'])
    q=d[d.n_cells>=20]
    x=q.pivot(index='donor',columns=['celltype','gene'],values='mean').reindex(index=donors,columns=index).values
    det=q.pivot(index='donor',columns=['celltype','gene'],values='detect').reindex(index=donors,columns=index).values
    mask=np.isfinite(x);den=mask.sum(0);s=np.nansum(x,0)
    mean=np.divide(s,den,out=np.full_like(s,np.nan),where=den>=3)
    md=np.divide(np.nansum(det,0),den,out=np.full_like(s,np.nan),where=den>=3)
    rng=np.random.default_rng(seed);weights=rng.multinomial(len(donors),np.full(len(donors),1/len(donors)),size=1000)
    bden=weights@mask.astype(float)
    boot=np.divide(weights@np.nan_to_num(x),bden,out=np.full(bden.shape,np.nan),where=bden>=3)
    boot[:,den<3]=np.nan
    lden=den[None,:]-mask
    loo=np.divide(s[None,:]-np.nan_to_num(x),lden,out=np.full(x.shape,np.nan),where=lden>=3)
    means=mean.reshape(len(types),117);detect=md.reshape(len(types),117)
    boots=boot.reshape(1000,len(types),117);loos=loo.reshape(len(donors),len(types),117)
    rows=[];dist=[]
    for j,gene in enumerate(genes):
        available=np.flatnonzero(np.isfinite(means[:,j]))
        usable=len(available)>=2 and np.nanmax(detect[:,j])>=.01
        if len(available):
            order=available[np.argsort(-means[available,j],kind='stable')];top=order[0]
            runner=order[1] if len(order)>1 else None
        else:top=runner=None
        row={'cancer':'BRCA','cohort':cohort,'partition':partition,'gene':gene,'unit':'source_donor_label',
             'n_source_labels':len(donors),'n_evaluable_lineages':len(available),
             'status':'DONE' if usable else 'NOT_EVALUABLE',
             'reason':'' if usable else 'low_detection_or_fewer_than_two_covered_lineages',
             'top_lineage':types[top] if top is not None else 'NA',
             'runner_lineage':types[runner] if runner is not None else 'NA',
             'max_detection':float(np.nanmax(detect[:,j])) if len(available) else np.nan,
             'bootstrap_top_frequency':np.nan,'loo_top_frequency':np.nan,
             'paired_source_labels':0,'paired_mean_delta':np.nan,'paired_median_delta':np.nan,
             'paired_positive_fraction':np.nan,'paired_resampling_lower':np.nan,'paired_resampling_upper':np.nan,
             'p_value':np.nan,'q_value':np.nan,'inference':'descriptive;postselection contrast;not population significance'}
        if usable:
            for name,arr in [('bootstrap',boots[:,:,j]),('loo',loos[:,:,j])]:
                good=np.isfinite(arr).sum(1)>=2
                ar=arr[good];mx=np.nanmax(ar,axis=1)
                ties=np.isclose(ar,mx[:,None],rtol=1e-10,atol=1e-12)
                row[name+'_top_frequency']=float(np.mean(ties[:,top]/ties.sum(1))) if len(ar) else np.nan
            ti=top*117+j;ri=runner*117+j
            pairs=np.isfinite(x[:,ti])&np.isfinite(x[:,ri]);delta=x[pairs,ti]-x[pairs,ri]
            row['paired_source_labels']=int(pairs.sum())
            if len(delta)>=3:
                v=weights[:,pairs];n=v.sum(1);bs=(v@delta)[n>=3]/n[n>=3]
                row.update(paired_mean_delta=float(delta.mean()),paired_median_delta=float(np.median(delta)),
                           paired_positive_fraction=float((delta>0).mean()),
                           paired_resampling_lower=float(np.quantile(bs,.025)),paired_resampling_upper=float(np.quantile(bs,.975)))
        rows.append(row)
        for k,t in enumerate(types):
            ix=k*117+j;v=x[:,ix];v=v[np.isfinite(v)];b=boot[:,ix];b=b[np.isfinite(b)]
            good=len(v)>=3
            dist.append({'cancer':'BRCA','cohort':cohort,'partition':partition,'gene':gene,'celltype':t,'n':len(v),
                'effect_type':'mean_log1p_counts_per10k','effect':float(v.mean()) if good else np.nan,
                'median':float(np.median(v)) if good else np.nan,'q25':float(np.quantile(v,.25)) if good else np.nan,
                'q75':float(np.quantile(v,.75)) if good else np.nan,
                'resampling_lower':float(np.quantile(b,.025)) if good and len(b) else np.nan,
                'resampling_upper':float(np.quantile(b,.975)) if good and len(b) else np.nan,
                'n_bootstrap_evaluable':len(b),'mean_detection':md[ix],
                'status':'DONE' if good else 'NOT_EVALUABLE','p_value':np.nan,'q_value':np.nan})
    return pd.DataFrame(rows),pd.DataFrame(dist)


def main(root,previous):
    out=root/'public';allrows=[];alldist=[];manifest=[]
    for cohort in ['Wu2021','Pal2021_reprocessed','Reed2024']:
        path=previous/'private'/(cohort+'_donor_profiles.tsv');d=pd.read_csv(path,sep='\t')
        manifest.append({'kind':'private_source_summary','path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
        partitions=[('ALL',d)]
        if cohort=='Wu2021':
            partitions += [('Naive',d[d.treatment=='Naïve'])]
            partitions += [('Naive_'+str(k),v) for k,v in d[d.treatment=='Naïve'].groupby('stratum')]
        for part,block in partitions:
            seed=int(hashlib.sha256((cohort+part+'20260919').encode()).hexdigest()[:8],16)
            r,z=one(block,cohort,part,seed);allrows.append(r);alldist.append(z)
            print(cohort,part,'genes',len(r),'interpretable',r.status.eq('DONE').sum(),flush=True)
    r=pd.concat(allrows,ignore_index=True);z=pd.concat(alldist,ignore_index=True)
    save(r,out/'all117_source_stability.tsv');save(z,out/'all117_source_distributions.tsv')
    save(pd.DataFrame(manifest),out/'sc_input_manifest.tsv')
    cohorts=['Wu2021','Pal2021_reprocessed','Reed2024'];genes=sorted(r.gene.unique())
    a=r[r.partition=='ALL'].pivot(index='gene',columns='cohort',values='bootstrap_top_frequency').reindex(index=genes,columns=cohorts)
    cmap=plt.get_cmap('YlGnBu').copy();cmap.set_bad('#dddddd')
    fig,ax=plt.subplots(figsize=(7,24));im=ax.imshow(np.ma.masked_invalid(a.values),vmin=0,vmax=1,aspect='auto',cmap=cmap)
    ax.set_yticks(range(117));ax.set_yticklabels(genes,fontsize=7);ax.set_xticks(range(3));ax.set_xticklabels(cohorts,rotation=25,ha='right')
    ax.set_title('All117: source-label bootstrap top-lineage frequency\nGray = low signal / insufficient coverage')
    fig.colorbar(im,ax=ax,label='Descriptive stability, not a P value',shrink=.2);fig.tight_layout()
    fig.savefig(out/'all117_bootstrap_source_stability.png',dpi=160);plt.close(fig)
    vals={}
    for c,g in r[r.partition=='ALL'].groupby('cohort'):
        v=g[g.status=='DONE'];vals[c]={'retained':len(g),'interpretable':len(v),
            'bootstrap_ge_080':int(v.bootstrap_top_frequency.ge(.8).sum()),
            'loo_ge_090':int(v.loo_top_frequency.ge(.9).sum())}
    report={'status':'DONE','scope':'ALL117 source-label descriptive resampling','cohorts':vals,
        'rows':len(r),'distribution_rows':len(z),'new_hypothesis_tests':0,'bootstraps':1000,
        'paired_contrasts_postselected':True,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'sc_stability_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
    a=p.parse_args();main(a.root,a.previous)
