"""All117 integration and aggregate plots; preserve every prior column verbatim."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')


def prefix(d,root,analysis_type,source):
    cols='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
    defaults={'cancer':'BRCA','stage_id':'04_ROBUSTNESS','run_id':root.name,'analysis_version':'all117_robustness_v2',
              'analysis_type':analysis_type,'source_id':source,'test_family':'NA_descriptive_no_tests'}
    for c in cols:
        if c not in d:d[c]=defaults.get(c,np.nan)
    return d[cols+[c for c in d if c not in cols]]


def main(root,previous):
    p=root/'public'
    stability=pd.read_csv(p/'all117_source_stability.tsv',sep='\t')
    distributions=pd.read_csv(p/'all117_source_distributions.tsv',sep='\t')
    relations=pd.read_csv(p/'all174_composition_sensitivity.tsv',sep='\t')
    base=stability[stability.partition=='ALL']
    assert len(base)==351 and base.groupby('cohort').gene.nunique().eq(117).all()
    oldpath=previous/'public/candidate_comparison_117_sc_appended.tsv'
    legacy=pd.read_csv(oldpath,sep='\t',dtype=str,keep_default_na=False)
    genes=legacy.gene.tolist();assert len(genes)==len(set(genes))==117
    assert relations.groupby('test_family').size().eq(174).all()
    rows=[]
    for gene in genes:
        row={'gene':gene}
        for cohort in ['Wu2021','Pal2021_reprocessed','Reed2024']:
            r=base[(base.gene==gene)&(base.cohort==cohort)].iloc[0]
            for col in ['status','top_lineage','runner_lineage','bootstrap_top_frequency','loo_top_frequency',
                        'paired_source_labels','paired_mean_delta','paired_positive_fraction']:
                row['rob_'+cohort+'_'+col]=r[col]
        for model in ['ER_PC12','ER_MonoEndoFib']:
            a=relations[(relations.gene==gene)&(relations.test_family==model)]
            valid=a[a.status=='DONE'].sort_values(['q_value','relation_id'])
            row['rob_'+model+'_n_evaluable']=len(valid)
            row['rob_'+model+'_n_q_lt005']=int(valid.q_value.lt(.05).sum())
            if len(valid):
                for col in ['relation_id','effect','ci_lower','ci_upper','p_value','q_value','ER_v3_rho','ER_v3_q','abs_effect_change','n']:
                    row['rob_'+model+'_best_'+col]=valid.iloc[0][col]
        old=legacy[legacy.gene==gene].iloc[0]
        w=row['rob_Wu2021_bootstrap_top_frequency'];t=row['rob_Pal2021_reprocessed_bootstrap_top_frequency']
        stable=(np.isfinite(w) and np.isfinite(t) and w>=.8 and t>=.8 and
                row['rob_Wu2021_top_lineage']==row['rob_Pal2021_reprocessed_top_lineage'])
        row['rob_two_tumor_source_stable_descriptive']=stable
        if row['rob_ER_PC12_n_q_lt005']>0:
            decision='患者关联在主代理模型中保留；结合适用功能背景继续'
        elif row['rob_ER_PC12_n_evaluable']==0:
            decision='患者项不可评估；保留功能和细胞来源证据'
        elif old['receiver_v4_genetic_report_count_basis']=='IMPORTED_GENETIC_BRCA_CLASS':
            decision='已有遗传干预报告；按原研究模型推进，不因本轮q不显著淘汰'
        elif stable:
            decision='细胞来源较稳定；功能或患者支持仍需补充'
        else:decision='保留；按检出覆盖、背景差异及功能缺口安排下一步'
        row['rob_下一步方向']=decision
        row['rob_证据边界']='组成相关表达代理，不是真实细胞比例；标本独立性未核验；不证明因果或介导'
        rows.append(row)
    summary=pd.DataFrame(rows)
    assert not ((set(summary)-{'gene'})&set(legacy))
    joined=legacy.merge(summary,on='gene',how='left',sort=False,validate='one_to_one')
    assert joined[legacy.columns].equals(legacy)
    save(summary,p/'all117_next_actions_CN.tsv');save(joined,p/'all117_comparison_robustness_appended.tsv')
    # All117 same-donor top-v-runner differences. Selection is descriptive, not a test.
    fig,axes=plt.subplots(1,3,figsize=(15,25),sharey=True)
    for ax,cohort in zip(axes,['Wu2021','Pal2021_reprocessed','Reed2024']):
        d=base[base.cohort==cohort].set_index('gene').reindex(genes)
        for j,(_,r) in enumerate(d.iterrows()):
            if np.isfinite(r.paired_mean_delta):
                ax.plot([r.paired_resampling_lower,r.paired_resampling_upper],[j,j],color='#789ab2',lw=.8)
                ax.plot(r.paired_mean_delta,j,'o',color='#145a7a',ms=2.8)
        ax.axvline(0,color='#888888',lw=.8,ls='--');ax.set_title(cohort)
        ax.set_yticks(range(117));ax.set_yticklabels(genes,fontsize=7);ax.invert_yaxis()
        ax.set_xlabel('Within-source mean log1p10k difference')
    fig.suptitle('All117: paired top-v-runner descriptive differences\nPostselected comparison; resampling ranges are not significance tests',y=.998)
    fig.tight_layout(rect=[0,0,1,.985]);fig.savefig(p/'all117_paired_source_distributions.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,5),sharex=True,sharey=True)
    for ax,(model,d) in zip(axes,relations.groupby('test_family',sort=True)):
        d=d[d.status=='DONE'];sig=d.q_value.lt(.05)
        ax.scatter(d.ER_v3_rho,d.effect,c=np.where(sig,'#c34a36','#417d9c'),s=20,alpha=.75)
        ax.plot([-1,1],[-1,1],ls='--',color='gray',lw=1)
        ax.axhline(0,c='gray',lw=.5);ax.axvline(0,c='gray',lw=.5)
        ax.set(xlim=(-.8,.8),ylim=(-.8,.8),xlabel='Same-subset ER-adjusted rho',ylabel='Composition-proxy-adjusted rho',title=model)
    fig.suptitle('All evaluable relations; red = new within-model q<0.05\nProxy scores are not measured cell fractions',fontsize=11)
    fig.tight_layout();fig.savefig(p/'all174_composition_before_after.png',dpi=160);plt.close(fig)
    # Use the common cross-cancer statistical prefix in final public numeric tables.
    stability['n']=stability.n_source_labels;stability['effect_type']='bootstrap_top_lineage_frequency'
    stability['effect']=stability.bootstrap_top_frequency
    save(prefix(stability,root,'source_label_rank_stability','previous_sc117_private_aggregates'),p/'all117_source_stability.tsv')
    distributions['unit']='source_donor_label';distributions['ci_lower']=np.nan;distributions['ci_upper']=np.nan
    save(prefix(distributions,root,'source_label_expression_distribution','previous_sc117_private_aggregates'),p/'all117_source_distributions.tsv')
    save(prefix(relations,root,'partial_rank_composition_proxy_sensitivity','CAMP_BRCA1'),p/'all174_composition_sensitivity.tsv')
    vals={'status':'DONE','all117_retained':True,'previous_columns_preserved':len(legacy.columns),
          'previous_values_exactly_unchanged':True,'n_relation_rows':len(relations),
          'two_tumor_source_stable_descriptive':int(summary.rob_two_tumor_source_stable_descriptive.sum()),
          'genes_primary_q_lt005':summary.loc[summary.rob_ER_PC12_n_q_lt005.gt(0),'gene'].tolist(),
          'genes_secondary_q_lt005':summary.loc[summary.rob_ER_MonoEndoFib_n_q_lt005.gt(0),'gene'].tolist(),
          'primary_family':'ER_PC12','secondary_family':'ER_MonoEndoFib','clinical_target_validation':False,
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (p/'integration_validation.json').write_text(json.dumps(vals,indent=2)+'\n');print(json.dumps(vals),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
    a=p.parse_args();main(a.root,a.previous)
