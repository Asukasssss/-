"""Finalize public aggregate-only BRCA sensitivity outputs; no source matrices."""
from pathlib import Path
import argparse,hashlib,json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--run',required=True);a=ap.parse_args()
    out=a.repo/'results/BRCA/04_ROBUSTNESS'/a.run
    d=pd.read_csv(out/'association_174_sensitivities.tsv',sep='\t');e=pd.read_csv(out/'paired_RNA117.tsv',sep='\t')
    oldpath=a.repo/'results/BRCA/07_INTEGRATION/20260921T093756Z_candidate_discussion_v1/all117_comparison_discussion_appended.tsv'
    old=pd.read_csv(oldpath,sep='\t',dtype=str,keep_default_na=False);new=old.copy();assert len(old)==117 and old.gene.is_unique
    p=e.set_index('gene')
    for col in ['status','n','effect','ci_lower','ci_upper','p_value','q_value']:
        new['identity_v1_paired_RNA_'+col]=new.gene.map(p[col])
    for fam,g in d.groupby('test_family'):
        count=g[g.p_value.notna()].groupby('gene').size();sig=g[g.q_value.lt(.05)].groupby('gene').size()
        new['identity_v1_'+fam+'_n_evaluable']=new.gene.map(count).fillna(0).astype(int)
        new['identity_v1_'+fam+'_n_q_lt005']=new.gene.map(sig).fillna(0).astype(int)
    new['identity_v1_scope']='60 author tumor cases after conflict exclusion;45 paired cases;not external validation;old columns retained as history'
    new['identity_v1_source']=out.relative_to(a.repo).as_posix()
    save(new,out/'all117_comparison_identity_appended.tsv')
    check=pd.read_csv(out/'all117_comparison_identity_appended.tsv',sep='\t',dtype=str,keep_default_na=False)
    assert check[list(old.columns)].equals(old)
    focus=['ASNS','GLS','SLC6A8','NNMT','GPCPD1','GPI','MDH1','PNP','PYCR1']
    save(d[d.gene.isin(focus)],out/'focus_relations_same_ID.tsv');save(e[e.gene.isin(focus)],out/'focus_paired_RNA.tsv')
    rows=[]
    for fam,g in d.groupby('test_family'):
        v=g[g.p_value.notna()&g.previous_p.notna()];changed=v[v.q_threshold_changed]
        rows.append(dict(family=fam,evaluable=len(v),old_q_lt005=int(v.previous_q.lt(.05).sum()),new_q_lt005=int(v.q_value.lt(.05).sum()),
             max_abs_rho_change=float(v.effect_change.abs().max()),median_abs_rho_change=float(v.effect_change.abs().median()),
             q_threshold_crossings=len(changed),effect_sign_changes=int((np.sign(v.effect)!=np.sign(v.previous_effect)).sum()),
             reused_p=int(g.p_origin.eq('REUSED_IDENTICAL_RETAINED_OBSERVATIONS').sum())))
    summary=pd.DataFrame(rows);save(summary,out/'family_before_after.tsv')
    save(d[d.q_threshold_changed],out/'q_threshold_changes.tsv')
    fig,ax=plt.subplots(1,3,figsize=(15,4.6))
    q=d[d.test_family.eq('processed_Spearman')&d.effect.notna()]
    ax[0].scatter(q.previous_effect,q.effect,s=16,c='#2878A0',alpha=.7);ax[0].plot([-.65,.65],[-.65,.65],ls='--',c='grey',lw=1)
    ax[0].set(xlabel='Original Spearman rho',ylabel='After excluding conflict',title='All evaluable relations',xlim=(-.65,.65),ylim=(-.65,.65))
    labels=['processed_Spearman','available_Spearman','processed_ER','available_ER','available_ER_PC12','available_ER_MonoEndoFib']
    ax[1].boxplot([d.loc[d.test_family.eq(k),'effect_change'].dropna() for k in labels],vert=False,tick_labels=['Processed','Available','ER / processed','ER / available','ER + PC1/2','ER + Mono/Endo/Fib'],showfliers=True)
    ax[1].axvline(0,c='grey',lw=1);ax[1].set(xlabel='New rho - original rho',title='Effect changes by fixed model')
    v=e[e.p_value.notna()];color=np.where(v.q_value.lt(.05),'#2878A0','#AAAAAA');ax[2].scatter(v.effect,-np.log10(v.q_value.clip(lower=1e-300)),c=color,s=16)
    ax[2].axvline(0,c='grey',lw=1);ax[2].axhline(-np.log10(.05),c='grey',ls='--',lw=1)
    ax[2].set(xlabel='Paired tumor - normal mean\n(author expression scale; not fold change)',ylabel='-log10 BH q',title='117-gene paired RNA analysis')
    for b in ax:b.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(out/'analysis_overview.png',dpi=180);plt.close(fig)
    val=json.loads((out/'validation.json').read_text());val['historical_columns_preserved']=len(old.columns);val['appended_master_columns']=len(new.columns);val['master_rows']=len(new)
    val['finalizer_input_master_sha256']=hashlib.sha256(oldpath.read_bytes()).hexdigest();(out/'validation.json').write_bytes((json.dumps(val,indent=2)+'\n').encode())
    n=int(e.q_value.lt(.05).sum());oldn=int(e.previous_unpaired_q.lt(.05).sum());newonly=e[e.q_value.lt(.05)&~e.previous_unpaired_q.lt(.05)].gene.tolist();oldonly=e[~e.q_value.lt(.05)&e.previous_unpaired_q.lt(.05)].gene.tolist()
    lines=['# BRCA 冲突标本排除与45对RNA补充分析','','## 本轮问题','排除来源标签冲突标本后，原174条关联及调整结果是否变化？在作者确认的45对组织中，117基因RNA的肿瘤—正常差异是什么？','','## 输入与范围','沿用身份核查的显式病例表，保留60个肿瘤病例和45对肿瘤—正常。冲突标本排除，不重标为正常。输入矩阵、完整映射和病例测量仅在server165。原代谢物效应、P/q及所有历史候选列保持不变。',
      '6个174项检验家族分别计算并各自BH；主关联与作者数据可用性子集各一套，ER两套，ER+PC12与ER+单核/内皮/成纤维各一套。作者数据可用性不等于原始检测掩码。PCA沿用9类标记集合，在60例上重新标准化拟合。相同保留观察且协变量未变的结果复用原P；新q依照整个固定家族计算。',
      '配对RNA：未另加log变换，计算每对肿瘤减正常的表达差，双侧配对t检验、t分布95%区间、4000次整对重采样均值区间；BH覆盖全部可评估候选。检验均值差依赖病例独立及均值推断近似；不把作者处理尺度差值称作倍数。',
      '', '## 实际结果','|模型|可评估|原q<0.05|新q<0.05|最大绝对rho变化|','|---|---:|---:|---:|---:|']
    for _,r in summary.iterrows():lines.append(f'|{r.family}|{r.evaluable}|{r.old_q_lt005}|{r.new_q_lt005}|{r.max_abs_rho_change:.6f}|')
    lines += ['重点解释：GPCPD1—GPC与GPI—G6P在组成主模型中的效应分别约-0.466、+0.447，与原来-0.460、+0.449非常接近，但本轮q均为0.0513。次模型分别约-0.529/q0.0171、+0.459/q0.0342。不能把主模型未过线写成关联消失。9999次置换的随机误差、样本删减及PCA重拟合均需考虑；不为过线增加置换或更换检验族。',
      f'配对RNA：117个保留，{e.p_value.notna().sum()}个可评估，{n}个本轮q<0.05；旧非配对版本为{oldn}个。两版本同时改变配对设计、样本子集、效应指标和检验，显著数差不能单独归因于某一因素。',
      '配对显著而原非配对未显著：'+('; '.join(newonly) or '无')+'。','原非配对显著而本次配对未显著：'+('; '.join(oldonly) or '无')+'。',
      '', '## 新手解释','关联检查回答：去掉身份有争议的标本后，基因与代谢物是否仍一起变化。配对RNA检查回答：同一个患者的肿瘤与正常相比，基因表达平均变化多少。这两个问题分别解释，RNA差异不作为功能研究的新门槛。',
      '', '## 限制/反证','本轮为同队列的来源纠错敏感性，不是独立复现。组织标签冲突仍未解决；排除只用于减少歧义。作者病例表确认不等于基因型鉴定。组成变量仍为修改版表达代理评分，并非细胞比例。重采样区间为逐项区间，不是同时覆盖区间；调整模型的500次区间精度有限。检验族间q不直接排序。原冻结CAMP代谢差异的标签风险继续登记，未在本轮重算。',
      '', '## 当前决定','全117候选保留。新结果追加到全候选表；旧患者列是历史版本，不代表它们已完成本次身份敏感性更新。调整结果按完整关系ID与同模型比较。缺测不当作阴性，阈值跨越不自动改变功能候选价值。',
      '', '## 下一步','以本轮身份明确的敏感性和配对RNA列解释CAMP证据，结合已有外部与功能结果继续判断；不因q阈值变化反复调模型。若需正式替换原发现阶段统计，另立有明确范围的修正版。',
      '', '## 复现命令',f'`python3 brca_camp_pair_sensitivity_v1.py --out {ROOTSTR}/{a.run}`',f'`python code/brca_camp_pair_finalize_v1.py --repo . --run {a.run}`','首次计算需新目录和本任务.running锁；需将两个已版本化依赖脚本camp_per_cancer_repro.py与brca_all174_composition_v2.py放在同目录。参数预先写入analysis_spec.json，哈希见source_manifest.tsv；Python结果复核OLS残差、BH及配对t公式。']
    (out/'README_CN.md').write_bytes(('\n'.join(lines)+'\n').encode())
    print(summary.to_string(index=False));print('PAIRED',n,'newonly',newonly,'oldonly',oldonly)

ROOTSTR='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A'
if __name__=='__main__':main()
