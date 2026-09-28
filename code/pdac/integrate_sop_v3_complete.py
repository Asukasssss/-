"""Full relation/gene integration. Distinct evidence dimensions; no composite score."""
from pathlib import Path
import json,itertools,subprocess
import numpy as np
import pandas as pd
from record_sop_v3_resume import REPO,ROOT,INTERNAL,SOURCE,DISCOVERY,INTEGRATION,MAP,read,write,write_json,checksum,stage,readme,digest,VERSION,SOP_SHA

PAT=ROOT/'03_PATIENT/20260922T112700Z_sop_v3_patient'
COHORTS=['GSE263733','GSE278688','GSE242230']
def evidence(p,q):return 'NOT_EVALUABLE' if pd.isna(p) else 'FDR_SUPPORTED' if q<.05 else 'NOMINAL_EXPLORATORY' if p<.05 else 'NOT_SUPPORTED_THIS_TEST'

def main():
    INTEGRATION.mkdir(exist_ok=True)
    rel=pd.concat([read(MAP/'direct_relations.tsv'),read(MAP/'conditional_relations.tsv')],ignore_index=True);pool=read(MAP/'gene_pool_history_union.tsv')
    m=read(DISCOVERY/'metabolite_paired.tsv').set_index('metabolite_key');ma=read(DISCOVERY/'metabolite_available_sensitivity.tsv').set_index('metabolite_key')
    a=pd.concat([read(PAT/'tumor_association.tsv'),read(PAT/'conditional_tumor_association.tsv')]).set_index('relation_id');av=pd.concat([read(PAT/'tumor_association_available.tsv'),read(PAT/'conditional_tumor_association_available.tsv')]).set_index('relation_id')
    RNA=pd.concat([read(INTERNAL/'paired_RNA.tsv'),read(INTERNAL/'paired_RNA_history_supplement.tsv')]).set_index('stable_gene_id')
    sc=read(SOURCE/'sc_source_stability.tsv');profiles=read(SOURCE/'sc_celltype_profiles.tsv');coverage=read(SOURCE/'sc_gene_coverage.tsv');cross=read(SOURCE/'sc_cross_study.tsv')
    assert len(rel)==715 and rel.relation_id.is_unique and len(pool)==687 and pool.stable_gene_id.is_unique and len(RNA)==687 and RNA.index.is_unique
    old=read(ROOT/'02_MAPPING/20260920T111500Z_direct_mapping_v1/direct_relations_v1.tsv')
    source_summary=[]
    for _,g in pool.iterrows():
        tt=sc[sc.stable_gene_id==g.stable_gene_id].set_index('cohort');cc=coverage[coverage.stable_gene_id==g.stable_gene_id].set_index('cohort');assert set(tt.index)==set(COHORTS)
        valid=tt.status.eq('DONE')&tt.tie_status.eq('UNIQUE');same=bool(valid.all() and tt.top_celltype.nunique()==1);stable=bool(same and tt.bootstrap_top_frequency.ge(.8).all())
        cat='ALL3_SAME_TOP_AND_BOOTSTRAP_GE080' if stable else 'ALL3_SAME_TOP' if same else 'STUDY_DEPENDENT_OR_NOT_EVALUABLE'
        rr={'cancer':'PDAC','gene':g.gene,'stable_gene_id':g.stable_gene_id,'all3_same_unique_top':same,'all3_top_bootstrap_ge080':stable,'cell_background':cat,'source_type':'Expression-source descriptive only;not unique action cell'}
        for c in COHORTS:
            t=tt.loc[c];ct=t.top_celltype
            for k in ['top_celltype','runner_celltype','top_gap','bootstrap_top_frequency','bootstrap_valid','max_detection','n_evaluable_celltypes','status']:rr[c+'_'+k]=t[k]
            rr[c+'_matrix_coverage_status']=cc.loc[c,'status']
            p=profiles[(profiles.stable_gene_id==g.stable_gene_id)&(profiles.cohort==c)&(profiles.celltype==ct)]
            rr[c+'_top_detection_fraction']=p.mean_detection_fraction.iloc[0] if len(p)==1 else np.nan
            rr[c+'_top_source_labels']=p.n.iloc[0] if len(p)==1 else np.nan
        source_summary.append(rr)
    sg=pd.DataFrame(source_summary).set_index('stable_gene_id');assert len(sg)==687
    write(INTEGRATION/'gene_source_three_study_summary.tsv',sg.reset_index())
    rows=[]
    for _,r in rel.iterrows():
        mm=m.loc[r.metabolite_key];ms=ma.loc[r.metabolite_key];aa=a.loc[r.relation_id];vv=av.loc[r.relation_id];rn=RNA.loc[r.stable_gene_id];cell=sg.loc[r.stable_gene_id]
        assert rn.gene==r.gene and cell.gene==r.gene and aa.gene==r.gene and aa.metabolite_key==r.metabolite_key
        limitations=['same_cohort_postselection','clinical_identity_and_aliquot_not_recertified','cell_source_not_external_relation_validation']
        if r.mapping_status!='DIRECT':limitations.append('CONDITIONAL_BIOCHEMISTRY_NOT_UPGRADED')
        if pd.isna(ms.p_value):limitations.append('PAIRED_AVAILABLE_N_LT8_OR_NOT_EVALUABLE')
        if pd.isna(aa.p_value):limitations.append('PRIMARY_ASSOCIATION_NOT_EVALUABLE')
        if pd.isna(vv.p_value):limitations.append('AVAILABLE_ASSOCIATION_NOT_EVALUABLE')
        elif pd.notna(aa.effect) and np.sign(vv.effect)!=np.sign(aa.effect):limitations.append('ASSOCIATION_DIRECTION_CHANGED_IN_AVAILABLE_SUBSET')
        if pd.isna(rn.p_value):limitations.append('RNA_NOT_EVALUABLE')
        if not cell.all3_same_unique_top:limitations.append('SOURCE_NOT_ALL3_CONCORDANT')
        rec={'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','relation_id':r.relation_id,'metabolite_key':r.metabolite_key,'stable_gene_id':r.stable_gene_id,'gene':r.gene,'mapping_version':r.mapping_version,'mapping_status':r.mapping_status,'metabolite_name':mm.metabolite_name,'metabolite_effect':mm.effect,'metabolite_p':mm.p_value,'metabolite_q':mm.q_value,'metabolite_n':mm.n,'metabolite_source_run':DISCOVERY.name,'association_rho':aa.effect,'association_p':aa.p_value,'association_q':aa.q_value,'association_n':aa.n,'association_source_run':PAT.name,'RNA_effect':rn.effect,'RNA_p':rn.p_value,'RNA_q':rn.q_value,'RNA_n':rn.n,'RNA_source_run':INTERNAL.name,'sc_source_run':SOURCE.name,'metabolite_evidence':evidence(mm.p_value,mm.q_value),'patient_association_evidence':evidence(aa.p_value,aa.q_value),'RNA_background':evidence(rn.p_value,rn.q_value),'cell_background':cell.cell_background,'data_limitations':';'.join(limitations),'next_action':'Discuss multi-layer evidence with source limitations;retain regardless of RNA or correlation P','status':'DONE','reason':'Integration completed;individual unevaluable measurements remain NA with reasons','RNA_test_family':rn.test_family,'association_test_family':aa.test_family,'association_status':aa.status,'association_available_status':vv.status,'RNA_status':rn.status,'relation_type':r.relation_type,'biochemical_source_url':r.source_url,'identity_restriction':r.compound_specificity,'metabolite_up_pairs':mm.n_up,'metabolite_down_pairs':mm.n_down,'metabolite_equal_pairs':mm.n_equal,'metabolite_median_delta':mm.median_delta,'metabolite_ci_lower':mm.ci_lower,'metabolite_ci_upper':mm.ci_upper,'metabolite_available_n':ms.n,'metabolite_available_p':ms.p_value,'metabolite_available_q':ms.q_value,'metabolite_available_effect':ms.effect,'association_ci_lower':aa.ci_lower,'association_ci_upper':aa.ci_upper,'association_available_n':vv.n,'association_available_rho':vv.effect,'association_available_p':vv.p_value,'association_available_q':vv.q_value,'association_available_ci_lower':vv.ci_lower,'association_available_ci_upper':vv.ci_upper,'same_association_input':aa.input_hash==vv.input_hash,'RNA_up_pairs':rn.n_up,'RNA_down_pairs':rn.n_down,'RNA_equal_pairs':rn.n_equal,'RNA_ci_lower':rn.ci_lower,'RNA_ci_upper':rn.ci_upper,'RNA_bootstrap_mean_lower':rn.bootstrap_mean_lower,'RNA_bootstrap_mean_upper':rn.bootstrap_mean_upper,'same_cohort_postselection':True}
        rec.update({k:v for k,v in cell.items() if k not in ['cancer','gene','cell_background']});rows.append(rec)
    result=pd.DataFrame(rows);write(INTEGRATION/'candidate_relations_integrated.tsv',result)
    # All historical relation IDs remain indexed separately; old relations never inherit new P/q.
    history=[]
    for _,r in old.iterrows():history.append({'gene':r.gene,'historical_relation_id':r.relation_id,'historical_mapping_version':r.mapping_version,'metabolite_key':r.metabolite_key,'metabolite_name':r.feature_name,'source_table':'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1/direct_relations_v1.tsv','new_statistics_inherited':False})
    write(INTEGRATION/'historical_relation_index.tsv',pd.DataFrame(history))
    reviewed=ROOT/'06_EXTERNAL/20260921T030642Z_historical_review_v1';histindex={}
    for name in ['historical_cell_source_label_review.tsv','historical_cptac_scope_review.tsv']:
        dd=read(reviewed/name)
        for gene in dd.gene.unique():histindex.setdefault(gene,[]).append((reviewed/name).relative_to(REPO).as_posix())
    genes=[]
    for _,g in pool.iterrows():
        sub=result[result.stable_gene_id==g.stable_gene_id];rn=RNA.loc[g.stable_gene_id];cell=sg.loc[g.stable_gene_id];legacy=old[old.gene==g.gene]
        gr={'cancer':'PDAC','stable_gene_id':g.stable_gene_id,'gene':g.gene,'in_current_pool':g.in_current_pool,'history_only':g.history_only,'relation_ids':';'.join(sorted(sub.relation_id)),'n_relations':len(sub),'RNA_background':evidence(rn.p_value,rn.q_value),'RNA_source_run':INTERNAL.name,'cell_background':cell.cell_background,'sc_source_run':SOURCE.name,'functional_source_ref':'NOT_RUN;outside source-only SOP scope','external_source_ref':';'.join(histindex.get(g.gene,[])) or 'NO_LINKED_LEGACY_EXTERNAL_GENE_RESULT','data_limitations':'Separate currentDIRECT and supplementalRNA families;no best-relation selection;clinical overlap remains unverified','next_action':'Discuss complete evidence;no mechanism claim from highest expression','status':'DONE','reason':'All current and historical genes retained;measurement status recorded separately','n_current_direct_relations':int((sub.mapping_status=='DIRECT').sum()),'n_conditional_relations':int((sub.mapping_status=='CONDITIONAL').sum()),'historical_relation_ids':';'.join(sorted(legacy.relation_id)),'n_historical_original42_relations':len(legacy),'historical_relation_index':'historical_relation_index.tsv','RNA_effect':rn.effect,'RNA_p':rn.p_value,'RNA_q':rn.q_value,'RNA_n':rn.n,'RNA_test_family':rn.test_family,'RNA_status':rn.status,'RNA_up_pairs':rn.n_up,'RNA_down_pairs':rn.n_down,'RNA_equal_pairs':rn.n_equal}
        gr.update({k:v for k,v in cell.items() if k not in ['cancer','gene','cell_background']});genes.append(gr)
    gene_result=pd.DataFrame(genes);write(INTEGRATION/'candidate_genes_integrated.tsv',gene_result)
    work=read(DISCOVERY/'metabolite_workpool.tsv');unmapped=work[~work.metabolite_key.isin(rel.metabolite_key)].copy();unmapped['mapping_status']='NO_PLANNED_RELATION';unmapped['reason']='Preserved51workpool;see full mapping feature dispositions';write(INTEGRATION/'unmapped_workpool_features.tsv',unmapped)
    assert len(result)==715 and len(gene_result)==687 and len(unmapped)==7
    assert set(work.metabolite_key)==set(result.metabolite_key)|set(unmapped.metabolite_key)
    # Verify all numeric joins by complete stable keys, not gene-only best hits.
    for col,ref,key in [('metabolite_p',m,'metabolite_key'),('metabolite_q',m,'metabolite_key'),('association_p',a,'relation_id'),('association_q',a,'relation_id'),('RNA_p',RNA,'stable_gene_id'),('RNA_q',RNA,'stable_gene_id')]:
        field='p_value' if col.endswith('_p') else 'q_value';np.testing.assert_allclose(result[col],ref.loc[result[key],field],rtol=0,atol=1e-12,equal_nan=True)
    summary={'relations':715,'genes':687,'current_direct_genes':250,'workpool':51,'unmapped_preserved':7,'all3_same_top_genes':int(sg.all3_same_unique_top.sum()),'all3_top_bootstrap_ge080_genes':int(sg.all3_top_bootstrap_ge080.sum()),'all3_top_bootstrap_ge080_current_genes':int(gene_result.loc[gene_result.in_current_pool,'all3_top_bootstrap_ge080'].sum()),'RNA_current':json.loads((INTERNAL/'summary.json').read_text())['RNA']['RNA_PAIRED_CURRENT'],'source_cohorts':json.loads((SOURCE/'summary.json').read_text())['cohorts'],'no_composite_score':True}
    write_json(INTEGRATION/'summary.json',summary);write_json(INTEGRATION/'validation.json',{'status':'PASS','all51_features_have_destination':True,'all715_relation_rows_unique':True,'all687_stable_gene_rows_unique':True,'current250_equals_DIRECT_gene_set':set(pool.loc[pool.in_current_pool,'stable_gene_id'])==set(rel.loc[rel.mapping_status=='DIRECT','stable_gene_id']),'no_RNA_or_correlation_gate_for_SC':True,'complete_key_numeric_joins_checked':True,'all362_original42_relation_ids_indexed':len(old)==362,'old_numbers_not_modified':True,'not_proven':['Clinical patient independent identity','Aliquot-level equality','Mechanism','External matched-omics replication']})
    srcfiles=[MAP/'direct_relations.tsv',MAP/'conditional_relations.tsv',MAP/'gene_pool_history_union.tsv',DISCOVERY/'metabolite_paired.tsv',DISCOVERY/'metabolite_available_sensitivity.tsv',INTERNAL/'paired_RNA.tsv',INTERNAL/'paired_RNA_history_supplement.tsv',SOURCE/'sc_celltype_profiles.tsv',SOURCE/'sc_source_stability.tsv',SOURCE/'sc_cross_study.tsv',SOURCE/'sc_gene_coverage.tsv',PAT/'tumor_association.tsv',PAT/'tumor_association_available.tsv',PAT/'conditional_tumor_association.tsv',PAT/'conditional_tumor_association_available.tsv',ROOT/'02_MAPPING/20260920T111500Z_direct_mapping_v1/direct_relations_v1.tsv']
    if not (SOURCE/'sc_celltype_profiles.tsv').exists():
        srcfiles=[p for p in srcfiles if p.name!='sc_celltype_profiles.tsv']+[SOURCE/'sc_celltype_profile_parts.tsv']+[SOURCE/f for f in read(SOURCE/'sc_celltype_profile_parts.tsv').file]
    write(INTEGRATION/'source_manifest.tsv',pd.DataFrame([{'source_id':p.stem,'path_or_url':p.relative_to(REPO).as_posix(),'sha256':digest(p),'access_scope':'PUBLIC_AGGREGATE','release':VERSION} for p in srcfiles]))
    write_json(INTEGRATION/'analysis_spec.json',{'analysis_version':VERSION,'status':'DONE_INTEGRATION_SCOPE','run_id':INTEGRATION.name,'SOP_commit':SOP_SHA,'unit':'Exact feature-stableGene-version relation and unique stable gene','selection':'All715 currentdirect/conditional and all687 historyunion;7unmapped metabolites separate','RNA_families':'250currentDIRECT and437history supplement,separateBH','SC':'cellwise log1p10k then source-label equal mean;all687;descriptive80%label;no P/q','labels':'Separate metabolite,association,RNA,cell background;no A/B/C/D or score','historical_numbers':'Immutable;old362relation ids indexed without inheriting new stats','integration_code_base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()})
    readme(INTEGRATION,'PDAC SOP v3：完整关系级与基因级比较表','51项配对P<0.05工作池；357直接、358条件关系；250当前直接基因和687历史并集。',f"715关系表、687基因表、7项无关系特征单列。三研究同一明确最高类别{summary['all3_same_top_genes']}基因，其中三边重采样保持率均≥80%为{summary['all3_top_bootstrap_ge080_genes']}基因（当前DIRECT池{summary['all3_top_bootstrap_ge080_current_genes']}）。RNA当前91项P<0.05、40项q<0.05。",'51项代谢物均未通过配对FDR；关联来自同一队列的选择后探索。80%为描述性排名稳定标签，不是患者符合率，不是机制证明。三研究按原来源标签分析，临床跨研究去重未完成。上皮/导管标签不等于恶性CNV确认。','本规范主线交付到来源与比较表；后续选题按用户另行要求，不自动机制扩展。','`python code/pdac/integrate_sop_v3_complete.py`；来源计算/独立检查及逐图输入见对应批次和figure_manifest。')
    stage(INTEGRATION,'07_INTEGRATION','715relations687genes7unmapped;separate evidence dimensions','Exact-key merges validated;legacy362relation IDs indexed;no significance intersection','User discussion;source-only SOP scope completed')
    checksum(INTEGRATION);print(json.dumps(summary))

if __name__=='__main__':main()
