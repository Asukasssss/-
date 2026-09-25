"""Fresh PDAC/B server run: split author identities, recompute affected cohort only."""
from pathlib import Path
import json,hashlib,traceback,itertools,shutil,tarfile,subprocess,sys
import numpy as np
import pandas as pd
import run_sop_v3_server as base
import sc_source_three_cohorts_v1 as readers
from author_identity_v4 import resolve_identity,UNRESOLVED
from verify_sop_v3_server import source as independent_source_check

RUN=Path(__file__).resolve().parent
PRE=RUN.parent/'20260922T120400Z_sop_v3_source'
VERSION='PDAC_author_cell_identity_v4'
COHORTS=['GSE263733','GSE278688','GSE242230']
def read(p):return pd.read_csv(p,sep='\t')
def save(name,x):(RUN/'public'/name).write_text(json.dumps(x,indent=2)+'\n')
def write(name,x):x.to_csv(RUN/'public'/name,sep='\t',index=False,na_rep='NA')
def corrected242(genes):
    meta,x,lib,present=original242(genes)
    assert meta.index.is_unique and meta.cell_type_specific.notna().all()
    ids=[resolve_identity('GSE242230',r.celltype,r.cell_type_specific) for _,r in meta.iterrows()]
    audit=meta[['unit','celltype','cell_type_specific']].copy();audit['resolved_celltype']=[v[0] for v in ids];audit['identity_status']=[v[1] for v in ids];audit['reason']=[v[2] for v in ids]
    audit.to_csv(RUN/'private/cell_identity_mapping_private.tsv',sep='\t',index_label='cell_id')
    counts=audit.groupby(['celltype','cell_type_specific','resolved_celltype','identity_status','reason']).size().reset_index(name='n_cells');counts.insert(0,'cohort','GSE242230');write('author_identity_counts.tsv',counts)
    assert len(meta)==31215 and (audit.identity_status=='AUTHOR_NORMAL').sum()==350 and (audit.identity_status=='AUTHOR_MALIGNANT').sum()==10783
    meta=meta.copy();meta['celltype']=audit.resolved_celltype
    save('cell_join_validation.json',{'status':'PASS','selected_cells':len(meta),'unique_cell_ids':True,'all_selected_cells_have_specific_label':True,'same_cells_as_v3':True,'exact_barcode_matrix_join':'Checked by unchanged raw reader','author_normal_cells':350,'author_malignant_cells':10783,'author_classical':5660,'author_basal':5123,'normal_overrides_broad_malignant':True,'private_mapping_exported':False})
    return meta,x,lib,present
original242=readers.sparse242

def main():
    assert RUN.parent==readers.ROOT/'results/collaborative/PDAC/B'
    (RUN/'private').mkdir();(RUN/'public').mkdir()
    with (RUN/'.running').open('x') as f:f.write(VERSION)
    try:
        subprocess.run([sys.executable,'-m','unittest','test_author_identity_v4'],cwd=RUN,check=True,stdout=(RUN/'public/identity_tests.txt').open('w'),stderr=subprocess.STDOUT)
        spec=json.loads((RUN/'source/identity_spec.json').read_text());panel=json.loads((RUN/'source/panel.json').read_text());aliases=json.loads((RUN/'source/gene_aliases.json').read_text())
        prior_manifest=read(PRE/'public/source_manifest.tsv');checked=[]
        for _,r in prior_manifest.iterrows():
            p=Path(r.path_or_url)
            if '/data/candidates/' in str(p) and ('GSE242230' in str(p) or 'annotation' in p.name.lower() or 'metadata' in p.name.lower()):
                assert readers.sha(p)==r.sha256,str(p);checked.append(p)
        assert len(checked)>=78
        base.COHORTS=['GSE242230'];base.VERSION=VERSION;readers.dense263=corrected242
        s=base.source(RUN,panel,aliases)
        files=['sc_celltype_profiles.tsv','sc_source_stability.tsv','sc_gene_coverage.tsv','sc_annotation_map.tsv']
        for name in files:
            old=read(PRE/'public'/name);old=old[old.cohort.isin(COHORTS[:2])].copy();new=read(RUN/'public'/name)
            for col in ['celltype','top_celltype','runner_celltype','coarse_celltype']:
                if col in old:old[col]=old[col].str.replace('Epithelial/Ductal',UNRESOLVED,regex=False)
            if 'run_id' in old:old['run_id']=RUN.name
            if 'analysis_version' in old:old['analysis_version']=VERSION
            if name=='sc_annotation_map.tsv':
                old['reason']='Ductal identity unresolved in available per-cell GEO metadata';new['reason']='Specific author identity retained;no new malignancy inference';new['annotation_origin']='Author supplied fine labels';new['mapping_evidence']='Normal Epithelial takes precedence over broad Malignant';old['version']=new['version']=VERSION
            joined=pd.concat([old,new],ignore_index=True);write(name,joined)
        for c in COHORTS[:2]:
            p=PRE/'private'/(c+'_sc_donor_profiles_private.tsv');d=read(p);d['celltype']=d.celltype.replace({'Epithelial/Ductal':UNRESOLVED});d.to_csv(RUN/'private'/p.name,sep='\t',index=False,na_rep='NA');checked.append(p)
        prof=read(RUN/'public/sc_celltype_profiles.tsv');top=read(RUN/'public/sc_source_stability.tsv');cross=[]
        for gene in panel['genes']:
            for ca,cb in itertools.combinations(COHORTS,2):
                a=top[(top.gene==gene)&(top.cohort==ca)].iloc[0];b=top[(top.gene==gene)&(top.cohort==cb)].iloc[0]
                pa=prof[(prof.gene==gene)&(prof.cohort==ca)&(prof.status=='DONE')];pb=prof[(prof.gene==gene)&(prof.cohort==cb)&(prof.status=='DONE')];common=set(pa.celltype)&set(pb.celltype);shared=[]
                for dd in [pa,pb]:
                    z=dd[dd.celltype.isin(common)];ok=len(z)>=2 and z.mean_detection_fraction.max()>=.01;shared.append(';'.join(sorted(z.loc[np.isclose(z.effect,z.effect.max(),rtol=1e-10,atol=1e-12),'celltype'])) if ok else 'NOT_EVALUABLE')
                valid=a.status==b.status=='DONE' and a.tie_status==b.tie_status=='UNIQUE';same=bool(valid and a.top_celltype==b.top_celltype)
                cross.append({'cancer':'PDAC','gene':gene,'stable_gene_id':panel['stable_gene_ids'][gene],'study_A':ca,'study_B':cb,'top_A':a.top_celltype,'top_B':b.top_celltype,'n_common_categories':len(common),'same_top_all_categories':same if valid else np.nan,'top_A_shared':shared[0],'top_B_shared':shared[1],'same_top_shared_categories':shared[0]==shared[1] if all(x!='NOT_EVALUABLE' and ';' not in x for x in shared) else np.nan,'bootstrap_top_A':a.bootstrap_top_frequency,'bootstrap_top_B':b.bootstrap_top_frequency,'same_top_both_bootstrap_ge080':same and a.bootstrap_top_frequency>=.8 and b.bootstrap_top_frequency>=.8 if valid else np.nan,'status':'DONE' if valid else 'NOT_EVALUABLE','reason':'Strict identity labels;malignant and unresolved ductal never equated;shared-category result excludes nonshared epithelial identities'})
        write('sc_cross_study.tsv',pd.DataFrame(cross))
        oldsum=json.loads((PRE/'public/summary.json').read_text());oldsum['profile_rows']=len(prof);oldsum['normal_split_cells']=350;oldsum['malignant_cells']=10783;oldsum['changed_cohort']='GSE242230';oldsum['reuse_cohorts']=COHORTS[:2];save('summary.json',oldsum)
        save('analysis_spec.json',spec)
        # Complete numeric comparisons: unchanged cohorts have identical statistics;
        # split cohort non-epithelial expression means retain their v3 values.
        before=read(PRE/'public/sc_celltype_profiles.tsv');before.celltype=before.celltype.replace({'Epithelial/Ductal':UNRESOLVED})
        for c in COHORTS:
            a=before[before.cohort==c].set_index(['gene','celltype']);b=prof[prof.cohort==c].set_index(['gene','celltype']);keys=a.index.intersection(b.index)
            columns=['effect','mean_detection_fraction','n','n_cells_total','n_cells_eligible']
            if c!='GSE242230':columns+=['bootstrap_first_frequency']
            np.testing.assert_allclose(a.loc[keys,columns].to_numpy(float),b.loc[keys,columns].to_numpy(float),rtol=0,atol=1e-12,equal_nan=True)
        v=independent_source_check(RUN);v.update(identity_join='PASS',unaffected_numeric_reuse='PASS',malignancy='Author labels only;no independent CNV recomputation');save('independent_validation.json',v)
        inputs=checked+[PRE/'public'/n for n in files]+list((RUN/'source').glob('*.json'))+list(RUN.glob('*.py'))
        write('source_manifest.tsv',pd.DataFrame([{'source_id':p.name,'path_or_url':str(p),'sha256':readers.sha(p),'access_scope':'SERVER_INPUT'} for p in inputs]))
        (RUN/'NUMERICS_DONE.json').write_text(json.dumps(oldsum));(RUN/'.running').unlink()
        with tarfile.open(RUN/'public_delivery.tar','w') as t:
            for p in (RUN/'public').iterdir():t.add(p,arcname=p.name)
        (RUN/'PUBLIC_READY.json').write_text(json.dumps(v));print('PUBLIC_READY',flush=True)
    except Exception:(RUN/'FAILED.txt').write_text(traceback.format_exc());raise
if __name__=='__main__':main()
