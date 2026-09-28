"""Focused integrity checks for public author-identity revision."""
import json,hashlib,subprocess
import numpy as np
from deliver_author_identity_v4 import SOURCE,INTEGRATION,PRE_SOURCE,PRE_INTEGRATION,REPO,MAP,read

def main():
    profile=read(SOURCE/'sc_celltype_profiles.tsv');assert len(profile)==18549 and profile.p_value.isna().all() and profile.q_value.isna().all()
    top=read(SOURCE/'sc_source_stability.tsv');assert len(top)==2061 and not top.duplicated(['cohort','stable_gene_id']).any()
    assert not profile.celltype.eq('Epithelial/Ductal').any()
    p=profile[profile.cohort=='GSE242230'];assert {'Malignant (author)','Normal epithelial (author)'}<=set(p.celltype)
    for ct,n in [('Malignant (author)',10783),('Normal epithelial (author)',350)]:assert p.loc[p.celltype==ct,'n_cells_total'].eq(n).all()
    for c in ['GSE263733','GSE278688']:
        assert not profile.loc[profile.cohort==c,'celltype'].eq('Malignant (author)').any()
        assert 'Ductal (unresolved)' in set(profile.loc[profile.cohort==c,'celltype'])
    figures=read(SOURCE/'figures/figure_manifest.tsv');pool=set(read(MAP/'gene_pool_history_union.tsv').gene)
    assert len(figures)==30
    for _,dd in figures.groupby('kind'):
        genes=';'.join(dd.genes).split(';');assert len(genes)==len(set(genes))==687 and set(genes)==pool
    n=0
    for folder in [SOURCE,INTEGRATION]:
        for _,row in read(folder/'checksums.tsv').iterrows():
            file=folder/row.file;assert hashlib.sha256(file.read_bytes()).hexdigest()==row.sha256
            raw=subprocess.check_output(['git','show',':'+file.relative_to(REPO).as_posix()],cwd=REPO);assert hashlib.sha256(raw).hexdigest()==row.sha256;n+=1
    catalog=read(INTEGRATION/'SOP_DELIVERY_STATUS.tsv');assert len(catalog)==31
    for _,row in catalog.iterrows():
        if row.data_scope!='SERVER_PRIVATE':assert all((REPO/f).is_file() for f in row.result_paths.split(';')),row.file
    subprocess.run(['git','diff','--cached','--check'],check=True,capture_output=True,cwd=REPO)
    print(json.dumps({'status':'PASS','staged_hashes':n,'profile_rows':len(profile),'all687_30plots':True,'strict_author_identity_labels':True,'31table_paths_present':True}))
if __name__=='__main__':main()
