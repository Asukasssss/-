"""Server-only identity audit; no patient statistics or historical files changed."""
from pathlib import Path
import argparse, csv, hashlib, json, subprocess
import pandas as pd

ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', type=Path, required=True)
    out = ap.parse_args().out
    assert (out/'.running').read_text() == 'camp_sample_identity_v1'
    src = ROOT/'data/candidates/camp_primary_tissue_multicancer'
    mp = src/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    gp = ROOT/'results/collaborative/BRCA/A/20260919_followup_ER_v2/source/geo_metadata_lines.txt'
    xp = out/'source/JCI71180sd2.xlsx'; pp = out/'source/JCI71180sd.pdf'
    subprocess.run(['pdftotext','-layout',str(pp),str(out/'source/supplement_text.txt')],check=True)
    text = (out/'source/supplement_text.txt').read_text()
    records = []
    for line in text.splitlines():
        z = line.split()
        if len(z)==8 and z[0].isdigit() and z[4].isdigit() and z[4].startswith('499'):
            records.append(z)
    pairs = pd.DataFrame(records,columns=['tumor_lhc','tumor_cel','methylation','myc','tumor_met','normal_lhc','normal_cel','normal_met'])
    assert len(pairs)==67 and pairs.tumor_met.is_unique
    assert pairs.loc[pairs.normal_met.ne('no'),'normal_met'].is_unique
    long=[]
    for case,row in pairs.iterrows():
        for typ in ['tumor','normal']:
            if row[typ+'_met']!='no':
                long.append(dict(case_row=case+1,MetabID=row[typ+'_met'],pdf_lhc=row[typ+'_lhc'],pdf_TN=typ.title(),pdf_cel=row[typ+'_cel']))
    long=pd.DataFrame(long);assert len(long)==132 and long.MetabID.is_unique
    d=pd.read_excel(xp,sheet_name='OrigData',header=None)
    author=pd.DataFrame({'MetabID':d.iloc[1,10:].map(lambda x:str(int(x))), 'xlsx_lhc':d.iloc[2,10:].map(lambda x:str(int(x))), 'xlsx_tissue':d.iloc[5,10:].astype(str)}).reset_index(drop=True)
    author['xlsx_TN']=author.xlsx_tissue.map(lambda x:'Normal' if 'NORMAL' in x else 'Tumor')
    all_source=long.merge(author,on='MetabID',validate='one_to_one')
    assert len(all_source)==132 and all_source.pdf_lhc.eq(all_source.xlsx_lhc).all() and all_source.pdf_TN.eq(all_source.xlsx_TN).all()
    rows=list(csv.reader(gp.read_text().splitlines(),delimiter='\t'))
    ids=next(r[1:] for r in rows if r[0]=='!Sample_geo_accession'); assert len(set(ids))==108
    md=pd.DataFrame(index=ids)
    for row in rows:
        if row[0]=='!Sample_characteristics_ch1':
            assert len(row)==109
            md[row[1].split(': ',1)[0].lower()]=[v.split(': ',1)[1] for v in row[1:]]
    m=pd.read_csv(mp,dtype=str);m=m[m.Dataset.eq('BRCA1')].copy()
    assert len(m)==108 and all(m[c].is_unique for c in ['CommonID','MetabID','RNAID'])
    m['GSM']=m.RNAID.str.replace('.CEL.gz','',regex=False)
    assert m.GSM.isin(md.index).all()
    m['geo_lhc']=m.GSM.map(md['tumor or normal lhc'])
    m['geo_TN']=m.GSM.map(md['tissue type']).map({'Tumor':'Tumor','Non-tumor':'Normal'})
    m=m.merge(all_source,on='MetabID',validate='one_to_one')
    assert len(m)==108 and m.geo_lhc.eq(m.xlsx_lhc).all() and m.TN.eq(m.geo_TN).all()
    assert m.groupby('TN').case_row.nunique().to_dict()=={'Normal':47,'Tumor':61}
    rp=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
    metp=src/'processed_metabolomics'/m.MetabFile.iloc[0]
    rcols=pd.read_csv(rp,nrows=0).columns.astype(str)
    assert m.RNAID.isin(rcols).all()
    for sh,sub in m.groupby('MetabFile_sheet'):
        cols=pd.read_excel(metp,sheet_name=sh,nrows=0).columns.astype(str)
        assert sub.MetabID.isin(cols).all()
    conflict=m[m.TN.ne(m.pdf_TN)];assert len(conflict)==1
    consistent=m[m.TN.eq(m.pdf_TN)]
    t=set(consistent.loc[consistent.TN.eq('Tumor'),'case_row']);n=set(consistent.loc[consistent.TN.eq('Normal'),'case_row'])
    summary=dict(mapped_specimens=108,camp_tumor=61,camp_normal=47,crossomics_LHC_matches=108,
                 original_author_case_rows=67,original_author_pairs=65,
                 matched_subset_author_cases=int(m.case_row.nunique()),consistent_pairs=len(t&n),
                 consistent_tumor_only=len(t-n),consistent_normal_only=len(n-t),tissue_conflicts=len(conflict),
                 pdf_RNA_absent_but_GEO_available=int(m.pdf_cel.eq('no').sum()),
                 patient_identity_basis='Author case-row reconstruction, not genotype verification',
                 new_patient_statistics=0,historical_values_modified=False)
    assert summary['consistent_pairs']==45 and summary['matched_subset_author_cases']==63
    m.to_csv(out/'private/audited_mapping.tsv',sep='\t',index=False)
    pairs.to_csv(out/'private/author_pair_table.tsv',sep='\t',index=False)
    (out/'public/validation.json').write_text(json.dumps(summary,indent=2))
    pd.DataFrame([dict(check=k,value=v) for k,v in summary.items()]).to_csv(out/'public/results.tsv',sep='\t',index=False)
    paths=[mp,gp,xp,pp,rp,metp,out/'source/conflict_GEO_full.soft',Path(__file__)]
    pd.DataFrame([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]).to_csv(out/'public/source_manifest.tsv',sep='\t',index=False)
    spec=dict(version='camp_sample_identity_v1',method='Exact CAMP MetabID to original XLSX SAMPLE_ID; exact RNA GSM to GEO LHC; author PDF case rows define pairs',
              inference_from_order_or_numeric_suffix=False,statistics_recomputed=False,
              sources=['https://www.jci.org/articles/view/71180/sd/1','https://www.jci.org/articles/view/71180/sd/2','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37751'],
              software=dict(pandas=pd.__version__),limitation='Original source tissue conflict unresolved; no new identity assay')
    (out/'public/analysis_spec.json').write_text(json.dumps(spec,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
