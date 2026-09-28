"""Server-only GEO clinical linkage; public cohort aggregates only."""
from pathlib import Path
import re,json,hashlib,tarfile,xml.etree.ElementTree as ET
import pandas as pd
R=Path(__file__).resolve().parent
S=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/GSE205013/metadata/GSE205013_family.xml.tgz')
ns={'m':'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'}
with tarfile.open(S) as t:
 files=[m for m in t.getmembers() if m.name.endswith('.xml')];assert len(files)==1
 root=ET.fromstring(t.extractfile(files[0]).read())
rows=[]
for node in root.findall('m:Sample',ns):
 title=node.findtext('m:Title',namespaces=ns).strip()
 hits=re.findall(r'\bP\d{2}\b',title);assert len(hits)==1,(title,'Unresolved sample label')
 row=dict(sample=hits[0],accession=node.attrib['iid'])
 for item in node.findall('.//m:Characteristics',ns):
  key=item.attrib['tag'];assert key not in row;row[key]=item.text.strip()
 rows.append(row)
d=pd.DataFrame(rows);assert len(d)==27 and d['sample'].is_unique and d.accession.is_unique and not d.isna().any().any()
inventory=pd.read_csv(R/'private/sample_input_inventory.tsv',sep='\t')
assert set(d['sample'])==set(inventory['sample'])
d=d.merge(inventory,on='sample',validate='one_to_one')
assert set(d.tissue)=={'Primary PDAC','PDAC liver met'} and set(d.treatment)=={'Untreated','Treated'}
d['procedure_original']=d.procedure;d['procedure']=d.procedure.replace({'Liver-BIopsy':'Liver-Biopsy'})
d.to_csv(R/'private/geo_clinical_input_linkage.tsv',sep='\t',index=False)
a=d.groupby(['tissue','treatment'],sort=True).agg(n_samples=('sample','size'),input_barcodes=('n_barcodes','sum')).reset_index()
a['eligible_malignant_cells']='NOT_YET_VERIFIED';a['same_donor_pairs']='NOT_YET_VERIFIED';a.to_csv(R/'public/clinical_stratum_summary.tsv',sep='\t',index=False)
v=dict(status='PASS',scope='GEO sample-level linkage only',n_samples=27,all27_input_samples_linked=True,unique_sample_and_accession=True,missing_characteristics=False,
 metadata_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),metadata_bytes=S.stat().st_size,
 limits=['Title identifies sample; patient uniqueness relies on original study description and is not independently clinically audited','Stage is not tissue; metastatic stage primary specimens stay in primary-tissue stratum','Input barcodes are unannotated cells, not final author QC or malignant cells','Malignancy and basal/classical labels unresolved'])
(R/'public/clinical_readiness_validation.json').write_text(json.dumps(v,indent=2));print(a.to_string(index=False));print(json.dumps(v))
