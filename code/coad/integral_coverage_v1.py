"""Coverage only: explicit sample map, no correlations or significance selection."""
import pathlib,json,csv,math,collections,openpyxl,hashlib
root=pathlib.Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
run=root/'results/collaborative/COAD/B/20260921T050501Z_integral_acquire_v1'
source=root/'data/candidates/coad_integral_20260921/ac4c04421_si_003.xlsx'
w=openpyxl.load_workbook(source,read_only=True,data_only=True)
targets=json.loads((run/'targets.json').read_text())
samples=list(w['A_sample_information'].values);meta=samples[1:]
assert len(meta)==12 and len(set(x[0]for x in meta))==12
types=collections.Counter(x[2]for x in meta)
tumor=[x[0]for x in meta if 'tumor' in x[2].lower()and 'normal' not in x[2].lower()]
assert len(tumor)==6 and len({x[1]for x in meta if x[0]in tumor})==6
aliases={
'C02918':('N-Methylnicotinamide','AMBIGUOUS_NAME'),
'C04501':('N-Acetylglucosamine 1-Phosphate','NAME_COMPATIBLE'),
'C01042':('N-Acetylaspartate','NAME_COMPATIBLE'),
'C00019':('S-Adenosyl-L-Methionine','NAME_COMPATIBLE'),
'C00147':('Adenine','NAME_COMPATIBLE'),
'C03794':('Adenylocuccinic Acid','TYPO_REQUIRES_CONFIRMATION'),
'C00152':('L-Asparagine Anhydrous','NAME_COMPATIBLE'),
'C00475':('Cytidine','NAME_COMPATIBLE'),
'C00051':('Glutathione Reducedform','NAME_COMPATIBLE'),
'C00670':('Sn-Glycero-3-Phosphocholine','NAME_COMPATIBLE'),
'C00242':('Guanine','NAME_COMPATIBLE'),
'C00387':('Guanosine','NAME_COMPATIBLE'),
'C00388':('Histamine','NAME_COMPATIBLE'),
'C00130':("Inosine 5'-monophosphate",'NAME_COMPATIBLE'),
'C00407':('L-Isoleucine','NAME_COMPATIBLE'),
'C00123':('DL-Leucine','STEREOCHEMISTRY_UNRESOLVED'),
'C00148':('L-Proline','NAME_COMPATIBLE'),
'C00245':('Taurine','NAME_COMPATIBLE'),
'C00299':('Uridine','NAME_COMPATIBLE'),
'C00015':('Uridine 5’‑Diphosphate'.replace('‑','-'),'NAME_COMPATIBLE')}
tables={};coverage={};structure=[]
for sheet,keycol,start in [('B_metabolome',0,1),('D_proteome_protein',1,2),('G_transcriptome_gene',0,2)]:
 data=list(w[sheet].values);header=data[0];assert set(header[start:])=={x[0]for x in meta}
 ti=[header.index(x)for x in tumor];group=collections.defaultdict(list)
 for row in data[1:]:group[str(row[keycol])].append(row)
 tables[sheet]=group
 def finite(x):return isinstance(x,(int,float))and math.isfinite(x)
 coverage[sheet]={k:dict(records=len(v),finite_tumor=min(sum(finite(row[i])for i in ti)for row in v),nonzero_tumor=min(sum(finite(row[i])and row[i]!=0 for i in ti)for row in v))for k,v in group.items()}
 structure.append(dict(sheet=sheet,features=len(data)-1,measurement_columns=len(header)-start,exact_sample_set_match=True))
def write(name,rows):
 with (run/name).open('w',newline='')as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
met=[]
for key,name in sorted(set((x['metabolite_key'],x['metabolite_name'])for x in targets)):
 label,status=aliases.get(key.split(':')[-1],('NA','NOT_FOUND_IN_NAME_SCREEN'))
 c=coverage['B_metabolome'].get(label)
 if label!='NA':assert c,repr(label)
 met.append(dict(metabolite_key=key,metabolite_name=name,source_name=label,match_status=status,source_rows=c['records']if c else 0,finite_tumor=c['finite_tumor']if c else 'NA',nonzero_tumor=c['nonzero_tumor']if c else 'NA',identity_certified=False))
genes=[]
for g in sorted(set(x['gene']for x in targets)):
 rec=dict(gene=g)
 for layer,sheet in [('rna','G_transcriptome_gene'),('protein','D_proteome_protein')]:
  c=coverage[sheet].get(g);rec.update({layer+'_exact_rows':c['records']if c else 0,layer+'_finite_tumor':c['finite_tumor']if c else 'NA',layer+'_nonzero_tumor':c['nonzero_tumor']if c else 'NA'})
 genes.append(rec)
m={x['metabolite_key']:x for x in met};gg={x['gene']:x for x in genes}
rels=[dict(**t,metabolite_match_status=m[t['metabolite_key']]['match_status'],source_metabolite=m[t['metabolite_key']]['source_name'],rna_rows=gg[t['gene']]['rna_exact_rows'],protein_rows=gg[t['gene']]['protein_exact_rows'],analysis_status='NOT_RUN',reason='Coverage only; chemical identity, expression scale and clinical eligibility require lock')for t in targets]
write('metabolite_coverage23.tsv',met);write('gene_coverage35.tsv',genes);write('relation_coverage39.tsv',rels)
summary=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),sample_types=dict(types),patients=len(set(x[1]for x in meta)),tumor_samples=len(tumor),unique_tumor_patients=6,structure=structure,metabolite_status=dict(collections.Counter(x['match_status']for x in met)),rna_genes=sum(x['rna_exact_rows']>0 for x in genes),protein_genes=sum(x['protein_exact_rows']>0 for x in genes),new_tests=0)
(run/'coverage_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
