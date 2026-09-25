"""Server-only coverage and donor-equal expression; never exports individual rows."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
VERSION = 'COAD_cell_origin_v1'
PREFIX = 'cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
KEYS = ['study', 'tissue', 'enrichment', 'technology', 'annotation_level', 'cell_type']

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def vector(path):
    # Author readDataRobj.m/fastTxtRead uses headerless dDvec files.
    with gzip.open(path, 'rt') as f:
        return np.array([x.rstrip('\r\n').strip('"') for x in f], dtype=object)

def load_metadata(d):
    a = pd.read_csv(d / 'lee_annotation.txt.gz', sep='\t', dtype=str)
    assert a.Index.is_unique and not a.isna().any().any()
    assert a.groupby('Sample')[['Patient', 'Class']].nunique().max().max() == 1
    lee = pd.DataFrame(dict(cell_id=a.Index, patient=a.Patient, specimen=a.Sample,
        library='NA', study='Lee2020_GSE132465', tissue=a.Class,
        enrichment='NO_CD45_REPORTED_FICOLL_PURIFIED', technology='10x_3prime_v2',
        lineage=a.Cell_type, subtype=a.Cell_subtype))
    m = pd.read_csv(d / 'pelka_author_metadata.csv.gz', dtype=str, keep_default_na=False)
    c = pd.read_csv(d / 'pelka_author_annot.csv.gz', dtype=str, keep_default_na=False)
    ids = vector(d / 'colon10x_default_dDvec_sampleID.csv.gz')
    libs = vector(d / 'colon10x_default_dDvec_batchID.csv.gz')
    assert len(m) == len(c) == len(ids) == len(libs) == 371223
    assert len(set(ids)) == len(ids)
    # Synchronized arrays in the original author's exported object, not an
    # inferred patient join. readDataRobj.m and Figure_1.m document this contract.
    pelka = pd.DataFrame(dict(cell_id=ids, patient=m.PID, specimen=m.PatientTypeID,
        library=libs, study='Pelka2021_GSE178341', tissue=m.SPECIMEN_TYPE,
        enrichment=m.PROCESSING_TYPE, technology=m.SINGLECELL_TYPE,
        lineage=c.clTopLevel, midway=c.clMidway, subtype=c.clFull))
    assert not pelka[['patient', 'specimen', 'enrichment', 'technology']].isin(['', 'NA']).any().any()
    assert pelka.groupby('library')[['patient', 'specimen', 'tissue', 'enrichment', 'technology']].nunique().max().max() == 1
    return {'Lee2020_GSE132465': lee, 'Pelka2021_GSE178341': pelka}

def levels(frame):
    for col in ['lineage', 'midway', 'subtype']:
        if col in frame:
            f = frame.copy()
            f['annotation_level'] = col
            f['cell_type'] = f[col]
            yield f

def base(run, key, kind, n, status='DONE', reason='Descriptive only; no hypothesis tests'):
    row = dict.fromkeys(PREFIX, 'NA')
    row.update(cancer='COAD', cohort=key[0], stage_id='06_EXTERNAL', run_id=run,
        analysis_version=VERSION, analysis_type=kind, unit='patient', n=n,
        status=status, reason=reason, source_id=key[0])
    row.update(zip(KEYS, key))
    return row

def write_table(rows, path):
    f = pd.DataFrame(rows)
    extra = [x for x in f if x not in PREFIX]
    f.reindex(columns=PREFIX + extra).to_csv(path, sep='\t', index=False, na_rep='NA')

def coverage(data, out):
    rows, audit = [], {}
    for study, original in data.items():
        audit[study] = {'cells': len(original), 'patients': original.patient.nunique(),
            'specimens': original.specimen.nunique(), 'libraries': (original.loc[original.library != 'NA', 'library'].nunique() if (original.library != 'NA').any() else None),
            'library_id_missing': int((original.library == 'NA').sum()),
            'tissue_cells': original.tissue.value_counts().to_dict(),
            'enrichment_cells': original.enrichment.value_counts().to_dict()}
        for f in levels(original):
            for key, group in f.groupby(KEYS, dropna=False, sort=True):
                v = group.groupby('patient').size()
                row = base(out.name, key, 'cell_coverage', len(v))
                row.update(n_cells=int(v.sum()), n_patients_ge20=int((v >= 20).sum()),
                    patient_cells_min=int(v.min()), patient_cells_q25=float(v.quantile(.25)),
                    patient_cells_median=float(v.median()), patient_cells_q75=float(v.quantile(.75)),
                    patient_cells_max=int(v.max()), n_specimens=group.specimen.nunique(),
                    composition_eligible=key[2] in ['unsorted', 'NO_CD45_REPORTED_FICOLL_PURIFIED'])
                rows.append(row)
    write_table(rows, out / 'coverage_by_lineage.tsv')
    (out / 'metadata_audit.json').write_text(json.dumps(audit, indent=2))
    return audit

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--mode', choices=['coverage'], required=True)
    a = p.parse_args()
    assert a.data_dir.resolve().parent == ROOT / 'data/candidates'
    assert a.out.resolve().parent == ROOT / 'results/collaborative/COAD/B'
    assert (a.out / '.running').exists()
    data = load_metadata(a.data_dir)
    if a.mode == 'coverage':
        assert not (a.out / 'coverage_by_lineage.tsv').exists()
        print(json.dumps(coverage(data, a.out)), flush=True)

if __name__ == '__main__':
    main()
