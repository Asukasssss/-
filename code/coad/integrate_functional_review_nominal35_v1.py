#!/usr/bin/env python3
"""Append COAD functional evidence without recalculating or editing source statistics.

Requires Python 3.9+. Uses only standard-library modules.
The source files must match the pinned Git blob IDs in analysis_spec.json.
No patient matrix, network access, git write, or DepMap calculation is performed.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, pathlib, sys

PREFIX = "func35_"
ANN_FIELDS = [
    "review_status", "annotation_level", "pair_scope", "evidence_class", "research_decision",
    "functional_conclusion", "metabolic_link", "counterevidence_or_limit",
    "research_increment", "next_action", "evidence_ids", "source_urls",
    "depmap_status", "depmap_score", "depmap_n_models"
]
MAP = {
    "evidence_class": "evidence_class_cn",
    "research_decision": "research_decision_cn",
    "functional_conclusion": "functional_conclusion_cn",
    "metabolic_link": "metabolic_link_cn",
    "counterevidence_or_limit": "counterevidence_or_limit_cn",
    "research_increment": "research_increment_cn",
    "next_action": "next_action_cn",
    "evidence_ids": "evidence_ids",
    "source_urls": "evidence_urls"
}

def read_tsv(path: pathlib.Path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        fields = list(rd.fieldnames or [])
        rows = list(rd)
    if not fields or any(None in r for r in rows):
        raise ValueError(f"Invalid TSV structure: {path}")
    return fields, rows

def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def git_blob_sha(path: pathlib.Path) -> str:
    data = path.read_bytes()
    # Text checkouts may use CRLF. Git's canonical blob here is LF.
    data = data.replace(b"\r\n", b"\n")
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()

def key(row):
    return row["metabolite_key"], row["gene"]

def enrich(fields, rows, review, selected_keys=None):
    appended = [PREFIX + x for x in ANN_FIELDS]
    if set(appended) & set(fields):
        raise ValueError("This evidence prefix already exists. Use a new version; do not overwrite.")
    out = []
    for old in rows:
        record = dict(old)
        annotation = {PREFIX + x: "" for x in ANN_FIELDS}
        r = review.get(old.get("gene", ""))
        annotation[PREFIX+"review_status"] = (
            "REVIEWED_IN_35_GENE_BATCH" if r else "NOT_REVIEWED_THIS_BATCH")
        annotation[PREFIX+"annotation_level"] = "GENE_LEVEL_LITERATURE"
        annotation[PREFIX+"pair_scope"] = ("INITIAL39_NOMINAL_PAIR" if selected_keys is not None and key(old) in selected_keys else "NOT_REVIEWED_AS_A_PAIR_IN_THIS_BATCH")
        annotation[PREFIX+"depmap_status"] = "NOT_RUN/DEFERRED"
        if r:
            for target, source in MAP.items():
                annotation[PREFIX+target] = r.get(source, "") or ""
        record.update(annotation)
        # Preserve EVERY original string field, not merely P/q.
        if any(record[c] != old[c] for c in fields):
            raise AssertionError("Source cell changed.")
        out.append(record)
    return fields + appended, out

def self_test():
    fields = ["metabolite_key","gene","effect","p_value","q_value","status"]
    original = [
        dict(zip(fields,["KEGG:X","G","0.123456789","0.04","0.9","DONE"])),
        dict(zip(fields,["KEGG:Y","UNREVIEWED","NA","NA","NA","NOT_EVALUABLE"]))
    ]
    reviews = {"G":{"evidence_class_cn":"test", "evidence_ids":"E-test"}}
    headers, output = enrich(fields, original, reviews)
    assert len(output)==len(original)
    assert [{c:r[c] for c in fields} for r in output]==original
    assert output[1][PREFIX+"review_status"]=="NOT_REVIEWED_THIS_BATCH"
    assert output[1][PREFIX+"depmap_score"]==""
    assert output[0][PREFIX+"evidence_ids"]=="E-test"
    print(json.dumps({"self_test":"PASS","real_repository_join":"NOT_RUN"}))
    return 0

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=pathlib.Path)
    p.add_argument("--evidence-root", type=pathlib.Path)
    p.add_argument("--out", type=pathlib.Path)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args(argv)
    if a.self_test:
        return self_test()
    if not all((a.repo, a.evidence_root, a.out)):
        p.error("--repo, --evidence-root and --out are required")
    repo, eroot, out = a.repo.resolve(), a.evidence_root.resolve(), a.out.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {out}")
    spec = json.loads((eroot/"analysis_spec.json").read_text(encoding="utf-8"))
    reviews = json.loads((eroot/"results"/"gene_function_review.json").read_text(encoding="utf-8"))
    review = {r["gene"]:r for r in reviews}
    if len(reviews)!=35 or len(review)!=35:
        raise ValueError("Expected 35 unique reviewed genes.")
    np = repo / spec["source_nominal_path"]
    fp = repo / spec["source_full_catalog_path"]
    for path, expected in [
        (np, spec["source_nominal_git_blob_sha1"]),
        (fp, spec["source_full_catalog_git_blob_sha1"])
    ]:
        got = git_blob_sha(path)
        if got != expected:
            raise ValueError(f"Source revision mismatch: {path}; got {got}, expected {expected}. "
                             "Do not silently apply the old review to a changed source.")
    nf, nominal = read_tsv(np)
    ff, catalog = read_tsv(fp)
    _, display = read_tsv(eroot/"source"/"nominal39_display_extract.tsv")
    if len(nominal)!=39 or len(catalog)!=974:
        raise ValueError("Expected 39 nominal rows and 974 full catalog rows.")
    if {r["gene"] for r in nominal} != set(review):
        raise ValueError("Nominal gene set differs from reviewed 35.")
    dk = {key(r): r for r in display}
    if len(dk)!=39 or {key(r) for r in nominal}!=set(dk):
        raise ValueError("Display extract and source pair identities differ.")
    for row in nominal:
        d = dk[key(row)]
        for source, shown in [("effect","rho_display"),("p_value","p_display"),("q_value","q_display")]:
            if abs(float(row[source])-float(d[shown]))>0.000051:
                raise ValueError(f"Rounded display mismatch for {key(row)}: {source}")
        if int(float(row["n"])) != int(float(d["n_display"])):
            raise ValueError(f"Sample count mismatch for {key(row)}")
    before = {str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (np,fp)}
    nk = {key(r) for r in nominal}
    nh, nr = enrich(nf, nominal, review, nk)
    fh, fr = enrich(ff, catalog, review, nk)
    out.mkdir(parents=True, exist_ok=False)
    write_tsv(out/"nominal39_annotated_exact.tsv", nh, nr)
    write_tsv(out/"catalog974_annotated_exact.tsv", fh, fr)
    write_tsv(out/"gene_function_review.tsv", list(reviews[0]), [
        {k:("" if v is None else v) for k,v in r.items()} for r in reviews])
    after = {str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (np,fp)}
    if before != after:
        raise AssertionError("Source file hashes changed.")
    validation = {
        "status":"PASS", "source_blob_ids_verified":True,
        "original_columns_preserved":True,"catalog_rows":len(fr),
        "nominal_rows":len(nr),"reviewed_genes":len(review),
        "patient_tests_recomputed":0,"depmap_computed":False,
        "source_sha256_before":before,"source_sha256_after":after,
        "literature_review":"Targeted evidence review, not systematic completeness certification",
        "git_upload":"NOT_PERFORMED_BY_SCRIPT"
    }
    (out/"join_validation.json").write_text(json.dumps(validation,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    # Emit an append proposal, never overwrite the active coordinator's index.
    rel_out = out.relative_to(repo).as_posix() if out.is_relative_to(repo) else str(out)
    rows=[]
    for stage in ("05_FUNCTION","07_INTEGRATION"):
        rows.append({
            "cancer":"COAD","stage_id":stage,"run_id":out.name,
            "analysis_version":spec["analysis_version"],"status":"PARTIAL",
            "scope":"35 genes/39 nominal relations functionally reviewed; full974 retained",
            "result_path":rel_out,"code_path":"code/coad/integrate_functional_review_nominal35_v1.py",
            "git_branch":"analysis/coad-functional-review-20260919",
            "reason":"Scoped batch complete; wider candidate review and external/causal validation incomplete",
            "next_action":"Review context-specific questions; no repeated significance selection"})
    write_tsv(out/"stage_index_append.tsv",list(rows[0]),rows)
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    return 0

if __name__=="__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(1)
