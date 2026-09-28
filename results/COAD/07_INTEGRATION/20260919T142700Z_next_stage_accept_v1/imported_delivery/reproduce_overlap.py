#!/usr/bin/env python3
"""Reproduce gene-set comparison from supplied identity extracts, NOT patient analysis."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.root
    coad = load_tsv(root / "COAD_35_gene_identities_extracted.tsv")
    brca = load_tsv(root / "BRCA_174_relation_identities_extracted.tsv")
    coad_genes = {r["gene"] for r in coad}
    brca_genes = {r["gene"] for r in brca}
    ids = [r["relation_id"] for r in brca]
    if len(coad) != len(coad_genes) or len(ids) != len(set(ids)):
        raise ValueError("Duplicate gene or relation identity in the supplied extracts.")
    if (len(coad_genes), len(brca_genes), len(ids)) != (35, 117, 174):
        raise ValueError("Extract scope differs from the fixed-version delivery.")
    common = sorted(coad_genes & brca_genes)
    expected = ["BCAT2", "NNMT", "SLC6A6", "UPP2"]
    if common != expected:
        raise ValueError(f"Unexpected overlap: {common!r}")
    priority = {r["gene"] for r in coad if r["existing_decision"] == "优先核查"}
    result = {
        "scope": "COAD reviewed nominal35 versus BRCA full117, not COAD all458",
        "coad_genes": len(coad_genes), "brca_genes": len(brca_genes),
        "brca_relations": len(ids), "shared_genes": common,
        "priority9_shared": sorted(priority & brca_genes),
        "patient_statistics_run": False,
        "limitation": "Validates supplied field extracts, not original source bytes or patient statistics."
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
