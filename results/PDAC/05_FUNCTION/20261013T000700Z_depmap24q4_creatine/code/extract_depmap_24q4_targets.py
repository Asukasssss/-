"""Extract exact DepMap 24Q4 ModelID-matched target columns from official files.

Run on the server beside the three original CSV files. No row-order inference.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
EXPECTED_MD5 = {
    "Model.csv": "675210d17675f3517b0ce39a3c274f16",
    "CRISPRGeneEffect.csv": "6edf7ade09b9b34199210b559d4745d3",
    "CRISPRGeneDependency.csv": "0d3bdadf0c59264e39f7fbadf232ccdb",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": "71794802b750ce77c422dad0720a40af",
}
TARGETS = {"NNMT", "SLC6A8", "CKB", "GATM", "GAMT",
           "KRAS", "RPA3", "SLC6A6"}
EXPRESSION_ONLY = {"CKMT1A"}  # Not present in 24Q4 CRISPRGeneEffect.


def verify(path: Path, expected: str) -> None:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    assert digest.hexdigest() == expected, (path.name, digest.hexdigest())


def load_target_matrix(path: Path, suffix: str, required: set[str] = TARGETS) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0).columns
    index = header[0]
    chosen = {name: name.split(" (")[0] for name in header[1:] if name.split(" (")[0] in required}
    assert set(chosen.values()) == required, (path.name, chosen)
    data = pd.read_csv(path, usecols=[index, *chosen], low_memory=False)
    data = data.rename(columns={index: "ModelID", **{col: f"{sym}_{suffix}" for col, sym in chosen.items()}})
    assert data.ModelID.notna().all() and data.ModelID.is_unique, path.name
    return data


def main() -> None:
    for filename, digest in EXPECTED_MD5.items():
        verify(ROOT / filename, digest)
    meta = pd.read_csv(ROOT / "Model.csv", low_memory=False)
    assert meta.ModelID.is_unique
    effect = load_target_matrix(ROOT / "CRISPRGeneEffect.csv", "gene_effect")
    dependency = load_target_matrix(ROOT / "CRISPRGeneDependency.csv", "dependency_probability")
    expression = load_target_matrix(ROOT / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", "rna_log2tpm1", TARGETS | EXPRESSION_ONLY)
    result = effect.merge(meta, on="ModelID", how="left", validate="one_to_one", indicator=True)
    assert (result._merge == "both").all()
    result = result.drop(columns="_merge")
    result = result.merge(dependency, on="ModelID", how="left", validate="one_to_one", indicator=True)
    assert (result._merge == "both").all()
    result = result.drop(columns="_merge")
    result = result.merge(expression, on="ModelID", how="left", validate="one_to_one", indicator=True)
    result = result.rename(columns={"_merge": "expression_join_status"})
    keep = ["ModelID", "CellLineName", "ModelType", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode", "PrimaryOrMetastasis", "expression_join_status"]
    keep += [f"{gene}_{suffix}" for gene in sorted(TARGETS) for suffix in ("gene_effect", "dependency_probability", "rna_log2tpm1")]
    keep += [f"{gene}_rna_log2tpm1" for gene in sorted(EXPRESSION_ONLY)]
    result[keep].to_csv(ROOT / "depmap_24q4_target_models.tsv", sep="\t", index=False)
    print({"effect_models": len(effect), "expression_models": len(expression), "exact_expression_matches": int(result.expression_join_status.eq("both").sum()), "pdac_primary_disease_models": int(result.OncotreePrimaryDisease.eq("Pancreatic Adenocarcinoma").sum())})


if __name__ == "__main__":
    main()
