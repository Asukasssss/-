"""Exploratory CKB/creatine-pathway context in PDAC cell lines, exact ModelID joins."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, spearmanr

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
data = pd.read_csv(TABLES / "depmap_24q4_pdac_cell_lines.tsv", sep="\t")
assert len(data) == 45 and data.ModelID.is_unique

summary = []
for gene in ("CKB", "CKMT1A", "GATM", "GAMT"):
    rna = data[f"{gene}_rna_log2tpm1"].dropna()
    row = {"gene": gene, "n_pdac_cell_lines": len(data), "n_rna": len(rna),
           "rna_median_log2tpm1": rna.median(),
           "rna_range_min": rna.min(), "rna_range_max": rna.max()}
    if f"{gene}_gene_effect" in data:
        effect = data[f"{gene}_gene_effect"].dropna()
        probability = data[f"{gene}_dependency_probability"].dropna()
        row.update(n_effect=len(effect), effect_median=effect.median(),
                   n_effect_le_minus_0p5=int((effect <= -.5).sum()),
                   n_probability_ge_0p9=int((probability >= .9).sum()),
                   effect_status="evaluated")
    else:
        row.update(n_effect=0, effect_median=np.nan,
                   n_effect_le_minus_0p5=np.nan,
                   n_probability_ge_0p9=np.nan,
                   effect_status="not_measured_in_crispr_matrix")
    summary.append(row)
pd.DataFrame(summary).to_csv(TABLES / "depmap_24q4_creatine_pathway_summary.tsv",
                             sep="\t", index=False)

comparisons = [
    ("CKB RNA", "SLC6A8 CRISPR effect", "CKB_rna_log2tpm1", "SLC6A8_gene_effect"),
    ("CKB RNA", "SLC6A8 RNA", "CKB_rna_log2tpm1", "SLC6A8_rna_log2tpm1"),
    ("CKMT1A RNA", "SLC6A8 CRISPR effect", "CKMT1A_rna_log2tpm1", "SLC6A8_gene_effect"),
    ("GATM RNA", "GAMT RNA", "GATM_rna_log2tpm1", "GAMT_rna_log2tpm1"),
]
rows = []
for label_x, label_y, col_x, col_y in comparisons:
    pair = data[[col_x, col_y]].dropna()
    rho, p = spearmanr(pair[col_x], pair[col_y])
    rows.append({"measure_x": label_x, "measure_y": label_y,
                 "n_exact_models": len(pair), "spearman_rho": rho,
                 "p_nominal": p})
q = false_discovery_control(np.array([x["p_nominal"] for x in rows]), method="bh")
for row, qi in zip(rows, q):
    row["q_bh_four_exploratory_tests"] = qi
pd.DataFrame(rows).to_csv(TABLES / "depmap_24q4_creatine_context_correlations.tsv",
                          sep="\t", index=False)
print(pd.DataFrame(summary).to_string(index=False))
print(pd.DataFrame(rows).to_string(index=False))
