"""Descriptive, exact-ID PDAC cell-line dependency check from official 24Q4 data."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
GREEN, CORAL, BLUE, GREY = "#0A5A43", "#CF765F", "#4B92A6", "#AEB9B6"


def main() -> None:
    data = pd.read_csv(TABLES / "depmap_24q4_target_models.tsv", sep="\t", low_memory=False)
    assert data.ModelID.is_unique and len(data) == 1178
    pdac = data.loc[
        data.OncotreePrimaryDisease.eq("Pancreatic Adenocarcinoma")
        & data.ModelType.eq("Cell Line")
    ].copy()
    assert len(pdac) == 45
    pdac.to_csv(TABLES / "depmap_24q4_pdac_cell_lines.tsv", sep="\t", index=False)

    summary = []
    for gene in ("NNMT", "SLC6A8", "KRAS", "RPA3", "SLC6A6"):
        effect = pdac[f"{gene}_gene_effect"].dropna()
        prob = pdac[f"{gene}_dependency_probability"].dropna()
        expr = pdac[f"{gene}_rna_log2tpm1"]
        pair = pdac[[f"{gene}_rna_log2tpm1", f"{gene}_gene_effect"]].dropna()
        rho, p = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1]) if len(pair) >= 3 else (np.nan, np.nan)
        summary.append({
            "gene": gene, "n_pdac_cell_lines": len(pdac),
            "n_gene_effect": len(effect), "effect_median": effect.median(),
            "effect_q25": effect.quantile(.25), "effect_q75": effect.quantile(.75),
            "effect_min": effect.min(), "effect_max": effect.max(),
            "n_effect_le_minus_0p5": int((effect <= -.5).sum()),
            "n_effect_le_minus_1": int((effect <= -1).sum()),
            "n_dependency_probability": len(prob),
            "n_probability_ge_0p5": int((prob >= .5).sum()),
            "n_probability_ge_0p9": int((prob >= .9).sum()),
            "n_expression": int(expr.notna().sum()),
            "expression_median_log2tpm1": expr.median(),
            "spearman_rna_effect": rho, "spearman_rna_effect_p": p,
            "n_rna_effect_pairs": len(pair),
        })
    pd.DataFrame(summary).to_csv(TABLES / "depmap_24q4_pdac_summary.tsv", sep="\t", index=False)

    all_lines = data.loc[data.ModelType.eq("Cell Line") & data.OncotreeLineage.notna()].copy()
    lineage = all_lines.groupby("OncotreeLineage", dropna=False).agg(
        n_models=("SLC6A8_gene_effect", "count"),
        slc6a8_median_effect=("SLC6A8_gene_effect", "median"),
        n_slc6a8_effect_le_minus_0p5=("SLC6A8_gene_effect", lambda x: int((x <= -.5).sum())),
        n_slc6a8_probability_ge_0p9=("SLC6A8_dependency_probability", lambda x: int((x >= .9).sum())),
    ).reset_index().sort_values("slc6a8_median_effect")
    lineage.to_csv(TABLES / "depmap_24q4_lineage_summary.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.7), gridspec_kw={"wspace": .28})
    rng = np.random.default_rng(20261013)
    for ax, genes, ylim, title in (
        (axes[0], ("NNMT", "SLC6A8"), (-.78, .44), "A  Candidate gene effects"),
        (axes[1], ("KRAS", "RPA3"), (-4.5, .45), "B  Positive-control effects"),
    ):
        ax.axhline(0, color="#7D8E89", lw=1)
        if ax is axes[0]:
            ax.axhline(-.5, color="#A8B5B0", lw=.9, ls="--")
        for i, gene in enumerate(genes):
            values = pdac[f"{gene}_gene_effect"].dropna().to_numpy()
            color = CORAL if gene == "SLC6A8" else GREEN if gene == "NNMT" else BLUE
            ax.scatter(i + rng.uniform(-.19, .19, len(values)), values,
                       s=25, c=color, alpha=.8, edgecolors="white", linewidths=.3)
            ax.plot([i-.26, i+.26], [np.median(values)]*2, color="#1C3831", lw=2.8, zorder=4)
        ax.set_xticks(range(len(genes)), genes)
        ax.set_xlim(-.48, len(genes)-.52)
        ax.set_ylim(*ylim)
        ax.set_title(title, color=GREEN, loc="left", fontweight="bold", fontsize=12)
        ax.set_ylabel("Chronos CRISPR gene effect")
        ax.grid(axis="y", alpha=.15)
    axes[0].text(.97, .06, "dashed: -0.5 triage guide", transform=axes[0].transAxes,
                 ha="right", va="bottom", fontsize=8, color="#6C7D78")

    ax = axes[2]
    pair = pdac[["CellLineName", "SLC6A8_rna_log2tpm1", "SLC6A8_gene_effect"]].dropna()
    ax.scatter(pair.SLC6A8_rna_log2tpm1, pair.SLC6A8_gene_effect,
               s=30, color=CORAL, alpha=.85, edgecolor="white", linewidth=.4)
    ax.axhline(0, color="#7D8E89", lw=1)
    ax.axhline(-.5, color="#A8B5B0", lw=.9, ls="--")
    ax.set(xlabel="SLC6A8 RNA, log2(TPM + 1)", ylabel="SLC6A8 CRISPR gene effect",
           ylim=(-.78, .25))
    ax.set_title("C  Expression vs knockout effect", color=GREEN, loc="left",
                 fontweight="bold", fontsize=12)
    ax.grid(axis="y", alpha=.15)
    for _, row in pair[pair.SLC6A8_gene_effect <= -.5].iterrows():
        ax.annotate(row.CellLineName, (row.SLC6A8_rna_log2tpm1,
                                      row.SLC6A8_gene_effect), xytext=(5, -8),
                    textcoords="offset points", fontsize=8, color="#374A43")
    fig.suptitle("PDAC cell lines: SLC6A8 is not a broad baseline growth dependency",
                 x=.045, ha="left", y=1.04, color=GREEN, fontsize=16, fontweight="bold")
    fig.text(.045, -.065,
             "DepMap 24Q4, 45 annotated pancreatic adenocarcinoma cell lines; each dot = one model, short bar = median."
             " More negative scores indicate stronger knockout effects.\n"
             "RNA matches for 43 models. A standard culture screen does not test CAF-derived metabolites or context-specific creatine dependence.",
             ha="left", fontsize=8.5, color="#61736C")
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"15_depmap_24q4_pdac_functional.{suffix}", dpi=250,
                    bbox_inches="tight")
    plt.close(fig)
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
