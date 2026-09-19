#!/usr/bin/env python3
"""Multi-cancer CPTAC support analysis for the frozen inosine hypothesis.

The analysis is intentionally separate from CAMP_Phase1_v1.0 and the existing
PDAC run.  It uses exact PDC aliquot/sample metadata and explicit GDC
aliquot-to-sample-type mapping; no sample order, imputation, purity inference,
or direct metabolite--gene correlation is used.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import html
import json
import math
import re
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import requests
import yaml
from scipy import stats


PDC_URL = "https://pdc.cancer.gov/graphql"
GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data"
ANALYSIS_NAME = "Inosine_pan_cancer_CPTAC_v0.1"
TUMOR = "Primary Tumor"
NORMAL = "Solid Tissue Normal"


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: Sequence[float]) -> List[float]:
    array = np.asarray(values, dtype=float)
    result = np.full(array.shape, np.nan, dtype=float)
    finite = np.isfinite(array)
    if not finite.any():
        return result.tolist()
    p = array[finite]
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    result[finite] = adjusted
    return result.tolist()


def pdc_query(query: str) -> Dict[str, Any]:
    response = requests.post(PDC_URL, json={"query": query}, timeout=300)
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError("PDC GraphQL error: " + json.dumps(payload["errors"], ensure_ascii=False)[:2000])
    return payload["data"]


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def text_or_na(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return str(value)


def gdc_case_id_from_refs(row: Dict[str, Any]) -> Optional[str]:
    for reference in row.get("externalReferences") or []:
        location = reference.get("reference_entity_location") or ""
        match = re.search(r"/cases/([0-9a-fA-F-]{36})", location)
        if match:
            return match.group(1).lower()
    return None


def fetch_biospecimen(study_id: str) -> pd.DataFrame:
    query = f'''query {{ biospecimenPerStudy(pdc_study_id: "{study_id}", acceptDUA: true) {{
      aliquot_id sample_id case_id aliquot_submitter_id sample_submitter_id case_submitter_id
      sample_type disease_type primary_site project_name externalReferences {{ reference_entity_location }}
    }} }}'''
    records = pdc_query(query)["biospecimenPerStudy"]
    return pd.DataFrame(records)


def matrix_frame(records: Any) -> pd.DataFrame:
    if not isinstance(records, list) or len(records) < 2:
        raise RuntimeError("PDC quantDataMatrix did not return a usable matrix.")
    return pd.DataFrame(records[1:], columns=records[0])


def fetch_matrix(study_id: str) -> pd.DataFrame:
    query = f'''query {{ quantDataMatrix(pdc_study_id: "{study_id}", data_type: "log2_ratio") }}'''
    data = pdc_query(query)
    return matrix_frame(data["quantDataMatrix"])


def extract_feature_table(matrix: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    first = matrix.columns[0]
    labels = matrix[first].astype(str)
    symbols = labels.str.split("|", regex=False).str[0].str.strip()
    rows: List[Dict[str, Any]] = []
    for gene in genes:
        matches = matrix.index[symbols == gene].tolist()
        if len(matches) == 1:
            rows.append({"gene": gene, "protein_feature_label": labels.loc[matches[0]], "feature_status": "UNIQUE"})
        elif len(matches) == 0:
            rows.append({"gene": gene, "protein_feature_label": "NOT_FOUND", "feature_status": "NOT_MEASURED_OR_NOT_MAPPED"})
        else:
            rows.append({"gene": gene, "protein_feature_label": "AMBIGUOUS_MULTIPLE_FEATURES", "feature_status": "AMBIGUOUS"})
    return pd.DataFrame(rows)


def make_protein_long(matrix: pd.DataFrame, biospecimen: pd.DataFrame, disease: str, genes: Sequence[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    first = matrix.columns[0]
    column_lookup = {str(column).split(":", 1)[1]: str(column) for column in matrix.columns[1:] if ":" in str(column)}
    allowed = biospecimen[(biospecimen["disease_type"] == disease) & (biospecimen["sample_type"].isin([TUMOR, NORMAL]))].copy()
    allowed["matrix_column"] = allowed["aliquot_submitter_id"].map(column_lookup)
    allowed["matrix_mapping_status"] = np.where(
        allowed["matrix_column"].notna(), "MAPPED_BY_EXACT_PDC_ALIQUOT_SUBMITTER_ID", "EXCLUDE_NOT_IN_LOG2_RATIO_MATRIX"
    )
    eligible = allowed[allowed["matrix_mapping_status"] == "MAPPED_BY_EXACT_PDC_ALIQUOT_SUBMITTER_ID"].copy()
    eligible["n_matrix_aliquots_per_case_type"] = eligible.groupby(["case_submitter_id", "sample_type"])["aliquot_id"].transform("size")
    unique = eligible[eligible["n_matrix_aliquots_per_case_type"] == 1].copy()
    tissue_sets = unique.groupby("case_submitter_id")["sample_type"].agg(lambda values: set(values))
    pair_cases = set(tissue_sets[tissue_sets.apply(lambda values: values == {TUMOR, NORMAL})].index)
    allowed["selection_status"] = "EXCLUDE"
    selected_ids = set(unique[unique["case_submitter_id"].isin(pair_cases)]["aliquot_id"])
    allowed.loc[allowed["aliquot_id"].isin(selected_ids), "selection_status"] = "SELECT_FOR_EXACT_PDC_CASE_PAIRED_PROTEIN"
    allowed["selection_reason"] = np.where(
        allowed["selection_status"] == "SELECT_FOR_EXACT_PDC_CASE_PAIRED_PROTEIN",
        "Exact aliquot ID in matrix; one aliquot per case and tissue type; both tissue types present.",
        "Not selected under exact PDC case/type uniqueness gate.",
    )
    selected = allowed[allowed["selection_status"] == "SELECT_FOR_EXACT_PDC_CASE_PAIRED_PROTEIN"].copy()
    feature_table = extract_feature_table(matrix, genes)
    value_rows: List[Dict[str, Any]] = []
    for _, feature in feature_table.iterrows():
        label = feature["protein_feature_label"]
        if feature["feature_status"] != "UNIQUE":
            continue
        original = matrix[matrix[first].astype(str) == label].iloc[0]
        for _, sample in selected.iterrows():
            value_rows.append(
                {
                    "gene": feature["gene"],
                    "protein_feature_label": label,
                    "case_submitter_id": sample["case_submitter_id"],
                    "sample_type": sample["sample_type"],
                    "aliquot_id": sample["aliquot_id"],
                    "sample_submitter_id": sample.get("sample_submitter_id"),
                    "value": pd.to_numeric(original.get(sample["matrix_column"], np.nan), errors="coerce"),
                }
            )
    value_columns = ["gene", "protein_feature_label", "case_submitter_id", "sample_type", "aliquot_id", "sample_submitter_id", "value"]
    return allowed, pd.DataFrame(value_rows, columns=value_columns)


def gdc_files_for_cases(case_ids: Sequence[str]) -> List[Dict[str, Any]]:
    if not case_ids:
        return []
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.case_id", "value": list(case_ids)}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["RNA-Seq"]}},
        ],
    }
    fields = (
        "file_id,file_name,file_size,md5sum,data_type,data_format,experimental_strategy,analysis.workflow_type,"
        "associated_entities.entity_id,associated_entities.entity_type,"
        "cases.case_id,cases.submitter_id,cases.primary_site,cases.diagnoses.primary_diagnosis,"
        "cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type,"
        "cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id"
    )
    response = requests.post(
        GDC_FILES_URL,
        json={"filters": json.dumps(filters), "format": "JSON", "fields": fields, "size": 10000},
        timeout=300,
    )
    response.raise_for_status()
    return response.json()["data"]["hits"]


def diagnosis_from_hit(hit: Dict[str, Any]) -> str:
    values: List[str] = []
    for case in hit.get("cases") or []:
        for diagnosis in case.get("diagnoses") or []:
            if diagnosis.get("primary_diagnosis"):
                values.append(str(diagnosis["primary_diagnosis"]))
    return "; ".join(sorted(set(values))) if values else "NA"


def nested_aliquot_sample_map(case: Dict[str, Any], associated_aliquot_id: str) -> List[Dict[str, Any]]:
    found: List[Dict[str, Any]] = []
    for sample in case.get("samples") or []:
        for portion in sample.get("portions") or []:
            for analyte in portion.get("analytes") or []:
                for aliquot in analyte.get("aliquots") or []:
                    if aliquot.get("aliquot_id") == associated_aliquot_id:
                        found.append(sample)
    return list({sample.get("sample_id", str(i)): sample for i, sample in enumerate(found)}.values())


def derive_rna_manifest(hits: Sequence[Dict[str, Any]]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for hit in hits:
        cases = hit.get("cases") or []
        if len(cases) != 1:
            continue
        case = cases[0]
        associated = [entity.get("entity_id") for entity in hit.get("associated_entities") or [] if entity.get("entity_type") == "aliquot"]
        samples: List[Dict[str, Any]] = []
        for aliquot_id in associated:
            samples.extend(nested_aliquot_sample_map(case, aliquot_id))
        unique_samples = {sample.get("sample_id", str(i)): sample for i, sample in enumerate(samples)}
        types = sorted({sample.get("sample_type") for sample in unique_samples.values() if sample.get("sample_type")})
        tissue_type = types[0] if len(types) == 1 and types[0] in {TUMOR, NORMAL} else "NA"
        mapping_status = "EXPLICIT_ALIQUOT_TO_SINGLE_TISSUE_TYPE" if tissue_type != "NA" else (
            "EXCLUDE_NO_EXPLICIT_SAMPLE_MAPPING" if not types else "EXCLUDE_AMBIGUOUS_TISSUE_TYPE"
        )
        rows.append(
            {
                "file_id": hit.get("file_id") or hit.get("id"),
                "file_name": hit.get("file_name"),
                "file_size": hit.get("file_size"),
                "md5sum": hit.get("md5sum"),
                "workflow_type": (hit.get("analysis") or {}).get("workflow_type"),
                "case_id": case.get("case_id"),
                "case_submitter_id": case.get("submitter_id"),
                "primary_diagnosis": diagnosis_from_hit(hit),
                "tissue_type": tissue_type,
                "mapping_status": mapping_status,
                "n_explicit_samples": len(unique_samples),
                "explicit_sample_submitter_ids": ";".join(sorted({str(sample.get("submitter_id")) for sample in unique_samples.values()})),
                "associated_aliquot_ids": ";".join(sorted(set(str(x) for x in associated if x))),
            }
        )
    manifest = pd.DataFrame(rows)
    if manifest.empty:
        return pd.DataFrame(columns=["file_id", "file_name", "file_size", "md5sum", "workflow_type", "case_id", "case_submitter_id", "primary_diagnosis", "tissue_type", "mapping_status", "n_explicit_samples", "explicit_sample_submitter_ids", "associated_aliquot_ids", "selection_status", "selection_reason"])
    manifest["selection_status"] = "EXCLUDE"
    eligible = manifest[(manifest["mapping_status"] == "EXPLICIT_ALIQUOT_TO_SINGLE_TISSUE_TYPE") & (manifest["workflow_type"] == "STAR - Counts")].copy()
    eligible["_n_case_type_files"] = eligible.groupby(["case_id", "tissue_type"])["file_id"].transform("size")
    unique = eligible[eligible["_n_case_type_files"] == 1]
    pair_sets = unique.groupby("case_id")["tissue_type"].agg(lambda values: set(values))
    pair_cases = set(pair_sets[pair_sets.apply(lambda values: values == {TUMOR, NORMAL})].index)
    chosen_ids = set(unique[unique["case_id"].isin(pair_cases)]["file_id"])
    manifest.loc[manifest["file_id"].isin(chosen_ids), "selection_status"] = "SELECT_FOR_CASE_LEVEL_PAIRED_RNA"
    manifest["selection_reason"] = np.where(
        manifest["selection_status"] == "SELECT_FOR_CASE_LEVEL_PAIRED_RNA",
        "One explicit STAR-count file per tissue type within a GDC case with both tissues.",
        "Not selected under explicit aliquot/tissue and case-level uniqueness gate.",
    )
    return manifest.sort_values(["case_submitter_id", "tissue_type", "file_id"]).reset_index(drop=True)


def download_gdc_file(record: pd.Series, destination: Path) -> Tuple[str, int]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected_md5 = str(record["md5sum"])
    if destination.exists() and md5(destination) == expected_md5:
        return "REUSED_VERIFIED", destination.stat().st_size
    if destination.exists():
        destination.unlink()
    for attempt in range(1, 4):
        temporary = destination.with_suffix(destination.suffix + ".part")
        if temporary.exists():
            temporary.unlink()
        try:
            response = requests.get(f"{GDC_DATA_URL}/{record['file_id']}", stream=True, timeout=300)
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
            if md5(temporary) != expected_md5:
                raise RuntimeError("MD5 mismatch after download")
            temporary.replace(destination)
            return f"DOWNLOADED_ATTEMPT_{attempt}", destination.stat().st_size
        except Exception:
            if temporary.exists():
                temporary.unlink()
            if attempt == 3:
                raise
            time.sleep(3 * attempt)
    raise AssertionError("unreachable")


def read_rna_candidate_values(path: Path, genes: Sequence[str]) -> Dict[str, float]:
    table = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    required = {"gene_name", "fpkm_uq_unstranded"}
    if not required.issubset(table.columns):
        raise RuntimeError(f"Unexpected GDC STAR-count columns in {path.name}: {list(table.columns)}")
    values = {gene: np.nan for gene in genes}
    for _, row in table[table["gene_name"].isin(genes)][["gene_name", "fpkm_uq_unstranded"]].iterrows():
        value = pd.to_numeric(row["fpkm_uq_unstranded"], errors="coerce")
        values[str(row["gene_name"])] = float(np.log2(value + 1.0)) if pd.notna(value) and value >= 0 else np.nan
    return values


def paired_effects(long_values: pd.DataFrame, genes: pd.DataFrame, prefix: str, value_scale: str) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for _, gene_row in genes.iterrows():
        gene = gene_row["gene"]
        subset = long_values[long_values["gene"] == gene].copy() if not long_values.empty else pd.DataFrame()
        base = {"gene": gene, "status": "NOT_MEASURED_OR_NOT_MAPPED", "value_scale": value_scale}
        if subset.empty:
            rows.append(base)
            continue
        pivot = subset.pivot_table(index="case_submitter_id", columns="sample_type", values="value", aggfunc="first")
        if not {TUMOR, NORMAL}.issubset(pivot.columns):
            base["status"] = "NO_COMPLETE_CASE_PAIRS"
            rows.append(base)
            continue
        paired = pivot[[TUMOR, NORMAL]].dropna()
        difference = paired[TUMOR] - paired[NORMAL]
        n = len(difference)
        if n < 5:
            base.update({"n_complete_pairs": n, "status": "INSUFFICIENT_COMPLETE_PAIRS_LT_5"})
            rows.append(base)
            continue
        mean_delta = float(difference.mean())
        sd_delta = float(difference.std(ddof=1))
        se = sd_delta / math.sqrt(n)
        tcritical = stats.t.ppf(0.975, df=n - 1)
        ttest = stats.ttest_rel(paired[TUMOR], paired[NORMAL], nan_policy="omit")
        try:
            wilcoxon_p = float(stats.wilcoxon(difference, alternative="two-sided", zero_method="wilcox").pvalue)
        except ValueError:
            wilcoxon_p = np.nan
        rows.append(
            {
                "gene": gene,
                "n_complete_pairs": n,
                "tumor_mean": float(paired[TUMOR].mean()),
                "solid_tissue_normal_mean": float(paired[NORMAL].mean()),
                "paired_mean_delta_tumor_minus_normal": mean_delta,
                "paired_delta_sd": sd_delta,
                "paired_delta_se": se,
                "ci_lower": mean_delta - tcritical * se,
                "ci_upper": mean_delta + tcritical * se,
                "paired_t_p_value": float(ttest.pvalue),
                "paired_wilcoxon_p_value": wilcoxon_p,
                "direction": "TUMOR_HIGHER" if mean_delta > 0 else "TUMOR_LOWER" if mean_delta < 0 else "NO_DIFFERENCE",
                "status": "ANALYSED",
                "value_scale": value_scale,
            }
        )
    result = pd.DataFrame(rows)
    if "paired_t_p_value" not in result.columns:
        result["paired_t_p_value"] = np.nan
    if "paired_wilcoxon_p_value" not in result.columns:
        result["paired_wilcoxon_p_value"] = np.nan
    result = result.add_prefix(f"{prefix}_").rename(columns={f"{prefix}_gene": "gene"})
    # Add FDR after prefixing so downstream columns are exactly
    # rna_paired_t_fdr / protein_paired_t_fdr (not double-prefixed).
    result[f"{prefix}_paired_t_fdr"] = bh_adjust(result[f"{prefix}_paired_t_p_value"].tolist())
    result[f"{prefix}_paired_wilcoxon_fdr"] = bh_adjust(result[f"{prefix}_paired_wilcoxon_p_value"].tolist())
    return result


def fmt(value: Any, digits: int = 3) -> str:
    try:
        number = float(value)
        if not np.isfinite(number):
            return "NA"
        if abs(number) < 0.001 and number != 0:
            return f"{number:.2e}"
        return f"{number:.{digits}f}"
    except (TypeError, ValueError):
        return text_or_na(value)


def integrated_effect(effect: pd.DataFrame) -> pd.DataFrame:
    effect = effect.copy()
    effect["integrated_direction"] = "NOT_EVALUABLE"
    rna_delta = pd.to_numeric(effect.get("rna_paired_mean_delta_tumor_minus_normal"), errors="coerce")
    protein_delta = pd.to_numeric(effect.get("protein_paired_mean_delta_tumor_minus_normal"), errors="coerce")
    rna_fdr = pd.to_numeric(effect.get("rna_paired_t_fdr"), errors="coerce")
    protein_fdr = pd.to_numeric(effect.get("protein_paired_t_fdr"), errors="coerce")
    same_up = (rna_delta > 0) & (protein_delta > 0)
    same_down = (rna_delta < 0) & (protein_delta < 0)
    effect.loc[same_up, "integrated_direction"] = "RNA_PROTEIN_TUMOR_HIGHER"
    effect.loc[same_down, "integrated_direction"] = "RNA_PROTEIN_TUMOR_LOWER"
    effect.loc[(rna_delta * protein_delta < 0), "integrated_direction"] = "RNA_PROTEIN_DISCORDANT"
    effect["cptac_support_interpretation"] = "UNRESOLVED_OR_NOT_EVALUABLE"
    effect.loc[(rna_fdr <= 0.05) & (protein_fdr <= 0.05) & same_up, "cptac_support_interpretation"] = "SUPPORTIVE_SAME_DIRECTION_TUMOR_HIGHER_NOT_CAUSAL"
    effect.loc[(rna_fdr <= 0.05) & (protein_fdr <= 0.05) & same_down, "cptac_support_interpretation"] = "SUPPORTIVE_SAME_DIRECTION_TUMOR_LOWER_NOT_CAUSAL"
    effect.loc[(rna_delta * protein_delta < 0), "cptac_support_interpretation"] = "CROSS_LAYER_DISCORDANT_DO_NOT_INTERPRET_AS_ABSENCE"
    effect.loc[effect["integrated_direction"] == "NOT_EVALUABLE", "cptac_support_interpretation"] = "ONE_OR_BOTH_OMIC_LAYERS_NOT_EVALUABLE"
    return effect


def file_inventory(paths: Iterable[Path], base: Path) -> pd.DataFrame:
    records = []
    for path in sorted({path for path in paths if path.exists()}):
        records.append({"relative_path": str(path.relative_to(base)).replace("\\", "/") if base in path.parents or path == base else str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    return pd.DataFrame(records)


def build_report(path: Path, effects: pd.DataFrame, coverage: pd.DataFrame, source_audit: Dict[str, Any], gene_count: Optional[int] = None) -> None:
    direct_genes = effects.loc[effects.get("panel_group", pd.Series(index=effects.index, dtype=object)).isin(["direct_inosine_core", "exploratory_inosine_context"]), "gene"].drop_duplicates().tolist() if "panel_group" in effects.columns else []
    direct_note = "本轮新增/探索性直接 inosine 候选：" + "、".join(direct_genes) + "。" if direct_genes else ""
    rows: List[str] = []
    for _, item in effects.iterrows():
        direction = item.get("integrated_direction", "NA")
        colour = "#b42318" if direction == "RNA_PROTEIN_TUMOR_HIGHER" else "#1d5f9c" if direction == "RNA_PROTEIN_TUMOR_LOWER" else "#6b7280"
        rows.append(
            "<tr>"
            f"<td><b>{html.escape(text_or_na(item.get('cancer')))}</b><br>{html.escape(text_or_na(item.get('gene')))}<br><small>{html.escape(text_or_na(item.get('tier')))} / {html.escape(text_or_na(item.get('panel_group')))}</small></td>"
            f"<td>{fmt(item.get('rna_n_complete_pairs'))}<br>{fmt(item.get('rna_paired_mean_delta_tumor_minus_normal'))}<br>FDR {fmt(item.get('rna_paired_t_fdr'))}</td>"
            f"<td>{fmt(item.get('protein_n_complete_pairs'))}<br>{fmt(item.get('protein_paired_mean_delta_tumor_minus_normal'))}<br>FDR {fmt(item.get('protein_paired_t_fdr'))}</td>"
            f"<td style='color:{colour};font-weight:700'>{html.escape(direction)}</td>"
            f"<td>{html.escape(text_or_na(item.get('cptac_support_interpretation')))}</td>"
            "</tr>"
        )
    coverage_rows = "".join(
        f"<tr><td>{html.escape(str(row.cancer))}</td><td>{html.escape(str(row.pdc_study_id))}</td><td>{html.escape(str(row.pdc_tissue_records))}</td><td>{html.escape(str(row.pdc_exact_case_pairs))}</td><td>{html.escape(str(row.gdc_selected_files))}</td><td>{html.escape(str(row.gdc_exact_case_pairs))}</td><td>{html.escape(str(row.camp_inosine_direction))}</td></tr>"
        for row in coverage.itertuples()
    )
    audit_lines = "".join(f"<li><b>{html.escape(str(k))}</b>: {html.escape(str(v))}</li>" for k, v in source_audit.items())
    gene_count = gene_count if gene_count is not None else int(effects["gene"].nunique())
    document = f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><title>Inosine pan-cancer CPTAC</title>
<style>body{{font-family:Arial,'Microsoft YaHei',sans-serif;line-height:1.55;color:#1f2937;margin:34px;max-width:1450px}}h1{{color:#064e3b}}h2{{color:#115e59;margin-top:28px}}.badge{{background:#fef3c7;color:#92400e;padding:5px 9px;border-radius:999px;font-weight:700}}.note{{background:#eff6ff;border-left:4px solid #2563eb;padding:12px 16px}}table{{border-collapse:collapse;width:100%;font-size:13px}}th{{background:#ecfdf5;color:#065f46}}td,th{{border:1px solid #d1d5db;padding:8px;vertical-align:top;text-align:left}}small{{color:#6b7280}}</style></head><body>
<h1>Inosine（肌苷）跨癌种 CPTAC RNA/蛋白支持性分析</h1>
<p><span class='badge'>EXPLORATORY_CPTAC_SUPPORT</span> 新分析癌种：BRCA、COAD、ccRCC；PDAC 保留既有独立结果作对照，不在本次重算。PRAD 在本轮未纳入，原因见 cancer inclusion audit。</p>
<div class='note'><b>问题边界：</b>冻结 CAMP 结果中的 inosine（KEGG C00294/HMDB00195）是组织肿瘤–正常差异背景；本报告检验锁定的核苷代谢相关基因在 CPTAC RNA/蛋白层面是否有同向差异。CPTAC 样本与 CAMP 不是逐患者配对，因此这不是 inosine 的外部代谢物验证，也不能建立基因–代谢物直接相关或因果通量。</div>
<h2>癌种与配对覆盖</h2><table><thead><tr><th>癌种</th><th>PDC</th><th>疾病记录</th><th>蛋白精确病例对</th><th>RNA入选文件</th><th>RNA精确病例对</th><th>CAMP inosine方向</th></tr></thead><tbody>{coverage_rows}</tbody></table>
<h2>候选基因层间结果</h2><p>效应均为肿瘤 − Solid Tissue Normal；RNA 为 log2(FPKM-UQ+1)，蛋白为 PDC 提供的 log2_ratio。FDR 在每个癌种、每个组学层内对 {gene_count} 个预设基因分别进行 BH 校正。表中依次显示配对数、均值效应、t 检验 FDR。{html.escape(direct_note)}</p>
<table><thead><tr><th>癌种/基因</th><th>RNA</th><th>蛋白</th><th>综合方向</th><th>解释状态</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<h2>审计摘要</h2><ul>{audit_lines}</ul>
<h2>解释限制</h2><ul><li>BRCA、COAD 若 RNA 没有明确病例内肿瘤–正常配对，RNA 状态保留为不可评估；不将缺失写成阴性。</li><li>蛋白结果只表示组织层面的表达/定量差异；不能替代酶活性、底物通量、细胞来源或纯度校正。</li><li>ccRCC 采用 PDC000127 中明确标注为 Clear Cell Renal Cell Carcinoma 的记录，未把 Non-Clear Cell 或 Other 推断为 ccRCC。</li><li>跨癌种方向比较是描述性的；未把不同癌种的原始表达值直接合并，也未进行跨癌种因果或机制排序。</li></ul>
<p>机器可读结果位于 <code>tables/</code>；含既有 PDAC 对照的方向矩阵为 <code>inosine_cptac_direction_matrix_including_pdac.tsv</code>。源矩阵只保留在 server165 的候选数据目录，具体文件与校验值见 <code>checksums/</code>。</p></body></html>"""
    path.write_text(document, encoding="utf-8")


def reuse_existing_pdac(root: Path, genes: pd.DataFrame, camp_direction: str) -> pd.DataFrame:
    """Load the already completed PDAC CPTAC table without re-running it."""
    path = root / "results" / "PDAC_mechanism_exploration_v0.1" / "tables" / "cptac_candidate_rna_protein_effects.tsv"
    if not path.exists():
        return pd.DataFrame()
    old = pd.read_csv(path, sep="\t")
    old = old[old["gene"].isin(genes["gene"].tolist())].copy()
    rows: List[Dict[str, Any]] = []
    for _, row in old.iterrows():
        rows.append(
            {
                "cancer": "PDAC",
                "gene": row.get("gene"),
                "tier": row.get("tier"),
                "axis": row.get("axis"),
                "relation": row.get("relation"),
                "rna_n_complete_pairs": row.get("rna_n_complete_pairs"),
                "rna_paired_mean_delta_tumor_minus_normal": row.get("rna_paired_mean_delta_tumor_minus_normal"),
                "rna_paired_t_fdr": row.get("rna_paired_t_fdr"),
                "rna_direction": row.get("rna_direction"),
                "rna_status": row.get("rna_status"),
                "protein_n_complete_pairs": row.get("protein_n_complete_pairs"),
                "protein_paired_mean_delta_tumor_minus_normal": row.get("protein_paired_mean_delta_tumor_minus_normal"),
                "protein_paired_t_fdr": row.get("protein_paired_t_fdr"),
                "protein_direction": row.get("protein_direction"),
                "protein_status": row.get("protein_status"),
                "integrated_direction": row.get("integrated_direction"),
                "cptac_support_interpretation": row.get("cptac_support_interpretation"),
                "pdc_study_id": "PDC000270",
                "pdc_study_name": "Existing PDAC CPTAC support run (reused; no re-computation)",
                "camp_inosine_direction": camp_direction,
                "evidence_stage": "REUSED_EXISTING_PDAC_CPTAC_SUPPORT",
                "causal_claim_allowed": False,
                "result_origin": "REUSED_EXISTING_PDAC_RUN",
            }
        )
    # Newly added genes were not in the frozen PDAC panel. Keep explicit rows
    # in the comparison table rather than silently treating them as absent or
    # negative in PDAC.
    old_genes = set(old["gene"].tolist())
    for _, gene_row in genes[~genes["gene"].isin(old_genes)].iterrows():
        rows.append(
            {
                "cancer": "PDAC",
                "gene": gene_row["gene"],
                "tier": gene_row.get("tier"),
                "axis": gene_row.get("axis"),
                "relation": gene_row.get("relation"),
                "rna_status": "NOT_EVALUABLE_NOT_IN_REUSED_PDAC_PANEL",
                "protein_status": "NOT_EVALUABLE_NOT_IN_REUSED_PDAC_PANEL",
                "integrated_direction": "NOT_EVALUABLE",
                "cptac_support_interpretation": "NOT_RECOMPUTED_NEW_GENE_PANEL_EXTENSION",
                "pdc_study_id": "PDC000270",
                "pdc_study_name": "Existing PDAC CPTAC support run (reused; no re-computation)",
                "camp_inosine_direction": camp_direction,
                "evidence_stage": "REUSED_EXISTING_PDAC_CPTAC_SUPPORT",
                "causal_claim_allowed": False,
                "result_origin": "NOT_IN_EXISTING_PDAC_PANEL",
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty and "panel_group" in genes.columns:
        panel_map = genes.set_index("gene")["panel_group"].to_dict()
        result["panel_group"] = result["gene"].map(panel_map)
    return result


def analyse_cancer(cancer: str, spec: Dict[str, Any], genes: pd.DataFrame, source_root: Path, output: Path, refresh: bool, max_download_bytes: int) -> Tuple[pd.DataFrame, Dict[str, Any], Dict[str, Any]]:
    study_id = spec["pdc_study_id"]
    disease = spec["pdc_disease_type"]
    cancer_source = source_root / cancer
    (cancer_source / "pdc").mkdir(parents=True, exist_ok=True)
    (cancer_source / "gdc_rna" / "files").mkdir(parents=True, exist_ok=True)
    biospecimen_path = cancer_source / "pdc" / f"{study_id.lower()}_biospecimen.json"
    matrix_path = cancer_source / "pdc" / f"{study_id.lower()}_log2_ratio_matrix.json"
    if refresh or not biospecimen_path.exists():
        write_json(biospecimen_path, fetch_biospecimen(study_id).to_dict(orient="records"))
    if refresh or not matrix_path.exists():
        matrix = fetch_matrix(study_id)
        write_json(matrix_path, [matrix.columns.tolist()] + matrix.astype(object).where(pd.notna(matrix), None).values.tolist())
    biospecimen = pd.DataFrame(json.loads(biospecimen_path.read_text(encoding="utf-8")))
    matrix = matrix_frame(json.loads(matrix_path.read_text(encoding="utf-8")))
    matrix_tsv = cancer_source / "pdc" / f"{study_id.lower()}_log2_ratio_matrix.tsv.gz"
    if not matrix_tsv.exists():
        matrix.to_csv(matrix_tsv, sep="\t", index=False, compression="gzip")
    protein_audit, protein_long = make_protein_long(matrix, biospecimen, disease, genes["gene"].tolist())
    pdac_like = biospecimen[(biospecimen["disease_type"] == disease) & (biospecimen["sample_type"].isin([TUMOR, NORMAL]))]
    case_ids = sorted({case_id for _, row in pdac_like.iterrows() if (case_id := gdc_case_id_from_refs(row))})
    gdc_manifest_path = cancer_source / "gdc_rna" / "gdc_star_counts_manifest.tsv"
    if gdc_manifest_path.exists() and not refresh:
        gdc_manifest = pd.read_csv(gdc_manifest_path, sep="\t", dtype=str)
    else:
        gdc_manifest = derive_rna_manifest(gdc_files_for_cases(case_ids))
        gdc_manifest.to_csv(gdc_manifest_path, sep="\t", index=False)
    selected = gdc_manifest[gdc_manifest.get("selection_status", pd.Series(dtype=str)) == "SELECT_FOR_CASE_LEVEL_PAIRED_RNA"].copy()
    if not selected.empty:
        selected["file_size_num"] = pd.to_numeric(selected["file_size"], errors="coerce").fillna(0)
    projected_bytes = int(selected["file_size_num"].sum()) if not selected.empty else 0
    if projected_bytes > max_download_bytes:
        raise RuntimeError(f"{cancer} selected RNA download is {projected_bytes} bytes, above configured limit {max_download_bytes}.")
    rna_rows: List[Dict[str, Any]] = []
    download_rows: List[Dict[str, Any]] = []
    for index, record in selected.reset_index(drop=True).iterrows():
        destination = cancer_source / "gdc_rna" / "files" / str(record["file_name"])
        state, byte_count = download_gdc_file(record, destination)
        values = read_rna_candidate_values(destination, genes["gene"].tolist())
        download_rows.append({"cancer": cancer, "file_id": record["file_id"], "file_name": record["file_name"], "case_submitter_id": record["case_submitter_id"], "tissue_type": record["tissue_type"], "download_status": state, "bytes": byte_count, "md5_expected": record["md5sum"], "sha256": sha256(destination)})
        for gene, value in values.items():
            rna_rows.append({"cancer": cancer, "gene": gene, "case_submitter_id": record["case_submitter_id"], "case_id": record["case_id"], "sample_type": record["tissue_type"], "file_id": record["file_id"], "file_name": record["file_name"], "value": value, "value_scale": "log2(fpkm_uq_unstranded+1)"})
        if (index + 1) % 20 == 0 or index + 1 == len(selected):
            print(f"{cancer} RNA files processed: {index + 1}/{len(selected)}", flush=True)
    rna_long = pd.DataFrame(rna_rows, columns=["cancer", "gene", "case_submitter_id", "case_id", "sample_type", "file_id", "file_name", "value", "value_scale"])
    protein_long["cancer"] = cancer
    rna_long.to_csv(output / "tables" / f"{cancer}_rna_candidate_values.tsv", sep="\t", index=False)
    protein_long.to_csv(output / "tables" / f"{cancer}_protein_candidate_values.tsv", sep="\t", index=False)
    protein_audit.to_csv(output / "tables" / f"{cancer}_protein_pair_audit.tsv", sep="\t", index=False)
    gdc_manifest.to_csv(output / "tables" / f"{cancer}_rna_file_manifest_audit.tsv", sep="\t", index=False)
    pd.DataFrame(download_rows).to_csv(output / "tables" / f"{cancer}_rna_file_download_audit.tsv", sep="\t", index=False)
    rna_effect = paired_effects(rna_long, genes, "rna", "log2(FPKM-UQ+1)")
    protein_effect = paired_effects(protein_long, genes, "protein", "PDC supplied log2_ratio")
    effect = integrated_effect(genes.merge(rna_effect, on="gene", how="left").merge(protein_effect, on="gene", how="left"))
    effect.insert(0, "cancer", cancer)
    effect["pdc_study_id"] = study_id
    effect["pdc_study_name"] = spec.get("pdc_study_name", "")
    effect["camp_inosine_direction"] = spec["camp_inosine_direction"]
    effect["evidence_stage"] = "EXPLORATORY_CPTAC_SUPPORT"
    effect["causal_claim_allowed"] = False
    coverage = {
        "cancer": cancer,
        "pdc_study_id": study_id,
        "pdc_study_name": spec.get("pdc_study_name", ""),
        "pdc_disease_type": disease,
        "pdc_tissue_records": int(len(pdac_like)),
        "pdc_cases": int(pdac_like["case_submitter_id"].nunique()),
        "pdc_exact_case_pairs": int(protein_audit.loc[protein_audit["selection_status"] == "SELECT_FOR_EXACT_PDC_CASE_PAIRED_PROTEIN", "case_submitter_id"].nunique()),
        "pdc_matrix_columns": int(max(0, matrix.shape[1] - 1)),
        "gdc_case_ids_from_pdc": int(len(case_ids)),
        "gdc_files_queried": int(len(gdc_manifest)),
        "gdc_selected_files": int(len(selected)),
        "gdc_exact_case_pairs": int(selected["case_submitter_id"].nunique()) if not selected.empty else 0,
        "gdc_projected_download_bytes": projected_bytes,
        "camp_inosine_direction": spec["camp_inosine_direction"],
    }
    audit = {f"{cancer}_PDC": f"{study_id}; {spec.get('pdc_study_name', '')}; disease_type={disease}", f"{cancer}_PDC_tissue_records": coverage["pdc_tissue_records"], f"{cancer}_PDC_exact_case_pairs": coverage["pdc_exact_case_pairs"], f"{cancer}_GDC_files_queried": coverage["gdc_files_queried"], f"{cancer}_GDC_selected_files": coverage["gdc_selected_files"], f"{cancer}_GDC_exact_case_pairs": coverage["gdc_exact_case_pairs"], f"{cancer}_RNA_projected_download": f"{projected_bytes / 1024**2:.1f} MiB", f"{cancer}_camp_inosine_direction": spec["camp_inosine_direction"]}
    return effect, coverage, audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--refresh-pdc", action="store_true")
    parser.add_argument("--max-rna-download-gb", type=float, default=2.0)
    parser.add_argument("--config-dir", default="config/inosine_pan_cancer_cptac_v0.1", help="Configuration directory relative to project root.")
    parser.add_argument("--analysis-name", default=ANALYSIS_NAME, help="Result directory/name; defaults to v0.1.")
    parser.add_argument("--source-dir-name", default="cptac_inosine_pan_cancer_v0.1", help="Remote source directory name under data/candidates.")
    args = parser.parse_args()
    root = args.project_root.resolve()
    analysis_name = args.analysis_name
    config_path = root / args.config_dir / "analysis_parameters.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    output = root / "results" / analysis_name
    source_root = root / "data" / "candidates" / args.source_dir_name
    for directory in [output / "parameters", output / "tables", output / "reports", output / "figures", output / "logs", output / "checksums", source_root]:
        directory.mkdir(parents=True, exist_ok=True)
    (output / "parameters" / config_path.name).parent.mkdir(parents=True, exist_ok=True)
    (output / "parameters" / config_path.name).write_text(config_path.read_text(encoding="utf-8"), encoding="utf-8")
    genes = pd.DataFrame([{"gene": gene, **details} for gene, details in config["genes"].items()])
    for cancer, spec in config["cancers"].items():
        spec["camp_inosine_direction"] = config["metabolite"]["camp_direction_by_cancer"][cancer]
    started = now_utc()
    effects: List[pd.DataFrame] = []
    coverage_rows: List[Dict[str, Any]] = []
    source_audit: Dict[str, Any] = {"analysis": analysis_name, "candidate_gene_count": len(genes), "PDC_join_rule": "exact aliquot_submitter_id", "GDC_pairing_rule": "explicit aliquot to single sample type then case-level pair", "missing_data_rule": "unevaluable_not_negative", "CAMP_context": "inosine C00294/HMDB00195 direction from frozen CAMP Phase1 table; no re-computation", "purity_stroma_adjustment": "NOT_PERFORMED", "metabolite_gene_direct_correlation": "PROHIBITED_DIFFERENT_PATIENTS"}
    for cancer, spec in config["cancers"].items():
        effect, coverage, audit = analyse_cancer(cancer, spec, genes, source_root, output, args.refresh_pdc, int(args.max_rna_download_gb * 1024**3))
        effects.append(effect)
        coverage_rows.append(coverage)
        source_audit.update(audit)
    effect_table = pd.concat(effects, ignore_index=True)
    coverage = pd.DataFrame(coverage_rows)
    effect_table.to_csv(output / "tables" / "inosine_cptac_effects.tsv", sep="\t", index=False)
    # Separate view of the newly added direct inosine candidates. The full
    # table retains the original context genes for comparability.
    direct_extension = effect_table[effect_table["panel_group"].isin(["direct_inosine_core", "exploratory_inosine_context"])].copy()
    direct_extension.to_csv(output / "tables" / "inosine_direct_extension_effects.tsv", sep="\t", index=False)
    genes[["gene", "tier", "panel_group", "axis", "relation"]].to_csv(output / "tables" / "inosine_direct_gene_panel.tsv", sep="\t", index=False)
    coverage.to_csv(output / "tables" / "inosine_cptac_coverage.tsv", sep="\t", index=False)
    # If PDAC is explicitly included in the current configuration, analyse it
    # with the current locked gene panel. Otherwise retain the historical
    # behaviour of reusing the prior PDAC result without re-computation.
    pdac_reuse = pd.DataFrame() if "PDAC" in config.get("cancers", {}) else reuse_existing_pdac(root, genes, config["metabolite"]["camp_direction_by_cancer"]["PDAC"])
    comparison = pd.concat([effect_table, pdac_reuse], ignore_index=True, sort=False) if not pdac_reuse.empty else effect_table.copy()
    comparison.to_csv(output / "tables" / "inosine_cptac_effects_including_pdac_reused.tsv", sep="\t", index=False)
    effect_table.pivot_table(index="gene", columns="cancer", values="integrated_direction", aggfunc="first").reset_index().to_csv(output / "tables" / "inosine_cptac_direction_matrix.tsv", sep="\t", index=False)
    comparison.pivot_table(index="gene", columns="cancer", values="integrated_direction", aggfunc="first").reset_index().to_csv(output / "tables" / "inosine_cptac_direction_matrix_including_pdac.tsv", sep="\t", index=False)
    comparison.pivot_table(index="gene", columns="cancer", values="protein_paired_mean_delta_tumor_minus_normal", aggfunc="first").reset_index().to_csv(output / "tables" / "inosine_cptac_protein_effect_matrix_including_pdac.tsv", sep="\t", index=False)
    inclusion = pd.DataFrame(
        [
            {"cancer": "BRCA", "status": "ANALYSED_NEW", "basis": "PDC000120 exact disease label and protein matrix"},
            {"cancer": "COAD", "status": "ANALYSED_NEW", "basis": "PDC000116 exact disease label and protein matrix"},
            {"cancer": "ccRCC", "status": "ANALYSED_NEW", "basis": "PDC000127 exact Clear Cell Renal Cell Carcinoma label and protein matrix"},
            {"cancer": "PDAC", "status": "REUSED_EXISTING_RESULT", "basis": "Existing PDAC_mechanism_exploration_v0.1; no re-computation"},
            {"cancer": "PRAD", "status": "NOT_INCLUDED_THIS_RUN", "basis": "No confirmed PDC protein study with the same strict paired workflow was selected"},
        ]
    )
    inclusion.to_csv(output / "tables" / "inosine_cptac_cancer_inclusion_audit.tsv", sep="\t", index=False)
    write_json(output / "tables" / "inosine_cptac_source_and_pairing_audit.json", source_audit)
    build_report(output / "reports" / "inosine_pan_cancer_cptac_report_zh.html", effect_table, coverage, source_audit, len(genes))
    source_manifest = []
    for cancer, spec in config["cancers"].items():
        source_manifest.append({"cancer": cancer, "source": "PDC GraphQL", "url": PDC_URL, "study_or_endpoint": spec["pdc_study_id"], "resource": "biospecimenPerStudy; quantDataMatrix(log2_ratio)", "retrieved_at_utc": now_utc()})
        source_manifest.append({"cancer": cancer, "source": "GDC API", "url": GDC_FILES_URL, "study_or_endpoint": "case IDs from PDC externalReferences", "resource": "Gene Expression Quantification; STAR - Counts", "retrieved_at_utc": now_utc()})
    pd.DataFrame(source_manifest).to_csv(output / "tables" / "inosine_cptac_source_manifest.tsv", sep="\t", index=False)
    analysis_files = list((output / "parameters").glob("*")) + list((output / "tables").glob("*")) + list((output / "reports").glob("*"))
    file_inventory(analysis_files, root).to_csv(output / "checksums" / "analysis_file_checksums.tsv", sep="\t", index=False)
    # Include every remote source artifact, including downloaded STAR-count
    # files, so the audit package can verify all source bytes without copying
    # those large files into the local workspace.
    source_files = [path for path in source_root.rglob("*") if path.is_file()]
    file_inventory(source_files, root).to_csv(output / "checksums" / "source_file_checksums.tsv", sep="\t", index=False)
    log_lines = [f"RUN_STARTED_UTC\t{started}", f"RUN_FINISHED_UTC\t{now_utc()}", f"ANALYSIS\t{analysis_name}", f"CANCERS\t{','.join(config['cancers'])}", f"GENES\t{len(genes)}", f"PDC_EXACT_CASE_PAIRS_TOTAL\t{int(coverage['pdc_exact_case_pairs'].sum())}", f"GDC_EXACT_CASE_PAIRS_TOTAL\t{int(coverage['gdc_exact_case_pairs'].sum())}", "STATUS\tPASS_CPTAC_MULTICANCER_EXPLORATORY_SUPPORT"]
    (output / "logs" / "inosine_cptac_run.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS_CPTAC_MULTICANCER_EXPLORATORY_SUPPORT", "output": str(output), "coverage": coverage.to_dict(orient="records")}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
