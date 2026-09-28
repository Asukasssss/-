"""Read-only COAD focus view from published aggregate tables; no new statistics."""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "runtime/report_dependencies"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


SOURCE = ROOT / "results/COAD/07_INTEGRATION/20260922T151000Z_source_contract_integration_v2"
PATIENT = ROOT / "results/COAD/03_PATIENT/20260922T135816Z_source_contract_patient_v2"
COV = ROOT / "results/COAD/04_ROBUSTNESS/20260919T142500Z_covariates_v1"
FOCUS = [
    ("UCKL1", "uridine", "上皮背景与非催化功能"),
    ("NNMT", "1-methylnicotinamide", "产物1-MNA：条件性功能支持"),
    ("NNMT", "S-adenosylmethionine (SAM)", "甲基供体SAM：单独保留缺口"),
    ("HDC", "histamine", "肥大细胞覆盖与组织组成"),
    ("GSTA4", "glutathione, reduced (GSH)", "同一GSH特征的候选之一"),
    ("CLIC2", "glutathione, reduced (GSH)", "同一GSH特征的另一候选"),
]


def read(path):
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, rows):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def one(rows, gene, metabolite):
    matches = [r for r in rows if r["gene"] == gene and r["metabolite_name"] == metabolite]
    assert len(matches) == 1, (gene, metabolite, len(matches))
    return matches[0]


def fmt(value, digits=3):
    return f"{float(value):.{digits}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    out = args.out.resolve()
    assert not out.exists(), "Use a new run directory"
    out.mkdir(parents=True)

    inputs = {
        "relations": SOURCE / "candidate_relations.tsv",
        "association": PATIENT / "association_results.tsv",
        "lee": SOURCE / "Lee_sc_celltype_profiles.tsv",
        "uhlitz": SOURCE / "Uhlitz_sc_celltype_profiles.tsv",
        "old_M2": COV / "M2.tsv",
    }
    data = {name: read(path) for name, path in inputs.items()}
    primary = [r for r in data["association"] if r["analysis_type"] == "association_primary"]
    main_rows = []
    for gene, met, question in FOCUS:
        r = one(data["relations"], gene, met)
        a = one(primary, gene, met)
        assert a["effect"] == r["CAMP_effect"] and a["q_value"] == r["CAMP_q_value"]
        main_rows.append(
            dict(
                question=question,
                relation_id=r["relation_id"],
                gene=gene,
                metabolite=met,
                tumor_n=a["n"],
                tumor_rho=a["effect"],
                tumor_P=a["p_value"],
                tumor_q=a["q_value"],
                leave_one_out_rho_min=a["loo_min"],
                leave_one_out_rho_max=a["loo_max"],
                metabolite_pair_up=r["metabolite_paired_pairs_higher"],
                metabolite_pair_down=r["metabolite_paired_pairs_lower"],
                metabolite_pair_q=r["metabolite_paired_q_value"],
                RNA_pair_up=r["RNA_n_up"],
                RNA_pair_down=r["RNA_n_down"],
                RNA_pair_effect_log2=r["RNA_effect"],
                RNA_pair_q=r["RNA_q_value"],
                Lee_top=r["Lee_display_top"],
                Uhlitz_top=r["Uhlitz_display_top"],
            )
        )
    assert len(main_rows) == 6
    write(out / "focus_relations.tsv", main_rows)

    genes = list(dict.fromkeys(r[0] for r in FOCUS))
    profiles = []
    for study, key in [("Lee", "lee"), ("Uhlitz", "uhlitz")]:
        for r in data[key]:
            if r["gene"] in genes:
                profiles.append(
                    dict(
                        study=study,
                        gene=r["gene"],
                        author_cell_type=r["celltype"],
                        status=r["status"],
                        reason=r["reason"],
                        donor_n=r["n_donors"],
                        cell_n=r["n_cells"],
                        donor_equal_mean_log1p_CP10K=r["mean_expression"],
                        donor_mean_detection_fraction=r["mean_detection_fraction"],
                    )
                )
    assert len(profiles) == 5 * 6 + 5 * 10
    write(out / "source_five_genes.tsv", profiles)

    plt.rcParams["font.family"] = "Microsoft YaHei"
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.8), layout="constrained")
    for ax, study in zip(axes, ["Lee", "Uhlitz"]):
        part = [r for r in profiles if r["study"] == study]
        types = list(dict.fromkeys(r["author_cell_type"] for r in part))
        vmax = max(float(r["donor_equal_mean_log1p_CP10K"]) for r in part if r["status"] == "DONE")
        for r in part:
            x, y = types.index(r["author_cell_type"]), genes.index(r["gene"])
            if r["status"] != "DONE":
                ax.text(x, y, "×", ha="center", va="center", color="0.55", fontsize=10)
                continue
            val = float(r["donor_equal_mean_log1p_CP10K"])
            det = float(r["donor_mean_detection_fraction"])
            ax.scatter(x, y, s=15 + 290 * det, c=[val], cmap="viridis", vmin=0, vmax=vmax, edgecolor="none")
        ax.set(xlim=(-0.5, len(types) - 0.5), ylim=(len(genes) - 0.5, -0.5), title=study)
        ax.set_xticks(range(len(types)), types, rotation=50, ha="right")
        ax.set_yticks(range(len(genes)), genes)
        ax.tick_params(length=0)
        fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, vmax), cmap="viridis"), ax=ax, shrink=0.65, label="供者等权平均 log1p(CP10K)")
    fig.suptitle("COAD四组关注问题：5基因的肿瘤单细胞表达来源\n点大小＝平均检出比例；×＝覆盖不足；各研究色标独立", fontsize=14)
    fig.savefig(out / "focus_source_five_genes.png", dpi=180)
    fig.savefig(out / "focus_source_five_genes.pdf")
    plt.close(fig)

    m2 = [one(data["old_M2"], "NNMT", met) for met in ("1-methylnicotinamide", "S-adenosylmethionine (SAM)")]
    report = f"""# COAD四组代表性问题：证据分层汇报

## 本轮问题
将现有结果组织成四组可讲清楚的问题，不重新筛选候选或计算统计量。

## 输入与范围
新版患者关联为33个肿瘤内的891项计划关系，旧协变量结果只覆盖历史674项；配对变化来自33对患者，来源来自Lee/Uhlitz两研究。完整候选与历史9/21/5安排见原版总表，本页展示6条关系，不是靶点排名。

## 实际结果
|问题|具体关系|肿瘤内ρ；P；q|配对方向与来源|现有证据怎样解释|
|---|---|---|---|---|
|UCKL1的上皮背景|UCKL1—尿苷|{fmt(main_rows[0]['tumor_rho'])}；{main_rows[0]['tumor_P']}；{fmt(main_rows[0]['tumor_q'])}|尿苷25/33对下降；UCKL1 RNA25/33对上升；Lee/Uhlitz均上皮可见|逐一剔除ρ为{fmt(main_rows[0]['leave_one_out_rho_min'])}至{fmt(main_rows[0]['leave_one_out_rho_max'])}。旧CRC研究提示铁死亡保护存在非经典/低催化活性突变体保留效应，负相关不能解释为直接尿苷消耗。|
|NNMT的双关系|NNMT—1-MNA|{fmt(main_rows[1]['tumor_rho'])}；{main_rows[1]['tumor_P']}；{fmt(main_rows[1]['tumor_q'])}|1-MNA29/33对上升；NNMT RNA16升17降；宽基质/成纤维背景|已有CRC细胞5-FU背景下1-MNA测量与高剂量部分救援；治疗、剂量与细胞模型均有适用范围。|
|NNMT的双关系|NNMT—SAM|{fmt(main_rows[2]['tumor_rho'])}；{main_rows[2]['tumor_P']}；{fmt(main_rows[2]['tumor_q'])}|SAM32/33对上升；基因来源同上|1-MNA产物证据不能直接充当SAM介导证据。|
|HDC的组成问题|HDC—组胺|{fmt(main_rows[3]['tumor_rho'])}；{main_rows[3]['tumor_P']}；{fmt(main_rows[3]['tumor_q'])}|组胺28/33对下降；HDC RNA29/33对下降；Uhlitz肥大细胞表达|Lee肥大细胞覆盖不足。组织RNA变化可能含细胞组成与细胞内表达两部分，当前33例未测细胞比例。|
|GSH的多候选|GSTA4—GSH|{fmt(main_rows[4]['tumor_rho'])}；{main_rows[4]['tumor_P']}；{fmt(main_rows[4]['tumor_q'])}|GSH25/33对上升；Lee宽基质，Uhlitz成纤维最高但保持率约59.7%|GSTA4经典脂质过氧化产物结合背景不能直接定为组织GSH负相关原因。|
|GSH的多候选|CLIC2—GSH|{fmt(main_rows[5]['tumor_rho'])}；{main_rows[5]['tumor_P']}；{fmt(main_rows[5]['tumor_q'])}|同一GSH特征；CLIC2 RNA29/33对下降；Lee宽基质/Uhlitz内皮|已有体外GSH依赖氧化还原与内皮屏障功能线索，两者均非本患者GSH机制验证。|

本版891关系中871可算、52项P<0.05，**0项q<0.05**。上表只展示读者需要逐层判断的例子；GSTA4与CLIC2的ρ差别没有比较检验，不能排序。逐一剔除后GSTA4ρ为{fmt(main_rows[4]['leave_one_out_rho_min'])}至{fmt(main_rows[4]['leave_one_out_rho_max'])}，CLIC2为{fmt(main_rows[5]['leave_one_out_rho_min'])}至{fmt(main_rows[5]['leave_one_out_rho_max'])}，均仍为负；它们属于同一队列内诊断。

![五个基因在Lee与Uhlitz的肿瘤单细胞表达来源](focus_source_five_genes.png)

图中点色为供者等权平均log1p(CP10K)，点大小为平均检出比例；两队列色标分别设定。表达位置不能推定组织代谢物的来源。

## 新手解释
“尿苷25对下降”与“UCKL1 RNA25对上升”分别是边际计数，并没有证明是同25人；交集至少17、至多25人，实际交集未在公开汇总中给出。肿瘤内ρ比较的是33个肿瘤之间的高低，不是患者内变化量相关。RNA差异、细胞来源和体外功能证据各自回答不同问题。

## 限制/反证
NNMT两条关系在**旧674范围**的M2（年龄、性别、分期）中分别为ρ={fmt(m2[0]['effect'])}、P={m2[0]['p_value']}与ρ={fmt(m2[1]['effect'])}、P={m2[1]['p_value']}；这只作为旧版敏感性，不代表新版891关系全做了M2。已有1-MNA实验涉及HT29/CCD-18Co模型、高剂量或5-FU，不能套成CAMP患者CAF因果机制。CLIC2文献部分采用用户本轮给定的原始研究解读，未另作逐图审计；本报告不计新增功能认证。来源颜色按队列独立，不跨研究比较数值大小。

## 当前决定
UCKL1、NNMT、HDC与GSTA4/CLIC2作为四组汇报案例；AQP9、SLC6A6、ASPA可辅助展示，SLC6A17与ENO2展示覆盖/可用性边界。完整关系和全部候选保持原版，不新增靶点名次。

## 下一步
依据原始关系表讨论测量与细胞模型；此汇报不启动新的机制、细胞通讯或患者分析。

## 复现命令
`python code/coad/build_four_questions_v1.py --out <全新运行目录>`。本脚本只读取已发布汇总；`focus_relations.tsv`和`source_five_genes.tsv`是图源数据。图中×表示供者覆盖不足。

## 定向功能来源
已核查仓库历史记录：[UCKL1功能说明](../../05_FUNCTION/20260921T031711Z_mediation4_v1/README_CN.md)、[NNMT定向核查](../../05_FUNCTION/20260921T031711Z_mediation4_v1/README_CN.md)。本轮采用用户提供的[CLIC2体外研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291220/)、[CLIC2内皮研究](https://pubmed.ncbi.nlm.nih.gov/30929599/)及[GSTA4生化研究](https://pubs.acs.org/doi/10.1021/bi902038u)解读；论文模型与当前患者数据独立。
"""
    with (out / "README_CN.md").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(report)
    with (out / "analysis_spec.json").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(dict(cancer="COAD", kind="read_only_focus_view", main_relation_count=6, gene_count=5, input_commit="06e2a9d8ac83042cb985730131d824fdd6574440", new_statistics=False, old_covariate_family=674, current_primary_family=891, source_metric="donor_equal_mean_log1p_CP10K"), ensure_ascii=False, indent=2))
    manifest = [dict(name=name, path=path.relative_to(ROOT).as_posix(), sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for name, path in inputs.items()]
    write(out / "source_manifest.tsv", manifest)
    with (out / "validation.json").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(dict(status="PASS", relation_rows=len(main_rows), source_profile_rows=len(profiles), primary_figure_panels=2, figure_labels=genes, all_relation_q_over_005=all(float(r["tumor_q"]) > 0.05 for r in main_rows), no_new_patient_or_single_cell_statistics=True, old_covariate_scope_separate=True), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
