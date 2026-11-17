"""Exploratory SASA controls for corrected-v2 between-PDB position IQR.

Receptor, not residue or PDB, is the inferential unit. No upstream writes.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/role_d/sasa_iqr"
ANALYSIS = ROOT / "data/rebuild_v2/analysis"
FUNCTIONAL = ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr.csv"
SEED = 20261029
N_BOOT = 10000
MIN_PDB = 5


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full(len(values), np.nan)
    finite = np.flatnonzero(np.isfinite(values))
    order = finite[np.argsort(values[finite])]
    if len(order):
        result[order] = np.minimum(1, np.minimum.accumulate((values[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1])
    return result


def summarize_effects(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"n_receptors": 0, "median_effect": np.nan, "ci95_low": np.nan,
                "ci95_high": np.nan, "wilcoxon_p": np.nan, "n_positive": 0}
    rng = np.random.default_rng(SEED)
    boot = np.median(values[rng.integers(len(values), size=(N_BOOT, len(values)))], axis=1)
    return {"n_receptors": len(values), "median_effect": np.median(values),
            "ci95_low": np.quantile(boot, .025), "ci95_high": np.quantile(boot, .975),
            "wilcoxon_p": wilcoxon(values).pvalue if np.any(values != 0) and len(values) >= 5 else np.nan,
            "n_positive": int((values > 0).sum())}


def context(data):
    element = data.canonical_elements.fillna("")
    return element.where(element.ne(""), "INTER:" + data.inter_helix_flank_before.fillna("NONE") + "-" + data.inter_helix_flank_after.fillna("NONE"))


def matched_effects(data, design):
    """Frozen observed strata; >=1 LBP and >=3 nonfunctional controls."""
    selected = data.loc[data.eligible & (data.is_lbp | ~data.is_any_functional)].copy()
    if design == "CONTEXT_SUPPORT_RSA_NO_REVIEW":
        selected = selected.loc[~selected.review_required_flag]
    selected["match_stratum"] = selected.context + "|" + selected.support_band
    if design != "CONTEXT_SUPPORT":
        selected["match_stratum"] += "|" + selected.rsa_band
    if design == "CONTEXT_SUPPORT_RSA_AMINOACID":
        selected["match_stratum"] += "|" + selected.canonical_residue_name
    rows = []
    for (uid, stratum), group in selected.groupby(["uniprot_id", "match_stratum"], sort=True):
        site, control = group.loc[group.is_lbp], group.loc[~group.is_lbp]
        if len(site) and len(control) >= 3:
            rows.append({"design": design, "uniprot_id": uid, "common_name": group.common_name.iloc[0],
                         "match_stratum": stratum, "n_site": len(site), "n_control": len(control),
                         "site_iqr": site.paired_iqr_z.median(), "control_iqr": control.paired_iqr_z.median(),
                         "effect": site.paired_iqr_z.median() - control.paired_iqr_z.median()})
    strata = pd.DataFrame(rows, columns=["design", "uniprot_id", "common_name", "match_stratum", "n_site", "n_control", "site_iqr", "control_iqr", "effect"])
    receptor = strata.groupby(["design", "uniprot_id", "common_name"], as_index=False).agg(
        effect=("effect", "median"), n_strata=("effect", "size"), n_site=("n_site", "sum"), n_control=("n_control", "sum"))
    return strata, receptor


def build_positions():
    sasa = pd.read_parquet(OUT / "pdb_position_sasa.parquet")
    bfactor = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")
    keys = ["uniprot_id", "pdb_id", "canonical_uniprot_position"]
    paired = bfactor.merge(sasa, on=keys, how="inner", validate="one_to_one")
    paired = paired.loc[np.isfinite(paired.b_factor_z) & np.isfinite(paired.rsa)]
    positions = paired.groupby(["uniprot_id", "canonical_uniprot_position"], as_index=False).agg(
        n_paired_pdb=("pdb_id", "nunique"), median_rsa=("rsa", "median"),
        median_sasa_a2=("sasa_a2", "median"), paired_median_z=("b_factor_z", "median"),
        paired_q25_z=("b_factor_z", lambda s: s.quantile(.25)), paired_q75_z=("b_factor_z", lambda s: s.quantile(.75)))
    positions["paired_iqr_z"] = positions.paired_q75_z - positions.paired_q25_z
    profile = pd.read_csv(ANALYSIS / "receptor_position_professor_profile_v2.csv")
    positions = profile.merge(positions, on=["uniprot_id", "canonical_uniprot_position"], how="left", validate="one_to_one")
    positions["n_paired_pdb"] = positions.n_paired_pdb.fillna(0).astype(int)
    positions["eligible"] = positions.n_paired_pdb.ge(MIN_PDB) & positions.paired_iqr_z.notna()
    positions["review_required_flag"] = positions.review_required_flag.astype(str).str.lower().eq("true")
    functional = pd.read_csv(FUNCTIONAL)
    all_keys = set(zip(functional.uniprot_id, functional.canonical_uniprot_position))
    pocket = functional.loc[functional.functional_site.eq("LIGAND_BINDING_POCKET")]
    lbp_keys = set(zip(pocket.uniprot_id, pocket.canonical_uniprot_position))
    h12 = functional.loc[functional.functional_site.eq("AF2_H12")]
    h12_keys = set(zip(h12.uniprot_id, h12.canonical_uniprot_position))
    keys2 = list(zip(positions.uniprot_id, positions.canonical_uniprot_position))
    positions["is_lbp"] = [k in lbp_keys for k in keys2]
    positions["is_af2_h12"] = [k in h12_keys for k in keys2]
    positions["is_any_functional"] = [k in all_keys for k in keys2]
    positions["context"] = context(positions)
    positions["support_band"] = pd.cut(positions.n_paired_pdb, [4, 9, 19, np.inf], labels=["5_9", "10_19", "20_PLUS"]).astype(str)
    positions["rsa_band"] = pd.cut(positions.median_rsa, [-np.inf, .05, .25, np.inf], labels=["BURIED_LE_005", "INTERMEDIATE_005_025", "EXPOSED_GT_025"]).astype(str)
    return positions, paired


def figures(positions, correlations, effects):
    plt.rcParams.update({"font.size": 11, "savefig.dpi": 180})
    fig, ax = plt.subplots(figsize=(9, 5))
    for i, (metric, group) in enumerate(correlations.groupby("endpoint", sort=True)):
        ax.scatter(np.linspace(i - .12, i + .12, len(group)), group.rho, alpha=.75)
        ax.plot([i - .2, i + .2], [group.rho.median()] * 2, color="black", lw=3)
    ax.set_xticks([0, 1], ["Position IQR(z) across PDB entries", "Median relative B-factor z"])
    ax.axhline(0, c="gray", ls="--")
    ax.set(ylabel="Within-receptor Spearman rho with median RSA", title="Solvent exposure: one point per receptor")
    fig.tight_layout(); fig.savefig(OUT / "01_sasa_correlations.png"); plt.close(fig)
    designs = effects.design.unique()
    fig, ax = plt.subplots(figsize=(10, 6))
    labels = []
    for i, design in enumerate(designs):
        group = effects.loc[effects.design.eq(design)]
        ax.scatter(group.effect, np.linspace(i - .12, i + .12, len(group)), alpha=.55, s=22)
        s = summarize_effects(group.effect)
        ax.plot([s['ci95_low'], s['ci95_high']], [i, i], color="black", lw=3)
        ax.scatter([s['median_effect']], [i], color="black", marker="D", zorder=4)
        labels.append(design.replace("CONTEXT_SUPPORT", "Context + support").replace("_", " ") + f" (n={len(group)})")
    ax.axvline(0, c="gray", ls="--")
    ax.set_yticks(range(len(designs)), labels)
    ax.set(xlabel="LBP − nonfunctional matched-control position IQR(z)", title="Exploratory contrasts; receptor medians and bootstrap 95% CIs")
    fig.tight_layout(); fig.savefig(OUT / "02_sasa_adjusted_lbp_iqr.png"); plt.close(fig)
    top = positions.loc[positions.eligible].groupby("common_name").size().sort_values(ascending=False).head(6).index
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharex=True, sharey=True)
    for name, ax in zip(top, axes.flat):
        d = positions.loc[positions.eligible & positions.common_name.eq(name)]
        for flag, label, color in [(~d.is_any_functional, "Other", "#b2bec3"), (d.is_lbp, "LBP", "#0984e3"), (d.is_af2_h12, "AF2/H12", "#d63031")]:
            ax.scatter(d.loc[flag, "median_rsa"], d.loc[flag, "paired_iqr_z"], s=13, alpha=.7, label=label, color=color)
        ax.set_title(f"{name} (positions n={len(d)})")
        ax.set(xlabel="Median isolated-chain RSA", ylabel="Between-PDB IQR(z)")
    axes.flat[0].legend(fontsize=8)
    fig.suptitle("Six receptors with most eligible positions (descriptive; no pooled residue test)")
    fig.tight_layout(); fig.savefig(OUT / "03_sasa_iqr_examples.png"); plt.close(fig)


def site_heatmap(sites):
    names = sorted(sites.common_name.unique())
    columns = ["AF2_H12", "LIGAND_BINDING_POCKET"]
    fig, axes = plt.subplots(1, 2, figsize=(9, max(8, len(names) * .24)), sharey=True)
    for ax, metric, title, cmap in zip(axes,
            ["median_position_iqr", "median_position_rsa"],
            ["Median position IQR(z)", "Median position RSA"], ["viridis", "YlOrRd"]):
        table = sites.pivot(index="common_name", columns="functional_site", values=metric).reindex(index=names, columns=columns)
        color = plt.get_cmap(cmap).copy(); color.set_bad("#dddddd")
        image = ax.imshow(table.to_numpy(), aspect="auto", cmap=color, vmin=0)
        ax.set_xticks([0, 1], ["AF2/H12", "LBP"])
        ax.set_yticks(range(len(names)), names, fontsize=9)
        ax.set_title(title)
        for i, row in enumerate(table.to_numpy()):
            for j, value in enumerate(row):
                if np.isfinite(value):
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8,
                            color="white" if value > np.nanmax(table.to_numpy()) * .6 else "black")
        fig.colorbar(image, ax=ax, shrink=.5)
    fig.suptitle("Functional sites: descriptive receptor summaries; gray = unavailable")
    fig.tight_layout(); fig.savefig(OUT / "04_functional_site_iqr_rsa_heatmap.png"); plt.close(fig)


def main():
    positions, paired = build_positions()
    positions.to_csv(OUT / "position_sasa_iqr.csv", index=False)
    correlations = []
    for uid, group in positions.loc[positions.eligible].groupby("uniprot_id"):
        if len(group) < 20:
            continue
        for endpoint, column in [("IQR_Z", "paired_iqr_z"), ("MEDIAN_Z", "paired_median_z")]:
            if group[column].nunique() > 1 and group.median_rsa.nunique() > 1:
                correlations.append({"uniprot_id": uid, "common_name": group.common_name.iloc[0],
                                     "endpoint": endpoint, "n_positions": len(group),
                                     "rho": spearmanr(group.median_rsa, group[column]).statistic})
    correlations = pd.DataFrame(correlations)
    correlations.to_csv(OUT / "receptor_sasa_correlations.csv", index=False)
    cs = pd.DataFrame([dict(endpoint=endpoint, **summarize_effects(g.rho)) for endpoint, g in correlations.groupby("endpoint")])
    cs["bh_q_two_correlations"] = bh(cs.wilcoxon_p)
    cs.to_csv(OUT / "sasa_correlation_summary.csv", index=False)
    strata, effects, summaries = [], [], []
    for design in ["CONTEXT_SUPPORT", "CONTEXT_SUPPORT_RSA", "CONTEXT_SUPPORT_RSA_NO_REVIEW", "CONTEXT_SUPPORT_RSA_AMINOACID"]:
        st, ef = matched_effects(positions, design)
        strata.append(st); effects.append(ef)
        summaries.append(dict(design=design, n_strata=len(st), n_site_positions=int(st.n_site.sum()),
                              n_control_positions=int(st.n_control.sum()), **summarize_effects(ef.effect)))
    strata = pd.concat(strata, ignore_index=True)
    effects = pd.concat(effects, ignore_index=True)
    summaries = pd.DataFrame(summaries)
    summaries["bh_q_four_designs"] = bh(summaries.wilcoxon_p)
    strata.to_csv(OUT / "lbp_sasa_matched_strata.csv", index=False)
    effects.to_csv(OUT / "lbp_sasa_receptor_effects.csv", index=False)
    summaries.to_csv(OUT / "lbp_sasa_summary.csv", index=False)
    site_rows = []
    for flag, site in [("is_lbp", "LIGAND_BINDING_POCKET"), ("is_af2_h12", "AF2_H12")]:
        selected = positions.loc[positions[flag] & positions.eligible]
        g = selected.groupby(["uniprot_id", "common_name"], as_index=False).agg(
            n_positions=("paired_iqr_z", "size"), median_position_iqr=("paired_iqr_z", "median"),
            median_position_rsa=("median_rsa", "median"), median_position_z=("paired_median_z", "median"))
        g["functional_site"] = site
        site_rows.append(g)
    sites = pd.concat(site_rows, ignore_index=True)
    sites.to_csv(OUT / "functional_site_sasa_iqr_summary.csv", index=False)
    figures(positions, correlations, effects)
    site_heatmap(sites)
    inputs = [ANALYSIS / "pdb_position_bfactor_profiles.parquet", ANALYSIS / "receptor_position_professor_profile_v2.csv", FUNCTIONAL, OUT / "pdb_position_sasa.parquet"]
    audit = {"random_seed": SEED, "bootstrap_replicates": N_BOOT, "min_paired_pdb_per_position": MIN_PDB,
             "n_paired_entry_positions": len(paired), "n_eligible_positions": int(positions.eligible.sum()),
             "n_receptors_with_eligible_positions": positions.loc[positions.eligible, "uniprot_id"].nunique(),
             "input_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
             "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / "analysis_manifest.json").write_text(json.dumps(audit, indent=2) + "\n")
    cache_audit = pd.read_csv(OUT / "coordinate_cache_audit.csv")
    excluded_pdbs = cache_audit.loc[~cache_audit.status.eq("VERIFIED"), "pdb_id"].tolist()
    calc_audit = pd.read_csv(OUT / "calculation_audit.csv")
    failed_pdbs = calc_audit.loc[~calc_audit.status.eq("OK"), "pdb_id"].tolist()
    lines = ["# SASA + IQR: результаты для команды (29 октября 2026)", "",
             "**Новый exploratory-анализ исправленного v2; старые результаты не заменяет.**", "",
             f"SASA и B-factor сопоставлены для {len(paired):,} PDB × residue наблюдений. "
             f"Основной порог: ≥5 PDB с полными тяжёлыми атомами; доступны {audit['n_eligible_positions']:,} позиций "
             f"в {audit['n_receptors_with_eligible_positions']} рецепторах.", "",
             f"Координаты: {int(cache_audit.status.eq('VERIFIED').sum())}/{len(cache_audit)} PDB совпали с frozen SHA256. "
             f"Без verified coordinates: {', '.join(excluded_pdbs) or 'нет'}. "
             f"Calculation failures: {', '.join(failed_pdbs) or 'нет'}. "
             "Эти записи сохранены в master, но не участвуют в новом SASA-анализе.", "",
             "## 1. Связь экспозиции с B-factor и IQR", "",
             "Корреляция Spearman рассчитана отдельно в каждом рецепторе (≥20 позиций). "
             "Доверительный интервал — bootstrap рецепторов, не остатков. "
             "Wilcoxon проверяет распределение receptor-level rho относительно нуля; BH для двух endpoints.", ""]
    for r in cs.itertuples():
        lines.append(f"- {r.endpoint}: n={r.n_receptors}, медианная rho={r.median_effect:.3f}, "
                     f"95% CI [{r.ci95_low:.3f}; {r.ci95_high:.3f}], q={r.bh_q_two_correlations:.4g}.")
    lines += ["", "## 2. Лигандный карман: IQR после matching", "",
              "Сравниваются LBP и позиции вне LBP/AF2-H12 одного рецептора: точный professor-контекст, "
              "support 5–9/10–19/≥20, затем RSA ≤0.05 / (0.05,0.25] / >0.25. "
              "≥1 site и ≥3 controls на страту; медиана разностей по стратам → один эффект на рецептор. "
              "Отрицательный эффект = меньший IQR у LBP. P — receptor-level Wilcoxon; BH для четырёх designs.", ""]
    for r in summaries.itertuples():
        lines.append(f"- **{r.design}**: n={r.n_receptors}; site positions={r.n_site_positions}; "
                     f"ΔIQR={r.median_effect:.4f}, 95% CI [{r.ci95_low:.4f}; {r.ci95_high:.4f}]; "
                     f"p={r.wilcoxon_p:.4g}; q={r.bh_q_four_designs:.4g}.")
    primary = summaries.loc[summaries.design.eq("CONTEXT_SUPPORT_RSA")].iloc[0]
    supported = (primary.ci95_low > 0 or primary.ci95_high < 0) and primary.bh_q_four_designs < .05
    lines += ["", "## Вывод для презентации", "",
              ("После coarse RSA/context matching наблюдается ненулевой LBP-контраст в этом exploratory-анализе. "
               "Нужно учитывать sensitivity designs и ограничения, прежде чем делать биологический вывод."
               if supported else "После учёта экспозиции и структурного контекста убедительного LBP-специфичного различия IQR не установлено "
               "по совместному критерию CI вне нуля и BH q<0.05. Это не доказывает отсутствия эффекта."), "",
              "**IQR(z) ≠ median(z).** IQR показывает неоднородность относительного кристаллографического disorder "
              "между PDB-записями, а не абсолютную гибкость или динамику в растворе. "
              "Корреляции с median(z) и с IQR(z) отвечают на разные вопросы.", "",
              "## Что нельзя утверждать", "",
              "- Нельзя объявлять ligand-induced эффект: это cross-sectional association, не APO/HOLO intervention.",
              "- AF2/H12 описан по IQR и RSA, но exact-context site-specific тест не выполнялся: весь H12 определяет сайт, нет независимого same-H12 background.",
              "- Это isolated-chain SASA, не SASA комплекса/биологической сборки/кристаллической решётки; ligand/partner occlusion не измерена.",
              "- Неполные side chains и неканонические residue types не получают RSA. Missingness и coarse bins могут влиять на результат.",
              "- Разрешение, refinement, packing, ligand composition и пространственная зависимость позиций полностью не контролируются; confidence intervals не снимают эти ограничения.",
              "- Matching меняет набор доступных рецепторов/позиций; сравнение designs не является формальной оценкой mediation/attenuation.", "",
              "## Графики", "", "1. `01_sasa_correlations.png` — receptor-level корреляции.",
              "2. `02_sasa_adjusted_lbp_iqr.png` — эффекты matching и CIs.",
              "3. `03_sasa_iqr_examples.png` — шесть рецепторов с максимальным числом eligible positions, без pooled-тестов.",
              "4. `04_functional_site_iqr_rsa_heatmap.png` — IQR/RSA функциональных сайтов по рецепторам, описательно.", "",
              "Методы, команды и дополнительные ограничения: `README.md`. Provenance: `sasa_manifest.json`, `analysis_manifest.json` и audit CSV."]
    (OUT / "RESULTS_RU.md").write_text("\n".join(lines) + "\n")
    slide = ["# Slide text — IQR and solvent exposure (29 October 2026)", "",
             "## Slide 1 — What was calculated", "",
             "- Corrected v2 numbering and professor structural-element annotations.",
             f"- {audit['n_eligible_positions']:,} positions in {audit['n_receptors_with_eligible_positions']} receptors; ≥5 paired PDB entries per position.",
             "- Between-PDB IQR of within-instance-normalized Cα B-factor z, not solution dynamics.",
             "- Isolated-chain SASA (Shrake–Rupley; 1.4 Å probe); complete residues only; Tien-normalized RSA.",
             "- Figure: 04_functional_site_iqr_rsa_heatmap.png (descriptive).", "",
             "## Slide 2 — Solvent exposure is a potential confounder", ""]
    for r in cs.itertuples():
        slide.append(f"- {r.endpoint}: median within-receptor Spearman rho {r.median_effect:.3f} "
                     f"(receptor-bootstrap 95% CI {r.ci95_low:.3f} to {r.ci95_high:.3f}; n={r.n_receptors}; BH q={r.bh_q_two_correlations:.3g}).")
    slide += ["- Figure: 01_sasa_correlations.png; each point is one receptor.", "",
              "## Slide 3 — Ligand-pocket contrast after exposure matching", "",
              f"- Context + support + RSA matched LBP − controls: median ΔIQR {primary.median_effect:.4f} "
              f"(95% CI {primary.ci95_low:.4f} to {primary.ci95_high:.4f}; n={primary.n_receptors} receptors; BH q={primary.bh_q_four_designs:.3g}).",
              ("- Exploratory nonzero contrast; inspect sensitivity designs before generalizing." if supported else
               "- No convincing LBP-specific contrast under the combined CI and BH criterion; not evidence of equivalence."),
              "- Figure: 02_sasa_adjusted_lbp_iqr.png; cohort/support changes across matching designs.", "",
              "## Slide 4 — Interpretation and limitations", "",
              "- IQR is heterogeneity across deposited structures; median z is relative disorder within structures.",
              "- Exposure matching is coarse and exploratory; packing, resolution, refinement and composition remain potential confounders.",
              "- AF2/H12: descriptive only here; no independent non-site same-H12 background.",
              "- No causal ligand-induced mechanism or solution-dynamics conclusion is established."]
    (OUT / "SLIDE_TEXT_EN.md").write_text("\n".join(slide) + "\n")
    # Small shareable package: no coordinate cache or large per-residue parquet.
    archive = OUT.parent / "Role_D_IQR_SASA_2026-10-29.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as package:
        for path in sorted(OUT.iterdir()):
            if path.suffix in {".md", ".png", ".csv", ".json", ".txt"}:
                package.write(path, "sasa_iqr/" + path.name)
        for path in [ROOT / "scripts/113_role_d_sasa.py", Path(__file__), ROOT / "tests/test_role_d_sasa_iqr.py"]:
            package.write(path, str(path.relative_to(ROOT)))
    print(cs.to_string(index=False)); print(summaries.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
