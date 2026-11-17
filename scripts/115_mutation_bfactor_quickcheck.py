"""Quick check: are engineered mutations associated with different normalized B-factors?

Observational and descriptive. Reads the corrected v2 catalogue tables and v2 profiles;
writes only reports/mutation_bfactor_quickcheck/. No analysis output is changed.

Units
- Entry = receptor x PDB entry with z-scores. Classes from catalogue_mutation_audit_v2.csv:
  ENGINEERED_MUTANT vs WT_NO_ENGINEERED_SUBSTITUTION; UNCERTAIN excluded (never WT).
- Per-position z = v2 entry-level b_factor_z (per-chain z, median over the entry's
  representative chains), canonical v2 UniProt numbering. Missing stays missing.
- Engineered sites = catalogue_sequence_differences_v2.csv rows with difference_class
  ENGINEERED (natural variants and modified residues are not used).

Q1 entry-level shape: median, IQR and 90th percentile of z per entry (the chain mean is
   0 by construction and is not compared). TRACE_LBD-only fragments are excluded.
Q2/Q3 local, position-matched: for each engineered site p in a mutant entry and offset d,
   diff = z_mutant(p+d) - median z of WT entries of the same receptor at p+d
   (>= 5 WT entries observing that position). Summarized per mutation observation by
   distance bin, then across observations.
Q4 ecological: Spearman of engineered-mutant fraction vs receptor median IQR(z).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "reports" / "catalogue_v2" / "tables"
OUT = ROOT / "reports" / "mutation_bfactor_quickcheck"
MIN_WT = 5            # WT entries needed at a position for a matched reference (project rule)
MIN_GROUP = 5         # mutant and WT entries needed for a within-receptor Q1 test
WINDOW = 10
FAR = (11, 20)        # distal comparison bin
BINS = [("0", 0, 0), ("1-2", 1, 2), ("3-5", 3, 5), ("6-10", 6, 10), ("11-20", *FAR)]
N_BOOT = 2000
RNG = np.random.default_rng(363)
MUT, WT = "ENGINEERED_MUTANT", "WT_NO_ENGINEERED_SUBSTITUTION"

INK, MUTED, GRID = "#1d2533", "#5f6876", "#e9ecf0"
DARK, MID, LIGHT = "#3c4a5e", "#8ea4c2", "#cdd5e0"
plt.rcParams.update({"font.family": "Helvetica Neue", "font.size": 13, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#9aa3af", "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.labelcolor": INK, "figure.facecolor": "white"})


def boot_ci(x, f=np.median):
    x = np.asarray(x, float)
    if len(x) < 3:
        return np.nan, np.nan
    b = f(RNG.choice(x, (N_BOOT, len(x)), replace=True), axis=1)
    return tuple(np.percentile(b, [2.5, 97.5]))


def signed_rank(x):
    """Wilcoxon signed-rank vs 0 and matched-pairs rank-biserial r."""
    x = np.asarray(x, float)
    x = x[x != 0]
    if len(x) < 5:
        return np.nan, np.nan
    r = stats.rankdata(np.abs(x))
    rb = (r[x > 0].sum() - r[x < 0].sum()) / r.sum()
    return stats.wilcoxon(x).pvalue, rb


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return u.pvalue, 2 * u.statistic / (len(a) * len(b)) - 1


def van_elteren(df, val):
    """Receptor-stratified rank-sum test (mutant vs WT), weights 1/(n+1)."""
    T = E = V = 0.0
    for _, g in df.groupby("uniprot_id"):
        m, n = (g["cls"] == MUT).sum(), (g["cls"] == WT).sum()
        if m == 0 or n == 0:
            continue
        N = m + n
        r = stats.rankdata(g[val])
        w = 1 / (N + 1)
        T += w * r[(g["cls"] == MUT).to_numpy()].sum()
        E += w * m * (N + 1) / 2
        ties = np.unique(r, return_counts=True)[1]
        V += w * w * m * n / 12 * ((N + 1) - (ties ** 3 - ties).sum() / (N * (N - 1)))
    z = (T - E) / np.sqrt(V)
    return 2 * stats.norm.sf(abs(z)), z


def load():
    audit = pd.read_csv(CAT / "catalogue_mutation_audit_v2.csv", dtype={"pdb_id": str}, keep_default_na=False)
    sd = pd.read_csv(CAT / "catalogue_sequence_differences_v2.csv", dtype={"pdb_id": str}, keep_default_na=False)
    S = pd.read_csv(CAT / "receptor_catalogue_summary_v2.csv")
    ctx = pd.read_csv(CAT / "catalogue_entry_context_v2.csv", dtype={"pdb_id": str}, keep_default_na=False,
                      usecols=["uniprot_id", "pdb_id", "catalogue_state", "resolution"])
    prof = pd.read_parquet(ROOT / "data" / "rebuild_v2" / "analysis" / "pdb_position_bfactor_profiles.parquet",
                           columns=["uniprot_id", "common_name", "pdb_id", "canonical_uniprot_position", "b_factor_z"])
    prof = prof.dropna(subset=["b_factor_z"])
    assert audit[["uniprot_id", "pdb_id"]].drop_duplicates().shape[0] == len(audit) == 1805
    return audit, sd, S, ctx, prof


def q1(audit, ctx, prof):
    m = prof.groupby(["uniprot_id", "pdb_id"])["b_factor_z"].agg(
        median_z="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75),
        p90_z=lambda x: x.quantile(0.9), n_pos="size").reset_index()
    m["iqr_z"] = m["q75"] - m["q25"]
    d = m.merge(audit[["uniprot_id", "receptor", "pdb_id", "final_classification"]], on=["uniprot_id", "pdb_id"])
    d = d.merge(ctx, on=["uniprot_id", "pdb_id"])
    assert len(d) == 1805
    n_trace = int((d["catalogue_state"] == "TRACE_LBD").sum())
    d = d[(d["catalogue_state"] != "TRACE_LBD") & d["final_classification"].isin([MUT, WT])].rename(
        columns={"final_classification": "cls"})
    rows = []
    for val in ["median_z", "iqr_z", "p90_z"]:
        a, b = d.loc[d["cls"] == MUT, val], d.loc[d["cls"] == WT, val]
        p, delta = cliffs(a, b)
        pv, zv = van_elteren(d, val)
        rows.append(dict(question="Q1", analysis="all receptors pooled", receptor="all", metric=val, n_mutant=len(a),
                         n_wt=len(b), mutant_median=a.median(), wt_median=b.median(), difference=a.median() - b.median(),
                         effect_size=delta, effect_size_type="Cliff's delta", p_value=p, test="Mann-Whitney U",
                         stratified_p=pv, stratified_test="van Elteren (receptor strata)"))
        for (uid, rec), g in d.groupby(["uniprot_id", "receptor"]):
            a, b = g.loc[g["cls"] == MUT, val], g.loc[g["cls"] == WT, val]
            if len(a) >= MIN_GROUP and len(b) >= MIN_GROUP:
                p, delta = cliffs(a, b)
                rows.append(dict(question="Q1", analysis="within receptor", receptor=rec, metric=val, n_mutant=len(a),
                                 n_wt=len(b), mutant_median=a.median(), wt_median=b.median(),
                                 difference=a.median() - b.median(), effect_size=delta,
                                 effect_size_type="Cliff's delta", p_value=p, test="Mann-Whitney U"))
    R = pd.DataFrame(rows)
    w = R["analysis"] == "within receptor"
    for val in R["metric"].unique():
        k = w & (R["metric"] == val)
        R.loc[k, "p_holm"] = holm(R.loc[k, "p_value"].to_numpy())
    res = {"n_trace_excluded": n_trace, "res_mut": d.loc[d["cls"] == MUT, "resolution"].astype(float).median(),
           "res_wt": d.loc[d["cls"] == WT, "resolution"].astype(float).median()}
    return R, d, res


def holm(p):
    o = np.argsort(p)
    adj = np.empty_like(p, dtype=float)
    run = 0
    for i, j in enumerate(o):
        run = max(run, min(1, (len(p) - i) * p[j]))
        adj[j] = run
    return adj


def local(audit, sd, prof):
    cls = audit.set_index(["uniprot_id", "pdb_id"])["final_classification"]
    eng = sd[sd["difference_class"] == "ENGINEERED"][["uniprot_id", "pdb_id", "position", "label"]].drop_duplicates(
        ["uniprot_id", "pdb_id", "position"])
    eng = eng[[cls.get((u, p)) == MUT for u, p in zip(eng["uniprot_id"], eng["pdb_id"])]]
    # WT reference: median z per receptor x position over WT entries, with >= MIN_WT entries
    key = list(zip(prof["uniprot_id"], prof["pdb_id"]))
    prof = prof.assign(cls=[cls.get(k) for k in key])
    wt = prof[prof["cls"] == WT].groupby(["uniprot_id", "canonical_uniprot_position"])["b_factor_z"].agg(["median", "size"])
    wt = wt[wt["size"] >= MIN_WT]["median"]
    z = prof.set_index(["uniprot_id", "pdb_id", "canonical_uniprot_position"])["b_factor_z"]
    sites = eng.groupby(["uniprot_id", "pdb_id"])["position"].apply(list).to_dict()
    names = prof.drop_duplicates("uniprot_id").set_index("uniprot_id")["common_name"]
    obs = []
    for _, r in eng.iterrows():
        u, pdb, p = r["uniprot_id"], r["pdb_id"], int(r["position"])
        others = [q for q in sites[(u, pdb)] if q != p]
        for dd in range(-FAR[1], FAR[1] + 1):
            q = p + dd
            zm, zw = z.get((u, pdb, q), np.nan), wt.get((u, q), np.nan)
            obs.append(dict(uniprot_id=u, receptor=names.get(u, u), pdb_id=pdb, site=p, label=r["label"], offset=dd,
                            position=q, z_mutant=zm, wt_median=zw, diff=zm - zw,
                            near_other_site=any(abs(q - o) <= WINDOW for o in others),
                            site_isolated=not any(abs(p - o) <= WINDOW for o in others)))
    return pd.DataFrame(obs), eng


def summarize_local(L):
    L = L.assign(dist=L["offset"].abs())
    rows = []
    variants = [("all mutation observations", L),
                ("isolated sites (no other engineered site within ±10)", L[L["site_isolated"]]),
                ("excluding ERα (ERα holds most observations)", L[L["receptor"] != "ERα"]),
                ("one value per receptor × mutation label", None)]
    for name, D in variants:
        for b, lo, hi in BINS + [("window ±2", 0, 2), ("window ±5", 0, 5), ("window ±10", 0, 10)]:
            src = L if D is None else D
            k = src[(src["dist"] >= lo) & (src["dist"] <= hi)].dropna(subset=["diff"])
            per = k.groupby(["uniprot_id", "pdb_id", "site", "label"])["diff"].mean().reset_index()
            if D is None:
                per = per.groupby(["uniprot_id", "label"])["diff"].median().reset_index()
            x = per["diff"].to_numpy()
            lo_ci, hi_ci = boot_ci(x)
            p, rb = signed_rank(x)
            rows.append(dict(question="Q2/Q3", analysis=name, distance=b, n_mutation_observations=len(x),
                             n_receptors=per["uniprot_id"].nunique(), median_difference=np.median(x) if len(x) else np.nan,
                             ci95_low=lo_ci, ci95_high=hi_ci, fraction_positive=(x > 0).mean() if len(x) else np.nan,
                             p_value=p, test="Wilcoxon signed-rank vs 0", effect_size=rb,
                             effect_size_type="matched-pairs rank-biserial"))
    return pd.DataFrame(rows)


def per_receptor_local(L):
    k = L[(L["offset"].abs() <= 2)].dropna(subset=["diff"])
    per = k.groupby(["receptor", "pdb_id", "site", "label"])["diff"].mean().reset_index()
    rows = []
    for rec, g in per.groupby("receptor"):
        if len(g) < 10:
            continue
        p, rb = signed_rank(g["diff"])
        lo, hi = boot_ci(g["diff"])
        top = g.groupby("label")["diff"].agg(["median", "size"]).sort_values("size", ascending=False).head(3)
        rows.append(dict(question="Q2 by receptor", analysis="window ±2, per receptor", receptor=rec,
                         n_mutation_observations=len(g), n_unique_mutations=g["label"].nunique(),
                         median_difference=g["diff"].median(), ci95_low=lo, ci95_high=hi, p_value=p,
                         effect_size=rb, effect_size_type="matched-pairs rank-biserial",
                         most_frequent_mutations="; ".join(f"{i} (n={int(r['size'])}, {r['median']:+.2f})"
                                                           for i, r in top.iterrows())))
    R = pd.DataFrame(rows)
    R["p_holm"] = holm(R["p_value"].to_numpy())
    return R


def q4(S):
    A = S[S["heterogeneity_class"] != "not assessable"]
    rho, p = stats.spearmanr(A["engineered_mutant_fraction"], A["median_iqr_z"])
    excl = A["n_pdb_engineered_mutant"] / (A["n_pdb_engineered_mutant"] + A["n_pdb_no_engineered_substitution"])
    rho2, p2 = stats.spearmanr(excl, A["median_iqr_z"])
    rows = [dict(question="Q4", analysis="receptors with assessable IQR", receptor=f"{len(A)} receptors",
                 metric="engineered-mutant fraction (of z-scored entries) vs median IQR(z)", effect_size=rho,
                 effect_size_type="Spearman rho", p_value=p, n_mutant=len(A)),
            dict(question="Q4", analysis="receptors with assessable IQR, UNCERTAIN removed from denominator",
                 receptor=f"{len(A)} receptors", metric="mutant / (mutant + WT) vs median IQR(z)", effect_size=rho2,
                 effect_size_type="Spearman rho", p_value=p2, n_mutant=len(A))]
    return pd.DataFrame(rows), A


def figures(L, Q1, A, d):
    # 1 offset profile
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.axhline(0, color=MUTED, lw=1, ls=(0, (4, 3)))
    for sub, c, lab in [(L, DARK, "all receptors"), (L[L["receptor"] != "ERα"], MID, "excluding ERα")]:
        k = sub.dropna(subset=["diff"])
        g = k.groupby("offset")["diff"]
        med = g.median()
        ci = np.array([boot_ci(v.to_numpy()) for _, v in g])
        ax.fill_between(med.index, ci[:, 0], ci[:, 1], color=c, alpha=0.18, lw=0)
        n0 = int(k.loc[k["offset"] == 0, "diff"].size)
        ax.plot(med.index, med.values, color=c, lw=2, marker="o", ms=4, label=f"{lab} ({n0} sites at 0)")
    ax.legend(frameon=False, loc="lower right")
    ax.axvspan(-0.5, 0.5, color=MID, alpha=0.25, lw=0)
    ax.set_xlabel("distance from engineered residue (residues, signed)")
    ax.set_ylabel("mutant z − WT median z\n(same receptor, same position)")
    ax.set_title("Position-matched difference around engineered sites (median, bootstrap 95% CI)", fontsize=13,
                 loc="left", color=INK)
    ax.grid(axis="y", color=GRID)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_local_effect_by_offset.png", dpi=150)
    plt.close(fig)
    # 2 within-receptor Q1 effect sizes
    W = Q1[Q1["analysis"] == "within receptor"]
    recs = W.groupby("receptor")["n_mutant"].max().sort_values().index.tolist()
    fig, axes = plt.subplots(1, 3, figsize=(12, 0.45 * len(recs) + 1.8), sharey=True)
    for ax, (val, lab) in zip(axes, [("median_z", "median z"), ("iqr_z", "IQR of z"), ("p90_z", "90th percentile z")]):
        w = W[W["metric"] == val].set_index("receptor").reindex(recs)
        y = np.arange(len(recs))
        ax.axvline(0, color=MUTED, lw=1, ls=(0, (4, 3)))
        ax.scatter(w["effect_size"], y, s=50, color=[DARK if p < 0.05 else MID for p in w["p_holm"]], zorder=3)
        pool = Q1[(Q1["analysis"] == "all receptors pooled") & (Q1["metric"] == val)]["effect_size"].iloc[0]
        ax.axvline(pool, color=DARK, lw=1, alpha=0.4)
        ax.set_xlim(-1, 1)
        ax.set_title(lab, fontsize=13, color=INK)
        ax.set_yticks(y, [f"{r}  ({int(w.loc[r, 'n_mutant'])} / {int(w.loc[r, 'n_wt'])})" for r in recs])
        ax.grid(axis="x", color=GRID)
        ax.set_xlabel("Cliff's delta (mutant vs WT)")
    fig.suptitle("Entry-level z distribution: mutant vs WT within receptor (dark = Holm p < 0.05)",
                 fontsize=13, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_entry_level_within_receptor.png", dpi=150)
    plt.close(fig)
    # 3 ecological scatter
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.scatter(A["engineered_mutant_fraction"] * 100, A["median_iqr_z"], s=40, color=DARK)
    for _, r in A.iterrows():
        if not (r["engineered_mutant_fraction"] >= 0.15 or r["median_iqr_z"] >= 0.7 or r["median_iqr_z"] <= 0.35):
            continue
        ax.annotate(r["common_name"], (r["engineered_mutant_fraction"] * 100, r["median_iqr_z"]), fontsize=9,
                    color=MUTED, xytext=(4, 2), textcoords="offset points")
    ax.set_xlabel("engineered-mutant entries (% of z-scored entries)")
    ax.set_ylabel("receptor median position IQR(z)")
    rho, p = stats.spearmanr(A["engineered_mutant_fraction"], A["median_iqr_z"])
    ax.set_title(f"Receptor level, {len(A)} receptors (descriptive only): Spearman rho = {rho:.2f}, p = {p:.2f}", fontsize=13, loc="left", color=INK)
    ax.grid(color=GRID)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_mutant_fraction_vs_iqr.png", dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    audit, sd, S, ctx, prof = load()
    Q1, d, res = q1(audit, ctx, prof)
    L, eng = local(audit, sd, prof)
    LS = summarize_local(L)
    LR = per_receptor_local(L)
    Q4, A = q4(S)
    figures(L, Q1, A, d)
    R = pd.concat([Q1, LS, LR, Q4], ignore_index=True)
    R.to_csv(OUT / "mutation_bfactor_quickcheck_results.csv", index=False, float_format="%.4g")
    site = L[L["offset"] == 0]
    meta = dict(n_engineered_sites=len(eng), n_sites_observed=int(site["z_mutant"].notna().sum()),
                n_sites_with_wt_ref=int(site["wt_median"].notna().sum()), n_sites_compared=int(site["diff"].notna().sum()),
                n_mutant_entries_with_site=eng[["uniprot_id", "pdb_id"]].drop_duplicates().shape[0],
                n_isolated_sites=int(site["site_isolated"].sum()), **res)
    pd.Series(meta).to_csv(OUT / "run_metadata.csv", header=["value"])
    print(pd.Series(meta).to_string())


if __name__ == "__main__":
    main()
