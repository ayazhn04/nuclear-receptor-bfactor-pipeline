"""Digest of the two "Dataset at a glance" catalogue slides: one question per page.

Reads only reports/catalogue_v2/tables/receptor_catalogue_summary_v2.csv (the table the
glance slides are drawn from) and recomputes nothing: every value on a page is a column
of that table, or a sum / share of its columns. Writes to reports/catalogue_digest/:

  figures/digest_<n>_<name>.png          full 1920x1080 pages (= PDF pages)
  figures/chart_only/digest_<n>_<name>.png the same page below the title band (deck images)
  catalogue_digest.pdf
  digest_values.csv                      every value shown, per receptor and page
  digest_page_text.csv                   each page's title and one-line message
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import matplotlib.lines

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from PIL import Image  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC = PROJECT_ROOT / "reports" / "catalogue_v2" / "tables" / "receptor_catalogue_summary_v2.csv"
OUT = PROJECT_ROOT / "reports" / "catalogue_digest"
FIG = OUT / "figures"
HEADER_Y = 0.815          # everything above this is the title band (cut off for the deck images)

INK, MUTED, FAINT, GRID, BG = "#1d2533", "#5f6876", "#9aa3af", "#e9ecf0", "#ffffff"
DARK, MID, LIGHT, PALE = "#3c4a5e", "#8ea4c2", "#cdd5e0", "#eef1f5"
plt.rcParams.update({"font.family": "Helvetica Neue", "font.size": 16, "axes.edgecolor": FAINT,
                     "xtick.color": MUTED, "ytick.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": BG, "axes.facecolor": BG,
                     "savefig.facecolor": BG})


def load() -> pd.DataFrame:
    S = pd.read_csv(SRC)
    assert len(S) == 48 and S["common_name"].is_unique
    return S.set_index("common_name", drop=False)


TEXT: list[dict] = []   # page titles and messages, reused verbatim by the editable deck


def page(title: str, message: str) -> plt.Figure:
    TEXT.append({"title": title, "message": message})
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.text(0.04, 0.925, title, fontsize=36, fontweight="bold", color=INK, va="baseline")
    fig.text(0.04, 0.868, message, fontsize=22, color=MUTED, va="baseline")
    return fig


def hbar_axes(fig, rect, names, bold=()):
    ax = fig.add_axes(rect)
    y = np.arange(len(names))[::-1]
    ax.set_ylim(-0.7, len(names) - 0.3)
    ax.set_yticks(y, names)
    for t in ax.get_yticklabels():
        t.set_fontsize(17 if len(names) <= 20 else 15)
        if t.get_text() in bold:
            t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.tick_params(axis="x", labelsize=14, length=0)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    return ax, y


def names_block(fig, x, y, heading, lines, size=17, width=None):
    fig.text(x, y, heading, fontsize=size, fontweight="bold", color=INK, va="top")
    fig.text(x, y - 0.042, "\n".join(lines), fontsize=size - 1, color=MUTED, va="top", linespacing=1.55)


# ---------------------------------------------------------------------------
# pages
# ---------------------------------------------------------------------------


def p1_entries(S, rec):
    z, t = S["n_entries_with_zscores"], S["n_entries_total"]
    top2 = z.sort_values(ascending=False).index[:2]
    share = z[top2].sum() / z.sum()
    fig = page("How much usable structural data is there per receptor?",
               f"Very uneven: {top2[1]} and {top2[0]} hold {share:.0%} of all {z.sum():,} z-scored entries.")
    big = S[z >= 20].sort_values("n_entries_with_zscores", ascending=False)
    small = S[(z < 20) & (t > 0)].sort_values(["n_entries_with_zscores", "n_entries_total"], ascending=False)
    none = S[t == 0].index.tolist()
    for rect, D, head in [([0.10, 0.14, 0.36, 0.58], big, f"≥ 20 z-scored entries  ({len(big)} receptors)"),
                          ([0.60, 0.14, 0.30, 0.58], small, f"< 20 z-scored entries  ({len(small)} receptors, zoomed scale)")]:
        ax, y = hbar_axes(fig, rect, D.index.tolist())
        ax.barh(y, D["n_entries_total"], color=PALE, edgecolor=LIGHT, linewidth=1, height=0.72)
        ax.barh(y, D["n_entries_with_zscores"], color=DARK, height=0.72)
        for yi, a, b in zip(y, D["n_entries_with_zscores"], D["n_entries_total"]):
            ax.text(b, yi, f"  {a} / {b}", va="center", fontsize=14 if len(D) > 20 else 16, color=INK)
        ax.set_xlim(0, D["n_entries_total"].max() * 1.28)
        fig.text(rect[0] - 0.06, 0.775, head, fontsize=19, fontweight="bold", color=INK)
    # direct key
    for x, lab, kw in [(0.04, "z-scored", dict(color=DARK)), (0.12, "all entries found", dict(facecolor=PALE, edgecolor=LIGHT))]:
        fig.add_artist(Rectangle((x, 0.733), 0.012, 0.018, transform=fig.transFigure, **kw))
        fig.text(x + 0.016, 0.742, lab, fontsize=15, color=INK, va="center")
    few = S[(t > 0) & (z <= 3)].index.tolist()
    frac = (z / t.replace(0, np.nan))
    lo = frac[(z >= 20)].sort_values()
    notes = [f"{len(few)} receptors have only 0–3 z-scored entries.",
             f"Share of found entries that are z-scored differs: {lo.index[0]} {lo.iloc[0]:.0%}, "
             f"{lo.index[1]} {lo.iloc[1]:.0%}, but {lo.index[-1]} {lo.iloc[-1]:.0%}.",
             f"No PDB entry found: {', '.join(none)}.  Dataset composition, not a biological result."]
    fig.text(0.04, 0.025, "     ·     ".join(notes[:2]) + "\n" + notes[2], fontsize=15, color=MUTED, linespacing=1.5)
    for n in S.index:
        rec.append(dict(page=1, receptor=n, n_entries_total=t[n], n_entries_with_zscores=z[n]))
    return fig


def p2_apo_holo(S, rec):
    zs = S[S["n_entries_with_zscores"] > 0]
    mixed = zs[(zs["n_apo"] > 0) & (zs["n_holo"] > 0)].sort_values(["n_apo", "n_holo"], ascending=False)
    holo = zs[(zs["n_apo"] == 0) & (zs["n_holo"] > 0)].sort_values("n_holo", ascending=False)
    apo = zs[(zs["n_apo"] > 0) & (zs["n_holo"] == 0)].sort_values("n_apo", ascending=False)
    neither = S[~S.index.isin(mixed.index.union(holo.index).union(apo.index))]
    assert len(mixed) + len(holo) + len(apo) + len(neither) == 48
    hd = int((mixed["n_holo"] / mixed["n_entries_with_zscores"] >= 0.8).sum())
    fig = page("Which receptors support an APO vs HOLO comparison?",
               f"Only {len(mixed)} receptors have both APO and HOLO entries, and {hd} of them are ≥ 80% HOLO.")
    ax, y = hbar_axes(fig, [0.11, 0.09, 0.40, 0.62], mixed.index.tolist())
    n = mixed["n_entries_with_zscores"]
    left = np.zeros(len(mixed))
    for k, c, h in [("n_holo", DARK, None), ("n_apo", MID, None), ("n_ambiguous", LIGHT, None), ("n_trace_lbd", BG, "///")]:
        w = mixed[k] / n * 100
        ax.barh(y, w, left=left, color=c, height=0.72, hatch=h, edgecolor=FAINT if h else c, linewidth=0.8 if h else 0)
        left += w
    for yi, (nm, r) in zip(y, mixed.iterrows()):
        tier = {"PRIMARY": DARK, "EXPLORATORY": BG}.get(r["apo_holo_tier"], "")
        ax.text(104, yi, f"  {r['n_apo']} APO · {r['n_holo']} HOLO" + (f" · {r['n_trace_lbd']} TRACE" if r["n_trace_lbd"] else ""),
                va="center", fontsize=14, color=INK)
        if tier:
            ax.plot(103, yi, "o", ms=9, mfc=tier, mec=DARK, mew=1.5, clip_on=False)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100], ["0%", "50%", "100%"])
    fig.text(0.05, 0.775, f"Both APO and HOLO  ({len(mixed)} receptors, % of z-scored entries)",
             fontsize=19, fontweight="bold", color=INK)
    for i, (lab, c, h) in enumerate([("HOLO", DARK, None), ("APO", MID, None), ("ambiguous", LIGHT, None), ("TRACE_LBD fragment", BG, "///")]):
        x = [0.05, 0.11, 0.16, 0.235][i]
        fig.add_artist(Rectangle((x, 0.738), 0.012, 0.018, facecolor=c, hatch=h, edgecolor=FAINT if h else c,
                                 transform=fig.transFigure))
        fig.text(x + 0.015, 0.747, lab, fontsize=14, color=INK, va="center")
    # side lists
    hl = [f"{nm}  {r['n_holo']}" + (f"  (+{r['n_trace_lbd']} TRACE)" if r["n_trace_lbd"] else "") for nm, r in holo.iterrows()]
    half = (len(hl) + 1) // 2
    names_block(fig, 0.73, 0.775, f"HOLO only  ({len(holo)})", hl[:half])
    fig.text(0.855, 0.733, "\n".join(hl[half:]), fontsize=16, color=MUTED, va="top", linespacing=1.55)
    names_block(fig, 0.73, 0.40, f"APO only  ({len(apo)})", [f"{nm}  {r['n_apo']}" for nm, r in apo.iterrows()][:3])
    fig.text(0.855, 0.358, "\n".join(f"{nm}  {r['n_apo']}" for nm, r in list(apo.iterrows())[3:]),
             fontsize=16, color=MUTED, va="top", linespacing=1.55)
    tr = neither[neither["n_trace_lbd"] > 0]
    nz = neither[(neither["n_trace_lbd"] == 0) & (neither["n_entries_total"] > 0)].index.tolist()
    nf = neither[neither["n_entries_total"] == 0].index.tolist()
    names_block(fig, 0.73, 0.20, f"No APO or HOLO entry  ({len(neither)})",
                [f"{nm}: {r['n_trace_lbd']} TRACE_LBD fragments only" for nm, r in tr.iterrows()]
                + [f"{', '.join(nz)}: no z-scored entry", f"{', '.join(nf)}: no entry found"])
    for x, fc, lab in [(0.36, DARK, "in the APO/HOLO analysis"), (0.50, BG, "exploratory")]:
        fig.add_artist(matplotlib.lines.Line2D([x], [0.747], marker="o", ms=9, mfc=fc, mec=DARK, mew=1.5,
                                               transform=fig.transFigure))
        fig.text(x + 0.008, 0.747, lab, fontsize=14, color=INK, va="center")
    fig.text(0.04, 0.025, "TRACE_LBD fragments (hatched) are short receptor peptides; they are never counted as APO or HOLO.",
             fontsize=15, color=MUTED)
    for nm, r in S.iterrows():
        g = ("both" if nm in mixed.index else "holo_only" if nm in holo.index else "apo_only" if nm in apo.index else "neither")
        rec.append(dict(page=2, receptor=nm, apo_holo_group=g, n_apo=r["n_apo"], n_holo=r["n_holo"], n_ambiguous=r["n_ambiguous"],
                        n_trace_lbd=r["n_trace_lbd"], n_entries_with_zscores=r["n_entries_with_zscores"],
                        apo_holo_tier=r["apo_holo_tier"]))
    return fig


def p3_ligands(S, rec):
    L = S[S["n_unique_ligands"] >= 5].sort_values("n_unique_ligands", ascending=False)
    few = S[(S["n_unique_ligands"] > 0) & (S["n_unique_ligands"] < 5)].sort_values("n_unique_ligands", ascending=False)
    zero = S[(S["n_unique_ligands"] == 0) & (S["n_entries_with_zscores"] > 0)].index.tolist()
    TOP = 50
    top = L[L["n_unique_ligands"] >= TOP].index.tolist()
    fig = page("Ligand diversity is very uneven",
               "Structure count and ligand diversity are related, but they are not the same thing.")
    ax, y = hbar_axes(fig, [0.10, 0.12, 0.60, 0.62], L.index.tolist(), bold=top)
    ax.barh(y, L["n_unique_ligands"], color=[DARK if n in top else MID for n in L.index], height=0.72)
    for yi, (nm, r) in zip(y, L.iterrows()):
        ax.text(r["n_unique_ligands"], yi, f"  {r['n_unique_ligands']}", va="center", fontsize=14, color=INK,
                fontweight="bold" if nm in top else "normal")
        ax.text(r["n_unique_ligands"] + L["n_unique_ligands"].max() * 0.055, yi, f"  in {r['n_entries_with_zscores']} entries",
                va="center", fontsize=12, color=FAINT)
    ax.set_xlim(0, L["n_unique_ligands"].max() * 1.22)
    ax.set_xlabel("unique functional ligands (CCD IDs) in z-scored entries", fontsize=15, color=MUTED)
    fig.text(0.05, 0.775, f"{len(L)} receptors with ≥ 5 unique ligands", fontsize=19, fontweight="bold", color=INK)
    names_block(fig, 0.78, 0.74, f"≥ {TOP} unique ligands  ({len(top)})", [", ".join(top[i:i + 3]) for i in range(0, len(top), 3)])
    fi = [f"{nm} {r['n_unique_ligands']}" for nm, r in few.iterrows()]
    names_block(fig, 0.78, 0.52, f"1–4 unique ligands  ({len(few)})", [", ".join(fi[i:i + 3]) for i in range(0, len(fi), 3)])
    names_block(fig, 0.78, 0.31, f"No bound functional ligand  ({len(zero)})", [", ".join(zero[i:i + 4]) for i in range(0, len(zero), 4)])
    fig.text(0.78, 0.14, "Receptors without z-scored entries\nare not shown.", fontsize=15, color=FAINT, va="top")
    for nm, r in S.iterrows():
        rec.append(dict(page=3, receptor=nm, n_unique_ligands=r["n_unique_ligands"], n_entries_with_zscores=r["n_entries_with_zscores"]))
    return fig


def p4_partners(S, rec):
    MIN_N = 5
    P = S[S["n_entries_with_zscores"] >= MIN_N].copy()
    n = P["n_entries_with_zscores"]
    P["nr"] = P["n_partner_nr_primary"] / n * 100
    P["co"] = P["n_partner_coregulator_peptide_primary"] / n * 100
    P["ot"] = (P["n_with_partner"] - P["n_partner_nr_primary"] - P["n_partner_coregulator_peptide_primary"]) / n * 100
    assert (P["ot"] >= -1e-9).all() and np.allclose(P[["nr", "co", "ot"]].sum(axis=1), P["frac_with_partner"] * 100, atol=0.06)
    P = P.sort_values(["frac_with_partner", "n_entries_with_zscores"], ascending=[False, False])
    rich = P[P["frac_with_partner"] >= 0.9].index.tolist()
    alone = P[P["frac_with_partner"] == 0].index.tolist()
    fig = page("Are structures mostly standalone or partner-bound?",
               "Some receptor datasets are almost always partner-bound; others are entirely standalone.")
    ax, y = hbar_axes(fig, [0.10, 0.13, 0.56, 0.585], P.index.tolist(), bold=rich + alone)
    left = np.zeros(len(P))
    for k, c in [("nr", DARK), ("co", MID), ("ot", LIGHT)]:
        ax.barh(y, P[k], left=left, color=c, height=0.72)
        left += P[k].to_numpy()
    for yi, (nm, r) in zip(y, P.iterrows()):
        ax.text(r["frac_with_partner"] * 100, yi, f"  {r['frac_with_partner']:.0%}", va="center", fontsize=14, color=INK)
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("% of z-scored entries with a partner chain contacting the LBD (5 Å)", fontsize=15, color=MUTED)
    fig.text(0.05, 0.775, f"{len(P)} receptors with ≥ {MIN_N} z-scored entries", fontsize=19, fontweight="bold", color=INK)
    for i, (lab, c) in enumerate([("NR / heterodimer partner", DARK), ("coregulator peptide", MID), ("other partner", LIGHT)]):
        x = 0.05 + i * 0.13
        fig.add_artist(Rectangle((x, 0.738), 0.012, 0.018, color=c, transform=fig.transFigure))
        fig.text(x + 0.015, 0.747, lab, fontsize=14, color=INK, va="center")
    names_block(fig, 0.73, 0.70, f"≥ 90% partner-bound  ({len(rich)})", [", ".join(rich[i:i + 3]) for i in range(0, len(rich), 3)])
    names_block(fig, 0.73, 0.50, f"Fully standalone (0%)  ({len(alone)})", [", ".join(alone[i:i + 3]) for i in range(0, len(alone), 3)])
    few = S[(S["n_entries_with_zscores"] > 0) & (S["n_entries_with_zscores"] < MIN_N)]
    fig.text(0.73, 0.30, f"Not ranked (< {MIN_N} z-scored entries):", fontsize=15, color=MUTED, va="top", fontweight="bold")
    items = [f"{nm} {r['n_with_partner']}/{r['n_entries_with_zscores']}" for nm, r in few.iterrows()]
    fig.text(0.73, 0.265, "\n".join(", ".join(items[i:i + 4]) for i in range(0, len(items), 4)),
             fontsize=14, color=FAINT, va="top", linespacing=1.5)
    fig.text(0.04, 0.025, "One class per entry: NR partner first, then coregulator peptide, then any other partner chain.",
             fontsize=15, color=MUTED)
    for nm, r in S.iterrows():
        rec.append(dict(page=4, receptor=nm, ranked=nm in P.index, frac_with_partner=r["frac_with_partner"],
                        n_with_partner=r["n_with_partner"], n_partner_nr_primary=r["n_partner_nr_primary"],
                        n_partner_coregulator_peptide_primary=r["n_partner_coregulator_peptide_primary"],
                        n_entries_with_zscores=r["n_entries_with_zscores"]))
    return fig


def p5_iqr(S, rec):
    A = S[S["heterogeneity_class"] != "not assessable"].sort_values("median_iqr_z", ascending=False)
    lo_cut, hi_cut = S["heterogeneity_cut_low_moderate"].iloc[0], S["heterogeneity_cut_moderate_high"].iloc[0]
    excl = S[S["heterogeneity_class"] == "not assessable"]
    col = {"high": DARK, "moderate": MID, "low": LIGHT}
    fig = page("Between-PDB flexibility heterogeneity differs by receptor",
               "Median IQR(z) of each canonical position across PDB entries, highest to lowest.")
    ax, y = hbar_axes(fig, [0.10, 0.12, 0.56, 0.60], A.index.tolist(),
                      bold=A.index[A["heterogeneity_class"] != "moderate"].tolist())
    ax.barh(y, A["median_iqr_z"], color=[col[c] for c in A["heterogeneity_class"]], height=0.72)
    for yi, v in zip(y, A["median_iqr_z"]):
        ax.text(v, yi, f"  {v:.2f}", va="center", fontsize=14, color=INK)
    ax.set_xlim(0, A["median_iqr_z"].max() * 1.15)
    ax.set_xlabel("median position IQR(z) across PDB entries", fontsize=15, color=MUTED)
    for cl in ["high", "moderate", "low"]:
        yy = y[(A["heterogeneity_class"] == cl).to_numpy()]
        ax.text(A["median_iqr_z"].max() * 1.15, yy.mean(), cl.upper(), ha="right", va="center", fontsize=17,
                fontweight="bold", color=col[cl] if cl != "low" else MID)
    fig.text(0.05, 0.775, f"{len(A)} receptors where the IQR is assessable", fontsize=19, fontweight="bold", color=INK)
    fig.text(0.05, 0.742, f"tertiles among them:  low ≤ {lo_cut:.3f} < moderate ≤ {hi_cut:.3f} < high",
             fontsize=14, color=MUTED, va="center")
    names_block(fig, 0.72, 0.70, "What this is", ["Spread of per-entry B-factor z between", "crystal structures of one receptor.",
                                                 "Not motion in solution, and not",
                                                 "controlled for burial / surface exposure."])
    nm_ex = excl[excl["n_entries_with_zscores"] > 0].index.tolist()
    names_block(fig, 0.72, 0.42, f"Not assessable  ({len(excl)})",
                ["< 20 positions seen in ≥ 5 entries:"] + [", ".join(nm_ex[i:i + 4]) for i in range(0, len(nm_ex), 4)]
                + ["no z-scored entry:", ", ".join(excl.index[excl["n_entries_with_zscores"] == 0][:4]),
                   ", ".join(excl.index[excl["n_entries_with_zscores"] == 0][4:])])
    for nm, r in S.iterrows():
        rec.append(dict(page=5, receptor=nm, median_iqr_z=r["median_iqr_z"], heterogeneity_class=r["heterogeneity_class"]))
    return fig


def p6_contrasts(S, rec):
    pairs = [("ERα", "PPARγ"), ("LXRα", "LXRβ"), ("VDR", "FXR"), ("GR", "PXR")]
    g = lambda n, k: S.loc[n, k]  # noqa: E731
    iqr = lambda n: f"{g(n, 'median_iqr_z'):.2f}"  # noqa: E731
    # each observation is checked against the table before it is drawn
    assert min(g(n, "n_entries_with_zscores") for n in pairs[0]) >= 300 and min(g(n, "n_unique_ligands") for n in pairs[0]) >= 200
    assert g("PPARγ", "median_iqr_z") > g("ERα", "median_iqr_z")
    assert all(g(n, "n_apo") == 0 and g(n, "n_holo") > 0 for n in pairs[1])
    assert g("VDR", "frac_with_partner") == 0 and g("FXR", "frac_with_partner") >= 0.9
    assert g("VDR", "heterogeneity_class") == "low" and g("FXR", "heterogeneity_class") == "high"
    assert g("GR", "heterogeneity_class") == "high" and g("PXR", "heterogeneity_class") == "low"
    obs = {
        0: f"Both data-rich and ligand-diverse; PPARγ has the higher median IQR ({iqr('PPARγ')} vs {iqr('ERα')}).",
        1: f"Closely related, both HOLO-only; median IQR {iqr('LXRα')} vs {iqr('LXRβ')} "
           f"(LXRα: only {g('LXRα', 'n_entries_with_zscores')} entries).",
        2: "Both rich in structures and ligands; VDR standalone, low IQR; FXR partner-bound, high IQR.",
        3: f"Both ligand-diverse; GR has high IQR ({iqr('GR')}), PXR low ({iqr('PXR')}).",
    }
    metrics = [("z-scored entries", "n_entries_with_zscores", "{:.0f}", S["n_entries_with_zscores"].max()),
               ("unique ligands", "n_unique_ligands", "{:.0f}", S["n_unique_ligands"].max()),
               ("with partner", "frac_with_partner", "{:.0%}", 1.0),
               ("median IQR(z)", "median_iqr_z", "{:.2f}", S["median_iqr_z"].max())]
    fig = page("What do these variables tell us together?",
               "Dataset size, ligand diversity and partner context do not map simply onto between-PDB heterogeneity.")
    cells = [(0.04, 0.43), (0.52, 0.43), (0.04, 0.06), (0.52, 0.06)]
    for i, ((a, b), (x0, y0)) in enumerate(zip(pairs, cells)):
        W, H = 0.44, 0.34
        fig.add_artist(Rectangle((x0, y0), W, H, facecolor=PALE, edgecolor="none", transform=fig.transFigure))
        fig.text(x0 + 0.015, y0 + H - 0.045, f"{a}  vs  {b}", fontsize=24, fontweight="bold", color=INK)
        fig.text(x0 + 0.015, y0 + 0.022, obs[i], fontsize=13.5, color=INK, wrap=False)
        for j, (lab, k, fmt, mx) in enumerate(metrics):
            yy = y0 + H - 0.105 - j * 0.052
            fig.text(x0 + 0.015, yy, lab, fontsize=15, color=MUTED, va="center")
            for kk, (nm, c) in enumerate([(a, DARK), (b, MID)]):
                v = g(nm, k)
                bx, by, bw = x0 + 0.16, yy + (0.012 if kk == 0 else -0.012), 0.20
                fig.add_artist(Rectangle((bx, by - 0.009), bw * (0 if pd.isna(v) else v / mx), 0.018, color=c,
                                         transform=fig.transFigure))
                cl = f"  ({g(nm, 'heterogeneity_class')})" if k == "median_iqr_z" else ""
                fig.text(bx + bw * (0 if pd.isna(v) else v / mx) + 0.004, by, f"{nm} {fmt.format(v)}{cl}",
                         fontsize=13, color=INK, va="center")
        for nm in (a, b):
            rec.append(dict(page=6, receptor=nm, pair=f"{a} vs {b}", **{k: g(nm, k) for _, k, _, _ in metrics},
                            n_apo=g(nm, "n_apo"), n_holo=g(nm, "n_holo")))
    fig.text(0.04, 0.025, "Associations only, no causal claim: each pair is a candidate for receptor-specific follow-up.",
             fontsize=15, color=MUTED)
    return fig


PAGES = [("1_entries", p1_entries), ("2_apo_holo", p2_apo_holo), ("3_ligands", p3_ligands),
         ("4_partners", p4_partners), ("5_iqr", p5_iqr), ("6_contrasts", p6_contrasts)]


def main():
    S = load()
    (FIG / "chart_only").mkdir(parents=True, exist_ok=True)
    rec: list[dict] = []
    with PdfPages(OUT / "catalogue_digest.pdf") as pdf:
        for name, fn in PAGES:
            fig = fn(S, rec)
            png = FIG / f"digest_{name}.png"
            fig.savefig(png, dpi=100)
            pdf.savefig(fig)
            plt.close(fig)
            im = Image.open(png)
            im.crop((0, int(round((1 - HEADER_Y) * im.height)), im.width, im.height)).save(FIG / "chart_only" / png.name)
    pd.DataFrame(rec).to_csv(OUT / "digest_values.csv", index=False)
    pd.DataFrame([dict(page=i + 1, figure=f"digest_{n}.png", **t) for i, ((n, _), t) in enumerate(zip(PAGES, TEXT))]
                 ).to_csv(OUT / "digest_page_text.csv", index=False)
    print(f"{len(PAGES)} pages -> {OUT}")


if __name__ == "__main__":
    main()
