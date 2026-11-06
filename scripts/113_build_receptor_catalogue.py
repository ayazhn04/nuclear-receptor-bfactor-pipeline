"""
Receptor catalogue (v2), step 2: one summary row and one 1920x1080 page per receptor,
two overview figures, and a PDF of all pages. Reads the corrected v2 outputs and the
catalogue entry table from scripts/112; changes nothing in the analysis pipeline.

Denominator of every composition panel: the receptor's entries that contribute
z-scores (the entries the v2 flexibility profiles are built from).

Per-position flexibility: the committed v2 receptor-position profile (median and IQR of
the per-entry z across PDB entries). A position is "supported" when >= 5 entries observe
it (the project's rule for IQR); positions seen in 1-4 entries are drawn faint and never
enter the IQR summaries. Heterogeneity class = tertile of the receptor's median IQR among
receptors with >= 20 supported positions (cut-offs written to the summary table).

Writes reports/catalogue_v2/:
  tables/receptor_catalogue_summary_v2.csv      one row per receptor (all page numbers)
  tables/catalogue_rules_v2.json                every threshold/rule used for labels and notes
  figures/receptor_<NRcode>.png                 one page per receptor (1920x1080)
  figures/overview_*.png                        dataset-at-a-glance and heterogeneity overview
  receptor_catalogue_v2.pdf                     dataset at a glance (2 pages), then the 48 receptor pages
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402

V2 = PROJECT_ROOT / "data" / "rebuild_v2"
V2T = PROJECT_ROOT / "reports" / "rebuild_v2" / "tables"
OUT = PROJECT_ROOT / "reports" / "catalogue_v2"
TAB, FIG = OUT / "tables", OUT / "figures"

MIN_ENTRIES = 5          # project rule: a position needs >= 5 entries for IQR
MIN_POS_FOR_CLASS = 20   # receptors with fewer supported positions get no heterogeneity class
Y_LIM = (-2.5, 5.0)      # fixed z axis on every page so receptors compare at a glance
RULES = {
    "denominator": "entries contributing z-scores (v2 profiles)",
    "zscored_label": "'<z-scored> / <total> entries z-scored (%)' beside the z-plot axis title; total = every "
                     "receptor x PDB entry found at discovery, z-scored = >= 1 finite b_factor_z in the v2 profiles; "
                     "TRACE_LBD-only entries that are z-scored are named ('incl. N TRACE_LBD')",
    "data_level": {"none": "0 entries", "data-poor": "1-4", "moderate data": "5-19", "data-rich": ">= 20"},
    "composition": {"holo only": "no APO, >= 1 HOLO", "mostly holo": "HOLO >= 90%", "mixed apo/holo": ">= 3 APO and >= 3 HOLO",
                    "mostly apo": "APO >= 60%", "fragment-heavy": "TRACE_LBD-only >= 50% of entries"},
    "partner": {"partner-rich": ">= 50% of entries with a contacting partner chain", "mostly standalone": "<= 10%",
                "mixed partner context": "otherwise"},
    "ligand": {"one dominant ligand": "most common functional ligand in >= 50% of HOLO entries",
               "highly diverse": "unique ligands >= 60% of HOLO entries with a ligand", "mixed": "otherwise"},
    "construct_diversity_high": "distinct construct boundary patterns >= 50% of entries (>= 10 entries); "
                                "construct = first/last UniProt residue of the entity's sequence alignment",
    "point_mutants_common": ">= 30% of entries carry >= 1 point substitution vs the reference sequence",
    "resolution_broad": "IQR >= 0.8 A or worst entry >= 3.5 A",
    "low_coverage": "< 50% of observed positions are supported (>= 5 entries)",
    "strong_holo": "HOLO >= 90% with >= 10 entries",
    "heterogeneity_class": f"tertile of median IQR among receptors with >= {MIN_POS_FOR_CLASS} supported positions",
    "construct_line": "page shows engineered mutants when >= 20% of entries, fusion/chimera when >= 5% (or < 10 entries), "
                      "internal deletions when >= 10% (or < 10 entries); every count stays in the summary table",
    "caveats": "at most 3, receptor-specific: no data / no professor map, TRACE_LBD fragments, weak coverage, "
               "internal deletions >= 50%, fusion >= 20%, ligand dominance >= 70%, crystal forms >= 50% of entries",
    "substitution_classes": "ENGINEERED / NATURAL_VARIANT / SEQUENCE_CONFLICT / DECLARED_OTHER from the depositor's "
                            "_struct_ref_seq_dif (or a UniProt natural variant); UNDECLARED = no annotation; "
                            "MODIFIED_RESIDUE = non-standard component on a reference residue",
    "crystal_form": "same space group and unit cell within 2% (lengths) and 2 degrees (angles)",
    "engineered_mutant": "entry with >= 1 residue the depositor marks 'engineered mutation' (_struct_ref_seq_dif, "
                         "reference residue identical to ours) or the entity's _entity.pdbx_mutation lists at the same "
                         "residue; UNCERTAIN = only conflict/unannotated differences, or a declared substitution that no "
                         "observed residue confirms; see catalogue_mutation_audit_v2.csv",
    "partner_identity": "normalized gene-level identity (scripts/utils/partner_identity.py); DNA and same-receptor "
                        "chains excluded; unidentified peptides count once",
    "contexts": "per entry with an LBD context (TRACE_LBD excluded): ligand state x partner category, and functional "
                "ligand CCD id(s) x partner identity(ies); several are joined in sorted order (deterministic)",
}

# palette (matches the report figures: slate on white)
INK, MUTED, GRID, BG = "#1d2533", "#6b7280", "#e6e9ee", "#fbfbf9"
SLATE, SLATE_LT, PALE = "#4a5568", "#a3b1c6", "#dfe5ee"
HOLO_C, APO_C, AMB_C, TRACE_C = "#3c4a5e", "#8ea4c2", "#d6dce5", "#e3a857"
POCKET_C, AF2_C, GROOVE_C, LINE5_C = "#2563eb", "#d97706", "#7c3aed", "#c53030"
FONT = "Helvetica Neue"
MONO = "Menlo"
plt.rcParams.update({"font.family": FONT, "font.size": 13, "axes.edgecolor": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG})

PARTNER_ORDER = [("COREGULATOR_PEPTIDE", "Coregulator peptide"), ("NR_PARTNER", "NR / heterodimer partner"),
                 ("SAME_RECEPTOR_CHAIN", "Same-receptor chain"), ("OTHER_PEPTIDE", "Other peptide"),
                 ("OTHER_PROTEIN", "Other protein"), ("NUCLEIC_ACID", "DNA / RNA")]


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------


def load() -> dict:
    d = {}
    d["meta"] = pd.DataFrame(NR_METADATA)[["uniprot_id", "nr_code", "common_name", "group"]]
    d["entries"] = pd.read_csv(TAB / "catalogue_entry_context_v2.csv", dtype={"pdb_id": str}, keep_default_na=False)
    for c in ["has_coordinate_data", "contributes_z", "trace_lbd_only", "in_primary_qc_set", "fusion_or_chimera"]:
        d["entries"][c] = d["entries"][c].astype(str) == "True"
    d["contacts"] = pd.read_csv(TAB / "catalogue_partner_contacts_v2.csv", dtype={"pdb_id": str}, keep_default_na=False)
    d["seqdiff"] = pd.read_csv(TAB / "catalogue_sequence_differences_v2.csv", dtype={"pdb_id": str}, keep_default_na=False)
    d["audit"] = pd.read_csv(TAB / "catalogue_mutation_audit_v2.csv", dtype={"pdb_id": str}, keep_default_na=False)
    d["profile"] = pd.read_csv(V2 / "analysis" / "receptor_position_professor_profile_v2.csv", low_memory=False)
    d["span"] = pd.read_csv(V2 / "inputs" / "analysis_span_v2.csv")
    d["segmap"] = pd.read_csv(V2 / "analysis" / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    d["sites"] = pd.read_parquet(V2 / "analysis" / "current_functional_site_definitions.parquet")
    d["elig"] = pd.read_csv(V2T / "apo_holo_receptor_eligibility.csv")
    d["cohort"] = pd.read_csv(V2 / "analysis" / "modeling_receptor_cohort.csv")
    d["instances"] = pd.read_csv(V2 / "analysis" / "instance_bfactor_summary.csv", dtype={"pdb_id": str},
                                 usecols=["uniprot_id", "pdb_id", "n_ca_usable", "normalization_status"])
    return d


def zscore_inventory(E: pd.DataFrame, inst: pd.DataFrame) -> dict:
    """Entry-level z-score accounting. total = every receptor x PDB entry found (discovery);
    z-scored = entries with >= 1 finite b_factor_z in the v2 profiles. Every entry not
    z-scored gets exactly one reason: no coordinate data (with the discovery reason), or
    coordinate data whose instances could not be normalized (classified from the v2
    instance normalization_status and usable CA count)."""
    n_total = len(E)
    n_coord = int(E["has_coordinate_data"].sum())
    n_z = int(E["contributes_z"].sum())
    assert not (E["contributes_z"] & ~E["has_coordinate_data"]).any()
    out = {"n_entries_total": n_total, "n_entries_with_coordinates": n_coord, "n_entries_with_zscores": n_z,
           "n_entries_without_zscores": n_total - n_z,
           "fraction_entries_zscored": round(n_z / n_total, 4) if n_total else np.nan}
    noc = E.loc[~E["has_coordinate_data"], "no_coordinate_reason"]
    out["n_not_zscored_receptor_chain_outside_lbd"] = int((noc == "receptor chain outside the LBD").sum())
    out["n_not_zscored_not_xray"] = int((noc == "not X-ray (NMR/EM)").sum())
    out["n_not_zscored_no_coordinates_other"] = len(noc) - out["n_not_zscored_receptor_chain_outside_lbd"] \
        - out["n_not_zscored_not_xray"]
    reasons = Counter()
    for pid in E.loc[E["has_coordinate_data"] & ~E["contributes_z"], "pdb_id"]:
        I = inst[inst["pdb_id"] == pid]
        if I.empty:
            reasons["other"] += 1
        elif (I["n_ca_usable"] == 0).all():
            reasons["no_usable_ca"] += 1
        elif (I["normalization_status"] == "NOT_NORMALIZABLE_FEWER_THAN_2_USABLE_B").all():
            reasons["too_few_values"] += 1
        elif I["normalization_status"].isin(["NOT_NORMALIZABLE_ZERO_SD", "NOT_NORMALIZABLE_FEWER_THAN_2_USABLE_B"]).all():
            reasons["zero_sd"] += 1
        else:
            reasons["other"] += 1
    for k in ["no_usable_ca", "too_few_values", "zero_sd", "other"]:
        out[f"n_not_zscored_{k}"] = reasons[k]
    # n_trace_lbd (set later) stays the composition count = TRACE_LBD among z-scored entries;
    # n_trace_lbd_only (inventory) = all TRACE_LBD-only entries found
    out["n_trace_lbd_not_zscored"] = int((E["trace_lbd_only"] & ~E["contributes_z"]).sum())
    out["n_trace_lbd_zscored"] = int((E["trace_lbd_only"] & E["contributes_z"]).sum())
    out["n_zscored_lbd_context"] = n_z - out["n_trace_lbd_zscored"]
    out["zscored_label"] = zscored_label(out)
    assert out["n_entries_without_zscores"] == sum(v for k, v in out.items() if k.startswith("n_not_zscored_"))
    return out


def zscored_label(z: dict) -> str:
    """The small label beside the z-plot axis title (denominator = all entries found)."""
    t, n = z["n_entries_total"], z["n_entries_with_zscores"]
    if t == 0:
        return "no PDB entries found"
    pct = f"{n / t:.0%}"
    if n < t and pct == "100%":
        pct = ">99%"   # never show 100% when an entry is missing
    s = f"{n} / {t} entr{'y' if t == 1 else 'ies'} z-scored ({pct})"
    if z["n_trace_lbd_zscored"]:
        s += " · all TRACE_LBD" if z["n_trace_lbd_zscored"] == n else f" · incl. {z['n_trace_lbd_zscored']} TRACE_LBD"
    return s


def helix_blocks(seg: pd.DataFrame) -> list[dict]:
    """Contiguous runs of one professor element (an element with two intervals gives two
    blocks; a dual position belongs to both neighbouring blocks)."""
    h = seg[seg["structural_segment"] != ""].sort_values("canonical_uniprot_position")
    blocks = []
    for el, g in h.groupby("structural_segment", sort=False):
        pos = sorted(g["canonical_uniprot_position"].astype(int))
        review = set(g["assignment_confidence"]) == {"review_required"}
        start = prev = pos[0]
        for p in pos[1:] + [None]:
            if p is None or p != prev + 1:
                blocks.append({"label": el.replace("'", "′"), "start": start, "end": prev, "review": review})
                start = p
            prev = p if p is not None else prev
    return sorted(blocks, key=lambda b: b["start"])


def summarise(uid: str, d: dict, het_cuts: tuple | None) -> dict:
    m = d["meta"].set_index("uniprot_id").loc[uid]
    E = d["entries"][d["entries"]["uniprot_id"] == uid]
    C = E[E["contributes_z"]]
    n_c = len(C)
    sp = d["span"].set_index("uniprot_id").loc[uid]
    prof = d["profile"][d["profile"]["uniprot_id"] == uid]
    seg = d["segmap"][d["segmap"]["uniprot_id"] == uid]
    r = {"uniprot_id": uid, "nr_code": m["nr_code"], "common_name": m["common_name"], "group": m["group"],
         "reference_accession": sp["reference_accession"], "professor_map": bool(sp["professor_mapping"]),
         "analysis_span": f"{int(sp['analysis_span_start'])}-{int(sp['analysis_span_end'])}",
         "n_professor_elements": int(seg.loc[seg["structural_segment"] != "", "structural_segment"].nunique())}
    # inventory
    r.update({"n_entries_discovered": len(E), "n_entries_coordinate": int(E["has_coordinate_data"].sum()),
              "n_entries_contributing": n_c, "n_trace_lbd_only": int(E["trace_lbd_only"].sum()),
              "n_trace_lbd_only_contributing": int(C["trace_lbd_only"].sum()),
              "n_coordinate_not_normalizable": int((E["has_coordinate_data"] & ~E["contributes_z"]).sum()),
              "n_no_coordinate_data": int((~E["has_coordinate_data"]).sum()),
              "no_coordinate_reasons": "; ".join(f"{k} ({v})" for k, v in Counter(
                  x for v in E.loc[~E["has_coordinate_data"], "no_coordinate_reason"] for x in v.split(";") if x).items()),
              "n_entries_primary_qc": int(E["in_primary_qc_set"].sum()),
              "n_entries_two_receptor_entities": int((E["n_receptor_entities"].astype(int) > 1).sum())})
    inst = d["instances"][d["instances"]["uniprot_id"] == uid]
    r.update(zscore_inventory(E, inst))
    assert r["n_entries_with_zscores"] == n_c and r["n_entries_with_coordinates"] == r["n_entries_coordinate"]
    # apo/holo
    st = C["catalogue_state"].value_counts()
    for k in ["APO", "HOLO", "AMBIGUOUS", "TRACE_LBD"]:
        r[f"n_{k.lower()}"] = int(st.get(k, 0))
    # ligands
    holo = C[C["catalogue_state"] == "HOLO"]
    ligs = [x for v in holo["functional_ligands"] for x in v.split(";") if x]
    names = {}
    for ids, nm in zip(holo["functional_ligands"], holo["functional_ligand_names"]):
        for i, n in zip(ids.split(";"), nm.split(";")):
            names.setdefault(i, n)
    cnt = Counter(ligs)
    n_holo_lig = int((holo["functional_ligands"] != "").sum())
    top = cnt.most_common(5)
    r.update({"n_holo_with_ligand_id": n_holo_lig, "n_unique_ligands": len(cnt),
              "top_ligand": top[0][0] if top else "", "top_ligand_name": names.get(top[0][0], "") if top else "",
              "top_ligand_entries": top[0][1] if top else 0,
              "top_ligand_share_of_holo": round(top[0][1] / n_holo_lig, 3) if top and n_holo_lig else np.nan,
              "top5_ligands": "; ".join(f"{k} ({v})" for k, v in top)})
    if not n_holo_lig:
        r["ligand_context"] = "no ligand-bound entries"
    elif r["top_ligand_share_of_holo"] >= 0.5:
        r["ligand_context"] = "one dominant ligand"
    elif len(cnt) / n_holo_lig >= 0.6:
        r["ligand_context"] = "highly diverse"
    else:
        r["ligand_context"] = "mixed"
    # partners
    cls = [set(v.split(";")) if v else set() for v in C["partner_classes"]]
    r["n_standalone"] = sum(1 for c in cls if not c)
    r["n_with_partner"] = n_c - r["n_standalone"]
    for k, _ in PARTNER_ORDER:
        r[f"n_partner_{k.lower()}"] = sum(1 for c in cls if k in c)
    r["frac_with_partner"] = round(r["n_with_partner"] / n_c, 3) if n_c else np.nan
    # construct and experiment
    cons = C.dropna(subset=["construct_first_position"])
    cons = cons[cons["construct_first_position"] != ""]
    pairs = set(zip(cons["construct_first_position"], cons["construct_last_position"]))
    res = pd.to_numeric(C["resolution"], errors="coerce").dropna()
    subst = pd.to_numeric(C["n_point_substitutions"], errors="coerce").fillna(0)
    nseg = pd.to_numeric(C["construct_n_aligned_segments"], errors="coerce")
    ofr = pd.to_numeric(C["observed_span_fraction"], errors="coerce")
    r.update({"n_construct_patterns": len(pairs), "n_fusion_or_chimera": int(C["fusion_or_chimera"].sum()),
              "n_internal_deletion_constructs": int((nseg > 1).sum()),
              "median_observed_span_fraction": round(ofr.median(), 3) if len(ofr.dropna()) else np.nan,
              "n_entries_with_point_substitution": int((subst > 0).sum()),
              "n_non_human_or_ambiguous_identity": int((C["identity_class"] != "HUMAN_OR_HUMAN_DERIVED").sum()),
              "resolution_median": round(res.median(), 2) if len(res) else np.nan,
              "resolution_q25": round(res.quantile(.25), 2) if len(res) else np.nan,
              "resolution_q75": round(res.quantile(.75), 2) if len(res) else np.nan,
              "resolution_min": round(res.min(), 2) if len(res) else np.nan,
              "resolution_max": round(res.max(), 2) if len(res) else np.nan})
    # sequence differences, crystal forms, partner single-class assignment
    for c in ["engineered", "natural_variant", "sequence_conflict", "declared_other", "undeclared", "modified_residue"]:
        col = pd.to_numeric(C.get(f"n_diff_{c}", pd.Series(dtype=float)), errors="coerce").fillna(0)
        r[f"n_entries_{c}" if c != "engineered" else "n_entries_engineered"] = int((col > 0).sum())
    r["frac_entries_engineered"] = round(r["n_entries_engineered"] / n_c, 3) if n_c else np.nan
    sd = d["seqdiff"]
    eng = sd[(sd["uniprot_id"] == uid) & (sd["difference_class"] == "ENGINEERED") &
             sd["pdb_id"].isin(set(C["pdb_id"]))].drop_duplicates(["pdb_id", "label"])
    r["top_engineered_substitutions"] = "; ".join(f"{k} ({v})" for k, v in Counter(eng["label"]).most_common(5))
    r["n_crystal_forms"] = int(C["crystal_form_id"].replace("", np.nan).nunique())
    r["n_space_groups"] = int(C["space_group"].replace("", np.nan).nunique())
    sg = Counter(x for x in C["space_group"] if x)
    r["top_space_groups"] = "; ".join(f"{k} ({v})" for k, v in sg.most_common(3))
    yrs = pd.to_numeric(C["deposition_year"], errors="coerce").dropna()
    r["deposition_year_min"] = int(yrs.min()) if len(yrs) else np.nan
    r["deposition_year_max"] = int(yrs.max()) if len(yrs) else np.nan
    r["median_polymer_chains_asu"] = float(pd.to_numeric(C["n_polymer_chains_asu"], errors="coerce").median()) if n_c else np.nan
    r["median_receptor_copies_asu"] = float(pd.to_numeric(C["n_receptor_copies_asu"], errors="coerce").median()) if n_c else np.nan
    # one class per entry for the partner bar: NR partner > coregulator > other
    r["n_partner_nr_primary"] = sum(1 for c in cls if "NR_PARTNER" in c)
    r["n_partner_coregulator_peptide_primary"] = sum(1 for c in cls if "NR_PARTNER" not in c and "COREGULATOR_PEPTIDE" in c)
    ct = d["contacts"]
    ct = ct[(ct["uniprot_id"] == uid) & ct["pdb_id"].isin(set(C["pdb_id"]))].drop_duplicates(["pdb_id", "partner_entity_id"])
    def short(x: str, k: int = 42) -> str:
        x = " ".join(str(x).split()).lower()
        return x if len(x) <= k else x[:k].rsplit(" ", 1)[0] + "…"
    desc = Counter(short(x) for x in ct["partner_description"] if str(x).strip())
    r["top_partner"] = desc.most_common(1)[0][0] if desc else ""
    r["top_partners"] = "; ".join(f"{k} ({v})" for k, v in desc.most_common(5))
    # audited mutation status (catalogue_mutation_audit_v2.csv; only depositor-confirmed engineered
    # substitutions count) -- frac_entries_engineered above is kept for continuity and equals it
    au = d["audit"][d["audit"]["uniprot_id"] == uid]
    r["n_pdb_engineered_mutant"] = int((au["final_classification"] == "ENGINEERED_MUTANT").sum())
    r["n_pdb_substitution_uncertain"] = int((au["final_classification"] == "UNCERTAIN").sum())
    r["n_pdb_no_engineered_substitution"] = int((au["final_classification"] == "WT_NO_ENGINEERED_SUBSTITUTION").sum())
    r["engineered_mutant_fraction"] = round(r["n_pdb_engineered_mutant"] / len(au), 3) if len(au) else np.nan
    r["substitution_uncertain_fraction"] = round(r["n_pdb_substitution_uncertain"] / len(au), 3) if len(au) else np.nan
    r["mutation_classification_source"] = au["classification_source"].iloc[0] if len(au) else ""
    # partner identities and biological contexts (entries with an LBD context; TRACE_LBD excluded)
    X = C[C["ligand_partner_context"] != ""]
    ids = Counter(i for v in X["partner_identities"] for i in set(v.split(";")) if i)
    r["n_unique_partner_identities"] = len(ids)
    r["partner_identities"] = "; ".join(f"{k} ({v})" for k, v in ids.most_common())
    r["dominant_partner_identity"] = ids.most_common(1)[0][0] if ids else ""
    r["dominant_partner_identity_entries"] = ids.most_common(1)[0][1] if ids else 0
    lname = {}
    for i_, n_ in zip(C["functional_ligands"], C["functional_ligand_names"]):
        for a, b in zip(i_.split(";"), n_.split(";")):
            if a:
                lname.setdefault(a, b)
    for kind, col in [("ligand_partner", "ligand_partner_context"), ("state_partner", "state_partner_context")]:
        vc = X[col].value_counts()
        r[f"n_unique_{kind}_contexts"] = len(vc)
        r[f"dominant_{kind}_context"] = vc.index[0] if len(vc) else ""
        r[f"dominant_{kind}_context_fraction"] = round(vc.iloc[0] / len(X), 3) if len(vc) else np.nan
    r["dominant_ligand_partner_context_display"] = context_display(r["dominant_ligand_partner_context"], lname)
    r["n_context_entries"] = len(X)
    # flexibility heterogeneity
    span_len = int(sp["analysis_span_end"] - sp["analysis_span_start"] + 1)
    obs = prof[prof["n_pdb_entries_observed_at_position"] >= 1]
    sup = prof[prof["n_pdb_entries_observed_at_position"] >= MIN_ENTRIES]
    r.update({"n_span_positions": span_len, "n_observed_positions": len(obs), "n_supported_positions": len(sup),
              "median_iqr_z": round(sup["iqr_z"].median(), 3) if len(sup) else np.nan,
              "p90_iqr_z": round(sup["iqr_z"].quantile(.9), 3) if len(sup) else np.nan,
              "max_iqr_z": round(sup["iqr_z"].max(), 3) if len(sup) else np.nan,
              "max_iqr_position": int(sup.loc[sup["iqr_z"].idxmax(), "canonical_uniprot_position"]) if len(sup) else np.nan})
    if het_cuts and len(sup) >= MIN_POS_FOR_CLASS:
        r["heterogeneity_class"] = ("low" if r["median_iqr_z"] <= het_cuts[0] else
                                    "high" if r["median_iqr_z"] > het_cuts[1] else "moderate")
    else:
        r["heterogeneity_class"] = "not assessable"
    # where the upper tail sits (top-decile IQR positions by professor element)
    if len(sup):
        top = sup[sup["iqr_z"] >= sup["iqr_z"].quantile(.9)]
        where = Counter((e if isinstance(e, str) and e else "between helices") for e in top["canonical_elements"])
        r["top_decile_iqr_location"] = "; ".join(f"{k.replace('|', '/').replace(chr(39), '′')} {v}"
                                                 for k, v in where.most_common(4))
    else:
        r["top_decile_iqr_location"] = ""
    # analysis membership
    el = d["elig"].set_index("uniprot_id")
    r["apo_holo_tier"] = ("PRIMARY" if uid in el.index and el.loc[uid, "eligible_primary_5_5"] else
                          "EXPLORATORY" if uid in el.index and el.loc[uid, "eligible_exploratory_3_3"] else "not testable")
    co = d["cohort"].set_index("uniprot_id")
    r["in_modeling_cohort"] = bool(uid in co.index and co.loc[uid, "included_primary"])
    r["interpretation"] = interpretation(r)
    r["notes"] = " | ".join(notes(r))
    return r


def context_display(ctx: str, lname: dict) -> str:
    """'EST | NCOA2' -> 'estradiol + NCOA2'; ligand CCD ids shown by their (short) chemical name."""
    if not ctx:
        return ""
    lig, part = ctx.split(" | ")
    def short(i: str) -> str:
        nm = (lname.get(i) or i).lower()
        return i if len(nm) > 22 else nm
    lig = " + ".join(short(x) for x in lig.split(" + ")) if lig not in ("APO", "AMBIGUOUS", "UNKNOWN") else lig
    return f"{lig} + {part}"


def interpretation(r: dict) -> str:
    """One-line context summary (data amount is shown separately in the header)."""
    n = r["n_entries_contributing"]
    if n == 0:
        return "no usable LBD structures"
    apo, holo, trace = r["n_apo"], r["n_holo"], r["n_trace_lbd"]
    if trace / n >= 0.5:
        comp = "fragment-only" if trace == n else "fragment-heavy"
    elif holo == 0 and apo > 0:
        comp = "apo only"
    elif apo == 0 and holo > 0:
        comp = "holo only"
    elif holo / n >= 0.9:
        comp = "mostly holo"
    elif apo / n >= 0.6:
        comp = "mostly apo"
    else:
        comp = "mixed apo/holo"
    lig = {"no ligand-bound entries": "no bound ligands", "one dominant ligand": "one dominant ligand",
           "highly diverse": "ligand-diverse", "mixed": "mixed ligands"}[r["ligand_context"]]
    f = r["frac_with_partner"]
    part = "partner-rich" if f >= 0.5 else "mostly standalone" if f <= 0.1 else "partly partner-bound"
    het = ("IQR not assessable" if r["heterogeneity_class"] == "not assessable"
           else f"{r['heterogeneity_class']} IQR heterogeneity")
    return " · ".join([comp, lig, part, het])


def notes(r: dict) -> list[str]:
    """At most three receptor-specific caveats that could change interpretation. Facts already
    obvious from a page graphic (0 APO, IQR not assessable) are not repeated here."""
    n = r["n_entries_contributing"]
    out = []
    if not r["professor_map"]:
        out.append("No professor helix map: no canonical elements.")
    if n == 0:
        return out[:3]
    if r["n_trace_lbd"]:
        out.append(f"{r['n_trace_lbd']} TRACE_LBD fragment{'s' if r['n_trace_lbd'] > 1 else ''} "
                   "(receptor peptides) excluded from apo/holo.")
    if r["n_observed_positions"] and n >= MIN_ENTRIES and r["n_supported_positions"] / r["n_observed_positions"] < 0.5:
        out.append(f"Weak coverage: {r['n_supported_positions']} of {r['n_observed_positions']} observed positions "
                   "seen in ≥5 structures.")
    if r["n_internal_deletion_constructs"] / n >= 0.5:
        out.append(f"{r['n_internal_deletion_constructs']} of {n} constructs carry an internal deletion.")
    if r["n_fusion_or_chimera"] / n >= 0.2:
        out.append(f"Fusion/chimera constructs in {r['n_fusion_or_chimera'] / n:.0%} of structures.")
    if r["ligand_context"] == "one dominant ligand" and r["n_holo_with_ligand_id"] >= 5 and r["top_ligand_share_of_holo"] >= 0.7:
        out.append(f"Strong ligand dominance: {r['top_ligand']} in {r['top_ligand_share_of_holo']:.0%} of HOLO structures.")
    if n >= 10 and r["n_crystal_forms"] / n >= 0.5:
        out.append(f"Very heterogeneous crystal forms ({r['n_crystal_forms']} forms in {n} structures).")
    return out[:3]


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------

LABEL_C = "#5b6472"   # small-caps block labels


def block_label(fig, x: float, y: float, text: str) -> None:
    fig.text(x, y, text.upper(), fontsize=14, color=LABEL_C, fontweight="bold", va="bottom")


def stacked_bar(fig, rect, parts: list[tuple[str, int, str]], total: int, text_dark=("#3c4a5e", "#4a5568")) -> None:
    """One horizontal 100% bar; a segment is labelled inside when wide enough, else below."""
    ax = fig.add_axes(rect)
    ax.set_xlim(0, total)
    ax.set_ylim(0, 1)
    ax.axis("off")
    left, below = 0, []
    for lab, v, col in parts:
        if not v:
            continue
        ax.barh(0.5, v, left=left, height=1.0, color=col, edgecolor=BG, linewidth=1.5)
        txt = f"{lab} {v}"
        dark_bg = col in text_dark
        if v / total * rect[2] * 1920 >= len(txt) * 10 + 16:
            ax.text(left + v / 2, 0.5, txt, ha="center", va="center", fontsize=15,
                    color="#f8f9fb" if dark_bg else INK, fontweight="bold" if dark_bg else "normal")
        else:
            below.append(txt)
        left += v
    if below:  # small segments: one line of labels directly under the bar
        fig.text(rect[0], rect[1] - 0.028, "   ".join(below), fontsize=15, color=INK, va="baseline")


def draw_profile(fig, rect, uid: str, d: dict, r: dict) -> None:
    x0, y0, w, h = rect
    sp = d["span"].set_index("uniprot_id").loc[uid]
    s0, s1 = int(sp["analysis_span_start"]), int(sp["analysis_span_end"])
    pos = np.arange(s0, s1 + 1)
    prof = d["profile"][d["profile"]["uniprot_id"] == uid].set_index("canonical_uniprot_position").reindex(pos)
    n = prof["n_pdb_entries_observed_at_position"].fillna(0).to_numpy()
    seg = d["segmap"][d["segmap"]["uniprot_id"] == uid]
    sites = d["sites"][d["sites"]["uniprot_id"] == uid]

    ax_h = fig.add_axes([x0, y0 + h * 0.905, w, h * 0.075])
    ax_s = fig.add_axes([x0, y0 + h * 0.855, w, h * 0.04], sharex=ax_h)
    ax_m = fig.add_axes([x0, y0 + h * 0.22, w, h * 0.62], sharex=ax_h)
    ax_c = fig.add_axes([x0, y0, w, h * 0.17], sharex=ax_h)
    for a in (ax_h, ax_s):
        a.set_ylim(0, 1)
        a.axis("off")
    ax_h.set_xlim(s0 - 0.5, s1 + 0.5)
    for a in (ax_m, ax_c):
        a.spines["left"].set_color(GRID)
        a.spines["bottom"].set_color(GRID)

    # unobserved stretches: grey bands in the profile and coverage panels
    gaps = n == 0
    if gaps.any():
        edges = np.flatnonzero(np.diff(np.r_[0, gaps.astype(int), 0]))
        for a0, a1 in zip(edges[::2], edges[1::2]):
            for a in (ax_m, ax_c):
                a.axvspan(pos[a0] - 0.5, pos[a1 - 1] + 0.5, color="#eceef2", linewidth=0, zorder=0)
            if a1 - a0 >= 12 and n.max() > 0:
                ax_m.text((pos[a0] + pos[a1 - 1]) / 2, Y_LIM[0] + 0.35, "not observed", ha="center", fontsize=13,
                          color=MUTED)

    # helix track
    blocks = helix_blocks(seg)
    px_per_res = w * 1920 / (s1 - s0 + 1)
    for b in blocks:
        h12 = b["label"] == "H12"
        ax_h.add_patch(Rectangle((b["start"] - 0.5, 0.05), b["end"] - b["start"] + 1, 0.9,
                                 facecolor="#f6dfbf" if h12 else PALE, edgecolor=AF2_C if h12 else SLATE_LT,
                                 linewidth=0.8, hatch="////" if b["review"] else None))
        if (b["end"] - b["start"] + 1) * px_per_res >= len(b["label"]) * 10 + 6:
            ax_h.text((b["start"] + b["end"]) / 2, 0.5, b["label"], ha="center", va="center", fontsize=14,
                      color=INK, fontweight="bold" if h12 else "normal")
    if not blocks:
        ax_h.text(s0, 0.5, "no professor helix map", fontsize=14, color=MUTED, va="center")
    # site ticks with a minimal key at the right end
    lbp = sites[(sites["functional_site"] == "LIGAND_BINDING_POCKET")]["canonical_uniprot_position"]
    groove = sites[(sites["functional_site"] == "COACTIVATOR_GROOVE") &
                   (sites["definition_source"] == "ROLE_D_CURATED")]["canonical_uniprot_position"]
    ax_s.vlines(lbp, 0.0, 1.0, color=POCKET_C, linewidth=1.6)
    ax_s.vlines(groove, 0.0, 1.0, color=GROOVE_C, linewidth=2.6)
    if len(lbp):
        fig.text(x0 + w + 0.006, y0 + h * 0.875, "pocket", fontsize=13, color=POCKET_C, va="center")
    if len(groove):
        fig.text(x0 + w + 0.006, y0 + h * 0.840, "groove", fontsize=13, color=GROOVE_C, va="center")

    # median + IQR
    med = prof["median_z"].to_numpy(dtype=float)
    q25, q75 = prof["q25_z"].to_numpy(dtype=float), prof["q75_z"].to_numpy(dtype=float)
    sup = n >= MIN_ENTRIES
    faint = (n >= 1) & ~sup
    ax_m.axhline(0, color=MUTED, linewidth=0.8, linestyle=(0, (4, 3)))
    ax_m.fill_between(pos, np.where(sup, q25, np.nan), np.where(sup, q75, np.nan), color="#cfd8e5", linewidth=0)
    ax_m.plot(pos, np.where(sup, med, np.nan), color=HOLO_C, linewidth=2.0)
    ax_m.plot(pos, np.where(faint | np.r_[False, sup[:-1] & faint[1:]] | np.r_[faint[1:] & sup[:-1], False], med, np.nan),
              color=SLATE_LT, linewidth=1.3, linestyle=(0, (2, 2)))
    clipped = np.where((n >= 1) & (med > Y_LIM[1]))[0]
    ax_m.scatter(pos[clipped], np.full(len(clipped), Y_LIM[1] - 0.12), marker="^", s=22, color=SLATE, zorder=3)
    ax_m.set_ylim(*Y_LIM)
    ax_m.set_ylabel("B-factor z  (line: median · band: IQR)", fontsize=14, color=INK)
    fig.canvas.draw()
    yl = ax_m.yaxis.label.get_window_extent()
    fig.text((yl.x0 - 4) / fig.bbox.width, (yl.y0 + yl.y1) / 2 / fig.bbox.height, r["zscored_label"],
             rotation=90, ha="right", va="center", fontsize=12, color=MUTED)
    ax_m.grid(axis="y", color=GRID, linewidth=0.8)
    ax_m.tick_params(labelbottom=False, labelsize=13, length=0)
    for b in blocks:
        if b["label"] == "H12":
            ax_m.axvspan(b["start"] - 0.5, b["end"] + 0.5, color=AF2_C, alpha=0.08, linewidth=0)
    if n.max() == 0:
        ax_m.text((s0 + s1) / 2, 1.2, "no structures", ha="center", fontsize=22, color=MUTED)

    # coverage
    ax_c.bar(pos, n, width=1.0, color=SLATE_LT, linewidth=0)
    ax_c.axhline(MIN_ENTRIES, color=LINE5_C, linewidth=1.0, linestyle=(0, (4, 3)))
    ax_c.set_ylabel("structures", fontsize=14, color=INK)
    ax_c.set_ylim(0, max(n.max() * 1.1, MIN_ENTRIES + 1))
    ax_c.set_yticks([0, int(max(n.max(), MIN_ENTRIES))])
    ax_c.tick_params(labelsize=13, length=0)
    ax_c.set_xlabel(f"canonical position ({r['reference_accession']})", fontsize=14, color=INK)
    step = 25 if s1 - s0 < 320 else 50
    ax_c.set_xticks(np.arange((s0 // step + 1) * step, s1 + 1, step))


def wrap(text: str, width: int) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width) or [""]


def page(uid: str, d: dict, r: dict) -> plt.Figure:
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    L = 0.04
    n = r["n_entries_contributing"]

    # 1 header
    name = fig.text(L, 0.935, r["common_name"], fontsize=44, color=INK, fontweight="bold", va="baseline")
    fig.canvas.draw()
    nx = L + name.get_window_extent().width / 1920 + 0.015
    fig.text(nx, 0.935, f"{r['nr_code']}  ·  {n} usable PDB{'s' if n != 1 else ''}", fontsize=26, color=MUTED, va="baseline")
    fig.text(L, 0.885, r["interpretation"], fontsize=22, color=INK, va="baseline")

    # 2 data composition | 3 biological context (top band)
    TOP = 0.815
    block_label(fig, L, TOP, "Structures")
    if n:
        stacked_bar(fig, [L, TOP - 0.065, 0.27, 0.045],
                    [("HOLO", r["n_holo"], HOLO_C), ("APO", r["n_apo"], APO_C), ("AMBIG.", r["n_ambiguous"], AMB_C),
                     ("TRACE", r["n_trace_lbd"], TRACE_C)], n)
        sub = [f"{r['n_entries_discovered']} found", f"{r['n_entries_coordinate']} with coordinates"]
        if not pd.isna(r["median_observed_span_fraction"]):
            sub.append(f"{r['median_observed_span_fraction']:.0%} of LBD observed (median)")
        fig.text(L, TOP - 0.13, "  ·  ".join(sub), fontsize=14, color=MUTED, va="baseline")
        if r["n_apo"] == 0 and r["n_trace_lbd"] < n:
            fig.text(L + 0.275, TOP - 0.043, "no APO", fontsize=16, color=AF2_C, fontweight="bold", va="center")
    else:
        fig.text(L, TOP - 0.05, "No usable LBD structures", fontsize=26, color=AF2_C, fontweight="bold", va="center")
        fig.text(L, TOP - 0.10, f"{r['n_entries_discovered']} entr{'y' if r['n_entries_discovered'] == 1 else 'ies'} found"
                 + (f" ({r['no_coordinate_reasons']})" if r["no_coordinate_reasons"] else ""), fontsize=14, color=MUTED)

    LX = 0.385
    block_label(fig, LX, TOP, "Ligands")
    if r["n_holo_with_ligand_id"]:
        big = fig.text(LX, TOP - 0.06, f"{r['n_unique_ligands']}", fontsize=40, color=INK, fontweight="bold", va="baseline")
        fig.canvas.draw()
        bx = LX + big.get_window_extent().width / 1920 + 0.008
        fig.text(bx, TOP - 0.06, "unique" + (f" · {r['top_ligand']} in {r['top_ligand_share_of_holo']:.0%}"
                                              if r["ligand_context"] == "one dominant ligand" else ""),
                 fontsize=18, color=INK, va="baseline")
        top3 = [t_ for t_ in r["top5_ligands"].split("; ") if t_][:3]
        fig.text(LX, TOP - 0.105, "  ·  ".join(t_.replace(" (", " ").rstrip(")") for t_ in top3), fontsize=16,
                 color=MUTED, va="baseline")
    else:
        fig.text(LX, TOP - 0.06, "no bound ligands" if n else "–", fontsize=22, color=MUTED, va="baseline")
    if r["n_unique_ligand_partner_contexts"]:
        k = r["n_unique_ligand_partner_contexts"]
        fig.text(LX, TOP - 0.135, f"{k} ligand–partner context{'s' if k != 1 else ''}", fontsize=15, color=INK, va="baseline")
        fig.text(LX, TOP - 0.16, f"dominant: {r['dominant_ligand_partner_context_display']} "
                 f"({r['dominant_ligand_partner_context_fraction']:.0%})", fontsize=14, color=MUTED, va="baseline")

    PX = 0.60
    block_label(fig, PX, TOP, "Partners / cofactors")
    if n:
        if r["n_standalone"] == n:
            fig.text(PX, TOP - 0.06, "Standalone: 100%", fontsize=24, color=INK, va="baseline")
        else:
            other = n - r["n_standalone"] - r["n_partner_coregulator_peptide_primary"] - r["n_partner_nr_primary"]
            stacked_bar(fig, [PX, TOP - 0.065, 0.36, 0.045],
                        [("coregulator", r["n_partner_coregulator_peptide_primary"], HOLO_C),
                         ("NR partner", r["n_partner_nr_primary"], "#4f6d9a"),
                         ("other", other, APO_C), ("standalone", r["n_standalone"], AMB_C)], n,
                        text_dark=(HOLO_C, "#4f6d9a"))
            if r["top_partner"]:
                fig.text(PX, TOP - 0.13, f"most common: {r['top_partner']}", fontsize=14, color=MUTED, va="baseline")
            if r["n_unique_partner_identities"]:
                k = r["n_unique_partner_identities"]
                fig.text(PX, TOP - 0.16, f"{k} unique partner identit{'ies' if k != 1 else 'y'}", fontsize=15,
                         color=INK, va="baseline")
    else:
        fig.text(PX, TOP - 0.06, "–", fontsize=22, color=MUTED, va="baseline")

    # 4 flexibility (main visual)
    draw_profile(fig, (L + 0.03, 0.065, 0.655, 0.555), uid, d, r)

    # heterogeneity, constructs and caveats beside the plot
    RX = 0.79
    block_label(fig, RX, 0.585, "Between-PDB heterogeneity")
    if r["heterogeneity_class"] != "not assessable":
        fig.text(RX, 0.52, r["heterogeneity_class"].upper(), fontsize=36, color=INK, fontweight="bold", va="baseline")
        for i, (k, v) in enumerate([("median IQR", f"{r['median_iqr_z']:.2f}"), ("90th percentile", f"{r['p90_iqr_z']:.2f}"),
                                    ("supported positions", f"{r['n_supported_positions']} / {r['n_span_positions']}")]):
            fig.text(RX, 0.475 - i * 0.034, k, fontsize=16, color=MUTED, va="baseline")
            fig.text(0.965, 0.475 - i * 0.034, v, fontsize=16, color=INK, va="baseline", ha="right")
    else:
        fig.text(RX, 0.52, "not assessable", fontsize=30, color=AF2_C, fontweight="bold", va="baseline")
        why = (f"Only {n} usable structure{'s' if n != 1 else ''}" if n < MIN_ENTRIES else
               f"Only {r['n_supported_positions']} positions in ≥5 structures")
        fig.text(RX, 0.475, why if n else "No structures", fontsize=18, color=INK, va="baseline")

    y = 0.335 if r["heterogeneity_class"] != "not assessable" else 0.415
    if n:
        block_label(fig, RX, y, "Constructs & crystals")
        y -= 0.045
        for line in construct_lines(r):
            fig.text(RX, y, line, fontsize=16, color=INK, va="baseline")
            y -= 0.034
        y -= 0.025
    cav = notes(r)
    if cav:
        block_label(fig, RX, y, "Caveats")
        y -= 0.04
        for c in cav:
            if y - 0.03 * len(wrap(c, 33)) < 0.02:
                break  # never draw past the page; every caveat stays in the summary table
            for j, line in enumerate(wrap(c, 33)):
                fig.text(RX + (0 if j == 0 else 0.012), y, ("•  " if j == 0 else "") + line, fontsize=15, color=INK, va="baseline")
                y -= 0.03
            y -= 0.006
    return fig


def construct_lines(r: dict) -> list[str]:
    """Compact construct/experiment summary: always the patterns, resolution and crystal forms;
    other features only when present."""
    n = r["n_entries_contributing"]
    lines = [f"{r['n_construct_patterns']} construct pattern{'s' if r['n_construct_patterns'] != 1 else ''}"
             + (f" · {r['n_crystal_forms']} crystal form{'s' if r['n_crystal_forms'] != 1 else ''}" if n > 1 else "")]
    if not pd.isna(r["resolution_median"]):
        lines.append(f"{r['resolution_median']:.2f} Å median" + (f" ({r['resolution_q25']:.2f}–{r['resolution_q75']:.2f})"
                                                                 if n >= 4 else ""))
    # unusual features only (thresholds in RULES["construct_line"])
    if r["n_pdb_engineered_mutant"] and r["engineered_mutant_fraction"] >= 0.2:  # audited, depositor-confirmed
        top = r["top_engineered_substitutions"].split(" (")[0] if r["top_engineered_substitutions"] else ""
        lines.append(f"engineered mutants: {r['engineered_mutant_fraction']:.0%}" + (f" (top {top})" if top else ""))
    if r["substitution_uncertain_fraction"] >= 0.1 and r["n_pdb_substitution_uncertain"] >= 2:
        lines.append(f"substitution status uncertain: {r['substitution_uncertain_fraction']:.0%}")
    extra = []
    if r["n_fusion_or_chimera"] and (r["n_fusion_or_chimera"] / n >= 0.05 or n < 10):
        extra.append(f"{r['n_fusion_or_chimera']} fusion/chimera")
    if r["n_internal_deletion_constructs"] and (r["n_internal_deletion_constructs"] / n >= 0.1 or n < 10):
        extra.append(f"{r['n_internal_deletion_constructs']} internal deletion")
    if extra:
        lines.append(" · ".join(extra))
    return lines


def overview_figures(S: pd.DataFrame) -> list[tuple[str, plt.Figure]]:
    figs = []
    S = S.sort_values("nr_code").reset_index(drop=True)
    for part, sub in enumerate([S.iloc[:24], S.iloc[24:]], 1):
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.text(0.035, 0.95, f"Dataset at a glance ({part}/2)", fontsize=30, color=INK, fontweight="bold", va="center")
        fig.text(0.035, 0.905, "One row per receptor. First panel: entries with z-scores (dark) of all PDB entries found (light). "
                 "Other panels use the z-scored entries.",
                 fontsize=16, color=MUTED, va="center")
        y = np.arange(len(sub))[::-1]
        names = [f"{a}  {b}" for a, b in zip(sub["common_name"], sub["nr_code"])]
        cols = [(0.13, 0.17, "Z-scored / total entries"), (0.33, 0.17, "Apo / holo / ambiguous / TRACE (%)"),
                (0.53, 0.14, "With a partner chain (%)"), (0.70, 0.12, "Unique ligands"),
                (0.85, 0.12, "Median IQR(z)")]
        axes = [fig.add_axes([x, 0.06, w, 0.78]) for x, w, _ in cols]
        for ax, (_, _, t) in zip(axes, cols):
            ax.set_title(t, fontsize=14, color=INK, loc="left")
            ax.set_ylim(-0.6, len(sub) - 0.4)
            ax.set_yticks([])
            ax.grid(axis="x", color=GRID)
            ax.tick_params(labelsize=11)
        axes[0].set_yticks(y, names, fontsize=13, color=INK)
        axes[0].tick_params(axis="y", length=0)
        axes[0].barh(y, sub["n_entries_total"], color=PALE, edgecolor=SLATE_LT, linewidth=0.6, height=0.6)
        axes[0].barh(y, sub["n_entries_with_zscores"], color=SLATE, height=0.6)
        for yi, v, t in zip(y, sub["n_entries_with_zscores"], sub["n_entries_total"]):
            axes[0].text(t, yi, f" {v} / {t}", va="center", fontsize=11, color=INK)
        axes[0].set_xlim(0, S["n_entries_total"].max() * 1.32)
        left = np.zeros(len(sub))
        for k, c in [("n_holo", HOLO_C), ("n_apo", APO_C), ("n_ambiguous", AMB_C), ("n_trace_lbd", TRACE_C)]:
            frac = np.where(sub["n_entries_contributing"] > 0, sub[k] / sub["n_entries_contributing"].replace(0, 1), 0) * 100
            axes[1].barh(y, frac, left=left, color=c, height=0.6)
            left += frac
        axes[1].set_xlim(0, 100)
        axes[2].barh(y, sub["frac_with_partner"].fillna(0) * 100, color=SLATE, height=0.6)
        axes[2].set_xlim(0, 100)
        axes[3].barh(y, sub["n_unique_ligands"], color=HOLO_C, height=0.6)
        axes[3].set_xlim(0, S["n_unique_ligands"].max() * 1.1)
        axes[4].barh(y, sub["median_iqr_z"].fillna(0), color=SLATE_LT, height=0.6)
        for yi, v, cl, nn in zip(y, sub["median_iqr_z"], sub["heterogeneity_class"], sub["n_entries_contributing"]):
            if nn == 0:
                continue
            axes[4].text((0 if pd.isna(v) else v), yi, " –" if pd.isna(v) else f" {v:.2f} {cl[0].upper() if cl != 'not assessable' else ''}",
                         va="center", fontsize=11, color=INK)
        axes[4].set_xlim(0, S["median_iqr_z"].max() * 1.45)
        for ax in axes[1:]:
            for yi, v in zip(y, sub["n_entries_contributing"]):
                if v == 0:
                    ax.text(ax.get_xlim()[1] * 0.02, yi, "no data", va="center", fontsize=11, color=MUTED)
        fig.text(0.035, 0.02, "Apo/holo colours: HOLO dark, APO mid, AMBIGUOUS light, TRACE_LBD orange.  "
                 "Heterogeneity letter: L/M/H = low/moderate/high tertile.", fontsize=12, color=MUTED)
        figs.append((f"overview_glance_{part}", fig))
    # heterogeneity overview
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.text(0.035, 0.95, "Where within-receptor heterogeneity is concentrated", fontsize=30, color=INK, fontweight="bold", va="center")
    fig.text(0.035, 0.905, "Each point = one receptor (positions seen in ≥5 entries). x = typical position IQR(z); "
             "y = upper tail (90th percentile); size = entries with z-scores.", fontsize=16, color=MUTED, va="center")
    ax = fig.add_axes([0.07, 0.09, 0.62, 0.76])
    T = S.dropna(subset=["median_iqr_z"])
    ax.scatter(T["median_iqr_z"], T["p90_iqr_z"], s=30 + T["n_entries_contributing"].clip(upper=400) * 1.2,
               color=SLATE, alpha=0.55, edgecolor=BG, linewidth=1)
    for _, rr in T.iterrows():
        ax.annotate(rr["common_name"], (rr["median_iqr_z"], rr["p90_iqr_z"]), xytext=(5, 4), textcoords="offset points",
                    fontsize=12, color=INK)
    ax.set_xlabel("median position IQR(z)", fontsize=14)
    ax.set_ylabel("90th-percentile position IQR(z)", fontsize=14)
    ax.grid(color=GRID)
    tx = fig.add_axes([0.73, 0.09, 0.24, 0.76])
    tx.axis("off")
    tx.text(0, 1, "Top-10% IQR positions, by element", fontsize=16, color=INK, fontweight="bold", va="top")
    agg = Counter()
    for v in T["top_decile_iqr_location"]:
        for part in v.split("; "):
            if part:
                k, c = part.rsplit(" ", 1)
                agg[k] += int(c)
    tot = sum(agg.values())
    for i, (k, c) in enumerate(agg.most_common(14)):
        tx.text(0, 0.93 - i * 0.058, f"{k}", fontsize=14, color=INK, va="top")
        tx.text(0.62, 0.93 - i * 0.058, f"{c}  ({c / tot:.0%})", fontsize=14, color=MUTED, va="top")
    tx.text(0, 0.06, "Counted over all receptors; 'between helices' = positions outside\nevery professor interval. "
            "Descriptive only (no burial/SASA control).", fontsize=12, color=MUTED, va="top")
    figs.append(("overview_heterogeneity", fig))
    return figs


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    d = load()
    uids = list(d["meta"].sort_values("nr_code")["uniprot_id"])
    first = pd.DataFrame([summarise(u, d, None) for u in uids])
    ok = first[first["n_supported_positions"] >= MIN_POS_FOR_CLASS]["median_iqr_z"]
    cuts = (round(ok.quantile(1 / 3), 3), round(ok.quantile(2 / 3), 3))
    S = pd.DataFrame([summarise(u, d, cuts) for u in uids])
    S["heterogeneity_cut_low_moderate"], S["heterogeneity_cut_moderate_high"] = cuts
    S.to_csv(TAB / "receptor_catalogue_summary_v2.csv", index=False, lineterminator="\n")
    # every context per receptor with its entry count (TRACE_LBD-only entries have no LBD context)
    E = d["entries"][d["entries"]["ligand_partner_context"] != ""]
    ctx = pd.concat([E.groupby(["uniprot_id", "nr_code", "common_name", c]).size().rename("n_entries").reset_index()
                     .rename(columns={c: "context"}).assign(context_type=k)
                     for c, k in [("ligand_partner_context", "ligand x partner identity"),
                                  ("state_partner_context", "ligand state x partner category")]])
    ctx["fraction_of_receptor_entries"] = (ctx["n_entries"] / ctx["uniprot_id"].map(E.groupby("uniprot_id").size())).round(4)
    ctx.sort_values(["nr_code", "context_type", "n_entries", "context"], ascending=[True, True, False, True]) \
        .to_csv(TAB / "catalogue_ligand_partner_contexts_v2.csv", index=False, lineterminator="\n")
    rules = dict(RULES, heterogeneity_cutoffs_median_iqr={"low_le": cuts[0], "high_gt": cuts[1],
                                                         "n_receptors_classed": int(len(ok))})
    (TAB / "catalogue_rules_v2.json").write_text(json.dumps(rules, indent=2, ensure_ascii=False) + "\n")
    with PdfPages(OUT / "receptor_catalogue_v2.pdf") as pdf:
        for name, fig in overview_figures(S):  # the PDF opens with the dataset-at-a-glance pages
            fig.savefig(FIG / f"{name}.png", dpi=100)
            if name.startswith("overview_glance"):
                pdf.savefig(fig)
            plt.close(fig)
        for _, r in S.iterrows():
            fig = page(r["uniprot_id"], d, r.to_dict())
            fig.savefig(FIG / f"receptor_{r['nr_code']}.png", dpi=100)
            pdf.savefig(fig)
            plt.close(fig)
    print(f"{len(S)} receptors; heterogeneity cut-offs {cuts}; pages in {FIG.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
