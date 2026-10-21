"""
Post-Phase-2 validation of the corrected v2 rerun, and the old (v1 frozen) vs v2
result comparison.

Checks (each PASS/FAIL, written to reports/rebuild_v2/phase2_validation_v2.json):
  * no v3 canonical label / v3 assignment vocabulary in any v2 Phase 2 output;
  * no hard-coded v1 conclusion text in the v2 summaries;
  * missing != 0: undefined professor elements stay NaN in the element tables; every
    imputed matrix cell is NaN before imputation and equals the cohort median after;
  * dual-label positions: exactly two memberships at PXR 261 / CAR 344, one residue-test
    row per position, no joined label in any element column;
  * derived cohorts and test counts agree across the tables that report them;
  * script 95 consistency audit passed with zero failures.

Comparison (reports/rebuild_v2/old_vs_v2_comparison.csv): every metric of the frozen
final_project_numbers.csv next to the v2 value, plus headline result lists.

Reads only: frozen tracked v1 outputs (read-only) and data/rebuild_v2, reports/rebuild_v2.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
V1A, V1T = PROJECT_ROOT / "data" / "analysis", PROJECT_ROOT / "reports" / "tables"
V2A, V2T = PROJECT_ROOT / "data" / "rebuild_v2" / "analysis", PROJECT_ROOT / "reports" / "rebuild_v2" / "tables"
OUT = PROJECT_ROOT / "reports" / "rebuild_v2"
V3_TOKENS = ["loop_", "beta_sheet", "LOOP_UNCERTAIN", "own_canonical_template", "structural_alignment_transfer",
             "complementary_structure_rescue", "NO_EVIDENCE"]
# v1-specific wording (a bare number such as "29 receptors" can be a legitimate v2 count)
V1_PHRASES = ["29-receptor", "~29", "29 primary-cohort", "8 segments", "8-feature", "loop_H3_H5", "35 segment",
              "~35 tests", "1-vs-28", "**0 segments**", "no robust evidence", "3 within-receptor", "0 global",
              "loop_H7_H10", "loop_H6_H10", "All 3 significant", "best 0.41"]
DUAL = {("O75469", 261): ("H3", "H4"), ("Q14994", 344): ("H11'", "H12")}

checks: list[dict] = []


def check(name: str, ok: bool, detail="") -> None:
    checks.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})


def phase2_files() -> list[Path]:
    m = json.loads((OUT / "phase2_export_manifest.json").read_text())
    return [PROJECT_ROOT / f for f in m["exported_files"]]


def scan_v3(files) -> None:
    hits = []
    for f in files:
        if f.suffix == ".csv":
            df = pd.read_csv(f, keep_default_na=False, dtype=str)
        elif f.suffix == ".parquet":
            df = pd.read_parquet(f)
        else:
            continue
        cols = [c for c in df.columns if df[c].dtype == object or str(df[c].dtype).startswith("str")]
        for c in cols:
            vals = df[c].astype(str)
            for t in V3_TOKENS:
                if vals.str.contains(t, regex=False).any():
                    hits.append(f"{f.name}:{c}:{t}")
        for t in V3_TOKENS:
            if any(t in str(c) for c in df.columns):
                hits.append(f"{f.name}:column:{t}")
    check("no v3 labels or v3 assignment vocabulary in Phase 2 outputs", not hits, hits[:20])


def scan_v1_phrases() -> None:
    hits = []
    for md in V2T.glob("*.md"):
        text = md.read_text()
        hits += [f"{md.name}: {p!r}" for p in V1_PHRASES if p in text]
    check("no hard-coded v1 conclusion phrases in v2 summaries", not hits, hits)


def missing_not_zero() -> None:
    seg = pd.read_csv(V2A / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    defined = seg[seg["structural_segment"] != ""].groupby("uniprot_id")["structural_segment"].agg(set).to_dict()
    med = pd.read_parquet(V2A / "receptor_segment_feature_matrix_median.parquet").set_index("uniprot_id")
    bad = [(u, c) for u in med.index for c in med.columns
           if c not in defined.get(u, set()) and not pd.isna(med.loc[u, c])]
    check("undefined professor elements are NaN (never a value) in the receptor x element matrix", not bad, bad[:10])
    prof = pd.read_parquet(V2A / "receptor_segment_flexibility_profiles.parquet")
    bad2 = prof[prof.apply(lambda r: r["structural_segment"] not in defined.get(r["uniprot_id"], set()), axis=1)]
    check("element profiles exist only for professor-defined elements", bad2.empty, len(bad2))
    mask = pd.read_parquet(V2A / "receptor_segment_matrix_primary_imputation_mask.parquet").set_index("uniprot_id")
    unimp = pd.read_parquet(V2A / "receptor_segment_matrix_primary_unscaled.parquet").set_index("uniprot_id") \
        .reindex(mask.index)[mask.columns]
    raw = med.reindex(mask.index)[mask.columns]
    m = mask.values.astype(bool)
    ok_nan = bool(np.isnan(raw.values[m].astype(float)).all())
    obs_med = raw.median()
    ok_med = all(np.isclose(unimp.loc[u, c], obs_med[c]) for u in mask.index for c in mask.columns if mask.loc[u, c])
    ok_obs = bool(np.allclose(unimp.values[~m].astype(float), raw.values[~m].astype(float)))
    zero_imputed = [(u, c) for u in mask.index for c in mask.columns if mask.loc[u, c] and unimp.loc[u, c] == 0]
    check("every imputed cell was missing before imputation", ok_nan, int(mask.values.sum()))
    check("imputed cells equal the cohort median of observed values (never 0 by construction)",
          ok_med and not zero_imputed, zero_imputed)
    check("observed cells unchanged by imputation", ok_obs)
    imp = pd.read_csv(V2T / "modeling_imputed_cells.csv")
    check("imputed-cell documentation lists exactly the imputed cells", len(imp) == int(mask.values.sum()),
          f"{len(imp)} documented / {int(mask.values.sum())} imputed")


def dual_labels() -> None:
    seg = pd.read_csv(V2A / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    c = seg.groupby(["uniprot_id", "canonical_uniprot_position"]).size()
    check("only PXR 261 and CAR 344 carry two professor elements", {k: int(v) for k, v in c[c > 1].items()} == {k: 2 for k in DUAL})
    res = pd.read_parquet(V2A / "apo_holo_residue_effects.parquet")
    check("apo/holo residue tests: one row per position", not res.duplicated(["uniprot_id", "canonical_uniprot_position"]).any())
    joined = res[res["structural_segment"].str.contains("|", regex=False)]
    check("joined labels appear only at dual positions (display only)",
          set(zip(joined["uniprot_id"], joined["canonical_uniprot_position"])) <= set(DUAL), joined[["nr_code", "canonical_uniprot_position"]].values.tolist())
    for f in ["pdb_segment_flexibility_features.parquet", "receptor_segment_flexibility_profiles.parquet",
              "superfamily_segment_flexibility_atlas.csv"]:
        d = pd.read_parquet(V2A / f) if f.endswith("parquet") else pd.read_csv(V2A / f)
        check(f"no joined label used as an element in {f}", not d["structural_segment"].astype(str).str.contains("|", regex=False).any())
    defs = pd.read_parquet(V2A / "current_functional_site_definitions.parquet")
    dup = defs[defs.duplicated(["uniprot_id", "functional_site", "canonical_uniprot_position"], keep=False)]
    check("site-definition duplicates are only the intended dual memberships",
          set(zip(dup["uniprot_id"], dup["canonical_uniprot_position"])) <= set(DUAL), dup[["nr_code", "functional_site", "canonical_uniprot_position", "structural_segment"]].values.tolist())
    tiers = pd.read_csv(V2T / "functional_site_matched_control_tiers.csv", keep_default_na=False)
    inter_leak = []
    for f in phase2_files():
        if f.suffix == ".csv" and f.name != "functional_site_matched_control_tiers.csv":
            if "INTER:" in f.read_text():
                inter_leak.append(f.name)
    check("derived inter-element strata appear only in the matched-control tier table", not inter_leak, inter_leak)
    check("matched-control tier table names strata as matching strata", "site_matching_strata" in tiers.columns)


def consistency() -> None:
    cs = json.loads((V2A / "modeling_cohort_matrix_summary.json").read_text())
    coh = pd.read_csv(V2A / "modeling_receptor_cohort.csv")
    sup = pd.read_csv(V2T / "modeling_feature_support.csv")
    feats = sorted(sup.loc[sup["fraction_supported"] >= cs["primary_support_threshold"], "structural_segment"])
    check("primary features == elements with >= 50% support", feats == sorted(cs["primary_features"]), feats)
    check("cohort size consistent (json, cohort table, matrix)",
          cs["n_cohort_included"] == int(coh["included_primary"].sum()) ==
          len(pd.read_parquet(V2A / "receptor_segment_matrix_primary_scaled.parquet")))
    fp = pd.read_csv(V2T / "final_project_numbers.csv").set_index("metric")["value"]
    kw = pd.read_csv(V2T / "subfamily_segment_kruskal.csv")
    seg = pd.read_csv(V2T / "apo_holo_segment_tests.csv")
    within = pd.read_csv(V2T / "current_functional_site_within_receptor_tests.csv")
    pairs = {
        "primary_clustering_receptor_cohort": int(coh["included_primary"].sum()),
        "modeling_primary_features": len(cs["primary_features"]),
        "total_omnibus_tests": int((kw["test_status"] == "TESTED").sum()),
        "fdr_significant": int(kw.loc[kw["test_status"] == "TESTED", "significant_fdr_05"].sum()),
        "n_tests": len(seg), "within_receptor_fdr_significant": int(seg["significant_within_receptor"].sum()),
        "original_within_receptor_tests": len(within),
    }
    bad = {k: (fp.get(k), v) for k, v in pairs.items() if str(fp.get(k)) != str(v)}
    check("final_project_numbers agree with the result tables", not bad, bad)
    audit = (V2T / "final_statistical_consistency_audit.md").read_text()
    m = re.search(r"\*\*(\d+)/(\d+) PASS, (\d+) FAIL\.\*\*", audit)
    check("script 95 consistency audit: all claims PASS", bool(m) and m.group(3) == "0",
          m.group(0) if m else "pattern not found")


def comparison() -> pd.DataFrame:
    v1 = pd.read_csv(V1T / "final_project_numbers.csv").set_index("metric")
    v2 = pd.read_csv(V2T / "final_project_numbers.csv").set_index("metric")
    rows = []
    for metric in list(dict.fromkeys(list(v1.index) + list(v2.index))):
        cat = (v2 if metric in v2.index else v1).loc[metric, "category"]
        a = v1.loc[metric, "value"] if metric in v1.index else ""
        b = v2.loc[metric, "value"] if metric in v2.index else ""
        rows.append({"category": cat, "metric": metric, "v1_frozen": a, "v2_corrected": b,
                     "same": str(a) == str(b)})

    def sig_list(tbl, filt, cols, root):
        d = pd.read_csv(root / tbl)
        return "; ".join(" ".join(str(x) for x in r) for r in d[d[filt]][cols].values)
    for name, tbl, filt, cols in [
        ("apo/holo within-FDR significant (receptor element)", "apo_holo_segment_tests.csv", "significant_within_receptor", ["nr_code", "structural_segment"]),
        ("subfamily FDR-significant (scheme element)", "subfamily_segment_kruskal.csv", "significant_fdr_05", ["grouping_scheme", "structural_segment"]),
        ("functional-site superfamily FDR-significant", "current_functional_site_superfamily_tests.csv", "significant_fdr_05", ["functional_site"]),
    ]:
        rows.append({"category": "headline_lists", "metric": name, "v1_frozen": sig_list(tbl, filt, cols, V1T),
                     "v2_corrected": sig_list(tbl, filt, cols, V2T), "same": None})
    c1 = pd.read_csv(V1T / "functional_site_interpretation_classification.csv").set_index("functional_site")["interpretation_category"]
    c2 = pd.read_csv(V2T / "functional_site_interpretation_classification.csv").set_index("functional_site")["interpretation_category"]
    for s in sorted(set(c1.index) | set(c2.index)):
        rows.append({"category": "headline_lists", "metric": f"specificity classification {s}",
                     "v1_frozen": c1.get(s, ""), "v2_corrected": c2.get(s, ""), "same": c1.get(s) == c2.get(s)})
    sup1 = pd.read_csv(V1T / "current_functional_site_superfamily_tests.csv").set_index("functional_site")["median_delta_vs_background"]
    sup2 = pd.read_csv(V2T / "current_functional_site_superfamily_tests.csv").set_index("functional_site")["median_delta_vs_background"]
    for s in sorted(set(sup1.index) | set(sup2.index)):
        rows.append({"category": "headline_lists", "metric": f"superfamily median site-minus-background delta {s}",
                     "v1_frozen": round(float(sup1.get(s, np.nan)), 3), "v2_corrected": round(float(sup2.get(s, np.nan)), 3),
                     "same": None})
    e1 = pd.read_csv(V1T / "apo_holo_receptor_eligibility.csv")
    e2 = pd.read_csv(V2T / "apo_holo_receptor_eligibility.csv")
    for tier in ["eligible_primary_5_5", "eligible_exploratory_3_3"]:
        a, b = ";".join(sorted(e1.loc[e1[tier], "nr_code"])), ";".join(sorted(e2.loc[e2[tier], "nr_code"]))
        rows.append({"category": "headline_lists", "metric": f"apo/holo {tier} receptors", "v1_frozen": a, "v2_corrected": b, "same": a == b})
    f1 = json.loads((V1A / "modeling_cohort_matrix_summary.json").read_text())["primary_features"]
    f2 = json.loads((V2A / "modeling_cohort_matrix_summary.json").read_text())["primary_features"]
    rows.append({"category": "headline_lists", "metric": "primary modeling features", "v1_frozen": ";".join(f1),
                 "v2_corrected": ";".join(f2), "same": f1 == f2})
    return pd.DataFrame(rows)


def main() -> None:
    files = phase2_files()
    scan_v3(files)
    scan_v1_phrases()
    missing_not_zero()
    dual_labels()
    consistency()
    comp = comparison()
    comp.to_csv(OUT / "old_vs_v2_comparison.csv", index=False)
    res = {"n_checks": len(checks), "n_fail": sum(c["status"] == "FAIL" for c in checks), "checks": checks}
    (OUT / "phase2_validation_v2.json").write_text(json.dumps(res, indent=2, default=str) + "\n")
    for c in checks:
        print(f"{c['status']}  {c['check']}" + (f"  -> {c['detail']}" if c["status"] == "FAIL" else ""))
    print(f"{res['n_checks'] - res['n_fail']}/{res['n_checks']} PASS")


if __name__ == "__main__":
    main()
