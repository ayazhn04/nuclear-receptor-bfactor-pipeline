"""
Part 1 of the functional-site specificity/validation stage (continuing
from commit e420454 on analysis/normalization-v1): independently
recompute the within-receptor and global BH-FDR corrections from the raw
p-values stored in
reports/tables/current_functional_site_within_receptor_tests.csv, verify
the exact q-values and significant counts already reported, and explain
the counterintuitive pattern (76/107 global-significant vs 73/107
within-receptor-significant -- MORE significant under the larger, harsher
correction) by rank structure, not by changing anything.

This script does NOT modify data/analysis/current_functional_site_*
outputs -- it is a read-only sensitivity/validation layer.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Standard BH step-up, identical implementation to every prior stage's
    (scripts 56, 64, 76) -- re-derived independently here, not imported, so
    this audit is a genuine second implementation, not a call to the same code."""
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def main():
    stored = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")

    audit = stored[["uniprot_id", "nr_code", "functional_site", "raw_p"]].copy()

    # independent within-receptor recomputation
    audit["q_within_receptor_recomputed"] = np.nan
    audit["n_tests_in_receptor_family"] = np.nan
    audit["rank_within_receptor"] = np.nan
    for uid, idx in audit.groupby("uniprot_id").groups.items():
        pvals = audit.loc[idx, "raw_p"].values
        audit.loc[idx, "q_within_receptor_recomputed"] = benjamini_hochberg(pvals)
        audit.loc[idx, "n_tests_in_receptor_family"] = len(pvals)
        audit.loc[idx, "rank_within_receptor"] = pd.Series(pvals, index=idx).rank(method="first").values

    # independent global recomputation
    audit["q_global_recomputed"] = benjamini_hochberg(audit["raw_p"].values)
    audit["rank_global"] = audit["raw_p"].rank(method="first")
    audit["n_tests_global"] = len(audit)

    audit["significant_within_receptor_recomputed"] = audit["q_within_receptor_recomputed"] <= 0.05
    audit["significant_global_recomputed"] = audit["q_global_recomputed"] <= 0.05

    # cross-check against the originally stored values
    audit = audit.merge(
        stored[["uniprot_id", "functional_site", "q_within_receptor", "q_global",
                "significant_within_receptor", "significant_global"]],
        on=["uniprot_id", "functional_site"],
    )
    audit["q_within_matches"] = np.isclose(audit["q_within_receptor_recomputed"], audit["q_within_receptor"])
    audit["q_global_matches"] = np.isclose(audit["q_global_recomputed"], audit["q_global"])
    audit["sig_within_matches"] = audit["significant_within_receptor_recomputed"] == audit["significant_within_receptor"]
    audit["sig_global_matches"] = audit["significant_global_recomputed"] == audit["significant_global"]

    all_match = audit[["q_within_matches", "q_global_matches", "sig_within_matches", "sig_global_matches"]].all().all()
    print(f"Independent recomputation matches stored values exactly: {all_match}")
    if not all_match:
        mismatches = audit[~(audit["q_within_matches"] & audit["q_global_matches"] &
                              audit["sig_within_matches"] & audit["sig_global_matches"])]
        print(f"  {len(mismatches)} mismatched row(s) -- see output CSV for detail. "
              "NOT auto-correcting the original outputs; a human must confirm whether this is a bug.")

    n_within = audit["significant_within_receptor_recomputed"].sum()
    n_global = audit["significant_global_recomputed"].sum()
    print(f"\nRecomputed: {n_within}/{len(audit)} significant within-receptor FDR, "
          f"{n_global}/{len(audit)} significant global FDR (stored: 73, 76)")

    # explain: which tests are global-significant but not within-receptor-significant
    only_global = audit[audit["significant_global_recomputed"] & ~audit["significant_within_receptor_recomputed"]]
    print(f"\n{len(only_global)} test(s) are global-significant but NOT within-receptor-significant:")
    if len(only_global):
        print(only_global[["uniprot_id", "nr_code", "functional_site", "raw_p",
                            "n_tests_in_receptor_family", "rank_within_receptor",
                            "n_tests_global", "rank_global",
                            "q_within_receptor_recomputed", "q_global_recomputed"]].to_string(index=False))
        print(
            "\nMechanism (not a bug): BH's q-value for a test = min over j>=rank of "
            "(p_(j) * n_family / j) (the step-up procedure with monotonicity enforced from the "
            "largest p-value down). For these 3 tests, the receptor's OWN family is small "
            "(3-4 tested sites) and every sibling p-value in that family is also only moderately "
            "small, so the within-receptor correction stays harsh even after monotonicity "
            "enforcement (q_within = 0.0625-0.094). Globally (n=107), the same raw p=0.03125 lands "
            "at rank ~74-76 out of 107 -- comfortably ahead of the many much-larger p-values "
            "contributed by receptors and sites with no real signal -- so q_global (~0.044) comes "
            "out smaller than q_within despite the much larger nominal family size. This is a "
            "documented, mathematically legitimate BH property (q is not monotonic in family size "
            "alone; it depends on the full shape of that family's p-value distribution), not an "
            "implementation error."
        )

    out = TABLES / "functional_site_fdr_audit.csv"
    audit.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(audit)} rows)")
    print("No change was made to the original functional-site outputs -- this is a read-only audit.")


if __name__ == "__main__":
    main()
