"""
Part 8 of the functional-site specificity/validation stage: classify the
evidence for each site class. This is a REPORTING classification driven
directly by the already-computed numbers (Parts 2, 6, 7), not an
arbitrary score:

SITE-SPECIFIC SUPPORTED:
    the confirmatory superfamily matched-background reanalysis (Part 7)
    is FDR-significant, AND at least 70% of individual receptor-level
    results that had both an original and a matched test preserve
    significance (Part 6) -- i.e. the effect is broad, not a few outliers.

STRUCTURAL-CONTEXT EXPLAINED:
    the original whole-LBD-background superfamily result was significant
    (prior stage) but the matched-background superfamily reanalysis
    (Part 7) is NOT significant -- the effect is attributable to
    structural context (e.g. helix vs. loop) rather than the site itself.

INSUFFICIENT RESIDUE-LEVEL MAPPING:
    Part 7 status is INSUFFICIENT_N (too few HIGH/MODERATE-confidence
    RESIDUE-level definitions to test at all).

MIXED:
    the confirmatory superfamily result IS significant, but fewer than
    70% of individual receptor-level results preserve significance under
    matching -- broad average signal with substantial per-receptor
    inconsistency.
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TABLES = PROJECT_ROOT / "reports" / "tables"

MIXED_THRESHOLD = 0.70  # fixed, documented threshold for "broad vs. a few outliers"


def main():
    superfamily_original = pd.read_csv(TABLES / "current_functional_site_superfamily_tests.csv")
    superfamily_matched = pd.read_csv(TABLES / "functional_site_superfamily_matched_reanalysis.csv")
    sensitivity = pd.read_csv(TABLES / "functional_site_background_sensitivity.csv")

    rows = []
    all_sites = sorted(set(superfamily_original["functional_site"]) | set(superfamily_matched["functional_site"]))
    for site in all_sites:
        orig_row = superfamily_original[superfamily_original["functional_site"] == site]
        matched_row = superfamily_matched[superfamily_matched["functional_site"] == site]
        orig_sig = bool(orig_row.iloc[0]["significant_fdr_05"]) if len(orig_row) else False
        matched_status = matched_row.iloc[0]["status"] if len(matched_row) else "NOT_TESTED"
        matched_sig = bool(matched_row.iloc[0]["significant_fdr_05"]) if len(matched_row) and matched_status == "TESTED" else False

        site_sensitivity = sensitivity[sensitivity["functional_site"] == site]
        n_with_both = len(site_sensitivity)
        n_preserved = int(site_sensitivity["significance_preserved"].sum()) if n_with_both else 0
        n_orig_sig_individual = int(site_sensitivity["significant_original"].sum()) if n_with_both else 0
        preservation_fraction = (n_preserved / n_orig_sig_individual) if n_orig_sig_individual else None

        if matched_status == "INSUFFICIENT_N":
            category = "INSUFFICIENT_RESIDUE_LEVEL_MAPPING"
            rationale = "Too few HIGH/MODERATE-confidence RESIDUE-level definitions to test the matched-background superfamily reanalysis at all."
        elif matched_sig and (preservation_fraction is None or preservation_fraction >= MIXED_THRESHOLD):
            category = "SITE_SPECIFIC_SUPPORTED"
            rationale = (
                f"Matched-background superfamily reanalysis remains FDR-significant "
                f"({'no individually-significant original results to check preservation against' if preservation_fraction is None else f'{preservation_fraction:.0%} of individually-significant original results preserve significance'})."
            )
        elif matched_sig and preservation_fraction is not None and preservation_fraction < MIXED_THRESHOLD:
            category = "MIXED"
            rationale = (
                f"Matched-background superfamily reanalysis is FDR-significant overall, but only "
                f"{preservation_fraction:.0%} of individually-significant original receptor-level "
                f"results preserve significance under matching -- broad average signal with "
                "substantial per-receptor inconsistency."
            )
        elif orig_sig and not matched_sig:
            category = "STRUCTURAL_CONTEXT_EXPLAINED"
            rationale = "Original whole-LBD-background superfamily effect was significant, but disappears once compared against a structurally matched (same-segment/same-class) background."
        else:
            category = "MIXED"
            rationale = "Neither the original nor the matched superfamily result reaches FDR significance, or the evidence is otherwise inconclusive."

        rows.append({
            "functional_site": site,
            "original_superfamily_significant": orig_sig,
            "matched_superfamily_status": matched_status,
            "matched_superfamily_significant": matched_sig,
            "n_individual_results_with_both_tests": n_with_both,
            "n_individually_significant_original": n_orig_sig_individual,
            "n_individually_significant_preserved": n_preserved,
            "preservation_fraction": preservation_fraction,
            "interpretation_category": category,
            "rationale": rationale,
        })

    classification = pd.DataFrame(rows)
    out = TABLES / "functional_site_interpretation_classification.csv"
    classification.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(classification)} site classes)")
    print(classification[["functional_site", "interpretation_category", "rationale"]].to_string(index=False))


if __name__ == "__main__":
    main()
