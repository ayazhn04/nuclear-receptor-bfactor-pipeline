# Final Statistical Consistency Audit

> 21 claims cross-checked: every value below was INDEPENDENTLY RECOMPUTED from its machine-readable source table (never hand-typed) and confirmed to appear correctly in the corresponding Markdown summary. **21/21 PASS, 0 FAIL.**

| Claim | Recomputed value | Target file | Status |
|---|---|---|---|
| PC1 explained variance | `28.4` | reports/tables/modeling_summary.md | ✅ PASS |
| DBSCAN max clusters found (PCA-primary) | `**1**` | reports/tables/modeling_summary.md | ✅ PASS |
| Best silhouette value (PCA-primary) | `silhouette=0.410` | reports/tables/modeling_summary.md | ✅ PASS |
| Subfamily PRIMARY-grouping FDR-significant segment count | `significant segments: 0 (PRIMARY)` | reports/tables/subfamily_flexibility_summary.md | ✅ PASS |
| Atlas most-flexible well-supported segment (rank 1) | `loop_H11_H12` | reports/tables/subfamily_flexibility_summary.md | ✅ PASS |
| apo/holo PRIMARY-eligible receptor count | `PRIMARY (>=5 APO and >=5 HOLO PDB entries): 7/42 receptors eligible` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo EXPLORATORY-eligible receptor count | `EXPLORATORY (>=3/>=3): 9/42 eligible` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo within-receptor FDR significant count | `**3** significant within-receptor FDR` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo global FDR significant count | `**0** significant under the SECONDARY global FDR` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3A1/H1 | `+0.656` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3A1/H3 | `-0.145` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3B3/loop_H6_H10 | `+0.138` | reports/tables/apo_holo_summary.md | ✅ PASS |
| Functional-site within-receptor FDR significant count | `**73**` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for AF2_H12 | `+0.264` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for COACTIVATOR_GROOVE | `-0.340` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for DIMERIZATION_INTERFACE | `+0.091` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for LIGAND_BINDING_POCKET | `-0.353` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Interpretation classification for AF2_H12 | `SITE_SPECIFIC_SUPPORTED` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
| Interpretation classification for COACTIVATOR_GROOVE | `STRUCTURAL_CONTEXT_EXPLAINED` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
| Interpretation classification for DIMERIZATION_INTERFACE | `INSUFFICIENT_RESIDUE_LEVEL_MAPPING` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
| Interpretation classification for LIGAND_BINDING_POCKET | `SITE_SPECIFIC_SUPPORTED` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
