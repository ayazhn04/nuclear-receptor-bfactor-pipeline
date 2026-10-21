# Final Statistical Consistency Audit

> 29 claims cross-checked: every value below was INDEPENDENTLY RECOMPUTED from its machine-readable source table (never hand-typed) and confirmed to appear correctly in the corresponding Markdown summary. **29/29 PASS, 0 FAIL.**

| Claim | Recomputed value | Target file | Status |
|---|---|---|---|
| PC1 explained variance | `31.4` | reports/tables/modeling_summary.md | ✅ PASS |
| DBSCAN max clusters found (PCA-primary) | `**1**` | reports/tables/modeling_summary.md | ✅ PASS |
| Best silhouette value (PCA-primary) | `silhouette=0.219` | reports/tables/modeling_summary.md | ✅ PASS |
| Subfamily PRIMARY-grouping FDR-significant segment count | `significant segments: 0 (PRIMARY)` | reports/tables/subfamily_flexibility_summary.md | ✅ PASS |
| Atlas most-flexible well-supported segment (rank 1) | `H2` | reports/tables/subfamily_flexibility_summary.md | ✅ PASS |
| apo/holo PRIMARY-eligible receptor count | `PRIMARY (>=5 APO and >=5 HOLO PDB entries): 7/42 receptors eligible` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo EXPLORATORY-eligible receptor count | `EXPLORATORY (>=3/>=3): 9/42 eligible` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo within-receptor FDR significant count | `**11** significant within-receptor FDR` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo global FDR significant count | `**0** significant under the SECONDARY global FDR` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3A1/H1 | `+0.656` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3A1/H3 | `-0.123` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3A1/H4 | `+0.410` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3A1/H9 | `+0.229` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR2B1/H12 | `-0.807` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR4A1/H7 | `-0.163` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR4A1/H3 | `-0.194` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR1C3/H9 | `+0.068` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3B3/H9 | `+0.210` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3B3/H8 | `+0.162` | reports/tables/apo_holo_summary.md | ✅ PASS |
| apo/holo delta for NR3B3/H6 | `-0.569` | reports/tables/apo_holo_summary.md | ✅ PASS |
| Functional-site within-receptor FDR significant count | `**82**` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for AF2_H12 | `+0.257` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for COACTIVATOR_GROOVE | `-0.369` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for DIMERIZATION_INTERFACE | `-0.091` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Superfamily site-class delta for LIGAND_BINDING_POCKET | `-0.386` | reports/tables/current_functional_site_summary.md | ✅ PASS |
| Interpretation classification for AF2_H12 | `SITE_SPECIFIC_SUPPORTED` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
| Interpretation classification for COACTIVATOR_GROOVE | `STRUCTURAL_CONTEXT_EXPLAINED` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
| Interpretation classification for DIMERIZATION_INTERFACE | `INSUFFICIENT_RESIDUE_LEVEL_MAPPING` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
| Interpretation classification for LIGAND_BINDING_POCKET | `SITE_SPECIFIC_SUPPORTED` | reports/tables/functional_site_specificity_summary.md | ✅ PASS |
