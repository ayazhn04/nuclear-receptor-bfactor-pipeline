# Slide text — IQR and solvent exposure (29 October 2026)

## Slide 1 — What was calculated

- Corrected v2 numbering and professor structural-element annotations.
- 7,065 positions in 31 receptors; ≥5 paired PDB entries per position.
- Between-PDB IQR of within-instance-normalized Cα B-factor z, not solution dynamics.
- Isolated-chain SASA (Shrake–Rupley; 1.4 Å probe); complete residues only; Tien-normalized RSA.
- Figure: 04_functional_site_iqr_rsa_heatmap.png (descriptive).

## Slide 2 — Solvent exposure is a potential confounder

- IQR_Z: median within-receptor Spearman rho 0.414 (receptor-bootstrap 95% CI 0.368 to 0.459; n=30; BH q=1.86e-09).
- MEDIAN_Z: median within-receptor Spearman rho 0.560 (receptor-bootstrap 95% CI 0.536 to 0.576; n=30; BH q=1.86e-09).
- Figure: 01_sasa_correlations.png; each point is one receptor.

## Slide 3 — Ligand-pocket contrast after exposure matching

- Context + support + RSA matched LBP − controls: median ΔIQR 0.0020 (95% CI -0.0523 to 0.0333; n=29 receptors; BH q=1).
- No convincing LBP-specific contrast under the combined CI and BH criterion; not evidence of equivalence.
- Figure: 02_sasa_adjusted_lbp_iqr.png; cohort/support changes across matching designs.

## Slide 4 — Interpretation and limitations

- IQR is heterogeneity across deposited structures; median z is relative disorder within structures.
- Exposure matching is coarse and exploratory; packing, resolution, refinement and composition remain potential confounders.
- AF2/H12: descriptive only here; no independent non-site same-H12 background.
- No causal ligand-induced mechanism or solution-dynamics conclusion is established.
