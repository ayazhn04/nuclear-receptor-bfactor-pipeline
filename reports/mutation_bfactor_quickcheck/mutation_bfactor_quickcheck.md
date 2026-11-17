# Engineered mutations and normalized B-factors: quick check (v2)

## Meeting answer

1. **Do engineered mutations affect B-factor? NO CLEAR EVIDENCE** of a general effect.
   Mutant and WT entries differ by about 0.02 z (negligible). A few receptors differ, but
   those differences are not centred on the mutated residues.
2. **Is any effect mainly local around the mutated position? NO.** Mutant z at the
   mutated residue, compared with WT at the same receptor position, is not higher. There
   is no peak at the site that fades with distance.
3. **Do mutant-heavy receptors tend to have higher between-PDB IQR? NO** (Spearman
   rho = 0.18, p = 0.35, 30 receptors). This is a receptor-level association only.

Observational crystallographic data: these are associations, not effects of the mutation.
z = relative crystallographic disorder within each chain.

## Data

- **Entries:** 1,805 z-scored entries, from `catalogue_mutation_audit_v2.csv`.
  - 537 `ENGINEERED_MUTANT`
  - 1,158 `WT_NO_ENGINEERED_SUBSTITUTION`
  - 110 `UNCERTAIN`: excluded, never counted as WT.
- **Mutation sites:** only depositor-confirmed `ENGINEERED` residues are used. Natural
  variants and modified residues are not used.
  - 955 engineered receptor × entry × position sites, in all 537 mutant entries.
  - 954 sites are observed in the mutant entry. The one missing is ERα 5UFX C530S, and it
    stays missing.
  - 809 sites have a WT reference at the same receptor position (≥ 5 WT entries observing
    that position).
  - 145 sites have no WT reference, so no matched comparison: GR 88, MR 37 (no WT entry at
    all), SF-1 8, TLX 4, ERRβ 3, TR4 2, and 1 each for ERRα, HNF4α and PNR.
- **z values:** the v2 entry-level `b_factor_z` (per-chain z, median over representative
  chains) at corrected canonical UniProt positions. Raw B-factors are never compared.

## 1. Entry-level distribution (mutant vs WT)

Chain mean z is 0 by construction, so it is not compared. Per entry: median, IQR and
90th percentile of z over the LBD positions. The 19 TRACE_LBD-only fragments are excluded.

| Measure | Mutant (537) | WT (1,139) | Difference | Cliff's delta | Mann-Whitney p | Receptor-stratified p (van Elteren) |
|---|---|---|---|---|---|---|
| Median z | −0.234 | −0.254 | +0.020 | 0.14 | 5e-6 | 8e-7 |
| IQR of z | 1.170 | 1.153 | +0.017 | 0.06 | 0.047 | 0.013 |
| 90th percentile z | 1.375 | 1.346 | +0.029 | 0.08 | 0.011 | 0.46 |

The differences are about 0.02–0.03 z, and every Cliff's delta is below 0.15
(negligible). The p-values are small only because of the large N. Median resolution is
similar (2.16 Å mutant vs 2.20 Å WT).

**Within receptor** (10 receptors with ≥ 5 mutant and ≥ 5 WT entries; Holm-corrected
within each measure; Figure 2). Only three results pass:
- **FXR, median z:** +0.11, delta 0.67, Holm p = 0.001 (21 mutant vs 25 WT).
- **ERα, 90th percentile z:** +0.22, delta 0.43, Holm p = 0.0007. ERα IQR also passes:
  +0.18, delta 0.35, Holm p = 0.012. ERα has 309 mutant vs only 31 WT entries.
- **ERβ, 90th percentile z:** −0.33, delta −0.83, Holm p = 0.004, from only 6 mutants,
  all Cys→Ser.

AR, RORγ, PPARγ, PXR, TRβ, LXRβ and NUR77 show no difference after correction. The
receptor-level shifts may reflect what else differs between mutant and WT entries
(ligand series, construct, crystal form), and they are not located at the mutated
residues (section 2).

## 2–3. Local, position-matched effect

For each engineered site p and offset d, the comparison is
**diff = z_mutant(p + d) − median z of WT entries of the same receptor at p + d.** Each
mutation observation is averaged within a distance bin, then the bins are summarized
across observations. The test is a Wilcoxon signed-rank vs 0, the effect size is the
matched-pairs rank-biserial r, and the CI is a bootstrap 95% CI of the median (2,000
resamples).

| Distance from site | N | Median diff | 95% CI | p | r |
|---|---|---|---|---|---|
| 0 (the residue itself) | 809 | −0.047 | −0.093 to −0.005 | 0.91 | 0.00 |
| 1–2 | 810 | +0.029 | −0.008 to 0.068 | 0.0008 | 0.14 |
| 3–5 | 811 | −0.010 | −0.042 to 0.017 | 0.69 | 0.02 |
| 6–10 | 811 | −0.061 | −0.085 to −0.027 | 0.005 | −0.11 |
| 11–20 | 811 | −0.019 | −0.039 to 0.001 | 0.64 | −0.02 |
| window ±2 | 810 | +0.015 | −0.033 to 0.052 | 0.009 | 0.11 |
| window ±5 | 811 | +0.004 | −0.029 to 0.041 | 0.10 | 0.07 |
| window ±10 | 811 | −0.015 | −0.053 to 0.021 | 0.83 | −0.01 |

The mutated residue itself is not more disordered than WT at that position. There is no
gradient that fades with distance: the 6–10 bin is as far from 0 as the 1–2 bin, in the
opposite direction. All medians are under 0.07 z.

The small positive 1–2 signal does not hold up across sensitivity analyses:
- **ERα dominates:** it has 544 of the 809 sites, 215 of them Y537S, so repeated
  mutations carry the pooled profile (Figure 1).
- **Isolated sites** (no other engineered site within ±10): 1–2 bin +0.069 (CI 0.029 to
  0.120), r = 0.26. This is still small.
- **Excluding ERα** (265 sites, 16 receptors): every bin's median is between −0.07 and
  +0.03, and every p ≥ 0.09. At the residue itself the median is −0.071 (CI −0.160 to
  −0.001), p = 0.09.
- **One value per receptor × mutation** (117 distinct mutations, 17 receptors): every
  bin's |median| ≤ 0.04 and every p ≥ 0.6.

**Receptor examples, window ±2** (receptors with ≥ 10 sites; Holm-corrected):
- **ERα:** +0.048 (CI 0.001 to 0.099), Holm p = 0.0009, r = 0.19. Y537S, the most
  frequent mutation (n = 215), gives +0.13. L536S gives 0.00 and C381S −0.17.
  Significant but small.
- **LXRβ:** −0.36 (CI −0.75 to 0.09), Holm p = 0.31. Its surface mutations (R263S,
  D261S, Q258A) sit lower than WT; not significant.
- **PPARγ:** −0.24, and **RORγ:** −0.31. Neither is significant after correction.
- **AR, FXR, ERβ:** no consistent local difference.

## 4. Mutant-heavy receptors and between-PDB IQR

Across the 30 receptors with an assessable IQR, the engineered-mutant fraction of
z-scored entries vs receptor median position IQR(z) gives Spearman rho = 0.18, p = 0.35.
With UNCERTAIN removed from the denominator, rho = 0.21, p = 0.26.
- **Mutant-heavy receptors** span the whole IQR range: ERα 88% mutant, IQR 0.46 (low end
  of moderate); GR 82%, 0.68 (high); MR 76%, 0.58 (high); SF-1 67%, 0.38 (low).
- **Several high-IQR receptors** are mostly WT: RXRα, RARα, LRH-1 (Figure 3).

This is ecological, receptor-level evidence only.

## Caveats

- **What mutants differ by:** mutant and WT entries also differ in ligand series,
  constructs and crystal forms. Position matching controls for the region of the
  receptor, not for those.
- **Imbalance:** the comparison is uneven. ERα has 309 mutants vs 31 WT, MR has no WT
  entry, and GR has 6.
- **Repeated mutations** (ERα Y537S, L536S, C381S) dominate any observation-level pooling.
  The per-mutation analysis is the fairer summary.

## Files

| File | Content |
|---|---|
| `mutation_bfactor_quickcheck_results.csv` | every test: Q1 pooled and per receptor, Q2/Q3 bins (4 variants), per-receptor local, Q4 |
| `run_metadata.csv` | site counts and exclusions |
| `fig1_local_effect_by_offset.png` | median position-matched difference by signed distance, all vs excluding ERα |
| `fig2_entry_level_within_receptor.png` | Cliff's delta per receptor for median, IQR, 90th percentile z |
| `fig3_mutant_fraction_vs_iqr.png` | receptor mutant fraction vs median IQR(z) |

Rebuild: `python scripts/115_mutation_bfactor_quickcheck.py` (reads the catalogue tables
and the v2 profile parquet; writes only this folder).
