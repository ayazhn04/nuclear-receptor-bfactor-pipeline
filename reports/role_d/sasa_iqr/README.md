# Role D: SASA control of between-PDB position IQR

New **exploratory** analysis, 29 October 2026. The available verified baseline
is documented in [IQR_VALIDATION_RU.md](IQR_VALIDATION_RU.md). The completed
SASA-adjusted results are in [RESULTS_RU.md](RESULTS_RU.md), with presentation
text in [SLIDE_TEXT_EN.md](SLIDE_TEXT_EN.md). This layer does not replace the
frozen v1 release, corrected v2 baseline, or prior contextual-IQR files.
Existing presentation files were not modified.

`../Role_D_IQR_SASA_2026-10-29.zip` is the small team handoff: four PNG figures,
Russian results, English slide text, methods, audit/summary CSVs, manifests and
new code/tests. Raw structures and large intermediate Parquets are intentionally
not included. Reproduction requires the original repository with its corrected-v2
inputs; the archive is not a standalone reconstruction of the whole project.

## Endpoints and hierarchy

- Source: **corrected v2 only**, including corrected numbering, analysis span,
  professor canonical elements, dual labels and explicit review flags.
- B-factor z is already normalized within instance. Use the approved representative
  instances, then receptor × PDB × canonical position median (never duplicate chains
  as independent observations).
- **IQR(z)** is q75−q25 across PDB-entry medians at a position. It describes
  between-deposit crystallographic-disorder heterogeneity, not solution dynamics.
  It is not equivalent to median(z), the relative-disorder location endpoint.
- New SASA control requires ≥5 PDB entries with **both** usable z and a complete
  residue for SASA. IQR is recalculated on this paired-entry subset; it is not
  silently taken from a different set of PDB structures. Original profile IQR and
  counts are retained alongside paired IQR/counts in `position_sasa_iqr.csv`.
- Inferential unit for cross-receptor summaries is **receptor**. No pooled
  residue-level significance test or independent-chain test is used.

## SASA method (SAS = solvent-accessible surface)

Shrake–Rupley, Biotite, 1.4 Å probe, 1,000 Fibonacci points per atom. ProtOr
radii; element-specific fallback for unclassified atoms (counts audited).
Each representative chain/model is calculated **in isolation**, with its entire
observed polymer construct, not a cut-out LBD. This avoids artificial exposure at
LBD boundaries. Other chains, ligands, water, symmetry mates and biological-assembly
operators are absent. This is intrinsic isolated-construct accessibility, **not**
accessibility in the ligand-bound complex or crystal lattice.

Only positive-occupancy heavy atoms are used. One nonblank altloc is chosen per
residue, preferring the established selected Cα altloc where available, otherwise
the highest mean occupancy with lexicographic tie-break; shared blank-altloc atoms
are retained. Duplicate atom names are resolved deterministically by occupancy.
No B-factor value is used for SASA atom/conformer selection.

RSA = summed residue SASA / theoretical maximum ASA for its **observed residue
type**, Tien et al. 2013. Residues lacking any expected standard heavy atom,
including backbone atoms, do not receive usable SASA/RSA. Nonstandard or mixed
residue identities are explicit unsupported rows. Such atoms still occlude
neighboring atoms. RSA is not clipped at 1; geometric/reference differences can
produce values above 1. Missing residues are not assigned zero exposure.

All restored coordinates live in a separate ignored cache. Before use, the
**decompressed SHA256 must match** the original coordinate-download manifest
(including the v2 additions). Revisions, missing hashes, download failures and
calculation failures are excluded from SASA only and reported explicitly; master
membership is never changed. See `coordinate_cache_audit.csv`,
`calculation_audit.csv`, `sasa_manifest.json`.

## Exploratory statistical design

1. Spearman RSA–median(z) and RSA–IQR(z) correlations separately within receptors
   with ≥20 eligible positions. Report median receptor rho and percentile bootstrap
   CI (10,000 receptor resamples, seed 20261029). Two-sided Wilcoxon against zero,
   BH correction across the two endpoints. Spatial dependence can affect each rho;
   receptor resampling does not make residues independent or remove confounding.
2. LBP site versus positions outside the **two analyzed site classes** (LBP and
   AF2/H12), within exact professor-context and PDB-support strata (5–9, 10–19,
   ≥20). “Nonfunctional” is shorthand for this exclusion, **not** proof that those
   control positions lack biological function. Inter-element strata are matching
   devices, not newly assigned secondary structure. Dual labels are retained as
   an exact composite context, conservatively rather than double-counted.
3. Add median-position RSA bands ≤0.05 (buried), (0.05, 0.25] (intermediate),
   >0.25 (exposed). These prespecified coarse operational bands are not universal
   biological categories. Require ≥1 site and ≥3 control positions per stratum.
   Per-stratum median IQR difference → median across strata → one receptor effect.
4. Sensitivities: exclude professor review-required positions; additionally match
   canonical amino-acid identity (a conservative, potentially sparse control).
   All four designs share a BH family for receptor-level Wilcoxon p-values;
   bootstrap CIs are pointwise, not multiplicity-adjusted. CIs and p-values may
   disagree because they summarize different aspects of the effect distribution.
5. AF2/H12 gets descriptive IQR/RSA summaries only: all of H12 defines that site,
   leaving no independent same-H12 non-site background. No circular test.

These are position-level association analyses, not causal ligand-effect tests,
not APO/HOLO comparisons, and not a substitute for resolution, refinement,
packing, mutation, partner or pharmacology controls. Matching can substantially
change the cohort; differences between designs are **not** a mediation estimate.
Coarse exposure bins and the use of median RSA do not capture all accessibility
variation. Missing sidechains/coordinates can be nonrandom.

Numerical pilot: `sphere_point_convergence_pilot.csv` compares 1,000 versus
4,000 sphere points on PR/1A28. For 236 complete residues the median absolute RSA
difference was 0.00101, maximum 0.00737, and all exposure-band assignments agreed.
This is a single-structure numerical check, not whole-cohort convergence proof.

## Reproduce

Python dependencies: numpy, pandas, scipy, matplotlib, pyarrow, biopython,
biotite, requests; pytest for tests. Exact scientific versions are recorded in
`sasa_manifest.json`. No change to the legacy conda environment is required.

```bash
python scripts/113_role_d_sasa.py --download --workers 8
python scripts/114_role_d_sasa_iqr_analysis.py
python -m pytest -q tests/test_role_d_sasa_iqr.py
```

Omit `--download` for offline reproduction once the coordinate cache is restored.
The scripts write only this report directory and the separate ignored raw cache.
`representative_instance_sasa.parquet` retains QC statuses including missing
residues; `pdb_position_sasa.parquet` contains usable PDB-level medians.
`analysis_manifest.json` records input hashes, thresholds and random seed.

## References

- Shrake A, Rupley JA (1973). Environment and exposure to solvent of protein
  atoms. *J Mol Biol* 79:351–371. DOI: 10.1016/0022-2836(73)90011-9.
- Tsai J et al. (1999). The packing density in proteins: standard radii and
  volumes. *J Mol Biol* 290:253–266. DOI: 10.1006/jmbi.1999.2829.
- Tien MZ et al. (2013). Maximum allowed solvent accessibilities of residues
  in proteins. *PLoS ONE* 8:e80635. DOI: 10.1371/journal.pone.0080635.
