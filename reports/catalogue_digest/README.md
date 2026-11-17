# Catalogue digest (meeting companion)

Six pages that break the two "Dataset at a glance" slides of the receptor catalogue into
one question per page. It is a separate document; the receptor catalogue
(`reports/catalogue_v2/`) is not changed.

## Source

One table, read as is: `reports/catalogue_v2/tables/receptor_catalogue_summary_v2.csv`
(commit `c0de7a5`, the same table the glance slides are drawn from). No definition is
recomputed. Every value on a page is a column of that table, or a share or sum of its
columns; `digest_values.csv` lists each value per receptor and page. The script checks the
statements on page 6 against the table before it draws them.

Rebuild: `python scripts/114_build_catalogue_digest.py`

## Pages

| Page | Question | Columns used |
|---|---|---|
| 1 | How much usable structural data is there per receptor? | `n_entries_total`, `n_entries_with_zscores` |
| 2 | Which receptors support an APO vs HOLO comparison? | `n_holo`, `n_apo`, `n_ambiguous`, `n_trace_lbd`, `n_entries_with_zscores`, `apo_holo_tier` |
| 3 | Ligand diversity is very uneven | `n_unique_ligands`, `n_entries_with_zscores` |
| 4 | Are structures mostly standalone or partner-bound? | `frac_with_partner`, `n_with_partner`, `n_partner_nr_primary`, `n_partner_coregulator_peptide_primary` |
| 5 | Between-PDB flexibility heterogeneity differs by receptor | `median_iqr_z`, `heterogeneity_class`, `heterogeneity_cut_*` |
| 6 | What do these variables tell us together? (4 receptor pairs) | the columns above |

## Definitions

- **Entry**: receptor × PDB entry. **Total** = every entry found; **z-scored** = at least one
  finite B-factor z-score (the catalogue's "usable PDBs").
- **APO / HOLO / ambiguous / TRACE_LBD**: catalogue state of each z-scored entry. TRACE_LBD
  fragments (short receptor peptides) are drawn hatched or listed separately and never
  count as APO or HOLO. "In the APO/HOLO analysis" = `apo_holo_tier` PRIMARY (filled dot)
  or EXPLORATORY (open dot).
- **Unique ligands**: distinct functional-ligand CCD IDs in the z-scored entries.
- **With a partner chain**: share of z-scored entries with any polymer chain within 5 Å of
  the LBD. Split into one class per entry: NR / heterodimer partner first, then
  coregulator peptide, then any other partner.
- **Median IQR(z)**: median over canonical positions seen in ≥ 5 entries of the IQR of
  per-entry B-factor z. It is between-PDB crystallographic heterogeneity: not motion in
  solution, and not controlled for burial. Low / moderate / high = the catalogue's tertiles
  among the 30 assessable receptors.

## Display choices (no data removed)

- Page 1: two panels, ≥ 20 and < 20 z-scored entries; the second has its own zoomed
  scale. RORβ, TR2 and GCNF have no entry found and are listed in the footnote.
- Page 2: bars only for receptors with both APO and HOLO entries; HOLO-only, APO-only and
  no-APO/HOLO receptors are side lists with their counts.
- Page 3: bars for receptors with ≥ 5 unique ligands (dark = ≥ 50); 1–4 and 0 are side
  lists. The 3 receptors without z-scored entries (DAX-1, COUP-TFI, NOR-1) and the 3 with
  no entry found are not shown.
- Page 4: bars for receptors with ≥ 5 z-scored entries; the 11 with 1–4 are listed as
  "partner / z-scored" counts without ranking.
- Page 5: the 30 receptors with an assessable IQR (≥ 20 positions seen in ≥ 5 entries);
  the other 18 are listed as not assessable.
- Page 6: values only; no statement implies cause. "Associated with" / "does not map
  simply onto" wording only.

## Files

| File | Content |
|---|---|
| `catalogue_digest.pdf` | the 6 pages |
| `figures/digest_<n>_<name>.png` | each page, 1920 × 1080 |
| `figures/chart_only/` | the same pages without the title band (used in the editable deck) |
| `digest_values.csv` | every value shown, per page and receptor |
| `digest_page_text.csv` | each page's title and one-line message |
