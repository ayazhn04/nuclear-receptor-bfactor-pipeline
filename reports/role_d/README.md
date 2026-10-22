# Role D — biological evidence and analysis

This is the single entry point for participant D's new work.

## Literature review

[`literature_review/`](literature_review/) contains the systematic biological
review, bibliography, PubMed search log, and the complete 48 receptor × 4 site
evidence audit. It keeps direct residue-level literature evidence separate from
project geometry and structural-proxy definitions.

## Contextual IQR analysis

[`contextual_iqr/`](contextual_iqr/) contains generated tables and the figure
for the v2 between-PDB IQR analysis and its LBP matched-context comparison. It
is computational analysis, not literature-review evidence.

## Code and validation

- `scripts/109_role_d_position_iqr.py`
- `scripts/110_role_d_position_iqr_figures.py`
- `scripts/111_role_d_lbp_matched_iqr.py`
- `tests/test_role_d_position_iqr.py`
- `tests/test_role_d_literature_audit.py`

Historical Role D inputs remain in `reports/tables/functional_site_definitions.csv`
and `reports/ROLE_D_PROGRESS_REPORT.md`, since the historical analysis script
uses those stable paths.
