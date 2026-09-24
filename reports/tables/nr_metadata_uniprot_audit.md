# Nuclear Receptor Metadata — Current UniProt Audit

## Methods

Every UniProt accession in the professor-supplied `NR_METADATA` snapshot (`config/nr_metadata_professor.py`) was queried against the official UniProt REST API (`https://rest.uniprot.org/uniprotkb`) and compared against the professor-supplied values. Comparisons were conservative: name differences alone were never treated as proof of a mismatch, and no professor-supplied value was altered. See `data/manifests/nr_metadata_uniprot_audit.csv` for the full per-record table and `data/manifests/uniprot_api_manifest.csv` for the raw fetch audit trail.

## Overall counts

- Professor records: 48
- Successful API responses: 48
- Failed API responses: 0
- PASS: 48
- REVIEW_NEEDED: 0
- FAIL: 0
- Human / taxon 9606: 48
- Reviewed (Swiss-Prot) entries: 48
- Primary accession differs from professor accession: 0

## PASS / REVIEW_NEEDED / FAIL summary

| Status | Count |
|---|---|
| PASS | 48 |
| REVIEW_NEEDED | 0 |
| FAIL | 0 |

## Non-PASS records

None.

## Unexpected findings

No accessions resolved to a different current primary accession.

## Data integrity statement

No professor-supplied value (`uniprot_id`, `nr_code`, `common_name`, `search_terms`, `group`) was modified, overwritten, or silently corrected during this audit. `config/nr_metadata_professor.py` and `config/nr_metadata.py` (still `NR_METADATA = []`) were not touched by this script beyond being read. All discrepancies are recorded side-by-side in the audit table for manual review.
