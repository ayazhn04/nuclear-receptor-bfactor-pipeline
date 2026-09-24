# `config/` — Nuclear Receptor Metadata & Validation Configuration

This directory holds the project's reference metadata for the human nuclear
receptor (NR) superfamily and the configuration that governs how that
metadata is validated in later stages.

## Provenance model

1. **Professor-supplied mapping is the starting reference.**
   The course supervisor supplied a mapping of 48 human nuclear receptors
   (names, NR nomenclature codes, UniProt accessions, search terms,
   project-level groupings). That mapping was preserved **immutably**, byte
   for byte, in [`nr_metadata_professor.py`](nr_metadata_professor.py)
   during Stage 1.

2. **The professor-supplied values are not ground truth to be silently
   overwritten.** They are the project's starting reference point. Stage 1
   will independently validate every entry against current authoritative
   sources (e.g. UniProt, HGNC, NR nomenclature committee resources).

3. **Both the original and validated values must remain traceable.**
   When validation occurs, the pipeline must retain the professor-supplied
   value *and* the database-confirmed value side by side (e.g. as separate
   fields or columns), not overwrite one with the other.

4. **Never silently change biological identifiers.**
   Receptor names, NR nomenclature codes (e.g. `NR3C1`), UniProt accessions,
   and structural/functional classifications must never be modified without
   an explicit, recorded justification. "Looks close enough" or fuzzy
   text-matching is not sufficient justification.

5. **Every discrepancy must be logged.**
   Any correction, ambiguity, or disagreement discovered while validating
   against authoritative databases must be written to an explicit audit /
   change log (introduced in Stage 1), including: what was supplied, what
   the database says, the source consulted, and the resolution.

6. **No inference from filenames or fuzzy text matching when an
   authoritative identifier is available.** If a UniProt accession, PDB ID,
   or other authoritative identifier exists for a receptor, it must be used
   directly rather than guessed from a gene name, filename, or free text.

## Files

- [`nr_metadata_professor.py`](nr_metadata_professor.py) — Immutable,
  verbatim snapshot of the professor-supplied 48-receptor list. Never
  hand-edited after creation; any future correction must go through a
  documented, re-run validation step, not a silent edit here.
- [`nr_metadata.py`](nr_metadata.py) — The project's validated working
  receptor master. **Generated**, not hand-typed, by
  `scripts/03_freeze_validated_metadata.py`. All 48 UniProt accessions were
  independently checked against the current UniProt REST API during Stage 1;
  all 48 resolved to current, reviewed (Swiss-Prot) *Homo sapiens* entries,
  and no primary-accession replacements were required. `nr_metadata.py`
  therefore now contains the full 48-record working master, with `uniprot_id`
  set to the validated current primary accession (identical to the
  professor-supplied accession for every record) and `nr_code`,
  `common_name`, `search_terms`, and `group` preserved exactly as supplied.

## Stage boundary

`NR_METADATA` in `nr_metadata.py` is **populated** (48 records) as of Stage
1.5. Detailed audit and provenance — the full per-record UniProt validation
results, checksums, and generation metadata — live under
[`../data/manifests/`](../data/manifests/), principally
`nr_metadata_uniprot_audit.csv`, `nr_metadata_uniprot_audit_summary.json`,
`nr_metadata_validated.csv`, and `nr_metadata_validated_provenance.json`.

This does not relax the provenance rules above: any future biological
correction (e.g. arising from RCSB/structural-stage discoveries) must still
never be applied silently. It must be documented as an explicit,
auditable change — never a direct hand-edit of `nr_metadata.py` or
`nr_metadata_professor.py`.
