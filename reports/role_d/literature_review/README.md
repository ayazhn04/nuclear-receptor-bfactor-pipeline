# Role D — literature review and curated functional-site evidence

This directory makes Role D’s manual biological evidence review visible in Git. It does **not** change the historical table consumed by analysis code.

## Deliverables

- [Systematic biological literature review](ROLE_D_SYSTEMATIC_LITERATURE_REVIEW.md) — functional-site biology, five anchor receptors, explicit cross-receptor limits, and the crystallographic B-factor/solution-dynamics distinction.
- [Machine-readable bibliography](bibliography.csv) — verified DOI/PMID/URL identifiers, source type, scope/site, and claim supported.
- [Evidence coverage and gaps](evidence_coverage_and_gaps.csv) — residue-level primary evidence versus curated-resource support, geometry-derived definitions, and segment proxies.
- [All-receptor × site evidence matrix](per_receptor_site_evidence_matrix.csv) — one auditable literature-evidence decision for each of the 48 project receptors and four current project site classes (192 unique cells). It separately records literature tier, residue-level coverage, verified citation key(s), the current project-definition status, and the exact gap/reason.
- [PubMed systematic-search log](pubmed_systematic_search_log.csv) — the 192 receptor × site discovery queries, date, verbatim query, retrieval count, and request status. Retrieval is a discovery record, **not** evidence adjudication.
- [Historical curated source table](../../tables/functional_site_definitions.csv) — stable 44-position input read by the historical Role D script; **read-only for this review**.

## Evidence boundary

The historical source table contains 44 canonical positions for five anchors: VDR, PPARγ, ERα, GR and NR4A1. It is a five-anchor curated set, **not** residue-level literature coverage for every receptor in the 48-receptor project. Outside those curated assignments, current v2 labels may be H12 mapping, ligand-contact geometry, or structural-segment proxies. See the coverage table before interpreting any project-wide label as a direct functional contact.

## All-48 audit completion (2026-10-22)

The evidence matrix is complete as an **audit**: 48 metadata receptors × AF2_H12, LIGAND_BINDING_POCKET, COACTIVATOR_GROOVE and DIMERIZATION_INTERFACE = **192/192 unique cells**. Its literature result is deliberately sparse rather than inferred: 9 cells have verified receptor-specific primary residue evidence, 2 have curated-resource-only residue/sequence support, 1 is the limited/non-standard NR4A1 cytosporone-B annotation, and 180 have no verified receptor-specific residue citation entered by this audit. The latter label does **not** prove that no paper exists.

Project-definition availability is a separate structural-data field: 64 cells are `AVAILABLE`, 74 `PARTIAL`, 30 `UNRESOLVED`, and 24 `NOT_ASSESSED_NO_COORDINATE_CANDIDATES`. The last category comprises all four cells for each of six metadata receptors absent from the current coordinate-candidate availability table; it is not a biological absence claim.

Nonblank `citation_keys` in the matrix resolve to `bibliography.csv`. The only sources elevated to residue-level support are the already verified, receptor-specific anchor records; no structural proxy, homologous residue, or search-result snippet was promoted to a literature annotation.

## Search and screening boundary

On 2026-10-22, a reproducible PubMed discovery search was executed for every matrix cell (192 successful requests). Queries combine the current gene symbol and recommended protein name with site-specific AF-2/H12, ligand-pocket, coactivator/coregulator, or dimerization concepts. Exact queries and counts are retained in `pubmed_systematic_search_log.csv`. A hit count is neither a verified residue claim nor evidence of absence. A record was added to the matrix only when the cited source had been checked for receptor identity, site context, and canonical-numbering/claim compatibility. Thus, the matrix makes the remaining source-level review and numbering reconciliation workload explicit instead of fabricating all-receptor residue curation.

## Related methods and limitations

- Historical curation/numbering validation: [`../../ROLE_D_PROGRESS_REPORT.md`](../../ROLE_D_PROGRESS_REPORT.md)
- Current functional-site construction: [`../../FINAL_METHODS_AUDIT.md`](../../FINAL_METHODS_AUDIT.md), §§18–20
- v2 functional-site and Role D endpoint context: [`../contextual_iqr/`](../contextual_iqr/)
- Interpretation boundaries: [`../../FINAL_LIMITATIONS.md`](../../FINAL_LIMITATIONS.md), especially functional-site and B-factor limitations

The literature review supplies biological context; it does not perform or revise statistical analyses.
