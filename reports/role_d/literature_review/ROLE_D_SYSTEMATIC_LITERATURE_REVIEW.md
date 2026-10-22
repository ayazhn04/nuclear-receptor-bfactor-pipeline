# Role D systematic biological literature review: nuclear-receptor LBD functional sites

**Scope and status.** This review documents the biological evidence underlying Role D functional-site interpretation. It is a literature and curation companion to, not a replacement for, the immutable historical table [`../tables/functional_site_definitions.csv`](../tables/functional_site_definitions.csv). The review covers (i) the ligand-binding pocket (LBP), AF-2/helix 12 (H12) and coactivator-facing surface, and LBD dimerization; (ii) the five manually curated anchor receptors (VDR, PPARγ, ERα/ESR1, GR/NR3C1 and NR4A1/NUR77); and (iii) the appropriately limited framework for transferring *structural* site labels across receptors. Citation metadata and the claim-to-source map are in [`bibliography.csv`](bibliography.csv); residue-level coverage and gaps are in [`evidence_coverage_and_gaps.csv`](evidence_coverage_and_gaps.csv). The completed all-48 audit is [`per_receptor_site_evidence_matrix.csv`](per_receptor_site_evidence_matrix.csv), with the discovery-search record in [`pubmed_systematic_search_log.csv`](pubmed_systematic_search_log.csv).

## Review protocol and evidence rules

### Question and units of evidence

The question was: *which LBD positions or regions have direct published evidence for ligand recognition, AF-2/coactivator recognition, or LBD dimerization, and what may instead be called only a structural proxy?* The unit of curation is a human UniProt canonical-sequence position when an assignment is present in the historical anchor table. Published construct/PDB numbering is not silently treated as canonical numbering; PPARγ positions in the table are explicitly converted from the PPARγ1 numbering reported by Nolte *et al.* to PPARγ2/UniProt P37231 (+28), as documented in [`../ROLE_D_PROGRESS_REPORT.md`](../ROLE_D_PROGRESS_REPORT.md).

Included evidence was peer-reviewed primary structural/functional work and a small number of foundational primary structural papers needed to delimit the common LBD/AF-2 architecture. Stable DOI and/or PMID metadata were verified before inclusion; curated UniProt records are identified as sequence/feature resources, not primary experimental literature. The review does **not** infer a residue-level function from sequence homology, helix name, proximity to a ligand in an unrelated receptor, or an automated 5 Å contact calculation.

### All-48 systematic discovery audit and source-screening boundary (2026-10-22)

The audit universe was the 48 unique human receptor accessions in `data/manifests/nr_metadata_validated.csv`, crossed with the four current project site classes: `AF2_H12`, `LIGAND_BINDING_POCKET`, `COACTIVATOR_GROOVE`, and `DIMERIZATION_INTERFACE` (192 cells). PubMed was queried once per cell on 2026-10-22 using the receptor's current primary gene symbol **or** recommended protein name in title/abstract plus a pre-specified site-concept block: AF-2/helix-12/activation-function-2; ligand-binding-pocket/domain; coactivator/coregulator plus LxxLL/charge-clamp/groove/binding; or dimerization/dimer-interface/homo-/heterodimer. The verbatim query, returned count and request outcome for every cell are retained in `pubmed_systematic_search_log.csv`; all 192 requests returned successfully.

Search retrieval was intentionally separated from evidence adjudication. A PubMed count is not a residue-level claim and zero hits do not prove biological absence. A citation was entered into the evidence matrix only after source-level checking of receptor identity, functional-site context, and canonical-numbering compatibility. In this pass, the verified record set remains the anchor literature already represented in the bibliography; no candidate result was converted into a residue-level annotation merely from abstract/search metadata. This is a complete coverage audit, not a claim of exhaustive full-text curation for all 192 cells.

## All-receptor × site completion and evidence counts

`per_receptor_site_evidence_matrix.csv` has exactly **192 unique (`uniprot_id`, `site_class`) rows**: 48 per site class. It links each cell to the project availability/definition data while keeping that structural evidence distinct from literature evidence.

| Literature-evidence tier | Cells | Meaning |
| --- | ---: | --- |
| `RESIDUE_LEVEL_PRIMARY` | 9 | Verified receptor-specific primary structural/functional evidence (the VDR, PPARγ, ERα and GR anchor assignments). |
| `RESIDUE_LEVEL_CURATED_RESOURCE` | 2 | Curated UniProt sequence/numbering support only (VDR coactivator annotation; NR4A1 H12 segment). |
| `RESIDUE_LEVEL_PRIMARY_LIMITED_NONSTANDARD` | 1 | NR4A1 Tyr453 cytosporone-B response; explicitly not a conventional deep-pocket map. |
| `NO_VERIFIED_RESIDUE_LEVEL_EVIDENCE_IN_AUDIT` | 180 | No receptor-specific residue citation was source-verified and entered. This is an explicit review gap, not evidence of absence. |

| Site class | Primary | Curated resource | Limited/non-standard | No verified residue evidence | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| AF2_H12 | 0 | 1 | 0 | 47 | 48 |
| Ligand-binding pocket | 4 | 0 | 1 | 43 | 48 |
| Coactivator groove | 3 | 1 | 0 | 44 | 48 |
| Dimerization interface | 2 | 0 | 0 | 46 | 48 |

| Current project availability (not literature evidence) | AF2/H12 | LBP | Coactivator | Dimerization |
| --- | ---: | ---: | ---: | ---: |
| `AVAILABLE` | 28 | 30 | 4 | 2 |
| `PARTIAL` | 0 | 0 | 37 | 37 |
| `UNRESOLVED` | 14 | 12 | 1 | 3 |
| `NOT_ASSESSED_NO_COORDINATE_CANDIDATES` | 6 | 6 | 6 | 6 |

The last row represents six metadata receptors absent from `current_functional_site_availability.csv`, because the current structural layer is limited to receptors with coordinate candidates. It must not be read as lack of a site or lack of literature.

### Evidence labels used here

- **Residue-level primary evidence:** a cited experiment/structure supports the named residue(s) in that receptor and site context.
- **Residue-level curated resource:** a canonical identity or sequence segment comes from UniProt; this establishes numbering/identity, not an experimental interaction mechanism by itself.
- **Structure/geometry-derived proxy:** an analysis label arises from the project’s current structural mapping, H12 mapping, or ligand-contact geometry. It can localize a region/contact in deposited structures, but is not a receptor-specific literature annotation.
- **Segment proxy:** a whole helix or set of helices is used because residue-level partner-contact evidence was not curated/computed. It must be interpreted as the region containing a functional surface.

These labels describe the *definition*, not the precision of a B-factor result.

## Biological basis for the three LBD site concepts

### Ligand-binding pocket

Nuclear-receptor LBDs share a helical fold but do not share a single invariant set of ligand-contact side chains or an identical cavity. Early RARγ and ERα structures established that ligand recognition and agonist/antagonist-dependent arrangements can be resolved in receptor-specific LBD structures [Renaud1995; Brzozowski1997]. PPARγ, VDR and GR structures subsequently identified their own ligand-contacting residues and pocket chemistry [Nolte1998; Rochel2000; Bledsoe2002]. Thus, “LBP” is biologically meaningful at the receptor level, whereas a pooled cross-receptor LBP is an operational class rather than one conserved residue set.

The current project-wide LBP definition is more limited and reproducible: in holo structures, a receptor heavy atom is marked when within 5.0 Å of a pre-classified functional LBD ligand; receptor-level contacts are retained by a specified contact-fraction threshold (methods audit §§18–19). This is **structure/geometry-derived** evidence of contact in the dataset, not proof that every marked residue contributes equally to affinity, specificity, activation, or solution-state dynamics. It should not overwrite the curated anchor residues.

### AF-2/H12 and the coactivator-facing surface

The C-terminal activation function 2 (AF-2) contains H12 and, in an agonist-compatible conformation, contributes to a hydrophobic coactivator-binding cleft. Surface mutagenesis of thyroid hormone receptor showed that ligand-dependent coactivator binding involves a hydrophobic cleft formed when the C-terminal helix packs against other helices [Feng1998]. ERα peptide-complex structures directly resolved coactivator recognition and antagonist-dependent disruption of that interaction [Shiau1998]. The often-used “charge clamp” terminology refers to receptor residues that help position an LxxLL-type coactivator motif; exact identities and relevance remain receptor- and construct-context dependent.

Therefore H12 is a defensible **structural element** across mapped receptors, but H12 alone is not synonymous with a complete coactivator surface, nor does it establish that every H12 residue contacts a coactivator. In the current v2 dataset AF-2/H12 is defined from the professor-supplied canonical H12 intervals, not from a literature curation of every receptor’s contact residues. This mapping supports region-level interpretation only outside anchors.

### LBD dimerization

LBD dimerization is not a universal, fixed interface. ERα structural work supports an LBD homodimer interface involving hydrophobic residues in the H10/H11 area [Brzozowski1997]. In contrast, the GR LBD structure reported a distinct dimer configuration with an intermolecular β-sheet and functional support for the observed interface [Bledsoe2002]. Heterodimeric receptors may use still other arrangements; for example, the PPARγ/RXRα complex documents an asymmetric heterodimer architecture [Gampe2000]. It is consequently not valid to call H9/H10/H11 residues a residue-level dimer interface for all nuclear receptors merely because they lie in a conventionally dimer-associated region.

The project’s cross-receptor dimerization label is consequently a segment correspondence/proxy except for ERα and GR curated residues. It indicates an LBD region potentially relevant to association, not a measured interface in every receptor or a statement about oligomeric state in solution.

## Five curated anchors

### VDR (NR1I1; P11473)

Rochel *et al.* solved a calcitriol-bound VDR LBD and provide the structural basis used for the six curated ligand-pocket residues: Tyr143, Ser237, Arg274, Ser278, His305 and His397 [Rochel2000]. These are residue-level ligand-contact annotations in the historical table. Lys246 and Glu420 are curated as AF-2/coactivator-groove charge-clamp residues from the UniProt canonical record; their sequence identities and numbering are well defined, but this review does not elevate that resource annotation to a receptor-specific coactivator-peptide contact experiment. No VDR LBD dimerization residue is manually curated in the source table.

### PPARγ (NR1C3; P37231)

Nolte *et al.* described ligand binding and coactivator assembly in PPARγ [Nolte1998]. The anchor table records Ser317, His351, His477 and Tyr501 as ligand-pocket residues and Lys329/Glu499 as AF-2 charge-clamp residues in PPARγ2/UniProt numbering. The paper’s PPARγ1 numbering is the reason for the explicit +28 conversion noted above. These six assignments are direct, receptor-specific structural evidence after documented numbering reconciliation; they should not be generalized to PPARα/δ or to other NR1 receptors residue-for-residue.

### ERα (NR3A1; P03372)

Brzozowski *et al.* resolved agonist and antagonist ERα LBD states and supports the curated estradiol-recognition residues Glu353, Arg394 and His524 as well as the hydrophobic LBD homodimer-interface residues Leu504, Leu508 and Leu511 [Brzozowski1997]. Shiau *et al.* resolved ERα/coactivator recognition and supports the curated charge-clamp residues Lys362 and Glu542 [Shiau1998]. Together these are the most complete anchor-level coverage in the historical table: LBP, coactivator surface and a defined LBD homodimerization interface. They do not establish the same interface or charge-clamp positions in ERβ or non-estrogen receptors.

### GR (NR3C1; P04150)

Bledsoe *et al.* solved dexamethasone- and coactivator-motif-bound GR LBD, reported a distinct ligand pocket, coactivator recognition, and a novel dimerization mode [Bledsoe2002]. The table’s Asn564, Gln570, Arg611, Gln642 and Thr739 are GR steroid-recognition residues; Lys579 and Glu755 define the curated coactivator charge clamp; Pro625 and Ile628 are curated D-loop dimerization residues. This is receptor-specific evidence for the observed experimental construct and conditions. It does not imply that an observed GR LBD dimer is the only biologically relevant GR oligomeric arrangement in all cellular contexts.

### NR4A1/NUR77 (NR4A1; P22736)

NR4A-family LBDs require special handling. Nurr1 structural work identified a class of ligand-independent receptors with an occluded/atypical ligand-binding cavity [Wang2003], so it is not scientifically safe to assign NR4A1 the conventional deep-pocket residue list used for lipid or steroid receptors. Zhan *et al.* reported cytosporone B binding/agonism for Nur77 and the source table retains Tyr453 as a **cytosporone-B response residue** [Zhan2008]. That designation should remain as written; it is not evidence that Tyr453 defines a generic deep LBP across NR4A1 structures.

The current table also contains the canonical H12 sequence segment Pro584–Thr595 from UniProt P22736. This is a sequence/structural-element definition, not a published residue-by-residue coactivator-contact map. The appropriate Role D interpretation is therefore “H12 segment” rather than inferred standard AF-2 charge-clamp contacts or conventional pocket chemistry.

## Bounded cross-receptor evidence framework

The project contains 48 receptors, but this review makes **no claim that all 48 have residue-level literature evidence**. Cross-receptor statements are restricted as follows.

| Functional label | What may be said across mapped receptors | What may not be said without receptor-specific evidence |
| --- | --- | --- |
| AF-2/H12 | The mapped C-terminal H12 region is an LBD structural element and a biologically motivated AF-2-region proxy. | Every residue directly binds a coactivator; H12 has the same conformational or transcriptional consequence in every receptor. |
| Ligand-binding pocket | The geometric workflow marks residues observed near classified LBD ligands under its explicit 5 Å/contact-fraction rule. | Each marked residue is a conserved literature-supported affinity/activation determinant; apo receptors have a demonstrable pocket contact. |
| Coactivator groove | H3/H4/H5 correspondence can identify a region that may contain a coactivator-facing surface. | A segment-level label identifies peptide-contact residues or a receptor-specific charge clamp. |
| Dimerization interface | H9/H10/H11 correspondence can identify a region plausibly associated with LBD association in some receptor classes. | All receptors dimerize through that region, have the same oligomeric state, or retain an ERα/GR-like residue interface. |

This framework corresponds to the current methods audit (§18): tier A is read-only curated Role-D residues; tier B (bound peptide/protein contact) was not computed project-wide; and tier C is a structural-segment proxy. The functional-site specificity report accordingly found residue-level mapped AF-2/H12 and LBP classes testable in a matched structural-context analysis, while the dimerization class lacked enough residue-level mappings; this is an evidence-coverage constraint, not a negative biological conclusion.

## Relationship to crystallographic B-factors

The literature reviewed here supplies a functional interpretation for selected positions or regions; it does **not** convert crystallographic B-factors into a solution-dynamics measurement. In this project, B values are normalized within coordinate instances and represent relative crystallographic mobility/disorder. They can reflect thermal motion, static/conformational disorder, occupancy/modeling and refinement effects, and crystal packing. A site’s relative B-factor or IQR across deposited structures is therefore an association with a structurally defined/curated region under heterogeneous crystallographic conditions—not direct evidence for molecular dynamics, ligand-induced kinetics, or altered receptor function in solution. This boundary is consistent with [`../FINAL_LIMITATIONS.md`](../FINAL_LIMITATIONS.md) and the v2 Role D endpoint definition in [`../rebuild_v2/README.md`](../rebuild_v2/README.md).

## Conclusions and priority gaps

1. The five anchors provide 44 manually curated canonical positions, but the evidence types differ: VDR/PPARγ/ERα/GR have primary-structure-supported residue assignments; NR4A1 deliberately has an atypical, limited definition; some charge-clamp/H12 entries rest on curated canonical-resource annotation rather than a direct peptide-contact experiment.
2. The strongest transferable conclusion is architectural, not residue-specific: many nuclear receptors have an LBD with ligand-responsive and AF-2-associated regions, but pocket chemistry, partner contacts and dimer interfaces are receptor-specific.
3. For project-wide functional-site results, report coactivator-groove and dimerization labels outside the anchor residues as **regional structural proxies**. Do not describe them as literature-curated contact residues.
4. Priority literature curation would be receptor-specific bound-peptide/protein contact evidence (tier B), experimentally supported dimer interfaces beyond ERα and GR, and explicit primary citations for any residue-level transfer proposed for non-anchor receptors. Such additions should enter a new curation file rather than modify the historical source table.

## References

Citation keys resolve to complete, machine-readable records in [`bibliography.csv`](bibliography.csv). Citations were checked against DOI landing metadata and/or PubMed on 2026-10-22.
