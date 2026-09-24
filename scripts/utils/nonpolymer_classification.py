"""
Final Data Release, Step 8: nonpolymer relevance classification.

Multi-evidence, documented cascade — never a single hand-picked list
assumed complete. Evidence combined, in priority order:

1. Monatomic components (atom_count_heavy == 1) are always ION. NR LBDs are
   never themselves metal-ion ligands (unlike e.g. zinc-finger DBDs), so a
   monatomic species near the LBD is treated as buffer/crystallization
   chemistry regardless of contact, per the explicit project instruction
   that simple ions must never by themselves create HOLO status.
2. A documented, explicitly-listed set of common crystallization
   additives/cryoprotectants/buffers/detergents/heavy-atom phasing reagents
   (CRYSTALLIZATION_REFERENCE_SET below, each with a one-line identity
   note) -> CRYSTALLIZATION_ADDITIVE_OR_BUFFER, UNLESS this specific
   entity instance also carries RCSB's own SUBJECT_OF_INVESTIGATION
   annotation (observed, in this dataset, on some instances of otherwise
   generic compounds like sulfate/glycerol/acetate — RCSB itself is
   telling us that particular occurrence was the analyte of interest, not
   incidental buffer chemistry) -> POSSIBLE_FUNCTIONAL_LIGAND instead.
3. A deposited covalent linkage to the polymer (`_struct_conn`,
   conn_type_id contains "coval") -> COVALENT_OR_MODIFIED_COMPONENT.
4. RCSB's own SUBJECT_OF_INVESTIGATION annotation AND geometric LBD
   contact (<=4.5 A) -> FUNCTIONAL_LBD_LIGAND.
5. A deposited `_struct_site` binding-site record naming this exact
   component AND LBD contact AND a drug/hormone-scale heavy-atom count
   (>=5) -> FUNCTIONAL_LBD_LIGAND.
6. LBD contact AND heavy-atom count >=10 (well above typical
   solvent/buffer/ion size; drug- and hormone-like ligands in this
   dataset are essentially all >=14 heavy atoms) and not on the
   crystallization reference set -> FUNCTIONAL_LBD_LIGAND.
7. LBD contact with 5-9 heavy atoms (ambiguous size band) and not on the
   reference set -> POSSIBLE_FUNCTIONAL_LIGAND.
8. LBD contact with <5 heavy atoms and not ion/additive ->
   POSSIBLE_FUNCTIONAL_LIGAND (small molecule genuinely in contact with
   the pocket; too small to confidently call a drug-like ligand, so
   flagged for review rather than assumed).
9. No LBD contact: on the reference set -> CRYSTALLIZATION_ADDITIVE_OR_BUFFER;
   heavy-atom count <=2 -> SOLVENT; otherwise -> OTHER_NONPOLYMER.
10. Anything not resolved by the above -> REVIEW_NEEDED (should not occur
    given 1-9 are exhaustive, but kept as an explicit fallback rather than
    a silent default).

Water (HOH/DOD) never reaches this function — it is enumerated and
classified WATER directly, without geometry, because per the project's
apo/holo policy water can never by itself indicate HOLO regardless of
contact distance.
"""

from __future__ import annotations

CONTACT_THRESHOLD_A = 4.5

# Common crystallization additives, cryoprotectants, buffers, detergents,
# and heavy-atom phasing reagents identified by reviewing every nonpolymer
# component that recurs >=2 times across this project's 1900 primary-set
# nonpolymer entities (see data/manifests/final_nonpolymer_inventory.csv).
# This list is used only as ONE input among several (see module docstring)
# — it never silently overrides RCSB's own per-instance
# SUBJECT_OF_INVESTIGATION annotation, and a component's ABSENCE from this
# list never by itself implies functional relevance (that still requires
# LBD contact and/or another positive signal).
CRYSTALLIZATION_REFERENCE_SET: dict[str, str] = {
    "GOL": "glycerol — cryoprotectant",
    "EDO": "1,2-ethanediol (ethylene glycol) — cryoprotectant",
    "PGO": "S-1,2-propanediol — cryoprotectant",
    "PGR": "R-1,2-propanediol — cryoprotectant",
    "MPD": "(4S)-2-methyl-2,4-pentanediol — cryoprotectant/precipitant",
    "MRD": "(4R)-2-methyl-2,4-pentanediol — cryoprotectant/precipitant",
    "IPA": "isopropyl alcohol — cryoprotectant/solvent",
    "DMS": "dimethyl sulfoxide — cosolvent",
    "PEG": "di(hydroxyethyl)ether (short PEG) — precipitant/cryoprotectant",
    "PG4": "tetraethylene glycol — precipitant/cryoprotectant",
    "PGE": "triethylene glycol — precipitant/cryoprotectant",
    "PE4": "octaethylene-glycol-like PEG chain — precipitant/cryoprotectant",
    "1PE": "pentaethylene glycol — precipitant/cryoprotectant",
    "P6G": "hexaethylene glycol — precipitant/cryoprotectant",
    "15P": "PEG-family chain — precipitant/cryoprotectant",
    "BU1": "1,4-butanediol — cryoprotectant",
    "SO4": "sulfate ion — crystallization precipitant",
    "PO4": "phosphate ion — buffer/precipitant",
    "ACT": "acetate ion — buffer",
    "FMT": "formic acid/formate — buffer/precipitant",
    "FLC": "citrate anion — buffer",
    "MLA": "malonic acid — buffer",
    "MLI": "malonate ion — buffer",
    "CAC": "cacodylate ion — buffer",
    "TRS": "Tris buffer (2-amino-2-hydroxymethyl-propane-1,3-diol)",
    "EPE": "HEPES buffer",
    "NHE": "CHES-like buffer (N-cyclohexylaminoethanesulfonic acid)",
    "IMD": "imidazole — buffer/additive",
    "BME": "beta-mercaptoethanol — reducing agent",
    "DTT": "dithiothreitol — reducing agent",
    "B7G": "heptyl beta-D-glucopyranoside — detergent",
    "BOG": "octyl beta-D-glucopyranoside — detergent",
    "JZR": "hexyl beta-D-glucopyranoside — detergent",
    "CPS": "CHAPS-like zwitterionic detergent",
    "TBY": "tributyl(chloro)stannane — heavy-atom phasing reagent",
    "T9T": "chlorotri(phenyl)stannane — heavy-atom phasing reagent",
}


def classify_nonpolymer_instance(
    comp_id: str,
    atom_count_heavy: int | None,
    contacts_lbd: bool,
    min_distance_to_lbd: float | None,
    has_subject_of_investigation_annotation: bool,
    has_struct_site_record: bool,
    is_covalently_linked: bool,
) -> tuple[str, str]:
    """Returns (classification, evidence_summary)."""
    heavy = atom_count_heavy if atom_count_heavy is not None else -1
    evidence_bits = [
        f"heavy_atoms={heavy}",
        f"contacts_lbd_4.5A={contacts_lbd}",
        f"min_dist={min_distance_to_lbd:.2f}A" if min_distance_to_lbd is not None else "min_dist=NA",
        f"rcsb_subject_of_investigation={has_subject_of_investigation_annotation}",
        f"struct_site_record={has_struct_site_record}",
        f"covalently_linked={is_covalently_linked}",
    ]
    evidence = "; ".join(evidence_bits)

    if heavy == 1:
        return "ION", f"Monatomic component (atom_count_heavy==1); {evidence}"

    on_reference_set = comp_id in CRYSTALLIZATION_REFERENCE_SET
    if on_reference_set and not has_subject_of_investigation_annotation:
        return (
            "CRYSTALLIZATION_ADDITIVE_OR_BUFFER",
            f"{CRYSTALLIZATION_REFERENCE_SET[comp_id]}; {evidence}",
        )

    if is_covalently_linked:
        return "COVALENT_OR_MODIFIED_COMPONENT", f"Deposited covalent linkage (_struct_conn); {evidence}"

    if has_subject_of_investigation_annotation and contacts_lbd:
        return "FUNCTIONAL_LBD_LIGAND", f"RCSB SUBJECT_OF_INVESTIGATION annotation and LBD contact; {evidence}"

    if has_struct_site_record and contacts_lbd and heavy >= 5:
        return "FUNCTIONAL_LBD_LIGAND", f"Deposited _struct_site binding-site record and LBD contact; {evidence}"

    if contacts_lbd and heavy >= 10 and not on_reference_set:
        return "FUNCTIONAL_LBD_LIGAND", f"LBD contact with drug/hormone-scale heavy-atom count; {evidence}"

    if contacts_lbd and 5 <= heavy < 10 and not on_reference_set:
        return "POSSIBLE_FUNCTIONAL_LIGAND", f"LBD contact, mid-size (ambiguous) heavy-atom count; {evidence}"

    if contacts_lbd and heavy < 5 and not on_reference_set:
        return "POSSIBLE_FUNCTIONAL_LIGAND", f"LBD contact, small heavy-atom count, not a recognized additive; {evidence}"

    if not contacts_lbd:
        if on_reference_set:
            return "CRYSTALLIZATION_ADDITIVE_OR_BUFFER", f"{CRYSTALLIZATION_REFERENCE_SET[comp_id]}; no LBD contact; {evidence}"
        if heavy <= 2:
            return "SOLVENT", f"No LBD contact, minimal heavy-atom count; {evidence}"
        return "OTHER_NONPOLYMER", f"No LBD contact, not a recognized additive/ion; {evidence}"

    return "REVIEW_NEEDED", f"Did not match any classification rule; {evidence}"
