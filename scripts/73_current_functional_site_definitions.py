"""
Parts 1, 2, 4, 5, 6 of the current functional-site flexibility stage:
assemble data/analysis/current_functional_site_definitions.parquet, a NEW
analysis layer entirely separate from the historical Role D
reports/tables/functional_site_definitions.csv (read-only, never modified,
never overwritten).

Four site classes (Part 1): AF-2/H12, ligand-binding pocket (LBP),
coactivator-binding groove, dimerization interface. Not every receptor has
every site in an equivalent (or any) form -- every receptor x site
combination gets an explicit AVAILABLE / PARTIAL / AMBIGUOUS /
NOT_APPLICABLE / UNRESOLVED status (data/analysis/
current_functional_site_availability.csv); unresolved sites are never
given a placeholder position or a zero value.

Part 2 (AF-2/H12): taken directly from the CURRENT validated v3
structural-segment map (H12 rows of nr_lbd_structural_segment_map.csv) --
never the old heuristic H12 definitions.

Part 4 (coactivator groove) and Part 5 (dimerization interface) use an
evidence hierarchy:
  A. existing curated Role-D functional_site_definitions.csv residues
     (read-only; RESIDUE-level, HIGH confidence) -- covers only the
     handful of receptors Role D hand-curated.
  B. bound coactivator/coregulator peptide contacts, or observed
     receptor-receptor protein-protein contacts for dimerization --
     NOT implemented in this pass: this pipeline does not track
     polymer-partner chain identity outside the historical strict-QC
     subset (`final_polymer_partner_annotations.csv` covers only that
     1382-structure subset, per TEAM_HANDOFF.md), so this tier is
     reported as unavailable rather than approximated or invented.
  C. validated structural correspondence to the canonical H3/H4/H5
     (coactivator) or H9/H10/H11 (dimerization) segments from the v3
     structural-segment map -- SEGMENT-level, confidence inherited from
     the segment mapping's own assignment_confidence, explicitly lower
     than tier A.
Never automatically treats the entire H3/H4/H5 or H9/H10/H11 span as the
functional surface without saying so -- SEGMENT-level rows are always
labeled as such, never presented as exact residue-level annotations.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
sys.path.insert(0, str(PROJECT_ROOT))

POCKET_PRIMARY_THRESHOLD = 0.25
POCKET_SENSITIVITY_THRESHOLDS = [0.10, 0.50]

COACTIVATOR_SEGMENTS = ["H3", "H4", "H5"]
DIMERIZATION_SEGMENTS = ["H9", "H10", "H11"]

SITE_AF2_H12 = "AF2_H12"
SITE_LBP = "LIGAND_BINDING_POCKET"
SITE_COACT = "COACTIVATOR_GROOVE"
SITE_DIMER = "DIMERIZATION_INTERFACE"
ALL_SITES = [SITE_AF2_H12, SITE_LBP, SITE_COACT, SITE_DIMER]


def load_role_d_definitions() -> pd.DataFrame:
    """Read-only reference to the historical curated definitions -- never modified."""
    return pd.read_csv(TABLES / "functional_site_definitions.csv")


def receptor_list() -> pd.DataFrame:
    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    return rfs[rfs["n_coordinate_candidates"] > 0][["uniprot_id", "nr_code", "common_name"]].reset_index(drop=True)


def build_af2_h12(seg_map: pd.DataFrame) -> pd.DataFrame:
    h12 = seg_map[seg_map["structural_segment"] == "H12"].copy()
    rows = []
    for _, r in h12.iterrows():
        rows.append({
            "uniprot_id": r["uniprot_id"], "nr_code": r["nr_code"], "functional_site": SITE_AF2_H12,
            "canonical_uniprot_position": r["canonical_uniprot_position"], "structural_segment": "H12",
            "definition_level": "RESIDUE", "assignment_status": "ASSIGNED",
            "assignment_confidence": r["assignment_confidence"], "evidence_type": "current_structural_segment_map",
            "source_or_method": f"nr_lbd_structural_segment_map.csv ({r['assignment_method']})",
            "notes": "AF-2/H12 taken directly from the current v3 canonical H12 mapping, never the old heuristic definition.",
        })
    return pd.DataFrame(rows)


def build_lbp() -> pd.DataFrame:
    consensus = pd.read_parquet(ANALYSIS / "ligand_pocket_consensus_residues.parquet")
    primary = consensus[consensus["contact_fraction"] >= POCKET_PRIMARY_THRESHOLD]
    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")[["uniprot_id", "nr_code"]]
    primary = primary.merge(rfs, on="uniprot_id", how="left")
    rows = []
    for _, r in primary.iterrows():
        confidence = "HIGH" if r["contact_fraction"] >= 0.50 else "MODERATE"
        rows.append({
            "uniprot_id": r["uniprot_id"], "nr_code": r["nr_code"], "functional_site": SITE_LBP,
            "canonical_uniprot_position": r["canonical_uniprot_position"], "structural_segment": "",
            "definition_level": "RESIDUE", "assignment_status": "ASSIGNED",
            "assignment_confidence": confidence, "evidence_type": "geometric_ligand_contact_consensus",
            "source_or_method": (
                f"5.0A heavy-atom contact to the classified FUNCTIONAL_LBD_LIGAND, aggregated across "
                f"{int(r['n_holo_pdbs_observed'])} holo PDB entries; contact_fraction="
                f"{r['contact_fraction']:.2f} (>= {POCKET_PRIMARY_THRESHOLD:.0%} primary threshold)"
            ),
            "notes": "Receptor-level consensus pocket from real atomic geometry, not a universal residue list.",
        })
    return pd.DataFrame(rows)


def _tier_a_residues(role_d: pd.DataFrame, site_names: list[str], functional_site: str, note_suffix: str) -> pd.DataFrame:
    sub = role_d[role_d["site_name"].isin(site_names)]
    rows = []
    for _, r in sub.iterrows():
        rows.append({
            "uniprot_id": r["uniprot_id"], "nr_code": r["nr_code"], "functional_site": functional_site,
            "canonical_uniprot_position": r["uniprot_position"], "structural_segment": "",
            "definition_level": "RESIDUE", "assignment_status": "ASSIGNED", "assignment_confidence": "HIGH",
            "evidence_type": "role_d_curated_definition",
            "source_or_method": f"reports/tables/functional_site_definitions.csv (site '{r['site_name']}', {r['source']}), read-only",
            "notes": f"Tier A (curated residue-level evidence). {note_suffix} {r['notes']}",
        })
    return pd.DataFrame(rows)


def _tier_c_segment(seg_map: pd.DataFrame, segments: list[str], functional_site: str, tier_a_uids: set, note: str) -> pd.DataFrame:
    """Segment-level fallback for receptors without tier-A residue evidence."""
    sub = seg_map[seg_map["structural_segment"].isin(segments) & (~seg_map["uniprot_id"].isin(tier_a_uids))]
    rows = []
    for _, r in sub.iterrows():
        rows.append({
            "uniprot_id": r["uniprot_id"], "nr_code": r["nr_code"], "functional_site": functional_site,
            "canonical_uniprot_position": r["canonical_uniprot_position"], "structural_segment": r["structural_segment"],
            "definition_level": "SEGMENT", "assignment_status": "ASSIGNED",
            "assignment_confidence": r["assignment_confidence"], "evidence_type": "structural_segment_correspondence",
            "source_or_method": f"nr_lbd_structural_segment_map.csv, segment {r['structural_segment']} ({r['assignment_method']})",
            "notes": (
                f"Tier C (segment-level structural correspondence, NOT an exact residue-level functional "
                f"annotation -- the entire {r['structural_segment']} segment is not assumed to be "
                f"functionally equivalent at every position). {note}"
            ),
        })
    return pd.DataFrame(rows)


def build_coactivator_groove(role_d: pd.DataFrame, seg_map: pd.DataFrame) -> pd.DataFrame:
    tier_a = _tier_a_residues(role_d, ["AF-2 coactivator groove"], SITE_COACT,
                               "Charge-clamp / coactivator-contact residue.")
    tier_a_uids = set(tier_a["uniprot_id"].unique())
    tier_c = _tier_c_segment(seg_map, COACTIVATOR_SEGMENTS, SITE_COACT, tier_a_uids,
                              "Tier B (bound-coactivator-peptide contact) was not computed -- this pipeline "
                              "does not track polymer-partner chain identity outside the historical strict-QC "
                              "subset; reported as unavailable rather than approximated.")
    return pd.concat([tier_a, tier_c], ignore_index=True)


def build_dimerization_interface(role_d: pd.DataFrame, seg_map: pd.DataFrame) -> pd.DataFrame:
    tier_a = _tier_a_residues(role_d, ["LBD homodimerization", "LBD dimerization"], SITE_DIMER,
                               "Dimerization-interface residue.")
    tier_a_uids = set(tier_a["uniprot_id"].unique())
    tier_c = _tier_c_segment(seg_map, DIMERIZATION_SEGMENTS, SITE_DIMER, tier_a_uids,
                              "Tier B (observed receptor-receptor protein contacts) was not computed -- "
                              "same polymer-partner-identity limitation as the coactivator groove.")
    return pd.concat([tier_a, tier_c], ignore_index=True)


def build_availability_table(definitions: pd.DataFrame, receptors: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, rec in receptors.iterrows():
        uid, nr_code = rec["uniprot_id"], rec["nr_code"]
        for site in ALL_SITES:
            sub = definitions[(definitions["uniprot_id"] == uid) & (definitions["functional_site"] == site)]
            if len(sub) == 0:
                status, level, confidence = "UNRESOLVED", "", ""
            elif (sub["definition_level"] == "RESIDUE").any():
                status, level = "AVAILABLE", "RESIDUE"
                confidence = sub.loc[sub["definition_level"] == "RESIDUE", "assignment_confidence"].mode().iloc[0]
            else:
                status, level = "PARTIAL", "SEGMENT"
                confidence = sub["assignment_confidence"].mode().iloc[0] if len(sub) else ""
            rows.append({
                "uniprot_id": uid, "nr_code": nr_code, "functional_site": site,
                "status": status, "definition_level": level, "representative_confidence": confidence,
                "n_positions": len(sub),
            })
    return pd.DataFrame(rows)


def main():
    role_d = load_role_d_definitions()
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    receptors = receptor_list()

    af2 = build_af2_h12(seg_map)
    print(f"AF-2/H12: {af2['uniprot_id'].nunique()} receptors, {len(af2)} positions (from current v3 mapping)")

    lbp = build_lbp()
    print(f"LBP: {lbp['uniprot_id'].nunique()} receptors, {len(lbp)} positions "
          f"(consensus >= {POCKET_PRIMARY_THRESHOLD:.0%} contact fraction)")

    coact = build_coactivator_groove(role_d, seg_map)
    n_tier_a_coact = (coact["evidence_type"] == "role_d_curated_definition").sum()
    n_tier_c_coact = (coact["evidence_type"] == "structural_segment_correspondence").sum()
    print(f"Coactivator groove: {coact['uniprot_id'].nunique()} receptors "
          f"({n_tier_a_coact} tier-A residue rows, {n_tier_c_coact} tier-C segment rows)")

    dimer = build_dimerization_interface(role_d, seg_map)
    n_tier_a_dimer = (dimer["evidence_type"] == "role_d_curated_definition").sum()
    n_tier_c_dimer = (dimer["evidence_type"] == "structural_segment_correspondence").sum()
    print(f"Dimerization interface: {dimer['uniprot_id'].nunique()} receptors "
          f"({n_tier_a_dimer} tier-A residue rows, {n_tier_c_dimer} tier-C segment rows)")

    definitions = pd.concat([af2, lbp, coact, dimer], ignore_index=True)
    definitions = definitions.sort_values(["functional_site", "uniprot_id", "canonical_uniprot_position"]).reset_index(drop=True)

    out = ANALYSIS / "current_functional_site_definitions.parquet"
    definitions.to_parquet(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(definitions)} rows) -- "
          "a NEW analysis layer, separate from the historical Role D definitions file.")

    availability = build_availability_table(definitions, receptors)
    avail_out = TABLES / "current_functional_site_availability.csv"
    availability.to_csv(avail_out, index=False)
    print(f"Wrote {avail_out.relative_to(PROJECT_ROOT)} ({len(availability)} receptor x site rows)")
    print(availability.groupby(["functional_site", "status"]).size().unstack(fill_value=0).to_string())

    # sensitivity pocket thresholds recorded for later reference (Part 3)
    consensus = pd.read_parquet(ANALYSIS / "ligand_pocket_consensus_residues.parquet")
    for t in POCKET_SENSITIVITY_THRESHOLDS:
        n = (consensus["contact_fraction"] >= t).sum()
        n_rec = consensus.loc[consensus["contact_fraction"] >= t, "uniprot_id"].nunique()
        print(f"Pocket sensitivity check: >= {t:.0%} -> {n} residues across {n_rec} receptors")


if __name__ == "__main__":
    main()
