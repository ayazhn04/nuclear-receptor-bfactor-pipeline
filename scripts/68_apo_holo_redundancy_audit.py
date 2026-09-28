"""
Part 12 of the APO vs HOLO flexibility stage: audit whether one receptor's
APO or HOLO group is dominated by a single repeated ligand or an
apparently identical construct series, using only metadata this pipeline
already tracks (data/manifests/final_ligand_annotations.csv for ligand
identity/names; structures.csv construct-boundary fields as a weak
construct-similarity proxy). Publication/deposition-series and mutation
metadata are NOT tracked anywhere in this pipeline, so that specific check
is reported as "not assessable with current metadata" rather than
fabricated. A PDB deposition is never assumed to be a fully independent
biological experiment -- dominance by a single ligand or construct is
flagged explicitly wherever it is found.
"""
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
TABLES = PROJECT_ROOT / "reports" / "tables"

DOMINANCE_THRESHOLD = 0.5  # a single ligand/construct value accounting for >50% of a
                            # group is flagged as dominant


def ligand_dominance(la: pd.DataFrame, uid: str, pdb_ids: set) -> dict:
    sub = la[(la["uniprot_id"] == uid) & (la["pdb_id"].isin(pdb_ids))]
    names = []
    for raw in sub["functional_ligand_names_json"]:
        try:
            parsed = json.loads(raw) if raw else []
        except (json.JSONDecodeError, TypeError):
            parsed = []
        if parsed:
            names.append(parsed[0])  # first-listed functional ligand name for this entry
        else:
            names.append("(none/unresolved)")
    if not names:
        return {"n_entries": 0, "n_distinct_ligands": 0, "top_ligand": None, "top_ligand_fraction": None}
    counts = pd.Series(names).value_counts()
    return {
        "n_entries": len(names), "n_distinct_ligands": len(counts),
        "top_ligand": counts.index[0], "top_ligand_fraction": float(counts.iloc[0] / len(names)),
    }


def main():
    eligibility = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    eligible = eligibility[eligibility["eligible_exploratory_3_3"]]

    structures = pd.read_csv(PROCESSED / "structures.csv", dtype={"pdb_id": str}, keep_default_na=False)
    coord = structures[structures["coordinate_data_status"] == "COORDINATE_DATA_AVAILABLE"]
    coord_pdbs = coord.groupby(["uniprot_id", "pdb_id"]).agg(ligand_state=("ligand_state", "first")).reset_index()

    la = pd.read_csv(MANIFESTS / "final_ligand_annotations.csv", keep_default_na=False)

    lines = []
    lines.append("# APO/HOLO Redundancy & Dominance Audit\n")
    lines.append(
        "> Checks whether a receptor's APO or HOLO group is dominated by a single repeated "
        "ligand. Publication/deposition-series and mutation metadata are **not tracked** "
        "anywhere in this pipeline -- that specific check is reported as not assessable "
        "with current metadata rather than fabricated. No PDB deposition is assumed to be a "
        f"fully independent biological experiment. Dominance threshold: a single ligand "
        f"accounting for >{DOMINANCE_THRESHOLD:.0%} of a group's entries is flagged.\n"
    )

    any_flagged = False
    for _, r in eligible.iterrows():
        uid, nr_code = r["uniprot_id"], r["nr_code"]
        g = coord_pdbs[coord_pdbs["uniprot_id"] == uid]
        apo_pdbs = set(g.loc[g["ligand_state"] == "APO", "pdb_id"])
        holo_pdbs = set(g.loc[g["ligand_state"] == "HOLO", "pdb_id"])

        apo_dom = ligand_dominance(la, uid, apo_pdbs)
        holo_dom = ligand_dominance(la, uid, holo_pdbs)

        lines.append(f"## {nr_code} ({uid})\n")
        lines.append(f"- APO: {apo_dom['n_entries']} entries, {apo_dom['n_distinct_ligands']} distinct "
                     f"ligand identities (APO structures have no functional ligand by definition; "
                     f"'{apo_dom['top_ligand']}' reflects crystallization additives/incidental "
                     f"components if any, not the APO/HOLO classification itself).")
        holo_flag = holo_dom["top_ligand_fraction"] is not None and holo_dom["top_ligand_fraction"] > DOMINANCE_THRESHOLD
        lines.append(
            f"- HOLO: {holo_dom['n_entries']} entries, {holo_dom['n_distinct_ligands']} distinct ligands; "
            f"most common = **{holo_dom['top_ligand']}** "
            f"({holo_dom['top_ligand_fraction']:.0%} of HOLO entries)"
            + (" -- **FLAGGED: HOLO group dominated by a single ligand.**" if holo_flag else "")
        )
        if holo_flag:
            any_flagged = True
        lines.append(
            "- Publication/deposition-series dominance: **not assessable** -- this pipeline does "
            "not track citation/deposition-group metadata.\n"
        )

    out = TABLES / "apo_holo_redundancy_audit.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out.relative_to(PROJECT_ROOT)}")
    print(f"\n{'Some' if any_flagged else 'No'} receptor(s) flagged for single-ligand HOLO dominance "
          f"(threshold >{DOMINANCE_THRESHOLD:.0%}).")


if __name__ == "__main__":
    main()
