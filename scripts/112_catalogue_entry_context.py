"""
Receptor catalogue (v2), step 1: one row per receptor x PDB entry with every context
dimension the per-receptor pages show. Catalogue output only: nothing in the analysis
pipeline is read back from here, and no v2 analysis table is changed.

Per entry:
  inventory    discovered / coordinate data / contributes z-scores / TRACE_LBD-only /
               primary-QC member / reason when there are no coordinate data
  ligand       analysis ligand state (APO/HOLO/AMBIGUOUS, TRACE_LBD-only kept apart) and
               the functional ligand component(s)
  partners     every other polymer chain within 5 A of the receptor's LBD atoms, classified
               (see below); "standalone" = no contacting partner chain
  construct    first/last UniProt residue of the entity's sequence-evidence alignment (the whole
               construct, not clipped to the analysis span), number of aligned segments
               (>1 = internal deletion), fraction of the analysis span observed in the
               representative chain, fusion/chimera flag, point substitutions against the
               reference sequence (standard residues only; modified residues are not counted)
  experiment   method, resolution, R-free; from the entry's own mmCIF: space group, unit cell,
               deposition date, polymer chains and receptor copies in the asymmetric unit,
               and a crystal-form id (same space group, cell within 2% / 2 degrees)
  sequence     every observed residue of the representative chain that differs from the
               reference: engineered / natural variant / sequence conflict / other declared
               (depositor's _struct_ref_seq_dif annotation or a UniProt natural variant),
               UNDECLARED when there is no annotation (no claim about its origin), and
               MODIFIED_RESIDUE for non-standard components on a reference residue

Partner context extends scripts/27 (which annotates only the 1,395 primary-QC entries)
to all coordinate entries, with the same contact rule (any heavy atom of the partner
chain within 5.0 A of an LBD heavy atom with positive occupancy) and the same
classifier (scripts/27.classify_partner). Two catalogue differences, both recorded:
  * entity metadata (description, length, UniProt) is read from the local mmCIF
    (_entity, _entity_poly_seq, _struct_ref) instead of the RCSB API cache;
  * another project receptor's entity in the same entry counts as a partner
    (scripts/27 skips every project-receptor entity, so e.g. PPARg-RXRa heterodimers
    show no partner there). Other chains of the receptor's own entity (homodimer copies)
    are not partners, as in scripts/27.
Catalogue classes (presentation grouping of scripts/27 categories + description):
  COREGULATOR_PEPTIDE   scripts/27 COREGULATOR_CANDIDATE (keyword evidence)
  NR_PARTNER            another nuclear receptor (project-receptor UniProt, or a
                        receptor-name description on a chain > 40 aa)
  SAME_RECEPTOR_CHAIN   a different entity of the same receptor (e.g. a receptor peptide)
  OTHER_PEPTIDE         <= 40 aa without coregulator evidence (inhibitors, unnamed peptides)
  NUCLEIC_ACID          DNA/RNA polymer
  OTHER_PROTEIN         anything else

Reads (all committed v2 outputs, plus the local mmCIF caches listed in
data/rebuild_v2/inputs/raw_cache_requirements_v2.csv):
  data/rebuild_v2/processed/{structures.csv, lbd_residue_bfactors_all.parquet}
  data/rebuild_v2/inputs/{rcsb_sequence_coordinate_alignments_v2.csv, analysis_span_v2.csv}
  data/rebuild_v2/analysis/pdb_position_bfactor_profiles.parquet
  data/rebuild_v2/manifests/{final_ligand_annotations.csv, final_polymer_partner_annotations.csv,
                             final_structure_selection_audit.csv}
Writes:
  reports/catalogue_v2/tables/catalogue_entry_context_v2.csv
  reports/catalogue_v2/tables/catalogue_partner_contacts_v2.csv
  reports/catalogue_v2/tables/catalogue_partner_validation_vs_script27.csv
  reports/catalogue_v2/tables/catalogue_sequence_differences_v2.csv  (audited class + evidence per residue)
  reports/catalogue_v2/tables/catalogue_mutation_audit_v2.csv        (one row per usable entry)
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402
from scripts.utils import partner_identity  # noqa: E402
from scripts.utils.mmcif_parse import as_list, get_entity_types, get_heavy_atoms, parse_mmcif_gz  # noqa: E402

_spec = importlib.util.spec_from_file_location("s27", PROJECT_ROOT / "scripts" / "27_polymer_partner_annotation.py")
s27 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s27)

V2 = PROJECT_ROOT / "data" / "rebuild_v2"
OUT = PROJECT_ROOT / "reports" / "catalogue_v2" / "tables"
MMCIF_CACHES = [PROJECT_ROOT / "data" / "raw" / "mmcif", PROJECT_ROOT / "data" / "interim" / "raw_mmcif_v2_additions"]
CONTACT_A = s27.CONTACT_THRESHOLD_A
PEPTIDE_MAX = s27.COREGULATOR_MAX_LENGTH
STANDARD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
               "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
PROJECT_UNIPROT = {r["uniprot_id"] for r in NR_METADATA}
NR_NAME = re.compile(r"(nuclear receptor subfamily|retinoi[cd] (acid )?x? ?receptor|retinoid x receptor|\brxr|"
                     r"estrogen receptor|estrogen-related receptor|androgen receptor|glucocorticoid receptor|"
                     r"mineralocorticoid receptor|progesterone receptor|thyroid hormone receptor|vitamin d3? receptor|"
                     r"peroxisome proliferator[- ]activated receptor|liver x receptor|oxysterols receptor|"
                     r"farnesoid|bile acid receptor|pregnane x receptor|constitutive androstane|rev-?erb|"
                     r"ror[- ]?(alpha|beta|gamma)|hepatocyte nuclear factor 4|nr4a|nurr1|nur77|steroidogenic factor|"
                     r"coup-tf|orphan nuclear receptor)", re.I)


def mmcif_path(pdb: str) -> Path:
    for c in MMCIF_CACHES:
        if (c / f"{pdb}.cif.gz").exists():
            return c / f"{pdb}.cif.gz"
    raise FileNotFoundError(f"{pdb}.cif.gz not in the local mmCIF caches (see scripts/108)")


def entity_meta(d: dict) -> dict[str, dict]:
    desc = dict(zip(as_list(d.get("_entity.id", [])), as_list(d.get("_entity.pdbx_description", []))))
    ptype = dict(zip(as_list(d.get("_entity_poly.entity_id", [])), as_list(d.get("_entity_poly.type", []))))
    length: dict[str, set] = {}
    for e, n in zip(as_list(d.get("_entity_poly_seq.entity_id", [])), as_list(d.get("_entity_poly_seq.num", []))):
        length.setdefault(e, set()).add(n)
    unp: dict[str, list] = {}
    if "_struct_ref.entity_id" in d:
        for e, db, acc in zip(as_list(d["_struct_ref.entity_id"]), as_list(d["_struct_ref.db_name"]),
                              as_list(d["_struct_ref.pdbx_db_accession"])):
            if db == "UNP" and acc not in unp.get(e, []):
                unp.setdefault(e, []).append(acc)
    return {e: {"description": "" if desc.get(e) in (None, "?", ".") else desc[e], "poly_type": ptype.get(e, ""),
                "length": len(length.get(e, ())), "uniprot": unp.get(e, [])} for e in desc}


def catalogue_class(meta: dict, own_uid: str, partner_uid: str | None, official: tuple | None = None) -> tuple[str, str, str]:
    """(catalogue class, scripts/27 category, evidence). `official` = the (category, evidence)
    scripts/27 recorded for this chain, used whenever it exists (primary-QC entries)."""
    if "polypeptide" not in meta["poly_type"]:
        return "NUCLEIC_ACID", official[0] if official else "", f"polymer type {meta['poly_type']}"
    cat, ev = official or s27.classify_partner(meta["description"], meta["length"] or None, meta["uniprot"], own_uid)
    if partner_uid == own_uid or (cat == "OTHER_RECEPTOR_PARTNER" and meta["uniprot"] == [own_uid]):
        return "SAME_RECEPTOR_CHAIN", cat, ev
    if partner_uid and partner_uid in PROJECT_UNIPROT:
        return "NR_PARTNER", cat, f"project receptor {partner_uid} in the same entry"
    if cat == "OTHER_RECEPTOR_PARTNER":
        return "NR_PARTNER", cat, ev
    if cat == "COREGULATOR_CANDIDATE":
        return "COREGULATOR_PEPTIDE", cat, ev
    if meta["length"] and meta["length"] > PEPTIDE_MAX and NR_NAME.search(meta["description"] or ""):
        return "NR_PARTNER", cat, f"nuclear-receptor name in description: '{meta['description']}'"
    if cat == "UNCLASSIFIED_PARTNER":
        return "OTHER_PEPTIDE", cat, ev
    return "OTHER_PROTEIN", cat, ev


def first(d: dict, key: str) -> str:
    v = as_list(d.get(key, []))
    return "" if not v or v[0] in ("?", ".") else v[0]


def crystal_meta(d: dict, pdb: str) -> dict:
    """Experimental metadata recorded in the entry's own mmCIF (no network)."""
    poly = set(as_list(d.get("_pdbx_poly_seq_scheme.asym_id", [])))
    ent_of_asym = dict(zip(as_list(d.get("_struct_asym.id", [])), as_list(d.get("_struct_asym.entity_id", []))))
    copies: dict[str, int] = {}
    for a in poly:
        copies[ent_of_asym.get(a, "")] = copies.get(ent_of_asym.get(a, ""), 0) + 1
    num = lambda k: pd.to_numeric(first(d, k), errors="coerce")  # noqa: E731
    seq_dif = []
    if "_struct_ref_seq_dif.seq_num" in d:
        k = len(as_list(d["_struct_ref_seq_dif.seq_num"]))
        col = lambda c: as_list(d.get(f"_struct_ref_seq_dif.{c}", ["?"] * k))  # noqa: E731
        for row in zip(col("pdbx_pdb_strand_id"), col("seq_num"), col("mon_id"), col("db_mon_id"), col("details"),
                       col("pdbx_seq_db_seq_num"), col("pdbx_seq_db_accession_code")):
            seq_dif.append(row)
    entity_mutation = {e: m for e, m in zip(as_list(d.get("_entity.id", [])), as_list(d.get("_entity.pdbx_mutation", [])))
                       if m not in ("?", ".", "")}
    return {"pdb_id": pdb,
            "space_group": first(d, "_symmetry.space_group_name_H-M") or first(d, "_space_group.name_H-M_alt"),
            "cell_a": num("_cell.length_a"), "cell_b": num("_cell.length_b"), "cell_c": num("_cell.length_c"),
            "cell_alpha": num("_cell.angle_alpha"), "cell_beta": num("_cell.angle_beta"), "cell_gamma": num("_cell.angle_gamma"),
            "deposition_date": first(d, "_pdbx_database_status.recvd_initial_deposition_date"),
            "n_polymer_chains_asu": len(poly), "copies_by_entity": copies, "seq_dif": seq_dif,
            "entity_mutation": entity_mutation}


def partners_for_pdb(job: tuple) -> dict:
    pdb, cands, lbd_by_instance, uid_by_entity, official = job
    d = parse_mmcif_gz(mmcif_path(pdb))
    cmeta = crystal_meta(d, pdb)
    etypes = get_entity_types(d)
    meta = entity_meta(d)
    atoms_by_asym: dict[str, list] = {}
    for a in get_heavy_atoms(d):
        if (a["occupancy"] or 0) > 0:
            atoms_by_asym.setdefault(a["label_asym_id"], []).append(a)
    rows = []
    for c in cands:
        lbd = lbd_by_instance.get(c["rep"])
        if not lbd:
            continue
        seq_ids, pos_by_seq = lbd
        lbd_atoms = [a for a in atoms_by_asym.get(c["asym"], []) if a["label_seq_id"] in seq_ids]
        if not lbd_atoms:
            continue
        lc = np.array([[a["x"], a["y"], a["z"]] for a in lbd_atoms])
        lp = np.array([pos_by_seq.get(a["label_seq_id"], -1) for a in lbd_atoms])
        for asym, atoms in atoms_by_asym.items():
            ent = atoms[0]["label_entity_id"]
            if ent == c["entity"] or etypes.get(ent) != "polymer":
                continue
            dm = cdist(np.array([[a["x"], a["y"], a["z"]] for a in atoms]), lc)
            if dm.min() > CONTACT_A:
                continue
            m = meta.get(ent, {"description": "", "poly_type": "", "length": 0, "uniprot": []})
            cls, cat27, ev = catalogue_class(m, c["uid"], uid_by_entity.get(ent),
                                             official.get((c["uid"], f"{pdb}_{ent}", asym)))
            rows.append({"uniprot_id": c["uid"], "pdb_id": pdb, "receptor_instance_id": c["rep"],
                         "partner_entity_id": f"{pdb}_{ent}", "partner_label_asym_id": asym,
                         "partner_description": m["description"], "partner_polymer_type": m["poly_type"],
                         "partner_length": m["length"], "partner_uniprot_json": json.dumps(m["uniprot"]),
                         "partner_project_receptor": uid_by_entity.get(ent, ""),
                         "min_distance_A": round(float(dm.min()), 3),
                         "n_contacted_lbd_positions": int(len({int(p) for p in lp[(dm <= CONTACT_A).any(axis=0)] if p >= 0})),
                         "catalogue_partner_class": cls, "script27_rule_category": cat27, "evidence": ev})
    return {"contacts": rows, "meta": cmeta}


DIFF_CLASSES = ["ENGINEERED", "NATURAL_VARIANT_OR_ISOFORM", "MODIFIED_RESIDUE", "MAPPING_REFERENCE_DIFFERENCE",
                "UNCERTAIN"]
MUTATION_SOURCE = ("per residue: depositor _struct_ref_seq_dif ('engineered mutation') consistent with the reference "
                   "residue, and/or the entity's _entity.pdbx_mutation list (same residues and position)")
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
                "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
                "TYR": "Y", "VAL": "V"}


def uniprot_features(uid: str) -> tuple[dict[int, set], dict[int, set]]:
    """UniProt single-residue 'Natural variant' and 'Alternative sequence' (isoform) features:
    position -> alternative residues (one-letter)."""
    path = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot" / f"{uid}.json"
    nv: dict[int, set] = {}
    iso: dict[int, set] = {}
    if not path.exists():
        return nv, iso
    for f in json.loads(path.read_text()).get("features", []):
        loc = f["location"]
        if loc["start"]["value"] != loc["end"]["value"]:
            continue
        alts = (f.get("alternativeSequence") or {}).get("alternativeSequences") or []
        target = nv if f.get("type") == "Natural variant" else iso if f.get("type") == "Alternative sequence" else None
        if target is not None:
            target.setdefault(int(loc["start"]["value"]), set()).update(a for a in alts if len(a) == 1)
    return nv, iso


MUT_TOKEN = re.compile(r"([A-Z])(\d+)([A-Z])")


def classify_substitutions(diff: pd.DataFrame, s: pd.DataFrame, cmeta: dict) -> pd.DataFrame:
    """Audited class of every observed residue that differs from the reference sequence in a
    representative chain (corrected v2 numbering, so position mismatches are already removed).

    ENGINEERED                    depositor marks the residue 'engineered mutation' in
                                  _struct_ref_seq_dif AND its reference residue (db_mon_id) is the
                                  same as ours; or the entity's _entity.pdbx_mutation lists the same
                                  reference/observed residues at our or the depositor's position
    NATURAL_VARIANT_OR_ISOFORM    depositor 'variant'/'polymorphism', or a UniProt natural variant /
                                  single-residue isoform difference at the position with that residue
    MODIFIED_RESIDUE              non-standard component (MSE, CSO, ...) on a reference residue
    MAPPING_REFERENCE_DIFFERENCE  depositor says 'engineered' but its reference residue differs from
                                  ours (a different reference sequence or numbering): not counted
    UNCERTAIN                     'conflict', other annotations, or no annotation at all
    Only ENGINEERED counts toward engineered-mutant fractions."""
    rep = s[s["representative_instance_id"] != ""]
    auth = dict(zip(rep["representative_instance_id"], rep["auth_asym_id"]))
    ent = dict(zip(rep["representative_instance_id"], rep["entity_id"].astype(str)))
    feat_cache: dict[str, tuple] = {}
    rows = []
    for r in diff.itertuples():
        pos = int(r.canonical_uniprot_position)
        obs_aa, ref_aa = r.observed_comp, r.canonical_residue_name
        ref1, obs1 = THREE_TO_ONE.get(ref_aa, "X"), THREE_TO_ONE.get(obs_aa, "X")
        meta = cmeta.get(r.pdb_id, {})
        det = db_mon = db_num = ""
        for strand, seq, mon, db, d, dbn, _acc in meta.get("seq_dif", []):
            if strand == auth.get(r.instance_id) and seq == str(int(r.entity_label_seq_id)) and mon == obs_aa:
                det, db_mon, db_num = (x if x not in ("?", ".") else "" for x in (d, db, dbn))
                break
        ent_list = meta.get("entity_mutation", {}).get(ent.get(r.instance_id, ""), "")
        ent_hit = any(a == ref1 and c == obs1 and int(n) in {pos, int(db_num) if db_num.isdigit() else -1}
                      for a, n, c in MUT_TOKEN.findall(ent_list.upper()))
        dl = det.lower()
        sd_eng = "engineer" in dl or "mutation" in dl or "mutagen" in dl
        sd_consistent = db_mon == ref_aa
        if obs_aa not in STANDARD_AA:
            cls, detail = "MODIFIED_RESIDUE", "non-standard component"
        elif (sd_eng and sd_consistent) or ent_hit:
            cls = "ENGINEERED"
            detail = "seq_dif+entity" if (sd_eng and sd_consistent and ent_hit) else "seq_dif" if sd_eng else "entity_mutation"
        elif sd_eng:
            cls, detail = "MAPPING_REFERENCE_DIFFERENCE", f"engineered vs depositor reference {db_mon or '?'}"
        else:
            nv, iso = feat_cache.setdefault(r.uniprot_id, uniprot_features(r.uniprot_id))
            ok_ref = not (r.uniprot_id == "P10827" and pos > 370)  # TRalpha1 numbering beyond the shared region
            if "variant" in dl or "polymorph" in dl:
                cls, detail = "NATURAL_VARIANT_OR_ISOFORM", "depositor variant"
            elif ok_ref and obs1 in nv.get(pos, set()):
                cls, detail = "NATURAL_VARIANT_OR_ISOFORM", "UniProt natural variant"
            elif ok_ref and obs1 in iso.get(pos, set()):
                cls, detail = "NATURAL_VARIANT_OR_ISOFORM", "UniProt isoform"
            else:
                cls = "UNCERTAIN"
                detail = "depositor conflict" if "conflict" in dl else f"depositor: {det}" if det else "no annotation"
        rows.append({"uniprot_id": r.uniprot_id, "pdb_id": r.pdb_id, "instance_id": r.instance_id, "position": pos,
                     "reference_residue": ref_aa, "observed_residue": obs_aa, "difference_class": cls,
                     "evidence": detail, "depositor_annotation": det, "depositor_reference_residue": db_mon,
                     "depositor_position": db_num, "entity_mutation_list": ent_list,
                     "label": f"{ref1}{pos}{obs1}" if cls != "MODIFIED_RESIDUE" else f"{ref_aa}{pos}{obs_aa}"})
    return pd.DataFrame(rows, columns=["uniprot_id", "pdb_id", "instance_id", "position", "reference_residue",
                                       "observed_residue", "difference_class", "evidence", "depositor_annotation",
                                       "depositor_reference_residue", "depositor_position", "entity_mutation_list",
                                       "label"])


def mutation_audit(subst: pd.DataFrame, entries: pd.DataFrame, cmeta: dict, s: pd.DataFrame) -> pd.DataFrame:
    """One row per usable entry. An entity mutation list naming a residue outside the observed
    chain cannot be checked residue by residue and is reported, not counted."""
    rep = s[s["representative_instance_id"] != ""]
    ent = rep.groupby(["uniprot_id", "pdb_id"])["entity_id"].agg(lambda e: sorted(set(map(str, e))))
    rows = []
    for e in entries[entries["contributes_z"]].itertuples():
        g = subst[(subst["uniprot_id"] == e.uniprot_id) & (subst["pdb_id"] == e.pdb_id)]
        c = g["difference_class"].value_counts()
        n_eng = int(c.get("ENGINEERED", 0))
        n_unc = int(c.get("UNCERTAIN", 0)) + int(c.get("MAPPING_REFERENCE_DIFFERENCE", 0))
        lists = [cmeta.get(e.pdb_id, {}).get("entity_mutation", {}).get(x, "") for x in ent.get((e.uniprot_id, e.pdb_id), [])]
        # substitutions the entity record declares but no observed residue confirms (residue
        # unobserved or outside the span); deletion notes ("DEL(165-215)") are construct deletions
        eng_labels = set(g.loc[g["difference_class"] == "ENGINEERED", "label"].str[0] +
                         g.loc[g["difference_class"] == "ENGINEERED", "label"].str[-1])
        declared = [f"{a}{n}{b}" for x in lists if not re.search(r"DEL|LACK|INSER", x.upper())
                    for a, n, b in MUT_TOKEN.findall(x.upper())]
        unmatched = [m for m in declared if m[0] + m[-1] not in eng_labels] if not n_eng else []
        status = ("ENGINEERED_MUTANT" if n_eng else "UNCERTAIN" if (n_unc or unmatched)
                  else "WT_NO_ENGINEERED_SUBSTITUTION")
        rows.append({"uniprot_id": e.uniprot_id, "receptor": e.common_name, "nr_code": e.nr_code, "pdb_id": e.pdb_id,
                     "n_observed_differences": len(g),
                     "n_substitutions": int(len(g) - c.get("MODIFIED_RESIDUE", 0)),
                     "n_engineered": n_eng, "n_modified_residue": int(c.get("MODIFIED_RESIDUE", 0)),
                     "n_natural_variant_or_isoform": int(c.get("NATURAL_VARIANT_OR_ISOFORM", 0)),
                     "n_mapping_reference_difference": int(c.get("MAPPING_REFERENCE_DIFFERENCE", 0)),
                     "n_uncertain": int(c.get("UNCERTAIN", 0)),
                     "engineered_substitutions": ";".join(g.loc[g["difference_class"] == "ENGINEERED", "label"]),
                     "engineered_evidence": ";".join(sorted(set(g.loc[g["difference_class"] == "ENGINEERED", "evidence"]))),
                     "uncertain_substitutions": ";".join(g.loc[g["difference_class"].isin(["UNCERTAIN", "MAPPING_REFERENCE_DIFFERENCE"]), "label"]),
                     "entity_mutation_annotation": " | ".join(x for x in lists if x),
                     "declared_substitutions_not_observed": ";".join(unmatched),
                     "final_classification": status, "classification_source": MUTATION_SOURCE})
    return pd.DataFrame(rows)


CELL_TOL, ANGLE_TOL = 0.02, 2.0


def crystal_forms(s: pd.DataFrame, cmeta: dict) -> dict:
    """Crystal form = same space group and unit cell within 2% (lengths) / 2 degrees (angles),
    grouped greedily per receptor in deposition order. Returns (uid, pdb) -> form id."""
    out = {}
    coord = s[s["has_coord"]].drop_duplicates(["uniprot_id", "pdb_id"])
    for uid, g in coord.groupby("uniprot_id"):
        reps: list[tuple] = []
        for pdb in sorted(g["pdb_id"], key=lambda p: (cmeta.get(p, {}).get("deposition_date", ""), p)):
            m = cmeta.get(pdb)
            if not m or not m["space_group"] or pd.isna(m["cell_a"]):
                continue
            cell = np.array([m["cell_a"], m["cell_b"], m["cell_c"]], dtype=float)
            ang = np.array([m["cell_alpha"], m["cell_beta"], m["cell_gamma"]], dtype=float)
            for i, (sg, c0, a0) in enumerate(reps):
                if sg == m["space_group"] and np.all(np.abs(cell - c0) / c0 <= CELL_TOL) and \
                        np.all(np.nan_to_num(np.abs(ang - a0)) <= ANGLE_TOL):
                    out[(uid, pdb)] = f"F{i + 1}"
                    break
            else:
                reps.append((m["space_group"], cell, ang))
                out[(uid, pdb)] = f"F{len(reps)}"
    return out


PARTNER_CATEGORY = {"COREGULATOR_PEPTIDE": "coregulator", "NR_PARTNER": "NR partner", "OTHER_PEPTIDE": "other partner",
                    "OTHER_PROTEIN": "other partner"}


def add_partner_identities(contacts: pd.DataFrame, s: pd.DataFrame) -> pd.DataFrame:
    """Normalized identity per partner chain (scripts/utils/partner_identity.py). DNA/RNA and other
    entities of the same receptor get no partner identity; neither does a fusion's own-receptor part."""
    names = dict(zip(s["uniprot_id"], s["common_name"]))
    first = [partner_identity.identities(d, a) for d, a in zip(contacts["partner_description"], contacts["partner_uniprot_json"])]
    learned: dict[str, set] = {}
    for ids, acc in zip(first, contacts["partner_uniprot_json"]):
        accs = json.loads(acc or "[]")
        if len(accs) == 1 and accs[0] not in partner_identity.ACCESSION_TO_GENE and len(ids) == 1 and ids[0][1] == "NAME":
            learned.setdefault(accs[0], set()).add(ids[0][0])
    learned1 = {a: next(iter(g)) for a, g in learned.items() if len(g) == 1}
    out_id, out_rule = [], []
    for r, d, a in zip(contacts.itertuples(), contacts["partner_description"], contacts["partner_uniprot_json"]):
        if r.catalogue_partner_class in ("NUCLEIC_ACID", "SAME_RECEPTOR_CHAIN"):
            out_id.append(""), out_rule.append("")
            continue
        ids = [(i, ru) for i, ru in partner_identity.identities(d, a, learned1) if i != names.get(r.uniprot_id)]
        out_id.append(";".join(i for i, _ in ids)), out_rule.append(";".join(ru for _, ru in ids))
    contacts = contacts.assign(partner_identity=out_id, partner_identity_rule=out_rule)
    # a fusion contributes several identities: one row per identity in the identity view
    return contacts.assign(partner_identity=contacts["partner_identity"].str.split(";")).explode("partner_identity") \
        .fillna({"partner_identity": ""}).reset_index(drop=True)


def add_contexts(entries: pd.DataFrame) -> pd.DataFrame:
    """Biological context per entry with z-scores (TRACE_LBD-only fragments are not an LBD
    context and get none). Several ligands or partners are joined in sorted order, so a
    combination is deterministic and synonyms cannot inflate counts (ligands are CCD ids,
    partners normalized identities)."""
    ctx_state, ctx_lig = [], []
    for e in entries.itertuples():
        if not e.contributes_z or e.trace_lbd_only:
            ctx_state.append(""), ctx_lig.append("")
            continue
        cats = sorted({PARTNER_CATEGORY[c] for c in e.partner_classes.split(";") if c in PARTNER_CATEGORY})
        partner_cat = " + ".join(cats) if cats else "standalone"
        partners = e.partner_identities.replace(";", " + ") if e.partner_identities else "standalone"
        state = e.ligand_state or "UNKNOWN"
        ligand = (e.functional_ligands.replace(";", " + ") if e.functional_ligands else "ligand not identified") \
            if state == "HOLO" else state
        ctx_state.append(f"{state} + {partner_cat}"), ctx_lig.append(f"{ligand} | {partners}")
    return entries.assign(state_partner_context=ctx_state, ligand_partner_context=ctx_lig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(V2 / "processed" / "structures.csv", dtype={"pdb_id": str}, keep_default_na=False, low_memory=False)
    s["has_coord"] = s["has_stage3b_coordinate_data"].astype(str) == "True"
    audit = pd.read_csv(V2 / "manifests" / "final_structure_selection_audit.csv", dtype={"pdb_id": str}, low_memory=False)
    res = pd.read_parquet(V2 / "processed" / "lbd_residue_bfactors_all.parquet",
                          columns=["uniprot_id", "pdb_id", "instance_id", "label_asym_id", "entity_label_seq_id",
                                   "canonical_uniprot_position", "is_construct_mapped", "ca_record_present",
                                   "is_representative_instance", "canonical_residue_name", "label_comp_ids_json"])
    rep = res[res["is_representative_instance"] & res["is_construct_mapped"]]

    # ---- partner contacts for every coordinate entry ----
    coord = s[s["has_coord"] & (s["representative_instance_id"] != "")]
    lbd_by_instance = {iid: (set(int(v) for v in g["entity_label_seq_id"].dropna()),
                             dict(zip(g["entity_label_seq_id"].astype("Int64"), g["canonical_uniprot_position"])))
                       for iid, g in rep.groupby("instance_id")}
    p27 = pd.read_csv(V2 / "manifests" / "final_polymer_partner_annotations.csv", dtype={"pdb_id": str}, keep_default_na=False)
    official_by_pdb: dict[str, dict] = {}
    for r in p27.itertuples():
        official_by_pdb.setdefault(r.pdb_id, {})[(r.receptor_uniprot_id, r.partner_polymer_entity_id, r.partner_chain)] = (
            r.partner_category, r.partner_category_evidence)
    jobs = []
    for pdb, g in coord.groupby("pdb_id"):
        cands = [{"uid": r.uniprot_id, "rep": r.representative_instance_id, "asym": r.label_asym_id,
                  "entity": str(r.entity_id)} for r in g.itertuples()]
        uid_by_entity = {str(r.entity_id): r.uniprot_id for r in s[s["pdb_id"] == pdb].itertuples()}
        jobs.append((pdb, cands, {c["rep"]: lbd_by_instance.get(c["rep"]) for c in cands}, uid_by_entity,
                     official_by_pdb.get(pdb, {})))
    with ProcessPoolExecutor() as ex:
        results = list(ex.map(partners_for_pdb, jobs, chunksize=8))
    contacts = pd.DataFrame([row for r in results for row in r["contacts"]])
    cmeta = {r["meta"]["pdb_id"]: r["meta"] for r in results}
    contacts = contacts.sort_values(["uniprot_id", "pdb_id", "partner_entity_id", "partner_label_asym_id"])
    contacts = add_partner_identities(contacts, s)
    contacts.to_csv(OUT / "catalogue_partner_contacts_v2.csv", index=False, lineterminator="\n")

    # ---- validation against scripts/27 on its own universe ----
    # every chain scripts/27 lists must be found by the catalogue's own geometry; the catalogue
    # additionally lists project-receptor entities (heterodimer partners) and non-primary entries
    key = lambda d, u, e, a: set(zip(d[u], d["pdb_id"], d[e], d[a]))  # noqa: E731
    a = key(contacts, "uniprot_id", "partner_entity_id", "partner_label_asym_id")
    b = key(p27, "receptor_uniprot_id", "partner_polymer_entity_id", "partner_chain")
    primary = set(zip(audit.loc[audit["in_primary_set"], "uniprot_id"], audit.loc[audit["in_primary_set"], "pdb_id"]))
    in_primary = {k for k in a if (k[0], k[1]) in primary}
    val = pd.DataFrame([
        {"check": "scripts/27 partner chains (primary-QC entries)", "n": len(b)},
        {"check": "of these, found by the catalogue contact computation", "n": len(a & b)},
        {"check": "catalogue chains in primary-QC entries", "n": len(in_primary)},
        {"check": "  of which project-receptor entities (skipped by scripts/27 by design)",
         "n": int(sum(1 for k in in_primary - b if k in set(key(contacts[contacts["partner_project_receptor"] != ""],
                                                               "uniprot_id", "partner_entity_id", "partner_label_asym_id"))))},
        {"check": "  other catalogue-only chains in primary-QC entries", "n": len(in_primary - b) - int(sum(
            1 for k in in_primary - b if k in set(key(contacts[contacts["partner_project_receptor"] != ""],
                                                      "uniprot_id", "partner_entity_id", "partner_label_asym_id"))))},
        {"check": "catalogue chains in non-primary coordinate entries (not covered by scripts/27)",
         "n": len(a) - len(in_primary)},
    ])
    val.to_csv(OUT / "catalogue_partner_validation_vs_script27.csv", index=False, lineterminator="\n")

    # ---- entry table ----
    prof = pd.read_parquet(V2 / "analysis" / "pdb_position_bfactor_profiles.parquet",
                           columns=["uniprot_id", "pdb_id", "b_factor_z", "ligand_state"])
    contrib = prof.dropna(subset=["b_factor_z"]).groupby(["uniprot_id", "pdb_id"])["ligand_state"].first()
    lig = pd.read_csv(V2 / "manifests" / "final_ligand_annotations.csv", dtype={"pdb_id": str}, keep_default_na=False)
    lig_ids = lig.groupby(["uniprot_id", "pdb_id"])["functional_ligand_component_ids_json"].agg(
        lambda c: sorted({x for v in c for x in json.loads(v or "[]")}))
    lig_names = {}
    for r in lig.itertuples():
        for i, n in zip(json.loads(r.functional_ligand_component_ids_json or "[]"), json.loads(r.functional_ligand_names_json or "[]")):
            lig_names.setdefault(i, n)
    obs = rep[rep["ca_record_present"]].copy()
    obs["observed_comp"] = obs["label_comp_ids_json"].map(lambda v: (json.loads(v) or [""])[0] if v else "")
    diff = obs[obs["observed_comp"] != obs["canonical_residue_name"]]
    subst = classify_substitutions(diff, s, cmeta)
    subst.to_csv(OUT / "catalogue_sequence_differences_v2.csv", index=False, lineterminator="\n")
    n_subst = subst[subst["difference_class"] != "MODIFIED_RESIDUE"].groupby(["uniprot_id", "pdb_id"]).size()
    by_class = subst.groupby(["uniprot_id", "pdb_id", "difference_class"]).size().unstack(fill_value=0)
    labels = subst[subst["difference_class"] != "MODIFIED_RESIDUE"].groupby(["uniprot_id", "pdb_id"])["label"].agg(";".join)
    forms = crystal_forms(s, cmeta)
    # construct = the entity's full sequence-evidence alignment to UniProt (not clipped to the span)
    aln = pd.read_csv(V2 / "inputs" / "rcsb_sequence_coordinate_alignments_v2.csv", dtype={"pdb_id": str},
                      keep_default_na=False)
    aln = aln[aln["polymer_entity_id"].isin(set(s.loc[s["has_coord"], "polymer_entity_id"]))]
    reg = aln.assign(r=aln["aligned_regions_json"].map(json.loads))
    reg = reg.assign(first=reg["r"].map(lambda x: min(i["query_begin"] for i in x)),
                     last=reg["r"].map(lambda x: max(i["query_end"] for i in x)), nseg=reg["r"].map(len))
    construct = reg.groupby(["uniprot_id", "pdb_id"]).agg(min=("first", "min"), max=("last", "max"), nseg=("nseg", "max"))
    span = pd.read_csv(V2 / "inputs" / "analysis_span_v2.csv").set_index("uniprot_id")
    obs_frac = obs.groupby(["uniprot_id", "pdb_id"])["canonical_uniprot_position"].nunique()
    in_primary = audit.groupby(["uniprot_id", "pdb_id"])["in_primary_set"].any()
    pc = contacts.groupby(["uniprot_id", "pdb_id"])["catalogue_partner_class"].agg(lambda c: sorted(set(c)))
    pid = contacts[contacts["partner_identity"] != ""].groupby(["uniprot_id", "pdb_id"])["partner_identity"].agg(
        lambda c: sorted(set(c)))

    rows = []
    for (uid, pdb), g in s.groupby(["uniprot_id", "pdb_id"]):
        key = (uid, pdb)
        cg = g[g["has_coord"]]
        trace_only = len(cg) > 0 and set(cg["lbd_mapping_category"]) == {"TRACE_LBD"}
        states = sorted(set(cg["ligand_state"])) if len(cg) else []
        state = contrib.get(key) if key in contrib.index else (states[0] if len(states) == 1 else ("AMBIGUOUS" if states else ""))
        classes = pc.get(key, []) if len(cg) else []
        rows.append({
            "uniprot_id": uid, "nr_code": g["nr_code"].iloc[0], "common_name": g["common_name"].iloc[0],
            "group": g["group"].iloc[0], "pdb_id": pdb, "n_receptor_entities": len(g),
            "experimental_method": g["experimental_method"].iloc[0],
            "has_coordinate_data": len(cg) > 0,
            "no_coordinate_reason": "" if len(cg) else (
                "not X-ray (NMR/EM)" if not (g["method_is_xray"].astype(str) == "True").any() else
                "receptor chain outside the LBD" if set(g["lbd_mapping_category"]) == {"NO_LBD_OVERLAP"} else
                "excluded at the coordinate stage"),
            "contributes_z": key in contrib.index, "trace_lbd_only": trace_only,
            "in_primary_qc_set": bool(in_primary.get(key, False)),
            "ligand_state": state or "", "catalogue_state": "TRACE_LBD" if trace_only else (state or ""),
            "functional_ligands": ";".join(lig_ids.get(key, [])),
            "functional_ligand_names": ";".join(lig_names.get(i, "") for i in lig_ids.get(key, [])),
            "resolution": pd.to_numeric(g["resolution"], errors="coerce").min(),
            "r_free": pd.to_numeric(g["r_free"], errors="coerce").min(),
            "fusion_or_chimera": bool((g["fusion_or_chimera"].astype(str) == "True").any()),
            "identity_class": ";".join(sorted(set(g["identity_class"]))),
            "construct_first_position": construct["min"].get(key, np.nan) if key in construct.index else np.nan,
            "construct_last_position": construct["max"].get(key, np.nan) if key in construct.index else np.nan,
            "construct_n_aligned_segments": construct["nseg"].get(key, np.nan) if key in construct.index else np.nan,
            "observed_span_fraction": round(obs_frac.get(key, 0) / (span.loc[uid, "analysis_span_end"]
                                                                     - span.loc[uid, "analysis_span_start"] + 1), 3)
            if len(cg) else np.nan,
            "n_point_substitutions": int(n_subst.get(key, 0)) if len(cg) else np.nan,
            **{f"n_diff_{c.lower()}": (int(by_class.loc[key, c]) if key in by_class.index and c in by_class.columns else 0)
               if len(cg) else np.nan for c in DIFF_CLASSES},
            "substitutions": labels.get(key, "") if len(cg) else "",
            **({k: v for k, v in cmeta[pdb].items() if k not in ("pdb_id", "copies_by_entity", "seq_dif", "entity_mutation")}
               if pdb in cmeta else {}),
            "deposition_year": int(cmeta[pdb]["deposition_date"][:4]) if pdb in cmeta and cmeta[pdb]["deposition_date"] else np.nan,
            "n_receptor_copies_asu": sum(cmeta[pdb]["copies_by_entity"].get(str(e), 0) for e in cg["entity_id"])
            if pdb in cmeta else np.nan,
            "crystal_form_id": forms.get(key, ""),
            "partner_classes": ";".join(classes),
            "partner_identities": ";".join(pid.get(key, [])) if len(cg) else "",
            "standalone": (len(classes) == 0) if len(cg) else np.nan,
        })
    entries = add_contexts(pd.DataFrame(rows).sort_values(["nr_code", "pdb_id"]))
    audit_tbl = mutation_audit(subst, entries, cmeta, s)
    audit_tbl.to_csv(OUT / "catalogue_mutation_audit_v2.csv", index=False, lineterminator="\n")
    entries = entries.merge(audit_tbl[["uniprot_id", "pdb_id", "final_classification"]].rename(
        columns={"final_classification": "mutation_status"}), on=["uniprot_id", "pdb_id"], how="left")
    entries.to_csv(OUT / "catalogue_entry_context_v2.csv", index=False, lineterminator="\n")
    print(f"entries {len(entries)}; coordinate {int(entries['has_coordinate_data'].sum())}; partner contacts {len(contacts)}")
    print(val.to_string(index=False))


if __name__ == "__main__":
    main()
