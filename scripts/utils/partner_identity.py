"""
Partner / cofactor identity normalization for the receptor catalogue.

The deposited entity description and UniProt accession of a partner chain are kept
unchanged in the catalogue tables; this module only adds a normalized identity
(gene-level name) so that naming variants of one protein are counted once.

Evidence order (first that applies wins; the rule used is recorded):
  1 ACCESSION   explicit UniProt accession -> gene (human entries and their orthologs,
                listed below; every entry was checked against its deposited descriptions)
  2 NAME        description pattern -> gene (synonyms listed below)
  3 MOTIF       synthetic NR-box peptides named only by sequence that contain the core
                ILHRLLQ are one identity ("HKILHRLLQ NR-box peptide"); no gene is claimed
  4 GENERIC     generic labels with no protein named ("co-activator peptide", "unknown
                peptide", "stapled peptide", ...) -> one identity, "unidentified peptide"
  5 LEARNED     an accession not in the table that every NAME-matched row maps to one gene
  6 FALLBACK    the cleaned description itself (uncertain identities are never merged)
A chain whose accessions name two proteins (a fusion, e.g. RXRa-NCOA1) contributes both.
Nucleic acids are never a partner identity; another entity of the same receptor is not
counted as a partner.
"""
from __future__ import annotations

import json
import re

ACCESSION_TO_GENE = {
    # coregulators
    "Q15788": "NCOA1", "P70365": "NCOA1",                  # SRC-1 (human, mouse)
    "Q15596": "NCOA2", "Q61026": "NCOA2",                  # SRC-2 / TIF2 / GRIP1 (human, mouse)
    "Q9Y6Q9": "NCOA3", "Q13772": "NCOA4", "Q9HCD5": "NCOA5",
    "Q9UBK2": "PPARGC1A", "Q86YN6": "PPARGC1B",
    "O75376": "NCOR1", "Q9Y618": "NCOR2",
    "Q15648": "MED1", "P48552": "NRIP1", "Q92793": "CREBBP",
    # nuclear receptors (project receptors and orthologs)
    "P19793": "RXRα", "P28700": "RXRα", "P28702": "RXRβ", "P37231": "PPARγ", "Q96RI1": "FXR",
    "Q13133": "LXRα", "Q9Z0Y9": "LXRα", "Q91X41": "LXRα", "P55055": "LXRβ", "F1D8P7": "LXRβ",
    "Q14994": "CAR", "O35627": "CAR", "Q15466": "SHP", "P97947": "SHP", "O00482": "LRH-1",
    "P10828": "TRβ", "P04625": "TRα", "P10826": "RARβ", "P10275": "AR", "O75469": "PXR",
    "F1D8P4": "DAX-1", "A0A1X8XLE9": "ancestral GR (AncGR2)",
    # other proteins
    "P63104": "14-3-3ζ", "P31947": "14-3-3σ", "P35222": "β-catenin", "P02751": "fibronectin",
    "Q15653": "NFKBIB", "Q9UQ80": "PA2G4", "P41223": "BUD31", "Q96L73": "NSD1", "C9JGY8": "JAZF1",
}
NAME_RULES = [  # order matters: PGC-1 before NCOA1 ("coactivator 1-alpha" vs "coactivator 1")
    (r"COACTIVATOR 1[- ]?BETA", "PPARGC1B"),
    (r"PGC-?1 ?(A|ALPHA|ALFA)\b|PGC1ALPHA|GAMMA COACTIVATOR[- ]1[- ]?(A\b|ALPHA)|COACTIVATOR-1A", "PPARGC1A"),
    (r"NUCLEAR RECEPTOR COACTIVATOR 2|NCOA-?2|\bSRC-?2|SCR2-2|GRIP-?1|GRIP PEPTIDE|\bTIF-?2|"
     r"TRANSCRIPTION(AL)? INTERMEDIARY FACTOR[- ]2|GLUCOCORTICOID RECEPTOR[- ]INTERACTING PROTEIN 1", "NCOA2"),
    (r"NUCLEAR RECEPTOR COACTIVATOR 1\b|NCOA-?1|\bSRC-?1|STEROID RECEPTOR COACTIVATOR[- ]1", "NCOA1"),
    (r"NUCLEAR RECEPTOR COACTIVATOR 3|\bSRC-?3|STEROID RECEPTOR COACTIVATOR-3|\bRAC3|\bACTR\b|\bAIB1", "NCOA3"),
    (r"NUCLEAR RECEPTOR CO-?REPRESSOR 2|\bSMRT", "NCOR2"),
    (r"NUCLEAR RECEPTOR CO-?REPRESSOR 1\b", "NCOR1"),
    (r"MEDIATOR OF RNA POLYMERASE II TRANSCRIPTION SUBUNIT 1|TRAP220", "MED1"),
    (r"NUCLEAR RECEPTOR-INTERACTING PROTEIN 1|RIP140", "NRIP1"),
    (r"ARA70", "NCOA4"),
]
GENERIC = re.compile(r"^(CO-?REGULATOR PEPTIDE|CO-?ACTIVATOR( PEPTIDE)?|COREPRESSOR PEPTIDE|ACTIVATOR PEPTIDE|"
                     r"LXXLL MOTIF COACTIVATOR|PEPTIDE FROM NUCLEAR RECEPTOR COACTIVATOR|UNKNOWN PEPTIDE|"
                     r"STAPLED PEPTIDE(-LIKE LIGAND| \d+\w*)?|PEPTIDE (INHIBITOR|ANTAGONIST)|\d+-MER PEPTIDE|"
                     r"SRC PEPTIDE( FRAGMENT)?|D30 PEPTIDE|AB5 PEPTIDE|PG\d+(-NL)?)$")
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}


def _sequence(desc: str) -> str | None:
    toks = re.findall(r"[A-Z0-9]{3}", desc.replace("PEPTIDE", ""))
    if desc.count("-") >= 3 and toks and sum(t in THREE for t in toks) >= 0.7 * len(toks):
        return "".join(THREE.get(t, "X") for t in toks)
    if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,}( PEPTIDE)?", desc):
        return desc.split()[0]
    return None


def _clean(desc: str) -> str:
    return re.sub(r"\s+", " ", desc.strip().upper())


def identities(description: str, accessions_json: str, learned: dict | None = None) -> list[tuple[str, str]]:
    """[(identity, rule)] for one partner chain."""
    desc = _clean(description)
    accs = json.loads(accessions_json or "[]")
    out = [(ACCESSION_TO_GENE[a], "ACCESSION") for a in accs if a in ACCESSION_TO_GENE]
    if out and len(out) == len(accs):
        return sorted(set(out))
    for pat, gene in NAME_RULES:
        if re.search(pat, desc):
            return sorted(set(out + [(gene, "NAME")]))
    seq = _sequence(desc)
    if seq and "ILHRLLQ" in seq:
        return sorted(set(out + [("HKILHRLLQ NR-box peptide", "MOTIF")]))
    if GENERIC.match(desc):
        return sorted(set(out + [("unidentified peptide", "GENERIC")]))
    for a in accs:
        if learned and a in learned:
            out.append((learned[a], "LEARNED"))
    return sorted(set(out)) if out else [(desc.lower(), "FALLBACK")]
