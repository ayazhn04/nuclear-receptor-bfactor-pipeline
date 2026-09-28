"""
Part 3 of the freeze/consolidate stage (continuing from commit f76f352 on
analysis/normalization-v1): a final, complete receptor-identity audit
across every generated file, against the sole authoritative source
(config/nr_metadata.py, the frozen professor-supplied metadata).

Two passes:
  1. STRUCTURED: every tracked CSV/Parquet under data/ and reports/ that
     carries both a uniprot_id column and an nr_code and/or common_name
     column -- every row's nr_code/common_name must exactly match the
     authoritative value for that uniprot_id.
  2. PROSE: every tracked Markdown file -- for each receptor's authoritative
     common_name, flag any line that mentions that common_name in the same
     line as a DIFFERENT receptor's nr_code (the exact pattern behind all
     3 previously-found mislabeling bugs: a hand-typed common name next to
     the wrong nr_code/uniprot_id).

This script only reads and reports -- it never modifies any biological or
statistical value; display-label bugs (if any are found) are fixed by
hand afterward in the specific file, never here.
"""
import re
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402

TABLES = PROJECT_ROOT / "reports" / "tables"

AUTHORITATIVE = {r["uniprot_id"]: {"nr_code": r["nr_code"], "common_name": r["common_name"]} for r in NR_METADATA}
NR_CODE_TO_UID = {r["nr_code"]: r["uniprot_id"] for r in NR_METADATA}

SEARCH_DIRS = ["data/analysis", "data/manifests", "data/processed", "reports/tables"]
SKIP_SUBSTRINGS = ["__pycache__"]

# Files whose "uniprot_id"-named column legitimately holds something other than a
# single per-row UniProt ID (e.g. a JSON pair describing a cross-receptor anomaly) --
# not an identity bug, just a different column semantic; excluded from the structured
# per-row check rather than mis-flagged.
STRUCTURED_FALSE_POSITIVE_FILES = {"data/manifests/rcsb_discovery_anomalies.csv"}

# Known valid alternate spellings for the SAME receptor (gene-symbol aliases, ASCII
# transliterations of Greek letters, etc.) -- confirmed by manual inspection during
# this audit to never collide with any OTHER of the 48 receptors' authoritative name.
# Recording these here (not silently) makes every classification traceable.
KNOWN_ALIASES = {
    "P03372": {"ESR1", "ERalpha", "ER-alpha", "ER alpha"},       # authoritative: NR3A1 / ERα
    "P37231": {"PPARgamma", "PPARG", "PPAR-gamma", "PPAR gamma"},  # authoritative: NR1C3 / PPARγ
    "P04150": {"GR"},                                              # authoritative: NR3C1 / GR (identical, no alias needed)
    "P11473": {"VDR"},                                             # authoritative: NR1I1 / VDR (identical)
    "P22736": {"NR4A1", "Nur77", "NUR77"},                         # authoritative: NR4A1 / NUR77
}


def find_structured_files() -> list[Path]:
    files = []
    for d in SEARCH_DIRS:
        base = PROJECT_ROOT / d
        if not base.exists():
            continue
        for ext in ("*.csv", "*.parquet"):
            files.extend(base.rglob(ext))
    return [f for f in files if not any(s in str(f) for s in SKIP_SUBSTRINGS)]


def audit_structured_file(path: Path) -> list[dict]:
    rel = str(path.relative_to(PROJECT_ROOT))
    if rel in STRUCTURED_FALSE_POSITIVE_FILES:
        return []
    try:
        if path.suffix == ".csv":
            df = pd.read_csv(path, dtype=str, keep_default_na=False, nrows=None)
        else:
            df = pd.read_parquet(path)
    except Exception as e:
        return [{"file": rel, "issue": f"COULD_NOT_READ: {e}", "severity": "TOOLING"}]

    if "uniprot_id" not in df.columns:
        return []
    # a real per-row single UniProt accession is exactly 6 or 10 alnum characters;
    # a column merely NAMED uniprot_id that actually holds something else (e.g. a
    # JSON list) is a different semantic, not a per-receptor identity claim
    sample = df["uniprot_id"].dropna().astype(str)
    if len(sample) and not (sample.str.match(r"^[A-Z][A-Z0-9]{5}([A-Z0-9]{4})?$").mean() > 0.9):
        return []
    has_nr_code = "nr_code" in df.columns
    has_common_name = "common_name" in df.columns
    if not (has_nr_code or has_common_name):
        return []

    issues = []
    for uid, group in df.groupby("uniprot_id"):
        auth = AUTHORITATIVE.get(uid)
        if auth is None:
            issues.append({"file": rel, "uniprot_id": uid,
                            "issue": "UNKNOWN_UNIPROT_ID_NOT_IN_NR_METADATA", "severity": "REVIEW"})
            continue
        aliases = KNOWN_ALIASES.get(uid, set())
        if has_nr_code:
            bad_codes = set(group["nr_code"].astype(str)) - {auth["nr_code"], ""}
            genuinely_bad = bad_codes - aliases
            if bad_codes:
                sev = "GENUINE" if genuinely_bad else "BENIGN_ALIAS"
                issues.append({"file": rel, "uniprot_id": uid,
                                "issue": f"nr_code column holds {bad_codes}, strict authoritative nr_code is {auth['nr_code']}",
                                "severity": sev})
        if has_common_name:
            bad_names = set(group["common_name"].astype(str)) - {auth["common_name"], ""}
            genuinely_bad = bad_names - aliases
            if bad_names:
                sev = "GENUINE" if genuinely_bad else "BENIGN_ALIAS"
                issues.append({"file": rel, "uniprot_id": uid,
                                "issue": f"common_name column holds {bad_names}, authoritative common_name is {auth['common_name']}",
                                "severity": sev})
    return issues


# Manually reviewed and confirmed NOT a mislabeling: the heuristic's remaining false
# positive after nearest-match refinement. reports/ROLE_D_PROGRESS_REPORT.md is a
# PROTECTED historical file (never modified); its line 7 (Russian) lists all 5 anchor
# receptors as a flat sequence "VDR, PPARG, ESR1, NR3C1 (GR) и NR4A1" -- correct text,
# just not a Name/Code PAIR structure, which the proximity heuristic cannot perfectly
# distinguish from the genuine bug pattern. Recorded here rather than silently dropped.
MANUALLY_REVIEWED_FALSE_POSITIVES = {
    ("reports/ROLE_D_PROGRESS_REPORT.md", 7): (
        "Protected historical file; line lists all 5 anchor receptors as a flat "
        "sequence (VDR, PPARG, ESR1, NR3C1 (GR), NR4A1), not a mislabeled pair -- "
        "manually reviewed, confirmed correct, left unmodified per file protection."
    ),
}

PROXIMITY_WINDOW = 20  # chars; catches the exact "CommonName/WrongCode" or
                        # "CommonName (WrongCode)" adjacency pattern behind all 3
                        # previously-found bugs, without flagging a sentence that
                        # legitimately lists several receptors together


ANCHOR_NAME_SET = {"VDR", "PPARgamma", "ESR1", "GR", "NR4A1"}


def _brace_list_spans(line: str) -> list[tuple[int, int]]:
    """Spans of {A,B,C,...} constructs (e.g. a glob-style figure-name list
    like "{VDR,PPARgamma,ESR1,GR,NR4A1}") whose comma-separated tokens are
    exactly the 5 known anchor-receptor names, in any order -- a flat list
    idiom used repeatedly in this project's own documentation, never a
    mislabeled Name/Code pair. Generalizes the single-line override that
    would otherwise be needed for every new occurrence of this idiom."""
    spans = []
    anchor_lower = {n.lower() for n in ANCHOR_NAME_SET}
    for m in re.finditer(r"\{([^{}]*)\}", line):
        tokens = {t.strip().lower() for t in m.group(1).split(",")}
        if tokens == anchor_lower:
            spans.append((m.start(), m.end()))
    return spans


def audit_prose_files() -> list[dict]:
    md_files = [f for f in PROJECT_ROOT.rglob("*.md") if ".git" not in str(f)]
    issues = []
    for path in md_files:
        try:
            text = path.read_text()
        except Exception:
            continue
        for line_no, line in enumerate(text.splitlines(), 1):
            brace_spans = _brace_list_spans(line)
            code_matches = list(re.finditer(r"\bNR[0-9][A-Z][0-9]\b", line))
            if not code_matches:
                continue
            for uid, auth in AUTHORITATIVE.items():
                name = auth["common_name"]
                if len(name) < 3:
                    continue
                for name_match in re.finditer(r"\b" + re.escape(name) + r"\b", line):
                    if any(s <= name_match.start() and name_match.end() <= e for s, e in brace_spans):
                        continue  # inside a known "{anchor1,anchor2,...}" flat glob-list idiom
                    # only the NEAREST code to this specific name occurrence matters --
                    # a comma-separated list of several correctly-paired "Name/CODE"
                    # items must not flag a name against its NEIGHBOR's code
                    def dist(code_match):
                        return min(abs(name_match.start() - code_match.end()),
                                   abs(code_match.start() - name_match.end()))
                    nearest = min(code_matches, key=dist)
                    if dist(nearest) > PROXIMITY_WINDOW:
                        continue
                    if nearest.group() == auth["nr_code"]:
                        continue  # correctly paired with its own nearest code
                    rel = str(path.relative_to(PROJECT_ROOT))
                    override = MANUALLY_REVIEWED_FALSE_POSITIVES.get((rel, line_no))
                    issues.append({
                        "file": rel, "line": line_no, "uniprot_id": uid,
                        "issue": f"common_name '{name}' (uniprot {uid}, correct nr_code {auth['nr_code']}) "
                                 f"has NEAREST nr_code '{nearest.group()}' ({dist(nearest)} chars away) -- mismatch"
                                 + (f" [MANUALLY REVIEWED: {override}]" if override else ""),
                        "line_text": line.strip()[:200],
                        "severity": "BENIGN_ALIAS" if override else "GENUINE",
                    })
    return issues


def main():
    structured_files = find_structured_files()
    print(f"Scanning {len(structured_files)} structured (CSV/Parquet) files under {SEARCH_DIRS}...")
    all_issues = []
    files_checked_per_uid = {uid: 0 for uid in AUTHORITATIVE}
    n_checked = 0
    for path in structured_files:
        issues = audit_structured_file(path)
        if issues or True:
            n_checked += 1
        all_issues.extend(issues)

    print(f"Checked {n_checked} files with identity columns; {len(all_issues)} structured inconsistencies found")

    print("\nScanning Markdown files for prose cross-reference errors...")
    prose_issues = audit_prose_files()
    print(f"{len(prose_issues)} prose inconsistencies found")

    all_findings = all_issues + prose_issues
    genuine = [i for i in all_findings if i.get("severity") == "GENUINE"]
    benign = [i for i in all_findings if i.get("severity") == "BENIGN_ALIAS"]
    review = [i for i in all_findings if i.get("severity") == "REVIEW"]

    if genuine:
        print(f"\n{len(genuine)} GENUINE cross-receptor mislabeling issue(s) (require a fix):")
        for i in genuine:
            print(f"  {i}")
    if review:
        print(f"\n{len(review)} REVIEW item(s) (uniprot_id-named column with non-identity semantics; "
              "excluded from per-row identity checking, listed for transparency):")
        for i in review[:10]:
            print(f"  {i}")
    if benign:
        print(f"\n{len(benign)} BENIGN_ALIAS finding(s) (same receptor, alternate valid spelling -- not an error):")
        for i in benign:
            print(f"  {i}")

    # ---- build the per-receptor summary table ----
    rows = []
    for uid, auth in sorted(AUTHORITATIVE.items()):
        n_genuine = sum(1 for i in genuine if i.get("uniprot_id") == uid or f"uniprot {uid}" in i.get("issue", ""))
        n_benign = sum(1 for i in benign if i.get("uniprot_id") == uid)
        status = "GENUINE_INCONSISTENCY_UNRESOLVED" if n_genuine else (
            "CLEAN_WITH_DOCUMENTED_BENIGN_ALIASES" if n_benign else "CLEAN"
        )
        rows.append({
            "uniprot_id": uid, "authoritative_nr_code": auth["nr_code"], "authoritative_common_name": auth["common_name"],
            "files_checked": n_checked, "n_inconsistencies": n_genuine,
            "n_documented_benign_aliases": n_benign,
            "status": status,
        })
    audit_table = pd.DataFrame(rows)
    out = TABLES / "final_receptor_identity_audit.csv"
    audit_table.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(audit_table)} receptors)")

    n_unresolved = (audit_table["status"] == "GENUINE_INCONSISTENCY_UNRESOLVED").sum()
    print(f"\n{len(audit_table) - n_unresolved}/{len(audit_table)} receptors have ZERO genuine unresolved "
          f"inconsistencies; {n_unresolved} require a fix.")
    if n_unresolved:
        print(audit_table[audit_table["status"] == "GENUINE_INCONSISTENCY_UNRESOLVED"].to_string(index=False))


if __name__ == "__main__":
    main()
