#!/usr/bin/env python3
"""Classify how two VEP Storable cache builds differ on the transcripts whose content differs.

Reads two dumps written by `dump_storable_transcripts.pl` for the same transcript ids and,
per transcript, names the field groups that differ once number quoting, `dbID` and
`_gene_version` are neutralised (the same normalisation as the content digest):
`peptide` (the stored translation), `predictions` (SIFT/PolyPhen matrices present in one
build only), `protein_features` (domain annotations), `other` (anything else, with the
key). Writes one TSV row per transcript and a summary JSON with the count per class.

With --json-cache, the vep-rs JSON cache (its `transcripts/<contig>/*.json[.gz]` shards) is
read too and each transcript's stored peptide is matched against build A's and build B's,
so the column `json_cache_peptide_matches` says which build the conversion came from.

Usage:
  classify_storable_transcript_differences.py --version-a 113 --dump-a a.txt --version-b 115 --dump-b b.txt --out-dir <dir>
      [--json-cache <cache dir with transcripts/> [--json-cache-label <name>]]
"""

from __future__ import annotations

import argparse
import collections
import difflib
import gzip
import json
import re
import sys
from pathlib import Path

PEPTIDE_KEYS = {"peptide", "translateable_seq"}
PREDICTION_KEYS = {"sift", "polyphen_humdiv", "polyphen_humvar", "protein_function_predictions", "matrix", "matrix_compressed", "peptide_length", "sub_analysis", "translation_md5"}
FEATURE_KEYS = {"hseqname", "_display_label", "analysis", "protein_features", "hstart", "hend", "score", "percent_id", "p_value", "logic_name", "db", "db_version", "interpro_ac", "idesc", "ilabel"}


def load(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for chunk in path.read_text(encoding="utf-8", errors="replace").split("### ")[1:]:
        tid, body = chunk.split("\n", 1)
        lines = []
        for line in body.split("\n"):
            line = re.sub(r"=> '(-?\d+)'(,?)$", r"=> \1\2", line)
            if "'dbID' =>" in line or "'_gene_version' =>" in line:
                continue
            lines.append(line)
        out[tid.strip()] = lines
    return out


def peptide_of(dump_lines: list[str]) -> str | None:
    for line in dump_lines:
        m = re.search(r"'peptide' => '([A-Z*]*)'", line)
        if m:
            return m.group(1)
    return None


def json_cache_peptides(cache: Path, ids: set[str]) -> dict[str, str | None]:
    found: dict[str, str | None] = {}
    shards = list((cache / "transcripts").glob("*/*.json")) + list((cache / "transcripts").glob("*/*.json.gz"))
    for shard in sorted(shards):
        opener = gzip.open if shard.name.endswith(".gz") else open
        with opener(shard, "rt", encoding="utf-8") as fh:
            data = json.load(fh)
        transcripts = data if isinstance(data, list) else [t for v in data.values() if isinstance(v, list) for t in v]
        for t in transcripts:
            if not isinstance(t, dict):
                continue
            sid = (t.get("stable_id") or "").split(".")[0]
            if sid in ids and sid not in found:
                found[sid] = (t.get("variation_effect_feature_cache") or {}).get("peptide")
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version-a", required=True)
    ap.add_argument("--dump-a", required=True, type=Path)
    ap.add_argument("--version-b", required=True)
    ap.add_argument("--dump-b", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--json-cache", type=Path)
    ap.add_argument("--json-cache-label", help="name recorded for the JSON cache (default: the directory's name)")
    a = ap.parse_args()
    da, db = load(a.dump_a), load(a.dump_b)
    ids = sorted(set(da) & set(db))
    missing = sorted((set(da) | set(db)) - set(ids))
    if missing:
        print(f"ERROR: [classify_storable_transcript_differences] {len(missing)} ids are in one dump only: {missing[:5]}", file=sys.stderr)
        return 1
    a.out_dir.mkdir(parents=True, exist_ok=True)
    per_class: collections.Counter[str] = collections.Counter()
    combos: collections.Counter[str] = collections.Counter()
    json_pep = json_cache_peptides(a.json_cache, set(ids)) if a.json_cache else {}
    json_match: collections.Counter[str] = collections.Counter()
    rows = []
    for tid in ids:
        diff = [l for l in difflib.unified_diff(da[tid], db[tid], lineterm="", n=0) if l[:1] in "+-" and not l.startswith(("+++", "---"))]
        classes: set[str] = set()
        other_keys: set[str] = set()
        # A differing value line without a key belongs to the key of the nearest preceding
        # keyed line in the same hunk; the peptide and matrix values are single lines, so
        # the keyed line itself carries them.
        for l in diff:
            m = re.search(r"'([A-Za-z_0-9]+)' =>", l)
            if not m:
                continue
            k = m.group(1)
            if k in PEPTIDE_KEYS:
                classes.add("peptide")
            elif k in PREDICTION_KEYS:
                classes.add("predictions")
            elif k in FEATURE_KEYS or k in ("start", "end"):
                classes.add("protein_features")
            else:
                classes.add("other")
                other_keys.add(k)
        if not diff:
            classes.add("identical_after_normalisation")
        for c in classes:
            per_class[c] += 1
        combos["+".join(sorted(classes))] += 1
        match = ""
        if a.json_cache:
            pa, pb, pj = peptide_of(da[tid]), peptide_of(db[tid]), json_pep.get(tid)
            if pj is None:
                match = "absent"
            elif pa == pb:
                match = "both" if pj == pa else "neither"
            elif pj == pa:
                match = a.version_a
            elif pj == pb:
                match = a.version_b
            else:
                match = "neither"
            json_match[match] += 1
        rows.append((tid, "+".join(sorted(classes)), ",".join(sorted(other_keys)), match))
    with open(a.out_dir / f"differing_{a.version_a}_{a.version_b}_classes.tsv", "w", encoding="utf-8") as fh:
        fh.write("transcript_id\tclasses\tother_keys\tjson_cache_peptide_matches\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    summary = {"version_a": a.version_a, "version_b": a.version_b, "transcripts": len(ids),
               "per_class": dict(per_class), "class_combinations": dict(combos)}
    if a.json_cache:
        summary["json_cache"] = a.json_cache_label or a.json_cache.name
        summary["json_cache_peptide_matches"] = dict(json_match)
    (a.out_dir / f"differing_{a.version_a}_{a.version_b}_classes.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[classify_storable_transcript_differences] {len(ids)} transcripts: " + ", ".join(f"{k} {v}" for k, v in sorted(per_class.items())) + "; combinations " + json.dumps(dict(combos))
          + (f"; JSON cache peptide matches {json.dumps(dict(json_match))}" if a.json_cache else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
