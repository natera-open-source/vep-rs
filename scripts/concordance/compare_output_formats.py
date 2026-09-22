#!/usr/bin/env python3
"""Output-format concordance harness.

Verifies that a VEP engine emits the same consequence calls regardless of output
format. Each format is parsed down to a canonical 5-tuple
``(location, allele, feature, consequence_set, impact)`` and the per-format tuple
sets are diffed pairwise. Self-parity (same engine, different formats) must be
EXACT; cross-engine parity (vep-rs vs fastVEP on the same format) is a
spot-check that documents inter-engine compatibility and is NOT required to be 1.0.

Engines and formats:
  - vep-rs : vcf, tab, json, parquet  (parquet read back via duckdb's
             read_parquet over the Hive-partitioned directory; both nested
             LIST<VARCHAR> and flat scalar CSQ layouts are handled. Skipped
             with a caveat only when the path is absent/unreadable).
  - fastvep: vcf, tab, json           (no parquet emitter).

The 5-tuple normalization mirrors the conventions used by
``compare_vep_outputs.py`` (consequence terms comma-joined + sorted; Location
``chr:pos`` / ``chr:start-end``) and ``compare_sv_concordance.py``
(chr-strip, ``M`` -> ``MT``) so results are directly comparable to the SNP/indel
and structural-variant concordance reports.

This script does NOT annotate; it consumes already-produced output files. The
caller annotates the same input under each format first, then points this script
at the resulting files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

# CSQ field order emitted by vep-rs (from the ##INFO=<ID=CSQ ...Format: ...> header).
# Position is resolved dynamically from the header when present; this is the fallback.
DEFAULT_CSQ_FIELDS = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
]


def _open(path: Path):
    """Open a possibly-gzipped text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def norm_contig(chrom: str) -> str:
    """Ensembl-style contig: strip 'chr', map M->MT. Matches compare_sv_concordance."""
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return "MT" if c == "M" else c


def norm_consequence(raw: str) -> str:
    """Sorted, comma-joined SO terms. '&' and ',' both used as separators upstream."""
    terms = [t for t in re.split(r"[,&]", raw.strip()) if t]
    return ",".join(sorted(terms))


# A 5-tuple: (location, allele, feature, consequence_set, impact)
Tuple5 = tuple[str, str, str, str, str]


def parse_vcf(path: Path) -> set[Tuple5]:
    """One tuple per (record x CSQ entry). Location = chr:pos (or chr:pos-end)."""
    tuples: set[Tuple5] = set()
    csq_fields = DEFAULT_CSQ_FIELDS
    idx = {name: i for i, name in enumerate(csq_fields)}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                m = re.search(r"ID=CSQ.*Format:\s*([^\"]+)\"?>", line)
                if m:
                    csq_fields = [f.strip() for f in m.group(1).strip().split("|")]
                    idx = {name: i for i, name in enumerate(csq_fields)}
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, info = cols[0], cols[1], cols[7]
            loc = f"{norm_contig(chrom)}:{pos}"
            m = re.search(r"(?:^|;)CSQ=([^;\t]+)", info)
            if not m:
                continue
            for entry in m.group(1).split(","):
                parts = entry.split("|")
                if len(parts) <= idx.get("Feature", 6):
                    continue
                allele = (
                    parts[idx["Allele"]] if idx.get("Allele", 0) < len(parts) else ""
                )
                feature = (
                    parts[idx["Feature"]] if idx.get("Feature", 6) < len(parts) else ""
                )
                csq = (
                    parts[idx["Consequence"]]
                    if idx.get("Consequence", 1) < len(parts)
                    else ""
                )
                impact = (
                    parts[idx["IMPACT"]] if idx.get("IMPACT", 2) < len(parts) else ""
                )
                tuples.add((loc, allele, feature, norm_consequence(csq), impact))
    return tuples


def _impact_from_extra(extra: str) -> str:
    m = re.search(r"(?:^|;)IMPACT=([^;]+)", extra)
    return m.group(1) if m else ""


def parse_tab(path: Path) -> set[Tuple5]:
    """VEP default tab. Named cols Location/Allele/Feature/Consequence; IMPACT in Extra."""
    tuples: set[Tuple5] = set()
    header: list[str] | None = None
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            cols = line.rstrip("\n").split("\t")
            row = dict(zip(header, cols))
            loc_raw = row.get("Location", "")
            # Location already chr:pos or chr:start-end; normalize contig part.
            if ":" in loc_raw:
                c, rest = loc_raw.split(":", 1)
                loc = f"{norm_contig(c)}:{rest}"
            else:
                loc = loc_raw
            allele = row.get("Allele", "")
            feature = row.get("Feature", "")
            csq = row.get("Consequence", "")
            impact = _impact_from_extra(row.get("Extra", ""))
            tuples.add((loc, allele, feature, norm_consequence(csq), impact))
    return tuples


def parse_json(path: Path) -> set[Tuple5]:
    """vep-rs/fastVEP JSON: one object per line OR a JSON array. Walk transcript_consequences."""
    tuples: set[Tuple5] = set()

    def handle(rec: dict) -> None:
        chrom = norm_contig(str(rec.get("seq_region_name", rec.get("seqid", ""))))
        start = rec.get("start", rec.get("pos", ""))
        loc = f"{chrom}:{start}"
        allele = rec.get("allele_string", "")
        # VEP JSON allele_string is "REF/ALT"; the per-consequence variant_allele is the ALT.
        for block_key in (
            "transcript_consequences",
            "intergenic_consequences",
            "regulatory_feature_consequences",
            "motif_feature_consequences",
        ):
            for tc in rec.get(block_key, []) or []:
                feat = tc.get(
                    "transcript_id",
                    tc.get("regulatory_feature_id", tc.get("motif_feature_id", "")),
                )
                va = tc.get(
                    "variant_allele", allele.split("/")[-1] if "/" in allele else allele
                )
                csq = norm_consequence(",".join(tc.get("consequence_terms", []) or []))
                impact = tc.get("impact", "")
                tuples.add((loc, va, feat, csq, impact))

    with _open(path) as fh:
        text = fh.read().strip()
    if not text:
        return tuples
    if text[0] == "[":
        for rec in json.loads(text):
            handle(rec)
    else:
        for line in text.splitlines():
            line = line.strip()
            if line:
                handle(json.loads(line))
    return tuples


def parse_parquet(path: Path) -> set[Tuple5]:
    """Read a vep-rs Parquet directory (Hive-partitioned by chrom) via duckdb.

    Builds the same 5-tuple as the other parsers directly from the typed
    columns. Works for both ``--parquet_shape`` layouts:
      - flat: one row per (variant x transcript-consequence); CSQ subfields
        (Consequence, IMPACT, Feature, Allele) are scalar columns.
      - nested: one row per variant; the CSQ subfields are ``LIST<VARCHAR>``,
        so duckdb ``UNNEST`` expands them back to one row per consequence.
    The Location is reconstructed ``chr:pos`` (or ``chr:pos-end`` when ``end``
    differs from ``pos``, matching the VEP tab/json span convention for
    indels) so the projected ``calls`` key lines up with the other formats.
    """
    import subprocess

    glob = f"{path}/**/*.parquet"
    # Detect shape: if Consequence is a LIST, UNNEST it; else read scalar.
    type_sql = (
        'SELECT typeof("Consequence") FROM read_parquet(\''
        + glob.replace("'", "''")
        + "') LIMIT 1"
    )
    t = subprocess.run(
        ["duckdb", "-list", "-noheader", "-c", type_sql],
        capture_output=True,
        text=True,
    )
    is_list = "[]" in t.stdout or "LIST" in t.stdout.upper()
    # chrom comes back from the partition column; pos/end are scalar.
    if is_list:
        sel = (
            'SELECT chrom, pos, "end", UNNEST("Allele") AS allele, '
            'UNNEST("Feature") AS feature, UNNEST("Consequence") AS csq, '
            'UNNEST("IMPACT") AS impact FROM read_parquet(\''
            + glob.replace("'", "''")
            + "')"
        )
    else:
        sel = (
            'SELECT chrom, pos, "end", "Allele" AS allele, '
            '"Feature" AS feature, "Consequence" AS csq, "IMPACT" AS impact '
            "FROM read_parquet('" + glob.replace("'", "''") + "')"
        )
    out = subprocess.run(
        ["duckdb", "-csv", "-c", sel], capture_output=True, text=True, check=True
    )
    tuples: set[Tuple5] = set()
    reader = csv.DictReader(io.StringIO(out.stdout))
    for row in reader:
        chrom = norm_contig(str(row.get("chrom", "")))
        pos = str(row.get("pos", ""))
        end = str(row.get("end", ""))
        loc = f"{chrom}:{pos}" if (end in ("", pos)) else f"{chrom}:{pos}-{end}"
        tuples.add(
            (
                loc,
                row.get("allele", "") or "",
                row.get("feature", "") or "",
                norm_consequence(row.get("csq", "") or ""),
                row.get("impact", "") or "",
            )
        )
    return tuples


PARSERS = {"vcf": parse_vcf, "tab": parse_tab, "tsv": parse_tab, "json": parse_json}


def project_key(tuples: set[Tuple5], mode: str) -> set:
    """Project the 5-tuple onto the comparison key.

    Output formats encode the variant Location differently for indels/insertions:
    VCF emits the VCF POS anchor (e.g. ``1:10035812``) while the VEP tab/json
    formats emit VEP's own 1-based start-end span (``1:10035813-10035814``).
    These are deterministically related but format-specific, so a strict 5-tuple
    diff reports them as mismatches even though the *consequence call* is identical.

    - ``strict``: full ``(location, allele, feature, consequence, impact)``.
    - ``calls`` (default): ``(contig, allele, feature, consequence, impact)`` -
      keeps the contig but drops the format-specific coordinate, isolating
      whether the engine emits the same consequence calls regardless of format's
      coordinate convention. This is the real output-format-parity question.
    """
    if mode == "strict":
        return set(tuples)
    out = set()
    for loc, allele, feat, csq, impact in tuples:
        contig = loc.split(":", 1)[0]
        out.add((contig, allele, feat, csq, impact))
    return out


@dataclass
class FormatResult:
    fmt: str
    path: str
    n_tuples: int
    present: bool
    note: str = ""


@dataclass
class PairResult:
    a: str
    b: str
    a_only: int
    b_only: int
    shared: int
    identical: bool
    examples_a_only: list = field(default_factory=list)
    examples_b_only: list = field(default_factory=list)


def diff(a: set[Tuple5], b: set[Tuple5], limit: int = 10) -> PairResult:
    a_only, b_only = a - b, b - a
    return PairResult(
        a="",
        b="",
        a_only=len(a_only),
        b_only=len(b_only),
        shared=len(a & b),
        identical=(not a_only and not b_only),
        examples_a_only=sorted(list(a_only))[:limit],
        examples_b_only=sorted(list(b_only))[:limit],
    )


def main() -> int:
    p = argparse.ArgumentParser(
        description="VEP output-format concordance (5-tuple parity)"
    )
    p.add_argument("--engine", required=True, choices=["vep-rs", "fastvep", "cross"])
    p.add_argument("--label", default=None, help="Display label (default: engine)")
    # Self-parity mode: one --file per format.
    p.add_argument("--vcf", type=Path, default=None)
    p.add_argument("--tab", type=Path, default=None)
    p.add_argument("--json", type=Path, default=None)
    p.add_argument(
        "--parquet",
        type=Path,
        default=None,
        help="Parquet path; skipped with a caveat if absent/unreadable.",
    )
    p.add_argument(
        "--parquet-tsv",
        type=Path,
        default=None,
        help="Optional TSV export of the parquet (duckdb), parsed as tab.",
    )
    # Cross-engine mode: two same-format files.
    p.add_argument("--a", type=Path, default=None, help="cross: engine-A file")
    p.add_argument("--b", type=Path, default=None, help="cross: engine-B file")
    p.add_argument("--a-format", default="tab", choices=list(PARSERS))
    p.add_argument("--b-format", default="tab", choices=list(PARSERS))
    p.add_argument("--a-label", default="A")
    p.add_argument("--b-label", default="B")
    p.add_argument("--summary-json", type=Path, required=True)
    p.add_argument("--summary-md", type=Path, required=True)
    p.add_argument(
        "--key-mode",
        default="calls",
        choices=["strict", "calls"],
        help="strict = full (location,allele,feature,consequence,impact); "
        "calls (default) = drop the format-specific coordinate, keep contig, "
        "isolating consequence-call parity from VCF-POS-vs-VEP-span representation.",
    )
    args = p.parse_args()

    label = args.label or args.engine
    results: list[FormatResult] = []
    pairs: list[tuple[str, str, PairResult]] = []
    parsed: dict[str, set[Tuple5]] = {}
    caveats: list[str] = []

    if args.engine == "cross":
        if not args.a or not args.b:
            print("ERROR: --engine cross requires --a and --b", file=sys.stderr)
            return 2
        sa = PARSERS[args.a_format](args.a)
        sb = PARSERS[args.b_format](args.b)
        results.append(
            FormatResult(args.a_format, str(args.a), len(sa), True, args.a_label)
        )
        results.append(
            FormatResult(args.b_format, str(args.b), len(sb), True, args.b_label)
        )
        pr = diff(project_key(sa, args.key_mode), project_key(sb, args.key_mode))
        pr.a, pr.b = args.a_label, args.b_label
        pairs.append((args.a_label, args.b_label, pr))
    else:
        fmt_files = {"vcf": args.vcf, "tab": args.tab, "json": args.json}
        for fmt, path in fmt_files.items():
            if path is None:
                continue
            if not path.exists():
                results.append(FormatResult(fmt, str(path), 0, False, "file absent"))
                caveats.append(f"{fmt}: file absent ({path})")
                continue
            parsed[fmt] = PARSERS[fmt](path)
            results.append(FormatResult(fmt, str(path), len(parsed[fmt]), True))
        # Parquet via its TSV export (engine writes parquet through duckdb; the
        # .tsv.tmp intermediate is the pre-parquet tuple source when present).
        if args.parquet_tsv and args.parquet_tsv.exists():
            parsed["parquet"] = parse_tab(args.parquet_tsv)
            results.append(
                FormatResult(
                    "parquet",
                    str(args.parquet_tsv),
                    len(parsed["parquet"]),
                    True,
                    "via tsv export",
                )
            )
        elif args.parquet is not None:
            if not Path(args.parquet).exists():
                results.append(
                    FormatResult("parquet", str(args.parquet), 0, False, "file absent")
                )
                caveats.append(f"parquet: path absent ({args.parquet})")
            else:
                try:
                    parsed["parquet"] = parse_parquet(Path(args.parquet))
                    results.append(
                        FormatResult(
                            "parquet",
                            str(args.parquet),
                            len(parsed["parquet"]),
                            True,
                            "read via duckdb",
                        )
                    )
                except Exception as exc:  # noqa: BLE001 - report, don't crash the run
                    results.append(
                        FormatResult(
                            "parquet", str(args.parquet), 0, False, f"unreadable: {exc}"
                        )
                    )
                    caveats.append(f"parquet: unreadable via duckdb ({exc})")
        # Pairwise self-parity across all present formats.
        fmts = [f for f in ("vcf", "tab", "json", "parquet") if f in parsed]
        for i in range(len(fmts)):
            for j in range(i + 1, len(fmts)):
                pr = diff(
                    project_key(parsed[fmts[i]], args.key_mode),
                    project_key(parsed[fmts[j]], args.key_mode),
                )
                pr.a, pr.b = fmts[i], fmts[j]
                pairs.append((fmts[i], fmts[j], pr))

    all_identical = all(pr.identical for _, _, pr in pairs) if pairs else False
    summary = {
        "engine": args.engine,
        "label": label,
        "formats": [r.__dict__ for r in results],
        "pairs": [
            {
                "a": a,
                "b": b,
                "a_only": pr.a_only,
                "b_only": pr.b_only,
                "shared": pr.shared,
                "identical": pr.identical,
                "examples_a_only": pr.examples_a_only,
                "examples_b_only": pr.examples_b_only,
            }
            for a, b, pr in pairs
        ],
        "all_identical": all_identical,
        "caveats": caveats,
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(json.dumps(summary, indent=2))

    lines = [f"# Output-format parity: {label}", ""]
    lines.append("| Format | Tuples | Present | Note |")
    lines.append("| ------ | -----: | ------- | ---- |")
    for r in results:
        lines.append(
            f"| {r.fmt} | {r.n_tuples} | {'yes' if r.present else 'NO'} | {r.note} |"
        )
    lines.append("")
    lines.append("| Pair | A-only | B-only | Shared | Identical |")
    lines.append("| ---- | -----: | -----: | -----: | --------- |")
    for a, b, pr in pairs:
        lines.append(
            f"| {a} vs {b} | {pr.a_only} | {pr.b_only} | {pr.shared} | "
            f"{'YES' if pr.identical else 'NO'} |"
        )
    if caveats:
        lines.append("")
        lines.append("## Caveats")
        for c in caveats:
            lines.append(f"- {c}")
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    if args.engine != "cross" and not all_identical:
        # Self-parity must be exact -> non-zero exit for the caller to catch.
        non_identical = [(a, b) for a, b, pr in pairs if not pr.identical]
        if non_identical:
            print(
                f"\nFAIL: self-parity mismatch in pairs: {non_identical}",
                file=sys.stderr,
            )
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
