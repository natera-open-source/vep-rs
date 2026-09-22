#!/usr/bin/env python3
"""Build a golden corpus: exemplar input records for every consequence-set combination a
VEP reference output contains, keyed by Ensembl release and assembly.

The corpus is the fixture behind `crates/vep-cli/tests/golden.rs` and
`crates/vep-cli/tests/format_parity.rs`: a small VCF whose records, annotated by VEP and by
vep-rs against the same transcript cache, must yield the same rows in every output format.

Subcommands:

  select    Stream one or more VEP default-format reference outputs, enumerate the distinct
            consequence sets they contain, choose up to K exemplar records per set
            (stratified by variant class and transcript strand, preferring records that
            exercise transcript flags, coordinate ranges and unmappable ends), resolve each
            exemplar back to its input VCF record, and write `variants.vcf` plus
            `manifest.json`.
  classify  Compare a VEP default-format output for `variants.vcf` with a vep-rs output for
            the same file and record, per exemplar row VEP emits that vep-rs does not, the
            documented divergence class that explains it (or `unexplained_residual`), so the
            golden test asserts the documented difference instead of equality.
  check     Recompute the observed combination space from the reference outputs and fail
            when a combination has no exemplar in the manifest.

Row keys follow `compare_vep_outputs.py`: (Location, Allele, Feature, Feature_type,
consequence set), consequence terms sorted and comma-joined. A record is identified by
VEP's (Uploaded_variation, Location) pair, which every row of one input record shares.

Selection is deterministic under `--seed`: reservoir sampling bounds the candidates kept per
combination and the same seed reproduces the same corpus from the same references.
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import re
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import IO, Iterable, Iterator

GENERATOR_VERSION = "2"
RESERVOIR = 64
SO_CLASS_WORDS = {
    "deletion", "duplication", "inversion", "insertion", "copy_number_variation",
    "tandem_duplication", "mobile_element_insertion", "tandem_repeat", "indel",
    "sequence_alteration", "chromosome_breakpoint", "SNV", "substitution",
    "copy_number_gain", "copy_number_loss", "translocation",
}
DEFAULT_COLUMNS = (
    "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "Extra",
)


def opener(path: Path):
    return gzip.open if str(path).endswith(".gz") else open


def normalize_consequence_set(raw: str) -> str:
    return ",".join(sorted({t for t in raw.strip().split(",") if t}))


def parse_extra(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not raw or raw == "-":
        return out
    for part in raw.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        elif part:
            out[part] = ""
    return out


def variant_class(location: str, allele: str) -> str:
    """Infer VEP's variant shape from the Location and Allele columns of one row."""
    if allele in SO_CLASS_WORDS or (allele and not re.fullmatch(r"[ACGTNacgtn*-]+", allele)):
        return "symbolic"
    span = location.split(":", 1)[1] if ":" in location else location
    if allele == "-":
        return "deletion"
    if "-" in span:
        a, b = span.split("-", 1)
        if a.isdigit() and b.isdigit() and int(b) == int(a) + 1:
            return "insertion"
        return "mnv_or_complex"
    return "snv" if len(allele) == 1 else "mnv_or_complex"


@dataclass
class Row:
    suite: str
    uploaded_variation: str
    location: str
    allele: str
    feature: str
    feature_type: str
    consequence_set: str
    strand: str
    flags: str
    positions: tuple[str, str, str]

    @property
    def record_key(self) -> tuple[str, str]:
        return (self.uploaded_variation, self.location)

    @property
    def klass(self) -> str:
        return variant_class(self.location, self.allele)

    def preference(self) -> int:
        score = 0
        if self.flags:
            score += 4
        if any("?" in p for p in self.positions):
            score += 3
        if any("-" in p and p != "-" for p in self.positions):
            score += 2
        return score

    @property
    def span(self) -> int:
        """Genomic width of the Location, in bases (0 for a single position)."""
        coords = self.location.split(":", 1)[-1]
        if "-" not in coords:
            return 0
        a, b = coords.split("-", 1)
        try:
            return abs(int(b) - int(a))
        except ValueError:
            return 0


def iter_rows(path: Path, suite: str) -> Iterator[Row]:
    with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            while len(f) < 14:
                f.append("-")
            extra = parse_extra(f[13])
            yield Row(
                suite=suite,
                uploaded_variation=f[0],
                location=f[1],
                allele=f[2],
                feature=f[4],
                feature_type=f[5],
                consequence_set=normalize_consequence_set(f[6]),
                strand=extra.get("STRAND", ""),
                flags=extra.get("FLAGS", ""),
                positions=(f[7], f[8], f[9]),
            )


def parse_suite_arg(arg: str) -> tuple[str, Path, Path]:
    """`sid=reference.txt=input.vcf[.gz]`."""
    parts = arg.split("=")
    if len(parts) != 3:
        raise SystemExit(f"ERROR: [build_golden_corpus] --suite wants sid=reference=input, got {arg!r}")
    return parts[0], Path(parts[1]), Path(parts[2])


# ----------------------------------------------------------------------------- select
def reservoir_pass(suites: list[tuple[str, Path, Path]], seed: int) -> tuple[dict[str, list[Row]], Counter[str], dict[str, Counter[str]]]:
    """One streaming pass: per consequence set, a bounded uniform sample of rows plus counts."""
    rng = random.Random(seed)
    reservoir: dict[str, list[Row]] = defaultdict(list)
    seen: Counter[str] = Counter()
    per_suite: dict[str, Counter[str]] = defaultdict(Counter)
    for sid, ref, _inp in suites:
        for row in iter_rows(ref, sid):
            c = row.consequence_set
            seen[c] += 1
            per_suite[sid][c] += 1
            bucket = reservoir[c]
            if len(bucket) < RESERVOIR:
                bucket.append(row)
            else:
                j = rng.randrange(seen[c])
                if j < RESERVOIR:
                    bucket[j] = row
    return reservoir, seen, per_suite


def choose_exemplars(reservoir: dict[str, list[Row]], k: int, rng: random.Random, max_span: int = 0) -> dict[str, list[Row]]:
    """Up to K rows per combination, distinct (class, strand) first, higher preference first.

    Rows wider than `max_span` bases (when > 0) come after every narrower row, so a
    combination takes a wide structural variant only when nothing narrower shows it:
    each wide record pulls every transcript of its span into the fixture cache.
    """
    chosen: dict[str, list[Row]] = {}
    for combo in sorted(reservoir):
        cands = list(reservoir[combo])
        rng.shuffle(cands)
        cands.sort(key=lambda r: (1 if max_span and r.span > max_span else 0, -r.preference()))
        picked: list[Row] = []
        strata_seen: set[tuple[str, str]] = set()
        for r in cands:
            stratum = (r.klass, r.strand)
            if stratum in strata_seen:
                continue
            picked.append(r)
            strata_seen.add(stratum)
            if len(picked) == k:
                break
        for r in cands:
            if len(picked) == k:
                break
            if r not in picked:
                picked.append(r)
        chosen[combo] = picked
    return chosen


def vep_auto_names(chrom: str, pos: int, ref: str, alts: list[str]) -> list[str]:
    """Reproduce VEP's auto-generated Uploaded_variation for an ID-less VCF record.

    The alleles are always the RAW REF and ALTs joined by `/` (VEP's
    `nontrimmed_allele_string`). The start is the variant's after parsing: unchanged
    when no allele differs in length from REF; for a bi-allelic indel the first shared
    base is chopped and the alleles are then trimmed from both ends
    (`trim_sequences`); for a multi-allelic indel the first base is chopped when every
    allele (ignoring `*`) shares it. Symbolic and breakend alleles: `chr_start_<ALT...>`
    with start POS or POS+1 (both are returned; which one VEP used depends on whether
    the REF base was padding).
    """
    if any(a.startswith("<") or "[" in a or "]" in a for a in alts):
        joined = "/".join(alts)
        return [f"{chrom}_{pos}_{joined}", f"{chrom}_{pos + 1}_{joined}"]
    raw = "/".join([ref] + alts)
    if not any(len(a) != len(ref) for a in alts):
        return [f"{chrom}_{pos}_{raw}"]
    if len(alts) == 1:
        r, a, start = ref, alts[0], pos
        if r[:1] == a[:1]:
            r, a, start = r[1:] or "-", a[1:] or "-", start + 1
        while r and a and r[0] == a[0]:
            r, a, start = r[1:], a[1:], start + 1
        while r and a and r[-1] == a[-1]:
            r, a = r[:-1], a[:-1]
        return [f"{chrom}_{start}_{raw}"]
    firsts = {x[:1] for x in [ref] + alts if "*" not in x}
    start = pos + 1 if len(firsts) == 1 else pos
    return [f"{chrom}_{start}_{raw}"]


MULTI_ALLELIC_STRATA = ("multi_allelic_snv", "multi_allelic_indel", "multi_allelic_indel_trimmable")


def multi_allelic_stratum(ref: str, alts: list[str]) -> str:
    """Which multi-allelic stratum a record belongs to: SNV-only alleles; alleles of
    differing length; or such alleles that VEP's bi-allelic minimisation WOULD trim
    further after the shared first base is chopped (a shared prefix or suffix left
    between REF and some ALT), the shape on which a per-allele minimiser and VEP's
    untouched multi-allelic record disagree on positions and distances."""
    if not any(len(a) != len(ref) for a in alts):
        return "multi_allelic_snv"
    alleles = [ref] + alts
    if len({x[:1] for x in alleles}) == 1:
        alleles = [x[1:] or "-" for x in alleles]
    r = alleles[0]
    for a in alleles[1:]:
        if r == "-" or a == "-" or len(r) == len(a):
            continue
        if r[0] == a[0] or r[-1] == a[-1]:
            return "multi_allelic_indel_trimmable"
    return "multi_allelic_indel"


def resolve_records(
    inputs: Iterable[tuple[str, Path]],
    wanted: dict[str, set[str]],
    multi_per_suite: int = 0,
    seed: int = 1,
) -> tuple[dict[str, dict[tuple[str, str], list[str]]], list[str], dict[str, dict[str, list[str]]]]:
    """Scan each input VCF once; return per suite the record lines whose Uploaded_variation
    (ID, or VEP's auto-name) is wanted, the VCF header of each input, and a uniform sample
    of `multi_per_suite` multi-allelic sequence records per suite and stratum (SNV-only
    alleles; alleles of differing length), which no consequence set can select for and
    which exercise VEP's treatment of a multi-allelic record as one variant."""
    rng = random.Random(seed)
    found: dict[str, dict[tuple[str, str], list[str]]] = defaultdict(lambda: defaultdict(list))
    multi: dict[str, dict[str, list[str]]] = defaultdict(lambda: {k: [] for k in MULTI_ALLELIC_STRATA})
    multi_seen: Counter[tuple[str, str]] = Counter()
    headers: list[str] = []
    for sid, path in inputs:
        want = wanted.get(sid, set())
        if not want and not multi_per_suite:
            continue
        with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    if not headers:
                        headers.append(line.rstrip("\n"))
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
                alts = alt.split(",")
                record = "\t".join(f[:8])
                if multi_per_suite and len(alts) > 1 and not any(a.startswith("<") or "[" in a or "]" in a or a == "*" for a in alts):
                    stratum = multi_allelic_stratum(ref, alts)
                    multi_seen[(sid, stratum)] += 1
                    bucket = multi[sid][stratum]
                    if len(bucket) < multi_per_suite:
                        bucket.append(record)
                    else:
                        j = rng.randrange(multi_seen[(sid, stratum)])
                        if j < multi_per_suite:
                            bucket[j] = record
                if not want:
                    continue
                names: list[str] = []
                if vid != ".":
                    names.append(vid.split(";")[0])
                names.extend(vep_auto_names(chrom, pos, ref, alts))
                for n in names:
                    if n in want:
                        found[sid][(n, chrom)].append(record)
                        break
    return found, headers, multi


def cmd_select(args: argparse.Namespace) -> int:
    suites = [parse_suite_arg(s) for s in args.suite]
    rng = random.Random(args.seed)
    reservoir, seen, per_suite = reservoir_pass(suites, args.seed)
    chosen = choose_exemplars(reservoir, args.k, rng, args.max_span)
    wanted: dict[str, set[str]] = defaultdict(set)
    for combo, rows in chosen.items():
        for r in rows:
            wanted[r.suite].add(r.uploaded_variation)
    found, headers, multi = resolve_records(
        [(sid, inp) for sid, _ref, inp in suites], wanted, args.multi_allelic_per_suite, args.seed
    )
    records: dict[str, dict] = {}
    missing: list[str] = []
    for combo, rows in chosen.items():
        for r in rows:
            chrom = r.location.split(":", 1)[0]
            lines = found.get(r.suite, {}).get((r.uploaded_variation, chrom), [])
            if not lines:
                missing.append(f"{r.suite}\t{r.uploaded_variation}\t{r.location}")
                continue
            for line in lines:
                key = f"{r.suite}\t{line}"
                rec = records.setdefault(key, {"suite": r.suite, "vcf": line, "combinations": [], "exemplar_rows": []})
                if combo not in rec["combinations"]:
                    rec["combinations"].append(combo)
                rec["exemplar_rows"].append({
                    "uploaded_variation": r.uploaded_variation, "location": r.location, "allele": r.allele,
                    "feature": r.feature, "variant_class": r.klass, "strand": r.strand,
                })
    for sid in sorted(multi):
        for stratum, lines in sorted(multi[sid].items()):
            for line in lines:
                key = f"{sid}\t{line}"
                rec = records.setdefault(key, {"suite": sid, "vcf": line, "combinations": [], "exemplar_rows": []})
                rec.setdefault("strata", []).append(stratum)
    if len(records) > args.max_records:
        # Keep every record that is the sole exemplar of some combination, then fill by
        # preference until the cap; deterministic because the record order is sorted.
        by_combo: dict[str, list[str]] = defaultdict(list)
        for key, rec in records.items():
            for c in rec["combinations"]:
                by_combo[c].append(key)
        keep = {keys[0] for keys in by_combo.values()}
        for key in sorted(records):
            if len(keep) >= args.max_records:
                break
            keep.add(key)
        records = {k: v for k, v in records.items() if k in keep}

    def sort_key(item):
        f = item[1]["vcf"].split("\t")
        chrom = f[0]
        chrom_rank = (0, int(chrom)) if chrom.isdigit() else (1, {"X": 1, "Y": 2, "MT": 3}.get(chrom, 9))
        return (chrom_rank, int(f[1]), f[3], f[4])

    ordered = sorted(records.items(), key=sort_key)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    header = headers[0] if headers else "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    with (out / "variants.vcf").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(header.split("\tFORMAT")[0] + "\n")
        for _key, rec in ordered:
            fh.write(rec["vcf"] + "\n")
    covered = {c for _k, rec in ordered for c in rec["combinations"]}
    manifest = {
        "release": args.release,
        "assembly": args.assembly,
        "generator_version": GENERATOR_VERSION,
        "seed": args.seed,
        "exemplars_per_combination": args.k,
        "max_span": args.max_span,
        # File names only: the directory each suite was staged in belongs to the build
        # host, not to the corpus, and the tests read neither field.
        "suites": {sid: {"reference": Path(ref).name, "input": Path(inp).name, "rows": sum(per_suite[sid].values()), "combinations": len(per_suite[sid])} for sid, ref, inp in suites},
        "combinations": {c: {"observed_rows": seen[c], "exemplar_records": sum(1 for _k, rec in ordered if c in rec["combinations"])} for c in sorted(seen)},
        "records": [
            {"suite": rec["suite"], "vcf": rec["vcf"], "combinations": sorted(rec["combinations"]), **({"strata": rec["strata"]} if rec.get("strata") else {})}
            for _k, rec in ordered
        ],
        "multi_allelic_per_suite": args.multi_allelic_per_suite,
        "uncovered_combinations": sorted(set(seen) - covered),
        "unresolved_exemplars": missing,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    print(f"combinations observed {len(seen)}, covered {len(covered)}, records {len(ordered)}, unresolved exemplars {len(missing)}")
    return 0 if not missing and not manifest["uncovered_combinations"] else 1


# --------------------------------------------------------------------------- classify
def load_tuples(path: Path) -> tuple[dict[tuple[str, str, str, str], set[str]], dict[tuple[str, str, str, str], set[str]]]:
    """(Location, Allele, Feature, Feature_type) -> consequence-set strings, and -> Uploaded_variation names."""
    out: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    names: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for r in iter_rows(path, "x"):
        key = (r.location, r.allele, r.feature, r.feature_type)
        out[key].add(r.consequence_set)
        names[key].add(r.uploaded_variation)
    return out, names


def record_name_index(records: list[dict]) -> dict[str, list[int]]:
    """Uploaded_variation (the ID, else VEP's auto-name) -> ordinals of the corpus records it names."""
    index: dict[str, list[int]] = defaultdict(list)
    for i, rec in enumerate(records):
        f = rec["vcf"].split("\t")
        chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
        names = [vid.split(";")[0]] if vid != "." else vep_auto_names(chrom, pos, ref, alt.split(","))
        for n in names:
            index[n].append(i)
    return index


def classify_pair(perl_csq: str, rust_csq: str, allele: str) -> str | None:
    """First matching rule of the key comparator's registry, by taxonomy class."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
    import compare_vep_outputs as cvo  # noqa: E402

    for rule in cvo.EXCLUSION_REGISTRY:
        if rule.matches(perl_csq, rust_csq, allele):
            return rule.taxonomy_class
    return None


def cmd_classify(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    perl, perl_names = load_tuples(Path(args.vep_default))
    rust, rust_names = load_tuples(Path(args.vep_rs_output))
    by_name = record_name_index(manifest["records"])

    def records_for(key: tuple[str, str, str, str]) -> list[int]:
        found: set[int] = set()
        for n in perl_names.get(key, set()) | rust_names.get(key, set()):
            found.update(by_name.get(n, []))
        return sorted(found)

    divergences: list[dict] = []
    for key, perl_sets in sorted(perl.items()):
        rust_sets = rust.get(key, set())
        for cs in sorted(perl_sets - rust_sets):
            klass = None
            for rs in sorted(rust_sets):
                klass = classify_pair(cs, rs, key[1])
                if klass:
                    break
            if klass is None and not rust_sets:
                klass = "vep_only_transcript" if variant_class(key[0], key[1]) == "symbolic" else None
            divergences.append({
                "location": key[0], "allele": key[1], "feature": key[2], "feature_type": key[3],
                "record_indices": records_for(key),
                "vep_consequence_set": cs, "vep_rs_consequence_sets": sorted(rust_sets),
                "expected_divergence": klass or "unexplained_residual",
            })
    extra_rust = [
        {"location": k[0], "allele": k[1], "feature": k[2], "feature_type": k[3],
         "record_indices": records_for(k), "vep_rs_consequence_sets": sorted(v)}
        for k, v in sorted(rust.items()) if k not in perl
    ]
    manifest["divergences"] = divergences
    manifest["vep_rs_only_tuples"] = extra_rust
    manifest["divergence_summary"] = dict(Counter(d["expected_divergence"] for d in divergences))
    (corpus / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    print(f"divergent VEP rows {len(divergences)} ({manifest['divergence_summary']}); vep-rs-only tuples {len(extra_rust)}")
    return 0


# ------------------------------------------------------------------------------ check
def cmd_check(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    suites = [parse_suite_arg(s) for s in args.suite]
    observed: set[str] = set()
    for sid, ref, _inp in suites:
        for r in iter_rows(ref, sid):
            observed.add(r.consequence_set)
    covered = {c for rec in manifest["records"] for c in rec["combinations"]}
    missing = sorted(observed - covered)
    print(f"observed {len(observed)}, covered by corpus {len(covered & observed)}, missing {len(missing)}")
    for m in missing[:50]:
        print("  MISSING", m)
    return 1 if missing else 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--release", required=True)
    s.add_argument("--assembly", required=True)
    s.add_argument("--suite", action="append", required=True, help="sid=reference.txt=input.vcf[.gz]")
    s.add_argument("--out", required=True)
    s.add_argument("--k", type=int, default=3)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--max-records", type=int, default=1500)
    s.add_argument("--multi-allelic-per-suite", type=int, default=4,
                   help="multi-allelic sequence records sampled per suite and stratum (SNV-only, indel)")
    s.add_argument("--max-span", type=int, default=1_000_000,
                   help="prefer exemplars no wider than this many bases; 0 disables the preference")
    s.set_defaults(fn=cmd_select)
    c = sub.add_parser("classify")
    c.add_argument("--corpus", required=True)
    c.add_argument("--vep-default", required=True, help="VEP default-format output for variants.vcf")
    c.add_argument("--vep-rs-output", required=True, help="vep-rs default-format output for variants.vcf")
    c.set_defaults(fn=cmd_classify)
    k = sub.add_parser("check")
    k.add_argument("--corpus", required=True)
    k.add_argument("--suite", action="append", required=True)
    k.set_defaults(fn=cmd_check)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
