#!/usr/bin/env python3
"""Render the GitHub release body for a vep-rs version from the repository's own files.

The page, every header `##`, in this order: the hand-written opener, then `## Highlights` and
`## Upgrade notes` when `.github/release-notes/vX.Y.Z.md` carries them (its lines are the only
hand-written ones); the version's CHANGELOG.md section with each `###` kind header raised to `##`;
`## Concordance and wall time` from the release and population provenance records; `## Install`
from SHA256SUMS, the build and the workspace rust-version; `## Cite` from CITATION.cff; and the
full-changelog footer. A page runs about 60 lines; one over 80 is refused.

usage: render_release_notes.py --repo DIR --version X.Y.Z --previous vX.Y.Z --commit SHA
                               --toolchain 1.97.1 --sha256sums FILE [--notes FILE]
                               [--changelog-section NAME] [--date YYYY-MM-DD] [--out FILE]

A dry run before the CHANGELOG has the version's dated section passes
`--changelog-section Unreleased --date <today>`.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
from reference_release import release_from_provenance  # noqa: E402  (one reading of a reference's release for every script)

REPO_URL = "https://github.com/natera-open-source/vep-rs"
CONCEPT_DOI = "10.5281/zenodo.22837896"
LINE_CAP = 80
SCORED = ("s01", "s02", "s03", "s04", "s05", "s06", "s07", "s08")
TARGETS = {
    "x86_64-unknown-linux-gnu": (
        "Linux x86_64, glibc 2.39 or later",
        "x86-64-v3 (AVX2, BMI2, FMA: Intel Haswell 2013 or later, AMD Zen or later)",
    ),
    "aarch64-unknown-linux-gnu": (
        "Linux aarch64, glibc 2.39 or later",
        "Armv8.4-A with SVE (Arm Neoverse V1 or later: AWS Graviton3 and later, not Graviton2, Ampere Altra or Raspberry Pi)",
    ),
    "aarch64-apple-darwin": ("macOS 11 or later on Apple silicon", "any Apple M-series"),
}
SUITE_LABELS = {
    "s01": ("ClinVar full", "GRCh37"),
    "s02": ("ClinVar full", "GRCh38"),
    "s03": ("gnomAD v2.1.1 chr21", "GRCh37"),
    "s04": ("gnomAD v4.1 chr21", "GRCh38"),
    "s05": ("1KG Phase 3 chr21", "GRCh37"),
    "s06": ("1KG high-cov chr21", "GRCh38"),
    "s07": ("SV per-VCF (16 files)", "GRCh37"),
    "s08": ("SV per-VCF (16 files)", "GRCh38"),
}
CLINVAR = {"s01": "ClinVar GRCh37", "s02": "ClinVar GRCh38"}
SV = {"s07": "GRCh37", "s08": "GRCh38"}
FAMILIES = ("snv_indel", "sv")
ARCH_LABELS = {"arm64": "ARM", "x86_64": "x86"}
SV_SETS = (("clinvar_sv_", "ClinVar SV"), ("gnomad_sv_", "gnomAD SV"), ("1kg_sv_", "1000 Genomes SV"))
SYNTHETIC = re.compile(r"^\d{2}_")
SYNTHETIC_FILES = 13
LIST_ITEM = re.compile(r"^\s*(?:[-*+]|\d+\.)\s")
BULLET = re.compile(r"[-*] ")
SENTENCE_END = re.compile(r"[.!?][\"'`)\]]*(?=\s|$)")


def fail(message: str) -> None:
    sys.exit(f"ERROR: [render_release_notes] {message}")


def half_up(x: float, places: int) -> str:
    return str(Decimal(str(x)).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))


def f1_9(intersection, a, b) -> str:
    """F1 as 2I/(A+B) from the integer tuple counts at nine decimal places; an exact match reads 1.000000000,
    and so do two empty sides, the figure the records' own per-file `f1` carries for that case. Computed from
    the counts rather than copied from a record's rounded field, so the page never prints a coarser rounding
    than the counts support."""
    denominator = int(a) + int(b)
    if denominator == 0:
        return "1.000000000"
    q = Decimal(2 * int(intersection)) / Decimal(denominator)
    return str(q.quantize(Decimal("0.000000001"), rounding=ROUND_HALF_UP))


def pooled_f1(cells: list[tuple[int, int, int]]) -> str:
    """F1 of several cells from their summed (intersection, A, B) counts."""
    return f1_9(*(sum(c) for c in zip(*cells)))


def whole_seconds(x: Decimal) -> str:
    """A wall time as whole seconds, rounded half up with thousands separators; one too large for the Decimal
    context to round to a whole number is refused."""
    try:
        return f"{int(x.quantize(Decimal(1), rounding=ROUND_HALF_UP)):,}"
    except InvalidOperation:
        fail(f"a wall time of {x.normalize()} seconds is too large to round to whole seconds")


def changelog_section(text: str, name: str) -> tuple[str, str | None]:
    """Return (body, date) of the `## [name]` section, body without the heading."""
    pattern = re.compile(rf"^## \[{re.escape(name)}\](?: - (\d{{4}}-\d{{2}}-\d{{2}}))?[^\n]*\n", re.MULTILINE)
    m = pattern.search(text)
    if not m:
        fail(f"no '## [{name}]' section in CHANGELOG.md")
    rest = text[m.end() :]
    nxt = re.search(r"^## \[", rest, re.MULTILINE)
    body = rest[: nxt.start()] if nxt else rest
    body = re.sub(r"\n{3,}", "\n\n", body).strip("\n")
    return body, m.group(1)


def unwrap_bullets(text: str) -> str:
    """Join each list item's wrapped continuation lines (indented, not list items themselves, outside a
    fenced block) into its line, so the page carries one physical line per bullet in the changelog's words."""
    out: list[str] = []
    fenced = False
    for line in text.split("\n"):
        if line.lstrip().startswith(("```", "~~~")):
            fenced = not fenced
        elif (
            not fenced
            and out
            and line.startswith("  ")
            and line.strip()
            and not LIST_ITEM.match(line)
            and LIST_ITEM.match(out[-1])
        ):
            out[-1] = out[-1] + " " + line.strip()
            continue
        out.append(line)
    return "\n".join(out)


def kind_blocks(section: str) -> list[str]:
    """The page blocks of a CHANGELOG section: text before the first kind header as it stands, then each
    `### Kind` header as `## Kind` followed by its content; a kind with no content is dropped."""
    groups: list[tuple[str | None, list[str]]] = [(None, [])]
    for line in unwrap_bullets(section).split("\n"):
        m = re.match(r"^### (.+)$", line)
        if m:
            groups.append((m.group(1).strip(), []))
        else:
            groups[-1][1].append(line)
    blocks: list[str] = []
    for kind, lines in groups:
        content = "\n".join(lines).strip("\n")
        if not content:
            continue
        if kind is not None:
            blocks.append(f"## {kind}")
        blocks.append(content)
    return blocks


def hand_written(path: Path) -> dict[str, str]:
    """The maintainer's lines: the opener paragraph, then the optional `## Highlights` and
    `## Upgrade notes` sections, each bullet on one line; any other section is an error."""
    out = {"opener": "", "highlights": "", "upgrade": ""}
    if not path.exists():
        return out
    text = re.sub(r"^# .*\n", "", path.read_text(), count=1)
    parts = re.split(r"^## ", text, flags=re.MULTILINE)
    out["opener"] = parts[0].strip()
    for part in parts[1:]:
        title, _, body = part.partition("\n")
        key = {"Highlights": "highlights", "Upgrade notes": "upgrade"}.get(title.strip())
        if key is None:
            fail(f"unknown section '## {title.strip()}' in {path}")
        out[key] = unwrap_bullets(body.strip("\n")).strip()
    return out


def bullets_only(name: str, body: str, notes: Path) -> list[str]:
    lines = [line for line in body.splitlines() if line.strip()]
    if any(not BULLET.match(line) for line in lines):
        fail(f"'## {name}' in {notes} is bullets only")
    return lines


def check_hand_written(hand: dict[str, str], section: str, notes: Path) -> None:
    """The opener is one paragraph of at most two sentences and 50 words; Highlights, when present, is two
    to four bullets over a CHANGELOG section of at least three; Upgrade notes, when present, is one to five
    bullets with every `**Breaking:**` bullet before the first that is not."""
    opener = hand["opener"]
    if not opener:
        fail("the release-notes file needs an opener paragraph")
    if "\n\n" in opener:
        fail(f"the opener in {notes} is more than one paragraph")
    words = len(opener.split())
    if words > 50:
        fail(f"the opener in {notes} is {words} words; at most 50")
    sentences = len(SENTENCE_END.findall(opener))
    if sentences > 2:
        fail(f"the opener in {notes} is {sentences} sentences; at most two")
    if hand["highlights"]:
        lines = bullets_only("Highlights", hand["highlights"], notes)
        if not 2 <= len(lines) <= 4:
            fail(f"'## Highlights' in {notes} has {len(lines)} bullets; two to four")
        changes = sum(1 for line in unwrap_bullets(section).splitlines() if BULLET.match(line))
        if changes < 3:
            fail(f"'## Highlights' needs at least three CHANGELOG bullets; the section has {changes}")
    if hand["upgrade"]:
        lines = bullets_only("Upgrade notes", hand["upgrade"], notes)
        if not 1 <= len(lines) <= 5:
            fail(f"'## Upgrade notes' in {notes} has {len(lines)} bullets; one to five")
        breaking = [line[2:].startswith("**Breaking:**") for line in lines]
        first_plain = breaking.index(False) if False in breaking else len(breaking)
        if any(breaking[first_plain:]):
            fail(f"'## Upgrade notes' in {notes} lists a **Breaking:** bullet after one that is not")


def read_doi(cff: Path, version: str) -> str | None:
    """This version's own DOI, or None.

    Read from an `identifiers` entry whose description names `version X.Y.Z`, or from the
    top-level `doi` when the file's `version` is this version and the DOI is not the concept
    DOI; the concept DOI never counts as a version's own."""
    text = cff.read_text()
    for m in re.finditer(r"value: (10\.5281/zenodo\.\d+)\n\s+description: (.*)", text):
        doi, desc = m.group(1), m.group(2)
        vm = re.search(r"version (\d+\.\d+\.\d+)", desc)
        if vm and doi != CONCEPT_DOI and vm.group(1) == version:
            return doi
    top_doi = re.search(r"^doi: (10\.5281/zenodo\.\d+)$", text, re.MULTILINE)
    top_version = re.search(r"^version: (\d+\.\d+\.\d+)$", text, re.MULTILINE)
    if top_doi and top_version and top_doi.group(1) != CONCEPT_DOI and top_version.group(1) == version:
        return top_doi.group(1)
    return None


def latest_record(repo: Path, kind: str, version: str) -> Path:
    """The last docs/concordance-provenance/<date>-<kind>-v<version>.json by name, which is the latest date."""
    pattern = f"*-{kind}-v{version}.json"
    records = sorted((repo / "docs" / "concordance-provenance").glob(pattern))
    if not records:
        fail(f"no docs/concordance-provenance/{pattern}")
    return records[-1]


def six_dataset_shape(release: dict, conc: dict) -> bool:
    """Whether a release record is of the six-dataset shape: its ClinVar GRCh37 cell carries the adjusted
    counts, or its `perl_reference` block is keyed by suite family (`snv_indel`, `sv`). The version number
    never decides, so a record of the earlier shape renders the earlier page whatever its version."""
    block = release.get("perl_reference")
    by_family = isinstance(block, dict) and any(isinstance(block.get(f), dict) for f in FAMILIES)
    s01 = conc.get("s01")
    return by_family or (isinstance(s01, dict) and "adj_intersection" in s01)


def load_release_record(repo: Path, version: str) -> tuple[dict, dict, dict, str]:
    """The release record, its per-suite concordance, its whole-genome concordance and the record directory's
    path; the record and its two blocks are refused by name when not objects, and the whole-genome block when
    it lacks `per_dataset` or `pooled`. A record of the six-dataset shape needs the two ClinVar cells and the
    two structural-variant cells; any other needs the eight scored suites."""
    record_path = latest_record(repo, "release", version)
    record = json.loads(record_path.read_text())
    if not isinstance(record, dict):
        fail(f"{record_path.name} is not an object")
    conc = record.get("concordance")
    if not isinstance(conc, dict):
        fail("the release record's concordance is not an object")
    needed = (*CLINVAR, *SV) if six_dataset_shape(record, conc) else SCORED
    missing = [s for s in needed if not isinstance(conc.get(s), dict)]
    if missing:
        fail(f"{record_path.name} lacks concordance rows for {missing}")
    pop = record.get("population_concordance")
    if pop is None:
        fail(f"{record_path.name} lacks population_concordance")
    if not isinstance(pop, dict):
        fail("the release record's population_concordance is not an object")
    per_dataset, _ = required(pop, ("per_dataset", "pooled"), "population_concordance")
    if not isinstance(per_dataset, list):
        fail("the release record's population_concordance.per_dataset is not a list")
    return record, conc, pop, record_path.relative_to(repo).parent.as_posix()


def load_population_record(repo: Path, version: str) -> dict:
    """The release's population record: the whole-genome wall time per chromosome and architecture; one that is
    not an object is refused by name."""
    path = latest_record(repo, "population", version)
    population = json.loads(path.read_text())
    if not isinstance(population, dict):
        fail(f"{path.name} is not an object")
    return population


def comparator_version(population: dict, pop: dict, release: dict | None = None) -> str:
    """The Ensembl VEP release the whole genome's concordance was scored against: the `perl_reference`
    block of the release record (top level or under `population_concordance`) or of the population
    record when one is present, read as the comparators read a reference's provenance (the image tag,
    else the archived prefix, else the cache version mapped to its release, so a block naming cache
    116 alone renders 116.2); else the population record's `ensembl_vep` engine version, else the
    release number in the release record's reference prefix. The wall-time table's version stays the
    engine's own `engine_version`, which a carried-forward comparator block keeps."""
    for block in ((release or {}).get("perl_reference"), pop.get("perl_reference"), population.get("perl_reference")):
        if isinstance(block, dict):
            try:
                return release_from_provenance(block)
            except ValueError as exc:
                fail(f"a perl_reference block names no release: {exc}")
    engine = (population.get("wall_time") or {}).get("ensembl_vep") or {}
    m = re.search(r"(\d+(?:\.\d+)?)$", str(engine.get("engine_version", "")))
    if m:
        return m.group(1)
    m = re.search(r"_r(\d+(?:\.\d+)?)_", str(pop.get("perl_reference_prefix", "")))
    if m:
        return m.group(1)
    fail("neither record names the Ensembl VEP version")


def named_releases(where: str, block) -> dict[str, str]:
    """The Ensembl VEP release each part of a `perl_reference` block names, keyed by its place: one part per
    suite family when the block is keyed by family, else the block itself; an absent block names nothing."""
    if not isinstance(block, dict):
        return {}
    parts = {f"{where}.{f}": block[f] for f in FAMILIES if isinstance(block.get(f), dict)} or {where: block}
    out: dict[str, str] = {}
    for name, part in parts.items():
        try:
            out[name] = release_from_provenance(part)
        except ValueError as exc:
            fail(f"the perl_reference block of {name} names no release: {exc}")
    return out


def reference_release(release: dict, pop: dict, population: dict) -> str:
    """The one Ensembl VEP release the records' `perl_reference` blocks name: the release record's, the one
    under its `population_concordance` and the population record's. A record without a block is refused
    rather than read from an engine version, which a carried-forward wall-time run keeps at an earlier
    release; blocks naming different releases are refused."""
    named = named_releases("the release record", release.get("perl_reference"))
    if not named:
        fail("the release record carries no perl_reference block")
    named.update(named_releases("population_concordance", pop.get("perl_reference")))
    from_population = named_releases("the population record", population.get("perl_reference"))
    if not from_population:
        fail("the population record carries no perl_reference block")
    named.update(from_population)
    if len(set(named.values())) != 1:
        fail("the records name different Ensembl VEP releases: " + ", ".join(f"{k} {v}" for k, v in named.items()))
    return next(iter(named.values()))


def wall_time_sums(population: dict) -> dict[tuple[str, str], Decimal]:
    """Per (dataset, arch), the sum of the per-chromosome wall-time medians of the released binary; a cell's
    dataset, arch or median is refused by name when absent, and its median when it is not a finite number."""
    sums: dict[tuple[str, str], Decimal] = {}
    for index, cell in enumerate(population["wall_time"]["vep_rs"]["per_cell"]):
        where = f"the population record's wall_time.vep_rs.per_cell[{index}]"
        dataset, arch, median = required(cell, ("dataset", "arch", "median_sec"), where)
        key = (dataset, arch)
        sums[key] = sums.get(key, Decimal(0)) + median_seconds(median, where)
    return sums


def peak_rss_gb(population: dict) -> str:
    """The largest per-chromosome median peak RSS of the released binary, the record's MB over 1000 at one
    decimal, the README's own reading of the same field."""
    peak = max(float(c["peak_rss_mb_median"]) for c in population["wall_time"]["vep_rs"]["per_cell"])
    return half_up(peak / 1000, 1)


def short_label(label: str) -> str:
    """A dataset label without its ', whole genome' suffix: the caption carries the scale."""
    return label.replace(", whole genome", "")


def population_counts(entry, where: str) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    """A population dataset's (intersection, Ensembl VEP tuples, vep-rs tuples), raw and then less the excluded
    counts, each count refused by name when it is absent or not an integer."""
    keys = ("intersection", "perl_tuples", "vep_rs_tuples", "excluded_perl", "excluded_rust")
    values = required(entry, keys, where)
    i, a, b, excluded_perl, excluded_rust = (integer(v, k, where) for v, k in zip(values, keys))
    return (i, a, b), (i, a - excluded_perl, b - excluded_rust)


def whole_genome_row(label: str, entry: dict, arm: Decimal, x86: Decimal, where: str) -> str:
    """A whole-genome row: raw F1 over every tuple and adjusted F1 over the tuples left after the excluded
    counts, both from the integer counts, then the wall time per architecture."""
    raw, adjusted = population_counts(entry, where)
    return f"| {label} | {f1_9(*raw)} | {f1_9(*adjusted)} | {whole_seconds(arm)} | {whole_seconds(x86)} |"


def whole_genome_setup(population: dict) -> tuple[dict, dict[str, str], dict[tuple[str, str], Decimal], list[str]]:
    """Where both pages' whole-genome tables start: the released binary's wall-time block, its instance type
    per architecture, the per-(dataset, arch) sums of its per-chromosome medians and the table's header rows."""
    vep_rs = population["wall_time"]["vep_rs"]
    header = ["| Dataset | Raw F1 | Adjusted F1 | ARM (s) | x86 (s) |", "| --- | --- | --- | --- | --- |"]
    return vep_rs, vep_rs["instance_types"], wall_time_sums(population), header


def dataset_wall_time(dataset: str, sums: dict[tuple[str, str], Decimal]) -> tuple[Decimal, Decimal]:
    """A population dataset's summed wall time on ARM, then on x86; one the population record did not time on
    both architectures is refused."""
    if any((dataset, arch) not in sums for arch in ARCH_LABELS):
        fail(f"the population record has no wall time for {dataset}")
    return sums[(dataset, "arm64")], sums[(dataset, "x86_64")]


def required(block, keys: tuple[str, ...], where: str) -> tuple:
    """The values under `keys` of a record block, refused by name when any is absent or the block is not an
    object."""
    values = block if isinstance(block, dict) else {}
    missing = [k for k in keys if values.get(k) is None]
    if missing:
        fail(f"{where} carries no {', '.join(missing)}")
    return tuple(values[k] for k in keys)


def integer(value, key: str, where: str) -> int:
    """A record's integer count, refused by name when it is not one: `int()` would read 19.6 as 19 and `true`
    as 1, and `bool` is a subclass of `int`."""
    if isinstance(value, bool) or not isinstance(value, int):
        fail(f"{where} carries {key} that is not an integer: {value!r}")
    return value


def median_seconds(value, where: str) -> Decimal:
    """A record's `median_sec` as a Decimal, refused by name when it is not a number or not finite."""
    try:
        seconds = Decimal(str(value))
    except InvalidOperation:
        fail(f"{where} carries a median_sec that is not a number: {value!r}")
    if not seconds.is_finite():
        fail(f"{where} carries a median_sec that is not finite: {value!r}")
    return seconds


def int_triple(block, keys: tuple[str, str, str], where: str) -> tuple[int, int, int]:
    """Three integer counts of a record block, each refused by name when absent or not an integer, so a record's
    rounded F1 field is never printed in their place."""
    i, a, b = required(block, keys, where)
    return integer(i, keys[0], where), integer(a, keys[1], where), integer(b, keys[2], where)


def whole_file_aggregate(release: dict, suite: str, arch: str) -> dict:
    """The release record's `per_cell_aggregates` entry of one ClinVar cell and architecture, the released
    binary's whole-file wall time; refused when absent or repeated, as is a `per_cell_aggregates` that is
    not a list of objects."""
    entries = release.get("per_cell_aggregates")
    if not isinstance(entries, list) or not all(isinstance(a, dict) for a in entries):
        fail("the release record's per_cell_aggregates is not a list of objects")
    hits = [
        a
        for a in entries
        if a.get("suite_id") == suite and a.get("arch") == arch and a.get("engine", "vep_rs") == "vep_rs"
    ]
    if len(hits) != 1:
        fail(f"the release record has {len(hits)} per_cell_aggregates entries for {suite} on {arch}; one is read")
    return hits[0]


def whole_file_wall_time(release: dict, suite: str, arch: str) -> tuple[Decimal, str, int]:
    """One ClinVar cell's aggregate as (median seconds, instance type, clone count), each refused by name when
    absent, the median when it is not a finite number and the clone count when it is not an integer."""
    where = f"the per_cell_aggregates entry for {suite} on {arch}"
    keys = ("median_sec", "instance_type", "n_clones")
    median, instance_type, n_clones = required(whole_file_aggregate(release, suite, arch), keys, where)
    return median_seconds(median, where), str(instance_type), integer(n_clones, "n_clones", where)


def one(values: set, what: str) -> str:
    """The one value the ClinVar cells agree on, refused when they name several."""
    if len(values) != 1:
        fail(f"the ClinVar cells name {len(values)} {what}: {sorted(values)}")
    return str(next(iter(values)))


def clone_count_clause(clones: dict[tuple[str, str], int]) -> str:
    """The clone count the ClinVar cells agree on, else each cell's own in the table's order, ClinVar GRCh37 then
    GRCh38 and ARM before x86 within each, whatever the record's order:
    `19 (ClinVar GRCh37, ARM), 20 (ClinVar GRCh37, x86), ... and 20 (ClinVar GRCh38, x86)`."""
    if len(set(clones.values())) == 1:
        return str(next(iter(clones.values())))
    parts = [f"{n} ({CLINVAR[suite]}, {ARCH_LABELS[arch]})" for (suite, arch), n in clones.items()]
    return f"{', '.join(parts[:-1])} and {parts[-1]}"


def file_counts(block, where: str) -> dict[str, tuple[int, int, int]]:
    """A `per_file` or `per_file_adjusted` block as {file: (intersection, Ensembl VEP tuples, vep-rs tuples)}."""
    if not isinstance(block, dict) or not block:
        fail(f"{where} is absent or empty")
    keys = ("intersection", "perl_tuples", "engine_tuples")
    return {name: int_triple(entry, keys, f"{where} {name}") for name, entry in block.items()}


def sv_sets(cell: dict, where: str) -> list[tuple[str, list, list]]:
    """A structural-variant cell's file sets as (label, raw cells, adjusted cells): the ClinVar, gnomAD and
    1000 Genomes sets, each one file found by the prefix of its name, then the synthetic feature files,
    named by a two-digit index, pooled. The two per-file blocks must list the same files, each real set
    one file and the synthetic set SYNTHETIC_FILES; any other file is refused."""
    listing = f"{where}.per_file"
    raw = file_counts(cell.get("per_file"), listing)
    adjusted = file_counts(cell.get("per_file_adjusted"), f"{where}.per_file_adjusted")
    if raw.keys() != adjusted.keys():
        fail(f"{where} per_file and per_file_adjusted list different files")
    real: dict[str, list[str]] = {}
    synthetic: list[str] = []
    for name in raw:
        label = next((label for prefix, label in SV_SETS if name.startswith(prefix)), None)
        if label is not None:
            real.setdefault(label, []).append(name)
        elif SYNTHETIC.match(name):
            synthetic.append(name)
        else:
            fail(f"{listing} lists {name}, a file of no known set")
    sets: list[tuple[str, list, list]] = []
    for _, label in SV_SETS:
        files = real.get(label, [])
        if len(files) != 1:
            fail(f"{listing} lists {len(files)} files of the {label} set; one is read")
        sets.append((label, [raw[files[0]]], [adjusted[files[0]]]))
    if len(synthetic) != SYNTHETIC_FILES:
        fail(f"{listing} lists {len(synthetic)} synthetic files; the set has {SYNTHETIC_FILES}")
    sets.append((f"Synthetic ({SYNTHETIC_FILES} files)", [raw[n] for n in synthetic], [adjusted[n] for n in synthetic]))
    return sets


def six_dataset_blocks(
    release: dict, conc: dict, pop: dict, population: dict, record_dir: str, version: str
) -> list[str]:
    """The page of the six-dataset shape: the whole-genome table (the four population datasets in the
    record's order, the two genome-wide ClinVar files, then `All six` pooled from the summed integer counts
    of all six) under a caption naming the reference release, both protocols, the instance types and the
    peak memory, then the structural-variant F1 table (per assembly the ClinVar, gnomAD and 1000 Genomes
    sets and the synthetic files pooled) under the caption that the paper's figures stay as published.
    `All six` carries no wall time: the ClinVar whole-file medians were measured on another instance class
    and protocol than the per-chromosome sums, so their sum is not a figure. A record's chromosome 21 cells
    render no row."""
    vep_rs, types, sums, rows = whole_genome_setup(population)
    raw: list[tuple[int, int, int]] = []
    adjusted: list[tuple[int, int, int]] = []
    datasets = pop["per_dataset"]
    if len(datasets) != 4:
        fail(
            f"the population concordance has {len(datasets)} datasets; "
            "the six-dataset table pools four with the two ClinVar files"
        )
    for index, r in enumerate(datasets):
        dataset, label = required(r, ("dataset", "label"), f"population_concordance.per_dataset[{index}]")
        where = f"population_concordance.per_dataset {dataset}"
        arm, x86 = dataset_wall_time(dataset, sums)
        raw_counts, adjusted_counts = population_counts(r, where)
        raw.append(raw_counts)
        adjusted.append(adjusted_counts)
        rows.append(whole_genome_row(short_label(str(label)), r, arm, x86, where))
    clinvar_types: dict[str, set[str]] = {"arm64": set(), "x86_64": set()}
    clones: dict[tuple[str, str], int] = {}
    for suite, label in CLINVAR.items():
        cell, where = conc[suite], f"concordance.{suite}"
        i, a, b = int_triple(cell, ("intersection", "perl", "rust"), where)
        adj = int_triple(cell, ("adj_intersection", "adj_perl", "adj_rust"), where)
        raw.append((i, a, b))
        adjusted.append(adj)
        medians: dict[str, str] = {}
        for arch, instance_types in clinvar_types.items():
            seconds, instance_type, n_clones = whole_file_wall_time(release, suite, arch)
            medians[arch] = whole_seconds(seconds)
            instance_types.add(instance_type)
            clones[(suite, arch)] = n_clones
        rows.append(f"| {label} | {f1_9(i, a, b)} | {f1_9(*adj)} | {medians['arm64']} | {medians['x86_64']} |")
    rows.append(f"| All six | {pooled_f1(raw)} | {pooled_f1(adjusted)} | | |")
    caption = (
        f"Whole genome against Ensembl VEP {reference_release(release, pop, population)}: the population datasets "
        "one run per chromosome, each ClinVar release one genome-wide file; wall time is the sum of the "
        f"per-chromosome medians of {int(vep_rs['n_clones_per_cell'])} machines at {int(population['fork'])} threads on "
        f"`{types['arm64']}` (ARM) and `{types['x86_64']}` (x86) for the population datasets, and the whole-file "
        f"median of {clone_count_clause(clones)} machines on `{one(clinvar_types['arm64'], 'ARM instance types')}` (ARM) "
        f"and `{one(clinvar_types['x86_64'], 'x86 instance types')}` (x86) for the ClinVar files, another instance class "
        "and protocol, so the `All six` row carries no wall time; peak memory at most "
        f"{peak_rss_gb(population)} GB on any chromosome; tuple counts, per-chromosome rows and per-file counts are in "
        f"the [provenance records]({REPO_URL}/tree/v{version}/{record_dir})."
    )
    sv_caption = (
        "Structural-variant file sets measured on the released binary against the same reference, each set's files "
        "pooled by summed counts; the paper's own figures stay as published in "
        f"[docs/published-figures.md]({REPO_URL}/blob/v{version}/docs/published-figures.md)."
    )
    sv = ["| File set | Assembly | Raw F1 | Adjusted F1 |", "| --- | --- | --- | --- |"]
    for suite, assembly in SV.items():
        for label, raw_files, adjusted_files in sv_sets(conc[suite], f"concordance.{suite}"):
            sv.append(f"| {label} | {assembly} | {pooled_f1(raw_files)} | {pooled_f1(adjusted_files)} |")
    return ["## Concordance and wall time", caption, "\n".join(rows), sv_caption, "\n".join(sv)]


def measurements_blocks(repo: Path, version: str) -> list[str]:
    """The `## Concordance and wall time` blocks for the record's shape: six_dataset_blocks for a record of
    the six-dataset shape, else the earlier page: the whole-genome table (one row per population dataset
    and the pooled row: raw and adjusted F1 from the release record's counts, wall time per architecture
    from the population record) under its caption, then the chromosome 21 table over the eight scored
    suites under its caption. A suite row's adjusted F1 is the record's `adj_f1`: the suite rows carry no
    excluded counts to derive it from."""
    release, conc, pop, record_dir = load_release_record(repo, version)
    population = load_population_record(repo, version)
    if six_dataset_shape(release, conc):
        return six_dataset_blocks(release, conc, pop, population, record_dir, version)
    vep_rs, types, sums, rows = whole_genome_setup(population)
    caption = (
        f"Whole genome against Ensembl VEP {comparator_version(population, pop, release)}, one run per chromosome; "
        f"wall time is the sum of the per-chromosome medians of {int(vep_rs['n_clones_per_cell'])} machines "
        f"at {int(population['fork'])} threads on `{types['arm64']}` (ARM) and `{types['x86_64']}` (x86); "
        f"peak memory at most {peak_rss_gb(population)} GB on any chromosome; tuple counts and "
        f"per-chromosome rows are in the [provenance records]({REPO_URL}/tree/v{version}/{record_dir})."
    )
    total_arm = total_x86 = Decimal(0)
    for index, r in enumerate(pop["per_dataset"]):
        dataset, label = required(r, ("dataset", "label"), f"population_concordance.per_dataset[{index}]")
        arm, x86 = dataset_wall_time(dataset, sums)
        total_arm += arm
        total_x86 += x86
        where = f"population_concordance.per_dataset {dataset}"
        rows.append(whole_genome_row(short_label(str(label)), r, arm, x86, where))
    p = pop["pooled"]
    (label,) = required(p, ("label",), "population_concordance.pooled")
    rows.append(whole_genome_row(str(label), p, total_arm, total_x86, "population_concordance.pooled"))
    chr21_caption = (
        "The paper's chromosome 21 cells measured on the released binary; the paper's own figures stay as "
        f"published in [docs/published-figures.md]({REPO_URL}/blob/v{version}/docs/published-figures.md)."
    )
    chr21 = ["| Dataset | Assembly | Raw F1 | Adjusted F1 |", "| --- | --- | --- | --- |"]
    for s in SCORED:
        r = conc[s]
        dataset, assembly = SUITE_LABELS[s]
        chr21.append(
            f"| {dataset} | {assembly} | {f1_9(r['intersection'], r['perl'], r['rust'])} | {half_up(float(r['adj_f1']), 9)} |"
        )
    return ["## Concordance and wall time", caption, "\n".join(rows), chr21_caption, "\n".join(chr21)]


def read_sums(path: Path) -> dict[str, str]:
    """SHA256SUMS as {file: digest}."""
    if not path.exists():
        fail(f"no SHA256SUMS at {path}")
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[1].lstrip("*")] = parts[0]
    return out


def install_bullets(
    version: str, sums: dict[str, str], commit: str, toolchain: str, rust_version: str, version_doi: str | None
) -> str:
    """One bullet per binary archive in the order of TARGETS, then the archive contents and the source build,
    the container image, the build and its attestation, and the source archive. Every target's archive and
    the source archive must be in SHA256SUMS, and every binary archive there must be a known target."""
    tag = f"v{version}"
    dl = f"{REPO_URL}/releases/download/{tag}"
    prefix, suffix = f"vep-{version}-", ".tar.gz"
    missing = [t for t in TARGETS if f"{prefix}{t}{suffix}" not in sums]
    if missing:
        fail(f"SHA256SUMS lacks the archive for {missing}")
    unknown = [
        n for n in sums if n.startswith(prefix) and n.endswith(suffix) and n[len(prefix) : -len(suffix)] not in TARGETS
    ]
    if unknown:
        fail(f"SHA256SUMS lists {unknown}, archives of no known target")
    src = f"vep-rs-{version}.tar.gz"
    if src not in sums:
        fail(f"SHA256SUMS lacks {src}")
    bullets = [
        f"- [`{prefix}{triple}{suffix}`]({dl}/{prefix}{triple}{suffix}): {platform}, {floor}; "
        f"sha256 `{sums[f'{prefix}{triple}{suffix}']}`"
        for triple, (platform, floor) in TARGETS.items()
    ]
    doi = f"; version DOI [{version_doi}](https://doi.org/{version_doi})" if version_doi else ""
    bullets += [
        f"- Each archive unpacks to `vep-{version}-<target>/` holding `vep`, `vep-cache-builder` and "
        "`vep-cache-converter`; a binary below its CPU floor aborts with an illegal-instruction fault, so on "
        f"such a host build from source (Rust {rust_version} or later), which compiles for the machine it runs "
        'on: `RUSTFLAGS="-C target-cpu=native" cargo build --release --locked`',
        f"- Container image: `docker pull ghcr.io/natera-open-source/vep-rs:{version}` (linux/amd64: the x86_64 "
        "binaries with the `duckdb` CLI that Parquet output needs)",
        f"- Built by the [release workflow]({REPO_URL}/actions/workflows/release.yml) from tag `{tag}` "
        f"(commit `{commit}`) with Rust {toolchain}; every file in [`SHA256SUMS`]({dl}/SHA256SUMS) carries a "
        "GitHub build-provenance attestation: `gh attestation verify <file> --repo natera-open-source/vep-rs`",
        f"- [`{src}`]({dl}/{src}) is `git archive --format=tar.gz --prefix=vep-rs-{version}/ {tag}`; "
        f"sha256 `{sums[src]}`{doi}",
    ]
    return "\n".join(bullets)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--version", required=True)
    ap.add_argument("--previous", required=True, help="previous release tag, e.g. v0.1.0")
    ap.add_argument("--commit", required=True)
    ap.add_argument("--toolchain", required=True)
    ap.add_argument("--sha256sums", required=True, type=Path)
    ap.add_argument("--notes", type=Path, help=".github/release-notes/vX.Y.Z.md")
    ap.add_argument("--changelog-section", help="defaults to the version")
    ap.add_argument("--date", help="release date; defaults to the CHANGELOG heading's")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()

    v, tag = a.version, f"v{a.version}"
    repo = a.repo
    section, cl_date = changelog_section((repo / "CHANGELOG.md").read_text(), a.changelog_section or v)
    date = a.date or cl_date
    if not date:
        fail("no date: pass --date or date the CHANGELOG heading")
    rust_version = re.search(r'^rust-version\s*=\s*"([^"]+)"', (repo / "Cargo.toml").read_text(), re.MULTILINE)
    if not rust_version:
        fail("no rust-version in Cargo.toml")
    notes = a.notes or repo / ".github" / "release-notes" / f"{tag}.md"
    hand = hand_written(notes)
    check_hand_written(hand, section, notes)
    sums = read_sums(a.sha256sums)
    version_doi = read_doi(repo / "CITATION.cff", v)
    cite_doi = version_doi or CONCEPT_DOI

    blocks = [hand["opener"]]
    if hand["highlights"]:
        blocks += ["## Highlights", hand["highlights"]]
    if hand["upgrade"]:
        blocks += ["## Upgrade notes", hand["upgrade"]]
    blocks += kind_blocks(section)
    blocks += measurements_blocks(repo, v)
    blocks += [
        "## Install",
        install_bullets(v, sums, a.commit, a.toolchain, rust_version.group(1), version_doi),
        "## Cite",
        f"Porter M, Borkowski R. vep-rs, version {v}. {date[:4]}. doi:[{cite_doi}](https://doi.org/{cite_doi})",
        f"**Full changelog**: [CHANGELOG]({REPO_URL}/blob/{tag}/CHANGELOG.md) · "
        f"[`{a.previous}...{tag}`]({REPO_URL}/compare/{a.previous}...{tag})",
    ]
    body = "\n\n".join(blocks) + "\n"
    lines = body.count("\n")
    if lines > LINE_CAP:
        fail(f"the page is {lines} lines; the cap is {LINE_CAP}")
    if a.out:
        a.out.write_text(body)
    else:
        sys.stdout.write(body)


if __name__ == "__main__":
    main()
