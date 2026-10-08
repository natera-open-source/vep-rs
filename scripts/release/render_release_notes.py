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
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

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
LIST_ITEM = re.compile(r"^\s*(?:[-*+]|\d+\.)\s")
BULLET = re.compile(r"[-*] ")
SENTENCE_END = re.compile(r"[.!?][\"'`)\]]*(?=\s|$)")


def fail(message: str) -> None:
    sys.exit(f"ERROR: [render_release_notes] {message}")


def half_up(x: float, places: int) -> str:
    return str(Decimal(str(x)).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))


def f1_9(intersection, a, b) -> str:
    """F1 as 2I/(A+B) from the integer tuple counts at nine decimal places; an exact match reads 1.000000000.
    Computed from the counts rather than copied from a record's rounded field, so the page never prints
    a coarser rounding than the counts support."""
    q = Decimal(2 * int(intersection)) / Decimal(int(a) + int(b))
    return str(q.quantize(Decimal("0.000000001"), rounding=ROUND_HALF_UP))


def whole_seconds(x: Decimal) -> str:
    return f"{int(x.quantize(Decimal(1), rounding=ROUND_HALF_UP)):,}"


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


def load_release_record(repo: Path, version: str) -> tuple[dict, dict, str]:
    """The release record's per-suite concordance, its whole-genome concordance and the record directory's path."""
    record_path = latest_record(repo, "release", version)
    record = json.loads(record_path.read_text())
    conc = record.get("concordance") or {}
    missing = [s for s in SCORED if s not in conc]
    if missing:
        fail(f"{record_path.name} lacks concordance rows for {missing}")
    pop = record.get("population_concordance")
    if not pop:
        fail(f"{record_path.name} lacks population_concordance")
    return conc, pop, record_path.relative_to(repo).parent.as_posix()


def load_population_record(repo: Path, version: str) -> dict:
    """The release's population record: the whole-genome wall time per chromosome and architecture."""
    return json.loads(latest_record(repo, "population", version).read_text())


def comparator_version(population: dict, pop: dict) -> str:
    """The Ensembl VEP release the whole genome was scored against: the population record's
    `ensembl_vep` engine version, else the release number in the release record's reference prefix."""
    engine = (population.get("wall_time") or {}).get("ensembl_vep") or {}
    m = re.search(r"(\d+(?:\.\d+)?)$", str(engine.get("engine_version", "")))
    if m:
        return m.group(1)
    m = re.search(r"_r(\d+(?:\.\d+)?)_", str(pop.get("perl_reference_prefix", "")))
    if m:
        return m.group(1)
    fail("neither record names the Ensembl VEP version")


def wall_time_sums(population: dict) -> dict[tuple[str, str], Decimal]:
    """Per (dataset, arch), the sum of the per-chromosome wall-time medians of the released binary."""
    sums: dict[tuple[str, str], Decimal] = {}
    for cell in population["wall_time"]["vep_rs"]["per_cell"]:
        key = (cell["dataset"], cell["arch"])
        sums[key] = sums.get(key, Decimal(0)) + Decimal(str(cell["median_sec"]))
    return sums


def peak_rss_gb(population: dict) -> str:
    """The largest per-chromosome median peak RSS of the released binary, the record's MB over 1000 at one
    decimal, the README's own reading of the same field."""
    peak = max(float(c["peak_rss_mb_median"]) for c in population["wall_time"]["vep_rs"]["per_cell"])
    return half_up(peak / 1000, 1)


def short_label(label: str) -> str:
    """A dataset label without its ', whole genome' suffix: the caption carries the scale."""
    return label.replace(", whole genome", "")


def whole_genome_row(label: str, entry: dict, arm: Decimal, x86: Decimal) -> str:
    """A whole-genome row: raw F1 over every tuple and adjusted F1 over the tuples left after the excluded
    counts, both from the integer counts, then the wall time per architecture."""
    i, a, b = int(entry["intersection"]), int(entry["perl_tuples"]), int(entry["vep_rs_tuples"])
    adjusted = f1_9(i, a - int(entry["excluded_perl"]), b - int(entry["excluded_rust"]))
    return f"| {label} | {f1_9(i, a, b)} | {adjusted} | {whole_seconds(arm)} | {whole_seconds(x86)} |"


def measurements_blocks(repo: Path, version: str) -> list[str]:
    """The whole-genome table (one row per population dataset and the pooled row: raw and adjusted F1
    from the release record's counts, wall time per architecture from the population record) under its
    caption, then the chromosome 21 table over the eight scored suites under its caption. A suite row's
    adjusted F1 is the record's `adj_f1`: the suite rows carry no excluded counts to derive it from."""
    conc, pop, record_dir = load_release_record(repo, version)
    population = load_population_record(repo, version)
    vep_rs = population["wall_time"]["vep_rs"]
    types = vep_rs["instance_types"]
    sums = wall_time_sums(population)
    caption = (
        f"Whole genome against Ensembl VEP {comparator_version(population, pop)}, one run per chromosome; "
        f"wall time is the sum of the per-chromosome medians of {int(vep_rs['n_clones_per_cell'])} machines "
        f"at {int(population['fork'])} threads on `{types['arm64']}` (ARM) and `{types['x86_64']}` (x86); "
        f"peak memory at most {peak_rss_gb(population)} GB on any chromosome; tuple counts and "
        f"per-chromosome rows are in the [provenance records]({REPO_URL}/tree/v{version}/{record_dir})."
    )
    rows = ["| Dataset | Raw F1 | Adjusted F1 | ARM (s) | x86 (s) |", "| --- | --- | --- | --- | --- |"]
    total = {"arm64": Decimal(0), "x86_64": Decimal(0)}
    for r in pop["per_dataset"]:
        dataset = r["dataset"]
        if any((dataset, arch) not in sums for arch in total):
            fail(f"the population record has no wall time for {dataset}")
        for arch in total:
            total[arch] += sums[(dataset, arch)]
        rows.append(whole_genome_row(short_label(r["label"]), r, sums[(dataset, "arm64")], sums[(dataset, "x86_64")]))
    p = pop["pooled"]
    rows.append(whole_genome_row(p["label"], p, total["arm64"], total["x86_64"]))
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
