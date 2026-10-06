#!/usr/bin/env python3
"""Render the GitHub release body for a vep-rs version from the repository's own files.

Every line of the body except the hand-written summary, upgrade notes and known issues is derived:
the CHANGELOG section for the version, the release's provenance record, SHA256SUMS, CITATION.cff's DOIs,
the workspace rust-version and the build's commit and toolchain.

usage: render_release_notes.py --repo DIR --version X.Y.Z --previous vX.Y.Z --commit SHA
                               --toolchain 1.97.1 [--sha256sums FILE] [--notes FILE]
                               [--changelog-section NAME] [--date YYYY-MM-DD] [--out FILE]

The release workflow runs it on a tag (the CHANGELOG section of the version, dated by the release
commit) and on a dry run (`--changelog-section Unreleased --date <today>`).
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
SCORED = ("s01", "s02", "s03", "s04", "s05", "s06", "s07", "s08")
SNP_INDEL = ("s01", "s02", "s03", "s04", "s05", "s06")
TARGETS = (
    (
        "x86_64-unknown-linux-gnu",
        "Linux x86_64 (glibc 2.39 or later)",
        "x86-64-v3: AVX2, BMI2, FMA (Intel Haswell 2013 or later, AMD Zen or later)",
    ),
    (
        "aarch64-unknown-linux-gnu",
        "Linux aarch64 (glibc 2.39 or later)",
        "Armv8.4-A with SVE (Arm Neoverse V1 or later: AWS Graviton3 and later; not Graviton2, Ampere Altra or Raspberry Pi)",
    ),
    (
        "aarch64-apple-darwin",
        "macOS 11 or later on Apple silicon",
        "any Apple M-series",
    ),
)
ABOUT = (
    "> vep-rs is a from-scratch Rust implementation of the Ensembl Variant Effect Predictor: the same\n"
    "> inputs and flags, VEP's consequence vocabulary in VEP's output formats plus Parquet, measured for\n"
    "> concordance against Ensembl VEP at population scale."
)


def half_up(x: float, places: int) -> str:
    return str(
        Decimal(str(x)).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP)
    )


def f1_9(intersection, a, b) -> str:
    """F1 as 2I/(A+B) from the integer tuple counts at nine decimal places; an exact match reads 1.000000000.
    Computed from the counts rather than copied from a record's rounded field, so the page never prints
    a coarser rounding than the counts support."""
    q = Decimal(2 * int(intersection)) / Decimal(int(a) + int(b))
    return str(q.quantize(Decimal("0.000000001"), rounding=ROUND_HALF_UP))


def changelog_section(text: str, name: str) -> tuple[str, str | None]:
    """Return (body, date) of the `## [name]` section, body without the heading."""
    pattern = re.compile(
        rf"^## \[{re.escape(name)}\](?: - (\d{{4}}-\d{{2}}-\d{{2}}))?[^\n]*\n", re.MULTILINE
    )
    m = pattern.search(text)
    if not m:
        sys.exit(
            f"ERROR: [render_release_notes] no '## [{name}]' section in CHANGELOG.md"
        )
    rest = text[m.end() :]
    nxt = re.search(r"^## \[", rest, re.MULTILINE)
    body = rest[: nxt.start()] if nxt else rest
    body = re.sub(r"\n{3,}", "\n\n", body).strip("\n")
    return body, m.group(1)


def hand_written(path: Path | None) -> dict[str, str]:
    """The maintainer's lines: the first paragraph is the summary; optional '## Upgrade notes' and
    '## Known issues' sections are copied through."""
    out = {"summary": "", "upgrade": "", "known": ""}
    if path is None or not path.exists():
        return out
    text = path.read_text()
    text = re.sub(r"^# .*\n", "", text, count=1)
    parts = re.split(r"^## ", text, flags=re.MULTILINE)
    out["summary"] = parts[0].strip()
    for part in parts[1:]:
        title, _, body = part.partition("\n")
        key = {"Upgrade notes": "upgrade", "Known issues": "known"}.get(title.strip())
        if key is None:
            sys.exit(
                f"ERROR: [render_release_notes] unknown section '## {title.strip()}' in {path}"
            )
        out[key] = body.strip()
    return out


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


def load_release_measurements(repo: Path, version: str):
    """The release's provenance record (docs/concordance-provenance/<date>-release-v<version>.json):
    its per-suite concordance and per-cell wall-time aggregates. The record, not a file under
    manuscript/data, carries a release's own measurements."""
    records = sorted((repo / "docs" / "concordance-provenance").glob(f"*-release-v{version}.json"))
    if not records:
        sys.exit(
            f"ERROR: [render_release_notes] no docs/concordance-provenance/*-release-v{version}.json"
        )
    record_path = records[-1]
    record = json.loads(record_path.read_text())
    conc = record.get("concordance") or {}
    walls = {
        (cell["arch"], cell["suite_id"]): (float(cell["median_sec"]), int(cell["n_clones"]))
        for cell in record.get("per_cell_aggregates") or []
    }
    missing = [s for s in SCORED if s not in conc]
    if missing:
        sys.exit(
            f"ERROR: [render_release_notes] {record_path.name} lacks concordance rows for {missing}"
        )
    return conc, walls, record_path.relative_to(repo).as_posix(), record.get("population_concordance")


def load_population_record(repo: Path, version: str) -> dict | None:
    """The release's population record (docs/concordance-provenance/<date>-population-v<version>.json):
    the whole-genome wall time per dataset and architecture. None when the release has none."""
    records = sorted((repo / "docs" / "concordance-provenance").glob(f"*-population-v{version}.json"))
    return json.loads(records[-1].read_text()) if records else None


def population_wall_time_lines(population: dict) -> list[str]:
    """The whole-genome wall-time line: per dataset, the sum of the per-chromosome medians on each
    architecture, with the instance types the record names."""
    wt = population["wall_time"]["vep_rs"]
    types = wt["instance_types"]
    n = int(wt["n_clones_per_cell"])
    labels = {d["dataset"]: d["label"].replace(", whole genome", "") for d in population["datasets"]}
    sums = {(c["dataset"], c["arch"]): float(c["sum_of_cell_medians_sec"]) for c in wt["per_dataset"]}
    cells = [
        f"{labels[d['dataset']]} {int(round(sums[(d['dataset'], 'arm64')])):,} s ARM / "
        f"{int(round(sums[(d['dataset'], 'x86_64')])):,} s x86"
        for d in population["datasets"]
    ]
    return [
        "",
        f"Whole-genome wall time as the sum of the per-chromosome medians, each the median of {n} independent",
        f"machines, ARM Graviton4 (`{types['arm64']}`) and x86 Intel (`{types['x86_64']}`), 16 threads, local NVMe,",
        "sites-only inputs, the discarded warmup's output deleted before the timed run: " + "; ".join(cells) + ".",
    ]


def population_lines(pop: dict) -> list[str]:
    """The whole-genome table: one row per population dataset in the record's order and the
    pooled row, from the record's population_concordance section, whose rows carry their own
    labels (variants are VCF records; tuples are output rows reduced to Location, Allele,
    Feature, Feature_type and Consequence set)."""
    lines = [
        "",
        "Whole-genome concordance of the released binary on four population datasets, one run per",
        "chromosome on ARM Graviton4, against Ensembl VEP 115.2's output on the same inputs; a variant is",
        "one VCF record, a tuple one output row reduced to its Location, Allele, Feature, Feature_type and",
        "Consequence set, and F1 pools the chromosomes by summed counts:",
        "",
        "| Dataset | Assembly | Chromosomes | Variants | VEP tuples | vep-rs tuples | Matched | Raw F1 | Adjusted F1 |",
        "| ------- | -------- | ----------- | -------- | ---------- | ------------- | ------- | ------ | ----------- |",
    ]
    for r in pop["per_dataset"]:
        lines.append(
            f"| {r['label']} | {r['assembly']} | {int(r['shards'])} | {int(r['variants']):,} | {int(r['perl_tuples']):,} | "
            f"{int(r['vep_rs_tuples']):,} | {int(r['intersection']):,} | {f1_9(r['intersection'], r['perl_tuples'], r['vep_rs_tuples'])} | {half_up(float(r['adjusted_f1']), 9)} |"
        )
    p = pop["pooled"]
    lines.append(
        f"| {p['label']} | GRCh37 and GRCh38 | {int(p['shards'])} | {int(p['variants']):,} | {int(p['perl_tuples']):,} | "
        f"{int(p['vep_rs_tuples']):,} | {int(p['intersection']):,} | {f1_9(p['intersection'], p['perl_tuples'], p['vep_rs_tuples'])} | {half_up(float(p['adjusted_f1']), 9)} |"
    )
    return lines


def measurements_section(repo: Path, version: str, date: str) -> str:
    """The whole-genome concordance table and wall-time line first (the release record's
    population_concordance and the population record), then the chromosome 21 suites: the eight-row
    table and their wall-time line."""
    conc, walls, record_rel, pop = load_release_measurements(repo, version)
    population = load_population_record(repo, version)
    record_dir = record_rel.rsplit("/", 1)[0]
    lines = [
        "## Concordance and wall time for this release",
        "",
        "The paper's published figures are unchanged; these are the released version's own, from the",
        f"release's provenance records under [`{record_dir}/`]({REPO_URL}/tree/v{version}/{record_dir})",
        "(method: `scripts/concordance/run_clone_measurement.sh`).",
    ]
    if pop:
        lines += population_lines(pop)
    if population:
        lines += population_wall_time_lines(population)
    if pop or population:
        lines += ["", "The chromosome 21 suites, the paper's cells measured on the released binary:"]
    lines += [
        "",
        "| Dataset | Assembly | Raw F1 | Adjusted F1 | VEP tuples | vep-rs tuples | Matched |",
        "| ------- | -------- | ------ | ----------- | ---------- | ------------- | ------- |",
    ]
    for s in SCORED:
        r = conc[s]
        dataset, assembly = SUITE_LABELS[s]
        lines.append(
            f"| {dataset} | {assembly} | {f1_9(r['intersection'], r['perl'], r['rust'])} | {half_up(float(r['adj_f1']), 9)} | "
            f"{int(r['perl']):,} | {int(r['rust']):,} | {int(r['intersection']):,} |"
        )
    n = next(iter(walls.values()))[1] if walls else 0
    cells = []
    for s in SNP_INDEL:
        if ("arm64", s) in walls and ("x86_64", s) in walls:
            dataset, assembly = SUITE_LABELS[s]
            label = dataset.replace(" full", f" {assembly}")
            cells.append(
                f"{label} {half_up(walls[('arm64', s)][0], 2)} s ARM / "
                f"{half_up(walls[('x86_64', s)][0], 2)} s x86"
            )
    lines += [
        "",
        f"Wall time, median of {n} independent machines per cell, ARM Graviton4 (`r8gd.8xlarge`)",
        "and x86 Intel (`r8id.8xlarge`), 16 threads, local NVMe, sites-only inputs, the discarded warmup's",
        "output deleted before the timed run: " + "; ".join(cells) + ".",
    ]
    return "\n".join(lines)

def read_sums(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[1].lstrip("*")] = parts[0]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--version", required=True)
    ap.add_argument(
        "--previous", required=True, help="previous release tag, e.g. v0.1.0"
    )
    ap.add_argument("--commit", required=True)
    ap.add_argument("--toolchain", required=True)
    ap.add_argument("--sha256sums", type=Path)
    ap.add_argument("--notes", type=Path, help=".github/release-notes/vX.Y.Z.md")
    ap.add_argument("--changelog-section", help="defaults to the version")
    ap.add_argument("--date", help="release date; defaults to the CHANGELOG heading's")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()

    v, tag = a.version, f"v{a.version}"
    repo = a.repo
    section, cl_date = changelog_section(
        (repo / "CHANGELOG.md").read_text(), a.changelog_section or v
    )
    date = a.date or cl_date
    if not date:
        sys.exit(
            "ERROR: [render_release_notes] no date: pass --date or date the CHANGELOG heading"
        )
    rust_version = re.search(
        r'^rust-version\s*=\s*"([^"]+)"', (repo / "Cargo.toml").read_text(), re.MULTILINE
    )
    if not rust_version:
        sys.exit("ERROR: [render_release_notes] no rust-version in Cargo.toml")
    hand = hand_written(
        a.notes if a.notes else repo / ".github" / "release-notes" / f"{tag}.md"
    )
    if not hand["summary"]:
        sys.exit(
            "ERROR: [render_release_notes] the release-notes file needs a summary paragraph"
        )
    sums = read_sums(a.sha256sums)
    version_doi = read_doi(repo / "CITATION.cff", v)
    dl = f"{REPO_URL}/releases/download/{tag}"
    ph = "`<sha256 from SHA256SUMS>`"

    out = [hand["summary"], "", ABOUT, ""]
    if hand["upgrade"]:
        out += ["## Upgrade notes", "", hand["upgrade"], ""]
    out += [
        "## Changes",
        "",
        section,
        "",
        f"Full entry: [CHANGELOG.md]({REPO_URL}/blob/{tag}/CHANGELOG.md) (this version is its first section) ·",
        f"**Full Changelog**: {REPO_URL}/compare/{a.previous}...{tag}",
        "",
        measurements_section(repo, v, date),
        "",
        "## Install",
        "",
        "| File | Platform | CPU floor | SHA256 |",
        "| ---- | -------- | --------- | ------ |",
    ]
    for triple, platform, floor in TARGETS:
        name = f"vep-{v}-{triple}.tar.gz"
        digest = f"`{sums[name]}`" if name in sums else ph
        out.append(f"| [{name}]({dl}/{name}) | {platform} | {floor} | {digest} |")
    x86 = f"vep-{v}-x86_64-unknown-linux-gnu.tar.gz"
    out += [
        "",
        f"Each archive unpacks to `vep-{v}-<target>/` holding `vep`, `vep-cache-builder` and",
        "`vep-cache-converter`. A binary run on a CPU below its floor aborts with an illegal-instruction",
        "fault rather than running slowly; on such a host build from source, which compiles for the machine",
        'it runs on when `RUSTFLAGS="-C target-cpu=native"` is set.',
        "",
        "```sh",
        f"curl -LO {dl}/{x86}",
        f"curl -LO {dl}/SHA256SUMS",
        "sha256sum -c --ignore-missing SHA256SUMS   # macOS: shasum -a 256 -c --ignore-missing SHA256SUMS",
        f"tar xzf {x86}",
        "```",
        "",
        f"Container image, pushed when this release is published: `docker pull ghcr.io/natera-open-source/vep-rs:{v}`",
        "(linux/amd64: the x86_64 binaries above and the `duckdb` CLI that Parquet output needs).",
        "",
        f"From source (Rust {rust_version.group(1)} or newer): `cargo build --release --locked`. The README covers",
        "the JSON transcript cache every run needs.",
        "",
        "## Verify",
        "",
        f"Built by the [release workflow]({REPO_URL}/actions/workflows/release.yml)",
        f"from tag `{tag}` (commit `{a.commit}`) with Rust {a.toolchain} and `cargo build --release --locked`.",
        "Every file in `SHA256SUMS` carries a GitHub build-provenance attestation:",
        "",
        "```sh",
        f"gh attestation verify {x86} --repo natera-open-source/vep-rs",
        "```",
        "",
        "## Source archive",
        "",
    ]
    src = f"vep-rs-{v}.tar.gz"
    src_digest = f"`{sums[src]}`" if src in sums else ph
    out.append(
        f"`{src}` is `git archive --format=tar.gz --prefix=vep-rs-{v}/ {tag}`, sha256 {src_digest}."
    )
    if version_doi:
        out.append(f"Version DOI [{version_doi}](https://doi.org/{version_doi}).")
    cite_doi = version_doi or CONCEPT_DOI
    out += [
        "",
        "## Citing",
        "",
        f"Porter M, Borkowski R. vep-rs, version {v}. {date[:4]}. doi:{cite_doi}.",
        "`CITATION.cff` carries the citation in machine-readable form; the README's Citing section names the",
        "accompanying paper.",
    ]
    if hand["known"]:
        out += ["", "## Known issues", "", hand["known"]]
    body = "\n".join(out).rstrip("\n") + "\n"
    if a.out:
        a.out.write_text(body)
    else:
        sys.stdout.write(body)


if __name__ == "__main__":
    main()
