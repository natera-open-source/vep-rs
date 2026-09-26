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
import csv
import json
import math
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
        "Linux x86_64 (glibc 2.17 or later)",
        "x86-64-v3: AVX2, BMI2, FMA (Intel Haswell 2013 or later, AMD Zen or later)",
    ),
    (
        "aarch64-unknown-linux-gnu",
        "Linux aarch64 (glibc 2.17 or later)",
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


def speedup(x: float) -> str:
    if x >= 100:
        return f"{half_up(x, 0)}×"
    if x >= 10:
        return f"{half_up(x, 1)}×"
    return f"{half_up(x, 2)}×"


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


def read_dois(cff: Path, version: str) -> tuple[str | None, str | None]:
    """(this version's DOI or None, "X.Y.Z (doi)" of the newest archived version or None).

    A version DOI is read from an `identifiers` entry whose description names `version X.Y.Z`, or
    from the top-level `doi` when the file's `version` is this version and the DOI is not the
    concept DOI; the concept DOI never counts as a version's own."""
    text = cff.read_text()
    this, latest = None, None
    for m in re.finditer(r"value: (10\.5281/zenodo\.\d+)\n\s+description: (.*)", text):
        doi, desc = m.group(1), m.group(2)
        vm = re.search(r"version (\d+\.\d+\.\d+)", desc)
        if vm and doi != CONCEPT_DOI:
            if vm.group(1) == version:
                this = doi
            latest = (vm.group(1), doi)
    top_doi = re.search(r"^doi: (10\.5281/zenodo\.\d+)$", text, re.MULTILINE)
    top_version = re.search(r"^version: (\d+\.\d+\.\d+)$", text, re.MULTILINE)
    if top_doi and top_version and top_doi.group(1) != CONCEPT_DOI:
        if top_version.group(1) == version:
            this = this or top_doi.group(1)
        if latest is None:
            latest = (top_version.group(1), top_doi.group(1))
    return this, (None if latest is None else f"{latest[0]} ({latest[1]})")


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
    its per-suite concordance and per-cell wall-time aggregates, plus the paper's Perl VEP medians
    from manuscript/data/wall_times.csv; the record, not a file under manuscript/data, carries a
    release's own measurements."""
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
    perl = {}
    for row in csv.DictReader(open(repo / "manuscript" / "data" / "wall_times.csv", newline="")):
        if row["engine"] == "perl" and row["run_index"].startswith("median_of"):
            notes = dict(
                kv.split("=", 1) for kv in row["notes"].split(";") if "=" in kv
            )
            perl[(notes["arch"], notes["suite_id"])] = float(row["wall_time_sec"])
    missing = [s for s in SCORED if s not in conc]
    if missing:
        sys.exit(
            f"ERROR: [render_release_notes] {record_path.name} lacks concordance rows for {missing}"
        )
    return conc, walls, perl, record_path.relative_to(repo).as_posix()


def measurements_section(repo: Path, version: str, date: str) -> str:
    conc, walls, perl, record_rel = load_release_measurements(repo, version)
    record_dir = record_rel.rsplit("/", 1)[0]
    lines = [
        "## Concordance and wall time for this release",
        "",
        "The paper's published figures are unchanged; these are the released binary's own, from the",
        f"release's provenance record under [`{record_dir}/`]({REPO_URL}/tree/v{version}/{record_dir})",
        "(method: `scripts/concordance/run_clone_measurement.sh`).",
        "",
        "| Dataset | Assembly | Raw F1 | Adjusted F1 | VEP tuples | vep-rs tuples | Matched |",
        "| ------- | -------- | ------ | ----------- | ---------- | ------------- | ------- |",
    ]
    for s in SCORED:
        r = conc[s]
        dataset, assembly = SUITE_LABELS[s]
        lines.append(
            f"| {dataset} | {assembly} | {float(r['raw_f1']):.6f} | {float(r['adj_f1']):.6f} | "
            f"{int(r['perl']):,} | {int(r['rust']):,} | {int(r['intersection']):,} |"
        )
    n = next(iter(walls.values()))[1] if walls else 0
    cells, geo = [], {}
    for arch, label in (("arm64", "ARM"), ("x86_64", "x86")):
        logs = []
        for s in SNP_INDEL:
            if (arch, s) in walls and (arch, s) in perl:
                logs.append(math.log(perl[(arch, s)] / walls[(arch, s)][0]))
        geo[label] = math.exp(sum(logs) / len(logs)) if logs else float("nan")
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
        "and x86 Intel (`r8id.8xlarge`), 16 threads, local NVMe, sites-only inputs: "
        + "; ".join(cells)
        + ".",
        f"Geomean speedup over Ensembl VEP 115.2 (the paper's medians): {speedup(geo['ARM'])} ARM, {speedup(geo['x86'])} x86.",
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
    version_doi, archived = read_dois(repo / "CITATION.cff", v)
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
        "## Source archive and DOI",
        "",
    ]
    src = f"vep-rs-{v}.tar.gz"
    src_digest = f"`{sums[src]}`" if src in sums else ph
    out.append(
        f"`{src}` is `git archive --format=tar.gz --prefix=vep-rs-{v}/ {tag}`, sha256 {src_digest}."
    )
    if version_doi:
        out.append(f"Version DOI [{version_doi}](https://doi.org/{version_doi}).")
    elif archived:
        av, adoi = archived.split(" (")
        adoi = adoi.rstrip(")")
        out.append(
            f"This version has no Zenodo deposit of its own; the archived version is {av},"
            f" [{adoi}](https://doi.org/{adoi})."
        )
    else:
        out.append("This version has no Zenodo deposit of its own.")
    out.append(f"All versions: [{CONCEPT_DOI}](https://doi.org/{CONCEPT_DOI}).")
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
