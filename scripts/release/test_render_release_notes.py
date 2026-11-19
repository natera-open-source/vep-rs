"""Tests for render_release_notes.py on three inputs: the repository's own 0.3.2 inputs (the release and
population records, the CHANGELOG section, CITATION.cff, Cargo.toml and the published SHA256SUMS digests),
copied into a scratch repository that each test may edit; release and population records of the
six-dataset shape that the tests build under the same version, on the instance classes and clone count a
test names; and the repository's own 0.4.0 inputs, whose rendering is the published 0.4.0 page. The scratch
release-notes file keeps the repository file's Upgrade notes under an opener within the renderer's budget; the
tests of the page by area write their own notes file over the same scratch repository."""

from __future__ import annotations

import copy
import hashlib
import json
import re
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import NamedTuple

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "render_release_notes.py"
ROOT = HERE.parent.parent
VERSION = "0.3.2"
TAG = f"v{VERSION}"
PREVIOUS = "v0.3.1"
COMMIT = "dd55b91219a9bf433608906950acbc29645493f1"
TOOLCHAIN = "1.97.1"
REPO_URL = "https://github.com/natera-open-source/vep-rs"
DOWNLOAD = f"{REPO_URL}/releases/download/{TAG}"
CONCEPT_DOI = "10.5281/zenodo.22837896"
OPENER = (
    "A performance release with no change in output: every gated input annotates byte-identically to 0.3.1, "
    "and the whole of gnomAD v4.1 takes 724 s on ARM against 974 s under 0.3.1. Hosts wider than 32 logical "
    "CPUs on the default `--fork` gain the most from upgrading."
)
DIGESTS = {
    "vep-0.3.2-aarch64-apple-darwin.tar.gz": "3b138436813751ebfc265f72605786dd9069cb51669d0389b9ba0623e8e8729e",
    "vep-0.3.2-aarch64-unknown-linux-gnu.tar.gz": "692c66dbb11af0056bdd1511ec4028aedf462a1b761576ff5c1dfc1d56458217",
    "vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz": "a52df441c299e7690c6732d81da13e54192801f577814b5db282820eee32baf4",
    "vep-rs-0.3.2.tar.gz": "cc61520bc8542e96289218e68af9e57f05fab313edbab9fb6cb13ea3cb2518b9",
}
PLATFORMS = {
    "vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz": (
        "Linux x86_64, glibc 2.39 or later, x86-64-v3 (AVX2, BMI2, FMA: Intel Haswell 2013 or later, AMD Zen or later)"
    ),
    "vep-0.3.2-aarch64-unknown-linux-gnu.tar.gz": (
        "Linux aarch64, glibc 2.39 or later, Armv8.4-A with SVE (Arm Neoverse V1 or later: "
        "AWS Graviton3 and later, not Graviton2, Ampere Altra or Raspberry Pi)"
    ),
    "vep-0.3.2-aarch64-apple-darwin.tar.gz": "macOS 11 or later on Apple silicon, any Apple M-series",
}
SUITES = (
    ("s01", "ClinVar full", "GRCh37"),
    ("s02", "ClinVar full", "GRCh38"),
    ("s03", "gnomAD v2.1.1 chr21", "GRCh37"),
    ("s04", "gnomAD v4.1 chr21", "GRCh38"),
    ("s05", "1KG Phase 3 chr21", "GRCh37"),
    ("s06", "1KG high-cov chr21", "GRCh38"),
    ("s07", "SV per-VCF (16 files)", "GRCh37"),
    ("s08", "SV per-VCF (16 files)", "GRCh38"),
)


def repository_upgrade_notes() -> str:
    text = (ROOT / ".github" / "release-notes" / f"{TAG}.md").read_text()
    return text.split("## Upgrade notes\n", 1)[1].strip()


def fixture_repo(tmp_path: Path) -> Path:
    repo = tmp_path / "repo"
    (repo / "docs" / "concordance-provenance").mkdir(parents=True)
    (repo / ".github" / "release-notes").mkdir(parents=True)
    for rel in ("CHANGELOG.md", "CITATION.cff", "Cargo.toml"):
        shutil.copy(ROOT / rel, repo / rel)
    # The repository's CITATION.cff names the current version, not the fixture release, and `read_doi`
    # reads the top-level DOI as a version's own only when the file's version is the one rendered.
    cff = repo / "CITATION.cff"
    cff.write_text(re.sub(r"^version: .*$", f"version: {VERSION}", cff.read_text(), flags=re.MULTILINE))
    for record in (ROOT / "docs" / "concordance-provenance").glob(f"*-v{VERSION}.json"):
        shutil.copy(record, repo / "docs" / "concordance-provenance" / record.name)
    write_notes(repo, OPENER, f"Upgrade notes\n\n{repository_upgrade_notes()}")
    (repo / "SHA256SUMS").write_text("".join(f"{d}  {n}\n" for n, d in DIGESTS.items()))
    write_floors(repo, {n: "2.39" for n in DIGESTS if "linux-gnu" in n})
    return repo


def write_floors(repo: Path, floors: dict[str, str]) -> None:
    """The measured glibc floor per Linux archive, the file check_glibc_floor.py --json writes."""
    (repo / "glibc-floors.json").write_text(json.dumps(floors))


def render(repo: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--repo", str(repo), "--version", VERSION, "--previous", PREVIOUS,
         "--commit", COMMIT, "--toolchain", TOOLCHAIN, "--sha256sums", str(repo / "SHA256SUMS"),
         "--glibc-floors", str(repo / "glibc-floors.json"), *extra],
        check=False, capture_output=True, text=True,
    )


def page(repo: Path, *extra: str) -> str:
    out = render(repo, *extra)
    assert out.returncode == 0, out.stderr
    return out.stdout


def error(repo: Path, *extra: str) -> str:
    out = render(repo, *extra)
    assert out.returncode == 1 and out.stderr.startswith("ERROR: [render_release_notes]"), out.stderr
    return out.stderr


def notes(repo: Path) -> Path:
    return repo / ".github" / "release-notes" / f"{TAG}.md"


def write_notes(repo: Path, opener: str = "A patch release.", *sections: str) -> None:
    notes(repo).write_text(f"# vep-rs {VERSION}\n\n{opener}\n" + "".join(f"\n## {s}\n" for s in sections))


def write_changelog(repo: Path, section: str) -> None:
    (repo / "CHANGELOG.md").write_text(
        f"# Changelog\n\n## [Unreleased]\n\n### Changed\n\n- Unreleased work.\n\n## [{VERSION}] - 2026-10-06\n\n{section}\n"
    )


def edit_record(repo: Path, kind: str, edit) -> None:
    path = next((repo / "docs" / "concordance-provenance").glob(f"*-{kind}-v{VERSION}.json"))
    record = json.loads(path.read_text())
    edit(record)
    path.write_text(json.dumps(record))


def half_up(q: Fraction, places: int) -> str:
    """ROUND_HALF_UP of a non-negative rational to `places` decimals, with the fractions module rather than
    the renderer's Decimal arithmetic."""
    scaled = q * 10**places
    n = scaled.numerator // scaled.denominator
    if scaled - n >= Fraction(1, 2):
        n += 1
    return f"{n // 10**places}.{n % 10**places:0{places}d}" if places else f"{n:,}"


def f1(intersection: int, a: int, b: int) -> str:
    return half_up(Fraction(2 * intersection, a + b), 9)


def release_record(repo: Path = ROOT) -> dict:
    return json.loads(next((repo / "docs" / "concordance-provenance").glob(f"*-release-v{VERSION}.json")).read_text())


def population_record(repo: Path = ROOT) -> dict:
    return json.loads(next((repo / "docs" / "concordance-provenance").glob(f"*-population-v{VERSION}.json")).read_text())


def declared_msrv(repo: Path) -> str:
    return re.search(r'^rust-version = "([^"]+)"', (repo / "Cargo.toml").read_text(), re.MULTILINE).group(1)


def expected_whole_genome_rows(repo: Path = ROOT) -> list[str]:
    pop = release_record(repo)["population_concordance"]
    cells = population_record(repo)["wall_time"]["vep_rs"]["per_cell"]
    sums: dict[tuple[str, str], Fraction] = {}
    for c in cells:
        key = (c["dataset"], c["arch"])
        sums[key] = sums.get(key, Fraction(0)) + Fraction(str(c["median_sec"]))

    def row(label: str, r: dict, arm: Fraction, x86: Fraction) -> str:
        raw = f1(r["intersection"], r["perl_tuples"], r["vep_rs_tuples"])
        adjusted = f1(r["intersection"], r["perl_tuples"] - r["excluded_perl"], r["vep_rs_tuples"] - r["excluded_rust"])
        return f"| {label} | {raw} | {adjusted} | {half_up(arm, 0)} | {half_up(x86, 0)} |"

    rows = [
        row(r["label"].replace(", whole genome", ""), r, sums[(r["dataset"], "arm64")], sums[(r["dataset"], "x86_64")])
        for r in pop["per_dataset"]
    ]
    arm = sum((sums[(r["dataset"], "arm64")] for r in pop["per_dataset"]), Fraction(0))
    x86 = sum((sums[(r["dataset"], "x86_64")] for r in pop["per_dataset"]), Fraction(0))
    rows.append(row("All four", pop["pooled"], arm, x86))
    return rows


def expected_chromosome_21_rows() -> list[str]:
    conc = release_record()["concordance"]
    return [
        f"| {dataset} | {assembly} | {f1(conc[s]['intersection'], conc[s]['perl'], conc[s]['rust'])} | 1.000000000 |"
        for s, dataset, assembly in SUITES
    ]


def changelog_bullets(repo: Path, version: str) -> list[str]:
    """The version's CHANGELOG bullets with each one's wrapped lines joined by single spaces."""
    text = (repo / "CHANGELOG.md").read_text()
    section = re.split(r"^## \[", text, flags=re.MULTILINE)
    body = next(s for s in section if s.startswith(f"{version}]"))
    body = body.split("\n", 1)[1].replace("### Changed\n", "")
    return ["- " + " ".join(b.split()) for b in body.split("\n- ")[1:]]


def expected_page(repo: Path) -> str:
    msrv = declared_msrv(repo)
    install = [f"- [`{n}`]({DOWNLOAD}/{n}): {PLATFORMS[n]}; sha256 `{DIGESTS[n]}`" for n in PLATFORMS]
    install += [
        f"- Each archive unpacks to `vep-{VERSION}-<target>/` holding `vep`, `vep-cache-builder` and `vep-cache-converter`; "
        "a binary below its CPU floor aborts with an illegal-instruction fault, so on such a host build from source "
        f"(Rust {msrv} or later), which compiles for the machine it runs on: "
        '`RUSTFLAGS="-C target-cpu=native" cargo build --release --locked`',
        f"- Container image: `docker pull ghcr.io/natera-open-source/vep-rs:{VERSION}` (linux/amd64: the x86_64 binaries "
        "with `bgzip`, `tabix` and `ps`)",
        f"- Built by the [release workflow]({REPO_URL}/actions/workflows/release.yml) from tag `{TAG}` (commit `{COMMIT}`) "
        f"with Rust {TOOLCHAIN}; every file in [`SHA256SUMS`]({DOWNLOAD}/SHA256SUMS) carries a GitHub build-provenance "
        "attestation: `gh attestation verify <file> --repo natera-open-source/vep-rs`",
        f"- [`vep-rs-{VERSION}.tar.gz`]({DOWNLOAD}/vep-rs-{VERSION}.tar.gz) is `git archive --format=tar.gz "
        f"--prefix=vep-rs-{VERSION}/ {TAG}`; sha256 `{DIGESTS[f'vep-rs-{VERSION}.tar.gz']}`",
    ]
    blocks = [
        OPENER,
        "## Upgrade notes",
        repository_upgrade_notes(),
        "## Changed",
        "\n".join(changelog_bullets(repo, VERSION)),
        "## Concordance and wall time",
        "Whole genome against Ensembl VEP 115.2, one run per chromosome; wall time is the sum of the per-chromosome "
        "medians of 20 machines at 16 threads on `c8gd.8xlarge` (ARM) and `c8id.8xlarge` (x86); peak memory at most "
        "2.5 GB on any chromosome; tuple counts and per-chromosome rows are in the "
        f"[provenance records]({REPO_URL}/tree/{TAG}/docs/concordance-provenance).",
        "\n".join(["| Dataset | Raw F1 | Adjusted F1 | ARM (s) | x86 (s) |", "| --- | --- | --- | --- | --- |", *expected_whole_genome_rows()]),
        "The paper's chromosome 21 cells measured on the released binary; the paper's own figures stay as published in "
        f"[docs/published-figures.md]({REPO_URL}/blob/{TAG}/docs/published-figures.md).",
        "\n".join(["| Dataset | Assembly | Raw F1 | Adjusted F1 |", "| --- | --- | --- | --- |", *expected_chromosome_21_rows()]),
        "## Install",
        "\n".join(install),
        "## Cite",
        f"Porter M, Borkowski R. vep-rs, version {VERSION}. 2026. doi:[{CONCEPT_DOI}](https://doi.org/{CONCEPT_DOI})",
        f"**Full changelog**: [CHANGELOG]({REPO_URL}/blob/{TAG}/CHANGELOG.md) · [`{PREVIOUS}...{TAG}`]({REPO_URL}/compare/{PREVIOUS}...{TAG})",
    ]
    return "\n\n".join(blocks) + "\n"


def test_the_0_3_2_page_renders_from_the_repository_files(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    assert len(OPENER.split()) < 50 and OPENER.count(". ") == 1
    assert page(repo) == expected_page(repo)


def test_the_page_runs_about_sixty_lines_and_is_refused_over_eight_hundred(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    lines = page(repo).splitlines()
    assert len(lines) <= 60
    assert [l for l in lines if l.startswith("#")] == [
        "## Upgrade notes", "## Changed", "## Concordance and wall time", "## Install", "## Cite",
    ]
    assert not [l for l in lines if l.startswith((">", "<"))]
    assert "**Full changelog**" not in "\n".join(lines[:-1])
    assert lines[-1].startswith("**Full changelog**: ")

    def changed_bullets(n: int) -> None:
        write_changelog(repo, "### Changed\n\n" + "".join(f"- Change number {i}.\n" for i in range(n)))

    lines_outside_changed = len(lines) - len(changelog_bullets(repo, VERSION))
    changed_bullets(800 - lines_outside_changed)
    assert page(repo).count("\n") == 800
    changed_bullets(801 - lines_outside_changed)
    assert error(repo) == f"{ERR}the page is 801 lines; the cap is 800\n"


def test_the_page_renders_at_120_000_characters_and_is_refused_one_character_over(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)

    def one_change_of(digits: int) -> None:
        write_changelog(repo, f"### Changed\n\n- Change {'0' * digits}.\n")

    one_change_of(1)
    characters_outside_the_digits = len(page(repo)) - 1
    one_change_of(120_000 - characters_outside_the_digits)
    at_the_cap = page(repo)
    assert len(at_the_cap) == 120_000
    # The footer's middle dot is two bytes of UTF-8, so a count of bytes would refuse this page.
    assert len(at_the_cap.encode()) > 120_000
    one_change_of(120_001 - characters_outside_the_digits)
    assert error(repo) == (
        f"{ERR}the page is 120,001 characters; the cap is 120,000, "
        "and GitHub refuses a release body over 125,000 characters\n"
    )


def test_whole_genome_table_derives_every_cell_from_the_records(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    body = page(repo)
    rows = expected_whole_genome_rows()
    assert "| Dataset | Raw F1 | Adjusted F1 | ARM (s) | x86 (s) |\n| --- | --- | --- | --- | --- |\n" + "\n".join(rows) in body
    assert rows[-1] == "| All four | 0.999998884 | 1.000000000 | 1,121 | 1,818 |"
    cells = [c.strip() for l in body.splitlines() if l.startswith("| ") for c in l.split("|")]
    assert cells.count("1.000000000") == 14
    assert "whole genome |" not in body and "VEP tuples" not in body and "Chromosomes" not in body

    def spoil(record: dict) -> None:
        for r in [*record["population_concordance"]["per_dataset"], record["population_concordance"]["pooled"]]:
            r["raw_f1"], r["adjusted_f1"] = 0.5, 0.5
        record["concordance"]["s01"]["raw_f1"] = 0.5

    edit_record(repo, "release", spoil)
    assert page(repo) == body

    def unmask(record: dict) -> None:
        record["population_concordance"]["per_dataset"][3]["excluded_rust"] = 0

    edit_record(repo, "release", unmask)
    r = release_record(repo)["population_concordance"]["per_dataset"][3]
    adjusted = f1(r["intersection"], r["perl_tuples"] - r["excluded_perl"], r["vep_rs_tuples"])
    assert adjusted != "1.000000000"
    assert f"| 1000 Genomes high-coverage | {rows[3].split(' | ')[1]} | {adjusted} | " in page(repo)


def test_chromosome_21_table_covers_the_eight_scored_suites_in_suite_order(tmp_path: Path) -> None:
    body = page(fixture_repo(tmp_path))
    rows = expected_chromosome_21_rows()
    assert "| Dataset | Assembly | Raw F1 | Adjusted F1 |\n| --- | --- | --- | --- |\n" + "\n".join(rows) + "\n" in body
    assert rows[0].startswith("| ClinVar full | GRCh37 | 0.999979065 |")
    assert rows[-1].startswith("| SV per-VCF (16 files) | GRCh38 | 0.902664038 |")
    assert body.index("| All four |") < body.index("published-figures.md") < body.index("| ClinVar full | GRCh37 |")
    assert "Wall time, median of" not in body and "r8gd" not in body


def test_captions_carry_the_comparator_protocol_memory_and_record_link(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    body = page(repo)
    assert (
        "Whole genome against Ensembl VEP 115.2, one run per chromosome; wall time is the sum of the per-chromosome "
        "medians of 20 machines at 16 threads on `c8gd.8xlarge` (ARM) and `c8id.8xlarge` (x86); peak memory at most "
        "2.5 GB on any chromosome; tuple counts and per-chromosome rows are in the "
        f"[provenance records]({REPO_URL}/tree/{TAG}/docs/concordance-provenance)."
    ) in body
    assert "-release-v" not in body and "-population-v" not in body

    def set_peak(mb: float):
        def edit(record: dict) -> None:
            for cell in record["wall_time"]["vep_rs"]["per_cell"]:
                cell["peak_rss_mb_median"] = 1_900.0
            record["wall_time"]["vep_rs"]["per_cell"][0]["peak_rss_mb_median"] = mb

        return edit

    # The bound rounds up, never down: 2,500.0 MB is exactly 2.5 GB, anything above it reads 2.6 GB.
    for mb, shown in ((2_500.0, "2.5"), (2_500.1, "2.6"), (2_502.7, "2.6"), (2_438.3, "2.5")):
        edit_record(repo, "population", set_peak(mb))
        assert f"peak memory at most {shown} GB on any chromosome" in page(repo), mb

    # The earlier page labels a pair of two generations as the six-dataset page does.
    edit_record(repo, "population", lambda r: r["wall_time"]["vep_rs"].update(instance_types=NEWEST, n_clones_per_cell=5))
    assert f"medians of 5 machines at 16 threads on {pair(NEWEST)}, {NEWEST_LABEL}; peak memory at most" in page(repo)


def test_the_concordance_caption_reads_the_perl_reference_block_when_a_record_carries_one(tmp_path: Path) -> None:
    """A release scored against a later Ensembl VEP names it in a `perl_reference` block (image, digest,
    cache_version, prefix); the caption follows the block and prints no trace of the earlier release.
    The wall-time columns stay the engine's own `engine_version`, so a comparator block carried forward
    from an earlier record keeps its version there and nowhere else."""
    repo = fixture_repo(tmp_path)
    block = {
        "image": "ensemblorg/ensembl-vep:release_116.2",
        "digest": "sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da",
        "cache_version": 116,
        "prefix": "population_r116.2_c116_fixture",
    }

    def name_the_reference(record: dict) -> None:
        record["perl_reference"] = block

    edit_record(repo, "release", name_the_reference)
    body = page(repo)
    assert "Whole genome against Ensembl VEP 116.2, one run per chromosome;" in body
    assert "115.2" not in body
    assert "5c57abdc" not in body and "population_r116.2" not in body

    # The block under population_concordance, and in the population record, are read the same way;
    # the block wins over the population record's ensembl_vep engine version.
    def move_under_population_concordance(record: dict) -> None:
        del record["perl_reference"]
        record["population_concordance"]["perl_reference"] = dict(block, image="ensemblorg/ensembl-vep@sha256:abc")

    edit_record(repo, "release", move_under_population_concordance)
    assert "Whole genome against Ensembl VEP 116.2," in page(repo)

    def prefix_only(record: dict) -> None:
        record["population_concordance"]["perl_reference"] = {"prefix": "population_r116.2_c116_fixture", "cache_version": 116}

    edit_record(repo, "release", prefix_only)
    assert "Whole genome against Ensembl VEP 116.2," in page(repo)

    def cache_version_only(record: dict) -> None:
        record["population_concordance"]["perl_reference"] = {"cache_version": 116}

    edit_record(repo, "release", cache_version_only)
    assert "Whole genome against Ensembl VEP 116.2," in page(repo), "cache 116 is release 116.2, as the comparators read it"

    def empty_block(record: dict) -> None:
        record["population_concordance"]["perl_reference"] = {"digest": "sha256:abc"}

    edit_record(repo, "release", empty_block)
    assert "a perl_reference block names no release" in error(repo)

    # Without any block the caption still follows the population record's engine version.
    def drop_block(record: dict) -> None:
        del record["population_concordance"]["perl_reference"]

    edit_record(repo, "release", drop_block)
    assert "Whole genome against Ensembl VEP 115.2," in page(repo)


def test_the_opener_is_one_paragraph_of_two_sentences_within_fifty_words(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    shutil.copy(ROOT / ".github" / "release-notes" / f"{TAG}.md", notes(repo))
    assert re.search(r"the opener in .* is 1\d\d words; at most 50", error(repo))
    write_notes(repo, "One. Two. Three.")
    assert "is 3 sentences; at most two" in error(repo)
    write_notes(repo, "Drop-in for 0.3.1: flags, formats and the cache format are unchanged. Upgrade freely!")
    assert page(repo).startswith("Drop-in for 0.3.1: flags, formats and the cache format are unchanged. Upgrade freely!\n\n## Changed\n")
    write_notes(repo, "One paragraph.\n\nA second paragraph.")
    assert "is more than one paragraph" in error(repo)
    write_notes(repo, "", "Upgrade notes\n\n- A note.")
    assert "needs an opener paragraph" in error(repo)
    notes(repo).unlink()
    assert "needs an opener paragraph" in error(repo)


def test_highlights_is_two_to_four_bullets_over_at_least_three_changes(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_notes(repo, "A release.", "Highlights\n\n- One.\n- Two.\n- Three.", "Upgrade notes\n\n- A note.")
    body = page(repo)
    assert "A release.\n\n## Highlights\n\n- One.\n- Two.\n- Three.\n\n## Upgrade notes\n\n- A note.\n\n## Changed\n" in body
    write_notes(repo, "A release.", "Highlights\n\n* One wrapped\n  over two lines.\n* Two.")
    assert "## Highlights\n\n* One wrapped over two lines.\n* Two.\n\n## Changed\n" in page(repo)
    write_notes(repo, "A release.", "Highlights\n\n- One.")
    assert "has 1 bullets; two to four" in error(repo)
    write_notes(repo, "A release.", "Highlights\n\n- One.\n- Two.\n- Three.\n- Four.\n- Five.")
    assert "has 5 bullets; two to four" in error(repo)
    write_notes(repo, "A release.", "Highlights\n\nProse, not a bullet.\n- Two.")
    assert "'## Highlights' in" in error(repo) and "is bullets only" in error(repo)
    write_notes(repo, "A release.", "Highlights\n\n- One.\n- Two.")
    write_changelog(repo, "### Fixed\n\n- A fix.\n- Another fix.")
    assert "needs at least three CHANGELOG bullets; the section has 2" in error(repo)
    write_changelog(repo, "### Fixed\n\n- A fix.\n- Another fix.\n\n### Added\n\n- A thing.")
    assert "## Highlights\n\n- One.\n- Two.\n\n## Fixed\n" in page(repo)


def test_upgrade_notes_are_bullets_with_breaking_ones_first(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_notes(repo, "A release with nothing to upgrade.")
    body = page(repo)
    assert "## Upgrade notes" not in body
    assert body.startswith("A release with nothing to upgrade.\n\n## Changed\n\n- ")
    write_notes(repo, "A release.", "Upgrade notes\n\n- **Breaking:** `--old` is gone; pass `--new`.\n* A cache converted\n  for 0.3.0 is valid.")
    assert "A release.\n\n## Upgrade notes\n\n- **Breaking:** `--old` is gone; pass `--new`.\n* A cache converted for 0.3.0 is valid.\n\n## Changed\n" in page(repo)
    write_notes(repo, "A release.", "Upgrade notes\n\nNothing breaks.\n- A note.")
    assert "'## Upgrade notes' in" in error(repo) and "is bullets only" in error(repo)
    # No count bounds the section: six bullets render, as would sixty.
    six_bullets = "".join(f"- Note {i}.\n" for i in range(6))
    write_notes(repo, "A release.", f"Upgrade notes\n\n{six_bullets}")
    assert f"## Upgrade notes\n\n{six_bullets}\n## Changed\n" in page(repo)
    write_notes(repo, "A release.", "Upgrade notes\n\n- A note.\n- **Breaking:** `--old` is gone.")
    assert "'## Upgrade notes' in" in error(repo) and "lists a Breaking: bullet after one that is not" in error(repo)
    write_notes(repo, "A release.", "Upgrade notes\n\n- A note.\n- Breaking: `--old` is gone.")
    assert "lists a Breaking: bullet after one that is not" in error(repo)
    write_notes(repo, "A release.", "Upgrade notes\n\n- `Breaking:` `--old` is gone.")
    assert error(repo) == f"{ERR}'## Upgrade notes' in {notes(repo)} marks a breaking change as '`Breaking:`'; the mark is Breaking:\n"
    # The mark is read after the bullet's own spacing, as a change under an area is.
    write_notes(repo, "A release.", "Upgrade notes\n\n- A note.\n-  **Breaking:** `--old` is gone.")
    assert error(repo) == f"{ERR}'## Upgrade notes' in {notes(repo)} lists a Breaking: bullet after one that is not\n"


AREAS = """\
- Input and formats
  - Breaking: standard input needs `--format`; a run without it exits 1 with Ensembl VEP's message.
  - Plain gzip and `.bgz` inputs are read.
  - `--format guess` detects the format from the first data line.
- Output
  - Breaking: `<output>_warnings.txt` is written beside the output; `--warning_file STDERR` keeps the old stream.
  - **Breaking:** `AF` prints four decimals.
  - `MAX_AF_POPS` lists every tied population.
- Caches
  - Converted caches carry the transcript attributes and the variation files.
  - Breaking the input into batches of 5000 records stays the default."""
SECURITY = """\
- Malformed inputs
  - A FASTA index whose line length is zero is refused with `ERROR: [fasta]`.
  - A structural-variant allele of an unknown type is skipped with a warning.
- Every input line is bounded in length."""
DEPENDENCIES = "- `rustls` 0.23.46 clears four advisories, none reachable from `vep`.\n- The container image carries `bgzip` and `tabix` from the base distribution."
REMOVED = "- `--transcript_index_impl` and the index implementations it selected.\n- The unused second VCF parser."
ONE_AREA = "What changed, by area\n\n- Output\n  - A change."


def test_a_page_by_area_renders_its_sections_as_written_in_place_of_the_kind_lists(tmp_path: Path) -> None:
    """A notes file with `## What changed, by area` renders that section and its three companions byte for
    byte, nested bullets and inline `Breaking:` marks included, where the flat page carries Upgrade notes and
    the CHANGELOG kinds; the rest of the page is the flat page's."""
    repo = fixture_repo(tmp_path)
    write_notes(
        repo, "A release.", "Highlights\n\n- One.\n- Two.", f"What changed, by area\n\n{AREAS}",
        f"Security and robustness\n\n{SECURITY}", f"Dependencies and tooling\n\n{DEPENDENCIES}", f"Removed\n\n{REMOVED}",
    )
    body = page(repo)
    assert body.startswith(
        f"A release.\n\n## Highlights\n\n- One.\n- Two.\n\n## What changed, by area\n\n{AREAS}\n\n"
        f"## Security and robustness\n\n{SECURITY}\n\n## Dependencies and tooling\n\n{DEPENDENCIES}\n\n"
        f"## Removed\n\n{REMOVED}\n\n## Concordance and wall time\n\n"
    )
    assert [l for l in body.splitlines() if l.startswith("#")] == [
        "## Highlights", "## What changed, by area", "## Security and robustness", "## Dependencies and tooling",
        "## Removed", "## Concordance and wall time", "## Install", "## Cite",
    ]
    assert "## Changed" not in body and "## Upgrade notes" not in body and "###" not in body
    flat = page(fixture_repo(tmp_path / "flat"))
    assert body.split("## Concordance and wall time", 1)[1] == flat.split("## Concordance and wall time", 1)[1]
    # The sections render in the page's order whatever the file's, the companions are each optional, an area
    # name may run to AREA_NAME_CAP characters, and a change wrapped over two lines is one line on the page.
    name = "A" * 60
    write_notes(repo, "A release.", f"Removed\n\n{REMOVED}", f"What changed, by area\n\n- {name}\n  - A change wrapped\n    over two lines.\n  - Another.")
    assert (
        f"A release.\n\n## What changed, by area\n\n- {name}\n  - A change wrapped over two lines.\n  - Another.\n\n"
        f"## Removed\n\n{REMOVED}\n\n## Concordance and wall time\n\n"
    ) in page(repo)


@pytest.mark.parametrize(
    ("sections", "expected"),
    [
        (("What changed, by area\n\n- Output\n  - A change.\n    - A detail.",), "'## What changed, by area' in {notes} nests a bullet three levels deep; two is the most"),
        (("What changed, by area\n\n- Output\n- Caches\n  - A change.",), "the area 'Output' in {notes} has no bullets"),
        (("What changed, by area\n\n- Output\n  - A change.\n  - Breaking: a late one.",), "the area 'Output' in {notes} lists a Breaking: bullet after one that is not"),
        (("What changed, by area\n\n- Output\n  - A change.\n  - **Breaking:** a late one.",), "the area 'Output' in {notes} lists a Breaking: bullet after one that is not"),
        ((f"What changed, by area\n\n- {'A' * 61}\n  - A change.",), f"the area '{'A' * 61}' in {{notes}} is 61 characters; at most 60"),
        (("What changed, by area\n\n  - A change with no area.",), "'## What changed, by area' in {notes} opens with a nested bullet; a nested bullet sits under one at the margin"),
        (("What changed, by area\n\n- Output\n - A change.",), "'## What changed, by area' in {notes} indents a bullet by 1; a nested bullet is indented two spaces"),
        (("What changed, by area\n\nProse first.\n- Output\n  - A change.",), "'## What changed, by area' in {notes} is bullets only, each a '-' bullet"),
        (("What changed, by area\n\n* Output\n  - A change.",), "'## What changed, by area' in {notes} is bullets only, each a '-' bullet"),
        (("What changed, by area\n\n- Output\n  * A change.",), "'## What changed, by area' in {notes} is bullets only, each a '-' bullet"),
        (("What changed, by area\n",), "'## What changed, by area' in {notes} has nothing under its header"),
        (("What changed, by area\n\n- \n  - A change.",), "'## What changed, by area' in {notes} has a bullet with no text"),
        (("What changed, by area\n\n- Output\n  - \n  - A change.",), "'## What changed, by area' in {notes} has a bullet with no text"),
        (("What changed, by area\n\n- Output\n  -\n  - A change.",), "'## What changed, by area' in {notes} has a bullet with no text"),
        (("What changed, by area\n\n- Output\n  - A change.\n  - ",), "'## What changed, by area' in {notes} has a bullet with no text"),
        (("What changed, by area\n\n- Output\n  - A change.\n  - BREAKING: a late one.",), "the area 'Output' in {notes} marks a breaking change as 'BREAKING:'; the mark is Breaking:"),
        (("What changed, by area\n\n- Output\n  - **Breaking**: a first one.",), "the area 'Output' in {notes} marks a breaking change as '**Breaking**:'; the mark is Breaking:"),
        (("What changed, by area\n\n- Output\n  - Breaking : a first one.",), "the area 'Output' in {notes} marks a breaking change as 'Breaking :'; the mark is Breaking:"),
        (("What changed, by area\n\n- Output\n  - A change.\n  - Breaking change: a late one.",), "the area 'Output' in {notes} marks a breaking change as 'Breaking change:'; the mark is Breaking:"),
        (("What changed, by area\n\n- Output\n  - A change.\n  - [Breaking]: a late one.",), "the area 'Output' in {notes} marks a breaking change as '[Breaking]:'; the mark is Breaking:"),
        (("What changed, by area\n\n- Output\n  - A change.\n  - ~~Breaking:~~ a late one.",), "the area 'Output' in {notes} marks a breaking change as '~~Breaking:~~'; the mark is Breaking:"),
        (("What changed, by area\n\n- Breaking: `AF` prints four decimals.\n  - A change.",), "the area 'Breaking: `AF` prints four decimals.' in {notes} carries a breaking mark; the mark goes on a change under it"),
        (("What changed, by area\n\n- BREAKING: output\n  - A change.",), "the area 'BREAKING: output' in {notes} carries a breaking mark; the mark goes on a change under it"),
        (("What changed, by area\n\n- Output\n  - **Breaking: an unclosed bold mark.",), "the area 'Output' in {notes} marks a breaking change as '**Breaking:'; the mark is Breaking:"),
        (("What changed, by area\n\n- Output\n  - Breaking:** an odd close.",), "the area 'Output' in {notes} marks a breaking change as 'Breaking:**'; the mark is Breaking:"),
        (("Removed\n\n- `--gone`.",), "'## Removed' in {notes} needs '## What changed, by area'; without it the page lists the CHANGELOG kinds"),
        (("Dependencies and tooling\n\n- A bump.", "Upgrade notes\n\n- A note."), "'## Dependencies and tooling' in {notes} needs '## What changed, by area'; without it the page lists the CHANGELOG kinds"),
        (("Highlights\n", "Upgrade notes\n\n- A note."), "'## Highlights' in {notes} has nothing under its header"),
        ((ONE_AREA, "Removed\n"), "'## Removed' in {notes} has nothing under its header"),
        ((ONE_AREA, "Upgrade notes\n\n- **Breaking:** a note."), "'## Upgrade notes' and '## What changed, by area' in {notes} are two shapes of one page"),
        ((ONE_AREA, "Removed\n\n- A.\n  - B.\n    - C."), "'## Removed' in {notes} nests a bullet three levels deep; two is the most"),
        ((ONE_AREA, "Security and robustness\n\n- A note.\n- Breaking: a late one."), "'## Security and robustness' in {notes} lists a Breaking: bullet after one that is not"),
        ((ONE_AREA, "Removed\n\n- A group.\n  - A plain one.\n  - Breaking: a late one."), "'## Removed' under 'A group.' in {notes} lists a Breaking: bullet after one that is not"),
        ((ONE_AREA, "Dependencies and tooling\n\nProse.\n- A bump."), "'## Dependencies and tooling' in {notes} is bullets only, each a '-' bullet"),
        ((ONE_AREA, "Removed\n\n- `--gone`.\n- breaking: a late one."), "'## Removed' in {notes} marks a breaking change as 'breaking:'; the mark is Breaking:"),
        ((ONE_AREA, ONE_AREA), "'## What changed, by area' appears twice in {notes}"),
    ],
)
def test_a_page_by_area_is_refused_for_each_fault_by_its_line(tmp_path: Path, sections: tuple[str, ...], expected: str) -> None:
    repo = fixture_repo(tmp_path)
    write_notes(repo, "A release.", *sections)
    assert error(repo) == f"{ERR}{expected.format(notes=notes(repo))}\n"


def test_a_nested_changelog_keeps_every_bullet_on_the_flat_page_and_counts_its_changes(tmp_path: Path) -> None:
    """A CHANGELOG kind whose bullets are grouped by area (`- <Area>` with `  - ` entries) renders every line
    on the flat page, each entry's wrapped lines joined; the Highlights rule counts the entries, not the groups."""
    repo = fixture_repo(tmp_path)
    write_notes(repo, "A release.", "Highlights\n\n- One.\n- Two.")
    write_changelog(
        repo,
        "### Changed\n\n- Input and formats\n  - `--format guess` detects the format\n    from the first data line. (#N)\n"
        "  - Plain gzip is read. (#N)\n  - `.bgz` is read. (#N)\n\n### Removed\n\n- `--gone`. (#N)",
    )
    assert (
        "## Changed\n\n- Input and formats\n  - `--format guess` detects the format from the first data line. (#N)\n"
        "  - Plain gzip is read. (#N)\n  - `.bgz` is read. (#N)\n\n## Removed\n\n- `--gone`. (#N)\n\n## Concordance"
    ) in page(repo)
    write_changelog(repo, "### Changed\n\n- Input and formats\n  - Plain gzip is read. (#N)\n- Output\n  - `AF` prints four decimals. (#N)")
    assert "needs at least three CHANGELOG bullets; the section has 2" in error(repo)
    # A list-looking line inside a fenced block under a bullet is code, not a change.
    write_changelog(repo, "### Changed\n\n- One:\n  ```\n  - a\n  - b\n  ```")
    assert error(repo) == f"{ERR}'## Highlights' needs at least three CHANGELOG bullets; the section has 1\n"


def test_kind_headers_rise_to_level_two_in_file_order_and_an_empty_kind_is_absent(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_changelog(repo, "### Removed\n\n- `--gone`.\n\n### Fixed\n\n### Added\n\n- `--new`.\n\n### Deprecated\n\n- `--old`.")
    body = page(repo)
    assert "## Removed\n\n- `--gone`.\n\n## Added\n\n- `--new`.\n\n## Deprecated\n\n- `--old`.\n\n## Concordance" in body
    assert "Fixed" not in body and "###" not in body


def test_each_bullet_is_one_physical_line_in_the_changelog_words(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_changelog(
        repo,
        "### Changed\n\n- A bullet wrapped over\n  three lines of the\n  changelog file.\n"
        "- A bullet with a nested item:\n  - the nested item stays\n    its own line\n"
        "- A fenced block right under its bullet:\n  ```\n  kept\n  as is\n  ```\n"
        "- A tilde fence too:\n  ~~~\n  kept\n  ~~~\n",
    )
    body = page(repo)
    assert (
        "## Changed\n\n- A bullet wrapped over three lines of the changelog file.\n"
        "- A bullet with a nested item:\n  - the nested item stays its own line\n"
        "- A fenced block right under its bullet:\n  ```\n  kept\n  as is\n  ```\n"
        "- A tilde fence too:\n  ~~~\n  kept\n  ~~~\n\n## Concordance" in body
    )
    original = fixture_repo(tmp_path / "original")
    bullets = [l for l in page(original).split("## Concordance")[0].splitlines() if l.startswith("- ")]
    assert len(bullets) == 5 and bullets[2:] == changelog_bullets(original, VERSION)


def test_install_lists_every_archive_in_target_order_with_its_digest(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    body = page(repo)
    install = body.split("## Install\n\n")[1].split("\n\n## Cite")[0].splitlines()
    assert len(install) == 7 and all(l.startswith("- ") for l in install)
    assert [l.split("`")[1] for l in install[:3]] == list(PLATFORMS)
    assert [l.rsplit("`", 2)[1] for l in install[:3]] == [DIGESTS[n] for n in PLATFORMS]
    msrv = declared_msrv(repo)
    assert install[3].startswith("- Each archive unpacks to `vep-0.3.2-<target>/`") and f"Rust {msrv} or later" in install[3]
    assert [l.split(": ", 1)[1].rsplit("; sha256", 1)[0] for l in install[:3]] == [PLATFORMS[n] for n in PLATFORMS]
    assert install[4] == (
        "- Container image: `docker pull ghcr.io/natera-open-source/vep-rs:0.3.2` (linux/amd64: the x86_64 binaries "
        "with `bgzip`, `tabix` and `ps`)"
    )
    assert f"from tag `{TAG}` (commit `{COMMIT}`) with Rust {TOOLCHAIN}; every file in [`SHA256SUMS`]({DOWNLOAD}/SHA256SUMS)" in install[5]
    assert install[6].endswith(f"sha256 `{DIGESTS['vep-rs-0.3.2.tar.gz']}`") and "DOI" not in install[6]
    reordered = fixture_repo(tmp_path / "reordered")
    (reordered / "SHA256SUMS").write_text("".join(f"{DIGESTS[n]}  {n}\n" for n in reversed(list(DIGESTS))))
    assert page(reordered) == body


def test_install_refuses_an_incomplete_or_unknown_sha256sums(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    sums = repo / "SHA256SUMS"
    sums.write_text("".join(f"{d}  {n}\n" for n, d in DIGESTS.items() if "apple" not in n))
    assert "SHA256SUMS lacks the archive for ['aarch64-apple-darwin']" in error(repo)
    sums.write_text("".join(f"{d}  {n}\n" for n, d in DIGESTS.items() if "vep-rs-" not in n))
    assert "SHA256SUMS lacks vep-rs-0.3.2.tar.gz" in error(repo)
    sums.write_text("".join(f"{d}  {n}\n" for n, d in DIGESTS.items()) + "e" * 64 + "  vep-0.3.2-riscv64gc-unknown-linux-gnu.tar.gz\n")
    assert "SHA256SUMS lists ['vep-0.3.2-riscv64gc-unknown-linux-gnu.tar.gz'], archives of no known target" in error(repo)
    sums.unlink()
    assert "no SHA256SUMS at" in error(repo)


GENERIC_DIGESTS = {
    "vep-0.3.2-x86_64-unknown-linux-gnu-v2.tar.gz": "1" * 64,
    "vep-0.3.2-aarch64-unknown-linux-gnu-generic.tar.gz": "2" * 64,
}
GENERIC_PLATFORMS = {
    "vep-0.3.2-x86_64-unknown-linux-gnu-v2.tar.gz": (
        "Linux x86_64, glibc 2.17 or later, x86-64-v2 (SSE4.2, POPCNT: Intel Nehalem 2008 or later, AMD Bulldozer 2011 "
        "or later, and VMs that hide AVX2)"
    ),
    "vep-0.3.2-aarch64-unknown-linux-gnu-generic.tar.gz": (
        "Linux aarch64, glibc 2.17 or later, Armv8.0-A with NEON (any 64-bit Arm: AWS Graviton2, Ampere Altra, "
        "Raspberry Pi, Apple silicon under Linux)"
    ),
}


def with_generic_archives(repo: Path, floors: dict[str, str]) -> None:
    """SHA256SUMS listing the two generic archives beside the three measured ones, with the given floors."""
    (repo / "SHA256SUMS").write_text("".join(f"{d}  {n}\n" for n, d in {**DIGESTS, **GENERIC_DIGESTS}.items()))
    write_floors(repo, floors)


def install_lines(body: str) -> list[str]:
    return body.split("## Install\n\n")[1].split("\n\n## Cite")[0].splitlines()


def test_install_states_the_glibc_floor_measured_on_each_linux_archive(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_floors(repo, {"vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz": "2.17", "vep-0.3.2-aarch64-unknown-linux-gnu.tar.gz": "2.28"})
    install = install_lines(page(repo))
    assert "Linux x86_64, glibc 2.17 or later, x86-64-v3" in install[0]
    assert "Linux aarch64, glibc 2.28 or later, Armv8.4-A with SVE" in install[1]
    assert "glibc" not in install[2] and install[2].split("`")[1] == "vep-0.3.2-aarch64-apple-darwin.tar.gz"


def test_install_refuses_floors_that_are_absent_incomplete_stale_or_not_versions(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    floors = repo / "glibc-floors.json"
    out = render(repo, "--glibc-floors", str(tmp_path / "absent.json"))
    assert out.returncode == 1 and "no glibc floors at" in out.stderr
    write_floors(repo, {"vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz": "2.17"})
    assert "has no glibc floor for ['vep-0.3.2-aarch64-unknown-linux-gnu.tar.gz']" in error(repo)
    write_floors(repo, {n: "2.17" for n in DIGESTS if "linux-gnu" in n} | {"vep-0.3.2-riscv64gc-unknown-linux-gnu.tar.gz": "2.17"})
    assert "names ['vep-0.3.2-riscv64gc-unknown-linux-gnu.tar.gz'], archives SHA256SUMS lacks" in error(repo)
    write_floors(repo, {n: "2.17" for n in DIGESTS if "linux-gnu" in n} | {"vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz": "glibc 2.17"})
    assert "carries a glibc floor for vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz that is not a version: 'glibc 2.17'" in error(repo)
    floors.write_text("[]")
    assert f"{floors} is not an object" in error(repo)
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "--repo", str(repo), "--version", VERSION, "--previous", PREVIOUS, "--commit", COMMIT,
         "--toolchain", TOOLCHAIN, "--sha256sums", str(repo / "SHA256SUMS")],
        check=False, capture_output=True, text=True,
    )
    assert out.returncode == 1 and "Linux archives need --glibc-floors" in out.stderr


def test_generic_archives_follow_their_measured_archive_with_their_own_cpu_floor(tmp_path: Path) -> None:
    """A release with the two generic builds lists five archives, x86_64 then aarch64 then macOS, each generic
    archive right after the measured archive of its architecture; a record without the generic builds' cells
    says their wall time was not measured, and the container bullet names both platforms and the fast tags."""
    repo = fixture_repo(tmp_path)
    with_generic_archives(repo, {n: "2.17" for n in {**DIGESTS, **GENERIC_DIGESTS} if "linux-gnu" in n})
    install = install_lines(page(repo))
    assert len(install) == 9
    names = [l.split("`")[1] for l in install[:5]]
    assert names == [
        "vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz", "vep-0.3.2-x86_64-unknown-linux-gnu-v2.tar.gz",
        "vep-0.3.2-aarch64-unknown-linux-gnu.tar.gz", "vep-0.3.2-aarch64-unknown-linux-gnu-generic.tar.gz",
        "vep-0.3.2-aarch64-apple-darwin.tar.gz",
    ]
    for line in install[:5]:
        name = line.split("`")[1]
        digest = {**DIGESTS, **GENERIC_DIGESTS}[name]
        assert line.endswith(f"; sha256 `{digest}`")
    assert install[1] == (
        f"- [`vep-0.3.2-x86_64-unknown-linux-gnu-v2.tar.gz`]({DOWNLOAD}/vep-0.3.2-x86_64-unknown-linux-gnu-v2.tar.gz): "
        f"{GENERIC_PLATFORMS['vep-0.3.2-x86_64-unknown-linux-gnu-v2.tar.gz']}, wall time against the x86-64-v3 build "
        f"not measured; sha256 `{'1' * 64}`"
    )
    assert install[3] == (
        f"- [`vep-0.3.2-aarch64-unknown-linux-gnu-generic.tar.gz`]({DOWNLOAD}/vep-0.3.2-aarch64-unknown-linux-gnu-generic.tar.gz): "
        f"{GENERIC_PLATFORMS['vep-0.3.2-aarch64-unknown-linux-gnu-generic.tar.gz']}, wall time against the neoverse-v1 build "
        f"not measured; sha256 `{'2' * 64}`"
    )
    assert install[6] == (
        "- Container image: `docker pull ghcr.io/natera-open-source/vep-rs:0.3.2` (linux/amd64 and linux/arm64: the "
        "`x86-64-v2` and the generic aarch64 binaries, which run on almost every 64-bit host; the single-platform tags "
        "`0.3.2-x86-64-v3` and `0.3.2-neoverse-v1` carry the measured binaries; every image adds `bgzip`, `tabix` and `ps`)"
    )
    assert all(l.startswith("- ") for l in install)


INTEL_MAC = "vep-0.3.2-x86_64-apple-darwin.tar.gz"
INTEL_MAC_PLATFORM = "macOS on Intel, x86-64-v3 (AVX2, BMI2, FMA: every Intel Mac macOS 26 supports)"


def test_x86_64_apple_darwin_archive_renders_its_bullet(tmp_path: Path) -> None:
    """A SHA256SUMS carrying the Intel Mac archive renders its bullet last among the archives, after the Apple
    silicon one, with its platform and CPU floor and no glibc clause, and the page changes by that line alone;
    a release without the archive renders without it, and the measured archives stay required."""
    repo = fixture_repo(tmp_path)
    without = page(repo)
    (repo / "SHA256SUMS").write_text("".join(f"{d}  {n}\n" for n, d in {**DIGESTS, INTEL_MAC: "3" * 64}.items()))
    body = page(repo)
    install = install_lines(body)
    assert len(install) == 8 and [l.split("`")[1] for l in install[:4]] == [*PLATFORMS, INTEL_MAC]
    assert install[3] == f"- [`{INTEL_MAC}`]({DOWNLOAD}/{INTEL_MAC}): {INTEL_MAC_PLATFORM}; sha256 `{'3' * 64}`"
    assert body.replace(f"{install[3]}\n", "") == without
    with_generic_archives(repo, {n: "2.17" for n in {**DIGESTS, **GENERIC_DIGESTS} if "linux-gnu" in n})
    (repo / "SHA256SUMS").write_text(f"{'3' * 64}  {INTEL_MAC}\n" + (repo / "SHA256SUMS").read_text())
    install = install_lines(page(repo))
    assert len(install) == 10 and [l.split("`")[1] for l in install[:6]] == [
        "vep-0.3.2-x86_64-unknown-linux-gnu.tar.gz", "vep-0.3.2-x86_64-unknown-linux-gnu-v2.tar.gz",
        "vep-0.3.2-aarch64-unknown-linux-gnu.tar.gz", "vep-0.3.2-aarch64-unknown-linux-gnu-generic.tar.gz",
        "vep-0.3.2-aarch64-apple-darwin.tar.gz", INTEL_MAC,
    ]
    assert install[5].split(": ", 1)[1].rsplit("; sha256", 1)[0] == INTEL_MAC_PLATFORM
    (repo / "SHA256SUMS").write_text("".join(f"{d}  {n}\n" for n, d in {**DIGESTS, INTEL_MAC: "3" * 64}.items() if "aarch64-apple" not in n))
    assert "SHA256SUMS lacks the archive for ['aarch64-apple-darwin']" in error(repo)


def generic_cells(record: dict) -> None:
    """The campaign's generic-build cells: ClinVar GRCh38 on each architecture, named by `target_cpu`, beside the
    measured build's cells; the measured build's own cells may name their CPU too."""
    s02 = {a["arch"]: a for a in record["per_cell_aggregates"] if a["suite_id"] == "s02"}
    record["per_cell_aggregates"] += [
        {**s02["arm64"], "target_cpu": "generic", "median_sec": 14.0},  # 14.0 / 12.5
        {**s02["x86_64"], "target_cpu": "x86-64-v2", "median_sec": 24.05},  # 24.05 / 21.2
    ]
    s02["arm64"]["target_cpu"] = "neoverse-v1"


def test_generic_archives_state_their_wall_time_against_the_measured_build_from_the_record(tmp_path: Path) -> None:
    repo = six_dataset_repo(tmp_path)
    with_generic_archives(repo, {n: "2.17" for n in {**DIGESTS, **GENERIC_DIGESTS} if "linux-gnu" in n})
    before = page(repo)
    edit_record(repo, "release", generic_cells)
    body = page(repo)
    install = install_lines(body)
    assert "x86-64-v2 (SSE4.2, POPCNT: Intel Nehalem 2008 or later, AMD Bulldozer 2011 or later, and VMs that hide AVX2), 1.13 times the x86-64-v3 build's wall time on ClinVar GRCh38; sha256" in install[1]
    assert "Apple silicon under Linux), 1.12 times the neoverse-v1 build's wall time on ClinVar GRCh38; sha256" in install[3]
    # The generic cells change nothing else: the measured build's ClinVar rows and caption read as before.
    assert concordance_section(body) == concordance_section(before)
    edit_record(repo, "release", lambda r: r["per_cell_aggregates"].append({**aggregate(r, "s02", "arm64"), "target_cpu": "generic"}))
    assert "the release record has 2 per_cell_aggregates entries for s02 on arm64 of the generic build; one is read" in error(repo)
    edit_record(repo, "release", lambda r: r["per_cell_aggregates"].pop())
    edit_record(repo, "release", lambda r: aggregate(r, "s02", "arm64").update(median_sec=0))
    assert "the measured build's s02 median on arm64 is zero" in error(repo)
    edit_record(repo, "release", lambda r: aggregate(r, "s02", "arm64").update(median_sec=12.5))
    edit_record(repo, "release", lambda r: next(a for a in r["per_cell_aggregates"] if a.get("target_cpu") == "generic").update(median_sec="slow"))
    assert "the per_cell_aggregates entry for s02 on arm64 of the generic build carries a median_sec that is not a number: 'slow'" in error(repo)


def test_cite_and_source_archive_carry_the_version_doi_when_minted(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    body = page(repo)
    assert f"## Cite\n\nPorter M, Borkowski R. vep-rs, version 0.3.2. 2026. doi:[{CONCEPT_DOI}](https://doi.org/{CONCEPT_DOI})\n\n**Full changelog**" in body
    assert "version DOI" not in body
    cff = repo / "CITATION.cff"
    cff.write_text(cff.read_text() + "  - type: doi\n    value: 10.5281/zenodo.99999999\n    description: version 0.3.2\n")
    body = page(repo)
    assert "doi:[10.5281/zenodo.99999999](https://doi.org/10.5281/zenodo.99999999)\n" in body
    assert f"{TAG}`; sha256 `{DIGESTS['vep-rs-0.3.2.tar.gz']}`; version DOI [10.5281/zenodo.99999999](https://doi.org/10.5281/zenodo.99999999)\n" in body
    assert CONCEPT_DOI not in body
    cff.write_text(re.sub(r"^doi: .*$", "doi: 10.5281/zenodo.88888888", cff.read_text().split("identifiers:")[0], flags=re.MULTILINE))
    assert "doi:[10.5281/zenodo.88888888](https://doi.org/10.5281/zenodo.88888888)" in page(repo)


def test_footer_is_the_changelog_and_compare_links(tmp_path: Path) -> None:
    body = page(fixture_repo(tmp_path))
    assert body.endswith(
        f"\n\n**Full changelog**: [CHANGELOG]({REPO_URL}/blob/{TAG}/CHANGELOG.md) · "
        f"[`{PREVIOUS}...{TAG}`]({REPO_URL}/compare/{PREVIOUS}...{TAG})\n"
    )
    assert "Full entry:" not in body and "**Full Changelog**: https" not in body


def test_dry_run_renders_the_unreleased_section_with_a_given_date(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_changelog(repo, "### Fixed\n\n- A released fix.")
    body = page(repo, "--changelog-section", "Unreleased", "--date", "2031-01-02")
    assert "## Changed\n\n- Unreleased work.\n\n## Concordance" in body
    assert "A released fix." not in body
    assert "vep-rs, version 0.3.2. 2031. doi:" in body


def test_the_notes_file_carries_only_the_known_sections(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_notes(repo, "A release.", "Known issues\n\n- An issue.")
    assert "unknown section '## Known issues'" in error(repo)
    write_notes(repo, "A release.", "Highlights\n\n- One.\n- Two.", "Highlights\n\n- Three.\n- Four.")
    assert error(repo) == f"{ERR}'## Highlights' appears twice in {notes(repo)}\n"


def test_missing_records_or_changelog_section_fail(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    assert "no '## [9.9.9]' section" in error(repo, "--changelog-section", "9.9.9")
    records = repo / "docs" / "concordance-provenance"
    population = next(records.glob(f"*-population-v{VERSION}.json"))
    population.unlink()
    assert f"no docs/concordance-provenance/*-population-v{VERSION}.json" in error(repo)
    release = next(records.glob(f"*-release-v{VERSION}.json"))
    record = json.loads(release.read_text())
    del record["population_concordance"]
    release.write_text(json.dumps(record))
    assert "lacks population_concordance" in error(repo)
    release.unlink()
    assert f"no docs/concordance-provenance/*-release-v{VERSION}.json" in error(repo)


def test_an_earlier_shape_record_whose_clinvar_grch37_cell_is_not_an_object_is_refused(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    edit_record(repo, "release", lambda r: r["concordance"].update(s01=7))
    assert "lacks concordance rows for ['s01']" in error(repo)


def test_a_record_that_is_not_an_object_is_refused_by_name(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    release = next((repo / "docs" / "concordance-provenance").glob(f"*-release-v{VERSION}.json"))
    release.write_text("[1, 2]")
    assert error(repo).rstrip("\n") == f"ERROR: [render_release_notes] {release.name} is not an object"
    repo = fixture_repo(tmp_path / "population")
    population = next((repo / "docs" / "concordance-provenance").glob(f"*-population-v{VERSION}.json"))
    population.write_text("7")
    assert error(repo).rstrip("\n") == f"ERROR: [render_release_notes] {population.name} is not an object"


def test_an_earlier_shape_record_whose_pooled_block_is_malformed_is_refused_by_name(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    edit_record(repo, "release", lambda r: r["population_concordance"]["pooled"].pop("excluded_perl"))
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] population_concordance.pooled carries no excluded_perl"
    edit_record(repo, "release", lambda r: r["population_concordance"].update(pooled=[1]))
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] population_concordance.pooled carries no label"


def test_a_flat_perl_reference_block_with_a_nested_object_keeps_the_earlier_shape(tmp_path: Path) -> None:
    """A flat block whose one nested value sits under a key of no suite family is not the six-dataset shape: the
    record's ClinVar GRCh37 cell carries no adjusted counts, so the earlier page renders, its caption from the block."""
    repo = fixture_repo(tmp_path)
    block = {"image": "ensemblorg/ensembl-vep:release_115.2", "cache_version": 115, "cells": {"s01": {"scored": 1}}}
    edit_record(repo, "release", lambda r: r.update(perl_reference=block))
    assert page(repo) == expected_page(repo)


IMAGE = "ensemblorg/ensembl-vep:release_116.2"
REFERENCE = {"image": IMAGE, "digest": "sha256:" + "5" * 64, "cache_version": 116}
POPULATION = (  # dataset, label, assembly, perl, vep-rs, intersection, excluded_perl, excluded_rust
    ("gnomad_v4.1_genomes", "gnomAD v4.1, whole genome", "GRCh38", 7_409_818_550, 7_409_818_550, 7_409_815_011, 3_539, 3_539),
    ("gnomad_v2.1.1_genomes", "gnomAD v2.1.1, whole genome", "GRCh37", 885_664_902, 885_664_902, 885_664_478, 424, 424),
    ("1000g_phase3", "1000 Genomes Phase 3, whole genome", "GRCh37", 286_649_478, 286_649_478, 286_649_396, 82, 82),
    ("1000g_highcov", "1000 Genomes high-coverage, whole genome", "GRCh38", 699_522_065, 699_534_118, 699_521_778, 287, 12_340),
)
CLINVAR = {  # suite: label, perl, vep-rs, intersection, adj_perl, adj_rust, adj_intersection
    "s01": ("ClinVar GRCh37", 33_483_881, 33_483_881, 33_483_180, 33_483_180, 33_483_180, 33_483_180),
    "s02": ("ClinVar GRCh38", 125_572_619, 125_572_619, 125_569_159, 125_569_159, 125_569_159, 125_569_159),
}
SYNTHETIC_FILES = (
    "01_snv", "02_mnp_complex", "03_small_indels", "04_large_explicit", "05_symbolic_del_ins", "06_symbolic_dup_inv",
    "07_cnv_repeat", "08_mobile_elements", "09_breakends", "10_special_alleles", "11_multi_allelic",
    "12_complex_imprecise", "13_vcf45_features",
)
SV_CELLS = {"s07": ("GRCh37", "gnomad_sv_v2.1_chr21", 7_000), "s08": ("GRCh38", "gnomad_sv_v4.1_chr21", 23_000)}
WALL = {  # dataset: per-chromosome medians on arm64, then on x86_64; the sums carry halves to pin ROUND_HALF_UP
    "gnomad_v4.1_genomes": ((400.125, 323.5), (650.25, 518.3)),
    "gnomad_v2.1.1_genomes": ((150.5, 131.0), (240.75, 234.75)),
    "1000g_phase3": ((44.2, 40.3), (70.5, 65.0)),
    "1000g_highcov": ((38.25, 35.25), (60.1, 57.4)),
}
CLINVAR_WALL = {("s01", "arm64"): 4.31, ("s01", "x86_64"): 6.9, ("s02", "arm64"): 12.5, ("s02", "x86_64"): 21.2}
COMPUTE = {"arm64": "c8gd.8xlarge", "x86_64": "c8id.8xlarge"}
MEMORY = {"arm64": "r8gd.8xlarge", "x86_64": "r8id.8xlarge"}
# A pair of two generations: a 9th-generation ARM class beside the 8th-generation x86 class.
NEWEST = {"arm64": "c9gd.8xlarge", "x86_64": "c8id.8xlarge"}
NEWEST_MEMORY = {"arm64": "r9gd.8xlarge", "x86_64": "r8id.8xlarge"}
NEWEST_LABEL = "the newest available class per architecture"


def pair(types: dict[str, str]) -> str:
    return f"`{types['arm64']}` (ARM) and `{types['x86_64']}` (x86)"


class FileCounts(NamedTuple):
    perl: int
    rust: int
    intersection: int
    adj_perl: int
    adj_rust: int
    adj_intersection: int


def sv_files(suite: str) -> dict[str, FileCounts]:
    """Per file of a structural-variant cell: perl, vep-rs and intersection counts, raw then adjusted. The
    integers differ per file and per cell, so no pooled figure equals a single file's."""
    _, gnomad, base = SV_CELLS[suite]
    out = {}
    for k, name in enumerate((*SYNTHETIC_FILES, "1kg_sv_chr21", "clinvar_sv_chr21", gnomad), 1):
        perl = base + 311 * k
        intersection = perl - 5 * k
        out[name] = FileCounts(perl, perl + 17 * k, intersection, intersection + k % 2, intersection + k % 3, intersection)
    return out


def six_dataset_records(compute: dict[str, str] = COMPUTE, clinvar: dict[str, str] = MEMORY, clones: int = 20) -> tuple[dict, dict]:
    """The release and population records of the six-dataset shape, every float field a coarse or wrong
    rounding of the counts beside it; the whole genome timed on `compute`, the ClinVar files on `clinvar`,
    both by `clones` machines per cell."""
    per_dataset = [
        {"dataset": ds, "label": label, "assembly": asm, "perl_tuples": p, "vep_rs_tuples": r, "intersection": i,
         "excluded_perl": xp, "excluded_rust": xr, "raw_f1": 0.5, "adjusted_f1": 0.5}
        for ds, label, asm, p, r, i, xp, xr in POPULATION
    ]
    conc = {
        s: {"perl": p, "rust": r, "intersection": i, "adj_perl": ap, "adj_rust": ar, "adj_intersection": ai, "raw_f1": 0.5, "adj_f1": 0.5}
        for s, (_, p, r, i, ap, ar, ai) in CLINVAR.items()
    }
    for s in SV_CELLS:
        files = sv_files(s)
        conc[s] = {
            "perl": sum(f.perl for f in files.values()), "rust": sum(f.rust for f in files.values()),
            "intersection": sum(f.intersection for f in files.values()), "adj_perl": sum(f.adj_perl for f in files.values()),
            "adj_rust": sum(f.adj_rust for f in files.values()), "adj_intersection": sum(f.adj_intersection for f in files.values()),
            "raw_f1": 0.5, "adj_f1": 0.5,
            "per_file": {n: {"perl_tuples": f.perl, "engine_tuples": f.rust, "intersection": f.intersection, "f1": 0.5} for n, f in files.items()},
            "per_file_adjusted": {
                n: {"perl_tuples": f.adj_perl, "engine_tuples": f.adj_rust, "intersection": f.adj_intersection, "f1": 0.5} for n, f in files.items()
            },
        }
    aggregates = [
        {"arch": arch, "suite_id": s, "median_sec": sec, "instance_type": clinvar[arch], "n_clones": clones}
        for (s, arch), sec in CLINVAR_WALL.items()
    ]
    aggregates += [{"arch": a, "suite_id": s, "median_sec": 1.5, "instance_type": clinvar[a], "n_clones": clones} for s in SV_CELLS for a in clinvar]
    release = {
        "release": TAG, "record": "release-addendum", "concordance": conc, "per_cell_aggregates": aggregates,
        "perl_reference": {
            "snv_indel": {**REFERENCE, "prefix": "canonical_116_fixture", "fork": 16},
            "sv": {**REFERENCE, "prefix": "sv_reference_116_fixture", "fork": 4},
        },
        "population_concordance": {
            "release": TAG, "per_dataset": per_dataset, "pooled": {"label": "All four", "raw_f1": 0.5, "adjusted_f1": 0.5},
            "perl_reference": {**REFERENCE, "prefix": "population_r116.2_c116_fixture"},
        },
    }
    per_cell = [
        {"dataset": ds, "label": label, "assembly": asm, "chromosome": str(c), "arch": arch, "instance_type": compute[arch],
         "n_clones": clones, "median_sec": sec, "peak_rss_mb_median": 2_501.0 if (ds, c, arch) == ("gnomad_v4.1_genomes", 2, "x86_64") else 1_900.0}
        for ds, label, asm, *_ in POPULATION
        for arch, medians in zip(("arm64", "x86_64"), WALL[ds])
        for c, sec in enumerate(medians, 1)
    ]
    population = {
        "release": TAG, "record": "population-measurements", "fork": 16,
        "wall_time": {
            "vep_rs": {"engine": "vep-rs", "engine_version": TAG, "instance_types": compute, "n_clones_per_cell": clones, "per_cell": per_cell},
            "ensembl_vep": {"engine": "ensembl_vep", "engine_version": "release_115.2", "instance_types": MEMORY},
        },
        "perl_reference": {**REFERENCE, "prefix": "population_r116.2_c116_fixture"},
    }
    return release, population


def six_dataset_repo(tmp_path: Path, **shape) -> Path:
    """The 0.3.2 scratch repository with its two records replaced by records of the six-dataset shape under
    the same version: the renderer reads the shape from the record, never from the version."""
    repo = fixture_repo(tmp_path)
    records = repo / "docs" / "concordance-provenance"
    for old in records.glob(f"*-v{VERSION}.json"):
        old.unlink()
    release, population = six_dataset_records(**shape)
    (records / f"six-dataset-release-v{VERSION}.json").write_text(json.dumps(release))
    (records / f"six-dataset-population-v{VERSION}.json").write_text(json.dumps(population))
    return repo


def pooled(cells: list[tuple[int, int, int]]) -> str:
    return f1(sum(c[0] for c in cells), sum(c[1] for c in cells), sum(c[2] for c in cells))


def expected_caption(clones: int, population_classes: str, clinvar_classes: str, single_file: str = "") -> str:
    """The six-dataset caption over the fixture records: `population_classes` is the pair the whole genome was
    timed on, `clinvar_classes` the pair the ClinVar files were timed on with what differs from the whole-genome
    measurement, `single_file` the clause a `single_file_whole_genome` block adds."""
    return (
        "Whole genome against Ensembl VEP 116.2: the population datasets one run per chromosome, each ClinVar release one "
        f"genome-wide file; wall time is the sum of the per-chromosome medians of {clones} machines at 16 threads on "
        f"{population_classes} for the population datasets, and the whole-file median of {clones} machines on "
        f"{clinvar_classes}, so the `All six` row carries no wall time; peak memory is each dataset's largest per-chromosome "
        "median, at most 2.6 GB on any chromosome, and a whole-file cell's own median where its campaign recorded one"
        f"{single_file}; tuple counts, per-chromosome rows and per-file counts are in the "
        f"[provenance records]({REPO_URL}/tree/{TAG}/docs/concordance-provenance)."
    )


def expected_six_dataset_section(caption: str | None = None, clinvar_peaks: dict[tuple[str, str], str] | None = None) -> list[str]:
    """The six-dataset section over the fixture records; `clinvar_peaks` is the peak RSS cell per (ClinVar cell,
    arch) of the ClinVar rows, empty where the fixture aggregate carries none."""
    clinvar_peaks = clinvar_peaks or {}
    raw: list[tuple[int, int, int]] = []
    adjusted: list[tuple[int, int, int]] = []
    rows = [
        "| Dataset | Raw F1 | Adjusted F1 | ARM (s) | x86 (s) | ARM peak (GB) | x86 peak (GB) |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for ds, label, _, p, r, i, xp, xr in POPULATION:
        raw.append((i, p, r))
        adjusted.append((i, p - xp, r - xr))
        arm, x86 = (sum(Fraction(str(m)) for m in medians) for medians in WALL[ds])
        # Every fixture cell peaks at 1,900 MB but gnomAD v4.1 chromosome 2 on x86, at 2,501 MB: 2.6 GB rounded up.
        x86_peak = "2.6" if ds == "gnomad_v4.1_genomes" else "1.9"
        rows.append(
            f"| {label.replace(', whole genome', '')} | {f1(i, p, r)} | {f1(i, p - xp, r - xr)} | {half_up(arm, 0)} | "
            f"{half_up(x86, 0)} | 1.9 | {x86_peak} |"
        )
    for s, (label, p, r, i, ap, ar, ai) in CLINVAR.items():
        raw.append((i, p, r))
        adjusted.append((ai, ap, ar))
        arm, x86 = (half_up(Fraction(str(CLINVAR_WALL[(s, arch)])), 0) for arch in ("arm64", "x86_64"))
        peaks = "".join(f" {clinvar_peaks[(s, a)]} |" if (s, a) in clinvar_peaks else " |" for a in ("arm64", "x86_64"))
        rows.append(f"| {label} | {f1(i, p, r)} | {f1(ai, ap, ar)} | {arm} | {x86} |{peaks}")
    rows.append(f"| All six | {pooled(raw)} | {pooled(adjusted)} | | | | |")
    sv = ["| File set | Assembly | Raw F1 | Adjusted F1 |", "| --- | --- | --- | --- |"]
    for s, (assembly, gnomad, _) in SV_CELLS.items():
        files = sv_files(s)
        for label, names in (("ClinVar SV", ["clinvar_sv_chr21"]), ("gnomAD SV", [gnomad]), ("1000 Genomes SV", ["1kg_sv_chr21"]),
                             ("Synthetic (13 files)", list(SYNTHETIC_FILES))):
            raw_files = [(files[n].intersection, files[n].perl, files[n].rust) for n in names]
            adjusted_files = [(files[n].adj_intersection, files[n].adj_perl, files[n].adj_rust) for n in names]
            sv.append(f"| {label} | {assembly} | {pooled(raw_files)} | {pooled(adjusted_files)} |")
    if caption is None:
        caption = expected_caption(20, pair(COMPUTE), f"{pair(MEMORY)} for the ClinVar files, another instance class and protocol")
    sv_caption = (
        "Structural-variant file sets measured on the released binary against the same reference, each set's files pooled by "
        "summed counts; the paper's own figures stay as published in "
        f"[docs/published-figures.md]({REPO_URL}/blob/{TAG}/docs/published-figures.md)."
    )
    return [caption, "\n".join(rows), sv_caption, "\n".join(sv)]


def concordance_section(body: str) -> str:
    return body.split("## Concordance and wall time\n\n", 1)[1].split("\n\n## Install", 1)[0]


def test_six_dataset_records_render_six_rows_all_six_and_the_eight_structural_variant_rows(tmp_path: Path) -> None:
    repo = six_dataset_repo(tmp_path)
    body = page(repo)
    section = concordance_section(body)
    caption, table, sv_caption, sv_table = expected_six_dataset_section()
    assert section == f"{caption}\n\n{table}\n\n{sv_caption}\n\n{sv_table}"
    rows = table.splitlines()[2:]
    assert [r.split(" | ")[0][2:] for r in rows] == [
        "gnomAD v4.1", "gnomAD v2.1.1", "1000 Genomes Phase 3", "1000 Genomes high-coverage", "ClinVar GRCh37", "ClinVar GRCh38", "All six",
    ]
    assert rows[0].startswith("| gnomAD v4.1 | 0.999999522 | 1.000000000 | 724 | 1,169 |")
    assert rows[4] == "| ClinVar GRCh37 | 0.999979065 | 1.000000000 | 4 | 7 | | |"
    assert rows[5] == "| ClinVar GRCh38 | 0.999972446 | 1.000000000 | 13 | 21 | | |", "12.5 s rounds half up, never to even"
    assert rows[6] == "| All six | 0.999998462 | 1.000000000 | | | | |"
    sv_rows = sv_table.splitlines()[2:]
    assert [r.split(" | ")[0][2:] for r in sv_rows] == ["ClinVar SV", "gnomAD SV", "1000 Genomes SV", "Synthetic (13 files)"] * 2
    assert [r.split(" | ")[1] for r in sv_rows] == ["GRCh37"] * 4 + ["GRCh38"] * 4
    assert sv_table.count("1.000000000") == 0 and len(set(sv_rows)) == 8
    assert "chromosome 21" not in section and "All four" not in body and "ClinVar full" not in body and "SV per-VCF" not in body
    assert "115.2" not in section and "r8gd" in caption and "c8gd" in caption
    assert len([l for l in section.splitlines() if l.strip()]) <= 23
    assert len(body.splitlines()) <= 60
    assert [l for l in body.splitlines() if l.startswith("#")] == [
        "## Upgrade notes", "## Changed", "## Concordance and wall time", "## Install", "## Cite",
    ]

    # A record carrying the chromosome 21 cells renders no row for them; the record's own rounded F1 fields are not read.
    def add_development_cells_and_spoil_floats(record: dict) -> None:
        for s in ("s03", "s04", "s05", "s06"):
            record["concordance"][s] = {"perl": 10, "rust": 10, "intersection": 9, "raw_f1": 0.9, "adj_f1": 1.0}
        for cell in record["concordance"].values():
            cell["raw_f1"], cell["adj_f1"] = 0.1, 0.2
            for block in ("per_file", "per_file_adjusted"):
                for entry in cell.get(block, {}).values():
                    entry["f1"] = 0.3
        for r in [*record["population_concordance"]["per_dataset"], record["population_concordance"]["pooled"]]:
            r["raw_f1"], r["adjusted_f1"] = 0.1, 0.2

    edit_record(repo, "release", add_development_cells_and_spoil_floats)
    assert page(repo) == body


def test_a_pair_of_two_generations_is_labelled_the_newest_available_class_per_architecture(tmp_path: Path) -> None:
    """Records timed on a 9th-generation ARM class beside the 8th-generation x86 class, five machines per cell,
    the ClinVar files on the whole genome's classes: the caption labels the pair, carries the record's clone
    count and says the ClinVar files differ in protocol alone; a same-generation pair carries no label."""
    repo = six_dataset_repo(tmp_path, compute=NEWEST, clinvar=NEWEST, clones=5)
    body = page(repo)
    section = concordance_section(body)
    caption, table, sv_caption, sv_table = expected_six_dataset_section(
        expected_caption(5, f"{pair(NEWEST)}, {NEWEST_LABEL},", "the same classes for the ClinVar files, another protocol")
    )
    assert section == f"{caption}\n\n{table}\n\n{sv_caption}\n\n{sv_table}"
    assert caption.count(NEWEST_LABEL) == 1 and caption.count("5 machines") == 2 and caption.count("c9gd.8xlarge") == 1
    assert "20 machines" not in section and "another instance class" not in section
    assert "c8gd" not in section and "r8gd" not in section and "r8id" not in section
    assert len(body.splitlines()) <= 60

    # The label follows each pair's own generations: the whole genome's pair alone, the ClinVar pair alone, neither.
    repo = six_dataset_repo(tmp_path / "clinvar-on-the-memory-class", compute=NEWEST, clinvar=MEMORY, clones=5)
    assert expected_caption(
        5, f"{pair(NEWEST)}, {NEWEST_LABEL},", f"{pair(MEMORY)} for the ClinVar files, another instance class and protocol"
    ) in page(repo)
    repo = six_dataset_repo(tmp_path / "clinvar-on-the-newer-memory-class", clinvar=NEWEST_MEMORY)
    assert expected_caption(
        20, pair(COMPUTE), f"{pair(NEWEST_MEMORY)}, {NEWEST_LABEL}, for the ClinVar files, another instance class and protocol"
    ) in page(repo)
    repo = six_dataset_repo(tmp_path / "one-generation-throughout", clinvar=COMPUTE)
    body = page(repo)
    assert expected_caption(20, pair(COMPUTE), "the same classes for the ClinVar files, another protocol") in body
    assert NEWEST_LABEL not in body


def test_six_dataset_adjusted_f1_is_computed_from_the_adjusted_counts(tmp_path: Path) -> None:
    repo = six_dataset_repo(tmp_path)
    body = page(repo)

    # 2 x 125,569,059 / (2 x 125,569,159) rounds to 0.999999204 at nine places; a rendering that copied the
    # record's six-place adj_f1 would print 0.999999000.
    def lower_the_adjusted_intersection(record: dict) -> None:
        record["concordance"]["s02"]["adj_intersection"] -= 100
        record["concordance"]["s02"]["adj_f1"] = 0.999999

    edit_record(repo, "release", lower_the_adjusted_intersection)
    moved = page(repo)
    assert "| ClinVar GRCh38 | 0.999972446 | 0.999999204 | 13 | 21 |" in moved
    assert f1(125_569_059, 125_569_159, 125_569_159) == "0.999999204" and "0.999999000" not in moved
    all_six = next(l for l in moved.splitlines() if l.startswith("| All six |"))
    assert all_six != next(l for l in body.splitlines() if l.startswith("| All six |")) and all_six.endswith("| | | | |")

    def drop_adjusted_counts(record: dict) -> None:
        del record["concordance"]["s01"]["adj_rust"]

    edit_record(repo, "release", drop_adjusted_counts)
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] concordance.s01 carries no adj_rust"

    def unmask_a_population_dataset(record: dict) -> None:
        record["population_concordance"]["per_dataset"][3]["excluded_rust"] = 0

    repo = six_dataset_repo(tmp_path / "unmasked")
    edit_record(repo, "release", unmask_a_population_dataset)
    _, _, _, p, r, i, xp, _ = POPULATION[3]
    assert f"| 1000 Genomes high-coverage | {f1(i, p, r)} | {f1(i, p - xp, r)} | " in page(repo) and f1(i, p - xp, r) != "1.000000000"

    def drop_a_population_dataset(record: dict) -> None:
        del record["population_concordance"]["per_dataset"][3]

    repo = six_dataset_repo(tmp_path / "three")
    edit_record(repo, "release", drop_a_population_dataset)
    assert error(repo).rstrip("\n") == (
        "ERROR: [render_release_notes] the population concordance has 3 datasets; the six-dataset table pools four with the two ClinVar files"
    )


def test_records_naming_different_reference_releases_are_refused_with_one_error_line(tmp_path: Path) -> None:
    repo = six_dataset_repo(tmp_path)
    earlier = {"image": "ensemblorg/ensembl-vep:release_115.2", "cache_version": 115, "prefix": "population_r115.2_c115_fixture"}

    def population_scored_against_the_earlier_release(record: dict) -> None:
        record["perl_reference"] = earlier

    edit_record(repo, "population", population_scored_against_the_earlier_release)
    err = error(repo)
    assert err.count("\n") == 1 and err == (
        "ERROR: [render_release_notes] the records name different Ensembl VEP releases: the release record.snv_indel 116.2, "
        "the release record.sv 116.2, population_concordance 116.2, the population record 115.2\n"
    )

    repo = six_dataset_repo(tmp_path / "families")

    def families_disagree(record: dict) -> None:
        record["perl_reference"]["sv"] = earlier

    edit_record(repo, "release", families_disagree)
    assert "the release record.sv 115.2" in error(repo)

    def no_block_on_the_release_record(record: dict) -> None:
        del record["perl_reference"]

    edit_record(repo, "release", no_block_on_the_release_record)
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] the release record carries no perl_reference block"

    repo = six_dataset_repo(tmp_path / "population")

    def no_block_on_the_population_record(record: dict) -> None:
        del record["perl_reference"]

    edit_record(repo, "population", no_block_on_the_population_record)
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] the population record carries no perl_reference block"
    assert "115.2" not in error(repo), "the engine version of a carried-forward wall-time run is never the reference"


def test_structural_variant_sets_follow_the_file_names_and_a_malformed_cell_is_refused(tmp_path: Path) -> None:
    repo = six_dataset_repo(tmp_path)

    def rename_a_file(record: dict) -> None:
        cell = record["concordance"]["s08"]
        for block in ("per_file", "per_file_adjusted"):
            cell[block]["extra_set_chr21"] = cell[block].pop("gnomad_sv_v4.1_chr21")

    edit_record(repo, "release", rename_a_file)
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] concordance.s08.per_file lists extra_set_chr21, a file of no known set"

    repo = six_dataset_repo(tmp_path / "twelve")

    def drop_a_synthetic_file(record: dict) -> None:
        for block in ("per_file", "per_file_adjusted"):
            del record["concordance"]["s07"][block]["13_vcf45_features"]

    edit_record(repo, "release", drop_a_synthetic_file)
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] concordance.s07.per_file lists 12 synthetic files; the set has 13"

    repo = six_dataset_repo(tmp_path / "adjusted")

    def drop_a_file_from_the_adjusted_block(record: dict) -> None:
        del record["concordance"]["s07"]["per_file_adjusted"]["clinvar_sv_chr21"]

    edit_record(repo, "release", drop_a_file_from_the_adjusted_block)
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] concordance.s07 per_file and per_file_adjusted list different files"

    repo = six_dataset_repo(tmp_path / "counts")

    def drop_a_count(record: dict) -> None:
        del record["concordance"]["s08"]["per_file"]["09_breakends"]["engine_tuples"]

    edit_record(repo, "release", drop_a_count)
    assert error(repo).rstrip("\n") == "ERROR: [render_release_notes] concordance.s08.per_file 09_breakends carries no engine_tuples"

    repo = six_dataset_repo(tmp_path / "aggregate")

    def drop_a_wall_time(record: dict) -> None:
        record["per_cell_aggregates"] = [a for a in record["per_cell_aggregates"] if (a["suite_id"], a["arch"]) != ("s02", "x86_64")]

    edit_record(repo, "release", drop_a_wall_time)
    assert error(repo).rstrip("\n") == (
        "ERROR: [render_release_notes] the release record has 0 per_cell_aggregates entries for s02 on x86_64; one is read"
    )

    repo = six_dataset_repo(tmp_path / "cells")

    def drop_the_grch38_sv_cell(record: dict) -> None:
        del record["concordance"]["s08"]

    edit_record(repo, "release", drop_the_grch38_sv_cell)
    assert "lacks concordance rows for ['s08']" in error(repo)


ERR = "ERROR: [render_release_notes] "
NO_RELEASE = (
    "names no release: provenance carries neither an image tag (perl_docker_image, perl_image, perl_vep_docker_image, "
    "vep_image_tag, image, perl_reference.image) nor a cache_version"
)


def aggregate(record: dict, suite: str, arch: str) -> dict:
    return next(a for a in record["per_cell_aggregates"] if (a["suite_id"], a["arch"]) == (suite, arch))


def sv_blocks(record: dict, suite: str) -> tuple[dict, dict]:
    cell = record["concordance"][suite]
    return cell["per_file"], cell["per_file_adjusted"]


def copy_file(record: dict, suite: str, source: str, name: str) -> None:
    for block in sv_blocks(record, suite):
        block[name] = dict(block[source])


def drop_file(record: dict, suite: str, name: str) -> None:
    for block in sv_blocks(record, suite):
        del block[name]


def zero_file(record: dict, suite: str, name: str) -> None:
    for block in sv_blocks(record, suite):
        block[name] = {"perl_tuples": 0, "engine_tuples": 0, "intersection": 0, "f1": 1.0}


def drop_keys(block: dict, *keys: str) -> None:
    for key in keys:
        del block[key]


def keep_timed_cells(record: dict, keep) -> None:
    record["wall_time"]["vep_rs"]["per_cell"] = [c for c in record["wall_time"]["vep_rs"]["per_cell"] if keep(c)]


def timed_cell(record: dict, index: int) -> dict:
    return record["wall_time"]["vep_rs"]["per_cell"][index]


def set_every_median(record: dict, dataset: str, arch: str, value) -> None:
    for c in record["wall_time"]["vep_rs"]["per_cell"]:
        if (c["dataset"], c["arch"]) == (dataset, arch):
            c["median_sec"] = value


def every_aggregate_engine(record: dict, engine: str) -> None:
    for a in record["per_cell_aggregates"]:
        a["engine"] = engine


def add_an_ensembl_vep_twin(record: dict) -> None:
    every_aggregate_engine(record, "vep_rs")
    record["per_cell_aggregates"].append({**aggregate(record, "s01", "arm64"), "engine": "ensembl_vep", "median_sec": 999})


def released(record: dict) -> dict:
    return record["wall_time"]["vep_rs"]


def retime_clinvar(record: dict, arch: str, instance_type: str) -> None:
    for a in record["per_cell_aggregates"]:
        if a["arch"] == arch:
            a["instance_type"] = instance_type


RELEASED = "the population record's wall_time.vep_rs"
SIX_DATASET_PROBES = [  # record kind, its edit, then the exact ERROR line or a line of the page that renders
    pytest.param(
        "population", lambda r: drop_keys(r, "wall_time"),
        ERR + "the population record carries no wall_time", id="no-wall-time",
    ),
    pytest.param(
        "population", lambda r: drop_keys(r["wall_time"], "vep_rs"),
        ERR + "the population record's wall_time carries no vep_rs", id="no-released-binary-under-wall-time",
    ),
    pytest.param(
        "population", lambda r: r["wall_time"].update(vep_rs=[1]),
        ERR + f"{RELEASED} is not an object", id="a-wall-time-block-that-is-a-list",
    ),
    pytest.param(
        "population", lambda r: drop_keys(released(r), "instance_types"),
        ERR + f"{RELEASED} carries no instance_types", id="no-instance-types",
    ),
    pytest.param(
        "population", lambda r: drop_keys(released(r), "n_clones_per_cell"),
        ERR + f"{RELEASED} carries no n_clones_per_cell", id="no-clone-count",
    ),
    pytest.param(
        "population", lambda r: drop_keys(released(r), "instance_types", "n_clones_per_cell"),
        ERR + f"{RELEASED} carries no instance_types, n_clones_per_cell", id="neither-instance-types-nor-clone-count",
    ),
    pytest.param(
        "population", lambda r: released(r).update(n_clones_per_cell=5.0),
        ERR + f"{RELEASED} carries n_clones_per_cell that is not an integer: 5.0", id="a-clone-count-that-is-not-an-integer",
    ),
    pytest.param(
        "population", lambda r: drop_keys(r, "fork"),
        ERR + "the population record carries no fork", id="no-fork",
    ),
    pytest.param(
        "population", lambda r: r.update(fork="16"),
        ERR + "the population record carries fork that is not an integer: '16'", id="a-fork-that-is-not-an-integer",
    ),
    pytest.param(
        "population", lambda r: drop_keys(released(r)["instance_types"], "x86_64"),
        ERR + f"{RELEASED}.instance_types carries no x86_64", id="no-x86-instance-type",
    ),
    pytest.param(
        "population", lambda r: released(r).update(instance_types="c8gd.8xlarge"),
        ERR + f"{RELEASED}.instance_types carries no arm64, x86_64", id="instance-types-that-are-one-string",
    ),
    pytest.param(
        "population", lambda r: released(r)["instance_types"].update(arm64=7),
        ERR + f"{RELEASED}.instance_types carries arm64 that is not a string: 7", id="an-instance-type-that-is-not-a-string",
    ),
    pytest.param(
        "population", lambda r: released(r)["instance_types"].update(arm64="metal"),
        ERR + f"{RELEASED} names the instance type 'metal', whose name carries no generation number",
        id="an-instance-type-without-a-generation",
    ),
    pytest.param(
        "release", lambda r: retime_clinvar(r, "x86_64", "metal"),
        ERR + "the release record names the instance type 'metal', whose name carries no generation number",
        id="a-clinvar-instance-type-without-a-generation",
    ),
    pytest.param(
        "population", lambda r: drop_keys(released(r), "per_cell"),
        ERR + f"{RELEASED} carries no per_cell", id="no-per-cell",
    ),
    pytest.param(
        "population", lambda r: released(r).update(per_cell={}),
        ERR + f"{RELEASED}.per_cell is not a list", id="a-per-cell-that-is-an-object",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s01", "arm64").update(n_clones=19),
        (
            "the whole-file median of 19 (ClinVar GRCh37, ARM), 20 (ClinVar GRCh37, x86), 20 (ClinVar GRCh38, ARM) and "
            "20 (ClinVar GRCh38, x86) machines on `r8gd.8xlarge` (ARM) and `r8id.8xlarge` (x86) for the ClinVar files"
        ),
        id="clone-counts-differ-and-the-caption-names-each-cell",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s01", "arm64").update(instance_type="r7gd.8xlarge"),
        ERR + "the ClinVar cells name 2 ARM instance types: ['r7gd.8xlarge', 'r8gd.8xlarge']",
        id="arm-instance-types-differ",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s02", "x86_64").update(instance_type="r7id.8xlarge"),
        ERR + "the ClinVar cells name 2 x86 instance types: ['r7id.8xlarge', 'r8id.8xlarge']",
        id="x86-instance-types-differ",
    ),
    pytest.param(
        "release", lambda r: r["per_cell_aggregates"].append(dict(aggregate(r, "s01", "arm64"))),
        ERR + "the release record has 2 per_cell_aggregates entries for s01 on arm64; one is read",
        id="a-duplicate-aggregate",
    ),
    pytest.param(
        "release", lambda r: every_aggregate_engine(r, "ensembl_vep"),
        ERR + "the release record has 0 per_cell_aggregates entries for s01 on arm64; one is read",
        id="every-aggregate-of-the-other-engine",
    ),
    pytest.param(
        "release", add_an_ensembl_vep_twin, "| ClinVar GRCh37 | 0.999979065 | 1.000000000 | 4 | 7 |",
        id="an-ensembl-vep-twin-aggregate-is-not-read",
    ),
    pytest.param(
        "release", lambda r: r["perl_reference"].update(sv={"digest": "sha256:abc"}),
        ERR + "the perl_reference block of the release record.sv " + NO_RELEASE,
        id="a-family-part-naming-no-release",
    ),
    pytest.param(
        "release", lambda r: r.update(perl_reference={"digest": "sha256:abc"}),
        ERR + "the perl_reference block of the release record " + NO_RELEASE,
        id="a-flat-block-naming-no-release",
    ),
    pytest.param(
        "release", lambda r: copy_file(r, "s07", "13_vcf45_features", "14_extra"),
        ERR + "concordance.s07.per_file lists 14 synthetic files; the set has 13",
        id="fourteen-synthetic-files",
    ),
    pytest.param(
        "release", lambda r: copy_file(r, "s07", "clinvar_sv_chr21", "clinvar_sv_extra"),
        ERR + "concordance.s07.per_file lists 2 files of the ClinVar SV set; one is read",
        id="two-files-of-one-set",
    ),
    pytest.param(
        "release", lambda r: drop_file(r, "s08", "gnomad_sv_v4.1_chr21"),
        ERR + "concordance.s08.per_file lists 0 files of the gnomAD SV set; one is read",
        id="no-file-of-one-set",
    ),
    pytest.param(
        "release", lambda r: r["population_concordance"]["per_dataset"].append(dict(r["population_concordance"]["per_dataset"][0])),
        ERR + "the population concordance has 5 datasets; the six-dataset table pools four with the two ClinVar files",
        id="five-datasets",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["concordance"]["s07"], "per_file"),
        ERR + "concordance.s07.per_file is absent or empty", id="per-file-absent",
    ),
    pytest.param(
        "release", lambda r: r["concordance"]["s07"].update(per_file={}),
        ERR + "concordance.s07.per_file is absent or empty", id="per-file-empty",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["concordance"]["s08"], "per_file_adjusted"),
        ERR + "concordance.s08.per_file_adjusted is absent or empty", id="per-file-adjusted-absent",
    ),
    pytest.param(
        "release", lambda r: r["concordance"]["s08"]["per_file"].update({"09_breakends": 7}),
        ERR + "concordance.s08.per_file 09_breakends carries no intersection, perl_tuples, engine_tuples",
        id="a-per-file-entry-that-is-not-an-object",
    ),
    pytest.param(
        "population", lambda r: keep_timed_cells(r, lambda c: c["dataset"] != "1000g_phase3"),
        ERR + "the population record has no wall time for 1000g_phase3", id="no-wall-time-for-a-dataset",
    ),
    pytest.param(
        "population", lambda r: keep_timed_cells(r, lambda c: (c["dataset"], c["arch"]) != ("1000g_highcov", "x86_64")),
        ERR + "the population record has no wall time for 1000g_highcov", id="no-wall-time-on-one-architecture",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["concordance"]["s02"], "intersection", "perl", "rust"),
        ERR + "concordance.s02 carries no intersection, perl, rust", id="a-clinvar-cell-without-raw-counts",
    ),
    pytest.param(
        "release", lambda r: drop_keys(aggregate(r, "s01", "arm64"), "median_sec"),
        ERR + "the per_cell_aggregates entry for s01 on arm64 carries no median_sec", id="an-aggregate-without-median-sec",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s01", "arm64").update(median_sec=None),
        ERR + "the per_cell_aggregates entry for s01 on arm64 carries no median_sec", id="an-aggregate-with-a-null-median-sec",
    ),
    pytest.param(
        "release", lambda r: drop_keys(aggregate(r, "s02", "x86_64"), "instance_type", "n_clones"),
        ERR + "the per_cell_aggregates entry for s02 on x86_64 carries no instance_type, n_clones",
        id="an-aggregate-without-instance-type-or-clone-count",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s01", "arm64").update(median_sec="n/a"),
        ERR + "the per_cell_aggregates entry for s01 on arm64 carries a median_sec that is not a number: 'n/a'",
        id="an-aggregate-whose-median-sec-is-not-a-number",
    ),
    pytest.param(
        "release", lambda r: zero_file(r, "s07", "clinvar_sv_chr21"), "| ClinVar SV | GRCh37 | 1.000000000 | 1.000000000 |",
        id="a-file-with-no-tuples-on-either-side-reads-one",
    ),
    pytest.param(
        "release", lambda r: r["concordance"]["s01"].update(perl=0, rust=0, intersection=0),
        "| ClinVar GRCh37 | 1.000000000 | 1.000000000 | 4 | 7 |", id="a-cell-with-no-raw-tuples-on-either-side-reads-one",
    ),
    pytest.param(
        "release", lambda r: r["concordance"]["s01"].update(adj_perl=0, adj_rust=0, adj_intersection=0),
        "| ClinVar GRCh37 | 0.999979065 | 1.000000000 | 4 | 7 |", id="a-cell-with-no-adjusted-tuples-on-either-side-reads-one",
    ),
    pytest.param(
        "release", lambda r: r["concordance"].update(s07=7),
        ERR + f"six-dataset-release-v{VERSION}.json lacks concordance rows for ['s07']",
        id="a-structural-variant-cell-that-is-not-an-object",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s02", "x86_64").update(n_clones="twenty"),
        ERR + "the per_cell_aggregates entry for s02 on x86_64 carries n_clones that is not an integer: 'twenty'",
        id="an-aggregate-whose-clone-count-is-not-an-integer",
    ),
    pytest.param(
        "release", lambda r: r["concordance"]["s08"]["per_file"]["09_breakends"].update(perl_tuples=7.9),
        ERR + "concordance.s08.per_file 09_breakends carries perl_tuples that is not an integer: 7.9",
        id="a-per-file-count-that-is-not-an-integer",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s01", "arm64").update(median_sec=float("nan")),
        ERR + "the per_cell_aggregates entry for s01 on arm64 carries a median_sec that is not finite: nan",
        id="an-aggregate-whose-median-sec-is-not-finite",
    ),
    pytest.param(
        "release", lambda r: r["per_cell_aggregates"].append(7),
        ERR + "the release record's per_cell_aggregates is not a list of objects",
        id="an-aggregate-that-is-not-an-object",
    ),
    pytest.param(
        "release", lambda r: r["population_concordance"]["per_dataset"][0].update(perl_tuples=19.6),
        ERR + "population_concordance.per_dataset gnomad_v4.1_genomes carries perl_tuples that is not an integer: 19.6",
        id="a-population-count-that-is-not-an-integer",
    ),
    pytest.param(
        "release", lambda r: r["population_concordance"]["per_dataset"][1].update(intersection=True),
        ERR + "population_concordance.per_dataset gnomad_v2.1.1_genomes carries intersection that is not an integer: True",
        id="a-population-count-that-is-a-boolean",
    ),
    pytest.param(
        "release", lambda r: r.update(concordance=[1, 2]),
        ERR + "the release record's concordance is not an object", id="a-concordance-that-is-not-an-object",
    ),
    pytest.param(
        "release", lambda r: r.update(population_concordance="abc"),
        ERR + "the release record's population_concordance is not an object",
        id="a-population-concordance-that-is-not-an-object",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["population_concordance"], "per_dataset"),
        ERR + "population_concordance carries no per_dataset", id="a-population-concordance-without-per-dataset",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["population_concordance"], "pooled"),
        ERR + "population_concordance carries no pooled", id="a-population-concordance-without-pooled",
    ),
    pytest.param(
        "release", lambda r: r["population_concordance"].update(per_dataset={}),
        ERR + "the release record's population_concordance.per_dataset is not a list",
        id="a-per-dataset-that-is-not-a-list",
    ),
    pytest.param(
        "release", lambda r: r["population_concordance"].update(per_dataset=[1, 2, 3, 4]),
        ERR + "population_concordance.per_dataset[0] carries no dataset, label",
        id="a-per-dataset-entry-that-is-not-an-object",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["population_concordance"]["per_dataset"][0], "dataset"),
        ERR + "population_concordance.per_dataset[0] carries no dataset", id="a-per-dataset-entry-without-its-dataset",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["population_concordance"]["per_dataset"][1], "label"),
        ERR + "population_concordance.per_dataset[1] carries no label", id="a-per-dataset-entry-without-its-label",
    ),
    pytest.param(
        "release", lambda r: drop_keys(r["population_concordance"]["per_dataset"][0], "excluded_perl"),
        ERR + "population_concordance.per_dataset gnomad_v4.1_genomes carries no excluded_perl",
        id="a-per-dataset-entry-without-a-count",
    ),
    pytest.param(
        "release", lambda r: r.update(per_cell_aggregates={}),
        ERR + "the release record's per_cell_aggregates is not a list of objects", id="an-aggregates-block-that-is-an-object",
    ),
    pytest.param(
        "release", lambda r: aggregate(r, "s01", "arm64").update(median_sec=1e30),
        ERR + "a wall time of 1E+30 seconds is too large to round to whole seconds",
        id="an-aggregate-whose-median-sec-is-too-large-to-round",
    ),
    pytest.param(
        "population", lambda r: timed_cell(r, 0).update(median_sec="n/a"),
        ERR + "the population record's wall_time.vep_rs.per_cell[0] carries a median_sec that is not a number: 'n/a'",
        id="a-chromosome-whose-median-sec-is-not-a-number",
    ),
    pytest.param(
        "population", lambda r: timed_cell(r, 0).update(median_sec=None),
        ERR + "the population record's wall_time.vep_rs.per_cell[0] carries no median_sec",
        id="a-chromosome-with-a-null-median-sec",
    ),
    pytest.param(
        "population", lambda r: drop_keys(timed_cell(r, 3), "median_sec"),
        ERR + "the population record's wall_time.vep_rs.per_cell[3] carries no median_sec",
        id="a-chromosome-without-median-sec",
    ),
    pytest.param(
        "population", lambda r: timed_cell(r, 0).update(median_sec=float("nan")),
        ERR + "the population record's wall_time.vep_rs.per_cell[0] carries a median_sec that is not finite: nan",
        id="a-chromosome-whose-median-sec-is-not-finite",
    ),
    pytest.param(
        "population", lambda r: r["wall_time"]["vep_rs"]["per_cell"].insert(2, 7),
        ERR + "the population record's wall_time.vep_rs.per_cell[2] carries no dataset, arch, median_sec",
        id="a-chromosome-cell-that-is-not-an-object",
    ),
    pytest.param(
        "population", lambda r: set_every_median(r, "gnomad_v4.1_genomes", "arm64", 1e30),
        ERR + "a wall time of 2E+30 seconds is too large to round to whole seconds",
        id="a-dataset-whose-summed-medians-are-too-large-to-round",
    ),
]


@pytest.mark.parametrize("kind, edit, expected", SIX_DATASET_PROBES)
def test_each_six_dataset_fault_is_refused_by_its_line_or_rendered_with_its_figure(tmp_path: Path, kind: str, edit, expected: str) -> None:
    repo = six_dataset_repo(tmp_path)
    edit_record(repo, kind, edit)
    if expected.startswith(ERR):
        assert error(repo) == f"{expected}\n"
    else:
        assert expected in page(repo)


def whole_genome_rows(body: str) -> dict[str, str]:
    """The six-dataset whole-genome table's rows by dataset label."""
    table = concordance_section(body).split("\n\n")[1]
    return {r.split(" | ")[0][2:]: r for r in table.splitlines()[2:]}


def test_peak_memory_renders_per_dataset_beside_wall_time(tmp_path: Path) -> None:
    """Each population row carries the dataset's largest per-chromosome median peak RSS per architecture, in GB
    rounded up: the fixture's gnomAD v4.1 chromosome 2 on x86 peaks at 2,501 MB and every other cell at 1,900 MB,
    so that row alone reads 2.6 on x86. A ClinVar row's cells are empty until its aggregate carries a peak, and
    the caption's bound is the largest cell of all."""
    repo = six_dataset_repo(tmp_path)
    body = page(repo)
    rows = whole_genome_rows(body)
    assert rows["gnomAD v4.1"].endswith("| 724 | 1,169 | 1.9 | 2.6 |")
    for label in ("gnomAD v2.1.1", "1000 Genomes Phase 3", "1000 Genomes high-coverage"):
        assert rows[label].endswith("| 1.9 | 1.9 |"), label
    assert rows["ClinVar GRCh37"].endswith("| 4 | 7 | | |") and rows["ClinVar GRCh38"].endswith("| 13 | 21 | | |")
    assert rows["All six"].endswith("| 1.000000000 | | | | |")
    assert "peak memory is each dataset's largest per-chromosome median, at most 2.6 GB on any chromosome" in body

    # A higher peak on one chromosome moves its dataset's cell and the caption's bound alone.
    def raise_a_phase3_cell(record: dict) -> None:
        cell = next(c for c in released(record)["per_cell"] if (c["dataset"], c["arch"]) == ("1000g_phase3", "arm64"))
        cell["peak_rss_mb_median"] = 3_000.1

    edit_record(repo, "population", raise_a_phase3_cell)
    body = page(repo)
    rows = whole_genome_rows(body)
    assert rows["1000 Genomes Phase 3"].endswith("| 3.1 | 1.9 |") and rows["gnomAD v4.1"].endswith("| 1.9 | 2.6 |")
    assert "at most 3.1 GB on any chromosome" in body

    # A ClinVar aggregate that records its peak renders it; 3,100.0 MB is exactly 3.1 GB and 3,100.1 rounds up.
    def record_clinvar_peaks(record: dict) -> None:
        aggregate(record, "s01", "arm64")["peak_rss_mb_median"] = 3_100.0
        aggregate(record, "s01", "x86_64")["peak_rss_mb_median"] = 3_100.1

    edit_record(repo, "release", record_clinvar_peaks)
    rows = whole_genome_rows(page(repo))
    assert rows["ClinVar GRCh37"].endswith("| 4 | 7 | 3.1 | 3.2 |") and rows["ClinVar GRCh38"].endswith("| 13 | 21 | | |")
    table = concordance_section(page(repo)).split("\n\n")[1]
    expected = expected_six_dataset_section(clinvar_peaks={("s01", "arm64"): "3.1", ("s01", "x86_64"): "3.2"})[1]
    expected = "\n".join(
        r[: -len("| 1.9 | 1.9 |")] + "| 3.1 | 1.9 |" if r.startswith("| 1000 Genomes Phase 3 |") else r
        for r in expected.splitlines()
    )
    assert table == expected

    edit_record(repo, "release", lambda r: aggregate(r, "s01", "arm64").update(peak_rss_mb_median="n/a"))
    assert error(repo) == (
        f"{ERR}the per_cell_aggregates entry for s01 on arm64 carries a peak_rss_mb_median that is not a number: 'n/a'\n"
    )


SINGLE_FILE = {
    "label": "gnomAD v4.1, whole genome",
    "assembly": "GRCh38",
    "fork": 16,
    "per_arch": [
        {"arch": "arm64", "instance_type": "c9gd.8xlarge", "n_clones": 5, "median_sec": 640.4, "peak_rss_mb_median": 15_321.0},
        {"arch": "x86_64", "instance_type": "c8id.8xlarge", "n_clones": 5, "median_sec": 1_049.6, "peak_rss_mb_median": 16_040.2},
    ],
}
SINGLE_FILE_CLAUSE = (
    "; the whole of gnomAD v4.1 GRCh38 as one input file ran in 640 s at a peak of 15.4 GB on ARM and 1,050 s at a "
    "peak of 16.1 GB on x86, the median of 5 machines at 16 threads"
)


def test_the_single_file_whole_genome_bound_renders_from_the_record_when_present(tmp_path: Path) -> None:
    """A release record without a `single_file_whole_genome` block renders no such clause; one with the block
    states the dataset as one input file with its wall time and peak RSS per architecture and the clone count,
    and a malformed block is refused by name."""
    repo = six_dataset_repo(tmp_path)
    assert "as one input file" not in page(repo)
    edit_record(repo, "release", lambda r: r.update(single_file_whole_genome=copy.deepcopy(SINGLE_FILE)))
    body = page(repo)
    caption = concordance_section(body).split("\n\n")[0]
    assert caption == expected_caption(
        20, pair(COMPUTE), f"{pair(MEMORY)} for the ClinVar files, another instance class and protocol", SINGLE_FILE_CLAUSE
    )
    assert len(body.splitlines()) <= 60

    def with_block(edit) -> None:
        block = copy.deepcopy(SINGLE_FILE)
        edit(block)
        edit_record(repo, "release", lambda r: r.update(single_file_whole_genome=block))

    with_block(lambda b: b["per_arch"][1].update(n_clones=3))
    assert "the median of 5 and 3 machines at 16 threads" in page(repo)
    where = f"{ERR}the release record's single_file_whole_genome"
    for edit, expected in (
        (lambda b: drop_keys(b, "assembly"), f"{where} carries no assembly"),
        (lambda b: b.update(fork="16"), f"{where} carries fork that is not an integer: '16'"),
        (lambda b: b.update(per_arch={}), f"{where}.per_arch is not a list"),
        (lambda b: drop_keys(b["per_arch"][0], "peak_rss_mb_median"), f"{where}.per_arch[0] carries no peak_rss_mb_median"),
        (lambda b: b["per_arch"][1].update(median_sec="slow"), f"{where}.per_arch[1] carries a median_sec that is not a number: 'slow'"),
        (lambda b: b["per_arch"][1].update(arch="riscv64"), f"{where}.per_arch[1] names the architecture 'riscv64'; one of arm64, x86_64 is read"),
        (lambda b: b["per_arch"][1].update(arch="arm64"), f"{where}.per_arch carries arm64 twice"),
        (lambda b: b["per_arch"].pop(), f"{where}.per_arch carries no x86_64"),
    ):
        with_block(edit)
        assert error(repo) == f"{expected}\n"


SAME_DAY_RATIOS = {  # dataset: the stored ratio on arm64, then on x86_64, at the record's three decimals
    "gnomad_v4.1_genomes": (0.926, 0.934),
    "gnomad_v2.1.1_genomes": (0.945, 1.004),
    "1000g_phase3": (0.871, 1.207),
    "1000g_highcov": (0.9, 0.995),
}
# Each stored ratio at two decimals, rounded by hand: 0.926 rounds up to 0.93 and 0.934 down to 0.93; 0.945 and 0.995
# are halves and round up to 0.95 and 1.00; 1.004 rounds down to 1.00, 0.871 down to 0.87 and 1.207 up to 1.21; 0.9
# gains its second decimal.
SAME_DAY_FIGURES = (
    "0.93 (ARM) and 0.93 (x86) for gnomAD v4.1, 0.95 (ARM) and 1.00 (x86) for gnomAD v2.1.1, 0.87 (ARM) and 1.21 (x86) "
    "for 1000 Genomes Phase 3, 0.90 (ARM) and 1.00 (x86) for 1000 Genomes high-coverage"
)
SAME_DAY = "the population record's previous_release_wall_time"


def same_day_sentence(figures: str, compared: str = PREVIOUS[1:]) -> str:
    return (
        f"The summed wall time of vep-rs {VERSION} over that of vep-rs {compared}, whose archived binaries were timed on "
        f"the same day, instance classes and clone count, is {figures}."
    )


def same_day_block(datasets: list[tuple[str, str, str, float, float]]) -> dict:
    """A `previous_release_wall_time` block over `datasets` (id, label, assembly, then the stored ratio on arm64
    and on x86_64), its entries in that order with arm64 before x86_64. Every entry's two sums divide to 0.5,
    which no stored ratio rounds to."""
    return {
        "release": PREVIOUS, "campaign_id": "fixture-retiming", "release_campaign_id": "fixture-campaign", "n_clones_per_cell": 20,
        "instance_types": dict(COMPUTE), "binary_sha256": {"arm64": "a" * 64, "x86_64": "b" * 64},
        "warmup_run_discarded_per_clone": True, "note": "A fixture.",
        "per_dataset": [
            {"dataset": ds, "label": label, "assembly": asm, "arch": arch, "n_cells": 2, "release_sum_of_cell_medians_sec": 100.0,
             "previous_release_sum_of_cell_medians_sec": 200.0, "ratio_release_over_previous": ratio}
            for ds, label, asm, *ratios in datasets
            for arch, ratio in zip(("arm64", "x86_64"), ratios)
        ],
    }


def six_dataset_same_day_block() -> dict:
    return same_day_block([(ds, label, asm, *SAME_DAY_RATIOS[ds]) for ds, label, asm, *_ in POPULATION])


def same_day(record: dict) -> dict:
    return record["previous_release_wall_time"]


def test_the_same_day_comparison_renders_as_one_caption_sentence_from_the_stored_ratios(tmp_path: Path) -> None:
    """A population record without a `previous_release_wall_time` block renders the caption as it stands; one with
    the block gains one sentence between the peak-memory clause and the records clause, which becomes a sentence
    of its own: the summed wall time of the version rendered over that of the release the block names, each
    dataset of the table under its label with the stored ratio per architecture at two decimals. Nothing else on
    the page moves, the tables included."""
    repo = six_dataset_repo(tmp_path)
    without = page(repo)
    assert without.count("; tuple counts") == 1 and "over that of" not in without
    edit_record(repo, "population", lambda r: r.update(previous_release_wall_time=six_dataset_same_day_block()))
    body = page(repo)
    sentence = same_day_sentence(SAME_DAY_FIGURES)
    assert body == without.replace("; tuple counts", f". {sentence} Tuple counts")
    caption = concordance_section(body).split("\n\n")[0]
    assert caption.endswith(
        "and a whole-file cell's own median where its campaign recorded one. "
        f"{sentence} Tuple counts, per-chromosome rows and per-file counts are in the "
        f"[provenance records]({REPO_URL}/tree/{TAG}/docs/concordance-provenance)."
    )
    assert "\n" not in caption and len(body.splitlines()) == len(without.splitlines())
    # The halves pin the rounding to half up on the stored decimal, where binary float formatting prints 0.94 and
    # 0.99; the two sums beside each ratio divide to 0.5 and are not read.
    assert (f"{0.945:.2f}", f"{0.995:.2f}") == ("0.94", "0.99") and "0.50" not in caption

    # The compared version is the block's own `release`: the footer's compare link still starts at v0.3.1.
    edit_record(repo, "population", lambda r: same_day(r).update(release="v0.5.0"))
    moved = page(repo)
    assert same_day_sentence(SAME_DAY_FIGURES, "0.5.0") in moved and f"compare/{PREVIOUS}...{TAG}" in moved
    edit_record(repo, "population", lambda r: same_day(r).update(release=PREVIOUS))

    # The sentence follows the table's dataset order and names ARM before x86, whatever the block's order.
    edit_record(repo, "population", lambda r: same_day(r)["per_dataset"].reverse())
    assert page(repo) == body

    # It stands after the single-file clause when the release record carries one.
    edit_record(repo, "release", lambda r: r.update(single_file_whole_genome=copy.deepcopy(SINGLE_FILE)))
    assert f"recorded one{SINGLE_FILE_CLAUSE}. {sentence} Tuple counts, per-chromosome rows" in page(repo)

    # A null block is an absent one.
    repo = six_dataset_repo(tmp_path / "null")
    edit_record(repo, "population", lambda r: r.update(previous_release_wall_time=None))
    assert page(repo) == without


def test_the_fixture_repo_page_carries_the_same_day_sentence_in_its_own_caption(tmp_path: Path) -> None:
    """The whole-genome caption of the page rendered from the `fixture_repo` repository takes the sentence at the same
    place, over its own table's datasets: 0.744 rounds down to 0.74 and the half 1.005 up to 1.01."""
    repo = fixture_repo(tmp_path)
    datasets = release_record(repo)["population_concordance"]["per_dataset"]
    block = same_day_block([(r["dataset"], r["label"], r["assembly"], 0.744, 1.005) for r in datasets])
    edit_record(repo, "population", lambda r: r.update(previous_release_wall_time=block))
    labels = ("gnomAD v4.1", "gnomAD v2.1.1", "1000 Genomes Phase 3", "1000 Genomes high-coverage")
    sentence = same_day_sentence(", ".join(f"0.74 (ARM) and 1.01 (x86) for {label}" for label in labels))
    expected = expected_page(repo)
    assert expected.count("; tuple counts and per-chromosome rows") == 1
    assert page(repo) == expected.replace("; tuple counts", f". {sentence} Tuple counts")


SAME_DAY_PROBES = [  # the population record's edit over the fixture block, then the exact ERROR line
    pytest.param(lambda r: drop_keys(same_day(r), "release"), f"{SAME_DAY} carries no release", id="no-release"),
    pytest.param(
        lambda r: drop_keys(same_day(r), "note", "binary_sha256"),
        f"{SAME_DAY} carries no binary_sha256, note", id="no-key-the-caption-does-not-read",
    ),
    pytest.param(
        lambda r: r.update(previous_release_wall_time=[1]),
        f"{SAME_DAY} carries no release, per_dataset, campaign_id, release_campaign_id, n_clones_per_cell, instance_types, "
        "binary_sha256, warmup_run_discarded_per_clone, note",
        id="a-block-that-is-not-an-object",
    ),
    pytest.param(
        lambda r: same_day(r).update(release="0.3.1"),
        f"{SAME_DAY} names the release '0.3.1'; a tag of the form vX.Y.Z is read", id="a-release-that-is-not-a-tag",
    ),
    pytest.param(
        lambda r: same_day(r).update(release="v0.3.1-rc1"),
        f"{SAME_DAY} names the release 'v0.3.1-rc1'; a tag of the form vX.Y.Z is read", id="a-release-tag-with-a-suffix",
    ),
    pytest.param(
        lambda r: same_day(r).update(release=31),
        f"{SAME_DAY} names the release 31; a tag of the form vX.Y.Z is read", id="a-release-that-is-not-a-string",
    ),
    pytest.param(
        lambda r: same_day(r).update(release=TAG),
        f"{SAME_DAY} names v0.3.2, the version rendered, as the release it is compared with",
        id="a-release-that-is-the-version-rendered",
    ),
    pytest.param(
        lambda r: same_day(r).update(per_dataset={}),
        f"{SAME_DAY}.per_dataset is not a list", id="a-per-dataset-that-is-not-a-list",
    ),
    pytest.param(
        lambda r: drop_keys(same_day(r)["per_dataset"][2], "ratio_release_over_previous"),
        f"{SAME_DAY}.per_dataset[2] carries no ratio_release_over_previous", id="an-entry-without-its-ratio",
    ),
    pytest.param(
        lambda r: drop_keys(same_day(r)["per_dataset"][4], "n_cells", "label"),
        f"{SAME_DAY}.per_dataset[4] carries no label, n_cells", id="an-entry-without-a-key-the-caption-does-not-read",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][3].update(dataset="clinvar_grch38"),
        f"{SAME_DAY}.per_dataset[3] names the dataset 'clinvar_grch38', no population dataset of the whole-genome table",
        id="a-dataset-the-table-does-not-carry",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][3].update(dataset=["1000g_phase3"]),
        f"{SAME_DAY}.per_dataset[3] names the dataset ['1000g_phase3'], no population dataset of the whole-genome table",
        id="a-dataset-that-is-not-a-string",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][1].update(arch="riscv64"),
        f"{SAME_DAY}.per_dataset[1] names the architecture 'riscv64'; one of arm64, x86_64 is read",
        id="an-architecture-other-than-the-two",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][1].update(arch={"name": "x86_64"}),
        f"{SAME_DAY}.per_dataset[1] names the architecture {{'name': 'x86_64'}}; one of arm64, x86_64 is read",
        id="an-architecture-that-is-not-a-string",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][1].update(arch="arm64"),
        f"{SAME_DAY}.per_dataset carries gnomad_v4.1_genomes on arm64 twice", id="an-entry-given-twice",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][0].update(ratio_release_over_previous=0),
        f"{SAME_DAY}.per_dataset[0] carries a ratio_release_over_previous that is not a positive number: 0",
        id="a-ratio-of-zero",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][7].update(ratio_release_over_previous=-0.926),
        f"{SAME_DAY}.per_dataset[7] carries a ratio_release_over_previous that is not a positive number: -0.926",
        id="a-ratio-below-zero",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][0].update(ratio_release_over_previous="n/a"),
        f"{SAME_DAY}.per_dataset[0] carries a ratio_release_over_previous that is not a positive number: 'n/a'",
        id="a-ratio-that-is-not-a-number",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][0].update(ratio_release_over_previous=float("nan")),
        f"{SAME_DAY}.per_dataset[0] carries a ratio_release_over_previous that is not a positive number: nan",
        id="a-ratio-that-is-not-finite",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"][0].update(ratio_release_over_previous=1e30),
        f"{SAME_DAY}.per_dataset[0] carries a ratio_release_over_previous too large to state at two decimals: 1e+30",
        id="a-ratio-too-large-to-round",
    ),
    pytest.param(
        lambda r: same_day(r).update(per_dataset=same_day(r)["per_dataset"][:6]),
        f"{SAME_DAY}.per_dataset carries no entry for 1000g_highcov on arm64, x86_64", id="a-table-dataset-with-no-entry",
    ),
    pytest.param(
        lambda r: same_day(r)["per_dataset"].pop(5),
        f"{SAME_DAY}.per_dataset carries no entry for 1000g_phase3 on x86_64", id="a-table-dataset-on-one-architecture-only",
    ),
]


@pytest.mark.parametrize("edit, expected", SAME_DAY_PROBES)
def test_a_malformed_same_day_block_is_refused_by_name_before_anything_is_printed(tmp_path: Path, edit, expected: str) -> None:
    repo = six_dataset_repo(tmp_path)
    edit_record(repo, "population", lambda r: r.update(previous_release_wall_time=six_dataset_same_day_block()))
    edit_record(repo, "population", edit)
    out = render(repo)
    assert (out.returncode, out.stdout, out.stderr) == (1, "", f"{ERR}{expected}\n")


def test_a_chromosome_peak_rss_that_is_absent_or_not_a_finite_number_is_refused_by_name(tmp_path: Path) -> None:
    repo = six_dataset_repo(tmp_path)
    where = f"{ERR}the population record's wall_time.vep_rs.per_cell[0] carries"
    for edit, expected in (
        (lambda r: drop_keys(timed_cell(r, 0), "peak_rss_mb_median"), f"{where} no peak_rss_mb_median"),
        (lambda r: timed_cell(r, 0).update(peak_rss_mb_median=None), f"{where} no peak_rss_mb_median"),
        (lambda r: timed_cell(r, 0).update(peak_rss_mb_median="n/a"), f"{where} a peak_rss_mb_median that is not a number: 'n/a'"),
        (lambda r: timed_cell(r, 0).update(peak_rss_mb_median=float("nan")), f"{where} a peak_rss_mb_median that is not finite: nan"),
        (lambda r: timed_cell(r, 0).update(peak_rss_mb_median=float("inf")), f"{where} a peak_rss_mb_median that is not finite: inf"),
    ):
        edit_record(repo, "population", edit)
        assert error(repo) == f"{expected}\n"
    # The earlier page builds its caption before its dataset rows, so an empty cell list reaches peak_rss_gb before
    # the rows refuse it; the six-dataset page refuses it as a dataset without wall time first.
    repo = fixture_repo(tmp_path / "empty")
    edit_record(repo, "population", lambda r: r["wall_time"]["vep_rs"]["per_cell"].clear())
    assert error(repo) == f"{ERR}the population record's wall_time.vep_rs.per_cell is empty\n"


# The published 0.4.0 release: the tag's commit, the SHA256SUMS digests of its assets, and the concordance figures
# of its page as the renderer states them from the records the repository carries.
COMMIT_040 = "f64272326dee8f813a57f9fa832539ccd8067dab"
DIGESTS_040 = {
    "vep-0.4.0-aarch64-apple-darwin.tar.gz": "b22283552d1e4c995090025c04dedbc884f1b1afc477ee9f88334f9bb4674fc2",
    "vep-0.4.0-aarch64-unknown-linux-gnu.tar.gz": "6057da61fb971e4ff822b1fe5c5b2e6c7f251fadedf2737f9b94ac2c195537d0",
    "vep-0.4.0-x86_64-unknown-linux-gnu.tar.gz": "d72ff534447dd38ceef576c7d51d479cb7af354c09bba42c47f0991343bcb3a8",
    "vep-rs-0.4.0.tar.gz": "0ff936860386d74d5f532ef86d0b19f4c31eadc90998062d0291c5bc402756ab",
}
PUBLISHED_040_CAPTION = (
    "Whole genome against Ensembl VEP 116.2: the population datasets one run per chromosome, each ClinVar release one "
    "genome-wide file; wall time is the sum of the per-chromosome medians of 20 machines at 16 threads on `c8gd.8xlarge` "
    "(ARM) and `c8id.8xlarge` (x86) for the population datasets, and the whole-file median of 20 machines on "
    "`r8gd.8xlarge` (ARM) and `r8id.8xlarge` (x86) for the ClinVar files, another instance class and protocol, so the "
    "`All six` row carries no wall time; peak memory is each dataset's largest per-chromosome median, at most 2.6 GB on "
    "any chromosome, and a whole-file cell's own median where its campaign recorded one; tuple counts, per-chromosome rows "
    f"and per-file counts are in the [provenance records]({REPO_URL}/tree/v0.4.0/docs/concordance-provenance)."
)
# The peak columns are the per-dataset maxima of the 0.4.0 population record's cells (gnomAD v4.1 2,502.7 and
# 2,316.9 MB; the ClinVar aggregates of that campaign carry no peak), rounded up to one decimal.
PUBLISHED_040_TABLES = """\
| Dataset | Raw F1 | Adjusted F1 | ARM (s) | x86 (s) | ARM peak (GB) | x86 peak (GB) |
| --- | --- | --- | --- | --- | --- | --- |
| gnomAD v4.1 | 0.999999353 | 1.000000000 | 881 | 1,365 | 2.6 | 2.4 |
| gnomAD v2.1.1 | 0.999999521 | 1.000000000 | 282 | 472 | 1.0 | 1.0 |
| 1000 Genomes Phase 3 | 0.999999714 | 1.000000000 | 49 | 76 | 1.0 | 1.0 |
| 1000 Genomes high-coverage | 0.999991204 | 1.000000000 | 85 | 116 | 1.8 | 1.8 |
| ClinVar GRCh37 | 0.999979065 | 1.000000000 | 5 | 6 | | |
| ClinVar GRCh38 | 0.999972446 | 1.000000000 | 14 | 14 | | |
| All six | 0.999998408 | 1.000000000 | | | | |

Structural-variant file sets measured on the released binary against the same reference, each set's files pooled by \
summed counts; the paper's own figures stay as published in \
[docs/published-figures.md](https://github.com/natera-open-source/vep-rs/blob/v0.4.0/docs/published-figures.md).

| File set | Assembly | Raw F1 | Adjusted F1 |
| --- | --- | --- | --- |
| ClinVar SV | GRCh37 | 0.999952183 | 1.000000000 |
| gnomAD SV | GRCh37 | 0.726719694 | 1.000000000 |
| 1000 Genomes SV | GRCh37 | 1.000000000 | 1.000000000 |
| Synthetic (13 files) | GRCh37 | 0.986694682 | 1.000000000 |
| ClinVar SV | GRCh38 | 1.000000000 | 1.000000000 |
| gnomAD SV | GRCh38 | 0.791428005 | 1.000000000 |
| 1000 Genomes SV | GRCh38 | 0.989455966 | 1.000000000 |
| Synthetic (13 files) | GRCh38 | 0.983528044 | 1.000000000 |"""


def release_040_repo(tmp_path: Path) -> Path:
    """The repository's own 0.4.0 inputs in a scratch repository: the release-notes file, the two records,
    CHANGELOG.md, CITATION.cff and Cargo.toml as they stand, and the published digests."""
    repo = tmp_path / "repo"
    (repo / "docs" / "concordance-provenance").mkdir(parents=True)
    (repo / ".github" / "release-notes").mkdir(parents=True)
    for rel in ("CHANGELOG.md", "CITATION.cff", "Cargo.toml", ".github/release-notes/v0.4.0.md"):
        shutil.copy(ROOT / rel, repo / rel)
    for record in (ROOT / "docs" / "concordance-provenance").glob("*-v0.4.0.json"):
        shutil.copy(record, repo / "docs" / "concordance-provenance" / record.name)
    (repo / "SHA256SUMS").write_text("".join(f"{d}  {n}\n" for n, d in DIGESTS_040.items()))
    write_floors(repo, {n: "2.39" for n in DIGESTS_040 if "linux-gnu" in n})
    return repo


def test_the_0_4_0_records_render_the_published_concordance_section(tmp_path: Path) -> None:
    """The repository's 0.4.0 records, timed on one generation per pair with the ClinVar files on another class
    than the whole genome, render the published page's concordance figures in the six-dataset page's shape: the
    caption without a label, then the two tables."""
    repo = release_040_repo(tmp_path)
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "--repo", str(repo), "--version", "0.4.0", "--previous", "v0.3.2", "--commit", COMMIT_040,
         "--toolchain", TOOLCHAIN, "--sha256sums", str(repo / "SHA256SUMS"), "--glibc-floors", str(repo / "glibc-floors.json")],
        capture_output=True, text=True, check=False,
    )
    assert out.returncode == 0, out.stderr
    assert concordance_section(out.stdout) == f"{PUBLISHED_040_CAPTION}\n\n{PUBLISHED_040_TABLES}"
    assert NEWEST_LABEL not in out.stdout
    # A flat notes file, Highlights and Upgrade notes without an area section, keeps the published page's shape.
    assert [l for l in out.stdout.splitlines() if l.startswith("#")] == [
        "## Highlights", "## Upgrade notes", "## Added", "## Changed", "## Fixed", "## Concordance and wall time", "## Install", "## Cite",
    ]


PLUGIN_RECORD = ROOT / "scripts" / "docs" / "testdata" / "plugin_record_fixture.json"
PLUGINS = "the release record's plugin_concordance"
PLUGIN_COUNTS = ("perl", "rust", "intersection", "adj_perl", "adj_rust", "adj_intersection")
PLUGIN_DATA = {"clinvar": "ClinVar", "synthetic_fixture": "synthetic fixture", "committed_table": "committed table"}
PLUGIN_CAPTION = (
    "Each supported plugin is measured on the released binary against the same Ensembl VEP reference, as raw and "
    "adjusted F1 over the plugin's populated fields. A plugin verified on real data is measured on the genome-wide "
    "ClinVar dataset of each assembly; a plugin whose data cannot be redistributed, on synthetic fixtures; a plugin "
    "that reads the user's own table, on a committed table. An empty cell is an assembly the plugin is not supported on."
)
NOT_PORTED = 56
PLUGIN_HEADER = "| Plugin | Data | GRCh37 raw F1 | GRCh37 adjusted F1 | GRCh38 raw F1 | GRCh38 adjusted F1 |"


def plugin_record() -> dict:
    return json.loads(PLUGIN_RECORD.read_text())


def plugin_block() -> dict:
    """The `plugin_concordance` block of a release record, built from the fixture plugin record: per gate row
    vep-rs's per-feature counts summed over the row's expected sets, vep-rs's support, and the labels and
    counts of the groups of plugins not ported."""
    record = plugin_record()
    entries = [
        e for name in ("real_data", "synthetic_fixtures", "committed_table") for e in record[name]
        if e["engine"] == "vep-rs" and e["projection"] == "feature"
    ]
    rows = []
    for row in record["rows"]:
        key = (row["plugin"], row["assembly"], row["dataset_or_corpus"])
        sets = [e for e in entries if (e["plugin"], e["assembly"], e.get("dataset", e.get("corpus"))) == key]
        assert sets, key
        masked: dict[str, int] = {}
        for e in sets:
            for name, n in e["masked"].items():
                masked[name] = masked.get(name, 0) + n
        rows.append({
            **row, "expected_sets": len(sets), **{k: sum(e[k] for e in sets) for k in PLUGIN_COUNTS},
            "masked": masked, "residual_rows": sum(e["residual_rows"] for e in sets),
        })
    support = next(s for s in record["support"] if s["engine"] == "vep-rs")
    groups = [{"label": g["label"], "count": g["count"]} for g in record["not_ported"]["groups"]]
    return {
        "plugins_record": {"file": PLUGIN_RECORD.name, "sha256": hashlib.sha256(PLUGIN_RECORD.read_bytes()).hexdigest()},
        "ensembl_vep": "116.2",
        "rows": rows,
        "support": {k: v for k, v in support.items() if k != "engine"},
        "not_ported": {"of": record["not_ported"]["of"], "groups": groups},
    }


def plugin_repo(tmp_path: Path, edit=None) -> Path:
    """The six-dataset scratch repository whose release record carries the fixture block, after `edit`."""
    repo = six_dataset_repo(tmp_path)
    block = plugin_block()
    if edit is not None:
        edit(block)
    edit_record(repo, "release", lambda r: r.update(plugin_concordance=block))
    return repo


def expected_plugin_table(block: dict) -> list[str]:
    """The plugin table over a block: one row per label and kind of data in the block's order, each figure from
    the row's integer counts with the fractions module."""
    cells: dict[tuple[str, str], dict[str, str]] = {}
    for row in block["rows"]:
        raw = f1(row["intersection"], row["perl"], row["rust"])
        adjusted = f1(row["adj_intersection"], row["adj_perl"], row["adj_rust"])
        cells.setdefault((row["plugin"], row["scored_on"]), {})[row["assembly"]] = f" {raw} | {adjusted} |"
    table = [PLUGIN_HEADER, "| --- | --- | --- | --- | --- | --- |"]
    for (label, scored_on), pairs in cells.items():
        table.append(f"| {label} | {PLUGIN_DATA[scored_on]} |" + pairs.get("GRCh37", " | |") + pairs.get("GRCh38", " | |"))
    return table


def test_a_plugin_block_renders_the_plugin_table_and_the_plugins_not_ported_after_the_structural_variant_table(tmp_path: Path) -> None:
    block = plugin_block()
    assert len(block["rows"]) == 57
    repo = plugin_repo(tmp_path)
    body = page(repo)
    section = concordance_section(body)
    six_dataset = "\n\n".join(expected_six_dataset_section())
    table = expected_plugin_table(block)
    labels = [g["label"] for g in block["not_ported"]["groups"]]
    reasons = [
        "| Reason | Count |", "| --- | --- |",
        *(f"| {label} | {n} |" for label, n in zip(labels, (16, 13, 7, 7, 5, 4, 4))),
        "| Total | 56 |",
    ]
    sentence = (
        "vep-rs supports 31 of the 87 Ensembl VEP plugins in its support table. The 56 it does not port are counted "
        "by reason below; the README names each one."
    )
    assert sum(g["count"] for g in block["not_ported"]["groups"]) == NOT_PORTED and "count" not in block["not_ported"]
    assert section == "\n\n".join([six_dataset, PLUGIN_CAPTION, "\n".join(table), sentence, "\n".join(reasons)])

    rows = table[2:]
    assert len(rows) == 32
    names = [r.split(" | ")[0][2:] for r in rows]
    assert len(set(names)) == 32 and len({n.partition(" (")[0] for n in names}) == 31
    assert [n for n in names if " (" in n] == ["CADD (structural variants)"]
    data = [r.split(" | ")[1] for r in rows]
    assert (data.count("ClinVar"), data.count("synthetic fixture"), data.count("committed table")) == (20, 11, 1)
    by_name = dict(zip(names, rows))
    assert by_name["AlphaMissense"] == "| AlphaMissense | ClinVar | 0.999888723 | 1.000000000 | 0.999997449 | 1.000000000 |"
    assert by_name["Conservation"].startswith("| Conservation | ClinVar | | | ")
    assert by_name["CADD (structural variants)"].startswith("| CADD (structural variants) | synthetic fixture | | | ")
    assert by_name["TranscriptAnnotator"].startswith("| TranscriptAnnotator | committed table | ")
    figures = [cell for r in rows for cell in r.strip("| ").split(" | ")[2:] if cell.strip("| ")]
    adjusted = [r.split("|")[i].strip() for r in rows for i in (4, 6) if r.split("|")[i].strip()]
    assert len(adjusted) == 57 and set(adjusted) == {"1.000000000"}
    assert all(re.fullmatch(r"[01]\.\d{9}", cell) for cell in figures)
    assert len({r.split("|")[i].strip() for r in rows for i in (3, 5)} - {"", "1.000000000"}) > 5, "raw F1 comes from the counts"

    assert [l for l in body.splitlines() if l.startswith("#")] == [
        "## Upgrade notes", "## Changed", "## Concordance and wall time", "## Install", "## Cite",
    ]

    # A block that carries more of the plugin record than the page reads renders the same page: a row's rounded
    # figures, sets and classes, the record's file, a stated total, each group's plugins and the other engines.
    record = plugin_record()
    unread_names = [p for g in record["not_ported"]["groups"] for p in g["plugins"]]
    unread_names += [text for e in record["engines"] if e["id"] != "vep-rs" for text in (e["id"], e["name"])]

    def carry_unread_fields(release: dict) -> None:
        carried = release["plugin_concordance"]
        for row in carried["rows"]:
            row.update(raw_f1=0.1, adjusted_f1=0.2, expected_sets=99, masked={"a_class": 3})
        carried["plugins_record"] = {"file": "another.json", "sha256": "0" * 64}
        carried["not_ported"]["count"] = 999
        for group, source in zip(carried["not_ported"]["groups"], record["not_ported"]["groups"]):
            group["plugins"] = source["plugins"]
        carried.update(engines=record["engines"], absent=record["absent"])

    edit_record(repo, "release", carry_unread_fields)
    carried = json.dumps(release_record(repo)["plugin_concordance"])
    assert len(unread_names) == NOT_PORTED + 2 and all(name in carried for name in unread_names)
    assert page(repo) == body
    assert not [name for name in unread_names if name in body]
    edit_record(repo, "release", lambda r: r.pop("plugin_concordance"))
    assert concordance_section(page(repo)) == six_dataset


def test_the_plugin_tables_add_fifty_lines_to_the_fixed_lines_of_the_page(tmp_path: Path) -> None:
    """The page's fixed lines run from the Concordance header through the footer: neither the notes file nor the
    CHANGELOG writes one. On the fixture repository, which lists three binary archives, they are 42 without the
    block and 92 with the 57 gate rows: the caption, the 34 lines of the plugin table, the sentence, the 10
    lines of the table of plugins not ported, and the blank line before each of the four."""

    def fixed_lines(repo: Path) -> list[str]:
        lines = page(repo).splitlines()
        return lines[lines.index("## Concordance and wall time") :]

    repo = plugin_repo(tmp_path)
    with_tables = fixed_lines(repo)
    edit_record(repo, "release", lambda r: r.pop("plugin_concordance"))
    without = fixed_lines(repo)
    assert (len(without), len(with_tables)) == (42, 92)
    install = with_tables.index("## Install")
    added = with_tables[install - 50 : install - 1]
    assert with_tables[install - 51] == ""
    assert [len(block.splitlines()) for block in "\n".join(added).split("\n\n")] == [1, 34, 1, 10]
    assert with_tables[: len(without) - 15] + with_tables[-15:] == without[:-15] + without[-15:]


def test_a_block_of_fewer_kinds_of_data_renders_its_own_rows_and_numbers_under_the_same_caption(tmp_path: Path) -> None:
    def keep(*kinds: str):
        def edit(block: dict) -> None:
            block["rows"] = [row for row in block["rows"] if row["scored_on"] in kinds]
            scored = len({row["plugin"].partition(" (")[0] for row in block["rows"]})
            block["support"]["read"] = scored
            block["not_ported"]["of"] = scored + NOT_PORTED
        return edit

    def plugin_part(repo: Path) -> list[str]:
        """The four blocks after the structural-variant table: the caption, the plugin table, the sentence and the reasons."""
        return concordance_section(page(repo)).split("\n\n")[4:]

    caption, table, sentence, _ = plugin_part(plugin_repo(tmp_path / "one", keep("clinvar")))
    assert caption == PLUGIN_CAPTION
    assert len(table.splitlines()) == 2 + 20 and {r.split(" | ")[1] for r in table.splitlines()[2:]} == {"ClinVar"}
    assert sentence.startswith("vep-rs supports 20 of the 76 Ensembl VEP plugins in its support table. The 56 it does not port ")
    caption, table, sentence, _ = plugin_part(plugin_repo(tmp_path / "two", keep("committed_table", "synthetic_fixture")))
    assert caption == PLUGIN_CAPTION
    assert len(table.splitlines()) == 2 + 12
    assert {r.split(" | ")[1] for r in table.splitlines()[2:]} == {"synthetic fixture", "committed table"}
    assert sentence.startswith("vep-rs supports 11 of the 67 Ensembl VEP plugins in its support table. The 56 it does not port ")


def test_an_earlier_shape_record_carrying_a_plugin_block_is_refused(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    edit_record(repo, "release", lambda r: r.update(plugin_concordance=plugin_block()))
    assert error(repo) == (
        f"{ERR}{PLUGINS} is rendered after the structural-variant table, which a record of the earlier shape lacks\n"
    )


def plugin_row(block: dict, plugin: str, assembly: str) -> dict:
    return next(row for row in block["rows"] if (row["plugin"], row["assembly"]) == (plugin, assembly))


SV_ROW = "rows[11] (CADD (structural variants) on GRCh38)"
PLUGIN_PROBES = [  # the block's edit, then the exact ERROR line
    pytest.param(
        lambda b: b["not_ported"]["groups"][0].update(count=15),
        f"{ERR}{PLUGINS} counts 31 plugins supported and 55 not ported, which are not the 87 it lists",
        id="groups-whose-sum-and-the-supported-count-do-not-give-the-plugins-listed",
    ),
    pytest.param(
        lambda b: b["not_ported"].update(of=88),
        f"{ERR}{PLUGINS} counts 31 plugins supported and 56 not ported, which are not the 88 it lists",
        id="supported-and-not-ported-that-do-not-sum-to-the-plugins-listed",
    ),
    pytest.param(
        lambda b: b["support"].update(read=32),
        f"{ERR}{PLUGINS} counts 32 plugins supported and 56 not ported, which are not the 87 it lists",
        id="a-supported-count-that-does-not-sum",
    ),
    pytest.param(
        lambda b: (b["support"].update(read=30), b["not_ported"].update(of=86)),
        f"{ERR}{PLUGINS}.support reads 30 plugins and the rows score 31",
        id="a-supported-count-that-is-not-the-plugins-the-rows-score",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").update(adj_rust=269569),
        f"{ERR}{PLUGINS}.rows[0] (AlphaMissense on GRCh37) is not a perfect adjusted match: adj_perl 269568, adj_rust 269569, "
        "adj_intersection 269568",
        id="adjusted-counts-that-differ-on-one-side",
    ),
    pytest.param(
        lambda b: plugin_row(b, "CADD (structural variants)", "GRCh38").update(adj_intersection=0),
        f"{ERR}{PLUGINS}.{SV_ROW} is not a perfect adjusted match: adj_perl 157, adj_rust 157, adj_intersection 0",
        id="an-adjusted-intersection-below-both-sides",
    ),
    pytest.param(
        lambda b: plugin_row(b, "CADD (structural variants)", "GRCh38").update(residual_rows=2),
        f"{ERR}{PLUGINS}.{SV_ROW} carries 2 residual rows; the page states a plugin with none",
        id="a-row-with-residual-rows",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").pop("residual_rows"),
        f"{ERR}{PLUGINS}.rows[0] carries no residual_rows",
        id="a-row-without-its-residual-count",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").update(perl=269598.0),
        f"{ERR}{PLUGINS}.rows[0] (AlphaMissense on GRCh37) carries perl that is not an integer: 269598.0",
        id="a-count-that-is-not-an-integer",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").pop("adj_perl"),
        f"{ERR}{PLUGINS}.rows[0] (AlphaMissense on GRCh37) carries no adj_perl",
        id="a-row-without-an-adjusted-count",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").update(scored_on="real_data"),
        f"{ERR}{PLUGINS}.rows[0] (AlphaMissense on GRCh37) is scored on 'real_data'; one of clinvar, synthetic_fixture, "
        "committed_table is read",
        id="a-kind-of-data-the-table-has-no-cell-for",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").update(scored_on=["clinvar"]),
        f"{ERR}{PLUGINS}.rows[0] (AlphaMissense on GRCh37) is scored on ['clinvar']; one of clinvar, synthetic_fixture, "
        "committed_table is read",
        id="a-kind-of-data-that-is-not-text",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").update(assembly="GRCm39"),
        f"{ERR}{PLUGINS}.rows[0] (AlphaMissense on GRCm39) names an assembly the table has no columns for; one of GRCh37, "
        "GRCh38 is read",
        id="an-assembly-the-table-has-no-columns-for",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").update(assembly="GRCh38"),
        f"{ERR}{PLUGINS}.rows carries AlphaMissense on GRCh38 twice",
        id="a-plugin-on-one-assembly-twice",
    ),
    pytest.param(
        lambda b: plugin_row(b, "AlphaMissense", "GRCh37").update(plugin="Alpha | Missense"),
        f"{ERR}{PLUGINS}.rows[0] carries plugin that is not the text of one table cell: 'Alpha | Missense'",
        id="a-label-that-would-break-its-row",
    ),
    pytest.param(
        lambda b: b["not_ported"]["groups"][2].update(label=" "),
        f"{ERR}{PLUGINS}.not_ported.groups[2] carries label that is not the text of one table cell: ' '",
        id="a-group-with-a-blank-label",
    ),
    pytest.param(
        lambda b: b["not_ported"]["groups"][2].update(count=-7),
        f"{ERR}{PLUGINS}.not_ported.groups[2] carries count that is below zero: -7",
        id="a-group-count-below-zero",
    ),
    pytest.param(
        lambda b: b["not_ported"]["groups"][2].pop("count"),
        f"{ERR}{PLUGINS}.not_ported.groups[2] carries no count",
        id="a-group-without-its-count",
    ),
    pytest.param(
        lambda b: b["not_ported"].update(groups={"label": "one", "count": 56}),
        f"{ERR}{PLUGINS}.not_ported.groups is not a list",
        id="groups-that-are-not-a-list",
    ),
    pytest.param(
        lambda b: b["not_ported"].pop("of"),
        f"{ERR}{PLUGINS}.not_ported carries no of",
        id="a-block-without-the-plugins-listed",
    ),
    pytest.param(
        lambda b: b["not_ported"].update(of=True),
        f"{ERR}{PLUGINS}.not_ported carries of that is not an integer: True",
        id="a-listed-count-that-is-not-an-integer",
    ),
    pytest.param(
        lambda b: b["support"].pop("read"),
        f"{ERR}{PLUGINS}.support carries no read",
        id="a-support-block-without-the-plugins-read",
    ),
    pytest.param(lambda b: b.pop("not_ported"), f"{ERR}{PLUGINS} carries no not_ported", id="a-block-without-not-ported"),
    pytest.param(lambda b: b.pop("rows"), f"{ERR}{PLUGINS} carries no rows", id="a-block-without-rows"),
    pytest.param(lambda b: b.update(rows=[]), f"{ERR}{PLUGINS}.rows is not a list of rows", id="a-block-with-no-rows"),
    pytest.param(
        lambda b: b.update(ensembl_vep="115.2"),
        f"{ERR}{PLUGINS} names Ensembl VEP '115.2'; the records' perl_reference blocks name 116.2",
        id="a-block-naming-another-reference-release",
    ),
    pytest.param(lambda b: b.pop("ensembl_vep"), f"{ERR}{PLUGINS} carries no ensembl_vep", id="a-block-without-its-reference-release"),
]


@pytest.mark.parametrize("edit, expected", PLUGIN_PROBES)
def test_a_malformed_plugin_block_is_refused_by_name_before_anything_is_printed(tmp_path: Path, edit, expected: str) -> None:
    out = render(plugin_repo(tmp_path, edit))
    assert (out.returncode, out.stdout, out.stderr) == (1, "", f"{expected}\n")


def test_a_plugin_block_that_is_not_an_object_is_refused(tmp_path: Path) -> None:
    repo = six_dataset_repo(tmp_path)
    edit_record(repo, "release", lambda r: r.update(plugin_concordance=[plugin_block()]))
    assert error(repo) == f"{ERR}{PLUGINS} is not an object\n"
