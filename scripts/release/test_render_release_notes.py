"""Tests for render_release_notes.py on two input shapes: the repository's own 0.3.2 inputs (the release and
population records, the CHANGELOG section, CITATION.cff, Cargo.toml and the published SHA256SUMS digests),
copied into a scratch repository that each test may edit, and release and population records of the
six-dataset shape that the tests build under the same version. The scratch release-notes file keeps the
repository file's Upgrade notes under an opener within the renderer's budget."""

from __future__ import annotations

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
    for record in (ROOT / "docs" / "concordance-provenance").glob(f"*-v{VERSION}.json"):
        shutil.copy(record, repo / "docs" / "concordance-provenance" / record.name)
    write_notes(repo, OPENER, f"Upgrade notes\n\n{repository_upgrade_notes()}")
    (repo / "SHA256SUMS").write_text("".join(f"{d}  {n}\n" for n, d in DIGESTS.items()))
    return repo


def render(repo: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--repo", str(repo), "--version", VERSION, "--previous", PREVIOUS,
         "--commit", COMMIT, "--toolchain", TOOLCHAIN, "--sha256sums", str(repo / "SHA256SUMS"), *extra],
        capture_output=True, text=True,
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
    msrv = re.search(r'^rust-version = "([^"]+)"', (repo / "Cargo.toml").read_text(), re.MULTILINE).group(1)
    install = [f"- [`{n}`]({DOWNLOAD}/{n}): {PLATFORMS[n]}; sha256 `{DIGESTS[n]}`" for n in PLATFORMS]
    install += [
        f"- Each archive unpacks to `vep-{VERSION}-<target>/` holding `vep`, `vep-cache-builder` and `vep-cache-converter`; "
        "a binary below its CPU floor aborts with an illegal-instruction fault, so on such a host build from source "
        f"(Rust {msrv} or later), which compiles for the machine it runs on: "
        '`RUSTFLAGS="-C target-cpu=native" cargo build --release --locked`',
        f"- Container image: `docker pull ghcr.io/natera-open-source/vep-rs:{VERSION}` (linux/amd64: the x86_64 binaries "
        "with the `duckdb` CLI that Parquet output needs)",
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
        "2.4 GB on any chromosome; tuple counts and per-chromosome rows are in the "
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


def test_the_page_runs_about_sixty_lines_and_is_refused_over_eighty(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    lines = page(repo).splitlines()
    assert len(lines) <= 60
    assert [l for l in lines if l.startswith("#")] == [
        "## Upgrade notes", "## Changed", "## Concordance and wall time", "## Install", "## Cite",
    ]
    assert not [l for l in lines if l.startswith((">", "<"))]
    assert "**Full changelog**" not in "\n".join(lines[:-1])
    assert lines[-1].startswith("**Full changelog**: ")
    write_changelog(repo, "### Changed\n\n" + "".join(f"- Change number {i}.\n" for i in range(45)))
    m = re.search(r"the page is (\d+) lines; the cap is 80", error(repo))
    assert m and int(m.group(1)) > 80


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
        "2.4 GB on any chromosome; tuple counts and per-chromosome rows are in the "
        f"[provenance records]({REPO_URL}/tree/{TAG}/docs/concordance-provenance)."
    ) in body
    assert "-release-v" not in body and "-population-v" not in body

    def raise_peak(record: dict) -> None:
        record["wall_time"]["vep_rs"]["per_cell"][0]["peak_rss_mb_median"] = 2500

    edit_record(repo, "population", raise_peak)
    assert "peak memory at most 2.5 GB on any chromosome" in page(repo)


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
        "prefix": "population_r116.2_c116_20261020",
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
        record["population_concordance"]["perl_reference"] = {"prefix": "population_r116.2_c116_20261020", "cache_version": 116}

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


def test_upgrade_notes_are_one_to_five_bullets_with_breaking_ones_first(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_notes(repo, "A release with nothing to upgrade.")
    body = page(repo)
    assert "## Upgrade notes" not in body
    assert body.startswith("A release with nothing to upgrade.\n\n## Changed\n\n- ")
    write_notes(repo, "A release.", "Upgrade notes\n\n- **Breaking:** `--old` is gone; pass `--new`.\n* A cache converted\n  for 0.3.0 is valid.")
    assert "A release.\n\n## Upgrade notes\n\n- **Breaking:** `--old` is gone; pass `--new`.\n* A cache converted for 0.3.0 is valid.\n\n## Changed\n" in page(repo)
    write_notes(repo, "A release.", "Upgrade notes\n\nNothing breaks.\n- A note.")
    assert "'## Upgrade notes' in" in error(repo) and "is bullets only" in error(repo)
    write_notes(repo, "A release.", "Upgrade notes\n\n" + "".join(f"- Note {i}.\n" for i in range(6)))
    assert "has 6 bullets; one to five" in error(repo)
    write_notes(repo, "A release.", "Upgrade notes\n\n- A note.\n- **Breaking:** `--old` is gone.")
    assert "lists a **Breaking:** bullet after one that is not" in error(repo)


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
    assert install[3].startswith("- Each archive unpacks to `vep-0.3.2-<target>/`") and "Rust 1.88 or later" in install[3]
    assert install[4] == (
        "- Container image: `docker pull ghcr.io/natera-open-source/vep-rs:0.3.2` (linux/amd64: the x86_64 binaries "
        "with the `duckdb` CLI that Parquet output needs)"
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


def test_the_notes_file_is_an_opener_plus_highlights_and_upgrade_notes_only(tmp_path: Path) -> None:
    repo = fixture_repo(tmp_path)
    write_notes(repo, "A release.", "Known issues\n\n- An issue.")
    assert "unknown section '## Known issues'" in error(repo)


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


def six_dataset_records() -> tuple[dict, dict]:
    """The release and population records of the six-dataset shape, every float field a coarse or wrong
    rounding of the counts beside it."""
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
        {"arch": arch, "suite_id": s, "median_sec": sec, "instance_type": MEMORY[arch], "n_clones": 20}
        for (s, arch), sec in CLINVAR_WALL.items()
    ]
    aggregates += [{"arch": a, "suite_id": s, "median_sec": 1.5, "instance_type": MEMORY[a], "n_clones": 20} for s in SV_CELLS for a in MEMORY]
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
        {"dataset": ds, "label": label, "assembly": asm, "chromosome": str(c), "arch": arch, "instance_type": COMPUTE[arch],
         "n_clones": 20, "median_sec": sec, "peak_rss_mb_median": 2_501.0 if (ds, c, arch) == ("gnomad_v4.1_genomes", 2, "x86_64") else 1_900.0}
        for ds, label, asm, *_ in POPULATION
        for arch, medians in zip(("arm64", "x86_64"), WALL[ds])
        for c, sec in enumerate(medians, 1)
    ]
    population = {
        "release": TAG, "record": "population-measurements", "fork": 16,
        "wall_time": {
            "vep_rs": {"engine": "vep-rs", "engine_version": TAG, "instance_types": COMPUTE, "n_clones_per_cell": 20, "per_cell": per_cell},
            "ensembl_vep": {"engine": "ensembl_vep", "engine_version": "release_115.2", "instance_types": MEMORY},
        },
        "perl_reference": {**REFERENCE, "prefix": "population_r116.2_c116_fixture"},
    }
    return release, population


def six_dataset_repo(tmp_path: Path) -> Path:
    """The 0.3.2 scratch repository with its two records replaced by records of the six-dataset shape under
    the same version: the renderer reads the shape from the record, never from the version."""
    repo = fixture_repo(tmp_path)
    records = repo / "docs" / "concordance-provenance"
    for old in records.glob(f"*-v{VERSION}.json"):
        old.unlink()
    release, population = six_dataset_records()
    (records / f"six-dataset-release-v{VERSION}.json").write_text(json.dumps(release))
    (records / f"six-dataset-population-v{VERSION}.json").write_text(json.dumps(population))
    return repo


def pooled(cells: list[tuple[int, int, int]]) -> str:
    return f1(sum(c[0] for c in cells), sum(c[1] for c in cells), sum(c[2] for c in cells))


def expected_six_dataset_section() -> list[str]:
    raw: list[tuple[int, int, int]] = []
    adjusted: list[tuple[int, int, int]] = []
    rows = ["| Dataset | Raw F1 | Adjusted F1 | ARM (s) | x86 (s) |", "| --- | --- | --- | --- | --- |"]
    for ds, label, _, p, r, i, xp, xr in POPULATION:
        raw.append((i, p, r))
        adjusted.append((i, p - xp, r - xr))
        arm, x86 = (sum(Fraction(str(m)) for m in medians) for medians in WALL[ds])
        rows.append(f"| {label.replace(', whole genome', '')} | {f1(i, p, r)} | {f1(i, p - xp, r - xr)} | {half_up(arm, 0)} | {half_up(x86, 0)} |")
    for s, (label, p, r, i, ap, ar, ai) in CLINVAR.items():
        raw.append((i, p, r))
        adjusted.append((ai, ap, ar))
        arm, x86 = (half_up(Fraction(str(CLINVAR_WALL[(s, arch)])), 0) for arch in ("arm64", "x86_64"))
        rows.append(f"| {label} | {f1(i, p, r)} | {f1(ai, ap, ar)} | {arm} | {x86} |")
    rows.append(f"| All six | {pooled(raw)} | {pooled(adjusted)} | | |")
    sv = ["| File set | Assembly | Raw F1 | Adjusted F1 |", "| --- | --- | --- | --- |"]
    for s, (assembly, gnomad, _) in SV_CELLS.items():
        files = sv_files(s)
        for label, names in (("ClinVar SV", ["clinvar_sv_chr21"]), ("gnomAD SV", [gnomad]), ("1000 Genomes SV", ["1kg_sv_chr21"]),
                             ("Synthetic (13 files)", list(SYNTHETIC_FILES))):
            raw_files = [(files[n].intersection, files[n].perl, files[n].rust) for n in names]
            adjusted_files = [(files[n].adj_intersection, files[n].adj_perl, files[n].adj_rust) for n in names]
            sv.append(f"| {label} | {assembly} | {pooled(raw_files)} | {pooled(adjusted_files)} |")
    caption = (
        "Whole genome against Ensembl VEP 116.2: the population datasets one run per chromosome, each ClinVar release one "
        "genome-wide file; wall time is the sum of the per-chromosome medians of 20 machines at 16 threads on `c8gd.8xlarge` "
        "(ARM) and `c8id.8xlarge` (x86) for the population datasets, and the whole-file median of 20 machines on "
        "`r8gd.8xlarge` (ARM) and `r8id.8xlarge` (x86) for the ClinVar files, another instance class and protocol, so the "
        "`All six` row carries no wall time; peak memory at most 2.5 GB on any chromosome; tuple counts, per-chromosome rows "
        f"and per-file counts are in the [provenance records]({REPO_URL}/tree/{TAG}/docs/concordance-provenance)."
    )
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
    assert rows[4] == "| ClinVar GRCh37 | 0.999979065 | 1.000000000 | 4 | 7 |"
    assert rows[5] == "| ClinVar GRCh38 | 0.999972446 | 1.000000000 | 13 | 21 |", "12.5 s rounds half up, never to even"
    assert rows[6] == "| All six | 0.999998462 | 1.000000000 | | |"
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
    assert all_six != next(l for l in body.splitlines() if l.startswith("| All six |")) and all_six.endswith("| | |")

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


SIX_DATASET_PROBES = [  # record kind, its edit, then the exact ERROR line or a line of the page that renders
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
