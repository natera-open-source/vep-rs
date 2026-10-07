"""Tests for render_release_notes.py on the repository's own 0.3.2 inputs: the release and population
records, the CHANGELOG section, CITATION.cff, Cargo.toml and the published SHA256SUMS digests, copied into
a scratch repository that each test may edit. The scratch release-notes file keeps the repository file's
Upgrade notes under an opener within the renderer's budget."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

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
