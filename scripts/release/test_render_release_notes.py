"""Tests for render_release_notes.py on a small fixture repository."""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "render_release_notes.py"
SUITES = ("s01", "s02", "s03", "s04", "s05", "s06", "s07", "s08")


def fixture_repo(tmp_path: Path, *, record: bool = True, cff_version_doi: bool = False) -> Path:
    repo = tmp_path / "repo"
    (repo / "manuscript" / "data").mkdir(parents=True)
    (repo / "docs" / "concordance-provenance").mkdir(parents=True)
    (repo / ".github" / "release-notes").mkdir(parents=True)
    (repo / "Cargo.toml").write_text('[workspace.package]\nversion = "9.9.9"\nrust-version = "1.88"\n')
    (repo / "CHANGELOG.md").write_text(
        "# Changelog\n\n## [Unreleased]\n\n### Changed\n\n- Pending.\n\n"
        "## [9.9.9] - 2031-01-01\n\n### Fixed\n\n- A fix.\n\n### Added\n\n- A thing.\n\n"
        "## [0.1.0] - 2026-09-22\n\n- First release.\n"
    )
    cff = "cff-version: 1.2.0\ntitle: vep-rs\nversion: 9.9.9\ndate-released: 2031-01-01\n"
    if cff_version_doi:
        cff += "doi: 10.5281/zenodo.99999999\n"
    else:
        cff += (
            "doi: 10.5281/zenodo.22837896\nidentifiers:\n"
            "  - type: doi\n    value: 10.5281/zenodo.22837896\n    description: all versions\n"
            "  - type: doi\n    value: 10.5281/zenodo.22837897\n    description: version 0.1.0\n"
        )
    (repo / "CITATION.cff").write_text(cff)
    rows = ["date,engine,engine_version,suite,assembly,instance_type,run_index,wall_time_sec,peak_rss_mb,pct_cpu,major_page_faults,minor_page_faults,cache_state,notes"]
    for arch in ("arm64", "x86_64"):
        for i, s in enumerate(SUITES, start=1):
            rows.append(f"2026-09-20,perl,115.2,x,GRCh37,r8gd.8xlarge,median_of_20,{100.0 * i},,,,,warm,arch={arch};suite_id={s}")
            rows.append(f"2026-09-20,vep-rs,0.1.0,x,GRCh37,r8gd.8xlarge,median_of_20,{2.0 * i},,,,,warm,arch={arch};suite_id={s}")
    (repo / "manuscript" / "data" / "wall_times.csv").write_text("\n".join(rows) + "\n")
    if record:
        conc = {s: {"raw_f1": 0.999979, "adj_f1": 1.0, "perl": 1000 + i, "rust": 1000 + i, "intersection": 999 + i} for i, s in enumerate(SUITES)}
        cells = [
            {"arch": arch, "suite_id": s, "median_sec": 1.0 * i, "p5": 0.9 * i, "p95": 1.1 * i, "n_clones": 20, "engine_version": "v9.9.9"}
            for arch in ("arm64", "x86_64")
            for i, s in enumerate(SUITES, start=1)
        ]
        (repo / "docs" / "concordance-provenance" / "2031-01-01-release-v9.9.9.json").write_text(
            json.dumps({"release": "v9.9.9", "concordance": conc, "per_cell_aggregates": cells})
        )
    (repo / ".github" / "release-notes" / "v9.9.9.md").write_text(
        "# vep-rs 9.9.9\n\nOne-sentence summary.\n\n## Upgrade notes\n\n- An upgrade note.\n\n## Known issues\n\n- A known issue.\n"
    )
    (repo / "SHA256SUMS").write_text(
        "a" * 64 + "  vep-9.9.9-x86_64-unknown-linux-gnu.tar.gz\n"
        + "b" * 64 + "  vep-9.9.9-aarch64-unknown-linux-gnu.tar.gz\n"
        + "c" * 64 + "  vep-9.9.9-aarch64-apple-darwin.tar.gz\n"
        + "d" * 64 + "  vep-rs-9.9.9.tar.gz\n"
    )
    return repo


def render(repo: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--repo", str(repo), "--version", "9.9.9", "--previous", "v0.1.0",
         "--commit", "0123abcd", "--toolchain", "1.97.1", "--sha256sums", str(repo / "SHA256SUMS"), *extra],
        capture_output=True, text=True,
    )


def test_body_follows_the_template(tmp_path: Path) -> None:
    out = render(fixture_repo(tmp_path))
    assert out.returncode == 0, out.stderr
    body = out.stdout
    lines = body.splitlines()
    assert lines[0] == "One-sentence summary."
    headings = [l for l in lines if l.startswith("## ")]
    assert headings == [
        "## Upgrade notes", "## Changes", "## Concordance and wall time for this release",
        "## Install", "## Verify", "## Source archive and DOI", "## Citing", "## Known issues",
    ]
    # the CHANGELOG section of the version, verbatim, without its heading
    assert "### Fixed\n\n- A fix.\n\n### Added\n\n- A thing." in body
    assert "- Pending." not in body
    assert "**Full Changelog**: https://github.com/natera-open-source/vep-rs/compare/v0.1.0...v9.9.9" in body
    # the record's directory is linked at the tag; the dated file name is not printed
    assert "(https://github.com/natera-open-source/vep-rs/tree/v9.9.9/docs/concordance-provenance)" in body
    assert "2031-01-01-release-v9.9.9.json" not in body
    # eight concordance rows with grouped tuple counts
    assert "| ClinVar full | GRCh37 | 0.999979 | 1.000000 | 1,000 | 1,000 | 999 |" in body
    assert body.count("| 1.000000 |") == 8
    # wall time: the medians, and the geomean of perl/release over the six SNP/indel suites
    assert "Wall time, median of 20 independent machines per cell" in body
    assert "ClinVar GRCh37 1.00 s ARM / 1.00 s x86" in body
    geo = math.exp(sum(math.log(100.0 * i / (1.0 * i)) for i in range(1, 7)) / 6)
    assert f"{round(geo):d}× ARM" in body  # 100×
    # install rows carry the digests of SHA256SUMS
    assert "| `" + "a" * 64 + "` |" in body
    assert "vep-rs-9.9.9.tar.gz` is `git archive --format=tar.gz --prefix=vep-rs-9.9.9/ v9.9.9`, sha256 `" + "d" * 64 + "`" in body
    assert "from tag `v9.9.9` (commit `0123abcd`) with Rust 1.97.1" in body
    assert "From source (Rust 1.88 or newer)" in body
    # no deposit of its own: the archived 0.1.0 deposit and the concept DOI are named
    assert "This version has no Zenodo deposit of its own; the archived version is 0.1.0, [10.5281/zenodo.22837897]" in body
    assert "Porter M, Borkowski R. vep-rs, version 9.9.9. 2031. doi:10.5281/zenodo.22837896." in body
    assert body.rstrip().endswith("- A known issue.")


def test_version_doi_from_the_top_level_doi(tmp_path: Path) -> None:
    out = render(fixture_repo(tmp_path, cff_version_doi=True))
    assert out.returncode == 0, out.stderr
    assert "Version DOI [10.5281/zenodo.99999999](https://doi.org/10.5281/zenodo.99999999)." in out.stdout
    assert "doi:10.5281/zenodo.99999999." in out.stdout


def test_dry_run_renders_the_unreleased_section_with_a_given_date(tmp_path: Path) -> None:
    out = render(fixture_repo(tmp_path), "--changelog-section", "Unreleased", "--date", "2031-01-02")
    assert out.returncode == 0, out.stderr
    assert "- Pending." in out.stdout and "- A fix." not in out.stdout


def test_missing_record_summary_or_section_fail(tmp_path: Path) -> None:
    assert render(fixture_repo(tmp_path, record=False)).returncode == 1
    repo = fixture_repo(tmp_path / "b")
    (repo / ".github" / "release-notes" / "v9.9.9.md").write_text("# vep-rs 9.9.9\n\n## Known issues\n\n- x\n")
    assert "needs a summary paragraph" in render(repo).stderr
    (repo / ".github" / "release-notes" / "v9.9.9.md").write_text("# vep-rs 9.9.9\n\nSummary.\n\n## Extra\n\n- x\n")
    assert "unknown section" in render(repo).stderr
    repo = fixture_repo(tmp_path / "c")
    assert "no '## [0.3.0]' section" in render(repo, "--changelog-section", "0.3.0").stderr


@pytest.mark.parametrize("name", ["vep-9.9.9-aarch64-apple-darwin.tar.gz"])
def test_install_table_lists_every_target(tmp_path: Path, name: str) -> None:
    body = render(fixture_repo(tmp_path)).stdout
    assert f"[{name}](https://github.com/natera-open-source/vep-rs/releases/download/v9.9.9/{name})" in body
