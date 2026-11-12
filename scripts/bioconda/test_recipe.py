"""The Bioconda recipe builds the three release binaries from the git-archive tarball."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT = HERE / "check_recipe.py"
sys.path.insert(0, str(HERE))

import check_recipe  # noqa: E402

# The v0.5.0 release asset vep-rs-0.5.0.tar.gz.sha256.
PUBLISHED_SHA256 = "3ac0ab1e9af04691fa0ed1c19b6a1512180755fc383f4e2fcbfd800b6530926b"


def shipped() -> tuple[str, str, str]:
    meta = (ROOT / "bioconda" / "vep-rs" / "meta.yaml").read_text(encoding="utf-8")
    build = (ROOT / "bioconda" / "vep-rs" / "build.sh").read_text(encoding="utf-8")
    version = check_recipe.cargo_version((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    assert version is not None
    return meta, build, version


def test_the_shipped_recipe_matches_the_release_contract() -> None:
    meta, build, version = shipped()
    assert check_recipe.recipe_errors(meta, build, version) == []
    assert check_recipe.check_root(ROOT) == []
    out = subprocess.run([sys.executable, str(SCRIPT), "--root", str(ROOT)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stderr == ""


def test_the_source_tarball_digest_is_the_published_release_asset() -> None:
    meta, _, version = shipped()
    assert version == "0.5.0"
    assert f"  sha256: {PUBLISHED_SHA256}" in meta.splitlines()


def test_a_version_other_than_cargo_toml_is_refused() -> None:
    meta, build, version = shipped()
    moved = meta.replace(f'version = "{version}"', 'version = "9.9.9"', 1)
    errors = check_recipe.recipe_errors(moved, build, version)
    assert errors == [f"meta.yaml version is 9.9.9, Cargo.toml is {version}"]


def test_a_github_archive_url_or_a_short_sha256_is_refused() -> None:
    meta, build, version = shipped()
    github = meta.replace(
        "releases/download/v{{ version }}/vep-rs-{{ version }}.tar.gz",
        "archive/refs/tags/v{{ version }}.tar.gz",
    )
    assert "source url is not the release git-archive tarball" in check_recipe.recipe_errors(github, build, version)
    short = meta.replace(PUBLISHED_SHA256, "abc")
    assert "source sha256 is not 64 hex characters" in check_recipe.recipe_errors(short, build, version)


def test_omitting_a_binary_or_installing_bench_plugins_is_refused() -> None:
    meta, build, version = shipped()
    dropped = build.replace("  --bin vep-cache-converter\n", "")
    dropped = dropped.replace("  target/release/vep-cache-converter \\\n", "")
    errors = check_recipe.recipe_errors(meta, dropped, version)
    assert "build.sh does not build vep-cache-converter" in errors
    assert "build.sh does not install vep-cache-converter" in errors
    extra = build.replace("--bin vep-cache-converter", "--bin vep-cache-converter --bin bench_plugins")
    assert "build.sh installs bench_plugins" in check_recipe.recipe_errors(meta, extra, version)


def test_an_unlocked_build_or_a_measured_cpu_floor_is_refused() -> None:
    meta, build, version = shipped()
    unlocked = build.replace("--release --locked", "--release")
    errors = check_recipe.recipe_errors(meta, unlocked, version)
    assert "build.sh does not cargo build --release --locked" in errors
    assert "RUSTFLAGS is not exported before cargo build" in errors
    measured = build.replace("cpu=x86-64-v2", "cpu=x86-64-v3")
    assert any("x86-64-v3" in error for error in check_recipe.recipe_errors(meta, measured, version))
    sve = build.replace("cpu=generic", "cpu=neoverse-v1")
    assert any("neoverse" in error for error in check_recipe.recipe_errors(meta, sve, version))


def test_openssl_tests_and_license_are_required() -> None:
    meta, build, version = shipped()
    no_ssl = meta.replace("    - openssl  # [linux]\n", "")
    assert "openssl is not a linux host requirement" in check_recipe.recipe_errors(no_ssl, build, version)
    no_test = meta.replace("    - vep --version\n", "")
    assert any("vep --version" in error for error in check_recipe.recipe_errors(no_test, build, version))
    other = meta.replace("license: Apache-2.0", "license: MIT")
    assert "license is not Apache-2.0" in check_recipe.recipe_errors(other, build, version)
    no_third = meta.replace("    - THIRDPARTY.yml\n", "")
    assert "license_file lacks THIRDPARTY.yml" in check_recipe.recipe_errors(no_third, build, version)
    host_ssl = build.replace(
        'if [ "$(uname -s)" = "Linux" ]; then\n  export OPENSSL_DIR="${PREFIX}"\n  export OPENSSL_NO_VENDOR=1\nfi\n',
        'export OPENSSL_DIR="/usr"\n',
    )
    assert "build.sh does not point Linux OpenSSL at the conda prefix" in check_recipe.recipe_errors(
        meta, host_ssl, version
    )


def test_recipe_maintainers_are_the_codeowners() -> None:
    meta, _, _ = shipped()
    owners = []
    for line in (ROOT / ".github" / "CODEOWNERS").read_text(encoding="utf-8").splitlines():
        parts = line.split("#", 1)[0].split()
        if parts and parts[0] == "*":
            owners = [part[1:] for part in parts[1:] if part.startswith("@")]
    block = meta.split("recipe-maintainers:", 1)[1]
    listed = [
        match.group(1)
        for line in block.splitlines()
        if (match := check_recipe.MAINTAINER_RE.match(line))
    ]
    assert listed == owners and owners
