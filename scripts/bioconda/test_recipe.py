"""The Bioconda recipe builds the three release binaries from the git-archive tarball."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT = HERE / "check_recipe.py"
sys.path.insert(0, str(HERE))

import check_recipe


def shipped() -> tuple[str, str, str]:
    meta = (ROOT / "bioconda" / "vep-rs" / "meta.yaml").read_text(encoding="utf-8")
    build = (ROOT / "bioconda" / "vep-rs" / "build.sh").read_text(encoding="utf-8")
    version = check_recipe.cargo_version((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    assert version is not None
    return meta, build, version


def recipe_version(meta: str) -> str:
    found = check_recipe.VERSION_SET_RE.search(meta)
    assert found is not None
    return found.group(1)


def scratch_tree(root: Path, workspace: str) -> Path:
    """The shipped recipe beside a Cargo.toml whose workspace version is `workspace`."""
    recipe = root / "bioconda" / "vep-rs"
    recipe.mkdir(parents=True)
    for name in ("meta.yaml", "build.sh"):
        shipped_text = (ROOT / "bioconda" / "vep-rs" / name).read_text(encoding="utf-8")
        (recipe / name).write_text(shipped_text, encoding="utf-8")
    (root / "Cargo.toml").write_text(f'[workspace.package]\nversion = "{workspace}"\n', encoding="utf-8")
    return root


def test_the_shipped_recipe_matches_the_release_contract() -> None:
    meta, build, version = shipped()
    assert check_recipe.recipe_errors(meta, build, version) == []
    assert check_recipe.check_root(ROOT) == []
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(ROOT)], check=False, capture_output=True, text=True
    )
    assert out.returncode == 0, out.stderr
    assert out.stderr == ""


def test_the_source_tarball_digest_is_the_published_release_asset() -> None:
    meta, _, _ = shipped()
    assert recipe_version(meta) == check_recipe.PUBLISHED_VERSION
    assert f"  sha256: {check_recipe.PUBLISHED_SHA256}" in meta.splitlines()


def test_a_workspace_at_or_above_the_published_release_is_accepted(tmp_path: Path) -> None:
    meta, build, _ = shipped()
    major, minor, patch = (int(part) for part in recipe_version(meta).split("."))
    # The published release's own tag, then what a version bump of main, a release
    # branch and the tag cut from it carry.
    for workspace in (
        f"{major}.{minor}.{patch}",
        f"{major}.{minor}.{patch + 1}",
        f"{major}.{minor + 1}.0",
        f"{major + 1}.0.0",
    ):
        assert check_recipe.recipe_errors(meta, build, workspace) == []
        root = scratch_tree(tmp_path / workspace, workspace)
        assert check_recipe.check_root(root) == []
        out = subprocess.run(
            [sys.executable, str(SCRIPT), "--root", str(root)], check=False, capture_output=True, text=True
        )
        assert (out.returncode, out.stderr) == (0, "")


def test_a_version_other_than_the_published_release_is_refused() -> None:
    meta, build, version = shipped()
    recipe = recipe_version(meta)
    older = meta.replace(f'version = "{recipe}"', 'version = "0.0.1"', 1)
    assert check_recipe.recipe_errors(older, build, version) == [
        f"meta.yaml version is 0.0.1, the published release is {recipe}"
    ]


def test_a_version_above_the_workspace_is_refused(tmp_path: Path) -> None:
    meta, build, version = shipped()
    recipe = recipe_version(meta)
    above = [f"meta.yaml version {recipe} is above Cargo.toml's 0.0.1"]
    assert check_recipe.recipe_errors(meta, build, "0.0.1") == above
    assert check_recipe.check_root(scratch_tree(tmp_path, "0.0.1")) == above
    ahead = meta.replace(f'version = "{recipe}"', 'version = "9.9.9"', 1)
    assert check_recipe.recipe_errors(ahead, build, version) == [
        f"meta.yaml version is 9.9.9, the published release is {recipe}",
        f"meta.yaml version 9.9.9 is above Cargo.toml's {version}",
    ]


def test_a_workspace_version_that_is_not_three_numbers_is_refused() -> None:
    meta, build, _ = shipped()
    assert check_recipe.recipe_errors(meta, build, "next") == ["Cargo.toml version next is not MAJOR.MINOR.PATCH"]


def test_a_github_archive_url_or_a_short_sha256_is_refused() -> None:
    meta, build, version = shipped()
    github = meta.replace(
        "releases/download/v{{ version }}/vep-rs-{{ version }}.tar.gz",
        "archive/refs/tags/v{{ version }}.tar.gz",
    )
    assert "source url is not the release git-archive tarball" in check_recipe.recipe_errors(github, build, version)
    short = meta.replace(check_recipe.PUBLISHED_SHA256, "abc")
    assert "source sha256 is not 64 hex characters" in check_recipe.recipe_errors(short, build, version)


def test_omitting_a_binary_or_installing_bench_plugins_is_refused() -> None:
    meta, build, version = shipped()
    dropped = build.replace("  --bin vep-cache-converter\n", "")
    dropped = dropped.replace('  "${release_dir}/vep-cache-converter" \\\n', "")
    errors = check_recipe.recipe_errors(meta, dropped, version)
    assert "build.sh does not build vep-cache-converter" in errors
    assert "build.sh does not install vep-cache-converter" in errors
    extra = build.replace("--bin vep-cache-converter", "--bin vep-cache-converter --bin bench_plugins")
    assert "build.sh installs bench_plugins" in check_recipe.recipe_errors(meta, extra, version)


def test_binaries_install_from_the_cargo_build_target_release_dir() -> None:
    meta, build, version = shipped()
    assert check_recipe.RELEASE_DIR in build.splitlines()
    for target, expected in (
        ("", "target/release"),
        ("x86_64-unknown-linux-gnu", "target/x86_64-unknown-linux-gnu/release"),
    ):
        out = subprocess.run(
            ["bash", "-c", f'{check_recipe.RELEASE_DIR}; printf %s "$release_dir"'],
            check=False,
            capture_output=True,
            text=True,
            env={"PATH": "/usr/bin:/bin", **({"CARGO_BUILD_TARGET": target} if target else {})},
        )
        assert out.returncode == 0, out.stderr
        assert out.stdout == expected
    fixed = build.replace(check_recipe.RELEASE_DIR + "\n", "")
    for binary in check_recipe.BINS:
        fixed = fixed.replace(f'"${{release_dir}}/{binary}"', f"target/release/{binary}")
    errors = check_recipe.recipe_errors(meta, fixed, version)
    assert "build.sh does not take release_dir from CARGO_BUILD_TARGET" in errors
    assert "build.sh names target/release, which CARGO_BUILD_TARGET moves" in errors
    assert all(f"build.sh does not install {binary}" in errors for binary in check_recipe.BINS)


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


def test_a_digest_other_than_the_published_one_is_refused() -> None:
    meta, build, version = shipped()
    digest = check_recipe.SHA256_RE.search(meta)
    assert digest is not None
    zeros = meta.replace(digest.group(1), "0" * 64)
    assert check_recipe.recipe_errors(zeros, build, version) == [
        f"source sha256 is not the published vep-rs-{recipe_version(meta)}.tar.gz digest"
    ]


def test_a_target_cpu_of_native_is_refused() -> None:
    meta, build, version = shipped()
    native = build.replace("target-cpu=${cpu}", "target-cpu=native")
    assert native != build
    errors = check_recipe.recipe_errors(meta, native, version)
    assert "build.sh sets target-cpu=native, not the CPU floor in ${cpu}" in errors
    assert "build.sh does not export RUSTFLAGS with -C target-cpu=${cpu}" in errors


def test_a_commented_out_requirement_is_refused() -> None:
    meta, build, version = shipped()
    for line, error in (
        ("    - {{ compiler('rust') }}\n", "build requirements lack {{ compiler('rust') }}"),
        (
            "    - {{ pin_subpackage('vep-rs', max_pin=\"x.x\") }}\n",
            'run_exports does not pin vep-rs with max_pin="x.x"',
        ),
    ):
        commented = meta.replace(line, "    # " + line[4:])
        assert commented != meta
        assert check_recipe.recipe_errors(commented, build, version) == [error]


def test_rustflags_unset_before_the_build_is_refused() -> None:
    meta, build, version = shipped()
    unset = build.replace("\ncargo build --release", "\nunset RUSTFLAGS\ncargo build --release")
    assert unset != build
    assert check_recipe.recipe_errors(meta, unset, version) == [
        "build.sh changes RUSTFLAGS after its export: unset RUSTFLAGS"
    ]


def test_the_notice_is_a_license_file() -> None:
    meta, build, version = shipped()
    dropped = meta.replace("    - NOTICE\n", "")
    assert dropped != meta
    assert check_recipe.recipe_errors(dropped, build, version) == ["license_file lacks NOTICE"]


def test_the_run_requirements_are_those_of_the_published_release(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    meta, build, version = shipped()
    recipe = recipe_version(meta)
    host = "  host:\n    - openssl  # [linux]\n"
    assert meta.count(host) == 1

    def requiring(*entries: str) -> str:
        """The shipped recipe with these lines as a `run:` list under its requirements."""
        return meta.replace(host, host + "  run:\n" + "".join(f"    {entry}\n" for entry in entries))

    # The published release starts no other program, so a recipe that requires one is refused.
    assert check_recipe.listed(check_recipe.section(meta, "requirements"), "run") == []
    for item in ("duckdb-cli >=1.3.2", "duckdb-cli", "samtools"):
        assert check_recipe.recipe_errors(requiring(f"- {item}"), build, version) == [
            f"run requirements name {item}, which {recipe} does not run"
        ]
    commented = requiring("# - duckdb-cli >=1.3.2")
    assert check_recipe.recipe_errors(commented, build, version) == []
    root = scratch_tree(tmp_path, version)
    (root / "bioconda" / "vep-rs" / "meta.yaml").write_text(requiring("- duckdb-cli >=1.3.2"), encoding="utf-8")
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(root)], check=False, capture_output=True, text=True
    )
    assert (out.returncode, out.stderr) == (
        1,
        f"ERROR: [bioconda] run requirements name duckdb-cli >=1.3.2, which {recipe} does not run\n",
    )
    # Each block holds a run requirement in a shape `listed` reads no entry from.
    unread_key = "run requirements are not written as one `  run:` block list under requirements"
    unread_entry = "run requirements hold {}, which is no `    - item` line"
    for block, error in (
        ("  run: [duckdb-cli]\n", unread_key),
        ("  run:  # programs\n    - duckdb-cli\n", unread_key),
        ('  "run":\n    - duckdb-cli\n', unread_key),
        ("{% if linux %}\n  run:\n    - duckdb-cli\n{% endif %}\n", unread_key),
        ("  run:\n      - duckdb-cli\n", unread_entry.format("- duckdb-cli")),
        ("  run:\n  - duckdb-cli\n", unread_entry.format("- duckdb-cli")),
        ("  run:\n    [duckdb-cli]\n", unread_entry.format("[duckdb-cli]")),
    ):
        unread = meta.replace(host, host + block)
        assert check_recipe.listed(check_recipe.section(unread, "requirements"), "run") == []
        assert check_recipe.recipe_errors(unread, build, version) == [error]
    # A release that starts another program records it, and its recipe lists that program alone.
    monkeypatch.setattr(check_recipe, "PUBLISHED_RUN", ("duckdb-cli >=1.3.2",))
    lack = ["run requirements lack duckdb-cli >=1.3.2"]
    assert check_recipe.recipe_errors(meta, build, version) == lack
    assert check_recipe.recipe_errors(commented, build, version) == lack
    assert check_recipe.recipe_errors(requiring("- duckdb-cli >=1.3.2"), build, version) == []
    above_host = meta.replace(host, "  run:\n    - duckdb-cli >=1.3.2\n" + host)
    assert check_recipe.recipe_errors(above_host, build, version) == []
    assert check_recipe.recipe_errors(requiring("- duckdb-cli >=1.3.2", "- samtools"), build, version) == [
        f"run requirements name samtools, which {recipe} does not run"
    ]


def test_a_cpu_assigned_outside_the_case_block_floors_is_refused() -> None:
    meta, build, version = shipped()
    for line, value in (
        ("cpu=native", "native"),
        ('cpu="native"', "native"),
        ("export cpu=generic", "generic"),
    ):
        late = build.replace("esac\n", f"esac\n{line}\n")
        assert late != build
        assert check_recipe.recipe_errors(meta, late, version) == [
            f"build.sh assigns cpu={value} outside the CPU floors of its case block"
        ]
    arm = build.replace("  *)\n", "  FreeBSD-amd64) cpu=native ;;\n  *)\n")
    assert arm != build
    assert check_recipe.recipe_errors(meta, arm, version) == [
        "build.sh assigns cpu=native outside the CPU floors of its case block"
    ]


def test_a_second_cpu_or_feature_flag_on_the_export_line_is_refused() -> None:
    meta, build, version = shipped()
    export = 'export RUSTFLAGS="-C target-cpu=${cpu}"'
    for flag, error in (
        ("-C target_cpu=native", "build.sh sets target_cpu=native, not the CPU floor in ${cpu}"),
        (
            "-C target-feature=+avx2",
            "build.sh sets target-feature=+avx2, which moves the build off the CPU floor in ${cpu}",
        ),
        (
            "-Ctarget_feature=+sve",
            "build.sh sets target_feature=+sve, which moves the build off the CPU floor in ${cpu}",
        ),
    ):
        added = build.replace(export, f'{export[:-1]} {flag}"')
        assert added != build
        assert check_recipe.recipe_errors(meta, added, version) == [error]


def test_a_cargo_build_under_a_cleared_environment_is_refused() -> None:
    meta, build, version = shipped()
    for prefix, through in (
        ('env -i PATH="${PATH}" ', 'env -i PATH="${PATH}"'),
        ("env -i \\\n", "env -i"),
        ("exec -c ", "exec -c"),
    ):
        cleared = build.replace("\ncargo build --release", f"\n{prefix}cargo build --release")
        assert cleared != build
        assert check_recipe.recipe_errors(meta, cleared, version) == [
            f"build.sh starts cargo build through {through}, not as its own command under the exported RUSTFLAGS"
        ]


def test_a_commented_out_build_command_is_refused() -> None:
    meta, build, version = shipped()
    licenses = "cargo-bundle-licenses --format yaml --output THIRDPARTY.yml\n"
    compile_lines = "cargo build --release --locked \\\n" + "".join(
        f"  --bin {binary}" + (" \\\n" if binary != check_recipe.BINS[-1] else "\n") for binary in check_recipe.BINS
    )
    install_lines = (
        "install -m 0755 \\\n"
        + "".join(f'  "${{release_dir}}/{binary}" \\\n' for binary in check_recipe.BINS)
        + '  "${PREFIX}/bin/"\n'
    )
    for command, errors in (
        (
            licenses,
            ["build.sh does not write THIRDPARTY.yml", "cargo-bundle-licenses does not run before cargo build"],
        ),
        (
            compile_lines,
            [
                "build.sh does not cargo build --release --locked",
                *(f"build.sh does not build {binary}" for binary in check_recipe.BINS),
                "RUSTFLAGS is not exported before cargo build",
                "cargo-bundle-licenses does not run before cargo build",
            ],
        ),
        (install_lines, [f"build.sh does not install {binary}" for binary in check_recipe.BINS]),
    ):
        assert command in build
        commented = build.replace(command, "".join(f"# {line}\n" for line in command.splitlines()))
        assert check_recipe.recipe_errors(meta, commented, version) == errors


def test_a_trailing_comment_that_names_rustflags_is_accepted() -> None:
    meta, build, version = shipped()
    licenses = "cargo-bundle-licenses --format yaml --output THIRDPARTY.yml"
    noted = build.replace(licenses, licenses + "  # RUSTFLAGS does not reach this tool")
    assert noted != build
    assert check_recipe.recipe_errors(meta, noted, version) == []
    # Inside quotes a `#` opens no comment, so this line still names the variable.
    directory = 'install -d "${PREFIX}/bin"'
    quoted = build.replace(directory, directory + '; echo "# RUSTFLAGS"')
    assert quoted != build
    assert check_recipe.recipe_errors(meta, quoted, version) == [
        f'build.sh changes RUSTFLAGS after its export: {directory}; echo "# RUSTFLAGS"'
    ]


def test_a_pre_release_workspace_version_is_read_as_a_version(tmp_path: Path) -> None:
    meta, build, _ = shipped()
    recipe = recipe_version(meta)
    major, minor, patch = (int(part) for part in recipe.split("."))
    for workspace in (f"{major}.{minor + 1}.0-rc.1", f"{major}.{minor}.{patch + 1}-alpha", f"{major + 1}.0.0-0"):
        assert check_recipe.recipe_errors(meta, build, workspace) == []
        out = subprocess.run(
            [sys.executable, str(SCRIPT), "--root", str(scratch_tree(tmp_path / workspace, workspace))],
            check=False,
            capture_output=True,
            text=True,
        )
        assert (out.returncode, out.stderr) == (0, "")
    # A pre-release sorts below its own release, so the recipe's release is above it.
    below = f"{recipe}-rc.1"
    assert check_recipe.recipe_errors(meta, build, below) == [
        f"meta.yaml version {recipe} is above Cargo.toml's {below}"
    ]
    ordered = [
        check_recipe.precedence(text)
        for text in (
            "1.0.0-alpha",
            "1.0.0-alpha.1",
            "1.0.0-alpha.beta",
            "1.0.0-beta",
            "1.0.0-beta.2",
            "1.0.0-beta.11",
            "1.0.0-rc.1",
            "1.0.0",
        )
    ]
    assert None not in ordered
    assert ordered == sorted(set(ordered))
    for text in ("1.2", "1.2.3.4", "v1.2.3", "1.2.3-", "1.2.3-rc..1", "1.2.3-01", "1.2.3-rc_1"):
        assert check_recipe.recipe_errors(meta, build, text) == [f"Cargo.toml version {text} is not MAJOR.MINOR.PATCH"]
