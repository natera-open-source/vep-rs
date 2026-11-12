#!/usr/bin/env python3
"""Hold bioconda/vep-rs to the release contract.

The recipe builds vep, vep-cache-builder and vep-cache-converter from the git
archive the release workflow publishes, with cargo --locked, at the generic CPU
floors. .cargo/config.toml would otherwise compile for x86-64-v3 and neoverse-v1.
Linux links conda-forge OpenSSL, which is built for glibc 2.17. The binaries
install from the release directory CARGO_BUILD_TARGET selects.

usage: check_recipe.py [--root DIR]

Exit 0 when the recipe matches. Exit 1 with one ERROR: [bioconda] line per miss.
No dependency beyond the standard library.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

SOURCE_URL = (
    "https://github.com/natera-open-source/vep-rs/releases/download/"
    "v{{ version }}/vep-rs-{{ version }}.tar.gz"
)
BINS = ("vep", "vep-cache-builder", "vep-cache-converter")
# conda-forge's rust activation exports CARGO_BUILD_TARGET, which moves the binaries to
# target/<triple>/release.
RELEASE_DIR = 'release_dir="target/${CARGO_BUILD_TARGET:+${CARGO_BUILD_TARGET}/}release"'
SHA256_RE = re.compile(r"^  sha256: ([0-9a-f]{64})$", re.MULTILINE)
VERSION_SET_RE = re.compile(r'^\{% set version = "([^"]+)" %\}$', re.MULTILINE)
CARGO_VERSION_RE = re.compile(r'^version = "([^"]+)"$', re.MULTILINE)
MAINTAINER_RE = re.compile(r"^    - ([A-Za-z0-9][A-Za-z0-9-]*)$")


def fail_lines(errors: list[str]) -> None:
    for error in errors:
        print(f"ERROR: [bioconda] {error}", file=sys.stderr)


def cargo_version(text: str) -> str | None:
    """The workspace package version, the first `version = "..."` in Cargo.toml."""
    match = CARGO_VERSION_RE.search(text)
    return match.group(1) if match else None


def section(text: str, heading: str) -> str:
    """The lines of a top-level mapping, from `heading:` to the next top-level key."""
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines) if line == f"{heading}:"), None)
    if start is None:
        return ""
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i] and not lines[i].startswith((" ", "#"))),
        len(lines),
    )
    return "\n".join(lines[start:end])


def recipe_errors(meta: str, build: str, version: str) -> list[str]:
    """What is wrong with the recipe texts, given the workspace version."""
    errors: list[str] = []
    if '{% set name = "vep-rs" %}' not in meta.splitlines():
        errors.append('meta.yaml does not set name to "vep-rs"')
    found = VERSION_SET_RE.search(meta)
    if found is None:
        errors.append("meta.yaml does not set version")
    elif found.group(1) != version:
        errors.append(f"meta.yaml version is {found.group(1)}, Cargo.toml is {version}")
    if "  name: {{ name }}" not in meta.splitlines() or "  version: {{ version }}" not in meta.splitlines():
        errors.append("package name or version is not the jinja variable")
    if f"  url: {SOURCE_URL}" not in meta.splitlines():
        errors.append("source url is not the release git-archive tarball")
    if not SHA256_RE.search(meta):
        errors.append("source sha256 is not 64 hex characters")
    if "  number: 0" not in meta.splitlines():
        errors.append("build number is not 0")
    if 'pin_subpackage(\'vep-rs\', max_pin="x.x")' not in meta:
        errors.append('run_exports does not pin vep-rs with max_pin="x.x"')
    for need in (
        "{{ compiler('rust') }}",
        "{{ compiler('c') }}",
        "{{ stdlib('c') }}",
        "cargo-bundle-licenses",
    ):
        if need not in meta:
            errors.append(f"build requirements lack {need}")
    if "    - pkg-config  # [linux]" not in meta.splitlines():
        errors.append("pkg-config is not a linux build requirement")
    # `host` sits under `requirements`, so it is not a top-level key.
    if "  host:\n    - openssl  # [linux]\n" not in meta:
        errors.append("openssl is not a linux host requirement")
    tests = section(meta, "test")
    for command in (
        "    - vep --version",
        "    - vep-cache-builder --help > /dev/null",
        "    - vep-cache-converter --help > /dev/null",
    ):
        if command not in tests.splitlines():
            errors.append(f"test commands lack {command.strip()}")
    about = section(meta, "about")
    if "  license: Apache-2.0" not in about.splitlines():
        errors.append("license is not Apache-2.0")
    for name in ("LICENSE", "THIRDPARTY.yml"):
        if f"    - {name}" not in about.splitlines():
            errors.append(f"license_file lacks {name}")
    if "  summary: Fast Rust reimplementation of Ensembl VEP." not in about.splitlines():
        errors.append("summary is missing")
    extra = section(meta, "extra")
    for platform in ("    - linux-aarch64", "    - osx-arm64"):
        if platform not in extra.splitlines():
            errors.append(f"additional-platforms lacks {platform.strip()}")
    block = extra.split("recipe-maintainers:", 1)
    maintainers: list[str] = []
    if len(block) == 2:
        for line in block[1].splitlines():
            match = MAINTAINER_RE.match(line)
            if match:
                maintainers.append(match.group(1))
            elif line.strip():
                break
    if not maintainers:
        errors.append("recipe-maintainers is empty")

    if not build.startswith("#!/usr/bin/env bash\n"):
        errors.append("build.sh has no bash shebang")
    if "set -euo pipefail" not in build.splitlines():
        errors.append("build.sh does not set -euo pipefail")
    # The comment names the floors .cargo/config.toml sets; only a setting is a miss.
    code = "\n".join(line for line in build.splitlines() if not line.lstrip().startswith("#"))
    for floor in ("cpu=x86-64-v2", "cpu=generic", "cpu=apple-m1"):
        if floor not in code:
            errors.append(f"build.sh does not set {floor}")
    for measured in ("x86-64-v3", "neoverse"):
        if measured in code:
            errors.append(f"build.sh names the measured CPU floor {measured}")
    if 'echo "ERROR: [bioconda] no CPU floor for $(uname -s)-$(uname -m)"' not in build:
        errors.append("build.sh does not refuse an unknown platform")
    openssl = (
        'if [ "$(uname -s)" = "Linux" ]; then\n'
        '  export OPENSSL_DIR="${PREFIX}"\n'
        "  export OPENSSL_NO_VENDOR=1\n"
        "fi"
    )
    if openssl not in build:
        errors.append("build.sh does not point Linux OpenSSL at the conda prefix")
    if "cargo-bundle-licenses --format yaml --output THIRDPARTY.yml" not in build:
        errors.append("build.sh does not write THIRDPARTY.yml")
    if "cargo build --release --locked" not in build:
        errors.append("build.sh does not cargo build --release --locked")
    if RELEASE_DIR not in build.splitlines():
        errors.append("build.sh does not take release_dir from CARGO_BUILD_TARGET")
    if "target/release" in code:
        errors.append("build.sh names target/release, which CARGO_BUILD_TARGET moves")
    for binary in BINS:
        if f"--bin {binary}" not in build:
            errors.append(f"build.sh does not build {binary}")
        if f'"${{release_dir}}/{binary}"' not in build:
            errors.append(f"build.sh does not install {binary}")
    if "bench_plugins" in build:
        errors.append("build.sh installs bench_plugins")
    if '"${PREFIX}/bin/"' not in build and '"${PREFIX}/bin"' not in build:
        errors.append("build.sh does not install into $PREFIX/bin")
    rustflags = build.find("export RUSTFLAGS=")
    cargo = build.find("cargo build --release --locked")
    licenses = build.find("cargo-bundle-licenses --format yaml --output THIRDPARTY.yml")
    if rustflags < 0 or cargo < 0 or rustflags > cargo:
        errors.append("RUSTFLAGS is not exported before cargo build")
    if licenses < 0 or cargo < 0 or licenses > cargo:
        errors.append("cargo-bundle-licenses does not run before cargo build")
    return errors


def check_root(root: Path) -> list[str]:
    """The recipe errors for a checkout, or why a required file cannot be read."""
    cargo = root / "Cargo.toml"
    meta = root / "bioconda" / "vep-rs" / "meta.yaml"
    build = root / "bioconda" / "vep-rs" / "build.sh"
    missing = [str(path.relative_to(root)) for path in (cargo, meta, build) if not path.is_file()]
    if missing:
        return [f"missing {', '.join(missing)}"]
    version = cargo_version(cargo.read_text(encoding="utf-8"))
    if version is None:
        return ["Cargo.toml has no version"]
    return recipe_errors(
        meta.read_text(encoding="utf-8"),
        build.read_text(encoding="utf-8"),
        version,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    errors = check_root(args.root)
    if errors:
        fail_lines(errors)
        sys.exit(1)


if __name__ == "__main__":
    main()
