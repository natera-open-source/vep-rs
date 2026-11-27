#!/usr/bin/env python3
"""Hold bioconda/vep-rs to the release contract.

The recipe names the newest published release, whose version, source digest and
run requirements are recorded below, and never a version above Cargo.toml's.

The recipe builds vep, vep-cache-builder and vep-cache-converter from the git
archive the release workflow publishes, with cargo --locked, at the generic CPU
floors.
Linux links conda-forge OpenSSL, which is built for glibc 2.17. The binaries
install from the release directory CARGO_BUILD_TARGET selects. build.sh is read
as bash runs it: a comment is no command.

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
VERSION_RE = re.compile(r"(\d+)\.(\d+)\.(\d+)(?:-([0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*))?")
CPU_FLOORS = ("x86-64-v2", "generic", "apple-m1")
# rustc reads `-C target_cpu` and `-C target_feature` as the hyphenated options.
TARGET_CPU_RE = re.compile(r"(target[-_]cpu)=([^\s\"']+)")
TARGET_FEATURE_RE = re.compile(r"(target[-_]feature)=([^\s\"']+)")
CPU_SET_RE = re.compile(r"(?<![\w$-])cpu=[\"']?([^\s;\"']*)")
ESAC_RE = re.compile(r"(?<![\w-])esac(?![\w-])")
CARGO_BUILD_RE = re.compile(r"(?<![\w-])cargo\s+build(?![\w-])")
RUN_KEY_RE = re.compile(r"""\s*["']?run["']?\s*:""")
REQUIREMENTS_KEY_RE = re.compile(r"  [^\s#-]")

# The newest published release: its version, the sha256 its SHA256SUMS states for
# vep-rs-<version>.tar.gz, and what its binaries run. Under crates/ of the 1.0.0
# tarball a process is started only in test code and in bench_plugins, which the
# recipe does not install, so its three binaries start no other program.
PUBLISHED_VERSION = "1.0.0"
PUBLISHED_SHA256 = "a216e0a9449801b18247764a004334f145b013113fca6ac3dbd2e6c2a2794d16"
PUBLISHED_RUN: tuple[str, ...] = ()


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


def listed(block: str, key: str) -> list[str]:
    """The `- item` entries under `  key:` of a section; a comment line is no entry."""
    lines = block.splitlines()
    start = next((i for i, line in enumerate(lines) if line == f"  {key}:"), None)
    if start is None:
        return []
    items: list[str] = []
    for line in lines[start + 1 :]:
        if line.startswith("    - "):
            items.append(line[6:])
        elif line.strip() and not line.lstrip().startswith("#"):
            break
    return items


def precedence(version: str) -> tuple[int, int, int, bool, tuple[tuple[bool, int, str], ...]] | None:
    """The semantic-version sort key of `X.Y.Z` or `X.Y.Z-PRE`, or None for any other shape.

    A pre-release sorts below its own release, a numeric identifier below a named
    one, and a shorter run of identifiers below the longer run it opens.
    """
    match = VERSION_RE.fullmatch(version)
    if match is None:
        return None
    major, minor, patch = (int(part) for part in match.group(1, 2, 3))
    identifiers: list[tuple[bool, int, str]] = []
    for part in match.group(4).split(".") if match.group(4) else ():
        if not part.isdigit():
            identifiers.append((True, 0, part))
        elif part == str(int(part)):
            identifiers.append((False, int(part), ""))
        else:
            return None
    return (major, minor, patch, match.group(4) is None, tuple(identifiers))


def shell_code(script: str) -> str:
    """The script without its comments, line for line.

    An unquoted `#` that starts a word opens a comment; a quote closes on its own line.
    """
    lines: list[str] = []
    for line in script.splitlines():
        quote, escaped, starts_word, end = "", False, True, len(line)
        for index, char in enumerate(line):
            if escaped:
                escaped, starts_word = False, False
            elif char == "\\" and quote != "'":
                escaped = True
            elif quote:
                quote = "" if char == quote else quote
            elif char == "#" and starts_word:
                end = index
                break
            else:
                quote = char if char in "'\"" else ""
                starts_word = char in " \t;&|()"
        lines.append(line[:end].rstrip() if end < len(line) else line)
    return "\n".join(lines)


def recipe_errors(meta: str, build: str, version: str) -> list[str]:
    """What is wrong with the recipe texts, given the workspace version."""
    errors: list[str] = []
    if '{% set name = "vep-rs" %}' not in meta.splitlines():
        errors.append('meta.yaml does not set name to "vep-rs"')
    found = VERSION_SET_RE.search(meta)
    if found is None:
        errors.append("meta.yaml does not set version")
    else:
        recipe = found.group(1)
        if recipe != PUBLISHED_VERSION:
            errors.append(f"meta.yaml version is {recipe}, the published release is {PUBLISHED_VERSION}")
        ours, workspace = precedence(recipe), precedence(version)
        if workspace is None:
            errors.append(f"Cargo.toml version {version} is not MAJOR.MINOR.PATCH")
        elif ours is not None and ours > workspace:
            errors.append(f"meta.yaml version {recipe} is above Cargo.toml's {version}")
    if "  name: {{ name }}" not in meta.splitlines() or "  version: {{ version }}" not in meta.splitlines():
        errors.append("package name or version is not the jinja variable")
    if f"  url: {SOURCE_URL}" not in meta.splitlines():
        errors.append("source url is not the release git-archive tarball")
    digest = SHA256_RE.search(meta)
    if digest is None:
        errors.append("source sha256 is not 64 hex characters")
    elif digest.group(1) != PUBLISHED_SHA256:
        errors.append(f"source sha256 is not the published vep-rs-{PUBLISHED_VERSION}.tar.gz digest")
    if "  number: 0" not in meta.splitlines():
        errors.append("build number is not 0")
    if "{{ pin_subpackage('vep-rs', max_pin=\"x.x\") }}" not in listed(section(meta, "build"), "run_exports"):
        errors.append('run_exports does not pin vep-rs with max_pin="x.x"')
    requirements = section(meta, "requirements")
    tools = listed(requirements, "build")
    for need in (
        "{{ compiler('rust') }}",
        "{{ compiler('c') }}",
        "{{ stdlib('c') }}",
        "cargo-bundle-licenses",
    ):
        if need not in tools:
            errors.append(f"build requirements lack {need}")
    if "    - pkg-config  # [linux]" not in meta.splitlines():
        errors.append("pkg-config is not a linux build requirement")
    # `host` sits under `requirements`, so it is not a top-level key.
    if "  host:\n    - openssl  # [linux]\n" not in meta:
        errors.append("openssl is not a linux host requirement")
    run = listed(requirements, "run")
    for need in PUBLISHED_RUN:
        if need not in run:
            errors.append(f"run requirements lack {need}")
    for item in run:
        if item not in PUBLISHED_RUN:
            errors.append(f"run requirements name {item}, which {PUBLISHED_VERSION} does not run")
    # `listed` reads one `  run:` line of the requirements and the `    - item` lines under it;
    # a run key or entry in any other shape is a requirement the two loops above never see.
    lines = requirements.splitlines()
    keys = [line for line in meta.splitlines() if RUN_KEY_RE.match(line)]
    if keys and (keys != ["  run:"] or "  run:" not in lines):
        errors.append("run requirements are not written as one `  run:` block list under requirements")
    elif keys:
        for line in lines[lines.index("  run:") + 1 :]:
            if REQUIREMENTS_KEY_RE.match(line):
                break
            if line.strip() and not line.lstrip().startswith("#") and not line.startswith("    - "):
                errors.append(f"run requirements hold {line.strip()}, which is no `    - item` line")
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
    for name in ("LICENSE", "NOTICE", "THIRDPARTY.yml"):
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
    # A comment may name any CPU floor; only a setting is a miss.
    code = shell_code(build)
    if "set -euo pipefail" not in code.splitlines():
        errors.append("build.sh does not set -euo pipefail")
    for floor in CPU_FLOORS:
        if f"cpu={floor}" not in code:
            errors.append(f"build.sh does not set cpu={floor}")
    for measured in ("x86-64-v3", "neoverse"):
        if measured in code:
            errors.append(f"build.sh names the measured CPU floor {measured}")
    # An assignment after the case block replaces the floor the block picked.
    esac = ESAC_RE.search(code)
    for assigned in CPU_SET_RE.finditer(code):
        if assigned.group(1) not in CPU_FLOORS or (esac is not None and assigned.start() > esac.start()):
            errors.append(f"build.sh assigns cpu={assigned.group(1)} outside the CPU floors of its case block")
    for flag, value in TARGET_CPU_RE.findall(code):
        if value != "${cpu}":
            errors.append(f"build.sh sets {flag}={value}, not the CPU floor in ${{cpu}}")
    for flag, value in TARGET_FEATURE_RE.findall(code):
        errors.append(f"build.sh sets {flag}={value}, which moves the build off the CPU floor in ${{cpu}}")
    if 'echo "ERROR: [bioconda] no CPU floor for $(uname -s)-$(uname -m)"' not in code:
        errors.append("build.sh does not refuse an unknown platform")
    openssl = (
        'if [ "$(uname -s)" = "Linux" ]; then\n'
        '  export OPENSSL_DIR="${PREFIX}"\n'
        "  export OPENSSL_NO_VENDOR=1\n"
        "fi"
    )
    if openssl not in code:
        errors.append("build.sh does not point Linux OpenSSL at the conda prefix")
    if "cargo-bundle-licenses --format yaml --output THIRDPARTY.yml" not in code:
        errors.append("build.sh does not write THIRDPARTY.yml")
    if "cargo build --release --locked" not in code:
        errors.append("build.sh does not cargo build --release --locked")
    # A command in front of cargo build can start it in an environment without the export.
    for command in code.replace("\\\n", " ").splitlines():
        found = CARGO_BUILD_RE.search(command)
        if found and command[: found.start()].strip():
            errors.append(
                f"build.sh starts cargo build through {command[: found.start()].strip()}, "
                "not as its own command under the exported RUSTFLAGS"
            )
    if RELEASE_DIR not in code.splitlines():
        errors.append("build.sh does not take release_dir from CARGO_BUILD_TARGET")
    if "target/release" in code:
        errors.append("build.sh names target/release, which CARGO_BUILD_TARGET moves")
    for binary in BINS:
        if f"--bin {binary}" not in code:
            errors.append(f"build.sh does not build {binary}")
        if f'"${{release_dir}}/{binary}"' not in code:
            errors.append(f"build.sh does not install {binary}")
    if "bench_plugins" in build:
        errors.append("build.sh installs bench_plugins")
    if '"${PREFIX}/bin/"' not in code and '"${PREFIX}/bin"' not in code:
        errors.append("build.sh does not install into $PREFIX/bin")
    rustflags = code.find("export RUSTFLAGS=")
    cargo = code.find("cargo build --release --locked")
    licenses = code.find("cargo-bundle-licenses --format yaml --output THIRDPARTY.yml")
    if rustflags < 0 or cargo < 0 or rustflags > cargo:
        errors.append("RUSTFLAGS is not exported before cargo build")
    # RUSTFLAGS replaces every rustflags setting, so one line names it: the export of the floor.
    flags = [line.strip() for line in code.splitlines() if "RUSTFLAGS" in line]
    if not (flags and flags[0].startswith("export RUSTFLAGS=") and "-C target-cpu=${cpu}" in flags[0]):
        errors.append("build.sh does not export RUSTFLAGS with -C target-cpu=${cpu}")
    for line in flags[1:]:
        errors.append(f"build.sh changes RUSTFLAGS after its export: {line}")
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
