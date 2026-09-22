"""Build manuscript/data/code_metrics.csv.

Per-crate lines of code and `#[test]` counts for the released crates, plus a
TOTAL row. No Cargo build is required: the script walks `crates/<name>/src`
and `crates/<name>/tests`.

LOC is "non-blank, non-comment Rust source lines", counted with a stripped-down
replication of tokei's logic so the script has no external dependencies. Test
count is "occurrences of `#[test]` in crate sources", which under-counts test
cases declared via macro expansion.

Every directory under crates/ is either in CRATE_PURPOSES (counted) or named
by `--exclude` (skipped). An unclassified directory raises, so a newly added
crate is classified deliberately instead of silently moving the totals.
"""

from __future__ import annotations

import csv
import argparse
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CRATES_DIR = REPO_ROOT / "crates"

CRATE_PURPOSES = {
    "vep-core": "Canonical types, SO terms, prediction decoder",
    "vep-effects": "Consequence engine + 6 SV submodules",
    "vep-fasta": "Indexed FASTA wrapper",
    "vep-io": "Input parsers + output formatters",
    "vep-builtins": "Built-in annotation plugins on shared tabix engine",
    "vep-plugin": "C-FFI dynamic library loader",
    "vep-cli": "Binary, runner, transcript index",
    "vep-cache-converter": "JSON cache validator + plugin data optimizer",
    "vep-cache-builder": "Native GFF3+FASTA+VCF cache builder",
}


_BLOCK_COMMENT_OPEN = re.compile(r"/\*")
_BLOCK_COMMENT_CLOSE = re.compile(r"\*/")
_LINE_COMMENT = re.compile(r"^\s*//")


def _count_loc(path: Path) -> int:
    """Approximate non-blank, non-comment Rust LOC of `path`; 0 when it is not a file.

    Block comments are tracked across lines; line comments and blank lines are
    skipped; strings containing `/*` or `*/` are not handled specially. Line comments
    are tested before the block-comment open: a `//` comment whose text contains `/*`
    (a doc note citing the Perl codon output ``Y/*``) is not a block-comment opener,
    and testing the block case first would latch the counter into a block comment
    and discard every remaining line of the file.
    """
    if not path.is_file():
        return 0
    in_block = False
    loc = 0
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if in_block:
            if _BLOCK_COMMENT_CLOSE.search(line):
                in_block = False
            continue
        if _LINE_COMMENT.match(line):
            continue
        if _BLOCK_COMMENT_OPEN.search(line) and not _BLOCK_COMMENT_CLOSE.search(line):
            in_block = True
            continue
        loc += 1
    return loc


def _count_tests(path: Path) -> int:
    if not path.is_file():
        return 0
    return path.read_text(errors="replace").count("#[test]")


def _walk_rust_sources(crate_root: Path):
    src = crate_root / "src"
    if not src.is_dir():
        return
    yield from src.rglob("*.rs")
    tests = crate_root / "tests"
    if tests.is_dir():
        yield from tests.rglob("*.rs")


def main() -> None:
    ap = argparse.ArgumentParser(description="Count Rust lines and tests per released crate into manuscript/data/code_metrics.csv")
    ap.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="CRATE",
        help="a workspace member to leave out of the table, skipped without being classified (repeatable)",
    )
    args = ap.parse_args()
    excluded = set(args.exclude)
    out_dir = Path(__file__).resolve().parent.parent / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "code_metrics.csv"

    # An absent workspace must fail here: otherwise the loop below iterates an
    # empty list, writes a TOTAL-only CSV of zeros, and exits 0.
    if not CRATES_DIR.is_dir():
        raise SystemExit(
            f"no crates/ directory at {CRATES_DIR}: run this from a checkout that "
            "contains crates/"
        )

    rows: list[dict[str, str]] = []
    total_loc = 0
    total_tests = 0
    for crate_path in sorted(CRATES_DIR.iterdir()):
        if not crate_path.is_dir():
            continue
        name = crate_path.name
        if name in excluded:
            continue
        if name not in CRATE_PURPOSES:
            raise SystemExit(
                f"unclassified workspace crate {name!r}: add it to CRATE_PURPOSES "
                f"(counted) or pass --exclude {name} (skipped)"
            )
        loc = sum(_count_loc(p) for p in _walk_rust_sources(crate_path))
        tests = sum(_count_tests(p) for p in _walk_rust_sources(crate_path))
        purpose = CRATE_PURPOSES.get(name, "")
        rows.append(
            {
                "crate": name,
                "loc": str(loc),
                "tests": str(tests),
                "purpose": purpose,
            }
        )
        total_loc += loc
        total_tests += tests

    rows.append(
        {
            "crate": "TOTAL",
            "loc": str(total_loc),
            "tests": str(total_tests),
            "purpose": "workspace aggregate",
        }
    )

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["crate", "loc", "tests", "purpose"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} rows -> {out_path}")
    print(f"total LOC: {total_loc:,}  tests: {total_tests:,}")


if __name__ == "__main__":
    main()
