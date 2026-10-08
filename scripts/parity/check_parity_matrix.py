#!/usr/bin/env python3
"""Check tests/parity/matrix.tsv against the tree.

usage: check_parity_matrix.py [--root DIR]

tests/parity/README.md describes the matrix's columns and how a row closes. The check fails a pull
request when:

- a long flag declared in crates/vep-cli/src/args.rs has no row. A flag `--foo-bar` or `--foo_bar`
  is the row `foo_bar`; aliases are not flags of their own; a flag this project adds beyond
  Ensembl VEP's needs a row too. The flags are read from the `#[arg(...)]` attributes of args.rs,
  and an attribute the reader cannot follow to its field, a `#[command(flatten)]`, a
  `#[command(subcommand)]` or a `rename_all` is an error, never a silent omission;
  crates/vep-cli/tests/parity_matrix.rs asks clap itself and is the arbiter;
- a row in state `landed:<target>` names no `test1`, or no `test2` when its target is
  `implemented` or `documented_subset` (the other targets are proven by `test1` alone), or names
  an `interim` (a landed row has none);
- a named `test1`, `test2` or `interim` entry does not resolve. `test1` and `interim` entries are
  `<path>.rs::<fn>` under `crates/` (no `..` segment) and resolve when the file declares `fn <fn>`
  under a test attribute (`#[test]`, `#[tokio::test]`, `#[rstest]`, any attribute whose name
  carries `test`); `test2` entries are golden corpora `tests/golden/<set>/<corpus>` and resolve
  when the directory exists and its manifest.json has a `covers` list naming the row. Several
  entries in one cell are separated by `;`;
- a row lands at a target other than its own. The check reads only `target` and `state`: a state
  is `open` or `landed:<target>`, the landed target must repeat the row's `target`, so a row whose
  work proved more or less than its target is changed by a target decision, never landed under
  another target; a bare `landed` is an error and an `excluded` row never lands;
- a `removed` row is landed while args.rs still declares its flag;
- a `perl_reference` is empty (`-` is the value for an item with none).

`interim` names a refusal or stop-gap shipped before the row lands, as the test that proves it.

Exit 0 with a one-line summary when clean; exit 1 with one `ERROR: [parity] ...` line per failure.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

COLUMNS = [
    "item",
    "perl_reference",
    "target",
    "state",
    "test1",
    "test2",
    "interim",
]
TARGETS = (
    "implemented",
    "documented_subset",
    "accepted_noop",
    "error",
    "removed",
    "excluded",
)
NEEDS_TEST2 = ("implemented", "documented_subset")
MATRIX = Path("tests/parity/matrix.tsv")
ARGS_RS = Path("crates/vep-cli/src/args.rs")
RUST_TEST = Path("crates/vep-cli/tests/parity_matrix.rs")
TEST_PATH_RE = re.compile(r"^(crates/[A-Za-z0-9_./-]+\.rs)::([A-Za-z_][A-Za-z0-9_]*)$")
CORPUS_RE = re.compile(r"^tests/golden/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
# The clap attributes that declare a flag, in the current and the legacy spelling.
ARG_ATTR_RE = re.compile(r"#\[\s*(?:arg|clap)\s*\(")
# The struct field after an attribute: visibility optional, raw identifiers allowed.
FIELD_RE = re.compile(r"(?:pub(?:\([^)]*\))?\s+)?(?:r#)?([A-Za-z_][A-Za-z0-9_]*)\s*:")
LONG_RE = re.compile(r'^long\s*=\s*r?"([^"]+)"$')
ID_RE = re.compile(r'^id\s*=\s*r?"([^"]+)"$')
# Shapes whose flags this reader cannot derive from args.rs alone.
UNREAD_RE = re.compile(
    r"#\[\s*(?:command|clap)\s*\(\s*(?:flatten|subcommand)\b|\brename_all\b"
)
# A fn's leading keywords, matched backwards from the `fn` token.
FN_KEYWORDS_RE = re.compile(r"(?:\b(?:pub(?:\([^)]*\))?|async|unsafe|const)\s*)+$")
TEST_ATTR_RE = re.compile(r"^\s*(?:\w+::)*\w*test\w*\s*(?:\(|$)")
LINE_COMMENT_RE = re.compile(r"//[^\n]*")
BLOCK_COMMENT_RE = re.compile(r"/\*.*?\*/", re.DOTALL)


def line_of(text: str, index: int) -> int:
    return text.count("\n", 0, index) + 1


def balanced_end(text: str, start: int, open_ch: str, close_ch: str) -> int:
    """Index just past the bracket closing the one before `start`; string and char literals are opaque."""
    depth = 1
    quote = ""
    i = start
    while i < len(text) and depth:
        ch = text[i]
        if quote:
            if ch == "\\":
                i += 1
            elif ch == quote:
                quote = ""
        elif ch in "\"'":
            quote = ch
        elif ch == open_ch:
            depth += 1
        elif ch == close_ch:
            depth -= 1
        i += 1
    return i


def split_top_level(text: str) -> list[str]:
    """Split a clap attribute list on the commas outside quotes, brackets and parentheses."""
    parts: list[str] = []
    depth = 0
    quote = ""
    current: list[str] = []
    i = 0
    while i < len(text):
        ch = text[i]
        if quote:
            current.append(ch)
            if ch == "\\" and i + 1 < len(text):
                current.append(text[i + 1])
                i += 1
            elif ch == quote:
                quote = ""
        elif ch in "\"'":
            quote = ch
            current.append(ch)
        elif ch in "([{":
            depth += 1
            current.append(ch)
        elif ch in ")]}":
            depth -= 1
            current.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
        i += 1
    tail = "".join(current).strip()
    if tail:
        parts.append(tail)
    return parts


def skip_trivia(text: str, i: int) -> int:
    """Index past the whitespace, line comments and further attributes that may precede a field."""
    while True:
        while i < len(text) and text[i].isspace():
            i += 1
        if text.startswith("//", i):
            end = text.find("\n", i)
            i = len(text) if end < 0 else end + 1
        elif text.startswith("#[", i):
            i = balanced_end(text, i + 2, "[", "]")
        else:
            return i


def scan_args(source: str) -> tuple[list[str], list[str]]:
    """(long names, errors) read from the `#[arg(...)]` attributes of args.rs.

    A bare `long` is the field name in kebab-case, or the attribute's `id` when one is given;
    `long = "name"` is the name. An attribute that is not followed by a field, and every shape
    whose flags come from elsewhere (`flatten`, `subcommand`, `rename_all`), is an error.
    """
    flags: list[str] = []
    errors: list[str] = []
    for unread in UNREAD_RE.finditer(source):
        errors.append(
            f"{ARGS_RS}:{line_of(source, unread.start())}: {unread.group(0).strip()} is not read by this check; "
            f"{RUST_TEST} covers the flags it brings"
        )
    pos = 0
    while (attr := ARG_ATTR_RE.search(source, pos)) is not None:
        line = line_of(source, attr.start())
        end = balanced_end(source, attr.end(), "(", ")")
        parts = split_top_level(source[attr.end() : end - 1])
        after = skip_trivia(source, end)
        if not source.startswith("]", after):
            errors.append(f"{ARGS_RS}:{line}: attribute is not closed by `]`")
            pos = end
            continue
        field = FIELD_RE.match(source, skip_trivia(source, after + 1))
        if field is None:
            errors.append(
                f"{ARGS_RS}:{line}: #[arg(...)] is not followed by a field this check can read"
            )
            pos = after + 1
            continue
        ident = next((m.group(1) for part in parts if (m := ID_RE.match(part))), None)
        for part in parts:
            if part == "long":
                flags.append(
                    ident if ident is not None else field.group(1).replace("_", "-")
                )
            elif (named := LONG_RE.match(part)) is not None:
                flags.append(named.group(1))
        pos = field.end()
    return flags, errors


def read_matrix(root: Path) -> tuple[list[dict], list[str]]:
    errors: list[str] = []
    with (root / MATRIX).open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != COLUMNS:
            return [], [f"{MATRIX}: header is {header}, expected {COLUMNS}"]
        rows = []
        for lineno, cells in enumerate(reader, start=2):
            if len(cells) != len(COLUMNS):
                errors.append(
                    f"{MATRIX}:{lineno}: {len(cells)} columns, expected {len(COLUMNS)}"
                )
                continue
            rows.append(dict(zip(COLUMNS, cells), lineno=lineno))
    seen: dict[str, int] = {}
    for row in rows:
        item, state = row["item"], row["state"]
        where = f"{MATRIX}:{row['lineno']}: {item}"
        if not item:
            errors.append(f"{MATRIX}:{row['lineno']}: empty item")
        elif item in seen:
            errors.append(
                f"{MATRIX}:{row['lineno']}: item {item!r} repeats line {seen[item]}"
            )
        else:
            seen[item] = row["lineno"]
        if not row["perl_reference"]:
            errors.append(f"{where}: empty perl_reference (`-` names none)")
        if row["target"] not in TARGETS:
            errors.append(
                f"{where}: target {row['target']!r} is not one of {', '.join(TARGETS)}"
            )
        if state == "landed":
            errors.append(
                f"{where}: state 'landed' must name the target it landed at: landed:<target>"
            )
        elif state != "open" and not state.startswith("landed:"):
            errors.append(f"{where}: state {state!r} is not open or landed:<target>")
    return rows, errors


def entries(cell: str) -> list[str]:
    return [entry.strip() for entry in cell.split(";") if entry.strip()]


def attribute_spans(text: str) -> dict[int, tuple[int, str]]:
    """Every `#[...]` attribute, keyed by its end index: (start index, attribute path)."""
    spans: dict[int, tuple[int, str]] = {}
    i = 0
    while (start := text.find("#[", i)) >= 0:
        end = balanced_end(text, start + 2, "[", "]")
        spans[end] = (start, text[start + 2 : end - 1])
        i = end
    return spans


def declares_test_fn(text: str, fn: str) -> bool:
    """Whether the source declares `fn <fn>` under a test attribute, comments stripped first."""
    code = LINE_COMMENT_RE.sub("", BLOCK_COMMENT_RE.sub("", text))
    spans = attribute_spans(code)
    for decl in re.finditer(r"\bfn\s+" + re.escape(fn) + r"\b", code):
        head = FN_KEYWORDS_RE.sub("", code[: decl.start()].rstrip()).rstrip()
        while head.endswith("]") and len(head) in spans:
            start, attribute = spans[len(head)]
            if TEST_ATTR_RE.match(attribute):
                return True
            head = code[:start].rstrip()
    return False


def resolve_test(root: Path, entry: str) -> str | None:
    """None when `<path>.rs::<fn>` names a file under crates/ declaring the test function, else the reason."""
    match = TEST_PATH_RE.match(entry)
    if not match:
        return "is not <path>.rs::<fn> under crates/"
    path, fn = match.groups()
    if ".." in path.split("/"):
        return "has a `..` segment"
    file = root / path
    if not file.is_file():
        return f"{path} does not exist"
    if not declares_test_fn(file.read_text(encoding="utf-8", errors="replace"), fn):
        return f"{path} declares no test fn {fn}"
    return None


def resolve_corpus(root: Path, entry: str, item: str) -> str | None:
    """None when the corpus directory exists and its manifest's `covers` names the item, else the reason."""
    if not CORPUS_RE.match(entry):
        return "is not tests/golden/<set>/<corpus>"
    corpus = root / entry
    if not corpus.is_dir():
        return f"{entry} does not exist"
    manifest = corpus / "manifest.json"
    if not manifest.is_file():
        return f"{entry} has no manifest.json"
    try:
        covers = json.loads(manifest.read_text(encoding="utf-8")).get("covers")
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        return f"{entry}/manifest.json is not JSON ({exc})"
    if not isinstance(covers, list):
        return f"{entry}/manifest.json has no covers list"
    if item not in covers:
        return f"{entry}/manifest.json covers does not name {item}"
    return None


def check(root: Path) -> tuple[list[str], str]:
    """Run every rule; return (errors, one-line summary)."""
    args_rs = root / ARGS_RS
    if not (root / MATRIX).is_file():
        return [f"{MATRIX} is missing"], ""
    if not args_rs.is_file():
        return [f"{ARGS_RS} is missing"], ""
    rows, errors = read_matrix(root)
    if errors:
        return errors, ""
    items = {row["item"] for row in rows}

    flags, errors = scan_args(args_rs.read_text(encoding="utf-8"))
    declared = {flag.replace("-", "_") for flag in flags}
    for flag in flags:
        name = flag.replace("-", "_")
        if name not in items:
            errors.append(f"{ARGS_RS}: flag --{flag} has no row {name!r} in {MATRIX}")

    landed = 0
    interim = 0
    for row in rows:
        item, target, state = row["item"], row["target"], row["state"]
        where = f"{MATRIX}:{row['lineno']}: {item}"
        is_landed = state != "open"
        if is_landed:
            landed += 1
            landed_at = state.split(":", 1)[1]
            if target == "excluded":
                errors.append(f"{where}: an excluded row never lands")
            elif landed_at != target:
                errors.append(
                    f"{where}: landed at {landed_at!r}, its target is {target!r}"
                )
            if target == "removed" and item in declared:
                errors.append(
                    f"{where}: landed as removed while {ARGS_RS} still declares the flag"
                )
            if not entries(row["test1"]):
                errors.append(f"{where}: landed with no test1")
            if target in NEEDS_TEST2 and not entries(row["test2"]):
                errors.append(f"{where}: landed {target} with no test2")
            if entries(row["interim"]):
                errors.append(f"{where}: a landed row has no interim")
        for entry in entries(row["test1"]):
            reason = resolve_test(root, entry)
            if reason:
                errors.append(f"{where}: test1 {entry} {reason}")
        for entry in entries(row["test2"]):
            reason = resolve_corpus(root, entry, item)
            if reason:
                errors.append(f"{where}: test2 {entry} {reason}")
        for entry in entries(row["interim"]):
            interim += 1
            reason = resolve_test(root, entry)
            if reason:
                errors.append(f"{where}: interim {entry} {reason}")

    summary = f"parity: {len(rows)} rows, {len(flags)} flags covered, {landed} landed, {interim} interim entries resolve"
    return errors, summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--root", default=".", help="repository root (default: the working directory)"
    )
    args = parser.parse_args(argv)
    errors, summary = check(Path(args.root))
    for error in errors:
        print(f"ERROR: [parity] {error}", file=sys.stderr)
    if errors:
        return 1
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
