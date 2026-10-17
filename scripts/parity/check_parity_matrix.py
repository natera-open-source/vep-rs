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
  carries `test`); a `test1` entry may instead be `scripts/<dir>/test_<name>.py::<fn>`, a pytest
  function, resolving when the module exists under `scripts/` (no `..` segment) and declares
  `def <fn>(` at module level or inside a class (comments and strings are not declarations);
  `test2` entries are golden corpora `tests/golden/<set>/<corpus>` and resolve when the directory
  exists and its manifest.json names the row in its `covers` list or in the `covers` of one of
  its `flag_sets` entries. Several entries in one cell are separated by `;`;
- a process row names an Ensembl VEP module. A row with a scripts `test1` entry is a process row:
  it registers a release procedure rather than a behaviour of Ensembl VEP's, so its
  `perl_reference` must be `-`; its `test2` resolves as any other row's;
- a `covers` list reached through a `test2` cell names something that is not a row of the matrix
  (the manifest's own list and every flag set's are read);
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
import ast
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
SCRIPT_TEST_RE = re.compile(
    r"^(scripts/[A-Za-z0-9_./-]+/test_[A-Za-z0-9_]+\.py)::(test_[A-Za-z0-9_]*)$"
)
CRATE_SHAPE = "<path>.rs::<fn> under crates/"
SCRIPT_SHAPE = "scripts/<dir>/test_<name>.py::<fn>"
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
IDENT_CHARS = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")


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


def strip_comments(text: str) -> str:
    """The source without its line and (nested) block comments, every string, raw-string and
    char literal reduced to its delimiters around blanks of the body's length (newlines kept):
    a `/*` inside a string (`"U/*"`, an amino-acid change) or a `//` inside one opens no
    comment, and a `fn` inside a comment or a literal declares nothing."""
    out: list[str] = []

    def blanked(literal: str, open_len: int, close_len: int) -> str:
        end = len(literal) - close_len
        return literal[:open_len] + re.sub(r"[^\n]", " ", literal[open_len:end]) + literal[end:]

    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            end = text.find("\n", i)
            i = n if end < 0 else end
        elif text.startswith("/*", i):
            depth, i = 1, i + 2
            while i < n and depth:
                if text.startswith("/*", i):
                    depth, i = depth + 1, i + 2
                elif text.startswith("*/", i):
                    depth, i = depth - 1, i + 2
                else:
                    i += 1
        elif c == "r" and (hashes := _raw_string_hashes(text, i)) is not None:
            close = '"' + "#" * hashes
            end = text.find(close, i + 1 + hashes + 1)
            end = n if end < 0 else end + len(close)
            out.append(blanked(text[i:end], 1 + hashes + 1, len(close) if text.endswith(close, 0, end) else 0))
            i = end
        elif c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            out.append(blanked(text[i : j + 1], 1, 1 if j < n else 0))
            i = j + 1
        elif c == "'" and i + 2 < n and (text[i + 1] == "\\" or text[i + 2] == "'"):
            # A char literal (`'x'`, `'\''`, `'\u{1F600}'`); a lifetime (`'a`) falls through.
            j = text.find("'", i + 3 if text[i + 1] == "\\" else i + 2)
            j = n - 1 if j < 0 else j
            out.append(blanked(text[i : j + 1], 1, 1))
            i = j + 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def _raw_string_hashes(text: str, i: int) -> int | None:
    """The hash count of a raw string literal starting at `text[i] == "r"` (`r"…"`, `r#"…"#`,
    `br"…"`), or None when the `r` is part of an identifier or no quote follows."""
    before = text[i - 1] if i > 0 else ""
    if before == "b":
        before = text[i - 2] if i > 1 else ""
    if before in IDENT_CHARS:
        return None
    j = i + 1
    while j < len(text) and text[j] == "#":
        j += 1
    return j - i - 1 if j < len(text) and text[j] == '"' else None


def declares_test_fn(text: str, fn: str) -> bool:
    """Whether the source declares `fn <fn>` under a test attribute, comments stripped first."""
    code = strip_comments(text)
    spans = attribute_spans(code)
    for decl in re.finditer(r"\bfn\s+" + re.escape(fn) + r"\b", code):
        head = FN_KEYWORDS_RE.sub("", code[: decl.start()].rstrip()).rstrip()
        while head.endswith("]") and len(head) in spans:
            start, attribute = spans[len(head)]
            if TEST_ATTR_RE.match(attribute):
                return True
            head = code[:start].rstrip()
    return False


def declares_pytest_fn(source: str, fn: str) -> tuple[bool, str | None]:
    """(declared, reason): whether the module declares `def <fn>(` at module level or inside a
    class; the source is parsed, so a `def` inside a comment or a string declares nothing, and a
    module that does not parse declares nothing with the parse error as the reason."""
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        return False, f"does not parse ({exc.msg}, line {exc.lineno})"
    functions = (ast.FunctionDef, ast.AsyncFunctionDef)
    for node in tree.body:
        if isinstance(node, functions) and node.name == fn:
            return True, None
        if isinstance(node, ast.ClassDef) and any(
            isinstance(inner, functions) and inner.name == fn for inner in node.body
        ):
            return True, None
    return False, None


def resolve_test(root: Path, entry: str, scripts_allowed: bool = False) -> str | None:
    """None when the entry names a test the tree declares, else the reason. The crate shape,
    `<path>.rs::<fn>` under crates/, resolves to a function under a test attribute; the scripts
    shape, `scripts/<dir>/test_<name>.py::<fn>`, accepted for `test1` only, resolves to a pytest
    function of the module."""
    script = SCRIPT_TEST_RE.match(entry) if scripts_allowed else None
    match = TEST_PATH_RE.match(entry)
    if not match and not script:
        shapes = f"{CRATE_SHAPE} or {SCRIPT_SHAPE}" if scripts_allowed else CRATE_SHAPE
        return f"is not {shapes}"
    path, fn = (script or match).groups()
    if ".." in path.split("/"):
        return "has a `..` segment"
    file = root / path
    if not file.is_file():
        return f"{path} does not exist"
    source = file.read_text(encoding="utf-8", errors="replace")
    if script:
        declared, reason = declares_pytest_fn(source, fn)
        if reason:
            return f"{path} {reason}"
        if not declared:
            return f"{path} declares no def {fn}"
    elif not declares_test_fn(source, fn):
        return f"{path} declares no test fn {fn}"
    return None


def is_process_row(row: dict) -> bool:
    """A row with a scripts `test1` entry registers a release procedure, not an Ensembl VEP behaviour."""
    return any(SCRIPT_TEST_RE.match(entry) for entry in entries(row["test1"]))


def manifest_covers(manifest: dict) -> tuple[list[tuple[str | None, list]], str | None]:
    """Every `covers` list of a manifest, (None, the manifest's own) then (set name, the set's)
    for each `flag_sets` entry, or the reason one of them cannot be read."""
    covers = manifest.get("covers")
    if not isinstance(covers, list):
        return [], "has no covers list"
    lists: list[tuple[str | None, list]] = [(None, covers)]
    flag_sets = manifest.get("flag_sets", [])
    if not isinstance(flag_sets, list):
        return [], "has a flag_sets that is not a list"
    for flag_set in flag_sets:
        name = flag_set.get("name") if isinstance(flag_set, dict) else None
        if not isinstance(name, str):
            return [], "has a flag_sets entry without a name"
        set_covers = flag_set.get("covers", [])
        if not isinstance(set_covers, list):
            return [], f"flag set {name} has a covers that is not a list"
        lists.append((name, set_covers))
    return lists, None


def read_corpus_manifest(root: Path, entry: str) -> tuple[dict | None, str | None]:
    """The corpus's manifest, or the reason it cannot be read."""
    if not CORPUS_RE.match(entry):
        return None, "is not tests/golden/<set>/<corpus>"
    corpus = root / entry
    if not corpus.is_dir():
        return None, f"{entry} does not exist"
    manifest = corpus / "manifest.json"
    if not manifest.is_file():
        return None, f"{entry} has no manifest.json"
    try:
        parsed = json.loads(manifest.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        return None, f"{entry}/manifest.json is not JSON ({exc})"
    if not isinstance(parsed, dict):
        return None, f"{entry}/manifest.json is not a JSON object"
    return parsed, None


def resolve_corpus(root: Path, entry: str, item: str, manifests: dict) -> str | None:
    """None when the corpus directory exists and its manifest names the item in its `covers` or in
    a flag set's, else the reason. `manifests` caches the parsed manifests by entry."""
    if entry not in manifests:
        manifests[entry] = read_corpus_manifest(root, entry)
    manifest, reason = manifests[entry]
    if reason:
        return reason
    lists, reason = manifest_covers(manifest)
    if reason:
        return f"{entry}/manifest.json {reason}"
    if not any(item in covers for _name, covers in lists):
        return f"{entry}/manifest.json covers does not name {item}"
    return None


def covers_outside_matrix(entry: str, manifest: dict, items: set[str]) -> list[str]:
    """One error per `covers` name of the manifest (its own list and every flag set's) that is
    not a row; a value that is not a string names no row."""
    errors: list[str] = []
    lists, _reason = manifest_covers(manifest)
    for name, covers in lists:
        where = f"{entry}/manifest.json" + (f" flag set {name}" if name else "")
        errors.extend(
            f"{where} covers names {covered}, which is not a row"
            for covered in covers
            if not isinstance(covered, str) or covered not in items
        )
    return errors


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
    manifests: dict[str, tuple[dict | None, str | None]] = {}
    for row in rows:
        item, target, state = row["item"], row["target"], row["state"]
        where = f"{MATRIX}:{row['lineno']}: {item}"
        if is_process_row(row) and row["perl_reference"] != "-":
            errors.append(
                f"{where}: a process row names no Ensembl VEP module; perl_reference must be -"
            )
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
            reason = resolve_test(root, entry, scripts_allowed=True)
            if reason:
                errors.append(f"{where}: test1 {entry} {reason}")
        for entry in entries(row["test2"]):
            reason = resolve_corpus(root, entry, item, manifests)
            if reason:
                errors.append(f"{where}: test2 {entry} {reason}")
        for entry in entries(row["interim"]):
            interim += 1
            reason = resolve_test(root, entry)
            if reason:
                errors.append(f"{where}: interim {entry} {reason}")
    for entry, (manifest, reason) in manifests.items():
        if manifest is not None and reason is None:
            errors.extend(covers_outside_matrix(entry, manifest, items))

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
