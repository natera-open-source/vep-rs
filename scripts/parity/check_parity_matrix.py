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
- a landed row whose target is `documented_subset` names no `test1` entry that reads a public text,
  the test that the documented exclusion is stated. A crate or tools entry reads one when the body of
  its test function holds a string literal that is a Markdown file's path, with no white space in
  it and the file's name at its end (`"cache-setup.md"`, the page or the README whose sentences
  it asserts; a message that names a page is no path), or asks clap for a help text (`get_help`,
  `get_long_help`, `render_help`, `render_long_help`), or calls, by its bare name, a function of
  its own file whose body does, through any number of such calls. The file a test stands in and
  the name it carries show nothing of what it reads: a behaviour test of
  crates/vep-cli/tests/documented_subsets.rs states no exclusion;
- a named `test1`, `test2` or `interim` entry does not resolve. `test1` and `interim` entries are
  `<path>.rs::<fn>` under `crates/` (no `..` segment) and resolve when the file declares `fn <fn>`
  under a test attribute (`#[test]`, `#[tokio::test]`, `#[rstest]`, any attribute whose name
  carries `test`); a `test1` entry may instead be `scripts/<dir>/test_<name>.py::<fn>`, a pytest
  function, resolving when the module exists under `scripts/` (no `..` segment) and declares
  `def <fn>(` at module level or inside a class (comments and strings are not declarations), or
  `scripts/[<dir>/]test_<name>.sh::<check>`, a named check of a shell suite, resolving when the
  suite exists under `scripts/` (no `..` segment) and calls a check function, as a command,
  with `<check>`, whole and quoted, as its first argument: `ok`, the function that prints a
  passed check's line, or a function of the suite's own that hands its first argument on to a
  check function (a comment, a here-document and a quoted text declare nothing, and a name
  holds no `$`, backtick or backslash, so a name the suite computes is never an entry), or
  `tools/<tool>/<path>.rs::<fn>`, a test of a tool's own Cargo root, resolving as a crate entry
  does;
  `test2` entries are golden corpora `tests/golden/<set>/<corpus>` and resolve when the directory
  exists and its manifest.json names the row in its `covers` list or in the `covers` of one of
  its `flag_sets` entries. Several entries in one cell are separated by `;`;
- a process row names an Ensembl VEP module. A row with a scripts or tools `test1` entry is a process row:
  it registers a release procedure rather than a behaviour of Ensembl VEP's, so its
  `perl_reference` must be `-`; its `test2` resolves as any other row's;
- a `covers` list reached through a `test2` cell names something that is not a row of the matrix
  (the manifest's own list and every flag set's are read);
- a row lands at a target other than its own. The check reads only `target` and `state`: a state
  is `open` or `landed:<target>`, the landed target must repeat the row's `target`, so a row whose
  work proved more or less than its target is changed by a target decision, never landed under
  another target; a bare `landed` is an error and an `excluded` row never lands;
- a `removed` row is landed while args.rs still declares its flag;
- a `perl_reference` is empty (`-` is the value for an item with none);
- a `perl_reference` is in none of the forms tests/parity/README.md gives. A cell is `-`; a line
  of tests/parity/reference_files.txt, which lists the scripts of ensembl-vep by name, its modules
  by their path below `modules/Bio/EnsEMBL/VEP/`, the files under `modules/Bio/EnsEMBL/Variation/`
  of ensembl-variation and under `modules/Bio/EnsEMBL/IO/` of ensembl-io by their path below
  `modules/Bio/EnsEMBL/`, and one BioPerl package; or a plugin's module, `<Name>.pm` for a plugin
  `<Name>` of the support table and `LoF.pm` for LoFTEE (a plugin row's own `<item>.pm` is held
  by the table rule below). A file may carry after `::` the consuming sub or the parameter the
  module sets, one identifier. A module named by its last components is therefore refused,
  except where those components are a line of their own (`Parser/VCF.pm`, also the end of
  `Haplo/Parser/VCF.pm`). The list itself is read whole: a missing list, a line that is not one
  name and a line that does not sort after the line before it, in byte order, stop the check;
- the plugin support table of docs/plugins.md disagrees with the matrix. The table is the six-column
  table under the `## Support table` heading, one row per Ensembl VEP plugin, its first cell the
  plugin name (linked or bare) and its fifth the status. Every table row names a matrix row (its
  `item`), and the status fixes the row's `target`: `supported, verified on real data` and
  `supported on the user's licensed file, verified on synthetic fixtures` are `implemented`, `not ported`
  is `error`; a supported status also needs the row landed, since the table's "verified" is a claim
  about shipped code, and `not ported` fits any open row, whatever its target. Every plugin row of
  the matrix, one whose `perl_reference` is `<item>.pm`, has a table row; a plugin row named
  otherwise (LoFTEE's `LoF.pm`) is checked through its table row alone. The rows keep alphabetical
  order of the plugin name, case ignored. A row without six cells, a status outside the three
  phrases or an empty or repeated name is malformed.

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
PLUGINS_PAGE = Path("docs/plugins.md")
REFERENCE_FILES = Path("tests/parity/reference_files.txt")
# The sub or the parameter a reference names after `::`.
REFERENCE_SUB_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
# The plugins whose module is not named after the plugin.
PLUGIN_MODULES_NAMED_OTHERWISE = {"LoFTEE": "LoF.pm"}
SUPPORT_TABLE_HEADING = "## Support table"
SUPPORT_TABLE_COLUMNS = 6
STATUS_TARGETS = {
    "supported, verified on real data": "implemented",
    "supported on the user's licensed file, verified on synthetic fixtures": "implemented",
    "not ported": "error",
}
# A table cell that is one Markdown link: the plugin name is its text.
LINK_CELL_RE = re.compile(r"^\[([^\]]+)\]\([^)]*\)$")
SEPARATOR_ROW_RE = re.compile(r"^\|[\s:|-]*\|$")
TEST_PATH_RE = re.compile(r"^(crates/[A-Za-z0-9_./-]+\.rs)::([A-Za-z_][A-Za-z0-9_]*)$")
SCRIPT_TEST_RE = re.compile(
    r"^(scripts/[A-Za-z0-9_./-]+/test_[A-Za-z0-9_]+\.py)::(test_[A-Za-z0-9_]*)$"
)
SHELL_TEST_RE = re.compile(
    r"^(scripts/(?:[A-Za-z0-9_./-]+/)?test_[A-Za-z0-9_]+\.sh)::([^$`\\\t]+)$"
)
TOOL_TEST_RE = re.compile(
    r"^(tools/[A-Za-z0-9_.-]+/[A-Za-z0-9_./-]+\.rs)::([A-Za-z_][A-Za-z0-9_]*)$"
)
CRATE_SHAPE = "<path>.rs::<fn> under crates/"
SCRIPT_SHAPE = "scripts/<dir>/test_<name>.py::<fn>"
SHELL_SHAPE = "scripts/[<dir>/]test_<name>.sh::<check>"
TOOL_SHAPE = "tools/<tool>/<path>.rs::<fn>"
# A shell function's definition line: `name() {` or `function name {`, with what precedes it.
SHELL_FN_RE = re.compile(
    r"^([ \t]*)(?:function[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*(?:\(\))?|([A-Za-z_][A-Za-z0-9_]*)[ \t]*\(\))[ \t]*\{"
)
# Where a command starts on a shell line: the line's start, a separator, the closing bracket of
# a case arm, or a word that a command follows; then the assignments (`NAME=value`) that stand
# before the command's name.
SHELL_COMMAND_START = (
    r"(?:^|[;&|({!)]|\b(?:then|else|elif|do|if|while|until)\b)[ \t]*"
    r"(?:[A-Za-z_][A-Za-z0-9_]*=(?:[^\s\"'\\;&|()<>]|\\.|\"(?:[^\"\\]|\\.)*\"|'[^']*')*[ \t]+)*"
)
# The opening of a here-document: `<<WORD`, `<<-WORD`, the word in quotes, after a backslash
# or bare; never `<<<`.
SHELL_HEREDOC_RE = re.compile(
    r"(?<!<)<<(?!<)(-?)[ \t]*(?:(['\"])([^'\"]+)\2|\\?([A-Za-z0-9_][A-Za-z0-9_.-]*))"
)
# A variable a shell function sets to its first argument: `name="$1"`, `name=$1`, `name=${1}`.
SHELL_FIRST_ARGUMENT_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)=\"?\$\{?1\}?(?![0-9])")
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
# The target whose landed rows name a test that reads a public text.
STATED_TARGET = "documented_subset"
# A Rust function's declaration, and the tokens that end its signature or nest inside it.
FN_DECL_RE = re.compile(r"\bfn\s+([A-Za-z_][A-Za-z0-9_]*)\b")
SIGNATURE_TOKEN_RE = re.compile(r"[()\[\]{;]")
BRACE_RE = re.compile(r"[{}]")
# A string literal that is a Markdown file's path, with no white space in it and the file's name
# at its end: `"cache-setup.md"`, `"{}/../../docs/{path}.md"`.
PAGE_NAME_RE = re.compile(r'"[^\s"]*[A-Za-z0-9_)}-]\.md"')
# The clap calls that give an argument's or a command's help text.
HELP_TEXT_RE = re.compile(r"\b(?:get_help|get_long_help|render_help|render_long_help)\b")
# A function called by its bare name: no `.` and no `::` before it.
BARE_CALL_RE = re.compile(r"(?<![.:A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)\s*\(")


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
    # A byte that is not UTF-8 is replaced and a cell past csv's default limit is read, so neither
    # stops the checker before it judges every row it can.
    csv.field_size_limit(1 << 24)
    with (root / MATRIX).open(newline="", encoding="utf-8", errors="replace") as fh:
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


def table_cells(line: str) -> list[str]:
    """The cells of a Markdown table row, split on the pipes not escaped as `\\|`, each stripped."""
    inner = line.strip().removeprefix("|")
    if inner.endswith("|") and not inner.endswith("\\|"):
        inner = inner[:-1]
    return [cell.replace("\\|", "|").strip() for cell in re.split(r"(?<!\\)\|", inner)]


def read_support_table(root: Path) -> tuple[list[dict], list[str]]:
    """The rows of the plugin support table, each `{"name", "status", "lineno"}`, or the reasons the
    table cannot be read: the page or the heading missing, no table under the heading, no separator
    row, a header, separator or row without six cells, a status outside STATUS_TARGETS, an empty or
    a repeated name."""
    page = root / PLUGINS_PAGE
    if not page.is_file():
        return [], [f"{PLUGINS_PAGE} is missing"]
    lines = page.read_text(encoding="utf-8").splitlines()
    heading = next(
        (i for i, line in enumerate(lines) if line.strip() == SUPPORT_TABLE_HEADING), None
    )
    if heading is None:
        return [], [f"{PLUGINS_PAGE}: no {SUPPORT_TABLE_HEADING!r} heading"]
    header = None
    for index in range(heading + 1, len(lines)):
        if lines[index].startswith("#"):
            break
        if lines[index].lstrip().startswith("|"):
            header = index
            break
    if header is None:
        return [], [f"{PLUGINS_PAGE}: no table under {SUPPORT_TABLE_HEADING!r}"]
    malformed = f"{PLUGINS_PAGE}: support table malformed at line"
    header_cells = len(table_cells(lines[header]))
    if header_cells != SUPPORT_TABLE_COLUMNS:
        return [], [
            f"{malformed} {header + 1}: the header has {header_cells} cells, not {SUPPORT_TABLE_COLUMNS}"
        ]
    if header + 1 >= len(lines) or not SEPARATOR_ROW_RE.match(lines[header + 1].strip()):
        return [], [f"{malformed} {header + 2}: no separator row"]
    separator_cells = len(table_cells(lines[header + 1]))
    if separator_cells != SUPPORT_TABLE_COLUMNS:
        return [], [
            f"{malformed} {header + 2}: the separator has {separator_cells} cells, not {SUPPORT_TABLE_COLUMNS}"
        ]
    rows: list[dict] = []
    errors: list[str] = []
    seen: dict[str, int] = {}
    for index in range(header + 2, len(lines)):
        line = lines[index]
        if not line.lstrip().startswith("|"):
            break
        lineno = index + 1
        cells = table_cells(line)
        if len(cells) != SUPPORT_TABLE_COLUMNS:
            errors.append(f"{malformed} {lineno}: {len(cells)} cells, not {SUPPORT_TABLE_COLUMNS}")
            continue
        link = LINK_CELL_RE.match(cells[0])
        name = link.group(1) if link else cells[0]
        status = cells[4]
        if not name:
            errors.append(f"{malformed} {lineno}: empty plugin name")
        elif name in seen:
            errors.append(f"{malformed} {lineno}: {name} repeats line {seen[name]}")
        else:
            seen[name] = lineno
        if status not in STATUS_TARGETS:
            phrases = ", ".join(repr(phrase) for phrase in STATUS_TARGETS)
            errors.append(f"{malformed} {lineno}: status {status!r} is not one of {phrases}")
        rows.append({"name": name, "status": status, "lineno": lineno})
    return rows, errors


def check_support_table(table: list[dict], rows: list[dict]) -> list[str]:
    """One error per table row without a matrix row, whose supported status names an open matrix
    row, or whose status contradicts a landed row's target (`not ported` fits any open row: vep-rs
    ships no such plugin until the row lands); one per row that sorts before the row above it, case
    ignored; and one per plugin row of the matrix (`perl_reference` equal to `<item>.pm`) the table
    lacks. A status outside STATUS_TARGETS asserts no target; `read_support_table` reports it."""
    errors: list[str] = []
    by_item = {row["item"]: row for row in rows}
    previous = ""
    for entry in table:
        name, status = entry["name"], entry["status"]
        if name.casefold() < previous.casefold():
            errors.append(
                f"{PLUGINS_PAGE}: support table row {name} is out of alphabetical order after {previous}"
            )
        previous = name
        row = by_item.get(name)
        target = STATUS_TARGETS.get(status)
        if row is None:
            errors.append(f"{PLUGINS_PAGE}: support table row {name} has no matrix row")
        elif row["state"] == "open":
            if target == "implemented":
                errors.append(
                    f"{PLUGINS_PAGE}: support table row {name} status {status} names an open matrix row"
                )
        elif target is not None and target != row["target"]:
            errors.append(
                f"{PLUGINS_PAGE}: support table row {name} status {status} "
                f"does not match target {row['target']}"
            )
    named = {entry["name"] for entry in table}
    for row in rows:
        item = row["item"]
        if row["perl_reference"] == f"{item}.pm" and item not in named:
            errors.append(
                f"{MATRIX}:{row['lineno']}: {item} is a plugin row (perl_reference {item}.pm) "
                "absent from the support table"
            )
    return errors


def read_reference_files(root: Path) -> tuple[set[str], list[str]]:
    """The lines of the reference list, or the reasons it cannot be read: the file missing, a line
    that is empty or holds a blank, or a line that does not sort after the line before it in byte
    order, which a repeated line does not."""
    path = root / REFERENCE_FILES
    if not path.is_file():
        return set(), [f"{REFERENCE_FILES} is missing"]
    errors: list[str] = []
    names: set[str] = set()
    previous: str | None = None
    lines = path.read_text(encoding="utf-8", errors="replace").split("\n")
    if lines[-1] == "":
        lines.pop()
    for lineno, line in enumerate(lines, start=1):
        if not line or any(ch.isspace() for ch in line):
            errors.append(f"{REFERENCE_FILES}:{lineno}: {line!r} is not one name")
            continue
        if previous is not None and line.encode() <= previous.encode():
            errors.append(f"{REFERENCE_FILES}:{lineno}: {line} does not sort after {previous}")
        previous = line
        names.add(line)
    return names, errors


def plugin_module(name: str) -> str:
    """The module of the plugin `name`: `<name>.pm`, or the one PLUGIN_MODULES_NAMED_OTHERWISE
    gives it."""
    return PLUGIN_MODULES_NAMED_OTHERWISE.get(name, f"{name}.pm")


def reference_error(row: dict, files: set[str], plugins: set[str]) -> str | None:
    """Why the row's `perl_reference` is in none of the forms the README gives, or None. A cell is
    `-`, a line of the reference list, or the module of a plugin of the support table, as
    `plugin_module` names it. Each but `-` may carry one identifier after `::`, a sub or a
    parameter. The module of the row's own item with nothing after it passes here whatever the
    table holds: `check_support_table` reads `<item>.pm` as a plugin row and reports the one the
    table lacks. An empty cell is `read_matrix`'s to report."""
    cell, item = row["perl_reference"], row["item"]
    modules = {plugin_module(name) for name in plugins}

    def names_a_file(name: str) -> bool:
        return name in files or name in modules

    if cell in ("", "-", plugin_module(item)) or names_a_file(cell):
        return None
    path, separator, sub = cell.rpartition("::")
    if separator and names_a_file(path):
        if REFERENCE_SUB_RE.fullmatch(sub):
            return None
        return f"perl_reference {cell}: {sub!r} after `::` is not one identifier"
    return (
        f"perl_reference {cell} names no file of {REFERENCE_FILES} "
        "and no module of a plugin of the support table"
    )


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


def strip_comments(text: str, keep_literals: bool = False) -> str:
    """The source without its line and (nested) block comments, every string, raw-string and
    char literal reduced to its delimiters around blanks of the body's length (newlines kept):
    a `/*` inside a string (`"U/*"`, an amino-acid change) or a `//` inside one opens no
    comment, and a `fn` inside a comment or a literal declares nothing. With `keep_literals`
    the literals stand as written, so the two readings of one source agree index for index
    and differ inside the literals alone."""
    out: list[str] = []

    def blanked(literal: str, open_len: int, close_len: int) -> str:
        if keep_literals:
            return literal
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


def under_test_attribute(code: str, spans: dict[int, tuple[int, str]], at: int) -> bool:
    """Whether the `fn` at index `at` of comment-stripped source stands under a test attribute,
    `spans` being the source's `attribute_spans`."""
    head = FN_KEYWORDS_RE.sub("", code[:at].rstrip()).rstrip()
    while head.endswith("]") and len(head) in spans:
        start, attribute = spans[len(head)]
        if TEST_ATTR_RE.match(attribute):
            return True
        head = code[:start].rstrip()
    return False


def declares_test_fn(text: str, fn: str) -> bool:
    """Whether the source declares `fn <fn>` under a test attribute, comments stripped first."""
    code = strip_comments(text)
    spans = attribute_spans(code)
    return any(
        under_test_attribute(code, spans, decl.start())
        for decl in re.finditer(r"\bfn\s+" + re.escape(fn) + r"\b", code)
    )


def function_bodies(code: str) -> dict[str, list[tuple[int, int, int]]]:
    """The body of every function of comment-stripped source, by the function's name and in the
    file's order: (the index of its `fn`, the index of the brace that opens the body, the index
    past the brace that closes it). The body opens at the first brace after the name that
    stands outside the parentheses and brackets of the signature; a declaration that ends at
    a `;` before one has no body. Literal bodies are blank in such source, so every brace is
    code."""
    bodies: dict[str, list[tuple[int, int, int]]] = {}
    for decl in FN_DECL_RE.finditer(code):
        nesting = 0
        opened = -1
        for token in SIGNATURE_TOKEN_RE.finditer(code, decl.end()):
            char = token.group(0)
            if char in "([":
                nesting += 1
            elif char in ")]":
                nesting -= 1
            elif nesting == 0:
                opened = token.start() if char == "{" else -1
                break
        if opened < 0:
            continue
        depth = 0
        for brace in BRACE_RE.finditer(code, opened):
            depth += 1 if brace.group(0) == "{" else -1
            if depth == 0:
                bodies.setdefault(decl.group(1), []).append((decl.start(), opened, brace.end()))
                break
    return bodies


def rust_source(text: str) -> dict:
    """A Rust file as the public-text rule reads it: its `code` (comments dropped, literal
    bodies blank), the same text with its `literals` kept, its attribute `spans` and the
    `bodies` of its functions."""
    code = strip_comments(text)
    return {
        "code": code,
        "literals": strip_comments(text, keep_literals=True),
        "spans": attribute_spans(code),
        "bodies": function_bodies(code),
    }


def reads_public_text(source: dict, fn: str) -> bool:
    """Whether a test function `fn` of a `rust_source` reads a public text: its body holds a
    string literal that is a Markdown file's path or asks clap for a help text, or calls by its
    bare name a function of the same file whose body does, through any number of such calls. A
    called function is found by its name alone, every declaration of that name being read; a
    declaration of `fn` that stands under no test attribute is no test and reads nothing."""
    code, literals, bodies = source["code"], source["literals"], source["bodies"]
    waiting = [
        (start, end)
        for at, start, end in bodies.get(fn, [])
        if under_test_attribute(code, source["spans"], at)
    ]
    followed: set[str] = set()
    while waiting:
        start, end = waiting.pop()
        body = code[start:end]
        if HELP_TEXT_RE.search(body):
            return True
        # A whole literal is its two quotes around blanks in the code and written in the kept text.
        if any(
            code[start + path.start() : start + path.end()] == '"' + " " * (len(path.group(0)) - 2) + '"'
            for path in PAGE_NAME_RE.finditer(literals[start:end])
        ):
            return True
        for call in BARE_CALL_RE.finditer(body):
            called = call.group(1)
            if called in bodies and called not in followed:
                followed.add(called)
                waiting.extend((opened, closed) for _at, opened, closed in bodies[called])
    return False


def entry_reads_public_text(root: Path, entry: str, sources: dict[str, dict | None]) -> bool:
    """Whether a `test1` entry of the crate or the tools shape names a test function that reads
    a public text. `sources` caches each file's `rust_source` by its path; an entry of another
    shape, or one whose file the tree lacks, reads none."""
    match = TEST_PATH_RE.match(entry) or TOOL_TEST_RE.match(entry)
    if not match:
        return False
    path, fn = match.groups()
    if ".." in path.split("/"):
        return False
    if path not in sources:
        file = root / path
        sources[path] = (
            rust_source(file.read_text(encoding="utf-8", errors="replace")) if file.is_file() else None
        )
    source = sources[path]
    return source is not None and reads_public_text(source, fn)


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


def shell_quotes(line: str) -> list[str]:
    """The quote that is open at each index of a shell line and, as the last item, the outermost
    one open at the line's end, read from the line's start: `"`, `'` or nothing. A backslash
    outside single quotes takes the next character with it. The text of a `$( )` inside
    double quotes is code with quotes of its own, down to the parenthesis that closes it."""
    quotes: list[str] = []
    # The quotes open at the index, innermost last, with a `(` for the `$(` that opens code
    # inside a double quote and for each parenthesis open in that code.
    opened: list[str] = []
    index = 0
    while index < len(line):
        char = line[index]
        quote = opened[-1] if opened and opened[-1] != "(" else ""
        quotes.append(quote)
        if char == "\\" and quote != "'":
            quotes.append(quote)
            index += 2
            continue
        if quote:
            if char == quote:
                opened.pop()
            elif quote == '"' and line.startswith("$(", index):
                opened.append("(")
                quotes.append(quote)
                index += 1
        elif char in "\"'" or (opened and char == "("):
            opened.append(char)
        elif opened and char == ")":
            opened.pop()
        index += 1
    return quotes[: len(line)] + [opened[0] if opened else ""]


def without_shell_comment(line: str) -> str:
    """A shell line without its comment: the text before a `#` that stands outside quotes at the
    line's start or after a blank or one of `;&|(`."""
    quotes = shell_quotes(line)
    for found in re.finditer(r"(?:^|(?<=[ \t;&|(]))#", line):
        if not quotes[found.start()]:
            return line[: found.start()]
    return line


def shell_code_lines(source: str, spanning: bool = True) -> list[str]:
    """A shell file's lines as code, one for each line of the file: a line without its comment,
    and an empty line in place of each line of a here-document, whose text is data. A
    here-document is read where its closing word stands on a later line. A quote that a line
    leaves open stays open on the lines after it: each of them is an empty line, down to the
    line that closes the quote, which keeps what follows the closing quote, after an empty
    quoted word. A file that this reading leaves inside a quote at its end is read with
    `spanning` off, every quote ending with its line."""
    lines = source.splitlines()
    code: list[str] = []
    closes_at = -1
    quote = ""
    for number, line in enumerate(lines):
        if number <= closes_at:
            code.append("")
            continue
        if quote:
            closed = 0
            while closed < len(line) and line[closed] != quote:
                closed += 2 if line[closed] == "\\" and quote == '"' else 1
            if closed >= len(line):
                code.append("")
                continue
            line = quote * 2 + line[closed + 1 :]
        bare = without_shell_comment(line)
        code.append(bare)
        quotes = shell_quotes(bare)
        quote = quotes[-1] if spanning else ""
        for opened in SHELL_HEREDOC_RE.finditer(bare):
            if quotes[opened.start()]:
                continue
            tabs, word = opened.group(1), opened.group(3) or opened.group(4)
            closing = next(
                (
                    later
                    for later in range(number + 1, len(lines))
                    if (lines[later].lstrip("\t") if tabs else lines[later]) == word
                ),
                -1,
            )
            closes_at = max(closes_at, closing)
    return shell_code_lines(source, spanning=False) if quote else code


def shell_functions(source: str) -> dict[str, str]:
    """Each function a shell file defines with its opening brace on the definition's line, and
    its body as code: the rest of that line when its braces outside quotes close on it, else
    that and the lines after it, down to the first line that opens with a closing brace at the
    definition's indentation (or to the file's end), a line that is no part of the body. The
    body of a function defined inside another is its own and no part of the outer one's. Of
    two functions of one name the first is read."""
    lines = shell_code_lines(source)
    spans: list[tuple[str, int, int, int]] = []
    for number, line in enumerate(lines):
        definition = SHELL_FN_RE.match(line)
        if not definition:
            continue
        last = number
        quotes = shell_quotes(line)
        braces = [char for at, char in enumerate(line) if char in "{}" and not quotes[at]]
        if braces.count("{") != braces.count("}"):
            closing = definition.group(1) + "}"
            last = next(
                (later for later in range(number + 1, len(lines)) if lines[later].startswith(closing)),
                len(lines),
            )
        spans.append((definition.group(2) or definition.group(3), number, last, definition.end()))
    functions: dict[str, str] = {}
    for name, first, last, start in spans:
        inner = [(begin, end) for _, begin, end, _ in spans if first < begin and end <= last]
        body = [lines[first][start:]]
        body += [
            lines[number]
            for number in range(first + 1, last)
            if not any(begin <= number <= end for begin, end in inner)
        ]
        functions.setdefault(name, "\n".join(body))
    return functions


def shell_calls(lines: list[str], called: set[str], argument: str) -> bool:
    """Whether one of the code lines calls one of the `called` functions, as a command and
    outside quotes, with a first argument that the pattern `argument` matches."""
    names = "|".join(re.escape(name) for name in sorted(called))
    call = re.compile(rf"{SHELL_COMMAND_START}(?P<name>{names})[ \t]+(?:{argument})")
    for line in lines:
        quotes: list[str] = []
        for found in call.finditer(line):
            quotes = quotes or shell_quotes(line)
            if not quotes[found.start("name")]:
                return True
    return False


def shell_check_functions(source: str) -> set[str]:
    """The check functions of a shell suite: `ok`, which prints a passed check's line, where
    the suite defines it, and every function of the suite whose body calls a check function
    with the function's own first argument, as `$1` or as a variable the body sets to `$1`, at
    the start of that call's first argument."""
    functions = shell_functions(source)
    checks = {"ok"} & functions.keys()
    grew = bool(checks)
    while grew:
        grew = False
        for name, body in functions.items():
            if name in checks:
                continue
            first = "|".join(["1", *map(re.escape, SHELL_FIRST_ARGUMENT_RE.findall(body))])
            handed_on = rf"\"?\$\{{?(?:{first})\}}?(?![A-Za-z0-9_])"
            if shell_calls(body.splitlines(), checks, handed_on):
                checks.add(name)
                grew = True
    return checks


def declares_shell_check(source: str, name: str) -> bool:
    """Whether a shell suite calls one of its check functions, as a command, with `name`, whole
    and quoted, as the first argument."""
    checks = shell_check_functions(source)
    quoted = [quote + re.escape(name) + quote for quote in "\"'" if quote not in name]
    if not checks or not quoted:
        return False
    return shell_calls(shell_code_lines(source), checks, rf"(?:{'|'.join(quoted)})(?=[\s;)|&]|$)")


def resolve_test(root: Path, entry: str, process_allowed: bool = False) -> str | None:
    """None when the entry names a test the tree declares, else the reason. The crate shape,
    `<path>.rs::<fn>` under crates/, resolves to a function under a test attribute. The three
    process shapes are accepted for `test1` only: the scripts shape,
    `scripts/<dir>/test_<name>.py::<fn>`, resolves to a pytest function of the module, the shell
    shape, `scripts/[<dir>/]test_<name>.sh::<check>`, to a check the suite names in a call of one
    of its check functions, and the tools shape, `tools/<tool>/<path>.rs::<fn>`, as the crate
    shape does."""
    script = SCRIPT_TEST_RE.match(entry) if process_allowed else None
    shell = SHELL_TEST_RE.match(entry) if process_allowed else None
    match = TEST_PATH_RE.match(entry) or (TOOL_TEST_RE.match(entry) if process_allowed else None)
    if not match and not script and not shell:
        shapes = (
            f"{CRATE_SHAPE}, {SCRIPT_SHAPE}, {SHELL_SHAPE} or {TOOL_SHAPE}"
            if process_allowed
            else CRATE_SHAPE
        )
        return f"is not {shapes}"
    path, fn = (script or shell or match).groups()
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
    elif shell:
        if not declares_shell_check(source, fn):
            return f"{path} names no check {fn!r} in a call of a check function"
    elif not declares_test_fn(source, fn):
        return f"{path} declares no test fn {fn}"
    return None


def is_process_row(row: dict) -> bool:
    """A row with a scripts or tools `test1` entry registers a release procedure, not an Ensembl VEP behaviour."""
    return any(
        SCRIPT_TEST_RE.match(entry) or SHELL_TEST_RE.match(entry) or TOOL_TEST_RE.match(entry)
        for entry in entries(row["test1"])
    )


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
        parsed = json.loads(manifest.read_text(encoding="utf-8", errors="replace"))
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
    table, errors = read_support_table(root)
    if errors:
        return errors, ""
    reference_files, errors = read_reference_files(root)
    if errors:
        return errors, ""
    plugins = {entry["name"] for entry in table}
    items = {row["item"] for row in rows}

    flags, errors = scan_args(args_rs.read_text(encoding="utf-8", errors="replace"))
    errors.extend(check_support_table(table, rows))
    declared = {flag.replace("-", "_") for flag in flags}
    for flag in flags:
        name = flag.replace("-", "_")
        if name not in items:
            errors.append(f"{ARGS_RS}: flag --{flag} has no row {name!r} in {MATRIX}")

    landed = 0
    interim = 0
    manifests: dict[str, tuple[dict | None, str | None]] = {}
    sources: dict[str, dict | None] = {}
    for row in rows:
        item, target, state = row["item"], row["target"], row["state"]
        where = f"{MATRIX}:{row['lineno']}: {item}"
        if is_process_row(row) and row["perl_reference"] != "-":
            errors.append(
                f"{where}: a process row names no Ensembl VEP module; perl_reference must be -"
            )
        reason = reference_error(row, reference_files, plugins)
        if reason:
            errors.append(f"{where}: {reason}")
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
            if (
                target == STATED_TARGET
                and entries(row["test1"])
                and not any(
                    entry_reads_public_text(root, entry, sources) for entry in entries(row["test1"])
                )
            ):
                errors.append(
                    f"{where}: landed {target} with no test1 that reads a public text: "
                    "a string literal that is a Markdown file's path or a help text asked of clap"
                )
        for entry in entries(row["test1"]):
            reason = resolve_test(root, entry, process_allowed=True)
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

    summary = (
        f"parity: {len(rows)} rows, {len(flags)} flags covered, {landed} landed, "
        f"{interim} interim entries resolve, {len(table)} support table rows"
    )
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
