#!/usr/bin/env python3
"""Hold a pull request to the admission policy and close one that misses it.

usage: admission.py --event <github event json>          the pull_request_target job
       admission.py --pull-request <number>              the same checks, run by hand (workflow_dispatch)
       admission.py --workflow-run <github event json>   the job that runs when the CI workflow completes
       admission.py --issue <github event json>          the issues job: the intended-divergences choice
       admission.py --local [--base origin/main] [--body-file <f>]
                                                        the static checks on the local branch before pushing

The static checks: the branch contains the tip of the base branch; a change to code under `crates/`,
`tools/` or `scripts/` (benchmarks and examples apart) comes with a test that gains a line (a file under
a crate's or tool's `tests/` directory, a line of a `#[cfg(test)]` module or a `#[test]` function of a
source file, a `test_*.py` or `test_*.sh` under `scripts/`, a golden `expected/` or `expected_*/` file
or a corpus `manifest.json`), unless a maintainer sets the `tests-not-needed` label; a commit declares
the model that wrote the change with a `Co-Authored-By` trailer (`No LLM` for none); the description
carries the template's three text sections and its `## Intended divergences` section with one of the
three statements checked (none of the index's entries touched, which the change may not claim while it
edits the index, a page the index links or `scripts/divergences/`; existing entries touched, one line
per entry naming it as the index's first column does and saying what the change does to it and why;
a new intended divergence proposed, with its description), the first statement excluding the other
two; a crate change adds a line under `## [Unreleased]` of CHANGELOG.md, or under a release heading
the change itself adds; a change from outside the maintainers leaves `.github/workflows/` alone; the
parity matrix resolves on the contributed tree; and every commit carries a `Signed-off-by` trailer
with its author's name and address, a miss the comment reports without closing, since the only remedy
on a pushed branch is a history rewrite after which a closed pull request cannot be reopened, and the
DCO check holds the merge anyway. The job that runs after CI adds one line per required check that
failed.

The entry names the divergences check accepts are read from `docs/intended-divergences.md` at the base
tip and at the pull request's head: the first cell of each table row that links a page (the literal
`none` and a number alone excepted), and each link target of the row with and without the `docs/`
prefix and any fragment, so an entry the pull request itself adds resolves and the check keeps working
when the index is regrouped.

The issues job reads the `Intended divergences` choice an issue form writes under its `### ` heading
and the entries or description under the next one, checks the named entries against the index on the
checked-out base branch, and asks for what is missing with the `needs-divergence-check` label and one
comment edited in place; it never closes an issue, an edit that fixes the body clears both, and an
issue opened before the rule landed is not judged.

The two Actions modes read the pull request through the REST API (commits, changed files, labels,
comments) and the contribution's files as git objects: a `git fetch` of the pull request's head ref, `git
show` for a file's text, and a `checkout-index` into a temporary directory for the parity check, whose
checker is the one of the checked-out base branch. Nothing of the contribution is checked out into the
workspace or executed. Every string a contributor controls (a path, a commit id, a matrix cell) is
rendered inside a code span, so it carries no mention, link or formatting.

One comment per pull request, edited in place on every run, states each miss with its fix. A miss closes
the pull request unless its author is a maintainer (a CODEOWNERS entry for `*`) or a bot, the
`admission-override` label is set, the pull request is a draft, or it was opened before this workflow
first landed on the base branch; the divergences miss alone is reported without a close on a pull
request opened before that rule landed, whose template had no such section; the author reopens a
closed pull request once it is fixed and the check runs again. The
newest commit that declares a model sets the `model: <name>` label (at most two names; none when
the commits declare more). Before the comment and the close are written the pull request's head and the comment
are read again, so a push in between or a run of the other job loses nothing. Exit 1 only when the run
closes the pull request, 0 otherwise, 2 on a usage or environment error.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import NamedTuple

MARKER = "<!-- admission -->"
OVERRIDE_LABEL = "admission-override"
NO_TESTS_LABEL = "tests-not-needed"
MODEL_LABEL_PREFIX = "model: "
NO_MODEL_DECLARATION = "No LLM"
LABEL_MAX = 50
LABEL_MODELS_MAX = 2
SECTIONS = ("Summary", "Verification", "Upgrade impact")
# The fourth section holds checklist statements alone, so it is read by its own check, not by the
# text test of `check_body_sections`, which would read it as empty.
DIVERGENCE_SECTION = "Intended divergences"
INDEX_PATH = "docs/intended-divergences.md"
DIVERGENCE_SCRIPTS = "scripts/divergences/"
# The rule's own landing date: the oldest base-branch commit touching a file that arrived with it.
DIVERGENCE_LANDING_PATH = "scripts/ci/fixtures/admission/divergences-none-clean.json"
# The three statements, told apart by a word each alone carries; the order breaks a tie in a
# reworded line.
STATEMENT_WORDS = (
    ("proposes", re.compile(r"\bpropos", re.IGNORECASE)),
    ("touches", re.compile(r"\bexisting\b", re.IGNORECASE)),
    ("none", re.compile(r"\bnone\b", re.IGNORECASE)),
)
NONE_STATEMENT = "This change touches none of the intended divergences."
ISSUE_MARKER = "<!-- admission:issue -->"
NEEDS_DIVERGENCE_LABEL = "needs-divergence-check"
ISSUE_CHOICE_LABEL = "Intended divergences"
ISSUE_ENTRIES_LABEL = "Entries or proposed divergence"
ISSUE_CHOICES = {
    "none": "Touches none of the intended divergences",
    "touches": "Touches existing intended divergences (named below, with the reasoning)",
    "proposes": "Proposes a new intended divergence (described below)",
}
NO_RESPONSE = "_No response_"
# The base branch ruleset names the required checks; this list stands in when it cannot be read.
REQUIRED_CHECKS = (
    "Test",
    "Clippy",
    "Format",
    "Module Docs",
    "Comparator Tests",
    "Harness Tests",
    "Advisory Reachability",
    "MSRV",
    "Check DCO",
)
FAILED_CONCLUSIONS = ("failure", "timed_out")
# A run with one of these conclusions ran no job to judge: held for a maintainer, cancelled, or
# never started.
UNJUDGED_CONCLUSIONS = (None, "cancelled", "skipped", "action_required", "startup_failure")
PARITY_CHECKER = Path("scripts/parity/check_parity_matrix.py")
CODEOWNERS = Path(".github/CODEOWNERS")
WORKFLOW_FILE = ".github/workflows/admission.yml"
API = "https://api.github.com"
PARITY_TIMEOUT_S = 600

TEST_PATH_RES = (
    re.compile(r"^(crates|tools)/[^/]+/tests/"),
    re.compile(r"^scripts/(.*/)?test_[^/]*\.(py|sh)$"),
    re.compile(r"^tests/golden/.*/(expected(_[^/]+)?/|manifest\.json$)"),
)
# A fixture or recorded-data file under `scripts/` is an input of a test, not a test.
NOT_TEST_PATH_RE = re.compile(r"^scripts/.*/(fixtures|testdata)/")
RUST_CODE_RE = re.compile(r"^(crates|tools)/.*\.rs$")
# Benchmarks and examples are code a crate change counts, with no test asked for them.
RUST_UNTESTED_RE = re.compile(r"^(crates|tools)/[^/]+/(benches|examples)/.*\.rs$")
WORKFLOW_PATH_RE = re.compile(r"^\.github/workflows/")
SCRIPT_CODE_RE = re.compile(r"^scripts/.*\.(py|sh)$")
CRATE_MANIFEST_RE = re.compile(r"^crates/[^/]+/Cargo\.toml$")
MOD_DECL_TEMPLATE = r"^\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+{name}\s*;"
HUNK_RE = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")
SIGNOFF_RE = re.compile(r"^Signed-off-by:\s*(.*?)\s*<([^<>]+)>\s*$", re.MULTILINE | re.IGNORECASE)
DECLARATION_RE = re.compile(r"^Co-Authored-By:\s*(.+?)\s*$", re.MULTILINE | re.IGNORECASE)
AUTHOR_RE = re.compile(r"^author (.*?) <([^>]*)> \d+ [+-]\d{4}$")
HEADING_RE = re.compile(r"^#{1,6}\s+(.+?)\s*#*\s*$")
# The template's checklist sits under its last heading; a checkbox line is never a section's text.
CHECKLIST_RE = re.compile(r"^\s*[-*]\s+\[[ xX]\]")
CHECKBOX_RE = re.compile(r"^\s*[-*]\s+\[([ xX])\]\s*(.*?)\s*$")
# An entry line: the entry in a code span of one to three backticks (a name holding a backtick is
# written inside a longer run; an unbounded run backtracks cubically on a line of backticks), with
# or without a list marker, then the rest of the line, read whole and stripped by the caller (a lazy
# tail with trailing whitespace is quadratic on a run of spaces).
ENTRY_LINE_RE = re.compile(r"^\s*(?:(?:[-*+]|\d+[.)])\s+)?(`{1,3})(.+?)\1(?!`)(.*)$")
LINK_RE = re.compile(r"\[([^\]]*)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
TABLE_ROW_RE = re.compile(r"^\s*\|")
FORM_HEADING_RE = re.compile(r"^### (.+?)\s*$")
FENCE_RE = re.compile(r"^(`{3,}|~{3,})")
BLOCK_RE = re.compile(
    r"<!-- admission:(static|ci) sha=([0-9a-f]+) closes=([01]) -->\n(.*?)\n<!-- /admission:\1 -->",
    re.DOTALL,
)
CFG_RE = re.compile(r"^cfg\s*\((.*)\)\s*$", re.DOTALL)
CFG_ALL_RE = re.compile(r"^all\s*\((.*)\)\s*$", re.DOTALL)
CFG_DEPTH_MAX = 16
TEST_ATTR_RE = re.compile(
    r"^(?:[A-Za-z_][A-Za-z0-9_]*::)*(?:test|rstest|test_case|quickcheck|proptest)\b"
)
STRING_PREFIXES = ("r", "br", "cr", "b", "c")
# The words that open an item after an attribute and its visibility; an attribute whose target opens
# with any other word, or sits inside the parentheses of its block, or spells a weak keyword as a
# name, is on a struct field, an enum variant, a struct literal's field or an argument, which ends
# at the first comma outside its own brackets.
ITEM_KEYWORDS = frozenset(
    (
        "fn",
        "mod",
        "impl",
        "struct",
        "enum",
        "union",
        "trait",
        "use",
        "const",
        "static",
        "type",
        "macro_rules",
        "extern",
        "unsafe",
        "async",
    )
)
VISIBILITY_KEYWORD = "pub"
# A weak keyword opens an item only with an item's shape after it: `union U`, `macro_rules!`.
WEAK_ITEM_RES = {"union": re.compile(r"\s+[^\W\d]"), "macro_rules": re.compile(r"\s*!")}
LABEL_CHARS_RE = re.compile(r"[^A-Za-z0-9 ._+/-]")


class AdmissionError(Exception):
    """An environment or API failure: the run stops, closes nothing and exits 2."""


class Finding(NamedTuple):
    check: str
    text: str
    closes: bool = True


def code(value: str, limit: int = 120) -> str:
    """A contributor-controlled value as a code span: one line, no backticks, bounded."""
    flat = re.sub(r"\s+", " ", str(value)).replace("`", "'").strip()
    if len(flat) > limit:
        flat = flat[: limit - 1] + "…"
    return f"`{flat or '-'}`"


# --- Rust test regions -----------------------------------------------------------------------------


def split_operands(predicate: str) -> list[str]:
    """The top-level comma-separated operands of a `cfg` predicate list, strings and nested
    parentheses kept whole."""
    operands: list[str] = []
    depth = 0
    current: list[str] = []
    in_string = False
    i = 0
    while i < len(predicate):
        ch = predicate[i]
        if in_string:
            current.append(ch)
            if ch == "\\" and i + 1 < len(predicate):
                current.append(predicate[i + 1])
                i += 1
            elif ch == '"':
                in_string = False
        elif ch == '"':
            in_string = True
            current.append(ch)
        elif ch == "(":
            depth += 1
            current.append(ch)
        elif ch == ")":
            depth -= 1
            current.append(ch)
        elif ch == "," and depth == 0:
            operands.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
        i += 1
    tail = "".join(current).strip()
    if tail:
        operands.append(tail)
    return operands


def cfg_is_test(attr: str) -> bool:
    """Whether an attribute's text is a `cfg` predicate true only in a test build: `cfg(test)`, or
    `cfg(all(...))` with `test` among its operands, directly or inside a nested `all(...)` up to
    `CFG_DEPTH_MAX` levels. `any(...)`, `not(...)` and a feature name containing the word are
    predicates a release build can satisfy, so they are code; a deeper nest is not parsed and
    counts as code too."""
    match = CFG_RE.match(attr.strip())
    if not match:
        return False
    pending = [(match.group(1).strip(), 0)]
    while pending:
        predicate, depth = pending.pop()
        if predicate == "test":
            return True
        inner = CFG_ALL_RE.match(predicate)
        if inner and depth < CFG_DEPTH_MAX:
            pending.extend((operand, depth + 1) for operand in split_operands(inner.group(1)))
    return False


def rust_test_regions(text: str) -> list[tuple[int, int]]:
    """Line ranges (1-based, inclusive) of the items a `#[cfg(test)]` or a test attribute (`#[test]`,
    `#[tokio::test]`, `#[rstest]`) applies to: the attribute line through the item's closing brace or
    semicolon, or, for a struct field, an enum variant, a struct literal's field or an argument,
    through its comma or its last token when no comma follows. Braces are matched outside strings,
    character literals and comments; an inner `#![cfg(test)]` makes the whole file a test region."""
    n = len(text)
    total = len(text.splitlines())
    regions: list[tuple[int, int]] = []
    i = 0
    line = 1
    depth = 0
    # The armed attribute, the one whose target has not ended: its line, and the `depth` and
    # `parens` counts where it sits.
    armed: tuple[int, int, int] | None = None
    # Whether the armed attribute's target is an item, decided by its first word with the
    # visibility skipped; None until that word arrives.
    item_target: bool | None = None
    # Generic brackets open inside the armed attribute's target: a comma inside them does not end
    # a field. `<` opens one only right after a path character, which a comparison written `a < b`
    # never has, and `>` closes one unless it ends `->` or `=>`.
    angle = 0
    parens = 0
    # The `parens` count at each open brace: an item opens only at its block's own count.
    block_parens: list[int] = []
    region_start: int | None = None
    region_depth = 0
    last_token_line = 1

    def skip_string(j: int) -> int:
        j += 1
        while j < n:
            if text[j] == "\\":
                j += 2
                continue
            if text[j] == '"':
                return j + 1
            j += 1
        return n

    def skip_raw_string(j: int) -> int:
        hashes = 0
        while j < n and text[j] == "#":
            hashes += 1
            j += 1
        if j >= n or text[j] != '"':
            return j
        end = text.find('"' + "#" * hashes, j + 1)
        return n if end < 0 else end + 1 + hashes

    while i < n:
        ch = text[i]
        is_token = not ch.isspace() and not text.startswith(("//", "/*"), i)
        if ch == "\n":
            line += 1
            i += 1
        elif text.startswith("//", i):
            end = text.find("\n", i)
            i = n if end < 0 else end
        elif text.startswith("/*", i):
            level = 1
            j = i + 2
            while j < n and level:
                if text.startswith("/*", j):
                    level += 1
                    j += 2
                elif text.startswith("*/", j):
                    level -= 1
                    j += 2
                else:
                    j += 1
            line += text.count("\n", i, j)
            i = j
        elif ch == '"':
            j = skip_string(i)
            line += text.count("\n", i, j)
            i = j
        elif ch == "'":
            if i + 1 < n and text[i + 1] == "\\":
                end = text.find("'", i + 2)
                i = n if end < 0 else end + 1
            elif i + 2 < n and text[i + 2] == "'":
                i += 3
            else:
                i += 1
        elif ch.isalpha() or ch == "_":
            j = i
            while j < n and (text[j].isalnum() or text[j] == "_"):
                j += 1
            word = text[i:j]
            if (
                armed is not None
                and region_start is None
                and item_target is None
                and parens == armed[2]
                and word != VISIBILITY_KEYWORD
            ):
                weak = WEAK_ITEM_RES.get(word)
                item_target = (
                    word in ITEM_KEYWORDS
                    and parens == (block_parens[-1] if block_parens else 0)
                    and (weak is None or weak.match(text, j) is not None)
                )
            if word in STRING_PREFIXES and j < n and text[j] in '"#':
                if word in ("b", "c"):
                    j = skip_string(j) if text[j] == '"' else j
                else:
                    j = skip_raw_string(j)
                line += text.count("\n", i, j)
            i = j
        elif ch == "#" and (text.startswith("#[", i) or text.startswith("#![", i)):
            inner = text[i + 1] == "!"
            start = i + (3 if inner else 2)
            j = start
            level = 1
            while j < n and level:
                if text[j] == '"':
                    j = skip_string(j)
                    continue
                if text[j] == "[":
                    level += 1
                elif text[j] == "]":
                    level -= 1
                j += 1
            attr = text[start : j - 1].strip()
            if inner and cfg_is_test(attr):
                return [(1, total)]
            if (
                not inner
                and armed is None
                and region_start is None
                and (cfg_is_test(attr) or TEST_ATTR_RE.match(attr))
            ):
                armed = (line, depth, parens)
                item_target = None
                angle = 0
            line += text.count("\n", i, j)
            i = j
        elif ch == "{":
            depth += 1
            block_parens.append(parens)
            if armed is not None and region_start is None and armed[1] == depth - 1:
                region_start = armed[0]
                region_depth = depth
            i += 1
        elif ch == "}":
            if block_parens:
                block_parens.pop()
            if region_start is not None and depth == region_depth:
                regions.append((region_start, line))
                region_start = None
                armed = None
            elif armed is not None and depth == armed[1] and not item_target:
                # The body's last field or variant, with no comma after it.
                regions.append((armed[0], last_token_line))
                armed = None
            depth = max(depth - 1, 0)
            if armed is not None and depth < armed[1]:
                armed = None
            i += 1
        elif ch == ";":
            if armed is not None and region_start is None and armed[1:] == (depth, parens):
                regions.append((armed[0], line))
                armed = None
            i += 1
        elif ch == ",":
            if (
                armed is not None
                and region_start is None
                and armed[1:] == (depth, parens)
                and angle == 0
                and not item_target
            ):
                regions.append((armed[0], line))
                armed = None
            i += 1
        elif ch in "([":
            parens += 1
            i += 1
        elif ch in ")]":
            parens = max(parens - 1, 0)
            if armed is not None and parens < armed[2]:
                # The attribute is on the group's last argument or element, with no comma after it.
                if region_start is None and not item_target:
                    regions.append((armed[0], last_token_line))
                armed = None
            i += 1
        elif ch == "<":
            before = text[i - 1] if i else " "
            if armed is not None and region_start is None and (before.isalnum() or before in "_:"):
                angle += 1
            i += 1
        elif ch == ">":
            if angle and text[i - 1] not in "-=":
                angle -= 1
            i += 1
        else:
            i += 1
        if is_token:
            last_token_line = line
    if region_start is not None:
        regions.append((region_start, line))
    return regions


def in_regions(lineno: int, regions: list[tuple[int, int]]) -> bool:
    return any(start <= lineno <= end for start, end in regions)


def module_parents(path: str) -> tuple[str, list[str]]:
    """The module name a Rust source file declares and the files that may declare it with `mod name;`:
    the directory's `mod.rs`, `lib.rs` or `main.rs`, or the file named after the directory."""
    parts = path.split("/")
    name = parts[-1].removesuffix(".rs")
    directory = parts[:-1]
    if name == "mod" and directory:
        name = directory[-1]
        directory = directory[:-1]
    prefix = "/".join(directory)
    candidates = [f"{prefix}/{leaf}" for leaf in ("mod.rs", "lib.rs", "main.rs")] + [f"{prefix}.rs"]
    return name, [c for c in candidates if c != path]


def declared_under_cfg_test(path: str, read: Callable[[str], str | None]) -> bool:
    """True when a parent module declares this file with `mod name;` inside a `#[cfg(test)]` region,
    which makes the whole file test code."""
    name, parents = module_parents(path)
    pattern = re.compile(MOD_DECL_TEMPLATE.format(name=re.escape(name)))
    for parent in parents:
        text = read(parent)
        if not text:
            continue
        regions = rust_test_regions(text)
        for number, line in enumerate(text.splitlines(), 1):
            if pattern.match(line) and in_regions(number, regions):
                return True
    return False


def tree_regions(read: Callable[[str], str | None], path: str) -> list[tuple[int, int]]:
    """The test regions of a Rust file in a tree: the whole file when a parent declares it under
    `#[cfg(test)]`, else the regions its own attributes mark."""
    text = read(path)
    if not text:
        return []
    if declared_under_cfg_test(path, read):
        return [(1, len(text.splitlines()))]
    return rust_test_regions(text)


# --- the changed files ---------------------------------------------------------------------------


def changed_lines(patch: str) -> tuple[list[tuple[int, str]], list[tuple[int, str]]]:
    """(removed, added) lines of a unified diff as (line number, text): removed lines numbered in the
    old file, added lines in the new one."""
    removed: list[tuple[int, str]] = []
    added: list[tuple[int, str]] = []
    old = new = None
    for raw in patch.splitlines():
        hunk = HUNK_RE.match(raw)
        if hunk:
            old, new = int(hunk.group(1)), int(hunk.group(3))
            continue
        if old is None or new is None or raw.startswith("\\"):
            continue
        tag = raw[:1]
        if tag == "+":
            added.append((new, raw[1:]))
            new += 1
        elif tag == "-":
            removed.append((old, raw[1:]))
            old += 1
        else:
            old += 1
            new += 1
    return removed, added


def is_test_path(path: str) -> bool:
    if NOT_TEST_PATH_RE.match(path):
        return False
    return any(pattern.search(path) for pattern in TEST_PATH_RES)


def asks_for_test(path: str) -> bool:
    """Code a change to which needs a test change beside it; a benchmark or an example does not."""
    return not RUST_UNTESTED_RE.match(path)


def file_changed(entry: dict) -> bool:
    return entry.get("status") in ("added", "removed") or bool(
        (entry.get("additions") or 0) + (entry.get("deletions") or 0)
    )


RegionsFor = Callable[[str, str], list[tuple[int, int]]]


def classify_file(entry: dict, regions_for: RegionsFor) -> tuple[bool, bool]:
    """(changes code, changes a test) for one entry of the pull request's file list. `regions_for(side,
    path)` gives the test regions of the file's text on the `old` or `new` side. A test change is an
    added line in a test file or a test region, or a file added or renamed into a test path; a
    deleted test, a removed test function alone, or a file moved out of a test path is none."""
    path = entry.get("filename", "")
    previous = entry.get("previous_filename") or path
    if is_test_path(path):
        arrived = entry.get("status") == "added" or (previous != path and not is_test_path(previous))
        return False, arrived or bool(entry.get("additions"))
    if not file_changed(entry):
        return False, False
    if RUST_UNTESTED_RE.match(path):
        return True, False
    if RUST_CODE_RE.match(path):
        patch = entry.get("patch")
        if patch is None:
            return True, False
        removed, added = changed_lines(patch)
        new_regions = regions_for("new", path) if added else []
        old_regions = regions_for("old", previous) if removed else []
        added_in_test = [in_regions(n, new_regions) for n, _ in added]
        removed_in_test = [in_regions(n, old_regions) for n, _ in removed]
        changes_code = not all(added_in_test) or not all(removed_in_test)
        return changes_code, any(added_in_test)
    if SCRIPT_CODE_RE.match(path):
        return True, False
    return False, False


def check_tests_included(
    files: list[dict], labels: set[str], regions_for: RegionsFor
) -> Finding | None:
    if NO_TESTS_LABEL in labels:
        return None
    changes_code = changes_test = False
    for entry in files:
        is_code, is_test = classify_file(entry, regions_for)
        changes_code = changes_code or (is_code and asks_for_test(entry.get("filename", "")))
        changes_test = changes_test or is_test
    if changes_code and not changes_test:
        return Finding(
            "tests",
            "Code under `crates/`, `tools/` or `scripts/` changes and no test gains a line: add or "
            "update a test for the change (a file under a `tests/` directory, a `#[test]` function or a "
            "`#[cfg(test)]` module, a `test_*.py` or `test_*.sh`, or the golden expectation the change "
            f"moves), or a maintainer adds the `{NO_TESTS_LABEL}` label.",
        )
    return None


def check_workflow_files(files: list[dict]) -> Finding | None:
    """A change under `.github/workflows/` from outside the maintainers: the required checks the
    after-CI job judges are produced by those files, so a contributor cannot edit, add, remove or
    rename them. The `admission-override` label holds the close, as for every finding."""
    touched = sorted(
        {
            path
            for entry in files
            for path in (entry.get("filename", ""), entry.get("previous_filename") or "")
            if WORKFLOW_PATH_RE.match(path)
        }
    )
    if not touched:
        return None
    listed = ", ".join(code(path) for path in touched[:5])
    more = f" and {len(touched) - 5} more" if len(touched) > 5 else ""
    return Finding(
        "workflows",
        f"The change touches {listed}{more} under `.github/workflows/`: workflow changes are the "
        "maintainers'; ask on an issue and leave the workflows out of this pull request.",
    )


# --- the commits ---------------------------------------------------------------------------------


def check_contains_base(contains: bool, base_ref: str) -> Finding | None:
    if contains:
        return None
    return Finding(
        "base",
        f"The branch does not contain the tip of `{base_ref}`: bring it in with "
        f"`git fetch origin && git merge --signoff origin/{base_ref}` (from a fork clone, add this "
        "repository as the `upstream` remote, then `git fetch upstream && git merge --signoff "
        f"upstream/{base_ref}`), then push.",
    )


def check_signoff(commits: list[dict]) -> Finding | None:
    """Merge commits are not judged, as the DCO check leaves them out. The miss never closes: the
    remedy on a pushed branch rewrites history, after which a closed pull request cannot be reopened."""
    unsigned = []
    for commit in commits:
        if len(commit.get("parents") or []) > 1:
            continue
        message = commit["commit"].get("message") or ""
        author = commit["commit"].get("author") or {}
        name = re.sub(r"\s+", " ", author.get("name") or "").strip().lower()
        email = (author.get("email") or "").strip().lower()
        signers = {
            (re.sub(r"\s+", " ", match.group(1)).strip().lower(), match.group(2).strip().lower())
            for match in SIGNOFF_RE.finditer(message)
        }
        if not email or (name, email) not in signers:
            unsigned.append(commit["sha"][:7])
    if not unsigned:
        return None
    listed = ", ".join(code(sha) for sha in unsigned[:5])
    more = f" and {len(unsigned) - 5} more" if len(unsigned) > 5 else ""
    return Finding(
        "signoff",
        f"Commit {listed}{more} carries no `Signed-off-by` trailer with its author's name and address: "
        "before the first push, sign every commit with `git rebase --signoff origin/main`; on a pushed "
        "branch, open a new pull request from a signed branch, since a force-push would stop this one "
        "from being reopened.",
        closes=False,
    )


def declared_models(commits: list[dict]) -> list[str]:
    """The distinct names of the `Co-Authored-By` trailers, their addresses dropped, `No LLM` normalised."""
    names: list[str] = []
    for commit in commits:
        for match in DECLARATION_RE.finditer(commit["commit"].get("message") or ""):
            name = re.sub(r"\s+", " ", re.sub(r"<[^>]*>", "", match.group(1))).strip()
            if name.lower() == NO_MODEL_DECLARATION.lower():
                name = NO_MODEL_DECLARATION
            if name and name not in names:
                names.append(name)
    return names


def check_declaration(models: list[str]) -> Finding | None:
    if models:
        return None
    return Finding(
        "declaration",
        "No commit declares the model that wrote the change: add one new signed commit whose message "
        "carries `Co-Authored-By: <model name> <email>` (or `Co-Authored-By: No LLM`); do not amend or "
        "force-push a pushed branch.",
    )


def label_models(commits: list[dict]) -> list[str] | None:
    """The names the `model:` labels are set from: the declarations of the newest commit that carries
    one (a merge of the base branch at the head declares nothing), at most `LABEL_MODELS_MAX`; None
    when more distinct names appear across the commits, so a run of unrelated trailers never
    creates repository labels."""
    if len(declared_models(commits)) > LABEL_MODELS_MAX:
        return None
    for commit in reversed(commits):
        names = declared_models([commit])
        if names:
            return names[:LABEL_MODELS_MAX]
    return []


def model_label(name: str) -> str | None:
    """`model: <name>` with the name reduced to letters, digits, spaces and `._+/-`, at most
    `LABEL_MAX` characters; None when nothing of the name survives."""
    clean = re.sub(r"\s+", " ", LABEL_CHARS_RE.sub("", name)).strip()
    if not clean:
        return None
    return (MODEL_LABEL_PREFIX + clean)[:LABEL_MAX].rstrip()


# --- the description and the changelog -----------------------------------------------------------


def opens_fence(stripped: str) -> str | None:
    """The fence a line opens, or None. A backtick run followed by text that holds another backtick
    is an inline code span, not a fence: a backtick fence's info string may not contain one."""
    opening = FENCE_RE.match(stripped)
    if opening is None:
        return None
    fence = opening.group(1)
    if fence[0] == "`" and "`" in stripped[len(fence) :]:
        return None
    return fence


def sectioned(text: str, heading_re: re.Pattern[str]) -> dict[str, list[str]]:
    """The text's sections by lower-cased heading, each the lines under it; a heading inside a
    fenced code block opens no section."""
    sections: dict[str, list[str]] = {}
    current: str | None = None
    fence: str | None = None
    for raw in text.splitlines():
        stripped = raw.strip()
        if fence is not None:
            if stripped.startswith(fence) and stripped.rstrip(fence[0]) == "":
                fence = None
        elif (opened := opens_fence(stripped)) is not None:
            fence = opened
        elif (heading := heading_re.match(stripped)) is not None:
            current = heading.group(1).strip().lower()
            sections.setdefault(current, [])
            continue
        if current is not None:
            sections[current].append(raw)
    return sections


def body_sections(body: str | None) -> dict[str, list[str]]:
    """The description's sections, HTML comments stripped first."""
    return sectioned(re.sub(r"<!--.*?-->", "", body or "", flags=re.DOTALL), HEADING_RE)


def check_body_sections(body: str | None) -> Finding | None:
    """Every template text section is present as a heading with at least one line of text under it;
    the template's HTML comments and checklist lines do not count as text."""
    sections = body_sections(body)
    missing = [
        name
        for name in SECTIONS
        if not any(
            line.strip() and not CHECKLIST_RE.match(line) for line in sections.get(name.lower(), [])
        )
    ]
    if not missing:
        return None
    listed = ", ".join(f"`## {name}`" for name in missing)
    return Finding(
        "body",
        f"The description lacks the {listed} section of the pull request template, or leaves it "
        "empty: edit the description and fill it in.",
    )


# --- the intended divergences --------------------------------------------------------------------


class Index(NamedTuple):
    """What the divergences check reads from `docs/intended-divergences.md`: the names a contributor
    may cite, normalised by `entry_name`, and the repository paths of the pages its rows link."""

    names: frozenset[str]
    pages: frozenset[str]

    def __or__(self, other: Index) -> Index:
        return Index(self.names | other.names, self.pages | other.pages)


def entry_name(text: str) -> str:
    """A cell or a cited token as the check compares it: link markup reduced to its text, code
    spans to their text, whitespace collapsed, lower-cased."""
    plain = LINK_RE.sub(lambda m: m.group(1), text).replace("`", "")
    return re.sub(r"\s+", " ", plain).strip().lower()


def table_cells(row: str) -> list[str]:
    """The cells of a table row, split on unescaped pipes, the outer pipes dropped."""
    text = row.strip().removeprefix("|")
    if text.endswith("|") and not text.endswith("\\|"):
        text = text[:-1]
    return [cell.replace("\\|", "|").strip() for cell in re.split(r"(?<!\\)\|", text)]


def is_separator_row(line: str) -> bool:
    return "-" in line and not re.sub(r"[|:\-\s]", "", line)


def page_links(row: str) -> list[str]:
    """The targets of a table row's links to pages: a fragment alone and an address are not pages."""
    targets = []
    for match in LINK_RE.finditer(row):
        target = match.group(2).removeprefix("./")
        if not target.startswith("#") and "://" not in target:
            targets.append(target)
    return targets


def read_index(text: str | None) -> Index:
    """The entry names and pages of one revision of the index. A table row that links a page is an
    entry: its names are its first cell's text (a defect letter, a difference's name, an identifier)
    and each link target with and without the `docs/` prefix and any fragment; a header row and a
    row without a page link name nothing, and neither a number alone (a class) nor the literal `none`
    a row without a defect carries is a name. A page is the repository path of a row's link target,
    relative to the `docs/` directory the index sits in."""
    names: set[str] = set()
    pages: set[str] = set()
    lines = (text or "").splitlines()
    for number, line in enumerate(lines):
        if not TABLE_ROW_RE.match(line) or is_separator_row(line):
            continue
        if number + 1 < len(lines) and is_separator_row(lines[number + 1]):
            continue
        targets = page_links(line)
        cells = table_cells(line)
        if not targets or not cells:
            continue
        name = entry_name(cells[0])
        if name and name != "none" and not name.isdigit():
            names.add(name)
        for target in targets:
            page = target.partition("#")[0]
            for variant in {target, page}:
                bare = variant.removeprefix("docs/")
                names.update({bare.lower(), f"docs/{bare}".lower()})
            if page:
                repo_path = page if page.startswith("docs/") else f"docs/{page}"
                pages.add(os.path.normpath(repo_path))
    return Index(frozenset(names), frozenset(pages))


def read_statements(lines: list[str]) -> dict[str, tuple[bool, list[str]]]:
    """The checklist statements of the divergences section by kind, each with whether it is checked
    and the lines that follow it up to the next checklist line."""
    statements: dict[str, tuple[bool, list[str]]] = {}
    current: str | None = None
    for line in lines:
        box = CHECKBOX_RE.match(line)
        if box is None:
            if current is not None:
                statements[current][1].append(line)
            continue
        current = next((kind for kind, word in STATEMENT_WORDS if word.search(box.group(2))), None)
        if current is not None:
            statements[current] = (box.group(1) != " ", [])
    return statements


def entry_lines(lines: list[str]) -> list[tuple[str, str]]:
    """(cited entry, the rest of the line) for every line of the entry shape."""
    entries = []
    for line in lines:
        match = ENTRY_LINE_RE.match(line)
        if match is not None:
            rest = re.sub(r"^[\s:,\-]+", "", match.group(3)).strip()
            entries.append((match.group(2).strip(), rest))
    return entries


def listed_codes(values: list[str]) -> str:
    listed = ", ".join(code(value) for value in values[:5])
    more = f" and {len(values) - 5} more" if len(values) > 5 else ""
    return listed + more


CHANGE_REMEDY = "what the change does to it (keeps, changes or removes it)"
REPORT_REMEDY = "what should change"


def entry_problems(
    entries: list[tuple[str, str]], index: Index, claim: str, remedy: str
) -> list[str]:
    """What is wrong with the cited entries: none named, a name the index lacks, a line that says
    nothing after the name. `claim` opens the first sentence (the statement or choice that names
    existing entries) and `remedy` is what the rest of an entry line must say: for a pull request
    what the change does to the entry, for an issue what should change."""
    if not entries:
        return [
            (
                f"{claim} and no entry line follows it: add one line per entry, ``- `<entry>`: "
                f"{remedy}, because ...``, the entry named as the first column of `{INDEX_PATH}` "
                "names it (a defect letter such as `B`, or a difference's name as written there, "
                "inner backticks dropped)."
            )
        ]
    problems = []
    unknown = [token for token, _rest in entries if entry_name(token) not in index.names]
    if unknown:
        problems.append(
            f"No entry of `{INDEX_PATH}` is named {listed_codes(unknown)}: name an entry as the "
            "index's first column does (a defect letter such as `B`, or a difference's name as "
            "written there, inner backticks dropped), or by the path of its page."
        )
    bare = [token for token, rest in entries if not rest]
    if bare:
        problems.append(
            f"The entry line for {listed_codes(bare)} says nothing after the entry: add {remedy} "
            "and why the index's reasoning is wrong for this case."
        )
    return problems


def divergence_artefacts(files: list[dict], index: Index) -> list[str]:
    """The changed paths that record the intended divergences: the index, a page one of its rows
    links, or a file under `scripts/divergences/`."""
    touched = {
        path
        for entry in files
        for path in (entry.get("filename", ""), entry.get("previous_filename") or "")
        if path
    }
    return sorted(
        path
        for path in touched
        if path == INDEX_PATH or path in index.pages or path.startswith(DIVERGENCE_SCRIPTS)
    )


def check_divergences(
    body: str | None,
    files: list[dict],
    base_index: str | None,
    head_index: str | None,
    closes: bool = True,
) -> Finding | None:
    """The `## Intended divergences` section carries one of its three statements checked: the first
    (none of the index's entries touched) excludes the other two and is refused while the change
    edits a divergence artefact; the second needs one line per entry, the entry resolving against
    the index at the base tip or the head and the line saying what the change does to it; the third
    needs its description. One finding states every problem; `closes` is False for a pull request
    opened before the rule landed, whose template had no such section."""
    index = read_index(base_index) | read_index(head_index)
    section_lines = body_sections(body).get(DIVERGENCE_SECTION.lower(), [])
    statements = read_statements(section_lines)
    checked = {kind for kind, (is_checked, _lines) in statements.items() if is_checked}
    problems: list[str] = []
    if not statements:
        problems.append(
            f"The description lacks the `## {DIVERGENCE_SECTION}` section of the pull request "
            "template, or the section carries none of its three statements: add the section and "
            "check the statement that applies (the change touches none of the index's entries; it "
            "touches existing entries, one line each; it proposes a new intended divergence)."
        )
    elif not checked:
        problems.append(
            f"No statement of the `## {DIVERGENCE_SECTION}` section is checked: check the one "
            "that applies (the change touches none of the index's entries; it touches existing "
            "entries, one line each; it proposes a new intended divergence)."
        )
    elif "none" in checked and checked != {"none"}:
        problems.append(
            f"The `## {DIVERGENCE_SECTION}` section checks the first statement (none touched) "
            "together with another: the first statement excludes the other two."
        )
    else:
        if "touches" in checked:
            problems.extend(
                entry_problems(
                    entry_lines(statements["touches"][1]),
                    index,
                    f"The second statement of the `## {DIVERGENCE_SECTION}` section is checked",
                    CHANGE_REMEDY,
                )
            )
        if "proposes" in checked and not any(line.strip() for line in statements["proposes"][1]):
            problems.append(
                f"The third statement of the `## {DIVERGENCE_SECTION}` section is checked and no "
                "description follows it: in a line or two, give the input, what Ensembl VEP and "
                "vep-rs each write, and why vep-rs's output is the one to keep."
            )
        if "none" in checked:
            artefacts = divergence_artefacts(files, index)
            if artefacts:
                problems.append(
                    f"The first statement of the `## {DIVERGENCE_SECTION}` section (none "
                    f"touched) is checked while the change edits {listed_codes(artefacts)}: name "
                    "the entries the change touches under the second statement, or describe the "
                    "proposed divergence under the third."
                )
    if not problems:
        return None
    return Finding("divergences", " ".join(problems), closes=closes)


def form_fields(body: str | None) -> dict[str, str]:
    """An issue form's fields as GitHub writes them: `### <label>` then the value up to the next
    heading, keyed by the lower-cased label; `_No response_` is an empty field."""
    fields: dict[str, str] = {}
    for label, lines in sectioned(body or "", FORM_HEADING_RE).items():
        value = "\n".join(lines).strip()
        fields[label] = "" if value == NO_RESPONSE else value
    return fields


def check_issue_divergences(body: str | None, index: Index) -> str | None:
    """What an issue is missing against the choice its form carries, or None: the choice itself
    (an issue written outside the forms), the entries the second choice names, resolved against the
    index, or the description the third choice promises."""
    fields = form_fields(body)
    choice = re.sub(r"\s+", " ", fields.get(ISSUE_CHOICE_LABEL.lower(), "")).strip().lower()
    kind = next((k for k, text in ISSUE_CHOICES.items() if text.lower() == choice), None)
    if kind is None:
        options = "; ".join(f"`{text}`" for text in ISSUE_CHOICES.values())
        return (
            f"The issue carries no `{ISSUE_CHOICE_LABEL}` choice: edit the description to add a "
            f"`### {ISSUE_CHOICE_LABEL}` heading followed by one of {options}, and under a "
            f"`### {ISSUE_ENTRIES_LABEL}` heading the entries or the description the choice calls for."
        )
    text = fields.get(ISSUE_ENTRIES_LABEL.lower(), "")
    if kind == "touches":
        problems = entry_problems(
            entry_lines(text.splitlines()),
            index,
            f"The `{ISSUE_CHOICE_LABEL}` choice names existing entries",
            REPORT_REMEDY,
        )
        if problems:
            return " ".join(problems)
    if kind == "proposes" and not text:
        return (
            f"The `{ISSUE_CHOICE_LABEL}` choice proposes a new intended divergence and the "
            f"`{ISSUE_ENTRIES_LABEL}` field is empty: in a line or two, give the input, what Ensembl "
            "VEP and vep-rs each write, and why vep-rs's output is the one to keep."
        )
    return None


def section_range(changelog: str, start: int) -> tuple[int, int]:
    """Line numbers (1-based, inclusive) of the section whose `## [` heading is line `start`, up to the
    next release heading."""
    lines = changelog.splitlines()
    end = len(lines)
    for n in range(start + 1, len(lines) + 1):
        if lines[n - 1].startswith("## ["):
            end = n - 1
            break
    return start, end


def unreleased_range(changelog: str) -> tuple[int, int] | None:
    """Line numbers of the `## [Unreleased]` section, up to the next release heading."""
    lines = changelog.splitlines()
    start = next((n for n, text in enumerate(lines, 1) if text.startswith("## [Unreleased]")), None)
    return None if start is None else section_range(changelog, start)


def is_crate_change(files: list[dict], regions_for: RegionsFor) -> bool:
    """A code change to a crate, or a change to a crate's manifest."""
    for entry in files:
        path = entry.get("filename", "")
        if CRATE_MANIFEST_RE.match(path) and file_changed(entry):
            return True
        if path.startswith("crates/") and classify_file(entry, regions_for)[0]:
            return True
    return False


def check_changelog(
    files: list[dict], head_changelog: str | None, regions_for: RegionsFor
) -> Finding | None:
    """A crate change needs an added, non-blank line inside the Unreleased section of the pull
    request's CHANGELOG.md, or a release heading the change itself adds with lines under it (a release
    moves the entries under its version, added or left in place). A CHANGELOG.md diff too large for
    the API to carry its hunks is taken as satisfying it."""
    if not is_crate_change(files, regions_for):
        return None
    finding = Finding(
        "changelog",
        "A crate changes and `CHANGELOG.md` gains no line under `## [Unreleased]`: add the entry in the "
        "grammar CONTRIBUTING.md gives for changelog entries.",
    )
    entry = next((f for f in files if f.get("filename") == "CHANGELOG.md"), None)
    if entry is None or entry.get("status") == "removed":
        return finding
    patch = entry.get("patch")
    if patch is None:
        return None if file_changed(entry) else finding
    changelog = head_changelog or ""
    lines = changelog.splitlines()
    _removed, added = changed_lines(patch)
    unreleased = unreleased_range(changelog)
    if unreleased and any(unreleased[0] < n <= unreleased[1] and text.strip() for n, text in added):
        return None
    for n, text in added:
        if text.startswith("## [") and n <= len(lines):
            start, end = section_range(changelog, n)
            if any(lines[k - 1].strip() for k in range(start + 1, end + 1)):
                return None
    return finding


# --- the parity matrix ---------------------------------------------------------------------------


def check_parity(tree: Path, checker: Path) -> Finding | None:
    """Run the base branch's checker over the contributed tree. The checker reads the tree and runs
    nothing from it. A checker that exits with a traceback was stopped by what it read, which is the
    contribution's miss; a timeout, a signal, no output at all or a failed import of the checker's
    own modules is an environment error."""
    try:
        proc = subprocess.run(
            [sys.executable, str(checker), "--root", str(tree)],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=PARITY_TIMEOUT_S,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise AdmissionError(f"the parity check did not finish in {PARITY_TIMEOUT_S} s") from exc
    if proc.returncode == 0:
        return None
    prefix = "ERROR: [parity] "
    errors = [
        line.removeprefix(prefix) for line in proc.stderr.splitlines() if line.startswith(prefix)
    ]
    if not errors:
        last = next((line for line in reversed(proc.stderr.splitlines()) if line.strip()), "")
        if proc.returncode < 0 or not last or last.startswith(("ModuleNotFoundError", "ImportError")):
            raise AdmissionError(
                f"the parity checker stopped without reading the tree (exit {proc.returncode}): "
                f"{last or 'no output'}"
            )
        return Finding(
            "parity",
            f"The parity checker could not read this change ({code(last, 200)}): run "
            f"`python3 {PARITY_CHECKER}` on the branch and fix what stops it.",
        )
    listed = "; ".join(code(line, 200) for line in errors[:5])
    more = f"; and {len(errors) - 5} more" if len(errors) > 5 else ""
    return Finding(
        "parity",
        f"The parity matrix does not resolve on this change ({listed}{more}): run "
        f"`python3 {PARITY_CHECKER}` and fix every line it prints.",
    )


# --- CI ------------------------------------------------------------------------------------------


def check_ci(jobs: list[dict], required: set[str]) -> list[Finding]:
    findings = []
    for job in jobs:
        name = job.get("name") or ""
        if name in required and job.get("conclusion") in FAILED_CONCLUSIONS:
            url = job.get("html_url") or ""
            findings.append(
                Finding("ci", f"The CI job {code(name)} failed ({url}): fix it and push.")
            )
    return findings


# --- the verdict ---------------------------------------------------------------------------------


def maintainers_from_codeowners(text: str) -> set[str]:
    """The logins owning `*` in CODEOWNERS, lower-cased; teams and addresses are not logins."""
    for raw in text.splitlines():
        parts = raw.split("#", 1)[0].split()
        if parts and parts[0] == "*":
            return {p[1:].lower() for p in parts[1:] if p.startswith("@") and "/" not in p}
    return set()


def read_maintainers(repo_dir: Path) -> set[str]:
    """The maintainers named by the base checkout's CODEOWNERS; none when it has no such file."""
    codeowners = repo_dir / CODEOWNERS
    return maintainers_from_codeowners(codeowners.read_text() if codeowners.is_file() else "")


def is_bot(user: dict) -> bool:
    """A GitHub App's user: type `Bot`, or a login of the `name[bot]` or `app/name` form."""
    login = user.get("login") or ""
    return user.get("type") == "Bot" or login.endswith("[bot]") or login.startswith("app/")


def is_maintainer(user: dict, maintainers: set[str]) -> bool:
    """By CODEOWNERS login: the author association of a maintainer can read COLLABORATOR."""
    return (user.get("login") or "").lower() in maintainers


def predates(created_at: str | None, landed_at: str | None) -> bool:
    """True when the pull request was opened before the workflow first landed on the base branch."""
    if not created_at or not landed_at:
        return False
    return datetime.fromisoformat(created_at) < datetime.fromisoformat(landed_at)


REPORTED_ONLY = (
    (
        "signoff",
        (
            "the sign-off miss is reported only, and the DCO check holds the merge until every "
            "commit is signed off"
        ),
    ),
    (
        "divergences",
        (
            f"the `## {DIVERGENCE_SECTION}` miss is reported only, since the pull request was opened "
            "before that section was asked"
        ),
    ),
)


def close_footer(
    findings: list[Finding],
    user: dict,
    labels: set[str],
    draft: bool,
    maintainers: set[str],
    grandfathered: bool = False,
) -> tuple[bool, str | None]:
    """(close it, the sentence that ends the comment). Findings on an exempt pull request, and
    findings that are reported only (a sign-off miss; the divergences miss on a pull request older
    than that rule), are stated with the reason the check does not close."""
    if not findings:
        return False, None
    if is_maintainer(user, maintainers):
        return False, "No automatic close: the author is a maintainer."
    if is_bot(user):
        return False, "No automatic close: the author is a bot."
    if OVERRIDE_LABEL in labels:
        return False, f"No automatic close: the `{OVERRIDE_LABEL}` label is set."
    if draft:
        return False, "No automatic close: the pull request is a draft."
    if grandfathered:
        return False, "No automatic close: the pull request was opened before the admission check."
    if not any(finding.closes for finding in findings):
        checks = {finding.check for finding in findings}
        reasons = [reason for check, reason in REPORTED_ONLY if check in checks] or [
            "the misses above are reported only"
        ]
        return False, "No automatic close: " + "; ".join(reasons) + "."
    return True, (
        "The pull request is closed until these are fixed: push the fixes to the same branch, then "
        "reopen it and the check runs again."
    )


# --- the comment ---------------------------------------------------------------------------------


Block = tuple[str, str, bool]


def parse_blocks(body: str | None) -> dict[str, Block]:
    """{kind: (head sha, text, closes)} for the blocks an earlier run wrote into the comment."""
    return {
        m.group(1): (m.group(2), m.group(4), m.group(3) == "1")
        for m in BLOCK_RE.finditer(body or "")
    }


def static_block(
    findings: list[Finding], head_sha: str, policy_url: str, note: str | None = None
) -> str:
    head = code(head_sha[:7])
    if not findings:
        lines = [f"The admission check ([policy]({policy_url})) passes on head {head}."]
    else:
        lines = [f"The admission check ([policy]({policy_url})) finds on head {head}:"]
        lines.extend(f"- {finding.text}" for finding in findings)
    if note:
        lines.append(note)
    return "\n".join(lines)


LABELS_SKIPPED_NOTE = (
    f"The `{MODEL_LABEL_PREFIX.strip()}` label is not set: the commits declare more than "
    f"{LABEL_MODELS_MAX} distinct models; the newest declaring commit is the one a label would name."
)


def ci_block(findings: list[Finding], head_sha: str) -> str:
    head = code(head_sha[:7])
    if not findings:
        return f"Every required CI check passes on head {head}."
    lines = [f"After CI on head {head}:"]
    lines.extend(f"- {finding.text}" for finding in findings)
    return "\n".join(lines)


def render_comment(blocks: dict[str, Block], head_sha: str, footer: str | None) -> str:
    """The comment body: the marker, the blocks written for this head, the closing sentence."""
    parts = [MARKER]
    for kind in ("static", "ci"):
        if kind in blocks and blocks[kind][0] == head_sha:
            sha, text, closes = blocks[kind]
            parts.append(
                f"<!-- admission:{kind} sha={sha} closes={int(closes)} -->\n"
                f"{text}\n"
                f"<!-- /admission:{kind} -->"
            )
    if footer:
        parts.append(footer)
    return "\n\n".join(parts)


def block_has_findings(block: Block | None) -> bool:
    return block is not None and "\n- " in block[1]


# --- GitHub --------------------------------------------------------------------------------------


class GitHub:
    """The REST API over urllib: no dependency to install in the write-token job."""

    def __init__(self, repo: str, token: str | None, api: str = API):
        self.repo = repo
        self.token = token
        self.api = api

    def request(
        self,
        method: str,
        path: str,
        body: dict | None = None,
        params: dict | None = None,
    ):
        url = path if path.startswith("http") else self.api + path
        if params:
            url += ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
        data = json.dumps(body).encode() if body is not None else None
        headers = {"Accept": "application/vnd.github+json"}
        if self.token:
            headers["Authorization"] = f"Bearer {self.token}"
        if data is not None:
            headers["Content-Type"] = "application/json"
        for attempt in (1, 2):
            # Only a GET is retried: a write the server applied before the error would run twice.
            retry = attempt == 1 and method == "GET"
            req = urllib.request.Request(url, data=data, method=method, headers=headers)
            try:
                with urllib.request.urlopen(req, timeout=60) as resp:
                    payload = resp.read()
                    return (
                        resp.status,
                        dict(resp.headers),
                        (json.loads(payload) if payload else None),
                    )
            except urllib.error.HTTPError as exc:
                if retry and exc.code >= 500:
                    time.sleep(5)
                    continue
                detail = exc.read().decode(errors="replace")[:300]
                raise AdmissionError(f"{method} {url} -> {exc.code}: {detail}") from exc
            except urllib.error.URLError as exc:
                if retry:
                    time.sleep(5)
                    continue
                raise AdmissionError(f"{method} {url}: {exc.reason}") from exc
        raise AdmissionError(f"{method} {url}: no response")

    def get(self, path: str, params: dict | None = None):
        return self.request("GET", path, params=params)[2]

    def paginate(self, path: str, params: dict | None = None) -> list:
        """Every item of a list endpoint. The jobs of a workflow run arrive wrapped as
        `{"total_count": n, "jobs": [...]}`; a page of that shape yields its list."""
        items: list = []
        params = dict(params or {}, per_page=100)
        url: str | None = path
        while url:
            _status, headers, page = self.request("GET", url, params=params)
            params = None
            if isinstance(page, dict):
                page = next((value for value in page.values() if isinstance(value, list)), [])
            items.extend(page or [])
            link = headers.get("Link") or headers.get("link") or ""
            next_page = re.search(r"<([^>]+)>;\s*rel=\"next\"", link)
            url = next_page.group(1) if next_page else None
        return items

    def post(self, path: str, body: dict):
        return self.request("POST", path, body=body)[2]

    def patch(self, path: str, body: dict):
        return self.request("PATCH", path, body=body)[2]

    def delete(self, path: str):
        return self.request("DELETE", path)[2]


def find_comment(gh: GitHub, number: int, marker: str = MARKER) -> dict | None:
    for comment in gh.paginate(f"/repos/{gh.repo}/issues/{number}/comments"):
        if (comment.get("body") or "").startswith(marker) and is_bot(comment.get("user") or {}):
            return comment
    return None


def write_comment(gh: GitHub, number: int, existing: dict | None, body: str) -> None:
    if existing is not None:
        if existing.get("body") != body:
            gh.patch(f"/repos/{gh.repo}/issues/comments/{existing['id']}", {"body": body})
    else:
        gh.post(f"/repos/{gh.repo}/issues/{number}/comments", {"body": body})


def sync_model_labels(gh: GitHub, number: int, labels: set[str], models: list[str]) -> None:
    """One `model: <name>` label per declared model; a label from an earlier declaration goes."""
    wanted = {label for label in map(model_label, models) if label is not None}
    present = {label for label in labels if label.startswith(MODEL_LABEL_PREFIX)}
    try:
        for stale in sorted(present - wanted):
            gh.delete(
                f"/repos/{gh.repo}/issues/{number}/labels/{urllib.parse.quote(stale, safe='')}"
            )
        if wanted - present:
            gh.post(
                f"/repos/{gh.repo}/issues/{number}/labels",
                {"labels": sorted(wanted - present)},
            )
    except AdmissionError as exc:
        print(
            f"WARNING: [admission] the model label was not written: {exc}",
            file=sys.stderr,
        )


def landed_at(gh: GitHub, base_ref: str, path: str = WORKFLOW_FILE) -> str | None:
    """The committer date of the first commit on the base branch that carries `path`: the workflow
    file for the check as a whole, a fixture of the divergences rule for that rule."""
    commits = gh.paginate(f"/repos/{gh.repo}/commits", {"path": path, "sha": base_ref})
    if not commits:
        return None
    # The API lists commits newest first, so the last one is the oldest.
    return ((commits[-1].get("commit") or {}).get("committer") or {}).get("date")


def opened_before(gh: GitHub, item: dict, path: str, base_ref: str | None = None) -> bool:
    """Whether the pull request or issue predates the landing of `path` on the base branch (a pull
    request's own base, or `base_ref` for an issue); a failed read of the landing date counts as
    yes, so an API fault never closes a pull request or asks an issue for anything."""
    ref = base_ref or item["base"]["ref"]
    try:
        return predates(item.get("created_at"), landed_at(gh, ref, path))
    except AdmissionError as exc:
        print(f"WARNING: [admission] the landing date was not read ({path}): {exc}", file=sys.stderr)
        return True


def required_checks(gh: GitHub, base_ref: str) -> set[str]:
    try:
        rules = gh.get(f"/repos/{gh.repo}/rules/branches/{urllib.parse.quote(base_ref, safe='')}")
    except AdmissionError as exc:
        print(f"WARNING: [admission] the ruleset was not read: {exc}", file=sys.stderr)
        return set(REQUIRED_CHECKS)
    names = {
        check.get("context")
        for rule in rules or []
        if rule.get("type") == "required_status_checks"
        for check in (rule.get("parameters") or {}).get("required_status_checks") or []
    }
    return {name for name in names if name} or set(REQUIRED_CHECKS)


# --- git -----------------------------------------------------------------------------------------


def parse_commit_object(sha: str, raw: bytes) -> dict:
    """A commit object's headers and message as the API's commit entry: the parents, the author's
    name and address, and the message with its trailers."""
    text = raw.decode("utf-8", errors="replace")
    header, _, message = text.partition("\n\n")
    parents: list[dict] = []
    name = email = ""
    for line in header.splitlines():
        if line.startswith("parent "):
            parents.append({"sha": line.split()[1]})
        elif line.startswith("author "):
            match = AUTHOR_RE.match(line)
            if match:
                name, email = match.group(1), match.group(2)
    return {
        "sha": sha,
        "parents": parents,
        "commit": {"message": message, "author": {"name": name, "email": email}},
    }


class Git:
    def __init__(self, repo_dir: Path):
        self.repo_dir = repo_dir

    def _git(self, *args: str, env: dict | None = None) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["git", "-C", str(self.repo_dir), *args],
            capture_output=True,
            text=True,
            errors="replace",
            env=env,
            check=False,
        )

    def run(self, *args: str, check: bool = True, env: dict | None = None) -> str:
        proc = self._git(*args, env=env)
        if check and proc.returncode != 0:
            raise AdmissionError(f"git {' '.join(args[:3])} failed: {proc.stderr.strip()[-300:]}")
        return proc.stdout

    def rev_parse(self, rev: str) -> str | None:
        out = self.run("rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}", check=False).strip()
        return out or None

    def is_ancestor(self, ancestor: str, descendant: str) -> bool:
        return self._git("merge-base", "--is-ancestor", ancestor, descendant).returncode == 0

    def merge_base(self, a: str, b: str) -> str | None:
        return self.run("merge-base", a, b, check=False).strip() or None

    def show(self, rev: str, path: str) -> str | None:
        proc = self._git("show", f"{rev}:{path}")
        return proc.stdout if proc.returncode == 0 else None

    def fetch_pull_head(self, number: int) -> str:
        """Fetch the pull request's head ref into the object store and return its commit id."""
        ref = f"refs/admission/pull/{number}/head"
        self.run("fetch", "--quiet", "origin", f"+refs/pull/{number}/head:{ref}")
        sha = self.rev_parse(ref)
        if sha is None:
            raise AdmissionError(f"refs/pull/{number}/head was not fetched")
        return sha

    def log_commits(self, base: str, head: str) -> list[dict]:
        """The commits reachable from `head` and not from `base`, parents before children, in the
        shape the API gives a pull request's commits (the fields the checks read). The objects come
        through `cat-file --batch`, whose size-prefixed records a message cannot break."""
        shas = self.run("rev-list", "--topo-order", "--reverse", f"{base}..{head}").split()
        if not shas:
            return []
        proc = subprocess.run(
            ["git", "-C", str(self.repo_dir), "cat-file", "--batch"],
            input=("\n".join(shas) + "\n").encode(),
            capture_output=True,
            check=False,
        )
        if proc.returncode != 0:
            raise AdmissionError(f"git cat-file failed: {proc.stderr.decode(errors='replace')[-300:]}")
        commits = []
        data = proc.stdout
        pos = 0
        for sha in shas:
            header_end = data.index(b"\n", pos)
            fields = data[pos:header_end].decode(errors="replace").split()
            if len(fields) != 3 or fields[1] != "commit":
                raise AdmissionError(f"{sha} is not a commit object: {' '.join(fields[1:])}")
            size = int(fields[2])
            body_start = header_end + 1
            body_end = body_start + size
            commits.append(parse_commit_object(sha, data[body_start:body_end]))
            pos = body_end + 1
        return commits

    def export_tree(self, rev: str, dest: Path) -> None:
        """Write the tree at `rev` under `dest` through a temporary index: no hook runs, the tree's
        `.gitattributes` can name a filter but not define its command, the runner's own git-lfs
        smudge is skipped, and a symbolic link is written as a plain file."""
        dest.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, GIT_INDEX_FILE=str(Path(tmp) / "index"), GIT_LFS_SKIP_SMUDGE="1")
            self.run("read-tree", rev, env=env)
            self.run(
                "-c",
                "core.symlinks=false",
                "-c",
                "filter.lfs.required=false",
                "checkout-index",
                "--all",
                "--force",
                f"--prefix={dest}/",
                env=env,
            )


def git_regions_for(git: Git, head_sha: str, merge_base: str | None) -> RegionsFor:
    def regions_for(side: str, path: str) -> list[tuple[int, int]]:
        rev = head_sha if side == "new" else merge_base
        if not rev:
            return []
        return tree_regions(lambda p: git.show(rev, p), path)

    return regions_for


# --- modes ---------------------------------------------------------------------------------------


def policy_url(repo: str, base_ref: str) -> str:
    return f"https://github.com/{repo}/blob/{base_ref}/CONTRIBUTING.md#pull-requests"


def static_findings(
    contains: bool,
    base_ref: str,
    commits: list[dict],
    files: list[dict],
    labels: set[str],
    body: str | None,
    regions_for: RegionsFor,
    head_changelog: str | None,
    parity: Finding | None,
    models: list[str],
    base_index: str | None = None,
    head_index: str | None = None,
    divergences_close: bool = True,
) -> list[Finding]:
    """The static findings in comment order. `base_index` and `head_index` are the text of
    `docs/intended-divergences.md` at the base tip and the head; `divergences_close` is False for a
    pull request opened before the divergences rule landed, whose miss is then reported only."""
    checks = [
        check_contains_base(contains, base_ref),
        check_tests_included(files, labels, regions_for),
        check_signoff(commits),
        check_declaration(models),
        check_body_sections(body),
        check_divergences(body, files, base_index, head_index, closes=divergences_close),
        check_changelog(files, head_changelog, regions_for),
        parity,
    ]
    return [finding for finding in checks if finding is not None]


def report(findings: list[Finding]) -> None:
    for finding in findings:
        print(f"ERROR: [admission] {finding.text}", file=sys.stderr)


def verdict(findings: list[Finding], closed: bool) -> int:
    """The job's exit status: 1 only when this run closed the pull request, so the check reads red
    for a close and never for a comment alone (an exempt author, the sign-off miss)."""
    report(findings)
    if not findings:
        print("admission: OK")
        return 0
    print(f"admission: {len(findings)} miss(es)" + (", closed" if closed else ", not closed"))
    return 1 if closed else 0


def label_names(pr: dict) -> set[str]:
    return {label.get("name", "") for label in pr.get("labels") or []}


class Publication(NamedTuple):
    """`closed`: this run closed the pull request. `current`: the pull request's head as re-read
    before the write matched the judged head, so the comment and the close were decided here and
    the labels may follow. `labels`: the pull request's labels from the same re-read (empty when
    the head moved)."""

    closed: bool
    current: bool
    labels: set[str]


def publish_findings(
    gh: GitHub,
    pr: dict,
    repo_dir: Path,
    kind: str,
    block: str,
    findings: list[Finding],
    noted: bool = False,
) -> Publication:
    """Write this run's block into the pull request's comment, close the pull request when the
    findings call for it and it is open, and return the `Publication`: nothing is written when the
    head moved since the judgement. The footer is the close or exemption sentence for the findings
    on this head: this run's, else the other block's when an earlier run on the same head left
    some. A clean run writes only to refresh a comment that exists, unless its block carries a note
    (`noted`)."""
    number = pr["number"]
    head_sha = pr["head"]["sha"]
    maintainers = read_maintainers(repo_dir)
    user = pr.get("user") or {}
    labels = label_names(pr)
    draft = bool(pr.get("draft"))
    older = opened_before(gh, pr, WORKFLOW_FILE)
    close, footer = close_footer(findings, user, labels, draft, maintainers, older)
    # The judgement is of one head: a push since the event moved it, and the run for the new head
    # writes the comment and decides the close; this run then writes nothing.
    current = gh.get(f"/repos/{gh.repo}/pulls/{number}") or {}
    if (current.get("head") or {}).get("sha", head_sha) != head_sha:
        print(f"admission: the head moved past {head_sha[:7]}; this run writes nothing.")
        return Publication(closed=False, current=False, labels=set())
    current_labels = label_names(current) if current.get("labels") is not None else labels
    # The comment is read right before it is written and the other job's block is kept, so two runs
    # that interleave lose neither block.
    existing = find_comment(gh, number)
    blocks = parse_blocks(existing.get("body") if existing else None)
    blocks[kind] = (head_sha, block, any(f.closes for f in findings))
    other = blocks.get("ci" if kind == "static" else "static")
    if not footer and block_has_findings(other) and other[0] == head_sha:
        _close, footer = close_footer(
            [Finding("other", "", closes=other[2])], user, labels, draft, maintainers, older
        )
    if existing is not None or findings or noted:
        write_comment(gh, number, existing, render_comment(blocks, head_sha, footer))
    if close and (current.get("state") or pr.get("state")) == "open":
        gh.patch(f"/repos/{gh.repo}/pulls/{number}", {"state": "closed"})
        return Publication(closed=True, current=True, labels=current_labels)
    return Publication(closed=False, current=True, labels=current_labels)


def run_event(event: dict, gh: GitHub, git: Git, repo_dir: Path) -> int:
    pr = event["pull_request"]
    number = pr["number"]
    head_sha = pr["head"]["sha"]
    base_ref = pr["base"]["ref"]
    fetched = git.fetch_pull_head(number)
    if fetched != head_sha:
        print(f"admission: the head moved to {fetched[:7]}; the run for that head reports instead.")
        return 0
    base_tip = git.rev_parse(f"origin/{base_ref}") or pr["base"]["sha"]
    contains = git.is_ancestor(base_tip, head_sha)
    merge_base = git.merge_base(base_tip, head_sha)
    # The commits come from the fetched objects: the API's list stops at 250.
    commits = git.log_commits(merge_base or base_tip, head_sha)
    files = gh.paginate(f"/repos/{gh.repo}/pulls/{number}/files")
    models = declared_models(commits)
    labels = label_names(pr)
    maintainers = read_maintainers(repo_dir)
    findings = static_findings(
        contains,
        base_ref,
        commits,
        files,
        labels,
        pr.get("body"),
        git_regions_for(git, head_sha, merge_base),
        git.show(head_sha, "CHANGELOG.md"),
        None,
        models,
        base_index=git.show(base_tip, INDEX_PATH),
        head_index=git.show(head_sha, INDEX_PATH),
        divergences_close=not opened_before(gh, pr, DIVERGENCE_LANDING_PATH),
    )
    if not is_maintainer(pr.get("user") or {}, maintainers):
        workflows = check_workflow_files(files)
        if workflows is not None:
            findings.append(workflows)
    # The parity finding is appended last rather than passed to static_findings, so the workflow
    # finding precedes it.
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "tree"
        git.export_tree(head_sha, tree)
        parity = check_parity(tree, repo_dir / PARITY_CHECKER)
    if parity is not None:
        findings.append(parity)
    labelled = label_models(commits)
    note = LABELS_SKIPPED_NOTE if labelled is None else None
    block = static_block(findings, head_sha, policy_url(gh.repo, base_ref), note)
    published = publish_findings(gh, pr, repo_dir, "static", block, findings, noted=note is not None)
    if published.current and labelled is not None:
        sync_model_labels(gh, number, published.labels, labelled)
    print(f"admission: {len(commits)} commits, {len(files)} files, head {head_sha[:7]}")
    return verdict(findings, published.closed)


def run_dispatch(number: int, gh: GitHub, git: Git, repo_dir: Path) -> int:
    """The static checks on one pull request, run by hand through `workflow_dispatch`."""
    pr = gh.get(f"/repos/{gh.repo}/pulls/{number}")
    return run_event({"action": "dispatch", "pull_request": pr}, gh, git, repo_dir)


def find_pull_request(gh: GitHub, run: dict) -> dict | None:
    """The open pull request whose head is the run's head commit. The run's own `pull_requests` list
    is empty for a fork, so the head repository and branch, then the commit, are asked."""
    head_sha = run.get("head_sha")
    for ref in run.get("pull_requests") or []:
        pr = gh.get(f"/repos/{gh.repo}/pulls/{ref['number']}")
        if pr.get("head", {}).get("sha") == head_sha and pr.get("state") == "open":
            return pr
    head_repo = run.get("head_repository") or {}
    owner = (head_repo.get("owner") or {}).get("login")
    branch = run.get("head_branch")
    if owner and branch:
        by_branch = gh.get(
            f"/repos/{gh.repo}/pulls",
            {"state": "open", "head": f"{owner}:{branch}", "per_page": 20},
        )
        for pr in by_branch or []:
            if pr.get("head", {}).get("sha") == head_sha:
                return pr
    by_commit = gh.get(f"/repos/{gh.repo}/commits/{head_sha}/pulls", {"per_page": 20})
    for pr in by_commit or []:
        if pr.get("head", {}).get("sha") == head_sha and pr.get("state") == "open":
            return pr
    return None


def resolve_pull_request(event: dict, gh: GitHub) -> int:
    """Print the number of the open pull request a completed CI run belongs to, or nothing: the
    first job of the after-CI workflow, whose output keys the judging job's concurrency group."""
    run = event["workflow_run"]
    if run.get("event") != "pull_request" or run.get("conclusion") in UNJUDGED_CONCLUSIONS:
        return 0
    pr = find_pull_request(gh, run)
    if pr is not None:
        print(pr["number"])
    return 0


def run_workflow_run(event: dict, gh: GitHub, repo_dir: Path) -> int:
    run = event["workflow_run"]
    if run.get("event") != "pull_request" or run.get("conclusion") in UNJUDGED_CONCLUSIONS:
        print("admission: not a completed pull request run; nothing to judge.")
        return 0
    head_sha = run["head_sha"]
    pr = find_pull_request(gh, run)
    if pr is None:
        print(f"admission: no open pull request has head {head_sha[:7]}; nothing to judge.")
        return 0
    jobs = gh.paginate(f"/repos/{gh.repo}/actions/runs/{run['id']}/jobs")
    findings = check_ci(jobs, required_checks(gh, pr["base"]["ref"]))
    published = publish_findings(gh, pr, repo_dir, "ci", ci_block(findings, head_sha), findings)
    print(f"admission: CI run {run['id']} judged on head {head_sha[:7]}")
    return verdict(findings, published.closed)


def run_issue(event: dict, gh: GitHub, repo_dir: Path) -> int:
    """The issues job: hold an open issue to the intended-divergences choice its form carries, with
    the `needs-divergence-check` label and one comment edited in place, never a close; an edit that
    clears the miss removes the label and rewrites the comment. An issue opened before the rule
    landed on the default branch is not judged. The index read is the checked-out base branch's;
    the labels and state are re-read before the writes, so a run queued behind another never acts
    on a stale payload. Exit 0 whatever the outcome: nothing here blocks a merge."""
    issue = event["issue"]
    number = issue["number"]
    if "pull_request" in issue or issue.get("state") != "open":
        print(f"admission: #{number} is not an open issue; nothing to judge.")
        return 0
    default_branch = (event.get("repository") or {}).get("default_branch") or "main"
    if opened_before(gh, issue, DIVERGENCE_LANDING_PATH, default_branch):
        print(f"admission: issue #{number} was opened before the check; nothing to judge.")
        return 0
    index_file = repo_dir / INDEX_PATH
    index = read_index(index_file.read_text(encoding="utf-8") if index_file.is_file() else None)
    problem = check_issue_divergences(issue.get("body"), index)
    current = gh.get(f"/repos/{gh.repo}/issues/{number}") or issue
    if current.get("state") != "open":
        print(f"admission: issue #{number} is closed; nothing to write.")
        return 0
    labels = label_names(current)
    existing = find_comment(gh, number, ISSUE_MARKER)
    if problem is not None:
        print(f"ERROR: [admission] {problem}", file=sys.stderr)
        body = f"{ISSUE_MARKER}\n{problem} Edit the issue and the check runs again."
        write_comment(gh, number, existing, body)
        if NEEDS_DIVERGENCE_LABEL not in labels:
            gh.post(f"/repos/{gh.repo}/issues/{number}/labels", {"labels": [NEEDS_DIVERGENCE_LABEL]})
        print(f"admission: issue #{number}: 1 miss, labelled `{NEEDS_DIVERGENCE_LABEL}`")
        return 0
    if existing is not None:
        write_comment(gh, number, existing, f"{ISSUE_MARKER}\nThe intended-divergences check passes.")
    if NEEDS_DIVERGENCE_LABEL in labels:
        quoted = urllib.parse.quote(NEEDS_DIVERGENCE_LABEL, safe="")
        gh.delete(f"/repos/{gh.repo}/issues/{number}/labels/{quoted}")
    print(f"admission: issue #{number}: OK")
    return 0


LOCAL_STATUS = {
    "A": "added",
    "M": "modified",
    "D": "removed",
    "R": "renamed",
    "C": "added",
    "T": "modified",
}


def local_files(git: Git, merge_base: str) -> list[dict]:
    """The pull request's file list as the API would give it, from `git diff` against the merge base."""
    fields = git.run("diff", "--name-status", "-M", "-z", merge_base, "HEAD").split("\0")
    files = []
    i = 0
    while i < len(fields) and fields[i]:
        status = fields[i]
        kind = LOCAL_STATUS.get(status[0], "modified")
        previous = None
        if status[0] in "RC":
            previous, path = fields[i + 1], fields[i + 2]
            i += 3
        else:
            path = fields[i + 1]
            i += 2
        paths = [p for p in (previous, path) if p]
        diff = git.run(
            "diff",
            "--no-color",
            "--no-ext-diff",
            "-M",
            merge_base,
            "HEAD",
            "--",
            *paths,
        )
        hunks = diff[diff.find("\n@@") + 1 :] if "\n@@" in diff else ""
        hunk_lines = hunks.splitlines()
        entry = {
            "filename": path,
            "status": kind,
            "additions": sum(1 for line in hunk_lines if line.startswith("+")),
            "deletions": sum(1 for line in hunk_lines if line.startswith("-")),
        }
        if hunks:
            entry["patch"] = hunks
        if previous:
            entry["previous_filename"] = previous
        files.append(entry)
    return files


def run_local(repo_dir: Path, base: str, body_file: Path | None) -> int:
    git = Git(repo_dir)
    head = git.rev_parse("HEAD")
    base_sha = git.rev_parse(base)
    if head is None or base_sha is None:
        remote = base.split("/", 1)[0] if "/" in base else "origin"
        raise AdmissionError(f"cannot resolve HEAD or {base}; run `git fetch {remote}` first")
    merge_base = git.merge_base(base_sha, head) or base_sha
    commits = git.log_commits(base, "HEAD")
    files = local_files(git, merge_base)
    models = declared_models(commits)
    # Without a description to check, the synthesised one fills the text sections and certifies that
    # the change touches no intended divergence, which the artefact test still holds it to.
    filled_template = "\n".join(f"## {name}\nx" for name in SECTIONS)
    filled_template += f"\n## {DIVERGENCE_SECTION}\n- [x] {NONE_STATEMENT}\n"
    body = body_file.read_text() if body_file else filled_template
    findings = static_findings(
        git.is_ancestor(base_sha, head),
        base.split("/", 1)[-1],
        commits,
        files,
        set(),
        body,
        git_regions_for(git, head, merge_base),
        git.show(head, "CHANGELOG.md"),
        check_parity(repo_dir, repo_dir / PARITY_CHECKER),
        models,
        base_index=git.show(base_sha, INDEX_PATH),
        head_index=git.show(head, INDEX_PATH),
    )
    report(findings)
    workflows = check_workflow_files(files)
    if workflows is not None:
        print(
            f"WARNING: [admission] {workflows.text} From outside the maintainers this closes the "
            f"pull request unless a maintainer sets the `{OVERRIDE_LABEL}` label.",
            file=sys.stderr,
        )
    if findings:
        print(f"admission: {len(findings)} miss(es) against {base}; fix them before pushing.")
        return 1
    note = "" if body_file else "; the description is checked on the pull request"
    print(
        f"admission: OK ({len(commits)} commits, {len(files)} files against {base}, model {', '.join(models)}{note})"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--event", metavar="JSON", help="the pull_request_target event payload")
    mode.add_argument("--pull-request", metavar="N", type=int, help="one pull request, by number")
    mode.add_argument("--workflow-run", metavar="JSON", help="the workflow_run event payload")
    mode.add_argument("--issue", metavar="JSON", help="the issues event payload")
    mode.add_argument("--local", action="store_true", help="check the local branch before pushing")
    parser.add_argument("--repo-dir", default=".", help="the base repository checkout (default: .)")
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"), help="owner/name")
    parser.add_argument(
        "--base",
        default="origin/main",
        help="--local: the base ref (default: origin/main)",
    )
    parser.add_argument("--body-file", help="--local: the pull request description to check")
    parser.add_argument(
        "--resolve-pull-request",
        action="store_true",
        help="--workflow-run: print the open pull request's number for the run's head, or nothing",
    )
    args = parser.parse_args(argv)
    repo_dir = Path(args.repo_dir).resolve()
    try:
        if args.local:
            return run_local(repo_dir, args.base, Path(args.body_file) if args.body_file else None)
        if not args.repo:
            raise AdmissionError("--repo or GITHUB_REPOSITORY is required")
        gh = GitHub(args.repo, os.environ.get("GITHUB_TOKEN"))
        if args.event:
            event = json.loads(Path(args.event).read_text())
            return run_event(event, gh, Git(repo_dir), repo_dir)
        if args.pull_request is not None:
            return run_dispatch(args.pull_request, gh, Git(repo_dir), repo_dir)
        if args.issue:
            return run_issue(json.loads(Path(args.issue).read_text()), gh, repo_dir)
        event = json.loads(Path(args.workflow_run).read_text())
        if args.resolve_pull_request:
            return resolve_pull_request(event, gh)
        return run_workflow_run(event, gh, repo_dir)
    except AdmissionError as exc:
        print(f"ERROR: [admission] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
