#!/usr/bin/env python3
"""Hold a pull request to the admission policy and close one that misses it.

usage: admission.py --event <github event json>          the pull_request_target job
       admission.py --pull-request <number>              the same checks, run by hand (workflow_dispatch)
       admission.py --workflow-run <github event json>   the job that runs when the CI workflow completes
       admission.py --local [--base origin/main] [--body-file <f>]
                                                        the static checks on the local branch before pushing

The static checks: the branch contains the tip of the base branch; a change to code under `crates/`,
`tools/` or `scripts/` comes with a new or changed test (a file under a crate's or tool's `tests/`
directory, a line of a `#[cfg(test)]` module or a `#[test]` function of a source file, a `test_*.py` or
`test_*.sh` under `scripts/`, a golden `expected/` or `expected_*/` file or a corpus `manifest.json`),
unless a maintainer sets the `tests-not-needed` label; a commit declares the model that wrote the change with a
`Co-Authored-By` trailer (`No LLM` for none); the description carries the template's three sections; a
crate change adds a line under `## [Unreleased]` of CHANGELOG.md, or under a release heading the change
itself adds; the parity matrix resolves on the contributed tree; and every commit carries a
`Signed-off-by` trailer with its author's address, a miss the comment reports without closing, since the
only remedy on a pushed branch is a history rewrite after which a closed pull request cannot be reopened,
and the DCO check holds the merge anyway. The job that runs after CI adds one line per required check
that failed.

The two Actions modes read the pull request through the REST API (commits, changed files, labels,
comments) and the contribution's files as git objects: a `git fetch` of the pull request's head ref, `git
show` for a file's text, and a `checkout-index` into a temporary directory for the parity check, whose
checker is the one of the checked-out base branch. Nothing of the contribution is checked out into the
workspace or executed. Every string a contributor controls (a path, a commit id, a matrix cell) is
rendered inside a code span, so it carries no mention, link or formatting.

One comment per pull request, edited in place on every run, states each miss with its fix. A miss closes
the pull request unless its author is a maintainer (a CODEOWNERS entry for `*`) or a bot, the
`admission-override` label is set, the pull request is a draft, or it was opened before this workflow
first landed on the base branch; the author reopens it once it is fixed and the check runs again. Exit 1
only when the run closes the pull request, 0 otherwise, 2 on a usage or environment error.
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
LABEL_MAX = 50
SECTIONS = ("Summary", "Verification", "Upgrade impact")
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
RUST_CODE_RE = re.compile(r"^(crates|tools)/.*\.rs$")
SCRIPT_CODE_RE = re.compile(r"^scripts/.*\.(py|sh)$")
CRATE_MANIFEST_RE = re.compile(r"^crates/[^/]+/Cargo\.toml$")
MOD_DECL_TEMPLATE = r"^\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+{name}\s*;"
HUNK_RE = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")
SIGNOFF_RE = re.compile(r"^Signed-off-by:\s*(.*?)\s*<([^<>]+)>\s*$", re.MULTILINE | re.IGNORECASE)
DECLARATION_RE = re.compile(r"^Co-Authored-By:\s*(.+?)\s*$", re.MULTILINE | re.IGNORECASE)
HEADING_RE = re.compile(r"^#{1,6}\s+(.+?)\s*#*\s*$")
BLOCK_RE = re.compile(
    r"<!-- admission:(static|ci) sha=([0-9a-f]+) closes=([01]) -->\n(.*?)\n<!-- /admission:\1 -->",
    re.DOTALL,
)
# `cfg(test)`, `cfg(all(test, ...))` and the like mark test-only items; `cfg(not(test))` does not.
CFG_TEST_RE = re.compile(r"^cfg\s*\((?!\s*not\s*\(\s*test\b).*\btest\b", re.DOTALL)
TEST_ATTR_RE = re.compile(
    r"^(?:[A-Za-z_][A-Za-z0-9_]*::)*(?:test|rstest|test_case|quickcheck|proptest)\b"
)
STRING_PREFIXES = ("r", "br", "cr", "b", "c")
# The words that open an item after an attribute; an attribute followed by a comma before one of them
# sits on a struct field or an enum variant, not on an item.
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
        "pub",
    )
)
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


def rust_test_regions(text: str) -> list[tuple[int, int]]:
    """Line ranges (1-based, inclusive) of the items a `#[cfg(test)]` or a test attribute (`#[test]`,
    `#[tokio::test]`, `#[rstest]`) applies to: the attribute line through the item's closing brace or
    semicolon. Braces are matched outside strings, character literals and comments; an inner
    `#![cfg(test)]` makes the whole file a test region."""
    n = len(text)
    total = len(text.splitlines())
    regions: list[tuple[int, int]] = []
    i = 0
    line = 1
    depth = 0
    armed: tuple[int, int] | None = None
    item_started = False
    parens = 0
    region_start: int | None = None
    region_depth = 0

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
            if armed is not None and region_start is None and parens == 0 and word in ITEM_KEYWORDS:
                item_started = True
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
            if inner and CFG_TEST_RE.match(attr):
                return [(1, total)]
            if (
                not inner
                and armed is None
                and region_start is None
                and (CFG_TEST_RE.match(attr) or TEST_ATTR_RE.match(attr))
            ):
                armed = (line, depth)
                item_started = False
            line += text.count("\n", i, j)
            i = j
        elif ch == "{":
            depth += 1
            if armed is not None and region_start is None and armed[1] == depth - 1:
                region_start = armed[0]
                region_depth = depth
            i += 1
        elif ch == "}":
            if region_start is not None and depth == region_depth:
                regions.append((region_start, line))
                region_start = None
                armed = None
            depth = max(depth - 1, 0)
            if armed is not None and depth < armed[1]:
                armed = None
            i += 1
        elif ch == ";":
            if armed is not None and region_start is None and armed[1] == depth and parens == 0:
                regions.append((armed[0], line))
                armed = None
            i += 1
        elif ch == ",":
            if (
                armed is not None
                and region_start is None
                and armed[1] == depth
                and parens == 0
                and not item_started
            ):
                armed = None
            i += 1
        elif ch in "([":
            parens += 1
            i += 1
        elif ch in ")]":
            parens = max(parens - 1, 0)
            i += 1
        else:
            i += 1
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
    return any(pattern.search(path) for pattern in TEST_PATH_RES)


def file_changed(entry: dict) -> bool:
    return entry.get("status") in ("added", "removed") or bool(
        (entry.get("additions") or 0) + (entry.get("deletions") or 0)
    )


RegionsFor = Callable[[str, str], list[tuple[int, int]]]


def classify_file(entry: dict, regions_for: RegionsFor) -> tuple[bool, bool]:
    """(changes code, changes a test) for one entry of the pull request's file list. `regions_for(side,
    path)` gives the test regions of the file's text on the `old` or `new` side."""
    path = entry.get("filename", "")
    previous = entry.get("previous_filename") or path
    if not file_changed(entry):
        return False, False
    if is_test_path(path):
        return False, True
    moved_test = previous != path and is_test_path(previous)
    if RUST_CODE_RE.match(path):
        patch = entry.get("patch")
        if patch is None:
            return True, moved_test
        removed, added = changed_lines(patch)
        new_regions = regions_for("new", path) if added else []
        old_regions = regions_for("old", previous) if removed else []
        in_test = [in_regions(n, new_regions) for n, _ in added] + [
            in_regions(n, old_regions) for n, _ in removed
        ]
        if not in_test:
            return True, moved_test
        return not all(in_test), any(in_test) or moved_test
    if SCRIPT_CODE_RE.match(path):
        return True, moved_test
    return False, moved_test


def check_tests_included(
    files: list[dict], labels: set[str], regions_for: RegionsFor
) -> Finding | None:
    if NO_TESTS_LABEL in labels:
        return None
    changes_code = changes_test = False
    for entry in files:
        is_code, is_test = classify_file(entry, regions_for)
        changes_code = changes_code or is_code
        changes_test = changes_test or is_test
    if changes_code and not changes_test:
        return Finding(
            "tests",
            "Code under `crates/`, `tools/` or `scripts/` changes and no test does: add or update a test "
            "for the change (a file under a `tests/` directory, a `#[test]` function or a `#[cfg(test)]` "
            "module, a `test_*.py` or `test_*.sh`, or the golden expectation the change moves), or a "
            f"maintainer adds the `{NO_TESTS_LABEL}` label.",
        )
    return None


# --- the commits ---------------------------------------------------------------------------------


def check_contains_base(contains: bool, base_ref: str) -> Finding | None:
    if contains:
        return None
    return Finding(
        "base",
        f"The branch does not contain the tip of `{base_ref}`: bring it in with "
        f"`git fetch origin && git merge --signoff origin/{base_ref}`, then push.",
    )


def check_signoff(commits: list[dict]) -> Finding | None:
    """Merge commits are not judged, as the DCO check leaves them out. The miss never closes: the
    remedy on a pushed branch rewrites history, after which a closed pull request cannot be reopened."""
    unsigned = []
    for commit in commits:
        if len(commit.get("parents") or []) > 1:
            continue
        message = commit["commit"].get("message") or ""
        author = ((commit["commit"].get("author") or {}).get("email") or "").lower()
        signers = {match.group(2).strip().lower() for match in SIGNOFF_RE.finditer(message)}
        if not author or author not in signers:
            unsigned.append(commit["sha"][:7])
    if not unsigned:
        return None
    listed = ", ".join(code(sha) for sha in unsigned[:5])
    more = f" and {len(unsigned) - 5} more" if len(unsigned) > 5 else ""
    return Finding(
        "signoff",
        f"Commit {listed}{more} carries no `Signed-off-by` trailer with its author's address: before "
        "the first push, sign every commit with `git rebase --signoff origin/main`; on a pushed branch, "
        "open a new pull request from a signed branch, since a force-push would stop this one from "
        "being reopened.",
        closes=False,
    )


def declared_models(commits: list[dict]) -> list[str]:
    """The distinct names of the `Co-Authored-By` trailers, their addresses dropped, `No LLM` normalised."""
    names: list[str] = []
    for commit in commits:
        for match in DECLARATION_RE.finditer(commit["commit"].get("message") or ""):
            name = re.sub(r"\s+", " ", re.sub(r"<[^>]*>", "", match.group(1))).strip()
            if name.lower() == "no llm":
                name = "No LLM"
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


def model_label(name: str) -> str | None:
    """`model: <name>` with the name reduced to letters, digits, spaces and `._+/-`, at most
    `LABEL_MAX` characters; None when nothing of the name survives."""
    clean = re.sub(r"\s+", " ", LABEL_CHARS_RE.sub("", name)).strip()
    if not clean:
        return None
    return (MODEL_LABEL_PREFIX + clean)[:LABEL_MAX].rstrip()


# --- the description and the changelog -----------------------------------------------------------


def check_body_sections(body: str | None) -> Finding | None:
    """Every template section is present as a heading with at least one line of text under it;
    the template's HTML comments do not count as text."""
    text = re.sub(r"<!--.*?-->", "", body or "", flags=re.DOTALL)
    sections: dict[str, list[str]] = {}
    current = None
    for raw in text.splitlines():
        heading = HEADING_RE.match(raw.strip())
        if heading:
            current = heading.group(1).strip().lower()
            sections.setdefault(current, [])
        elif current is not None:
            sections[current].append(raw)
    missing = [
        name
        for name in SECTIONS
        if not any(line.strip() for line in sections.get(name.lower(), []))
    ]
    if not missing:
        return None
    listed = ", ".join(f"`## {name}`" for name in missing)
    return Finding(
        "body",
        f"The description lacks the {listed} section of the pull request template, or leaves it "
        "empty: edit the description and fill it in.",
    )


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
    nothing from it; a crash of the checker is an environment error, not a miss."""
    try:
        proc = subprocess.run(
            [sys.executable, str(checker), "--root", str(tree)],
            capture_output=True,
            text=True,
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
        raise AdmissionError(
            f"the parity checker failed without a finding: {proc.stderr.strip()[-500:]}"
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


def close_footer(
    findings: list[Finding],
    user: dict,
    labels: set[str],
    draft: bool,
    maintainers: set[str],
    grandfathered: bool = False,
) -> tuple[bool, str | None]:
    """(close it, the sentence that ends the comment). Findings on an exempt pull request, and a
    sign-off miss alone, are stated with the reason the check does not close."""
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
        return False, (
            "No automatic close: the sign-off miss is reported only, and the DCO check holds the merge "
            "until every commit is signed off."
        )
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


def static_block(findings: list[Finding], head_sha: str, policy_url: str) -> str:
    head = code(head_sha[:7])
    if not findings:
        return f"The admission check ([policy]({policy_url})) passes on head {head}."
    lines = [f"The admission check ([policy]({policy_url})) finds on head {head}:"]
    lines.extend(f"- {finding.text}" for finding in findings)
    return "\n".join(lines)


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


def find_comment(gh: GitHub, number: int) -> dict | None:
    for comment in gh.paginate(f"/repos/{gh.repo}/issues/{number}/comments"):
        if (comment.get("body") or "").startswith(MARKER) and is_bot(comment.get("user") or {}):
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


def landed_at(gh: GitHub, base_ref: str) -> str | None:
    """The committer date of the first commit on the base branch that carries this workflow file."""
    commits = gh.paginate(f"/repos/{gh.repo}/commits", {"path": WORKFLOW_FILE, "sha": base_ref})
    if not commits:
        return None
    # The API lists commits newest first, so the last one is the oldest.
    return ((commits[-1].get("commit") or {}).get("committer") or {}).get("date")


def predates_workflow(gh: GitHub, pr: dict) -> bool:
    """Whether the pull request predates the workflow's landing; a failed read of the landing date
    counts as yes, so an API fault never closes a pull request."""
    try:
        return predates(pr.get("created_at"), landed_at(gh, pr["base"]["ref"]))
    except AdmissionError as exc:
        print(f"WARNING: [admission] the landing date was not read: {exc}", file=sys.stderr)
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
) -> list[Finding]:
    checks = [
        check_contains_base(contains, base_ref),
        check_tests_included(files, labels, regions_for),
        check_signoff(commits),
        check_declaration(models),
        check_body_sections(body),
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


def publish_findings(
    gh: GitHub,
    pr: dict,
    repo_dir: Path,
    kind: str,
    block: str,
    findings: list[Finding],
) -> bool:
    """Write this run's block into the pull request's comment, close the pull request when the
    findings call for it and it is open, and return whether this run closed it. The footer is the
    close or exemption sentence for the findings on this head: this run's, else the other block's
    when an earlier run on the same head left some. A clean run writes only to refresh a comment
    that exists."""
    number = pr["number"]
    head_sha = pr["head"]["sha"]
    codeowners = repo_dir / CODEOWNERS
    maintainers = maintainers_from_codeowners(
        codeowners.read_text() if codeowners.is_file() else ""
    )
    user = pr.get("user") or {}
    labels = label_names(pr)
    draft = bool(pr.get("draft"))
    older = predates_workflow(gh, pr)
    close, footer = close_footer(findings, user, labels, draft, maintainers, older)
    existing = find_comment(gh, number)
    blocks = parse_blocks(existing.get("body") if existing else None)
    blocks[kind] = (head_sha, block, any(f.closes for f in findings))
    other = blocks.get("ci" if kind == "static" else "static")
    if not footer and block_has_findings(other) and other[0] == head_sha:
        _close, footer = close_footer(
            [Finding("other", "", closes=other[2])], user, labels, draft, maintainers, older
        )
    if existing is not None or findings:
        write_comment(gh, number, existing, render_comment(blocks, head_sha, footer))
    if close and pr.get("state") == "open":
        gh.patch(f"/repos/{gh.repo}/pulls/{number}", {"state": "closed"})
        return True
    return False


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
    commits = gh.paginate(f"/repos/{gh.repo}/pulls/{number}/commits")
    files = gh.paginate(f"/repos/{gh.repo}/pulls/{number}/files")
    models = declared_models(commits)
    labels = label_names(pr)
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "tree"
        git.export_tree(head_sha, tree)
        parity = check_parity(tree, repo_dir / PARITY_CHECKER)
    findings = static_findings(
        contains,
        base_ref,
        commits,
        files,
        labels,
        pr.get("body"),
        git_regions_for(git, head_sha, merge_base),
        git.show(head_sha, "CHANGELOG.md"),
        parity,
        models,
    )
    block = static_block(findings, head_sha, policy_url(gh.repo, base_ref))
    close = publish_findings(gh, pr, repo_dir, "static", block, findings)
    sync_model_labels(gh, number, labels, models)
    print(f"admission: {len(commits)} commits, {len(files)} files, head {head_sha[:7]}")
    return verdict(findings, close)


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
    close = publish_findings(gh, pr, repo_dir, "ci", ci_block(findings, head_sha), findings)
    print(f"admission: CI run {run['id']} judged on head {head_sha[:7]}")
    return verdict(findings, close)


def local_commits(git: Git, base: str) -> list[dict]:
    out = git.run("log", "--format=%H%x00%P%x00%an%x00%ae%x00%B%x1e", f"{base}..HEAD")
    commits = []
    for record in out.split("\x1e"):
        if not record.strip():
            continue
        sha, parents, name, email, message = record.lstrip("\n").split("\x00", 4)
        commits.append(
            {
                "sha": sha,
                "parents": [{"sha": p} for p in parents.split()],
                "commit": {
                    "message": message,
                    "author": {"name": name, "email": email},
                },
            }
        )
    return commits


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
        raise AdmissionError(f"cannot resolve HEAD or {base}; run `git fetch origin` first")
    merge_base = git.merge_base(base_sha, head) or base_sha
    commits = local_commits(git, base)
    files = local_files(git, merge_base)
    models = declared_models(commits)
    filled_template = "\n".join(f"## {name}\nx" for name in SECTIONS)
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
    )
    report(findings)
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
    mode.add_argument("--local", action="store_true", help="check the local branch before pushing")
    parser.add_argument("--repo-dir", default=".", help="the base repository checkout (default: .)")
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"), help="owner/name")
    parser.add_argument(
        "--base",
        default="origin/main",
        help="--local: the base ref (default: origin/main)",
    )
    parser.add_argument("--body-file", help="--local: the pull request description to check")
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
        event = json.loads(Path(args.workflow_run).read_text())
        return run_workflow_run(event, gh, repo_dir)
    except AdmissionError as exc:
        print(f"ERROR: [admission] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
