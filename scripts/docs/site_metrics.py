#!/usr/bin/env python3
"""Measure the legibility of the documentation pages and gate on the limits.

The pages are the `(path.md)` links of `docs/SUMMARY.md`, or the paths given on
the command line. Each page is measured as the site renders it: a chapter whose
content is a repository file (`changelog.md`, `contributing.md`) is read from
that file, and a `<!-- readme-section: <heading prefix> -->` line is replaced by
the README section whose `## ` heading starts with the prefix, as the
preprocessor does. Per page: prose words (whitespace-separated tokens outside
table rows and fenced code, as `wc -w` counts them), headings per 300 prose
words, paragraphs over 120 words, sentences over 45 words, tables over 4
columns, and lines over 80 characters inside fenced code blocks. Headings,
tables and fenced code are not prose; each list item and each callout is a
paragraph of its own. A sentence that a fenced block interrupts is counted
once: when the paragraph before the block ends without a sentence terminator
and the paragraph after it opens with a continuing word (`and`, `or`, `but`,
`with`, ...), the two parts are one sentence, counted at the first part's line.

With `--gate` the exit status is 1 when any page breaks a limit it is not exempt
from: a page over 4000 prose words, a paragraph over 120 words, a sentence over
45 words, a table over 4 columns, or a code line over 80 characters outside a
block that follows an `N.2`, `N.3` or `N.6` heading (a divergence class's
record, engine output and reproducible commands, which are verbatim and cannot
wrap). `EXEMPT` lists every exemption with its reason; its `README.md` entry
applies to the lines a readme-section directive pulls in. Warnings never change the
exit status: heading density outside 0.6 to 1.5 per 300 words (pages of at
least 300 words), more than a quarter of the sentences over 25 words, a
paragraph over 5 sentences, more than 2 callouts on a page, two callouts in a
row.

Usage::

    python3 scripts/docs/site_metrics.py [--docs docs] [--gate] [--json FILE] [PAGE ...]

Prints a Markdown table (page, prose words, paragraphs over 120 words,
sentences over 45 words, tables over 4 columns, code lines over 80 characters,
warnings), then every violation with its line number, then the verdict. `--json`
writes the same per page, plus the total word count. Stdlib only.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdbook_repo_links import REPO_FILE_CHAPTERS

# repository path of a page -> rule -> why the page is exempt from it; the rules
# are words (the page cap), paragraphs, sentences, tables and code (the per-item
# caps below). The "README.md" entry covers the lines a readme-section directive
# pulls into a page.
EXEMPT = {
    "README.md": {
        "paragraphs": "the README sections pages include are written by the release renderer and kept byte-identical",
        "sentences": "the README sections pages include are written by the release renderer and kept byte-identical",
    },
    "docs/changelog.md": {
        "words": "append-only release history",
        "paragraphs": "released entries are never edited",
        "sentences": "released entries are never edited",
    },
    "docs/concordance-provenance/116.2-divergence-classes.md": {
        "words": "one catalogue of the comparator's report classes",
    },
    "docs/published-figures.md": {
        "paragraphs": "the paper's text stays as published",
        "sentences": "the paper's text stays as published",
        "tables": "the paper's tables stay as published",
    },
    "docs/intended-divergences.md": {
        "tables": "the defect and class tables, kept whole with a sticky first column",
        "code": "the output column-header line and two Ensembl VEP error messages, verbatim",
    },
}

WORD_CAP = 4000
PARAGRAPH_CAP = 120
SENTENCE_CAP = 45
SENTENCE_SOFT_CAP = 25
LONG_SENTENCE_SHARE = 0.25
TABLE_CAP = 4
CODE_LINE_CAP = 80
HEADING_DENSITY = (0.6, 1.5)
HEADING_DENSITY_MIN_WORDS = 300
PARAGRAPH_SENTENCE_CAP = 5
CALLOUTS_PER_PAGE = 2

SUMMARY_LINK = re.compile(r"\]\(([^)#\s]+\.md)\)")
DIRECTIVE = re.compile(r"^\s*<!--\s*readme-section:\s*(.*?)\s*-->\s*$")
HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
FENCE = re.compile(r"^\s*(`{3,}|~{3,})(.*)$")
LIST_ITEM = re.compile(r"^\s*(?:[-*+]|\d+[.)])\s+")
QUOTE_MARKS = re.compile(r"^\s*(?:>\s?)+")
CALLOUT_TAG = re.compile(r"^\[![A-Za-z]+\]\s*")
HTML_LINE = re.compile(r"^\s*(?:<[A-Za-z/!][^>]*>\s*)+$")
INCLUDE = re.compile(r"^\s*\{\{#\w+[^}]*\}\}\s*$")
THEMATIC_BREAK = re.compile(r"^\s*(?:(?:-\s*){3,}|(?:\*\s*){3,}|(?:_\s*){3,})$")
CODE_SPAN = re.compile(r"(`+)(.+?)\1")
# The record (N.2), engine output (N.3) and commands (N.6) parts of a divergence class.
VERBATIM_HEADING = re.compile(r"^\d+\.[236]\b")
SENTENCE_END = re.compile(r"[.!?]+[\"'”’)\]*_`]*(?:\s+|$)")
OPEN_END = re.compile(r"[.!?][\"'”’)\]*_`]*$")
# Lowercase words that open a paragraph only to carry on a sentence; "vep-rs" is
# lowercase too but opens sentences of its own, so a lowercase start is not enough.
CONTINUATION = re.compile(r"^(?:and|or|but|with|then|while|whereas|where|which|so|plus|followed)\b")
ABBREVIATIONS = {"e.g.", "i.e.", "vs.", "cf.", "etc.", "fig.", "no.", "approx."}
CLOSERS = "\"'”’()[]*_`"


def summary_pages(summary_path):
    with open(summary_path, encoding="utf-8") as fh:
        text = fh.read()
    pages = []
    for match in SUMMARY_LINK.finditer(text):
        page = match.group(1)
        if page not in pages:
            pages.append(page)
    return pages


def readme_section(readme_lines, prefix):
    """Return (start, end) such that readme_lines[start:end] is the README section the prefix names."""
    headings, in_fence = [], False
    for i, line in enumerate(readme_lines):
        if FENCE.match(line):
            in_fence = not in_fence
        elif not in_fence and line.startswith("## "):
            headings.append(i)
    starts = [i for i in headings if readme_lines[i][3:].lstrip().startswith(prefix)]
    if len(starts) != 1:
        sys.exit(
            f"ERROR: site_metrics: readme-section {prefix!r} matches {len(starts)} "
            "README headings, not one"
        )
    following = [i for i in headings if i > starts[0]]
    return starts[0], following[0] if following else len(readme_lines)


def load_page(docs_dir, page):
    """Return (source path, lines, origins) of a page as the site renders it.

    origins[i] is the line number of lines[i] in the source, or `README.md:N`
    for a line a readme-section directive pulled in.
    """
    root = os.path.dirname(os.path.abspath(docs_dir))
    source = page if os.path.isabs(page) else os.path.join(docs_dir, page)
    if page in REPO_FILE_CHAPTERS and os.path.exists(os.path.join(root, REPO_FILE_CHAPTERS[page])):
        source = os.path.join(root, REPO_FILE_CHAPTERS[page])
    if not os.path.exists(source):
        sys.exit(f"ERROR: site_metrics: {source} does not exist")
    with open(source, encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    out, origins = [], []
    readme = None
    for lineno, line in enumerate(lines, 1):
        directive = DIRECTIVE.match(line)
        if not directive:
            out.append(line)
            origins.append(lineno)
            continue
        if readme is None:
            with open(os.path.join(root, "README.md"), encoding="utf-8") as fh:
                readme = fh.read().split("\n")
        prefix = directive.group(1).split("|")[0].strip()
        start, end = readme_section(readme, prefix)
        out.extend(readme[start:end])
        origins.extend(f"README.md:{i + 1}" for i in range(start, end))
    return os.path.relpath(source), out, origins


def mask_code_spans(text):
    """Replace sentence punctuation inside inline code so it never ends a sentence."""

    def mask(match):
        return re.sub(r"[.!?]", "_", match.group(0))

    return CODE_SPAN.sub(mask, text)


def split_sentences(text):
    sentences = []
    start = 0
    for match in SENTENCE_END.finditer(text):
        last_token = text[start : match.end()].split()[-1]
        if last_token.lower().strip(CLOSERS) in ABBREVIATIONS:
            continue
        sentences.append(text[start : match.end()].strip())
        start = match.end()
    tail = text[start:].strip()
    if tail:
        sentences.append(tail)
    return sentences


def table_columns(row):
    cells = re.split(r"(?<!\\)\|", row.strip())
    if cells and cells[0].strip() == "":
        cells = cells[1:]
    if cells and cells[-1].strip() == "":
        cells = cells[:-1]
    return len(cells)


class Measurement:
    def __init__(self):
        self.headings = []
        self.paragraphs = []
        self.tables = []
        self.code_lines = []
        self.callouts = []
        self.blocks = []
        self.prose_words = 0
        self._para = None

    @property
    def kind(self):
        return self._para["kind"] if self._para else None

    def begin(self, lineno, text, kind, continues=None):
        """Open a paragraph; `continues` is the index of the paragraph whose open last
        sentence this one carries on past a fenced block."""
        self.flush()
        self._para = {"line": lineno, "lines": [text], "kind": kind, "continues": continues}

    def extend(self, text):
        self._para["lines"].append(text)

    def open_sentence(self):
        """Index of the paragraph holding the last paragraph's unfinished sentence, or None."""
        if not self.paragraphs or not self.blocks or self.blocks[-1][0] != "prose":
            return None
        last = self.paragraphs[-1]
        return last["tail_owner"] if last["open_end"] else None

    def flush(self):
        para = self._para
        self._para = None
        if para is None or not any(s.strip() for s in para["lines"]):
            return
        text = " ".join(s.strip() for s in para["lines"] if s.strip())
        masked = mask_code_spans(text)
        sentences = [len(s.split()) for s in split_sentences(masked)]
        owner = len(self.paragraphs)
        if para["continues"] is not None and sentences:
            self.paragraphs[para["continues"]]["sentences"][-1] += sentences.pop(0)
            if not sentences:
                owner = self.paragraphs[para["continues"]]["tail_owner"]
        self.paragraphs.append(
            {
                "line": para["line"],
                "words": len(masked.split()),
                "sentences": sentences,
                "text": text[:60],
                "open_end": not OPEN_END.search(masked.rstrip()),
                "tail_owner": owner,
            }
        )
        self.blocks.append((para["kind"], para["line"]))


def measure(lines):
    m = Measurement()
    fence = None
    in_verbatim = False
    last_heading_verbatim = False
    table_open = False
    interrupted = None
    for lineno, raw in enumerate(lines, 1):
        if fence is not None:
            closing = raw.strip()
            if closing.startswith(fence) and set(closing) == {fence[0]}:
                fence = None
                continue
            if len(raw) > CODE_LINE_CAP:
                m.code_lines.append({"line": lineno, "length": len(raw), "exempt": in_verbatim})
            continue
        opened = FENCE.match(raw)
        if opened:
            m.flush()
            table_open = False
            fence = opened.group(1)
            in_verbatim = last_heading_verbatim
            if not (m.blocks and m.blocks[-1][0] == "code"):
                interrupted = m.open_sentence()
            m.blocks.append(("code", lineno))
            continue
        if not raw.strip():
            m.flush()
            table_open = False
            continue
        if not (m.blocks and m.blocks[-1][0] == "code") or m.kind is not None:
            interrupted = None
        heading = HEADING.match(raw)
        if heading:
            m.flush()
            table_open = False
            text = heading.group(2)
            m.headings.append({"line": lineno, "level": len(heading.group(1)), "text": text})
            last_heading_verbatim = bool(VERBATIM_HEADING.match(text))
            m.blocks.append(("heading", lineno))
            continue
        if raw.lstrip().startswith("|"):
            m.flush()
            if not table_open:
                table_open = True
                m.tables.append({"line": lineno, "columns": table_columns(raw)})
                m.blocks.append(("table", lineno))
            continue
        table_open = False
        if HTML_LINE.match(raw) or INCLUDE.match(raw) or THEMATIC_BREAK.match(raw):
            m.flush()
            interrupted = None
            continue
        m.prose_words += len(raw.split())
        text = raw
        kind = "prose"
        quoted = QUOTE_MARKS.match(raw)
        if quoted:
            kind = "quote"
            text = raw[quoted.end() :]
            if not text.strip():
                m.flush()
                continue
            tagged = CALLOUT_TAG.match(text)
            if tagged:
                m.flush()
                m.callouts.append(lineno)
                m.blocks.append(("callout", lineno))
                text = text[tagged.end() :]
        item = LIST_ITEM.match(text)
        if item:
            m.begin(lineno, text[item.end() :], kind)
        elif m.kind != kind:
            carries_on = interrupted is not None and kind == "prose" and CONTINUATION.match(text.lstrip())
            m.begin(lineno, text, kind, interrupted if carries_on else None)
        else:
            m.extend(text)
        interrupted = None
    m.flush()
    return m


def where(page, line):
    return f"{page}:{line}" if isinstance(line, int) else f"{page} ({line})"


def cell(count, rule, exempt):
    return f"{count} (exempt)" if rule in exempt and count else str(count)


def report_page(page, source, lines, origins, exempt):
    m = measure(lines)
    for record in m.headings + m.paragraphs + m.tables + m.code_lines:
        record["line"] = origins[record["line"] - 1]
    m.callouts = [origins[line - 1] for line in m.callouts]
    m.blocks = [(kind, origins[line - 1]) for kind, line in m.blocks]
    words = m.prose_words
    words_total = sum(len(line.split()) for line in lines)
    sentences = [n for p in m.paragraphs for n in p["sentences"]]
    long_paragraphs = [p for p in m.paragraphs if p["words"] > PARAGRAPH_CAP]
    long_sentences = []
    for p in m.paragraphs:
        for n in p["sentences"]:
            if n > SENTENCE_CAP:
                long_sentences.append({"line": p["line"], "words": n, "text": p["text"]})
    wide_tables = [t for t in m.tables if t["columns"] > TABLE_CAP]
    long_code = [c for c in m.code_lines if not c["exempt"]]
    exempt_code = [c for c in m.code_lines if c["exempt"]]

    included = EXEMPT.get("README.md", {})

    def exempt_from(rule, line):
        return rule in exempt or (rule in included and isinstance(line, str))

    violations = []
    if words > WORD_CAP and "words" not in exempt:
        violations.append(f"{page}: {words} prose words, over {WORD_CAP}")
    for p in long_paragraphs:
        if not exempt_from("paragraphs", p["line"]):
            violations.append(f'{where(page, p["line"])} paragraph of {p["words"]} words: "{p["text"]}"')
    for s in long_sentences:
        if not exempt_from("sentences", s["line"]):
            violations.append(f'{where(page, s["line"])} sentence of {s["words"]} words in "{s["text"]}"')
    if "tables" not in exempt:
        for t in wide_tables:
            violations.append(f'{where(page, t["line"])} table of {t["columns"]} columns')
    if "code" not in exempt:
        for c in long_code:
            violations.append(f'{where(page, c["line"])} code line of {c["length"]} characters')

    warnings = []
    density = len(m.headings) * 300 / words if words else 0.0
    if words >= HEADING_DENSITY_MIN_WORDS and not (HEADING_DENSITY[0] <= density <= HEADING_DENSITY[1]):
        warnings.append(f"headings per 300 words {density:.1f}")
    over_soft = sum(1 for n in sentences if n > SENTENCE_SOFT_CAP)
    share = over_soft / len(sentences) if sentences else 0.0
    if share > LONG_SENTENCE_SHARE:
        warnings.append(f"sentences over {SENTENCE_SOFT_CAP} words {share:.0%} ({over_soft} of {len(sentences)})")
    long_runs = [p for p in m.paragraphs if len(p["sentences"]) > PARAGRAPH_SENTENCE_CAP]
    if long_runs:
        lines_text = ", ".join(f'L{p["line"]}' for p in long_runs)
        warnings.append(f"paragraphs over {PARAGRAPH_SENTENCE_CAP} sentences: {lines_text}")
    if len(m.callouts) > CALLOUTS_PER_PAGE:
        warnings.append(f"{len(m.callouts)} callouts")
    callout_blocks = [(i, line) for i, (kind, line) in enumerate(m.blocks) if kind == "callout"]
    for k in range(1, len(callout_blocks)):
        (i, line_a), (j, line_b) = callout_blocks[k - 1], callout_blocks[k]
        if all(kind == "quote" for kind, _ in m.blocks[i + 1 : j]):
            warnings.append(f"consecutive callouts L{line_a} and L{line_b}")

    return {
        "page": page,
        "source": source,
        "words": words,
        "words_total": words_total,
        "headings": len(m.headings),
        "headings_per_300_words": round(density, 2),
        "paragraphs": len(m.paragraphs),
        "sentences": len(sentences),
        "sentences_over_25": over_soft,
        "paragraphs_over_120": [{k: p[k] for k in ("line", "words", "text")} for p in long_paragraphs],
        "sentences_over_45": long_sentences,
        "included_over_caps": {
            "paragraphs": sum(1 for p in long_paragraphs if isinstance(p["line"], str)),
            "sentences": sum(1 for s in long_sentences if isinstance(s["line"], str)),
        },
        "wide_tables": wide_tables,
        "long_code_lines": long_code,
        "exempt_code_lines": exempt_code,
        "exempt": exempt,
        "violations": violations,
        "warnings": warnings,
    }


def print_table(reports):
    print("| Page | Prose words | Paragraphs over 120 words | Sentences over 45 words | Tables over 4 columns | Code lines over 80 chars | Warnings |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | --- |")
    for r in reports:
        exempt = r["exempt"]
        words = f'{r["words"]:,}' + (" (exempt)" if "words" in exempt and r["words"] > WORD_CAP else "")
        code = cell(len(r["long_code_lines"]), "code", exempt)
        if r["exempt_code_lines"]:
            code += f' (+{len(r["exempt_code_lines"])} in verbatim blocks)'
        included = r["included_over_caps"]
        paragraphs = cell(len(r["paragraphs_over_120"]) - included["paragraphs"], "paragraphs", exempt)
        sentences = cell(len(r["sentences_over_45"]) - included["sentences"], "sentences", exempt)
        if included["paragraphs"]:
            paragraphs += f' (+{included["paragraphs"]} in README sections)'
        if included["sentences"]:
            sentences += f' (+{included["sentences"]} in README sections)'
        print(
            f'| {r["page"]} | {words} | {paragraphs} | {sentences} '
            f'| {cell(len(r["wide_tables"]), "tables", exempt)} | {code} | {"; ".join(r["warnings"])} |'
        )


def main():
    ap = argparse.ArgumentParser(description="Measure the legibility of the documentation pages.")
    ap.add_argument("pages", nargs="*", help="pages to measure; default every page docs/SUMMARY.md lists")
    ap.add_argument("--docs", default="docs", help="the book's source directory")
    ap.add_argument("--gate", action="store_true", help="exit 1 on any violation")
    ap.add_argument("--json", metavar="FILE", help="write the per-page measurements as JSON")
    args = ap.parse_args()

    docs_dir = args.docs
    if args.pages:
        pages = []
        for path in args.pages:
            rel = os.path.relpath(os.path.abspath(path), os.path.abspath(docs_dir))
            pages.append(rel.replace(os.sep, "/") if not rel.startswith("..") else path)
    else:
        pages = summary_pages(os.path.join(docs_dir, "SUMMARY.md"))

    reports = []
    for page in pages:
        source, lines, origins = load_page(docs_dir, page)
        reports.append(report_page(page, source, lines, origins, EXEMPT.get("docs/" + page, {})))

    print_table(reports)
    violations = [v for r in reports for v in r["violations"]]
    if violations:
        print()
        for v in violations:
            print(f"- {v}")
    total_sentences = sum(r["sentences"] for r in reports)
    total_soft = sum(r["sentences_over_25"] for r in reports)
    share = f"{total_soft / total_sentences:.0%}" if total_sentences else "0%"
    print()
    print(f"Sentences over {SENTENCE_SOFT_CAP} words: {total_soft} of {total_sentences} ({share}).")
    failing = sum(1 for r in reports if r["violations"])
    verdict = "PASS" if not violations else f"FAIL ({len(violations)} violations on {failing} of {len(reports)} pages)"
    print(f"Gate: {verdict}")

    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump({"verdict": verdict, "violations": len(violations), "pages": reports}, fh, indent=1)
            fh.write("\n")
    if args.gate and violations:
        sys.exit(1)


if __name__ == "__main__":
    main()
