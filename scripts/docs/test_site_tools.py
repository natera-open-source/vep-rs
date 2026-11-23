"""Tests for the site tooling in this directory, on fixtures the tests write: the link check, the mdBook preprocessor driven
through its stdin protocol, the page metrics and their gate, the post-build finish with its two checks, the
home page's writer, the screenshot script's pure helpers, the site build's refusal of a builder warning on a
stand-in for the builder, and the theme script's version selector run under node on a stand-in for its DOM calls.
No test runs mdbook or a browser."""

from __future__ import annotations

import io
import json
import os
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import check_site_links as link_check
import mdbook_repo_links as links
import render_home_page as home
import render_release_figure as figure
import site_build
import site_metrics as metrics
import site_postbuild as postbuild
import site_screenshots as shots
import test_render_release_figure as figure_records

REPO_URL = "https://github.com/natera-open-source/vep-rs"
EDIT_URL = f"{REPO_URL}/edit/main/"
RELEASE = "v1.2.3"
BASE = "https://docs.example/vep-rs/"
STUB = "# vep-rs documentation\n\nThe home page is the README.\n"
OVERVIEW = "# Overview\n\nWelcome.\n"

# The two fixture pages sit in a subdirectory of the book, so link resolution across
# directories is exercised.
SUB = "guide"
GUIDE_PAGE = f"{SUB}/page.md"
GUIDE_SETUP = f"{SUB}/setup.md"

# One README for the preprocessor and the metrics: a section with links and a fenced `## ` line, a
# section that a shorter prefix also matches, and a last section that runs to the end of the file.
README = f"""# Project

Intro.

## Install

Download the [archive](dist/archive.tar.gz), read [Notes](#notes) and [`docs/{SUB}/setup.md`](docs/{SUB}/setup.md).

```text
## Not a heading
```

## Notes

A note.

## Notes on output

Columns.
"""
CHANGELOG = f"""# Changelog

## [1.0.0] - first release

First.

[1.0.0]: {REPO_URL}/releases/tag/v1.0.0
"""


@pytest.fixture(autouse=True)
def _fresh_repository_files():
    """The preprocessor caches each repository file it reads by path, for the life of the process."""
    links.read_repo_file.cache_clear()


def git(root, *args):
    subprocess.run(
        ["git", "-C", str(root), "-c", "commit.gpgsign=false", "-c", "user.name=Tester",
         "-c", "user.email=tester@docs.example", *args],
        check=True,
        capture_output=True,
    )


def repository(tmp_path, readme=README):
    root = tmp_path / "repo"
    (root / "docs" / SUB).mkdir(parents=True)
    (root / "scripts").mkdir()
    (root / "README.md").write_text(readme, encoding="utf-8")
    (root / "CHANGELOG.md").write_text(CHANGELOG, encoding="utf-8")
    (root / "CONTRIBUTING.md").write_text("# Contributing\n\nOpen a pull request.\n", encoding="utf-8")
    (root / "scripts" / "README.md").write_text("# Scripts\n", encoding="utf-8")
    return root


def chapter(source, content):
    return {"Chapter": {"name": source, "content": content, "source_path": source, "sub_items": []}}


def site_book(page):
    """A book whose `GUIDE_PAGE` chapter holds `page`, beside a home page, a sibling and the two
    chapters that mirror repository files."""
    return [
        chapter("index.md", OVERVIEW),
        chapter(GUIDE_PAGE, page),
        chapter(GUIDE_SETUP, "# Setup\n\n## Step 2\n"),
        chapter("changelog.md", ""),
        chapter("contributing.md", ""),
    ]


def context(root, repo_url=REPO_URL):
    html = {"git-repository-url": repo_url} if repo_url else {}
    return {"root": str(root), "config": {"book": {"src": "docs"}, "output": {"html": html}}}


def set_variable(monkeypatch, name, value):
    """Set the environment variable, or remove it when `value` is None."""
    if value is None:
        monkeypatch.delenv(name, raising=False)
    else:
        monkeypatch.setenv(name, value)


def run_preprocessor(monkeypatch, root, items, env=None, repo_url=REPO_URL):
    settings = {"VEP_RS_DOCS_REF": None, "VEP_RS_DOCS_RELEASE": RELEASE, "VEP_RS_DOCS_README_HOME": None}
    settings.update(env or {})
    for name, value in settings.items():
        set_variable(monkeypatch, name, value)
    monkeypatch.setattr(sys, "stdin", io.StringIO(json.dumps([context(root, repo_url), {"items": items}])))
    monkeypatch.setattr(sys, "argv", ["mdbook_repo_links.py"])
    links.main()


def preprocess(monkeypatch, capsys, root, items, env=None):
    """The chapters' content by source path after the preprocessor has rewritten the book."""
    run_preprocessor(monkeypatch, root, items, env)
    book = json.loads(capsys.readouterr().out)
    return {c["source_path"]: c["content"] for c in links.walk(book["items"])}


def preprocess_error(monkeypatch, root, items, env=None, repo_url=REPO_URL):
    with pytest.raises(SystemExit) as excinfo:
        run_preprocessor(monkeypatch, root, items, env, repo_url)
    return excinfo.value.code


def test_supports_exits_zero(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["mdbook_repo_links.py", "supports", "html"])
    with pytest.raises(SystemExit) as excinfo:
        links.main()
    assert excinfo.value.code == 0


@pytest.mark.parametrize(
    ("setting", "content", "home"),
    [
        ("1", STUB, True),
        ("1", OVERVIEW, True),
        ("0", STUB, False),
        (None, STUB, True),
        (None, OVERVIEW, False),
        (None, STUB + "\n", False),
        (None, STUB.replace("\n", "\r\n"), False),
        (None, None, False),
    ],
    ids=["set-stub", "set-landing", "unset-by-zero", "stub", "landing", "stub-trailing-blank", "crlf-stub", "no-index"],
)
def test_readme_is_home_follows_the_variable_else_the_three_line_stub(monkeypatch, setting, content, home):
    set_variable(monkeypatch, "VEP_RS_DOCS_README_HOME", setting)
    items = [chapter(GUIDE_PAGE, "Text.\n")]
    if content is not None:
        items.append(chapter("index.md", content))
    assert links.readme_is_home({"items": items}) is home


def test_stub_home_page_shows_the_readme_with_its_links_rewritten(monkeypatch, capsys, tmp_path):
    readme = "# Project\n\nSee [the scripts](scripts/README.md) and [the changelog](CHANGELOG.md).\n"
    root = repository(tmp_path, readme)
    pages = preprocess(monkeypatch, capsys, root, [chapter("index.md", STUB), chapter("changelog.md", "")])
    assert pages["index.md"] == (
        f"# Project\n\nSee [the scripts]({REPO_URL}/blob/main/scripts/README.md) and [the changelog](changelog.md).\n"
    )


def test_home_variable_zero_keeps_the_stub_and_sends_its_readme_link_to_github(monkeypatch, capsys, tmp_path):
    stub = "# vep-rs documentation\n\nThe home page is the [README](../README.md).\n"
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), [chapter("index.md", stub)],
                       env={"VEP_RS_DOCS_README_HOME": "0"})
    assert pages["index.md"] == f"# vep-rs documentation\n\nThe home page is the [README]({REPO_URL}/blob/main/README.md).\n"


def test_directive_pastes_the_one_matching_section_up_to_the_next_heading(monkeypatch, capsys, tmp_path):
    page = "Before.\n\n<!-- readme-section: Install -->\n\nAfter.\n"
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book(page))
    assert pages[GUIDE_PAGE] == (
        "Before.\n\n## Install\n\n"
        f"Download the [archive]({REPO_URL}/blob/main/dist/archive.tar.gz), read [Notes](#notes) and [`setup.md`](setup.md).\n\n"
        "```text\n## Not a heading\n```\n\n\nAfter.\n"
    )


def test_directive_for_the_last_section_runs_to_the_end_of_the_readme(monkeypatch, capsys, tmp_path):
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book("<!-- readme-section: Notes on -->\n"))
    assert pages[GUIDE_PAGE] == "## Notes on output\n\nColumns.\n"


def test_noheading_option_drops_the_heading_line(monkeypatch, capsys, tmp_path):
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book("<!-- readme-section: Notes on | noheading -->\n"))
    assert pages[GUIDE_PAGE] == "\nColumns.\n"


def test_fragments_option_points_bare_fragment_links_at_the_named_chapter(monkeypatch, capsys, tmp_path):
    page = "<!-- readme-section: Install | noheading | fragments: index.md -->\n"
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book(page))
    assert pages[GUIDE_PAGE] == (
        f"\nDownload the [archive]({REPO_URL}/blob/main/dist/archive.tar.gz), read [Notes](../index.md#notes) and "
        "[`setup.md`](setup.md).\n\n```text\n## Not a heading\n```\n\n"
    )


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ("Notes", "readme-section 'Notes' must match one `## ` heading of README.md; matched: ## Notes, ## Notes on output"),
        ("Not a", "readme-section 'Not a' must match one `## ` heading of README.md; matched: none"),
        ("Install | fragments: nowhere.md", "readme-section fragments: 'nowhere.md' is not a chapter"),
        ("Install | bogus", "readme-section option 'bogus' is not `noheading` or `fragments: <chapter>`"),
    ],
    ids=["two-matches", "fenced-heading-is-no-match", "fragments-not-a-chapter", "unknown-option"],
)
def test_directive_errors_name_the_cause(monkeypatch, tmp_path, spec, message):
    page = f"<!-- readme-section: {spec} -->\n"
    assert preprocess_error(monkeypatch, repository(tmp_path), site_book(page)) == f"mdbook_repo_links: {message}"


def test_directive_failure_exits_one_with_the_message_on_stderr(tmp_path):
    root = repository(tmp_path)
    env = {k: v for k, v in os.environ.items() if not k.startswith("VEP_RS_DOCS_")}
    env["VEP_RS_DOCS_RELEASE"] = RELEASE
    out = subprocess.run(
        [sys.executable, str(HERE / "mdbook_repo_links.py")],
        input=json.dumps([context(root), {"items": site_book("<!-- readme-section: Not a -->\n")}]),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert out.returncode == 1
    assert out.stdout == ""
    assert out.stderr == "mdbook_repo_links: readme-section 'Not a' must match one `## ` heading of README.md; matched: none\n"


def test_tokens_take_the_ref_and_release_from_the_environment(monkeypatch, capsys, tmp_path):
    page = "Ref {{vep-rs-ref}}, release {{vep-rs-release}}, [scripts](../../scripts/README.md).\n"
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book(page),
                       env={"VEP_RS_DOCS_REF": "v2.0.0", "VEP_RS_DOCS_RELEASE": "v2.0.0"})
    assert pages[GUIDE_PAGE] == f"Ref v2.0.0, release v2.0.0, [scripts]({REPO_URL}/blob/v2.0.0/scripts/README.md).\n"


def test_tokens_default_to_main_and_a_link_to_the_latest_release(monkeypatch, capsys, tmp_path):
    root = repository(tmp_path)
    git(root, "init", "-q")
    pages = preprocess(monkeypatch, capsys, root, site_book("{{vep-rs-ref}} {{vep-rs-release}}\n"),
                       env={"VEP_RS_DOCS_RELEASE": None})
    assert pages[GUIDE_PAGE] == f"main [on GitHub]({REPO_URL}/releases/latest)\n"


def test_release_token_defaults_to_the_newest_tag_reachable_from_head(monkeypatch, capsys, tmp_path):
    root = repository(tmp_path)
    git(root, "init", "-q")
    git(root, "commit", "-q", "--allow-empty", "-m", "first")
    git(root, "tag", "v9.9.9")
    pages = preprocess(monkeypatch, capsys, root, site_book("{{vep-rs-ref}} {{vep-rs-release}}\n"),
                       env={"VEP_RS_DOCS_RELEASE": None})
    assert pages[GUIDE_PAGE] == "main v9.9.9\n"


def test_links_directives_and_tokens_inside_fenced_code_are_untouched(monkeypatch, capsys, tmp_path):
    page = "```md\n[scripts](../../scripts/README.md) {{vep-rs-ref}}\n<!-- readme-section: Install -->\n```\n"
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book(page))
    assert pages[GUIDE_PAGE] == page


@pytest.mark.parametrize(
    ("target", "expected"),
    [
        ("setup.md", "setup.md"),
        ("setup.md#step-2", "setup.md#step-2"),
        ("../index.md", "../index.md"),
        ("../../CHANGELOG.md", "../changelog.md"),
        ("../../CONTRIBUTING.md", "../contributing.md"),
        ("../../README.md", f"{REPO_URL}/blob/main/README.md"),
        ("../../scripts/README.md", f"{REPO_URL}/blob/main/scripts/README.md"),
        ("../../scripts", f"{REPO_URL}/tree/main/scripts"),
        ("missing.md", f"{REPO_URL}/blob/main/docs/{SUB}/missing.md"),
        ("https://example.org/page", "https://example.org/page"),
        ("#top", "#top"),
        ("../../../outside.md", "../../../outside.md"),
    ],
    ids=["page", "page-fragment", "home", "changelog", "contributing", "readme-not-home", "repository-file",
         "repository-directory", "page-not-in-book", "absolute-url", "bare-fragment", "above-the-repository"],
)
def test_link_targets_are_rewritten_for_the_site(monkeypatch, capsys, tmp_path, target, expected):
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book(f"See [it]({target}).\n"))
    assert pages[GUIDE_PAGE] == f"See [it]({expected}).\n"


@pytest.mark.parametrize(
    ("target", "expected"),
    [
        ("../assets/figure.png", "../assets/figure.png"),
        ("figure.png", "figure.png"),
        ("../assets/figure.png#top", "../assets/figure.png#top"),
        ("../../manuscript/figure.png", f"{REPO_URL}/blob/main/manuscript/figure.png"),
        ("https://example.org/figure.png", "https://example.org/figure.png"),
    ],
    ids=["book-file", "beside-the-chapter", "book-file-fragment", "outside-the-book", "absolute-url"],
)
def test_image_under_the_book_source_stays_relative_to_its_chapter(monkeypatch, capsys, tmp_path, target, expected):
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book(f"![Figure]({target})\n"))
    assert pages[GUIDE_PAGE] == f"![Figure]({expected})\n"


def test_link_to_a_book_file_that_is_no_chapter_opens_its_github_view_beside_an_image_of_it(monkeypatch, capsys, tmp_path):
    page = 'See [the figure](../assets/figure.png) and ![Figure](../assets/figure.png "A title").\n'
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book(page))
    assert pages[GUIDE_PAGE] == (
        f'See [the figure]({REPO_URL}/blob/main/docs/assets/figure.png) and ![Figure](../assets/figure.png "A title").\n'
    )


def test_image_of_a_readme_section_points_at_the_site_copy_from_every_page_that_shows_it(monkeypatch, capsys, tmp_path):
    readme = (
        "# Project\n\n## Results\n\n![The figure](docs/assets/figure.png)\n\n"
        "A [linked ![badge](docs/assets/badge.svg)](docs/assets/figure.png) image.\n"
    )
    book = [chapter("index.md", STUB), chapter(GUIDE_PAGE, "<!-- readme-section: Results -->\n")]
    pages = preprocess(monkeypatch, capsys, repository(tmp_path, readme), book)
    view = f"{REPO_URL}/blob/main/docs/assets/figure.png"
    assert pages[GUIDE_PAGE] == (
        f"## Results\n\n![The figure](../assets/figure.png)\n\nA [linked ![badge](../assets/badge.svg)]({view}) image.\n"
    )
    assert pages["index.md"] == (
        f"# Project\n\n## Results\n\n![The figure](assets/figure.png)\n\nA [linked ![badge](assets/badge.svg)]({view}) image.\n"
    )


def test_readme_link_points_at_the_home_page_when_the_readme_is_home(monkeypatch, capsys, tmp_path):
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book("See [it](../../README.md).\n"),
                       env={"VEP_RS_DOCS_README_HOME": "1"})
    assert pages[GUIDE_PAGE] == "See [it](../index.md).\n"


def test_changelog_chapter_shows_the_repository_changelog_with_its_heading_link_moved(monkeypatch, capsys, tmp_path):
    pages = preprocess(monkeypatch, capsys, repository(tmp_path), site_book("Text.\n"))
    assert pages["changelog.md"] == (
        "# Changelog\n\n## 1.0.0 - first release\n\n"
        f"[Release v1.0.0]({REPO_URL}/releases/tag/v1.0.0)\n\n\nFirst.\n\n[1.0.0]: {REPO_URL}/releases/tag/v1.0.0\n"
    )


def test_missing_repository_url_is_refused(monkeypatch, tmp_path):
    code = preprocess_error(monkeypatch, repository(tmp_path), site_book("Text.\n"), repo_url=None)
    assert code == "mdbook_repo_links: output.html.git-repository-url is not set in book.toml"


def test_missing_repository_file_behind_a_chapter_is_refused(monkeypatch, tmp_path):
    root = repository(tmp_path)
    (root / "CHANGELOG.md").unlink()
    assert preprocess_error(monkeypatch, root, site_book("Text.\n")) == f"mdbook_repo_links: CHANGELOG.md is missing from {root}"


def test_a_repository_file_that_is_not_utf8_is_refused_in_one_line(monkeypatch, tmp_path):
    root = repository(tmp_path)
    readable = b"# Project\n\n## Install\n\n"
    (root / "README.md").write_bytes(readable + b"\xff\n")
    code = preprocess_error(monkeypatch, root, site_book("<!-- readme-section: Install -->\n"))
    assert code == f"mdbook_repo_links: README.md in {root} is not UTF-8: byte {len(readable)} does not decode"


LONG_PARAGRAPH = " ".join(["One two three four five."] * 25)
LONG_SENTENCE = " ".join(["word"] * 45) + " end."
PAGE = f"""# Title

Short opening paragraph here.

{LONG_PARAGRAPH}

{LONG_SENTENCE}

| A | B | C | D | E |
| --- | --- | --- | --- | --- |
| 1 | 2 | 3 | 4 | 5 |

```text
{"x" * 81}
```
"""


def docs_tree(tmp_path, pages):
    docs = tmp_path / "repo" / "docs"
    docs.mkdir(parents=True)
    (docs.parent / "README.md").write_text(README, encoding="utf-8")
    for name, content in pages.items():
        (docs / name).write_text(content, encoding="utf-8")
    (docs / "SUMMARY.md").write_text("# Summary\n\n" + "".join(f"- [{n}]({n})\n" for n in pages), encoding="utf-8")
    return docs


def measure_site(docs, *flags):
    return subprocess.run(
        [sys.executable, str(HERE / "site_metrics.py"), "--docs", str(docs), *flags],
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize(("prefix", "span"), [("Install", (4, 12)), ("Notes on", (16, 20))])
def test_readme_section_spans_its_heading_to_the_next_heading_or_the_end(prefix, span):
    assert metrics.readme_section(README.split("\n"), prefix) == span


@pytest.mark.parametrize(("prefix", "count"), [("Notes", 2), ("Not a", 0)])
def test_readme_section_needs_exactly_one_matching_heading(prefix, count):
    with pytest.raises(SystemExit) as excinfo:
        metrics.readme_section(README.split("\n"), prefix)
    assert excinfo.value.code == f"ERROR: site_metrics: readme-section {prefix!r} matches {count} README headings, not one"


def test_page_measures_find_each_kind_of_excess(tmp_path):
    docs = docs_tree(tmp_path, {"page.md": PAGE})
    source, lines, origins = metrics.load_page(str(docs), "page.md")
    report = metrics.report_page("page.md", source, lines, origins, {})
    assert report["words"] == 175
    assert report["paragraphs_over_120"] == [{"line": 5, "words": 125, "text": LONG_PARAGRAPH[:60]}]
    assert report["sentences_over_45"] == [{"line": 7, "words": 46, "text": LONG_SENTENCE[:60]}]
    assert report["wide_tables"] == [{"line": 9, "columns": 5}]
    assert report["long_code_lines"] == [{"line": 14, "length": 81, "exempt": False}]
    assert report["violations"] == [
        f'page.md:5 paragraph of 125 words: "{LONG_PARAGRAPH[:60]}"',
        f'page.md:7 sentence of 46 words in "{LONG_SENTENCE[:60]}"',
        "page.md:9 table of 5 columns",
        "page.md:14 code line of 81 characters",
    ]


def test_gate_passes_a_tree_within_every_limit(tmp_path):
    docs = docs_tree(tmp_path, {"index.md": "# Home\n\nWelcome.\n\n<!-- readme-section: Notes on -->\n"})
    out = measure_site(docs, "--gate")
    assert out.returncode == 0, out.stderr
    assert "| index.md | 2 | 0 | 0 | 0 | 0 |  |\n" in out.stdout
    assert out.stdout.endswith("Gate: PASS\n")


def test_gate_fails_a_tree_with_a_violation_and_lists_it(tmp_path):
    docs = docs_tree(tmp_path, {"index.md": "# Home\n\nWelcome.\n", "page.md": PAGE})
    out = measure_site(docs, "--gate")
    assert out.returncode == 1
    assert "- page.md:9 table of 5 columns\n" in out.stdout
    assert out.stdout.endswith("Gate: FAIL (4 violations on 1 of 2 pages)\n")


def test_without_the_gate_flag_a_violation_leaves_the_exit_status_zero(tmp_path):
    docs = docs_tree(tmp_path, {"page.md": PAGE})
    out = measure_site(docs)
    assert out.returncode == 0
    assert out.stdout.endswith("Gate: FAIL (4 violations on 1 of 1 pages)\n")


def test_the_index_of_intended_divergences_has_a_ceiling_of_its_own(tmp_path):
    """Every page is held to `WORD_CAP` prose words but the index of intended divergences, whose ceiling `WORD_CAPS` gives."""
    assert metrics.WORD_CAP == 4000
    assert metrics.WORD_CAPS == {"docs/intended-divergences.md": 4500}

    def page(words):
        """A page of `words` prose words: paragraphs of four 25-word sentences, then the rest as one sentence."""
        sentence = " ".join(["word"] * 24) + " word."
        paragraphs = [" ".join([sentence] * 4)] * (words // 100)
        if words % 100:
            paragraphs.append(" ".join(["word"] * (words % 100 - 1) + ["word."]))
        return "# Page\n\n" + "\n\n".join(paragraphs) + "\n"

    docs = docs_tree(tmp_path, {"intended-divergences.md": page(4500), "other.md": page(4000)})
    out = measure_site(docs, "--gate")
    assert out.returncode == 0, out.stdout
    assert "| intended-divergences.md | 4,500 |" in out.stdout and "| other.md | 4,000 |" in out.stdout
    (docs / "intended-divergences.md").write_text(page(4501), encoding="utf-8")
    (docs / "other.md").write_text(page(4001), encoding="utf-8")
    out = measure_site(docs, "--gate")
    assert out.returncode == 1
    assert "- intended-divergences.md: 4501 prose words, over 4500\n" in out.stdout
    assert "- other.md: 4001 prose words, over 4000\n" in out.stdout


VERSIONS = {
    "release": RELEASE,
    "channels": [{"name": RELEASE, "path": "/vep-rs/"}, {"name": "main", "path": "/vep-rs/main/"}],
}


def page_html(title, body, book=f"vep-rs {RELEASE}"):
    return (
        "<!DOCTYPE html>\n<html>\n<head>\n"
        f"<title>{title} - {book}</title>\n</head>\n<body>\n"
        f'<h1 class="menu-title">{book}</h1>\n{body}\n</body>\n</html>\n'
    )


def redirect_html(destination):
    return (
        "<!DOCTYPE html>\n<html>\n<head>\n"
        f'<meta http-equiv="refresh" content="0; URL={destination}">\n</head>\n<body></body>\n</html>\n'
    )


HOME_PROSE = "The home page of the fixture site opens with a paragraph long enough to stand as its description."
SETUP_PROSE = "The setup page of the fixture site opens with a paragraph of its own, as long as a description must be."
OTHER_PROSE = "A third page of the fixture site opens with a paragraph that no other page of the site begins with."


def built_site(tmp_path, redirect="guide/setup.html#step-2", book=f"vep-rs {RELEASE}"):
    """A built site of two reader pages, each with a title and an opening paragraph of its own, two redirect
    stubs, and the two pages mdBook generates for printing and not-found."""
    site = tmp_path / "book"
    (site / "guide").mkdir(parents=True)
    edit = f'<a href="{EDIT_URL}docs/index.md">Edit</a>\n<h1 id="overview">Overview</h1>\n<p>{HOME_PROSE}</p>'
    (site / "index.html").write_text(page_html("vep-rs", edit, book), encoding="utf-8")
    (site / "guide" / "setup.html").write_text(
        page_html("Setup", f'<h1 id="setup">Setup</h1>\n<p>{SETUP_PROSE}</p>\n<h2 id="step-2">Step 2</h2>', book),
        encoding="utf-8",
    )
    (site / "old.html").write_text(redirect_html(redirect), encoding="utf-8")
    moved = redirect_html("guide/setup.html").replace(
        "</head>",
        '<script>\nlet redirectUrl = "guide/setup.html";\n'
        'let fragmentMap = {"#old": "guide/setup.html#step-2", "#ext": "https://example.org/"};\n</script>\n</head>',
    )
    (site / "moved.html").write_text(moved, encoding="utf-8")
    (site / "print.html").write_text(page_html("Print", "<p>Every page.</p>"), encoding="utf-8")
    (site / "404.html").write_text(page_html("Not found", "<p>No such page.</p>"), encoding="utf-8")
    return site


def finish(monkeypatch, capsys, site, *extra, base=BASE, readme_home=None):
    """Exit status and output of the post-build step on `site`."""
    set_variable(monkeypatch, "VEP_RS_DOCS_README_HOME", readme_home)
    monkeypatch.setattr(sys, "argv", ["site_postbuild.py", str(site), "--base", base, "--release", RELEASE, *extra])
    with pytest.raises(SystemExit) as excinfo:
        postbuild.main()
    return excinfo.value.code, capsys.readouterr().out


def test_canonical_and_og_url_are_inserted_once_and_a_second_run_adds_nothing(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert index.count(f'<link rel="canonical" href="{BASE}">') == 1
    assert index.count(f'<meta property="og:url" content="{BASE}">') == 1
    setup = (site / "guide" / "setup.html").read_text(encoding="utf-8")
    assert setup.count(f'<link rel="canonical" href="{BASE}guide/setup.html">') == 1
    for page in ("old.html", "moved.html", "print.html", "404.html"):
        assert 'rel="canonical"' not in (site / page).read_text(encoding="utf-8")
    assert "6 pages finished: menu title demoted, head tags (canonical, og:url) inserted on 2, 0 already carried a canonical link\n" in out
    assert "redirect destinations: 3 checked, 1 external not fetched, 0 broken\n" in out
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert (site / "index.html").read_text(encoding="utf-8") == index
    assert "inserted on 0, 2 already carried a canonical link\n" in out


def test_menu_title_is_demoted_and_a_title_repeating_the_book_name_collapses(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    finish(monkeypatch, capsys, site)
    index = (site / "index.html").read_text(encoding="utf-8")
    assert postbuild.HOME_TAGLINE == "fast Ensembl VEP (Variant Effect Predictor) in Rust"
    assert f'<p class="menu-title">vep-rs {RELEASE}</p>' in index
    assert "<h1 class=" not in index
    assert f"<title>vep-rs: {postbuild.HOME_TAGLINE}</title>" in index
    assert "<title>Setup - vep-rs documentation</title>" in (site / "guide" / "setup.html").read_text(encoding="utf-8")
    assert "<title>Not found - vep-rs</title>" in (site / "404.html").read_text(encoding="utf-8")


def test_preview_channel_marks_its_pages_and_writes_no_sitemap(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f'<link rel="canonical" href="{BASE}main/">' in index
    assert '<meta name="vep-rs-channel" content="main">' in index
    assert '<meta name="robots" content="noindex">' in index
    assert not (site / "sitemap.xml").exists()
    assert json.loads((site / "versions.json").read_text(encoding="utf-8")) == VERSIONS


def test_release_form_adds_no_channel_metas_and_writes_versions_and_the_sitemap(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert "vep-rs-channel" not in index
    assert "noindex" not in index
    assert (site / "sitemap.xml").read_text(encoding="utf-8") == (
        '<?xml version="1.0" encoding="UTF-8"?>\n<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
        f"  <url><loc>{BASE}guide/setup.html</loc></url>\n  <url><loc>{BASE}</loc></url>\n</urlset>\n"
    )
    assert json.loads((site / "versions.json").read_text(encoding="utf-8")) == VERSIONS


BOOK_DESCRIPTION = "The book-wide description."
PROSE = (
    "vep-rs takes <a href=\"x.html\">Ensembl VEP&#39;s</a> inputs and flags and writes its consequence vocabulary "
    "in its output formats, two orders of magnitude faster than the Perl implementation on every dataset measured."
)


def described_page(title, body, book=f"vep-rs {RELEASE}"):
    metas = "".join(
        f'<meta {attr} content="{BOOK_DESCRIPTION}">\n'
        for attr in ('name="description"', 'property="og:description"', 'name="twitter:description"')
    )
    metas += f'<meta property="og:site_name" content="{book}">\n<meta property="og:title" content="vep-rs">\n'
    metas += '<meta name="twitter:title" content="vep-rs">\n'
    return page_html(title, body, book).replace("</head>", metas + "</head>")


def test_each_page_gets_its_first_paragraph_of_prose_as_its_description(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    (site / "index.html").write_text(
        described_page("vep-rs", f"<main>\n<p>Documentation: short.</p>\n<p>{PROSE}</p>\n</main>"), encoding="utf-8"
    )
    (site / "guide" / "setup.html").write_text(
        described_page("Setup", "<main>\n<h2 id=\"step-2\">Step 2</h2>\n<p>Too short to describe anything.</p>\n</main>"), encoding="utf-8"
    )
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    expected = postbuild.page_description(index)
    assert expected.startswith("vep-rs takes Ensembl VEP's inputs and flags")
    assert expected.endswith("\u2026") and len(expected) <= postbuild.DESCRIPTION_MAX
    assert BOOK_DESCRIPTION not in index
    assert index.count(postbuild.html.escape(expected)) == 3
    assert f'<meta property="og:title" content="vep-rs: {postbuild.HOME_TAGLINE}">' in index
    assert f'<meta name="twitter:title" content="vep-rs: {postbuild.HOME_TAGLINE}">' in index
    assert '<meta property="og:site_name" content="vep-rs">' in index
    record = json.loads(index.split('<script type="application/ld+json">')[1].split("</script>")[0])
    assert record["@type"] == "SoftwareApplication"
    assert record["softwareVersion"] == RELEASE.lstrip("v")
    assert record["description"] == expected
    setup = (site / "guide" / "setup.html").read_text(encoding="utf-8")
    assert setup.count(BOOK_DESCRIPTION) == 3
    assert '<meta property="og:site_name" content="vep-rs">' in setup
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert (site / "index.html").read_text(encoding="utf-8") == index


def test_preview_home_gets_the_tagline_but_no_json_ld(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    (site / "index.html").write_text(
        described_page("vep-rs", f"<main><p>{PROSE}</p></main>", "vep-rs main"), encoding="utf-8"
    )
    code, out = finish(monkeypatch, capsys, site, "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f"<title>vep-rs main: {postbuild.HOME_TAGLINE}</title>" in index
    assert '<meta property="og:site_name" content="vep-rs main">' in index
    assert "application/ld+json" not in index


def test_page_heads_line_counts_the_pages_a_reader_lands_on(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert "page heads: 2 pages, 0 titles shared, 0 descriptions shared, 0 on the book's description\n" in out
    assert "WARNING" not in out


def third_page_finding(monkeypatch, capsys, tmp_path, title, prose, counts, finding):
    """Holds the post-build step, on each channel, to its finding over the fixture site plus `guide/other.html`:
    a warning and exit status 0 on the release, a failure on the preview. `{name}` in `finding` is what a short
    title ends with on the channel."""
    for channel, name in ((None, "vep-rs documentation"), ("main", "vep-rs main")):
        book = name if channel else f"vep-rs {RELEASE}"
        site = built_site(tmp_path / (channel or "release"), book=book)
        (site / "guide" / "other.html").write_text(page_html(title, f"<p>{prose}</p>", book), encoding="utf-8")
        extra = ("--channel", channel) if channel else ()
        base = BASE + (f"{channel}/" if channel else "")
        code, out = finish(monkeypatch, capsys, site, *extra, base=base)
        assert f"page heads: 3 pages, {counts}\n" in out
        assert f"{'ERROR' if channel else 'WARNING'}: [site_postbuild] {finding.format(name=name)}\n" in out
        assert code == (1 if channel else 0)


def test_a_title_two_pages_share_fails_the_preview_and_is_a_warning_on_the_release(monkeypatch, capsys, tmp_path):
    third_page_finding(
        monkeypatch, capsys, tmp_path, "Setup", OTHER_PROSE,
        "1 titles shared, 0 descriptions shared, 0 on the book's description",
        "title 'Setup - {name}' is on 2 pages: guide/other.html, guide/setup.html",
    )


def test_a_description_two_pages_share_fails_the_preview_and_is_a_warning_on_the_release(monkeypatch, capsys, tmp_path):
    third_page_finding(
        monkeypatch, capsys, tmp_path, "Other", SETUP_PROSE,
        "0 titles shared, 1 descriptions shared, 0 on the book's description",
        f"description {SETUP_PROSE!r} is on 2 pages: guide/other.html, guide/setup.html",
    )


def test_a_page_left_on_the_book_description_fails_the_preview_and_is_a_warning_on_the_release(
    monkeypatch, capsys, tmp_path
):
    third_page_finding(
        monkeypatch, capsys, tmp_path, "Other", "Too short to describe anything.",
        "0 titles shared, 0 descriptions shared, 1 on the book's description",
        "guide/other.html: no paragraph of 80 characters or more, so its description is the book's",
    )


def test_a_sitemap_that_does_not_list_as_many_pages_as_a_reader_lands_on_fails(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    reader_pages = {"index.html", "guide/setup.html"}
    assert postbuild.check_page_heads(str(site), reader_pages, None) == []
    sitemap = site / "sitemap.xml"
    listed = sitemap.read_text(encoding="utf-8")
    home = f"  <url><loc>{BASE}</loc></url>\n"
    sitemap.write_text(listed.replace(home, ""), encoding="utf-8")
    assert postbuild.check_page_heads(str(site), reader_pages, None) == ["sitemap.xml lists 1 pages, a reader lands on 2"]
    sitemap.write_text(listed.replace(home, home + home), encoding="utf-8")
    assert postbuild.check_page_heads(str(site), reader_pages, None) == ["sitemap.xml lists 3 pages, a reader lands on 2"]
    sitemap.unlink()
    assert postbuild.check_page_heads(str(site), reader_pages, None) == ["sitemap.xml lists 0 pages, a reader lands on 2"]
    assert postbuild.check_page_heads(str(site), reader_pages, "main") == []


ENDING, ENDING_LONG = " - vep-rs", " - vep-rs documentation"
HOME_PAGE_TITLE = "vep-rs: fast Ensembl VEP (Variant Effect Predictor) in Rust"
CUT_A = "Co-located variants, frequencies and known structural variants"
CUT_B = "Co-located variants, frequencies and known structural rearrangements"
CUT_TITLE = "Co-located variants, frequencies and known structural\u2026"
TITLES_LINE = 'page titles: {} pages, {} end " - vep-rs", {} end " - vep-rs documentation", {} without an ending, {} cut, {} outside 25 to 60 characters\n'


@pytest.mark.parametrize(
    ("heading", "title", "form"),
    [
        ("Quick start", "Quick start" + ENDING_LONG, "lengthened"),
        ("x" * 15, "x" * 15 + ENDING_LONG, "lengthened"),
        ("x" * 16, "x" * 16 + ENDING, "kept"),
        ("x" * 51, "x" * 51 + ENDING, "kept"),
        ("x" * 52, "x" * 52, "bare"),
        ("x" * 60, "x" * 60, "bare"),
        ("x" * 61, "x" * 59 + "\u2026", "cut"),
    ],
    ids=["20-with-the-ending", "24-with-the-ending", "25-with-the-ending", "60-with-the-ending", "61-with-the-ending", "60-alone", "61-alone"],
)
def test_a_title_takes_the_form_its_length_calls_for(heading, title, form):
    assert postbuild.fitted_title(heading) == title
    assert postbuild.title_form(title) == form
    assert 25 <= len(title) <= 60


@pytest.mark.parametrize(
    ("heading", "title"),
    [
        (CUT_A, CUT_TITLE),
        ("The fields a flag switches on in each output format, and the fields it leaves", "The fields a flag switches on in each output format\u2026"),
        ("The warning file, the statistics report and the fields of a plugin's own file", "The warning file, the statistics report and the fields\u2026"),
        ("x" * 21 + " of " + "y" * 50, "x" * 21 + " of\u2026"),
        ("x" * 23 + ", " + "y" * 50, "x" * 23 + ",\u2026"),
        ("ab cd " + "x" * 70, "ab cd " + "x" * 53 + "\u2026"),
    ],
    ids=["last-word-that-fits", "comma-and-a-small-word-dropped", "two-small-words-dropped", "small-word-kept-at-25", "comma-kept-at-25", "no-space-past-24"],
)
def test_a_heading_over_60_characters_is_cut_at_a_word_and_ends_in_an_ellipsis(heading, title):
    assert len(heading) > 60
    assert postbuild.fitted_title(heading) == title
    assert postbuild.title_form(title) == "cut"
    assert 25 <= len(title) <= 60


def test_the_home_title_is_within_the_limits_and_the_fit_leaves_it():
    assert f"{postbuild.BOOK_NAME}: {postbuild.HOME_TAGLINE}" == HOME_PAGE_TITLE
    assert len(HOME_PAGE_TITLE) == 59 == len("vep-rs: ") + len(postbuild.HOME_TAGLINE)
    page = f"<head>\n<title>{HOME_PAGE_TITLE}</title>\n</head>\n"
    assert postbuild.fit_title(page) == page


def chapter_html(name, n, book=f"vep-rs {RELEASE}"):
    """A page as the theme writes it, the chapter name escaped before the book's in its title and whole in its
    `og:title`, opening with a paragraph that `n` makes its own."""
    escaped = postbuild.html.escape(name)
    prose = f"Page {n} of the fixture site opens with a paragraph of its own, as long as a description must be."
    meta = f'<meta property="og:title" content="{escaped}">\n'
    return page_html(escaped, f"<p>{prose}</p>", book).replace("</head>", meta + "</head>")


def add_chapters(site, names, book=f"vep-rs {RELEASE}"):
    """Writes `guide/page<n>.html` for each chapter name; the pages' paths under the site."""
    pages = [f"guide/page{n}.html" for n in range(len(names))]
    for n, name in enumerate(names):
        (site / pages[n]).write_text(chapter_html(name, n, book), encoding="utf-8")
    return pages


def title_of(site, page):
    return re.search(r"<title>(.*)</title>", (site / page).read_text(encoding="utf-8"))[1]


def test_the_release_fits_each_title_and_a_second_run_changes_no_byte(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    kept, bare, cut = add_chapters(site, ["Setting up a cache", "x" * 52, CUT_A])
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert {page: title_of(site, page) for page in ("index.html", "guide/setup.html", kept, bare, cut)} == {
        "index.html": HOME_PAGE_TITLE,
        "guide/setup.html": "Setup" + ENDING_LONG,
        kept: "Setting up a cache" + ENDING,
        bare: "x" * 52,
        cut: CUT_TITLE,
    }
    assert f'<meta property="og:title" content="{CUT_A}">' in (site / cut).read_text(encoding="utf-8")
    assert TITLES_LINE.format(5, 1, 1, 2, 1, 0) in out
    finished = {path: path.read_bytes() for path in site.rglob("*") if path.is_file()}
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert {path: path.read_bytes() for path in site.rglob("*") if path.is_file()} == finished
    assert TITLES_LINE.format(5, 1, 1, 2, 1, 0) in out


def test_the_preview_keeps_every_title_whole_and_prints_no_titles_line(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path, book="vep-rs main")
    names = ["A", "Setting up a cache", "x" * 52, CUT_A]
    pages = add_chapters(site, names, "vep-rs main")
    code, out = finish(monkeypatch, capsys, site, "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    assert "page titles:" not in out
    assert [title_of(site, page) for page in pages] == [f"{name} - vep-rs main" for name in names]
    assert title_of(site, "guide/setup.html") == "Setup - vep-rs main"


def test_a_title_is_measured_unescaped_and_written_back_escaped(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    kept, lengthened = add_chapters(site, ["R&D " + "x" * 47, "A <B> & C's"])
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert title_of(site, kept) == "R&amp;D " + "x" * 47 + ENDING
    assert title_of(site, lengthened) == "A &lt;B&gt; &amp; C's" + ENDING_LONG
    assert TITLES_LINE.format(4, 1, 2, 1, 0, 0) in out


def test_a_title_under_25_characters_fails_the_release_and_passes_the_preview(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path / "release")
    (page,) = add_chapters(site, ["A"])
    code, out = finish(monkeypatch, capsys, site)
    assert code == 1
    assert f"ERROR: [site_postbuild] {page}: title 'A - vep-rs documentation' is 24 characters, under 25\n" in out
    assert TITLES_LINE.format(3, 0, 2, 1, 0, 1) in out
    site = built_site(tmp_path / "main", book="vep-rs main")
    add_chapters(site, ["A"], "vep-rs main")
    code, out = finish(monkeypatch, capsys, site, "--channel", "main", base=BASE + "main/")
    assert code == 0, out


def test_a_title_over_60_characters_fails_the_release(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(postbuild, "HOME_TAGLINE", postbuild.HOME_TAGLINE + ", every flag")
    site = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site)
    assert code == 1
    assert f"ERROR: [site_postbuild] index.html: title '{HOME_PAGE_TITLE}, every flag' is 71 characters, over 60\n" in out
    assert TITLES_LINE.format(2, 0, 1, 1, 0, 1) in out


def test_two_headings_cut_to_one_title_fail_the_release_and_pass_the_preview(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path / "release")
    first, second = add_chapters(site, [CUT_A, CUT_B])
    code, out = finish(monkeypatch, capsys, site)
    assert code == 1
    assert "page heads: 4 pages, 1 titles shared, 0 descriptions shared, 0 on the book's description\n" in out
    assert f"ERROR: [site_postbuild] title {CUT_TITLE!r} is the fit of 2 headings: {first}, {second}\n" in out
    assert "WARNING" not in out
    site = built_site(tmp_path / "main", book="vep-rs main")
    add_chapters(site, [CUT_A, CUT_B], "vep-rs main")
    code, out = finish(monkeypatch, capsys, site, "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    assert "page heads: 4 pages, 0 titles shared, 0 descriptions shared, 0 on the book's description\n" in out


def test_one_heading_cut_on_two_pages_is_a_warning_on_the_release(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    first, second = add_chapters(site, [CUT_A, CUT_A])
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert f"WARNING: [site_postbuild] title {CUT_TITLE!r} is on 2 pages: {first}, {second}\n" in out
    assert "is the fit of" not in out


def test_a_chapter_name_that_ends_in_the_release_tag_keeps_it_on_a_second_run(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    kept, bare = add_chapters(site, [f"Upgrading to {RELEASE}", f"{'x' * 45} {RELEASE}"])
    for _run in range(2):
        code, out = finish(monkeypatch, capsys, site)
        assert code == 0, out
        assert title_of(site, kept) == f"Upgrading to {RELEASE}" + ENDING
        assert title_of(site, bare) == f"{'x' * 45} {RELEASE}"


def test_a_home_chapter_of_another_name_is_titled_with_the_book_name_alone(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    (site / "index.html").write_text(described_page("Overview", f"<p>{HOME_PROSE}</p>"), encoding="utf-8")
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f"<title>{HOME_PAGE_TITLE}</title>" in index
    assert f'<meta property="og:title" content="{HOME_PAGE_TITLE}">' in index


FRAME = '<noscript>\n<iframe class="sidebar-iframe-outer" src="toc.html"></iframe>\n</noscript>'
CHAPTERS = (
    '<ol class="chapter"><li><a href="index.html" target="_parent">vep-rs</a></li>'
    '<li><a href="guide/setup.html" target="_parent">Setup</a></li></ol>'
)


def test_home_page_lists_every_chapter_where_its_noscript_block_framed_the_list(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    for page in (site / "index.html", site / "guide" / "setup.html"):
        framed = page.read_text(encoding="utf-8").replace("<body>\n", f"<body>\n<nav>{FRAME}</nav>\n")
        page.write_text(framed, encoding="utf-8")
    (site / "toc.html").write_text(
        f'<!DOCTYPE HTML>\n<html>\n<head></head>\n<body class="sidebar-iframe-inner">\n{CHAPTERS}\n</body>\n</html>\n',
        encoding="utf-8",
    )
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f'<nav><noscript>\n<div class="sidebar-scrollbox">{CHAPTERS}</div>\n</noscript></nav>' in index
    assert "<iframe" not in index
    assert f"<nav>{FRAME}</nav>" in (site / "guide" / "setup.html").read_text(encoding="utf-8")
    assert "home page: 2 chapter links listed in its <noscript> block\n" in out
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert (site / "index.html").read_text(encoding="utf-8") == index


SUMMARY_CHAPTER = re.compile(r"\[(.+)\]\(([^()\s]+\.md)\)\s*$")


def test_every_chapter_of_the_summary_has_a_title_of_its_own():
    lines = (ROOT / "docs" / "SUMMARY.md").read_text(encoding="utf-8").splitlines()
    names = [m[1] for m in map(SUMMARY_CHAPTER.search, lines) if m]
    assert names
    assert sorted({name for name in names if names.count(name) > 1}) == []


def test_every_chapter_of_the_summary_keeps_a_title_of_its_own_within_the_limits_once_fitted():
    """mdBook drops the backticks of inline code in a chapter name, so a title is fitted from the name without them."""
    lines = (ROOT / "docs" / "SUMMARY.md").read_text(encoding="utf-8").splitlines()
    chapters = [(m[2], m[1].replace("`", "")) for m in map(SUMMARY_CHAPTER.search, lines) if m]
    titles = {}
    for path, name in chapters:
        title = HOME_PAGE_TITLE if path == "index.md" else postbuild.fitted_title(name)
        titles.setdefault(title, []).append(name)
    assert HOME_PAGE_TITLE in titles
    assert {title: names for title, names in titles.items() if len(names) > 1} == {}
    assert sorted(title for title in titles if not 25 <= len(title) <= 60) == []


def test_the_readme_opening_the_book_description_and_the_home_title_name_the_same_things():
    paragraphs = (ROOT / "README.md").read_text(encoding="utf-8").split("\n\n")
    opening = next(p for p in paragraphs if not p.startswith(("#", "[![")) and len(p) >= postbuild.DESCRIPTION_MIN)
    book = re.search(r'^description = "(.*)"$', (ROOT / "book.toml").read_text(encoding="utf-8"), re.MULTILINE)
    assert book, "book.toml has no description line"
    texts = {
        "the README's opening": opening[: postbuild.DESCRIPTION_MAX],
        "book.toml's description": book[1],
        "the home page's title": postbuild.HOME_TAGLINE,
    }
    for where, text in texts.items():
        for name in ("Rust", "Ensembl", "Variant Effect Predictor"):
            assert name in text, f"{where} does not name {name}: {text}"


def test_the_header_example_of_the_default_format_names_the_workspace_version():
    """The engine line `## vep-rs <version>` of the header example on `docs/output/vep-and-tab.md` names the version
    `[workspace.package]` of the root `Cargo.toml` declares, the version the binary writes on that line."""
    page = (ROOT / "docs" / "output" / "vep-and-tab.md").read_text(encoding="utf-8")
    workspace = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))["workspace"]
    assert re.findall(r"^## vep-rs (\S+)$", page, re.MULTILINE) == [workspace["package"]["version"]]


@pytest.mark.parametrize(
    ("destination", "message"),
    [
        ("gone.html", "old.html: redirect to gone.html: no such file in the site"),
        ("guide/setup.html#nowhere", "old.html: redirect to guide/setup.html#nowhere: no id 'nowhere' on guide/setup.html"),
        ("/elsewhere/page.html", "old.html: redirect to /elsewhere/page.html leaves the site"),
        ("/vep-rs-old/page.html", "old.html: redirect to /vep-rs-old/page.html leaves the site"),
        ("guide/", "old.html: redirect to guide/: no such file in the site"),
    ],
    ids=["missing-page", "missing-id", "outside-the-site", "longer-name-than-the-site-root", "directory-without-an-index-page"],
)
def test_a_redirect_without_its_destination_fails_by_name(monkeypatch, capsys, tmp_path, destination, message):
    site = built_site(tmp_path, redirect=destination)
    code, out = finish(monkeypatch, capsys, site)
    assert code == 1
    assert f"ERROR: [site_postbuild] {message}\n" in out


@pytest.mark.parametrize(
    "destination", ["/vep-rs", "/vep-rs/", "guide/"], ids=["site-root-without-the-slash", "site-root", "directory"]
)
def test_a_redirect_to_a_directory_is_checked_against_its_index_page(monkeypatch, capsys, tmp_path, destination):
    site = built_site(tmp_path, redirect=destination)
    (site / "guide" / "index.html").write_text(page_html("Guide", '<h1 id="guide">Guide</h1>'), encoding="utf-8")
    code, out = finish(monkeypatch, capsys, site)
    assert code == 0, out
    assert "redirect destinations: 3 checked, 1 external not fetched, 0 broken\n" in out


LINK_CHECK_BOOK = """[book]
title = "vep-rs"
src = "docs"

[output.html]
site-url = "/vep-rs/"
{table}
[output.html.search]
limit-results = 20
"""


def linked_site(tmp_path, entries, home_link="guide/setup.html#step-2"):
    """A checkout holding a built site of three reader pages, their sources beside the chapter list under `docs`, and
    the site's book configuration, whose redirect table holds `entries`; `None` writes a configuration without the
    table."""
    root = tmp_path / "checkout"
    site = root / "target" / "book"
    (site / "guide").mkdir(parents=True)
    (root / "docs" / "guide").mkdir(parents=True)
    for source in ("SUMMARY.md", "index.md", "guide/setup.md", "guide/page.md"):
        (root / "docs" / source).write_text("# A page\n", encoding="utf-8")
    home = f'<h1 id="overview">Overview</h1>\n<a href="{home_link}">Setup</a>\n<a href="https://example.org/">Elsewhere</a>'
    setup = '<h1 id="setup">Setup</h1>\n<h2 id="step-2">Step 2</h2>\n<a href="../index.html#overview">Home</a>'
    (site / "index.html").write_text(page_html("vep-rs", home), encoding="utf-8")
    (site / "guide" / "setup.html").write_text(page_html("Setup", setup), encoding="utf-8")
    (site / "guide" / "page.html").write_text(page_html("Page", '<h1 id="page">Page</h1>'), encoding="utf-8")
    (site / "print.html").write_text(page_html("Print", '<a href="gone.html">Gone</a>'), encoding="utf-8")
    table = "" if entries is None else "\n[output.html.redirect]\n# Pages and headings that moved.\n\n" + "".join(
        f"{entry}\n" for entry in entries
    )
    (root / "book.toml").write_text(LINK_CHECK_BOOK.format(table=table), encoding="utf-8")
    return root


def check_links(cwd, *flags):
    """The finished run of the link check started in `cwd`."""
    return subprocess.run(
        [sys.executable, str(HERE / "check_site_links.py"), *flags], cwd=cwd, capture_output=True, text=True, check=False
    )


LANDING_REDIRECTS = (
    '"/old.html" = "guide/setup.html#step-2"',
    '"/guide/page.html#was-here" = "setup.html#step-2"',
    '"/guide/page.html#up" = "../index.html#overview"',
    '"/index.html#top" = "#overview"',
    f'"/index.html#source" = "{REPO_URL}/blob/main/book.toml"',
    '"/index.html#readme" = "https://example.org/readme" # another host',
    '"/index.html#root" = "/vep-rs/guide/setup.html"',
)
NOT_FOLLOWED = "not followed (a URL or a path from the server root)"
ALL_BUILT = "pages under docs: 3 checked, 0 not in the site\n"
NOT_BUILT = "in the site: mdBook builds the chapters {}/SUMMARY.md lists\n"


def test_the_link_check_passes_a_site_whose_links_and_redirects_all_land(tmp_path):
    out = check_links(linked_site(tmp_path, LANDING_REDIRECTS))
    assert (out.returncode, out.stderr) == (0, "")
    assert out.stdout == (
        "checked 2 links on 3 pages, 0 broken\n"
        f"redirects of book.toml: 5 checked, 2 {NOT_FOLLOWED}, 0 broken\n"
        + ALL_BUILT
    )


def test_the_link_check_lists_each_redirect_whose_destination_the_site_or_the_repository_lacks(tmp_path):
    lost = (
        '"/old.html" = "guide/gone.html"',
        '"/guide/page.html#was-here" = "setup.html#step-9"',
        '"/guide/page.html#up" = "../index.html#overview"',
        '"/index.html#top" = "#summary"',
        f'"/index.html#source" = "{REPO_URL}/blob/main/docs/gone.md"',
    )
    out = check_links(linked_site(tmp_path, lost))
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "book.toml: /old.html -> guide/gone.html: no such file in the site\n"
        "book.toml: /guide/page.html#was-here -> setup.html#step-9: no id 'step-9' on guide/setup.html\n"
        "book.toml: /index.html#top -> #summary: no id 'summary' on index.html\n"
        f"book.toml: /index.html#source -> {REPO_URL}/blob/main/docs/gone.md: no such path in the repository\n"
        "checked 2 links on 3 pages, 0 broken\n"
        f"redirects of book.toml: 5 checked, 0 {NOT_FOLLOWED}, 4 broken\n"
        + ALL_BUILT
    )


def test_the_link_check_holds_a_redirect_to_a_directory_to_its_index_page(tmp_path):
    entries = ('"/old.html" = "guide/"', '"/guide/page.html#up" = "../#overview"', '"/index.html#top" = "guide/#step-2"')
    root = linked_site(tmp_path, entries)
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "book.toml: /old.html -> guide/: no such file in the site\n"
        "book.toml: /index.html#top -> guide/#step-2: no such file in the site\n"
        "checked 2 links on 3 pages, 0 broken\n"
        f"redirects of book.toml: 3 checked, 0 {NOT_FOLLOWED}, 2 broken\n"
        + ALL_BUILT
    )
    (root / "target" / "book" / "guide" / "index.html").write_text(page_html("Guide", '<h1 id="guide">Guide</h1>'), encoding="utf-8")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "book.toml: /index.html#top -> guide/#step-2: no id 'step-2' on guide/index.html\n"
        "checked 2 links on 4 pages, 0 broken\n"
        f"redirects of book.toml: 3 checked, 0 {NOT_FOLLOWED}, 1 broken\n"
        "pages under docs: 3 checked, 0 not in the site\n"
    )


def test_the_link_check_lists_a_broken_link_before_a_lost_redirect_and_counts_each_apart(tmp_path):
    root = linked_site(tmp_path, ['"/old.html" = "guide/setup.html#step-9"'], home_link="guide/setup.html#step-9")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "index.html: guide/setup.html#step-9: no id 'step-9' on guide/setup.html\n"
        "book.toml: /old.html -> guide/setup.html#step-9: no id 'step-9' on guide/setup.html\n"
        "checked 2 links on 3 pages, 1 broken\n"
        f"redirects of book.toml: 1 checked, 0 {NOT_FOLLOWED}, 1 broken\n"
        + ALL_BUILT
    )


def test_the_link_check_reads_the_book_its_option_names(tmp_path):
    root = linked_site(tmp_path, ['"/old.html" = "guide/setup.html"', '"/guide/page.html#was-here" = "setup.html#step-9"'])
    book = root / "book.toml"
    out = check_links(tmp_path, "--site", str(root / "target" / "book"), "--book", str(book))
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        f"{book}: /guide/page.html#was-here -> setup.html#step-9: no id 'step-9' on guide/setup.html\n"
        "checked 2 links on 3 pages, 0 broken\n"
        f"redirects of {book}: 2 checked, 0 {NOT_FOLLOWED}, 1 broken\n"
        f"pages under {root / 'docs'}: 3 checked, 0 not in the site\n"
    )


def test_the_link_check_passes_a_book_without_a_redirect_table(tmp_path):
    out = check_links(linked_site(tmp_path, None))
    assert (out.returncode, out.stderr) == (0, "")
    assert out.stdout == (
        f"checked 2 links on 3 pages, 0 broken\nredirects of book.toml: 0 checked, 0 {NOT_FOLLOWED}, 0 broken\n"
        + ALL_BUILT
    )


def test_the_link_check_refuses_a_site_whose_book_it_cannot_find(tmp_path):
    root = linked_site(tmp_path, LANDING_REDIRECTS)
    (root / "book.toml").unlink()
    out = check_links(root)
    assert (out.returncode, out.stdout) == (1, "")
    assert out.stderr == (
        "ERROR: [check_site_links] no book configuration at book.toml: name the one the site was built from with --book\n"
    )


@pytest.mark.parametrize(
    "entry",
    ["'/old.html' = 'guide/setup.html'", '"/old.html" = { to = "guide/setup.html" }', '"/old.html" ='],
    ids=["literal-strings", "inline-table", "no-destination"],
)
def test_the_link_check_refuses_a_redirect_entry_it_cannot_read(tmp_path, entry):
    out = check_links(linked_site(tmp_path, ['"/guide/page.html#up" = "../index.html"', entry]))
    assert (out.returncode, out.stdout) == (1, "")
    assert out.stderr == f"ERROR: [check_site_links] unrecognised redirect entry in book.toml: {entry}\n"


def test_the_link_check_reads_every_entry_of_the_committed_redirect_table():
    book = ROOT / "book.toml"
    sources = [line.split('"')[1] for line in book.read_text(encoding="utf-8").splitlines() if line.startswith('"/')]
    assert sources
    assert [source for source, _ in link_check.read_redirects(str(book))] == sources


def unlisted(root, *sources, docs="docs"):
    """Writes each of `sources` under the fixture's source directory as a page no chapter lists and the site lacks."""
    for source in sources:
        (root / docs / source).write_text("# A page no chapter lists\n", encoding="utf-8")


def test_the_link_check_lists_each_source_page_the_site_lacks_after_the_lost_redirects(tmp_path):
    root = linked_site(tmp_path, ['"/old.html" = "guide/setup.html#step-9"'])
    unlisted(root, "notes.md", "guide/extra.md")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "book.toml: /old.html -> guide/setup.html#step-9: no id 'step-9' on guide/setup.html\n"
        "docs/guide/extra.md: no page guide/extra.html " + NOT_BUILT.format("docs")
        + "docs/notes.md: no page notes.html " + NOT_BUILT.format("docs")
        + "checked 2 links on 3 pages, 0 broken\n"
        f"redirects of book.toml: 1 checked, 0 {NOT_FOLLOWED}, 1 broken\n"
        "pages under docs: 5 checked, 2 not in the site\n"
    )


def test_a_missing_page_alone_fails_the_link_check(tmp_path):
    root = linked_site(tmp_path, LANDING_REDIRECTS)
    unlisted(root, "guide/extra.md")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "docs/guide/extra.md: no page guide/extra.html " + NOT_BUILT.format("docs")
        + "checked 2 links on 3 pages, 0 broken\n"
        f"redirects of book.toml: 5 checked, 2 {NOT_FOLLOWED}, 0 broken\n"
        "pages under docs: 4 checked, 1 not in the site\n"
    )


def test_the_link_check_holds_a_readme_to_the_index_of_its_directory_and_the_not_found_source_to_its_page(tmp_path):
    root = linked_site(tmp_path, None)
    unlisted(root, "404.md", "guide/README.md")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "docs/404.md: no page 404.html " + NOT_BUILT.format("docs")
        + "docs/guide/README.md: no page guide/index.html " + NOT_BUILT.format("docs")
        + f"checked 2 links on 3 pages, 0 broken\nredirects of book.toml: 0 checked, 0 {NOT_FOLLOWED}, 0 broken\n"
        "pages under docs: 5 checked, 2 not in the site\n"
    )
    site = root / "target" / "book"
    (site / "404.html").write_text(page_html("Not found", "<p>No such page.</p>"), encoding="utf-8")
    (site / "guide" / "index.html").write_text(page_html("Guide", '<h1 id="guide">Guide</h1>'), encoding="utf-8")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (0, "")
    assert out.stdout == (
        f"checked 2 links on 4 pages, 0 broken\nredirects of book.toml: 0 checked, 0 {NOT_FOLLOWED}, 0 broken\n"
        "pages under docs: 5 checked, 0 not in the site\n"
    )


def test_the_link_check_reads_the_source_pages_of_the_tree_its_option_names(tmp_path):
    root = linked_site(tmp_path, None)
    unlisted(root, "guide/extra.md")
    git(root, "init", "-q")
    git(root, "add", "book.toml", "docs")
    git(root, "commit", "-q", "-m", "The pages")
    (root / "docs" / "guide" / "extra.md").unlink()
    unlisted(root, "notes.md")
    links_and_redirects = (
        f"checked 2 links on 3 pages, 0 broken\nredirects of book.toml: 0 checked, 0 {NOT_FOLLOWED}, 0 broken\n"
    )
    out = check_links(root, "--tree", "HEAD")
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "docs/guide/extra.md: no page guide/extra.html " + NOT_BUILT.format("docs")
        + links_and_redirects
        + "pages under docs in HEAD: 4 checked, 1 not in the site\n"
    )
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "docs/notes.md: no page notes.html " + NOT_BUILT.format("docs")
        + links_and_redirects
        + "pages under docs: 4 checked, 1 not in the site\n"
    )


def test_the_link_check_takes_no_directory_of_the_tree_for_a_page(tmp_path):
    root = linked_site(tmp_path, None)
    (root / "docs" / "figures.md").mkdir()
    (root / "docs" / "figures.md" / "plot.txt").write_text("A plot.\n", encoding="utf-8")
    git(root, "init", "-q")
    git(root, "add", "book.toml", "docs")
    git(root, "commit", "-q", "-m", "The pages")
    out = check_links(root, "--tree", "HEAD")
    assert (out.returncode, out.stderr) == (0, "")
    assert out.stdout.endswith("pages under docs in HEAD: 3 checked, 0 not in the site\n")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (0, "")
    assert out.stdout.endswith(ALL_BUILT)


def test_under_a_tree_the_link_check_reads_the_source_directory_from_the_root_of_that_tree(tmp_path):
    root = linked_site(tmp_path, None)
    unlisted(root, "guide/extra.md")
    git(root, "init", "-q")
    git(root, "add", "book.toml", "docs")
    git(root, "commit", "-q", "-m", "The pages")
    book = tmp_path / "elsewhere" / "book.toml"
    book.parent.mkdir()
    shutil.copyfile(root / "book.toml", book)
    out = check_links(root, "--book", str(book), "--tree", "HEAD")
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        "docs/guide/extra.md: no page guide/extra.html " + NOT_BUILT.format("docs")
        + f"checked 2 links on 3 pages, 0 broken\nredirects of {book}: 0 checked, 0 {NOT_FOLLOWED}, 0 broken\n"
        "pages under docs in HEAD: 4 checked, 1 not in the site\n"
    )


def test_the_link_check_reads_a_readme_in_another_letter_case_as_an_index_page(tmp_path):
    root = linked_site(tmp_path, None)
    unlisted(root, "guide/ReadMe.md")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout.startswith("docs/guide/ReadMe.md: no page guide/index.html " + NOT_BUILT.format("docs"))


def test_a_directory_at_the_address_of_a_page_is_not_that_page(tmp_path):
    root = linked_site(tmp_path, None)
    unlisted(root, "guide/extra.md")
    (root / "target" / "book" / "guide" / "extra.html").mkdir()
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout.startswith("docs/guide/extra.md: no page guide/extra.html " + NOT_BUILT.format("docs"))


def test_the_link_check_reads_the_source_directory_from_the_book_table_alone(tmp_path):
    root = linked_site(tmp_path, None)
    book = root / "book.toml"
    named = book.read_text(encoding="utf-8")
    site_url = 'site-url = "/vep-rs/"\n'
    assert 'src = "docs"\n' in named and site_url in named
    book.write_text(named.replace('src = "docs"\n', "").replace(site_url, site_url + 'src = "elsewhere"\n'), encoding="utf-8")
    (root / "docs").rename(root / "src")
    out = check_links(root)
    assert (out.returncode, out.stderr) == (0, "")
    assert out.stdout.endswith("pages under src: 3 checked, 0 not in the site\n")


@pytest.mark.parametrize(("key", "directory"), [('src = "pages"\n', "pages"), ("", "src")], ids=["named", "left-out"])
def test_the_link_check_reads_the_source_directory_from_the_book(tmp_path, key, directory):
    root = linked_site(tmp_path, None)
    book = root / "book.toml"
    named = book.read_text(encoding="utf-8")
    assert 'src = "docs"\n' in named
    book.write_text(named.replace('src = "docs"\n', key), encoding="utf-8")
    (root / "docs").rename(root / directory)
    unlisted(root, "notes.md", docs=directory)
    out = check_links(root)
    assert (out.returncode, out.stderr) == (1, "")
    assert out.stdout == (
        f"{directory}/notes.md: no page notes.html " + NOT_BUILT.format(directory)
        + f"checked 2 links on 3 pages, 0 broken\nredirects of book.toml: 0 checked, 0 {NOT_FOLLOWED}, 0 broken\n"
        f"pages under {directory}: 4 checked, 1 not in the site\n"
    )


def test_the_link_check_refuses_a_source_directory_without_a_page(tmp_path):
    root = linked_site(tmp_path, LANDING_REDIRECTS)
    git(root, "init", "-q")
    git(root, "add", "book.toml")
    git(root, "commit", "-q", "-m", "The book alone")
    refusal = "ERROR: [check_site_links] no Markdown page under {}, the source directory of the book book.toml\n"
    out = check_links(root, "--tree", "HEAD")
    assert (out.returncode, out.stdout, out.stderr) == (1, "", refusal.format("docs in HEAD"))
    for source in ("index.md", "guide/setup.md", "guide/page.md"):
        (root / "docs" / source).unlink()
    out = check_links(root)
    assert (out.returncode, out.stdout, out.stderr) == (1, "", refusal.format("docs"))
    shutil.rmtree(root / "docs")
    out = check_links(root)
    assert (out.returncode, out.stdout, out.stderr) == (1, "", refusal.format("docs"))


def test_the_link_check_reads_the_source_directory_of_the_committed_book():
    book = ROOT / "book.toml"
    assert link_check.read_source(str(book)) == tomllib.loads(book.read_text(encoding="utf-8"))["book"]["src"]


def test_the_release_site_script_has_the_link_check_read_the_table_it_built_from():
    script = (HERE / "build_release_site.sh").read_text(encoding="utf-8")
    runs = [line for line in script.splitlines() if "check_site_links.py" in line and not line.startswith("#")]
    arguments = '--site "$out/target/book" --book "$out/book.toml" --ref "$docs_ref" --tree "$ref"'
    assert runs == [f'python3 "$checkout/scripts/docs/check_site_links.py" {arguments}']


QUIET_BUILD = " INFO Book building has started\n INFO Running the html backend\n"
INDEX_WARNING = " WARN search index is very large (16413305 bytes)\n"
# What mdBook 0.5.4 writes for a code span that wraps so that `<dir>` opens a line: in a paragraph the
# first line alone, in a list item both.
UNCLOSED_IN_PARAGRAPH = " WARN unclosed HTML tag `<dir>` found in `developers/writing-a-plugin.md`\n"
UNCLOSED_IN_ITEM = (
    " WARN unclosed HTML tag `<dir>` found in `plugins.md` while exiting Item\n"
    "HTML tags must be closed before exiting a markdown element.\n"
)


def stand_in_builder(tmp_path, output, status=0):
    """A directory whose `mdbook` writes its arguments, one a line, to `arguments.txt` beside it and
    the value of `MDBOOK_LOG`, or `unset`, to `log_variable.txt`, then `output` to standard error,
    and ends with `status`. It calls shell built-ins alone, since the directory is the whole PATH of
    the run."""
    directory = tmp_path / "bin"
    directory.mkdir()
    (directory / "output.txt").write_text(output, encoding="utf-8")
    builder = directory / "mdbook"
    builder.write_text(
        "#!/bin/sh\n"
        "here=${0%/*}\n"
        'printf \'%s\\n\' "$@" > "$here/arguments.txt"\n'
        'printf \'%s\\n\' "${MDBOOK_LOG-unset}" > "$here/log_variable.txt"\n'
        'while IFS= read -r line; do printf \'%s\\n\' "$line" >&2; done < "$here/output.txt"\n'
        f"exit {status}\n",
        encoding="utf-8",
    )
    builder.chmod(0o755)
    return directory


def build_site(path, *arguments):
    """The finished run of the site build with `path` as its whole PATH."""
    return subprocess.run(
        [sys.executable, str(HERE / "site_build.py"), *arguments],
        env={**os.environ, "PATH": str(path)},
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize(
    ("warning", "message"),
    [
        (UNCLOSED_IN_PARAGRAPH, "unclosed HTML tag `<dir>` found in `developers/writing-a-plugin.md`"),
        (UNCLOSED_IN_ITEM, "unclosed HTML tag `<dir>` found in `plugins.md` while exiting Item"),
    ],
    ids=["in a paragraph", "in a list item"],
)
def test_the_site_build_refuses_a_page_that_leaves_a_tag_unclosed(tmp_path, warning, message):
    """The builder ends with status 0 under the warning; the site build names it and ends with 1."""
    run = build_site(stand_in_builder(tmp_path, QUIET_BUILD + warning + INDEX_WARNING))
    assert run.returncode == 1
    assert run.stderr == f"ERROR: [site_build] the builder warned: {message}\n"
    assert run.stdout == QUIET_BUILD + warning + INDEX_WARNING


@pytest.mark.parametrize("output", [QUIET_BUILD, QUIET_BUILD + INDEX_WARNING], ids=["none", "the admitted one"])
def test_the_site_build_passes_a_build_with_no_warning_but_an_admitted_one(tmp_path, output):
    run = build_site(stand_in_builder(tmp_path, output))
    assert (run.returncode, run.stderr, run.stdout) == (0, "", output)


def test_only_a_line_at_the_warning_level_is_a_warning_and_each_is_named_once():
    lines = [
        " INFO a page says WARN in its text\n",
        "WARN first word of a line, no padding\n",
        " WARN search index is very large (1 bytes)\n",
        "HTML tags must be closed before exiting a markdown element.\n",
        "  WARN  padded twice  \n",
    ]
    assert site_build.refused(lines) == ["first word of a line, no padding", "padded twice"]
    assert all(reason for reason in site_build.ADMITTED.values())


def test_the_site_build_ends_with_the_builders_status_when_the_build_fails(tmp_path):
    failed = "ERROR Rendering failed\n"
    run = build_site(stand_in_builder(tmp_path, QUIET_BUILD + UNCLOSED_IN_PARAGRAPH + failed, status=101))
    assert run.returncode == 101
    assert run.stderr == "ERROR: [site_build] mdbook build ended with status 101\n"
    assert run.stdout.endswith(failed)


def test_the_site_build_hands_its_arguments_to_the_builder(tmp_path):
    directory = stand_in_builder(tmp_path, QUIET_BUILD)
    assert build_site(directory, "--dest-dir", "target/other book").returncode == 0
    assert (directory / "arguments.txt").read_text(encoding="utf-8") == "build\n--dest-dir\ntarget/other book\n"


def test_the_site_build_runs_the_builder_without_its_log_variable(monkeypatch, tmp_path):
    """mdBook reads `MDBOOK_LOG`: `error` silences its warnings, and any level puts each line's
    source before its message, which then opens with no entry of `ADMITTED`. The builder runs
    without the variable."""
    monkeypatch.setenv("MDBOOK_LOG", "error")
    directory = stand_in_builder(tmp_path, QUIET_BUILD + UNCLOSED_IN_PARAGRAPH)
    assert build_site(directory).returncode == 1
    assert (directory / "log_variable.txt").read_text(encoding="utf-8") == "unset\n"


def test_the_site_build_writes_a_character_its_standard_output_cannot_encode_as_its_escape(monkeypatch, tmp_path):
    """A line of the builder that the encoding of standard output cannot hold is passed through with
    the character escaped, and the warning on it is still refused."""
    monkeypatch.setenv("PYTHONIOENCODING", "ascii")
    message = "unclosed HTML tag `<dir>` found in `{}`"
    named, escaped = message.format("caf" + chr(0xE9) + ".md"), message.format("caf\\xe9.md")
    run = build_site(stand_in_builder(tmp_path, QUIET_BUILD + f" WARN {named}\n"))
    assert run.returncode == 1
    assert run.stdout == QUIET_BUILD + f" WARN {escaped}\n"
    assert run.stderr == f"ERROR: [site_build] the builder warned: {escaped}\n"


def test_the_site_build_refuses_to_run_without_the_builder(tmp_path):
    run = build_site(tmp_path)
    assert (run.returncode, run.stdout) == (2, "")
    assert run.stderr == "ERROR: [site_build] mdbook is not on PATH\n"


def test_the_pull_request_step_of_the_pages_job_builds_through_the_site_build():
    """The step builds the pages with the script that refuses a warning, never with the builder alone,
    and builds before the checks that read the built site. The builder the job installs is of the
    series whose warning lines the script reads."""
    workflow = (ROOT / ".github" / "workflows" / "pages.yml").read_text(encoding="utf-8")
    installed = re.findall(r"^\s*MDBOOK_VERSION=(\S+)$", workflow, re.MULTILINE)
    assert len(installed) == 1, installed
    assert installed[0].startswith(site_build.MDBOOK_SERIES + "."), installed[0]
    step = workflow[workflow.index("      - name: Build the pull request's pages\n") :]
    commands = [line.strip() for line in step[: step.index("\n      - name: ", 1)].splitlines()]
    assert "python3 scripts/docs/site_build.py" in commands
    assert not [command for command in commands if command.split()[:1] == ["mdbook"]]
    built = commands.index("python3 scripts/docs/site_build.py")
    assert built < min(i for i, command in enumerate(commands) if "check_site_links.py" in command)


def test_readme_home_points_the_edit_link_at_the_readme(monkeypatch, capsys, tmp_path):
    site = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, readme_home="1")
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f'<a href="{EDIT_URL}README.md">Edit</a>' in index
    assert "docs/index.md" not in index


def probe_record(**overrides):
    rec = {
        "tiles": ["index-375-light-top.png", "index-375-light-s1.png"],
        "overflow": False,
        "overflowing_total": 0,
        "wrapped_inline_code": [],
        "long_span_wraps": [{"text": "a long span"}],
        "scrollers": [{"kind": "pre", "scrollbar_visible": False}],
        "min_font_px": 14.4,
        "min_font_page_px": 14.4,
        "measure_px": 343,
        "chars_per_line": 68.2,
        "chars_per_line_text": 72.5,
    }
    rec.update(overrides)
    return rec


def wide_record(**overrides):
    return probe_record(
        tiles=["index-1280-light-top.png"], overflow=True, overflowing_total=3,
        wrapped_inline_code=[{"text": "a"}, {"text": "b"}], long_span_wraps=[], scrollers=[],
        min_font_px=16, min_font_page_px=12, chars_per_line_text=95, **overrides,
    )


def probe(narrow=None, wide=None):
    return {
        "site": "/site", "base": "/vep-rs/", "themes": ["light", "navy"], "widths": [375, 1280],
        "full_widths": [375, 1280], "axe": "skipped",
        "pages": {"index": {"slug": "index", "widths": {"375": narrow or probe_record(), "1280": wide or wide_record()}}},
    }


SUMMARY_HEAD = (
    "Probe of /site at /vep-rs/, themes light, navy.\n\n"
    "| Page | Width | Tiles | Overflow | Wrapped code | Long-span wraps | Scrollers without a bar | Min font px | Chars per line |\n"
    "| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |\n"
)


def test_tile_names_join_the_stem_and_the_tile():
    assert shots.tile_name("guide-setup-375-light", "top") == "guide-setup-375-light-top.png"
    assert shots.tile_path("out", "guide-setup-375-light", "s2") == os.path.join("out", "guide-setup-375-light-s2.png")


def test_summary_md_has_a_row_per_page_and_width_then_the_totals():
    assert shots.summary_md(probe()) == SUMMARY_HEAD + (
        "| index | 375 | 2 | no | 0 | 1 | 1 | 14.4 | 72.5 |\n"
        "| index | 1280 | 1 | yes | 2 | 0 | 0 | 16 | 95 |\n\n"
        "- overflow: 1 of 2 cells\n"
        "- wrapped inline code that fits its block: 2 in 1 of 2 cells\n"
        "- long-span wraps: 1 in 1 of 2 cells\n"
        "- scrollers without a visible bar: 1 in 1 of 2 cells\n"
        "- minimum font: 14.4 px (index at 375 px)\n"
        "- characters per line outside 60 to 80 at 768 px and wider: 1 of 2 cells\n"
    )


def test_summary_md_shows_long_span_wraps_as_not_available_when_a_probe_lacks_them():
    narrow, wide = probe_record(), wide_record()
    del narrow["long_span_wraps"], wide["long_span_wraps"]
    assert shots.summary_md(probe(narrow, wide)) == SUMMARY_HEAD + (
        "| index | 375 | 2 | no | 0 | n/a | 1 | 14.4 | 72.5 |\n"
        "| index | 1280 | 1 | yes | 2 | n/a | 0 | 16 | 95 |\n\n"
        "- overflow: 1 of 2 cells\n"
        "- inline code on two or more lines, fit not recorded: 2 in 1 of 2 cells\n"
        "- long-span wraps: not recorded in this probe\n"
        "- scrollers without a visible bar: 1 in 1 of 2 cells\n"
        "- minimum font: 14.4 px (index at 375 px)\n"
        "- characters per line outside 60 to 80 at 768 px and wider: 1 of 2 cells\n"
    )


def test_summary_html_marks_the_bad_measures():
    assert shots.summary_html(probe_record()) == (
        "overflow no | 68.2 n/line, 72.5 chars/line (343 px) | min font 14.4 px | "
        '<b class="bad">1 pre scroll (1 without scrollbar)</b> | 0 table scroll (0 without scrollbar) | '
        "0 inline code wrapped | 1 long-span wraps"
    )
    axe = {"light": {"count": 1, "violations": [{"id": "color-contrast", "nodes": 2}]}}
    assert shots.summary_html(wide_record(axe=axe)) == (
        '<b class="bad">overflow yes (3 elements)</b> | 68.2 n/line, 95 chars/line (343 px) | '
        "min font 16 px (page 12 px) | 0 pre scroll (0 without scrollbar) | 0 table scroll (0 without scrollbar) | "
        '<b class="bad">2 inline code wrapped</b> | 0 long-span wraps | '
        '<b class="bad">axe light: 1 violations (color-contrast 2)</b>'
    )


def test_summary_of_probe_files_prints_their_tables_without_playwright(monkeypatch, capsys, tmp_path):
    path = tmp_path / "probe.json"
    path.write_text(json.dumps(probe()), encoding="utf-8")
    monkeypatch.setitem(sys.modules, "playwright", None)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", None)
    monkeypatch.setattr(sys, "argv", ["site_screenshots.py", "--summary", str(path)])
    shots.main()
    assert capsys.readouterr().out == f"## {path}\n\n{shots.summary_md(probe())}\n"


HOME_README = f"# Project\n\n## Current release: {RELEASE}\n"
PROVENANCE = "docs/concordance-provenance"
POPULATION_RECORD = "first-population-v1.2.3.json"
RELEASE_RECORD = "first-release-v1.2.3.json"
PLUGIN_RECORD = "first-plugins-v1.2.3.json"
HOME_IMAGE = "docs/assets/comparison.png"
HOME_VALUES = "docs/assets/comparison.values.json"
HOME_NOTE = "Measured on made-up machines."
RELEASED_TOOL = ("vep-rs", "vep-rs", RELEASE)
COMPARATOR_TOOL = ("fastvep", "fastVEP", "v9.8.7")
REFERENCE_TOOL = ("ensembl-vep", "Ensembl VEP", "100.1")
# The page's sentence on a reference whose wall time stands at another version, here at 100.0.
HOME_TIMED = "Ensembl VEP's wall times are Ensembl VEP 100.0's figures, carried forward from the earlier reference generation."
PANEL_A = "raw F1 per dataset"
PANEL_B = "wall time per dataset and architecture"
PANEL_C = "raw F1 per plugin"
PANEL_D = "the plugins each tool reads"
PANEL_E = "the ready-made Ensembl VEP caches that are exact, of those checked"
FOUR_PANELS = f"{PANEL_A}, {PANEL_B}, {PANEL_C} and {PANEL_D}"
FIVE_PANELS = f"{PANEL_A}, {PANEL_B}, {PANEL_C}, {PANEL_D} and {PANEL_E}"
CHECK_REFUSAL = (
    "ERROR: [render_home_page] docs/index.md line {} is not the line the records of v1.2.3 give; "
    "scripts/docs/render_home_page.py writes the page"
)


def home_dataset(key="gnomad_v4.1_genomes", label="Cohort A v7, whole genome", variants=4321987, arm=(600, 20), x86=(900, 45)):
    """A dataset of a population record, by default the one the page's first reason reads. `arm` and `x86` give
    the wall-time sums of Ensembl VEP and of vep-rs on that architecture, or None where the record times neither."""
    return {"dataset": key, "label": label, "variants": variants, "arm64": arm, "x86_64": x86}


COHORT_B = home_dataset("cohort_b_genomes", "Cohort B, whole genome", 1234, arm=(500, 25), x86=(770, 35))
HOME_DATASETS = (home_dataset(), COHORT_B)


def home_population(datasets=HOME_DATASETS, release=RELEASE):
    """A population record of the datasets, with one wall-time row per engine, dataset and architecture."""
    def rows(engine):
        return [
            {"dataset": entry["dataset"], "arch": arch, "sum_of_cell_medians_sec": entry[arch][engine]}
            for entry in datasets for arch in ("arm64", "x86_64") if entry[arch] is not None
        ]

    return {
        "record": "population-measurements",
        "release": release,
        "datasets": [{key: entry[key] for key in ("dataset", "label", "variants")} for entry in datasets],
        "wall_time": {"ensembl_vep": {"per_dataset": rows(0)}, "vep_rs": {"per_dataset": rows(1)}},
    }


def home_release(references=("100.1",), release=RELEASE):
    """A release record with the two genome-wide ClinVar cells and one reference image per entry of
    `references`, tagged `release_<entry>`."""
    return {
        "record": "release-addendum",
        "release": release,
        "perl_reference": {
            f"family_{index}": {"image": f"registry.example/ensembl-vep:release_{reference}"}
            for index, reference in enumerate(references, start=1)
        },
        "concordance": {"s01": {}, "s02": {}},
    }


def home_plugins(read, release=RELEASE):
    """A plugin record whose support block of the released tool reads `read` plugins."""
    return {
        "record_type": "plugins",
        "release": release,
        "support": [{"engine": "vep-rs", "read": read}, {"engine": "fastvep", "read": 3}],
    }


def home_values(tools=(RELEASED_TOOL, COMPARATOR_TOOL, REFERENCE_TOOL), notes=(HOME_NOTE,), panels="abcd", image=HOME_IMAGE,
                panel_notes=None):
    """The values file of a release figure under every key the figure script writes: its tools as (id, name,
    version) and, for a tool whose wall time stands at another version, that version as a fourth item, its
    notes, its panel ids and, for each panel id of `panel_notes`, that panel's note, beside the caption and the
    alternative text of the README, which the page does not copy."""
    panel_notes = panel_notes or {}
    return {
        "image": image,
        "canvas": {"width_px": 100, "height_px": 125},
        "alt_text": "The alternative text the README shows.",
        "caption": "The caption the README shows.",
        "notes": list(notes),
        "sources": [],
        "tools": [
            {"id": tool_id, "name": name, "version": version} | ({"wall_time_version": timed[0]} if timed else {}) | {"colour": "#000000"}
            for tool_id, name, version, *timed in tools
        ],
        "records": [],
        "panels": [
            {"id": panel, "title": f"Panel {panel}", "groups": []} | ({"note": panel_notes[panel]} if panel in panel_notes else {})
            for panel in panels
        ],
    }


def figure_tree(tmp_path, genomes=False, timed_apart=False):
    """A made-up tree of one release whose records the figure script and the page's script both read, with the
    figure's pair: the values file as the figure script's values layer writes it from those records, and made-up
    image bytes beside it. With `genomes` the release record holds the block of the caches checked, which gives
    the figure its fifth panel; with `timed_apart` the timing records name the reference at another version
    than its concordance records. Returns the tree and the values file's document."""
    made_up = figure_records.with_genomes() if genomes else figure_records.record_set(**(figure_records.timed_at() if timed_apart else {}))
    records = json.loads(json.dumps(made_up).replace('"pop_a"', '"gnomad_v4.1_genomes"'))
    for record in records.values():
        if record.get("record") == "population-measurements":
            for entry in record["datasets"]:
                entry["variants"] = 4321987
    repo = tmp_path / "repo"
    figure_records.write_records(repo / PROVENANCE, records)
    (repo / "README.md").write_text(f"# Project\n\n## Current release: {figure_records.RELEASE}\n", encoding="utf-8")
    figure.main(["--repo", str(repo), "--values-only"])
    (repo / HOME_IMAGE).write_bytes(b"made-up image bytes")
    return repo, json.loads((repo / HOME_VALUES).read_text(encoding="utf-8"))


def home_tree(tmp_path, *, readme=HOME_README, datasets=HOME_DATASETS, references=("100.1",), plugins=None,
              values=None, png=None, records=None):
    """A made-up tree of the release `RELEASE` names: the README, a population record of `datasets` and a
    release record of `references`; the plugin record when `plugins` gives the released tool's count; the
    figure's values file when `values` gives it, with its image beside it unless `png` says otherwise.
    `records` adds files to the provenance directory, and leaves out a file it gives as None."""
    repo = tmp_path / "repo"
    (repo / PROVENANCE).mkdir(parents=True)
    (repo / "docs" / "assets").mkdir()
    (repo / "README.md").write_text(readme, encoding="utf-8")
    files = {POPULATION_RECORD: home_population(datasets), RELEASE_RECORD: home_release(references)}
    if plugins is not None:
        files[PLUGIN_RECORD] = home_plugins(plugins)
    files.update(records or {})
    for name, record in files.items():
        if record is not None:
            (repo / PROVENANCE / name).write_text(json.dumps(record), encoding="utf-8")
    if values is not None:
        (repo / HOME_VALUES).write_text(json.dumps(values), encoding="utf-8")
    with_image = values is not None if png is None else png
    if with_image:
        (repo / HOME_IMAGE).write_bytes(b"made-up image bytes")
    return repo


def home_build(repo, template=None, registry=None):
    """What the script builds for the tree, from the committed template and registry unless a test passes
    its own."""
    text = home.TEMPLATE.read_text(encoding="utf-8") if template is None else template
    return home.build(repo, text, registry or home.REGISTRY)


def home_page(repo, template=None, registry=None):
    return home_build(repo, template, registry)[0]


def home_reasons(page):
    """The reasons of a page, one line each."""
    return [line for line in page.splitlines() if line.startswith("- **")]


def figure_paragraphs(alt, caption):
    """The figure block as the paragraphs of a page: the image with its alternative text inside a link to the
    site's own copy of it, then the caption."""
    return [
        '<figure class="home-figure">\n<a href="assets/comparison.png">',
        f"![{alt}](assets/comparison.png)",
        "</a>\n<figcaption>",
        caption,
        "</figcaption>\n</figure>",
    ]


def files_of(repo):
    """Every file under the tree, with its bytes and the time of its last write."""
    return {
        path.relative_to(repo).as_posix(): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in repo.rglob("*") if path.is_file()
    }


def readme_row(readme, name):
    """The second cell of the README table row whose first cell opens with `name` in bold, or None."""
    found = re.search(rf"^\| \*\*{re.escape(name)}[^|\n]*\|([^|\n]*)\|", readme, re.MULTILINE)
    return found.group(1) if found else None


def test_the_committed_home_page_is_the_page_its_records_give(capsys):
    assert (ROOT / "docs" / "index.md").read_bytes() == home_page(ROOT).encode("utf-8")
    home.main(["--repo", str(ROOT), "--check"])
    version = home.current_release(ROOT / "README.md")
    assert capsys.readouterr().out == f"docs/index.md is the page the records of v{version} give\n"


def test_the_committed_home_page_states_the_figures_of_the_readme_rows():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    page = (ROOT / "docs" / "index.md").read_text(encoding="utf-8")
    speed, genome = readme_row(readme, "Speed vs Ensembl VEP"), readme_row(readme, "All of gnomAD v4.1")
    if speed is None and genome is None:
        pytest.skip("the README has neither the speed row nor the whole-genome row")
    rows, stated = {}, {}
    if speed is not None:
        upper = re.findall(r" to (\d[\d,.]*)x per (?:population )?dataset", speed)
        rows["speed number"] = [figure.replace(",", "") for figure in upper]
        stated["speed number"] = [figure.replace(",", "") for figure in re.findall(r"(\d[\d,]*) times faster than Ensembl VEP", page)]
    if genome is not None:
        figures = r"(\d[\d,]*) variants\b.*?\b(\d[\d,]*) s on one ARM machine"
        rows["first reason"] = re.findall(figures, genome)
        stated["first reason"] = re.findall(figures, page)
    assert all(rows.values())
    assert stated == rows


def test_check_passes_on_the_page_the_script_wrote_and_writes_nothing(capsys, tmp_path):
    repo = home_tree(tmp_path, plugins=17, values=home_values())
    home.main(["--repo", str(repo)])
    written = files_of(repo)
    assert "docs/index.md" in written
    capsys.readouterr()
    home.main(["--repo", str(repo), "--check"])
    assert capsys.readouterr().out == "docs/index.md is the page the records of v1.2.3 give\n"
    assert files_of(repo) == written


def test_check_names_the_line_of_a_changed_digit_and_leaves_the_page_as_it_was(tmp_path):
    repo = home_tree(tmp_path)
    home.main(["--repo", str(repo)])
    page = repo / "docs" / "index.md"
    lines = page.read_text(encoding="utf-8").splitlines(keepends=True)
    number = next(index for index, line in enumerate(lines, start=1) if "4,321,987 variants" in line)
    lines[number - 1] = lines[number - 1].replace("4,321,987", "4,321,988")
    changed = "".join(lines)
    page.write_text(changed, encoding="utf-8")
    with pytest.raises(SystemExit) as excinfo:
        home.main(["--repo", str(repo), "--check"])
    assert excinfo.value.code == CHECK_REFUSAL.format(number)
    assert page.read_text(encoding="utf-8") == changed


def test_check_names_the_last_line_of_a_page_without_it_or_without_its_newline(tmp_path):
    repo = home_tree(tmp_path)
    home.main(["--repo", str(repo)])
    page = repo / "docs" / "index.md"
    lines = page.read_text(encoding="utf-8").splitlines(keepends=True)
    for held in ("".join(lines[:-1]), "".join(lines).removesuffix("\n")):
        page.write_text(held, encoding="utf-8")
        with pytest.raises(SystemExit) as excinfo:
            home.main(["--repo", str(repo), "--check"])
        assert excinfo.value.code == CHECK_REFUSAL.format(len(lines))
        assert page.read_text(encoding="utf-8") == held


def test_the_speed_number_is_the_largest_quotient_over_the_datasets_and_both_architectures(tmp_path):
    repo = home_tree(tmp_path, datasets=(home_dataset(), {**COHORT_B, "x86_64": (2460, 60)}))
    assert " runs up to 41 times faster than Ensembl VEP " in home_page(repo)
    row = "wall_time.{}.per_dataset[dataset=cohort_b_genomes, arch=x86_64].sum_of_cell_medians_sec"
    assert home_build(repo)[2]["speed-up"][2] == f"{row.format('ensembl_vep')} / {row.format('vep_rs')}"


def test_a_half_rounds_up_in_the_speed_number_and_the_seconds_from_the_records_decimal_text(tmp_path):
    # 1174.32 over 10.08 is 116.5 in decimal, and the quotient of the two doubles lies below it.
    assert 1174.32 / 10.08 < 116.5
    datasets = (home_dataset(arm=(1000, 12.5)), {**COHORT_B, "arm64": (1174.32, 10.08)})
    page = home_page(home_tree(tmp_path, datasets=datasets))
    assert " runs up to 117 times faster than Ensembl VEP " in page
    assert home_reasons(page)[0].endswith(" annotates in 13 s on one ARM machine.")


def test_the_variants_the_seconds_and_the_speed_number_carry_thousands_separators(tmp_path):
    page = home_page(home_tree(tmp_path, datasets=(home_dataset(arm=(3086000, 1234.4)), COHORT_B)))
    assert " runs up to 2,500 times faster than Ensembl VEP " in page
    assert home_reasons(page)[0].endswith(", 4,321,987 variants, annotates in 1,234 s on one ARM machine.")


def test_the_count_word_counts_the_population_datasets_and_the_two_clinvar_cells(tmp_path):
    third = home_dataset("cohort_c_genomes", "Cohort C, whole genome", 56, arm=(300, 30), x86=(400, 40))
    page = home_page(home_tree(tmp_path, datasets=(*HOME_DATASETS, third)))
    assert home_reasons(page)[1].endswith(", on five whole-genome datasets.")


def test_the_label_loses_its_whole_genome_ending(tmp_path):
    page = home_page(home_tree(tmp_path))
    assert home_reasons(page)[0].endswith(" All of Cohort A v7, 4,321,987 variants, annotates in 20 s on one ARM machine.")
    assert ", whole genome" not in page


def test_the_reference_release_is_the_image_tag_without_its_prefix(tmp_path):
    page = home_page(home_tree(tmp_path, references=("100.1", "100.1")))
    assert home_reasons(page)[1].endswith(" Release 100.1, on four whole-genome datasets.")
    assert "release_" not in page


@pytest.mark.parametrize(
    ("values", "alt", "caption"),
    [
        (
            {},
            f"{HOME_NOTE} Four panels comparing vep-rs v1.2.3 and fastVEP v9.8.7 against Ensembl VEP 100.1: {FOUR_PANELS}.",
            f"{HOME_NOTE} vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 as comparator.",
        ),
        (
            {"tools": (RELEASED_TOOL, REFERENCE_TOOL)},
            f"{HOME_NOTE} Four panels comparing vep-rs v1.2.3 against Ensembl VEP 100.1: {FOUR_PANELS}.",
            f"{HOME_NOTE} vep-rs v1.2.3 against Ensembl VEP 100.1.",
        ),
        (
            {"notes": ()},
            f"Four panels comparing vep-rs v1.2.3 and fastVEP v9.8.7 against Ensembl VEP 100.1: {FOUR_PANELS}.",
            "vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 as comparator.",
        ),
        (
            {"notes": ("First.", "Second.")},
            f"First. Second. Four panels comparing vep-rs v1.2.3 and fastVEP v9.8.7 against Ensembl VEP 100.1: {FOUR_PANELS}.",
            "First. Second. vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 as comparator.",
        ),
        (
            {"panel_notes": {"d": "Under the last panel.", "b": "Under the second panel."}},
            (
                f"{HOME_NOTE} Under the second panel. Under the last panel. Four panels comparing vep-rs v1.2.3 and "
                f"fastVEP v9.8.7 against Ensembl VEP 100.1: {FOUR_PANELS}."
            ),
            (
                f"{HOME_NOTE} Under the second panel. Under the last panel. vep-rs v1.2.3 against Ensembl VEP 100.1, "
                "with fastVEP v9.8.7 as comparator."
            ),
        ),
        (
            {"panels": "b"},
            f"{HOME_NOTE} One panel comparing vep-rs v1.2.3 and fastVEP v9.8.7 against Ensembl VEP 100.1: {PANEL_B}.",
            f"{HOME_NOTE} vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 as comparator.",
        ),
        (
            {"panels": "da"},
            f"{HOME_NOTE} Two panels comparing vep-rs v1.2.3 and fastVEP v9.8.7 against Ensembl VEP 100.1: {PANEL_D} and {PANEL_A}.",
            f"{HOME_NOTE} vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 as comparator.",
        ),
        (
            {"panels": "abcde"},
            f"{HOME_NOTE} Five panels comparing vep-rs v1.2.3 and fastVEP v9.8.7 against Ensembl VEP 100.1: {FIVE_PANELS}.",
            f"{HOME_NOTE} vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 as comparator.",
        ),
        (
            {"notes": ("On *one* machine [made up].",), "tools": (RELEASED_TOOL, ("fastvep", "fast_VEP*", "v9.8.7"), REFERENCE_TOOL)},
            (
                "On \\*one\\* machine \\[made up\\]. Four panels comparing vep-rs v1.2.3 and fast\\_VEP\\* v9.8.7 against "
                f"Ensembl VEP 100.1: {FOUR_PANELS}."
            ),
            "On \\*one\\* machine \\[made up\\]. vep-rs v1.2.3 against Ensembl VEP 100.1, with fast\\_VEP\\* v9.8.7 as comparator.",
        ),
        (
            {"tools": (RELEASED_TOOL, COMPARATOR_TOOL, (*REFERENCE_TOOL, "100.0"))},
            f"{HOME_NOTE} Four panels comparing vep-rs v1.2.3 and fastVEP v9.8.7 against Ensembl VEP 100.1: {FOUR_PANELS}. {HOME_TIMED}",
            f"{HOME_NOTE} vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 as comparator. {HOME_TIMED}",
        ),
    ],
    ids=[
        "one-comparator", "no-comparator", "no-notes", "two-notes", "panel-notes-after-the-notes", "one-panel",
        "panels-in-the-file-order", "five-panels", "markup-escaped", "reference-timed-at-another-version",
    ],
)
def test_the_figure_block_is_composed_from_the_values_file(tmp_path, values, alt, caption):
    page = home_page(home_tree(tmp_path, values=home_values(**values)))
    assert page.split("\n\n")[3:9] == [*figure_paragraphs(alt, caption), "## Why vep-rs"]


def test_the_page_is_built_from_the_values_file_the_figure_script_writes(tmp_path):
    repo, document = figure_tree(tmp_path)
    made_up = home_values(panel_notes={"d": "A note."})
    assert set(made_up) == set(document)
    assert [set(tool) for tool in made_up["tools"]] == [set(tool) for tool in document["tools"]]
    assert [set(panel) for panel in made_up["panels"]] == [set(panel) for panel in document["panels"]]
    released, comparator, reference = (f"{tool['name']} {tool['version']}" for tool in document["tools"])
    note = next(panel["note"] for panel in document["panels"] if panel["id"] == "d")
    lead = f"{released} against {reference}, with {comparator} as comparator."
    alt = f"{note} Four panels comparing {released} and {comparator} against {reference}: {FOUR_PANELS}."
    assert home_page(repo).split("\n\n")[3:9] == [*figure_paragraphs(alt, f"{note} {lead}"), "## Why vep-rs"]
    # The page's own sentences on the figure open the sentences the figure script writes for the README.
    assert document["notes"] == []
    assert document["caption"].startswith(f"**{lead}** ")
    assert document["alt_text"].startswith(f"Four panels comparing {released} and {comparator} against {reference}. ")
    for panel_id, clause in zip("abcd", (PANEL_A, PANEL_B, PANEL_C, PANEL_D)):
        assert f" ({panel_id}) {clause[0].upper()}{clause[1:]}" in document["caption"]


def test_the_page_counts_five_panels_from_a_values_file_whose_release_record_holds_the_caches_checked(tmp_path):
    repo, document = figure_tree(tmp_path, genomes=True)
    assert [panel["id"] for panel in document["panels"]] == list("abcde")
    released, comparator, reference = (f"{tool['name']} {tool['version']}" for tool in document["tools"])
    note = next(panel["note"] for panel in document["panels"] if panel["id"] == "d")
    lead = f"{released} against {reference}, with {comparator} as comparator."
    alt = f"{note} Five panels comparing {released} and {comparator} against {reference}: {FIVE_PANELS}."
    assert home_page(repo).split("\n\n")[3:9] == [*figure_paragraphs(alt, f"{note} {lead}"), "## Why vep-rs"]
    # Each clause of the page opens its panel's part of the caption the figure script writes, the fifth too,
    # and the page's clause for the fifth panel names no tool.
    assert document["alt_text"].startswith(f"Five panels comparing {released} and {comparator} against {reference}. ")
    for panel_id, clause in zip("abcde", (PANEL_A, PANEL_B, PANEL_C, PANEL_D, PANEL_E)):
        assert f" ({panel_id}) {clause[0].upper()}{clause[1:]}" in document["caption"]
    assert PANEL_E == figure.GENOMES_CLAUSE and not any(tool["name"] in PANEL_E for tool in document["tools"] if tool["id"] != "ensembl-vep")


def test_the_page_names_the_version_the_references_wall_time_stands_at_in_the_figure_scripts_words(tmp_path):
    repo, document = figure_tree(tmp_path, timed_apart=True)
    released, comparator, reference = (f"{tool['name']} {tool['version']}" for tool in document["tools"])
    stated = document["tools"][2]["wall_time_version"]
    assert stated == figure_records.TIMED_VERSION != document["tools"][2]["version"]
    timed = f"Ensembl VEP's wall times are Ensembl VEP {stated}'s figures, carried forward from the earlier reference generation."
    note = next(panel["note"] for panel in document["panels"] if panel["id"] == "d")
    lead = f"{released} against {reference}, with {comparator} as comparator. {timed}"
    alt = f"{note} Four panels comparing {released} and {comparator} against {reference}: {FOUR_PANELS}. {timed}"
    assert home_page(repo).split("\n\n")[3:9] == [*figure_paragraphs(alt, f"{note} {lead}"), "## Why vep-rs"]
    # The page's sentence is the first of the two the figure script writes after the wall-time panel's part of
    # its caption and of its alternative text, and the page's lead names the version of the second.
    for text in (document["caption"], document["alt_text"]):
        assert text.count(f". {timed} The concordance panels compare against {reference}'s output. (c) ") == 1
    _, parts, places, _ = home_build(repo)
    assert parts["timed-apart"] is True
    assert places["reference-timed"][:3] == (f"Ensembl VEP {stated}", HOME_VALUES, "tools, by the role of each id in comparators.tsv: wall_time_version")
    # The same records at one version give no such sentence: the part is left out with its place.
    same, _ = figure_tree(tmp_path / "same")
    _, parts, places, _ = home_build(same)
    assert parts["timed-apart"] is False and places["reference-timed"] is None and "wall times are" not in home_page(same)


def test_two_comparators_are_named_together_in_the_figure_block(tmp_path):
    registry = tmp_path / "comparators.tsv"
    registry.write_text(
        "id\trole\trecord_keys\nvep-rs\treleased\tvep_rs\nfastvep\tcomparator\tfastvep\n"
        "othervep\tcomparator\tothervep\nensembl-vep\treference\tensembl_vep\n",
        encoding="utf-8",
    )
    tools = (RELEASED_TOOL, COMPARATOR_TOOL, ("othervep", "otherVEP", "v4.5.6"), REFERENCE_TOOL)
    page = home_page(home_tree(tmp_path, values=home_values(tools=tools)), registry=registry)
    assert page.split("\n\n")[3:9] == [
        *figure_paragraphs(
            (
                f"{HOME_NOTE} Four panels comparing vep-rs v1.2.3, fastVEP v9.8.7 and otherVEP v4.5.6 against "
                f"Ensembl VEP 100.1: {FOUR_PANELS}."
            ),
            f"{HOME_NOTE} vep-rs v1.2.3 against Ensembl VEP 100.1, with fastVEP v9.8.7 and otherVEP v4.5.6 as comparators.",
        ),
        "## Why vep-rs",
    ]


def test_a_tree_without_the_figure_pair_has_no_figure_block(tmp_path):
    page = home_page(home_tree(tmp_path))
    assert "<figure" not in page
    assert "comparison.png" not in page
    title, opening, second, following = page.split("\n\n")[:4]
    assert title == "# vep-rs"
    assert " runs up to 30 times faster than Ensembl VEP " in opening
    assert second.startswith("It does the work of Ensembl VEP, ")
    assert following == "## Why vep-rs"


def test_the_paragraph_under_the_opening_says_what_the_tool_is_and_links_the_install_lines(tmp_path):
    page = home_page(home_tree(tmp_path, plugins=17, values=home_values()))
    second = page.split("\n\n")[2]
    assert second.startswith("It does the work of Ensembl VEP, the Variant Effect Predictor: for each variant of a VCF file ")
    assert "[Install](#install)" in second
    assert "\n## Install\n" in page
    assert page.count("Variant Effect Predictor") == 1
    assert page.count("VCF") == 1


@pytest.mark.parametrize(
    ("tree", "message"),
    [
        (
            {"png": True},
            f"the tree holds {HOME_IMAGE} without {HOME_VALUES}; the page shows the figure from both",
        ),
        (
            {"values": home_values(), "png": False},
            f"the tree holds {HOME_VALUES} without {HOME_IMAGE}; the page shows the figure from both",
        ),
        (
            {"values": home_values(tools=(("vep-rs", "vep-rs", "v1.2.2"), COMPARATOR_TOOL, REFERENCE_TOOL))},
            f"{HOME_VALUES} draws vep-rs v1.2.2; the README names v1.2.3",
        ),
        (
            {"values": home_values(image="docs/assets/other.png")},
            f"{HOME_VALUES} names the image 'docs/assets/other.png'; the page shows {HOME_IMAGE}",
        ),
        (
            {"values": home_values(panels="az")},
            f"{HOME_VALUES} holds the panel z, and home_page.template.md has no '<<panel z>>' clause",
        ),
        (
            {"values": home_values(tools=(RELEASED_TOOL, ("othervep", "otherVEP", "v4.5.6"), REFERENCE_TOOL))},
            f"{HOME_VALUES} tools[1] names othervep, which comparators.tsv does not list",
        ),
        (
            {"records": {"second-population-v1.2.3.json": home_population()}},
            (
                f"2 files of {PROVENANCE} are the population record (record: population-measurements) of v1.2.3: "
                "first-population-v1.2.3.json, second-population-v1.2.3.json"
            ),
        ),
        (
            {"records": {RELEASE_RECORD: None}},
            f"no file of {PROVENANCE} is the release record (record: release-addendum) of v1.2.3",
        ),
        (
            {"plugins": 17, "records": {"second-plugins-v1.2.3.json": home_plugins(18)}},
            (
                f"2 files of {PROVENANCE} are the plugin record (record_type: plugins) of v1.2.3: "
                "first-plugins-v1.2.3.json, second-plugins-v1.2.3.json"
            ),
        ),
        ({"readme": "# Project\n"}, "README.md has 0 '## Current release: vX.Y.Z' headings; one is read"),
        (
            {"readme": HOME_README + "\n## Current release: v1.2.2\n"},
            "README.md has 2 '## Current release: vX.Y.Z' headings; one is read",
        ),
        (
            {"datasets": (COHORT_B,)},
            f"{PROVENANCE}/{POPULATION_RECORD} datasets holds no gnomad_v4.1_genomes",
        ),
        (
            {"datasets": (home_dataset(), {**COHORT_B, "x86_64": None})},
            f"{PROVENANCE}/{POPULATION_RECORD} wall_time.vep_rs.per_dataset has no row of cohort_b_genomes on x86_64",
        ),
        (
            {"references": ("100.1", "100.2")},
            f"{PROVENANCE}/{RELEASE_RECORD} perl_reference names the releases 100.1 and 100.2; the page names one",
        ),
    ],
    ids=[
        "image-without-values", "values-without-image", "figure-of-another-release", "another-image-path",
        "panel-without-a-clause", "tool-outside-the-registry", "two-population-records", "no-release-record",
        "two-plugin-records", "no-heading", "two-headings", "no-gnomad-dataset", "dataset-untimed-on-x86",
        "two-reference-releases",
    ],
)
def test_a_tree_the_page_cannot_be_written_from_is_refused_in_one_line(tmp_path, tree, message):
    with pytest.raises(SystemExit) as excinfo:
        home_page(home_tree(tmp_path, **tree))
    assert excinfo.value.code == f"ERROR: [render_home_page] {message}"


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (
            "Reads <<plugins-read>> plugins<<if plugins>> of the support table<<end plugins>>.\n",
            "home_page.template.md holds the place '<<plugins-read>>' outside its part, and this tree gives it no value",
        ),
        (
            "Up to <<no-such-place>> times.\n",
            "home_page.template.md holds the mark '<<no-such-place>>', which the script does not fill",
        ),
        (
            "Up to <<speed-up>> times<<if no-such-part>> or more<<end no-such-part>>.\n",
            "home_page.template.md names the part '<<if no-such-part>>', which the script does not know",
        ),
    ],
    ids=["place-outside-its-part", "unknown-place", "unknown-part"],
)
def test_a_template_mark_the_tree_gives_no_value_is_refused_in_one_line(tmp_path, text, message):
    with pytest.raises(SystemExit) as excinfo:
        home_page(home_tree(tmp_path), template=f"<<page>>\n{text}")
    assert excinfo.value.code == f"ERROR: [render_home_page] {message}"


def test_a_template_of_another_text_is_filled_and_the_site_tokens_pass_through(tmp_path):
    template = "<<page>>\nUp to <<speed-up>> times<<if plugins>> with plugins<<end plugins>> at {{vep-rs-ref}}.\n"
    assert home_page(home_tree(tmp_path), template=template) == "Up to 30 times at {{vep-rs-ref}}.\n"


def test_the_plugins_clause_states_the_plugins_the_released_tool_reads(tmp_path):
    page = home_page(home_tree(tmp_path, plugins=17))
    assert home_reasons(page)[1].endswith(", on four whole-genome datasets and [17 plugins](plugins.md#support-table).")


def test_a_tree_without_a_plugin_record_of_the_release_has_no_plugins_clause(tmp_path):
    without = home_page(home_tree(tmp_path / "without"))
    assert home_reasons(without)[1].endswith(", on four whole-genome datasets.")
    assert "plugins.md" not in without
    other = home_tree(tmp_path / "other", records={"first-plugins-v1.2.2.json": home_plugins(9, release="v1.2.2")})
    assert home_page(other) == without


def test_values_prints_each_filled_place_with_its_file_and_key_and_writes_nothing(capsys, tmp_path):
    repo = home_tree(tmp_path, plugins=17, values=home_values())
    before = files_of(repo)
    home.main(["--repo", str(repo), "--values"])
    rows = [line.split("\t") for line in capsys.readouterr().out.splitlines()]
    population, release, plugins = (f"{PROVENANCE}/{name}" for name in (POPULATION_RECORD, RELEASE_RECORD, PLUGIN_RECORD))
    values = HOME_VALUES
    seconds = "wall_time.{}.per_dataset[dataset=gnomad_v4.1_genomes, arch=arm64].sum_of_cell_medians_sec"
    roles = "tools, by the role of each id in comparators.tsv"
    assert rows[0] == ["place", "value", "file", "key", "rule"]
    assert all(len(row) == 5 for row in rows)
    assert {row[0]: row[1:4] for row in rows[1:]} == {
        "if plugins": ["holds", "", ""],
        "if figure": ["holds", "", ""],
        "if one-comparator": ["holds", "", ""],
        "if several-comparators": ["left out", "", ""],
        "if timed-apart": ["left out", "", ""],
        "speed-up": ["30", population, f"{seconds.format('ensembl_vep')} / {seconds.format('vep_rs')}"],
        "gnomad-label": ["Cohort A v7", population, "datasets[dataset=gnomad_v4.1_genomes].label"],
        "gnomad-variants": ["4,321,987", population, "datasets[dataset=gnomad_v4.1_genomes].variants"],
        "gnomad-arm-seconds": ["20", population, seconds.format("vep_rs")],
        "reference-release": ["100.1", release, "perl_reference.*.image"],
        "whole-genome-datasets": ["four", release, "concordance.s01, s02"],
        "plugins-read": ["17", plugins, "support[engine=vep-rs].read"],
        "figure-notes": [f"{HOME_NOTE} ", values, "notes, panels[].note"],
        "released": ["vep-rs v1.2.3", values, roles],
        "reference": ["Ensembl VEP 100.1", values, roles],
        "comparators": ["fastVEP v9.8.7", values, roles],
        "compared": ["vep-rs v1.2.3 and fastVEP v9.8.7", values, roles],
        "panel-count": ["Four panels", values, "panels"],
        "panels": [FOUR_PANELS, values, "panels[].id"],
    }
    assert files_of(repo) == before


def test_values_prints_no_row_for_a_place_the_tree_leaves_out(capsys, tmp_path):
    home.main(["--repo", str(home_tree(tmp_path)), "--values"])
    assert [line.split("\t")[0] for line in capsys.readouterr().out.splitlines()[1:]] == [
        "if plugins", "if figure", "if one-comparator", "if several-comparators", "if timed-apart",
        "speed-up", "gnomad-label", "gnomad-variants", "gnomad-arm-seconds", "reference-release", "whole-genome-datasets",
    ]


def test_two_runs_on_one_tree_write_the_same_bytes(tmp_path):
    repo = home_tree(tmp_path, plugins=17, values=home_values())
    page = repo / "docs" / "index.md"
    written = []
    for seed in ("1", "2"):
        subprocess.run(
            [sys.executable, str(HERE / "render_home_page.py"), "--repo", str(repo)],
            env={**os.environ, "PYTHONHASHSEED": seed},
            check=True,
            capture_output=True,
        )
        written.append(page.read_bytes())
        page.unlink()
    assert written[0] == written[1] == home_page(repo).encode("utf-8")


def test_check_refuses_a_page_with_other_line_endings(tmp_path):
    repo = home_tree(tmp_path)
    home.main(["--repo", str(repo)])
    page = repo / "docs" / "index.md"
    page.write_bytes(page.read_bytes().replace(b"\n", b"\r\n"))
    with pytest.raises(SystemExit) as excinfo:
        home.main(["--repo", str(repo), "--check"])
    assert excinfo.value.code == CHECK_REFUSAL.format(1)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (
            "<<if figure>>Shows <<no-such-place>>.<<end figure>>\n",
            "holds the mark '<<no-such-place>>', which the script does not fill",
        ),
        ("Up to <<speed-up>> times<<if plugins>> or more.\n", "holds '<<if plugins>>' without its pair"),
    ],
    ids=["unknown-place-in-a-part-left-out", "part-without-its-end"],
)
def test_a_template_mark_is_read_wherever_it_stands(tmp_path, text, message):
    with pytest.raises(SystemExit) as excinfo:
        home_page(home_tree(tmp_path), template=f"<<page>>\n{text}")
    assert excinfo.value.code == f"ERROR: [render_home_page] home_page.template.md {message}"


def test_a_note_that_is_not_one_line_is_refused(tmp_path):
    repo = home_tree(tmp_path, values=home_values(notes=("Two\nlines.",)))
    with pytest.raises(SystemExit) as excinfo:
        home_page(repo)
    assert excinfo.value.code == f"ERROR: [render_home_page] {HOME_VALUES} notes is not a list of lines"


@pytest.mark.parametrize("note", ["Two\nlines.", 7], ids=["two-lines", "a-number"])
def test_a_panel_note_that_is_not_one_line_is_refused(tmp_path, note):
    repo = home_tree(tmp_path, values=home_values(panel_notes={"c": note}))
    with pytest.raises(SystemExit) as excinfo:
        home_page(repo)
    assert excinfo.value.code == f"ERROR: [render_home_page] {HOME_VALUES} panels[2] note is not a line"


def redirect_filter():
    """The redirect filter of `build_release_site.sh`: the program the script feeds its third Python run."""
    script = (HERE / "build_release_site.sh").read_text(encoding="utf-8")
    body = script[script.index("\n", script.index('"$docs_ref" "$site_url" <<\'EOF\'')) + 1:]
    return body[:body.index("\nEOF\n") + 1]


def filter_redirects(tmp_path, built, table):
    """Run the redirect filter over `table` on a site whose pages, each with its Markdown source, carry the
    ids `built` gives them; returns the finished run and the table as the filter left it."""
    docs, site = tmp_path / "docs", tmp_path / "site"
    for page, ids in built.items():
        headings = "".join(f'<h2 id="{anchor}">Heading</h2>' for anchor in ids)
        for root, suffix, text in ((docs, ".md", "# Page\n"), (site, ".html", headings)):
            path = root / f"{page}{suffix}"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
    toml = tmp_path / "book.toml"
    toml.write_text(table, encoding="utf-8")
    done = subprocess.run(
        [sys.executable, "-", str(toml), str(docs), str(site), RELEASE, "/vep-rs/"],
        input=redirect_filter(), capture_output=True, text=True, check=False,
    )
    return done, toml.read_text(encoding="utf-8")


def test_a_redirect_to_an_anchor_the_built_page_lacks_is_sent_to_the_page(tmp_path):
    table = (
        "[output.html.redirect]\n"
        f'"/index.html#first" = "{SUB}/setup.html#step-1"\n'
        f'"/index.html#second" = "{SUB}/setup.html#step-2"\n'
    )
    done, left = filter_redirects(tmp_path, {"index": ["page"], f"{SUB}/setup": ["step-1"]}, table)
    assert done.returncode == 0, done.stderr
    assert done.stdout == (
        f"redirect sent to the page, anchor not on the {RELEASE} page: /index.html#second -> {SUB}/setup.html#step-2\n"
        "redirects kept: 2, dropped: 0, sent to the page: 1, sent to the source's anchor: 0\n"
    )
    assert left == table.replace("setup.html#step-2", "setup.html")


def test_a_redirect_whose_destination_page_carries_the_sources_anchor_name_is_sent_to_that_anchor(tmp_path):
    """The setup page carries `step-1` and `step-2`. Only a redirect whose destination anchor the page lacks
    moves: to the source's anchor name where the page carries it, else to the page."""
    built = {"index": ["page"], f"{SUB}/page": ["page"], f"{SUB}/setup": ["step-1", "step-2"]}
    table = (
        "[output.html.redirect]\n"
        f'"/index.html#step-1" = "{SUB}/setup.html#first-step"\n'
        f'"/index.html#step-2" = "{SUB}/setup.html#step-1"\n'
        f'"/index.html#step-3" = "{SUB}/setup.html#third-step"\n'
        f'"/moved.html" = "{SUB}/setup.html#first-step"\n'
        f'"/{SUB}/page.html#step-2" = "setup.html#second-step"\n'
    )
    done, left = filter_redirects(tmp_path, built, table)
    sent = f"redirect sent to the source's anchor, the destination's anchor not on the {RELEASE} page"
    to_page = f"redirect sent to the page, anchor not on the {RELEASE} page"
    assert done.returncode == 0, done.stderr
    assert done.stdout == (
        f"{sent}: /index.html#step-1 -> {SUB}/setup.html#first-step\n"
        f"{to_page}: /index.html#step-3 -> {SUB}/setup.html#third-step\n"
        f"{to_page}: /moved.html -> {SUB}/setup.html#first-step\n"
        f"{sent}: /{SUB}/page.html#step-2 -> setup.html#second-step\n"
        "redirects kept: 5, dropped: 0, sent to the page: 2, sent to the source's anchor: 2\n"
    )
    assert left == (
        "[output.html.redirect]\n"
        f'"/index.html#step-1" = "{SUB}/setup.html#step-1"\n'
        f'"/index.html#step-2" = "{SUB}/setup.html#step-1"\n'
        f'"/index.html#step-3" = "{SUB}/setup.html"\n'
        f'"/moved.html" = "{SUB}/setup.html"\n'
        f'"/{SUB}/page.html#step-2" = "setup.html#step-2"\n'
    )


THEME_SCRIPT = ROOT / "theme" / "vep-rs.js"
SITE = "https://docs.example"
ATTACKER = "attacker.example"

# The theme script runs whole under node against this stand-in for the DOM calls it makes on a page with nothing to
# fit: one content element with its main, the menu bar and its buttons, and elements that keep their children and
# listeners. The select's value answers as a browser's does, the chosen option's own value or else its text, unless the
# scenario plants one, so a script that navigated to text read back from the page would be seen doing it. Each choice
# reports the addresses the script navigated to and the requests it made; the report also carries the address the
# preview banner links, when the scenario names the page's channel and the script shows the banner.
SELECTOR_HARNESS = r"""
const fs = require('fs');
const vm = require('vm');
const [script, scenarioText] = process.argv.slice(2);
const scenario = JSON.parse(scenarioText);

function element(tag) {
    const listeners = {};
    return {
        tagName: String(tag).toUpperCase(),
        children: [],
        className: '',
        textContent: '',
        style: { setProperty() {}, removeProperty() {} },
        classList: { add() {}, remove() {}, toggle() {} },
        setAttribute() {},
        removeAttribute() {},
        getAttribute() { return null; },
        hasAttribute() { return false; },
        querySelector() { return null; },
        querySelectorAll() { return []; },
        appendChild(child) { this.children.push(child); return child; },
        append(...nodes) { this.children.push(...nodes); },
        insertBefore(node) { this.children.unshift(node); return node; },
        get firstChild() { return this.children[0] || null; },
        addEventListener(type, listener) { (listeners[type] = listeners[type] || []).push(listener); },
        fire(type) { (listeners[type] || []).forEach(listener => listener({ type, target: this })); },
        focus() {},
        blur() {},
    };
}

function select() {
    const node = element('select');
    let planted;
    node.selectedIndex = -1;
    Object.defineProperty(node, 'value', {
        get() {
            if (planted !== undefined) {
                return planted;
            }
            const option = node.children[node.selectedIndex];
            if (!option) {
                return '';
            }
            return option.value === undefined ? option.textContent : option.value;
        },
        set(text) { planted = text; },
    });
    return node;
}

const main = element('main');
const content = element('div');
content.querySelector = selector => (selector === 'main' ? main : null);
const menuBar = element('div');
const buttons = element('div');
const document = {
    documentElement: element('html'),
    body: element('body'),
    activeElement: null,
    getElementById: id => ({ 'mdbook-content': content, 'mdbook-menu-bar': menuBar })[id] || null,
    querySelector(selector) {
        if (selector === 'meta[name="vep-rs-channel"]' && scenario.channel) {
            return { getAttribute: () => scenario.channel };
        }
        return selector === '#mdbook-menu-bar .right-buttons' ? buttons : null;
    },
    querySelectorAll: () => [],
    createElement: tag => (tag === 'select' ? select() : element(tag)),
    addEventListener() {},
};

const navigations = [];
const requests = [];
const location = {
    origin: scenario.origin,
    pathname: scenario.pathname,
    search: '',
    hash: scenario.hash || '',
    get href() { return this.origin + this.pathname + this.hash; },
    set href(address) { navigations.push(String(address)); },
    assign(address) { navigations.push(String(address)); },
    replace(address) { navigations.push(String(address)); },
};
const window = {
    location,
    matchMedia: () => ({ matches: false, addEventListener() {} }),
    addEventListener() {},
    open(address) { navigations.push(String(address)); },
};

function fetch(address, options) {
    const text = String(address);
    if (text.endsWith('versions.json')) {
        return Promise.resolve({ ok: true, json: () => Promise.resolve(scenario.versions) });
    }
    requests.push(((options && options.method) || 'GET') + ' ' + text);
    if (scenario.offline) {
        return Promise.reject(new Error('offline'));
    }
    return Promise.resolve({ ok: scenario.pages.includes(text) });
}

const settle = () => new Promise(resolve => setImmediate(resolve));

(async () => {
    vm.runInNewContext(fs.readFileSync(script, 'utf8'), {
        window,
        document,
        fetch,
        URL,
        localStorage: { getItem: () => null },
        getComputedStyle: () => ({}),
        requestAnimationFrame: callback => callback(),
    });
    await settle();
    const label = buttons.children[0];
    const chooser = label ? label.children.find(child => child.tagName === 'SELECT') : null;
    const note = content.children.find(child => child.className === 'vep-banner');
    const link = note ? note.children.find(child => child.tagName === 'A') : null;
    const report = {
        options: chooser ? chooser.children.map(option => option.textContent) : null,
        banner: link ? link.href : null,
        choices: [],
    };
    for (const choice of scenario.choices) {
        navigations.length = 0;
        requests.length = 0;
        chooser.selectedIndex = choice.index;
        if (choice.planted !== undefined) {
            chooser.value = choice.planted;
        }
        chooser.fire('change');
        await settle();
        await settle();
        report.choices.push({ navigations: [...navigations], requests: [...requests] });
    }
    process.stdout.write(JSON.stringify(report));
})().catch(error => {
    process.stderr.write(String(error && error.stack ? error.stack : error) + '\n');
    process.exit(1);
});
"""

# A list as the site writes it, the release at the site's path and the preview under it, then entries no site tool
# writes: a script address, another site by its full address, by a scheme-relative one and by a backslash a browser
# reads as a slash, and a path that resolves outside every listed one.
LISTED = [{"name": RELEASE, "path": "/vep-rs/"}, {"name": "main", "path": "/vep-rs/main/"}]
UNLISTED = [
    {"name": "script", "path": "javascript:alert(document.domain)"},
    {"name": "elsewhere", "path": f"https://{ATTACKER}/vep-rs/"},
    {"name": "scheme-relative", "path": f"//{ATTACKER}/vep-rs/"},
    {"name": "backslash", "path": f"/\\{ATTACKER}/vep-rs/"},
    {"name": "another scheme", "path": "ftp://docs.example/vep-rs/"},
    {"name": "outside", "path": "/vep-rs/main/../../private/"},
]
PREVIEW_PAGE = f"/vep-rs/main/{GUIDE_PAGE[:-3]}.html"
RELEASE_PAGE = f"{SITE}/vep-rs/{GUIDE_PAGE[:-3]}.html"


def choose(
    tmp_path, choices, *, script=THEME_SCRIPT, channels=None, pathname=PREVIEW_PAGE, pages=(), offline=False, channel=None,
):
    """What the theme script at `script` does for each choice in the version selector of the page at `pathname`: the
    option names it listed, the address its preview banner links on a page of `channel`, and per choice the addresses
    it navigated to and the requests it made."""
    node = shutil.which("node")
    if node is None:
        pytest.skip("no node on PATH: the version selector's behaviour is tested by running theme/vep-rs.js under node")
    harness = tmp_path / "selector_harness.js"
    harness.write_text(SELECTOR_HARNESS, encoding="utf-8")
    scenario = {
        "origin": SITE,
        "pathname": pathname,
        "hash": "#section",
        "channel": channel,
        "versions": {"release": RELEASE, "channels": LISTED + UNLISTED if channels is None else channels},
        "pages": list(pages),
        "offline": offline,
        "choices": choices,
    }
    done = subprocess.run(
        [node, str(harness), str(script), json.dumps(scenario)], capture_output=True, text=True, check=False,
    )
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


def test_the_version_selector_lists_every_entry_in_the_lists_order(tmp_path):
    report = choose(tmp_path, [])
    assert report["options"] == [entry["name"] for entry in LISTED + UNLISTED]


def test_choosing_a_version_opens_the_same_page_there_with_the_fragment(tmp_path):
    report = choose(tmp_path, [{"index": 0}], pages=[RELEASE_PAGE])
    assert report["choices"] == [{"navigations": [f"{RELEASE_PAGE}#section"], "requests": [f"HEAD {RELEASE_PAGE}"]}]


@pytest.mark.parametrize("offline", [False, True], ids=["the page is not there", "the request fails"])
def test_choosing_a_version_without_the_page_opens_the_versions_home(tmp_path, offline):
    report = choose(tmp_path, [{"index": 0}], offline=offline)
    assert report["choices"] == [{"navigations": [f"{SITE}/vep-rs/"], "requests": [f"HEAD {RELEASE_PAGE}"]}]


def test_a_site_served_at_an_unlisted_path_opens_the_chosen_versions_home(tmp_path):
    report = choose(tmp_path, [{"index": 1}, {"index": 2}], pathname="/local-copy/index.html", pages=[RELEASE_PAGE])
    assert report["choices"] == [
        {"navigations": [f"{SITE}/vep-rs/main/"], "requests": []},
        {"navigations": [], "requests": []},
    ]


@pytest.mark.parametrize("entry", UNLISTED, ids=[entry["name"] for entry in UNLISTED])
def test_an_entry_off_the_site_or_off_the_list_is_never_navigated_to(tmp_path, entry):
    index = (LISTED + UNLISTED).index(entry)
    report = choose(tmp_path, [{"index": index}, {"index": 0}], pages=[RELEASE_PAGE])
    assert report["choices"] == [
        {"navigations": [], "requests": []},
        {"navigations": [f"{RELEASE_PAGE}#section"], "requests": [f"HEAD {RELEASE_PAGE}"]},
    ]


@pytest.mark.parametrize(
    "choice, navigations",
    [
        ({"index": 0, "planted": "javascript:alert(document.domain)"}, [f"{RELEASE_PAGE}#section"]),
        ({"index": 0, "planted": f"https://{ATTACKER}/"}, [f"{RELEASE_PAGE}#section"]),
        ({"index": -1, "planted": f"https://{ATTACKER}/"}, []),
        ({"index": len(LISTED + UNLISTED), "planted": "/vep-rs/"}, []),
    ],
    ids=["a script address", "another site", "no option chosen", "an option past the list"],
)
def test_the_destination_is_the_lists_entry_whatever_the_pages_select_holds(tmp_path, choice, navigations):
    report = choose(tmp_path, [choice], pages=[RELEASE_PAGE])
    assert report["choices"][0]["navigations"] == navigations


def test_the_preview_banner_links_the_listed_release(tmp_path):
    assert choose(tmp_path, [], channel="main")["banner"] == f"{SITE}/vep-rs/"
    assert choose(tmp_path, [])["banner"] is None


@pytest.mark.parametrize("entry", UNLISTED, ids=[entry["name"] for entry in UNLISTED])
def test_the_preview_banner_links_the_parent_site_when_the_release_path_is_off_the_site_or_the_list(tmp_path, entry):
    channels = [{"name": RELEASE, "path": entry["path"]}, LISTED[1]]
    assert choose(tmp_path, [], channel="main", channels=channels)["banner"] == "../"


def test_the_theme_script_navigates_in_one_place_and_never_to_a_value_read_from_the_page():
    text = THEME_SCRIPT.read_text(encoding="utf-8")
    assert len(re.findall(r"location\.href\s*=(?!=)", text)) == 1
    assert not re.search(r"location\.(href\s*=(?!=)|assign\(|replace\().*\.value\b", text)
    assert not re.search(r"\bselect\.value\b", text)
