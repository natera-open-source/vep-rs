"""Tests for the site tooling in this directory, on fixtures the tests write: the mdBook preprocessor driven
through its stdin protocol, the page metrics and their gate, the post-build finish with its three checks, the
benchmark chart renderer, the home page's writer, and the screenshot script's pure helpers. No test runs mdbook or a
browser."""

from __future__ import annotations

import io
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import mdbook_repo_links as links
import render_benchmark_chart as chart
import render_home_page as home
import render_release_figure as figure
import site_metrics as metrics
import site_postbuild as postbuild
import site_screenshots as shots
import test_render_release_figure as figure_records

REPO_URL = "https://github.com/natera-open-source/vep-rs"
EDIT_URL = f"{REPO_URL}/edit/main/"
RELEASE = "v1.2.3"
VERSION = "1.2.3"
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


def built_site(tmp_path, chart_release=RELEASE, redirect="guide/setup.html#step-2", book=f"vep-rs {RELEASE}"):
    """A built site of two reader pages, each with a title and an opening paragraph of its own, two redirect
    stubs, the two pages mdBook generates for printing and not-found, and the benchmark chart; beside it the
    source tree whose README names the release."""
    site = tmp_path / "book"
    (site / "guide").mkdir(parents=True)
    (site / "assets").mkdir()
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
    (site / "assets" / "benchmark.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" class="benchmark-chart" data-release="{chart_release}">\n'
        f"<desc>Release {chart_release}</desc>\n</svg>\n",
        encoding="utf-8",
    )
    source = tmp_path / "source"
    source.mkdir()
    (source / "README.md").write_text(f"# Project\n\n## Current release: {RELEASE}\n\nText.\n", encoding="utf-8")
    return site, source


def finish(monkeypatch, capsys, site, *extra, base=BASE, readme_home=None):
    """Exit status and output of the post-build step on `site`."""
    set_variable(monkeypatch, "VEP_RS_DOCS_README_HOME", readme_home)
    monkeypatch.setattr(sys, "argv", ["site_postbuild.py", str(site), "--base", base, "--release", RELEASE, *extra])
    with pytest.raises(SystemExit) as excinfo:
        postbuild.main()
    return excinfo.value.code, capsys.readouterr().out


def test_canonical_and_og_url_are_inserted_once_and_a_second_run_adds_nothing(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
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
    assert "benchmark chart: v1.2.3 matches the README\n" in out
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    assert (site / "index.html").read_text(encoding="utf-8") == index
    assert "inserted on 0, 2 already carried a canonical link\n" in out


def test_menu_title_is_demoted_and_a_title_repeating_the_book_name_collapses(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    finish(monkeypatch, capsys, site, "--source", str(source))
    index = (site / "index.html").read_text(encoding="utf-8")
    assert postbuild.HOME_TAGLINE == "fast Ensembl VEP (Variant Effect Predictor) in Rust"
    assert f'<p class="menu-title">vep-rs {RELEASE}</p>' in index
    assert "<h1 class=" not in index
    assert f"<title>vep-rs: {postbuild.HOME_TAGLINE}</title>" in index
    assert "<title>Setup - vep-rs documentation</title>" in (site / "guide" / "setup.html").read_text(encoding="utf-8")
    assert "<title>Not found - vep-rs</title>" in (site / "404.html").read_text(encoding="utf-8")


def test_preview_channel_marks_its_pages_and_writes_no_sitemap(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, "--source", str(source), "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f'<link rel="canonical" href="{BASE}main/">' in index
    assert '<meta name="vep-rs-channel" content="main">' in index
    assert '<meta name="robots" content="noindex">' in index
    assert not (site / "sitemap.xml").exists()
    assert json.loads((site / "versions.json").read_text(encoding="utf-8")) == VERSIONS


def test_release_form_adds_no_channel_metas_and_writes_versions_and_the_sitemap(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
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
    site, source = built_site(tmp_path)
    (site / "index.html").write_text(
        described_page("vep-rs", f"<main>\n<p>Documentation: short.</p>\n<p>{PROSE}</p>\n</main>"), encoding="utf-8"
    )
    (site / "guide" / "setup.html").write_text(
        described_page("Setup", "<main>\n<h2 id=\"step-2\">Step 2</h2>\n<p>Too short to describe anything.</p>\n</main>"), encoding="utf-8"
    )
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
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
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    assert (site / "index.html").read_text(encoding="utf-8") == index


def test_preview_home_gets_the_tagline_but_no_json_ld(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    (site / "index.html").write_text(
        described_page("vep-rs", f"<main><p>{PROSE}</p></main>", "vep-rs main"), encoding="utf-8"
    )
    code, out = finish(monkeypatch, capsys, site, "--source", str(source), "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f"<title>vep-rs main: {postbuild.HOME_TAGLINE}</title>" in index
    assert '<meta property="og:site_name" content="vep-rs main">' in index
    assert "application/ld+json" not in index


def test_page_heads_line_counts_the_pages_a_reader_lands_on(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    assert "page heads: 2 pages, 0 titles shared, 0 descriptions shared, 0 on the book's description\n" in out
    assert "WARNING" not in out


def third_page_finding(monkeypatch, capsys, tmp_path, title, prose, counts, finding):
    """Holds the post-build step, on each channel, to its finding over the fixture site plus `guide/other.html`:
    a warning and exit status 0 on the release, a failure on the preview. `{name}` in `finding` is what a short
    title ends with on the channel."""
    for channel, name in ((None, "vep-rs documentation"), ("main", "vep-rs main")):
        book = name if channel else f"vep-rs {RELEASE}"
        site, source = built_site(tmp_path / (channel or "release"), book=book)
        (site / "guide" / "other.html").write_text(page_html(title, f"<p>{prose}</p>", book), encoding="utf-8")
        extra = ("--channel", channel) if channel else ()
        base = BASE + (f"{channel}/" if channel else "")
        code, out = finish(monkeypatch, capsys, site, "--source", str(source), *extra, base=base)
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
    site, source = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
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
    site, source = built_site(tmp_path)
    kept, bare, cut = add_chapters(site, ["Setting up a cache", "x" * 52, CUT_A])
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
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
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    assert {path: path.read_bytes() for path in site.rglob("*") if path.is_file()} == finished
    assert TITLES_LINE.format(5, 1, 1, 2, 1, 0) in out


def test_the_preview_keeps_every_title_whole_and_prints_no_titles_line(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path, book="vep-rs main")
    names = ["A", "Setting up a cache", "x" * 52, CUT_A]
    pages = add_chapters(site, names, "vep-rs main")
    code, out = finish(monkeypatch, capsys, site, "--source", str(source), "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    assert "page titles:" not in out
    assert [title_of(site, page) for page in pages] == [f"{name} - vep-rs main" for name in names]
    assert title_of(site, "guide/setup.html") == "Setup - vep-rs main"


def test_a_title_is_measured_unescaped_and_written_back_escaped(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    kept, lengthened = add_chapters(site, ["R&D " + "x" * 47, "A <B> & C's"])
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    assert title_of(site, kept) == "R&amp;D " + "x" * 47 + ENDING
    assert title_of(site, lengthened) == "A &lt;B&gt; &amp; C's" + ENDING_LONG
    assert TITLES_LINE.format(4, 1, 2, 1, 0, 0) in out


def test_a_title_under_25_characters_fails_the_release_and_passes_the_preview(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path / "release")
    (page,) = add_chapters(site, ["A"])
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 1
    assert f"ERROR: [site_postbuild] {page}: title 'A - vep-rs documentation' is 24 characters, under 25\n" in out
    assert TITLES_LINE.format(3, 0, 2, 1, 0, 1) in out
    site, source = built_site(tmp_path / "main", book="vep-rs main")
    add_chapters(site, ["A"], "vep-rs main")
    code, out = finish(monkeypatch, capsys, site, "--source", str(source), "--channel", "main", base=BASE + "main/")
    assert code == 0, out


def test_a_title_over_60_characters_fails_the_release(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(postbuild, "HOME_TAGLINE", postbuild.HOME_TAGLINE + ", every flag")
    site, source = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 1
    assert f"ERROR: [site_postbuild] index.html: title '{HOME_PAGE_TITLE}, every flag' is 71 characters, over 60\n" in out
    assert TITLES_LINE.format(2, 0, 1, 1, 0, 1) in out


def test_two_headings_cut_to_one_title_fail_the_release_and_pass_the_preview(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path / "release")
    first, second = add_chapters(site, [CUT_A, CUT_B])
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 1
    assert "page heads: 4 pages, 1 titles shared, 0 descriptions shared, 0 on the book's description\n" in out
    assert f"ERROR: [site_postbuild] title {CUT_TITLE!r} is the fit of 2 headings: {first}, {second}\n" in out
    assert "WARNING" not in out
    site, source = built_site(tmp_path / "main", book="vep-rs main")
    add_chapters(site, [CUT_A, CUT_B], "vep-rs main")
    code, out = finish(monkeypatch, capsys, site, "--source", str(source), "--channel", "main", base=BASE + "main/")
    assert code == 0, out
    assert "page heads: 4 pages, 0 titles shared, 0 descriptions shared, 0 on the book's description\n" in out


def test_one_heading_cut_on_two_pages_is_a_warning_on_the_release(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    first, second = add_chapters(site, [CUT_A, CUT_A])
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    assert f"WARNING: [site_postbuild] title {CUT_TITLE!r} is on 2 pages: {first}, {second}\n" in out
    assert "is the fit of" not in out


def test_a_chapter_name_that_ends_in_the_release_tag_keeps_it_on_a_second_run(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    kept, bare = add_chapters(site, [f"Upgrading to {RELEASE}", f"{'x' * 45} {RELEASE}"])
    for _run in range(2):
        code, out = finish(monkeypatch, capsys, site, "--source", str(source))
        assert code == 0, out
        assert title_of(site, kept) == f"Upgrading to {RELEASE}" + ENDING
        assert title_of(site, bare) == f"{'x' * 45} {RELEASE}"


def test_a_home_chapter_of_another_name_is_titled_with_the_book_name_alone(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path)
    (site / "index.html").write_text(described_page("Overview", f"<p>{HOME_PROSE}</p>"), encoding="utf-8")
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
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
    site, source = built_site(tmp_path)
    for page in (site / "index.html", site / "guide" / "setup.html"):
        framed = page.read_text(encoding="utf-8").replace("<body>\n", f"<body>\n<nav>{FRAME}</nav>\n")
        page.write_text(framed, encoding="utf-8")
    (site / "toc.html").write_text(
        f'<!DOCTYPE HTML>\n<html>\n<head></head>\n<body class="sidebar-iframe-inner">\n{CHAPTERS}\n</body>\n</html>\n',
        encoding="utf-8",
    )
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f'<nav><noscript>\n<div class="sidebar-scrollbox">{CHAPTERS}</div>\n</noscript></nav>' in index
    assert "<iframe" not in index
    assert f"<nav>{FRAME}</nav>" in (site / "guide" / "setup.html").read_text(encoding="utf-8")
    assert "home page: 2 chapter links listed in its <noscript> block\n" in out
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
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
    site, source = built_site(tmp_path, redirect=destination)
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 1
    assert f"ERROR: [site_postbuild] {message}\n" in out


@pytest.mark.parametrize(
    "destination", ["/vep-rs", "/vep-rs/", "guide/"], ids=["site-root-without-the-slash", "site-root", "directory"]
)
def test_a_redirect_to_a_directory_is_checked_against_its_index_page(monkeypatch, capsys, tmp_path, destination):
    site, source = built_site(tmp_path, redirect=destination)
    (site / "guide" / "index.html").write_text(page_html("Guide", '<h1 id="guide">Guide</h1>'), encoding="utf-8")
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 0, out
    assert "redirect destinations: 3 checked, 1 external not fetched, 0 broken\n" in out


def test_chart_release_differing_from_the_readme_heading_fails(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path, chart_release="v1.2.2")
    code, out = finish(monkeypatch, capsys, site, "--source", str(source))
    assert code == 1
    assert "ERROR: [site_postbuild] benchmark chart says v1.2.2, the README's Current release heading says v1.2.3\n" in out


def test_chart_without_a_source_tree_cannot_be_checked(monkeypatch, capsys, tmp_path):
    site, _source = built_site(tmp_path)
    code, out = finish(monkeypatch, capsys, site)
    assert code == 1
    assert "ERROR: [site_postbuild] assets/benchmark.svg says v1.2.3; --source is needed to compare it with the README\n" in out


def test_readme_home_skips_the_chart_check_and_points_the_edit_link_at_the_readme(monkeypatch, capsys, tmp_path):
    site, source = built_site(tmp_path, chart_release="v1.2.2")
    code, out = finish(monkeypatch, capsys, site, "--source", str(source), readme_home="1")
    assert code == 0, out
    assert "benchmark chart: the README is the home page; release check skipped\n" in out
    index = (site / "index.html").read_text(encoding="utf-8")
    assert f'<a href="{EDIT_URL}README.md">Edit</a>' in index
    assert "docs/index.md" not in index


DATASETS = [
    ("alpha", "Alpha, whole genome", "GRCh38"),
    ("beta", "Beta, whole genome", "GRCh37"),
    ("gamma", "Gamma, whole genome", "GRCh37"),
    ("delta", "Delta, whole genome", "GRCh38"),
]
# Seconds per chromosome on ARM; every dataset has chromosomes 1 and 2, so a dataset's time is twice its entry.
MEDIANS = {
    "vep_rs": {"alpha": 10, "beta": 20, "gamma": 30, "delta": 40},
    "fastvep": {"alpha": 100, "beta": 200, "gamma": 300, "delta": 600},
    "ensembl_vep": {"alpha": 2000, "beta": 4000, "gamma": 6000, "delta": 8000},
}
ENGINE_VERSIONS = {"vep_rs": "v1.2.3", "fastvep": "v2.0.0", "ensembl_vep": "release_115.2"}


def chart_record(engines=("vep_rs", "fastvep", "ensembl_vep"), release=RELEASE):
    """A population record of four two-chromosome datasets; the x86 cells carry other times, so a sum that
    reads them is wrong."""
    wall_time = {}
    for engine in engines:
        cells = []
        for key, _label, _assembly in DATASETS:
            for chromosome in ("1", "2"):
                for arch, extra in (("arm64", 0), ("x86_64", 1)):
                    cells.append({"dataset": key, "arch": arch, "chromosome": chromosome,
                                  "median_sec": MEDIANS[engine][key] + extra})
        wall_time[engine] = {"engine_version": ENGINE_VERSIONS[engine], "n_clones_per_cell": 5, "per_cell": cells}
    return {
        "release": release,
        "fork": 16,
        "datasets": [
            {"dataset": key, "label": label, "assembly": assembly, "chromosomes": ["1", "2"]}
            for key, label, assembly in DATASETS
        ],
        "wall_time": wall_time,
    }


def chart_repo(tmp_path, records):
    repo = tmp_path / "repo"
    provenance = repo / "docs" / "concordance-provenance"
    provenance.mkdir(parents=True)
    (repo / "README.md").write_text(f"# Project\n\n## Current release: {RELEASE}\n", encoding="utf-8")
    for name, record in records.items():
        (provenance / name).write_text(json.dumps(record), encoding="utf-8")
    return repo


def render_chart(monkeypatch, repo, out):
    monkeypatch.setattr(sys, "argv", ["render_benchmark_chart.py", "--repo", str(repo), "--out", str(out)])
    chart.main()


def test_population_record_is_the_one_file_of_the_version(tmp_path):
    repo = chart_repo(tmp_path, {
        "first-population-v1.2.3.json": chart_record(),
        "other-population-v1.2.2.json": chart_record(release="v1.2.2"),
    })
    assert chart.population_record(repo, VERSION) == repo / "docs" / "concordance-provenance" / "first-population-v1.2.3.json"


def test_no_population_record_fails_by_pattern(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        chart.population_record(chart_repo(tmp_path, {}), VERSION)
    assert excinfo.value.code == "ERROR: [render_benchmark_chart] no docs/concordance-provenance/*-population-v1.2.3.json"


def test_two_population_records_fail_naming_both(tmp_path):
    repo = chart_repo(tmp_path, {
        "first-population-v1.2.3.json": chart_record(),
        "second-population-v1.2.3.json": chart_record(),
    })
    with pytest.raises(SystemExit) as excinfo:
        chart.population_record(repo, VERSION)
    assert excinfo.value.code == (
        "ERROR: [render_benchmark_chart] 2 records match docs/concordance-provenance/*-population-v1.2.3.json: "
        "first-population-v1.2.3.json, second-population-v1.2.3.json"
    )


def test_load_record_refuses_a_release_other_than_the_readme_version(tmp_path):
    path = tmp_path / "first-population-v1.2.3.json"
    path.write_text(json.dumps(chart_record(release="v1.2.2")), encoding="utf-8")
    with pytest.raises(SystemExit) as excinfo:
        chart.load_record(path, VERSION)
    assert excinfo.value.code == (
        "ERROR: [render_benchmark_chart] first-population-v1.2.3.json records release 'v1.2.2'; the README names v1.2.3"
    )


def test_the_committed_chart_is_the_render_of_the_current_record(monkeypatch, capsys, tmp_path):
    out = tmp_path / "benchmark.svg"
    render_chart(monkeypatch, ROOT, out)
    assert out.read_bytes() == (ROOT / "docs" / "assets" / "benchmark.svg").read_bytes()
    version = chart.current_release(ROOT / "README.md")
    record = chart.population_record(ROOT, version)
    assert capsys.readouterr().out.endswith(f" bytes) for v{version} from {record.name}\n")


def test_a_record_without_a_required_engine_fails_by_name():
    with pytest.raises(SystemExit) as excinfo:
        chart.read_chart(chart_record(engines=("vep_rs", "fastvep")), "record.json", VERSION)
    assert excinfo.value.code == "ERROR: [render_benchmark_chart] record.json carries no wall_time.ensembl_vep"


def test_a_record_without_fastvep_omits_its_series_with_a_note(monkeypatch, capsys, tmp_path):
    repo = chart_repo(tmp_path, {"first-population-v1.2.3.json": chart_record(engines=("vep_rs", "ensembl_vep"))})
    out = tmp_path / "chart.svg"
    render_chart(monkeypatch, repo, out)
    printed = capsys.readouterr().out
    assert printed.startswith("note: first-population-v1.2.3.json has no wall_time.fastvep; the chart omits fastVEP\n")
    svg = out.read_text(encoding="utf-8")
    assert "bar-fastvep" not in svg
    assert "bar-ensembl-vep" in svg
    assert "bar-vep-rs" in svg


def test_svg_carries_the_release_stamp_and_no_blank_line():
    svg = chart.render(chart.read_chart(chart_record(), "record.json", VERSION))
    first = svg.split("\n", 1)[0]
    assert first.startswith("<svg ")
    assert first.endswith(f'data-release="{RELEASE}">')
    assert (
        f"<desc>Release {RELEASE}: whole-genome wall time per dataset on ARM for vep-rs v1.2.3, fastVEP v2.0.0, "
        "Ensembl VEP 115.2, from the release's population record.</desc>"
    ) in svg
    assert "\n\n" not in svg
    assert svg.endswith("</svg>\n")


def test_headline_states_the_speed_up_range_against_each_other_series():
    assert chart.headline(chart.read_chart(chart_record(), "record.json", VERSION)) == (
        "vep-rs annotates the four whole-genome datasets 10.0 to 15.0 times faster than fastVEP v2.0.0 and "
        "200 times faster than Ensembl VEP 115.2 on ARM."
    )


def test_title_names_the_base_machine_count_and_every_series_measured_on_another():
    assert chart.title_text(chart.read_chart(chart_record(), "record.json", VERSION)) == (
        "Wall time, sum of per-chromosome medians over 5 machines at 16 threads, ARM"
    )
    record = chart_record()
    record["wall_time"]["ensembl_vep"]["n_clones_per_cell"] = 1
    assert chart.title_text(chart.read_chart(record, "record.json", VERSION)) == (
        "Wall time, sum of per-chromosome medians over 5 machines (Ensembl VEP over 1) at 16 threads, ARM"
    )


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
PANEL_A = "raw F1 per dataset"
PANEL_B = "wall time per dataset and architecture"
PANEL_C = "raw F1 per plugin"
PANEL_D = "the plugins each tool reads"
FOUR_PANELS = f"{PANEL_A}, {PANEL_B}, {PANEL_C} and {PANEL_D}"
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
    version), its notes, its panel ids and, for each panel id of `panel_notes`, that panel's note, beside the
    caption and the alternative text of the README, which the page does not copy."""
    panel_notes = panel_notes or {}
    return {
        "image": image,
        "canvas": {"width_px": 100, "height_px": 125},
        "alt_text": "The alternative text the README shows.",
        "caption": "The caption the README shows.",
        "notes": list(notes),
        "sources": [],
        "tools": [{"id": tool_id, "name": name, "version": version, "colour": "#000000"} for tool_id, name, version in tools],
        "records": [],
        "panels": [
            {"id": panel, "title": f"Panel {panel}", "groups": []} | ({"note": panel_notes[panel]} if panel in panel_notes else {})
            for panel in panels
        ],
    }


def figure_tree(tmp_path):
    """A made-up tree of one release whose records the figure script and the page's script both read, with the
    figure's pair: the values file as the figure script's values layer writes it from those records, and made-up
    image bytes beside it. Returns the tree and the values file's document."""
    records = json.loads(json.dumps(figure_records.record_set()).replace('"pop_a"', '"gnomad_v4.1_genomes"'))
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
            {"notes": ("On *one* machine [made up].",), "tools": (RELEASED_TOOL, ("fastvep", "fast_VEP*", "v9.8.7"), REFERENCE_TOOL)},
            (
                "On \\*one\\* machine \\[made up\\]. Four panels comparing vep-rs v1.2.3 and fast\\_VEP\\* v9.8.7 against "
                f"Ensembl VEP 100.1: {FOUR_PANELS}."
            ),
            "On \\*one\\* machine \\[made up\\]. vep-rs v1.2.3 against Ensembl VEP 100.1, with fast\\_VEP\\* v9.8.7 as comparator.",
        ),
    ],
    ids=[
        "one-comparator", "no-comparator", "no-notes", "two-notes", "panel-notes-after-the-notes", "one-panel",
        "panels-in-the-file-order", "markup-escaped",
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
            {"values": home_values(panels="ae")},
            f"{HOME_VALUES} holds the panel e, and home_page.template.md has no '<<panel e>>' clause",
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
        "if plugins", "if figure", "if one-comparator", "if several-comparators",
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


def test_a_redirect_to_an_anchor_the_built_page_lacks_is_sent_to_the_page(tmp_path):
    docs, site = tmp_path / "docs", tmp_path / "site"
    (docs / SUB).mkdir(parents=True)
    (site / SUB).mkdir(parents=True)
    for page in ("index", f"{SUB}/setup"):
        (docs / f"{page}.md").write_text("# Page\n", encoding="utf-8")
    (site / "index.html").write_text('<h1 id="page">Page</h1>', encoding="utf-8")
    (site / SUB / "setup.html").write_text('<h2 id="step-1">Step 1</h2>', encoding="utf-8")
    toml = tmp_path / "book.toml"
    table = (
        "[output.html.redirect]\n"
        f'"/index.html#first" = "{SUB}/setup.html#step-1"\n'
        f'"/index.html#second" = "{SUB}/setup.html#step-2"\n'
    )
    toml.write_text(table, encoding="utf-8")
    done = subprocess.run(
        [sys.executable, "-", str(toml), str(docs), str(site), RELEASE, "/vep-rs/"],
        input=redirect_filter(), capture_output=True, text=True, check=False,
    )
    assert done.returncode == 0, done.stderr
    assert done.stdout == (
        f"redirect sent to the page, anchor not on the {RELEASE} page: /index.html#second -> {SUB}/setup.html#step-2\n"
        "redirects kept: 2, dropped: 0, sent to the page: 1\n"
    )
    assert toml.read_text(encoding="utf-8") == table.replace("setup.html#step-2", "setup.html")
