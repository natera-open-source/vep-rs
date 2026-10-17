"""Tests for the site tooling in this directory, on fixtures the tests write: the mdBook preprocessor driven
through its stdin protocol, the page metrics and their gate, the post-build finish with its two checks, the
benchmark chart renderer, and the screenshot script's pure helpers. No test runs mdbook or a browser."""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import mdbook_repo_links as links
import render_benchmark_chart as chart
import site_metrics as metrics
import site_postbuild as postbuild
import site_screenshots as shots

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


VERSIONS = {
    "release": RELEASE,
    "channels": [{"name": RELEASE, "path": "/vep-rs/"}, {"name": "main", "path": "/vep-rs/main/"}],
}


def page_html(title, body):
    return (
        "<!DOCTYPE html>\n<html>\n<head>\n"
        f"<title>{title} - vep-rs {RELEASE}</title>\n</head>\n<body>\n"
        f'<h1 class="menu-title">vep-rs {RELEASE}</h1>\n{body}\n</body>\n</html>\n'
    )


def redirect_html(destination):
    return (
        "<!DOCTYPE html>\n<html>\n<head>\n"
        f'<meta http-equiv="refresh" content="0; URL={destination}">\n</head>\n<body></body>\n</html>\n'
    )


def built_site(tmp_path, chart_release=RELEASE, redirect="guide/setup.html#step-2"):
    """A built site of two reader pages, two redirect stubs, the two pages mdBook generates for printing and
    not-found, and the benchmark chart; beside it the source tree whose README names the release."""
    site = tmp_path / "book"
    (site / "guide").mkdir(parents=True)
    (site / "assets").mkdir()
    edit = f'<a href="{EDIT_URL}docs/index.md">Edit</a>\n<h1 id="overview">Overview</h1>'
    (site / "index.html").write_text(page_html("vep-rs", edit), encoding="utf-8")
    (site / "guide" / "setup.html").write_text(
        page_html("Setup", '<h1 id="setup">Setup</h1>\n<h2 id="step-2">Step 2</h2>'), encoding="utf-8"
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
    assert f'<p class="menu-title">vep-rs {RELEASE}</p>' in index
    assert "<h1 class=" not in index
    assert f"<title>vep-rs {RELEASE}</title>" in index
    assert f"<title>Setup - vep-rs {RELEASE}</title>" in (site / "guide" / "setup.html").read_text(encoding="utf-8")


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
    assert "benchmark chart: not shown, the README is the home page; release check skipped\n" in out
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
