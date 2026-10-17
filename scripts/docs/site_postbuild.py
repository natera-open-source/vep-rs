#!/usr/bin/env python3
"""Finish one channel of the built documentation site and check what the build cannot.

Runs after `mdbook build` on the release channel at the site root or on the
development preview of main under `main/`. On every page it rewrites the menu
bar's `<h1 class="menu-title">` as a `<p>`, so the page's own heading is its
one `h1`, and collapses a `<title>` whose chapter name repeats the book title.
On every page except `print.html`, `toc.html`, `404.html` and the stub mdBook
writes for a redirected path, it inserts `<link rel="canonical">` and `og:url`
for the page's address under `--base` and, on the preview,
`<meta name="vep-rs-channel" content="main">` and a `noindex` robots meta, so a
search engine lists the release and not the preview. It writes `versions.json`,
which the header's channel switcher reads, and for a release `sitemap.xml`.

When `VEP_RS_DOCS_README_HOME` is `1` the home page renders the README, so its
edit link is pointed at `README.md` and the benchmark chart, which that page
does not show, is not held to the README.

Then two checks. Every redirect destination the built pages name must exist:
the page as a file and, when the destination carries one, the id on that page.
The benchmark chart's `data-release` must equal the version in the README's
`## Current release: vX.Y.Z` heading of the tree the site was built from
(`--source`). Every step prints one line; any failure exits 1 after all of
them have run.

Usage::

    python3 scripts/docs/site_postbuild.py <site> --base <url> --release <tag>
        [--source <dir>] [--channel main]

Stdlib only.
"""

import argparse
import html
import html.parser
import json
import os
import posixpath
import re
import sys
import urllib.parse

# Pages mdBook generates whose address is not a reader's destination.
SKIP_PAGES = {"print.html", "toc.html", "404.html"}
PREVIEW_CHANNEL_NAME = "main"
MENU_TITLE = re.compile(r'<h1 class="menu-title">(.*?)</h1>', re.DOTALL)
TITLE = re.compile(r"<title>(.*?) - (.*?)</title>")
EDIT_HOME = re.compile(r'(href="[^"]*/edit/[^"]*/)[^"/]+/index\.md"')
REFRESH = re.compile(r'http-equiv="refresh"\s+content="0;\s*URL=([^"]+)"', re.IGNORECASE)
FRAGMENT_MAP = re.compile(r"fragmentMap\s*=\s*(\{.*?\})\s*;", re.DOTALL)
REDIRECT_URL = re.compile(r'let redirectUrl = "([^"]+)"')
SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
SVG_RELEASE = re.compile(r'<svg[^>]*\sdata-release="([^"]+)"')
README_RELEASE = re.compile(r"^## Current release: (\S+)\s*$", re.MULTILINE)


class Ids(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name == "id" and value:
                self.ids.add(value)


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def write(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def site_pages(site):
    pages = []
    for dirpath, _, files in os.walk(site):
        for name in files:
            if name.endswith(".html"):
                rel = os.path.relpath(os.path.join(dirpath, name), site)
                pages.append(rel.replace(os.sep, "/"))
    return sorted(pages)


def page_url(base, page):
    if posixpath.basename(page) == "index.html":
        directory = posixpath.dirname(page)
        return base + (directory + "/" if directory else "")
    return base + page


def collapse_title(m):
    """`<chapter> - <book>` becomes `<book>` when the chapter is the book's name or its first word."""
    chapter, book = m[1], m[2]
    return f"<title>{book}</title>" if chapter in (book, book.split(" ")[0]) else m[0]


def finish_pages(site, pages, reader_pages, base, channel, readme_home):
    problems, inserted, present = [], 0, 0
    for page in pages:
        path = os.path.join(site, page)
        text = MENU_TITLE.sub(r'<p class="menu-title">\1</p>', read(path))
        text = TITLE.sub(collapse_title, text, 1)
        if readme_home and page == "index.html":
            text = EDIT_HOME.sub(r'\1README.md"', text)
        if page in reader_pages:
            if 'rel="canonical"' in text:
                present += 1
            elif "</head>" not in text:
                problems.append(f"{page}: no </head> to insert the canonical link before")
            else:
                url = html.escape(page_url(base, page))
                tags = [f'<link rel="canonical" href="{url}">', f'<meta property="og:url" content="{url}">']
                if channel:
                    tags.append(f'<meta name="vep-rs-channel" content="{channel}">')
                    tags.append('<meta name="robots" content="noindex">')
                text = text.replace("</head>", "\n        ".join(tags) + "\n    </head>", 1)
                inserted += 1
        write(path, text)
    what = "canonical, og:url" + (", channel, noindex" if channel else "")
    print(
        f"{len(pages)} pages finished: menu title demoted, head tags ({what}) inserted on {inserted},"
        f" {present} already carried a canonical link"
    )
    return problems


def write_versions(site, base, release, channel):
    root = urllib.parse.urlsplit(base).path
    if channel:
        suffix = channel + "/"
        if not root.endswith(suffix):
            return [f"--base path {root} does not end in {suffix}, so the release path cannot be derived"]
        root = root[: -len(suffix)]
    preview_path = root + PREVIEW_CHANNEL_NAME + "/"
    versions = {
        "release": release,
        "channels": [
            {"name": release, "path": root},
            {"name": PREVIEW_CHANNEL_NAME, "path": preview_path},
        ],
    }
    write(os.path.join(site, "versions.json"), json.dumps(versions, indent=2) + "\n")
    print(f"versions.json: release {release}, channels at {root} and {preview_path}")
    return []


def write_sitemap(site, pages, base):
    lines = ['<?xml version="1.0" encoding="UTF-8"?>']
    lines.append('<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">')
    for page in pages:
        lines.append(f"  <url><loc>{html.escape(page_url(base, page))}</loc></url>")
    lines.append("</urlset>")
    write(os.path.join(site, "sitemap.xml"), "\n".join(lines) + "\n")
    print(f"sitemap.xml: {len(pages)} pages")


def redirect_destinations(text):
    """Every destination a page's redirect script or refresh tag names, in page order."""
    destinations = []
    for m in FRAGMENT_MAP.finditer(text):
        destinations.extend(json.loads(m[1]).values())
    destinations.extend(m[1] for m in REDIRECT_URL.finditer(text))
    destinations.extend(m[1] for m in REFRESH.finditer(text))
    return list(dict.fromkeys(destinations))


def check_redirects(site, pages, base):
    base_path = urllib.parse.urlsplit(base).path
    ids = {}
    broken, checked, external = [], 0, 0
    for page in pages:
        for dest in redirect_destinations(read(os.path.join(site, page))):
            if SCHEME.match(dest):
                external += 1
                continue
            checked += 1
            path, _, fragment = dest.partition("#")
            if path.endswith("/"):
                path += "index.html"
            if not path:
                target = page
            elif path.startswith("/"):
                if not path.startswith(base_path) and path != base_path.rstrip("/"):
                    broken.append(f"{page}: redirect to {dest} leaves the site")
                    continue
                target = path[len(base_path):]
            else:
                target = posixpath.normpath(posixpath.join(posixpath.dirname(page), path))
            if not target:
                target = "index.html"
            if not os.path.isfile(os.path.join(site, target)):
                broken.append(f"{page}: redirect to {dest}: no such file in the site")
                continue
            if fragment and target.endswith(".html"):
                if target not in ids:
                    parser = Ids()
                    parser.feed(read(os.path.join(site, target)))
                    ids[target] = parser.ids
                if fragment not in ids[target]:
                    broken.append(f"{page}: redirect to {dest}: no id {fragment!r} on {target}")
    print(f"redirect destinations: {checked} checked, {external} external not fetched, {len(broken)} broken")
    return broken


def check_chart_release(site, source, readme_home):
    svg = os.path.join(site, "assets", "benchmark.svg")
    if not os.path.isfile(svg):
        print("benchmark chart: not in the site, release check skipped")
        return []
    if readme_home:
        print("benchmark chart: not shown, the README is the home page; release check skipped")
        return []
    m = SVG_RELEASE.search(read(svg))
    if not m:
        return ["assets/benchmark.svg carries no data-release attribute on its root element"]
    chart_release = m[1]
    if not source:
        return [f"assets/benchmark.svg says {chart_release}; --source is needed to compare it with the README"]
    readme = os.path.join(source, "README.md")
    if not os.path.isfile(readme):
        return [f"{readme} does not exist, so the chart's release cannot be checked"]
    heading = README_RELEASE.search(read(readme))
    if not heading:
        return [f"{readme} has no '## Current release: vX.Y.Z' heading, so the chart's release cannot be checked"]
    if heading[1] != chart_release:
        return [f"benchmark chart says {chart_release}, the README's Current release heading says {heading[1]}"]
    print(f"benchmark chart: {chart_release} matches the README")
    return []


def main():
    ap = argparse.ArgumentParser(description="Finish one channel of the built documentation site and check it.")
    ap.add_argument("site", help="the built site, mdBook's output directory")
    ap.add_argument("--base", required=True, help="the URL the site is served at, ending in a slash")
    ap.add_argument("--release", required=True, help="the release tag the site root serves (vX.Y.Z)")
    ap.add_argument("--source", help="the tree the site was built from; holds the README the chart is held to")
    ap.add_argument("--channel", help="the channel built under <base>; omitted for the release at the root")
    args = ap.parse_args()
    base = args.base if args.base.endswith("/") else args.base + "/"
    readme_home = os.environ.get("VEP_RS_DOCS_README_HOME") == "1"

    pages = site_pages(args.site)
    if not pages:
        sys.exit(f"ERROR: [site_postbuild] no .html page under {args.site}")
    reader_pages = {
        p for p in pages if p not in SKIP_PAGES and not REFRESH.search(read(os.path.join(args.site, p)))
    }

    failures = []
    failures += finish_pages(args.site, pages, reader_pages, base, args.channel, readme_home)
    failures += write_versions(args.site, base, args.release, args.channel)
    if not args.channel:
        write_sitemap(args.site, sorted(reader_pages), base)
    failures += check_redirects(args.site, pages, base)
    failures += check_chart_release(args.site, args.source, readme_home)

    for failure in failures:
        print(f"ERROR: [site_postbuild] {failure}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
