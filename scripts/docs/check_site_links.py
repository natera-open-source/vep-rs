#!/usr/bin/env python3
"""Check the links of the built documentation site.

Every relative link in the site must name a file the build wrote and, when it
carries a fragment, an `id` on that page. Every link into this repository on
GitHub (`<repo>/blob/<ref>/<path>` or `<repo>/tree/<ref>/<path>`) must name a
path that exists in the tree the site was built from. Every destination of the
redirect table of the book the site was built from (`[output.html.redirect]`
in `book.toml`) is held to the same rules, read as a link on the page it
redirects from: mdBook writes a redirect without checking where it lands. A
destination that ends in a slash names the `index.html` of that directory.
External links are not fetched.

Every Markdown file under the book's source directory (`src` of `[book]` in
`book.toml`), its `SUMMARY.md` excepted, must have its page in the site: mdBook
builds a page for each chapter `SUMMARY.md` lists and passes over every other
file without a message. A `README.md` is held to the `index.html` of its
directory. The files are those of the tree `--tree` names, read from that
tree's root, or else those beside the book configuration.

Usage::

    python3 scripts/docs/check_site_links.py [--site target/book] [--book book.toml] [--ref main] [--tree TREE-ISH]

Exits 1 and lists each broken link, each broken redirect and each page the
site lacks, if any. Stdlib only.
"""

import argparse
import html.parser
import os
import posixpath
import re
import subprocess
import sys
import urllib.parse

REPO_URL = "https://github.com/natera-open-source/vep-rs"
# Pages mdBook generates whose links all point at the chapters themselves.
SKIP_PAGES = {"print.html", "toc.html", "404.html"}
REDIRECT_TABLE = "[output.html.redirect]"
# The one entry shape scripts/docs/build_release_site.sh reads when it filters the table for a site.
REDIRECT_ENTRY = re.compile(r'^\s*"([^"]+)"\s*=\s*"([^"]*)"\s*(#.*)?$')
BOOK_TABLE = "[book]"
SOURCE_KEY = re.compile(r'^\s*src\s*=\s*"([^"]*)"\s*(#.*)?$')
# mdBook reads a book without the key from this directory.
DEFAULT_SOURCE = "src"


class Page(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids, self.links = set(), []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("href"):
            self.links.append(attrs["href"])


def parse(path):
    page = Page()
    with open(path, encoding="utf-8") as f:
        page.feed(f.read())
    return page


def read_redirects(book):
    """The (source, destination) pairs of the redirect table in the book configuration `book`, in file order.
    Exits with one line when the file is absent or an entry of the table is not `"source" = "destination"`."""
    if not os.path.isfile(book):
        sys.exit(
            f"ERROR: [check_site_links] no book configuration at {book}: "
            "name the one the site was built from with --book"
        )
    entries, in_table = [], False
    with open(book, encoding="utf-8") as f:
        for line in f:
            if line.lstrip().startswith("["):
                in_table = line.split("#", 1)[0].strip() == REDIRECT_TABLE
            elif in_table and line.strip() and not line.lstrip().startswith("#"):
                entry = REDIRECT_ENTRY.match(line)
                if not entry:
                    sys.exit(f"ERROR: [check_site_links] unrecognised redirect entry in {book}: {line.strip()}")
                entries.append((entry[1], entry[2]))
    return entries


def read_source(book):
    """The source directory the book configuration `book` names (`src` of its `[book]` table), as written
    there, or mdBook's default when the table has no such key."""
    in_table = False
    with open(book, encoding="utf-8") as f:
        for line in f:
            if line.lstrip().startswith("["):
                in_table = line.split("#", 1)[0].strip() == BOOK_TABLE
            elif in_table:
                key = SOURCE_KEY.match(line)
                if key:
                    return key[1]
    return DEFAULT_SOURCE


def source_pages(directory, listing):
    """The Markdown files under the source directory `directory` as sorted paths from it, the chapter list
    `SUMMARY.md` at its top left out: the paths of `listing` (a tree's paths from its root, its directories
    among them) when that is a list, else the files on disk."""
    if listing is None:
        found = [
            os.path.relpath(os.path.join(dirpath, name), directory).replace(os.sep, "/")
            for dirpath, _, files in os.walk(directory)
            for name in files
        ]
    else:
        top = directory + "/"
        folders = {posixpath.dirname(path) for path in listing}
        found = [path[len(top) :] for path in listing if path.startswith(top) and path not in folders]
    return sorted(page for page in found if page.endswith(".md") and page != "SUMMARY.md")


def built_page(source_page):
    """The site file mdBook writes for the source page `source_page`: the same path ending in `.html`, and the
    `index.html` of its directory for a `README.md` in any letter case."""
    directory, name = posixpath.split(source_page)
    stem = name[: -len(".md")]
    return posixpath.join(directory, ("index" if stem.lower() == "readme" else stem) + ".html")


def link_fault(href, page_path, site, pages, prefixes, in_repo):
    """Why `href`, read as a link on the page `page_path`, lands nowhere: an empty string when it lands, and
    None when it is not followed (a URL outside this repository, or a path from the server root)."""
    repo_prefix = next((p for p in prefixes if href.startswith(p)), None)
    if repo_prefix:
        target = urllib.parse.unquote(href[len(repo_prefix) :].split("#")[0])
        return "" if in_repo(target.rstrip("/")) else "no such path in the repository"
    url = urllib.parse.urlsplit(href)
    if url.scheme or url.netloc or href.startswith("/"):
        return None
    if url.path:
        target = os.path.normpath(
            os.path.join(os.path.dirname(page_path), urllib.parse.unquote(url.path))
        ).replace(os.sep, "/")
    else:
        target = page_path
    if not os.path.exists(os.path.join(site, target)):
        return "no such file in the site"
    if url.fragment and target in pages and url.fragment not in pages[target].ids:
        return f"no id {url.fragment!r} on {target}"
    return ""


def served_page(destination):
    """The redirect destination `destination` as the page a server answers it with: a path that ends in a
    slash names the `index.html` of that directory."""
    url = urllib.parse.urlsplit(destination)
    if url.scheme or url.netloc or not url.path.endswith("/"):
        return destination
    return urllib.parse.urlunsplit(url._replace(path=url.path + "index.html"))


def main():
    ap = argparse.ArgumentParser(description="Check the links of the built documentation site.")
    ap.add_argument("--site", default="target/book")
    ap.add_argument("--book", default="book.toml", help="the book configuration the site was built from")
    ap.add_argument("--repo-url", default=REPO_URL)
    ap.add_argument("--ref", default="main", help="the ref the site's GitHub links name")
    ap.add_argument("--tree", help="git tree-ish the site was built from; default the working tree")
    args = ap.parse_args()

    redirects = read_redirects(args.book)
    book_source = read_source(args.book)
    pages = {}
    for dirpath, _, files in os.walk(args.site):
        for name in files:
            if name.endswith(".html"):
                rel = os.path.relpath(os.path.join(dirpath, name), args.site)
                pages[rel.replace(os.sep, "/")] = parse(os.path.join(dirpath, name))

    if args.tree:
        listing = subprocess.run(
            ["git", "ls-tree", "-r", "-t", "--name-only", args.tree],
            check=True, capture_output=True, text=True,
        ).stdout.splitlines()
        in_repo = set(listing).__contains__
        source_dir = posixpath.normpath(book_source)
        held = f"{source_dir} in {args.tree}"
    else:
        listing = None
        in_repo = os.path.exists
        source_dir = held = os.path.normpath(os.path.join(os.path.dirname(args.book), book_source))
    sources = source_pages(source_dir, listing)
    if not sources:
        sys.exit(
            f"ERROR: [check_site_links] no Markdown page under {held}, the source directory of the book {args.book}"
        )

    prefixes = [args.repo_url.rstrip("/") + f"/{kind}/{args.ref}/" for kind in ("blob", "tree")]
    broken, checked = [], 0
    for page_path, page in sorted(pages.items()):
        if page_path in SKIP_PAGES:
            continue
        for href in page.links:
            fault = link_fault(href, page_path, args.site, pages, prefixes, in_repo)
            if fault is None:
                continue
            checked += 1
            if fault:
                broken.append((page_path, href, fault))

    lost, followed = [], 0
    for source, destination in redirects:
        # mdBook resolves a destination against the address of the page the source names.
        fault = link_fault(
            served_page(destination), source.lstrip("/").split("#")[0], args.site, pages, prefixes, in_repo
        )
        if fault is None:
            continue
        followed += 1
        if fault:
            lost.append((source, destination, fault))

    unbuilt = [page for page in sources if not os.path.isfile(os.path.join(args.site, built_page(page)))]

    for page_path, href, why in broken:
        print(f"{page_path}: {href}: {why}")
    for source, destination, why in lost:
        print(f"{args.book}: {source} -> {destination}: {why}")
    shown = source_dir.replace(os.sep, "/")
    for page in unbuilt:
        print(
            f"{shown}/{page}: no page {built_page(page)} in the site: "
            f"mdBook builds the chapters {shown}/SUMMARY.md lists"
        )
    print(f"checked {checked} links on {len(pages) - len(SKIP_PAGES & pages.keys())} pages, {len(broken)} broken")
    print(
        f"redirects of {args.book}: {followed} checked, {len(redirects) - followed} not followed "
        f"(a URL or a path from the server root), {len(lost)} broken"
    )
    print(f"pages under {held}: {len(sources)} checked, {len(unbuilt)} not in the site")
    sys.exit(1 if broken or lost or unbuilt else 0)


if __name__ == "__main__":
    main()
