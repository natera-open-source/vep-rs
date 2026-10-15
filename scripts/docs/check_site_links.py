#!/usr/bin/env python3
"""Check the links of the built documentation site.

Every relative link in the site must name a file the build wrote and, when it
carries a fragment, an `id` on that page. Every link into this repository on
GitHub (`<repo>/blob/<ref>/<path>` or `<repo>/tree/<ref>/<path>`) must name a
path that exists in the tree the site was built from. External links are not
fetched.

Usage::

    python3 scripts/docs/check_site_links.py [--site target/book] [--ref main] [--tree TREE-ISH]

Exits 1 and lists each broken link if any. Stdlib only.
"""

import argparse
import html.parser
import os
import subprocess
import sys
import urllib.parse

REPO_URL = "https://github.com/natera-open-source/vep-rs"
# Pages mdBook generates whose links all point at the chapters themselves.
SKIP_PAGES = {"print.html", "toc.html", "404.html"}


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


def main():
    ap = argparse.ArgumentParser(description="Check the links of the built documentation site.")
    ap.add_argument("--site", default="target/book")
    ap.add_argument("--repo-url", default=REPO_URL)
    ap.add_argument("--ref", default="main", help="the ref the site's GitHub links name")
    ap.add_argument("--tree", help="git tree-ish the site was built from; default the working tree")
    args = ap.parse_args()

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
    else:
        in_repo = os.path.exists

    prefixes = [args.repo_url.rstrip("/") + f"/{kind}/{args.ref}/" for kind in ("blob", "tree")]
    broken, checked = [], 0
    for page_path, page in sorted(pages.items()):
        if page_path in SKIP_PAGES:
            continue
        for href in page.links:
            url = urllib.parse.urlsplit(href)
            repo_prefix = next((p for p in prefixes if href.startswith(p)), None)
            if repo_prefix:
                checked += 1
                target = urllib.parse.unquote(href[len(repo_prefix) :].split("#")[0])
                if not in_repo(target.rstrip("/")):
                    broken.append((page_path, href, "no such path in the repository"))
                continue
            if url.scheme or url.netloc or href.startswith("/"):
                continue
            checked += 1
            if url.path:
                target = os.path.normpath(
                    os.path.join(os.path.dirname(page_path), urllib.parse.unquote(url.path))
                ).replace(os.sep, "/")
            else:
                target = page_path
            if not os.path.exists(os.path.join(args.site, target)):
                broken.append((page_path, href, "no such file in the site"))
            elif url.fragment and target in pages and url.fragment not in pages[target].ids:
                broken.append((page_path, href, f"no id {url.fragment!r} on {target}"))

    for page_path, href, why in broken:
        print(f"{page_path}: {href}: {why}")
    print(f"checked {checked} links on {len(pages) - len(SKIP_PAGES & pages.keys())} pages, {len(broken)} broken")
    sys.exit(1 if broken else 0)


if __name__ == "__main__":
    main()
