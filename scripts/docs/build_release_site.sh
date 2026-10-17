#!/usr/bin/env bash
# Builds one channel of the documentation site and checks it: the site of a
# release, or the development preview of the checkout.
#
# Release form. The pages (README.md, CHANGELOG.md, CONTRIBUTING.md, docs/)
# come from the release tag; the site tooling (book.toml, theme/, scripts/docs/)
# comes from the current checkout, so a tooling fix applies to the release
# already published. A tag cut before a page of the site's frame existed
# (docs/SUMMARY.md, index.md, changelog.md, contributing.md, 404.md, docs/assets/)
# takes it from the checkout. A SUMMARY.md entry naming a page the tag does not
# have is dropped. So is a redirect whose destination page the tag does not
# have, or whose source anchor is still a heading of the tag's own page. Each
# drop prints a line; a redirect entry the filter cannot read (anything but
# "source" = "destination" on one line), or whose destination leaves the site,
# fails the build. A tag whose
# docs/index.md is absent, or is the stub that points at the README, keeps the
# README as its home page (VEP_RS_DOCS_README_HOME=1).
#
# Preview form (--main). The pages come from HEAD; the site is titled
# "vep-rs main", addressed under /vep-rs/main/ and marked with its channel, which
# the theme shows as the development-version banner and the header switcher.
# VEP_RS_DOCS_RELEASE names the latest release: the environment's
# value when set, else the newest tag reachable from origin/main, else the
# highest tag.
#
# Both forms end with scripts/docs/site_postbuild.py, which writes the canonical
# links, versions.json and, for a release, sitemap.xml, then checks the redirect
# destinations and the benchmark chart's release. The site's GitHub links name
# the ref built, its title carries it, and the link check holds repository paths
# to that tree. The site is written to <out-dir>/target/book.
#
# Usage: scripts/docs/build_release_site.sh <tag> <out-dir>
#        scripts/docs/build_release_site.sh --main <out-dir>
# Needs git, python3 and mdbook on PATH, and the tag fetched.
set -euo pipefail

usage="usage: build_release_site.sh <tag> <out-dir> | build_release_site.sh --main <out-dir>"
pages_origin="https://natera-open-source.github.io"

fail() {
    echo "ERROR: [build_release_site] $*" >&2
    exit 1
}

latest_release() {
    git describe --tags --abbrev=0 --match 'v[0-9]*' origin/main 2>/dev/null && return
    git tag --list 'v[0-9]*' --sort=-v:refname | head -n 1
}

if [ "${1:-}" = --main ]; then
    out="${2:?$usage}"
    ref=HEAD
    docs_ref=main
    title="vep-rs main"
    site_url=/vep-rs/main/
    release="${VEP_RS_DOCS_RELEASE:-$(latest_release)}"
    [ -n "$release" ] || fail "no release tag found; set VEP_RS_DOCS_RELEASE"
    channel=main
else
    ref="${1:?$usage}"
    out="${2:?$usage}"
    git rev-parse --verify --quiet "${ref}^{commit}" >/dev/null || fail "no such tag: $ref"
    docs_ref="$ref"
    title="vep-rs ${ref}"
    site_url=/vep-rs/
    release="$ref"
    channel=
fi
[ ! -e "$out" ] || [ -z "$(ls -A "$out")" ] || fail "$out exists and is not empty"
checkout="$(git rev-parse --show-toplevel)"

in_ref() {
    [ -n "$(git ls-tree --name-only "$ref" -- "$1")" ]
}

mkdir -p "$out"
pages=(README.md CHANGELOG.md docs)
if in_ref CONTRIBUTING.md; then
    pages+=(CONTRIBUTING.md)
fi
git archive "$ref" "${pages[@]}" | tar -x -C "$out"

cp "$checkout/book.toml" "$out/"
[ ! -d "$checkout/theme" ] || cp -R "$checkout/theme" "$out/"
mkdir -p "$out/scripts"
cp -R "$checkout/scripts/docs" "$out/scripts/"

index="${out}/docs/index.md"
if [ ! -e "$index" ] || {
    [ "$(head -n 1 "$index")" = "# vep-rs documentation" ] && [ "$(awk 'END { print NR }' "$index")" -eq 3 ]
}; then
    export VEP_RS_DOCS_README_HOME=1
    echo "$docs_ref has no landing page; the README is the home page"
fi

summary_from_checkout=
for page in SUMMARY.md index.md changelog.md contributing.md 404.md; do
    [ ! -e "$out/docs/$page" ] || continue
    if [ "$page" = contributing.md ] && [ ! -e "$out/CONTRIBUTING.md" ]; then
        echo "$docs_ref has no CONTRIBUTING.md; docs/$page is left out so its SUMMARY.md entry drops"
        continue
    fi
    if [ ! -e "$checkout/docs/$page" ]; then
        echo "neither $docs_ref nor the checkout has docs/$page"
        continue
    fi
    cp "$checkout/docs/$page" "$out/docs/$page"
    echo "$docs_ref has no docs/$page; using the checkout's"
    [ "$page" != SUMMARY.md ] || summary_from_checkout=1
done
if [ -n "$summary_from_checkout" ]; then
    python3 - "$out/docs" <<'EOF'
import os, re, sys

docs = sys.argv[1]
path = os.path.join(docs, "SUMMARY.md")
kept = []
for line in open(path, encoding="utf-8"):
    m = re.search(r"\]\(([^)#]+)\)", line)
    if m and not os.path.exists(os.path.join(docs, m[1])):
        print(f"dropped from SUMMARY.md, not in the tag: {m[1]}")
        continue
    kept.append(line)
open(path, "w", encoding="utf-8").writelines(kept)
EOF
fi

if [ ! -e "$out/docs/assets" ] && [ -d "$checkout/docs/assets" ]; then
    cp -R "$checkout/docs/assets" "$out/docs/"
    echo "$docs_ref has no docs/assets; using the checkout's"
fi

build_site() {
    (
        cd "$out"
        MDBOOK_BOOK__TITLE="$title" MDBOOK_OUTPUT__HTML__SITE_URL="$site_url" \
            VEP_RS_DOCS_REF="$docs_ref" VEP_RS_DOCS_RELEASE="$release" mdbook build
    )
}

# mdBook writes a redirect without checking that its destination exists, and it
# redirects a fragment even when the source page still has that heading, so an
# unfiltered table would send a reader to a 404 page or away from the heading
# they asked for. The ids a page has come from a first build: mdBook slugifies
# the headings and the preprocessor replaces whole pages, so the Markdown alone
# does not show them.
if grep -q '^\[output\.html\.redirect\]' "$out/book.toml"; then
    build_site
    python3 - "$out/book.toml" "$out/docs" "$out/target/book" "$docs_ref" "$site_url" <<'EOF'
import html.parser, os, posixpath, re, sys

toml, docs, site, ref, site_url = sys.argv[1:6]
ENTRY = re.compile(r'^\s*"([^"]+)"\s*=\s*"([^"]*)"\s*(#.*)?$')
SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")


class Ids(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        self.ids.update(value for name, value in attrs if name == "id" and value)


page_ids = {}


def ids_of(page):
    if page not in page_ids:
        parser = Ids()
        path = os.path.join(site, page)
        if os.path.isfile(path):
            parser.feed(open(path, encoding="utf-8").read())
        page_ids[page] = parser.ids
    return page_ids[page]


def destination_exists(key, dest):
    page = dest.split("#")[0]
    if not page or SCHEME.match(page):
        return True
    if page.startswith("/"):
        # A leading slash is a browser path, as mdBook writes it, so it must
        # start with the site's own URL.
        if page == site_url.rstrip("/"):
            page = site_url
        if not page.startswith(site_url):
            sys.exit(f"ERROR: [build_release_site] redirect destination leaves the site: {key} -> {dest}")
        # An empty segment after the site URL (//) leaves an absolute path,
        # which os.path.join takes in place of the docs root.
        resolved = page[len(site_url):].lstrip("/")
    else:
        resolved = posixpath.normpath(posixpath.join(posixpath.dirname(key.lstrip("/")), page))
    if page.endswith("/"):
        resolved = posixpath.join(resolved, "index.md")
    return os.path.exists(os.path.join(docs, re.sub(r"\.html$", ".md", resolved)))


kept, in_table, seen, dropped = [], False, 0, 0
for line in open(toml, encoding="utf-8"):
    if line.lstrip().startswith("["):
        in_table = line.strip() == "[output.html.redirect]"
    elif in_table and line.strip() and not line.lstrip().startswith("#"):
        m = ENTRY.match(line)
        if not m:
            sys.exit(f"ERROR: [build_release_site] unrecognised redirect entry: {line.rstrip()}")
        seen += 1
        key, dest = m[1], m[2]
        source, _, fragment = key.lstrip("/").partition("#")
        if not destination_exists(key, dest):
            print(f"dropped redirect, destination not in {ref}: {key} -> {dest}")
            dropped += 1
            continue
        if fragment and fragment in ids_of(source):
            print(f"dropped redirect, anchor still on the {ref} page: {key} -> {dest}")
            dropped += 1
            continue
    kept.append(line)
open(toml, "w", encoding="utf-8").writelines(kept)
print(f"redirects kept: {seen - dropped}, dropped: {dropped}")
EOF
fi

build_site
python3 "$checkout/scripts/docs/check_site_links.py" --site "$out/target/book" --ref "$docs_ref" --tree "$ref"
postbuild=("$out/target/book" --base "${pages_origin}${site_url}" --release "$release" --source "$out")
[ -z "$channel" ] || postbuild+=(--channel "$channel")
python3 "$checkout/scripts/docs/site_postbuild.py" "${postbuild[@]}"
echo "site for $docs_ref: $out/target/book"
