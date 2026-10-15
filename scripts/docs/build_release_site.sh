#!/usr/bin/env bash
# Builds the documentation site of one release and checks its links.
#
# The pages (README.md, CHANGELOG.md, docs/) come from the release tag; the
# site tooling (book.toml, scripts/docs/) comes from the current checkout, so a
# tooling fix applies to the release already published. A tag cut before the
# site existed has no docs/SUMMARY.md, docs/index.md or docs/changelog.md: those
# are taken from the checkout, and a SUMMARY.md entry naming a page the tag does
# not have is dropped, with a line saying so.
#
# The site's GitHub links name the tag, its title carries the tag, and the link
# check holds repository paths to the tag's tree. The site is written to
# <out-dir>/target/book.
#
# Usage: scripts/docs/build_release_site.sh <tag> <out-dir>
# Needs git, python3 and mdbook on PATH, and the tag fetched.
set -euo pipefail

usage="usage: build_release_site.sh <tag> <out-dir>"
tag="${1:?$usage}"
out="${2:?$usage}"

fail() {
    echo "ERROR: [build_release_site] $*" >&2
    exit 1
}

git rev-parse --verify --quiet "${tag}^{commit}" >/dev/null || fail "no such tag: $tag"
[ ! -e "$out" ] || [ -z "$(ls -A "$out")" ] || fail "$out exists and is not empty"
checkout="$(git rev-parse --show-toplevel)"

mkdir -p "$out"
git archive "$tag" README.md CHANGELOG.md docs | tar -x -C "$out"

cp "$checkout/book.toml" "$out/"
mkdir -p "$out/scripts"
cp -R "$checkout/scripts/docs" "$out/scripts/"

summary_from_checkout=
for page in SUMMARY.md index.md changelog.md; do
    if [ ! -e "$out/docs/$page" ]; then
        cp "$checkout/docs/$page" "$out/docs/$page"
        echo "$tag has no docs/$page; using the checkout's"
        [ "$page" != SUMMARY.md ] || summary_from_checkout=1
    fi
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

(
    cd "$out"
    MDBOOK_BOOK__TITLE="vep-rs ${tag}" VEP_RS_DOCS_REF="$tag" mdbook build
)
python3 "$checkout/scripts/docs/check_site_links.py" --site "$out/target/book" --ref "$tag" --tree "$tag"
echo "site for $tag: $out/target/book"
