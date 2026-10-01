#!/usr/bin/env bash
# Prints the image references to push for one release, one per line:
# <image>:X.Y.Z always; <image>:X.Y when the tag is the newest published
# release in its X.Y line; <image>:latest when it is the newest of all.
# Published release tags are read from stdin, one per line, so a rebuild of an
# old release never moves X.Y or latest backwards.
#
# Usage: gh release list --exclude-drafts --exclude-pre-releases --json tagName \
#            --jq '.[].tagName' | scripts/docker/image_tags.sh <image> <tag>
set -euo pipefail

usage="usage: image_tags.sh <image> <tag> < published-tags"
image="${1:?$usage}"
tag="${2:?$usage}"

fail() {
    echo "ERROR: [image_tags] $*" >&2
    exit 1
}

semver='^v([0-9]+)\.([0-9]+)\.([0-9]+)$'
[[ "$tag" =~ $semver ]] || fail "tag '$tag' is not vX.Y.Z"
minor="${BASH_REMATCH[1]}.${BASH_REMATCH[2]}"

published="$(grep -E "$semver" || true)"
grep -qxF "$tag" <<<"$published" || fail "$tag is not a published release"
newest="$(sort -V <<<"$published" | tail -1)"
newest_in_minor="$(grep -F "v$minor." <<<"$published" | sort -V | tail -1)"

echo "$image:${tag#v}"
[[ "$tag" == "$newest_in_minor" ]] && echo "$image:$minor"
[[ "$tag" == "$newest" ]] && echo "$image:latest"
true
