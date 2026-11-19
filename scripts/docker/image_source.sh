#!/usr/bin/env bash
# Prints what one run of .github/workflows/docker.yml packages and whether it
# pushes, as the `key=value` lines of a step's outputs:
#   source   `release`: a published release's attested archives
#            (stage_release.sh); `tree`: the binaries compiled from the
#            checked-out tree (stage_tree.sh)
#   tag      the release's tag; empty for `tree`
#   version  the version the image is labelled with and `vep --version` must
#            print: the tag's, or the workspace manifest's
#   push     `true` or `false`
# A `release` event packages its own tag and pushes. `workflow_dispatch`
# packages the tag it is given and pushes only when asked, from main; a tag
# older than FIRST_FITTING is refused by name. A `pull_request` packages its
# own tree and never pushes, so every pushed image holds attested release
# binaries. Any other event is refused, and a refused run prints no output
# line.
#
# Usage: EVENT=<event> [GITHUB_REF=<ref>] [RELEASE_TAG=<tag>] [INPUT_TAG=<tag>]
#            [INPUT_PUSH=true] scripts/docker/image_source.sh
# MANIFEST overrides the workspace manifest the tree's version is read from.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MANIFEST="${MANIFEST:-$ROOT/Cargo.toml}"
event="${EVENT:?EVENT names the event that started the run}"
ref="${GITHUB_REF:-}"

fail() {
    echo "ERROR: [image_source] $*" >&2
    exit 1
}

# The first release the tree's image fits. Its `vep` writes Parquet output
# itself; an older `vep` runs the `duckdb` CLI for it, docker/Dockerfile
# installs none, and smoke_test.sh checks that output.
FIRST_FITTING=v0.6.0

case "$event" in
    release)
        source=release
        tag="${RELEASE_TAG:-}"
        push=true
        ;;
    workflow_dispatch)
        source=release
        tag="${INPUT_TAG:-}"
        push=false
        if [[ "${INPUT_PUSH:-}" == true ]]; then
            [[ "$ref" == refs/heads/main ]] || fail "push=true is allowed only from main, not $ref"
            push=true
        fi
        ;;
    pull_request)
        source=tree
        tag=""
        push=false
        ;;
    *) fail "event '$event' is not one the workflow packages an image for" ;;
esac

if [[ "$source" == release ]]; then
    [[ "$tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || fail "tag '$tag' is not vX.Y.Z"
    # A release event packages its own tag from that tag's tree, whose image
    # fits it; only a dispatch pairs this tree's image with another release.
    if [[ "$event" == workflow_dispatch && "$tag" != "$FIRST_FITTING" &&
        "$(printf '%s\n' "$tag" "$FIRST_FITTING" | sort -V | head -n 1)" == "$tag" ]]; then
        fail "tag $tag is older than $FIRST_FITTING, the first release this image fits: its vep needs the duckdb CLI for Parquet output, and the image carries none"
    fi
    version="${tag#v}"
else
    # The first `version` line of the manifest, the one release.yml names a
    # dry run's archives by.
    version="$(awk -F'"' '/^version = "[^"]+"$/ { print $2; exit }' "$MANIFEST")"
    [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || fail "$MANIFEST states no X.Y.Z version: '$version'"
fi

printf '%s\n' "source=$source" "tag=$tag" "version=$version" "push=$push"
echo "packaging ${tag:-the tree at $version} (source=$source, push=$push)" >&2
