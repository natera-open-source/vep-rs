#!/usr/bin/env bash
# Tests for image_source.sh: what each event of .github/workflows/docker.yml
# packages and whether it pushes, the runs it refuses, a dispatch of a tag
# older than the first release the image fits among them, and the workflow
# steps that take its answer.
#
# Pure bash. No network, no GitHub.
# Run: bash scripts/docker/test_image_source.sh
# shellcheck disable=SC2016  # the workflow's text is matched as written, unexpanded
set -uo pipefail

PASS=0
FAIL=0
ok() {
    PASS=$((PASS + 1))
    echo "  ok   $1"
}
bad() {
    FAIL=$((FAIL + 1))
    echo "  FAIL $1"
}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE="$HERE/image_source.sh"
ROOT=$(cd "$HERE/../.." && pwd)
WORKFLOW="$ROOT/.github/workflows/docker.yml"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# A workspace manifest whose first `version` line is the workspace's own.
printf '%s\n' '[workspace.package]' 'version = "9.9.9"' 'edition = "2021"' '' '[workspace.dependencies.dep]' 'version = "1.0.0"' >"$TMP/Cargo.toml"
printf '%s\n' '[workspace.package]' 'edition = "2021"' >"$TMP/no-version.toml"
printf '%s\n' '[workspace.package]' 'version = "9.9.9-rc1"' >"$TMP/pre-release.toml"

run() { # run [VAR=value ...]: image_source.sh in an environment of those variables alone
    env -i PATH="$PATH" MANIFEST="$TMP/Cargo.toml" "$@" bash "$SOURCE"
}

expect() { # expect <name> <expected stdout> [VAR=value ...]
    local name="$1" want="$2" got rc
    shift 2
    got=$(run "$@" 2>/dev/null)
    rc=$?
    if [[ "$rc" -eq 0 && "$got" == "$want" ]]; then ok "$name"; else bad "$name: rc=$rc, got [$got], want [$want]"; fi
}

expect_fail() { # expect_fail <name> <stderr substring> [VAR=value ...]: exits non-zero with no output line
    local name="$1" want="$2" out err rc
    shift 2
    out=$(run "$@" 2>"$TMP/stderr")
    rc=$?
    err=$(<"$TMP/stderr")
    if [[ "$rc" -ne 0 && "$err" == *"$want"* && -z "$out" ]]; then ok "$name"; else bad "$name (rc=$rc, stdout: $out, stderr: $err)"; fi
}

bash -n "$SOURCE" && ok "image_source.sh parses" || bad "image_source.sh has a syntax error"

RELEASE=$'source=release\ntag=v1.2.3\nversion=1.2.3\npush=true'
UNPUSHED=$'source=release\ntag=v1.2.3\nversion=1.2.3\npush=false'
TREE=$'source=tree\ntag=\nversion=9.9.9\npush=false'

expect "a release event packages its own tag and pushes" "$RELEASE" \
    EVENT=release RELEASE_TAG=v1.2.3 GITHUB_REF=refs/tags/v1.2.3
expect "a dispatch packages the tag it is given and does not push unasked" "$UNPUSHED" \
    EVENT=workflow_dispatch INPUT_TAG=v1.2.3 INPUT_PUSH=false GITHUB_REF=refs/heads/topic
expect "a dispatch from main pushes when asked" "$RELEASE" \
    EVENT=workflow_dispatch INPUT_TAG=v1.2.3 INPUT_PUSH=true GITHUB_REF=refs/heads/main
expect "a pull request packages its own tree at the manifest's version and does not push" "$TREE" \
    EVENT=pull_request GITHUB_REF=refs/pull/7/merge
expect "a pull request cannot be made to push or to package a release" "$TREE" \
    EVENT=pull_request GITHUB_REF=refs/heads/main RELEASE_TAG=v1.2.3 INPUT_TAG=v1.2.3 INPUT_PUSH=true

# Every answer over the events, push inputs and refs: one that pushes packages
# a release, and a pull request always packages the tree.
answers=0
unsound=""
for event in release workflow_dispatch pull_request; do
    for push in true false ""; do
        for ref in refs/heads/main refs/heads/topic refs/pull/7/merge ""; do
            out=$(run EVENT="$event" RELEASE_TAG=v1.2.3 INPUT_TAG=v1.2.3 INPUT_PUSH="$push" GITHUB_REF="$ref" 2>/dev/null) || continue
            answers=$((answers + 1))
            [[ "$out" == *"push=true"* && "$out" != *"source=release"* ]] && unsound+=" $event/$push/$ref"
            [[ "$event" == pull_request && "$out" != "$TREE" ]] && unsound+=" $event/$push/$ref"
        done
    done
done
[[ "$answers" -eq 33 && -z "$unsound" ]] && ok "no run that pushes packages the tree" ||
    bad "$answers answers of 33; a pushing run packages the tree or a pull request a release:$unsound"

# A dispatch pairs this tree's image with the release it names; a release
# event runs from its own tag's tree. v0.10.0 sorts after v0.6.0 by version
# and before it as text.
FIRST=$(sed -nE 's/^FIRST_FITTING=(v[0-9]+\.[0-9]+\.[0-9]+)$/\1/p' "$SOURCE")
REFUSAL="ERROR: [image_source] tag v0.5.0 is older than v0.6.0, the first release this image fits: its vep needs the duckdb CLI for Parquet output, and the image carries none"
for push in false true; do
    out=$(run EVENT=workflow_dispatch INPUT_TAG=v0.5.0 INPUT_PUSH="$push" GITHUB_REF=refs/heads/main 2>"$TMP/stderr")
    rc=$?
    [[ "$rc" -eq 1 && -z "$out" && "$(<"$TMP/stderr")" == "$REFUSAL" ]] &&
        ok "a dispatch of a tag older than the first fitting release is refused by one line that names both and the tool (push=$push)" ||
        bad "a dispatch of v0.5.0 with push=$push: rc=$rc, stdout [$out], stderr [$(<"$TMP/stderr")]"
done
expect_fail "a dispatch of the first release's own line is refused below it" "tag v0.5.9 is older than v0.6.0" \
    EVENT=workflow_dispatch INPUT_TAG=v0.5.9 INPUT_PUSH=false GITHUB_REF=refs/heads/topic
# The tags below are at and after the recorded first fitting release.
fitting() { [[ "$FIRST" == v0.6.0 && "$(grep -c 'FIRST_FITTING=' "$SOURCE")" -eq 1 ]]; }
fitting && expect "a dispatch of the first fitting release resolves" $'source=release\ntag=v0.6.0\nversion=0.6.0\npush=true' \
    EVENT=workflow_dispatch INPUT_TAG="$FIRST" INPUT_PUSH=true GITHUB_REF=refs/heads/main ||
    bad "a dispatch of the first fitting release: image_source.sh records [$FIRST], not v0.6.0 once"
fitting && expect "a dispatch of a later patch release resolves" $'source=release\ntag=v0.6.1\nversion=0.6.1\npush=false' \
    EVENT=workflow_dispatch INPUT_TAG=v0.6.1 INPUT_PUSH=false GITHUB_REF=refs/heads/main ||
    bad "a dispatch of a later patch release: image_source.sh records [$FIRST], not v0.6.0 once"
fitting && expect "a dispatch of a later release is compared by version, not as text" $'source=release\ntag=v0.10.0\nversion=0.10.0\npush=false' \
    EVENT=workflow_dispatch INPUT_TAG=v0.10.0 INPUT_PUSH=false GITHUB_REF=refs/heads/main ||
    bad "a dispatch of a later release: image_source.sh records [$FIRST], not v0.6.0 once"
fitting && expect "a release event is never refused for the age of its tag" $'source=release\ntag=v0.5.0\nversion=0.5.0\npush=true' \
    EVENT=release RELEASE_TAG=v0.5.0 GITHUB_REF=refs/tags/v0.5.0 ||
    bad "a release event of an older tag: image_source.sh records [$FIRST], not v0.6.0 once"

expect_fail "a dispatch that asks to push from another branch fails" "push=true is allowed only from main, not refs/heads/topic" \
    EVENT=workflow_dispatch INPUT_TAG=v1.2.3 INPUT_PUSH=true GITHUB_REF=refs/heads/topic
expect_fail "a dispatch that asks to push with no ref fails" "push=true is allowed only from main" \
    EVENT=workflow_dispatch INPUT_TAG=v1.2.3 INPUT_PUSH=true
expect_fail "a dispatch without a tag fails" "tag '' is not vX.Y.Z" EVENT=workflow_dispatch INPUT_PUSH=false
expect_fail "a release with a two-part tag fails" "tag 'v1.2' is not vX.Y.Z" EVENT=release RELEASE_TAG=v1.2
expect_fail "a release with a pre-release tag fails" "tag 'v1.2.3-rc1' is not vX.Y.Z" EVENT=release RELEASE_TAG=v1.2.3-rc1
expect_fail "an event the workflow does not package for fails" "event 'push' is not one the workflow packages an image for" \
    EVENT=push GITHUB_REF=refs/heads/main RELEASE_TAG=v1.2.3 INPUT_TAG=v1.2.3
expect_fail "a run that names no event fails" "EVENT names the event" GITHUB_REF=refs/heads/main
expect_fail "a manifest without a version fails" "states no X.Y.Z version" EVENT=pull_request MANIFEST="$TMP/no-version.toml"
expect_fail "a manifest with a pre-release version fails" "states no X.Y.Z version: '9.9.9-rc1'" \
    EVENT=pull_request MANIFEST="$TMP/pre-release.toml"

# The tree's own manifest, read as release.yml reads it for a dry run.
want=$(sed -nE 's/^version = "([^"]+)"$/\1/p' "$ROOT/Cargo.toml" | head -1)
got=$(env -i PATH="$PATH" EVENT=pull_request bash "$SOURCE" 2>/dev/null)
[[ -n "$want" && "$got" == *$'\n'"version=$want"$'\n'* ]] && ok "the tree's version is the one release.yml names a dry run's archives by" ||
    bad "the tree's version: got [$got], want version=$want"

# The image is built for the tree's own version, so the first fitting release
# is never recorded past it.
[[ -n "$FIRST" && "$(printf '%s\n' "$FIRST" "v$want" | sort -V | head -n 1)" == "$FIRST" ]] &&
    expect "the first fitting release is no later than the tree's version, whose dispatch resolves" \
        "source=release"$'\n'"tag=v$want"$'\n'"version=$want"$'\n'"push=false" \
        EVENT=workflow_dispatch INPUT_TAG="v$want" INPUT_PUSH=false GITHUB_REF=refs/heads/main ||
    bad "the first fitting release [$FIRST] is later than the tree's version v$want"

[[ "$(grep -cxF '        run: scripts/docker/image_source.sh >> "$GITHUB_OUTPUT"' "$WORKFLOW")" -eq 2 ]] &&
    ok "both image jobs of the workflow take their answer from image_source.sh" ||
    bad "the workflow does not write image_source.sh's lines to a step's outputs in both image jobs"

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
