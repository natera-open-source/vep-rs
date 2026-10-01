#!/usr/bin/env bash
# Tests for image_tags.sh: which tags a release pushes, given the published
# releases, and the inputs it refuses.
#
# Pure bash. Run: bash scripts/docker/test_image_tags.sh
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
TAGS="$HERE/image_tags.sh"
IMG=ghcr.io/o/r
PUBLISHED=$'v0.1.0\nv0.2.0\nv0.2.1\nv0.10.0\nv0.10.1-rc1\nnot-a-version'

expect() { # expect <name> <tag> <expected lines>
    local got
    got=$(bash "$TAGS" "$IMG" "$2" <<<"$PUBLISHED" 2>&1)
    [[ "$got" == "$3" ]] && ok "$1" || bad "$1: got [$got], want [$3]"
}
expect_fail() { # expect_fail <name> <tag> <stderr substring>
    local err rc
    err=$(bash "$TAGS" "$IMG" "$2" <<<"$PUBLISHED" 2>&1 >/dev/null)
    rc=$?
    [[ "$rc" -ne 0 && "$err" == *"$3"* ]] && ok "$1" || bad "$1 (rc=$rc, stderr: $err)"
}

bash -n "$TAGS" && ok "image_tags.sh parses" || bad "image_tags.sh has a syntax error"
expect "the newest release gets X.Y.Z, X.Y and latest" v0.10.0 $'ghcr.io/o/r:0.10.0\nghcr.io/o/r:0.10\nghcr.io/o/r:latest'
expect "the newest of an older line gets X.Y.Z and X.Y" v0.2.1 $'ghcr.io/o/r:0.2.1\nghcr.io/o/r:0.2'
expect "an older patch gets only X.Y.Z" v0.2.0 'ghcr.io/o/r:0.2.0'
expect "0.10 sorts above 0.2 (version order, not text order)" v0.1.0 $'ghcr.io/o/r:0.1.0\nghcr.io/o/r:0.1'
expect_fail "an unpublished tag fails" v0.3.0 "is not a published release"
expect_fail "a pre-release tag fails" v0.10.1-rc1 "is not vX.Y.Z"
expect_fail "a tag without the v prefix fails" 0.2.0 "is not vX.Y.Z"

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
