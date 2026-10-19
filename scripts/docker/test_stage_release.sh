#!/usr/bin/env bash
# Tests for stage_release.sh: a stub `gh` serves a fabricated release with two
# Linux archives, the target argument picks the archive, and each fail-closed
# gate (tag shape, target shape, existing destination, download, checksum line,
# checksum value, attestation) is shown to stop the run with nothing staged.
#
# Pure bash. No network, no GitHub.
# Run: bash scripts/docker/test_stage_release.sh
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
STAGE="$HERE/stage_release.sh"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

ARCHIVE=vep-9.9.9-x86_64-unknown-linux-gnu.tar.gz
GENERIC=vep-9.9.9-aarch64-unknown-linux-gnu-generic.tar.gz
RELEASE="$TMP/release"
STUBS="$TMP/bin"
mkdir -p "$RELEASE" "$STUBS" "$TMP/src/vep-9.9.9-x86_64-unknown-linux-gnu" "$TMP/src/vep-9.9.9-aarch64-unknown-linux-gnu-generic"
for bin in vep vep-cache-builder vep-cache-converter; do
    echo "$bin" >"$TMP/src/vep-9.9.9-x86_64-unknown-linux-gnu/$bin"
    echo "$bin generic" >"$TMP/src/vep-9.9.9-aarch64-unknown-linux-gnu-generic/$bin"
done
tar -czf "$RELEASE/$ARCHIVE" -C "$TMP/src" vep-9.9.9-x86_64-unknown-linux-gnu
tar -czf "$RELEASE/$GENERIC" -C "$TMP/src" vep-9.9.9-aarch64-unknown-linux-gnu-generic
sum() { (cd "$RELEASE" && sha256sum "$ARCHIVE" "$GENERIC"); }
good_sums() { { echo "$(printf '0%.0s' {1..64})  vep-rs-9.9.9.tar.gz"; sum; } >"$RELEASE/SHA256SUMS"; }
good_sums

cat >"$STUBS/gh" <<'STUB'
#!/usr/bin/env bash
echo "$*" >>"$STUB_LOG"
case "$1 $2" in
    "release download")
        [[ "${STUB_DOWNLOAD_RC:-0}" -eq 0 ]] || exit "$STUB_DOWNLOAD_RC"
        dir="" pats=()
        shift 3
        while [[ $# -gt 0 ]]; do
            case "$1" in
                --dir) dir="$2"; shift ;;
                --pattern) pats+=("$2"); shift ;;
            esac
            shift
        done
        for p in "${pats[@]}"; do cp "$STUB_RELEASE/$p" "$dir/"; done ;;
    "attestation verify") exit "${STUB_ATTEST_RC:-0}" ;;
    *) exit 2 ;;
esac
STUB
chmod +x "$STUBS/gh"

run() { # run <dest> [VAR=value ...]; TAG and TARGET in the caller's env override v9.9.9 and the x86_64 archive
    local dest="$1"
    shift
    env PATH="$STUBS:$PATH" STUB_RELEASE="$RELEASE" STUB_LOG="$TMP/gh.log" "$@" \
        bash "$STAGE" "${TAG:-v9.9.9}" "${TARGET:-x86_64-unknown-linux-gnu}" "$dest"
}

expect_fail() { # expect_fail <name> <stderr substring> [VAR=value ...]
    local name="$1" want="$2" dest="$TMP/ctx-$RANDOM" err rc
    shift 2
    err=$(run "$dest" "$@" 2>&1 >/dev/null)
    rc=$?
    if [[ "$rc" -ne 0 && "$err" == *"$want"* && ! -e "$dest/vep" ]]; then
        ok "$name"
    else
        bad "$name (rc=$rc, staged=$([[ -e "$dest/vep" ]] && echo yes || echo no), stderr: $err)"
    fi
}

bash -n "$STAGE" && ok "stage_release.sh parses" || bad "stage_release.sh has a syntax error"

: >"$TMP/gh.log"
dest="$TMP/ctx-good"
if run "$dest" >/dev/null 2>&1 &&
    [[ -x "$dest/vep" && -x "$dest/vep-cache-builder" && -x "$dest/vep-cache-converter" &&
        -s "$dest/LICENSE" && -s "$dest/NOTICE" && "$(cat "$dest/vep")" == vep ]]; then
    ok "a verified release stages three executables, LICENSE and NOTICE"
else
    bad "a verified release did not stage the context: $(ls -l "$dest" 2>&1)"
fi
grep -q -- "--signer-workflow natera-open-source/vep-rs/.github/workflows/release.yml" "$TMP/gh.log" &&
    ok "the attestation must come from release.yml" ||
    bad "gh attestation verify was not pinned to release.yml: $(cat "$TMP/gh.log")"
grep -q -- "--source-ref refs/tags/v9.9.9" "$TMP/gh.log" &&
    ok "the attestation must come from this release's tag" ||
    bad "gh attestation verify was not pinned to refs/tags/v9.9.9: $(cat "$TMP/gh.log")"

dest="$TMP/ctx-generic"
if TARGET=aarch64-unknown-linux-gnu-generic run "$dest" >/dev/null 2>&1 &&
    [[ "$(cat "$dest/vep")" == "vep generic" && "$(cat "$dest/vep-cache-converter")" == "vep-cache-converter generic" ]]; then
    ok "the target argument picks that archive of the release"
else
    bad "the generic target did not stage the generic archive: $(cat "$dest/vep" 2>&1)"
fi
grep -q -- "--pattern $GENERIC" "$TMP/gh.log" &&
    ok "only the named archive and SHA256SUMS are downloaded" ||
    bad "the generic archive was not the download pattern: $(cat "$TMP/gh.log")"

TAG=0.2.0 expect_fail "a tag without the v prefix fails" "is not vX.Y.Z"
TAG=v1.0 expect_fail "a two-part tag fails" "is not vX.Y.Z"
TARGET=aarch64-apple-darwin expect_fail "a non-Linux target fails" "is not a Linux archive name"
TARGET=x86_64-unknown-linux-gnu/../x expect_fail "a target with a path in it fails" "is not a Linux archive name"
TARGET=aarch64-unknown-linux-gnu-sve2 expect_fail "a target whose archive the release lacks fails" "SHA256SUMS has no line for vep-9.9.9-aarch64-unknown-linux-gnu-sve2.tar.gz"
mkdir -p "$TMP/exists"
err=$(run "$TMP/exists" 2>&1 >/dev/null) && bad "an existing destination was overwritten" ||
    { [[ "$err" == *"already exists"* ]] && ok "an existing destination is never overwritten" || bad "existing destination: $err"; }
expect_fail "a failed download fails" "could not download" STUB_DOWNLOAD_RC=1

echo "$(printf '0%.0s' {1..64})  vep-rs-9.9.9.tar.gz" >"$RELEASE/SHA256SUMS"
expect_fail "SHA256SUMS without the archive's line fails" "has no line for $ARCHIVE"
echo "$(printf 'a%.0s' {1..64})  $ARCHIVE" >"$RELEASE/SHA256SUMS"
expect_fail "a checksum mismatch fails" "does not match SHA256SUMS"
good_sums
expect_fail "a missing attestation fails" "has no attestation" STUB_ATTEST_RC=1

rm "$TMP/src/vep-9.9.9-x86_64-unknown-linux-gnu/vep-cache-converter"
tar -czf "$RELEASE/$ARCHIVE" -C "$TMP/src" vep-9.9.9-x86_64-unknown-linux-gnu
good_sums
dest="$TMP/ctx-partial"
err=$(run "$dest" 2>&1 >/dev/null) && bad "an archive without vep-cache-converter was staged" ||
    { [[ "$err" == *"has no vep-cache-converter"* && ! -e "$dest" ]] && ok "an incomplete archive leaves nothing staged" ||
        bad "incomplete archive: staged=$([[ -e "$dest" ]] && echo yes || echo no), stderr: $err"; }

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
