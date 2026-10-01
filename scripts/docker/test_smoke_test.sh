#!/usr/bin/env bash
# Tests for smoke_test.sh: a passing run against stub binaries, then one
# injected fault per check proving that check fails the run.
#
# Pure bash. No image, no Docker, no cache.
# Run: bash scripts/docker/test_smoke_test.sh
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
SMOKE="$HERE/smoke_test.sh"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

STUBS="$TMP/bin"
CORPUS="$TMP/corpus"
CA="$TMP/ca.crt"
mkdir -p "$STUBS" "$CORPUS/json_cache"
echo "cert" >"$CA"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >"$CORPUS/variants.vcf"
for i in 1 2 3; do printf '21\t%s\t.\tA\tG\t.\t.\t.\n' "$i" >>"$CORPUS/variants.vcf"; done

cat >"$STUBS/vep" <<'STUB'
#!/usr/bin/env bash
case "$1" in
    --version) echo "vep ${STUB_VERSION:-9.9.9}"; exit 0 ;;
    --help) exit 0 ;;
esac
out="" fmt=vep
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) out="$2"; shift ;;
        --vcf) fmt=vcf ;;
        --tab) fmt=tab ;;
        --output_format) fmt="$2"; shift ;;
    esac
    shift
done
case "$fmt" in
    vcf) {
        echo '##fileformat=VCFv4.2'
        for ((i = 0; i < ${STUB_VCF_RECORDS:-3}; i++)); do
            if ((i < ${STUB_VCF_CSQ:-3})); then echo "21 $i CSQ=G|missense_variant"; else echo "21 $i ."; fi
        done
    } >"$out" ;;
    tab) { echo '#Uploaded_variation'; for ((i = 0; i < ${STUB_TAB_ROWS:-5}; i++)); do echo "r$i"; done; } >"$out" ;;
    parquet) mkdir -p "$out" ;;
esac
STUB
cat >"$STUBS/duckdb" <<'STUB'
#!/usr/bin/env bash
[[ "$1" == --version ]] && { echo v1.5.5; exit "${STUB_DUCKDB_RC:-0}"; }
[[ "${STUB_PQ_RC:-0}" -eq 0 ]] || { echo "IO Error: No files found" >&2; exit 1; }
echo "${STUB_PQ_ROWS:-5}"
STUB
cat >"$STUBS/vep-cache-builder" <<'STUB'
#!/usr/bin/env bash
exit "${STUB_BUILDER_RC:-0}"
STUB
cat >"$STUBS/vep-cache-converter" <<'STUB'
#!/usr/bin/env bash
exit "${STUB_CONVERTER_RC:-0}"
STUB
cat >"$STUBS/id" <<'STUB'
#!/usr/bin/env bash
[[ "$1" == -u ]] && { echo "${STUB_UID:-1000}"; exit 0; }
exec /usr/bin/id "$@"
STUB
chmod +x "$STUBS"/*

smoke() { # smoke [VAR=value ...]: runs smoke_test.sh against the stubs
    env PATH="$STUBS:$PATH" CORPUS="$CORPUS" CA_BUNDLE="$CA" "$@" bash "$SMOKE" 9.9.9
}

expect_fail() { # expect_fail <name> <stderr substring> [VAR=value ...]
    local name="$1" want="$2" err rc
    shift 2
    err=$(smoke "$@" 2>&1 >/dev/null)
    rc=$?
    if [[ "$rc" -ne 0 && "$err" == *"$want"* ]]; then ok "$name"; else bad "$name (rc=$rc, stderr: $err)"; fi
}

bash -n "$SMOKE" && ok "smoke_test.sh parses" || bad "smoke_test.sh has a syntax error"

out=$(smoke 2>&1)
rc=$?
n=$(printf '%s\n' "$out" | grep -c '^  ok ') || n=0
[[ "$rc" -eq 0 && "$n" -eq 10 ]] && ok "a healthy image passes all 10 checks" || bad "healthy run: rc=$rc, $n ok lines: $out"

expect_fail "a wrong vep version fails" "vep --version printed 'vep 1.0.0'" STUB_VERSION=1.0.0
expect_fail "a broken vep-cache-builder fails" "vep-cache-builder --help exited non-zero" STUB_BUILDER_RC=1
expect_fail "a broken vep-cache-converter fails" "vep-cache-converter --help exited non-zero" STUB_CONVERTER_RC=1
expect_fail "a broken duckdb fails" "duckdb --version exited non-zero" STUB_DUCKDB_RC=1
expect_fail "a missing CA bundle fails" "no CA bundle" CA_BUNDLE="$TMP/absent.crt"
expect_fail "running as root fails" "runs as root" STUB_UID=0
expect_fail "a dropped VCF record fails" "--vcf wrote 2 records for 3" STUB_VCF_RECORDS=2
expect_fail "VCF records without CSQ fail" "2 of 3 --vcf records carry CSQ" STUB_VCF_CSQ=2
expect_fail "an empty tab output fails" "--tab wrote no rows" STUB_TAB_ROWS=0
expect_fail "a Parquet row-count mismatch fails" "Parquet has 4 rows, --tab has 5" STUB_PQ_ROWS=4
expect_fail "an unreadable Parquet output fails" "duckdb could not read the Parquet output" STUB_PQ_RC=1

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
