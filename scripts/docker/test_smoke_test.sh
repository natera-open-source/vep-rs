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
EXAMPLE="$TMP/examples/homo_sapiens_GRCh37.vcf"
CA="$TMP/ca.crt"
CONF="$TMP/ca-certificates.conf"
mkdir -p "$STUBS" "$CORPUS/json_cache" "$TMP/examples"
echo "cert" >"$CA"
# The certificate list of a healthy store: the two roots deactivated; a list
# with one of them active, and one with both.
printf '%s\n' 'mozilla/ISRG_Root_X1.crt' '!mozilla/emSign_ECC_Root_CA_-_G3.crt' '!mozilla/Comodo_AAA_Services_root.crt' >"$CONF"
printf '%s\n' 'mozilla/ISRG_Root_X1.crt' '!mozilla/emSign_ECC_Root_CA_-_G3.crt' 'mozilla/Comodo_AAA_Services_root.crt' >"$TMP/one-active.conf"
printf '%s\n' 'mozilla/emSign_ECC_Root_CA_-_G3.crt' 'mozilla/Comodo_AAA_Services_root.crt' >"$TMP/both-active.conf"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >"$CORPUS/variants.vcf"
for i in 1 2 3; do printf '21\t%s\t.\tA\tG\t.\t.\t.\n' "$i" >>"$CORPUS/variants.vcf"; done
# The example: two records with identifiers of their own.
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >"$EXAMPLE"
for i in 1 2; do printf '21\t%s\tid%s\tA\tG\t.\t.\t.\n' "$i" "$i" >>"$EXAMPLE"; done

cat >"$STUBS/vep" <<'STUB'
#!/usr/bin/env bash
case "$1" in
    --version) echo "vep ${STUB_VERSION:-9.9.9}"; exit 0 ;;
    --help) exit 0 ;;
esac
out="" fmt=vep input="" warnings=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) input="$2"; shift ;;
        -o) out="$2"; shift ;;
        --warning_file) warnings="$2"; shift ;;
        --vcf) fmt=vcf ;;
        --tab) fmt=tab ;;
        --output_format) fmt="$2"; shift ;;
    esac
    shift
done
case "$fmt" in
    vep)
        # The default format: one row per input record keyed by its ID, the
        # last record dropped under STUB_EXAMPLE_DROP=1; STUB_EXAMPLE_WARN=1
        # writes a warnings file beside the output.
        {
            echo '#Uploaded_variation'
            grep -v '^#' "$input" | cut -f3 | { if [[ "${STUB_EXAMPLE_DROP:-0}" -eq 1 ]]; then sed '$d'; else cat; fi; }
        } >"$out"
        if [[ "${STUB_EXAMPLE_WARN:-0}" -eq 1 && -n "$warnings" ]]; then
            echo "WARNING: line 2 skipped" >"$warnings"
        fi
        ;;
    vcf) {
        echo '##fileformat=VCFv4.2'
        for ((i = 0; i < ${STUB_VCF_RECORDS:-3}; i++)); do
            if ((i < ${STUB_VCF_CSQ:-3})); then echo "21 $i CSQ=G|missense_variant"; else echo "21 $i ."; fi
        done
    } >"$out" ;;
    tab) { echo '#Uploaded_variation'; for ((i = 0; i < ${STUB_TAB_ROWS:-5}; i++)); do echo "r$i"; done; } >"$out" ;;
    parquet)
        # A Parquet-shaped file: the magic bytes, a body, a footer naming the writer's
        # version key, the magic bytes again. STUB_PQ_PARTS part files in a directory
        # (default 2; 0 leaves the directory empty), under chrom=21/ with STUB_PQ_HIVE=1,
        # or one file at the output path with STUB_PQ_SINGLE_FILE=1. STUB_PQ_TRUNCATE=1
        # drops the closing magic of the last part, STUB_PQ_NO_VERSION=1 the footer key.
        part() {
            local footer="vep_rs_version=9.9.9"
            [[ "${STUB_PQ_NO_VERSION:-0}" -eq 1 ]] && footer="created_by=stub"
            if [[ "${STUB_PQ_TRUNCATE:-0}" -eq 1 && "$2" -eq 1 ]]; then
                printf 'PAR1rows%s' "$footer" >"$1"
            else
                printf 'PAR1rows%sPAR1' "$footer" >"$1"
            fi
        }
        if [[ "${STUB_PQ_SINGLE_FILE:-0}" -eq 1 ]]; then
            part "$out" 1
        else
            dir="$out"
            [[ "${STUB_PQ_HIVE:-0}" -eq 1 ]] && dir="$out/chrom=21"
            mkdir -p "$dir"
            n=${STUB_PQ_PARTS:-2}
            for ((i = 0; i < n; i++)); do part "$dir/part-0000$i.parquet" $((i == n - 1)); done
        fi
        ;;
esac
STUB
cat >"$STUBS/bgzip" <<'STUB'
#!/usr/bin/env bash
[[ "$1" == --version ]] && { echo "bgzip (htslib) 1.21"; exit "${STUB_BGZIP_RC:-0}"; }
STUB
cat >"$STUBS/tabix" <<'STUB'
#!/usr/bin/env bash
[[ "$1" == --version ]] && { echo "tabix (htslib) 1.21"; exit "${STUB_TABIX_RC:-0}"; }
STUB
cat >"$STUBS/ps" <<'STUB'
#!/usr/bin/env bash
[[ "$*" == "-e -o pid= -o ppid=" ]] && { echo "1 0"; exit "${STUB_PS_RC:-0}"; }
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
    env PATH="$STUBS:$PATH" CORPUS="$CORPUS" EXAMPLE="$EXAMPLE" CA_BUNDLE="$CA" CA_CONF="$CONF" "$@" bash "$SMOKE" 9.9.9
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
[[ "$rc" -eq 0 && "$n" -eq 14 ]] && ok "a healthy image passes all 14 checks" || bad "healthy run: rc=$rc, $n ok lines: $out"
[[ "$out" == *"Parquet wrote 2 well-formed file(s)"* ]] && ok "every part file of a Parquet directory is checked" || bad "part count: $out"
[[ "$out" == *"the bundled example annotates every record (2) with no warning"* ]] && ok "the bundled example's records are counted" || bad "example count: $out"
[[ "$out" == *"the emSign ECC G3 and Comodo AAA Services roots are inactive"* ]] && ok "a store with both roots deactivated passes" || bad "roots check: $out"

# A release tree without the example mounts none: the check is skipped by name
# and the other 13 pass.
out=$(smoke EXAMPLE="$TMP/absent/homo_sapiens_GRCh37.vcf" 2>&1)
rc=$?
n=$(printf '%s\n' "$out" | grep -c '^  ok ') || n=0
[[ "$rc" -eq 0 && "$n" -eq 13 && "$out" == *"skip no bundled example at $TMP/absent/homo_sapiens_GRCh37.vcf"* ]] &&
    ok "a missing example is skipped by name" || bad "missing example: rc=$rc, $n ok lines: $out"

expect_pass() { # expect_pass <name> <stdout substring> [VAR=value ...]
    local name="$1" want="$2" out rc
    shift 2
    out=$(smoke "$@" 2>&1)
    rc=$?
    if [[ "$rc" -eq 0 && "$out" == *"$want"* ]]; then ok "$name"; else bad "$name (rc=$rc, output: $out)"; fi
}
expect_pass "a single-file Parquet output passes" "Parquet wrote 1 well-formed file(s)" STUB_PQ_SINGLE_FILE=1
expect_pass "a Parquet directory partitioned by chromosome passes" "Parquet wrote 2 well-formed file(s)" STUB_PQ_HIVE=1

expect_fail "a wrong vep version fails" "vep --version printed 'vep 1.0.0'" STUB_VERSION=1.0.0
expect_fail "a broken vep-cache-builder fails" "vep-cache-builder --help exited non-zero" STUB_BUILDER_RC=1
expect_fail "a broken vep-cache-converter fails" "vep-cache-converter --help exited non-zero" STUB_CONVERTER_RC=1
expect_fail "a broken bgzip fails" "bgzip --version exited non-zero" STUB_BGZIP_RC=1
expect_fail "a broken tabix fails" "tabix --version exited non-zero" STUB_TABIX_RC=1
expect_fail "a broken ps fails" "ps -e -o pid= -o ppid= exited non-zero" STUB_PS_RC=1
# An image without ps: a PATH of every stub but ps, plus the host commands
# smoke_test.sh runs before that check and on exit.
NO_PS="$TMP/no-ps"
mkdir "$NO_PS"
for tool in "$STUBS"/* "$(command -v bash)" "$(command -v mktemp)" "$(command -v rm)"; do
    [[ "${tool##*/}" == ps ]] || ln -s "$tool" "$NO_PS/"
done
expect_fail "a missing ps fails" "ps -e -o pid= -o ppid= exited non-zero" PATH="$NO_PS"
expect_fail "a missing CA bundle fails" "no CA bundle" CA_BUNDLE="$TMP/absent.crt"
expect_fail "a missing certificate list fails" "no certificate list at $TMP/absent.conf" CA_CONF="$TMP/absent.conf"
expect_fail "one active root fails" "activates 1 of the emSign ECC G3 and Comodo AAA Services roots" CA_CONF="$TMP/one-active.conf"
expect_fail "both active roots fail" "activates 2 of the emSign ECC G3 and Comodo AAA Services roots" CA_CONF="$TMP/both-active.conf"
expect_fail "running as root fails" "runs as root" STUB_UID=0
expect_fail "a dropped VCF record fails" "--vcf wrote 2 records for 3" STUB_VCF_RECORDS=2
expect_fail "VCF records without CSQ fail" "2 of 3 --vcf records carry CSQ" STUB_VCF_CSQ=2
expect_fail "an example record without a row fails" "the bundled example's 2 records got rows for 1" STUB_EXAMPLE_DROP=1
expect_fail "a warning while annotating the example fails" "annotating the bundled example warned: WARNING: line 2 skipped" STUB_EXAMPLE_WARN=1
expect_fail "an empty tab output fails" "--tab wrote no rows" STUB_TAB_ROWS=0
expect_fail "a Parquet directory without part files fails" "wrote no Parquet file" STUB_PQ_PARTS=0
expect_fail "a truncated Parquet part fails" "part-00001.parquet does not open and close with the Parquet magic bytes" STUB_PQ_TRUNCATE=1
expect_fail "a Parquet part without the writer's version fails" "part-00000.parquet carries no vep_rs_version in its footer" STUB_PQ_NO_VERSION=1

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
