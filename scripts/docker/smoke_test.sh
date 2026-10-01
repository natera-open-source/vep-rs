#!/usr/bin/env bash
# Runs inside the vep-rs image and exits 1 on the first failed check: every
# binary starts, the CA bundle is present, the default user is not root, and
# the bundled GRCh37 corpus annotates in VCF, tab and Parquet.
#
# From the repository root:
#   docker run --rm -v "$PWD/scripts/docker:/smoke:ro" \
#     -v "$PWD/tests/golden/115/GRCh37:/corpus:ro" IMAGE bash /smoke/smoke_test.sh X.Y.Z
set -euo pipefail

expected_version="${1:?usage: smoke_test.sh <expected version>}"
CORPUS="${CORPUS:-/corpus}"
CA_BUNDLE="${CA_BUNDLE:-/etc/ssl/certs/ca-certificates.crt}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() {
    echo "ERROR: [smoke_test] $*" >&2
    exit 1
}
pass() { echo "  ok   $*"; }

got="$(vep --version)"
[[ "$got" == "vep $expected_version" ]] || fail "vep --version printed '$got', expected 'vep $expected_version'"
pass "vep --version is $expected_version"

for bin in vep-cache-builder vep-cache-converter; do
    "$bin" --help >/dev/null || fail "$bin --help exited non-zero"
    pass "$bin starts"
done

duckdb --version >/dev/null || fail "duckdb --version exited non-zero"
pass "duckdb starts"

[[ -s "$CA_BUNDLE" ]] || fail "no CA bundle at $CA_BUNDLE"
pass "CA bundle present"

[[ "$(id -u)" != 0 ]] || fail "the image runs as root"
pass "default user is not root"

common=(-i "$CORPUS/variants.vcf" --offline --json_cache "$CORPUS/json_cache"
    --assembly GRCh37 --force --no_stats --quiet)

vep "${common[@]}" --vcf -o "$WORK/out.vcf"
records_in=$(grep -vc '^#' "$CORPUS/variants.vcf") || records_in=0
records_out=$(grep -vc '^#' "$WORK/out.vcf") || records_out=0
[[ "$records_in" -gt 0 && "$records_out" -eq "$records_in" ]] ||
    fail "--vcf wrote $records_out records for $records_in input records"
pass "--vcf writes every input record ($records_in)"

vep "${common[@]}" --tab -o "$WORK/out.tab"
tab_rows=$(grep -vc '^#' "$WORK/out.tab") || tab_rows=0
[[ "$tab_rows" -gt 0 ]] || fail "--tab wrote no rows"
pass "--tab writes $tab_rows rows"

vep "${common[@]}" --output_format parquet --parquet_shape flat -o "$WORK/out.parquet"
pq_rows=$(duckdb -noheader -list -c \
    "SELECT count(*) FROM read_parquet('$WORK/out.parquet/**/*.parquet', hive_partitioning=false)")
[[ "$pq_rows" == "$tab_rows" ]] || fail "Parquet has $pq_rows rows, --tab has $tab_rows"
pass "Parquet row count matches --tab"
