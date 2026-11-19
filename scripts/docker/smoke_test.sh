#!/usr/bin/env bash
# Runs inside the vep-rs image and exits 1 on the first failed check: every
# binary starts, ps lists the processes, the CA bundle is present, the two
# roots the store must not activate are inactive, the default user is not root,
# the bundled GRCh37 corpus annotates in VCF, tab and Parquet, and the bundled
# GRCh37 example annotates every record against the corpus cache when the
# release's examples directory is mounted.
#
# From the repository root:
#   docker run --rm -v "$PWD/scripts/docker:/smoke:ro" \
#     -v "$PWD/tests/golden/116/GRCh37:/corpus:ro" \
#     -v "$PWD/examples:/examples:ro" IMAGE bash /smoke/smoke_test.sh X.Y.Z
set -euo pipefail

expected_version="${1:?usage: smoke_test.sh <expected version>}"
CORPUS="${CORPUS:-/corpus}"
EXAMPLE="${EXAMPLE:-/examples/homo_sapiens_GRCh37.vcf}"
CA_BUNDLE="${CA_BUNDLE:-/etc/ssl/certs/ca-certificates.crt}"
CA_CONF="${CA_CONF:-/etc/ca-certificates.conf}"
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

for bin in bgzip tabix; do
    "$bin" --version >/dev/null || fail "$bin --version exited non-zero"
    pass "$bin starts"
done

# The call Nextflow's task wrapper makes to walk a task's process tree.
ps -e -o pid= -o ppid= >/dev/null || fail "ps -e -o pid= -o ppid= exited non-zero"
pass "ps lists the processes"

[[ -s "$CA_BUNDLE" ]] || fail "no CA bundle at $CA_BUNDLE"
pass "CA bundle present"

# ca-certificates.conf lists the store's roots, a deactivated one prefixed
# with `!`. Neither of these two roots may be active, whichever ca-certificates
# version the base image installs; a missing list fails instead of reading as
# no active root.
[[ -s "$CA_CONF" ]] || fail "no certificate list at $CA_CONF"
active=$(grep -cE '^(mozilla/emSign_ECC_Root_CA_-_G3|mozilla/Comodo_AAA_Services_root)\.crt$' "$CA_CONF") || active=0
[[ "$active" -eq 0 ]] || fail "$CA_CONF activates $active of the emSign ECC G3 and Comodo AAA Services roots"
pass "the emSign ECC G3 and Comodo AAA Services roots are inactive"

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

# 1,002 of the corpus's 1,008 records carry a consequence; a cache that failed
# to load echoes the records back without one.
annotated=$(grep -v '^#' "$WORK/out.vcf" | grep -c 'CSQ=') || annotated=0
((annotated * 10 >= records_in * 9)) ||
    fail "$annotated of $records_in --vcf records carry CSQ, expected at least 90%"
pass "--vcf annotates $annotated of $records_in records"

# The example is a subset of the corpus's records, so the corpus cache covers
# every one: each record gets at least one row and the run warns of nothing. A
# release up to 0.5.0 has no examples directory and mounts no example, and the
# check says so instead of failing the image.
if [[ -f "$EXAMPLE" ]]; then
    vep -i "$EXAMPLE" --offline --json_cache "$CORPUS/json_cache" --assembly GRCh37 \
        --force --no_stats --quiet --warning_file "$WORK/example_warnings.txt" \
        -o "$WORK/example.txt"
    example_in=$(grep -vc '^#' "$EXAMPLE") || example_in=0
    example_annotated=$(grep -v '^#' "$WORK/example.txt" | cut -f1 | sort -u | grep -c .) || example_annotated=0
    [[ "$example_in" -gt 0 && "$example_annotated" -eq "$example_in" ]] ||
        fail "the bundled example's $example_in records got rows for $example_annotated"
    [[ ! -s "$WORK/example_warnings.txt" ]] ||
        fail "annotating the bundled example warned: $(head -n 1 "$WORK/example_warnings.txt")"
    pass "the bundled example annotates every record ($example_in) with no warning"
else
    echo "  skip no bundled example at $EXAMPLE"
fi

vep "${common[@]}" --tab -o "$WORK/out.tab"
tab_rows=$(grep -vc '^#' "$WORK/out.tab") || tab_rows=0
[[ "$tab_rows" -gt 0 ]] || fail "--tab wrote no rows"
pass "--tab writes $tab_rows rows"

vep "${common[@]}" --output_format parquet -o "$WORK/out.parquet"
# The output is one file or a directory of part files. The image carries no
# Parquet reader, so each file is checked for its shape: it opens and closes
# with the format's magic bytes and its footer names the writer's version key;
# the rows themselves are compared with the --tab rows by the repository's
# tests.
parts=()
if [[ -f "$WORK/out.parquet" ]]; then
    parts=("$WORK/out.parquet")
elif [[ -d "$WORK/out.parquet" ]]; then
    while IFS= read -r -d '' part; do parts+=("$part"); done \
        < <(find "$WORK/out.parquet" -type f -name 'part-*.parquet' -print0 | sort -z)
fi
[[ ${#parts[@]} -gt 0 ]] || fail "--output_format parquet wrote no Parquet file at $WORK/out.parquet"
for part in "${parts[@]}"; do
    [[ "$(head -c 4 "$part")" == PAR1 && "$(tail -c 4 "$part")" == PAR1 ]] ||
        fail "$part does not open and close with the Parquet magic bytes"
    versions=$(tail -c 65536 "$part" | grep -a -c vep_rs_version) || versions=0
    [[ "$versions" -ge 1 ]] || fail "$part carries no vep_rs_version in its footer"
done
pass "Parquet wrote ${#parts[@]} well-formed file(s)"
