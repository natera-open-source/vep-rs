#!/usr/bin/env bash
# Tests for feature_oracle.sh: the case-table contract, the dry-run command shapes, every
# fail-closed path, and the verdicts on the real path through a stub `docker` and a stub
# `vep`, both of which answer with the committed Ensembl VEP outputs of
# tests/golden/116/GRCh38-release116 so the real comparator decides.
#
# Pure bash. No Docker, no VEP, no cache.
#
# Run: bash scripts/concordance/test_feature_oracle.sh

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
ORACLE="$HERE/feature_oracle.sh"
REPO=$(cd "$HERE/../.." && pwd)
CORPUS="$REPO/tests/golden/116/GRCh38-release116"
IMAGE="ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [[ ! -f "$ORACLE" ]]; then
    bad "feature_oracle.sh not found beside this test"
    echo "passed: $PASS   failed: $FAIL"
    exit 1
fi
bash -n "$ORACLE" && ok "feature_oracle.sh parses" || bad "feature_oracle.sh has a syntax error"

# The committed reference outputs, decompressed once: what both stubs answer with.
REF="$TMP/ref"
mkdir -p "$REF"
for f in default.txt tab.txt vcf.vcf json.jsonl; do
    gzip -dc "$CORPUS/expected/$f.gz" >"$REF/$f"
done

# A fabricated data tree: an Ensembl cache slice, a candidate JSON cache per assembly.
DATA="$TMP/data"
mkdir -p "$DATA/caches/perl/vep-cache/homo_sapiens/116_GRCh38" "$DATA/caches/perl/vep-cache/homo_sapiens/116_GRCh37" \
    "$DATA/caches/vep-rs/grch38/transcripts" "$DATA/caches/vep-rs/grch37/transcripts"

# Stubs. `docker run` and `vep` copy the reference output named by -o from $STUB_REF_DIR,
# write a warnings companion beside a default run when STUB_PERL_WARNINGS / STUB_CAND_WARNINGS
# is set, mutate the candidate's default output per STUB_CAND_MUTATE, exit STUB_PERL_RC /
# STUB_CAND_RC, and append their argv to STUB_DOCKER_LOG / STUB_CAND_LOG.
STUBDIR="$TMP/stubbin"
mkdir -p "$STUBDIR"
cat >"$STUBDIR/docker" <<'EOF'
#!/usr/bin/env bash
[[ -n "${STUB_DOCKER_LOG:-}" ]] && printf '%s\n' "$*" >>"$STUB_DOCKER_LOG"
case "$1" in
image)
    # STUB_PULL_MARKER names a file `docker pull` creates: until then the image is absent.
    [[ -n "${STUB_PULL_MARKER:-}" && ! -e "$STUB_PULL_MARKER" ]] && exit 0
    echo "ensemblorg/ensembl-vep@sha256:stubdigest"; exit 0 ;;
pull) [[ -z "${STUB_PULL_MARKER:-}" ]] || : >"$STUB_PULL_MARKER"; exit 0 ;;
run) ;;
*) echo "stub docker: unexpected subcommand $1" >&2; exit 99 ;;
esac
host_out=""
out=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    -v) [[ "$2" == *:/work/out ]] && host_out="${2%%:*}"; shift 2 ;;
    -o) out="$2"; shift 2 ;;
    *) shift ;;
    esac
done
[[ -n "$host_out" && -n "$out" ]] || { echo "stub docker: no output mapping" >&2; exit 98; }
name=$(basename "$out")
target="$host_out/$name"
if [[ -n "${STUB_PERL_SLEEP:-}" ]]; then
    [[ -z "${STUB_PERL_PIDFILE:-}" ]] || echo "$$" >>"$STUB_PERL_PIDFILE"
    if [[ -n "${STUB_PERL_TIMES:-}" ]]; then
        start=$(date +%s.%N)
        sleep "$STUB_PERL_SLEEP"
        echo "$start $(date +%s.%N)" >>"$STUB_PERL_TIMES"
        exit 0
    fi
    exec sleep "$STUB_PERL_SLEEP"
fi
[[ "${STUB_PERL_RC:-0}" -eq 0 ]] || exit "${STUB_PERL_RC}"
[[ -z "${STUB_PERL_NOFILE:-}" ]] || exit 0
if [[ -n "${STUB_PERL_EMPTY:-}" ]]; then
    grep '^#' "$STUB_REF_DIR/$name" >"$target"
    exit 0
fi
cp "$STUB_REF_DIR/$name" "$target"
if [[ "$name" == "default.txt" && -n "${STUB_PERL_WARNINGS:-}" ]]; then
    printf '%s\n' "$STUB_PERL_WARNINGS" >"${target}_warnings.txt"
fi
exit 0
EOF
cat >"$STUBDIR/vep" <<'EOF'
#!/usr/bin/env bash
[[ -n "${STUB_CAND_LOG:-}" ]] && printf '%s\n' "$*" >>"$STUB_CAND_LOG"
[[ "${1:-}" == "--version" ]] && { echo "vep 0.0.0-stub"; exit 0; }
out=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    -o) out="$2"; shift 2 ;;
    *) shift ;;
    esac
done
[[ -n "$out" ]] || { echo "stub vep: no -o" >&2; exit 98; }
[[ "${STUB_CAND_RC:-0}" -eq 0 ]] || { echo "stub vep: failing on request" >&2; exit "${STUB_CAND_RC}"; }
name=$(basename "$out")
if [[ -n "${STUB_CAND_EMPTY:-}" ]]; then
    grep '^#' "$STUB_REF_DIR/$name" >"$out"
    exit 0
fi
cp "$STUB_REF_DIR/$name" "$out"
[[ "$name" == "json.jsonl" && "${STUB_CAND_MUTATE:-}" == "corrupt_json" ]] && echo '{not json' >"$out"
if [[ "$name" == "default.txt" ]]; then
    case "${STUB_CAND_MUTATE:-}" in
    change_impact) awk '!done && sub(/IMPACT=MODIFIER/, "IMPACT=LOW") {done=1} 1' "$out" >"$out.tmp" && mv "$out.tmp" "$out" ;;
    drop_last) sed -i.bak '$d' "$out" && rm -f "$out.bak" ;;
    add_row) grep -v '^#' "$out" | tail -1 >>"$out" ;;
    esac
    [[ -n "${STUB_CAND_WARNINGS:-}" ]] && printf '%s\n' "$STUB_CAND_WARNINGS" >"${out}_warnings.txt"
fi
exit 0
EOF
chmod +x "$STUBDIR/docker" "$STUBDIR/vep"
cp "$STUBDIR/vep" "$STUBDIR/vep2"
# The stored-reference checks replace the docker stub with a failing one; this copy restores it.
cp "$STUBDIR/docker" "$TMP/docker-answering"
export STUB_REF_DIR="$REF"

# row <cells...>: one tab-separated line.
row() {
    local IFS=$'\t'
    printf '%s\n' "$*"
}
HEADER=$(row case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs)
CORPUS_FLAGS="--protein --protein_version --mask_header_cache_path --max_sv_size -1"
CORPUS_CACHE="116_GRCh38=tests/golden/116/GRCh38-release116/json_cache"
CORPUS_INPUT="tests/golden/116/GRCh38-release116/variants.vcf"

# run_oracle <table> <extra args...>: the stubs on PATH, every location explicit.
run_oracle() {
    local table="$1"
    shift
    PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$table" --data-dir "$DATA" "$@"
}
dry() { # dry <table> <out file> <extra args...>
    local table="$1" out="$2"
    shift 2
    run_oracle "$table" --out "$TMP/never" --candidates "$STUBDIR/vep" --candidate-pins stub --dry-run "$@" >"$out" 2>"$out.err"
}
# cmd_after <marker> <file>: the command line that follows a DRY-RUN marker line, and
# cmd_second the one after it (the comparator command), trailing blanks dropped.
cmd_after() {
    awk -v m="$1" '$0 == m {getline; sub(/ +$/, ""); print; exit}' "$2"
}
cmd_second() {
    awk -v m="$1" '$0 == m {getline; getline; sub(/ +$/, ""); print; exit}' "$2"
}
count_lines() {
    local n
    n=$(grep -c -- "$1" "$2") || n=0
    echo "$n"
}
# reasons_of <report.json>: the case's reasons on one line, joined by " | ".
reasons_of() {
    python3 -c 'import json,sys; print(" | ".join(json.load(open(sys.argv[1]))["reasons"]))' "$1"
}

echo "DRY-RUN: the engine and comparator commands of every case"
T1="$TMP/t1.tsv"
{
    echo "$HEADER"
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json,warnings "$IMAGE" offline 116.2 - - -
    row c-db fixture "$CORPUS_INPUT" 116_GRCh37 - "--format hgvs" default "$IMAGE" database 116.2 breakend_mate_local_read IMPACT=3,Feature=2 "IMPACT=LOW=>HIGH;Consequence=a=>b"
    row c-pq fixture "$CORPUS_INPUT" 116_GRCh38 - - tab,parquet "$IMAGE" offline 116.2 - - -
    row c-near fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest symbol" default,tab "$IMAGE" offline 116.2 - - -
    row GRCh38-release116 population "$CORPUS_INPUT" 116_GRCh38 - - default "$IMAGE" offline 116.2 - - -
} >"$T1"
dry "$T1" "$TMP/d1.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "dry run over a complete table exits 0" || {
    bad "dry run exited $rc"
    head -5 "$TMP/d1.out.err" | sed 's/^/         /'
}
n=$(count_lines '^DRY-RUN c-four perl ' "$TMP/d1.out")
[[ "$n" -eq 4 ]] && ok "one Ensembl VEP run per output format (4 for default,tab,vcf,json; none for warnings)" || bad "$n Ensembl VEP runs planned for c-four, expected 4"
n=$(grep -A1 '^DRY-RUN c-four perl ' "$TMP/d1.out" | grep '^docker run ' | grep -c -- '--offline --cache --dir_cache /work/cache --cache_version 116 --species homo_sapiens --assembly GRCh38 --format vcf --fork 1 --no_stats --quiet --force_overwrite') || n=0
[[ "$n" -eq 4 ]] && ok "every offline Ensembl VEP run carries --offline --cache --dir_cache --cache_version 116 and the fixed flags" || bad "$n of 4 Ensembl VEP commands carry the offline flag set"
pc=$(cmd_after 'DRY-RUN c-four perl tab' "$TMP/d1.out")
[[ "$pc" == *"-o /work/out/tab.txt"*"--tab" ]] && ok "the tab run writes /work/out/tab.txt with --tab last" || bad "tab run: $pc"
[[ "$pc" == *"--mask_header_cache_path --max_sv_size -1 --tab"* ]] && ok "the case flags are appended verbatim before the format flag" || bad "flags not verbatim: $pc"
[[ "$pc" == *"-v $CORPUS/variants.vcf:/work/input/variants.vcf:ro"* && "$pc" == *"-i /work/input/variants.vcf"* ]] && ok "a repository-relative input resolves and is mounted read-only under /work/input" || bad "input mount: $pc"
[[ "$pc" == *"-v $DATA/caches/perl/vep-cache:/work/cache:ro"* ]] && ok "the Ensembl cache root is mounted read-only at /work/cache" || bad "cache mount: $pc"
[[ "$pc" == *"--user $(id -u):$(id -g) "* ]] && ok "the container runs as the invoking user" || bad "no --user: $pc"
[[ "$pc" == *" -e PERL_HASH_SEED=0 -e PERL_PERTURB_KEYS=0 "* ]] && ok "the container runs Perl under PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0, the seed of the committed references" || bad "no Perl hash seed: $pc"
dc=$(cmd_after 'DRY-RUN c-db perl default' "$TMP/d1.out")
[[ "$dc" == *" -e PERL_HASH_SEED=0 -e PERL_PERTURB_KEYS=0 "* ]] && ok "a database-mode run carries the same Perl hash seed" || bad "no Perl hash seed in database mode: $dc"
[[ "$dc" == *"--database --db_version 116 "* ]] && ok "database mode passes --database --db_version 116" || bad "database mode: $dc"
[[ "$dc" != *"--offline"* && "$dc" != *"--cache "* && "$dc" != *"--dir_cache"* && "$dc" != *"--cache_version"* ]] && ok "database mode passes no --offline, --cache, --dir_cache or --cache_version" || bad "database mode carries a cache flag: $dc"
[[ "$dc" != *":/work/cache:ro"* ]] && ok "database mode mounts no cache" || bad "database mode mounts the cache: $dc"
[[ "$dc" == *"--format hgvs"* && "$dc" != *"--format vcf"* ]] && ok "a --format in the flags replaces the runner's --format vcf" || bad "format handling: $dc"
[[ "$dc" == *"--assembly GRCh37"* ]] && ok "the assembly comes from the cache cell (116_GRCh37)" || bad "assembly: $dc"
cc=$(cmd_after 'DRY-RUN c-four stub default' "$TMP/d1.out")
[[ "$cc" == "$STUBDIR/vep -i $CORPUS/variants.vcf -o $TMP/never/work/cases/c-four/candidates/stub/default.txt --offline --json_cache $CORPUS/json_cache --species homo_sapiens --assembly GRCh38 --format vcf --no_stats --quiet --force_overwrite --protein --protein_version --mask_header_cache_path --max_sv_size -1" ]] &&
    ok "the candidate runs the same flags with --offline --json_cache <the cache cell's =dir>" || bad "candidate default: $cc"
cc=$(cmd_after 'DRY-RUN c-pq stub tab' "$TMP/d1.out")
[[ "$cc" == *"--json_cache $DATA/caches/vep-rs/grch38 "* && "$cc" == *" --tab" ]] && ok "a bare cache cell selects <json-cache-dir>/<assembly> for the candidate" || bad "candidate cache: $cc"
cdb=$(cmd_after 'DRY-RUN c-db stub default' "$TMP/d1.out")
[[ "$cdb" == *"--format hgvs"* && "$cdb" != *"--format vcf"* && "$cdb" == *"--offline --json_cache"* ]] && ok "the candidate of a database case still runs offline on its JSON cache" || bad "candidate in database case: $cdb"
pq=$(cmd_after 'DRY-RUN c-pq stub parquet' "$TMP/d1.out")
[[ "$pq" == *"-o $TMP/never/work/cases/c-pq/candidates/stub/parquet.parquet "*"--output_format parquet" ]] && ok "the parquet run passes --output_format parquet" || bad "parquet run: $pq"
pqc=$(cmd_second 'DRY-RUN c-pq stub parquet' "$TMP/d1.out")
[[ "$pqc" == *"--format parquet --expected $TMP/never/work/cases/c-pq/perl/tab.txt --actual $TMP/never/work/cases/c-pq/candidates/stub/parquet.parquet "* ]] && ok "the parquet comparison reads the Ensembl VEP tab output" || bad "parquet comparison: $pqc"
n=$(count_lines '^DRY-RUN c-pq perl ' "$TMP/d1.out")
[[ "$n" -eq 1 ]] && ok "tab,parquet plans one Ensembl VEP run (tab)" || bad "$n Ensembl VEP runs for tab,parquet"
vc=$(cmd_second 'DRY-RUN c-four stub vcf' "$TMP/d1.out")
jc=$(cmd_second 'DRY-RUN c-four stub json' "$TMP/d1.out")
tc=$(cmd_second 'DRY-RUN c-four stub tab' "$TMP/d1.out")
[[ "$vc" == *"--format vcf --expected $TMP/never/work/cases/c-four/perl/vcf.vcf --actual $TMP/never/work/cases/c-four/candidates/stub/vcf.vcf --report $TMP/never/stub/c-four/vcf.json --input $CORPUS/variants.vcf"* && "$jc" == *"--input "* && "$tc" != *"--input "* ]] &&
    ok "the vcf and json comparisons carry --input; default and tab do not" || bad "comparator --input: vcf=$vc tab=$tc"
[[ "$tc" == *"--manifest $CORPUS/manifest.json"* ]] && ok "a manifest.json beside the input is passed to the comparator" || bad "manifest: $tc"
[[ "$tc" != *"--masks"* && "$tc" != *"--floors"* && "$tc" != *"--field-pairs"* ]] && ok "a - cell passes no --masks, --floors or --field-pairs" || bad "dash cells leaked: $tc"
dcc=$(cmd_second 'DRY-RUN c-db stub default' "$TMP/d1.out")
# printf %q escapes the comma, the semicolon and the > of the printed command.
[[ "$dcc" == *'--masks breakend_mate_local_read --floors IMPACT=3\,Feature=2 --field-pairs IMPACT=LOW=\>HIGH\;Consequence=a=\>b' ]] && ok "masks, floors and field_pairs reach the comparator" || bad "comparator cells: $dcc"
wc_=$(cmd_after 'DRY-RUN c-four stub warnings' "$TMP/d1.out")
[[ "$wc_" == *"--warnings --expected $TMP/never/work/cases/c-four/perl/default.txt_warnings.txt --actual $TMP/never/work/cases/c-four/candidates/stub/default.txt_warnings.txt --report $TMP/never/stub/c-four/warnings.json" ]] &&
    ok "warnings compares the default runs' _warnings.txt files" || bad "warnings command: $wc_"
nc=$(cmd_after 'DRY-RUN c-near perl default' "$TMP/d1.out")
[[ "$nc" == *"-v $DATA/caches/perl/vep-cache:/work/cache "* && "$nc" != *":/work/cache:ro"* ]] && ok "a --nearest case mounts the Ensembl cache writable" || bad "nearest cache mount: $nc"
[[ "$pc" == *":/work/cache:ro"* && "$tc" != *"--nearest-sites"* ]] && ok "a case without --nearest keeps the read-only mount and passes no --nearest-sites" || bad "nearest leaked: $pc / $tc"
ncc=$(cmd_second 'DRY-RUN c-near stub tab' "$TMP/d1.out")
[[ "$ncc" == *"--nearest-sites $CORPUS/json_cache" ]] && ok "a --nearest case passes the candidate's JSON cache as --nearest-sites" || bad "nearest sites: $ncc"
[[ ! -e "$TMP/never" ]] && ok "a dry run writes nothing" || bad "a dry run created the output directory"

dry "$T1" "$TMP/d2.out" --perl-fork 3 --fork 2
pc=$(cmd_after 'DRY-RUN c-four perl default' "$TMP/d2.out")
cc=$(cmd_after 'DRY-RUN c-four stub default' "$TMP/d2.out")
[[ "$pc" == *"--fork 3 "* && "$cc" == *"--fork 2 "* ]] && ok "--perl-fork and --fork set each engine's own --fork" || bad "fork options: perl=$pc cand=$cc"
cc=$(cmd_after 'DRY-RUN c-four stub default' "$TMP/d1.out")
[[ "$cc" != *"--fork"* ]] && ok "without --fork the candidate's own default stands" || bad "candidate got --fork unasked: $cc"

echo
echo "FASTA: a .gz FASTA is decompressed under the work directory for both engines"
T2="$TMP/t2.tsv"
{
    echo "$HEADER"
    row c-fa fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - -
} >"$T2"
dry "$T2" "$TMP/d3.out"
pc=$(cmd_after 'DRY-RUN c-fa perl default' "$TMP/d3.out")
cc=$(cmd_after 'DRY-RUN c-fa stub default' "$TMP/d3.out")
[[ "$pc" == *"-v $TMP/never/work/fasta/1/reference.fa:/work/ref.fa:ro -v $TMP/never/work/fasta/1/reference.fa.fai:/work/ref.fa.fai:ro"* && "$pc" == *"--hgvs --fasta /work/ref.fa"* ]] &&
    ok "Ensembl VEP mounts the decompressed FASTA and its index and gets --fasta" || bad "Ensembl VEP FASTA: $pc"
[[ "$cc" == *"--hgvs --fasta $TMP/never/work/fasta/1/reference.fa" ]] && ok "the candidate gets --fasta with the decompressed FASTA" || bad "candidate FASTA: $cc"

echo
echo "SELECTION: --tier and --only"
dry "$T1" "$TMP/d4.out" --tier population
n=$(count_lines '^DRY-RUN GRCh38-release116 ' "$TMP/d4.out")
m=$(count_lines '^DRY-RUN c-four ' "$TMP/d4.out")
[[ "$n" -gt 0 && "$m" -eq 0 ]] && ok "--tier population selects the population rows only" || bad "--tier population: GRCh38-release116=$n c-four=$m"
n=$(count_lines '^DRY-RUN GRCh38-release116 ' "$TMP/d1.out")
[[ "$n" -eq 0 ]] && ok "the default tier is fixture, so a population row is not planned" || bad "a population row was planned under the default tier"
dry "$T1" "$TMP/d5.out" --only c-db
n=$(count_lines '^DRY-RUN c-db ' "$TMP/d5.out")
m=$(count_lines '^DRY-RUN c-four ' "$TMP/d5.out")
[[ "$n" -gt 0 && "$m" -eq 0 ]] && ok "--only selects the named case" || bad "--only c-db: c-db=$n c-four=$m"
dry "$T1" "$TMP/d6.out" --only c-nope
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "names 'c-nope', which is not a case" "$TMP/d6.out.err" && ok "--only with an unknown case is an error" || bad "--only c-nope exited $rc"
dry "$T2" "$TMP/d7.out" --tier population
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'no case of tier population' "$TMP/d7.out.err" && ok "a tier with no rows is an error, not an empty success" || bad "empty tier exited $rc"

echo
echo "FAIL-CLOSED: table and option errors stop the run before anything is planned"
expect_error() { # expect_error <label> <pattern> <table> [extra args]
    local label="$1" pattern="$2" table="$3"
    shift 3
    dry "$table" "$TMP/e.out" "$@"
    local rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q -- "$pattern" "$TMP/e.out.err" && ! grep -q '^DRY-RUN' "$TMP/e.out"; then
        ok "$label"
    else
        bad "$label (rc=$rc)"
        head -3 "$TMP/e.out.err" | sed 's/^/         /'
    fi
}
T3="$TMP/t3.tsv"
{
    echo "$HEADER"
    row c-miss fixture "$DATA/absent.vcf" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a missing case input exits non-zero and names it" "c-miss: input not found: $DATA/absent.vcf" "$T3"
grep -q 'preflight failed; nothing ran' "$TMP/e.out.err" && ok "the abort is reported as a preflight failure (the other case was fine)" || bad "no preflight-failure line"
{
    echo "$HEADER"
    row c-fa fixture "$CORPUS_INPUT" "$CORPUS_CACHE" "$DATA/absent.fa" - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a missing FASTA is an error" "FASTA not found" "$T3"
{
    echo "$HEADER"
    row c-cache fixture "$CORPUS_INPUT" 116_GRCh38 - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a missing Ensembl cache slice is an error in offline mode" "Ensembl VEP cache slice not found: $DATA/x/homo_sapiens/116_GRCh38" "$T3" --perl-cache-dir "$DATA/x"
expect_error "a candidate JSON cache without transcripts/ is an error" "candidate JSON cache has no transcripts/" "$T3" --json-cache-dir "$DATA/x"
{
    echo "$HEADER"
    row c-fmt fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default,xml "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "an unknown format is an error" "formats names 'xml'" "$T3"
{
    echo "$HEADER"
    row c-w fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab,warnings "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "warnings without default is an error" "warnings without default" "$T3"
{
    echo "$HEADER"
    row c-flag fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--offline --hgvs" default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a flags cell naming a runner-owned flag is an error" "flags names --offline" "$T3"
{
    echo "$HEADER"
    row c-mode fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" online 116.2 - - -
} >"$T3"
expect_error "an unknown perl_mode is an error" "perl_mode 'online'" "$T3"
{
    echo "$HEADER"
    row c-tier nightly "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "an unknown tier is an error" "tier 'nightly'" "$T3"
{
    echo "$HEADER"
    row c-dup fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
    row c-dup fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a duplicate case is an error" "case 'c-dup' appears twice" "$T3"
{
    echo "$HEADER"
    row "../escape" fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a case that is not a directory-safe name is an error (it names the work and report directories)" "case '../escape' is not a directory-safe name" "$T3"
{
    echo "$HEADER"
    row c-empty fixture "$CORPUS_INPUT" "$CORPUS_CACHE" "" - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "an empty cell is an error (write -)" "an empty cell" "$T3"
{
    echo "$HEADER"
    row c-short fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - -
} >"$T3"
expect_error "a short row is an error" "12 columns, expected 13" "$T3"
{
    row case input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs
    row c-x "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a header without the tier column is an error" "header lacks the column 'tier'" "$T3"
{
    echo "$HEADER"
    row c-masks fixture "$DATA/caches/vep-rs/grch38/transcripts" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 some_class - -
} >"$T3"
expect_error "masks without a manifest.json beside the input is an error" "masks names 'some_class', which applies from a manifest.json beside the input, and none is there (the built-in classes are regulatory_feature_rows order_only individual_allele_order)" "$T3"
mkdir -p "$DATA/nomanifest"
cp "$CORPUS/variants.vcf" "$DATA/nomanifest/variants.vcf"
{
    echo "$HEADER"
    row c-builtin fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,regulatory_feature_rows,individual_allele_order - -
    row c-mixed fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,some_class - -
} >"$T3"
dry "$T3" "$TMP/d-builtin.out" --only c-builtin
rc=$?
bc=$(cmd_second 'DRY-RUN c-builtin stub default' "$TMP/d-builtin.out")
[[ "$rc" -eq 0 && "$bc" == *'--masks order_only\,regulatory_feature_rows\,individual_allele_order' && "$bc" != *"--manifest"* ]] && ok "the comparator's built-in classes need no manifest beside the input" || bad "built-in masks: rc=$rc $bc"
expect_error "a manifest class beside a built-in one still needs the manifest" "masks names 'some_class', which applies from a manifest.json" "$T3" --only c-mixed
{
    echo "$HEADER"
    row c-dash fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 -,order_only - -
} >"$T3"
expect_error "a - beside a class in the masks cell is an error (the comparator would read it as an option)" "masks '-,order_only' has an item beginning with -" "$T3"
{
    echo "$HEADER"
    row c-dashfloor fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - -,IMPACT=1 -
    row c-dashpair fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - "-;IMPACT=LOW=>HIGH"
    row c-dashclass fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 -foo - -
    row c-dashafter fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,-foo - -
} >"$T3"
expect_error "a - item in the floors cell is an error" "floors '-,IMPACT=1' has an item beginning with -" "$T3" --only c-dashfloor
expect_error "a - item in the field_pairs cell is an error" "field_pairs '-;IMPACT=LOW=>HIGH' has an item beginning with -" "$T3" --only c-dashpair
expect_error "a masks class beginning with - is an error even beside a manifest" "masks '-foo' has an item beginning with -" "$T3" --only c-dashclass
expect_error "a masks class beginning with - after a separator is an error" "masks 'order_only,-foo' has an item beginning with -" "$T3" --only c-dashafter
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_option_error() { # expect_option_error <label> <pattern> <candidate and tier options>
    local label="$1" pattern="$2"
    shift 2
    run_oracle "$T3" --out "$TMP/never" "$@" --dry-run >/dev/null 2>"$TMP/e.err"
    local rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q -- "$pattern" "$TMP/e.err"; then
        ok "$label"
    else
        bad "$label (rc=$rc)"
    fi
}
expect_option_error "a pin count that differs from the candidate count is an error" 'names 2 binaries and --candidate-pins 1 pins' --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins one
expect_option_error "a duplicate pin is an error" "pin 'one' is given twice" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins one,one
expect_option_error "a candidate that is not an executable file is an error" 'candidate is not an executable file' --candidates "$TMP/no-such-vep" --candidate-pins one
expect_option_error "--candidate-pins is required" '--candidate-pins required' --candidates "$STUBDIR/vep"
expect_option_error "an unknown --tier is refused" 'invalid --tier: smoke' --candidates "$STUBDIR/vep" --candidate-pins one --tier smoke
out=$(bash "$ORACLE" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--cases required' <<<"$out" && ok "no arguments: --cases required, non-zero" || bad "no arguments exited $rc"

echo
echo "REAL PATH: identical outputs are green"
T4="$TMP/t4.tsv"
{
    echo "$HEADER"
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json,warnings "$IMAGE" offline 116.2 - - -
} >"$T4"
OUT="$TMP/out1"
export STUB_DOCKER_LOG="$TMP/docker1.log" STUB_CAND_LOG="$TMP/cand1.log"
run_oracle "$T4" --out "$OUT" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >"$TMP/r1.out" 2>"$TMP/r1.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "identical outputs in every format exit 0" || {
    bad "identical outputs exited $rc"
    tail -5 "$TMP/r1.err" | sed 's/^/         /'
}
[[ -f "$OUT/status.tsv" ]] && ok "status.tsv is written" || bad "no status.tsv"
[[ "$(head -1 "$OUT/status.tsv")" == "$(row case candidate_pin verdict report_digest run_utc)" ]] && ok "status.tsv header is case, candidate_pin, verdict, report_digest, run_utc" || bad "status.tsv header: $(head -1 "$OUT/status.tsv")"
rowline=$(sed -n 2p "$OUT/status.tsv")
[[ "$rowline" == "c-four	abc1234	green	"* ]] && ok "the status row carries the case name, the pin and green" || bad "status row: $rowline"
digest=$(printf '%s' "$rowline" | cut -f4)
expected_digest=$(sha256sum "$OUT/abc1234/c-four/report.json" 2>/dev/null | awk '{print $1}' || shasum -a 256 "$OUT/abc1234/c-four/report.json" | awk '{print $1}')
[[ -n "$digest" && "$digest" == "$expected_digest" ]] && ok "report_digest is the sha256 of report.json" || bad "report_digest $digest != $expected_digest"
run_utc=$(printf '%s' "$rowline" | cut -f5)
[[ "$run_utc" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$ ]] && ok "run_utc is an ISO 8601 UTC timestamp" || bad "run_utc: $run_utc"
[[ "$(wc -l <"$OUT/status.tsv" | tr -d ' ')" -eq 2 ]] && ok "one candidate and one case give one status row" || bad "$(wc -l <"$OUT/status.tsv") status lines"
REPORT="$OUT/abc1234/c-four/report.json"
python3 - "$REPORT" <<'PY' && ok "report.json: verdict green, every requested format green, the image digest, mode, pin and flags recorded" || bad "report.json content"
import json, sys
r = json.load(open(sys.argv[1]))
assert r["verdict"] == "green", r["verdict"]
assert r["reasons"] == []
assert sorted(r["formats"]) == ["default", "json", "tab", "vcf", "warnings"], sorted(r["formats"])
assert all(v["verdict"] == "green" for v in r["formats"].values())
assert r["formats"]["default"]["one_sided"] == {"expected_only": 0, "actual_only": 0}
assert r["perl_image_digest"] == "ensemblorg/ensembl-vep@sha256:stubdigest"
assert r["perl_mode"] == "offline" and r["perl_server"] is None
assert r["candidate_pin"] == "abc1234" and r["candidate_version"] == "vep 0.0.0-stub"
assert r["flags"] == ["--protein", "--protein_version", "--mask_header_cache_path", "--max_sv_size", "-1"]
assert r["formats_requested"] == ["default", "tab", "vcf", "json", "warnings"]
assert r["reference_release"] == "116.2" and r["cache"] == "116_GRCh38"
assert r["manifest"].endswith("GRCh38-release116/manifest.json")
PY
for f in default tab vcf json warnings; do
    [[ -f "$OUT/abc1234/c-four/$f.json" ]] || bad "comparator report $f.json missing"
done
[[ -f "$OUT/abc1234/c-four/warnings.json" ]] && ok "the comparator's per-format reports sit beside report.json" || bad "comparator reports missing"
n=$(count_lines '^run ' "$STUB_DOCKER_LOG")
[[ "$n" -eq 4 ]] && ok "Ensembl VEP ran once per output format (4 docker run calls)" || bad "$n docker run calls, expected 4"
n=$(count_lines ' -o ' "$STUB_CAND_LOG")
[[ "$n" -eq 4 ]] && ok "the candidate ran once per output format" || bad "$n candidate runs, expected 4"
[[ -d "$OUT/work/cases/c-four" && ! -d "$OUT/work/cases/c-four/perl" && ! -d "$OUT/work/cases/c-four/candidates/abc1234" ]] && ok "the engines' outputs are deleted after the comparison" || bad "engine outputs were kept without --keep-outputs"
[[ -f "$OUT/run.json" ]] && grep -q '"abc1234"' "$OUT/run.json" && grep -q 'stubdigest' "$OUT/run.json" && ok "run.json records the pins and the image digest" || bad "run.json missing or incomplete"
grep -q '"verdict": "green"' "$OUT/abc1234/c-four/warnings.json" && ok "no warnings file on either side compares green (an absent file is no warnings)" || bad "warnings verdict without files"
T4F="$TMP/t4f.tsv"
{
    echo "$HEADER"
    row c-floors fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,json "$IMAGE" offline 116.2 - Uploaded_variation=1,IMPACT=1 -
} >"$T4F"
run_oracle "$T4F" --out "$TMP/out-floors" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >/dev/null 2>&1
rc=$?
python3 - "$TMP/out-floors/abc1234/c-floors/report.json" "$rc" <<'PY2' && ok "a floor on a record-level column is met in the default format and listed as not applicable on the json entry of report.json" || bad "json floors: rc=$rc $(reasons_of "$TMP/out-floors/abc1234/c-floors/report.json" 2>/dev/null)"
import json, sys
r = json.load(open(sys.argv[1]))
assert sys.argv[2] == "0", sys.argv[2]
assert r["verdict"] == "green", r["reasons"]
assert r["formats"]["json"]["floors_not_applicable"] == ["Uploaded_variation"], r["formats"]["json"]
assert "floors_not_applicable" not in r["formats"]["default"], r["formats"]["default"]
sub = json.load(open(r["formats"]["json"]["report"]))
assert sub["columns"]["IMPACT"]["floor_fields"] == ["impact"] and sub["columns"]["IMPACT"]["floor_met"] is True, sub["columns"]["IMPACT"]
PY2
{
    echo "$HEADER"
    row c-jsonfloor fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" json "$IMAGE" offline 116.2 - Uploaded_variation=1,IMPACT=1 -
} >"$T4F"
run_oracle "$T4F" --out "$TMP/out-jsonfloor" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(reasons_of "$TMP/out-jsonfloor/abc1234/c-jsonfloor/report.json")" == "json: floor Uploaded_variation is not applicable to json and no text format of the case checked it" ]] &&
    ok "a json-only case with a floor the json format cannot count is red: no format checked the floor" || bad "json-only floor: rc=$rc $(reasons_of "$TMP/out-jsonfloor/abc1234/c-jsonfloor/report.json" 2>/dev/null)"

echo
echo "REAL PATH: every difference the comparator sees makes the case red"
red_run() { # red_run <label> <mutate> <expected reason pattern> [extra oracle args]
    local label="$1" mutate="$2" pattern="$3"
    shift 3
    local out="$TMP/out-$mutate"
    STUB_CAND_MUTATE="$mutate" run_oracle "$T4" --out "$out" --candidates "$STUBDIR/vep" --candidate-pins p1 "$@" >"$TMP/rr.out" 2>"$TMP/rr.err"
    local rc=$?
    local verdict reasons
    verdict=$(sed -n 2p "$out/status.tsv" | cut -f3)
    reasons=$(python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); print("\n".join(r["reasons"])); print("TAB:" + r["formats"]["tab"]["verdict"])' "$out/p1/c-four/report.json" 2>/dev/null)
    if [[ "$rc" -eq 1 && "$verdict" == "red" && "$reasons" == *"$pattern"* && "$reasons" == *"TAB:green"* ]]; then
        ok "$label"
    else
        bad "$label (rc=$rc verdict=$verdict reasons=$reasons)"
    fi
}
red_run "one changed field: red, exit 1, only the default format red" change_impact "default: field IMPACT: 1 differing entries"
red_run "a row missing from the candidate: red (expected_only 1)" drop_last "default: missing in the candidate: 1 entries"
red_run "a row the candidate adds: red (actual_only 1)" add_row "default: unexpected in the candidate: 1 entries"
T5="$TMP/t5.tsv"
{
    echo "$HEADER"
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - IMPACT=100000 -
} >"$T5"
run_oracle "$T5" --out "$TMP/out-floor" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-floor/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: floor IMPACT=100000"* && "$reasons" == *"tab: floor IMPACT=100000"* ]] && ok "an unmet floor is red in every format it applies to" || bad "floor: rc=$rc reasons=$reasons"
STUB_PERL_RC=3 run_oracle "$T4" --out "$TMP/out-perlrc" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>"$TMP/perlrc.err"
rc=$?
reasons=$(reasons_of "$TMP/out-perlrc/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: Ensembl VEP exited 3"* && "$reasons" == *"warnings: Ensembl VEP exited 3"* && "$(sed -n 2p "$TMP/out-perlrc/status.tsv" | cut -f3)" == "red" ]] &&
    ok "an Ensembl VEP run that exits non-zero is red for that format and for warnings" || bad "Ensembl VEP exit 3: rc=$rc reasons=$reasons"
STUB_CAND_RC=4 run_oracle "$T4" --out "$TMP/out-candrc" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-candrc/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: candidate exited 4: stub vep: failing on request"* ]] && ok "a candidate run that exits non-zero is red with its stderr" || bad "candidate exit 4: rc=$rc reasons=$reasons"

echo
echo "WARNINGS: the (line, description) sets of the two default runs"
STUB_PERL_WARNINGS="WARNING: line 12 skipped (21 30165800 del ...): Skipped record" run_oracle "$T4" --out "$TMP/out-w1" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-w1/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == "warnings: warnings: 1 entries only in the expected file" ]] && ok "an Ensembl VEP warning the candidate lacks is red on warnings alone" || bad "warnings one-sided: rc=$rc reasons=$reasons"
W="WARNING: line 12 skipped (21 30165800 del ...): Skipped record"
STUB_PERL_WARNINGS="$W" STUB_CAND_WARNINGS="$W" run_oracle "$T4" --out "$TMP/out-w2" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 ]] && ok "the same warning on both sides is green" || bad "same warnings exited $rc"
STUB_CAND_WARNINGS="WARNING: line 3 skipped: extra" run_oracle "$T4" --out "$TMP/out-w3" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-w3/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == "warnings: warnings: 1 entries only in the actual file" ]] && ok "a candidate warning Ensembl VEP lacks is red" || bad "candidate-only warning: rc=$rc reasons=$reasons"

echo
echo "CANDIDATES: several per run"
run_oracle "$T4" --out "$TMP/out-two" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins aaa,bbb >/dev/null 2>"$TMP/two.err"
rc=$?
[[ "$rc" -eq 0 && -f "$TMP/out-two/aaa/c-four/report.json" && -f "$TMP/out-two/bbb/c-four/report.json" ]] && ok "two candidates produce two report trees" || bad "two candidates: rc=$rc"
[[ "$(cut -f1,2,3 "$TMP/out-two/status.tsv" | tail -n +2 | tr '\n\t' '  ')" == "c-four aaa green c-four bbb green " ]] && ok "status.tsv carries one block per candidate" || bad "status rows: $(tail -n +2 "$TMP/out-two/status.tsv" | tr '\n' ' ')"
: >"$STUB_DOCKER_LOG"
run_oracle "$T4" --out "$TMP/out-two2" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins aaa,bbb >/dev/null 2>&1
n=$(count_lines '^run ' "$STUB_DOCKER_LOG")
[[ "$n" -eq 4 ]] && ok "Ensembl VEP runs once per case and format, not once per candidate" || bad "$n docker run calls for two candidates, expected 4"
STUB_CAND_MUTATE=change_impact run_oracle "$T4" --out "$TMP/out-two3" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins aaa,bbb >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(cut -f3 "$TMP/out-two3/status.tsv" | tail -n +2 | tr '\n' ' ')" == "red red " ]] && ok "a red for any candidate exits 1" || bad "two red candidates: rc=$rc"

echo
echo "DATABASE MODE on the real path"
T6="$TMP/t6.tsv"
{
    echo "$HEADER"
    row c-db fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" database 116.2 - - -
} >"$T6"
: >"$STUB_DOCKER_LOG"
run_oracle "$T6" --out "$TMP/out-db" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
runline=$(grep '^run ' "$STUB_DOCKER_LOG")
[[ "$rc" -eq 0 && "$runline" == *"--database --db_version 116 --species homo_sapiens --assembly GRCh38"* && "$runline" != *"--offline"* && "$runline" != *"/work/cache"* ]] &&
    ok "a database case runs --database --db_version 116 with no --offline and no cache mount" || bad "database run: rc=$rc $runline"
server=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["perl_server"])' "$TMP/out-db/p1/c-db/report.json")
[[ "$server" == "ensembldb.ensembl.org:3306" ]] && ok "the report records the GRCh38 server ensembldb.ensembl.org:3306" || bad "perl_server: $server"

echo
echo "STORED REFERENCE: --reference-dir runs no container"
cat >"$STUBDIR/docker" <<'EOF'
#!/usr/bin/env bash
: >"${STUB_DOCKER_CALLED:?}"
exit 99
EOF
export STUB_DOCKER_CALLED="$TMP/docker-called"
T7="$TMP/t7.tsv"
{
    echo "$HEADER"
    row GRCh38-release116 fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json "$IMAGE" offline 116.2 - - -
} >"$T7"
run_oracle "$T7" --out "$TMP/out-ref" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >"$TMP/ref.out" 2>"$TMP/ref.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "a golden corpus directory is a stored reference: green against the committed expected outputs" || {
    bad "stored reference exited $rc"
    tail -5 "$TMP/ref.err" | sed 's/^/         /'
}
[[ ! -e "$STUB_DOCKER_CALLED" ]] && ok "--reference-dir runs no container" || bad "docker was called with --reference-dir"
grep -q "Ensembl VEP: stored reference $REPO/tests/golden/116/GRCh38-release116/expected" "$TMP/ref.out" && ok "the log names the stored reference" || bad "stored reference not logged"
python3 - "$TMP/out-ref/p1/GRCh38-release116/report.json" <<'PY' && ok "the report records the stored reference, its provenance.json and the image the provenance names" || bad "stored-reference report"
import json, sys
r = json.load(open(sys.argv[1]))
assert r["reference_dir"].endswith("tests/golden/116/GRCh38-release116"), r["reference_dir"]
assert r["reference_provenance"].endswith("tests/golden/116/GRCh38-release116/provenance.json"), r["reference_provenance"]
assert r["perl_image_digest"] == "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da", r["perl_image_digest"]
assert r["reference_provenance_image"] == r["perl_image_digest"]
assert r["verdict"] == "green"
PY
{
    echo "$HEADER"
    row GRCh38-release116 fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,warnings "$IMAGE" offline 116.2 - - -
} >"$T7"
run_oracle "$T7" --out "$TMP/out-ref2" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >/dev/null 2>"$TMP/ref2.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'stored reference lacks default.txt_warnings.txt' "$TMP/ref2.err" && ok "a stored reference without a warnings file fails a warnings case before anything runs" || bad "stored reference without warnings exited $rc"
{
    echo "$HEADER"
    row no-such-corpus fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T7"
run_oracle "$T7" --out "$TMP/out-ref3" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >/dev/null 2>"$TMP/ref3.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'stored reference lacks default.txt\[.gz\]' "$TMP/ref3.err" && ok "a case with no stored reference fails closed rather than running Ensembl VEP" || bad "missing stored reference exited $rc"
[[ ! -e "$STUB_DOCKER_CALLED" ]] && ok "no container ran for the failed stored-reference cases either" || bad "docker was called"
cp "$TMP/docker-answering" "$STUBDIR/docker"

echo
echo "OPTIONS: --keep-outputs, --jobs and a .gz FASTA on the real path"
run_oracle "$T4" --out "$TMP/out-keep" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs --work-dir "$TMP/work-keep" >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && -f "$TMP/work-keep/cases/c-four/perl/default.txt" && -f "$TMP/work-keep/cases/c-four/candidates/p1/json.jsonl" ]] && ok "--keep-outputs keeps both engines' outputs under --work-dir" || bad "--keep-outputs: rc=$rc"
T8="$TMP/t8.tsv"
{
    echo "$HEADER"
    row c-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - -
    row c-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" vcf,json "$IMAGE" offline 116.2 - - -
} >"$T8"
run_oracle "$T8" --out "$TMP/out-jobs" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs 2 >"$TMP/jobs.out" 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(cut -f1,3 "$TMP/out-jobs/status.tsv" | tail -n +2 | tr '\n\t' '  ')" == "c-a green c-b green " ]] && ok "--jobs 2 runs both cases and reports both verdicts in table order" || bad "--jobs 2: rc=$rc $(tr '\n' ' ' <"$TMP/out-jobs/status.tsv")"
grep -q '^\[.*\] \[c-a\] ' "$TMP/jobs.out" && grep -q '^\[.*\] \[c-b\] ' "$TMP/jobs.out" && ok "each case's log lines carry its case name" || bad "case names missing from the log"
# Two --nearest cases on a slice without transcript_gene_tss.txt run Ensembl VEP one after
# the other (the stub records each run's start and end); a slice that has the file takes no lock.
T_NEAR="$TMP/t-near.tsv"
{
    echo "$HEADER"
    row c-near-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest symbol" default "$IMAGE" offline 116.2 - - -
    row c-near-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest gene" default "$IMAGE" offline 116.2 - - -
} >"$T_NEAR"
TSS_FILE="$DATA/caches/perl/vep-cache/homo_sapiens/116_GRCh38/transcript_gene_tss.txt"
if command -v flock >/dev/null; then
    STUB_PERL_SLEEP=1.5 STUB_PERL_TIMES="$TMP/near.times" run_oracle "$T_NEAR" --out "$TMP/out-near" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs 2 --keep-outputs >/dev/null 2>&1 || true
    [[ -e "$TMP/out-near/work/locks/homo_sapiens_116_GRCh38.tss.lock" ]] && ok "two --nearest cases on a slice without the tree file take the slice's lock" || bad "lock file missing under $TMP/out-near/work/locks"
    overlap=$(sort -n "$TMP/near.times" | awk 'NR > 1 && $1 < prev_end {o = 1} {prev_end = $2} END {print o + 0}')
    [[ "$(wc -l <"$TMP/near.times")" -eq 2 && "$overlap" -eq 0 ]] && ok "their Ensembl VEP runs do not overlap" || bad "nearest runs overlapped: $(tr '\n' ' ' <"$TMP/near.times")"
    : >"$TSS_FILE"
    STUB_PERL_SLEEP=0.1 STUB_PERL_TIMES="$TMP/near2.times" run_oracle "$T_NEAR" --out "$TMP/out-near2" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs 2 --keep-outputs >/dev/null 2>&1 || true
    [[ "$(wc -l <"$TMP/near2.times")" -eq 2 && ! -e "$TMP/out-near2/work/locks" ]] && ok "a slice that has the tree file runs its --nearest cases without the lock" || bad "lock taken with the tree file present"
    rm -f "$TSS_FILE"
else
    run_oracle "$T_NEAR" --out "$TMP/out-near" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs 2 >/dev/null 2>"$TMP/near.err"
    rc=$?
    [[ "$rc" -ne 0 && ! -e "$TMP/out-near" ]] && grep -q 'flock is required to run the --nearest cases of .*116_GRCh38 under --jobs 2' "$TMP/near.err" && ok "without flock, two --nearest cases on a slice without the tree file refuse --jobs 2 before anything runs" || bad "without flock at --jobs 2: rc=$rc $(tr '\n' ' ' <"$TMP/near.err")"
    run_oracle "$T_NEAR" --out "$TMP/out-near1" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs >/dev/null 2>&1
    rc=$?
    [[ "$rc" -eq 0 && "$(cut -f3 "$TMP/out-near1/status.tsv" | tail -n +2 | tr '\n' ' ')" == "green green " && ! -e "$TMP/out-near1/work/locks" ]] && ok "without flock, --jobs 1 runs the same cases one after the other and takes no lock" || bad "without flock at --jobs 1: rc=$rc $(tr '\n' ' ' <"$TMP/out-near1/status.tsv")"
fi
run_oracle "$T2" --out "$TMP/out-fa" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs >/dev/null 2>"$TMP/fa.err"
rc=$?
[[ -f "$TMP/out-fa/work/fasta/1/reference.fa" && -f "$TMP/out-fa/work/fasta/1/reference.fa.fai" && "$(head -c 1 "$TMP/out-fa/work/fasta/1/reference.fa")" == ">" ]] &&
    ok "a .gz FASTA is decompressed into <work>/fasta/<n>/ with its .fai beside it" || bad "FASTA decompression (rc=$rc)"

echo
echo "MORE FAIL-CLOSED PATHS: cells, names and option values the first sections did not reach"
{
    echo "$HEADER"
    row c-cache fixture "$CORPUS_INPUT" "GRCh38=tests/golden/116/GRCh38-release116/json_cache" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "a cache cell that is not <version>_<assembly> is an error" "cache 'GRCh38' is not <version>_<assembly>" "$T3"
{
    echo "$HEADER"
    row c-img fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default - offline 116.2 - - -
} >"$T3"
expect_error "a - perl_image is an error" "perl_image is -" "$T3"
{
    echo "$HEADER"
    row c-rel fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline - - - -
} >"$T3"
expect_error "a - reference_release is an error" "reference_release is -" "$T3"
mkdir -p "$DATA/other"
cp "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.gz" "$DATA/other/reference.fa.gz"
cp "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.fai" "$DATA/other/reference.fa.fai"
{
    echo "$HEADER"
    row c-fa1 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - -
    row c-fa2 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache other/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - -
} >"$T3"
dry "$T3" "$TMP/d-fa2.out"
rc=$?
fa1=$(cmd_after 'DRY-RUN c-fa1 stub default' "$TMP/d-fa2.out")
fa2=$(cmd_after 'DRY-RUN c-fa2 stub default' "$TMP/d-fa2.out")
[[ "$rc" -eq 0 && "$fa1" == *"--fasta $TMP/never/work/fasta/1/reference.fa" && "$fa2" == *"--fasta $TMP/never/work/fasta/2/reference.fa" ]] &&
    ok "two .gz FASTAs that decompress to one file name each get a directory of their own" || bad "two FASTAs of one name: rc=$rc fa1=$fa1 fa2=$fa2"
{
    echo "$HEADER"
    row c-fa1 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - -
    row c-fa1b fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - -
} >"$T3"
dry "$T3" "$TMP/d-fa3.out"
fa1=$(cmd_after 'DRY-RUN c-fa1 stub default' "$TMP/d-fa3.out")
fa1b=$(cmd_after 'DRY-RUN c-fa1b stub default' "$TMP/d-fa3.out")
[[ "$fa1" == *"--fasta $TMP/never/work/fasta/1/reference.fa" && "$fa1b" == *"--fasta $TMP/never/work/fasta/1/reference.fa" ]] &&
    ok "two cases on one .gz FASTA share its one decompressed copy" || bad "shared FASTA: fa1=$fa1 fa1b=$fa1b"
{
    echo "$HEADER"
    row c-refseq fixture "$CORPUS_INPUT" 116_GRCh38 - --refseq default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_error "--refseq in the flags selects the homo_sapiens_refseq slice" "Ensembl VEP cache slice not found: $DATA/caches/perl/vep-cache/homo_sapiens_refseq/116_GRCh38" "$T3"
dry "$T1" "$TMP/d8.out" --only GRCh38-release116
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "names 'GRCh38-release116', which is not a case of tier fixture" "$TMP/d8.out.err" && ok "--only naming a case of the other tier is an error, not a silent drop" || bad "--only across tiers exited $rc"
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
expect_option_error "a pin named like the work directory is refused" "candidate pin 'work' names a file or directory the run writes under --out" --candidates "$STUBDIR/vep" --candidate-pins work
expect_option_error "a pin named status.tsv is refused" "candidate pin 'status.tsv' names a file or directory" --candidates "$STUBDIR/vep" --candidate-pins status.tsv
expect_option_error "a pin named run.json is refused" "candidate pin 'run.json' names a file or directory" --candidates "$STUBDIR/vep" --candidate-pins run.json
PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$T3" --data-dir "$DATA" --out "$TMP/never" --work-dir "$TMP/elsewhere" --candidates "$STUBDIR/vep" --candidate-pins work --dry-run >"$TMP/e.out" 2>"$TMP/e.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "the pin work is free when --work-dir lies outside --out" || bad "pin work with an outside work dir exited $rc"
dry "$T1" "$TMP/d9.out" --tier population
grep -q 'WARN: \[feature_oracle\] population tier without --reference-dir' "$TMP/d9.out.err" && ok "a population run without --reference-dir warns that Ensembl VEP runs live" || bad "no population warning"
dry "$T1" "$TMP/d10.out" --tier population --reference-dir "$REPO/tests/golden/116"
rc=$?
if [[ "$rc" -eq 0 ]] && grep -q '^DRY-RUN GRCh38-release116 perl stored' "$TMP/d10.out" && ! grep -q 'population tier without' "$TMP/d10.out.err"; then
    ok "the population warning is silent with --reference-dir (the case was planned against its stored reference)"
else
    bad "population with --reference-dir: rc=$rc"
fi
{
    echo "$HEADER"
    row c-glob fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--filter *" default "$IMAGE" offline 116.2 - - -
} >"$T3"
dry "$T3" "$TMP/d11.out"
gc=$(cmd_after 'DRY-RUN c-glob perl default' "$TMP/d11.out")
[[ "$gc" == *'--filter \*' ]] && ok "a glob character in a flags cell reaches the engines as written" || bad "glob expanded: $gc"

echo
echo "NAMES THAT COLLIDE WITH THE LAYOUT: a pin or case may be named perl or fasta"
STUB_CAND_MUTATE=change_impact run_oracle "$T4" --out "$TMP/out-perlpin" --candidates "$STUBDIR/vep" --candidate-pins perl >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(sed -n 2p "$TMP/out-perlpin/status.tsv" | cut -f3)" == "red" ]] && ok "a candidate pinned perl with a differing output is red (its outputs never reach the Ensembl VEP directory)" || bad "pin perl: rc=$rc $(sed -n 2p "$TMP/out-perlpin/status.tsv")"
{
    echo "$HEADER"
    row fasta fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-fastacase" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && -f "$TMP/out-fastacase/work/fasta/1/reference.fa" ]] && ok "a case named fasta leaves the run's decompressed FASTA in place" || bad "case fasta: rc=$rc"
run_oracle "$T2" --out "$TMP/out-refasta" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs >/dev/null 2>&1
echo "stale" >"$TMP/out-refasta/work/fasta/1/reference.fa"
run_oracle "$T2" --out "$TMP/out-refasta" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(head -c 1 "$TMP/out-refasta/work/fasta/1/reference.fa")" == ">" ]] && ok "a .gz FASTA is decompressed again on every run, never reused from an earlier one" || bad "FASTA reuse: rc=$rc"

echo
echo "EMPTY AND BROKEN OUTPUTS"
T9="$TMP/t9.tsv"
{
    echo "$HEADER"
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - -
} >"$T9"
STUB_PERL_EMPTY=1 STUB_CAND_EMPTY=1 run_oracle "$T9" --out "$TMP/out-empty" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-empty/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: no entry compared: both outputs carry no annotation"* && "$reasons" == *"tab: no entry compared"* ]] && ok "two header-only outputs are red, not a vacuous green" || bad "empty outputs: rc=$rc reasons=$reasons"
STUB_PERL_NOFILE=1 run_oracle "$T9" --out "$TMP/out-nofile" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-nofile/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: Ensembl VEP exited 0 and wrote no default.txt"* ]] && ok "an Ensembl VEP run that exits 0 without its output is red" || bad "no output file: rc=$rc reasons=$reasons"
STUB_CAND_MUTATE=corrupt_json run_oracle "$T4" --out "$TMP/out-corrupt" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-corrupt/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"json: comparator exited 2: ERROR: [compare_vep_formats]"* && "$reasons" != *"default:"* ]] && ok "a comparator error (exit 2) is red for that format with its stderr" || bad "comparator exit 2: rc=$rc reasons=$reasons"

echo
echo "IMAGE PULL, GRCh37 SERVER, INTERRUPTION"
export STUB_PULL_MARKER="$TMP/pulled"
rm -f "$STUB_PULL_MARKER"
run_oracle "$T9" --out "$TMP/out-pull" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/pull.out" 2>&1
rc=$?
[[ "$rc" -eq 0 && -e "$STUB_PULL_MARKER" ]] && grep -q "pulling $IMAGE" "$TMP/pull.out" && grep -q '"perl_image_digest": "ensemblorg/ensembl-vep@sha256:stubdigest"' "$TMP/out-pull/p1/c-four/report.json" &&
    ok "an absent image is pulled once and its digest recorded" || bad "pull path: rc=$rc"
unset STUB_PULL_MARKER
{
    echo "$HEADER"
    row c-db37 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache - - default "$IMAGE" database 116.2 - - -
} >"$T6"
run_oracle "$T6" --out "$TMP/out-db37" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
server=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["perl_server"])' "$TMP/out-db37/p1/c-db37/report.json")
[[ "$server" == "ensembldb.ensembl.org:3337" ]] && ok "a GRCh37 database case records ensembldb.ensembl.org:3337" || bad "GRCh37 perl_server: $server"
# An interrupted run stops its engines at once, under --jobs 1 (the engine is a waited
# job of the main shell) and under --jobs 2 (the engines are children of case
# subshells), and leaves no status table. The stub engine sleeps 30 s; a run that only
# died after the sleep would take 30 s, so the elapsed time is asserted too.
interrupt_run() { # interrupt_run <jobs> [table] [expected engines]: ok lines for one --jobs value
    local jobs="$1" table="${2:-$T8}" expected_engines="${3:-$1}" runner started elapsed alive=0 spid n
    rm -rf "$TMP/out-int" "$TMP/sleepers"
    mkdir -p "$TMP/out-int"
    echo "stale" >"$TMP/out-int/status.tsv"
    started=$SECONDS
    STUB_PERL_SLEEP=30 STUB_PERL_PIDFILE="$TMP/sleepers" PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$table" --data-dir "$DATA" --out "$TMP/out-int" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs "$jobs" >"$TMP/int.out" 2>&1 &
    runner=$!
    for _ in 1 2 3 4 5 6 7 8 9 10; do
        [[ -s "$TMP/sleepers" && "$(wc -l <"$TMP/sleepers" | tr -d ' ')" -ge "$expected_engines" ]] && break
        sleep 1
    done
    # The engines running when the signal lands; a lock holder's death lets a waiting
    # case's engine start for an instant before it too is stopped, so the count precedes the kill.
    n=$(wc -l <"$TMP/sleepers" | tr -d ' ')
    kill -TERM "$runner" 2>/dev/null
    wait "$runner" 2>/dev/null
    elapsed=$((SECONDS - started))
    sleep 1
    while read -r spid; do kill -0 "$spid" 2>/dev/null && alive=$((alive + 1)); done <"$TMP/sleepers"
    [[ "$n" -eq "$expected_engines" && "$alive" -eq 0 && "$elapsed" -lt 15 ]] && ok "SIGTERM to a --jobs $jobs run stops its $n running stub engine(s) at once (${elapsed}s, not after the 30 s sleep)" || bad "--jobs $jobs: $n engine(s) running at the signal (expected $expected_engines), $alive alive after it, exited after ${elapsed}s"
    [[ ! -e "$TMP/out-int/status.tsv" ]] && ok "an interrupted --jobs $jobs run leaves no status.tsv (the stale one was removed at start)" || bad "a status.tsv survived the interrupted --jobs $jobs run"
}
interrupt_run 1
interrupt_run 2
# A first --nearest case holds its slice's lock while it runs, so under --jobs 2 the second
# waits and one engine is running when the signal lands; the lock holder is the engine itself.
if command -v flock >/dev/null; then interrupt_run 2 "$T_NEAR" 1; fi

echo
echo "BOUNDARIES: path spellings, masked-only cases, a broken provenance, a preflight failure"
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
for spelling in "--out $TMP/never/" "--out $TMP/never --work-dir $TMP/never/work/" "--out $TMP/never --work-dir $TMP/never/./work"; do
    # shellcheck disable=SC2086
    PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$T3" --data-dir "$DATA" $spelling --candidates "$STUBDIR/vep" --candidate-pins work --dry-run >/dev/null 2>"$TMP/e.err"
    rc=$?
    [[ "$rc" -ne 0 ]] && grep -q "candidate pin 'work' names a file or directory" "$TMP/e.err" && ok "the pin work is refused under the spelling '$spelling'" || bad "pin work accepted under '$spelling' (rc=$rc)"
done
# a*b is a directory named a*b: were the segment matched against the working directory's
# a1b and a2b, --out would normalize to a1b/a2b and the work directory's parent to
# a1b/a1b/a2b, and the pin work would pass.
mkdir -p "$TMP/globcwd/a1b" "$TMP/globcwd/a2b"
(cd "$TMP/globcwd" && PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$T3" --data-dir "$DATA" --out "$TMP/globcwd/a*b" --work-dir "$TMP/globcwd/a*b/../a*b/work" --candidates "$STUBDIR/vep" --candidate-pins work --dry-run) >/dev/null 2>"$TMP/e.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "candidate pin 'work' names a file or directory" "$TMP/e.err" && ok "a glob character in --out is taken literally by the reserved-pin check (a*b/../a*b/work lies under a*b)" || bad "pin work accepted under a glob --out (rc=$rc)"
# Every entry kept from the breakend records is a documented divergence: with the mask
# on, the comparator compares nothing and skips them all. Skipped entries still count as entries,
# so the empty-output guard stays quiet; the red is the comparator's own, since the stub
# candidate copies the Ensembl VEP output and so never shows the documented divergence.
mkdir -p "$TMP/masked"
{
    grep '^#' "$REF/default.txt"
    grep -v '^#' "$REF/default.txt" | grep -E '^bnd_intergenic_(pair_a|with_genic_mate)' | grep 'ENST00000334068'
} >"$TMP/masked/default.txt"
cp "$CORPUS/manifest.json" "$TMP/masked/manifest.json"
grep -E 'bnd_intergenic_(pair_a|with_genic_mate)' "$CORPUS/variants.vcf" >"$TMP/masked/variants.vcf"
n=$(grep -vc '^#' "$TMP/masked/default.txt")
{
    echo "$HEADER"
    row c-masked fixture "$TMP/masked/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 breakend_mate_local_read - -
} >"$T3"
STUB_REF_DIR="$TMP/masked" run_oracle "$T3" --out "$TMP/out-masked" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/masked.out" 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-masked/p1/c-masked/report.json")
skipped=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["counts"]["documented_divergences_skipped"])' "$TMP/out-masked/p1/c-masked/default.json")
[[ "$n" -ge 1 && "$skipped" -ge 1 && "$rc" -eq 1 && "$reasons" == *"documented divergences not observed"* && "$reasons" != *"no entry compared"* ]] &&
    ok "a case whose every entry is a documented divergence is judged by the comparator, not read as two empty outputs ($skipped skipped)" || bad "masked-only case: rc=$rc skipped=$skipped reasons=$reasons"
mkdir -p "$TMP/refbad/c-badprov/expected"
cp "$REF/default.txt" "$TMP/refbad/c-badprov/expected/default.txt"
echo '{not json' >"$TMP/refbad/c-badprov/provenance.json"
{
    echo "$HEADER"
    row c-badprov fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-badprov" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$TMP/refbad" >"$TMP/badprov.out" 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-badprov/p1/c-badprov/report.json" 2>/dev/null)
[[ "$rc" -eq 1 && "$reasons" == *"stored reference provenance.json unreadable"* && "$(sed -n 2p "$TMP/out-badprov/status.tsv" | cut -f3)" == "red" ]] && ! grep -q Traceback "$TMP/badprov.out" &&
    ok "an unreadable provenance.json beside a stored reference is a red reason, not a traceback" || bad "bad provenance: rc=$rc reasons=$reasons"
mkdir -p "$TMP/out-pre"
echo "stale" >"$TMP/out-pre/status.tsv"
{
    echo "$HEADER"
    row c-miss fixture "$DATA/absent.vcf" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-pre" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
[[ "$rc" -ne 0 && ! -e "$TMP/out-pre/status.tsv" ]] && ok "a run that fails preflight removes the previous run's status.tsv" || bad "preflight failure left status.tsv (rc=$rc)"
mkdir -p "$TMP/out-dry"
echo "stale" >"$TMP/out-dry/status.tsv"
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-dry" --candidates "$STUBDIR/vep" --candidate-pins p1 --dry-run >/dev/null 2>&1
[[ -e "$TMP/out-dry/status.tsv" ]] && ok "a dry run leaves an existing status.tsv alone" || bad "a dry run removed status.tsv"

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
