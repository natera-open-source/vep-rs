#!/usr/bin/env bash
# Tests for feature_oracle.sh: the case-table contract, the dry-run command shapes, every
# fail-closed path, and the verdicts on the real path through a stub `docker` and a stub
# `vep`, both of which answer with the committed Ensembl VEP outputs of
# tests/golden/116/GRCh38-release116 so the real comparator decides.
#
# Bash and python3. No Docker, no VEP, no cache.
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
# STUB_CAND_RC (printing STUB_PERL_STDERR / STUB_CAND_STDERR to stderr), sleep in place of a
# run under STUB_PERL_SLEEP / STUB_CAND_SLEEP (recording their pids in STUB_PERL_PIDFILE /
# STUB_CAND_PIDFILE), and append their argv to STUB_DOCKER_LOG / STUB_CAND_LOG.
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
if [[ "${STUB_PERL_RC:-0}" -ne 0 ]]; then
    [[ -z "${STUB_PERL_STDERR:-}" ]] || echo "$STUB_PERL_STDERR" >&2
    exit "${STUB_PERL_RC}"
fi
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
if [[ -n "${STUB_CAND_SLEEP:-}" ]]; then
    [[ -z "${STUB_CAND_PIDFILE:-}" ]] || echo "$$" >>"$STUB_CAND_PIDFILE"
    exec sleep "$STUB_CAND_SLEEP"
fi
# STUB_CAND_KILL_JOB kills the shell that runs the engine, a format job under --format-jobs,
# after writing STUB_CAND_JOB_STDERR, when set, to that job's own stderr file beside the output.
if [[ -n "${STUB_CAND_KILL_JOB:-}" ]]; then
    fmt=$(basename "$out")
    [[ -z "${STUB_CAND_JOB_STDERR:-}" ]] || printf '%s\n' "$STUB_CAND_JOB_STDERR" >>"$(dirname "$out")/${fmt%%.*}.job.err"
    kill -KILL "$PPID"
    exit 0
fi
[[ "${STUB_CAND_RC:-0}" -eq 0 ]] || { echo "${STUB_CAND_STDERR:-stub vep: failing on request}" >&2; exit "${STUB_CAND_RC}"; }
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
    add_splice_region) awk -F'\t' -v OFS='\t' '!done && !/^#/ && $7 == "intron_variant" {$7 = "splice_region_variant,intron_variant"; done=1} 1' "$out" >"$out.tmp" && mv "$out.tmp" "$out" ;;
    own_span_protein) sed -i.bak 's/DownstreamProtein=VK;/DownstreamProtein=PKVK;/' "$out" && rm -f "$out.bak" ;;
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
HEADER=$(row case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs species unordered reference refusal)
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
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json,warnings "$IMAGE" offline 116.2 - - - - - - -
    row c-db fixture "$CORPUS_INPUT" 116_GRCh37 - "--format hgvs" default "$IMAGE" database 116.2 breakend_mate_local_read IMPACT=3,Feature=2 "IMPACT=LOW=>HIGH;Consequence=a=>b" - - - -
    row c-pq fixture "$CORPUS_INPUT" 116_GRCh38 - - tab,parquet "$IMAGE" offline 116.2 - - - - - - -
    row c-near fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest symbol" default,tab "$IMAGE" offline 116.2 - - - - - - -
    row c-skip fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 reference_skipped_record=4 - - - - - -
    row GRCh38-release116 population "$CORPUS_INPUT" 116_GRCh38 - - default "$IMAGE" offline 116.2 - - - - - - -
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
[[ "$vc" == *"--format vcf --expected $TMP/never/work/cases/c-four/perl/vcf.vcf --actual $TMP/never/work/cases/c-four/candidates/stub/vcf.vcf --report $TMP/never/stub/c-four/vcf.json --input $CORPUS/variants.vcf"* && "$jc" == *"--input "* ]] &&
    ok "the vcf and json comparisons carry --input" || bad "comparator --input: vcf=$vc json=$jc"
[[ "$tc" == *"--manifest $CORPUS/manifest.json"* ]] && ok "a manifest.json beside the input is passed to the comparator" || bad "manifest: $tc"
# The manifest's flag set over the input is selected by the input's name, so the text
# formats name it too whenever a manifest is read.
dfc=$(cmd_second 'DRY-RUN c-four stub default' "$TMP/d1.out")
[[ "$tc" == *"--format tab "*" --input $CORPUS/variants.vcf --manifest $CORPUS/manifest.json"* && "$dfc" == *"--format default "*" --input $CORPUS/variants.vcf --manifest $CORPUS/manifest.json"* ]] &&
    ok "a default or tab comparison reading a manifest names the input before the manifest" || bad "text-format --input beside a manifest: tab=$tc default=$dfc"
[[ "$tc" != *"--masks"* && "$tc" != *"--floors"* && "$tc" != *"--field-pairs"* ]] && ok "a - cell passes no --masks, --floors or --field-pairs" || bad "dash cells leaked: $tc"
dcc=$(cmd_second 'DRY-RUN c-db stub default' "$TMP/d1.out")
# printf %q escapes the comma, the semicolon and the > of the printed command.
[[ "$dcc" == *'--masks breakend_mate_local_read --floors IMPACT=3\,Feature=2 --field-pairs IMPACT=LOW=\>HIGH\;Consequence=a=\>b --json-cache '"$DATA/caches/vep-rs/grch37" ]] && ok "masks, floors and field_pairs reach the comparator; a cache-reading class brings --json-cache" || bad "comparator cells: $dcc"
[[ "$dcc" == *" --input $CORPUS/variants.vcf --manifest "* ]] && ok "a default comparison naming no input-reading class still names the input for its manifest" || bad "input beside a manifest: $dcc"
wc_=$(cmd_after 'DRY-RUN c-four stub warnings' "$TMP/d1.out")
[[ "$wc_" == *"--warnings --expected $TMP/never/work/cases/c-four/perl/default.txt_warnings.txt --actual $TMP/never/work/cases/c-four/candidates/stub/default.txt_warnings.txt --report $TMP/never/stub/c-four/warnings.json" ]] &&
    ok "warnings compares the default runs' _warnings.txt files" || bad "warnings command: $wc_"
nc=$(cmd_after 'DRY-RUN c-near perl default' "$TMP/d1.out")
[[ "$nc" == *"-v $DATA/caches/perl/vep-cache:/work/cache "* && "$nc" != *":/work/cache:ro"* ]] && ok "a --nearest case mounts the Ensembl cache writable" || bad "nearest cache mount: $nc"
[[ "$pc" == *":/work/cache:ro"* && "$tc" != *"--nearest-sites"* ]] && ok "a case without --nearest keeps the read-only mount and passes no --nearest-sites" || bad "nearest leaked: $pc / $tc"
ncc=$(cmd_second 'DRY-RUN c-near stub tab' "$TMP/d1.out")
[[ "$ncc" == *"--nearest-sites $CORPUS/json_cache" ]] && ok "a --nearest case passes the candidate's JSON cache as --nearest-sites" || bad "nearest sites: $ncc"
skc=$(cmd_second 'DRY-RUN c-skip stub default' "$TMP/d1.out")
[[ "$skc" == *"--format default "*" --input $CORPUS/variants.vcf "*"--masks reference_skipped_record" && "$skc" != *"--json-cache"* ]] &&
    ok "a default comparison naming reference_skipped_record carries --input, the count stripped, and no --json-cache" || bad "skipped-record input: $skc"
[[ ! -e "$TMP/never" ]] && ok "a dry run writes nothing" || bad "a dry run created the output directory"
# vep reads its options with Getopt::Long's defaults (vep:27, Config.pm:62 at release/116.2):
# a single dash, any unambiguous abbreviation and any letter case name the option, so a case
# spelled -nearest, --near or --NEAREST is a --nearest case, its flag word passed as written.
T_SP="$TMP/t-spellings.tsv"
{
    echo "$HEADER"
    row c-near-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-nearest symbol" default,tab "$IMAGE" offline 116.2 - - - - - - -
    row c-near-abbr fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--near symbol" default,tab "$IMAGE" offline 116.2 - - - - - - -
    row c-near-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--NEAREST symbol" default,tab "$IMAGE" offline 116.2 - - - - - - -
} >"$T_SP"
dry "$T_SP" "$TMP/dsp.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "a table spelling --nearest as -nearest, --near and --NEAREST plans" || bad "spellings dry run exited $rc: $(head -3 "$TMP/dsp.out.err" | tr '\n' ' ')"
for case_spelling in "c-near-dash:-nearest" "c-near-abbr:--near" "c-near-upper:--NEAREST"; do
    c="${case_spelling%%:*}"
    sp="${case_spelling#*:}"
    nc=$(cmd_after "DRY-RUN $c perl default" "$TMP/dsp.out")
    ncc=$(cmd_second "DRY-RUN $c stub tab" "$TMP/dsp.out")
    [[ "$nc" == *"-v $DATA/caches/perl/vep-cache:/work/cache "* && "$nc" != *":/work/cache:ro"* && "$nc" == *"--force_overwrite $sp symbol" && "$ncc" == *"--nearest-sites $CORPUS/json_cache" ]] &&
        ok "the $sp spelling of --nearest mounts the cache writable, is passed as written and settles NEAREST by start site" || bad "$sp: $nc / $ncc"
done

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
    row c-fa fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - -
} >"$T2"
dry "$T2" "$TMP/d3.out"
pc=$(cmd_after 'DRY-RUN c-fa perl default' "$TMP/d3.out")
cc=$(cmd_after 'DRY-RUN c-fa stub default' "$TMP/d3.out")
[[ "$pc" == *"-v $TMP/never/work/fasta/1/reference.fa:/work/ref.fa:ro -v $TMP/never/work/fasta/1/reference.fa.fai:/work/ref.fa.fai:ro"* && "$pc" == *"--hgvs --fasta /work/ref.fa"* ]] &&
    ok "Ensembl VEP mounts the decompressed FASTA and its index and gets --fasta" || bad "Ensembl VEP FASTA: $pc"
[[ "$cc" == *"--hgvs --fasta $TMP/never/work/fasta/1/reference.fa" ]] && ok "the candidate gets --fasta with the decompressed FASTA" || bad "candidate FASTA: $cc"

echo
echo "SPECIES: the species cell names both engines' --species and the Ensembl cache slice"
# feature_oracle_dry_run_names_the_case_species_for_both_engines: a yeast row plans both
# commands under --species saccharomyces_cerevisiae with its slice under the species directory;
# a - cell is homo_sapiens.
YEAST_SLICE="$DATA/caches/perl/vep-cache/saccharomyces_cerevisiae/116_R64-1-1"
mkdir -p "$YEAST_SLICE" "$DATA/caches/vep-rs/yeast/transcripts"
TS="$TMP/ts.tsv"
{
    echo "$HEADER"
    row c-yeast fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - - default "$IMAGE" offline 116.2 - - - saccharomyces_cerevisiae - - -
    row c-human fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
    row c-yeast-db fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - --genomes default "$IMAGE" database 116.2 - - - saccharomyces_cerevisiae - - -
} >"$TS"
dry "$TS" "$TMP/ds.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "feature_oracle_dry_run_names_the_case_species_for_both_engines: a table with a yeast row plans" || {
    bad "species dry run exited $rc"
    head -5 "$TMP/ds.out.err" | sed 's/^/         /'
}
pc=$(cmd_after 'DRY-RUN c-yeast perl default' "$TMP/ds.out")
cc=$(cmd_after 'DRY-RUN c-yeast stub default' "$TMP/ds.out")
[[ "$pc" == *"--cache_version 116 --species saccharomyces_cerevisiae --assembly R64-1-1 --format vcf"* ]] && ok "Ensembl VEP runs under the row's --species and the cache cell's assembly" || bad "Ensembl VEP species: $pc"
[[ "$cc" == *"--json_cache $DATA/caches/vep-rs/yeast --species saccharomyces_cerevisiae --assembly R64-1-1 --format vcf"* ]] && ok "the candidate runs under the same --species and --assembly" || bad "candidate species: $cc"
grep -q "^DRY-RUN c-yeast cache species=saccharomyces_cerevisiae perl=$YEAST_SLICE candidate=$DATA/caches/vep-rs/yeast$" "$TMP/ds.out" && ok "the dry run names the slice under the species directory, <perl-cache-dir>/saccharomyces_cerevisiae/116_R64-1-1" || bad "species slice line: $(grep '^DRY-RUN c-yeast cache' "$TMP/ds.out")"
hc=$(cmd_after 'DRY-RUN c-human perl default' "$TMP/ds.out")
[[ "$hc" == *"--species homo_sapiens --assembly GRCh38"* ]] && grep -q "^DRY-RUN c-human cache species=homo_sapiens perl=$DATA/caches/perl/vep-cache/homo_sapiens/116_GRCh38 " "$TMP/ds.out" && ok "a - species cell is homo_sapiens for both the command and the slice" || bad "default species: $hc"
dc=$(cmd_after 'DRY-RUN c-yeast-db perl default' "$TMP/ds.out")
[[ "$dc" == *"--database --db_version 116 --species saccharomyces_cerevisiae --assembly R64-1-1 "*"--genomes" ]] && ok "a database case passes the species beside --db_version and keeps --genomes in its flags" || bad "database species: $dc"
dcc=$(cmd_after 'DRY-RUN c-yeast-db stub default' "$TMP/ds.out")
[[ "$dcc" == *"--species saccharomyces_cerevisiae --assembly R64-1-1 "* && "$dcc" != *"--genomes"* ]] && ok "the candidate of a --genomes case runs without --genomes, which selects Ensembl VEP's server alone" || bad "candidate under --genomes: $dcc"
rmdir "$YEAST_SLICE"
dry "$TS" "$TMP/ds2.out" --only c-yeast
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "c-yeast: Ensembl VEP cache slice not found: $YEAST_SLICE" "$TMP/ds2.out.err" && ok "a missing slice is named under the species directory" || bad "missing species slice exited $rc"
mkdir -p "$YEAST_SLICE"
{
    echo "$HEADER"
    row c-yeast-db fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - --genomes default "$IMAGE" database 116.2 - - - saccharomyces_cerevisiae - - -
    row c-yeast fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - - default "$IMAGE" offline 116.2 - - - saccharomyces_cerevisiae - - -
} >"$TS"
run_oracle "$TS" --out "$TMP/out-species" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>"$TMP/species.err"
rc=$?
[[ "$rc" -eq 0 ]] && python3 - "$TMP/out-species/p1" <<'PY' && ok "the real path records the species and the Ensembl Genomes server of a --genomes database case in report.json" || bad "species real path: rc=$rc $(tail -2 "$TMP/species.err" | tr '\n' ' ')"
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
db = json.loads((out / "c-yeast-db" / "report.json").read_text())
assert db["species"] == "saccharomyces_cerevisiae" and db["cache"] == "116_R64-1-1", (db["species"], db["cache"])
assert db["perl_server"] == "mysql-eg-publicsql.ebi.ac.uk:4157", db["perl_server"]
offline = json.loads((out / "c-yeast" / "report.json").read_text())
assert offline["species"] == "saccharomyces_cerevisiae" and offline["perl_server"] is None, (offline["species"], offline["perl_server"])
PY

echo
echo "PLACEHOLDERS: {corpus} and {sources} name the input's directory and --sources-dir"
# corpus_and_sources_placeholders_mount_read_only: Ensembl VEP mounts both directories
# read-only and reads the words rewritten to the mounts; the candidate reads the host paths.
mkdir -p "$DATA/sources/shared"
TP="$TMP/tp.tsv"
PH_FLAGS="--custom file={corpus}/custom/sites.vcf.gz,short_name=CE,format=vcf,type=exact --custom {sources}/shared/clinvar.vcf.gz,CV,vcf,exact"
{
    echo "$HEADER"
    row c-both fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$PH_FLAGS" default "$IMAGE" offline 116.2 - - - - - - -
    row c-corpus fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--custom {corpus}/custom/regions.bed.gz,RE,bed,overlap" default "$IMAGE" offline 116.2 - - - - - - -
} >"$TP"
dry "$TP" "$TMP/dp.out" --sources-dir "$DATA/sources"
rc=$?
[[ "$rc" -eq 0 ]] && ok "corpus_and_sources_placeholders_mount_read_only: a table whose flags carry {corpus} and {sources} plans under --sources-dir" || {
    bad "placeholder dry run exited $rc"
    head -5 "$TMP/dp.out.err" | sed 's/^/         /'
}
pc=$(cmd_after 'DRY-RUN c-both perl default' "$TMP/dp.out")
cc=$(cmd_after 'DRY-RUN c-both stub default' "$TMP/dp.out")
[[ "$pc" == *"-v $CORPUS:/work/corpus:ro -v $DATA/sources:/work/sources:ro $IMAGE vep "* ]] && ok "Ensembl VEP mounts the input's directory at /work/corpus and --sources-dir at /work/sources, both read-only" || bad "placeholder mounts: $pc"
perl_words=$(printf '%q ' --custom "file=/work/corpus/custom/sites.vcf.gz,short_name=CE,format=vcf,type=exact" --custom "/work/sources/shared/clinvar.vcf.gz,CV,vcf,exact")
[[ "$pc" == *"--force_overwrite ${perl_words% }" ]] && ok "Ensembl VEP reads the words with {corpus} and {sources} rewritten to the mounts, keyed or positional" || bad "Ensembl VEP words: $pc"
cand_words=$(printf '%q ' --custom "file=$CORPUS/custom/sites.vcf.gz,short_name=CE,format=vcf,type=exact" --custom "$DATA/sources/shared/clinvar.vcf.gz,CV,vcf,exact")
[[ "$cc" == *"--force_overwrite ${cand_words% }" ]] && ok "the candidate reads the words with the host directories in place of the placeholders" || bad "candidate words: $cc"
grep '^docker run \|^/' "$TMP/dp.out" | grep -q 'corpus}\|sources}' && bad "a placeholder reached an engine command unrewritten" || ok "no engine command carries a placeholder"
pc=$(cmd_after 'DRY-RUN c-corpus perl default' "$TMP/dp.out")
[[ "$pc" == *":/work/corpus:ro"* && "$pc" != *"/work/sources"* ]] && ok "a case carrying {corpus} alone mounts the corpus and no sources directory" || bad "corpus-only mounts: $pc"
pc=$(cmd_after 'DRY-RUN c-four perl default' "$TMP/d1.out")
[[ "$pc" != *"/work/corpus"* && "$pc" != *"/work/sources"* ]] && ok "a case without a placeholder mounts neither directory" || bad "a placeholder mount leaked into a plain case: $pc"
dry "$TP" "$TMP/dp2.out"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'c-both: flags carry {sources} and no --sources-dir' "$TMP/dp2.out.err" && ! grep -q '^DRY-RUN' "$TMP/dp2.out" &&
    ok "a {sources} word without --sources-dir fails preflight by case, and nothing is planned" || bad "{sources} without --sources-dir exited $rc"
dry "$TP" "$TMP/dp3.out" --sources-dir "$TMP/no-such-sources"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'sources-dir not a directory' "$TMP/dp3.out.err" && ok "a --sources-dir that is not a directory is refused by name" || bad "bad --sources-dir exited $rc"
run_oracle "$TP" --out "$TMP/out-ph" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-both --sources-dir "$DATA/sources" >/dev/null 2>"$TMP/ph.err"
rc=$?
[[ "$rc" -eq 0 ]] && python3 - "$TMP/out-ph" "$DATA/sources" <<'PY' && ok "the real path records --sources-dir in report.json and run.json and the case is green" || bad "placeholder real path: rc=$rc"
import json, sys
from pathlib import Path
out, sources = Path(sys.argv[1]), sys.argv[2]
report = json.loads((out / "p1" / "c-both" / "report.json").read_text())
assert report["verdict"] == "green", report["reasons"]
assert report["sources_dir"] == sources, report["sources_dir"]
assert report["flags"][1].startswith("file={corpus}/"), report["flags"]
assert json.loads((out / "run.json").read_text())["sources_dir"] == sources
PY

echo
echo "DIR: a flags cell names --dir through {corpus} or {sources}"
# dir_flag_names_its_directory_through_a_placeholder: both engines get --dir as they get
# every placeholder word, Ensembl VEP at the read-only mount and the candidate at the host
# path, and each command keeps the runner's own cache word.
mkdir -p "$DATA/sources/gwas/vepdir" "$DATA/sources/ini"
: >"$DATA/sources/ini/vep.ini"
TD="$TMP/td.tsv"
{
    echo "$HEADER"
    row c-dir fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {sources}/gwas/vepdir --plugin StubPlugin" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-eq fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir={corpus} --plugin StubPlugin" default "$IMAGE" offline 116.2 - - - - - - -
    row c-word-free fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-HGVS +Symbol" default "$IMAGE" offline 116.2 - - - - - - -
} >"$TD"
dry "$TD" "$TMP/dd.out" --sources-dir "$DATA/sources"
rc=$?
[[ "$rc" -eq 0 ]] && ok "dir_flag_names_its_directory_through_a_placeholder: a table whose flags name --dir {sources}/<directory> and --dir={corpus} plans" || {
    bad "--dir dry run exited $rc"
    head -5 "$TMP/dd.out.err" | sed 's/^/         /'
}
pc=$(cmd_after 'DRY-RUN c-dir perl default' "$TMP/dd.out")
cc=$(cmd_after 'DRY-RUN c-dir stub default' "$TMP/dd.out")
[[ "$pc" == *"-v $DATA/sources:/work/sources:ro $IMAGE vep "*"--force_overwrite --dir /work/sources/gwas/vepdir --plugin StubPlugin" ]] && ok "Ensembl VEP reads --dir at the read-only /work/sources mount" || bad "Ensembl VEP --dir: $pc"
[[ "$pc" == *" --offline --cache --dir_cache /work/cache --cache_version 116 "*" --dir /work/sources/"* ]] && ok "the Ensembl VEP command of a --dir case keeps the runner's --dir_cache /work/cache" || bad "Ensembl VEP cache words beside --dir: $pc"
[[ "$cc" == *" --offline --json_cache $CORPUS/json_cache "*"--force_overwrite --dir $DATA/sources/gwas/vepdir --plugin StubPlugin" ]] && ok "the candidate reads --dir at the host path beside the runner's --json_cache" || bad "candidate --dir: $cc"
pc=$(cmd_after 'DRY-RUN c-dir-eq perl default' "$TMP/dd.out")
cc=$(cmd_after 'DRY-RUN c-dir-eq stub default' "$TMP/dd.out")
[[ "$pc" == *"-v $CORPUS:/work/corpus:ro $IMAGE vep "*" --dir_cache /work/cache "*"--force_overwrite --dir=/work/corpus --plugin StubPlugin" && "$cc" == *" --json_cache $CORPUS/json_cache "*"--force_overwrite --dir=$CORPUS --plugin StubPlugin" ]] &&
    ok "the --dir=<directory> spelling reaches both engines the same two ways, {corpus} being the input's directory" || bad "--dir= spelling: perl=$pc cand=$cc"
STUB_DOCKER_LOG="$TMP/dir-docker.log" STUB_CAND_LOG="$TMP/dir-cand.log" run_oracle "$TD" --out "$TMP/out-dir" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-dir --sources-dir "$DATA/sources" >/dev/null 2>"$TMP/dir.err"
rc=$?
[[ "$rc" -eq 0 ]] && grep -q -- '--force_overwrite --dir /work/sources/gwas/vepdir --plugin StubPlugin$' "$TMP/dir-docker.log" && grep -q -- "--force_overwrite --dir $DATA/sources/gwas/vepdir --plugin StubPlugin\$" "$TMP/dir-cand.log" &&
    python3 - "$TMP/out-dir/p1/c-dir/report.json" <<'PY' && ok "the real path runs each engine with its own --dir directory, and the case is green with the cell's words in report.json" || bad "--dir real path: rc=$rc $(tail -2 "$TMP/dir.err" | tr '\n' ' ')"
import json, sys
report = json.load(open(sys.argv[1]))
assert report["verdict"] == "green", report["reasons"]
assert report["flags"] == ["--dir", "{sources}/gwas/vepdir", "--plugin", "StubPlugin"], report["flags"]
PY
dry "$TD" "$TMP/dd0.out"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'c-dir: flags carry {sources} and no --sources-dir' "$TMP/dd0.out.err" && ! grep -q 'flags give --dir\|flags names --dir' "$TMP/dd0.out.err" &&
    ok "a --dir {sources} word without --sources-dir is refused for the missing option alone" || bad "--dir {sources} without --sources-dir exited $rc: $(head -3 "$TMP/dd0.out.err" | tr '\n' ' ')"
# dir_flag_outside_the_placeholders_or_over_a_vep_ini_is_refused: Ensembl VEP's container
# sees the two mounts alone, and a vep.ini under --dir is read as options by both engines
# (Config.pm:704-705 at release/116.2), so each such cell is refused by case and by name;
# the value is read after `=` or from the next word, in any spelling of the flag.
# runner_word_is_refused_in_every_spelling: vep reads a single dash or a plus sign and any
# letter case as the option (vep:27, Config.pm:62), so a word the runner sets is refused
# in those spellings too, and a word it does not set passes in them.
{
    echo "$HEADER"
    row c-dir-abs fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir $DATA/sources/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-rel fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir=gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-up fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {corpus}/../GRCh37-hgvs" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-bare fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--plugin StubPlugin --dir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-dir /work/sources/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-ini fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {sources}/ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-ini-eq fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--DIR={sources}/ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-cache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir_cache {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-config fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--config {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-dir-two fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {corpus}/{sources}/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-word-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-config {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-word-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--CONFIG={sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-word-dircache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-dir_cache {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
} >"$TD"
dry "$TD" "$TMP/dd2.out" --sources-dir "$DATA/sources"
rc=$?
not_a_mount="which is not {corpus}, {sources} or a path under one without a .. segment"
[[ "$rc" -ne 0 ]] && ! grep -q '^DRY-RUN' "$TMP/dd2.out" && grep -qF "c-dir-abs: flags give --dir '$DATA/sources/gwas/vepdir', $not_a_mount" "$TMP/dd2.out.err" &&
    ok "an absolute --dir directory is refused by case and by name, and nothing is planned" || bad "absolute --dir exited $rc: $(grep -F 'c-dir-abs:' "$TMP/dd2.out.err")"
grep -qF "c-dir-rel: flags give --dir 'gwas/vepdir', $not_a_mount" "$TMP/dd2.out.err" && ok "a relative --dir=<directory> is refused by case and by name" || bad "relative --dir: $(grep -F 'c-dir-rel:' "$TMP/dd2.out.err")"
grep -qF "c-dir-up: flags give --dir '{corpus}/../GRCh37-hgvs', $not_a_mount" "$TMP/dd2.out.err" && ok "a --dir path that leaves its placeholder through a .. segment is refused" || bad "--dir with ..: $(grep -F 'c-dir-up:' "$TMP/dd2.out.err")"
grep -qF "c-dir-bare: flags give --dir '', $not_a_mount" "$TMP/dd2.out.err" && ok "a --dir with no directory after it is refused" || bad "bare --dir: $(grep -F 'c-dir-bare:' "$TMP/dd2.out.err")"
grep -qF "c-dir-dash: flags give --dir '/work/sources/gwas/vepdir', $not_a_mount" "$TMP/dd2.out.err" && ok "the -dir spelling is held to the same rule: a container path is not a placeholder" || bad "-dir spelling: $(grep -F 'c-dir-dash:' "$TMP/dd2.out.err")"
grep -qF "c-dir-ini: flags give --dir '{sources}/ini', which holds a vep.ini ($DATA/sources/ini/vep.ini)" "$TMP/dd2.out.err" && ok "a --dir directory holding a vep.ini is refused by case, naming the file" || bad "vep.ini under --dir: $(grep -F 'c-dir-ini:' "$TMP/dd2.out.err")"
grep -qF "c-dir-ini-eq: flags give --dir '{sources}/ini', which holds a vep.ini" "$TMP/dd2.out.err" && ok "the --DIR=<directory> spelling over a vep.ini is refused the same way" || bad "vep.ini under --DIR=: $(grep -F 'c-dir-ini-eq:' "$TMP/dd2.out.err")"
grep -qF "c-dir-cache: flags names --dir_cache, which the runner sets" "$TMP/dd2.out.err" && ok "control: --dir_cache in a flags cell stays refused as the runner's own word" || bad "control --dir_cache: $(grep -F 'c-dir-cache:' "$TMP/dd2.out.err")"
grep -qF "c-dir-config: flags names --config, which the runner sets" "$TMP/dd2.out.err" && ok "control: --config in a flags cell stays refused as the runner's own word" || bad "control --config: $(grep -F 'c-dir-config:' "$TMP/dd2.out.err")"
grep -qF "c-dir-two: flags give --dir '{corpus}/{sources}/vepdir', $not_a_mount" "$TMP/dd2.out.err" && ok "a --dir directory carrying a second placeholder is refused: the engines would resolve two directories" || bad "--dir with two placeholders: $(grep -F 'c-dir-two:' "$TMP/dd2.out.err")"
grep -qF "c-word-dash: flags names -config, which the runner sets" "$TMP/dd2.out.err" && ok "runner_word_is_refused_in_every_spelling: -config is refused as --config is" || bad "-config: $(grep -F 'c-word-dash:' "$TMP/dd2.out.err")"
grep -qF "c-word-upper: flags names --CONFIG={sources}/ini/vep.ini, which the runner sets" "$TMP/dd2.out.err" && ok "--CONFIG=<file> is refused as --config is" || bad "--CONFIG=: $(grep -F 'c-word-upper:' "$TMP/dd2.out.err")"
grep -qF "c-word-dircache: flags names -dir_cache, which the runner sets" "$TMP/dd2.out.err" && ok "-dir_cache is refused as --dir_cache is" || bad "-dir_cache: $(grep -F 'c-word-dircache:' "$TMP/dd2.out.err")"
fc=$(cmd_after 'DRY-RUN c-word-free stub default' "$TMP/dd.out")
[[ "$fc" == *"--force_overwrite -HGVS +Symbol" ]] && ok "control: a word the runner does not set passes to the candidate as written in a single-dash upper-case or plus-sign spelling" || bad "control unlisted word: $fc"
# abbreviation_of_a_runner_word_is_refused: an engine reads an abbreviation as the one option
# whose names begin with it, and another name of an option as that option, so a word either
# engine's option list resolves to a word the runner sets is refused under that word: a word
# with a value, a switch, a word another option's name begins with (--json begins the
# candidate's --json_cache, so --js is --json by Ensembl VEP's list alone), a word of the
# candidate alone (--threads) and two other names (-u, --full_cache_dir).
{
    echo "$HEADER"
    row c-abbr-value fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--conf {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-value-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-conf {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-value-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--CONF={sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-switch fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--offl" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-switch-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-offl" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-switch-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--OFFL" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-begun fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--js" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-begun-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-js" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-begun-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--JS" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-dircache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir_c {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-candidate fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--thr 2" default "$IMAGE" offline 116.2 - - - - - - -
    row c-name-user fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-u anonymous" default "$IMAGE" offline 116.2 - - - - - - -
    row c-name-fullcache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--full_cache_dir {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
} >"$TD"
dry "$TD" "$TMP/da.out" --sources-dir "$DATA/sources"
rc=$?
[[ "$rc" -ne 0 ]] && ! grep -q '^DRY-RUN' "$TMP/da.out" && ok "abbreviation_of_a_runner_word_is_refused: a table of abbreviated runner words fails preflight, and nothing is planned" || bad "abbreviated runner words exited $rc with $(count_lines '^DRY-RUN' "$TMP/da.out") planned line(s)"
for case_word in "c-abbr-value|--conf|--config|a word with a value" "c-abbr-value-dash|-conf|--config|a word with a value, one dash" "c-abbr-value-upper|--CONF={sources}/ini/vep.ini|--config|a word with a value, upper case" \
    "c-abbr-switch|--offl|--offline|a switch" "c-abbr-switch-dash|-offl|--offline|a switch, one dash" "c-abbr-switch-upper|--OFFL|--offline|a switch, upper case" \
    "c-abbr-begun|--js|--json|a word another option's name begins with" "c-abbr-begun-dash|-js|--json|a word another option's name begins with, one dash" "c-abbr-begun-upper|--JS|--json|a word another option's name begins with, upper case" \
    "c-abbr-dircache|--dir_c|--dir_cache|the cache root" "c-abbr-candidate|--thr|--threads|a word of the candidate alone" \
    "c-name-user|-u|--user|another name of the option" "c-name-fullcache|--full_cache_dir|--json_cache|the candidate's other name of the option"; do
    IFS='|' read -r c word runner_word kind <<<"$case_word"
    grep -qF "$c: flags names $word, which vep reads as $runner_word, a word the runner sets" "$TMP/da.out.err" && ok "$word is refused as $runner_word ($kind)" || bad "$word: $(grep -F "$c:" "$TMP/da.out.err")"
done
# A word that abbreviates an option the runner does not set passes as written, and so does a
# prefix several options begin with, which neither engine reads as any of them: --vc begins
# --vcf and --vcf_info_field, --cach begins --cache and --cache_version.
{
    echo "$HEADER"
    row c-abbr-free fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--symb --dir_p {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-shared fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--vc" default "$IMAGE" offline 116.2 - - - - - - -
    row c-abbr-shared-two fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--cach" default "$IMAGE" offline 116.2 - - - - - - -
} >"$TD"
dry "$TD" "$TMP/df.out" --sources-dir "$DATA/sources"
rc=$?
pc=$(cmd_after 'DRY-RUN c-abbr-free perl default' "$TMP/df.out")
cc=$(cmd_after 'DRY-RUN c-abbr-free stub default' "$TMP/df.out")
[[ "$rc" -eq 0 && "$pc" == *"--force_overwrite --symb --dir_p /work/sources/gwas/vepdir" && "$cc" == *"--force_overwrite --symb --dir_p $DATA/sources/gwas/vepdir" ]] &&
    ok "control: an abbreviation of a word the runner does not set (--symb, --dir_p) plans and reaches both engines as written" || bad "control free abbreviation exited $rc: perl=$pc cand=$cc"
pc=$(cmd_after 'DRY-RUN c-abbr-shared perl default' "$TMP/df.out")
cc=$(cmd_after 'DRY-RUN c-abbr-shared stub default' "$TMP/df.out")
[[ "$pc" == *"--force_overwrite --vc" && "$cc" == *"--force_overwrite --vc" ]] &&
    ok "control: a prefix a runner word shares with another option (--vc: --vcf, --vcf_info_field) is not refused as the runner word and reaches both engines as written" || bad "control shared prefix: perl=$pc cand=$cc"
pc=$(cmd_after 'DRY-RUN c-abbr-shared-two perl default' "$TMP/df.out")
cc=$(cmd_after 'DRY-RUN c-abbr-shared-two stub default' "$TMP/df.out")
[[ "$pc" == *"--force_overwrite --cach" && "$cc" == *"--force_overwrite --cach" ]] &&
    ok "control: a prefix two runner words share (--cach: --cache, --cache_version) names neither and reaches both engines as written" || bad "control prefix of two runner words: perl=$pc cand=$cc"
# option_lists_are_the_engines_own: the two lists the refusal reads are held to their sources,
# Ensembl VEP's `@VEP_PARAMS` as crates/vep-cli/src/argv.rs carries it and the candidate's
# parser as the reference pages generated from it list it (the group pages the landing page
# names: each flag row's names, and a flag the error table catches before parsing), so an
# option either engine gains fails here until the list carries it.
perl_list=$(sed -n 's/^PERL_OPTIONS="\(.*\)"$/\1/p' "$ORACLE")
perl_source=$(sed -n '/^pub const VEP_PARAMS: &\[&str\] = &\[$/,/^\];$/p' "$REPO/crates/vep-cli/src/argv.rs" | sed -n 's/^    "\([^"=:]*\)[^"]*",$/\1/p' | tr '\n' ' ')
[[ -n "$perl_list" && "$perl_list" == "${perl_source% }" ]] && ok "option_lists_are_the_engines_own: PERL_OPTIONS is VEP_PARAMS of argv.rs, name for name and in its order" ||
    bad "PERL_OPTIONS differs from VEP_PARAMS of argv.rs: $(diff <(tr ' ' '\n' <<<"$perl_list") <(tr ' ' '\n' <<<"${perl_source% }") | grep '^[<>]' | head -5 | tr '\n' ' ')"
candidate_list=$(sed -n 's/^CANDIDATE_OPTIONS="\(.*\)"$/\1/p' "$ORACLE")
candidate_source=$(sed -n 's/^| \[[^]]*\](\(cli\/[a-z-]*\.md\)) | .*/\1/p' "$REPO/docs/cli-reference.md" | while IFS= read -r page; do
    awk '
        /^\| Flag \| Type \| Default \| Description \|$/ { kind = 1; next }
        /^\| Flag \| / { kind = 2; next }
        /^\| --- / { next }
        !/^\|/ { kind = 0; next }
        kind == 0 || (kind == 2 && !/^\| `--[A-Za-z0-9_]+` \(not declared; caught before parsing\) \| /) { next }
        {
            cell = $0
            sub(/ \| .*$/, "", cell)
            names = ""
            while (match(cell, /`-+[A-Za-z0-9_]+`/)) {
                name = substr(cell, RSTART + 1, RLENGTH - 2)
                sub(/^-+/, "", name)
                names = names (names == "" ? "" : "|") name
                cell = substr(cell, RSTART + RLENGTH)
            }
            if (names != "") print names
        }' "$REPO/docs/$page"
done | LC_ALL=C sort -u | tr '\n' ' ')
[[ -n "$candidate_list" && "$candidate_list" == "${candidate_source% }" ]] && ok "CANDIDATE_OPTIONS is the parser's options as the reference pages under docs/cli/ list them, sorted" ||
    bad "CANDIDATE_OPTIONS differs from the reference pages: $(diff <(tr ' ' '\n' <<<"$candidate_list") <(tr ' ' '\n' <<<"${candidate_source% }") | grep '^[<>]' | head -5 | tr '\n' ' ')"
unlisted=""
runner_list=$(sed -n 's/^RUNNER_FLAGS="\(.*\)"$/\1/p' "$ORACLE")
for word in $runner_list; do
    name="${word#-}"
    [[ " $(tr '|' ' ' <<<"$perl_list $candidate_list") " == *" ${name#-} "* ]] || unlisted="$unlisted $word"
done
[[ -n "$perl_list" && -z "$unlisted" ]] && ok "every word of RUNNER_FLAGS is a name of an option in one of the two lists" || bad "RUNNER_FLAGS words neither list names:${unlisted:- (PERL_OPTIONS is empty)}"

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
    row c-miss fixture "$DATA/absent.vcf" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a missing case input exits non-zero and names it" "c-miss: input not found: $DATA/absent.vcf" "$T3"
grep -q 'preflight failed; nothing ran' "$TMP/e.out.err" && ok "the abort is reported as a preflight failure (the other case was fine)" || bad "no preflight-failure line"
{
    echo "$HEADER"
    row c-fa fixture "$CORPUS_INPUT" "$CORPUS_CACHE" "$DATA/absent.fa" - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a missing FASTA is an error" "FASTA not found" "$T3"
{
    echo "$HEADER"
    row c-cache fixture "$CORPUS_INPUT" 116_GRCh38 - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a missing Ensembl cache slice is an error in offline mode" "Ensembl VEP cache slice not found: $DATA/x/homo_sapiens/116_GRCh38" "$T3" --perl-cache-dir "$DATA/x"
expect_error "a candidate JSON cache without transcripts/ is an error" "candidate JSON cache has no transcripts/" "$T3" --json-cache-dir "$DATA/x"
{
    echo "$HEADER"
    row c-fmt fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default,xml "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "an unknown format is an error" "formats names 'xml'" "$T3"
{
    echo "$HEADER"
    row c-w fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab,warnings "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "warnings without default is an error" "warnings without default" "$T3"
{
    echo "$HEADER"
    row c-flag fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--offline --hgvs" default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a flags cell naming a runner-owned flag is an error" "flags names --offline" "$T3"
{
    echo "$HEADER"
    row c-mode fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" online 116.2 - - - - - - -
} >"$T3"
expect_error "an unknown perl_mode is an error" "perl_mode 'online'" "$T3"
{
    echo "$HEADER"
    row c-tier nightly "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "an unknown tier is an error" "tier 'nightly'" "$T3"
{
    echo "$HEADER"
    row c-dup fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
    row c-dup fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a duplicate case is an error" "case 'c-dup' appears twice" "$T3"
{
    echo "$HEADER"
    row "../escape" fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a case that is not a directory-safe name is an error (it names the work and report directories)" "case '../escape' is not a directory-safe name" "$T3"
{
    echo "$HEADER"
    row c-empty fixture "$CORPUS_INPUT" "$CORPUS_CACHE" "" - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "an empty cell is an error (write -)" "an empty cell" "$T3"
{
    echo "$HEADER"
    row c-short fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - -
} >"$T3"
expect_error "a short row is an error" "16 columns, expected 17" "$T3"
{
    row case input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs species unordered reference refusal
    row c-x "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a header without the tier column is an error" "header lacks the column 'tier'" "$T3"
{
    row case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs unordered reference refusal
    row c-x fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - -
} >"$T3"
expect_error "a header without the species column is an error" "header lacks the column 'species'" "$T3"
{
    echo "$HEADER"
    row c-suffixed fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - homo_sapiens_refseq - - -
} >"$TS"
expect_error "a species carrying the slice suffix is refused, since Ensembl VEP refuses it as --species" "species 'homo_sapiens_refseq' carries the slice suffix; write the species and the --refseq flag" "$TS"
{
    echo "$HEADER"
    row c-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - Homo_sapiens - - -
} >"$TS"
expect_error "a species that is not a lower-case production name is refused" "species 'Homo_sapiens' is not an Ensembl species name" "$TS"
{
    echo "$HEADER"
    row c-masks fixture "$DATA/caches/vep-rs/grch38/transcripts" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 some_class - - - - - -
} >"$T3"
expect_error "masks without a manifest.json beside the input is an error" "masks names 'some_class', which applies from a manifest.json beside the input, and none is there (the built-in classes are regulatory_feature_rows order_only individual_allele_order conservation_short_name_spelling pc_description_text motif_score_reverse_strand_mnv colocated_list_cut_at_cache_cap start_cooccurrence_swap splice_family_swap gxf_reverse_strand_utr_sequence organelle_peptide_paths_translation_table conservation_max_text_order gwas_haplotype_cross_product gwas_interaction_first_snp_allele_whitespace alphamissense_two_digit_transcript_version cadd_structural_allele_match splicevault_sibling_allele_rows splicevault_json_later_rows utrannotator_annotation_pair_order utrannotator_gained_start_window utrannotator_consequence_precedence utrannotator_structural_allele_text utrannotator_lost_start_after_insertion utrannotator_gained_start_remade_in_place utrannotator_frameshift_at_stop_boundary utrannotator_lost_stop_written_before_stop_codon nearest_exon_jb_item_order loeuf_zero_minimum maxentscan_overlapping_windows alphamissense_sibling_allele_rows revel_sibling_allele_rows revel_other_reference_amino_acid_rows gnomadmt_info_value_second_equals loftee_anc_allele_multiallelic_record loftee_first_alt_donor_filter loftee_gerp_zoom_level_summary nearest_exon_jb_stretch_priming nearest_exon_jb_shifted_span nearest_exon_jb_id_less_exon_fold downstream_reverse_strand_rotation spliceai_chopped_multiallelic spliceai_multi_entry_cell plugin_contig_mitochondrial_spelling refseq_edited_model_hgvsc_offset refseq_edited_model_shifted_insertion_bases refseq_edited_model_seam_insertion_hgvsc refseq_edited_model_unminimised_allele_hgvsc refseq_edited_model_hgvs_offset refseq_edited_model_deleted_span_term refseq_edited_model_tail_term mature_mirna_unshifted_span stop_codon_unshifted_span stop_retained_without_stop_codon_test stop_codon_overlap_genomic_span stop_retained_three_base_stop_codon_test stop_codon_removed_without_stop_codon_test stop_codon_overlap_insertion_range splice_term_past_the_record_span refseq_edited_model_deleted_span_codons hgvs_intron_offset_abutting_insertion hgvs_shifted_insertion_clamped_end shift_intron_term_unshifted_span hgvs_shifted_indel_written_span_gate hgvsp_shifted_indel_written_span_gate hgvs_shifted_deletion_clamped_span hgvs_clipped_reference_against_whole_alternate hgvs_added_bases_left_at_the_span_end hgvs_insertion_left_where_the_shift_stops hgvs_added_bases_written_outside_the_transcript hgvs_shift_sequence_edge_window hgvs_shift_shared_across_alleles hgvsp_stop_codon_protein_comparison hgvsp_stop_named_on_incomplete_coding_sequence hgvsp_intron_end_insertion hgvsp_intron_start_insertion frequency_filter_zero_value merged_same_stable_id_collapse hgnc_id_copy downstream_utr_frame_past_coding_end downstream_start_past_own_codon plugin_hgvs_shifted_cds_coordinates plugin_value_follows_splice_term_past_the_record_span plugin_value_follows_stop_codon_term_pair plugin_value_on_reference_skipped_record blosum62_hgvs_shift_shared_across_alleles reference_skipped_record breakend_mate_local_read bnd_own_chromosome_rows_lost dosage_sensitivity_regulatory_row_warnings fork_warning_pid_prefix)" "$T3"
mkdir -p "$DATA/nomanifest"
cp "$CORPUS/variants.vcf" "$DATA/nomanifest/variants.vcf"
{
    echo "$HEADER"
    row c-builtin fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,regulatory_feature_rows,individual_allele_order,conservation_short_name_spelling,pc_description_text - - - - - -
    row c-mixed fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,some_class - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-builtin.out" --only c-builtin
rc=$?
bc=$(cmd_second 'DRY-RUN c-builtin stub default' "$TMP/d-builtin.out")
[[ "$rc" -eq 0 && "$bc" == *'--masks order_only\,regulatory_feature_rows\,individual_allele_order\,conservation_short_name_spelling\,pc_description_text' && "$bc" != *"--manifest"* ]] && ok "the comparator's built-in classes need no manifest beside the input" || bad "built-in masks: rc=$rc $bc"
[[ "$bc" != *" --input "* ]] && ok "a default comparison without a manifest or an input-reading class passes no --input" || bad "input leaked without a manifest: $bc"
expect_error "a manifest class beside a built-in one still needs the manifest" "masks names 'some_class', which applies from a manifest.json" "$T3" --only c-mixed
# The runner's built-in list is the comparator's: the row and text classes, the registry pair
# classes, the field classes, the transcript row classes and the warnings-mode classes, no more
# and no fewer.
builtin_listed=$(awk -F'"' '/^BUILTIN_MASKS="/ {print $2; exit}' "$ORACLE" | tr ' ' '\n' | sort)
builtin_implemented=$(python3 - "$HERE" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
import compare_vep_formats as c
names = {c.REGULATORY_FEATURE_ROWS, c.ORDER_ONLY, c.INDIVIDUAL_ALLELE_ORDER, c.CONSERVATION_SHORT_NAME_SPELLING,
         c.PC_DESCRIPTION_TEXT, c.MOTIF_SCORE_REVERSE_STRAND_MNV, c.MERGED_SAME_STABLE_ID_COLLAPSE,
         c.COLOCATED_LIST_CUT_AT_CACHE_CAP,
         *c.TERM_PAIR_CLASSES, *c.REGISTRY_CLASSES, *c.FIELD_CLASSES, *c.STRUCTURAL_VARIANT_CLASSES,
         *c.PLUGIN_WARNING_CLASSES, *c.WARNING_TEXT_CLASSES}
print("\n".join(sorted(names)))
PY
)
[[ "$builtin_listed" == "$builtin_implemented" ]] && ok "BUILTIN_MASKS names exactly the comparator's built-in classes" || bad "BUILTIN_MASKS drift: listed [$(tr '\n' ' ' <<<"$builtin_listed")] implemented [$(tr '\n' ' ' <<<"$builtin_implemented")]"
# The classes that read beside the two outputs. A class that reads a record's alleles, or the
# reference bases a row's allele stands in for, brings --input to a text-format comparison
# without a manifest. The two GWAS classes set nothing
# aside without the catalog, so their comparison names the file the case's flags hand the GWAS
# plugin, the host path, in either spelling of the value. hgvsp_intron_end_insertion and
# hgvsp_intron_start_insertion set nothing aside without the candidate's JSON cache, so their
# comparison names the case's cache.
TB="$TMP/tb.tsv"
{
    echo "$HEADER"
    row c-alleles fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 revel_sibling_allele_rows - - - - - -
    row c-span fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 splice_term_past_the_record_span - - - - - -
    row c-gwas fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv --plugin GWAS,file={corpus}/catalog.tsv" default "$IMAGE" offline 116.2 gwas_haplotype_cross_product=2 - - - - - -
    row c-gwas-joined fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin=GWAS,type=curated,file={corpus}/catalog.tsv" tab "$IMAGE" offline 116.2 gwas_interaction_first_snp_allele_whitespace - - - - - -
    row c-gwas-unmasked fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin GWAS,file={corpus}/catalog.tsv" default "$IMAGE" offline 116.2 order_only - - - - - -
    row c-gwas-fileless fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv" default "$IMAGE" offline 116.2 gwas_haplotype_cross_product - - - - - -
    row c-clipped fixture nomanifest/variants.vcf "$CORPUS_CACHE" - --hgvs default "$IMAGE" offline 116.2 hgvs_clipped_reference_against_whole_alternate - - - - - -
    row c-downstream fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin Downstream" default,tab "$IMAGE" offline 116.2 downstream_utr_frame_past_coding_end,downstream_start_past_own_codon - - - - - -
    row c-organelle fixture nomanifest/variants.vcf "$CORPUS_CACHE" - --hgvs default "$IMAGE" offline 116.2 organelle_peptide_paths_translation_table - - - - - -
    row c-edge fixture nomanifest/variants.vcf 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 hgvs_shift_sequence_edge_window - - - - - -
    row c-edge-manifest fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 hgvs_shift_sequence_edge_window - - - - - -
    row c-intron-end fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--hgvs" default "$IMAGE" offline 116.2 hgvsp_intron_end_insertion - - - - - -
    row c-intron-start fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--hgvs" default "$IMAGE" offline 116.2 hgvsp_intron_start_insertion - - - - - -
    row c-two-windows fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 maxentscan_overlapping_windows - - - - - -
} >"$TB"
dry "$TB" "$TMP/d-beside.out"
rc=$?
input_listed=$(awk -F'"' '/^INPUT_READING_MASKS="/ {print $2; exit}' "$ORACLE" | tr ' ' '\n' | sort)
input_implemented=$(python3 - "$HERE" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
import compare_vep_formats as c
print("\n".join(sorted({c.REFERENCE_SKIPPED_RECORD, *c.ALLELE_READING_CLASSES, *c.RECORD_REFERENCE_CLASSES})))
PY
)
ac=$(cmd_second 'DRY-RUN c-alleles stub default' "$TMP/d-beside.out")
sc=$(cmd_second 'DRY-RUN c-span stub default' "$TMP/d-beside.out")
[[ "$input_listed" == "$input_implemented" && "$ac" == *" --input $DATA/nomanifest/variants.vcf --masks revel_sibling_allele_rows" && "$sc" == *" --input $DATA/nomanifest/variants.vcf --masks splice_term_past_the_record_span" ]] &&
    ok "INPUT_READING_MASKS names the skipped-record class and the comparator's allele-reading and reference-reading classes, whose text-format comparison carries --input" || bad "INPUT_READING_MASKS: listed [$(tr '\n' ' ' <<<"$input_listed")] implemented [$(tr '\n' ' ' <<<"$input_implemented")] default=$ac span=$sc"
gc=$(cmd_second 'DRY-RUN c-gwas stub default' "$TMP/d-beside.out")
gj=$(cmd_second 'DRY-RUN c-gwas-joined stub tab' "$TMP/d-beside.out")
[[ "$rc" -eq 0 && "$gc" == *"--masks gwas_haplotype_cross_product --gwas-catalog $DATA/nomanifest/catalog.tsv" && "$gj" == *"--masks gwas_interaction_first_snp_allele_whitespace --gwas-catalog $DATA/nomanifest/catalog.tsv" ]] &&
    ok "a comparison naming a GWAS class carries the GWAS plugin's file as --gwas-catalog" || bad "GWAS catalog: rc=$rc default=$gc tab=$gj"
gu=$(cmd_second 'DRY-RUN c-gwas-unmasked stub default' "$TMP/d-beside.out")
gf=$(cmd_second 'DRY-RUN c-gwas-fileless stub default' "$TMP/d-beside.out")
[[ "$gu" == *"--masks order_only" && "$gf" == *"--masks gwas_haplotype_cross_product" ]] &&
    ok "a comparison naming no GWAS class, or a GWAS class beside no GWAS plugin file, carries no --gwas-catalog" || bad "GWAS catalog leaked: unmasked=$gu fileless=$gf"
# The sequence-edge class is decided from the reference both engines read and the candidate's
# cache: its comparison names the case's FASTA, decompressed as the engines read it, and the
# cache. A case with a manifest and no FASTA leaves the class to the manifest, and a case with
# neither is refused before anything runs.
reference_listed=$(awk -F'"' '/^REFERENCE_READING_MASKS="/ {print $2; exit}' "$ORACLE")
ec=$(cmd_second 'DRY-RUN c-edge stub default' "$TMP/d-beside.out")
[[ "$reference_listed" == "hgvs_shift_sequence_edge_window" && "$ec" == *"--masks hgvs_shift_sequence_edge_window --json-cache $REPO/tests/golden/116/GRCh37-hgvs/json_cache --reference-fasta $TMP/never/work/fasta/1/reference.fa" ]] &&
    ok "a comparison naming the sequence-edge class carries the candidate's cache and the case's FASTA as --reference-fasta" || bad "sequence-edge reference: listed [$reference_listed] default=$ec"
em=$(cmd_second 'DRY-RUN c-edge-manifest stub default' "$TMP/d-beside.out")
[[ "$em" == *"--manifest $CORPUS/manifest.json --reference-sequence absent --masks hgvs_shift_sequence_edge_window" ]] &&
    ok "a case that names the sequence-edge class beside a manifest and no FASTA carries neither --reference-fasta nor --json-cache" || bad "sequence-edge beside a manifest: $em"
# Beside a manifest the comparator is told the candidate's run: its run of the manifest, and
# whether it has a reference sequence, which the candidate takes from --fasta or from an
# indexed file inside its cache directory. A fabricated corpus with its own run, two flag
# sets and a manifest side that names --grch37 and a --fasta word, and candidate caches that
# hold each form of reference file.
RS="$DATA/runsetting"
mkdir -p "$RS/corpus" "$RS/bare/transcripts" "$RS/plain/transcripts" "$RS/unindexed/transcripts" "$RS/bgzf/transcripts" \
    "$RS/bgzf-no-gzi/transcripts" "$RS/gzip/transcripts" "$RS/other-name/transcripts"
printf '21\t100\trs1\tA\tG\t.\t.\t.\n' >"$RS/corpus/variants.vcf"
cp "$RS/corpus/variants.vcf" "$RS/corpus/second.vcf"
cat >"$RS/corpus/manifest.json" <<'JSON'
{"flags": ["--grch37", "--symbol", "--hgvs"],
 "flag_sets": [
  {"name": "picked", "flags": ["--grch37", "--symbol", "--pick", "--fasta", "{corpus}/ref.fa"]},
  {"name": "second-input", "input": "second.vcf", "flags": ["--symbol", "--hgvs"]},
  {"name": "twin-a", "flags": ["--twin"]},
  {"name": "twin-b", "flags": ["--twin", "--fasta=x.fa"]}]}
JSON
printf '>21\nA\n' >"$RS/ref.fa"
printf '21\t1\t4\t1\t2\n' >"$RS/ref.fa.fai"
bgzf_head() { printf '\037\213\010\004\000\000\000\000\000\377\006\000BC\002\000\033\000'; }
printf '>21\nA\n' >"$RS/plain/ref.fa"
cp "$RS/ref.fa.fai" "$RS/plain/ref.fa.fai"
printf '>21\nA\n' >"$RS/unindexed/ref.fa"
bgzf_head >"$RS/bgzf/ref.fa.gz"
cp "$RS/ref.fa.fai" "$RS/bgzf/ref.fa.gz.fai"
: >"$RS/bgzf/ref.fa.gz.gzi"
bgzf_head >"$RS/bgzf-no-gzi/ref.fa.gz"
cp "$RS/ref.fa.fai" "$RS/bgzf-no-gzi/ref.fa.gz.fai"
printf '>21\nA\n' | gzip -c >"$RS/gzip/ref.fa.gz"
cp "$RS/ref.fa.fai" "$RS/gzip/ref.fa.gz.fai"
: >"$RS/gzip/ref.fa.gz.gzi"
printf '>21\nA\n' >"$RS/other-name/ref.fasta"
cp "$RS/ref.fa.fai" "$RS/other-name/ref.fasta.fai"
T_RS="$TMP/t-runsetting.tsv"
rs_row() { # rs_row <case> <input> <cache dir under runsetting/> <fasta cell> <flags cell>
    row "$1" fixture "runsetting/corpus/$2" "116_GRCh38=runsetting/$3" "$4" "$5" default "$IMAGE" offline 116.2 - - - - - - -
}
{
    echo "$HEADER"
    rs_row c-rs-own variants.vcf bare - "--hgvs --symbol"
    rs_row c-rs-set variants.vcf bare - "--pick --symbol"
    rs_row c-rs-second second.vcf bare - "--symbol --hgvs"
    rs_row c-rs-none variants.vcf bare - "--symbol"
    rs_row c-rs-twin variants.vcf bare - "--twin"
    rs_row c-rs-cell variants.vcf bare runsetting/ref.fa "--hgvs --symbol"
    rs_row c-rs-cell-no-fasta variants.vcf plain runsetting/ref.fa "--symbol --no_fasta"
    rs_row c-rs-plain variants.vcf plain - "--hgvs --symbol"
    rs_row c-rs-unindexed variants.vcf unindexed - "--hgvs --symbol"
    rs_row c-rs-bgzf variants.vcf bgzf - "--hgvs --symbol"
    rs_row c-rs-bgzf-no-gzi variants.vcf bgzf-no-gzi - "--hgvs --symbol"
    rs_row c-rs-gzip variants.vcf gzip - "--hgvs --symbol"
    rs_row c-rs-other-name variants.vcf other-name - "--hgvs --symbol"
    rs_row c-rs-no-fasta variants.vcf plain - "--symbol -NO_FASTA"
    rs_row c-rs-fasta-dir variants.vcf plain - "--symbol --fasta_dir {corpus}"
    row c-rs-no-manifest fixture nomanifest/variants.vcf "116_GRCh38=runsetting/plain" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T_RS"
dry "$T_RS" "$TMP/d-rs.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "the run and reference setting table plans" || bad "run setting dry run exited $rc: $(head -3 "$TMP/d-rs.out.err" | tr '\n' ' ')"
rs_cmd() { cmd_second "DRY-RUN $1 stub default" "$TMP/d-rs.out"; }
rs_told() { # rs_told <case> <words after --manifest, or nothing>: the comparator command names exactly those words there
    local c
    c=$(rs_cmd "$1")
    if [[ -n "$2" ]]; then
        [[ "$c" == *" --manifest $RS/corpus/manifest.json $2" ]]
    else
        [[ "$c" == *" --manifest $RS/corpus/manifest.json" ]]
    fi
}
rs_told c-rs-own "--run default --reference-sequence absent" &&
    ok "a case whose flags are the manifest's own run's, in another order and without its --grch37, is told --run default" || bad "own run: $(rs_cmd c-rs-own)"
rs_told c-rs-set "--run picked --reference-sequence absent" &&
    ok "a case whose flags are a flag set's, the set's --fasta word and its file set aside, is told that set's name" || bad "flag set run: $(rs_cmd c-rs-set)"
[[ "$(rs_cmd c-rs-second)" == *" --input $RS/corpus/second.vcf --manifest $RS/corpus/manifest.json --run second-input --reference-sequence absent" ]] &&
    ok "a run is one whose input has the case input's file name: the flag set over another input, not the own run of the same flags" || bad "run by input: $(rs_cmd c-rs-second)"
rs_told c-rs-none "--reference-sequence absent" &&
    ok "a case that is no run of its manifest is told no run" || bad "no run: $(rs_cmd c-rs-none)"
rs_told c-rs-twin "--reference-sequence absent" &&
    ok "a case that fits two runs of its manifest is told no run" || bad "two runs: $(rs_cmd c-rs-twin)"
rs_told c-rs-cell "--run default --reference-sequence present" &&
    ok "a case with a FASTA cell is told the run has a reference sequence" || bad "FASTA cell: $(rs_cmd c-rs-cell)"
rs_told c-rs-cell-no-fasta "--reference-sequence present" &&
    ok "--no_fasta beside a FASTA cell leaves the reference sequence present: the candidate is given --fasta" || bad "FASTA cell under --no_fasta: $(rs_cmd c-rs-cell-no-fasta)"
rs_told c-rs-plain "--run default --reference-sequence present" &&
    ok "a .fa with its .fai inside the candidate's cache directory is a reference sequence" || bad "cache .fa: $(rs_cmd c-rs-plain)"
rs_told c-rs-unindexed "--run default --reference-sequence absent" &&
    ok "a .fa without its .fai inside the cache directory is none" || bad "cache .fa without .fai: $(rs_cmd c-rs-unindexed)"
rs_told c-rs-bgzf "--run default --reference-sequence present" &&
    ok "a BGZF .fa.gz with its .fai and .gzi inside the cache directory is a reference sequence" || bad "cache BGZF: $(rs_cmd c-rs-bgzf)"
rs_told c-rs-bgzf-no-gzi "--run default --reference-sequence absent" &&
    ok "a BGZF .fa.gz without its .gzi is none" || bad "cache BGZF without .gzi: $(rs_cmd c-rs-bgzf-no-gzi)"
rs_told c-rs-gzip "--run default --reference-sequence absent" &&
    ok "a gzip .fa.gz that is not BGZF is none, whatever index files lie beside it" || bad "cache gzip: $(rs_cmd c-rs-gzip)"
rs_told c-rs-other-name "--run default --reference-sequence absent" &&
    ok "a file of another name than .fa or .fa.gz is none" || bad "cache other name: $(rs_cmd c-rs-other-name)"
rs_told c-rs-no-fasta "--reference-sequence absent" &&
    ok "--no_fasta, in any spelling vep reads, keeps the candidate from the cache directory's file: absent" || bad "--no_fasta: $(rs_cmd c-rs-no-fasta)"
rs_told c-rs-fasta-dir "" &&
    ok "a flags cell naming --fasta_dir tells the comparator no reference setting" || bad "--fasta_dir: $(rs_cmd c-rs-fasta-dir)"
nm=$(rs_cmd c-rs-no-manifest)
[[ "$nm" != *"--run"* && "$nm" != *"--reference-sequence"* && "$nm" != *"--manifest"* ]] &&
    ok "a case without a manifest beside its input tells the comparator neither" || bad "no manifest: $nm"
# The candidate's own rule, name for name: the .fa and .fa.gz names, the .fai, and the .gzi of BGZF bytes.
fasta_rule="$REPO/crates/vep-fasta/src/lib.rs"
grep -q 'name.ends_with(".fa") || name.ends_with(".fa.gz")' "$fasta_rule" && grep -q 'sibling(path, "fai").is_file()' "$fasta_rule" &&
    grep -q 'Ok(Encoding::Bgzf) => sibling(path, "gzi").is_file(),' "$fasta_rule" && grep -q 'Ok(Encoding::Gzip) | Err(_) => false,' "$fasta_rule" &&
    grep -q 'let fasta = if config.fasta.is_some() || config.no_fasta {' "$REPO/crates/vep-cli/src/cache_resolve.rs" &&
    ok "the candidate's detection reads the names and index files the runner reads" || bad "the candidate's FASTA detection moved: read is_fasta_name, is_indexed and resolve again"
{
    echo "$HEADER"
    row c-edge-bare fixture nomanifest/variants.vcf "$CORPUS_CACHE" - --hgvs default "$IMAGE" offline 116.2 hgvs_shift_sequence_edge_window - - - - - -
} >"$T3"
expect_error "the sequence-edge class beside neither a FASTA nor a manifest is an error" "masks names 'hgvs_shift_sequence_edge_window', which is decided from the case's FASTA or applies from a manifest.json beside the input, and the case has neither" "$T3"
# The organelle peptide class reads its transcripts' codon tables from the candidate's cache.
oc=$(cmd_second 'DRY-RUN c-organelle stub default' "$TMP/d-beside.out")
[[ "$oc" == *"--masks organelle_peptide_paths_translation_table --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming the organelle peptide class carries the candidate's cache as --json-cache" || bad "organelle peptide cache: $oc"
# The two Downstream classes recompute both values from the candidate's cache.
dc=$(cmd_second 'DRY-RUN c-downstream stub default' "$TMP/d-beside.out")
dt=$(cmd_second 'DRY-RUN c-downstream stub tab' "$TMP/d-beside.out")
[[ "$dc" == *"--masks downstream_utr_frame_past_coding_end\,downstream_start_past_own_codon --json-cache $CORPUS/json_cache" && "$dt" == *"--json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming the Downstream classes carries the candidate's cache as --json-cache" || bad "Downstream cache: default=$dc tab=$dt"
# The class of an allele over a transcript's edge counts the exon spans of the candidate's cache.
kc=$(cmd_second 'DRY-RUN c-clipped stub default' "$TMP/d-beside.out")
[[ "$kc" == *"--masks hgvs_clipped_reference_against_whole_alternate --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming the class of an allele over a transcript's edge carries the candidate's cache as --json-cache" || bad "transcript-edge cache: $kc"
ic=$(cmd_second 'DRY-RUN c-intron-end stub default' "$TMP/d-beside.out")
[[ "$ic" == *"--masks hgvsp_intron_end_insertion --json-cache $CORPUS/json_cache" && "$ac" != *"--json-cache"* ]] &&
    ok "a comparison naming hgvsp_intron_end_insertion carries the candidate's JSON cache as --json-cache, one naming no cache-reading class none" || bad "intron-end cache: named=$ic unnamed=$ac"
isc=$(cmd_second 'DRY-RUN c-intron-start stub default' "$TMP/d-beside.out")
[[ "$isc" == *"--masks hgvsp_intron_start_insertion --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming hgvsp_intron_start_insertion carries the candidate's JSON cache as --json-cache" || bad "intron-start cache: $isc"
# maxentscan_overlapping_windows reads a transcript's introns from the candidate's cache.
tw=$(cmd_second 'DRY-RUN c-two-windows stub default' "$TMP/d-beside.out")
[[ "$tw" == *"--masks maxentscan_overlapping_windows --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming maxentscan_overlapping_windows carries the candidate's JSON cache as --json-cache" || bad "two-window cache: $tw"
# A class may carry its expected applied count: the count stays with the runner's report check
# and the comparator receives the class names alone; an item that is neither form is refused.
{
    echo "$HEADER"
    row c-counted fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 start_cooccurrence_swap=127,splice_family_swap,conservation_max_text_order=3 - - - - - -
    row c-badcount fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=x - - - - - -
    row c-doubleeq fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=1=2 - - - - - -
    row c-zero fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=0 - - - - - -
    row c-empty fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,,regulatory_feature_rows - - - - - -
    row c-trail fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=1, - - - - - -
    row c-blank fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 "order_only, regulatory_feature_rows" - - - - - -
    row c-glob fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 "*" - - - - - -
    row c-twice fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=1,regulatory_feature_rows,order_only=3 - - - - - -
    row c-countedw fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default,warnings "$IMAGE" offline 116.2 start_cooccurrence_swap=127,dosage_sensitivity_regulatory_row_warnings - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-counted.out" --only c-counted
rc=$?
bc=$(cmd_second 'DRY-RUN c-counted stub default' "$TMP/d-counted.out")
[[ "$rc" -eq 0 && "$bc" == *'--masks start_cooccurrence_swap\,splice_family_swap\,conservation_max_text_order' && "$bc" != *"=127"* ]] && ok "a masks item's expected count is stripped before the comparator sees the classes" || bad "counted masks: rc=$rc $bc"
expect_error "a masks count that is not a number is an error" "masks item 'order_only=x' is not <class> or <class>=<expected applied count>" "$T3" --only c-badcount
expect_error "a masks item with two = is an error" "masks item 'order_only=1=2' is not <class> or <class>=<expected applied count>" "$T3" --only c-doubleeq
expect_error "a masks count of zero is an error (a class applied to nothing is red as unused)" "masks item 'order_only=0' is not <class> or <class>=<expected applied count>" "$T3" --only c-zero
expect_error "an empty masks item between two commas is an error" "masks 'order_only,,regulatory_feature_rows' has an empty or blank item" "$T3" --only c-empty
! grep -q "twice\|masks names ''" "$TMP/e.out.err" && ok "an empty item draws its one refusal, not a duplicate or built-in-list message for a class with no name" || bad "empty item noise: $(grep -c '' "$TMP/e.out.err") stderr lines"
expect_error "a trailing comma in the masks cell is an error" "masks 'order_only=1,' has an empty or blank item" "$T3" --only c-trail
expect_error "a blank inside the masks cell is an error (the report splits the cell on commas alone)" "masks 'order_only, regulatory_feature_rows' has an empty or blank item" "$T3" --only c-blank
# The runner reads the table from the directory holding a file a glob item would name.
mkdir -p "$TMP/globdir" && touch "$TMP/globdir/order_only"
pushd "$TMP/globdir" >/dev/null || exit 1
expect_error "a glob character in a masks item is refused as written, not read as a file name of the working directory" "masks item '[*]' is not <class> or <class>=<expected applied count>" "$T3" --only c-glob
popd >/dev/null || exit 1
expect_error "a class named twice in the masks cell is an error (two counts could not both hold)" "masks names 'order_only' twice" "$T3" --only c-twice
dry "$T3" "$TMP/d-countedw.out" --only c-countedw
rc=$?
wcm=$(cmd_after 'DRY-RUN c-countedw stub warnings' "$TMP/d-countedw.out")
[[ "$rc" -eq 0 && "$wcm" == *'--warnings --expected '* && "$wcm" == *'--masks start_cooccurrence_swap\,dosage_sensitivity_regulatory_row_warnings' ]] && ok "the warnings comparison receives the case's classes, counts stripped, so a warnings-mode class can apply there" || bad "warnings masks: rc=$rc $wcm"
{
    echo "$HEADER"
    row c-dash fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 -,order_only - - - - - -
} >"$T3"
expect_error "a - beside a class in the masks cell is an error (the comparator would read it as an option)" "masks '-,order_only' has an item beginning with -" "$T3"
{
    echo "$HEADER"
    row c-dashfloor fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - -,IMPACT=1 - - - - -
    row c-dashpair fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - "-;IMPACT=LOW=>HIGH" - - - -
    row c-dashclass fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 -foo - - - - - -
    row c-dashafter fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,-foo - - - - - -
} >"$T3"
expect_error "a - item in the floors cell is an error" "floors '-,IMPACT=1' has an item beginning with -" "$T3" --only c-dashfloor
expect_error "a - item in the field_pairs cell is an error" "field_pairs '-;IMPACT=LOW=>HIGH' has an item beginning with -" "$T3" --only c-dashpair
expect_error "a masks class beginning with - is an error even beside a manifest" "masks '-foo' has an item beginning with -" "$T3" --only c-dashclass
expect_error "a masks class beginning with - after a separator is an error" "masks 'order_only,-foo' has an item beginning with -" "$T3" --only c-dashafter
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
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
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json,warnings "$IMAGE" offline 116.2 - - - - - - -
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
    row c-floors fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,json "$IMAGE" offline 116.2 - Uploaded_variation=1,IMPACT=1 - - - - -
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
    row c-jsonfloor fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" json "$IMAGE" offline 116.2 - Uploaded_variation=1,IMPACT=1 - - - - -
} >"$T4F"
run_oracle "$T4F" --out "$TMP/out-jsonfloor" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(reasons_of "$TMP/out-jsonfloor/abc1234/c-jsonfloor/report.json")" == "json: floor Uploaded_variation is not applicable to json and no other format of the case counted it" ]] &&
    ok "a json-only case with a floor the json format cannot count is red: no format checked the floor" || bad "json-only floor: rc=$rc $(reasons_of "$TMP/out-jsonfloor/abc1234/c-jsonfloor/report.json" 2>/dev/null)"
# Extra is the default format's column alone: the tab header does not declare it and the json
# format has no row count for it, so a tab,json case sets the floor aside in both formats.
{
    echo "$HEADER"
    row c-twoaside fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" tab,json "$IMAGE" offline 116.2 - Extra=1,IMPACT=1 - - - - -
} >"$T4F"
run_oracle "$T4F" --out "$TMP/out-twoaside" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(reasons_of "$TMP/out-twoaside/abc1234/c-twoaside/report.json")" == "json: floor Extra is not applicable to json and no other format of the case counted it | tab: floor Extra is not applicable to tab and no other format of the case counted it" ]] &&
    ok "a floor every compared format sets aside is red even when a text format was compared" || bad "two formats set the floor aside: rc=$rc $(reasons_of "$TMP/out-twoaside/abc1234/c-twoaside/report.json" 2>/dev/null)"

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
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - IMPACT=100000 - - - - -
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
    row c-db fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" database 116.2 - - - - - - -
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
    row GRCh38-release116 fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json "$IMAGE" offline 116.2 - - - - - - -
} >"$T7"
run_oracle "$T7" --out "$TMP/out-ref" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >"$TMP/ref.out" 2>"$TMP/ref.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "a golden corpus directory is a stored reference: green against the committed expected outputs" || {
    bad "stored reference exited $rc"
    tail -5 "$TMP/ref.err" | sed 's/^/         /'
}
[[ ! -e "$STUB_DOCKER_CALLED" ]] && ok "--reference-dir runs no container" || bad "docker was called with --reference-dir"
grep -q "Ensembl VEP: stored reference $REPO/tests/golden/116/GRCh38-release116/expected" "$TMP/ref.out" && ok "the log names the stored reference" || bad "stored reference not logged"
python3 - "$TMP/out-ref/p1/GRCh38-release116/report.json" <<'PY' && ok "the report records the stored reference directory, the provenance.json beside it and the image the provenance names" || bad "stored-reference report"
import json, sys
r = json.load(open(sys.argv[1]))
assert r["reference_dir"].endswith("tests/golden/116/GRCh38-release116/expected"), r["reference_dir"]
assert r["reference_provenance"].endswith("tests/golden/116/GRCh38-release116/provenance.json"), r["reference_provenance"]
assert r["perl_image_digest"] == "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da", r["perl_image_digest"]
assert r["reference_provenance_image"] == r["perl_image_digest"]
assert r["refusal_dir"] is None and r["candidate_dir_cache"] is None and r["unordered"] is False
assert r["verdict"] == "green"
PY
{
    echo "$HEADER"
    row GRCh38-release116 fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,warnings "$IMAGE" offline 116.2 - - - - - - -
} >"$T7"
run_oracle "$T7" --out "$TMP/out-ref2" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >/dev/null 2>"$TMP/ref2.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'stored reference lacks default.txt_warnings.txt' "$TMP/ref2.err" && ok "a stored reference without a warnings file fails a warnings case before anything runs" || bad "stored reference without warnings exited $rc"
{
    echo "$HEADER"
    row no-such-corpus fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T7"
run_oracle "$T7" --out "$TMP/out-ref3" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >/dev/null 2>"$TMP/ref3.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "stored reference directory not found: $REPO/tests/golden/116/no-such-corpus/expected" "$TMP/ref3.err" && ok "a case with no stored reference fails closed rather than running Ensembl VEP" || bad "missing stored reference exited $rc"
[[ ! -e "$STUB_DOCKER_CALLED" ]] && ok "no container ran for the failed stored-reference cases either" || bad "docker was called"

echo
echo "STORED REFERENCE PER CASE: the reference cell names a case's own directory"
# The docker stub still fails every call: a case whose reference cell names a directory of
# stored outputs runs no container, with or without --reference-dir.
REF_CELL="tests/golden/116/GRCh38-release116/expected"
TR="$TMP/tr.tsv"
{
    echo "$HEADER"
    row c-stored fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json "$IMAGE" offline 116.2 - - - - - "$REF_CELL" -
} >"$TR"
run_oracle "$TR" --out "$TMP/out-cell" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/cell.out" 2>"$TMP/cell.err"
rc=$?
[[ "$rc" -eq 0 && ! -e "$STUB_DOCKER_CALLED" ]] && ok "a reference cell compares the case against its stored outputs and runs no container, without --reference-dir" || {
    bad "reference cell exited $rc (docker called: $([[ -e "$STUB_DOCKER_CALLED" ]] && echo yes || echo no))"
    tail -5 "$TMP/cell.err" | sed 's/^/         /'
}
grep -q "Ensembl VEP: stored reference $CORPUS/expected" "$TMP/cell.out" && ok "the log names the cell's directory, resolved like the input" || bad "reference cell not logged: $(grep 'stored reference' "$TMP/cell.out")"
python3 - "$TMP/out-cell/p1/c-stored/report.json" "$CORPUS/expected" <<'PY' && ok "the report records the cell's directory, the provenance.json beside it and the image it names" || bad "reference-cell report"
import json, sys
r = json.load(open(sys.argv[1]))
assert r["reference_dir"] == sys.argv[2], r["reference_dir"]
assert r["reference_provenance"].endswith("GRCh38-release116/provenance.json"), r["reference_provenance"]
assert r["perl_image_digest"] == "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da", r["perl_image_digest"]
assert r["verdict"] == "green", r["reasons"]
PY
dry "$TR" "$TMP/d-cell.out"
grep -q "^DRY-RUN c-stored perl stored: $CORPUS/expected/$" "$TMP/d-cell.out" && ! grep -q '^DRY-RUN c-stored perl default' "$TMP/d-cell.out" && ok "a dry run plans the case against its stored directory and no Ensembl VEP run" || bad "dry run of a reference cell: $(grep 'c-stored perl' "$TMP/d-cell.out")"
cc=$(cmd_second 'DRY-RUN c-stored stub vcf' "$TMP/d-cell.out")
[[ "$cc" == *"--format vcf --expected $CORPUS/expected/vcf.vcf.gz --actual "* ]] && ok "the comparison reads the stored file, gzipped as the corpus commits it" || bad "comparison of a reference cell: $cc"
{
    echo "$HEADER"
    row no-such-corpus fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - "$REF_CELL" -
} >"$TR"
run_oracle "$TR" --out "$TMP/out-cell2" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >/dev/null 2>"$TMP/cell2.err"
rc=$?
[[ "$rc" -eq 0 && ! -e "$STUB_DOCKER_CALLED" ]] && ok "beside --reference-dir, a case's reference cell is read before the directory under it" || bad "reference cell beside --reference-dir exited $rc: $(tail -2 "$TMP/cell2.err" | tr '\n' ' ')"
mkdir -p "$DATA/stored/partial"
cp "$REF/default.txt" "$DATA/stored/partial/default.txt"
{
    echo "$HEADER"
    row c-partial fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,vcf "$IMAGE" offline 116.2 - - - - - stored/partial -
} >"$TR"
expect_error "a reference directory lacking a requested format fails preflight by name" "c-partial: stored reference lacks vcf.vcf\[.gz\] under $DATA/stored/partial/" "$TR"
{
    echo "$HEADER"
    row c-plain fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - stored/partial -
} >"$TR"
dry "$TR" "$TMP/d-plain.out"
rc=$?
[[ "$rc" -eq 0 ]] && grep -q "^DRY-RUN c-plain perl stored: $DATA/stored/partial/$" "$TMP/d-plain.out" && ok "a relative reference cell resolves under --data-dir first, and a plain file serves as the .gz does" || bad "relative reference cell: rc=$rc"
{
    echo "$HEADER"
    row c-noref fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - "$DATA/stored/absent" -
} >"$TR"
expect_error "a reference cell naming no directory is an error" "c-noref: stored reference directory not found: $DATA/stored/absent" "$TR"
[[ "$(grep -c 'c-noref' "$TMP/e.out.err")" -eq 1 ]] && ok "a missing reference directory is reported once, without a line per format" || bad "missing reference directory reported $(grep -c 'c-noref' "$TMP/e.out.err") times"
{
    echo "$HEADER"
    row c-badcell fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - yes - -
} >"$TR"
expect_error "an unordered cell other than - or 1 is an error" "c-badcell: unordered 'yes' is neither - nor 1" "$TR"
[[ ! -e "$STUB_DOCKER_CALLED" ]] && ok "no container ran in the reference-cell checks" || bad "docker was called in the reference-cell checks"
cp "$TMP/docker-answering" "$STUBDIR/docker"
{
    echo "$HEADER"
    row c-stored fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - "$REF_CELL" -
    row c-live fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - - -
} >"$TR"
: >"$STUB_DOCKER_LOG"
run_oracle "$TR" --out "$TMP/out-mixed" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/mixed.out" 2>&1
rc=$?
n=$(count_lines '^run ' "$STUB_DOCKER_LOG")
[[ "$rc" -eq 0 && "$n" -eq 2 && "$(cut -f1,3 "$TMP/out-mixed/status.tsv" | tail -n +2 | tr '\n\t' '  ')" == "c-stored green c-live green " ]] && ok "a table mixing a stored case and a live one runs Ensembl VEP for the live case alone (2 docker run calls for its two formats)" || bad "mixed table: rc=$rc docker runs=$n $(tr '\n' ' ' <"$TMP/out-mixed/status.tsv")"
grep -q 'live=1 ' "$TMP/mixed.out" && ok "the run log counts the cases that run Ensembl VEP live" || bad "live count missing from the log"

echo
echo "UNORDERED AND LAYOUT: the comparator's --unordered and a cache root in Ensembl VEP's layout"
mkdir -p "$DATA/cache_root/homo_sapiens/116_GRCh38/transcripts" "$DATA/cache_root/homo_sapiens_refseq/116_GRCh38/transcripts"
TU="$TMP/tu.tsv"
{
    echo "$HEADER"
    row c-unordered fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - --pick_allele default,json "$IMAGE" offline 116.2 - - - - 1 - -
    row c-layout fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - - default "$IMAGE" offline 116.2 - - - - - - -
    row c-layout-near fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - "--nearest symbol" default "$IMAGE" offline 116.2 - - - - - - -
    row c-layout-refseq fixture "$CORPUS_INPUT" "116_GRCh38=layout:$DATA/cache_root" - --refseq default "$IMAGE" offline 116.2 - - - - - "$REF_CELL" -
} >"$TU"
dry "$TU" "$TMP/du.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "a table with an unordered case and layout caches plans" || {
    bad "unordered and layout dry run exited $rc"
    head -5 "$TMP/du.out.err" | sed 's/^/         /'
}
uc=$(cmd_second 'DRY-RUN c-unordered stub default' "$TMP/du.out")
uj=$(cmd_second 'DRY-RUN c-unordered stub json' "$TMP/du.out")
[[ "$uc" == *" --unordered"* && "$uj" == *" --unordered"* ]] && ok "an unordered cell of 1 passes --unordered to the comparator in every format" || bad "unordered: default=$uc json=$uj"
tc=$(cmd_second 'DRY-RUN c-four stub tab' "$TMP/d1.out")
[[ "$tc" != *"--unordered"* ]] && ok "a - unordered cell passes no --unordered" || bad "unordered leaked: $tc"
lc=$(cmd_after 'DRY-RUN c-layout stub default' "$TMP/du.out")
[[ "$lc" == "$STUBDIR/vep -i $CORPUS/variants.vcf -o $TMP/never/work/cases/c-layout/candidates/stub/default.txt --offline --dir_cache $DATA/cache_root --species homo_sapiens --cache_version 116 --assembly GRCh38 --format vcf --no_stats --quiet --force_overwrite" ]] &&
    ok "a layout: cache cell gives the candidate --dir_cache <root> --species --cache_version --assembly and no --json_cache" || bad "layout candidate: $lc"
lp=$(cmd_after 'DRY-RUN c-layout perl default' "$TMP/du.out")
[[ "$lp" == *"--offline --cache --dir_cache /work/cache --cache_version 116 --species homo_sapiens --assembly GRCh38 "* ]] && ok "the Ensembl VEP command of a layout case is unchanged" || bad "layout perl: $lp"
ln_=$(cmd_second 'DRY-RUN c-layout-near stub default' "$TMP/du.out")
[[ "$ln_" == *"--nearest-sites $DATA/cache_root/homo_sapiens/116_GRCh38" ]] && ok "a --nearest layout case passes the leaf the candidate resolves as --nearest-sites" || bad "layout nearest sites: $ln_"
lr=$(cmd_after 'DRY-RUN c-layout-refseq stub default' "$TMP/du.out")
[[ "$lr" == *"--dir_cache $DATA/cache_root --species homo_sapiens --cache_version 116 --assembly GRCh38 "*" --refseq" ]] && ok "an absolute layout root is taken as given, and --refseq stays a flag the candidate resolves (the refseq leaf is checked)" || bad "layout refseq: $lr"
{
    echo "$HEADER"
    row c-layout-merged fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - --merged default "$IMAGE" offline 116.2 - - - - - - -
} >"$TU"
expect_error "a layout root lacking the leaf the flags select fails preflight by its layout path" "c-layout-merged: cache root $DATA/cache_root has no homo_sapiens_merged/116_GRCh38/transcripts/" "$TU"
{
    echo "$HEADER"
    row c-layout fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - -
} >"$TU"
run_oracle "$TU" --out "$TMP/out-layout" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
python3 - "$TMP/out-layout/p1/c-layout/report.json" "$DATA/cache_root" "$rc" <<'PY' && ok "the report records the layout root as candidate_dir_cache and the resolved leaf as candidate_json_cache" || bad "layout report (rc=$rc)"
import json, sys
r = json.load(open(sys.argv[1]))
assert sys.argv[3] == "0" and r["verdict"] == "green", r["reasons"]
assert r["candidate_dir_cache"] == sys.argv[2], r["candidate_dir_cache"]
assert r["candidate_json_cache"] == sys.argv[2] + "/homo_sapiens/116_GRCh38", r["candidate_json_cache"]
PY

echo
echo "REFUSAL: a case both engines are expected to refuse"
mkdir -p "$DATA/refusals/mixed" "$DATA/refusals/incomplete" "$DATA/refusals/blank"
printf 'ERROR: [plugin] LOEUF: parameters mix key=value and positional styles\n' >"$DATA/refusals/mixed/perl.txt"
printf 'failing on request\n' >"$DATA/refusals/mixed/candidate.txt"
cp "$DATA/refusals/mixed/perl.txt" "$DATA/refusals/incomplete/perl.txt"
cp "$DATA/refusals/mixed/perl.txt" "$DATA/refusals/blank/perl.txt"
printf '\n' >"$DATA/refusals/blank/candidate.txt"
TF="$TMP/tf.tsv"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv.gz,gene" default,tab "$IMAGE" offline 116.2 - - - - - - refusals/mixed
} >"$TF"
PERL_TEXT="ERROR: [plugin] LOEUF: parameters mix key=value and positional styles"
: >"$STUB_DOCKER_LOG"
STUB_PERL_RC=2 STUB_PERL_STDERR="$PERL_TEXT (Plugin.pm)" STUB_CAND_RC=1 run_oracle "$TF" --out "$TMP/out-refuse" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/refuse.out" 2>"$TMP/refuse.err"
rc=$?
[[ "$rc" -eq 0 && "$(sed -n 2p "$TMP/out-refuse/status.tsv" | cut -f3)" == "green" ]] && ok "both engines refusing every format with the expected texts is green, exit 0" || {
    bad "refusal green path: rc=$rc $(reasons_of "$TMP/out-refuse/p1/c-refuse/report.json" 2>/dev/null)"
    tail -3 "$TMP/refuse.err" | sed 's/^/         /'
}
! grep -q 'WARN: \[feature_oracle\] c-refuse: Ensembl VEP' "$TMP/refuse.err" && ok "an expected Ensembl VEP refusal raises no warning" || bad "a refusal was warned about"
n=$(count_lines '^run ' "$STUB_DOCKER_LOG")
[[ "$n" -eq 2 ]] && ok "Ensembl VEP ran once per requested format of the refusal case" || bad "$n docker run calls for the refusal case, expected 2"
python3 - "$TMP/out-refuse/p1/c-refuse/report.json" "$DATA/refusals/mixed" <<'PY' && ok "the report records the refusal directory and, per format, each engine's exit status and text match, with no comparator report" || bad "refusal report"
import json, os, sys
r = json.load(open(sys.argv[1]))
assert r["refusal_dir"] == sys.argv[2], r["refusal_dir"]
assert sorted(r["formats"]) == ["default", "tab"], sorted(r["formats"])
for entry in r["formats"].values():
    assert entry["verdict"] == "green" and entry["reasons"] == [], entry
    assert entry["refusal"] == {"perl": {"exit": 2, "text_matched": True}, "candidate": {"exit": 1, "text_matched": True}}, entry["refusal"]
    assert "report" not in entry
assert not os.path.exists(os.path.join(os.path.dirname(sys.argv[1]), "default.json"))
PY
STUB_CAND_RC=1 run_oracle "$TF" --out "$TMP/out-refuse2" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-refuse2/p1/c-refuse/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: Ensembl VEP exited 0 where a refusal was expected | tab: Ensembl VEP exited 0 where a refusal was expected" ]] && ok "an Ensembl VEP run that exits 0 is red in every format of a refusal case" || bad "refusal with Ensembl VEP exit 0: rc=$rc reasons=$reasons"
STUB_PERL_RC=2 STUB_PERL_STDERR="$PERL_TEXT" STUB_CAND_RC=1 STUB_CAND_STDERR="ERROR: [cli] something else" run_oracle "$TF" --out "$TMP/out-refuse3" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-refuse3/p1/c-refuse/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: candidate exited 1 without the expected refusal text: ERROR: [cli] something else"* && "$reasons" != *"Ensembl VEP"* ]] && ok "a candidate refusing with another text is red with its stderr; the matching Ensembl VEP side stays green" || bad "refusal text mismatch: rc=$rc reasons=$reasons"
STUB_PERL_RC=2 STUB_PERL_STDERR="$PERL_TEXT" run_oracle "$TF" --out "$TMP/out-refuse4" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-refuse4/p1/c-refuse/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: candidate exited 0 where a refusal was expected"* ]] && ok "a candidate that runs the case is red" || bad "candidate exit 0 in a refusal case: rc=$rc reasons=$reasons"
dry "$TF" "$TMP/df.out"
rc=$?
[[ "$rc" -eq 0 ]] && grep -q "^DRY-RUN c-refuse refusal: $DATA/refusals/mixed (perl.txt, candidate.txt)$" "$TMP/df.out" && [[ "$(count_lines '^DRY-RUN c-refuse perl ' "$TMP/df.out")" -eq 2 && "$(count_lines '^DRY-RUN c-refuse stub ' "$TMP/df.out")" -eq 2 ]] && ! grep -q compare_vep_formats "$TMP/df.out" &&
    ok "a dry run names the refusal directory and plans both engines' runs and no comparison" || bad "refusal dry run: rc=$rc"
# Under --reference-dir a refusal case still runs both engines live: the stored outputs of
# other cases are not a refusal.
: >"$STUB_DOCKER_LOG"
STUB_PERL_RC=2 STUB_PERL_STDERR="$PERL_TEXT" STUB_CAND_RC=1 run_oracle "$TF" --out "$TMP/out-refuse5" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >"$TMP/refuse5.out" 2>&1
rc=$?
n=$(count_lines '^run ' "$STUB_DOCKER_LOG")
[[ "$rc" -eq 0 && "$n" -eq 2 && "$(sed -n 2p "$TMP/out-refuse5/status.tsv" | cut -f3)" == "green" ]] && grep -q 'live=1 ' "$TMP/refuse5.out" && ok "a refusal case under --reference-dir runs both engines live (2 docker runs for its two formats) and is green on matching refusals" || bad "refusal under --reference-dir: rc=$rc docker runs=$n"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - refusals/incomplete
} >"$TF"
expect_error "a refusal directory lacking candidate.txt is an error" "c-refuse: refusal directory lacks candidate.txt: $DATA/refusals/incomplete" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - refusals/blank
} >"$TF"
expect_error "a blank refusal text is an error (it would match every stderr)" "c-refuse: refusal text $DATA/refusals/blank/candidate.txt is empty" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - "$DATA/refusals/absent"
} >"$TF"
expect_error "a refusal cell naming no directory is an error" "c-refuse: refusal directory not found: $DATA/refusals/absent" "$TF"
[[ "$(grep -c 'c-refuse' "$TMP/e.out.err")" -eq 1 ]] && ok "a missing refusal directory is reported once, without a line per text file" || bad "missing refusal directory reported $(grep -c 'c-refuse' "$TMP/e.out.err") times"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - "$REF_CELL" refusals/mixed
} >"$TF"
expect_error "refusal beside reference is an error" "c-refuse: refusal and reference are both set" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - IMPACT=1 - - - - refusals/mixed
} >"$TF"
expect_error "refusal beside a floor is an error" "c-refuse: refusal with floors 'IMPACT=1'" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default,warnings "$IMAGE" offline 116.2 - - - - - - refusals/mixed
} >"$TF"
expect_error "refusal beside warnings is an error" "c-refuse: refusal with warnings in formats" "$TF"

echo
echo "OPTIONS: --keep-outputs, --jobs and a .gz FASTA on the real path"
run_oracle "$T4" --out "$TMP/out-keep" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs --work-dir "$TMP/work-keep" >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && -f "$TMP/work-keep/cases/c-four/perl/default.txt" && -f "$TMP/work-keep/cases/c-four/candidates/p1/json.jsonl" ]] && ok "--keep-outputs keeps both engines' outputs under --work-dir" || bad "--keep-outputs: rc=$rc"
T8="$TMP/t8.tsv"
{
    echo "$HEADER"
    row c-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - - -
    row c-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" vcf,json "$IMAGE" offline 116.2 - - - - - - -
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
    row c-near-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest symbol" default "$IMAGE" offline 116.2 - - - - - - -
    row c-near-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest gene" default "$IMAGE" offline 116.2 - - - - - - -
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
# Preflight's TSS_SLICES and run_perl's lock read one predicate, so they agree on a case
# whatever its spelling: two cases spelled -nearest and --NEAREST on a slice without the tree
# file take the lock, and without flock on PATH preflight refuses --jobs 2 for the same table.
T_NEAR_SP="$TMP/t-near-sp.tsv"
{
    echo "$HEADER"
    row c-near-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-nearest symbol" default "$IMAGE" offline 116.2 - - - - - - -
    row c-near-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--NEAREST gene" default "$IMAGE" offline 116.2 - - - - - - -
} >"$T_NEAR_SP"
# noflock_path: PATH with every directory that carries flock replaced by a directory of links
# to its other entries, so the runner's `command -v flock` finds none.
noflock_path() {
    local dir out="" farm n=0
    local IFS=':'
    for dir in $PATH; do
        if [[ -x "$dir/flock" ]]; then
            farm="$TMP/noflock/$((n += 1))"
            mkdir -p "$farm"
            ln -s "$dir"/* "$farm/" 2>/dev/null
            rm -f "$farm/flock"
            dir="$farm"
        fi
        out="$out${out:+:}$dir"
    done
    printf '%s' "$out"
}
if command -v flock >/dev/null; then
    STUB_PERL_SLEEP=1.5 STUB_PERL_TIMES="$TMP/near-sp.times" run_oracle "$T_NEAR_SP" --out "$TMP/out-near-sp" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs 2 --keep-outputs >/dev/null 2>&1 || true
    overlap=$(sort -n "$TMP/near-sp.times" | awk 'NR > 1 && $1 < prev_end {o = 1} {prev_end = $2} END {print o + 0}')
    [[ -e "$TMP/out-near-sp/work/locks/homo_sapiens_116_GRCh38.tss.lock" && "$(wc -l <"$TMP/near-sp.times")" -eq 2 && "$overlap" -eq 0 ]] && ok "run_perl locks the -nearest and --NEAREST cases of a slice without the tree file, and their runs do not overlap" || bad "spelled nearest cases: lock $([[ -e "$TMP/out-near-sp/work/locks/homo_sapiens_116_GRCh38.tss.lock" ]] && echo taken || echo missing), runs $(tr '\n' ' ' <"$TMP/near-sp.times")"
fi
PATH="$(noflock_path)" run_oracle "$T_NEAR_SP" --out "$TMP/out-near-sp2" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs 2 >/dev/null 2>"$TMP/near-sp.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$TMP/out-near-sp2" ]] && grep -q 'flock is required to run the --nearest cases of .*116_GRCh38 under --jobs 2' "$TMP/near-sp.err" && ok "preflight's TSS_SLICES and run_perl's lock agree on the same case: without flock, the -nearest and --NEAREST cases refuse --jobs 2 before anything runs" || bad "spelled nearest cases without flock: rc=$rc $(tr '\n' ' ' <"$TMP/near-sp.err")"
run_oracle "$T2" --out "$TMP/out-fa" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs >/dev/null 2>"$TMP/fa.err"
rc=$?
[[ -f "$TMP/out-fa/work/fasta/1/reference.fa" && -f "$TMP/out-fa/work/fasta/1/reference.fa.fai" && "$(head -c 1 "$TMP/out-fa/work/fasta/1/reference.fa")" == ">" ]] &&
    ok "a .gz FASTA is decompressed into <work>/fasta/<n>/ with its .fai beside it" || bad "FASTA decompression (rc=$rc)"

echo
echo "MORE FAIL-CLOSED PATHS: cells, names and option values the first sections did not reach"
{
    echo "$HEADER"
    row c-cache fixture "$CORPUS_INPUT" "GRCh38=tests/golden/116/GRCh38-release116/json_cache" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a cache cell that is not <version>_<assembly> is an error" "cache 'GRCh38' is not <version>_<assembly>" "$T3"
{
    echo "$HEADER"
    row c-img fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default - offline 116.2 - - - - - - -
} >"$T3"
expect_error "a - perl_image is an error" "perl_image is -" "$T3"
{
    echo "$HEADER"
    row c-rel fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline - - - - - - - -
} >"$T3"
expect_error "a - reference_release is an error" "reference_release is -" "$T3"
mkdir -p "$DATA/other"
cp "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.gz" "$DATA/other/reference.fa.gz"
cp "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.fai" "$DATA/other/reference.fa.fai"
{
    echo "$HEADER"
    row c-fa1 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - -
    row c-fa2 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache other/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-fa2.out"
rc=$?
fa1=$(cmd_after 'DRY-RUN c-fa1 stub default' "$TMP/d-fa2.out")
fa2=$(cmd_after 'DRY-RUN c-fa2 stub default' "$TMP/d-fa2.out")
[[ "$rc" -eq 0 && "$fa1" == *"--fasta $TMP/never/work/fasta/1/reference.fa" && "$fa2" == *"--fasta $TMP/never/work/fasta/2/reference.fa" ]] &&
    ok "two .gz FASTAs that decompress to one file name each get a directory of their own" || bad "two FASTAs of one name: rc=$rc fa1=$fa1 fa2=$fa2"
{
    echo "$HEADER"
    row c-fa1 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - -
    row c-fa1b fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-fa3.out"
fa1=$(cmd_after 'DRY-RUN c-fa1 stub default' "$TMP/d-fa3.out")
fa1b=$(cmd_after 'DRY-RUN c-fa1b stub default' "$TMP/d-fa3.out")
[[ "$fa1" == *"--fasta $TMP/never/work/fasta/1/reference.fa" && "$fa1b" == *"--fasta $TMP/never/work/fasta/1/reference.fa" ]] &&
    ok "two cases on one .gz FASTA share its one decompressed copy" || bad "shared FASTA: fa1=$fa1 fa1b=$fa1b"
{
    echo "$HEADER"
    row c-refseq fixture "$CORPUS_INPUT" 116_GRCh38 - --refseq default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "--refseq in the flags selects the homo_sapiens_refseq slice" "Ensembl VEP cache slice not found: $DATA/caches/perl/vep-cache/homo_sapiens_refseq/116_GRCh38" "$T3"
dry "$T1" "$TMP/d8.out" --only GRCh38-release116
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "names 'GRCh38-release116', which is not a case of tier fixture" "$TMP/d8.out.err" && ok "--only naming a case of the other tier is an error, not a silent drop" || bad "--only across tiers exited $rc"
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_option_error "a pin named like the work directory is refused" "candidate pin 'work' names a file or directory the run writes under --out" --candidates "$STUBDIR/vep" --candidate-pins work
expect_option_error "a pin named status.tsv is refused" "candidate pin 'status.tsv' names a file or directory" --candidates "$STUBDIR/vep" --candidate-pins status.tsv
expect_option_error "a pin named run.json is refused" "candidate pin 'run.json' names a file or directory" --candidates "$STUBDIR/vep" --candidate-pins run.json
PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$T3" --data-dir "$DATA" --out "$TMP/never" --work-dir "$TMP/elsewhere" --candidates "$STUBDIR/vep" --candidate-pins work --dry-run >"$TMP/e.out" 2>"$TMP/e.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "the pin work is free when --work-dir lies outside --out" || bad "pin work with an outside work dir exited $rc"
dry "$T1" "$TMP/d9.out" --tier population
grep -q 'WARN: \[feature_oracle\] population tier: Ensembl VEP runs live on 1 case(s) without a stored reference' "$TMP/d9.out.err" && ok "a population run without a stored reference warns that Ensembl VEP runs live" || bad "no population warning"
dry "$T1" "$TMP/d10.out" --tier population --reference-dir "$REPO/tests/golden/116"
rc=$?
if [[ "$rc" -eq 0 ]] && grep -q '^DRY-RUN GRCh38-release116 perl stored' "$TMP/d10.out" && ! grep -q 'population tier:' "$TMP/d10.out.err"; then
    ok "the population warning is silent with --reference-dir (the case was planned against its stored reference)"
else
    bad "population with --reference-dir: rc=$rc"
fi
{
    echo "$HEADER"
    row c-glob fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--filter *" default "$IMAGE" offline 116.2 - - - - - - -
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
    row fasta fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - -
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
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - - -
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
    row c-db37 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache - - default "$IMAGE" database 116.2 - - - - - - -
} >"$T6"
run_oracle "$T6" --out "$TMP/out-db37" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
server=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["perl_server"])' "$TMP/out-db37/p1/c-db37/report.json")
[[ "$server" == "ensembldb.ensembl.org:3337" ]] && ok "a GRCh37 database case records ensembldb.ensembl.org:3337" || bad "GRCh37 perl_server: $server"
# An interrupted run stops its engines at once, under --jobs 1 (the engine is a waited
# job of the main shell) and under --jobs 2 (the engines are children of case
# subshells), and leaves no status table. The elapsed time is bounded too: the stub engine
# sleeps far longer than the bound, so a run that only died after the sleep exceeds it.
interrupt_run() { # interrupt_run <jobs> [table] [expected engines]: ok lines for one --jobs value
    local jobs="$1" table="${2:-$T8}" expected_engines="${3:-$1}" runner started elapsed alive=0 spid n
    rm -rf "$TMP/out-int" "$TMP/sleepers"
    mkdir -p "$TMP/out-int"
    echo "stale" >"$TMP/out-int/status.tsv"
    STUB_PERL_SLEEP=120 STUB_PERL_PIDFILE="$TMP/sleepers" PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$table" --data-dir "$DATA" --out "$TMP/out-int" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs "$jobs" >"$TMP/int.out" 2>&1 &
    runner=$!
    # A loaded host can take tens of seconds to start the engines; the stubs sleep 120 s so
    # every one is still running when the signal lands.
    for _ in $(seq 1 60); do
        [[ -s "$TMP/sleepers" && "$(wc -l <"$TMP/sleepers" | tr -d ' ')" -ge "$expected_engines" ]] && break
        sleep 1
    done
    # The engines running when the signal lands; a lock holder's death lets a waiting
    # case's engine start for an instant before it too is stopped, so the count precedes the kill.
    n=$(wc -l <"$TMP/sleepers" 2>/dev/null | tr -d ' '); n=${n:-0}
    started=$SECONDS
    kill -TERM "$runner" 2>/dev/null
    wait "$runner" 2>/dev/null
    elapsed=$((SECONDS - started))
    sleep 1
    while read -r spid; do kill -0 "$spid" 2>/dev/null && alive=$((alive + 1)); done <"$TMP/sleepers"
    [[ "$n" -eq "$expected_engines" && "$alive" -eq 0 && "$elapsed" -lt 15 ]] && ok "SIGTERM to a --jobs $jobs run stops its $n running stub engine(s) at once (${elapsed}s after the signal, not after the 120 s sleep)" || bad "--jobs $jobs: $n engine(s) running at the signal (expected $expected_engines), $alive alive after it, exited after ${elapsed}s"
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
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
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
    row c-masked fixture "$TMP/masked/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 breakend_mate_local_read - - - - - -
} >"$T3"
STUB_REF_DIR="$TMP/masked" run_oracle "$T3" --out "$TMP/out-masked" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/masked.out" 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-masked/p1/c-masked/report.json")
skipped=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["counts"]["documented_divergences_skipped"])' "$TMP/out-masked/p1/c-masked/default.json")
[[ "$n" -ge 1 && "$skipped" -ge 1 && "$rc" -eq 1 && "$reasons" == *"documented divergences not observed"* && "$reasons" != *"no entry compared"* ]] &&
    ok "a case whose every entry is a documented divergence is judged by the comparator, not read as two empty outputs ($skipped skipped)" || bad "masked-only case: rc=$rc skipped=$skipped reasons=$reasons"
# A requested mask must set something aside, and a mask with an expected count must apply to
# that many entries in every output format: the stub candidate gains splice_region_variant on
# one intron_variant row of its default output, the registry pair class splice_family_swap's
# shape; the tab copy is the reference's, so the class applies there to nothing.
{
    echo "$HEADER"
    row c-swap fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 splice_family_swap=1 - - - - - -
    row c-swapcount fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 splice_family_swap=2 - - - - - -
    row c-swapidle fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 splice_family_swap - - - - - -
    row c-swapboth fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 splice_family_swap=1 - - - - - -
} >"$T3"
STUB_CAND_MUTATE=add_splice_region run_oracle "$T3" --out "$TMP/out-swap" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-swap >"$TMP/swap.out" 2>&1
rc=$?
applied=$(python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); print(r["formats"]["default"]["masks_applied"].get("splice_family_swap"), r["masks"])' "$TMP/out-swap/p1/c-swap/report.json" 2>/dev/null)
[[ "$rc" -eq 0 && "$(sed -n 2p "$TMP/out-swap/status.tsv" | cut -f3)" == "green" && "$applied" == "1 splice_family_swap=1" ]] &&
    ok "a registry pair class with its expected count is green when the comparator applies it that often; the report carries the count per format" || bad "counted mask green path: rc=$rc applied=$applied $(reasons_of "$TMP/out-swap/p1/c-swap/report.json" 2>/dev/null)"
STUB_CAND_MUTATE=add_splice_region run_oracle "$T3" --out "$TMP/out-swapcount" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-swapcount >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-swapcount/p1/c-swapcount/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: mask splice_family_swap applied to 1 entries, 2 expected" ]] &&
    ok "a mask applied to another number of entries than its count is red, the format and both numbers named" || bad "counted mask red path: rc=$rc reasons=$reasons"
run_oracle "$T3" --out "$TMP/out-swapidle" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-swapidle >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-swapidle/p1/c-swapidle/report.json")
[[ "$rc" -eq 1 && "$reasons" == "mask splice_family_swap is requested and applied in no compared format" ]] &&
    ok "a requested mask the comparator applies in no format is red" || bad "idle mask: rc=$rc reasons=$reasons"
STUB_CAND_MUTATE=add_splice_region run_oracle "$T3" --out "$TMP/out-swapboth" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-swapboth >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-swapboth/p1/c-swapboth/report.json")
[[ "$rc" -eq 1 && "$reasons" == "tab: mask splice_family_swap applied to 0 entries, 1 expected" ]] &&
    ok "an expected count is held in every output format compared, a format where the class applied nowhere included" || bad "count per format: rc=$rc reasons=$reasons"
# A class the comparator's report names among those asked for and leaves out of its unused
# list is used though nothing is counted under it. The corpus is one frameshift row whose
# DownstreamProtein the manifest documents under plugin_hgvs_shifted_cds_coordinates with
# downstream_start_past_own_codon as its second class: the comparator sets the row aside with
# both named, counts it under the first alone and lists neither unused.
NEW_CLASS=plugin_hgvs_shifted_cds_coordinates
START_CLASS=downstream_start_past_own_codon
END_CLASS=downstream_utr_frame_past_coding_end
two_cause_entry() { # two_cause_entry <second class or empty>: the manifest of the one-row corpus
    python3 -c 'import json, sys
entry = {"location": "21:6006", "allele": "-", "feature": "ENST00000000205", "feature_type": "Transcript",
         "consequence_set": "frameshift_variant", "field": "DownstreamProtein", "vep": "VK", "vep_rs": "PKVK",
         "expected_divergence": sys.argv[1]}
if sys.argv[2]:
    entry["second_expected_divergence"] = sys.argv[2]
print(json.dumps({"field_divergences": [entry]}))' "$NEW_CLASS" "$1"
}
for corpus in needed plain; do
    mkdir -p "$TMP/$corpus"
    {
        grep '^#' "$REF/default.txt"
        printf 'v\t21:6006\t-\t-\tENST00000000205\tTranscript\tframeshift_variant\t-\t6\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVS_OFFSET=3;DownstreamProtein=VK;ProteinLengthChange=1\n'
    } >"$TMP/$corpus/default.txt"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t6005\tv\tCA\tC\t.\t.\t.\n' >"$TMP/$corpus/variants.vcf"
done
two_cause_entry "$START_CLASS" >"$TMP/needed/manifest.json"
two_cause_entry "" >"$TMP/plain/manifest.json"
{
    echo "$HEADER"
    row c-needed fixture "$TMP/needed/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS" - - - - - -
    row c-needed-count fixture "$TMP/needed/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS=1" - - - - - -
    row c-needed-idle fixture "$TMP/needed/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS,$END_CLASS" - - - - - -
    row c-plain fixture "$TMP/plain/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS" - - - - - -
} >"$T3"
needed_oracle() { # needed_oracle <case> <corpus>: the one case over the one-row corpus, its candidate's protein the own span's
    STUB_REF_DIR="$TMP/$2" STUB_CAND_MUTATE=own_span_protein run_oracle "$T3" --out "$TMP/out-$1" --candidates "$STUBDIR/vep" --candidate-pins p1 --only "$1" >"$TMP/$1.out" 2>&1
}
needed_oracle c-needed needed
rc=$?
masks_read=$(python3 -c 'import json,sys
r = json.load(open(sys.argv[1] + "/report.json")); m = json.load(open(sys.argv[1] + "/default.json"))["masks"]
print(r["verdict"], r["formats"]["default"]["masks_applied"], m["requested"], m["unused"])' "$TMP/out-c-needed/p1/c-needed" 2>/dev/null)
[[ "$rc" -eq 0 && "$masks_read" == "green {'$NEW_CLASS': 1} ['$START_CLASS', '$NEW_CLASS'] []" ]] &&
    ok "a class the report names and leaves out of its unused list is used: the case is green with nothing counted under it" || bad "needed class: rc=$rc read=$masks_read $(reasons_of "$TMP/out-c-needed/p1/c-needed/report.json" 2>/dev/null)"
needed_oracle c-needed-count needed
rc=$?
reasons=$(reasons_of "$TMP/out-c-needed-count/p1/c-needed-count/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: mask $START_CLASS applied to 0 entries, 1 expected" ]] &&
    ok "a count expected under such a class is red: the value is counted under its own class alone" || bad "needed class with a count: rc=$rc reasons=$reasons"
needed_oracle c-needed-idle needed
rc=$?
reasons=$(reasons_of "$TMP/out-c-needed-idle/p1/c-needed-idle/report.json")
[[ "$rc" -eq 1 && "$reasons" == "mask $END_CLASS is requested and applied in no compared format" ]] &&
    ok "a class beside it that nothing needs and nothing is counted under is red, and it alone" || bad "idle class beside a needed one: rc=$rc reasons=$reasons"
needed_oracle c-plain plain
rc=$?
reasons=$(reasons_of "$TMP/out-c-plain/p1/c-plain/report.json")
[[ "$rc" -eq 1 && "$reasons" == "mask $START_CLASS is requested and applied in no compared format" ]] &&
    ok "the same class is red where the row's entry names no second class: the report lists it unused" || bad "class no entry needs: rc=$rc reasons=$reasons"
# The same comparison through an interpreter that, after the comparator, lists every class of
# the report without a count as unused: the oracle reads the list, so the case is red.
cat >"$TMP/listing-python" <<'EOF'
#!/usr/bin/env bash
python3 "$@"
rc=$?
[[ "${1:-}" == */compare_vep_formats.py ]] || exit "$rc"
report=""
while [[ $# -gt 0 ]]; do
    [[ "$1" == "--report" ]] && report="${2:-}"
    shift
done
[[ -f "$report" ]] && python3 -c 'import json, sys
report = json.load(open(sys.argv[1]))
masks = report["masks"]
masks["unused"] = sorted(klass for klass in masks["requested"] if not masks["applied"].get(klass))
json.dump(report, open(sys.argv[1], "w"))' "$report"
exit "$rc"
EOF
chmod +x "$TMP/listing-python"
PYTHON="$TMP/listing-python" STUB_REF_DIR="$TMP/needed" STUB_CAND_MUTATE=own_span_protein run_oracle "$T3" --out "$TMP/out-listed" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-needed >"$TMP/listed.out" 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-listed/p1/c-needed/report.json")
[[ "$rc" -eq 1 && "$reasons" == "mask $START_CLASS is requested and applied in no compared format" ]] &&
    ok "the same comparison is red where its report lists the class unused" || bad "class listed unused: rc=$rc reasons=$reasons"
mkdir -p "$TMP/refbad/c-badprov/expected"
cp "$REF/default.txt" "$TMP/refbad/c-badprov/expected/default.txt"
echo '{not json' >"$TMP/refbad/c-badprov/provenance.json"
{
    echo "$HEADER"
    row c-badprov fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - -
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
    row c-miss fixture "$DATA/absent.vcf" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-pre" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
[[ "$rc" -ne 0 && ! -e "$TMP/out-pre/status.tsv" ]] && ok "a run that fails preflight removes the previous run's status.tsv" || bad "preflight failure left status.tsv (rc=$rc)"
mkdir -p "$TMP/out-dry"
echo "stale" >"$TMP/out-dry/status.tsv"
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-dry" --candidates "$STUBDIR/vep" --candidate-pins p1 --dry-run >/dev/null 2>&1
[[ -e "$TMP/out-dry/status.tsv" ]] && ok "a dry run leaves an existing status.tsv alone" || bad "a dry run removed status.tsv"

echo
echo "FORMAT JOBS: --format-jobs runs a case's formats at once and reports as the serial run does"
expect_option_error "--format-jobs 0 is refused" "--fork, --perl-fork, --jobs and --format-jobs take a positive integer (got '0')" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 0
expect_option_error "a non-numeric --format-jobs is refused" "take a positive integer (got 'two')" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs two
FJ="$TMP/fj"
mkdir -p "$FJ/serial" "$FJ/serial-rc" "$FJ/serial-refuse"
# same_report.py <serial report.json> <parallel report.json>: exit 0 when the two records agree
# byte for byte once the three clock fields are dropped (write_report dumps sorted keys at indent 1).
cat >"$FJ/same_report.py" <<'PY'
import json, sys
records = [json.load(open(p)) for p in sys.argv[1:3]]
for record in records:
    for clock in ("run_utc", "started_utc", "finished_utc"):
        del record[clock]
dumps = [json.dumps(r, indent=1, sort_keys=True) for r in records]
sys.exit(0 if dumps[0] == dumps[1] else 1)
PY
# Both runs write the same --out, so the per-format reports name the same paths; the serial
# run's reports are copied aside before the parallel run overwrites them.
run_oracle "$T4" --out "$FJ/out" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$FJ/serial.out" 2>"$FJ/serial.err"
rc1=$?
cp -R "$FJ/out/p1/c-four/." "$FJ/serial/"
run_oracle "$T4" --out "$FJ/out" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 3 --keep-outputs >"$FJ/parallel.out" 2>"$FJ/parallel.err"
rc2=$?
[[ "$rc1" -eq 0 && "$rc2" -eq 0 && "$(sed -n 2p "$FJ/out/status.tsv" | cut -f1-3)" == "$(row c-four p1 green)" ]] && ok "--format-jobs 3 over default,tab,vcf,json,warnings exits 0 with a green status row" || bad "--format-jobs 3: rc=$rc1/$rc2 $(sed -n 2p "$FJ/out/status.tsv")"
same=1
for f in default tab vcf json warnings; do
    cmp -s "$FJ/serial/$f.json" "$FJ/out/p1/c-four/$f.json" || same=0
done
[[ "$same" -eq 1 ]] && ok "every per-format comparator report is byte-identical to the serial run's" || bad "a per-format report differs between the serial and the parallel run"
python3 "$FJ/same_report.py" "$FJ/serial/report.json" "$FJ/out/p1/c-four/report.json" && ok "report.json is byte-identical to the serial run's once the three clock fields are dropped" || bad "report.json differs from the serial run's beyond the clock fields"
order=$(grep -oE '\[c-four\] p1 (default|tab|vcf|json)$' "$FJ/parallel.out" | awk '{print $3}' | paste -sd, -)
[[ "$order" == "default,tab,vcf,json" ]] && ok "the per-format start lines reach the run's log in format order as the jobs are dispatched" || bad "log order: $order"
cand="$FJ/out/work/cases/c-four/candidates/p1"
if [[ -f "$cand/default.job.out" && -f "$cand/default.job.err" && -f "$cand/json.job.out" && -f "$cand/json.job.err" && ! -s "$cand/default.job.out" && ! -s "$cand/default.job.err" ]]; then
    ok "each format job keeps its own log pair beside the candidate's outputs (<fmt>.job.out, <fmt>.job.err), empty for a quiet job"
else
    bad "format job log pair missing or wrong under $cand"
fi
[[ "$(<"$cand/default.result")" == "format_report=default=$FJ/out/p1/c-four/default.json" ]] && ok "a format's result file carries its write_report argument" || bad "result file: $(<"$cand/default.result" 2>/dev/null)"
run_oracle "$T4" --out "$FJ/out" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 9 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 ]] && python3 "$FJ/same_report.py" "$FJ/serial/report.json" "$FJ/out/p1/c-four/report.json" && ok "a --format-jobs above the format count runs every format at once and reports the same" || bad "--format-jobs 9: rc=$rc"
red_run "a changed field under --format-jobs 2 is red on the default format alone" change_impact "default: field IMPACT: 1 differing entries" --format-jobs 2
STUB_CAND_RC=4 run_oracle "$T4" --out "$FJ/out-rc" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
cp -R "$FJ/out-rc/p1/c-four/." "$FJ/serial-rc/"
STUB_CAND_RC=4 run_oracle "$T4" --out "$FJ/out-rc" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 4 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$FJ/out-rc/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: candidate exited 4: stub vep: failing on request"*"| tab: candidate exited 4"*"| warnings: the candidate's default run did not complete" ]] && python3 "$FJ/same_report.py" "$FJ/serial-rc/report.json" "$FJ/out-rc/p1/c-four/report.json" &&
    ok "candidate failures in every format job are reported per format in format order, then warnings, as the serial run reports them" || bad "candidate exit 4 under --format-jobs 4: rc=$rc reasons=$reasons"
# The REFUSAL section's green case, both engines refusing both formats, under format jobs.
TFR="$FJ/tf.tsv"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv.gz,gene" default,tab "$IMAGE" offline 116.2 - - - - - - refusals/mixed
} >"$TFR"
STUB_PERL_RC=2 STUB_PERL_STDERR="$PERL_TEXT (Plugin.pm)" STUB_CAND_RC=1 run_oracle "$TFR" --out "$FJ/out-refuse" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
cp -R "$FJ/out-refuse/p1/c-refuse/." "$FJ/serial-refuse/"
STUB_PERL_RC=2 STUB_PERL_STDERR="$PERL_TEXT (Plugin.pm)" STUB_CAND_RC=1 run_oracle "$TFR" --out "$FJ/out-refuse" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 2 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(sed -n 2p "$FJ/out-refuse/status.tsv" | cut -f3)" == "green" ]] && python3 "$FJ/same_report.py" "$FJ/serial-refuse/report.json" "$FJ/out-refuse/p1/c-refuse/report.json" && ok "a refusal case under --format-jobs 2 is green with the serial run's report" || bad "refusal under --format-jobs 2: rc=$rc $(reasons_of "$FJ/out-refuse/p1/c-refuse/report.json" 2>/dev/null)"
STUB_CAND_KILL_JOB=1 run_oracle "$T4" --out "$FJ/out-killed" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 2 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$FJ/out-killed/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: the format's job ended without a result: no stderr | tab: the format's job ended without a result: no stderr"* ]] &&
    ok "a format job that dies before writing its result is red with its stderr's tail (none for a killed job)" || bad "killed format job: rc=$rc reasons=$reasons"
STUB_CAND_KILL_JOB=1 STUB_CAND_JOB_STDERR="the engine lost its cache" run_oracle "$T4" --out "$FJ/out-killed2" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 2 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$FJ/out-killed2/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: the format's job ended without a result: the engine lost its cache  | tab: the format's job ended without a result: the engine lost its cache "* ]] &&
    ok "a format job that wrote to its stderr before dying is red with that text in its reason" || bad "killed format job with stderr: rc=$rc reasons=$reasons"
# An interrupted run stops every running format job's candidate engine at once: the stub
# engine sleeps 120 s and records its pid; the cases name a stored reference, so no container
# runs and the engines are the candidate's. Under --jobs 1 --format-jobs 2 one case runs two
# of its three formats; under --jobs 2 --format-jobs 2 both cases run two formats each.
TFJ="$FJ/t.tsv"
{
    echo "$HEADER"
    row c-fj-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf "$IMAGE" offline 116.2 - - - - - "$REF_CELL" -
    row c-fj-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - "$REF_CELL" -
} >"$TFJ"
for spec in "1 2 2" "2 2 4"; do
    read -r jobs fjobs expected <<<"$spec"
    rm -rf "$FJ/out-int" "$FJ/sleepers"
    mkdir -p "$FJ/out-int"
    echo "stale" >"$FJ/out-int/status.tsv"
    STUB_CAND_SLEEP=120 STUB_CAND_PIDFILE="$FJ/sleepers" PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$TFJ" --data-dir "$DATA" --out "$FJ/out-int" --candidates "$STUBDIR/vep" --candidate-pins p1 --jobs "$jobs" --format-jobs "$fjobs" >"$FJ/int.out" 2>&1 &
    runner=$!
    for _ in $(seq 1 60); do
        [[ -s "$FJ/sleepers" && "$(wc -l <"$FJ/sleepers" | tr -d ' ')" -ge "$expected" ]] && break
        sleep 1
    done
    n=$(wc -l <"$FJ/sleepers" 2>/dev/null | tr -d ' '); n=${n:-0}
    started=$SECONDS
    kill -TERM "$runner" 2>/dev/null
    wait "$runner" 2>/dev/null
    elapsed=$((SECONDS - started))
    sleep 1
    alive=0
    while read -r spid; do kill -0 "$spid" 2>/dev/null && alive=$((alive + 1)); done <"$FJ/sleepers"
    [[ "$n" -eq "$expected" && "$alive" -eq 0 && "$elapsed" -lt 15 ]] && ok "SIGTERM to a --jobs $jobs --format-jobs $fjobs run stops its $n running candidate engine(s) at once (${elapsed}s after the signal, not after the 120 s sleep)" || bad "--jobs $jobs --format-jobs $fjobs: $n engine(s) running at the signal (expected $expected), $alive alive after it, exited after ${elapsed}s"
    [[ ! -e "$FJ/out-int/status.tsv" ]] && ok "an interrupted --jobs $jobs --format-jobs $fjobs run leaves no status.tsv" || bad "a status.tsv survived the interrupted --jobs $jobs --format-jobs $fjobs run"
done
# A format's start line is logged once its job has a slot: with two slots over three formats and
# engines that sleep, the log names two formats while two engines run, and all three at the end.
rm -rf "$FJ/out-slot" "$FJ/slot-sleepers"
STUB_CAND_SLEEP=8 STUB_CAND_PIDFILE="$FJ/slot-sleepers" PATH="$STUBDIR:$PATH" bash "$ORACLE" --cases "$TFJ" --data-dir "$DATA" --out "$FJ/out-slot" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-fj-a --format-jobs 2 >"$FJ/slot.out" 2>&1 &
runner=$!
# The start lines are read while two engines run: two pids recorded and both alive, seen before
# and after the read, so no engine ended and no third started while the log was read.
two_slot_engines_run() {
    local spid
    [[ -s "$FJ/slot-sleepers" && "$(wc -l <"$FJ/slot-sleepers" | tr -d ' ')" -eq 2 ]] || return 1
    while read -r spid; do kill -0 "$spid" 2>/dev/null || return 1; done <"$FJ/slot-sleepers"
}
started_lines=0
two_seen=0
for _ in $(seq 1 60); do
    if two_slot_engines_run; then
        seen=$(grep -c '\[c-fj-a\] p1 \(default\|tab\|vcf\)$' "$FJ/slot.out") || seen=0
        if two_slot_engines_run; then
            started_lines=$seen
            two_seen=1
            break
        fi
    fi
    kill -0 "$runner" 2>/dev/null || break
    sleep 1
done
# The sleeping stub writes no output, so the case ends red; the start lines are what is read.
wait "$runner" || true
[[ "$two_seen" -eq 1 ]] || bad "slot logging: two running engines were not seen before the run ended"
[[ "$started_lines" -eq 2 && "$(grep -c '\[c-fj-a\] p1 \(default\|tab\|vcf\)$' "$FJ/slot.out")" -eq 3 ]] && ok "a format's start line is logged once its job has a slot: two of three formats named while two engines ran, three at the end" || bad "slot logging: $started_lines start lines while two engines ran (expected 2), $(grep -c '\[c-fj-a\] p1 ' "$FJ/slot.out") at the end"

# front_end_switch_reaches_the_comparator_through_the_runner: the runner starts the comparator
# with its own environment, so VEP_FORMATS_REMAINDER set for the runner hands the default and
# tab comparisons of a case over to the front end. One case, its candidate differing in one
# field, is run with the variable unset and with it naming a release build of
# tools/vep-formats-remainder, both runs at the same --out and --work-dir: every format's
# report and the comparator's standard output are the same bytes, and the comparator's
# standard error says what each comparison did. The checks need the build and are skipped,
# with one line that is no `ok` line, where there is none.
echo "FRONT END: the switch reaches the comparator through the runner"
FRONT_END=""
for root in ${CARGO_TARGET_DIR:+"$CARGO_TARGET_DIR"} "$REPO/tools/vep-formats-remainder/target"; do
    if [[ -z "$FRONT_END" && -x "$root/release/vep-formats-remainder" ]]; then
        FRONT_END="$root/release/vep-formats-remainder"
    fi
done
if [[ -z "$FRONT_END" ]]; then
    echo "  skip the front end checks: no release build of tools/vep-formats-remainder"
else
    TFE="$TMP/tfe.tsv"
    {
        echo "$HEADER"
        row c-front fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json "$IMAGE" offline 116.2 - - - - - - -
    } >"$TFE"
    FE="$TMP/fe"
    mkdir -p "$FE"
    for side in whole handed; do
        (
            unset VEP_FORMATS_REMAINDER
            [[ "$side" == whole ]] || export VEP_FORMATS_REMAINDER="$FRONT_END"
            STUB_CAND_MUTATE=change_impact run_oracle "$TFE" --out "$FE/out" --work-dir "$FE/work" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs >"$FE/$side.out" 2>"$FE/$side.err"
        )
        echo "$?" >"$FE/$side.rc"
        mv "$FE/out" "$FE/out-$side"
        mv "$FE/work" "$FE/work-$side"
    done
    [[ "$(<"$FE/whole.rc")" -eq 1 && "$(<"$FE/handed.rc")" -eq 1 ]] && ok "a case differing in one field is red with the variable unset and with it set" || bad "the runner exited $(<"$FE/whole.rc") unset and $(<"$FE/handed.rc") set, expected 1 and 1"
    same=0
    for f in default tab vcf json; do
        cmp -s "$FE/out-whole/p1/c-front/$f.json" "$FE/out-handed/p1/c-front/$f.json" && same=$((same + 1))
    done
    [[ "$same" -eq 4 ]] && ok "every format's report is the same bytes on both sides of the switch" || bad "$same of 4 format reports are the same bytes on both sides of the switch"
    same=0
    for f in default tab vcf json; do
        cmp -s "$FE/work-whole/cases/c-front/candidates/p1/$f.compare.out" "$FE/work-handed/cases/c-front/candidates/p1/$f.compare.out" && same=$((same + 1))
    done
    [[ "$same" -eq 4 ]] && ok "the comparator's standard output is the same bytes on both sides of the switch" || bad "$same of 4 comparator standard outputs are the same bytes on both sides of the switch"
    reasons=$(reasons_of "$FE/out-handed/p1/c-front/report.json")
    [[ "$reasons" == "default: field Extra: 1 differing entries | default: field IMPACT: 1 differing entries" ]] && ok "the handed-over comparison reports the one differing field" || bad "handed-over reasons: $reasons"
    kept=0
    for f in default tab; do
        grep -q "^\[compare_vep_formats\] [0-9]* of 977 reference rows were kept back by $FRONT_END, " "$FE/work-handed/cases/c-front/candidates/p1/$f.compare.err" && kept=$((kept + 1))
    done
    [[ "$kept" -eq 2 ]] && ok "the default and tab comparisons say how many rows the front end kept back" || bad "$kept of 2 text comparisons name the front end on the standard error"
    outside=0
    for f in vcf json; do
        [[ "$(<"$FE/work-handed/cases/c-front/candidates/p1/$f.compare.err")" == "[compare_vep_formats] the whole outputs were compared: the $f format is outside the front end" ]] && outside=$((outside + 1))
    done
    [[ "$outside" -eq 2 ]] && ok "the vcf and json comparisons say they compared the whole outputs" || bad "$outside of 2 comparisons outside the front end say so"
    quiet=0
    for f in default tab vcf json; do
        [[ ! -s "$FE/work-whole/cases/c-front/candidates/p1/$f.compare.err" ]] && quiet=$((quiet + 1))
    done
    [[ "$quiet" -eq 4 ]] && ok "with the variable unset the comparator writes nothing on the standard error" || bad "$quiet of 4 comparisons are silent on the standard error with the variable unset"
fi

echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
