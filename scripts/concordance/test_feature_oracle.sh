#!/usr/bin/env bash
# Tests for feature_oracle.sh: the case-table contract, the dry-run command shapes, every
# fail-closed path, and the verdicts on the real path through a stub `docker` and a stub
# `vep`, both of which answer with the committed Ensembl VEP outputs of
# tests/golden/116/GRCh38-release116 so the real comparator decides. The statistics checks
# compare through a stub named in VEP_STATS_COMPARATOR, never a build of
# tools/vep-stats-compare.
#
# Bash, python3 and a release build of tools/vep-formats-compare, the comparator: the one
# VEP_FORMATS_COMPARATOR names, else the one under CARGO_TARGET_DIR, else the tree's own.
# No Docker, no VEP, no cache.
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
IMAGE="ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
# The suite's cases borrow one corpus as an input under flags of their own, which are no run
# of its manifest, and the runner refuses to compare such a case unnamed where the manifest
# lists several runs over its input. The data tree therefore carries the corpus with its
# manifest cut to its own run and every other entry linked as it is; a relative input
# resolves under --data-dir first. The checks of the named run build corpora of several
# runs of their own.
DATA="$TMP/data"
CORPUS="$DATA/tests/golden/116/GRCh38-release116"
mkdir -p "$CORPUS"
for entry in "$REPO"/tests/golden/116/GRCh38-release116/*; do
    [[ "$(basename "$entry")" == manifest.json ]] || ln -s "$entry" "$CORPUS/"
done
python3 -c 'import json, sys
manifest = json.load(open(sys.argv[1], encoding="utf-8"))
manifest.pop("flag_sets", None)
json.dump(manifest, open(sys.argv[2], "w", encoding="utf-8"), indent=2)' "$REPO/tests/golden/116/GRCh38-release116/manifest.json" "$CORPUS/manifest.json"

if [[ ! -f "$ORACLE" ]]; then
    bad "feature_oracle.sh not found beside this test"
    echo "passed: $PASS   failed: $FAIL"
    exit 1
fi
bash -n "$ORACLE" && ok "feature_oracle.sh parses" || bad "feature_oracle.sh has a syntax error"

# The comparator of every real-path check, exported for the runner.
TOOL="$REPO/tools/vep-formats-compare"
if [[ -z "${VEP_FORMATS_COMPARATOR:-}" ]]; then
    for root in ${CARGO_TARGET_DIR:+"$CARGO_TARGET_DIR"} "$TOOL/target"; do
        if [[ -x "$root/release/vep-formats-compare" ]]; then
            export VEP_FORMATS_COMPARATOR="$root/release/vep-formats-compare"
            break
        fi
    done
fi
if [[ -z "${VEP_FORMATS_COMPARATOR:-}" ]]; then
    bad "no release build of tools/vep-formats-compare: cargo build --release --locked --manifest-path $TOOL/Cargo.toml"
    echo "passed: $PASS   failed: $FAIL"
    exit 1
fi
COMPARATOR="$VEP_FORMATS_COMPARATOR"
# class_table_members <table>...: the classes of the comparator's class tables named, one a
# line. A line of data/classes.tsv is a class, and its cell under a table's column is the
# class's place in that table, or nothing (src/classes.rs of the tool reads it so). A class
# whose CORRECTED_IN cell names an Ensembl VEP release is left out: the comparator applies it
# under no release from that one on, so the runner lists it under no built-in class.
class_table_members() {
    awk -F'\t' -v wanted="$*" '
        NR == 1 {
            n = split(wanted, names, " ")
            for (j = 1; j <= n; j++) {
                for (i = 2; i <= NF; i++) if ($i == names[j]) columns[i] = 1
            }
            for (i = 2; i <= NF; i++) if ($i == "CORRECTED_IN") corrected = i
            next
        }
        corrected && $corrected != "" { next }
        { for (i in columns) if ($i ~ /^[0-9]+$/) { print $1; break } }
    ' "$TOOL/data/classes.tsv"
}

# The committed reference outputs, decompressed once: what both stubs answer with.
REF="$TMP/ref"
mkdir -p "$REF"
for f in default.txt tab.txt vcf.vcf json.jsonl; do
    gzip -dc "$CORPUS/expected/$f.gz" >"$REF/$f"
done

# The fabricated data tree: an Ensembl cache slice, a candidate JSON cache per assembly.
mkdir -p "$DATA/caches/perl/vep-cache/homo_sapiens/116_GRCh38" "$DATA/caches/perl/vep-cache/homo_sapiens/116_GRCh37" \
    "$DATA/caches/vep-rs/grch38/transcripts" "$DATA/caches/vep-rs/grch37/transcripts"

# Stubs. `docker run` and `vep` copy the reference output named by -o from $STUB_REF_DIR,
# write a warnings companion beside a default run when STUB_PERL_WARNINGS / STUB_CAND_WARNINGS
# is set (beside the run of the output file STUB_WARNINGS_OF names, when that is set), copy
# the statistics report STUB_PERL_STATS / STUB_CAND_STATS names beside a default
# run that carries --stats_text (both engines write `<out>_summary.txt` under that flag),
# mutate the candidate's default output per STUB_CAND_MUTATE, exit STUB_PERL_RC /
# STUB_CAND_RC before writing anything (printing STUB_PERL_STDERR / STUB_CAND_STDERR to
# stderr) or STUB_CAND_EXIT_AFTER once the candidate's files are written, write the count
# alone for a --count run (STUB_PERL_COUNT / STUB_CAND_COUNT, 7 by default), sleep in place of a
# run under STUB_PERL_SLEEP / STUB_CAND_SLEEP (recording their pids in STUB_PERL_PIDFILE /
# STUB_CAND_PIDFILE), and append their argv to STUB_DOCKER_LOG / STUB_CAND_LOG. Given a FASTA
# (Ensembl VEP through a directory mounted at /work/fasta, the candidate by --fasta), each
# hands the file to stub-fasta beside it when STUB_FASTA_LOG or its own STUB_PERL_FASTA /
# STUB_CAND_FASTA is set, before it exits or writes anything else.
STUBDIR="$TMP/stubbin"
mkdir -p "$STUBDIR"
# stub-fasta <engine> <FASTA> <words>: what a stub engine does beside the FASTA it is given.
# It appends the engine, the mode of the FASTA's directory and what the directory holds to
# STUB_FASTA_LOG, then writes what the words name: fai and gzi are the two index files,
# other is a .fai of other bytes, lock is <FASTA>.vep.lock and tmp is <FASTA>.fai.tmp; force
# gives the directory its write bit first.
cat >"$STUBDIR/stub-fasta" <<'EOF'
#!/usr/bin/env bash
engine="$1" fasta="$2" words=" $3 "
dir=$(dirname "$fasta")
if [[ -n "${STUB_FASTA_LOG:-}" ]]; then
    mode=$(stat -c %a "$dir" 2>/dev/null || stat -f %Lp "$dir")
    held=$(ls -A "$dir" | tr '\n' ' ')
    printf '%s\t%s\t%s\n' "$engine" "$mode" "${held% }" >>"$STUB_FASTA_LOG"
fi
[[ "$words" != *" force "* ]] || chmod u+w "$dir"
[[ "$words" != *" fai "* ]] || printf 'an index a stub wrote\n' >"$fasta.fai"
[[ "$words" != *" other "* ]] || printf 'an index of other bytes\n' >"$fasta.fai"
[[ "$words" != *" gzi "* ]] || printf 'a block index a stub wrote\n' >"$fasta.gzi"
[[ "$words" != *" lock "* ]] || : >"$fasta.vep.lock"
[[ "$words" != *" tmp "* ]] || : >"$fasta.fai.tmp"
exit 0
EOF
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
stats_text=""
count=""
fasta_dir=""
fasta=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    -v)
        [[ "$2" == *:/work/out ]] && host_out="${2%%:*}"
        # The host directory of a staged FASTA is the one mounted at /work/fasta.
        case "$2" in
        *:/work/fasta | *:/work/fasta:ro) fasta_dir="${2%%:*}" ;;
        esac
        shift 2 ;;
    --fasta) fasta="$2"; shift 2 ;;
    -o) out="$2"; shift 2 ;;
    --stats_text) stats_text=1; shift ;;
    --count) count=1; shift ;;
    *) shift ;;
    esac
done
[[ -n "$host_out" && -n "$out" ]] || { echo "stub docker: no output mapping" >&2; exit 98; }
name=$(basename "$out")
target="$host_out/$name"
if [[ -n "$fasta_dir" && -n "${STUB_FASTA_LOG:-}${STUB_PERL_FASTA:-}" ]]; then
    "${0%/*}/stub-fasta" perl "$fasta_dir/${fasta##*/}" "${STUB_PERL_FASTA:-}"
fi
if [[ -n "$count" ]]; then
    printf '%s\n' "${STUB_PERL_COUNT-7}" >"$target"
    exit 0
fi
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
if [[ "$name" == "${STUB_WARNINGS_OF:-default.txt}" && -n "${STUB_PERL_WARNINGS:-}" ]]; then
    printf '%s\n' "$STUB_PERL_WARNINGS" >"${target}_warnings.txt"
fi
if [[ "$name" == "default.txt" && -n "$stats_text" && -n "${STUB_PERL_STATS:-}" ]]; then
    cp "$STUB_PERL_STATS" "${target}_summary.txt"
fi
exit 0
EOF
cat >"$STUBDIR/vep" <<'EOF'
#!/usr/bin/env bash
[[ -n "${STUB_CAND_LOG:-}" ]] && printf '%s\n' "$*" >>"$STUB_CAND_LOG"
[[ "${1:-}" == "--version" ]] && { echo "vep 0.0.0-stub"; exit 0; }
out=""
stats_text=""
count=""
fasta=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    -o) out="$2"; shift 2 ;;
    --fasta) fasta="$2"; shift 2 ;;
    --stats_text) stats_text=1; shift ;;
    --count) count=1; shift ;;
    *) shift ;;
    esac
done
[[ -n "$out" ]] || { echo "stub vep: no -o" >&2; exit 98; }
if [[ -n "$fasta" && -n "${STUB_FASTA_LOG:-}${STUB_CAND_FASTA:-}" ]]; then
    "${0%/*}/stub-fasta" candidate "$fasta" "${STUB_CAND_FASTA:-}"
fi
if [[ -n "${STUB_CAND_SLEEP:-}" ]]; then
    [[ -z "${STUB_CAND_PIDFILE:-}" ]] || echo "$$" >>"$STUB_CAND_PIDFILE"
    exec sleep "$STUB_CAND_SLEEP"
fi
# STUB_CAND_KILL_JOB kills the shell that runs the engine, a format job under --format-jobs,
# after writing STUB_CAND_JOB_STDERR, when set, to that job's own stderr file beside the output.
# The engine's parent is the runner's step timer, so the job's shell is the timer's parent.
if [[ -n "${STUB_CAND_KILL_JOB:-}" ]]; then
    fmt=$(basename "$out")
    [[ -z "${STUB_CAND_JOB_STDERR:-}" ]] || printf '%s\n' "$STUB_CAND_JOB_STDERR" >>"$(dirname "$out")/${fmt%%.*}.job.err"
    kill -KILL "$(ps -o ppid= -p "$PPID" | tr -d ' ')"
    exit 0
fi
[[ "${STUB_CAND_RC:-0}" -eq 0 ]] || { echo "${STUB_CAND_STDERR:-stub vep: failing on request}" >&2; exit "${STUB_CAND_RC}"; }
name=$(basename "$out")
if [[ -n "$count" ]]; then
    printf '%s\n' "${STUB_CAND_COUNT-7}" >"$out"
    exit 0
fi
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
    [[ -z "$stats_text" || -z "${STUB_CAND_STATS:-}" ]] || cp "$STUB_CAND_STATS" "${out}_summary.txt"
fi
if [[ "$name" == "${STUB_WARNINGS_OF:-default.txt}" && -n "${STUB_CAND_WARNINGS:-}" ]]; then
    printf '%s\n' "$STUB_CAND_WARNINGS" >"${out}_warnings.txt"
fi
exit "${STUB_CAND_EXIT_AFTER:-0}"
EOF
# The statistics comparator of every stats check. It takes the tool's three options and no
# other, appends its arguments to STUB_STATS_LOG, exits STUB_STATS_RC without a report when
# that is set, and otherwise writes a report with the tool's `mode`, `format`, `verdict` and
# `reasons`: green and exit 0 where the two files hold the same text (a .gz read as its
# text), red and exit 1 where they differ.
cat >"$STUBDIR/stats-compare" <<'EOF'
#!/usr/bin/env bash
[[ -n "${STUB_STATS_LOG:-}" ]] && printf '%s\n' "$*" >>"$STUB_STATS_LOG"
expected="" actual="" report=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    --expected) expected="$2"; shift 2 ;;
    --actual) actual="$2"; shift 2 ;;
    --report) report="$2"; shift 2 ;;
    *) echo "stub statistics comparator: unexpected argument $1" >&2; exit 2 ;;
    esac
done
if [[ -n "${STUB_STATS_RC:-}" ]]; then
    echo "ERROR: [compare_vep_stats] failing on request" >&2
    exit "$STUB_STATS_RC"
fi
text_of() {
    case "$1" in
    *.gz) gzip -dc "$1" ;;
    *) cat "$1" ;;
    esac
}
if [[ "$(text_of "$expected")" == "$(text_of "$actual")" ]]; then
    printf '{"mode": "stats", "format": "stats", "verdict": "green", "reasons": []}\n' >"$report"
    exit 0
fi
printf '{"mode": "stats", "format": "stats", "verdict": "red", "reasons": ["[General statistics]: 1 difference"]}\n' >"$report"
exit 1
EOF
chmod +x "$STUBDIR/docker" "$STUBDIR/vep" "$STUBDIR/stats-compare" "$STUBDIR/stub-fasta"
export VEP_STATS_COMPARATOR="$STUBDIR/stats-compare"
cp "$STUBDIR/vep" "$STUBDIR/vep2"
# The stored-reference checks replace the docker stub with a failing one; this copy restores it.
cp "$STUBDIR/docker" "$TMP/docker-answering"
export STUB_REF_DIR="$REF"

# row <cells...>: one tab-separated line.
row() {
    local IFS=$'\t'
    printf '%s\n' "$*"
}
HEADER=$(row case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs species unordered reference refusal tool)
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
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json,warnings "$IMAGE" offline 116.2 - - - - - - - -
    row c-db fixture "$CORPUS_INPUT" 116_GRCh37 - "--format hgvs" default "$IMAGE" database 116.2 breakend_mate_local_read IMPACT=3,Feature=2 "IMPACT=LOW=>HIGH;Consequence=a=>b" - - - - -
    row c-pq fixture "$CORPUS_INPUT" 116_GRCh38 - - tab,parquet "$IMAGE" offline 116.2 - - - - - - - -
    row c-near fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest symbol" default,tab "$IMAGE" offline 116.2 - - - - - - - -
    row c-skip fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 reference_skipped_record=4 - - - - - - -
    row GRCh38-release116 population "$CORPUS_INPUT" 116_GRCh38 - - default "$IMAGE" offline 116.2 - - - - - - - -
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
# vep reads its options with Getopt::Long's defaults (vep:27, Config.pm:62 at release/116.3):
# a single dash, any unambiguous abbreviation and any letter case name the option, so a case
# spelled -nearest, --near or --NEAREST is a --nearest case, its flag word passed as written.
T_SP="$TMP/t-spellings.tsv"
{
    echo "$HEADER"
    row c-near-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-nearest symbol" default,tab "$IMAGE" offline 116.2 - - - - - - - -
    row c-near-abbr fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--near symbol" default,tab "$IMAGE" offline 116.2 - - - - - - - -
    row c-near-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--NEAREST symbol" default,tab "$IMAGE" offline 116.2 - - - - - - - -
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

echo
echo "DRY-RUN STATS: the default run under --stats_text and the statistics comparator"
mkdir -p "$DATA/stored/stats"
cp "$CORPUS/expected/default.txt.gz" "$DATA/stored/stats/default.txt.gz"
printf '[VEP run statistics]\nSpecies\thomo_sapiens\n' >"$DATA/stored/stats/default.summary.log"
T1S="$TMP/t1s.tsv"
{
    echo "$HEADER"
    row c-stats fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,stats "$IMAGE" offline 116.2 order_only - - - - - - -
    row c-stats-flag fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--stats_text --sift b" default,stats "$IMAGE" offline 116.2 - - - - - - - -
    row c-stats-stored fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,stats "$IMAGE" offline 116.2 - - - - - stored/stats - -
} >"$T1S"
dry "$T1S" "$TMP/d1s.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "a dry run over stats cases exits 0" || {
    bad "stats dry run exited $rc"
    head -5 "$TMP/d1s.out.err" | sed 's/^/         /'
}
spc=$(cmd_after 'DRY-RUN c-stats perl default' "$TMP/d1s.out")
[[ "$spc" == *"--format vcf --fork 1 --stats_text --quiet --force_overwrite "* && "$spc" != *"--no_stats"* ]] && ok "the Ensembl VEP default run of a stats case carries --stats_text where --no_stats stood" || bad "stats default run: $spc"
spt=$(cmd_after 'DRY-RUN c-stats perl tab' "$TMP/d1s.out")
[[ "$spt" == *"--format vcf --fork 1 --no_stats --quiet --force_overwrite "* && "$spt" != *"--stats_text"* ]] && ok "the other formats' Ensembl VEP runs of a stats case keep --no_stats" || bad "stats tab run: $spt"
scc=$(cmd_after 'DRY-RUN c-stats stub default' "$TMP/d1s.out")
sct=$(cmd_after 'DRY-RUN c-stats stub tab' "$TMP/d1s.out")
[[ "$scc" == *"--format vcf --stats_text --quiet --force_overwrite "* && "$scc" != *"--no_stats"* && "$sct" == *"--format vcf --no_stats --quiet --force_overwrite "* && "$sct" != *"--stats_text"* ]] &&
    ok "the candidate's default run carries --stats_text in place of --no_stats, and its other runs --no_stats" || bad "stats candidate runs: default=$scc tab=$sct"
ssc=$(cmd_after 'DRY-RUN c-stats stub stats' "$TMP/d1s.out")
[[ "$ssc" == "$STUBDIR/stats-compare --expected $TMP/never/work/cases/c-stats/perl/default.txt_summary.txt --actual $TMP/never/work/cases/c-stats/candidates/stub/default.txt_summary.txt --report $TMP/never/stub/c-stats/stats.json" ]] &&
    ok "stats compares the default runs' _summary.txt reports through the statistics comparator, which is given no masks" || bad "stats command: $ssc"
n=$(count_lines '^DRY-RUN c-stats perl ' "$TMP/d1s.out")
[[ "$n" -eq 2 ]] && ok "stats adds no Ensembl VEP run of its own (2 for default,tab)" || bad "$n Ensembl VEP runs planned for c-stats, expected 2"
sfp=$(cmd_after 'DRY-RUN c-stats-flag perl default' "$TMP/d1s.out")
sfc=$(cmd_after 'DRY-RUN c-stats-flag stub default' "$TMP/d1s.out")
[[ "$sfp" == *"--fork 1 --quiet --force_overwrite --stats_text --sift b" && "$sfc" == *"--format vcf --quiet --force_overwrite --stats_text --sift b" ]] && ok "a case whose flags carry --stats_text is given it once, by its flags" || bad "stats flag once: perl=$sfp cand=$sfc"
sts=$(cmd_after 'DRY-RUN c-stats-stored stub stats' "$TMP/d1s.out")
[[ "$sts" == "$STUBDIR/stats-compare --expected $DATA/stored/stats/default.summary.log --actual "* ]] && ! grep -q '^DRY-RUN c-stats-stored perl default' "$TMP/d1s.out" &&
    ok "a stats case with a stored reference compares against its default.summary.log and plans no Ensembl VEP run" || bad "stored stats command: $sts"
# A table without a stats case plans what it planned before, and needs no statistics comparator.
VEP_STATS_COMPARATOR="$TMP/no-such-statistics-comparator" dry "$T1" "$TMP/d1-nostats.out"
rc=$?
[[ "$rc" -eq 0 ]] && cmp -s <(grep -v '^\[' "$TMP/d1.out") <(grep -v '^\[' "$TMP/d1-nostats.out") && ok "a table without a stats case plans the same commands under a VEP_STATS_COMPARATOR that names nothing" || bad "no stats case, no statistics comparator: rc=$rc $(head -2 "$TMP/d1-nostats.out.err" | tr '\n' ' ')"
VEP_STATS_COMPARATOR="$TMP/no-such-statistics-comparator" dry "$T1S" "$TMP/d1s-none.out"
rc=$?
[[ "$rc" -eq 1 && "$(<"$TMP/d1s-none.out.err")" == "ERROR: [feature_oracle] the statistics comparator $TMP/no-such-statistics-comparator is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO/tools/vep-stats-compare/Cargo.toml' or name a built one in VEP_STATS_COMPARATOR" ]] && ! grep -q '^DRY-RUN' "$TMP/d1s-none.out" &&
    ok "a stats case under a VEP_STATS_COMPARATOR that is not an executable stops the run on one line before any command is printed" || bad "stats case without a statistics comparator: rc=$rc $(head -2 "$TMP/d1s-none.out.err" | tr '\n' ' ')"

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
    row c-fa fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-yeast fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - - default "$IMAGE" offline 116.2 - - - saccharomyces_cerevisiae - - - -
    row c-human fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
    row c-yeast-db fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - --genomes default "$IMAGE" database 116.2 - - - saccharomyces_cerevisiae - - - -
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
    row c-yeast-db fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - --genomes default "$IMAGE" database 116.2 - - - saccharomyces_cerevisiae - - - -
    row c-yeast fixture "$CORPUS_INPUT" "116_R64-1-1=caches/vep-rs/yeast" - - default "$IMAGE" offline 116.2 - - - saccharomyces_cerevisiae - - - -
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
    row c-both fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$PH_FLAGS" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-corpus fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--custom {corpus}/custom/regions.bed.gz,RE,bed,overlap" default "$IMAGE" offline 116.2 - - - - - - - -
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
# A printed command writes a brace behind a backslash (printf %q): a word left as it is reads \{corpus\}.
grep '^docker run \|^/' "$TMP/dp.out" | grep -E '(corpus|sources)\\?\}' >/dev/null && bad "a placeholder reached an engine command unrewritten" || ok "no engine command carries a placeholder"
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
# manifest_words_take_the_corpus_directories: beside a manifest.json a flags cell is written as
# one of the manifest's runs writes it, {plugins} standing for <corpus>/plugin_data and
# {generated} for <corpus>/generated; each engine reads the directory on the path {corpus}
# takes for it, and the cell is that run of the manifest.
PG="$DATA/manifestwords"
mkdir -p "$PG/plugin_data" "$PG/generated/fordownload"
cp "$CORPUS/variants.vcf" "$PG/variants.vcf"
cat >"$PG/manifest.json" <<'JSON'
{"flags": ["--dir_plugins", "{plugins}", "--plugin", "StubScores,snv={corpus}/plugin_data/snv.vcf.gz"],
 "flag_sets": [
  {"name": "made", "flags": ["--dir_plugins", "{corpus}", "--plugin", "StubMatrices,{generated}/fordownload,SWA"], "input": "variants.vcf"}]}
JSON
TPG="$TMP/tpg.tsv"
{
    echo "$HEADER"
    row c-plugins fixture manifestwords/variants.vcf "$CORPUS_CACHE" - "--dir_plugins {plugins} --plugin StubScores,snv={corpus}/plugin_data/snv.vcf.gz" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-generated fixture manifestwords/variants.vcf "$CORPUS_CACHE" - "--dir_plugins {corpus} --plugin StubMatrices,{generated}/fordownload,SWA" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TPG"
dry "$TPG" "$TMP/dpg.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "manifest_words_take_the_corpus_directories: a table whose cells are written as a manifest's runs write {plugins} and {generated} plans" || {
    bad "manifest words dry run exited $rc"
    head -5 "$TMP/dpg.out.err" | sed 's/^/         /'
}
pc=$(cmd_after 'DRY-RUN c-plugins perl default' "$TMP/dpg.out")
cc=$(cmd_after 'DRY-RUN c-plugins stub default' "$TMP/dpg.out")
perl_words=$(printf '%q ' --dir_plugins /work/corpus/plugin_data --plugin "StubScores,snv=/work/corpus/plugin_data/snv.vcf.gz")
[[ "$pc" == *"-v $PG:/work/corpus:ro $IMAGE vep "* && "$pc" == *"--force_overwrite ${perl_words% }" ]] && ok "Ensembl VEP reads {plugins} as the plugin_data directory of the corpus mount, the corpus mounted read-only" || bad "{plugins}, Ensembl VEP: $pc"
cand_words=$(printf '%q ' --dir_plugins "$PG/plugin_data" --plugin "StubScores,snv=$PG/plugin_data/snv.vcf.gz")
[[ "$cc" == *"--force_overwrite ${cand_words% }" ]] && ok "the candidate reads {plugins} as the corpus's plugin_data directory" || bad "{plugins}, candidate: $cc"
pc=$(cmd_after 'DRY-RUN c-generated perl default' "$TMP/dpg.out")
cc=$(cmd_after 'DRY-RUN c-generated stub default' "$TMP/dpg.out")
perl_words=$(printf '%q ' --dir_plugins /work/corpus --plugin "StubMatrices,/work/corpus/generated/fordownload,SWA")
[[ "$pc" == *"-v $PG:/work/corpus:ro $IMAGE vep "* && "$pc" == *"--force_overwrite ${perl_words% }" ]] && ok "Ensembl VEP reads {generated} as the generated directory of the corpus mount" || bad "{generated}, Ensembl VEP: $pc"
cand_words=$(printf '%q ' --dir_plugins "$PG" --plugin "StubMatrices,$PG/generated/fordownload,SWA")
[[ "$cc" == *"--force_overwrite ${cand_words% }" ]] && ok "the candidate reads {generated} as the corpus's generated directory" || bad "{generated}, candidate: $cc"
grep -E '^(docker run |/).*(plugins|generated)\\?\}' "$TMP/dpg.out" >/dev/null && bad "a manifest word reached an engine command unrewritten" || ok "no engine command carries {plugins} or {generated}"
mcp=$(cmd_second 'DRY-RUN c-plugins stub default' "$TMP/dpg.out")
mcg=$(cmd_second 'DRY-RUN c-generated stub default' "$TMP/dpg.out")
[[ "$mcp" == *"--manifest $PG/manifest.json --run default"* && "$mcg" == *"--manifest $PG/manifest.json --run made"* ]] && ok "each cell is the run of the manifest that writes those words: the manifest's own, and the flag set" || bad "manifest word runs: $mcp | $mcg"
# Without a manifest.json beside the input the two words stand for nothing.
mkdir -p "$DATA/barewords"
cp "$CORPUS/variants.vcf" "$DATA/barewords/variants.vcf"
TPB="$TMP/tpb.tsv"
{
    echo "$HEADER"
    row c-plugins-bare fixture barewords/variants.vcf "$CORPUS_CACHE" - "--dir_plugins {plugins}" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-generated-bare fixture barewords/variants.vcf "$CORPUS_CACHE" - "--plugin StubMatrices,{generated}/fordownload" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TPB"
dry "$TPB" "$TMP/dpb.out"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'c-plugins-bare: flags carry {plugins} or {generated}, .* and no manifest.json sits beside the input' "$TMP/dpb.out.err" && grep -q 'c-generated-bare: flags carry {plugins} or {generated}' "$TMP/dpb.out.err" && ! grep -q '^DRY-RUN' "$TMP/dpb.out" &&
    ok "a {plugins} or {generated} word with no manifest.json beside the input fails preflight by case, and nothing is planned" || bad "manifest words without a manifest exited $rc"

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
    row c-dir fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {sources}/gwas/vepdir --plugin StubPlugin" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-eq fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir={corpus} --plugin StubPlugin" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-word-free fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-HGVS +Symbol" default "$IMAGE" offline 116.2 - - - - - - - -
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
# (Config.pm:704-705 at release/116.3), so each such cell is refused by case and by name;
# the value is read after `=` or from the next word, in any spelling of the flag.
# runner_word_is_refused_in_every_spelling: vep reads a single dash or a plus sign and any
# letter case as the option (vep:27, Config.pm:62), so a word the runner sets is refused
# in those spellings too, and a word it does not set passes in them.
{
    echo "$HEADER"
    row c-dir-abs fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir $DATA/sources/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-rel fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir=gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-up fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {corpus}/../GRCh37-hgvs" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-bare fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--plugin StubPlugin --dir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-dir /work/sources/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-ini fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {sources}/ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-ini-eq fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--DIR={sources}/ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-cache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir_cache {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-config fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--config {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-two fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {corpus}/{sources}/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-word fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir {corpus}/{plugins}" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dir-word-eq fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir={corpus}/made/{generated}" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-word-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-config {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-word-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--CONFIG={sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-word-dircache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-dir_cache {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
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
# Beside a manifest.json a manifest's word is read as a path under {corpus}, so one after the
# directory's placeholder is a second placeholder; the input's directory carries a manifest.json.
grep -qF "c-dir-word: flags give --dir '{corpus}/{plugins}', $not_a_mount" "$TMP/dd2.out.err" && grep -qF "c-dir-word-eq: flags give --dir '{corpus}/made/{generated}', $not_a_mount" "$TMP/dd2.out.err" &&
    ok "a --dir directory carrying {plugins} or {generated} after its placeholder is refused as a second placeholder" || bad "--dir with a manifest's word: $(grep -F 'c-dir-word' "$TMP/dd2.out.err")"
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
    row c-abbr-value fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--conf {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-value-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-conf {sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-value-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--CONF={sources}/ini/vep.ini" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-switch fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--offl" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-switch-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-offl" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-switch-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--OFFL" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-begun fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--js" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-begun-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-js" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-begun-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--JS" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-dircache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--dir_c {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-candidate fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--thr 2" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-name-user fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-u anonymous" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-name-fullcache fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--full_cache_dir {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-abbr-free fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--symb --dir_p {sources}/gwas/vepdir" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-shared fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--vc" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-abbr-shared-two fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--cach" default "$IMAGE" offline 116.2 - - - - - - - -
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
echo "ON DEMAND: a case of oracle_cases_on_demand.tsv beside the table runs only when --only names it"
# The list is read beside the case table, so each reading has a table in a directory of its own.
mkdir -p "$TMP/od" "$TMP/od-unknown" "$TMP/od-fixture" "$TMP/od-header"
{
    cat "$T1"
    row pop-asked population "$CORPUS_INPUT" 116_GRCh38 - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TMP/od/cases.tsv"
printf '# a comment line\n\ncase\treason\npop-asked\ta run of its own\n' >"$TMP/od/oracle_cases_on_demand.tsv"
dry "$TMP/od/cases.tsv" "$TMP/od.out" --tier population
rc=$?
n=$(count_lines '^DRY-RUN GRCh38-release116 ' "$TMP/od.out")
m=$(count_lines '^DRY-RUN pop-asked ' "$TMP/od.out")
[[ "$rc" -eq 0 && "$n" -gt 0 && "$m" -eq 0 ]] && ok "a run without --only leaves an on-demand case of its tier out and plans the tier's other cases" || bad "on demand without --only: rc=$rc GRCh38-release116=$n pop-asked=$m"
n=$(count_lines 'left out, on demand (--only names a case to run it): pop-asked$' "$TMP/od.out")
[[ "$n" -eq 1 ]] && ok "the run logs one line that names the on-demand cases it left out" || bad "the left-out line was logged $n times"
dry "$TMP/od/cases.tsv" "$TMP/od-named.out" --tier population --only pop-asked
rc=$?
n=$(count_lines '^DRY-RUN pop-asked ' "$TMP/od-named.out")
m=$(count_lines 'left out, on demand' "$TMP/od-named.out")
[[ "$rc" -eq 0 && "$n" -gt 0 && "$m" -eq 0 ]] && ok "--only runs the on-demand case it names, and no case is logged as left out" || bad "on demand under --only: rc=$rc pop-asked=$n left-out lines=$m"
dry "$TMP/od/cases.tsv" "$TMP/od-other.out"
rc=$?
n=$(count_lines '^DRY-RUN c-four ' "$TMP/od-other.out")
m=$(count_lines 'left out, on demand' "$TMP/od-other.out")
[[ "$rc" -eq 0 && "$n" -gt 0 && "$m" -eq 0 ]] && ok "a run of another tier plans its cases and logs no case as left out" || bad "on demand, the fixture tier: rc=$rc c-four=$n left-out lines=$m"
cp "$T1" "$TMP/od-unknown/cases.tsv"
printf 'case\treason\npop-nowhere\ta case the table lacks\n' >"$TMP/od-unknown/oracle_cases_on_demand.tsv"
dry "$TMP/od-unknown/cases.tsv" "$TMP/od-unknown.out" --tier population
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "oracle_cases_on_demand.tsv names 'pop-nowhere', which is not a case of " "$TMP/od-unknown.out.err" && ! grep -q '^DRY-RUN' "$TMP/od-unknown.out" &&
    ok "an on-demand row for a case the table lacks is an error before anything is planned" || bad "an on-demand row for an unknown case exited $rc: $(head -1 "$TMP/od-unknown.out.err")"
cp "$T1" "$TMP/od-fixture/cases.tsv"
printf 'case\treason\nc-four\ta fixture case\n' >"$TMP/od-fixture/oracle_cases_on_demand.tsv"
dry "$TMP/od-fixture/cases.tsv" "$TMP/od-fixture.out" --tier population
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "oracle_cases_on_demand.tsv names 'c-four', a fixture case; the list holds population cases alone" "$TMP/od-fixture.out.err" && ! grep -q '^DRY-RUN' "$TMP/od-fixture.out" &&
    ok "an on-demand row for a fixture case is an error, whatever tier runs" || bad "an on-demand row for a fixture case exited $rc: $(head -1 "$TMP/od-fixture.out.err")"
cp "$T1" "$TMP/od-header/cases.tsv"
printf 'case\tnote\nGRCh38-release116\tanother header\n' >"$TMP/od-header/oracle_cases_on_demand.tsv"
dry "$TMP/od-header/cases.tsv" "$TMP/od-header.out" --tier population
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "oracle_cases_on_demand.tsv:1: the header is not case, reason" "$TMP/od-header.out.err" && ok "an on-demand list under another header is an error" || bad "an on-demand list under another header exited $rc: $(head -1 "$TMP/od-header.out.err")"
printf 'case\treason\nGRCh38-release116\n' >"$TMP/od-header/oracle_cases_on_demand.tsv"
dry "$TMP/od-header/cases.tsv" "$TMP/od-row.out" --tier population
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "oracle_cases_on_demand.tsv:2: a row is case, reason" "$TMP/od-row.out.err" && ok "an on-demand row without its reason is an error" || bad "an on-demand row without its reason exited $rc: $(head -1 "$TMP/od-row.out.err")"
printf 'case\treason\nGRCh38-release116 \ta blank after the name\n' >"$TMP/od-header/oracle_cases_on_demand.tsv"
dry "$TMP/od-header/cases.tsv" "$TMP/od-id.out" --tier population
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "oracle_cases_on_demand.tsv:2: 'GRCh38-release116 ' is not a case's name" "$TMP/od-id.out.err" && ! grep -q '^DRY-RUN' "$TMP/od-id.out" &&
    ok "an on-demand row whose case cell is no case's name is an error, not a case left out" || bad "an on-demand row whose case cell is no case's name exited $rc: $(head -1 "$TMP/od-id.out.err")"
printf 'case\treason\nGRCh38-release116\tthe tier holds no other case\n' >"$TMP/od-header/oracle_cases_on_demand.tsv"
dry "$TMP/od-header/cases.tsv" "$TMP/od-empty.out" --tier population
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'no case of tier population selected from .*; left out, on demand (--only names a case to run it): GRCh38-release116$' "$TMP/od-empty.out.err" &&
    ok "a tier whose every case is on demand is an error without --only that names the cases left out, not an empty success" || bad "a tier of on-demand cases alone exited $rc: $(head -1 "$TMP/od-empty.out.err")"

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
    row c-miss fixture "$DATA/absent.vcf" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a missing case input exits non-zero and names it" "c-miss: input not found: $DATA/absent.vcf" "$T3"
grep -q 'preflight failed; nothing ran' "$TMP/e.out.err" && ok "the abort is reported as a preflight failure (the other case was fine)" || bad "no preflight-failure line"
{
    echo "$HEADER"
    row c-fa fixture "$CORPUS_INPUT" "$CORPUS_CACHE" "$DATA/absent.fa" - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a missing FASTA is an error" "FASTA not found" "$T3"
{
    echo "$HEADER"
    row c-cache fixture "$CORPUS_INPUT" 116_GRCh38 - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a missing Ensembl cache slice is an error in offline mode" "Ensembl VEP cache slice not found: $DATA/x/homo_sapiens/116_GRCh38" "$T3" --perl-cache-dir "$DATA/x"
expect_error "a candidate JSON cache without transcripts/ is an error" "candidate JSON cache has no transcripts/" "$T3" --json-cache-dir "$DATA/x"
{
    echo "$HEADER"
    row c-fmt fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default,xml "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "an unknown format is an error" "formats names 'xml'; the formats are default,tab,vcf,json,parquet, warnings and stats" "$T3"
{
    echo "$HEADER"
    row c-w fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab,vcf,warnings "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "warnings without default beside two output formats is an error" "c-w: formats names warnings without default beside 'tab vcf'" "$T3"
{
    echo "$HEADER"
    row c-s fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab,stats "$IMAGE" offline 116.2 - - - - - - - -
    row c-sf fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--stats_file other.txt" default,stats "$IMAGE" offline 116.2 - - - - - - - -
    row c-sf-name fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-sf other.txt" default,stats "$IMAGE" offline 116.2 - - - - - - - -
    row c-sf-free fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--stats_file other.txt" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-s-only fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - stats "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "stats without default is an error" "c-s: formats names stats without default; the statistics compared are the default run's" "$T3" --only c-s
expect_error "a --stats_file in the flags of a stats case is an error" "c-sf: flags names --stats_file with stats in formats; the statistics compared are the default run's <out>_summary.txt" "$T3" --only c-sf
expect_error "the option's other name, -sf, is refused the same way" "c-sf-name: flags names --stats_file with stats in formats" "$T3" --only c-sf-name
dry "$T3" "$TMP/d-sf-free.out" --only c-sf-free
rc=$?
[[ "$rc" -eq 0 ]] && ok "control: --stats_file in the flags of a case without stats plans" || bad "--stats_file without stats exited $rc: $(head -2 "$TMP/d-sf-free.out.err" | tr '\n' ' ')"
expect_error "stats alone names no output format, and the message lists it beside warnings" "c-s-only: formats 'stats' names no output format (default,tab,vcf,json,parquet, and warnings or stats)" "$T3" --only c-s-only
{
    echo "$HEADER"
    row c-flag fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--offline --hgvs" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a flags cell naming a runner-owned flag is an error" "flags names --offline" "$T3"
{
    echo "$HEADER"
    row c-mode fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" online 116.2 - - - - - - - -
} >"$T3"
expect_error "an unknown perl_mode is an error" "perl_mode 'online'" "$T3"
{
    echo "$HEADER"
    row c-tier nightly "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "an unknown tier is an error" "tier 'nightly'" "$T3"
{
    echo "$HEADER"
    row c-dup fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
    row c-dup fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a duplicate case is an error" "case 'c-dup' appears twice" "$T3"
{
    echo "$HEADER"
    row "../escape" fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a case that is not a directory-safe name is an error (it names the work and report directories)" "case '../escape' is not a directory-safe name" "$T3"
{
    echo "$HEADER"
    row c-empty fixture "$CORPUS_INPUT" "$CORPUS_CACHE" "" - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "an empty cell is an error (write -)" "an empty cell" "$T3"
{
    echo "$HEADER"
    row c-short fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a short row is an error" "17 columns, expected 18" "$T3"
{
    row case input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs species unordered reference refusal
    row c-x "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a header without the tier column is an error" "header lacks the column 'tier'" "$T3"
{
    row case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs unordered reference refusal
    row c-x fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$T3"
expect_error "a header without the species column is an error" "header lacks the column 'species'" "$T3"
{
    echo "$HEADER"
    row c-suffixed fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - homo_sapiens_refseq - - - -
} >"$TS"
expect_error "a species carrying the slice suffix is refused, since Ensembl VEP refuses it as --species" "species 'homo_sapiens_refseq' carries the slice suffix; write the species and the --refseq flag" "$TS"
{
    echo "$HEADER"
    row c-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - Homo_sapiens - - - -
} >"$TS"
expect_error "a species that is not a lower-case production name is refused" "species 'Homo_sapiens' is not an Ensembl species name" "$TS"
{
    echo "$HEADER"
    row c-masks fixture "$DATA/caches/vep-rs/grch38/transcripts" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 some_class - - - - - - -
} >"$T3"
expect_error "masks without a manifest.json beside the input is an error" "masks names 'some_class', which applies from a manifest.json beside the input, and none is there (the built-in classes are regulatory_feature_rows order_only individual_allele_order conservation_short_name_spelling pc_description_text motif_score_reverse_strand_mnv colocated_list_cut_at_cache_cap start_cooccurrence_swap splice_family_swap gxf_reverse_strand_utr_sequence organelle_peptide_paths_translation_table conservation_max_text_order gwas_haplotype_cross_product gwas_interaction_first_snp_allele_whitespace alphamissense_two_digit_transcript_version cadd_structural_allele_match splicevault_sibling_allele_rows splicevault_json_later_rows utrannotator_annotation_pair_order utrannotator_gained_start_window utrannotator_consequence_precedence utrannotator_structural_allele_text utrannotator_lost_start_after_insertion utrannotator_gained_start_remade_in_place utrannotator_frameshift_at_stop_boundary utrannotator_lost_stop_written_before_stop_codon utrannotator_without_reference nearest_exon_jb_item_order loeuf_zero_minimum maxentscan_overlapping_windows alphamissense_sibling_allele_rows revel_sibling_allele_rows revel_other_reference_amino_acid_rows gnomadmt_info_value_second_equals loftee_anc_allele_multiallelic_record loftee_first_alt_donor_filter loftee_gerp_zoom_level_summary loftee_gerp_walk_without_conservation_database loftee_nagnag_site_first_transcript_strand nearest_exon_jb_stretch_priming nearest_exon_jb_shifted_span nearest_exon_jb_id_less_exon_fold nearest_start_site_one_base_above downstream_reverse_strand_rotation shifted_allele_rotation_by_strand spliceai_chopped_multiallelic spliceai_multi_entry_cell plugin_contig_mitochondrial_spelling tabix_allele_case_insensitive_match refseq_edited_model_hgvsc_offset refseq_edited_model_shifted_insertion_bases refseq_edited_model_seam_insertion_hgvsc refseq_edited_model_unminimised_allele_hgvsc refseq_edited_model_hgvs_offset refseq_edited_model_deleted_span_term refseq_edited_model_tail_term refseq_edited_model_equal_span_insertion_term mature_mirna_unshifted_span stop_codon_unshifted_span structural_deletion_stop_lost_without_stop_codon stop_codon_completed_by_polyadenylation stop_codon_overlap_without_stop_codon stop_codon_behind_last_residue_without_stop_codon start_codon_overlap_without_start_codon stop_retained_without_stop_codon_test stop_codon_overlap_genomic_span stop_retained_three_base_stop_codon_test stop_codon_removed_without_stop_codon_test stop_codon_overlap_insertion_range stop_codon_settled_before_unknown_bases frameshift_inside_a_held_stop_codon start_codon_removed_without_five_prime_utr start_codon_kept_without_five_prime_utr cnv_tr_symbolic_loss_named_a_gain splice_term_past_the_record_span refseq_edited_model_deleted_span_codons alternate_codon_past_the_coding_end refseq_edited_model_intron_span_coding_reference refseq_edited_model_inserted_span_codons hgvs_intron_offset_abutting_insertion hgvs_shifted_insertion_clamped_end shift_intron_term_unshifted_span hgvs_shifted_indel_written_span_gate hgvsp_shifted_indel_written_span_gate hgvsp_shifted_indel_written_span_start_test hgvsp_start_codon_test_past_the_shifted_allele hgvs_shifted_deletion_clamped_span hgvs_clipped_reference_against_whole_alternate hgvs_added_bases_left_at_the_span_end hgvs_insertion_left_where_the_shift_stops hgvs_added_bases_written_outside_the_transcript hgvs_shift_sequence_edge_window hgvs_shift_shared_across_alleles hgvsp_stop_codon_protein_comparison hgvsp_stop_named_on_incomplete_coding_sequence hgvsp_frameshift_count_past_a_stop_read_as_a_residue hgvsp_frameshift_with_no_changed_residue hgvsp_frameshift_count_at_the_first_residue allele_number_repeated_or_reference_equal_alt hgvsp_intron_end_insertion hgvsp_intron_start_insertion merged_same_stable_id_collapse use_transcript_ref_cache_cdna refseq_failed_edit_model_hgvsc_offset hgnc_id_copy downstream_utr_frame_past_coding_end downstream_start_past_own_codon plugin_hgvs_shifted_cds_coordinates plugin_value_follows_splice_term_past_the_record_span plugin_value_follows_stop_codon_term_pair total_length_json_range_end plugin_value_on_reference_skipped_record blosum62_hgvs_shift_shared_across_alleles reference_skipped_record breakend_mate_local_read bnd_own_chromosome_rows_lost dosage_sensitivity_regulatory_row_warnings fork_warning_pid_prefix json_mismatch_warning_names_no_record)" "$T3"
mkdir -p "$DATA/nomanifest"
cp "$CORPUS/variants.vcf" "$DATA/nomanifest/variants.vcf"
{
    echo "$HEADER"
    row c-builtin fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,regulatory_feature_rows,individual_allele_order,conservation_short_name_spelling,pc_description_text - - - - - - -
    row c-mixed fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,some_class - - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-builtin.out" --only c-builtin
rc=$?
bc=$(cmd_second 'DRY-RUN c-builtin stub default' "$TMP/d-builtin.out")
[[ "$rc" -eq 0 && "$bc" == *'--masks order_only\,regulatory_feature_rows\,individual_allele_order\,conservation_short_name_spelling\,pc_description_text' && "$bc" != *"--manifest"* ]] && ok "the comparator's built-in classes need no manifest beside the input" || bad "built-in masks: rc=$rc $bc"
[[ "$bc" != *" --input "* ]] && ok "a default comparison without a manifest or an input-reading class passes no --input" || bad "input leaked without a manifest: $bc"
expect_error "a manifest class beside a built-in one still needs the manifest" "masks names 'some_class', which applies from a manifest.json" "$T3" --only c-mixed
# The runner's built-in list is the comparator's: the row and text classes, the registry pair
# classes, the field classes, the transcript row classes and the warnings-mode classes, no more
# and no fewer, but for one field class. dbnsfp_mitochondrial_record_named_m reads the assembly
# of a manifest and accepts nothing without one, so it applies beside a manifest.json alone and
# stands outside the list: the runner refuses a case that names it with no manifest beside its
# input.
builtin_columns="TERM_PAIR_CLASSES REGISTRY_CLASSES FIELD_CLASSES STRUCTURAL_VARIANT_CLASSES PLUGIN_WARNING_CLASSES WARNING_TEXT_CLASSES"
builtin_listed=$(awk -F'"' '/^BUILTIN_MASKS="/ {print $2; exit}' "$ORACLE" | tr ' ' '\n' | sort)
builtin_implemented=$({
    printf '%s\n' regulatory_feature_rows order_only individual_allele_order conservation_short_name_spelling \
        pc_description_text motif_score_reverse_strand_mnv merged_same_stable_id_collapse \
        colocated_list_cut_at_cache_cap
    class_table_members "$builtin_columns" | grep -vx dbnsfp_mitochondrial_record_named_m
} | sort -u)
[[ "$builtin_listed" == "$builtin_implemented" ]] && ok "BUILTIN_MASKS names exactly the comparator's built-in classes" || bad "BUILTIN_MASKS drift: listed [$(tr '\n' ' ' <<<"$builtin_listed")] implemented [$(tr '\n' ' ' <<<"$builtin_implemented")]"
# The comment above the runner's list says what the list holds: it cites the columns of the
# class file the list is held to, and it names each class of those columns that the list
# leaves out. A cited column is a whole word that ends _CLASSES; a word of another shape is
# not read as one, so the comment may name a column of the class file that holds no class list.
builtin_comment=$(awk '/^BUILTIN_MASKS="/ {print block; exit} /^#/ {block = block " " $0; next} {block = ""}' "$ORACLE")
comment_columns=$(grep -owE '[A-Z]+(_[A-Z]+)*_CLASSES' <<<"$builtin_comment" | sort -u)
comment_words=$(tr -c 'A-Za-z0-9_' '\n' <<<"$builtin_comment" | sort -u)
comment_unnamed=""
while read -r table_class; do
    [[ -n "$table_class" ]] || continue
    grep -x "$table_class" <<<"$builtin_listed" >/dev/null && continue
    grep -x "$table_class" <<<"$comment_words" >/dev/null || comment_unnamed="$comment_unnamed $table_class"
done < <(class_table_members "$builtin_columns")
[[ "$comment_columns" == "$(tr ' ' '\n' <<<"$builtin_columns" | sort -u)" && -z "$comment_unnamed" ]] &&
    ok "the comment above BUILTIN_MASKS cites the list's columns and names each class of them the list leaves out" || bad "the comment above BUILTIN_MASKS: columns cited [$(tr '\n' ' ' <<<"$comment_columns")], classes outside the list it does not name [$comment_unnamed ]"
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
    row c-alleles fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 revel_sibling_allele_rows - - - - - - -
    row c-span fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 splice_term_past_the_record_span - - - - - - -
    row c-gwas fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv --plugin GWAS,file={corpus}/catalog.tsv" default "$IMAGE" offline 116.2 gwas_haplotype_cross_product=2 - - - - - - -
    row c-gwas-joined fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin=GWAS,type=curated,file={corpus}/catalog.tsv" tab "$IMAGE" offline 116.2 gwas_interaction_first_snp_allele_whitespace - - - - - - -
    row c-gwas-unmasked fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin GWAS,file={corpus}/catalog.tsv" default "$IMAGE" offline 116.2 order_only - - - - - - -
    row c-gwas-fileless fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv" default "$IMAGE" offline 116.2 gwas_haplotype_cross_product - - - - - - -
    row c-clipped fixture nomanifest/variants.vcf "$CORPUS_CACHE" - --hgvs default "$IMAGE" offline 116.2 hgvs_clipped_reference_against_whole_alternate - - - - - - -
    row c-shifted-protein fixture nomanifest/variants.vcf "$CORPUS_CACHE" - --hgvs default "$IMAGE" offline 116.2 hgvsp_shifted_indel_written_span_gate - - - - - - -
    row c-downstream fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin Downstream" default,tab "$IMAGE" offline 116.2 downstream_utr_frame_past_coding_end,downstream_start_past_own_codon - - - - - - -
    row c-organelle fixture nomanifest/variants.vcf "$CORPUS_CACHE" - --hgvs default "$IMAGE" offline 116.2 organelle_peptide_paths_translation_table - - - - - - -
    row c-edge fixture nomanifest/variants.vcf 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 hgvs_shift_sequence_edge_window - - - - - - -
    row c-edge-manifest fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 hgvs_shift_sequence_edge_window - - - - - - -
    row c-intron-end fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--hgvs" default "$IMAGE" offline 116.2 hgvsp_intron_end_insertion - - - - - - -
    row c-intron-start fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--hgvs" default "$IMAGE" offline 116.2 hgvsp_intron_start_insertion - - - - - - -
    row c-two-windows fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 maxentscan_overlapping_windows - - - - - - -
    row c-inserted-span fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 refseq_edited_model_inserted_span_codons - - - - - - -
} >"$TB"
dry "$TB" "$TMP/d-beside.out"
rc=$?
input_listed=$(awk -F'"' '/^INPUT_READING_MASKS="/ {print $2; exit}' "$ORACLE" | tr ' ' '\n' | sort)
input_implemented=$({
    echo reference_skipped_record
    class_table_members ALLELE_READING_CLASSES RECORD_REFERENCE_CLASSES
} | sort -u)
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
# The comparator builds a side table from the candidate's JSON cache for the classes its
# data/json_cache_classes.tsv names, a class a line; of the other classes its masks name it
# opens a cache it is given for reference_skipped_record alone.
# The runner hands the cache to the classes of CACHE_READING_MASKS and, beside the FASTA, to
# those of REFERENCE_READING_MASKS, each list read as the runner sets it, a list a later line
# adds to with it. Every class of the file stands in one of the two lists, and every name of
# the two is a class of the file, but for a class of a score-reading list where the runner
# has one: the comparator's score table reads the cache on --spliceai-table alone and its
# class has no line in the file.
# The lists are read from the runner's text: a line that sets one of the three is read where
# the whole line is NAME="..." or NAME+="..." and its value holds class names and the three
# lists alone. A line that sets one of them with = or += in any other form is not read and
# fails the check.
cache_implemented=$(awk -F'\t' 'NR > 1 {print $1}' "$TOOL/data/json_cache_classes.tsv" | sort)
cache_line='(CACHE|REFERENCE|SCORE)_READING_MASKS\+?="([a-z0-9_ ]|\$\{?(CACHE|REFERENCE|SCORE)_READING_MASKS\}?)*"'
cache_lists=$(grep -E "^$cache_line\$" "$ORACLE")
cache_other_form=$(grep -E '(^|[^A-Za-z0-9_])(CACHE|REFERENCE|SCORE)_READING_MASKS\+?=' "$ORACLE" | grep -cvE "^$cache_line\$")
cache_handed=$(
    eval "$cache_lists"
    printf '%s\n' ${CACHE_READING_MASKS:-} ${REFERENCE_READING_MASKS:-} | sort -u
)
cache_score=$(
    eval "$cache_lists"
    echo "${SCORE_READING_MASKS:-}"
)
cache_unlisted=""
for cache_class in $cache_implemented; do
    grep -Fx "$cache_class" <<<"$cache_handed" >/dev/null || cache_unlisted="$cache_unlisted $cache_class"
done
cache_unread=""
for cache_class in $cache_handed; do
    grep -Fx "$cache_class" <<<"$cache_implemented" >/dev/null && continue
    [[ " $cache_score " == *" $cache_class "* ]] || cache_unread="$cache_unread $cache_class"
done
[[ -n "$cache_implemented" && -n "$cache_handed" && "$cache_other_form" == 0 && -z "$cache_unlisted" && -z "$cache_unread" ]] &&
    ok "CACHE_READING_MASKS and REFERENCE_READING_MASKS name exactly the classes the comparator builds a side table from the candidate's JSON cache for" || bad "cache-reading lists: classes of data/json_cache_classes.tsv outside both lists [$cache_unlisted ], listed names the file lacks [$cache_unread ], lines that set a list in a form the check does not read: $cache_other_form"
# Beside a manifest the comparator is told the candidate's run: its run of the manifest, and
# whether it has a reference sequence, which the candidate takes from --fasta or from an
# indexed file inside its cache directory. A fabricated corpus with its own run, two flag
# sets and a manifest side that names --grch37 and a --fasta word, and candidate caches that
# hold each form of reference file.
RS="$DATA/runsetting"
mkdir -p "$RS/corpus" "$RS/bare/transcripts" "$RS/plain/transcripts" "$RS/unindexed/transcripts" "$RS/bgzf/transcripts" \
    "$RS/bgzf-no-gzi/transcripts" "$RS/gzip/transcripts" "$RS/other-name/transcripts" "$RS/dot-name/transcripts" "$RS/dot-name-bgzf/transcripts"
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
# Names that open with a dot, which a listing of the directory's entries holds: a .fa behind
# a longer name, and a BGZF file whose whole name is .fa.gz.
printf '>21\nA\n' >"$RS/dot-name/.ref.fa"
cp "$RS/ref.fa.fai" "$RS/dot-name/.ref.fa.fai"
bgzf_head >"$RS/dot-name-bgzf/.fa.gz"
cp "$RS/ref.fa.fai" "$RS/dot-name-bgzf/.fa.gz.fai"
: >"$RS/dot-name-bgzf/.fa.gz.gzi"
T_RS="$TMP/t-runsetting.tsv"
rs_row() { # rs_row <case> <input> <cache dir under runsetting/> <fasta cell> <flags cell>
    row "$1" fixture "runsetting/corpus/$2" "116_GRCh38=runsetting/$3" "$4" "$5" default "$IMAGE" offline 116.2 - - - - - - - -
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
    rs_row c-rs-dot-name variants.vcf dot-name - "--hgvs --symbol"
    rs_row c-rs-dot-name-bgzf variants.vcf dot-name-bgzf - "--hgvs --symbol"
    rs_row c-rs-no-fasta variants.vcf plain - "--symbol -NO_FASTA"
    rs_row c-rs-fasta-dir variants.vcf plain - "--symbol --fasta_dir {corpus}"
    row c-rs-no-manifest fixture nomanifest/variants.vcf "116_GRCh38=runsetting/plain" - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T_RS"
# The manifest lists four runs over variants.vcf, so a case whose flags are not exactly one of
# them has a row of the named-runs table: a named run where the flags fit two, pooled where
# no run carries them.
RS_RUNS="$TMP/rs-runs.tsv"
{
    row case run reason
    row c-rs-none pooled "no run of the manifest carries --symbol alone"
    row c-rs-twin twin-a "the flags fit twin-a and twin-b alike, and the case gives no FASTA"
    row c-rs-cell-no-fasta pooled "no run of the manifest carries --no_fasta"
    row c-rs-no-fasta pooled "no run of the manifest carries --no_fasta"
    row c-rs-fasta-dir pooled "no run of the manifest carries --fasta_dir"
} >"$RS_RUNS"
dry "$T_RS" "$TMP/d-rs.out" --case-runs "$RS_RUNS"
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
    ok "a case that is no run of its manifest and has a pooled row is told no run" || bad "no run: $(rs_cmd c-rs-none)"
rs_told c-rs-twin "--run twin-a --reference-sequence absent" &&
    ok "a case that fits two runs of its manifest is told the run its row names" || bad "two runs: $(rs_cmd c-rs-twin)"
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
rs_told c-rs-dot-name "--run default --reference-sequence present" && rs_told c-rs-dot-name-bgzf "--run default --reference-sequence present" &&
    ok "a .fa or .fa.gz whose name opens with a dot is a reference sequence with its index files beside it: the candidate lists every entry of its cache directory" || bad "cache dot name: $(rs_cmd c-rs-dot-name) | BGZF: $(rs_cmd c-rs-dot-name-bgzf)"
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
    row c-edge-bare fixture nomanifest/variants.vcf "$CORPUS_CACHE" - --hgvs default "$IMAGE" offline 116.2 hgvs_shift_sequence_edge_window - - - - - - -
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
# hgvs_insertion_left_where_the_shift_stops holds a duplication to the exon spans of the
# candidate's cache. Its case is the transcript-edge case's row under its own name and class,
# whatever columns the table has.
TLW="$TMP/tlw.tsv"
awk -F'\t' -v OFS='\t' 'NR == 1 {print} $1 == "c-clipped" {$1 = "c-left-where"; sub(/hgvs_clipped_reference_against_whole_alternate/, "hgvs_insertion_left_where_the_shift_stops"); print}' "$TB" >"$TLW"
dry "$TLW" "$TMP/d-left-where.out"
lwc=$(cmd_second 'DRY-RUN c-left-where stub default' "$TMP/d-left-where.out")
[[ "$lwc" == *"--masks hgvs_insertion_left_where_the_shift_stops --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming hgvs_insertion_left_where_the_shift_stops carries the candidate's JSON cache as --json-cache" || bad "left-where cache: $lwc"
# The class of an allele over a transcript's edge counts the exon spans of the candidate's cache.
kc=$(cmd_second 'DRY-RUN c-clipped stub default' "$TMP/d-beside.out")
[[ "$kc" == *"--masks hgvs_clipped_reference_against_whole_alternate --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming the class of an allele over a transcript's edge carries the candidate's cache as --json-cache" || bad "transcript-edge cache: $kc"
# The protein notation of a shifted indel is held to the exon spans of the candidate's cache
# where the shift is longer than the coding coordinate the notation gives.
sp=$(cmd_second 'DRY-RUN c-shifted-protein stub default' "$TMP/d-beside.out")
[[ "$sp" == *"--masks hgvsp_shifted_indel_written_span_gate --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming the protein notation of a shifted indel carries the candidate's cache as --json-cache" || bad "shifted-indel protein cache: $sp"
# hgvsp_stop_codon_protein_comparison reads a transcript's coding sequence from the candidate's
# cache where more than one residue stands before an inserted stop: the two-window case, named
# under that class, carries the cache.
TK="$TMP/tk.tsv"
{
    head -n 1 "$TB"
    awk -F'\t' -v OFS='\t' '$1 == "c-two-windows" {
        $1 = "c-kept-protein"
        for (i = 2; i <= NF; i++) if ($i == "maxentscan_overlapping_windows") $i = "hgvsp_stop_codon_protein_comparison"
        print
    }' "$TB"
} >"$TK"
dry "$TK" "$TMP/d-kept.out"
kp=$(cmd_second 'DRY-RUN c-kept-protein stub default' "$TMP/d-kept.out")
[[ "$kp" == *"--masks hgvsp_stop_codon_protein_comparison --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming hgvsp_stop_codon_protein_comparison carries the candidate's JSON cache as --json-cache" || bad "kept-protein cache: $kp"
# hgvsp_frameshift_count_past_a_stop_read_as_a_residue derives each count from the candidate's
# cache: the two-window case, named under that class, carries the cache.
TF="$TMP/tf.tsv"
{
    head -n 1 "$TB"
    awk -F'\t' -v OFS='\t' '$1 == "c-two-windows" {
        $1 = "c-counted-frameshift"
        for (i = 2; i <= NF; i++) if ($i == "maxentscan_overlapping_windows") $i = "hgvsp_frameshift_count_past_a_stop_read_as_a_residue"
        print
    }' "$TB"
} >"$TF"
dry "$TF" "$TMP/d-counted.out"
cf=$(cmd_second 'DRY-RUN c-counted-frameshift stub default' "$TMP/d-counted.out")
[[ "$cf" == *"--masks hgvsp_frameshift_count_past_a_stop_read_as_a_residue --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming hgvsp_frameshift_count_past_a_stop_read_as_a_residue carries the candidate's JSON cache as --json-cache" || bad "counted-frameshift cache: $cf"
# refseq_failed_edit_model_hgvsc_offset reads a transcript's edits and coding positions from
# the candidate's cache: the two-window case, named under that class, carries the cache.
TO="$TMP/to.tsv"
{
    head -n 1 "$TB"
    awk -F'\t' -v OFS='\t' '$1 == "c-two-windows" {
        $1 = "c-position-past-an-edit"
        for (i = 2; i <= NF; i++) if ($i == "maxentscan_overlapping_windows") $i = "refseq_failed_edit_model_hgvsc_offset"
        print
    }' "$TB"
} >"$TO"
dry "$TO" "$TMP/d-offset.out"
po=$(cmd_second 'DRY-RUN c-position-past-an-edit stub default' "$TMP/d-offset.out")
[[ "$po" == *"--masks refseq_failed_edit_model_hgvsc_offset --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming refseq_failed_edit_model_hgvsc_offset carries the candidate's JSON cache as --json-cache" || bad "position-past-an-edit cache: $po"
# hgvsp_frameshift_count_at_the_first_residue derives each count from the candidate's cache
# and reads the length of the record's reference bases from the input: the two-window case,
# named under that class, carries both.
TFR="$TMP/tfr.tsv"
{
    head -n 1 "$TB"
    awk -F'\t' -v OFS='\t' '$1 == "c-two-windows" {
        $1 = "c-counted-at-the-first-residue"
        for (i = 2; i <= NF; i++) if ($i == "maxentscan_overlapping_windows") $i = "hgvsp_frameshift_count_at_the_first_residue"
        print
    }' "$TB"
} >"$TFR"
dry "$TFR" "$TMP/d-counted-first.out"
cfr=$(cmd_second 'DRY-RUN c-counted-at-the-first-residue stub default' "$TMP/d-counted-first.out")
[[ "$cfr" == *" --input $DATA/nomanifest/variants.vcf --masks hgvsp_frameshift_count_at_the_first_residue --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming hgvsp_frameshift_count_at_the_first_residue carries the input and the candidate's JSON cache as --json-cache" || bad "first-residue count cache and input: $cfr"
# stop_codon_settled_before_unknown_bases reads a transcript's model and coding sequence from
# the candidate's cache for an insertion carried out of an intron: the two-window case, named
# under that class, carries the cache.
TSS="$TMP/tss.tsv"
{
    head -n 1 "$TB"
    awk -F'\t' -v OFS='\t' '$1 == "c-two-windows" {
        $1 = "c-settled-stop"
        for (i = 2; i <= NF; i++) if ($i == "maxentscan_overlapping_windows") $i = "stop_codon_settled_before_unknown_bases"
        print
    }' "$TB"
} >"$TSS"
dry "$TSS" "$TMP/d-settled-stop.out"
ss=$(cmd_second 'DRY-RUN c-settled-stop stub default' "$TMP/d-settled-stop.out")
[[ "$ss" == *"--masks stop_codon_settled_before_unknown_bases --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming stop_codon_settled_before_unknown_bases carries the candidate's JSON cache as --json-cache" || bad "settled-stop cache: $ss"
ic=$(cmd_second 'DRY-RUN c-intron-end stub default' "$TMP/d-beside.out")
[[ "$ic" == *"--masks hgvsp_intron_end_insertion --json-cache $CORPUS/json_cache" && "$ac" != *"--json-cache"* ]] &&
    ok "a comparison naming hgvsp_intron_end_insertion carries the candidate's JSON cache as --json-cache, one naming no cache-reading class none" || bad "intron-end cache: named=$ic unnamed=$ac"
isc=$(cmd_second 'DRY-RUN c-intron-start stub default' "$TMP/d-beside.out")
[[ "$isc" == *"--masks hgvsp_intron_start_insertion --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming hgvsp_intron_start_insertion carries the candidate's JSON cache as --json-cache" || bad "intron-start cache: $isc"
# The class of a stop named on an incomplete coding sequence reads from the candidate's cache
# the transcripts whose coding sequence is no whole number of codons. Its case is the organelle
# case's row under its own name and class, whatever columns the table has.
TSI="$TMP/tsi.tsv"
awk -F'\t' -v OFS='\t' 'NR == 1 {print} $1 == "c-organelle" {$1 = "c-stop-incomplete"; sub(/organelle_peptide_paths_translation_table/, "hgvsp_stop_named_on_incomplete_coding_sequence"); print}' "$TB" >"$TSI"
dry "$TSI" "$TMP/d-stop-incomplete.out"
snc=$(cmd_second 'DRY-RUN c-stop-incomplete stub default' "$TMP/d-stop-incomplete.out")
[[ "$snc" == *"--masks hgvsp_stop_named_on_incomplete_coding_sequence --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming hgvsp_stop_named_on_incomplete_coding_sequence carries the candidate's JSON cache as --json-cache" || bad "incomplete coding sequence cache: $snc"
# maxentscan_overlapping_windows reads a transcript's introns from the candidate's cache.
tw=$(cmd_second 'DRY-RUN c-two-windows stub default' "$TMP/d-beside.out")
[[ "$tw" == *"--masks maxentscan_overlapping_windows --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming maxentscan_overlapping_windows carries the candidate's JSON cache as --json-cache" || bad "two-window cache: $tw"
# refseq_edited_model_inserted_span_codons reads a transcript's coding sequence from the candidate's cache.
isp=$(cmd_second 'DRY-RUN c-inserted-span stub default' "$TMP/d-beside.out")
[[ "$isp" == *"--masks refseq_edited_model_inserted_span_codons --json-cache $CORPUS/json_cache" ]] &&
    ok "a comparison naming refseq_edited_model_inserted_span_codons carries the candidate's JSON cache as --json-cache" || bad "inserted-span cache: $isp"
# The five classes of a coding end read from the candidate's cache what stands at a
# transcript's coding ends and take nothing without it. The case of each is the organelle
# case's row under the class's own name: it names no FASTA and stands beside no manifest.
for coding_end_class in structural_deletion_stop_lost_without_stop_codon \
    stop_codon_completed_by_polyadenylation stop_codon_overlap_without_stop_codon \
    stop_codon_behind_last_residue_without_stop_codon start_codon_overlap_without_start_codon; do
    TCE="$TMP/tce.tsv"
    awk -F'\t' -v OFS='\t' -v class="$coding_end_class" 'NR == 1 {print} $1 == "c-organelle" {$1 = "c-coding-end"; sub(/organelle_peptide_paths_translation_table/, class); print}' "$TB" >"$TCE"
    dry "$TCE" "$TMP/d-coding-end.out"
    cec=$(cmd_second 'DRY-RUN c-coding-end stub default' "$TMP/d-coding-end.out")
    [[ "$cec" == *"--masks $coding_end_class --json-cache $CORPUS/json_cache" ]] &&
        ok "a comparison naming $coding_end_class carries the candidate's JSON cache as --json-cache" || bad "coding end cache ($coding_end_class): $cec"
done
# The class of a stop codon that polyadenylation completes reads the reference bases past a
# transcript's end for one of its shapes: its comparison carries the case's FASTA where the
# case names one, beside the cache and the input record it reads in every case.
TCF="$TMP/tcf.tsv"
awk -F'\t' -v OFS='\t' 'NR == 1 {print} $1 == "c-edge" {$1 = "c-completed-stop"; sub(/hgvs_shift_sequence_edge_window/, "stop_codon_completed_by_polyadenylation"); print}' "$TB" >"$TCF"
dry "$TCF" "$TMP/d-completed-stop.out"
csf=$(cmd_second 'DRY-RUN c-completed-stop stub default' "$TMP/d-completed-stop.out")
reference_all=$(awk -F'"' '/^REFERENCE_READING_MASKS(\+)?="/ {print $2}' "$ORACLE" | tr -s ' \n' ' ')
[[ "$reference_all" == "hgvs_shift_sequence_edge_window stop_codon_completed_by_polyadenylation " && "$csf" == *" --input $DATA/nomanifest/variants.vcf --masks stop_codon_completed_by_polyadenylation --json-cache $REPO/tests/golden/116/GRCh37-hgvs/json_cache --reference-fasta $TMP/never/work/fasta/1/reference.fa" ]] &&
    ok "a comparison naming stop_codon_completed_by_polyadenylation beside a FASTA carries the input, the candidate's cache and the FASTA as --reference-fasta" || bad "completed stop codon reference: listed [$reference_all] default=$csf"
# A class of a stop codon or a start codon that stands is held, beside the cache, to a
# transcript that shows the codon: a case naming one carries the cache.
landed_without_cache=""
for landed_class in start_cooccurrence_swap stop_codon_unshifted_span stop_retained_without_stop_codon_test \
    stop_codon_overlap_genomic_span stop_retained_three_base_stop_codon_test \
    stop_codon_removed_without_stop_codon_test stop_codon_overlap_insertion_range; do
    TLC="$TMP/tlc.tsv"
    awk -F'\t' -v OFS='\t' -v class="$landed_class" 'NR == 1 {print} $1 == "c-organelle" {$1 = "c-landed"; sub(/organelle_peptide_paths_translation_table/, class); print}' "$TB" >"$TLC"
    dry "$TLC" "$TMP/d-landed.out"
    lcc=$(cmd_second 'DRY-RUN c-landed stub default' "$TMP/d-landed.out")
    [[ "$lcc" == *"--masks $landed_class --json-cache $CORPUS/json_cache" ]] || landed_without_cache="$landed_without_cache $landed_class"
done
[[ -z "$landed_without_cache" ]] &&
    ok "a comparison naming a class of a stop codon or a start codon that stands carries the candidate's JSON cache as --json-cache" || bad "no cache beside:$landed_without_cache"
# A class may carry its expected applied count: the count stays with the runner's report check
# and the comparator receives the class names alone; an item that is neither form is refused.
{
    echo "$HEADER"
    row c-counted fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 start_cooccurrence_swap=127,splice_family_swap,conservation_max_text_order=3 - - - - - - -
    row c-badcount fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=x - - - - - - -
    row c-doubleeq fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=1=2 - - - - - - -
    row c-zero fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=0 - - - - - - -
    row c-empty fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,,regulatory_feature_rows - - - - - - -
    row c-trail fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=1, - - - - - - -
    row c-blank fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 "order_only, regulatory_feature_rows" - - - - - - -
    row c-glob fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 "*" - - - - - - -
    row c-twice fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only=1,regulatory_feature_rows,order_only=3 - - - - - - -
    row c-countedw fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default,warnings "$IMAGE" offline 116.2 start_cooccurrence_swap=127,dosage_sensitivity_regulatory_row_warnings - - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-counted.out" --only c-counted
rc=$?
bc=$(cmd_second 'DRY-RUN c-counted stub default' "$TMP/d-counted.out")
[[ "$rc" -eq 0 && "$bc" == *'--masks start_cooccurrence_swap\,splice_family_swap\,conservation_max_text_order'" --json-cache $CORPUS/json_cache" && "$bc" != *"=127"* ]] && ok "a masks item's expected count is stripped before the comparator sees the classes, and the cache follows a class that reads it" || bad "counted masks: rc=$rc $bc"
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
    row c-dash fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 -,order_only - - - - - - -
} >"$T3"
expect_error "a - beside a class in the masks cell is an error (the comparator would read it as an option)" "masks '-,order_only' has an item beginning with -" "$T3"
{
    echo "$HEADER"
    row c-dashfloor fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - -,IMPACT=1 - - - - - -
    row c-dashpair fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - "-;IMPACT=LOW=>HIGH" - - - - -
    row c-dashclass fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 -foo - - - - - - -
    row c-dashafter fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,-foo - - - - - - -
} >"$T3"
expect_error "a - item in the floors cell is an error" "floors '-,IMPACT=1' has an item beginning with -" "$T3" --only c-dashfloor
expect_error "a - item in the field_pairs cell is an error" "field_pairs '-;IMPACT=LOW=>HIGH' has an item beginning with -" "$T3" --only c-dashpair
expect_error "a masks class beginning with - is an error even beside a manifest" "masks '-foo' has an item beginning with -" "$T3" --only c-dashclass
expect_error "a masks class beginning with - after a separator is an error" "masks 'order_only,-foo' has an item beginning with -" "$T3" --only c-dashafter
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json,warnings "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-floors fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,json "$IMAGE" offline 116.2 - Uploaded_variation=1,IMPACT=1 - - - - - -
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
    row c-jsonfloor fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" json "$IMAGE" offline 116.2 - Uploaded_variation=1,IMPACT=1 - - - - - -
} >"$T4F"
run_oracle "$T4F" --out "$TMP/out-jsonfloor" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(reasons_of "$TMP/out-jsonfloor/abc1234/c-jsonfloor/report.json")" == "json: floor Uploaded_variation is not applicable to json and no other format of the case counted it" ]] &&
    ok "a json-only case with a floor the json format cannot count is red: no format checked the floor" || bad "json-only floor: rc=$rc $(reasons_of "$TMP/out-jsonfloor/abc1234/c-jsonfloor/report.json" 2>/dev/null)"
# Extra is the default format's column alone: the tab header does not declare it and the json
# format has no row count for it, so a tab,json case sets the floor aside in both formats.
{
    echo "$HEADER"
    row c-twoaside fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" tab,json "$IMAGE" offline 116.2 - Extra=1,IMPACT=1 - - - - - -
} >"$T4F"
run_oracle "$T4F" --out "$TMP/out-twoaside" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(reasons_of "$TMP/out-twoaside/abc1234/c-twoaside/report.json")" == "json: floor Extra is not applicable to json and no other format of the case counted it | tab: floor Extra is not applicable to tab and no other format of the case counted it" ]] &&
    ok "a floor every compared format sets aside is red even when a text format was compared" || bad "two formats set the floor aside: rc=$rc $(reasons_of "$TMP/out-twoaside/abc1234/c-twoaside/report.json" 2>/dev/null)"

echo
echo "TERMINAL: a run whose standard error is a terminal reads the case table as a run without one does"
# From bash 5.0 a script's shell sets LINES and COLUMNS to the size of the terminal on its
# standard error after each command it waits for, so a list the runner kept under either
# name would be a number by the line that reads it.
# on_terminal.py <rows> <columns> <capture file> <command...>: the command with its standard
# error on a pseudo-terminal of that size; what it writes there goes to the capture file, and
# its exit code is the helper's.
cat >"$TMP/on_terminal.py" <<'PY'
import fcntl
import os
import struct
import subprocess
import sys
import termios

rows, columns, capture = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
master, slave = os.openpty()
fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", rows, columns, 0, 0))
child = subprocess.Popen(sys.argv[4:], stderr=slave)
os.close(slave)
with open(capture, "wb") as out:
    while True:
        try:
            data = os.read(master, 65536)
        except OSError:
            break
        if not data:
            break
        out.write(data)
sys.exit(child.wait())
PY
# Where no command starts on a pseudo-terminal of that size, as on a machine without the device,
# no run on a terminal can be read: the checks are left out behind one skip line.
size=$(python3 "$TMP/on_terminal.py" 40 120 "$TMP/term-size.err" python3 -c 'import os; size = os.get_terminal_size(2); print(size.lines, size.columns)') || size=""
if [[ "$size" == "40 120" ]]; then
    TT="$TMP/t-term.tsv"
    echo "$HEADER" >"$TT"
    STUB_DOCKER_LOG="" STUB_CAND_LOG="" run_oracle "$TT" --out "$TMP/out-term-none" --candidates "$STUBDIR/vep" --candidate-pins stub >"$TMP/term-none.out" 2>"$TMP/term-none.err"
    rc_none=$?
    STUB_DOCKER_LOG="" STUB_CAND_LOG="" PATH="$STUBDIR:$PATH" python3 "$TMP/on_terminal.py" 40 120 "$TMP/term-on.err" \
        bash "$ORACLE" --cases "$TT" --data-dir "$DATA" --out "$TMP/out-term-on" --candidates "$STUBDIR/vep" --candidate-pins stub >"$TMP/term-on.out"
    rc_on=$?
    [[ "$rc_none" -eq 1 && "$rc_on" -eq 1 ]] && grep -qF "no case of tier fixture selected from $TT" "$TMP/term-none.err" && grep -qF "no case of tier fixture selected from $TT" "$TMP/term-on.err" &&
        ! grep -qF "header lacks the column" "$TMP/term-on.err" && cmp -s "$TMP/term-none.out" "$TMP/term-on.out" &&
        ok "a real run over a table of a header alone passes the header check on a terminal of 40 rows and 120 columns and ends on the line and the exit code of a run without one" ||
        bad "a header-only table on a terminal: rc $rc_on (without one: $rc_none): $(tail -1 "$TMP/term-on.err")"
    STUB_DOCKER_LOG="" STUB_CAND_LOG="" PATH="$STUBDIR:$PATH" python3 "$TMP/on_terminal.py" 40 120 "$TMP/term-green.err" \
        bash "$ORACLE" --cases "$T4" --data-dir "$DATA" --out "$TMP/out-term-green" --candidates "$STUBDIR/vep" --candidate-pins abc1234 >"$TMP/term-green.out"
    rc=$?
    [[ "$rc" -eq 0 && "$(sed -n 2p "$TMP/out-term-green/status.tsv" 2>/dev/null | cut -f1-3)" == "$(row c-four abc1234 green)" ]] &&
        ok "a real run of a case whose outputs are identical is green with its standard error on that terminal" ||
        bad "a green case on a terminal: rc $rc: $(tail -1 "$TMP/term-green.err")"
else
    echo "  skip the runs on a terminal: no command starts on a pseudo-terminal of 40 rows and 120 columns on this machine"
fi

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
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - IMPACT=100000 - - - - - -
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
STUB_CAND_EXIT_AFTER=3 STUB_PERL_WARNINGS="$W" STUB_CAND_WARNINGS="$W" run_oracle "$T4" --out "$TMP/out-w4" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-w4/p1/c-four/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"default: candidate exited 3"* && "$reasons" == *"warnings: the candidate's default run exited 3" ]] && ok "a candidate default run that fails after writing its warnings is red on warnings with its exit status" || bad "warnings under candidate exit 3 after writing: rc=$rc reasons=$reasons"

echo
echo "WARNINGS OF A NAMED RUN: a cell naming warnings without default compares the runs of its one format"
JW="WARNING: line 12 skipped (21 30165800 del ...): Skipped record"
TJW="$TMP/tjw.tsv"
{
    echo "$HEADER"
    row c-jw fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" json,warnings "$IMAGE" offline 116.2 - - - - - - - -
    row c-vw fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" vcf,warnings "$IMAGE" offline 116.2 - - - - - - - -
    row c-djw fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,json,warnings "$IMAGE" offline 116.2 - - - - - - - -
} >"$TJW"
dry "$TJW" "$TMP/d-jw.out"
rc=$?
wc_=$(cmd_after 'DRY-RUN c-jw stub warnings' "$TMP/d-jw.out")
[[ "$rc" -eq 0 && "$wc_" == *"--warnings --expected $TMP/never/work/cases/c-jw/perl/json.jsonl_warnings.txt --actual $TMP/never/work/cases/c-jw/candidates/stub/json.jsonl_warnings.txt --report $TMP/never/stub/c-jw/warnings.json" ]] &&
    ok "json,warnings compares the JSON runs' _warnings.txt files" || bad "json,warnings: rc=$rc command: $wc_ $(head -2 "$TMP/d-jw.out.err" | tr '\n' ' ')"
n=$(count_lines '^DRY-RUN c-jw perl ' "$TMP/d-jw.out")
[[ "$n" -eq 1 ]] && ok "json,warnings plans one Ensembl VEP run, the JSON run" || bad "$n Ensembl VEP runs planned for c-jw, expected 1"
wc_=$(cmd_after 'DRY-RUN c-vw stub warnings' "$TMP/d-jw.out")
[[ "$wc_" == *"--warnings --expected $TMP/never/work/cases/c-vw/perl/vcf.vcf_warnings.txt --actual $TMP/never/work/cases/c-vw/candidates/stub/vcf.vcf_warnings.txt --report $TMP/never/stub/c-vw/warnings.json" ]] &&
    ok "vcf,warnings compares the VCF runs' _warnings.txt files" || bad "vcf,warnings command: $wc_"
wc_=$(cmd_after 'DRY-RUN c-djw stub warnings' "$TMP/d-jw.out")
[[ "$wc_" == *"--warnings --expected $TMP/never/work/cases/c-djw/perl/default.txt_warnings.txt --actual $TMP/never/work/cases/c-djw/candidates/stub/default.txt_warnings.txt --report $TMP/never/stub/c-djw/warnings.json" ]] &&
    ok "a cell naming default beside json compares the default runs' _warnings.txt files" || bad "default,json,warnings command: $wc_"
JW_RUN=(--candidates "$STUBDIR/vep" --candidate-pins p1)
STUB_WARNINGS_OF=json.jsonl STUB_PERL_WARNINGS="$JW" run_oracle "$TJW" --out "$TMP/out-jw1" "${JW_RUN[@]}" --only c-jw >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-jw1/p1/c-jw/report.json")
[[ "$rc" -eq 1 && "$reasons" == "warnings: warnings: 1 entries only in the expected file" ]] && ok "a warning of Ensembl VEP's JSON run that the candidate's JSON run lacks is red on warnings alone" || bad "json,warnings one-sided: rc=$rc reasons=$reasons"
STUB_WARNINGS_OF=json.jsonl STUB_PERL_WARNINGS="$JW" STUB_CAND_WARNINGS="$JW" run_oracle "$TJW" --out "$TMP/out-jw2" "${JW_RUN[@]}" --only c-jw >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(python3 -c 'import json, sys; print(",".join(sorted(json.load(open(sys.argv[1]))["formats"])))' "$TMP/out-jw2/p1/c-jw/report.json")" == "json,warnings" ]] &&
    ok "the same warning in both JSON runs is green, and the report holds json and warnings" || bad "json,warnings same warning: rc=$rc"
STUB_WARNINGS_OF=json.jsonl STUB_PERL_WARNINGS="$JW" run_oracle "$TJW" --out "$TMP/out-jw3" "${JW_RUN[@]}" --only c-djw >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 ]] && ok "a cell naming default reads the default runs' warnings alone: a warning of the JSON run is not compared" || bad "default,json,warnings with a JSON-run warning exited $rc: $(reasons_of "$TMP/out-jw3/p1/c-djw/report.json")"
STUB_WARNINGS_OF=json.jsonl STUB_CAND_EXIT_AFTER=3 STUB_PERL_WARNINGS="$JW" STUB_CAND_WARNINGS="$JW" run_oracle "$TJW" --out "$TMP/out-jw4" "${JW_RUN[@]}" --only c-jw >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-jw4/p1/c-jw/report.json")
[[ "$rc" -eq 1 && "$reasons" == "json: candidate exited 3"*" | warnings: the candidate's json run exited 3" ]] && ok "a candidate JSON run that fails after writing its warnings is red on warnings with its exit status" || bad "json,warnings under candidate exit 3 after writing: rc=$rc reasons=$reasons"
STUB_PERL_RC=3 run_oracle "$TJW" --out "$TMP/out-jw5" "${JW_RUN[@]}" --only c-jw >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-jw5/p1/c-jw/report.json")
[[ "$rc" -eq 1 && "$reasons" == "json: Ensembl VEP exited 3"*" | warnings: Ensembl VEP exited 3"* ]] && ok "an Ensembl VEP JSON run that exits non-zero is red for json and for warnings" || bad "json,warnings under Ensembl VEP exit 3: rc=$rc reasons=$reasons"
# The class of the warnings file of a JSON run on a RefSeq cache, on the committed logs of
# tests/golden/116/GRCh38-refseq: Ensembl VEP's JSON run writes the mismatch warning with no
# record name, and the candidate's lines are the named ones both default runs write.
REFSEQ_LOGS="$REPO/tests/golden/116/GRCh38-refseq/expected"
UNNAMED=$(<"$REFSEQ_LOGS/json.warnings.log")
NAMED=$(sed -E 's/^WARNING: [0-9]+ : //' "$REFSEQ_LOGS/default.warnings.log")
{
    echo "$HEADER"
    row c-jw-class fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" json,warnings "$IMAGE" offline 116.2 fork_warning_pid_prefix,json_mismatch_warning_names_no_record - - - - - - -
    row c-jw-fork fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" json,warnings "$IMAGE" offline 116.2 fork_warning_pid_prefix - - - - - - -
    row c-dw-class fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,warnings "$IMAGE" offline 116.2 fork_warning_pid_prefix,json_mismatch_warning_names_no_record - - - - - - -
} >"$TJW"
STUB_WARNINGS_OF=json.jsonl STUB_PERL_WARNINGS="$UNNAMED" STUB_CAND_WARNINGS="$NAMED" run_oracle "$TJW" --out "$TMP/out-jw6" "${JW_RUN[@]}" --only c-jw-class,c-jw-fork >/dev/null 2>&1
rc=$?
python3 - "$TMP/out-jw6/p1" <<'PY' && [[ "$rc" -eq 1 ]] && ok "the JSON runs' warnings of a RefSeq corpus are green under json_mismatch_warning_names_no_record and red without it" || bad "the class over the JSON runs' warnings: rc=$rc"
import json, sys
root = sys.argv[1]
report = json.load(open(f"{root}/c-jw-class/report.json"))
assert report["verdict"] == "green", report["reasons"]
applied = json.load(open(f"{root}/c-jw-class/warnings.json"))["warnings"]["masks_applied"]
assert applied.get("json_mismatch_warning_names_no_record", 0) > 0, applied
assert applied.get("fork_warning_pid_prefix", 0) > 0, applied
report = json.load(open(f"{root}/c-jw-fork/report.json"))
assert report["verdict"] == "red", report
assert [r for r in report["reasons"] if r.startswith("warnings: ")], report["reasons"]
PY
STUB_PERL_WARNINGS="$(<"$REFSEQ_LOGS/default.warnings.log")" STUB_CAND_WARNINGS="$NAMED" run_oracle "$TJW" --out "$TMP/out-jw7" "${JW_RUN[@]}" --only c-dw-class >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-jw7/p1/c-dw-class/report.json")
[[ "$rc" -eq 1 && "$reasons" == *"json_mismatch_warning_names_no_record"* && "$reasons" != *"entries only in"* ]] && ok "control: over the default runs' warnings, which name each record, the class sets nothing aside and the case is red for the unused class alone" || bad "the class over the default runs' warnings: rc=$rc reasons=$reasons"
{
    echo "$HEADER"
    row c-two fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - tab,json,warnings "$IMAGE" offline 116.2 - - - - - - - -
    row c-pqw fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - parquet,warnings "$IMAGE" offline 116.2 - - - - - - - -
    row c-jw fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" json,warnings "$IMAGE" offline 116.2 - - - - - - - -
} >"$TJW"
expect_error "warnings beside two output formats and no default is an error" "c-two: formats names warnings without default beside 'tab json'" "$TJW" --only c-two
expect_error "warnings beside parquet alone is an error" "c-pqw: formats names warnings without default beside 'parquet'" "$TJW" --only c-pqw
# A stored reference of a json,warnings case holds the JSON run's warnings file.
REFJW="$TMP/refjw"
mkdir -p "$REFJW/c-jw/expected"
cp "$REF/json.jsonl" "$REFJW/c-jw/expected/json.jsonl"
run_oracle "$TJW" --out "$TMP/out-jw8" "${JW_RUN[@]}" --only c-jw --reference-dir "$REFJW" >/dev/null 2>"$TMP/jw8.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -qF "stored reference lacks json.jsonl_warnings.txt[.gz] under $REFJW/c-jw/expected/" "$TMP/jw8.err" && ok "a stored reference without the JSON run's warnings file fails a json,warnings case before anything runs" || bad "stored reference without json.jsonl_warnings.txt exited $rc: $(head -2 "$TMP/jw8.err" | tr '\n' ' ')"
printf '%s\n' "$JW" | gzip -c >"$REFJW/c-jw/expected/json.jsonl_warnings.txt.gz"
: >"$STUB_DOCKER_LOG"
STUB_WARNINGS_OF=json.jsonl STUB_CAND_WARNINGS="$JW" run_oracle "$TJW" --out "$TMP/out-jw9" "${JW_RUN[@]}" --only c-jw --reference-dir "$REFJW" >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(count_lines '^run ' "$STUB_DOCKER_LOG")" -eq 0 ]] && ok "a stored json.jsonl_warnings.txt.gz is the reference of a json,warnings case: green against the same warning, no container" || bad "stored JSON-run warnings: rc=$rc"
run_oracle "$TJW" --out "$TMP/out-jw10" "${JW_RUN[@]}" --only c-jw --reference-dir "$REFJW" >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-jw10/p1/c-jw/report.json")
[[ "$rc" -eq 1 && "$reasons" == "warnings: warnings: 1 entries only in the expected file" ]] && ok "a stored warning of the JSON run that the candidate's JSON run lacks is red" || bad "stored JSON-run warning one-sided: rc=$rc reasons=$reasons"

echo
echo "STATS: the two default runs' statistics reports through the statistics comparator"
# Two made-up reports of the grammar the tool reads, the second differing in one count. The
# stub comparator calls two files of one text green and any other pair red.
STATS_DIR="$TMP/stats"
mkdir -p "$STATS_DIR"
printf '%s\n' '[VEP run statistics]' "Species	homo_sapiens" '' '[General statistics]' "Lines of input read	36" "Overlapped genes	12" '' '[Variant classes]' "SNV	19" >"$STATS_DIR/perl.txt"
sed 's/^Overlapped genes	12$/Overlapped genes	11/' "$STATS_DIR/perl.txt" >"$STATS_DIR/cand-genes.txt"
cmp -s "$STATS_DIR/perl.txt" "$STATS_DIR/cand-genes.txt" && bad "the two statistics fixtures are equal" || ok "the second statistics fixture differs from the first in one count"
T4S="$TMP/t4s.tsv"
{
    echo "$HEADER"
    row c-stats fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,stats "$IMAGE" offline 116.2 - - - - - - - -
    row c-stats-w fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,warnings,stats "$IMAGE" offline 116.2 - - - - - - - -
} >"$T4S"
# stats_run <out name> <extra oracle args...>: the case c-stats on the real path; the caller sets the stubs' variables.
stats_run() {
    local out="$TMP/out-$1"
    shift
    run_oracle "$T4S" --out "$out" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-stats "$@" >"$out.log" 2>&1
}
: >"$TMP/stats-args.log"
STUB_STATS_LOG="$TMP/stats-args.log" STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/perl.txt" stats_run s1 --keep-outputs
rc=$?
[[ "$rc" -eq 0 ]] && ok "two default runs whose reports the comparator calls equal are green" || {
    bad "stats green run exited $rc: $(reasons_of "$TMP/out-s1/p1/c-stats/report.json" 2>/dev/null)"
    tail -3 "$TMP/out-s1.log" | sed 's/^/         /'
}
python3 - "$TMP/out-s1/p1/c-stats" <<'PY' && ok "report.json records stats as a format with the comparator's report, stats.json beside it" || bad "report.json stats entry"
import json, sys
r = json.load(open(sys.argv[1] + "/report.json"))
assert sorted(r["formats"]) == ["default", "stats"], sorted(r["formats"])
assert r["formats"]["stats"] == {"verdict": "green", "report": sys.argv[1] + "/stats.json", "reasons": []}, r["formats"]["stats"]
assert r["formats_requested"] == ["default", "stats"], r["formats_requested"]
sub = json.load(open(sys.argv[1] + "/stats.json"))
assert sub["mode"] == "stats" and sub["verdict"] == "green", sub
PY
[[ "$(<"$TMP/stats-args.log")" == "--expected $TMP/out-s1/work/cases/c-stats/perl/default.txt_summary.txt --actual $TMP/out-s1/work/cases/c-stats/candidates/p1/default.txt_summary.txt --report $TMP/out-s1/p1/c-stats/stats.json" ]] &&
    ok "the statistics comparator is started once, on the two reports the default runs left beside default.txt" || bad "statistics comparator arguments: $(<"$TMP/stats-args.log")"
grep -q -- ' --stats_text ' "$STUB_DOCKER_LOG" && ok "the Ensembl VEP default run of a stats case is started with --stats_text" || bad "no --stats_text reached the Ensembl VEP run"
STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/cand-genes.txt" stats_run s2
rc=$?
reasons=$(reasons_of "$TMP/out-s2/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "stats: [General statistics]: 1 difference" ]] && ok "a report the comparator calls red is red on stats alone, with the comparator's reason" || bad "stats red: rc=$rc reasons=$reasons"
STUB_STATS_RC=2 STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/perl.txt" stats_run s2e
rc=$?
reasons=$(reasons_of "$TMP/out-s2e/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "stats: comparator exited 2: ERROR: [compare_vep_stats] failing on request"* && ! -e "$TMP/out-s2e/p1/c-stats/stats.json" ]] && ok "a statistics comparator that exits 2 is a format error for stats, with its stderr" || bad "stats comparator exit 2: rc=$rc reasons=$reasons"
STUB_PERL_STATS="$STATS_DIR/perl.txt" stats_run s3
rc=$?
reasons=$(reasons_of "$TMP/out-s3/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "stats: the candidate exited 0 and wrote no default.txt_summary.txt" ]] && ok "a candidate that writes no statistics report is red on stats" || bad "stats missing candidate report: rc=$rc reasons=$reasons"
STUB_CAND_STATS="$STATS_DIR/perl.txt" stats_run s4
rc=$?
reasons=$(reasons_of "$TMP/out-s4/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "stats: Ensembl VEP exited 0 and wrote no default.txt_summary.txt" ]] && ok "an Ensembl VEP run that writes no statistics report is red on stats" || bad "stats missing Ensembl VEP report: rc=$rc reasons=$reasons"
STUB_PERL_RC=3 STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/perl.txt" stats_run s5
rc=$?
reasons=$(reasons_of "$TMP/out-s5/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: Ensembl VEP exited 3"*" | stats: Ensembl VEP exited 3"* ]] && ok "an Ensembl VEP default run that exits non-zero is red for stats with the run's error" || bad "stats under Ensembl VEP exit 3: rc=$rc reasons=$reasons"
STUB_CAND_EXIT_AFTER=3 STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/perl.txt" stats_run s6
rc=$?
reasons=$(reasons_of "$TMP/out-s6/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: candidate exited 3"*" | stats: the candidate's default run exited 3" ]] && ok "a candidate default run that fails after writing its report is red on stats with its exit status" || bad "stats under candidate exit 3 after writing: rc=$rc reasons=$reasons"
STUB_CAND_RC=4 STUB_PERL_STATS="$STATS_DIR/perl.txt" stats_run s7
rc=$?
reasons=$(reasons_of "$TMP/out-s7/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: candidate exited 4"*" | stats: the candidate's default run exited 4" ]] && ok "a candidate default run that fails before writing anything is red on stats with its exit status" || bad "stats under candidate exit 4: rc=$rc reasons=$reasons"
STUB_CAND_KILL_JOB=1 STUB_PERL_STATS="$STATS_DIR/perl.txt" stats_run s8 --format-jobs 2
rc=$?
reasons=$(reasons_of "$TMP/out-s8/p1/c-stats/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: the format's job ended without a result: no stderr | stats: the candidate's default run did not complete" ]] && ok "a default run whose format job died without an exit status did not complete, for stats" || bad "stats under a killed format job: rc=$rc reasons=$reasons"
STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/perl.txt" run_oracle "$T4S" --out "$TMP/out-s9" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-stats-w --format-jobs 3 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(python3 -c 'import json, sys; print(",".join(sorted(json.load(open(sys.argv[1]))["formats"])))' "$TMP/out-s9/p1/c-stats-w/report.json")" == "default,stats,tab,warnings" ]] && ok "stats beside warnings and a second format is green under --format-jobs 3" || bad "stats beside warnings: rc=$rc $(reasons_of "$TMP/out-s9/p1/c-stats-w/report.json" 2>/dev/null)"

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
# A pin is read as text: awk reads the pins 1 and 1.0 as one number, and they name two candidates.
run_oracle "$T4" --out "$TMP/out-numerals" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins 1,1.0 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(cut -f1,2,3 "$TMP/out-numerals/status.tsv" | tail -n +2 | tr '\n\t' '  ')" == "c-four 1 green c-four 1.0 green " ]] && ok "two candidates pinned 1 and 1.0 each take their own status row, and two green rows exit 0" || bad "numeral pins: rc=$rc $(cut -f1,2,3 "$TMP/out-numerals/status.tsv" | tail -n +2 | tr '\n\t' '  ')"

echo
echo "DATABASE MODE on the real path"
T6="$TMP/t6.tsv"
{
    echo "$HEADER"
    row c-db fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" database 116.2 - - - - - - - -
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
    row GRCh38-release116 fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json "$IMAGE" offline 116.2 - - - - - - - -
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
assert r["perl_image_digest"] == "ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3", r["perl_image_digest"]
assert r["reference_provenance_image"] == r["perl_image_digest"]
assert r["refusal_dir"] is None and r["candidate_dir_cache"] is None and r["unordered"] is False
assert r["verdict"] == "green"
PY
{
    echo "$HEADER"
    row GRCh38-release116 fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,warnings "$IMAGE" offline 116.2 - - - - - - - -
} >"$T7"
run_oracle "$T7" --out "$TMP/out-ref2" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >/dev/null 2>"$TMP/ref2.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'stored reference lacks default.txt_warnings.txt' "$TMP/ref2.err" && ok "a stored reference without a warnings file fails a warnings case before anything runs" || bad "stored reference without warnings exited $rc"
{
    echo "$HEADER"
    row no-such-corpus fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-stored fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json "$IMAGE" offline 116.2 - - - - - "$REF_CELL" - -
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
assert r["perl_image_digest"] == "ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3", r["perl_image_digest"]
assert r["verdict"] == "green", r["reasons"]
PY
dry "$TR" "$TMP/d-cell.out"
grep -q "^DRY-RUN c-stored perl stored: $CORPUS/expected/$" "$TMP/d-cell.out" && ! grep -q '^DRY-RUN c-stored perl default' "$TMP/d-cell.out" && ok "a dry run plans the case against its stored directory and no Ensembl VEP run" || bad "dry run of a reference cell: $(grep 'c-stored perl' "$TMP/d-cell.out")"
cc=$(cmd_second 'DRY-RUN c-stored stub vcf' "$TMP/d-cell.out")
[[ "$cc" == *"--format vcf --expected $CORPUS/expected/vcf.vcf.gz --actual "* ]] && ok "the comparison reads the stored file, gzipped as the corpus commits it" || bad "comparison of a reference cell: $cc"
{
    echo "$HEADER"
    row no-such-corpus fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - "$REF_CELL" - -
} >"$TR"
run_oracle "$TR" --out "$TMP/out-cell2" --candidates "$STUBDIR/vep" --candidate-pins p1 --reference-dir "$REPO/tests/golden/116" >/dev/null 2>"$TMP/cell2.err"
rc=$?
[[ "$rc" -eq 0 && ! -e "$STUB_DOCKER_CALLED" ]] && ok "beside --reference-dir, a case's reference cell is read before the directory under it" || bad "reference cell beside --reference-dir exited $rc: $(tail -2 "$TMP/cell2.err" | tr '\n' ' ')"
mkdir -p "$DATA/stored/partial"
cp "$REF/default.txt" "$DATA/stored/partial/default.txt"
{
    echo "$HEADER"
    row c-partial fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,vcf "$IMAGE" offline 116.2 - - - - - stored/partial - -
} >"$TR"
expect_error "a reference directory lacking a requested format fails preflight by name" "c-partial: stored reference lacks vcf.vcf\[.gz\] under $DATA/stored/partial/" "$TR"
{
    echo "$HEADER"
    row c-partial-stats fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,stats "$IMAGE" offline 116.2 - - - - - stored/partial - -
} >"$TR"
expect_error "a reference directory lacking the statistics report of a stats case fails preflight by name" "c-partial-stats: stored reference lacks default.summary.log\[.gz\] under $DATA/stored/partial/" "$TR"
# A stored default.summary.log, plain or gzipped, is the expected report of a stats case.
cp "$STATS_DIR/perl.txt" "$DATA/stored/stats/default.summary.log"
{
    echo "$HEADER"
    row c-stored-stats fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,stats "$IMAGE" offline 116.2 - - - - - stored/stats - -
} >"$TR"
STUB_CAND_STATS="$STATS_DIR/perl.txt" run_oracle "$TR" --out "$TMP/out-cell-stats" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>"$TMP/cell-stats.err"
rc=$?
[[ "$rc" -eq 0 && ! -e "$STUB_DOCKER_CALLED" ]] && grep -q '"verdict": "green"' "$TMP/out-cell-stats/p1/c-stored-stats/stats.json" && ok "a stored default.summary.log is the expected report of a stats case, with no container run" || {
    bad "stored stats reference exited $rc"
    tail -3 "$TMP/cell-stats.err" | sed 's/^/         /'
}
gzip -c "$DATA/stored/stats/default.summary.log" >"$DATA/stored/stats/default.summary.log.gz" && rm -f "$DATA/stored/stats/default.summary.log"
: >"$TMP/stats-args.log"
STUB_STATS_LOG="$TMP/stats-args.log" STUB_CAND_STATS="$STATS_DIR/perl.txt" run_oracle "$TR" --out "$TMP/out-cell-stats-gz" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(<"$TMP/stats-args.log")" == "--expected $DATA/stored/stats/default.summary.log.gz --actual "* ]] && ok "a gzipped stored report is handed to the comparator as the plain one is" || bad "gzipped stored stats reference exited $rc: $(reasons_of "$TMP/out-cell-stats-gz/p1/c-stored-stats/report.json" 2>/dev/null)"
{
    echo "$HEADER"
    row c-plain fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - stored/partial - -
} >"$TR"
dry "$TR" "$TMP/d-plain.out"
rc=$?
[[ "$rc" -eq 0 ]] && grep -q "^DRY-RUN c-plain perl stored: $DATA/stored/partial/$" "$TMP/d-plain.out" && ok "a relative reference cell resolves under --data-dir first, and a plain file serves as the .gz does" || bad "relative reference cell: rc=$rc"
{
    echo "$HEADER"
    row c-noref fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - "$DATA/stored/absent" - -
} >"$TR"
expect_error "a reference cell naming no directory is an error" "c-noref: stored reference directory not found: $DATA/stored/absent" "$TR"
[[ "$(grep -c 'c-noref' "$TMP/e.out.err")" -eq 1 ]] && ok "a missing reference directory is reported once, without a line per format" || bad "missing reference directory reported $(grep -c 'c-noref' "$TMP/e.out.err") times"
{
    echo "$HEADER"
    row c-badcell fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - yes - - -
} >"$TR"
expect_error "an unordered cell other than - or 1 is an error" "c-badcell: unordered 'yes' is neither - nor 1" "$TR"
[[ ! -e "$STUB_DOCKER_CALLED" ]] && ok "no container ran in the reference-cell checks" || bad "docker was called in the reference-cell checks"
cp "$TMP/docker-answering" "$STUBDIR/docker"
{
    echo "$HEADER"
    row c-stored fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - "$REF_CELL" - -
    row c-live fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-unordered fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - --pick_allele default,json "$IMAGE" offline 116.2 - - - - 1 - - -
    row c-layout fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - - default "$IMAGE" offline 116.2 - - - - - - - -
    row c-layout-near fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - "--nearest symbol" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-layout-refseq fixture "$CORPUS_INPUT" "116_GRCh38=layout:$DATA/cache_root" - --refseq default "$IMAGE" offline 116.2 - - - - - "$REF_CELL" - -
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
    row c-layout-merged fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - --merged default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TU"
expect_error "a layout root lacking the leaf the flags select fails preflight by its layout path" "c-layout-merged: cache root $DATA/cache_root has no homo_sapiens_merged/116_GRCh38/transcripts/" "$TU"
{
    echo "$HEADER"
    row c-layout fixture "$CORPUS_INPUT" 116_GRCh38=layout:cache_root - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv.gz,gene" default,tab "$IMAGE" offline 116.2 - - - - - - refusals/mixed -
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
[[ "$rc" -eq 0 ]] && grep -q "^DRY-RUN c-refuse refusal: $DATA/refusals/mixed (perl.txt, candidate.txt)$" "$TMP/df.out" && [[ "$(count_lines '^DRY-RUN c-refuse perl ' "$TMP/df.out")" -eq 2 && "$(count_lines '^DRY-RUN c-refuse stub ' "$TMP/df.out")" -eq 2 ]] && ! grep -q -- ' --report ' "$TMP/df.out" &&
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
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - refusals/incomplete -
} >"$TF"
expect_error "a refusal directory lacking candidate.txt is an error" "c-refuse: refusal directory lacks candidate.txt: $DATA/refusals/incomplete" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - refusals/blank -
} >"$TF"
expect_error "a blank refusal text is an error (it would match every stderr)" "c-refuse: refusal text $DATA/refusals/blank/candidate.txt is empty" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - "$DATA/refusals/absent" -
} >"$TF"
expect_error "a refusal cell naming no directory is an error" "c-refuse: refusal directory not found: $DATA/refusals/absent" "$TF"
[[ "$(grep -c 'c-refuse' "$TMP/e.out.err")" -eq 1 ]] && ok "a missing refusal directory is reported once, without a line per text file" || bad "missing refusal directory reported $(grep -c 'c-refuse' "$TMP/e.out.err") times"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - "$REF_CELL" refusals/mixed -
} >"$TF"
expect_error "refusal beside reference is an error" "c-refuse: refusal and reference are both set" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - IMPACT=1 - - - - refusals/mixed -
} >"$TF"
expect_error "refusal beside a floor is an error" "c-refuse: refusal with floors 'IMPACT=1'" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default,warnings "$IMAGE" offline 116.2 - - - - - - refusals/mixed -
} >"$TF"
expect_error "refusal beside warnings is an error" "c-refuse: refusal with warnings in formats" "$TF"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default,stats "$IMAGE" offline 116.2 - - - - - - refusals/mixed -
} >"$TF"
expect_error "refusal beside stats is an error" "c-refuse: refusal with stats in formats; a refused run writes no statistics report" "$TF"

echo
echo "UNINDEXED FASTA: a fasta cell marked unindexed: or unindexed-readonly: gives every engine run a copy of its own, without an index"
# A plain FASTA, a copy of it whose own name ends .tmp and a file that opens as a BGZF file
# does, none with an index beside it;
# the text both engines refuse a directory they cannot write with; a stand-in `id` that
# prints STUB_ID_U as the user's number and hands every other call to /usr/bin/id; and a
# stand-in `setpriv` that logs its arguments and runs the command after its one option.
FI="$DATA/unindexed"
FI_BIN="$TMP/fi-bin"
mkdir -p "$FI/refusal" "$FI_BIN/id" "$FI_BIN/setpriv"
cp "$RS/ref.fa" "$FI/plain.fa"
cp "$RS/ref.fa" "$FI/named.fa.tmp"
bgzf_head >"$FI/blocks.fa.gz"
FI_TEXT="Could not write to FASTA lock file"
printf '%s\n' "$FI_TEXT" >"$FI/refusal/perl.txt"
printf '%s\n' "$FI_TEXT" >"$FI/refusal/candidate.txt"
cat >"$FI_BIN/id/id" <<'EOF'
#!/usr/bin/env bash
if [[ "$*" == "-u" && -n "${STUB_ID_U:-}" ]]; then
    echo "$STUB_ID_U"
    exit 0
fi
exec /usr/bin/id "$@"
EOF
cat >"$FI_BIN/setpriv/setpriv" <<'EOF'
#!/usr/bin/env bash
[[ -z "${STUB_SETPRIV_LOG:-}" ]] || printf '%s\n' "$*" >>"$STUB_SETPRIV_LOG"
shift
exec "$@"
EOF
chmod +x "$FI_BIN/id/id" "$FI_BIN/setpriv/setpriv"
# fi_path_without <command>: PATH with every directory that carries the command replaced by a
# directory of links to its other entries, so the runner's `command -v` finds none.
fi_path_without() {
    local name="$1" dir out="" farm n=0
    local IFS=':'
    for dir in $PATH; do
        if [[ -x "$dir/$name" ]]; then
            farm="$TMP/fi-without-$name/$((n += 1))"
            mkdir -p "$farm"
            ln -s "$dir"/* "$farm/" 2>/dev/null
            rm -f "$farm/$name"
            dir="$farm"
        fi
        out="$out${out:+:}$dir"
    done
    printf '%s' "$out"
}
fi_row() { # fi_row <case> <fasta cell> <formats> [masks cell] [reference cell] [refusal cell]
    row "$1" fixture "$CORPUS_INPUT" "$CORPUS_CACHE" "$2" "$CORPUS_FLAGS" "$3" "$IMAGE" offline 116.2 "${4:--}" - - - - "${5:--}" "${6:--}" -
}
T_FI="$TMP/fi.tsv"
{
    echo "$HEADER"
    fi_row u-plain unindexed:unindexed/plain.fa default,tab,vcf,json
    fi_row u-one unindexed:unindexed/plain.fa default
    fi_row u-blocks unindexed:unindexed/blocks.fa.gz default
    fi_row u-pq unindexed:unindexed/plain.fa tab,parquet
    fi_row u-ro unindexed-readonly:unindexed/plain.fa default,tab - - unindexed/refusal
    fi_row u-tmp unindexed:unindexed/named.fa.tmp default
} >"$T_FI"
# Preflight. Every run of the read-only case has the stand-in `id` first on PATH, so that the
# checks read the same whatever user runs them.
STUB_ID_U=1000 PATH="$FI_BIN/id:$PATH" dry "$T_FI" "$TMP/fi-dry.out"
rc=$?
[[ "$rc" -eq 0 && ! -e "$FI/plain.fa.fai" && ! -e "$FI/blocks.fa.gz.fai" ]] && ok "a marked case passes preflight with no .fai beside its source" || bad "marked cases: the dry run exited $rc: $(head -n 2 "$TMP/fi-dry.out.err" | cut -c1-200 | tr '\n' ' ')"
{
    echo "$HEADER"
    fi_row u-plain unindexed/plain.fa default
} >"$T3"
expect_error "the same case without the marker still fails on its missing index" "u-plain: FASTA index not found: $FI/plain.fa.fai" "$T3"
{
    echo "$HEADER"
    fi_row u-ro unindexed-readonly:unindexed/plain.fa default
    fi_row u-ref unindexed:unindexed/plain.fa default - "$REF_CELL"
    fi_row u-mask unindexed:unindexed/plain.fa default hgvs_shift_sequence_edge_window
} >"$T3"
expect_error "the read-only marker on a row that names no refusal is an error" "u-ro: fasta 'unindexed-readonly:unindexed/plain.fa' on a row whose refusal cell is -; a directory no engine can write is a refusal row's" "$T3"
expect_error "a marked case that names a stored reference is an error" "u-ref: fasta 'unindexed:unindexed/plain.fa' with reference '$REF_CELL'; each engine indexes a copy of its own in a live run, which a stored output does not hold" "$T3"
expect_error "a marked case whose masks name a class the comparator decides from the FASTA is an error" "u-mask: masks names 'hgvs_shift_sequence_edge_window' beside fasta 'unindexed:unindexed/plain.fa'; the comparator reads the case's FASTA through its index" "$T3"
FI_NO_SETPRIV=$(fi_path_without setpriv)
STUB_ID_U=0 PATH="$FI_BIN/id:$FI_NO_SETPRIV" expect_error "the read-only marker as user 0 with no setpriv on PATH is an error" "u-ro: fasta 'unindexed-readonly:unindexed/plain.fa' as user 0 without setpriv on PATH" "$T_FI"
STUB_ID_U=0 PATH="$FI_BIN/id:$FI_BIN/setpriv:$PATH" dry "$T_FI" "$TMP/fi-dry-root.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "as user 0 with setpriv on PATH the read-only case passes preflight" || bad "user 0 with setpriv: the dry run exited $rc: $(head -n 2 "$TMP/fi-dry-root.out.err" | cut -c1-200 | tr '\n' ' ')"

# The dry run: one directory per engine and output format, under the case's work directory.
FI_WORK="$TMP/never/work/cases"
grep -qxF "DRY-RUN u-plain fasta: a copy without its index for every engine run (unindexed)" "$TMP/fi-dry.out" && grep -qxF "DRY-RUN u-ro fasta: a copy without its index for every engine run (unindexed-readonly)" "$TMP/fi-dry.out" &&
    ok "a dry run says of each marked case that every engine run gets a copy without its index, with the case's marker" || bad "the dry run's fasta lines: $(grep ' fasta: ' "$TMP/fi-dry.out" | tr '\n' ' ')"
pc=$(cmd_after 'DRY-RUN u-plain perl default' "$TMP/fi-dry.out")
cc=$(cmd_after 'DRY-RUN u-plain stub default' "$TMP/fi-dry.out")
[[ "$pc" == *" -v $FI_WORK/u-plain/fasta/perl/default:/work/fasta "* && "$pc" != *"ref.fa"* && "$pc" == *" --fasta /work/fasta/plain.fa" ]] &&
    ok "Ensembl VEP gets its staged directory as one writable mount, no index mount, and --fasta inside it" || bad "Ensembl VEP's staged FASTA: $pc"
[[ "$cc" == *" --fasta $FI_WORK/u-plain/fasta/candidates/stub/default/plain.fa" ]] && ok "the candidate gets --fasta with its own staged copy" || bad "the candidate's staged FASTA: $cc"
n=0
for f in default tab vcf json; do
    pc=$(cmd_after "DRY-RUN u-plain perl $f" "$TMP/fi-dry.out")
    cc=$(cmd_after "DRY-RUN u-plain stub $f" "$TMP/fi-dry.out")
    [[ "$pc" == *" -v $FI_WORK/u-plain/fasta/perl/$f:/work/fasta "* && "$cc" == *" --fasta $FI_WORK/u-plain/fasta/candidates/stub/$f/plain.fa"* ]] && n=$((n + 1))
done
[[ "$n" -eq 4 ]] && ok "every engine run of a marked case is planned on a directory of its own, one per engine and output format" || bad "staged directories: $n of 4 formats planned on their own"
pc=$(cmd_after 'DRY-RUN u-blocks perl default' "$TMP/fi-dry.out")
cc=$(cmd_after 'DRY-RUN u-blocks stub default' "$TMP/fi-dry.out")
[[ "$pc" == *" --fasta /work/fasta/blocks.fa.gz" && "$cc" == *" --fasta $FI_WORK/u-blocks/fasta/candidates/stub/default/blocks.fa.gz" ]] && ! grep -q "$TMP/never/work/fasta/" "$TMP/fi-dry.out" &&
    ok "a marked .gz FASTA is planned as it is, with no decompressed copy" || bad "a marked .gz FASTA: perl=$pc candidate=$cc"
pc=$(cmd_after 'DRY-RUN u-pq perl tab' "$TMP/fi-dry.out")
cc=$(cmd_after 'DRY-RUN u-pq stub parquet' "$TMP/fi-dry.out")
[[ "$(count_lines '^DRY-RUN u-pq perl ' "$TMP/fi-dry.out")" -eq 1 && "$pc" == *" -v $FI_WORK/u-pq/fasta/perl/tab:/work/fasta "* && "$cc" == *" --fasta $FI_WORK/u-pq/fasta/candidates/stub/parquet/plain.fa "* ]] &&
    ok "a parquet run gets a staged copy of its own, and Ensembl VEP runs the tab format alone" || bad "a parquet run's staged copy: perl=$pc candidate=$cc"
pc=$(cmd_after 'DRY-RUN u-ro perl default' "$TMP/fi-dry.out")
cc=$(cmd_after 'DRY-RUN u-ro stub default' "$TMP/fi-dry.out")
[[ "$pc" == *" -v $FI_WORK/u-ro/fasta/perl/default:/work/fasta:ro "* && "$cc" == "$STUBDIR/vep "* ]] &&
    ok "the read-only form mounts Ensembl VEP's directory read-only, and the candidate of another user than 0 runs as it is" || bad "the read-only form's staged directory: perl=$pc candidate=$cc"
cc=$(cmd_after 'DRY-RUN u-ro stub default' "$TMP/fi-dry-root.out")
[[ "$cc" == "setpriv --bounding-set=-dac_override $STUBDIR/vep "* && "$(cmd_after 'DRY-RUN u-plain stub default' "$TMP/fi-dry-root.out")" == "$STUBDIR/vep "* ]] &&
    ok "as user 0 the candidate's command of a read-only case is planned through setpriv --bounding-set=-dac_override, and no other case's is" || bad "user 0, the candidate's command: $cc"

# The real path. Each stub hands the FASTA it is given to stub-fasta, which logs the mode of
# the FASTA's directory and what it holds, then writes what STUB_PERL_FASTA or STUB_CAND_FASTA
# names.
FI_LOG="$TMP/fi-fasta.log"
fi_run() { # fi_run <out name> <case> <runner arguments...>: one case of the table, the logs emptied first
    local out="$1" only="$2"
    shift 2
    : >"$FI_LOG"
    : >"$STUB_DOCKER_LOG"
    : >"$STUB_CAND_LOG"
    STUB_FASTA_LOG="$FI_LOG" STUB_ID_U="${STUB_ID_U:-1000}" PATH="$FI_BIN/id:$PATH" run_oracle "$T_FI" --out "$TMP/$out" --candidates "$STUBDIR/vep" --candidate-pins p1 --only "$only" "$@" >"$TMP/$out.out" 2>"$TMP/$out.err"
}
# fi_started <engine> [mode]: the engine's runs in the log that started on a directory holding
# plain.fa alone, of that mode when one is given.
fi_started() {
    awk -F'\t' -v e="$1" -v m="${2:-}" '$1 == e && $3 == "plain.fa" && (m == "" || $2 == m) {n++} END {print n + 0}' "$FI_LOG"
}
fi_reasons() { # fi_reasons <out name> <case>
    reasons_of "$TMP/$1/p1/$2/report.json" 2>/dev/null
}
STUB_PERL_FASTA=fai STUB_CAND_FASTA=fai fi_run out-fi-green u-plain
rc=$?
[[ "$rc" -eq 0 && "$(sed -n 2p "$TMP/out-fi-green/status.tsv" 2>/dev/null | cut -f3)" == "green" ]] && ok "a marked case whose engines write equal indexes beside their copies is green" || bad "equal indexes: rc=$rc $(fi_reasons out-fi-green u-plain) $(tail -n 1 "$TMP/out-fi-green.err" | cut -c1-200)"
[[ "$(wc -l <"$FI_LOG" | tr -d ' ')" -eq 8 && "$(fi_started perl)" -eq 4 && "$(fi_started candidate)" -eq 4 ]] && ok "each engine's run of each format started on a directory that held the source alone" || bad "the staged directories at each run's start: $(tr '\t\n' ' |' <"$FI_LOG")"
[[ "$rc" -eq 0 && ! -e "$TMP/out-fi-green/work/cases/u-plain/fasta" && ! -e "$FI/plain.fa.fai" ]] && ok "after the run no staged directory is left under the case's work directory, and the source has no index beside it" || bad "after the run: rc=$rc, the case's work directory holds $(ls "$TMP/out-fi-green/work/cases/u-plain" 2>/dev/null | tr '\n' ' ')"
python3 - "$TMP/out-fi-green/p1/u-plain/report.json" "$FI/plain.fa" "$TMP/out1/abc1234/c-four/report.json" <<'PY' && ok "the report records the source FASTA and the marker, with no index reason; a case without a marker records neither" || bad "a marked case's report"
import json, sys
r = json.load(open(sys.argv[1]))
assert r["fasta"] == sys.argv[2], r["fasta"]
assert r["fasta_staged"] == "unindexed", r["fasta_staged"]
assert r["fasta_index_reasons"] == [] and r["reasons"] == [] and r["verdict"] == "green", r
plain = json.load(open(sys.argv[3]))
assert plain["fasta_staged"] is None and plain["fasta_index_reasons"] == [], plain
PY
STUB_PERL_FASTA=fai STUB_CAND_FASTA=fai fi_run out-fi-jobs u-plain --format-jobs 4
rc=$?
[[ "$rc" -eq 0 && "$(wc -l <"$FI_LOG" | tr -d ' ')" -eq 8 && "$(fi_started perl)" -eq 4 && "$(fi_started candidate)" -eq 4 ]] && ok "under --format-jobs each format's run starts on a directory of its own and the case is green" || bad "--format-jobs over a marked case: rc=$rc $(fi_reasons out-fi-jobs u-plain)"
STUB_PERL_FASTA=fai STUB_CAND_FASTA=fai fi_run out-fi-keep u-plain --keep-outputs
n=0
for d in perl/default perl/tab perl/vcf perl/json candidates/p1/default candidates/p1/tab candidates/p1/vcf candidates/p1/json; do
    d="$TMP/out-fi-keep/work/cases/u-plain/fasta/$d"
    [[ -w "$d" && -f "$d/plain.fa" && -f "$d/plain.fa.fai" ]] && n=$((n + 1))
done
[[ "$n" -eq 8 ]] && ok "--keep-outputs keeps every staged directory, writable, with the index its engine wrote" || bad "--keep-outputs: $n of 8 staged directories kept"
# A pin may be named perl: the candidate's staged copies stand under candidates/<pin>/, apart
# from Ensembl VEP's under perl/.
: >"$FI_LOG"
STUB_FASTA_LOG="$FI_LOG" STUB_PERL_FASTA=fai STUB_CAND_FASTA=fai run_oracle "$T_FI" --out "$TMP/out-fi-perlpin" --candidates "$STUBDIR/vep" --candidate-pins perl --only u-one >/dev/null 2>"$TMP/out-fi-perlpin.err"
rc=$?
[[ "$rc" -eq 0 && "$(fi_started perl)" -eq 1 && "$(fi_started candidate)" -eq 1 ]] && ok "a candidate pinned perl starts on a staged copy of its own, apart from Ensembl VEP's" || bad "a marked case under the pin perl: rc=$rc $(tr '\t\n' ' |' <"$FI_LOG") $(tail -n 1 "$TMP/out-fi-perlpin.err" | cut -c1-200)"
STUB_PERL_FASTA="fai gzi" STUB_CAND_FASTA="fai gzi" fi_run out-fi-blocks u-blocks
rc=$?
[[ "$rc" -eq 0 && "$(awk -F'\t' '$3 == "blocks.fa.gz"' "$FI_LOG" | wc -l | tr -d ' ')" -eq 2 && ! -e "$TMP/out-fi-blocks/work/fasta" ]] && ok "a marked BGZF file reaches both engines as it is, nothing is decompressed, and the case is green where both write its .fai and its .gzi" || bad "a marked BGZF file: rc=$rc $(fi_reasons out-fi-blocks u-blocks)"
# fi_red <label> <out name> <case> <the case's reasons>: the last run was red for exactly those reasons.
fi_red() {
    local rc=$? label="$1" reasons
    reasons=$(fi_reasons "$2" "$3")
    [[ "$rc" -eq 1 && "$reasons" == "$4" ]] && ok "$label" || bad "$label: rc=$rc reasons=$reasons"
}
STUB_PERL_FASTA="fai gzi" STUB_CAND_FASTA=fai fi_run out-fi-nogzi u-blocks
fi_red "a BGZF file whose .gzi an engine did not write is red, naming the engine and the file" out-fi-nogzi u-blocks "default: fasta index: the candidate wrote no blocks.fa.gz.gzi beside a BGZF file"
STUB_PERL_FASTA=fai STUB_CAND_FASTA="" fi_run out-fi-nofai u-one
fi_red "a .fai an engine did not write is red" out-fi-nofai u-one "default: fasta index: the candidate wrote no plain.fa.fai"
STUB_PERL_FASTA="fai lock" STUB_CAND_FASTA=fai fi_run out-fi-lock u-one
fi_red "a lock file an engine left is red" out-fi-lock u-one "default: fasta index: Ensembl VEP left plain.fa.vep.lock"
python3 - "$TMP/out-fi-lock/p1/u-one/report.json" <<'PY' && ok "the format's comparison is reported green beside the index reason, which the record lists under fasta_index_reasons" || bad "the comparison's verdict beside an index reason"
import json, sys
r = json.load(open(sys.argv[1]))
entry = r["formats"]["default"]
assert entry["verdict"] == "green" and entry["reasons"] == [] and "report" in entry, entry
assert r["verdict"] == "red" and r["fasta_index_reasons"] == r["reasons"] == ["default: fasta index: Ensembl VEP left plain.fa.vep.lock"], r["reasons"]
PY
STUB_PERL_FASTA=fai STUB_CAND_FASTA="fai tmp" fi_run out-fi-tmp u-one
fi_red "a .tmp file an engine left is red" out-fi-tmp u-one "default: fasta index: the candidate left plain.fa.fai.tmp"
STUB_PERL_FASTA=fai STUB_CAND_FASTA=fai fi_run out-fi-tmpname u-tmp
rc=$?
[[ "$rc" -eq 0 && "$(sed -n 2p "$TMP/out-fi-tmpname/status.tsv" 2>/dev/null | cut -f3)" == "green" ]] && ok "a source whose own name ends .tmp is green where both engines index it: the copy is no file an engine left" || bad "a source named .tmp: rc=$rc $(fi_reasons out-fi-tmpname u-tmp)"
STUB_PERL_FASTA=fai STUB_CAND_FASTA="fai other" fi_run out-fi-other u-one
fi_red "an index whose bytes differ between the engines is red" out-fi-other u-one "default: fasta index: plain.fa.fai differs between Ensembl VEP and the candidate"
STUB_CAND_MUTATE=change_impact STUB_PERL_FASTA="fai lock" STUB_CAND_FASTA=fai fi_run out-fi-both u-one
rc=$?
reasons=$(fi_reasons out-fi-both u-one)
[[ "$rc" -eq 1 && "$reasons" == *"default: field IMPACT: 1 differing entries | default: fasta index: Ensembl VEP left plain.fa.vep.lock" ]] &&
    ok "a difference the comparator sees is reported beside the index reason" || bad "a comparator difference beside an index reason: rc=$rc reasons=$reasons"
STUB_CAND_RC=4 STUB_PERL_FASTA=fai fi_run out-fi-candrc u-one
rc=$?
reasons=$(fi_reasons out-fi-candrc u-one)
[[ "$rc" -eq 1 && "$reasons" == "default: candidate exited 4: stub vep: failing on request"* && "$reasons" != *"fasta index"* ]] &&
    ok "a format one engine's run failed is red for the run alone: no index is asked of it" || bad "a failed run of a marked case: rc=$rc reasons=$reasons"
STUB_PERL_FASTA="fai lock" STUB_CAND_FASTA=fai fi_run out-fi-pq u-pq
reasons=$(fi_reasons out-fi-pq u-pq)
[[ "$reasons" == *"tab: fasta index: Ensembl VEP left plain.fa.vep.lock"* && "$reasons" == *"parquet: fasta index: Ensembl VEP left plain.fa.vep.lock"* && "$(fi_started perl)" -eq 1 && "$(fi_started candidate)" -eq 2 ]] &&
    ok "a parquet format's checks read the staged directory of Ensembl VEP's tab run" || bad "parquet index checks: $reasons"
STUB_PERL_FASTA=fai STUB_CAND_FASTA=fai fi_run out-fi-refdir u-one --reference-dir "$REPO/tests/golden/116"
rc=$?
[[ "$rc" -eq 0 && "$(count_lines '^run ' "$STUB_DOCKER_LOG")" -eq 1 ]] && grep -q 'live=1 ' "$TMP/out-fi-refdir.out" "$TMP/out-fi-refdir.err" &&
    ok "a marked case under --reference-dir still starts the container, once per format, and is green" || bad "a marked case under --reference-dir: rc=$rc docker runs=$(count_lines '^run ' "$STUB_DOCKER_LOG") $(tail -n 1 "$TMP/out-fi-refdir.err" | cut -c1-200)"
STUB_PERL_FASTA="fai lock" STUB_CAND_FASTA=fai fi_run out-fi-disc u-one --discover
rc=$?
facts=$(python3 -c 'import json, sys
r = json.load(open(sys.argv[1]))
d = r["discovery"]
print(r["status"], d["errors"], [reason for reason in d["residual"]["reasons"] if "fasta index" in reason], r["fasta_index_reasons"])' "$TMP/out-fi-disc/p1/u-one/report.json" 2>&1)
[[ "$rc" -eq 0 && "$facts" == "discovered [] ['default: fasta index: Ensembl VEP left plain.fa.vep.lock'] ['default: fasta index: Ensembl VEP left plain.fa.vep.lock']" ]] &&
    ok "in a discovery pass an index finding is a residual difference, not a format that did not run" || bad "an index finding in a discovery pass: rc=$rc $facts"

# The read-only form: both engines refuse, and neither writes beside its copy.
FI_SETPRIV_LOG="$TMP/fi-setpriv.log"
fi_refused() { # fi_refused <out name> <runner arguments...>: the read-only case with both engines refusing by their texts, a stand-in setpriv on PATH
    local out="$1"
    shift
    : >"$FI_SETPRIV_LOG"
    STUB_PERL_RC=2 STUB_PERL_STDERR="MSG: ERROR: $FI_TEXT /work/fasta/plain.fa.vep.lock" STUB_CAND_RC=1 STUB_CAND_STDERR="ERROR: [fasta] $FI_TEXT" \
        STUB_SETPRIV_LOG="$FI_SETPRIV_LOG" PATH="$FI_BIN/setpriv:$PATH" fi_run "$out" u-ro "$@"
}
fi_refused out-fi-ro
rc=$?
[[ "$rc" -eq 0 && "$(sed -n 2p "$TMP/out-fi-ro/status.tsv" 2>/dev/null | cut -f3)" == "green" ]] && ok "a read-only case is green where both engines refuse with their texts and write nothing" || bad "the read-only case: rc=$rc $(fi_reasons out-fi-ro u-ro) $(tail -n 1 "$TMP/out-fi-ro.err" | cut -c1-200)"
[[ "$(wc -l <"$FI_LOG" | tr -d ' ')" -eq 4 && "$(fi_started perl 555)" -eq 2 && "$(fi_started candidate 555)" -eq 2 && "$(count_lines ':/work/fasta:ro ' "$STUB_DOCKER_LOG")" -eq 2 ]] &&
    ok "each engine's run of a read-only case starts on a directory of mode 555 that holds the source alone, Ensembl VEP's mounted read-only" || bad "the read-only directories at each run's start: $(tr '\t\n' ' |' <"$FI_LOG")"
[[ "$rc" -eq 0 && ! -s "$FI_SETPRIV_LOG" && "$(count_lines '^-i ' "$STUB_CAND_LOG")" -eq 2 ]] && ok "the candidate of another user than 0 does not run through setpriv, which is on PATH" || bad "setpriv as another user than 0: rc=$rc $(head -n 1 "$FI_SETPRIV_LOG" | cut -c1-200)"
[[ "$rc" -eq 0 && ! -e "$TMP/out-fi-ro/work/cases/u-ro/fasta" ]] && ok "the read-only directories go with the case's outputs" || bad "read-only directories left: rc=$rc"
STUB_CAND_FASTA="force fai" fi_refused out-fi-ro-wrote
fi_red "a file an engine wrote beside the FASTA of a refused run is red in each format" out-fi-ro-wrote u-ro "default: fasta index: the candidate wrote plain.fa.fai into a directory it was given read-only | tab: fasta index: the candidate wrote plain.fa.fai into a directory it was given read-only"
fi_refused out-fi-ro-keep --keep-outputs
n=0
for d in perl/default perl/tab candidates/p1/default candidates/p1/tab; do
    d="$TMP/out-fi-ro-keep/work/cases/u-ro/fasta/$d"
    [[ -w "$d" && "$(ls -A "$d" 2>/dev/null)" == "plain.fa" ]] && n=$((n + 1))
done
[[ "$n" -eq 4 ]] && ok "--keep-outputs keeps a read-only case's staged directories with their write bit given back" || bad "--keep-outputs over the read-only case: $n of 4 staged directories kept writable"
# A staged directory an interrupted run left in mode 555 would stop the removal of the case's
# work directory at the next run's start.
chmod 0555 "$TMP/out-fi-ro-keep/work/cases/u-ro/fasta/perl/default" 2>/dev/null
fi_refused out-fi-ro-keep
rc=$?
[[ "$rc" -eq 0 && ! -e "$TMP/out-fi-ro-keep/work/cases/u-ro/fasta" ]] && ok "a run over a work tree that holds a staged directory of mode 555 removes it and is green" || bad "a run over a left-over read-only directory: rc=$rc $(tail -n 1 "$TMP/out-fi-ro-keep.err" | cut -c1-200)"
chmod -R u+w "$TMP/out-fi-ro-keep" 2>/dev/null
STUB_ID_U=0 fi_refused out-fi-root
rc=$?
[[ "$rc" -eq 0 && "$(wc -l <"$FI_SETPRIV_LOG" | tr -d ' ')" -eq 2 && "$(count_lines "^--bounding-set=-dac_override $STUBDIR/vep -i " "$FI_SETPRIV_LOG")" -eq 2 ]] &&
    ok "as user 0 the candidate's runs of a read-only case go through setpriv --bounding-set=-dac_override, and the case is green" || bad "setpriv as user 0: rc=$rc $(head -n 1 "$FI_SETPRIV_LOG" | cut -c1-200) $(tail -n 1 "$TMP/out-fi-root.err" | cut -c1-200)"
# A format job that dies while its engine runs has not given its directory the write bit back.
STUB_CAND_KILL_JOB=1 fi_refused out-fi-ro-killed --format-jobs 2
rc=$?
[[ "$rc" -eq 1 && "$(sed -n 2p "$TMP/out-fi-ro-killed/status.tsv" 2>/dev/null | cut -f3)" == "red" && "$(fi_reasons out-fi-ro-killed u-ro)" == "default: the format's job ended without a result"* && ! -e "$TMP/out-fi-ro-killed/work/cases/u-ro/fasta" ]] &&
    ok "a read-only case whose format jobs died is red for the jobs, and its read-only directories go with its outputs" || bad "format jobs killed in a read-only case: rc=$rc $(tail -n 1 "$TMP/out-fi-ro-killed.err" | cut -c1-200)"
chmod -R u+w "$TMP/out-fi-ro-killed" 2>/dev/null
# An interrupted run: the candidate sleeps in its read-only directory when the signal lands.
rm -f "$TMP/fi-sleepers"
STUB_ID_U=1000 STUB_PERL_RC=2 STUB_CAND_SLEEP=120 STUB_CAND_PIDFILE="$TMP/fi-sleepers" PATH="$FI_BIN/id:$STUBDIR:$PATH" bash "$ORACLE" --cases "$T_FI" --data-dir "$DATA" --out "$TMP/out-fi-int" --candidates "$STUBDIR/vep" --candidate-pins p1 --only u-ro >/dev/null 2>&1 &
runner=$!
for _ in $(seq 1 60); do
    [[ -s "$TMP/fi-sleepers" ]] && break
    sleep 1
done
kill -TERM "$runner" 2>/dev/null
wait "$runner" 2>/dev/null
d="$TMP/out-fi-int/work/cases/u-ro/fasta/candidates/p1/default"
[[ -s "$TMP/fi-sleepers" && -d "$d" && -w "$d" ]] && rm -rf "$TMP/out-fi-int" 2>/dev/null &&
    ok "an interrupted run of a read-only case leaves each staged directory with its write bit, so its work tree can be removed" || bad "an interrupted run of a read-only case: $(ls -ld "$d" 2>&1 | cut -c1-120)"
chmod -R u+w "$TMP/out-fi-int" 2>/dev/null

echo
echo "OPTIONS: --keep-outputs, --jobs and a .gz FASTA on the real path"
run_oracle "$T4" --out "$TMP/out-keep" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs --work-dir "$TMP/work-keep" >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && -f "$TMP/work-keep/cases/c-four/perl/default.txt" && -f "$TMP/work-keep/cases/c-four/candidates/p1/json.jsonl" ]] && ok "--keep-outputs keeps both engines' outputs under --work-dir" || bad "--keep-outputs: rc=$rc"
T8="$TMP/t8.tsv"
{
    echo "$HEADER"
    row c-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - - - -
    row c-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" vcf,json "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-near-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest symbol" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-near-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--nearest gene" default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-near-dash fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "-nearest symbol" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-near-upper fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--NEAREST gene" default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-cache fixture "$CORPUS_INPUT" "GRCh38=tests/golden/116/GRCh38-release116/json_cache" - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a cache cell that is not <version>_<assembly> is an error" "cache 'GRCh38' is not <version>_<assembly>" "$T3"
{
    echo "$HEADER"
    row c-img fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default - offline 116.2 - - - - - - - -
} >"$T3"
expect_error "a - perl_image is an error" "perl_image is -" "$T3"
{
    echo "$HEADER"
    row c-rel fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline - - - - - - - - -
} >"$T3"
expect_error "a - reference_release is an error" "reference_release is -" "$T3"
mkdir -p "$DATA/other"
cp "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.gz" "$DATA/other/reference.fa.gz"
cp "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.fai" "$DATA/other/reference.fa.fai"
{
    echo "$HEADER"
    row c-fa1 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - - -
    row c-fa2 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache other/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-fa2.out"
rc=$?
fa1=$(cmd_after 'DRY-RUN c-fa1 stub default' "$TMP/d-fa2.out")
fa2=$(cmd_after 'DRY-RUN c-fa2 stub default' "$TMP/d-fa2.out")
[[ "$rc" -eq 0 && "$fa1" == *"--fasta $TMP/never/work/fasta/1/reference.fa" && "$fa2" == *"--fasta $TMP/never/work/fasta/2/reference.fa" ]] &&
    ok "two .gz FASTAs that decompress to one file name each get a directory of their own" || bad "two FASTAs of one name: rc=$rc fa1=$fa1 fa2=$fa2"
{
    echo "$HEADER"
    row c-fa1 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - - -
    row c-fa1b fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
dry "$T3" "$TMP/d-fa3.out"
fa1=$(cmd_after 'DRY-RUN c-fa1 stub default' "$TMP/d-fa3.out")
fa1b=$(cmd_after 'DRY-RUN c-fa1b stub default' "$TMP/d-fa3.out")
[[ "$fa1" == *"--fasta $TMP/never/work/fasta/1/reference.fa" && "$fa1b" == *"--fasta $TMP/never/work/fasta/1/reference.fa" ]] &&
    ok "two cases on one .gz FASTA share its one decompressed copy" || bad "shared FASTA: fa1=$fa1 fa1b=$fa1b"
{
    echo "$HEADER"
    row c-refseq fixture "$CORPUS_INPUT" 116_GRCh38 - --refseq default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
expect_error "--refseq in the flags selects the homo_sapiens_refseq slice" "Ensembl VEP cache slice not found: $DATA/caches/perl/vep-cache/homo_sapiens_refseq/116_GRCh38" "$T3"
dry "$T1" "$TMP/d8.out" --only GRCh38-release116
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "names 'GRCh38-release116', which is not a case of tier fixture" "$TMP/d8.out.err" && ok "--only naming a case of the other tier is an error, not a silent drop" || bad "--only across tiers exited $rc"
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-glob fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--filter *" default "$IMAGE" offline 116.2 - - - - - - - -
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
    row fasta fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - - -
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
echo "LOOKUPS: a table read by key is read whole, whatever its size"
# Two cases on two images. The stand-in answers the second image with a digest longer than a
# pipe holds (65,536 bytes on Linux), so an awk that takes its input a block at a time, as
# gawk does, finds the first image's row while the table is still being written: a reader
# that left at that row would end the writer, and the run with it. The digest stays under
# the 131,072 bytes one argument may hold there, since the runner hands the table to its
# record as one argument; mawk reads that much ahead, and the check of the runner's text
# at the end of this block holds the form for every awk.
LK_IMAGE="ensemblorg/ensembl-vep@sha256:$(printf '%064d' 1)"
LK_ZEROS=120000
mkdir -p "$TMP/lk-bin"
cat >"$TMP/lk-bin/docker" <<STUB
#!/usr/bin/env bash
if [[ "\$1" == image && "\${*: -1}" == "$LK_IMAGE" ]]; then
    printf 'ensemblorg/ensembl-vep@sha256:%0${LK_ZEROS}d\n' 0
    exit 0
fi
exec "$STUBDIR/docker" "\$@"
STUB
chmod +x "$TMP/lk-bin/docker"
T_LK="$TMP/lk.tsv"
{
    echo "$HEADER"
    row c-first fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-second fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$LK_IMAGE" offline 116.2 - - - - - - - -
} >"$T_LK"
PATH="$TMP/lk-bin:$STUBDIR:$PATH" bash "$ORACLE" --cases "$T_LK" --data-dir "$DATA" --out "$TMP/out-lk" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/lk.out" 2>"$TMP/lk.err"
rc=$?
[[ "$rc" -eq 0 && "$(cut -f1,3 "$TMP/out-lk/status.tsv" 2>/dev/null | tr '\n\t' '  ')" == "case verdict c-first green c-second green " ]] &&
    ok "a run whose table of image digests is longer than a pipe holds ends with both cases green" || bad "long digest table: rc=$rc $(tail -n 2 "$TMP/lk.err" | cut -c1-200 | tr '\n' ' ')"
lk_digests=$(python3 -c 'import json, sys
for path in sys.argv[1:]:
    digest = json.load(open(path, encoding="utf-8"))["perl_image_digest"]
    print(digest if len(digest) < 100 else len(digest))' "$TMP/out-lk/p1/c-first/report.json" "$TMP/out-lk/p1/c-second/report.json" 2>/dev/null | tr '\n' ' ')
[[ "$lk_digests" == "ensemblorg/ensembl-vep@sha256:stubdigest $((30 + LK_ZEROS)) " ]] &&
    ok "each case's report carries the digest of its own image's row" || bad "digests by image: $lk_digests"
# The run facts of a manifest one of whose runs has a name of 700,000 bytes, more than a
# pipe holds and more than mawk reads ahead: the run the flags fit, the run a row names and
# the refusal of a case that fits no run each read the facts by key before their last line
# is written.
LK="$DATA/lookups"
mkdir -p "$LK/corpus"
cp "$CORPUS/variants.vcf" "$LK/corpus/variants.vcf"
python3 -c 'import json, sys
runs = [{"name": "picked", "flags": ["--symbol", "--pick"]}, {"name": "n" * 700000, "input": "second.vcf", "flags": ["--symbol"]}]
json.dump({"flags": ["--symbol"], "flag_sets": runs}, open(sys.argv[1], "w", encoding="utf-8"))' "$LK/corpus/manifest.json"
{
    echo "$HEADER"
    row l-own fixture lookups/corpus/variants.vcf "$CORPUS_CACHE" - "--symbol" default "$IMAGE" offline 116.2 - - - - - - - -
    row l-named fixture lookups/corpus/variants.vcf "$CORPUS_CACHE" - "--symbol --canonical" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TMP/lk-runs.tsv"
{
    row case run reason
    row l-named picked "the row names the second of the two runs over the input"
} >"$TMP/lk-named.tsv"
dry "$TMP/lk-runs.tsv" "$TMP/lk-runs.out" --case-runs "$TMP/lk-named.tsv"
rc=$?
lk_cmd=$(cmd_second 'DRY-RUN l-own stub default' "$TMP/lk-runs.out")
[[ "$rc" -eq 0 && "$lk_cmd" == *" --manifest $LK/corpus/manifest.json --run default "* ]] &&
    ok "a case whose flags fit a run is named by it where the manifest's run facts are longer than a pipe holds" || bad "long run facts, the run the flags fit: rc=$rc $lk_cmd"
lk_cmd=$(cmd_second 'DRY-RUN l-named stub default' "$TMP/lk-runs.out")
[[ "$rc" -eq 0 && "$lk_cmd" == *" --manifest $LK/corpus/manifest.json --run picked "* ]] &&
    ok "a row's run is found among run facts longer than a pipe holds" || bad "long run facts, the run a row names: rc=$rc $lk_cmd"
{
    echo "$HEADER"
    row l-unfit fixture lookups/corpus/variants.vcf "$CORPUS_CACHE" - "--symbol --canonical" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TMP/lk-unfit.tsv"
dry "$TMP/lk-unfit.tsv" "$TMP/lk-unfit.out"
rc=$?
[[ "$rc" -eq 1 ]] && grep -q "^ERROR: \[feature_oracle\] l-unfit: $LK/corpus/manifest.json lists 2 runs over variants.vcf (default, picked) and the flags are not exactly one of them" "$TMP/lk-unfit.out.err" &&
    ok "a case that fits no run of such a manifest is refused by name, with the runs over its input" || bad "long run facts, a case that fits no run: rc=$rc $(head -n 2 "$TMP/lk-unfit.out.err" | cut -c1-200 | tr '\n' ' ')"
# The runner's own text. A writer of several lines is ended at random as well, between two
# of its writes, so the form is held on every line, whatever the size of its table.
lk_early=$(awk '!/^[ \t]*#/ && /\|[ \t]*awk[ \t].*exit[ \t]*[;}]/ {printf "%s%d", sep, NR; sep = ", "}' "$ORACLE")
[[ -z "$lk_early" ]] && ok "no line of the runner pipes a writer into an awk that leaves at its first match" || bad "a writer piped into an awk that leaves at its first match: feature_oracle.sh line(s) $lk_early"

echo
echo "EMPTY AND BROKEN OUTPUTS"
T9="$TMP/t9.tsv"
{
    echo "$HEADER"
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-db37 fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache - --hgvs default "$IMAGE" database 116.2 - - - - - - - -
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
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-masked fixture "$TMP/masked/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 breakend_mate_local_read - - - - - - -
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
    row c-swap fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 splice_family_swap=1 - - - - - - -
    row c-swapcount fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 splice_family_swap=2 - - - - - - -
    row c-swapidle fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 splice_family_swap - - - - - - -
    row c-swapboth fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 splice_family_swap=1 - - - - - - -
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
    row c-needed fixture "$TMP/needed/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS" - - - - - - -
    row c-needed-count fixture "$TMP/needed/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS=1" - - - - - - -
    row c-needed-idle fixture "$TMP/needed/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS,$END_CLASS" - - - - - - -
    row c-plain fixture "$TMP/plain/variants.vcf" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 "$NEW_CLASS=1,$START_CLASS" - - - - - - -
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
# The same comparison through a comparator that, after the tool, lists every class of the
# report without a count as unused: the oracle reads the list, so the case is red.
cat >"$TMP/listing-comparator" <<'EOF'
#!/usr/bin/env bash
"$LISTED_COMPARATOR" "$@"
rc=$?
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
chmod +x "$TMP/listing-comparator"
LISTED_COMPARATOR="$COMPARATOR" VEP_FORMATS_COMPARATOR="$TMP/listing-comparator" STUB_REF_DIR="$TMP/needed" STUB_CAND_MUTATE=own_span_protein run_oracle "$T3" --out "$TMP/out-listed" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-needed >"$TMP/listed.out" 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-listed/p1/c-needed/report.json")
[[ "$rc" -eq 1 && "$reasons" == "mask $START_CLASS is requested and applied in no compared format" ]] &&
    ok "the same comparison is red where its report lists the class unused" || bad "class listed unused: rc=$rc reasons=$reasons"
mkdir -p "$TMP/refbad/c-badprov/expected"
cp "$REF/default.txt" "$TMP/refbad/c-badprov/expected/default.txt"
echo '{not json' >"$TMP/refbad/c-badprov/provenance.json"
{
    echo "$HEADER"
    row c-badprov fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - - -
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
    row c-miss fixture "$DATA/absent.vcf" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-pre" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
[[ "$rc" -ne 0 && ! -e "$TMP/out-pre/status.tsv" ]] && ok "a run that fails preflight removes the previous run's status.tsv" || bad "preflight failure left status.tsv (rc=$rc)"
mkdir -p "$TMP/out-dry"
echo "stale" >"$TMP/out-dry/status.tsv"
{
    echo "$HEADER"
    row c-ok fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$T3"
run_oracle "$T3" --out "$TMP/out-dry" --candidates "$STUBDIR/vep" --candidate-pins p1 --dry-run >/dev/null 2>&1
[[ -e "$TMP/out-dry/status.tsv" ]] && ok "a dry run leaves an existing status.tsv alone" || bad "a dry run removed status.tsv"

echo
echo "TOOL ROWS: a filter_vep cell runs the companion filter tools over an annotated input"
# The annotated inputs: the committed Ensembl VEP VCF, tab and default outputs, which the
# stubs also answer with, so the "filtered" outputs are identical and the real comparator
# decides.
TOOLDIR="$DATA/tooled"
mkdir -p "$TOOLDIR"
cp "$REF/vcf.vcf" "$TOOLDIR/annotated.vcf"
cp "$REF/tab.txt" "$TOOLDIR/annotated.tab"
cp "$REF/default.txt" "$TOOLDIR/annotated.txt"
TOOL_VCF="tooled/annotated.vcf"
MISSENSE='filter_vep:--filter "Consequence is missense_variant and IMPACT is MODERATE"'
TT="$TMP/tt.tsv"
{
    echo "$HEADER"
    row t-vcf fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-tab fixture tooled/annotated.tab "$CORPUS_CACHE" - - tab "$IMAGE" offline 116.2 - - - - - - - 'filter_vep:--filter "SYMBOL match ^TPTE"'
    row t-count fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - 'filter_vep:--count --filter "IMPACT is HIGH"'
    row t-db fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" database 116.2 - - - - - - - 'filter_vep:--ontology --filter "Consequence is coding_sequence_variant"'
    row t-default fixture tooled/annotated.txt "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - 'filter_vep:--filter "IMPACT is HIGH"'
    row c-four fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json,warnings "$IMAGE" offline 116.2 - - - - - - - -
} >"$TT"
# A tool row needs no cache: the table's cache cell is recorded alone, so a data tree without
# the slice or the candidate cache still plans the tool rows.
dry "$TT" "$TMP/dt-nocache.out" --perl-cache-dir "$DATA/x" --json-cache-dir "$DATA/x" --only t-vcf,t-tab,t-count,t-db,t-default
rc=$?
[[ "$rc" -eq 0 ]] && ok "tool rows plan without the Ensembl cache slice or the candidate cache" || {
    bad "tool rows without caches exited $rc"
    head -3 "$TMP/dt-nocache.out.err" | sed 's/^/         /'
}
# The cache cell's own =<dir> and =layout:<root> forms name no directory the tool row reads either.
TTE="$TMP/tte.tsv"
{
    echo "$HEADER"
    row t-nocorpus fixture "$TOOL_VCF" 116_GRCh38=tests/golden/116/no-such-corpus/json_cache - - vcf "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-noroot fixture "$TOOL_VCF" 116_GRCh38=layout:tests/golden/116/no-such-root - - vcf "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row c-nocorpus fixture "$CORPUS_INPUT" 116_GRCh38=tests/golden/116/no-such-corpus/json_cache - - default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TTE"
dry "$TTE" "$TMP/dt-nocorpus.out" --only t-nocorpus,t-noroot
rc=$?
[[ "$rc" -eq 0 ]] && grep -q "^DRY-RUN t-noroot cache species=homo_sapiens perl=.* candidate=.*/tests/golden/116/no-such-root/homo_sapiens/116_GRCh38$" "$TMP/dt-nocorpus.out" &&
    ok "a tool row's absent =<dir> or =layout:<root> cache is recorded and not required" || bad "tool rows with absent cache cells exited $rc: $(head -2 "$TMP/dt-nocorpus.out.err" | tr '\n' ' ')"
expect_error "an engine row's absent =<dir> cache is still refused" "c-nocorpus: candidate JSON cache has no transcripts/" "$TTE" --only c-nocorpus
dry "$TT" "$TMP/dt.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "a table with tool rows plans" || {
    bad "tool dry run exited $rc"
    head -5 "$TMP/dt.out.err" | sed 's/^/         /'
}
n=$(count_lines '^DRY-RUN t-vcf perl ' "$TMP/dt.out")
m=$(count_lines '^DRY-RUN t-vcf stub ' "$TMP/dt.out")
[[ "$n" -eq 1 && "$m" -eq 1 ]] && ok "a tool row plans one Ensembl VEP tool run and one candidate tool run, its one format" || bad "tool row planned $n Ensembl VEP and $m candidate runs"
tp=$(cmd_after 'DRY-RUN t-vcf perl vcf' "$TMP/dt.out")
[[ "$tp" == "docker run --rm --user $(id -u):$(id -g) -e PERL_HASH_SEED=0 -e PERL_PERTURB_KEYS=0 -v $TOOLDIR/annotated.vcf:/work/input/annotated.vcf:ro -v $TMP/never/work/cases/t-vcf/perl:/work/out $IMAGE filter_vep -i /work/input/annotated.vcf -o /work/out/vcf.vcf --force_overwrite --filter Consequence\\ is\\ missense_variant\\ and\\ IMPACT\\ is\\ MODERATE" ]] &&
    ok "Ensembl VEP's side runs filter_vep in the image over the mounted input with -i, -o, --force_overwrite and the row's words, the filter string one word" || bad "filter_vep command: $tp"
[[ "$tp" != *"/work/cache"* && "$tp" != *"--offline"* && "$tp" != *"--format"* && "$tp" != *"--fork"* && "$tp" != *"--species"* && "$tp" != *"_stats"* ]] && ok "the tool run mounts no cache and carries none of the engine flags" || bad "engine flags leaked into the tool run: $tp"
tc=$(cmd_after 'DRY-RUN t-vcf stub vcf' "$TMP/dt.out")
[[ "$tc" == "$STUBDIR/vep filter -i $TOOLDIR/annotated.vcf -o $TMP/never/work/cases/t-vcf/candidates/stub/vcf.vcf --force_overwrite --filter Consequence\\ is\\ missense_variant\\ and\\ IMPACT\\ is\\ MODERATE" ]] &&
    ok "the candidate's side runs vep filter over the same input with the same words and no cache" || bad "vep filter command: $tc"
tcc=$(cmd_second 'DRY-RUN t-vcf stub vcf' "$TMP/dt.out")
[[ "$tcc" == *"--format vcf --expected $TMP/never/work/cases/t-vcf/perl/vcf.vcf --actual $TMP/never/work/cases/t-vcf/candidates/stub/vcf.vcf --report $TMP/never/stub/t-vcf/vcf.json --input $TOOLDIR/annotated.vcf" ]] &&
    ok "the two filtered outputs are compared in the row's format, the annotated input keying the VCF records" || bad "tool comparison: $tcc"
tt=$(cmd_after 'DRY-RUN t-tab perl tab' "$TMP/dt.out")
[[ "$tt" == *"filter_vep -i /work/input/annotated.tab -o /work/out/tab.txt --force_overwrite --filter SYMBOL\\ match\\ \\^TPTE" ]] && ok "a tab row filters the tab input into tab.txt" || bad "tab tool command: $tt"
td=$(cmd_after 'DRY-RUN t-db perl vcf' "$TMP/dt.out")
[[ "$td" == *"filter_vep -i /work/input/annotated.vcf -o /work/out/vcf.vcf --force_overwrite --host ensembldb.ensembl.org --port 3306 --user anonymous --version 116 --ontology --filter "* ]] &&
    ok "in database mode filter_vep gets --host, --port, --user and --version <cache version> before the row's words" || bad "database tool command: $td"
tdc=$(cmd_after 'DRY-RUN t-db stub vcf' "$TMP/dt.out")
[[ "$tdc" == *" filter -i "* && "$tdc" != *"--host"* && "$tdc" != *"--version"* && "$tdc" == *"--ontology --filter "* ]] && ok "the candidate of a database tool row gets the row's words alone" || bad "candidate under database tool row: $tdc"
tx=$(cmd_second 'DRY-RUN t-count stub vcf' "$TMP/dt.out")
[[ "$tx" == "compared as text: $TMP/never/work/cases/t-count/perl/vcf.vcf $TMP/never/work/cases/t-count/candidates/stub/vcf.vcf" ]] && ok "a --count row plans a text comparison of the two counts, not the comparator" || bad "count comparison line: $tx"
cf=$(cmd_after 'DRY-RUN c-four perl tab' "$TMP/dt.out")
[[ "$cf" == *" vep -i /work/input/variants.vcf -o /work/out/tab.txt --offline --cache --dir_cache /work/cache --cache_version 116 "* ]] && ok "a - tool cell leaves an engine row as it was" || bad "engine row changed: $cf"

{
    echo "$HEADER"
    row t-tool fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - haplo:--x
    row t-notool fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - ':--filter "IMPACT is HIGH"'
    row t-colon fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - filter_vep
    row t-empty fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - 'filter_vep: '
    row t-quote fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - 'filter_vep:--filter "IMPACT is HIGH'
    row t-out fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - 'filter_vep:-o x --filter "IMPACT is HIGH"'
    row t-list fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - filter_vep:--list
    row t-version fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" database 116.2 - - - - - - - 'filter_vep:--version 116 --ontology --filter "Consequence is coding_sequence_variant"'
    row t-onto fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - - - 'filter_vep:--ontology --filter "Consequence is coding_sequence_variant"'
    row t-two fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf,tab "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-json fixture "$TOOL_VCF" "$CORPUS_CACHE" - - json "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-warn fixture "$TOOL_VCF" "$CORPUS_CACHE" - - default,warnings "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-stats fixture "$TOOL_VCF" "$CORPUS_CACHE" - - default,stats "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-flags fixture "$TOOL_VCF" "$CORPUS_CACHE" - --hgvs vcf "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-fasta fixture "$TOOL_VCF" "$CORPUS_CACHE" tests/golden/116/GRCh37-hgvs/reference.fa.gz - vcf "$IMAGE" offline 116.2 - - - - - - - "$MISSENSE"
    row t-floor fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - IMPACT=1 - - - - - 'filter_vep:--count --filter "IMPACT is HIGH"'
} >"$TTE"
expect_error "a tool cell naming another tool is refused by name" "t-tool: tool 'haplo:--x' is not - or <tool>:<args> with a tool of: filter_vep" "$TTE" --only t-tool
expect_error "a tool cell with an empty tool name is refused, not run as an engine row" "t-notool: tool ':--filter \"IMPACT is HIGH\"' is not - or <tool>:<args> with a tool of: filter_vep" "$TTE" --only t-notool
expect_error "a tool cell without the colon is refused" "t-colon: tool 'filter_vep' is not - or <tool>:<args>" "$TTE" --only t-colon
expect_error "a tool cell with no argument is refused" "t-empty: tool 'filter_vep: ' carries no argument" "$TTE" --only t-empty
expect_error "a tool cell whose quotes do not close is refused as not shell words" "t-quote: tool arguments '--filter \"IMPACT is HIGH' are not shell words: No closing quotation" "$TTE" --only t-quote
expect_error "a tool argument the runner sets (-o) is refused" "t-out: tool arguments name -o, which the runner sets" "$TTE" --only t-out
expect_error "--list is refused: its output never reaches the output file" "t-list: tool arguments name --list, which the runner sets" "$TTE" --only t-list
expect_error "the database connection (--version) is the runner's in database mode and refused in the cell" "t-version: tool arguments name --version, which the runner sets" "$TTE" --only t-version
expect_error "--ontology under perl_mode offline is refused: filter_vep answers it from the database" "t-onto: tool arguments carry --ontology under perl_mode offline" "$TTE" --only t-onto
expect_error "two formats on a tool row are refused" "t-two: formats 'vcf,tab' on a tool row; the tool reads one annotated file and writes one, in one of: default tab vcf" "$TTE" --only t-two
expect_error "json on a tool row is refused" "t-json: formats 'json' on a tool row" "$TTE" --only t-json
expect_error "warnings on a tool row is refused" "t-warn: formats 'default,warnings' on a tool row" "$TTE" --only t-warn
expect_error "stats on a tool row is refused" "t-stats: formats 'default,stats' on a tool row" "$TTE" --only t-stats
expect_error "flags on a tool row are refused: no engine runs" "t-flags: flags '--hgvs' on a tool row; the row runs no engine, so the flags cell is -" "$TTE" --only t-flags
expect_error "a FASTA on a tool row is refused" "t-fasta: fasta 'tests/golden/116/GRCh37-hgvs/reference.fa.gz' on a tool row; the tool reads no FASTA" "$TTE" --only t-fasta
expect_error "a comparator cell on a --count row is refused" "t-floor: a --count tool row with floors 'IMPACT=1'; the two counts are compared as text, so the cell is -" "$TTE" --only t-floor
{
    row case run reason
    row t-count default "a count row compares against no manifest"
} >"$TMP/tt-runs.tsv"
expect_error "a named run on a --count tool row is refused" "t-count: $TMP/tt-runs.tsv names the run 'default' for the case, which compares against no manifest (no manifest.json beside its input, a refusal row or a --count tool row)" "$TT" --case-runs "$TMP/tt-runs.tsv" --only t-count
{
    row case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs species unordered reference refusal
    row c-x fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - -
} >"$TTE"
expect_error "a header without the tool column is an error" "header lacks the column 'tool'" "$TTE"

run_oracle "$TT" --out "$TMP/out-tool" --candidates "$STUBDIR/vep" --candidate-pins p1 --only t-vcf,t-tab,t-count,t-db,t-default >"$TMP/tool.out" 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(grep -c 'green' "$TMP/out-tool/status.tsv")" -eq 5 ]] && ok "on the real path five tool rows (vcf, tab, default, a count, a database row) with identical stub outputs are green" || {
    bad "tool real path: rc=$rc"
    tail -5 "$TMP/tool.out" | sed 's/^/         /'
}
python3 - "$TMP/out-tool/p1" "$MISSENSE" <<'PY' && ok "report.json records the tool cell, the text comparison's two counts and the database tool row's server" || bad "tool report.json: $(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["formats"])' "$TMP/out-tool/p1/t-count/report.json" 2>&1 | head -3)"
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
vcf = json.loads((out / "t-vcf" / "report.json").read_text())
assert vcf["tool"] == sys.argv[2] and vcf["perl_server"] is None, (vcf["tool"], vcf["perl_server"])
assert vcf["formats"]["vcf"]["verdict"] == "green" and "report" in vcf["formats"]["vcf"]
count = json.loads((out / "t-count" / "report.json").read_text())
entry = count["formats"]["vcf"]
assert entry == {"verdict": "green", "compared": "text", "text": {"expected": "7", "actual": "7"}, "reasons": []}, entry
db = json.loads((out / "t-db" / "report.json").read_text())
assert db["perl_server"] == "ensembldb.ensembl.org:3306" and db["perl_mode"] == "database", db["perl_server"]
PY
python3 -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1]))["tool"] is None else 1)' "$OUT/abc1234/c-four/report.json" && ok "an engine row's report records no tool" || bad "engine row's tool: $(grep '"tool"' "$OUT/abc1234/c-four/report.json")"
STUB_CAND_COUNT=6 run_oracle "$TT" --out "$TMP/out-tool6" --candidates "$STUBDIR/vep" --candidate-pins p1 --only t-count >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-tool6/p1/t-count/report.json" 2>/dev/null)
[[ "$rc" -eq 1 && "$reasons" == "vcf: count differs: Ensembl VEP '7', candidate '6'" ]] && ok "a count the candidate gets wrong is red with both counts in the reason" || bad "count mismatch: rc=$rc reasons=$reasons"
STUB_CAND_COUNT="" run_oracle "$TT" --out "$TMP/out-tool0" --candidates "$STUBDIR/vep" --candidate-pins p1 --only t-count >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-tool0/p1/t-count/report.json" 2>/dev/null)
[[ "$rc" -eq 1 && "$reasons" == "vcf: candidate output is empty where the tool writes its count" ]] && ok "an empty count output is red, never an agreement on nothing" || bad "empty count: rc=$rc reasons=$reasons"

echo
echo "TOOL ROWS: the reference cell may name the stored filtered output file"
# The stored file stands in for the Ensembl VEP tool run, its provenance the corpus's beside the
# input: here the tooled directory gets a provenance.json naming an image.
cp "$REF/vcf.vcf" "$TOOLDIR/ontology.vcf"
echo '{"vep_image": "ensemblorg/ensembl-vep@sha256:toolprov"}' >"$TOOLDIR/provenance.json"
TTS="$TMP/tts.tsv"
{
    echo "$HEADER"
    row t-stored fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" database 116.2 - - - - - tooled/ontology.vcf - 'filter_vep:--ontology --filter "Consequence is coding_sequence_variant"'
    row t-nofile fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" database 116.2 - - - - - "$TOOLDIR/absent.vcf" - 'filter_vep:--ontology --filter "Consequence is coding_sequence_variant"'
    row c-fileref fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - tooled/ontology.vcf - -
} >"$TTS"
dry "$TTS" "$TMP/dts.out" --only t-stored
rc=$?
[[ "$rc" -eq 0 ]] && grep -q "^DRY-RUN t-stored perl stored: $TOOLDIR/ontology.vcf$" "$TMP/dts.out" && ! grep -q '^docker run' "$TMP/dts.out" &&
    ok "a tool row whose reference names a file plans no container and reads the file" || bad "stored file dry run: rc=$rc $(grep 'perl stored' "$TMP/dts.out")"
sc=$(cmd_second 'DRY-RUN t-stored stub vcf' "$TMP/dts.out")
[[ "$sc" == *"--format vcf --expected $TOOLDIR/ontology.vcf --actual "* ]] && ok "the comparison's expected output is the stored file" || bad "stored file comparison: $sc"
expect_error "a tool row whose reference is neither a directory nor a file is refused with the file hint" "t-nofile: stored reference directory not found: $TOOLDIR/absent.vcf; a tool row may name the stored filtered output file instead" "$TTS" --only t-nofile
expect_error "an engine row's reference must still be a directory" "c-fileref: stored reference directory not found: $TOOLDIR/ontology.vcf" "$TTS" --only c-fileref
grep -q "c-fileref: stored reference directory not found: $TOOLDIR/ontology.vcf$" "$TMP/e.out.err" && ok "the file hint is given to tool rows alone" || bad "file hint on an engine row: $(grep c-fileref "$TMP/e.out.err")"
cat >"$STUBDIR/docker" <<'EOF'
#!/usr/bin/env bash
: >"${STUB_DOCKER_CALLED:?}"
exit 99
EOF
rm -f "$STUB_DOCKER_CALLED"
run_oracle "$TTS" --out "$TMP/out-stored" --candidates "$STUBDIR/vep" --candidate-pins p1 --only t-stored >"$TMP/stored.out" 2>&1
rc=$?
cp "$TMP/docker-answering" "$STUBDIR/docker"
[[ "$rc" -eq 0 && ! -e "$STUB_DOCKER_CALLED" ]] && python3 - "$TMP/out-stored/p1/t-stored/report.json" "$TOOLDIR" <<'PY' && ok "on the real path a stored filtered file is compared with no container, its provenance the corpus's beside the input" || bad "stored file real path: rc=$rc $(tail -2 "$TMP/stored.out" | tr '\n' ' ')"
import json, sys
r = json.load(open(sys.argv[1]))
tooled = sys.argv[2]
assert r["verdict"] == "green", r["reasons"]
assert r["reference_dir"] == f"{tooled}/ontology.vcf" and r["reference_provenance"] == f"{tooled}/provenance.json", (r["reference_dir"], r["reference_provenance"])
assert r["perl_image_digest"] == "ensemblorg/ensembl-vep@sha256:toolprov" and r["perl_server"] == "ensembldb.ensembl.org:3306", (r["perl_image_digest"], r["perl_server"])
PY
printf '7\n' | gzip -n >"$TOOLDIR/count.txt.gz"
{
    echo "$HEADER"
    row t-gzcount fixture "$TOOL_VCF" "$CORPUS_CACHE" - - vcf "$IMAGE" offline 116.2 - - - - - tooled/count.txt.gz - 'filter_vep:--count --filter "IMPACT is HIGH"'
} >"$TTS"
run_oracle "$TTS" --out "$TMP/out-gzcount" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/gzcount.out" 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["formats"]["vcf"]["text"]["expected"])' "$TMP/out-gzcount/p1/t-gzcount/report.json" 2>/dev/null)" == "7" ]] &&
    ok "a gzipped stored count file is read through gzip and compared as text" || bad "gz stored count: rc=$rc $(tail -2 "$TMP/gzcount.out" | tr '\n' ' ')"

echo
echo "SCORE FILES: a score-reading class is handed the SpliceAI plugin's score files and the cache"
# spliceai_multi_entry_cell reads the score files the case's flags hand the SpliceAI plugin
# (snv=, indel=), each the host path as its own --spliceai-table in the cell's order, and a
# row's gene from the candidate's cache; --plugin is read in every spelling vep accepts.
TSP="$TMP/tsp.tsv"
{
    echo "$HEADER"
    row c-spliceai fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv --plugin SpliceAI,snv={corpus}/scores.snv.vcf.gz,indel={corpus}/scores.indel.vcf.gz,cutoff=0.5" default "$IMAGE" offline 116.2 spliceai_multi_entry_cell=2 - - - - - - -
    row c-spliceai-joined fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin=SpliceAI,indel={corpus}/scores.indel.vcf.gz,split_output=1,snv={corpus}/scores.snv.vcf.gz" tab "$IMAGE" offline 116.2 spliceai_multi_entry_cell - - - - - - -
    row c-spliceai-one fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin SpliceAI,snv={corpus}/scores.snv.vcf.gz" default "$IMAGE" offline 116.2 spliceai_multi_entry_cell - - - - - - -
    row c-spliceai-spelled fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "-PLUGIN SpliceAI,indel={corpus}/scores.indel.vcf.gz" default "$IMAGE" offline 116.2 spliceai_multi_entry_cell - - - - - - -
    row c-spliceai-unmasked fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin SpliceAI,snv={corpus}/scores.snv.vcf.gz,indel={corpus}/scores.indel.vcf.gz" default "$IMAGE" offline 116.2 order_only - - - - - - -
} >"$TSP"
dry "$TSP" "$TMP/s.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "the table of the score-reading class dry-runs with exit 0" || bad "score-reading class rc=$rc: $(head -c 300 "$TMP/s.out.err")"
ssa=$(cmd_second 'DRY-RUN c-spliceai stub default' "$TMP/s.out")
ssj=$(cmd_second 'DRY-RUN c-spliceai-joined stub tab' "$TMP/s.out")
[[ "$ssa" == *"--masks spliceai_multi_entry_cell --json-cache $CORPUS/json_cache --spliceai-table $DATA/nomanifest/scores.snv.vcf.gz --spliceai-table $DATA/nomanifest/scores.indel.vcf.gz" &&
    "$ssj" == *"--masks spliceai_multi_entry_cell --json-cache $CORPUS/json_cache --spliceai-table $DATA/nomanifest/scores.indel.vcf.gz --spliceai-table $DATA/nomanifest/scores.snv.vcf.gz" ]] &&
    ok "a comparison naming spliceai_multi_entry_cell carries the cache and each of two score files as its own --spliceai-table, in the cell's order" || bad "SpliceAI tables: default=$ssa tab=$ssj"
sso=$(cmd_second 'DRY-RUN c-spliceai-one stub default' "$TMP/s.out")
[[ "$sso" == *"--masks spliceai_multi_entry_cell --json-cache $CORPUS/json_cache --spliceai-table $DATA/nomanifest/scores.snv.vcf.gz" ]] &&
    ok "a comparison naming spliceai_multi_entry_cell beside one score file carries --spliceai-table once" || bad "SpliceAI table, one file: $sso"
sss=$(cmd_second 'DRY-RUN c-spliceai-spelled stub default' "$TMP/s.out")
[[ "$sss" == *"--masks spliceai_multi_entry_cell --json-cache $CORPUS/json_cache --spliceai-table $DATA/nomanifest/scores.indel.vcf.gz" ]] &&
    ok "a score file is read behind --plugin in another spelling vep accepts" || bad "SpliceAI table, -PLUGIN: $sss"
ssu=$(cmd_second 'DRY-RUN c-spliceai-unmasked stub default' "$TMP/s.out")
[[ "$ssu" == *"--masks order_only" ]] &&
    ok "a comparison naming no score-reading class carries no --json-cache and no --spliceai-table beside the same flags" || bad "score options leaked: $ssu"
# The class sets nothing aside without a score file, so a cell naming it beside flags that
# hand the SpliceAI plugin none is refused before anything runs.
TSF="$TMP/tsf.tsv"
{
    echo "$HEADER"
    row c-spliceai-fileless fixture nomanifest/variants.vcf "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv --plugin SpliceAI,cutoff=0.5" default "$IMAGE" offline 116.2 spliceai_multi_entry_cell - - - - - - -
} >"$TSF"
expect_error "spliceai_multi_entry_cell beside a SpliceAI plugin without snv= or indel= is an error" "c-spliceai-fileless: masks names 'spliceai_multi_entry_cell', which reads the score files the flags hand the SpliceAI plugin (snv=, indel=), and the flags name none" "$TSF"
{
    echo "$HEADER"
    row c-spliceai-flagless fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only,spliceai_multi_entry_cell=3 - - - - - - -
} >"$TSF"
expect_error "spliceai_multi_entry_cell beside no flags is an error" "c-spliceai-flagless: masks names 'spliceai_multi_entry_cell', which reads the score files the flags hand the SpliceAI plugin (snv=, indel=), and the flags name none" "$TSF"
# A reading class outside the built-in list would be refused in a case without a manifest
# before its file could be handed over. The five lists are read as the runner sets them, a
# list a later line adds to with it.
builtin_listed=$(awk -F'"' '/^BUILTIN_MASKS="/ {print $2; exit}' "$ORACLE")
reading_listed=$(
    eval "$(grep -E '^(INPUT|CACHE|REFERENCE|CATALOG|SCORE)_READING_MASKS\+?="' "$ORACLE")"
    # The awk's last pattern names the lines the eval reads; a line that sets a list in
    # any other form is reported as a class outside the built-in list.
    awk '/(^|[^A-Za-z0-9_])(INPUT|CACHE|REFERENCE|CATALOG|SCORE)_READING_MASKS[+:]?=/ && !/^[ \t]*#/ && !/^(INPUT|CACHE|REFERENCE|CATALOG|SCORE)_READING_MASKS[+]?="/ {print "a-list-set-in-a-form-not-read:line-" NR}' "$ORACLE"
    for reading_list in INPUT CACHE REFERENCE CATALOG SCORE; do
        reading_list="${reading_list}_READING_MASKS"
        [[ -n "${!reading_list:-}" ]] || echo "a-list-the-runner-lacks:$reading_list"
        printf '%s\n' ${!reading_list:-}
    done
)
reading_outside=""
for reading_class in $reading_listed; do
    [[ " $builtin_listed " == *" $reading_class "* ]] || reading_outside="$reading_outside $reading_class"
done
[[ -n "$reading_listed" && -z "$reading_outside" && " $reading_listed " == *spliceai_multi_entry_cell*spliceai_multi_entry_cell* ]] &&
    ok "every class of the input-, cache-, reference-, catalog- and score-reading lists is a built-in class, the score-reading class in the cache-reading list as well" || bad "reading classes outside BUILTIN_MASKS:$reading_outside"

echo
echo "FASTA INDEX SPELLING: a compressed FASTA indexed as <name>.fa.gz.fai alone"
TGZ="$TMP/tgz.tsv"
mkdir -p "$DATA/gzfasta"
cp "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.gz" "$REPO/tests/golden/116/GRCh37-hgvs/reference.fa.gz.fai" "$DATA/gzfasta/"
{
    echo "$HEADER"
    row c-gzfai fixture "$REPO/tests/golden/116/GRCh37-hgvs/variants.vcf" 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache gzfasta/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TGZ"
run_oracle "$TGZ" --out "$TMP/out-gzfai" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$TMP/gzfai.out" 2>&1
rc=$?
[[ "$rc" -eq 0 && -f "$TMP/out-gzfai/work/fasta/1/reference.fa.fai" ]] && cmp -s "$DATA/gzfasta/reference.fa.gz.fai" "$TMP/out-gzfai/work/fasta/1/reference.fa.fai" &&
    ok "preflight accepts <name>.fa.gz.fai and the decompressed copy gets it as its .fai" || bad "gz.fai spelling: rc=$rc $(tail -2 "$TMP/gzfai.out" | tr '\n' ' ')"
rm -f "$DATA/gzfasta/reference.fa.gz.fai"
expect_error "a compressed FASTA with neither index spelling is refused naming the plain one" "c-gzfai: FASTA index not found: $DATA/gzfasta/reference.fa.fai" "$TGZ"


echo
echo "SIZING: each step is timed, and the case's report and run.json carry what sizing a machine needs"
SZ="$TMP/sizing"
mkdir -p "$SZ"
# sizing_check.py <check> <out dir> [value...]: exit 0 when the named statement holds of the run
# under <out dir>, whose one case is c-four under candidate p1.
cat >"$SZ/sizing_check.py" <<'PY'
import hashlib, json, os, sys
check, out = sys.argv[1], sys.argv[2]
case_dir = os.path.join(out, "p1", "c-four")
report = json.load(open(os.path.join(case_dir, "report.json")))
sizing = report["sizing"]
steps = sizing["steps"]
formats = ("default", "tab", "vcf", "json")
def number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)
if check == "program":
    comparator = sys.argv[3]
    digest = hashlib.sha256(open(comparator, "rb").read()).hexdigest()[:12]
    held = sizing["comparison_program"] == "vep-formats-compare" and sizing["comparison_program_version"] == "sha256:" + digest
elif check == "rows":
    held = set(sizing["rows_compared"]) == set(formats)
    for fmt in formats:
        counts = json.load(open(os.path.join(case_dir, fmt + ".json")))["counts"]
        held = held and sizing["rows_compared"][fmt] == max(counts["expected_entries"], counts["actual_entries"]) > 0
elif check == "steps":
    wanted = {f"candidate.{fmt}" for fmt in formats} | {f"comparison.{fmt}" for fmt in formats + ("warnings",)}
    held = wanted <= set(steps)
    # These steps stay far under 4 GiB: a larger peak is a figure in bytes, as Darwin's ru_maxrss is.
    for name in wanted:
        entry = steps[name]
        held = held and number(entry["wall_seconds"]) and isinstance(entry["peak_rss_kib"], int) and 0 < entry["peak_rss_kib"] < 1 << 22 and entry["exit"] == 0
elif check == "engine":
    live = [name for name in steps if name.startswith("ensembl-vep.")]
    held = bool(live) and all(steps[name]["peak_rss_kib"] == "not measured" and number(steps[name]["wall_seconds"]) for name in live)
elif check == "totals":
    measured = [entry["peak_rss_kib"] for name, entry in steps.items() if not name.startswith("ensembl-vep.")]
    held = (
        sizing["peak_memory_kib"] == max(measured)
        and isinstance(sizing["work_disk_peak_bytes"], int) and sizing["work_disk_peak_bytes"] > 0
        and isinstance(sizing["work_disk_growth_bytes"], int) and 0 <= sizing["work_disk_growth_bytes"] <= sizing["work_disk_peak_bytes"]
        and number(sizing["wall_minutes"]) and sizing["wall_minutes"] >= 0
    )
elif check == "run":
    run = json.load(open(os.path.join(out, "run.json")))
    held = (
        run["comparison_program"] == sizing["comparison_program"]
        and run["comparison_program_version"] == sizing["comparison_program_version"]
        and run["case_sizes"] == [{"case": "c-four", "candidate_pin": "p1", "sizing": sizing}]
    )
elif check == "unmeasured":
    held = (
        sizing["rows_compared"] == "not measured"
        and not [name for name in steps if name.startswith("comparison.")]
        and all(steps[f"candidate.{fmt}"]["exit"] == 4 for fmt in formats)
        and isinstance(sizing["peak_memory_kib"], int)
    )
elif check == "no-records":
    held = sizing["steps"] == {} and all(
        sizing[key] == "not measured" for key in ("peak_memory_kib", "work_disk_peak_bytes", "work_disk_growth_bytes")
    ) and number(sizing["wall_minutes"]) and report["verdict"] == "green"
elif check == "signalled":
    held = all(steps[f"candidate.{fmt}"]["exit"] == 143 for fmt in formats)
else:
    held = False
sys.exit(0 if held else 1)
PY
run_oracle "$T4" --out "$SZ/out" --candidates "$STUBDIR/vep" --candidate-pins p1 >"$SZ/run.out" 2>"$SZ/run.err"
rc=$?
[[ "$rc" -eq 0 ]] && python3 "$SZ/sizing_check.py" program "$SZ/out" "$COMPARATOR" &&
    ok "a case's report names the comparison program, vep-formats-compare, and its version, the sha256 of the file that ran as it" || bad "sizing program: rc=$rc $(tail -2 "$SZ/run.err" | tr '\n' ' ')"
python3 "$SZ/sizing_check.py" rows "$SZ/out" && ok "rows_compared holds, per compared format, the larger entry count of the comparator's report" || bad "sizing rows_compared differs from the format reports' counts"
python3 "$SZ/sizing_check.py" steps "$SZ/out" && ok "every candidate run and every comparison, the warnings one included, leaves its wall seconds, its peak resident set in KiB and its exit status" || bad "sizing steps lack a candidate run's or a comparison's figures, or carry a peak that is not in KiB"
python3 "$SZ/sizing_check.py" engine "$SZ/out" && ok "a live Ensembl VEP step is timed and its memory reads not measured, the engine being a container" || bad "sizing: the Ensembl VEP steps are absent or carry a memory figure"
python3 "$SZ/sizing_check.py" totals "$SZ/out" && ok "peak_memory_kib is the largest peak of the candidate's runs and the comparisons, beside the work disk's high-water mark, its growth and the case's wall minutes" || bad "sizing totals: $(python3 -c 'import json,sys; s=json.load(open(sys.argv[1]))["sizing"]; s.pop("steps"); print(s)' "$SZ/out/p1/c-four/report.json" 2>&1)"
python3 "$SZ/sizing_check.py" run "$SZ/out" && ok "run.json names the program and its version and lists the case's sizing block under case_sizes" || bad "run.json case_sizes differs from the report's sizing block"
STUB_CAND_RC=4 run_oracle "$T4" --out "$SZ/out-rc" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 ]] && python3 "$SZ/sizing_check.py" unmeasured "$SZ/out-rc" &&
    ok "a case whose candidate failed compares nothing: rows_compared reads not measured, the failed runs keep their figures and the run is red as before" || bad "sizing of a failed candidate: rc=$rc"
# The program keeps the tool's name whatever file VEP_FORMATS_COMPARATOR names, as a path or
# as a command on PATH, and the version is that file's digest: here a script that runs the tool.
mkdir -p "$SZ/bin"
printf '%s\n' '#!/bin/bash' "exec \"$COMPARATOR\" \"\$@\"" >"$SZ/bin/wrapped-comparator"
chmod +x "$SZ/bin/wrapped-comparator"
VEP_FORMATS_COMPARATOR="$SZ/bin/wrapped-comparator" run_oracle "$T4" --out "$SZ/out-named" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
PATH="$SZ/bin:$PATH" VEP_FORMATS_COMPARATOR=wrapped-comparator run_oracle "$T4" --out "$SZ/out-named-path" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
python3 "$SZ/sizing_check.py" program "$SZ/out-named" "$SZ/bin/wrapped-comparator" && python3 "$SZ/sizing_check.py" program "$SZ/out-named-path" "$SZ/bin/wrapped-comparator" &&
    ok "with another file named by VEP_FORMATS_COMPARATOR, as a path or as a command on PATH, the program keeps the tool's name and the version carries that file's digest" || bad "sizing: the program's name or version under a comparator named by one of its two spellings"
# A run whose steps leave no record (here a python3 that runs each step without the timer) fails
# nothing: the figures read not measured and the verdict stands.
cat >"$SZ/python-no-timer" <<SH
#!/bin/bash
# The step timer is the one program handed over as -c text that imports threading.
if [[ "\${1:-}" == "-c" && "\${2:-}" == *"import json, os, signal, subprocess, sys, threading, time"* ]]; then
    shift 4
    exec "\$@"
fi
exec python3 "\$@"
SH
chmod +x "$SZ/python-no-timer"
PYTHON="$SZ/python-no-timer" run_oracle "$T4" --out "$SZ/out-none" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>"$SZ/none.err"
rc=$?
[[ "$rc" -eq 0 ]] && python3 "$SZ/sizing_check.py" no-records "$SZ/out-none" &&
    ok "a run whose steps left no record is green as before, with every measured figure reading not measured" || bad "sizing without records: rc=$rc $(tail -2 "$SZ/none.err" | tr '\n' ' ')"
# A SIGTERM that reaches a step timer while it starts its command is handed to the command. Here
# a python3 whose Popen signals the timer itself as the child starts, over a stub engine that
# sleeps: the engine is stopped and its step still leaves its record.
cat >"$SZ/python-signalled-start" <<'SH'
#!/bin/bash
if [[ "${1:-}" == "-c" && "${2:-}" == *"import json, os, signal, subprocess, sys, threading, time"* ]]; then
    program="import os, signal, subprocess
class SignalledPopen(subprocess.Popen):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        os.kill(os.getpid(), signal.SIGTERM)
subprocess.Popen = SignalledPopen
$2"
    shift 2
    exec python3 -c "$program" "$@"
fi
exec python3 "$@"
SH
chmod +x "$SZ/python-signalled-start"
STUB_CAND_SLEEP=30 STUB_CAND_PIDFILE="$SZ/signalled.pids" PYTHON="$SZ/python-signalled-start" run_oracle "$T4" --out "$SZ/out-signalled" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>&1
alive=0
[[ ! -e "$SZ/signalled.pids" ]] || while read -r spid; do kill -0 "$spid" 2>/dev/null && alive=$((alive + 1)); done <"$SZ/signalled.pids"
[[ "$alive" -eq 0 ]] && python3 "$SZ/sizing_check.py" signalled "$SZ/out-signalled" &&
    ok "a SIGTERM that reaches a step timer as it starts its command stops the command: no stub engine is left running, and each candidate step records exit 143" || bad "a signal at a step's start: $alive stub engine(s) left running, or a candidate step without its exit 143"
# A case that left no report for a candidate reads not measured under case_sizes, whatever an
# earlier run left at the report's path. The first candidate's engine kills the case's shell
# (under --jobs 2 the case runs in one), so neither candidate reports, and the second
# candidate's report in the directory is the earlier run's.
run_oracle "$T4" --out "$SZ/out-two" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins p1,p2 >/dev/null 2>&1
STUB_CAND_KILL_JOB=1 run_oracle "$T4" --out "$SZ/out-two" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins p1,p2 --jobs 2 >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 1 && -f "$SZ/out-two/p2/c-four/report.json" ]] && python3 -c 'import json, sys
sizes = json.load(open(sys.argv[1]))["case_sizes"]
sys.exit(0 if sizes == [{"case": "c-four", "candidate_pin": pin, "sizing": "not measured"} for pin in ("p1", "p2")] else 1)' "$SZ/out-two/run.json" &&
    ok "a case that left no report for a candidate reads not measured under case_sizes, an earlier run's report at the same path unread" || bad "case_sizes of a case that did not complete: rc=$rc $(python3 -c 'import json, sys
print([(e["candidate_pin"], e["sizing"] if isinstance(e["sizing"], str) else "a sizing block") for e in json.load(open(sys.argv[1]))["case_sizes"]])' "$SZ/out-two/run.json" 2>&1)"


# The statistics comparison of a case naming stats is a step like any other comparison.
STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/perl.txt" stats_run sizing-stats
rc=$?
[[ "$rc" -eq 0 ]] && python3 -c 'import json, sys
steps = json.load(open(sys.argv[1]))["sizing"]["steps"]
step = steps["comparison.stats"]
sys.exit(0 if set(steps) >= {"candidate.default", "comparison.default", "comparison.stats"} and step["exit"] == 0 and isinstance(step["peak_rss_kib"], int) else 1)' "$TMP/out-sizing-stats/p1/c-stats/report.json" &&
    ok "a case naming stats leaves the statistics comparison's figures under comparison.stats" || bad "sizing of a statistics comparison: rc=$rc"

echo
echo "FORMAT JOBS: --format-jobs runs a case's formats at once and reports as the serial run does"
expect_option_error "--format-jobs 0 is refused" "--fork, --perl-fork, --jobs and --format-jobs take a positive integer (got '0')" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs 0
expect_option_error "a non-numeric --format-jobs is refused" "take a positive integer (got 'two')" --candidates "$STUBDIR/vep" --candidate-pins p1 --format-jobs two
FJ="$TMP/fj"
mkdir -p "$FJ/serial" "$FJ/serial-rc" "$FJ/serial-refuse"
# same_report.py <serial report.json> <parallel report.json>: exit 0 when the two records agree
# byte for byte once the three clock fields and the sizing block, each run's own measurement, are
# dropped (write_report dumps sorted keys at indent 1).
cat >"$FJ/same_report.py" <<'PY'
import json, sys
records = [json.load(open(p)) for p in sys.argv[1:3]]
for record in records:
    for clock in ("run_utc", "started_utc", "finished_utc", "sizing"):
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
python3 "$FJ/same_report.py" "$FJ/serial/report.json" "$FJ/out/p1/c-four/report.json" && ok "report.json is byte-identical to the serial run's once the three clock fields and the sizing block are dropped" || bad "report.json differs from the serial run's beyond the clock fields"
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
[[ "$rc" -eq 1 && "$reasons" == "default: candidate exited 4: stub vep: failing on request"*"| tab: candidate exited 4"*"| warnings: the candidate's default run exited 4" ]] && python3 "$FJ/same_report.py" "$FJ/serial-rc/report.json" "$FJ/out-rc/p1/c-four/report.json" &&
    ok "candidate failures in every format job are reported per format in format order, then warnings, as the serial run reports them" || bad "candidate exit 4 under --format-jobs 4: rc=$rc reasons=$reasons"
# The REFUSAL section's green case, both engines refusing both formats, under format jobs.
TFR="$FJ/tf.tsv"
{
    echo "$HEADER"
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "--plugin LOEUF,file={corpus}/loeuf.tsv.gz,gene" default,tab "$IMAGE" offline 116.2 - - - - - - refusals/mixed -
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
    row c-fj-a fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf "$IMAGE" offline 116.2 - - - - - "$REF_CELL" - -
    row c-fj-b fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab "$IMAGE" offline 116.2 - - - - - "$REF_CELL" - -
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

echo
echo "COMPARATOR: VEP_FORMATS_COMPARATOR names the format comparator; unset or empty, the tree's release build runs"
# A copy of the runner in a tree of its own reads that tree's release build. Each stub
# comparator appends its name and its arguments to STUB_COMPARATOR_LOG, one line a start.
CT="$TMP/comparator-trees"
TOOL_REL="tools/vep-formats-compare/target/release/vep-formats-compare"
mkdir -p "$CT/built/scripts/concordance" "$CT/built/$(dirname "$TOOL_REL")" "$CT/unbuilt/scripts/concordance" "$CT/bin"
cp "$ORACLE" "$CT/built/scripts/concordance/feature_oracle.sh"
cp "$ORACLE" "$CT/unbuilt/scripts/concordance/feature_oracle.sh"
stub_comparator() { # stub_comparator <file> <name>
    cat >"$1" <<EOF
#!/usr/bin/env bash
printf '%s\n' "$2 \$*" >>"\$STUB_COMPARATOR_LOG"
exit 2
EOF
    chmod +x "$1"
}
stub_comparator "$CT/built/$TOOL_REL" tool
stub_comparator "$CT/bin/named-comparator" named
: >"$CT/not-executable"
: >"$CT/bin/inert-comparator"
TC="$TMP/tc.tsv"
{
    echo "$HEADER"
    row c-one fixture "$CORPUS/variants.vcf" "116_GRCh38=$CORPUS/json_cache" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TC"
export STUB_COMPARATOR_LOG="$CT/started"
# comparator_run <tree> <value> [--dry-run]: the one case through the runner of <tree> with
# the variable unset (value ""), empty (value `empty`) or set; stdout in $CT/run.out, stderr
# in $CT/run.err, the engines' starts in $CT/docker.log and $CT/cand.log.
comparator_run() {
    local tree="$1" value="$2"
    shift 2
    rm -rf "$CT/out" "$CT/docker.log" "$CT/cand.log"
    : >"$STUB_COMPARATOR_LOG"
    (
        case "$value" in
        "") unset VEP_FORMATS_COMPARATOR ;;
        empty) export VEP_FORMATS_COMPARATOR="" ;;
        *) export VEP_FORMATS_COMPARATOR="$value" ;;
        esac
        STUB_DOCKER_LOG="$CT/docker.log" STUB_CAND_LOG="$CT/cand.log" PATH="$CT/bin:$STUBDIR:$PATH" bash "$CT/$tree/scripts/concordance/feature_oracle.sh" \
            --cases "$TC" --data-dir "$DATA" --out "$CT/out" --candidates "$STUBDIR/vep" --candidate-pins p1 "$@" >"$CT/run.out" 2>"$CT/run.err"
    )
}
# compares_through <tree> <value> <file> <stub name>: a dry run prints <file> as the first
# word of the comparator command, and a real run starts that stub once, with the words the
# dry run printed after it.
compares_through() {
    local planned rc
    comparator_run "$1" "$2" --dry-run || return 1
    planned=$(cmd_second 'DRY-RUN c-one p1 default' "$CT/run.out")
    [[ "$planned" == "$3 --format default --expected "*" --report $CT/out/p1/c-one/default.json"* ]] || return 1
    comparator_run "$1" "$2"
    rc=$?
    [[ "$rc" -eq 1 && "$(<"$STUB_COMPARATOR_LOG")" == "$4 ${planned#* }" ]]
}
for value in "" empty; do
    compares_through built "$value" "$CT/built/$TOOL_REL" tool &&
        ok "VEP_FORMATS_COMPARATOR ${value:-unset}: the release build of the runner's own tree compares, started with the comparison's arguments alone" ||
        bad "VEP_FORMATS_COMPARATOR ${value:-unset}: started [$(<"$STUB_COMPARATOR_LOG")] $(head -2 "$CT/run.err")"
done
for value in "$CT/bin/named-comparator" named-comparator; do
    compares_through unbuilt "$value" "$CT/bin/named-comparator" named &&
        ok "VEP_FORMATS_COMPARATOR=${value##*/} ($([[ "$value" == */* ]] && echo a path || echo a bare name on PATH)): the comparator named compares, in a tree that has no build" ||
        bad "VEP_FORMATS_COMPARATOR=$value: started [$(<"$STUB_COMPARATOR_LOG")] $(head -2 "$CT/run.err")"
done
reasons=$(reasons_of "$CT/out/p1/c-one/report.json")
[[ "$reasons" == "default: comparator exited 2: no stderr" ]] && ok "a comparator that exits 2 without a word is red for its format" || bad "stub comparator exit 2: $reasons"
comparator_refusal() { # comparator_refusal <tree> <the comparator the line names>
    printf '%s' "ERROR: [feature_oracle] the format comparator $2 is not an executable: build the tool with 'cargo build --release --locked --manifest-path $CT/$1/tools/vep-formats-compare/Cargo.toml' or name a built one in VEP_FORMATS_COMPARATOR"
}
refused() { # refused <tree> <value> <the comparator the line names>: the run stopped on that line before anything started
    comparator_run "$1" "$2"
    local rc=$?
    [[ "$rc" -eq 1 && "$(<"$CT/run.err")" == "$(comparator_refusal "$1" "$3")" && ! -s "$STUB_COMPARATOR_LOG" && ! -e "$CT/docker.log" && ! -e "$CT/cand.log" && ! -e "$CT/out" ]]
}
for value in "" empty; do
    refused unbuilt "$value" "$CT/unbuilt/$TOOL_REL" &&
        ok "VEP_FORMATS_COMPARATOR ${value:-unset} in a tree without the build: one line names the build command and the variable, exit 1, no engine runs and nothing is written" ||
        bad "unbuilt tree, VEP_FORMATS_COMPARATOR ${value:-unset}: $(head -2 "$CT/run.err")"
done
for value in "$CT/not-executable" "$CT" inert-comparator "$CT/absent"; do
    refused built "$value" "$value" &&
        ok "VEP_FORMATS_COMPARATOR=${value##*/} (not an executable file): refused with the same line, exit 1, and the tree's own build does not run in its place" ||
        bad "VEP_FORMATS_COMPARATOR=$value: $(head -2 "$CT/run.err")"
done
# The statistics comparator is found the same ways, where a selected case names stats: the
# one case with stats beside default, through a format comparator that is named.
TCS="$TMP/tcs.tsv"
{
    echo "$HEADER"
    row c-one fixture "$CORPUS/variants.vcf" "116_GRCh38=$CORPUS/json_cache" - "$CORPUS_FLAGS" default,stats "$IMAGE" offline 116.2 - - - - - - - -
} >"$TCS"
STATS_REL="tools/vep-stats-compare/target/release/vep-stats-compare"
mkdir -p "$CT/built/$(dirname "$STATS_REL")"
cp "$STUBDIR/stats-compare" "$CT/built/$STATS_REL"
cp "$STUBDIR/stats-compare" "$CT/bin/named-statistics"
: >"$CT/bin/inert-statistics"
# stats_comparator_run <tree> <value> [--dry-run]: the stats case through the runner of <tree>
# with VEP_STATS_COMPARATOR unset (value ""), empty (value `empty`) or set.
stats_comparator_run() {
    local tree="$1" value="$2"
    shift 2
    rm -rf "$CT/out" "$CT/docker.log" "$CT/cand.log"
    (
        case "$value" in
        "") unset VEP_STATS_COMPARATOR ;;
        empty) export VEP_STATS_COMPARATOR="" ;;
        *) export VEP_STATS_COMPARATOR="$value" ;;
        esac
        VEP_FORMATS_COMPARATOR="$CT/bin/named-comparator" STUB_DOCKER_LOG="$CT/docker.log" STUB_CAND_LOG="$CT/cand.log" PATH="$CT/bin:$STUBDIR:$PATH" bash "$CT/$tree/scripts/concordance/feature_oracle.sh" \
            --cases "$TCS" --data-dir "$DATA" --out "$CT/out" --candidates "$STUBDIR/vep" --candidate-pins p1 "$@" >"$CT/run.out" 2>"$CT/run.err"
    )
}
for value in "" empty; do
    stats_comparator_run built "$value" --dry-run && [[ "$(cmd_after 'DRY-RUN c-one p1 stats' "$CT/run.out")" == "$CT/built/$STATS_REL --expected "*" --report $CT/out/p1/c-one/stats.json" ]] &&
        ok "VEP_STATS_COMPARATOR ${value:-unset}: the release build of the runner's own tree compares the statistics" ||
        bad "VEP_STATS_COMPARATOR ${value:-unset}: $(cmd_after 'DRY-RUN c-one p1 stats' "$CT/run.out") $(head -2 "$CT/run.err")"
done
for value in "$CT/bin/named-statistics" named-statistics; do
    stats_comparator_run unbuilt "$value" --dry-run && [[ "$(cmd_after 'DRY-RUN c-one p1 stats' "$CT/run.out")" == "$CT/bin/named-statistics --expected "* ]] &&
        ok "VEP_STATS_COMPARATOR=${value##*/} ($([[ "$value" == */* ]] && echo a path || echo a bare name on PATH)): the comparator named compares the statistics, in a tree that has no build" ||
        bad "VEP_STATS_COMPARATOR=$value: $(cmd_after 'DRY-RUN c-one p1 stats' "$CT/run.out") $(head -2 "$CT/run.err")"
done
# stats_refused <tree> <value> <the comparator the line names>: the run stopped on that line before anything started.
stats_refused() {
    stats_comparator_run "$1" "$2"
    local rc=$?
    [[ "$rc" -eq 1 && "$(<"$CT/run.err")" == "ERROR: [feature_oracle] the statistics comparator $3 is not an executable: build the tool with 'cargo build --release --locked --manifest-path $CT/$1/tools/vep-stats-compare/Cargo.toml' or name a built one in VEP_STATS_COMPARATOR" && ! -s "$STUB_COMPARATOR_LOG" && ! -e "$CT/docker.log" && ! -e "$CT/cand.log" && ! -e "$CT/out" ]]
}
: >"$STUB_COMPARATOR_LOG"
for value in "" empty; do
    stats_refused unbuilt "$value" "$CT/unbuilt/$STATS_REL" &&
        ok "VEP_STATS_COMPARATOR ${value:-unset} in a tree without the build: one line names the build command and the variable, exit 1, no engine runs and nothing is written" ||
        bad "unbuilt tree, VEP_STATS_COMPARATOR ${value:-unset}: $(head -2 "$CT/run.err")"
done
for value in "$CT/not-executable" "$CT" inert-statistics "$CT/absent"; do
    stats_refused built "$value" "$value" &&
        ok "VEP_STATS_COMPARATOR=${value##*/} (not an executable file): refused with the same line, exit 1, and the tree's own build does not run in its place" ||
        bad "VEP_STATS_COMPARATOR=$value: $(head -2 "$CT/run.err")"
done
unset STUB_COMPARATOR_LOG
# The tool through the runner, every format of one case whose candidate differs in one field
# of its default output: the case is red on that field alone, and a comparison that ran
# leaves its standard error empty.
TQ="$TMP/tq.tsv"
{
    echo "$HEADER"
    row c-quiet fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,tab,vcf,json "$IMAGE" offline 116.2 - - - - - - - -
} >"$TQ"
STUB_CAND_MUTATE=change_impact run_oracle "$TQ" --out "$TMP/out-quiet" --work-dir "$TMP/work-quiet" --candidates "$STUBDIR/vep" --candidate-pins p1 --keep-outputs >/dev/null 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-quiet/p1/c-quiet/report.json")
[[ "$rc" -eq 1 && "$reasons" == "default: field Extra: 1 differing entries | default: field IMPACT: 1 differing entries" ]] && ok "a case of four formats differing in one field of its default output is red on that field alone" || bad "four formats, one field: rc=$rc reasons=$reasons"
quiet=0
for f in default tab vcf json; do
    [[ -e "$TMP/work-quiet/cases/c-quiet/candidates/p1/$f.compare.err" && ! -s "$TMP/work-quiet/cases/c-quiet/candidates/p1/$f.compare.err" ]] && quiet=$((quiet + 1))
done
[[ "$quiet" -eq 4 ]] && ok "the comparator writes nothing on the standard error of a comparison that ran" || bad "$quiet of 4 comparisons are silent on the standard error"

echo
echo "DISCOVER: a discovery pass reports the classes that applied and gives no verdict"
# A discovery pass (--discover) reports which classes apply and gives no verdict. The case
# compares the default format and the warnings over a copy of the corpus whose manifest keeps
# the four candidate-only keys it documents under reference_skipped_record, which applies in
# the default format, and the reference's one warning line carries a forked child's pid
# prefix, which fork_warning_pid_prefix removes. Its cell names fork_warning_pid_prefix,
# which applies, and order_only, which applies nowhere, with a count no format would meet.
TD="$TMP/discover.tsv"
{
    echo "$HEADER"
    row c-disc fixture discover/variants.vcf "$CORPUS_CACHE" - "$CORPUS_FLAGS" default,warnings "$IMAGE" offline 116.2 order_only=7,fork_warning_pid_prefix - - - - - - -
    row c-plain fixture nomanifest/variants.vcf "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 - - - - - - - -
    row c-refuse fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - "$TMP/refusal-disc" -
    row c-tool fixture "$CORPUS_INPUT" "$CORPUS_CACHE" - - default "$IMAGE" offline 116.2 order_only - - - - - - "filter_vep:--filter SYMBOL"
    row c-odd fixture discover-odd/variants.vcf "$CORPUS_CACHE" - "$CORPUS_FLAGS" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TD"
mkdir -p "$DATA/discover" "$TMP/refusal-disc"
cp "$CORPUS/variants.vcf" "$DATA/discover/variants.vcf"
python3 -c 'import json, sys
manifest = json.load(open(sys.argv[1]))
kept = {key: manifest[key] for key in ("release", "assembly", "flags", "vep_rs_only_tuples")}
json.dump(kept, open(sys.argv[2], "w"))' "$CORPUS/manifest.json" "$DATA/discover/manifest.json"
echo "refused by design" >"$TMP/refusal-disc/perl.txt"
echo "refused by design" >"$TMP/refusal-disc/candidate.txt"
DISC_WARNING="WARNING: line 12 skipped (21 30165800 del ...): Skipped record"
# Every run of c-disc below reads the two warnings files and takes these arguments.
export STUB_PERL_WARNINGS="WARNING: 4242 : $DISC_WARNING" STUB_CAND_WARNINGS="$DISC_WARNING"
DISC_RUN=(run_oracle "$TD" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-disc --out)
"${DISC_RUN[@]}" "$TMP/out-disc" --discover >"$TMP/disc.out" 2>&1
rc=$?
[[ "$rc" -eq 0 ]] && ok "a discovery pass exits 0 though a class its cell names applied nowhere" || bad "discovery exit: rc=$rc $(tail -3 "$TMP/disc.out")"
[[ "$(sed -n 1p "$TMP/out-disc/status.tsv")" == "case	candidate_pin	verdict	report_digest	run_utc" && "$(sed -n 2p "$TMP/out-disc/status.tsv" | cut -f1-3)" == "c-disc	p1	discovered" ]] && ! grep -qw -e green -e red "$TMP/out-disc/status.tsv" &&
    ok "a discovered case's status is the word discovered, and the table carries no verdict word" || bad "discovery status: $(cat "$TMP/out-disc/status.tsv" 2>/dev/null)"
want=$(printf 'format\tclass\trows\tnamed\tapplied\nwarnings\tfork_warning_pid_prefix\t1\tyes\tyes\ndefault\treference_skipped_record\t4\tno\tyes\n-\torder_only\t0\tyes\tno')
[[ "$(cat "$TMP/out-disc/p1/c-disc/discovery.tsv" 2>/dev/null)" == "$want" ]] &&
    ok "the case's discovery.tsv names each class that applied with its format and rows, whether the cell names it, and the named class that applied nowhere" || bad "per-case discovery.tsv: $(cat "$TMP/out-disc/p1/c-disc/discovery.tsv" 2>/dev/null)"
want=$(printf 'case\tcandidate_pin\tclass\tapplied_formats\trows\nc-disc\tp1\tfork_warning_pid_prefix\twarnings\t1\nc-disc\tp1\treference_skipped_record\tdefault\t4')
[[ "$(cat "$TMP/out-disc/discovery.tsv" 2>/dev/null)" == "$want" ]] &&
    ok "the run's discovery.tsv carries one line per class that applied, with its formats and rows" || bad "run discovery.tsv: $(cat "$TMP/out-disc/discovery.tsv" 2>/dev/null)"
grep -q 'discovery c-disc p1: applied: fork_warning_pid_prefix,reference_skipped_record; named and applied nowhere: order_only; applied and not named: reference_skipped_record; residual differences: 0$' "$TMP/disc.out" &&
    ok "the summary line states the three facts: the classes applied, the named class that applied nowhere, the applied class the cell does not name" || bad "discovery summary: $(grep 'discovery c-disc' "$TMP/disc.out")"
facts=$(python3 -c 'import json, sys
r = json.load(open(sys.argv[1]))
d = r["discovery"]
print(r.get("mode"), r.get("status"), "verdict" in r, any("verdict" in e for e in r["formats"].values()), sorted(d["applied"]), d["applied"]["reference_skipped_record"], d["named_applied_nowhere"], d["applied_not_named"], d["residual"]["count"])' "$TMP/out-disc/p1/c-disc/report.json" 2>&1)
[[ "$facts" == "discover discovered False False ['fork_warning_pid_prefix', 'reference_skipped_record'] {'default': 4} ['order_only'] ['reference_skipped_record'] 0" ]] &&
    ok "report.json carries the mode, the status and the three facts, and no verdict for the case or a format" || bad "discovery report.json: $facts"
run_mode=$(python3 -c 'import json, sys
r = json.load(open(sys.argv[1]))
print(r.get("mode"), r.get("discovered"), r.get("incomplete"), r.get("discover_classes"), "red" in r)' "$TMP/out-disc/run.json" 2>&1)
[[ "$run_mode" == "discover 1 0 None False" ]] && ok "run.json carries the mode discover, the discovered and incomplete counts, and no red count" || bad "discovery run.json: $run_mode"
asked=$(python3 -c 'import json, sys
print(",".join(json.load(open(sys.argv[1]))["masks"]["requested"]))' "$TMP/out-disc/p1/c-disc/default.json" 2>&1)
builtin_n=$(awk -F'"' '/^BUILTIN_MASKS="/ {print split($2, w, " "); exit}' "$ORACLE")
opt_in=$(awk -F'"' '/^DISCOVERY_OPT_IN_MASKS="/ {print $2; exit}' "$ORACLE")
asked_n=$(printf '%s' "$asked" | tr ',' '\n' | grep -c .)
# The registry asked for: every built-in class but the five asked for by name alone, of which
# the case's manifest documents one (reference_skipped_record).
[[ "$asked_n" -eq $((builtin_n - 4)) && ",$asked," == *",reference_skipped_record,"* && ",$asked," != *",breakend_mate_local_read,"* && ",$asked," == *",order_only,"* && ",$asked," != *",regulatory_feature_rows,"* && ",$asked," != *",individual_allele_order,"* && ",$asked," != *",bnd_own_chromosome_rows_lost,"* && "$(printf '%s' "$opt_in" | wc -w | tr -d ' ')" -eq 5 && "$asked" != *"="* ]] &&
    ok "the comparator is asked for the registry in place of the cell: every built-in class but those asked for by name alone, the manifest's own classes, the cell's without their counts ($asked_n classes)" || bad "discovery request: $asked_n of $builtin_n classes: $asked"
# The same case as a verdict run: red, for the class named and applied in no compared format.
"${DISC_RUN[@]}" "$TMP/out-disc-verdict" >"$TMP/disc-verdict.out" 2>&1
rc=$?
reasons=$(reasons_of "$TMP/out-disc-verdict/p1/c-disc/report.json")
[[ "$rc" -eq 1 && "$reasons" == "mask order_only is requested and applied in no compared format" && "$(sed -n 2p "$TMP/out-disc-verdict/status.tsv" | cut -f3)" == "red" && ! -e "$TMP/out-disc-verdict/discovery.tsv" && ! -e "$TMP/out-disc-verdict/p1/c-disc/discovery.tsv" ]] &&
    [[ "$(python3 -c 'import json, sys; r = json.load(open(sys.argv[1])); print(r["mode"], r["red"])' "$TMP/out-disc-verdict/run.json")" == "verdict 1" ]] &&
    ok "the same case without --discover is red for the named class that applied nowhere, writes no discovery table, and its run.json reads mode verdict" || bad "verdict run of the discovery case: rc=$rc reasons=$reasons"
# --discover-classes narrows the registry to the listed classes; the cell's own stay.
"${DISC_RUN[@]}" "$TMP/out-disc-narrow" --discover --discover-classes reference_skipped_record,splice_family_swap >"$TMP/disc-narrow.out" 2>&1
rc=$?
asked=$(python3 -c 'import json, sys
print(",".join(json.load(open(sys.argv[1]))["masks"]["requested"]))' "$TMP/out-disc-narrow/p1/c-disc/default.json" 2>&1)
[[ "$rc" -eq 0 && "$asked" == "fork_warning_pid_prefix,order_only,reference_skipped_record,splice_family_swap" && "$(sed -n 3p "$TMP/out-disc-narrow/discovery.tsv")" == "c-disc	p1	fork_warning_pid_prefix	warnings	1" ]] &&
    [[ "$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["discover_classes"])' "$TMP/out-disc-narrow/run.json")" == "['reference_skipped_record', 'splice_family_swap']" ]] &&
    ok "--discover-classes narrows the request to the listed classes beside the cell's own, and run.json records the list" || bad "narrowed discovery: rc=$rc asked=$asked $(cat "$TMP/out-disc-narrow/discovery.tsv" 2>/dev/null)"
"${DISC_RUN[@]}" "$TMP/out-disc-one" --discover --discover-classes splice_family_swap >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(sed -n '2,$p' "$TMP/out-disc-one/discovery.tsv")" == "c-disc	p1	fork_warning_pid_prefix	warnings	1" ]] &&
    ok "a class the narrowed list leaves out is not asked for, so the pass does not report it" || bad "narrowed to one class: rc=$rc $(cat "$TMP/out-disc-one/discovery.tsv" 2>/dev/null)"
dry "$TD" "$TMP/disc-dry.out" --only c-disc --discover --discover-classes reference_skipped_record
rc=$?
[[ "$rc" -eq 0 && "$(cmd_second "DRY-RUN c-disc stub default" "$TMP/disc-dry.out")" == *'--masks reference_skipped_record\,order_only\,fork_warning_pid_prefix'* && "$(cmd_after "DRY-RUN c-disc stub warnings" "$TMP/disc-dry.out")" == *'--masks reference_skipped_record\,order_only\,fork_warning_pid_prefix' ]] &&
    ok "a dry run of a discovery pass prints the comparator commands with the classes asked for, the list's first" || bad "discovery dry run: rc=$rc $(cmd_second "DRY-RUN c-disc stub default" "$TMP/disc-dry.out")"
expect_error "--discover-classes without --discover is an error" "--discover-classes names the classes a discovery pass asks for; pass --discover with it" "$TD" --only c-disc --discover-classes order_only
expect_error "an empty --discover-classes is an error, not a pass over the registry" "--discover-classes '' is not <class>\[,<class>...\]" "$TD" --only c-disc --discover --discover-classes ""
expect_error "a --discover-classes item with a count is an error" "is not <class>\[,<class>...\]" "$TD" --only c-disc --discover --discover-classes order_only=3
expect_error "a --discover-classes name that is neither built-in nor a selected case's manifest class is an error" "--discover-classes names 'no_such_class', which is neither a built-in class nor a class a selected case's manifest.json documents" "$TD" --only c-disc --discover --discover-classes order_only,no_such_class
expect_error "discovery.tsv is no candidate pin: the run writes that name under --out" "candidate pin 'discovery.tsv' names a file or directory the run writes under --out" "$TD" --only c-disc --candidate-pins discovery.tsv
# A residual difference no class covers is counted and listed, and does not fail the pass.
STUB_CAND_MUTATE=change_impact "${DISC_RUN[@]}" "$TMP/out-disc-residual" --discover >"$TMP/disc-residual.out" 2>&1
rc=$?
residual=$(python3 -c 'import json, sys
r = json.load(open(sys.argv[1]))
print(r["status"], r["discovery"]["residual"]["count"], r["discovery"]["residual"]["reasons"], r["formats"]["default"]["residual"])' "$TMP/out-disc-residual/p1/c-disc/report.json" 2>&1)
[[ "$rc" -eq 0 && "$residual" == "discovered 2 ['default: field Extra: 1 differing entries', 'default: field IMPACT: 1 differing entries'] ['field Extra: 1 differing entries', 'field IMPACT: 1 differing entries']" ]] && grep -q 'discovery c-disc p1: .*residual differences: 2$' "$TMP/disc-residual.out" &&
    ok "a difference no class covers is counted and listed in the report and the summary, and the pass still exits 0" || bad "residual difference: rc=$rc $residual"
# A case whose engine fails did not run: its status is incomplete, never a verdict word, it
# leaves no line in the run's table, and the pass exits 1.
STUB_CAND_RC=3 "${DISC_RUN[@]}" "$TMP/out-disc-failed" --discover >"$TMP/disc-failed.out" 2>&1
rc=$?
[[ "$rc" -eq 1 && "$(sed -n 2p "$TMP/out-disc-failed/status.tsv" | cut -f1-3)" == "c-disc	p1	incomplete" && "$(wc -l <"$TMP/out-disc-failed/discovery.tsv" | tr -d ' ')" -eq 1 && ! -e "$TMP/out-disc-failed/p1/c-disc/discovery.tsv" ]] && ! grep -qw -e green -e red "$TMP/out-disc-failed/status.tsv" &&
    grep -q 'discovery c-disc p1: did not run: default: candidate exited 3' "$TMP/disc-failed.out" &&
    [[ "$(python3 -c 'import json, sys; r = json.load(open(sys.argv[1])); print(r["discovered"], r["incomplete"])' "$TMP/out-disc-failed/run.json")" == "0 1" ]] &&
    ok "a case whose engine fails is incomplete, leaves no line in the run's table, and the pass exits 1" || bad "incomplete discovery: rc=$rc $(cat "$TMP/out-disc-failed/status.tsv" 2>/dev/null) $(tail -2 "$TMP/disc-failed.out")"
unset STUB_PERL_WARNINGS STUB_CAND_WARNINGS
# A case no class applies to is one line of the run's table whose class is -; a refusal row
# and a tool row are asked for their own cell alone.
run_oracle "$TD" --out "$TMP/out-disc-rest" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-plain --discover >"$TMP/disc-rest.out" 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(sed -n '2,$p' "$TMP/out-disc-rest/discovery.tsv")" == "c-plain	p1	-	-	-" && "$(cat "$TMP/out-disc-rest/p1/c-plain/discovery.tsv")" == "format	class	rows	named	applied" ]] &&
    ok "a case no class applies to is one line of the run's table whose class is -, and its own table is the header alone" || bad "discovery of a plain case: rc=$rc $(cat "$TMP/out-disc-rest/discovery.tsv" 2>/dev/null)"
# The run's table and the summary lines read a pin as text too: the pins 1 and 1.0 are two candidates.
run_oracle "$TD" --out "$TMP/out-disc-pins" --candidates "$STUBDIR/vep,$STUBDIR/vep2" --candidate-pins 1,1.0 --only c-plain --discover >"$TMP/disc-pins.out" 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(sed -n '2,$p' "$TMP/out-disc-pins/discovery.tsv" | tr '\n\t' '; ')" == "c-plain 1 - - -;c-plain 1.0 - - -;" && "$(grep -c '^applied: ' "$TMP/disc-pins.out")" -eq 0 ]] &&
    ok "a discovery pass over two candidates pinned 1 and 1.0 gives each its own line of the run's table and its own summary line, and exits 0" || bad "discovery under numeral pins: rc=$rc $(sed -n '2,$p' "$TMP/out-disc-pins/discovery.tsv" 2>/dev/null | tr '\n\t' '; ')"
dry "$TD" "$TMP/disc-dry2.out" --only c-refuse,c-tool --discover
rc=$?
[[ "$rc" -eq 0 && "$(grep -c -- '--masks' "$TMP/disc-dry2.out")" -eq 1 && "$(cmd_second "DRY-RUN c-tool stub default" "$TMP/disc-dry2.out")" == *"--masks order_only" ]] &&
    ok "a refusal row is asked for nothing and a tool row for its own cell alone" || bad "discovery of refusal and tool rows: rc=$rc $(grep -- '--masks' "$TMP/disc-dry2.out")"
# A manifest item whose class is written as a list names no class: it is passed over, and the
# manifest's other classes are still asked for (reference_skipped_record is asked for by name
# alone, so the request holds it only because the manifest documents it).
mkdir -p "$DATA/discover-odd"
cp "$CORPUS/variants.vcf" "$DATA/discover-odd/variants.vcf"
python3 -c 'import json, sys
manifest = json.load(open(sys.argv[1]))
manifest["field_divergences"] = [{"expected_divergence": ["a_list"], "second_expected_divergence": "class_of_the_odd_item"}]
json.dump(manifest, open(sys.argv[2], "w"))' "$DATA/discover/manifest.json" "$DATA/discover-odd/manifest.json"
dry "$TD" "$TMP/disc-odd.out" --only c-odd --discover
rc=$?
odd=$(cmd_second "DRY-RUN c-odd stub default" "$TMP/disc-odd.out")
# The request is read class by class: plugin_value_on_reference_skipped_record, a built-in
# class every registry carries, ends in the name of the class looked for.
odd_asked=$(sed -n 's/.* --masks \([^ ]*\).*/\1/p' <<<"$odd" | tr -d '\\')
[[ "$rc" -eq 0 && ",$odd_asked," == *",reference_skipped_record,"* && "$odd" == *"class_of_the_odd_item"* && "$odd" != *"a_list"* ]] && ! grep -q Traceback "$TMP/disc-odd.out.err" &&
    ok "a manifest class written as a list is passed over and the manifest's other classes are still asked for" || bad "odd manifest class: rc=$rc $(head -3 "$TMP/disc-odd.out.err")"
# A verdict run into the directory of an earlier discovery pass leaves no discovery table.
export STUB_PERL_WARNINGS="WARNING: 4242 : $DISC_WARNING" STUB_CAND_WARNINGS="$DISC_WARNING"
"${DISC_RUN[@]}" "$TMP/out-disc" >/dev/null 2>&1
[[ ! -e "$TMP/out-disc/discovery.tsv" && "$(sed -n 2p "$TMP/out-disc/status.tsv" | cut -f3)" == "red" ]] &&
    ok "a verdict run into an earlier discovery pass's directory removes the run's discovery table" || bad "stale discovery.tsv survived a verdict run"
unset STUB_PERL_WARNINGS STUB_CAND_WARNINGS
# The sequence-edge class is decided from the reference both engines read: a discovery pass
# asks for it with the registry, so its comparison of a case that gives a FASTA names that
# FASTA, decompressed as the engines read it, and the candidate's cache once.
TDE="$TMP/discover-edge.tsv"
{
    echo "$HEADER"
    row c-edge-unnamed fixture nomanifest/variants.vcf 116_GRCh37=tests/golden/116/GRCh37-hgvs/json_cache tests/golden/116/GRCh37-hgvs/reference.fa.gz --hgvs default "$IMAGE" offline 116.2 order_only - - - - - - -
} >"$TDE"
dry "$TDE" "$TMP/e-edge-disc.out" --discover
rc=$?
bed=$(cmd_second 'DRY-RUN c-edge-unnamed stub default' "$TMP/e-edge-disc.out")
[[ "$rc" -eq 0 && "$bed" == *"--masks "*"hgvs_shift_sequence_edge_window"* && "$bed" == *" --json-cache $REPO/tests/golden/116/GRCh37-hgvs/json_cache "* && "$bed" == *" --reference-fasta $TMP/never/work/fasta/1/reference.fa"* ]] &&
    [[ "$(grep -o -- '--json-cache ' <<<"$bed" | wc -l | tr -d ' ')" -eq 1 ]] &&
    ok "a discovery pass, which asks for the sequence-edge class with the registry, hands the comparison the case's FASTA and names the cache once" || bad "sequence-edge reference in a discovery pass: rc=$rc $bed"
dry "$TDE" "$TMP/e-edge-verdict.out"
beu=$(cmd_second 'DRY-RUN c-edge-unnamed stub default' "$TMP/e-edge-verdict.out")
[[ "$beu" == *"--masks order_only" ]] &&
    ok "the same case as a verdict run asks for its cell alone and carries no --reference-fasta" || bad "sequence-edge reference leaked into a verdict run: $beu"

# A pooled case compares unnamed, and the comparator refuses an unnamed comparison asked for a
# class the lists of several runs over the input document. The manifest here lists two runs
# over the input, its own documenting reference_skipped_record: a discovery pass asks the
# pooled case for the registry less that class, and for the class where the cell names it.
mkdir -p "$DATA/discover-pooled"
cp "$CORPUS/variants.vcf" "$DATA/discover-pooled/variants.vcf"
python3 -c 'import json, sys
manifest = json.load(open(sys.argv[1]))
manifest["flag_sets"] = [{"name": "second", "flags": manifest["flags"] + ["--canonical"], "input": "variants.vcf"}]
json.dump(manifest, open(sys.argv[2], "w"))' "$DATA/discover/manifest.json" "$DATA/discover-pooled/manifest.json"
TDP="$TMP/discover-pooled.tsv"
{
    echo "$HEADER"
    row c-pooled fixture discover-pooled/variants.vcf "$CORPUS_CACHE" - "$CORPUS_FLAGS --symbol" default "$IMAGE" offline 116.2 - - - - - - - -
    row c-pooled-named fixture discover-pooled/variants.vcf "$CORPUS_CACHE" - "$CORPUS_FLAGS --symbol" default "$IMAGE" offline 116.2 reference_skipped_record - - - - - - -
    row c-pooled-one fixture discover/variants.vcf "$CORPUS_CACHE" - "$CORPUS_FLAGS --symbol" default "$IMAGE" offline 116.2 - - - - - - - -
} >"$TDP"
printf 'case\trun\treason\nc-pooled\tpooled\tno run of the manifest carries --symbol\nc-pooled-named\tpooled\tno run of the manifest carries --symbol\nc-pooled-one\tpooled\tthe one run of the manifest carries no --symbol\n' >"$TMP/discover-pooled-runs.tsv"
run_oracle "$TDP" --case-runs "$TMP/discover-pooled-runs.tsv" --out "$TMP/out-disc-pooled" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-pooled --discover >"$TMP/disc-pooled.out" 2>&1
rc=$?
asked=$(python3 -c 'import json, sys
print(",".join(json.load(open(sys.argv[1]))["masks"]["requested"]))' "$TMP/out-disc-pooled/p1/c-pooled/default.json" 2>&1)
[[ "$rc" -eq 0 && "$(sed -n 2p "$TMP/out-disc-pooled/status.tsv" | cut -f1-3)" == "c-pooled	p1	discovered" && ",$asked," == *",order_only,"* && ",$asked," != *",reference_skipped_record,"* ]] &&
    ok "a discovery pass asks a pooled case for no class its manifest documents, so the unnamed comparison runs and the case is discovered" || bad "discovery of a pooled case: rc=$rc $(sed -n 2p "$TMP/out-disc-pooled/status.tsv" 2>/dev/null | cut -f1-3) $(tail -2 "$TMP/disc-pooled.out" | tr '\n' ' ')"
dry "$TDP" "$TMP/disc-pooled-dry.out" --case-runs "$TMP/discover-pooled-runs.tsv" --only c-pooled-named --discover
rc=$?
asked=$(cmd_second "DRY-RUN c-pooled-named stub default" "$TMP/disc-pooled-dry.out" | sed -n 's/.* --masks \([^ ]*\).*/\1/p' | tr -d '\\')
[[ "$rc" -eq 0 && ",$asked," == *",reference_skipped_record,"* ]] &&
    ok "a pooled case whose cell names a class its manifest documents is still asked for it" || bad "discovery of a pooled case that names the class: rc=$rc asked=$asked"
# On a manifest of one run over the input the comparator refuses no unnamed comparison: a
# pooled case there is asked for the manifest's classes, and the pass reports the one that
# applies.
run_oracle "$TDP" --case-runs "$TMP/discover-pooled-runs.tsv" --out "$TMP/out-disc-pooled-one" --candidates "$STUBDIR/vep" --candidate-pins p1 --only c-pooled-one --discover >"$TMP/disc-pooled-one.out" 2>&1
rc=$?
[[ "$rc" -eq 0 ]] && grep -q 'c-pooled-one: compares under no named run of its manifest, pooled by ' "$TMP/disc-pooled-one.out" &&
    grep -q 'discovery c-pooled-one p1: applied: reference_skipped_record; named and applied nowhere: -; applied and not named: reference_skipped_record; residual differences: 0$' "$TMP/disc-pooled-one.out" &&
    ok "a pooled case on a manifest of one run over its input is asked for the manifest's classes, and the pass reports the class that applies" || bad "discovery of a pooled case on one run: rc=$rc $(grep 'c-pooled-one' "$TMP/disc-pooled-one.out" | tail -2 | tr '\n' ' ')"
# The statistics comparison of a case naming stats runs in a discovery pass as in a verdict
# run, with no class asked of it: a difference between the two reports is a residual one.
STUB_PERL_STATS="$STATS_DIR/perl.txt" STUB_CAND_STATS="$STATS_DIR/cand-genes.txt" stats_run disc-stats --discover
rc=$?
residual=$(python3 -c 'import json, sys
r = json.load(open(sys.argv[1]))
print(r["status"], "stats: [General statistics]: 1 difference" in r["discovery"]["residual"]["reasons"], r["formats"]["stats"].get("residual"), "verdict" in r["formats"]["stats"])' "$TMP/out-disc-stats/p1/c-stats/report.json" 2>&1)
[[ "$rc" -eq 0 && "$residual" == "discovered True ['[General statistics]: 1 difference'] False" ]] &&
    ok "a discovery pass compares the statistics reports and counts their difference among the residual ones" || bad "statistics in a discovery pass: rc=$rc $residual"

echo
echo "THE NAMED RUN: a case on a manifest of several runs compares under one of them"
# The pooled rows of the run-setting table: the dry run says why no run is named.
grep -qxF "DRY-RUN c-rs-none run pooled, no run named by $RS_RUNS: no run of the manifest carries --symbol alone" "$TMP/d-rs.out" &&
    ok "a pooled row keeps the unnamed comparison and the dry run prints its reason" || bad "pooled line: $(grep '^DRY-RUN c-rs-none run' "$TMP/d-rs.out")"
grep -qxF "DRY-RUN c-rs-twin run twin-a, named by $RS_RUNS: the flags fit twin-a and twin-b alike, and the case gives no FASTA" "$TMP/d-rs.out" &&
    ok "a named row's run and reason are printed" || bad "named line: $(grep '^DRY-RUN c-rs-twin run' "$TMP/d-rs.out")"
grep -q 'runs of a manifest: 12 case(s) named by their flags, 1 by '"$RS_RUNS"', 4 pooled by a row of it' "$TMP/d-rs.out" &&
    ok "the run log counts the cases named by their flags, by the table and pooled" || bad "tally: $(grep 'runs of a manifest' "$TMP/d-rs.out")"
# Without the table the five cases are refused by name before anything is planned.
dry "$T_RS" "$TMP/d-rs-bare.out"
rc=$?
n=0
for rs_case in c-rs-none c-rs-twin c-rs-cell-no-fasta c-rs-no-fasta c-rs-fasta-dir; do
    grep -q "^ERROR: \[feature_oracle\] $rs_case: $RS/corpus/manifest.json lists 4 runs over variants.vcf (default, picked, twin-a, twin-b) and the flags are not exactly one of them" "$TMP/d-rs-bare.out.err" && n=$((n + 1))
done
[[ "$rc" -ne 0 && "$n" -eq 5 ]] && ! grep -q '^DRY-RUN' "$TMP/d-rs-bare.out" && [[ "$(grep -c '^ERROR' "$TMP/d-rs-bare.out.err")" -eq 6 ]] &&
    ok "without a row, each case that is no run or fits two is refused by name before any comparator is planned, and no other case is" || bad "unnamed cases: rc=$rc, $n of 5 named: $(head -2 "$TMP/d-rs-bare.out.err")"

# A corpus whose manifest writes a plugin directory, a path from the corpus directory, the
# two directory words and options under their first names; the cases write the same runs with
# another plugin directory, the corpus placeholder and other spellings of the options.
NR="$DATA/namedruns"
mkdir -p "$NR/corpus" "$TMP/nr" "$TMP/nr-default"
cp "$RS/corpus/variants.vcf" "$NR/corpus/variants.vcf"
cp "$RS/corpus/variants.vcf" "$NR/corpus/second.vcf"
cat >"$NR/corpus/manifest.json" <<'JSON'
{"flags": ["--symbol"],
 "flag_sets": [
  {"name": "plugin", "flags": ["--symbol", "--dir_plugins", "{plugins}", "--plugin", "Score,{plugins}/scores.tsv"]},
  {"name": "relative", "flags": ["--symbol", "--plugin", "Table,file=data/table.tsv"]},
  {"name": "generated", "flags": ["--symbol", "--custom", "file={generated}/sites.vcf.gz,short_name=sites"]},
  {"name": "versioned", "flags": ["--symbol", "--protein_version"]},
  {"name": "everything", "flags": ["--everything"]},
  {"name": "other-input", "input": "second.vcf", "flags": ["--symbol", "--pick"]}]}
JSON
nr_row() { # nr_row <case> <input under --data-dir> <flags cell> [masks cell]: a case over the bare candidate cache
    row "$1" fixture "$2" "116_GRCh38=runsetting/bare" - "$3" default "$IMAGE" offline 116.2 "${4:--}" - - - - - - -
}
T_NR="$TMP/nr/matched.tsv"
{
    echo "$HEADER"
    nr_row n-own namedruns/corpus/variants.vcf "--symbol"
    nr_row n-dir namedruns/corpus/variants.vcf "--symbol --dir_plugins /anywhere"
    nr_row n-dir-joined namedruns/corpus/variants.vcf "--dir_plugins=/anywhere --symbol"
    nr_row n-plugin namedruns/corpus/variants.vcf "--dir_plugins {corpus}/plugin_data --plugin Score,{corpus}/plugin_data/scores.tsv --symbol"
    nr_row n-relative namedruns/corpus/variants.vcf "--symbol --plugin Table,file={corpus}/data/table.tsv"
    nr_row n-generated namedruns/corpus/variants.vcf "--symbol --custom file={corpus}/generated/sites.vcf.gz,short_name=sites"
    nr_row n-spelled namedruns/corpus/variants.vcf "-SYMBOL --protein_v"
    nr_row n-other-name namedruns/corpus/variants.vcf "-e"
    nr_row n-manifest-only namedruns/corpus/variants.vcf "--symbol" transcript_annotator_insertion_window
} >"$T_NR"
dry "$T_NR" "$TMP/nr-matched.out"
rc=$?
[[ "$rc" -eq 0 ]] && ok "every case that differs from a run by words the match sees through plans without a row" || bad "matched table rc=$rc: $(head -3 "$TMP/nr-matched.out.err" | tr '\n' ' ')"
while IFS='|' read -r nr_case nr_run nr_label; do
    nr_cmd=$(cmd_second "DRY-RUN $nr_case stub default" "$TMP/nr-matched.out")
    [[ "$nr_cmd" == *" --manifest $NR/corpus/manifest.json --run $nr_run --reference-sequence absent"* ]] && ok "$nr_label" || bad "$nr_label: $nr_cmd"
done <<'CHECKS'
n-own|default|a case whose flags are the manifest's own run's compares under default
n-dir|default|a case whose flags add a plugin-directory word to a run's flags compares under that run
n-dir-joined|default|a plugin-directory word joined to its option by = is seen through as well
n-plugin|plugin|a plugin directory written under the corpus placeholder is the run that writes {plugins}
n-relative|relative|a file under the corpus placeholder is the run that writes the path from the corpus directory
n-generated|generated|a path under the corpus's generated directory is the run that writes {generated}
n-spelled|versioned|an option under one dash, in upper case or abbreviated is the run that writes its full name
n-other-name|everything|an option spelled by its other name is the run that writes its first name
n-manifest-only|default|a manifest-only class beside a manifest is planned
CHECKS
! grep -q '^DRY-RUN .* run ' "$TMP/nr-matched.out" && grep -q 'runs of a manifest: 9 case(s) named by their flags, 0 by no named-runs table, 0 pooled by a row of it' "$TMP/nr-matched.out" &&
    ok "a matched case needs no row and prints no row line" || bad "matched tally: $(grep 'runs of a manifest' "$TMP/nr-matched.out")"

# An added flag: no run of the manifest, so the case is refused by name without a row, compares
# under the run a row names, and keeps the unnamed comparison under a pooled row.
T_NA="$TMP/nr/added.tsv"
{
    echo "$HEADER"
    nr_row n-own namedruns/corpus/variants.vcf "--symbol"
    nr_row n-added namedruns/corpus/variants.vcf "--symbol --canonical"
    nr_row n-second namedruns/corpus/second.vcf "--symbol --canonical"
    nr_row n-bare nomanifest/variants.vcf "--symbol --canonical"
} >"$T_NA"
dry "$T_NA" "$TMP/nr-added.out"
rc=$?
[[ "$rc" -ne 0 ]] && ! grep -q '^DRY-RUN' "$TMP/nr-added.out" && [[ "$(grep -c '^ERROR' "$TMP/nr-added.out.err")" -eq 2 ]] &&
    grep -q "^ERROR: \[feature_oracle\] n-added: $NR/corpus/manifest.json lists 6 runs over variants.vcf (default, plugin, relative, generated, versioned, everything) and the flags are not exactly one of them, so the comparison would read every run's documented items as one list: name the case's run in a row of a named-runs table (--case-runs) (case, run, reason), or mark it pooled there with the reason no run is right" "$TMP/nr-added.out.err" &&
    ok "a case with an added flag and no row fails preflight by name; a case over an input one run annotates and a case without a manifest do not" || bad "added flag without a row: rc=$rc $(head -3 "$TMP/nr-added.out.err" | tr '\n' ' ')"
{
    row case run reason
    row n-added default "the own run plus a flag that changes no documented item"
} >"$TMP/nr/named.tsv"
dry "$T_NA" "$TMP/nr-named.out" --case-runs "$TMP/nr/named.tsv"
rc=$?
nr_cmd=$(cmd_second 'DRY-RUN n-added stub default' "$TMP/nr-named.out")
[[ "$rc" -eq 0 && "$nr_cmd" == *" --manifest $NR/corpus/manifest.json --run default --reference-sequence absent"* ]] &&
    grep -qxF "DRY-RUN n-added run default, named by $TMP/nr/named.tsv: the own run plus a flag that changes no documented item" "$TMP/nr-named.out" &&
    ok "a case with an added flag and a row compares under the run the row names" || bad "named row: rc=$rc $nr_cmd"
grep -q "runs of a manifest: 1 case(s) named by their flags, 1 by $TMP/nr/named.tsv, 0 pooled by a row of it" "$TMP/nr-named.out" &&
    ok "the tally tells the row from the match" || bad "named tally: $(grep 'runs of a manifest' "$TMP/nr-named.out")"
nr_cmd=$(cmd_second 'DRY-RUN n-second stub default' "$TMP/nr-named.out")
[[ "$nr_cmd" == *" --manifest $NR/corpus/manifest.json --reference-sequence absent"* && "$nr_cmd" != *"--run "* ]] &&
    ok "a case over an input that one run annotates, and whose flags are not that run's, compares unnamed: one list is read" || bad "one run over the input: $nr_cmd"
{
    row case run reason
    row n-added pooled "no run of the manifest carries --canonical"
    row n-bare pooled "no manifest stands beside its input"
} >"$TMP/nr/pooled.tsv"
dry "$T_NA" "$TMP/nr-pooled.out" --case-runs "$TMP/nr/pooled.tsv"
rc=$?
nr_cmd=$(cmd_second 'DRY-RUN n-added stub default' "$TMP/nr-pooled.out")
[[ "$rc" -eq 0 && "$nr_cmd" == *" --manifest $NR/corpus/manifest.json --reference-sequence absent"* && "$nr_cmd" != *"--run "* ]] &&
    grep -qxF "DRY-RUN n-added run pooled, no run named by $TMP/nr/pooled.tsv: no run of the manifest carries --canonical" "$TMP/nr-pooled.out" &&
    ok "a pooled row keeps the unnamed call and prints the reason" || bad "pooled row: rc=$rc $nr_cmd"
! grep -q '^DRY-RUN n-bare run ' "$TMP/nr-pooled.out" && grep -q "runs of a manifest: 1 case(s) named by their flags, 0 by $TMP/nr/pooled.tsv, 1 pooled by a row of it" "$TMP/nr-pooled.out" &&
    ok "a pooled row on a case without a manifest keeps nothing unnamed and is not read" || bad "pooled row without a manifest: $(grep 'runs of a manifest' "$TMP/nr-pooled.out")"
{
    row case run reason
    row n-own versioned "a named row is the case's run whatever the cell fits"
    row n-added default "the own run plus a flag"
} >"$TMP/nr/outranks.tsv"
dry "$T_NA" "$TMP/nr-outranks.out" --case-runs "$TMP/nr/outranks.tsv"
nr_cmd=$(cmd_second 'DRY-RUN n-own stub default' "$TMP/nr-outranks.out")
[[ "$nr_cmd" == *" --manifest $NR/corpus/manifest.json --run versioned "* ]] && ok "the table is read first: a named row outranks the run the flags fit" || bad "row over match: $nr_cmd"
# The default table: a file oracle_case_runs.tsv beside the case table, read without the option.
cp "$T_NA" "$TMP/nr-default/cases.tsv"
cp "$TMP/nr/named.tsv" "$TMP/nr-default/oracle_case_runs.tsv"
dry "$TMP/nr-default/cases.tsv" "$TMP/nr-default.out"
rc=$?
[[ "$rc" -eq 0 ]] && grep -qF "DRY-RUN n-added run default, named by $TMP/nr-default/oracle_case_runs.tsv: " "$TMP/nr-default.out" &&
    ok "a file oracle_case_runs.tsv beside the case table is read without --case-runs" || bad "default table: rc=$rc $(head -2 "$TMP/nr-default.out.err" | tr '\n' ' ')"

# Each refusal below stops the run before anything is planned (expect_error).
{
    row case run reason
    row n-added no-such-run "a run the manifest lacks"
} >"$TMP/nr/e.tsv"
expect_error "a row naming a run the manifest lacks is refused, with the manifest's runs" "n-added: $TMP/nr/e.tsv names the run 'no-such-run', and $NR/corpus/manifest.json lists no run of that name (its runs: default, plugin, relative, generated, versioned, everything, other-input)" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row n-added other-input "a run over another input"
} >"$TMP/nr/e.tsv"
expect_error "a row naming a run over another input is refused" "n-added: $TMP/nr/e.tsv names the run 'other-input', which annotates second.vcf, and the case's input is variants.vcf" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row n-own pooled "a matched case"
    row n-added default "the own run plus a flag"
} >"$TMP/nr/e.tsv"
expect_error "a pooled row on a case whose flags are a run of its manifest is refused" "n-own: $TMP/nr/e.tsv marks the case pooled, and its flags are the run 'default' of $NR/corpus/manifest.json" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row n-added default "the own run plus a flag"
    row n-bare default "a run where no manifest is"
} >"$TMP/nr/e.tsv"
expect_error "a named run on a case without a manifest beside its input is refused" "n-bare: $TMP/nr/e.tsv names the run 'default' for the case, which compares against no manifest" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row n-added default "the own run plus a flag"
    row n-added pooled "a second row"
} >"$TMP/nr/e.tsv"
expect_error "a second row for a case is refused" "$TMP/nr/e.tsv:3: case 'n-added' has a second row" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row n-added default "the own run plus a flag"
    row n-gone default "a case the table lacks"
} >"$TMP/nr/e.tsv"
expect_error "a row for a case the case table lacks is refused" "$TMP/nr/e.tsv:3: 'n-gone' is not a case of $T_NA" "$T_NA" --case-runs "$TMP/nr/e.tsv"
# A case cell that holds a blank names no case: the table's first case behind a blank, and two
# cases that stand side by side in the table written in one cell.
{
    row case run reason
    row " n-own" versioned "the first case behind a blank"
    row n-added default "the own run plus a flag"
} >"$TMP/nr/e.tsv"
expect_error "a case cell that is the first case's name behind a blank is refused as no case of the table" "$TMP/nr/e.tsv:2: ' n-own' is not a case of $T_NA" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row "n-own n-added" default "two cases in one cell"
} >"$TMP/nr/e.tsv"
expect_error "a case cell that holds two neighbouring cases is refused as no case of the table" "$TMP/nr/e.tsv:2: 'n-own n-added' is not a case of $T_NA" "$T_NA" --case-runs "$TMP/nr/e.tsv"
printf 'case\trun\n' >"$TMP/nr/e.tsv"
expect_error "a table whose header is not case, run, reason is refused" "$TMP/nr/e.tsv:1: the header is not case, run, reason" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row n-added default
} >"$TMP/nr/e.tsv"
expect_error "a row without a reason is refused" "$TMP/nr/e.tsv:2: a row is case, run, reason: three tab-separated cells, none empty" "$T_NA" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row n-added "two words" "a run's name has no blank"
} >"$TMP/nr/e.tsv"
expect_error "a run cell that is no name is refused" "$TMP/nr/e.tsv:2: run 'two words' is not a run's name" "$T_NA" --case-runs "$TMP/nr/e.tsv"
expect_error "a --case-runs file that is not there is refused" "named-runs table not found: $TMP/nr/absent.tsv" "$T_NA" --case-runs "$TMP/nr/absent.tsv"
# A real run stops at the same refusal before any engine is called.
{
    row case run reason
    row n-added no-such-run "a run the manifest lacks"
} >"$TMP/nr/e.tsv"
STUB_DOCKER_LOG="$TMP/nr-docker.log" STUB_CAND_LOG="$TMP/nr-cand.log" run_oracle "$T_NA" --case-runs "$TMP/nr/e.tsv" --out "$TMP/out-nr" --candidates "$STUBDIR/vep" --candidate-pins p1 >/dev/null 2>"$TMP/nr-real.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$TMP/nr-docker.log" && ! -e "$TMP/nr-cand.log" && ! -e "$TMP/out-nr/status.tsv" ]] && grep -q "names the run 'no-such-run'" "$TMP/nr-real.err" &&
    ok "a row naming an unknown run stops a real run before any engine call" || bad "unknown run on a real run: rc=$rc"

# A case and a run are found by their names as text, and the table is read whole: awk reads
# the numerals 7 and 7.0, or 1 and 1.0, as one number.
mkdir -p "$NR/numerals" "$NR/unprintable"
cp "$RS/corpus/variants.vcf" "$NR/numerals/variants.vcf"
cat >"$NR/numerals/manifest.json" <<'JSON'
{"flags": ["--symbol"],
 "flag_sets": [
  {"name": "1", "flags": ["--symbol", "--pick"]},
  {"name": "1.0", "input": "second.vcf", "flags": ["--symbol", "--pick"]}]}
JSON
T_NN="$TMP/nr/numerals.tsv"
{
    echo "$HEADER"
    nr_row 7 namedruns/numerals/variants.vcf "--symbol --canonical"
    nr_row 7.0 namedruns/numerals/variants.vcf "--symbol --canonical"
} >"$T_NN"
{
    row case run reason
    row 7.0 pooled "no run of the manifest carries --canonical"
    row 7 1 "the row of the case 7"
} >"$TMP/nr/numerals-runs.tsv"
dry "$T_NN" "$TMP/nr-numerals.out" --case-runs "$TMP/nr/numerals-runs.tsv"
rc=$?
[[ "$rc" -eq 0 ]] && grep -qxF "DRY-RUN 7 run 1, named by $TMP/nr/numerals-runs.tsv: the row of the case 7" "$TMP/nr-numerals.out" && grep -qxF "DRY-RUN 7.0 run pooled, no run named by $TMP/nr/numerals-runs.tsv: no run of the manifest carries --canonical" "$TMP/nr-numerals.out" &&
    ok "a case's row is found by its name as text: the cases 7 and 7.0 each take their own row" || bad "numeral case names: rc=$rc $(grep '^DRY-RUN 7[.0]* run ' "$TMP/nr-numerals.out" | tr '\n' ' ')"
{
    row case run reason
    row 7 1.0 "a run over another input"
    row 7.0 pooled "no run of the manifest carries --canonical"
} >"$TMP/nr/e.tsv"
expect_error "a run is found by its name as text: a row naming the run 1.0 over another input is refused beside the run 1 over the case's" "7: $TMP/nr/e.tsv names the run '1.0', which annotates second.vcf, and the case's input is variants.vcf" "$T_NN" --case-runs "$TMP/nr/e.tsv"
{
    row case run reason
    row 7 1 "the first row of a long table"
    row 7.0 pooled "$(printf '%0300000d' 0)"
} >"$TMP/nr/long-runs.tsv"
dry "$T_NN" "$TMP/nr-long.out" --case-runs "$TMP/nr/long-runs.tsv"
rc=$?
[[ "$rc" -eq 0 ]] && grep -qxF "DRY-RUN 7 run 1, named by $TMP/nr/long-runs.tsv: the first row of a long table" "$TMP/nr-long.out" &&
    ok "a table longer than a pipe holds is read whole: the case whose row comes first still plans" || bad "long named-runs table: rc=$rc $(head -2 "$TMP/nr-long.out.err" | cut -c1-200 | tr '\n' ' ')"
# A run match that fails says nothing of the runs over the input: a run's name that standard
# output cannot carry ends the match, and the case is refused where it would compare unnamed.
cp "$RS/corpus/variants.vcf" "$NR/unprintable/variants.vcf"
printf '%s\n' '{"flags": ["--symbol"], "flag_sets": [{"name": "\ud800", "flags": ["--symbol", "--pick"]}]}' >"$NR/unprintable/manifest.json"
{
    echo "$HEADER"
    nr_row n-unprintable namedruns/unprintable/variants.vcf "--symbol --canonical"
} >"$TMP/nr/unprintable.tsv"
expect_error "a case whose run match fails is refused by name, not compared unnamed" "n-unprintable: the run match over $NR/unprintable/manifest.json failed" "$TMP/nr/unprintable.tsv"

# A manifest-only class: named where no manifest stands beside the input, it is refused as
# such, from the runner's list.
T_MO="$TMP/nr/manifest-only.tsv"
{
    echo "$HEADER"
    nr_row n-only nomanifest/variants.vcf "--symbol" transcript_annotator_insertion_window
} >"$T_MO"
expect_error "a case naming a manifest-only class with no manifest beside its input fails preflight as such" "n-only: masks names 'transcript_annotator_insertion_window', a manifest-only class (the comparator's MANIFEST_CLASSES: no rule over the two outputs decides it), and no manifest.json stands beside the input $DATA/nomanifest/variants.vcf" "$T_MO"
{
    echo "$HEADER"
    nr_row n-only nomanifest/variants.vcf "--symbol" class_of_no_list
} >"$T_MO"
expect_error "a class of neither list keeps the manifest refusal" "n-only: masks names 'class_of_no_list', which applies from a manifest.json beside the input, and none is there" "$T_MO"
grep -q 'manifest-only' "$TMP/e.out.err" && bad "a class outside the runner's list is called manifest-only" || ok "a class outside the runner's list is not called manifest-only"
# The runner's list is held to the comparator's class table: the manifest-only classes are
# the table's MANIFEST_CLASSES less the built-in classes. The column is read whole, a class
# whose CORRECTED_IN cell names an Ensembl VEP release with the others: the comparator keeps
# such a class among its manifest classes and applies it from a manifest's list under an
# earlier release, so the runner still calls it manifest-only where no manifest stands
# beside the input.
builtin_line=$(awk -F'"' '/^BUILTIN_MASKS="/ {print $2; exit}' "$ORACLE")
only_listed=$(awk -F'"' '/^MANIFEST_ONLY_MASKS="/ {print $2; exit}' "$ORACLE" | tr ' ' '\n' | sort)
only_table=$(awk -F'\t' '
    NR == 1 { for (i = 2; i <= NF; i++) if ($i == "MANIFEST_CLASSES") column = i; next }
    column && $column ~ /^[0-9]+$/ { print $1 }
' "$TOOL/data/classes.tsv" | while read -r table_class; do
    [[ " $builtin_line " == *" $table_class "* ]] || echo "$table_class"
done | sort)
[[ -n "$only_listed" && "$only_listed" == "$only_table" ]] &&
    ok "the runner's manifest-only classes are the comparator's MANIFEST_CLASSES less the built-in classes" || bad "MANIFEST_ONLY_MASKS [$(tr '\n' ' ' <<<"$only_listed")] against the class table [$(tr '\n' ' ' <<<"$only_table")]"
# The list of the classes a discovery pass asks for by name alone follows the class table
# too: regulatory_feature_rows and the classes under which the comparator compares the whole
# outputs, its column FRONT_END_WHOLE_CLASSES.
opt_in_listed=$(awk -F'"' '/^DISCOVERY_OPT_IN_MASKS="/ {print $2; exit}' "$ORACLE" | tr ' ' '\n' | sort)
opt_in_table=$({
    echo regulatory_feature_rows
    class_table_members FRONT_END_WHOLE_CLASSES
} | sort)
[[ -n "$opt_in_listed" && "$opt_in_listed" == "$opt_in_table" ]] &&
    ok "a discovery pass asks by name alone for regulatory_feature_rows and the classes the comparator reads whole outputs under" || bad "DISCOVERY_OPT_IN_MASKS [$(tr '\n' ' ' <<<"$opt_in_listed")] against the class table [$(tr '\n' ' ' <<<"$opt_in_table")]"

echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
