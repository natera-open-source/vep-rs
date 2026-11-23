#!/usr/bin/env bash
# Tests for check_genome.sh: each mode's command lines, each exit code, each exclusion,
# the stored mode's refusal of records that differ, --stop-after at each stage, the
# option both engines are given on a cache edited from an alignment file, and that
# nothing is written outside --work and --out.
#
# Bash, tar, gzip and python3 with made-up caches. No Docker, no Ensembl VEP, no vep-rs
# build: a stub `docker` stands first on PATH, a stub build of the converter stands under
# a made-up CARGO_TARGET_DIR, and stub engines and tools are passed by path. Each stub
# appends its arguments to a log the checks read.
#
# Run: bash scripts/genomes/test_check_genome.sh
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
JOB="$HERE/check_genome.sh"
TABLE="$HERE/../../tests/genomes"

for tool in tar gzip python3; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: [test_check_genome] $tool is not on PATH" >&2
        exit 2
    }
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
BIN="$TMP/bin"
TOOLS="$TMP/tools"
CTL="$TMP/ctl"
mkdir -p "$BIN" "$TOOLS" "$CTL"
export CTL

# The info.json a stub conversion writes, by the converter's rule for the transcript
# set: the suffix of the directory above <cache version>[_<assembly>] in the path the
# converter is given, and no set for a path outside that layout. The second argument is
# the directory that holds the cache's info.txt, whose `bam` line, the alignment file
# whose edits the cache's transcripts carry, the conversion carries over.
INFO_JSON_STUB="$TMP/info_json.stub"
export INFO_JSON_STUB
cat >"$INFO_JSON_STUB" <<'STUB'
#!/usr/bin/env bash
layout='(^|/)([^/]+)/[0-9]+(_[^/]+)?/?$'
if [[ ! $1 =~ $layout ]]; then
    echo '{}'
    exit 0
fi
bam=""
if [ -f "${2:-}/info.txt" ]; then
    bam=$(awk -F '\t' '$1 == "bam" && $2 != "-" { value = $2 } END { print value }' "$2/info.txt")
fi
[ -z "$bam" ] || bam="\"bam\": \"$bam\", "
case "${BASH_REMATCH[2]}" in
*_refseq) echo "{$bam\"source_type\": \"refseq\"}" ;;
*_merged) echo "{$bam\"source_type\": \"merged\"}" ;;
*) echo "{$bam\"source_type\": \"ensembl\"}" ;;
esac
STUB

# The stub docker: `image inspect` answers for the image unless ctl/no_image stands;
# `run` appends its arguments to ctl/docker.log and plays the converter or Ensembl VEP
# by what the arguments name, writing into the host directory mounted at /output. The
# converter declares the transcript set of the path it is given. Ensembl VEP ends with
# its own line and exit 2 on a cache whose info.txt names an alignment file where the
# command holds no --use_given_ref (CacheDir.pm:412-414 and Runner.pm:762 at 116.3).
cat >"$BIN/docker" <<'STUB'
#!/usr/bin/env bash
mode() { if [ -f "$CTL/$1" ]; then cat "$CTL/$1"; else echo ok; fi; }
case "$1" in
image)
    [ ! -f "$CTL/no_image" ] || exit 1
    case "$*" in *'{{.Id}}'*) echo "sha256:0000000000000000000000000000000000000000000000000000000000000001" ;; esac
    exit 0
    ;;
inspect) echo 0; exit 0 ;;
run) ;;
*) exit 64 ;;
esac
echo "$*" >>"$CTL/docker.log"
output=""
cache_path=""
cache_host=""
cache_root=""
previous=""
for word in "$@"; do
    if [ "$previous" = -v ]; then
        case "$word" in
        *:/output) output=${word%:/output} ;;
        *:/input/*:ro) cache_host=${word%%:/input/*} ;;
        *:/cache:ro) cache_root=${word%:/cache:ro} ;;
        esac
    fi
    case "$previous" in */vep-storable-to-json) cache_path=$word ;; esac
    previous=$word
done
case "$*" in
*vep-storable-to-json*)
    case "$(mode convert_mode)" in
    ok)
        mkdir -p "$output/transcripts/1"
        bash "$INFO_JSON_STUB" "$cache_path" "$cache_host" >"$output/info.json"
        printf 'name\tlength\tcircular\n1\t1000\t0\n' >"$output/sequence_regions.tsv"
        echo "Done. Converted: 1 / 1 files, 5 total transcripts, 0 variation files copied, 0 regulatory files"
        ;;
    zero)
        echo '{}' >"$output/info.json"
        echo "Done. Converted: 1 / 1 files, 0 total transcripts, 0 variation files copied, 0 regulatory files"
        ;;
    sereal) echo "ERROR: [storable_to_json] /input/x uses the sereal serialiser; only a Storable cache converts" >&2; exit 1 ;;
    noshards) echo "ERROR: [storable_to_json] no transcript shards (<start>-<end>.gz) under /input/x" >&2; exit 1 ;;
    shards) echo "ERROR: [storable_to_json] 2 of 3 shards failed to convert (see the WARN lines above)" >&2; exit 1 ;;
    crash) echo "docker: the daemon stopped" >&2; exit 125 ;;
    esac
    ;;
*' vep '*)
    edited=$(awk -F '\t' '$1 == "bam" && $2 != "" && $2 != "0" && $2 != "-" { print "yes"; exit }' \
        "$cache_root"/*/*/info.txt 2>/dev/null)
    if [ -n "$edited" ]; then
        case " $* " in
        *' --use_given_ref '*) ;;
        *)
            echo "MSG: ERROR: Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref" >&2
            exit 2
            ;;
        esac
    fi
    case "$(mode ensembl_mode)" in
    ok) printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n%s\n%s\n' "$(printf 'r1\t1:100\tG')" "$(printf 'r2\t1:200\tT')" >"$output/variants.txt" ;;
    warn)
        printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n%s\n' "$(printf 'r1\t1:100\tG')" >"$output/variants.txt"
        echo "WARNING: a made-up warning" >"$output/variants.warnings"
        ;;
    norows) printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n' >"$output/variants.txt" ;;
    fail) echo "a made-up failure" >&2; exit 2 ;;
    esac
    ;;
esac
exit 0
STUB

# The stub build of the converter, for a conversion without an image: it appends its
# arguments to the log STUB_CONVERTER_LOG names, ctl/converter.log by default, and plays
# the converter into its second argument, declaring the transcript set of the path in
# its first. The copy under the made-up CARGO_TARGET_DIR is the build every case finds;
# the stub docker plays the converter itself, so a conversion in the image mounts this
# file and never starts it.
cat >"$TMP/converter.stub" <<'STUB'
#!/usr/bin/env bash
echo "$*" >>"$CTL/${STUB_CONVERTER_LOG:-converter.log}"
bash "$INFO_JSON_STUB" "$1" "$1" >"$2/info.json"
echo "Done. Converted: 1 / 1 files, 5 total transcripts, 0 variation files copied, 0 regulatory files"
STUB
STUB_TARGET="$TMP/cargo-target"
STUB_BUILD="$STUB_TARGET/release/vep-storable-to-json"
mkdir -p "$STUB_TARGET/release"
cp "$TMP/converter.stub" "$STUB_BUILD"
export CARGO_TARGET_DIR="$STUB_TARGET"

cat >"$TOOLS/vep-seeded-variants" <<'STUB'
#!/usr/bin/env bash
if [ "$1" = --version ]; then echo "vep-seeded-variants 0.0.0 (generator version 1)"; exit 0; fi
echo "$*" >>"$CTL/generator.log"
output="" summary="" seed=""
while [ $# -gt 0 ]; do
    case "$1" in
    --output) output=$2; shift ;;
    --summary) summary=$2; shift ;;
    --seed) seed=$2; shift ;;
    esac
    shift
done
version=$(if [ -f "$CTL/generator_version" ]; then cat "$CTL/generator_version"; else echo 1; fi)
case "$(if [ -f "$CTL/generator_mode" ]; then cat "$CTL/generator_mode"; else echo ok; fi)" in
fail) echo "ERROR: [vep-seeded-variants] a made-up failure" >&2; exit 1 ;;
empty) printf '##fileformat=VCFv4.2\n##source=vep-seeded-variants generator_version=%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' "$version" >"$output" ;;
null_byte) printf '##fileformat=VCFv4.2\n##source=vep-seeded-variants generator_version=%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\tfirst_base_%s\0\tA\tG\t.\t.\t.\n1\t200\tintron_%s\tC\tT\t.\t.\t.\n' "$version" "$seed" "$seed" >"$output" ;;
*)
    {
        printf '##fileformat=VCFv4.2\n##source=vep-seeded-variants generator_version=%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' "$version"
        printf '1\t100\tfirst_base_%s\tA\tG\t.\t.\t.\n' "$seed"
        printf '1\t200\tintron_%s\tC\tT\t.\t.\t.\n' "$seed"
        printf '1\t300\tupstream_%s\tG\tA\t.\t.\t.\n' "$seed"
    } >"$output"
    ;;
esac
printf 'category\twanted\twritten\tpassed_over\treason\n' >"$summary"
STUB

cat >"$TOOLS/vep" <<'STUB'
#!/usr/bin/env bash
if [ "$1" = --version ]; then echo "vep 0.0.0"; exit 0; fi
echo "$*" >>"$CTL/vep_rs.log"
output="" json_cache="" selected=ensembl given_ref=false
while [ $# -gt 0 ]; do
    case "$1" in
    -o) output=$2 ;;
    --json_cache) json_cache=$2 ;;
    --refseq) selected=refseq ;;
    --merged) selected=merged ;;
    --use_given_ref) given_ref=true ;;
    esac
    shift
done
# vep-rs refuses a converted cache that declares another transcript set than the one
# its flags select.
declared=$(sed -n 's/.*"source_type": "\([a-z]*\)".*/\1/p' "$json_cache/info.json" 2>/dev/null)
if [ -n "$declared" ] && [ "$declared" != "$selected" ]; then
    echo "ERROR: [cache] the flags select the $selected transcript set but $json_cache/info.json declares source_type $declared" >&2
    exit 1
fi
# vep-rs refuses a converted cache that names an alignment file, the `bam` string of
# info.json's top-level object, where its flags hold neither --use_given_ref nor a
# FASTA file.
bam=$(python3 -c 'import json, sys
value = json.load(open(sys.argv[1])).get("bam")
print(value if isinstance(value, str) else "")' "$json_cache/info.json" 2>/dev/null)
if [ -n "$bam" ] && [ "$bam" != 0 ] && [ "$given_ref" = false ]; then
    echo "ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref" >&2
    exit 1
fi
case "$(if [ -f "$CTL/vep_rs_mode" ]; then cat "$CTL/vep_rs_mode"; else echo ok; fi)" in
fail) echo "ERROR: [vep] a made-up failure" >&2; exit 1 ;;
norows) printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n' >"$output" ;;
*) printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n%s\n%s\n' "$(printf 'r1\t1:100\tG')" "$(printf 'r2\t1:200\tT')" >"$output" ;;
esac
STUB

cat >"$TOOLS/vep-outputs-compare" <<'STUB'
#!/usr/bin/env bash
if [ "$1" = --version ]; then echo "usage: vep-outputs-compare [options]" >&2; exit 2; fi
echo "$*" >>"$CTL/comparator.log"
report=""
while [ $# -gt 0 ]; do
    if [ "$1" = --report-dir ]; then report=$2; fi
    shift
done
mode=$(if [ -f "$CTL/compare_mode" ]; then cat "$CTL/compare_mode"; else echo same; fi)
[ "$mode" != fail ] || { echo "ERROR: [vep-outputs-compare] a made-up failure" >&2; exit 1; }
mkdir -p "$report"
echo '{}' >"$report/summary.json"
echo '# summary' >"$report/summary.md"
printf 'file\tside\tkey\n' >"$report/discordant.tsv"
printf 'file\tside\tkey\n' >"$report/discordant_open.tsv"
if [ "$mode" = differ ]; then
    printf 'variants.txt\tmissing\tr1\nvariants.txt\textra\tr9\n' >>"$report/discordant_open.tsv"
fi
STUB

# The stub table tool: `digest` prints ctl/digest or a fixed digest; `row` prints the
# field ctl/row.<field> holds, and refuses the key while ctl/row.missing stands.
cat >"$TOOLS/vep-genome-table" <<'STUB'
#!/usr/bin/env bash
echo "$*" >>"$CTL/table.log"
case "$1" in
digest)
    if [ -f "$CTL/digest" ]; then cat "$CTL/digest"; else echo "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"; fi
    ;;
row)
    if [ -f "$CTL/row.missing" ]; then echo "ERROR: [vep-genome-table] no row has the key $3" >&2; exit 1; fi
    cat "$CTL/row.$5"
    ;;
*) echo "ERROR: [vep-genome-table] no subcommand $1" >&2; exit 2 ;;
esac
STUB
chmod +x "$BIN/docker" "$STUB_BUILD" "$TOOLS"/*

if [[ ! -f "$JOB" ]]; then
    bad "check_genome.sh not found beside this test"
    echo "passed: $PASS   failed: $FAIL"
    exit 1
fi
if bash -n "$JOB"; then ok "check_genome.sh parses"; else bad "check_genome.sh has a syntax error"; fi

# Made-up caches: an archive as Ensembl ships one, an archive of a genome kept in a
# collection, a RefSeq-set cache directory, a cache whose info.txt names sereal and a
# file that is no archive.
SRC="$TMP/src"
mkdir -p "$SRC/plain/made_up_genome/63_ASM1/1" "$SRC/collection/made_up_collection/made_up_genome/63_ASM1/1" \
    "$SRC/refseq/made_up_genome_refseq/116_ASM1/1" "$SRC/merged/made_up_genome_merged/63_ASM1/1" \
    "$SRC/sereal/made_up_genome/63_ASM1/1" "$SRC/other/another_genome/63_ASM1/1"
printf 'species\tmade_up_genome\nassembly\tASM1\n' >"$SRC/plain/made_up_genome/63_ASM1/info.txt"
printf 'species\tmade_up_genome\nserialiser_type\tsereal\n' >"$SRC/sereal/made_up_genome/63_ASM1/info.txt"
tar -czf "$TMP/plain.tar.gz" -C "$SRC/plain" made_up_genome
tar -czf "$TMP/collection.tar.gz" -C "$SRC/collection" made_up_collection
tar -czf "$TMP/other.tar.gz" -C "$SRC/other" another_genome
echo "not an archive" >"$TMP/broken.tar.gz"
mkdir -p "$TMP/ready_json"
echo '{}' >"$TMP/ready_json/info.json"

IDENTITY=(--genome made_up_genome --assembly ASM1 --transcripts ensembl --set fungi --cache-version 63)
BINARIES=(--vep-rs "$TOOLS/vep" --generator "$TOOLS/vep-seeded-variants"
    --comparator "$TOOLS/vep-outputs-compare" --table-tool "$TOOLS/vep-genome-table")
IMAGE=(--vep-image made-up/image:tag)
KEY="fungi__made_up_genome__ASM1__ensembl"

# run_case <name> <check_genome.sh arguments...>: a fresh --work and --out under
# $TMP/case/<name>; the exit code lands in CODE, the output in $CASE/stdout and stderr.
# The control files and logs are cleared, so each case sets the modes it wants first
# through `set_mode`, which run_case does not undo.
PENDING_MODES=()
set_mode() { PENDING_MODES+=("$1=$2"); }
run_case() {
    local name=$1 pair
    shift
    CASE="$TMP/case/$name"
    rm -rf "$CASE" "$CTL"
    mkdir -p "$CASE" "$CTL"
    for pair in ${PENDING_MODES[@]+"${PENDING_MODES[@]}"}; do
        if [ -n "${pair#*=}" ]; then printf '%s\n' "${pair#*=}" >"$CTL/${pair%%=*}"; else : >"$CTL/${pair%%=*}"; fi
    done
    PENDING_MODES=()
    PATH="$BIN:$PATH" bash "$JOB" "$@" --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr"
    CODE=$?
}
expect_code() { # expect_code <name> <code>
    if [[ $CODE -eq $2 ]]; then ok "$1: exit $2"; else bad "$1: exit $CODE, wanted $2 ($(tail -n 1 "$CASE/stderr"))"; fi
}
expect_file() {
    if [[ -e "$CASE/$2" ]]; then ok "$1: $2 written"; else bad "$1: $2 is missing"; fi
}
expect_no_file() {
    if [[ ! -e "$CASE/$2" ]]; then ok "$1: no $2"; else bad "$1: $2 was written"; fi
}
expect_row() { # expect_row <name> <whole row, tab-separated>
    local row
    row=$(sed -n '2p' "$CASE/out/row.tsv" 2>/dev/null)
    if [[ $row == "$2" ]]; then ok "$1: row $(printf '%s' "$2" | cut -f8,9 | tr '\t' ' ')"; else bad "$1: row is [$row]"; fi
}
expect_log() { # expect_log <name> <log file under ctl> <fixed text of one line>
    if grep -q -F -- "$3" "$CTL/$2" 2>/dev/null; then ok "$1"; else bad "$1 ($(cat "$CTL/$2" 2>/dev/null))"; fi
}
expect_no_log() {
    if ! grep -q -F -- "$3" "$CTL/$2" 2>/dev/null; then ok "$1"; else bad "$1 ($(grep -F -- "$3" "$CTL/$2"))"; fi
}
expect_stderr() {
    if grep -q -F -- "$2" "$CASE/stderr"; then ok "$1"; else bad "$1 ($(cat "$CASE/stderr"))"; fi
}
expect_declared() { # expect_declared <name> <set>: the set the converted cache's info.json declares
    if grep -q -F -- "\"source_type\": \"$2\"" "$CASE/work/json_cache/info.json" 2>/dev/null; then
        ok "$1: the converted cache declares the $2 set"
    else
        bad "$1: the converted cache's info.json is [$(cat "$CASE/work/json_cache/info.json" 2>/dev/null)]"
    fi
}
row_of() { # row_of <records> <outcome> <sha256> [seed]
    printf 'made_up_genome\tASM1\tensembl\tfungi\t63\t%s\t%s\t%s\t%s' "${4:-1}" "$1" "$2" "$3"
}
DIGEST="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"

# 1. Parity on an archive: every step's command line, the files left, exit 0.
run_case parity "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code parity 0
expect_row parity "$(row_of 3 exact "$DIGEST")"
if [[ $(head -n 1 "$CASE/out/row.tsv") == "$(head -n 1 "$TABLE/genomes.tsv")" ]]; then
    ok "parity: row.tsv opens with the table's header line"
else
    bad "parity: row.tsv header is [$(head -n 1 "$CASE/out/row.tsv")]"
fi
expect_log "parity: the conversion runs in the image as the caller, the converter's build mounted read-only" docker.log \
    "--user $(id -u):$(id -g) -v $CASE/work/ensembl_cache/made_up_genome/63_ASM1:/input/made_up_genome/63_ASM1:ro -v $CASE/work/json_cache:/output -v $STUB_BUILD:/tools/vep-storable-to-json:ro made-up/image:tag /tools/vep-storable-to-json /input/made_up_genome/63_ASM1 /output"
expect_declared parity ensembl
if [[ ! -e "$CTL/converter.log" ]]; then
    ok "parity: the build is not started on the host when the conversion runs in the image"
else
    bad "parity: the build ran on the host ($(cat "$CTL/converter.log"))"
fi
expect_log "parity: Ensembl VEP runs offline on the cache under default flags" docker.log \
    "-v $CASE/work/ensembl_cache:/cache:ro -v $CASE/out:/input:ro -v $CASE/work/ensembl_vep:/output made-up/image:tag vep --offline --cache --dir_cache /cache --species made_up_genome --assembly ASM1 --cache_version 63 -i /input/variants.vcf -o /output/variants.txt --warning_file /output/variants.warnings --no_stats --force_overwrite"
expect_no_log "parity: Ensembl VEP is given no --fork" docker.log "--fork"
expect_log "parity: docker never pulls" docker.log "run --rm --pull never"
expect_log "parity: the generator reads the converted cache under seed 1" generator.log \
    "--cache $CASE/work/json_cache --seed 1 --output $CASE/out/variants.vcf --summary $CASE/out/variants.summary.tsv"
expect_log "parity: vep-rs runs offline on the converted cache under default flags" vep_rs.log \
    "-i $CASE/out/variants.vcf -o $CASE/work/vep_rs/variants.txt --offline --json_cache $CASE/work/json_cache --species made_up_genome --assembly ASM1 --warning_file $CASE/work/logs/vep_rs.warnings --no_stats"
expect_log "parity: the comparator reads both outputs with the input records, no contig filter and the Ensembl cache" comparator.log \
    "--perl-dir $CASE/work/ensembl_vep --rust-dir $CASE/work/vep_rs --report-dir $CASE/out/comparison --input-vcf $CASE/out/variants.vcf --no-canonical-contigs --assembly ASM1 --ensembl-cache $CASE/work/ensembl_cache/made_up_genome/63_ASM1"
expect_no_log "parity: the comparator is given no converted cache for a cache of the Ensembl transcript set" comparator.log "--json-cache"
expect_log "parity: the digest is the table tool's over vep-rs's output" table.log "digest $CASE/work/vep_rs/variants.txt"
for left in out/variants.vcf out/variants.summary.tsv out/ensembl_vep.txt.gz out/ensembl_vep.warnings.txt \
    out/comparison/summary.json out/job.json work/layout.tsv work/steps.tsv; do
    expect_file parity "$left"
done
if [[ $(gzip -dc "$CASE/out/ensembl_vep.txt.gz" | grep -c -v '^#') -eq 2 ]]; then
    ok "parity: ensembl_vep.txt.gz holds Ensembl VEP's rows"
else
    bad "parity: ensembl_vep.txt.gz does not hold Ensembl VEP's rows"
fi
if [[ ! -s "$CASE/out/ensembl_vep.warnings.txt" ]]; then
    ok "parity: the warnings file is empty where Ensembl VEP wrote none"
else
    bad "parity: the warnings file is not empty"
fi
if python3 - "$CASE/out/job.json" "$KEY" <<'PY'; then
import json, sys
job = json.load(open(sys.argv[1]))
assert job["genome_key"] == sys.argv[2] and job["mode"] == "parity" and job["seed"] == 1, job
assert job["outcome"] == "exact" and job["exit_code"] == 0 and job["records"] == 3, job
assert len(job["records_sha256"]) == 64 and job["generator_version"] == "1", job
assert job["vep_image"] == "made-up/image:tag" and job["vep_image_id"].startswith("sha256:"), job
assert job["ensembl_vep_command"].startswith("vep --offline --cache --dir_cache /cache"), job
assert job["tarball_bytes"] > 0 and job["converted_bytes"] > 0, job
binaries = job["binaries"]
assert sorted(binaries) == ["comparator", "generator", "table_tool", "vep_rs"], binaries
assert all(len(binary["sha256"]) == 64 for binary in binaries.values()), binaries
assert binaries["vep_rs"]["version"] == "vep 0.0.0" and binaries["comparator"]["version"] is None, binaries
steps = [step["step"] for step in job["steps"]]
assert steps == ["extract", "convert", "generate", "ensembl_vep", "vep_rs", "digest", "compare"], steps
assert all(step["seconds"] >= 0 for step in job["steps"]), job["steps"]
assert all(step["peak_rss_kb"] is None or step["peak_rss_kb"] > 0 for step in job["steps"]), job["steps"]
PY
    ok "parity: job.json holds the steps, the image, the binaries and the records' digest"
else
    bad "parity: job.json ($(cat "$CASE/out/job.json" 2>/dev/null | head -c 400))"
fi
layout=$(cat "$CASE/work/layout.tsv")
wanted=$(printf 'ensembl_cache_root\t%s\nensembl_cache_dir\t%s\njson_cache\t%s\nsequence_regions\t%s' \
    "$CASE/work/ensembl_cache" "$CASE/work/ensembl_cache/made_up_genome/63_ASM1" "$CASE/work/json_cache" \
    "$CASE/work/json_cache/sequence_regions.tsv")
if [[ $layout == "$wanted" ]]; then ok "parity: layout.tsv names the four paths"; else bad "parity: layout.tsv is [$layout]"; fi
if [[ $(tail -n 1 "$CASE/stdout") == "check_genome: $KEY exact" ]]; then
    ok "parity: the last line names the key and the outcome"
else
    bad "parity: the last line is [$(tail -n 1 "$CASE/stdout")]"
fi

# 2. The set flag and the cache version come from the row, and the converter is given
# the cache under its species directory, whose suffix is the set the conversion declares.
run_case merged --genome made_up_genome --assembly ASM1 --transcripts merged --set plants --cache-version 63 \
    --mode parity --cache-dir "$SRC/merged" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code merged 0
expect_log "merged: the Ensembl VEP line carries --cache_version 63 and --merged" docker.log \
    "vep --offline --cache --dir_cache /cache --species made_up_genome --assembly ASM1 --cache_version 63 --merged -i /input/variants.vcf"
expect_log "merged: the converter is given the cache under its _merged species directory" docker.log \
    "-v $SRC/merged/made_up_genome_merged/63_ASM1:/input/made_up_genome_merged/63_ASM1:ro -v $CASE/work/json_cache:/output -v $STUB_BUILD:/tools/vep-storable-to-json:ro made-up/image:tag /tools/vep-storable-to-json /input/made_up_genome_merged/63_ASM1 /output"
expect_declared merged merged
expect_log "merged: vep-rs is given --merged" vep_rs.log "--species made_up_genome --assembly ASM1 --merged --warning_file"
expect_log "merged: the comparator is given the converted cache, which the classes of a shared dbID and of a shared stable id are told from" comparator.log \
    "--input-vcf $CASE/out/variants.vcf --no-canonical-contigs --assembly ASM1 --json-cache $CASE/work/json_cache"
expect_log "merged: the comparator is given the Ensembl cache after the converted cache" comparator.log \
    "--assembly ASM1 --json-cache $CASE/work/json_cache --ensembl-cache $SRC/merged/made_up_genome_merged/63_ASM1"
expect_row merged "$(printf 'made_up_genome\tASM1\tmerged\tplants\t63\t1\t3\texact\t%s' "$DIGEST")"
run_case refseq --genome made_up_genome --assembly ASM1 --transcripts refseq --set vertebrates --cache-version 116 \
    --mode parity --cache-dir "$SRC/refseq" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code refseq 0
expect_log "refseq: the Ensembl VEP line carries --cache_version 116 and --refseq" docker.log \
    "--assembly ASM1 --cache_version 116 --refseq -i /input/variants.vcf"
expect_log "refseq: the converter is given the cache under its _refseq species directory" docker.log \
    "-v $SRC/refseq/made_up_genome_refseq/116_ASM1:/input/made_up_genome_refseq/116_ASM1:ro -v $CASE/work/json_cache:/output -v $STUB_BUILD:/tools/vep-storable-to-json:ro made-up/image:tag /tools/vep-storable-to-json /input/made_up_genome_refseq/116_ASM1 /output"
expect_declared refseq refseq
expect_log "refseq: vep-rs is given --refseq" vep_rs.log "--assembly ASM1 --refseq --warning_file"
expect_log "refseq: the comparator is given the converted cache, which the class of a shared stable id is told from" comparator.log \
    "--input-vcf $CASE/out/variants.vcf --no-canonical-contigs --assembly ASM1 --json-cache $CASE/work/json_cache"
expect_log "refseq: the comparator is given the Ensembl cache after the converted cache" comparator.log \
    "--assembly ASM1 --json-cache $CASE/work/json_cache --ensembl-cache $SRC/refseq/made_up_genome_refseq/116_ASM1"
expect_no_log "refseq: a cache that names no alignment file leaves Ensembl VEP without --use_given_ref" docker.log "--use_given_ref"
expect_no_log "refseq: a cache that names no alignment file leaves vep-rs without --use_given_ref" vep_rs.log "--use_given_ref"
expect_row refseq "$(printf 'made_up_genome\tASM1\trefseq\tvertebrates\t116\t1\t3\texact\t%s' "$DIGEST")"
# A RefSeq-set archive, as most caches arrive.
tar -czf "$TMP/refseq.tar.gz" -C "$SRC/refseq" made_up_genome_refseq
run_case refseq-archive --genome made_up_genome --assembly ASM1 --transcripts refseq --set vertebrates --cache-version 116 \
    --mode parity --cache-tarball "$TMP/refseq.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code refseq-archive 0
expect_declared refseq-archive refseq
expect_row refseq-archive "$(printf 'made_up_genome\tASM1\trefseq\tvertebrates\t116\t1\t3\texact\t%s' "$DIGEST")"
# Without an image the build on the host is given the same directory by its own path.
run_case refseq-host --genome made_up_genome --assembly ASM1 --transcripts refseq --set vertebrates \
    --cache-version 116 --cache-dir "$SRC/refseq" --stop-after convert
expect_code refseq-host 0
expect_log "refseq-host: the build on the host is given the cache under its _refseq species directory" converter.log \
    "$SRC/refseq/made_up_genome_refseq/116_ASM1 $CASE/work/json_cache"
expect_declared refseq-host refseq
# A converted cache that declares another set than the row's: vep-rs refuses it, and the
# row says the outputs differ.
mkdir -p "$TMP/ensembl_json"
echo '{"source_type": "ensembl"}' >"$TMP/ensembl_json/info.json"
run_case refseq-other-set --genome made_up_genome --assembly ASM1 --transcripts refseq --set vertebrates --cache-version 116 \
    --mode parity --cache-dir "$SRC/refseq" --json-cache "$TMP/ensembl_json" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code refseq-other-set 1
expect_row refseq-other-set "$(printf 'made_up_genome\tASM1\trefseq\tvertebrates\t116\t1\t3\tdiffers:unclassified\t-')"
expect_stderr "refseq-other-set: the error says vep-rs wrote no row" "vep-rs wrote no row (exit 1;"

# 3. A difference: exit 1, the row says so and carries no digest.
set_mode compare_mode differ
run_case differs "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code differs 1
expect_row differs "$(row_of 3 differs:unclassified -)"
expect_file differs out/comparison/discordant_open.tsv
set_mode vep_rs_mode fail
run_case vep-rs-fails "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code vep-rs-fails 1
expect_row vep-rs-fails "$(row_of 3 differs:unclassified -)"
expect_file vep-rs-fails out/ensembl_vep.txt.gz
set_mode vep_rs_mode norows
run_case vep-rs-no-rows "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code vep-rs-no-rows 1
expect_row vep-rs-no-rows "$(row_of 3 differs:unclassified -)"
set_mode ensembl_mode warn
run_case warnings "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
if grep -q -F 'a made-up warning' "$CASE/out/ensembl_vep.warnings.txt" 2>/dev/null; then
    ok "warnings: Ensembl VEP's warnings are kept"
else
    bad "warnings: Ensembl VEP's warnings are not in the output directory"
fi

# 4. The four exclusions: exit 3 and a row that names the reason.
run_case unreadable-archive "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/broken.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code unreadable-archive 3
expect_row unreadable-archive "$(row_of 0 excluded:cache-unreadable -)"
run_case other-genome "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/other.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code other-genome 3
expect_row other-genome "$(row_of 0 excluded:cache-unreadable -)"
expect_stderr "other-genome: the error names the directory the cache lacks" "holds no directory made_up_genome/63_ASM1"
set_mode convert_mode shards
run_case failed-shards "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code failed-shards 3
expect_row failed-shards "$(row_of 0 excluded:cache-unreadable -)"
run_case sereal-info "${IDENTITY[@]}" --mode parity --cache-dir "$SRC/sereal" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code sereal-info 3
expect_row sereal-info "$(row_of 0 excluded:serialiser-not-storable -)"
expect_no_log "sereal-info: the conversion is not started" docker.log "vep-storable-to-json"
set_mode convert_mode sereal
run_case sereal-converter "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code sereal-converter 3
expect_row sereal-converter "$(row_of 0 excluded:serialiser-not-storable -)"
set_mode convert_mode noshards
run_case no-shards "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code no-shards 3
expect_row no-shards "$(row_of 0 excluded:no-transcripts -)"
set_mode convert_mode zero
run_case no-transcripts "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code no-transcripts 3
expect_row no-transcripts "$(row_of 0 excluded:no-transcripts -)"
set_mode ensembl_mode fail
run_case ensembl-fails "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code ensembl-fails 3
expect_row ensembl-fails "$(row_of 3 excluded:ensembl-vep-failed -)"
expect_no_file ensembl-fails out/ensembl_vep.txt.gz
expect_no_log "ensembl-fails: vep-rs is not run" vep_rs.log "--offline"
set_mode ensembl_mode norows
run_case ensembl-no-rows "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code ensembl-no-rows 3
expect_row ensembl-no-rows "$(row_of 3 excluded:ensembl-vep-failed -)"
declared=true
exclusions=0
while IFS= read -r name; do
    exclusions=$((exclusions + 1))
    grep -q "^excluded:$name"$'\t' "$TABLE/outcomes.tsv" || declared=false
done < <(sed -n 's/.*exclude \([a-z][a-z-]*\) .*/\1/p' "$JOB" | sort -u)
if [[ $declared == true && $exclusions -eq 4 ]]; then
    ok "the four exclusions the script writes are declared in tests/genomes/outcomes.tsv"
else
    bad "an exclusion the script writes is not declared in tests/genomes/outcomes.tsv"
fi

# 5. A failed step that is no exclusion: exit 3 and no row.
set_mode convert_mode crash
run_case converter-crash "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code converter-crash 3
expect_no_file converter-crash out/row.tsv
expect_file converter-crash out/job.json
set_mode generator_mode fail
run_case generator-fails "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code generator-fails 3
expect_no_file generator-fails out/row.tsv
set_mode generator_mode empty
run_case no-records "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code no-records 3
expect_no_file no-records out/row.tsv
set_mode compare_mode fail
run_case comparator-fails "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code comparator-fails 3
expect_no_file comparator-fails out/row.tsv

# 6. The skips and the usage errors: exit 2, one line, nothing run.
set_mode no_image ""
run_case no-image "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code no-image 2
expect_stderr "no-image: the skip names the image" "skip: no image made-up/image:tag in docker"
expect_no_file no-image out/row.tsv
CASE="$TMP/case/no-docker"
rm -rf "$CASE" "$CTL"
mkdir -p "$CASE" "$CTL" "$TMP/empty-path"
for tool in bash dirname cat mkdir rm awk; do ln -sf "$(command -v "$tool")" "$TMP/empty-path/$tool"; done
PATH="$TMP/empty-path" bash "$JOB" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}" \
    --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr"
CODE=$?
expect_code no-docker 2
expect_stderr "no-docker: the skip says docker is missing" "skip: no docker on PATH"
# A copy of the script under a made-up repository root, which holds no build of the
# converter at first. run_in_tree <name> <CARGO_TARGET_DIR or empty> <arguments...> runs
# the copy with that variable, or without it, in a fresh case directory.
TREE="$TMP/tree"
mkdir -p "$TREE/scripts/genomes" "$TMP/no-build-target/release"
cp "$JOB" "$TREE/scripts/genomes/check_genome.sh"
run_in_tree() {
    local name=$1 target=$2
    shift 2
    CASE="$TMP/case/$name"
    rm -rf "$CASE" "$CTL"
    mkdir -p "$CASE" "$CTL"
    if [ -n "$target" ]; then
        PATH="$BIN:$PATH" CARGO_TARGET_DIR="$target" bash "$TREE/scripts/genomes/check_genome.sh" "$@" \
            --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr"
    else
        PATH="$BIN:$PATH" env -u CARGO_TARGET_DIR bash "$TREE/scripts/genomes/check_genome.sh" "$@" \
            --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr"
    fi
    CODE=$?
}
NO_BUILD="ERROR: [check_genome] skip: no build of the converter at $TREE/tools/vep-storable-to-json/target/release/vep-storable-to-json: build it with 'cargo build --release --locked --manifest-path $TREE/tools/vep-storable-to-json/Cargo.toml'"
# expect_no_build <name>: exit 2 on the one line, nothing written and nothing started.
expect_no_build() {
    if [[ $CODE -eq 2 && $(cat "$CASE/stderr") == "$NO_BUILD" && ! -e "$CASE/work" && ! -e "$CASE/out" && ! -e "$CTL/docker.log" ]]; then
        ok "$1: exit 2 on the one line that names the build command, before anything is written or started"
    else
        bad "$1: exit $CODE ($(cat "$CASE/stderr"))"
    fi
}
run_in_tree no-build "" "${IDENTITY[@]}" --mode stored --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}" --reference-dir "$TMP"
expect_stderr "no-build: a missing reference file is a usage error before any skip" "no ensembl_vep.txt.gz under"
mkdir -p "$TMP/reference"
printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n%s\n' "$(printf 'r1\t1:100\tG')" | gzip -n -c >"$TMP/reference/ensembl_vep.txt.gz"
printf '##fileformat=VCFv4.2\n##source=vep-seeded-variants generator_version=1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\tfirst_base_1\tA\tG\t.\t.\t.\n1\t200\tintron_1\tC\tT\t.\t.\t.\n1\t300\tupstream_1\tG\tA\t.\t.\t.\n' |
    gzip -n -c >"$TMP/reference/variants.vcf.gz"
run_in_tree no-build "" "${IDENTITY[@]}" --mode stored --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}" --reference-dir "$TMP/reference"
expect_code no-build 2
expect_stderr "no-build: the skip says the converter has no build" "skip: no build of the converter at $TREE/tools/vep-storable-to-json/target/release/vep-storable-to-json"
expect_no_build "no-build: a conversion on the host"
run_in_tree no-build-image "" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_no_build "no-build: a conversion in the image, which supplies no converter"
run_in_tree no-build-target "$TMP/no-build-target" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_no_build "no-build: CARGO_TARGET_DIR names a directory without a build, and the line names the tool's own target directory"
: >"$TMP/no-build-target/release/vep-storable-to-json"
run_in_tree no-build-not-executable "$TMP/no-build-target" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_no_build "no-build: the file under CARGO_TARGET_DIR is not executable"
run_in_tree no-build-converted "" "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code no-build-converted 0
run_in_tree no-build-extract "" "${IDENTITY[@]}" --stop-after extract --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}"
expect_code no-build-extract 0
if [[ ! -s "$CASE/stderr" ]]; then
    ok "no-build: a run that converts nothing asks for no build"
else
    bad "no-build: a run that converts nothing printed [$(cat "$CASE/stderr")]"
fi

# The build is found under CARGO_TARGET_DIR first, then under the tool's own target
# directory: the tree is given a build of its own, which logs to ctl/own.log.
mkdir -p "$TREE/tools/vep-storable-to-json/target/release"
OWN_BUILD="$TREE/tools/vep-storable-to-json/target/release/vep-storable-to-json"
: >"$OWN_BUILD"
run_in_tree no-build-own-not-executable "" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_no_build "no-build: the file under the tool's own target directory is not executable"
printf '#!/usr/bin/env bash\nSTUB_CONVERTER_LOG=own.log exec "%s" "$@"\n' "$TMP/converter.stub" >"$OWN_BUILD"
chmod +x "$OWN_BUILD" "$TMP/converter.stub"
run_in_tree own-build "" "${IDENTITY[@]}" --mode stored --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}" --reference-dir "$TMP/reference"
expect_code own-build 0
expect_log "own-build: without CARGO_TARGET_DIR the build under the tool's own target directory converts" own.log \
    "$CASE/work/ensembl_cache/made_up_genome/63_ASM1 $CASE/work/json_cache"
run_in_tree own-build-image "" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code own-build-image 0
expect_log "own-build-image: the image mounts the build under the tool's own target directory" docker.log \
    "-v $OWN_BUILD:/tools/vep-storable-to-json:ro made-up/image:tag /tools/vep-storable-to-json /input/made_up_genome/63_ASM1 /output"
run_in_tree both-builds "$STUB_TARGET" "${IDENTITY[@]}" --mode stored --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}" --reference-dir "$TMP/reference"
expect_code both-builds 0
if grep -q -F -- "$CASE/work/json_cache" "$CTL/converter.log" 2>/dev/null && [[ ! -e "$CTL/own.log" ]]; then
    ok "both-builds: the build under CARGO_TARGET_DIR converts, and the tool's own is not started"
else
    bad "both-builds: converter.log [$(cat "$CTL/converter.log" 2>/dev/null)] own.log [$(cat "$CTL/own.log" 2>/dev/null)]"
fi
run_in_tree own-build-over-unusable-target "$TMP/no-build-target" "${IDENTITY[@]}" --mode stored --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}" --reference-dir "$TMP/reference"
expect_code own-build-over-unusable-target 0
expect_log "own-build-over-unusable-target: a file under CARGO_TARGET_DIR that is not executable is passed over" own.log \
    "$CASE/work/ensembl_cache/made_up_genome/63_ASM1 $CASE/work/json_cache"
CASE="$TMP/case/relative-target"
rm -rf "$CASE" "$CTL"
mkdir -p "$CASE" "$CTL"
(cd "$TMP" && PATH="$BIN:$PATH" CARGO_TARGET_DIR=cargo-target bash "$JOB" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" \
    "${IMAGE[@]}" "${BINARIES[@]}" --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr")
CODE=$?
expect_code relative-target 0
expect_log "relative-target: a relative CARGO_TARGET_DIR is mounted by its absolute path" docker.log \
    "-v $STUB_BUILD:/tools/vep-storable-to-json:ro made-up/image:tag"

usage_case() { # usage_case <name> <text of the error> <arguments...>
    local name=$1 text=$2
    shift 2
    run_case "$name" "$@"
    if [[ $CODE -eq 2 ]] && grep -q -F -- "$text" "$CASE/stderr" && [[ ! -e "$CASE/out/row.tsv" ]]; then
        ok "usage: $name"
    else
        bad "usage: $name (exit $CODE: $(cat "$CASE/stderr"))"
    fi
}
usage_case "a mode is needed" "--mode parity, stored or recheck is needed" "${IDENTITY[@]}" --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "an unknown mode" "--mode is parity, stored or recheck" "${IDENTITY[@]}" --mode compare --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}"
usage_case "no cache" "one of --cache-tarball, --cache-dir and --json-cache is needed" "${IDENTITY[@]}" --mode parity "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "two caches" "name one cache" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" --cache-dir "$SRC/plain" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "a converted cache beside an archive" "--json-cache accompanies --cache-dir alone" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" --json-cache "$TMP/ready_json" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "parity without the Ensembl cache" "--mode parity needs --cache-tarball or --cache-dir" "${IDENTITY[@]}" --mode parity --json-cache "$TMP/ready_json" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "parity without an image" "--mode parity needs --vep-image" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}"
usage_case "stored without a reference" "--mode stored needs --reference-dir" "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" "${BINARIES[@]}"
usage_case "recheck without a table" "--mode recheck needs --table-dir" "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" "${BINARIES[@]}"
usage_case "an unknown transcript set" "--transcripts is ensembl, refseq or merged" --genome made_up_genome --assembly ASM1 --transcripts gencode --set fungi --cache-version 63 --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "a cache version that is no number" "--cache-version is a whole number" --genome made_up_genome --assembly ASM1 --transcripts ensembl --set fungi --cache-version latest --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "a seed with a leading zero" "--seed is a whole number without a leading zero" "${IDENTITY[@]}" --mode parity --seed 007 --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "a genome with a slash" "--genome is needed, without a blank or a slash" --genome made/up --assembly ASM1 --transcripts ensembl --set fungi --cache-version 63 --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "an unknown stage" "--stop-after is extract, convert or generate" "${IDENTITY[@]}" --mode parity --stop-after compare --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "a missing archive" "no file $TMP/absent.tar.gz" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/absent.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
usage_case "a binary that is not there" "--vep-rs: $TOOLS/absent is not an executable file" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" --vep-rs "$TOOLS/absent" --generator "$TOOLS/vep-seeded-variants" --comparator "$TOOLS/vep-outputs-compare" --table-tool "$TOOLS/vep-genome-table"
usage_case "a missing comparator" "--comparator is needed" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" --vep-rs "$TOOLS/vep" --generator "$TOOLS/vep-seeded-variants" --table-tool "$TOOLS/vep-genome-table"
usage_case "an unknown argument" "unknown argument --fork" "${IDENTITY[@]}" --mode parity --fork 4 --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
CASE="$TMP/case/no-value"
mkdir -p "$CASE"
bash "$JOB" "${IDENTITY[@]}" --mode parity --work "$CASE/work" --out "$CASE/out" --seed >"$CASE/stdout" 2>"$CASE/stderr"
CODE=$?
if [[ $CODE -eq 2 ]] && grep -q -F -- "--seed takes a value" "$CASE/stderr"; then
    ok "usage: an option without its value"
else
    bad "usage: an option without its value (exit $CODE: $(cat "$CASE/stderr"))"
fi
if bash "$JOB" --help | grep -F -- '--mode parity|stored|recheck' >/dev/null; then
    ok "usage: --help names the three modes"
else
    bad "usage: --help does not name the three modes"
fi

# 7. Stored mode: vep-rs alone against the kept output; records that differ stop it.
run_case stored "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code stored 0
expect_row stored "$(row_of 3 exact "$DIGEST")"
if [[ ! -e "$CTL/docker.log" ]]; then ok "stored: neither the conversion nor Ensembl VEP runs"; else bad "stored: docker ran ($(cat "$CTL/docker.log"))"; fi
expect_log "stored: vep-rs reads the converted cache as given" vep_rs.log "--json_cache $TMP/ready_json --species made_up_genome"
expect_log "stored: the comparator reads the kept output" comparator.log "--perl-dir $CASE/work/ensembl_vep --rust-dir $CASE/work/vep_rs"
expect_no_log "stored: a run that holds no Ensembl cache names none to the comparator" comparator.log "--ensembl-cache"
if [[ $(grep -c -v '^#' "$CASE/work/ensembl_vep/variants.txt") -eq 1 ]]; then
    ok "stored: the kept Ensembl VEP output is the one compared"
else
    bad "stored: the compared Ensembl VEP output is not the kept one"
fi
expect_no_file stored out/ensembl_vep.txt.gz
run_case stored-other-seed "${IDENTITY[@]}" --mode stored --seed 2 --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code stored-other-seed 3
expect_stderr "stored-other-seed: the error says the stored input is not reproduced" "the generator does not reproduce the stored input"
expect_no_file stored-other-seed out/row.tsv
expect_no_log "stored-other-seed: vep-rs is not run" vep_rs.log "--offline"
# A stored input of another generator version is held by its records, never by the
# header line that names the version.
if ! grep -q -F "generator version" "$TMP/case/stored/stdout" "$TMP/case/stored/stderr"; then
    ok "stored: a run whose generator is the stored input's version says nothing of versions"
else
    bad "stored: a run whose generator is the stored input's version names a version"
fi
set_mode generator_version 2
run_case stored-other-version "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
if [[ $(gzip -dc "$TMP/reference/variants.vcf.gz" | diff - "$CASE/out/variants.vcf" | grep '^[<>]') == "$(printf '< ##source=vep-seeded-variants generator_version=1\n> ##source=vep-seeded-variants generator_version=2')" ]]; then
    ok "stored-other-version: the generated input differs from the stored one in the version line alone"
else
    bad "stored-other-version: the generated input does not differ from the stored one in the version line alone"
fi
expect_code stored-other-version 0
expect_row stored-other-version "$(row_of 3 exact "$DIGEST")"
expect_log "stored-other-version: vep-rs annotates the generated records" vep_rs.log "-i $CASE/out/variants.vcf -o $CASE/work/vep_rs/variants.txt"
if grep -q -x -F "check_genome: the generated records equal the stored ones; the stored input names generator version 1 and this run's generator is version 2" "$CASE/stdout"; then
    ok "stored-other-version: one line names both generator versions"
else
    bad "stored-other-version: no line names both generator versions ($(cat "$CASE/stdout"))"
fi
if python3 -c 'import json, sys; job = json.load(open(sys.argv[1])); assert job["generator_version"] == "2" and job["outcome"] == "exact", job' "$CASE/out/job.json" 2>/dev/null; then
    ok "stored-other-version: job.json holds this run's generator version"
else
    bad "stored-other-version: job.json ($(head -c 300 "$CASE/out/job.json" 2>/dev/null))"
fi
set_mode generator_version 2
run_case stored-other-version-other-records "${IDENTITY[@]}" --mode stored --seed 2 --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code stored-other-version-other-records 3
expect_stderr "stored-other-version-other-records: the error names both generator versions and says the stored input is not reproduced" \
    "the generated records differ from $TMP/reference/variants.vcf.gz; the stored input names generator version 1 and this run's generator is version 2: the generator does not reproduce the stored input"
expect_no_file stored-other-version-other-records out/row.tsv
expect_no_log "stored-other-version-other-records: vep-rs is not run" vep_rs.log "--offline"
mkdir -p "$TMP/reference-other-header"
cp "$TMP/reference/ensembl_vep.txt.gz" "$TMP/reference-other-header/"
printf '##fileformat=VCFv4.2\n##made_up_line=1\n##another_made_up_line=2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\tfirst_base_1\tA\tG\t.\t.\t.\n1\t200\tintron_1\tC\tT\t.\t.\t.\n1\t300\tupstream_1\tG\tA\t.\t.\t.\n' |
    gzip -n -c >"$TMP/reference-other-header/variants.vcf.gz"
run_case stored-other-header "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference-other-header" "${BINARIES[@]}"
expect_code stored-other-header 0
expect_row stored-other-header "$(row_of 3 exact "$DIGEST")"
if grep -q -x -F "check_genome: the generated records equal the stored ones; the stored input names no generator version and this run's generator is version 1" "$CASE/stdout"; then
    ok "stored-other-header: the line says the stored input names no generator version"
else
    bad "stored-other-header: no line says the stored input names no generator version ($(cat "$CASE/stdout"))"
fi
printf '##fileformat=VCFv4.2\n##source=vep-seeded-variants generator_version=1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\tfirst_base_1\tA\tG\t.\t.\t.\n1\t200\tintron_1\tC\tT\t.\t.\t.\n' |
    gzip -n -c >"$TMP/reference-other-header/variants.vcf.gz"
run_case stored-fewer-records "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference-other-header" "${BINARIES[@]}"
expect_code stored-fewer-records 3
expect_stderr "stored-fewer-records: a stored input that lacks a generated record is not reproduced" \
    "the generated records differ from $TMP/reference-other-header/variants.vcf.gz: the generator does not reproduce the stored input"
echo "not a gzip file" >"$TMP/reference-other-header/variants.vcf.gz"
run_case stored-unreadable "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference-other-header" "${BINARIES[@]}"
expect_code stored-unreadable 3
expect_stderr "stored-unreadable: the error says the stored input cannot be read" "cannot read $TMP/reference-other-header/variants.vcf.gz"
expect_no_file stored-unreadable out/row.tsv
# A null byte in the first record of both inputs: the second records, which differ,
# are still compared.
mkdir -p "$TMP/reference-null-byte"
cp "$TMP/reference/ensembl_vep.txt.gz" "$TMP/reference-null-byte/"
printf '##fileformat=VCFv4.2\n##source=vep-seeded-variants generator_version=1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\tfirst_base_1\0\tA\tG\t.\t.\t.\n1\t200\tintron_2\tC\tT\t.\t.\t.\n' |
    gzip -n -c >"$TMP/reference-null-byte/variants.vcf.gz"
set_mode generator_mode null_byte
run_case stored-null-byte "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference-null-byte" "${BINARIES[@]}"
expect_code stored-null-byte 3
expect_stderr "stored-null-byte: records that differ after a null byte in both inputs are not reproduced" \
    "the generated records differ from $TMP/reference-null-byte/variants.vcf.gz: the generator does not reproduce the stored input"
set_mode compare_mode differ
run_case stored-differs "${IDENTITY[@]}" --mode stored --json-cache "$TMP/ready_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code stored-differs 1
expect_row stored-differs "$(row_of 3 differs:unclassified -)"
run_case stored-converts "${IDENTITY[@]}" --mode stored --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code stored-converts 0
expect_log "stored-converts: an archive is converted in the image" docker.log "made-up/image:tag /tools/vep-storable-to-json /input/made_up_genome/63_ASM1 /output"
expect_no_log "stored-converts: Ensembl VEP is not run" docker.log " vep --offline"
run_case stored-host "${IDENTITY[@]}" --mode stored --cache-tarball "$TMP/plain.tar.gz" "${BINARIES[@]}" --reference-dir "$TMP/reference"
expect_code stored-host 0
expect_log "stored-host: without an image the conversion runs the build on the host" converter.log \
    "$CASE/work/ensembl_cache/made_up_genome/63_ASM1 $CASE/work/json_cache"
if [[ ! -e "$CTL/docker.log" ]]; then
    ok "stored-host: docker is not started"
else
    bad "stored-host: docker ran ($(cat "$CTL/docker.log"))"
fi

# 8. Recheck: the digest and the record count against the genome's row.
recheck_row() { # recheck_row <outcome> <sha256> <records> <seed>
    set_mode row.outcome "$1"
    set_mode row.sha256 "$2"
    set_mode row.records "$3"
    set_mode row.seed "$4"
}
recheck_row exact "$DIGEST" 3 1
run_case recheck-match "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-match 0
expect_row recheck-match "$(row_of 3 exact "$DIGEST")"
expect_log "recheck: the row is read by the genome key" table.log "row $TABLE $KEY --field outcome"
if [[ ! -e "$CTL/docker.log" && ! -e "$CTL/comparator.log" ]]; then
    ok "recheck: neither Ensembl VEP nor the comparator runs"
else
    bad "recheck: docker or the comparator ran"
fi
recheck_row exact "$DIGEST" 3 7
run_case recheck-row-seed "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-row-seed 0
expect_log "recheck-row-seed: the records are generated under the row's seed" generator.log "--seed 7 --output"
expect_row recheck-row-seed "$(row_of 3 exact "$DIGEST" 7)"
recheck_row exact "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb" 3 1
run_case recheck-other-digest "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-other-digest 1
expect_no_file recheck-other-digest out/row.tsv
expect_stderr "recheck-other-digest: the error gives both digests" "digest bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb; this run made 3 records and digest $DIGEST"
recheck_row exact "$DIGEST" 30 1
run_case recheck-other-count "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-other-count 1
recheck_row exact "$DIGEST" 3 1
set_mode vep_rs_mode fail
run_case recheck-vep-rs-fails "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-vep-rs-fails 1
expect_no_file recheck-vep-rs-fails out/row.tsv
recheck_row differs:unclassified - 3 1
run_case recheck-not-exact "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-not-exact 2
expect_stderr "recheck-not-exact: the skip names the row's outcome" "skip: the row of $KEY is differs:unclassified"
set_mode row.missing ""
run_case recheck-no-row "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-no-row 3
expect_no_file recheck-no-row out/row.tsv
recheck_row exact "$DIGEST" 3 ""
run_case recheck-no-seed "${IDENTITY[@]}" --mode recheck --json-cache "$TMP/ready_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-no-seed 3
if python3 -c 'import json, sys; job = json.load(open(sys.argv[1])); assert job["seed"] is None and job["outcome"] == "failed", job' "$CASE/out/job.json" 2>/dev/null; then
    ok "recheck-no-seed: job.json parses and holds no seed"
else
    bad "recheck-no-seed: job.json ($(head -c 200 "$CASE/out/job.json" 2>/dev/null))"
fi

# 9. --stop-after at each stage, and the run that carries on from it.
run_case stop-extract "${IDENTITY[@]}" --stop-after extract --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}"
expect_code stop-extract 0
expect_no_file stop-extract out/row.tsv
expect_no_file stop-extract work/json_cache
if [[ ! -e "$CTL/docker.log" && ! -e "$CTL/generator.log" ]]; then ok "stop-extract: nothing runs past the extraction"; else bad "stop-extract: a later step ran"; fi
if [[ $(cat "$CASE/work/layout.tsv") == "$(printf 'ensembl_cache_root\t%s\nensembl_cache_dir\t%s' "$CASE/work/ensembl_cache" "$CASE/work/ensembl_cache/made_up_genome/63_ASM1")" ]]; then
    ok "stop-extract: layout.tsv names the extracted cache"
else
    bad "stop-extract: layout.tsv is [$(cat "$CASE/work/layout.tsv")]"
fi
if [[ $(tail -n 1 "$CASE/stdout") == "check_genome: $KEY stage extract" ]]; then ok "stop-extract: the last line names the stage"; else bad "stop-extract: the last line is [$(tail -n 1 "$CASE/stdout")]"; fi
run_case stop-extract-absent-binary "${IDENTITY[@]}" --stop-after extract --cache-tarball "$TMP/plain.tar.gz" --vep-rs "$TOOLS/absent"
expect_code stop-extract-absent-binary 0
if [[ ! -s "$CASE/stderr" ]]; then
    ok "stop-extract-absent-binary: a binary the stage does not run, named and not there, prints nothing"
else
    bad "stop-extract-absent-binary: standard error holds [$(cat "$CASE/stderr")]"
fi
run_case stop-convert "${IDENTITY[@]}" --stop-after convert --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}"
expect_code stop-convert 0
expect_file stop-convert work/json_cache/info.json
expect_no_file stop-convert out/variants.vcf
expect_no_file stop-convert out/row.tsv
# The second run carries on in the first one's directories, so it is not a run_case.
rm -rf "$CTL"
mkdir -p "$CTL"
PATH="$BIN:$PATH" bash "$JOB" "${IDENTITY[@]}" --mode parity --cache-dir "$CASE/work/ensembl_cache" --json-cache "$CASE/work/json_cache" \
    "${IMAGE[@]}" "${BINARIES[@]}" --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr"
CODE=$?
expect_code carry-on 0
expect_row carry-on "$(row_of 3 exact "$DIGEST")"
expect_no_log "carry-on: a converted cache beside the Ensembl cache is not converted again" docker.log "vep-storable-to-json"
expect_log "carry-on: Ensembl VEP reads the cache directory as given" docker.log "-v $CASE/work/ensembl_cache:/cache:ro"
if python3 - "$CASE/out/job.json" <<'PY'; then
import json, sys
steps = [step["step"] for step in json.load(open(sys.argv[1]))["steps"]]
assert steps == ["extract", "convert", "generate", "ensembl_vep", "vep_rs", "digest", "compare"], steps
PY
    ok "carry-on: job.json keeps the earlier run's steps"
else
    bad "carry-on: job.json lost the earlier run's steps"
fi
run_case stop-generate "${IDENTITY[@]}" --mode parity --stop-after generate --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" --generator "$TOOLS/vep-seeded-variants"
expect_code stop-generate 0
expect_file stop-generate out/variants.vcf
expect_no_file stop-generate out/row.tsv
expect_no_log "stop-generate: Ensembl VEP is not run" docker.log " vep --offline"

# 10. A genome kept in a collection: the species directory one level down.
run_case collection "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/collection.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code collection 0
expect_log "collection: Ensembl VEP's cache root is the collection directory" docker.log \
    "-v $CASE/work/ensembl_cache/made_up_collection:/cache:ro"

# 11. The seed is passed on and written to the row.
run_case seed "${IDENTITY[@]}" --mode parity --seed 9 --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}"
expect_log "seed: the generator is given --seed 9" generator.log "--seed 9 --output"
expect_row seed "$(row_of 3 exact "$DIGEST" 9)"

# 12. A run writes under --work and --out alone: the listing of its home, TMPDIR and
# working directories is the same before and after.
SANDBOX="$TMP/sandbox"
rm -rf "$CTL" "$SANDBOX"
mkdir -p "$CTL" "$SANDBOX/home" "$SANDBOX/tmp" "$SANDBOX/cwd" "$SANDBOX/cache"
cp "$TMP/plain.tar.gz" "$SANDBOX/cache/"
listing() { (cd "$SANDBOX" && find . -path ./work -prune -o -path ./out -prune -o -print | sort); }
before=$(listing)
(cd "$SANDBOX/cwd" && HOME="$SANDBOX/home" TMPDIR="$SANDBOX/tmp" PATH="$BIN:$PATH" bash "$JOB" "${IDENTITY[@]}" --mode parity \
    --cache-tarball "$SANDBOX/cache/plain.tar.gz" "${IMAGE[@]}" "${BINARIES[@]}" --work "$SANDBOX/work" --out "$SANDBOX/out" >/dev/null 2>&1)
CODE=$?
if [[ $CODE -eq 0 && $(listing) == "$before" ]]; then
    ok "a run writes under --work and --out alone"
else
    bad "a run wrote outside --work and --out (exit $CODE)"
fi

# 13. A binary named without a slash is the file in the working directory, whatever
# program PATH holds under that name.
CASE="$TMP/case/bare-name"
rm -rf "$CASE" "$CTL"
mkdir -p "$CASE" "$CTL" "$TMP/other-vep"
cat >"$TMP/other-vep/vep" <<'STUB'
#!/usr/bin/env bash
echo "$*" >>"$CTL/path_vep.log"
STUB
chmod +x "$TMP/other-vep/vep"
(cd "$TOOLS" && PATH="$TMP/other-vep:$BIN:$PATH" bash "$JOB" "${IDENTITY[@]}" --mode parity --cache-tarball "$TMP/plain.tar.gz" "${IMAGE[@]}" \
    --vep-rs vep --generator vep-seeded-variants --comparator vep-outputs-compare --table-tool vep-genome-table \
    --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr")
CODE=$?
expect_code bare-name 0
expect_row bare-name "$(row_of 3 exact "$DIGEST")"
if [[ ! -e "$CTL/path_vep.log" ]]; then
    ok "bare-name: the vep on PATH is not run"
else
    bad "bare-name: the vep on PATH ran ($(cat "$CTL/path_vep.log"))"
fi

# 14. A name that holds a pattern character reaches Ensembl VEP as written, whatever
# the working directory holds.
mkdir -p "$SRC/pattern/made_up_genome/63_ASM?/1" "$TMP/pattern-cwd"
: >"$TMP/pattern-cwd/ASM1"
: >"$TMP/pattern-cwd/ASM2"
CASE="$TMP/case/pattern"
rm -rf "$CASE" "$CTL"
mkdir -p "$CASE" "$CTL"
(cd "$TMP/pattern-cwd" && PATH="$BIN:$PATH" bash "$JOB" --genome made_up_genome --assembly 'ASM?' --transcripts ensembl --set fungi --cache-version 63 \
    --mode parity --cache-dir "$SRC/pattern" "${IMAGE[@]}" "${BINARIES[@]}" --work "$CASE/work" --out "$CASE/out" >"$CASE/stdout" 2>"$CASE/stderr")
CODE=$?
expect_code pattern 0
expect_log "pattern: the assembly reaches Ensembl VEP as one word, as written" docker.log "--species made_up_genome --assembly ASM? --cache_version 63"

# 15. A cache that names an alignment file, whose edits its transcripts carry: both
# engines are given --use_given_ref, in every mode, and on no other cache. Without the
# option the stub engines refuse such a cache, as Ensembl VEP and vep-rs do where the
# command names no FASTA file.
EDITED_IDENTITY=(--genome made_up_genome --assembly ASM1 --transcripts refseq --set vertebrates --cache-version 116)
mkdir -p "$SRC/edited/made_up_genome_refseq/116_ASM1/1" "$SRC/edited-merged/made_up_genome_merged/116_ASM1/1" \
    "$TMP/edited_json" "$TMP/unedited_json"
printf 'species\tmade_up_genome\nassembly\tASM1\nbam\tmade_up_alignments.bam\n' >"$SRC/edited/made_up_genome_refseq/116_ASM1/info.txt"
cp "$SRC/edited/made_up_genome_refseq/116_ASM1/info.txt" "$SRC/edited-merged/made_up_genome_merged/116_ASM1/info.txt"
# info.json as the conversion lays it out: one key a line, a blank each side of the colon.
printf '{\n   "assembly" : "ASM1",\n   "bam" : "made_up_alignments.bam",\n   "source_type" : "refseq"\n}\n' >"$TMP/edited_json/info.json"
# A `bam` value of 0 is false to both engines and names no alignment file.
printf '{\n   "assembly" : "ASM1",\n   "bam" : "0",\n   "source_type" : "refseq"\n}\n' >"$TMP/unedited_json/info.json"
run_case edited "${EDITED_IDENTITY[@]}" --mode parity --cache-dir "$SRC/edited" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code edited 0
expect_row edited "$(printf 'made_up_genome\tASM1\trefseq\tvertebrates\t116\t1\t3\texact\t%s' "$DIGEST")"
expect_log "edited: Ensembl VEP is given --use_given_ref after the set flag" docker.log \
    "--assembly ASM1 --cache_version 116 --refseq --use_given_ref -i /input/variants.vcf"
expect_log "edited: vep-rs is given --use_given_ref after the set flag" vep_rs.log \
    "--assembly ASM1 --refseq --use_given_ref --warning_file"
expect_file edited out/ensembl_vep.txt.gz
if python3 - "$CASE/out/job.json" <<'PY'; then
import json, sys
job = json.load(open(sys.argv[1]))
assert job["outcome"] == "exact", job
assert " --refseq --use_given_ref -i /input/variants.vcf " in job["ensembl_vep_command"], job
assert job["vep_rs_command"].endswith(" --refseq --use_given_ref --warning_file <log> --no_stats"), job
PY
    ok "edited: job.json records the option in both engines' commands"
else
    bad "edited: job.json ($(head -c 600 "$CASE/out/job.json" 2>/dev/null))"
fi
run_case edited-merged --genome made_up_genome --assembly ASM1 --transcripts merged --set vertebrates --cache-version 116 \
    --mode parity --cache-dir "$SRC/edited-merged" "${IMAGE[@]}" "${BINARIES[@]}"
expect_code edited-merged 0
expect_log "edited-merged: Ensembl VEP is given --use_given_ref after --merged" docker.log \
    "--assembly ASM1 --cache_version 116 --merged --use_given_ref -i /input/variants.vcf"
expect_log "edited-merged: vep-rs is given --use_given_ref after --merged" vep_rs.log \
    "--assembly ASM1 --merged --use_given_ref --warning_file"
run_case edited-stored "${EDITED_IDENTITY[@]}" --mode stored --json-cache "$TMP/edited_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code edited-stored 0
expect_row edited-stored "$(printf 'made_up_genome\tASM1\trefseq\tvertebrates\t116\t1\t3\texact\t%s' "$DIGEST")"
expect_log "edited-stored: the comparator is given the converted cache a stored run on a RefSeq cache names" comparator.log \
    "--no-canonical-contigs --assembly ASM1 --json-cache $TMP/edited_json"
expect_log "edited-stored: vep-rs is given --use_given_ref on a converted cache that names an alignment file" vep_rs.log \
    "--assembly ASM1 --refseq --use_given_ref --warning_file"
recheck_row exact "$DIGEST" 3 1
run_case recheck-edited "${EDITED_IDENTITY[@]}" --mode recheck --json-cache "$TMP/edited_json" --table-dir "$TABLE" "${BINARIES[@]}"
expect_code recheck-edited 0
expect_log "recheck-edited: vep-rs is given --use_given_ref" vep_rs.log "--assembly ASM1 --refseq --use_given_ref --warning_file"
run_case unedited-stored "${EDITED_IDENTITY[@]}" --mode stored --json-cache "$TMP/unedited_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code unedited-stored 0
expect_no_log "unedited-stored: a bam value of 0 leaves vep-rs without --use_given_ref" vep_rs.log "--use_given_ref"
# A `bam` key of an inner object names no alignment file: under source_versions it is
# the version of a source of that name.
mkdir -p "$TMP/inner_json" "$TMP/inner_first_json"
printf '{\n   "assembly" : "ASM1",\n   "source_type" : "refseq",\n   "source_versions" : {\n      "bam" : "1"\n   }\n}\n' >"$TMP/inner_json/info.json"
run_case inner-stored "${EDITED_IDENTITY[@]}" --mode stored --json-cache "$TMP/inner_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code inner-stored 0
expect_no_log "inner-stored: a bam key under source_versions leaves vep-rs without --use_given_ref" vep_rs.log "--use_given_ref"
# The top-level key is read past an inner one that stands before it, here on one line.
printf '{"source_versions": {"bam": "0"}, "bam": "made_up_alignments.bam", "source_type": "refseq"}\n' >"$TMP/inner_first_json/info.json"
run_case inner-first-stored "${EDITED_IDENTITY[@]}" --mode stored --json-cache "$TMP/inner_first_json" --reference-dir "$TMP/reference" "${BINARIES[@]}"
expect_code inner-first-stored 0
expect_log "inner-first-stored: vep-rs is given --use_given_ref by the top-level bam key" vep_rs.log \
    "--assembly ASM1 --refseq --use_given_ref --warning_file"

echo "passed: $PASS   failed: $FAIL"
[[ $FAIL -eq 0 ]]
