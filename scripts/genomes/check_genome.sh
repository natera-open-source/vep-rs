#!/usr/bin/env bash
# Checks one genome's ready-made Ensembl VEP cache: extracts it, converts it with the
# release build of tools/vep-storable-to-json, generates the genome's seeded input
# records, annotates them with Ensembl VEP and with vep-rs under default flags in the
# default output format, compares the two outputs and writes the genome's row of
# tests/genomes/. On a cache that names an alignment file, whose edits its transcripts
# carry, both engines are given --use_given_ref. docs/genome-checks.md describes the
# row and the files left behind.
#
# Every path is local: the script reads no remote store and writes under --work and
# --out alone. It runs under bash 3.2.
#
# Run: bash scripts/genomes/check_genome.sh --help
set -uo pipefail
export LC_ALL=C

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "$HERE/../.." && pwd)
CONVERTER_NAME=vep-storable-to-json
ROW_HEADER=$(printf 'genome\tassembly\ttranscripts\tset\tcache_version\tseed\trecords\toutcome\tsha256')

usage() {
    cat <<'USAGE'
Usage: check_genome.sh [--mode parity|stored|recheck]
           --cache-tarball <file> | --cache-dir <dir> | --json-cache <dir>
           --genome <name> --assembly <name> --transcripts ensembl|refseq|merged
           --set <set> --cache-version <n> [--seed <n>]
           [--vep-image <reference>]
           [--vep-rs <vep>] [--generator <vep-seeded-variants>]
           [--comparator <vep-outputs-compare>] [--table-tool <vep-genome-table>]
           [--reference-dir <dir>] [--table-dir <dir>]
           [--stop-after extract|convert|generate]
           --work <dir> --out <dir>

The cache, one of:
  --cache-tarball <file>  a cache archive as Ensembl ships it
  --cache-dir <dir>       an extracted cache: the directory that holds
                          <species>/<cache version>_<assembly>, the species directory
                          being the genome, suffixed _refseq or _merged by --transcripts
  --json-cache <dir>      a converted cache, the directory that holds info.json; it may
                          accompany --cache-dir, and the conversion is then not run

Modes:
  parity   runs Ensembl VEP and vep-rs and compares their outputs; needs the Ensembl
           cache and --vep-image; leaves ensembl_vep.txt.gz and
           ensembl_vep.warnings.txt in --out
  stored   runs vep-rs alone and compares its output with ensembl_vep.txt.gz of
           --reference-dir, after checking that the generated records equal those of
           variants.vcf.gz there; the header lines, which name the generator's
           version, are not compared
  recheck  runs vep-rs alone and compares its digest and record count with the row
           --table-dir holds for the genome; the seed is the row's unless --seed is given

--stop-after ends the run with exit 0 once the named stage is done; --mode may be left
out with --stop-after extract or convert. A later run on the same --work and --out,
given --cache-dir <ensembl_cache_root of layout.tsv> and, after a conversion,
--json-cache <json_cache of layout.tsv>, carries on without repeating the stage.

A cache whose info.txt names an alignment file under `bam` holds transcripts edited
from that file; the converted cache's info.json keeps the line. On such a cache both
engines are given --use_given_ref, without which each refuses a run that names no
FASTA file, and no other cache changes either command.

The conversion and Ensembl VEP run in --vep-image under docker; without --vep-image the
conversion runs on the host. Either way the converter is the release build of
tools/vep-storable-to-json, which is mounted read-only into the image: the build under
$CARGO_TARGET_DIR/release when that holds one, else the one under the tool's own
target/release. This script builds nothing; the build command is
  cargo build --release --locked --manifest-path tools/vep-storable-to-json/Cargo.toml

Binaries: parity and stored take all four; recheck takes --vep-rs, --generator and
--table-tool; --stop-after generate takes --generator.

Left in --out: row.tsv (the table's header line and the genome's row), variants.vcf,
variants.summary.tsv, comparison/ (the report of vep-outputs-compare), job.json.
Left in --work: ensembl_cache/, json_cache/, logs/, steps.tsv and layout.tsv, the
names and absolute paths of ensembl_cache_root, ensembl_cache_dir, json_cache and
sequence_regions.

Exit codes:
  0  exact, a recheck match, or the --stop-after stage reached
  1  the outputs differ (the row says differs:unclassified), or a recheck mismatch
  2  a usage error, or a skip: no build of the converter, no docker, no image, or
     a recheck of a row that is not exact
  3  a step failed; the row is written when the failure is a named exclusion
     (cache-unreadable, serialiser-not-storable, no-transcripts, ensembl-vep-failed)
USAGE
}

usage_error() {
    echo "ERROR: [check_genome] $1" >&2
    echo "Run with --help for the usage." >&2
    exit 2
}

skip() {
    echo "ERROR: [check_genome] skip: $1" >&2
    exit 2
}

MODE=""
CACHE_TARBALL=""
CACHE_DIR=""
JSON_CACHE_GIVEN=""
GENOME=""
ASSEMBLY=""
TRANSCRIPTS=""
SET_NAME=""
CACHE_VERSION=""
SEED=""
VEP_IMAGE=""
VEP_RS=""
GENERATOR=""
COMPARATOR=""
TABLE_TOOL=""
REFERENCE_DIR=""
TABLE_DIR=""
STOP_AFTER=""
WORK=""
OUT=""

while [ $# -gt 0 ]; do
    case "$1" in
    --help | -h)
        usage
        exit 0
        ;;
    --mode | --cache-tarball | --cache-dir | --json-cache | --genome | --assembly | \
        --transcripts | --set | --cache-version | --seed | --vep-image | --vep-rs | \
        --generator | --comparator | --table-tool | --reference-dir | --table-dir | \
        --stop-after | --work | --out)
        [ $# -ge 2 ] || usage_error "$1 takes a value"
        case "$1" in
        --mode) MODE=$2 ;;
        --cache-tarball) CACHE_TARBALL=$2 ;;
        --cache-dir) CACHE_DIR=$2 ;;
        --json-cache) JSON_CACHE_GIVEN=$2 ;;
        --genome) GENOME=$2 ;;
        --assembly) ASSEMBLY=$2 ;;
        --transcripts) TRANSCRIPTS=$2 ;;
        --set) SET_NAME=$2 ;;
        --cache-version) CACHE_VERSION=$2 ;;
        --seed) SEED=$2 ;;
        --vep-image) VEP_IMAGE=$2 ;;
        --vep-rs) VEP_RS=$2 ;;
        --generator) GENERATOR=$2 ;;
        --comparator) COMPARATOR=$2 ;;
        --table-tool) TABLE_TOOL=$2 ;;
        --reference-dir) REFERENCE_DIR=$2 ;;
        --table-dir) TABLE_DIR=$2 ;;
        --stop-after) STOP_AFTER=$2 ;;
        --work) WORK=$2 ;;
        --out) OUT=$2 ;;
        esac
        shift 2
        ;;
    *) usage_error "unknown argument $1" ;;
    esac
done

# Digits without a leading zero, the one spelling the table tool takes in a row and a
# JSON file takes as a number.
is_whole_number() {
    case "$1" in
    '' | *[!0-9]* | 0?*) return 1 ;;
    *) return 0 ;;
    esac
}

# A name that goes into a path, a key and a command line: no blank, slash or control
# character.
is_plain_name() {
    case "$1" in
    '' | *[[:space:]/[:cntrl:]]*) return 1 ;;
    *) return 0 ;;
    esac
}

case "$MODE" in
parity | stored | recheck) ;;
'')
    case "$STOP_AFTER" in
    extract | convert) ;;
    *) usage_error "--mode parity, stored or recheck is needed" ;;
    esac
    ;;
*) usage_error "--mode is parity, stored or recheck, not $MODE" ;;
esac
case "$STOP_AFTER" in
'' | extract | convert | generate) ;;
*) usage_error "--stop-after is extract, convert or generate, not $STOP_AFTER" ;;
esac
case "$TRANSCRIPTS" in
ensembl) SPECIES_SUFFIX="" SET_FLAG="" ;;
refseq) SPECIES_SUFFIX="_refseq" SET_FLAG="--refseq" ;;
merged) SPECIES_SUFFIX="_merged" SET_FLAG="--merged" ;;
*) usage_error "--transcripts is ensembl, refseq or merged" ;;
esac
is_plain_name "$GENOME" || usage_error "--genome is needed, without a blank or a slash"
is_plain_name "$ASSEMBLY" || usage_error "--assembly is needed, without a blank or a slash"
is_plain_name "$SET_NAME" || usage_error "--set is needed, without a blank or a slash"
is_whole_number "$CACHE_VERSION" || usage_error "--cache-version is a whole number without a leading zero"
SEED_GIVEN=true
if [ -z "$SEED" ]; then
    SEED=1
    SEED_GIVEN=false
fi
is_whole_number "$SEED" || usage_error "--seed is a whole number without a leading zero"
[ -n "$WORK" ] || usage_error "--work is needed"
[ -n "$OUT" ] || usage_error "--out is needed"

if [ -n "$CACHE_TARBALL" ] && [ -n "$CACHE_DIR" ]; then
    usage_error "--cache-tarball and --cache-dir are given together; name one cache"
fi
if [ -n "$CACHE_TARBALL" ] && [ -n "$JSON_CACHE_GIVEN" ]; then
    usage_error "--json-cache accompanies --cache-dir alone, not --cache-tarball"
fi
if [ -z "$CACHE_TARBALL" ] && [ -z "$CACHE_DIR" ] && [ -z "$JSON_CACHE_GIVEN" ]; then
    usage_error "one of --cache-tarball, --cache-dir and --json-cache is needed"
fi
[ -z "$CACHE_TARBALL" ] || [ -f "$CACHE_TARBALL" ] || usage_error "no file $CACHE_TARBALL"
[ -z "$CACHE_DIR" ] || [ -d "$CACHE_DIR" ] || usage_error "no directory $CACHE_DIR"
if [ -n "$JSON_CACHE_GIVEN" ] && [ ! -f "$JSON_CACHE_GIVEN/info.json" ]; then
    usage_error "no info.json under $JSON_CACHE_GIVEN"
fi

HAS_ENSEMBL_CACHE=false
if [ -n "$CACHE_TARBALL" ] || [ -n "$CACHE_DIR" ]; then HAS_ENSEMBL_CACHE=true; fi
NEEDS_CONVERSION=false
if [ -z "$JSON_CACHE_GIVEN" ]; then NEEDS_CONVERSION=true; fi

case "$MODE" in
parity)
    [ "$HAS_ENSEMBL_CACHE" = true ] || usage_error "--mode parity needs --cache-tarball or --cache-dir"
    [ -n "$VEP_IMAGE" ] || usage_error "--mode parity needs --vep-image"
    ;;
stored)
    [ -n "$REFERENCE_DIR" ] || usage_error "--mode stored needs --reference-dir"
    for stored_file in ensembl_vep.txt.gz variants.vcf.gz; do
        [ -f "$REFERENCE_DIR/$stored_file" ] || usage_error "no $stored_file under $REFERENCE_DIR"
    done
    ;;
recheck)
    [ -n "$TABLE_DIR" ] || usage_error "--mode recheck needs --table-dir"
    [ -d "$TABLE_DIR" ] || usage_error "no directory $TABLE_DIR"
    ;;
esac

# need_binary <option> <path>: the binaries a run uses are named and executable.
need_binary() {
    [ -n "$2" ] || usage_error "$1 is needed"
    if [ ! -f "$2" ] || [ ! -x "$2" ]; then usage_error "$1: $2 is not an executable file"; fi
}
RUNS_GENERATE=true
RUNS_ENGINES=true
case "$STOP_AFTER" in
extract | convert) RUNS_GENERATE=false RUNS_ENGINES=false ;;
generate) RUNS_ENGINES=false ;;
esac
if [ "$RUNS_GENERATE" = true ]; then need_binary --generator "$GENERATOR"; fi
if [ "$RUNS_ENGINES" = true ]; then
    need_binary --vep-rs "$VEP_RS"
    need_binary --table-tool "$TABLE_TOOL"
    case "$MODE" in
    parity | stored) need_binary --comparator "$COMPARATOR" ;;
    esac
elif [ "$MODE" = recheck ]; then
    need_binary --table-tool "$TABLE_TOOL"
fi
# A command word without a slash is looked up on PATH, so a binary named without one is
# run as ./<name>: the file of that name in the working directory, which the checks
# above read.
in_working_directory() {
    case "$1" in
    '' | */*) printf '%s\n' "$1" ;;
    *) printf './%s\n' "$1" ;;
    esac
}
VEP_RS=$(in_working_directory "$VEP_RS")
GENERATOR=$(in_working_directory "$GENERATOR")
COMPARATOR=$(in_working_directory "$COMPARATOR")
TABLE_TOOL=$(in_working_directory "$TABLE_TOOL")

# The skips, before anything is written: what the run needs of the host.
RUNS_ENSEMBL_VEP=false
if [ "$MODE" = parity ] && [ "$RUNS_ENGINES" = true ]; then RUNS_ENSEMBL_VEP=true; fi
RUNS_CONVERSION=false
if [ "$NEEDS_CONVERSION" = true ] && [ "$STOP_AFTER" != extract ]; then RUNS_CONVERSION=true; fi
# The converter is the release build of tools/vep-storable-to-json, which nothing here
# builds: the one under CARGO_TARGET_DIR when that is set and holds one (cargo writes a
# build there, not under the tool's own target directory), else the one under the tool's
# own target directory. The image supplies none, so a conversion in the image needs it
# as much as one on the host.
CONVERTER="$REPO_DIR/tools/$CONVERTER_NAME/target/release/$CONVERTER_NAME"
if [ -n "${CARGO_TARGET_DIR:-}" ] && [ -f "$CARGO_TARGET_DIR/release/$CONVERTER_NAME" ] &&
    [ -x "$CARGO_TARGET_DIR/release/$CONVERTER_NAME" ]; then
    CONVERTER="$CARGO_TARGET_DIR/release/$CONVERTER_NAME"
fi
if [ "$RUNS_CONVERSION" = true ]; then
    if [ ! -f "$CONVERTER" ] || [ ! -x "$CONVERTER" ]; then
        skip "no build of the converter at $CONVERTER: build it with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/$CONVERTER_NAME/Cargo.toml'"
    fi
    # docker takes a mount's host path as a volume name unless it is absolute.
    CONVERTER="$(cd "$(dirname "$CONVERTER")" && pwd)/$CONVERTER_NAME"
fi
USES_IMAGE=false
if [ "$RUNS_ENSEMBL_VEP" = true ]; then USES_IMAGE=true; fi
if [ "$RUNS_CONVERSION" = true ] && [ -n "$VEP_IMAGE" ]; then USES_IMAGE=true; fi
if [ "$USES_IMAGE" = true ]; then
    command -v docker >/dev/null 2>&1 || skip "no docker on PATH"
    docker image inspect "$VEP_IMAGE" >/dev/null 2>&1 || skip "no image $VEP_IMAGE in docker"
fi

mkdir -p "$WORK/logs" "$OUT" || usage_error "cannot create $WORK or $OUT"
WORK=$(cd "$WORK" && pwd)
OUT=$(cd "$OUT" && pwd)
[ -z "$CACHE_DIR" ] || CACHE_DIR=$(cd "$CACHE_DIR" && pwd)
[ -z "$JSON_CACHE_GIVEN" ] || JSON_CACHE_GIVEN=$(cd "$JSON_CACHE_GIVEN" && pwd)

SPECIES_DIR="$GENOME$SPECIES_SUFFIX"
VERSION_DIR="${CACHE_VERSION}_${ASSEMBLY}"
GENOME_KEY="${SET_NAME}__${GENOME}__${ASSEMBLY}__${TRANSCRIPTS}"
STEPS_FILE="$WORK/steps.tsv"
[ -f "$STEPS_FILE" ] || : >"$STEPS_FILE"

ENSEMBL_ROOT=""
ENSEMBL_DIR=""
JSON_CACHE=""
RECORDS=0
RECORDS_SHA=""
GENERATOR_VERSION=""
ENSEMBL_VEP_COMMAND=""
VEP_RS_COMMAND=""

TIME_BIN=$(type -P time 2>/dev/null || true)
if [ -n "$TIME_BIN" ] && ! "$TIME_BIN" -f '%M' -o /dev/null true >/dev/null 2>&1; then
    TIME_BIN=""
fi

now() {
    if [ -n "${EPOCHREALTIME:-}" ]; then
        printf '%s\n' "$EPOCHREALTIME"
    else
        date +%s
    fi
}

# record_step <step> <start> <end> <peak resident kB or -> <exit status>
record_step() {
    local seconds
    seconds=$(awk -v a="$2" -v b="$3" 'BEGIN { printf "%.3f", b - a }')
    printf '%s\t%s\t%s\t%s\n' "$1" "$seconds" "$4" "$5" >>"$STEPS_FILE"
    echo "check_genome: $1 ${seconds}s exit $5"
}

# run_native <step> <command...>: the command's output goes to logs/<step>.log; the
# peak resident memory is GNU time's where the host has it.
run_native() {
    local step=$1
    shift
    local log="$WORK/logs/$step.log" rss_file="$WORK/logs/$step.rss" rss="-" start end status
    start=$(now)
    if [ -n "$TIME_BIN" ]; then
        "$TIME_BIN" -f '%M' -o "$rss_file" "$@" >"$log" 2>&1
        status=$?
        if [ -s "$rss_file" ]; then rss=$(tail -n 1 "$rss_file"); fi
        is_whole_number "$rss" || rss="-"
    else
        "$@" >"$log" 2>&1
        status=$?
    fi
    end=$(now)
    record_step "$step" "$start" "$end" "$rss" "$status"
    return "$status"
}

# run_in_image <step> <docker run arguments...>: runs as the calling user, never pulls.
# The peak resident memory is the high-water mark of the container's first process,
# read from the host's /proc while the container runs; a host that does not show the
# process there records none.
run_in_image() {
    local step=$1
    shift
    local log="$WORK/logs/$step.log" cid_file="$WORK/logs/$step.cid"
    local rss=0 pid="" mark start end status docker_pid
    rm -f "$cid_file"
    start=$(now)
    docker run --rm --pull never --cidfile "$cid_file" --user "$(id -u):$(id -g)" "$@" >"$log" 2>&1 &
    docker_pid=$!
    while kill -0 "$docker_pid" 2>/dev/null; do
        if [ -z "$pid" ] && [ -s "$cid_file" ]; then
            pid=$(docker inspect --format '{{.State.Pid}}' "$(cat "$cid_file")" 2>/dev/null) || pid=""
            if ! is_whole_number "$pid" || [ "$pid" -eq 0 ]; then pid=""; fi
        fi
        if [ -n "$pid" ] && [ -r "/proc/$pid/status" ]; then
            mark=$(awk '/^VmHWM:/ { print $2 }' "/proc/$pid/status" 2>/dev/null) || mark=""
            if is_whole_number "$mark" && [ "$mark" -gt "$rss" ]; then rss=$mark; fi
        fi
        sleep 0.2
    done
    wait "$docker_pid"
    status=$?
    end=$(now)
    rm -f "$cid_file"
    if [ "$rss" -eq 0 ]; then rss="-"; fi
    record_step "$step" "$start" "$end" "$rss" "$status"
    return "$status"
}

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{ print $1 }'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{ print $1 }'
    fi
}

json_string() {
    printf '%s' "$1" | awk 'BEGIN { ORS = "" } {
        gsub(/\\/, "\\\\"); gsub(/"/, "\\\""); gsub(/\t/, "\\t"); gsub(/[[:cntrl:]]/, "")
        if (NR > 1) print "\\n"
        print
    }'
}

# json_text <value>: a JSON string, or null for an empty value.
json_text() {
    if [ -n "$1" ]; then
        printf '"%s"' "$(json_string "$1")"
    else
        printf 'null'
    fi
}

json_number() {
    if is_whole_number "$1"; then printf '%s' "$1"; else printf 'null'; fi
}

# json_binary <path>: the binary's digest and the first line its --version prints, null
# where the binary has no --version; null alone for a path that is empty or names no file.
json_binary() {
    local version=""
    if [ -z "$1" ] || [ ! -f "$1" ]; then
        printf 'null'
        return
    fi
    version=$("$1" --version 2>/dev/null | head -n 1) || version=""
    printf '{"sha256": %s, "version": %s}' "$(json_text "$(sha256_of "$1")")" "$(json_text "$version")"
}

write_layout() {
    {
        [ -z "$ENSEMBL_ROOT" ] || printf 'ensembl_cache_root\t%s\n' "$ENSEMBL_ROOT"
        [ -z "$ENSEMBL_DIR" ] || printf 'ensembl_cache_dir\t%s\n' "$ENSEMBL_DIR"
        [ -z "$JSON_CACHE" ] || printf 'json_cache\t%s\n' "$JSON_CACHE"
        if [ -n "$JSON_CACHE" ] && [ -f "$JSON_CACHE/sequence_regions.tsv" ]; then
            printf 'sequence_regions\t%s\n' "$JSON_CACHE/sequence_regions.tsv"
        fi
    } >"$WORK/layout.tsv"
}

# write_job <exit code> <outcome or stage>: steps.tsv keeps every run on this --work, and
# a step's last line is the one reported.
write_job() {
    local image_id="" tarball_bytes="" converted_bytes="" steps
    if [ -n "$VEP_IMAGE" ] && command -v docker >/dev/null 2>&1; then
        image_id=$(docker image inspect --format '{{.Id}}' "$VEP_IMAGE" 2>/dev/null | head -n 1) || image_id=""
    fi
    if [ -n "$CACHE_TARBALL" ]; then tarball_bytes=$(wc -c <"$CACHE_TARBALL" | tr -d ' '); fi
    if [ -n "$JSON_CACHE" ] && [ -d "$JSON_CACHE" ]; then
        converted_bytes=$(du -sk "$JSON_CACHE" 2>/dev/null | awk '{ print $1 * 1024 }')
    fi
    steps=$(awk -F '\t' '
        !($1 in seen) { order[++n] = $1 }
        { seen[$1] = 1; seconds[$1] = $2; rss[$1] = $3 }
        END {
            for (i = 1; i <= n; i++) {
                step = order[i]
                printf "%s    {\"step\": \"%s\", \"seconds\": %s, \"peak_rss_kb\": %s}", (i > 1 ? ",\n" : ""), step, seconds[step], (rss[step] == "-" ? "null" : rss[step])
            }
        }' "$STEPS_FILE")
    {
        printf '{\n'
        printf '  "genome_key": %s,\n' "$(json_text "$GENOME_KEY")"
        printf '  "mode": %s,\n' "$(json_text "$MODE")"
        printf '  "seed": %s,\n' "$(json_number "$SEED")"
        printf '  "outcome": %s,\n' "$(json_text "$2")"
        printf '  "exit_code": %s,\n' "$1"
        printf '  "records": %s,\n' "$RECORDS"
        printf '  "records_sha256": %s,\n' "$(json_text "$RECORDS_SHA")"
        printf '  "generator_version": %s,\n' "$(json_text "$GENERATOR_VERSION")"
        printf '  "vep_image": %s,\n' "$(json_text "$VEP_IMAGE")"
        printf '  "vep_image_id": %s,\n' "$(json_text "$image_id")"
        printf '  "ensembl_vep_command": %s,\n' "$(json_text "$ENSEMBL_VEP_COMMAND")"
        printf '  "vep_rs_command": %s,\n' "$(json_text "$VEP_RS_COMMAND")"
        printf '  "tarball_bytes": %s,\n' "$(json_number "$tarball_bytes")"
        printf '  "converted_bytes": %s,\n' "$(json_number "$converted_bytes")"
        printf '  "binaries": {\n'
        printf '    "vep_rs": %s,\n' "$(json_binary "$VEP_RS")"
        printf '    "generator": %s,\n' "$(json_binary "$GENERATOR")"
        printf '    "comparator": %s,\n' "$(json_binary "$COMPARATOR")"
        printf '    "table_tool": %s\n' "$(json_binary "$TABLE_TOOL")"
        printf '  },\n'
        printf '  "steps": [\n%s\n  ]\n' "$steps"
        printf '}\n'
    } >"$OUT/job.json"
}

# finish <exit code> <outcome or stage>
finish() {
    write_layout
    write_job "$1" "$2"
    echo "check_genome: $GENOME_KEY $2"
    exit "$1"
}

# write_row <outcome> <sha256 or ->
write_row() {
    printf '%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$ROW_HEADER" \
        "$GENOME" "$ASSEMBLY" "$TRANSCRIPTS" "$SET_NAME" "$CACHE_VERSION" "$SEED" "$RECORDS" "$1" "$2" \
        >"$OUT/row.tsv"
}

# exclude <name> <what happened>: the genome is not compared; its row says why.
exclude() {
    echo "ERROR: [check_genome] $2" >&2
    write_row "excluded:$1" "-"
    finish 3 "excluded:$1"
}

# step_failed <what happened>: no row is written.
step_failed() {
    echo "ERROR: [check_genome] $1" >&2
    finish 3 "failed"
}

stage_reached() {
    if [ "$STOP_AFTER" = "$1" ]; then finish 0 "stage $1"; fi
}

# count_rows <output file>: the lines that are not header lines.
count_rows() {
    local rows
    rows=$(grep -c -v '^#' "$1" 2>/dev/null) || rows=0
    printf '%s\n' "$rows"
}

# alignment_file_of <info.json>: the `bam` string of the file's top-level object, the
# alignment file whose edits the cache's transcripts carry; nothing for a cache that
# names none. A `bam` key of an inner object is passed over: under source_versions it
# is the version of a source of that name, which neither engine reads as a file.
alignment_file_of() {
    awk '
        {
            n = length($0)
            for (i = 1; i <= n; i++) {
                c = substr($0, i, 1)
                if (in_string) {
                    if (c == "\\") i++
                    else if (c == "\"") in_string = 0
                } else if (c == "\"") {
                    in_string = 1
                    if (depth == 1 && match(substr($0, i), /^"bam"[[:space:]]*:[[:space:]]*"[^"]*"/)) {
                        value = substr($0, i, RLENGTH)
                        sub(/^"bam"[[:space:]]*:[[:space:]]*"/, "", value)
                        sub(/"$/, "", value)
                        print value
                        exit
                    }
                } else if (c == "{" || c == "[") depth++
                else if (c == "}" || c == "]") depth--
            }
        }' "$1" 2>/dev/null
}

rm -f "$OUT/row.tsv" "$OUT/job.json"

# Recheck reads the genome's row first: a row that is not exact has no digest to hold
# vep-rs's output against.
ROW_SHA=""
ROW_RECORDS=""
if [ "$MODE" = recheck ]; then
    row_outcome=$("$TABLE_TOOL" row "$TABLE_DIR" "$GENOME_KEY" --field outcome 2>"$WORK/logs/row.log") ||
        step_failed "$TABLE_DIR holds no readable row for $GENOME_KEY ($(tail -n 1 "$WORK/logs/row.log"))"
    if [ "$row_outcome" != exact ]; then
        echo "ERROR: [check_genome] skip: the row of $GENOME_KEY is $row_outcome, so there is no digest to recheck" >&2
        exit 2
    fi
    ROW_SHA=$("$TABLE_TOOL" row "$TABLE_DIR" "$GENOME_KEY" --field sha256 2>>"$WORK/logs/row.log") || ROW_SHA=""
    ROW_RECORDS=$("$TABLE_TOOL" row "$TABLE_DIR" "$GENOME_KEY" --field records 2>>"$WORK/logs/row.log") || ROW_RECORDS=""
    if [ "$SEED_GIVEN" = false ]; then
        SEED=$("$TABLE_TOOL" row "$TABLE_DIR" "$GENOME_KEY" --field seed 2>>"$WORK/logs/row.log") || SEED=""
        is_whole_number "$SEED" || step_failed "the row of $GENOME_KEY carries no seed"
    fi
fi

# Extract.
if [ -n "$CACHE_TARBALL" ]; then
    rm -rf "$WORK/ensembl_cache"
    mkdir -p "$WORK/ensembl_cache"
    run_native extract tar -xzf "$CACHE_TARBALL" -C "$WORK/ensembl_cache" ||
        exclude cache-unreadable "$CACHE_TARBALL does not extract ($(tail -n 1 "$WORK/logs/extract.log"))"
    CACHE_DIR="$WORK/ensembl_cache"
fi
if [ -n "$CACHE_DIR" ]; then
    # An archive of a genome kept in a collection holds the species directory one level
    # down, under the collection's name.
    for candidate in "$CACHE_DIR" "$CACHE_DIR"/*/ "$CACHE_DIR"/*/*/; do
        candidate=${candidate%/}
        if [ -d "$candidate/$SPECIES_DIR/$VERSION_DIR" ]; then
            ENSEMBL_ROOT=$candidate
            ENSEMBL_DIR="$candidate/$SPECIES_DIR/$VERSION_DIR"
            break
        fi
    done
    [ -n "$ENSEMBL_DIR" ] ||
        exclude cache-unreadable "the cache holds no directory $SPECIES_DIR/$VERSION_DIR"
    if [ -f "$ENSEMBL_DIR/info.txt" ] &&
        [ "$(awk -F '\t' '$1 == "serialiser_type" { print $2 }' "$ENSEMBL_DIR/info.txt" | tail -n 1)" = sereal ]; then
        exclude serialiser-not-storable "$SPECIES_DIR/$VERSION_DIR uses the sereal serialiser"
    fi
fi
stage_reached extract

# Convert.
if [ -n "$JSON_CACHE_GIVEN" ]; then
    JSON_CACHE=$JSON_CACHE_GIVEN
else
    JSON_CACHE="$WORK/json_cache"
    rm -rf "$JSON_CACHE"
    mkdir -p "$JSON_CACHE"
    if [ -n "$VEP_IMAGE" ]; then
        # The converter reads the transcript set from the name of the directory above
        # the version directory, so the cache is mounted under its species directory.
        run_in_image convert \
            -v "$ENSEMBL_DIR:/input/$SPECIES_DIR/$VERSION_DIR:ro" -v "$JSON_CACHE:/output" \
            -v "$CONVERTER:/tools/$CONVERTER_NAME:ro" \
            "$VEP_IMAGE" "/tools/$CONVERTER_NAME" "/input/$SPECIES_DIR/$VERSION_DIR" /output
    else
        run_native convert "$CONVERTER" "$ENSEMBL_DIR" "$JSON_CACHE"
    fi
    convert_status=$?
    convert_error=$(grep -F 'ERROR: [storable_to_json]' "$WORK/logs/convert.log" | tail -n 1)
    if [ "$convert_status" -ne 0 ] || [ ! -f "$JSON_CACHE/info.json" ]; then
        case "$convert_error" in
        *sereal*) exclude serialiser-not-storable "the conversion refused the cache: $convert_error" ;;
        *'no transcript shards'*) exclude no-transcripts "the conversion found no transcript: $convert_error" ;;
        '') step_failed "the conversion did not run to its end (exit $convert_status; $WORK/logs/convert.log)" ;;
        *) exclude cache-unreadable "the conversion failed: $convert_error" ;;
        esac
    fi
    if grep -q -E '^Done\. .* 0 total transcripts' "$WORK/logs/convert.log"; then
        exclude no-transcripts "the cache's shards hold no transcript"
    fi
fi
stage_reached convert

# On a cache that names an alignment file both engines switch --use_transcript_ref on
# (CacheDir.pm:412-414) and then refuse a run that names no FASTA file (Runner.pm:762),
# and this job names none. --use_given_ref leaves the switch off (CacheDir.pm:412). The
# value is tested as both engines test it: empty and 0 name no file.
REFERENCE_FLAG=""
case "$(alignment_file_of "$JSON_CACHE/info.json")" in
'' | 0) ;;
*) REFERENCE_FLAG="--use_given_ref" ;;
esac

# Generate.
run_native generate "$GENERATOR" --cache "$JSON_CACHE" --seed "$SEED" \
    --output "$OUT/variants.vcf" --summary "$OUT/variants.summary.tsv" ||
    step_failed "the generator failed ($(tail -n 1 "$WORK/logs/generate.log"))"
[ -f "$OUT/variants.vcf" ] || step_failed "the generator wrote no $OUT/variants.vcf"
RECORDS=$(count_rows "$OUT/variants.vcf")
RECORDS_SHA=$(sha256_of "$OUT/variants.vcf")
GENERATOR_VERSION=$(sed -n '2s/^##source=.* generator_version=\([0-9][0-9]*\).*$/\1/p' "$OUT/variants.vcf")
[ "$RECORDS" -gt 0 ] || step_failed "the generator wrote no record"
# The stored Ensembl VEP output stands for the stored records, so the records are held
# and the header is not: it names the generator version, and two versions give a cache
# the same records wherever the change between them does not reach it.
if [ "$MODE" = stored ]; then
    gzip -t "$REFERENCE_DIR/variants.vcf.gz" 2>/dev/null ||
        step_failed "cannot read $REFERENCE_DIR/variants.vcf.gz"
    stored_version=$(gzip -dc "$REFERENCE_DIR/variants.vcf.gz" 2>/dev/null |
        sed -n '2s/^##source=.* generator_version=\([0-9][0-9]*\).*$/\1/p')
    versions=""
    if [ "$stored_version" != "$GENERATOR_VERSION" ]; then
        if [ -n "$stored_version" ]; then
            versions="; the stored input names generator version $stored_version"
        else
            versions="; the stored input names no generator version"
        fi
        versions="$versions and this run's generator is version ${GENERATOR_VERSION:--}"
    fi
    # grep -a: without it GNU grep prints no line of an input that holds a null byte,
    # and two such inputs compare equal whatever their records.
    if ! gzip -dc "$REFERENCE_DIR/variants.vcf.gz" 2>/dev/null | grep -a -v '^#' |
        cmp -s - <(grep -a -v '^#' "$OUT/variants.vcf"); then
        step_failed "the generated records differ from $REFERENCE_DIR/variants.vcf.gz$versions: the generator does not reproduce the stored input"
    fi
    [ -z "$versions" ] || echo "check_genome: the generated records equal the stored ones$versions"
fi
stage_reached generate

# Ensembl VEP's output: generated in parity mode, read back in stored mode. The
# comparator pairs the files of two directories by name, so each engine's output is
# variants.txt in a directory of its own and no other .txt file is written beside it.
ENSEMBL_OUT="$WORK/ensembl_vep"
VEP_RS_OUT="$WORK/vep_rs"
rm -rf "$ENSEMBL_OUT" "$VEP_RS_OUT"
mkdir -p "$ENSEMBL_OUT" "$VEP_RS_OUT"
case "$MODE" in
parity)
    rm -f "$OUT/ensembl_vep.txt.gz" "$OUT/ensembl_vep.warnings.txt"
    ENSEMBL_VEP_COMMAND="vep --offline --cache --dir_cache /cache --species $GENOME --assembly $ASSEMBLY --cache_version $CACHE_VERSION"
    [ -z "$SET_FLAG" ] || ENSEMBL_VEP_COMMAND="$ENSEMBL_VEP_COMMAND $SET_FLAG"
    [ -z "$REFERENCE_FLAG" ] || ENSEMBL_VEP_COMMAND="$ENSEMBL_VEP_COMMAND $REFERENCE_FLAG"
    ENSEMBL_VEP_COMMAND="$ENSEMBL_VEP_COMMAND -i /input/variants.vcf -o /output/variants.txt --warning_file /output/variants.warnings --no_stats --force_overwrite"
    # Every word of the command line is free of blanks, which the name checks above
    # hold, and pathname expansion is off for the call, so the line splits back into
    # the arguments it records whatever *, ? or [ a name holds.
    set -f
    # shellcheck disable=SC2086
    run_in_image ensembl_vep \
        -v "$ENSEMBL_ROOT:/cache:ro" -v "$OUT:/input:ro" -v "$ENSEMBL_OUT:/output" \
        "$VEP_IMAGE" $ENSEMBL_VEP_COMMAND
    ensembl_status=$?
    set +f
    if [ "$ensembl_status" -ne 0 ] || [ ! -f "$ENSEMBL_OUT/variants.txt" ]; then
        exclude ensembl-vep-failed "Ensembl VEP failed (exit $ensembl_status; $WORK/logs/ensembl_vep.log)"
    fi
    [ "$(count_rows "$ENSEMBL_OUT/variants.txt")" -gt 0 ] ||
        exclude ensembl-vep-failed "Ensembl VEP wrote no row for $RECORDS records"
    gzip -n -c "$ENSEMBL_OUT/variants.txt" >"$OUT/ensembl_vep.txt.gz" ||
        step_failed "cannot write $OUT/ensembl_vep.txt.gz"
    if [ -f "$ENSEMBL_OUT/variants.warnings" ]; then
        cp "$ENSEMBL_OUT/variants.warnings" "$OUT/ensembl_vep.warnings.txt"
    else
        : >"$OUT/ensembl_vep.warnings.txt"
    fi
    ;;
stored)
    gzip -dc "$REFERENCE_DIR/ensembl_vep.txt.gz" >"$ENSEMBL_OUT/variants.txt" ||
        step_failed "cannot read $REFERENCE_DIR/ensembl_vep.txt.gz"
    ;;
esac

# vep-rs.
VEP_RS_COMMAND="vep -i variants.vcf -o variants.txt --offline --json_cache <converted cache> --species $GENOME --assembly $ASSEMBLY"
[ -z "$SET_FLAG" ] || VEP_RS_COMMAND="$VEP_RS_COMMAND $SET_FLAG"
[ -z "$REFERENCE_FLAG" ] || VEP_RS_COMMAND="$VEP_RS_COMMAND $REFERENCE_FLAG"
VEP_RS_COMMAND="$VEP_RS_COMMAND --warning_file <log> --no_stats"
# shellcheck disable=SC2086
run_native vep_rs "$VEP_RS" -i "$OUT/variants.vcf" -o "$VEP_RS_OUT/variants.txt" --offline \
    --json_cache "$JSON_CACHE" --species "$GENOME" --assembly "$ASSEMBLY" $SET_FLAG $REFERENCE_FLAG \
    --warning_file "$WORK/logs/vep_rs.warnings" --no_stats
vep_rs_status=$?
VEP_RS_WROTE=true
if [ "$vep_rs_status" -ne 0 ] || [ ! -f "$VEP_RS_OUT/variants.txt" ] ||
    [ "$(count_rows "$VEP_RS_OUT/variants.txt")" -eq 0 ]; then
    VEP_RS_WROTE=false
    echo "ERROR: [check_genome] vep-rs wrote no row (exit $vep_rs_status; $WORK/logs/vep_rs.log)" >&2
fi

DIGEST="-"
if [ "$VEP_RS_WROTE" = true ]; then
    run_native digest "$TABLE_TOOL" digest "$VEP_RS_OUT/variants.txt" ||
        step_failed "the digest failed ($(tail -n 1 "$WORK/logs/digest.log"))"
    DIGEST=$(tail -n 1 "$WORK/logs/digest.log")
    case "$DIGEST" in
    *[!0-9a-f]* | '') step_failed "the digest is not 64 hex digits: $DIGEST" ;;
    esac
    [ "${#DIGEST}" -eq 64 ] || step_failed "the digest is not 64 hex digits: $DIGEST"
fi

if [ "$MODE" = recheck ]; then
    if [ "$VEP_RS_WROTE" = true ] && [ "$DIGEST" = "$ROW_SHA" ] && [ "$RECORDS" = "$ROW_RECORDS" ]; then
        write_row exact "$DIGEST"
        finish 0 "exact"
    fi
    echo "ERROR: [check_genome] $GENOME_KEY: the row holds $ROW_RECORDS records and digest $ROW_SHA; this run made $RECORDS records and digest $DIGEST" >&2
    finish 1 "recheck-mismatch"
fi

# Compare. The outputs are exact when the comparator leaves no row open, the
# documented divergences set aside. Where the run holds the Ensembl cache, the comparator
# reads it for the names its synonyms file gives to several of its sequences.
if [ "$VEP_RS_WROTE" = false ]; then
    write_row "differs:unclassified" "-"
    finish 1 "differs:unclassified"
fi
rm -rf "$OUT/comparison"
REFERENCE_CACHE=()
[ -z "$ENSEMBL_DIR" ] || REFERENCE_CACHE=(--ensembl-cache "$ENSEMBL_DIR")
set -- --perl-dir "$ENSEMBL_OUT" --rust-dir "$VEP_RS_OUT" \
    --report-dir "$OUT/comparison" --input-vcf "$OUT/variants.vcf" \
    --no-canonical-contigs --assembly "$ASSEMBLY"
# On a RefSeq and on a merged cache Ensembl VEP keeps one transcript a stable id among
# those an input buffer loads, and on a merged cache, whose two transcript sets each number
# their own transcripts, one a number before that. The comparator tells the rows of the
# transcripts dropped so from the converted cache alone.
case "$TRANSCRIPTS" in
refseq | merged) set -- "$@" --json-cache "$JSON_CACHE" ;;
esac
run_native compare "$COMPARATOR" "$@" ${REFERENCE_CACHE[@]+"${REFERENCE_CACHE[@]}"} ||
    step_failed "the comparison failed ($(tail -n 1 "$WORK/logs/compare.log"))"
OPEN_ROWS_FILE="$OUT/comparison/discordant_open.tsv"
if [ ! -f "$OPEN_ROWS_FILE" ] || [ ! -f "$OUT/comparison/summary.json" ]; then
    step_failed "the comparison left no report under $OUT/comparison"
fi
open_rows=$(awk 'NR > 1 && NF { n++ } END { print n + 0 }' "$OPEN_ROWS_FILE")
if [ "$open_rows" -eq 0 ]; then
    write_row exact "$DIGEST"
    finish 0 "exact"
fi
echo "check_genome: $open_rows open rows in $OPEN_ROWS_FILE"
write_row "differs:unclassified" "-"
finish 1 "differs:unclassified"
