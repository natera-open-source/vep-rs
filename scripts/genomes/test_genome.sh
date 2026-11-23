#!/usr/bin/env bash
# Checks genomes of the table under tests/genomes, one at a time: reads the genome's row,
# takes its ready-made Ensembl VEP cache from --cache-root or from Ensembl's public site,
# holds a download to the directory's CHECKSUMS file, and runs check_genome.sh on it.
# docs/genome-checks.md describes the two modes and the line each genome prints.
#
# Standard output carries one line a genome: `ok <key> <mode>`, `FAIL <key> <mode>: <what>`
# or `skip <key>: <why>`. The script runs under bash 3.2.
#
# Run: bash scripts/genomes/test_genome.sh --help
set -uo pipefail
export LC_ALL=C

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
JOB="$HERE/check_genome.sh"
CONVERTER_NAME=vep-storable-to-json
IMAGE_SOURCE="$ROOT/scripts/concordance/generate_reference_output.sh"
ENSEMBL_SITE="https://ftp.ensembl.org"
GENOMES_SITE="https://ftp.ensemblgenomes.ebi.ac.uk"

usage() {
    cat <<'USAGE'
Usage: test_genome.sh --key <genome key>... | --keys-file <file> |
                      --set <set> --sample <n> [--seed <n>]
           [--mode recheck|parity] [--cache-root <dir>] [--max-bytes <n>]
           [--table-dir <dir>] [--vep-image <reference>]
           [--vep-rs <vep>] [--generator <vep-seeded-variants>]
           [--comparator <vep-outputs-compare>] [--table-tool <vep-genome-table>]
           [--work <dir>]

The genomes, one of:
  --key <genome key>      a key of the table, <set>__<genome>__<assembly>__<transcripts>;
                          may be given more than once
  --keys-file <file>      one key a line; blank lines and lines starting with # are left out
  --set <set> --sample <n> [--seed <n>]
                          n keys of the set, picked by the seed (default 1); the same
                          table, set, n and seed pick the same keys

Modes:
  recheck  (default) runs vep-rs alone on the genome's generated records and holds its
           digest and record count to the row
  parity   runs Ensembl VEP and vep-rs and holds the outcome to the row's; a row that
           says differs: stands when the two outputs still differ

The cache:
  --cache-root <dir>   the cache archives are files of this directory, under the names
                       Ensembl ships them; nothing is downloaded
  without it           each archive is downloaded from Ensembl's public site and held to
                       the length the site states and to the CHECKSUMS file of its
                       directory; an archive that file does not list, and one whose
                       directory has no such file, as every human pangenome archive,
                       is held to gzip's own check
  --max-bytes <n>      a genome whose archive is larger is skipped by name

--table-dir is the table directory (default tests/genomes). --vep-image is the Ensembl VEP
image (default: the reference scripts/concordance/generate_reference_output.sh pins).
The four binaries default to the release builds of this checkout: target/release/vep and
tools/<name>/target/release/<name> for vep-seeded-variants, vep-outputs-compare and
vep-genome-table.

The conversion runs the release build of tools/vep-storable-to-json that check_genome.sh
finds; its --help says where. recheck converts in the image when docker holds it, and on
the host otherwise; parity needs the image, and so does a run given --vep-image. This
script builds nothing; the build command is
  cargo build --release --locked --manifest-path tools/vep-storable-to-json/Cargo.toml

--work keeps each genome's files under <dir>/<genome key>/: out/ holds row.tsv, job.json,
the input records and the comparison report, and job.log the job's own lines. A genome
that passes leaves out/ and job.log alone. Without --work a new directory is used, and
removed when no genome fails.

A row whose outcome is excluded:<name> is skipped by that name and never run.

Exit codes:
  0  every genome that ran is ok
  1  a genome failed
  2  a usage error, or no genome ran: every one was skipped, or the host has no docker,
     no image, or no build of the converter
USAGE
}

usage_error() {
    echo "ERROR: [test_genome] $1" >&2
    echo "Run with --help for the usage." >&2
    exit 2
}

# Digits without a leading zero, the spelling the table takes.
is_whole_number() {
    case "$1" in
    '' | *[!0-9]* | 0?*) return 1 ;;
    *) return 0 ;;
    esac
}

MODE=recheck
KEYS=()
KEYS_FILE=""
SET_NAME=""
SAMPLE=""
PICK_SEED=""
CACHE_ROOT=""
MAX_BYTES=""
TABLE_DIR="$ROOT/tests/genomes"
VEP_IMAGE=""
VEP_RS="$ROOT/target/release/vep"
GENERATOR="$ROOT/tools/vep-seeded-variants/target/release/vep-seeded-variants"
COMPARATOR="$ROOT/tools/vep-outputs-compare/target/release/vep-outputs-compare"
TABLE_TOOL="$ROOT/tools/vep-genome-table/target/release/vep-genome-table"
WORK=""

while [ $# -gt 0 ]; do
    case "$1" in
    --help | -h)
        usage
        exit 0
        ;;
    --key | --keys-file | --set | --sample | --seed | --mode | --cache-root | --max-bytes | \
        --table-dir | --vep-image | --vep-rs | --generator | --comparator | --table-tool | --work)
        [ $# -ge 2 ] || usage_error "$1 takes a value"
        case "$1" in
        --key) KEYS+=("$2") ;;
        --keys-file) KEYS_FILE=$2 ;;
        --set) SET_NAME=$2 ;;
        --sample) SAMPLE=$2 ;;
        --seed) PICK_SEED=$2 ;;
        --mode) MODE=$2 ;;
        --cache-root) CACHE_ROOT=$2 ;;
        --max-bytes) MAX_BYTES=$2 ;;
        --table-dir) TABLE_DIR=$2 ;;
        --vep-image) VEP_IMAGE=$2 ;;
        --vep-rs) VEP_RS=$2 ;;
        --generator) GENERATOR=$2 ;;
        --comparator) COMPARATOR=$2 ;;
        --table-tool) TABLE_TOOL=$2 ;;
        --work) WORK=$2 ;;
        esac
        shift 2
        ;;
    *) usage_error "unknown argument $1" ;;
    esac
done

case "$MODE" in
recheck | parity) ;;
*) usage_error "--mode is recheck or parity, not $MODE" ;;
esac
[ -d "$TABLE_DIR" ] || usage_error "no directory $TABLE_DIR"
[ -f "$TABLE_DIR/sets.tsv" ] || usage_error "no sets.tsv under $TABLE_DIR"
[ -z "$CACHE_ROOT" ] || [ -d "$CACHE_ROOT" ] || usage_error "no directory $CACHE_ROOT"
[ -z "$MAX_BYTES" ] || is_whole_number "$MAX_BYTES" || usage_error "--max-bytes is a whole number of bytes"
[ -z "$KEYS_FILE" ] || [ -f "$KEYS_FILE" ] || usage_error "no file $KEYS_FILE"
[ -f "$JOB" ] || usage_error "no check_genome.sh beside this script"

# need_binary <option> <path>
need_binary() {
    if [ ! -f "$2" ] || [ ! -x "$2" ]; then
        usage_error "$1: $2 is not an executable file; build it or name another"
    fi
}
need_binary --table-tool "$TABLE_TOOL"
need_binary --vep-rs "$VEP_RS"
need_binary --generator "$GENERATOR"
[ "$MODE" != parity ] || need_binary --comparator "$COMPARATOR"
if [ -z "$CACHE_ROOT" ]; then
    command -v curl >/dev/null 2>&1 || usage_error "curl is not on PATH; name the archives' directory with --cache-root"
fi

# The genomes. --set with --sample picks among the set's keys in the table's order.
# pick_sample <n> <seed>: the keys on standard input, n of them, in key order. The
# generator is the Lehmer one of modulus 2^31 - 1: the seed is reduced digit by digit and
# every product stays exact in the floating point every awk computes in, so a seed of
# any length picks the same keys on every host.
pick_sample() {
    awk -v want="$1" -v seed="$2" '
        { key[NR] = $0 }
        END {
            state = 0
            for (i = 1; i <= length(seed); i++) state = (state * 10 + substr(seed, i, 1)) % 2147483646
            state += 1
            for (i = 0; i < 3; i++) state = (state * 48271) % 2147483647
            if (want > NR) want = NR
            for (i = 1; i <= want; i++) {
                state = (state * 48271) % 2147483647
                j = i + state % (NR - i + 1)
                held = key[i]; key[i] = key[j]; key[j] = held
                print key[i]
            }
        }' | sort
}

if [ -n "$SET_NAME" ]; then
    if [ "${#KEYS[@]}" -gt 0 ] || [ -n "$KEYS_FILE" ]; then
        usage_error "--set accompanies neither --key nor --keys-file; name the genomes one way"
    fi
    is_whole_number "$SAMPLE" && [ "$SAMPLE" -gt 0 ] || usage_error "--set needs --sample <n>, a whole number above zero"
    [ -n "$PICK_SEED" ] || PICK_SEED=1
    is_whole_number "$PICK_SEED" || usage_error "--seed is a whole number without a leading zero"
    set_keys=$("$TABLE_TOOL" keys "$TABLE_DIR" --set "$SET_NAME") ||
        usage_error "the table under $TABLE_DIR does not list its keys"
    [ -n "$set_keys" ] || usage_error "the table under $TABLE_DIR holds no row of the set $SET_NAME"
    while IFS= read -r key; do
        KEYS+=("$key")
    done <<<"$(printf '%s\n' "$set_keys" | pick_sample "$SAMPLE" "$PICK_SEED")"
else
    [ -z "$SAMPLE" ] && [ -z "$PICK_SEED" ] || usage_error "--sample and --seed accompany --set"
    if [ -n "$KEYS_FILE" ]; then
        while IFS= read -r key || [ -n "$key" ]; do
            case "$key" in
            '' | '#'*) ;;
            *) KEYS+=("$key") ;;
            esac
        done <"$KEYS_FILE"
    fi
fi
[ "${#KEYS[@]}" -gt 0 ] || usage_error "no genome is named: give --key, --keys-file or --set with --sample"

# What the host lacks for the mode, found once and before anything is downloaded. Parity
# runs Ensembl VEP in the image. Recheck needs the conversion alone, which the job runs
# in the image when docker holds it and on the host otherwise. Either way the job runs
# the release build of the converter, looked for here where the job looks for it: under
# CARGO_TARGET_DIR when that holds one, else under the tool's own target directory.
CONVERTER="$ROOT/tools/$CONVERTER_NAME/target/release/$CONVERTER_NAME"
if [ -n "${CARGO_TARGET_DIR:-}" ] && [ -f "$CARGO_TARGET_DIR/release/$CONVERTER_NAME" ] &&
    [ -x "$CARGO_TARGET_DIR/release/$CONVERTER_NAME" ]; then
    CONVERTER="$CARGO_TARGET_DIR/release/$CONVERTER_NAME"
fi
NO_CONVERTER=""
if [ ! -f "$CONVERTER" ] || [ ! -x "$CONVERTER" ]; then
    NO_CONVERTER="no build of the converter at $CONVERTER: build it with 'cargo build --release --locked --manifest-path $ROOT/tools/$CONVERTER_NAME/Cargo.toml'"
fi
pinned_image() {
    sed -n 's/^PERL_IMAGE="\([^"]*\)"$/\1/p' "$IMAGE_SOURCE" 2>/dev/null | head -n 1
}
IMAGE_GIVEN=true
if [ -z "$VEP_IMAGE" ]; then
    IMAGE_GIVEN=false
    VEP_IMAGE=$(pinned_image)
    [ -n "$VEP_IMAGE" ] || usage_error "$IMAGE_SOURCE pins no image; name one with --vep-image"
fi
has_image() {
    command -v docker >/dev/null 2>&1 && docker image inspect "$VEP_IMAGE" >/dev/null 2>&1
}
JOB_IMAGE=""
HOST_LACKS=""
if has_image; then
    JOB_IMAGE=$VEP_IMAGE
    HOST_LACKS=$NO_CONVERTER
elif [ "$MODE" = parity ] || [ "$IMAGE_GIVEN" = true ]; then
    if command -v docker >/dev/null 2>&1; then
        HOST_LACKS="no image $VEP_IMAGE in docker"
    else
        HOST_LACKS="no docker on PATH"
    fi
    [ -z "$NO_CONVERTER" ] || HOST_LACKS="$HOST_LACKS and $NO_CONVERTER"
elif [ -n "$NO_CONVERTER" ]; then
    HOST_LACKS="no image $VEP_IMAGE in docker and $NO_CONVERTER"
fi

WORK_GIVEN=true
if [ -z "$WORK" ]; then
    WORK_GIVEN=false
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/test_genome.XXXXXX") || usage_error "cannot create a work directory"
fi
mkdir -p "$WORK" || usage_error "cannot create $WORK"
WORK=$(cd "$WORK" && pwd)

OKS=0
FAILS=0
SKIPS=0
passed() {
    OKS=$((OKS + 1))
    echo "ok $1 $MODE"
}
failed() {
    FAILS=$((FAILS + 1))
    echo "FAIL $1 $MODE: $2"
}
skipped() {
    SKIPS=$((SKIPS + 1))
    echo "skip $1: $2"
}

# remote_bytes <url>: the length the site states for the file, or a failure when the
# site does not answer for it.
remote_bytes() {
    curl -fsSIL --connect-timeout 30 --retry 3 "$1" 2>/dev/null | tr -d '\r' |
        awk 'tolower($1) == "content-length:" { bytes = $2 } END { if (bytes == "") exit 1; print bytes }'
}

# sum_of <file>: the checksum and the block count a CHECKSUMS line carries, as numbers.
sum_of() {
    sum "$1" | awk '{ print $1 + 0, $2 + 0 }'
}

# last_error <log>: the job's last error line, without its component.
last_error() {
    grep -F 'ERROR: [check_genome] ' "$1" 2>/dev/null | tail -n 1 | sed 's/^ERROR: \[check_genome\] //'
}

# check_one <genome key>: prints the genome's one line.
check_one() {
    local key=$1 dir="$WORK/$1"
    local row genome assembly transcripts set_name cache_version seed outcome
    local release accession="" suffix="" names url name tarball="" bytes="" listed unlisted answer found status
    local job_outcome what

    # The key names the genome's directory under the work directory and is an argument of
    # the table tool, so a key that is empty, a path or an option fails with no file touched.
    case "$key" in
    '' | . | .. | */* | -*)
        failed "$key" "the key is not <set>__<genome>__<assembly>__<transcripts>: it is empty, a path or an option"
        return
        ;;
    esac

    rm -rf "$dir"
    mkdir -p "$dir/download" || {
        failed "$key" "cannot create $dir"
        return
    }
    row=$("$TABLE_TOOL" row "$TABLE_DIR" "$key" 2>"$dir/row.log") || {
        failed "$key" "$(tail -n 1 "$dir/row.log" | sed 's/^ERROR: \[[^]]*\] //')"
        return
    }
    IFS=$'\t' read -r genome assembly transcripts set_name cache_version seed _ outcome _ \
        <<<"$(printf '%s\n' "$row" | sed -n 2p)"
    case "$outcome" in
    excluded:*)
        skipped "$key" "the row is $outcome"
        return
        ;;
    exact) ;;
    *)
        if [ "$MODE" = recheck ]; then
            skipped "$key" "the row is $outcome, so there is no digest to recheck"
            return
        fi
        ;;
    esac
    if [ -n "$HOST_LACKS" ]; then
        skipped "$key" "$HOST_LACKS"
        return
    fi

    # The archive's name and the directory Ensembl serves it from.
    release=$(awk -F '\t' -v set="$set_name" 'NR > 1 && $1 == set { print $2 }' "$TABLE_DIR/sets.tsv" | head -n 1)
    if [ -z "$release" ]; then
        failed "$key" "sets.tsv has no line for the set $set_name"
        return
    fi
    case "$transcripts" in
    refseq) suffix="_refseq" ;;
    merged) suffix="_merged" ;;
    esac
    names="${genome}${suffix}_vep_${cache_version}_${assembly}.tar.gz"
    case "$set_name" in
    vertebrates) url="$ENSEMBL_SITE/pub/release-$release/variation/indexed_vep_cache" ;;
    metazoa | plants | fungi | protists)
        url="$GENOMES_SITE/pub/$set_name/release-$release/variation/indexed_vep_cache"
        ;;
    bacteria_*_collection)
        url="$GENOMES_SITE/pub/bacteria/release-$release/variation/indexed_vep_cache/$set_name"
        ;;
    fungi_*_collection)
        url="$GENOMES_SITE/pub/fungi/release-$release/variation/indexed_vep_cache/$set_name"
        ;;
    protists_*_collection)
        url="$GENOMES_SITE/pub/protists/release-$release/variation/indexed_vep_cache/$set_name"
        ;;
    pangenome)
        # The genome name ends in the assembly accession, written gca<digits>v<version>.
        # The archive is named by the accession and the dated directory, the date with
        # an underscore or with a hyphen.
        accession=$(printf '%s\n' "$genome" | sed -n 's/^.*_gca\([0-9][0-9]*\)v\([0-9][0-9]*\)$/GCA_\1.\2/p')
        if [ -z "$accession" ]; then
            failed "$key" "the genome name $genome ends in no assembly accession"
            return
        fi
        url="$ENSEMBL_SITE/pub/rapid-release/species/Homo_sapiens/$accession/ensembl/variation/$release/indexed_vep_cache"
        names="Homo_sapiens-$accession-$release.tar.gz"
        case "$release" in
        *_*) names="$names Homo_sapiens-$accession-$(printf '%s' "$release" | tr _ -).tar.gz" ;;
        esac
        ;;
    *)
        failed "$key" "no address is known for the set $set_name"
        return
        ;;
    esac

    if [ -n "$CACHE_ROOT" ]; then
        for name in $names; do
            if [ -f "$CACHE_ROOT/$name" ]; then
                tarball="$CACHE_ROOT/$name"
                break
            fi
        done
        if [ -z "$tarball" ]; then
            failed "$key" "no ${names// / or } under $CACHE_ROOT"
            return
        fi
        bytes=$(wc -c <"$tarball" | tr -d ' ')
    else
        for name in $names; do
            if bytes=$(remote_bytes "$url/$name"); then
                tarball="$dir/download/$name"
                break
            fi
        done
        if [ -z "$tarball" ]; then
            failed "$key" "$url does not answer for ${names// / or }"
            return
        fi
    fi
    name=${tarball##*/}
    if [ -n "$MAX_BYTES" ] && is_whole_number "$bytes" && [ "$bytes" -gt "$MAX_BYTES" ]; then
        skipped "$key" "$name is $bytes bytes, over --max-bytes $MAX_BYTES"
        return
    fi

    if [ -z "$CACHE_ROOT" ]; then
        if ! curl -fsSL --connect-timeout 30 --retry 3 -o "$tarball" "$url/$name" 2>"$dir/download.log"; then
            failed "$key" "cannot download $url/$name ($(tail -n 1 "$dir/download.log"))"
            return
        fi
        found=$(wc -c <"$tarball" | tr -d ' ')
        if [ "$found" != "$bytes" ]; then
            failed "$key" "$name is $found bytes where the site states $bytes"
            return
        fi
        # An archive that the directory's CHECKSUMS file does not list, or whose directory
        # has no such file, is held to gzip's own check. curl exits 22 on every HTTP error
        # status, so the status 404 alone says the directory has no such file.
        listed=""
        unlisted="no CHECKSUMS file beside $name"
        if [ "$set_name" != pangenome ]; then
            answer=$(curl -fsSL --connect-timeout 30 --retry 3 -o "$dir/download/CHECKSUMS" \
                -w '%{http_code}' "$url/CHECKSUMS" 2>>"$dir/download.log")
            case "$?:$answer" in
            0:*)
                listed=$(awk -v name="$name" '$3 == name { print $1 + 0, $2 + 0 }' "$dir/download/CHECKSUMS" | head -n 1)
                unlisted="no CHECKSUMS line for $name"
                ;;
            22:404) ;;
            *)
                failed "$key" "cannot read $url/CHECKSUMS, so $name is not verified ($(tail -n 1 "$dir/download.log"))"
                return
                ;;
            esac
        fi
        if [ -n "$listed" ]; then
            found=$(sum_of "$tarball")
            if [ "$found" != "$listed" ]; then
                failed "$key" "$name sums to $found where CHECKSUMS says $listed"
                return
            fi
        elif gzip -t "$tarball" 2>/dev/null; then
            echo "test_genome: $key: $unlisted; held to the stated length and gzip's check" >&2
        else
            failed "$key" "$name is not a whole gzip archive"
            return
        fi
    fi

    # The job, with the row's seed.
    local job_args
    job_args=(--mode "$MODE" --cache-tarball "$tarball"
        --genome "$genome" --assembly "$assembly" --transcripts "$transcripts"
        --set "$set_name" --cache-version "$cache_version" --seed "$seed"
        --vep-rs "$VEP_RS" --generator "$GENERATOR" --table-tool "$TABLE_TOOL"
        --work "$dir/work" --out "$dir/out")
    [ -z "$JOB_IMAGE" ] || job_args+=(--vep-image "$JOB_IMAGE")
    case "$MODE" in
    recheck) job_args+=(--table-dir "$TABLE_DIR") ;;
    parity) job_args+=(--comparator "$COMPARATOR") ;;
    esac
    bash "$JOB" "${job_args[@]}" >"$dir/job.log" 2>&1
    status=$?

    what=$(last_error "$dir/job.log")
    if [ "$status" -eq 2 ]; then
        case "$what" in
        skip:*) skipped "$key" "${what#skip: }" ;;
        *) failed "$key" "the job refused its command line: $what" ;;
        esac
        return
    fi
    if [ "$MODE" = recheck ]; then
        case "$status" in
        0) ;;
        1)
            failed "$key" "${what#"$key": }"
            return
            ;;
        *)
            failed "$key" "the job ended with exit $status: $what"
            return
            ;;
        esac
    else
        job_outcome=""
        if [ -f "$dir/out/row.tsv" ]; then
            job_outcome=$(awk -F '\t' 'NR == 2 { print $8 }' "$dir/out/row.tsv")
        fi
        if [ -z "$job_outcome" ]; then
            failed "$key" "the job ended with exit $status and wrote no row: $what"
            return
        fi
        # The job names no difference: every row that says differs: stands when the
        # outputs differ.
        case "$outcome" in
        differs:*) [ "${job_outcome%%:*}" = differs ] || job_outcome="other:$job_outcome" ;;
        *) [ "$job_outcome" = "$outcome" ] || job_outcome="other:$job_outcome" ;;
        esac
        case "$job_outcome" in
        other:*)
            failed "$key" "the row says $outcome; this run gave ${job_outcome#other:}"
            return
            ;;
        esac
    fi
    rm -rf "$dir/work" "$dir/download" "$dir/row.log" "$dir/download.log"
    passed "$key"
}

for key in "${KEYS[@]}"; do
    check_one "$key"
done

echo "test_genome: $OKS ok, $FAILS failed, $SKIPS skipped" >&2
if [ "$FAILS" -gt 0 ]; then
    echo "test_genome: each failed genome's files are under $WORK" >&2
    exit 1
fi
[ "$WORK_GIVEN" = true ] || rm -rf "$WORK"
if [ "$OKS" -eq 0 ]; then
    echo "ERROR: [test_genome] skip: no genome ran" >&2
    exit 2
fi
exit 0
