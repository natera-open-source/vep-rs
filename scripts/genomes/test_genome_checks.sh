#!/usr/bin/env bash
# Tests for test_genome.sh, and for .github/workflows/genome-check.yml read as text: the
# address formed for a row of each kind of set, the refusal of a download that fails its
# check, each mode's verdicts, the skips, where the converter's build is looked for, the
# exit codes, the seeded pick, and the rules the workflow is held to.
#
# Bash, tar and gzip over a made-up table and made-up cache archives. No network, no
# Docker, no Ensembl VEP, no vep-rs build: a stub `curl` and a stub `docker` stand first
# on PATH, a stub build of the converter stands under a made-up CARGO_TARGET_DIR, and
# stub engines and tools are passed by path. check_genome.sh itself runs under
# test_genome.sh.
#
# Run: bash scripts/genomes/test_genome_checks.sh
# shellcheck disable=SC2016  # the workflow's text is matched as written, unexpanded
set -uo pipefail
export LC_ALL=C

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
same() { # same <name> <expected> <found>
    if [ "$2" = "$3" ]; then ok "$1"; else bad "$1: expected [$2], found [$3]"; fi
}
has() { # has <name> <text> <where>
    case "$3" in
    *"$2"*) ok "$1" ;;
    *) bad "$1: [$2] is not in [$3]" ;;
    esac
}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
SCRIPT="$HERE/test_genome.sh"
WORKFLOW="$ROOT/.github/workflows/genome-check.yml"
CI="$ROOT/.github/workflows/ci.yml"
IMAGE_SOURCE="$ROOT/scripts/concordance/generate_reference_output.sh"

for tool in tar gzip sum awk; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: [test_genome_checks] $tool is not on PATH" >&2
        exit 2
    }
done
for file in "$SCRIPT" "$HERE/check_genome.sh" "$WORKFLOW" "$CI" "$IMAGE_SOURCE"; do
    [ -f "$file" ] || {
        bad "no file $file"
        echo "passed: $PASS   failed: $FAIL"
        exit 1
    }
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
BIN="$TMP/bin"
TOOLS="$TMP/tools"
CTL="$TMP/ctl"
SITE="$TMP/site"
TABLE="$TMP/table"
STUB_TARGET="$TMP/cargo-target"
mkdir -p "$BIN" "$TOOLS" "$CTL" "$SITE" "$TABLE" "$STUB_TARGET/release"
export CTL SITE

# The stub curl serves $SITE/<host>/<path> for https://<host>/<path>. A request for the
# headers states the file's length, or the one ctl/stated_bytes holds. A request for a
# CHECKSUMS file fails with the exit code ctl/checksums_exit holds, and a request with
# -w prints the HTTP status, the one ctl/checksums_status holds for such a failure. Each
# request is a line of ctl/curl.log.
cat >"$BIN/curl" <<'STUB'
#!/usr/bin/env bash
headers=false output="" url="" stated=false
state() { [ "$stated" = false ] || printf '%s' "$1"; }
while [ $# -gt 0 ]; do
    case "$1" in
    -o) output=$2; shift ;;
    -w) stated=true; shift ;;
    --connect-timeout | --retry) shift ;;
    -*I*) headers=true ;;
    -*) ;;
    *) url=$1 ;;
    esac
    shift
done
file="$SITE/${url#https://}"
if [ "$headers" = true ]; then
    echo "HEAD $url" >>"$CTL/curl.log"
    [ -f "$file" ] || { echo "curl: (22) the server answered 404" >&2; exit 22; }
    bytes=$(wc -c <"$file" | tr -d ' ')
    [ ! -f "$CTL/stated_bytes" ] || bytes=$(cat "$CTL/stated_bytes")
    printf 'HTTP/1.1 200 OK\r\nContent-Type: application/gzip\r\nContent-Length: %s\r\n\r\n' "$bytes"
    exit 0
fi
echo "GET $url" >>"$CTL/curl.log"
case "$url" in
*/CHECKSUMS)
    if [ -f "$CTL/checksums_exit" ]; then
        echo "curl: a made-up failure" >&2
        state "$(cat "$CTL/checksums_status")"
        exit "$(cat "$CTL/checksums_exit")"
    fi
    ;;
esac
[ -f "$file" ] || { echo "curl: (22) the server answered 404" >&2; state 404; exit 22; }
cp "$file" "$output"
state 200
STUB

# The stub docker: `image inspect` answers for the image unless ctl/no_image stands;
# `run` plays the converter or Ensembl VEP by what its arguments name, writing into the
# host directory mounted at /output. Each call is a line of ctl/docker.log.
cat >"$BIN/docker" <<'STUB'
#!/usr/bin/env bash
echo "$*" >>"$CTL/docker.log"
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
output=""
previous=""
for word in "$@"; do
    if [ "$previous" = -v ]; then
        case "$word" in *:/output) output=${word%:/output} ;; esac
    fi
    previous=$word
done
case "$*" in
*vep-storable-to-json*)
    echo '{}' >"$output/info.json"
    echo "Done. Converted: 1 / 1 files, 5 total transcripts, 0 variation files copied, 0 regulatory files"
    ;;
*' vep '*)
    [ ! -f "$CTL/ensembl_fails" ] || { echo "a made-up failure" >&2; exit 2; }
    printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n%s\n' "$(printf 'r1\t1:100\tG')" >"$output/variants.txt"
    ;;
esac
exit 0
STUB

# The stub build of the converter, for a conversion without an image: it appends its
# arguments to the log STUB_CONVERTER_LOG names, ctl/converter.log by default, and plays
# the converter into its second argument. The copy under the made-up CARGO_TARGET_DIR is
# the build every run finds unless a check names another; the stub docker plays the
# converter itself, so a conversion in the image mounts this file and never starts it.
cat >"$TMP/converter.stub" <<'STUB'
#!/usr/bin/env bash
echo "$*" >>"$CTL/${STUB_CONVERTER_LOG:-converter.log}"
echo '{}' >"$2/info.json"
echo "Done. Converted: 1 / 1 files, 5 total transcripts, 0 variation files copied, 0 regulatory files"
STUB
cp "$TMP/converter.stub" "$STUB_TARGET/release/vep-storable-to-json"
export CARGO_TARGET_DIR="$STUB_TARGET"

# The stub generator writes three records, or as many as ctl/records holds.
cat >"$TOOLS/vep-seeded-variants" <<'STUB'
#!/usr/bin/env bash
if [ "$1" = --version ]; then echo "vep-seeded-variants 0.0.0 (generator version 1)"; exit 0; fi
echo "$*" >>"$CTL/generator.log"
output="" summary=""
while [ $# -gt 0 ]; do
    case "$1" in
    --output) output=$2; shift ;;
    --summary) summary=$2; shift ;;
    esac
    shift
done
records=3
[ ! -f "$CTL/records" ] || records=$(cat "$CTL/records")
{
    printf '##fileformat=VCFv4.2\n##source=vep-seeded-variants generator_version=1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    i=0
    while [ "$i" -lt "$records" ]; do
        i=$((i + 1))
        printf '1\t%s00\trecord_%s\tA\tG\t.\t.\t.\n' "$i" "$i"
    done
} >"$output"
printf 'category\twanted\twritten\tpassed_over\treason\n' >"$summary"
STUB

cat >"$TOOLS/vep" <<'STUB'
#!/usr/bin/env bash
if [ "$1" = --version ]; then echo "vep 0.0.0"; exit 0; fi
echo "$*" >>"$CTL/vep_rs.log"
output=""
while [ $# -gt 0 ]; do
    if [ "$1" = -o ]; then output=$2; fi
    shift
done
printf '## made-up header\n#Uploaded_variation\tLocation\tAllele\n%s\n' "$(printf 'r1\t1:100\tG')" >"$output"
STUB

# The stub comparator leaves no open row, or two while ctl/outputs_differ stands.
cat >"$TOOLS/vep-outputs-compare" <<'STUB'
#!/usr/bin/env bash
if [ "$1" = --version ]; then exit 2; fi
echo "$*" >>"$CTL/comparator.log"
report=""
while [ $# -gt 0 ]; do
    if [ "$1" = --report-dir ]; then report=$2; fi
    shift
done
mkdir -p "$report"
echo '{}' >"$report/summary.json"
printf 'file\tside\tkey\n' >"$report/discordant_open.tsv"
if [ -f "$CTL/outputs_differ" ]; then
    printf 'variants.txt\tmissing\tr1\nvariants.txt\textra\tr9\n' >>"$report/discordant_open.tsv"
fi
STUB

# The stub table tool reads the two row files of the directory it is given: `row` prints
# the header and the row of a key, or one field of it; `keys` prints the keys, of one
# set with --set; `digest` prints the digest every exact row of the made-up table holds.
# Each call is a line of ctl/table.log.
cat >"$TOOLS/vep-genome-table" <<'STUB'
#!/usr/bin/env bash
echo "$*" >>"$CTL/table.log"
header=$(printf 'genome\tassembly\ttranscripts\tset\tcache_version\tseed\trecords\toutcome\tsha256')
rows() {
    tail -n +2 "$1/genomes.tsv"
    gzip -dc "$1/bacteria.tsv.gz" | tail -n +2
}
case "$1" in
digest) printf '%s\n' aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa ;;
row)
    line=$(rows "$2" | awk -F '\t' -v key="$3" '$4 "__" $1 "__" $2 "__" $3 == key')
    [ -n "$line" ] || { echo "ERROR: [vep-genome-table] $2: no row has the key $3" >&2; exit 1; }
    if [ "${4:-}" = --field ]; then
        printf '%s\n%s\n' "$header" "$line" | awk -F '\t' -v field="$5" '
            NR == 1 { for (i = 1; i <= NF; i++) if ($i == field) column = i; next }
            { print $column }'
    else
        printf '%s\n%s\n' "$header" "$line"
    fi
    ;;
keys)
    set_name=""
    [ "${3:-}" != --set ] || set_name=$4
    rows "$2" | awk -F '\t' -v set="$set_name" 'set == "" || $4 == set { print $4 "__" $1 "__" $2 "__" $3 }' | sort
    ;;
*) echo "ERROR: [vep-genome-table] no subcommand $1" >&2; exit 2 ;;
esac
STUB
chmod +x "$BIN"/* "$TOOLS"/* "$TMP/converter.stub" "$STUB_TARGET/release/vep-storable-to-json"
STUB_PATH="$BIN:$PATH"

# The made-up table. Every exact row holds the stub tool's digest and three records,
# apart from the two rows named for what they hold in their place.
DIGEST=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
OTHER_DIGEST=bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
HEADER=$(printf 'genome\tassembly\ttranscripts\tset\tcache_version\tseed\trecords\toutcome\tsha256')
row() { # row <genome> <assembly> <transcripts> <set> <cache version> <records> <outcome> <sha256>
    printf '%s\t%s\t%s\t%s\t%s\t7\t%s\t%s\t%s\n' "$@"
}
{
    echo "$HEADER"
    row example_fungus_gca_000000002 ASM4v1 ensembl fungi_ascomycota1_collection 63 3 exact "$DIGEST"
    row example_differs ASM6v1 ensembl metazoa 63 3 differs:unclassified -
    row example_excluded ASM7v1 ensembl metazoa 63 0 excluded:no-transcripts -
    row homo_sapiens_gca000000004v2 EXAMPLE.pri ensembl pangenome 107 3 exact "$DIGEST"
    row homo_sapiens_gca000000005v1 EXAMPLE.alt ensembl pangenome 107 3 exact "$DIGEST"
    row example_plant ASM2v1 ensembl plants 63 3 exact "$DIGEST"
    row example_protist_gca_000000003 ASM5v1 ensembl protists_alveolata1_collection 63 3 exact "$DIGEST"
    row example_absent ASM13v1 ensembl vertebrates 116 3 exact "$DIGEST"
    row example_badsum ASM11v1 ensembl vertebrates 116 3 exact "$DIGEST"
    row example_miscount ASM10v1 ensembl vertebrates 116 4 exact "$DIGEST"
    row example_mismatch ASM9v1 ensembl vertebrates 116 3 exact "$OTHER_DIGEST"
    row example_torn ASM12v1 ensembl vertebrates 116 3 exact "$DIGEST"
    row example_unlisted ASM8v1 ensembl vertebrates 116 3 exact "$DIGEST"
    row example_unplaced ASM14v1 ensembl fungi 63 3 exact "$DIGEST"
    row example_vertebrate ASM1v1 ensembl vertebrates 116 3 exact "$DIGEST"
    row example_vertebrate ASM1v1 merged vertebrates 116 3 exact "$DIGEST"
    row example_vertebrate ASM1v1 refseq vertebrates 116 3 exact "$DIGEST"
} >"$TABLE/genomes.tsv"
{
    echo "$HEADER"
    row example_bacterium_gca_000000001 ASM3v1 ensembl bacteria_7_collection 63 3 exact "$DIGEST"
} | gzip -c >"$TABLE/bacteria.tsv.gz"
printf 'outcome\tmeaning\ndiffers:unclassified\tThe outputs differ.\nexcluded:no-transcripts\tNo transcript.\n' >"$TABLE/outcomes.tsv"
{
    printf 'set\trelease\tcache_version\tcaches\n'
    printf 'vertebrates\t116\t116\t9\n'
    printf 'metazoa\t63\t63\t2\n'
    printf 'plants\t63\t63\t1\n'
    printf 'pangenome\t2022_10\t107\t2\n'
    printf 'bacteria_7_collection\t63\t63\t1\n'
    printf 'fungi_ascomycota1_collection\t63\t63\t1\n'
    printf 'protists_alveolata1_collection\t63\t63\t1\n'
} >"$TABLE/sets.tsv"

VERTEBRATE=vertebrates__example_vertebrate__ASM1v1__ensembl
PLANT=plants__example_plant__ASM2v1__ensembl
DIFFERS=metazoa__example_differs__ASM6v1__ensembl
EXCLUDED=metazoa__example_excluded__ASM7v1__ensembl
MISMATCH=vertebrates__example_mismatch__ASM9v1__ensembl

# The made-up site. serve <directory> <archive> <path inside the archive> [unlisted]: an
# archive that holds the cache directory, and its line in the directory's CHECKSUMS.
ENSEMBL="ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache"
GENOMES="ftp.ensemblgenomes.ebi.ac.uk/pub"
PANGENOME="ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens"
serve() {
    local dir="$SITE/$1" stage="$TMP/stage"
    rm -rf "$stage"
    mkdir -p "$stage/$3" "$dir"
    printf 'species\tmade-up\n' >"$stage/$3/info.txt"
    tar -czf "$dir/$2" -C "$stage" "${3%%/*}"
    [ "${4:-}" = unlisted ] || printf '%s %s\n' "$(sum "$dir/$2" | awk '{ print $1, $2 }')" "$2" >>"$dir/CHECKSUMS"
}
serve "$ENSEMBL" example_vertebrate_vep_116_ASM1v1.tar.gz example_vertebrate/116_ASM1v1
serve "$ENSEMBL" example_vertebrate_refseq_vep_116_ASM1v1.tar.gz example_vertebrate_refseq/116_ASM1v1
serve "$ENSEMBL" example_vertebrate_merged_vep_116_ASM1v1.tar.gz example_vertebrate_merged/116_ASM1v1
serve "$ENSEMBL" example_mismatch_vep_116_ASM9v1.tar.gz example_mismatch/116_ASM9v1
serve "$ENSEMBL" example_miscount_vep_116_ASM10v1.tar.gz example_miscount/116_ASM10v1
serve "$ENSEMBL" example_unlisted_vep_116_ASM8v1.tar.gz example_unlisted/116_ASM8v1 unlisted
serve "$ENSEMBL" example_badsum_vep_116_ASM11v1.tar.gz example_badsum/116_ASM11v1 unlisted
printf '1 1 example_badsum_vep_116_ASM11v1.tar.gz\n' >>"$SITE/$ENSEMBL/CHECKSUMS"
printf 'not an archive\n' >"$SITE/$ENSEMBL/example_torn_vep_116_ASM12v1.tar.gz"
serve "$GENOMES/plants/release-63/variation/indexed_vep_cache" example_plant_vep_63_ASM2v1.tar.gz example_plant/63_ASM2v1
serve "$GENOMES/metazoa/release-63/variation/indexed_vep_cache" example_differs_vep_63_ASM6v1.tar.gz example_differs/63_ASM6v1
serve "$GENOMES/bacteria/release-63/variation/indexed_vep_cache/bacteria_7_collection" \
    example_bacterium_gca_000000001_vep_63_ASM3v1.tar.gz bacteria_7_collection/example_bacterium_gca_000000001/63_ASM3v1
serve "$GENOMES/fungi/release-63/variation/indexed_vep_cache/fungi_ascomycota1_collection" \
    example_fungus_gca_000000002_vep_63_ASM4v1.tar.gz fungi_ascomycota1_collection/example_fungus_gca_000000002/63_ASM4v1
serve "$GENOMES/protists/release-63/variation/indexed_vep_cache/protists_alveolata1_collection" \
    example_protist_gca_000000003_vep_63_ASM5v1.tar.gz protists_alveolata1_collection/example_protist_gca_000000003/63_ASM5v1
serve "$PANGENOME/GCA_000000004.2/ensembl/variation/2022_10/indexed_vep_cache" \
    Homo_sapiens-GCA_000000004.2-2022_10.tar.gz homo_sapiens_gca000000004v2/107_EXAMPLE.pri unlisted
serve "$PANGENOME/GCA_000000005.1/ensembl/variation/2022_10/indexed_vep_cache" \
    Homo_sapiens-GCA_000000005.1-2022-10.tar.gz homo_sapiens_gca000000005v1/107_EXAMPLE.alt unlisted

# run <arguments...>: test_genome.sh over the made-up table with the stub tools, or the
# copy of it RUN_SCRIPT names. OUT and ERR hold its two streams and CODE its exit code.
OUT="" ERR="" CODE=0
RUN_SCRIPT=$SCRIPT
run() {
    PATH="$STUB_PATH" bash "$RUN_SCRIPT" --table-dir "$TABLE" --work "$TMP/work" \
        --vep-rs "$TOOLS/vep" --generator "$TOOLS/vep-seeded-variants" \
        --comparator "$TOOLS/vep-outputs-compare" --table-tool "$TOOLS/vep-genome-table" \
        "$@" >"$TMP/stdout" 2>"$TMP/stderr"
    CODE=$?
    OUT=$(cat "$TMP/stdout")
    ERR=$(cat "$TMP/stderr")
}
fresh() {
    rm -rf "$CTL" "$TMP/work"
    mkdir -p "$CTL"
}
logged() { # logged <log of ctl>: its lines on one line, or nothing when there is none
    if [ -f "$CTL/$1" ]; then tr '\n' ' ' <"$CTL/$1"; fi
}
containers() { # how many containers the stub docker was asked to run
    local count=0
    if [ -f "$CTL/docker.log" ]; then count=$(grep -c '^run ' "$CTL/docker.log") || count=0; fi
    echo "$count"
}

# The address of each kind of set.
address() { # address <name> <genome key> <address of the archive>
    fresh
    run --key "$2"
    if [ "$OUT" = "ok $2 recheck" ] && [ "$CODE" -eq 0 ] && grep -qxF "GET https://$3" "$CTL/curl.log"; then
        ok "$1"
    else
        bad "$1: exit $CODE [$OUT] [$(logged curl.log)] [$ERR]"
    fi
}
address "a vertebrate's archive comes from the Ensembl release directory" "$VERTEBRATE" \
    "$ENSEMBL/example_vertebrate_vep_116_ASM1v1.tar.gz"
has "the job is given the row's seed" "--seed 7 " "$(logged generator.log)"
has "the job is given the row's genome and assembly" "--species example_vertebrate --assembly ASM1v1 " "$(logged vep_rs.log)"
address "a RefSeq cache's archive carries the suffix in its name" vertebrates__example_vertebrate__ASM1v1__refseq \
    "$ENSEMBL/example_vertebrate_refseq_vep_116_ASM1v1.tar.gz"
has "a RefSeq cache is annotated with --refseq" " --refseq " "$(logged vep_rs.log)"
address "a merged cache's archive carries the suffix in its name" vertebrates__example_vertebrate__ASM1v1__merged \
    "$ENSEMBL/example_vertebrate_merged_vep_116_ASM1v1.tar.gz"
address "an Ensembl Genomes division's archive comes from the division's release directory" "$PLANT" \
    "$GENOMES/plants/release-63/variation/indexed_vep_cache/example_plant_vep_63_ASM2v1.tar.gz"
address "a bacterial collection's archive comes from the collection's folder" \
    bacteria_7_collection__example_bacterium_gca_000000001__ASM3v1__ensembl \
    "$GENOMES/bacteria/release-63/variation/indexed_vep_cache/bacteria_7_collection/example_bacterium_gca_000000001_vep_63_ASM3v1.tar.gz"
address "a fungal collection's archive comes from the collection's folder" \
    fungi_ascomycota1_collection__example_fungus_gca_000000002__ASM4v1__ensembl \
    "$GENOMES/fungi/release-63/variation/indexed_vep_cache/fungi_ascomycota1_collection/example_fungus_gca_000000002_vep_63_ASM4v1.tar.gz"
address "a protist collection's archive comes from the collection's folder" \
    protists_alveolata1_collection__example_protist_gca_000000003__ASM5v1__ensembl \
    "$GENOMES/protists/release-63/variation/indexed_vep_cache/protists_alveolata1_collection/example_protist_gca_000000003_vep_63_ASM5v1.tar.gz"
address "a pangenome assembly's archive is named by the accession its genome name ends in" \
    pangenome__homo_sapiens_gca000000004v2__EXAMPLE.pri__ensembl \
    "$PANGENOME/GCA_000000004.2/ensembl/variation/2022_10/indexed_vep_cache/Homo_sapiens-GCA_000000004.2-2022_10.tar.gz"
case "$(logged curl.log)" in
*CHECKSUMS*) bad "a pangenome directory was asked for a CHECKSUMS file: $(logged curl.log)" ;;
*) ok "a pangenome directory is asked for no CHECKSUMS file" ;;
esac
address "a pangenome archive dated with a hyphen is found at the second name" \
    pangenome__homo_sapiens_gca000000005v1__EXAMPLE.alt__ensembl \
    "$PANGENOME/GCA_000000005.1/ensembl/variation/2022_10/indexed_vep_cache/Homo_sapiens-GCA_000000005.1-2022-10.tar.gz"
fresh
run --key fungi__example_unplaced__ASM14v1__ensembl
same "a row whose set sets.tsv lacks fails by the set's name" \
    "FAIL fungi__example_unplaced__ASM14v1__ensembl recheck: sets.tsv has no line for the set fungi" "$OUT"

# A download is held to its length and its checksum.
fresh
run --key vertebrates__example_badsum__ASM11v1__ensembl
has "a download whose checksum differs is refused with both sums" "where CHECKSUMS says 1 1" "$OUT"
has "the refusal is the genome's FAIL line" "FAIL vertebrates__example_badsum__ASM11v1__ensembl recheck: example_badsum_vep_116_ASM11v1.tar.gz sums to " "$OUT"
same "a refused download exits 1" 1 "$CODE"
same "nothing runs on a refused download" 0 "$(containers)$(logged generator.log)$(logged vep_rs.log)"
fresh
run --key vertebrates__example_unlisted__ASM8v1__ensembl
same "an archive CHECKSUMS does not list runs once gzip's check passes" "ok vertebrates__example_unlisted__ASM8v1__ensembl recheck" "$OUT"
has "the archive without a line is named on standard error" "no CHECKSUMS line for example_unlisted_vep_116_ASM8v1.tar.gz" "$ERR"
fresh
echo 22 >"$CTL/checksums_exit"
echo 404 >"$CTL/checksums_status"
run --key "$VERTEBRATE"
same "an archive whose directory has no CHECKSUMS file runs once gzip's check passes" "ok $VERTEBRATE recheck" "$OUT"
has "the directory without the file is named on standard error" "no CHECKSUMS file beside example_vertebrate_vep_116_ASM1v1.tar.gz" "$ERR"
fresh
echo 22 >"$CTL/checksums_exit"
echo 503 >"$CTL/checksums_status"
run --key "$VERTEBRATE"
same "a CHECKSUMS file the site answers another error status for leaves the archive unverified, a FAIL" \
    "FAIL $VERTEBRATE recheck: cannot read https://$ENSEMBL/CHECKSUMS, so example_vertebrate_vep_116_ASM1v1.tar.gz is not verified (curl: a made-up failure)" "$OUT"
same "nothing runs on an archive whose CHECKSUMS file the site fails to serve" 0 "$(containers)$(logged generator.log)"
fresh
echo 7 >"$CTL/checksums_exit"
echo 000 >"$CTL/checksums_status"
run --key "$VERTEBRATE"
same "a CHECKSUMS file that cannot be read leaves the archive unverified, a FAIL" \
    "FAIL $VERTEBRATE recheck: cannot read https://$ENSEMBL/CHECKSUMS, so example_vertebrate_vep_116_ASM1v1.tar.gz is not verified (curl: a made-up failure)" "$OUT"
same "nothing runs on an unverified archive" 0 "$(containers)$(logged generator.log)"
fresh
run --key vertebrates__example_torn__ASM12v1__ensembl
same "an archive without a line that fails gzip's check is refused" \
    "FAIL vertebrates__example_torn__ASM12v1__ensembl recheck: example_torn_vep_116_ASM12v1.tar.gz is not a whole gzip archive" "$OUT"
fresh
echo 999999 >"$CTL/stated_bytes"
run --key "$VERTEBRATE"
has "a download of another length than the site states is refused" "bytes where the site states 999999" "$OUT"
same "nothing runs on a download of another length" "" "$(logged generator.log)"
fresh
run --key vertebrates__example_absent__ASM13v1__ensembl
same "an archive the site does not serve fails by its address" \
    "FAIL vertebrates__example_absent__ASM13v1__ensembl recheck: https://$ENSEMBL does not answer for example_absent_vep_116_ASM13v1.tar.gz" "$OUT"
fresh
run --key vertebrates__example_nowhere__ASM1v1__ensembl
same "a key the table lacks fails with the table tool's line" \
    "FAIL vertebrates__example_nowhere__ASM1v1__ensembl recheck: $TABLE: no row has the key vertebrates__example_nowhere__ASM1v1__ensembl" "$OUT"
# A key names its directory under the work directory, so the key ../beside names $TMP/beside.
NO_KEY="the key is not <set>__<genome>__<assembly>__<transcripts>: it is empty, a path or an option"
fresh
mkdir -p "$TMP/beside"
echo kept >"$TMP/beside/file"
run --key "$VERTEBRATE" --key ../beside --key '' --key --field
same "a key that is a path, empty or an option fails by the form a key takes" \
    "ok $VERTEBRATE recheck|FAIL ../beside recheck: $NO_KEY|FAIL  recheck: $NO_KEY|FAIL --field recheck: $NO_KEY" \
    "$(printf '%s' "$OUT" | tr '\n' '|')"
if [ -f "$TMP/beside/file" ] && [ ! -e "$TMP/beside/download" ] && [ -f "$TMP/work/$VERTEBRATE/job.log" ]; then
    ok "such a key touches no file outside its own directory or of another genome"
else
    bad "a key that is no key left [$(find "$TMP/beside" "$TMP/work" -maxdepth 2 | tr '\n' ' ')]"
fi

# Recheck.
fresh
run --key "$VERTEBRATE" --mode recheck
same "a matching digest and record count is ok" "ok $VERTEBRATE recheck" "$OUT"
same "a run of ok lines alone exits 0" 0 "$CODE"
has "recheck hands the job the table directory" "row $TABLE $VERTEBRATE --field sha256 " "$(logged table.log)"
if [ -f "$TMP/work/$VERTEBRATE/out/row.tsv" ] && [ -f "$TMP/work/$VERTEBRATE/out/job.json" ] &&
    [ -f "$TMP/work/$VERTEBRATE/job.log" ] && [ ! -e "$TMP/work/$VERTEBRATE/work" ] &&
    [ ! -e "$TMP/work/$VERTEBRATE/download" ]; then
    ok "a genome that passes keeps its row, job record and log and drops its archive and caches"
else
    bad "a passing genome left [$(find "$TMP/work/$VERTEBRATE" -maxdepth 2 | tr '\n' ' ')]"
fi
fresh
run --key "$MISMATCH"
same "another digest is a FAIL with both digests" \
    "FAIL $MISMATCH recheck: the row holds 3 records and digest $OTHER_DIGEST; this run made 3 records and digest $DIGEST" "$OUT"
same "a FAIL exits 1" 1 "$CODE"
has "a failed run names the work directory on standard error" "each failed genome's files are under $TMP/work" "$ERR"
if [ -d "$TMP/work/$MISMATCH/work" ] && [ -d "$TMP/work/$MISMATCH/download" ]; then
    ok "a genome that fails keeps its archive and caches"
else
    bad "a failed genome left [$(find "$TMP/work/$MISMATCH" -maxdepth 1 | tr '\n' ' ')]"
fi
fresh
run --key vertebrates__example_miscount__ASM10v1__ensembl
same "another record count is a FAIL with both counts" \
    "FAIL vertebrates__example_miscount__ASM10v1__ensembl recheck: the row holds 4 records and digest $DIGEST; this run made 3 records and digest $DIGEST" "$OUT"
fresh
run --key "$VERTEBRATE" --key "$MISMATCH" --key "$EXCLUDED"
same "each genome prints its one line, in the order given" \
    "ok $VERTEBRATE recheck|FAIL $MISMATCH recheck: the row holds 3 records and digest $OTHER_DIGEST; this run made 3 records and digest $DIGEST|skip $EXCLUDED: the row is excluded:no-transcripts" \
    "$(printf '%s' "$OUT" | tr '\n' '|')"
same "one FAIL among ok and skip lines exits 1" 1 "$CODE"
has "standard error counts the three kinds of line" "test_genome: 1 ok, 1 failed, 1 skipped" "$ERR"

# Parity.
fresh
run --key "$VERTEBRATE" --mode parity
same "an outcome equal to the row's is ok" "ok $VERTEBRATE parity" "$OUT"
has "parity runs Ensembl VEP in the image" " vep --offline --cache " "$(logged docker.log)"
has "parity hands the job the comparator" "--report-dir $TMP/work/$VERTEBRATE/out/comparison" "$(logged comparator.log)"
if [ -f "$TMP/work/$VERTEBRATE/out/comparison/summary.json" ]; then
    ok "a parity genome keeps its comparison report"
else
    bad "no comparison report under $TMP/work/$VERTEBRATE/out"
fi
fresh
touch "$CTL/outputs_differ"
run --key "$VERTEBRATE" --mode parity
same "an outcome that differs from the row's is a FAIL with both" \
    "FAIL $VERTEBRATE parity: the row says exact; this run gave differs:unclassified" "$OUT"
same "a parity FAIL exits 1" 1 "$CODE"
run --key "$DIFFERS" --mode parity
same "a row that says differs: is ok while the outputs differ" "ok $DIFFERS parity" "$OUT"
fresh
run --key "$DIFFERS" --mode parity
same "a row that says differs: is a FAIL once the outputs agree" \
    "FAIL $DIFFERS parity: the row says differs:unclassified; this run gave exact" "$OUT"
fresh
touch "$CTL/ensembl_fails"
run --key "$VERTEBRATE" --mode parity
same "an exclusion the job writes for an exact row is a FAIL" \
    "FAIL $VERTEBRATE parity: the row says exact; this run gave excluded:ensembl-vep-failed" "$OUT"

# Skips and exit codes.
fresh
run --key "$EXCLUDED" --mode parity
same "an excluded row is a skip that names the exclusion" "skip $EXCLUDED: the row is excluded:no-transcripts" "$OUT"
same "an excluded row is never fetched or run" 0 "$(containers)$(logged curl.log)$(logged generator.log)"
same "a run of skips alone exits 2" 2 "$CODE"
has "a run where no genome ran says so" "ERROR: [test_genome] skip: no genome ran" "$ERR"
fresh
run --key "$DIFFERS" --mode recheck
same "recheck skips a row that is not exact" "skip $DIFFERS: the row is differs:unclassified, so there is no digest to recheck" "$OUT"
fresh
run --key "$VERTEBRATE" --key "$EXCLUDED"
same "an ok line beside a skip exits 0" 0 "$CODE"
fresh
run --key "$VERTEBRATE" --max-bytes 10
has "an archive over --max-bytes is a skip by its name" "skip $VERTEBRATE: example_vertebrate_vep_116_ASM1v1.tar.gz is " "$OUT"
has "the skip states the bound" " bytes, over --max-bytes 10" "$OUT"
case "$(logged curl.log)" in
*GET*) bad "an archive over --max-bytes was downloaded: $(logged curl.log)" ;;
*HEAD*) ok "an archive over --max-bytes is not downloaded" ;;
*) bad "the archive's length was never asked for" ;;
esac
same "a run that skipped its one genome by size exits 2" 2 "$CODE"
fresh
run --key "$VERTEBRATE" --max-bytes 100000000
same "an archive under --max-bytes runs" "ok $VERTEBRATE recheck" "$OUT"
fresh
touch "$CTL/no_image"
run --key "$VERTEBRATE" --key "$PLANT" --mode parity --vep-image example/image:1
same "a missing image skips every genome of a parity run by the image's name" \
    "skip $VERTEBRATE: no image example/image:1 in docker|skip $PLANT: no image example/image:1 in docker" \
    "$(printf '%s' "$OUT" | tr '\n' '|')"
same "a missing image exits 2" 2 "$CODE"
same "nothing is fetched without the image" "" "$(logged curl.log)"
fresh
touch "$CTL/no_image"
run --key "$VERTEBRATE" --mode recheck --vep-image example/image:1
same "recheck skips when the image it was given is missing" "skip $VERTEBRATE: no image example/image:1 in docker" "$OUT"
fresh
touch "$CTL/no_image"
run --key "$VERTEBRATE" --mode recheck
same "recheck converts on the host with the build when docker lacks the pinned image" "ok $VERTEBRATE recheck" "$OUT"
has "the conversion ran the build on the host" \
    "/ensembl_cache/example_vertebrate/116_ASM1v1 $TMP/work/$VERTEBRATE/work/json_cache" "$(logged converter.log)"
same "the job is given no image, so no container runs" 0 "$(containers)"

# The converter's build. A copy of the two scripts stands under a made-up repository
# root, which holds no build at first; the script that pins the image is one line there.
TREE="$TMP/tree"
TREE_BUILD="$TREE/tools/vep-storable-to-json/target/release/vep-storable-to-json"
mkdir -p "$TREE/scripts/genomes" "$TREE/scripts/concordance" "$TMP/no-build-target/release"
cp "$SCRIPT" "$HERE/check_genome.sh" "$TREE/scripts/genomes/"
printf 'PERL_IMAGE="made-up/image:tag"\n' >"$TREE/scripts/concordance/generate_reference_output.sh"
NO_BUILD="no build of the converter at $TREE_BUILD: build it with 'cargo build --release --locked --manifest-path $TREE/tools/vep-storable-to-json/Cargo.toml'"
RUN_SCRIPT="$TREE/scripts/genomes/test_genome.sh"
unset CARGO_TARGET_DIR
fresh
touch "$CTL/no_image"
run --key "$VERTEBRATE" --mode recheck
same "recheck skips with neither the image nor a build" \
    "skip $VERTEBRATE: no image made-up/image:tag in docker and $NO_BUILD" "$OUT"
same "a host that can convert nothing exits 2" 2 "$CODE"
same "nothing is fetched on a host that can convert nothing" "" "$(logged curl.log)"
fresh
touch "$CTL/no_image"
run --key "$VERTEBRATE" --mode parity
same "parity skips with neither the image nor a build on a line that names both" \
    "skip $VERTEBRATE: no image made-up/image:tag in docker and $NO_BUILD" "$OUT"
fresh
run --key "$VERTEBRATE" --key "$PLANT" --mode parity
same "a host with the image and no build skips every genome on the line that gives the build command" \
    "skip $VERTEBRATE: $NO_BUILD|skip $PLANT: $NO_BUILD" "$(printf '%s' "$OUT" | tr '\n' '|')"
same "nothing is fetched or run without a build" 0 "$(containers)$(logged curl.log)"
JOB_SKIP=$(PATH="$STUB_PATH" bash "$TREE/scripts/genomes/check_genome.sh" --stop-after convert \
    --cache-tarball "$SITE/$ENSEMBL/example_vertebrate_vep_116_ASM1v1.tar.gz" \
    --genome example_vertebrate --assembly ASM1v1 --transcripts ensembl --set vertebrates \
    --cache-version 116 --work "$TMP/job/work" --out "$TMP/job/out" 2>&1 >/dev/null)
same "the job's own skip names the same build and the same command" "ERROR: [check_genome] skip: $NO_BUILD" "$JOB_SKIP"
export CARGO_TARGET_DIR="$TMP/no-build-target"
fresh
run --key "$VERTEBRATE" --mode parity
same "a CARGO_TARGET_DIR that holds no build leaves the skip naming the tool's own target directory" \
    "skip $VERTEBRATE: $NO_BUILD" "$OUT"
unset CARGO_TARGET_DIR
mkdir -p "${TREE_BUILD%/*}"
printf '#!/usr/bin/env bash\nSTUB_CONVERTER_LOG=own.log exec "%s" "$@"\n' "$TMP/converter.stub" >"$TREE_BUILD"
chmod +x "$TREE_BUILD"
fresh
touch "$CTL/no_image"
run --key "$VERTEBRATE" --mode recheck
same "without CARGO_TARGET_DIR the build under the tool's own target directory is found" "ok $VERTEBRATE recheck" "$OUT"
has "the tool's own build converts on the host" "$TMP/work/$VERTEBRATE/work/json_cache" "$(logged own.log)"
RUN_SCRIPT=$SCRIPT
export CARGO_TARGET_DIR="$STUB_TARGET"

fresh
run --key "$VERTEBRATE" --mode parity
PINNED=$(sed -n 's/^PERL_IMAGE="\([^"]*\)"$/\1/p' "$IMAGE_SOURCE")
if [ -n "$PINNED" ] && grep -qxF "image inspect $PINNED" "$CTL/docker.log"; then
    ok "the image is the reference generate_reference_output.sh pins unless one is given"
else
    bad "the pinned reference [$PINNED] is not the image asked of docker: $(logged docker.log)"
fi

# The cache root.
fresh
mkdir -p "$TMP/caches"
cp "$SITE/$ENSEMBL/example_vertebrate_vep_116_ASM1v1.tar.gz" "$TMP/caches/"
run --key "$VERTEBRATE" --cache-root "$TMP/caches"
same "an archive under --cache-root runs" "ok $VERTEBRATE recheck" "$OUT"
same "nothing is requested with --cache-root" "" "$(logged curl.log)"
if [ -f "$TMP/caches/example_vertebrate_vep_116_ASM1v1.tar.gz" ]; then
    ok "an archive under --cache-root stays where it is"
else
    bad "the archive under --cache-root is gone"
fi
run --key "$PLANT" --cache-root "$TMP/caches"
same "an archive --cache-root lacks fails by its name" \
    "FAIL $PLANT recheck: no example_plant_vep_63_ASM2v1.tar.gz under $TMP/caches" "$OUT"
run --key "$VERTEBRATE" --cache-root "$TMP/caches" --max-bytes 10
has "--max-bytes holds an archive under --cache-root too" " bytes, over --max-bytes 10" "$OUT"

# The genomes named.
fresh
printf '# two genomes\n%s\n\n%s\n' "$VERTEBRATE" "$PLANT" >"$TMP/keys.txt"
run --keys-file "$TMP/keys.txt"
same "--keys-file takes one key a line and leaves out blank and # lines" \
    "ok $VERTEBRATE recheck|ok $PLANT recheck" "$(printf '%s' "$OUT" | tr '\n' '|')"
# A table of twelve genomes in one set, for the pick. Each run is a parity run without
# the image, so every picked key prints its skip line and nothing is fetched.
PICKS="$TMP/pick-table"
mkdir -p "$PICKS"
{
    echo "$HEADER"
    for n in 01 02 03 04 05 06 07 08 09 10 11 12; do
        row "example_plant_$n" ASM1v1 ensembl plants 63 3 exact "$DIGEST"
    done
} >"$PICKS/genomes.tsv"
echo "$HEADER" | gzip -c >"$PICKS/bacteria.tsv.gz"
printf 'set\trelease\tcache_version\tcaches\nplants\t63\t63\t12\n' >"$PICKS/sets.tsv"
picked() { # picked <arguments...>: the keys a run picks, on one line
    fresh
    touch "$CTL/no_image"
    run --table-dir "$PICKS" --mode parity --set plants "$@"
    printf '%s\n' "$OUT" | sed -n 's/^skip \([^:]*\): .*$/\1/p' | tr '\n' ' '
}
first=$(picked --sample 4 --seed 5)
second=$(picked --sample 4 --seed 5)
same "--set with --sample picks the same keys for the same seed" "$first" "$second"
same "the pick holds as many keys as --sample names" 4 "$(printf '%s' "$first" | wc -w | tr -d ' ')"
same "the picked keys are distinct and in key order" "$first" "$(printf '%s' "$first" | tr ' ' '\n' | sort -u | tr '\n' ' ')"
distinct=$(for seed in 1 2 3 4 5 6; do picked --sample 4 --seed "$seed"; echo; done | sort -u | wc -l | tr -d ' ')
if [ "$distinct" -gt 1 ]; then ok "other seeds pick other keys"; else bad "six seeds picked one sample"; fi
same "the seed is 1 unless one is given" "$(picked --sample 4 --seed 1)" "$(picked --sample 4)"
same "the generator's pick for seed 1 is the same on every host" \
    "plants__example_plant_01__ASM1v1__ensembl plants__example_plant_08__ASM1v1__ensembl plants__example_plant_09__ASM1v1__ensembl " \
    "$(picked --sample 3 --seed 1)"
same "a seed of 400 digits picks the keys whole-number arithmetic gives" \
    "plants__example_plant_01__ASM1v1__ensembl plants__example_plant_07__ASM1v1__ensembl plants__example_plant_10__ASM1v1__ensembl " \
    "$(picked --sample 3 --seed "$(printf '%400s' '' | tr ' ' 9)")"
same "a sample larger than the set takes every key" 12 "$(picked --sample 50 | wc -w | tr -d ' ')"

# The command line.
usage_refused() { # usage_refused <name> <text of the error> <arguments...>
    local name=$1 text=$2
    shift 2
    fresh
    run "$@"
    if [ "$CODE" -eq 2 ] && [ -z "$OUT" ]; then
        has "$name" "ERROR: [test_genome] $text" "$ERR"
    else
        bad "$name: exit $CODE [$OUT] [$ERR]"
    fi
}
usage_refused "no genome named exits 2" "no genome is named"
usage_refused "--set without --sample exits 2" "--set needs --sample" --set plants
usage_refused "--set beside --key exits 2" "--set accompanies neither --key nor --keys-file" --set plants --sample 2 --key "$PLANT"
usage_refused "--sample without --set exits 2" "--sample and --seed accompany --set" --key "$PLANT" --sample 2
usage_refused "a set without a row exits 2" "the table under $TABLE holds no row of the set fungi_none_collection" --set fungi_none_collection --sample 2
usage_refused "another mode exits 2" "--mode is recheck or parity, not stored" --key "$PLANT" --mode stored
usage_refused "an unknown argument exits 2" "unknown argument --keep" --key "$PLANT" --keep
usage_refused "a --max-bytes that is no number exits 2" "--max-bytes is a whole number of bytes" --key "$PLANT" --max-bytes 1GB
usage_refused "a binary that is no file exits 2" "--vep-rs: $TMP/none is not an executable file" --key "$PLANT" --vep-rs "$TMP/none"
usage_refused "a table directory that is none exits 2" "no directory $TMP/none" --key "$PLANT" --table-dir "$TMP/none"
fresh
run --help
if [ "$CODE" -eq 0 ] && [ -z "$ERR" ]; then
    has "--help prints the usage and exits 0" "Usage: test_genome.sh --key <genome key>..." "$OUT"
else
    bad "--help: exit $CODE [$ERR]"
fi

# The workflow, read as text.
plain() { # the workflow's lines, comment lines left out
    grep -v '^ *#' "$WORKFLOW"
}
block() { # block <top-level key>: the lines under the key, their two-space indent taken off
    plain | awk -v key="$1:" '
        $0 == key { inside = 1; next }
        inside && /^[^ ]/ { inside = 0 }
        inside && NF { sub(/^  /, ""); print }'
}
same "the workflow's only trigger is workflow_dispatch" "workflow_dispatch:" "$(block on | grep '^[^ ]' | tr '\n' ' ' | sed 's/ $//')"
found=""
for trigger in schedule push pull_request pull_request_target release workflow_run workflow_call; do
    if plain | grep -E "^  $trigger:|^on:.*[^_a-z]$trigger([^_a-z]|\$)" >/dev/null; then found="$found $trigger"; fi
done
if plain | grep 'cron:' >/dev/null; then found="$found cron"; fi
same "the workflow has no schedule, push, pull-request or release trigger" "" "$found"
same "the workflow's permissions are contents: read and nothing else" "contents: read" "$(block permissions | tr '\n' ' ' | sed 's/ $//')"
same "no job of the workflow states permissions of its own" 1 "$(plain | grep -c 'permissions:')"
uses=$(plain | grep -E '^ *(- )?uses: ')
unpinned=$(printf '%s\n' "$uses" | grep -vE 'uses: [A-Za-z0-9_./-]+@[0-9a-f]{40}( |$)')
if [ -n "$uses" ] && [ -z "$unpinned" ]; then
    ok "every action of the workflow is pinned by a 40-character commit"
else
    bad "the workflow has an action without a commit pin: [$unpinned]"
fi
strangers=""
for action in $(printf '%s\n' "$uses" | sed 's/^.*uses: \([^ ]*\).*$/\1/'); do
    grep -qF "uses: $action" "$CI" || strangers="$strangers $action"
done
same "every action is pinned at the commit ci.yml pins it at" "" "$strangers"
same "the workflow's toolchain is the one most jobs of ci.yml name" \
    "$(grep -E '^ +toolchain: ' "$CI" | awk '{ print $2 }' | sort | uniq -c | sort -rn | awk 'NR == 1 { print $2 }')" \
    "$(plain | grep -E '^ +toolchain: ' | awk '{ print $2 }' | sort -u | tr '\n' ' ' | sed 's/ $//')"
same "one step of the workflow calls test_genome.sh" 1 "$(plain | grep -c 'bash scripts/genomes/test_genome.sh ')"
same "no step runs the per-genome job, an engine or a container of its own" 0 \
    "$(plain | grep -cE 'check_genome\.sh|docker run|target/release/vep |--offline')"
has "the workflow refuses inputs that name no genome before it builds anything" "name: Refuse inputs that name no genome" \
    "$(plain | grep -nE 'name: Refuse inputs that name no genome|cargo build' | head -n 1)"
built=$(plain | grep -n -F 'run: cargo build --release --locked --manifest-path tools/vep-storable-to-json/Cargo.toml' | cut -d : -f 1)
called=$(plain | grep -n 'bash scripts/genomes/test_genome.sh ' | cut -d : -f 1)
if [ "$(printf '%s\n' "$built" | grep -c .)" -eq 1 ] && [ "$built" -lt "$called" ] 2>/dev/null; then
    ok "one step builds the cache converter before the step that calls test_genome.sh"
else
    bad "the converter is built at the workflow's lines [$built] and test_genome.sh is called at [$called]"
fi
same "the workflow has one job" 1 "$(plain | awk '/^jobs:$/ { inside = 1; next } inside && /^  [a-z_-]+:$/' | wc -l | tr -d ' ')"
same "the job runs on ubuntu-latest" 1 "$(plain | grep -c '^    runs-on: ubuntu-latest$')"
same "the job has a time limit" 1 "$(plain | grep -cE '^    timeout-minutes: [0-9]+$')"
same "the workflow types no image reference" 0 "$(grep -c 'ensemblorg/' "$WORKFLOW")"
REFERENCE_READ='sed -n '"'"'s/^PERL_IMAGE="\([^"]*\)"$/\1/p'"'"
DIGEST_READ='sed -n '"'"'s/^PERL_IMAGE_DIGEST=.*:-\(sha256:[0-9a-f]\{64\}\)}"$/\1/p'"'"
if grep -qF "$REFERENCE_READ" "$WORKFLOW" && grep -qF "$DIGEST_READ" "$WORKFLOW"; then
    ok "the workflow reads the image reference and its digest from the tree"
else
    bad "the workflow does not read the reference and the digest with the expressions this test runs"
fi
reference=$(sed -n 's/^PERL_IMAGE="\([^"]*\)"$/\1/p' "$IMAGE_SOURCE")
digest=$(sed -n 's/^PERL_IMAGE_DIGEST=.*:-\(sha256:[0-9a-f]\{64\}\)}"$/\1/p' "$IMAGE_SOURCE")
if [ "$(printf '%s\n' "$reference" | grep -c .)" -eq 1 ] && [ "$(printf '%s\n' "$digest" | grep -c .)" -eq 1 ]; then
    ok "those expressions find one reference and one digest in the tree"
else
    bad "the expressions found the reference [$reference] and the digest [$digest]"
fi
inputs=$(plain | awk '/^    inputs:$/ { inside = 1; next } inside && /^    [^ ]|^[^ ]/ { inside = 0 } inside && /^      [a-z_]+:$/ { sub(/^ +/, ""); sub(/:$/, ""); print }' | tr '\n' ' ')
same "the workflow takes keys, set, sample, mode and max_bytes" "keys set sample mode max_bytes " "$inputs"
typed=$(plain | grep -F '${{ inputs.' | grep -cvE '^ +[A-Z_]+: \$\{\{ inputs\.[a-z_]+ \}\}$') || typed=0
same "an input reaches a step through its environment alone" 0 "$typed"
same "the artifact holds each genome's row, job record and comparison report" 3 \
    "$(plain | grep -cE '/genomes/\*/out/(row\.tsv|job\.json|comparison/)$')"

echo "passed: $PASS   failed: $FAIL"
[ "$FAIL" -eq 0 ]
