#!/usr/bin/env bash
# Tests that capture_baseline.sh computes its metrics and its cache checksums with the file
# census of its own tree, tools/vep-file-census, and with nothing else. The script is run
# whole in trees of its own with a stub vep-rs: where the tool is not built (nothing at the
# build's path, a file that is not executable, a directory), with a stub tool that records
# the command lines and the standard input it was started with, and once on this tree's own
# release build, a check that is skipped where the tool is not built.
#
# Run: bash scripts/validation/test_capture_baseline.sh

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

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAPTURE="$HERE/capture_baseline.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

CENSUS_REL="tools/vep-file-census/target/release/vep-file-census"

bash -n "$CAPTURE" && ok "capture_baseline.sh parses" || bad "capture_baseline.sh has a syntax error"

echo "FILE-CENSUS: capture_baseline.sh finds the tree's release build before vep-rs runs and reads the output and the cache through it"

grep -q "^census_bin=\"\\\$REPO_DIR/$CENSUS_REL\"\$" "$CAPTURE" &&
    ok "the tree's release build of tools/vep-file-census is the script's one candidate" ||
    bad "the script does not name the tree's release build of tools/vep-file-census"
_block_line=$(grep -n -e '^census_bin=' "$CAPTURE" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^# Step 1:' "$CAPTURE" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "the file census is found before vep-rs runs (line $_block_line before $_first_run)" ||
    bad "the file census block must stand before Step 1 (block at ${_block_line:-?}, step at ${_first_run:-?})"
_named=$(grep -vE '^[[:space:]]*#' "$CAPTURE" | grep -c -E -e 'python3|md5sum|md5 -r') || _named=0
[[ "$_named" -eq 0 ]] &&
    ok "no live line of the script names an interpreter or a digest program" ||
    bad "$_named live line(s) of the script name python3, md5sum or md5"

# A cache of four JSON files, one of them under a name with a space, and one file that is
# no JSON file; an input VCF; a stub vep-rs that marks its start and writes three rows.
CACHE="$TMP/cache"
mkdir -p "$CACHE/transcripts/21" "$CACHE/with space" "$TMP/bin"
printf '{"a": 1}\n' >"$CACHE/info.json"
printf '{"b": 2}\n' >"$CACHE/transcripts/21/1-1000000.json"
: >"$CACHE/transcripts/21/empty.json"
printf '{"c": 3}' >"$CACHE/with space/c d.json"
printf 'not a cache file\n' >"$CACHE/chr_synonyms.txt"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1000\tv1\tA\tG\t.\t.\t.\n' >"$TMP/input.vcf"
HEADER='#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n'
ROWS='v1\t21:1000\tG\tENSG01\tENST01\tTranscript\tmissense_variant,splice_region_variant\t10\t10\t4\tK/R\taAg/aGg\trs1\tIMPACT=MODERATE;CLIN_SIG=benign;AF=0.5;HGVSc=ENST01:c.10A>G\nv1\t21:1000\tG\tENSG02\tENST02\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs1\tIMPACT=MODIFIER;AF=0.5;gnomADe_AF=0.25\nv2\t21:2000\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n'
# The row formats are the output's bytes.
# shellcheck disable=SC2059
{
    printf '## ENSEMBL VARIANT EFFECT PREDICTOR\n' && printf "$HEADER" && printf "$ROWS"
} >"$TMP/engine_output.txt"
cat >"$TMP/bin/vep" <<STUB
#!/usr/bin/env bash
out=""; prev=""
for a in "\$@"; do [[ "\$prev" == -o ]] && out="\$a"; prev="\$a"; done
: >"$TMP/engine.ran"
cp "$TMP/engine_output.txt" "\$out"
STUB
chmod +x "$TMP/bin/vep"

# make_tree <name>: a tree that holds a copy of the script and nothing at the build's path.
make_tree() {
    mkdir -p "$TMP/$1/scripts/validation" "$TMP/$1/$(dirname "$CENSUS_REL")"
    cp "$CAPTURE" "$TMP/$1/scripts/validation/"
}
# capture_in <tree> <output dir>: the copy of the script in the tree, whole.
capture_in() {
    rm -rf "$2" "$TMP/engine.ran"
    bash "$TMP/$1/scripts/validation/capture_baseline.sh" --cache-dir "$CACHE" --vcf "$TMP/input.vcf" --vep-binary "$TMP/bin/vep" --output-dir "$2"
}

for _shape in unbuilt inert directory; do
    make_tree "tree-$_shape"
    case "$_shape" in
    inert) : >"$TMP/tree-$_shape/$CENSUS_REL" ;;
    directory) mkdir -p "$TMP/tree-$_shape/$CENSUS_REL" ;;
    esac
    _out=$(capture_in "tree-$_shape" "$TMP/out-$_shape" 2>&1)
    rc=$?
    [[ "$rc" -eq 1 && "$_out" == "ERROR: [capture_baseline] the file census $TMP/tree-$_shape/$CENSUS_REL is not an executable: build the tool with 'cargo build --release --locked --manifest-path $TMP/tree-$_shape/tools/vep-file-census/Cargo.toml'" && ! -e "$TMP/engine.ran" && ! -e "$TMP/out-$_shape" ]] &&
        ok "a run in a tree whose build is $_shape: exit 1 on one line naming the build command, vep-rs never started, no output directory made" ||
        bad "a run in a tree whose build is $_shape was not refused before vep-rs (rc=$rc): $_out"
done

# A tool that records each start, its arguments on one line and its standard input after.
make_tree tree-recorded
cat >"$TMP/tree-recorded/$CENSUS_REL" <<STUB
#!/usr/bin/env bash
printf '%s\\n' "\$*" >>"$TMP/started"
[[ "\$1" == digest ]] && cat >"$TMP/names"
exit 0
STUB
chmod +x "$TMP/tree-recorded/$CENSUS_REL"
: >"$TMP/started"
_out=$(capture_in tree-recorded "$TMP/out-recorded" 2>&1)
rc=$?
_metrics_start=$(sed -n '1p' "$TMP/started")
[[ "$rc" -eq 0 && -e "$TMP/engine.ran" && "$(wc -l <"$TMP/started" | tr -d ' ')" -eq 2 && "$_metrics_start" =~ ^metrics\ $TMP/out-recorded/baseline_output\.txt\ $TMP/out-recorded/baseline_metrics\.json\ [0-9]+$ ]] &&
    ok "after vep-rs, Step 2 starts the tool once on the output, the report's path and the elapsed seconds" ||
    bad "Step 2 did not start the tool as expected (rc=$rc): $(tr '\n' '|' <"$TMP/started")"
[[ "$(sed -n '2p' "$TMP/started")" == "digest --algorithm md5 --names read" && "$(cat "$TMP/names")" == "$(find "$CACHE" -type f -name '*.json' | sort)" ]] &&
    ok "Step 3 starts the tool once, on the sorted names of the cache's JSON files, one a line on its standard input" ||
    bad "Step 3 did not start the tool on the cache's JSON files: $(sed -n '2p' "$TMP/started") with $(tr '\n' '|' <"$TMP/names")"
grep -q "^  Wrote 0 checksums to $TMP/out-recorded/cache_checksums.txt\$" <<<"$_out" &&
    ok "the checksum count is the line count of what the tool wrote" ||
    bad "the checksum count line is missing: $(grep 'checksums' <<<"$_out")"

# The whole script on this tree's own release build.
REAL_CENSUS="$(cd "$HERE/../.." && pwd)/$CENSUS_REL"
if [[ -x "$REAL_CENSUS" ]]; then
    make_tree tree-real
    ln -s "$REAL_CENSUS" "$TMP/tree-real/$CENSUS_REL"
    _out=$(capture_in tree-real "$TMP/out-real" 2>&1)
    rc=$?
    [[ "$rc" -eq 0 ]] && ok "a run on the tree's own build exits 0" || bad "a run on the tree's own build exited $rc: $(tail -3 <<<"$_out" | tr '\n' ' ')"
    _summary=$(sed -n '/^Step 2:/,/^Step 3:/p' <<<"$_out" | sed '1d;$d')
    [[ "$_summary" == "  Total annotation lines: 3
  Unique consequence types: 4
  Existing_variation non-empty: 2
  CLIN_SIG non-empty: 1
  HGVSc non-empty: 1
  HGVSp non-empty: 0
  Max frequency field count: 2" ]] &&
        ok "Step 2 prints the seven summary lines of the output's three rows" ||
        bad "Step 2 printed: $(tr '\n' '|' <<<"$_summary")"
    _elapsed=$(sed -n 's/^  "elapsed_seconds": \([0-9]*\)$/\1/p' "$TMP/out-real/baseline_metrics.json")
    _want='{
  "total_annotation_lines": 3,
  "consequence_counts": {
    "missense_variant": 1,
    "splice_region_variant": 1,
    "intron_variant": 1,
    "intergenic_variant": 1
  },
  "existing_variation_count": 2,
  "clin_sig_count": 1,
  "hgvsc_count": 1,
  "hgvsp_count": 0,
  "frequency_field_counts": {
    "AF": 2,
    "gnomADe_AF": 1,
    "gnomADg_AF": 0,
    "MAX_AF": 0,
    "EUR_AF": 0,
    "AFR_AF": 0,
    "AMR_AF": 0,
    "EAS_AF": 0,
    "SAS_AF": 0
  },
  "elapsed_seconds": '"$_elapsed"'
}'
    [[ -n "$_elapsed" && "$(cat "$TMP/out-real/baseline_metrics.json")" == "$_want" ]] &&
        ok "baseline_metrics.json holds the row count, the terms in the order first met, the filled keys and the elapsed seconds" ||
        bad "baseline_metrics.json is: $(tr '\n' ' ' <"$TMP/out-real/baseline_metrics.json")"
    if command -v md5sum >/dev/null 2>&1; then
        _want_sums=$(find "$CACHE" -type f -name '*.json' | sort | while read -r f; do md5sum "$f"; done)
        [[ "$(cat "$TMP/out-real/cache_checksums.txt")" == "$_want_sums" ]] && grep -q "^  Wrote 4 checksums to $TMP/out-real/cache_checksums.txt\$" <<<"$_out" &&
            ok "cache_checksums.txt is one md5sum line a JSON file of the cache, in sorted order, the file under a spaced name among them" ||
            bad "cache_checksums.txt is not the md5sum lines of the cache's JSON files: $(tr '\n' '|' <"$TMP/out-real/cache_checksums.txt")"
    else
        echo "  skip the checksum lines against md5sum: no md5sum on this host"
    fi
else
    echo "  skip the run on the tree's own build: no release build of tools/vep-file-census"
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
