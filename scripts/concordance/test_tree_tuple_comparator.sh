#!/usr/bin/env bash
# Tests that run_concordance.sh and scripts/validation/compare_cache_outputs.sh score tuples
# with the tuple comparator of their own tree, tools/vep-outputs-compare, and with nothing
# else. Each script's comparator block and its comparator call are extracted and run against
# stub outputs, with a stub cargo, a stub python3 and stub comparators that record the
# command line they were started with; each script is also run whole where its comparator
# cannot run, and the comparison step of compare_cache_outputs.sh runs once on this tree's own
# release build, a check that is skipped where the tool is not built. The input preparation
# of run_concordance.sh is held the same way: its block and its step are extracted and run
# with the stub cargo, which shows the tree's tools/vep-prepare-benchmark-vcfs built and then
# run, and a failing build refused.
# So is its file census: its block and the two functions that digest a directory and a file
# are extracted and run with the stub cargo, which shows the tree's tools/vep-file-census
# built and then run, and a failing build refused; both functions run once on this tree's own
# release build, a check that is skipped where the tool is not built.
#
# The cases, for each script:
#   TREE      the Rust tool of the tree runs: built and run through cargo by
#             run_concordance.sh, the tree's release build by compare_cache_outputs.sh.
#   EXPORTED  an exported VEP_OUTPUTS_COMPARATOR or VEP_FIELDS_COMPARATOR changes nothing,
#             whatever it names (the word `python`, an executable, a file that is not one):
#             the same commands start and nothing else does.
#   NO-PYTHON no line of the script names compare_vep_outputs.py, compare_vep_fields.py or
#             either variable.
#   REFUSAL   a comparator that cannot run (a build that fails, a tool that is not built or
#             is not an executable file) stops the script with one line naming the build
#             command, before any engine run, and nothing runs in its place.
#
# And for run_concordance.sh alone, the release it hands the comparator:
#   RELEASE   a --reference-release outside the accepted releases stops the script on one
#             line naming them; an accepted one goes on to the next check.
#
# Run: bash scripts/concordance/test_tree_tuple_comparator.sh

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
CONCORDANCE="$HERE/run_concordance.sh"
CACHE_OUTPUTS="$HERE/../validation/compare_cache_outputs.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

TOOL_REL="tools/vep-outputs-compare/target/release/vep-outputs-compare"
LOG="$TMP/started"
mkdir -p "$TMP/bin" "$TMP/built/$(dirname "$TOOL_REL")" "$TMP/unbuilt" "$TMP/inert_build/$(dirname "$TOOL_REL")" "$TMP/dir_build/$TOOL_REL"
# Each stub appends its label and arguments, one line a start. The cargo stub fails a
# `build` while CARGO_BUILD_FAILS is set.
stub() {
    printf '#!/bin/bash\nprintf "%%s\\n" "%s $*" >>"%s"\n%s\n' "$2" "$LOG" "${3:-}" >"$1"
    chmod +x "$1"
}
stub "$TMP/bin/cargo" cargo '[[ "$1" == build && -n "${CARGO_BUILD_FAILS:-}" ]] && exit 101; exit 0'
stub "$TMP/bin/python3" python3
stub "$TMP/bin/comparator" comparator
stub "$TMP/bin/vep" vep
stub "$TMP/built/$TOOL_REL" tool
: >"$TMP/not_executable"
: >"$TMP/bin/inert"
: >"$TMP/inert_build/$TOOL_REL"

# What a caller may have exported, and what each form is. "" leaves both variables unset,
# `empty` exports them empty, any other word is exported as the value of both.
EXPORTED=(empty python "$TMP/bin/comparator" comparator "$TMP/not_executable" "$TMP" inert)
EXPORTED_AS=("the empty string" "the word python" "an executable file" "a bare name on PATH" "a file that is not executable" "a directory" "a bare name whose PATH match is not executable")
export_comparators() {
    case "$1" in
    "") unset VEP_OUTPUTS_COMPARATOR VEP_FIELDS_COMPARATOR ;;
    empty) export VEP_OUTPUTS_COMPARATOR="" VEP_FIELDS_COMPARATOR="" ;;
    *) export VEP_OUTPUTS_COMPARATOR="$1" VEP_FIELDS_COMPARATOR="$1" ;;
    esac
}
# $1 is a script; prints how many of its lines name a Python tuple or field comparator or
# one of the two variables.
names_python_or_variable() {
    local n
    n=$(grep -c -E -e 'compare_vep_(outputs|fields)\.py|VEP_(OUTPUTS|FIELDS)_COMPARATOR' "$1") || n=0
    echo "$n"
}

echo "RUN-CONCORDANCE: run_concordance.sh builds the tree's tuple comparator before any engine runs and compares through it"

_block=$(awk '
    !p && /^OUTPUTS_COMPARE_MANIFEST=/ { p = 1 }
    p { print }
    p && /^OUTPUTS_COMPARE=\(/ { exit }
' "$CONCORDANCE")
[[ "$_block" == *'if ! cargo build -q --release --locked --manifest-path "${OUTPUTS_COMPARE_MANIFEST}"; then'* && "$_block" == *$'\nOUTPUTS_COMPARE=(cargo run -q --release --locked --manifest-path "${OUTPUTS_COMPARE_MANIFEST}" --)' ]] &&
    ok "run_concordance.sh: the comparator block is extracted whole, a cargo build then the cargo run command" ||
    bad "run_concordance.sh: could not extract the comparator block"
_call=$(grep -e '^run_timed "compare_outputs" ' "$CONCORDANCE")
[[ "$_call" == 'run_timed "compare_outputs" "${OUTPUTS_COMPARE[@]}" "${COMPARE_ARGS[@]}"' ]] &&
    ok "run_concordance.sh: the one compare step runs the tree's command under its timer" ||
    bad "run_concordance.sh: expected one compare step through OUTPUTS_COMPARE, got '$_call'"
_named=$(names_python_or_variable "$CONCORDANCE")
[[ "$_named" -eq 0 ]] &&
    ok "run_concordance.sh: no line names a Python tuple or field comparator, VEP_OUTPUTS_COMPARATOR or VEP_FIELDS_COMPARATOR" ||
    bad "run_concordance.sh: $_named line(s) name a Python comparator or a comparator variable"
_block_line=$(grep -n -e '^OUTPUTS_COMPARE_MANIFEST=' "$CONCORDANCE" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^[[:space:]]*run_timed "' "$CONCORDANCE" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "run_concordance.sh: the comparator is built before the first timed step (line $_block_line before $_first_run)" ||
    bad "run_concordance.sh: the comparator block must stand before the first timed step (block at ${_block_line:-?}, step at ${_first_run:-?})"

_args="--perl-dir $TMP/work/perl --rust-dir $TMP/work/rust --report-dir $TMP/work/reports"
_manifest="$TMP/built/tools/vep-outputs-compare/Cargo.toml"
# $1 is what the caller exported; runs the block and the compare step; prints stderr, then
# `rc=`.
run_concordance_compare() {
    : >"$LOG"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034,SC2329
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$TMP/built"
        export_comparators "$1"
        run_timed() {
            shift
            "$@"
        }
        # The compare step's arguments, split on spaces as the stub lines are.
        # shellcheck disable=SC2206
        COMPARE_ARGS=($_args)
        eval "$_block"
        eval "$_call"
    ) 2>&1
    echo "rc=$?"
}
_tree_build="cargo build -q --release --locked --manifest-path $_manifest"
_tree_run="cargo run -q --release --locked --manifest-path $_manifest -- $_args"
_build_refusal="ERROR: [run_concordance] the tuple comparator tools/vep-outputs-compare did not build: repair 'cargo build --release --locked --manifest-path $_manifest'"

_out=$(run_concordance_compare "")
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "$_tree_build"$'\n'"$_tree_run" ]] &&
    ok "run_concordance.sh: the Rust tool of the tree is built, then run on the outputs, through cargo" ||
    bad "run_concordance.sh: expected a cargo build then a cargo run of the tool, got '$_out' and '$(cat "$LOG")'"
for _i in "${!EXPORTED[@]}"; do
    _out=$(run_concordance_compare "${EXPORTED[$_i]}")
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "$_tree_build"$'\n'"$_tree_run" ]] &&
        ok "run_concordance.sh, both comparator variables exported as ${EXPORTED_AS[$_i]}: the same cargo build and run of the tree's tool, and nothing else starts" ||
        bad "run_concordance.sh with ${EXPORTED[$_i]} exported: expected the cargo build and run alone, got '$_out' and '$(cat "$LOG")'"
done
_out=$(CARGO_BUILD_FAILS=1 run_concordance_compare "")
[[ "$_out" == "$_build_refusal"$'\n'"rc=1" && "$(cat "$LOG")" == "$_tree_build" ]] &&
    ok "run_concordance.sh, the tool's build fails: one line names the build command, exit 1, nothing is compared and nothing runs in the tool's place" ||
    bad "run_concordance.sh with a failing build: expected the refusal, rc=1 and the build alone, got '$_out' and '$(cat "$LOG")'"
_out=$(CARGO_BUILD_FAILS=1 run_concordance_compare python)
[[ "$_out" == "$_build_refusal"$'\n'"rc=1" && "$(cat "$LOG")" == "$_tree_build" ]] &&
    ok "run_concordance.sh, the tool's build fails with the word python exported: the same refusal, and python3 does not start" ||
    bad "run_concordance.sh with a failing build and python exported: expected the refusal, rc=1 and the build alone, got '$_out' and '$(cat "$LOG")'"
# The whole script, not the block: a run whose comparator does not build stops on the
# refusal with the build as the one command started, so before any engine runs.
_real_manifest="$(cd "$HERE/../.." && pwd)/tools/vep-outputs-compare/Cargo.toml"
mkdir -p "$TMP/perl_cache"
: >"$LOG"
_out=$(CARGO_BUILD_FAILS=1 PATH="$TMP/bin:$PATH" bash "$CONCORDANCE" --perl-cache-dir "$TMP/perl_cache" --no-fasta 2>&1)
rc=$?
[[ "$rc" -eq 1 && "$(printf '%s\n' "$_out" | grep -c '^ERROR')" -eq 1 && "$(printf '%s\n' "$_out" | tail -1)" == "ERROR: [run_concordance] the tuple comparator tools/vep-outputs-compare did not build: repair 'cargo build --release --locked --manifest-path $_real_manifest'" && "$(cat "$LOG")" == "cargo build -q --release --locked --manifest-path $_real_manifest" ]] &&
    ok "a run_concordance.sh run whose comparator does not build: exit 1 on the one refusal line, the build the only command started" ||
    bad "run_concordance.sh whole with a failing build: expected rc=1, the refusal as the one ERROR and last line and the build alone, got rc=$rc '$(printf '%s\n' "$_out" | tail -2)' and '$(cat "$LOG")'"

echo "RUN-CONCORDANCE-PREPARATION: run_concordance.sh builds the tree's input preparation tool before any engine runs and prepares through it"

_prep_block=$(awk '
    !p && /^PREPARE_VCFS_MANIFEST=/ { p = 1 }
    p { print }
    p && /^PREPARE_VCFS=\(/ { exit }
' "$CONCORDANCE")
[[ "$_prep_block" == *'if ! cargo build -q --release --locked --manifest-path "${PREPARE_VCFS_MANIFEST}"; then'* && "$_prep_block" == *$'\nPREPARE_VCFS=(cargo run -q --release --locked --manifest-path "${PREPARE_VCFS_MANIFEST}" --)' ]] &&
    ok "run_concordance.sh: the preparation block is extracted whole, a cargo build then the cargo run command" ||
    bad "run_concordance.sh: could not extract the preparation block"
_prep_call=$(grep -e '^run_timed "prepare_benchmark_vcfs" ' "$CONCORDANCE")
[[ "$_prep_call" == 'run_timed "prepare_benchmark_vcfs" "${PREPARE_VCFS[@]}" "${PREP_ARGS[@]}"' ]] &&
    ok "run_concordance.sh: the one preparation step runs the tree's command under its timer" ||
    bad "run_concordance.sh: expected one preparation step through PREPARE_VCFS, got '$_prep_call'"
_prep_line=$(grep -n -e '^PREPARE_VCFS_MANIFEST=' "$CONCORDANCE" | head -1 | cut -d: -f1)
_first_timed=$(grep -n -e '^[[:space:]]*run_timed "' "$CONCORDANCE" | head -1 | cut -d: -f1)
[[ -n "$_prep_line" && -n "$_first_timed" && "$_prep_line" -lt "$_first_timed" ]] &&
    ok "run_concordance.sh: the preparation tool is built before the first timed step (line $_prep_line before $_first_timed)" ||
    bad "run_concordance.sh: the preparation block must stand before the first timed step (block at ${_prep_line:-?}, step at ${_first_timed:-?})"

_prep_args="--input-dir $TMP/work/benchmark --output-dir $TMP/work/inputs --chrom-mode auto --canonical-contigs"
_prep_manifest="$TMP/built/tools/vep-prepare-benchmark-vcfs/Cargo.toml"
# Runs the block and the preparation step; prints stderr, then `rc=`.
run_concordance_prepare() {
    : >"$LOG"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034,SC2329
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$TMP/built"
        run_timed() {
            shift
            "$@"
        }
        # The preparation step's arguments, split on spaces as the stub lines are.
        # shellcheck disable=SC2206
        PREP_ARGS=($_prep_args)
        eval "$_prep_block"
        eval "$_prep_call"
    ) 2>&1
    echo "rc=$?"
}
_prep_build="cargo build -q --release --locked --manifest-path $_prep_manifest"
_out=$(run_concordance_prepare)
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "$_prep_build"$'\n'"cargo run -q --release --locked --manifest-path $_prep_manifest -- $_prep_args" ]] &&
    ok "run_concordance.sh: the preparation tool of the tree is built, then run on the benchmark directory, through cargo, and nothing else starts" ||
    bad "run_concordance.sh: expected a cargo build then a cargo run of the preparation tool, got '$_out' and '$(cat "$LOG")'"
_out=$(CARGO_BUILD_FAILS=1 run_concordance_prepare)
[[ "$_out" == "ERROR: [run_concordance] the input preparation tool tools/vep-prepare-benchmark-vcfs did not build: repair 'cargo build --release --locked --manifest-path $_prep_manifest'"$'\n'"rc=1" && "$(cat "$LOG")" == "$_prep_build" ]] &&
    ok "run_concordance.sh, the preparation tool's build fails: one line names the build command, exit 1, nothing is prepared and nothing runs in the tool's place" ||
    bad "run_concordance.sh with a failing preparation build: expected the refusal, rc=1 and the build alone, got '$_out' and '$(cat "$LOG")'"

echo "RUN-CONCORDANCE-CENSUS: run_concordance.sh builds the tree's file census before any engine runs and digests the caches, the FASTA and the inputs through it"

_census_block=$(awk '
    !p && /^FILE_CENSUS_MANIFEST=/ { p = 1 }
    p { print }
    p && /^FILE_CENSUS=\(/ { exit }
' "$CONCORDANCE")
[[ "$_census_block" == *'if ! cargo build -q --release --locked --manifest-path "${FILE_CENSUS_MANIFEST}"; then'* && "$_census_block" == *$'\nFILE_CENSUS=(cargo run -q --release --locked --manifest-path "${FILE_CENSUS_MANIFEST}" --)' ]] &&
    ok "run_concordance.sh: the file census block is extracted whole, a cargo build then the cargo run command" ||
    bad "run_concordance.sh: could not extract the file census block"
_census_line=$(grep -n -e '^FILE_CENSUS_MANIFEST=' "$CONCORDANCE" | head -1 | cut -d: -f1)
[[ -n "$_census_line" && -n "$_first_timed" && "$_census_line" -lt "$_first_timed" && "$_census_line" -gt "$_prep_line" ]] &&
    ok "run_concordance.sh: the file census is built after the comparator and the preparation tool and before the first timed step (line $_census_line before $_first_timed)" ||
    bad "run_concordance.sh: the file census block must stand after the preparation block and before the first timed step (block at ${_census_line:-?}, step at ${_first_timed:-?})"
_digests=$(grep -c -F -e '"${FILE_CENSUS[@]}" digest --algorithm sha256 ' "$CONCORDANCE") || _digests=0
_whole_file=$(grep -vE '^[[:space:]]*#' "$CONCORDANCE" | grep -c -E -e 'xargs -0 sha256sum|sha256sum "\$|hashlib\.sha256\(\)') || _whole_file=0
[[ "$_digests" -eq 4 && "$_whole_file" -eq 0 ]] &&
    ok "run_concordance.sh: the four digests of whole files (a cache tree, the FASTA, a prepared input twice) run the tree's command, and no line digests a whole file otherwise" ||
    bad "run_concordance.sh: expected 4 digests through FILE_CENSUS and none through sha256sum or hashlib, found $_digests and $_whole_file"

_dir_sha=$(sed -n '/^dir_sha256() {/,/^}/p' "$CONCORDANCE")
_file_sha=$(sed -n '/^sha256_file() {/,/^}/p' "$CONCORDANCE")
_census_manifest="$TMP/built/tools/vep-file-census/Cargo.toml"
mkdir -p "$TMP/census_dir/sub"
printf 'one\n' >"$TMP/census_dir/a file"
printf 'two\n' >"$TMP/census_dir/sub/b"
: >"$TMP/census_dir/sub/empty"
# Runs the block, then the digest of a directory and of one file; prints stderr, then `rc=`.
run_concordance_census() {
    : >"$LOG"
    # The extracted block and functions read the names this subshell sets.
    # shellcheck disable=SC2034
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$TMP/built"
        eval "$_census_block"
        eval "$_dir_sha"
        eval "$_file_sha"
        dir_sha256 "$TMP/census_dir" >/dev/null
        sha256_file "$TMP/census_dir/a file" >/dev/null
    ) 2>&1
    echo "rc=$?"
}
_census_build="cargo build -q --release --locked --manifest-path $_census_manifest"
_census_run="cargo run -q --release --locked --manifest-path $_census_manifest --"
_out=$(run_concordance_census)
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "$_census_build"$'\n'"$_census_run digest --algorithm sha256 --names nul --directory $TMP/census_dir"$'\n'"$_census_run digest --algorithm sha256 -" ]] &&
    ok "run_concordance.sh: the file census of the tree is built, then run through cargo on a directory's sorted names and on one file's bytes, and nothing else starts" ||
    bad "run_concordance.sh: expected a cargo build then two cargo runs of the file census, got '$_out' and '$(cat "$LOG")'"
_out=$(CARGO_BUILD_FAILS=1 run_concordance_census)
[[ "$_out" == "ERROR: [run_concordance] the file census tools/vep-file-census did not build: repair 'cargo build --release --locked --manifest-path $_census_manifest'"$'\n'"rc=1" && "$(cat "$LOG")" == "$_census_build" ]] &&
    ok "run_concordance.sh, the file census's build fails: one line names the build command, exit 1, nothing is digested and nothing runs in the tool's place" ||
    bad "run_concordance.sh with a failing file census build: expected the refusal, rc=1 and the build alone, got '$_out' and '$(cat "$LOG")'"
# The two functions on this tree's own release build: the values are those of the
# pipelines the tool stands for.
_real_census="$(cd "$HERE/../.." && pwd)/tools/vep-file-census/target/release/vep-file-census"
if [[ -x "$_real_census" ]] && command -v sha256sum >/dev/null 2>&1; then
    # The extracted functions read the name this subshell sets.
    # shellcheck disable=SC2034
    _got=$(
        FILE_CENSUS=("$_real_census")
        eval "$_dir_sha"
        eval "$_file_sha"
        dir_sha256 "$TMP/census_dir"
        sha256_file "$TMP/census_dir/a file"
    )
    _want_dir=$(cd "$TMP/census_dir" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}')
    _want_file=$(sha256sum <"$TMP/census_dir/a file" | awk '{print $1}')
    [[ "$_got" == "$_want_dir"$'\n'"$_want_file" && "${#_want_dir}" -eq 64 ]] &&
        ok "run_concordance.sh on the tree's own release build: a directory's digest is that of find, sort, xargs sha256sum and sha256sum, a file's that of sha256sum" ||
        bad "run_concordance.sh on the tree's own build: expected '$_want_dir' and '$_want_file', got '$_got'"
else
    echo "  skip the digests on the tree's own build: no release build of tools/vep-file-census, or no sha256sum"
fi

echo "REFERENCE-RELEASE: run_concordance.sh holds the comparator's --reference-release to the accepted releases"

# The release is checked before the FASTA: with a cache directory and no FASTA, an accepted
# release stops the script on the FASTA's line and any other on the release's own.
for _release in 115.2 116.2 116.3; do
    bash "$CONCORDANCE" --perl-cache-dir "$TMP" --reference-release "$_release" >"$TMP/release.log" 2>&1
    [[ "$(head -1 "$TMP/release.log")" == "ERROR: --fasta is required (or pass --no-fasta to opt out explicitly)" ]] &&
        ok "run_concordance.sh, --reference-release $_release: accepted, the script goes on to its FASTA check" ||
        bad "run_concordance.sh --reference-release $_release: stopped on '$(head -1 "$TMP/release.log")'"
done
bash "$CONCORDANCE" --perl-cache-dir "$TMP" --reference-release 116.1 >"$TMP/release.log" 2>&1
_rc=$?
[[ "$_rc" -eq 1 && "$(cat "$TMP/release.log")" == "ERROR: --reference-release must be 115.2, 116.2 or 116.3 (got 116.1)" ]] &&
    ok "run_concordance.sh, --reference-release 116.1: refused on one line naming the three accepted releases, exit 1" ||
    bad "run_concordance.sh --reference-release 116.1: expected the refusal and rc=1, got rc=$_rc and '$(cat "$TMP/release.log")'"

echo "CACHE-OUTPUTS: compare_cache_outputs.sh finds the tree's release build before vep-rs runs and compares through it"

_cache_block=$(awk '
    !p && /^OUTPUTS_COMPARE=\(\)$/ { p = 1 }
    p { print }
    p && /^fi$/ { exit }
' "$CACHE_OUTPUTS")
[[ "$_cache_block" == *"outputs_bin=\"\${VEP_RS_DIR}/$TOOL_REL\""* && "$_cache_block" == *$'\nfi' ]] &&
    ok "compare_cache_outputs.sh: the comparator block is extracted whole, the tree's release build its one candidate" ||
    bad "compare_cache_outputs.sh: could not extract the comparator block"
_step=$(awk '
    /^# Step 5:/ { exit }
    /^# Step 4:/ { p = 1 }
    p { print }
' "$CACHE_OUTPUTS")
[[ "$_step" == *'"${OUTPUTS_COMPARE[@]}" \'* && "$_step" != *python3* ]] &&
    ok "compare_cache_outputs.sh: the comparison step runs the tree's command and names no interpreter" ||
    bad "compare_cache_outputs.sh: the comparison step must run OUTPUTS_COMPARE alone"
_named=$(names_python_or_variable "$CACHE_OUTPUTS")
[[ "$_named" -eq 0 ]] &&
    ok "compare_cache_outputs.sh: no line names a Python tuple or field comparator, VEP_OUTPUTS_COMPARATOR or VEP_FIELDS_COMPARATOR" ||
    bad "compare_cache_outputs.sh: $_named line(s) name a Python comparator or a comparator variable"
_block_line=$(grep -n -e '^OUTPUTS_COMPARE=()$' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^# Step 1:' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "compare_cache_outputs.sh: the comparator is found before the first vep-rs run (line $_block_line before $_first_run)" ||
    bad "compare_cache_outputs.sh: the comparator block must stand before Step 1 (block at ${_block_line:-?}, step at ${_first_run:-?})"

_work="$TMP/cache_work"
mkdir -p "$_work/baseline" "$_work/builder"
printf '#Uploaded_variation\tLocation\n' >"$_work/baseline/baseline_output.txt"
cp "$_work/baseline/baseline_output.txt" "$_work/builder/baseline_output.txt"
_cache_args="--perl-dir $_work/baseline_out --rust-dir $_work/builder_out --report-dir $_work/concordance_report --summary-json concordance.json"
# $1 is the repository root, $2 what the caller exported; runs the block and the comparison
# step; prints stderr and the step's own lines, then `rc=`.
run_cache_compare() {
    : >"$LOG"
    rm -rf "$_work/baseline_out" "$_work/builder_out" "$_work/concordance_report"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$1"
        WORK_DIR="$_work"
        export_comparators "$2"
        eval "$_cache_block"
        eval "$_step"
        [[ "$CONCORDANCE_REPORT" == "$_work/concordance_report/concordance.json" ]] || echo "report path: $CONCORDANCE_REPORT"
    ) 2>&1 | grep -v -e '^=== Step 4' -e '^$'
    echo "rc=${PIPESTATUS[0]}"
}
# $1 is the repository root: the one line a tree whose build cannot run is refused with.
_cache_refusal() {
    printf '%s' "ERROR: [compare_cache_outputs] the tuple comparator $1/$TOOL_REL is not an executable: build the tool with 'cargo build --release --locked --manifest-path $1/tools/vep-outputs-compare/Cargo.toml'"
}

_out=$(run_cache_compare "$TMP/built" "")
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "tool $_cache_args" && "$(readlink "$_work/baseline_out/output.txt")" == "$_work/baseline/baseline_output.txt" && "$(readlink "$_work/builder_out/output.txt")" == "$_work/builder/baseline_output.txt" ]] &&
    ok "compare_cache_outputs.sh: the tree's release build of the Rust tool compares the two outputs" ||
    bad "compare_cache_outputs.sh: expected the tool alone on the linked outputs, got '$_out' and '$(cat "$LOG")'"
for _i in "${!EXPORTED[@]}"; do
    _out=$(run_cache_compare "$TMP/built" "${EXPORTED[$_i]}")
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "tool $_cache_args" ]] &&
        ok "compare_cache_outputs.sh, both comparator variables exported as ${EXPORTED_AS[$_i]}: the tree's release build compares the two outputs, and nothing else starts" ||
        bad "compare_cache_outputs.sh with ${EXPORTED[$_i]} exported: expected the tool alone, got '$_out' and '$(cat "$LOG")'"
done
_out=$(run_cache_compare "$TMP/unbuilt" "")
[[ "$_out" == "$(_cache_refusal "$TMP/unbuilt")"$'\n'"rc=1" && ! -s "$LOG" && ! -e "$_work/concordance_report" ]] &&
    ok "compare_cache_outputs.sh, the tool not built: one line names the build command, exit 1, nothing is compared and nothing runs in the tool's place" ||
    bad "compare_cache_outputs.sh in an unbuilt tree: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
for _i in 1 2 3; do
    _out=$(run_cache_compare "$TMP/unbuilt" "${EXPORTED[$_i]}")
    [[ "$_out" == "$(_cache_refusal "$TMP/unbuilt")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
        ok "compare_cache_outputs.sh, the tool not built and both comparator variables exported as ${EXPORTED_AS[$_i]}: the same refusal, and what the variables name does not start" ||
        bad "compare_cache_outputs.sh in an unbuilt tree with ${EXPORTED[$_i]} exported: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
done
_out=$(run_cache_compare "$TMP/inert_build" "")
[[ "$_out" == "$(_cache_refusal "$TMP/inert_build")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
    ok "compare_cache_outputs.sh, the tool's release path a file that is not executable: refused with the build command, exit 1, nothing runs" ||
    bad "compare_cache_outputs.sh with a non-executable build: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
_out=$(run_cache_compare "$TMP/dir_build" "")
[[ "$_out" == "$(_cache_refusal "$TMP/dir_build")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
    ok "compare_cache_outputs.sh, the tool's release path a directory: refused with the build command, exit 1, nothing runs" ||
    bad "compare_cache_outputs.sh with a directory at the build's path: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
# The whole script, not the block: a run in a tree without the tool stops on the refusal
# before it starts vep-rs.
mkdir -p "$TMP/tree/scripts/validation"
cp "$CACHE_OUTPUTS" "$TMP/tree/scripts/validation/"
: >"$LOG"
_out=$(PATH="$TMP/bin:$PATH" bash "$TMP/tree/scripts/validation/compare_cache_outputs.sh" --baseline-cache "$TMP/perl_cache" --builder-cache "$TMP/perl_cache" --vcf "$TMP/input.vcf" --vep-binary "$TMP/bin/vep" 2>&1)
rc=$?
[[ "$rc" -eq 1 && "$_out" == "$(_cache_refusal "$TMP/tree")" && ! -s "$LOG" ]] &&
    ok "a compare_cache_outputs.sh run without the tool: exit 1 on the one refusal line, vep-rs never started" ||
    bad "compare_cache_outputs.sh whole without the tool: expected rc=1, the refusal alone and nothing started, got rc=$rc '$_out' and '$(cat "$LOG")'"
# The block and the step on this tree's own release build: the binary takes the step's
# command line over two outputs of three rows that differ in one row's Consequence, and
# writes the report the gates of Step 5 read.
_real_root="$(cd "$HERE/../.." && pwd)"
if [[ -x "$_real_root/$TOOL_REL" ]]; then
    _header='#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n'
    _rows='v1\t21:1000\tG\tENSG01\tENST01\tTranscript\tmissense_variant\t10\t10\t4\tK/R\taAg/aGg\t-\tIMPACT=MODERATE\nv2\t21:2000\tT\tENSG01\tENST01\tTranscript\t%s\t20\t20\t7\tL\tctC/ctT\t-\tIMPACT=LOW\nv3\t21:3000\tA\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n'
    # The row format carries the one Consequence that differs between the two outputs.
    # shellcheck disable=SC2059
    {
        printf "$_header" && printf "$_rows" synonymous_variant
    } >"$_work/baseline/baseline_output.txt"
    # shellcheck disable=SC2059
    {
        printf "$_header" && printf "$_rows" stop_gained
    } >"$_work/builder/baseline_output.txt"
    _out=$(run_cache_compare "$_real_root" "")
    if [[ "$_out" == *$'\n'"rc=0" && ! -s "$LOG" ]] && python3 - "$_work/concordance_report/concordance.json" <<'PY'; then
import json, sys
aggregate = json.load(open(sys.argv[1]))["aggregate"]
counts = (aggregate["perl_tuple_count"], aggregate["rust_tuple_count"], aggregate["intersection_count"])
assert counts == (3, 3, 2), counts
PY
        ok "compare_cache_outputs.sh on the tree's own release build: the comparison step writes concordance_report/concordance.json, three tuples a side, two shared"
    else
        bad "compare_cache_outputs.sh on the tree's own build: expected rc=0, no stub started and 3, 3 and 2 tuples in concordance.json, got '$_out' and '$(cat "$LOG")'"
    fi
else
    echo "  skip the comparison step on the tree's own build: no release build of tools/vep-outputs-compare"
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
