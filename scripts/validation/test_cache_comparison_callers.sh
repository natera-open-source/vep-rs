#!/usr/bin/env bash
# Tests that compare_cache_outputs.sh and validate_builder_cache.sh compare two caches with the
# cache comparison of their own tree, tools/vep-cache-diff, and with nothing else. Each
# script's block for the tool and its Step 3 are extracted and run with a stub cargo, a stub
# python3 and a stub tool that record the command line they were started with;
# compare_cache_outputs.sh is also run whole where the tool cannot run, and its Step 3 runs
# once on this tree's own release build, a check that is skipped where the tool is not built.
#
# The cases:
#   TREE      the tool of the tree runs: its release build by compare_cache_outputs.sh, built
#             and run through cargo by validate_builder_cache.sh.
#   NO-PYTHON no line of either script names a Python file of the scripts directory.
#   REFUSAL   a tool that cannot run (not built, not an executable file, a build that fails)
#             stops the script with one line naming the build command, and nothing runs in
#             its place.
#
# Run: bash scripts/validation/test_cache_comparison_callers.sh

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
CACHE_OUTPUTS="$HERE/compare_cache_outputs.sh"
VALIDATE_BUILDER="$HERE/validate_builder_cache.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

TOOL_REL="tools/vep-cache-diff/target/release/vep-cache-diff"
OUTPUTS_REL="tools/vep-outputs-compare/target/release/vep-outputs-compare"
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
stub "$TMP/bin/vep" vep
stub "$TMP/built/$TOOL_REL" tool
: >"$TMP/inert_build/$TOOL_REL"

# $1 is a script; prints how many of its lines name a Python file of the scripts directory.
names_python() {
    local n
    n=$(grep -c -E -e 'SCRIPT_DIR\}?/[A-Za-z0-9_/]+\.py' "$1") || n=0
    echo "$n"
}

echo "CACHE-OUTPUTS: compare_cache_outputs.sh finds the tree's release build before vep-rs runs and compares the caches through it"

_block=$(awk '
    !p && /^CACHE_DIFF=\(\)$/ { p = 1 }
    p { print }
    p && /^fi$/ { exit }
' "$CACHE_OUTPUTS")
[[ "$_block" == *"cache_diff_bin=\"\${VEP_RS_DIR}/$TOOL_REL\""* && "$_block" == *$'\nfi' ]] &&
    ok "compare_cache_outputs.sh: the tool's block is extracted whole, the tree's release build its one candidate" ||
    bad "compare_cache_outputs.sh: could not extract the tool's block"
_step=$(awk '
    /^# Step 4:/ { exit }
    /^# Step 3:/ { p = 1 }
    p { print }
' "$CACHE_OUTPUTS")
[[ "$_step" == *'"${CACHE_DIFF[@]}" \'* && "$_step" != *python3* ]] &&
    ok "compare_cache_outputs.sh: Step 3 runs the tree's command and names no interpreter" ||
    bad "compare_cache_outputs.sh: Step 3 must run CACHE_DIFF alone"
_named=$(names_python "$CACHE_OUTPUTS")
[[ "$_named" -eq 0 ]] &&
    ok "compare_cache_outputs.sh: no line names a Python file of the scripts directory" ||
    bad "compare_cache_outputs.sh: $_named line(s) name a Python file of the scripts directory"
_block_line=$(grep -n -e '^CACHE_DIFF=()$' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^# Step 1:' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "compare_cache_outputs.sh: the tool is found before the first vep-rs run (line $_block_line before $_first_run)" ||
    bad "compare_cache_outputs.sh: the tool's block must stand before Step 1 (block at ${_block_line:-?}, step at ${_first_run:-?})"

_work="$TMP/cache_work"
mkdir -p "$_work"
_args="--baseline-dir $TMP/baseline_cache --builder-dir $TMP/builder_cache --output-report $_work/cache_diff.json"
# $1 is the repository root; runs the block and Step 3; prints stderr and the step's own
# lines, then `rc=`.
run_cache_step() {
    : >"$LOG"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$1"
        WORK_DIR="$_work"
        BASELINE_CACHE="$TMP/baseline_cache"
        BUILDER_CACHE="$TMP/builder_cache"
        eval "$_block"
        eval "$_step"
    ) 2>&1 | grep -v -e '^=== Step 3' -e '^$'
    echo "rc=${PIPESTATUS[0]}"
}
# $1 is the repository root: the one line a tree whose build cannot run is refused with.
_refusal() {
    printf '%s' "ERROR: [compare_cache_outputs] the cache comparison $1/$TOOL_REL is not an executable: build the tool with 'cargo build --release --locked --manifest-path $1/tools/vep-cache-diff/Cargo.toml'"
}

_out=$(run_cache_step "$TMP/built")
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "tool $_args" ]] &&
    ok "compare_cache_outputs.sh: the tree's release build compares the two caches into cache_diff.json" ||
    bad "compare_cache_outputs.sh: expected the tool alone on the two caches, got '$_out' and '$(cat "$LOG")'"
for _tree in unbuilt inert_build dir_build; do
    _out=$(run_cache_step "$TMP/$_tree")
    [[ "$_out" == "$(_refusal "$TMP/$_tree")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
        ok "compare_cache_outputs.sh, a tree whose tool is $_tree: one line names the build command, exit 1, nothing runs in the tool's place" ||
        bad "compare_cache_outputs.sh in the $_tree tree: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
done
# The whole script, not the block: a tree with the tuple comparator and without this tool
# stops on the one refusal before it starts vep-rs.
mkdir -p "$TMP/tree/scripts/validation" "$TMP/tree/$(dirname "$OUTPUTS_REL")"
cp "$CACHE_OUTPUTS" "$TMP/tree/scripts/validation/"
stub "$TMP/tree/$OUTPUTS_REL" outputs
: >"$LOG"
_out=$(PATH="$TMP/bin:$PATH" bash "$TMP/tree/scripts/validation/compare_cache_outputs.sh" --baseline-cache "$TMP/baseline_cache" --builder-cache "$TMP/builder_cache" --vcf "$TMP/input.vcf" --vep-binary "$TMP/bin/vep" 2>&1)
rc=$?
[[ "$rc" -eq 1 && "$_out" == "$(_refusal "$TMP/tree")" && ! -s "$LOG" ]] &&
    ok "a compare_cache_outputs.sh run without the tool: exit 1 on the one refusal line, vep-rs never started" ||
    bad "compare_cache_outputs.sh whole without the tool: expected rc=1, the refusal alone and nothing started, got rc=$rc '$_out' and '$(cat "$LOG")'"

echo "BUILDER-CACHE: validate_builder_cache.sh builds the tree's tool before a cache is generated and compares through it"

_builder_block=$(awk '
    !p && /^CACHE_DIFF_MANIFEST=/ { p = 1 }
    p { print }
    p && /^CACHE_DIFF=\(/ { exit }
' "$VALIDATE_BUILDER")
[[ "$_builder_block" == *'if ! cargo build -q --release --locked --manifest-path "${CACHE_DIFF_MANIFEST}"; then'* && "$_builder_block" == *$'\nCACHE_DIFF=(cargo run -q --release --locked --manifest-path "${CACHE_DIFF_MANIFEST}" --)' ]] &&
    ok "validate_builder_cache.sh: the tool's block is extracted whole, a cargo build then the cargo run command" ||
    bad "validate_builder_cache.sh: could not extract the tool's block"
_call=$(awk '
    /^# Step 4:/ { exit }
    /^"\$\{CACHE_DIFF\[@\]\}" \\$/ { p = 1 }
    p { print }
' "$VALIDATE_BUILDER")
[[ "$_call" == '"${CACHE_DIFF[@]}" \'* && "$_call" == *'--output-report "${WORK_DIR}/validation_diff.json"'* ]] &&
    ok "validate_builder_cache.sh: Step 3 runs the tree's command into validation_diff.json" ||
    bad "validate_builder_cache.sh: expected one comparison through CACHE_DIFF, got '$_call'"
_named=$(names_python "$VALIDATE_BUILDER")
[[ "$_named" -eq 0 ]] &&
    ok "validate_builder_cache.sh: no line names a Python file of the scripts directory" ||
    bad "validate_builder_cache.sh: $_named line(s) name a Python file of the scripts directory"
_block_line=$(grep -n -e '^CACHE_DIFF_MANIFEST=' "$VALIDATE_BUILDER" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^# Step 2:' "$VALIDATE_BUILDER" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "validate_builder_cache.sh: the tool is built before the cache is generated (line $_block_line before $_first_run)" ||
    bad "validate_builder_cache.sh: the tool's block must stand before Step 2 (block at ${_block_line:-?}, step at ${_first_run:-?})"

_manifest="$TMP/built/tools/vep-cache-diff/Cargo.toml"
_builder_args="--baseline-dir $TMP/baseline_cache --builder-dir $_work/builder_cache --output-report $_work/validation_diff.json"
# Runs the block and the comparison; prints stderr, then `rc=`.
run_builder_step() {
    : >"$LOG"
    # The extracted block and call read the names this subshell sets.
    # shellcheck disable=SC2034
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$TMP/built"
        WORK_DIR="$_work"
        REFERENCE_DIR="$TMP/baseline_cache"
        BUILDER_OUTPUT="$_work/builder_cache"
        eval "$_builder_block"
        eval "$_call"
    ) 2>&1
    echo "rc=$?"
}
_tree_build="cargo build -q --release --locked --manifest-path $_manifest"
_tree_run="cargo run -q --release --locked --manifest-path $_manifest -- $_builder_args"
_out=$(run_builder_step)
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "$_tree_build"$'\n'"$_tree_run" ]] &&
    ok "validate_builder_cache.sh: the tool of the tree is built, then run on the two caches, through cargo" ||
    bad "validate_builder_cache.sh: expected a cargo build then a cargo run of the tool, got '$_out' and '$(cat "$LOG")'"
_out=$(CARGO_BUILD_FAILS=1 run_builder_step)
[[ "$_out" == "ERROR: [validate_builder_cache] the cache comparison tools/vep-cache-diff did not build: repair 'cargo build --release --locked --manifest-path $_manifest'"$'\n'"rc=1" && "$(cat "$LOG")" == "$_tree_build" ]] &&
    ok "validate_builder_cache.sh, the tool's build fails: one line names the build command, exit 1, nothing is compared and nothing runs in the tool's place" ||
    bad "validate_builder_cache.sh with a failing build: expected the refusal, rc=1 and the build alone, got '$_out' and '$(cat "$LOG")'"

# The block and Step 3 of compare_cache_outputs.sh on this tree's own release build: the
# binary takes the step's command line over one fixture cache on both sides and writes the
# report the gates of Step 5 read.
_real_root="$(cd "$HERE/../.." && pwd)"
if [[ -x "$_real_root/$TOOL_REL" ]]; then
    _fixture="$_real_root/tools/vep-cache-diff/tests/fixtures/caches/identical/cache"
    rm -rf "$TMP/baseline_cache" "$TMP/builder_cache"
    cp -r "$_fixture" "$TMP/baseline_cache"
    cp -r "$_fixture" "$TMP/builder_cache"
    _out=$(run_cache_step "$_real_root")
    if [[ "$_out" == *$'\n'"rc=0" && ! -s "$LOG" ]] &&
        grep -q -e '^    "shared": 3,$' "$_work/cache_diff.json" &&
        [[ "$(grep -c -e '^    "match_rate": 1.0' "$_work/cache_diff.json")" -eq 2 ]] &&
        grep -q -e '^    "matches": true,$' "$_work/cache_diff.json"; then
        ok "compare_cache_outputs.sh on the tree's own release build: Step 3 writes cache_diff.json, three transcripts shared, both rates 1.0, info.json matching"
    else
        bad "compare_cache_outputs.sh on the tree's own build: expected rc=0, no stub started and the report of two equal caches, got '$_out' and '$(cat "$LOG")'"
    fi
else
    echo "  skip Step 3 on the tree's own build: no release build of tools/vep-cache-diff"
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
