#!/usr/bin/env bash
# Tests of the plugin comparator switch, VEP_PLUGIN_COMPARATOR, in run_concordance.sh. The
# script's resolution block and its plugin comparison step are extracted and run with a stub
# cargo, a stub python3 and a stub comparator that record the command line they were started
# with.
#
# The cases:
#   DEFAULT   unset or empty under --plugins, the Rust tool of the tree, tools/vep-plugin-pairs,
#             is built and its `compare` command run, both through cargo.
#   NO-PLUGINS  a run without --plugins builds and resolves nothing.
#   BINARY    an executable file, or a bare name on PATH, runs its `compare` command as given.
#   REFUSAL   a comparator that cannot run (a build that fails, a value that names no
#             executable file) stops the script with one line naming the remedy, before any
#             engine run.
#   No case starts an interpreter: the step has no other program to run.
#
# Run: bash scripts/concordance/test_plugin_comparator_switch.sh

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
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

LOG="$TMP/started"
mkdir -p "$TMP/bin" "$TMP/tree"
# Each stub appends its label and arguments, one line a start. The cargo stub fails a
# `build` while CARGO_BUILD_FAILS is set.
stub() {
    printf '#!/bin/bash\nprintf "%%s\\n" "%s $*" >>"%s"\n%s\n' "$2" "$LOG" "${3:-}" >"$1"
    chmod +x "$1"
}
stub "$TMP/bin/cargo" cargo '[[ "$1" == build && -n "${CARGO_BUILD_FAILS:-}" ]] && exit 101; exit 0'
stub "$TMP/bin/python3" python3
stub "$TMP/bin/comparator" comparator
: >"$TMP/not_executable"
: >"$TMP/bin/inert"

echo "RUN-CONCORDANCE: run_concordance.sh resolves the plugin comparator before any engine runs and compares through it"

_block=$(awk '
    !p && /^PLUGIN_COMPARE_MANIFEST=/ { p = 1 }
    p { print }
    p && /^fi$/ { exit }
' "$CONCORDANCE")
[[ "$_block" == *VEP_PLUGIN_COMPARATOR* && "$_block" == *$'\nfi' ]] &&
    ok "the resolution block is extracted whole" ||
    bad "could not extract the resolution block"
_call=$(awk '
    !p && /"\$\{PLUGIN_COMPARE\[@\]\}"/ { p = 1 }
    p { print }
    p && !/\\$/ { exit }
' "$CONCORDANCE")
_runs=$(grep -c -e '"\${PLUGIN_COMPARE\[@\]}"' "$CONCORDANCE") || _runs=0
[[ "$_runs" -eq 1 && "$_call" == *'run_timed "'*'"${PLUGIN_COMPARE[@]}" \'* && "$_call" == *'|| PLUGIN_COMPARE_STATUS=$?' ]] &&
    ok "the one plugin comparison step runs the resolved command under its timer and keeps its exit status" ||
    bad "expected one plugin comparison step through PLUGIN_COMPARE, found $_runs: '$_call'"
[[ "$_call" != *python3* && "$_block" != *python* ]] &&
    ok "neither the resolution block nor the step names an interpreter" ||
    bad "the plugin comparison must run the resolved command alone"
_block_line=$(grep -n -e '^PLUGIN_COMPARE_MANIFEST=' "$CONCORDANCE" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^[[:space:]]*run_timed "' "$CONCORDANCE" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "the comparator is resolved before the first timed step (line $_block_line before $_first_run)" ||
    bad "the resolution block must stand before the first timed step (block at ${_block_line:-?}, step at ${_first_run:-?})"

_args="--perl-dir $TMP/work/perl --rust-dir $TMP/work/rust --plugins REVEL,LoFtool --min-annotated 1 --json-output $TMP/work/reports/plugin_concordance.json"
_manifest="$TMP/tree/tools/vep-plugin-pairs/Cargo.toml"
# $1 is the switch value ("" leaves the variable unset, `empty` sets it empty), $2 the value
# of --plugins; runs the block and, under plugins, the comparison step; prints stderr, the
# step's exit status where it is not 0, then `rc=`.
run_plugin_compare() {
    : >"$LOG"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034,SC2329
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$TMP/tree"
        WORK_DIR="$TMP/work"
        PLUGINS="$2"
        PLUGIN_LIST_CSV="$2"
        PLUGIN_MIN_ANNOTATED=1
        case "$1" in
        "") unset VEP_PLUGIN_COMPARATOR ;;
        empty) export VEP_PLUGIN_COMPARATOR="" ;;
        *) export VEP_PLUGIN_COMPARATOR="$1" ;;
        esac
        run_timed() {
            shift
            "$@"
        }
        eval "$_block"
        if [[ -n "$PLUGINS" ]]; then
            PLUGIN_COMPARE_STATUS=0
            eval "$_call"
            [[ "$PLUGIN_COMPARE_STATUS" -eq 0 ]] || echo "status=$PLUGIN_COMPARE_STATUS"
        fi
    ) 2>&1
    echo "rc=$?"
}
_refusal() {
    printf '%s' "ERROR: [run_concordance] the plugin comparator $1 is not an executable: name a built vep-plugin-pairs in VEP_PLUGIN_COMPARATOR or unset it to build tools/vep-plugin-pairs"
}

for _unset in "" empty; do
    _out=$(run_plugin_compare "$_unset" REVEL,LoFtool)
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "cargo build -q --release --locked --manifest-path $_manifest"$'\n'"cargo run -q --release --locked --manifest-path $_manifest -- compare $_args" ]] &&
        ok "VEP_PLUGIN_COMPARATOR ${_unset:-unset}: the Rust tool is built, then its compare command run on the outputs, through cargo" ||
        bad "${_unset:-unset}: expected a cargo build then a cargo run of the tool's compare command, got '$_out' and '$(cat "$LOG")'"
done
for _value in "" "$TMP/not_executable"; do
    _out=$(run_plugin_compare "$_value" "")
    [[ "$_out" == "rc=0" && ! -s "$LOG" ]] &&
        ok "no --plugins, VEP_PLUGIN_COMPARATOR $([[ -n "$_value" ]] && echo "naming no executable" || echo unset): nothing is built and the value is not read" ||
        bad "no --plugins with '${_value}': expected nothing started and rc=0, got '$_out' and '$(cat "$LOG")'"
done
_out=$(CARGO_BUILD_FAILS=1 run_plugin_compare "" REVEL,LoFtool)
[[ "$_out" == "ERROR: [run_concordance] the plugin comparator tools/vep-plugin-pairs did not build: repair 'cargo build --release --locked --manifest-path $_manifest' or name a built one in VEP_PLUGIN_COMPARATOR"$'\n'"rc=1" && "$(cat "$LOG")" == "cargo build -q --release --locked --manifest-path $_manifest" ]] &&
    ok "the tool's build fails: one line names both remedies, exit 1, nothing is compared" ||
    bad "with a failing build: expected the refusal, rc=1 and the build alone, got '$_out' and '$(cat "$LOG")'"
for _value in "$TMP/bin/comparator" comparator; do
    _out=$(run_plugin_compare "$_value" REVEL,LoFtool)
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "comparator compare $_args" ]] &&
        ok "VEP_PLUGIN_COMPARATOR=${_value##*/} ($([[ "$_value" == */* ]] && echo a path || echo a bare name on PATH)): the binary's compare command runs on the outputs, and cargo builds nothing" ||
        bad "=$_value: expected the comparator stub's compare command alone, got '$_out' and '$(cat "$LOG")'"
done
for _value in "$TMP/not_executable" "$TMP" inert; do
    _out=$(run_plugin_compare "$_value" REVEL,LoFtool)
    [[ "$_out" == "$(_refusal "$_value")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
        ok "VEP_PLUGIN_COMPARATOR=${_value##*/} (not an executable file): refused with both remedies, exit 1, nothing runs" ||
        bad "=$_value: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
done
stub "$TMP/bin/failing" failing 'exit 1'
_out=$(run_plugin_compare failing REVEL,LoFtool)
[[ "$_out" == "status=1"$'\n'"rc=0" && "$(cat "$LOG")" == "failing compare $_args" ]] &&
    ok "a comparison that exits 1 leaves its status for the run's exit, the step's own line not ended by it" ||
    bad "a failing comparison: expected status=1 kept and the step survived, got '$_out' and '$(cat "$LOG")'"

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
