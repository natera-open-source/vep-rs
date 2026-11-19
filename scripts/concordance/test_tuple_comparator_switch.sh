#!/usr/bin/env bash
# Tests of the tuple comparator switch, VEP_OUTPUTS_COMPARATOR, in run_concordance.sh and
# scripts/validation/compare_cache_outputs.sh. Each script's resolution block and its
# comparator call are extracted and run against stub outputs, with a stub cargo, a stub
# python3 and stub comparators that record the command line they were started with.
#
# The cases, for each script:
#   DEFAULT   unset or empty, the Rust tool of the tree runs: built and run through cargo by
#             run_concordance.sh, the tree's release build by compare_cache_outputs.sh.
#   PYTHON    the value `python`, and no other, runs compare_vep_outputs.py.
#   BINARY    an executable file, or a bare name on PATH, runs as given.
#   REFUSAL   a comparator that cannot run (a build that fails, a tool that is not built, a
#             value that names no executable file) stops the script with one line naming
#             both remedies, before any engine run, and the Python does not run in its place.
#
# Run: bash scripts/concordance/test_tuple_comparator_switch.sh

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
mkdir -p "$TMP/bin" "$TMP/built/$(dirname "$TOOL_REL")" "$TMP/unbuilt"
# Each stub appends its label and arguments, one line a start. The cargo stub fails a
# `build` while CARGO_BUILD_FAILS is set.
stub() {
    printf '#!/bin/bash\nprintf "%%s\\n" "%s $*" >>"%s"\n%s\n' "$2" "$LOG" "${3:-}" >"$1"
    chmod +x "$1"
}
stub "$TMP/bin/cargo" cargo '[[ "$1" == build && -n "${CARGO_BUILD_FAILS:-}" ]] && exit 101; exit 0'
stub "$TMP/bin/python3" python3
stub "$TMP/bin/comparator" comparator
stub "$TMP/built/$TOOL_REL" tool
: >"$TMP/not_executable"
: >"$TMP/bin/inert"

# The value forms the two scripts share: "" leaves the variable unset, `empty` sets it empty.
set_switch() {
    case "$1" in
    "") unset VEP_OUTPUTS_COMPARATOR ;;
    empty) export VEP_OUTPUTS_COMPARATOR="" ;;
    *) export VEP_OUTPUTS_COMPARATOR="$1" ;;
    esac
}

echo "RUN-CONCORDANCE: run_concordance.sh resolves the switch before any engine runs and compares through it"

_block=$(awk '
    !p && /^OUTPUTS_COMPARE_MANIFEST=/ { p = 1 }
    p { print }
    p && /^esac$/ { exit }
' "$CONCORDANCE")
[[ "$_block" == *VEP_OUTPUTS_COMPARATOR* && "$_block" == *$'\nesac' ]] &&
    ok "run_concordance.sh: the resolution block is extracted whole" ||
    bad "run_concordance.sh: could not extract the resolution block"
_call=$(grep -e '^run_timed "compare_outputs" ' "$CONCORDANCE")
[[ "$_call" == 'run_timed "compare_outputs" "${OUTPUTS_COMPARE[@]}" "${COMPARE_ARGS[@]}"' ]] &&
    ok "run_concordance.sh: the one compare step runs the resolved command under its timer" ||
    bad "run_concordance.sh: expected one compare step through OUTPUTS_COMPARE, got '$_call'"
_hard=$(grep -c -e 'python3 .*compare_vep_outputs\.py' "$CONCORDANCE") || _hard=0
[[ "$_hard" -eq 1 ]] &&
    ok "run_concordance.sh: the Python is run from the one line that selects it" ||
    bad "run_concordance.sh: expected 1 line that runs the Python, found $_hard"
_block_line=$(grep -n -e '^OUTPUTS_COMPARE_MANIFEST=' "$CONCORDANCE" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^[[:space:]]*run_timed "' "$CONCORDANCE" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "run_concordance.sh: the comparator is resolved before the first timed step (line $_block_line before $_first_run)" ||
    bad "run_concordance.sh: the resolution block must stand before the first timed step (block at ${_block_line:-?}, step at ${_first_run:-?})"

_args="--perl-dir $TMP/work/perl --rust-dir $TMP/work/rust --report-dir $TMP/work/reports"
_manifest="$TMP/built/tools/vep-outputs-compare/Cargo.toml"
# $1 is the switch value; runs the block and the compare step; prints stderr, then `rc=`.
run_concordance_compare() {
    : >"$LOG"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034,SC2329
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$TMP/built"
        set_switch "$1"
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
_concordance_refusal() {
    printf '%s' "ERROR: [run_concordance] the tuple comparator $1 is not an executable: name a built one in VEP_OUTPUTS_COMPARATOR or unset it to build tools/vep-outputs-compare, or set VEP_OUTPUTS_COMPARATOR=python to score with the Python reference implementation"
}

for _unset in "" empty; do
    _out=$(run_concordance_compare "$_unset")
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "cargo build -q --release --locked --manifest-path $_manifest"$'\n'"cargo run -q --release --locked --manifest-path $_manifest -- $_args" ]] &&
        ok "run_concordance.sh, VEP_OUTPUTS_COMPARATOR ${_unset:-unset}: the Rust tool is built, then run on the outputs, through cargo" ||
        bad "run_concordance.sh ${_unset:-unset}: expected a cargo build then a cargo run of the tool, got '$_out' and '$(cat "$LOG")'"
done
_out=$(CARGO_BUILD_FAILS=1 run_concordance_compare "")
[[ "$_out" == "ERROR: [run_concordance] the tuple comparator tools/vep-outputs-compare did not build: repair 'cargo build --release --locked --manifest-path $_manifest' or name a built one in VEP_OUTPUTS_COMPARATOR, or set VEP_OUTPUTS_COMPARATOR=python to score with the Python reference implementation"$'\n'"rc=1" && "$(cat "$LOG")" == "cargo build -q --release --locked --manifest-path $_manifest" ]] &&
    ok "run_concordance.sh, the tool's build fails: one line names both remedies, exit 1, nothing is compared and the Python does not run" ||
    bad "run_concordance.sh with a failing build: expected the refusal, rc=1 and the build alone, got '$_out' and '$(cat "$LOG")'"
_out=$(run_concordance_compare python)
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "python3 $TMP/built/scripts/concordance/compare_vep_outputs.py $_args" ]] &&
    ok "run_concordance.sh, VEP_OUTPUTS_COMPARATOR=python: the Python compares the outputs and cargo builds nothing" ||
    bad "run_concordance.sh =python: expected the Python alone, got '$_out' and '$(cat "$LOG")'"
for _value in "$TMP/bin/comparator" comparator; do
    _out=$(run_concordance_compare "$_value")
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "comparator $_args" ]] &&
        ok "run_concordance.sh, VEP_OUTPUTS_COMPARATOR=${_value##*/} ($([[ "$_value" == */* ]] && echo a path || echo a bare name on PATH)): the binary compares the outputs, run as given" ||
        bad "run_concordance.sh =$_value: expected the comparator stub alone, got '$_out' and '$(cat "$LOG")'"
done
for _value in "$TMP/not_executable" "$TMP" inert; do
    _out=$(run_concordance_compare "$_value")
    [[ "$_out" == "$(_concordance_refusal "$_value")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
        ok "run_concordance.sh, VEP_OUTPUTS_COMPARATOR=${_value##*/} (not an executable file): refused with both remedies, exit 1, nothing runs" ||
        bad "run_concordance.sh =$_value: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
done

echo "CACHE-OUTPUTS: compare_cache_outputs.sh resolves the switch before vep-rs runs and compares through it"

_cache_block=$(awk '
    !p && /^case "\$\{VEP_OUTPUTS_COMPARATOR:-\}" in$/ { p = 1 }
    p { print }
    p && /^esac$/ { exit }
' "$CACHE_OUTPUTS")
[[ "$_cache_block" == *OUTPUTS_COMPARE=* && "$_cache_block" == *$'\nesac' ]] &&
    ok "compare_cache_outputs.sh: the resolution block is extracted whole" ||
    bad "compare_cache_outputs.sh: could not extract the resolution block"
_step=$(awk '
    /^# Step 5:/ { exit }
    /^# Step 4:/ { p = 1 }
    p { print }
' "$CACHE_OUTPUTS")
[[ "$_step" == *'"${OUTPUTS_COMPARE[@]}" \'* && "$_step" != *python3* ]] &&
    ok "compare_cache_outputs.sh: the comparison step runs the resolved command and names no interpreter" ||
    bad "compare_cache_outputs.sh: the comparison step must run OUTPUTS_COMPARE alone"
_block_line=$(grep -n -e '^case "\${VEP_OUTPUTS_COMPARATOR:-}" in$' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^# Step 1:' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "compare_cache_outputs.sh: the comparator is resolved before the first vep-rs run (line $_block_line before $_first_run)" ||
    bad "compare_cache_outputs.sh: the resolution block must stand before Step 1 (block at ${_block_line:-?}, step at ${_first_run:-?})"

_work="$TMP/cache_work"
mkdir -p "$_work/baseline" "$_work/builder"
printf '#Uploaded_variation\tLocation\n' >"$_work/baseline/baseline_output.txt"
cp "$_work/baseline/baseline_output.txt" "$_work/builder/baseline_output.txt"
_cache_args="--perl-dir $_work/baseline_out --rust-dir $_work/builder_out --report-dir $_work/concordance_report --summary-json concordance.json"
# $1 is the repository root, $2 the switch value; runs the block and the comparison step;
# prints stderr and the step's own lines, then `rc=`.
run_cache_compare() {
    : >"$LOG"
    rm -rf "$_work/baseline_out" "$_work/builder_out" "$_work/concordance_report"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034
    (
        PATH="$TMP/bin:$PATH"
        VEP_RS_DIR="$1"
        WORK_DIR="$_work"
        set_switch "$2"
        eval "$_cache_block"
        eval "$_step"
        [[ "$CONCORDANCE_REPORT" == "$_work/concordance_report/concordance.json" ]] || echo "report path: $CONCORDANCE_REPORT"
    ) 2>&1 | grep -v -e '^=== Step 4' -e '^$'
    echo "rc=${PIPESTATUS[0]}"
}
_cache_refusal() {
    printf '%s' "ERROR: [compare_cache_outputs] the tuple comparator $1 is not an executable: build the tool with 'cargo build --release --locked --manifest-path $2/tools/vep-outputs-compare/Cargo.toml' or name a built one in VEP_OUTPUTS_COMPARATOR, or set VEP_OUTPUTS_COMPARATOR=python to score with the Python reference implementation"
}

for _unset in "" empty; do
    _out=$(run_cache_compare "$TMP/built" "$_unset")
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "tool $_cache_args" && "$(readlink "$_work/baseline_out/output.txt")" == "$_work/baseline/baseline_output.txt" && "$(readlink "$_work/builder_out/output.txt")" == "$_work/builder/baseline_output.txt" ]] &&
        ok "compare_cache_outputs.sh, VEP_OUTPUTS_COMPARATOR ${_unset:-unset}: the tree's release build of the Rust tool compares the two outputs" ||
        bad "compare_cache_outputs.sh ${_unset:-unset}: expected the tool alone on the linked outputs, got '$_out' and '$(cat "$LOG")'"
done
_out=$(run_cache_compare "$TMP/unbuilt" "")
[[ "$_out" == "$(_cache_refusal "$TMP/unbuilt/$TOOL_REL" "$TMP/unbuilt")"$'\n'"rc=1" && ! -s "$LOG" && ! -e "$_work/concordance_report" ]] &&
    ok "compare_cache_outputs.sh, the tool not built: one line names both remedies, exit 1, nothing is compared and the Python does not run" ||
    bad "compare_cache_outputs.sh in an unbuilt tree: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
_out=$(run_cache_compare "$TMP/unbuilt" python)
[[ "$_out" == "rc=0" && "$(cat "$LOG")" == "python3 $TMP/unbuilt/scripts/concordance/compare_vep_outputs.py $_cache_args" ]] &&
    ok "compare_cache_outputs.sh, VEP_OUTPUTS_COMPARATOR=python: the Python compares the two outputs, the tool unbuilt" ||
    bad "compare_cache_outputs.sh =python: expected the Python alone, got '$_out' and '$(cat "$LOG")'"
for _value in "$TMP/bin/comparator" comparator; do
    _out=$(run_cache_compare "$TMP/unbuilt" "$_value")
    [[ "$_out" == "rc=0" && "$(cat "$LOG")" == "comparator $_cache_args" ]] &&
        ok "compare_cache_outputs.sh, VEP_OUTPUTS_COMPARATOR=${_value##*/} ($([[ "$_value" == */* ]] && echo a path || echo a bare name on PATH)): the binary compares the two outputs, run as given" ||
        bad "compare_cache_outputs.sh =$_value: expected the comparator stub alone, got '$_out' and '$(cat "$LOG")'"
done
for _value in "$TMP/not_executable" "$TMP" inert; do
    _out=$(run_cache_compare "$TMP/built" "$_value")
    [[ "$_out" == "$(_cache_refusal "$_value" "$TMP/built")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
        ok "compare_cache_outputs.sh, VEP_OUTPUTS_COMPARATOR=${_value##*/} (not an executable file): refused with both remedies, exit 1, nothing runs" ||
        bad "compare_cache_outputs.sh =$_value: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
done

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
