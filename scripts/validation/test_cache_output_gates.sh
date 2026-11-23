#!/usr/bin/env bash
# Tests that compare_cache_outputs.sh scores its gates with the gate scorer of its own tree,
# tools/vep-cache-output-gates, and with nothing else. The script's block for the tool and its
# Step 5 are extracted and run with a stub tool that records the command line it was started
# with and leaves with a chosen status; the script is run whole on stub tools and a stub
# capture step, and whole where the tool cannot run; and its Step 5 runs on this tree's own
# release build over three cases of the tool's fixtures, checks that are skipped where the
# tool is not built.
#
# The cases:
#   TREE      the tool of the tree runs, its release build, on ten operands in the tool's
#             order; its exit code is the verdict the script prints and leaves with.
#   NO-PYTHON no line of the script names an interpreter.
#   REFUSAL   a tool that cannot run (not built, not an executable file) stops the script
#             with one line naming the build command, before vep-rs runs.
#   SUMMARY   a tuple comparison that fails is scored on a summary that is missing, never
#             on the summary an earlier run left in the work directory.
#
# Run: bash scripts/validation/test_cache_output_gates.sh

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
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

TOOL_REL="tools/vep-cache-output-gates/target/release/vep-cache-output-gates"
OUTPUTS_REL="tools/vep-outputs-compare/target/release/vep-outputs-compare"
CACHE_DIFF_REL="tools/vep-cache-diff/target/release/vep-cache-diff"
LOG="$TMP/started"
mkdir -p "$TMP/bin" "$TMP/built/$(dirname "$TOOL_REL")" "$TMP/unbuilt" "$TMP/inert_build/$(dirname "$TOOL_REL")" "$TMP/dir_build/$TOOL_REL"
# Each stub appends its label and arguments, one line a start, then runs its third argument.
stub() {
    printf '#!/bin/bash\nprintf "%%s\\n" "%s $*" >>"%s"\n%s\n' "$2" "$LOG" "${3:-}" >"$1"
    chmod +x "$1"
}
# The stub tool leaves with GATES_EXIT, as the tool leaves with its verdict.
stub "$TMP/bin/vep" vep
stub "$TMP/built/$TOOL_REL" gates 'exit "${GATES_EXIT:-0}"'
: >"$TMP/inert_build/$TOOL_REL"

BANNER="============================================"

echo "CACHE-OUTPUT-GATES: compare_cache_outputs.sh finds the tree's release build before vep-rs runs and scores its gates through it"

_block=$(awk '
    !p && /^OUTPUT_GATES=\(\)$/ { p = 1 }
    p { print }
    p && /^fi$/ { exit }
' "$CACHE_OUTPUTS")
[[ "$_block" == *"output_gates_bin=\"\${VEP_RS_DIR}/$TOOL_REL\""* && "$_block" == *$'\nfi' ]] &&
    ok "compare_cache_outputs.sh: the tool's block is extracted whole, the tree's release build its one candidate" ||
    bad "compare_cache_outputs.sh: could not extract the tool's block"
_step=$(awk '
    /^# Step 5:/ { p = 1 }
    p { print }
' "$CACHE_OUTPUTS")
[[ "$_step" == *'"${OUTPUT_GATES[@]}" \'* && "$_step" == *$'\nfi' ]] &&
    ok "compare_cache_outputs.sh: Step 5 runs the tree's command and is the script's last step" ||
    bad "compare_cache_outputs.sh: Step 5 must run OUTPUT_GATES"
_named=$(grep -c -i -e 'python' "$CACHE_OUTPUTS") || _named=0
[[ "$_named" -eq 0 ]] &&
    ok "compare_cache_outputs.sh: no line names an interpreter" ||
    bad "compare_cache_outputs.sh: $_named line(s) name an interpreter"
_block_line=$(grep -n -e '^OUTPUT_GATES=()$' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
_first_run=$(grep -n -e '^# Step 1:' "$CACHE_OUTPUTS" | head -1 | cut -d: -f1)
[[ -n "$_block_line" && -n "$_first_run" && "$_block_line" -lt "$_first_run" ]] &&
    ok "compare_cache_outputs.sh: the tool is found before the first vep-rs run (line $_block_line before $_first_run)" ||
    bad "compare_cache_outputs.sh: the tool's block must stand before Step 1 (block at ${_block_line:-?}, step at ${_first_run:-?})"

_work="$TMP/gate_work"
mkdir -p "$_work"
_operands="$_work/baseline/baseline_metrics.json $_work/builder/baseline_metrics.json $_work/cache_diff.json $_work/concordance_report/concordance.json $_work/verdict.json"
# $1 is the repository root; runs the block and Step 5 with the script's default gates;
# prints both streams, then `rc=`.
run_gate_step() {
    : >"$LOG"
    # The extracted block and step read the names this subshell sets.
    # shellcheck disable=SC2034
    (
        VEP_RS_DIR="$1"
        WORK_DIR="$_work"
        CONCORDANCE_REPORT="$_work/concordance_report/concordance.json"
        GATE_F1=0.999
        GATE_TOTAL_DELTA=0.001
        GATE_COLOC_DELTA=0.01
        GATE_FREQ_DELTA=0.01
        GATE_PER_CSQ_DELTA=0.01
        eval "$_block"
        eval "$_step"
    ) 2>&1
    echo "rc=$?"
}
# $1 is the verdict: what Step 5 prints around a tool that prints nothing, with its status.
_verdict_lines() {
    printf '%s\n' "=== Step 5: Validation gates ===" "" "$BANNER" "  VERDICT: $1" "$BANNER"
}
# $1 is the repository root: the one line a tree whose build cannot run is refused with.
_refusal() {
    printf '%s' "ERROR: [compare_cache_outputs] the gate scorer $1/$TOOL_REL is not an executable: build the tool with 'cargo build --release --locked --manifest-path $1/tools/vep-cache-output-gates/Cargo.toml'"
}

_out=$(run_gate_step "$TMP/built")
[[ "$_out" == "$(_verdict_lines PASS)"$'\n'"rc=0" && "$(cat "$LOG")" == "gates $_operands 0.999 0.001 0.01 0.01 0.01" ]] &&
    ok "compare_cache_outputs.sh, the tool leaves with 0: ten operands in the tool's order, VERDICT: PASS, exit 0" ||
    bad "compare_cache_outputs.sh with a passing tool: expected the tool on its ten operands, the PASS lines and rc=0, got '$_out' and '$(cat "$LOG")'"
_out=$(GATES_EXIT=1 run_gate_step "$TMP/built")
[[ "$_out" == "$(_verdict_lines FAIL)"$'\n'"rc=1" && "$(cat "$LOG")" == "gates $_operands 0.999 0.001 0.01 0.01 0.01" ]] &&
    ok "compare_cache_outputs.sh, the tool leaves with 1: VERDICT: FAIL, exit 1" ||
    bad "compare_cache_outputs.sh with a failing tool: expected the FAIL lines and rc=1, got '$_out' and '$(cat "$LOG")'"
for _status in 2 101 134; do
    _out=$(GATES_EXIT=$_status run_gate_step "$TMP/built")
    [[ "$_out" == "=== Step 5: Validation gates ==="$'\n'"rc=1" ]] &&
        ok "compare_cache_outputs.sh, the tool leaves with $_status: no verdict is printed, exit 1" ||
        bad "compare_cache_outputs.sh with a tool that leaves with $_status: expected the heading alone and rc=1, got '$_out'"
done
for _tree in unbuilt inert_build dir_build; do
    _out=$(run_gate_step "$TMP/$_tree")
    [[ "$_out" == "$(_refusal "$TMP/$_tree")"$'\n'"rc=1" && ! -s "$LOG" ]] &&
        ok "compare_cache_outputs.sh, a tree whose tool is $_tree: one line names the build command, exit 1, nothing runs in the tool's place" ||
        bad "compare_cache_outputs.sh in the $_tree tree: expected the refusal, rc=1 and nothing started, got '$_out' and '$(cat "$LOG")'"
done

# The whole script on stub tools and a stub capture step: its five gate options reach the
# tool in the tool's order, behind the five reports of its work directory.
mkdir -p "$TMP/whole/scripts/validation" "$TMP/whole/$(dirname "$OUTPUTS_REL")" "$TMP/whole/$(dirname "$CACHE_DIFF_REL")" "$TMP/whole/$(dirname "$TOOL_REL")" "$TMP/cache_a" "$TMP/cache_b"
cp "$CACHE_OUTPUTS" "$TMP/whole/scripts/validation/"
stub "$TMP/whole/scripts/validation/capture_baseline.sh" capture
stub "$TMP/whole/$OUTPUTS_REL" outputs
stub "$TMP/whole/$CACHE_DIFF_REL" cache_diff
stub "$TMP/whole/$TOOL_REL" gates 'exit "${GATES_EXIT:-0}"'
: >"$TMP/input.vcf"
_whole_work="$TMP/whole_work"
_whole_operands="$_whole_work/baseline/baseline_metrics.json $_whole_work/builder/baseline_metrics.json $_whole_work/cache_diff.json $_whole_work/concordance_report/concordance.json $_whole_work/verdict.json"
# Runs the whole script with five gate values of its own; prints stdout and stderr.
run_whole() {
    : >"$LOG"
    bash "$TMP/whole/scripts/validation/compare_cache_outputs.sh" \
        --baseline-cache "$TMP/cache_a" --builder-cache "$TMP/cache_b" --vcf "$TMP/input.vcf" --vep-binary "$TMP/bin/vep" \
        --work-dir "$_whole_work" \
        --gate-f1 0.5 --gate-total-delta 0.25 --gate-coloc-delta 1e-3 --gate-freq-delta 2 --gate-per-csq-delta 0.125 2>&1
}
_out=$(run_whole)
rc=$?
[[ "$rc" -eq 0 && "$_out" == *$'\n'"$(_verdict_lines PASS)" && "$(grep -e '^gates ' "$LOG")" == "gates $_whole_operands 0.5 0.25 1e-3 2 0.125" && "$(grep -c -e '^capture ' "$LOG")" -eq 2 && "$(tail -1 "$LOG")" == gates\ * ]] &&
    ok "a whole compare_cache_outputs.sh run: the five --gate options reach the tool in its order, last of all, VERDICT: PASS, exit 0" ||
    bad "compare_cache_outputs.sh whole on stubs: expected rc=0, the PASS lines last and the tool on its ten operands, got rc=$rc '$_out' and '$(cat "$LOG")'"
_out=$(GATES_EXIT=1 run_whole)
rc=$?
[[ "$rc" -eq 1 && "$_out" == *$'\n'"$(_verdict_lines FAIL)" ]] &&
    ok "a whole compare_cache_outputs.sh run whose tool leaves with 1: VERDICT: FAIL, exit 1" ||
    bad "compare_cache_outputs.sh whole with a failing tool: expected rc=1 and the FAIL lines last, got rc=$rc '$_out'"
_out=$(GATES_EXIT=2 run_whole)
rc=$?
[[ "$rc" -eq 1 && "$_out" == *$'\n'"=== Step 5: Validation gates ===" ]] &&
    ok "a whole compare_cache_outputs.sh run whose tool leaves with 2: exit 1 after the step's heading, no verdict printed" ||
    bad "compare_cache_outputs.sh whole with a tool that reached no verdict: expected rc=1 and the heading last, got rc=$rc '$_out'"
# A tuple comparison that fails writes no summary, and the summary of an earlier run in the
# same work directory must not be scored in its place: the tool starts on a path that names
# nothing. The comparator stub fails here and the tool stub records a summary it finds.
mkdir -p "$_whole_work/concordance_report"
printf '{"aggregate": {"f1": 1.0}}\n' >"$_whole_work/concordance_report/concordance.json"
stub "$TMP/whole/$OUTPUTS_REL" outputs 'exit 1'
stub "$TMP/whole/$TOOL_REL" gates '[[ -e "$4" ]] && echo "summary found" >>"'"$LOG"'"; exit 1'
_out=$(run_whole)
rc=$?
[[ "$rc" -eq 1 && "$_out" == *$'\n'"$(_verdict_lines FAIL)" && "$(grep -c -e '^outputs ' "$LOG")" -eq 1 && "$(grep -e '^gates ' "$LOG")" == "gates $_whole_operands 0.5 0.25 1e-3 2 0.125" && "$(grep -c -e '^summary found$' "$LOG")" -eq 0 && ! -e "$_whole_work/concordance_report/concordance.json" ]] &&
    ok "a whole compare_cache_outputs.sh run whose tuple comparison fails over an earlier run's summary: the summary is gone when the tool starts, VERDICT: FAIL, exit 1" ||
    bad "compare_cache_outputs.sh whole with a failing comparison over an earlier summary: expected rc=1, the FAIL lines last and the tool started on no summary, got rc=$rc '$_out' and '$(cat "$LOG")'"
# The whole script in a tree with the two comparison tools and without this one stops on
# the one refusal before it starts vep-rs.
rm "$TMP/whole/$TOOL_REL"
_out=$(run_whole)
rc=$?
[[ "$rc" -eq 1 && "$_out" == "$(_refusal "$TMP/whole")" && ! -s "$LOG" ]] &&
    ok "a compare_cache_outputs.sh run without the tool: exit 1 on the one refusal line, vep-rs never started" ||
    bad "compare_cache_outputs.sh whole without the tool: expected rc=1, the refusal alone and nothing started, got rc=$rc '$_out' and '$(cat "$LOG")'"

# The block and Step 5 on this tree's own release build, over cases of the tool's fixtures:
# the lines the tool's golden holds, the script's verdict lines behind them, the verdict
# report of the golden, and the exit code of the case.
_real_root="$(cd "$HERE/../.." && pwd)"
_fixtures="$_real_root/tools/vep-cache-output-gates/tests/fixtures"
# $1 is a case of cases.tsv: lays its four reports out as the script's work directory holds
# them and prints the case's exit code.
lay_case() {
    local row
    row=$(awk -F'\t' -v name="$1" '$1 == name' "$_fixtures/cases.tsv")
    rm -rf "$_work"
    mkdir -p "$_work/baseline" "$_work/builder" "$_work/concordance_report"
    cp "$_fixtures/reports/$(cut -f2 <<<"$row")" "$_work/baseline/baseline_metrics.json"
    cp "$_fixtures/reports/$(cut -f3 <<<"$row")" "$_work/builder/baseline_metrics.json"
    cp "$_fixtures/reports/$(cut -f4 <<<"$row")" "$_work/cache_diff.json"
    cp "$_fixtures/reports/$(cut -f5 <<<"$row")" "$_work/concordance_report/concordance.json"
    cut -f11 <<<"$row"
}
# $1 is a case: the tool's lines of its golden, the report named as the script names it.
golden_lines() {
    sed -e "s|^  Report:  verdict.json\$|  Report:  $_work/verdict.json|" "$_fixtures/expected/$1.stdout"
}
if [[ -x "$_real_root/$TOOL_REL" ]]; then
    _status=$(lay_case equal_pass)
    _out=$(run_gate_step "$_real_root")
    [[ "$_status" == 0 && "$_out" == "=== Step 5: Validation gates ==="$'\n'"$(golden_lines equal_pass)"$'\n\n'"$BANNER"$'\n'"  VERDICT: PASS"$'\n'"$BANNER"$'\n'"rc=0" ]] &&
        cmp -s "$_work/verdict.json" "$_fixtures/expected/equal_pass.verdict.json" &&
        ok "compare_cache_outputs.sh on the tree's own release build, equal metrics: the golden's lines, VERDICT: PASS, the golden's verdict.json, exit 0" ||
        bad "compare_cache_outputs.sh on the tree's own build over equal_pass: got '$_out'"
    _status=$(lay_case far_fail)
    _out=$(run_gate_step "$_real_root")
    [[ "$_status" == 1 && "$_out" == "=== Step 5: Validation gates ==="$'\n'"$(golden_lines far_fail)"$'\n\n'"$BANNER"$'\n'"  VERDICT: FAIL"$'\n'"$BANNER"$'\n'"rc=1" ]] &&
        cmp -s "$_work/verdict.json" "$_fixtures/expected/far_fail.verdict.json" &&
        ok "compare_cache_outputs.sh on the tree's own release build, metrics apart: the golden's lines, VERDICT: FAIL, the golden's verdict.json, exit 1" ||
        bad "compare_cache_outputs.sh on the tree's own build over far_fail: got '$_out'"
    _status=$(lay_case cache_diff_list)
    _out=$(run_gate_step "$_real_root" | grep -v -e '^ERROR: \[vep-cache-output-gates\] ')
    [[ "$_status" == 2 && "$_out" == "=== Step 5: Validation gates ==="$'\n'"$(cat "$_fixtures/expected/cache_diff_list.stdout")"$'\n'"rc=1" && ! -e "$_work/verdict.json" ]] &&
        ok "compare_cache_outputs.sh on the tree's own release build, a cache comparison report that is no object: the five gate lines, no verdict printed or written, exit 1" ||
        bad "compare_cache_outputs.sh on the tree's own build over cache_diff_list: got '$_out'"
else
    echo "  skip Step 5 on the tree's own build: no release build of tools/vep-cache-output-gates"
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
