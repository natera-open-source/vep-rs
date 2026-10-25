#!/usr/bin/env bash
# Adversarial tests for the measurement harness's failure-propagation gates.
#
# These exercise the exact loop and fallback SHAPES used by
# run_clone_measurement.sh, with a stub engine and a stub bcftools, so each test
# proves a gate fires on an injected fault rather than asserting the happy path.
#
# The injected faults:
#   LOOP-EXIT      a mid-loop per-VCF crash that records rc=0 plus a SHORTER wall
#                  time when the loop's exit status is its last statement
#                  (`true`).
#   FAIL-CLOSED    a missing bcftools that degrades to `cp` of the raw input with
#                  only a WARN, silently voiding the canonical-contigs guarantee
#                  every reported figure rests on.
#   CELL-ISOLATION one failing cell that aborts the whole run instead of being
#                  recorded and skipped.
#   GROUND-TRUTH   a suite with no Perl ground truth that leaves a dangling
#                  symlink, which a copy of the output tree reports as an error.
#   WARMUP-OUTPUT  a timed run that starts with the discarded warmup's output
#                  still on disk, whose pages the timed run then reclaims while
#                  writing its own.
#   EXTRA-ARGS     --vep-extra-args that reach the SNP/indel call but not the SV
#                  loop's subshell, or a provenance record that names other
#                  arguments than the ones the binary ran with, hashes a file the
#                  run did not execute, or sits beside an earlier run's wrapper;
#                  both loops run the wrapper, an empty value appends nothing, and
#                  a host without python3 or a binary that is not a readable file
#                  is refused before any cell.
#
# Run: bash scripts/concordance/test_run_clone_measurement.sh

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

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# The per-VCF loop shape, extracted verbatim from the vep-rs SV cell.
run_loop() {
    # $1 = engine stub, $2 = input dir, $3 = output dir
    bash -c '
		set -uo pipefail
		bin="$1"; in_dir="$2"; out_dir="$3"
		shopt -s nullglob
		failed=0
		for vcf in "$in_dir"/*.vcf; do
			b=$(basename "$vcf"); b=${b%.vcf}
			[[ "$b" == "manifest" ]] && continue
			"$bin" "$vcf" "$out_dir/$b.txt" 2>"$out_dir/$b.stderr" || {
				failed=$((failed+1))
				echo "ERROR: per-VCF invocation failed: $b" >&2
			}
		done
		[[ "$failed" -eq 0 ]] || {
			echo "ERROR: $failed per-VCF invocation(s) failed" >&2
			exit "$failed"
		}
	' _ "$1" "$2" "$3"
}

echo "LOOP-EXIT: mid-loop failure must produce a non-zero loop exit"

IN="$TMP/in"
OUT="$TMP/out"
mkdir -p "$IN" "$OUT"
for i in 1 2 3 4 5; do echo "##fileformat=VCFv4.2" >"$IN/v$i.vcf"; done

# Stub engine that succeeds on everything.
cat >"$TMP/engine_ok.sh" <<'EOF'
#!/usr/bin/env bash
echo "annotated" >"$2"
EOF
chmod +x "$TMP/engine_ok.sh"

# Stub engine that fails on exactly one input, mid-loop.
cat >"$TMP/engine_flaky.sh" <<'EOF'
#!/usr/bin/env bash
case "$1" in
*v3.vcf) echo "boom" >&2; exit 42 ;;
esac
echo "annotated" >"$2"
EOF
chmod +x "$TMP/engine_flaky.sh"

run_loop "$TMP/engine_ok.sh" "$IN" "$OUT" 2>/dev/null
rc_ok=$?
[[ "$rc_ok" -eq 0 ]] && ok "all-succeed loop exits 0" || bad "all-succeed loop exited $rc_ok, expected 0"

rm -f "$OUT"/*
run_loop "$TMP/engine_flaky.sh" "$IN" "$OUT" 2>"$TMP/flaky.err"
rc_flaky=$?
[[ "$rc_flaky" -ne 0 ]] &&
    ok "mid-loop failure exits non-zero (rc=$rc_flaky)" ||
    bad "mid-loop failure exited 0 -- a partial cell would record as a fast success"

[[ "$rc_flaky" -eq 1 ]] &&
    ok "exit status equals the failure count (1)" ||
    bad "exit status $rc_flaky does not equal the 1 injected failure"

grep -q "per-VCF invocation failed: v3" "$TMP/flaky.err" &&
    ok "the failing input is named on stderr" ||
    bad "stderr does not name the failing input"

# The loop must still process every OTHER input: a fail-open harness exists so
# one bad VCF does not truncate the cell's wall time.
produced=$(find "$OUT" -name '*.txt' -size +0c | wc -l | tr -d ' ')
[[ "$produced" -eq 4 ]] &&
    ok "the other 4 inputs still annotated (loop did not abort)" ||
    bad "only $produced of 4 surviving inputs annotated"

rm -f "$OUT"/*
cat >"$TMP/engine_dead.sh" <<'EOF'
#!/usr/bin/env bash
exit 7
EOF
chmod +x "$TMP/engine_dead.sh"
run_loop "$TMP/engine_dead.sh" "$IN" "$OUT" 2>/dev/null
rc_dead=$?
[[ "$rc_dead" -eq 5 ]] &&
    ok "all-fail loop exit status equals 5 failures" ||
    bad "all-fail loop exited $rc_dead, expected 5"

echo
echo "FAIL-CLOSED: canonicalize must fail closed when bcftools is unavailable"

# The fallback under test: no bcftools means the canonical-contigs guarantee cannot be honored,
# so the cell fails rather than measuring a different (unfiltered,
# genotype-carrying) input than the one the ground truth was built from.
canonicalize_hardened() {
    local in_vcf="$1" out_vcf="$2"
    if ! command -v bcftools_absent_stub >/dev/null 2>&1; then
        echo "ERROR: [run_clone_measurement] bcftools unavailable; cannot produce a canonical-contig input" >&2
        return 1
    fi
    cp "$in_vcf" "$out_vcf"
}

echo "##fileformat=VCFv4.2" >"$TMP/raw.vcf"
canonicalize_hardened "$TMP/raw.vcf" "$TMP/canon.vcf" 2>"$TMP/canon.err"
rc_canon=$?
[[ "$rc_canon" -ne 0 ]] &&
    ok "missing bcftools fails the cell instead of cp-ing raw input" ||
    bad "missing bcftools returned 0 -- canonical-contigs guarantee silently voided"

[[ ! -f "$TMP/canon.vcf" ]] &&
    ok "no unfiltered input is left behind as if it were canonical" ||
    bad "an unfiltered VCF was written to the canonical path"

grep -q "cannot produce a canonical-contig input" "$TMP/canon.err" &&
    ok "the reason is stated on stderr" ||
    bad "stderr does not explain the failure"

echo
echo "CELL-ISOLATION: a failed cell is recorded and SKIPPED, never fatal"

# The two gates above make measure_p1/measure_sv return non-zero. The main suite
# loop runs under `set -euo pipefail`, so calling them bare makes the first bad
# cell kill the script, and a caller that ALSO runs under `set -e` inherits that
# abort, skipping every remaining step of the run. This asserts the guard shape
# in the shipped script rather than a copy of it: the two call sites must
# tolerate a non-zero return, and the run must still reach its completion line.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"

if [[ -f "$HARNESS" ]]; then
    # COMMENT LINES EXCLUDED. Without the exclusion the measure_p1 arm is unfalsifiable:
    # the harness's own explanatory comment reads
    #   `# site is guarded (\`measure_p1 ... || failed_suites+=(...)\`). Exiting`
    # which the pattern matches, so deleting the REAL guard would still report ok.
    for _fn in measure_p1 measure_sv; do
        if grep -e "$_fn .*||" "$HARNESS" | grep -qv '^[[:space:]]*#'; then
            ok "$_fn's call site tolerates a non-zero return"
        else
            bad "$_fn is called bare under set -e -- one bad cell aborts the run"
        fi
    done

    # The guard above is only load-bearing if the function RETURNS on failure. An
    # `exit` inside measure_p1 bypasses the `||` entirely and kills the run, so
    # the fastvep-normalizer fail-closed path must use `return 1`. Scoped to the
    # function body: the `exit 1` calls outside it, in argument parsing and setup,
    # are where aborting is correct. Comments are stripped so the rationale
    # comment beside the guard cannot itself trip this.
    p1_exits=$(sed -n '/^measure_p1() {/,/^}/p' "$HARNESS" |
        grep -vE '^[[:space:]]*#' | grep -c 'exit 1') || p1_exits=0
    [[ "$p1_exits" -eq 0 ]] &&
        ok "measure_p1 never calls exit -- its failures return to the guarded loop" ||
        bad "measure_p1 contains $p1_exits exit call(s) -- a failed cell kills the run"

    # Simulate the guarded loop: cell 2 of 3 fails, and all three must be visited.
    visited=$(bash -c '
		set -euo pipefail
		failed_suites=()
		cell() { [[ "$1" == "s02" ]] && return 1; return 0; }
		for s in s01 s02 s03; do
			echo -n "$s "
			cell "$s" || failed_suites+=("$s")
		done
		echo "| failed=${#failed_suites[@]}"
	')
    [[ "$visited" == "s01 s02 s03 | failed=1" ]] &&
        ok "the loop visits every cell and tallies the failure ($visited)" ||
        bad "guarded loop did not visit all cells: '$visited'"

    # And the unguarded shape must NOT survive -- proving the guard is load-bearing.
    unguarded=$(bash -c '
		set -euo pipefail
		cell() { [[ "$1" == "s02" ]] && return 1; return 0; }
		for s in s01 s02 s03; do echo -n "$s "; cell "$s"; done
		echo "| reached-end"
	' 2>/dev/null)
    [[ "$unguarded" != *"reached-end"* ]] &&
        ok "the unguarded shape dies mid-loop, which is what the guard prevents" ||
        bad "unguarded loop reached the end -- the test cannot detect the defect"
else
    bad "run_clone_measurement.sh not found next to this test"
fi

echo
echo "GROUND-TRUTH: a suite with no Perl ground truth must not leave a DANGLING symlink"

# The fault: the sample-scaling suites s05m/s06m have no Perl ground truth of their own
# (the genotype matrix's cost is measured separately, and genotypes change no
# consequence, so the sites-only truth is the truth for both forms), so a concordance
# branch that symlinks the nonexistent target anyway leaves a dangling link, which a
# copy of the output tree is liable to treat as an error: a complete result set
# reported as FAILED.
#
# Both halves are asserted, because a guard keyed on existence alone would fix this fault
# and introduce its mirror image: a genuinely unstaged s01 reading as a clean skip.

if [[ -f "$HARNESS" ]]; then
    # One guard per engine branch: the vep-rs cell and the fastvep cell.
    n=$(grep -cE 'WALLTIME_ONLY.*!=.*true.*&&.*ground_truth_ready' "$HARNESS") || n=0
    [[ "$n" -eq 2 ]] &&
        ok "both concordance branches are guarded by ground_truth_ready ($n of 2)" ||
        bad "expected 2 ground_truth_ready guards in the shipped harness, found $n"

    # Every ground-truth symlink must sit inside a guarded branch. Counting guards alone
    # would pass if a third symlink site were added outside one, so this walks the file and
    # ties each `ln -sf` to the WALLTIME_ONLY condition that encloses it. The two
    # WALLTIME_ONLY branches in measure_sv are deliberately unguarded and create no
    # symlink: s07/s08 always have ground truth and there is no sample-scaling SV suite.
    unguarded_n=$(awk '
        /if \[\[ "\$WALLTIME_ONLY" != "true" \]\]/ { cond = $0 }
        /ln -sf .*perl\/output\.txt/ { if (cond !~ /ground_truth_ready/) n++ }
        END { print n + 0 }
    ' "$HARNESS")
    [[ "$unguarded_n" -eq 0 ]] &&
        ok "every ground-truth symlink sits inside a ground_truth_ready branch" ||
        bad "$unguarded_n symlink site(s) still stage ground truth unconditionally"

    # Behavioural halves, with the predicate extracted from the shipped script so the test
    # exercises the real function body rather than a restatement of it.
    eval "$(sed -n '/^ground_truth_ready() {/,/^}/p' "$HARNESS")"
    csv_row() { echo "$*" >>"$TMP/rows.txt"; }
    : >"$TMP/rows.txt"
    mkdir -p "$TMP/gt_present" "$TMP/gt_absent" "$TMP/link"
    echo "ground truth" >"$TMP/gt_present/output.txt"

    # (a) a sample-scaling suite: skipped, quietly, with no csv failure row.
    if ground_truth_ready s05m 1KG_multisample GRCh37 "$TMP/gt_absent" >/dev/null; then
        bad "ground_truth_ready accepted s05m, which has no ground truth by design"
    else
        [[ ! -s "$TMP/rows.txt" ]] &&
            ok "s05m is skipped without recording a failure row" ||
            bad "s05m recorded a failure row: $(cat "$TMP/rows.txt")"
    fi

    # (b) a headline suite with the ground truth genuinely missing: loud, with a row.
    if ground_truth_ready s01 ClinVar_GRCh37 GRCh37 "$TMP/gt_absent" 2>/dev/null; then
        bad "ground_truth_ready accepted s01 with no ground truth staged"
    else
        grep -q 'ground_truth_missing' "$TMP/rows.txt" &&
            ok "an unstaged s01 fails the cell and records ground_truth_missing" ||
            bad "an unstaged s01 was skipped silently, which hides a staging failure"
    fi

    # (c) the happy path still passes.
    ground_truth_ready s01 ClinVar_GRCh37 GRCh37 "$TMP/gt_present" >/dev/null &&
        ok "a staged suite still reaches the comparator" ||
        bad "ground_truth_ready rejected a suite whose ground truth exists"

    # (d) the mechanism itself: a symlink to a missing target is a file a copy of the
    #     output tree is liable to skip or error on, which turns a complete run into a
    #     FAILED one.
    ln -sf "$TMP/gt_absent/output.txt" "$TMP/link/output.txt"
    if [[ -L "$TMP/link/output.txt" && ! -e "$TMP/link/output.txt" ]]; then
        ok "a symlink to a missing target dangles, which is the error this guard removes"
    else
        bad "could not reproduce the dangling symlink; the test cannot detect the defect"
    fi
else
    bad "run_clone_measurement.sh not found next to this test"
fi

echo
echo "SV-FORK-PARITY: every engine's SV cell must take the SAME --fork as its SNP/indel cell"

# The fault: a vep-rs SV branch that hardcodes `--fork 1` while the Perl and fastVEP SV
# branches pass `$FORK` (16). On the two SV cells vep-rs is then the only engine of the
# three annotating serially, and every SV ratio is a floor on its advantage rather than a
# measurement of it. The 1 is not needed to match the per-VCF SV ground truth: vep-rs
# output is fork-invariant (FORK-INVARIANCE below) and the ground truth is a fixed
# archived artifact.
#
# Asserted structurally rather than behaviourally: an engine invocation's --fork argument is
# not observable without a real cache and a real 30-second run, and the defect is a literal
# in the source. Any hardcoded numeric --fork inside measure_sv fails this.

if [[ -f "$HARNESS" ]]; then
    sv_body=$(sed -n '/^measure_sv() {/,/^}/p' "$HARNESS" | grep -vE '^[[:space:]]*#')

    hardcoded=$(printf '%s\n' "$sv_body" | grep -cE -- '--fork ("?[0-9]+"?)') || hardcoded=0
    [[ "$hardcoded" -eq 0 ]] &&
        ok "measure_sv hardcodes no numeric --fork" ||
        bad "measure_sv hardcodes a numeric --fork on $hardcoded line(s) -- that engine is handicapped against the others"

    # All three engine branches must forward the harness-level FORK. Counting only the
    # absence of a literal would pass a branch that omitted --fork entirely and silently
    # took the CLI default.
    forwarded=$(printf '%s\n' "$sv_body" | grep -cE -- '--fork "\$(fork|FORK)"|RAYON_NUM_THREADS="\$threads"') || forwarded=0
    [[ "$forwarded" -ge 3 ]] &&
        ok "all three SV engine branches forward the harness --fork ($forwarded sites)" ||
        bad "only $forwarded of 3 SV engine branches forward --fork; the others take a default"

    # And each branch must actually be PASSED $FORK from the enclosing scope, not read a
    # variable that is never set inside the `bash -c` subshell (which would expand empty).
    # Line continuations are collapsed first: the Perl branch's argument list wraps across
    # two lines, so a per-line match reports 2 of 3 and fails a correct harness.
    passed_in=$(printf '%s\n' "$sv_body" | sed -e :a -e '/\\$/N; s/\\\n[[:space:]]*/ /; ta' |
        grep -cE "' _ .*\"\\\$FORK\"") || passed_in=0
    [[ "$passed_in" -eq 3 ]] &&
        ok "all three SV subshells receive \$FORK in their argument list ($passed_in of 3)" ||
        bad "$passed_in of 3 SV subshells receive \$FORK -- an unset fork expands empty"
fi

echo
echo "FORK-INVARIANCE: vep-rs output must be independent of --fork"

# The premise SV-FORK-PARITY rests on, asserted at the source rather than in prose: the
# parallel and serial annotation paths must run the SAME per-variant function over the SAME
# buffer, and output must be emitted in input order. A parallel path with its own body
# would make --fork more than a performance knob and every fork-crossing comparison in
# the harness unsound.
#
# The behavioural counterpart is the `test_annotate_batch_fork_1_matches_fork_4` unit test in
# crates/vep-cli/src/runner.rs, which annotates one batch at both fork settings and asserts
# positional equality of the emitted consequences. This check is the cheap always-on half: it
# needs no cache and runs in a bash-only image.

RUNNER="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/crates/vep-cli/src/runner.rs"
if [[ -f "$RUNNER" ]]; then
    # Both arms of the use_parallel branch call annotate_one, and neither inlines a
    # second implementation.
    par_calls=$(grep -cE 'par_iter_mut\(\)' "$RUNNER") || par_calls=0
    ser_calls=$(grep -cE 'iter_mut\(\)' "$RUNNER") || ser_calls=0
    [[ "$par_calls" -ge 1 && "$ser_calls" -gt "$par_calls" ]] &&
        ok "runner.rs keeps both a parallel and a serial arm over the same batch" ||
        bad "runner.rs has $par_calls parallel and $((ser_calls - par_calls)) serial batch arms; the fork paths may have diverged"

    annotate_one_defs=$(grep -cE '^[[:space:]]*let annotate_one =( |$)' "$RUNNER") || annotate_one_defs=0
    [[ "$annotate_one_defs" -eq 1 ]] &&
        ok "one annotate_one closure serves both fork paths" ||
        bad "$annotate_one_defs annotate_one definitions found; expected exactly 1 shared by both paths"

    # The behavioural test the comment above points at must actually exist. A dangling
    # pointer here is worse than none: it tells a maintainer the expensive half is already
    # covered, so nobody writes it.
    fork_test=$(grep -cE 'fn test_annotate_batch_fork_1_matches_fork_4\(\)' "$RUNNER") || fork_test=0
    [[ "$fork_test" -eq 1 ]] &&
        ok "the behavioural fork-invariance test exists in runner.rs" ||
        bad "test_annotate_batch_fork_1_matches_fork_4 not found in runner.rs; this gate's behavioural counterpart is missing"
else
    bad "crates/vep-cli/src/runner.rs not found; cannot assert the fork-invariance premise"
fi

echo "CACHE-COMPLETE: a chromosome-subset Perl cache must fail the cell before it is timed"

# This gate covers the case assert_perl_annotated provably cannot. Over the 16 GRCh37 SV VCFs,
# whose records are chr21 by design, a chr21-only cache clears that gate's rows>=records
# floor, while `09_breakends` silently loses every off-chromosome mate.
# So the assertion has to be made against the cache, not the output volume.
#
# The real function is extracted from the harness and evaluated here, rather than reimplemented,
# so this exercises the shipped code and cannot drift from it.
HARNESS="$(dirname "$0")/run_clone_measurement.sh"
if [[ -f "$HARNESS" ]]; then
    guard_src=$(awk '/^assert_perl_cache_complete\(\) \{$/,/^}$/' "$HARNESS")
    if [[ -n "$guard_src" ]]; then
        ok "assert_perl_cache_complete is defined in the harness"
        eval "$guard_src"

        # A chr21-only slice.
        PERL_CACHE_DIR="$TMP/cache_chr21only"
        mkdir -p "$PERL_CACHE_DIR/homo_sapiens/116_GRCh37/21"
        : >"$PERL_CACHE_DIR/homo_sapiens/116_GRCh37/info.txt"
        if assert_perl_cache_complete GRCh37 >/dev/null 2>&1; then
            bad "a chr21-only GRCh37 cache PASSED the completeness guard"
        else
            ok "a chr21-only GRCh37 cache fails the guard"
        fi

        # All 25 primary contigs, no scaffolds. Ensembl's full tree carries 699 directories, so a
        # guard demanding all of them would fail this legitimately pruned cache.
        PERL_CACHE_DIR="$TMP/cache_primary"
        for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
            mkdir -p "$PERL_CACHE_DIR/homo_sapiens/116_GRCh37/$c"
        done
        if assert_perl_cache_complete GRCh37 >/dev/null 2>&1; then
            ok "a cache with all 25 primary contigs and no scaffolds passes"
        else
            bad "a 25-primary-contig cache FAILED the guard; it must not require scaffolds"
        fi

        # One contig short is the interesting near-miss: 24 of 25 is what an incomplete
        # cache extraction looks like, and it must not pass.
        rmdir "$PERL_CACHE_DIR/homo_sapiens/116_GRCh37/MT"
        if assert_perl_cache_complete GRCh37 >/dev/null 2>&1; then
            bad "a cache missing only MT PASSED the guard"
        else
            ok "a cache missing a single primary contig fails the guard"
        fi

        # A missing slice must fail loudly rather than pass by absence: a check that measured
        # nothing reports clean.
        PERL_CACHE_DIR="$TMP/cache_absent"
        if assert_perl_cache_complete GRCh37 >/dev/null 2>&1; then
            bad "an ABSENT cache slice PASSED the guard"
        else
            ok "an absent cache slice fails the guard rather than passing by omission"
        fi
    else
        bad "assert_perl_cache_complete not found in run_clone_measurement.sh"
    fi

    # Both Perl branches must call it. A guard defined and never invoked is the failure mode
    # this assertion exists for.
    guard_calls=$(grep -cE '^[[:space:]]*if ! assert_perl_cache_complete ' "$HARNESS") || guard_calls=0
    [[ "$guard_calls" -eq 2 ]] &&
        ok "both Perl branches (SNP/indel and SV) call the cache guard (2 of 2)" ||
        bad "$guard_calls Perl branches call assert_perl_cache_complete; expected 2"

    # It has to run before the timing loop, or a failing cell is timed anyway -- and the
    # output-volume assert has to run after it and before the scratch delete, or the
    # evidence is gone before anything looks at it. Assert the ORDER, not the presence.
    #
    # WHY ORDER AND NOT TEXT. A scanner that tracks `case` arms by comparing each stripped
    # line to `perl)` and then matches `*assert_perl_cache_complete*` is satisfied by a
    # COMMENT, and gives three false verdicts:
    #   * false PASS -- delete both real guard calls, leave a comment naming the guard in
    #     each perl) arm, put syntactically-matching calls in the vep-rs) arms (never
    #     entered when ENGINE=perl).
    #   * false PASS -- move the guard AFTER the timing loop, which is the exact defect
    #     the comment above says this assertion prevents.
    #   * false FAIL -- append a trailing comment to a `perl)` label, or reindent it. Both
    #     are no-ops `bash -n` accepts with both guards still executing.
    # One that `break`s on the first hit in any arm never verifies the second one.
    #
    # NO awk: a literal \t is read differently by BSD and busybox awk. `grep -n` plus
    # `cut` are POSIX and present in busybox, and neither needs a character class.
    _guard_lines=$(grep -n -e '^[[:space:]]*if ! assert_perl_cache_complete ' "$HARNESS" | cut -d: -f1)
    _annot_lines=$(grep -n -e '^[[:space:]]*if ! assert_perl_annotated ' "$HARNESS" | cut -d: -f1)
    _loop_lines=$(grep -n -e 'for run_kind in warmup timed' "$HARNESS" | cut -d: -f1)
    _rm_lines=$(grep -n -e 'rm -rf "\$suite_dir/perl_scratch"' "$HARNESS" | cut -d: -f1)

    _n_annot=$(printf '%s\n' $_annot_lines | grep -c . ) || _n_annot=0
    [[ "$_n_annot" -eq 2 ]] &&
        ok "both Perl branches call the output-volume gate (2 of 2)" ||
        bad "$_n_annot Perl branches call assert_perl_annotated; expected 2 (SNP/indel and SV)"

    if grep -q '^assert_perl_annotated()' "$HARNESS"; then
        ok "assert_perl_annotated is defined in the harness"
    else
        bad "assert_perl_annotated not found in run_clone_measurement.sh"
    fi

    # Pair each guard with the next output-volume assert after it; that pair brackets one
    # Perl branch's timed work. A timing loop must sit BETWEEN them, and a scratch delete
    # must NOT sit before the assert.
    _ordered=0
    _branches=0
    for _g in $_guard_lines; do
        _a=""
        for _cand in $_annot_lines; do
            if [[ "$_cand" -gt "$_g" ]]; then _a="$_cand"; break; fi
        done
        [[ -n "$_a" ]] || continue
        _branches=$((_branches + 1))
        _loop_between=0
        for _l in $_loop_lines; do
            if [[ "$_l" -gt "$_g" && "$_l" -lt "$_a" ]]; then _loop_between=1; fi
        done
        _rm_before_assert=0
        for _r in $_rm_lines; do
            if [[ "$_r" -gt "$_g" && "$_r" -lt "$_a" ]]; then _rm_before_assert=1; fi
        done
        if [[ "$_loop_between" -eq 1 && "$_rm_before_assert" -eq 0 ]]; then
            _ordered=$((_ordered + 1))
        fi
    done
    if [[ "$_branches" -eq 2 && "$_ordered" -eq 2 ]]; then
        ok "on both Perl branches the cache guard precedes the timed loop, which precedes the output-volume gate, which precedes the scratch delete"
    else
        bad "only $_ordered of $_branches Perl branches order guard -> timed loop -> output-volume gate -> scratch delete"
    fi

    # THE FLAG. `--full_cache_dir` is not a VEP flag: VEP ignores it, finds no annotation
    # source, emits one warning per record and exits 0, so a Perl wall time measured with
    # it times a VEP that annotated nothing. A live occurrence at any site passes
    # `bash -n`, so this assertion is the only check on it.
    _bad_flag=$(grep -c -e '--full_cache_dir' "$HARNESS") || _bad_flag=0
    _bad_flag_live=$(grep -e '--full_cache_dir' "$HARNESS" | grep -vc '^[[:space:]]*#') || _bad_flag_live=0
    [[ "$_bad_flag_live" -eq 0 ]] &&
        ok "no live --full_cache_dir in the harness ($_bad_flag comment reference(s))" ||
        bad "$_bad_flag_live non-comment line(s) pass --full_cache_dir, which VEP ignores; the cell would time a VEP that annotates nothing"

    _dir_cache_live=$(grep -e '--dir_cache' "$HARNESS" | grep -vc '^[[:space:]]*#') || _dir_cache_live=0
    [[ "$_dir_cache_live" -ge 2 ]] &&
        ok "the Perl invocations pass --dir_cache ($_dir_cache_live live site(s))" ||
        bad "only $_dir_cache_live live --dir_cache site(s); expected at least 2 (SNP/indel and SV Perl timing)"
else
    bad "run_clone_measurement.sh not found; cannot assert the cache guard"
fi

echo "SV-INPUT-GATE: an SV cell whose input glob matches nothing must fail, not time an empty loop"

# The per-VCF loop runs under `nullglob`. A data tree whose SV root does not match the
# harness's SV_GT_DIR times an empty loop: wall 0 s, rc 0, no output, no comparison, and a
# clean exit. Both ends of the loop are gated; this exercises the shipped gates.
if [[ -f "$HARNESS" ]]; then
    for fn in sv_input_count assert_sv_inputs_present assert_sv_outputs_present; do
        _src=$(awk -v f="$fn" '$0 ~ "^"f"\\(\\) \\{$",/^}$/' "$HARNESS")
        if [[ -n "$_src" ]]; then
            eval "$_src"
        else
            bad "$fn not found in run_clone_measurement.sh"
        fi
    done
    if declare -f assert_sv_inputs_present >/dev/null && declare -f assert_sv_outputs_present >/dev/null; then
        ok "the SV input and output gates are defined in the harness"
        SV_GT_DIR="$TMP/sv_gt_probe"
        # absent directory
        if assert_sv_inputs_present "$TMP/sv_gt_probe/canonical_inputs/grch37" s07 >/dev/null 2>&1; then
            bad "an ABSENT SV input directory passed the input gate"
        else
            ok "an absent SV input directory fails the input gate"
        fi
        # present but holding only the manifest
        mkdir -p "$TMP/sv_only_manifest"
        : >"$TMP/sv_only_manifest/manifest.json"
        if assert_sv_inputs_present "$TMP/sv_only_manifest" s07 >/dev/null 2>&1; then
            bad "a directory holding only manifest.json passed the input gate"
        else
            ok "a directory holding only manifest.json fails the input gate"
        fi
        # one real input, plain and gzipped spellings both count
        mkdir -p "$TMP/sv_two"
        : >"$TMP/sv_two/01_snv.vcf"
        : >"$TMP/sv_two/02_mnp.vcf.gz"
        : >"$TMP/sv_two/manifest.json"
        if assert_sv_inputs_present "$TMP/sv_two" s07 >/dev/null 2>&1; then
            ok "a directory with input VCFs passes the input gate"
        else
            bad "a directory with two input VCFs FAILED the input gate"
        fi
        [[ "$(sv_input_count "$TMP/sv_two")" == "2" ]] &&
            ok "sv_input_count counts .vcf and .vcf.gz and skips manifest (2)" ||
            bad "sv_input_count returned $(sv_input_count "$TMP/sv_two"), expected 2"
        # outputs: nothing written -> fail; header-only files -> fail; one data row each -> pass
        mkdir -p "$TMP/sv_out"
        if assert_sv_outputs_present "$TMP/sv_out" 2 s07 >/dev/null 2>&1; then
            bad "an EMPTY output directory passed the output gate"
        else
            ok "an empty output directory fails the output gate"
        fi
        printf '## header only\n' >"$TMP/sv_out/01_snv.txt"
        printf '## header only\n' >"$TMP/sv_out/02_mnp.txt"
        if assert_sv_outputs_present "$TMP/sv_out" 2 s07 >/dev/null 2>&1; then
            bad "header-only outputs passed the output gate"
        else
            ok "header-only outputs fail the output gate"
        fi
        printf '## h\nrow1\n' >"$TMP/sv_out/01_snv.txt"
        if assert_sv_outputs_present "$TMP/sv_out" 2 s07 >/dev/null 2>&1; then
            bad "one non-empty output of two required passed the output gate"
        else
            ok "fewer non-empty outputs than inputs fails the output gate"
        fi
        printf '## h\nrow1\n' >"$TMP/sv_out/02_mnp.txt"
        if assert_sv_outputs_present "$TMP/sv_out" 2 s07 >/dev/null 2>&1; then
            ok "one non-empty output per input passes the output gate"
        else
            bad "complete outputs FAILED the output gate"
        fi
    fi
    # The gates have to be CALLED: the input gate before the pointer/engine dispatch, the
    # output gate inside both timing loops (vep-rs and fastVEP), each followed by a failed
    # csv_row and `return 1`.
    _in_calls=$(grep -c -e '^[[:space:]]*if ! assert_sv_inputs_present ' "$HARNESS") || _in_calls=0
    [[ "$_in_calls" -eq 1 ]] &&
        ok "measure_sv calls the input gate once, before dispatch" ||
        bad "$_in_calls call(s) to assert_sv_inputs_present; expected 1"
    _out_calls=$(grep -c -e '^[[:space:]]*if ! assert_sv_outputs_present ' "$HARNESS") || _out_calls=0
    [[ "$_out_calls" -eq 2 ]] &&
        ok "both engine timing loops (vep-rs, fastVEP) call the output gate (2 of 2)" ||
        bad "$_out_calls call(s) to assert_sv_outputs_present; expected 2"
    _in_line=$(grep -n -e '^[[:space:]]*if ! assert_sv_inputs_present ' "$HARNESS" | head -1 | cut -d: -f1)
    # The SNP/indel loop uses the same time_run shape earlier in the file, so anchor on the SV
    # function's own definition and take the first timing call after it.
    _sv_def=$(grep -n -e '^measure_sv() {' "$HARNESS" | head -1 | cut -d: -f1)
    _first_time_run=$(grep -n -e 'res=\$(time_run "\${suite_id}_\${run_kind}"' "$HARNESS" | cut -d: -f1 | awk -v d="${_sv_def:-0}" '$1 > d' | head -1)
    if [[ -n "$_in_line" && -n "$_first_time_run" && "$_in_line" -lt "$_first_time_run" ]]; then
        ok "the input gate (line $_in_line) precedes the first SV timing call (line $_first_time_run)"
    else
        bad "the input gate does not precede the SV timing loop (gate line ${_in_line:-none}, time_run line ${_first_time_run:-none})"
    fi
    # The SV root must be overridable by env so a caller can point at its own tree.
    if grep -q -e '^SV_GT_DIR="\${VEP_SV_GT_DIR:-\$DATA_DIR/ground_truth/perl/sv_per_vcf}"$' "$HARNESS"; then
        ok "SV_GT_DIR honours the VEP_SV_GT_DIR override"
    else
        bad "SV_GT_DIR does not honour a VEP_SV_GT_DIR override"
    fi
fi

echo "COMPARE-FAIL-CLOSED: a comparator that runs and fails must fail the cell, never warn and continue"

# A comparison ending in `|| echo "[WARN] ... continuing"` lets a run whose comparator
# crashed (a missing mask cache, a traceback) record its wall times, write no report, and
# exit 0. The comparator calls must fail the cell by name.
if [[ -f "$HARNESS" ]]; then
    _warn_continue=$(grep -c 'compare failed.*continuing' "$HARNESS") || _warn_continue=0
    [[ "$_warn_continue" -eq 0 ]] &&
        ok "no comparator call warns and continues" ||
        bad "$_warn_continue comparator call(s) still warn and continue on failure"
    _fail_rows=$(grep -c '"compare_failed" "" "1"' "$HARNESS") || _fail_rows=0
    [[ "$_fail_rows" -ge 4 ]] &&
        ok "the four comparator calls record compare_failed and return 1 ($_fail_rows)" ||
        bad "only $_fail_rows comparator call(s) record compare_failed; expected 4"
    grep -q '"mask_cache_missing" "" "1"' "$HARNESS" &&
        ok "the fastVEP SV branch asserts the vep-rs mask cache before comparing" ||
        bad "the fastVEP SV branch does not assert the vep-rs mask cache"
fi

# ---------------------------------------------------------------------------
# Field-level agreement rides beside the F1 report on every vep-rs SNP/indel cell: the
# field comparator runs inside the same ground_truth_ready branch, writes
# report/fields.json, and on failure records a CSV row without failing the cell
# (its report is not F1-bearing).
if [[ -f "$HARNESS" ]]; then
    _fields_in_guard=$(awk '
        /if \[\[ "\$WALLTIME_ONLY" != "true" \]\]/ { cond = $0 }
        /compare_vep_fields\.py/ { if (cond ~ /ground_truth_ready/) n++ }
        END { print n + 0 }
    ' "$HARNESS")
    [[ "$_fields_in_guard" -eq 1 ]] &&
        ok "compare_vep_fields.py runs once, inside the guarded vep-rs SNP/indel branch" ||
        bad "expected compare_vep_fields.py once inside a ground_truth_ready branch, found $_fields_in_guard"
    grep -q -- '--out "$suite_dir/report/fields.json"' "$HARNESS" &&
        ok "the field comparator writes report/fields.json beside the F1 report" ||
        bad "the field comparator does not write report/fields.json"
    grep -q '"fields_compare_failed" "" "1"' "$HARNESS" &&
        ok "a field-comparison failure records fields_compare_failed and lets the cell complete" ||
        bad "a field-comparison failure is not recorded as fields_compare_failed"
fi

# ---------------------------------------------------------------------------
# --skip-fields drops ONLY the field comparator: the flag is parsed, it defaults
# off (a default run gets the field report by omission), the comparator
# call sits inside its guard, and the F1 compare stays outside it.
if [[ -f "$HARNESS" ]]; then
    grep -q '^SKIP_FIELDS="false"$' "$HARNESS" &&
        ok "--skip-fields defaults off (SKIP_FIELDS=false)" ||
        bad "SKIP_FIELDS does not default to false"
    grep -A 2 -e '^[[:space:]]*--skip-fields)$' "$HARNESS" | grep -q 'SKIP_FIELDS="true"' &&
        ok "--skip-fields is parsed and sets SKIP_FIELDS=true" ||
        bad "--skip-fields is not parsed into SKIP_FIELDS=true"
    _fields_skip_guarded=$(awk '
        /if \[\[ "\$SKIP_FIELDS" != "true" \]\]/ { g = 1 }
        /compare_vep_fields\.py/ { if (g) n++; g = 0 }
        END { print n + 0 }
    ' "$HARNESS")
    [[ "$_fields_skip_guarded" -eq 1 ]] &&
        ok "compare_vep_fields.py runs only while SKIP_FIELDS is not true" ||
        bad "expected compare_vep_fields.py once under the SKIP_FIELDS guard, found $_fields_skip_guarded"
    _skip_line=$(grep -n -e 'if \[\[ "\$SKIP_FIELDS" != "true" \]\]' "$HARNESS" | head -1 | cut -d: -f1)
    _f1_line=$(grep -n -e 'compare_vep_outputs\.py' "$HARNESS" | grep -v '^[0-9]*:[[:space:]]*#' | head -1 | cut -d: -f1)
    [[ -n "$_skip_line" && -n "$_f1_line" && "$_f1_line" -lt "$_skip_line" ]] &&
        ok "the F1 compare is outside the SKIP_FIELDS guard (line $_f1_line before $_skip_line)" ||
        bad "the F1 compare must not depend on SKIP_FIELDS (F1 at ${_f1_line:-?}, guard at ${_skip_line:-?})"
fi

# VEP_FIELDS_COMPARATOR swaps the field comparator's command and nothing else: the
# SKIP_FIELDS block is extracted and run with a python3 stub and a comparator stub
# that record their argv. An executable file runs in the Python's place with the
# same flags, a bare name resolving through PATH; a set value that is not an
# executable file (a directory included) warns and runs the Python; unset runs the
# Python.
if [[ -f "$HARNESS" ]]; then
    _fields_block=$(awk '
        !p && /^[[:space:]]*if \[\[ "\$SKIP_FIELDS" != "true" \]\]; then$/ { match($0, /^[[:space:]]*/); ind = substr($0, 1, RLENGTH); p = 1 }
        p { print }
        p && $0 == ind "fi" { exit }
    ' "$HARNESS")
    eval "_fields_block_fn() {
$_fields_block
}"
    _fdir="$TMP/fields_dispatch"
    mkdir -p "$_fdir/bin" "$_fdir/suite/report"
    for _stub in python3 comparator; do
        cat >"$_fdir/bin/$_stub" <<EOF
#!/bin/bash
printf '%s\n' "\$0" "\$@" >"$_fdir/argv"
EOF
        chmod +x "$_fdir/bin/$_stub"
    done
    : >"$_fdir/not_executable"
    # $1 is the VEP_FIELDS_COMPARATOR value, "" for unset; clears the stubs' argv record,
    # then prints the block's stdout and stderr.
    _run_fields_block() {
        rm -f "$_fdir/argv"
        (
            PATH="$_fdir/bin:$PATH"
            REPO_DIR="$_fdir/repo" suite_dir="$_fdir/suite" suite_id=s06 suite_name=stub assembly=GRCh38 SKIP_FIELDS=false
            csv_row() { echo "csv_row $*"; }
            if [[ -n "$1" ]]; then export VEP_FIELDS_COMPARATOR="$1"; else unset VEP_FIELDS_COMPARATOR; fi
            _fields_block_fn
        ) 2>&1
    }
    # The stub the last run invoked: the first line of the argv record, "" when none ran.
    _fields_invoked() { head -n 1 "$_fdir/argv" 2>/dev/null; }
    _out=$(_run_fields_block "")
    [[ "$(_fields_invoked)" == "$_fdir/bin/python3" && "$_out" != *WARNING* ]] &&
        ok "VEP_FIELDS_COMPARATOR unset: the Python comparator runs" ||
        bad "VEP_FIELDS_COMPARATOR unset: expected the python3 stub without a warning, got argv '$(_fields_invoked)' output '$_out'"
    # The Python's flags start after its script path, the second argv line.
    _py_flags=$(tail -n +3 "$_fdir/argv" 2>/dev/null | paste -sd' ' -)
    _out=$(_run_fields_block "$_fdir/bin/comparator")
    [[ "$(_fields_invoked)" == "$_fdir/bin/comparator" && "$_out" != *WARNING* ]] &&
        ok "VEP_FIELDS_COMPARATOR executable: the binary runs in the Python's place" ||
        bad "VEP_FIELDS_COMPARATOR executable: expected the comparator stub without a warning, got argv '$(_fields_invoked)' output '$_out'"
    # The binary's flags start after its own path, the first argv line.
    _bin_flags=$(tail -n +2 "$_fdir/argv" 2>/dev/null | paste -sd' ' -)
    [[ -n "$_py_flags" && "$_py_flags" == "$_bin_flags" && "$_py_flags" == *"--out $_fdir/suite/report/fields.json --tmp-dir $_fdir/suite"* ]] &&
        ok "VEP_FIELDS_COMPARATOR: the binary receives the Python's flags (--perl --rust --out --tmp-dir)" ||
        bad "VEP_FIELDS_COMPARATOR: flag lists differ; python: $_py_flags binary: $_bin_flags"
    for _bad_value in "$_fdir/not_executable" "$_fdir"; do
        _out=$(_run_fields_block "$_bad_value")
        [[ "$(_fields_invoked)" == "$_fdir/bin/python3" && "$_out" == *"WARNING: [run_clone_measurement] VEP_FIELDS_COMPARATOR=$_bad_value is not executable"* ]] &&
            ok "VEP_FIELDS_COMPARATOR=${_bad_value##*/} (not an executable file): warns and runs the Python comparator" ||
            bad "VEP_FIELDS_COMPARATOR=$_bad_value: expected a warning and the python3 stub, got argv '$(_fields_invoked)' output '$_out'"
    done
    # A bare name resolves through PATH, as the exec resolves it; the test's cwd holds no
    # such file, so a guard that tests the name against the cwd warns and runs the Python.
    _out=$(_run_fields_block comparator)
    [[ "$(_fields_invoked)" == "$_fdir/bin/comparator" && "$_out" != *WARNING* ]] &&
        ok "VEP_FIELDS_COMPARATOR=comparator (a bare name on PATH): the binary runs, resolved through PATH" ||
        bad "VEP_FIELDS_COMPARATOR=comparator: expected the comparator stub through PATH without a warning, got argv '$(_fields_invoked)' output '$_out'"
    # bash's PATH search returns a non-executable match when no executable one exists.
    : >"$_fdir/bin/inert"
    _out=$(_run_fields_block inert)
    [[ "$(_fields_invoked)" == "$_fdir/bin/python3" && "$_out" == *"WARNING: [run_clone_measurement] VEP_FIELDS_COMPARATOR=inert is not executable"* ]] &&
        ok "VEP_FIELDS_COMPARATOR=inert (a bare name whose only PATH match is not executable): warns and runs the Python comparator" ||
        bad "VEP_FIELDS_COMPARATOR=inert: expected a warning and the python3 stub, got argv '$(_fields_invoked)' output '$_out'"
fi

echo "CONTAINER-PROBE: the container start-up probe times bare image starts, Perl only, after the suites"

# Every timed Perl invocation pays one container start; --container-probe N measures that
# cost on its own (N bare `docker run --rm <image> true` starts under /usr/bin/time -v, rows
# in <output-dir>/container_start.csv). The flag is Perl-only, defaults to 0, runs after the
# suite loop, and a non-zero probe exit is recorded rather than hidden.
if [[ -f "$HARNESS" ]]; then
    grep -q '^CONTAINER_PROBE=0$' "$HARNESS" &&
        ok "--container-probe defaults to 0 (no probe)" ||
        bad "CONTAINER_PROBE does not default to 0"
    grep -A 2 -e '^[[:space:]]*--container-probe)$' "$HARNESS" | grep -q 'CONTAINER_PROBE="$2"' &&
        ok "--container-probe is parsed into CONTAINER_PROBE" ||
        bad "--container-probe is not parsed"
    grep -q -- '--container-probe applies to --engine perl only' "$HARNESS" &&
        ok "a non-Perl engine with a probe count is rejected" ||
        bad "the Perl-only guard on --container-probe is missing"
    grep -q 'docker run --rm --user "$(id -u):$(id -g)" "$PERL_IMAGE" true' "$HARNESS" &&
        ok "the probe starts the image with the timed cells' --user flag and runs its true" ||
        bad "the probe does not start the bare image with the timed cells' --user flag"
    grep -q '^[[:space:]]*echo "engine,run_id,host_id,replicate,image,probe_index,wall_seconds,exit_code" >"$csv"$' "$HARNESS" &&
        ok "container_start.csv carries the wall_times.csv provenance columns plus image and probe index" ||
        bad "container_start.csv header differs from the documented columns"
    _probe_def=$(grep -n -e '^container_probe() {' "$HARNESS" | head -1 | cut -d: -f1)
    _loop_end=$(grep -n -e '^done$' "$HARNESS" | awk -F: -v d="${_probe_def:-0}" '$1 < d' | tail -1 | cut -d: -f1)
    _last_measure=$(grep -n -e 'measure_sv "$suite_id" "$suite_name" "$assembly"' "$HARNESS" | tail -1 | cut -d: -f1)
    if [[ -n "$_probe_def" && -n "$_last_measure" && "$_last_measure" -lt "$_probe_def" ]]; then
        ok "the probe is defined and run after the suite loop (loop measure at line $_last_measure, probe at $_probe_def)"
    else
        bad "the probe does not follow the suite loop (measure at ${_last_measure:-none}, probe at ${_probe_def:-none})"
    fi
    grep -q 'container probe $i exited $rc' "$HARNESS" &&
        ok "a failing probe start is reported by index" ||
        bad "a failing probe start is not reported"
    # Functional: the probe loop shape against a stub docker, where GNU time is available.
    if /usr/bin/time -v true >/dev/null 2>&1; then
        _pd="$TMP/probe"; mkdir -p "$_pd/bin"
        cat >"$_pd/bin/docker" <<'STUB'
#!/usr/bin/env bash
# stub: the third start fails
[[ -f "$PROBE_STATE" ]] && n=$(cat "$PROBE_STATE") || n=0
n=$((n + 1)); echo "$n" >"$PROBE_STATE"
[[ "$n" -eq 3 ]] && exit 7
exit 0
STUB
        chmod +x "$_pd/bin/docker"
        _probe_out=$(PROBE_STATE="$_pd/state" PATH="$_pd/bin:$PATH" bash -c '
            set -uo pipefail
            ENGINE=perl RUN_ID=run-test HOST_ID=host-test REPLICATE=1 PERL_IMAGE=stub:image OUTPUT_DIR="$1"
            source <(sed -n "/^time_run() {/,/^}/p" "$2")
            source <(sed -n "/^container_probe() {/,/^}/p" "$2")
            container_probe 4 2>&1
        ' _ "$_pd" "$HARNESS")
        _rows=$(grep -c '^perl,run-test,host-test,1,stub:image,' "$_pd/container_start.csv") || _rows=0
        [[ "$_rows" -eq 4 ]] &&
            ok "four probe starts write four rows" ||
            bad "expected 4 probe rows, found $_rows"
        grep -q '^perl,run-test,host-test,1,stub:image,3,[^,]*,7$' "$_pd/container_start.csv" &&
            ok "the failing third start records exit_code 7 in its row" ||
            bad "the failing start's exit code is not in its row"
        echo "$_probe_out" | grep -q 'container probe 3 exited 7' &&
            ok "the failing start is reported on stderr" ||
            bad "the failing start was not reported"
        [[ -f "$_pd/container_probe/probe-4.time.txt" ]] &&
            ok "each start keeps its /usr/bin/time -v block" ||
            bad "the per-start timing block is missing"
    else
        echo "  skip functional probe run: /usr/bin/time -v is not GNU time here"
    fi
fi

echo
echo "WARMUP-OUTPUT: the timed run must start with the warmup's output gone"

# The function under test is the harness's own, extracted by name so the test
# proves the shipped definition rather than a copy of it.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"
eval "$(sed -n '/^discard_warmup_output() {/,/^}/p' "$HARNESS")"
declare -F discard_warmup_output >/dev/null &&
    ok "discard_warmup_output is defined in the harness" ||
    bad "discard_warmup_output is not defined in the harness"

# Stub engine: writes its output like a run, and on the timed run records
# whether the warmup's output and stderr were still present when it started.
SUITE="$TMP/suite"
mkdir -p "$SUITE/rust" "$SUITE/perl_scratch" "$SUITE/perf"
cat >"$TMP/engine_warm.sh" <<'EOF'
#!/usr/bin/env bash
# $1 = run kind, $2 = suite dir
if [[ "$1" == timed ]]; then
    n=$(find "$2/rust" "$2/perl_scratch" -maxdepth 1 -type f \( -name '*.txt' -o -name '*.stderr' \) | wc -l | tr -d ' ')
    echo "$n" >"$2/perf/leftover_at_timed_start"
fi
head -c 4096 /dev/zero >"$2/rust/output.txt"
echo "log" >"$2/rust/output.stderr"
echo "scratch" >"$2/perl_scratch/output.txt"
EOF
chmod +x "$TMP/engine_warm.sh"

# The harness's warmup-then-timed loop shape, with the discard step.
for run_kind in warmup timed; do
    "$TMP/engine_warm.sh" "$run_kind" "$SUITE"
    [[ "$run_kind" == warmup ]] && discard_warmup_output "$SUITE"
done
leftover=$(cat "$SUITE/perf/leftover_at_timed_start")
[[ "$leftover" -eq 0 ]] &&
    ok "the timed run started with no warmup output on disk" ||
    bad "the timed run found $leftover warmup file(s) still on disk"
[[ -s "$SUITE/rust/output.txt" ]] &&
    ok "the timed run's own output survives the step" ||
    bad "the timed run's output is missing after the step"

# The same loop without the step is the shape the test exists to catch.
rm -f "$SUITE"/rust/* "$SUITE"/perl_scratch/* "$SUITE/perf/leftover_at_timed_start"
for run_kind in warmup timed; do
    "$TMP/engine_warm.sh" "$run_kind" "$SUITE"
done
leftover=$(cat "$SUITE/perf/leftover_at_timed_start")
[[ "$leftover" -eq 3 ]] &&
    ok "without the step the timed run starts over the warmup's 3 files (the fault the step removes)" ||
    bad "without the step the timed run found $leftover file(s), expected 3"

echo
echo "PERL-SV-CACHE-VERSION: the Perl SV cell's child shell must receive the cache version it passes to vep"

# The fault: the per-VCF Perl loop runs in a single-quoted `bash -c` child under `set -u`.
# A harness variable the parent never exported is unbound there, so a child that reads
# `$PERL_CACHE_VERSION` aborts on its first record with "unbound variable", every Perl SV
# cell records a non-zero exit, and the campaign has no Perl SV wall time. The value must
# travel as a positional argument like the other nine.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"
if [[ -f "$HARNESS" ]]; then
    # The child body: the single-quoted `bash -c` block of the Perl branch of measure_sv,
    # extracted from the shipped script so the test runs the code that ships.
    perl_child=$(sed -n '/^measure_sv() {/,/^}/p' "$HARNESS" |
        awk '/^[[:space:]]*perl\)$/ { p = 1 } p' |
        sed -n "/bash -c '/,/^[[:space:]]*' _ \"\$PERL_IMAGE\"/p" |
        sed '1d;$d')
    if [[ -n "$perl_child" ]]; then
        ok "the Perl SV child body was extracted from the harness ($(printf '%s\n' "$perl_child" | wc -l | tr -d ' ') lines)"
    else
        bad "the Perl SV child body could not be extracted from the harness"
    fi

    # Structural: the child names no harness-level variable, and the argument line hands
    # it the cache version.
    _unbound=$(printf '%s\n' "$perl_child" | grep -c 'PERL_CACHE_VERSION') || _unbound=0
    [[ "$_unbound" -eq 0 ]] &&
        ok "the child reads no \$PERL_CACHE_VERSION of its own" ||
        bad "the child reads \$PERL_CACHE_VERSION on $_unbound line(s); unexported, it is unbound under set -u"
    _perl_args=$(sed -n '/^measure_sv() {/,/^}/p' "$HARNESS" | sed -e :a -e '/\\$/N; s/\\\n[[:space:]]*/ /; ta' |
        grep -e "' _ \"\$PERL_IMAGE\"")
    printf '%s\n' "$_perl_args" | grep -q '"\$PERL_CACHE_VERSION"' &&
        ok "the Perl SV subshell receives \$PERL_CACHE_VERSION in its argument list" ||
        bad "the Perl SV subshell's argument list does not pass \$PERL_CACHE_VERSION"

    # Behavioural: run the child under set -u with a stub docker, one input, and read the
    # --cache_version the stub was handed.
    STUBBIN="$TMP/stubbin"
    mkdir -p "$STUBBIN" "$TMP/sv_in" "$TMP/sv_scratch"
    cat >"$STUBBIN/docker" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$@" >"${STUB_DOCKER_LOG:?}"
EOF
    chmod +x "$STUBBIN/docker"
    printf '##fileformat=VCFv4.2\n' | gzip -c >"$TMP/sv_in/one.vcf.gz"
    : >"$TMP/sv_in/manifest.json"
    (
        unset PERL_CACHE_VERSION
        export STUB_DOCKER_LOG="$TMP/docker_args.txt"
        PATH="$STUBBIN:$PATH" bash -c "$perl_child" _ "ensemblorg/ensembl-vep:release_116.2" "$TMP/cache" "$TMP/ref.fa" \
            GRCh37 "$TMP/sv_in" "$TMP/sv_scratch" "0:0" 5000 4 116
    ) >"$TMP/perl_child.out" 2>"$TMP/perl_child.err"
    rc_child=$?
    [[ "$rc_child" -eq 0 ]] &&
        ok "the child ran to completion under set -u (rc=0)" ||
        bad "the child exited $rc_child under set -u: $(head -1 "$TMP/perl_child.err")"
    grep -q 'unbound variable' "$TMP/perl_child.err" &&
        bad "the child aborted on an unbound variable: $(grep 'unbound variable' "$TMP/perl_child.err" | head -1)" ||
        ok "no unbound-variable abort"
    if [[ -f "$TMP/docker_args.txt" ]]; then
        _cv=$(grep -A 1 -x -- '--cache_version' "$TMP/docker_args.txt" | tail -1)
        [[ "$_cv" == "116" ]] &&
            ok "vep received --cache_version 116, the value passed as the tenth argument" ||
            bad "vep received --cache_version '$_cv'; expected 116"
    else
        bad "the stub docker was never invoked"
    fi
fi

echo
echo "REFERENCE-RELEASE: every comparator call scores against the release the ground truth names"

# The comparators default to the 116.2 registry. A cell scored against a 115.2 ground truth
# under it, or the reverse, applies the wrong divergence classes and reports an adjusted F1
# that means nothing. The release travels from the ground-truth directory's provenance (or
# the run's --reference-release) into every comparator call, and a ground truth that names
# neither an image tag nor a cache version fails the cell.
if [[ -f "$HARNESS" ]]; then
    _rr_src=$(awk '/^reference_release_for\(\) \{$/,/^}$/' "$HARNESS")
    if [[ -n "$_rr_src" ]]; then
        eval "$_rr_src"
        ok "reference_release_for is defined in the harness"
    else
        bad "reference_release_for not found in run_clone_measurement.sh"
    fi
    SCRIPT_DIR="$(dirname "$HARNESS")"
    GT116="$TMP/gt116"
    GT115="$TMP/gt115"
    GTNONE="$TMP/gtnone"
    mkdir -p "$GT116" "$GT115" "$GTNONE"
    printf '{"perl_image": "ensemblorg/ensembl-vep:release_116.2", "cache_version": 116}\n' >"$GT116/output.provenance.json"
    printf '{"perl_image": "ensemblorg/ensembl-vep@sha256:0000", "cache_version": 115}\n' >"$GT115/provenance.json"
    printf '{"fork": 4, "buffer_size": 5000}\n' >"$GTNONE/output.provenance.json"
    if declare -f reference_release_for >/dev/null; then
        REFERENCE_RELEASE=""
        got=$(reference_release_for "$GT116" 2>/dev/null)
        [[ "$got" == "116.2" ]] &&
            ok "a ground truth whose sidecar names release_116.2 scores under 116.2" ||
            bad "release_116.2 sidecar resolved to '$got'"
        got=$(reference_release_for "$GT115" 2>/dev/null)
        [[ "$got" == "115.2" ]] &&
            ok "a digest-only image with cache_version 115 scores under 115.2" ||
            bad "cache_version 115 resolved to '$got'"
        if reference_release_for "$GTNONE" >"$TMP/rr_none.out" 2>"$TMP/rr_none.err"; then
            bad "a ground truth naming no release passed; it would be scored under a guessed registry"
        else
            ok "a ground truth naming no release fails the cell"
        fi
        grep -q '^ERROR: \[run_clone_measurement\] no reference release for' "$TMP/rr_none.err" &&
            ok "the failure names the ground-truth directory on stderr" ||
            bad "the failure does not name the ground-truth directory: $(head -1 "$TMP/rr_none.err")"
        REFERENCE_RELEASE="115.2"
        got=$(reference_release_for "$GT116" 2>/dev/null)
        [[ "$got" == "115.2" ]] &&
            ok "--reference-release overrides the ground truth's own provenance" ||
            bad "the override was ignored (got '$got')"
        REFERENCE_RELEASE=""
        # A host without python3 cannot run the derivation (nor the comparators): the
        # failure names the interpreter, not the ground truth's provenance. PATH is set
        # as its own statement: a `PATH=... func` prefix leaves bash 3.2's hash table
        # intact, so `command -v python3` still finds the interpreter this shell hashed.
        mkdir -p "$TMP/nopython"
        if (PATH="$TMP/nopython"; reference_release_for "$GT116") >"$TMP/rr_nopy.out" 2>"$TMP/rr_nopy.err"; then
            bad "without python3 the derivation returned '$(cat "$TMP/rr_nopy.out")' instead of failing"
        elif grep -q '^ERROR: \[run_clone_measurement\] python3 is not on PATH' "$TMP/rr_nopy.err"; then
            ok "without python3 the cell fails naming the interpreter, not the provenance"
        else
            bad "without python3 the failure blames something else: $(head -1 "$TMP/rr_nopy.err")"
        fi
    fi
    # The three layouts a reference set is written in: one record at the set root over
    # suite directories that hold only outputs (the shape of the published SNP/indel set),
    # a sidecar beside each output, and neither.
    SET="$TMP/canonical_set"
    mkdir -p "$SET/s01" "$SET/s03" "$TMP/bare_set/s01"
    : >"$SET/s01/output.txt"
    : >"$SET/s01/output.txt_warnings.txt"
    : >"$SET/s03/output.txt"
    : >"$TMP/bare_set/s01/output.txt"
    printf '{"perl_image": "ensemblorg/ensembl-vep:release_115.2", "cache_version": 115}\n' >"$SET/provenance.json"
    printf '{"perl_image": "ensemblorg/ensembl-vep:release_116.2", "cache_version": 116}\n' >"$SET/s03/output.provenance.json"
    if declare -f reference_release_for >/dev/null; then
        REFERENCE_RELEASE=""
        got=$(reference_release_for "$SET/s01" 2>/dev/null)
        [[ "$got" == "115.2" ]] &&
            ok "layout 1, one record at the set root over output-only suite directories, resolves (115.2)" ||
            bad "the set-root record did not resolve for s01 (got '$got')"
        got=$(reference_release_for "$SET/s03" 2>/dev/null)
        [[ "$got" == "116.2" ]] &&
            ok "layout 2, a sidecar beside the suite's output, resolves and wins over the set record (116.2)" ||
            bad "the suite sidecar did not win (got '$got')"
        if reference_release_for "$TMP/bare_set/s01" >/dev/null 2>"$TMP/rr_bare.err"; then
            bad "layout 3, neither record, resolved to a release"
        else
            ok "layout 3, neither a set record nor a sidecar, fails the cell"
        fi
    fi
    _rr_calls=$(grep -c -e '--reference-release "\$reference_release"' "$HARNESS") || _rr_calls=0
    [[ "$_rr_calls" -eq 4 ]] &&
        ok "all four comparator calls pass --reference-release (4 of 4)" ||
        bad "$_rr_calls of 4 comparator calls pass --reference-release; the others score under the default registry"
    # The release is derived once per cell, BEFORE the warmup and timed runs, so a nameless
    # ground truth costs no annotation time: one derivation in each measure_* function,
    # each on a line before that function's first time_run call.
    for fn in measure_p1 measure_sv; do
        _fn_start=$(grep -n -e "^${fn}() {" "$HARNESS" | head -1 | cut -d: -f1)
        _fn_body=$(sed -n "/^${fn}() {/,/^}/p" "$HARNESS")
        _derive_rel=$(printf '%s\n' "$_fn_body" | grep -n -e 'reference_release=$(reference_release_for ' | cut -d: -f1)
        _time_rel=$(printf '%s\n' "$_fn_body" | grep -n -e 'res=$(time_run ' | head -1 | cut -d: -f1)
        _n_derive=$(printf '%s\n' "$_derive_rel" | grep -c .) || _n_derive=0
        if [[ "$_n_derive" -eq 1 && -n "$_time_rel" && "$_derive_rel" -lt "$_time_rel" ]]; then
            ok "$fn derives the release once, before its first timed run (line $((_fn_start + _derive_rel - 1)) < $((_fn_start + _time_rel - 1)))"
        else
            bad "$fn derives the release $_n_derive time(s) (relative line ${_derive_rel:-none}) against its first time_run at relative line ${_time_rel:-none}"
        fi
    done
fi

echo
echo "EXTRA-ARGS: --vep-extra-args reaches every vep-rs invocation through the wrapper and is recorded in provenance.json; an empty value appends nothing"

# The extra arguments travel through <output-dir>/vep-wrapper.sh, which VEP_BINARY is pointed
# at, so the SNP/indel call and the SV loop's `bash -c` subshell carry them without a change to
# either command line. A stub binary records the argument vector it was run with and writes one
# data line, enough for the SV output gate; the run is --walltime-only on one SNP/indel cell and
# one SV cell with two inputs, six invocations in all.
if [[ -f "$HARNESS" ]]; then
    grep -q '^VEP_EXTRA_ARGS=""$' "$HARNESS" &&
        ok "--vep-extra-args defaults to nothing (VEP_EXTRA_ARGS empty)" ||
        bad "VEP_EXTRA_ARGS does not default to empty"
    grep -A 2 -e '^[[:space:]]*--vep-extra-args)$' "$HARNESS" | grep -q 'VEP_EXTRA_ARGS="$2"' &&
        ok "--vep-extra-args is parsed into VEP_EXTRA_ARGS" ||
        bad "--vep-extra-args is not parsed"
    grep -q -- '--vep-extra-args applies to --engine vep-rs only' "$HARNESS" &&
        ok "a non-vep-rs engine with extra arguments is refused" ||
        bad "the vep-rs-only guard on --vep-extra-args is missing"
    # The two engine command lines are not edited for the pass-through: the extra arguments
    # appear in the wrapper the harness writes, never in the SNP/indel call or the SV subshell.
    _extra_in_cmds=$(sed -n '/^measure_p1() {/,/^}/p; /^measure_sv() {/,/^}/p' "$HARNESS" |
        grep -cE 'VEP_EXTRA_ARG[SV]') || _extra_in_cmds=0
    [[ "$_extra_in_cmds" -eq 0 ]] &&
        ok "neither measure_p1 nor measure_sv names the extra arguments (they travel through VEP_BINARY)" ||
        bad "measure_p1 or measure_sv names VEP_EXTRA_ARGS on $_extra_in_cmds line(s); the pass-through must not touch the command lines"
    _ea="$TMP/extra_args"
    _ea_sv_inputs="$_ea/data/ground_truth/perl/sv_per_vcf/canonical_inputs/grch37"
    mkdir -p "$_ea/bin" "$_ea/data/vcf/grch37/all_variants_canonical" "$_ea_sv_inputs"
    : >"$_ea/data/vcf/grch37/all_variants_canonical/clinvar_grch37.canonical.vcf.gz"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t1\n' >"$_ea_sv_inputs/01_a.vcf"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t2\n' >"$_ea_sv_inputs/02_b.vcf"
    cat >"$_ea/bin/vep" <<'STUB'
#!/usr/bin/env bash
# stub engine: one line per invocation, the argument vector joined by spaces; one data line at -o
out=""; prev=""
for a in "$@"; do [[ "$prev" == -o ]] && out="$a"; prev="$a"; done
printf '%s\n' "$*" >>"$STUB_ARGV"
printf '#Uploaded_variation\tLocation\nv\t21:1\n' >"$out"
STUB
    chmod +x "$_ea/bin/vep"
    run_extra_args_harness() { # run_extra_args_harness <run name> [harness args...]
        local run_name="$1"; shift
        : >"$_ea/$run_name.argv"
        STUB_ARGV="$_ea/$run_name.argv" bash "$HARNESS" --engine vep-rs --vep-binary "$_ea/bin/vep" --data-dir "$_ea/data" \
            --output-dir "$_ea/$run_name" --suites "s01 s07" --walltime-only "$@" >"$_ea/$run_name.log" 2>&1
    }
    if /usr/bin/time -v true >/dev/null 2>&1; then
        run_extra_args_harness with --vep-extra-args "--transcript_index_impl coitrees --fork 4"; rc=$?
        _n_inv=$(grep -c . "$_ea/with.argv") || _n_inv=0
        _n_tail=$(grep -c -- ' --quiet --transcript_index_impl coitrees --fork 4$' "$_ea/with.argv") || _n_tail=0
        _n_sv=$(grep -c -- '--species homo_sapiens --assembly GRCh37 --quiet --transcript_index_impl coitrees --fork 4$' "$_ea/with.argv") || _n_sv=0
        _n_p1=$(grep -c -- '/s01/rust/output.txt --offline .* --force --no_stats --quiet --transcript_index_impl coitrees --fork 4$' "$_ea/with.argv") || _n_p1=0
        if [[ "$rc" -eq 0 && "$_n_inv" -eq 6 && "$_n_tail" -eq 6 ]]; then
            ok "all six vep-rs invocations (SNP/indel warmup and timed, SV loop twice over two inputs) end with the extra arguments after the harness's own"
        else
            bad "rc=$rc; $_n_inv invocation(s), $_n_tail of them ending with the extra arguments (expected 6 and 6): $(head -2 "$_ea/with.argv" | cut -c1-160 | tr '\n' '|') $(tail -3 "$_ea/with.log" | tr '\n' '|')"
        fi
        [[ "$_n_p1" -eq 2 && "$_n_sv" -eq 4 ]] &&
            ok "the SNP/indel call (2) and the SV subshell's per-VCF calls (4) both ran through the wrapper" ||
            bad "$_n_p1 SNP/indel and $_n_sv SV invocations carried the extra arguments (expected 2 and 4)"
        _bad_rows=$(awk -F, 'NR>1 && ($10!="0" || ($8!="warmup" && $8!="timed"))' "$_ea/with/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        _rows=$(awk -F, 'NR>1' "$_ea/with/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$_rows" -eq 4 && "$_bad_rows" -eq 0 ]] &&
            ok "both cells record their warmup and timed rows with exit 0 through the wrapper" ||
            bad "wall_times.csv has $_rows rows, $_bad_rows of them failed (expected 4 and 0)"
        [[ -x "$_ea/with/vep-wrapper.sh" ]] && grep -q "^exec $_ea/bin/vep \"\$@\" --transcript_index_impl coitrees --fork 4\$" "$_ea/with/vep-wrapper.sh" &&
            ok "the wrapper execs the real binary with the harness's arguments first and the extra ones after" ||
            bad "vep-wrapper.sh is missing or does not exec the real binary with the extra arguments appended: $(cat "$_ea/with/vep-wrapper.sh" 2>/dev/null | tail -1)"
        _sha=$( (sha256sum "$_ea/bin/vep" 2>/dev/null || shasum -a 256 "$_ea/bin/vep") | awk '{print $1}')
        if python3 - "$_ea/with/provenance.json" "$_ea/bin/vep" "$_sha" "$_ea/with/vep-wrapper.sh" <<'PY'
import json, sys
p, binary, sha, wrapper = sys.argv[1:]
d = json.load(open(p))
assert d["engine"] == "vep-rs", d
assert d["vep_binary"] == binary and d["vep_binary_sha256"] == sha, (d["vep_binary"], d["vep_binary_sha256"])
assert d["vep_extra_args"] == "--transcript_index_impl coitrees --fork 4", d["vep_extra_args"]
assert d["vep_extra_argv"] == ["--transcript_index_impl", "coitrees", "--fork", "4"], d["vep_extra_argv"]
assert d["wrapper"] == wrapper, d["wrapper"]
assert len(d["executed_argv"]) == 6, len(d["executed_argv"])
for argv in d["executed_argv"]:
    assert argv[0] == binary and argv[-4:] == ["--transcript_index_impl", "coitrees", "--fork", "4"], argv
    assert "--fork" in argv[:-4] and "--offline" in argv, argv
assert sum(1 for a in d["executed_argv"] if "--species" in a) == 4, "SV invocations"
PY
        then
            ok "provenance.json names the real binary, its sha256, the extra arguments as given and split, the wrapper, and the six executed argument vectors"
        else
            bad "provenance.json does not carry the expected record: $(cat "$_ea/with/provenance.json" 2>/dev/null | tr -d '\n' | cut -c1-300)"
        fi
        run_extra_args_harness empty --vep-extra-args ""; rc=$?
        _n_inv=$(grep -c . "$_ea/empty.argv") || _n_inv=0
        _n_plain=$(grep -cE -- ' --quiet$' "$_ea/empty.argv") || _n_plain=0
        if [[ "$rc" -eq 0 && "$_n_inv" -eq 6 && "$_n_plain" -eq 6 && ! -e "$_ea/empty/vep-wrapper.sh" ]]; then
            ok "an empty value appends nothing: six invocations end with the harness's own --quiet and no wrapper is written"
        else
            bad "empty value: rc=$rc, $_n_inv invocations, $_n_plain ending with --quiet, wrapper present=$([[ -e "$_ea/empty/vep-wrapper.sh" ]] && echo yes || echo no)"
        fi
        if python3 - "$_ea/empty/provenance.json" "$_ea/bin/vep" 2>/dev/null <<'PY'
import json, sys
p, binary = sys.argv[1:]
d = json.load(open(p))
assert d["vep_binary"] == binary, d
assert d["wrapper"] is None, d
assert d["vep_extra_args"] == "" and d["vep_extra_argv"] == [], d
assert d["executed_argv"] == [], d
PY
        then
            ok "provenance.json is still written for a run without extra arguments, naming the binary and no wrapper"
        else
            bad "provenance.json for the empty value is missing or names a wrapper: $(cat "$_ea/empty/provenance.json" 2>/dev/null | tr -d '\n' | cut -c1-200)"
        fi
        # "Split on whitespace" includes a newline: a value read line by line would keep only
        # its first line and drop the rest without a word.
        run_extra_args_harness lines --vep-extra-args $'--transcript_index_impl coitrees\n--fork 4'; rc=$?
        _n_tail=$(grep -c -- ' --quiet --transcript_index_impl coitrees --fork 4$' "$_ea/lines.argv") || _n_tail=0
        if [[ "$rc" -eq 0 && "$_n_tail" -eq 6 ]] && python3 - "$_ea/lines/provenance.json" 2>/dev/null <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
assert d["vep_extra_argv"] == ["--transcript_index_impl", "coitrees", "--fork", "4"], d
PY
        then
            ok "a newline inside the value separates arguments like a space: all six invocations carry both flags"
        else
            bad "newline-separated value: rc=$rc, $_n_tail of 6 invocations carry both flags: $(head -1 "$_ea/lines.argv" 2>/dev/null | cut -c1-160)"
        fi
        # A binary named without a slash: the executable test reads the working directory
        # while a bare command word resolves through PATH, so a decoy first on PATH would run
        # and the record would hash the file that did not. The harness makes the operand
        # absolute, so the working directory's file runs and is the one recorded.
        mkdir -p "$_ea/decoy"
        printf '#!/usr/bin/env bash\n: >"$DECOY_MARK"\n' >"$_ea/decoy/vep"
        chmod +x "$_ea/decoy/vep"
        : >"$_ea/rel.argv"
        (cd "$_ea/bin" && DECOY_MARK="$_ea/decoy.ran" STUB_ARGV="$_ea/rel.argv" PATH="$_ea/decoy:$PATH" bash "$HARNESS" --engine vep-rs --vep-binary vep \
            --data-dir "$_ea/data" --output-dir "$_ea/rel" --suites s01 --walltime-only >"$_ea/rel.log" 2>&1); rc=$?
        _n_inv=$(grep -c . "$_ea/rel.argv") || _n_inv=0
        if [[ "$rc" -eq 0 && "$_n_inv" -eq 2 && ! -e "$_ea/decoy.ran" ]] && python3 - "$_ea/rel/provenance.json" "$_ea/bin/vep" 2>/dev/null <<'PY'
import json, sys
p, binary = sys.argv[1:]
d = json.load(open(p))
assert d["vep_binary"] == binary, d["vep_binary"]
PY
        then
            ok "a slash-less --vep-binary runs the working directory's file, never a PATH decoy, and is recorded by its absolute path"
        else
            bad "slash-less --vep-binary: rc=$rc, $_n_inv stub invocation(s), decoy ran=$([[ -e "$_ea/decoy.ran" ]] && echo yes || echo no), recorded as $(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["vep_binary"])' "$_ea/rel/provenance.json" 2>/dev/null)"
        fi
        # A run without extra arguments into the directory of a run with them leaves no
        # wrapper or log behind to contradict its own record.
        run_extra_args_harness rerun --vep-extra-args "--fork 4"; rc1=$?
        run_extra_args_harness rerun; rc2=$?
        if [[ "$rc1" -eq 0 && "$rc2" -eq 0 && ! -e "$_ea/rerun/vep-wrapper.sh" && ! -e "$_ea/rerun/vep-wrapper.argv" ]] && python3 - "$_ea/rerun/provenance.json" 2>/dev/null <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
assert d["wrapper"] is None and d["executed_argv"] == [], d
PY
        then
            ok "a re-run without extra arguments removes the earlier run's wrapper and log from the output directory"
        else
            bad "re-run into the same output directory: rc=$rc1/$rc2, left behind: $(ls "$_ea/rerun" 2>/dev/null | grep vep-wrapper | tr '\n' ' ')"
        fi
        # Only a vep-rs run writes the record, so a run under another engine into the same
        # directory must not leave a vep-rs record beside its own wall_times.csv. The record's
        # presence before the run is part of the check: its absence alone also holds when the
        # vep-rs run above wrote none.
        _record_before=$([[ -e "$_ea/rerun/provenance.json" ]] && echo yes || echo no)
        STUB_ARGV="$_ea/rerun-fastvep.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_ea/bin/vep" --data-dir "$_ea/data" --output-dir "$_ea/rerun" --suites s01 --walltime-only >"$_ea/rerun-fastvep.log" 2>&1; rc=$?
        if [[ "$_record_before" == yes && "$rc" -eq 0 && ! -e "$_ea/rerun/provenance.json" ]] && grep -q '^fastvep,' "$_ea/rerun/wall_times.csv"; then
            ok "a run under another engine into the same output directory removes the earlier vep-rs record"
        else
            bad "fastvep re-run into a vep-rs output directory: rc=$rc, provenance.json before=$_record_before, after=$([[ -e "$_ea/rerun/provenance.json" ]] && echo yes || echo no): $(tail -2 "$_ea/rerun-fastvep.log" | tr '\n' '|' | cut -c1-200)"
        fi
    else
        echo "  skip functional extra-args run: /usr/bin/time -v is not GNU time here"
    fi
    # A directory passes `-x`; the harness refuses it with its own line before anything is
    # written, instead of failing inside the record writer.
    bash "$HARNESS" --engine vep-rs --vep-binary "$_ea/bin" --data-dir "$_ea/data" --output-dir "$_ea/dirbin" --suites s01 --walltime-only >"$_ea/dirbin.log" 2>&1; rc=$?
    if [[ "$rc" -eq 1 && ! -e "$_ea/dirbin" ]] && grep -q '^ERROR: vep-rs binary not a readable executable file: ' "$_ea/dirbin.log" && ! grep -q 'Traceback' "$_ea/dirbin.log"; then
        ok "a directory given as --vep-binary is refused before any cell with the harness's own line"
    else
        bad "directory as --vep-binary: rc=$rc, output dir exists=$([[ -e "$_ea/dirbin" ]] && echo yes || echo no): $(head -2 "$_ea/dirbin.log" | tr '\n' '|' | cut -c1-200)"
    fi
    # An execute-only file passes `-x` and `-f` and cannot be hashed; it is refused like the
    # directory. Root reads every file, so the check is meaningful only for another user.
    if [[ "$EUID" -ne 0 ]]; then
        cp "$_ea/bin/vep" "$_ea/bin/vep-xonly"
        chmod 0111 "$_ea/bin/vep-xonly"
        bash "$HARNESS" --engine vep-rs --vep-binary "$_ea/bin/vep-xonly" --data-dir "$_ea/data" --output-dir "$_ea/xonly" --suites s01 --walltime-only >"$_ea/xonly.log" 2>&1; rc=$?
        if [[ "$rc" -eq 1 && ! -e "$_ea/xonly" ]] && grep -q '^ERROR: vep-rs binary not a readable executable file: ' "$_ea/xonly.log" && ! grep -q 'Traceback' "$_ea/xonly.log"; then
            ok "an execute-only --vep-binary is refused before any cell with the harness's own line"
        else
            bad "execute-only --vep-binary: rc=$rc, output dir exists=$([[ -e "$_ea/xonly" ]] && echo yes || echo no): $(head -2 "$_ea/xonly.log" | tr '\n' '|' | cut -c1-200)"
        fi
        chmod 0755 "$_ea/bin/vep-xonly"
    else
        echo "  skip execute-only check: root reads every file"
    fi
    # Every vep-rs run writes provenance.json through python3, so a host without it is refused
    # at startup by the harness's own line. The PATH holds the tools the harness runs before
    # the provenance write and no interpreter.
    _np="$_ea/nopython/bin"
    mkdir -p "$_np"
    for _t in bash dirname date mkdir; do ln -sf "$(type -P "$_t")" "$_np/$_t"; done
    PATH="$_np" "$_np/bash" "$HARNESS" --engine vep-rs --vep-binary "$_ea/bin/vep" --data-dir "$_ea/data" --output-dir "$_ea/nopython/out" --suites s01 --walltime-only >"$_ea/nopython.log" 2>&1; rc=$?
    _n_err=$(grep -c '^ERROR: \[run_clone_measurement\] python3 is not on PATH' "$_ea/nopython.log") || _n_err=0
    if [[ "$rc" -eq 1 && "$_n_err" -eq 1 ]] && ! grep -q -e 'command not found' -e ' starting$' "$_ea/nopython.log"; then
        ok "without python3 a vep-rs run stops before any cell with one line naming the interpreter"
    else
        bad "without python3: rc=$rc, $_n_err interpreter line(s): $(head -3 "$_ea/nopython.log" | tr '\n' '|')"
    fi
    bash "$HARNESS" --engine perl --data-dir "$_ea/data" --output-dir "$_ea/perl" --suites s01 --vep-extra-args "--transcript_index_impl coitrees" >"$_ea/perl.log" 2>&1; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: --vep-extra-args applies to --engine vep-rs only (got engine 'perl')" "$_ea/perl.log" && ! grep -q "run_clone_measurement engine=perl" "$_ea/perl.log"; then
        ok "extra arguments under another engine stop the run before any cell"
    else
        bad "extra arguments under --engine perl were not refused before the run (rc=$rc): $(head -2 "$_ea/perl.log" | tr '\n' ' ')"
    fi
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
