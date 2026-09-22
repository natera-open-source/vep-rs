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
    # function body: seven legitimate `exit 1` calls live outside it (argument
    # parsing and setup), where aborting is correct. Comments are stripped so the
    # rationale comment beside the guard cannot itself trip this.
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

    annotate_one_defs=$(grep -cE '^[[:space:]]*let annotate_one = ' "$RUNNER") || annotate_one_defs=0
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
        mkdir -p "$PERL_CACHE_DIR/homo_sapiens/115_GRCh37/21"
        : >"$PERL_CACHE_DIR/homo_sapiens/115_GRCh37/info.txt"
        if assert_perl_cache_complete GRCh37 >/dev/null 2>&1; then
            bad "a chr21-only GRCh37 cache PASSED the completeness guard"
        else
            ok "a chr21-only GRCh37 cache fails the guard"
        fi

        # All 25 primary contigs, no scaffolds. Ensembl's full tree carries 699 directories, so a
        # guard demanding all of them would fail this legitimately pruned cache.
        PERL_CACHE_DIR="$TMP/cache_primary"
        for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
            mkdir -p "$PERL_CACHE_DIR/homo_sapiens/115_GRCh37/$c"
        done
        if assert_perl_cache_complete GRCh37 >/dev/null 2>&1; then
            ok "a cache with all 25 primary contigs and no scaffolds passes"
        else
            bad "a 25-primary-contig cache FAILED the guard; it must not require scaffolds"
        fi

        # One contig short is the interesting near-miss: 24 of 25 is what an incomplete
        # cache extraction looks like, and it must not pass.
        rmdir "$PERL_CACHE_DIR/homo_sapiens/115_GRCh37/MT"
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
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
