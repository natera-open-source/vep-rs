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
#   IMAGE-DIGEST   a Perl image whose pulled digest is not the pinned one, which
#                  must stop the run before any cell is timed under that
#                  engine_version.
#   CELL-LIST      a cell id the suite table does not carry, which would match
#                  nothing, time nothing and exit 0, or a repeated one; --cells
#                  (a comma-separated list) and --suites refuse both before any
#                  cell runs.
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
#   FASTVEP-EXTRA-ARGS  --fastvep-extra-args that miss one of the fastVEP `annotate`
#                  runs (the SNP/indel warmup or timed run, an input of the SV loop),
#                  a word the rule does not admit reaching the binary, the words
#                  under another engine, a run without them that writes a wrapper or
#                  a record, a --sa-dir directory run over cells of two assemblies,
#                  or a path holding `=` handed to `env`, which reads it as an
#                  assignment: the wrapper's output directory with words, the
#                  binary's own path without them; or a binary named without a
#                  slash tested in the working directory while the cells start
#                  the file PATH resolves.
#   REFUSED-RUN    a run refused for its cell list or its --sa-dir cells after the
#                  output directory is made, which resets the wall_times.csv an
#                  earlier run left there and removes its record, wrapper and log.
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
        if assert_perl_cache_complete homo_sapiens GRCh37 >/dev/null 2>&1; then
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
        if assert_perl_cache_complete homo_sapiens GRCh37 >/dev/null 2>&1; then
            ok "a cache with all 25 primary contigs and no scaffolds passes"
        else
            bad "a 25-primary-contig cache FAILED the guard; it must not require scaffolds"
        fi

        # One contig short is the interesting near-miss: 24 of 25 is what an incomplete
        # cache extraction looks like, and it must not pass.
        rmdir "$PERL_CACHE_DIR/homo_sapiens/116_GRCh37/MT"
        if assert_perl_cache_complete homo_sapiens GRCh37 >/dev/null 2>&1; then
            bad "a cache missing only MT PASSED the guard"
        else
            ok "a cache missing a single primary contig fails the guard"
        fi

        # A missing slice must fail loudly rather than pass by absence: a check that measured
        # nothing reports clean.
        PERL_CACHE_DIR="$TMP/cache_absent"
        if assert_perl_cache_complete homo_sapiens GRCh37 >/dev/null 2>&1; then
            bad "an ABSENT cache slice PASSED the guard"
        else
            ok "an absent cache slice fails the guard rather than passing by omission"
        fi

        # Another species' slice lives under its own directory; its primary set is not known,
        # so presence of a chromosome directory is what the guard asserts there.
        PERL_CACHE_DIR="$TMP/cache_yeast"
        mkdir -p "$PERL_CACHE_DIR/saccharomyces_cerevisiae/116_R64-1-1"
        if assert_perl_cache_complete saccharomyces_cerevisiae R64-1-1 >/dev/null 2>&1; then
            bad "an EMPTY yeast slice PASSED the guard"
        else
            ok "a slice of another species with no chromosome directory fails the guard"
        fi
        mkdir -p "$PERL_CACHE_DIR/saccharomyces_cerevisiae/116_R64-1-1/I"
        if assert_perl_cache_complete saccharomyces_cerevisiae R64-1-1 >/dev/null 2>&1; then
            ok "a slice of another species holding a chromosome directory passes under <species>/<version>_<assembly>"
        else
            bad "a yeast slice with a chromosome directory FAILED the guard"
        fi
        if assert_perl_cache_complete homo_sapiens R64-1-1 >/dev/null 2>&1; then
            bad "the yeast slice PASSED the guard under homo_sapiens"
        else
            ok "the slice is looked up under the species given, not under homo_sapiens"
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
    # A scanner that tracks `case` arms by text is satisfied by a comment naming the guard,
    # passes with the guard moved after the timing loop, and fails on a reindented `perl)`
    # label, so the order is read from the line numbers of the anchored calls.
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
    # The F1 compare runs the resolved tuple comparator command, `"${OUTPUTS_COMPARE[@]}"`.
    _f1_line=$(grep -n -e '"\${OUTPUTS_COMPARE\[@\]}"' "$HARNESS" | head -1 | cut -d: -f1)
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

# VEP_OUTPUTS_COMPARATOR selects the tuple comparator and nothing else: the script-level
# resolution block (OUTPUTS_COMPARE) is extracted and run with the same stubs. Unset or empty,
# the command is the tree's own release build of tools/vep-outputs-compare; `python` is the
# only value that runs the Python; an executable file is the command on its own, a bare name
# resolving through PATH; a comparator that is not an executable file, the unbuilt default
# included, stops the run with one line naming both remedies and never runs the Python; a run
# that scores no F1 (--walltime-only, the Perl engine) asks for none; the block stands before
# the first engine run and both F1 call sites run the resolved command.
if [[ -f "$HARNESS" ]]; then
    _outputs_block=$(awk '
        !p && /^OUTPUTS_COMPARE=\(\)$/ { p = 1 }
        p { print }
        p && /^fi$/ { exit }
    ' "$HARNESS")
    [[ "$_outputs_block" == *VEP_OUTPUTS_COMPARATOR* && "$_outputs_block" == *$'\nfi' ]] &&
        ok "the OUTPUTS_COMPARE resolution block is extracted whole" ||
        bad "could not extract the OUTPUTS_COMPARE resolution block from the harness"
    _tool_rel="tools/vep-outputs-compare/target/release/vep-outputs-compare"
    _built="$_fdir/repo/$_tool_rel"
    _unbuilt="$_fdir/repo-unbuilt/$_tool_rel"
    mkdir -p "$(dirname "$_built")"
    cp "$_fdir/bin/comparator" "$_built"
    # $1 is the repository root, $2 the VEP_OUTPUTS_COMPARATOR value ("" for unset, "empty"
    # for the empty string), $3 the engine and $4 the WALLTIME_ONLY value; prints the block's
    # stderr, then the resolved command one word per line, then `rc=<exit status>`.
    _run_outputs_block() {
        (
            PATH="$_fdir/bin:$PATH"
            REPO_DIR="$1"
            ENGINE="${3:-vep-rs}"
            WALLTIME_ONLY="${4:-false}"
            case "$2" in
            "") unset VEP_OUTPUTS_COMPARATOR ;;
            empty) export VEP_OUTPUTS_COMPARATOR="" ;;
            *) export VEP_OUTPUTS_COMPARATOR="$2" ;;
            esac
            eval "$_outputs_block"
            [[ ${#OUTPUTS_COMPARE[@]} -eq 0 ]] || printf '%s\n' "${OUTPUTS_COMPARE[@]}"
        ) 2>&1
        echo "rc=$?"
    }
    for _unset in "" empty; do
        _out=$(_run_outputs_block "$_fdir/repo" "$_unset")
        [[ "$_out" == "$_built"$'\n'"rc=0" ]] &&
            ok "VEP_OUTPUTS_COMPARATOR ${_unset:-unset}: the tree's release build of the Rust tool is the command" ||
            bad "VEP_OUTPUTS_COMPARATOR ${_unset:-unset}: expected $_built alone, got '$_out'"
    done
    _out=$(_run_outputs_block "$_fdir/repo" python)
    [[ "$_out" == "python3"$'\n'"$_fdir/repo/scripts/concordance/compare_vep_outputs.py"$'\n'"rc=0" ]] &&
        ok "VEP_OUTPUTS_COMPARATOR=python: the Python comparator is the command" ||
        bad "VEP_OUTPUTS_COMPARATOR=python: expected python3 and the script path, got '$_out'"
    _out=$(_run_outputs_block "$_fdir/repo-unbuilt" python)
    [[ "$_out" == "python3"$'\n'"$_fdir/repo-unbuilt/scripts/concordance/compare_vep_outputs.py"$'\n'"rc=0" ]] &&
        ok "VEP_OUTPUTS_COMPARATOR=python: the Python runs in a tree whose Rust tool is not built" ||
        bad "VEP_OUTPUTS_COMPARATOR=python in an unbuilt tree: expected the Python command, got '$_out'"
    _out=$(_run_outputs_block "$_fdir/repo-unbuilt" "$_fdir/bin/comparator")
    [[ "$_out" == "$_fdir/bin/comparator"$'\n'"rc=0" ]] &&
        ok "VEP_OUTPUTS_COMPARATOR executable: the binary is the command, run as given" ||
        bad "VEP_OUTPUTS_COMPARATOR executable: expected the comparator stub alone, got '$_out'"
    _out=$(_run_outputs_block "$_fdir/repo-unbuilt" comparator)
    [[ "$_out" == "$_fdir/bin/comparator"$'\n'"rc=0" ]] &&
        ok "VEP_OUTPUTS_COMPARATOR=comparator (a bare name on PATH): the binary is the command, resolved through PATH" ||
        bad "VEP_OUTPUTS_COMPARATOR=comparator: expected the comparator stub through PATH, got '$_out'"
    # The one line a comparator that cannot run is refused with: the value it was asked for,
    # then both remedies.
    _refusal() {
        printf '%s' "ERROR: [run_clone_measurement] the tuple comparator $1 is not an executable: build the tool with 'cargo build --release --locked --manifest-path $2/tools/vep-outputs-compare/Cargo.toml' or name a built one in VEP_OUTPUTS_COMPARATOR, or set VEP_OUTPUTS_COMPARATOR=python to score with the Python reference implementation"
    }
    _out=$(_run_outputs_block "$_fdir/repo-unbuilt" "")
    [[ "$_out" == "$(_refusal "$_unbuilt" "$_fdir/repo-unbuilt")"$'\n'"rc=1" ]] &&
        ok "VEP_OUTPUTS_COMPARATOR unset, the tool not built: one line names both remedies, exit 1, the Python does not run" ||
        bad "unbuilt default: expected the refusal alone and rc=1, got '$_out'"
    for _bad_value in "$_fdir/not_executable" "$_fdir" inert; do
        _out=$(_run_outputs_block "$_fdir/repo" "$_bad_value")
        [[ "$_out" == "$(_refusal "$_bad_value" "$_fdir/repo")"$'\n'"rc=1" ]] &&
            ok "VEP_OUTPUTS_COMPARATOR=${_bad_value##*/} (not an executable file): refused with both remedies, exit 1, the Python does not run" ||
            bad "VEP_OUTPUTS_COMPARATOR=$_bad_value: expected the refusal alone and rc=1, got '$_out'"
    done
    _out=$(_run_outputs_block "$_fdir/repo-unbuilt" "" vep-rs true)
    [[ "$_out" == "rc=0" ]] &&
        ok "--walltime-only: no comparator is asked for, so an unbuilt tool does not stop the run" ||
        bad "--walltime-only with an unbuilt tool: expected no command and rc=0, got '$_out'"
    _out=$(_run_outputs_block "$_fdir/repo-unbuilt" "" perl false)
    [[ "$_out" == "rc=0" ]] &&
        ok "the Perl engine scores no F1: no comparator is asked for" ||
        bad "the Perl engine with an unbuilt tool: expected no command and rc=0, got '$_out'"
    _out=$(_run_outputs_block "$_fdir/repo-unbuilt" "" fastvep false)
    [[ "$_out" == "$(_refusal "$_unbuilt" "$_fdir/repo-unbuilt")"$'\n'"rc=1" ]] &&
        ok "the fastVEP engine scores F1: an unbuilt tool stops it as it stops vep-rs" ||
        bad "the fastVEP engine with an unbuilt tool: expected the refusal and rc=1, got '$_out'"
    _block_line=$(grep -n -e '^OUTPUTS_COMPARE=()$' "$HARNESS" | head -1 | cut -d: -f1)
    _cells_line=$(grep -n -e '^for entry in "\${SUITES\[@\]}"; do$' "$HARNESS" | head -1 | cut -d: -f1)
    [[ -n "$_block_line" && -n "$_cells_line" && "$_block_line" -lt "$_cells_line" ]] &&
        ok "the comparator is resolved before the loop that runs the cells (line $_block_line before $_cells_line)" ||
        bad "the resolution block must stand before the cell loop (block at ${_block_line:-?}, loop at ${_cells_line:-?})"
    _outputs_calls=$(grep -c -e '^[[:space:]]*"\${OUTPUTS_COMPARE\[@\]}" \\$' "$HARNESS") || _outputs_calls=0
    _python_runs=$(grep -c -e 'python3 "\$REPO_DIR/scripts/concordance/compare_vep_outputs\.py"' "$HARNESS") || _python_runs=0
    [[ "$_outputs_calls" -eq 2 && "$_python_runs" -eq 1 ]] &&
        ok "both F1 call sites (vep-rs, fastvep) run the resolved command and the Python is run from the one line that selects it" ||
        bad "expected 2 OUTPUTS_COMPARE call sites and 1 line that runs the Python, found $_outputs_calls and $_python_runs"
    # The whole harness, not the block: a vep-rs run that would score F1 in a tree without
    # the tool stops on the refusal before it makes its output directory or runs the engine.
    _od="$TMP/outputs_default"
    mkdir -p "$_od/tree/scripts/concordance" "$_od/bin" "$_od/data"
    cp "$HARNESS" "$_od/tree/scripts/concordance/"
    printf '#!/bin/bash\n: >"%s"\n' "$_od/engine.ran" >"$_od/bin/vep"
    chmod +x "$_od/bin/vep"
    (unset VEP_OUTPUTS_COMPARATOR; bash "$_od/tree/scripts/concordance/run_clone_measurement.sh" --engine vep-rs --vep-binary "$_od/bin/vep" --data-dir "$_od/data" --output-dir "$_od/out" --suites s01) >"$_od/log" 2>&1
    rc=$?
    [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_od/log")" -eq 1 && "$(tail -1 "$_od/log")" == "$(_refusal "$_od/tree/$_tool_rel" "$_od/tree")" && ! -e "$_od/engine.ran" ]] &&
        ok "a harness run without the tool: exit 1 on the one refusal line, the engine never started" ||
        bad "harness without the tool: expected rc=1, the refusal as the one ERROR and last line and no engine run, got rc=$rc '$(tail -2 "$_od/log")'"
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
echo "IMAGE-DIGEST: the Perl image must carry the pinned digest before any cell runs"

# The engine_version every Perl row carries is the image tag, and a tag can move. The harness
# pulls the image when absent and refuses one whose RepoDigests lack the expected digest, so a
# moved tag cannot be timed under the pinned version. Run against a stub docker whose inspect
# answers with a configurable digest; an empty data tree is enough, because the check precedes
# every cell.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"
if [[ -f "$HARNESS" ]]; then
    PINNED_DIGEST="sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da"
    OTHER_DIGEST="sha256:0000000000000000000000000000000000000000000000000000000000000000"
    _dd="$TMP/digest"; mkdir -p "$_dd/bin" "$_dd/data"
    cat >"$_dd/bin/docker" <<'STUB'
#!/usr/bin/env bash
echo "docker $*" >>"$DIGEST_CALLS"
case "$1" in
image) echo "ensemblorg/ensembl-vep@$STUB_DIGEST"; exit 0 ;;
pull) exit 0 ;;
*) exit 0 ;;
esac
STUB
    chmod +x "$_dd/bin/docker"
    run_perl_harness() { # run_perl_harness <stub digest> <out dir> [harness args...]
        local d="$1" o="$2"; shift 2
        : >"$_dd/calls"
        DIGEST_CALLS="$_dd/calls" STUB_DIGEST="$d" PATH="$_dd/bin:$PATH" \
            bash "$HARNESS" --engine perl --data-dir "$_dd/data" --output-dir "$o" --suites s01 "$@" >"$o.log" 2>&1
    }
    grep -q '^PERL_IMAGE="ensemblorg/ensembl-vep:release_116.2"$' "$HARNESS" &&
        ok "the default Perl image is the release_116.2 tag" ||
        bad "the default Perl image is not ensemblorg/ensembl-vep:release_116.2"
    grep -q "^PERL_IMAGE_DIGEST=\"\${VEP_PERL_IMAGE_DIGEST:-$PINNED_DIGEST}\"$" "$HARNESS" &&
        ok "the expected digest is one constant, overridable by VEP_PERL_IMAGE_DIGEST" ||
        bad "PERL_IMAGE_DIGEST is not the pinned manifest-list digest with the env override"
    run_perl_harness "$OTHER_DIGEST" "$_dd/out1"; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: \[perl-image\] digest mismatch: ensemblorg/ensembl-vep@$OTHER_DIGEST (expected $PINNED_DIGEST)$" "$_dd/out1.log"; then
        ok "a pulled image with another digest stops the run with the digest found and the digest expected"
    else
        bad "a mismatched digest did not stop the run (rc=$rc): $(head -3 "$_dd/out1.log" | tr '\n' ' ')"
    fi
    if ! grep -q "run_clone_measurement engine=perl" "$_dd/out1.log" && ! grep -q "^docker run" "$_dd/calls"; then
        ok "the mismatch is caught before the first suite and before any docker run"
    else
        bad "the harness went on past the digest mismatch"
    fi
    run_perl_harness "$PINNED_DIGEST" "$_dd/out2"; rc=$?
    if grep -q "run_clone_measurement engine=perl" "$_dd/out2.log" && ! grep -q "digest mismatch" "$_dd/out2.log"; then
        ok "the pinned digest passes the check and the suites start"
    else
        bad "the pinned digest did not pass the check: $(head -3 "$_dd/out2.log" | tr '\n' ' ')"
    fi
    grep -q '^docker image inspect --format {{join .RepoDigests ","}} ensemblorg/ensembl-vep:release_116.2$' "$_dd/calls" &&
        ok "the check reads every repo digest of the pinned tag" ||
        bad "the check did not inspect the pinned tag's RepoDigests"
    run_perl_harness "$OTHER_DIGEST" "$_dd/out3" --perl-image-digest "$OTHER_DIGEST"; rc=$?
    if grep -q "run_clone_measurement engine=perl" "$_dd/out3.log" && ! grep -q "digest mismatch" "$_dd/out3.log"; then
        ok "--perl-image-digest names the digest another image must carry"
    else
        bad "--perl-image-digest did not override the expected digest: $(head -3 "$_dd/out3.log" | tr '\n' ' ')"
    fi
    run_perl_harness "$PINNED_DIGEST" "$_dd/out4" --perl-image-digest 5c57abdc; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q -- "--perl-image-digest must be sha256:<64 hex>" "$_dd/out4.log"; then
        ok "a digest that is not sha256:<64 hex> is refused, so the check cannot be switched off"
    else
        bad "a malformed --perl-image-digest was accepted (rc=$rc)"
    fi
fi

echo
echo "CELL-LIST: --cells takes a comma-separated cell list; an unknown or repeated cell id is refused before any cell runs"

# --cells takes a comma-separated cell list (s01,s02,s07,s08,s09,s10), the comma form of --suites;
# the same stub docker and empty data tree as IMAGE-DIGEST, since the filter is applied before
# any cell's inputs are read.
if [[ -f "$HARNESS" ]]; then
    _cd="$TMP/cells"; mkdir -p "$_cd/data"
    run_cells_harness() { # run_cells_harness <out name> [harness args...]
        local o="$1"; shift
        : >"$_dd/calls"
        DIGEST_CALLS="$_dd/calls" STUB_DIGEST="$PINNED_DIGEST" PATH="$_dd/bin:$PATH" \
            bash "$HARNESS" --engine perl --data-dir "$_cd/data" --output-dir "$_cd/$o" "$@" >"$_cd/$o.log" 2>&1
    }
    run_cells_harness out1 --cells s07,s08,s09,s10
    if grep -q "run_clone_measurement engine=perl suites=s07 s08 s09 s10 walltime_only" "$_cd/out1.log" &&
        grep -q "suite s01 skipped (suites=s07 s08 s09 s10)" "$_cd/out1.log" && grep -q "suite s06 skipped" "$_cd/out1.log" &&
        grep -q "suite s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out1.log" && ! grep -q "suite s07 skipped" "$_cd/out1.log"; then
        ok "--cells s07,s08,s09,s10 selects those four cells and skips s01 to s06"
    else
        bad "--cells did not select the named cells: $(grep -E 'suites=|skipped|starting' "$_cd/out1.log" | head -4 | tr '\n' ' ')"
    fi
    run_cells_harness out2 --cells s07,s99; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: \[run_clone_measurement\] unknown cell 's99' in --suites/--cells; the cells are: s01 s02 s03 s04 s05 s06 s07 s08 s09 s10 s05m s06m" "$_cd/out2.log" &&
        ! grep -q "suite s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out2.log" && ! grep -q "^docker run" "$_dd/calls"; then
        ok "an unknown cell id stops the run before any cell starts and names the cells the table carries"
    else
        bad "an unknown cell id was not refused (rc=$rc): $(head -3 "$_cd/out2.log" | tr '\n' ' ')"
    fi
    run_cells_harness out3 --suites "s01 s0x"; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "unknown cell 's0x'" "$_cd/out3.log"; then
        ok "--suites with an unknown id is refused the same way (the id would have matched nothing and exited 0)"
    else
        bad "--suites with an unknown id was accepted (rc=$rc)"
    fi
    grep -q 'SUITES_FILTER="${2//,/ }"' "$HARNESS" && ! grep -qE '^cells[a-z_]*\(\) \{' "$HARNESS" &&
        ok "--cells is an inline rewrite of the comma list into the --suites filter, with no column-0 function of its own" ||
        bad "--cells is not the inline comma-to-space rewrite, or it has a column-0 function of its own"
    run_cells_harness out4 --cells ""; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: \[run_clone_measurement\] empty cell list in --suites/--cells" "$_cd/out4.log" &&
        ! grep -q "skipped" "$_cd/out4.log"; then
        ok "an empty cell list is refused (it would have skipped every cell and exited 0)"
    else
        bad "an empty cell list was accepted (rc=$rc): $(head -3 "$_cd/out4.log" | tr '\n' ' ')"
    fi
    # A glob-shaped id is validated as text: a file named like a cell in the working directory must not satisfy it.
    mkdir -p "$_cd/cwd" && : >"$_cd/cwd/s01"
    DIGEST_CALLS="$_dd/calls" STUB_DIGEST="$PINNED_DIGEST" PATH="$_dd/bin:$PATH" \
        bash -c 'cd "$1" && bash "$2" --engine perl --data-dir "$3" --output-dir "$4" --suites "s0*"' _ "$_cd/cwd" "$HARNESS" "$_cd/data" "$_cd/out5" >"$_cd/out5.log" 2>&1; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "unknown cell 's0\*'" "$_cd/out5.log"; then
        ok "a glob-shaped cell id is refused as text, not expanded against the working directory"
    else
        bad "a glob-shaped cell id was accepted (rc=$rc): $(head -3 "$_cd/out5.log" | tr '\n' ' ')"
    fi
    # The validation loop splits on any whitespace while should_run matches between single spaces:
    # a tab-separated list must run its cells, not validate and then skip every one with exit 0.
    run_cells_harness out6 --suites $'s07\ts08'
    if grep -q "run_clone_measurement engine=perl suites=s07 s08 walltime_only" "$_cd/out6.log" &&
        grep -q "suite s01 skipped (suites=s07 s08)" "$_cd/out6.log" &&
        grep -q "suite s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out6.log" && ! grep -q "suite s07 skipped" "$_cd/out6.log"; then
        ok "a tab-separated --suites list is normalised to one space per gap, so its cells run"
    else
        bad "a tab-separated --suites list validated but did not run its cells: $(grep -E 'suites=|skipped|starting' "$_cd/out6.log" | head -3 | tr '\n' ' ')"
    fi
    # A repeated id would pass the known-id check, and should_run's substring match would run
    # its cell once under a list that names it twice.
    run_cells_harness out7 --cells s07,s07,s08; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: \[run_clone_measurement\] repeated cell 's07' in --suites/--cells" "$_cd/out7.log" &&
        ! grep -q "suite s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out7.log" && ! grep -q "^docker run" "$_dd/calls"; then
        ok "a repeated cell id is refused before any cell starts, naming the id"
    else
        bad "a repeated cell id was accepted (rc=$rc): $(head -3 "$_cd/out7.log" | tr '\n' ' ')"
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
            GRCh37 "$TMP/sv_in" "$TMP/sv_scratch" "0:0" 5000 4 116 homo_sapiens
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
        _sp=$(grep -A 1 -x -- '--species' "$TMP/docker_args.txt" | tail -1)
        [[ "$_sp" == "homo_sapiens" ]] &&
            ok "vep received --species homo_sapiens, the value passed as the eleventh argument" ||
            bad "vep received --species '$_sp'; expected homo_sapiens"
    else
        bad "the stub docker was never invoked"
    fi
    _species_literal=$(grep -c -- '--species homo_sapiens' "$HARNESS") || _species_literal=0
    [[ "$_species_literal" -eq 0 ]] &&
        ok "no command line of the harness passes a literal --species homo_sapiens; the suite row's species travels as an argument" ||
        bad "$_species_literal command line(s) pass a literal --species homo_sapiens"
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
    for _t in bash dirname date mkdir tr; do ln -sf "$(type -P "$_t")" "$_np/$_t"; done
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
echo "FASTVEP-EXTRA-ARGS: --fastvep-extra-args reaches every fastVEP annotate run through the wrapper and is recorded in provenance.json; a run without the words writes neither; --sa-dir takes cells of one assembly"

# The words travel through <output-dir>/fastvep-wrapper.sh, which FASTVEP_BINARY is pointed at, so
# the SNP/indel call and the SV loop's `bash -c` subshell carry them with neither command line
# edited. The stub binary records its thread count and argument vector and writes one data line,
# enough for the SV output gate; a run is --walltime-only on one SNP/indel cell and one SV cell
# with two inputs, six invocations in all.
if [[ -f "$HARNESS" ]]; then
    grep -q '^FASTVEP_EXTRA_ARGS=""$' "$HARNESS" &&
        ok "--fastvep-extra-args defaults to nothing (FASTVEP_EXTRA_ARGS empty)" ||
        bad "FASTVEP_EXTRA_ARGS does not default to empty"
    grep -A 2 -e '^[[:space:]]*--fastvep-extra-args)$' "$HARNESS" | grep -q 'FASTVEP_EXTRA_ARGS="$2"' &&
        ok "--fastvep-extra-args is parsed into FASTVEP_EXTRA_ARGS" ||
        bad "--fastvep-extra-args is not parsed"
    _fextra_in_cmds=$(grep -cE 'FASTVEP_(EXTRA|WRAPPER|ARGV)' <(sed -n '/^measure_p1() {/,/^}/p; /^measure_sv() {/,/^}/p' "$HARNESS")) || _fextra_in_cmds=0
    [[ "$_fextra_in_cmds" -eq 0 ]] &&
        ok "neither measure_p1 nor measure_sv names the fastVEP extra arguments (they travel through FASTVEP_BINARY)" ||
        bad "measure_p1 or measure_sv names the fastVEP extra arguments on $_fextra_in_cmds line(s); the pass-through must not touch the command lines"
    _help=$(bash "$HARNESS" --help 2>&1)
    grep -q -e '^  --fastvep-extra-args "<flags>"$' <<<"$_help" &&
        ok "the usage text lists --fastvep-extra-args" ||
        bad "the usage text does not list --fastvep-extra-args"

    _fa="$TMP/fastvep_extra_args"
    _fa_sv="$_fa/data/ground_truth/perl/sv_per_vcf/canonical_inputs/grch37"
    mkdir -p "$_fa/bin" "$_fa/data/vcf/grch37/all_variants_canonical" "$_fa/data/vcf/grch38/all_variants_canonical" "$_fa_sv"
    : >"$_fa/data/vcf/grch37/all_variants_canonical/clinvar_grch37.canonical.vcf.gz"
    : >"$_fa/data/vcf/grch38/all_variants_canonical/clinvar_grch38.canonical.vcf.gz"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t1\n' >"$_fa_sv/01_a.vcf"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t2\n' >"$_fa_sv/02_b.vcf"
    cat >"$_fa/bin/fastvep" <<'STUB'
#!/usr/bin/env bash
# stub engine: one line per invocation, the thread count then the argument vector joined by spaces; one data line at -o
out=""; prev=""
for a in "$@"; do [[ "$prev" == -o ]] && out="$a"; prev="$a"; done
printf 'RAYON_NUM_THREADS=%s %s\n' "${RAYON_NUM_THREADS:-unset}" "$*" >>"$STUB_ARGV"
printf '#Uploaded_variation\tLocation\nv\t21:1\n' >"$out"
STUB
    chmod +x "$_fa/bin/fastvep"
    mkdir -p "$_fa/eq=bin"
    cp "$_fa/bin/fastvep" "$_fa/eq=bin/fastvep"
    # For the runs that name the binary without a slash: an executable copy of the stub in one
    # directory, a file of the same name that is not executable in another, and a directory
    # without one. Each is the working directory of one run or a PATH entry of another.
    mkdir -p "$_fa/bare-exec" "$_fa/bare-inert" "$_fa/bare-none"
    cp "$_fa/bin/fastvep" "$_fa/bare-exec/fastvep-bare"
    : >"$_fa/bare-inert/fastvep-bare"
    run_fastvep_extra() { # run_fastvep_extra <out name> [harness args...]; the cells are s01 and s07 unless --suites or --cells follows
        local o="$1"; shift
        : >"$_fa/$o.argv"
        STUB_ARGV="$_fa/$o.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" \
            --output-dir "$_fa/$o" --suites "s01 s07" --walltime-only "$@" >"$_fa/$o.log" 2>&1
    }
    # The wrapper, the log and the record a run left in its output directory, by name.
    fastvep_extra_files() {
        local f
        for f in fastvep-wrapper.argv fastvep-wrapper.sh provenance.json; do
            [[ ! -e "$_fa/$1/$f" ]] || printf '%s ' "$f"
        done
    }
    _fwords="--exclude-readthrough --exclude-artifact"
    if /usr/bin/time -v true >/dev/null 2>&1; then
        run_fastvep_extra with --fastvep-extra-args "$_fwords"; rc=$?
        _n_inv=$(grep -c . "$_fa/with.argv") || _n_inv=0
        _n_tail=$(grep -c -e "^RAYON_NUM_THREADS=16 annotate --gff3 .* --output-format tab $_fwords\$" "$_fa/with.argv") || _n_tail=0
        _n_p1=$(grep -c -e " -o $_fa/with/s01/rust/raw.txt --output-format tab $_fwords\$" "$_fa/with.argv") || _n_p1=0
        _n_sv=$(grep -c -e " -i $_fa_sv/0[12]_[ab].vcf -o $_fa/with/s07/rust/0[12]_[ab].txt --output-format tab $_fwords\$" "$_fa/with.argv") || _n_sv=0
        if [[ "$rc" -eq 0 && "$_n_inv" -eq 6 && "$_n_tail" -eq 6 ]]; then
            ok "all six fastVEP annotate runs end with the words after the harness's own arguments, at the harness's thread count"
        else
            bad "rc=$rc; $_n_inv invocation(s), $_n_tail of them ending with the words (expected 6 and 6): $(head -2 "$_fa/with.argv" | cut -c1-200 | tr '\n' '|') $(tail -3 "$_fa/with.log" | tr '\n' '|')"
        fi
        [[ "$_n_p1" -eq 2 && "$_n_sv" -eq 4 ]] &&
            ok "the SNP/indel cell's warmup and timed runs (2) and each input of the SV loop on both passes (4) ran through the wrapper" ||
            bad "$_n_p1 SNP/indel and $_n_sv SV invocations carried the words (expected 2 and 4)"
        _bad_rows=$(awk -F, 'NR>1 && ($1!="fastvep" || $10!="0" || ($8!="warmup" && $8!="timed"))' "$_fa/with/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        _rows=$(awk -F, 'NR>1' "$_fa/with/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$_rows" -eq 4 && "$_bad_rows" -eq 0 ]] &&
            ok "both cells record their warmup and timed rows with exit 0 through the fastVEP wrapper" ||
            bad "wall_times.csv has $_rows rows, $_bad_rows of them failed (expected 4 and 0)"
        [[ -x "$_fa/with/fastvep-wrapper.sh" ]] && grep -q "^exec $_fa/bin/fastvep \"\$@\" $_fwords\$" "$_fa/with/fastvep-wrapper.sh" &&
            ok "the fastVEP wrapper execs the real binary with the harness's arguments first and the words after" ||
            bad "fastvep-wrapper.sh is missing or does not exec the real binary with the words appended: $(tail -1 "$_fa/with/fastvep-wrapper.sh" 2>/dev/null)"
        if python3 - "$_fa/with/provenance.json" "$_fa/bin/fastvep" "$_fa/with/fastvep-wrapper.sh" <<'PY'
import json, sys
p, binary, wrapper = sys.argv[1:]
d = json.load(open(p))
words = ["--exclude-readthrough", "--exclude-artifact"]
assert sorted(d) == ["engine", "fastvep_binary", "fastvep_executed_argv", "fastvep_extra_args", "fastvep_extra_argv", "fastvep_wrapper"], sorted(d)
assert d["engine"] == "fastvep" and d["fastvep_binary"] == binary, d
assert d["fastvep_extra_args"] == " ".join(words) and d["fastvep_extra_argv"] == words, d
assert d["fastvep_wrapper"] == wrapper, d["fastvep_wrapper"]
assert len(d["fastvep_executed_argv"]) == 6, len(d["fastvep_executed_argv"])
for argv in d["fastvep_executed_argv"]:
    assert argv[0] == binary and argv[1] == "annotate" and argv[-4:] == ["--output-format", "tab"] + words, argv
assert sum(1 for a in d["fastvep_executed_argv"] if "/s07/rust/" in a[a.index("-o") + 1]) == 4, "SV invocations"
PY
        then
            ok "provenance.json names the engine, the real binary, the words joined and split, the wrapper and the six executed argument vectors, under fastVEP's own keys alone"
        else
            bad "provenance.json does not carry the expected record: $(tr -d '\n' <"$_fa/with/provenance.json" 2>/dev/null | cut -c1-300)"
        fi
        run_fastvep_extra absent; rc=$?
        _n_inv=$(grep -c . "$_fa/absent.argv") || _n_inv=0
        _n_plain=$(grep -c -e '^RAYON_NUM_THREADS=16 annotate --gff3 .* --output-format tab$' "$_fa/absent.argv") || _n_plain=0
        if [[ "$rc" -eq 0 && "$_n_inv" -eq 6 && "$_n_plain" -eq 6 && -z "$(fastvep_extra_files absent)" ]]; then
            ok "without the parameter the six runs end with the harness's own --output-format tab, and no wrapper, log or record is written"
        else
            bad "parameter absent: rc=$rc, $_n_inv invocations, $_n_plain ending with --output-format tab, left behind: $(fastvep_extra_files absent)"
        fi
        run_fastvep_extra empty --fastvep-extra-args ""; rc=$?
        run_fastvep_extra blank --fastvep-extra-args $' \n\t '; rc2=$?
        _n_plain=$(cat "$_fa/empty.argv" "$_fa/blank.argv" | grep -c -e ' --output-format tab$') || _n_plain=0
        if [[ "$rc" -eq 0 && "$rc2" -eq 0 && "$_n_plain" -eq 12 && -z "$(fastvep_extra_files empty)$(fastvep_extra_files blank)" ]]; then
            ok "an empty value and a value of whitespace alone append nothing and write no wrapper, log or record"
        else
            bad "empty value: rc=$rc/$rc2, $_n_plain of 12 invocations ending with --output-format tab, left behind: $(fastvep_extra_files empty)$(fastvep_extra_files blank)"
        fi
        _files_absent=$(cd "$_fa/absent" && find . | sort)
        [[ -n "$_files_absent" && "$_files_absent" == "$(cd "$_fa/empty" && find . | sort)" && "$_files_absent" == "$(cd "$_fa/blank" && find . | sort)" ]] &&
            ok "the empty and the whitespace value leave exactly the files a run without the parameter leaves" ||
            bad "the output directories of the absent, empty and whitespace runs hold different file lists"
        # "Split on whitespace" includes a newline, and the record's string is the words joined
        # by single spaces however the value separated them.
        run_fastvep_extra lines --fastvep-extra-args $'--exclude-readthrough\n   --exclude-artifact'; rc=$?
        _n_tail=$(grep -c -e " --output-format tab $_fwords\$" "$_fa/lines.argv") || _n_tail=0
        if [[ "$rc" -eq 0 && "$_n_tail" -eq 6 ]] && python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["fastvep_extra_args"] == sys.argv[2] and d["fastvep_extra_argv"] == sys.argv[2].split(), d' "$_fa/lines/provenance.json" "$_fwords" 2>/dev/null; then
            ok "a newline and repeated spaces inside the value separate words like one space: all six runs carry both, and the record joins them by single spaces"
        else
            bad "newline-separated value: rc=$rc, $_n_tail of 6 invocations carry both words: $(head -1 "$_fa/lines.argv" 2>/dev/null | cut -c1-200)"
        fi
        run_fastvep_extra values --fastvep-extra-args "--buffer-size 8 --distance=5000 -q a/b:c,d"; rc=$?
        _n_tail=$(grep -c -e ' --output-format tab --buffer-size 8 --distance=5000 -q a/b:c,d$' "$_fa/values.argv") || _n_tail=0
        [[ "$rc" -eq 0 && "$_n_tail" -eq 6 ]] &&
            ok "an option, a value after it, the --name=value form, a short option and a path-shaped value are admitted and passed in the order given" ||
            bad "option and value words: rc=$rc, $_n_tail of 6 invocations end with them in order: $(tail -2 "$_fa/values.log" | tr '\n' '|' | cut -c1-200)"
        # A --sa-dir directory over cells of one assembly is two more words (or one, in the
        # --name=value form), passed to every run and recorded like any other.
        _sa_words="$_fwords --sa-dir /data/sa/grch37"
        run_fastvep_extra sadir --fastvep-extra-args "$_sa_words"; rc=$?
        _n_tail=$(grep -c -e " --output-format tab $_sa_words\$" "$_fa/sadir.argv") || _n_tail=0
        if [[ "$rc" -eq 0 && "$_n_tail" -eq 6 ]] && python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["fastvep_extra_args"] == sys.argv[2] and d["fastvep_extra_argv"][-2:] == ["--sa-dir", "/data/sa/grch37"] and len(d["fastvep_executed_argv"]) == 6, d' "$_fa/sadir/provenance.json" "$_sa_words" 2>/dev/null; then
            ok "--sa-dir and its directory over two cells of one assembly reach all six runs and the record"
        else
            bad "--sa-dir over cells of one assembly: rc=$rc, $_n_tail of 6 invocations end with the words: $(tail -2 "$_fa/sadir.log" | tr '\n' '|' | cut -c1-240)"
        fi
        run_fastvep_extra sadireq --fastvep-extra-args "--sa-dir=/data/sa/grch37"; rc=$?
        _n_tail=$(grep -c -e ' --output-format tab --sa-dir=/data/sa/grch37$' "$_fa/sadireq.argv") || _n_tail=0
        [[ "$rc" -eq 0 && "$_n_tail" -eq 6 ]] &&
            ok "the --sa-dir=<directory> form over cells of one assembly is admitted the same way" ||
            bad "--sa-dir=<directory> over cells of one assembly: rc=$rc, $_n_tail of 6 invocations end with it: $(tail -2 "$_fa/sadireq.log" | tr '\n' '|' | cut -c1-240)"
        # The one-assembly rule binds only words that carry --sa-dir: the same cells of two
        # assemblies run under words without it.
        run_fastvep_extra twoasm --suites "s01 s02" --fastvep-extra-args "$_fwords"; rc=$?
        _n_37=$(grep -c -e "/vcf/grch37/.* --output-format tab $_fwords\$" "$_fa/twoasm.argv") || _n_37=0
        _n_38=$(grep -c -e "/vcf/grch38/.* --output-format tab $_fwords\$" "$_fa/twoasm.argv") || _n_38=0
        [[ "$rc" -eq 0 && "$_n_37" -eq 2 && "$_n_38" -eq 2 ]] &&
            ok "words without --sa-dir run cells of two assemblies, as a run without the parameter does" ||
            bad "words without --sa-dir over two assemblies: rc=$rc, $_n_37 GRCh37 and $_n_38 GRCh38 invocations (expected 2 and 2): $(tail -2 "$_fa/twoasm.log" | tr '\n' '|' | cut -c1-240)"
        # A run without the words into the directory of a run with them leaves no wrapper, log
        # or record behind to describe a wall_times.csv they did not produce.
        run_fastvep_extra rerun --fastvep-extra-args "$_fwords"; rc1=$?
        _left_before=$(fastvep_extra_files rerun)
        run_fastvep_extra rerun; rc2=$?
        if [[ "$rc1" -eq 0 && "$rc2" -eq 0 && "$_left_before" == "fastvep-wrapper.argv fastvep-wrapper.sh provenance.json " && -z "$(fastvep_extra_files rerun)" ]]; then
            ok "a re-run without the words removes the earlier run's wrapper, log and record from the output directory"
        else
            bad "re-run into the same output directory: rc=$rc1/$rc2, before '$_left_before', left behind '$(fastvep_extra_files rerun)'"
        fi
        # The refusal of an output directory whose path holds `=` binds a run with words alone:
        # without them no wrapper is written there, and the same directory runs.
        run_fastvep_extra eqplain --output-dir "$_fa/eq=plain/out" --suites s01; rc=$?
        _n_inv=$(grep -c -e " -o $_fa/eq=plain/out/s01/rust/raw.txt --output-format tab\$" "$_fa/eqplain.argv") || _n_inv=0
        _bad_rows=$(awk -F, 'NR>1 && $10!="0"' "$_fa/eq=plain/out/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$rc" -eq 0 && "$_n_inv" -eq 2 && "$_bad_rows" -eq 0 ]] &&
            ok "without the words an output directory whose path holds = runs its cell, both rows exit 0" ||
            bad "no words, = in the output directory: rc=$rc, $_n_inv of 2 invocations, $_bad_rows failed row(s): $(tail -2 "$_fa/eqplain.log" | tr '\n' '|' | cut -c1-240)"
        # With words the wrapper is the operand `env` reads and it execs the binary, so a binary
        # whose own path holds `=` runs.
        run_fastvep_extra eqbinwords --fastvep-binary "$_fa/eq=bin/fastvep" --suites s01 --fastvep-extra-args "$_fwords"; rc=$?
        _n_tail=$(grep -c -e " --output-format tab $_fwords\$" "$_fa/eqbinwords.argv") || _n_tail=0
        _bad_rows=$(awk -F, 'NR>1 && $10!="0"' "$_fa/eqbinwords/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$rc" -eq 0 && "$_n_tail" -eq 2 && "$_bad_rows" -eq 0 ]] &&
            ok "with the words a fastVEP binary whose path holds = runs through the wrapper, both rows exit 0" ||
            bad "words, = in the binary's path: rc=$rc, $_n_tail of 2 invocations, $_bad_rows failed row(s): $(tail -2 "$_fa/eqbinwords.log" | tr '\n' '|' | cut -c1-240)"
        run_fastvep_extra handover --fastvep-extra-args "$_fwords"; rc1=$?
        STUB_ARGV="$_fa/handover-veprs.argv" bash "$HARNESS" --engine vep-rs --vep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/handover" --suites s01 --walltime-only >"$_fa/handover-veprs.log" 2>&1; rc2=$?
        if [[ "$rc1" -eq 0 && "$rc2" -eq 0 && "$(fastvep_extra_files handover)" == "provenance.json " ]] && python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["engine"] == "vep-rs" and not [k for k in d if k.startswith("fastvep")], d' "$_fa/handover/provenance.json" 2>/dev/null; then
            ok "a vep-rs run into the directory of a fastVEP run with the words leaves the vep-rs record alone, no fastVEP wrapper or log"
        else
            bad "vep-rs run into a fastVEP output directory: rc=$rc1/$rc2, left behind '$(fastvep_extra_files handover)'"
        fi
        # The cells start a binary named without a slash through PATH, so the file PATH resolves
        # is the one tested: a file of that name in the working directory that is not executable
        # does not refuse the run.
        : >"$_fa/bare.argv"
        (cd "$_fa/bare-inert" && STUB_ARGV="$_fa/bare.argv" PATH="$_fa/bare-exec:$PATH" bash "$HARNESS" --engine fastvep --fastvep-binary fastvep-bare \
            --data-dir "$_fa/data" --output-dir "$_fa/bare" --suites s01 --walltime-only >"$_fa/bare.log" 2>&1); rc=$?
        _n_inv=$(grep -c -e " -o $_fa/bare/s01/rust/raw.txt --output-format tab\$" "$_fa/bare.argv") || _n_inv=0
        _bad_rows=$(awk -F, 'NR>1 && $10!="0"' "$_fa/bare/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$rc" -eq 0 && "$_n_inv" -eq 2 && "$_bad_rows" -eq 0 ]] &&
            ok "a slash-less --fastvep-binary is tested and started where PATH resolves it, whatever file of that name the working directory holds" ||
            bad "slash-less --fastvep-binary on PATH, a file of that name in the working directory: rc=$rc, $_n_inv of 2 invocations, $_bad_rows failed row(s): $(tail -2 "$_fa/bare.log" | tr '\n' '|' | cut -c1-240)"
    else
        echo "  skip functional fastVEP extra-args runs: /usr/bin/time -v is not GNU time here"
    fi
    # Each refusal comes before the output directory is made and before the run's opening line.
    for _eng in vep-rs perl; do
        bash "$HARNESS" --engine "$_eng" --vep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/refused-$_eng" --suites s01 --walltime-only --fastvep-extra-args "$_fwords" >"$_fa/refused-$_eng.log" 2>&1; rc=$?
        if [[ "$rc" -ne 0 && ! -e "$_fa/refused-$_eng" ]] && grep -q "^ERROR: --fastvep-extra-args applies to --engine fastvep only (got engine '$_eng')" "$_fa/refused-$_eng.log" && ! grep -q "run_clone_measurement engine=" "$_fa/refused-$_eng.log"; then
            ok "the words under --engine $_eng are refused by name before the run"
        else
            bad "the words under --engine $_eng were not refused before the run (rc=$rc): $(head -2 "$_fa/refused-$_eng.log" | tr '\n' ' ' | cut -c1-200)"
        fi
    done
    # value<TAB>the word of it the rule refuses: a shell metacharacter in an option and after a
    # value's first character, a letter outside ASCII (which the rule's ranges match under a
    # collating locale unless they are read in the C locale), a command substitution, a third
    # dash, a dash before a digit, a glob, and a quoted word after an admitted one. The harness
    # reads /dev/null, so no run can consume the rows the loop has yet to read.
    #
    # The rows run under a collating locale the host holds, the first that matches a letter
    # outside ASCII to [a-z], so the row with such a letter fails on a match made outside the C
    # locale whatever locale this file runs under. A host that holds none runs them under C.
    _word_locale=C
    while IFS= read -r _loc; do
        if LC_ALL="$_loc" bash -c '[[ é =~ ^[a-z]$ ]]' </dev/null 2>/dev/null; then
            _word_locale="$_loc"
            break
        fi
    done < <(echo en_US.UTF-8; locale -a 2>/dev/null)
    [[ "$_word_locale" != C ]] || echo "  skip the word rows' run under a collating locale: this host holds none that matches a letter outside ASCII to [a-z]"
    _wi=0
    while IFS=$'\t' read -r _value _word; do
        _wi=$((_wi + 1))
        : >"$_fa/word-$_wi.argv"
        LC_ALL="$_word_locale" STUB_ARGV="$_fa/word-$_wi.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/word-$_wi" --suites s01 --walltime-only --fastvep-extra-args "$_value" </dev/null >"$_fa/word-$_wi.log" 2>&1; rc=$?
        if [[ "$rc" -eq 1 && ! -e "$_fa/word-$_wi" && ! -s "$_fa/word-$_wi.argv" ]] && grep -qF -e "ERROR: [run_clone_measurement] --fastvep-extra-args word '$_word' is neither an option" "$_fa/word-$_wi.log" && ! grep -q "run_clone_measurement engine=" "$_fa/word-$_wi.log"; then
            ok "the word $_word is refused by name before any cell runs and before anything is written"
        else
            bad "the value '$_value' was not refused on the word '$_word' before the run (rc=$rc): $(head -2 "$_fa/word-$_wi.log" | tr '\n' ' ' | cut -c1-200)"
        fi
    done <<'WORDS'
--a;b	--a;b
a;b	a;b
café	café
$(id)	$(id)
---x	---x
-1	-1
a*	a*
--exclude-readthrough "--exclude-artifact"	"--exclude-artifact"
WORDS
    # On a host without a collating locale the rows above cannot tell a match made in the C locale
    # from one made in the caller's, so the match is read from the harness's source as well.
    grep -q -e '^[[:space:]]*(LC_ALL=C && \[\[ "\$word" =~ \$FASTVEP_EXTRA_OPTION || "\$word" =~ \$FASTVEP_EXTRA_VALUE \]\]) || {$' "$HARNESS" &&
        ok "the word rule is matched under LC_ALL=C, set inside the match's own subshell" ||
        bad "the word rule is not matched under LC_ALL=C; under a collating locale its ranges admit letters and digits outside ASCII"
    # With words the SNP/indel cells run the wrapper as the command operand of `env`, which reads
    # an operand holding `=` as an assignment and runs the next word instead; an output directory
    # whose path holds one is refused by name before it is made.
    : >"$_fa/eqdir.argv"
    STUB_ARGV="$_fa/eqdir.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/eq=refused/out" --suites s01 --walltime-only --fastvep-extra-args "$_fwords" >"$_fa/eqdir.log" 2>&1; rc=$?
    _n_err=$(grep -c '^ERROR: ' "$_fa/eqdir.log") || _n_err=0
    if [[ "$rc" -eq 1 && "$_n_err" -eq 1 && ! -e "$_fa/eq=refused" && ! -s "$_fa/eqdir.argv" ]] && grep -qF -e "ERROR: [run_clone_measurement] --output-dir '$_fa/eq=refused/out' holds '=': with --fastvep-extra-args" "$_fa/eqdir.log" && ! grep -q "run_clone_measurement engine=" "$_fa/eqdir.log"; then
        ok "with the words an output directory whose path holds = is refused by name with one ERROR line, before it is made"
    else
        bad "words and = in the output directory: rc=$rc, $_n_err ERROR line(s), directory made=$([[ -e "$_fa/eq=refused" ]] && echo yes || echo no): $(head -2 "$_fa/eqdir.log" | tr '\n' '|' | cut -c1-240)"
    fi
    # Without words the operand `env` reads is the binary's own path, so a binary whose path
    # holds `=` is refused by name the same way, before the output directory is made.
    : >"$_fa/eqbin.argv"
    STUB_ARGV="$_fa/eqbin.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/eq=bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/eqbin" --suites s01 --walltime-only >"$_fa/eqbin.log" 2>&1; rc=$?
    _n_err=$(grep -c '^ERROR: ' "$_fa/eqbin.log") || _n_err=0
    if [[ "$rc" -eq 1 && "$_n_err" -eq 1 && ! -e "$_fa/eqbin" && ! -s "$_fa/eqbin.argv" ]] && grep -qF -e "ERROR: [run_clone_measurement] --fastvep-binary '$_fa/eq=bin/fastvep' holds '=': " "$_fa/eqbin.log" && ! grep -q "run_clone_measurement engine=" "$_fa/eqbin.log"; then
        ok "without the words a fastVEP binary whose path holds = is refused by name with one ERROR line, before the output directory is made"
    else
        bad "no words and = in the binary's path: rc=$rc, $_n_err ERROR line(s), directory made=$([[ -e "$_fa/eqbin" ]] && echo yes || echo no), $(grep -c . "$_fa/eqbin.argv") invocation(s): $(head -2 "$_fa/eqbin.log" | tr '\n' '|' | cut -c1-240)"
    fi
    # A name without a slash that PATH resolves to no executable is refused by name before the
    # output directory is made. The working directory holds an executable of that name in both
    # runs, which the cells would never start.
    refused_bare_binary() { # refused_bare_binary <PATH entry under $_fa> <label>
        local on_path="$1" label="$2" n_err
        : >"$_fa/$on_path.argv"
        (cd "$_fa/bare-exec" && STUB_ARGV="$_fa/$on_path.argv" PATH="$_fa/$on_path:$PATH" bash "$HARNESS" --engine fastvep --fastvep-binary fastvep-bare \
            --data-dir "$_fa/data" --output-dir "$_fa/$on_path-out" --suites s01 --walltime-only >"$_fa/$on_path.log" 2>&1); rc=$?
        n_err=$(grep -c '^ERROR: ' "$_fa/$on_path.log") || n_err=0
        if [[ "$rc" -eq 1 && "$n_err" -eq 1 && ! -e "$_fa/$on_path-out" && ! -s "$_fa/$on_path.argv" ]] && grep -qF -e "ERROR: [run_clone_measurement] --fastvep-binary 'fastvep-bare' names no executable on PATH: " "$_fa/$on_path.log" && ! grep -q "run_clone_measurement engine=" "$_fa/$on_path.log"; then
            ok "a slash-less --fastvep-binary $label is refused by name with one ERROR line, before the output directory is made, an executable of that name in the working directory notwithstanding"
        else
            bad "slash-less --fastvep-binary $label: rc=$rc, $n_err ERROR line(s), directory made=$([[ -e "$_fa/$on_path-out" ]] && echo yes || echo no), $(grep -c . "$_fa/$on_path.argv") invocation(s): $(head -2 "$_fa/$on_path.log" | tr '\n' '|' | cut -c1-240)"
        fi
    }
    refused_bare_binary bare-none "that PATH does not resolve"
    refused_bare_binary bare-inert "whose only PATH match is not executable"
    # A --sa-dir directory holds one assembly's databases, so cells of two assemblies under it
    # stop the run with one ERROR line before any cell starts: a list of two cells, the comma
    # form with the --name=value word, and the default selection of every cell.
    refused_sa_dir() { # refused_sa_dir <out name> <label> <words> [cell selection...]
        local o="$1" label="$2" words="$3" n_err
        shift 3
        : >"$_fa/$o.argv"
        STUB_ARGV="$_fa/$o.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/$o" --walltime-only "$@" --fastvep-extra-args "$words" >"$_fa/$o.log" 2>&1; rc=$?
        n_err=$(grep -c '^ERROR: ' "$_fa/$o.log") || n_err=0
        if [[ "$rc" -eq 1 && "$n_err" -eq 1 && ! -s "$_fa/$o.argv" ]] && grep -q -e "^ERROR: \[run_clone_measurement\] --fastvep-extra-args carries --sa-dir, a directory of one assembly's databases, and the selected cells are of 2 assemblies: GRCh37 GRCh38; " "$_fa/$o.log" && ! grep -q ' starting$' "$_fa/$o.log"; then
            ok "--sa-dir over $label is refused with one ERROR line naming both assemblies, before any cell starts"
        else
            bad "--sa-dir over $label was not refused before any cell (rc=$rc, $n_err ERROR line(s), $(grep -c . "$_fa/$o.argv") invocation(s)): $(tail -2 "$_fa/$o.log" | tr '\n' '|' | cut -c1-300)"
        fi
    }
    refused_sa_dir sa-list "a list of one GRCh37 and one GRCh38 cell" "$_fwords --sa-dir /data/sa/grch38" --suites "s01 s02"
    refused_sa_dir sa-cells "the comma form naming both SV cells, in the --sa-dir=<directory> form" "--sa-dir=/data/sa/grch38" --cells s07,s08
    refused_sa_dir sa-all "the default selection of every cell" "--sa-dir /data/sa/grch38"
    # The record of a run with the words is written through python3, so a host without it is
    # refused at startup by the harness's own line. The PATH holds the tools the harness runs
    # before that check (tr normalises the cell list) and no interpreter.
    _fnp="$_fa/nopython/bin"
    mkdir -p "$_fnp"
    for _t in bash dirname date mkdir tr; do ln -sf "$(type -P "$_t")" "$_fnp/$_t"; done
    PATH="$_fnp" "$_fnp/bash" "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/nopython/out" --suites s01 --walltime-only --fastvep-extra-args "$_fwords" >"$_fa/nopython.log" 2>&1; rc=$?
    _n_err=$(grep -c '^ERROR: \[run_clone_measurement\] python3 is not on PATH: a fastVEP run with --fastvep-extra-args' "$_fa/nopython.log") || _n_err=0
    if [[ "$rc" -eq 1 && "$_n_err" -eq 1 && ! -e "$_fa/nopython/out" ]] && ! grep -q -e 'command not found' -e ' starting$' "$_fa/nopython.log"; then
        ok "without python3 a fastVEP run with the words stops before any cell with one line naming the interpreter"
    else
        bad "without python3: rc=$rc, $_n_err interpreter line(s): $(head -3 "$_fa/nopython.log" | tr '\n' '|')"
    fi
fi

echo
echo "REFUSED-RUN: a run refused on its arguments leaves the output directory as it found it"

# Each run below is refused for its cell selection, the last refusal the harness makes, twice:
# into a directory that holds an earlier run's wall_times.csv, record, wrapper and log, and
# into one that does not exist. The stub binary and the data tree are FASTVEP-EXTRA-ARGS's.
if [[ -f "$HARNESS" ]]; then
    _rr="$TMP/refused_run"
    mkdir -p "$_rr"
    refused_run_leaves() { # refused_run_leaves <out name> <label> <error text> [harness args...]
        local o="$1" label="$2" error="$3" f before rc_kept rc_absent
        shift 3
        mkdir -p "$_rr/$o"
        printf 'engine,run_id,host_id,replicate,suite_id,suite_name,assembly,run_kind,wall_seconds,exit_code\nfastvep,manual,local,local,s01,ClinVar_GRCh37,GRCh37,timed,1.5,0\n' >"$_rr/$o/wall_times.csv"
        for f in provenance.json fastvep-wrapper.sh fastvep-wrapper.argv; do
            echo "an earlier run's $f" >"$_rr/$o/$f"
        done
        before=$(cd "$_rr/$o" && cksum ./* 2>&1)
        bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_rr/$o" --walltime-only "$@" >"$_rr/$o.log" 2>&1; rc_kept=$?
        bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_rr/$o-absent" --walltime-only "$@" >"$_rr/$o-absent.log" 2>&1; rc_absent=$?
        if [[ "$rc_kept" -eq 1 && "$rc_absent" -eq 1 && "$(cd "$_rr/$o" && cksum ./* 2>&1)" == "$before" && ! -e "$_rr/$o-absent" ]] && grep -qF -e "ERROR: [run_clone_measurement] $error" "$_rr/$o.log"; then
            ok "$label is refused with the output directory untouched: an earlier run's wall_times.csv, record, wrapper and log keep their bytes, and a directory that does not exist is not made"
        else
            bad "$label: rc=$rc_kept/$rc_absent, the earlier run's files $([[ "$(cd "$_rr/$o" && cksum ./* 2>&1)" == "$before" ]] && echo kept || echo changed) (it holds $(cd "$_rr/$o" && echo ./*)), missing directory made=$([[ -e "$_rr/$o-absent" ]] && echo yes || echo no): $(head -2 "$_rr/$o.log" | tr '\n' '|' | cut -c1-200)"
        fi
    }
    refused_run_leaves sa-dir "--sa-dir over cells of two assemblies" "--fastvep-extra-args carries --sa-dir, a directory of one assembly's databases" --suites "s01 s02" --fastvep-extra-args "--sa-dir /data/sa/grch38"
    refused_run_leaves unknown "an unknown cell id" "unknown cell 's99' in --suites/--cells" --cells s07,s99
    refused_run_leaves empty "an empty cell list" "empty cell list in --suites/--cells" --cells ""
    refused_run_leaves repeated "a repeated cell id" "repeated cell 's07' in --suites/--cells" --cells s07,s07,s08
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
