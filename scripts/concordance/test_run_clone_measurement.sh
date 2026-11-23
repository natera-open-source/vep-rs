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
#   GROUND-TRUTH   a dataset with no Ensembl VEP ground truth that leaves a dangling
#                  symlink, which a copy of the output tree reports as an error.
#   IMAGE-DIGEST   an Ensembl VEP image whose pulled digest is not the pinned one, which
#                  must stop the run before any cell is timed under that
#                  engine_version.
#   CELL-LIST      a cell id the dataset table does not carry, which would match
#                  nothing, time nothing and exit 0, or a repeated one; --cells
#                  (a comma-separated list) and --datasets refuse both before any
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
#   FASTVEP-CACHE  a fastVEP cell started at a --transcript-cache path that names
#                  no file, where fastVEP parses the GFF3 beside it,
#                  saves a cache of its own at the path and exits 0, so the cell
#                  is timed and scored on the GFF3's gene models with nothing
#                  recorded.
#   TREE-COMPARATORS  a tuple or field comparison handed to anything but this tree's
#                  release builds (a Python comparator, a program an exported
#                  VEP_OUTPUTS_COMPARATOR or VEP_FIELDS_COMPARATOR names), or an
#                  engine started before a comparator the run would start is
#                  found not built.
#   TREE-NORMALISER  a fastVEP cell whose raw output is rewritten by anything but
#                  this tree's release build of tools/vep-fastvep-normalize (an
#                  interpreter, a program an exported FASTVEP_NORMALIZE names), an
#                  engine started before that build is found not built, or a raw
#                  output scored where the rewrite failed.
#   SV-COMPARATOR  a structural-variant cell scored by anything but the release build of
#                  tools/vep-sv-compare, a build under CARGO_TARGET_DIR passed over for
#                  the one under the tool's own target directory, or an engine started
#                  where neither holds a build and a structural-variant cell is to be
#                  scored.
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

echo "PIPES: no shell file of the tree pipes a writer into a quiet grep"
# Under pipefail a writer that has output left when a quiet grep leaves at its first match is
# ended by the closed pipe, and the pipeline's status is the writer's, whatever grep found: a
# check that must pass fails, a negated one passes and a script takes its other branch. A shell
# file that names pipefail reads a writer through a pipe with a grep that reads to the end.
PIPED_QUIET='(^|[^|])[|]&?[[:space:]]*([A-Za-z_][A-Za-z0-9_]*=[^[:space:]]*[[:space:]]+)*(command[[:space:]]+)?[ef]?grep[[:space:]]+(-[A-Za-z0-9]+[[:space:]]+)*(-[A-Za-z]*q|--quiet|--silent)'
_pipes_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
_piped=""
_piped_read=0
while IFS= read -r _pipes_file; do
    grep -q pipefail "$_pipes_file" || continue
    _piped_read=$((_piped_read + 1))
    _pipes_n=$(grep -E "$PIPED_QUIET" "$_pipes_file" | grep -cvE '^[[:space:]]*#')
    [[ $_pipes_n == 0 ]] || _piped+=" ${_pipes_file#"$_pipes_root"/}:$_pipes_n"
done < <(find "$_pipes_root/scripts" "$_pipes_root/bioconda" -type f -name '*.sh' 2>/dev/null | sort)
# The floor on the files read holds the walk: a walk that finds nothing reads no line of the form.
if [[ $_piped_read -ge 30 && -z $_piped ]]; then
    ok "no shell file that names pipefail pipes a writer into a quiet grep ($_piped_read files read)"
else
    bad "shell files that name pipefail pipe a writer into a quiet grep, a file and its count of lines ($_piped_read files read):$_piped"
fi

echo
echo "CACHE-TOOLS: the harness's comment on the vep-rs caches names the tool that converts a Storable cache"
# tools/vep-storable-to-json reads an Ensembl VEP Storable cache. vep-cache-converter reads a flat
# JSON directory or a plugin table and no Storable file. The comment is read as one line without
# its code-span marks, so the sentence is held however it is wrapped or marked.
_cache_tools_comment=$(awk '/^# vep-rs JSON caches\./ {on = 1} on && !/^#/ {on = 0} on {sub(/^#/, ""); for (i = 1; i <= NF; i++) printf "%s ", $i}' \
    "$(dirname "${BASH_SOURCE[0]}")/run_clone_measurement.sh")
_cache_tools_comment=${_cache_tools_comment//\`/}
if [[ $_cache_tools_comment == *'Storable cache with tools/vep-storable-to-json'* && $_cache_tools_comment != *'Storable cache with vep-cache-converter'* ]]; then
    ok "the comment gives tools/vep-storable-to-json as the converter of a Storable cache"
else
    bad "the comment on the vep-rs caches does not give tools/vep-storable-to-json as the converter of a Storable cache: '$_cache_tools_comment'"
fi

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

# The two gates above make measure_p1/measure_sv return non-zero. The main dataset
# loop runs under `set -euo pipefail`, so calling them bare makes the first bad
# cell kill the script, and a caller that ALSO runs under `set -e` inherits that
# abort, skipping every remaining step of the run. This asserts the guard shape
# in the shipped script rather than a copy of it: the two call sites must
# tolerate a non-zero return, and the run must still reach its completion line.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"

if [[ -f "$HARNESS" ]]; then
    # COMMENT LINES EXCLUDED. Without the exclusion the measure_p1 arm is unfalsifiable:
    # the harness's own explanatory comment reads
    #   `# site is guarded (\`measure_p1 ... || failed_datasets+=(...)\`). Exiting`
    # which the pattern matches, so deleting the REAL guard would still report ok.
    for _fn in measure_p1 measure_sv; do
        if grep -e "$_fn .*||" "$HARNESS" | grep -v '^[[:space:]]*#' >/dev/null; then
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
		failed_datasets=()
		cell() { [[ "$1" == "s02" ]] && return 1; return 0; }
		for s in s01 s02 s03; do
			echo -n "$s "
			cell "$s" || failed_datasets+=("$s")
		done
		echo "| failed=${#failed_datasets[@]}"
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
echo "GROUND-TRUTH: a dataset with no Perl ground truth must not leave a DANGLING symlink"

# The fault: the sample-scaling datasets s05m/s06m have no Ensembl VEP ground truth of their own
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
    # comparison branches in measure_sv are deliberately unguarded and create no symlink:
    # an SV dataset (s07 to s10) reads its assembly's per-VCF ground truth in place, or
    # records a pointer to the dataset of that assembly that did, and none of the four
    # is sample-scaling.
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

    # (a) a sample-scaling dataset: skipped, quietly, with no csv failure row.
    if ground_truth_ready s05m 1KG_multisample GRCh37 "$TMP/gt_absent" >/dev/null; then
        bad "ground_truth_ready accepted s05m, which has no ground truth by design"
    else
        [[ ! -s "$TMP/rows.txt" ]] &&
            ok "s05m is skipped without recording a failure row" ||
            bad "s05m recorded a failure row: $(cat "$TMP/rows.txt")"
    fi

    # (b) a headline dataset with the ground truth genuinely missing: loud, with a row.
    if ground_truth_ready s01 ClinVar_GRCh37 GRCh37 "$TMP/gt_absent" 2>/dev/null; then
        bad "ground_truth_ready accepted s01 with no ground truth staged"
    else
        grep -q 'ground_truth_missing' "$TMP/rows.txt" &&
            ok "an unstaged s01 fails the cell and records ground_truth_missing" ||
            bad "an unstaged s01 was skipped silently, which hides a staging failure"
    fi

    # (c) the happy path still passes.
    ground_truth_ready s01 ClinVar_GRCh37 GRCh37 "$TMP/gt_present" >/dev/null &&
        ok "a staged dataset still reaches the comparator" ||
        bad "ground_truth_ready rejected a dataset whose ground truth exists"

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

# The fault: a vep-rs SV branch that hardcodes `--fork 1` while the Ensembl VEP and fastVEP SV
# branches pass `$FORK` (16). On the four SV cells (s07 to s10, annotated once an assembly)
# vep-rs is then the only engine of the three annotating serially, and every SV ratio is a
# floor on its advantage rather than a measurement of it. The 1 is not needed to match the
# per-VCF SV ground truth: vep-rs output is fork-invariant (FORK-INVARIANCE below) and the
# ground truth is a fixed archived artifact.
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
    # Line continuations are collapsed first: the Ensembl VEP branch's argument list wraps across
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

    # Both Ensembl VEP branches must call it. A guard defined and never invoked is the failure mode
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
    _rm_lines=$(grep -n -e 'rm -rf "\$dataset_dir/perl_scratch"' "$HARNESS" | cut -d: -f1)

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
    # Ensembl VEP branch's timed work. A timing loop must sit BETWEEN them, and a scratch delete
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
    # source, emits one warning per record and exits 0, so an Ensembl VEP wall time measured with
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
    _first_time_run=$(grep -n -e 'res=\$(time_run "\${dataset_id}_\${run_kind}"' "$HARNESS" | cut -d: -f1 | awk -v d="${_sv_def:-0}" '$1 > d' | head -1)
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
        /^[[:space:]]*"\$[{]FIELDS_COMPARE\[@\][}]" \\$/ { if (cond ~ /ground_truth_ready/) n++ }
        END { print n + 0 }
    ' "$HARNESS")
    [[ "$_fields_in_guard" -eq 1 ]] &&
        ok "the field comparator runs once, inside the guarded vep-rs SNP/indel branch" ||
        bad "expected the field comparator once inside a ground_truth_ready branch, found $_fields_in_guard"
    grep -q -- '--out "$dataset_dir/report/fields.json"' "$HARNESS" &&
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
    grep -A 2 -e '^[[:space:]]*--skip-fields)$' "$HARNESS" | grep 'SKIP_FIELDS="true"' >/dev/null &&
        ok "--skip-fields is parsed and sets SKIP_FIELDS=true" ||
        bad "--skip-fields is not parsed into SKIP_FIELDS=true"
    _fields_skip_guarded=$(awk '
        /if \[\[ "\$SKIP_FIELDS" != "true" \]\]/ { g = 1 }
        /^[[:space:]]*"\$[{]FIELDS_COMPARE\[@\][}]" \\$/ { if (g) n++; g = 0 }
        END { print n + 0 }
    ' "$HARNESS")
    [[ "$_fields_skip_guarded" -eq 1 ]] &&
        ok "the field comparator runs only while SKIP_FIELDS is not true" ||
        bad "expected the field comparator once under the SKIP_FIELDS guard, found $_fields_skip_guarded"
    _skip_line=$(grep -n -e 'if \[\[ "\$SKIP_FIELDS" != "true" \]\]' "$HARNESS" | head -1 | cut -d: -f1)
    # The F1 compare runs the resolved tuple comparator command, `"${OUTPUTS_COMPARE[@]}"`.
    _f1_line=$(grep -n -e '"\${OUTPUTS_COMPARE\[@\]}"' "$HARNESS" | head -1 | cut -d: -f1)
    [[ -n "$_skip_line" && -n "$_f1_line" && "$_f1_line" -lt "$_skip_line" ]] &&
        ok "the F1 compare is outside the SKIP_FIELDS guard (line $_f1_line before $_skip_line)" ||
        bad "the F1 compare must not depend on SKIP_FIELDS (F1 at ${_f1_line:-?}, guard at ${_skip_line:-?})"
fi

echo "TREE-COMPARATORS: a cell that scores starts this tree's release builds of the two comparators and nothing else"

# The field comparator call of a vep-rs cell starts the command the run resolved,
# FIELDS_COMPARE, and nothing else: the SKIP_FIELDS block is extracted and run with a stub at
# the tree's field comparator path, a python3 stub and a stub a caller's variable could name,
# each recording its start. The tree's binary runs with --perl, --rust, --out and --tmp-dir;
# an exported VEP_FIELDS_COMPARATOR or VEP_OUTPUTS_COMPARATOR changes nothing, whatever it
# names; and a comparison that fails is recorded with nothing run in its place.
if [[ -f "$HARNESS" ]]; then
    _fields_block=$(awk '
        !p && /^[[:space:]]*if \[\[ "\$SKIP_FIELDS" != "true" \]\]; then$/ { match($0, /^[[:space:]]*/); ind = substr($0, 1, RLENGTH); p = 1 }
        p { print }
        p && $0 == ind "fi" { exit }
    ' "$HARNESS")
    eval "_fields_block_fn() {
$_fields_block
}"
    _fdir="$TMP/tree_comparators"
    _outputs_rel="tools/vep-outputs-compare/target/release/vep-outputs-compare"
    _fields_rel="tools/vep-fields-compare/target/release/vep-fields-compare"
    _normalise_rel="tools/vep-fastvep-normalize/target/release/vep-fastvep-normalize"
    _started="$_fdir/started"
    mkdir -p "$_fdir/bin" "$_fdir/dataset/report"
    # Each stub appends its own path and arguments, one line a start; `failing` then exits 3.
    for _stub in python3 comparator failing; do
        printf '#!/bin/bash\nprintf "%%s\\n" "$0 $*" >>"%s"\n' "$_started" >"$_fdir/bin/$_stub"
        chmod +x "$_fdir/bin/$_stub"
    done
    echo "exit 3" >>"$_fdir/bin/failing"
    : >"$_fdir/not_executable"
    # bash's PATH search returns a non-executable match when no executable one exists.
    : >"$_fdir/bin/inert"
    # _tree <name> <tuple> <field>: a repository root under $_fdir whose two release paths
    # are each `built` (a stub), `absent`, `inert` (a file that is not executable) or `dir`.
    _tree() {
        local root="$_fdir/$1" slot rel
        mkdir -p "$root"
        for slot in "$2 $_outputs_rel" "$3 $_fields_rel"; do
            rel="${slot#* }"
            case "${slot%% *}" in
            built) mkdir -p "$root/$(dirname "$rel")" && cp "$_fdir/bin/comparator" "$root/$rel" ;;
            inert) mkdir -p "$root/$(dirname "$rel")" && : >"$root/$rel" ;;
            dir) mkdir -p "$root/$rel" ;;
            esac
        done
    }
    _tree repo built built
    _tree repo-unbuilt absent absent
    _tree repo-tuple-only built absent
    _tree repo-fields-only absent built
    _tree repo-tuple-inert inert built
    _tree repo-tuple-dir dir built
    _tree repo-fields-inert built inert
    _tree repo-fields-dir built dir
    # _normaliser <tree> <state>: the tree's release path of the fastVEP normaliser, `built`
    # (a stub), `failing` (a stub that exits 3), `inert` (a file that is not executable) or
    # `dir`. The trees above hold none.
    _normaliser() {
        local path="$_fdir/$1/$_normalise_rel"
        mkdir -p "$(dirname "$path")"
        case "$2" in
        built) cp "$_fdir/bin/comparator" "$path" ;;
        failing) cp "$_fdir/bin/failing" "$path" ;;
        inert) : >"$path" ;;
        dir) mkdir -p "$path" ;;
        esac
    }
    for _state in built failing inert dir; do
        _tree "repo-normaliser-$_state" built absent
        _normaliser "repo-normaliser-$_state" "$_state"
    done
    _tree repo-all built built
    _normaliser repo-all built
    # What a caller may have exported, and what each form is: "" leaves both variables
    # unset, `empty` exports them empty, any other word is exported as the value of both.
    _exported=(empty python "$_fdir/bin/comparator" comparator "$_fdir/not_executable" "$_fdir" inert)
    _exported_as=("the empty string" "the word python" "an executable file" "a bare name on PATH" "a file that is not executable" "a directory" "a bare name whose PATH match is not executable")
    _export_comparators() {
        case "$1" in
        "") unset VEP_OUTPUTS_COMPARATOR VEP_FIELDS_COMPARATOR ;;
        empty) export VEP_OUTPUTS_COMPARATOR="" VEP_FIELDS_COMPARATOR="" ;;
        *) export VEP_OUTPUTS_COMPARATOR="$1" VEP_FIELDS_COMPARATOR="$1" ;;
        esac
    }
    # $1 is what the caller exported, $2 the resolved field comparator (the tree's by
    # default); clears the record of starts, then prints the block's stdout and stderr.
    _run_fields_block() {
        : >"$_started"
        # The extracted block reads the names this subshell sets.
        # shellcheck disable=SC2034,SC2329
        (
            PATH="$_fdir/bin:$PATH"
            dataset_dir="$_fdir/dataset" dataset_id=s06 dataset_name=stub assembly=GRCh38 SKIP_FIELDS=false
            FIELDS_COMPARE=("${2:-$_fdir/repo/$_fields_rel}")
            csv_row() { echo "csv_row $*"; }
            _export_comparators "$1"
            _fields_block_fn
        ) 2>&1
    }
    _fields_flags="--perl $_fdir/dataset/perl/output.txt --rust $_fdir/dataset/rust/output.txt --out $_fdir/dataset/report/fields.json --tmp-dir $_fdir/dataset"
    _out=$(_run_fields_block "")
    [[ "$(grep -c . "$_started")" -eq 1 && "$(cut -d' ' -f1 "$_started")" == "$_fdir/repo/$_fields_rel" && -z "$_out" ]] &&
        ok "the field comparator call starts the tree's release build of tools/vep-fields-compare, once and with no warning" ||
        bad "field comparator call: expected one start of the tree's binary and no output, got '$(cat "$_started")' output '$_out'"
    [[ "$(cut -d' ' -f2- "$_started")" == "$_fields_flags" ]] &&
        ok "the field comparator receives --perl, --rust, --out report/fields.json and --tmp-dir" ||
        bad "field comparator flags: expected '$_fields_flags', got '$(cut -d' ' -f2- "$_started")'"
    for _i in "${!_exported[@]}"; do
        _out=$(_run_fields_block "${_exported[$_i]}")
        [[ "$(cat "$_started")" == "$_fdir/repo/$_fields_rel $_fields_flags" && -z "$_out" ]] &&
            ok "field comparator call, both comparator variables exported as ${_exported_as[$_i]}: the tree's binary runs with no warning, and nothing else starts" ||
            bad "field comparator call with ${_exported[$_i]} exported: expected the tree's binary alone and no output, got '$(cat "$_started")' output '$_out'"
    done
    _out=$(_run_fields_block python "$_fdir/bin/failing")
    [[ "$(cat "$_started")" == "$_fdir/bin/failing $_fields_flags" && "$_out" == "WARNING: [run_clone_measurement] s06 field comparison failed; report/fields.json is absent"$'\n'"csv_row s06 stub GRCh38 fields_compare_failed  1" ]] &&
        ok "a field comparison that fails is warned about and recorded as fields_compare_failed, and nothing runs in its place" ||
        bad "failing field comparator: expected its one start, the warning and the CSV row, got '$(cat "$_started")' output '$_out'"
fi

# The two comparators are this tree's release builds and nothing else: the script-level block
# that resolves them (OUTPUTS_COMPARE, FIELDS_COMPARE) is extracted and run in trees whose
# release paths hold a stub, nothing, a file that is not executable or a directory. A run
# that scores resolves the tree's two builds, whatever VEP_OUTPUTS_COMPARATOR or
# VEP_FIELDS_COMPARATOR a caller exported; a comparator the run would start that is not an
# executable file stops it with one line naming the build command, and nothing is taken in
# its place; a run that scores no F1 (--walltime-only, the Ensembl VEP engine) asks for neither and
# one that compares no fields (--skip-fields, the fastVEP engine) for the tuple comparator
# alone; the block stands before the first engine run and the call sites run what it resolved.
if [[ -f "$HARNESS" ]]; then
    _comparators_block=$(awk '
        !p && /^OUTPUTS_COMPARE=\(\)$/ { p = 1 }
        p { print }
        p && /^fi$/ { exit }
    ' "$HARNESS")
    [[ "$_comparators_block" == *"outputs_bin=\"\$REPO_DIR/$_outputs_rel\""* && "$_comparators_block" == *"fields_bin=\"\$REPO_DIR/$_fields_rel\""* && "$_comparators_block" == *"normalize_bin=\"\$REPO_DIR/$_normalise_rel\""* && "$_comparators_block" == *$'\nfi' ]] &&
        ok "the block that resolves the two comparators and the fastVEP normaliser is extracted whole, the tree's release builds its only candidates" ||
        bad "could not extract the block that resolves the two comparators and the fastVEP normaliser from the harness"
    # $1 is the repository root, $2 what the caller exported, $3 the engine, $4 the
    # WALLTIME_ONLY value and $5 the SKIP_FIELDS value; prints the block's stderr, then
    # `tuple: <command>`, `fields: <command>` and `normaliser: <command>` for each command
    # resolved, then `rc=<exit status>`.
    _run_comparators_block() {
        # The extracted block reads the names this subshell sets.
        # shellcheck disable=SC2034
        (
            PATH="$_fdir/bin:$PATH"
            REPO_DIR="$_fdir/$1"
            ENGINE="${3:-vep-rs}"
            WALLTIME_ONLY="${4:-false}"
            SKIP_FIELDS="${5:-false}"
            _export_comparators "$2"
            eval "$_comparators_block"
            [[ ${#OUTPUTS_COMPARE[@]} -eq 0 ]] || echo "tuple: ${OUTPUTS_COMPARE[*]}"
            [[ ${#FIELDS_COMPARE[@]} -eq 0 ]] || echo "fields: ${FIELDS_COMPARE[*]}"
            [[ ${#FASTVEP_NORMALIZE[@]} -eq 0 ]] || echo "normaliser: ${FASTVEP_NORMALIZE[*]}"
        ) 2>&1
        echo "rc=$?"
    }
    # The one line a comparator that cannot run is refused with: $1 its role (tuple or
    # field), $2 its tool, $3 the tree.
    _refusal() {
        printf '%s' "ERROR: [run_clone_measurement] the $1 comparator $_fdir/$3/tools/$2/target/release/$2 is not an executable: build the tool with 'cargo build --release --locked --manifest-path $_fdir/$3/tools/$2/Cargo.toml'"
    }
    _both="tuple: $_fdir/repo/$_outputs_rel"$'\n'"fields: $_fdir/repo/$_fields_rel"$'\n'"rc=0"
    _out=$(_run_comparators_block repo "")
    [[ "$_out" == "$_both" ]] &&
        ok "a vep-rs run that scores: the tree's release builds of tools/vep-outputs-compare and tools/vep-fields-compare are the two commands" ||
        bad "both tools built: expected the tree's two builds, got '$_out'"
    for _i in "${!_exported[@]}"; do
        _out=$(_run_comparators_block repo "${_exported[$_i]}")
        [[ "$_out" == "$_both" ]] &&
            ok "both comparator variables exported as ${_exported_as[$_i]}: the same two commands, with no warning and no refusal" ||
            bad "${_exported[$_i]} exported: expected the tree's two builds alone, got '$_out'"
    done
    _out=$(_run_comparators_block repo-unbuilt "")
    [[ "$_out" == "$(_refusal tuple vep-outputs-compare repo-unbuilt)"$'\n'"rc=1" ]] &&
        ok "neither tool built: one line names the tuple comparator's build command, exit 1, no command is resolved" ||
        bad "neither tool built: expected the tuple comparator's refusal alone and rc=1, got '$_out'"
    _out=$(_run_comparators_block repo-fields-only "")
    [[ "$_out" == "$(_refusal tuple vep-outputs-compare repo-fields-only)"$'\n'"rc=1" ]] &&
        ok "the tuple comparator not built beside a built field comparator: the tuple comparator's refusal, exit 1" ||
        bad "tuple comparator not built: expected its refusal alone and rc=1, got '$_out'"
    _out=$(_run_comparators_block repo-tuple-only "")
    [[ "$_out" == "$(_refusal field vep-fields-compare repo-tuple-only)"$'\n'"rc=1" ]] &&
        ok "the field comparator not built: one line names its build command, exit 1, and no interpreter is taken in its place" ||
        bad "field comparator not built: expected its refusal alone and rc=1, got '$_out'"
    for _i in 1 2 3; do
        _out=$(_run_comparators_block repo-unbuilt "${_exported[$_i]}")
        [[ "$_out" == "$(_refusal tuple vep-outputs-compare repo-unbuilt)"$'\n'"rc=1" ]] &&
            ok "neither tool built and both comparator variables exported as ${_exported_as[$_i]}: the same refusal, and what the variables name is not taken" ||
            bad "neither tool built with ${_exported[$_i]} exported: expected the tuple comparator's refusal alone and rc=1, got '$_out'"
        _out=$(_run_comparators_block repo-tuple-only "${_exported[$_i]}")
        [[ "$_out" == "$(_refusal field vep-fields-compare repo-tuple-only)"$'\n'"rc=1" ]] &&
            ok "the field comparator not built and both comparator variables exported as ${_exported_as[$_i]}: the same refusal, and what the variables name is not taken" ||
            bad "field comparator not built with ${_exported[$_i]} exported: expected its refusal alone and rc=1, got '$_out'"
    done
    for _shape in inert dir; do
        _out=$(_run_comparators_block "repo-tuple-$_shape" "")
        [[ "$_out" == "$(_refusal tuple vep-outputs-compare "repo-tuple-$_shape")"$'\n'"rc=1" ]] &&
            ok "the tuple comparator's release path $([[ "$_shape" == dir ]] && echo "a directory" || echo "a file that is not executable"): refused with its build command, exit 1" ||
            bad "tuple comparator path $_shape: expected its refusal alone and rc=1, got '$_out'"
        _out=$(_run_comparators_block "repo-fields-$_shape" "")
        [[ "$_out" == "$(_refusal field vep-fields-compare "repo-fields-$_shape")"$'\n'"rc=1" ]] &&
            ok "the field comparator's release path $([[ "$_shape" == dir ]] && echo "a directory" || echo "a file that is not executable"): refused with its build command, exit 1" ||
            bad "field comparator path $_shape: expected its refusal alone and rc=1, got '$_out'"
    done
    _out=$(_run_comparators_block repo-unbuilt "" vep-rs true)
    [[ "$_out" == "rc=0" ]] &&
        ok "--walltime-only: no comparator is asked for, so a tree without the tools does not stop the run" ||
        bad "--walltime-only in a tree without the tools: expected no command and rc=0, got '$_out'"
    _out=$(_run_comparators_block repo-unbuilt "" perl false)
    [[ "$_out" == "rc=0" ]] &&
        ok "the Perl engine scores no F1: no comparator is asked for" ||
        bad "the Perl engine in a tree without the tools: expected no command and rc=0, got '$_out'"
    _out=$(_run_comparators_block repo-unbuilt "" fastvep false)
    [[ "$_out" == "$(_refusal tuple vep-outputs-compare repo-unbuilt)"$'\n'"rc=1" ]] &&
        ok "the fastVEP engine scores F1: a tuple comparator that is not built stops it as it stops vep-rs" ||
        bad "the fastVEP engine in a tree without the tools: expected the tuple comparator's refusal and rc=1, got '$_out'"
    _out=$(_run_comparators_block repo-normaliser-built "" fastvep false)
    [[ "$_out" == "tuple: $_fdir/repo-normaliser-built/$_outputs_rel"$'\n'"normaliser: $_fdir/repo-normaliser-built/$_normalise_rel"$'\n'"rc=0" ]] &&
        ok "the fastVEP engine compares no fields: the tuple comparator and the normaliser are resolved, and a field comparator that is not built does not stop the run" ||
        bad "the fastVEP engine without the field comparator: expected the tuple comparator and the normaliser and rc=0, got '$_out'"
    _out=$(_run_comparators_block repo-tuple-only "" vep-rs false true)
    [[ "$_out" == "tuple: $_fdir/repo-tuple-only/$_outputs_rel"$'\n'"rc=0" ]] &&
        ok "--skip-fields: the field comparator is not asked for, so a tree without it still scores F1" ||
        bad "--skip-fields without the field comparator: expected the tuple comparator alone and rc=0, got '$_out'"
    _out=$(_run_comparators_block repo-unbuilt "" vep-rs false true)
    [[ "$_out" == "$(_refusal tuple vep-outputs-compare repo-unbuilt)"$'\n'"rc=1" ]] &&
        ok "--skip-fields still scores F1: a tuple comparator that is not built stops the run" ||
        bad "--skip-fields in a tree without the tools: expected the tuple comparator's refusal and rc=1, got '$_out'"
    _block_line=$(grep -n -e '^OUTPUTS_COMPARE=()$' "$HARNESS" | head -1 | cut -d: -f1)
    _cells_line=$(grep -n -e '^for entry in "\${DATASETS\[@\]}"; do$' "$HARNESS" | head -1 | cut -d: -f1)
    [[ -n "$_block_line" && -n "$_cells_line" && "$_block_line" -lt "$_cells_line" ]] &&
        ok "the comparators are resolved before the loop that runs the cells (line $_block_line before $_cells_line)" ||
        bad "the block that resolves the comparators must stand before the cell loop (block at ${_block_line:-?}, loop at ${_cells_line:-?})"
    _outputs_calls=$(grep -c -e '^[[:space:]]*"\${OUTPUTS_COMPARE\[@\]}" \\$' "$HARNESS") || _outputs_calls=0
    _fields_calls=$(grep -c -e '^[[:space:]]*"\${FIELDS_COMPARE\[@\]}" \\$' "$HARNESS") || _fields_calls=0
    [[ "$_outputs_calls" -eq 2 && "$_fields_calls" -eq 1 ]] &&
        ok "both F1 call sites (vep-rs, fastvep) run the resolved tuple comparator and the one field call site the resolved field comparator" ||
        bad "expected 2 OUTPUTS_COMPARE call sites and 1 FIELDS_COMPARE call site, found $_outputs_calls and $_fields_calls"
    _named=$(grep -c -E -e 'compare_vep_(outputs|fields)\.py|VEP_(OUTPUTS|FIELDS)_COMPARATOR' "$HARNESS") || _named=0
    [[ "$_named" -eq 0 ]] &&
        ok "no line of the harness names a Python tuple or field comparator, VEP_OUTPUTS_COMPARATOR or VEP_FIELDS_COMPARATOR" ||
        bad "$_named line(s) of the harness name a Python tuple or field comparator or a comparator variable"

    # The whole harness, not the blocks: a copy of it in each tree runs one vep-rs cell, s01,
    # over a three-record input with a stub engine. The engine's output differs from the
    # ground truth in one row's Consequence and in another row's Existing_variation.
    _od="$TMP/tree_comparators_whole"
    _od_vcf="$_od/data/vcf/grch37/all_variants_canonical/clinvar_grch37.canonical.vcf.gz"
    _od_gt="$_od/data/ground_truth/perl/snp_indel/s01"
    mkdir -p "$_od/bin" "$(dirname "$_od_vcf")" "$_od_gt" "$_od/data/caches/vep-rs/grch37/transcripts"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1000\tv1\tA\tG\t.\t.\t.\n21\t2000\tv2\tC\tT\t.\t.\t.\n21\t3000\tv3\tG\tA\t.\t.\t.\n' | gzip -c >"$_od_vcf"
    _od_header='#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n'
    _od_row1='v1\t21:1000\tG\tENSG01\tENST01\tTranscript\tmissense_variant\t10\t10\t4\tK/R\taAg/aGg\t%s\tIMPACT=MODERATE\n'
    _od_row2='v2\t21:2000\tT\tENSG01\tENST01\tTranscript\t%s\t20\t20\t7\tL\tctC/ctT\t-\tIMPACT=LOW\n'
    _od_row3='v3\t21:3000\tA\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n'
    # The row formats carry the one field that differs between the two outputs.
    # shellcheck disable=SC2059
    {
        printf "$_od_header" && printf "$_od_row1" rs1 && printf "$_od_row2" synonymous_variant && printf "$_od_row3"
    } >"$_od_gt/output.txt"
    # shellcheck disable=SC2059
    {
        printf "$_od_header" && printf "$_od_row1" - && printf "$_od_row2" stop_gained && printf "$_od_row3"
    } >"$_od/engine_output.txt"
    cat >"$_od/bin/vep" <<STUB
#!/usr/bin/env bash
# stub engine: marks its start and writes the prepared output at -o
out=""; prev=""
for a in "\$@"; do [[ "\$prev" == -o ]] && out="\$a"; prev="\$a"; done
: >"$_od/engine.ran"
cp "$_od/engine_output.txt" "\$out"
STUB
    chmod +x "$_od/bin/vep"
    # _run_in_tree <tree> [harness args...]: the copy of the harness in the tree runs s01
    # into a fresh $_od/out; the log is $_od/log and the stubs' starts are in $_started.
    _run_in_tree() {
        local tree="$_fdir/$1"
        shift
        mkdir -p "$tree/scripts/concordance"
        cp "$HARNESS" "$(dirname "$HARNESS")/harness_genomes.sh" "$(dirname "$HARNESS")/harness_genomes.tsv" "$tree/scripts/concordance/"
        rm -rf "$_od/out" "$_od/engine.ran"
        : >"$_started"
        bash "$tree/scripts/concordance/run_clone_measurement.sh" --engine vep-rs --vep-binary "$_od/bin/vep" --data-dir "$_od/data" \
            --output-dir "$_od/out" --datasets s01 --reference-release 116.2 "$@" >"$_od/log" 2>&1
    }
    _tuple_start="$_outputs_rel --perl-dir $_od/out/s01/perl --rust-dir $_od/out/s01/rust --report-dir $_od/out/s01/report --reference-release 116.2 --input-vcf $_od_vcf --assembly GRCh37"
    # The vep-rs cell's call ends with the cache the cell annotated with; the fastVEP cell's does not.
    _veprs_tuple_start="$_tuple_start --json-cache $_od/data/caches/vep-rs/grch37"
    _fields_start="$_fields_rel --perl $_od/out/s01/perl/output.txt --rust $_od/out/s01/rust/output.txt --out $_od/out/s01/report/fields.json --tmp-dir $_od/out/s01"
    _run_in_tree repo-unbuilt
    rc=$?
    [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_od/log")" -eq 1 && "$(tail -1 "$_od/log")" == "$(_refusal tuple vep-outputs-compare repo-unbuilt)" && ! -e "$_od/engine.ran" && ! -e "$_od/out" ]] &&
        ok "a harness run in a tree without the tools: exit 1 on the tuple comparator's refusal line, no output directory made, the engine never started" ||
        bad "harness without the tools: expected rc=1, the refusal as the one ERROR and last line and no engine run, got rc=$rc '$(tail -2 "$_od/log")'"
    _run_in_tree repo-tuple-only
    rc=$?
    [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_od/log")" -eq 1 && "$(tail -1 "$_od/log")" == "$(_refusal field vep-fields-compare repo-tuple-only)" && ! -e "$_od/engine.ran" && ! -e "$_od/out" && ! -s "$_started" ]] &&
        ok "a harness run in a tree without the field comparator: exit 1 on its refusal line, no output directory made, neither the engine nor the tuple comparator started" ||
        bad "harness without the field comparator: expected rc=1, the refusal as the one ERROR and last line and nothing started, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
    (
        export VEP_OUTPUTS_COMPARATOR="$_fdir/bin/comparator" VEP_FIELDS_COMPARATOR="$_fdir/bin/comparator"
        _run_in_tree repo-unbuilt
    )
    rc=$?
    [[ "$rc" -eq 1 && "$(tail -1 "$_od/log")" == "$(_refusal tuple vep-outputs-compare repo-unbuilt)" && ! -e "$_od/engine.ran" && ! -s "$_started" ]] &&
        ok "a harness run in a tree without the tools, a built comparator exported in both variables: the same refusal, and the exported comparator never starts" ||
        bad "harness without the tools and a comparator exported: expected rc=1, the refusal and nothing started, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
    if /usr/bin/time -v true >/dev/null 2>&1; then
        _run_in_tree repo-unbuilt --walltime-only
        rc=$?
        [[ "$rc" -eq 0 && -e "$_od/engine.ran" && "$(grep -c 'is not an executable' "$_od/log")" -eq 0 ]] &&
            ok "a --walltime-only harness run in a tree without the tools: no refusal, the engine runs" ||
            bad "--walltime-only harness without the tools: expected rc=0 and an engine run, got rc=$rc '$(tail -2 "$_od/log")'"
        _run_in_tree repo-tuple-only --skip-fields
        rc=$?
        [[ "$rc" -eq 0 && -e "$_od/engine.ran" && "$(cat "$_started")" == "$_fdir/repo-tuple-only/$_veprs_tuple_start" ]] &&
            ok "a --skip-fields harness run in a tree without the field comparator: the cell runs and the tree's tuple comparator is the one comparator started" ||
            bad "--skip-fields harness without the field comparator: expected rc=0 and the tuple comparator's one start, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
        _run_in_tree repo
        rc=$?
        [[ "$rc" -eq 0 && "$(cat "$_started")" == "$_fdir/repo/$_veprs_tuple_start"$'\n'"$_fdir/repo/$_fields_start" ]] &&
            ok "a harness run that scores: the tree's tuple comparator starts on the cell's two output directories, then its field comparator on the two outputs, and nothing else" ||
            bad "scoring harness run: expected rc=0 and the two tree comparators' starts, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
        (
            export VEP_OUTPUTS_COMPARATOR=python VEP_FIELDS_COMPARATOR="$_fdir/bin/comparator"
            _run_in_tree repo
        )
        rc=$?
        [[ "$rc" -eq 0 && "$(cat "$_started")" == "$_fdir/repo/$_veprs_tuple_start"$'\n'"$_fdir/repo/$_fields_start" ]] &&
            ok "a harness run that scores with VEP_OUTPUTS_COMPARATOR=python and a built VEP_FIELDS_COMPARATOR exported: the same two starts of the tree's comparators" ||
            bad "scoring harness run with the variables exported: expected rc=0 and the two tree comparators' starts, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
        # The same cell by the harness of this tree, on the tree's own release builds: the
        # two binaries take the command lines above and write the cell's reports.
        _real_root="$(cd "$(dirname "$HARNESS")/../.." && pwd)"
        if [[ -x "$_real_root/$_outputs_rel" && -x "$_real_root/$_fields_rel" ]]; then
            rm -rf "$_od/real"
            bash "$HARNESS" --engine vep-rs --vep-binary "$_od/bin/vep" --data-dir "$_od/data" --output-dir "$_od/real" --datasets s01 --reference-release 116.2 >"$_od/real.log" 2>&1
            rc=$?
            _run_kinds=$(awk -F, 'NR>1 {print $8 ":" $10}' "$_od/real/wall_times.csv" 2>/dev/null | paste -sd' ' -)
            [[ "$rc" -eq 0 && "$_run_kinds" == "warmup:0 timed:0" ]] &&
                ok "the cell scored on the tree's own release builds: exit 0, and wall_times.csv carries the warmup and timed rows and no failed comparison" ||
                bad "cell on the tree's builds: expected rc=0 and rows 'warmup:0 timed:0', got rc=$rc rows '$_run_kinds' '$(tail -2 "$_od/real.log")'"
            if python3 - "$_od/real/s01/report" <<'PY'
import json, sys
report = sys.argv[1]
aggregate = json.load(open(report + "/summary.json"))["aggregate"]
counts = (aggregate["perl_tuple_count"], aggregate["rust_tuple_count"], aggregate["intersection_count"])
assert counts == (3, 3, 2), counts
PY
            then
                ok "the tree's tuple comparator writes report/summary.json: three tuples a side, two shared"
            else
                bad "report/summary.json of the cell on the tree's builds does not carry 3, 3 and 2 tuples: $(tail -2 "$_od/real.log")"
            fi
            if python3 - "$_od/real/s01/report" <<'PY'
import json, sys
report = sys.argv[1]
fields = json.load(open(report + "/fields.json"))
totals = fields["totals"]
assert (totals["intersection"], totals["reference_only"], totals["scored_only"]) == (2, 1, 1), totals
existing = fields["fields"]["Existing_variation"]
assert (existing["compared"], existing["differ"]) == (2, 1), existing
PY
            then
                ok "the tree's field comparator writes report/fields.json: two shared rows, one differing in Existing_variation"
            else
                bad "report/fields.json of the cell on the tree's builds does not carry 2 shared rows with 1 Existing_variation difference: $(tail -2 "$_od/real.log")"
            fi
        else
            echo "  skip the cell scored on the tree's own builds: no release build of tools/vep-outputs-compare and tools/vep-fields-compare"
        fi
    else
        echo "  skip the harness runs that reach the engine: /usr/bin/time -v is not GNU time here"
    fi
fi

echo "TREE-NORMALISER: a fastVEP cell that scores rewrites its raw output with this tree's release build of the normaliser and nothing else"

# fastVEP through v0.4.0 writes `Transcript` in the Feature_type of an intergenic row and
# from v0.5.0 `-`, and the raw output of a fastVEP cell of any version is rewritten before
# the tuple comparator reads it. The rewrite is this tree's release build of
# tools/vep-fastvep-normalize, resolved with the comparators before the first engine run: a
# build that is missing, not executable or a directory stops the run with one line naming
# the build command, whatever FASTVEP_NORMALIZE a caller exported; a run that scores nothing
# and another engine ask for none; the one call site runs what was resolved on the timed
# run's raw output; and a rewrite that fails is recorded with nothing scored in its place.
if [[ -f "$HARNESS" ]]; then
    # The one line a normaliser that cannot run is refused with: $1 the tree.
    _normaliser_refusal() {
        printf '%s' "ERROR: [run_clone_measurement] the fastVEP output normaliser $_fdir/$1/$_normalise_rel is not an executable: build the tool with 'cargo build --release --locked --manifest-path $_fdir/$1/tools/vep-fastvep-normalize/Cargo.toml'"
    }
    _out=$(_run_comparators_block repo-tuple-only "" fastvep false)
    [[ "$_out" == "$(_normaliser_refusal repo-tuple-only)"$'\n'"rc=1" ]] &&
        ok "the fastVEP engine with the normaliser not built: one line names its build command, exit 1, and no interpreter is taken in its place" ||
        bad "normaliser not built: expected its refusal alone and rc=1, got '$_out'"
    for _shape in inert dir; do
        _out=$(_run_comparators_block "repo-normaliser-$_shape" "" fastvep false)
        [[ "$_out" == "$(_normaliser_refusal "repo-normaliser-$_shape")"$'\n'"rc=1" ]] &&
            ok "the normaliser's release path $([[ "$_shape" == dir ]] && echo "a directory" || echo "a file that is not executable"): refused with its build command, exit 1" ||
            bad "normaliser path $_shape: expected its refusal alone and rc=1, got '$_out'"
    done
    _out=$(export FASTVEP_NORMALIZE="$_fdir/bin/comparator" && _run_comparators_block repo-tuple-only "" fastvep false)
    [[ "$_out" == "$(_normaliser_refusal repo-tuple-only)"$'\n'"rc=1" ]] &&
        ok "the normaliser not built and FASTVEP_NORMALIZE exported as an executable file: the same refusal, and what the variable names is not taken" ||
        bad "normaliser not built with FASTVEP_NORMALIZE exported: expected its refusal alone and rc=1, got '$_out'"
    _out=$(export FASTVEP_NORMALIZE="$_fdir/bin/python3" && _run_comparators_block repo-normaliser-built "" fastvep false)
    [[ "$_out" == "tuple: $_fdir/repo-normaliser-built/$_outputs_rel"$'\n'"normaliser: $_fdir/repo-normaliser-built/$_normalise_rel"$'\n'"rc=0" ]] &&
        ok "the normaliser built and FASTVEP_NORMALIZE exported as another program: the tree's build is the command, with no warning" ||
        bad "normaliser built with FASTVEP_NORMALIZE exported: expected the tree's two builds alone, got '$_out'"
    _out=$(_run_comparators_block repo-unbuilt "" fastvep true)
    [[ "$_out" == "rc=0" ]] &&
        ok "the fastVEP engine under --walltime-only: neither the tuple comparator nor the normaliser is asked for" ||
        bad "fastVEP --walltime-only in a tree without the tools: expected no command and rc=0, got '$_out'"
    _out=$(_run_comparators_block repo-all "")
    [[ "$_out" == "tuple: $_fdir/repo-all/$_outputs_rel"$'\n'"fields: $_fdir/repo-all/$_fields_rel"$'\n'"rc=0" ]] &&
        ok "the vep-rs engine asks for no normaliser, in a tree that holds its build as in one that does not" ||
        bad "vep-rs in a tree with the normaliser built: expected the two comparators alone, got '$_out'"
    _normalise_calls=$(grep -c -e '^[[:space:]]*"\${FASTVEP_NORMALIZE\[@\]}" "\$raw" "\$dataset_dir/rust/output.txt" || {$' "$HARNESS") || _normalise_calls=0
    _normalise_set=$(grep -c -e 'FASTVEP_NORMALIZE=' "$HARNESS") || _normalise_set=0
    _normalise_other=$(grep -e 'FASTVEP_NORMALIZE' "$HARNESS" | grep -c -v -e '^FASTVEP_NORMALIZE=()$' -e '^[[:space:]]*FASTVEP_NORMALIZE=("\$normalize_bin")$' -e '^[[:space:]]*"\${FASTVEP_NORMALIZE\[@\]}" "\$raw" ') || _normalise_other=0
    [[ "$_normalise_calls" -eq 1 && "$_normalise_set" -eq 2 && "$_normalise_other" -eq 0 ]] &&
        ok "the one fastVEP call site runs the resolved normaliser on the raw output, and no other line sets the command or hands the rewrite to an interpreter" ||
        bad "expected 1 call site of the resolved normaliser, 2 lines that set it and no other line that names it, found $_normalise_calls, $_normalise_set and $_normalise_other"

    # The whole harness: a copy of it in each tree runs one fastVEP cell, s01, over the input
    # and the ground truth of TREE-COMPARATORS. The stub engine writes the ground truth's three
    # rows with `Transcript` in the Feature_type of the intergenic one, as fastVEP through
    # v0.4.0 does.
    _od_row3_fastvep='v3\t21:3000\tA\t-\t-\tTranscript\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n'
    # shellcheck disable=SC2059
    {
        printf "$_od_header" && printf "$_od_row1" rs1 && printf "$_od_row2" synonymous_variant && printf "$_od_row3_fastvep"
    } >"$_od/fastvep_output.txt"
    cat >"$_od/bin/fastvep" <<STUB
#!/usr/bin/env bash
# stub engine: marks its start and writes the prepared output at -o
out=""; prev=""
for a in "\$@"; do [[ "\$prev" == -o ]] && out="\$a"; prev="\$a"; done
: >"$_od/engine.ran"
cp "$_od/fastvep_output.txt" "\$out"
STUB
    chmod +x "$_od/bin/fastvep"
    # A fastVEP cell is refused without its transcript cache file (FASTVEP-CACHE below).
    mkdir -p "$_od/data/caches/fastvep"
    echo "stub transcript cache" >"$_od/data/caches/fastvep/grch37.bin"
    # _run_fastvep_in_tree <tree> [harness args...]: as _run_in_tree, the engine fastVEP.
    _run_fastvep_in_tree() {
        local tree="$_fdir/$1"
        shift
        mkdir -p "$tree/scripts/concordance"
        cp "$HARNESS" "$(dirname "$HARNESS")/harness_genomes.sh" "$(dirname "$HARNESS")/harness_genomes.tsv" "$tree/scripts/concordance/"
        rm -rf "$_od/out" "$_od/engine.ran"
        : >"$_started"
        bash "$tree/scripts/concordance/run_clone_measurement.sh" --engine fastvep --fastvep-binary "$_od/bin/fastvep" --data-dir "$_od/data" \
            --output-dir "$_od/out" --datasets s01 --reference-release 116.2 "$@" >"$_od/log" 2>&1
    }
    _normalise_start="$_normalise_rel $_od/out/s01/rust/raw.txt $_od/out/s01/rust/output.txt"
    _fastvep_tuple_start="$_tuple_start --glob output.txt"
    _run_fastvep_in_tree repo-tuple-only
    rc=$?
    [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_od/log")" -eq 1 && "$(tail -1 "$_od/log")" == "$(_normaliser_refusal repo-tuple-only)" && ! -e "$_od/engine.ran" && ! -e "$_od/out" && ! -s "$_started" ]] &&
        ok "a fastVEP harness run in a tree without the normaliser: exit 1 on its refusal line, no output directory made, neither the engine nor the tuple comparator started" ||
        bad "fastVEP harness without the normaliser: expected rc=1, the refusal as the one ERROR and last line and nothing started, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
    if /usr/bin/time -v true >/dev/null 2>&1; then
        _run_fastvep_in_tree repo-tuple-only --walltime-only
        rc=$?
        [[ "$rc" -eq 0 && -e "$_od/engine.ran" && "$(grep -c 'is not an executable' "$_od/log")" -eq 0 && ! -s "$_started" ]] &&
            ok "a --walltime-only fastVEP harness run in a tree without the normaliser: no refusal, the engine runs and nothing rewrites or scores" ||
            bad "--walltime-only fastVEP harness without the normaliser: expected rc=0, an engine run and nothing else started, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
        _run_fastvep_in_tree repo-normaliser-built
        rc=$?
        [[ "$rc" -eq 0 && "$(cat "$_started")" == "$_fdir/repo-normaliser-built/$_normalise_start"$'\n'"$_fdir/repo-normaliser-built/$_fastvep_tuple_start" ]] &&
            ok "a fastVEP harness run that scores: the tree's normaliser starts on the timed run's raw output and the cell's output.txt, then its tuple comparator on the cell's two output directories, and nothing else" ||
            bad "scoring fastVEP harness run: expected rc=0 and the starts of the tree's normaliser and tuple comparator, got rc=$rc '$(tail -2 "$_od/log")' and '$(cat "$_started")'"
        _run_fastvep_in_tree repo-normaliser-failing
        _failed_rows=$(awk -F, 'NR>1 && $8 == "normalizer_failed" && $10 == "1"' "$_od/out/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$(cat "$_started")" == "$_fdir/repo-normaliser-failing/$_normalise_start" && "$_failed_rows" -eq 1 ]] && grep -q -F -e "[ERROR] fastvep normalizer failed on s01; refusing to compare un-normalised output" "$_od/log" &&
            ok "a normaliser that fails: the cell is recorded as normalizer_failed with one line naming it, and the tuple comparator never starts on the raw output" ||
            bad "failing normaliser: expected its one start, the line naming the cell and one normalizer_failed row, got '$(cat "$_started")', $_failed_rows row(s) and '$(tail -3 "$_od/log" | tr '\n' '|' | cut -c1-300)'"
        # The same cell by the harness of this tree, on the tree's own release builds: the
        # normaliser writes the output the tuple comparator scores.
        if [[ -x "$_real_root/$_outputs_rel" && -x "$_real_root/$_normalise_rel" ]]; then
            rm -rf "$_od/real-fastvep"
            bash "$HARNESS" --engine fastvep --fastvep-binary "$_od/bin/fastvep" --data-dir "$_od/data" --output-dir "$_od/real-fastvep" --datasets s01 --reference-release 116.2 >"$_od/real-fastvep.log" 2>&1
            rc=$?
            _run_kinds=$(awk -F, 'NR>1 {print $8 ":" $10}' "$_od/real-fastvep/wall_times.csv" 2>/dev/null | paste -sd' ' -)
            [[ "$rc" -eq 0 && "$_run_kinds" == "warmup:0 timed:0" ]] &&
                ok "the fastVEP cell scored on the tree's own release builds: exit 0, and wall_times.csv carries the warmup and timed rows and no failed step" ||
                bad "fastVEP cell on the tree's builds: expected rc=0 and rows 'warmup:0 timed:0', got rc=$rc rows '$_run_kinds' '$(tail -2 "$_od/real-fastvep.log")'"
            cmp -s "$_od/real-fastvep/s01/rust/output.txt" "$_od_gt/output.txt" && ! cmp -s "$_od/real-fastvep/s01/rust/raw.txt" "$_od_gt/output.txt" &&
                grep -q -x -F -e "normalize: rewrote 1/3 intergenic Feature_type=Transcript → -" "$_od/real-fastvep.log" &&
                ok "the tree's normaliser writes rust/output.txt: one row of three rewritten, the intergenic row's Feature_type the ground truth's" ||
                bad "rust/output.txt of the fastVEP cell on the tree's builds is not the raw output with its intergenic row rewritten: $(tail -2 "$_od/real-fastvep.log")"
            if python3 - "$_od/real-fastvep/s01/report" <<'PY'
import json, sys
report = sys.argv[1]
aggregate = json.load(open(report + "/summary.json"))["aggregate"]
counts = (aggregate["perl_tuple_count"], aggregate["rust_tuple_count"], aggregate["intersection_count"])
assert counts == (3, 3, 3), counts
PY
            then
                ok "the tree's tuple comparator scores the rewritten output: three tuples a side, all three shared"
            else
                bad "report/summary.json of the fastVEP cell on the tree's builds does not carry 3, 3 and 3 tuples: $(tail -2 "$_od/real-fastvep.log")"
            fi
        else
            echo "  skip the fastVEP cell scored on the tree's own builds: no release build of tools/vep-outputs-compare and tools/vep-fastvep-normalize"
        fi
    else
        echo "  skip the fastVEP harness runs that reach the engine: /usr/bin/time -v is not GNU time here"
    fi
fi

echo "CONTAINER-PROBE: the container start-up probe times bare image starts, Perl only, after the datasets"

# Every timed Ensembl VEP invocation pays one container start; --container-probe N measures that
# cost on its own (N bare `docker run --rm <image> true` starts under /usr/bin/time -v, rows
# in <output-dir>/container_start.csv). The flag is Ensembl VEP-only, defaults to 0, runs after the
# dataset loop, and a non-zero probe exit is recorded rather than hidden.
if [[ -f "$HARNESS" ]]; then
    grep -q '^CONTAINER_PROBE=0$' "$HARNESS" &&
        ok "--container-probe defaults to 0 (no probe)" ||
        bad "CONTAINER_PROBE does not default to 0"
    grep -A 2 -e '^[[:space:]]*--container-probe)$' "$HARNESS" | grep 'CONTAINER_PROBE="$2"' >/dev/null &&
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
    _last_measure=$(grep -n -e 'measure_sv "$dataset_id" "$dataset_name" "$assembly"' "$HARNESS" | tail -1 | cut -d: -f1)
    if [[ -n "$_probe_def" && -n "$_last_measure" && "$_last_measure" -lt "$_probe_def" ]]; then
        ok "the probe is defined and run after the dataset loop (loop measure at line $_last_measure, probe at $_probe_def)"
    else
        bad "the probe does not follow the dataset loop (measure at ${_last_measure:-none}, probe at ${_probe_def:-none})"
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
        echo "$_probe_out" | grep 'container probe 3 exited 7' >/dev/null &&
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

# The engine_version every Ensembl VEP row carries is the image tag, and a tag can move. The harness
# pulls the image when absent and refuses one whose RepoDigests lack the expected digest, so a
# moved tag cannot be timed under the pinned version. Run against a stub docker whose inspect
# answers with a configurable digest; an empty data tree is enough, because the check precedes
# every cell.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"
if [[ -f "$HARNESS" ]]; then
    PINNED_DIGEST="sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
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
    # An Ensembl VEP run asks for the release build of tools/vep-file-census in its tree before any
    # cell runs. These runs reach no pass over a whole file, so they run a copy of the
    # harness in a tree of its own, whose build records a start and counts nothing.
    _pt="$TMP/perl_tree"
    _census_rel="tools/vep-file-census/target/release/vep-file-census"
    mkdir -p "$_pt/scripts/concordance" "$_pt/$(dirname "$_census_rel")"
    cp "$HARNESS" "$(dirname "$HARNESS")/harness_genomes.sh" "$(dirname "$HARNESS")/harness_genomes.tsv" "$_pt/scripts/concordance/"
    printf '#!/usr/bin/env bash\nprintf "%%s\\n" "$*" >>"%s"\n' "$_pt/census.started" >"$_pt/$_census_rel"
    chmod +x "$_pt/$_census_rel"
    PERL_HARNESS="$_pt/scripts/concordance/run_clone_measurement.sh"
    run_perl_harness() { # run_perl_harness <stub digest> <out dir> [harness args...]
        local d="$1" o="$2"; shift 2
        : >"$_dd/calls"
        DIGEST_CALLS="$_dd/calls" STUB_DIGEST="$d" PATH="$_dd/bin:$PATH" \
            bash "$PERL_HARNESS" --engine perl --data-dir "$_dd/data" --output-dir "$o" --datasets s01 "$@" >"$o.log" 2>&1
    }
    grep -q '^PERL_IMAGE="ensemblorg/ensembl-vep:release_116.3"$' "$HARNESS" &&
        ok "the default Perl image is the release_116.3 tag" ||
        bad "the default Perl image is not ensemblorg/ensembl-vep:release_116.3"
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
        ok "the mismatch is caught before the first dataset and before any docker run"
    else
        bad "the harness went on past the digest mismatch"
    fi
    run_perl_harness "$PINNED_DIGEST" "$_dd/out2"; rc=$?
    if grep -q "run_clone_measurement engine=perl" "$_dd/out2.log" && ! grep -q "digest mismatch" "$_dd/out2.log"; then
        ok "the pinned digest passes the check and the datasets start"
    else
        bad "the pinned digest did not pass the check: $(head -3 "$_dd/out2.log" | tr '\n' ' ')"
    fi
    grep -q '^docker image inspect --format {{join .RepoDigests ","}} ensemblorg/ensembl-vep:release_116.3$' "$_dd/calls" &&
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

# --cells takes a comma-separated cell list (s01,s02,s07,s08,s09,s10), the comma form of --datasets;
# the same stub docker and empty data tree as IMAGE-DIGEST, since the filter is applied before
# any cell's inputs are read.
if [[ -f "$HARNESS" ]]; then
    _cd="$TMP/cells"; mkdir -p "$_cd/data"
    run_cells_harness() { # run_cells_harness <out name> [harness args...]
        local o="$1"; shift
        : >"$_dd/calls"
        DIGEST_CALLS="$_dd/calls" STUB_DIGEST="$PINNED_DIGEST" PATH="$_dd/bin:$PATH" \
            bash "$PERL_HARNESS" --engine perl --data-dir "$_cd/data" --output-dir "$_cd/$o" "$@" >"$_cd/$o.log" 2>&1
    }
    run_cells_harness out1 --cells s07,s08,s09,s10
    if grep -q "run_clone_measurement engine=perl datasets=s07 s08 s09 s10 walltime_only" "$_cd/out1.log" &&
        grep -q "dataset s01 skipped (datasets=s07 s08 s09 s10)" "$_cd/out1.log" && grep -q "dataset s06 skipped" "$_cd/out1.log" &&
        grep -q "dataset s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out1.log" && ! grep -q "dataset s07 skipped" "$_cd/out1.log"; then
        ok "--cells s07,s08,s09,s10 selects those four cells and skips s01 to s06"
    else
        bad "--cells did not select the named cells: $(grep -E 'datasets=|skipped|starting' "$_cd/out1.log" | head -4 | tr '\n' ' ')"
    fi
    run_cells_harness out2 --cells s07,s99; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: \[run_clone_measurement\] unknown cell 's99' in --datasets/--cells; the cells are: s01 s02 s03 s04 s05 s06 s07 s08 s09 s10 s05m s06m" "$_cd/out2.log" &&
        ! grep -q "dataset s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out2.log" && ! grep -q "^docker run" "$_dd/calls"; then
        ok "an unknown cell id stops the run before any cell starts and names the cells the table carries"
    else
        bad "an unknown cell id was not refused (rc=$rc): $(head -3 "$_cd/out2.log" | tr '\n' ' ')"
    fi
    run_cells_harness out3 --datasets "s01 s0x"; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "unknown cell 's0x'" "$_cd/out3.log"; then
        ok "--datasets with an unknown id is refused the same way (the id would have matched nothing and exited 0)"
    else
        bad "--datasets with an unknown id was accepted (rc=$rc)"
    fi
    grep -q 'DATASETS_FILTER="${2//,/ }"' "$HARNESS" && ! grep -qE '^cells[a-z_]*\(\) \{' "$HARNESS" &&
        ok "--cells is an inline rewrite of the comma list into the --datasets filter, with no column-0 function of its own" ||
        bad "--cells is not the inline comma-to-space rewrite, or it has a column-0 function of its own"
    run_cells_harness out4 --cells ""; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: \[run_clone_measurement\] empty cell list in --datasets/--cells" "$_cd/out4.log" &&
        ! grep -q "skipped" "$_cd/out4.log"; then
        ok "an empty cell list is refused (it would have skipped every cell and exited 0)"
    else
        bad "an empty cell list was accepted (rc=$rc): $(head -3 "$_cd/out4.log" | tr '\n' ' ')"
    fi
    # A glob-shaped id is validated as text: a file named like a cell in the working directory must not satisfy it.
    mkdir -p "$_cd/cwd" && : >"$_cd/cwd/s01"
    DIGEST_CALLS="$_dd/calls" STUB_DIGEST="$PINNED_DIGEST" PATH="$_dd/bin:$PATH" \
        bash -c 'cd "$1" && bash "$2" --engine perl --data-dir "$3" --output-dir "$4" --datasets "s0*"' _ "$_cd/cwd" "$PERL_HARNESS" "$_cd/data" "$_cd/out5" >"$_cd/out5.log" 2>&1; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "unknown cell 's0\*'" "$_cd/out5.log"; then
        ok "a glob-shaped cell id is refused as text, not expanded against the working directory"
    else
        bad "a glob-shaped cell id was accepted (rc=$rc): $(head -3 "$_cd/out5.log" | tr '\n' ' ')"
    fi
    # The validation loop splits on any whitespace while should_run matches between single spaces:
    # a tab-separated list must run its cells, not validate and then skip every one with exit 0.
    run_cells_harness out6 --datasets $'s07\ts08'
    if grep -q "run_clone_measurement engine=perl datasets=s07 s08 walltime_only" "$_cd/out6.log" &&
        grep -q "dataset s01 skipped (datasets=s07 s08)" "$_cd/out6.log" &&
        grep -q "dataset s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out6.log" && ! grep -q "dataset s07 skipped" "$_cd/out6.log"; then
        ok "a tab-separated --datasets list is normalised to one space per gap, so its cells run"
    else
        bad "a tab-separated --datasets list validated but did not run its cells: $(grep -E 'datasets=|skipped|starting' "$_cd/out6.log" | head -3 | tr '\n' ' ')"
    fi
    # A repeated id would pass the known-id check, and should_run's substring match would run
    # its cell once under a list that names it twice.
    run_cells_harness out7 --cells s07,s07,s08; rc=$?
    if [[ "$rc" -ne 0 ]] && grep -q "^ERROR: \[run_clone_measurement\] repeated cell 's07' in --datasets/--cells" "$_cd/out7.log" &&
        ! grep -q "dataset s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_cd/out7.log" && ! grep -q "^docker run" "$_dd/calls"; then
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
DATASET="$TMP/dataset"
mkdir -p "$DATASET/rust" "$DATASET/perl_scratch" "$DATASET/perf"
cat >"$TMP/engine_warm.sh" <<'EOF'
#!/usr/bin/env bash
# $1 = run kind, $2 = dataset dir
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
    "$TMP/engine_warm.sh" "$run_kind" "$DATASET"
    [[ "$run_kind" == warmup ]] && discard_warmup_output "$DATASET"
done
leftover=$(cat "$DATASET/perf/leftover_at_timed_start")
[[ "$leftover" -eq 0 ]] &&
    ok "the timed run started with no warmup output on disk" ||
    bad "the timed run found $leftover warmup file(s) still on disk"
[[ -s "$DATASET/rust/output.txt" ]] &&
    ok "the timed run's own output survives the step" ||
    bad "the timed run's output is missing after the step"

# The same loop without the step is the shape the test exists to catch.
rm -f "$DATASET"/rust/* "$DATASET"/perl_scratch/* "$DATASET/perf/leftover_at_timed_start"
for run_kind in warmup timed; do
    "$TMP/engine_warm.sh" "$run_kind" "$DATASET"
done
leftover=$(cat "$DATASET/perf/leftover_at_timed_start")
[[ "$leftover" -eq 3 ]] &&
    ok "without the step the timed run starts over the warmup's 3 files (the fault the step removes)" ||
    bad "without the step the timed run found $leftover file(s), expected 3"

echo
echo "PERL-SV-CACHE-VERSION: the Perl SV cell's child shell must receive the cache version it passes to vep"

# The fault: the per-VCF Ensembl VEP loop runs in a single-quoted `bash -c` child under `set -u`.
# A harness variable the parent never exported is unbound there, so a child that reads
# `$PERL_CACHE_VERSION` aborts on its first record with "unbound variable", every Ensembl VEP SV
# cell records a non-zero exit, and the campaign has no Ensembl VEP SV wall time. The value must
# travel as a positional argument like the other ten.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"
if [[ -f "$HARNESS" ]]; then
    # The child body: the single-quoted `bash -c` block of the Ensembl VEP branch of measure_sv,
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
    printf '%s\n' "$_perl_args" | grep '"\$FORK" "\$cache_version" "\$species")' >/dev/null &&
        ok "the Perl SV subshell receives the genome's cache version as its tenth argument" ||
        bad "the Perl SV subshell's argument list does not pass \$cache_version between \$FORK and \$species"
    _cv_assigned=$(sed -n '/^measure_sv() {/,/^}/p' "$HARNESS" | grep -c -e '^[[:space:]]*cache_version=\$(cache_version_for "\$species" "\$assembly")$') || _cv_assigned=0
    [[ "$_cv_assigned" -eq 1 ]] &&
        ok "measure_sv reads that cache version from cache_version_for, under the cell's species and assembly" ||
        bad "measure_sv assigns cache_version from cache_version_for on $_cv_assigned line(s); expected 1"

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
        ok "no command line of the harness passes a literal --species homo_sapiens; the dataset row's species travels as an argument" ||
        bad "$_species_literal command line(s) pass a literal --species homo_sapiens"
fi

echo
echo "REFERENCE-RELEASE: every comparator call scores against the release the ground truth names"

# The comparators default to a 116 release's registry. A cell scored against a 115.2 ground truth
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
    GT116_3="$TMP/gt116_3"
    GT115="$TMP/gt115"
    GTNONE="$TMP/gtnone"
    mkdir -p "$GT116" "$GT116_3" "$GT115" "$GTNONE"
    printf '{"perl_image": "ensemblorg/ensembl-vep:release_116.2", "cache_version": 116}\n' >"$GT116/output.provenance.json"
    printf '{"perl_image": "ensemblorg/ensembl-vep:release_116.3", "cache_version": 116}\n' >"$GT116_3/output.provenance.json"
    printf '{"perl_image": "ensemblorg/ensembl-vep@sha256:0000", "cache_version": 115}\n' >"$GT115/provenance.json"
    printf '{"fork": 4, "buffer_size": 5000}\n' >"$GTNONE/output.provenance.json"
    if declare -f reference_release_for >/dev/null; then
        REFERENCE_RELEASE=""
        got=$(reference_release_for "$GT116" 2>/dev/null)
        [[ "$got" == "116.2" ]] &&
            ok "a ground truth whose sidecar names release_116.2 scores under 116.2" ||
            bad "release_116.2 sidecar resolved to '$got'"
        got=$(reference_release_for "$GT116_3" 2>/dev/null)
        [[ "$got" == "116.3" ]] &&
            ok "a ground truth whose sidecar names release_116.3 scores under 116.3, beside the same cache version" ||
            bad "release_116.3 sidecar resolved to '$got'"
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
        # A host without python3 cannot run the derivation: the
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
    # dataset directories that hold only outputs (the shape of the published SNP/indel set),
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
            ok "layout 1, one record at the set root over output-only dataset directories, resolves (115.2)" ||
            bad "the set-root record did not resolve for s01 (got '$got')"
        got=$(reference_release_for "$SET/s03" 2>/dev/null)
        [[ "$got" == "116.2" ]] &&
            ok "layout 2, a sidecar beside the dataset's output, resolves and wins over the set record (116.2)" ||
            bad "the dataset sidecar did not win (got '$got')"
        if reference_release_for "$TMP/bare_set/s01" >/dev/null 2>"$TMP/rr_bare.err"; then
            bad "layout 3, neither record, resolved to a release"
        else
            ok "layout 3, neither a set record nor a sidecar, fails the cell"
        fi
    fi
    # The run's own --reference-release is checked first of the arguments: an accepted release
    # reaches the next check, the missing engine, and any other stops the run on its own line.
    for _rel in 115.2 116.2 116.3; do
        bash "$HARNESS" --reference-release "$_rel" >"$TMP/rr_arg.log" 2>&1
        [[ "$(head -1 "$TMP/rr_arg.log")" == "ERROR: --engine required" ]] &&
            ok "--reference-release $_rel is accepted" ||
            bad "--reference-release $_rel stopped the run on: $(head -1 "$TMP/rr_arg.log")"
    done
    bash "$HARNESS" --reference-release 116.1 >"$TMP/rr_arg.log" 2>&1; rc=$?
    if [[ "$rc" -ne 0 && "$(head -1 "$TMP/rr_arg.log")" == "ERROR: --reference-release must be 115.2, 116.2 or 116.3 (got '116.1')" ]]; then
        ok "--reference-release 116.1 stops the run on a line naming the three accepted releases"
    else
        bad "--reference-release 116.1: rc=$rc, first line: $(head -1 "$TMP/rr_arg.log")"
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
    grep -A 2 -e '^[[:space:]]*--vep-extra-args)$' "$HARNESS" | grep 'VEP_EXTRA_ARGS="$2"' >/dev/null &&
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
    # A fastVEP cell is refused without its transcript cache file; one run below is fastVEP's.
    mkdir -p "$_ea/data/caches/fastvep"
    echo "stub transcript cache" >"$_ea/data/caches/fastvep/grch37.bin"
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
            --output-dir "$_ea/$run_name" --datasets "s01 s07" --walltime-only "$@" >"$_ea/$run_name.log" 2>&1
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
            --data-dir "$_ea/data" --output-dir "$_ea/rel" --datasets s01 --walltime-only >"$_ea/rel.log" 2>&1); rc=$?
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
        STUB_ARGV="$_ea/rerun-fastvep.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_ea/bin/vep" --data-dir "$_ea/data" --output-dir "$_ea/rerun" --datasets s01 --walltime-only >"$_ea/rerun-fastvep.log" 2>&1; rc=$?
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
    bash "$HARNESS" --engine vep-rs --vep-binary "$_ea/bin" --data-dir "$_ea/data" --output-dir "$_ea/dirbin" --datasets s01 --walltime-only >"$_ea/dirbin.log" 2>&1; rc=$?
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
        bash "$HARNESS" --engine vep-rs --vep-binary "$_ea/bin/vep-xonly" --data-dir "$_ea/data" --output-dir "$_ea/xonly" --datasets s01 --walltime-only >"$_ea/xonly.log" 2>&1; rc=$?
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
    PATH="$_np" "$_np/bash" "$HARNESS" --engine vep-rs --vep-binary "$_ea/bin/vep" --data-dir "$_ea/data" --output-dir "$_ea/nopython/out" --datasets s01 --walltime-only >"$_ea/nopython.log" 2>&1; rc=$?
    _n_err=$(grep -c '^ERROR: \[run_clone_measurement\] python3 is not on PATH' "$_ea/nopython.log") || _n_err=0
    if [[ "$rc" -eq 1 && "$_n_err" -eq 1 ]] && ! grep -q -e 'command not found' -e ' starting$' "$_ea/nopython.log"; then
        ok "without python3 a vep-rs run stops before any cell with one line naming the interpreter"
    else
        bad "without python3: rc=$rc, $_n_err interpreter line(s): $(head -3 "$_ea/nopython.log" | tr '\n' '|')"
    fi
    bash "$HARNESS" --engine perl --data-dir "$_ea/data" --output-dir "$_ea/perl" --datasets s01 --vep-extra-args "--transcript_index_impl coitrees" >"$_ea/perl.log" 2>&1; rc=$?
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
    grep -A 2 -e '^[[:space:]]*--fastvep-extra-args)$' "$HARNESS" | grep 'FASTVEP_EXTRA_ARGS="$2"' >/dev/null &&
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
    # A fastVEP cell is refused without its transcript cache file (FASTVEP-CACHE below).
    mkdir -p "$_fa/data/caches/fastvep"
    echo "stub transcript cache" >"$_fa/data/caches/fastvep/grch37.bin"
    echo "stub transcript cache" >"$_fa/data/caches/fastvep/grch38.bin"
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
    run_fastvep_extra() { # run_fastvep_extra <out name> [harness args...]; the cells are s01 and s07 unless --datasets or --cells follows
        local o="$1"; shift
        : >"$_fa/$o.argv"
        STUB_ARGV="$_fa/$o.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" \
            --output-dir "$_fa/$o" --datasets "s01 s07" --walltime-only "$@" >"$_fa/$o.log" 2>&1
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
        run_fastvep_extra twoasm --datasets "s01 s02" --fastvep-extra-args "$_fwords"; rc=$?
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
        run_fastvep_extra eqplain --output-dir "$_fa/eq=plain/out" --datasets s01; rc=$?
        _n_inv=$(grep -c -e " -o $_fa/eq=plain/out/s01/rust/raw.txt --output-format tab\$" "$_fa/eqplain.argv") || _n_inv=0
        _bad_rows=$(awk -F, 'NR>1 && $10!="0"' "$_fa/eq=plain/out/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$rc" -eq 0 && "$_n_inv" -eq 2 && "$_bad_rows" -eq 0 ]] &&
            ok "without the words an output directory whose path holds = runs its cell, both rows exit 0" ||
            bad "no words, = in the output directory: rc=$rc, $_n_inv of 2 invocations, $_bad_rows failed row(s): $(tail -2 "$_fa/eqplain.log" | tr '\n' '|' | cut -c1-240)"
        # With words the wrapper is the operand `env` reads and it execs the binary, so a binary
        # whose own path holds `=` runs.
        run_fastvep_extra eqbinwords --fastvep-binary "$_fa/eq=bin/fastvep" --datasets s01 --fastvep-extra-args "$_fwords"; rc=$?
        _n_tail=$(grep -c -e " --output-format tab $_fwords\$" "$_fa/eqbinwords.argv") || _n_tail=0
        _bad_rows=$(awk -F, 'NR>1 && $10!="0"' "$_fa/eqbinwords/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$rc" -eq 0 && "$_n_tail" -eq 2 && "$_bad_rows" -eq 0 ]] &&
            ok "with the words a fastVEP binary whose path holds = runs through the wrapper, both rows exit 0" ||
            bad "words, = in the binary's path: rc=$rc, $_n_tail of 2 invocations, $_bad_rows failed row(s): $(tail -2 "$_fa/eqbinwords.log" | tr '\n' '|' | cut -c1-240)"
        run_fastvep_extra handover --fastvep-extra-args "$_fwords"; rc1=$?
        STUB_ARGV="$_fa/handover-veprs.argv" bash "$HARNESS" --engine vep-rs --vep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/handover" --datasets s01 --walltime-only >"$_fa/handover-veprs.log" 2>&1; rc2=$?
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
            --data-dir "$_fa/data" --output-dir "$_fa/bare" --datasets s01 --walltime-only >"$_fa/bare.log" 2>&1); rc=$?
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
        bash "$HARNESS" --engine "$_eng" --vep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/refused-$_eng" --datasets s01 --walltime-only --fastvep-extra-args "$_fwords" >"$_fa/refused-$_eng.log" 2>&1; rc=$?
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
        LC_ALL="$_word_locale" STUB_ARGV="$_fa/word-$_wi.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/word-$_wi" --datasets s01 --walltime-only --fastvep-extra-args "$_value" </dev/null >"$_fa/word-$_wi.log" 2>&1; rc=$?
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
    STUB_ARGV="$_fa/eqdir.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/eq=refused/out" --datasets s01 --walltime-only --fastvep-extra-args "$_fwords" >"$_fa/eqdir.log" 2>&1; rc=$?
    _n_err=$(grep -c '^ERROR: ' "$_fa/eqdir.log") || _n_err=0
    if [[ "$rc" -eq 1 && "$_n_err" -eq 1 && ! -e "$_fa/eq=refused" && ! -s "$_fa/eqdir.argv" ]] && grep -qF -e "ERROR: [run_clone_measurement] --output-dir '$_fa/eq=refused/out' holds '=': with --fastvep-extra-args" "$_fa/eqdir.log" && ! grep -q "run_clone_measurement engine=" "$_fa/eqdir.log"; then
        ok "with the words an output directory whose path holds = is refused by name with one ERROR line, before it is made"
    else
        bad "words and = in the output directory: rc=$rc, $_n_err ERROR line(s), directory made=$([[ -e "$_fa/eq=refused" ]] && echo yes || echo no): $(head -2 "$_fa/eqdir.log" | tr '\n' '|' | cut -c1-240)"
    fi
    # Without words the operand `env` reads is the binary's own path, so a binary whose path
    # holds `=` is refused by name the same way, before the output directory is made.
    : >"$_fa/eqbin.argv"
    STUB_ARGV="$_fa/eqbin.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/eq=bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/eqbin" --datasets s01 --walltime-only >"$_fa/eqbin.log" 2>&1; rc=$?
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
            --data-dir "$_fa/data" --output-dir "$_fa/$on_path-out" --datasets s01 --walltime-only >"$_fa/$on_path.log" 2>&1); rc=$?
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
    refused_sa_dir sa-list "a list of one GRCh37 and one GRCh38 cell" "$_fwords --sa-dir /data/sa/grch38" --datasets "s01 s02"
    refused_sa_dir sa-cells "the comma form naming both SV cells, in the --sa-dir=<directory> form" "--sa-dir=/data/sa/grch38" --cells s07,s08
    refused_sa_dir sa-all "the default selection of every cell" "--sa-dir /data/sa/grch38"
    # The record of a run with the words is written through python3, so a host without it is
    # refused at startup by the harness's own line. The PATH holds the tools the harness runs
    # before that check (tr normalises the cell list) and no interpreter.
    _fnp="$_fa/nopython/bin"
    mkdir -p "$_fnp"
    for _t in bash dirname date mkdir tr; do ln -sf "$(type -P "$_t")" "$_fnp/$_t"; done
    PATH="$_fnp" "$_fnp/bash" "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fa/data" --output-dir "$_fa/nopython/out" --datasets s01 --walltime-only --fastvep-extra-args "$_fwords" >"$_fa/nopython.log" 2>&1; rc=$?
    _n_err=$(grep -c '^ERROR: \[run_clone_measurement\] python3 is not on PATH: a fastVEP run with --fastvep-extra-args' "$_fa/nopython.log") || _n_err=0
    if [[ "$rc" -eq 1 && "$_n_err" -eq 1 && ! -e "$_fa/nopython/out" ]] && ! grep -q -e 'command not found' -e ' starting$' "$_fa/nopython.log"; then
        ok "without python3 a fastVEP run with the words stops before any cell with one line naming the interpreter"
    else
        bad "without python3: rc=$rc, $_n_err interpreter line(s): $(head -3 "$_fa/nopython.log" | tr '\n' '|')"
    fi
fi

echo
echo "REFUSED-RUN: a run refused on its arguments leaves the output directory as it found it"

# Each run below is refused for its cell selection, the last refusal before SV-COMPARATOR's, twice:
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
    refused_run_leaves sa-dir "--sa-dir over cells of two assemblies" "--fastvep-extra-args carries --sa-dir, a directory of one assembly's databases" --datasets "s01 s02" --fastvep-extra-args "--sa-dir /data/sa/grch38"
    refused_run_leaves unknown "an unknown cell id" "unknown cell 's99' in --datasets/--cells" --cells s07,s99
    refused_run_leaves empty "an empty cell list" "empty cell list in --datasets/--cells" --cells ""
    refused_run_leaves repeated "a repeated cell id" "repeated cell 's07' in --datasets/--cells" --cells s07,s07,s08
fi

echo
echo "GENOME-TABLE: a cell's genome reads its FASTA, vep-rs cache and cache version from its row of the genome table"

# The table is harness_genomes.tsv beside the harness, or the file VEP_HARNESS_GENOMES names.
# A genome other than the two human ones needs a dataset row to be a cell at all, and the
# shipped dataset table is human, so the copy under test is the shipped harness with two dataset
# rows appended for a made-up genome (one SNP/indel cell, one SV cell) and nothing else
# changed. A stub engine and a stub docker record the argument vectors they are started with.
if [[ -f "$HARNESS" ]]; then
    _gt="$TMP/genome_table"
    _gt_data="$_gt/data"
    _gt_table_header='species	assembly	cache_version	fasta	json_cache'
    _gt_human_rows='homo_sapiens	GRCh37	116	reference/grch37/genome.fa	caches/vep-rs/grch37
homo_sapiens	GRCh38	116	reference/grch38/genome.fa	caches/vep-rs/grch38'
    _gt_made_up_row='made_up_species	ASM1	63	reference/asm1/made_up.fa	caches/vep-rs/asm1'
    _gt_committed="$(dirname "$HARNESS")/harness_genomes.tsv"
    # The digest the stub docker answers with and the Ensembl VEP runs ask for: any one value
    # serves.
    _gt_digest="sha256:$(printf '%064d' 0)"
    mkdir -p "$_gt/copy/scripts/concordance" "$_gt/bin" "$_gt_data/vcf/asm1/all_variants_canonical" \
        "$_gt_data/vcf/grch37/all_variants_canonical" "$_gt_data/ground_truth/perl/sv_per_vcf/canonical_inputs/asm1" \
        "$_gt_data/caches/perl/vep-cache/made_up_species/63_ASM1/I" "$_gt_data/caches/perl/vep-cache/made_up_species/115_ASM1/I"
    _gt_copy="$_gt/copy/scripts/concordance/run_clone_measurement.sh"
    awk '{ print } /^\t"s06m:1KG_highcov_chr21_multisample:GRCh38:snp_indel:1kg_highcov_chr21_multisample:homo_sapiens"$/ {
        print "\t\"s90:Made_up:ASM1:snp_indel:made_up:made_up_species\""
        print "\t\"s91:Made_up_SV:ASM1:sv:-:made_up_species\""
    }' "$HARNESS" >"$_gt_copy"
    cp "$(dirname "$HARNESS")/harness_genomes.sh" "$_gt/copy/scripts/concordance/"
    # An Ensembl VEP run asks for the release build of tools/vep-file-census in its tree and counts
    # its rows and records through it, so the copy's tree holds a stand-in made of the
    # pipelines those two counts stand for.
    _gt_census="$_gt/copy/tools/vep-file-census/target/release/vep-file-census"
    mkdir -p "$(dirname "$_gt_census")"
    cat >"$_gt_census" <<'STUB'
#!/usr/bin/env bash
case "$1 $2" in
"count --gunzip") gzip -dc "$5" | grep -vc '^#' ;;
"count --not-starting") grep -vc '^#' "$4" ;;
*) echo "stand-in vep-file-census: unexpected arguments: $*" >&2; exit 2 ;;
esac
STUB
    chmod +x "$_gt_census"
    # diff exits 1 on files that differ, so its output goes to a file and not through a pipe.
    diff "$HARNESS" "$_gt_copy" >"$_gt/copy.diff"
    _added=$(grep -c '^[<>]' "$_gt/copy.diff") || _added=0
    _added_rows=$(grep -c -e '^> 	"s9[01]:Made_up' "$_gt/copy.diff") || _added_rows=0
    [[ "$_added" -eq 2 && "$_added_rows" -eq 2 ]] &&
        ok "the copy under test is the shipped harness plus two dataset rows for the made-up genome" ||
        bad "the copy differs from the shipped harness by $_added line(s), $_added_rows of them the made-up dataset rows; expected 2 and 2"

    : >"$_gt_data/vcf/asm1/all_variants_canonical/made_up.canonical.vcf.gz"
    : >"$_gt_data/vcf/grch37/all_variants_canonical/clinvar_grch37.canonical.vcf.gz"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\nI\t1\n' | gzip -c >"$_gt_data/ground_truth/perl/sv_per_vcf/canonical_inputs/asm1/01_a.vcf.gz"
    cat >"$_gt/bin/vep" <<'STUB'
#!/usr/bin/env bash
# stub engine: one line per invocation, the argument vector joined by spaces; one data line at -o
out=""; prev=""
for a in "$@"; do [[ "$prev" == -o ]] && out="$a"; prev="$a"; done
printf '%s\n' "$*" >>"$STUB_ARGV"
printf '#Uploaded_variation\tLocation\nv\tI:1\n' >"$out"
STUB
    cat >"$_gt/bin/docker" <<'STUB'
#!/usr/bin/env bash
# stub docker: answers the image inspection with the pinned digest, records every `run`
# and writes two annotated rows at the output the run names
case "$1" in
image) echo "ensemblorg/ensembl-vep@$STUB_DIGEST"; exit 0 ;;
run)
    printf '%s\n' "$*" >>"$STUB_ARGV"
    out_dir=""; out_name=""; prev=""
    for a in "$@"; do
        case "$a" in *:/work/perl_out) out_dir="${a%:/work/perl_out}" ;; esac
        [[ "$prev" == -o ]] && out_name="${a#/work/perl_out/}"
        prev="$a"
    done
    printf 'v\tI:1\nv\tI:1\n' >"$out_dir/$out_name"
    ;;
esac
exit 0
STUB
    chmod +x "$_gt/bin/vep" "$_gt/bin/docker"
    _gt_write_table() { # _gt_write_table <path> <rows...>: the header, then each argument as its lines
        local path="$1"
        shift
        {
            printf '%s\n' "$_gt_table_header"
            [[ $# -eq 0 ]] || printf '%s\n' "$@"
        } >"$path"
    }
    # _gt_run <name> <script> <table> <harness args...>: one run, its stub log at $_gt/<name>.argv,
    # its own log at $_gt/<name>.log and its output directory at $_gt/<name>.
    _gt_run() {
        local name="$1" script="$2" table="$3"
        shift 3
        : >"$_gt/$name.argv"
        STUB_ARGV="$_gt/$name.argv" STUB_DIGEST="$_gt_digest" VEP_HARNESS_GENOMES="$table" PATH="$_gt/bin:$PATH" \
            bash "$script" --data-dir "$_gt_data" --output-dir "$_gt/$name" "$@" >"$_gt/$name.log" 2>&1
    }
    _gt_write_table "$_gt/made_up.tsv" "$_gt_human_rows" "$_gt_made_up_row"

    # vep-rs: both the SNP/indel cell and the SV cell read the row's FASTA and JSON cache.
    _gt_run row "$_gt_copy" "$_gt/made_up.tsv" --engine vep-rs --vep-binary "$_gt/bin/vep" --datasets "s90 s91" --walltime-only
    rc=$?
    _n=$(grep -c . "$_gt/row.argv") || _n=0
    [[ "$rc" -eq 0 && "$_n" -eq 4 ]] &&
        ok "a made-up genome with a table row runs its two cells (warmup and timed run each)" ||
        bad "made-up genome with a row: rc=$rc invocations=$_n $(grep -m1 ERROR "$_gt/row.log")"
    _n=$(grep -c -- " --offline --json_cache $_gt_data/caches/vep-rs/asm1 --fasta $_gt_data/reference/asm1/made_up.fa --fork " "$_gt/row.argv") || _n=0
    [[ "$_n" -eq 4 ]] &&
        ok "every vep-rs invocation of its cells reads the row's json_cache and fasta, as paths under --data-dir (4 of 4)" ||
        bad "$_n of 4 vep-rs invocations carry the row's json_cache and fasta: $(head -1 "$_gt/row.argv")"
    _n=$(grep -c -- ' --species made_up_species --assembly ASM1 --quiet$' "$_gt/row.argv") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "the SV invocations carry the dataset row's species and assembly" ||
        bad "$_n of 2 SV invocations carry --species made_up_species --assembly ASM1"

    # No row: both cells fail before the engine starts, and the refusal names the genome and
    # the table that lacks it.
    _gt_run norow "$_gt_copy" "$_gt_committed" --engine vep-rs --vep-binary "$_gt/bin/vep" --datasets "s90 s91" --walltime-only
    _n=$(grep -c . "$_gt/norow.argv") || _n=0
    _rows=$(grep -c -e ',genome_resources_missing,,1$' "$_gt/norow/wall_times.csv") || _rows=0
    [[ "$_n" -eq 0 && "$_rows" -eq 2 ]] &&
        ok "a genome with no row fails both cells as genome_resources_missing before any engine runs" ||
        bad "genome with no row: $_n engine invocation(s), $_rows genome_resources_missing row(s); expected 0 and 2"
    _n=$(grep -c -x -F "ERROR: [run_clone_measurement] no FASTA is named for the genome made_up_species:ASM1; add its row to $_gt_committed" "$_gt/norow.log") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "each refusal names the genome and the table that has no row for it" ||
        bad "$_n of 2 refusals name made_up_species:ASM1 and $_gt_committed: $(grep -m1 ERROR "$_gt/norow.log")"
    # A row whose json_cache field is empty names a FASTA and no vep-rs cache.
    _gt_write_table "$_gt/nocache.tsv" "$_gt_human_rows" 'made_up_species	ASM1	63	reference/asm1/made_up.fa	'
    _gt_run nocache "$_gt_copy" "$_gt/nocache.tsv" --engine vep-rs --vep-binary "$_gt/bin/vep" --datasets s90 --walltime-only
    _n=$(grep -c . "$_gt/nocache.argv") || _n=0
    [[ "$_n" -eq 0 ]] && grep -q -x -F "ERROR: [run_clone_measurement] no vep-rs JSON cache is named for the genome made_up_species:ASM1; add its row to $_gt/nocache.tsv" "$_gt/nocache.log" &&
        ok "a row with an empty json_cache field names no vep-rs cache: the cell fails before the engine runs" ||
        bad "row with an empty json_cache: $_n invocation(s), $(grep -m1 ERROR "$_gt/nocache.log")"

    # The human rows, and what replaces them: a row below them, then the VEP_* variables.
    _gt_run human "$HARNESS" "$_gt_committed" --engine vep-rs --vep-binary "$_gt/bin/vep" --datasets s01 --walltime-only
    _n=$(grep -c -- " --offline --json_cache $_gt_data/caches/vep-rs/grch37 --fasta $_gt_data/reference/grch37/genome.fa --fork " "$_gt/human.argv") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "the committed table gives a human cell <data>/caches/vep-rs/grch37 and <data>/reference/grch37/genome.fa" ||
        bad "$_n of 2 s01 invocations carry the committed GRCh37 row's paths: $(head -1 "$_gt/human.argv")"
    _gt_write_table "$_gt/moved.tsv" "$_gt_human_rows" 'homo_sapiens	GRCh37	116	alt/other.fa	alt/cache'
    _gt_run moved "$HARNESS" "$_gt/moved.tsv" --engine vep-rs --vep-binary "$_gt/bin/vep" --datasets s01 --walltime-only
    _n=$(grep -c -- " --offline --json_cache $_gt_data/alt/cache --fasta $_gt_data/alt/other.fa --fork " "$_gt/moved.argv") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "the last row for a genome replaces the rows above it" ||
        bad "$_n of 2 s01 invocations carry the last GRCh37 row's paths: $(head -1 "$_gt/moved.argv")"
    VEP_FASTA_GRCH37="$_gt/env.fa" VEP_RS_CACHE_GRCH37="$_gt/env-cache" \
        _gt_run env "$HARNESS" "$_gt/moved.tsv" --engine vep-rs --vep-binary "$_gt/bin/vep" --datasets s01 --walltime-only
    _n=$(grep -c -- " --offline --json_cache $_gt/env-cache --fasta $_gt/env.fa --fork " "$_gt/env.argv") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "VEP_FASTA_GRCH37 and VEP_RS_CACHE_GRCH37 replace the GRCh37 row's fasta and json_cache, as given" ||
        bad "$_n of 2 s01 invocations carry the VEP_FASTA_GRCH37 and VEP_RS_CACHE_GRCH37 paths: $(head -1 "$_gt/env.argv")"

    # Ensembl VEP: the SV cell of the made-up genome reads the slice and passes --cache_version
    # under the row's cache_version; VEP_PERL_CACHE_VERSION names one version for every genome.
    _gt_run perl "$_gt_copy" "$_gt/made_up.tsv" --engine perl --perl-image-digest "$_gt_digest" --datasets s91
    rc=$?
    _n=$(grep -c -- " -v $_gt_data/reference/asm1/made_up.fa:/work/ref.fa:ro .* --cache_version 63 --species made_up_species --assembly ASM1 " "$_gt/perl.argv") || _n=0
    [[ "$rc" -eq 0 && "$_n" -eq 2 ]] &&
        ok "the Perl cell of the made-up genome mounts the row's fasta and passes the row's cache_version, 63" ||
        bad "Perl cell of the made-up genome: rc=$rc, $_n of 2 runs carry the row's fasta and --cache_version 63: $(grep -m1 ERROR "$_gt/perl.log")"
    grep -q -F "perl cache check (made_up_species ASM1): 1 chromosome directories; a slice of a species other than homo_sapiens is checked to hold at least one" "$_gt/perl.log" &&
        ok "the cache check of another species says what it checks" ||
        bad "the cache check line of the made-up genome is absent or worded otherwise"
    VEP_PERL_CACHE_VERSION=115 _gt_run perl115 "$_gt_copy" "$_gt/made_up.tsv" --engine perl --perl-image-digest "$_gt_digest" --datasets s91
    _n=$(grep -c -- ' --cache_version 115 --species made_up_species --assembly ASM1 ' "$_gt/perl115.argv") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "VEP_PERL_CACHE_VERSION names the cache version of every genome, whatever its row says" ||
        bad "$_n of 2 Perl runs pass --cache_version 115 under VEP_PERL_CACHE_VERSION=115"
    rm -rf "$_gt_data/caches/perl/vep-cache/made_up_species/63_ASM1"
    _gt_run perlslice "$_gt_copy" "$_gt/made_up.tsv" --engine perl --perl-image-digest "$_gt_digest" --datasets s91
    _n=$(grep -c . "$_gt/perlslice.argv") || _n=0
    [[ "$_n" -eq 0 ]] && grep -q -F "Perl cache slice not found: $_gt_data/caches/perl/vep-cache/made_up_species/63_ASM1" "$_gt/perlslice.log" &&
        ok "the slice the Perl cell checks is the one under the row's cache version" ||
        bad "a missing 63_ASM1 slice: $_n docker run(s), $(grep -m1 ERROR "$_gt/perlslice.log")"
    # A row whose cache_version field is empty names none: the genome reads the default, 116.
    _gt_write_table "$_gt/noversion.tsv" "$_gt_human_rows" 'made_up_species	ASM1		reference/asm1/made_up.fa	caches/vep-rs/asm1'
    mkdir -p "$_gt_data/caches/perl/vep-cache/made_up_species/116_ASM1/I"
    _gt_run perlnoversion "$_gt_copy" "$_gt/noversion.tsv" --engine perl --perl-image-digest "$_gt_digest" --datasets s91
    _n=$(grep -c -- ' --cache_version 116 --species made_up_species --assembly ASM1 ' "$_gt/perlnoversion.argv") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "a row with an empty cache_version field reads the default cache version, 116" ||
        bad "$_n of 2 Perl runs pass --cache_version 116 for a row with an empty cache_version: $(grep -m1 ERROR "$_gt/perlnoversion.log")"
    # The same holds where no row of the table names a cache_version.
    _gt_write_table "$_gt/onerow.tsv" 'made_up_species	ASM1		reference/asm1/made_up.fa	caches/vep-rs/asm1'
    _gt_run perlonerow "$_gt_copy" "$_gt/onerow.tsv" --engine perl --perl-image-digest "$_gt_digest" --datasets s91
    _n=$(grep -c -- ' --cache_version 116 --species made_up_species --assembly ASM1 ' "$_gt/perlonerow.argv") || _n=0
    [[ "$_n" -eq 2 ]] &&
        ok "a table whose one row leaves cache_version empty is read, and its genome reads the default cache version, 116" ||
        bad "$_n of 2 Perl runs pass --cache_version 116 under a one-row table with an empty cache_version: $(grep -m1 ERROR "$_gt/perlonerow.log")"

    # The Ensembl VEP SNP/indel cell reads the same per-genome version: a row that gives GRCh37
    # another cache_version moves s01. bcftools and tabix are stubs, so the cell's canonical
    # copy is its input as given.
    cat >"$_gt/bin/bcftools" <<'STUB'
#!/usr/bin/env bash
# stub bcftools: `view ... -o <out> <src>` copies <src>, the last argument, to <out>
out=""; prev=""; src=""
for a in "$@"; do [[ "$prev" == -o ]] && out="$a"; prev="$a"; src="$a"; done
cp "$src" "$out"
STUB
    printf '#!/usr/bin/env bash\nexit 0\n' >"$_gt/bin/tabix"
    chmod +x "$_gt/bin/bcftools" "$_gt/bin/tabix"
    mkdir -p "$_gt_data/vcf/grch37/all_variants"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t1\n' | gzip -c >"$_gt_data/vcf/grch37/all_variants/clinvar_grch37.vcf.gz"
    for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
        mkdir -p "$_gt_data/caches/perl/vep-cache/homo_sapiens/115_GRCh37/$c"
    done
    _gt_write_table "$_gt/v115.tsv" "$_gt_human_rows" 'homo_sapiens	GRCh37	115	reference/grch37/genome.fa	caches/vep-rs/grch37'
    _gt_run perlp1 "$_gt_copy" "$_gt/v115.tsv" --engine perl --perl-image-digest "$_gt_digest" --datasets s01
    rc=$?
    _n=$(grep -c -- ' --cache_version 115 --species homo_sapiens --assembly GRCh37 ' "$_gt/perlp1.argv") || _n=0
    [[ "$rc" -eq 0 && "$_n" -eq 2 ]] && grep -q -F "perl cache check (GRCh37): all 25 primary contigs present" "$_gt/perlp1.log" &&
        ok "the Perl SNP/indel cell checks the slice and passes --cache_version under its genome's row (115 for GRCh37)" ||
        bad "Perl s01 under a GRCh37 row at 115: rc=$rc, $_n of 2 runs pass --cache_version 115: $(grep -m1 ERROR "$_gt/perlp1.log")"

    # fastVEP is measured on the two human assemblies alone.
    cp "$_gt/bin/vep" "$_gt/bin/fastvep"
    _gt_run fastvep "$_gt_copy" "$_gt/made_up.tsv" --engine fastvep --fastvep-binary "$_gt/bin/fastvep" --datasets "s90 s91" --walltime-only
    _n=$(grep -c . "$_gt/fastvep.argv") || _n=0
    _rows=$(grep -c -x -F "ERROR: [run_clone_measurement] no fastVEP cache is named for the genome made_up_species:ASM1; fastVEP is measured on the two human assemblies alone" "$_gt/fastvep.log") || _rows=0
    [[ "$_n" -eq 0 && "$_rows" -eq 2 ]] &&
        ok "a fastVEP cell of another genome fails before the engine runs and says fastVEP is measured on the two human assemblies alone" ||
        bad "fastVEP cells of the made-up genome: $_n invocation(s), $_rows of 2 refusals"
    _two=$(grep -c -e '^FASTVEP_CACHE_TABLE="homo_sapiens:GRCh37	\$FASTVEP_CACHE_GRCH37$' -e '^homo_sapiens:GRCh38	\$FASTVEP_CACHE_GRCH38"$' -e '^GFF3_TABLE="homo_sapiens:GRCh37	\$GFF3_GRCH37$' -e '^homo_sapiens:GRCh38	\$GFF3_GRCH38"$' "$HARNESS") || _two=0
    [[ "$_two" -eq 4 ]] &&
        ok "the fastVEP cache and GFF3 tables are the two human lines each" ||
        bad "the fastVEP cache and GFF3 tables hold $_two of their 4 human lines"

    # A table the reader refuses stops the run before the output directory is made.
    _gt_refused() { # _gt_refused <name> <table> <expected last line>
        local name="$1" table="$2" expected="$3" rc
        _gt_run "refused" "$HARNESS" "$table" --engine vep-rs --vep-binary "$_gt/bin/vep" --datasets s01 --walltime-only
        rc=$?
        if [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_gt/refused.log")" -eq 1 && "$(tail -1 "$_gt/refused.log")" == "$expected" && ! -e "$_gt/refused" && ! -s "$_gt/refused.argv" ]]; then
            ok "$name"
        else
            bad "$name: rc=$rc, output directory $([[ -e "$_gt/refused" ]] && echo made || echo absent), last line '$(tail -1 "$_gt/refused.log")'"
        fi
    }
    _gt_refused "a table that is no file stops the run before the output directory is made" "$_gt/absent.tsv" \
        "ERROR: [run_clone_measurement] the genome table cannot be read: $_gt/absent.tsv"
    printf 'genome\tfasta\nhomo_sapiens\tgenome.fa\n' >"$_gt/nokey.tsv"
    _gt_refused "a table with no species and assembly columns stops the run" "$_gt/nokey.tsv" \
        "ERROR: [run_clone_measurement] the genome table names no genome under a species and an assembly column: $_gt/nokey.tsv"
    _gt_write_table "$_gt/short.tsv" "$_gt_human_rows" 'made_up_species	ASM1	63	reference/asm1/made_up.fa'
    _gt_refused "a short row is refused with the reader's own line" "$_gt/short.tsv" \
        "ERROR: [harness_genomes] $_gt/short.tsv: line 4 has 4 fields where the header has 5"
    printf 'species\tassembly\tcache_version\tfasta\nhomo_sapiens\tGRCh37\t116\treference/grch37/genome.fa\n' >"$_gt/nojson.tsv"
    _gt_refused "a table without a json_cache column is refused with the reader's own line" "$_gt/nojson.tsv" \
        "ERROR: [harness_genomes] $_gt/nojson.tsv: no column json_cache; the columns are species, assembly, cache_version, fasta"
    _gt_write_table "$_gt/badversion.tsv" "$_gt_human_rows" 'made_up_species	ASM1	sixty-three	reference/asm1/made_up.fa	caches/vep-rs/asm1'
    _gt_refused "a cache_version that is no positive integer stops the run and names its genome" "$_gt/badversion.tsv" \
        "ERROR: [run_clone_measurement] the genome table gives made_up_species:ASM1 the cache_version 'sixty-three', which is no positive integer: $_gt/badversion.tsv"
    # The same genome above 200 rows whose species names run to 2,000 characters, so that
    # the cache versions are more than a pipe holds (65,536 bytes on Linux) and more than
    # gawk or mawk reads ahead of its first line: the check leaves at that genome's line,
    # and a writer left behind would be ended by the closed pipe, and the run with it,
    # before the refusal is printed.
    _gt_long_rows=$(awk 'BEGIN { pad = sprintf("%2000s", ""); gsub(/ /, "x", pad); for (i = 1; i <= 200; i++) printf "made_up_species_%s_%03d\tASM1\t63\treference/asm1/made_up.fa\tcaches/vep-rs/asm1\n", pad, i }')
    _gt_write_table "$_gt/badversion_long.tsv" 'made_up_species	ASM1	sixty-three	reference/asm1/made_up.fa	caches/vep-rs/asm1' "$_gt_long_rows"
    _gt_refused "the same refusal is printed above rows whose cache versions are more than a pipe holds" "$_gt/badversion_long.tsv" \
        "ERROR: [run_clone_measurement] the genome table gives made_up_species:ASM1 the cache_version 'sixty-three', which is no positive integer: $_gt/badversion_long.tsv"
    # The harness's own text. A writer of several lines is ended at random as well, between
    # two of its writes, so the form is held on every line, whatever the size of its value:
    # no line whose status is a pipeline's (an assignment from it, the pipeline as a
    # statement, a one-line function's body) pipes printf or echo into an awk with a bare
    # exit, a head, a grep -m or a sed that quits.
    _gt_status_of_early_reader='^[[:space:]]*([A-Za-z_][A-Za-z0-9_]*=("[^"]*)?[$][(]|([A-Za-z_][A-Za-z0-9_]*[(][)][[:space:]]*[{][[:space:]]*)?)(printf|echo)[[:space:]]([^|$]|[$][^(])*[|][[:space:]]*(awk[[:space:]].*exit[[:space:]]*[;}]|head([[:space:]]|$)|[ef]?grep[[:space:]][^|]*-m[[:space:]]*[0-9]|sed[[:space:]][^|]*[0-9]q)'
    _gt_early=$(grep -nE "$_gt_status_of_early_reader" "$HARNESS" | cut -d: -f1 | paste -sd, -)
    [[ -z "$_gt_early" ]] &&
        ok "no line of the harness reads the status of printf or echo piped into a reader that leaves before the end" ||
        bad "printf or echo piped into a reader that leaves before the end, its status read: run_clone_measurement.sh line(s) $_gt_early"
fi

echo
echo "FASTVEP-CACHE: a fastVEP cell whose transcript cache path names no file with bytes is refused before fastVEP starts"

# fastVEP, given a --transcript-cache path that names no file, parses the --gff3 file, saves a
# cache of its own at the path and exits 0. The stub binary and the inputs are
# FASTVEP-EXTRA-ARGS's; each tree below is a copy of that data tree that differs in its cache.
if [[ -f "$HARNESS" ]]; then
    _fc="$TMP/fastvep_cache"
    _fc_cells=$(sed -n '/^measure_p1() {/,/^}/p; /^measure_sv() {/,/^}/p' "$HARNESS")
    _fc_calls=$(grep -c -F -e '[[ -f "$cache" && -s "$cache" ]] || {' <<<"$_fc_cells") || _fc_calls=0
    _fc_rows=$(grep -A 2 -F -e '[[ -f "$cache" && -s "$cache" ]] || {' <<<"$_fc_cells" | grep -c -e '"fastvep_cache_missing" "" "1"$') || _fc_rows=0
    [[ "$_fc_calls" -eq 2 && "$_fc_rows" -eq 2 ]] &&
        ok "both fastVEP cells (SNP/indel, SV) test the cache file and record fastvep_cache_missing on its refusal (2 of 2)" ||
        bad "$_fc_calls fastVEP cell(s) test the cache file and $_fc_rows record fastvep_cache_missing on its refusal (expected 2 and 2)"
    fastvep_cache_tree() { # fastvep_cache_tree <name>: copies the data tree and prints its GRCh37 cache path
        mkdir -p "$_fc/$1"
        cp -R "$_fa/data" "$_fc/$1/data"
        echo "$_fc/$1/data/caches/fastvep/grch37.bin"
    }
    run_fastvep_cache() { # run_fastvep_cache <name>: both GRCh37 cells on the tree, into <name>/out
        : >"$_fc/$1.argv"
        STUB_ARGV="$_fc/$1.argv" bash "$HARNESS" --engine fastvep --fastvep-binary "$_fa/bin/fastvep" --data-dir "$_fc/$1/data" \
            --output-dir "$_fc/$1/out" --datasets "s01 s07" --walltime-only >"$_fc/$1.log" 2>&1
    }
    fastvep_cache_refused() { # fastvep_cache_refused <name> <label>, once the tree's cache path is set up
        local o="$1" label="$2" cache="$_fc/$1/data/caches/fastvep/grch37.bin" rc n_inv rows refused n_err
        run_fastvep_cache "$o"; rc=$?
        n_inv=$(grep -c . "$_fc/$o.argv") || n_inv=0
        rows=$(awk -F, 'NR>1' "$_fc/$o/out/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        refused=$(awk -F, 'NR>1 && $1=="fastvep" && ($5=="s01" || $5=="s07") && $8=="fastvep_cache_missing" && $9=="" && $10=="1"' "$_fc/$o/out/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        n_err=$(grep -c -F -e "ERROR: [run_clone_measurement] s01: no fastVEP transcript cache file with bytes at $cache;" -e "ERROR: [run_clone_measurement] s07: no fastVEP transcript cache file with bytes at $cache;" "$_fc/$o.log") || n_err=0
        if [[ "$rc" -eq 0 && "$n_inv" -eq 0 && "$rows" -eq 2 && "$refused" -eq 2 && "$n_err" -eq 2 ]] &&
            grep -q -F -e "ERROR: [run_clone_measurement] 2 dataset(s) failed: s01 s07" "$_fc/$o.log"; then
            ok "$label: both cells are refused by name before fastVEP starts, each with one row of run kind fastvep_cache_missing and exit code 1"
        else
            bad "$label: rc=$rc, $n_inv fastVEP start(s), $rows row(s) of which $refused refused, $n_err ERROR line(s) naming the path: $(grep -e ERROR "$_fc/$o.log" | head -2 | tr '\n' '|' | cut -c1-300)"
        fi
    }
    if /usr/bin/time -v true >/dev/null 2>&1; then
        _cache=$(fastvep_cache_tree absent)
        rm -rf "$(dirname "$_cache")"
        fastvep_cache_refused absent "no file at the cache path"
        _cache=$(fastvep_cache_tree empty)
        : >"$_cache"
        fastvep_cache_refused empty "an empty file at the cache path"
        _cache=$(fastvep_cache_tree dangling)
        rm -f "$_cache"
        ln -s "$_fc/dangling/nowhere.bin" "$_cache"
        fastvep_cache_refused dangling "a link at the cache path that names no file"
        # The control: the same tree with its cache file runs every cell on that file.
        _cache=$(fastvep_cache_tree present)
        run_fastvep_cache present; rc=$?
        _n_inv=$(grep -c -F -e " --transcript-cache $_cache " "$_fc/present.argv") || _n_inv=0
        _rows=$(awk -F, 'NR>1' "$_fc/present/out/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        _bad_rows=$(awk -F, 'NR>1 && ($1!="fastvep" || $10!="0" || ($8!="warmup" && $8!="timed"))' "$_fc/present/out/wall_times.csv" 2>/dev/null | wc -l | tr -d ' ')
        [[ "$rc" -eq 0 && "$_n_inv" -eq 6 && "$_rows" -eq 4 && "$_bad_rows" -eq 0 ]] &&
            ok "the same tree with its cache file starts all six fastVEP annotate runs on that file and records four rows with exit 0" ||
            bad "a tree with its cache file: rc=$rc, $_n_inv of 6 starts on the file, $_rows row(s), $_bad_rows of them failed: $(grep -e ERROR "$_fc/present.log" | head -2 | tr '\n' '|' | cut -c1-300)"
    else
        echo "  skip functional fastVEP cache runs: /usr/bin/time -v is not GNU time here"
    fi
fi

echo
echo "FILE-CENSUS: the Perl cells' passes over whole files run the tree's release build of tools/vep-file-census"

# The Ensembl VEP cells strip the contig prefix of a raw input and count the rows and records of
# whole outputs and inputs. Those passes run tools/vep-file-census: the block that finds
# it and the two functions that call it are extracted and run against stub trees and a
# stub tool that records its command lines, the harness is run whole in a tree without the
# build, and both functions run once on this tree's own release build, a check that is
# skipped where the tool is not built.
HARNESS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_clone_measurement.sh"
if [[ -f "$HARNESS" ]]; then
    _fc="$TMP/file_census"
    _census_rel="tools/vep-file-census/target/release/vep-file-census"
    mkdir -p "$_fc/bin" "$_fc/data" "$_fc/repo-built/$(dirname "$_census_rel")" "$_fc/repo-unbuilt" "$_fc/repo-inert/$(dirname "$_census_rel")" "$_fc/repo-directory/$_census_rel"
    : >"$_fc/repo-inert/$_census_rel"
    # The stub tool: one line a start, then STUB_COUNT as a count where one is set; a
    # count of 0 leaves with status 1, as the tool does.
    cat >"$_fc/repo-built/$_census_rel" <<STUB
#!/usr/bin/env bash
printf '%s\\n' "\$*" >>"$_fc/started"
if [[ "\$1" == strip-chr ]]; then printf 'stripped\\n'; exit 0; fi
printf '%s\\n' "\${STUB_COUNT:-0}"
[[ "\${STUB_COUNT:-0}" -gt 0 ]]
STUB
    chmod +x "$_fc/repo-built/$_census_rel"

    _census_block=$(awk '
        !p && /^FILE_CENSUS=\(\)$/ { p = 1 }
        p { print }
        p && /^fi$/ { exit }
    ' "$HARNESS")
    [[ "$_census_block" == *'if [[ "$ENGINE" == "perl" ]]; then'* && "$_census_block" == *"census_bin=\"\$REPO_DIR/$_census_rel\""* && "$_census_block" == *$'\nfi' ]] &&
        ok "the file census block is extracted whole: the tree's release build its one candidate, asked for by the Perl engine alone" ||
        bad "could not extract the file census block of the harness"
    _block_line=$(grep -n -e '^FILE_CENSUS=()$' "$HARNESS" | head -1 | cut -d: -f1)
    _loop_line=$(grep -n -e '^for entry in "\${DATASETS\[@\]}"; do$' "$HARNESS" | head -1 | cut -d: -f1)
    [[ -n "$_block_line" && -n "$_loop_line" && "$_block_line" -lt "$_loop_line" ]] &&
        ok "the file census is found before the dataset loop (line $_block_line before $_loop_line)" ||
        bad "the file census block must stand before the dataset loop (block at ${_block_line:-?}, loop at ${_loop_line:-?})"
    # $1 is the engine, $2 a tree: runs the block; prints stderr, the command found, then `rc=`.
    _census_found() {
        # The extracted block reads the names this subshell sets.
        # shellcheck disable=SC2034
        (
            ENGINE="$1"
            REPO_DIR="$_fc/$2"
            eval "$_census_block"
            printf 'found=%s\n' "${FILE_CENSUS[*]:-}"
        ) 2>&1
        echo "rc=$?"
    }
    _census_refusal() {
        printf '%s' "ERROR: [run_clone_measurement] the file census $_fc/$1/$_census_rel is not an executable: build the tool with 'cargo build --release --locked --manifest-path $_fc/$1/tools/vep-file-census/Cargo.toml'"
    }
    [[ "$(_census_found perl repo-built)" == "found=$_fc/repo-built/$_census_rel"$'\n'"rc=0" ]] &&
        ok "a Perl run: the tree's release build of tools/vep-file-census is the command" ||
        bad "a Perl run in a built tree did not find the file census: $(_census_found perl repo-built)"
    for _engine in vep-rs fastvep; do
        [[ "$(_census_found "$_engine" repo-unbuilt)" == "found="$'\n'"rc=0" ]] &&
            ok "a $_engine run asks for no file census: a tree without the build is not refused" ||
            bad "a $_engine run was refused for the file census: $(_census_found "$_engine" repo-unbuilt)"
    done
    for _shape in unbuilt inert directory; do
        [[ "$(_census_found perl "repo-$_shape")" == "$(_census_refusal "repo-$_shape")"$'\n'"rc=1" ]] &&
            ok "a Perl run in a tree whose build is $_shape: one line names the build command, exit 1" ||
            bad "a Perl run in a tree whose build is $_shape was not refused: $(_census_found perl "repo-$_shape")"
    done
    # The whole harness in a tree without the build: the refusal before any dataset starts.
    mkdir -p "$_fc/repo-unbuilt/scripts/concordance"
    cp "$HARNESS" "$(dirname "$HARNESS")/harness_genomes.sh" "$(dirname "$HARNESS")/harness_genomes.tsv" "$_fc/repo-unbuilt/scripts/concordance/"
    : >"$_dd/calls"
    DIGEST_CALLS="$_dd/calls" STUB_DIGEST="$PINNED_DIGEST" PATH="$_dd/bin:$PATH" \
        bash "$_fc/repo-unbuilt/scripts/concordance/run_clone_measurement.sh" --engine perl --data-dir "$_fc/data" --output-dir "$_fc/out" --datasets s07 >"$_fc/whole.log" 2>&1
    rc=$?
    [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_fc/whole.log")" -eq 1 && "$(tail -1 "$_fc/whole.log")" == "$(_census_refusal repo-unbuilt)" ]] &&
        ! grep -q "dataset s07 (SV_GRCh37 homo_sapiens GRCh37 sv) starting" "$_fc/whole.log" && ! grep -q "^docker run" "$_dd/calls" &&
        ok "a Perl harness run in a tree without the build: exit 1 on the one refusal line, before any dataset starts and before any docker run" ||
        bad "a Perl harness run without the file census: expected rc=1 and the refusal as the one ERROR and last line, got rc=$rc '$(tail -2 "$_fc/whole.log" | tr '\n' ' ')'"

    # assert_perl_annotated, extracted: which files it counts, and with which rule.
    _annotated=$(sed -n '/^assert_perl_annotated() {/,/^}/p' "$HARNESS")
    mkdir -p "$_fc/perl_out" "$_fc/in"
    : >"$_fc/perl_out/a.txt"; : >"$_fc/perl_out/b.txt"; : >"$_fc/perl_out/b.txt_warnings.txt"
    : >"$_fc/in/x.vcf"; : >"$_fc/in/y.vcf.gz"; : >"$_fc/in/notes.txt"
    # $1 is the count the stub prints; the rest are the function's arguments.
    _annotated_run() {
        local count="$1"
        shift
        : >"$_fc/started"
        (
            FILE_CENSUS=("$_fc/repo-built/$_census_rel")
            export STUB_COUNT="$count"
            eval "$_annotated"
            assert_perl_annotated "$@"
        ) >"$_fc/annotated.log" 2>&1
    }
    _annotated_run 3 "$_fc/perl_out" "$_fc/in" s07
    rc=$?
    _want="count --not-starting # $_fc/perl_out/a.txt
count --not-starting # $_fc/perl_out/b.txt
count --not-starting # $_fc/in/x.vcf
count --gunzip --not-starting # $_fc/in/y.vcf.gz"
    [[ "$rc" -eq 0 && "$(cat "$_fc/started")" == "$_want" ]] && grep -q 's07 perl annotation check: 6 rows / 6 records$' "$_fc/annotated.log" &&
        ok "a directory of outputs and of inputs: one count a file, the warnings companions left out, a gzipped input inflated by the tool, the counts summed" ||
        bad "assert_perl_annotated over directories: expected four counts and 6 rows / 6 records, got rc=$rc '$(tr '\n' '|' <"$_fc/started")' and '$(tail -1 "$_fc/annotated.log")'"
    _annotated_run 5 "$_fc/perl_out/a.txt" "$_fc/in/y.vcf.gz" s01
    rc=$?
    [[ "$rc" -eq 0 && "$(cat "$_fc/started")" == "count --not-starting # $_fc/perl_out/a.txt"$'\n'"count --gunzip --not-starting # $_fc/in/y.vcf.gz" ]] && grep -q 's01 perl annotation check: 5 rows / 5 records$' "$_fc/annotated.log" &&
        ok "one output and one gzipped input: two counts, by the same rule" ||
        bad "assert_perl_annotated over one file each: got rc=$rc '$(tr '\n' '|' <"$_fc/started")' and '$(tail -1 "$_fc/annotated.log")'"
    _annotated_run 0 "$_fc/perl_out/a.txt" "$_fc/in/y.vcf.gz" s01
    rc=$?
    [[ "$rc" -eq 1 ]] && grep -q 's01 perl annotation check: 0 rows / 0 records$' "$_fc/annotated.log" && grep -q 'ERROR: \[run_clone_measurement\] s01 counted 0 input records' "$_fc/annotated.log" &&
        ok "a count of 0, which the tool prints with status 1, is read as 0 and fails the gate" ||
        bad "a zero count did not fail the gate as 0 records (rc=$rc): $(tail -2 "$_fc/annotated.log" | tr '\n' ' ')"

    # canonicalize, extracted, with stub bcftools, tabix and bgzip: the prefix pass is the tool's.
    _canonicalize=$(sed -n '/^CANONICAL_CONTIGS_LIST=/,/^}/p' "$HARNESS")
    printf '#!/usr/bin/env bash\nout=""; prev=""; for a in "$@"; do [[ "$prev" == -o ]] && out="$a"; prev="$a"; done\ncp "${@: -1}" "$out"\n' >"$_fc/bin/bcftools"
    printf '#!/usr/bin/env bash\nexit 0\n' >"$_fc/bin/tabix"
    printf '#!/usr/bin/env bash\nexec gzip -c\n' >"$_fc/bin/bgzip"
    chmod +x "$_fc/bin/bcftools" "$_fc/bin/tabix" "$_fc/bin/bgzip"
    printf '##contig=<ID=chr21>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr21\t1000\tv1\tA\tG\t.\t.\tAF=0.5\nchr21\t2000\tv2\tC\tT\t.\t.\t.\n' | gzip -c >"$_fc/prefixed.vcf.gz"
    printf '##contig=<ID=21>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1000\tv1\tA\tG\t.\t.\t.\n' | gzip -c >"$_fc/plain.vcf.gz"
    # $1 is the tool, $2 the raw input, $3 the canonical output.
    _canonicalize_run() {
        : >"$_fc/started"
        (
            PATH="$_fc/bin:$PATH"
            FILE_CENSUS=("$1")
            eval "$_canonicalize"
            canonicalize "$2" "$3"
        ) >"$_fc/canonicalize.log" 2>&1
    }
    _canonicalize_run "$_fc/repo-built/$_census_rel" "$_fc/prefixed.vcf.gz" "$_fc/canon1.vcf.gz"
    rc=$?
    [[ "$rc" -eq 0 && "$(cat "$_fc/started")" == "strip-chr $_fc/prefixed.vcf.gz" && "$(gzip -dc "$_fc/canon1.vcf.gz")" == "stripped" ]] &&
        ok "an input on chr-prefixed contigs: one strip-chr start of the tool on the raw input, its output block-gzipped and filtered" ||
        bad "canonicalize of a prefixed input: expected one strip-chr start, got rc=$rc '$(tr '\n' '|' <"$_fc/started")'"
    _canonicalize_run "$_fc/repo-built/$_census_rel" "$_fc/plain.vcf.gz" "$_fc/canon2.vcf.gz"
    rc=$?
    [[ "$rc" -eq 0 && ! -s "$_fc/started" ]] && cmp -s "$_fc/plain.vcf.gz" "$_fc/canon2.vcf.gz" &&
        ok "an input on unprefixed contigs: the tool is not started and the raw input is filtered as it is" ||
        bad "canonicalize of an unprefixed input started the tool or changed the input (rc=$rc): '$(tr '\n' '|' <"$_fc/started")'"
    # A prefix pass that stops part-way has written the lines before the stop, which the
    # bcftools stub would copy as a whole input: the pass's own status must fail the cell.
    printf '#!/usr/bin/env bash\nprintf "##contig=<ID=21>\\n"\nexit 1\n' >"$_fc/bin/stopping-census"
    chmod +x "$_fc/bin/stopping-census"
    _canonicalize_run "$_fc/bin/stopping-census" "$_fc/prefixed.vcf.gz" "$_fc/canon_stopped.vcf.gz"
    rc=$?
    [[ "$rc" -eq 1 && ! -e "$_fc/canon_stopped.vcf.gz" && ! -e "$_fc/canon_stopped.stripped.vcf.gz" && "$(cat "$_fc/canonicalize.log")" == "ERROR: [run_clone_measurement] contig-prefix pass failed for $_fc/prefixed.vcf.gz" ]] &&
        ok "a prefix pass that stops part-way: canonicalize fails on one line naming the input, and neither its partial copy nor a canonical input is left" ||
        bad "canonicalize over a prefix pass that fails: expected rc=1, the one line and no file left, got rc=$rc '$(tr '\n' '|' <"$_fc/canonicalize.log")' and $(find "$_fc" -maxdepth 1 -name 'canon_stopped*' | wc -l | tr -d ' ') file(s)"

    # Both functions on this tree's own release build.
    _real_census="$(cd "$(dirname "$HARNESS")/../.." && pwd)/$_census_rel"
    if [[ -x "$_real_census" ]]; then
        _canonicalize_run "$_real_census" "$_fc/prefixed.vcf.gz" "$_fc/canon3.vcf.gz"
        rc=$?
        [[ "$rc" -eq 0 && "$(gzip -dc "$_fc/canon3.vcf.gz")" == "$(printf '##contig=<ID=21>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1000\tv1\tA\tG\t.\t.\tAF=0.5\n21\t2000\tv2\tC\tT\t.\t.\t.')" ]] &&
            ok "canonicalize on the tree's own build: the contig line and both records lose their chr prefix" ||
            bad "canonicalize on the tree's own build wrote '$(gzip -dc "$_fc/canon3.vcf.gz" | tr '\n' '|')' (rc=$rc)"
        printf '## header\n#Uploaded_variation\tLocation\nv1\t21:1000\nv1\t21:1000\nv2\t21:2000\n' >"$_fc/real_out.txt"
        (
            FILE_CENSUS=("$_real_census")
            eval "$_annotated"
            assert_perl_annotated "$_fc/real_out.txt" "$_fc/prefixed.vcf.gz" s01
        ) >"$_fc/real1.log" 2>&1
        rc=$?
        [[ "$rc" -eq 0 ]] && grep -q 's01 perl annotation check: 3 rows / 2 records$' "$_fc/real1.log" &&
            ok "assert_perl_annotated on the tree's own build: three rows for two records pass the floor" ||
            bad "assert_perl_annotated on the tree's own build: expected 3 rows / 2 records and rc=0, got rc=$rc '$(tail -1 "$_fc/real1.log")'"
        printf 'v1\t21:1000\n' >"$_fc/real_short.txt"
        (
            FILE_CENSUS=("$_real_census")
            eval "$_annotated"
            assert_perl_annotated "$_fc/real_short.txt" "$_fc/prefixed.vcf.gz" s01
        ) >"$_fc/real2.log" 2>&1
        rc=$?
        [[ "$rc" -eq 1 ]] && grep -q 'ERROR: \[run_clone_measurement\] s01 Perl emitted 1 rows for 2 input records' "$_fc/real2.log" &&
            ok "assert_perl_annotated on the tree's own build: one row for two records fails the floor" ||
            bad "assert_perl_annotated on the tree's own build did not fail one row for two records (rc=$rc): $(tail -2 "$_fc/real2.log" | tr '\n' ' ')"
    else
        echo "  skip the two functions on the tree's own build: no release build of tools/vep-file-census"
    fi
else
    bad "run_clone_measurement.sh not found next to this test"
fi

echo
echo "SV-COMPARATOR: a structural-variant cell is scored by the release build of tools/vep-sv-compare, the one under CARGO_TARGET_DIR before the one under the tool's own target directory, and a run that finds neither is refused before any engine runs"

# The harness is copied into two trees of its own, so REPO_DIR is the copy and the tools it
# finds are the stubs placed there: a tuple comparator and a field comparator in both, which
# a vep-rs run that scores F1 resolves first, and in one the structural-variant tool. Three
# directories stand in for CARGO_TARGET_DIR: one holds a build of the tool, one holds
# nothing, and one holds a file at the build's path that is not executable. Each stub of the
# tool appends its own path and its arguments, one line each. The stub engine and the data
# tree (two inputs under the GRCh37 SV directory) are EXTRA-ARGS's.
if [[ -f "$HARNESS" ]]; then
    _sc="$TMP/sv_comparator"
    _sc_tool="tools/vep-sv-compare/target/release/vep-sv-compare"
    _sc_target="$_sc/cargo-target"
    _sv_cargo_target=""
    for _tree in built unbuilt; do
        mkdir -p "$_sc/$_tree/scripts/concordance"
        cp "$HARNESS" "$(dirname "$HARNESS")/harness_genomes.sh" "$(dirname "$HARNESS")/harness_genomes.tsv" "$_sc/$_tree/scripts/concordance/"
        for _stub in vep-outputs-compare vep-fields-compare; do
            mkdir -p "$_sc/$_tree/tools/$_stub/target/release"
            printf '#!/bin/bash\n' >"$_sc/$_tree/tools/$_stub/target/release/$_stub"
            chmod +x "$_sc/$_tree/tools/$_stub/target/release/$_stub"
        done
    done
    mkdir -p "$(dirname "$_sc/built/$_sc_tool")" "$_sc_target/release" "$_sc/cargo-target-unbuilt" "$_sc/cargo-target-inert/release"
    for _stub in "$_sc/built/$_sc_tool" "$_sc_target/release/vep-sv-compare"; do
        printf '#!/bin/bash\nprintf "%%s\\n" "$0" "$@" >>"%s"\n' "$_sc/tool.argv" >"$_stub"
        chmod +x "$_stub"
    done
    : >"$_sc/cargo-target-inert/release/vep-sv-compare"
    # run_sv_comparator_harness <tree> <run name> [harness args...]: the run's CARGO_TARGET_DIR
    # is $_sv_cargo_target where the call sets it and unset elsewhere, whatever the shell that
    # runs this file exports.
    run_sv_comparator_harness() {
        local tree="$1" run_name="$2"
        shift 2
        (
            unset VEP_OUTPUTS_COMPARATOR CARGO_TARGET_DIR
            [[ -z "$_sv_cargo_target" ]] || export CARGO_TARGET_DIR="$_sv_cargo_target"
            STUB_ARGV="$_sc/$run_name.argv" bash "$_sc/$tree/scripts/concordance/run_clone_measurement.sh" --engine vep-rs --vep-binary "$_ea/bin/vep" \
                --data-dir "$_ea/data" --output-dir "$_sc/$run_name" --reference-release 116.2 "$@"
        ) >"$_sc/$run_name.log" 2>&1
    }
    # _sv_scored_by <tree> <run name> <binary>: a run that scores s07 in <tree> exits 0, records
    # no failed comparison, and starts <binary> and no other build of the tool, once, with the
    # cell's directories, cache and release. _sv_seen is what the run left, for a failing check.
    _sv_scored_by() {
        local tree="$1" run_name="$2" binary="$3" rc expected started
        rm -f "$_sc/tool.argv"
        run_sv_comparator_harness "$tree" "$run_name" --datasets s07
        rc=$?
        expected=$(printf '%s\n' "$binary" --assembly grch37 --input-dir "$_sv_gt/canonical_inputs/grch37" --perl-dir "$_sv_gt/grch37" \
            --rust-dir "$_sc/$run_name/s07/rust" --vep-rs-cache "$_ea/data/caches/vep-rs/grch37" --reference-release 116.2 --output-dir "$_sc/$run_name/s07/report")
        started=$(cat "$_sc/tool.argv" 2>/dev/null)
        _sv_seen="rc=$rc; the tool's argument vector: $(printf '%s' "$started" | tr '\n' ' ' | cut -c1-300) $(tail -2 "$_sc/$run_name.log" | tr '\n' '|' | cut -c1-300)"
        [[ "$rc" -eq 0 && "$started" == "$expected" ]] && ! grep -q 'compare_failed' "$_sc/$run_name/wall_times.csv" 2>/dev/null
    }
    if /usr/bin/time -v true >/dev/null 2>&1; then
        _sv_gt="$_ea/data/ground_truth/perl/sv_per_vcf"
        _sv_calls=$(grep -c -e '^[[:space:]]*"\$SV_COMPARE" \\$' "$HARNESS") || _sv_calls=0
        if _sv_scored_by built scored "$_sc/built/$_sc_tool" && [[ "$_sv_calls" -eq 2 ]]; then
            ok "the structural-variant cell runs the tree's release build of tools/vep-sv-compare with the cell's directories, cache and release, and both call sites (vep-rs, fastvep) run that command"
        else
            bad "$_sv_calls call site(s) (expected 2); $_sv_seen"
        fi
        _sv_cargo_target="$_sc_target" _sv_scored_by unbuilt target-only "$_sc_target/release/vep-sv-compare" &&
            ok "a tree whose only build of the tool stands under CARGO_TARGET_DIR: the structural-variant cell is scored by that build" ||
            bad "a build under CARGO_TARGET_DIR alone: $_sv_seen"
        _sv_cargo_target="$_sc_target" _sv_scored_by built target-first "$_sc_target/release/vep-sv-compare" &&
            ok "builds under CARGO_TARGET_DIR and under the tool's own target directory: the one under CARGO_TARGET_DIR is the command, and the other is not started" ||
            bad "both builds present: $_sv_seen"
        _sv_cargo_target="$_sc/cargo-target-unbuilt" _sv_scored_by built target-unbuilt "$_sc/built/$_sc_tool" &&
            _sv_cargo_target="$_sc/cargo-target-inert" _sv_scored_by built target-inert "$_sc/built/$_sc_tool" &&
            ok "CARGO_TARGET_DIR naming a directory that holds no build, or one whose file at the build's path is not executable: the build under the tool's own target directory is the command" ||
            bad "CARGO_TARGET_DIR without an executable build: $_sv_seen"
    else
        echo "  skip functional structural-variant comparator runs: /usr/bin/time -v is not GNU time here"
    fi
    # Without the tool: a run that would score a structural-variant cell stops on one line
    # before it makes its output directory or starts the engine; a run that scores SNP/indel
    # cells alone, and one that scores nothing, ask for no such tool.
    _sv_refusal="ERROR: [run_clone_measurement] the structural-variant comparator $_sc/unbuilt/$_sc_tool is not an executable: build the tool with 'cargo build --release --locked --manifest-path $_sc/unbuilt/tools/vep-sv-compare/Cargo.toml'"
    run_sv_comparator_harness unbuilt refused --datasets "s01 s07"; rc=$?
    run_sv_comparator_harness unbuilt snp-indel --datasets s01
    run_sv_comparator_harness unbuilt walltime --datasets s07 --walltime-only
    if [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_sc/refused.log")" -eq 1 && "$(tail -1 "$_sc/refused.log")" == "$_sv_refusal" && ! -e "$_sc/refused" && ! -s "$_sc/refused.argv" ]] &&
        ! grep -q 'structural-variant comparator' "$_sc/snp-indel.log" "$_sc/walltime.log"; then
        ok "a tree without the tool: a run that scores a structural-variant cell exits 1 on the one refusal line, the engine never started and no output directory made; an SNP/indel run and a --walltime-only run are not refused"
    else
        bad "tree without the tool: rc=$rc, expected 1 and the refusal as the one ERROR and last line, got '$(tail -2 "$_sc/refused.log" | cut -c1-300)'; refusal lines in the SNP/indel and --walltime-only logs: $(cat "$_sc/snp-indel.log" "$_sc/walltime.log" 2>/dev/null | grep -c 'structural-variant comparator')"
    fi
    # With CARGO_TARGET_DIR naming a directory that holds no build either, the refusal is the
    # same line: the path it names is the one under the tool's own target directory.
    _sv_cargo_target="$_sc/cargo-target-unbuilt" run_sv_comparator_harness unbuilt refused-target --datasets "s01 s07"; rc=$?
    if [[ "$rc" -eq 1 && "$(grep -c '^ERROR' "$_sc/refused-target.log")" -eq 1 && "$(tail -1 "$_sc/refused-target.log")" == "$_sv_refusal" && ! -e "$_sc/refused-target" && ! -s "$_sc/refused-target.argv" ]]; then
        ok "no build under CARGO_TARGET_DIR nor under the tool's own target directory: exit 1 on the same one refusal line, which names the path under the tool's own target directory, the engine never started and no output directory made"
    else
        bad "no build in either place: rc=$rc, expected 1 and the refusal as the one ERROR and last line, got '$(tail -2 "$_sc/refused-target.log" | cut -c1-300)'"
    fi
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
