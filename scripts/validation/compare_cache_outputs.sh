#!/usr/bin/env bash
set -euo pipefail

# compare_cache_outputs.sh
#
# Full side-by-side comparison of two JSON caches: runs vep-rs with each cache,
# structurally diffs the caches, compares annotation outputs, and checks against
# hard validation gates.
#
# Usage:
#   scripts/validation/compare_cache_outputs.sh \
#     --baseline-cache /path/to/perl_cache \
#     --builder-cache /path/to/builder_cache \
#     --vcf /path/to/clinvar.vcf.gz \
#     --vep-binary target/release/vep

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VEP_RS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Defaults
BASELINE_CACHE=""
BUILDER_CACHE=""
VCF=""
VEP_BINARY=""
WORK_DIR=""
FORK=1
BUFFER_SIZE=5000
FASTA=""
EXTRA_VEP_ARGS=""

# Hard gates
GATE_F1=0.999
GATE_TOTAL_DELTA=0.001  # 0.1%
GATE_COLOC_DELTA=0.01   # 1%
GATE_FREQ_DELTA=0.01    # 1%
GATE_PER_CSQ_DELTA=0.01 # 1%

# Usage
usage() {
	cat <<'EOF'
Run full side-by-side comparison of baseline vs builder JSON caches.

Usage:
  scripts/validation/compare_cache_outputs.sh --baseline-cache <path> --builder-cache <path> --vcf <path> --vep-binary <path> [options]

Required:
  --baseline-cache <path>   Path to the baseline (Perl-derived) JSON cache
  --builder-cache <path>    Path to the builder-generated JSON cache
  --vcf <path>              Input VCF file (plain or bgzipped)
  --vep-binary <path>       Path to the vep-rs binary

Optional:
  --work-dir <path>         Working directory (default: vep-rs/tmp/cache_comparison)
  --fork <int>              Number of threads for vep-rs (default: 1)
  --buffer-size <int>       Buffer size for vep-rs (default: 5000)
  --fasta <path>            Reference FASTA for HGVS
  --extra-args <args>       Additional vep-rs arguments (quoted string)
  --gate-f1 <float>         Consequence F1 gate (default: 0.999)
  --gate-total-delta <float>  Total annotation delta gate (default: 0.001)
  --gate-coloc-delta <float>  Co-located variant delta gate (default: 0.01)
  --gate-freq-delta <float>   Frequency field delta gate (default: 0.01)
  --gate-per-csq-delta <float>  Per-consequence-type delta gate (default: 0.01)
  --help                    Print this help

Comparator:
  The annotation output comparison runs this tree's release build of the tuple
  comparator, tools/vep-outputs-compare/target/release/vep-outputs-compare
  (cargo build --release --locked --manifest-path tools/vep-outputs-compare/Cargo.toml).
  A build that is not an executable stops the run before vep-rs runs; nothing is
  built here. The structural cache comparison runs the release build of
  tools/vep-cache-diff, found and refused the same way
  (cargo build --release --locked --manifest-path tools/vep-cache-diff/Cargo.toml),
  and the gates are scored by the release build of tools/vep-cache-output-gates
  (cargo build --release --locked --manifest-path tools/vep-cache-output-gates/Cargo.toml).

Outputs:
  <work-dir>/baseline/         Baseline VEP output and metrics
  <work-dir>/builder/          Builder VEP output and metrics
  <work-dir>/cache_diff.json   Structural cache comparison
  <work-dir>/verdict.json      Final PASS/FAIL verdict with metrics
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
	case "$1" in
	--baseline-cache)
		BASELINE_CACHE="$2"
		shift 2
		;;
	--builder-cache)
		BUILDER_CACHE="$2"
		shift 2
		;;
	--vcf)
		VCF="$2"
		shift 2
		;;
	--vep-binary)
		VEP_BINARY="$2"
		shift 2
		;;
	--work-dir)
		WORK_DIR="$2"
		shift 2
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--buffer-size)
		BUFFER_SIZE="$2"
		shift 2
		;;
	--fasta)
		FASTA="$2"
		shift 2
		;;
	--extra-args)
		EXTRA_VEP_ARGS="$2"
		shift 2
		;;
	--gate-f1)
		GATE_F1="$2"
		shift 2
		;;
	--gate-total-delta)
		GATE_TOTAL_DELTA="$2"
		shift 2
		;;
	--gate-coloc-delta)
		GATE_COLOC_DELTA="$2"
		shift 2
		;;
	--gate-freq-delta)
		GATE_FREQ_DELTA="$2"
		shift 2
		;;
	--gate-per-csq-delta)
		GATE_PER_CSQ_DELTA="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERROR: [compare_cache_outputs] Unknown argument: $1" >&2
		exit 1
		;;
	esac
done

# Validate required arguments
if [[ -z "$BASELINE_CACHE" ]]; then
	echo "ERROR: [compare_cache_outputs] --baseline-cache is required" >&2
	exit 1
fi
if [[ -z "$BUILDER_CACHE" ]]; then
	echo "ERROR: [compare_cache_outputs] --builder-cache is required" >&2
	exit 1
fi
if [[ -z "$VCF" ]]; then
	echo "ERROR: [compare_cache_outputs] --vcf is required" >&2
	exit 1
fi
if [[ -z "$VEP_BINARY" ]]; then
	echo "ERROR: [compare_cache_outputs] --vep-binary is required" >&2
	exit 1
fi

# The tuple comparator of Step 4 is this tree's release build of tools/vep-outputs-compare.
# Nothing is built here and no other program runs in its place: a build that is not an
# executable file stops the run here, before vep-rs runs.
OUTPUTS_COMPARE=()
outputs_bin="${VEP_RS_DIR}/tools/vep-outputs-compare/target/release/vep-outputs-compare"
if [[ -f "$outputs_bin" && -x "$outputs_bin" ]]; then
	OUTPUTS_COMPARE=("$outputs_bin")
else
	echo "ERROR: [compare_cache_outputs] the tuple comparator $outputs_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path ${VEP_RS_DIR}/tools/vep-outputs-compare/Cargo.toml'" >&2
	exit 1
fi

# The structural cache comparison of Step 3 is this tree's release build of
# tools/vep-cache-diff, held to the same rule: found here, before vep-rs runs, or refused.
CACHE_DIFF=()
cache_diff_bin="${VEP_RS_DIR}/tools/vep-cache-diff/target/release/vep-cache-diff"
if [[ -f "$cache_diff_bin" && -x "$cache_diff_bin" ]]; then
	CACHE_DIFF=("$cache_diff_bin")
else
	echo "ERROR: [compare_cache_outputs] the cache comparison $cache_diff_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path ${VEP_RS_DIR}/tools/vep-cache-diff/Cargo.toml'" >&2
	exit 1
fi

# The gates of Step 5 are scored by this tree's release build of
# tools/vep-cache-output-gates, held to the same rule: found here, before vep-rs runs, or
# refused.
OUTPUT_GATES=()
output_gates_bin="${VEP_RS_DIR}/tools/vep-cache-output-gates/target/release/vep-cache-output-gates"
if [[ -f "$output_gates_bin" && -x "$output_gates_bin" ]]; then
	OUTPUT_GATES=("$output_gates_bin")
else
	echo "ERROR: [compare_cache_outputs] the gate scorer $output_gates_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path ${VEP_RS_DIR}/tools/vep-cache-output-gates/Cargo.toml'" >&2
	exit 1
fi

if [[ -z "$WORK_DIR" ]]; then
	WORK_DIR="${VEP_RS_DIR}/tmp/cache_comparison"
fi

mkdir -p "$WORK_DIR"

echo "============================================"
echo "  Cache Output Comparison"
echo "============================================"
echo "  Baseline cache: $BASELINE_CACHE"
echo "  Builder cache:  $BUILDER_CACHE"
echo "  VCF input:      $VCF"
echo "  VEP binary:     $VEP_BINARY"
echo "  Work dir:       $WORK_DIR"
echo ""

# Step 1: Run vep-rs with baseline cache
echo "=== Step 1: Annotate with baseline cache ==="

CAPTURE_ARGS=(
	--cache-dir "$BASELINE_CACHE"
	--vcf "$VCF"
	--vep-binary "$VEP_BINARY"
	--output-dir "${WORK_DIR}/baseline"
	--fork "$FORK"
	--buffer-size "$BUFFER_SIZE"
)
if [[ -n "$FASTA" ]]; then
	CAPTURE_ARGS+=(--fasta "$FASTA")
fi
if [[ -n "$EXTRA_VEP_ARGS" ]]; then
	CAPTURE_ARGS+=(--extra-args "$EXTRA_VEP_ARGS")
fi

bash "${SCRIPT_DIR}/capture_baseline.sh" "${CAPTURE_ARGS[@]}"
echo ""

# Step 2: Run vep-rs with builder cache
echo "=== Step 2: Annotate with builder cache ==="

CAPTURE_ARGS=(
	--cache-dir "$BUILDER_CACHE"
	--vcf "$VCF"
	--vep-binary "$VEP_BINARY"
	--output-dir "${WORK_DIR}/builder"
	--fork "$FORK"
	--buffer-size "$BUFFER_SIZE"
)
if [[ -n "$FASTA" ]]; then
	CAPTURE_ARGS+=(--fasta "$FASTA")
fi
if [[ -n "$EXTRA_VEP_ARGS" ]]; then
	CAPTURE_ARGS+=(--extra-args "$EXTRA_VEP_ARGS")
fi

bash "${SCRIPT_DIR}/capture_baseline.sh" "${CAPTURE_ARGS[@]}"
echo ""

# Step 3: Structural cache diff
echo "=== Step 3: Structural cache comparison ==="

"${CACHE_DIFF[@]}" \
	--baseline-dir "$BASELINE_CACHE" \
	--builder-dir "$BUILDER_CACHE" \
	--output-report "${WORK_DIR}/cache_diff.json"
echo ""

# Step 4: Annotation output comparison (the tuple comparator resolved above)
echo "=== Step 4: Annotation output comparison ==="

# The tuple comparator expects directories with paired .txt files, so each single
# output file is symlinked into a directory of its own.
BASELINE_OUT_DIR="${WORK_DIR}/baseline_out"
BUILDER_OUT_DIR="${WORK_DIR}/builder_out"
REPORT_DIR="${WORK_DIR}/concordance_report"
mkdir -p "$BASELINE_OUT_DIR" "$BUILDER_OUT_DIR" "$REPORT_DIR"

ln -sf "${WORK_DIR}/baseline/baseline_output.txt" "${BASELINE_OUT_DIR}/output.txt"
ln -sf "${WORK_DIR}/builder/baseline_output.txt" "${BUILDER_OUT_DIR}/output.txt"

# A comparison that fails writes no summary. A summary an earlier run left in the work
# directory is removed first, so that Step 5 scores a missing summary and never that one.
rm -f "${REPORT_DIR}/concordance.json"

"${OUTPUTS_COMPARE[@]}" \
	--perl-dir "$BASELINE_OUT_DIR" \
	--rust-dir "$BUILDER_OUT_DIR" \
	--report-dir "$REPORT_DIR" \
	--summary-json "concordance.json" || true

CONCORDANCE_REPORT="${REPORT_DIR}/concordance.json"
echo ""

# Step 5: Gate checks and verdict
echo "=== Step 5: Validation gates ==="

# The tool scores the five gates over the two runs' metrics, the cache comparison and the
# concordance summary, prints a line a gate and writes verdict.json. It leaves with the
# verdict: 0 for PASS, 1 for FAIL; any other status is a run that reached no verdict and
# stops here.
GATES_STATUS=0
"${OUTPUT_GATES[@]}" \
	"${WORK_DIR}/baseline/baseline_metrics.json" \
	"${WORK_DIR}/builder/baseline_metrics.json" \
	"${WORK_DIR}/cache_diff.json" \
	"${CONCORDANCE_REPORT:-}" \
	"${WORK_DIR}/verdict.json" \
	"$GATE_F1" \
	"$GATE_TOTAL_DELTA" \
	"$GATE_COLOC_DELTA" \
	"$GATE_FREQ_DELTA" \
	"$GATE_PER_CSQ_DELTA" || GATES_STATUS=$?

case "$GATES_STATUS" in
0) VERDICT="PASS" ;;
1) VERDICT="FAIL" ;;
*) exit 1 ;;
esac

echo ""
echo "============================================"
echo "  VERDICT: $VERDICT"
echo "============================================"

if [[ "$VERDICT" == "FAIL" ]]; then
	exit 1
fi
