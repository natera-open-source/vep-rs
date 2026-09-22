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

python3 "${SCRIPT_DIR}/diff_caches.py" \
	--baseline-dir "$BASELINE_CACHE" \
	--builder-dir "$BUILDER_CACHE" \
	--output-report "${WORK_DIR}/cache_diff.json"
echo ""

# Step 4: Annotation output comparison (via compare_vep_outputs.py if available)
echo "=== Step 4: Annotation output comparison ==="

CONCORDANCE_SCRIPT="${VEP_RS_DIR}/scripts/concordance/compare_vep_outputs.py"
CONCORDANCE_REPORT=""

if [[ -f "$CONCORDANCE_SCRIPT" ]]; then
	# compare_vep_outputs.py expects directories with paired .txt files, so each
	# single output file is symlinked into a directory of its own.
	BASELINE_OUT_DIR="${WORK_DIR}/baseline_out"
	BUILDER_OUT_DIR="${WORK_DIR}/builder_out"
	REPORT_DIR="${WORK_DIR}/concordance_report"
	mkdir -p "$BASELINE_OUT_DIR" "$BUILDER_OUT_DIR" "$REPORT_DIR"

	ln -sf "${WORK_DIR}/baseline/baseline_output.txt" "${BASELINE_OUT_DIR}/output.txt"
	ln -sf "${WORK_DIR}/builder/baseline_output.txt" "${BUILDER_OUT_DIR}/output.txt"

	python3 "$CONCORDANCE_SCRIPT" \
		--perl-dir "$BASELINE_OUT_DIR" \
		--rust-dir "$BUILDER_OUT_DIR" \
		--report-dir "$REPORT_DIR" \
		--summary-json "concordance.json" || true

	CONCORDANCE_REPORT="${REPORT_DIR}/concordance.json"
	echo ""
else
	echo "  WARNING: compare_vep_outputs.py not found, skipping semantic comparison"
	echo ""
fi

# Step 5: Gate checks and verdict
echo "=== Step 5: Validation gates ==="

python3 - \
	"${WORK_DIR}/baseline/baseline_metrics.json" \
	"${WORK_DIR}/builder/baseline_metrics.json" \
	"${WORK_DIR}/cache_diff.json" \
	"${CONCORDANCE_REPORT:-}" \
	"${WORK_DIR}/verdict.json" \
	"$GATE_F1" \
	"$GATE_TOTAL_DELTA" \
	"$GATE_COLOC_DELTA" \
	"$GATE_FREQ_DELTA" \
	"$GATE_PER_CSQ_DELTA" <<'PYEOF'
import json
import sys
from pathlib import Path

baseline_metrics_path = Path(sys.argv[1])
builder_metrics_path = Path(sys.argv[2])
cache_diff_path = Path(sys.argv[3])
concordance_path = sys.argv[4]  # may be empty string
verdict_path = Path(sys.argv[5])
gate_f1 = float(sys.argv[6])
gate_total_delta = float(sys.argv[7])
gate_coloc_delta = float(sys.argv[8])
gate_freq_delta = float(sys.argv[9])
gate_per_csq_delta = float(sys.argv[10])

with open(baseline_metrics_path) as fh:
    baseline = json.load(fh)
with open(builder_metrics_path) as fh:
    builder = json.load(fh)
with open(cache_diff_path) as fh:
    cache_diff = json.load(fh)

concordance = None
if concordance_path and Path(concordance_path).exists():
    with open(concordance_path) as fh:
        concordance = json.load(fh)

gates = []
all_pass = True

# Gate 1: Consequence F1 (from concordance report)
if concordance and "aggregate" in concordance:
    f1 = concordance["aggregate"].get("f1", 0.0)
    passed = f1 >= gate_f1
    gates.append({
        "gate": "consequence_f1",
        "threshold": gate_f1,
        "actual": round(f1, 6),
        "passed": passed,
    })
    if not passed:
        all_pass = False
    print(f"  Consequence F1: {f1:.6f} (gate >= {gate_f1}) {'PASS' if passed else 'FAIL'}")
else:
    gates.append({
        "gate": "consequence_f1",
        "threshold": gate_f1,
        "actual": None,
        "passed": False,
        "note": "concordance report not available",
    })
    all_pass = False
    print(f"  Consequence F1: N/A (concordance report missing) FAIL")

# Gate 2: Total annotation line delta
b_total = baseline.get("total_annotation_lines", 0)
c_total = builder.get("total_annotation_lines", 0)
if b_total > 0:
    delta = abs(c_total - b_total) / b_total
else:
    delta = 0.0 if c_total == 0 else 1.0
passed = delta < gate_total_delta
gates.append({
    "gate": "total_annotations_delta",
    "threshold": gate_total_delta,
    "baseline": b_total,
    "builder": c_total,
    "actual": round(delta, 6),
    "passed": passed,
})
if not passed:
    all_pass = False
print(f"  Total annotations delta: {delta:.6f} ({b_total} vs {c_total}) (gate < {gate_total_delta}) {'PASS' if passed else 'FAIL'}")

# Gate 3: Co-located variant match delta
b_coloc = baseline.get("existing_variation_count", 0)
c_coloc = builder.get("existing_variation_count", 0)
if b_coloc > 0:
    delta = abs(c_coloc - b_coloc) / b_coloc
else:
    delta = 0.0 if c_coloc == 0 else 1.0
passed = delta < gate_coloc_delta
gates.append({
    "gate": "colocated_variant_delta",
    "threshold": gate_coloc_delta,
    "baseline": b_coloc,
    "builder": c_coloc,
    "actual": round(delta, 6),
    "passed": passed,
})
if not passed:
    all_pass = False
print(f"  Co-located variant delta: {delta:.6f} ({b_coloc} vs {c_coloc}) (gate < {gate_coloc_delta}) {'PASS' if passed else 'FAIL'}")

# Gate 4: Frequency field population delta
b_freqs = baseline.get("frequency_field_counts", {})
c_freqs = builder.get("frequency_field_counts", {})
max_freq_delta = 0.0
freq_details = {}
for field in set(b_freqs) | set(c_freqs):
    bv = b_freqs.get(field, 0)
    cv = c_freqs.get(field, 0)
    if bv > 0:
        d = abs(cv - bv) / bv
    else:
        d = 0.0 if cv == 0 else 1.0
    freq_details[field] = {"baseline": bv, "builder": cv, "delta": round(d, 6)}
    max_freq_delta = max(max_freq_delta, d)

passed = max_freq_delta < gate_freq_delta
gates.append({
    "gate": "frequency_field_delta",
    "threshold": gate_freq_delta,
    "actual": round(max_freq_delta, 6),
    "details": freq_details,
    "passed": passed,
})
if not passed:
    all_pass = False
print(f"  Frequency field max delta: {max_freq_delta:.6f} (gate < {gate_freq_delta}) {'PASS' if passed else 'FAIL'}")

# Gate 5: Per-consequence-type max delta
b_csq = baseline.get("consequence_counts", {})
c_csq = builder.get("consequence_counts", {})
max_csq_delta = 0.0
csq_details = {}
for csq_type in set(b_csq) | set(c_csq):
    bv = b_csq.get(csq_type, 0)
    cv = c_csq.get(csq_type, 0)
    if bv > 0:
        d = abs(cv - bv) / bv
    else:
        d = 0.0 if cv == 0 else 1.0
    csq_details[csq_type] = {"baseline": bv, "builder": cv, "delta": round(d, 6)}
    max_csq_delta = max(max_csq_delta, d)

passed = max_csq_delta < gate_per_csq_delta
gates.append({
    "gate": "per_consequence_type_delta",
    "threshold": gate_per_csq_delta,
    "actual": round(max_csq_delta, 6),
    "details": csq_details,
    "passed": passed,
})
if not passed:
    all_pass = False
print(f"  Per-consequence max delta: {max_csq_delta:.6f} (gate < {gate_per_csq_delta}) {'PASS' if passed else 'FAIL'}")

# Verdict
verdict = {
    "verdict": "PASS" if all_pass else "FAIL",
    "gates": gates,
    "cache_diff_summary": {
        "transcript_match_rate": cache_diff.get("transcripts", {}).get("match_rate"),
        "variation_match_rate": cache_diff.get("variations", {}).get("match_rate"),
        "info_json_matches": cache_diff.get("info_json", {}).get("matches"),
    },
}

with open(verdict_path, "w") as fh:
    json.dump(verdict, fh, indent=2)
    fh.write("\n")

print()
print(f"  Verdict: {'PASS' if all_pass else 'FAIL'}")
print(f"  Report:  {verdict_path}")
PYEOF

VERDICT=$(python3 -c "import json; print(json.load(open('${WORK_DIR}/verdict.json'))['verdict'])")

echo ""
echo "============================================"
echo "  VERDICT: $VERDICT"
echo "============================================"

if [[ "$VERDICT" == "FAIL" ]]; then
	exit 1
fi
