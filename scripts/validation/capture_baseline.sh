#!/usr/bin/env bash
set -euo pipefail

# capture_baseline.sh
#
# Captures baseline metrics from a Perl-derived JSON cache by running vep-rs
# and computing annotation statistics from the output.
#
# Usage:
#   scripts/validation/capture_baseline.sh \
#     --cache-dir /path/to/perl_cache \
#     --vcf /path/to/clinvar.vcf.gz \
#     --vep-binary target/release/vep \
#     --output-dir /tmp/baseline

# Defaults
CACHE_DIR=""
VCF=""
VEP_BINARY=""
OUTPUT_DIR=""
FORK=1
BUFFER_SIZE=5000
FASTA=""
EXTRA_VEP_ARGS=""

# Usage
usage() {
	cat <<'EOF'
Capture baseline annotation metrics from a JSON cache.

Usage:
  scripts/validation/capture_baseline.sh --cache-dir <path> --vcf <path> --vep-binary <path> --output-dir <path> [options]

Required:
  --cache-dir <path>     Path to the JSON cache directory
  --vcf <path>           Input VCF file (plain or bgzipped)
  --vep-binary <path>    Path to the vep-rs binary
  --output-dir <path>    Output directory for metrics and VEP output

Optional:
  --fork <int>           Number of threads for vep-rs (default: 1)
  --buffer-size <int>    Buffer size for vep-rs (default: 5000)
  --fasta <path>         Reference FASTA for HGVS (enables --hgvs)
  --extra-args <args>    Additional vep-rs arguments (quoted string)
  --help                 Print this help

Outputs:
  <output-dir>/baseline_metrics.json   Annotation statistics
  <output-dir>/baseline_output.txt     Full VEP output
  <output-dir>/cache_checksums.txt     MD5 checksums of all cache files
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
	case "$1" in
	--cache-dir)
		CACHE_DIR="$2"
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
	--output-dir)
		OUTPUT_DIR="$2"
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
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERROR: [capture_baseline] Unknown argument: $1" >&2
		exit 1
		;;
	esac
done

# Validate required arguments
if [[ -z "$CACHE_DIR" ]]; then
	echo "ERROR: [capture_baseline] --cache-dir is required" >&2
	exit 1
fi
if [[ -z "$VCF" ]]; then
	echo "ERROR: [capture_baseline] --vcf is required" >&2
	exit 1
fi
if [[ -z "$VEP_BINARY" ]]; then
	echo "ERROR: [capture_baseline] --vep-binary is required" >&2
	exit 1
fi
if [[ -z "$OUTPUT_DIR" ]]; then
	echo "ERROR: [capture_baseline] --output-dir is required" >&2
	exit 1
fi

if [[ ! -d "$CACHE_DIR" ]]; then
	echo "ERROR: [capture_baseline] cache directory not found: $CACHE_DIR" >&2
	exit 1
fi
if [[ ! -f "$VCF" ]]; then
	echo "ERROR: [capture_baseline] VCF file not found: $VCF" >&2
	exit 1
fi
if [[ ! -x "$VEP_BINARY" ]]; then
	echo "ERROR: [capture_baseline] vep binary not found or not executable: $VEP_BINARY" >&2
	exit 1
fi

# Setup
mkdir -p "$OUTPUT_DIR"
OUTPUT_FILE="${OUTPUT_DIR}/baseline_output.txt"
METRICS_FILE="${OUTPUT_DIR}/baseline_metrics.json"
CHECKSUMS_FILE="${OUTPUT_DIR}/cache_checksums.txt"

echo "=== Capture Baseline ==="
echo "  Cache:  $CACHE_DIR"
echo "  VCF:    $VCF"
echo "  Binary: $VEP_BINARY"
echo "  Output: $OUTPUT_DIR"
echo ""

# Step 1: Run vep-rs
echo "Step 1: Running vep-rs..."

VEP_CMD=(
	"$VEP_BINARY"
	-i "$VCF"
	-o "$OUTPUT_FILE"
	--offline
	--json_cache "$CACHE_DIR"
	--fork "$FORK"
	--buffer_size "$BUFFER_SIZE"
	--force
	--no_stats
	--quiet
)

if [[ -n "$FASTA" ]]; then
	VEP_CMD+=(--fasta "$FASTA" --hgvs)
fi

if [[ -n "$EXTRA_VEP_ARGS" ]]; then
	# shellcheck disable=SC2206
	VEP_CMD+=($EXTRA_VEP_ARGS)
fi

START_TIME=$(date +%s)
"${VEP_CMD[@]}"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "  vep-rs completed in ${ELAPSED}s"

# Step 2: Compute metrics from VEP output
echo "Step 2: Computing metrics..."

python3 - "$OUTPUT_FILE" "$METRICS_FILE" "$ELAPSED" <<'PYEOF'
import json
import sys
from collections import Counter
from pathlib import Path

output_file = Path(sys.argv[1])
metrics_file = Path(sys.argv[2])
elapsed_seconds = int(sys.argv[3])

total_lines = 0
consequence_counts = Counter()
existing_variation_count = 0
clin_sig_count = 0
hgvsc_count = 0
hgvsp_count = 0
frequency_fields = {
    "AF": 0,
    "gnomADe_AF": 0,
    "gnomADg_AF": 0,
    "MAX_AF": 0,
    "EUR_AF": 0,
    "AFR_AF": 0,
    "AMR_AF": 0,
    "EAS_AF": 0,
    "SAS_AF": 0,
}

# Parse header to find Extra field index
header_cols = []

with open(output_file) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header_cols = line.lstrip("#").split("\t")
            continue
        if not line:
            continue

        total_lines += 1
        cols = line.split("\t")
        if len(cols) < 7:
            continue

        # Column 6 is Consequence
        consequences = cols[6].strip()
        for csq in consequences.split(","):
            csq = csq.strip()
            if csq:
                consequence_counts[csq] += 1

        # Extra field is the last column (col 13 in standard VEP output)
        extra_idx = -1
        if header_cols:
            for i, h in enumerate(header_cols):
                if h.strip() == "Extra":
                    extra_idx = i
                    break

        extra = ""
        if extra_idx >= 0 and extra_idx < len(cols):
            extra = cols[extra_idx]
        elif len(cols) > 13:
            extra = cols[13]

        extra_fields = {}
        if extra and extra != "-":
            for pair in extra.split(";"):
                if "=" in pair:
                    k, v = pair.split("=", 1)
                    extra_fields[k] = v

        # Existing_variation is col 12 in standard output
        existing_var = ""
        if len(cols) > 12:
            existing_var = cols[12].strip()
        if existing_var and existing_var != "-":
            existing_variation_count += 1

        if extra_fields.get("CLIN_SIG"):
            clin_sig_count += 1
        if extra_fields.get("HGVSc"):
            hgvsc_count += 1
        if extra_fields.get("HGVSp"):
            hgvsp_count += 1

        for freq_field in frequency_fields:
            if extra_fields.get(freq_field):
                frequency_fields[freq_field] += 1

metrics = {
    "total_annotation_lines": total_lines,
    "consequence_counts": dict(consequence_counts.most_common()),
    "existing_variation_count": existing_variation_count,
    "clin_sig_count": clin_sig_count,
    "hgvsc_count": hgvsc_count,
    "hgvsp_count": hgvsp_count,
    "frequency_field_counts": frequency_fields,
    "elapsed_seconds": elapsed_seconds,
}

with open(metrics_file, "w") as fh:
    json.dump(metrics, fh, indent=2)
    fh.write("\n")

print(f"  Total annotation lines: {total_lines}")
print(f"  Unique consequence types: {len(consequence_counts)}")
print(f"  Existing_variation non-empty: {existing_variation_count}")
print(f"  CLIN_SIG non-empty: {clin_sig_count}")
print(f"  HGVSc non-empty: {hgvsc_count}")
print(f"  HGVSp non-empty: {hgvsp_count}")
top_freq = max(frequency_fields.values()) if frequency_fields else 0
print(f"  Max frequency field count: {top_freq}")
PYEOF

# Step 3: Compute cache file checksums
echo "Step 3: Computing cache checksums..."

if command -v md5sum &>/dev/null; then
	MD5CMD="md5sum"
elif command -v md5 &>/dev/null; then
	MD5CMD="md5 -r"
else
	echo "  WARNING: No md5sum/md5 found, skipping checksums"
	MD5CMD=""
fi

if [[ -n "$MD5CMD" ]]; then
	find "$CACHE_DIR" -type f -name "*.json" | sort | while read -r f; do
		$MD5CMD "$f"
	done >"$CHECKSUMS_FILE"
	CHECKSUM_COUNT=$(wc -l <"$CHECKSUMS_FILE" | tr -d ' ')
	echo "  Wrote $CHECKSUM_COUNT checksums to $CHECKSUMS_FILE"
fi

echo ""
echo "=== Baseline capture complete ==="
echo "  Metrics: $METRICS_FILE"
echo "  Output:  $OUTPUT_FILE"
echo "  Sums:    $CHECKSUMS_FILE"
