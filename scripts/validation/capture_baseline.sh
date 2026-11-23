#!/usr/bin/env bash
set -euo pipefail

# capture_baseline.sh
#
# Captures baseline metrics from an Ensembl VEP-derived JSON cache by running vep-rs
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

# The metrics of Step 2 and the checksums of Step 3 read the whole output and every cache
# file: they are made by this tree's release build of tools/vep-file-census, across the
# host's cores. Nothing is built here and no other program runs in its place: a tree
# without the build stops here, before vep-rs runs.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
census_bin="$REPO_DIR/tools/vep-file-census/target/release/vep-file-census"
if [[ -f "$census_bin" && -x "$census_bin" ]]; then
	FILE_CENSUS=("$census_bin")
else
	echo "ERROR: [capture_baseline] the file census $census_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-file-census/Cargo.toml'" >&2
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

"${FILE_CENSUS[@]}" metrics "$OUTPUT_FILE" "$METRICS_FILE" "$ELAPSED"

# Step 3: Compute cache file checksums
echo "Step 3: Computing cache checksums..."

find "$CACHE_DIR" -type f -name "*.json" | sort |
	"${FILE_CENSUS[@]}" digest --algorithm md5 --names read >"$CHECKSUMS_FILE"
CHECKSUM_COUNT=$(wc -l <"$CHECKSUMS_FILE" | tr -d ' ')
echo "  Wrote $CHECKSUM_COUNT checksums to $CHECKSUMS_FILE"

echo ""
echo "=== Baseline capture complete ==="
echo "  Metrics: $METRICS_FILE"
echo "  Output:  $OUTPUT_FILE"
echo "  Sums:    $CHECKSUMS_FILE"
