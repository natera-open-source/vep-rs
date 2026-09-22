#!/usr/bin/env bash
set -euo pipefail

# validate_builder_cache.sh
#
# Validates the vep-cache-builder output against either captured fixtures or a
# Perl-derived cache.
#
# Usage:
#   scripts/validation/validate_builder_cache.sh --fixture-mode
#   scripts/validation/validate_builder_cache.sh --baseline-cache /path/to/perl_cache

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VEP_RS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Defaults
FIXTURE_MODE=0
BASELINE_CACHE=""
RELEASE=115
ASSEMBLY="GRCh37"
SPECIES="homo_sapiens"
REGION="21"
WORK_DIR=""
FIXTURE_DIR="${VEP_RS_DIR}/tests/cache_builder_validation"
GFF3=""
FASTA=""

# Usage
usage() {
	cat <<'EOF'
Validate vep-cache-builder output against a reference.

Usage:
  scripts/validation/validate_builder_cache.sh [options]

Modes (choose one):
  --fixture-mode             Compare against fixtures in tests/cache_builder_validation/
                             (a build product of this script, not committed; see
                             --fixture-dir. Exits 3 when they are absent.)
  --baseline-cache <path>    Compare against a Perl-derived JSON cache

Options:
  --release <int>            Ensembl release number (default: 115)
  --assembly <name>          Assembly (default: GRCh37)
  --species <name>           Species (default: homo_sapiens)
  --region <chr[,chr...]>    Chromosome(s) to build, passed to the builder as
                             --chromosomes (default: 21)
  --work-dir <path>          Working directory (default: vep-rs/tmp/builder_validation)
  --fixture-dir <path>       Fixture directory (default: tests/cache_builder_validation/,
                             which is not committed and must be created first)
  --gff3 <path>              Local GFF3 file (skip download)
  --fasta <path>             Local protein FASTA file (skip download)
  --help                     Print this help

Exit codes:
  0    Validation passed
  1    Validation failed or error
  3    Could not compare: --fixture-mode with no fixture directory
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
	case "$1" in
	--fixture-mode)
		FIXTURE_MODE=1
		shift
		;;
	--baseline-cache)
		BASELINE_CACHE="$2"
		shift 2
		;;
	--release)
		RELEASE="$2"
		shift 2
		;;
	--assembly)
		ASSEMBLY="$2"
		shift 2
		;;
	--species)
		SPECIES="$2"
		shift 2
		;;
	--region)
		REGION="$2"
		shift 2
		;;
	--work-dir)
		WORK_DIR="$2"
		shift 2
		;;
	--fixture-dir)
		FIXTURE_DIR="$2"
		shift 2
		;;
	--gff3)
		GFF3="$2"
		shift 2
		;;
	--fasta)
		FASTA="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERROR: [validate_builder_cache] Unknown argument: $1" >&2
		exit 1
		;;
	esac
done

# Validate arguments
if [[ "$FIXTURE_MODE" -eq 0 && -z "$BASELINE_CACHE" ]]; then
	echo "ERROR: [validate_builder_cache] Either --fixture-mode or --baseline-cache is required" >&2
	exit 1
fi

if [[ -z "$WORK_DIR" ]]; then
	WORK_DIR="${VEP_RS_DIR}/tmp/builder_validation"
fi

mkdir -p "$WORK_DIR"

BUILDER_OUTPUT="${WORK_DIR}/builder_cache"

echo "============================================"
echo "  Cache Builder Validation"
echo "============================================"
echo "  Mode:     $(if [[ $FIXTURE_MODE -eq 1 ]]; then echo 'fixture'; else echo 'baseline comparison'; fi)"
echo "  Release:  $RELEASE"
echo "  Assembly: $ASSEMBLY"
echo "  Region:   $REGION"
echo "  Work dir: $WORK_DIR"
echo ""

# Step 1: Build the vep-cache-builder binary
echo "=== Step 1: Building vep-cache-builder ==="

cd "$VEP_RS_DIR"

if ! cargo build --release --bin vep-cache-builder 2>"${WORK_DIR}/build.log"; then
	echo "ERROR: [validate_builder_cache] cargo build failed. See ${WORK_DIR}/build.log" >&2
	tail -20 "${WORK_DIR}/build.log" >&2
	exit 1
fi

BUILDER_BINARY="${VEP_RS_DIR}/target/release/vep-cache-builder"

if [[ ! -x "$BUILDER_BINARY" ]]; then
	echo "ERROR: [validate_builder_cache] builder binary not found at $BUILDER_BINARY" >&2
	exit 1
fi
echo "  Binary: $BUILDER_BINARY"
echo ""

# Step 2: Generate cache for test region
echo "=== Step 2: Generating cache for chr${REGION} ==="

rm -rf "$BUILDER_OUTPUT"
mkdir -p "$BUILDER_OUTPUT"

BUILDER_ARGS=(
	"$BUILDER_BINARY"
	--release "$RELEASE"
	--assembly "$ASSEMBLY"
	--species "$SPECIES"
	--chromosomes "$REGION"
	--output-dir "$BUILDER_OUTPUT"
)

if [[ -n "$GFF3" ]]; then
	BUILDER_ARGS+=(--gff3 "$GFF3")
fi
if [[ -n "$FASTA" ]]; then
	BUILDER_ARGS+=(--fasta "$FASTA")
fi

START_TIME=$(date +%s)
if ! "${BUILDER_ARGS[@]}" 2>"${WORK_DIR}/builder.log"; then
	echo "ERROR: [validate_builder_cache] cache builder failed. See ${WORK_DIR}/builder.log" >&2
	tail -30 "${WORK_DIR}/builder.log" >&2
	exit 1
fi
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "  Cache generated in ${ELAPSED}s at $BUILDER_OUTPUT"

# Verify output has content
TRANSCRIPT_COUNT=$(find "$BUILDER_OUTPUT" -path "*/transcripts/*.json" -type f 2>/dev/null | wc -l | tr -d ' ')
echo "  Transcript region files: $TRANSCRIPT_COUNT"

if [[ "$TRANSCRIPT_COUNT" -eq 0 ]]; then
	echo "ERROR: [validate_builder_cache] builder produced no transcript files" >&2
	exit 1
fi
echo ""

# Step 3: Compare against reference
echo "=== Step 3: Comparing against reference ==="

if [[ "$FIXTURE_MODE" -eq 1 ]]; then
	# No fixtures means there is nothing to compare against, which is NOT a pass.
	# Fixtures are not committed: they are a build product of this script. Exit 3
	# (distinct from a real mismatch) so a caller can tell "could not compare" from
	# "compared and differed".
	if [[ ! -d "$FIXTURE_DIR" ]]; then
		echo "ERROR: [validate_builder_cache] no fixtures at $FIXTURE_DIR" >&2
		echo "  Nothing to compare against, so this run proves nothing." >&2
		echo "" >&2
		echo "  Create them once from a build you trust:" >&2
		echo "    cp -r ${BUILDER_OUTPUT} ${FIXTURE_DIR}" >&2
		echo "" >&2
		echo "  Or compare against a Perl-derived cache instead:" >&2
		echo "    --baseline-cache <path>" >&2
		echo "============================================" >&2
		echo "  VERDICT: COULD NOT COMPARE (no fixtures)" >&2
		echo "============================================" >&2
		exit 3
	fi

	REFERENCE_DIR="$FIXTURE_DIR"
else
	# Baseline mode: compare against Perl-derived cache
	if [[ ! -d "$BASELINE_CACHE" ]]; then
		echo "ERROR: [validate_builder_cache] baseline cache not found: $BASELINE_CACHE" >&2
		exit 1
	fi
	REFERENCE_DIR="$BASELINE_CACHE"
fi

echo "  Reference: $REFERENCE_DIR"
echo "  Builder:   $BUILDER_OUTPUT"

python3 "${SCRIPT_DIR}/diff_caches.py" \
	--baseline-dir "$REFERENCE_DIR" \
	--builder-dir "$BUILDER_OUTPUT" \
	--output-report "${WORK_DIR}/validation_diff.json"

# Step 4: Check results
echo ""
echo "=== Step 4: Checking results ==="

RESULT=$(python3 -c "
import json, sys
with open('${WORK_DIR}/validation_diff.json') as fh:
    report = json.load(fh)

tr = report.get('transcripts', {})
var = report.get('variations', {})
info = report.get('info_json', {})

tr_rate = tr.get('match_rate', 0.0)
var_rate = var.get('match_rate', 0.0)
tr_only_b = tr.get('only_baseline', 0)
tr_only_c = tr.get('only_builder', 0)

print(f'Transcript match rate: {tr_rate}')
print(f'Variation match rate:  {var_rate}')
print(f'Transcripts only in baseline: {tr_only_b}')
print(f'Transcripts only in builder:  {tr_only_c}')
print(f'info.json matches: {info.get(\"matches\", \"N/A\")}')

# In fixture mode, require exact match
# In baseline mode, allow some tolerance (builder may add/improve)
if ${FIXTURE_MODE} == 1:
    if tr_rate == 1.0 and tr_only_b == 0 and tr_only_c == 0:
        print('VERDICT:PASS')
    else:
        print('VERDICT:FAIL')
else:
    # Baseline comparison: structural differences are expected, so the
    # transcript match rate has a 0.95 floor instead of an exact-match test.
    if tr_rate >= 0.95:
        print('VERDICT:PASS')
    else:
        print('VERDICT:FAIL')
")

echo "$RESULT" | grep -v '^VERDICT:'

VERDICT=$(echo "$RESULT" | grep '^VERDICT:' | cut -d: -f2)

echo ""
echo "============================================"
echo "  VERDICT: $VERDICT"
echo "============================================"

if [[ "$VERDICT" == "FAIL" ]]; then
	echo ""
	echo "Full diff report: ${WORK_DIR}/validation_diff.json"
	exit 1
fi
