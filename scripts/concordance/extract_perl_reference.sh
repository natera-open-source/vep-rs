#!/usr/bin/env bash
# Extract five Perl VEP modules from the Docker image into a local directory, for
# reading and tracing beside the Rust port. Re-run when the image version changes.
#
# Usage: scripts/concordance/extract_perl_reference.sh [--output-dir DIR] [--image REF]
set -euo pipefail

# Defaults
OUTPUT_DIR="${VEP_PERL_REF_DIR:-./.vep/perl_reference}"
DOCKER_IMAGE="ensemblorg/ensembl-vep:release_115.2"
PERL_BASE="/opt/vep/src/ensembl-vep"

usage() {
    cat <<'EOF'
Extract five Perl VEP modules from the Docker image for reference and debugging.

Optional:
  --output-dir <path>   Where to write extracted files (default: ${VEP_PERL_REF_DIR:-./.vep/perl_reference})
  --image <ref>         Docker image (default: ensemblorg/ensembl-vep:release_115.2)
  -h, --help            Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    --output-dir)
        OUTPUT_DIR="$2"
        shift 2
        ;;
    --image)
        DOCKER_IMAGE="$2"
        shift 2
        ;;
    -h | --help) usage ;;
    *)
        echo "ERROR: Unknown argument: $1"
        usage
        ;;
    esac
done

# Files to extract, as two parallel arrays (bash 3.2, which macOS ships, has no
# associative arrays): path inside the image relative to PERL_BASE, and local name.
# The ensembl-variation modules sit directly under Bio/; VEP's own modules sit
# under modules/Bio/.
DOCKER_RELS=(
    "Bio/EnsEMBL/Variation/Utils/VariationEffect.pm"
    "Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm"
    "Bio/EnsEMBL/Variation/StructuralVariationOverlapAllele.pm"
    "Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm"
    "modules/Bio/EnsEMBL/VEP/Parser/VCF.pm"
)
LOCAL_NAMES=(
    "VariationEffect.pm"
    "TranscriptVariationAllele.pm"
    "StructuralVariationOverlapAllele.pm"
    "BaseTranscriptVariation.pm"
    "VCF.pm"
)

echo "=== Extract Perl VEP Reference Modules ==="
echo "Docker image: $DOCKER_IMAGE"
echo "Output dir:   $OUTPUT_DIR"
echo "Started:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ""

# Check Docker is available
command -v docker >/dev/null 2>&1 || {
    echo "ERROR: docker not found in PATH"
    exit 1
}

# Verify image exists (pull if needed)
if ! docker image inspect "$DOCKER_IMAGE" >/dev/null 2>&1; then
    echo "Pulling Docker image (this may take a few minutes)..."
    docker pull "$DOCKER_IMAGE"
fi

mkdir -p "$OUTPUT_DIR"

echo "--- Extracting files ---"
FAILURES=0
for i in "${!DOCKER_RELS[@]}"; do
    docker_rel="${DOCKER_RELS[$i]}"
    local_name="${LOCAL_NAMES[$i]}"
    docker_path="${PERL_BASE}/${docker_rel}"
    output_path="${OUTPUT_DIR}/${local_name}"

    if docker run --rm "$DOCKER_IMAGE" cat "$docker_path" >"$output_path" 2>/dev/null; then
        echo "  OK   ${docker_rel} -> ${local_name}"
    else
        echo "  FAIL ${docker_rel} (not found in image)"
        rm -f "$output_path"
        FAILURES=$((FAILURES + 1))
    fi
done

if [[ $FAILURES -gt 0 ]]; then
    echo ""
    echo "ERROR: $FAILURES file(s) failed to extract"
    exit 1
fi

# Make all extracted files read-only
chmod 444 "$OUTPUT_DIR"/*.pm

# md5sum is GNU coreutils; macOS ships md5 instead.
file_md5() {
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "$1" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then
        md5 -q "$1"
    else
        echo "n/a"
    fi
}

echo ""
echo "--- Verification ---"
printf "%-45s %10s %10s\n" "File" "Lines" "MD5"
for local_name in $(printf '%s\n' "${LOCAL_NAMES[@]}" | sort); do
    fpath="${OUTPUT_DIR}/${local_name}"
    lines=$(wc -l <"$fpath" | tr -d ' ')
    md5=$(file_md5 "$fpath")
    printf "%-45s %10d %10s\n" "$local_name" "$lines" "${md5:0:12}..."
done

echo ""
echo "=== Complete: $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "Files at: $OUTPUT_DIR/"
