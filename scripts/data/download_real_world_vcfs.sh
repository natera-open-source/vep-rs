#!/usr/bin/env bash
set -euo pipefail

# download_real_world_vcfs.sh
#
# Idempotent seeder for the 14 real-world SOURCE VCFs used in concordance testing.
# Ensures every VCF exists locally (and, if an S3 bucket is configured, in S3).
# For each file:
#   1. Already on disk?  -> skip
#   2. In S3 bucket?     -> aws s3 cp to local        (only if a bucket is set)
#   3. Neither?          -> download from upstream, filter/index, copy to S3
#
# The S3 copy is optional: name the bucket with --s3-bucket, $VEP_BENCHMARK_S3
# (an s3:// URI) or $VEP_BENCHMARK_BUCKET (a bare bucket name, the variable
# setup_plugin_data.sh reads); a bare name is accepted everywhere and prefixed
# with s3://. With no bucket configured the script uses local disk + public
# upstream sources only (all aws calls are skipped).
#
# Covers both structural variant (SV) and all-variant (SNP/indel/multi-allelic)
# test suites for GRCh37 and GRCh38 assemblies, all filtered to chr21 except the
# two genome-wide ClinVar files.
#
# NOT produced here, and derived from these files by the caller:
#   <asm>/all_variants_canonical/<base>.canonical.vcf.gz   Ensembl-style contigs
#       (prepare_benchmark_vcfs.py, or bcftools view -t), the inputs the
#       measurement harness and the reference generator read;
#   the sites-only and genotype-carrying forms of the two 1000 Genomes files
#       (bcftools view -G drops genotypes) under the names the harness's suite
#       table gives them: 1kg_integrated_chr21_sitesonly and 1kg_highcov_chr21
#       (sites-only), 1kg_integrated_chr21 and 1kg_highcov_chr21_multisample
#       (with genotypes).
#
# Target layout (under OUTPUT_DIR):
#   grch37/sv/
#     clinvar_sv_chr21.vcf.gz          ClinVar GRCh37, chr21 SVs
#     gnomad_sv_v2.1_chr21.vcf.gz      gnomAD-SV v2.1 GRCh37, chr21
#     1kg_sv_chr21.vcf.gz              1000G Phase 3 SV GRCh37, chr21
#   grch37/all_variants/
#     clinvar_grch37.vcf.gz               ClinVar GRCh37, all chr, all variants
#     clinvar_chr21_grch37.vcf.gz         ClinVar GRCh37, chr21
#     gnomad_genomes_v2.1.1_chr21.vcf.gz  gnomAD Genomes v2.1.1 GRCh37, chr21
#     1kg_integrated_chr21.vcf.gz          1000G Phase 3 integrated GRCh37, chr21
#   grch38/sv/
#     clinvar_sv_chr21.vcf.gz          ClinVar GRCh38, chr21 SVs
#     gnomad_sv_v4.1_chr21.vcf.gz      gnomAD-SV v4.1 GRCh38, chr21
#     1kg_sv_chr21.vcf.gz              1000G high-cov SV GRCh38, chr21
#   grch38/all_variants/
#     clinvar_grch38.vcf.gz            ClinVar GRCh38, all chr, all variants
#     clinvar_chr21_grch38.vcf.gz      ClinVar GRCh38, chr21
#     gnomad_genomes_v4.1_chr21.vcf.gz gnomAD Genomes v4.1 GRCh38, chr21
#     1kg_highcov_chr21.vcf.gz         1000G high-cov integrated GRCh38, chr21
#
# Usage:
#   scripts/data/download_real_world_vcfs.sh [--output-dir DIR] [--s3-bucket URI|NAME] [--force]
#   plus --profile, which selects the AWS credentials for the S3 copy.

# Defaults
OUTPUT_DIR="${VEP_VCF_DIR:-./.vep/vcf}"
S3_BUCKET="${VEP_BENCHMARK_S3:-${VEP_BENCHMARK_BUCKET:-}}"
AWS_PROFILE="${AWS_PROFILE:-default}"
FORCE=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
    --output-dir)
        OUTPUT_DIR="$2"
        shift 2
        ;;
    --s3-bucket)
        S3_BUCKET="$2"
        shift 2
        ;;
    --profile)
        AWS_PROFILE="$2"
        shift 2
        ;;
    --force)
        FORCE=1
        shift
        ;;
    -h | --help)
        echo "Usage: $0 [--output-dir DIR] [--s3-bucket URI|NAME] [--profile PROFILE] [--force]"
        echo ""
        echo "  --output-dir DIR     Local output directory (default: \${VEP_VCF_DIR:-./.vep/vcf})"
        echo "  --s3-bucket URI|NAME Optional S3 bucket holding a copy of the downloads (default: none;"
        echo "                       set --s3-bucket, \$VEP_BENCHMARK_S3 or \$VEP_BENCHMARK_BUCKET)"
        echo "  --profile PROFILE    AWS CLI profile (default: default)"
        echo "  --force              Re-download even if target files already exist"
        echo ""
        echo "Idempotent: downloads 14 real-world VCFs (SV + all-variants, both assemblies)."
        echo "Checks local disk first, then (if an S3 bucket is set) S3, only hitting upstream when missing everywhere."
        echo "When an S3 bucket is set, backfills it for any file downloaded from upstream."
        exit 0
        ;;
    *)
        echo "ERROR: Unknown argument: $1" >&2
        exit 1
        ;;
    esac
done

# A bare bucket name becomes an s3:// URI, so both spellings work.
if [[ -n "$S3_BUCKET" && "$S3_BUCKET" != s3://* ]]; then
    S3_BUCKET="s3://${S3_BUCKET}"
fi

# Verify dependencies
REQUIRED_CMDS=(curl bcftools bgzip tabix)
# aws is only needed when mirroring to/from an S3 bucket
[[ -n "$S3_BUCKET" ]] && REQUIRED_CMDS+=(aws)
for cmd in "${REQUIRED_CMDS[@]}"; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: Required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# Directory setup
GRCH37_SV_DIR="${OUTPUT_DIR}/grch37/sv"
GRCH37_AV_DIR="${OUTPUT_DIR}/grch37/all_variants"
GRCH38_SV_DIR="${OUTPUT_DIR}/grch38/sv"
GRCH38_AV_DIR="${OUTPUT_DIR}/grch38/all_variants"
TMPDIR_BASE=$(mktemp -d "${TMPDIR:-/tmp}/real_world_vcfs.XXXXXX")

mkdir -p "$GRCH37_SV_DIR" "$GRCH37_AV_DIR" "$GRCH38_SV_DIR" "$GRCH38_AV_DIR"

trap 'rm -rf "$TMPDIR_BASE"' EXIT

echo "=== Real-World VCF Download Script ==="
echo "Output dir : $OUTPUT_DIR"
if [[ -n "$S3_BUCKET" ]]; then
    echo "S3 bucket  : $S3_BUCKET"
    echo "AWS profile: $AWS_PROFILE"
else
    echo "S3 bucket  : (none; local + upstream only)"
fi
echo "Force      : $FORCE"
echo ""

# Counters
TOTAL=0
SKIPPED=0
FROM_S3=0
FROM_UPSTREAM=0

# Helper: count variants in a VCF
count_variants() {
    local vcf="$1"
    bcftools view -H "$vcf" 2>/dev/null | wc -l | tr -d ' '
}

# Helper: detect chromosome naming convention in a VCF
detect_chr21_name() {
    local vcf="$1"
    if bcftools view -h "$vcf" 2>/dev/null | grep -qE '##contig=<ID=chr21'; then
        echo "chr21"
    elif bcftools view -h "$vcf" 2>/dev/null | grep -qE '##contig=<ID=21[,>]'; then
        echo "21"
    else
        local first_chr
        first_chr=$(bcftools view -H "$vcf" 2>/dev/null | head -1 | cut -f1)
        if [[ "$first_chr" == chr* ]]; then
            echo "chr21"
        else
            echo "21"
        fi
    fi
}

# Helper: download a URL to a local path (with progress)
download_url() {
    local url="$1"
    local dest="$2"
    echo "    Downloading: $url"
    curl -L --progress-bar -o "$dest" "$url"
}

# Helper: check if file exists in S3
s3_exists() {
    local s3_key="$1"
    # No bucket configured: behave as "not in S3" so callers fall through to upstream.
    [[ -z "$S3_BUCKET" ]] && return 1
    aws s3 ls "${S3_BUCKET}/${s3_key}" --profile "$AWS_PROFILE" >/dev/null 2>&1
}

# Helper: copy file to S3 (backfill)
s3_upload() {
    local local_path="$1"
    local s3_key="$2"
    # No bucket configured: nothing to mirror.
    [[ -z "$S3_BUCKET" ]] && return 0
    echo "    Uploading to S3: ${S3_BUCKET}/${s3_key}"
    aws s3 cp "$local_path" "${S3_BUCKET}/${s3_key}" --profile "$AWS_PROFILE" --quiet
    # Also upload index if it exists
    if [[ -f "${local_path}.tbi" ]]; then
        aws s3 cp "${local_path}.tbi" "${S3_BUCKET}/${s3_key}.tbi" --profile "$AWS_PROFILE" --quiet
    fi
}

# Helper: copy file from S3 to local
s3_download() {
    local s3_key="$1"
    local local_path="$2"
    [[ -z "$S3_BUCKET" ]] && return 1
    echo "    Fetching from S3: ${S3_BUCKET}/${s3_key}"
    aws s3 cp "${S3_BUCKET}/${s3_key}" "$local_path" --profile "$AWS_PROFILE" --quiet
    # Also fetch index if it exists
    if s3_exists "${s3_key}.tbi"; then
        aws s3 cp "${S3_BUCKET}/${s3_key}.tbi" "${local_path}.tbi" --profile "$AWS_PROFILE" --quiet
    fi
}

# Core: ensure a simple VCF exists (direct download, no filtering needed)
#   ensure_simple_vcf TARGET_PATH S3_KEY UPSTREAM_URL [INDEX_URL]
ensure_simple_vcf() {
    local target="$1"
    local s3_key="$2"
    local upstream_url="$3"
    local index_url="${4:-${upstream_url}.tbi}"

    TOTAL=$((TOTAL + 1))

    # 1. Already on disk?
    if [[ -f "$target" && "$FORCE" -eq 0 ]]; then
        echo "  SKIP (local): $(basename "$target") ($(count_variants "$target") variants)"
        SKIPPED=$((SKIPPED + 1))
        # Backfill S3 if missing
        if ! s3_exists "$s3_key"; then
            s3_upload "$target" "$s3_key"
        fi
        return 0
    fi

    # 2. In S3?
    if s3_exists "$s3_key" && [[ "$FORCE" -eq 0 ]]; then
        s3_download "$s3_key" "$target"
        echo "  FROM S3: $(basename "$target") ($(count_variants "$target") variants)"
        FROM_S3=$((FROM_S3 + 1))
        return 0
    fi

    # 3. Download from upstream
    local tmp_vcf="${TMPDIR_BASE}/$(basename "$target")"
    download_url "$upstream_url" "$tmp_vcf"
    if [[ -n "$index_url" ]]; then
        download_url "$index_url" "${tmp_vcf}.tbi"
    fi
    cp "$tmp_vcf" "$target"
    [[ -f "${tmp_vcf}.tbi" ]] && cp "${tmp_vcf}.tbi" "${target}.tbi"

    # Index if no .tbi yet
    if [[ ! -f "${target}.tbi" ]]; then
        tabix -p vcf "$target"
    fi

    echo "  FROM UPSTREAM: $(basename "$target") ($(count_variants "$target") variants)"
    FROM_UPSTREAM=$((FROM_UPSTREAM + 1))

    # Backfill S3
    s3_upload "$target" "$s3_key"
}

# Core: ensure a filtered VCF exists (download genome-wide, filter to chr21)
#   ensure_filtered_vcf TARGET_PATH S3_KEY UPSTREAM_URL CHR_REGION [FILTER_EXPR]
ensure_filtered_vcf() {
    local target="$1"
    local s3_key="$2"
    local upstream_url="$3"
    local chr_region="$4"
    local filter_expr="${5:-}"

    TOTAL=$((TOTAL + 1))

    # 1. Already on disk?
    if [[ -f "$target" && "$FORCE" -eq 0 ]]; then
        echo "  SKIP (local): $(basename "$target") ($(count_variants "$target") variants)"
        SKIPPED=$((SKIPPED + 1))
        if ! s3_exists "$s3_key"; then
            s3_upload "$target" "$s3_key"
        fi
        return 0
    fi

    # 2. In S3?
    if s3_exists "$s3_key" && [[ "$FORCE" -eq 0 ]]; then
        s3_download "$s3_key" "$target"
        echo "  FROM S3: $(basename "$target") ($(count_variants "$target") variants)"
        FROM_S3=$((FROM_S3 + 1))
        return 0
    fi

    # 3. Download from upstream, filter, index
    local tmp_vcf="${TMPDIR_BASE}/$(basename "$target" .vcf.gz)_full.vcf.gz"
    download_url "$upstream_url" "$tmp_vcf"
    download_url "${upstream_url}.tbi" "${tmp_vcf}.tbi"

    # Auto-detect chr21 naming if chr_region is "auto"
    if [[ "$chr_region" == "auto" ]]; then
        chr_region=$(detect_chr21_name "$tmp_vcf")
        echo "    Detected chr21 name: $chr_region"
    fi

    echo "    Filtering to ${chr_region}..."
    if [[ -n "$filter_expr" ]]; then
        bcftools view -r "$chr_region" "$tmp_vcf" |
            bcftools view -i "$filter_expr" - |
            bgzip -c >"$target"
    else
        bcftools view -r "$chr_region" "$tmp_vcf" |
            bgzip -c >"$target"
    fi
    tabix -p vcf "$target"

    echo "  FROM UPSTREAM: $(basename "$target") ($(count_variants "$target") variants)"
    FROM_UPSTREAM=$((FROM_UPSTREAM + 1))

    # Backfill S3
    s3_upload "$target" "$s3_key"
}

# ClinVar release pins. The rolling vcf_GRCh3x/clinvar.vcf.gz changes weekly, so
# every ClinVar-derived input (genome-wide, chr21 subset, chr21 SV subset) reads
# NCBI's dated archive copy, vcf_GRCh3x/archive_2.0/<year>/clinvar_<release>.vcf.gz,
# of the release the benchmark inputs were built from; <release> is that file's own
# ##fileDate. The GRCh38 structural-variant subset was cut from release 20260302 and
# the GRCh38 genome-wide and chr21 inputs from release 20260321, so the two carry
# separate pins. A release moves here and in the matching source_release row of
# manuscript/data/dataset_inventory.csv together.
CLINVAR_GRCH37_RELEASE="20260302"
CLINVAR_GRCH38_RELEASE="20260321"
CLINVAR_GRCH38_SV_RELEASE="20260302"
CLINVAR_GRCH37_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/${CLINVAR_GRCH37_RELEASE:0:4}/clinvar_${CLINVAR_GRCH37_RELEASE}.vcf.gz"
CLINVAR_GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/${CLINVAR_GRCH38_RELEASE:0:4}/clinvar_${CLINVAR_GRCH38_RELEASE}.vcf.gz"
CLINVAR_GRCH38_SV_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/${CLINVAR_GRCH38_SV_RELEASE:0:4}/clinvar_${CLINVAR_GRCH38_SV_RELEASE}.vcf.gz"

# ClinVar SV filter expression (shared by both assemblies)
CLNVC_SV_FILTER='CLNVC="Deletion" || CLNVC="Duplication" || CLNVC="Insertion" || CLNVC="Inversion" || CLNVC="copy_number_gain" || CLNVC="copy_number_loss" || CLNVC="copy number gain" || CLNVC="copy number loss"'

# GRCh37 SV VCFs (3 files)
echo "--- [1/14] ClinVar GRCh37 SV chr21 ---"
ensure_filtered_vcf \
    "${GRCH37_SV_DIR}/clinvar_sv_chr21.vcf.gz" \
    "vcf/grch37/sv/clinvar_sv_chr21.vcf.gz" \
    "$CLINVAR_GRCH37_URL" \
    "21" \
    "$CLNVC_SV_FILTER"

echo ""
echo "--- [2/14] gnomAD-SV v2.1 GRCh37 chr21 ---"
ensure_filtered_vcf \
    "${GRCH37_SV_DIR}/gnomad_sv_v2.1_chr21.vcf.gz" \
    "vcf/grch37/sv/gnomad_sv_v2.1_chr21.vcf.gz" \
    "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz" \
    "21"

echo ""
echo "--- [3/14] 1000G SV GRCh37 chr21 ---"
# 1000G Phase 3 SV: this is a smaller callset; download and filter
ensure_filtered_vcf \
    "${GRCH37_SV_DIR}/1kg_sv_chr21.vcf.gz" \
    "vcf/grch37/sv/1kg_sv_chr21.vcf.gz" \
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz" \
    "21"

# GRCh37 All-Variants VCFs (4 files)
echo ""
echo "--- [4/14] ClinVar GRCh37 Full (all chr) ---"
# ClinVar GRCh37, full genome, all variant types. This is the primary concordance VCF.
# Downloads genome-wide ClinVar, no chr filtering. Stored alongside the chr21 subset.
ensure_simple_vcf \
    "${GRCH37_AV_DIR}/clinvar_grch37.vcf.gz" \
    "vcf/grch37/all_variants/clinvar_grch37.vcf.gz" \
    "$CLINVAR_GRCH37_URL" \
    "${CLINVAR_GRCH37_URL}.tbi"

echo ""
echo "--- [5/14] ClinVar GRCh37 chr21 ---"
ensure_filtered_vcf \
    "${GRCH37_AV_DIR}/clinvar_chr21_grch37.vcf.gz" \
    "vcf/grch37/all_variants/clinvar_chr21_grch37.vcf.gz" \
    "$CLINVAR_GRCH37_URL" \
    "21"

echo ""
echo "--- [6/14] gnomAD Genomes v2.1.1 GRCh37 chr21 ---"
# Already a chr21-only file: direct download, no filtering
ensure_simple_vcf \
    "${GRCH37_AV_DIR}/gnomad_genomes_v2.1.1_chr21.vcf.gz" \
    "vcf/grch37/all_variants/gnomad_genomes_v2.1.1_chr21.vcf.gz" \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.21.vcf.bgz" \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.21.vcf.bgz.tbi"

echo ""
echo "--- [7/14] 1000G Phase 3 Integrated GRCh37 chr21 ---"
ensure_simple_vcf \
    "${GRCH37_AV_DIR}/1kg_integrated_chr21.vcf.gz" \
    "vcf/grch37/all_variants/1kg_integrated_chr21.vcf.gz" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"

# GRCh38 SV VCFs (3 files)
echo ""
echo "--- [8/14] ClinVar GRCh38 SV chr21 ---"
ensure_filtered_vcf \
    "${GRCH38_SV_DIR}/clinvar_sv_chr21.vcf.gz" \
    "vcf/grch38/sv/clinvar_sv_chr21.vcf.gz" \
    "$CLINVAR_GRCH38_SV_URL" \
    "auto" \
    "$CLNVC_SV_FILTER"

echo ""
echo "--- [9/14] gnomAD-SV v4.1 GRCh38 chr21 ---"
# gnomAD-SV v4.1 genome-wide SV callset: filter to chr21
ensure_filtered_vcf \
    "${GRCH38_SV_DIR}/gnomad_sv_v4.1_chr21.vcf.gz" \
    "vcf/grch38/sv/gnomad_sv_v4.1_chr21.vcf.gz" \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz" \
    "auto"

echo ""
echo "--- [10/14] 1000G SV GRCh38 chr21 ---"
ensure_filtered_vcf \
    "${GRCH38_SV_DIR}/1kg_sv_chr21.vcf.gz" \
    "vcf/grch38/sv/1kg_sv_chr21.vcf.gz" \
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz" \
    "auto"

# GRCh38 All-Variants VCFs (4 files)
echo ""
echo "--- [11/14] ClinVar GRCh38 Full (all chr) ---"
ensure_simple_vcf \
    "${GRCH38_AV_DIR}/clinvar_grch38.vcf.gz" \
    "vcf/grch38/all_variants/clinvar_grch38.vcf.gz" \
    "$CLINVAR_GRCH38_URL" \
    "${CLINVAR_GRCH38_URL}.tbi"

echo ""
echo "--- [12/14] ClinVar GRCh38 chr21 ---"
ensure_filtered_vcf \
    "${GRCH38_AV_DIR}/clinvar_chr21_grch38.vcf.gz" \
    "vcf/grch38/all_variants/clinvar_chr21_grch38.vcf.gz" \
    "$CLINVAR_GRCH38_URL" \
    "auto"

echo ""
echo "--- [13/14] gnomAD Genomes v4.1 GRCh38 chr21 ---"
ensure_simple_vcf \
    "${GRCH38_AV_DIR}/gnomad_genomes_v4.1_chr21.vcf.gz" \
    "vcf/grch38/all_variants/gnomad_genomes_v4.1_chr21.vcf.gz" \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz" \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz.tbi"

echo ""
echo "--- [14/14] 1000G High-Cov GRCh38 chr21 ---"
ensure_simple_vcf \
    "${GRCH38_AV_DIR}/1kg_highcov_chr21.vcf.gz" \
    "vcf/grch38/all_variants/1kg_highcov_chr21.vcf.gz" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"

# Summary
echo ""
echo "=== Summary ==="
echo ""
echo "Sources: $TOTAL total; $SKIPPED skipped (local), $FROM_S3 from S3, $FROM_UPSTREAM from upstream"
echo ""

for dir_label in "GRCh37 SV:${GRCH37_SV_DIR}" "GRCh37 All-Variants:${GRCH37_AV_DIR}" \
    "GRCh38 SV:${GRCH38_SV_DIR}" "GRCh38 All-Variants:${GRCH38_AV_DIR}"; do
    label="${dir_label%%:*}"
    dir="${dir_label##*:}"
    echo "$label ($dir):"
    for f in "$dir"/*.vcf.gz; do
        if [[ -f "$f" ]]; then
            n=$(count_variants "$f")
            size=$(du -h "$f" | cut -f1)
            tbi="${f}.tbi"
            idx_status="NO INDEX"
            [[ -f "$tbi" ]] && idx_status="indexed"
            printf "  %-45s %8s variants  %6s  [%s]\n" "$(basename "$f")" "$n" "$size" "$idx_status"
        fi
    done
    echo ""
done

echo "Done."
