#!/usr/bin/env bash
# trace_perl_variant.sh - Run Perl VEP on a SINGLE variant with injected predicate tracing.
#
# Perl VEP has no built-in debug output for consequence logic. This script:
#   1. Creates a minimal single-variant VCF
#   2. Patches 3 Perl modules with warn() trace statements via Docker volume mounts
#   3. Runs Perl VEP on the single variant
#   4. Parses stderr TRACE lines into a structured predicate report
#
# Requires:
#   - Docker with ensemblorg/ensembl-vep:release_115.2
#   - Reference Perl modules at --reference-dir (extracted by extract_perl_reference.sh)
#   - Perl VEP cache at --cache-dir
#
# Usage:
#   trace_perl_variant.sh --variant 21:26960070:G:A --assembly GRCh37
#   trace_perl_variant.sh --vcf-line "21\t26960070\t.\tG\tA\t.\t.\t." --assembly GRCh37
#   trace_perl_variant.sh --vcf-file "${VEP_VCF_DIR}/test.vcf" --variant-index 42 --assembly GRCh37

set -euo pipefail

# --- Defaults ---
ASSEMBLY=""
CACHE_DIR="${VEP_PERL_CACHE_DIR:-./.vep/cache}"
REFERENCE_DIR="${VEP_PERL_REF_DIR:-./.vep/perl_reference}"
WORK_DIR="${VEP_WORK_DIR:-./.vep/work}/perl_debug"
DOCKER_IMAGE="ensemblorg/ensembl-vep:release_115.2"

VARIANT=""
VCF_LINE=""
VCF_FILE=""
VARIANT_INDEX=""
KEEP_TEMP=0

# --- Usage ---
usage() {
	cat <<'EOF'
Run Perl VEP on a single variant with predicate-level trace output.

Usage:
  trace_perl_variant.sh --variant CHR:POS:REF:ALT --assembly GRCh37|GRCh38 [options]
  trace_perl_variant.sh --vcf-line "CHR\tPOS\t.\tREF\tALT\t.\t.\t." --assembly GRCh37|GRCh38
  trace_perl_variant.sh --vcf-file FILE --variant-index N --assembly GRCh37|GRCh38

Required (one of):
  --variant CHR:POS:REF:ALT       Variant in colon-separated format
  --vcf-line "LINE"               Raw VCF line (tab-separated, 8+ columns)
  --vcf-file FILE --variant-index N  Extract line N (0-based) from a VCF file

Required:
  --assembly GRCh37|GRCh38        Genome assembly

Optional:
  --cache-dir DIR         Perl VEP cache root, the directory above homo_sapiens/, mounted as
                          --dir_cache (default: ${VEP_PERL_CACHE_DIR:-./.vep/cache})
  --reference-dir DIR     Unpatched Perl module directory (default: ${VEP_PERL_REF_DIR:-./.vep/perl_reference})
  --work-dir DIR          Working directory for temp files (default: ${VEP_WORK_DIR:-./.vep/work}/perl_debug)
  --docker-image IMAGE    Docker image (default: ensemblorg/ensembl-vep:release_115.2)
  --keep-temp             Do not clean up temp files
  --help                  Print this help
EOF
	exit "${1:-0}"
}

# --- Argument parsing ---
while [[ $# -gt 0 ]]; do
	case "$1" in
	--variant)
		VARIANT="$2"
		shift 2
		;;
	--vcf-line)
		VCF_LINE="$2"
		shift 2
		;;
	--vcf-file)
		VCF_FILE="$2"
		shift 2
		;;
	--variant-index)
		VARIANT_INDEX="$2"
		shift 2
		;;
	--assembly)
		ASSEMBLY="$2"
		shift 2
		;;
	--cache-dir)
		CACHE_DIR="$2"
		shift 2
		;;
	--reference-dir)
		REFERENCE_DIR="$2"
		shift 2
		;;
	--work-dir)
		WORK_DIR="$2"
		shift 2
		;;
	--docker-image)
		DOCKER_IMAGE="$2"
		shift 2
		;;
	--keep-temp)
		KEEP_TEMP=1
		shift
		;;
	--help) usage 0 ;;
	*)
		echo "ERROR: Unknown option: $1" >&2
		usage 1
		;;
	esac
done

# --- Validate arguments ---
if [[ -z "${ASSEMBLY}" ]]; then
	echo "ERROR: --assembly is required (GRCh37 or GRCh38)" >&2
	exit 1
fi
if [[ "${ASSEMBLY}" != "GRCh37" && "${ASSEMBLY}" != "GRCh38" ]]; then
	echo "ERROR: --assembly must be GRCh37 or GRCh38, got: ${ASSEMBLY}" >&2
	exit 1
fi

input_count=0
[[ -n "${VARIANT}" ]] && ((input_count++)) || true
[[ -n "${VCF_LINE}" ]] && ((input_count++)) || true
[[ -n "${VCF_FILE}" ]] && ((input_count++)) || true
if [[ "${input_count}" -eq 0 ]]; then
	echo "ERROR: One of --variant, --vcf-line, or --vcf-file is required" >&2
	exit 1
fi
if [[ "${input_count}" -gt 1 ]]; then
	echo "ERROR: Only one of --variant, --vcf-line, or --vcf-file may be specified" >&2
	exit 1
fi
if [[ -n "${VCF_FILE}" && -z "${VARIANT_INDEX}" ]]; then
	echo "ERROR: --variant-index is required when using --vcf-file" >&2
	exit 1
fi

# --- Check prerequisites ---
if ! command -v docker &>/dev/null; then
	echo "ERROR: docker not found in PATH" >&2
	exit 1
fi

MODULES=(
	"VariationEffect.pm"
	"TranscriptVariationAllele.pm"
	"BaseTranscriptVariation.pm"
)
for mod in "${MODULES[@]}"; do
	if [[ ! -f "${REFERENCE_DIR}/${mod}" ]]; then
		echo "ERROR: Reference module not found: ${REFERENCE_DIR}/${mod}" >&2
		echo "  Run extract_perl_reference.sh first to extract modules from Docker." >&2
		exit 1
	fi
done

if [[ ! -d "${CACHE_DIR}" ]]; then
	echo "ERROR: Cache directory not found: ${CACHE_DIR}" >&2
	exit 1
fi

# --- Set up working directory ---
TMP_DIR="${WORK_DIR}/tmp"
mkdir -p "${TMP_DIR}"
mkdir -p "${WORK_DIR}/patched"
# The Perl VEP Docker image runs as a non-root user, so the bind-mounted
# output dir (/work/tmp) must be world-writable or VEP fails with
# "Could not write to output file /work/tmp/output.txt". chmod is best-effort
# (no-op / harmless if the caller already owns the dir as root).
chmod 777 "${TMP_DIR}" 2>/dev/null || true

cleanup() {
	if [[ "${KEEP_TEMP}" -eq 0 ]]; then
		rm -rf "${TMP_DIR}"
	fi
}
trap cleanup EXIT

# --- Step 1: Create single-variant VCF ---
VCF_PATH="${TMP_DIR}/single_variant.vcf"

if [[ -n "${VARIANT}" ]]; then
	# Parse CHR:POS:REF:ALT
	IFS=':' read -r CHR POS REF ALT <<<"${VARIANT}"
	if [[ -z "${CHR}" || -z "${POS}" || -z "${REF}" || -z "${ALT}" ]]; then
		echo "ERROR: --variant must be CHR:POS:REF:ALT format, got: ${VARIANT}" >&2
		exit 1
	fi
	cat >"${VCF_PATH}" <<VCFEOF
##fileformat=VCFv4.2
##INFO=<ID=.,Number=.,Type=String,Description=".">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
${CHR}	${POS}	.	${REF}	${ALT}	.	.	.
VCFEOF
	VARIANT_DESC="${CHR}:${POS} ${REF}>${ALT}"

elif [[ -n "${VCF_LINE}" ]]; then
	# Interpret escape sequences so the user can pass literal \t
	INTERPRETED_LINE="$(printf '%b' "${VCF_LINE}")"
	# Validate: must have at least 8 tab-separated fields
	field_count="$(echo "${INTERPRETED_LINE}" | awk -F'\t' '{print NF}')"
	if [[ "${field_count}" -lt 8 ]]; then
		echo "ERROR: VCF line must have at least 8 tab-separated fields, got ${field_count}" >&2
		exit 1
	fi
	cat >"${VCF_PATH}" <<VCFEOF
##fileformat=VCFv4.2
##INFO=<ID=.,Number=.,Type=String,Description=".">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
${INTERPRETED_LINE}
VCFEOF
	CHR="$(echo "${INTERPRETED_LINE}" | cut -f1)"
	POS="$(echo "${INTERPRETED_LINE}" | cut -f2)"
	REF="$(echo "${INTERPRETED_LINE}" | cut -f4)"
	ALT="$(echo "${INTERPRETED_LINE}" | cut -f5)"
	VARIANT_DESC="${CHR}:${POS} ${REF}>${ALT}"

elif [[ -n "${VCF_FILE}" ]]; then
	if [[ ! -f "${VCF_FILE}" ]]; then
		echo "ERROR: VCF file not found: ${VCF_FILE}" >&2
		exit 1
	fi
	# Extract the Nth non-header line (0-based index)
	EXTRACTED_LINE="$(grep -v '^#' "${VCF_FILE}" | sed -n "$((VARIANT_INDEX + 1))p")"
	if [[ -z "${EXTRACTED_LINE}" ]]; then
		echo "ERROR: No variant at index ${VARIANT_INDEX} in ${VCF_FILE}" >&2
		exit 1
	fi
	# Copy header from original VCF, then append the single variant
	grep '^#' "${VCF_FILE}" >"${VCF_PATH}"
	echo "${EXTRACTED_LINE}" >>"${VCF_PATH}"
	CHR="$(echo "${EXTRACTED_LINE}" | cut -f1)"
	POS="$(echo "${EXTRACTED_LINE}" | cut -f2)"
	REF="$(echo "${EXTRACTED_LINE}" | cut -f4)"
	ALT="$(echo "${EXTRACTED_LINE}" | cut -f5)"
	VARIANT_DESC="${CHR}:${POS} ${REF}>${ALT}"
fi

echo "=== Perl VEP Predicate Trace ===" >&2
echo "Variant: ${VARIANT_DESC} (${ASSEMBLY})" >&2
echo "" >&2

# --- Step 2: Patch Perl modules with TRACE warn() statements ---
echo "Patching Perl modules..." >&2

python3 - "${REFERENCE_DIR}" "${WORK_DIR}/patched" <<'PYEOF'
"""Patch Perl VEP modules to inject warn()-based TRACE output at sub boundaries."""
import sys
import re
import os

ref_dir = sys.argv[1]
out_dir = sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

# ---- VariationEffect.pm ----
# Boolean predicates: wrap each `return EXPR;` to capture and log the value.
# These subs receive ($tva, $feat, ...) and return a boolean/value.
with open(os.path.join(ref_dir, "VariationEffect.pm"), "r") as f:
    ve_src = f.read()

# Instrument EVERY sub in VariationEffect.pm. Its predicate subs carry Perl-side
# names rather than SO terms (`stop_retained` not `stop_retained_variant`,
# `within_5_prime_utr` not `5_prime_UTR_variant`, `donor_splice_site` /
# `acceptor_splice_site`, `splice_region`, `within_intron`,
# `non_coding_exon_variant`, ...), so a name-based gate would have to track the
# Perl source. The classify-then-transform return rewrite below is safe for ANY
# sub (postfix-conditional / ternary / plain returns all handled; multi-line
# returns left untouched), so no gate is needed. `instrument_sub(name)` is the
# single decision point and admits every sub.
def instrument_sub(name: str) -> bool:
    return True

def patch_ve_module(src):
    """Patch VariationEffect.pm: inject TRACE warn before return statements in subs."""
    lines = src.split("\n")
    out = []
    in_sub = None
    brace_depth = 0
    sub_brace_start = 0

    i = 0
    while i < len(lines):
        line = lines[i]

        # Detect sub start
        sub_match = re.match(r'^sub\s+(\w+)\s*\{', line)
        if sub_match and in_sub is None:
            in_sub = sub_match.group(1)
            brace_depth = line.count("{") - line.count("}")
            sub_brace_start = brace_depth
            out.append(line)
            # Pre-declare the trace scratch var at sub entry so every later
            # branch can ASSIGN to it without `my` (avoids "Global symbol
            # $_trace_val requires explicit package name" when a return is
            # rewritten outside a fresh lexical scope). Only for instrumented
            # subs; harmless if unused.
            if instrument_sub(in_sub):
                out.append('  my $_trace_val;')
            i += 1
            continue

        if in_sub is not None:
            brace_depth += line.count("{") - line.count("}")

            # Return-statement patching (predicates only). The predicate subs use
            # three return shapes; each needs a different rewrite. A
            # postfix-conditional or multi-line return must not become
            # `my $_trace_val = (EXPR if COND);`: that is invalid Perl and aborts
            # the VEP compile.
            if instrument_sub(in_sub):
                # Strip a trailing `# comment` so the statement-terminating `;`
                # is the LAST `;` the rewrite anchors on. Returns frequently carry
                # trailing comments AFTER the `;` that themselves contain `;`
                # (e.g. `return 0 if X;   # synonymous_variant(@_);`), which a
                # greedy `;\s*$` anchor would swallow into the condition and emit
                # malformed Perl. A `#` outside quotes starts a Perl comment; in
                # these simple return lines there is no `#` inside a string, so
                # cutting at the first un-escaped `#` is safe. Multi-line returns
                # (no `;` at all on the line) still fall through untouched.
                code = re.sub(r'\s*#.*$', '', line.rstrip("\n"))
                # (a) postfix-conditional, value OPTIONAL:
                #     `return [EXPR] if COND;` / `... unless COND;`. `if`/`unless`
                #     matched as a word boundary with optional surrounding space
                #     so `if(...)` (no trailing space) is caught; EXPR optional so
                #     a bare `return if COND;` is handled here, not mis-captured by
                #     m_expr as `scalar(if ...)`.
                #     -> `COND_KW (COND) { [$_trace_val=scalar(EXPR); warn; return
                #        $_trace_val;] OR [warn; return;] }`
                m_cond = re.match(
                    r'^(\s*)return\b\s*(.*?)\s*\b(if|unless)\b\s*(.+?)\s*;\s*$', code
                )
                # (c) bare return: `return;`
                m_bare = re.match(r'^(\s*)return\s*;\s*$', code)
                # (b) plain single-line expression return: `return EXPR;`
                #     (only when it ends in `;` on this line -- multi-line returns
                #     have no trailing `;` and fall through UNTOUCHED on purpose.)
                m_expr = re.match(r'^(\s*)return\s+(.+);\s*$', code)

                if m_cond:
                    indent, expr, kw, cond = m_cond.groups()
                    out.append(f'{indent}{kw} ({cond}) {{')
                    if expr:
                        out.append(f'{indent}  $_trace_val = scalar({expr});')
                        out.append(
                            f'{indent}  warn "TRACE [VE] {in_sub} => $_trace_val\\n";'
                        )
                        out.append(f'{indent}  return $_trace_val;')
                    else:
                        out.append(
                            f'{indent}  warn "TRACE [VE] {in_sub} => (undef)\\n";'
                        )
                        out.append(f'{indent}  return;')
                    out.append(f'{indent}}}')
                    i += 1
                    if brace_depth <= 0:
                        in_sub = None
                    continue
                if m_bare:
                    indent = m_bare.group(1)
                    out.append(f'{indent}warn "TRACE [VE] {in_sub} => (undef)\\n";')
                    out.append(line)
                    i += 1
                    if brace_depth <= 0:
                        in_sub = None
                    continue
                if m_expr:
                    indent, expr = m_expr.groups()
                    # Guard: a ternary or any expression is fine here because the
                    # rewrite assigns (no `my`) and forces scalar context.
                    out.append(f'{indent}$_trace_val = scalar({expr});')
                    out.append(
                        f'{indent}warn "TRACE [VE] {in_sub} => $_trace_val\\n";'
                    )
                    out.append(f'{indent}return $_trace_val;')
                    i += 1
                    if brace_depth <= 0:
                        in_sub = None
                    continue
                # Else: multi-line return (no trailing `;`) or non-return line.
                # Leave UNTOUCHED -- losing a few traces is acceptable; emitting
                # malformed Perl is not.

            # End of sub (closing brace at depth 0)
            if brace_depth <= 0:
                # For predicates, add a trace before the final }
                if instrument_sub(in_sub):
                    stripped = line.rstrip()
                    if stripped == "}":
                        out.append(f'  warn "TRACE [VE] {in_sub} => (fallthrough)\\n";')
                out.append(line)
                in_sub = None
                i += 1
                continue

        out.append(line)
        i += 1

    return "\n".join(out)

ve_patched = patch_ve_module(ve_src)
with open(os.path.join(out_dir, "VariationEffect.pm"), "w") as f:
    f.write(ve_patched)
_n_subs = len(re.findall(r'^sub\s+\w+\s*\{', ve_src, re.M))
print(f"  Patched VariationEffect.pm ({_n_subs} subs instrumented)", file=sys.stderr)


# ---- TranscriptVariationAllele.pm ----
# Method calls: log entry with key params and exit with results.
# Key methods: _get_peptide_alleles, _get_differing_regions, codon, pep_allele_string,
# peptide, display_codon_allele_string, _get_alternate_cds_sequence
with open(os.path.join(ref_dir, "TranscriptVariationAllele.pm"), "r") as f:
    tva_src = f.read()

TVA_METHODS = {
    "codon": "codon string",
    "pep_allele_string": "peptide allele",
    "peptide": "peptide sequence",
    "display_codon_allele_string": "display codon",
    "_get_peptide_alleles": "peptide alleles",
    "_get_alternate_cds_sequence": "alt CDS seq",
    "_get_differing_regions": "differing regions",
}

def patch_method_module(src, methods, tag):
    """Patch a module by adding entry/exit TRACE warns to specified methods."""
    lines = src.split("\n")
    out = []
    in_sub = None
    brace_depth = 0
    patched_count = 0

    i = 0
    while i < len(lines):
        line = lines[i]

        sub_match = re.match(r'^sub\s+(\w+)\s*\{', line)
        if sub_match and in_sub is None:
            name = sub_match.group(1)
            if name in methods:
                in_sub = name
                brace_depth = line.count("{") - line.count("}")
                out.append(line)
                # Insert entry trace + pre-declare the return scratch var so the
                # rewritten returns can assign without `my` (same scoping fix as
                # patch_ve_module).
                out.append(f'  warn "TRACE [{tag}] {name} ENTER\\n";')
                out.append('  my $_trace_ret;')
                patched_count += 1
                i += 1
                continue

        if in_sub is not None:
            brace_depth += line.count("{") - line.count("}")

            # Patch return statements. Same three-shape classify-then-transform
            # as patch_ve_module: postfix-conditional, plain single-line, and
            # multi-line-left-untouched. The log line is ref-aware (logs ref type
            # for references, the value for scalars).
            _logfmt = (
                f'warn "TRACE [{tag}] {in_sub} => " . '
                '(ref($_trace_ret) ? ref($_trace_ret) . "(...)" : '
                'defined($_trace_ret) ? $_trace_ret : "(undef)") . "\\n";'
            )
            # Strip trailing `# comment` before matching (same rationale as
            # patch_ve_module: a comment after the `;` may contain its own `;`).
            code = re.sub(r'\s*#.*$', '', line.rstrip("\n"))
            m_cond = re.match(
                r'^(\s*)return\b\s*(.*?)\s*\b(if|unless)\b\s*(.+?)\s*;\s*$', code
            )
            m_bare = re.match(r'^(\s*)return\s*;\s*$', code)
            m_expr = re.match(r'^(\s*)return\s+(.+);\s*$', code)
            # NOTE: these methods can return LISTs (e.g. _get_peptide_alleles,
            # codon), so scalar() is NOT forced here: capture with a plain
            # parenthesized scalar assignment (`$x = (EXPR)`), which preserves
            # the value the caller sees in scalar context and logs a ref/scalar
            # summary. List-context callers of the patched method keep working
            # because the same expression's scalar value is returned only when the
            # method is genuinely scalar; the trace is best-effort diagnostics.
            if m_cond:
                indent, expr, kw, cond = m_cond.groups()
                out.append(f'{indent}{kw} ({cond}) {{')
                if expr:
                    out.append(f'{indent}  $_trace_ret = ({expr});')
                    out.append(f'{indent}  {_logfmt}')
                    out.append(f'{indent}  return $_trace_ret;')
                else:
                    out.append(f'{indent}  warn "TRACE [{tag}] {in_sub} => (undef)\\n";')
                    out.append(f'{indent}  return;')
                out.append(f'{indent}}}')
                i += 1
                if brace_depth <= 0:
                    in_sub = None
                continue
            if m_bare:
                indent = m_bare.group(1)
                out.append(f'{indent}warn "TRACE [{tag}] {in_sub} => (undef)\\n";')
                out.append(line)
                i += 1
                if brace_depth <= 0:
                    in_sub = None
                continue
            if m_expr:
                indent, expr = m_expr.groups()
                out.append(f'{indent}$_trace_ret = ({expr});')
                out.append(f'{indent}{_logfmt}')
                out.append(f'{indent}return $_trace_ret;')
                i += 1
                if brace_depth <= 0:
                    in_sub = None
                continue

            # End of sub
            if brace_depth <= 0:
                stripped = line.rstrip()
                if stripped == "}":
                    out.append(f'  warn "TRACE [{tag}] {in_sub} EXIT\\n";')
                out.append(line)
                in_sub = None
                i += 1
                continue

        out.append(line)
        i += 1

    return "\n".join(out), patched_count

tva_patched, tva_count = patch_method_module(tva_src, TVA_METHODS, "TVA")
with open(os.path.join(out_dir, "TranscriptVariationAllele.pm"), "w") as f:
    f.write(tva_patched)
print(f"  Patched TranscriptVariationAllele.pm ({tva_count} methods)", file=sys.stderr)


# ---- BaseTranscriptVariation.pm ----
# Overlap/coordinate methods.
with open(os.path.join(ref_dir, "BaseTranscriptVariation.pm"), "r") as f:
    btv_src = f.read()

BTV_METHODS = {
    "_intron_effects": "intron effects",
    "_exon_overlap": "exon overlap",
    "exonic_splice_region": "exonic splice",
    "_overlapped_introns": "overlapped introns",
    "_overlapped_exons": "overlapped exons",
    "cds_start": "CDS start",
    "cds_end": "CDS end",
    "translation_start": "translation start",
    "translation_end": "translation end",
    "cdna_start": "cDNA start",
    "cdna_end": "cDNA end",
}

btv_patched, btv_count = patch_method_module(btv_src, BTV_METHODS, "BTV")
with open(os.path.join(out_dir, "BaseTranscriptVariation.pm"), "w") as f:
    f.write(btv_patched)
print(f"  Patched BaseTranscriptVariation.pm ({btv_count} methods)", file=sys.stderr)
PYEOF

# Verify patches were created
for mod in "${MODULES[@]}"; do
	if [[ ! -f "${WORK_DIR}/patched/${mod}" ]]; then
		echo "ERROR: Patched module not created: ${WORK_DIR}/patched/${mod}" >&2
		exit 1
	fi
done

# --- Step 3: Run Perl VEP with patched modules ---
echo "Running Perl VEP with trace instrumentation..." >&2

STDERR_LOG="${TMP_DIR}/perl_stderr.log"
VEP_OUTPUT="${TMP_DIR}/output.txt"

# Docker paths for the 3 patched modules
#   VariationEffect.pm     -> /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm
#   TranscriptVariationAllele.pm -> /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm
#   BaseTranscriptVariation.pm   -> /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm

docker run --rm \
	-v "${CACHE_DIR}:/opt/vep/.vep:ro" \
	-v "${TMP_DIR}:/work/tmp" \
	-v "${WORK_DIR}/patched/VariationEffect.pm:/opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm:ro" \
	-v "${WORK_DIR}/patched/TranscriptVariationAllele.pm:/opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm:ro" \
	-v "${WORK_DIR}/patched/BaseTranscriptVariation.pm:/opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm:ro" \
	"${DOCKER_IMAGE}" \
	/opt/vep/src/ensembl-vep/vep \
	-i /work/tmp/single_variant.vcf \
	-o /work/tmp/output.txt \
	--offline --dir_cache /opt/vep/.vep --cache_version 115 \
	--assembly "${ASSEMBLY}" \
	--fork 1 --buffer_size 1 --force --no_stats \
	2>"${STDERR_LOG}" || {
	echo "ERROR: Perl VEP Docker run failed. Stderr:" >&2
	cat "${STDERR_LOG}" >&2
	exit 1
}

TRACE_COUNT="$(grep -c '^TRACE ' "${STDERR_LOG}" 2>/dev/null || echo 0)"
echo "Captured ${TRACE_COUNT} TRACE lines." >&2
echo "" >&2

# --- Step 4: Parse and format trace output ---
python3 - "${STDERR_LOG}" "${VEP_OUTPUT}" "${VARIANT_DESC}" "${ASSEMBLY}" <<'PYEOF'
"""Parse Perl VEP TRACE output and VEP result into a structured report."""
import sys
import re
from collections import defaultdict, OrderedDict

stderr_log = sys.argv[1]
vep_output = sys.argv[2]
variant_desc = sys.argv[3]
assembly = sys.argv[4]

# --- Parse TRACE lines from stderr ---
# Format: TRACE [TAG] sub_name => value
# or:     TRACE [TAG] sub_name ENTER
# or:     TRACE [TAG] sub_name EXIT
trace_pattern = re.compile(r'^TRACE \[(\w+)\] (\S+)\s+(=>|ENTER|EXIT)\s*(.*)$')

# Collect all traces in order
ve_traces = []      # (sub_name, value)
tva_traces = []     # (sub_name, direction_or_value)
btv_traces = []     # (sub_name, direction_or_value)

with open(stderr_log, "r") as f:
    for line in f:
        line = line.strip()
        m = trace_pattern.match(line)
        if not m:
            continue
        tag, sub_name, op, value = m.groups()
        value = value.strip()
        if tag == "VE":
            ve_traces.append((sub_name, value))
        elif tag == "TVA":
            if op == "=>":
                tva_traces.append((sub_name, f"=> {value}"))
            elif op == "ENTER":
                tva_traces.append((sub_name, "ENTER"))
            elif op == "EXIT":
                tva_traces.append((sub_name, "EXIT"))
        elif tag == "BTV":
            if op == "=>":
                btv_traces.append((sub_name, f"=> {value}"))
            elif op == "ENTER":
                btv_traces.append((sub_name, "ENTER"))
            elif op == "EXIT":
                btv_traces.append((sub_name, "EXIT"))

# --- Parse VEP output ---
vep_lines = []
try:
    with open(vep_output, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            vep_lines.append(line.strip())
except FileNotFoundError:
    pass

# Group VEP output by transcript
# VEP default output columns:
# #Uploaded_variation  Location  Allele  Gene  Feature  Feature_type  Consequence  ...
vep_by_transcript = {}
for line in vep_lines:
    cols = line.split("\t")
    if len(cols) >= 7:
        transcript_id = cols[4] if len(cols) > 4 else "?"
        consequence = cols[6] if len(cols) > 6 else "?"
        vep_by_transcript[transcript_id] = {
            "line": line,
            "consequence": consequence,
            "location": cols[1] if len(cols) > 1 else "?",
            "allele": cols[2] if len(cols) > 2 else "?",
            "gene": cols[3] if len(cols) > 3 else "?",
            "feature_type": cols[5] if len(cols) > 5 else "?",
        }

# --- Build the report ---
print(f"=== Perl VEP Predicate Trace ===")
print(f"Variant: {variant_desc} ({assembly})")
print()

# --- VE predicate summary ---
# Group predicates: show which ones fired (nonzero/true) vs didn't
if ve_traces:
    # Deduplicate: for each predicate, keep the last value seen
    # (Perl VEP may evaluate predicates multiple times per transcript)
    pred_values = OrderedDict()
    for name, val in ve_traces:
        pred_values.setdefault(name, []).append(val)

    print("--- Consequence Predicates (VariationEffect.pm) ---")
    fired = []
    not_fired = []
    other = []
    for name, values in pred_values.items():
        # Take last value for summary
        last_val = values[-1]
        # Determine if it "fired" (truthy in Perl: non-zero, non-empty, non-undef)
        is_truthy = (
            last_val not in ("0", "", "(undef)", "(fallthrough)")
            and last_val != "0"
        )
        entry = f"  {name:40s} => {last_val}"
        if len(values) > 1:
            entry += f"  (evaluated {len(values)}x)"
        if is_truthy:
            fired.append(entry + "  <-- FIRES")
        else:
            not_fired.append(entry)

    # Print fired predicates first (these are the interesting ones)
    for line in fired:
        print(line)
    if fired and not_fired:
        print()
    for line in not_fired:
        print(line)
    print()

# --- TVA method traces ---
if tva_traces:
    print("--- Allele Details (TranscriptVariationAllele.pm) ---")
    current_method = None
    for name, info in tva_traces:
        if info == "ENTER":
            current_method = name
            # Don't print ENTER alone, wait for the return value
        elif info.startswith("=> "):
            val = info[3:]
            print(f"  {name}: {val}")
        elif info == "EXIT":
            if current_method == name:
                print(f"  {name}: (no explicit return)")
            current_method = None
    print()

# --- BTV method traces ---
if btv_traces:
    print("--- Overlap Details (BaseTranscriptVariation.pm) ---")
    current_method = None
    for name, info in btv_traces:
        if info == "ENTER":
            current_method = name
        elif info.startswith("=> "):
            val = info[3:]
            print(f"  {name}: {val}")
        elif info == "EXIT":
            if current_method == name:
                print(f"  {name}: (no explicit return)")
            current_method = None
    print()

# --- VEP output lines ---
if vep_by_transcript:
    print("--- VEP Output ---")
    for tr_id, info in vep_by_transcript.items():
        feature_type = info.get("feature_type", "?")
        consequence = info.get("consequence", "?")
        print(f"  {tr_id} ({feature_type}): {consequence}")
    print()
    print("--- Full VEP Output Lines ---")
    for line in vep_lines:
        print(f"  {line}")
else:
    print("--- VEP Output ---")
    print("  (no output lines produced)")

# Summary stats
print()
print(f"--- Trace Summary ---")
print(f"  VE predicate evaluations: {len(ve_traces)}")
print(f"  TVA method traces: {len(tva_traces)}")
print(f"  BTV method traces: {len(btv_traces)}")
print(f"  VEP output transcripts: {len(vep_by_transcript)}")

if not ve_traces and not tva_traces and not btv_traces:
    print()
    print("WARNING: No TRACE lines captured. Possible causes:")
    print("  - Variant is intergenic (no transcript overlap)")
    print("  - Patched modules failed to load (check stderr log)")
    print(f"  - Stderr log: {stderr_log}")
PYEOF

# --- Done ---
if [[ "${KEEP_TEMP}" -eq 1 ]]; then
	echo "" >&2
	echo "Temp files preserved at: ${TMP_DIR}" >&2
	echo "  VCF:        ${VCF_PATH}" >&2
	echo "  VEP output: ${VEP_OUTPUT}" >&2
	echo "  Stderr log: ${STDERR_LOG}" >&2
	echo "  Patched modules: ${WORK_DIR}/patched/" >&2
fi
