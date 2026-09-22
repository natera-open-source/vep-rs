#!/usr/bin/env bash
set -euo pipefail
#
# Stage plugin data files for vep-rs plugin testing, per assembly.
#
# TWO SOURCES. Neither downloads from the upstream databases unattended.
#
# `--source upstream` (the default) PRINTS, per applicable plugin, the upstream
# URL and the preparation each file needs (header fixes, per-assembly tabix
# indexes, decompression), then reports what is on disk. You run those steps;
# the script does not fetch tens of GB on its own.
#
# `--source s3` pulls an already-prepared copy from an S3 bucket you supply via
# --bucket or $VEP_BENCHMARK_BUCKET, and errors out when neither names one. That
# is worth setting up when several machines need the same tens of GB, or when a
# machine's local storage does not survive a reboot. Seed the bucket once from
# the `--source upstream` steps, then upload the prepared tree laid out as
#   <bucket>/plugin_data/<scope>/<relative path>
#   <bucket>/plugin_data/MANIFEST.json
# where <scope> is grch37, grch38, or assembly_agnostic (LoFtool, pLI, GWAS),
# and <relative path> is the path plugin_data_layout.sh assigns to the file, the
# same path it takes under --data-dir. MANIFEST.json is optional. When present, a
# staged file whose byte size differs from its entry is re-fetched as a truncated
# transfer. Its schema, one entry per prepared file:
#   {"scopes": {"grch37": [{"path": "<relative path>", "bytes": <int>,
#                           "sha256": "<hex>"}, ...],
#               "grch38": [...], "assembly_agnostic": [...]}}
# `bytes` is the field this script reads; `sha256` is recorded for auditing.
#
# Assemblies: each database defines which builds it supports. See the SUPPORT
# MATRIX below; this script refuses combinations the upstream data does not offer
# rather than staging something that will silently annotate nothing.
#
# Usage:
#   scripts/data/setup_plugin_data.sh --assembly GRCh38
#   scripts/data/setup_plugin_data.sh --assembly GRCh37 --plugins CADD,REVEL
#   scripts/data/setup_plugin_data.sh --assembly GRCh38 --tier full
#   scripts/data/setup_plugin_data.sh --assembly GRCh38 --source s3 --bucket <name>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/data/plugin_data_layout.sh
source "${SCRIPT_DIR}/plugin_data_layout.sh"
VEP_RS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

ASSEMBLY=""
DATA_DIR=""
TIER="quick"
PLUGINS="" # empty = all supported for the assembly
SOURCE="upstream"
S3_BUCKET="${VEP_BENCHMARK_BUCKET:-}"
S3_PREFIX="plugin_data"
AWS_PROFILE_ARG=()
DRY_RUN=0

# SUPPORT MATRIX: which assemblies each database actually publishes.
#
# A combination absent here is NOT a vep-rs gap -- the database does not offer
# it. Staging a substitute would annotate from the wrong genome build.
# CADD           GRCh37 + GRCh38   separate per-assembly trees
# dbscSNV        GRCh37 + GRCh38   pre-split per assembly
# AlphaMissense  GRCh37 + GRCh38   hg19 + hg38 files
# REVEL          GRCh37 + GRCh38   ONE file, two coordinate columns; per-assembly index
# gnomADc        GRCh37 + GRCh38   v2.1 genomes (37) / v4.0 exomes (38)
# LoFtool        both              gene-symbol keyed, assembly-independent
# pLI            both              gene keyed, assembly-independent
# LoFTEE         GRCh37 + GRCh38   STRUCTURALLY different bundles (tabix TSV vs bigWig)
# GWAS           GRCh38 only       the catalog publishes GRCh38 coordinates
# SpliceAI       GRCh38 SNV only   GRCh37 + all indels are Illumina-login-gated
# dbNSFP         UNSUPPORTED       commercial license required for clinical use

ALL_PLUGINS=(CADD REVEL SpliceAI gnomADc AlphaMissense dbscSNV GWAS LoFtool LoFTEE pLI)

# Plugins that only support one assembly.
GRCH38_ONLY_PLUGINS=(GWAS SpliceAI)

# Explicitly unsupported, with the reason surfaced to the user.
UNSUPPORTED_PLUGINS=(dbNSFP)

usage() {
	cat <<'EOF'
Stage plugin data files for vep-rs plugin testing, per assembly.

Usage:
  scripts/data/setup_plugin_data.sh --assembly <GRCh37|GRCh38> [options]

Required:
  --assembly <name>     GRCh37 or GRCh38. Determines which files are staged and
                        which plugins are applicable.

Options:
  --data-dir <path>     Destination (default: ~/.vep/plugin_data/<assembly>/<tier>)
  --source <upstream|s3>  Where the files come from (default: upstream)
                          upstream: print each database's URL and the
                                    preparation its file needs, for you to run;
                                    nothing is downloaded
                          s3:       sync an already-prepared copy from a bucket
                                    you supply (fast, same-region); needs
                                    --bucket or $VEP_BENCHMARK_BUCKET. The bucket
                                    layout is in this script's header.
  --tier <quick|full>   quick: chr21 subsets  full: genome-wide (default: quick)
  --plugins <list>      Comma-separated (default: all supported for the assembly)
                          CADD,REVEL,SpliceAI,gnomADc,AlphaMissense,dbscSNV,
                          GWAS,LoFtool,LoFTEE,pLI
  --bucket <name>       Your S3 bucket holding an already-prepared copy. Required
                        for --source s3; defaults to $VEP_BENCHMARK_BUCKET.
  --profile <name>      AWS profile for S3 access
  --dry-run             Print what would be staged, change nothing
  --help                Print this help

Not supported:
  dbNSFP    Its license requires a paid commercial license for clinical or
            commercial use, so it is neither distributed nor validated against.
  SpliceAI  GRCh38 SNV scores only (open, from the Ensembl MANE release).
            GRCh37 scores and the indel files for both assemblies are
            distributed solely via Illumina BaseSpace, which requires an
            account signup.
  GWAS      GRCh38 only, because the NHGRI-EBI catalog publishes GRCh38
            coordinates. Not a vep-rs limitation.

The script is idempotent: existing files are skipped, and tabix indexes are
verified for every file that needs one.
EOF
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--assembly)
		ASSEMBLY="$2"
		shift 2
		;;
	--data-dir)
		DATA_DIR="$2"
		shift 2
		;;
	--source)
		SOURCE="$2"
		shift 2
		;;
	--tier)
		TIER="$2"
		shift 2
		;;
	--plugins)
		PLUGINS="$2"
		shift 2
		;;
	--bucket)
		S3_BUCKET="$2"
		shift 2
		;;
	--profile)
		AWS_PROFILE_ARG=(--profile "$2")
		shift 2
		;;
	--dry-run)
		DRY_RUN=1
		shift 1
		;;
	--help | -h)
		usage
		exit 0
		;;
	*)
		echo "Unknown argument: $1" >&2
		usage >&2
		exit 1
		;;
	esac
done

# Validate arguments

if [[ -z "${ASSEMBLY}" ]]; then
	echo "ERROR: --assembly is required (GRCh37 or GRCh38)." >&2
	echo "  It must be explicit: the plugin data directory and the concordance" >&2
	echo "  harness pick their assembly independently, so an unstated assembly can" >&2
	echo "  read one assembly's data while annotating against the other, which" >&2
	echo "  annotates nothing." >&2
	exit 1
fi

ASSEMBLY_LC="$(echo "${ASSEMBLY}" | tr '[:upper:]' '[:lower:]')"
case "${ASSEMBLY_LC}" in
grch37 | hg19)
	ASSEMBLY="GRCh37"
	ASSEMBLY_LC="grch37"
	UCSC_BUILD="hg19"
	;;
grch38 | hg38)
	ASSEMBLY="GRCh38"
	ASSEMBLY_LC="grch38"
	UCSC_BUILD="hg38"
	;;
*)
	echo "ERROR: --assembly must be GRCh37 or GRCh38 (got '${ASSEMBLY}')" >&2
	exit 1
	;;
esac

if [[ "${TIER}" != "quick" && "${TIER}" != "full" ]]; then
	echo "ERROR: --tier must be quick or full" >&2
	exit 1
fi

if [[ "${SOURCE}" != "s3" && "${SOURCE}" != "upstream" ]]; then
	echo "ERROR: --source must be upstream or s3" >&2
	exit 1
fi

if [[ "${SOURCE}" == "s3" && "${DRY_RUN}" -eq 0 && -z "${S3_BUCKET}" ]]; then
	echo "ERROR: --source s3 needs a bucket: pass --bucket <name> or set \$VEP_BENCHMARK_BUCKET." >&2
	exit 1
fi

if [[ -z "${DATA_DIR}" ]]; then
	DATA_DIR="${HOME}/.vep/plugin_data/${ASSEMBLY_LC}/${TIER}"
fi

# Resolve the plugin list

if [[ -n "${PLUGINS}" ]]; then
	IFS=',' read -ra SELECTED_PLUGINS <<<"${PLUGINS}"
else
	SELECTED_PLUGINS=("${ALL_PLUGINS[@]}")
fi

in_list() {
	local needle="$1"
	shift
	local item
	for item in "$@"; do
		[[ "${item}" == "${needle}" ]] && return 0
	done
	return 1
}

# Drop (or reject) plugins the requested assembly cannot support.
APPLICABLE_PLUGINS=()
for plugin in "${SELECTED_PLUGINS[@]}"; do
	if in_list "${plugin}" "${UNSUPPORTED_PLUGINS[@]}"; then
		echo "SKIP  ${plugin}: not supported (commercial license required for clinical/commercial use)" >&2
		continue
	fi
	if ! in_list "${plugin}" "${ALL_PLUGINS[@]}"; then
		echo "ERROR: unknown plugin '${plugin}'. Known: ${ALL_PLUGINS[*]}" >&2
		exit 1
	fi
	if in_list "${plugin}" "${GRCH38_ONLY_PLUGINS[@]}" && [[ "${ASSEMBLY_LC}" != "grch38" ]]; then
		case "${plugin}" in
		GWAS)
			echo "SKIP  GWAS: GRCh38 only (the NHGRI-EBI catalog publishes GRCh38 coordinates)" >&2
			;;
		SpliceAI)
			echo "SKIP  SpliceAI: GRCh38 SNV only (GRCh37 + indels are behind an Illumina BaseSpace signup)" >&2
			;;
		esac
		continue
	fi
	APPLICABLE_PLUGINS+=("${plugin}")
done

if [[ ${#APPLICABLE_PLUGINS[@]} -eq 0 ]]; then
	echo "ERROR: no applicable plugins for ${ASSEMBLY}" >&2
	exit 1
fi

echo "Plugin data staging"
echo "  assembly:  ${ASSEMBLY}"
echo "  tier:      ${TIER}"
echo "  source:    ${SOURCE}"
echo "  data dir:  ${DATA_DIR}"
echo "  plugins:   ${APPLICABLE_PLUGINS[*]}"
echo ""

if [[ "${DRY_RUN}" -eq 1 ]]; then
	echo "(dry run: nothing will be downloaded or written)"
	exit 0
fi

mkdir -p "${DATA_DIR}"

# Helpers

# Every file this script stages, as "<plugin>|<relative path>".
#
# Paths are IDENTICAL in the bucket and on disk, so a sync needs no path
# translation.
staged_files_for() {
	# One table, shared with run_concordance.sh: scripts/data/plugin_data_layout.sh.
	plugin_staged_files "$1" "${ASSEMBLY_LC}"
}

# Which S3 prefix holds a plugin's data: per-assembly, or shared.
s3_scope_for() {
	case "$1" in
	LoFtool | pLI | GWAS) echo "assembly_agnostic" ;;
	*) echo "${ASSEMBLY_LC}" ;;
	esac
}

# Local copy of the bucket's MANIFEST.json, fetched once per run when staging
# from S3. Empty when unavailable (e.g. --source upstream, or no aws CLI).
MANIFEST_FILE=""

# Fetch the bucket's manifest (schema in the header): bytes + sha256 per file.
fetch_manifest() {
	local dest="${DATA_DIR}/.plugin_data_manifest.json"
	if aws s3 cp "${AWS_PROFILE_ARG[@]+"${AWS_PROFILE_ARG[@]}"}" \
		"s3://${S3_BUCKET}/${S3_PREFIX}/MANIFEST.json" "${dest}" --only-show-errors; then
		MANIFEST_FILE="${dest}"
	else
		echo "    [warn] no MANIFEST.json in the bucket; staged sizes cannot be verified" >&2
	fi
}

# Expected byte size for a staged file, from the manifest. Empty when unknown.
manifest_bytes_for() {
	local scope="$1" rel="$2"
	[[ -z "${MANIFEST_FILE}" ]] && return 0
	MANIFEST_PATH="${MANIFEST_FILE}" MANIFEST_SCOPE="${scope}" MANIFEST_REL="${rel}" python3 -c '
import json, os, sys
try:
    with open(os.environ["MANIFEST_PATH"], encoding="utf-8") as fh:
        manifest = json.load(fh)
except (OSError, ValueError):
    sys.exit(0)
entries = manifest.get("scopes", {}).get(os.environ["MANIFEST_SCOPE"], [])
want = os.environ["MANIFEST_REL"]
for entry in entries:
    if entry.get("path") == want and isinstance(entry.get("bytes"), int):
        print(entry["bytes"])
        break
' 2>/dev/null
}

# A staged file whose size does not match the manifest is a TRUNCATED transfer,
# not a staged file. Existence alone cannot tell the two apart, and an interrupted
# transfer is routine on databases this large: a half-written multi-gigabyte file
# otherwise reports READY at exit 0 and the plugin then annotates a fraction of
# its rows. Compare against the recorded size, never against zero.
verify_staged_size() {
	local scope="$1" rel="$2" file="$3" expected actual
	expected="$(manifest_bytes_for "${scope}" "${rel}")"
	[[ -z "${expected}" ]] && return 0
	actual="$(wc -c <"${file}" | tr -d '[:space:]')"
	if [[ "${actual}" != "${expected}" ]]; then
		echo "    [warn] size mismatch: ${rel} is ${actual} bytes, manifest says ${expected}" >&2
		return 1
	fi
	return 0
}

verify_tabix_index() {
	local file="$1"
	if [[ ! -f "${file}" ]]; then
		return 1
	fi
	# Only index-bearing formats need a .tbi; flat gene-keyed files do not.
	case "${file}" in
	*.txt | *.tsv | *.fa | *.fai | *.bw) return 0 ;;
	esac
	if [[ -f "${file}.tbi" ]]; then
		return 0
	fi
	echo "    [warn] missing tabix index: ${file}.tbi" >&2
	return 1
}

# S3 staging

stage_from_s3() {
	if ! command -v aws >/dev/null 2>&1; then
		echo "ERROR: aws CLI not found; needed for --source s3" >&2
		exit 1
	fi

	fetch_manifest

	local plugin scope rel dest src
	for plugin in "${APPLICABLE_PLUGINS[@]}"; do
		echo "[${plugin}]"
		scope="$(s3_scope_for "${plugin}")"
		while IFS= read -r rel; do
			[[ -z "${rel}" ]] && continue
			dest="${DATA_DIR}/${rel}"
			# Re-fetch a file whose size does not match the manifest: it is a
			# truncated leftover from an interrupted run, not a staged file.
			if [[ -f "${dest}" ]] && verify_staged_size "${scope}" "${rel}" "${dest}"; then
				echo "    [skip] ${rel}"
				continue
			fi
			if [[ -f "${dest}" ]]; then
				echo "    [re-sync] ${rel} (size mismatch)"
			fi
			mkdir -p "$(dirname "${dest}")"
			src="s3://${S3_BUCKET}/${S3_PREFIX}/${scope}/${rel}"
			echo "    [sync] ${rel}"
			if ! aws s3 cp "${AWS_PROFILE_ARG[@]+"${AWS_PROFILE_ARG[@]}"}" \
				"${src}" "${dest}" --only-show-errors; then
				echo "    [ERROR] not in the bucket: ${src}" >&2
				echo "            Re-seed it with --source upstream, then upload." >&2
				return 1
			fi
		done < <(staged_files_for "${plugin}")
		echo ""
	done
}

# Upstream staging (the default; also how a bucket is seeded)

CADD_BASE="https://kircherlab.bihealth.org/download/CADD/v1.7/${ASSEMBLY}"
# The kircherlab host is the one the CADD site advertises; the older
# krishna.gs.washington.edu host also serves these files.

upstream_note() {
	cat <<EOF
    [manual] ${1}
EOF
}

stage_from_upstream() {
	echo "Re-seeding from upstream hosts. Prepared output goes to ${DATA_DIR}."
	if [[ -n "${S3_BUCKET}" ]]; then
		echo "To reuse it with --source s3, upload it to s3://${S3_BUCKET}/${S3_PREFIX}/ when done."
	fi
	echo ""

	for cmd in curl bgzip tabix; do
		if ! command -v "${cmd}" >/dev/null 2>&1; then
			echo "ERROR: missing required tool: ${cmd} (install htslib for bgzip/tabix)" >&2
			exit 1
		fi
	done

	local plugin
	for plugin in "${APPLICABLE_PLUGINS[@]}"; do
		echo "[${plugin}]"
		case "${plugin}" in
		CADD)
			echo "    SNV:   ${CADD_BASE}/whole_genome_SNVs.tsv.gz (+ .tbi)"
			if [[ "${ASSEMBLY_LC}" == "grch37" ]]; then
				echo "    indel: ${CADD_BASE}/gnomad.genomes-exomes.r4.0.indel.tsv.gz (+ .tbi)"
				upstream_note "GRCh37 note: the SNV file's REMOTE tabix index is corrupt (bgzf read fails at every offset), so a streaming 'tabix <url> 21' slice does not work. Download the full file, or resolve chr21's offsets from the .tbi and range-GET them."
			else
				echo "    indel: ${CADD_BASE}/gnomad.genomes.r4.0.indel.tsv.gz (+ .tbi)"
				upstream_note "GRCh38: remote tabix works, so 'tabix <url> 21 | bgzip' yields a chr21 subset without the full download."
			fi
			;;
		dbscSNV)
			echo "    ${CADD_BASE}/dbscSNV1.1_${ASSEMBLY}.txt.gz (+ .tbi)"
			upstream_note "PREPARATION REQUIRED: the published header has NO leading '#' (GRCh37 'chr pos ref alt ada_score rf_score'; GRCh38 'hg38_chr hg38_pos ...'). vep-rs tolerates that, but Perl's dbscSNV.pm does not, so store it with '#' prepended: zcat f | sed '1s/^/#/' | bgzip > out; tabix -s 1 -b 2 -e 2 out. Files are also CRLF."
			;;
		AlphaMissense)
			echo "    https://zenodo.org/records/8208688/files/AlphaMissense_${UCSC_BUILD}.tsv.gz"
			upstream_note "No .tbi is shipped: build one with 'tabix -s 1 -b 2 -e 2'. Contigs are UCSC-style ('chr21') on BOTH builds. License CC-BY-NC-SA-4.0."
			;;
		REVEL)
			echo "    https://zenodo.org/records/7072866 (revel-v1.3_all_chromosomes.zip,"
			echo "      or revel-v1.3_segments_chrom_21.zip for the chr21 tier -- 6.8 MB vs 667 MB)"
			if [[ "${ASSEMBLY_LC}" == "grch37" ]]; then
				upstream_note "Index on the hg19 column: tabix -f -s 1 -b 2 -e 2"
			else
				upstream_note "Index on the grch38 column, which needs a re-sort first: keep the header, drop rows whose grch38_pos is '.', sort -k1,1 -k3,3n, then tabix -f -s 1 -b 3 -e 3. One REVEL file serves both assemblies via different indexes."
			fi
			;;
		gnomADc)
			if [[ "${ASSEMBLY_LC}" == "grch37" ]]; then
				echo "    s3://gnomad-public-us-east-1/release/2.1/coverage/genomes/gnomad.genomes.coverage.summary.tsv.bgz"
			else
				echo "    s3://gnomad-public-us-east-1/release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.summary.tsv.bgz"
			fi
			upstream_note "Public AWS Open Data bucket: fetch with 'aws s3 cp --no-sign-request'. No .tbi is published, so build one. Perl's gnomADc.pm needs a '#'-prefixed header: sed '1s/.*/#&/'."
			;;
		SpliceAI)
			echo "    https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/"
			echo "      spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz (+ .tbi)"
			upstream_note "GRCh38 SNV only, and MANE-restricted, so coverage is narrower than Illumina's full callset. Remote tabix works. Indel scores and all GRCh37 scores require an Illumina BaseSpace account and are NOT supported."
			;;
		GWAS)
			echo "    https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
			echo "      gwas-catalog-associations_ontology-annotated-full.zip"
			upstream_note "Unzip to a plain TSV: Perl's GWAS.pm needs it uncompressed, and the header row starts with 'DATE ADDED TO CATALOG'. GRCh38 coordinates only."
			;;
		LoFtool)
			echo "    https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/115/LoFtool_scores.txt"
			upstream_note "Fetched directly from GitHub, not from a VEP_plugins clone: vep-rs must not depend on the Perl plugin repo being checked out. Gene-symbol keyed, so it serves both assemblies."
			;;
		pLI)
			echo "    https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/115/pLI_values.txt"
			upstream_note "Same GitHub-direct rationale as LoFtool. Gene keyed, serves both assemblies."
			;;
		LoFTEE)
			echo "    https://personal.broadinstitute.org/konradk/loftee_data/${ASSEMBLY}/"
			if [[ "${ASSEMBLY_LC}" == "grch37" ]]; then
				echo "      GERP_scores.final.sorted.txt.gz (+ .tbi), human_ancestor.fa.gz (+ .fai/.gzi)"
				upstream_note "GRCh37 (master branch) ships GERP as a tabix per-base TSV, NOT a bigWig. Pass it as gerp_tabix=. The END_TRUNC cutoff is +180 here, vs -58 on grch38."
			else
				echo "      gerp_conservation_scores.homo_sapiens.GRCh38.bw, human_ancestor.fa.gz (+ .fai/.gzi)"
				upstream_note "GRCh38 (grch38 branch) ships GERP as a bigWig. Pass it as gerp_bigwig=."
			fi
			upstream_note "Stage both forms: the Perl plugin reads human_ancestor.fa.gz with its .fai and .gzi, vep-rs the plain human_ancestor.fa, so keep the upstream .gz trio and decompress a copy, re-indexed with samtools faidx."
			;;
		esac
		echo ""
	done

	echo "NOTE: upstream staging prints the sources and required preparation rather"
	echo "      than fetching tens of GB unattended. Run the steps and verify. To"
	echo "      reuse the result, upload it to a bucket of your own and pass"
	echo "      --bucket <name> --source s3 on later runs."
}

case "${SOURCE}" in
s3) stage_from_s3 ;;
upstream) stage_from_upstream ;;
esac

# Summary: report per-plugin readiness from what is actually on disk

echo "=== Staging summary (${ASSEMBLY}, ${TIER}) ==="
echo ""
MISSING_ANY=0
for plugin in "${APPLICABLE_PLUGINS[@]}"; do
	all_present=1
	summary_scope="$(s3_scope_for "${plugin}")"
	while IFS= read -r rel; do
		[[ -z "${rel}" ]] && continue
		if [[ ! -f "${DATA_DIR}/${rel}" ]]; then
			all_present=0
			break
		fi
		verify_tabix_index "${DATA_DIR}/${rel}" || all_present=0
		verify_staged_size "${summary_scope}" "${rel}" "${DATA_DIR}/${rel}" || all_present=0
	done < <(staged_files_for "${plugin}")

	if [[ "${all_present}" -eq 1 ]]; then
		echo "  ${plugin}: READY"
	else
		echo "  ${plugin}: MISSING"
		MISSING_ANY=1
	fi
done
echo ""
echo "Data directory: ${DATA_DIR}"
echo ""
echo "Run plugin concordance with (the harness requires the Perl cache root, a vep-rs JSON cache and a FASTA):"
echo "  scripts/concordance/run_concordance.sh --assembly ${ASSEMBLY} \\"
echo "    --perl-cache-dir <dir above homo_sapiens/> --json-cache-dir <json cache dir> --fasta <genome.fa> \\"
echo "    --plugins $(
	IFS=','
	echo "${APPLICABLE_PLUGINS[*]}"
) \\"
echo "    --plugin-data-dir ${DATA_DIR}"

# A partially-staged directory must not look like success: every later gate
# assumes these files exist.
if [[ "${MISSING_ANY}" -eq 1 ]]; then
	echo ""
	echo "ERROR: one or more plugins are not fully staged (see MISSING above)." >&2
	exit 1
fi
