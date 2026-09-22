#!/usr/bin/env bash
#
# run_clone_measurement.sh -- the measurement entrypoint for a single engine.
#
# Runs the full 10-suite measurement (6 SNP/indel + 4 structural-variant) for ONE
# engine (vep-rs | fastvep | perl) against an already-populated data tree,
# emitting per-suite concordance reports (compare_vep_outputs.py for the
# SNP/indel suites, compare_sv_concordance.py for the SV suites) and a
# machine-readable wall_times.csv.
#
# Scope: measurement only. It does not fetch inputs, build binaries, or manage
# hosts. Point it at a data tree you have already populated, and give it the
# engine binary path or the Perl Docker image reference.
#
# The exact command lines this script embeds:
#
#   vep-rs   SNP/indel: vep -i <canonical.vcf.gz> -o <out> --offline \
#                  --json_cache <cache dir> --fasta <fa> --fork 16 \
#                  --buffer_size 5000 --force --no_stats --quiet
#   vep-rs   SV: same, per-VCF loop, same --fork as SNP/indel,
#                  plus --species homo_sapiens --assembly <ASM>
#   fastvep  SNP/indel: RAYON_NUM_THREADS=16 fastvep annotate --gff3 <g> \
#                  --fasta <fa> --transcript-cache <cache.bin> \
#                  -i <vcf> -o <raw> --output-format tab
#                then fastvep_normalize.py <raw> <output.txt>
#                then compare_vep_outputs.py ... --glob 'output.txt'
#   fastvep  SV: same annotate form, per-VCF loop; concordance is best-effort
#                  and is skipped when the loop produces no non-empty output.
#   perl   both: docker run ... vep -i <in> -o <out> --offline --cache \
#                  --dir_cache /work/cache --cache_version <ver> --fork 16 \
#                  --fasta /work/ref.fa --buffer_size 5000 --no_headers \
#                  --no_stats --quiet --force_overwrite
#                Perl is the ground truth, so it has no concordance branch and
#                its scratch output is deleted. It contributes wall time only.
#
# Expected --data-dir layout. Every path below is overridable by the matching
# VEP_* environment variable; these are only the defaults relative to --data-dir.
#
#   vcf/<asm>/all_variants_canonical/<base>.canonical.vcf.gz   SNP/indel inputs
#   ground_truth/perl/snp_indel/<suite>/output.txt             SNP/indel truth
#   ground_truth/perl/sv_per_vcf/canonical_inputs/<asm>/       SV inputs
#   ground_truth/perl/sv_per_vcf/<asm>/                        SV truth
#   caches/vep-rs/<asm>/                                       vep-rs JSON cache
#   caches/fastvep/<asm>.bin                                   fastVEP cache
#   caches/perl/vep-cache/                                     Perl Storable cache
#   reference/<asm>/<genome>.fa[.fai]                           reference FASTA
#   reference/<asm>/<annotation>.gff3                           GFF3 annotation
#
# "canonical" means Ensembl-style contig names. Raw all_variants/ inputs score
# F1=0 on the GRCh38 chr21 cells, whose contigs are `chr21` while the ground
# truth and the caches expect `21`; canonical inputs fix that symmetrically for
# every engine. Generate them with `prepare_benchmark_vcfs.py`; the vep-rs and
# fastVEP engines read them as given, and the Perl engine derives its own
# canonical, sites-only copy from all_variants/ with bcftools at run time.
#
# The two ground-truth trees are WRITTEN by generate_reference_output.sh in this
# directory, which runs the Perl invocation above over the same inputs (SV at the
# ground truth's own --fork 4). This script only reads them.
#
# Usage:
#   scripts/concordance/run_clone_measurement.sh \
#     --engine vep-rs --vep-binary /path/to/vep \
#     --data-dir <data> --output-dir <data>/work \
#     [--suites all|snp-indel|sv|"s01 s04"] [--walltime-only] [--skip-fields] [--fork 16] \
#     [--run-id <id>] [--host-id <id>]
#
#   scripts/concordance/run_clone_measurement.sh --engine fastvep \
#     --fastvep-binary /path/to/fastvep --data-dir <data> --output-dir <data>/work
#
#   scripts/concordance/run_clone_measurement.sh --engine perl \
#     --perl-image ensemblorg/ensembl-vep:release_115.2 \
#     --data-dir <data> --output-dir <data>/work
set -euo pipefail

# Locate the repo so the compare_*.py scripts resolve relative to THIS file
# (this script lives at scripts/concordance/, so the repo root is ../..).
# The harness scripts ship beside this one in the same checkout, so nothing else
# has to be located.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Defaults
ENGINE=""
VEP_BINARY=""
FASTVEP_BINARY=""
PERL_IMAGE="ensemblorg/ensembl-vep:release_115.2"
DATA_DIR=""
OUTPUT_DIR=""
SUITES_FILTER="all"
WALLTIME_ONLY="false"
SKIP_FIELDS="false"
# Default parallelism. 16 matches the --fork value the reported measurements
# used; override with --fork to match your own core count.
FORK="${VEP_FORK:-16}"
BUFFER_SIZE="${VEP_BUFFER_SIZE:-5000}"
RUN_ID="manual"
HOST_ID="local"
REPLICATE=""
CONTAINER_PROBE=0

usage() {
	cat <<'EOF'
Run the full 10-suite measurement for ONE engine on a populated data tree.

Required:
  --engine {vep-rs|fastvep|perl}
  --data-dir <path>      Populated data tree (vcf/, caches/, reference/, ground_truth/)
  --output-dir <path>    Where to write per-suite results + wall_times.csv

Engine binary (one, matching --engine):
  --vep-binary <path>      vep-rs binary             (engine vep-rs)
  --fastvep-binary <path>  fastVEP binary            (engine fastvep)
  --perl-image <ref>       Perl VEP Docker image      (engine perl; default
                           ensemblorg/ensembl-vep:release_115.2)

Optional:
  --suites <filter>      all (default) | snp-indel | sv | space-separated "s01 s04"
  --walltime-only        Skip concordance compare / GT-sync / normalize; time
                         the annotation only (annotation output still written so
                         the real I/O happens). Default: full concordance + timing.
  --skip-fields          Skip the per-column field comparator (report/fields.json)
                         on the vep-rs SNP/indel cells; F1 is computed either way.
                         For runs that need only F1 and wall time. Default: run it.
  --fork <N>             SNP/indel parallelism (default 16; fastVEP uses
                         RAYON_NUM_THREADS=N)
  --buffer-size <N>      Engine buffer size (default 5000)
  --run-id <id>          Provenance: identifier for this measurement campaign,
                         written to wall_times.csv
  --host-id <id>         Provenance: identifier for the machine running this
                         invocation, written to wall_times.csv
  --replicate <n>        Provenance: 1..N index of this run within the campaign
                         (defaults to the host id on a single run)
  --container-probe <N>  Perl engine only: after the suites, time N bare starts of
                         the Perl image (`docker run --rm <image> true`, no
                         annotation) under /usr/bin/time -v and write them to
                         <output-dir>/container_start.csv (default 0: none). Every
                         Perl cell pays one container start per invocation; the
                         probe measures that cost on its own.
EOF
	exit 1
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--engine)
		ENGINE="$2"
		shift 2
		;;
	--vep-binary)
		VEP_BINARY="$2"
		shift 2
		;;
	--fastvep-binary)
		FASTVEP_BINARY="$2"
		shift 2
		;;
	--perl-image)
		PERL_IMAGE="$2"
		shift 2
		;;
	--data-dir)
		DATA_DIR="$2"
		shift 2
		;;
	--output-dir)
		OUTPUT_DIR="$2"
		shift 2
		;;
	--suites)
		SUITES_FILTER="$2"
		shift 2
		;;
	--walltime-only)
		WALLTIME_ONLY="true"
		shift
		;;
	--skip-fields)
		SKIP_FIELDS="true"
		shift
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--buffer-size)
		BUFFER_SIZE="$2"
		shift 2
		;;
	--run-id)
		RUN_ID="$2"
		shift 2
		;;
	--host-id)
		HOST_ID="$2"
		shift 2
		;;
	--replicate)
		REPLICATE="$2"
		shift 2
		;;
	--container-probe)
		CONTAINER_PROBE="$2"
		shift 2
		;;
	-h | --help) usage ;;
	*)
		echo "ERROR: unknown argument: $1" >&2
		usage
		;;
	esac
done

# Validate args
[[ "$CONTAINER_PROBE" =~ ^[0-9]+$ ]] || {
	echo "ERROR: --container-probe takes a non-negative integer (got '$CONTAINER_PROBE')" >&2
	usage
}
[[ "$CONTAINER_PROBE" -eq 0 || "$ENGINE" == "perl" ]] || {
	echo "ERROR: --container-probe applies to --engine perl only (got engine '$ENGINE')" >&2
	usage
}
case "$ENGINE" in
vep-rs | fastvep | perl) ;;
"")
	echo "ERROR: --engine required" >&2
	usage
	;;
*)
	echo "ERROR: invalid --engine: $ENGINE (vep-rs|fastvep|perl)" >&2
	usage
	;;
esac
[[ -z "$DATA_DIR" ]] && {
	echo "ERROR: --data-dir required" >&2
	usage
}
[[ -z "$OUTPUT_DIR" ]] && {
	echo "ERROR: --output-dir required" >&2
	usage
}
[[ -d "$DATA_DIR" ]] || {
	echo "ERROR: data-dir not a directory: $DATA_DIR" >&2
	exit 1
}
case "$ENGINE" in
vep-rs)
	[[ -n "$VEP_BINARY" ]] || {
		echo "ERROR: --vep-binary required for engine vep-rs" >&2
		exit 1
	}
	[[ -x "$VEP_BINARY" ]] || {
		echo "ERROR: vep-rs binary not executable: $VEP_BINARY" >&2
		exit 1
	}
	;;
fastvep)
	[[ -n "$FASTVEP_BINARY" ]] || {
		echo "ERROR: --fastvep-binary required for engine fastvep" >&2
		exit 1
	}
	[[ -x "$FASTVEP_BINARY" ]] || {
		echo "ERROR: fastVEP binary not executable: $FASTVEP_BINARY" >&2
		exit 1
	}
	;;
perl)
	command -v docker >/dev/null 2>&1 || {
		echo "ERROR: docker required for engine perl" >&2
		exit 1
	}
	;;
esac

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] run_clone_measurement engine=$ENGINE suites=$SUITES_FILTER walltime_only=$WALLTIME_ONLY skip_fields=$SKIP_FIELDS fork=$FORK"
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] data_dir=$DATA_DIR output_dir=$OUTPUT_DIR repo=$REPO_DIR"

# Resource locations, all relative to DATA_DIR by default and all overridable by
# the matching environment variable. Override any one of them when your data
# tree is laid out differently; see the layout note in the file header.
#
# The ground-truth trees must hold Perl VEP output produced with the SAME cache
# release and reference FASTA the engines are given here, or the concordance
# figures compare two different annotations.
SNP_INDEL_GT_DIR="${VEP_SNP_INDEL_GT_DIR:-$DATA_DIR/ground_truth/perl/snp_indel}"
SV_GT_DIR="${VEP_SV_GT_DIR:-$DATA_DIR/ground_truth/perl/sv_per_vcf}"

# Per-assembly reference FASTA (shared across all engines), indexed (.fa + .fai).
FASTA_GRCH37="${VEP_FASTA_GRCH37:-$DATA_DIR/reference/grch37/genome.fa}"
FASTA_GRCH38="${VEP_FASTA_GRCH38:-$DATA_DIR/reference/grch38/genome.fa}"
fasta_for() { [[ "$1" == "GRCh37" ]] && echo "$FASTA_GRCH37" || echo "$FASTA_GRCH38"; }

# Engine-specific resources.
# vep-rs JSON caches. Build them with `vep-cache-builder`, or convert an existing
# Perl Storable cache with `vep-cache-converter`. The cache must carry populated
# 3' UTR sequences; a cache built without them scores measurably lower.
VEPRS_CACHE_GRCH37="${VEP_RS_CACHE_GRCH37:-$DATA_DIR/caches/vep-rs/grch37}"
VEPRS_CACHE_GRCH38="${VEP_RS_CACHE_GRCH38:-$DATA_DIR/caches/vep-rs/grch38}"
# fastVEP transcript cache (.bin) plus the GFF3 it was built from.
FASTVEP_CACHE_GRCH37="${VEP_FASTVEP_CACHE_GRCH37:-$DATA_DIR/caches/fastvep/grch37.bin}"
FASTVEP_CACHE_GRCH38="${VEP_FASTVEP_CACHE_GRCH38:-$DATA_DIR/caches/fastvep/grch38.bin}"
GFF3_GRCH37="${VEP_GFF3_GRCH37:-$DATA_DIR/reference/grch37/annotation.gff3}"
GFF3_GRCH38="${VEP_GFF3_GRCH38:-$DATA_DIR/reference/grch38/annotation.gff3}"
# Committed alongside this script. It corrects fastVEP's
# `Feature_type=Transcript` on intergenic rows, where Perl and vep-rs both emit
# `-`; feature type is part of the comparison key, so without it every
# intergenic row of a fastVEP run scores discordant.
FASTVEP_NORMALIZE="$SCRIPT_DIR/fastvep_normalize.py"
# Perl Storable cache (the --dir_cache mount; contains homo_sapiens/...). VEP resolves
# homo_sapiens/<cache_version>_<assembly> underneath it from --cache_version and --assembly.
#
# The flag is `--dir_cache`. `--full_cache_dir` is not a VEP flag, and VEP does not
# error on it: it ignores it, falls back to
# a default cache directory that holds nothing, finds no annotation source, and emits one
# warning per record instead of annotations, so the run produces zero data lines and exits 0.
#
# Nothing downstream would notice: the timed Perl run's output is written to a scratch
# directory and deleted, and concordance is scored against a pinned ground truth rather than
# this run's output, so a Perl wall time can measure a VEP that annotated nothing. The
# `assert_perl_annotated` gate below exists to catch that, and `assert_perl_cache_complete`
# fails any cell whose cache is a chromosome subset.
PERL_CACHE_DIR="${VEP_PERL_CACHE_DIR:-$DATA_DIR/caches/perl/vep-cache}"
# Perl VEP cache release. Must match the release the ground truth was built with.
PERL_CACHE_VERSION="${VEP_PERL_CACHE_VERSION:-115}"

veprs_cache_for() { [[ "$1" == "GRCh37" ]] && echo "$VEPRS_CACHE_GRCH37" || echo "$VEPRS_CACHE_GRCH38"; }
fastvep_cache_for() { [[ "$1" == "GRCh37" ]] && echo "$FASTVEP_CACHE_GRCH37" || echo "$FASTVEP_CACHE_GRCH38"; }
gff3_for() { [[ "$1" == "GRCh37" ]] && echo "$GFF3_GRCH37" || echo "$GFF3_GRCH38"; }

# Output dir + wall_times.csv. The schema is identical across all three engines,
# so runs can be concatenated into one uniform CSV.
mkdir -p "$OUTPUT_DIR"
WALL_CSV="$OUTPUT_DIR/wall_times.csv"
echo "engine,run_id,host_id,replicate,suite_id,suite_name,assembly,run_kind,wall_seconds,exit_code" >"$WALL_CSV"
# replicate is a provenance column: this run's 1..N index across a set of repeated
# runs, passed in by the caller. When unset (a single run) it defaults to the host
# id so the column is never empty.
[[ -n "$REPLICATE" ]] || REPLICATE="$HOST_ID"

# time_run <label> <out_time_txt> -- <cmd...>
# Runs the command under GNU /usr/bin/time -v when available, capturing the
# timing block to <out_time_txt>. Echoes "<wall_seconds>|<exit_code>". Never
# aborts the run on a non-zero engine exit (recorded as exit_code in the CSV
# instead).
#
# GNU time is not present on macOS or the BSDs, where /usr/bin/time takes no -v.
# Without it the wall time still comes from the surrounding `date` delta, so the
# measurement holds; only the peak-RSS block in <out_time_txt> is lost. Install
# GNU time (`brew install gnu-time`) and set VEP_TIME_BIN=gtime to recover it.
time_run() {
	local label="$1" out_txt="$2"
	shift 2
	[[ "$1" == "--" ]] && shift
	local start end rc
	start=$(date +%s.%N)
	set +e
	# GNU time propagates the command's exit status, so the fallback to a bare run
	# is decided by whether GNU time is usable, never by the command failing.
	if "${VEP_TIME_BIN:-/usr/bin/time}" -v true >/dev/null 2>&1; then
		"${VEP_TIME_BIN:-/usr/bin/time}" -v "$@" 2>"$out_txt"
	else
		"$@" 2>"$out_txt"
	fi
	rc=$?
	set -e
	end=$(date +%s.%N)
	# Prefer GNU time's "Elapsed (wall clock)"; fall back to the date delta.
	local wall
	wall=$(awk '
		/Elapsed \(wall clock\)/ {
			n=split($NF, a, ":");
			if (n==3) { print a[1]*3600 + a[2]*60 + a[3] }
			else if (n==2) { print a[1]*60 + a[2] }
			else { print $NF }
			found=1; exit
		}
		END { if (!found) print "" }
	' "$out_txt")
	if [[ -z "$wall" ]]; then
		wall=$(awk "BEGIN{printf \"%.3f\", $end - $start}")
	fi
	echo "$wall|$rc"
}

# csv_row <suite_id> <suite_name> <assembly> <run_kind> <wall> <rc>
csv_row() {
	echo "$ENGINE,$RUN_ID,$HOST_ID,$REPLICATE,$1,$2,$3,$4,$5,$6" >>"$WALL_CSV"
}

# Perl canonical-contig filter (Perl engine only). The cached SNP/indel ground
# truth was built
# from canonical-filtered inputs, so the timed Perl run filters the same way.
#
# ALSO drops per-sample genotypes (bcftools view -G -> sites-only). VEP annotates
# per (position, allele) and its tab output carries NO genotype columns, so this
# leaves the annotation output (and therefore concordance F1, which keys on
# location/allele/feature/consequence) byte-identical, and it is what makes the
# multi-sample suites' timed cells comparable across engines: on s06 (1KG_highcov,
# 3,202 samples) the genotype matrix is most of the input bytes, and its cost is
# measured separately by the `*m` suites. The other 9 suites are already sites-only,
# so -G is a harmless no-op for them.
#
# FAILS CLOSED at every step. There is no cp fallback: an unfiltered or
# genotype-carrying VCF written to the canonical path would be measured and
# reported as if the canonical filter had been applied, and on the multi-sample
# suites would time the genotype-carrying input under a sites-only label. A missing bcftools,
# a failed filter, or a failed index each abort the cell instead.
CANONICAL_CONTIGS_LIST="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT,M"
canonicalize() {
	local in_vcf="$1" out_vcf="$2"
	if ! command -v bcftools >/dev/null 2>&1; then
		# FAIL CLOSED. Every canonical-contigs guarantee in the reported results
		# rests on this filter, and the cached Perl ground truth was
		# generated from canonical-filtered inputs, so annotating the raw input
		# instead measures a DIFFERENT input than the one the comparison assumes.
		# A bare cp with only a WARN here would silently void the policy claim and
		# (on the multi-sample suites) time the genotype matrix, whose cost is
		# measured separately; genotypes change no consequence.
		echo "ERROR: [run_clone_measurement] bcftools unavailable; cannot produce a canonical-contig input for $in_vcf" >&2
		return 1
	fi
	[[ -f "${in_vcf}.tbi" ]] || tabix -p vcf "$in_vcf" 2>/dev/null || true
	local first_chrom src
	first_chrom=$(zcat -f "$in_vcf" | grep -v '^#' | head -1 | cut -f1 || echo "")
	src="$in_vcf"
	if [[ "$first_chrom" =~ ^chr ]]; then
		local stripped="${out_vcf%.vcf.gz}.stripped.vcf.gz"
		zcat -f "$in_vcf" |
			awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}' |
			bgzip >"$stripped"
		tabix -p vcf "$stripped" 2>/dev/null || true
		src="$stripped"
	fi
	# -G drops genotypes (sites-only); -t applies the canonical-contig filter.
	#
	# Both steps FAIL CLOSED. A bare `cp` of the raw input on bcftools failure
	# would leave an unfiltered, genotype-carrying VCF at the canonical path --
	# measured and reported as if the filter had been applied. A canonical
	# filter that can silently not apply is not a filter.
	if ! bcftools view -G -t "$CANONICAL_CONTIGS_LIST" -O z -o "$out_vcf" "$src" 2>/dev/null; then
		echo "ERROR: [run_clone_measurement] canonical-contig + sites-only filter failed for $src" >&2
		rm -f "$out_vcf"
		return 1
	fi
	# The index is required: the record-count assertions and Perl both read it.
	if ! tabix -p vcf "$out_vcf" 2>/dev/null; then
		echo "ERROR: [run_clone_measurement] tabix index failed for $out_vcf" >&2
		return 1
	fi
}

# Suite table: id:name:assembly:type:vcf_base
declare -a SUITES=(
	"s01:ClinVar_GRCh37:GRCh37:snp_indel:clinvar_grch37"
	"s02:ClinVar_GRCh38:GRCh38:snp_indel:clinvar_grch38"
	"s03:gnomAD_v2.1.1_chr21:GRCh37:snp_indel:gnomad_genomes_v2.1.1_chr21"
	"s04:gnomAD_v4.1_chr21:GRCh38:snp_indel:gnomad_genomes_v4.1_chr21"
	# s05 + s06 are the ONLY suites whose SOURCE inputs carry per-sample genotypes
	# (1000 Genomes; the high-coverage set has 3,202 samples, FORMAT=GT). Genotypes
	# do NOT affect VEP annotation output (VEP keys on position/allele; tab output has
	# no genotype columns), so the apples-to-apples cells (s05/s06) use the SITES-ONLY
	# inputs (bcftools view -G; the 4 other SNP/indel suites s01-s04 are sites-only at
	# source). The `*m` rows below measure the SAME suites WITH genotypes for the
	# sample-scaling measurement (vep-rs + fastVEP; this script's Perl branch always
	# strips genotypes, and Perl VEP's one genotype-carrying measurement is the probe
	# row of manuscript/data/perl_bench.csv).
	# Concordance F1 is identical between the sites-only and multi-sample forms.
	"s05:1KG_integrated_chr21:GRCh37:snp_indel:1kg_integrated_chr21_sitesonly"
	"s06:1KG_highcov_chr21:GRCh38:snp_indel:1kg_highcov_chr21"
	"s07:SV_GRCh37:GRCh37:sv:-"
	"s08:SV_GRCh38:GRCh38:sv:-"
	"s09:SV_real_world_GRCh37:GRCh37:sv:-"
	"s10:SV_real_world_GRCh38:GRCh38:sv:-"
	# Sample-scaling rows (with-genotypes; NOT part of --suites all / the headline table).
	# Run explicitly, e.g. --suites "s05m s06m". vcf_base names the genotype-carrying input.
	"s05m:1KG_integrated_chr21_multisample:GRCh37:snp_indel:1kg_integrated_chr21"
	"s06m:1KG_highcov_chr21_multisample:GRCh38:snp_indel:1kg_highcov_chr21_multisample"
)

# SV per-assembly dedupe: the canonical per-VCF SV GT carries all 16 VCFs per
# assembly in one directory (13 synthetic + 3 real-world). s07/s09 (GRCh37) and
# s08/s10 (GRCh38) measure the SAME per-assembly SV set, so the annotation +
# comparison is done ONCE per assembly under s07 (GRCh37) / s08 (GRCh38);
# s09/s10 record a pointer marker. compare_sv_concordance.py reports overall +
# per-file F1, so the synthetic vs real-world split is recoverable.
# bash 3.2 has no associative arrays, and macOS ships 3.2 while `bash -n` and
# shellcheck both PASS on `declare -A` -- it fails at RUNTIME, mid-run. A
# delimited string is exact here: at most two keys (the two assemblies) and the
# values are suite ids with no colons.
SV_DONE=""

# should_run <suite_id> <stype> : honor --suites filter.
should_run() {
	local sid="$1" stype="$2"
	# Sample-scaling rows (ids ending in `m`, e.g. s05m/s06m) are OPT-IN: they run ONLY
	# when named explicitly (--suites "s05m s06m"), never under all/snp-indel, so a
	# default run stays the sites-only apples-to-apples set and never times the
	# genotype matrix, whose cost is measured separately.
	case "$SUITES_FILTER" in
	all) [[ "$sid" != *m ]] ;;
	snp-indel) [[ "$stype" == "snp_indel" && "$sid" != *m ]] ;;
	sv) [[ "$stype" == "sv" ]] ;;
	*) [[ " $SUITES_FILTER " == *" $sid "* ]] ;;
	esac
}

# True when this suite has Perl ground truth to compare against. Distinguishes an
# EXPECTED absence from a missing file, because conflating the two makes a run
# report failure while producing perfect data.
#
# The sample-scaling suites (`s05m`/`s06m`) have no Perl ground truth: the ground truth is
# generated from the sites-only inputs, and genotypes change no consequence, so the
# sites-only truth is the truth for both forms (see `should_run` above). Symlinking a nonexistent
# target anyway leaves a DANGLING link, which a subsequent copy of the output tree is
# liable to treat as an error, turning a complete run into a reported failure.
#
# A missing ground truth on any OTHER suite is a data-tree error and fails the cell loudly.
# Keying on existence alone would let a suite whose ground truth is absent read as a
# clean skip, which is the same conflation in the opposite direction.
ground_truth_ready() {
	local sid="$1" sname="$2" asm="$3" gt="$4"
	if [[ "$sid" == *m ]]; then
		echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $sid is sample-scaling: wall time only, no Perl ground truth exists"
		return 1
	fi
	if [[ ! -f "$gt/output.txt" ]]; then
		echo "[ERROR] $sid: no Perl ground truth at $gt/output.txt" >&2
		csv_row "$sid" "$sname" "$asm" "ground_truth_missing" "" "1"
		return 1
	fi
	return 0
}

# Per-engine SNP/indel measurement
# measure_p1 covers the SNP/indel suites (suite type `snp_indel`).
measure_p1() {
	local suite_id="$1" suite_name="$2" assembly="$3" vcf_base="$4"
	local suite_dir="$OUTPUT_DIR/$suite_id"
	mkdir -p "$suite_dir/rust" "$suite_dir/report" "$suite_dir/perf"
	local fasta
	fasta=$(fasta_for "$assembly")
	# Canonical input (Ensembl-style contigs) -- the set the
	# ground_truth/perl/snp_indel/ Perl GT was generated from.
	local vcf="$DATA_DIR/vcf/$(printf %s "$assembly" | tr "[:upper:]" "[:lower:]")/all_variants_canonical/${vcf_base}.canonical.vcf.gz"
	local gt_dir="$SNP_INDEL_GT_DIR/$suite_id"

	case "$ENGINE" in
	vep-rs)
		local cache
		cache=$(veprs_cache_for "$assembly")
		local out="$suite_dir/rust/output.txt"
		local run_kind perf_txt res wall rc
		for run_kind in warmup timed; do
			perf_txt="$suite_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${suite_id}_${run_kind}" "$perf_txt" -- \
				"$VEP_BINARY" -i "$vcf" -o "$out" \
				--offline --json_cache "$cache" \
				--fasta "$fasta" \
				--fork "$FORK" --buffer_size "$BUFFER_SIZE" \
				--force --no_stats --quiet)
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $suite_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$suite_id" "$suite_name" "$assembly" "$run_kind" "$wall" "$rc"
		done
		if [[ "$WALLTIME_ONLY" != "true" ]] && ground_truth_ready "$suite_id" "$suite_name" "$assembly" "$gt_dir"; then
			mkdir -p "$suite_dir/perl"
			ln -sf "$gt_dir/output.txt" "$suite_dir/perl/output.txt"
			python3 "$REPO_DIR/scripts/concordance/compare_vep_outputs.py" \
				--perl-dir "$suite_dir/perl" \
				--rust-dir "$suite_dir/rust" \
				--report-dir "$suite_dir/report" \
				--assembly "$assembly" || {
				echo "ERROR: [run_clone_measurement] $suite_id SNP/indel compare failed; the cell has no F1 and this run is incomplete" >&2
				csv_row "$suite_id" "$suite_name" "$assembly" "compare_failed" "" "1"
				return 1
			}
			# Per-column and Extra-key agreement, beside the F1 report. Not F1-bearing:
			# a failure is logged in the CSV and the cell still completes.
			if [[ "$SKIP_FIELDS" != "true" ]]; then
				python3 "$REPO_DIR/scripts/concordance/compare_vep_fields.py" \
					--perl "$suite_dir/perl/output.txt" \
					--rust "$suite_dir/rust/output.txt" \
					--out "$suite_dir/report/fields.json" \
					--tmp-dir "$suite_dir" >/dev/null || {
					echo "WARNING: [run_clone_measurement] $suite_id field comparison failed; report/fields.json is absent" >&2
					csv_row "$suite_id" "$suite_name" "$assembly" "fields_compare_failed" "" "1"
				}
			else
				echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $suite_id field comparison skipped (--skip-fields)"
			fi
		fi
		;;
	fastvep)
		local cache gff3
		cache=$(fastvep_cache_for "$assembly")
		gff3=$(gff3_for "$assembly")
		local raw="$suite_dir/rust/raw.txt"
		local run_kind perf_txt res wall rc
		for run_kind in warmup timed; do
			perf_txt="$suite_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${suite_id}_${run_kind}" "$perf_txt" -- \
				env RAYON_NUM_THREADS="$FORK" "$FASTVEP_BINARY" annotate \
				--gff3 "$gff3" --fasta "$fasta" --transcript-cache "$cache" \
				-i "$vcf" -o "$raw" --output-format tab)
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $suite_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$suite_id" "$suite_name" "$assembly" "$run_kind" "$wall" "$rc"
		done
		if [[ "$WALLTIME_ONLY" != "true" ]] && ground_truth_ready "$suite_id" "$suite_name" "$assembly" "$gt_dir"; then
			# Normalize the timed raw output (intergenic Feature_type=Transcript -> -).
			# FAIL CLOSED. A `cp` fallback here would substitute un-normalised output
			# and depress fastVEP's F1 for a reason unrelated to annotation quality,
			# with no error surfaced. A missing or broken normaliser must abort the
			# cell instead.
			if [[ -s "$raw" ]]; then
				# `return 1`, never `exit 1`: this runs inside measure_p1, whose call
				# site is guarded (`measure_p1 ... || failed_suites+=(...)`). Exiting
				# here would abort the whole run, skipping every later suite and any
				# cleanup the caller performs.
				if [[ ! -f "$FASTVEP_NORMALIZE" ]]; then
					echo "[ERROR] fastvep normalizer not found at $FASTVEP_NORMALIZE" >&2
					csv_row "$suite_id" "$suite_name" "$assembly" "normalizer_failed" "" "1"
					return 1
				fi
				python3 "$FASTVEP_NORMALIZE" "$raw" "$suite_dir/rust/output.txt" || {
					echo "[ERROR] fastvep normalizer failed on $suite_id; refusing to compare un-normalised output" >&2
					csv_row "$suite_id" "$suite_name" "$assembly" "normalizer_failed" "" "1"
					return 1
				}
			fi
			mkdir -p "$suite_dir/perl"
			ln -sf "$gt_dir/output.txt" "$suite_dir/perl/output.txt"
			python3 "$REPO_DIR/scripts/concordance/compare_vep_outputs.py" \
				--perl-dir "$suite_dir/perl" \
				--rust-dir "$suite_dir/rust" \
				--report-dir "$suite_dir/report" \
				--assembly "$assembly" --glob 'output.txt' || {
				echo "ERROR: [run_clone_measurement] $suite_id SNP/indel compare failed; the cell has no F1 and this run is incomplete" >&2
				csv_row "$suite_id" "$suite_name" "$assembly" "compare_failed" "" "1"
				return 1
			}
		fi
		;;
	perl)
		# Perl is the ground-truth generator: it has no concordance branch. The
		# scratch annotation is timed then deleted (the cached GT is the source
		# of truth; this engine contributes wall-time only).
		#
		# Cache completeness is asserted BEFORE any work: a chromosome-subset cache
		# produces a fast, wrong wall time and exits 0.
		if ! assert_perl_cache_complete "$assembly"; then
			csv_row "$suite_id" "$suite_name" "$assembly" "perl_cache_incomplete" "" "1"
			return 1
		fi
		local canon_vcf="$suite_dir/inputs/${vcf_base}.canonical.vcf.gz"
		local raw_vcf="$DATA_DIR/vcf/$(printf %s "$assembly" | tr "[:upper:]" "[:lower:]")/all_variants/${vcf_base}.vcf.gz"
		mkdir -p "$suite_dir/inputs" "$suite_dir/perl_scratch"
		# canonicalize fails closed (see its definition). A cell whose input is not
		# canonical cannot be compared against the canonical-filtered ground truth,
		# so record the failure and skip rather than time a different input.
		if ! canonicalize "$raw_vcf" "$canon_vcf"; then
			echo "ERROR: [run_clone_measurement] $suite_id skipped: canonical input unavailable" >&2
			csv_row "$suite_id" "$suite_name" "$assembly" "canonicalize_failed" "" "1"
			return 1
		fi
		local run_kind perf_txt res wall rc
		for run_kind in warmup timed; do
			perf_txt="$suite_dir/perf/${run_kind}.time.txt"
			res=$(perl_annotate_timed "${suite_id}_${run_kind}" "$perf_txt" \
				"$canon_vcf" "$suite_dir/perl_scratch/output.txt" "$assembly")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $suite_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$suite_id" "$suite_name" "$assembly" "$run_kind" "$wall" "$rc"
		done
		# Before the scratch directory goes, confirm the timed run actually annotated. This is
		# the only point at which the evidence exists.
		if ! assert_perl_annotated "$suite_dir/perl_scratch/output.txt" "$canon_vcf" "$suite_id"; then
			csv_row "$suite_id" "$suite_name" "$assembly" "perl_not_annotating" "" "1"
			rm -rf "$suite_dir/perl_scratch" "$suite_dir/inputs"
			return 1
		fi
		rm -rf "$suite_dir/perl_scratch" "$suite_dir/inputs"
		;;
	esac
}

# perl_annotate_timed <label> <perf_txt> <in.vcf.gz> <out.txt> <assembly>
# Runs Perl VEP in Docker under /usr/bin/time. --fork matches the value the other
# engines are given, so the comparison is like for like; the SV ground truth's own
# fork setting belongs to generate_reference_output.sh, not to a timed cell.
perl_annotate_timed() {
	local label="$1" perf_txt="$2" in_vcf="$3" out_txt="$4" assembly="$5"
	local fasta out_dir out_base
	fasta=$(fasta_for "$assembly")
	out_dir=$(dirname "$out_txt")
	out_base=$(basename "$out_txt")
	mkdir -p "$out_dir"
	time_run "$label" "$perf_txt" -- \
		docker run --rm \
		--user "$(id -u):$(id -g)" \
		-v "${in_vcf}:/work/input.vcf.gz:ro" \
		-v "${out_dir}:/work/perl_out" \
		-v "${PERL_CACHE_DIR}:/work/cache:ro" \
		-v "${fasta}:/work/ref.fa:ro" \
		-v "${fasta}.fai:/work/ref.fa.fai:ro" \
		"$PERL_IMAGE" \
		vep \
		-i /work/input.vcf.gz \
		-o "/work/perl_out/${out_base}" \
		--offline --cache --dir_cache /work/cache --cache_version "$PERL_CACHE_VERSION" \
		--species homo_sapiens --assembly "$assembly" \
		--format vcf --buffer_size "$BUFFER_SIZE" --fork "$FORK" \
		--fasta /work/ref.fa \
		--no_headers --no_stats --quiet --force_overwrite
}

# assert_perl_annotated <perl_out_file_or_dir> <input_file_or_dir> <suite_id>
#
# A Perl invocation that cannot reach its cache does not fail: it emits one warning per record,
# zero annotation rows, and exit 0. Because the timed run's output is deleted immediately after
# the loop, nothing downstream ever sees it, so without this gate a Perl wall time can measure a
# VEP that annotated nothing.
#
# THE FLOOR IS REFERENCE-FREE ON PURPOSE. It needs no ground truth staged, so it cannot be
# skipped on a data tree that carries only inputs. A working VEP emits at least one row per input
# record (one per overlapping transcript, or a single intergenic row) and several per record on
# genome-wide and structural-variant inputs alike. So `rows >= records` fails a broken run, which
# sits at zero, and passes a legitimately sparse one with an order of magnitude spare.
#
# It is not specific to one flag: a wrong cache mount, a wrong --cache_version, a wrong
# --assembly and a missing cache all present identically as "almost no output".
#
# `n=$(...) || n=0`, never `$(... || echo 0)`: grep -c prints 0 AND exits 1 when nothing matches,
# so the latter form makes n the two-line string "0\n0", and every later numeric test on it dies
# with "integer expression expected" and falls through to its else branch -- reporting a pass
# because the comparison ERRORED rather than because it succeeded.
assert_perl_annotated() {
	local perl_out="$1" inputs="$2" suite_id="$3"
	local out_rows=0 in_records=0 f n

	if [[ -d "$perl_out" ]]; then
		for f in "$perl_out"/*.txt; do
			[[ -e "$f" ]] || continue
			# VEP writes <out>_warnings.txt beside its output; counting it would mask the
			# defect exactly, since the broken run puts one line per record in there.
			case "$f" in *_warnings.txt) continue ;; esac
			n=$(grep -vc '^#' "$f") || n=0
			out_rows=$((out_rows + n))
		done
	elif [[ -e "$perl_out" ]]; then
		out_rows=$(grep -vc '^#' "$perl_out") || out_rows=0
	fi

	for f in "$inputs" "$inputs"/*.vcf "$inputs"/*.vcf.gz; do
		[[ -f "$f" ]] || continue
		case "$f" in
		*.vcf.gz) n=$(zcat "$f" | grep -vc '^#') || n=0 ;;
		*.vcf) n=$(grep -vc '^#' "$f") || n=0 ;;
		*) continue ;;
		esac
		in_records=$((in_records + n))
	done

	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $suite_id perl annotation check: ${out_rows} rows / ${in_records} records"
	# ZERO INPUT RECORDS IS A FAILURE, NOT A PASS. The floor below is `out_rows >=
	# in_records`, so `in_records == 0` would make the whole gate vacuous: an absent or empty
	# input directory yields 0 rows against 0 records, the predicate is false, and the cell
	# is recorded as a clean, very fast measurement. `shopt -s nullglob` on the SV branch
	# makes an empty directory iterate zero times and exit 0, so nothing upstream objects
	# either, and the cache guard passes because the cache is fine. That is a wrong
	# number reached with both gates green.
	if [[ "$in_records" -eq 0 ]]; then
		echo "ERROR: [run_clone_measurement] $suite_id counted 0 input records" >&2
		echo "  The inputs are missing or empty, so the output-volume floor cannot bound anything." >&2
		echo "  A cell timed in this state records a near-zero wall time as a success." >&2
		return 1
	fi
	if [[ "$out_rows" -lt "$in_records" ]]; then
		echo "ERROR: [run_clone_measurement] $suite_id Perl emitted $out_rows rows for $in_records input records" >&2
		echo "  A Perl VEP run that cannot reach its cache emits warnings rather than annotations and still exits 0." >&2
		echo "  Check --dir_cache, the cache mount, --cache_version and --assembly before trusting this cell." >&2
		return 1
	fi
	return 0
}

# assert_perl_cache_complete <assembly>
#
# A chromosome-subset cache (a `115_GRCh37/` holding only `21/` against a genome-wide
# `115_GRCh38/`) is a second way for Perl to annotate almost nothing and exit 0, independent
# of the flag.
#
# `assert_perl_annotated` does NOT subsume this, which is the whole reason the function exists.
# Its reference-free `rows >= records` floor catches a chr21-only cache only where the input is
# genome-wide. Over the 16 GRCh37 SV VCFs, whose records are chr21 by design, the same broken
# cache passes the floor with an order of magnitude spare while `09_breakends` silently loses
# every mate that points off chromosome 21. An output-volume gate cannot see a cache gap the
# input barely samples, so completeness is asserted on the cache itself, before any timing.
#
# Primary contigs only. Ensembl's full GRCh37 cache carries 699 directories, the other 674 being
# unplaced scaffolds and haplotype contigs, and no input here reaches them -- so requiring all
# 699 would fail a legitimately pruned cache. Requiring the 25 primary contigs fails every
# chromosome-subset cache, which is the defect class.
assert_perl_cache_complete() {
	local assembly="$1"
	# `:-115` is load-bearing rather than defensive. The CACHE-COMPLETE gate in
	# test_run_clone_measurement.sh extracts this function out of the harness with awk and evals
	# it, deliberately, so the gate exercises the shipped code instead of a reimplementation
	# that can drift from it. In that context the harness's own
	# `PERL_CACHE_VERSION` assignment has not run, and the test suite is `set -u`, so a bare
	# reference here aborts the whole gate suite on an unbound variable -- exiting 1 with no
	# failing assertion printed, which reads as a harness defect rather than a missing default.
	local ver="${PERL_CACHE_VERSION:-115}"
	local slice="$PERL_CACHE_DIR/homo_sapiens/${ver}_${assembly}"
	local primaries="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
	local c missing=""

	if [[ ! -d "$slice" ]]; then
		echo "ERROR: [run_clone_measurement] Perl cache slice not found: $slice" >&2
		echo "  Expected homo_sapiens/${ver}_${assembly} under \$PERL_CACHE_DIR ($PERL_CACHE_DIR)." >&2
		return 1
	fi

	for c in $primaries; do
		[[ -d "$slice/$c" ]] || missing="$missing $c"
	done

	if [[ -n "$missing" ]]; then
		echo "ERROR: [run_clone_measurement] Perl cache for $assembly is missing primary contigs:$missing" >&2
		echo "  Slice: $slice" >&2
		echo "  A chromosome-subset cache does not fail a VEP run. It under-annotates every record" >&2
		echo "  outside the contigs it holds, exits 0, and yields a wall time that is too fast." >&2
		echo "  Do NOT time this cell against it." >&2
		return 1
	fi

	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   perl cache check ($assembly): all 25 primary contigs present"
	return 0
}

# Per-engine SV measurement (per-VCF over the canonical 16-VCF set per assembly)
# assert_sv_inputs_present <inputs_dir> <suite_id>
# assert_sv_outputs_present <outputs_dir> <n_inputs> <suite_id>
#
# The per-VCF loop runs under `nullglob`, so an input directory that is absent, empty or
# holds only manifest.json is not an error to bash: the loop body never runs, `time_run`
# reports a 0 s wall time with rc=0, the cell writes no output, and the comparison that
# follows finds nothing to score. A cell that annotated nothing must FAIL, so both ends
# of the loop are asserted: at least one input VCF before it runs, and one non-empty
# output per input after it. Reference-free on purpose: neither needs the ground truth
# staged.
sv_input_count() {
	local dir="$1" n=0 f b
	shopt -s nullglob
	for f in "$dir"/*.vcf "$dir"/*.vcf.gz; do
		b=$(basename "$f"); b=${b%.vcf.gz}; b=${b%.vcf}
		[[ "$b" == "manifest" ]] && continue
		n=$((n + 1))
	done
	shopt -u nullglob
	echo "$n"
}
assert_sv_inputs_present() {
	local dir="$1" sid="$2" n
	n=$(sv_input_count "$dir")
	if [[ "$n" -lt 1 ]]; then
		echo "ERROR: [run_clone_measurement] $sid: no SV input VCF under $dir (SV_GT_DIR=$SV_GT_DIR); the cell would time an empty loop" >&2
		return 1
	fi
	return 0
}
assert_sv_outputs_present() {
	local dir="$1" want="$2" sid="$3" n=0 f
	shopt -s nullglob
	for f in "$dir"/*.txt; do
		if grep -q -v '^#' "$f" 2>/dev/null; then n=$((n + 1)); fi
	done
	shopt -u nullglob
	if [[ "$n" -lt "$want" ]]; then
		echo "ERROR: [run_clone_measurement] $sid: $n non-empty SV output file(s) under $dir for $want input VCF(s); the engine annotated nothing on at least one file" >&2
		return 1
	fi
	return 0
}

measure_sv() {
	local suite_id="$1" suite_name="$2" assembly="$3"
	local suite_dir="$OUTPUT_DIR/$suite_id"
	mkdir -p "$suite_dir/rust" "$suite_dir/report" "$suite_dir/perf"
	# `${var,,}` is bash 4+ and fails at runtime on macOS bash 3.2 while passing
	# `bash -n` and shellcheck.
	local asm_lc
	asm_lc=$(printf %s "$assembly" | tr "[:upper:]" "[:lower:]")
	local sv_inputs="$SV_GT_DIR/canonical_inputs/$asm_lc"
	local sv_gt="$SV_GT_DIR/$asm_lc"
	if ! assert_sv_inputs_present "$sv_inputs" "$suite_id"; then
		csv_row "$suite_id" "$suite_name" "$assembly" "sv_inputs_missing" "" "1"
		return 1
	fi
	local sv_n_inputs
	sv_n_inputs=$(sv_input_count "$sv_inputs")

	# Per-assembly dedupe: annotate once under s07/s08; s09/s10 are pointers.
	if [[ "$SV_DONE" == *":${asm_lc}="* ]]; then
		sv_done_suite="${SV_DONE#*:${asm_lc}=}"
		sv_done_suite="${sv_done_suite%%:*}"
		echo "SV $assembly annotated under ${sv_done_suite}; see that suite's report (per-file F1 covers synthetic + real-world)" \
			>"$suite_dir/POINTER"
		csv_row "$suite_id" "$suite_name" "$assembly" "pointer" "0" "0"
		return 0
	fi
	SV_DONE="${SV_DONE}:${asm_lc}=${suite_id}:"

	local fasta
	fasta=$(fasta_for "$assembly")

	case "$ENGINE" in
	vep-rs)
		local cache
		cache=$(veprs_cache_for "$assembly")
		local run_kind rust_out perf_txt res wall rc
		for run_kind in warmup timed; do
			rust_out="$suite_dir/rust"
			perf_txt="$suite_dir/perf/${run_kind}.time.txt"
			# Wrap the whole per-VCF loop in one /usr/bin/time so wall_seconds is
			# the full SV-cell annotation time.
			#
			# `$FORK` here, NOT a hardcoded 1. Every engine takes the same --fork on
			# every cell: the Perl branch below and the fastVEP branch above both pass
			# `$FORK` on the SV cells too, so hardcoding 1 here would make vep-rs the
			# only engine of the three annotating serially and turn every SV ratio
			# into a floor rather than a measurement. The 1 is not needed "to match
			# the per-VCF SV GT": the GT is a fixed archived artifact, and vep-rs
			# output is fork-invariant by construction -- `annotate_batch` dispatches
			# the same `annotate_one` over the same `&mut` batch through either
			# `par_iter_mut` or `iter_mut` and emits in input order (`runner.rs`), and
			# the per-chromosome cache parse runs on the global rayon pool regardless
			# (`json_cache.rs::load_transcripts_for_chr`). What --fork gates is only
			# the per-variant consequence fan-out, which is why a serial run is slowest
			# on the one 32,630-variant GRCh38 file.
			res=$(time_run "${suite_id}_${run_kind}" "$perf_txt" -- \
				bash -c '
					set -uo pipefail
					bin="$1"; cache="$2"; fasta="$3"; asm="$4"; buf="$5"; in_dir="$6"; out_dir="$7"; fork="$8"
					shopt -s nullglob
					failed=0
					for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
						b=$(basename "$vcf"); b=${b%.vcf.gz}; b=${b%.vcf}
						[[ "$b" == "manifest" ]] && continue
						"$bin" -i "$vcf" -o "$out_dir/$b.txt" \
							--offline --json_cache "$cache" --fasta "$fasta" \
							--fork "$fork" --buffer_size "$buf" --force --no_stats \
							--species homo_sapiens --assembly "$asm" --quiet \
							2>"$out_dir/$b.stderr" || {
							failed=$((failed+1))
							echo "ERROR: [run_clone_measurement] per-VCF invocation failed: $b" >&2
						}
					done
					[[ "$failed" -eq 0 ]] || {
						echo "ERROR: [run_clone_measurement] $failed per-VCF invocation(s) failed" >&2
						exit "$failed"
					}
				' _ "$VEP_BINARY" "$cache" "$fasta" "$assembly" "$BUFFER_SIZE" "$sv_inputs" "$rust_out" "$FORK")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $suite_id $run_kind wall=${wall}s rc=$rc"
			# One non-empty output per input, or the cell is a failure however fast it
			# ran: a 0 s SV cell with rc=0 is the signature of an empty input glob.
			if ! assert_sv_outputs_present "$rust_out" "$sv_n_inputs" "$suite_id"; then
				csv_row "$suite_id" "$suite_name" "$assembly" "${run_kind}_no_output" "$wall" "1"
				return 1
			fi
			csv_row "$suite_id" "$suite_name" "$assembly" "$run_kind" "$wall" "$rc"
		done
		if [[ "$WALLTIME_ONLY" != "true" ]]; then
			# `--vep-rs-cache` is what activates the cross-chromosome mask, which removes
			# Perl tuples naming a transcript the variant's chromosome does not hold. The
			# cache is the same one vep-rs just annotated against, so it is the right
			# authority. Passing it is NOT optional for a reported run: without it the
			# mask is inactive and the report says so.
			python3 "$REPO_DIR/scripts/validation/compare_sv_concordance.py" \
				--assembly "$asm_lc" \
				--input-dir "$sv_inputs" \
				--perl-dir "$sv_gt" \
				--rust-dir "$suite_dir/rust" \
				--vep-rs-cache "$cache" \
				--output-dir "$suite_dir/report" || {
				echo "ERROR: [run_clone_measurement] $suite_id SV compare failed; the cell has no F1 and this run is incomplete" >&2
				csv_row "$suite_id" "$suite_name" "$assembly" "compare_failed" "" "1"
				return 1
			}
		fi
		;;
	fastvep)
		local cache gff3
		cache=$(fastvep_cache_for "$assembly")
		gff3=$(gff3_for "$assembly")
		local run_kind rust_out perf_txt res wall rc
		for run_kind in warmup timed; do
			rust_out="$suite_dir/rust"
			perf_txt="$suite_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${suite_id}_${run_kind}" "$perf_txt" -- \
				bash -c '
					set -uo pipefail
					bin="$1"; cache="$2"; fasta="$3"; gff3="$4"; in_dir="$5"; out_dir="$6"; threads="$7"
					shopt -s nullglob
					failed=0
					for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
						b=$(basename "$vcf"); b=${b%.vcf.gz}; b=${b%.vcf}
						[[ "$b" == "manifest" ]] && continue
						RAYON_NUM_THREADS="$threads" "$bin" annotate \
							--gff3 "$gff3" --fasta "$fasta" --transcript-cache "$cache" \
							-i "$vcf" -o "$out_dir/$b.txt" --output-format tab \
							2>"$out_dir/$b.stderr" || {
							failed=$((failed+1))
							echo "ERROR: [run_clone_measurement] per-VCF invocation failed: $b" >&2
						}
					done
					[[ "$failed" -eq 0 ]] || {
						echo "ERROR: [run_clone_measurement] $failed per-VCF invocation(s) failed" >&2
						exit "$failed"
					}
				' _ "$FASTVEP_BINARY" "$cache" "$fasta" "$gff3" "$sv_inputs" "$rust_out" "$FORK")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $suite_id $run_kind wall=${wall}s rc=$rc"
			# Same gate as the vep-rs arm, at a floor of ONE non-empty output rather than
			# one per input: fastVEP's SV CONCORDANCE is best-effort below (it may decline
			# individual synthetic files), but a timing cell whose whole loop produced
			# nothing is the empty-glob signature, not a measurement.
			if ! assert_sv_outputs_present "$rust_out" 1 "$suite_id"; then
				csv_row "$suite_id" "$suite_name" "$assembly" "${run_kind}_no_output" "$wall" "1"
				return 1
			fi
			csv_row "$suite_id" "$suite_name" "$assembly" "$run_kind" "$wall" "$rc"
		done
		# SV concordance is best-effort for fastVEP: a loop that produces no non-empty
		# output is skipped, never hard-failed.
		if [[ "$WALLTIME_ONLY" != "true" ]]; then
			local produced
			produced=$(find "$suite_dir/rust" -name '*.txt' -size +0c 2>/dev/null | wc -l | tr -d ' ')
			if [[ "$produced" -eq 0 ]]; then
				echo "fastVEP produced no non-empty SV output; SV concordance skipped" \
					>"$suite_dir/report/SV_CONCORDANCE_SKIPPED.txt"
				echo "[WARN] $suite_id ($assembly) fastVEP SV concordance skipped: no usable output"
			else
				# `--vep-rs-cache` on the fastVEP arm too, so this arm's adjusted F1 applies
				# the SAME exclusion set as the vep-rs arm's. Without it the cross-chromosome
				# class is inactive here and the two engines' adjusted SV columns sit on
				# different masks.
				# `--scored-engine fastvep` is required with it: the filter's sanity assertion
				# that the scored engine cannot name a transcript the cache lacks for the
				# variant's chromosome holds for vep-rs and not for an engine carrying its own
				# annotation release, and it aborts the whole comparison when it fires.
				#
				# `veprs_cache_for`, NOT the local `$cache`, which in this branch is fastVEP's
				# OWN cache. The flag wants the vep-rs JSON cache: it is read as the authority
				# on which chromosome a transcript sits on, in the layout that loader parses,
				# and it is not the transcript set the scored engine annotated against.
				# The mask authority is the vep-rs JSON cache, which fastVEP itself never
				# reads, so a fastVEP data tree can be complete for annotation and still
				# lack it; assert it by name before the comparator does with a traceback.
				if [[ ! -d "$(veprs_cache_for "$assembly")" ]]; then
					echo "ERROR: [run_clone_measurement] $suite_id: vep-rs JSON cache absent at $(veprs_cache_for "$assembly"); the fastVEP tree must carry it for the cross-chromosome mask (--vep-rs-cache)" >&2
					csv_row "$suite_id" "$suite_name" "$assembly" "mask_cache_missing" "" "1"
					return 1
				fi
				python3 "$REPO_DIR/scripts/validation/compare_sv_concordance.py" \
					--assembly "$asm_lc" \
					--input-dir "$sv_inputs" \
					--perl-dir "$sv_gt" \
					--rust-dir "$suite_dir/rust" \
					--vep-rs-cache "$(veprs_cache_for "$assembly")" \
					--scored-engine fastvep \
					--output-dir "$suite_dir/report" || {
					# A compare that RAN and failed is not the documented no-output case
					# handled above: it means a missing input (the vep-rs cache the mask
					# reads, a ground-truth file) or a comparator error, and the cell has
					# no F1.
					echo "ERROR: [run_clone_measurement] $suite_id fastVEP SV compare failed; the cell has no F1 and this run is incomplete" >&2
					csv_row "$suite_id" "$suite_name" "$assembly" "compare_failed" "" "1"
					return 1
				}
			fi
		fi
		;;
	perl)
		# Perl SV: time over the canonical per-VCF SV inputs (already canonical at
		# GT-gen time). Scratch output deleted (Perl is GT; wall-time only).
		#
		# This is the branch `assert_perl_annotated` cannot protect: these inputs are
		# chr21 by design, so a chr21-only cache clears its rows/records floor while
		# `09_breakends` loses every off-chromosome mate. See the guard's own comment.
		if ! assert_perl_cache_complete "$assembly"; then
			csv_row "$suite_id" "$suite_name" "$assembly" "perl_cache_incomplete" "" "1"
			return 1
		fi
		local run_kind perf_txt res wall rc
		mkdir -p "$suite_dir/perl_scratch"
		for run_kind in warmup timed; do
			perf_txt="$suite_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${suite_id}_${run_kind}" "$perf_txt" -- \
				bash -c '
					set -uo pipefail
					image="$1"; cache="$2"; fasta="$3"; asm="$4"; in_dir="$5"; scratch="$6"; uidgid="$7"; buf="$8"; fork="$9"
					mkdir -p "$scratch"
					shopt -s nullglob
					failed=0
					for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
						b=$(basename "$vcf"); b=${b%.vcf.gz}; b=${b%.vcf}
						[[ "$b" == "manifest" ]] && continue
						if [[ "$vcf" == *.vcf ]]; then
							bgzip -c "$vcf" >"$scratch/$b.vcf.gz" 2>/dev/null || cp "$vcf" "$scratch/$b.vcf.gz"
							invcf="$scratch/$b.vcf.gz"
						else
							invcf="$vcf"
						fi
						docker run --rm --user "$uidgid" \
							-v "$invcf:/work/input.vcf.gz:ro" \
							-v "$scratch:/work/perl_out" \
							-v "$cache:/work/cache:ro" \
							-v "$fasta:/work/ref.fa:ro" \
							-v "$fasta.fai:/work/ref.fa.fai:ro" \
							"$image" \
							vep -i /work/input.vcf.gz -o "/work/perl_out/$b.txt" \
							--offline --cache --dir_cache /work/cache --cache_version "$PERL_CACHE_VERSION" \
							--species homo_sapiens --assembly "$asm" \
							--format vcf --buffer_size "$buf" --fork "$fork" \
							--fasta /work/ref.fa \
							--no_headers --no_stats --quiet --force_overwrite || {
							failed=$((failed+1))
							echo "ERROR: [run_clone_measurement] per-VCF invocation failed: $b" >&2
						}
					done
					[[ "$failed" -eq 0 ]] || {
						echo "ERROR: [run_clone_measurement] $failed per-VCF invocation(s) failed" >&2
						exit "$failed"
					}
				' _ "$PERL_IMAGE" "$PERL_CACHE_DIR" "$fasta" "$assembly" \
				"$sv_inputs" "$suite_dir/perl_scratch" "$(id -u):$(id -g)" "$BUFFER_SIZE" "$FORK")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $suite_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$suite_id" "$suite_name" "$assembly" "$run_kind" "$wall" "$rc"
		done
		# Same gate as the SNP/indel path: the per-VCF Perl output is about to be discarded.
		if ! assert_perl_annotated "$suite_dir/perl_scratch" "$sv_inputs" "$suite_id"; then
			csv_row "$suite_id" "$suite_name" "$assembly" "perl_not_annotating" "" "1"
			rm -rf "$suite_dir/perl_scratch"
			return 1
		fi
		rm -rf "$suite_dir/perl_scratch"
		;;
	esac
}

# Main suite loop
# RECORD, DON'T ABORT. A failing cell is recorded (the measure_* helpers already
# write a canonicalize_failed or per-VCF-failure row via csv_row) and the loop
# moves on; a caller rejects the run by reading the CSV rows whose exit_code is
# non-zero.
#
# The `|| failed_suites+=(...)` guard is load-bearing under `set -euo pipefail`.
# Both functions fail closed and return non-zero, so calling them bare would let
# the FIRST bad cell kill the script -- contradicting the "record the failure and
# skip" contract they document. A caller that also runs under `set -e` with no
# trap inherits that abort, so a single fail-closed cell would skip every
# remaining step of the surrounding run and produce no results at all.
failed_suites=()
for entry in "${SUITES[@]}"; do
	IFS=':' read -r suite_id suite_name assembly stype vcf_base <<<"$entry"
	if ! should_run "$suite_id" "$stype"; then
		echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] suite $suite_id skipped (suites=$SUITES_FILTER)"
		continue
	fi
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] suite $suite_id ($suite_name $assembly $stype) starting"
	if [[ "$stype" == "snp_indel" ]]; then
		measure_p1 "$suite_id" "$suite_name" "$assembly" "$vcf_base" ||
			failed_suites+=("$suite_id")
	else
		measure_sv "$suite_id" "$suite_name" "$assembly" ||
			failed_suites+=("$suite_id")
	fi
done

# container_probe <n>: n bare starts of the Perl image, the same --user flag as the
# timed cells and no mounts, so the wall time is the container start-up alone. Every
# probe is timed: the image is already pulled and warm from the cells above, which is
# the state every timed Perl invocation starts from. The rows carry the same
# provenance columns as wall_times.csv, keyed by probe index.
container_probe() {
	local n="$1" csv="$OUTPUT_DIR/container_start.csv" perf_dir="$OUTPUT_DIR/container_probe"
	mkdir -p "$perf_dir"
	echo "engine,run_id,host_id,replicate,image,probe_index,wall_seconds,exit_code" >"$csv"
	local i res wall rc
	for ((i = 1; i <= n; i++)); do
		res=$(time_run "container_probe_$i" "$perf_dir/probe-$i.time.txt" -- 			docker run --rm --user "$(id -u):$(id -g)" "$PERL_IMAGE" true)
		wall=${res%|*}
		rc=${res#*|}
		echo "$ENGINE,$RUN_ID,$HOST_ID,$REPLICATE,$PERL_IMAGE,$i,$wall,$rc" >>"$csv"
		if [[ "$rc" -ne 0 ]]; then
			echo "ERROR: [run_clone_measurement] container probe $i exited $rc; a probe with a non-zero exit must not enter an aggregate" >&2
		fi
	done
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] container probe: $n starts of $PERL_IMAGE written to $csv"
}

if [[ "$CONTAINER_PROBE" -gt 0 ]]; then
	container_probe "$CONTAINER_PROBE"
fi

if [[ ${#failed_suites[@]} -gt 0 ]]; then
	echo "ERROR: [run_clone_measurement] ${#failed_suites[@]} suite(s) failed: ${failed_suites[*]}" >&2
	echo "  Their rows in $WALL_CSV carry a non-zero rc; treat this run as failed." >&2
fi

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] run_clone_measurement complete (engine=$ENGINE); wall_times.csv at $WALL_CSV"
