#!/usr/bin/env bash
#
# generate_reference_output.sh -- write the Perl VEP reference outputs that
# run_clone_measurement.sh scores every engine against.
#
# The measurement harness READS an already-populated data tree (its header documents
# the layout) and nothing else in this repository writes the Perl ground truth into
# it. This script is that writer. It runs Ensembl VEP in the `ensemblorg/ensembl-vep`
# Docker image with the exact invocation the harness embeds for its Perl engine, and
# writes the outputs where the harness and both comparators expect them:
#
#   SNP/indel (s01-s06), one file per suite:
#     <data>/ground_truth/perl/snp_indel/<suite>/output.txt
#     from <data>/vcf/<asm>/all_variants_canonical/<base>.canonical.vcf.gz
#   Structural variants (s07/s08; s09/s10 are the harness's pointer suites), one
#   file per input VCF, 16 per assembly:
#     <data>/ground_truth/perl/sv_per_vcf/<asm>/<vcf-basename>.txt
#     from <data>/ground_truth/perl/sv_per_vcf/canonical_inputs/<asm>/*.vcf[.gz]
#
# Beside every output it writes <name>.provenance.json (image digest, flags, input and
# output digests, record and row counts, wall time) and renames VEP's own
# `<name>.txt_warnings.txt` companion to <name>.warnings.log, so no `*.txt` glob can
# ever pair a warnings file as an engine output.
#
# FORK AND BUFFER ARE PART OF THE GROUND TRUTH'S IDENTITY, NOT A TUNING KNOB.
#   SNP/indel: --fork 16 --buffer_size 5000, the settings the reported SNP/indel
#              reference set was generated with.
#   SV:        --fork 4 --buffer_size 5000, the settings the reported per-VCF SV
#              reference set was generated with (the `invocation` field of the SV
#              reference's provenance record under docs/concordance-provenance/).
#              VEP's transcript selection above --max_sv_size is batch-dependent, so
#              the fork count changes which variants share an input buffer and
#              therefore the row count. A regeneration must reuse the ground truth's
#              OWN settings, never the timing harness's.
#
# THE FLAG IS `--dir_cache`. `--full_cache_dir` is not a VEP flag: VEP ignores it,
# finds no annotation source, emits one warning per record and exits 0. The
# output-volume gate below (rows >= input records per suite or per SV set, the same
# floor the harness's `assert_perl_annotated` applies) makes that failure mode loud
# here as well.
#
# THE CACHE MUST BE THE FULL-GENOME ENSEMBL CACHE. A chromosome-subset cache does not
# fail a VEP run; it under-annotates every record whose transcripts lie outside the
# contigs it holds and exits 0. Completeness is asserted on the cache tree
# itself (all 25 primary contigs under homo_sapiens/<version>_<assembly>/), before
# anything runs, the same check the harness's `assert_perl_cache_complete` makes.
#
# Obtaining the cache (about 24 GB per assembly, gzip tar), with <v> the Ensembl
# release, 115:
#   GRCh38: https://ftp.ensembl.org/pub/release-<v>/variation/indexed_vep_cache/homo_sapiens_vep_<v>_GRCh38.tar.gz
#   GRCh37: https://ftp.ensembl.org/pub/grch37/release-<v>/variation/indexed_vep_cache/homo_sapiens_vep_<v>_GRCh37.tar.gz
#   Extract both under ONE directory, so that it contains
#   homo_sapiens/115_GRCh38/ and homo_sapiens/115_GRCh37/; that directory is
#   --perl-cache-dir (default <data>/caches/perl/vep-cache). VEP resolves
#   homo_sapiens/<cache_version>_<assembly> beneath the --dir_cache mount from
#   --cache_version and --assembly. Readiness is an artifact test, not a log line:
#   count the contig directories and confirm no tar process is still writing.
#
# The reference FASTA must be the primary assembly the cache was built against
# (Ensembl's Homo_sapiens.GRCh37.75.dna.primary_assembly.fa and
# Homo_sapiens.GRCh38.dna.primary_assembly.fa), bgzipped or plain, with a .fai
# beside it; the default location is <data>/reference/<asm>/genome.fa. VEP receives --fasta on every call: the comparators score against
# outputs generated WITH --fasta, and an asymmetric run hides UTR-level divergences.
#
# This script never fetches anything, needs no cloud credentials, and refuses to
# overwrite an existing reference output unless --force is passed: a reference is a
# pinned artifact whose digests the provenance records carry.
#
# Usage:
#   scripts/concordance/generate_reference_output.sh --data-dir <data> \
#     [--assembly GRCh37|GRCh38|both] [--suites all|snp-indel|sv|"s01 s07"] \
#     [--perl-image ensemblorg/ensembl-vep:release_115.2] \
#     [--perl-cache-dir <dir>] [--cache-version 115] \
#     [--fasta <fa>] [--fasta-grch37 <fa>] [--fasta-grch38 <fa>] \
#     [--fork 16] [--sv-fork 4] [--buffer-size 5000] \
#     [--snp-indel-gt-dir <dir>] [--sv-gt-dir <dir>] [--sv-inputs-dir <dir>] \
#     [--force] [--dry-run]
#
#   --dry-run runs every preflight check (inputs, FASTA, cache completeness, output
#   collisions) and prints the docker command for each output instead of running it.
#   A missing input fails a dry run exactly as it fails a real one.
set -euo pipefail

# Defaults. The data-tree locations and the VEP_* environment overrides are the same
# ones run_clone_measurement.sh reads, so a tree generated here is a tree the harness
# measures without further configuration.
DATA_DIR=""
ASSEMBLY="both"
SUITES_FILTER="all"
PERL_IMAGE="ensemblorg/ensembl-vep:release_115.2"
PERL_CACHE_DIR="${VEP_PERL_CACHE_DIR:-}"
PERL_CACHE_VERSION="${VEP_PERL_CACHE_VERSION:-115}"
FASTA_ANY=""
FASTA_GRCH37="${VEP_FASTA_GRCH37:-}"
FASTA_GRCH38="${VEP_FASTA_GRCH38:-}"
FORK=16
SV_FORK=4
BUFFER_SIZE=5000
SNP_INDEL_GT_DIR="${VEP_SNP_INDEL_GT_DIR:-}"
SV_GT_DIR="${VEP_SV_GT_DIR:-}"
SV_INPUTS_DIR=""
FORCE="false"
DRY_RUN="false"

usage() {
	cat <<'EOF'
Generate the Perl VEP reference outputs run_clone_measurement.sh scores against.

Required:
  --data-dir <path>        Data tree (vcf/, caches/, reference/, ground_truth/)

Optional:
  --assembly <A>           GRCh37 | GRCh38 | both (default both)
  --suites <filter>        all (default) | snp-indel | sv | space-separated "s01 s07"
                           (s09/s10 are aliases of s07/s08: the real-world SV VCFs are
                           part of the 16-file per-assembly set)
  --perl-image <ref>       Docker image (default ensemblorg/ensembl-vep:release_115.2)
  --perl-cache-dir <dir>   Directory holding homo_sapiens/<ver>_<asm>/
                           (default <data>/caches/perl/vep-cache)
  --cache-version <N>      Ensembl cache release (default 115)
  --fasta <fa>             Reference FASTA; valid only with a single --assembly
  --fasta-grch37 <fa>      Per-assembly FASTA overrides (defaults under
  --fasta-grch38 <fa>      <data>/reference/<asm>/)
  --fork <N>               SNP/indel --fork (default 16)
  --sv-fork <N>            SV --fork (default 4; the shipped SV ground truth's setting)
  --buffer-size <N>        --buffer_size for both (default 5000)
  --snp-indel-gt-dir <dir> Where <suite>/output.txt is written
  --sv-gt-dir <dir>        Where <asm>/<vcf-basename>.txt is written
  --sv-inputs-dir <dir>    Parent of the per-assembly SV input directories
                           (default <sv-gt-dir>/canonical_inputs)
  --force                  Overwrite an existing reference output
  --dry-run                Preflight everything, print the docker commands, run nothing
EOF
	exit 1
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--data-dir)
		DATA_DIR="$2"
		shift 2
		;;
	--assembly)
		ASSEMBLY="$2"
		shift 2
		;;
	--suites)
		SUITES_FILTER="$2"
		shift 2
		;;
	--perl-image)
		PERL_IMAGE="$2"
		shift 2
		;;
	--perl-cache-dir)
		PERL_CACHE_DIR="$2"
		shift 2
		;;
	--cache-version)
		PERL_CACHE_VERSION="$2"
		shift 2
		;;
	--fasta)
		FASTA_ANY="$2"
		shift 2
		;;
	--fasta-grch37)
		FASTA_GRCH37="$2"
		shift 2
		;;
	--fasta-grch38)
		FASTA_GRCH38="$2"
		shift 2
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--sv-fork)
		SV_FORK="$2"
		shift 2
		;;
	--buffer-size)
		BUFFER_SIZE="$2"
		shift 2
		;;
	--snp-indel-gt-dir)
		SNP_INDEL_GT_DIR="$2"
		shift 2
		;;
	--sv-gt-dir)
		SV_GT_DIR="$2"
		shift 2
		;;
	--sv-inputs-dir)
		SV_INPUTS_DIR="$2"
		shift 2
		;;
	--force)
		FORCE="true"
		shift
		;;
	--dry-run)
		DRY_RUN="true"
		shift
		;;
	-h | --help) usage ;;
	*)
		echo "ERROR: [generate_reference_output] unknown argument: $1" >&2
		usage
		;;
	esac
done

die() {
	echo "ERROR: [generate_reference_output] $*" >&2
	exit 1
}
log() {
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
}

# Validate arguments
[[ -n "$DATA_DIR" ]] || {
	echo "ERROR: [generate_reference_output] --data-dir required" >&2
	usage
}
[[ -d "$DATA_DIR" ]] || die "data-dir not a directory: $DATA_DIR"
case "$ASSEMBLY" in
GRCh37 | GRCh38 | both) ;;
*) die "invalid --assembly: $ASSEMBLY (GRCh37|GRCh38|both)" ;;
esac
for _n in "$FORK" "$SV_FORK" "$BUFFER_SIZE" "$PERL_CACHE_VERSION"; do
	[[ "$_n" =~ ^[0-9]+$ && "$_n" -gt 0 ]] || die "--fork, --sv-fork, --buffer-size and --cache-version take a positive integer (got '$_n')"
done
if [[ -n "$FASTA_ANY" ]]; then
	[[ "$ASSEMBLY" != "both" ]] || die "--fasta needs a single --assembly; use --fasta-grch37 / --fasta-grch38 with --assembly both"
	if [[ "$ASSEMBLY" == "GRCh37" ]]; then FASTA_GRCH37="$FASTA_ANY"; else FASTA_GRCH38="$FASTA_ANY"; fi
fi
if [[ "$DRY_RUN" != "true" ]]; then
	command -v docker >/dev/null 2>&1 || die "docker is required (pass --dry-run to preflight and print the commands without it)"
fi

[[ -n "$PERL_CACHE_DIR" ]] || PERL_CACHE_DIR="$DATA_DIR/caches/perl/vep-cache"
[[ -n "$FASTA_GRCH37" ]] || FASTA_GRCH37="$DATA_DIR/reference/grch37/genome.fa"
[[ -n "$FASTA_GRCH38" ]] || FASTA_GRCH38="$DATA_DIR/reference/grch38/genome.fa"
[[ -n "$SNP_INDEL_GT_DIR" ]] || SNP_INDEL_GT_DIR="$DATA_DIR/ground_truth/perl/snp_indel"
[[ -n "$SV_GT_DIR" ]] || SV_GT_DIR="$DATA_DIR/ground_truth/perl/sv_per_vcf"
[[ -n "$SV_INPUTS_DIR" ]] || SV_INPUTS_DIR="$SV_GT_DIR/canonical_inputs"

fasta_for() { [[ "$1" == "GRCh37" ]] && echo "$FASTA_GRCH37" || echo "$FASTA_GRCH38"; }
# `${var,,}` is bash 4+ and fails at runtime on macOS bash 3.2 while passing `bash -n`.
lower() { printf %s "$1" | tr '[:upper:]' '[:lower:]'; }

# Suite table: id:name:assembly:type:vcf_base -- the SAME rows as run_clone_measurement.sh,
# minus the sample-scaling rows, which need no Perl ground truth of their own: genotypes
# change no consequence, so the sites-only reference is the reference for both forms,
# and the genotype matrix's cost is measured separately.
SUITES="
s01:ClinVar_GRCh37:GRCh37:snp_indel:clinvar_grch37
s02:ClinVar_GRCh38:GRCh38:snp_indel:clinvar_grch38
s03:gnomAD_v2.1.1_chr21:GRCh37:snp_indel:gnomad_genomes_v2.1.1_chr21
s04:gnomAD_v4.1_chr21:GRCh38:snp_indel:gnomad_genomes_v4.1_chr21
s05:1KG_integrated_chr21:GRCh37:snp_indel:1kg_integrated_chr21_sitesonly
s06:1KG_highcov_chr21:GRCh38:snp_indel:1kg_highcov_chr21
s07:SV_GRCh37:GRCh37:sv:-
s08:SV_GRCh38:GRCh38:sv:-
"

# should_run <suite_id> <stype> <assembly>: honor --suites and --assembly. s09/s10 select
# the per-assembly SV set they point at.
should_run() {
	local sid="$1" stype="$2" asm="$3" f="$SUITES_FILTER"
	[[ "$ASSEMBLY" == "both" || "$ASSEMBLY" == "$asm" ]] || return 1
	f="${f// s09/ s07}"
	f="${f// s10/ s08}"
	f="${f/#s09/s07}"
	f="${f/#s10/s08}"
	case "$f" in
	all) return 0 ;;
	snp-indel) [[ "$stype" == "snp_indel" ]] ;;
	sv) [[ "$stype" == "sv" ]] ;;
	*) [[ " $f " == *" $sid "* ]] ;;
	esac
}

# Helpers
sha256_of() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	elif command -v shasum >/dev/null 2>&1; then
		shasum -a 256 "$1" | awk '{print $1}'
	else
		echo "unavailable"
	fi
}

# count_records <vcf|vcf.gz>: non-header lines. `n=$(...) || n=0`, never `$(... || echo 0)`:
# grep -c prints 0 AND exits 1 on no match, so the latter yields the two-line string "0\n0"
# and every later numeric test on it errors into its else branch.
count_records() {
	local n
	case "$1" in
	*.gz) n=$(gzip -dc "$1" | grep -vc '^#') || n=0 ;;
	*) n=$(grep -vc '^#' "$1") || n=0 ;;
	esac
	echo "$n"
}

count_data_lines() {
	local n
	n=$(grep -vc '^#' "$1") || n=0
	echo "$n"
}

# assert_cache_complete <assembly>: the full-genome check. Primary contigs only: Ensembl's
# full GRCh37 cache carries 699 directories, the other 674 being unplaced scaffolds and
# haplotype contigs no benchmark input reaches, so requiring all 699 would fail a
# legitimately pruned cache. Requiring the 25 primaries fails every chromosome-subset
# cache, which is the defect class.
assert_cache_complete() {
	local assembly="$1"
	local slice="$PERL_CACHE_DIR/homo_sapiens/${PERL_CACHE_VERSION}_${assembly}"
	local primaries="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
	local c missing=""
	if [[ ! -d "$slice" ]]; then
		echo "ERROR: [generate_reference_output] Perl cache slice not found: $slice" >&2
		echo "  Expected homo_sapiens/${PERL_CACHE_VERSION}_${assembly} under --perl-cache-dir ($PERL_CACHE_DIR)." >&2
		echo "  See the header of this script for the Ensembl FTP location of the full-genome cache." >&2
		return 1
	fi
	for c in $primaries; do
		[[ -d "$slice/$c" ]] || missing="$missing $c"
	done
	if [[ -n "$missing" ]]; then
		echo "ERROR: [generate_reference_output] Perl cache for $assembly is missing primary contigs:$missing" >&2
		echo "  Slice: $slice" >&2
		echo "  A chromosome-subset cache does not fail a VEP run; it under-annotates every record" >&2
		echo "  outside the contigs it holds and exits 0. A reference generated against it is wrong" >&2
		echo "  in a way no comparator can see. Extract the full-genome Ensembl cache first." >&2
		return 1
	fi
	log "perl cache check ($assembly): all 25 primary contigs present under $slice"
	return 0
}

assert_fasta() {
	local assembly="$1" fasta
	fasta=$(fasta_for "$assembly")
	[[ -f "$fasta" ]] || {
		echo "ERROR: [generate_reference_output] FASTA for $assembly not found: $fasta" >&2
		echo "  Pass --fasta-grch37/--fasta-grch38 (or --fasta with a single --assembly)." >&2
		return 1
	}
	[[ -f "${fasta}.fai" ]] || {
		echo "ERROR: [generate_reference_output] FASTA index not found: ${fasta}.fai" >&2
		echo "  VEP reads the index; build it with \`samtools faidx\` before generating." >&2
		return 1
	}
	return 0
}

# docker_cmd <assembly> <fork> <in.vcf[.gz]> <out_dir> <out_name>
# Prints the exact docker invocation as one shell-quoted line. This is the harness's
# Perl invocation, flag for flag, so the reference and the timed runs cannot diverge.
docker_cmd() {
	local assembly="$1" fork="$2" in_vcf="$3" out_dir="$4" out_name="$5"
	local fasta in_mount
	fasta=$(fasta_for "$assembly")
	case "$in_vcf" in
	*.gz) in_mount="/work/input.vcf.gz" ;;
	*) in_mount="/work/input.vcf" ;;
	esac
	printf '%q ' docker run --rm \
		--user "$(id -u):$(id -g)" \
		-v "${in_vcf}:${in_mount}:ro" \
		-v "${out_dir}:/work/perl_out" \
		-v "${PERL_CACHE_DIR}:/work/cache:ro" \
		-v "${fasta}:/work/ref.fa:ro" \
		-v "${fasta}.fai:/work/ref.fa.fai:ro" \
		"$PERL_IMAGE" \
		vep \
		-i "$in_mount" \
		-o "/work/perl_out/${out_name}" \
		--offline --cache --dir_cache /work/cache --cache_version "$PERL_CACHE_VERSION" \
		--species homo_sapiens --assembly "$assembly" \
		--format vcf --buffer_size "$BUFFER_SIZE" --fork "$fork" \
		--fasta /work/ref.fa \
		--no_headers --no_stats --quiet --force_overwrite
	echo
}

# write_provenance <json_path> <assembly> <fork> <in> <in_records> <out> <out_lines> <wall> <image_digest>
write_provenance() {
	local json="$1" assembly="$2" fork="$3" in_vcf="$4" in_records="$5" out="$6" out_lines="$7" wall="$8" digest="$9"
	local digest_json
	if [[ -n "$digest" ]]; then digest_json="\"$digest\""; else digest_json="null"; fi
	cat >"$json" <<EOF
{
  "generated_at_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "generator": "scripts/concordance/generate_reference_output.sh",
  "perl_image": "$PERL_IMAGE",
  "perl_image_digest": $digest_json,
  "assembly": "$assembly",
  "cache_version": $PERL_CACHE_VERSION,
  "cache_dir": "$PERL_CACHE_DIR",
  "fasta": "$(fasta_for "$assembly")",
  "fork": $fork,
  "buffer_size": $BUFFER_SIZE,
  "vep_flags": "--offline --cache --dir_cache <cache> --cache_version $PERL_CACHE_VERSION --species homo_sapiens --assembly $assembly --format vcf --buffer_size $BUFFER_SIZE --fork $fork --fasta <fa> --no_headers --no_stats --quiet --force_overwrite",
  "input": "$in_vcf",
  "input_sha256": "$(sha256_of "$in_vcf")",
  "input_records": $in_records,
  "output": "$out",
  "output_sha256": "$(sha256_of "$out")",
  "output_data_lines": $out_lines,
  "wall_seconds": $wall
}
EOF
}

# Plan: one line per output, "kind|assembly|fork|out_dir|input|out_name". Built in full
# BEFORE anything runs, so every missing input, FASTA, cache slice and output collision is
# reported together and nothing is generated against a half-valid tree. Lines sharing
# their first four fields form a GROUP (one SNP/indel suite, or one assembly's 16-file SV
# set), which is the unit the output-volume gate below is applied to. A newline-delimited
# string rather than an array of records keeps this bash-3.2 safe.
JOBS=""
PREFLIGHT_FAILED=0
CHECKED_ASM=""

preflight_assembly() {
	local asm="$1"
	[[ " $CHECKED_ASM " == *" $asm "* ]] && return 0
	CHECKED_ASM="$CHECKED_ASM $asm"
	assert_fasta "$asm" || PREFLIGHT_FAILED=1
	assert_cache_complete "$asm" || PREFLIGHT_FAILED=1
}

# check_output_collision <path>
check_output_collision() {
	if [[ -e "$1" && "$FORCE" != "true" ]]; then
		echo "ERROR: [generate_reference_output] reference output already exists: $1" >&2
		echo "  A reference is a pinned artifact. Pass --force to overwrite it." >&2
		return 1
	fi
	return 0
}

n_selected=0
while IFS= read -r entry; do
	[[ -n "$entry" ]] || continue
	IFS=':' read -r suite_id suite_name asm stype vcf_base <<<"$entry"
	should_run "$suite_id" "$stype" "$asm" || continue
	n_selected=$((n_selected + 1))
	preflight_assembly "$asm"
	asm_lc=$(lower "$asm")
	if [[ "$stype" == "snp_indel" ]]; then
		in_vcf="$DATA_DIR/vcf/$asm_lc/all_variants_canonical/${vcf_base}.canonical.vcf.gz"
		if [[ ! -f "$in_vcf" ]]; then
			echo "ERROR: [generate_reference_output] $suite_id ($suite_name): input not found: $in_vcf" >&2
			echo "  Canonical inputs carry Ensembl-style contigs; derive them with prepare_benchmark_vcfs.py." >&2
			PREFLIGHT_FAILED=1
			continue
		fi
		check_output_collision "$SNP_INDEL_GT_DIR/$suite_id/output.txt" || PREFLIGHT_FAILED=1
		JOBS="${JOBS}snp_indel|$asm|$FORK|$SNP_INDEL_GT_DIR/$suite_id|$in_vcf|output.txt
"
	else
		in_dir="$SV_INPUTS_DIR/$asm_lc"
		if [[ ! -d "$in_dir" ]]; then
			echo "ERROR: [generate_reference_output] $suite_id ($suite_name): SV input directory not found: $in_dir" >&2
			PREFLIGHT_FAILED=1
			continue
		fi
		n_sv=0
		for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
			[[ -f "$vcf" ]] || continue
			b=$(basename "$vcf")
			b=${b%.vcf.gz}
			b=${b%.vcf}
			[[ "$b" == "manifest" ]] && continue
			n_sv=$((n_sv + 1))
			check_output_collision "$SV_GT_DIR/$asm_lc/$b.txt" || PREFLIGHT_FAILED=1
			JOBS="${JOBS}sv|$asm|$SV_FORK|$SV_GT_DIR/$asm_lc|$vcf|$b.txt
"
		done
		if [[ "$n_sv" -eq 0 ]]; then
			echo "ERROR: [generate_reference_output] $suite_id ($suite_name): no .vcf or .vcf.gz under $in_dir" >&2
			echo "  An empty input directory would produce an empty reference that scores as a clean cell." >&2
			PREFLIGHT_FAILED=1
		fi
	fi
done <<<"$SUITES"

[[ "$n_selected" -gt 0 ]] || die "no suite matched --suites '$SUITES_FILTER' with --assembly $ASSEMBLY"
[[ "$PREFLIGHT_FAILED" -eq 0 ]] || die "preflight failed; nothing was generated"
[[ -n "$JOBS" ]] || die "preflight selected $n_selected suite(s) but planned no output"

n_jobs=$(printf '%s' "$JOBS" | grep -c '|') || n_jobs=0
log "generate_reference_output image=$PERL_IMAGE cache_version=$PERL_CACHE_VERSION fork=$FORK sv_fork=$SV_FORK buffer_size=$BUFFER_SIZE"
log "data_dir=$DATA_DIR snp_indel_gt=$SNP_INDEL_GT_DIR sv_gt=$SV_GT_DIR planned_outputs=$n_jobs dry_run=$DRY_RUN"

# Dry run: print every command and stop.
if [[ "$DRY_RUN" == "true" ]]; then
	while IFS='|' read -r kind asm fork out_dir in_vcf out_name; do
		[[ -n "$kind" ]] || continue
		echo "DRY-RUN [$kind $asm] -> $out_dir/$out_name"
		docker_cmd "$asm" "$fork" "$in_vcf" "$out_dir" "$out_name"
	done <<<"$JOBS"
	log "dry run complete: $n_jobs command(s) printed, nothing written"
	exit 0
fi

# Run, one group at a time. Every output of a group is produced in one scratch directory
# beside its destination, the group is gated on annotation volume, and only then is
# anything moved into place. An interrupted or non-annotating run therefore never leaves
# a partial file at a path the harness would trust.
#
# THE OUTPUT-VOLUME GATE is the harness's own floor, applied at the harness's own grain:
# summed rows >= summed input records, per SNP/indel suite and per assembly's SV set. A
# VEP that cannot reach its cache emits one warning per record, zero annotation rows and
# exit 0; a working VEP emits at least one row per input record (one per overlapping
# transcript, or a single intergenic row) and several per record on genome-wide and SV
# inputs alike. The floor is NOT applied per SV file: some synthetic
# files consist of allele classes VEP declines to annotate (`<*>`, `<NON_REF>`, `.`), so a
# legitimately sparse file would fail a per-file floor while the aggregate still bounds
# the failure mode that matters. A file with zero rows is reported, never silently kept.
IMAGE_DIGEST=$(docker image inspect --format '{{index .RepoDigests 0}}' "$PERL_IMAGE" 2>/dev/null || true)
if [[ -z "$IMAGE_DIGEST" ]]; then
	log "pulling $PERL_IMAGE"
	docker pull "$PERL_IMAGE" >/dev/null
	IMAGE_DIGEST=$(docker image inspect --format '{{index .RepoDigests 0}}' "$PERL_IMAGE" 2>/dev/null || true)
fi
[[ -n "$IMAGE_DIGEST" ]] && log "image digest $IMAGE_DIGEST"

done_jobs=0
GROUP_KEYS=$(printf '%s' "$JOBS" | cut -d'|' -f1-4 | uniq)
while IFS='|' read -r kind asm fork out_dir; do
	[[ -n "$kind" ]] || continue
	group_prefix="$kind|$asm|$fork|$out_dir|"
	mkdir -p "$out_dir"
	scratch=$(mktemp -d "$out_dir/.generate.XXXXXX")
	group_records=0
	group_rows=0
	group_files=0
	# Per-file facts for the provenance sidecars, recorded before the gate so a failure
	# discards them with the scratch directory. "in_vcf|records|rows|wall|out_name".
	group_facts=""
	log "group [$kind $asm --fork $fork] -> $out_dir/"
	while IFS='|' read -r _k _a _f _o in_vcf out_name; do
		[[ -n "$_k" ]] || continue
		in_records=$(count_records "$in_vcf")
		if [[ "$in_records" -eq 0 ]]; then
			rm -rf "$scratch"
			die "$out_dir/$out_name: input $in_vcf has 0 records; an empty reference would score as a clean cell"
		fi
		log "  $in_vcf ($in_records records) -> $out_name"
		start=$(date +%s)
		set +e
		eval "$(docker_cmd "$asm" "$fork" "$in_vcf" "$scratch" "$out_name")"
		rc=$?
		set -e
		wall=$(($(date +%s) - start))
		if [[ "$rc" -ne 0 ]]; then
			rm -rf "$scratch"
			die "$out_dir/$out_name: vep exited $rc after ${wall}s; the group's partial output was discarded"
		fi
		[[ -f "$scratch/$out_name" ]] || {
			rm -rf "$scratch"
			die "$out_dir/$out_name: vep exited 0 but wrote no output; the group's partial output was discarded"
		}
		out_lines=$(count_data_lines "$scratch/$out_name")
		[[ "$out_lines" -gt 0 ]] || echo "WARN: [generate_reference_output] $out_dir/$out_name has 0 annotation rows for $in_records records" >&2
		group_records=$((group_records + in_records))
		group_rows=$((group_rows + out_lines))
		group_files=$((group_files + 1))
		group_facts="${group_facts}${in_vcf}|${in_records}|${out_lines}|${wall}|${out_name}
"
		log "    $out_lines rows in ${wall}s"
	done <<<"$(printf '%s' "$JOBS" | awk -v p="$group_prefix" 'index($0, p) == 1')"

	if [[ "$group_rows" -lt "$group_records" ]]; then
		echo "ERROR: [generate_reference_output] $out_dir: $group_rows annotation rows for $group_records input records across $group_files file(s)" >&2
		echo "  A VEP run that cannot reach its cache emits warnings rather than annotations and still exits 0." >&2
		echo "  Check --perl-cache-dir, --cache-version and --assembly. The group's output was discarded." >&2
		for w in "$scratch"/*_warnings.txt; do
			[[ -f "$w" ]] || continue
			head -3 "$w" | sed 's/^/  warnings: /' >&2
			break
		done
		rm -rf "$scratch"
		exit 1
	fi

	while IFS='|' read -r in_vcf in_records out_lines wall out_name; do
		[[ -n "$out_name" ]] || continue
		mv -f "$scratch/$out_name" "$out_dir/$out_name"
		if [[ -f "$scratch/${out_name}_warnings.txt" ]]; then
			mv -f "$scratch/${out_name}_warnings.txt" "$out_dir/${out_name%.txt}.warnings.log"
		fi
		write_provenance "$out_dir/${out_name%.txt}.provenance.json" "$asm" "$fork" "$in_vcf" "$in_records" \
			"$out_dir/$out_name" "$out_lines" "$wall" "$IMAGE_DIGEST"
		done_jobs=$((done_jobs + 1))
		log "  wrote $out_dir/$out_name"
	done <<<"$group_facts"
	rm -rf "$scratch"
	log "group [$kind $asm] complete: $group_rows rows / $group_records records over $group_files file(s)"
done <<<"$GROUP_KEYS"

log "generate_reference_output complete: $done_jobs/$n_jobs output(s) written"
