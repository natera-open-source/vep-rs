#!/usr/bin/env bash
#
# generate_reference_output.sh -- write the Ensembl VEP reference outputs that
# run_clone_measurement.sh scores every engine against, a flag family's reference
# in any output format, or a golden corpus's expected outputs.
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
# OUTPUT FORMATS. Ensembl VEP builds one OutputFactory per run from its output_format
# (Runner.pm:871 get_OutputFactory, OutputFactory.pm:97-102 FORMAT_MAP), so a reference
# in several formats is one container run per format: --format default|tab|vcf|json|all.
# In the suite layout the default format keeps the names above and the other formats
# take their own extension on the same stem (output.tsv, output.vcf, output.jsonl;
# <vcf-basename>.tsv, .vcf, .jsonl), so a `*.txt` glob still pairs only the default
# output. --no_headers is passed with the default format alone: the tab output's column
# line and the VCF output's `##INFO=<ID=CSQ...>` line are what name the columns of those
# two formats, and the JSON format has no header to drop.
#
# FLAG FAMILIES. --flags "<words>" appends the words, split on whitespace, after the
# harness's flags and before the format flag, so a reference for one Ensembl VEP flag
# family is the harness's invocation plus that family. A word this script already sets
# (-i, -o, --fork, --fasta, --tab, --no_stats, ...) is refused rather than passed twice,
# in the single-dash and abbreviated spellings vep's option parser accepts as well. A
# family carrying --nearest mounts the Ensembl cache writable: Ensembl VEP writes
# transcript_gene_tss.txt into the cache directory on its first --nearest run and keeps
# the file, built under that run's transcript filters, for every later run; a run on a
# cache that lacks the file scans every region file once, so the first --nearest run on
# a cache goes alone.
#
# STATS. --stats keeps Ensembl VEP's statistics report beside each output as
# <stem>.summary.log instead of passing --no_stats. The run passes --stats_text, because
# a run that names neither --stats_text nor --stats_html writes the HTML report alone
# (Runner.pm:787-789); the text file is <output>_summary.txt (BaseRunner.pm:254-271) and
# is renamed out of the `*.txt` namespace for the reason the warnings file is.
#
# GOLDEN CORPORA. --corpus <dir> writes a corpus's expected outputs in the layout
# crates/vep-cli/tests/common/mod.rs reads: expected/{default.txt,tab.txt,vcf.vcf,
# json.jsonl}.gz from <dir>/variants.vcf, with the manifest's `assembly`, `flags` and
# `fasta`, --fork 1, headers kept, and <dir>/provenance.json carrying the image digest,
# the vep command and the sha256 of every expected file's uncompressed bytes. A key the
# script does not own (the divergence classification, the pruned-cache summary) survives
# a regeneration, and `input` and `fasta` are filled only when absent, so a record's own
# description of them stands. --format must be `all` for a corpus: the four files are one
# contract.
#
# THE IMAGE DIGEST. Before anything runs, the pulled image's RepoDigests must carry
# --perl-image-digest (default the release_116.2 tag's manifest-list digest, or
# VEP_PERL_IMAGE_DIGEST), so a moved tag cannot produce a reference under this name.
# --perl-image and --perl-image-digest move together.
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
#   Corpus:    --fork 1 --buffer_size 5000, the settings the committed corpora carry.
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
# release, 116:
#   GRCh38: https://ftp.ensembl.org/pub/release-<v>/variation/indexed_vep_cache/homo_sapiens_vep_<v>_GRCh38.tar.gz
#   GRCh37: https://ftp.ensembl.org/pub/grch37/release-<v>/variation/indexed_vep_cache/homo_sapiens_vep_<v>_GRCh37.tar.gz
#   Extract both under ONE directory, so that it contains
#   homo_sapiens/116_GRCh38/ and homo_sapiens/116_GRCh37/; that directory is
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
# A corpus passes --fasta only when its manifest names one, as the golden tests do.
#
# This script never fetches anything, needs no cloud credentials, and refuses to
# overwrite an existing reference output unless --force is passed: a reference is a
# pinned artifact whose digests the provenance records carry.
#
# Usage:
#   scripts/concordance/generate_reference_output.sh --data-dir <data> \
#     [--assembly GRCh37|GRCh38|both] [--suites all|snp-indel|sv|"s01 s07"] \
#     [--format default|tab|vcf|json|all] [--flags "<vep flags>"] [--stats] \
#     [--volume-floor records|any] \
#     [--perl-image ensemblorg/ensembl-vep:release_116.2] [--perl-image-digest sha256:<64 hex>] \
#     [--perl-cache-dir <dir>] [--cache-version 116] \
#     [--fasta <fa>] [--fasta-grch37 <fa>] [--fasta-grch38 <fa>] \
#     [--fork 16] [--sv-fork 4] [--buffer-size 5000] \
#     [--snp-indel-gt-dir <dir>] [--sv-gt-dir <dir>] [--sv-inputs-dir <dir>] \
#     [--force] [--dry-run]
#   scripts/concordance/generate_reference_output.sh --corpus tests/golden/116/GRCh38 \
#     --perl-cache-dir <dir> [--stats] [--volume-floor records|any] [--force] [--dry-run]
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
FORMAT=""
EXTRA_FLAGS=""
KEEP_STATS="false"
VOLUME_FLOOR="records"
CORPUS_DIR=""
PERL_IMAGE="ensemblorg/ensembl-vep:release_116.2"
# The tag's manifest-list digest on Docker Hub; the pulled image's RepoDigests must carry it
# before anything is generated, so a moved tag cannot produce a reference under this name.
PERL_IMAGE_DIGEST="${VEP_PERL_IMAGE_DIGEST:-sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da}"
PERL_CACHE_DIR="${VEP_PERL_CACHE_DIR:-}"
PERL_CACHE_VERSION="${VEP_PERL_CACHE_VERSION:-116}"
FASTA_ANY=""
FASTA_GRCH37="${VEP_FASTA_GRCH37:-}"
FASTA_GRCH38="${VEP_FASTA_GRCH38:-}"
FORK=""
SV_FORK=4
BUFFER_SIZE=5000
SNP_INDEL_GT_DIR="${VEP_SNP_INDEL_GT_DIR:-}"
SV_GT_DIR="${VEP_SV_GT_DIR:-}"
SV_INPUTS_DIR=""
FORCE="false"
DRY_RUN="false"
# Options that select or place suite outputs, recorded as seen so --corpus can refuse them.
SUITE_OPTS_SEEN=""

usage() {
	cat <<'EOF'
Generate the Ensembl VEP reference outputs run_clone_measurement.sh scores against, a
flag family's reference in any output format, or a golden corpus's expected outputs.

Required (one of):
  --data-dir <path>        Data tree (vcf/, caches/, reference/, ground_truth/)
  --corpus <dir>           A golden corpus (variants.vcf, manifest.json): writes its
                           expected/{default.txt,tab.txt,vcf.vcf,json.jsonl}.gz and
                           provenance.json; needs --perl-cache-dir or --data-dir

Optional:
  --assembly <A>           GRCh37 | GRCh38 | both (default both)
  --suites <filter>        all (default) | snp-indel | sv | space-separated "s01 s07"
                           (s09/s10 are aliases of s07/s08: the real-world SV VCFs are
                           part of the 16-file per-assembly set)
  --format <F>             Output format: default | tab | vcf | json | all, one
                           container run per format (default: default; a corpus is
                           always all). Ensembl VEP's own --format vcf names its input.
  --flags "<words>"        Ensembl VEP flags appended after the harness's, split on
                           whitespace; a flag this script sets itself is refused
  --stats                  Keep the statistics report as <stem>.summary.log (passes
                           --stats_text) instead of --no_stats
  --volume-floor <mode>    records (default): annotated rows >= input records per
                           suite or SV set, and every written VCF record or JSON object
                           annotated; any: at least one annotated row, for flags that
                           drop records or rows (--chr, --coding_only, --no_intergenic,
                           --filter_common, --allow_non_variant, --individual,
                           --individual_zyg, --dont_skip) or write one row per record
                           (--pick, --most_severe, --summary) on an input with records
                           vep skips
  --perl-image <ref>       Docker image (default ensemblorg/ensembl-vep:release_116.2)
  --perl-image-digest <sha256:hex>  Digest the pulled image must carry (default the
                           pinned tag's manifest-list digest, or VEP_PERL_IMAGE_DIGEST)
  --perl-cache-dir <dir>   Directory holding homo_sapiens/<ver>_<asm>/
                           (default <data>/caches/perl/vep-cache)
  --cache-version <N>      Ensembl cache release (default 116)
  --fasta <fa>             Reference FASTA; valid only with a single --assembly
  --fasta-grch37 <fa>      Per-assembly FASTA overrides (defaults under
  --fasta-grch38 <fa>      <data>/reference/<asm>/)
  --fork <N>               SNP/indel --fork (default 16; a corpus defaults to 1)
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
	--corpus)
		CORPUS_DIR="$2"
		shift 2
		;;
	--assembly)
		ASSEMBLY="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--suites)
		SUITES_FILTER="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--format)
		FORMAT="$2"
		shift 2
		;;
	--flags)
		EXTRA_FLAGS="$2"
		shift 2
		;;
	--stats)
		KEEP_STATS="true"
		shift
		;;
	--volume-floor)
		VOLUME_FLOOR="$2"
		shift 2
		;;
	--perl-image)
		PERL_IMAGE="$2"
		shift 2
		;;
	--perl-image-digest)
		PERL_IMAGE_DIGEST="$2"
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
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--fasta-grch37)
		FASTA_GRCH37="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--fasta-grch38)
		FASTA_GRCH38="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--sv-fork)
		SV_FORK="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--buffer-size)
		BUFFER_SIZE="$2"
		shift 2
		;;
	--snp-indel-gt-dir)
		SNP_INDEL_GT_DIR="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--sv-gt-dir)
		SV_GT_DIR="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--sv-inputs-dir)
		SV_INPUTS_DIR="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
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
if [[ -n "$CORPUS_DIR" ]]; then
	[[ -d "$CORPUS_DIR" ]] || die "corpus not a directory: $CORPUS_DIR"
	# Docker bind mounts take absolute host paths; the documented call names the corpus
	# relative to the repository root.
	CORPUS_DIR=$(cd "$CORPUS_DIR" && pwd)
	[[ -f "$CORPUS_DIR/variants.vcf" && -f "$CORPUS_DIR/manifest.json" ]] || die "corpus needs variants.vcf and manifest.json: $CORPUS_DIR"
	[[ -z "$SUITE_OPTS_SEEN" ]] || die "--corpus takes its assembly, flags, FASTA and output location from the corpus; drop$SUITE_OPTS_SEEN"
	[[ -z "$EXTRA_FLAGS" ]] || die "--corpus takes its flags from manifest.json; drop --flags"
	[[ -z "$FORMAT" || "$FORMAT" == "all" ]] || die "--corpus writes all four formats (--format all); a corpus with one format regenerated is not a corpus"
	FORMAT="all"
	[[ -n "$FORK" ]] || FORK=1
	[[ -n "$DATA_DIR" || -n "$PERL_CACHE_DIR" ]] || die "--corpus needs --perl-cache-dir (or --data-dir for its default cache location)"
	command -v python3 >/dev/null 2>&1 || die "python3 is required to read manifest.json and write provenance.json"
else
	[[ -n "$DATA_DIR" ]] || {
		echo "ERROR: [generate_reference_output] --data-dir or --corpus required" >&2
		usage
	}
	[[ -n "$FORMAT" ]] || FORMAT="default"
	[[ -n "$FORK" ]] || FORK=16
fi
[[ -z "$DATA_DIR" || -d "$DATA_DIR" ]] || die "data-dir not a directory: $DATA_DIR"
case "$ASSEMBLY" in
GRCh37 | GRCh38 | both) ;;
*) die "invalid --assembly: $ASSEMBLY (GRCh37|GRCh38|both)" ;;
esac
case "$FORMAT" in
default | tab | vcf | json | all) ;;
*) die "invalid --format: $FORMAT (default|tab|vcf|json|all)" ;;
esac
case "$VOLUME_FLOOR" in
records | any) ;;
*) die "invalid --volume-floor: $VOLUME_FLOOR (records|any)" ;;
esac
for _n in "$FORK" "$SV_FORK" "$BUFFER_SIZE" "$PERL_CACHE_VERSION"; do
	[[ "$_n" =~ ^[0-9]+$ && "$_n" -gt 0 ]] || die "--fork, --sv-fork, --buffer-size and --cache-version take a positive integer (got '$_n')"
done
[[ "$PERL_IMAGE_DIGEST" =~ ^sha256:[0-9a-f]{64}$ ]] || die "--perl-image-digest must be sha256:<64 hex> (got '$PERL_IMAGE_DIGEST')"
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

# The flags this script passes itself. A second copy inside --flags would win or lose
# by position in Ensembl VEP's option parsing, so it is refused instead. vep reads its
# options with Getopt::Long's defaults (vep:27, Config.pm:62 at release/116.2), which
# take a single dash or a plus sign, any unambiguous abbreviation and any letter case as
# the full option, so a word is refused when the lower-cased option name it carries is a
# prefix of one of these.
OWN_FLAGS="--input_file --output_file --offline --cache --dir_cache --cache_version --species --assembly --format --buffer_size --fork --fasta --no_headers --no_stats --stats_text --stats_html --stats_file --quiet --force_overwrite --tab --vcf --json --output_format"
# The INFO key Ensembl VEP writes its consequences under in the VCF format, the
# `--vcf_info_field` of the flag family or CSQ (OutputFactory/VCF.pm:312).
VCF_INFO_FIELD="CSQ"
# option_name <word>: the option named by a vep argument word, in its lower-case `--name`
# spelling with a value suffix removed, or nothing for a value word. `sf` is the one alias
# of an owned option (stats_file|sf) that is not a prefix of the option's own name.
option_name() {
	local name="${1%%=*}"
	case "$name" in
	--[A-Za-z]*) name=$(lower "$name") ;;
	-[A-Za-z]* | +[A-Za-z]*) name=$(lower "--${name#[-+]}") ;;
	*) return 0 ;;
	esac
	[[ "$name" != "--sf" ]] || name="--stats_file"
	printf '%s' "$name"
}
check_extra_flags() {
	local word name own take_field=""
	set -f
	# shellcheck disable=SC2086
	set -- $EXTRA_FLAGS
	set +f
	# One space between words from here on: docker_cmd splits the rendered flags on spaces
	# alone, so a tab or newline inside --flags would reach vep as part of a word or end
	# the command early, and the provenance records would carry it verbatim.
	EXTRA_FLAGS="$*"
	for word in "$@"; do
		if [[ -n "$take_field" ]]; then
			VCF_INFO_FIELD="$word"
			take_field=""
			continue
		fi
		name=$(option_name "$word")
		[[ -n "$name" ]] || continue
		# An exact option or alias name outranks an abbreviation in Getopt::Long, and three
		# vep names are prefixes of owned flags: `verbose|v` (Config.pm:159), `terms|t`
		# (:176) and `dir` (:241); vep reads them as themselves, so they pass.
		case "$name" in
		--v | --t | --dir) continue ;;
		esac
		for own in $OWN_FLAGS; do
			if [[ "$own" == "$name"* ]]; then
				die "--flags carries $word, which vep reads as $own, a flag this script sets; use its own option (--format, --fork, --stats, --fasta, --buffer-size, --cache-version) or drop it"
			fi
		done
		# `--vcf` and its prefixes were refused above, so a word that is still a prefix of
		# --vcf_info_field names that option alone.
		if [[ "--vcf_info_field" == "$name"* ]]; then
			if [[ "$word" == *=* ]]; then VCF_INFO_FIELD="${word#*=}"; else take_field=1; fi
		fi
	done
	[[ -z "$take_field" ]] || die "--flags ends with --vcf_info_field and no field name"
}

# A corpus takes its assembly, flags and FASTA from manifest.json. The flags are a list
# of words, each one argument to vep, so a word holding whitespace cannot be expressed.
# An existing provenance.json is parsed here too, so a damaged record stops the run
# before the container runs rather than after them.
CORPUS_FASTA_NAME=""
if [[ -n "$CORPUS_DIR" ]]; then
	manifest_fields=$(python3 - "$CORPUS_DIR/manifest.json" 2>&1 <<'EOF'
import json, sys
from pathlib import Path
manifest_path = Path(sys.argv[1])


def load(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except ValueError as e:
        sys.exit(f"{path}: {e}")


m = load(manifest_path)
record_path = manifest_path.with_name("provenance.json")
if record_path.is_file() and not isinstance(load(record_path), dict):
    sys.exit(f"{record_path}: the record is not a JSON object")
flags = m.get("flags") or []
if not isinstance(flags, list):
    sys.exit(f"manifest flags must be a list of words, got {flags!r}")
bad = [f for f in flags if not isinstance(f, str) or f != f.strip() or any(c.isspace() for c in f)]
if bad:
    sys.exit(f"manifest flags must be single words: {bad!r}")
assembly = m.get("assembly") or ""
if assembly not in ("GRCh37", "GRCh38"):
    sys.exit(f"manifest assembly must be GRCh37 or GRCh38, got {assembly!r}")
print(assembly)
print(m.get("fasta") or "")
print(" ".join(flags))
EOF
	) || die "corpus: $manifest_fields"
	ASSEMBLY=$(printf '%s\n' "$manifest_fields" | sed -n 1p)
	CORPUS_FASTA_NAME=$(printf '%s\n' "$manifest_fields" | sed -n 2p)
	EXTRA_FLAGS=$(printf '%s\n' "$manifest_fields" | sed -n 3p)
fi
check_extra_flags

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

# The formats --format selects, in the order they run.
formats_selected() {
	case "$FORMAT" in
	all) echo "default tab vcf json" ;;
	*) echo "$FORMAT" ;;
	esac
}

# out_name_for <kind> <base> <format>: the output file name. A suite output keeps the
# harness's name for the default format and takes the format's extension otherwise; a
# corpus output takes the names the golden tests read.
out_name_for() {
	local kind="$1" base="$2" format="$3"
	if [[ "$kind" == "corpus" ]]; then
		case "$format" in
		default) echo "default.txt" ;;
		tab) echo "tab.txt" ;;
		vcf) echo "vcf.vcf" ;;
		json) echo "json.jsonl" ;;
		esac
		return
	fi
	case "$format" in
	default) echo "$base.txt" ;;
	tab) echo "$base.tsv" ;;
	vcf) echo "$base.vcf" ;;
	json) echo "$base.jsonl" ;;
	esac
}

# sidecar_stem <kind> <format> <out_name>: what the provenance, warnings and summary
# sidecars are named after: the output name less a default-format `.txt`, or the format
# name for a corpus output.
sidecar_stem() {
	if [[ "$1" == "corpus" ]]; then echo "$2"; else echo "${3%.txt}"; fi
}

# headers_for <kind> <format>: whether the run keeps Ensembl VEP's header lines. The
# harness reads the default suite output without them; every other output keeps them,
# because the tab column line and the VCF `##INFO` line name the columns.
headers_for() {
	if [[ "$1" != "corpus" && "$2" == "default" ]]; then echo "no"; else echo "yes"; fi
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

json_escape() { printf '%s' "$1" | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g'; }

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

# annotated_units <format> <file>: what the output-volume gate counts. The default and tab
# formats write one row per consequence. The VCF format writes every buffered record as a
# line and appends the CSQ field only to a record with consequences (OutputFactory/VCF.pm:
# 299-361), and the JSON format writes one object per record with a `*_consequences` array
# only on an annotated one (OutputFactory/JSON.pm:213-260, :384), so for those two the unit
# is an annotated record.
annotated_units() {
	local n
	case "$1" in
	default | tab) n=$(count_data_lines "$2") ;;
	vcf) n=$(FIELD="$VCF_INFO_FIELD" awk -F'\t' '!/^#/ && index(";" $8, ";" ENVIRON["FIELD"] "=") {c++} END {print c + 0}' "$2") ;;
	json) n=$(grep -c '_consequences":' "$2") || n=0 ;;
	esac
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

# vep_args <assembly> <fork> <format> <cache> <fasta|""> <headers yes|no>
# Everything after `vep -i <in> -o <out>`, one shell word each, space-joined: the
# harness's invocation, --fasta when one is given, --no_headers when the headers go, the
# stats flag, --quiet --force_overwrite, the --flags words and the output-format flag
# (none for the default format). docker_cmd and the provenance records both render this
# one function, so the command run and the command recorded cannot drift.
vep_args() {
	local assembly="$1" fork="$2" format="$3" cache="$4" fasta="$5" headers="$6"
	local a="--offline --cache --dir_cache $cache --cache_version $PERL_CACHE_VERSION --species homo_sapiens --assembly $assembly --format vcf --buffer_size $BUFFER_SIZE --fork $fork"
	if [[ -n "$fasta" ]]; then a="$a --fasta $fasta"; fi
	if [[ "$headers" == "no" ]]; then a="$a --no_headers"; fi
	if [[ "$KEEP_STATS" == "true" ]]; then a="$a --stats_text"; else a="$a --no_stats"; fi
	a="$a --quiet --force_overwrite"
	if [[ -n "$EXTRA_FLAGS" ]]; then a="$a $EXTRA_FLAGS"; fi
	case "$format" in
	tab | vcf | json) a="$a --$format" ;;
	esac
	printf '%s' "$a"
}

# docker_cmd <assembly> <fork> <format> <in.vcf[.gz]> <out_dir> <out_name> <fasta|""> <headers>
# Prints the exact docker invocation as one shell-quoted line. This is the harness's
# Perl invocation, flag for flag, so the reference and the timed runs cannot diverge.
docker_cmd() {
	local assembly="$1" fork="$2" format="$3" in_vcf="$4" out_dir="$5" out_name="$6" fasta="$7" headers="$8"
	local in_mount fasta_in=""
	local -a argv
	case "$in_vcf" in
	*.gz) in_mount="/work/input.vcf.gz" ;;
	*) in_mount="/work/input.vcf" ;;
	esac
	argv=(docker run --rm
		--user "$(id -u):$(id -g)"
		-v "${in_vcf}:${in_mount}:ro"
		-v "${out_dir}:/work/perl_out")
	# Ensembl VEP's first --nearest run writes transcript_gene_tss.txt into the cache
	# directory and dies on a read-only one (AnnotationSource/Cache/Transcript.pm:288-289).
	if [[ " $EXTRA_FLAGS " == *" --nearest "* || " $EXTRA_FLAGS " == *" --nearest="* ]]; then
		argv+=(-v "${PERL_CACHE_DIR}:/work/cache")
	else
		argv+=(-v "${PERL_CACHE_DIR}:/work/cache:ro")
	fi
	if [[ -n "$fasta" ]]; then
		argv+=(-v "${fasta}:/work/ref.fa:ro" -v "${fasta}.fai:/work/ref.fa.fai:ro")
		fasta_in="/work/ref.fa"
	fi
	argv+=("$PERL_IMAGE" vep -i "$in_mount" -o "/work/perl_out/${out_name}")
	local -a words
	IFS=' ' read -r -a words <<<"$(vep_args "$assembly" "$fork" "$format" /work/cache "$fasta_in" "$headers")"
	argv+=("${words[@]}")
	printf '%q ' "${argv[@]}"
	echo
}

# write_provenance <json_path> <assembly> <fork> <format> <headers> <in> <in_records> <out> <out_units> <out_lines> <wall> <image_digest>
write_provenance() {
	local json="$1" assembly="$2" fork="$3" format="$4" headers="$5" in_vcf="$6" in_records="$7" out="$8" out_units="$9" out_lines="${10}" wall="${11}" digest="${12}"
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
  "output_format": "$format",
  "vep_flags": "$(json_escape "$(vep_args "$assembly" "$fork" "$format" '<cache>' '<fa>' "$headers")")",
  "input": "$in_vcf",
  "input_sha256": "$(sha256_of "$in_vcf")",
  "input_records": $in_records,
  "output": "$out",
  "output_sha256": "$(sha256_of "$out")",
  "output_annotated_units": $out_units,
  "output_data_lines": $out_lines,
  "wall_seconds": $wall
}
EOF
}

# Plan: one line per output, "kind|assembly|fork|format|out_dir|input|out_name". Built in
# full BEFORE anything runs, so every missing input, FASTA, cache slice and output collision
# is reported together and nothing is generated against a half-valid tree. Lines sharing
# their first five fields form a GROUP (one SNP/indel suite, or one assembly's 16-file SV
# set, in one format), which is the unit the output-volume gate below is applied to. A
# newline-delimited string rather than an array of records keeps this bash-3.2 safe.
JOBS=""
PREFLIGHT_FAILED=0
CHECKED_ASM=""

preflight_assembly() {
	local asm="$1"
	[[ " $CHECKED_ASM " == *" $asm "* ]] && return 0
	CHECKED_ASM="$CHECKED_ASM $asm"
	if [[ -z "$CORPUS_DIR" ]]; then
		assert_fasta "$asm" || PREFLIGHT_FAILED=1
	fi
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
CORPUS_STAGE=""
CORPUS_FASTA=""
if [[ -n "$CORPUS_DIR" ]]; then
	# The corpus's four outputs are staged in one scratch directory and moved into
	# expected/ together once every format has passed its gate, so a failing format
	# leaves the committed expected files untouched.
	if [[ "$DRY_RUN" == "true" ]]; then
		CORPUS_STAGE="$CORPUS_DIR/.generate.dry-run"
	else
		CORPUS_STAGE=$(mktemp -d "$CORPUS_DIR/.generate.XXXXXX")
		trap 'rm -rf "$CORPUS_STAGE"' EXIT
	fi
	n_selected=1
	preflight_assembly "$ASSEMBLY"
	if [[ -n "$CORPUS_FASTA_NAME" ]]; then
		fasta_plain="${CORPUS_FASTA_NAME%.gz}"
		if [[ ! -f "$CORPUS_DIR/$CORPUS_FASTA_NAME" || ! -f "$CORPUS_DIR/$fasta_plain.fai" ]]; then
			echo "ERROR: [generate_reference_output] corpus FASTA or its index not found: $CORPUS_DIR/$CORPUS_FASTA_NAME, $CORPUS_DIR/$fasta_plain.fai" >&2
			PREFLIGHT_FAILED=1
		elif [[ "$CORPUS_FASTA_NAME" == *.gz ]]; then
			# The golden tests decompress the FASTA once per run and read it with the
			# plain file's index; the container does the same.
			CORPUS_FASTA="$CORPUS_STAGE/$fasta_plain"
		else
			CORPUS_FASTA="$CORPUS_DIR/$CORPUS_FASTA_NAME"
		fi
	fi
	for fmt in $(formats_selected); do
		out_name=$(out_name_for corpus "" "$fmt")
		check_output_collision "$CORPUS_DIR/expected/$out_name.gz" || PREFLIGHT_FAILED=1
		check_output_collision "$CORPUS_DIR/expected/$out_name" || PREFLIGHT_FAILED=1
		JOBS="${JOBS}corpus|$ASSEMBLY|$FORK|$fmt|$CORPUS_STAGE|$CORPUS_DIR/variants.vcf|$out_name
"
	done
fi
while IFS= read -r entry; do
	[[ -n "$entry" && -z "$CORPUS_DIR" ]] || continue
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
		for fmt in $(formats_selected); do
			out_name=$(out_name_for snp_indel output "$fmt")
			check_output_collision "$SNP_INDEL_GT_DIR/$suite_id/$out_name" || PREFLIGHT_FAILED=1
			JOBS="${JOBS}snp_indel|$asm|$FORK|$fmt|$SNP_INDEL_GT_DIR/$suite_id|$in_vcf|$out_name
"
		done
	else
		in_dir="$SV_INPUTS_DIR/$asm_lc"
		if [[ ! -d "$in_dir" ]]; then
			echo "ERROR: [generate_reference_output] $suite_id ($suite_name): SV input directory not found: $in_dir" >&2
			PREFLIGHT_FAILED=1
			continue
		fi
		n_sv=0
		for fmt in $(formats_selected); do
			for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
				[[ -f "$vcf" ]] || continue
				b=$(basename "$vcf")
				b=${b%.vcf.gz}
				b=${b%.vcf}
				[[ "$b" == "manifest" ]] && continue
				n_sv=$((n_sv + 1))
				out_name=$(out_name_for sv "$b" "$fmt")
				check_output_collision "$SV_GT_DIR/$asm_lc/$out_name" || PREFLIGHT_FAILED=1
				JOBS="${JOBS}sv|$asm|$SV_FORK|$fmt|$SV_GT_DIR/$asm_lc|$vcf|$out_name
"
			done
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
bad_record=$(printf '%s' "$JOBS" | awk -F'|' 'NF != 7 && !seen { print; seen = 1 }')
[[ -z "$bad_record" ]] || die "a planned path holds '|' or a newline, which the plan's record format cannot carry: $bad_record"

# fasta_for_job <kind> <assembly>: the host FASTA a job mounts, empty for a corpus without one.
fasta_for_job() {
	if [[ "$1" == "corpus" ]]; then echo "$CORPUS_FASTA"; else fasta_for "$2"; fi
}

n_jobs=$(printf '%s' "$JOBS" | grep -c '|') || n_jobs=0
log "generate_reference_output image=$PERL_IMAGE cache_version=$PERL_CACHE_VERSION fork=$FORK sv_fork=$SV_FORK buffer_size=$BUFFER_SIZE format=$FORMAT stats=$KEEP_STATS volume_floor=$VOLUME_FLOOR"
[[ -z "$EXTRA_FLAGS" ]] || log "flags: $EXTRA_FLAGS"
if [[ -n "$CORPUS_DIR" ]]; then
	log "corpus=$CORPUS_DIR assembly=$ASSEMBLY planned_outputs=$n_jobs dry_run=$DRY_RUN"
else
	log "data_dir=$DATA_DIR snp_indel_gt=$SNP_INDEL_GT_DIR sv_gt=$SV_GT_DIR planned_outputs=$n_jobs dry_run=$DRY_RUN"
fi

# Dry run: print every command and stop.
if [[ "$DRY_RUN" == "true" ]]; then
	while IFS='|' read -r kind asm fork format out_dir in_vcf out_name; do
		[[ -n "$kind" ]] || continue
		if [[ "$kind" == "corpus" ]]; then
			echo "DRY-RUN [$kind $asm $format] -> $CORPUS_DIR/expected/$out_name.gz"
		else
			echo "DRY-RUN [$kind $asm $format] -> $out_dir/$out_name"
		fi
		docker_cmd "$asm" "$fork" "$format" "$in_vcf" "$out_dir" "$out_name" "$(fasta_for_job "$kind" "$asm")" "$(headers_for "$kind" "$format")"
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
# The VCF and JSON formats write one line per kept record whatever its consequences, so
# their floor is that every written record is annotated: a record VEP skipped before its
# input buffer is absent from both, while a record it could not annotate is present and
# bare. --volume-floor any lowers every format to one annotated unit per group, for a
# flag family that drops records or rows.
image_digests=$(docker image inspect --format '{{join .RepoDigests ","}}' "$PERL_IMAGE" 2>/dev/null || true)
if [[ -z "$image_digests" ]]; then
	log "pulling $PERL_IMAGE"
	docker pull "$PERL_IMAGE" >/dev/null
	image_digests=$(docker image inspect --format '{{join .RepoDigests ","}}' "$PERL_IMAGE" 2>/dev/null || true)
fi
[[ "$image_digests" == *"@$PERL_IMAGE_DIGEST"* ]] ||
	die "image digest mismatch for $PERL_IMAGE: ${image_digests:-none} (expected $PERL_IMAGE_DIGEST); pass --perl-image-digest for another image, nothing was generated"
IMAGE_DIGEST=$(printf '%s\n' "${image_digests//,/$'\n'}" | grep -F "@$PERL_IMAGE_DIGEST" | head -1)
log "image digest $IMAGE_DIGEST"

if [[ -n "$CORPUS_FASTA" && "$CORPUS_FASTA_NAME" == *.gz ]]; then
	gzip -dc "$CORPUS_DIR/$CORPUS_FASTA_NAME" >"$CORPUS_FASTA"
	cp "$CORPUS_DIR/${CORPUS_FASTA_NAME%.gz}.fai" "$CORPUS_FASTA.fai"
fi

done_jobs=0
GROUP_KEYS=$(printf '%s' "$JOBS" | cut -d'|' -f1-5 | uniq)
while IFS='|' read -r kind asm fork format out_dir; do
	[[ -n "$kind" ]] || continue
	group_prefix="$kind|$asm|$fork|$format|$out_dir|"
	fasta=$(fasta_for_job "$kind" "$asm")
	headers=$(headers_for "$kind" "$format")
	mkdir -p "$out_dir"
	scratch=$(mktemp -d "$out_dir/.generate.XXXXXX")
	group_records=0
	group_rows=0
	group_written=0
	group_files=0
	# Per-file facts for the provenance sidecars, recorded before the gate so a failure
	# discards them with the scratch directory. "in_vcf|records|units|lines|wall|out_name".
	group_facts=""
	log "group [$kind $asm $format --fork $fork] -> $out_dir/"
	while IFS='|' read -r _k _a _f _fmt _o in_vcf out_name; do
		[[ -n "$_k" ]] || continue
		in_records=$(count_records "$in_vcf")
		if [[ "$in_records" -eq 0 ]]; then
			rm -rf "$scratch"
			die "$out_dir/$out_name: input $in_vcf has 0 records; an empty reference would score as a clean cell"
		fi
		log "  $in_vcf ($in_records records) -> $out_name"
		start=$(date +%s)
		set +e
		eval "$(docker_cmd "$asm" "$fork" "$format" "$in_vcf" "$scratch" "$out_name" "$fasta" "$headers")"
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
		if [[ "$KEEP_STATS" == "true" && ! -f "$scratch/${out_name}_summary.txt" ]]; then
			rm -rf "$scratch"
			die "$out_dir/$out_name: --stats asked for the statistics report and vep wrote no ${out_name}_summary.txt; the group's partial output was discarded"
		fi
		out_lines=$(annotated_units "$format" "$scratch/$out_name")
		written=$(count_data_lines "$scratch/$out_name")
		[[ "$out_lines" -gt 0 ]] || echo "WARN: [generate_reference_output] $out_dir/$out_name has 0 annotated $format units for $in_records records" >&2
		group_records=$((group_records + in_records))
		group_rows=$((group_rows + out_lines))
		group_written=$((group_written + written))
		group_files=$((group_files + 1))
		group_facts="${group_facts}${in_vcf}|${in_records}|${out_lines}|${written}|${wall}|${out_name}
"
		log "    $out_lines annotated units in ${wall}s"
	# The prefix reaches awk through the environment: `-v` would read a backslash in a path
	# as an escape sequence and the group would match none of its own records.
	done <<<"$(printf '%s' "$JOBS" | PREFIX="$group_prefix" awk 'index($0, ENVIRON["PREFIX"]) == 1')"

	gate_failed=""
	if [[ "$VOLUME_FLOOR" == "any" ]]; then
		[[ "$group_rows" -ge 1 ]] || gate_failed="$out_dir ($format): no annotated unit for $group_records input records across $group_files file(s)"
	else
		case "$format" in
		default | tab)
			[[ "$group_rows" -ge "$group_records" ]] || gate_failed="$out_dir ($format): $group_rows annotation rows for $group_records input records across $group_files file(s)"
			;;
		vcf)
			[[ "$group_rows" -ge 1 && "$group_written" -eq "$group_rows" ]] || gate_failed="$out_dir ($format): $((group_written - group_rows)) of $group_written written records carry no $VCF_INFO_FIELD field ($group_records input records, $group_files file(s))"
			;;
		json)
			[[ "$group_rows" -ge 1 && "$group_written" -eq "$group_rows" ]] || gate_failed="$out_dir ($format): $((group_written - group_rows)) of $group_written written objects carry no consequences array ($group_records input records, $group_files file(s))"
			;;
		esac
	fi
	if [[ -n "$gate_failed" ]]; then
		echo "ERROR: [generate_reference_output] $gate_failed" >&2
		echo "  A VEP run that cannot reach its cache emits warnings rather than annotations and still exits 0." >&2
		echo "  Check --perl-cache-dir, --cache-version and --assembly; a flag family that drops records or rows" >&2
		echo "  (--chr, --coding_only, --no_intergenic, --filter_common, --allow_non_variant, --individual," >&2
		echo "  --individual_zyg, --dont_skip)" >&2
		echo "  or writes one row per record (--pick, --most_severe, --summary) on an input with records vep skips" >&2
		echo "  takes --volume-floor any. The group's output was discarded." >&2
		for w in "$scratch"/*_warnings.txt; do
			[[ -f "$w" ]] || continue
			head -3 "$w" | sed 's/^/  warnings: /' >&2
			break
		done
		rm -rf "$scratch"
		exit 1
	fi

	while IFS='|' read -r in_vcf in_records out_units out_lines wall out_name; do
		[[ -n "$out_name" ]] || continue
		stem=$(sidecar_stem "$kind" "$format" "$out_name")
		mv -f "$scratch/$out_name" "$out_dir/$out_name"
		if [[ -f "$scratch/${out_name}_warnings.txt" ]]; then
			mv -f "$scratch/${out_name}_warnings.txt" "$out_dir/$stem.warnings.log"
		fi
		if [[ -f "$scratch/${out_name}_summary.txt" ]]; then
			mv -f "$scratch/${out_name}_summary.txt" "$out_dir/$stem.summary.log"
		fi
		if [[ "$kind" != "corpus" ]]; then
			write_provenance "$out_dir/$stem.provenance.json" "$asm" "$fork" "$format" "$headers" "$in_vcf" "$in_records" \
				"$out_dir/$out_name" "$out_units" "$out_lines" "$wall" "$IMAGE_DIGEST"
		fi
		done_jobs=$((done_jobs + 1))
		log "  wrote $out_dir/$out_name"
	done <<<"$group_facts"
	rm -rf "$scratch"
	log "group [$kind $asm $format] complete: $group_rows annotated units / $group_records records over $group_files file(s)"
done <<<"$GROUP_KEYS"
[[ "$done_jobs" -eq "$n_jobs" ]] || die "wrote $done_jobs of $n_jobs planned output(s); a group matched none of its records"

# A corpus: gzip the four outputs and write provenance.json inside the scratch directory,
# then move everything into place, so expected/ changes only once all of it exists (the
# golden tests read `<name>` before `<name>.gz`, so a plain copy of an earlier generation
# goes). The keys the script owns are rewritten; `input` and `fasta` are filled only when
# absent, since a record may describe them in prose and the owned `vep_command` carries
# the file names; any other key of an existing record is kept.
if [[ -n "$CORPUS_DIR" ]]; then
	for fmt in $(formats_selected); do
		out_name=$(out_name_for corpus "" "$fmt")
		gzip -n -9 -c "$CORPUS_STAGE/$out_name" >"$CORPUS_STAGE/$out_name.gz"
	done
	fasta_text="none on either side"
	fasta_arg=""
	if [[ -n "$CORPUS_FASTA_NAME" ]]; then
		fasta_text="$CORPUS_FASTA_NAME"
		fasta_arg="${CORPUS_FASTA_NAME%.gz}"
	fi
	vep_command="vep -i variants.vcf -o <out> $(vep_args "$ASSEMBLY" "$FORK" "" '<cache>' "$fasta_arg" yes) [--tab|--vcf|--json]"
	python3 - "$CORPUS_DIR" "$CORPUS_STAGE" "$IMAGE_DIGEST" "$PERL_IMAGE" "$PERL_CACHE_VERSION" "$vep_command" \
		"homo_sapiens/${PERL_CACHE_VERSION}_${ASSEMBLY}" "$fasta_text" <<'EOF'
import hashlib, json, re, sys
from pathlib import Path

corpus, stage, image, tag, cache_version, command, perl_cache, fasta = sys.argv[1:9]
corpus, stage = Path(corpus), Path(stage)
path = corpus / "provenance.json"
record = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
if not isinstance(record, dict):
    sys.exit(f"{path}: the record is not a JSON object")
names = ("default.txt", "tab.txt", "vcf.vcf", "json.jsonl")
head = (stage / "default.txt").read_text(encoding="utf-8", errors="replace").splitlines()[:1]
version = re.match(r"## ENSEMBL VARIANT EFFECT PREDICTOR v(\S+)", head[0] if head else "")
# The version in the shape of vep's own version line, which reference_release.py parses.
reported = f"ensembl-vep          : {version.group(1)}" if version else None
owned = {
    "vep_image": image,
    "vep_image_tag": tag,
    "vep_version_reported": reported,
    "cache_version": int(cache_version),
    "vep_command": command,
    "perl_cache": perl_cache,
    "fasta": fasta,
    "input": "variants.vcf",
    "sha256": {n: hashlib.sha256((stage / n).read_bytes()).hexdigest() for n in names},
    "variants_vcf_sha256": hashlib.sha256((corpus / "variants.vcf").read_bytes()).hexdigest(),
}
filled_when_absent = ("fasta", "input")
for key, value in owned.items():
    if key in filled_when_absent:
        record.setdefault(key, value)
    else:
        record[key] = value
(stage / "provenance.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
EOF
	mkdir -p "$CORPUS_DIR/expected"
	for fmt in $(formats_selected); do
		out_name=$(out_name_for corpus "" "$fmt")
		rm -f "$CORPUS_DIR/expected/$out_name"
		mv -f "$CORPUS_STAGE/$out_name.gz" "$CORPUS_DIR/expected/$out_name.gz"
		for side in warnings summary; do
			[[ -f "$CORPUS_STAGE/$fmt.$side.log" ]] && mv -f "$CORPUS_STAGE/$fmt.$side.log" "$CORPUS_DIR/expected/$fmt.$side.log"
		done
		log "  wrote $CORPUS_DIR/expected/$out_name.gz"
	done
	mv -f "$CORPUS_STAGE/provenance.json" "$CORPUS_DIR/provenance.json"
	log "  wrote $CORPUS_DIR/provenance.json"
fi

log "generate_reference_output complete: $done_jobs/$n_jobs output(s) written"
