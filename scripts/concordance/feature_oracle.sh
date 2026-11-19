#!/usr/bin/env bash
#
# feature_oracle.sh -- the comparison oracle: the harness that runs Ensembl VEP and vep-rs
# on the same inputs and compares every column of every row in every output format. It runs
# the comparison cases of a case table: Ensembl VEP against one or more candidate `vep`
# binaries on the same input and flags, every requested output format compared column by
# column by compare_vep_formats.py, one report per candidate and case, one status table per
# run.
#
# THE CASE TABLE is a tab-separated file whose header names these columns, in any order.
# A cell with no value is `-`; blank lines and `#` lines are skipped.
#
#   case               unique name: the report directory and the status.tsv key
#   tier               fixture | population; --tier selects the rows of one tier
#   input              the input file: absolute, else relative to --data-dir, else to the
#                      repository root
#   cache              <version>_<assembly>, the Ensembl cache slice both engines annotate
#                      from: <species>/<version>_<assembly> under --perl-cache-dir
#                      (<species>_refseq or _merged when the flags say so) and
#                      <assembly, lower-cased> under --json-cache-dir; an optional
#                      =<dir> suffix names the candidate's JSON cache instead, and an
#                      =layout:<root> suffix a cache root in Ensembl VEP's layout
#                      (<root>/<species>[_refseq|_merged]/<version>_<assembly>), which
#                      the candidate resolves itself from --dir_cache <root> --species
#                      <species> --cache_version <version> --assembly <assembly> in
#                      place of --json_cache; <dir> and <root> resolve like input
#   species            the Ensembl species both engines run under (--species, the cache
#                      slice's directory: CacheDir.pm:306-315 at release/116.2 reads
#                      <dir_cache>/<species>[_refseq|_merged]/<version>_<assembly>), or -
#                      for homo_sapiens
#   fasta              - or an indexed FASTA given to both engines (its .fai beside it; a
#                      .gz is decompressed once per run, as the golden tests do, into a
#                      directory of its own under the work tree, so two sources that
#                      decompress to one file name never share a copy)
#   flags              whitespace-separated flags appended verbatim to both engines
#                      (--genomes to Ensembl VEP alone), or -;
#                      a flag is read in every spelling vep accepts (-nearest, --near,
#                      --NEAREST name --nearest); a --nearest case mounts the Ensembl
#                      cache writable, since Ensembl VEP writes transcript_gene_tss.txt
#                      there on its first run (a scan of every region file) and keeps the
#                      file, built under that run's transcript filters, for every later
#                      run; while a slice lacks the file, its --nearest runs take turns
#                      under one lock, because the temporary file is named after the
#                      process id, which is the same in every container; its comparisons
#                      settle NEAREST by start site through the candidate's JSON cache
#                      (compare_vep_formats.py --nearest-sites); a word may carry
#                      `{corpus}`, the input's directory, or `{sources}`, the
#                      --sources-dir directory: the candidate reads the host path and
#                      Ensembl VEP a read-only mount of it (/work/corpus, /work/sources),
#                      so a flag names a file beside the input or a shared source without
#                      an absolute path; a cell may name --dir, as `--dir <directory>` or
#                      `--dir=<directory>`, with the directory written as one of the two
#                      placeholders or a path under one, the placeholder standing in the
#                      first segment alone and no segment being `..`, since a flags cell
#                      reaches no other directory of Ensembl VEP's container; the mount is
#                      read-only, so a file a plugin keeps under --dir is there before
#                      the run (the GWAS plugin's processed catalogue, GWAS.pm:106-118,
#                      which a first run would write); both commands keep the runner's
#                      own cache words, which outrank --dir as the cache root
#                      (AnnotationSourceAdaptor.pm:127 at release/116.2), and a
#                      directory holding a vep.ini is refused, as --config is, since
#                      both engines read that file as options (Config.pm:704-705)
#   formats            comma-separated subset of default,tab,vcf,json,parquet,warnings;
#                      warnings compares the default run's `_warnings.txt` files and so
#                      needs default
#   perl_image         the Ensembl VEP image, pinned by digest (name@sha256:...)
#   perl_mode          offline | database
#   reference_release  the Ensembl VEP release the image carries, recorded per report
#   masks              compare_vep_formats.py --masks, applied from the manifest.json
#                      beside the input (whose flag set over the input is selected by the
#                      input's name, which every format's comparison passes), or -; the
#                      comparator's built-in classes
#                      (BUILTIN_MASKS below: regulatory_feature_rows, order_only,
#                      individual_allele_order, conservation_short_name_spelling,
#                      pc_description_text, motif_score_reverse_strand_mnv, the registry
#                      pair classes start_cooccurrence_swap and splice_family_swap, the
#                      field, row and warnings classes of the corrected behaviours, and the
#                      structural-variant comparator's classes) need no manifest. A class
#                      may carry its expected applied count, `<class>=<n>` (n at least 1):
#                      the case is red when an output format's report (not the warnings
#                      report, whose count is lines) applies the class to another number of
#                      entries, and a class with or without a count is red when every
#                      compared format, warnings included, leaves it unused, since a mask
#                      that sets nothing aside documents nothing; a class an output
#                      format's report names among those asked for and leaves out of its
#                      unused list is used though nothing is counted under it (the class
#                      a value set aside under another class needs named beside it).
#                      Beside a manifest.json the comparator is told two things of the
#                      candidate's run. Its run of the manifest (compare_vep_formats.py
#                      --run), whose lists alone are then read: the one run, the manifest's
#                      own (`default`) or a flag set, whose input has the case input's file
#                      name and whose flags are the flags cell's words in any order, the
#                      manifest's side read without --grch37 and without a --fasta word and
#                      its file, which the cache and fasta cells give; a case that is no run
#                      of its manifest, or fits several, names none, and every run over its
#                      input is read as one list. And whether the run has a reference
#                      sequence (--reference-sequence), as the candidate settles it: present
#                      under a fasta cell, which it is given as --fasta; else present where
#                      its cache directory holds a .fa or .fa.gz file with the .fai beside
#                      it, and the .gzi as well where the file's bytes are BGZF (a gzip file
#                      that is not BGZF is never taken), unless the flags name --no_fasta,
#                      which keeps the candidate from that file; else absent. A flags cell
#                      naming --fasta_dir says neither: the candidate takes that
#                      directory's file by the species and assembly in its name, and the
#                      cache directory's where none qualifies
#   floors             compare_vep_formats.py --floors (COL=MIN,...), or -; a floor on a
#                      column a format cannot count (the json format's record-level columns,
#                      Feature_type, Extra and a custom source's columns; a tab or vcf column
#                      neither header declares) is not applicable there: the report.json entry
#                      of that format lists it under floors_not_applicable, and a case whose
#                      formats give such a floor no format that counted it is red
#   field_pairs        compare_vep_formats.py --field-pairs (COL=EXPECTED=>ACTUAL;...), or -
#   unordered          - or 1: the comparator takes --unordered, so the shared rows agree in
#                      any order, for a run whose rows Ensembl VEP emits in hash order (the
#                      allele buckets of --pick_allele and --pick_allele_gene, the gene
#                      buckets of --per_gene)
#   reference          - or a directory of stored Ensembl VEP outputs for this case, resolved
#                      like input: {default.txt,tab.txt,vcf.vcf,json.jsonl}[.gz] for the
#                      case's formats (the layout of a golden corpus's expected directory),
#                      plus default.txt_warnings.txt[.gz] when the case names warnings (an
#                      empty file states no warnings); the case runs no container, a
#                      provenance.json beside the directory is recorded, and preflight
#                      refuses a directory lacking a format the case names
#   refusal            - or a directory, resolved like input, holding perl.txt and
#                      candidate.txt: the case expects both engines to refuse the run, and
#                      is green when, in every requested format, each engine exits non-zero
#                      with its file's text (without the trailing newline) in its stderr,
#                      red otherwise; nothing is compared, so the masks, floors,
#                      field_pairs, unordered and reference cells are -
#
# THE ENGINES. Ensembl VEP runs in Docker, once per case and output format, since a run
# writes one format (Runner.pm:871 at release/116.2 builds one OutputFactory per run):
#
#   vep -i <input> -o <out> --offline --cache --dir_cache <cache> --cache_version <v>
#       --species <species> --assembly <asm> --format vcf --fork 1 --no_stats --quiet
#       --force_overwrite <flags> [--fasta <fa>] [--tab|--vcf|--json]
#
# In database mode `--offline --cache --dir_cache --cache_version` become
# `--database --db_version <v>`: Config.pm:618 refuses --database beside --offline or
# --cache, and the HGVS, ID and SPDI input parsers (Parser/HGVS.pm:104, Parser/ID.pm:100,
# Parser/SPDI.pm:96) refuse to run offline. The server is Ensembl's public one, the
# Config.pm defaults (host ensembldb.ensembl.org, port 3306; Config.pm:711 moves GRCh37
# to 3337; a --genomes flag moves every species to Ensembl Genomes,
# mysql-eg-publicsql.ebi.ac.uk:4157, Config.pm:391-394), recorded in the report. The
# candidate runs with the same flags, less --genomes (it reads a local cache and refuses
# the flag by name), and
# `--offline --json_cache <cache>`; its parquet run is `--output_format parquet`,
# compared with the Ensembl VEP tab output. `--format vcf` is given to both engines
# unless the flags name an input format. The warnings compared are the default runs'
# `<out>_warnings.txt`, which Ensembl VEP opens on the first warning
# (BaseVEP.pm:878-897), so an absent file on either side means no warnings. A warning an
# Ensembl VEP worker raises while annotating is written with the worker's pid under
# --fork (Runner.pm:563 prefixes the captured worker stderr with it), so it differs from
# run to run; a record skipped while parsing is written by the parent without one.
#
# --reference-dir replaces the Ensembl VEP runs of every selected case with a stored
# reference: <reference-dir>/<case>/expected/, in the layout the reference column
# describes, with the provenance.json beside expected/ recorded. A case's own reference
# cell names its directory instead. No container runs for a case with a stored reference,
# and docker is needed only while a selected case has none; a refusal case runs both
# engines live under either.
#
# A CASE IS GREEN for a candidate when every requested format's comparison is green:
# no entry on one side only, every compared field equal after the masks and field
# pairs, every floor met, and for warnings the same (line, description) set. A format
# whose engine or comparison fails to run is red, and so is one whose two outputs both
# carry no annotation row: two empty outputs agree on nothing. A refusal case is green
# when both engines refuse every requested format with the expected texts. The run exits
# 1 on any red.
#
# OUTPUT
#   <out>/status.tsv                           case, candidate_pin, verdict,
#                                              report_digest (sha256 of report.json),
#                                              run_utc; one block per candidate
#   <out>/<candidate_pin>/<case>/report.json   the case record with its verdict
#   <out>/<candidate_pin>/<case>/<format>.json compare_vep_formats.py's report
#   <out>/run.json                             the run's options and image digests
#   <work-dir>/cases/<case>/                   the engines' outputs (perl/ and
#                                              candidates/<pin>/), deleted per case
#                                              unless --keep-outputs; fasta/ and logs/
#                                              beside cases/
#
# status.tsv and run.json describe the run that wrote them: both are removed when a run
# starts, so an interrupted run leaves no table, and a run over some cases of a table
# indexes those cases alone.
#
# Usage:
#   scripts/concordance/feature_oracle.sh --cases <table.tsv> --out <dir> \
#     --candidates <vep>[,<vep>...] --candidate-pins <pin>[,<pin>...] \
#     [--tier fixture|population] [--only <case>[,...]] [--data-dir <dir>] \
#     [--perl-cache-dir <dir>] [--json-cache-dir <dir>] [--reference-dir <dir>] \
#     [--sources-dir <dir>] [--work-dir <dir>] [--keep-outputs] [--fork <N>] \
#     [--perl-fork <N>] [--jobs <N>] [--format-jobs <N>] [--dry-run]
#
#   --format-jobs N runs a candidate's output formats of one case (its run and its
#   comparison) as N background jobs waited on, each logging to <fmt>.job.out and
#   <fmt>.job.err beside its outputs; the per-format verdicts are merged as under the
#   default, one format after another, so report.json reads the same either way.
#
#   --dry-run runs every preflight check and prints each engine and comparator command
#   instead of running them. A missing input fails a dry run as it fails a real one.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
COMPARATOR="$SCRIPT_DIR/compare_vep_formats.py"
PYTHON="${PYTHON:-python3}"

CASES=""
OUT_DIR=""
CANDIDATES_SPEC=""
PINS_SPEC=""
TIER="fixture"
ONLY_SPEC=""
DATA_DIR=""
PERL_CACHE_DIR=""
JSON_CACHE_DIR=""
REFERENCE_DIR=""
SOURCES_DIR=""
WORK_DIR=""
KEEP_OUTPUTS="false"
FORK=""
PERL_FORK=1
JOBS=1
FORMAT_JOBS=1
DRY_RUN="false"

usage() {
	cat <<'EOF'
Run a case table's comparison cases: Ensembl VEP against candidate vep binaries.

Required:
  --cases <tsv>             The case table (columns in the header of this script)
  --out <dir>               Where status.tsv and the per-case reports are written
  --candidates <bin>[,...]  Candidate vep binaries, each run on every selected case
  --candidate-pins <p>[,...] One pin per candidate (the commit or tag it was built
                            from): the report directory and status.tsv's candidate_pin

Optional:
  --tier <fixture|population>  The tier to run (default fixture)
  --only <case>[,...]    Restrict the run to these cases of the tier
  --data-dir <dir>          Root for relative inputs and the cache defaults
                            (default: the repository root)
  --perl-cache-dir <dir>    Directory holding <species>/<version>_<assembly>/
                            (default <data>/caches/perl/vep-cache)
  --json-cache-dir <dir>    Directory holding the candidate's <assembly>/ JSON caches
                            (default <data>/caches/vep-rs)
  --reference-dir <dir>     Stored Ensembl VEP output per case
                            (<dir>/<case>/expected/); a case's reference cell names
                            its own; no container runs for a case with one
  --sources-dir <dir>       The directory a {sources} flag word names: the candidate
                            reads it there, Ensembl VEP at a read-only mount
  --work-dir <dir>          Where the engines write (default <out>/work)
  --keep-outputs            Keep the engines' outputs after the comparison
  --fork <N>                The candidate's --fork (default: the binary's own)
  --perl-fork <N>           Ensembl VEP's --fork (default 1)
  --jobs <N>                Cases run at once (default 1)
  --format-jobs <N>         A candidate's output formats of one case run and compared
                            at once (default 1)
  --dry-run                 Preflight everything, print the commands, run nothing
EOF
	exit 1
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--cases)
		CASES="$2"
		shift 2
		;;
	--out)
		OUT_DIR="$2"
		shift 2
		;;
	--candidates)
		CANDIDATES_SPEC="$2"
		shift 2
		;;
	--candidate-pins)
		PINS_SPEC="$2"
		shift 2
		;;
	--tier)
		TIER="$2"
		shift 2
		;;
	--only)
		ONLY_SPEC="$2"
		shift 2
		;;
	--data-dir)
		DATA_DIR="$2"
		shift 2
		;;
	--perl-cache-dir)
		PERL_CACHE_DIR="$2"
		shift 2
		;;
	--json-cache-dir)
		JSON_CACHE_DIR="$2"
		shift 2
		;;
	--reference-dir)
		REFERENCE_DIR="$2"
		shift 2
		;;
	--sources-dir)
		SOURCES_DIR="$2"
		shift 2
		;;
	--work-dir)
		WORK_DIR="$2"
		shift 2
		;;
	--keep-outputs)
		KEEP_OUTPUTS="true"
		shift
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--perl-fork)
		PERL_FORK="$2"
		shift 2
		;;
	--jobs)
		JOBS="$2"
		shift 2
		;;
	--format-jobs)
		FORMAT_JOBS="$2"
		shift 2
		;;
	--dry-run)
		DRY_RUN="true"
		shift
		;;
	-h | --help) usage ;;
	*)
		echo "ERROR: [feature_oracle] unknown argument: $1" >&2
		usage
		;;
	esac
done

die() {
	echo "ERROR: [feature_oracle] $*" >&2
	exit 1
}
log() {
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
}
# `${var,,}` is bash 4+ and fails at runtime on bash 3.2 while passing `bash -n`.
lower() { printf %s "$1" | tr '[:upper:]' '[:lower:]'; }
abs_path() {
	local d
	d=$(cd "$(dirname "$1")" && pwd) || return 1
	printf '%s/%s' "$d" "$(basename "$1")"
}
# norm_path <path>: absolute, without `.`, `..`, repeated or trailing slashes; the path
# need not exist (symbolic links are not resolved).
norm_path() {
	local p="$1" out="" seg
	[[ "$p" == /* ]] || p="$PWD/$p"
	# The unquoted $p splits on IFS alone: set -f keeps a `*` or `?` in a segment from
	# matching the working directory's entries.
	local IFS='/'
	set -f
	for seg in $p; do
		case "$seg" in
		'' | '.') ;;
		'..') out="${out%/*}" ;;
		*) out="$out/$seg" ;;
		esac
	done
	set +f
	printf '%s' "${out:-/}"
}
sha256_of() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	elif command -v shasum >/dev/null 2>&1; then
		shasum -a 256 "$1" | awk '{print $1}'
	else
		echo "unavailable"
	fi
}
# Whether the whitespace-separated list $1 contains the word $2.
has_word() { [[ " $1 " == *" $2 "* ]]; }
# in_lines <list> <item>: whether the newline-separated <list> has a line equal to <item>.
in_lines() { printf '%s' "$1" | grep -Fx -- "$2" >/dev/null; }
# tsv_lookup <table> <key>: the second column of the first line of the two-column
# tab-separated <table> whose first column is <key>, or nothing.
tsv_lookup() {
	printf '%s' "$1" | awk -F'\t' -v k="$2" '$1 == k {print $2; exit}'
}

# Validate the options
[[ -n "$CASES" ]] || {
	echo "ERROR: [feature_oracle] --cases required" >&2
	usage
}
[[ -n "$OUT_DIR" ]] || {
	echo "ERROR: [feature_oracle] --out required" >&2
	usage
}
[[ -n "$CANDIDATES_SPEC" ]] || {
	echo "ERROR: [feature_oracle] --candidates required" >&2
	usage
}
[[ -n "$PINS_SPEC" ]] || {
	echo "ERROR: [feature_oracle] --candidate-pins required (one pin per candidate)" >&2
	usage
}
[[ -f "$CASES" ]] || die "case table not found: $CASES"
case "$TIER" in
fixture | population) ;;
*) die "invalid --tier: $TIER (fixture|population)" ;;
esac
for _n in "$PERL_FORK" "$JOBS" "$FORMAT_JOBS" ${FORK:+"$FORK"}; do
	[[ "$_n" =~ ^[0-9]+$ && "$_n" -gt 0 ]] || die "--fork, --perl-fork, --jobs and --format-jobs take a positive integer (got '$_n')"
done
[[ -n "$DATA_DIR" ]] || DATA_DIR="$REPO_DIR"
[[ -d "$DATA_DIR" ]] || die "data-dir not a directory: $DATA_DIR"
DATA_DIR=$(abs_path "$DATA_DIR")
[[ -n "$PERL_CACHE_DIR" ]] || PERL_CACHE_DIR="$DATA_DIR/caches/perl/vep-cache"
[[ -n "$JSON_CACHE_DIR" ]] || JSON_CACHE_DIR="$DATA_DIR/caches/vep-rs"
# Docker binds the Ensembl cache root by host path, which must be absolute; the report
# records both roots.
[[ ! -d "$PERL_CACHE_DIR" ]] || PERL_CACHE_DIR=$(abs_path "$PERL_CACHE_DIR")
[[ ! -d "$JSON_CACHE_DIR" ]] || JSON_CACHE_DIR=$(abs_path "$JSON_CACHE_DIR")
if [[ -n "$REFERENCE_DIR" ]]; then
	[[ -d "$REFERENCE_DIR" ]] || die "reference-dir not a directory: $REFERENCE_DIR"
	REFERENCE_DIR=$(abs_path "$REFERENCE_DIR")
fi
if [[ -n "$SOURCES_DIR" ]]; then
	[[ -d "$SOURCES_DIR" ]] || die "sources-dir not a directory: $SOURCES_DIR"
	SOURCES_DIR=$(abs_path "$SOURCES_DIR")
fi
# The flag-word placeholders and the container paths they take in the Ensembl VEP
# command; the candidate takes the host directories.
CORPUS_PLACEHOLDER='{corpus}'
SOURCES_PLACEHOLDER='{sources}'
CORPUS_MOUNT="/work/corpus"
SOURCES_MOUNT="/work/sources"
[[ -f "$COMPARATOR" ]] || die "comparator not found beside this script: $COMPARATOR"
command -v "$PYTHON" >/dev/null 2>&1 || die "$PYTHON is required (the comparator and the report writer run on it)"

# Candidates and pins: one pin per binary, pins unique and usable as directory names.
CANDIDATES=()
PINS=()
IFS=',' read -r -a CANDIDATES <<<"$CANDIDATES_SPEC"
IFS=',' read -r -a PINS <<<"$PINS_SPEC"
[[ "${#CANDIDATES[@]}" -eq "${#PINS[@]}" ]] || die "--candidates names ${#CANDIDATES[@]} binaries and --candidate-pins ${#PINS[@]} pins; they go one to one"
PIN_RE='^[A-Za-z0-9][A-Za-z0-9._-]*$'
[[ -n "$WORK_DIR" ]] || WORK_DIR="$OUT_DIR/work"
# The names the run itself writes under --out; a pin is a directory there.
RESERVED_PINS="status.tsv run.json"
NORM_WORK_DIR=$(norm_path "$WORK_DIR")
if [[ "${NORM_WORK_DIR%/*}" == "$(norm_path "$OUT_DIR")" ]]; then
	RESERVED_PINS="$RESERVED_PINS ${NORM_WORK_DIR##*/}"
fi
for ((i = 0; i < ${#CANDIDATES[@]}; i++)); do
	bin="${CANDIDATES[$i]}"
	pin="${PINS[$i]}"
	[[ -n "$bin" && -f "$bin" && -x "$bin" ]] || die "candidate is not an executable file: '$bin'"
	CANDIDATES[i]=$(abs_path "$bin")
	[[ "$pin" =~ $PIN_RE ]] || die "candidate pin '$pin' is not a directory-safe name ([A-Za-z0-9._-], not starting with . or -)"
	has_word "$RESERVED_PINS" "$pin" && die "candidate pin '$pin' names a file or directory the run writes under --out ($RESERVED_PINS)"
	for ((j = 0; j < i; j++)); do
		[[ "${PINS[$j]}" != "$pin" ]] || die "candidate pin '$pin' is given twice; every candidate needs its own pin"
	done
done

# The case table
COLUMNS="case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs species unordered reference refusal"
N_TABLE_COLUMNS=17
DEFAULT_SPECIES="homo_sapiens"
OUTPUT_FORMATS="default tab vcf json parquet"
# The cache cell's suffix that names a cache root in Ensembl VEP's layout.
LAYOUT_PREFIX="layout:"
# The comparator's own divergence classes, applied without a manifest: the row and text
# classes, the registry pair classes, the field classes of the corrected behaviours, the
# transcript row classes of the merged and RefSeq caches, the structural-variant comparator's
# classes and the warnings-mode classes (compare_vep_formats.py REGISTRY_CLASSES,
# FIELD_CLASSES, MERGED_SAME_STABLE_ID_COLLAPSE, TERM_PAIR_CLASSES, STRUCTURAL_VARIANT_CLASSES,
# PLUGIN_WARNING_CLASSES and WARNING_TEXT_CLASSES).
BUILTIN_MASKS="regulatory_feature_rows order_only individual_allele_order conservation_short_name_spelling pc_description_text motif_score_reverse_strand_mnv colocated_list_cut_at_cache_cap start_cooccurrence_swap splice_family_swap gxf_reverse_strand_utr_sequence organelle_peptide_paths_translation_table conservation_max_text_order gwas_haplotype_cross_product gwas_interaction_first_snp_allele_whitespace alphamissense_two_digit_transcript_version cadd_structural_allele_match splicevault_sibling_allele_rows splicevault_json_later_rows utrannotator_annotation_pair_order utrannotator_gained_start_window utrannotator_consequence_precedence utrannotator_structural_allele_text utrannotator_lost_start_after_insertion utrannotator_gained_start_remade_in_place utrannotator_frameshift_at_stop_boundary utrannotator_lost_stop_written_before_stop_codon nearest_exon_jb_item_order loeuf_zero_minimum maxentscan_overlapping_windows alphamissense_sibling_allele_rows revel_sibling_allele_rows revel_other_reference_amino_acid_rows gnomadmt_info_value_second_equals loftee_anc_allele_multiallelic_record loftee_first_alt_donor_filter loftee_gerp_zoom_level_summary nearest_exon_jb_stretch_priming nearest_exon_jb_shifted_span nearest_exon_jb_id_less_exon_fold downstream_reverse_strand_rotation spliceai_chopped_multiallelic spliceai_multi_entry_cell plugin_contig_mitochondrial_spelling refseq_edited_model_hgvsc_offset refseq_edited_model_shifted_insertion_bases refseq_edited_model_seam_insertion_hgvsc refseq_edited_model_unminimised_allele_hgvsc refseq_edited_model_hgvs_offset refseq_edited_model_deleted_span_term refseq_edited_model_tail_term mature_mirna_unshifted_span stop_codon_unshifted_span stop_retained_without_stop_codon_test stop_codon_overlap_genomic_span stop_retained_three_base_stop_codon_test stop_codon_removed_without_stop_codon_test stop_codon_overlap_insertion_range splice_term_past_the_record_span refseq_edited_model_deleted_span_codons hgvs_intron_offset_abutting_insertion hgvs_shifted_insertion_clamped_end shift_intron_term_unshifted_span hgvs_shifted_indel_written_span_gate hgvsp_shifted_indel_written_span_gate hgvs_shifted_deletion_clamped_span hgvs_clipped_reference_against_whole_alternate hgvs_added_bases_left_at_the_span_end hgvs_insertion_left_where_the_shift_stops hgvs_added_bases_written_outside_the_transcript hgvs_shift_sequence_edge_window hgvs_shift_shared_across_alleles hgvsp_stop_codon_protein_comparison hgvsp_stop_named_on_incomplete_coding_sequence hgvsp_intron_end_insertion hgvsp_intron_start_insertion frequency_filter_zero_value merged_same_stable_id_collapse hgnc_id_copy downstream_utr_frame_past_coding_end downstream_start_past_own_codon plugin_hgvs_shifted_cds_coordinates plugin_value_follows_splice_term_past_the_record_span plugin_value_follows_stop_codon_term_pair plugin_value_on_reference_skipped_record blosum62_hgvs_shift_shared_across_alleles reference_skipped_record breakend_mate_local_read bnd_own_chromosome_rows_lost dosage_sensitivity_regulatory_row_warnings fork_warning_pid_prefix"
# The built-in classes that read beside the two outputs: the skipped-record class, the
# classes that read a record's alleles and the two classes of splice_term_past_the_record_span,
# which read the reference bases of a row under an identifier, the input VCF in every format
# (the comparator reads it for the VCF and JSON formats in any case), the breakend classes,
# hgvsp_intron_end_insertion and hgvsp_intron_start_insertion the candidate's JSON cache
# (compare_vep_formats.py --json-cache).
# organelle_peptide_paths_translation_table reads its transcripts' codon tables from that cache.
# The two Downstream classes read their transcripts' coding sequences, 3' UTRs and codon tables from it.
# plugin_hgvs_shifted_cds_coordinates reads the transcripts' strands, exons, coding bounds and mapper
# pairs from it as well, and the record's alleles for a later allele's row.
# hgvs_shift_shared_across_alleles decides three columns of a later allele's row from the record's
# alleles and the cache's mapper pairs, coding sequences and codon tables; its other rows are the
# manifest's.
# hgvs_clipped_reference_against_whole_alternate counts a transcript's exon spans from it.
# maxentscan_overlapping_windows reads a transcript's introns from it and sets nothing aside without it.
INPUT_READING_MASKS="reference_skipped_record splicevault_sibling_allele_rows alphamissense_sibling_allele_rows revel_sibling_allele_rows revel_other_reference_amino_acid_rows loftee_anc_allele_multiallelic_record loftee_first_alt_donor_filter spliceai_chopped_multiallelic utrannotator_gained_start_window hgvs_shift_shared_across_alleles plugin_hgvs_shifted_cds_coordinates splice_term_past_the_record_span plugin_value_follows_splice_term_past_the_record_span blosum62_hgvs_shift_shared_across_alleles"
CACHE_READING_MASKS="breakend_mate_local_read bnd_own_chromosome_rows_lost organelle_peptide_paths_translation_table downstream_utr_frame_past_coding_end downstream_start_past_own_codon plugin_hgvs_shifted_cds_coordinates hgvs_shift_shared_across_alleles hgvs_clipped_reference_against_whole_alternate hgvsp_intron_end_insertion hgvsp_intron_start_insertion maxentscan_overlapping_windows"
# The built-in class the comparator decides from the FASTA the case gives both engines and
# the candidate's JSON cache (compare_vep_formats.py --reference-fasta with --json-cache); a
# case that names no FASTA leaves it to the manifest.json beside its input.
REFERENCE_READING_MASKS="hgvs_shift_sequence_edge_window"
# The built-in classes that read the GWAS catalog file the case's flags hand the GWAS plugin
# (compare_vep_formats.py --gwas-catalog) and set nothing aside without it.
CATALOG_READING_MASKS="gwas_haplotype_cross_product gwas_interaction_first_snp_allele_whitespace"
# Flags the runner sets from the table's own columns and its options, with the aliases
# Ensembl VEP's Config.pm and the candidate accept for them; a flags cell naming one would
# give an engine the flag twice, or move the input, the server or the assembly behind the
# report's back. A cell names one in any spelling an engine reads as it: a single dash or a
# plus sign, any letter case, another name of the same option (-u, --full_cache_dir) and an
# abbreviation, each resolved against the two option lists below.
# --dir is not among them: the runner sets none, a command's
# own cache words outrank it (AnnotationSourceAdaptor.pm:127 at release/116.2 reads
# `dir_cache || dir`), and preflight holds its directory to a read-only placeholder mount
# without a vep.ini, the file Config.pm:704-705 reads as options from under --dir.
RUNNER_FLAGS="-i --input_file --input_data --id -o --output_file --offline --cache --database --db_version --host --port --user --password --pass --dir_cache --cache_version --json_cache -s --species -a --assembly --grch37 --fasta --fa --tab --vcf --json --output_format --force_overwrite --force --no_stats -q --quiet --fork --threads --warning_file --skipped_variants_file --config"
# The options each engine declares, one word per option, its names joined by `|`.
# PERL_OPTIONS is Ensembl VEP's `@VEP_PARAMS` (Config.pm:73-276 at release/116.2) in the
# source's order, which VEP_PARAMS of crates/vep-cli/src/argv.rs carries; CANDIDATE_OPTIONS
# is the candidate's parser, sorted: the flag rows of the reference pages under docs/cli/
# that are generated from it, and the flags its error table catches before parsing.
# test_feature_oracle.sh holds each list to that source.
PERL_OPTIONS="help config input_file|i input_data|id format output_format delimiter no_check_variants_order species|s registry host port user|u password|pass db_version assembly|a grch37 genomes refseq merged all_refseq gencode_basic gencode_primary is_multispecies transcript_filter exclude_predicted minimal most_severe summary pick pick_allele per_gene pick_allele_gene flag_pick flag_pick_allele flag_pick_allele_gene flag_gencode_primary pick_order buffer_size failed gp chr check_ref lookup_ref check_existing check_svs no_check_alleles exclude_null_alleles check_frequency af af_1kg af_gnomade af_gnomad af_gnomadg old_maf max_af pubmed freq_filter freq_freq freq_gt_lt freq_pop filter_common allow_non_variant process_ref_homs individual individual_zyg phased fork dont_skip nearest distance clin_sig_allele clinvar_somatic_classification overlaps max_sv_size remove_hgvsp_version verbose|v quiet|q no_progress everything|e output_file|o compress_output no_headers mask_header_cache_path stats_file|sf stats_text stats_html no_stats warning_file skipped_variants_file force_overwrite|force terms|t coding_only canonical mane mane_select tsl appris ccds xref_refseq uniprot protein protein_version biotype hgnc symbol transcript_version gene_version gene_phenotype mirna spdi ga4gh_vrs hgvs hgvsg hgvsg_use_accession hgvsp_use_prediction shift_hgvs ambiguous_hgvs sift polyphen humdiv condel variant_class regulatory cell_type regulatory_gff extended_promoters convert no_intergenic vcf solr json tab vcf_info_field keep_csq keep_ann lrg fields domains numbers total_length allele_number show_ref_allele no_escape ambiguity var_synonyms shift_3prime shift_genomic shift_length uploaded_allele database cache cache_version show_cache_info dir dir_cache full_cache_dir dir_plugins offline fasta|fa fasta_dir no_fasta sereal synonyms bam use_transcript_ref use_given_ref custom tmpdir gff gtf bigwig phyloP phastCons ucsc_assembly ucsc_data_root custom_multi_allelic custom_suppress_filter plugin safe debug"
CANDIDATE_OPTIONS="af af_1kg af_exac af_gnomad af_gnomade af_gnomadg all_refseq allele_number allow_non_variant ambiguity ambiguous_hgvs appris assembly|a bam bigwig biotype buffer_size cache cache_version canonical ccds cell_type check_existing check_frequency check_ref check_svs chr clin_sig_allele clinvar_somatic_classification coding_only compress_output compression_threads condel config convert custom custom_multi_allelic custom_suppress_filter database db_version debug decompression_threads delimiter dir dir_cache dir_plugins distance domains dont_skip everything|e exclude_null_alleles exclude_predicted extended_promoters failed fasta_dir fasta|fa fields filter_common flag_gencode_primary flag_pick flag_pick_allele flag_pick_allele_gene force_overwrite|force fork|threads format freq_filter freq_freq freq_gt_lt freq_pop ga4gh_vrs gencode_basic gencode_primary gene_phenotype gene_version genomes gff gp grch37 gtf hgnc hgvs hgvsg hgvsg_use_accession hgvsp_use_prediction host humdiv individual individual_zyg input_data|id input_file|i is_multispecies json json_cache|full_cache_dir keep_ann keep_csq lookup_ref lrg mane mane_select mask_header_cache_path max_af max_sv_size merged minimal mirna most_severe nearest no_check_alleles no_check_variants_order no_escape no_fasta no_headers no_intergenic no_progress no_stats numbers offline old_maf output_buffer_size output_file|o output_format overlaps parquet_partition parse_threads password|pass per_gene phased phastCons phyloP pick pick_allele pick_allele_gene pick_order pipeline_depth pipeline_mode plugin polyphen port prewarm_ahead process_ref_homs protein protein_version pubmed quiet|q reader_core refseq registry regulatory regulatory_gff remove_hgvsp_version render_task_records resolved_config safe sereal shift_3prime shift_genomic shift_hgvs shift_length show_cache_info show_ref_allele sift skipped_variants_file solr spdi species|s stats_file|sf stats_html stats_text summary symbol synonyms tab terms|t tmpdir total_length transcript_filter transcript_version tsl ucsc_assembly ucsc_data_root uniprot uploaded_allele use_given_ref use_transcript_ref user|u var_synonyms variant_class vcf vcf_info_field verbose|v warning_file writer_strategy xref_refseq"
# The word of RUNNER_FLAGS that the option name in OPTION_NAME reaches in either list, or
# nothing. An engine reads a name as the option that has it among its names, else as the
# one option whose names begin with it (Getopt::Long's auto_abbrev, which Config.pm:62
# leaves on; resolve_spelling of argv.rs reads the candidate's options, then Ensembl
# VEP's), and a name that several options begin with reaches none.
# The run of a manifest.json a case is: the name of the one run, the manifest's own
# (`default`) or a flag set, whose input has the file name of argument 2 and whose flags are
# the words after it in any order, the manifest's side read without --grch37 and without a
# --fasta word and its file; nothing where no run or more than one fits, or the manifest
# cannot be read. The names are compare_vep_formats.py's OWN_RUN and CORPUS_INPUT.
CASE_RUN_PROGRAM='
import json, sys
try:
    with open(sys.argv[1], encoding="utf-8") as handle:
        manifest = json.load(handle)
except (OSError, ValueError):
    sys.exit(0)
def words(flags):
    kept, value = [], False
    for word in flags if isinstance(flags, list) else []:
        if value or word == "--grch37" or str(word).startswith("--fasta="):
            value = False
        elif word == "--fasta":
            value = True
        else:
            kept.append(str(word))
    return sorted(kept)
sets = manifest.get("flag_sets") if isinstance(manifest, dict) else None
runs = [("default", manifest.get("flags"), "variants.vcf")] if isinstance(manifest, dict) else []
runs += [(s.get("name"), s.get("flags"), s.get("input", "variants.vcf")) for s in sets or [] if isinstance(s, dict)]
fits = [name for name, flags, run_input in runs if run_input == sys.argv[2] and words(flags) == sorted(sys.argv[3:])]
if len(fits) == 1 and isinstance(fits[0], str):
    print(fits[0])
'
RUNNER_WORD_PROGRAM='
function option_of(name, list,    options, names, n, m, i, j, hit, hits) {
	n = split(list, options, " ")
	hit = hits = 0
	for (i = 1; i <= n; i++) {
		m = split(tolower(options[i]), names, "|")
		for (j = 1; j <= m; j++) {
			if (names[j] == name) return i
			if (index(names[j], name) == 1 && hit != i) {
				hit = i
				hits++
			}
		}
	}
	return hits == 1 ? hit : 0
}
function runner_word(name, list,    options, names, m, i, j, pass) {
	i = option_of(name, list)
	if (!i) return ""
	split(list, options, " ")
	m = split(options[i], names, "|")
	for (pass = 1; pass <= 2; pass++) {
		for (j = 1; j <= m; j++) {
			if (pass == 1 && index(tolower(names[j]), name) != 1) continue
			if (index(runner, " --" names[j] " ")) return "--" names[j]
			if (index(runner, " -" names[j] " ")) return "-" names[j]
		}
	}
	return ""
}
BEGIN {
	name = ENVIRON["OPTION_NAME"]
	word = runner_word(name, perl)
	if (word == "") word = runner_word(name, candidate)
	printf "%s", word
}'

# A run's table describes that run: the previous one goes before anything is checked, so
# a run that fails preflight or is interrupted leaves no table behind.
[[ "$DRY_RUN" == "true" ]] || rm -f "$OUT_DIR/status.tsv" "$OUT_DIR/run.json"

# Header: every column present once; the column order is read from it.
HEADER=""
while IFS= read -r line; do
	[[ -n "$line" && "$line" != \#* ]] || continue
	HEADER="$line"
	break
done <"$CASES"
[[ -n "$HEADER" ]] || die "$CASES: no header line"
N_COLUMNS=$(printf '%s\n' "$HEADER" | awk -F'\t' '{print NF}')
# COL_<name>=<1-based index>
for name in $COLUMNS; do
	idx=$(printf '%s\n' "$HEADER" | awk -F'\t' -v n="$name" '{for (i = 1; i <= NF; i++) if ($i == n) {print i; exit}}')
	[[ -n "$idx" ]] || die "$CASES: header lacks the column '$name' (needs: $COLUMNS)"
	eval "COL_$name=$idx"
done
[[ "$N_COLUMNS" -eq "$N_TABLE_COLUMNS" ]] || die "$CASES: header has $N_COLUMNS columns; the $N_TABLE_COLUMNS are: $COLUMNS"

# cell <line> <column name>: the cell's text.
cell() {
	local idx
	eval "idx=\$COL_$2"
	printf '%s\n' "$1" | awk -F'\t' -v c="$idx" '{print $c}'
}

# Row selection: the tier, then --only.
ONLY_NAMES=""
if [[ -n "$ONLY_SPEC" ]]; then
	ONLY_NAMES=$(printf '%s' "$ONLY_SPEC" | tr ',' ' ')
fi
SELECTED=""
ALL_NAMES=""
TIER_NAMES=""
line_no=0
seen_header="false"
while IFS= read -r line; do
	line_no=$((line_no + 1))
	[[ -n "$line" && "$line" != \#* ]] || continue
	if [[ "$seen_header" != "true" ]]; then
		seen_header="true"
		continue
	fi
	nf=$(printf '%s\n' "$line" | awk -F'\t' '{print NF}')
	[[ "$nf" -eq "$N_TABLE_COLUMNS" ]] || die "$CASES:$line_no: $nf columns, expected $N_TABLE_COLUMNS (an empty cell is written -)"
	printf '%s\n' "$line" | awk -F'\t' '{for (i = 1; i <= NF; i++) if ($i == "") exit 1}' ||
		die "$CASES:$line_no: an empty cell; write - for no value"
	cid=$(cell "$line" case)
	[[ "$cid" =~ $PIN_RE ]] || die "$CASES:$line_no: case '$cid' is not a directory-safe name ([A-Za-z0-9._-], not starting with . or -)"
	has_word "$ALL_NAMES" "$cid" && die "$CASES:$line_no: case '$cid' appears twice"
	ALL_NAMES="$ALL_NAMES $cid"
	tier=$(cell "$line" tier)
	case "$tier" in
	fixture | population) ;;
	*) die "$CASES:$line_no: tier '$tier' is not fixture or population" ;;
	esac
	[[ "$tier" == "$TIER" ]] || continue
	TIER_NAMES="$TIER_NAMES $cid"
	if [[ -n "$ONLY_NAMES" ]] && ! has_word "$ONLY_NAMES" "$cid"; then continue; fi
	SELECTED="${SELECTED}${line}
"
done <"$CASES"
for cid in $ONLY_NAMES; do
	has_word "$ALL_NAMES" "$cid" || die "--only names '$cid', which is not a case of $CASES"
	has_word "$TIER_NAMES" "$cid" || die "--only names '$cid', which is not a case of tier $TIER"
done
N_SELECTED=$(printf '%s' "$SELECTED" | grep -c . 2>/dev/null) || N_SELECTED=0
[[ "$N_SELECTED" -gt 0 ]] || die "no case of tier $TIER selected from $CASES${ONLY_SPEC:+ by --only $ONLY_SPEC}"

# One case: fields and resolved paths

# resolve_path <cell>: absolute as given, else under --data-dir, else under the repository.
resolve_path() {
	local p="$1"
	case "$p" in
	/*) printf '%s' "$p" ;;
	*)
		if [[ -e "$DATA_DIR/$p" ]]; then
			printf '%s/%s' "$DATA_DIR" "$p"
		else
			printf '%s/%s' "$REPO_DIR" "$p"
		fi
		;;
	esac
}

# option_name <word>: the option named by a vep argument word, in its lower-case `--name`
# spelling with a value suffix removed, or nothing for a value word. vep reads its options
# with Getopt::Long's defaults (vep:27, Config.pm:62 at release/116.2), which take a single
# dash or a plus sign, any unambiguous abbreviation and any letter case as the full option.
# An abbreviation is returned as written: preflight resolves it against the option lists.
option_name() {
	local name="${1%%=*}"
	case "$name" in
	--[A-Za-z]*) name=$(lower "$name") ;;
	-[A-Za-z]* | +[A-Za-z]*) name=$(lower "--${name#[-+]}") ;;
	*) return 0 ;;
	esac
	printf '%s' "$name"
}
# has_option <--name>: whether a flag word of the loaded case names the option <--name> in
# a spelling vep reads as it, that is, whether the word's option_name is a prefix of <--name>.
has_option() {
	local word name
	for word in ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}; do
		name=$(option_name "$word")
		[[ -n "$name" && "$1" == "$name"* ]] && return 0
	done
	return 1
}
# gwas_catalog_file: the file the loaded case's flags hand the GWAS plugin as the candidate
# reads it (`--plugin GWAS,file=<path>`, the value its own word or joined by `=`), or nothing
# where no flag word names one.
gwas_catalog_file() {
	local word name value parameter plugin_value="false"
	local parameters=()
	for word in ${C_CAND_FLAG_WORDS[@]+"${C_CAND_FLAG_WORDS[@]}"}; do
		name=$(option_name "$word")
		value=""
		if [[ -z "$name" ]]; then
			[[ "$plugin_value" != "true" ]] || value="$word"
			plugin_value="false"
		elif [[ "--plugin" == "$name"* && "$word" == *=* ]]; then
			value="${word#*=}"
			plugin_value="false"
		elif [[ "--plugin" == "$name"* ]]; then
			plugin_value="true"
		else
			plugin_value="false"
		fi
		[[ "$value" == GWAS,* ]] || continue
		IFS=',' read -r -a parameters <<<"$value"
		for parameter in "${parameters[@]}"; do
			if [[ "$parameter" == file=* ]]; then
				printf '%s' "${parameter#file=}"
				return 0
			fi
		done
	done
	return 0
}
# tss_first_run: whether the loaded case's first Ensembl VEP run writes
# transcript_gene_tss.txt: an offline --nearest case whose cache slice lacks the file.
tss_first_run() {
	[[ "$C_PERL_MODE" == "offline" && ! -e "$C_PERL_SLICE/transcript_gene_tss.txt" ]] && has_option --nearest
}

# load_case <line>: the C_* globals of one row.
load_case() {
	local line="$1" cache_spec species_dir w engine rewritten head corpus sources
	C_ID=$(cell "$line" case)
	C_TIER=$(cell "$line" tier)
	C_INPUT=$(cell "$line" input)
	C_FASTA=$(cell "$line" fasta)
	C_FLAGS=$(cell "$line" flags)
	C_FORMATS=$(cell "$line" formats)
	C_PERL_IMAGE=$(cell "$line" perl_image)
	C_PERL_MODE=$(cell "$line" perl_mode)
	C_REFERENCE_RELEASE=$(cell "$line" reference_release)
	C_MASKS=$(cell "$line" masks)
	# The cell's classes without their expected counts: the comparator's --masks value
	# (write_report reads the counts from the cell itself).
	C_MASK_CLASSES="$C_MASKS"
	[[ "$C_MASKS" == "-" ]] || C_MASK_CLASSES=$(printf '%s' "$C_MASKS" | sed 's/=[^,]*//g')
	C_FLOORS=$(cell "$line" floors)
	C_FIELD_PAIRS=$(cell "$line" field_pairs)
	C_SPECIES=$(cell "$line" species)
	[[ "$C_SPECIES" != "-" ]] || C_SPECIES="$DEFAULT_SPECIES"
	C_UNORDERED=$(cell "$line" unordered)
	C_REFERENCE=$(cell "$line" reference)
	C_REFUSAL=$(cell "$line" refusal)
	cache_spec=$(cell "$line" cache)
	[[ "$C_FLAGS" != "-" ]] || C_FLAGS=""
	# The flags as words: `read -a` splits on blanks without expanding a glob, so a `*` in
	# a flag value reaches the engines as written.
	C_FLAG_WORDS=()
	[[ -z "$C_FLAGS" ]] || read -r -a C_FLAG_WORDS <<<"$C_FLAGS"

	C_INPUT_PATH=$(resolve_path "$C_INPUT")
	C_CORPUS_DIR=$(dirname "$C_INPUT_PATH")
	C_MANIFEST="$C_CORPUS_DIR/manifest.json"
	[[ -f "$C_MANIFEST" ]] || C_MANIFEST=""
	# The flag words as each engine reads them: the two placeholders replaced by the
	# container mounts for Ensembl VEP and by the host directories for the candidate. Not
	# `${w//"$placeholder"/"$dir"}`: bash 3.2 keeps the quotes of a quoted replacement and
	# bash 5.2 reads `&` in an unquoted one as the matched text, so each occurrence is cut
	# out by prefix and suffix removal, literal on both.
	C_PERL_FLAG_WORDS=()
	C_CAND_FLAG_WORDS=()
	for w in ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}; do
		for engine in perl candidate; do
			if [[ "$engine" == perl ]]; then
				corpus="$CORPUS_MOUNT" sources="$SOURCES_MOUNT"
			else
				corpus="$C_CORPUS_DIR" sources="$SOURCES_DIR"
			fi
			rewritten="$w"
			head=""
			while [[ "$rewritten" == *"$CORPUS_PLACEHOLDER"* ]]; do
				head="$head${rewritten%%"$CORPUS_PLACEHOLDER"*}$corpus"
				rewritten="${rewritten#*"$CORPUS_PLACEHOLDER"}"
			done
			rewritten="$head$rewritten"
			head=""
			while [[ "$rewritten" == *"$SOURCES_PLACEHOLDER"* ]]; do
				head="$head${rewritten%%"$SOURCES_PLACEHOLDER"*}$sources"
				rewritten="${rewritten#*"$SOURCES_PLACEHOLDER"}"
			done
			if [[ "$engine" == perl ]]; then
				C_PERL_FLAG_WORDS+=("$head$rewritten")
			else
				C_CAND_FLAG_WORDS+=("$head$rewritten")
			fi
		done
	done
	C_FASTA_PATH=""
	[[ "$C_FASTA" == "-" ]] || C_FASTA_PATH=$(resolve_path "$C_FASTA")
	C_REFUSAL_DIR=""
	[[ "$C_REFUSAL" == "-" ]] || C_REFUSAL_DIR=$(resolve_path "$C_REFUSAL")
	# The directory the case's Ensembl VEP outputs are read from instead of a run: its own
	# reference cell, else its directory under --reference-dir; a refusal case runs live.
	C_REFERENCE_DIR=""
	if [[ -n "$C_REFUSAL_DIR" ]]; then
		:
	elif [[ "$C_REFERENCE" != "-" ]]; then
		C_REFERENCE_DIR=$(resolve_path "$C_REFERENCE")
	elif [[ -n "$REFERENCE_DIR" ]]; then
		C_REFERENCE_DIR="$REFERENCE_DIR/$C_ID/expected"
	fi

	C_CACHE="${cache_spec%%=*}"
	C_CACHE_VERSION="${C_CACHE%%_*}"
	C_ASSEMBLY="${C_CACHE#*_}"
	# CacheDir.pm:312 at release/116.2 appends _refseq, then _merged, for each flag set;
	# Config.pm:621 refuses the two together, so one suffix follows the species.
	species_dir="$C_SPECIES"
	if has_option --merged; then
		species_dir="${C_SPECIES}_merged"
	elif has_option --refseq; then
		species_dir="${C_SPECIES}_refseq"
	fi
	C_PERL_SLICE="$PERL_CACHE_DIR/$species_dir/$C_CACHE"
	# A layout root is what the candidate is given; the leaf it resolves is checked and recorded.
	C_DIR_CACHE=""
	if [[ "$cache_spec" == *="$LAYOUT_PREFIX"* ]]; then
		C_DIR_CACHE=$(resolve_path "${cache_spec#*="$LAYOUT_PREFIX"}")
		C_JSON_CACHE="$C_DIR_CACHE/$species_dir/$C_CACHE"
	elif [[ "$cache_spec" == *=* ]]; then
		C_JSON_CACHE=$(resolve_path "${cache_spec#*=}")
	else
		C_JSON_CACHE="$JSON_CACHE_DIR/$(lower "$C_ASSEMBLY")"
	fi

	C_FORMAT_LIST=$(printf '%s' "$C_FORMATS" | tr ',' ' ')
	C_OUTPUT_FORMATS=""
	C_WARNINGS="false"
	C_PERL_FORMATS=""
	local f
	for f in $C_FORMAT_LIST; do
		if [[ "$f" == "warnings" ]]; then
			C_WARNINGS="true"
		else
			C_OUTPUT_FORMATS="$C_OUTPUT_FORMATS $f"
		fi
	done
	C_OUTPUT_FORMATS="${C_OUTPUT_FORMATS# }"
	for f in $C_OUTPUT_FORMATS; do
		f=$(perl_format "$f")
		has_word "$C_PERL_FORMATS" "$f" || C_PERL_FORMATS="$C_PERL_FORMATS $f"
	done
	C_PERL_FORMATS="${C_PERL_FORMATS# }"
	# The server a database-mode run reaches: Config.pm:287-289 at release/116.2 default
	# host ensembldb.ensembl.org and port 3306; Config.pm:711-713 set port 3337 for GRCh37;
	# Config.pm:391-394 set the Ensembl Genomes server under --genomes, whatever the species.
	C_PERL_SERVER=""
	if [[ "$C_PERL_MODE" == "database" ]]; then
		if has_option --genomes; then
			C_PERL_SERVER="mysql-eg-publicsql.ebi.ac.uk:4157"
		elif [[ "$(lower "$C_ASSEMBLY")" == "grch37" ]]; then
			C_PERL_SERVER="ensembldb.ensembl.org:3337"
		else
			C_PERL_SERVER="ensembldb.ensembl.org:3306"
		fi
	fi
}

# The output file of a format inside an engine's output directory; parquet is a directory.
format_file() {
	case "$1" in
	default) echo "default.txt" ;;
	tab) echo "tab.txt" ;;
	vcf) echo "vcf.vcf" ;;
	json) echo "json.jsonl" ;;
	parquet) echo "parquet.parquet" ;;
	esac
}

# perl_format <format>: the Ensembl VEP output a candidate format is compared with; the
# parquet comparison reads the tab output.
perl_format() {
	if [[ "$1" == "parquet" ]]; then
		echo "tab"
	else
		echo "$1"
	fi
}

# plain_or_gz <path>: <path> if it is a file, else <path>.gz if that is, else nothing.
plain_or_gz() {
	local p
	for p in "$1" "$1.gz"; do
		if [[ -f "$p" ]]; then
			printf '%s' "$p"
			return 0
		fi
	done
	return 1
}
# stored_reference_file <format>: the loaded case's stored output, plain or .gz, or nothing.
stored_reference_file() {
	plain_or_gz "$C_REFERENCE_DIR/$(format_file "$1")"
}
# stored_warnings_file: the loaded case's stored default-run warnings file, plain or .gz, or nothing.
stored_warnings_file() {
	plain_or_gz "$C_REFERENCE_DIR/default.txt_warnings.txt"
}
# case_work_dir: the loaded case's directory under the work tree. Cases live under cases/,
# candidates under <case>/candidates/, so no case or pin name can reach the Ensembl VEP
# outputs (<case>/perl/), the decompressed FASTAs (fasta/) or the logs (logs/).
case_work_dir() { printf '%s/cases/%s' "$WORK_DIR" "$C_ID"; }

# expected_output <perl format>: the Ensembl VEP output a comparison of the loaded case
# reads: the stored reference, or the file its run writes under the work tree.
expected_output() {
	if [[ -n "$C_REFERENCE_DIR" ]]; then
		stored_reference_file "$1"
	else
		printf '%s/perl/%s' "$(case_work_dir)" "$(format_file "$1")"
	fi
}

# Preflight: every selected case is checked before anything runs
PREFLIGHT_FAILED=0
pf_error() {
	echo "ERROR: [feature_oracle] $C_ID: $*" >&2
	PREFLIGHT_FAILED=1
}
# FASTA_SOURCES: every FASTA path, one per line, in first-use order; a .gz source is
# decompressed into fasta/<its line number>/ under the work tree. TSS_SLICES: the cache
# slice of every offline --nearest case whose first Ensembl VEP run writes
# transcript_gene_tss.txt, one line per case. IMAGES: the image of every case that runs
# Ensembl VEP live, and N_LIVE their count.
FASTA_SOURCES=""
TSS_SLICES=""
IMAGES=""
N_LIVE=0
# An assembly name may carry a dot or a dash (BDGP6.54, R64-1-1).
CACHE_RE='^[0-9]+_[A-Za-z0-9][A-Za-z0-9._-]*$'
# An Ensembl species production name, the cache directory's name.
SPECIES_RE='^[a-z][a-z0-9_]*$'
FORMAT_SPEC_RE='^[a-z,]+$'

while IFS= read -r line; do
	[[ -n "$line" ]] || continue
	load_case "$line"
	[[ -f "$C_INPUT_PATH" ]] || pf_error "input not found: $C_INPUT_PATH (from '$C_INPUT')"
	if [[ -n "$C_FASTA_PATH" ]]; then
		[[ -f "$C_FASTA_PATH" ]] || pf_error "FASTA not found: $C_FASTA_PATH"
		fai="${C_FASTA_PATH%.gz}.fai"
		[[ -f "$fai" ]] || pf_error "FASTA index not found: $fai"
		in_lines "$FASTA_SOURCES" "$C_FASTA_PATH" || FASTA_SOURCES="${FASTA_SOURCES}${C_FASTA_PATH}
"
	fi
	[[ "$C_CACHE" =~ $CACHE_RE ]] || pf_error "cache '$C_CACHE' is not <version>_<assembly>"
	[[ "$C_SPECIES" =~ $SPECIES_RE ]] || pf_error "species '$C_SPECIES' is not an Ensembl species name (lower-case letters, digits and underscores)"
	# CacheDir.pm:310-311 at release/116.2 refuses a suffixed species: the flag selects the slice.
	case "$C_SPECIES" in
	*_refseq | *_merged) pf_error "species '$C_SPECIES' carries the slice suffix; write the species and the --${C_SPECIES##*_} flag" ;;
	esac
	if [[ -n "$C_DIR_CACHE" && ! -d "$C_JSON_CACHE/transcripts" ]]; then
		pf_error "cache root $C_DIR_CACHE has no $(basename "$(dirname "$C_JSON_CACHE")")/$C_CACHE/transcripts/ (the layout the candidate resolves from --dir_cache)"
	elif [[ ! -d "$C_JSON_CACHE/transcripts" ]]; then
		pf_error "candidate JSON cache has no transcripts/: $C_JSON_CACHE"
	fi
	case "$C_PERL_MODE" in
	offline)
		if [[ -z "$C_REFERENCE_DIR" && ! -d "$C_PERL_SLICE" ]]; then
			pf_error "Ensembl VEP cache slice not found: $C_PERL_SLICE"
		elif [[ -z "$C_REFERENCE_DIR" ]] && tss_first_run; then
			TSS_SLICES="${TSS_SLICES}${C_PERL_SLICE}
"
		fi
		;;
	database) ;;
	*) pf_error "perl_mode '$C_PERL_MODE' is not offline or database" ;;
	esac
	[[ "$C_PERL_IMAGE" != "-" ]] || pf_error "perl_image is -"
	if [[ -z "$C_REFERENCE_DIR" ]]; then
		N_LIVE=$((N_LIVE + 1))
		has_word "$IMAGES" "$C_PERL_IMAGE" || IMAGES="$IMAGES $C_PERL_IMAGE"
	fi
	[[ "$C_REFERENCE_RELEASE" != "-" ]] || pf_error "reference_release is -"
	[[ "$C_UNORDERED" == "-" || "$C_UNORDERED" == "1" ]] || pf_error "unordered '$C_UNORDERED' is neither - nor 1"
	if [[ -n "$C_REFUSAL_DIR" ]]; then
		[[ "$C_REFERENCE" == "-" ]] || pf_error "refusal and reference are both set; a refused run has no output to store"
		[[ "$C_WARNINGS" != "true" ]] || pf_error "refusal with warnings in formats; a refused run writes no warnings file"
		for column_cell in "masks:$C_MASKS" "floors:$C_FLOORS" "field_pairs:$C_FIELD_PAIRS" "unordered:$C_UNORDERED"; do
			[[ "${column_cell#*:}" == "-" ]] || pf_error "refusal with ${column_cell%%:*} '${column_cell#*:}'; a refused run compares nothing, so the cell is -"
		done
		if [[ ! -d "$C_REFUSAL_DIR" ]]; then
			pf_error "refusal directory not found: $C_REFUSAL_DIR (from '$C_REFUSAL')"
		else
			for f in perl.txt candidate.txt; do
				if [[ ! -f "$C_REFUSAL_DIR/$f" ]]; then
					pf_error "refusal directory lacks $f: $C_REFUSAL_DIR"
				elif [[ -z "$(tr -d '[:space:]' <"$C_REFUSAL_DIR/$f")" ]]; then
					pf_error "refusal text $C_REFUSAL_DIR/$f is empty, which every stderr would match"
				fi
			done
		fi
	fi
	if [[ ! "$C_FORMATS" =~ $FORMAT_SPEC_RE ]] || [[ -z "$C_OUTPUT_FORMATS" ]]; then
		pf_error "formats '$C_FORMATS' names no output format (default,tab,vcf,json,parquet, and warnings)"
	fi
	seen=""
	for f in $C_FORMAT_LIST; do
		if ! has_word "$OUTPUT_FORMATS warnings" "$f"; then
			pf_error "formats names '$f'; the formats are default,tab,vcf,json,parquet and warnings"
		elif has_word "$seen" "$f"; then
			pf_error "formats names '$f' twice"
		fi
		seen="$seen $f"
	done
	if [[ "$C_WARNINGS" == "true" ]] && ! has_word "$C_OUTPUT_FORMATS" "default"; then
		pf_error "formats names warnings without default; the warnings compared are the default run's"
	fi
	for ((k = 0; k < ${#C_FLAG_WORDS[@]}; k++)); do
		f="${C_FLAG_WORDS[$k]}"
		# The option the word names, in any spelling vep reads; the list spells a short flag
		# with one dash.
		flag_name=$(option_name "$f")
		if [[ -n "$flag_name" ]] && { has_word "$RUNNER_FLAGS" "$flag_name" || has_word "$RUNNER_FLAGS" "${flag_name#-}"; }; then
			pf_error "flags names $f, which the runner sets from the table's columns and its own options"
		elif [[ -n "$flag_name" ]]; then
			# A word the list does not spell is refused when an engine resolves it to a listed
			# option, as another of the option's names or as an abbreviation.
			runner_word=$(OPTION_NAME="${flag_name#--}" awk -v runner=" $RUNNER_FLAGS " -v perl="$PERL_OPTIONS" -v candidate="$CANDIDATE_OPTIONS" "$RUNNER_WORD_PROGRAM")
			[[ -z "$runner_word" ]] || pf_error "flags names $f, which vep reads as $runner_word, a word the runner sets from the table's columns and its own options"
		fi
		# The value of --dir is the rest of the word after `=`, else the next word whatever it
		# begins with (`dir=s`, Config.pm:241 at release/116.2); the candidate's word at the
		# same index carries the host directory.
		[[ "$flag_name" == "--dir" ]] || continue
		if [[ "$f" == *=* ]]; then
			dir_value="${f#*=}"
			dir_host="${C_CAND_FLAG_WORDS[$k]#*=}"
		else
			dir_value="${C_FLAG_WORDS[$((k + 1))]-}"
			dir_host="${C_CAND_FLAG_WORDS[$((k + 1))]-}"
		fi
		dir_root="${dir_value%%/*}"
		dir_rest="${dir_value#"$dir_root"}"
		if [[ ("$dir_root" != "$CORPUS_PLACEHOLDER" && "$dir_root" != "$SOURCES_PLACEHOLDER") || "/$dir_value/" == */../* || "$dir_rest" == *"$CORPUS_PLACEHOLDER"* || "$dir_rest" == *"$SOURCES_PLACEHOLDER"* ]]; then
			pf_error "flags give --dir '$dir_value', which is not $CORPUS_PLACEHOLDER, $SOURCES_PLACEHOLDER or a path under one without a .. segment or a second placeholder: a flags cell reaches no other directory of Ensembl VEP's container"
		# A {sources} directory has a host path only under --sources-dir, whose absence is
		# reported below.
		elif [[ "$dir_root" == "$CORPUS_PLACEHOLDER" || -n "$SOURCES_DIR" ]] && [[ -e "$dir_host/vep.ini" || -L "$dir_host/vep.ini" ]]; then
			pf_error "flags give --dir '$dir_value', which holds a vep.ini ($dir_host/vep.ini): both engines read that file as options (Config.pm:704-705 at release/116.2), so it is refused as --config is"
		fi
	done
	if [[ "$C_FLAGS" == *"$SOURCES_PLACEHOLDER"* && -z "$SOURCES_DIR" ]]; then
		pf_error "flags carry $SOURCES_PLACEHOLDER and no --sources-dir names the directory it stands for"
	fi
	if [[ -n "$C_REFERENCE_DIR" && ! -d "$C_REFERENCE_DIR" ]]; then
		pf_error "stored reference directory not found: $C_REFERENCE_DIR"
	elif [[ -n "$C_REFERENCE_DIR" ]]; then
		for f in $C_PERL_FORMATS; do
			stored_reference_file "$f" >/dev/null || pf_error "stored reference lacks $(format_file "$f")[.gz] under $C_REFERENCE_DIR/"
		done
		if [[ "$C_WARNINGS" == "true" ]]; then
			stored_warnings_file >/dev/null || pf_error "stored reference lacks default.txt_warnings.txt[.gz] under $C_REFERENCE_DIR/ (an empty file states no warnings)"
		fi
	fi
	# In a comparator cell `-` is the whole cell or nothing: no mask class, floor or field pair
	# begins with a dash, and a cell that does is read by the comparator's option parser as
	# another option.
	for column_cell in "masks:$C_MASKS" "floors:$C_FLOORS" "field_pairs:$C_FIELD_PAIRS"; do
		column="${column_cell%%:*}"
		value="${column_cell#*:}"
		items=$(printf '%s' "$value" | tr ',;' '\n\n')
		if [[ "$value" != "-" && "$items" =~ (^|$'\n')- ]]; then
			pf_error "$column '$value' has an item beginning with -; write - alone for no value"
		fi
	done
	if [[ "$C_MASKS" != "-" && ( "$C_MASKS" == ,* || "$C_MASKS" == *, || "$C_MASKS" == *,,* || "$C_MASKS" == *[[:space:]]* ) ]]; then
		pf_error "masks '$C_MASKS' has an empty or blank item"
	fi
	# `read -a` splits the cell at its commas without expanding a glob, so a `*` item is
	# refused as written instead of being read as the names of files in the working directory.
	IFS=',' read -r -a mask_items <<<"$C_MASKS"
	seen_mask_classes=""
	for f in "${mask_items[@]}"; do
		if [[ "$f" != "-" && ! "$f" =~ ^[a-z0-9_]+(=[1-9][0-9]*)?$ ]]; then
			pf_error "masks item '$f' is not <class> or <class>=<expected applied count>"
		fi
		f="${f%%=*}"
		# An item with no class name is refused above; has_word tests the padded list for two
		# spaces when the word is empty, so both checks below would report a class with no name.
		[[ -n "$f" ]] || continue
		if has_word "$seen_mask_classes" "$f"; then
			pf_error "masks names '$f' twice"
		fi
		seen_mask_classes="$seen_mask_classes $f"
		if [[ "$f" != "-" && -z "$C_MANIFEST" ]] && ! has_word "$BUILTIN_MASKS" "$f"; then
			pf_error "masks names '$f', which applies from a manifest.json beside the input, and none is there (the built-in classes are $BUILTIN_MASKS)"
		fi
		if [[ -z "$C_MANIFEST" && -z "$C_FASTA_PATH" ]] && has_word "$REFERENCE_READING_MASKS" "$f"; then
			pf_error "masks names '$f', which is decided from the case's FASTA or applies from a manifest.json beside the input, and the case has neither"
		fi
	done
done <<<"$SELECTED"
[[ "$PREFLIGHT_FAILED" -eq 0 ]] || die "preflight failed; nothing ran"
if [[ "$DRY_RUN" != "true" && "$N_LIVE" -gt 0 ]]; then
	command -v docker >/dev/null 2>&1 || die "docker is required to run Ensembl VEP on the $N_LIVE selected case(s) without a stored reference (pass --reference-dir, or --dry-run)"
fi
# Without flock, run_perl cannot make the --nearest cases of one slice take turns, and two
# first runs at once would tear transcript_gene_tss.txt; one such case per slice needs no turn.
if [[ "$JOBS" -gt 1 && "$DRY_RUN" != "true" ]] && ! command -v flock >/dev/null 2>&1; then
	slice=$(printf '%s' "$TSS_SLICES" | awk 'seen[$0]++ {print; exit}')
	[[ -z "$slice" ]] || die "flock is required to run the --nearest cases of $slice under --jobs $JOBS while the slice lacks transcript_gene_tss.txt (install flock, or pass --jobs 1)"
fi

# Commands

# perl_cmd <format> <host out dir> <host fasta or empty>: CMD becomes the docker command.
perl_cmd() {
	local fmt="$1" out_dir="$2" fasta="$3"
	local in_name
	in_name=$(basename "$C_INPUT_PATH")
	# Parser/VCF.pm:639 at release/116.2 orders an individual's alternate alleles, and
	# OutputFactory.pm:983 the ZYG list of --individual_zyg, by Perl hash order; the committed
	# --individual and --individual_zyg references are runs under this seed.
	CMD=(docker run --rm --user "$(id -u):$(id -g)" -e PERL_HASH_SEED=0 -e PERL_PERTURB_KEYS=0
		-v "${C_INPUT_PATH}:/work/input/${in_name}:ro"
		-v "${out_dir}:/work/out")
	if [[ "$C_PERL_MODE" == "offline" ]]; then
		# Ensembl VEP's first --nearest run writes transcript_gene_tss.txt into the cache
		# directory and dies on a read-only one (AnnotationSource/Cache/Transcript.pm:288-289).
		if has_option --nearest; then
			CMD+=(-v "${PERL_CACHE_DIR}:/work/cache")
		else
			CMD+=(-v "${PERL_CACHE_DIR}:/work/cache:ro")
		fi
	fi
	if [[ -n "$fasta" ]]; then
		CMD+=(-v "${fasta}:/work/ref.fa:ro" -v "${fasta}.fai:/work/ref.fa.fai:ro")
	fi
	if [[ "$C_FLAGS" == *"$CORPUS_PLACEHOLDER"* ]]; then
		CMD+=(-v "${C_CORPUS_DIR}:${CORPUS_MOUNT}:ro")
	fi
	if [[ "$C_FLAGS" == *"$SOURCES_PLACEHOLDER"* ]]; then
		CMD+=(-v "${SOURCES_DIR}:${SOURCES_MOUNT}:ro")
	fi
	CMD+=("$C_PERL_IMAGE" vep -i "/work/input/${in_name}" -o "/work/out/$(format_file "$fmt")")
	if [[ "$C_PERL_MODE" == "offline" ]]; then
		CMD+=(--offline --cache --dir_cache /work/cache --cache_version "$C_CACHE_VERSION")
	else
		CMD+=(--database --db_version "$C_CACHE_VERSION")
	fi
	CMD+=(--species "$C_SPECIES" --assembly "$C_ASSEMBLY")
	has_option --format || CMD+=(--format vcf)
	CMD+=(--fork "$PERL_FORK" --no_stats --quiet --force_overwrite)
	CMD+=(${C_PERL_FLAG_WORDS[@]+"${C_PERL_FLAG_WORDS[@]}"})
	[[ -z "$fasta" ]] || CMD+=(--fasta /work/ref.fa)
	case "$fmt" in
	default) ;;
	*) CMD+=("--$fmt") ;;
	esac
}

# candidate_cmd <binary> <format> <host out dir> <host fasta or empty>
candidate_cmd() {
	local bin="$1" fmt="$2" out_dir="$3" fasta="$4"
	CMD=("$bin" -i "$C_INPUT_PATH" -o "$out_dir/$(format_file "$fmt")" --offline)
	if [[ -n "$C_DIR_CACHE" ]]; then
		CMD+=(--dir_cache "$C_DIR_CACHE" --species "$C_SPECIES" --cache_version "$C_CACHE_VERSION" --assembly "$C_ASSEMBLY")
	else
		CMD+=(--json_cache "$C_JSON_CACHE" --species "$C_SPECIES" --assembly "$C_ASSEMBLY")
	fi
	has_option --format || CMD+=(--format vcf)
	[[ -z "$FORK" ]] || CMD+=(--fork "$FORK")
	CMD+=(--no_stats --quiet --force_overwrite)
	# --genomes selects the Ensembl VEP run's server, as --database and --db_version do; the
	# candidate reads a local cache and refuses the flag by name, so a word naming it in any
	# spelling vep reads is not given it.
	local w name
	for w in ${C_CAND_FLAG_WORDS[@]+"${C_CAND_FLAG_WORDS[@]}"}; do
		name=$(option_name "$w")
		[[ -n "$name" && "--genomes" == "$name"* ]] || CMD+=("$w")
	done
	[[ -z "$fasta" ]] || CMD+=(--fasta "$fasta")
	case "$fmt" in
	default) ;;
	parquet) CMD+=(--output_format parquet) ;;
	*) CMD+=("--$fmt") ;;
	esac
}

# compare_cmd <format> <expected> <actual> <report>
compare_cmd() {
	local fmt="$1" mask_class mask_reads_input="false" mask_reads_cache="false" mask_reads_catalog="false" catalog
	local mask_classes=()
	local mask_reads_reference="false" run reference_sequence fa head
	IFS=',' read -r -a mask_classes <<<"$C_MASK_CLASSES"
	for mask_class in "${mask_classes[@]}"; do
		has_word "$REFERENCE_READING_MASKS" "$mask_class" && mask_reads_reference="true"
		has_word "$INPUT_READING_MASKS" "$mask_class" && mask_reads_input="true"
		has_word "$CACHE_READING_MASKS" "$mask_class" && mask_reads_cache="true"
		has_word "$CATALOG_READING_MASKS" "$mask_class" && mask_reads_catalog="true"
	done
	CMD=("$PYTHON" "$COMPARATOR" --format "$fmt" --expected "$2" --actual "$3" --report "$4")
	# --input where the comparator reads the records (the vcf and json formats, an
	# input-reading class) and beside a manifest, which documents each flag set under its
	# input's name and is read whole when the input is unnamed; a text-format comparison
	# reads none of the input's lines for the manifest.
	if [[ "$fmt" == vcf || "$fmt" == json || "$mask_reads_input" == "true" || -n "$C_MANIFEST" ]]; then
		CMD+=(--input "$C_INPUT_PATH")
	fi
	[[ -z "$C_MANIFEST" ]] || CMD+=(--manifest "$C_MANIFEST")
	# Beside a manifest the comparator is told the candidate's run (the header's masks column
	# states both rules): its run of the manifest, and whether it has a reference sequence, as
	# vep settles it (fasta_source of crates/vep-cli/src/runner.rs, detect_fasta of
	# cache_resolve.rs, first_indexed_fasta of crates/vep-fasta/src/lib.rs).
	if [[ -n "$C_MANIFEST" ]]; then
		run=$("$PYTHON" -c "$CASE_RUN_PROGRAM" "$C_MANIFEST" "$(basename "$C_INPUT_PATH")" ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}) || run=""
		[[ -z "$run" ]] || CMD+=(--run "$run")
		reference_sequence=""
		if [[ -n "$C_FASTA_PATH" ]]; then
			reference_sequence="present"
		elif ! has_option --fasta_dir; then
			reference_sequence="absent"
			if ! has_option --no_fasta; then
				for fa in "$C_JSON_CACHE"/*.fa "$C_JSON_CACHE"/*.fa.gz; do
					[[ -f "$fa" && -f "$fa.fai" ]] || continue
					# The first 18 bytes, as hex pairs between single blanks, tell a BGZF file
					# (htslib check_header, bgzf.c: the gzip magic and deflate method, FLG.FEXTRA,
					# XLEN 6 and the BC subfield of SLEN 2) from any other gzip file and from
					# plain text.
					head=$(head -c 18 "$fa" | od -An -v -tx1 | tr -s ' \n' ' ')
					head=" ${head# }"
					head="${head% } "
					case "$head" in
					" 1f 8b 08 "?[4567cdef]" "??" "??" "??" "??" "??" "??" 06 00 42 43 02 00 "??" "??" ") [[ -f "$fa.gzi" ]] || continue ;;
					" 1f 8b "*) continue ;;
					esac
					reference_sequence="present"
					break
				done
			fi
		fi
		[[ -z "$reference_sequence" ]] || CMD+=(--reference-sequence "$reference_sequence")
	fi
	[[ "$C_MASKS" == "-" ]] || CMD+=(--masks "$C_MASK_CLASSES")
	[[ "$C_FLOORS" == "-" ]] || CMD+=(--floors "$C_FLOORS")
	[[ "$C_FIELD_PAIRS" == "-" ]] || CMD+=(--field-pairs "$C_FIELD_PAIRS")
	[[ "$C_UNORDERED" != "1" ]] || CMD+=(--unordered)
	# NEAREST compares by transcription start site: Ensembl VEP names one transcript at a
	# winning site and the candidate every transcript there (compare_vep_formats.py).
	if has_option --nearest; then
		CMD+=(--nearest-sites "$C_JSON_CACHE")
	fi
	[[ "$mask_reads_cache" != "true" ]] || CMD+=(--json-cache "$C_JSON_CACHE")
	if [[ "$mask_reads_catalog" == "true" ]]; then
		catalog=$(gwas_catalog_file)
		[[ -z "$catalog" ]] || CMD+=(--gwas-catalog "$catalog")
	fi
	# A reference-reading class is decided from the FASTA the engines read and the
	# candidate's cache; without a FASTA the comparator applies it from the manifest.
	if [[ "$mask_reads_reference" == "true" && -n "$C_FASTA_PATH" ]]; then
		[[ "$mask_reads_cache" == "true" ]] || CMD+=(--json-cache "$C_JSON_CACHE")
		CMD+=(--reference-fasta "$(case_fasta)")
	fi
}
warnings_cmd() {
	CMD=("$PYTHON" "$COMPARATOR" --warnings --expected "$1" --actual "$2" --report "$3")
	[[ "$C_MASKS" == "-" ]] || CMD+=(--masks "$C_MASK_CLASSES")
}

print_cmd() {
	printf '%q ' "${CMD[@]}"
	echo
}

# Run setup
RUN_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
N_CANDIDATES="${#CANDIDATES[@]}"
if [[ "$TIER" == "population" && "$N_LIVE" -gt 0 ]]; then
	echo "WARN: [feature_oracle] population tier: Ensembl VEP runs live on $N_LIVE case(s) without a stored reference" >&2
fi
log "feature_oracle tier=$TIER cases=$N_SELECTED live=$N_LIVE candidates=$N_CANDIDATES table=$CASES dry_run=$DRY_RUN"
log "out=$OUT_DIR work=$WORK_DIR reference_dir=${REFERENCE_DIR:-none} perl_cache=$PERL_CACHE_DIR json_cache=$JSON_CACHE_DIR"

# decompressed_fasta <.gz source>: where the run decompresses the source: a directory of
# its own, numbered by the source's line in FASTA_SOURCES, so two sources of one file name
# never share a copy.
decompressed_fasta() {
	local n
	n=$(printf '%s' "$FASTA_SOURCES" | grep -nFx -- "$1" | cut -d: -f1)
	printf '%s/fasta/%s/%s' "$WORK_DIR" "$n" "$(basename "${1%.gz}")"
}
# The FASTA an engine reads: a .gz source decompressed into the work tree at the start of
# every run (never a copy an earlier run left) with its index beside it, otherwise the
# source itself.
case_fasta() {
	case "$C_FASTA_PATH" in
	"") ;;
	*.gz) decompressed_fasta "$C_FASTA_PATH" ;;
	*) printf '%s' "$C_FASTA_PATH" ;;
	esac
}

if [[ "$DRY_RUN" == "true" ]]; then
	while IFS= read -r line; do
		[[ -n "$line" ]] || continue
		load_case "$line"
		fasta=$(case_fasta)
		case_work=$(case_work_dir)
		echo "DRY-RUN $C_ID cache species=$C_SPECIES perl=$C_PERL_SLICE candidate=$C_JSON_CACHE"
		[[ -z "$C_REFUSAL_DIR" ]] || echo "DRY-RUN $C_ID refusal: $C_REFUSAL_DIR (perl.txt, candidate.txt)"
		if [[ -z "$C_REFERENCE_DIR" ]]; then
			for f in $C_PERL_FORMATS; do
				echo "DRY-RUN $C_ID perl $f"
				perl_cmd "$f" "$case_work/perl" "$fasta"
				print_cmd
			done
		else
			echo "DRY-RUN $C_ID perl stored: $C_REFERENCE_DIR/"
		fi
		for ((i = 0; i < N_CANDIDATES; i++)); do
			pin="${PINS[$i]}"
			cand_dir="$case_work/candidates/$pin"
			for f in $C_OUTPUT_FORMATS; do
				echo "DRY-RUN $C_ID $pin $f"
				candidate_cmd "${CANDIDATES[$i]}" "$f" "$cand_dir" "$fasta"
				print_cmd
				[[ -z "$C_REFUSAL_DIR" ]] || continue
				expected=$(expected_output "$(perl_format "$f")")
				compare_cmd "$f" "$expected" "$cand_dir/$(format_file "$f")" "$OUT_DIR/$pin/$C_ID/$f.json"
				print_cmd
			done
			if [[ "$C_WARNINGS" == "true" && -z "$C_REFUSAL_DIR" ]]; then
				echo "DRY-RUN $C_ID $pin warnings"
				if [[ -n "$C_REFERENCE_DIR" ]]; then
					expected=$(stored_warnings_file)
				else
					expected="$case_work/perl/default.txt_warnings.txt"
				fi
				warnings_cmd "$expected" "$cand_dir/default.txt_warnings.txt" "$OUT_DIR/$pin/$C_ID/warnings.json"
				print_cmd
			fi
		done
	done <<<"$SELECTED"
	log "dry run complete: $N_SELECTED case(s) x $N_CANDIDATES candidate(s) planned, nothing written"
	exit 0
fi

mkdir -p "$OUT_DIR" "$WORK_DIR"
OUT_DIR=$(abs_path "$OUT_DIR")
WORK_DIR=$(abs_path "$WORK_DIR")

# An interrupted run stops its cases and their engines. Every engine and comparator runs
# as a background job waited on, since bash runs a trap only once a foreground command
# has returned; under --jobs the cases are background subshells, which bash leaves running
# when the parent dies, and each engine is a child of its case's subshell.
on_signal() {
	trap - INT TERM
	local pid
	for pid in $(jobs -p); do
		pkill -TERM -P "$pid" 2>/dev/null || true
		kill -TERM "$pid" 2>/dev/null || true
	done
	exit 130
}
trap on_signal INT TERM

# Decompress every .gz FASTA once per run.
while IFS= read -r src; do
	case "$src" in
	*.gz)
		target=$(decompressed_fasta "$src")
		mkdir -p "$(dirname "$target")"
		gzip -dc "$src" >"$target.tmp" || die "cannot decompress $src"
		mv -f "$target.tmp" "$target"
		cp -f "${src%.gz}.fai" "$target.fai"
		;;
	esac
done <<<"$FASTA_SOURCES"

# repo_digest <image>: the first RepoDigest of a local image, or nothing.
repo_digest() {
	docker image inspect --format '{{index .RepoDigests 0}}' "$1" 2>/dev/null || true
}

# The digest of every image a live case uses, pulled when absent; recorded per report.
IMAGE_DIGESTS=""
for image in $IMAGES; do
	digest=$(repo_digest "$image")
	if [[ -z "$digest" ]]; then
		log "pulling $image"
		docker pull "$image" >/dev/null || die "docker pull $image failed"
		digest=$(repo_digest "$image")
	fi
	[[ -n "$digest" ]] || die "docker image inspect reports no digest for $image"
	log "image $image digest $digest"
	IMAGE_DIGESTS="${IMAGE_DIGESTS}${image}	${digest}
"
done

# One case: Ensembl VEP once per format, then every candidate

# run_engine <log file> <cmd...>: runs, exit status returned, stderr kept for the report.
run_engine() {
	local errlog="$1"
	shift
	run_waited "$errlog.out" "$errlog" "$@"
}
# run_perl <log file>: runs CMD, the Ensembl VEP command. Its first --nearest run on a
# cache slice writes transcript_gene_tss.txt through a temporary file named after the
# process id (AnnotationSource/Cache/Transcript.pm:288, 317), which is the same in every
# container, so two first runs would tear one file: while the slice lacks it, the
# --nearest runs of this invocation hold one lock per slice. flock execs the engine
# (--no-fork) so that on_signal's pkill of a case's children reaches it; without flock
# the run is unlocked, which preflight allows only under --jobs 1 or where no two such
# cases share a slice.
run_perl() {
	local errlog="$1" slice_name lock
	if tss_first_run && command -v flock >/dev/null 2>&1; then
		mkdir -p "$WORK_DIR/locks"
		slice_name="$(basename "$(dirname "$C_PERL_SLICE")")_$(basename "$C_PERL_SLICE")"
		lock="$WORK_DIR/locks/$slice_name.tss.lock"
		run_engine "$errlog" flock --no-fork "$lock" "${CMD[@]}"
	else
		run_engine "$errlog" "${CMD[@]}"
	fi
}
# run_waited <stdout file> <stderr file> <cmd...>: runs the command as a job and waits for
# it, so a signal reaches the trap at once; returns the command's status.
run_waited() {
	local out="$1" err="$2" pid rc
	shift 2
	"$@" >"$out" 2>"$err" &
	pid=$!
	set +e
	wait "$pid"
	rc=$?
	set -e
	return "$rc"
}
error_text() { # error_text <prefix> <rc> <stderr file>
	local tail
	tail=$(tail -c 400 "$3" 2>/dev/null | tr '\n\t' '  ')
	printf '%s exited %s: %s' "$1" "$2" "${tail:-no stderr}"
}
# run_compare <name> <report> <log prefix>: runs CMD, the comparator, and prints the
# write_report argument for the result: format_report=<name>=<report>, or
# format_error=<name>=<text> when the comparator exited 2 or wrote no report.
run_compare() {
	local name="$1" report="$2" logs="$3" rc=0
	run_waited "$logs.compare.out" "$logs.compare.err" "${CMD[@]}" || rc=$?
	if [[ "$rc" -eq 2 || ! -f "$report" ]]; then
		printf 'format_error=%s=%s' "$name" "$(error_text "comparator" "$rc" "$logs.compare.err")"
	else
		printf 'format_report=%s=%s' "$name" "$report"
	fi
}

# write_report <report path> <key=value ...>: composes report.json and prints the verdict.
write_report() {
	"$PYTHON" - "$@" <<'PY'
import hashlib
import json
import sys

report_path = sys.argv[1]
record = {"formats": {}}
reasons = []
# Per compared format, the entries each requested mask class set aside (the warnings report
# counts its lines), read against the case's masks cell below, and per output format the
# classes its report names among those asked for and leaves out of its unused list.
applied_by_format = {}
used_by_format = {}
args = sys.argv[2:]
i = 0
while i < len(args):
    key, _, value = args[i].partition("=")
    i += 1
    if key == "format_report":
        fmt, _, path = value.partition("=")
        with open(path, encoding="utf-8") as fh:
            sub = json.load(fh)
        sub_reasons = list(sub.get("reasons", []))
        # Two outputs with nothing to compare agree vacuously; an engine that annotated
        # nothing must not pass as one that annotated everything alike.
        counts = sub.get("counts", {})
        if sub.get("mode") == "compare" and counts.get("expected_entries", 0) == 0 and counts.get("actual_entries", 0) == 0:
            sub_reasons.append("no entry compared: both outputs carry no annotation")
        if sub.get("verdict") != "green" and not sub_reasons:
            sub_reasons.append(f"comparator verdict {sub.get('verdict')!r} without a stated reason")
        entry = {"verdict": "red" if sub_reasons else "green", "report": path, "reasons": sub_reasons}
        if "one_sided" in sub:
            entry["one_sided"] = sub["one_sided"]
        if sub.get("floors_not_applicable"):
            entry["floors_not_applicable"] = sub["floors_not_applicable"]
        if sub.get("mode") == "compare":
            masks_report = sub.get("masks", {})
            applied_by_format[fmt] = dict(masks_report.get("applied", {}))
            if "unused" in masks_report:
                used_by_format[fmt] = set(masks_report.get("requested", [])) - set(masks_report["unused"])
        elif sub.get("mode") == "warnings":
            applied_by_format[fmt] = dict(sub.get("warnings", {}).get("masks_applied", {}))
        if fmt in applied_by_format:
            entry["masks_applied"] = applied_by_format[fmt]
        record["formats"][fmt] = entry
        reasons.extend(f"{fmt}: {r}" for r in sub_reasons)
    elif key == "format_error":
        fmt, _, message = value.partition("=")
        record["formats"][fmt] = {"verdict": "red", "error": message}
        reasons.append(f"{fmt}: {message}")
    elif key == "format_refusal":
        # <format>=<Ensembl VEP exit>\t<its text matched>\t<candidate exit>\t<matched>\t<Ensembl
        # VEP stderr tail>\t<candidate stderr tail>: both engines must exit non-zero with the
        # expected text in their stderr.
        fmt, _, fields = value.partition("=")
        perl_rc, perl_matched, cand_rc, cand_matched, perl_tail, cand_tail = fields.split("\t")
        sub_reasons = []
        entry = {"refusal": {}}
        for engine, label, rc, matched, tail in (
            ("perl", "Ensembl VEP", perl_rc, perl_matched, perl_tail),
            ("candidate", "candidate", cand_rc, cand_matched, cand_tail),
        ):
            entry["refusal"][engine] = {"exit": int(rc), "text_matched": matched == "1"}
            if rc == "0":
                sub_reasons.append(f"{label} exited 0 where a refusal was expected")
            elif matched != "1":
                sub_reasons.append(f"{label} exited {rc} without the expected refusal text: {tail or 'no stderr'}")
        entry["verdict"] = "red" if sub_reasons else "green"
        entry["reasons"] = sub_reasons
        record["formats"][fmt] = entry
        reasons.extend(f"{fmt}: {r}" for r in sub_reasons)
    elif key == "flags":
        record[key] = value.split()
    elif key == "masks":
        # The cell as written; its classes and their expected counts are checked once every
        # format is read.
        record[key] = value
    elif key == "unordered":
        record[key] = value == "1"
    elif key == "formats_requested":
        record["formats_requested"] = value.split(",")
    elif key == "reference_provenance":
        record[key] = value or None
        if value:
            try:
                with open(value, encoding="utf-8") as fh:
                    provenance = json.load(fh)
                if not isinstance(provenance, dict):
                    raise ValueError("not a JSON object")
            except (OSError, ValueError) as exc:
                reasons.append(f"stored reference provenance.json unreadable: {exc}")
                provenance = {}
            record["reference_provenance_image"] = provenance.get("vep_image")
    elif key in ("perl_image_digest", "perl_server", "reference_dir", "refusal_dir", "sources_dir", "fasta", "manifest", "candidate_dir_cache", "candidate_version"):
        record[key] = value or None
    else:
        record[key] = value
# A stored reference carries its image in its provenance; a live run records the digest
# docker reports.
if not record.get("perl_image_digest") and record.get("reference_provenance_image"):
    record["perl_image_digest"] = record["reference_provenance_image"]
# A floor a format could not count is met only where another output format of the same case
# counted it, one whose report does not set that floor aside too; a case that gives the floor
# no such format leaves it unchecked.
floors_set_aside = {
    fmt: set(entry.get("floors_not_applicable", []))
    for fmt, entry in record["formats"].items()
    if fmt in ("default", "tab", "vcf", "parquet", "json") and "report" in entry
}
for fmt, entry in sorted(record["formats"].items()):
    for name in entry.get("floors_not_applicable", []):
        counted_elsewhere = any(name not in aside for other, aside in floors_set_aside.items() if other != fmt)
        if not counted_elsewhere:
            reasons.append(f"{fmt}: floor {name} is not applicable to {fmt} and no other format of the case counted it")
# A class the cell names must be used in some compared format: it sets something aside
# there, or the format's report leaves it out of its unused list (a class a value set aside
# under another class needs, counted under that other class alone). A class with an expected
# count must apply to that many entries in every output format compared (the warnings report
# counts lines and is not held to it).
masks_cell = record.get("masks", "-")
if masks_cell != "-" and applied_by_format:
    for item in masks_cell.split(","):
        klass, _, expected_count = item.partition("=")
        counted = any(applied.get(klass) for applied in applied_by_format.values())
        if not counted and not any(klass in used for used in used_by_format.values()):
            reasons.append(f"mask {klass} is requested and applied in no compared format")
            continue
        if not expected_count:
            continue
        for fmt, applied in sorted(applied_by_format.items()):
            if fmt == "warnings":
                continue
            count = applied.get(klass, 0)
            if count != int(expected_count):
                reasons.append(f"{fmt}: mask {klass} applied to {count} entries, {expected_count} expected")
record["verdict"] = "red" if reasons else "green"
record["reasons"] = reasons
with open(report_path, "w", encoding="utf-8") as fh:
    json.dump(record, fh, indent=1, sort_keys=True)
    fh.write("\n")
print(record["verdict"])
PY
}

run_case() {
	load_case "$1"
	local case_work fasta perl_dir f rc started
	case_work=$(case_work_dir)
	local perl_image_digest="" perl_errors=""
	rm -rf "$case_work"
	mkdir -p "$case_work/perl"
	: >"$case_work/status.rows"
	fasta=$(case_fasta)
	started=$(date -u +%Y-%m-%dT%H:%M:%SZ)

	# A refusal case keeps, per Ensembl VEP format, the exit status, whether the stderr carried
	# perl.txt's text, and the stderr's tail, for the candidate's formats to be judged against.
	local reference_provenance="" perl_text="" cand_text="" perl_rcs="" perl_matches="" perl_tails="" matched
	if [[ -n "$C_REFERENCE_DIR" ]]; then
		perl_dir="$C_REFERENCE_DIR"
		[[ ! -f "$(dirname "$C_REFERENCE_DIR")/provenance.json" ]] || reference_provenance="$(dirname "$C_REFERENCE_DIR")/provenance.json"
		log "[$C_ID] Ensembl VEP: stored reference $perl_dir"
	else
		perl_dir="$case_work/perl"
		perl_image_digest=$(tsv_lookup "$IMAGE_DIGESTS" "$C_PERL_IMAGE")
		if [[ -n "$C_REFUSAL_DIR" ]]; then
			perl_text=$(<"$C_REFUSAL_DIR/perl.txt")
			cand_text=$(<"$C_REFUSAL_DIR/candidate.txt")
		fi
		for f in $C_PERL_FORMATS; do
			perl_cmd "$f" "$perl_dir" "$fasta"
			log "[$C_ID] Ensembl VEP $f ($C_PERL_MODE${C_PERL_SERVER:+ $C_PERL_SERVER})"
			rc=0
			run_perl "$case_work/perl.$f.stderr" || rc=$?
			if [[ -n "$C_REFUSAL_DIR" ]]; then
				matched=0
				if [[ "$rc" -ne 0 && "$(<"$case_work/perl.$f.stderr")" == *"$perl_text"* ]]; then
					matched=1
				fi
				perl_rcs="${perl_rcs}${f}	${rc}
"
				perl_matches="${perl_matches}${f}	${matched}
"
				perl_tails="${perl_tails}${f}	$(tail -c 400 "$case_work/perl.$f.stderr" 2>/dev/null | tr '\n\t' '  ')
"
				log "[$C_ID] Ensembl VEP $f exited $rc; refusal text matched: $matched"
			elif [[ "$rc" -ne 0 ]]; then
				perl_errors="${perl_errors}${f}	$(error_text "Ensembl VEP" "$rc" "$case_work/perl.$f.stderr")
"
				echo "WARN: [feature_oracle] $C_ID: Ensembl VEP $f exited $rc" >&2
			elif [[ ! -e "$perl_dir/$(format_file "$f")" ]]; then
				perl_errors="${perl_errors}${f}	Ensembl VEP exited 0 and wrote no $(format_file "$f")
"
			fi
		done
	fi

	# format_result <candidate> <format> <candidate dir> <report dir> <fasta or empty>: the
	# candidate's run of one output format of the loaded case, then its comparison; prints the
	# write_report argument for the result (format_refusal, format_error or run_compare's
	# line). Reads run_case's locals (perl_errors and, for a refusal case, cand_text,
	# perl_rcs, perl_matches and perl_tails) through dynamic scoping, and is defined inside
	# run_case so that the script's top-level functions stay its interface.
	format_result() {
		local bin="$1" f="$2" cand_dir="$3" report_dir="$4" fasta="$5"
		local rc=0 pf perl_error expected actual matched perl_rc perl_matched perl_tail cand_tail
		candidate_cmd "$bin" "$f" "$cand_dir" "$fasta"
		run_engine "$cand_dir/$f.stderr" "${CMD[@]}" || rc=$?
		pf=$(perl_format "$f")
		if [[ -n "$C_REFUSAL_DIR" ]]; then
			matched=0
			if [[ "$rc" -ne 0 && "$(<"$cand_dir/$f.stderr")" == *"$cand_text"* ]]; then
				matched=1
			fi
			perl_rc=$(tsv_lookup "$perl_rcs" "$pf")
			perl_matched=$(tsv_lookup "$perl_matches" "$pf")
			perl_tail=$(tsv_lookup "$perl_tails" "$pf")
			cand_tail=$(tail -c 400 "$cand_dir/$f.stderr" 2>/dev/null | tr '\n\t' '  ')
			printf 'format_refusal=%s=%s\t%s\t%s\t%s\t%s\t%s' \
				"$f" "$perl_rc" "$perl_matched" "$rc" "$matched" "$perl_tail" "$cand_tail"
			return
		fi
		if [[ "$rc" -ne 0 ]]; then
			printf 'format_error=%s=%s' "$f" "$(error_text "candidate" "$rc" "$cand_dir/$f.stderr")"
			return
		fi
		perl_error=$(tsv_lookup "$perl_errors" "$pf")
		if [[ -n "$perl_error" ]]; then
			printf 'format_error=%s=%s' "$f" "$perl_error"
			return
		fi
		expected=$(expected_output "$pf")
		actual="$cand_dir/$(format_file "$f")"
		compare_cmd "$f" "$expected" "$actual" "$report_dir/$f.json"
		run_compare "$f" "$report_dir/$f.json" "$cand_dir/$f"
	}
	# An interrupted format job stops the engine or comparator it waits on.
	on_format_job_signal() {
		trap - INT TERM
		local pid
		for pid in $(jobs -p); do
			kill -TERM "$pid" 2>/dev/null || true
		done
		exit 130
	}
	local i pin bin cand_dir report_dir perl_error expected actual report_args job_tail verdict digest
	for ((i = 0; i < N_CANDIDATES; i++)); do
		pin="${PINS[$i]}"
		bin="${CANDIDATES[$i]}"
		cand_dir="$case_work/candidates/$pin"
		report_dir="$OUT_DIR/$pin/$C_ID"
		rm -rf "$report_dir"
		mkdir -p "$cand_dir" "$report_dir"
		# A format's result lands in <cand_dir>/<format>.result. The format's start is logged
		# by the run once the job has a slot, not inside the job, so the log shows the formats
		# in flight while --format-jobs runs them; what a job itself writes goes to
		# <format>.job.out and .job.err beside its outputs and is replayed in format order once
		# every job has ended.
		for f in $C_OUTPUT_FORMATS; do
			if [[ "$FORMAT_JOBS" -eq 1 ]]; then
				log "[$C_ID] $pin $f"
				format_result "$bin" "$f" "$cand_dir" "$report_dir" "$fasta" >"$cand_dir/$f.result"
			else
				while [[ "$(running_jobs)" -ge "$FORMAT_JOBS" ]]; do sleep 1; done
				log "[$C_ID] $pin $f"
				(
					trap on_format_job_signal INT TERM
					format_result "$bin" "$f" "$cand_dir" "$report_dir" "$fasta" >"$cand_dir/$f.result"
				) >"$cand_dir/$f.job.out" 2>"$cand_dir/$f.job.err" &
			fi
		done
		if [[ "$FORMAT_JOBS" -gt 1 ]]; then
			wait
			for f in $C_OUTPUT_FORMATS; do
				cat "$cand_dir/$f.job.out"
				cat "$cand_dir/$f.job.err" >&2
			done
		fi
		# A job that died leaves its result file empty (the redirection opened it first); its
		# stderr's tail travels in the reason, since the candidate directory goes with the case.
		report_args=()
		for f in $C_OUTPUT_FORMATS; do
			if [[ -s "$cand_dir/$f.result" ]]; then
				report_args+=("$(<"$cand_dir/$f.result")")
			else
				job_tail=$(tail -c 400 "$cand_dir/$f.job.err" 2>/dev/null | tr '\n\t' '  ')
				report_args+=("format_error=$f=the format's job ended without a result: ${job_tail:-no stderr}")
			fi
		done
		if [[ "$C_WARNINGS" == "true" ]]; then
			perl_error=$(tsv_lookup "$perl_errors" default)
			if [[ -n "$perl_error" ]]; then
				report_args+=("format_error=warnings=$perl_error")
			elif [[ ! -e "$cand_dir/default.txt" ]]; then
				report_args+=("format_error=warnings=the candidate's default run did not complete")
			else
				if [[ -n "$C_REFERENCE_DIR" ]]; then
					expected=$(stored_warnings_file)
				else
					expected="$perl_dir/default.txt_warnings.txt"
					[[ -e "$expected" ]] || : >"$expected"
				fi
				actual="$cand_dir/default.txt_warnings.txt"
				[[ -e "$actual" ]] || : >"$actual"
				warnings_cmd "$expected" "$actual" "$report_dir/warnings.json"
				report_args+=("$(run_compare warnings "$report_dir/warnings.json" "$cand_dir/warnings")")
			fi
		fi
		verdict=$(write_report "$report_dir/report.json" \
			"case=$C_ID" "tier=$C_TIER" "input=$C_INPUT_PATH" "cache=$C_CACHE" "species=$C_SPECIES" \
			"candidate_json_cache=$C_JSON_CACHE" "candidate_dir_cache=$C_DIR_CACHE" "fasta=$fasta" "flags=$C_FLAGS" \
			"formats_requested=$C_FORMATS" "manifest=$C_MANIFEST" "masks=$C_MASKS" \
			"floors=$C_FLOORS" "field_pairs=$C_FIELD_PAIRS" "unordered=$C_UNORDERED" \
			"perl_image=$C_PERL_IMAGE" "perl_image_digest=$perl_image_digest" "perl_mode=$C_PERL_MODE" \
			"perl_server=$C_PERL_SERVER" "reference_release=$C_REFERENCE_RELEASE" \
			"reference_dir=$C_REFERENCE_DIR" "reference_provenance=$reference_provenance" \
			"refusal_dir=$C_REFUSAL_DIR" "sources_dir=$SOURCES_DIR" \
			"candidate=$bin" "candidate_pin=$pin" "candidate_version=$("$bin" --version 2>/dev/null | head -1 || true)" \
			"run_utc=$RUN_UTC" "started_utc=$started" "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
			"${report_args[@]}")
		digest=$(sha256_of "$report_dir/report.json")
		printf '%s\t%s\t%s\t%s\t%s\n' "$C_ID" "$pin" "$verdict" "$digest" "$RUN_UTC" >>"$case_work/status.rows"
		log "[$C_ID] $pin: $verdict"
		[[ "$KEEP_OUTPUTS" == "true" ]] || rm -rf "$cand_dir"
	done
	[[ "$KEEP_OUTPUTS" == "true" ]] || rm -rf "$case_work/perl"
}

# All cases, --jobs at a time
running_jobs() {
	local n
	n=$(jobs -rp | grep -c .) || n=0
	echo "$n"
}
CASE_NAMES=""
while IFS= read -r line; do
	[[ -n "$line" ]] || continue
	cid=$(cell "$line" case)
	CASE_NAMES="$CASE_NAMES $cid"
	if [[ "$JOBS" -eq 1 ]]; then
		run_case "$line"
	else
		while [[ "$(running_jobs)" -ge "$JOBS" ]]; do sleep 1; done
		mkdir -p "$WORK_DIR/logs"
		(run_case "$line") >"$WORK_DIR/logs/$cid.log" 2>&1 &
	fi
done <<<"$SELECTED"
if [[ "$JOBS" -gt 1 ]]; then
	wait
	for cid in $CASE_NAMES; do
		cat "$WORK_DIR/logs/$cid.log"
	done
fi

# status.tsv: one block per candidate, the cases in table order; a case that left no
# row for a candidate (the runner stopped inside it) is red.
STATUS="$OUT_DIR/status.tsv"
printf 'case\tcandidate_pin\tverdict\treport_digest\trun_utc\n' >"$STATUS.tmp"
N_RED=0
N_ROWS=0
for ((i = 0; i < N_CANDIDATES; i++)); do
	pin="${PINS[$i]}"
	for cid in $CASE_NAMES; do
		row=$(awk -F'\t' -v p="$pin" '$2 == p' "$WORK_DIR/cases/$cid/status.rows" 2>/dev/null || true)
		if [[ -z "$row" ]]; then
			row=$(printf '%s\t%s\tred\t-\t%s' "$cid" "$pin" "$RUN_UTC")
			echo "ERROR: [feature_oracle] $cid: no verdict for $pin; the case did not complete" >&2
		fi
		printf '%s\n' "$row" >>"$STATUS.tmp"
		N_ROWS=$((N_ROWS + 1))
		[[ "$row" == "$cid	$pin	green	"* ]] || N_RED=$((N_RED + 1))
	done
done
mv -f "$STATUS.tmp" "$STATUS"

"$PYTHON" - "$OUT_DIR/run.json" "$RUN_UTC" "$TIER" "$CASES" "$REFERENCE_DIR" "$SOURCES_DIR" "$PERL_FORK" "$FORK" "$IMAGE_DIGESTS" "$N_ROWS" "$N_RED" "${PINS[@]}" <<'PY'
import json
import sys

(path, run_utc, tier, cases, reference_dir, sources_dir, perl_fork, fork, digests, rows, red) = sys.argv[1:12]
record = {
    "run_utc": run_utc,
    "tier": tier,
    "case_table": cases,
    "reference_dir": reference_dir or None,
    "sources_dir": sources_dir or None,
    "perl_fork": int(perl_fork),
    "candidate_fork": int(fork) if fork else None,
    "image_digests": dict(line.split("\t", 1) for line in digests.splitlines() if line),
    "candidate_pins": sys.argv[12:],
    "status_rows": int(rows),
    "red": int(red),
}
with open(path, "w", encoding="utf-8") as fh:
    json.dump(record, fh, indent=1, sort_keys=True)
    fh.write("\n")
PY

log "feature_oracle complete: $N_ROWS verdict(s), $N_RED red; $STATUS"
[[ "$N_RED" -eq 0 ]] || exit 1
