#!/usr/bin/env bash
#
# feature_oracle.sh -- run the comparison cases of a case table: Ensembl VEP against one
# or more candidate `vep` binaries on the same input and flags, every requested output
# format compared column by column by compare_vep_formats.py, one report per candidate
# and case, one status table per run.
#
# THE CASE TABLE is a tab-separated file whose header names these columns, in any order.
# A cell with no value is `-`; blank lines and `#` lines are skipped.
#
#   case               unique name: the report directory and the status.tsv key
#   tier               fixture | population; --tier selects the rows of one tier
#   input              the input file: absolute, else relative to --data-dir, else to the
#                      repository root
#   cache              <version>_<assembly>, the Ensembl cache slice both engines annotate
#                      from: homo_sapiens/<version>_<assembly> under --perl-cache-dir
#                      (homo_sapiens_refseq or _merged when the flags say so) and
#                      <assembly, lower-cased> under --json-cache-dir; an optional
#                      =<dir> suffix names the candidate's JSON cache instead
#   fasta              - or an indexed FASTA given to both engines (its .fai beside it; a
#                      .gz is decompressed once per run, as the golden tests do, into a
#                      directory of its own under the work tree, so two sources that
#                      decompress to one file name never share a copy)
#   flags              whitespace-separated flags appended verbatim to both engines, or -;
#                      a --nearest case mounts the Ensembl cache writable, since Ensembl
#                      VEP writes transcript_gene_tss.txt there on its first run (a scan
#                      of every region file) and keeps the file, built under that run's
#                      transcript filters, for every later run; while a slice lacks the
#                      file, its --nearest runs take turns under one lock, because the
#                      temporary file is named after the process id, which is the same
#                      in every container; its comparisons settle NEAREST by start site
#                      through the candidate's JSON cache (compare_vep_formats.py
#                      --nearest-sites)
#   formats            comma-separated subset of default,tab,vcf,json,parquet,warnings;
#                      warnings compares the default run's `_warnings.txt` files and so
#                      needs default
#   perl_image         the Ensembl VEP image, pinned by digest (name@sha256:...)
#   perl_mode          offline | database
#   reference_release  the Ensembl VEP release the image carries, recorded per report
#   masks              compare_vep_formats.py --masks, applied from the manifest.json
#                      beside the input, or -; the comparator's built-in classes
#                      (regulatory_feature_rows, order_only, individual_allele_order) need
#                      no manifest
#   floors             compare_vep_formats.py --floors (COL=MIN,...), or -; a floor on a
#                      column the json format has no row count for (a record-level column,
#                      Feature_type, Extra) is not applicable there: the report.json entry
#                      of that format lists it under floors_not_applicable, and a case whose
#                      formats give such a floor no text format to be checked in is red
#   field_pairs        compare_vep_formats.py --field-pairs (COL=EXPECTED=>ACTUAL;...), or -
#
# THE ENGINES. Ensembl VEP runs in Docker, once per case and output format, since a run
# writes one format (Runner.pm:871 at release/116.2 builds one OutputFactory per run):
#
#   vep -i <input> -o <out> --offline --cache --dir_cache <cache> --cache_version <v>
#       --species homo_sapiens --assembly <asm> --format vcf --fork 1 --no_stats --quiet
#       --force_overwrite <flags> [--fasta <fa>] [--tab|--vcf|--json]
#
# In database mode `--offline --cache --dir_cache --cache_version` become
# `--database --db_version <v>`: Config.pm:618 refuses --database beside --offline or
# --cache, and the HGVS, ID and SPDI input parsers (Parser/HGVS.pm:104, Parser/ID.pm:100,
# Parser/SPDI.pm:96) refuse to run offline. The server is Ensembl's public one, the
# Config.pm defaults (host ensembldb.ensembl.org, port 3306; Config.pm:711 moves GRCh37
# to 3337), recorded in the report. The candidate runs with the same flags and
# `--offline --json_cache <cache>`; its parquet run is `--output_format parquet`,
# compared with the Ensembl VEP tab output. `--format vcf` is given to both engines
# unless the flags name an input format. The warnings compared are the default runs'
# `<out>_warnings.txt`, which Ensembl VEP opens on the first warning
# (BaseVEP.pm:878-897), so an absent file on either side means no warnings. A warning an
# Ensembl VEP worker raises while annotating is written with the worker's pid under
# --fork (Runner.pm:563 prefixes the captured worker stderr with it), so it differs from
# run to run; a record skipped while parsing is written by the parent without one.
#
# --reference-dir replaces the Ensembl VEP runs with a stored reference per case:
# <reference-dir>/<case>/expected/{default.txt,tab.txt,vcf.vcf,json.jsonl}[.gz],
# the layout of a golden corpus, plus expected/default.txt_warnings.txt[.gz] when the
# case names warnings (a stored reference states "no warnings" with an empty file);
# a provenance.json beside expected/ is recorded. No container runs in that mode.
#
# A CASE IS GREEN for a candidate when every requested format's comparison is green:
# no entry on one side only, every compared field equal after the masks and field
# pairs, every floor met, and for warnings the same (line, description) set. A format
# whose engine or comparison fails to run is red, and so is one whose two outputs both
# carry no annotation row: two empty outputs agree on nothing. The run exits 1 on any red.
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
#     [--work-dir <dir>] [--keep-outputs] [--fork <N>] [--perl-fork <N>] \
#     [--jobs <N>] [--dry-run]
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
WORK_DIR=""
KEEP_OUTPUTS="false"
FORK=""
PERL_FORK=1
JOBS=1
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
  --perl-cache-dir <dir>    Directory holding homo_sapiens/<version>_<assembly>/
                            (default <data>/caches/perl/vep-cache)
  --json-cache-dir <dir>    Directory holding the candidate's <assembly>/ JSON caches
                            (default <data>/caches/vep-rs)
  --reference-dir <dir>     Stored Ensembl VEP output per case
                            (<dir>/<case>/expected/); runs no container
  --work-dir <dir>          Where the engines write (default <out>/work)
  --keep-outputs            Keep the engines' outputs after the comparison
  --fork <N>                The candidate's --fork (default: the binary's own)
  --perl-fork <N>           Ensembl VEP's --fork (default 1)
  --jobs <N>                Cases run at once (default 1)
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
for _n in "$PERL_FORK" "$JOBS" ${FORK:+"$FORK"}; do
	[[ "$_n" =~ ^[0-9]+$ && "$_n" -gt 0 ]] || die "--fork, --perl-fork and --jobs take a positive integer (got '$_n')"
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
[[ -f "$COMPARATOR" ]] || die "comparator not found beside this script: $COMPARATOR"
command -v "$PYTHON" >/dev/null 2>&1 || die "$PYTHON is required (the comparator and the report writer run on it)"
if [[ "$DRY_RUN" != "true" && -z "$REFERENCE_DIR" ]]; then
	command -v docker >/dev/null 2>&1 || die "docker is required to run Ensembl VEP (pass --reference-dir for a stored reference, or --dry-run)"
fi

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
COLUMNS="case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs"
OUTPUT_FORMATS="default tab vcf json parquet"
# The comparator's own divergence classes, applied without a manifest.
BUILTIN_MASKS="regulatory_feature_rows order_only individual_allele_order"
# Flags the runner sets from the table's own columns and its options, with the aliases
# Ensembl VEP's Config.pm and the candidate accept for them; a flags cell naming one would
# give an engine the flag twice, or move the input, the server or the assembly behind the
# report's back.
RUNNER_FLAGS="-i --input_file --input_data --id -o --output_file --offline --cache --database --db_version --host --port --user --password --pass --dir --dir_cache --cache_version --json_cache -s --species -a --assembly --grch37 --fasta --fa --tab --vcf --json --output_format --force_overwrite --force --no_stats -q --quiet --fork --threads --warning_file --skipped_variants_file --config"

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
[[ "$N_COLUMNS" -eq 13 ]] || die "$CASES: header has $N_COLUMNS columns; the 13 are: $COLUMNS"

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
	[[ "$nf" -eq 13 ]] || die "$CASES:$line_no: $nf columns, expected 13 (an empty cell is written -)"
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

# has_flag <flag>: whether the loaded case's flags carry <flag>, as a word or as <flag>=value.
has_flag() {
	local f
	for f in ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}; do
		[[ "$f" == "$1" || "$f" == "$1="* ]] && return 0
	done
	return 1
}

# load_case <line>: the C_* globals of one row.
load_case() {
	local line="$1" cache_spec species_dir
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
	C_FLOORS=$(cell "$line" floors)
	C_FIELD_PAIRS=$(cell "$line" field_pairs)
	cache_spec=$(cell "$line" cache)
	[[ "$C_FLAGS" != "-" ]] || C_FLAGS=""
	# The flags as words: `read -a` splits on blanks without expanding a glob, so a `*` in
	# a flag value reaches the engines as written.
	C_FLAG_WORDS=()
	[[ -z "$C_FLAGS" ]] || read -r -a C_FLAG_WORDS <<<"$C_FLAGS"

	C_INPUT_PATH=$(resolve_path "$C_INPUT")
	C_MANIFEST="$(dirname "$C_INPUT_PATH")/manifest.json"
	[[ -f "$C_MANIFEST" ]] || C_MANIFEST=""
	C_FASTA_PATH=""
	[[ "$C_FASTA" == "-" ]] || C_FASTA_PATH=$(resolve_path "$C_FASTA")

	C_CACHE="${cache_spec%%=*}"
	C_CACHE_VERSION="${C_CACHE%%_*}"
	C_ASSEMBLY="${C_CACHE#*_}"
	species_dir="homo_sapiens"
	if has_flag "--merged"; then
		species_dir="homo_sapiens_merged"
	elif has_flag "--refseq"; then
		species_dir="homo_sapiens_refseq"
	fi
	C_PERL_SLICE="$PERL_CACHE_DIR/$species_dir/$C_CACHE"
	if [[ "$cache_spec" == *=* ]]; then
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
	# host ensembldb.ensembl.org and port 3306; Config.pm:711-713 set port 3337 for GRCh37.
	C_PERL_SERVER=""
	if [[ "$C_PERL_MODE" == "database" ]]; then
		if [[ "$(lower "$C_ASSEMBLY")" == "grch37" ]]; then
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
# stored_reference_file <case> <format>: the stored output, plain or .gz, or nothing.
stored_reference_file() {
	plain_or_gz "$REFERENCE_DIR/$1/expected/$(format_file "$2")"
}
# stored_warnings_file <case>: the stored default run's warnings file, plain or .gz, or nothing.
stored_warnings_file() {
	plain_or_gz "$REFERENCE_DIR/$1/expected/default.txt_warnings.txt"
}
# case_work_dir: the loaded case's directory under the work tree. Cases live under cases/,
# candidates under <case>/candidates/, so no case or pin name can reach the Ensembl VEP
# outputs (<case>/perl/), the decompressed FASTAs (fasta/) or the logs (logs/).
case_work_dir() { printf '%s/cases/%s' "$WORK_DIR" "$C_ID"; }

# expected_output <perl format>: the Ensembl VEP output a comparison of the loaded case
# reads: the stored reference, or the file its run writes under the work tree.
expected_output() {
	if [[ -n "$REFERENCE_DIR" ]]; then
		stored_reference_file "$C_ID" "$1"
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
# transcript_gene_tss.txt, one line per case.
FASTA_SOURCES=""
TSS_SLICES=""
IMAGES=""
CACHE_RE='^[0-9]+_[A-Za-z0-9]+$'
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
	[[ -d "$C_JSON_CACHE/transcripts" ]] || pf_error "candidate JSON cache has no transcripts/: $C_JSON_CACHE"
	case "$C_PERL_MODE" in
	offline)
		if [[ -z "$REFERENCE_DIR" && ! -d "$C_PERL_SLICE" ]]; then
			pf_error "Ensembl VEP cache slice not found: $C_PERL_SLICE"
		elif [[ -z "$REFERENCE_DIR" && ! -e "$C_PERL_SLICE/transcript_gene_tss.txt" ]] && has_flag --nearest; then
			TSS_SLICES="${TSS_SLICES}${C_PERL_SLICE}
"
		fi
		;;
	database) ;;
	*) pf_error "perl_mode '$C_PERL_MODE' is not offline or database" ;;
	esac
	[[ "$C_PERL_IMAGE" != "-" ]] || pf_error "perl_image is -"
	has_word "$IMAGES" "$C_PERL_IMAGE" || IMAGES="$IMAGES $C_PERL_IMAGE"
	[[ "$C_REFERENCE_RELEASE" != "-" ]] || pf_error "reference_release is -"
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
	for f in ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}; do
		has_word "$RUNNER_FLAGS" "${f%%=*}" && pf_error "flags names $f, which the runner sets from the table's columns and its own options"
	done
	if [[ -n "$REFERENCE_DIR" ]]; then
		for f in $C_PERL_FORMATS; do
			stored_reference_file "$C_ID" "$f" >/dev/null || pf_error "stored reference lacks $(format_file "$f")[.gz] under $REFERENCE_DIR/$C_ID/expected/"
		done
		if [[ "$C_WARNINGS" == "true" ]]; then
			stored_warnings_file "$C_ID" >/dev/null || pf_error "stored reference lacks default.txt_warnings.txt[.gz] under $REFERENCE_DIR/$C_ID/expected/ (an empty file states no warnings)"
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
	for f in $(printf '%s' "$C_MASKS" | tr ',' ' '); do
		if [[ "$f" != "-" && -z "$C_MANIFEST" ]] && ! has_word "$BUILTIN_MASKS" "$f"; then
			pf_error "masks names '$f', which applies from a manifest.json beside the input, and none is there (the built-in classes are $BUILTIN_MASKS)"
		fi
	done
done <<<"$SELECTED"
[[ "$PREFLIGHT_FAILED" -eq 0 ]] || die "preflight failed; nothing ran"
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
		if has_flag --nearest; then
			CMD+=(-v "${PERL_CACHE_DIR}:/work/cache")
		else
			CMD+=(-v "${PERL_CACHE_DIR}:/work/cache:ro")
		fi
	fi
	if [[ -n "$fasta" ]]; then
		CMD+=(-v "${fasta}:/work/ref.fa:ro" -v "${fasta}.fai:/work/ref.fa.fai:ro")
	fi
	CMD+=("$C_PERL_IMAGE" vep -i "/work/input/${in_name}" -o "/work/out/$(format_file "$fmt")")
	if [[ "$C_PERL_MODE" == "offline" ]]; then
		CMD+=(--offline --cache --dir_cache /work/cache --cache_version "$C_CACHE_VERSION")
	else
		CMD+=(--database --db_version "$C_CACHE_VERSION")
	fi
	CMD+=(--species homo_sapiens --assembly "$C_ASSEMBLY")
	has_flag --format || CMD+=(--format vcf)
	CMD+=(--fork "$PERL_FORK" --no_stats --quiet --force_overwrite)
	CMD+=(${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"})
	[[ -z "$fasta" ]] || CMD+=(--fasta /work/ref.fa)
	case "$fmt" in
	default) ;;
	*) CMD+=("--$fmt") ;;
	esac
}

# candidate_cmd <binary> <format> <host out dir> <host fasta or empty>
candidate_cmd() {
	local bin="$1" fmt="$2" out_dir="$3" fasta="$4"
	CMD=("$bin" -i "$C_INPUT_PATH" -o "$out_dir/$(format_file "$fmt")"
		--offline --json_cache "$C_JSON_CACHE" --species homo_sapiens --assembly "$C_ASSEMBLY")
	has_flag --format || CMD+=(--format vcf)
	[[ -z "$FORK" ]] || CMD+=(--fork "$FORK")
	CMD+=(--no_stats --quiet --force_overwrite)
	CMD+=(${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"})
	[[ -z "$fasta" ]] || CMD+=(--fasta "$fasta")
	case "$fmt" in
	default) ;;
	parquet) CMD+=(--output_format parquet) ;;
	*) CMD+=("--$fmt") ;;
	esac
}

# compare_cmd <format> <expected> <actual> <report>
compare_cmd() {
	local fmt="$1"
	CMD=("$PYTHON" "$COMPARATOR" --format "$fmt" --expected "$2" --actual "$3" --report "$4")
	case "$fmt" in
	vcf | json) CMD+=(--input "$C_INPUT_PATH") ;;
	esac
	[[ -z "$C_MANIFEST" ]] || CMD+=(--manifest "$C_MANIFEST")
	[[ "$C_MASKS" == "-" ]] || CMD+=(--masks "$C_MASKS")
	[[ "$C_FLOORS" == "-" ]] || CMD+=(--floors "$C_FLOORS")
	[[ "$C_FIELD_PAIRS" == "-" ]] || CMD+=(--field-pairs "$C_FIELD_PAIRS")
	# NEAREST compares by transcription start site: Ensembl VEP names one transcript at a
	# winning site and the candidate every transcript there (compare_vep_formats.py).
	if has_flag --nearest; then
		CMD+=(--nearest-sites "$C_JSON_CACHE")
	fi
}
warnings_cmd() {
	CMD=("$PYTHON" "$COMPARATOR" --warnings --expected "$1" --actual "$2" --report "$3")
}

print_cmd() {
	printf '%q ' "${CMD[@]}"
	echo
}

# Run setup
RUN_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
N_CANDIDATES="${#CANDIDATES[@]}"
if [[ "$TIER" == "population" && -z "$REFERENCE_DIR" ]]; then
	echo "WARN: [feature_oracle] population tier without --reference-dir: Ensembl VEP runs live on every case" >&2
fi
log "feature_oracle tier=$TIER cases=$N_SELECTED candidates=$N_CANDIDATES table=$CASES dry_run=$DRY_RUN"
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
		if [[ -z "$REFERENCE_DIR" ]]; then
			for f in $C_PERL_FORMATS; do
				echo "DRY-RUN $C_ID perl $f"
				perl_cmd "$f" "$case_work/perl" "$fasta"
				print_cmd
			done
		else
			echo "DRY-RUN $C_ID perl stored: $REFERENCE_DIR/$C_ID/expected/"
		fi
		for ((i = 0; i < N_CANDIDATES; i++)); do
			pin="${PINS[$i]}"
			cand_dir="$case_work/candidates/$pin"
			for f in $C_OUTPUT_FORMATS; do
				echo "DRY-RUN $C_ID $pin $f"
				candidate_cmd "${CANDIDATES[$i]}" "$f" "$cand_dir" "$fasta"
				print_cmd
				expected=$(expected_output "$(perl_format "$f")")
				compare_cmd "$f" "$expected" "$cand_dir/$(format_file "$f")" "$OUT_DIR/$pin/$C_ID/$f.json"
				print_cmd
			done
			if [[ "$C_WARNINGS" == "true" ]]; then
				echo "DRY-RUN $C_ID $pin warnings"
				if [[ -n "$REFERENCE_DIR" ]]; then
					expected=$(stored_warnings_file "$C_ID")
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

# The digest of every image used, pulled when absent; recorded per report.
IMAGE_DIGESTS=""
if [[ -z "$REFERENCE_DIR" ]]; then
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
fi

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
	if has_flag --nearest && [[ "$C_PERL_MODE" == "offline" && ! -e "$C_PERL_SLICE/transcript_gene_tss.txt" ]] &&
		command -v flock >/dev/null 2>&1; then
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
        record["formats"][fmt] = entry
        reasons.extend(f"{fmt}: {r}" for r in sub_reasons)
    elif key == "format_error":
        fmt, _, message = value.partition("=")
        record["formats"][fmt] = {"verdict": "red", "error": message}
        reasons.append(f"{fmt}: {message}")
    elif key == "flags":
        record[key] = value.split()
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
    elif key in ("perl_image_digest", "perl_server", "reference_dir", "fasta", "manifest", "candidate_version"):
        record[key] = value or None
    else:
        record[key] = value
# A stored reference carries its image in its provenance; a live run records the digest
# docker reports.
if not record.get("perl_image_digest") and record.get("reference_provenance_image"):
    record["perl_image_digest"] = record["reference_provenance_image"]
# A floor a format could not count is met only where a text format of the same case
# counted it; a case that compared no text format leaves it unchecked.
text_compared = any(
    fmt in ("default", "tab", "vcf", "parquet") and "report" in entry for fmt, entry in record["formats"].items()
)
if not text_compared:
    for fmt, entry in sorted(record["formats"].items()):
        for name in entry.get("floors_not_applicable", []):
            reasons.append(f"{fmt}: floor {name} is not applicable to {fmt} and no text format of the case checked it")
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

	local reference_provenance=""
	if [[ -n "$REFERENCE_DIR" ]]; then
		perl_dir="$REFERENCE_DIR/$C_ID/expected"
		[[ ! -f "$REFERENCE_DIR/$C_ID/provenance.json" ]] || reference_provenance="$REFERENCE_DIR/$C_ID/provenance.json"
		log "[$C_ID] Ensembl VEP: stored reference $perl_dir"
	else
		perl_dir="$case_work/perl"
		perl_image_digest=$(tsv_lookup "$IMAGE_DIGESTS" "$C_PERL_IMAGE")
		for f in $C_PERL_FORMATS; do
			perl_cmd "$f" "$perl_dir" "$fasta"
			log "[$C_ID] Ensembl VEP $f ($C_PERL_MODE${C_PERL_SERVER:+ $C_PERL_SERVER})"
			rc=0
			run_perl "$case_work/perl.$f.stderr" || rc=$?
			if [[ "$rc" -ne 0 ]]; then
				perl_errors="${perl_errors}${f}	$(error_text "Ensembl VEP" "$rc" "$case_work/perl.$f.stderr")
"
				echo "WARN: [feature_oracle] $C_ID: Ensembl VEP $f exited $rc" >&2
			elif [[ ! -e "$perl_dir/$(format_file "$f")" ]]; then
				perl_errors="${perl_errors}${f}	Ensembl VEP exited 0 and wrote no $(format_file "$f")
"
			fi
		done
	fi

	local i pin bin cand_dir report_dir pf perl_error expected actual report_args verdict digest
	for ((i = 0; i < N_CANDIDATES; i++)); do
		pin="${PINS[$i]}"
		bin="${CANDIDATES[$i]}"
		cand_dir="$case_work/candidates/$pin"
		report_dir="$OUT_DIR/$pin/$C_ID"
		rm -rf "$report_dir"
		mkdir -p "$cand_dir" "$report_dir"
		report_args=()
		for f in $C_OUTPUT_FORMATS; do
			candidate_cmd "$bin" "$f" "$cand_dir" "$fasta"
			log "[$C_ID] $pin $f"
			rc=0
			run_engine "$cand_dir/$f.stderr" "${CMD[@]}" || rc=$?
			if [[ "$rc" -ne 0 ]]; then
				report_args+=("format_error=$f=$(error_text "candidate" "$rc" "$cand_dir/$f.stderr")")
				continue
			fi
			pf=$(perl_format "$f")
			perl_error=$(tsv_lookup "$perl_errors" "$pf")
			if [[ -n "$perl_error" ]]; then
				report_args+=("format_error=$f=$perl_error")
				continue
			fi
			expected=$(expected_output "$pf")
			actual="$cand_dir/$(format_file "$f")"
			compare_cmd "$f" "$expected" "$actual" "$report_dir/$f.json"
			report_args+=("$(run_compare "$f" "$report_dir/$f.json" "$cand_dir/$f")")
		done
		if [[ "$C_WARNINGS" == "true" ]]; then
			perl_error=$(tsv_lookup "$perl_errors" default)
			if [[ -n "$perl_error" ]]; then
				report_args+=("format_error=warnings=$perl_error")
			elif [[ ! -e "$cand_dir/default.txt" ]]; then
				report_args+=("format_error=warnings=the candidate's default run did not complete")
			else
				if [[ -n "$REFERENCE_DIR" ]]; then
					expected=$(stored_warnings_file "$C_ID")
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
			"case=$C_ID" "tier=$C_TIER" "input=$C_INPUT_PATH" "cache=$C_CACHE" \
			"candidate_json_cache=$C_JSON_CACHE" "fasta=$fasta" "flags=$C_FLAGS" \
			"formats_requested=$C_FORMATS" "manifest=$C_MANIFEST" "masks=$C_MASKS" \
			"floors=$C_FLOORS" "field_pairs=$C_FIELD_PAIRS" \
			"perl_image=$C_PERL_IMAGE" "perl_image_digest=$perl_image_digest" "perl_mode=$C_PERL_MODE" \
			"perl_server=$C_PERL_SERVER" "reference_release=$C_REFERENCE_RELEASE" \
			"reference_dir=${REFERENCE_DIR:+$REFERENCE_DIR/$C_ID}" "reference_provenance=$reference_provenance" \
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

"$PYTHON" - "$OUT_DIR/run.json" "$RUN_UTC" "$TIER" "$CASES" "$REFERENCE_DIR" "$PERL_FORK" "$FORK" "$IMAGE_DIGESTS" "$N_ROWS" "$N_RED" "${PINS[@]}" <<'PY'
import json
import sys

(path, run_utc, tier, cases, reference_dir, perl_fork, fork, digests, rows, red) = sys.argv[1:11]
record = {
    "run_utc": run_utc,
    "tier": tier,
    "case_table": cases,
    "reference_dir": reference_dir or None,
    "perl_fork": int(perl_fork),
    "candidate_fork": int(fork) if fork else None,
    "image_digests": dict(line.split("\t", 1) for line in digests.splitlines() if line),
    "candidate_pins": sys.argv[11:],
    "status_rows": int(rows),
    "red": int(red),
}
with open(path, "w", encoding="utf-8") as fh:
    json.dump(record, fh, indent=1, sort_keys=True)
    fh.write("\n")
PY

log "feature_oracle complete: $N_ROWS verdict(s), $N_RED red; $STATUS"
[[ "$N_RED" -eq 0 ]] || exit 1
