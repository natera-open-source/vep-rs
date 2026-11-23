#!/usr/bin/env bash
# Tests of harness_genomes.sh over the committed harness_genomes.tsv and over tables
# written here: every column of the two committed rows, a genome and a column the table
# lacks, a later row that replaces an earlier one, and the tables the function refuses (a
# row of another width than the header, a header without the two key columns, an empty
# file, a path that is no file).
#
# A passing run prints `ok` lines alone; a failing check prints a `FAIL` line and the run
# exits 1.
#
# Run: bash scripts/concordance/test_harness_genomes.sh
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=scripts/concordance/harness_genomes.sh
source "$HERE/harness_genomes.sh"
TABLE="$HERE/harness_genomes.tsv"

PASS=0
FAIL=0
ok() {
	PASS=$((PASS + 1))
	echo "  ok   $1"
}
bad() {
	FAIL=$((FAIL + 1))
	echo "  FAIL $1"
}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# field_is <name> <expected> <table> <species> <assembly> <column>: the function prints
# exactly <expected> on standard output, nothing on standard error, and returns 0.
field_is() {
	local name="$1" expected="$2" got status
	shift 2
	got=$(harness_genome_field "$@" 2>"$WORK/stderr")
	status=$?
	if [[ $status -eq 0 && "$got" == "$expected" && ! -s "$WORK/stderr" ]]; then
		ok "$name"
	else
		bad "$name: status $status, printed '$got', stderr '$(cat "$WORK/stderr")'"
	fi
}

# refused <name> <expected stderr> <arguments...>: the function prints nothing on standard
# output, exactly <expected stderr> on standard error, and returns 1.
refused() {
	local name="$1" expected="$2" got status
	shift 2
	got=$(harness_genome_field "$@" 2>"$WORK/stderr")
	status=$?
	if [[ $status -eq 1 && -z "$got" && "$(cat "$WORK/stderr")" == "$expected" ]]; then
		ok "$name"
	else
		bad "$name: status $status, printed '$got', stderr '$(cat "$WORK/stderr")'"
	fi
}

field_is "grch37_species" homo_sapiens "$TABLE" homo_sapiens GRCh37 species
field_is "grch37_assembly" GRCh37 "$TABLE" homo_sapiens GRCh37 assembly
field_is "grch37_cache_version" 116 "$TABLE" homo_sapiens GRCh37 cache_version
field_is "grch37_fasta" reference/grch37/genome.fa "$TABLE" homo_sapiens GRCh37 fasta
field_is "grch37_json_cache" caches/vep-rs/grch37 "$TABLE" homo_sapiens GRCh37 json_cache
field_is "grch38_species" homo_sapiens "$TABLE" homo_sapiens GRCh38 species
field_is "grch38_assembly" GRCh38 "$TABLE" homo_sapiens GRCh38 assembly
field_is "grch38_cache_version" 116 "$TABLE" homo_sapiens GRCh38 cache_version
field_is "grch38_fasta" reference/grch38/genome.fa "$TABLE" homo_sapiens GRCh38 fasta
field_is "grch38_json_cache" caches/vep-rs/grch38 "$TABLE" homo_sapiens GRCh38 json_cache

refused "a_genome_the_table_lacks" \
	"ERROR: [harness_genomes] $TABLE: no row for the genome mus_musculus GRCm39" \
	"$TABLE" mus_musculus GRCm39 fasta
refused "an_assembly_the_table_lacks" \
	"ERROR: [harness_genomes] $TABLE: no row for the genome homo_sapiens GRCm39" \
	"$TABLE" homo_sapiens GRCm39 fasta
refused "a_column_the_table_lacks" \
	"ERROR: [harness_genomes] $TABLE: no column gff3; the columns are species, assembly, cache_version, fasta, json_cache" \
	"$TABLE" homo_sapiens GRCh38 gff3

# A later row for a genome replaces an earlier one, and leaves the other genomes alone.
LATER="$WORK/later.tsv"
cat "$TABLE" >"$LATER"
printf 'mus_musculus\tGRCm39\t116\treference/grcm39/genome.fa\tcaches/vep-rs/grcm39\n' >>"$LATER"
printf 'homo_sapiens\tGRCh38\t117\tother/genome.fa\tother/cache\n' >>"$LATER"
field_is "a_later_row_replaces_the_fasta" other/genome.fa "$LATER" homo_sapiens GRCh38 fasta
field_is "a_later_row_replaces_the_cache_version" 117 "$LATER" homo_sapiens GRCh38 cache_version
field_is "a_later_row_leaves_another_assembly_alone" reference/grch37/genome.fa "$LATER" homo_sapiens GRCh37 fasta
field_is "an_added_genome_is_read" caches/vep-rs/grcm39 "$LATER" mus_musculus GRCm39 json_cache

# The columns are found by name, wherever the header puts them.
MOVED="$WORK/moved.tsv"
printf 'json_cache\tassembly\tspecies\nsome/cache\tR64-1-1\tsaccharomyces_cerevisiae\n' >"$MOVED"
field_is "columns_are_found_by_name" some/cache "$MOVED" saccharomyces_cerevisiae R64-1-1 json_cache

# An empty field is a field: the row has the header's width.
BLANK="$WORK/blank.tsv"
printf 'species\tassembly\tcache_version\tfasta\tjson_cache\nhomo_sapiens\tGRCh38\t116\t\tcaches/vep-rs/grch38\n' >"$BLANK"
field_is "an_empty_field_prints_an_empty_line" "" "$BLANK" homo_sapiens GRCh38 fasta

# A name that looks like a number is compared as text: 1.0 is neither 1.00 nor 1.
NUMERIC="$WORK/numeric.tsv"
printf 'species\tassembly\tfasta\nexample\t1.0\tone.fa\nexample\t1.00\ttwo.fa\n' >"$NUMERIC"
field_is "a_numeric_assembly_is_compared_as_text" one.fa "$NUMERIC" example 1.0 fasta
refused "a_numeric_assembly_spelled_another_way_is_no_match" \
	"ERROR: [harness_genomes] $NUMERIC: no row for the genome example 1" \
	"$NUMERIC" example 1 fasta

# A short row is refused wherever it stands, for every genome.
SHORT="$WORK/short.tsv"
cat "$TABLE" >"$SHORT"
printf 'mus_musculus\tGRCm39\t116\treference/grcm39/genome.fa\n' >>"$SHORT"
refused "a_short_row_is_refused_for_its_own_genome" \
	"ERROR: [harness_genomes] $SHORT: line 4 has 4 fields where the header has 5" \
	"$SHORT" mus_musculus GRCm39 fasta
refused "a_short_row_is_refused_for_another_genome" \
	"ERROR: [harness_genomes] $SHORT: line 4 has 4 fields where the header has 5" \
	"$SHORT" homo_sapiens GRCh37 fasta

LONG="$WORK/long.tsv"
cat "$TABLE" >"$LONG"
printf 'mus_musculus\tGRCm39\t116\ta\tb\tc\n' >>"$LONG"
refused "a_long_row_is_refused" \
	"ERROR: [harness_genomes] $LONG: line 4 has 6 fields where the header has 5" \
	"$LONG" homo_sapiens GRCh37 fasta

NOKEY="$WORK/nokey.tsv"
printf 'genome\tfasta\nhomo_sapiens\tgenome.fa\n' >"$NOKEY"
refused "a_header_without_the_key_columns_is_refused" \
	"ERROR: [harness_genomes] $NOKEY: the header line names no species and assembly columns" \
	"$NOKEY" homo_sapiens GRCh38 fasta

EMPTY="$WORK/empty.tsv"
: >"$EMPTY"
refused "an_empty_table_is_refused" \
	"ERROR: [harness_genomes] $EMPTY: the table is empty" \
	"$EMPTY" homo_sapiens GRCh38 fasta

refused "a_missing_table_is_refused" \
	"ERROR: [harness_genomes] $WORK/absent.tsv: the table cannot be read" \
	"$WORK/absent.tsv" homo_sapiens GRCh38 fasta
refused "a_directory_is_refused" \
	"ERROR: [harness_genomes] $WORK: the table cannot be read" \
	"$WORK" homo_sapiens GRCh38 fasta
refused "three_arguments_are_refused" \
	"ERROR: [harness_genomes] harness_genome_field takes <table> <species> <assembly> <column>, not 3 arguments" \
	"$TABLE" homo_sapiens GRCh38

# A name with a backslash reaches the lookup as it is written.
refused "a_backslash_in_a_name_is_kept" \
	"ERROR: [harness_genomes] $TABLE: no row for the genome homo\\tsapiens GRCh38" \
	"$TABLE" 'homo\tsapiens' GRCh38 fasta

# The error line keeps the backslash where the shell's echo expands escapes.
shopt -s xpg_echo
refused "a_backslash_in_a_name_is_kept_under_xpg_echo" \
	"ERROR: [harness_genomes] $TABLE: no row for the genome homo\\tsapiens GRCh38" \
	"$TABLE" 'homo\tsapiens' GRCh38 fasta
shopt -u xpg_echo

# A relative path that awk would read as an assignment is still the table.
cp "$TABLE" "$WORK/genomes=2.tsv"
cd "$WORK" || exit 1
field_is "a_relative_path_with_an_equals_sign_is_the_table" 116 "genomes=2.tsv" homo_sapiens GRCh38 cache_version </dev/null
cd - >/dev/null || exit 1

if [[ $FAIL -gt 0 ]]; then
	echo "ERROR: [test_harness_genomes] $FAIL of $((PASS + FAIL)) checks failed" >&2
	exit 1
fi
