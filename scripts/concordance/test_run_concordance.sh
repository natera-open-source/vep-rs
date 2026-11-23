#!/usr/bin/env bash
# Tests of run_concordance.sh on the script's own text. Each block under test is extracted
# by its first and last line and evaluated in a subshell that holds the variables the
# script holds at that point, so no engine, cache build or comparator runs.
#
# The cases:
#   CONTIGS   the canonical-contigs flag: bare under GRCh37 and GRCh38, the contig_order of
#             the JSON cache's info.json joined by commas under any other assembly, a
#             refusal where that cache lists none, and --no-canonical-contigs under any
#             assembly when the filter is off. The argument lists of the input filter and
#             of the tuple comparator both carry it.
#   FASTA     the two messages that name a reference FASTA, for a missing --json-cache-dir
#             and a missing --fasta: the human file under homo_sapiens, the species' own
#             under any other, where the cache builder's example also names the species
#             and leaves out the human mitochondrial codon table.
#
# A passing run prints `ok` lines alone; a failing check prints a `FAIL` line and the run
# exits 1.
#
# Run: bash scripts/concordance/test_run_concordance.sh

set -uo pipefail

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

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARNESS="$HERE/run_concordance.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# harness_block <first line> <last line>: the harness's lines from the first that equals
# <first line> to the next that equals <last line>, both included.
harness_block() {
    FIRST="$1" LAST="$2" awk '
        !p && $0 == ENVIRON["FIRST"] { p = 1 }
        p { print }
        p && $0 == ENVIRON["LAST"] { exit }
    ' "$HARNESS"
}

# extracted <what> <block> <last line>: the block holds lines and ends on its last line.
extracted() {
    if [[ -n "$2" && "$2" == *$'\n'"$3" ]]; then
        ok "run_concordance.sh: $1 is extracted whole"
    else
        bad "run_concordance.sh: could not extract $1"
    fi
}

# has_lines <text> <lines>: <lines> stand in <text> as whole consecutive lines.
has_lines() {
    [[ $'\n'"$1"$'\n' == *$'\n'"$2"$'\n'* ]]
}

# line_of <line>: the number of the harness's first line that equals <line>.
line_of() {
    LINE="$1" awk '$0 == ENVIRON["LINE"] { print NR; exit }' "$HARNESS"
}

# The harness's own lines, matched verbatim.
# shellcheck disable=SC2016
{
    FASTA_GATE_FIRST='if [[ -z "${FASTA}" && "${NO_FASTA}" -ne 1 ]]; then'
    CACHE_GATE_FIRST='if [[ -z "${JSON_CACHE_DIR}" ]]; then'
    CACHE_ARGS_FIRST='COMPARE_CACHE_ARGS=(--json-cache "${JSON_CACHE_DIR}")'
}
_flag_block=$(harness_block 'CANONICAL_CONTIGS_FLAG=(--no-canonical-contigs)' 'fi')
_prep_args=$(harness_block 'PREP_ARGS=(' ')')
_cache_args=$(harness_block "$CACHE_ARGS_FIRST" 'fi')
_compare_args=$(harness_block 'COMPARE_ARGS=(' ')')
_fasta_name=$(harness_block 'ensembl_fasta_name() {' '}')
_fasta_gate=$(harness_block "$FASTA_GATE_FIRST" 'fi')
_cache_gate=$(harness_block "$CACHE_GATE_FIRST" 'fi')

echo "CONTIGS: the canonical-contigs flag follows the assembly, the cache and the filter switch"

extracted "the canonical-contigs flag block" "$_flag_block" 'fi'
extracted "the input filter's argument list" "$_prep_args" ')'
extracted "the tuple comparator's argument list" "$_compare_args" ')'

# Three fabricated caches of a made-up genome: one lists its contigs out of sort order, one
# has no contig_order key, one lists none.
LISTED="$TMP/cache_listed"
UNLISTED="$TMP/cache_unlisted"
EMPTY_LIST="$TMP/cache_empty_list"
mkdir -p "$LISTED" "$UNLISTED" "$EMPTY_LIST"
printf '%s\n' '{"species": "madeup_species", "assembly": "MadeUp1", "contig_order": ["2L", "X", "mitochondrion_genome", "4"]}' >"$LISTED/info.json"
printf '%s\n' '{"species": "madeup_species", "assembly": "MadeUp1"}' >"$UNLISTED/info.json"
printf '%s\n' '{"species": "madeup_species", "assembly": "MadeUp1", "contig_order": []}' >"$EMPTY_LIST/info.json"
CONTIG_ORDER="2L,X,mitochondrion_genome,4"

# contigs_flag <assembly> <filter: 1 or 0> <cache dir>: runs the flag block under the
# harness's shell options; prints its standard error, the flag's words one a line, then `rc=`.
contigs_flag() {
    # The extracted block reads the names this subshell sets.
    # shellcheck disable=SC2034
    (
        set -euo pipefail
        ASSEMBLY="$1" CANONICAL_CONTIGS="$2" JSON_CACHE_DIR="$3"
        eval "$_flag_block"
        printf '%s\n' "${CANONICAL_CONTIGS_FLAG[@]}"
    ) 2>&1
    echo "rc=$?"
}
_no_contig_order() {
    printf '%s' "ERROR: [run_concordance] $1/info.json lists no contig_order, so the canonical contigs of assembly MadeUp1 are unknown; pass --no-canonical-contigs to annotate every contig"
}

_out=$(contigs_flag MadeUp1 1 "$LISTED")
if [[ "$_out" == "--canonical-contigs"$'\n'"$CONTIG_ORDER"$'\n'"rc=0" ]]; then
    ok "an assembly that is neither human one, filter on: the flag's value is the cache's contig_order, in the cache's order, joined by commas"
else
    bad "MadeUp1 with the filter on: expected --canonical-contigs and $CONTIG_ORDER, got '$_out'"
fi
for _cache in "$UNLISTED" "$EMPTY_LIST"; do
    _out=$(contigs_flag MadeUp1 1 "$_cache")
    if [[ "$_out" == "$(_no_contig_order "$_cache")"$'\n'"rc=1" ]]; then
        ok "an assembly that is neither human one, ${_cache##*/cache_}: refused with one line that names the cache and the opt-out, exit 1"
    else
        bad "MadeUp1 on ${_cache##*/}: expected the refusal and rc=1, got '$_out'"
    fi
done
for _assembly in GRCh37 GRCh38; do
    _out=$(contigs_flag "$_assembly" 1 "$LISTED")
    if [[ "$_out" == "--canonical-contigs"$'\n'"rc=0" ]]; then
        ok "$_assembly, filter on: the flag is bare and the cache's contig_order is not its value"
    else
        bad "$_assembly with the filter on: expected a bare --canonical-contigs, got '$_out'"
    fi
done
for _case in "GRCh37 $LISTED" "GRCh38 $LISTED" "MadeUp1 $LISTED" "MadeUp1 $UNLISTED"; do
    _out=$(contigs_flag "${_case%% *}" 0 "${_case#* }")
    if [[ "$_out" == "--no-canonical-contigs"$'\n'"rc=0" ]]; then
        ok "${_case%% *}, filter off, ${_case##*/cache_}: the flag is --no-canonical-contigs and no cache is refused"
    else
        bad "${_case%% *} with the filter off on ${_case##*/}: expected --no-canonical-contigs, got '$_out'"
    fi
done

# Evaluates the flag block, then both argument lists, and writes each list's words one a
# line to a file of $TMP.
# shellcheck disable=SC2034
(
    set -euo pipefail
    ASSEMBLY=MadeUp1 CANONICAL_CONTIGS=1 JSON_CACHE_DIR="$LISTED"
    BENCHMARK_DIR=/bench WORK_DIR=/work REFERENCE_RELEASE=116.2
    eval "$_flag_block"
    eval "$_prep_args"
    eval "$_cache_args"
    eval "$_compare_args"
    printf '%s\n' "${PREP_ARGS[@]}" >"$TMP/prepare.words"
    printf '%s\n' "${COMPARE_ARGS[@]}" >"$TMP/compare.words"
) 2>"$TMP/lists.stderr" || echo "the argument lists did not evaluate" >>"$TMP/lists.stderr"
for _list in prepare compare; do
    if [[ ! -s "$TMP/lists.stderr" ]] && has_lines "$(cat "$TMP/$_list.words" 2>/dev/null)" "--canonical-contigs"$'\n'"$CONTIG_ORDER"; then
        ok "the $_list step's arguments carry the flag and the cache's contig_order as its value"
    else
        bad "the $_list step's arguments lack --canonical-contigs $CONTIG_ORDER: '$(tr '\n' ' ' <"$TMP/$_list.words" 2>/dev/null)' $(cat "$TMP/lists.stderr")"
    fi
done
if [[ ! -s "$TMP/lists.stderr" ]] && has_lines "$(cat "$TMP/compare.words" 2>/dev/null)" "--json-cache"$'\n'"$LISTED"; then
    ok "the compare step's arguments name the cache the vep-rs run reads as --json-cache"
else
    bad "the compare step's arguments lack --json-cache $LISTED: '$(tr '\n' ' ' <"$TMP/compare.words" 2>/dev/null)' $(cat "$TMP/lists.stderr")"
fi

extracted "the block that hands the comparator the cache" "$_cache_args" 'fi'

# compare_words <assembly> <cache dir>: evaluates, with the filter off and under the harness's
# shell options, the flag block, the block that hands the comparator the cache and the
# comparator's argument list; prints their standard error, the list's words on one line,
# then `rc=`.
compare_words() {
    # The extracted blocks read the names this subshell sets.
    # shellcheck disable=SC2034
    (
        set -euo pipefail
        ASSEMBLY="$1" CANONICAL_CONTIGS=0 JSON_CACHE_DIR="$2"
        WORK_DIR=/work REFERENCE_RELEASE=116.2
        eval "$_flag_block"
        eval "$_cache_args"
        eval "$_compare_args"
        echo "${COMPARE_ARGS[*]}"
    ) 2>&1
    echo "rc=$?"
}
# The comparator ends a run under an assembly that is neither human one on a cache that lists
# no contig_order, whatever the filter switch.
for _cache in "$UNLISTED" "$EMPTY_LIST"; do
    _out=$(compare_words MadeUp1 "$_cache")
    if [[ "$_out" == *" --no-canonical-contigs"$'\n'"rc=0" && "$_out" != *--json-cache* ]]; then
        ok "an assembly that is neither human one, filter off, ${_cache##*/cache_}: the compare step's arguments name no cache"
    else
        bad "MadeUp1 with the filter off on ${_cache##*/}: expected arguments without --json-cache, got '$_out'"
    fi
done
for _case in "MadeUp1 $LISTED" "GRCh37 $UNLISTED" "GRCh38 $UNLISTED"; do
    _out=$(compare_words "${_case%% *}" "${_case#* }")
    if [[ "$_out" == *" --no-canonical-contigs --json-cache ${_case#* }"$'\n'"rc=0" ]]; then
        ok "${_case%% *}, filter off, ${_case##*/cache_}: the compare step's arguments name the cache"
    else
        bad "${_case%% *} with the filter off on ${_case##*/}: expected arguments ending --json-cache ${_case#* }, got '$_out'"
    fi
done

echo "FASTA: a message that names a reference FASTA names the species' own"

extracted "ensembl_fasta_name" "$_fasta_name" '}'
extracted "the --fasta gate" "$_fasta_gate" 'fi'
extracted "the --json-cache-dir gate" "$_cache_gate" 'fi'

# A function defined below its first call, or an option parsed below the message that reads
# it, leaves that message with an empty word and no error.
_species_parsed=$(line_of '    --species)')
_name_defined=$(line_of 'ensembl_fasta_name() {')
_fasta_gate_line=$(line_of "$FASTA_GATE_FIRST")
_cache_gate_line=$(line_of "$CACHE_GATE_FIRST")
if [[ -n "$_species_parsed" && -n "$_name_defined" && -n "$_fasta_gate_line" && -n "$_cache_gate_line" && "$_species_parsed" -lt "$_fasta_gate_line" && "$_name_defined" -lt "$_fasta_gate_line" && "$_fasta_gate_line" -lt "$_cache_gate_line" ]]; then
    ok "--species is parsed and ensembl_fasta_name defined before the two messages that read them (lines $_species_parsed and $_name_defined before $_fasta_gate_line and $_cache_gate_line)"
else
    bad "--species (line ${_species_parsed:-?}) and ensembl_fasta_name (line ${_name_defined:-?}) must stand before the --fasta gate (line ${_fasta_gate_line:-?}) and the --json-cache-dir gate (line ${_cache_gate_line:-?})"
fi

# Ensembl's own file names for the mouse and fly genomes.
for _case in "mus_musculus GRCm39 Mus_musculus.GRCm39.dna.toplevel.fa" "drosophila_melanogaster BDGP6.54 Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa"; do
    read -r _species _assembly _want <<<"$_case"
    _out=$(
        eval "$_fasta_name"
        ensembl_fasta_name "$_species" "$_assembly" 2>&1
    )
    if [[ "$_out" == "$_want" ]]; then
        ok "ensembl_fasta_name $_species $_assembly: $_want"
    else
        bad "ensembl_fasta_name $_species $_assembly: expected $_want, got '$_out'"
    fi
done

# gate <block> <species> <assembly>: runs a gate with neither a FASTA nor a JSON cache given,
# under the harness's shell options; prints the gate's message, then `rc=`.
gate() {
    # The extracted block reads the names this subshell sets.
    # shellcheck disable=SC2034
    (
        set -euo pipefail
        SPECIES="$2" ASSEMBLY="$3" WORK_DIR=/work FASTA="" NO_FASTA=0 JSON_CACHE_DIR=""
        ASSEMBLY_LC="$(echo "$3" | tr '[:upper:]' '[:lower:]')"
        eval "$_fasta_name"
        eval "$1"
    ) 2>&1
    echo "rc=$?"
}

cat >"$TMP/human_cache_message" <<'EOF'
ERROR: --json-cache-dir is required; no JSON cache was given (for example /work/cache_json).
  Build one with vep-cache-builder (no Perl required):
    vep-cache-builder --release 116 --assembly GRCh38 \
      --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
      --core-dumps core_dumps_grch38 \
      --include-predictions --output-dir /work/cache_json
  --core-dumps names a directory of the release's table dumps (a missing table is
  downloaded into it): they supply what the GFF3 lacks and each region's codon table.
  A build without them needs --codon-table MT=2 where the annotation holds MT.
  Or provide an existing cache with --json-cache-dir.
rc=1
EOF
cat >"$TMP/other_cache_message" <<'EOF'
ERROR: --json-cache-dir is required; no JSON cache was given (for example /work/cache_json).
  Build one with vep-cache-builder (no Perl required):
    vep-cache-builder --release 116 --species madeup_species --assembly MadeUp1 \
      --genome-fasta Madeup_species.MadeUp1.dna.toplevel.fa \
      --core-dumps core_dumps_madeup1 \
      --include-predictions --output-dir /work/cache_json
  --core-dumps names a directory of the release's table dumps (a missing table is
  downloaded into it): they supply what the GFF3 lacks and each region's codon table.
  Or provide an existing cache with --json-cache-dir.
rc=1
EOF

_out=$(gate "$_cache_gate" homo_sapiens GRCh38)
if [[ "$_out" == "$(cat "$TMP/human_cache_message")" ]]; then
    ok "no JSON cache, homo_sapiens on GRCh38: the example names the primary-assembly FASTA and no species, the MT=2 sentence stands, exit 1"
else
    bad "no JSON cache under homo_sapiens on GRCh38: the message differs from the expected one: '$_out'"
fi
_out=$(gate "$_cache_gate" homo_sapiens GRCh37)
if [[ "$_out" == *$'\nrc=1' ]] && has_lines "$_out" "    vep-cache-builder --release 116 --assembly GRCh37 \\"$'\n'"      --genome-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \\" && has_lines "$_out" "  A build without them needs --codon-table MT=2 where the annotation holds MT."; then
    ok "no JSON cache, homo_sapiens on GRCh37: the example names that assembly's primary-assembly FASTA, the MT=2 sentence stands, exit 1"
else
    bad "no JSON cache under homo_sapiens on GRCh37: expected the GRCh37 example, the MT=2 sentence and rc=1, got '$_out'"
fi
_out=$(gate "$_cache_gate" madeup_species MadeUp1)
if [[ "$_out" == "$(cat "$TMP/other_cache_message")" && "$_out" != *MT=2* && "$_out" != *Homo_sapiens* ]]; then
    ok "no JSON cache, a made-up species: the example names --species and the species' toplevel FASTA, no MT=2 sentence, exit 1"
else
    bad "no JSON cache under a made-up species: the message differs from the expected one: '$_out'"
fi

_out=$(gate "$_fasta_gate" homo_sapiens GRCh37)
if [[ "$_out" == "ERROR: --fasta is required"*$'\nrc=1' ]] && has_lines "$_out" "  Fix: pass --fasta /path/to/Homo_sapiens.GRCh3{7,8}.dna.primary_assembly.fa"; then
    ok "no FASTA, homo_sapiens: the message names the human primary-assembly FASTA of either assembly, exit 1"
else
    bad "no FASTA under homo_sapiens: expected the human file name and rc=1, got '$_out'"
fi
_out=$(gate "$_fasta_gate" madeup_species MadeUp1)
if [[ "$_out" == "ERROR: --fasta is required"*$'\nrc=1' && "$_out" != *Homo_sapiens* ]] && has_lines "$_out" "  Fix: pass --fasta /path/to/Madeup_species.MadeUp1.dna.toplevel.fa"; then
    ok "no FASTA, a made-up species: the message names the species' toplevel FASTA and no human file, exit 1"
else
    bad "no FASTA under a made-up species: expected the species' file name and rc=1, got '$_out'"
fi

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
