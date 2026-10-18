#!/usr/bin/env bash
# Tests for generate_reference_output.sh: the dry-run contract and every fail-closed
# path, plus the output-volume gate on the real-run path via a stub `docker`.
#
# Pure bash with a fabricated data tree. No Docker, no VEP, no cache. Every location the
# script reads or writes is passed explicitly, so the tests do not depend on either
# tree's default layout.
#
# Run: bash scripts/concordance/test_generate_reference_output.sh

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

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GEN="$HERE/generate_reference_output.sh"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [[ ! -f "$GEN" ]]; then
    bad "generate_reference_output.sh not found beside this test"
    echo "passed: $PASS   failed: $FAIL"
    exit 1
fi

bash -n "$GEN" && ok "generate_reference_output.sh parses" || bad "generate_reference_output.sh has a syntax error"

# A fabricated data tree. Three records per SNP/indel input, two SV inputs per
# assembly plus a `manifest.vcf` the script must skip, a full 25-contig cache slice
# per assembly, and an indexed FASTA per assembly.
DATA="$TMP/data"
CACHE="$DATA/cache"
SNP_GT="$DATA/gt/snp_indel"
SV_GT="$DATA/gt/sv"
SV_IN="$SV_GT/canonical_inputs"
PINNED_DIGEST="sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da"

write_vcf() { # write_vcf <path> <records>
    local path="$1" n="$2" i
    {
        echo "##fileformat=VCFv4.2"
        echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
        for ((i = 1; i <= n; i++)); do
            echo "21	$((10000000 + i))	.	A	G	.	PASS	."
        done
    } >"$path"
}

build_tree() {
    rm -rf "$DATA"
    local asm asm_lc c base
    for asm in GRCh37 GRCh38; do
        asm_lc=$(printf %s "$asm" | tr '[:upper:]' '[:lower:]')
        for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
            mkdir -p "$CACHE/homo_sapiens/116_$asm/$c"
        done
        mkdir -p "$DATA/reference/$asm_lc" "$DATA/vcf/$asm_lc/all_variants_canonical" "$SV_IN/$asm_lc"
        printf '>21\nACGT\n' >"$DATA/reference/$asm_lc/genome.fa"
        printf '21\t4\t4\t4\t5\n' >"$DATA/reference/$asm_lc/genome.fa.fai"
        for base in clinvar_$asm_lc gnomad_genomes_v2.1.1_chr21 gnomad_genomes_v4.1_chr21 1kg_integrated_chr21_sitesonly 1kg_highcov_chr21; do
            write_vcf "$TMP/plain.vcf" 3
            gzip -c "$TMP/plain.vcf" >"$DATA/vcf/$asm_lc/all_variants_canonical/$base.canonical.vcf.gz"
        done
        write_vcf "$SV_IN/$asm_lc/01_snv.vcf" 4
        write_vcf "$TMP/plain.vcf" 2
        gzip -c "$TMP/plain.vcf" >"$SV_IN/$asm_lc/09_breakends.vcf.gz"
        write_vcf "$SV_IN/$asm_lc/manifest.vcf" 1
    done
}

# run_gen <extra args...>: dry-run against the fabricated tree with every location explicit.
run_gen() {
    bash "$GEN" --data-dir "$DATA" \
        --perl-cache-dir "$CACHE" \
        --fasta-grch37 "$DATA/reference/grch37/genome.fa" \
        --fasta-grch38 "$DATA/reference/grch38/genome.fa" \
        --snp-indel-gt-dir "$SNP_GT" --sv-gt-dir "$SV_GT" \
        "$@"
}

count_docker_cmds() { # count_docker_cmds <output-file>
    local n
    n=$(grep -c '^docker run ' "$1") || n=0
    echo "$n"
}

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'; else shasum -a 256 "$1" | awk '{print $1}'; fi
}

echo "DRY-RUN: a complete tree prints one docker command per output and writes nothing"
build_tree
run_gen --dry-run >"$TMP/dry.out" 2>"$TMP/dry.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "dry run over a complete tree exits 0" || {
    bad "dry run exited $rc"
    head -5 "$TMP/dry.err" | sed 's/^/         /'
}
n=$(count_docker_cmds "$TMP/dry.out")
# 6 SNP/indel suites + 2 SV files per assembly x 2 assemblies (manifest skipped) = 10.
[[ "$n" -eq 10 ]] && ok "10 docker commands printed (6 SNP/indel + 2x2 SV, manifest skipped)" || bad "$n docker commands printed, expected 10"
grep -q 'manifest' "$TMP/dry.out" && bad "the manifest file was planned as an SV input" || ok "manifest.vcf is skipped"
[[ ! -e "$SNP_GT" && ! -e "$SV_GT/grch37" ]] && ok "dry run wrote no output directory" || bad "dry run created an output directory"

n=$(grep -c -- '--dir_cache /work/cache' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "every command passes --dir_cache /work/cache" || bad "$n of 10 commands pass --dir_cache"
grep -q -- '--full_cache_dir' "$TMP/dry.out" && bad "a command passes --full_cache_dir, which VEP ignores" || ok "no command passes --full_cache_dir"
n=$(grep -c -- '--fasta /work/ref.fa' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "every command passes --fasta" || bad "$n of 10 commands pass --fasta"
n=$(grep -c -- '--buffer_size 5000' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "every command passes --buffer_size 5000" || bad "$n of 10 commands pass --buffer_size 5000"
n=$(grep -c -- '--cache_version 116' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "every command passes --cache_version 116" || bad "$n of 10 commands pass --cache_version 116"
n=$(grep -c -- '--no_headers --no_stats --quiet --force_overwrite' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "every command carries the harness's output flags" || bad "$n of 10 commands carry --no_headers --no_stats --quiet --force_overwrite"
n=$(grep -cE -- '--force_overwrite ?$' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "the default format passes no output-format flag (Ensembl VEP's output_format default is vep)" || bad "$n of 10 default-format commands end at --force_overwrite"

echo
echo "FORK-IDENTITY: SNP/indel commands run --fork 16 and SV commands --fork 4"
# VEP's transcript selection above --max_sv_size is batch-dependent, so the SV fork count
# is part of the ground truth's identity: a reference generated at another fork count is
# a different reference.
p1_fork16=$(grep -E '^docker run .*/output\.txt ' "$TMP/dry.out" | grep -c -- '--fork 16') || p1_fork16=0
[[ "$p1_fork16" -eq 6 ]] && ok "all 6 SNP/indel commands pass --fork 16" || bad "$p1_fork16 of 6 SNP/indel commands pass --fork 16"
sv_fork4=$(grep -E '^docker run .*/(01_snv|09_breakends)\.txt ' "$TMP/dry.out" | grep -c -- '--fork 4 ') || sv_fork4=0
[[ "$sv_fork4" -eq 4 ]] && ok "all 4 SV commands pass --fork 4" || bad "$sv_fork4 of 4 SV commands pass --fork 4"
run_gen --dry-run --sv-fork 7 --fork 3 >"$TMP/dry2.out" 2>/dev/null
sv7=$(grep -E '^docker run .*/(01_snv|09_breakends)\.txt ' "$TMP/dry2.out" | grep -c -- '--fork 7 ') || sv7=0
p3=$(grep -E '^docker run .*/output\.txt ' "$TMP/dry2.out" | grep -c -- '--fork 3 ') || p3=0
[[ "$sv7" -eq 4 && "$p3" -eq 6 ]] && ok "--fork and --sv-fork override their own class only" || bad "--fork/--sv-fork overrides leaked across classes (sv=$sv7 snp-indel=$p3)"

echo
echo "LAYOUT: outputs land where the harness and the comparators read them"
grep -q "DRY-RUN \[snp_indel GRCh37 default\] -> $SNP_GT/s01/output.txt" "$TMP/dry.out" && ok "s01 -> <snp-indel-gt>/s01/output.txt" || bad "s01 output path is not <snp-indel-gt>/s01/output.txt"
grep -q "DRY-RUN \[sv GRCh37 default\] -> $SV_GT/grch37/01_snv.txt" "$TMP/dry.out" && ok "01_snv.vcf -> <sv-gt>/grch37/01_snv.txt" || bad "SV output path is not <sv-gt>/<asm>/<basename>.txt"
grep -q "DRY-RUN \[sv GRCh38 default\] -> $SV_GT/grch38/09_breakends.txt" "$TMP/dry.out" && ok "09_breakends.vcf.gz -> <sv-gt>/grch38/09_breakends.txt (both suffixes stripped)" || bad ".vcf.gz basename not stripped to <basename>.txt"
grep -q -- "-v $SV_IN/grch37/01_snv.vcf:/work/input.vcf:ro" "$TMP/dry.out" && ok "a plain .vcf is mounted as /work/input.vcf" || bad "plain .vcf input mount is wrong"
grep -q -- "-v $SV_IN/grch38/09_breakends.vcf.gz:/work/input.vcf.gz:ro" "$TMP/dry.out" && ok "a .vcf.gz is mounted as /work/input.vcf.gz" || bad ".vcf.gz input mount is wrong"
grep -q -- "-v $CACHE:/work/cache:ro" "$TMP/dry.out" && ok "the cache is mounted read-only at /work/cache" || bad "cache mount is wrong"
grep -q -- "-v $DATA/reference/grch38/genome.fa.fai:/work/ref.fa.fai:ro" "$TMP/dry.out" && ok "the FASTA index is mounted beside the FASTA" || bad "FASTA index mount is wrong"

echo
echo "FILTERS: --suites and --assembly select exactly what they name"
run_gen --dry-run --suites sv --assembly GRCh38 >"$TMP/f1.out" 2>/dev/null
n=$(count_docker_cmds "$TMP/f1.out")
[[ "$n" -eq 2 ]] && grep -q 'sv GRCh38' "$TMP/f1.out" && ! grep -q 'GRCh37' "$TMP/f1.out" &&
    ok "--suites sv --assembly GRCh38 plans only the GRCh38 SV set" || bad "--suites sv --assembly GRCh38 planned $n command(s)"
run_gen --dry-run --suites "s01 s07" >"$TMP/f2.out" 2>/dev/null
n=$(count_docker_cmds "$TMP/f2.out")
[[ "$n" -eq 3 ]] && grep -q '/s01/output.txt' "$TMP/f2.out" && ! grep -q '/s02/' "$TMP/f2.out" &&
    ok "--suites \"s01 s07\" plans s01 plus the GRCh37 SV set" || bad "--suites \"s01 s07\" planned $n command(s)"
run_gen --dry-run --suites s09 >"$TMP/f3.out" 2>/dev/null
n=$(count_docker_cmds "$TMP/f3.out")
[[ "$n" -eq 2 ]] && grep -q 'sv GRCh37' "$TMP/f3.out" && ok "--suites s09 resolves to the GRCh37 SV set it points at" || bad "--suites s09 planned $n command(s)"
run_gen --dry-run --suites snp-indel >"$TMP/f4.out" 2>/dev/null
n=$(count_docker_cmds "$TMP/f4.out")
[[ "$n" -eq 6 ]] && ok "--suites snp-indel plans the six SNP/indel suites" || bad "--suites snp-indel planned $n command(s)"
run_gen --dry-run --suites p1 >"$TMP/f5.out" 2>"$TMP/f5.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'no suite matched' "$TMP/f5.err" && ok "an unknown filter word is an error, not an empty success" || bad "--suites p1 exited $rc"
run_gen --dry-run --suites s99 >"$TMP/f6.out" 2>"$TMP/f6.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'no suite matched' "$TMP/f6.err" && ok "an unknown suite id is an error, not an empty success" || bad "--suites s99 exited $rc"

echo
echo "FORMATS: one container run per output format (Runner.pm:871 builds one OutputFactory per run)"
# Ensembl VEP's --tab, --vcf and --json each set output_format (Config.pm:486-514) and
# Runner.pm:871 builds the one OutputFactory the run writes through, so a reference in
# several formats is several runs.
run_gen --dry-run --format tab >"$TMP/fmt_tab.out" 2>"$TMP/fmt_tab.err"
rc=$?
n=$(count_docker_cmds "$TMP/fmt_tab.out")
[[ "$rc" -eq 0 && "$n" -eq 10 ]] && ok "--format tab plans one command per output (10)" || bad "--format tab: rc=$rc commands=$n"
n=$(grep -cE -- '--force_overwrite --tab ?$' "$TMP/fmt_tab.out") || n=0
[[ "$n" -eq 10 ]] && ok "every tab command ends with --tab after the harness's flags" || bad "$n of 10 tab commands end with --tab"
grep -q -- '--no_headers' "$TMP/fmt_tab.out" && bad "a tab command passes --no_headers, dropping the column line the comparators read" || ok "tab commands keep the header (the column line names the columns)"
grep -q "DRY-RUN \[snp_indel GRCh37 tab\] -> $SNP_GT/s01/output.tsv" "$TMP/fmt_tab.out" && ok "the tab output is <suite>/output.tsv" || bad "tab output name is not output.tsv"
grep -q "DRY-RUN \[sv GRCh37 tab\] -> $SV_GT/grch37/01_snv.tsv" "$TMP/fmt_tab.out" && ok "the SV tab output is <asm>/<basename>.tsv" || bad "SV tab output name is not <basename>.tsv"
run_gen --dry-run --format all --suites s01 >"$TMP/fmt_all.out" 2>/dev/null
rc=$?
n=$(count_docker_cmds "$TMP/fmt_all.out")
[[ "$rc" -eq 0 && "$n" -eq 4 ]] && ok "--format all on one suite plans four commands, one per format" || bad "--format all --suites s01: rc=$rc commands=$n"
d=$(grep -cE -- '/output\.txt --offline .* --no_headers --no_stats --quiet --force_overwrite ?$' "$TMP/fmt_all.out") || d=0
t=$(grep -cE -- '/output\.tsv --offline .* --no_stats --quiet --force_overwrite --tab ?$' "$TMP/fmt_all.out") || t=0
v=$(grep -cE -- '/output\.vcf --offline .* --no_stats --quiet --force_overwrite --vcf ?$' "$TMP/fmt_all.out") || v=0
j=$(grep -cE -- '/output\.jsonl --offline .* --no_stats --quiet --force_overwrite --json ?$' "$TMP/fmt_all.out") || j=0
[[ "$d" -eq 1 && "$t" -eq 1 && "$v" -eq 1 && "$j" -eq 1 ]] && ok "output.txt, output.tsv, output.vcf and output.jsonl each take their own format flag" || bad "format-to-name mapping: default=$d tab=$t vcf=$v json=$j"
n=$(grep -c -- '--no_headers' "$TMP/fmt_all.out") || n=0
[[ "$n" -eq 1 ]] && ok "--no_headers goes with the default format alone" || bad "$n of 4 commands pass --no_headers, expected 1"
run_gen --dry-run --format all >"$TMP/fmt_all10.out" 2>/dev/null
n=$(count_docker_cmds "$TMP/fmt_all10.out")
[[ "$n" -eq 40 ]] && ok "--format all over the whole tree plans 40 commands (10 outputs x 4 formats)" || bad "--format all planned $n commands, expected 40"
out=$(run_gen --dry-run --format parquet 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'invalid --format' <<<"$out" && ok "an unknown --format is refused" || bad "--format parquet exited $rc"

echo
echo "FLAGS: --flags words are appended verbatim after the harness's flags, before the format flag"
run_gen --dry-run --flags "--everything --max_sv_size -1" >"$TMP/flags.out" 2>"$TMP/flags.err"
rc=$?
n=$(grep -cE -- '--force_overwrite --everything --max_sv_size -1 ?$' "$TMP/flags.out") || n=0
[[ "$rc" -eq 0 && "$n" -eq 10 ]] && ok "every command carries the --flags words after --force_overwrite (10 of 10)" || bad "--flags: rc=$rc commands carrying the words=$n"
run_gen --dry-run --flags "--everything --max_sv_size -1" --format vcf --suites s01 >"$TMP/flags_vcf.out" 2>/dev/null
grep -qE -- '--force_overwrite --everything --max_sv_size -1 --vcf ?$' "$TMP/flags_vcf.out" && ok "the format flag follows the --flags words" || bad "the format flag does not follow the --flags words"
run_gen --dry-run --flags $'--everything\n--max_sv_size\t-1' --suites s01 >"$TMP/flags_ws.out" 2>/dev/null
grep -qE -- '--force_overwrite --everything --max_sv_size -1 ?$' "$TMP/flags_ws.out" && ok "a newline or tab between --flags words is one word boundary like a space" || bad "a newline or tab inside --flags lost or merged a word"
out=$(run_gen --dry-run --flags "--everything --fork 2" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--flags carries --fork' <<<"$out" && ok "a flag this script sets (--fork) is refused inside --flags" || bad "--flags with --fork exited $rc"
out=$(run_gen --dry-run --flags "--no_stats" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--flags carries --no_stats' <<<"$out" && ok "--no_stats inside --flags is refused (--stats owns it)" || bad "--flags with --no_stats exited $rc"
out=$(run_gen --dry-run --flags "--fasta=/x.fa" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--flags carries --fasta=/x.fa' <<<"$out" && ok "the --flag=value shape of an owned flag is refused too" || bad "--flags with --fasta=/x.fa exited $rc"
run_gen --dry-run --flags "--nearest symbol" --suites s01 >"$TMP/flags_near.out" 2>/dev/null
grep -q -- "-v $CACHE:/work/cache " "$TMP/flags_near.out" && ! grep -q -- "-v $CACHE:/work/cache:ro" "$TMP/flags_near.out" && ok "a --nearest family mounts the cache writable" || bad "--nearest family cache mount: $(grep -o -- "-v $CACHE:/work/cache[^ ]*" "$TMP/flags_near.out" | head -1)"
grep -q -- "-v $CACHE:/work/cache:ro" "$TMP/flags_vcf.out" && ok "a family without --nearest keeps the read-only cache mount" || bad "a family without --nearest lost the read-only mount"
out=$(run_gen --dry-run --flags "--vcf_info_field" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'no field name' <<<"$out" && ok "a trailing --vcf_info_field without its name is refused" || bad "--flags \"--vcf_info_field\" exited $rc"
# vep reads its options with Getopt::Long's defaults (vep:27, Config.pm:62 at release/116.2),
# which take a single dash or a plus sign, any unambiguous abbreviation and any letter case
# as the full option, and an option's alias (stats_file|sf) as the option.
out=$(run_gen --dry-run --flags "-fork 2" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- 'carries -fork, which vep reads as --fork' <<<"$out" && ok "the single-dash spelling of an owned flag (-fork) is refused" || bad "--flags \"-fork 2\" exited $rc"
out=$(run_gen --dry-run --flags "-tab" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- 'carries -tab, which vep reads as --tab' <<<"$out" && ok "-tab is refused (it would switch the default-format file to tab output)" || bad "--flags \"-tab\" exited $rc"
out=$(run_gen --dry-run --format tab --flags "--no_head" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- 'carries --no_head, which vep reads as --no_headers' <<<"$out" && ok "an abbreviation of an owned flag (--no_head) is refused" || bad "--flags \"--no_head\" exited $rc"
out=$(run_gen --dry-run --flags "--FORK 2" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- 'carries --FORK, which vep reads as --fork' <<<"$out" && ok "an upper-case spelling of an owned flag (--FORK) is refused" || bad "--flags \"--FORK 2\" exited $rc"
out=$(run_gen --dry-run --flags "+fork 2" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- 'carries +fork, which vep reads as --fork' <<<"$out" && ok "the plus-sign spelling of an owned flag (+fork) is refused" || bad "--flags \"+fork 2\" exited $rc"
out=$(run_gen --dry-run --flags "-sf x" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- 'carries -sf, which vep reads as --stats_file' <<<"$out" && ok "the sf alias for --stats_file is refused and named" || bad "--flags \"-sf x\" exited $rc"
run_gen --dry-run --flags "--max_sv_size -1 --individual all" --suites s01 >"$TMP/flags_ok.out" 2>/dev/null
rc=$?
[[ "$rc" -eq 0 ]] && grep -qE -- '--force_overwrite --max_sv_size -1 --individual all ?$' "$TMP/flags_ok.out" && ok "words naming no owned flag pass, a negative value word included" || bad "--flags with --max_sv_size -1 --individual all exited $rc"
# An exact option or alias name outranks an abbreviation in Getopt::Long: vep's `verbose|v`
# (Config.pm:159), `terms|t` (:176) and `dir` (:241) are prefixes of owned flags and pass.
run_gen --dry-run --flags "-v -t SO --dir /x" --suites s01 >"$TMP/flags_alias.out" 2>/dev/null
rc=$?
[[ "$rc" -eq 0 ]] && grep -qE -- '--force_overwrite -v -t SO --dir /x ?$' "$TMP/flags_alias.out" && ok "-v, -t and --dir, exact vep names that prefix owned flags, pass through" || bad "--flags \"-v -t SO --dir /x\" exited $rc"

echo
echo "STATS: --stats keeps the text statistics report instead of passing --no_stats"
# A run naming neither --stats_text nor --stats_html writes the HTML report alone
# (Runner.pm:787-789), so the text report takes --stats_text.
run_gen --dry-run --stats --suites s01 >"$TMP/stats.out" 2>/dev/null
grep -q -- '--no_stats' "$TMP/stats.out" && bad "--stats still passes --no_stats" || ok "--no_stats is absent under --stats"
grep -qE -- '--no_headers --stats_text --quiet --force_overwrite' "$TMP/stats.out" && ok "--stats passes --stats_text in --no_stats's place" || bad "--stats does not pass --stats_text"
grep -q -- '--stats_text' "$TMP/dry.out" && bad "--stats_text is passed without --stats" || ok "--stats_text is absent without --stats"

echo
echo "FAIL-CLOSED: a missing input aborts the whole run before any command is planned"
build_tree
rm "$DATA/vcf/grch38/all_variants_canonical/clinvar_grch38.canonical.vcf.gz"
run_gen --dry-run >"$TMP/m1.out" 2>"$TMP/m1.err"
rc=$?
[[ "$rc" -ne 0 ]] && ok "a missing SNP/indel input exits non-zero (rc=$rc)" || bad "a missing SNP/indel input exited 0"
grep -q 'clinvar_grch38.canonical.vcf.gz' "$TMP/m1.err" && ok "the missing input is named on stderr" || bad "stderr does not name the missing input"
grep -q 'preflight failed; nothing was generated' "$TMP/m1.err" && ok "the abort is reported as a preflight failure" || bad "no preflight-failure line"
n=$(count_docker_cmds "$TMP/m1.out")
[[ "$n" -eq 0 ]] && ok "no docker command is printed when preflight fails (the other 9 inputs were fine)" || bad "$n docker command(s) printed despite a missing input"

build_tree
rm -rf "$SV_IN/grch37"
run_gen --dry-run >/dev/null 2>"$TMP/m2.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'SV input directory not found' "$TMP/m2.err" && ok "a missing SV input directory exits non-zero" || bad "a missing SV input directory exited $rc"

build_tree
rm -f "$SV_IN/grch37"/*
write_vcf "$SV_IN/grch37/manifest.vcf" 1
run_gen --dry-run >/dev/null 2>"$TMP/m3.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'no .vcf or .vcf.gz under' "$TMP/m3.err" && ok "an SV directory holding only a manifest exits non-zero" || bad "a manifest-only SV directory exited $rc"

echo
echo "CACHE-COMPLETE: a chromosome-subset cache fails before anything runs"
build_tree
rm -rf "$CACHE/homo_sapiens/116_GRCh37/5" "$CACHE/homo_sapiens/116_GRCh37/X"
run_gen --dry-run >"$TMP/c1.out" 2>"$TMP/c1.err"
rc=$?
[[ "$rc" -ne 0 ]] && ok "a cache missing primary contigs exits non-zero" || bad "a chromosome-subset cache exited 0"
grep -q 'missing primary contigs: 5 X' "$TMP/c1.err" && ok "the missing contigs are named (5 X)" || bad "stderr does not name the missing contigs"
n=$(count_docker_cmds "$TMP/c1.out")
[[ "$n" -eq 0 ]] && ok "nothing is planned against an incomplete cache" || bad "$n command(s) planned against an incomplete cache"
run_gen --dry-run --assembly GRCh38 >"$TMP/c2.out" 2>/dev/null
rc=$?
n=$(count_docker_cmds "$TMP/c2.out")
[[ "$rc" -eq 0 && "$n" -eq 5 ]] && ok "the intact GRCh38 slice still generates when GRCh37 is excluded (5 commands)" || bad "--assembly GRCh38 with a broken GRCh37 slice: rc=$rc commands=$n"

build_tree
rm -rf "$CACHE/homo_sapiens/116_GRCh38"
run_gen --dry-run >/dev/null 2>"$TMP/c3.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'Perl cache slice not found' "$TMP/c3.err" && ok "an absent cache slice exits non-zero and says so" || bad "an absent cache slice exited $rc"

build_tree
run_gen --dry-run --cache-version 114 >/dev/null 2>"$TMP/c4.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q '114_GRCh37' "$TMP/c4.err" && ok "--cache-version selects the slice that is checked" || bad "--cache-version 114 against a 116 cache exited $rc"

echo
echo "FASTA: a missing FASTA or index fails the assembly"
build_tree
rm "$DATA/reference/grch37/genome.fa.fai"
run_gen --dry-run >/dev/null 2>"$TMP/fa1.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'FASTA index not found' "$TMP/fa1.err" && ok "a missing .fai exits non-zero" || bad "a missing .fai exited $rc"
build_tree
rm "$DATA/reference/grch38/genome.fa"
run_gen --dry-run >/dev/null 2>"$TMP/fa2.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'FASTA for GRCh38 not found' "$TMP/fa2.err" && ok "a missing FASTA exits non-zero" || bad "a missing FASTA exited $rc"
build_tree
run_gen --dry-run --fasta "$DATA/reference/grch37/genome.fa" >/dev/null 2>"$TMP/fa3.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--fasta needs a single --assembly' "$TMP/fa3.err" && ok "--fasta with --assembly both is refused" || bad "--fasta with both assemblies exited $rc"
run_gen --dry-run --assembly GRCh37 --fasta "$DATA/reference/grch37/genome.fa" >"$TMP/fa4.out" 2>/dev/null
rc=$?
n=$(count_docker_cmds "$TMP/fa4.out")
[[ "$rc" -eq 0 && "$n" -eq 5 ]] && ok "--fasta with a single assembly is accepted (5 commands)" || bad "--fasta with --assembly GRCh37: rc=$rc commands=$n"

echo
echo "COLLISION: an existing reference output is never overwritten without --force"
build_tree
mkdir -p "$SNP_GT/s03"
echo "pinned" >"$SNP_GT/s03/output.txt"
run_gen --dry-run >/dev/null 2>"$TMP/col.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "reference output already exists: $SNP_GT/s03/output.txt" "$TMP/col.err" &&
    ok "an existing output fails preflight and is named" || bad "an existing output exited $rc"
run_gen --dry-run --force >"$TMP/col2.out" 2>/dev/null
rc=$?
n=$(count_docker_cmds "$TMP/col2.out")
[[ "$rc" -eq 0 && "$n" -eq 10 ]] && ok "--force plans the collision as an overwrite" || bad "--force: rc=$rc commands=$n"
[[ "$(cat "$SNP_GT/s03/output.txt")" == "pinned" ]] && ok "a dry run with --force still writes nothing" || bad "dry run with --force modified the existing output"
run_gen --dry-run --format tab >"$TMP/col3.out" 2>/dev/null
rc=$?
[[ "$rc" -eq 0 ]] && ok "an existing default output does not collide with a tab generation (output.tsv is another file)" || bad "--format tab collided with output.txt (rc=$rc)"

echo
echo "ARGUMENTS: bad arguments are refused with a usage or named error"
out=$(bash "$GEN" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--data-dir or --corpus required' <<<"$out" && ok "no arguments: --data-dir or --corpus required, non-zero" || bad "no arguments exited $rc"
out=$(bash "$GEN" --data-dir "$TMP/nope" --dry-run 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'data-dir not a directory' <<<"$out" && ok "a non-directory --data-dir is a named error" || bad "non-directory data-dir exited $rc"
out=$(run_gen --dry-run --assembly hg19 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'invalid --assembly' <<<"$out" && ok "an unknown assembly is refused" || bad "--assembly hg19 exited $rc"
out=$(run_gen --dry-run --fork 0 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'positive integer' <<<"$out" && ok "--fork 0 is refused" || bad "--fork 0 exited $rc"
out=$(run_gen --dry-run --bogus 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'unknown argument: --bogus' <<<"$out" && ok "an unknown flag is refused" || bad "--bogus exited $rc"
out=$(run_gen --dry-run --volume-floor half 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'invalid --volume-floor' <<<"$out" && ok "an unknown --volume-floor is refused" || bad "--volume-floor half exited $rc"
out=$(run_gen --dry-run --perl-image-digest sha256:abc 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--perl-image-digest must be sha256:<64 hex>' <<<"$out" && ok "a malformed --perl-image-digest is refused" || bad "--perl-image-digest sha256:abc exited $rc"
out=$(run_gen --dry-run --suites s01 --snp-indel-gt-dir "$TMP/a|b" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "holds '|' or a newline" <<<"$out" && ! grep -q '^docker run ' <<<"$out" &&
    ok "a path holding '|' (the plan's record separator) is refused before any command is printed" || bad "a '|' in an output path exited $rc"

echo
echo "OUTPUT-VOLUME: on the real path, a VEP that annotates nothing leaves no output behind"
# A stub `docker` on PATH. `docker image inspect` prints the pinned digest (STUB_DIGEST
# overrides it), `docker pull` succeeds, and `docker run` resolves the /work/perl_out
# mount back to the host path and writes what VEP would at the -o target for the format
# the flags select: STUB_ROWS consequence rows for the default and tab formats (with a
# header block unless --no_headers); one line per input record for the VCF format, the
# last STUB_VCF_BARE of them without the CSQ field, as VEP writes a buffered record it
# could not annotate (OutputFactory/VCF.pm:299-361); one object per record for the JSON
# format, the last STUB_JSON_BARE without a consequences array (OutputFactory/JSON.pm:
# 213-260). Every run writes a warnings companion, and a _summary.txt under --stats_text
# unless STUB_NO_SUMMARY is set. STUB_ROWS=0 is the --full_cache_dir failure shape: exit
# 0, a warnings file, and no annotation rows. Each run's arguments are appended to
# STUB_LOG when it is set.
STUBDIR="$TMP/stubbin"
mkdir -p "$STUBDIR"
cat >"$STUBDIR/docker" <<'EOF'
#!/usr/bin/env bash
case "$1" in
image) echo "${STUB_DIGEST:-ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da}"; exit 0 ;;
pull) exit 0 ;;
run) ;;
*) echo "stub docker: unexpected subcommand $1" >&2; exit 99 ;;
esac
[[ -n "${STUB_LOG:-}" ]] && printf '%s\n' "$*" >>"$STUB_LOG"
host_out=""; host_in=""; out=""; fmt="default"; headers=1; stats=0; field="CSQ"
while [[ $# -gt 0 ]]; do
    case "$1" in
    -v)
        case "$2" in
        *:/work/perl_out) host_out="${2%%:*}" ;;
        *:/work/input.vcf:ro | *:/work/input.vcf.gz:ro) host_in="${2%%:*}" ;;
        esac
        shift 2 ;;
    -o) out="$2"; shift 2 ;;
    --tab) fmt="tab"; shift ;;
    --vcf) fmt="vcf"; shift ;;
    --json) fmt="json"; shift ;;
    --no_headers) headers=0; shift ;;
    --stats_text) stats=1; shift ;;
    --vcf_info_field | --vcf_info) field="$2"; shift 2 ;;
    --vcf_info_field=*) field="${1#*=}"; shift ;;
    *) shift ;;
    esac
done
[[ -n "$host_out" && -n "$out" ]] || { echo "stub docker: no output mapping" >&2; exit 98; }
case "$host_in" in
*.gz) n=$(gzip -dc "$host_in" | grep -vc '^#') || n=0 ;;
*) n=$(grep -vc '^#' "$host_in") || n=0 ;;
esac
target="$host_out/${out#/work/perl_out/}"
: >"$target"
case "$fmt" in
default | tab)
    [[ "$headers" -eq 1 ]] && printf '## ENSEMBL VARIANT EFFECT PREDICTOR v116.2\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n' >>"$target"
    for ((i = 1; i <= ${STUB_ROWS:-0}; i++)); do
        printf 'var%d\t21:%d\tG\tENSG\tENST\tTranscript\tmissense_variant\n' "$i" "$i" >>"$target"
    done ;;
vcf)
    [[ "$headers" -eq 1 ]] && printf '##fileformat=VCFv4.2\n##INFO=<ID=%s,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' "$field" >>"$target"
    for ((i = 1; i <= n; i++)); do
        if [[ "$i" -le $((n - ${STUB_VCF_BARE:-0})) ]]; then info="$field=G|missense_variant"; else info="."; fi
        printf '21\t%d\tvar%d\tA\tG\t.\tPASS\t%s\n' "$((10000000 + i))" "$i" "$info" >>"$target"
    done ;;
json)
    for ((i = 1; i <= n; i++)); do
        if [[ "$i" -le $((n - ${STUB_JSON_BARE:-0})) ]]; then
            printf '{"id":"var%d","transcript_consequences":[{"consequence_terms":["missense_variant"]}]}\n' "$i" >>"$target"
        else
            printf '{"id":"var%d"}\n' "$i" >>"$target"
        fi
    done ;;
esac
echo "WARNING: stub" >"${target}_warnings.txt"
if [[ "$stats" -eq 1 && -z "${STUB_NO_SUMMARY:-}" ]]; then printf '[VEP run statistics]\nVEP version (API)\t116.2 (116)\n' >"${target}_summary.txt"; fi
exit "${STUB_RC:-0}"
EOF
chmod +x "$STUBDIR/docker"

build_tree
STUB_ROWS=0 PATH="$STUBDIR:$PATH" run_gen --suites s01 >"$TMP/v0.out" 2>"$TMP/v0.err"
rc=$?
[[ "$rc" -ne 0 ]] && ok "zero annotation rows exits non-zero (rc=$rc)" || bad "zero annotation rows exited 0"
grep -q '0 annotation rows for 3 input records' "$TMP/v0.err" && ok "the row/record deficit is stated" || bad "stderr does not state the row/record deficit"
[[ ! -e "$SNP_GT/s01/output.txt" ]] && ok "no output.txt was left at the reference path" || bad "a non-annotating run left output.txt in place"
[[ -z "$(find "$SNP_GT" -name '.generate.*' 2>/dev/null)" ]] && ok "the scratch directory was removed" || bad "a scratch directory survived the failure"

build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 >"$TMP/v3.out" 2>"$TMP/v3.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "rows == records passes the gate and exits 0" || {
    bad "rows == records exited $rc"
    head -5 "$TMP/v3.err" | sed 's/^/         /'
}
[[ -f "$SNP_GT/s01/output.txt" ]] && ok "output.txt is in place" || bad "output.txt is missing after a passing run"
[[ -f "$SNP_GT/s01/output.warnings.log" && ! -e "$SNP_GT/s01/output.txt_warnings.txt" ]] &&
    ok "VEP's warnings companion is renamed to output.warnings.log" || bad "the warnings companion was not renamed"
[[ -f "$SNP_GT/s01/output.provenance.json" ]] && ok "output.provenance.json is written beside the output" || bad "no provenance sidecar"
grep -q '"input_records": 3' "$SNP_GT/s01/output.provenance.json" && grep -q '"output_data_lines": 3' "$SNP_GT/s01/output.provenance.json" &&
    grep -q '"fork": 16' "$SNP_GT/s01/output.provenance.json" && grep -q "ensemblorg/ensembl-vep@$PINNED_DIGEST" "$SNP_GT/s01/output.provenance.json" &&
    ok "the sidecar records records, rows, fork and the image digest" || bad "the sidecar is missing a field"
grep -q '"output_format": "default"' "$SNP_GT/s01/output.provenance.json" && grep -q '"vep_flags": "--offline .* --fasta <fa> --no_headers --no_stats --quiet --force_overwrite"' "$SNP_GT/s01/output.provenance.json" &&
    ok "the sidecar records the output format and the flags as run" || bad "the sidecar's output_format or vep_flags is wrong"
[[ -z "$(find "$SNP_GT" -name '.generate.*' 2>/dev/null)" ]] && ok "no scratch directory survives a passing run" || bad "a scratch directory survived a passing run"
n=$(find "$SNP_GT/s01" -name '*.txt' | wc -l | tr -d ' ')
[[ "$n" -eq 1 ]] && ok "exactly one *.txt in the suite directory, so a *.txt glob pairs only the annotation" || bad "$n *.txt files in the suite directory"
[[ ! -e "$SNP_GT/s01/output.summary.log" ]] && ok "no statistics report is kept without --stats" || bad "a summary was kept without --stats"

# The SV set is gated in AGGREGATE: two files, 4 + 2 records; a stub writing 3 rows per
# file gives 6 >= 6 and passes, even though the 4-record file alone would fail a per-file
# floor. Some synthetic SV files carry allele classes VEP declines to annotate, so a
# per-file floor would fail a legitimately sparse file.
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s07 >"$TMP/sv3.out" 2>"$TMP/sv3.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SV_GT/grch37/01_snv.txt" && -f "$SV_GT/grch37/09_breakends.txt" ]] &&
    ok "the SV set passes on aggregate rows >= records (6 >= 6) and both files land" || bad "SV aggregate gate: rc=$rc"
[[ -f "$SV_GT/grch37/09_breakends.provenance.json" && -f "$SV_GT/grch37/09_breakends.warnings.log" ]] &&
    ok "SV sidecars use the VCF basename (09_breakends.provenance.json / .warnings.log)" || bad "SV sidecar naming is wrong"
build_tree
STUB_ROWS=2 PATH="$STUBDIR:$PATH" run_gen --suites s07 >/dev/null 2>"$TMP/sv2.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$SV_GT/grch37/01_snv.txt" && ! -e "$SV_GT/grch37/09_breakends.txt" ]] &&
    ok "an SV set below the aggregate floor (4 < 6) leaves neither file in place" || bad "SV set below the floor: rc=$rc, files left behind"

build_tree
STUB_ROWS=3 STUB_RC=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 >/dev/null 2>"$TMP/rc.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'vep exited 3' "$TMP/rc.err" && [[ ! -e "$SNP_GT/s01/output.txt" ]] &&
    ok "a non-zero vep exit is fatal and leaves no output" || bad "vep exit 3: rc=$rc"

build_tree
mkdir -p "$SNP_GT/s01"
echo "pinned" >"$SNP_GT/s01/output.txt"
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 >/dev/null 2>"$TMP/keep.err"
rc=$?
[[ "$rc" -ne 0 && "$(cat "$SNP_GT/s01/output.txt")" == "pinned" ]] && ok "the real path also refuses to overwrite without --force" || bad "real path overwrote a pinned output (rc=$rc)"
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --force >/dev/null 2>&1
rc=$?
[[ "$rc" -eq 0 && "$(cat "$SNP_GT/s01/output.txt")" != "pinned" ]] && ok "--force replaces the pinned output" || bad "--force did not replace the output (rc=$rc)"

# A backslash sequence in an output path is a literal to the plan's group matching, so the
# group still finds its records and the output lands at the path as written.
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --snp-indel-gt-dir "$TMP/a\\nb" >/dev/null 2>"$TMP/bs.err"
rc=$?
[[ "$rc" -eq 0 && -f "$TMP/a\\nb/s01/output.txt" ]] && ok "an output path holding a backslash sequence is matched literally and written" || {
    bad "backslash in an output path: rc=$rc"
    head -3 "$TMP/bs.err" | sed 's/^/         /'
}

echo
echo "VOLUME-FLOOR: --volume-floor any keeps the gate at one annotated unit for flags that drop rows"
build_tree
STUB_ROWS=1 PATH="$STUBDIR:$PATH" run_gen --suites s01 --flags "--coding_only" >/dev/null 2>"$TMP/vf1.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$SNP_GT/s01/output.txt" ]] && grep -q -- '--volume-floor any' "$TMP/vf1.err" &&
    ok "1 row for 3 records fails the records floor and the error names --volume-floor any" || bad "records floor with 1 < 3: rc=$rc"
build_tree
STUB_ROWS=1 PATH="$STUBDIR:$PATH" run_gen --suites s01 --flags "--coding_only" --volume-floor any >/dev/null 2>"$TMP/vf2.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.txt" ]] && ok "--volume-floor any passes 1 row for 3 records" || bad "--volume-floor any with 1 row: rc=$rc"
build_tree
STUB_ROWS=0 PATH="$STUBDIR:$PATH" run_gen --suites s01 --volume-floor any >/dev/null 2>"$TMP/vf3.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$SNP_GT/s01/output.txt" ]] && grep -q 'no annotated unit' "$TMP/vf3.err" && ok "--volume-floor any still fails zero rows" || bad "--volume-floor any passed zero rows (rc=$rc)"

echo
echo "FORMAT-GATES: the VCF and JSON formats are gated on every written record being annotated"
# VCF.pm:299-361 writes every buffered record and appends the CSQ field only to one with
# consequences; JSON.pm:213-260 and :384 write one object per record and the consequence
# arrays only on an annotated one. A record VEP skipped before its input buffer is absent
# from both, so the records floor for these formats is "no bare record", not a row count.
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format vcf >"$TMP/fv.out" 2>"$TMP/fv.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.vcf" ]] && ok "--format vcf lands output.vcf (3 of 3 records carry CSQ)" || {
    bad "--format vcf: rc=$rc"
    head -5 "$TMP/fv.err" | sed 's/^/         /'
}
[[ -f "$SNP_GT/s01/output.vcf.provenance.json" && -f "$SNP_GT/s01/output.vcf.warnings.log" && ! -e "$SNP_GT/s01/output.txt" ]] &&
    ok "the VCF sidecars are output.vcf.provenance.json and output.vcf.warnings.log, and no output.txt is written" || bad "VCF sidecar naming is wrong"
grep -q '"output_format": "vcf"' "$SNP_GT/s01/output.vcf.provenance.json" && grep -q -- '--force_overwrite --vcf"' "$SNP_GT/s01/output.vcf.provenance.json" &&
    grep -q '"output_annotated_units": 3' "$SNP_GT/s01/output.vcf.provenance.json" &&
    ok "the VCF sidecar records the format, the --vcf flag and the annotated-record count" || bad "the VCF sidecar's fields are wrong"
grep -q '^#CHROM' "$SNP_GT/s01/output.vcf" && ok "the VCF reference keeps its header lines" || bad "the VCF reference lost its header"
build_tree
STUB_ROWS=3 STUB_VCF_BARE=1 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format vcf >/dev/null 2>"$TMP/fvb.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$SNP_GT/s01/output.vcf" ]] && grep -q '1 of 3 written records carry no CSQ field' "$TMP/fvb.err" &&
    ok "a written VCF record without CSQ fails the records floor and is counted" || bad "bare VCF record: rc=$rc"
build_tree
STUB_ROWS=3 STUB_VCF_BARE=1 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format vcf --volume-floor any >/dev/null 2>"$TMP/fvb2.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.vcf" ]] && ok "--volume-floor any accepts a bare VCF record (a kept record with no row, as under --dont_skip)" || bad "bare VCF record under --volume-floor any: rc=$rc"
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format vcf --flags "--vcf_info_field ANN" >/dev/null 2>"$TMP/fva.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.vcf" ]] && grep -q 'ANN=' "$SNP_GT/s01/output.vcf" && ok "--vcf_info_field in --flags moves the gate to that INFO key" || bad "--vcf_info_field ANN: rc=$rc"
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format vcf --flags "--vcf_info_field=ANN" >/dev/null 2>"$TMP/fva2.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.vcf" ]] && grep -q 'ANN=' "$SNP_GT/s01/output.vcf" && ok "the --vcf_info_field=NAME spelling moves the gate's key too" || bad "--vcf_info_field=ANN: rc=$rc"
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format vcf --flags "--vcf_info ANN" >/dev/null 2>"$TMP/fva3.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.vcf" ]] && grep -q 'ANN=' "$SNP_GT/s01/output.vcf" && ok "an abbreviation vep accepts for --vcf_info_field (--vcf_info) moves the gate's key" || bad "--vcf_info ANN: rc=$rc"
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s07 --format json >/dev/null 2>"$TMP/fj.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SV_GT/grch37/01_snv.jsonl" && -f "$SV_GT/grch37/09_breakends.jsonl" && -f "$SV_GT/grch37/01_snv.jsonl.provenance.json" ]] &&
    ok "--format json lands <basename>.jsonl per SV file with its sidecar" || bad "--format json on the SV set: rc=$rc"
n=$(grep -c '_consequences' "$SV_GT/grch37/01_snv.jsonl") || n=0
[[ "$n" -eq 4 ]] && ok "the JSON reference carries one annotated object per record (4)" || bad "$n annotated JSON objects for 4 records"
build_tree
STUB_ROWS=3 STUB_JSON_BARE=2 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format json >/dev/null 2>"$TMP/fjb.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$SNP_GT/s01/output.jsonl" ]] && grep -q '2 of 3 written objects carry no consequences array' "$TMP/fjb.err" &&
    ok "a JSON object without consequences fails the records floor and is counted" || bad "bare JSON object: rc=$rc"
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --format all >/dev/null 2>"$TMP/fall.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.txt" && -f "$SNP_GT/s01/output.tsv" && -f "$SNP_GT/s01/output.vcf" && -f "$SNP_GT/s01/output.jsonl" ]] &&
    ok "--format all lands the four outputs of one suite side by side" || bad "--format all on s01: rc=$rc"
n=$(find "$SNP_GT/s01" -name '*.txt' | wc -l | tr -d ' ')
[[ "$n" -eq 1 ]] && ok "with four formats present, still exactly one *.txt (the default output)" || bad "$n *.txt files beside the four formats"
grep -q '^#Uploaded_variation' "$SNP_GT/s01/output.tsv" && ! grep -q '^#' "$SNP_GT/s01/output.txt" &&
    ok "the tab reference keeps its column line and the default reference has none" || bad "header handling differs from the contract"
n=$(find "$SNP_GT/s01" -name '*.provenance.json' | wc -l | tr -d ' ')
[[ "$n" -eq 4 ]] && ok "each format has its own provenance sidecar (4)" || bad "$n provenance sidecars for 4 formats"

echo
echo "STATS-RUN: --stats keeps the report as <stem>.summary.log and fails closed without one"
build_tree
STUB_ROWS=3 PATH="$STUBDIR:$PATH" run_gen --suites s01 --stats >/dev/null 2>"$TMP/st.err"
rc=$?
[[ "$rc" -eq 0 && -f "$SNP_GT/s01/output.summary.log" && ! -e "$SNP_GT/s01/output.txt_summary.txt" ]] &&
    ok "the statistics report lands as output.summary.log" || bad "--stats run: rc=$rc"
n=$(find "$SNP_GT/s01" -name '*.txt' | wc -l | tr -d ' ')
[[ "$n" -eq 1 ]] && ok "the kept report leaves exactly one *.txt in the suite directory" || bad "$n *.txt files with --stats"
grep -q '"vep_flags": ".*--stats_text --quiet --force_overwrite"' "$SNP_GT/s01/output.provenance.json" && ok "the sidecar records --stats_text" || bad "the sidecar does not record --stats_text"
build_tree
STUB_ROWS=3 STUB_NO_SUMMARY=1 PATH="$STUBDIR:$PATH" run_gen --suites s01 --stats >/dev/null 2>"$TMP/st2.err"
rc=$?
[[ "$rc" -ne 0 && ! -e "$SNP_GT/s01/output.txt" ]] && grep -q 'wrote no output.txt_summary.txt' "$TMP/st2.err" &&
    ok "--stats with no report written is fatal and leaves no output" || bad "--stats without a report: rc=$rc"

echo
echo "DIGEST: the pulled image must carry the pinned digest before anything runs"
build_tree
other="sha256:0000000000000000000000000000000000000000000000000000000000000000"
: >"$TMP/stub.log"
STUB_ROWS=3 STUB_LOG="$TMP/stub.log" PATH="$STUBDIR:$PATH" run_gen --suites s01 --perl-image-digest "$other" >/dev/null 2>"$TMP/dg.err"
rc=$?
[[ "$rc" -eq 1 ]] && grep -q 'image digest mismatch' "$TMP/dg.err" && grep -q "$other" "$TMP/dg.err" &&
    ok "a digest mismatch exits 1 and names the expected digest" || bad "digest mismatch: rc=$rc"
[[ ! -s "$TMP/stub.log" && ! -e "$SNP_GT/s01" ]] && ok "no container ran and nothing was written before the digest check" || bad "a container ran despite the digest mismatch"
STUB_ROWS=3 STUB_DIGEST="ensemblorg/ensembl-vep@$other" PATH="$STUBDIR:$PATH" run_gen --suites s01 >/dev/null 2>"$TMP/dg2.err"
rc=$?
[[ "$rc" -eq 1 ]] && grep -q 'image digest mismatch' "$TMP/dg2.err" && ok "an image whose RepoDigests lack the pinned digest is refused" || bad "moved-tag digest: rc=$rc"
STUB_ROWS=3 STUB_DIGEST="ensemblorg/ensembl-vep@$other" PATH="$STUBDIR:$PATH" run_gen --suites s01 --perl-image-digest "$other" >/dev/null 2>"$TMP/dg3.err"
rc=$?
[[ "$rc" -eq 0 ]] && grep -q "ensemblorg/ensembl-vep@$other" "$SNP_GT/s01/output.provenance.json" && ok "--perl-image-digest names another image's digest and the sidecar records it" || bad "--perl-image-digest override: rc=$rc"
build_tree
STUB_ROWS=3 VEP_PERL_IMAGE_DIGEST="$other" PATH="$STUBDIR:$PATH" run_gen --suites s01 >/dev/null 2>"$TMP/dg4.err"
rc=$?
[[ "$rc" -eq 1 ]] && grep -q 'image digest mismatch' "$TMP/dg4.err" && ok "VEP_PERL_IMAGE_DIGEST sets the expected digest" || bad "VEP_PERL_IMAGE_DIGEST: rc=$rc"

echo
echo "CORPUS: --corpus writes a golden corpus's four expected outputs and its provenance"
CORPUS="$TMP/corpus"
build_corpus() { # build_corpus [with-fasta]
    rm -rf "$CORPUS"
    mkdir -p "$CORPUS"
    write_vcf "$CORPUS/variants.vcf" 3
    {
        echo '{'
        echo ' "release": "116",'
        echo ' "assembly": "GRCh38",'
        echo ' "flags": ['
        echo '  "--protein",'
        echo '  "--max_sv_size",'
        echo '  "-1"'
        if [[ "${1:-}" == "with-fasta" ]]; then
            echo ' ],'
            echo ' "fasta": "reference.fa.gz"'
        else
            echo ' ]'
        fi
        echo '}'
    } >"$CORPUS/manifest.json"
    if [[ "${1:-}" == "with-fasta" ]]; then
        printf '>21\nACGT\n' | gzip -c >"$CORPUS/reference.fa.gz"
        printf '21\t4\t4\t4\t5\n' >"$CORPUS/reference.fa.fai"
    fi
}
run_corpus() { bash "$GEN" --corpus "$CORPUS" --perl-cache-dir "$CACHE" "$@"; }

build_tree
build_corpus
run_corpus --dry-run >"$TMP/cd.out" 2>"$TMP/cd.err"
rc=$?
n=$(count_docker_cmds "$TMP/cd.out")
[[ "$rc" -eq 0 && "$n" -eq 4 ]] && ok "a corpus dry run plans four commands, one per format" || {
    bad "corpus dry run: rc=$rc commands=$n"
    head -5 "$TMP/cd.err" | sed 's/^/         /'
}
n=$(grep -c -- '--fork 1 ' "$TMP/cd.out") || n=0
[[ "$n" -eq 4 ]] && ok "a corpus runs --fork 1" || bad "$n of 4 corpus commands pass --fork 1"
grep -q -- '--no_headers' "$TMP/cd.out" && bad "a corpus command drops the headers the golden tests compare" || ok "a corpus keeps every format's headers"
n=$(grep -cE -- '--force_overwrite --protein --max_sv_size -1( --(tab|vcf|json))? ?$' "$TMP/cd.out") || n=0
[[ "$n" -eq 4 ]] && ok "the manifest's flags are appended to every command" || bad "$n of 4 corpus commands carry the manifest flags"
grep -q -- '--fasta' "$TMP/cd.out" && bad "a corpus without a manifest fasta passes --fasta" || ok "no --fasta without a manifest fasta"
d=$(grep -c -- '-o /work/perl_out/default.txt ' "$TMP/cd.out") || d=0
t=$(grep -c -- '-o /work/perl_out/tab.txt ' "$TMP/cd.out") || t=0
v=$(grep -c -- '-o /work/perl_out/vcf.vcf ' "$TMP/cd.out") || v=0
j=$(grep -c -- '-o /work/perl_out/json.jsonl ' "$TMP/cd.out") || j=0
[[ "$d" -eq 1 && "$t" -eq 1 && "$v" -eq 1 && "$j" -eq 1 ]] && ok "the outputs are default.txt, tab.txt, vcf.vcf and json.jsonl" || bad "corpus output names: default=$d tab=$t vcf=$v json=$j"
grep -q "\[corpus GRCh38 json\] -> $CORPUS/expected/json.jsonl.gz" "$TMP/cd.out" && ok "the dry run names the gzipped destination under expected/" || bad "the corpus destination is not expected/<name>.gz"
[[ ! -e "$CORPUS/expected" && -z "$(find "$CORPUS" -name '.generate.*')" ]] && ok "a corpus dry run writes nothing" || bad "a corpus dry run wrote into the corpus"

build_corpus with-fasta
run_corpus --dry-run >"$TMP/cdf.out" 2>/dev/null
n=$(grep -cE -- "-v $CORPUS/\.generate\.[^/]*/reference\.fa:/work/ref\.fa:ro .* --fasta /work/ref\.fa" "$TMP/cdf.out") || n=0
[[ "$n" -eq 4 ]] && ok "a manifest fasta is staged decompressed and passed as --fasta on every command" || bad "$n of 4 corpus commands mount the staged FASTA"
rm "$CORPUS/reference.fa.fai"
run_corpus --dry-run >/dev/null 2>"$TMP/cdf2.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'corpus FASTA or its index not found' "$TMP/cdf2.err" && ok "a manifest fasta without its .fai fails preflight" || bad "corpus fasta without .fai: rc=$rc"

build_corpus with-fasta
STUB_ROWS=5 PATH="$STUBDIR:$PATH" run_corpus --stats >"$TMP/cr.out" 2>"$TMP/cr.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "a corpus run exits 0" || {
    bad "corpus run exited $rc"
    head -5 "$TMP/cr.err" | sed 's/^/         /'
}
[[ -f "$CORPUS/expected/default.txt.gz" && -f "$CORPUS/expected/tab.txt.gz" && -f "$CORPUS/expected/vcf.vcf.gz" && -f "$CORPUS/expected/json.jsonl.gz" ]] &&
    ok "expected/{default.txt,tab.txt,vcf.vcf,json.jsonl}.gz are written" || bad "the four expected files are not all present"
n=$(find "$CORPUS/expected" -name '*.txt' -o -name '*.vcf' -o -name '*.jsonl' | wc -l | tr -d ' ')
[[ "$n" -eq 0 ]] && ok "no uncompressed expected output is left beside the gzipped ones" || bad "$n plain expected file(s) beside the .gz"
[[ -f "$CORPUS/provenance.json" ]] && ok "provenance.json is written" || bad "no provenance.json"
python3 - "$CORPUS" <<'EOF' >"$TMP/cr.check" 2>&1
import gzip, hashlib, json, sys
from pathlib import Path
c = Path(sys.argv[1])
p = json.loads((c / "provenance.json").read_text())
names = ["default.txt", "tab.txt", "vcf.vcf", "json.jsonl"]
assert list(p["sha256"]) == names, p["sha256"]
for n in names:
    assert hashlib.sha256(gzip.open(c / "expected" / (n + ".gz"), "rb").read()).hexdigest() == p["sha256"][n], n
assert p["variants_vcf_sha256"] == hashlib.sha256((c / "variants.vcf").read_bytes()).hexdigest()
assert p["vep_image"] == "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da", p["vep_image"]
assert p["vep_image_tag"] == "ensemblorg/ensembl-vep:release_116.2"
assert p["vep_version_reported"] == "ensembl-vep          : 116.2", p["vep_version_reported"]
assert p["cache_version"] == 116
assert p["perl_cache"] == "homo_sapiens/116_GRCh38"
assert p["fasta"] == "reference.fa.gz"
assert p["input"] == "variants.vcf"
cmd = p["vep_command"]
assert cmd.startswith("vep -i variants.vcf -o <out> --offline --cache --dir_cache <cache> --cache_version 116 --species homo_sapiens --assembly GRCh38 --format vcf --buffer_size 5000 --fork 1 --fasta reference.fa --stats_text --quiet --force_overwrite --protein --max_sv_size -1 [--tab|--vcf|--json]"), cmd
print("OK")
EOF
grep -q '^OK$' "$TMP/cr.check" && ok "provenance.json carries the image digest, the vep command and a matching sha256 per expected file" || {
    bad "provenance.json content is wrong"
    tail -3 "$TMP/cr.check" | sed 's/^/         /'
}
[[ -f "$CORPUS/expected/default.summary.log" && -f "$CORPUS/expected/json.summary.log" && -f "$CORPUS/expected/vcf.warnings.log" ]] &&
    ok "the statistics reports and warnings land beside the expected files as <format>.summary.log and <format>.warnings.log" || bad "corpus sidecars are missing"
[[ -z "$(find "$CORPUS" -name '.generate.*')" ]] && ok "the corpus scratch directory is removed" || bad "a scratch directory survived the corpus run"
n=$(gzip -dc "$CORPUS/expected/default.txt.gz" | grep -c '^##') || n=0
[[ "$n" -eq 1 ]] && ok "the default expected output keeps its header lines" || bad "the default expected output has $n header lines"

# Regeneration: refused without --force; with --force the keys the script does not own survive.
STUB_ROWS=5 PATH="$STUBDIR:$PATH" run_corpus >/dev/null 2>"$TMP/cr2.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q "reference output already exists: $CORPUS/expected/default.txt.gz" "$TMP/cr2.err" &&
    ok "an existing expected file refuses a regeneration without --force" || bad "corpus regeneration without --force: rc=$rc"
python3 - "$CORPUS" <<'EOF'
import json, sys
from pathlib import Path
p = Path(sys.argv[1]) / "provenance.json"
d = json.loads(p.read_text())
d["divergences_classified_against"] = {"engine": "vep-rs", "binary_md5": "abc"}
d["input"] = "hand-built records"
d["fasta"] = "chromosome 21 of the primary assembly"
d["sha256"]["default.txt"] = "stale"
p.write_text(json.dumps(d, indent=1) + "\n")
EOF
echo "stale plain copy" >"$CORPUS/expected/default.txt"
STUB_ROWS=6 PATH="$STUBDIR:$PATH" run_corpus --force >/dev/null 2>"$TMP/cr3.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "--force regenerates an existing corpus" || bad "corpus --force: rc=$rc"
python3 - "$CORPUS" <<'EOF' >"$TMP/cr3.check" 2>&1
import gzip, hashlib, json, sys
from pathlib import Path
c = Path(sys.argv[1])
p = json.loads((c / "provenance.json").read_text())
assert p["divergences_classified_against"] == {"engine": "vep-rs", "binary_md5": "abc"}, p.keys()
assert p["input"] == "hand-built records" and p["fasta"] == "chromosome 21 of the primary assembly", (p["input"], p["fasta"])
assert p["sha256"]["default.txt"] == hashlib.sha256(gzip.open(c / "expected" / "default.txt.gz", "rb").read()).hexdigest()
assert list(p)[:3] == ["vep_image", "vep_image_tag", "vep_version_reported"], list(p)
print("OK")
EOF
grep -q '^OK$' "$TMP/cr3.check" && ok "a regeneration keeps a key it does not own and the input and fasta prose, and rewrites the sha256 map" || bad "regeneration lost or kept the wrong keys"
[[ ! -e "$CORPUS/expected/default.txt" ]] && ok "--force removes a stale plain expected file, which the golden tests would read before the .gz" || bad "a stale plain expected file survived --force"

# A failing format leaves the committed expected files untouched.
before=$(sha256_file "$CORPUS/expected/json.jsonl.gz")
STUB_ROWS=6 STUB_JSON_BARE=1 PATH="$STUBDIR:$PATH" run_corpus --force >/dev/null 2>"$TMP/cr4.err"
rc=$?
[[ "$rc" -ne 0 && "$(sha256_file "$CORPUS/expected/json.jsonl.gz")" == "$before" ]] && grep -q 'carry no consequences array' "$TMP/cr4.err" &&
    ok "a bare JSON object fails the corpus and the committed expected files are untouched" || bad "corpus with a bare JSON object: rc=$rc"
# A failure between the container runs and the record write leaves expected/ untouched:
# a `python3` on PATH that fails the provenance step alone (the one script importing
# hashlib) and hands every other call to the real interpreter.
PYSTUB="$TMP/pystub"
mkdir -p "$PYSTUB"
cat >"$PYSTUB/python3" <<EOF
#!/usr/bin/env bash
script=\$(cat)
if [[ -n "\${STUB_PY_FAIL_PROVENANCE:-}" && "\$script" == *hashlib* ]]; then echo "stub python3: provenance step failed" >&2; exit 1; fi
exec $(command -v python3) "\$@" <<<"\$script"
EOF
chmod +x "$PYSTUB/python3"
before_default=$(sha256_file "$CORPUS/expected/default.txt.gz")
STUB_ROWS=7 STUB_PY_FAIL_PROVENANCE=1 PATH="$PYSTUB:$STUBDIR:$PATH" run_corpus --force >/dev/null 2>"$TMP/cr6.err"
rc=$?
[[ "$rc" -ne 0 && "$(sha256_file "$CORPUS/expected/default.txt.gz")" == "$before_default" && "$(sha256_file "$CORPUS/expected/json.jsonl.gz")" == "$before" ]] &&
    grep -q 'provenance step failed' "$TMP/cr6.err" && [[ -z "$(find "$CORPUS" -name '.generate.*')" ]] &&
    ok "a provenance step that fails after the runs leaves expected/ untouched and no scratch directory" || bad "failing provenance step: rc=$rc"
STUB_ROWS=6 STUB_JSON_BARE=1 PATH="$STUBDIR:$PATH" run_corpus --force --volume-floor any >/dev/null 2>"$TMP/cr5.err"
rc=$?
[[ "$rc" -eq 0 ]] && ok "--volume-floor any applies to a corpus whose flags keep a record without rows" || bad "corpus --volume-floor any: rc=$rc"

build_corpus
out=$(run_corpus --dry-run --suites s01 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- 'drop --suites' <<<"$out" && ok "--corpus refuses a suite option (--suites)" || bad "--corpus --suites exited $rc"
out=$(run_corpus --dry-run --flags "--everything" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'takes its flags from manifest.json' <<<"$out" && ok "--corpus refuses --flags (the manifest's flags are the corpus's)" || bad "--corpus --flags exited $rc"
out=$(run_corpus --dry-run --format tab 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'writes all four formats' <<<"$out" && ok "--corpus refuses a single --format" || bad "--corpus --format tab exited $rc"
mkdir -p "$TMP/corpus_str"
cp "$CORPUS/variants.vcf" "$TMP/corpus_str/variants.vcf"
printf '{"assembly": "GRCh38", "flags": "--protein --max_sv_size -1"}\n' >"$TMP/corpus_str/manifest.json"
out=$(bash "$GEN" --corpus "$TMP/corpus_str" --perl-cache-dir "$CACHE" --dry-run 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'manifest flags must be a list of words' <<<"$out" && ! grep -q '^docker run ' <<<"$out" &&
    ok "a manifest whose flags is one string rather than a list is refused (a string would be split into characters)" || bad "string-valued manifest flags exited $rc"
out=$(bash "$GEN" --corpus "$CORPUS" --dry-run 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--corpus needs --perl-cache-dir' <<<"$out" && ok "--corpus without a cache location is refused" || bad "--corpus without --perl-cache-dir exited $rc"
out=$(bash "$GEN" --corpus "$TMP/nope" --perl-cache-dir "$CACHE" --dry-run 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'corpus not a directory' <<<"$out" && ok "a non-directory --corpus is a named error" || bad "non-directory corpus exited $rc"
build_corpus
printf '{"vep_image": "x",\n' >"$CORPUS/provenance.json"
: >"$TMP/stub.log"
STUB_ROWS=3 STUB_LOG="$TMP/stub.log" PATH="$STUBDIR:$PATH" run_corpus >/dev/null 2>"$TMP/cpv.err"
rc=$?
[[ "$rc" -ne 0 && ! -s "$TMP/stub.log" && ! -e "$CORPUS/expected" ]] && grep -q 'provenance.json' "$TMP/cpv.err" &&
    ok "a damaged existing provenance.json stops the corpus run before any container runs" || bad "damaged provenance.json: rc=$rc"
build_corpus
echo '[1, 2]' >"$CORPUS/provenance.json"
out=$(run_corpus --dry-run 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'not a JSON object' <<<"$out" && ok "a provenance.json that is not an object is refused by name" || bad "array provenance.json exited $rc"
build_corpus
rm -rf "$CACHE/homo_sapiens/116_GRCh38/MT"
out=$(run_corpus --dry-run 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'missing primary contigs: MT' <<<"$out" && ok "a corpus run checks the cache slice of the manifest's assembly" || bad "corpus against an incomplete cache exited $rc"

echo
echo "FLAG: the generator itself never passes --full_cache_dir"
_bad_live=$(grep -e '--full_cache_dir' "$GEN" | grep -vc '^[[:space:]]*#') || _bad_live=0
[[ "$_bad_live" -eq 0 ]] && ok "no live --full_cache_dir in generate_reference_output.sh" || bad "$_bad_live live line(s) pass --full_cache_dir"
_dir_cache=$(grep -e '--dir_cache' "$GEN" | grep -vc '^[[:space:]]*#') || _dir_cache=0
[[ "$_dir_cache" -ge 1 ]] && ok "the invocation passes --dir_cache" || bad "no live --dir_cache in the generator"

echo
echo "----------------------------------------"
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
