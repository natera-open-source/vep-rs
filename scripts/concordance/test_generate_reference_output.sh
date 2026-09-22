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
            mkdir -p "$CACHE/homo_sapiens/115_$asm/$c"
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
n=$(grep -c -- '--cache_version 115' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "every command passes --cache_version 115" || bad "$n of 10 commands pass --cache_version 115"
n=$(grep -c -- '--no_headers --no_stats --quiet --force_overwrite' "$TMP/dry.out") || n=0
[[ "$n" -eq 10 ]] && ok "every command carries the harness's output flags" || bad "$n of 10 commands carry --no_headers --no_stats --quiet --force_overwrite"

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
grep -q "DRY-RUN \[snp_indel GRCh37\] -> $SNP_GT/s01/output.txt" "$TMP/dry.out" && ok "s01 -> <snp-indel-gt>/s01/output.txt" || bad "s01 output path is not <snp-indel-gt>/s01/output.txt"
grep -q "DRY-RUN \[sv GRCh37\] -> $SV_GT/grch37/01_snv.txt" "$TMP/dry.out" && ok "01_snv.vcf -> <sv-gt>/grch37/01_snv.txt" || bad "SV output path is not <sv-gt>/<asm>/<basename>.txt"
grep -q "DRY-RUN \[sv GRCh38\] -> $SV_GT/grch38/09_breakends.txt" "$TMP/dry.out" && ok "09_breakends.vcf.gz -> <sv-gt>/grch38/09_breakends.txt (both suffixes stripped)" || bad ".vcf.gz basename not stripped to <basename>.txt"
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
rm -rf "$CACHE/homo_sapiens/115_GRCh37/5" "$CACHE/homo_sapiens/115_GRCh37/X"
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
rm -rf "$CACHE/homo_sapiens/115_GRCh38"
run_gen --dry-run >/dev/null 2>"$TMP/c3.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q 'Perl cache slice not found' "$TMP/c3.err" && ok "an absent cache slice exits non-zero and says so" || bad "an absent cache slice exited $rc"

build_tree
run_gen --dry-run --cache-version 114 >/dev/null 2>"$TMP/c4.err"
rc=$?
[[ "$rc" -ne 0 ]] && grep -q '114_GRCh37' "$TMP/c4.err" && ok "--cache-version selects the slice that is checked" || bad "--cache-version 114 against a 115 cache exited $rc"

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

echo
echo "ARGUMENTS: bad arguments are refused with a usage or named error"
out=$(bash "$GEN" 2>&1)
rc=$?
[[ "$rc" -ne 0 ]] && grep -q -- '--data-dir required' <<<"$out" && ok "no arguments: --data-dir required, non-zero" || bad "no arguments exited $rc"
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

echo
echo "OUTPUT-VOLUME: on the real path, a VEP that annotates nothing leaves no output behind"
# A stub `docker` on PATH. `docker image inspect` prints a digest, `docker pull` succeeds,
# and `docker run` resolves the /work/perl_out mount back to the host path and writes
# STUB_ROWS data lines (plus a warnings companion) to the -o target, exactly where VEP
# would. STUB_ROWS=0 is the --full_cache_dir failure shape: exit 0, a warnings file, and
# no annotation rows.
STUBDIR="$TMP/stubbin"
mkdir -p "$STUBDIR"
cat >"$STUBDIR/docker" <<'EOF'
#!/usr/bin/env bash
case "$1" in
image) echo "ensemblorg/ensembl-vep@sha256:stubdigest"; exit 0 ;;
pull) exit 0 ;;
run) ;;
*) echo "stub docker: unexpected subcommand $1" >&2; exit 99 ;;
esac
host_out=""
out=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    -v)
        case "$2" in
        *:/work/perl_out) host_out="${2%%:*}" ;;
        esac
        shift 2 ;;
    -o) out="$2"; shift 2 ;;
    *) shift ;;
    esac
done
[[ -n "$host_out" && -n "$out" ]] || { echo "stub docker: no output mapping" >&2; exit 98; }
target="$host_out/${out#/work/perl_out/}"
: >"$target"
for ((i = 1; i <= ${STUB_ROWS:-0}; i++)); do
    printf 'var%d\t21:%d\tG\tENSG\tENST\tTranscript\tmissense_variant\n' "$i" "$i" >>"$target"
done
echo "WARNING: stub" >"${target}_warnings.txt"
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
    grep -q '"fork": 16' "$SNP_GT/s01/output.provenance.json" && grep -q 'stubdigest' "$SNP_GT/s01/output.provenance.json" &&
    ok "the sidecar records records, rows, fork and the image digest" || bad "the sidecar is missing a field"
[[ -z "$(find "$SNP_GT" -name '.generate.*' 2>/dev/null)" ]] && ok "no scratch directory survives a passing run" || bad "a scratch directory survived a passing run"
n=$(find "$SNP_GT/s01" -name '*.txt' | wc -l | tr -d ' ')
[[ "$n" -eq 1 ]] && ok "exactly one *.txt in the suite directory, so a *.txt glob pairs only the annotation" || bad "$n *.txt files in the suite directory"

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
