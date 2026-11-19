#!/usr/bin/env bash
# Every data path run_concordance.sh passes to an engine must be a file
# setup_plugin_data.sh stages, on both assemblies, or the user's own file where a plugin
# reads one. Both scripts read the table in
# plugin_data_layout.sh; this test extracts the harness's build_plugin_args function and
# checks its output against plugin_staged_files, so a name edited in one place without the
# other fails here rather than at run time as an engine reading a missing file. It also
# evaluates the harness's PERL_VEP_SIGNATURE line, the Perl output cache's key, which must
# change with the plugin argument strings.
#
# Run: bash scripts/data/test_plugin_data_layout.sh
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
# shellcheck source=scripts/data/plugin_data_layout.sh
source "$HERE/plugin_data_layout.sh"

PASS=0; FAIL=0
ok() { PASS=$((PASS + 1)); echo "  ok   $1"; }
bad() { FAIL=$((FAIL + 1)); echo "  FAIL $1"; }

# The harness's function, extracted by name so the script's argument parsing never runs.
harness="$ROOT/scripts/concordance/run_concordance.sh"
fn_start=$(grep -n '^build_plugin_args()' "$harness" | cut -d: -f1)
fn_end=$(awk -v s="$fn_start" 'NR>s && /^}/{print NR; exit}' "$harness")
[[ -n "$fn_start" && -n "$fn_end" ]] || { echo "ERROR: [test_plugin_data_layout] build_plugin_args not found in $harness" >&2; exit 1; }
eval "$(sed -n "${fn_start},${fn_end}p" "$harness")"
runs="$ROOT/scripts/concordance/plugin_gate_runs.toml"

grep -q 'source "${SCRIPT_DIR}/../data/plugin_data_layout.sh"' "$harness" && ok "run_concordance.sh sources the shared layout" || bad "run_concordance.sh does not source plugin_data_layout.sh"
grep -q 'plugin_staged_files "$1" "${ASSEMBLY_LC}"' "$ROOT/scripts/data/setup_plugin_data.sh" && ok "setup_plugin_data.sh stages from the shared layout" || bad "setup_plugin_data.sh has its own file table"

PERL_VEP_MODE=native
for ASSEMBLY_LC in grch37 grch38; do
    ASSEMBLY=$(echo "$ASSEMBLY_LC" | sed 's/grch/GRCh/')
    export ASSEMBLY ASSEMBLY_LC
    for plugin in REVEL gnomADc AlphaMissense LoFtool pLI GWAS LoFTEE AncestralAllele OpenTargets Conservation LOEUF DosageSensitivity; do
        [[ ( "$plugin" == GWAS || "$plugin" == OpenTargets || "$plugin" == Conservation ) && "$ASSEMBLY_LC" == grch37 ]] && continue
        staged=$(plugin_staged_files "$plugin" "$ASSEMBLY_LC")
        for target in rust perl; do
            arg=$(LOFTEE_PERL_PATH=/lp build_plugin_args "$plugin" /pd "$target") || { bad "$plugin/$ASSEMBLY_LC/$target: build_plugin_args failed"; continue; }
            if [[ "$plugin" == LoFTEE ]]; then
                # Perl's master branch defaults run_splice_predictions to 1 (LoF.pm:98) and vep-rs has
                # no splice layer and refuses a GRCh37 run that leaves it on: the Perl string pins it
                # off once, the vep-rs string once on GRCh37 and never names it on GRCh38.
                n=$(grep -o 'run_splice_predictions' <<<"$arg" | grep -c .) || n=0
                if [[ "$target" == perl ]]; then
                    if [[ "$arg" == *,run_splice_predictions:0 && "$n" -eq 1 ]]; then ok "LoFTEE on $ASSEMBLY (perl): the splice-prediction layer is pinned off"; else bad "LoFTEE on $ASSEMBLY (perl): run_splice_predictions:0 is not the single last parameter: $arg"; fi
                elif [[ "$ASSEMBLY_LC" == grch37 ]]; then
                    if [[ "$arg" == *,run_splice_predictions=0 && "$n" -eq 1 ]]; then ok "LoFTEE on $ASSEMBLY (rust): the splice-prediction layer is pinned off"; else bad "LoFTEE on $ASSEMBLY (rust): run_splice_predictions=0 is not the single last parameter: $arg"; fi
                else
                    if [[ "$n" -eq 0 ]]; then ok "LoFTEE on $ASSEMBLY (rust): no splice-prediction parameter"; else bad "LoFTEE on $ASSEMBLY (rust): carries run_splice_predictions: $arg"; fi
                fi
            fi
            if [[ "$plugin" == LOEUF || "$plugin" == DosageSensitivity ]]; then
                # Each engine takes the word a run of plugin_gate_runs.toml gives it, the data root
                # apart: the runs name the file under {plugin_data}/<scope>/, the harness under /pd/.
                scope="$ASSEMBLY_LC"; [[ "$plugin" == DosageSensitivity ]] && scope=assembly_agnostic
                want="${arg%%/pd/*}{plugin_data}/$scope/${arg#*/pd/}"
                if grep -qF "\"$want\"" "$runs"; then ok "$plugin on $ASSEMBLY ($target): the word of the plugin gate's runs"; else bad "$plugin on $ASSEMBLY ($target): no run of plugin_gate_runs.toml carries '$want'"; fi
            fi
            missing=""
            # every /pd/<path> token the harness emits must be a file the layout stages
            for tok in $(echo "$arg" | tr ',' '\n' | grep -o '/pd/[^,]*' | sed 's#^/pd/##; s#/loftee_data/#/#'); do
                grep -qx "$tok" <<<"$staged" || missing="$missing $tok"
            done
            if [[ -z "$missing" ]]; then ok "$plugin on $ASSEMBLY ($target): every path is a staged file"; else bad "$plugin on $ASSEMBLY ($target): not staged:$missing"; fi
        done
    done
done

# A plugin without a data file on an assembly is a failure of plugin_staged_files, not an empty
# list that reads as nothing to stage.
if ! out=$(plugin_staged_files Conservation grch37 2>&1) && [[ "$out" == "ERROR: [plugin_data_layout] Conservation has no data file on grch37" ]]; then
    ok "plugin_staged_files Conservation on GRCh37 prints the layout error and exits 1"
else
    bad "plugin_staged_files Conservation on GRCh37 did not fail with the layout error: $out"
fi

# LOEUF reads one tabix table per assembly and DosageSensitivity one whole-read table for both,
# under the names the plugin gate's runs give the two engines below their data root: a name
# edited in the layout alone would stage a file that no run of plugin_gate_runs.toml reads.
for asm in grch37 grch38; do
    data=$(plugin_data_file LOEUF file "$asm" 2>/dev/null) || data=""
    staged=$(plugin_staged_files LOEUF "$asm" 2>&1)
    if [[ -n "$data" && "$staged" == "$(printf '%s\n%s.tbi' "$data" "$data")" ]]; then ok "LOEUF on $asm: stages $data and its index"; else bad "LOEUF on $asm: stages '$(tr '\n' ' ' <<<"$staged")', expected the table and its .tbi"; fi
    if [[ -n "$data" ]] && grep -qF "\"LOEUF,file={plugin_data}/$asm/$data,match_by=gene\"" "$runs"; then ok "LOEUF on $asm: the plugin gate's runs name $data"; else bad "LOEUF on $asm: no run of plugin_gate_runs.toml names '$data' under {plugin_data}/$asm/"; fi
    data=$(plugin_data_file DosageSensitivity file "$asm" 2>/dev/null) || data=""
    staged=$(plugin_staged_files DosageSensitivity "$asm" 2>&1)
    if [[ -n "$data" && "$staged" == "$data" ]]; then ok "DosageSensitivity on $asm: stages $data alone"; else bad "DosageSensitivity on $asm: stages '$(tr '\n' ' ' <<<"$staged")', expected the table alone"; fi
    if [[ -n "$data" ]] && grep -qF "\"DosageSensitivity,file={plugin_data}/assembly_agnostic/$data\"" "$runs"; then ok "DosageSensitivity on $asm: the plugin gate's runs name $data"; else bad "DosageSensitivity on $asm: no run of plugin_gate_runs.toml names '$data' under {plugin_data}/assembly_agnostic/"; fi
done
if a=$(plugin_data_file LOEUF file grch37 2>/dev/null) && b=$(plugin_data_file LOEUF file grch38 2>/dev/null) && [[ "$a" != "$b" ]]; then ok "LOEUF: one table per assembly"; else bad "LOEUF: the two assemblies do not have a table each"; fi
if a=$(plugin_data_file DosageSensitivity file grch37 2>/dev/null) && b=$(plugin_data_file DosageSensitivity file grch38 2>/dev/null) && [[ "$a" == "$b" ]]; then ok "DosageSensitivity: one table for both assemblies"; else bad "DosageSensitivity: the two assemblies do not share one table"; fi

# The data-preparation summary reads every listed file through verify_tabix_index: a fully
# prepared tree reads READY for every plugin and exits 0, and a tree missing one tabix index reads
# MISSING for that plugin and exits 1, so the index check neither demands an index of an
# index file nor overlooks a missing one. The files are empty: the summary checks presence
# and the index pairing alone under --source upstream, which fetches nothing but requires
# curl, bgzip and tabix on PATH, so stubs stand in for them.
setup="$ROOT/scripts/data/setup_plugin_data.sh"
stubs=$(mktemp -d "${TMPDIR:-/tmp}/plugin_data_layout_tools.XXXXXX")
for tool in curl bgzip tabix; do printf '#!/bin/sh\nexit 0\n' >"$stubs/$tool"; chmod +x "$stubs/$tool"; done
for ASSEMBLY_LC in grch37 grch38; do
    ASSEMBLY=$(echo "$ASSEMBLY_LC" | sed 's/grch/GRCh/')
    tree=$(mktemp -d "${TMPDIR:-/tmp}/plugin_data_layout.XXXXXX")
    for plugin in REVEL gnomADc AlphaMissense LoFtool pLI GWAS LoFTEE AncestralAllele OpenTargets Conservation LOEUF DosageSensitivity; do
        [[ ( "$plugin" == GWAS || "$plugin" == OpenTargets || "$plugin" == Conservation ) && "$ASSEMBLY_LC" == grch37 ]] && continue
        while IFS= read -r rel; do
            [[ -n "$rel" ]] || continue
            mkdir -p "$tree/$(dirname "$rel")" && : >"$tree/$rel"
        done < <(plugin_staged_files "$plugin" "$ASSEMBLY_LC")
    done
    if out=$(PATH="$stubs:$PATH" bash "$setup" --assembly "$ASSEMBLY" --data-dir "$tree" --source upstream 2>&1) && ! grep -q 'MISSING' <<<"$out"; then
        ok "data-preparation summary on $ASSEMBLY: every plugin READY on a fully prepared tree"
    else
        bad "data-preparation summary on $ASSEMBLY: a fully prepared tree is not READY: $(grep -E 'MISSING|\[warn\]|ERROR' <<<"$out" | head -3 | tr '\n' ' ')"
    fi
    # The DosageSensitivity table is read whole and staged with no index beside it; the LOEUF
    # table is read through tabix, so a tree without its index is not READY.
    if grep -q 'DosageSensitivity: READY' <<<"$out" && [[ ! -e "$tree/$(plugin_data_file DosageSensitivity file "$ASSEMBLY_LC").tbi" ]]; then
        ok "data-preparation summary on $ASSEMBLY: DosageSensitivity is READY on its table alone"
    else
        bad "data-preparation summary on $ASSEMBLY: DosageSensitivity is not READY on its table alone: $(grep -E 'DosageSensitivity|dosagesensitivity' <<<"$out" | head -2 | tr '\n' ' ')"
    fi
    rm -f "$tree/$(plugin_data_file LOEUF file "$ASSEMBLY_LC").tbi"
    if ! out=$(PATH="$stubs:$PATH" bash "$setup" --assembly "$ASSEMBLY" --data-dir "$tree" --source upstream --plugins LOEUF,DosageSensitivity 2>&1) && grep -q 'LOEUF: MISSING' <<<"$out" && grep -q 'DosageSensitivity: READY' <<<"$out"; then
        ok "data-preparation summary on $ASSEMBLY: LOEUF without its tabix index reads MISSING beside a READY DosageSensitivity"
    else
        bad "data-preparation summary on $ASSEMBLY: LOEUF without its tabix index was not reported: $(grep -E 'LOEUF|DosageSensitivity' <<<"$out" | tail -2 | tr '\n' ' ')"
    fi
    : >"$tree/$(plugin_data_file LOEUF file "$ASSEMBLY_LC").tbi"
    rm -f "$tree/$(plugin_data_file REVEL file "$ASSEMBLY_LC").tbi"
    if ! PATH="$stubs:$PATH" bash "$setup" --assembly "$ASSEMBLY" --data-dir "$tree" --source upstream >"$tree/summary.out" 2>&1 && grep -q 'REVEL: MISSING' "$tree/summary.out"; then
        ok "data-preparation summary on $ASSEMBLY: a missing tabix index reads MISSING and exits 1"
    else
        bad "data-preparation summary on $ASSEMBLY: a missing tabix index was not reported"
    fi
    rm -rf "$tree"
done
rm -rf "$stubs"

# A bucket copy holds the DosageSensitivity table once, in the shared scope, and a LOEUF table
# in each assembly's scope: setup_plugin_data.sh's s3_scope_for, extracted by name.
scope_start=$(grep -n '^s3_scope_for()' "$setup" | cut -d: -f1)
scope_end=$(awk -v s="$scope_start" 'NR>s && /^}/{print NR; exit}' "$setup")
[[ -n "$scope_start" && -n "$scope_end" ]] || { echo "ERROR: [test_plugin_data_layout] s3_scope_for not found in $setup" >&2; exit 1; }
eval "$(sed -n "${scope_start},${scope_end}p" "$setup")"
for ASSEMBLY_LC in grch37 grch38; do
    if [[ "$(s3_scope_for DosageSensitivity)" == assembly_agnostic ]]; then ok "DosageSensitivity on $ASSEMBLY_LC: the shared scope"; else bad "DosageSensitivity on $ASSEMBLY_LC: scope '$(s3_scope_for DosageSensitivity)', expected assembly_agnostic"; fi
    if [[ "$(s3_scope_for LOEUF)" == "$ASSEMBLY_LC" ]]; then ok "LOEUF on $ASSEMBLY_LC: the assembly's own scope"; else bad "LOEUF on $ASSEMBLY_LC: scope '$(s3_scope_for LOEUF)', expected $ASSEMBLY_LC"; fi
done

# The harness names no flat plugin data file; every path of the plugin data directory comes
# from the layout table.
if grep -qE '/revel\.tsv\.gz|/gnomad_coverage\.tsv\.gz' "$harness"; then bad "run_concordance.sh names a flat plugin data file"; else ok "run_concordance.sh names no flat plugin data file"; fi

# The note setup_plugin_data.sh prints for the GRCh37 LoFTEE bundle names the key the harness
# passes to vep-rs there, without which vep-rs refuses a GRCh37 run.
if grep -q 'Pass it as gerp_tabix= with run_splice_predictions=0' "$setup"; then ok "setup_plugin_data.sh's GRCh37 LoFTEE note names run_splice_predictions=0"; else bad "setup_plugin_data.sh's GRCh37 LoFTEE note does not name run_splice_predictions=0"; fi

# CADD, SpliceAI, AVI, dbNSFP and dbscSNV read the user's own file: the harness passes both
# engines the path its variable names, under Docker the mount of that file's directory, and
# names the variables a plugin cannot run without.
user_file_case() {
    local plugin="$1" want="$2" docker_want="$3" got
    for target in rust perl; do
        got=$(build_plugin_args "$plugin" /pd "$target") || got="(failed)"
        if [[ "$got" == "$want" ]]; then ok "$plugin ($target): the user's own file"; else bad "$plugin ($target): '$got', expected '$want'"; fi
    done
    got=$(PERL_VEP_MODE=docker build_plugin_args "$plugin" /pd perl) || got="(failed)"
    if [[ "$got" == "$docker_want" ]]; then ok "$plugin (perl, docker): the mount of the user's directory"; else bad "$plugin (perl, docker): '$got', expected '$docker_want'"; fi
}
vars_start=$(grep -n '^user_plugin_file_vars()' "$harness" | cut -d: -f1)
vars_end=$(awk -v s="$vars_start" 'NR>s && /^}/{print NR; exit}' "$harness")
[[ -n "$vars_start" && -n "$vars_end" ]] || { echo "ERROR: [test_plugin_data_layout] user_plugin_file_vars not found in $harness" >&2; exit 1; }
eval "$(sed -n "${vars_start},${vars_end}p" "$harness")"
ASSEMBLY=GRCh38 ASSEMBLY_LC=grch38
# shellcheck disable=SC2034  # build_plugin_args reads them by name
CADD_SNV_FILE=/u/c/snv_synthetic.tsv.gz CADD_INDEL_FILE=/u/c/indel_synthetic.tsv.gz SPLICEAI_SNV_FILE=/u/s/snv.synthetic.vcf.gz AVI_FILE=/u/a/avi_synthetic.tsv.gz DBNSFP_FILE=/u/d/dbNSFP4.9a_grch38_synthetic.gz DBSCSNV_FILE=/u/b/dbscsnv_synthetic.txt.gz
user_file_case CADD "CADD,snv=/u/c/snv_synthetic.tsv.gz,indels=/u/c/indel_synthetic.tsv.gz" "CADD,snv=/work/user_plugin_data/CADD_SNV_FILE/snv_synthetic.tsv.gz,indels=/work/user_plugin_data/CADD_INDEL_FILE/indel_synthetic.tsv.gz"
user_file_case SpliceAI "SpliceAI,snv=/u/s/snv.synthetic.vcf.gz" "SpliceAI,snv=/work/user_plugin_data/SPLICEAI_SNV_FILE/snv.synthetic.vcf.gz"
SPLICEAI_INDEL_FILE=/u/s/indel.synthetic.vcf.gz user_file_case SpliceAI "SpliceAI,snv=/u/s/snv.synthetic.vcf.gz,indel=/u/s/indel.synthetic.vcf.gz" "SpliceAI,snv=/work/user_plugin_data/SPLICEAI_SNV_FILE/snv.synthetic.vcf.gz,indel=/work/user_plugin_data/SPLICEAI_INDEL_FILE/indel.synthetic.vcf.gz"
user_file_case AVI "AVI,file=/u/a/avi_synthetic.tsv.gz,raw=1,links=1" "AVI,file=/work/user_plugin_data/AVI_FILE/avi_synthetic.tsv.gz,raw=1,links=1"
user_file_case dbNSFP "dbNSFP,/u/d/dbNSFP4.9a_grch38_synthetic.gz,SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred" "dbNSFP,/work/user_plugin_data/DBNSFP_FILE/dbNSFP4.9a_grch38_synthetic.gz,SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred"
user_file_case dbscSNV "dbscSNV,/u/b/dbscsnv_synthetic.txt.gz" "dbscSNV,/work/user_plugin_data/DBSCSNV_FILE/dbscsnv_synthetic.txt.gz"
if [[ "$(user_plugin_file_vars CADD)" == "CADD_SNV_FILE CADD_INDEL_FILE" && "$(user_plugin_file_vars SpliceAI)" == "SPLICEAI_SNV_FILE" && -z "$(user_plugin_file_vars REVEL)" ]]; then ok "the variables a plugin cannot run without: both CADD files, the SpliceAI SNV file, none for a plugin of the layout"; else bad "user_plugin_file_vars names other variables"; fi
unset CADD_SNV_FILE CADD_INDEL_FILE SPLICEAI_SNV_FILE AVI_FILE DBNSFP_FILE DBSCSNV_FILE

# The Perl output cache is keyed on PERL_VEP_SIGNATURE, evaluated here as the harness writes
# it: the same plugin names with a changed parameter must not be served an older output.
sig_line=$(grep -m1 '^PERL_VEP_SIGNATURE=' "$harness")
[[ -n "$sig_line" ]] || { echo "ERROR: [test_plugin_data_layout] PERL_VEP_SIGNATURE= not found in $harness" >&2; exit 1; }
# shellcheck disable=SC2034  # the eval of the harness line reads them
perl_signature() {
    (
        set -u
        PERL_VEP_MODE=native PERL_VEP_DOCKER_IMAGE=img PERL_CACHE_DIR=/pc SPECIES=homo_sapiens ASSEMBLY=GRCh38
        PLUGINS=LoFTEE PLUGIN_DATA_DIR=/pd PLUGIN_DATA_FINGERPRINT=fp FASTA=/ref.fa
        PERL_PLUGIN_ARGS=("$@")
        eval "$sig_line"
        echo "$PERL_VEP_SIGNATURE"
    )
}
if [[ "$(perl_signature --plugin 'LoF,run_splice_predictions:0')" != "$(perl_signature --plugin 'LoF,run_splice_predictions:1')" ]]; then ok "the Perl output-cache signature keys on the plugin argument strings"; else bad "the Perl output-cache signature ignores the plugin argument strings"; fi

# provenance.json records the same plugin argument strings the signature keys on, evaluated
# here as the harness writes it: the heredoc that writes reports/provenance.json, run whole
# against a scratch work directory.
prov_out=$(grep -n 'out = Path("${WORK_DIR}/reports/provenance.json")' "$harness" | cut -d: -f1)
[[ -n "$prov_out" ]] || { echo "ERROR: [test_plugin_data_layout] the provenance.json writer was not found in $harness" >&2; exit 1; }
prov_start=$(awk -v e="$prov_out" 'NR<e && /^python3 - <<PY$/{s=NR} END{print s}' "$harness")
prov_end=$(awk -v s="$prov_out" 'NR>s && /^PY$/{print NR; exit}' "$harness")
[[ -n "$prov_start" && -n "$prov_end" ]] || { echo "ERROR: [test_plugin_data_layout] the provenance.json heredoc was not delimited in $harness" >&2; exit 1; }
prov_block=$(sed -n "${prov_start},${prov_end}p" "$harness")
# shellcheck disable=SC2034  # the eval of the harness heredoc reads them
provenance_plugin_args() {
    (
        # A label the heredoc interpolates only inside a string expands to the empty string.
        set +u
        WORK_DIR=$(mktemp -d)
        mkdir -p "$WORK_DIR/reports" && : >"$WORK_DIR/reports/_input_shas.tsv"
        BENCHMARK_DIR=/bd ENSEMBL_VEP_DIR=/ev VEP_RS_DIR=/vr PERL_CACHE_DIR=/pc JSON_CACHE_DIR=/jc FASTA=/ref.fa
        PERL_VEP_MODE=native PERL_VEP_DOCKER_IMAGE=img FIXTURES_DIR=/fx MODE=smoke SMOKE_VARIANTS=5
        SPECIES=homo_sapiens ASSEMBLY=GRCh38 RUST_VERSION=116 RUST_RELEASE=1 PERL_CACHE_VERSION=116
        REFERENCE_RELEASE=116.2 PLUGINS="$1" PLUGIN_DATA_DIR=/pd PLUGIN_TIER=tier1 NO_FASTA=0
        CANONICAL_CONTIGS=1 NO_CANONICAL_CONTIGS=0 RUN_ID="" HOST_ID=""
        shift
        PERL_PLUGIN_ARGS=("$@")
        eval "$prov_block" >/dev/null || { echo "ERROR: [test_plugin_data_layout] the provenance.json heredoc failed to run" >&2; exit 1; }
        python3 -c 'import json, sys; print(json.dumps(json.load(open(sys.argv[1]))["perl_plugin_args"]))' "$WORK_DIR/reports/provenance.json"
        rm -rf "$WORK_DIR"
    )
}
sig_args="$(perl_signature --plugin 'LoF,run_splice_predictions:0' --plugin 'REVEL,file=/pd/revel.tsv.gz')"
sig_args="${sig_args##*;perl_plugin_args=}"
prov_args="$(provenance_plugin_args LoFTEE,REVEL --plugin 'LoF,run_splice_predictions:0' --plugin 'REVEL,file=/pd/revel.tsv.gz')"
if [[ "$prov_args" == "\"$sig_args\"" ]]; then ok "provenance.json records the plugin argument strings the signature keys on"; else bad "provenance.json perl_plugin_args $prov_args differs from the signature's $sig_args"; fi
prov_none="$(provenance_plugin_args "")"
if [[ "$prov_none" == "null" ]]; then ok "provenance.json records null plugin argument strings on a run without plugins"; else bad "provenance.json perl_plugin_args on a run without plugins: $prov_none"; fi

echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
