#!/usr/bin/env bash
# Every data path run_concordance.sh passes to an engine must be a file
# setup_plugin_data.sh stages, on both assemblies. Both scripts read the table in
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

grep -q 'source "${SCRIPT_DIR}/../data/plugin_data_layout.sh"' "$harness" && ok "run_concordance.sh sources the shared layout" || bad "run_concordance.sh does not source plugin_data_layout.sh"
grep -q 'plugin_staged_files "$1" "${ASSEMBLY_LC}"' "$ROOT/scripts/data/setup_plugin_data.sh" && ok "setup_plugin_data.sh stages from the shared layout" || bad "setup_plugin_data.sh has its own file table"

PERL_VEP_MODE=native
for ASSEMBLY_LC in grch37 grch38; do
    ASSEMBLY=$(echo "$ASSEMBLY_LC" | sed 's/grch/GRCh/')
    export ASSEMBLY ASSEMBLY_LC
    for plugin in CADD REVEL SpliceAI gnomADc AlphaMissense dbscSNV LoFtool pLI GWAS LoFTEE; do
        [[ "$plugin" == GWAS && "$ASSEMBLY_LC" == grch37 ]] && continue
        staged=$(plugin_staged_files "$plugin" "$ASSEMBLY_LC")
        for target in rust perl; do
            arg=$(LOFTEE_PERL_PATH=/lp build_plugin_args "$plugin" /pd "$target") || { bad "$plugin/$ASSEMBLY_LC/$target: build_plugin_args failed"; continue; }
            if [[ "$plugin" == LoFTEE ]]; then
                # Perl's master branch defaults run_splice_predictions to 1 (LoF.pm:98) and vep-rs has
                # no splice layer: the Perl string pins it off once, the vep-rs string never names it.
                n=$(grep -o 'run_splice_predictions' <<<"$arg" | grep -c .) || n=0
                if [[ "$target" == perl ]]; then
                    if [[ "$arg" == *,run_splice_predictions:0 && "$n" -eq 1 ]]; then ok "LoFTEE on $ASSEMBLY (perl): the splice-prediction layer is pinned off"; else bad "LoFTEE on $ASSEMBLY (perl): run_splice_predictions:0 is not the single last parameter: $arg"; fi
                else
                    if [[ "$n" -eq 0 ]]; then ok "LoFTEE on $ASSEMBLY (rust): no splice-prediction parameter"; else bad "LoFTEE on $ASSEMBLY (rust): carries run_splice_predictions: $arg"; fi
                fi
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

# The harness names no flat plugin data file; every path comes from the layout table.
if grep -qE 'cadd_snvs\.tsv\.gz|spliceai_snv\.vcf\.gz|/revel\.tsv\.gz|/dbnsfp\.gz' "$harness"; then bad "run_concordance.sh still names a flat plugin data file"; else ok "no flat plugin data name survives in run_concordance.sh"; fi

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

echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
