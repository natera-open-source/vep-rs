#!/usr/bin/env bash
# Every data path run_concordance.sh passes to an engine must be a file
# setup_plugin_data.sh stages, on both assemblies. Both scripts read the table in
# plugin_data_layout.sh; this test extracts the harness's build_plugin_args function and
# checks its output against plugin_staged_files, so a name edited in one place without the
# other fails here rather than at run time as an engine reading a missing file.
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

echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
