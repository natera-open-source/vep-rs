#!/usr/bin/env bash
# Every data path run_concordance.sh passes to an engine must be a file
# setup_plugin_data.sh stages, on both assemblies, or the user's own file where a plugin
# reads one. Both scripts read the table in
# plugin_data_layout.sh; this test extracts the harness's build_plugin_args function and
# checks its output against plugin_staged_files, so a name edited in one place without the
# other fails here rather than at run time as an engine reading a missing file. Every plugin
# data path a whole run of plugin_gate_runs.toml names must be a staged file too. It also
# evaluates the harness's PERL_VEP_SIGNATURE line, the Ensembl VEP output cache's key, which must
# change with the plugin argument strings. An assembly that is neither human one must be
# refused by the layout, by the harness's argument builder and by setup_plugin_data.sh. The
# harness's plugin_docker_mounts must mount the directory of a user's own file for a run of
# the plugin that reads it alone.
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
                # no splice layer and refuses a GRCh37 run that leaves it on: the Ensembl VEP string
                # pins it off once, the vep-rs string once on GRCh37 and never names it on GRCh38.
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

# refused <name> <expected stderr> <command...>: the command prints nothing on standard
# output, exactly <expected stderr> on standard error, and returns 1.
refused() {
    local name="$1" want="$2" out err status
    shift 2
    err=$(mktemp "${TMPDIR:-/tmp}/plugin_data_layout_stderr.XXXXXX")
    out=$("$@" 2>"$err")
    status=$?
    if [[ "$status" -eq 1 && -z "$out" && "$(cat "$err")" == "$want" ]]; then ok "$name"; else bad "$name: status $status, printed '$out', stderr '$(cat "$err")'"; fi
    rm -f "$err"
}
# gives <name> <expected stdout> <command...>: the command prints exactly <expected stdout>,
# nothing on standard error, and returns 0.
gives() {
    local name="$1" want="$2" out err status
    shift 2
    err=$(mktemp "${TMPDIR:-/tmp}/plugin_data_layout_stderr.XXXXXX")
    out=$("$@" 2>"$err")
    status=$?
    if [[ "$status" -eq 0 && "$out" == "$want" && ! -s "$err" ]]; then ok "$name"; else bad "$name: status $status, printed '$out', stderr '$(cat "$err")'"; fi
    rm -f "$err"
}
# args_on <assembly> <plugin> <data dir> <target>: build_plugin_args under an assembly, lower-cased
# as the harness lower-cases --assembly, with LoFTEE's data under the plugin data directory.
args_on() {
    # shellcheck disable=SC2034  # build_plugin_args reads them by name
    local ASSEMBLY="$1" ASSEMBLY_LC
    ASSEMBLY_LC=$(echo "$1" | tr '[:upper:]' '[:lower:]')
    shift
    LOFTEE_DATA_DIR="" LOFTEE_PERL_PATH=/lp build_plugin_args "$@"
}

# LOEUF and LoFTEE have a data bundle per human assembly and none for another genome: an
# assembly that is neither is refused by the layout and by the harness's argument builder, where
# reading it as GRCh38 would hand an engine the GRCh38 file as that genome's.
refused "LOEUF on a made-up assembly: the layout names no table" "ERROR: [plugin_data_layout] LOEUF has no data file on xyz1" plugin_data_file LOEUF file xyz1
for role in ancestor ancestor_bgzf gerp; do
    refused "LoFTEE $role on a made-up assembly: the layout names no file" "ERROR: [plugin_data_layout] LoFTEE has no data file on xyz1" plugin_data_file LoFTEE "$role" xyz1
done
for plugin in LOEUF LoFTEE; do
    refused "plugin_staged_files $plugin on a made-up assembly: the layout error, no file listed" "ERROR: [plugin_data_layout] $plugin has no data file on xyz1" plugin_staged_files "$plugin" xyz1
    for target in rust perl; do
        refused "$plugin on a made-up assembly ($target): the harness builds no argument string" "ERROR: [plugin_data_layout] $plugin has no data file on xyz1" args_on Xyz1 "$plugin" /pd "$target"
    done
done
# The layout's assembly is lower case. GRCh37 in the --assembly spelling is not grch37 and is
# refused; reading it as "not grch37" would name the GRCh38 file.
refused "LOEUF on GRCh37 in the --assembly spelling: refused, not given the GRCh38 table" "ERROR: [plugin_data_layout] LOEUF has no data file on GRCh37" plugin_data_file LOEUF file GRCh37
refused "LoFTEE gerp on GRCh37 in the --assembly spelling: refused, not given the GRCh38 bigWig" "ERROR: [plugin_data_layout] LoFTEE has no data file on GRCh37" plugin_data_file LoFTEE gerp GRCh37
# The plugin gate's GRCh37 run of LoFTEE reads made-up conservation data, which the layout names
# on grch37 alone: on any other assembly each of its three roles is refused by name.
for role in synthetic_gerp synthetic_conservation synthetic_ancestor_bgzf; do
    refused "LoFTEE $role on grch38: the layout names no synthetic conservation file" "ERROR: [plugin_data_layout] LoFTEE has no synthetic conservation file on grch38" plugin_data_file LoFTEE "$role" grch38
    refused "LoFTEE $role on a made-up assembly: the layout names no synthetic conservation file" "ERROR: [plugin_data_layout] LoFTEE has no synthetic conservation file on xyz1" plugin_data_file LoFTEE "$role" xyz1
done

# On the two human assemblies the layout names LOEUF's table per assembly and LoFTEE's GERP file
# as each branch of the plugin reads it, a table with its tabix index on GRCh37 (master's LoF.pm,
# lines 138 to 141) and a bigWig on GRCh38, and the harness passes both engines those files.
gives "LOEUF on grch37: the gnomAD v2.1.1 table" "loeuf/loeuf_dataset.tsv.gz" plugin_data_file LOEUF file grch37
gives "LOEUF on grch38: Ensembl's table" "loeuf/loeuf_dataset_grch38.tsv.gz" plugin_data_file LOEUF file grch38
gives "LoFTEE gerp on grch37: the tabix table" "loftee/grch37/GERP_scores.final.sorted.txt.gz" plugin_data_file LoFTEE gerp grch37
gives "LoFTEE gerp on grch38: the bigWig" "loftee/grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw" plugin_data_file LoFTEE gerp grch38
for asm in grch37 grch38; do
    gives "LoFTEE ancestor on $asm: the plain FASTA" "loftee/$asm/human_ancestor.fa" plugin_data_file LoFTEE ancestor "$asm"
    gives "LoFTEE ancestor_bgzf on $asm: the bgzipped FASTA" "loftee/$asm/human_ancestor.fa.gz" plugin_data_file LoFTEE ancestor_bgzf "$asm"
done
ancestors() { printf 'loftee/%s/human_ancestor.fa\nloftee/%s/human_ancestor.fa.fai\nloftee/%s/human_ancestor.fa.gz\nloftee/%s/human_ancestor.fa.gz.fai\nloftee/%s/human_ancestor.fa.gz.gzi' "$1" "$1" "$1" "$1" "$1"; }
# The six files of the plugin gate's GRCh37 run of LoFTEE, in one directory: the made-up GERP
# table with its .tbi, the made-up conservation database, and the ancestral FASTA with its
# .fai and .gzi.
synthetic_run() { printf 'loftee_synthetic_conservation/GERP_scores.synthetic.txt.gz\nloftee_synthetic_conservation/GERP_scores.synthetic.txt.gz.tbi\nloftee_synthetic_conservation/phylocsf_gerp.synthetic.sql\nloftee_synthetic_conservation/human_ancestor.fa.gz\nloftee_synthetic_conservation/human_ancestor.fa.gz.fai\nloftee_synthetic_conservation/human_ancestor.fa.gz.gzi'; }
gives "plugin_staged_files LoFTEE on grch37: both FASTA files with their indexes, the GERP table and its .tbi, then the six files of the synthetic conservation run" "$(ancestors grch37)"$'\n'"loftee/grch37/GERP_scores.final.sorted.txt.gz"$'\n'"loftee/grch37/GERP_scores.final.sorted.txt.gz.tbi"$'\n'"$(synthetic_run)" plugin_staged_files LoFTEE grch37
gives "plugin_staged_files LoFTEE on grch38: both FASTA files with their indexes and the bigWig alone" "$(ancestors grch38)"$'\n'"loftee/grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw" plugin_staged_files LoFTEE grch38
gives "LoFTEE synthetic_gerp on grch37: the made-up GERP table" "loftee_synthetic_conservation/GERP_scores.synthetic.txt.gz" plugin_data_file LoFTEE synthetic_gerp grch37
gives "LoFTEE synthetic_conservation on grch37: the made-up conservation database" "loftee_synthetic_conservation/phylocsf_gerp.synthetic.sql" plugin_data_file LoFTEE synthetic_conservation grch37
gives "LoFTEE synthetic_ancestor_bgzf on grch37: the run's ancestral FASTA" "loftee_synthetic_conservation/human_ancestor.fa.gz" plugin_data_file LoFTEE synthetic_ancestor_bgzf grch37
gives "plugin_synthetic_run_files LoFTEE on grch37: the six files of the synthetic conservation run" "$(synthetic_run)" plugin_synthetic_run_files LoFTEE grch37
gives "plugin_synthetic_run_files LoFTEE on grch38: no file" "" plugin_synthetic_run_files LoFTEE grch38
gives "plugin_synthetic_run_files REVEL on grch37: no file" "" plugin_synthetic_run_files REVEL grch37
for target in rust perl; do
    gives "LOEUF on GRCh37 ($target): the GRCh37 table by gene" "LOEUF,file=/pd/loeuf/loeuf_dataset.tsv.gz,match_by=gene" args_on GRCh37 LOEUF /pd "$target"
    gives "LOEUF on GRCh38 ($target): the GRCh38 table by gene" "LOEUF,file=/pd/loeuf/loeuf_dataset_grch38.tsv.gz,match_by=gene" args_on GRCh38 LOEUF /pd "$target"
done
gives "LoFTEE on GRCh37 (rust): the plain FASTA and the tabix table" "LoFTEE,human_ancestor_fa=/pd/loftee/grch37/human_ancestor.fa,gerp_tabix=/pd/loftee/grch37/GERP_scores.final.sorted.txt.gz,run_splice_predictions=0" args_on GRCh37 LoFTEE /pd rust
gives "LoFTEE on GRCh37 (perl): the bgzipped FASTA and the tabix table" "LoF,loftee_path:/lp,human_ancestor_fa:/pd/loftee/grch37/human_ancestor.fa.gz,gerp_file:/pd/loftee/grch37/GERP_scores.final.sorted.txt.gz,run_splice_predictions:0" args_on GRCh37 LoFTEE /pd perl
gives "LoFTEE on GRCh38 (rust): the plain FASTA and the bigWig" "LoFTEE,human_ancestor_fa=/pd/loftee/grch38/human_ancestor.fa,gerp_bigwig=/pd/loftee/grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw" args_on GRCh38 LoFTEE /pd rust
gives "LoFTEE on GRCh38 (perl): the bgzipped FASTA and the bigWig" "LoF,loftee_path:/lp,human_ancestor_fa:/pd/loftee/grch38/human_ancestor.fa.gz,gerp_bigwig:/pd/loftee/grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw,run_splice_predictions:0" args_on GRCh38 LoFTEE /pd perl

# gnomADMt reads the mitochondrial sites VCF with its tabix index, a file its own documentation
# calls GRCh38 based (gnomADMt.pm:50), and UTRAnnotator one uORF table per assembly, read
# whole (UTRAnnotator.pm:42-43): the layout names each under the path the plugin gate's runs
# give both engines, and refuses an assembly the file is not published for.
gives "gnomADMt on grch38: the mitochondrial sites VCF" "gnomadmt/gnomad.genomes.v3.1.sites.chrM.vcf.bgz" plugin_data_file gnomADMt file grch38
gives "plugin_staged_files gnomADMt on grch38: the VCF and its .tbi" "gnomadmt/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"$'\n'"gnomadmt/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi" plugin_staged_files gnomADMt grch38
refused "gnomADMt on grch37: the layout names no file" "ERROR: [plugin_data_layout] gnomADMt has no data file on grch37" plugin_data_file gnomADMt file grch37
refused "plugin_staged_files gnomADMt on grch37: the layout error, no file listed" "ERROR: [plugin_data_layout] gnomADMt has no data file on grch37" plugin_staged_files gnomADMt grch37
gives "UTRAnnotator on grch37: the GRCh37 uORF table" "utrannotator/uORF_5UTR_GRCh37_PUBLIC.txt" plugin_data_file UTRAnnotator file grch37
gives "UTRAnnotator on grch38: the GRCh38 uORF table" "utrannotator/uORF_5UTR_GRCh38_PUBLIC.txt" plugin_data_file UTRAnnotator file grch38
gives "plugin_staged_files UTRAnnotator on grch37: the table alone" "utrannotator/uORF_5UTR_GRCh37_PUBLIC.txt" plugin_staged_files UTRAnnotator grch37
gives "plugin_staged_files UTRAnnotator on grch38: the table alone" "utrannotator/uORF_5UTR_GRCh38_PUBLIC.txt" plugin_staged_files UTRAnnotator grch38
refused "UTRAnnotator on a made-up assembly: the layout names no table" "ERROR: [plugin_data_layout] UTRAnnotator has no data file on xyz1" plugin_data_file UTRAnnotator file xyz1
refused "UTRAnnotator on GRCh37 in the --assembly spelling: refused, not given the GRCh38 table" "ERROR: [plugin_data_layout] UTRAnnotator has no data file on GRCh37" plugin_data_file UTRAnnotator file GRCh37
for word in "gnomADMt,file={plugin_data}/assembly_agnostic/$(plugin_data_file gnomADMt file grch38)" "UTRAnnotator,file={plugin_data}/grch37/$(plugin_data_file UTRAnnotator file grch37)" "UTRAnnotator,file={plugin_data}/grch38/$(plugin_data_file UTRAnnotator file grch38)"; do
    if [[ "$(grep -cF "\"$word\"" "$runs")" -ge 1 && "$(grep -oF "\"$word\"" "$runs" | grep -c .)" -eq 2 ]]; then ok "${word%%,*}: the plugin gate's runs give both engines $word"; else bad "${word%%,*}: the runs of plugin_gate_runs.toml do not give both engines '$word'"; fi
done

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
    for plugin in REVEL gnomADc AlphaMissense LoFtool pLI GWAS LoFTEE AncestralAllele OpenTargets Conservation LOEUF DosageSensitivity UTRAnnotator gnomADMt; do
        [[ ( "$plugin" == GWAS || "$plugin" == OpenTargets || "$plugin" == Conservation || "$plugin" == gnomADMt ) && "$ASSEMBLY_LC" == grch37 ]] && continue
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
    # The note printed for LoFTEE's GERP file on GRCh37, whole: how the plugin's master branch
    # reads the file (its LoF.pm, lines 138 to 141) and the key vep-rs takes it by.
    if [[ "$ASSEMBLY_LC" == grch37 ]]; then
        if grep -qxF "    [manual] On GRCh37 the plugin's master branch reads its GERP file through tabix (LoF.pm lines 138 to 141), not from a bigWig. Pass it as gerp_tabix= with run_splice_predictions=0, without which vep-rs refuses a GRCh37 run. The END_TRUNC cutoff is +180 here, vs -58 on grch38." <<<"$out"; then
            ok "setup_plugin_data.sh on GRCh37: the LoFTEE GERP note says how the plugin reads the file"
        else
            bad "setup_plugin_data.sh on GRCh37: the LoFTEE GERP note reads: $(grep -F 'gerp_tabix=' <<<"$out" | head -1)"
        fi
    fi
    # No source publishes the made-up files of LoFTEE's synthetic conservation run, so the summary
    # reads the run's files in a line of their own, on GRCh37 alone: READY when the six stand, and
    # ABSENT, with LoFTEE itself READY and exit status 0, without the directory or one index of it.
    # The note that names the run's files and their generator is printed on GRCh37 alone.
    run_note="    [manual] The plugin gate's GRCh37 run of LoFTEE reads a made-up GERP table and a made-up conservation database, which tools/vep-loftee-synthetic-conservation writes, with the ancestral FASTA beside them, in loftee_synthetic_conservation/. docs/plugins/reproducing-concordance.md gives the commands."
    if [[ "$ASSEMBLY_LC" == grch37 ]]; then
        if grep -qxF "  LoFTEE: READY" <<<"$out" && grep -qxF "  LoFTEE synthetic conservation run: READY" <<<"$out"; then
            ok "data-preparation summary on $ASSEMBLY: the synthetic conservation run of LoFTEE is READY on its six files"
        else
            bad "data-preparation summary on $ASSEMBLY: the synthetic conservation run is not READY: $(grep -F 'LoFTEE' <<<"$out" | tail -2 | tr '\n' ' ')"
        fi
        if grep -qxF "$run_note" <<<"$out"; then
            ok "setup_plugin_data.sh on $ASSEMBLY: the note on the synthetic conservation run names its files and their generator"
        else
            bad "setup_plugin_data.sh on $ASSEMBLY: the note on the synthetic conservation run reads: $(grep -F 'made-up' <<<"$out" | head -1)"
        fi
        mv "$tree/loftee_synthetic_conservation" "$tree.synthetic_run"
        if syn_out=$(PATH="$stubs:$PATH" bash "$setup" --assembly "$ASSEMBLY" --data-dir "$tree" --source upstream 2>&1) && ! grep -q 'MISSING' <<<"$syn_out" && grep -qxF "  LoFTEE: READY" <<<"$syn_out" && grep -qxF "  LoFTEE synthetic conservation run: ABSENT" <<<"$syn_out"; then
            ok "data-preparation summary on $ASSEMBLY: LoFTEE is READY on its seven published files, the synthetic conservation run ABSENT, exit 0"
        else
            bad "data-preparation summary on $ASSEMBLY: without the synthetic conservation run: $(grep -E 'LoFTEE|MISSING|ERROR' <<<"$syn_out" | tail -3 | tr '\n' ' ')"
        fi
        mv "$tree.synthetic_run" "$tree/loftee_synthetic_conservation"
        rm -f "$tree/$(plugin_data_file LoFTEE synthetic_gerp "$ASSEMBLY_LC").tbi"
        if syn_out=$(PATH="$stubs:$PATH" bash "$setup" --assembly "$ASSEMBLY" --data-dir "$tree" --source upstream 2>&1) && ! grep -q 'MISSING' <<<"$syn_out" && grep -qxF "  LoFTEE synthetic conservation run: ABSENT" <<<"$syn_out"; then
            ok "data-preparation summary on $ASSEMBLY: the synthetic conservation run without its table's tabix index is ABSENT, exit 0"
        else
            bad "data-preparation summary on $ASSEMBLY: the synthetic conservation run without its table's tabix index: $(grep -E 'LoFTEE|MISSING|ERROR' <<<"$syn_out" | tail -3 | tr '\n' ' ')"
        fi
        : >"$tree/$(plugin_data_file LoFTEE synthetic_gerp "$ASSEMBLY_LC").tbi"
    else
        if ! grep -qF 'synthetic conservation run' <<<"$out"; then
            ok "data-preparation summary on $ASSEMBLY: no line on a synthetic conservation run"
        else
            bad "data-preparation summary on $ASSEMBLY: names a synthetic conservation run: $(grep -F 'synthetic conservation run' <<<"$out" | head -1)"
        fi
        if ! grep -qF 'made-up' <<<"$out"; then
            ok "setup_plugin_data.sh on $ASSEMBLY: no note on a synthetic conservation run"
        else
            bad "setup_plugin_data.sh on $ASSEMBLY: prints a note on made-up data: $(grep -F 'made-up' <<<"$out" | head -1)"
        fi
    fi
    # UTRAnnotator's table and gnomADMt's VCF: the source line of each, whole, from the place
    # the plugin's own documentation names (UTRAnnotator.pm:42-43, gnomADMt.pm:51-52), and
    # gnomADMt skipped by name on GRCh37.
    if grep -qxF "    https://raw.githubusercontent.com/Ensembl/UTRannotator/40e30db003d72435eb3744ee747029294abbcf3a/uORF_5UTR_${ASSEMBLY}_PUBLIC.txt" <<<"$out" && grep -qxF "  UTRAnnotator: READY" <<<"$out"; then
        ok "setup_plugin_data.sh on $ASSEMBLY: prints the source of UTRAnnotator's table and reads it READY"
    else
        bad "setup_plugin_data.sh on $ASSEMBLY: UTRAnnotator's source line or READY line is missing: $(grep -E 'UTRannotator|UTRAnnotator' <<<"$out" | head -3 | tr '\n' ' ')"
    fi
    if [[ "$ASSEMBLY_LC" == grch38 ]]; then
        if grep -qxF "    https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz (+ .tbi)" <<<"$out" && grep -qxF "  gnomADMt: READY" <<<"$out"; then
            ok "setup_plugin_data.sh on $ASSEMBLY: prints the source of gnomADMt's VCF and reads it READY"
        else
            bad "setup_plugin_data.sh on $ASSEMBLY: gnomADMt's source line or READY line is missing: $(grep -E 'gnomADMt|chrM' <<<"$out" | head -3 | tr '\n' ' ')"
        fi
    else
        if grep -qxF "SKIP  gnomADMt: GRCh38 only (the gnomAD v3.1 mitochondrial sites file is GRCh38 based)" <<<"$out" && ! grep -qE '^\[gnomADMt\]$|^  gnomADMt: ' <<<"$out"; then
            ok "setup_plugin_data.sh on $ASSEMBLY: gnomADMt is skipped by name and stages nothing"
        else
            bad "setup_plugin_data.sh on $ASSEMBLY: gnomADMt is not skipped: $(grep -F 'gnomADMt' <<<"$out" | head -3 | tr '\n' ' ')"
        fi
    fi
    # The run_concordance.sh command the summary prints names only plugins that harness builds an
    # argument string for; the two plugins the plugin gate's runs alone read are named in the
    # line after it.
    hint=$(grep -m1 -E '^    --plugins [A-Za-z0-9,]+ \\$' <<<"$out" | sed -E 's/^    --plugins //; s/ \\$//')
    unbuilt=""
    for plugin in ${hint//,/ }; do
        args_on "$ASSEMBLY" "$plugin" /pd rust >/dev/null 2>&1 || unbuilt="$unbuilt $plugin"
    done
    gate_only="UTRAnnotator,gnomADMt"; [[ "$ASSEMBLY_LC" == grch37 ]] && gate_only="UTRAnnotator"
    if [[ -n "$hint" && -z "$unbuilt" ]] && grep -qxF "run_concordance.sh builds no argument string for ${gate_only}: the runs of scripts/concordance/plugin_gate_runs.toml read their files (docs/plugins/reproducing-concordance.md)." <<<"$out"; then
        ok "setup_plugin_data.sh on $ASSEMBLY: the harness command names only plugins the harness runs"
    else
        bad "setup_plugin_data.sh on $ASSEMBLY: the harness command names '$hint' (no argument string for:$unbuilt) or the line on $gate_only is missing"
    fi
    # gnomADMt's VCF is read through tabix, so a tree without its index is not READY, while
    # UTRAnnotator's table is READY with no index beside it.
    if [[ "$ASSEMBLY_LC" == grch38 ]]; then
        rm -f "$tree/$(plugin_data_file gnomADMt file "$ASSEMBLY_LC").tbi"
        if ! mt_out=$(PATH="$stubs:$PATH" bash "$setup" --assembly "$ASSEMBLY" --data-dir "$tree" --source upstream --plugins gnomADMt,UTRAnnotator 2>&1) && grep -q 'gnomADMt: MISSING' <<<"$mt_out" && grep -q 'UTRAnnotator: READY' <<<"$mt_out" && ! grep -qF 'run_concordance.sh --assembly' <<<"$mt_out"; then
            ok "data-preparation summary on $ASSEMBLY: gnomADMt without its tabix index reads MISSING beside a READY UTRAnnotator, with no harness command"
        else
            bad "data-preparation summary on $ASSEMBLY: gnomADMt without its tabix index was not reported: $(grep -E 'gnomADMt|UTRAnnotator|run_concordance' <<<"$mt_out" | tail -3 | tr '\n' ' ')"
        fi
        : >"$tree/$(plugin_data_file gnomADMt file "$ASSEMBLY_LC").tbi"
    else
        mt_out=$(PATH="$stubs:$PATH" bash "$setup" --assembly "$ASSEMBLY" --data-dir "$tree" --source upstream --plugins gnomADMt 2>&1)
        status=$?
        if [[ "$status" -eq 1 ]] && grep -qxF "ERROR: no applicable plugins for $ASSEMBLY" <<<"$mt_out"; then
            ok "setup_plugin_data.sh on $ASSEMBLY: gnomADMt alone leaves no applicable plugin and exits 1"
        else
            bad "setup_plugin_data.sh on $ASSEMBLY: gnomADMt alone: status $status, printed $(tail -2 <<<"$mt_out" | tr '\n' ' ')"
        fi
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

# Under --source s3 a file of the synthetic conservation run that the copy does not bring is
# named in one line and the run goes on, to the plugins after LoFTEE and to the summary; any
# other file the copy does not bring ends the run there with status 1. The aws stub writes
# every object but MANIFEST.json and the keys of its pattern.
aws_stub() {
    printf '#!/bin/sh\ncase "$3" in\n*/MANIFEST.json | %s) exit 1 ;;\nesac\n: >"$4"\n' "$2" >"$1/aws"
    chmod +x "$1/aws"
}
stubs=$(mktemp -d "${TMPDIR:-/tmp}/plugin_data_layout_tools.XXXXXX")
tree=$(mktemp -d "${TMPDIR:-/tmp}/plugin_data_layout.XXXXXX")
aws_stub "$stubs" '*/loftee_synthetic_conservation/*'
want_absent=""
while IFS= read -r rel; do
    origin="written by tools/vep-loftee-synthetic-conservation"
    [[ "$rel" == *human_ancestor.fa.gz* ]] && origin="made from the ancestral FASTA"
    want_absent="$want_absent    [absent] $rel: $origin, see docs/plugins/reproducing-concordance.md"$'\n'
done < <(plugin_synthetic_run_files LoFTEE grch37)
if out=$(PATH="$stubs:$PATH" bash "$setup" --assembly GRCh37 --data-dir "$tree" --source s3 --bucket made-up-bucket --plugins LoFTEE,pLI 2>&1) && [[ "$(grep -F '[absent]' <<<"$out")"$'\n' == "$want_absent" ]] && grep -qxF "    [sync] $(plugin_data_file pLI file grch37)" <<<"$out" && grep -qxF "  LoFTEE: READY" <<<"$out" && grep -qxF "  LoFTEE synthetic conservation run: ABSENT" <<<"$out" && grep -qxF "  pLI: READY" <<<"$out"; then
    ok "setup_plugin_data.sh --source s3: the six files of the synthetic conservation run the copy lacks are named and the run goes on"
else
    bad "setup_plugin_data.sh --source s3 without the synthetic conservation run: $(grep -E 'absent|ERROR|LoFTEE|pLI' <<<"$out" | tail -4 | tr '\n' ' ')"
fi
rm -rf "$tree"
tree=$(mktemp -d "${TMPDIR:-/tmp}/plugin_data_layout.XXXXXX")
aws_stub "$stubs" "*/$(plugin_data_file LoFTEE gerp grch37)"
out=$(PATH="$stubs:$PATH" bash "$setup" --assembly GRCh37 --data-dir "$tree" --source s3 --bucket made-up-bucket --plugins LoFTEE,pLI 2>&1)
status=$?
if [[ "$status" -eq 1 ]] && grep -qx "    \[ERROR\] not in the bucket: .*/made-up-bucket/plugin_data/grch37/$(plugin_data_file LoFTEE gerp grch37)" <<<"$out" && ! grep -qE '\[absent\]|^\[pLI\]$|^=== ' <<<"$out"; then
    ok "setup_plugin_data.sh --source s3: a published file the copy lacks ends the run with status 1"
else
    bad "setup_plugin_data.sh --source s3 without LoFTEE's GERP table: status $status, printed $(tail -3 <<<"$out" | tr '\n' ' ')"
fi
rm -rf "$tree" "$stubs"

# Every bundle the script stages is human: its help says so in one line, and an assembly that
# is neither GRCh37 nor GRCh38 is refused before the data directory is made.
if out=$(bash "$setup" --help 2>&1) && grep -qxF 'The plugin data it stages is for the two human assemblies, GRCh37 and GRCh38.' <<<"$out"; then
    ok "setup_plugin_data.sh --help: the staged plugin data is for the two human assemblies"
else
    bad "setup_plugin_data.sh --help does not say the staged plugin data is for GRCh37 and GRCh38"
fi
tree=$(mktemp -d "${TMPDIR:-/tmp}/plugin_data_layout.XXXXXX")
out=$(bash "$setup" --assembly Xyz1 --data-dir "$tree/data" 2>"$tree/stderr")
status=$?
if [[ "$status" -eq 1 && -z "$out" && "$(cat "$tree/stderr")" == "ERROR: --assembly must be GRCh37 or GRCh38 (got 'Xyz1'): every plugin bundle is human" && ! -e "$tree/data" ]]; then
    ok "setup_plugin_data.sh on a made-up assembly: one ERROR line, exit 1, nothing staged"
else
    bad "setup_plugin_data.sh on a made-up assembly: status $status, printed '$out', stderr '$(cat "$tree/stderr")'"
fi
rm -rf "$tree"

# A bucket copy holds the DosageSensitivity table once, in the shared scope, and a LOEUF table
# in each assembly's scope: setup_plugin_data.sh's s3_scope_for, extracted by name.
scope_start=$(grep -n '^s3_scope_for()' "$setup" | cut -d: -f1)
scope_end=$(awk -v s="$scope_start" 'NR>s && /^}/{print NR; exit}' "$setup")
[[ -n "$scope_start" && -n "$scope_end" ]] || { echo "ERROR: [test_plugin_data_layout] s3_scope_for not found in $setup" >&2; exit 1; }
eval "$(sed -n "${scope_start},${scope_end}p" "$setup")"
for ASSEMBLY_LC in grch37 grch38; do
    if [[ "$(s3_scope_for DosageSensitivity)" == assembly_agnostic ]]; then ok "DosageSensitivity on $ASSEMBLY_LC: the shared scope"; else bad "DosageSensitivity on $ASSEMBLY_LC: scope '$(s3_scope_for DosageSensitivity)', expected assembly_agnostic"; fi
    if [[ "$(s3_scope_for LOEUF)" == "$ASSEMBLY_LC" ]]; then ok "LOEUF on $ASSEMBLY_LC: the assembly's own scope"; else bad "LOEUF on $ASSEMBLY_LC: scope '$(s3_scope_for LOEUF)', expected $ASSEMBLY_LC"; fi
    if [[ "$(s3_scope_for gnomADMt)" == assembly_agnostic ]]; then ok "gnomADMt on $ASSEMBLY_LC: the shared scope"; else bad "gnomADMt on $ASSEMBLY_LC: scope '$(s3_scope_for gnomADMt)', expected assembly_agnostic"; fi
    if [[ "$(s3_scope_for UTRAnnotator)" == "$ASSEMBLY_LC" ]]; then ok "UTRAnnotator on $ASSEMBLY_LC: the assembly's own scope"; else bad "UTRAnnotator on $ASSEMBLY_LC: scope '$(s3_scope_for UTRAnnotator)', expected $ASSEMBLY_LC"; fi
done

# Every plugin data path a whole run of plugin_gate_runs.toml names is <scope>/<a file the layout
# stages for a plugin setup_plugin_data.sh keeps under that scope>, as that file's header says. A
# path the layout lacks is a file setup_plugin_data.sh prints no source for and does not read
# before it reports READY.
all_line=$(grep -m1 '^ALL_PLUGINS=(' "$setup")
[[ -n "$all_line" ]] || { echo "ERROR: [test_plugin_data_layout] ALL_PLUGINS not found in $setup" >&2; exit 1; }
eval "$all_line"
unstaged=""
named=0
while IFS= read -r path; do
    named=$((named + 1))
    scope=${path%%/*}
    rel=${path#*/}
    found=0
    for ASSEMBLY_LC in grch37 grch38; do
        [[ "$scope" == assembly_agnostic || "$scope" == "$ASSEMBLY_LC" ]] || continue
        for plugin in "${ALL_PLUGINS[@]}"; do
            [[ "$(s3_scope_for "$plugin")" == "$scope" ]] || continue
            staged=$(plugin_staged_files "$plugin" "$ASSEMBLY_LC" 2>/dev/null) || continue
            grep -qxF "$rel" <<<"$staged" && found=1
        done
    done
    [[ "$found" -eq 1 ]] || unstaged="$unstaged $path"
done < <(grep -o '{plugin_data}/[^",]*' "$runs" | sed 's#^{plugin_data}/##' | sort -u)
if [[ "$named" -ge 30 && -z "$unstaged" ]]; then
    ok "every plugin data path of plugin_gate_runs.toml ($named) is a file the layout stages under its scope"
else
    bad "plugin_gate_runs.toml names $named plugin data paths; the layout stages no file for:$unstaged"
fi

# The harness names no flat plugin data file; every path of the plugin data directory comes
# from the layout table.
if grep -qE '/revel\.tsv\.gz|/gnomad_coverage\.tsv\.gz' "$harness"; then bad "run_concordance.sh names a flat plugin data file"; else ok "run_concordance.sh names no flat plugin data file"; fi

# The note setup_plugin_data.sh prints for the GRCh37 LoFTEE bundle names the key the harness
# passes to vep-rs there, without which vep-rs refuses a GRCh37 run.
if grep -q 'Pass it as gerp_tabix= with run_splice_predictions=0' "$setup"; then ok "setup_plugin_data.sh's GRCh37 LoFTEE note names run_splice_predictions=0"; else bad "setup_plugin_data.sh's GRCh37 LoFTEE note does not name run_splice_predictions=0"; fi

# The layout, the setup script and the harness word the GERP file LoFTEE reads on GRCh37 as
# the plugin's master branch reads it; none says in what form a download of it comes.
if grep -qE 'GRCh37 \(master branch\) ships GERP|GRCh37 as a tabix per-base TSV|tabix TSV vs bigWig' "$HERE/plugin_data_layout.sh" "$setup" "$harness"; then bad "a script says in what form the GRCh37 GERP file of LoFTEE is published: $(grep -nE 'GRCh37 \(master branch\) ships GERP|GRCh37 as a tabix per-base TSV|tabix TSV vs bigWig' "$HERE/plugin_data_layout.sh" "$setup" "$harness" | cut -d: -f1,2 | tr '\n' ' ')"; else ok "the layout, the setup script and the harness say how the plugin reads its GRCh37 GERP file"; fi

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

# Under Docker Ensembl VEP reads a user's own file through a mount of its directory. The
# harness mounts the directory of a variable only for a run of the plugin that reads it, and
# ends with one line naming the variable where that directory does not exist:
# plugin_docker_mounts and user_plugin_file_reader, extracted by name, with the list of the
# variables. Every file path is made up: the checks make directories and open no file.
for fn in plugin_docker_mounts user_plugin_file_reader; do
    fn_start=$(grep -n "^${fn}()" "$harness" | cut -d: -f1)
    fn_end=$(awk -v s="$fn_start" 'NR>s && /^}/{print NR; exit}' "$harness")
    [[ -n "$fn_start" && -n "$fn_end" ]] || { echo "ERROR: [test_plugin_data_layout] $fn not found in $harness" >&2; exit 1; }
    eval "$(sed -n "${fn_start},${fn_end}p" "$harness")"
done
vars_line=$(grep -m1 '^USER_PLUGIN_FILE_VARS=(' "$harness")
[[ -n "$vars_line" ]] || { echo "ERROR: [test_plugin_data_layout] USER_PLUGIN_FILE_VARS not found in $harness" >&2; exit 1; }
eval "$vars_line"
user=$(mktemp -d "${TMPDIR:-/tmp}/plugin_data_layout_user.XXXXXX")
user=$(cd "$user" && pwd)
mkdir "$user/plugin_data" "$user/cadd" "$user/spliceai" "$user/avi" "$user/dbnsfp"
# mounts_of <plugin>...: the mounts of a Docker run of those plugins over a plugin data
# directory and no plugins checkout: the data directory's line, then one for each user's file.
# shellcheck disable=SC2034  # plugin_docker_mounts reads them by name
mounts_of() {
    local PLUGIN_DATA_DIR="$user/plugin_data" VEP_PLUGINS_DIR="$user/no_plugins_checkout" name
    local -a PLUGIN_NAMES=("$@") PERL_PLUGIN_ARGS=()
    for name in "$@"; do PERL_PLUGIN_ARGS+=(--plugin "$name"); done
    plugin_docker_mounts
}
data_mount="-v"$'\n'"$user/plugin_data:/work/plugin_data:ro"
# shellcheck disable=SC2034  # plugin_docker_mounts and build_plugin_args read them by name
CADD_SNV_FILE="$user/cadd/snv_synthetic.tsv.gz" CADD_INDEL_FILE="$user/cadd/indel_synthetic.tsv.gz" SPLICEAI_SNV_FILE="$user/spliceai/snv.synthetic.vcf.gz" SPLICEAI_INDEL_FILE="$user/spliceai/indel.synthetic.vcf.gz" AVI_FILE="$user/avi/avi_synthetic.tsv.gz" DBNSFP_FILE="$user/dbnsfp/dbNSFP4.9a_grch38_synthetic.gz" DBSCSNV_FILE="$user/absent/dbscsnv_synthetic.txt.gz"
gives "a run of a layout plugin beside every user variable: no user directory is mounted, and one that does not exist is not read" "$data_mount" mounts_of REVEL
gives "a run of a layout plugin and dbNSFP beside the other plugins' variables: the directory of dbNSFP's file alone, under its variable's name" "$data_mount"$'\n'"-v"$'\n'"$user/dbnsfp:/work/user_plugin_data/DBNSFP_FILE:ro" mounts_of REVEL dbNSFP
gives "a run of SpliceAI beside the other plugins' variables: the directory of its SNV file and of its indel file, each under its variable's name" "$data_mount"$'\n'"-v"$'\n'"$user/spliceai:/work/user_plugin_data/SPLICEAI_SNV_FILE:ro"$'\n'"-v"$'\n'"$user/spliceai:/work/user_plugin_data/SPLICEAI_INDEL_FILE:ro" mounts_of SpliceAI
refused "a run of dbscSNV whose file's directory does not exist: one line names the variable, and no mount is printed" "ERROR: [run_concordance] DBSCSNV_FILE names a file whose directory does not exist: $user/absent" mounts_of dbscSNV
SPLICEAI_INDEL_FILE="$user/absent/indel.synthetic.vcf.gz" refused "a run of SpliceAI whose indel file's directory does not exist: one line names the variable, and no mount is printed" "ERROR: [run_concordance] SPLICEAI_INDEL_FILE names a file whose directory does not exist: $user/absent" mounts_of SpliceAI
# A directory closed to the user has no path the harness can hand Docker. The superuser enters
# a directory of any mode, and under it the same run mounts the directory.
mkdir "$user/closed" && chmod 000 "$user/closed"
if (cd "$user/closed") 2>/dev/null; then
    SPLICEAI_INDEL_FILE="$user/closed/indel.synthetic.vcf.gz" gives "a run of SpliceAI whose indel file's directory the user enters whatever its mode: mounted under its variable's name" "$data_mount"$'\n'"-v"$'\n'"$user/spliceai:/work/user_plugin_data/SPLICEAI_SNV_FILE:ro"$'\n'"-v"$'\n'"$user/closed:/work/user_plugin_data/SPLICEAI_INDEL_FILE:ro" mounts_of SpliceAI
else
    SPLICEAI_INDEL_FILE="$user/closed/indel.synthetic.vcf.gz" refused "a run of SpliceAI whose indel file's directory cannot be entered: one line names the variable, and no mount is printed" "ERROR: [run_concordance] SPLICEAI_INDEL_FILE names a file whose directory cannot be entered: $user/closed" mounts_of SpliceAI
fi
chmod 700 "$user/closed"
# The reader of a variable is the plugin whose Docker argument string names the variable's
# mount, and a plugin cannot run without a variable only where it is that variable's reader: a
# variable whose reader the harness misnames is a path Ensembl VEP is given and its container
# does not hold.
unread=""
for var in "${USER_PLUGIN_FILE_VARS[@]}"; do
    reader=$(user_plugin_file_reader "$var")
    arg=""
    [[ -z "$reader" ]] || arg=$(PERL_VEP_MODE=docker build_plugin_args "$reader" /pd perl 2>/dev/null) || arg=""
    [[ "$arg" == *"/work/user_plugin_data/$var/"* ]] || unread="$unread $var"
done
for plugin in CADD SpliceAI AVI dbNSFP dbscSNV; do
    for var in $(user_plugin_file_vars "$plugin"); do
        [[ "$(user_plugin_file_reader "$var")" == "$plugin" ]] || unread="$unread $plugin:$var"
    done
done
if [[ "${#USER_PLUGIN_FILE_VARS[@]}" -gt 0 && -z "$unread" ]]; then ok "every variable of a user's own file names the plugin whose Docker argument string reads its mount"; else bad "a variable of a user's own file names no plugin, or one whose Docker argument string does not read its mount:$unread"; fi
unset CADD_SNV_FILE CADD_INDEL_FILE SPLICEAI_SNV_FILE SPLICEAI_INDEL_FILE AVI_FILE DBNSFP_FILE DBSCSNV_FILE
rm -rf "$user"

# The Ensembl VEP output cache is keyed on PERL_VEP_SIGNATURE, evaluated here as the harness writes
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

# The page that lists the files of a genome-wide row gives every path as the layout's and sends
# the reader to setup_plugin_data.sh for every file; it names no plugin outside the layout.
page="$ROOT/docs/plugins/reproducing-concordance.md"
if grep -qF 'the path of the table below, the one `scripts/data/plugin_data_layout.sh`' "$page" && ! grep -qE 'outside that layout|preparation of the other two' "$page"; then
    ok "reproducing-concordance.md: every plugin data path of a genome-wide row is the layout's"
else
    bad "reproducing-concordance.md says a plugin's data file is outside the layout: $(grep -nE 'outside that layout|preparation of the other two' "$page" | cut -d: -f1 | tr '\n' ' ')"
fi

echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
