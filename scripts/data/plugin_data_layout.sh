#!/usr/bin/env bash
# The plugin data layout, shared by setup_plugin_data.sh (which stages the files) and
# run_concordance.sh (which reads them). One table: a name that appears here is the name
# the harness passes to both engines and the name the staging step writes, so the two
# cannot disagree on where a plugin's data lives. Paths are relative to the plugin data
# directory; the assembly is lower case (grch37, grch38).
#
# Sourced, never executed. Functions:
#   plugin_data_file <plugin> <role> <assembly_lc>   one data file's relative path
#   plugin_staged_files <plugin> <assembly_lc>       every file staged for the plugin,
#                                                    data files first, then their indexes

plugin_data_file() {
	local plugin="$1" role="$2" asm="$3"
	case "${plugin}:${role}" in
	CADD:snv) echo "cadd/whole_genome_SNVs.tsv.gz" ;;
	CADD:indel)
		# CADD v1.7 names the indel file per assembly.
		if [[ "${asm}" == "grch37" ]]; then
			echo "cadd/gnomad.genomes-exomes.r4.0.indel.tsv.gz"
		else
			echo "cadd/gnomad.genomes.r4.0.indel.tsv.gz"
		fi
		;;
	REVEL:file) echo "revel/revel.tsv.gz" ;;
	SpliceAI:snv) echo "spliceai/spliceai_scores.raw.snv.vcf.gz" ;;
	SpliceAI:indel) echo "spliceai/spliceai_scores.raw.indel.vcf.gz" ;;
	gnomADc:file) echo "gnomadc/gnomad_coverage.tsv.gz" ;;
	AlphaMissense:file) echo "alphamissense/alphamissense.tsv.gz" ;;
	dbNSFP:file) echo "dbnsfp/dbNSFP.gz" ;;
	dbscSNV:file) echo "dbscsnv/dbscsnv.txt.gz" ;;
	LoFtool:file) echo "loftool/LoFtool_scores.txt" ;;
	pLI:file) echo "pli/pLI_values.txt" ;;
	GWAS:file) echo "gwas/gwas_catalog_associations.tsv" ;;
	LoFTEE:ancestor) echo "loftee/${asm}/human_ancestor.fa" ;;
	LoFTEE:ancestor_bgzf) echo "loftee/${asm}/human_ancestor.fa.gz" ;;
	LoFTEE:gerp)
		# The two upstream LoFTEE bundles differ: GRCh37 (master branch) ships GERP as a
		# tabix per-base TSV, GRCh38 (grch38 branch) as a bigWig.
		if [[ "${asm}" == "grch37" ]]; then
			echo "loftee/${asm}/GERP_scores.final.sorted.txt.gz"
		else
			echo "loftee/${asm}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
		fi
		;;
	*)
		echo "ERROR: [plugin_data_layout] no data file for ${plugin} role ${role}" >&2
		return 1
		;;
	esac
}

plugin_staged_files() {
	local plugin="$1" asm="$2" f
	case "${plugin}" in
	CADD)
		for role in snv indel; do
			f="$(plugin_data_file CADD "${role}" "${asm}")"
			echo "${f}"
			echo "${f}.tbi"
		done
		;;
	REVEL | gnomADc | AlphaMissense | dbscSNV)
		f="$(plugin_data_file "${plugin}" file "${asm}")"
		echo "${f}"
		echo "${f}.tbi"
		;;
	SpliceAI)
		# The indel scores are not freely distributed, so only the SNV file is staged;
		# the harness passes an indel file only when one is present.
		f="$(plugin_data_file SpliceAI snv "${asm}")"
		echo "${f}"
		echo "${f}.tbi"
		;;
	GWAS | LoFtool | pLI)
		plugin_data_file "${plugin}" file "${asm}"
		;;
	LoFTEE)
		f="$(plugin_data_file LoFTEE ancestor "${asm}")"
		echo "${f}"
		echo "${f}.fai"
		f="$(plugin_data_file LoFTEE ancestor_bgzf "${asm}")"
		echo "${f}"
		echo "${f}.fai"
		echo "${f}.gzi"
		f="$(plugin_data_file LoFTEE gerp "${asm}")"
		echo "${f}"
		[[ "${asm}" == "grch37" ]] && echo "${f}.tbi"
		;;
	dbNSFP)
		# Unsupported (licence); nothing is staged.
		;;
	*)
		echo "ERROR: [plugin_data_layout] unknown plugin ${plugin}" >&2
		return 1
		;;
	esac
	return 0
}
