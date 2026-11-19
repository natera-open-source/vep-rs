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
	REVEL:file) echo "revel/revel.tsv.gz" ;;
	gnomADc:file) echo "gnomadc/gnomad_coverage.tsv.gz" ;;
	AlphaMissense:file) echo "alphamissense/alphamissense.tsv.gz" ;;
	LoFtool:file) echo "loftool/LoFtool_scores.txt" ;;
	pLI:file) echo "pli/pLI_values.txt" ;;
	LOEUF:file)
		# gnomAD publishes constraint per release: the GRCh37 file is the v2.1.1
		# per-transcript table prepared as LOEUF.pm:54-59 says, the GRCh38 file Ensembl's
		# cut of the v4.1.1 metrics with gene coordinates appended (LOEUF.pm:61-67).
		if [[ "${asm}" == "grch37" ]]; then
			echo "loeuf/loeuf_dataset.tsv.gz"
		else
			echo "loeuf/loeuf_dataset_grch38.tsv.gz"
		fi
		;;
	DosageSensitivity:file) echo "dosagesensitivity/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz" ;;
	GWAS:file) echo "gwas/gwas_catalog_associations.tsv" ;;
	# Ensembl's per-chromosome ancestral FASTA files concatenated and bgzipped, with the
	# .fai and .gzi samtools faidx writes; the assembly is the prefix's.
	AncestralAllele:file) echo "ancestral_allele/homo_sapiens_ancestor.fa.gz" ;;
	# The Platform's VEP export of every association, with its header line commented
	# and indexed with tabix -s 1 -b 2 -e 2; study_types= selects a family inside it.
	OpenTargets:file) echo "opentargets/open_targets_vep_full.tsv.gz" ;;
	# Ensembl Compara publishes the GERP bigWig for GRCh38 alone, and LoFTEE's GRCh38
	# bundle ships a copy of it (an earlier Compara release), so Conservation reads the
	# file staged for LoFTEE.
	Conservation:file)
		if [[ "${asm}" != "grch38" ]]; then
			echo "ERROR: [plugin_data_layout] Conservation has no data file on ${asm}" >&2
			return 1
		fi
		plugin_data_file LoFTEE gerp "${asm}"
		;;
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
	REVEL | gnomADc | AlphaMissense | LOEUF | OpenTargets)
		f="$(plugin_data_file "${plugin}" file "${asm}")"
		echo "${f}"
		echo "${f}.tbi"
		;;
	GWAS | LoFtool | pLI | DosageSensitivity | Conservation)
		plugin_data_file "${plugin}" file "${asm}" || return 1
		;;
	AncestralAllele)
		f="$(plugin_data_file AncestralAllele file "${asm}")"
		echo "${f}"
		echo "${f}.fai"
		echo "${f}.gzi"
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
	*)
		echo "ERROR: [plugin_data_layout] unknown plugin ${plugin}" >&2
		return 1
		;;
	esac
	return 0
}
