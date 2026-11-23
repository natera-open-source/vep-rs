#!/usr/bin/env bash
# The plugin data layout, shared by setup_plugin_data.sh (which stages the files) and
# run_concordance.sh (which reads them). One table: a name that appears here is the name
# the harness passes to both engines and the name the staging step writes, so the two
# cannot disagree on where a plugin's data lives. Paths are relative to the plugin data
# directory; the assembly is lower case (grch37, grch38).
#
# The runs of scripts/concordance/plugin_gate_runs.toml name every plugin data file as
# <scope>/<a path of this table>. gnomADMt and UTRAnnotator are read by those runs alone:
# run_concordance.sh builds no argument string for either.
#
# Sourced, never executed. Functions:
#   plugin_data_file <plugin> <role> <assembly_lc>   one data file's relative path
#   plugin_staged_files <plugin> <assembly_lc>       every file staged for the plugin,
#                                                    data files first, then their indexes,
#                                                    then the plugin's synthetic run files
#   plugin_synthetic_run_files <plugin> <assembly_lc>
#                                                    the files of a run on made-up data,
#                                                    each data file before its indexes
# A row that names one file per assembly (LOEUF, LoFTEE, Conservation, UTRAnnotator,
# gnomADMt) has none on an assembly its source does not publish: there both functions
# print one ERROR line on standard error, nothing on standard output, and return 1.
#
# loftee_synthetic_conservation/ holds the files of the plugin gate's GRCh37 run of LoFTEE:
# a made-up GERP table and a made-up conservation database, both written by
# tools/vep-loftee-synthetic-conservation, and the ancestral FASTA. No source publishes a
# made-up file, so plugin_synthetic_run_files lists the directory's files apart from the
# published ones, and prints nothing with status 0 for any other plugin or assembly.

# plugin_layout_require_human_assembly <plugin> <assembly_lc>: returns 0 on grch37 and
# grch38, the assemblies a per-assembly human bundle is published for; on any other
# spelling or assembly it prints the layout's refusal for <plugin> and returns 1.
plugin_layout_require_human_assembly() {
	local plugin="$1" asm="$2"
	case "${asm}" in
	grch37 | grch38) return 0 ;;
	esac
	echo "ERROR: [plugin_data_layout] ${plugin} has no data file on ${asm}" >&2
	return 1
}

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
		plugin_layout_require_human_assembly LOEUF "${asm}" || return 1
		if [[ "${asm}" == "grch37" ]]; then
			echo "loeuf/loeuf_dataset.tsv.gz"
		else
			echo "loeuf/loeuf_dataset_grch38.tsv.gz"
		fi
		;;
	DosageSensitivity:file) echo "dosagesensitivity/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz" ;;
	UTRAnnotator:file)
		# One uORF table per assembly, under the names UTRAnnotator.pm's documentation
		# gives the two files (UTRAnnotator.pm:42-43).
		plugin_layout_require_human_assembly UTRAnnotator "${asm}" || return 1
		if [[ "${asm}" == "grch37" ]]; then
			echo "utrannotator/uORF_5UTR_GRCh37_PUBLIC.txt"
		else
			echo "utrannotator/uORF_5UTR_GRCh38_PUBLIC.txt"
		fi
		;;
	gnomADMt:file)
		# gnomADMt.pm's documentation calls the mitochondrial sites VCF GRCh38 based
		# (gnomADMt.pm:50).
		if [[ "${asm}" != "grch38" ]]; then
			echo "ERROR: [plugin_data_layout] gnomADMt has no data file on ${asm}" >&2
			return 1
		fi
		echo "gnomadmt/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
		;;
	GWAS:file) echo "gwas/gwas_catalog_associations.tsv" ;;
	# Ensembl's per-chromosome ancestral FASTA files concatenated and bgzipped, with the
	# .fai and .gzi samtools faidx writes; the assembly is the prefix's.
	AncestralAllele:file) echo "ancestral_allele/homo_sapiens_ancestor.fa.gz" ;;
	# The Platform's VEP export of every association, with its header line commented
	# and indexed with tabix -s 1 -b 2 -e 2; study_types= selects a family inside it.
	OpenTargets:file) echo "opentargets/open_targets_vep_full.tsv.gz" ;;
	# Ensembl Compara publishes the GERP bigWig for GRCh38 alone, and LoFTEE's GRCh38
	# bundle ships a copy of it (another Compara release), so Conservation reads the
	# file staged for LoFTEE.
	Conservation:file)
		if [[ "${asm}" != "grch38" ]]; then
			echo "ERROR: [plugin_data_layout] Conservation has no data file on ${asm}" >&2
			return 1
		fi
		plugin_data_file LoFTEE gerp "${asm}"
		;;
	LoFTEE:ancestor)
		plugin_layout_require_human_assembly LoFTEE "${asm}" || return 1
		echo "loftee/${asm}/human_ancestor.fa"
		;;
	LoFTEE:ancestor_bgzf)
		plugin_layout_require_human_assembly LoFTEE "${asm}" || return 1
		echo "loftee/${asm}/human_ancestor.fa.gz"
		;;
	LoFTEE:gerp)
		# The plugin's two branches read GERP from two kinds of file: master (GRCh37) a
		# file it reads with tabix, under the default name of its gerp_file (LoF.pm lines
		# 132 and 138 to 141), and grch38 a bigWig.
		plugin_layout_require_human_assembly LoFTEE "${asm}" || return 1
		if [[ "${asm}" == "grch37" ]]; then
			echo "loftee/${asm}/GERP_scores.final.sorted.txt.gz"
		else
			echo "loftee/${asm}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
		fi
		;;
	LoFTEE:synthetic_gerp | LoFTEE:synthetic_conservation | LoFTEE:synthetic_ancestor_bgzf)
		# The plugin's master branch reads the three under gerp_file (LoF.pm lines 132 and 138
		# to 141), conservation_file (lines 290 to 297 and 567 to 570; gerp_dist.pl lines 58 to
		# 62 and 129 to 133) and human_ancestor_fa, by bare sequence name (lines 552 to 553).
		if [[ "${asm}" != "grch37" ]]; then
			echo "ERROR: [plugin_data_layout] LoFTEE has no synthetic conservation file on ${asm}" >&2
			return 1
		fi
		case "${role}" in
		synthetic_gerp) echo "loftee_synthetic_conservation/GERP_scores.synthetic.txt.gz" ;;
		synthetic_conservation) echo "loftee_synthetic_conservation/phylocsf_gerp.synthetic.sql" ;;
		synthetic_ancestor_bgzf) echo "loftee_synthetic_conservation/human_ancestor.fa.gz" ;;
		esac
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
	REVEL | gnomADc | AlphaMissense | LOEUF | OpenTargets | gnomADMt)
		f="$(plugin_data_file "${plugin}" file "${asm}")" || return 1
		echo "${f}"
		echo "${f}.tbi"
		;;
	GWAS | LoFtool | pLI | DosageSensitivity | Conservation | UTRAnnotator)
		plugin_data_file "${plugin}" file "${asm}" || return 1
		;;
	AncestralAllele)
		f="$(plugin_data_file AncestralAllele file "${asm}")"
		echo "${f}"
		echo "${f}.fai"
		echo "${f}.gzi"
		;;
	LoFTEE)
		f="$(plugin_data_file LoFTEE ancestor "${asm}")" || return 1
		echo "${f}"
		echo "${f}.fai"
		f="$(plugin_data_file LoFTEE ancestor_bgzf "${asm}")" || return 1
		echo "${f}"
		echo "${f}.fai"
		echo "${f}.gzi"
		f="$(plugin_data_file LoFTEE gerp "${asm}")" || return 1
		echo "${f}"
		[[ "${asm}" == "grch37" ]] && echo "${f}.tbi"
		plugin_synthetic_run_files LoFTEE "${asm}" || return 1
		;;
	*)
		echo "ERROR: [plugin_data_layout] unknown plugin ${plugin}" >&2
		return 1
		;;
	esac
	return 0
}

plugin_synthetic_run_files() {
	local plugin="$1" asm="$2" f
	[[ "${plugin}" == "LoFTEE" && "${asm}" == "grch37" ]] || return 0
	f="$(plugin_data_file LoFTEE synthetic_gerp "${asm}")" || return 1
	echo "${f}"
	echo "${f}.tbi"
	plugin_data_file LoFTEE synthetic_conservation "${asm}" || return 1
	f="$(plugin_data_file LoFTEE synthetic_ancestor_bgzf "${asm}")" || return 1
	echo "${f}"
	echo "${f}.fai"
	echo "${f}.gzi"
}
