#!/usr/bin/env bash
# The harness's genomes as data. harness_genomes.tsv beside this file holds one row a
# genome under a header line: species, assembly, cache_version (the number Ensembl VEP's
# --cache_version takes), and fasta and json_cache, two paths under the harness's data
# root.
#
# Sourced, never executed. One function:
#   harness_genome_field <table> <species> <assembly> <column>   one field of a genome's row
#
# The last row for a genome wins, so a line appended to a table of defaults replaces a
# default. A genome or a column the table lacks, a table that cannot be
# read and a row whose field count is not the header's each print one ERROR line on
# standard error and return 1. The columns are found by name in the header line.
#
# No associative array is used: bash 3.2 has none.

harness_genome_field() {
	if [[ $# -ne 4 ]]; then
		echo "ERROR: [harness_genomes] harness_genome_field takes <table> <species> <assembly> <column>, not $# arguments" >&2
		return 1
	fi
	if [[ ! -f "$1" || ! -r "$1" ]]; then
		printf 'ERROR: [harness_genomes] %s: the table cannot be read\n' "$1" >&2
		return 1
	fi
	local answer
	# awk prints the field and exits 0, or prints what the table lacks and exits 1. The
	# three names reach it through the environment, which keeps a backslash as it is, and
	# the table on standard input, since awk reads a name=value operand as an assignment.
	if answer=$(HARNESS_GENOME_SPECIES="$2" HARNESS_GENOME_ASSEMBLY="$3" HARNESS_GENOME_COLUMN="$4" awk -F'\t' '
		BEGIN {
			# The "" makes each name a string: awk compares two values that both look
			# like numbers as numbers, and 1.0 would match a row of 1.00.
			species = ENVIRON["HARNESS_GENOME_SPECIES"] ""
			assembly = ENVIRON["HARNESS_GENOME_ASSEMBLY"] ""
			column = ENVIRON["HARNESS_GENOME_COLUMN"] ""
		}
		NR == 1 {
			width = NF
			for (i = 1; i <= NF; i++) {
				if ($i == "species") species_at = i
				if ($i == "assembly") assembly_at = i
				if ($i == column) column_at = i
			}
			if (!species_at || !assembly_at) {
				print "the header line names no species and assembly columns"
				failed = 1
				exit 1
			}
			if (!column_at) {
				gsub(/\t/, ", ")
				printf "no column %s; the columns are %s\n", column, $0
				failed = 1
				exit 1
			}
			next
		}
		NF != width {
			printf "line %d has %d fields where the header has %d\n", NR, NF, width
			failed = 1
			exit 1
		}
		$species_at == species && $assembly_at == assembly {
			value = $column_at
			found = 1
		}
		END {
			if (failed) exit 1
			if (NR == 0) {
				print "the table is empty"
				exit 1
			}
			if (!found) {
				printf "no row for the genome %s %s\n", species, assembly
				exit 1
			}
			print value
		}
	' <"$1"); then
		printf '%s\n' "$answer"
	else
		printf 'ERROR: [harness_genomes] %s: %s\n' "$1" "$answer" >&2
		return 1
	fi
}
