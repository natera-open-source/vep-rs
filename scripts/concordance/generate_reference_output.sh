#!/usr/bin/env bash
#
# generate_reference_output.sh -- write the Ensembl VEP reference outputs that
# run_clone_measurement.sh scores every engine against, a flag family's reference
# in any output format, or a golden corpus's expected outputs.
#
# The measurement harness READS an already-populated data tree (its header documents
# the layout) and nothing else in this repository writes the Perl ground truth into
# it. This script is that writer. It runs Ensembl VEP in the `ensemblorg/ensembl-vep`
# Docker image with the exact invocation the harness embeds for its Perl engine, and
# writes the outputs where the harness and both comparators expect them:
#
#   SNP/indel (s01-s06), one file per suite:
#     <data>/ground_truth/perl/snp_indel/<suite>/output.txt
#     from <data>/vcf/<asm>/all_variants_canonical/<base>.canonical.vcf.gz
#   Structural variants (s07/s08; s09/s10 are the harness's pointer suites), one
#   file per input VCF, 16 per assembly:
#     <data>/ground_truth/perl/sv_per_vcf/<asm>/<vcf-basename>.txt
#     from <data>/ground_truth/perl/sv_per_vcf/canonical_inputs/<asm>/*.vcf[.gz]
#
# Beside every output it writes <name>.provenance.json (image digest, flags, input and
# output digests, record and row counts, wall time) and renames VEP's own
# `<name>.txt_warnings.txt` companion to <name>.warnings.log, so no `*.txt` glob can
# ever pair a warnings file as an engine output.
#
# OUTPUT FORMATS. Ensembl VEP builds one OutputFactory per run from its output_format
# (Runner.pm:871 get_OutputFactory, OutputFactory.pm:97-102 FORMAT_MAP), so a reference
# in several formats is one container run per format: --format default|tab|vcf|json|all.
# In the suite layout the default format keeps the names above and the other formats
# take their own extension on the same stem (output.tsv, output.vcf, output.jsonl;
# <vcf-basename>.tsv, .vcf, .jsonl), so a `*.txt` glob still pairs only the default
# output. --no_headers is passed with the default format alone: the tab output's column
# line and the VCF output's `##INFO=<ID=CSQ...>` line are what name the columns of those
# two formats, and the JSON format has no header to drop.
#
# FLAG FAMILIES. --flags "<words>" appends the words, split on whitespace, after the
# harness's flags and before the format flag, so a reference for one Ensembl VEP flag
# family is the harness's invocation plus that family. A word this script already sets
# (-i, -o, --fork, --fasta, --tab, --no_stats, ...) is refused rather than passed twice,
# in the single-dash and abbreviated spellings vep's option parser accepts as well; the
# one exception is `--format <input format>`, which names the input the words read (an
# Ensembl-format, region or HGVS file) and stands in place of the script's own `--format
# vcf`, as a golden corpus flag set over such an input needs. A
# family carrying --nearest, in any of those spellings, mounts the Ensembl cache
# writable: Ensembl VEP writes transcript_gene_tss.txt into the cache directory on its
# first --nearest run and keeps the file, built under that run's transcript filters, for
# every later run; a run on a cache that lacks the file scans every region file once, so
# the first --nearest run on a cache goes alone.
#
# STATS. --stats keeps Ensembl VEP's statistics report beside each output as
# <stem>.summary.log instead of passing --no_stats; a suite run drops it unless --stats
# is passed, a corpus run keeps it unless --no-stats is passed, so every corpus carries
# the report its expected outputs were written with. The run passes --stats_text, because
# a run that names neither --stats_text nor --stats_html writes the HTML report alone
# (Runner.pm:787-789); the text file is <output>_summary.txt (BaseRunner.pm:254-271) and
# is renamed out of the `*.txt` namespace for the reason the warnings file is.
#
# GOLDEN CORPORA. --corpus <dir> writes a corpus's expected outputs in the layout
# crates/vep-cli/tests/common/mod.rs reads: expected/{default.txt,tab.txt,vcf.vcf,
# json.jsonl}.gz from <dir>/variants.vcf, with the manifest's `species` (homo_sapiens when
# absent), `assembly`, `release` (the cache version, as the golden tests pass it under the
# Ensembl cache layout; --cache-version when the manifest names none), `flags` and
# `fasta`, --fork 1, headers kept, the statistics report of each format beside them as
# expected/<format>.summary.log, and <dir>/provenance.json carrying the image digest,
# the vep command and the sha256 of every expected file's uncompressed bytes, the
# statistics reports included. A key the
# script does not own (the divergence classification, the pruned-cache summary) survives
# a regeneration, and `input` and `fasta` are filled only when absent, so a record's own
# description of them stands. --format must be `all` for a corpus: the four files are one
# contract. A manifest flag word may carry `{corpus}`, the corpus directory, so a flag
# can name a file beside the manifest (`--custom file={corpus}/custom/sites.vcf.gz,...`)
# without an absolute path: the container mounts the corpus read-only at /work/corpus
# and reads the word with `{corpus}` rewritten to that mount, and the golden tests expand
# the same word to the corpus directory on the host; the recorded `vep_command` keeps
# the word as the manifest spells it. A corpus container also runs with the mount as its
# working directory, as the golden tests run vep-rs in the corpus directory, so a flag
# word that names a corpus file by a relative path (`file=data/scores.tsv.gz`) reads the
# same file in both engines and spells the same path: a plugin whose header text repeats
# the path it was given then writes the same header line on both sides.
#
# --flag-set <name> runs the manifest's `flag_sets` entry of that name instead of the
# manifest's own run: the entry's `flags` are the whole flag list of the run, its `input`
# (`variants.vcf` by default) the file annotated, its `formats` the formats written, and
# the outputs land under its `expected_dir` (`expected_<name>` by default) with the
# statistics reports beside them; provenance.json gains the record of that run under
# `flag_sets.<name>` (the image, the vep command, the flags, the input, the directory and
# the sha256 of every file written) and keeps every other key. An entry whose `stdin` is
# true or whose `perl_mode` is `database` is refused: this script mounts its input as a
# file and runs offline.
#
# A manifest whose flags carry --refseq or --merged annotates against the RefSeq or the
# merged cache: Ensembl VEP appends `_refseq` and `_merged` to the species directory name
# (CacheDir.pm:312 at release/116.2), so the slice checked for completeness and recorded
# as `perl_cache` is <species>_refseq/<version>_<assembly> or <species>_merged/..., and
# --species stays the species. A suite run passing the flags through --flags reads the
# same directory.
#
# A corpus that runs an Ensembl VEP plugin names the data files the plugin reads under the
# manifest's `plugin_files` (plugin name to a list of paths inside the corpus, which the
# flags name through `{corpus}`): each must exist, the corpus is mounted, and
# provenance.json records the sha256 of each under `plugin_files`, so a reference is tied
# to the fixture bytes it was generated over. A flag word may carry `{plugins}`, the plugin
# directory each engine reads: for Ensembl VEP the directory holding the plugin modules
# (the VEP_plugins checkout, or a plugin's own checkout), given as --dir-plugins <dir> and
# mounted read-only at /work/plugins, and for vep-rs the corpus's `plugin_data` directory,
# which the golden tests substitute; so `--dir_plugins {plugins}` names each engine's
# directory and `loftee_path:{plugins}` LoFTEE's own checkout. The recorded `vep_command`
# keeps `{plugins}` as the manifest spells it.
#
# COMPANION TOOLS. --tool filter_vep runs Ensembl VEP's filter_vep in the same image over a
# corpus's `tool_cases`, the manifest list of {name, input, args, expected} entries the
# `filter_tool` test of crates/vep-cli reads: each case runs
#   filter_vep -i <corpus>/<input> <args...> [<perl_args...>] [<--tool-args words>]
# in the container, its standard output is the expected file (so a --count or --list case
# is captured like any other), written to <corpus>/<expected> only when every case ran,
# and provenance.json gains a `tool_cases` map with each case's command, expected file,
# sha256, data-line count and wall time, plus the digests under its `sha256` map keyed by
# the file's path below expected/, as every corpus records its expected files. A case's
# `perl_args` are words filter_vep alone
# takes (a database connection for --ontology: filter_vep names it --host, --port, --user
# and --version, not vep's --db_version); --tool-args "<words>" appends words to every
# case's command and is recorded with it; --tool-case <name> runs one case (repeatable).
# filter_vep reads no cache, so --tool needs no --perl-cache-dir; a case's --filter words
# are the case's and --flags, --format and the suite options are refused.
# --tool variant_recoder runs Ensembl VEP's variant_recoder the same way over a corpus's
# `tool_cases` (the `recode_tool` test of crates/vep-cli reads them): each case runs
#   variant_recoder -i <corpus>/<input> <args...> [<perl_args...>] [<--tool-args words>]
# and its standard output, the JSON array, is the expected file. variant_recoder reads
# the Ensembl database by construction (VariantRecoder.pm sets `database`), so every case's
# `perl_args` pin the release (`--db_version 116`) and the container needs the network; it
# reads no cache either. A manifest's `tool` names the tool its cases run, and a case may
# name its own; both must be the --tool given.
#
# THE CONTAINER'S ENVIRONMENT. --perl-env NAME=VALUE (repeatable) sets a variable in the
# container for every vep run, recorded under `perl_env` in each provenance record. Its
# use is Perl's hash order: a plugin that joins a hash's keys (NearestExonJB's exon items)
# prints them in the order Perl's per-process random hash seed gives, so a reference of
# such a run is reproduced only under `--perl-env PERL_HASH_SEED=0 --perl-env
# PERL_PERTURB_KEYS=0`, which fixes that order; the golden tests still compare the items
# as a multiset. A word without `=`, with an empty name or with whitespace is refused, and
# a --tool run, which runs no vep, refuses the option.
#
# THE IMAGE DIGEST. Before anything runs, the pulled image's RepoDigests must carry
# --perl-image-digest (default the release_116.2 tag's manifest-list digest, or
# VEP_PERL_IMAGE_DIGEST), so a moved tag cannot produce a reference under this name.
# --perl-image and --perl-image-digest move together.
#
# FORK AND BUFFER ARE PART OF THE GROUND TRUTH'S IDENTITY, NOT A TUNING KNOB.
#   SNP/indel: --fork 16 --buffer_size 5000, the settings the reported SNP/indel
#              reference set was generated with.
#   SV:        --fork 4 --buffer_size 5000, the settings the reported per-VCF SV
#              reference set was generated with (the `invocation` field of the SV
#              reference's provenance record under docs/concordance-provenance/).
#              VEP's transcript selection above --max_sv_size is batch-dependent, so
#              the fork count changes which variants share an input buffer and
#              therefore the row count. A regeneration must reuse the ground truth's
#              OWN settings, never the timing harness's.
#   Corpus:    --fork 1 --buffer_size 5000, the settings the committed corpora carry.
#
# THE FLAG IS `--dir_cache`. `--full_cache_dir` is not a VEP flag: VEP ignores it,
# finds no annotation source, emits one warning per record and exits 0. The
# output-volume gate below (rows >= input records per suite or per SV set, the same
# floor the harness's `assert_perl_annotated` applies) makes that failure mode loud
# here as well.
#
# THE CACHE MUST BE THE FULL-GENOME ENSEMBL CACHE. A chromosome-subset cache does not
# fail a VEP run; it under-annotates every record whose transcripts lie outside the
# contigs it holds and exits 0. Completeness is asserted on the cache tree
# itself (all 25 primary contigs under homo_sapiens/<version>_<assembly>/; another
# species' slice must exist and hold at least one chromosome directory, since its
# primary set is not known here), before anything runs, the same check the harness's
# `assert_perl_cache_complete` makes.
#
# Obtaining the cache (about 24 GB per assembly, gzip tar), with <v> the Ensembl
# release, 116:
#   GRCh38: https://ftp.ensembl.org/pub/release-<v>/variation/indexed_vep_cache/homo_sapiens_vep_<v>_GRCh38.tar.gz
#   GRCh37: https://ftp.ensembl.org/pub/grch37/release-<v>/variation/indexed_vep_cache/homo_sapiens_vep_<v>_GRCh37.tar.gz
#   Extract both under ONE directory, so that it contains
#   homo_sapiens/116_GRCh38/ and homo_sapiens/116_GRCh37/; that directory is
#   --perl-cache-dir (default <data>/caches/perl/vep-cache). VEP resolves
#   <species>/<cache_version>_<assembly> beneath the --dir_cache mount from --species,
#   --cache_version and --assembly (CacheDir.pm:306-329 at release/116.2); another
#   species' cache tarball extracts beside the human ones. Readiness is an artifact test,
#   not a log line: count the contig directories and confirm no tar process is still
#   writing.
#
# The reference FASTA must be the primary assembly the cache was built against
# (Ensembl's Homo_sapiens.GRCh37.75.dna.primary_assembly.fa and
# Homo_sapiens.GRCh38.dna.primary_assembly.fa), bgzipped or plain, with a .fai
# beside it. --fasta-for <species>:<assembly>=<path>, repeatable, names the FASTA of one
# genome; homo_sapiens:GRCh37 and homo_sapiens:GRCh38 default to
# <data>/reference/<asm>/genome.fa (VEP_FASTA_GRCH37, VEP_FASTA_GRCH38). VEP receives
# --fasta on every call: the comparators score against outputs generated WITH --fasta, and
# an asymmetric run hides UTR-level divergences.
# A corpus passes --fasta only when its manifest names one, as the golden tests do.
#
# This script downloads no cache, FASTA or input (it pulls the image when the host lacks
# it, and a tool case in database mode reads the Ensembl database from its container),
# needs no cloud credentials, and refuses to overwrite an existing reference output unless
# --force is passed: a reference is a pinned artifact whose digests the provenance records
# carry.
#
# Usage:
#   scripts/concordance/generate_reference_output.sh --data-dir <data> \
#     [--assembly <assembly>|both] [--suites all|snp-indel|sv|"s01 s07"] \
#     [--format default|tab|vcf|json|all] [--flags "<vep flags>"] [--stats] \
#     [--volume-floor records|any] \
#     [--perl-image ensemblorg/ensembl-vep:release_116.2] [--perl-image-digest sha256:<64 hex>] \
#     [--perl-cache-dir <dir>] [--cache-version 116] \
#     [--fasta <fa>] [--fasta-for <species>:<assembly>=<fa>]... \
#     [--fork 16] [--sv-fork 4] [--buffer-size 5000] \
#     [--snp-indel-gt-dir <dir>] [--sv-gt-dir <dir>] [--sv-inputs-dir <dir>] \
#     [--force] [--dry-run]
#   scripts/concordance/generate_reference_output.sh --corpus tests/golden/116/GRCh38 \
#     --perl-cache-dir <dir> [--flag-set <name>] [--no-stats] [--volume-floor records|any] \
#     [--force] [--dry-run]
#   scripts/concordance/generate_reference_output.sh --corpus tests/golden/116/GRCh38-filter-tool \
#     --tool filter_vep [--tool-case <name>]... [--tool-args "<words>"] [--force] [--dry-run]
#   scripts/concordance/generate_reference_output.sh --corpus tests/golden/116/GRCh38-recoder \
#     --tool variant_recoder [--tool-case <name>]... [--tool-args "<words>"] [--force] [--dry-run]
#
#   --dry-run runs every preflight check (inputs, FASTA, cache completeness, output
#   collisions) and prints the docker command for each output instead of running it.
#   A missing input fails a dry run exactly as it fails a real one.
set -euo pipefail

# Defaults. The data-tree locations and the VEP_* environment overrides are the same
# ones run_clone_measurement.sh reads, so a tree generated here is a tree the harness
# measures without further configuration.
DATA_DIR=""
ASSEMBLY="both"
SPECIES="homo_sapiens"
SUITES_FILTER="all"
FORMAT=""
EXTRA_FLAGS=""
# true under --stats, false under --no-stats; a suite run defaults to false, a corpus run to true.
KEEP_STATS=""
VOLUME_FLOOR="records"
CORPUS_DIR=""
# The plugin modules directory of a corpus run (--dir-plugins), its mount in the container
# and the flag-word placeholder that names it.
DIR_PLUGINS=""
PLUGINS_MOUNT="/work/plugins"
PLUGINS_PLACEHOLDER='{plugins}'
# The corpus manifest's plugin data files (`plugin_files`), space-separated, relative to the corpus.
CORPUS_PLUGIN_FILES=""
# The container's environment (--perl-env NAME=VALUE words), space-separated.
PERL_ENV=""
# The `flag_sets` entry a corpus run generates (--flag-set), its input and its output
# directory relative to the corpus, and the formats it writes (space-separated; every
# format when empty). The manifest's own run reads variants.vcf and writes expected/.
FLAG_SET=""
CORPUS_INPUT="variants.vcf"
CORPUS_EXPECTED_DIR="expected"
CORPUS_FORMATS=""
PERL_IMAGE="ensemblorg/ensembl-vep:release_116.2"
# The tag's manifest-list digest on Docker Hub; the pulled image's RepoDigests must carry it
# before anything is generated, so a moved tag cannot produce a reference under this name.
PERL_IMAGE_DIGEST="${VEP_PERL_IMAGE_DIGEST:-sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da}"
PERL_CACHE_DIR="${VEP_PERL_CACHE_DIR:-}"
PERL_CACHE_VERSION="${VEP_PERL_CACHE_VERSION:-116}"
CACHE_VERSION_SEEN=""
FASTA_ANY=""
# The FASTA of each genome, "<species>:<assembly>	<path>" one per line (--fasta-for); the
# two human genomes default under <data>/reference/ once --data-dir is known.
FASTA_TABLE=""
FASTA_GRCH37="${VEP_FASTA_GRCH37:-}"
FASTA_GRCH38="${VEP_FASTA_GRCH38:-}"
FORK=""
SV_FORK=4
BUFFER_SIZE=5000
SNP_INDEL_GT_DIR="${VEP_SNP_INDEL_GT_DIR:-}"
SV_GT_DIR="${VEP_SV_GT_DIR:-}"
SV_INPUTS_DIR=""
FORCE="false"
DRY_RUN="false"
TOOL=""
TOOL_ARGS=""
TOOL_CASES=""
# Options that select or place suite outputs, recorded as seen so --corpus can refuse them.
SUITE_OPTS_SEEN=""

usage() {
	cat <<'EOF'
Generate the Ensembl VEP reference outputs run_clone_measurement.sh scores against, a
flag family's reference in any output format, or a golden corpus's expected outputs.

Required (one of):
  --data-dir <path>        Data tree (vcf/, caches/, reference/, ground_truth/)
  --corpus <dir>           A golden corpus (variants.vcf, manifest.json): writes its
                           expected/{default.txt,tab.txt,vcf.vcf,json.jsonl}.gz and
                           provenance.json; needs --perl-cache-dir or --data-dir

Optional:
  --assembly <A>           The assembly to generate for, or both (default both); every
                           input is human, GRCh37 or GRCh38
  --suites <filter>        all (default) | snp-indel | sv | space-separated "s01 s07"
                           (s09/s10 are aliases of s07/s08: the real-world SV VCFs are
                           part of the 16-file per-assembly set)
  --format <F>             Output format: default | tab | vcf | json | all, one
                           container run per format (default: default; a corpus is
                           always all). Ensembl VEP's own --format vcf names its input.
  --flags "<words>"        Ensembl VEP flags appended after the harness's, split on
                           whitespace; a flag this script sets itself is refused, except
                           --format <input format>, which replaces its --format vcf
  --stats                  Keep the statistics report as <stem>.summary.log (passes
                           --stats_text) instead of --no_stats; a corpus keeps it by default
  --no-stats               Pass --no_stats and keep no report (a corpus's opt-out)
  --volume-floor <mode>    records (default): annotated rows >= input records per
                           suite or SV set, and every written VCF record or JSON object
                           annotated; any: at least one annotated row, for flags that
                           drop records or rows (--chr, --coding_only, --no_intergenic,
                           --filter_common, --allow_non_variant, --individual,
                           --individual_zyg, --dont_skip) or write one row per record
                           (--pick, --most_severe, --summary) on an input with records
                           vep skips
  --perl-image <ref>       Docker image (default ensemblorg/ensembl-vep:release_116.2)
  --perl-image-digest <sha256:hex>  Digest the pulled image must carry (default the
                           pinned tag's manifest-list digest, or VEP_PERL_IMAGE_DIGEST)
  --perl-cache-dir <dir>   Directory holding <species>/<ver>_<asm>/
                           (default <data>/caches/perl/vep-cache)
  --cache-version <N>      Ensembl cache release (default 116)
  --fasta <fa>             Reference FASTA; valid only with a single --assembly
  --fasta-for <species>:<assembly>=<fa>
                           The FASTA of one genome, repeatable; homo_sapiens:GRCh37 and
                           homo_sapiens:GRCh38 default to <data>/reference/<asm>/genome.fa
                           (VEP_FASTA_GRCH37, VEP_FASTA_GRCH38)
  --fork <N>               SNP/indel --fork (default 16; a corpus defaults to 1)
  --sv-fork <N>            SV --fork (default 4; the shipped SV ground truth's setting)
  --buffer-size <N>        --buffer_size for both (default 5000)
  --snp-indel-gt-dir <dir> Where <suite>/output.txt is written
  --sv-gt-dir <dir>        Where <asm>/<vcf-basename>.txt is written
  --sv-inputs-dir <dir>    Parent of the per-assembly SV input directories
                           (default <sv-gt-dir>/canonical_inputs)
  --dir-plugins <dir>      With --corpus: the Ensembl VEP plugin modules directory (the
                           VEP_plugins checkout, or a plugin's own), mounted read-only at
                           /work/plugins, the value every {plugins} flag word takes;
                           required when the manifest names plugin_files or {plugins}
  --perl-env NAME=VALUE    Set a variable in the container for every vep run, recorded under
                           perl_env (repeatable); PERL_HASH_SEED=0 with PERL_PERTURB_KEYS=0
                           fixes the order a plugin prints a hash's keys in
  --flag-set <name>        With --corpus: run the manifest's flag_sets entry <name> (its
                           flags, input and formats) into its expected_dir and record the
                           run under provenance.json's flag_sets.<name>
  --tool <tool>            With --corpus: run Ensembl VEP's filter_vep or variant_recoder over
                           the manifest's tool_cases and write each case's expected file
  --tool-case <name>       Run one tool case (repeatable; default every case)
  --tool-args "<words>"    Words appended to every tool case's command, recorded with it
  --force                  Overwrite an existing reference output
  --dry-run                Preflight everything, print the docker commands, run nothing
EOF
	exit 1
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--data-dir)
		DATA_DIR="$2"
		shift 2
		;;
	--corpus)
		CORPUS_DIR="$2"
		shift 2
		;;
	--assembly)
		ASSEMBLY="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--suites)
		SUITES_FILTER="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--format)
		FORMAT="$2"
		shift 2
		;;
	--flags)
		EXTRA_FLAGS="$2"
		shift 2
		;;
	--stats)
		KEEP_STATS="true"
		shift
		;;
	--no-stats)
		KEEP_STATS="false"
		shift
		;;
	--volume-floor)
		VOLUME_FLOOR="$2"
		shift 2
		;;
	--perl-image)
		PERL_IMAGE="$2"
		shift 2
		;;
	--perl-image-digest)
		PERL_IMAGE_DIGEST="$2"
		shift 2
		;;
	--perl-cache-dir)
		PERL_CACHE_DIR="$2"
		shift 2
		;;
	--cache-version)
		PERL_CACHE_VERSION="$2"
		CACHE_VERSION_SEEN="$1"
		shift 2
		;;
	--fasta)
		FASTA_ANY="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--fasta-for)
		[[ "$2" == *=* && "${2%%=*}" == *:* ]] || {
			echo "ERROR: [generate_reference_output] --fasta-for takes <species>:<assembly>=<path> (got '$2')" >&2
			exit 1
		}
		FASTA_TABLE="${FASTA_TABLE}${2%%=*}	${2#*=}
"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--sv-fork)
		SV_FORK="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--buffer-size)
		BUFFER_SIZE="$2"
		shift 2
		;;
	--snp-indel-gt-dir)
		SNP_INDEL_GT_DIR="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--sv-gt-dir)
		SV_GT_DIR="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--sv-inputs-dir)
		SV_INPUTS_DIR="$2"
		SUITE_OPTS_SEEN="$SUITE_OPTS_SEEN $1"
		shift 2
		;;
	--dir-plugins)
		DIR_PLUGINS="$2"
		shift 2
		;;
	--perl-env)
		if [[ "$2" != *=* || -z "${2%%=*}" || "$2" == *[[:space:]]* ]]; then
			echo "ERROR: [generate_reference_output] --perl-env takes NAME=VALUE without whitespace, got '$2'" >&2
			exit 1
		fi
		PERL_ENV="${PERL_ENV:+$PERL_ENV }$2"
		shift 2
		;;
	--flag-set)
		FLAG_SET="$2"
		shift 2
		;;
	--tool)
		TOOL="$2"
		shift 2
		;;
	--tool-case)
		TOOL_CASES="$TOOL_CASES $2"
		shift 2
		;;
	--tool-args)
		TOOL_ARGS="$2"
		shift 2
		;;
	--force)
		FORCE="true"
		shift
		;;
	--dry-run)
		DRY_RUN="true"
		shift
		;;
	-h | --help) usage ;;
	*)
		echo "ERROR: [generate_reference_output] unknown argument: $1" >&2
		usage
		;;
	esac
done

die() {
	echo "ERROR: [generate_reference_output] $*" >&2
	exit 1
}
log() {
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
}

# Validate arguments
if [[ -n "$TOOL" ]]; then
	[[ "$TOOL" == "filter_vep" || "$TOOL" == "variant_recoder" ]] || die "invalid --tool: $TOOL (filter_vep, variant_recoder)"
	[[ -n "$CORPUS_DIR" ]] || die "--tool runs a corpus's tool_cases; pass --corpus"
	[[ -z "$SUITE_OPTS_SEEN" ]] || die "--tool runs $TOOL, which takes no suite option; drop$SUITE_OPTS_SEEN"
	[[ -z "$EXTRA_FLAGS" && -z "$FORMAT" ]] || die "--tool takes each case's words from manifest.json; drop --flags and --format"
	[[ -z "$KEEP_STATS" ]] || die "--tool runs $TOOL, which writes no statistics; drop --stats and --no-stats"
	[[ -z "$DIR_PLUGINS" ]] || die "--tool runs $TOOL, which loads no plugin; drop --dir-plugins"
	[[ -z "$PERL_ENV" ]] || die "--tool runs $TOOL, which reads no vep run's environment; drop --perl-env"
	[[ -z "$FLAG_SET" ]] || die "--tool runs $TOOL over the manifest's tool_cases; a flag set is a vep run, drop --flag-set"
else
	[[ -z "$TOOL_ARGS" && -z "$TOOL_CASES" ]] || die "--tool-args and --tool-case need --tool"
fi
[[ -z "$FLAG_SET" || -n "$CORPUS_DIR" ]] || die "--flag-set names a flag_sets entry of a corpus manifest; pass --corpus"
if [[ -n "$CORPUS_DIR" && -n "$TOOL" ]]; then
	[[ -d "$CORPUS_DIR" ]] || die "corpus not a directory: $CORPUS_DIR"
	CORPUS_DIR=$(cd "$CORPUS_DIR" && pwd)
	[[ -f "$CORPUS_DIR/manifest.json" ]] || die "corpus needs manifest.json: $CORPUS_DIR"
	command -v python3 >/dev/null 2>&1 || die "python3 is required to read manifest.json and write provenance.json"
	FORMAT="all"
	[[ -n "$FORK" ]] || FORK=1
elif [[ -n "$CORPUS_DIR" ]]; then
	[[ -d "$CORPUS_DIR" ]] || die "corpus not a directory: $CORPUS_DIR"
	# Docker bind mounts take absolute host paths; the documented call names the corpus
	# relative to the repository root.
	CORPUS_DIR=$(cd "$CORPUS_DIR" && pwd)
	[[ -f "$CORPUS_DIR/variants.vcf" && -f "$CORPUS_DIR/manifest.json" ]] || die "corpus needs variants.vcf and manifest.json: $CORPUS_DIR"
	[[ -z "$SUITE_OPTS_SEEN" ]] || die "--corpus takes its assembly, flags, FASTA and output location from the corpus; drop$SUITE_OPTS_SEEN"
	[[ -z "$EXTRA_FLAGS" ]] || die "--corpus takes its flags from manifest.json; drop --flags"
	if [[ -n "$DIR_PLUGINS" ]]; then
		[[ -d "$DIR_PLUGINS" ]] || die "--dir-plugins is not a directory: $DIR_PLUGINS"
		DIR_PLUGINS=$(cd "$DIR_PLUGINS" && pwd)
	fi
	[[ -z "$FORMAT" || "$FORMAT" == "all" ]] || die "--corpus writes all four formats (--format all); a corpus with one format regenerated is not a corpus"
	FORMAT="all"
	[[ -n "$FORK" ]] || FORK=1
	[[ -n "$KEEP_STATS" ]] || KEEP_STATS="true"
	[[ -n "$DATA_DIR" || -n "$PERL_CACHE_DIR" ]] || die "--corpus needs --perl-cache-dir (or --data-dir for its default cache location)"
	command -v python3 >/dev/null 2>&1 || die "python3 is required to read manifest.json and write provenance.json"
else
	[[ -n "$DATA_DIR" ]] || {
		echo "ERROR: [generate_reference_output] --data-dir or --corpus required" >&2
		usage
	}
	[[ -z "$DIR_PLUGINS" ]] || die "--dir-plugins mounts a corpus run's plugin modules; a suite run takes none"
	[[ -n "$FORMAT" ]] || FORMAT="default"
	[[ -n "$FORK" ]] || FORK=16
	[[ -n "$KEEP_STATS" ]] || KEEP_STATS="false"
fi
[[ -z "$DATA_DIR" || -d "$DATA_DIR" ]] || die "data-dir not a directory: $DATA_DIR"
# An assembly is a directory name under the species directory (CacheDir.pm:329 at
# release/116.2 matches it after the version prefix), so it takes the characters a name does.
ASSEMBLY_RE='^[A-Za-z0-9][A-Za-z0-9._-]*$'
[[ "$ASSEMBLY" == both || "$ASSEMBLY" =~ $ASSEMBLY_RE ]] || die "invalid --assembly: '$ASSEMBLY' (an assembly name, or both)"
case "$FORMAT" in
default | tab | vcf | json | all) ;;
*) die "invalid --format: $FORMAT (default|tab|vcf|json|all)" ;;
esac
case "$VOLUME_FLOOR" in
records | any) ;;
*) die "invalid --volume-floor: $VOLUME_FLOOR (records|any)" ;;
esac
for _n in "$FORK" "$SV_FORK" "$BUFFER_SIZE" "$PERL_CACHE_VERSION"; do
	[[ "$_n" =~ ^[0-9]+$ && "$_n" -gt 0 ]] || die "--fork, --sv-fork, --buffer-size and --cache-version take a positive integer (got '$_n')"
done
[[ "$PERL_IMAGE_DIGEST" =~ ^sha256:[0-9a-f]{64}$ ]] || die "--perl-image-digest must be sha256:<64 hex> (got '$PERL_IMAGE_DIGEST')"
if [[ -n "$FASTA_ANY" ]]; then
	[[ "$ASSEMBLY" != "both" ]] || die "--fasta needs a single --assembly; use --fasta-for <species>:<assembly>=<path> with --assembly both"
	FASTA_TABLE="${FASTA_TABLE}${SPECIES}:${ASSEMBLY}	${FASTA_ANY}
"
fi
if [[ "$DRY_RUN" != "true" ]]; then
	command -v docker >/dev/null 2>&1 || die "docker is required (pass --dry-run to preflight and print the commands without it)"
fi

[[ -n "$PERL_CACHE_DIR" ]] || PERL_CACHE_DIR="$DATA_DIR/caches/perl/vep-cache"
[[ -n "$FASTA_GRCH37" ]] || FASTA_GRCH37="$DATA_DIR/reference/grch37/genome.fa"
[[ -n "$FASTA_GRCH38" ]] || FASTA_GRCH38="$DATA_DIR/reference/grch38/genome.fa"
# The two human defaults go first and the last line for a genome wins the lookup, so a
# --fasta-for line replaces a default and a repeated one its predecessor.
FASTA_TABLE="homo_sapiens:GRCh37	${FASTA_GRCH37}
homo_sapiens:GRCh38	${FASTA_GRCH38}
${FASTA_TABLE}"
[[ -n "$SNP_INDEL_GT_DIR" ]] || SNP_INDEL_GT_DIR="$DATA_DIR/ground_truth/perl/snp_indel"
[[ -n "$SV_GT_DIR" ]] || SV_GT_DIR="$DATA_DIR/ground_truth/perl/sv_per_vcf"
[[ -n "$SV_INPUTS_DIR" ]] || SV_INPUTS_DIR="$SV_GT_DIR/canonical_inputs"

# fasta_for <species> <assembly>: the genome's FASTA from the table, or nothing.
fasta_for() { printf '%s' "$FASTA_TABLE" | awk -F'\t' -v k="$1:$2" '$1 == k {v = $2} END {if (v != "") print v}'; }
# `${var,,}` is bash 4+ and fails at runtime on macOS bash 3.2 while passing `bash -n`.
lower() { printf %s "$1" | tr '[:upper:]' '[:lower:]'; }

# The flags this script passes itself. A second copy inside --flags would win or lose
# by position in Ensembl VEP's option parsing, so it is refused instead. vep reads its
# options with Getopt::Long's defaults (vep:27, Config.pm:62 at release/116.2), which
# take a single dash or a plus sign, any unambiguous abbreviation and any letter case as
# the full option, so a word is refused when the lower-cased option name it carries is a
# prefix of one of these.
OWN_FLAGS="--input_file --output_file --offline --cache --dir_cache --cache_version --species --assembly --format --buffer_size --fork --fasta --no_headers --no_stats --stats_text --stats_html --stats_file --quiet --force_overwrite --tab --vcf --json --output_format"
# The INFO key Ensembl VEP writes its consequences under in the VCF format, the
# `--vcf_info_field` of the flag family or CSQ (OutputFactory/VCF.pm:312).
VCF_INFO_FIELD="CSQ"
# option_name <word>: the option named by a vep argument word, in its lower-case `--name`
# spelling with a value suffix removed, or nothing for a value word. `sf` is the one alias
# of an owned option (stats_file|sf) that is not a prefix of the option's own name.
option_name() {
	local name="${1%%=*}"
	case "$name" in
	--[A-Za-z]*) name=$(lower "$name") ;;
	-[A-Za-z]* | +[A-Za-z]*) name=$(lower "--${name#[-+]}") ;;
	*) return 0 ;;
	esac
	[[ "$name" != "--sf" ]] || name="--stats_file"
	printf '%s' "$name"
}
# has_option <--name>: whether a --flags word names the option <--name> in a spelling vep
# reads as it, that is, whether the word's option_name is a prefix of <--name>.
has_option() {
	local word name
	local -a words
	read -r -a words <<<"$EXTRA_FLAGS"
	for word in ${words[@]+"${words[@]}"}; do
		name=$(option_name "$word")
		[[ -n "$name" && "$1" == "$name"* ]] && return 0
	done
	return 1
}
# owned_flags_prefixed <--name>: how many owned flags the option name is a prefix of.
owned_flags_prefixed() {
	local own n=0
	for own in $OWN_FLAGS; do
		[[ "$own" == "$1"* ]] && n=$((n + 1))
	done
	printf '%s' "$n"
}
check_extra_flags() {
	local word name own take_field="" take_format=""
	set -f
	# shellcheck disable=SC2086
	set -- $EXTRA_FLAGS
	set +f
	# One space between words from here on: docker_cmd splits the rendered flags on spaces
	# alone, so a tab or newline inside --flags would reach vep as part of a word or end
	# the command early, and the provenance records would carry it verbatim.
	EXTRA_FLAGS="$*"
	for word in "$@"; do
		if [[ -n "$take_field" ]]; then
			VCF_INFO_FIELD="$word"
			take_field=""
			continue
		fi
		if [[ -n "$take_format" ]]; then
			take_format=""
			continue
		fi
		name=$(option_name "$word")
		[[ -n "$name" ]] || continue
		# An exact option or alias name outranks an abbreviation in Getopt::Long, and three
		# vep names are prefixes of owned flags: `verbose|v` (Config.pm:159), `terms|t`
		# (:176) and `dir` (:241); vep reads them as themselves, so they pass.
		case "$name" in
		--v | --t | --dir) continue ;;
		esac
		# `--format <input format>` names the input the words read (an Ensembl-format,
		# region or HGVS file) and stands in place of this script's own `--format vcf`
		# (vep_args); a spelling that is also a prefix of another owned flag (`--for`:
		# --fork, --force_overwrite) is ambiguous to vep and refused with it.
		if [[ "--format" == "$name"* && $(owned_flags_prefixed "$name") -eq 1 ]]; then
			[[ "$word" == *=* ]] || take_format=1
			continue
		fi
		for own in $OWN_FLAGS; do
			if [[ "$own" == "$name"* ]]; then
				die "--flags carries $word, which vep reads as $own, a flag this script sets; use its own option (--format, --fork, --stats, --no-stats, --fasta, --buffer-size, --cache-version) or drop it"
			fi
		done
		# `--vcf` and its prefixes were refused above, so a word that is still a prefix of
		# --vcf_info_field names that option alone.
		if [[ "--vcf_info_field" == "$name"* ]]; then
			if [[ "$word" == *=* ]]; then VCF_INFO_FIELD="${word#*=}"; else take_field=1; fi
		fi
	done
	[[ -z "$take_field" ]] || die "--flags ends with --vcf_info_field and no field name"
	[[ -z "$take_format" ]] || die "--flags ends with --format and no format name"
}

# A corpus takes its species, assembly, flags and FASTA from manifest.json; under
# --flag-set the named `flag_sets` entry supplies the flags, the input, the formats and the
# output directory in place of the manifest's own run. The flags are a list of words, each
# one argument to vep, so a word holding whitespace cannot be expressed. An existing
# provenance.json is parsed here too, so a damaged record stops the run before the
# container runs rather than after them.
CORPUS_FASTA_NAME=""
if [[ -n "$CORPUS_DIR" && -z "$TOOL" ]]; then
	manifest_fields=$(python3 - "$CORPUS_DIR/manifest.json" "$SPECIES" "$FLAG_SET" 2>&1 <<'EOF'
import json, re, sys
from pathlib import Path
manifest_path = Path(sys.argv[1])
default_species = sys.argv[2]
flag_set = sys.argv[3]
FORMATS = ("default", "tab", "vcf", "json")


def load(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except ValueError as e:
        sys.exit(f"{path}: {e}")


def word_list(flags, owner):
    if not isinstance(flags, list):
        sys.exit(f"{owner} flags must be a list of words, got {flags!r}")
    bad = [f for f in flags if not isinstance(f, str) or f != f.strip() or any(c.isspace() for c in f)]
    if bad:
        sys.exit(f"{owner} flags must be single words: {bad!r}")
    return flags


def corpus_path(value, owner):
    if not isinstance(value, str) or not value or value.startswith("/") or ".." in value.split("/") or any(c.isspace() for c in value):
        sys.exit(f"{owner} must be a path inside the corpus, got {value!r}")
    return value


m = load(manifest_path)
record_path = manifest_path.with_name("provenance.json")
record = load(record_path) if record_path.is_file() else {}
if not isinstance(record, dict):
    sys.exit(f"{record_path}: the record is not a JSON object")
if flag_set and "flag_sets" in record and not isinstance(record["flag_sets"], dict):
    sys.exit(f"{record_path}: flag_sets is not a JSON object")
flags = word_list(m.get("flags") or [], "manifest")
run_input = "variants.vcf"
expected_dir = "expected"
formats = list(FORMATS)
if flag_set:
    sets = m.get("flag_sets") or []
    if not isinstance(sets, list) or not all(isinstance(fs, dict) for fs in sets):
        sys.exit(f"manifest flag_sets must be a list of entries, got {sets!r}")
    found = [fs for fs in sets if fs.get("name") == flag_set]
    if not found:
        sys.exit(f"manifest flag_sets has no entry named {flag_set!r}")
    if len(found) > 1:
        sys.exit(f"manifest flag_sets has more than one entry named {flag_set!r}")
    entry = found[0]
    owner = f"flag set {flag_set}"
    flags = word_list(entry.get("flags") or [], owner)
    if entry.get("stdin"):
        sys.exit(f"{owner} reads its input on standard input, which this script does not run; it mounts the input as a file")
    if entry.get("perl_mode", "offline") != "offline":
        sys.exit(f"{owner} was generated in perl_mode {entry['perl_mode']!r}; this script runs offline")
    run_input = corpus_path(entry.get("input", "variants.vcf"), f"{owner} input")
    if not (manifest_path.parent / run_input).is_file():
        sys.exit(f"{owner} names an input the corpus lacks: {run_input}")
    expected_dir = corpus_path(entry.get("expected_dir", f"expected_{flag_set}"), f"{owner} expected_dir")
    formats = entry.get("formats", list(FORMATS))
    if not isinstance(formats, list) or not formats or any(f not in FORMATS for f in formats) or len(set(formats)) != len(formats):
        sys.exit(f"{owner} formats must list distinct names among {', '.join(FORMATS)}, got {formats!r}")
assembly = m.get("assembly") or ""
if not (isinstance(assembly, str) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", assembly)):
    sys.exit(f"manifest assembly must be an assembly name, got {assembly!r}")
species = m.get("species", default_species)
if not (isinstance(species, str) and re.fullmatch(r"[a-z][a-z0-9_]*", species)):
    sys.exit(f"manifest species must be an Ensembl species name (lower-case letters, digits and underscores), got {species!r}")
release = m.get("release", "")
if not (isinstance(release, str) and re.fullmatch(r"[0-9]*", release)):
    sys.exit(f"manifest release must be the cache version's digits, got {release!r}")
plugin_files = m.get("plugin_files") or {}
if not isinstance(plugin_files, dict) or not all(isinstance(v, list) for v in plugin_files.values()):
    sys.exit(f"manifest plugin_files must map each plugin to a list of paths inside the corpus, got {plugin_files!r}")
paths = []
for plugin, files in plugin_files.items():
    for f in files:
        if not isinstance(f, str) or not f or f.startswith("/") or ".." in f.split("/") or any(c.isspace() for c in f):
            sys.exit(f"manifest plugin_files for {plugin} must be paths inside the corpus, got {f!r}")
        if not (manifest_path.parent / f).is_file():
            sys.exit(f"manifest plugin_files for {plugin} names a file the corpus lacks: {f}")
        paths.append(f)
print(species)
print(assembly)
print(release)
print(m.get("fasta") or "")
print(" ".join(flags))
print(" ".join(paths))
print(run_input)
print(expected_dir)
print(" ".join(formats))
EOF
	) || die "corpus: $manifest_fields"
	SPECIES=$(printf '%s\n' "$manifest_fields" | sed -n 1p)
	ASSEMBLY=$(printf '%s\n' "$manifest_fields" | sed -n 2p)
	corpus_release=$(printf '%s\n' "$manifest_fields" | sed -n 3p)
	CORPUS_FASTA_NAME=$(printf '%s\n' "$manifest_fields" | sed -n 4p)
	EXTRA_FLAGS=$(printf '%s\n' "$manifest_fields" | sed -n 5p)
	CORPUS_PLUGIN_FILES=$(printf '%s\n' "$manifest_fields" | sed -n 6p)
	CORPUS_INPUT=$(printf '%s\n' "$manifest_fields" | sed -n 7p)
	CORPUS_EXPECTED_DIR=$(printf '%s\n' "$manifest_fields" | sed -n 8p)
	CORPUS_FORMATS=$(printf '%s\n' "$manifest_fields" | sed -n 9p)
	# A plugin's module directory is not part of the corpus: the data files are, the
	# code that reads them is the user's checkout.
	if [[ -n "$CORPUS_PLUGIN_FILES" || "$EXTRA_FLAGS" == *"$PLUGINS_PLACEHOLDER"* ]] && [[ -z "$DIR_PLUGINS" ]]; then
		die "the manifest names plugin_files or $PLUGINS_PLACEHOLDER, so the run loads an Ensembl VEP plugin; pass --dir-plugins <dir>, the directory holding its module (the VEP_plugins checkout, or the plugin's own)"
	fi
	# The manifest's release is the cache the corpus was cut from: the golden tests pass it as
	# --cache_version under the Ensembl layout, so the reference is generated under the same one.
	if [[ -n "$corpus_release" ]]; then
		[[ "$corpus_release" =~ ^[1-9][0-9]*$ ]] || die "manifest release must be a positive integer, the cache version (got '$corpus_release')"
		[[ -z "$CACHE_VERSION_SEEN" || "$PERL_CACHE_VERSION" == "$corpus_release" ]] || die "--cache-version $PERL_CACHE_VERSION differs from the manifest's release $corpus_release; drop $CACHE_VERSION_SEEN"
		PERL_CACHE_VERSION="$corpus_release"
	fi
fi
check_extra_flags
# The placeholder names the directory a corpus run mounts; a suite run has none.
CORPUS_PLACEHOLDER='{corpus}'
CORPUS_MOUNT="/work/corpus"
[[ -n "$CORPUS_DIR" || "$EXTRA_FLAGS" != *"$CORPUS_PLACEHOLDER"* ]] || die "--flags carries $CORPUS_PLACEHOLDER, which names the corpus directory of a --corpus run; a suite run has no corpus"
[[ -n "$CORPUS_DIR" || "$EXTRA_FLAGS" != *"$PLUGINS_PLACEHOLDER"* ]] || die "--flags carries $PLUGINS_PLACEHOLDER, which names the plugin directory of a --corpus run; a suite run has no plugin directory"

# Suite table: id:name:assembly:type:vcf_base:species -- the SAME rows as
# run_clone_measurement.sh, minus the pointer rows s09 and s10, which should_run maps onto
# s07 and s08, and minus the sample-scaling rows, which need no Ensembl VEP ground truth of
# their own: genotypes change no consequence, so the sites-only reference is the reference
# for both forms, and the genotype matrix's cost is measured separately.
SUITES="
s01:ClinVar_GRCh37:GRCh37:snp_indel:clinvar_grch37:homo_sapiens
s02:ClinVar_GRCh38:GRCh38:snp_indel:clinvar_grch38:homo_sapiens
s03:gnomAD_v2.1.1_chr21:GRCh37:snp_indel:gnomad_genomes_v2.1.1_chr21:homo_sapiens
s04:gnomAD_v4.1_chr21:GRCh38:snp_indel:gnomad_genomes_v4.1_chr21:homo_sapiens
s05:1KG_integrated_chr21:GRCh37:snp_indel:1kg_integrated_chr21_sitesonly:homo_sapiens
s06:1KG_highcov_chr21:GRCh38:snp_indel:1kg_highcov_chr21:homo_sapiens
s07:SV_GRCh37:GRCh37:sv:-:homo_sapiens
s08:SV_GRCh38:GRCh38:sv:-:homo_sapiens
"

# should_run <suite_id> <stype> <assembly>: honor --suites and --assembly. s09/s10 select
# the per-assembly SV set they point at.
should_run() {
	local sid="$1" stype="$2" asm="$3" f="$SUITES_FILTER"
	[[ "$ASSEMBLY" == "both" || "$ASSEMBLY" == "$asm" ]] || return 1
	f="${f// s09/ s07}"
	f="${f// s10/ s08}"
	f="${f/#s09/s07}"
	f="${f/#s10/s08}"
	case "$f" in
	all) return 0 ;;
	snp-indel) [[ "$stype" == "snp_indel" ]] ;;
	sv) [[ "$stype" == "sv" ]] ;;
	*) [[ " $f " == *" $sid "* ]] ;;
	esac
}

# The formats --format selects, in the order they run.
formats_selected() {
	if [[ -n "$CORPUS_FORMATS" ]]; then
		echo "$CORPUS_FORMATS"
		return
	fi
	case "$FORMAT" in
	all) echo "default tab vcf json" ;;
	*) echo "$FORMAT" ;;
	esac
}

# out_name_for <kind> <base> <format>: the output file name. A suite output keeps the
# harness's name for the default format and takes the format's extension otherwise; a
# corpus output takes the names the golden tests read.
out_name_for() {
	local kind="$1" base="$2" format="$3"
	if [[ "$kind" == "corpus" ]]; then
		case "$format" in
		default) echo "default.txt" ;;
		tab) echo "tab.txt" ;;
		vcf) echo "vcf.vcf" ;;
		json) echo "json.jsonl" ;;
		esac
		return
	fi
	case "$format" in
	default) echo "$base.txt" ;;
	tab) echo "$base.tsv" ;;
	vcf) echo "$base.vcf" ;;
	json) echo "$base.jsonl" ;;
	esac
}

# sidecar_stem <kind> <format> <out_name>: what the provenance, warnings and summary
# sidecars are named after: the output name less a default-format `.txt`, or the format
# name for a corpus output.
sidecar_stem() {
	if [[ "$1" == "corpus" ]]; then echo "$2"; else echo "${3%.txt}"; fi
}

# headers_for <kind> <format>: whether the run keeps Ensembl VEP's header lines. The
# harness reads the default suite output without them; every other output keeps them,
# because the tab column line and the VCF `##INFO` line name the columns.
headers_for() {
	if [[ "$1" != "corpus" && "$2" == "default" ]]; then echo "no"; else echo "yes"; fi
}

# Helpers
sha256_of() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	elif command -v shasum >/dev/null 2>&1; then
		shasum -a 256 "$1" | awk '{print $1}'
	else
		echo "unavailable"
	fi
}

json_escape() { printf '%s' "$1" | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g'; }

# count_records <vcf|vcf.gz>: non-header lines. `n=$(...) || n=0`, never `$(... || echo 0)`:
# grep -c prints 0 AND exits 1 on no match, so the latter yields the two-line string "0\n0"
# and every later numeric test on it errors into its else branch.
count_records() {
	local n
	case "$1" in
	*.gz) n=$(gzip -dc "$1" | grep -vc '^#') || n=0 ;;
	*) n=$(grep -vc '^#' "$1") || n=0 ;;
	esac
	echo "$n"
}

count_data_lines() {
	local n
	n=$(grep -vc '^#' "$1") || n=0
	echo "$n"
}

# annotated_units <format> <file>: what the output-volume gate counts. The default and tab
# formats write one row per consequence. The VCF format writes every buffered record as a
# line and appends the CSQ field only to a record with consequences (OutputFactory/VCF.pm:
# 299-361), and the JSON format writes one object per record with a `*_consequences` array
# only on an annotated one (OutputFactory/JSON.pm:213-260, :384), so for those two the unit
# is an annotated record.
annotated_units() {
	local n
	case "$1" in
	default | tab) n=$(count_data_lines "$2") ;;
	vcf) n=$(FIELD="$VCF_INFO_FIELD" awk -F'\t' '!/^#/ && index(";" $8, ";" ENVIRON["FIELD"] "=") {c++} END {print c + 0}' "$2") ;;
	json) n=$(grep -c '_consequences":' "$2") || n=0 ;;
	esac
	echo "$n"
}

# cache_species_dir <species>: the directory Ensembl VEP reads under --dir_cache for the
# species and the run's flags: the species, with `_refseq` and then `_merged` appended for
# each of --refseq and --merged the flags carry (CacheDir.pm:312 at release/116.2).
cache_species_dir() {
	local dir="$1"
	has_option --refseq && dir="${dir}_refseq"
	has_option --merged && dir="${dir}_merged"
	printf '%s' "$dir"
}

# assert_cache_complete <species> <assembly>: the full-genome check. Primary contigs only:
# Ensembl's full GRCh37 cache carries 699 directories, the other 674 being unplaced scaffolds
# and haplotype contigs no benchmark input reaches, so requiring all 699 would fail a
# legitimately pruned cache. Requiring the 25 primaries fails every chromosome-subset
# cache, which is the defect class. The primary set is known for homo_sapiens alone; another
# species' slice must hold at least one chromosome directory.
assert_cache_complete() {
	local species="$1" assembly="$2"
	local species_dir slice
	species_dir=$(cache_species_dir "$species")
	slice="$PERL_CACHE_DIR/$species_dir/${PERL_CACHE_VERSION}_${assembly}"
	local primaries="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
	local c missing="" n_chrom=0
	if [[ ! -d "$slice" ]]; then
		echo "ERROR: [generate_reference_output] Perl cache slice not found: $slice" >&2
		echo "  Expected $species_dir/${PERL_CACHE_VERSION}_${assembly} under --perl-cache-dir ($PERL_CACHE_DIR)." >&2
		echo "  See the header of this script for the Ensembl FTP location of the full-genome cache." >&2
		return 1
	fi
	if [[ "$species" != "homo_sapiens" ]]; then
		for c in "$slice"/*/; do
			[[ -d "$c" ]] && n_chrom=$((n_chrom + 1))
		done
		if [[ "$n_chrom" -eq 0 ]]; then
			echo "ERROR: [generate_reference_output] Perl cache slice holds no chromosome directory: $slice" >&2
			return 1
		fi
		log "perl cache check ($species $assembly): $n_chrom chromosome directories under $slice; completeness is asserted for homo_sapiens alone"
		return 0
	fi
	for c in $primaries; do
		[[ -d "$slice/$c" ]] || missing="$missing $c"
	done
	if [[ -n "$missing" ]]; then
		echo "ERROR: [generate_reference_output] Perl cache for $assembly is missing primary contigs:$missing" >&2
		echo "  Slice: $slice" >&2
		echo "  A chromosome-subset cache does not fail a VEP run; it under-annotates every record" >&2
		echo "  outside the contigs it holds and exits 0. A reference generated against it is wrong" >&2
		echo "  in a way no comparator can see. Extract the full-genome Ensembl cache first." >&2
		return 1
	fi
	log "perl cache check ($assembly): all 25 primary contigs present under $slice"
	return 0
}

# assert_fasta <species> <assembly>: the genome has a FASTA in the table, on disk, with its
# .fai beside it.
assert_fasta() {
	local species="$1" assembly="$2" fasta
	fasta=$(fasta_for "$species" "$assembly")
	[[ -n "$fasta" ]] || {
		echo "ERROR: [generate_reference_output] no FASTA named for $species:$assembly" >&2
		echo "  Pass --fasta-for $species:$assembly=<path> (or --fasta with a single --assembly)." >&2
		return 1
	}
	[[ -f "$fasta" ]] || {
		echo "ERROR: [generate_reference_output] FASTA for $assembly not found: $fasta" >&2
		echo "  Pass --fasta-for $species:$assembly=<path> (or --fasta with a single --assembly)." >&2
		return 1
	}
	[[ -f "${fasta}.fai" ]] || {
		echo "ERROR: [generate_reference_output] FASTA index not found: ${fasta}.fai" >&2
		echo "  VEP reads the index; build it with \`samtools faidx\` before generating." >&2
		return 1
	}
	return 0
}

# vep_args <species> <assembly> <fork> <format> <cache> <fasta|""> <headers yes|no>
# Everything after `vep -i <in> -o <out>`, one shell word each, space-joined: the
# harness's invocation, --fasta when one is given, --no_headers when the headers go, the
# stats flag, --quiet --force_overwrite, the --flags words and the output-format flag
# (none for the default format). docker_cmd and the provenance records both render this
# one function, so the command run and the command recorded cannot drift.
vep_args() {
	local species="$1" assembly="$2" fork="$3" format="$4" cache="$5" fasta="$6" headers="$7"
	local a="--offline --cache --dir_cache $cache --cache_version $PERL_CACHE_VERSION --species $species --assembly $assembly"
	# The --flags words name the input format themselves for a non-VCF input.
	if ! has_option --format; then a="$a --format vcf"; fi
	a="$a --buffer_size $BUFFER_SIZE --fork $fork"
	if [[ -n "$fasta" ]]; then a="$a --fasta $fasta"; fi
	if [[ "$headers" == "no" ]]; then a="$a --no_headers"; fi
	if [[ "$KEEP_STATS" == "true" ]]; then a="$a --stats_text"; else a="$a --no_stats"; fi
	a="$a --quiet --force_overwrite"
	if [[ -n "$EXTRA_FLAGS" ]]; then a="$a $EXTRA_FLAGS"; fi
	case "$format" in
	tab | vcf | json) a="$a --$format" ;;
	esac
	printf '%s' "$a"
}

# docker_cmd <species> <assembly> <fork> <format> <in.vcf[.gz]> <out_dir> <out_name> <fasta|""> <headers>
# Prints the exact docker invocation as one shell-quoted line. This is the harness's
# Perl invocation, flag for flag, so the reference and the timed runs cannot diverge.
docker_cmd() {
	local species="$1" assembly="$2" fork="$3" format="$4" in_vcf="$5" out_dir="$6" out_name="$7" fasta="$8" headers="$9"
	local in_mount fasta_in=""
	local -a argv
	case "$in_vcf" in
	*.gz) in_mount="/work/input.vcf.gz" ;;
	*) in_mount="/work/input.vcf" ;;
	esac
	argv=(docker run --rm
		--user "$(id -u):$(id -g)"
		-v "${in_vcf}:${in_mount}:ro"
		-v "${out_dir}:/work/perl_out")
	# Ensembl VEP's first --nearest run writes transcript_gene_tss.txt into the cache
	# directory and dies on a read-only one (AnnotationSource/Cache/Transcript.pm:288-289).
	if has_option --nearest; then
		argv+=(-v "${PERL_CACHE_DIR}:/work/cache")
	else
		argv+=(-v "${PERL_CACHE_DIR}:/work/cache:ro")
	fi
	if [[ -n "$fasta" ]]; then
		argv+=(-v "${fasta}:/work/ref.fa:ro" -v "${fasta}.fai:/work/ref.fa.fai:ro")
		fasta_in="/work/ref.fa"
	fi
	# A corpus run reads its flags' files from the corpus, by `{corpus}` or by a
	# relative path, so the corpus is mounted and is the working directory.
	if [[ -n "$CORPUS_DIR" ]]; then
		argv+=(-v "${CORPUS_DIR}:${CORPUS_MOUNT}:ro" -w "$CORPUS_MOUNT")
	fi
	if [[ -n "$DIR_PLUGINS" ]]; then
		argv+=(-v "${DIR_PLUGINS}:${PLUGINS_MOUNT}:ro")
	fi
	local env_word
	# The words are split on whitespace alone: a glob character in a value must
	# not expand against the working directory.
	set -f
	for env_word in $PERL_ENV; do
		argv+=(-e "$env_word")
	done
	set +f
	argv+=("$PERL_IMAGE" vep -i "$in_mount" -o "/work/perl_out/${out_name}")
	local -a words
	IFS=' ' read -r -a words <<<"$(vep_args "$species" "$assembly" "$fork" "$format" /work/cache "$fasta_in" "$headers")"
	words=("${words[@]//"$CORPUS_PLACEHOLDER"/$CORPUS_MOUNT}")
	words=("${words[@]//"$PLUGINS_PLACEHOLDER"/$PLUGINS_MOUNT}")
	argv+=("${words[@]}")
	printf '%q ' "${argv[@]}"
	echo
}

# perl_env_json: the --perl-env words as a JSON array of strings.
perl_env_json() {
	local word out="" sep=""
	set -f
	for word in $PERL_ENV; do
		out="$out$sep\"$(json_escape "$word")\""
		sep=", "
	done
	set +f
	printf '[%s]' "$out"
}

# write_provenance <json_path> <species> <assembly> <fork> <format> <headers> <in> <in_records> <out> <out_units> <out_lines> <wall> <image_digest>
write_provenance() {
	local json="$1" species="$2" assembly="$3" fork="$4" format="$5" headers="$6" in_vcf="$7" in_records="$8" out="$9" out_units="${10}" out_lines="${11}" wall="${12}" digest="${13}"
	local digest_json
	if [[ -n "$digest" ]]; then digest_json="\"$digest\""; else digest_json="null"; fi
	cat >"$json" <<EOF
{
  "generated_at_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "generator": "scripts/concordance/generate_reference_output.sh",
  "perl_image": "$PERL_IMAGE",
  "perl_image_digest": $digest_json,
  "species": "$species",
  "assembly": "$assembly",
  "cache_version": $PERL_CACHE_VERSION,
  "cache_dir": "$PERL_CACHE_DIR",
  "fasta": "$(fasta_for "$species" "$assembly")",
  "fork": $fork,
  "buffer_size": $BUFFER_SIZE,
  "output_format": "$format",
  "vep_flags": "$(json_escape "$(vep_args "$species" "$assembly" "$fork" "$format" '<cache>' '<fa>' "$headers")")",
  "perl_env": $(perl_env_json),
  "input": "$in_vcf",
  "input_sha256": "$(sha256_of "$in_vcf")",
  "input_records": $in_records,
  "output": "$out",
  "output_sha256": "$(sha256_of "$out")",
  "output_annotated_units": $out_units,
  "output_data_lines": $out_lines,
  "wall_seconds": $wall
}
EOF
}

# Plan: one line per output, "kind|species|assembly|fork|format|out_dir|input|out_name".
# Built in full BEFORE anything runs, so every missing input, FASTA, cache slice and output
# collision is reported together and nothing is generated against a half-valid tree. Lines
# sharing their first six fields form a GROUP (one SNP/indel suite, or one assembly's 16-file
# SV set, in one format), which is the unit the output-volume gate below is applied to. A
# newline-delimited string rather than an array of records keeps this bash-3.2 safe.
JOBS=""
PREFLIGHT_FAILED=0
CHECKED_GENOMES=""

# preflight_assembly <species> <assembly>: the FASTA and the cache slice, once per genome.
preflight_assembly() {
	local species="$1" asm="$2"
	[[ " $CHECKED_GENOMES " == *" $species:$asm "* ]] && return 0
	CHECKED_GENOMES="$CHECKED_GENOMES $species:$asm"
	if [[ -z "$CORPUS_DIR" ]]; then
		assert_fasta "$species" "$asm" || PREFLIGHT_FAILED=1
	fi
	assert_cache_complete "$species" "$asm" || PREFLIGHT_FAILED=1
}

# check_output_collision <path>
check_output_collision() {
	if [[ -e "$1" && "$FORCE" != "true" ]]; then
		echo "ERROR: [generate_reference_output] reference output already exists: $1" >&2
		echo "  A reference is a pinned artifact. Pass --force to overwrite it." >&2
		return 1
	fi
	return 0
}

# The pulled image must carry the pinned digest before anything runs; IMAGE_DIGEST is the
# `name@sha256:...` the provenance records.
IMAGE_DIGEST=""
check_image_digest() {
	local image_digests
	image_digests=$(docker image inspect --format '{{join .RepoDigests ","}}' "$PERL_IMAGE" 2>/dev/null || true)
	if [[ -z "$image_digests" ]]; then
		log "pulling $PERL_IMAGE"
		docker pull "$PERL_IMAGE" >/dev/null
		image_digests=$(docker image inspect --format '{{join .RepoDigests ","}}' "$PERL_IMAGE" 2>/dev/null || true)
	fi
	[[ "$image_digests" == *"@$PERL_IMAGE_DIGEST"* ]] ||
		die "image digest mismatch for $PERL_IMAGE: ${image_digests:-none} (expected $PERL_IMAGE_DIGEST); pass --perl-image-digest for another image, nothing was generated"
	IMAGE_DIGEST=$(printf '%s\n' "${image_digests//,/$'\n'}" | grep -F "@$PERL_IMAGE_DIGEST" | head -1)
	log "image digest $IMAGE_DIGEST"
}

# A tool corpus: every `tool_cases` entry of the manifest (or the --tool-case names) runs
# the tool in the container with its input mounted read-only under /work/corpus, the
# container's standard output captured as the case's expected file. The cases are planned
# in full first, run into one scratch directory, and moved into place together once every
# case exited 0 with output, so a failing case leaves the committed files untouched.
if [[ -n "$TOOL" ]]; then
	tool_plan=$(python3 - "$CORPUS_DIR" "$TOOL_CASES" "$TOOL" 2>&1 <<'EOF'
import json, shlex, sys
from pathlib import Path

corpus = Path(sys.argv[1])
selected = sys.argv[2].split()
tool_given = sys.argv[3]
path = corpus / "manifest.json"
try:
    manifest = json.loads(path.read_text(encoding="utf-8"))
except ValueError as e:
    sys.exit(f"{path}: {e}")
record_path = corpus / "provenance.json"
if record_path.is_file():
    try:
        record = json.loads(record_path.read_text(encoding="utf-8"))
    except ValueError as e:
        sys.exit(f"{record_path}: {e}")
    if not isinstance(record, dict):
        sys.exit(f"{record_path}: the record is not a JSON object")
cases = manifest.get("tool_cases")
if not isinstance(cases, list) or not cases:
    sys.exit("manifest tool_cases must be a non-empty list of cases")
manifest_tool = manifest.get("tool", "filter_vep")
if manifest_tool != tool_given:
    sys.exit(f"manifest tool is {manifest_tool!r}, not the --tool given ({tool_given})")
names = []
for case in cases:
    if not isinstance(case, dict):
        sys.exit(f"tool case is not an object: {case!r}")
    name = case.get("name")
    if not isinstance(name, str) or not name or any(c.isspace() or c in "/|" for c in name):
        sys.exit(f"tool case name must be one word without / or |: {name!r}")
    if name in names:
        sys.exit(f"tool case named twice: {name}")
    names.append(name)
    tool = case.get("tool", manifest_tool)
    if tool != tool_given:
        sys.exit(f"tool case {name}: tool must be {tool_given}, got {tool!r}")
    source = case.get("input")
    if not isinstance(source, str) or not source or "/" in source or source.startswith("."):
        sys.exit(f"tool case {name}: input must name a file of the corpus, got {source!r}")
    if not (corpus / source).is_file():
        sys.exit(f"tool case {name}: input not found: {corpus / source}")
    expected = case.get("expected")
    if not isinstance(expected, str) or not expected.startswith("expected/") or ".." in expected.split("/") or expected.endswith("/"):
        sys.exit(f"tool case {name}: expected must be a file path under expected/, got {expected!r}")
    for key, value in (("input", source), ("expected", expected)):
        if "\n" in value or "\t" in value:
            sys.exit(f"tool case {name}: {key} holds a tab or newline, which the plan's record format cannot carry")
    for key in ("args", "perl_args"):
        words = case.get(key, [])
        if not isinstance(words, list) or not all(isinstance(w, str) and w for w in words):
            sys.exit(f"tool case {name}: {key} must be a list of non-empty strings, got {words!r}")
        if any("\n" in w or "\t" in w for w in words):
            sys.exit(f"tool case {name}: {key} holds a tab or newline, which the command cannot carry")
    if selected and name not in selected:
        continue
    words = list(case.get("args", [])) + list(case.get("perl_args", []))
    print("\t".join([name, source, expected, " ".join(shlex.quote(w) for w in words)]))
missing = [n for n in selected if n not in names]
if missing:
    sys.exit(f"--tool-case names no case of the manifest: {' '.join(missing)}")
EOF
	) || die "tool corpus: $tool_plan"
	[[ -n "$tool_plan" ]] || die "tool corpus: no case selected"
	set -f
	# shellcheck disable=SC2086
	set -- $TOOL_ARGS
	set +f
	TOOL_ARGS="$*"
	tool_extra_quoted=""
	for word in "$@"; do tool_extra_quoted="$tool_extra_quoted $(printf '%q' "$word")"; done
	tool_failed=0
	while IFS=$'\t' read -r name source expected words; do
		[[ -n "$name" ]] || continue
		check_output_collision "$CORPUS_DIR/$expected" || tool_failed=1
	done <<<"$tool_plan"
	[[ "$tool_failed" -eq 0 ]] || die "preflight failed; nothing was generated"
	n_tool=$(printf '%s\n' "$tool_plan" | grep -c .) || n_tool=0
	log "generate_reference_output image=$PERL_IMAGE tool=$TOOL corpus=$CORPUS_DIR planned_cases=$n_tool dry_run=$DRY_RUN"
	[[ -z "$TOOL_ARGS" ]] || log "tool args: $TOOL_ARGS"
	# tool_docker_cmd <input> <quoted words>: the container invocation as one shell-quoted line,
	# its standard output the case's expected file.
	tool_docker_cmd() {
		printf '%q ' docker run --rm --user "$(id -u):$(id -g)" -v "${CORPUS_DIR}:${CORPUS_MOUNT}:ro" "$PERL_IMAGE" "$TOOL" -i "${CORPUS_MOUNT}/$1"
		printf '%s%s\n' "$2" "$tool_extra_quoted"
	}
	if [[ "$DRY_RUN" == "true" ]]; then
		while IFS=$'\t' read -r name source expected words; do
			[[ -n "$name" ]] || continue
			echo "DRY-RUN [$TOOL $name] -> $CORPUS_DIR/$expected"
			tool_docker_cmd "$source" "$words"
		done <<<"$tool_plan"
		log "dry run complete: $n_tool command(s) printed, nothing written"
		exit 0
	fi
	check_image_digest
	tool_scratch=$(mktemp -d "$CORPUS_DIR/.generate.XXXXXX")
	trap 'rm -rf "$tool_scratch"' EXIT
	tool_facts=""
	while IFS=$'\t' read -r name source expected words; do
		[[ -n "$name" ]] || continue
		log "  $TOOL $name: $source $words$tool_extra_quoted -> $expected"
		start=$(date +%s)
		set +e
		eval "$(tool_docker_cmd "$source" "$words")" >"$tool_scratch/$name.out" 2>"$tool_scratch/$name.err"
		rc=$?
		set -e
		wall=$(($(date +%s) - start))
		if [[ "$rc" -ne 0 ]]; then
			echo "ERROR: [generate_reference_output] $TOOL $name exited $rc after ${wall}s; every case's output was discarded" >&2
			head -5 "$tool_scratch/$name.err" | sed 's/^/  stderr: /' >&2
			exit 1
		fi
		[[ -s "$tool_scratch/$name.out" ]] || die "$TOOL $name exited 0 and wrote nothing; every case's output was discarded"
		lines=$(count_data_lines "$tool_scratch/$name.out")
		tool_facts="${tool_facts}${name}	${source}	${expected}	${words}	${lines}	${wall}
"
		log "  $lines data line(s) in ${wall}s"
	done <<<"$tool_plan"
	# The record and the expected files are written from the scratch directory in one pass.
	python3 - "$CORPUS_DIR" "$tool_scratch" "$IMAGE_DIGEST" "$PERL_IMAGE" "$TOOL" "$TOOL_ARGS" "$tool_facts" <<'EOF'
import hashlib, json, sys
from pathlib import Path

corpus, scratch, image, tag, tool, extra, facts = sys.argv[1:8]
corpus, scratch = Path(corpus), Path(scratch)
path = corpus / "provenance.json"
record = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
cases = record.get("tool_cases") if isinstance(record.get("tool_cases"), dict) else {}
# The digests also go under `sha256`, keyed by the file's path below expected/, the map
# every corpus carries for its expected files and the golden tests read.
digests = record.get("sha256") if isinstance(record.get("sha256"), dict) else {}
for line in facts.splitlines():
    if not line:
        continue
    name, source, expected, words, lines, wall = line.split("\t")
    out = scratch / f"{name}.out"
    command = f"{tool} -i {source} {words}".rstrip()
    if extra:
        command = f"{command} {extra}"
    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    cases[name] = {
        "input": source,
        "command": command,
        "expected": expected,
        "sha256": digest,
        "data_lines": int(lines),
        "wall_seconds": int(wall),
    }
    digests[expected[len("expected/"):]] = digest
    target = corpus / expected
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(out.read_bytes())
record["sha256"] = digests
record["tool"] = tool
record["tool_image"] = image
record["tool_image_tag"] = tag
record["tool_cases"] = cases
path.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
EOF
	while IFS=$'\t' read -r name source expected words lines wall; do
		[[ -n "$name" ]] || continue
		log "  wrote $CORPUS_DIR/$expected"
	done <<<"$tool_facts"
	log "  wrote $CORPUS_DIR/provenance.json"
	log "generate_reference_output complete: $n_tool/$n_tool case(s) written"
	exit 0
fi

n_selected=0
CORPUS_STAGE=""
CORPUS_FASTA=""
if [[ -n "$CORPUS_DIR" ]]; then
	# The corpus's four outputs are staged in one scratch directory and moved into
	# expected/ together once every format has passed its gate, so a failing format
	# leaves the committed expected files untouched.
	if [[ "$DRY_RUN" == "true" ]]; then
		CORPUS_STAGE="$CORPUS_DIR/.generate.dry-run"
	else
		CORPUS_STAGE=$(mktemp -d "$CORPUS_DIR/.generate.XXXXXX")
		trap 'rm -rf "$CORPUS_STAGE"' EXIT
	fi
	n_selected=1
	preflight_assembly "$SPECIES" "$ASSEMBLY"
	if [[ -n "$CORPUS_FASTA_NAME" ]]; then
		fasta_plain="${CORPUS_FASTA_NAME%.gz}"
		if [[ ! -f "$CORPUS_DIR/$CORPUS_FASTA_NAME" || ! -f "$CORPUS_DIR/$fasta_plain.fai" ]]; then
			echo "ERROR: [generate_reference_output] corpus FASTA or its index not found: $CORPUS_DIR/$CORPUS_FASTA_NAME, $CORPUS_DIR/$fasta_plain.fai" >&2
			PREFLIGHT_FAILED=1
		elif [[ "$CORPUS_FASTA_NAME" == *.gz ]]; then
			# The golden tests decompress the FASTA once per run and read it with the
			# plain file's index; the container does the same.
			CORPUS_FASTA="$CORPUS_STAGE/$fasta_plain"
		else
			CORPUS_FASTA="$CORPUS_DIR/$CORPUS_FASTA_NAME"
		fi
	fi
	for fmt in $(formats_selected); do
		out_name=$(out_name_for corpus "" "$fmt")
		check_output_collision "$CORPUS_DIR/$CORPUS_EXPECTED_DIR/$out_name.gz" || PREFLIGHT_FAILED=1
		check_output_collision "$CORPUS_DIR/$CORPUS_EXPECTED_DIR/$out_name" || PREFLIGHT_FAILED=1
		JOBS="${JOBS}corpus|$SPECIES|$ASSEMBLY|$FORK|$fmt|$CORPUS_STAGE|$CORPUS_DIR/$CORPUS_INPUT|$out_name
"
	done
fi
while IFS= read -r entry; do
	[[ -n "$entry" && -z "$CORPUS_DIR" ]] || continue
	IFS=':' read -r suite_id suite_name asm stype vcf_base species <<<"$entry"
	should_run "$suite_id" "$stype" "$asm" || continue
	n_selected=$((n_selected + 1))
	preflight_assembly "$species" "$asm"
	asm_lc=$(lower "$asm")
	if [[ "$stype" == "snp_indel" ]]; then
		in_vcf="$DATA_DIR/vcf/$asm_lc/all_variants_canonical/${vcf_base}.canonical.vcf.gz"
		if [[ ! -f "$in_vcf" ]]; then
			echo "ERROR: [generate_reference_output] $suite_id ($suite_name): input not found: $in_vcf" >&2
			echo "  Canonical inputs carry Ensembl-style contigs; derive them with prepare_benchmark_vcfs.py." >&2
			PREFLIGHT_FAILED=1
			continue
		fi
		for fmt in $(formats_selected); do
			out_name=$(out_name_for snp_indel output "$fmt")
			check_output_collision "$SNP_INDEL_GT_DIR/$suite_id/$out_name" || PREFLIGHT_FAILED=1
			JOBS="${JOBS}snp_indel|$species|$asm|$FORK|$fmt|$SNP_INDEL_GT_DIR/$suite_id|$in_vcf|$out_name
"
		done
	else
		in_dir="$SV_INPUTS_DIR/$asm_lc"
		if [[ ! -d "$in_dir" ]]; then
			echo "ERROR: [generate_reference_output] $suite_id ($suite_name): SV input directory not found: $in_dir" >&2
			PREFLIGHT_FAILED=1
			continue
		fi
		n_sv=0
		for fmt in $(formats_selected); do
			for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
				[[ -f "$vcf" ]] || continue
				b=$(basename "$vcf")
				b=${b%.vcf.gz}
				b=${b%.vcf}
				[[ "$b" == "manifest" ]] && continue
				n_sv=$((n_sv + 1))
				out_name=$(out_name_for sv "$b" "$fmt")
				check_output_collision "$SV_GT_DIR/$asm_lc/$out_name" || PREFLIGHT_FAILED=1
				JOBS="${JOBS}sv|$species|$asm|$SV_FORK|$fmt|$SV_GT_DIR/$asm_lc|$vcf|$out_name
"
			done
		done
		if [[ "$n_sv" -eq 0 ]]; then
			echo "ERROR: [generate_reference_output] $suite_id ($suite_name): no .vcf or .vcf.gz under $in_dir" >&2
			echo "  An empty input directory would produce an empty reference that scores as a clean cell." >&2
			PREFLIGHT_FAILED=1
		fi
	fi
done <<<"$SUITES"

[[ "$n_selected" -gt 0 ]] || die "no suite matched --suites '$SUITES_FILTER' with --assembly $ASSEMBLY"
[[ "$PREFLIGHT_FAILED" -eq 0 ]] || die "preflight failed; nothing was generated"
[[ -n "$JOBS" ]] || die "preflight selected $n_selected suite(s) but planned no output"
bad_record=$(printf '%s' "$JOBS" | awk -F'|' 'NF != 8 && !seen { print; seen = 1 }')
[[ -z "$bad_record" ]] || die "a planned path holds '|' or a newline, which the plan's record format cannot carry: $bad_record"

# fasta_for_job <kind> <species> <assembly>: the host FASTA a job mounts, empty for a corpus
# without one.
fasta_for_job() {
	if [[ "$1" == "corpus" ]]; then echo "$CORPUS_FASTA"; else fasta_for "$2" "$3"; fi
}

n_jobs=$(printf '%s' "$JOBS" | grep -c '|') || n_jobs=0
log "generate_reference_output image=$PERL_IMAGE cache_version=$PERL_CACHE_VERSION fork=$FORK sv_fork=$SV_FORK buffer_size=$BUFFER_SIZE format=$FORMAT stats=$KEEP_STATS volume_floor=$VOLUME_FLOOR"
[[ -z "$EXTRA_FLAGS" ]] || log "flags: $EXTRA_FLAGS"
if [[ -n "$CORPUS_DIR" ]]; then
	log "corpus=$CORPUS_DIR species=$SPECIES assembly=$ASSEMBLY${FLAG_SET:+ flag_set=$FLAG_SET} input=$CORPUS_INPUT expected_dir=$CORPUS_EXPECTED_DIR planned_outputs=$n_jobs dry_run=$DRY_RUN"
else
	log "data_dir=$DATA_DIR snp_indel_gt=$SNP_INDEL_GT_DIR sv_gt=$SV_GT_DIR planned_outputs=$n_jobs dry_run=$DRY_RUN"
fi

# Dry run: print every command and stop.
if [[ "$DRY_RUN" == "true" ]]; then
	while IFS='|' read -r kind species asm fork format out_dir in_vcf out_name; do
		[[ -n "$kind" ]] || continue
		if [[ "$kind" == "corpus" ]]; then
			echo "DRY-RUN [$kind $asm $format] -> $CORPUS_DIR/$CORPUS_EXPECTED_DIR/$out_name.gz"
		else
			echo "DRY-RUN [$kind $asm $format] -> $out_dir/$out_name"
		fi
		docker_cmd "$species" "$asm" "$fork" "$format" "$in_vcf" "$out_dir" "$out_name" "$(fasta_for_job "$kind" "$species" "$asm")" "$(headers_for "$kind" "$format")"
	done <<<"$JOBS"
	log "dry run complete: $n_jobs command(s) printed, nothing written"
	exit 0
fi

# Run, one group at a time. Every output of a group is produced in one scratch directory
# beside its destination, the group is gated on annotation volume, and only then is
# anything moved into place. An interrupted or non-annotating run therefore never leaves
# a partial file at a path the harness would trust.
#
# THE OUTPUT-VOLUME GATE is the harness's own floor, applied at the harness's own grain:
# summed rows >= summed input records, per SNP/indel suite and per assembly's SV set. A
# VEP that cannot reach its cache emits one warning per record, zero annotation rows and
# exit 0; a working VEP emits at least one row per input record (one per overlapping
# transcript, or a single intergenic row) and several per record on genome-wide and SV
# inputs alike. The floor is NOT applied per SV file: some synthetic
# files consist of allele classes VEP declines to annotate (`<*>`, `<NON_REF>`, `.`), so a
# legitimately sparse file would fail a per-file floor while the aggregate still bounds
# the failure mode that matters. A file with zero rows is reported, never silently kept.
# The VCF and JSON formats write one line per kept record whatever its consequences, so
# their floor is that every written record is annotated: a record VEP skipped before its
# input buffer is absent from both, while a record it could not annotate is present and
# bare. --volume-floor any lowers every format to one annotated unit per group, for a
# flag family that drops records or rows.
check_image_digest

if [[ -n "$CORPUS_FASTA" && "$CORPUS_FASTA_NAME" == *.gz ]]; then
	gzip -dc "$CORPUS_DIR/$CORPUS_FASTA_NAME" >"$CORPUS_FASTA"
	cp "$CORPUS_DIR/${CORPUS_FASTA_NAME%.gz}.fai" "$CORPUS_FASTA.fai"
fi

done_jobs=0
GROUP_KEYS=$(printf '%s' "$JOBS" | cut -d'|' -f1-6 | uniq)
while IFS='|' read -r kind species asm fork format out_dir; do
	[[ -n "$kind" ]] || continue
	group_prefix="$kind|$species|$asm|$fork|$format|$out_dir|"
	fasta=$(fasta_for_job "$kind" "$species" "$asm")
	headers=$(headers_for "$kind" "$format")
	mkdir -p "$out_dir"
	scratch=$(mktemp -d "$out_dir/.generate.XXXXXX")
	group_records=0
	group_rows=0
	group_written=0
	group_files=0
	# Per-file facts for the provenance sidecars, recorded before the gate so a failure
	# discards them with the scratch directory. "in_vcf|records|units|lines|wall|out_name".
	group_facts=""
	log "group [$kind $asm $format --fork $fork] -> $out_dir/"
	while IFS='|' read -r _k _s _a _f _fmt _o in_vcf out_name; do
		[[ -n "$_k" ]] || continue
		in_records=$(count_records "$in_vcf")
		if [[ "$in_records" -eq 0 ]]; then
			rm -rf "$scratch"
			die "$out_dir/$out_name: input $in_vcf has 0 records; an empty reference would score as a clean cell"
		fi
		log "  $in_vcf ($in_records records) -> $out_name"
		start=$(date +%s)
		set +e
		eval "$(docker_cmd "$species" "$asm" "$fork" "$format" "$in_vcf" "$scratch" "$out_name" "$fasta" "$headers")"
		rc=$?
		set -e
		wall=$(($(date +%s) - start))
		if [[ "$rc" -ne 0 ]]; then
			rm -rf "$scratch"
			die "$out_dir/$out_name: vep exited $rc after ${wall}s; the group's partial output was discarded"
		fi
		[[ -f "$scratch/$out_name" ]] || {
			rm -rf "$scratch"
			die "$out_dir/$out_name: vep exited 0 but wrote no output; the group's partial output was discarded"
		}
		if [[ "$KEEP_STATS" == "true" && ! -f "$scratch/${out_name}_summary.txt" ]]; then
			rm -rf "$scratch"
			die "$out_dir/$out_name: the statistics report was asked for (--stats, or a corpus run without --no-stats) and vep wrote no ${out_name}_summary.txt; the group's partial output was discarded"
		fi
		out_lines=$(annotated_units "$format" "$scratch/$out_name")
		written=$(count_data_lines "$scratch/$out_name")
		[[ "$out_lines" -gt 0 ]] || echo "WARN: [generate_reference_output] $out_dir/$out_name has 0 annotated $format units for $in_records records" >&2
		group_records=$((group_records + in_records))
		group_rows=$((group_rows + out_lines))
		group_written=$((group_written + written))
		group_files=$((group_files + 1))
		group_facts="${group_facts}${in_vcf}|${in_records}|${out_lines}|${written}|${wall}|${out_name}
"
		log "    $out_lines annotated units in ${wall}s"
	# The prefix reaches awk through the environment: `-v` would read a backslash in a path
	# as an escape sequence and the group would match none of its own records.
	done <<<"$(printf '%s' "$JOBS" | PREFIX="$group_prefix" awk 'index($0, ENVIRON["PREFIX"]) == 1')"

	gate_failed=""
	if [[ "$VOLUME_FLOOR" == "any" ]]; then
		[[ "$group_rows" -ge 1 ]] || gate_failed="$out_dir ($format): no annotated unit for $group_records input records across $group_files file(s)"
	else
		case "$format" in
		default | tab)
			[[ "$group_rows" -ge "$group_records" ]] || gate_failed="$out_dir ($format): $group_rows annotation rows for $group_records input records across $group_files file(s)"
			;;
		vcf)
			[[ "$group_rows" -ge 1 && "$group_written" -eq "$group_rows" ]] || gate_failed="$out_dir ($format): $((group_written - group_rows)) of $group_written written records carry no $VCF_INFO_FIELD field ($group_records input records, $group_files file(s))"
			;;
		json)
			[[ "$group_rows" -ge 1 && "$group_written" -eq "$group_rows" ]] || gate_failed="$out_dir ($format): $((group_written - group_rows)) of $group_written written objects carry no consequences array ($group_records input records, $group_files file(s))"
			;;
		esac
	fi
	if [[ -n "$gate_failed" ]]; then
		echo "ERROR: [generate_reference_output] $gate_failed" >&2
		echo "  A VEP run that cannot reach its cache emits warnings rather than annotations and still exits 0." >&2
		echo "  Check --perl-cache-dir, --cache-version and --assembly; a flag family that drops records or rows" >&2
		echo "  (--chr, --coding_only, --no_intergenic, --filter_common, --allow_non_variant, --individual," >&2
		echo "  --individual_zyg, --dont_skip)" >&2
		echo "  or writes one row per record (--pick, --most_severe, --summary) on an input with records vep skips" >&2
		echo "  takes --volume-floor any. The group's output was discarded." >&2
		for w in "$scratch"/*_warnings.txt; do
			[[ -f "$w" ]] || continue
			head -3 "$w" | sed 's/^/  warnings: /' >&2
			break
		done
		rm -rf "$scratch"
		exit 1
	fi

	while IFS='|' read -r in_vcf in_records out_units out_lines wall out_name; do
		[[ -n "$out_name" ]] || continue
		stem=$(sidecar_stem "$kind" "$format" "$out_name")
		mv -f "$scratch/$out_name" "$out_dir/$out_name"
		if [[ -f "$scratch/${out_name}_warnings.txt" ]]; then
			mv -f "$scratch/${out_name}_warnings.txt" "$out_dir/$stem.warnings.log"
		fi
		if [[ -f "$scratch/${out_name}_summary.txt" ]]; then
			mv -f "$scratch/${out_name}_summary.txt" "$out_dir/$stem.summary.log"
		fi
		if [[ "$kind" != "corpus" ]]; then
			write_provenance "$out_dir/$stem.provenance.json" "$species" "$asm" "$fork" "$format" "$headers" "$in_vcf" "$in_records" \
				"$out_dir/$out_name" "$out_units" "$out_lines" "$wall" "$IMAGE_DIGEST"
		fi
		done_jobs=$((done_jobs + 1))
		log "  wrote $out_dir/$out_name"
	done <<<"$group_facts"
	rm -rf "$scratch"
	log "group [$kind $asm $format] complete: $group_rows annotated units / $group_records records over $group_files file(s)"
done <<<"$GROUP_KEYS"
[[ "$done_jobs" -eq "$n_jobs" ]] || die "wrote $done_jobs of $n_jobs planned output(s); a group matched none of its records"

# A corpus: gzip the outputs and write provenance.json inside the scratch directory, then
# move everything into place, so the expected directory changes only once all of it exists
# (the golden tests read `<name>` before `<name>.gz`, so a plain copy of an earlier
# generation goes). The keys the script owns are rewritten; `input` and `fasta` are filled
# only when absent, since a record may describe them in prose and the owned `vep_command`
# carries the file names; any other key of an existing record is kept. The `sha256` map
# covers the outputs and each statistics report the run kept, so a committed report is
# pinned as the outputs are, and a run that kept none removes a report an earlier run left.
# Under --flag-set the same record is written under `flag_sets.<name>` and the top-level
# keys stand.
if [[ -n "$CORPUS_DIR" ]]; then
	for fmt in $(formats_selected); do
		out_name=$(out_name_for corpus "" "$fmt")
		gzip -n -9 -c "$CORPUS_STAGE/$out_name" >"$CORPUS_STAGE/$out_name.gz"
	done
	fasta_text="none on either side"
	fasta_arg=""
	if [[ -n "$CORPUS_FASTA_NAME" ]]; then
		fasta_text="$CORPUS_FASTA_NAME"
		fasta_arg="${CORPUS_FASTA_NAME%.gz}"
	fi
	vep_command="vep -i $CORPUS_INPUT -o <out> $(vep_args "$SPECIES" "$ASSEMBLY" "$FORK" "" '<cache>' "$fasta_arg" yes) [--tab|--vcf|--json]"
	python3 - "$CORPUS_DIR" "$CORPUS_STAGE" "$IMAGE_DIGEST" "$PERL_IMAGE" "$PERL_CACHE_VERSION" "$vep_command" \
		"$(cache_species_dir "$SPECIES")/${PERL_CACHE_VERSION}_${ASSEMBLY}" "$fasta_text" "$CORPUS_PLUGIN_FILES" "$PERL_ENV" \
		"$FLAG_SET" "$CORPUS_INPUT" "$CORPUS_EXPECTED_DIR" "$EXTRA_FLAGS" "$(formats_selected)" <<'EOF'
import hashlib, json, re, sys
from pathlib import Path

corpus, stage, image, tag, cache_version, command, perl_cache, fasta, plugin_files, perl_env = sys.argv[1:11]
flag_set, run_input, expected_dir, flags, formats = sys.argv[11:16]
corpus, stage = Path(corpus), Path(stage)
path = corpus / "provenance.json"
record = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
if not isinstance(record, dict):
    sys.exit(f"{path}: the record is not a JSON object")
outputs = {"default": "default.txt", "tab": "tab.txt", "vcf": "vcf.vcf", "json": "json.jsonl"}
names = [outputs[fmt] for fmt in formats.split()]
for fmt in formats.split():
    summary_log = f"{fmt}.summary.log"
    if (stage / summary_log).is_file():
        names.append(summary_log)
# The version in the shape of vep's own version line, which reference_release.py parses,
# read from the first header line of the default or the tab output.
reported = None
for name in ("default.txt", "tab.txt"):
    if (stage / name).is_file():
        head = (stage / name).read_text(encoding="utf-8", errors="replace").splitlines()[:1]
        version = re.match(r"## ENSEMBL VARIANT EFFECT PREDICTOR v(\S+)", head[0] if head else "")
        if version:
            reported = f"ensembl-vep          : {version.group(1)}"
            break
owned = {
    "vep_image": image,
    "vep_image_tag": tag,
    "vep_version_reported": reported,
    "cache_version": int(cache_version),
    "vep_command": command,
    "perl_cache": perl_cache,
    "fasta": fasta,
    "input": run_input,
    "sha256": {n: hashlib.sha256((stage / n).read_bytes()).hexdigest() for n in names},
}
# The container's environment, so a reference whose row text depends on Perl's hash
# order names the seed that reproduces it; a run that set none drops an earlier run's.
if perl_env:
    owned["perl_env"] = perl_env.split()
if flag_set:
    # A flag set's record: the manifest's own run keeps its keys, the set's run gets its own.
    sets = record.setdefault("flag_sets", {})
    if not isinstance(sets, dict):
        sys.exit(f"{path}: flag_sets is not a JSON object")
    sets[flag_set] = {**owned, "flags": flags.split(), "expected_dir": expected_dir}
    (stage / "provenance.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    sys.exit(0)
owned["variants_vcf_sha256"] = hashlib.sha256((corpus / "variants.vcf").read_bytes()).hexdigest()
# The plugin data files the reference was generated over, so a changed fixture shows;
# a manifest that names none drops an earlier run's digests.
if plugin_files:
    owned["plugin_files"] = {f: hashlib.sha256((corpus / f).read_bytes()).hexdigest() for f in plugin_files.split()}
else:
    record.pop("plugin_files", None)
if not perl_env:
    record.pop("perl_env", None)
filled_when_absent = ("fasta", "input")
for key, value in owned.items():
    if key in filled_when_absent:
        record.setdefault(key, value)
    else:
        record[key] = value
(stage / "provenance.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
EOF
	mkdir -p "$CORPUS_DIR/$CORPUS_EXPECTED_DIR"
	for fmt in $(formats_selected); do
		out_name=$(out_name_for corpus "" "$fmt")
		rm -f "$CORPUS_DIR/$CORPUS_EXPECTED_DIR/$out_name"
		mv -f "$CORPUS_STAGE/$out_name.gz" "$CORPUS_DIR/$CORPUS_EXPECTED_DIR/$out_name.gz"
		[[ "$KEEP_STATS" == "true" ]] || rm -f "$CORPUS_DIR/$CORPUS_EXPECTED_DIR/$fmt.summary.log"
		for side in warnings summary; do
			[[ -f "$CORPUS_STAGE/$fmt.$side.log" ]] && mv -f "$CORPUS_STAGE/$fmt.$side.log" "$CORPUS_DIR/$CORPUS_EXPECTED_DIR/$fmt.$side.log"
		done
		log "  wrote $CORPUS_DIR/$CORPUS_EXPECTED_DIR/$out_name.gz"
	done
	mv -f "$CORPUS_STAGE/provenance.json" "$CORPUS_DIR/provenance.json"
	log "  wrote $CORPUS_DIR/provenance.json"
fi

log "generate_reference_output complete: $done_jobs/$n_jobs output(s) written"
