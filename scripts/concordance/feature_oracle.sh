#!/usr/bin/env bash
#
# feature_oracle.sh -- the comparison oracle: the harness that runs Ensembl VEP and vep-rs
# on the same inputs and compares every column of every row in every output format. It runs
# the comparison cases of a case table: Ensembl VEP against one or more candidate `vep`
# binaries on the same input and flags, every requested output format compared column by
# column by vep-formats-compare, one report per candidate and case, one status table per
# run.
#
# THE CASE TABLE is a tab-separated file whose header names these columns, in any order.
# A cell with no value is `-`; blank lines and `#` lines are skipped.
#
#   case               unique name: the report directory and the status.tsv key
#   tier               fixture | population; --tier selects the rows of one tier
#   input              the input file: absolute, else relative to --data-dir, else to the
#                      repository root
#   cache              <version>_<assembly>, the Ensembl cache slice both engines annotate
#                      from: <species>/<version>_<assembly> under --perl-cache-dir
#                      (<species>_refseq or _merged when the flags say so) and
#                      <assembly, lower-cased> under --json-cache-dir; an optional
#                      =<dir> suffix names the candidate's JSON cache instead, and an
#                      =layout:<root> suffix a cache root in Ensembl VEP's layout
#                      (<root>/<species>[_refseq|_merged]/<version>_<assembly>), which
#                      the candidate resolves itself from --dir_cache <root> --species
#                      <species> --cache_version <version> --assembly <assembly> in
#                      place of --json_cache; <dir> and <root> resolve like input
#   species            the Ensembl species both engines run under (--species, the cache
#                      slice's directory: CacheDir.pm:306-315 at release/116.3 reads
#                      <dir_cache>/<species>[_refseq|_merged]/<version>_<assembly>), or -
#                      for homo_sapiens
#   fasta              - or an indexed FASTA given to both engines (its .fai beside it, a
#                      compressed file's as <name>.fa.gz.fai or <name>.fa.fai; a
#                      .gz is decompressed once per run, as the golden tests do, into a
#                      directory of its own under the work tree, so two sources that
#                      decompress to one file name never share a copy), or one of two
#                      marked forms for a FASTA that comes without its index, which each
#                      engine writes itself (FastaSequence.pm:307-359).
#                      unindexed:<path> gives every engine run, one per engine and output
#                      format, a copy of the file alone in a directory of its own under the
#                      case's work directory, made right before the run: the candidate is
#                      given the copy as --fasta, Ensembl VEP its directory as a writable
#                      mount and --fasta inside it. A .gz is copied as it is and never
#                      decompressed. For each output format whose two runs exited 0, each
#                      engine's directory must hold <name>.fai, and <name>.gzi beside a
#                      BGZF file, neither <name>.vep.lock (FastaSequence.pm:312-313) nor a
#                      *.tmp file, and each index must be byte for byte the other
#                      engine's; a finding is a red reason of its own, `<format>: fasta
#                      index: <finding>`, beside the format's comparison.
#                      unindexed-readonly:<path> gives the same copies in directories no
#                      engine can write (mode 0555, Ensembl VEP's mounted read-only), where
#                      Ensembl VEP refuses to run (FastaSequence.pm:334), so the row names
#                      a refusal directory; after the runs each directory must hold the
#                      copy alone. As user 0 the candidate of such a case runs through
#                      `setpriv --bounding-set=-dac_override`, since a directory's mode
#                      binds no process that holds the capability. A marked case runs
#                      Ensembl VEP live under --reference-dir as well and its reference
#                      cell is -: a stored output holds no index, so none is read for the
#                      case. Its masks name no class the comparator decides from the
#                      FASTA, which it reads through an index.
#   flags              whitespace-separated flags appended verbatim to both engines
#                      (--genomes to Ensembl VEP alone), or -;
#                      a flag is read in every spelling vep accepts (-nearest, --near,
#                      --NEAREST name --nearest); a --nearest case mounts the Ensembl
#                      cache writable, since Ensembl VEP writes transcript_gene_tss.txt
#                      there on its first run (a scan of every region file) and keeps the
#                      file, built under that run's transcript filters, for every later
#                      run; while a slice lacks the file, its --nearest runs take turns
#                      under one lock, because the temporary file is named after the
#                      process id, which is the same in every container; its comparisons
#                      settle NEAREST by start site through the candidate's JSON cache
#                      (vep-formats-compare --nearest-sites); a word may carry
#                      `{corpus}`, the input's directory, or `{sources}`, the
#                      --sources-dir directory: the candidate reads the host path and
#                      Ensembl VEP a read-only mount of it (/work/corpus, /work/sources),
#                      so a flag names a file beside the input or a shared source without
#                      an absolute path; beside a manifest.json a word may carry
#                      `{plugins}` or `{generated}` as the manifest's runs write them,
#                      read as `{corpus}/plugin_data` and `{corpus}/generated`, so the
#                      cell is that run of the manifest, and without a manifest.json
#                      beside the input either word is refused; a cell may name --dir, as
#                      `--dir <directory>` or `--dir=<directory>`, with the directory
#                      written as `{corpus}` or `{sources}` or a path under one, the
#                      placeholder standing in the first segment alone and no segment
#                      being `..`, since a flags cell reaches no other directory of
#                      Ensembl VEP's container; the mount is
#                      read-only, so a file a plugin keeps under --dir is there before
#                      the run (the GWAS plugin's processed catalogue, GWAS.pm:106-118,
#                      which a first run would write); both commands keep the runner's
#                      own cache words, which outrank --dir as the cache root
#                      (AnnotationSourceAdaptor.pm:127 at release/116.3), and a
#                      directory holding a vep.ini is refused, as --config is, since
#                      both engines read that file as options (Config.pm:704-705)
#   formats            comma-separated subset of default,tab,vcf,json,parquet,warnings,stats;
#                      warnings compares the `_warnings.txt` files of one run of each
#                      engine: the default run's where the cell names default, else those
#                      of the one output format the cell names, which is tab, vcf or json
#                      (`json,warnings` compares the JSON runs' files; a cell naming
#                      neither default nor one such format beside warnings is refused);
#                      stats compares the default runs' statistics reports
#                      through vep-stats-compare, so it needs default too, and a flags
#                      cell naming --stats_file, which moves the report, is refused
#   perl_image         the Ensembl VEP image, pinned by digest (name@sha256:...)
#   perl_mode          offline | database
#   reference_release  the Ensembl VEP release the image carries, recorded per report
#   masks              vep-formats-compare --masks, applied from the manifest.json
#                      beside the input (whose flag set over the input is selected by the
#                      input's name, which every format's comparison passes), or -; the
#                      comparator's built-in classes
#                      (BUILTIN_MASKS below: regulatory_feature_rows, order_only,
#                      individual_allele_order, conservation_short_name_spelling,
#                      pc_description_text, motif_score_reverse_strand_mnv, the registry
#                      pair classes start_cooccurrence_swap and splice_family_swap, the
#                      field, row and warnings classes of the corrected behaviours, and the
#                      structural-variant comparator's classes) need no manifest, but for
#                      one field class: dbnsfp_mitochondrial_record_named_m stands outside
#                      that list and applies beside a manifest.json alone; a class of
#                      MANIFEST_ONLY_MASKS applies from a manifest's list alone, and a case
#                      naming one with no manifest.json beside its input fails preflight
#                      as such. A class
#                      may carry its expected applied count, `<class>=<n>` (n at least 1):
#                      the case is red when an output format's report (not the warnings
#                      report, whose count is lines) applies the class to another number of
#                      entries, and a class with or without a count is red when every
#                      compared format, warnings included, leaves it unused, since a mask
#                      that sets nothing aside documents nothing; a class an output
#                      format's report names among those asked for and leaves out of its
#                      unused list is used though nothing is counted under it (the class
#                      a value set aside under another class needs named beside it).
#                      Beside a manifest.json the comparator is told two things of the
#                      candidate's run. Its run of the manifest (vep-formats-compare
#                      --run), whose lists alone are then read: the run the case's row of
#                      the named-runs table gives (THE NAMED-RUNS TABLE below), else the
#                      one run, the manifest's own (`default`) or a flag set, whose input
#                      has the case input's file name and whose flags are the flags cell's
#                      words in any order, the manifest's side read without --grch37 and
#                      without a --fasta word and its file, which the cache and fasta cells
#                      give. Both sides are read without --dir_plugins and its directory,
#                      with an option under any name an engine reads as it (PERL_OPTIONS,
#                      CANDIDATE_OPTIONS) and with a path under the corpus in either
#                      spelling of RUN_SPELLINGS. A case that is no run of its manifest,
#                      or fits several, names none, and an unnamed comparison reads the
#                      lists of every run over its input as one: where the manifest lists
#                      more than one run over the input and no row of the table names
#                      the case, it fails preflight by name before any comparator is
#                      started, as does a case whose match fails.
#                      And whether the run has a reference
#                      sequence (--reference-sequence), as the candidate settles it: present
#                      under a fasta cell, which it is given as --fasta; else present where
#                      its cache directory holds a .fa or .fa.gz file with the .fai beside
#                      it, and the .gzi as well where the file's bytes are BGZF (a gzip file
#                      that is not BGZF is never taken), unless the flags name --no_fasta,
#                      which keeps the candidate from that file; else absent. A flags cell
#                      naming --fasta_dir says neither: the candidate takes that
#                      directory's file by the species and assembly in its name, and the
#                      cache directory's where none qualifies
#   floors             vep-formats-compare --floors (COL=MIN,...), or -; a floor on a
#                      column a format cannot count (the json format's record-level columns,
#                      Feature_type, Extra and a custom source's columns; a tab or vcf column
#                      neither header declares) is not applicable there: the report.json entry
#                      of that format lists it under floors_not_applicable, and a case whose
#                      formats give such a floor no format that counted it is red
#   field_pairs        vep-formats-compare --field-pairs (COL=EXPECTED=>ACTUAL;...), or -
#   unordered          - or 1: the comparator takes --unordered, so the shared rows agree in
#                      any order, for a run whose rows Ensembl VEP emits in hash order (the
#                      allele buckets of --pick_allele and --pick_allele_gene, the gene
#                      buckets of --per_gene)
#   reference          - or a directory of stored Ensembl VEP outputs for this case, resolved
#                      like input: {default.txt,tab.txt,vcf.vcf,json.jsonl}[.gz] for the
#                      case's formats (the layout of a golden corpus's expected directory),
#                      plus, when the case names warnings, the `_warnings.txt[.gz]` file of
#                      the run whose warnings it compares (default.txt_warnings.txt[.gz],
#                      or json.jsonl_warnings.txt[.gz] for `json,warnings`; an
#                      empty file states no warnings) and default.summary.log[.gz] when it
#                      names stats; the case runs no container, a
#                      provenance.json beside the directory is recorded, and preflight
#                      refuses a directory lacking a format the case names; a tool row's
#                      cell may name the stored filtered output file itself (a corpus's
#                      expected/filter/<case>.<ext>), whose provenance.json is the
#                      corpus's, beside the input
#   refusal            - or a directory, resolved like input, holding perl.txt and
#                      candidate.txt: the case expects both engines to refuse the run, and
#                      is green when, in every requested format, each engine exits non-zero
#                      with its file's text (without the trailing newline) in its stderr,
#                      red otherwise; nothing is compared, so the masks, floors,
#                      field_pairs, unordered and reference cells are -
#   tool               - or filter_vep:<args>: the case compares the engines' companion
#                      filter tools in place of the engines. The input is an annotated
#                      output (a VCF, tab or default-format file), and each side filters it
#                      with <args>, shell words (a filter string is quoted): Ensembl VEP's
#                      filter_vep runs in the image, the candidate's `vep filter` beside
#                      it, and the two filtered outputs are compared in the row's one
#                      format (default, tab or vcf). A --count run writes the count alone,
#                      so the two counts are compared as text and the comparator cells
#                      are -. The runner sets -i, -o and --force_overwrite, and in database
#                      mode filter_vep's connection (--host, --port, --user, --version
#                      <cache version>; filter_vep answers --ontology from Ensembl's
#                      ontology database, so an --ontology row runs in database mode), so
#                      <args> names none of them and not --list, whose field list goes to
#                      standard output. The flags and fasta cells are -: no engine runs,
#                      and the cache cell is recorded and otherwise unread.
#
# THE NAMED-RUNS TABLE (--case-runs, else a file oracle_case_runs.tsv beside the case table
# where one is there) names the run of its manifest a case compares under where the flags
# cell is no run of it: a tab-separated file whose header is `case`, `run`, `reason`, one
# row per case, blank lines and `#` lines skipped. `run` is a run of the manifest beside the
# case's input whose input is the case's (`default`, or a flag set's name), or `pooled`,
# which keeps the unnamed comparison on purpose (the comparator refuses it where the runs'
# lists document a class the masks cell names); `reason` says why. A named row outranks the
# flags cell. Refused before anything runs: a second row for a case, a row for a case the
# case table lacks and, for a selected case, a run the manifest lacks, a run over another
# input, a `pooled` row on a case whose flags are a run of its manifest, and a named run on
# a case that compares against no manifest (none beside its input, a refusal row or a
# --count tool row); a `pooled` row there is not read.
#
# THE ON-DEMAND LIST (oracle_cases_on_demand.tsv beside the case table, where that file is
# there) names the population cases that run only when --only names them: a tab-separated
# file whose header is `case`, `reason`, blank lines and `#` lines skipped. A run without
# --only leaves the listed cases of its tier out and logs one line naming them; --only runs
# a listed case like any other. Refused before anything runs: a row for a case the case
# table lacks, and a row for a fixture case.
#
# THE ENGINES. Ensembl VEP runs in Docker, once per case and output format, since a run
# writes one format (Runner.pm:871 at release/116.3 builds one OutputFactory per run):
#
#   vep -i <input> -o <out> --offline --cache --dir_cache <cache> --cache_version <v>
#       --species <species> --assembly <asm> --format vcf --fork 1 --no_stats --quiet
#       --force_overwrite <flags> [--fasta <fa>] [--tab|--vcf|--json]
#
# In database mode `--offline --cache --dir_cache --cache_version` become
# `--database --db_version <v>`: Config.pm:618 refuses --database beside --offline or
# --cache, and the HGVS, ID and SPDI input parsers (Parser/HGVS.pm:104, Parser/ID.pm:100,
# Parser/SPDI.pm:96) refuse to run offline. The server is Ensembl's public one, the
# Config.pm defaults (host ensembldb.ensembl.org, port 3306; Config.pm:711 moves GRCh37
# to 3337; a --genomes flag moves every species to Ensembl Genomes,
# mysql-eg-publicsql.ebi.ac.uk:4157, Config.pm:391-394), recorded in the report. The
# candidate runs with the same flags, less --genomes (it reads a local cache and refuses
# the flag by name), and
# `--offline --json_cache <cache>`; its parquet run is `--output_format parquet`,
# compared with the Ensembl VEP tab output. `--format vcf` is given to both engines
# unless the flags name an input format. The warnings compared are the default runs'
# `<out>_warnings.txt`, or those of the one format a cell names beside warnings without
# default: a run's warnings follow its output format (the class
# json_mismatch_warning_names_no_record stands in a JSON run's file alone). Ensembl VEP
# opens the file on the first warning
# (BaseVEP.pm:878-897), so an absent file on either side means no warnings. A warning an
# Ensembl VEP worker raises while annotating is written with the worker's pid under
# --fork (Runner.pm:563 prefixes the captured worker stderr with it), so it differs from
# run to run; a record skipped while parsing is written by the parent without one. The
# statistics compared are the default runs' `<out>_summary.txt`: on a case naming stats each
# engine's default run takes --stats_text in place of --no_stats (BaseRunner.pm:254-271
# names the file; Runner.pm:787-789 adds the HTML report only where neither format is
# asked for), every other run keeps --no_stats, and a default run that leaves no report is
# red for stats.
#
# --reference-dir replaces the Ensembl VEP runs of every selected case with a stored
# reference: <reference-dir>/<case>/expected/, in the layout the reference column
# describes, with the provenance.json beside expected/ recorded. A case's own reference
# cell names its directory instead. No container runs for a case with a stored reference,
# and docker is needed only while a selected case has none; a refusal case runs both
# engines live under either, and so does a case whose fasta cell is marked.
#
# A CASE IS GREEN for a candidate when every requested format's comparison is green:
# no entry on one side only, every compared field equal after the masks and field
# pairs, every floor met, for warnings the same (line, description) set, and for stats
# the same sections and rows by vep-stats-compare's rules, which set aside the rows no
# second run reproduces; the masks cell is not given to that tool. A format
# whose engine or comparison fails to run is red, and so is one whose two outputs both
# carry no annotation row: two empty outputs agree on nothing. A refusal case is green
# when both engines refuse every requested format with the expected texts. A case whose
# fasta cell is marked is red as well for any finding of the checks the fasta column
# describes. The run exits 1 on any red.
#
# A DISCOVERY PASS (--discover) reports which classes apply to each selected case and gives
# no verdict. In place of the case's masks cell the comparator is asked for the registry:
# every class of BUILTIN_MASKS but those of DISCOVERY_OPT_IN_MASKS, the classes the case's
# own manifest.json documents, and the classes the cell names, which are always asked, so
# that the rows a named class covers are never read as differences. --discover-classes
# <class>[,<class>...] narrows the registry to the listed classes; the cell's own stay. A
# tool row and a refusal row are asked for their own cell alone: a tool row compares two
# filtered copies of one annotated output, which no class of the engines' annotation
# reaches, and a refusal row compares nothing. A case its row of the named-runs table keeps
# pooled, on a manifest of more than one run over its input, is asked for no class the
# manifest documents but those its cell names: the comparator refuses an unnamed comparison
# asked for a class the lists of several runs over the input document, as it refuses such
# a cell in a verdict run. The cell's expected counts are not held, and the statistics
# comparison of a case naming stats runs as in a verdict run, its differences counted
# among the residual ones, as the findings of a marked fasta cell's checks are.
#
# Per candidate and case the pass writes discovery.tsv beside report.json: format, class,
# rows (the entries the class set aside in that format), named (whether the cell names the
# class) and applied (whether the class applied in any compared format); a class the cell
# names that applied nowhere is one line whose format is `-`. The run writes
# <out>/discovery.tsv (case, candidate_pin, class, applied_formats, rows: one line per class
# that applied, the two lists in step; a case no class applied to is one line whose class is
# `-`) and one summary line per case: the classes applied, the classes its cell names that
# applied nowhere, the classes that applied and its cell does not name, and the count of
# residual differences, those no class asked for covers, which report.json lists.
#
# A discovery pass lands no verdict. A case's status in status.tsv is `discovered`, or
# `incomplete` where an engine or a comparison of the case did not run, never green or red;
# report.json carries a `status` and no `verdict`; run.json carries `"mode": "discover"`;
# and the run exits 0 when every selected case ran, whatever applied and whatever differed.
# A residual difference does not fail the pass: the pass says which classes a case's cell
# should name, and the verdict pass after it says whether the case is exact.
#
# OUTPUT
#   <out>/status.tsv                           case, candidate_pin, verdict,
#                                              report_digest (sha256 of report.json),
#                                              run_utc; one block per candidate
#   <out>/<candidate_pin>/<case>/report.json   the case record with its verdict
#   <out>/<candidate_pin>/<case>/<format>.json vep-formats-compare's report
#   <out>/<candidate_pin>/<case>/stats.json    vep-stats-compare's report, on a stats case
#   <out>/run.json                             the run's options and image digests,
#                                              its mode, verdict or discover, and what
#                                              sizing a machine for the run needs (below)
#   <out>/discovery.tsv                        a discovery pass's table of the classes
#   <out>/<candidate_pin>/<case>/discovery.tsv that applied, run-wide and per case
#   <work-dir>/cases/<case>/                   the engines' outputs (perl/ and
#                                              candidates/<pin>/) and the staged
#                                              copies of a marked fasta cell
#                                              (fasta/perl/<format>/ and
#                                              fasta/candidates/<pin>/<format>/),
#                                              deleted per case unless --keep-outputs;
#                                              fasta/ and logs/ beside cases/
#
# SIZING. Every engine run and every comparison is timed by a helper that waits for it and
# records its wall seconds, the peak resident set of the largest process it ran (wait4's
# ru_maxrss, the figure `/usr/bin/time -v` prints) and the used bytes of the work directory's
# filesystem, read every two seconds while the step runs. report.json carries them as
# `sizing`: comparison_program and comparison_program_version (the format comparator's name,
# vep-formats-compare, and the first 12 hex of the sha256 of the file that ran as it);
# rows_compared, per compared format the larger of the two outputs' entry counts
# in the comparator's report; peak_memory_kib, the largest peak of the candidate's runs and
# the comparisons; work_disk_peak_bytes, the filesystem's high-water mark over the case's
# steps, and work_disk_growth_bytes, that mark less the fewest used bytes a step of the
# case started at; wall_minutes, the case's start to its report; and steps, each step's own
# figures under ensembl-vep.<format>, candidate.<format> and comparison.<format>, the
# warnings and the statistics comparisons under comparison.warnings and comparison.stats.
# An Ensembl VEP step's memory reads `not measured`: the engine runs in a container, outside
# the timed process. A figure no step left reads `not measured` too, and fails nothing.
# run.json lists the same block per case and candidate under `case_sizes`, beside the
# program's name and version, so the machine of a later run can be sized from this one's.
#
# status.tsv, run.json and the run's discovery.tsv describe the run that wrote them: all
# three are removed when a run starts, so an interrupted run leaves no table, a verdict run
# leaves no discovery table of an earlier pass, and a run over some cases of a table
# indexes those cases alone.
#
# Usage:
#   scripts/concordance/feature_oracle.sh --cases <table.tsv> --out <dir> \
#     --candidates <vep>[,<vep>...] --candidate-pins <pin>[,<pin>...] \
#     [--tier fixture|population] [--only <case>[,...]] [--data-dir <dir>] \
#     [--perl-cache-dir <dir>] [--json-cache-dir <dir>] [--reference-dir <dir>] \
#     [--sources-dir <dir>] [--work-dir <dir>] [--keep-outputs] [--fork <N>] \
#     [--perl-fork <N>] [--jobs <N>] [--format-jobs <N>] [--dry-run] \
#     [--discover [--discover-classes <class>[,<class>...]]] [--case-runs <tsv>]
#
#   --format-jobs N runs a candidate's output formats of one case (its run and its
#   comparison) as N background jobs waited on, each logging to <fmt>.job.out and
#   <fmt>.job.err beside its outputs; the per-format verdicts are merged as under the
#   default, one format after another, so report.json reads the same either way.
#
#   --dry-run runs every preflight check and prints each engine and comparator command
#   instead of running them. A missing input fails a dry run as it fails a real one.
#
#   VEP_FORMATS_COMPARATOR names the format comparator, a build of
#   tools/vep-formats-compare, as a path or a command on PATH. Unset or empty: this tree's
#   release build, tools/vep-formats-compare/target/release/vep-formats-compare (cargo
#   build --release --locked --manifest-path tools/vep-formats-compare/Cargo.toml). The
#   runner builds nothing: a comparator that is not an executable stops the run, a dry run
#   as well, before the first engine run.
#
#   VEP_STATS_COMPARATOR names the statistics comparator, a build of
#   tools/vep-stats-compare, the same two ways; unset or empty, this tree's release build.
#   It is looked for only where a selected case names stats, once preflight has passed.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# The program a sizing figure belongs to.
COMPARISON_PROGRAM="vep-formats-compare"
PYTHON="${PYTHON:-python3}"

CASES=""
OUT_DIR=""
CANDIDATES_SPEC=""
PINS_SPEC=""
TIER="fixture"
ONLY_SPEC=""
DATA_DIR=""
PERL_CACHE_DIR=""
JSON_CACHE_DIR=""
REFERENCE_DIR=""
SOURCES_DIR=""
WORK_DIR=""
KEEP_OUTPUTS="false"
FORK=""
PERL_FORK=1
JOBS=1
FORMAT_JOBS=1
DRY_RUN="false"
DISCOVER="false"
DISCOVER_CLASSES=""
DISCOVER_CLASSES_GIVEN="false"
CASE_RUNS=""

usage() {
	cat <<'EOF'
Run a case table's comparison cases: Ensembl VEP against candidate vep binaries.

Required:
  --cases <tsv>             The case table (columns in the header of this script)
  --out <dir>               Where status.tsv and the per-case reports are written
  --candidates <bin>[,...]  Candidate vep binaries, each run on every selected case
  --candidate-pins <p>[,...] One pin per candidate (the commit or tag it was built
                            from): the report directory and status.tsv's candidate_pin

Optional:
  --tier <fixture|population>  The tier to run (default fixture)
  --only <case>[,...]    Restrict the run to these cases of the tier; a case of
                            oracle_cases_on_demand.tsv beside the case table runs only
                            when named here
  --data-dir <dir>          Root for relative inputs and the cache defaults
                            (default: the repository root)
  --perl-cache-dir <dir>    Directory holding <species>/<version>_<assembly>/
                            (default <data>/caches/perl/vep-cache)
  --json-cache-dir <dir>    Directory holding the candidate's <assembly>/ JSON caches
                            (default <data>/caches/vep-rs)
  --reference-dir <dir>     Stored Ensembl VEP output per case
                            (<dir>/<case>/expected/); a case's reference cell names
                            its own; no container runs for a case with one
  --sources-dir <dir>       The directory a {sources} flag word names: the candidate
                            reads it there, Ensembl VEP at a read-only mount
  --work-dir <dir>          Where the engines write (default <out>/work)
  --keep-outputs            Keep the engines' outputs after the comparison
  --fork <N>                The candidate's --fork (default: the binary's own)
  --perl-fork <N>           Ensembl VEP's --fork (default 1)
  --jobs <N>                Cases run at once (default 1)
  --format-jobs <N>         A candidate's output formats of one case run and compared
                            at once (default 1)
  --dry-run                 Preflight everything, print the commands, run nothing
  --discover                Report which classes apply to each case and give no verdict:
                            the comparator is asked for the registry in place of each
                            masks cell, a case's status is `discovered`, and the run
                            exits 0 when every selected case ran
  --discover-classes <c>[,...]  With --discover: the classes asked for, beside each
                            cell's own, in place of the registry
  --case-runs <tsv>         The named-runs table: case, run, reason (default
                            oracle_case_runs.tsv beside the case table, where it is there)
EOF
	exit 1
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--cases)
		CASES="$2"
		shift 2
		;;
	--out)
		OUT_DIR="$2"
		shift 2
		;;
	--candidates)
		CANDIDATES_SPEC="$2"
		shift 2
		;;
	--candidate-pins)
		PINS_SPEC="$2"
		shift 2
		;;
	--tier)
		TIER="$2"
		shift 2
		;;
	--only)
		ONLY_SPEC="$2"
		shift 2
		;;
	--data-dir)
		DATA_DIR="$2"
		shift 2
		;;
	--perl-cache-dir)
		PERL_CACHE_DIR="$2"
		shift 2
		;;
	--json-cache-dir)
		JSON_CACHE_DIR="$2"
		shift 2
		;;
	--reference-dir)
		REFERENCE_DIR="$2"
		shift 2
		;;
	--sources-dir)
		SOURCES_DIR="$2"
		shift 2
		;;
	--work-dir)
		WORK_DIR="$2"
		shift 2
		;;
	--keep-outputs)
		KEEP_OUTPUTS="true"
		shift
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--perl-fork)
		PERL_FORK="$2"
		shift 2
		;;
	--jobs)
		JOBS="$2"
		shift 2
		;;
	--format-jobs)
		FORMAT_JOBS="$2"
		shift 2
		;;
	--dry-run)
		DRY_RUN="true"
		shift
		;;
	--discover)
		DISCOVER="true"
		shift
		;;
	--discover-classes)
		DISCOVER_CLASSES="$2"
		DISCOVER_CLASSES_GIVEN="true"
		shift 2
		;;
	--case-runs)
		CASE_RUNS="$2"
		shift 2
		;;
	-h | --help) usage ;;
	*)
		echo "ERROR: [feature_oracle] unknown argument: $1" >&2
		usage
		;;
	esac
done

die() {
	echo "ERROR: [feature_oracle] $*" >&2
	exit 1
}
log() {
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
}
# `${var,,}` is bash 4+ and fails at runtime on bash 3.2 while passing `bash -n`.
lower() { printf %s "$1" | tr '[:upper:]' '[:lower:]'; }
abs_path() {
	local d
	d=$(cd "$(dirname "$1")" && pwd) || return 1
	printf '%s/%s' "$d" "$(basename "$1")"
}
# norm_path <path>: absolute, without `.`, `..`, repeated or trailing slashes; the path
# need not exist (symbolic links are not resolved).
norm_path() {
	local p="$1" out="" seg
	[[ "$p" == /* ]] || p="$PWD/$p"
	# The unquoted $p splits on IFS alone: set -f keeps a `*` or `?` in a segment from
	# matching the working directory's entries.
	local IFS='/'
	set -f
	for seg in $p; do
		case "$seg" in
		'' | '.') ;;
		'..') out="${out%/*}" ;;
		*) out="$out/$seg" ;;
		esac
	done
	set +f
	printf '%s' "${out:-/}"
}
sha256_of() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	elif command -v shasum >/dev/null 2>&1; then
		shasum -a 256 "$1" | awk '{print $1}'
	else
		echo "unavailable"
	fi
}
# Whether the whitespace-separated list $1 contains the word $2.
has_word() { [[ " $1 " == *" $2 "* ]]; }
# in_lines <list> <item>: whether the newline-separated <list> has a line equal to <item>.
in_lines() { printf '%s' "$1" | grep -Fx -- "$2" >/dev/null; }
# tsv_lookup <table> <key>: the second column of the first line of the two-column
# tab-separated <table> whose first column is <key>, or nothing. The table is awk's
# input as a here-string, here and wherever an awk leaves at its first match: a pipe's
# writer that is still writing when its reader leaves is ended by the closed pipe, which
# pipefail makes fatal.
tsv_lookup() {
	awk -F'\t' -v k="$2" '$1 == k {print $2; exit}' <<<"$1"
}

# Validate the options
[[ -n "$CASES" ]] || {
	echo "ERROR: [feature_oracle] --cases required" >&2
	usage
}
[[ -n "$OUT_DIR" ]] || {
	echo "ERROR: [feature_oracle] --out required" >&2
	usage
}
[[ -n "$CANDIDATES_SPEC" ]] || {
	echo "ERROR: [feature_oracle] --candidates required" >&2
	usage
}
[[ -n "$PINS_SPEC" ]] || {
	echo "ERROR: [feature_oracle] --candidate-pins required (one pin per candidate)" >&2
	usage
}
[[ -f "$CASES" ]] || die "case table not found: $CASES"
case "$TIER" in
fixture | population) ;;
*) die "invalid --tier: $TIER (fixture|population)" ;;
esac
for _n in "$PERL_FORK" "$JOBS" "$FORMAT_JOBS" ${FORK:+"$FORK"}; do
	[[ "$_n" =~ ^[0-9]+$ && "$_n" -gt 0 ]] || die "--fork, --perl-fork, --jobs and --format-jobs take a positive integer (got '$_n')"
done
# A list of classes is a discovery pass's: a verdict run asks for each case's own cell. An
# empty list is refused, not read as no list: it would ask every case for the registry.
if [[ "$DISCOVER_CLASSES_GIVEN" == "true" ]]; then
	[[ "$DISCOVER" == "true" ]] || die "--discover-classes names the classes a discovery pass asks for; pass --discover with it"
	[[ "$DISCOVER_CLASSES" =~ ^[a-z0-9_]+(,[a-z0-9_]+)*$ ]] || die "--discover-classes '$DISCOVER_CLASSES' is not <class>[,<class>...] (lower-case letters, digits and underscores; a list carries no expected count)"
fi
[[ -n "$DATA_DIR" ]] || DATA_DIR="$REPO_DIR"
[[ -d "$DATA_DIR" ]] || die "data-dir not a directory: $DATA_DIR"
DATA_DIR=$(abs_path "$DATA_DIR")
[[ -n "$PERL_CACHE_DIR" ]] || PERL_CACHE_DIR="$DATA_DIR/caches/perl/vep-cache"
[[ -n "$JSON_CACHE_DIR" ]] || JSON_CACHE_DIR="$DATA_DIR/caches/vep-rs"
# Docker binds the Ensembl cache root by host path, which must be absolute; the report
# records both roots.
[[ ! -d "$PERL_CACHE_DIR" ]] || PERL_CACHE_DIR=$(abs_path "$PERL_CACHE_DIR")
[[ ! -d "$JSON_CACHE_DIR" ]] || JSON_CACHE_DIR=$(abs_path "$JSON_CACHE_DIR")
if [[ -n "$REFERENCE_DIR" ]]; then
	[[ -d "$REFERENCE_DIR" ]] || die "reference-dir not a directory: $REFERENCE_DIR"
	REFERENCE_DIR=$(abs_path "$REFERENCE_DIR")
fi
if [[ -n "$SOURCES_DIR" ]]; then
	[[ -d "$SOURCES_DIR" ]] || die "sources-dir not a directory: $SOURCES_DIR"
	SOURCES_DIR=$(abs_path "$SOURCES_DIR")
fi
# The flag-word placeholders and the container paths they take in the Ensembl VEP
# command; the candidate takes the host directories.
CORPUS_PLACEHOLDER='{corpus}'
SOURCES_PLACEHOLDER='{sources}'
CORPUS_MOUNT="/work/corpus"
SOURCES_MOUNT="/work/sources"
# The words a corpus manifest's runs write for two directories of the corpus: {plugins}, its
# bundled plugin files (<corpus>/plugin_data), and {generated}, what its `generated` block
# writes (<corpus>/generated, which stands there before the run: the runner writes neither
# directory). A flags cell written as such a run writes it is that run of the manifest, so
# each word is read, beside a manifest.json, on the path {corpus} takes for the engine.
PLUGINS_PLACEHOLDER='{plugins}'
GENERATED_PLACEHOLDER='{generated}'
# The format comparator: the build of tools/vep-formats-compare that VEP_FORMATS_COMPARATOR
# names, else this tree's release build. `type -P` resolves a bare name through PATH as the
# exec does (`-x` alone would test it against the cwd) and returns a non-executable PATH
# match, hence the file tests on its result.
COMPARATOR_NAMED="${VEP_FORMATS_COMPARATOR:-$REPO_DIR/tools/vep-formats-compare/target/release/vep-formats-compare}"
if COMPARATOR=$(type -P -- "$COMPARATOR_NAMED") && [[ -f "$COMPARATOR" && -x "$COMPARATOR" ]]; then
	COMPARATOR=$(abs_path "$COMPARATOR")
else
	die "the format comparator $COMPARATOR_NAMED is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-formats-compare/Cargo.toml' or name a built one in VEP_FORMATS_COMPARATOR"
fi
command -v "$PYTHON" >/dev/null 2>&1 || die "$PYTHON is required (the report writer runs on it)"

# Candidates and pins: one pin per binary, pins unique and usable as directory names.
CANDIDATES=()
PINS=()
IFS=',' read -r -a CANDIDATES <<<"$CANDIDATES_SPEC"
IFS=',' read -r -a PINS <<<"$PINS_SPEC"
[[ "${#CANDIDATES[@]}" -eq "${#PINS[@]}" ]] || die "--candidates names ${#CANDIDATES[@]} binaries and --candidate-pins ${#PINS[@]} pins; they go one to one"
PIN_RE='^[A-Za-z0-9][A-Za-z0-9._-]*$'
[[ -n "$WORK_DIR" ]] || WORK_DIR="$OUT_DIR/work"
# The names the run itself writes under --out; a pin is a directory there.
RESERVED_PINS="status.tsv run.json discovery.tsv"
NORM_WORK_DIR=$(norm_path "$WORK_DIR")
if [[ "${NORM_WORK_DIR%/*}" == "$(norm_path "$OUT_DIR")" ]]; then
	RESERVED_PINS="$RESERVED_PINS ${NORM_WORK_DIR##*/}"
fi
for ((i = 0; i < ${#CANDIDATES[@]}; i++)); do
	bin="${CANDIDATES[$i]}"
	pin="${PINS[$i]}"
	[[ -n "$bin" && -f "$bin" && -x "$bin" ]] || die "candidate is not an executable file: '$bin'"
	CANDIDATES[i]=$(abs_path "$bin")
	[[ "$pin" =~ $PIN_RE ]] || die "candidate pin '$pin' is not a directory-safe name ([A-Za-z0-9._-], not starting with . or -)"
	has_word "$RESERVED_PINS" "$pin" && die "candidate pin '$pin' names a file or directory the run writes under --out ($RESERVED_PINS)"
	for ((j = 0; j < i; j++)); do
		[[ "${PINS[$j]}" != "$pin" ]] || die "candidate pin '$pin' is given twice; every candidate needs its own pin"
	done
done

# The case table. Its column names stand under a name of the runner's own: bash sets COLUMNS
# and LINES to the size of a terminal on standard error after each command a script waits for.
TABLE_COLUMNS="case tier input cache fasta flags formats perl_image perl_mode reference_release masks floors field_pairs species unordered reference refusal tool"
N_TABLE_COLUMNS=18
DEFAULT_SPECIES="homo_sapiens"
OUTPUT_FORMATS="default tab vcf json parquet"
# The companion tools a tool cell may name (candidate_cmd maps each to the candidate's
# subcommand) and the formats a tool reads (detect_format, filter_vep:423-438 at
# release/116.3, reads a VCF or a tab file, the default format among the latter).
TOOLS="filter_vep"
TOOL_FORMATS="default tab vcf"
# Tool words the runner sets (the input, the output, its overwrite, the database connection
# of filter_vep's --ontology) or whose output never reaches the output file (--list prints
# the field list on standard output and exits, filter_vep:219-223).
TOOL_FLAGS="-i --input_file -o --output_file --force_overwrite -l --list -h --help --host --user --pass --password --port --version --registry"
# The server filter_vep loads the Ensembl registry from for --ontology, its own defaults
# (filter_vep:67-69: host ensembldb.ensembl.org, port 3306, user anonymous).
TOOL_SERVER="ensembldb.ensembl.org:3306"
# The cache cell's suffix that names a cache root in Ensembl VEP's layout.
LAYOUT_PREFIX="layout:"
# The comparator's own divergence classes, applied without a manifest: the row and text
# classes, the registry pair classes, the field classes of the corrected behaviours but for
# dbnsfp_mitochondrial_record_named_m (it reads a manifest's assembly and accepts nothing
# without one, so a case that names it needs a manifest.json beside its input), the
# transcript row classes of the merged and RefSeq caches, the structural-variant comparator's
# classes and the warnings-mode classes (tools/vep-formats-compare/data/classes.tsv: the
# columns REGISTRY_CLASSES, FIELD_CLASSES, TERM_PAIR_CLASSES, STRUCTURAL_VARIANT_CLASSES,
# PLUGIN_WARNING_CLASSES and WARNING_TEXT_CLASSES, and the class merged_same_stable_id_collapse),
# less the classes whose CORRECTED_IN cell names an Ensembl VEP release: the comparator
# applies none of them under that release or a later one.
BUILTIN_MASKS="regulatory_feature_rows order_only individual_allele_order conservation_short_name_spelling pc_description_text motif_score_reverse_strand_mnv colocated_list_cut_at_cache_cap start_cooccurrence_swap splice_family_swap gxf_reverse_strand_utr_sequence organelle_peptide_paths_translation_table conservation_max_text_order gwas_haplotype_cross_product gwas_interaction_first_snp_allele_whitespace alphamissense_two_digit_transcript_version cadd_structural_allele_match splicevault_sibling_allele_rows splicevault_json_later_rows utrannotator_annotation_pair_order utrannotator_gained_start_window utrannotator_consequence_precedence utrannotator_structural_allele_text utrannotator_lost_start_after_insertion utrannotator_gained_start_remade_in_place utrannotator_frameshift_at_stop_boundary utrannotator_lost_stop_written_before_stop_codon utrannotator_without_reference nearest_exon_jb_item_order loeuf_zero_minimum maxentscan_overlapping_windows alphamissense_sibling_allele_rows revel_sibling_allele_rows revel_other_reference_amino_acid_rows gnomadmt_info_value_second_equals loftee_anc_allele_multiallelic_record loftee_first_alt_donor_filter loftee_gerp_zoom_level_summary loftee_gerp_walk_without_conservation_database loftee_nagnag_site_first_transcript_strand nearest_exon_jb_stretch_priming nearest_exon_jb_shifted_span nearest_exon_jb_id_less_exon_fold nearest_start_site_one_base_above downstream_reverse_strand_rotation spliceai_chopped_multiallelic spliceai_multi_entry_cell plugin_contig_mitochondrial_spelling tabix_allele_case_insensitive_match refseq_edited_model_hgvsc_offset refseq_edited_model_shifted_insertion_bases refseq_edited_model_seam_insertion_hgvsc refseq_edited_model_unminimised_allele_hgvsc refseq_edited_model_hgvs_offset refseq_edited_model_deleted_span_term refseq_edited_model_tail_term refseq_edited_model_equal_span_insertion_term mature_mirna_unshifted_span stop_codon_unshifted_span structural_deletion_stop_lost_without_stop_codon stop_codon_completed_by_polyadenylation stop_codon_overlap_without_stop_codon stop_codon_behind_last_residue_without_stop_codon start_codon_overlap_without_start_codon stop_retained_without_stop_codon_test stop_codon_overlap_genomic_span stop_retained_three_base_stop_codon_test stop_codon_removed_without_stop_codon_test stop_codon_overlap_insertion_range stop_codon_settled_before_unknown_bases frameshift_inside_a_held_stop_codon start_codon_removed_without_five_prime_utr start_codon_kept_without_five_prime_utr cnv_tr_symbolic_loss_named_a_gain splice_term_past_the_record_span refseq_edited_model_deleted_span_codons alternate_codon_past_the_coding_end refseq_edited_model_intron_span_coding_reference refseq_edited_model_inserted_span_codons hgvs_intron_offset_abutting_insertion hgvs_shifted_insertion_clamped_end shift_intron_term_unshifted_span hgvs_shifted_indel_written_span_gate hgvsp_shifted_indel_written_span_gate hgvsp_shifted_indel_written_span_start_test hgvsp_start_codon_test_past_the_shifted_allele hgvs_shifted_deletion_clamped_span hgvs_clipped_reference_against_whole_alternate hgvs_added_bases_left_at_the_span_end hgvs_insertion_left_where_the_shift_stops hgvs_added_bases_written_outside_the_transcript hgvs_shift_sequence_edge_window hgvs_shift_shared_across_alleles hgvsp_stop_codon_protein_comparison hgvsp_stop_named_on_incomplete_coding_sequence hgvsp_frameshift_count_past_a_stop_read_as_a_residue hgvsp_frameshift_with_no_changed_residue hgvsp_frameshift_count_at_the_first_residue allele_number_repeated_or_reference_equal_alt hgvsp_intron_end_insertion hgvsp_intron_start_insertion merged_same_stable_id_collapse use_transcript_ref_cache_cdna refseq_failed_edit_model_hgvsc_offset hgnc_id_copy downstream_utr_frame_past_coding_end downstream_start_past_own_codon plugin_hgvs_shifted_cds_coordinates plugin_value_follows_splice_term_past_the_record_span plugin_value_follows_stop_codon_term_pair total_length_json_range_end plugin_value_on_reference_skipped_record blosum62_hgvs_shift_shared_across_alleles reference_skipped_record breakend_mate_local_read bnd_own_chromosome_rows_lost dosage_sensitivity_regulatory_row_warnings fork_warning_pid_prefix json_mismatch_warning_names_no_record"
# The built-in classes that read beside the two outputs: the skipped-record class, the
# classes that read a record's alleles and the classes that read the reference bases of a row
# under an identifier (the two classes of splice_term_past_the_record_span,
# stop_codon_overlap_insertion_range, plugin_value_follows_stop_codon_term_pair,
# refseq_edited_model_equal_span_insertion_term and
# hgvsp_frameshift_count_at_the_first_residue), the input
# VCF in every format (the comparator reads it for the VCF and JSON formats in any case), the
# breakend classes, hgvsp_intron_end_insertion and hgvsp_intron_start_insertion the candidate's
# JSON cache (vep-formats-compare --json-cache).
# organelle_peptide_paths_translation_table reads its transcripts' codon tables from that cache.
# The two Downstream classes read their transcripts' coding sequences, 3' UTRs and codon tables from it.
# plugin_hgvs_shifted_cds_coordinates reads the transcripts' strands, exons, coding bounds and mapper
# pairs from it as well, and the record's alleles for a later allele's row.
# hgvs_shift_shared_across_alleles decides three columns of a later allele's row from the record's
# alleles and the cache's mapper pairs, coding sequences and codon tables; its other rows are the
# manifest's.
# hgvs_clipped_reference_against_whole_alternate counts a transcript's exon spans from it.
# hgvsp_shifted_indel_written_span_gate counts them for a row whose shift is longer than the
# coding coordinate its notation gives, and leaves such a row a difference without it.
# maxentscan_overlapping_windows reads a transcript's introns from it and sets nothing aside without it.
# refseq_edited_model_inserted_span_codons reads a transcript's coding sequence and codon table from it
# and sets nothing aside without it.
INPUT_READING_MASKS="reference_skipped_record splicevault_sibling_allele_rows alphamissense_sibling_allele_rows revel_sibling_allele_rows revel_other_reference_amino_acid_rows loftee_anc_allele_multiallelic_record loftee_first_alt_donor_filter loftee_nagnag_site_first_transcript_strand spliceai_chopped_multiallelic utrannotator_gained_start_window hgvs_shift_shared_across_alleles plugin_hgvs_shifted_cds_coordinates splice_term_past_the_record_span plugin_value_follows_splice_term_past_the_record_span stop_codon_overlap_insertion_range stop_codon_completed_by_polyadenylation stop_codon_overlap_without_stop_codon stop_codon_behind_last_residue_without_stop_codon start_codon_overlap_without_start_codon plugin_value_follows_stop_codon_term_pair blosum62_hgvs_shift_shared_across_alleles refseq_edited_model_equal_span_insertion_term hgvsp_frameshift_count_at_the_first_residue"
CACHE_READING_MASKS="breakend_mate_local_read bnd_own_chromosome_rows_lost organelle_peptide_paths_translation_table downstream_utr_frame_past_coding_end downstream_start_past_own_codon plugin_hgvs_shifted_cds_coordinates hgvs_shift_shared_across_alleles hgvs_clipped_reference_against_whole_alternate hgvsp_intron_end_insertion hgvsp_intron_start_insertion maxentscan_overlapping_windows refseq_edited_model_inserted_span_codons hgvsp_shifted_indel_written_span_gate"
# hgvsp_stop_named_on_incomplete_coding_sequence reads from the candidate's JSON cache the
# transcripts whose coding sequence is no whole number of codons, for a row whose FLAGS and
# terms do not show the incomplete coding sequence.
CACHE_READING_MASKS="$CACHE_READING_MASKS hgvsp_stop_named_on_incomplete_coding_sequence"
# The built-in class the comparator decides from the FASTA the case gives both engines and
# the candidate's JSON cache (vep-formats-compare --reference-fasta with --json-cache); a
# case that names no FASTA leaves it to the manifest.json beside its input.
REFERENCE_READING_MASKS="hgvs_shift_sequence_edge_window"
# hgvsp_stop_codon_protein_comparison reads a transcript's coding sequence and codon table from
# the candidate's JSON cache for one of its shapes, a pair of more than one residue before an
# inserted stop, which stays a difference without it.
CACHE_READING_MASKS+=" hgvsp_stop_codon_protein_comparison"
# hgvsp_frameshift_count_past_a_stop_read_as_a_residue derives the count of every row it
# takes from the coding sequence and 3' UTR of the row's transcript in the candidate's JSON
# cache, and leaves such a row a difference without it.
CACHE_READING_MASKS+=" hgvsp_frameshift_count_past_a_stop_read_as_a_residue"
# refseq_failed_edit_model_hgvsc_offset reads, from the candidate's JSON cache, the edits of the
# row's transcript that write cDNA positions as a text of another length and its coding
# positions, and leaves such a row a difference without it.
CACHE_READING_MASKS+=" refseq_failed_edit_model_hgvsc_offset"
# hgvsp_frameshift_count_at_the_first_residue derives the count of every row it takes from
# what the candidate's JSON cache holds of the row's transcript: that no start codon stands
# at its coding start (the cds_start_NF attribute, the start phase, the peptide, the
# consensus coding sequence id), its coding sequence, 3' UTR and codon table. It sets
# nothing aside without the cache, and stands in INPUT_READING_MASKS for the length of the
# record's reference bases at the row.
CACHE_READING_MASKS+=" hgvsp_frameshift_count_at_the_first_residue"
# The built-in classes that read the GWAS catalog file the case's flags hand the GWAS plugin
# (vep-formats-compare --gwas-catalog) and set nothing aside without it.
CATALOG_READING_MASKS="gwas_haplotype_cross_product gwas_interaction_first_snp_allele_whitespace"
# hgvs_insertion_left_where_the_shift_stops holds a duplication to its transcript's exon spans
# in the candidate's JSON cache; without the cache it reads Ensembl VEP's duplication by its
# count of bases alone and leaves the candidate's around an intron a difference.
CACHE_READING_MASKS+=" hgvs_insertion_left_where_the_shift_stops"
# stop_codon_settled_before_unknown_bases reads, from the candidate's JSON cache, the model and the
# coding sequence of the row's transcript for an insertion written in an intron that the HGVS 3'
# shift carries into the coding sequence; without it the class takes its rows with a Codons cell
# alone.
CACHE_READING_MASKS+=" stop_codon_settled_before_unknown_bases"
# The built-in class that reads the score files the case's flags hand the SpliceAI plugin
# (vep-formats-compare --spliceai-table, once a file) and sets nothing aside without them.
# It reads the candidate's JSON cache as well, for the gene of a row that writes no SYMBOL.
SCORE_READING_MASKS="spliceai_multi_entry_cell"
CACHE_READING_MASKS="$CACHE_READING_MASKS $SCORE_READING_MASKS"
# The built-in classes a discovery pass asks for only where the case's cell, its manifest or
# --discover-classes names them. regulatory_feature_rows drops the reference's regulatory
# rows whether or not the candidate writes any, so asking for it changes what is compared.
# Under individual_allele_order and the structural-variant classes the comparator compares
# the whole outputs, where it otherwise reads the two as streams and parses the record
# groups that differ alone (the column FRONT_END_WHOLE_CLASSES of its data/classes.tsv),
# and the structural-variant classes refuse an input that is no VCF.
DISCOVERY_OPT_IN_MASKS="regulatory_feature_rows individual_allele_order reference_skipped_record breakend_mate_local_read bnd_own_chromosome_rows_lost"
# The classes of a coding end at which Ensembl VEP reads a stop codon or a start codon that
# the transcript does not show. Each derives the terms vep-rs owes on a row from what the
# candidate's JSON cache holds of the row's transcript and sets nothing aside without the
# cache. The four that read a sequence allele's row stand in INPUT_READING_MASKS too, for
# the length of the record's reference bases at the row.
# structural_deletion_stop_lost_without_stop_codon reads whether a stop codon closes the
# coding sequence and whether a start codon opens it (the coding sequence, the codon table,
# the peptide, the start phase, the cds_end_NF and cds_start_NF attributes, the consensus
# coding sequence id) and the genomic bases of both coding ends (the strand, the mapper
# pairs); it sets nothing aside without the cache.
CACHE_READING_MASKS+=" structural_deletion_stop_lost_without_stop_codon"
# stop_codon_completed_by_polyadenylation reads those and, of a transcript whose stop codon
# adenines complete, the coding sequence, the peptide and the codon table whole, from which
# it derives the row's terms, residues, codons and protein notation; it sets nothing aside
# without the cache.
CACHE_READING_MASKS+=" stop_codon_completed_by_polyadenylation"
# stop_codon_overlap_without_stop_codon reads what the structural class reads and the
# transcript's coding bounds, mapper pairs and introns, by which it places the row's allele;
# it sets nothing aside without the cache.
CACHE_READING_MASKS+=" stop_codon_overlap_without_stop_codon"
# stop_codon_behind_last_residue_without_stop_codon reads what the overlap class reads and
# sets nothing aside without the cache.
CACHE_READING_MASKS+=" stop_codon_behind_last_residue_without_stop_codon"
# start_codon_overlap_without_start_codon reads what the overlap class reads and sets
# nothing aside without the cache.
CACHE_READING_MASKS+=" start_codon_overlap_without_start_codon"
# The classes of a stop codon or a start codon that stands. Beside the cache each takes a
# pair only on a transcript whose coding sequence shows that codon, read as the five classes
# above read it; given no cache the comparator applies each by the two rows alone, so a
# case that names one hands the cache.
CACHE_READING_MASKS+=" start_cooccurrence_swap stop_codon_unshifted_span stop_retained_without_stop_codon_test stop_codon_overlap_genomic_span stop_retained_three_base_stop_codon_test stop_codon_removed_without_stop_codon_test stop_codon_overlap_insertion_range"
# The two classes of a plugin value that follows a term pair take a value only beside a pair
# a row class takes, which beside the cache is read on the coding ends of the row's
# transcript as that class reads it; given no cache the comparator reads the two rows alone,
# so a case that names one hands the cache.
CACHE_READING_MASKS+=" plugin_value_follows_splice_term_past_the_record_span plugin_value_follows_stop_codon_term_pair"
# stop_codon_completed_by_polyadenylation reads the reference bases past a transcript's end
# for one shape, a deletion that runs past that end, from the FASTA the case gives both
# engines (vep-formats-compare --reference-fasta); without a FASTA that shape stays a
# difference and every other row of the class is decided from the cache.
REFERENCE_READING_MASKS+=" stop_codon_completed_by_polyadenylation"
# The manifest-only classes: the comparator's MANIFEST_CLASSES (the column of that name in
# its data/classes.tsv) less the classes of BUILTIN_MASKS. No rule over the two outputs
# decides one, so a case that names one needs a manifest.json beside its input.
MANIFEST_ONLY_MASKS="transcript_annotator_sibling_allele_rows transcript_annotator_insertion_window frequency_filter_allele_key_order minimal_symbolic_alleles shift_3prime_shared_predicates shift_start_codon_test"
# Flags the runner sets from the table's own columns and its options, with the aliases
# Ensembl VEP's Config.pm and the candidate accept for them; a flags cell naming one would
# give an engine the flag twice, or move the input, the server or the assembly behind the
# report's back. A cell names one in any spelling an engine reads as it: a single dash or a
# plus sign, any letter case, another name of the same option (-u, --full_cache_dir) and an
# abbreviation, each resolved against the two option lists below.
# --dir is not among them: the runner sets none, a command's
# own cache words outrank it (AnnotationSourceAdaptor.pm:127 at release/116.3 reads
# `dir_cache || dir`), and preflight holds its directory to a read-only placeholder mount
# without a vep.ini, the file Config.pm:704-705 reads as options from under --dir.
RUNNER_FLAGS="-i --input_file --input_data --id -o --output_file --offline --cache --database --db_version --host --port --user --password --pass --dir_cache --cache_version --json_cache -s --species -a --assembly --grch37 --fasta --fa --tab --vcf --json --output_format --force_overwrite --force --no_stats -q --quiet --fork --threads --warning_file --skipped_variants_file --config"
# The options each engine declares, one word per option, its names joined by `|`.
# PERL_OPTIONS is Ensembl VEP's `@VEP_PARAMS` (Config.pm:73-276 at release/116.3) in the
# source's order, which VEP_PARAMS of crates/vep-cli/src/argv.rs carries; CANDIDATE_OPTIONS
# is the candidate's parser, sorted: the flag rows of the reference pages under docs/cli/
# that are generated from it, and the flags its error table catches before parsing.
# test_feature_oracle.sh holds each list to that source.
PERL_OPTIONS="help config input_file|i input_data|id format output_format delimiter no_check_variants_order species|s registry host port user|u password|pass db_version assembly|a grch37 genomes refseq merged all_refseq gencode_basic gencode_primary is_multispecies transcript_filter exclude_predicted minimal most_severe summary pick pick_allele per_gene pick_allele_gene flag_pick flag_pick_allele flag_pick_allele_gene flag_gencode_primary pick_order buffer_size failed gp chr check_ref lookup_ref check_existing check_svs no_check_alleles exclude_null_alleles check_frequency af af_1kg af_gnomade af_gnomad af_gnomadg old_maf max_af pubmed freq_filter freq_freq freq_gt_lt freq_pop filter_common allow_non_variant process_ref_homs individual individual_zyg phased fork dont_skip nearest distance clin_sig_allele clinvar_somatic_classification overlaps max_sv_size remove_hgvsp_version verbose|v quiet|q no_progress everything|e output_file|o compress_output no_headers mask_header_cache_path stats_file|sf stats_text stats_html no_stats warning_file skipped_variants_file force_overwrite|force terms|t coding_only canonical mane mane_select tsl appris ccds xref_refseq uniprot protein protein_version biotype hgnc symbol transcript_version gene_version gene_phenotype mirna spdi ga4gh_vrs hgvs hgvsg hgvsg_use_accession hgvsp_use_prediction shift_hgvs ambiguous_hgvs sift polyphen humdiv condel variant_class regulatory cell_type regulatory_gff extended_promoters convert no_intergenic vcf solr json tab vcf_info_field keep_csq keep_ann lrg fields domains numbers total_length allele_number show_ref_allele no_escape ambiguity var_synonyms shift_3prime shift_genomic shift_length uploaded_allele database cache cache_version show_cache_info dir dir_cache full_cache_dir dir_plugins offline fasta|fa fasta_dir no_fasta sereal synonyms bam use_transcript_ref use_given_ref custom tmpdir gff gtf bigwig phyloP phastCons ucsc_assembly ucsc_data_root custom_multi_allelic custom_suppress_filter plugin safe debug"
CANDIDATE_OPTIONS="af af_1kg af_exac af_gnomad af_gnomade af_gnomadg all_refseq allele_number allow_non_variant ambiguity ambiguous_hgvs appris assembly|a bam bigwig biotype buffer_size cache cache_version canonical ccds cell_type check_existing check_frequency check_ref check_svs chr clin_sig_allele clinvar_somatic_classification coding_only compress_output compression_threads condel config convert custom custom_multi_allelic custom_suppress_filter database db_version debug decompression_threads delimiter dir dir_cache dir_plugins distance domains dont_skip everything|e exclude_null_alleles exclude_predicted extended_promoters failed fasta_dir fasta|fa fields filter_common flag_gencode_primary flag_pick flag_pick_allele flag_pick_allele_gene force_overwrite|force fork|threads format freq_filter freq_freq freq_gt_lt freq_pop ga4gh_vrs gencode_basic gencode_primary gene_phenotype gene_version genomes gff gp grch37 gtf hgnc hgvs hgvsg hgvsg_use_accession hgvsp_use_prediction host humdiv individual individual_zyg input_data|id input_file|i is_multispecies json json_cache|full_cache_dir keep_ann keep_csq lookup_ref lrg mane mane_select mask_header_cache_path max_af max_sv_size merged minimal mirna most_severe nearest no_check_alleles no_check_variants_order no_escape no_fasta no_headers no_intergenic no_progress no_stats numbers offline old_maf output_buffer_size output_file|o output_format overlaps parquet_partition parse_threads password|pass per_gene phased phastCons phyloP pick pick_allele pick_allele_gene pick_order pipeline_depth pipeline_mode plugin polyphen port prewarm_ahead process_ref_homs protein protein_version pubmed quiet|q reader_core refseq registry regulatory regulatory_gff remove_hgvsp_version render_task_records resolved_config safe sereal shift_3prime shift_genomic shift_hgvs shift_length show_cache_info show_ref_allele sift skipped_variants_file solr spdi species|s stats_file|sf stats_html stats_text summary symbol synonyms tab terms|t tmpdir total_length transcript_filter transcript_version tsl ucsc_assembly ucsc_data_root uniprot uploaded_allele use_given_ref use_transcript_ref user|u var_synonyms variant_class vcf vcf_info_field verbose|v warning_file writer_strategy xref_refseq"
# The word of RUNNER_FLAGS that the option name in OPTION_NAME reaches in either list, or
# nothing. An engine reads a name as the option that has it among its names, else as the
# one option whose names begin with it (Getopt::Long's auto_abbrev, which Config.pm:62
# leaves on; resolve_spelling of argv.rs reads the candidate's options, then Ensembl
# VEP's), and a name that several options begin with reaches none.
# The spellings of a path under the corpus that the run match reads as one, a line per
# spelling as <written><tab><read as>: the two words a manifest's runs write for the corpus's
# plugin_data and generated directories, and the corpus placeholder before a path, which a
# manifest leaves out where it writes the path from the corpus directory.
RUN_SPELLINGS=$(printf '%s\t%s\n' "$PLUGINS_PLACEHOLDER" "$CORPUS_PLACEHOLDER/plugin_data" \
	"$GENERATED_PLACEHOLDER" "$CORPUS_PLACEHOLDER/generated" "$CORPUS_PLACEHOLDER/" "")
# The runs of a manifest.json against a case, as tab-separated lines. `fit`: the name of the
# one run, the manifest's own (`default`) or a flag set, whose input has the file name of
# argument 2 and whose flags are the words after it in any order, or nothing where no run or
# more than one fits. `over`: the number of runs whose input has that file name, and their
# names. `run`: every run's name and input. Nothing where the manifest cannot be read. The
# manifest's side is read without --grch37 and without a --fasta word and its file; both
# sides without --dir_plugins and its directory, each word through RUN_SPELLINGS and each
# option under the first name of the option an engine reads the word as (PERL_OPTIONS, then
# CANDIDATE_OPTIONS). The names are vep-formats-compare's OWN_RUN and CORPUS_INPUT.
CASE_RUN_PROGRAM='
import json, os, re, sys
try:
    with open(sys.argv[1], encoding="utf-8") as handle:
        manifest = json.load(handle)
except (OSError, ValueError):
    sys.exit(0)
if not isinstance(manifest, dict):
    sys.exit(0)
spellings = [line.split("\t") for line in os.environ.get("RUN_SPELLINGS", "").split("\n") if "\t" in line]
lists = [[names.split("|") for names in os.environ.get(key, "").split()] for key in ("PERL_OPTIONS", "CANDIDATE_OPTIONS")]
def option(name):
    for options in lists:
        named = [names for names in options if name in [n.lower() for n in names]]
        begun = [names for names in options if any(n.lower().startswith(name) for n in names)]
        if named or len(begun) == 1:
            return (named or begun)[0][0]
    return None
def spelled(word):
    word = str(word)
    for written, read in spellings:
        word = word.replace(written, read)
    found = re.fullmatch(r"(?:--|-|\+)([A-Za-z][A-Za-z0-9_]*)(=.*)?", word, re.S)
    name = option(found.group(1).lower()) if found else None
    return "--" + name + (found.group(2) or "") if name else word
def words(flags, of_manifest):
    kept, value = [], False
    for word in map(spelled, flags if isinstance(flags, list) else []):
        name = word.split("=", 1)[0]
        if value:
            value = False
        elif name == "--dir_plugins" or (of_manifest and name == "--fasta"):
            value = "=" not in word
        elif not (of_manifest and word == "--grch37"):
            kept.append(word)
    return sorted(kept)
sets = manifest.get("flag_sets")
runs = [("default", manifest.get("flags"), "variants.vcf")]
runs += [(s.get("name"), s.get("flags"), s.get("input", "variants.vcf")) for s in sets if isinstance(s, dict)] if isinstance(sets, list) else []
case = words(sys.argv[3:], False)
fits = [name for name, flags, run_input in runs if run_input == sys.argv[2] and words(flags, True) == case]
over = [str(name) for name, _, run_input in runs if run_input == sys.argv[2]]
print("fit\t" + (fits[0] if len(fits) == 1 and isinstance(fits[0], str) else ""))
print("over\t" + str(len(over)) + "\t" + ", ".join(over))
for name, _, run_input in runs:
    if isinstance(name, str):
        print("run\t" + name + "\t" + str(run_input))
'
# The classes a manifest.json documents, as one line of words: the expected_divergence and
# the second_expected_divergence of every item of its own run and of each of its flag sets;
# nothing where the manifest cannot be read.
MANIFEST_CLASSES_PROGRAM='
import json, re, sys
try:
    with open(sys.argv[1], encoding="utf-8") as handle:
        manifest = json.load(handle)
except (OSError, ValueError):
    sys.exit(0)
runs = [manifest] if isinstance(manifest, dict) else []
sets = manifest.get("flag_sets") if runs else None
runs += [s for s in sets if isinstance(s, dict)] if isinstance(sets, list) else []
names = set()
for run in runs:
    for section in ("divergences", "vep_rs_only_tuples", "field_divergences"):
        items = run.get(section)
        for item in items if isinstance(items, list) else []:
            if isinstance(item, dict):
                values = (item.get(key) for key in ("expected_divergence", "second_expected_divergence"))
                names.update(value for value in values if isinstance(value, str))
print(" ".join(sorted(n for n in names if re.fullmatch("[a-z0-9_]+", n))))
'
RUNNER_WORD_PROGRAM='
function option_of(name, list,    options, names, n, m, i, j, hit, hits) {
	n = split(list, options, " ")
	hit = hits = 0
	for (i = 1; i <= n; i++) {
		m = split(tolower(options[i]), names, "|")
		for (j = 1; j <= m; j++) {
			if (names[j] == name) return i
			if (index(names[j], name) == 1 && hit != i) {
				hit = i
				hits++
			}
		}
	}
	return hits == 1 ? hit : 0
}
function runner_word(name, list,    options, names, m, i, j, pass) {
	i = option_of(name, list)
	if (!i) return ""
	split(list, options, " ")
	m = split(options[i], names, "|")
	for (pass = 1; pass <= 2; pass++) {
		for (j = 1; j <= m; j++) {
			if (pass == 1 && index(tolower(names[j]), name) != 1) continue
			if (index(runner, " --" names[j] " ")) return "--" names[j]
			if (index(runner, " -" names[j] " ")) return "-" names[j]
		}
	}
	return ""
}
BEGIN {
	name = ENVIRON["OPTION_NAME"]
	word = runner_word(name, perl)
	if (word == "") word = runner_word(name, candidate)
	printf "%s", word
}'

# A run's table describes that run: the previous one goes before anything is checked, so
# a run that fails preflight or is interrupted leaves no table behind.
[[ "$DRY_RUN" == "true" ]] || rm -f "$OUT_DIR/status.tsv" "$OUT_DIR/run.json" "$OUT_DIR/discovery.tsv"

# Header: every column present once; the column order is read from it.
HEADER=""
while IFS= read -r line; do
	[[ -n "$line" && "$line" != \#* ]] || continue
	HEADER="$line"
	break
done <"$CASES"
[[ -n "$HEADER" ]] || die "$CASES: no header line"
N_COLUMNS=$(printf '%s\n' "$HEADER" | awk -F'\t' '{print NF}')
# COL_<name>=<1-based index>
for name in $TABLE_COLUMNS; do
	idx=$(awk -F'\t' -v n="$name" '{for (i = 1; i <= NF; i++) if ($i == n) {print i; exit}}' <<<"$HEADER")
	[[ -n "$idx" ]] || die "$CASES: header lacks the column '$name' (needs: $TABLE_COLUMNS)"
	eval "COL_$name=$idx"
done
[[ "$N_COLUMNS" -eq "$N_TABLE_COLUMNS" ]] || die "$CASES: header has $N_COLUMNS columns; the $N_TABLE_COLUMNS are: $TABLE_COLUMNS"

# cell <line> <column name>: the cell's text.
cell() {
	local idx
	eval "idx=\$COL_$2"
	printf '%s\n' "$1" | awk -F'\t' -v c="$idx" '{print $c}'
}

# The on-demand list (the header's ON-DEMAND LIST): ON_DEMAND_NAMES holds the names of its cases.
ON_DEMAND="$(dirname "$CASES")/oracle_cases_on_demand.tsv"
ON_DEMAND_NAMES=""
if [[ -f "$ON_DEMAND" ]]; then
	line_no=0
	seen_header="false"
	while IFS= read -r line || [[ -n "$line" ]]; do
		line_no=$((line_no + 1))
		[[ -n "$line" && "$line" != \#* ]] || continue
		if [[ "$seen_header" != "true" ]]; then
			[[ "$line" == $'case\treason' ]] || die "$ON_DEMAND:$line_no: the header is not case, reason (tab-separated)"
			seen_header="true"
			continue
		fi
		awk -F'\t' 'NF != 2 {exit 1} {for (i = 1; i <= NF; i++) if ($i ~ /^[[:space:]]*$/) exit 1}' <<<"$line" ||
			die "$ON_DEMAND:$line_no: a row is case, reason: two tab-separated cells, neither empty"
		cid="${line%%$'\t'*}"
		# The names are read below as the words of one list: a cell holding a blank would name two
		# cases there, and one holding a glob the files of the working directory.
		[[ "$cid" =~ $PIN_RE ]] || die "$ON_DEMAND:$line_no: '$cid' is not a case's name ([A-Za-z0-9._-], not starting with . or -)"
		ON_DEMAND_NAMES="$ON_DEMAND_NAMES $cid"
	done <"$ON_DEMAND"
	[[ "$seen_header" == "true" ]] || die "$ON_DEMAND: no header line"
fi

# Row selection: the tier, then --only; without --only, the tier's cases less the on-demand ones.
ONLY_NAMES=""
if [[ -n "$ONLY_SPEC" ]]; then
	ONLY_NAMES=$(printf '%s' "$ONLY_SPEC" | tr ',' ' ')
fi
ON_DEMAND_LEFT_OUT=""
SELECTED=""
ALL_NAMES=""
TIER_NAMES=""
line_no=0
seen_header="false"
while IFS= read -r line; do
	line_no=$((line_no + 1))
	[[ -n "$line" && "$line" != \#* ]] || continue
	if [[ "$seen_header" != "true" ]]; then
		seen_header="true"
		continue
	fi
	nf=$(printf '%s\n' "$line" | awk -F'\t' '{print NF}')
	[[ "$nf" -eq "$N_TABLE_COLUMNS" ]] || die "$CASES:$line_no: $nf columns, expected $N_TABLE_COLUMNS (an empty cell is written -)"
	printf '%s\n' "$line" | awk -F'\t' '{for (i = 1; i <= NF; i++) if ($i == "") exit 1}' ||
		die "$CASES:$line_no: an empty cell; write - for no value"
	cid=$(cell "$line" case)
	[[ "$cid" =~ $PIN_RE ]] || die "$CASES:$line_no: case '$cid' is not a directory-safe name ([A-Za-z0-9._-], not starting with . or -)"
	has_word "$ALL_NAMES" "$cid" && die "$CASES:$line_no: case '$cid' appears twice"
	ALL_NAMES="$ALL_NAMES $cid"
	tier=$(cell "$line" tier)
	case "$tier" in
	fixture | population) ;;
	*) die "$CASES:$line_no: tier '$tier' is not fixture or population" ;;
	esac
	# A caller that sets a fixture case aside hands the runner the rest of the tier as an
	# --only list, which would name an on-demand fixture case.
	if [[ "$tier" != "population" ]] && has_word "$ON_DEMAND_NAMES" "$cid"; then
		die "$ON_DEMAND names '$cid', a $tier case; the list holds population cases alone"
	fi
	[[ "$tier" == "$TIER" ]] || continue
	TIER_NAMES="$TIER_NAMES $cid"
	if [[ -n "$ONLY_NAMES" ]]; then
		has_word "$ONLY_NAMES" "$cid" || continue
	elif has_word "$ON_DEMAND_NAMES" "$cid"; then
		ON_DEMAND_LEFT_OUT="$ON_DEMAND_LEFT_OUT $cid"
		continue
	fi
	SELECTED="${SELECTED}${line}
"
done <"$CASES"
for cid in $ONLY_NAMES; do
	has_word "$ALL_NAMES" "$cid" || die "--only names '$cid', which is not a case of $CASES"
	has_word "$TIER_NAMES" "$cid" || die "--only names '$cid', which is not a case of tier $TIER"
done
for cid in $ON_DEMAND_NAMES; do
	has_word "$ALL_NAMES" "$cid" || die "$ON_DEMAND names '$cid', which is not a case of $CASES"
done
N_SELECTED=$(printf '%s' "$SELECTED" | grep -c . 2>/dev/null) || N_SELECTED=0
[[ "$N_SELECTED" -gt 0 ]] || die "no case of tier $TIER selected from $CASES${ONLY_SPEC:+ by --only $ONLY_SPEC}${ON_DEMAND_LEFT_OUT:+; left out, on demand (--only names a case to run it):$ON_DEMAND_LEFT_OUT}"
[[ -z "$ON_DEMAND_LEFT_OUT" ]] || log "left out, on demand (--only names a case to run it):$ON_DEMAND_LEFT_OUT"

# The named-runs table (the header's NAMED-RUNS TABLE): CASE_RUN_ROWS holds its rows as
# they are written, one `case<tab>run<tab>reason` line each.
RUN_NAME_RE='^[A-Za-z0-9][A-Za-z0-9._-]*$'
if [[ -z "$CASE_RUNS" && -f "$(dirname "$CASES")/oracle_case_runs.tsv" ]]; then
	CASE_RUNS="$(dirname "$CASES")/oracle_case_runs.tsv"
fi
CASE_RUN_ROWS=""
if [[ -n "$CASE_RUNS" ]]; then
	[[ -f "$CASE_RUNS" ]] || die "named-runs table not found: $CASE_RUNS"
	line_no=0
	seen_header="false"
	seen_run_cases=""
	while IFS= read -r line || [[ -n "$line" ]]; do
		line_no=$((line_no + 1))
		[[ -n "$line" && "$line" != \#* ]] || continue
		if [[ "$seen_header" != "true" ]]; then
			[[ "$line" == $'case\trun\treason' ]] || die "$CASE_RUNS:$line_no: the header is not case, run, reason (tab-separated)"
			seen_header="true"
			continue
		fi
		awk -F'\t' 'NF != 3 {exit 1} {for (i = 1; i <= NF; i++) if ($i ~ /^[[:space:]]*$/) exit 1}' <<<"$line" ||
			die "$CASE_RUNS:$line_no: a row is case, run, reason: three tab-separated cells, none empty"
		cid="${line%%$'\t'*}"
		run_name="${line#*$'\t'}"
		run_name="${run_name%%$'\t'*}"
		# has_word pads its word with blanks, so a cell that holds one can match where no case
		# has that name: the cell is held to the form of a case's name.
		{ [[ "$cid" =~ $PIN_RE ]] && has_word "$ALL_NAMES" "$cid"; } || die "$CASE_RUNS:$line_no: '$cid' is not a case of $CASES"
		has_word "$seen_run_cases" "$cid" && die "$CASE_RUNS:$line_no: case '$cid' has a second row; a case compares under one run"
		[[ "$run_name" =~ $RUN_NAME_RE ]] || die "$CASE_RUNS:$line_no: run '$run_name' is not a run's name (default, a flag set's name) or pooled"
		seen_run_cases="$seen_run_cases $cid"
		CASE_RUN_ROWS="${CASE_RUN_ROWS}${line}
"
	done <"$CASE_RUNS"
	[[ "$seen_header" == "true" ]] || die "$CASE_RUNS: no header line"
fi

# One case: fields and resolved paths

# resolve_path <cell>: absolute as given, else under --data-dir, else under the repository.
resolve_path() {
	local p="$1"
	case "$p" in
	/*) printf '%s' "$p" ;;
	*)
		if [[ -e "$DATA_DIR/$p" ]]; then
			printf '%s/%s' "$DATA_DIR" "$p"
		else
			printf '%s/%s' "$REPO_DIR" "$p"
		fi
		;;
	esac
}

# option_name <word>: the option named by a vep argument word, in its lower-case `--name`
# spelling with a value suffix removed, or nothing for a value word. vep reads its options
# with Getopt::Long's defaults (vep:27, Config.pm:62 at release/116.3), which take a single
# dash or a plus sign, any unambiguous abbreviation and any letter case as the full option.
# An abbreviation is returned as written: preflight resolves it against the option lists.
option_name() {
	local name="${1%%=*}"
	case "$name" in
	--[A-Za-z]*) name=$(lower "$name") ;;
	-[A-Za-z]* | +[A-Za-z]*) name=$(lower "--${name#[-+]}") ;;
	*) return 0 ;;
	esac
	printf '%s' "$name"
}
# has_option <--name>: whether a flag word of the loaded case names the option <--name> in
# a spelling vep reads as it, that is, whether the word's option_name is a prefix of <--name>.
has_option() {
	local word name
	for word in ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}; do
		name=$(option_name "$word")
		[[ -n "$name" && "$1" == "$name"* ]] && return 0
	done
	return 1
}
# gwas_catalog_file: the file the loaded case's flags hand the GWAS plugin as the candidate
# reads it (`--plugin GWAS,file=<path>`, the value its own word or joined by `=`), or nothing
# where no flag word names one.
gwas_catalog_file() {
	local word name value parameter plugin_value="false"
	local parameters=()
	for word in ${C_CAND_FLAG_WORDS[@]+"${C_CAND_FLAG_WORDS[@]}"}; do
		name=$(option_name "$word")
		value=""
		if [[ -z "$name" ]]; then
			[[ "$plugin_value" != "true" ]] || value="$word"
			plugin_value="false"
		elif [[ "--plugin" == "$name"* && "$word" == *=* ]]; then
			value="${word#*=}"
			plugin_value="false"
		elif [[ "--plugin" == "$name"* ]]; then
			plugin_value="true"
		else
			plugin_value="false"
		fi
		[[ "$value" == GWAS,* ]] || continue
		IFS=',' read -r -a parameters <<<"$value"
		for parameter in "${parameters[@]}"; do
			if [[ "$parameter" == file=* ]]; then
				printf '%s' "${parameter#file=}"
				return 0
			fi
		done
	done
	return 0
}
# tss_first_run: whether the loaded case's first Ensembl VEP run writes
# transcript_gene_tss.txt: an offline --nearest case whose cache slice lacks the file.
tss_first_run() {
	[[ "$C_PERL_MODE" == "offline" && ! -e "$C_PERL_SLICE/transcript_gene_tss.txt" ]] && has_option --nearest
}

# load_case <line>: the C_* globals of one row.
load_case() {
	local line="$1" cache_spec species_dir w engine rewritten head corpus sources
	local name value plugin_value parameter parameters=() reads_scores="false"
	local registry="" asked=""
	local run_row run_name run_facts="" run_fit="" run_over="" run_input="" input_name words
	C_ID=$(cell "$line" case)
	C_TIER=$(cell "$line" tier)
	C_INPUT=$(cell "$line" input)
	C_FASTA=$(cell "$line" fasta)
	C_FLAGS=$(cell "$line" flags)
	C_FORMATS=$(cell "$line" formats)
	C_PERL_IMAGE=$(cell "$line" perl_image)
	C_PERL_MODE=$(cell "$line" perl_mode)
	C_REFERENCE_RELEASE=$(cell "$line" reference_release)
	C_MASKS=$(cell "$line" masks)
	# The cell's classes without their expected counts: the comparator's --masks value
	# (write_report reads the counts from the cell itself).
	C_MASK_CLASSES="$C_MASKS"
	[[ "$C_MASKS" == "-" ]] || C_MASK_CLASSES=$(printf '%s' "$C_MASKS" | sed 's/=[^,]*//g')
	C_FLOORS=$(cell "$line" floors)
	C_FIELD_PAIRS=$(cell "$line" field_pairs)
	C_SPECIES=$(cell "$line" species)
	[[ "$C_SPECIES" != "-" ]] || C_SPECIES="$DEFAULT_SPECIES"
	C_UNORDERED=$(cell "$line" unordered)
	C_REFERENCE=$(cell "$line" reference)
	C_REFUSAL=$(cell "$line" refusal)
	C_TOOL_CELL=$(cell "$line" tool)
	# The tool cell: <tool>:<args>, the args read as shell words so a filter string keeps its
	# spaces; C_TOOL_WORDS_ERROR carries the reason the words could not be read, for preflight.
	C_TOOL=""
	C_TOOL_WORDS=()
	C_TOOL_WORDS_ERROR=""
	C_TOOL_COUNT="false"
	if [[ "$C_TOOL_CELL" != "-" ]]; then
		C_TOOL="${C_TOOL_CELL%%:*}"
		# A cell carries no newline, so one word per line is a faithful transport; the closing
		# `end` line keeps a trailing empty word, which the substitution would otherwise strip.
		words=$("$PYTHON" -c 'import shlex, sys
try:
    sys.stdout.write("".join(w + "\n" for w in shlex.split(sys.argv[1])) + "end")
except ValueError as exc:
    sys.exit(f"{exc}")' "${C_TOOL_CELL#*:}" 2>&1) || C_TOOL_WORDS_ERROR="$words"
		if [[ -z "$C_TOOL_WORDS_ERROR" ]]; then
			while IFS= read -r w; do C_TOOL_WORDS+=("$w"); done <<<"$words"
			unset "C_TOOL_WORDS[$((${#C_TOOL_WORDS[@]} - 1))]"
		fi
		for w in ${C_TOOL_WORDS[@]+"${C_TOOL_WORDS[@]}"}; do
			[[ "$w" != "--count" && "$w" != "-c" ]] || C_TOOL_COUNT="true"
		done
	fi
	cache_spec=$(cell "$line" cache)
	[[ "$C_FLAGS" != "-" ]] || C_FLAGS=""
	# The flags as words: `read -a` splits on blanks without expanding a glob, so a `*` in
	# a flag value reaches the engines as written.
	C_FLAG_WORDS=()
	[[ -z "$C_FLAGS" ]] || read -r -a C_FLAG_WORDS <<<"$C_FLAGS"

	C_INPUT_PATH=$(resolve_path "$C_INPUT")
	C_CORPUS_DIR=$(dirname "$C_INPUT_PATH")
	C_MANIFEST="$C_CORPUS_DIR/manifest.json"
	[[ -f "$C_MANIFEST" ]] || C_MANIFEST=""
	# The classes the comparator is asked for, comma-separated, or -: the cell's own, or in a
	# discovery pass the registry beside them (the header's DISCOVERY PASS section), each class
	# once, the registry's first. A tool row and a refusal row are asked for their cell alone.
	C_MANIFEST_CLASSES=""
	C_REQUEST_CLASSES="$C_MASK_CLASSES"
	if [[ "$DISCOVER" == "true" ]]; then
		if [[ "$C_TOOL_CELL" == "-" && "$C_REFUSAL" == "-" ]]; then
			[[ -z "$C_MANIFEST" ]] || C_MANIFEST_CLASSES=$("$PYTHON" -c "$MANIFEST_CLASSES_PROGRAM" "$C_MANIFEST") || C_MANIFEST_CLASSES=""
			if [[ -n "$DISCOVER_CLASSES" ]]; then
				for w in $(printf '%s' "$DISCOVER_CLASSES" | tr ',' ' '); do
					if has_word "$BUILTIN_MASKS" "$w" || has_word "$C_MANIFEST_CLASSES" "$w"; then registry="$registry $w"; fi
				done
			else
				for w in $BUILTIN_MASKS; do
					has_word "$DISCOVERY_OPT_IN_MASKS" "$w" || registry="$registry $w"
				done
				registry="$registry $C_MANIFEST_CLASSES"
			fi
		fi
		[[ "$C_MASKS" == "-" ]] || registry="$registry $(printf '%s' "$C_MASK_CLASSES" | tr ',' ' ')"
		for w in $registry; do
			has_word "$asked" "$w" || asked="$asked $w"
		done
		C_REQUEST_CLASSES="-"
		[[ -z "$asked" ]] || C_REQUEST_CLASSES=$(printf '%s' "${asked# }" | tr ' ' ',')
	fi
	# The flag words as each engine reads them: the two placeholders replaced by the
	# container mounts for Ensembl VEP and by the host directories for the candidate. Not
	# `${w//"$placeholder"/"$dir"}`: bash 3.2 keeps the quotes of a quoted replacement and
	# bash 5.2 reads `&` in an unquoted one as the matched text, so each occurrence is cut
	# out by prefix and suffix removal, literal on both.
	C_PERL_FLAG_WORDS=()
	C_CAND_FLAG_WORDS=()
	for w in ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}; do
		for engine in perl candidate; do
			if [[ "$engine" == perl ]]; then
				corpus="$CORPUS_MOUNT" sources="$SOURCES_MOUNT"
			else
				corpus="$C_CORPUS_DIR" sources="$SOURCES_DIR"
			fi
			rewritten="$w"
			if [[ -n "$C_MANIFEST" ]]; then
				head=""
				while [[ "$rewritten" == *"$PLUGINS_PLACEHOLDER"* ]]; do
					head="$head${rewritten%%"$PLUGINS_PLACEHOLDER"*}$CORPUS_PLACEHOLDER/plugin_data"
					rewritten="${rewritten#*"$PLUGINS_PLACEHOLDER"}"
				done
				rewritten="$head$rewritten"
				head=""
				while [[ "$rewritten" == *"$GENERATED_PLACEHOLDER"* ]]; do
					head="$head${rewritten%%"$GENERATED_PLACEHOLDER"*}$CORPUS_PLACEHOLDER/generated"
					rewritten="${rewritten#*"$GENERATED_PLACEHOLDER"}"
				done
				rewritten="$head$rewritten"
			fi
			head=""
			while [[ "$rewritten" == *"$CORPUS_PLACEHOLDER"* ]]; do
				head="$head${rewritten%%"$CORPUS_PLACEHOLDER"*}$corpus"
				rewritten="${rewritten#*"$CORPUS_PLACEHOLDER"}"
			done
			rewritten="$head$rewritten"
			head=""
			while [[ "$rewritten" == *"$SOURCES_PLACEHOLDER"* ]]; do
				head="$head${rewritten%%"$SOURCES_PLACEHOLDER"*}$sources"
				rewritten="${rewritten#*"$SOURCES_PLACEHOLDER"}"
			done
			if [[ "$engine" == perl ]]; then
				C_PERL_FLAG_WORDS+=("$head$rewritten")
			else
				C_CAND_FLAG_WORDS+=("$head$rewritten")
			fi
		done
	done
	# The score files the flags hand the SpliceAI plugin as the candidate reads them, for a
	# case asked for a score-reading class: the values after snv= and indel= of every
	# `--plugin SpliceAI,...` value, in the cell's order, --plugin read as gwas_catalog_file
	# reads it.
	C_SPLICEAI_TABLES=()
	IFS=',' read -r -a parameters <<<"$C_REQUEST_CLASSES"
	for parameter in "${parameters[@]}"; do
		has_word "$SCORE_READING_MASKS" "$parameter" && reads_scores="true"
	done
	plugin_value="false"
	for w in ${C_CAND_FLAG_WORDS[@]+"${C_CAND_FLAG_WORDS[@]}"}; do
		[[ "$reads_scores" == "true" ]] || break
		name=$(option_name "$w")
		value=""
		if [[ -z "$name" ]]; then
			[[ "$plugin_value" != "true" ]] || value="$w"
			plugin_value="false"
		elif [[ "--plugin" == "$name"* && "$w" == *=* ]]; then
			value="${w#*=}"
			plugin_value="false"
		elif [[ "--plugin" == "$name"* ]]; then
			plugin_value="true"
		else
			plugin_value="false"
		fi
		[[ "$value" == SpliceAI,* ]] || continue
		IFS=',' read -r -a parameters <<<"$value"
		for parameter in "${parameters[@]}"; do
			case "$parameter" in
			snv=?* | indel=?*) C_SPLICEAI_TABLES+=("${parameter#*=}") ;;
			esac
		done
	done
	# The fasta cell's marker (the header's fasta column): C_FASTA_STAGED is `writable` for
	# unindexed:<path>, `readonly` for unindexed-readonly:<path>, else nothing; C_FASTA_PATH is
	# the path after the marker.
	C_FASTA_PATH=""
	C_FASTA_STAGED=""
	case "$C_FASTA" in
	-) ;;
	unindexed:*)
		C_FASTA_STAGED="writable"
		C_FASTA_PATH=$(resolve_path "${C_FASTA#unindexed:}")
		;;
	unindexed-readonly:*)
		C_FASTA_STAGED="readonly"
		C_FASTA_PATH=$(resolve_path "${C_FASTA#unindexed-readonly:}")
		;;
	*) C_FASTA_PATH=$(resolve_path "$C_FASTA") ;;
	esac
	C_REFUSAL_DIR=""
	[[ "$C_REFUSAL" == "-" ]] || C_REFUSAL_DIR=$(resolve_path "$C_REFUSAL")
	# The run of its manifest the case compares under (the header's NAMED-RUNS TABLE): C_RUN is
	# the name the comparator is told, or nothing; C_RUN_SOURCE is `table` for a named row,
	# `match` for the run the flags fit, `pooled` for a pooled row, else nothing; C_RUN_REASON
	# is a row's reason; C_RUN_ERROR is what preflight refuses.
	C_RUN=""
	C_RUN_SOURCE=""
	C_RUN_REASON=""
	C_RUN_ERROR=""
	input_name=$(basename "$C_INPUT_PATH")
	# The case and the run are compared as text (`k ""`): two numerals of one number (7, 7.0)
	# are equal to awk as numbers.
	run_row=$(awk -F'\t' -v k="$C_ID" '$1 == k "" {print $2 "\t" $3}' <<<"$CASE_RUN_ROWS")
	run_name="${run_row%%$'\t'*}"
	if [[ -n "$C_MANIFEST" && -z "$C_REFUSAL_DIR" && "$C_TOOL_COUNT" != "true" ]]; then
		run_facts=$(RUN_SPELLINGS="$RUN_SPELLINGS" PERL_OPTIONS="$PERL_OPTIONS" CANDIDATE_OPTIONS="$CANDIDATE_OPTIONS" "$PYTHON" -c "$CASE_RUN_PROGRAM" "$C_MANIFEST" "$input_name" ${C_FLAG_WORDS[@]+"${C_FLAG_WORDS[@]}"}) || {
			# A match that fails tells nothing of the runs over the input, so the case is refused.
			run_facts=""
			C_RUN_ERROR="the run match over $C_MANIFEST failed (its error is above), so the runs over $input_name are not known and the case's run cannot be settled"
		}
		run_fit=$(tsv_lookup "$run_facts" fit)
		run_over=$(tsv_lookup "$run_facts" over)
		C_RUN_REASON="${run_row#*$'\t'}"
		if [[ -n "$C_RUN_ERROR" ]]; then
			:
		elif [[ "$run_name" == "pooled" ]]; then
			if [[ -n "$run_fit" ]]; then
				C_RUN_ERROR="$CASE_RUNS marks the case pooled, and its flags are the run '$run_fit' of $C_MANIFEST: a case that is a run of its manifest compares under it: take the row out"
			else
				C_RUN_SOURCE="pooled"
			fi
		elif [[ -n "$run_row" ]]; then
			run_input=$(awk -F'\t' -v r="$run_name" '$1 == "run" && $2 == r "" {print $3; exit}' <<<"$run_facts")
			if [[ -z "$run_input" ]]; then
				C_RUN_ERROR="$CASE_RUNS names the run '$run_name', and $C_MANIFEST lists no run of that name (its runs: $(awk -F'\t' '$1 == "run" {printf "%s%s", sep, $2; sep = ", "}' <<<"$run_facts"))"
			elif [[ "$run_input" != "$input_name" ]]; then
				C_RUN_ERROR="$CASE_RUNS names the run '$run_name', which annotates $run_input, and the case's input is $input_name: a pair is of one run over one input"
			else
				C_RUN="$run_name"
				C_RUN_SOURCE="table"
			fi
		elif [[ -n "$run_fit" ]]; then
			C_RUN="$run_fit"
			C_RUN_SOURCE="match"
		elif [[ "${run_over:-0}" -gt 1 ]]; then
			C_RUN_ERROR="$C_MANIFEST lists $run_over runs over $input_name ($(awk -F'\t' '$1 == "over" {print $3; exit}' <<<"$run_facts")) and the flags are not exactly one of them, so the comparison would read every run's documented items as one list: name the case's run in a row of ${CASE_RUNS:-a named-runs table (--case-runs)} (case, run, reason), or mark it pooled there with the reason no run is right"
		fi
	elif [[ -n "$run_row" && "$run_name" != "pooled" ]]; then
		C_RUN_ERROR="$CASE_RUNS names the run '$run_name' for the case, which compares against no manifest (no manifest.json beside its input, a refusal row or a --count tool row)"
	fi
	# A pooled case compares unnamed, which the comparator refuses where it is asked for a
	# class the lists of several runs over the input document: on a manifest of several such
	# runs a discovery pass asks the case for no class the manifest documents but those its
	# cell names. On one run or none the comparator refuses nothing, and the request stays.
	if [[ "$DISCOVER" == "true" && "$C_RUN_SOURCE" == "pooled" && "${run_over:-0}" -gt 1 && -n "$C_MANIFEST_CLASSES" ]]; then
		asked=""
		for w in $(printf '%s' "$C_REQUEST_CLASSES" | tr ',' ' '); do
			if has_word "$C_MANIFEST_CLASSES" "$w" && ! has_word "$(printf '%s' "$C_MASK_CLASSES" | tr ',' ' ')" "$w"; then
				continue
			fi
			asked="$asked $w"
		done
		C_REQUEST_CLASSES="-"
		[[ -z "$asked" ]] || C_REQUEST_CLASSES=$(printf '%s' "${asked# }" | tr ' ' ',')
	fi
	# The directory the case's Ensembl VEP outputs are read from instead of a run: its own
	# reference cell, else its directory under --reference-dir; a refusal case runs live, and
	# so does a case whose fasta cell is marked, since a stored output holds no index an engine
	# wrote.
	C_REFERENCE_DIR=""
	C_REFERENCE_FILE=""
	if [[ -n "$C_REFUSAL_DIR" || -n "$C_FASTA_STAGED" ]]; then
		:
	elif [[ "$C_REFERENCE" != "-" ]]; then
		C_REFERENCE_DIR=$(resolve_path "$C_REFERENCE")
		# A tool row's stored reference may be the filtered output file itself.
		[[ -z "$C_TOOL" || ! -f "$C_REFERENCE_DIR" ]] || C_REFERENCE_FILE="$C_REFERENCE_DIR"
	elif [[ -n "$REFERENCE_DIR" ]]; then
		C_REFERENCE_DIR="$REFERENCE_DIR/$C_ID/expected"
	fi

	C_CACHE="${cache_spec%%=*}"
	C_CACHE_VERSION="${C_CACHE%%_*}"
	C_ASSEMBLY="${C_CACHE#*_}"
	# CacheDir.pm:312 at release/116.3 appends _refseq, then _merged, for each flag set;
	# Config.pm:621 refuses the two together, so one suffix follows the species.
	species_dir="$C_SPECIES"
	if has_option --merged; then
		species_dir="${C_SPECIES}_merged"
	elif has_option --refseq; then
		species_dir="${C_SPECIES}_refseq"
	fi
	C_PERL_SLICE="$PERL_CACHE_DIR/$species_dir/$C_CACHE"
	# A layout root is what the candidate is given; the leaf it resolves is checked and recorded.
	C_DIR_CACHE=""
	if [[ "$cache_spec" == *="$LAYOUT_PREFIX"* ]]; then
		C_DIR_CACHE=$(resolve_path "${cache_spec#*="$LAYOUT_PREFIX"}")
		C_JSON_CACHE="$C_DIR_CACHE/$species_dir/$C_CACHE"
	elif [[ "$cache_spec" == *=* ]]; then
		C_JSON_CACHE=$(resolve_path "${cache_spec#*=}")
	else
		C_JSON_CACHE="$JSON_CACHE_DIR/$(lower "$C_ASSEMBLY")"
	fi

	C_FORMAT_LIST=$(printf '%s' "$C_FORMATS" | tr ',' ' ')
	C_OUTPUT_FORMATS=""
	C_WARNINGS="false"
	C_STATS="false"
	C_PERL_FORMATS=""
	local f
	for f in $C_FORMAT_LIST; do
		if [[ "$f" == "warnings" ]]; then
			C_WARNINGS="true"
		elif [[ "$f" == "stats" ]]; then
			C_STATS="true"
		else
			C_OUTPUT_FORMATS="$C_OUTPUT_FORMATS $f"
		fi
	done
	C_OUTPUT_FORMATS="${C_OUTPUT_FORMATS# }"
	# The output format whose runs' warnings files a warnings case compares: default where
	# the cell names it, else the one output format the cell names where that is tab, vcf or
	# json; empty for any other cell, which preflight refuses.
	C_WARNINGS_FORMAT=""
	if [[ "$C_WARNINGS" == "true" ]]; then
		if has_word "$C_OUTPUT_FORMATS" "default"; then
			C_WARNINGS_FORMAT="default"
		elif [[ "$C_OUTPUT_FORMATS" == "tab" || "$C_OUTPUT_FORMATS" == "vcf" || "$C_OUTPUT_FORMATS" == "json" ]]; then
			C_WARNINGS_FORMAT="$C_OUTPUT_FORMATS"
		fi
	fi
	for f in $C_OUTPUT_FORMATS; do
		f=$(perl_format "$f")
		has_word "$C_PERL_FORMATS" "$f" || C_PERL_FORMATS="$C_PERL_FORMATS $f"
	done
	C_PERL_FORMATS="${C_PERL_FORMATS# }"
	# The server a database-mode run reaches: Config.pm:287-289 at release/116.3 default
	# host ensembldb.ensembl.org and port 3306; Config.pm:711-713 set port 3337 for GRCh37;
	# Config.pm:391-394 set the Ensembl Genomes server under --genomes, whatever the species.
	C_PERL_SERVER=""
	if [[ "$C_PERL_MODE" == "database" && -n "$C_TOOL" ]]; then
		C_PERL_SERVER="$TOOL_SERVER"
	elif [[ "$C_PERL_MODE" == "database" ]]; then
		if has_option --genomes; then
			C_PERL_SERVER="mysql-eg-publicsql.ebi.ac.uk:4157"
		elif [[ "$(lower "$C_ASSEMBLY")" == "grch37" ]]; then
			C_PERL_SERVER="ensembldb.ensembl.org:3337"
		else
			C_PERL_SERVER="ensembldb.ensembl.org:3306"
		fi
	fi
}

# The output file of a format inside an engine's output directory; parquet is a directory.
format_file() {
	case "$1" in
	default) echo "default.txt" ;;
	tab) echo "tab.txt" ;;
	vcf) echo "vcf.vcf" ;;
	json) echo "json.jsonl" ;;
	parquet) echo "parquet.parquet" ;;
	esac
}

# perl_format <format>: the Ensembl VEP output a candidate format is compared with; the
# parquet comparison reads the tab output.
perl_format() {
	if [[ "$1" == "parquet" ]]; then
		echo "tab"
	else
		echo "$1"
	fi
}

# plain_or_gz <path>: <path> if it is a file, else <path>.gz if that is, else nothing.
plain_or_gz() {
	local p
	for p in "$1" "$1.gz"; do
		if [[ -f "$p" ]]; then
			printf '%s' "$p"
			return 0
		fi
	done
	return 1
}
# stored_reference_file <format>: the loaded case's stored output, plain or .gz, or nothing.
stored_reference_file() {
	if [[ -n "$C_REFERENCE_FILE" ]]; then
		printf '%s' "$C_REFERENCE_FILE"
		return 0
	fi
	plain_or_gz "$C_REFERENCE_DIR/$(format_file "$1")"
}
# stored_warnings_file: the stored warnings file of the run whose warnings the loaded case
# compares, plain or .gz, or nothing.
stored_warnings_file() {
	plain_or_gz "$C_REFERENCE_DIR/$(format_file "$C_WARNINGS_FORMAT")_warnings.txt"
}
# stored_stats_file: the loaded case's stored default-run statistics report, plain or .gz, or
# nothing; a golden corpus keeps Ensembl VEP's `<out>_summary.txt` as `<format>.summary.log`.
stored_stats_file() {
	plain_or_gz "$C_REFERENCE_DIR/default.summary.log"
}
# case_work_dir: the loaded case's directory under the work tree. Cases live under cases/,
# candidates under <case>/candidates/, so no case or pin name can reach the Ensembl VEP
# outputs (<case>/perl/), the decompressed FASTAs (fasta/) or the logs (logs/). The staged
# copies of a marked fasta cell keep the same split under <case>/fasta/: perl/<format>/ and
# candidates/<pin>/<format>/.
case_work_dir() { printf '%s/cases/%s' "$WORK_DIR" "$C_ID"; }

# expected_output <perl format>: the Ensembl VEP output a comparison of the loaded case
# reads: the stored reference, or the file its run writes under the work tree.
expected_output() {
	if [[ -n "$C_REFERENCE_DIR" ]]; then
		stored_reference_file "$1"
	else
		printf '%s/perl/%s' "$(case_work_dir)" "$(format_file "$1")"
	fi
}

# Preflight: every selected case is checked before anything runs
PREFLIGHT_FAILED=0
pf_error() {
	echo "ERROR: [feature_oracle] $C_ID: $*" >&2
	PREFLIGHT_FAILED=1
}
# FASTA_SOURCES: the FASTA path of every fasta cell without a marker, one per line, in
# first-use order; a .gz source is decompressed into fasta/<its line number>/ under the
# work tree. TSS_SLICES: the cache
# slice of every offline --nearest case whose first Ensembl VEP run writes
# transcript_gene_tss.txt, one line per case. IMAGES: the image of every case that runs
# Ensembl VEP live, and N_LIVE their count.
FASTA_SOURCES=""
TSS_SLICES=""
IMAGES=""
N_LIVE=0
# The selected cases that name stats.
N_STATS=0
# The selected cases told a run by the named-runs table and by their flags, and those a row
# of the table keeps unnamed.
N_RUN_TABLE=0
N_RUN_MATCH=0
N_RUN_POOLED=0
# An assembly name may carry a dot or a dash (BDGP6.54, R64-1-1).
CACHE_RE='^[0-9]+_[A-Za-z0-9][A-Za-z0-9._-]*$'
# An Ensembl species production name, the cache directory's name.
SPECIES_RE='^[a-z][a-z0-9_]*$'
FORMAT_SPEC_RE='^[a-z,]+$'
# The classes the selected cases' manifests document, which --discover-classes may name
# beside the built-in ones.
DISCOVER_MANIFEST_CLASSES=""

while IFS= read -r line; do
	[[ -n "$line" ]] || continue
	load_case "$line"
	DISCOVER_MANIFEST_CLASSES="$DISCOVER_MANIFEST_CLASSES $C_MANIFEST_CLASSES"
	[[ -z "$C_RUN_ERROR" ]] || pf_error "$C_RUN_ERROR"
	case "$C_RUN_SOURCE" in
	table) N_RUN_TABLE=$((N_RUN_TABLE + 1)) ;;
	match) N_RUN_MATCH=$((N_RUN_MATCH + 1)) ;;
	pooled) N_RUN_POOLED=$((N_RUN_POOLED + 1)) ;;
	esac
	[[ -f "$C_INPUT_PATH" ]] || pf_error "input not found: $C_INPUT_PATH (from '$C_INPUT')"
	if [[ -n "$C_FASTA_STAGED" ]]; then
		# A marked case needs its source alone: every engine run reads a copy of its own, so no
		# index is looked for and the source joins no shared decompression.
		[[ -f "$C_FASTA_PATH" ]] || pf_error "FASTA not found: $C_FASTA_PATH"
		if [[ "$C_FASTA_STAGED" == "readonly" && "$C_REFUSAL" == "-" ]]; then
			pf_error "fasta '$C_FASTA' on a row whose refusal cell is -; a directory no engine can write is a refusal row's, so the row names its refusal directory"
		fi
		[[ "$C_REFERENCE" == "-" ]] || pf_error "fasta '$C_FASTA' with reference '$C_REFERENCE'; each engine indexes a copy of its own in a live run, which a stored output does not hold, so the reference cell is -"
		# A directory's mode binds no process that holds CAP_DAC_OVERRIDE, as user 0 does;
		# candidate_cmd takes the capability from the candidate through setpriv.
		if [[ "$C_FASTA_STAGED" == "readonly" && "$(id -u)" == "0" ]] && ! command -v setpriv >/dev/null 2>&1; then
			pf_error "fasta '$C_FASTA' as user 0 without setpriv on PATH; the candidate of a read-only case runs through setpriv --bounding-set=-dac_override there, since a directory's mode binds no process that holds the capability"
		fi
	elif [[ -n "$C_FASTA_PATH" ]]; then
		[[ -f "$C_FASTA_PATH" ]] || pf_error "FASTA not found: $C_FASTA_PATH"
		# samtools faidx names a compressed FASTA's index <name>.fa.gz.fai; a committed corpus
		# also carries the plain spelling, <name>.fa.fai, which the decompressed copy reads.
		fai="$C_FASTA_PATH.fai"
		[[ -f "$fai" ]] || fai="${C_FASTA_PATH%.gz}.fai"
		[[ -f "$fai" ]] || pf_error "FASTA index not found: $fai"
		in_lines "$FASTA_SOURCES" "$C_FASTA_PATH" || FASTA_SOURCES="${FASTA_SOURCES}${C_FASTA_PATH}
"
	fi
	[[ "$C_CACHE" =~ $CACHE_RE ]] || pf_error "cache '$C_CACHE' is not <version>_<assembly>"
	[[ "$C_SPECIES" =~ $SPECIES_RE ]] || pf_error "species '$C_SPECIES' is not an Ensembl species name (lower-case letters, digits and underscores)"
	# CacheDir.pm:310-311 at release/116.3 refuses a suffixed species: the flag selects the slice.
	case "$C_SPECIES" in
	*_refseq | *_merged) pf_error "species '$C_SPECIES' carries the slice suffix; write the species and the --${C_SPECIES##*_} flag" ;;
	esac
	# A tool row reads neither cache: its cache cell is recorded alone.
	if [[ -n "$C_TOOL" ]]; then
		:
	elif [[ -n "$C_DIR_CACHE" && ! -d "$C_JSON_CACHE/transcripts" ]]; then
		pf_error "cache root $C_DIR_CACHE has no $(basename "$(dirname "$C_JSON_CACHE")")/$C_CACHE/transcripts/ (the layout the candidate resolves from --dir_cache)"
	elif [[ ! -d "$C_JSON_CACHE/transcripts" ]]; then
		pf_error "candidate JSON cache has no transcripts/: $C_JSON_CACHE"
	fi
	case "$C_PERL_MODE" in
	offline)
		if [[ -n "$C_TOOL" ]]; then
			:
		elif [[ -z "$C_REFERENCE_DIR" && ! -d "$C_PERL_SLICE" ]]; then
			pf_error "Ensembl VEP cache slice not found: $C_PERL_SLICE"
		elif [[ -z "$C_REFERENCE_DIR" ]] && tss_first_run; then
			TSS_SLICES="${TSS_SLICES}${C_PERL_SLICE}
"
		fi
		;;
	database) ;;
	*) pf_error "perl_mode '$C_PERL_MODE' is not offline or database" ;;
	esac
	[[ "$C_PERL_IMAGE" != "-" ]] || pf_error "perl_image is -"
	if [[ -z "$C_REFERENCE_DIR" ]]; then
		N_LIVE=$((N_LIVE + 1))
		has_word "$IMAGES" "$C_PERL_IMAGE" || IMAGES="$IMAGES $C_PERL_IMAGE"
	fi
	[[ "$C_REFERENCE_RELEASE" != "-" ]] || pf_error "reference_release is -"
	[[ "$C_UNORDERED" == "-" || "$C_UNORDERED" == "1" ]] || pf_error "unordered '$C_UNORDERED' is neither - nor 1"
	if [[ -n "$C_REFUSAL_DIR" ]]; then
		[[ "$C_REFERENCE" == "-" ]] || pf_error "refusal and reference are both set; a refused run has no output to store"
		[[ "$C_WARNINGS" != "true" ]] || pf_error "refusal with warnings in formats; a refused run writes no warnings file"
		[[ "$C_STATS" != "true" ]] || pf_error "refusal with stats in formats; a refused run writes no statistics report"
		for column_cell in "masks:$C_MASKS" "floors:$C_FLOORS" "field_pairs:$C_FIELD_PAIRS" "unordered:$C_UNORDERED"; do
			[[ "${column_cell#*:}" == "-" ]] || pf_error "refusal with ${column_cell%%:*} '${column_cell#*:}'; a refused run compares nothing, so the cell is -"
		done
		if [[ ! -d "$C_REFUSAL_DIR" ]]; then
			pf_error "refusal directory not found: $C_REFUSAL_DIR (from '$C_REFUSAL')"
		else
			for f in perl.txt candidate.txt; do
				if [[ ! -f "$C_REFUSAL_DIR/$f" ]]; then
					pf_error "refusal directory lacks $f: $C_REFUSAL_DIR"
				elif [[ -z "$(tr -d '[:space:]' <"$C_REFUSAL_DIR/$f")" ]]; then
					pf_error "refusal text $C_REFUSAL_DIR/$f is empty, which every stderr would match"
				fi
			done
		fi
	fi
	if [[ ! "$C_FORMATS" =~ $FORMAT_SPEC_RE ]] || [[ -z "$C_OUTPUT_FORMATS" ]]; then
		pf_error "formats '$C_FORMATS' names no output format (default,tab,vcf,json,parquet, and warnings or stats)"
	fi
	seen=""
	for f in $C_FORMAT_LIST; do
		if ! has_word "$OUTPUT_FORMATS warnings stats" "$f"; then
			pf_error "formats names '$f'; the formats are default,tab,vcf,json,parquet, warnings and stats"
		elif has_word "$seen" "$f"; then
			pf_error "formats names '$f' twice"
		fi
		seen="$seen $f"
	done
	if [[ "$C_WARNINGS" == "true" && -z "$C_WARNINGS_FORMAT" && -n "$C_OUTPUT_FORMATS" ]]; then
		pf_error "formats names warnings without default beside '$C_OUTPUT_FORMATS'; the warnings compared are the default run's, or those of the one run of tab, vcf or json the cell names"
	fi
	if [[ "$C_STATS" == "true" ]]; then
		N_STATS=$((N_STATS + 1))
		has_word "$C_OUTPUT_FORMATS" "default" || pf_error "formats names stats without default; the statistics compared are the default run's"
		# -sf is the option's other name (`stats_file|sf=s`, Config.pm:169 at release/116.3).
		if has_option --stats_file || has_option --sf; then
			pf_error "flags names --stats_file with stats in formats; the statistics compared are the default run's <out>_summary.txt"
		fi
	fi
	for ((k = 0; k < ${#C_FLAG_WORDS[@]}; k++)); do
		f="${C_FLAG_WORDS[$k]}"
		# The option the word names, in any spelling vep reads; the list spells a short flag
		# with one dash.
		flag_name=$(option_name "$f")
		if [[ -n "$flag_name" ]] && { has_word "$RUNNER_FLAGS" "$flag_name" || has_word "$RUNNER_FLAGS" "${flag_name#-}"; }; then
			pf_error "flags names $f, which the runner sets from the table's columns and its own options"
		elif [[ -n "$flag_name" ]]; then
			# A word the list does not spell is refused when an engine resolves it to a listed
			# option, as another of the option's names or as an abbreviation.
			runner_word=$(OPTION_NAME="${flag_name#--}" awk -v runner=" $RUNNER_FLAGS " -v perl="$PERL_OPTIONS" -v candidate="$CANDIDATE_OPTIONS" "$RUNNER_WORD_PROGRAM")
			[[ -z "$runner_word" ]] || pf_error "flags names $f, which vep reads as $runner_word, a word the runner sets from the table's columns and its own options"
		fi
		# The value of --dir is the rest of the word after `=`, else the next word whatever it
		# begins with (`dir=s`, Config.pm:241 at release/116.3); the candidate's word at the
		# same index carries the host directory.
		[[ "$flag_name" == "--dir" ]] || continue
		if [[ "$f" == *=* ]]; then
			dir_value="${f#*=}"
			dir_host="${C_CAND_FLAG_WORDS[$k]#*=}"
		else
			dir_value="${C_FLAG_WORDS[$((k + 1))]-}"
			dir_host="${C_CAND_FLAG_WORDS[$((k + 1))]-}"
		fi
		dir_root="${dir_value%%/*}"
		dir_rest="${dir_value#"$dir_root"}"
		# A manifest's word after the first segment is read as a path under {corpus}: a second
		# placeholder, refused as one.
		if [[ ("$dir_root" != "$CORPUS_PLACEHOLDER" && "$dir_root" != "$SOURCES_PLACEHOLDER") || "/$dir_value/" == */../* || "$dir_rest" == *"$CORPUS_PLACEHOLDER"* || "$dir_rest" == *"$SOURCES_PLACEHOLDER"* || "$dir_rest" == *"$PLUGINS_PLACEHOLDER"* || "$dir_rest" == *"$GENERATED_PLACEHOLDER"* ]]; then
			pf_error "flags give --dir '$dir_value', which is not $CORPUS_PLACEHOLDER, $SOURCES_PLACEHOLDER or a path under one without a .. segment or a second placeholder: a flags cell reaches no other directory of Ensembl VEP's container"
		# A {sources} directory has a host path only under --sources-dir, whose absence is
		# reported below.
		elif [[ "$dir_root" == "$CORPUS_PLACEHOLDER" || -n "$SOURCES_DIR" ]] && [[ -e "$dir_host/vep.ini" || -L "$dir_host/vep.ini" ]]; then
			pf_error "flags give --dir '$dir_value', which holds a vep.ini ($dir_host/vep.ini): both engines read that file as options (Config.pm:704-705 at release/116.3), so it is refused as --config is"
		fi
	done
	if [[ "$C_FLAGS" == *"$SOURCES_PLACEHOLDER"* && -z "$SOURCES_DIR" ]]; then
		pf_error "flags carry $SOURCES_PLACEHOLDER and no --sources-dir names the directory it stands for"
	fi
	if [[ -z "$C_MANIFEST" && ("$C_FLAGS" == *"$PLUGINS_PLACEHOLDER"* || "$C_FLAGS" == *"$GENERATED_PLACEHOLDER"*) ]]; then
		pf_error "flags carry $PLUGINS_PLACEHOLDER or $GENERATED_PLACEHOLDER, a corpus manifest's words for two directories of its corpus, and no manifest.json sits beside the input"
	fi
	# On the cell, not C_TOOL: a cell whose tool name is empty (`:<args>`) loads as no tool
	# and would otherwise run as an engine row.
	if [[ "$C_TOOL_CELL" != "-" ]]; then
		if ! has_word "$TOOLS" "$C_TOOL" || [[ "$C_TOOL_CELL" != *:* ]]; then
			pf_error "tool '$C_TOOL_CELL' is not - or <tool>:<args> with a tool of: $TOOLS"
		elif [[ -n "$C_TOOL_WORDS_ERROR" ]]; then
			pf_error "tool arguments '${C_TOOL_CELL#*:}' are not shell words: $C_TOOL_WORDS_ERROR"
		elif [[ "${#C_TOOL_WORDS[@]}" -eq 0 ]]; then
			pf_error "tool '$C_TOOL_CELL' carries no argument"
		fi
		for f in ${C_TOOL_WORDS[@]+"${C_TOOL_WORDS[@]}"}; do
			has_word "$TOOL_FLAGS" "${f%%=*}" && pf_error "tool arguments name $f, which the runner sets (the input, the output and its overwrite, the database connection) or whose output never reaches the output file (--list)"
		done
		if [[ "$C_PERL_MODE" == "offline" ]]; then
			for f in ${C_TOOL_WORDS[@]+"${C_TOOL_WORDS[@]}"}; do
				[[ "$f" != "--ontology" && "$f" != "-y" ]] || pf_error "tool arguments carry $f under perl_mode offline; filter_vep answers it from Ensembl's ontology database, so the row runs in database mode"
			done
		fi
		if [[ "$C_FORMATS" == *,* ]] || ! has_word "$TOOL_FORMATS" "$C_FORMATS"; then
			pf_error "formats '$C_FORMATS' on a tool row; the tool reads one annotated file and writes one, in one of: $TOOL_FORMATS"
		fi
		[[ -z "$C_FLAGS" ]] || pf_error "flags '$C_FLAGS' on a tool row; the row runs no engine, so the flags cell is -"
		[[ "$C_FASTA" == "-" ]] || pf_error "fasta '$C_FASTA' on a tool row; the tool reads no FASTA, so the cell is -"
		if [[ "$C_TOOL_COUNT" == "true" ]]; then
			for column_cell in "masks:$C_MASKS" "floors:$C_FLOORS" "field_pairs:$C_FIELD_PAIRS" "unordered:$C_UNORDERED"; do
				[[ "${column_cell#*:}" == "-" ]] || pf_error "a --count tool row with ${column_cell%%:*} '${column_cell#*:}'; the two counts are compared as text, so the cell is -"
			done
		fi
	fi
	if [[ -n "$C_REFERENCE_DIR" && -z "$C_REFERENCE_FILE" && ! -d "$C_REFERENCE_DIR" ]]; then
		pf_error "stored reference directory not found: $C_REFERENCE_DIR${C_TOOL:+; a tool row may name the stored filtered output file instead}"
	elif [[ -n "$C_REFERENCE_DIR" ]]; then
		for f in $C_PERL_FORMATS; do
			stored_reference_file "$f" >/dev/null || pf_error "stored reference lacks $(format_file "$f")[.gz] under $C_REFERENCE_DIR/"
		done
		if [[ -n "$C_WARNINGS_FORMAT" ]]; then
			stored_warnings_file >/dev/null || pf_error "stored reference lacks $(format_file "$C_WARNINGS_FORMAT")_warnings.txt[.gz] under $C_REFERENCE_DIR/ (an empty file states no warnings)"
		fi
		if [[ "$C_STATS" == "true" ]]; then
			stored_stats_file >/dev/null || pf_error "stored reference lacks default.summary.log[.gz] under $C_REFERENCE_DIR/"
		fi
	fi
	# In a comparator cell `-` is the whole cell or nothing: no mask class, floor or field pair
	# begins with a dash, and a cell that does is read by the comparator's option parser as
	# another option.
	for column_cell in "masks:$C_MASKS" "floors:$C_FLOORS" "field_pairs:$C_FIELD_PAIRS"; do
		column="${column_cell%%:*}"
		value="${column_cell#*:}"
		items=$(printf '%s' "$value" | tr ',;' '\n\n')
		if [[ "$value" != "-" && "$items" =~ (^|$'\n')- ]]; then
			pf_error "$column '$value' has an item beginning with -; write - alone for no value"
		fi
	done
	if [[ "$C_MASKS" != "-" && ( "$C_MASKS" == ,* || "$C_MASKS" == *, || "$C_MASKS" == *,,* || "$C_MASKS" == *[[:space:]]* ) ]]; then
		pf_error "masks '$C_MASKS' has an empty or blank item"
	fi
	# `read -a` splits the cell at its commas without expanding a glob, so a `*` item is
	# refused as written instead of being read as the names of files in the working directory.
	IFS=',' read -r -a mask_items <<<"$C_MASKS"
	seen_mask_classes=""
	for f in "${mask_items[@]}"; do
		if [[ "$f" != "-" && ! "$f" =~ ^[a-z0-9_]+(=[1-9][0-9]*)?$ ]]; then
			pf_error "masks item '$f' is not <class> or <class>=<expected applied count>"
		fi
		f="${f%%=*}"
		# An item with no class name is refused above; has_word tests the padded list for two
		# spaces when the word is empty, so both checks below would report a class with no name.
		[[ -n "$f" ]] || continue
		if has_word "$seen_mask_classes" "$f"; then
			pf_error "masks names '$f' twice"
		fi
		seen_mask_classes="$seen_mask_classes $f"
		if [[ "$f" != "-" && -z "$C_MANIFEST" ]] && ! has_word "$BUILTIN_MASKS" "$f"; then
			if has_word "$MANIFEST_ONLY_MASKS" "$f"; then
				pf_error "masks names '$f', a manifest-only class (the comparator's MANIFEST_CLASSES: no rule over the two outputs decides it), and no manifest.json stands beside the input $C_INPUT_PATH to list the records it covers"
			else
				pf_error "masks names '$f', which applies from a manifest.json beside the input, and none is there (the built-in classes are $BUILTIN_MASKS)"
			fi
		fi
		# A reference-reading class of CACHE_READING_MASKS is decided from the cache where the
		# case names no FASTA, and leaves only its rows that read the FASTA a difference.
		if [[ -z "$C_MANIFEST" && -z "$C_FASTA_PATH" ]] && has_word "$REFERENCE_READING_MASKS" "$f" && ! has_word "$CACHE_READING_MASKS" "$f"; then
			pf_error "masks names '$f', which is decided from the case's FASTA or applies from a manifest.json beside the input, and the case has neither"
		fi
		if [[ -n "$C_FASTA_STAGED" ]] && has_word "$REFERENCE_READING_MASKS" "$f"; then
			pf_error "masks names '$f' beside fasta '$C_FASTA'; the comparator reads the case's FASTA through its index, and a marked cell's FASTA has none beside it"
		fi
		if has_word "$SCORE_READING_MASKS" "$f" && [[ ${#C_SPLICEAI_TABLES[@]} -eq 0 ]]; then
			pf_error "masks names '$f', which reads the score files the flags hand the SpliceAI plugin (snv=, indel=), and the flags name none"
		fi
	done
done <<<"$SELECTED"
for f in $(printf '%s' "$DISCOVER_CLASSES" | tr ',' ' '); do
	if ! has_word "$BUILTIN_MASKS $DISCOVER_MANIFEST_CLASSES" "$f"; then
		echo "ERROR: [feature_oracle] --discover-classes names '$f', which is neither a built-in class nor a class a selected case's manifest.json documents" >&2
		PREFLIGHT_FAILED=1
	fi
done
[[ "$PREFLIGHT_FAILED" -eq 0 ]] || die "preflight failed; nothing ran"
if [[ "$DRY_RUN" != "true" && "$N_LIVE" -gt 0 ]]; then
	command -v docker >/dev/null 2>&1 || die "docker is required to run Ensembl VEP on the $N_LIVE selected case(s) without a stored reference (pass --reference-dir, or --dry-run)"
fi
# Without flock, run_perl cannot make the --nearest cases of one slice take turns, and two
# first runs at once would tear transcript_gene_tss.txt; one such case per slice needs no turn.
if [[ "$JOBS" -gt 1 && "$DRY_RUN" != "true" ]] && ! command -v flock >/dev/null 2>&1; then
	slice=$(awk 'seen[$0]++ {print; exit}' <<<"$TSS_SLICES")
	[[ -z "$slice" ]] || die "flock is required to run the --nearest cases of $slice under --jobs $JOBS while the slice lacks transcript_gene_tss.txt (install flock, or pass --jobs 1)"
fi

# The statistics comparator, looked for where a selected case names stats: the build of
# tools/vep-stats-compare that VEP_STATS_COMPARATOR names, else this tree's release build,
# resolved as the format comparator is.
STATS_COMPARATOR=""
if [[ "$N_STATS" -gt 0 ]]; then
	STATS_COMPARATOR_NAMED="${VEP_STATS_COMPARATOR:-$REPO_DIR/tools/vep-stats-compare/target/release/vep-stats-compare}"
	if STATS_COMPARATOR=$(type -P -- "$STATS_COMPARATOR_NAMED") && [[ -f "$STATS_COMPARATOR" && -x "$STATS_COMPARATOR" ]]; then
		STATS_COMPARATOR=$(abs_path "$STATS_COMPARATOR")
	else
		die "the statistics comparator $STATS_COMPARATOR_NAMED is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-stats-compare/Cargo.toml' or name a built one in VEP_STATS_COMPARATOR"
	fi
fi

# Commands

# add_stats_flag <format>: appends the run's statistics flag to CMD: --stats_text on the
# default run of a case comparing stats (nothing when the case's flags carry it already, so
# the flag words add it once), --no_stats on every other run.
add_stats_flag() {
	if [[ "$1" == "default" && "$C_STATS" == "true" ]]; then
		has_option --stats_text || CMD+=(--stats_text)
	else
		CMD+=(--no_stats)
	fi
}

# perl_cmd <format> <host out dir> <host fasta or empty>: CMD becomes the docker command.
# The FASTA of a case whose fasta cell is marked is the run's staged copy.
perl_cmd() {
	local fmt="$1" out_dir="$2" fasta="$3"
	local in_name fasta_mount
	in_name=$(basename "$C_INPUT_PATH")
	# Parser/VCF.pm:639 at release/116.3 orders an individual's alternate alleles, and
	# OutputFactory.pm:983 the ZYG list of --individual_zyg, by Perl hash order; the committed
	# --individual and --individual_zyg references are runs under this seed.
	CMD=(docker run --rm --user "$(id -u):$(id -g)" -e PERL_HASH_SEED=0 -e PERL_PERTURB_KEYS=0
		-v "${C_INPUT_PATH}:/work/input/${in_name}:ro"
		-v "${out_dir}:/work/out")
	# A tool row runs the companion tool over the input in place of vep: no cache, no FASTA,
	# the connection flags filter_vep names (--version is the database version) in database
	# mode, then the row's words.
	if [[ -n "$C_TOOL" ]]; then
		CMD+=("$C_PERL_IMAGE" "$C_TOOL" -i "/work/input/${in_name}" -o "/work/out/$(format_file "$fmt")" --force_overwrite)
		if [[ "$C_PERL_MODE" == "database" ]]; then
			CMD+=(--host "${C_PERL_SERVER%%:*}" --port "${C_PERL_SERVER##*:}" --user anonymous --version "$C_CACHE_VERSION")
		fi
		CMD+=(${C_TOOL_WORDS[@]+"${C_TOOL_WORDS[@]}"})
		return
	fi
	if [[ "$C_PERL_MODE" == "offline" ]]; then
		# Ensembl VEP's first --nearest run writes transcript_gene_tss.txt into the cache
		# directory and dies on a read-only one (AnnotationSource/Cache/Transcript.pm:288-289).
		if has_option --nearest; then
			CMD+=(-v "${PERL_CACHE_DIR}:/work/cache")
		else
			CMD+=(-v "${PERL_CACHE_DIR}:/work/cache:ro")
		fi
	fi
	# A staged copy is alone in its directory, which is mounted whole so that Ensembl VEP
	# writes its lock file and its index beside the file (FastaSequence.pm:334-346); the
	# read-only form mounts the directory read-only.
	if [[ -n "$fasta" && -n "$C_FASTA_STAGED" ]]; then
		fasta_mount="$(dirname "$fasta"):/work/fasta"
		[[ "$C_FASTA_STAGED" != "readonly" ]] || fasta_mount="$fasta_mount:ro"
		CMD+=(-v "$fasta_mount")
	elif [[ -n "$fasta" ]]; then
		CMD+=(-v "${fasta}:/work/ref.fa:ro" -v "${fasta}.fai:/work/ref.fa.fai:ro")
	fi
	if [[ "$C_FLAGS" == *"$CORPUS_PLACEHOLDER"* || (-n "$C_MANIFEST" && ("$C_FLAGS" == *"$PLUGINS_PLACEHOLDER"* || "$C_FLAGS" == *"$GENERATED_PLACEHOLDER"*)) ]]; then
		CMD+=(-v "${C_CORPUS_DIR}:${CORPUS_MOUNT}:ro")
	fi
	if [[ "$C_FLAGS" == *"$SOURCES_PLACEHOLDER"* ]]; then
		CMD+=(-v "${SOURCES_DIR}:${SOURCES_MOUNT}:ro")
	fi
	CMD+=("$C_PERL_IMAGE" vep -i "/work/input/${in_name}" -o "/work/out/$(format_file "$fmt")")
	if [[ "$C_PERL_MODE" == "offline" ]]; then
		CMD+=(--offline --cache --dir_cache /work/cache --cache_version "$C_CACHE_VERSION")
	else
		CMD+=(--database --db_version "$C_CACHE_VERSION")
	fi
	CMD+=(--species "$C_SPECIES" --assembly "$C_ASSEMBLY")
	has_option --format || CMD+=(--format vcf)
	CMD+=(--fork "$PERL_FORK")
	add_stats_flag "$fmt"
	CMD+=(--quiet --force_overwrite)
	CMD+=(${C_PERL_FLAG_WORDS[@]+"${C_PERL_FLAG_WORDS[@]}"})
	if [[ -n "$fasta" && -n "$C_FASTA_STAGED" ]]; then
		CMD+=(--fasta "/work/fasta/${fasta##*/}")
	elif [[ -n "$fasta" ]]; then
		CMD+=(--fasta /work/ref.fa)
	fi
	case "$fmt" in
	default) ;;
	*) CMD+=("--$fmt") ;;
	esac
}

# candidate_cmd <binary> <format> <host out dir> <host fasta or empty>
candidate_cmd() {
	local bin="$1" fmt="$2" out_dir="$3" fasta="$4"
	# The candidate's companion of a tool is its subcommand: `vep filter` for filter_vep. It
	# reads no database, so the row's words alone follow.
	if [[ -n "$C_TOOL" ]]; then
		case "$C_TOOL" in
		filter_vep) CMD=("$bin" filter) ;;
		*) die "$C_ID: tool '$C_TOOL' is in TOOLS and candidate_cmd names no subcommand for it" ;;
		esac
		CMD+=(-i "$C_INPUT_PATH" -o "$out_dir/$(format_file "$fmt")" --force_overwrite)
		CMD+=(${C_TOOL_WORDS[@]+"${C_TOOL_WORDS[@]}"})
		return
	fi
	# A directory's mode binds no process that holds CAP_DAC_OVERRIDE, as user 0 does: there
	# the candidate of a read-only case runs with the capability out of its bounding set.
	CMD=()
	if [[ "$C_FASTA_STAGED" == "readonly" && "$(id -u)" == "0" ]]; then
		CMD=(setpriv --bounding-set=-dac_override)
	fi
	CMD+=("$bin" -i "$C_INPUT_PATH" -o "$out_dir/$(format_file "$fmt")" --offline)
	if [[ -n "$C_DIR_CACHE" ]]; then
		CMD+=(--dir_cache "$C_DIR_CACHE" --species "$C_SPECIES" --cache_version "$C_CACHE_VERSION" --assembly "$C_ASSEMBLY")
	else
		CMD+=(--json_cache "$C_JSON_CACHE" --species "$C_SPECIES" --assembly "$C_ASSEMBLY")
	fi
	has_option --format || CMD+=(--format vcf)
	[[ -z "$FORK" ]] || CMD+=(--fork "$FORK")
	add_stats_flag "$fmt"
	CMD+=(--quiet --force_overwrite)
	# --genomes selects the Ensembl VEP run's server, as --database and --db_version do; the
	# candidate reads a local cache and refuses the flag by name, so a word naming it in any
	# spelling vep reads is not given it.
	local w name
	for w in ${C_CAND_FLAG_WORDS[@]+"${C_CAND_FLAG_WORDS[@]}"}; do
		name=$(option_name "$w")
		[[ -n "$name" && "--genomes" == "$name"* ]] || CMD+=("$w")
	done
	[[ -z "$fasta" ]] || CMD+=(--fasta "$fasta")
	case "$fmt" in
	default) ;;
	parquet) CMD+=(--output_format parquet) ;;
	*) CMD+=("--$fmt") ;;
	esac
}

# compare_cmd <format> <expected> <actual> <report>
compare_cmd() {
	local fmt="$1" mask_class mask_reads_input="false" mask_reads_cache="false" mask_reads_catalog="false" catalog
	local mask_classes=()
	local mask_reads_reference="false" mask_reads_scores="false" reference_sequence fa head table
	IFS=',' read -r -a mask_classes <<<"$C_REQUEST_CLASSES"
	for mask_class in "${mask_classes[@]}"; do
		has_word "$REFERENCE_READING_MASKS" "$mask_class" && mask_reads_reference="true"
		has_word "$INPUT_READING_MASKS" "$mask_class" && mask_reads_input="true"
		has_word "$CACHE_READING_MASKS" "$mask_class" && mask_reads_cache="true"
		has_word "$CATALOG_READING_MASKS" "$mask_class" && mask_reads_catalog="true"
		has_word "$SCORE_READING_MASKS" "$mask_class" && mask_reads_scores="true"
	done
	CMD=("$COMPARATOR" --format "$fmt" --expected "$2" --actual "$3" --report "$4")
	# --input where the comparator reads the records (the vcf and json formats, an
	# input-reading class) and beside a manifest, which documents each flag set under its
	# input's name and is read whole when the input is unnamed; a text-format comparison
	# reads none of the input's lines for the manifest.
	if [[ "$fmt" == vcf || "$fmt" == json || "$mask_reads_input" == "true" || -n "$C_MANIFEST" ]]; then
		CMD+=(--input "$C_INPUT_PATH")
	fi
	[[ -z "$C_MANIFEST" ]] || CMD+=(--manifest "$C_MANIFEST")
	# Beside a manifest the comparator is told the candidate's run (the header's masks column
	# states both rules): its run of the manifest, and whether it has a reference sequence, as
	# vep settles it (fasta_source of crates/vep-cli/src/runner.rs, detect_fasta of
	# cache_resolve.rs, first_indexed_fasta of crates/vep-fasta/src/lib.rs).
	if [[ -n "$C_MANIFEST" ]]; then
		[[ -z "$C_RUN" ]] || CMD+=(--run "$C_RUN")
		reference_sequence=""
		if [[ -n "$C_FASTA_PATH" ]]; then
			reference_sequence="present"
		elif ! has_option --fasta_dir; then
			reference_sequence="absent"
			if ! has_option --no_fasta; then
				# Every entry of the directory is listed, as the candidate reads it
				# (indexed_fasta_files of crates/vep-fasta/src/lib.rs), and the name is tested
				# after: a glob's `*` matches no name that opens with a dot.
				for fa in "$C_JSON_CACHE"/* "$C_JSON_CACHE"/.*; do
					case "${fa##*/}" in
					*.fa | *.fa.gz) ;;
					*) continue ;;
					esac
					[[ -f "$fa" && -f "$fa.fai" ]] || continue
					# The first 18 bytes, as hex pairs between single blanks, tell a BGZF file
					# (htslib check_header, bgzf.c: the gzip magic and deflate method, FLG.FEXTRA,
					# XLEN 6 and the BC subfield of SLEN 2) from any other gzip file and from
					# plain text.
					head=$(head -c 18 "$fa" | od -An -v -tx1 | tr -s ' \n' ' ')
					head=" ${head# }"
					head="${head% } "
					case "$head" in
					" 1f 8b 08 "?[4567cdef]" "??" "??" "??" "??" "??" "??" 06 00 42 43 02 00 "??" "??" ") [[ -f "$fa.gzi" ]] || continue ;;
					" 1f 8b "*) continue ;;
					esac
					reference_sequence="present"
					break
				done
			fi
		fi
		[[ -z "$reference_sequence" ]] || CMD+=(--reference-sequence "$reference_sequence")
	fi
	[[ "$C_REQUEST_CLASSES" == "-" ]] || CMD+=(--masks "$C_REQUEST_CLASSES")
	[[ "$C_FLOORS" == "-" ]] || CMD+=(--floors "$C_FLOORS")
	[[ "$C_FIELD_PAIRS" == "-" ]] || CMD+=(--field-pairs "$C_FIELD_PAIRS")
	[[ "$C_UNORDERED" != "1" ]] || CMD+=(--unordered)
	# NEAREST compares by transcription start site: Ensembl VEP names one transcript at a
	# winning site and the candidate every transcript there (vep-formats-compare).
	if has_option --nearest; then
		CMD+=(--nearest-sites "$C_JSON_CACHE")
	fi
	[[ "$mask_reads_cache" != "true" ]] || CMD+=(--json-cache "$C_JSON_CACHE")
	if [[ "$mask_reads_catalog" == "true" ]]; then
		catalog=$(gwas_catalog_file)
		[[ -z "$catalog" ]] || CMD+=(--gwas-catalog "$catalog")
	fi
	# A reference-reading class is decided from the FASTA the engines read and the
	# candidate's cache; without a FASTA the comparator applies it from the manifest, or
	# from the cache alone where the class stands in CACHE_READING_MASKS, as it does for a
	# case whose fasta cell is marked, where a discovery pass asks for the class and the
	# FASTA has no index for the comparator to read.
	if [[ "$mask_reads_reference" == "true" && -n "$C_FASTA_PATH" && -z "$C_FASTA_STAGED" ]]; then
		[[ "$mask_reads_cache" == "true" ]] || CMD+=(--json-cache "$C_JSON_CACHE")
		CMD+=(--reference-fasta "$(case_fasta)")
	fi
	# The score files a class reads, each as the option's own occurrence.
	if [[ "$mask_reads_scores" == "true" ]]; then
		for table in ${C_SPLICEAI_TABLES[@]+"${C_SPLICEAI_TABLES[@]}"}; do
			CMD+=(--spliceai-table "$table")
		done
	fi
}
warnings_cmd() {
	CMD=("$COMPARATOR" --warnings --expected "$1" --actual "$2" --report "$3")
	[[ "$C_REQUEST_CLASSES" == "-" ]] || CMD+=(--masks "$C_REQUEST_CLASSES")
}
# stats_cmd <expected> <actual> <report>: CMD becomes the statistics comparator command.
stats_cmd() {
	CMD=("$STATS_COMPARATOR" --expected "$1" --actual "$2" --report "$3")
}
# stats_report_file <engine output dir>: the statistics report a default run under
# --stats_text leaves beside default.txt (`<output_file>_summary` plus `.txt`).
stats_report_file() {
	printf '%s/default.txt_summary.txt' "$1"
}

print_cmd() {
	printf '%q ' "${CMD[@]}"
	echo
}

# The step timer: runs the command after its two arguments (the record's path, a path on
# the filesystem to watch) as a child, hands it SIGTERM and SIGINT, one that lands before
# the child has started included, and exits with its status (128 plus the signal for one a
# signal ended), as the shell reports a waited job. The record is written once the child has
# ended; a record that cannot be written is no failure.
STEP_USAGE_PROGRAM='
import json, os, signal, subprocess, sys, threading, time
import errno
record_path, disk_path, command = sys.argv[1], sys.argv[2], sys.argv[3:]
def used_bytes():
    try:
        stat = os.statvfs(disk_path)
    except OSError:
        return None
    return (stat.f_blocks - stat.f_bfree) * stat.f_frsize
child = None
early = []
def forward(signum, frame):
    if child is None:
        early.append(signum)
        return
    try:
        os.kill(child.pid, signum)
    except OSError:
        pass
signal.signal(signal.SIGTERM, forward)
if signal.getsignal(signal.SIGINT) is not signal.SIG_IGN:
    signal.signal(signal.SIGINT, forward)
started = time.monotonic()
disk_start = used_bytes()
try:
    try:
        child = subprocess.Popen(command)
    except OSError as exc:
        # A script with no interpreter line runs under sh, as a shell runs it.
        if exc.errno != errno.ENOEXEC:
            raise
        child = subprocess.Popen(["/bin/sh"] + command)
except OSError as exc:
    sys.stderr.write(f"{command[0]}: {exc.strerror}\n")
    sys.exit(126 if isinstance(exc, PermissionError) else 127)
for signum in early:
    forward(signum, None)
peak = [disk_start]
ended = threading.Event()
def watch():
    while not ended.wait(2):
        used = used_bytes()
        if used is not None and (peak[0] is None or used > peak[0]):
            peak[0] = used
watcher = threading.Thread(target=watch, daemon=True)
watcher.start()
_, status, usage = os.wait4(child.pid, 0)
ended.set()
used = used_bytes()
if used is not None and (peak[0] is None or used > peak[0]):
    peak[0] = used
code = 128 + os.WTERMSIG(status) if os.WIFSIGNALED(status) else os.WEXITSTATUS(status)
# ru_maxrss is in bytes on Darwin and in kibibytes on Linux.
peak_rss = usage.ru_maxrss // 1024 if sys.platform == "darwin" else usage.ru_maxrss
record = {
    "wall_seconds": round(time.monotonic() - started, 3),
    "peak_rss_kib": peak_rss,
    "exit": code,
    "disk_used_start_bytes": disk_start,
    "disk_used_peak_bytes": peak[0],
}
try:
    with open(record_path + ".tmp", "w", encoding="utf-8") as handle:
        json.dump(record, handle)
    os.replace(record_path + ".tmp", record_path)
except OSError:
    pass
sys.exit(code)
'
# The version of the comparison program a sizing figure belongs to: the first 12 hex of the
# sha256 of the file that runs as the comparator.
COMPARISON_PROGRAM_VERSION="sha256:$(sha256_of "$COMPARATOR" | cut -c1-12)"

# Run setup
RUN_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
N_CANDIDATES="${#CANDIDATES[@]}"
if [[ "$TIER" == "population" && "$N_LIVE" -gt 0 ]]; then
	echo "WARN: [feature_oracle] population tier: Ensembl VEP runs live on $N_LIVE case(s) without a stored reference" >&2
fi
log "feature_oracle tier=$TIER cases=$N_SELECTED live=$N_LIVE candidates=$N_CANDIDATES table=$CASES dry_run=$DRY_RUN"
[[ "$DISCOVER" != "true" ]] || log "discovery pass: each case is asked for ${DISCOVER_CLASSES:-the registry} beside its own masks cell; no verdict is given"
if [[ $((N_RUN_TABLE + N_RUN_MATCH + N_RUN_POOLED)) -gt 0 ]]; then
	log "runs of a manifest: $N_RUN_MATCH case(s) named by their flags, $N_RUN_TABLE by ${CASE_RUNS:-no named-runs table}, $N_RUN_POOLED pooled by a row of it"
fi
log "out=$OUT_DIR work=$WORK_DIR reference_dir=${REFERENCE_DIR:-none} perl_cache=$PERL_CACHE_DIR json_cache=$JSON_CACHE_DIR"

# decompressed_fasta <.gz source>: where the run decompresses the source: a directory of
# its own, numbered by the source's line in FASTA_SOURCES, so two sources of one file name
# never share a copy.
decompressed_fasta() {
	local n
	n=$(printf '%s' "$FASTA_SOURCES" | grep -nFx -- "$1" | cut -d: -f1)
	printf '%s/fasta/%s/%s' "$WORK_DIR" "$n" "$(basename "${1%.gz}")"
}
# The FASTA an engine reads: a .gz source decompressed into the work tree at the start of
# every run (never a copy an earlier run left) with its index beside it, otherwise the
# source itself; nothing for a case whose fasta cell is marked, each of whose engine runs
# reads a staged copy.
case_fasta() {
	[[ -z "$C_FASTA_STAGED" ]] || return 0
	case "$C_FASTA_PATH" in
	"") ;;
	*.gz) decompressed_fasta "$C_FASTA_PATH" ;;
	*) printf '%s' "$C_FASTA_PATH" ;;
	esac
}

if [[ "$DRY_RUN" == "true" ]]; then
	while IFS= read -r line; do
		[[ -n "$line" ]] || continue
		load_case "$line"
		fasta=$(case_fasta)
		case_work=$(case_work_dir)
		echo "DRY-RUN $C_ID cache species=$C_SPECIES perl=$C_PERL_SLICE candidate=$C_JSON_CACHE"
		# A marked fasta cell's staged copies are printed at the paths run_case makes them.
		[[ -z "$C_FASTA_STAGED" ]] || echo "DRY-RUN $C_ID fasta: a copy without its index for every engine run (${C_FASTA%%:*})"
		case "$C_RUN_SOURCE" in
		table) echo "DRY-RUN $C_ID run $C_RUN, named by $CASE_RUNS: $C_RUN_REASON" ;;
		pooled) echo "DRY-RUN $C_ID run pooled, no run named by $CASE_RUNS: $C_RUN_REASON" ;;
		esac
		[[ -z "$C_REFUSAL_DIR" ]] || echo "DRY-RUN $C_ID refusal: $C_REFUSAL_DIR (perl.txt, candidate.txt)"
		if [[ -z "$C_REFERENCE_DIR" ]]; then
			for f in $C_PERL_FORMATS; do
				echo "DRY-RUN $C_ID perl $f"
				[[ -z "$C_FASTA_STAGED" ]] || fasta="$case_work/fasta/perl/$f/${C_FASTA_PATH##*/}"
				perl_cmd "$f" "$case_work/perl" "$fasta"
				print_cmd
			done
		else
			echo "DRY-RUN $C_ID perl stored: ${C_REFERENCE_FILE:-$C_REFERENCE_DIR/}"
		fi
		for ((i = 0; i < N_CANDIDATES; i++)); do
			pin="${PINS[$i]}"
			cand_dir="$case_work/candidates/$pin"
			for f in $C_OUTPUT_FORMATS; do
				echo "DRY-RUN $C_ID $pin $f"
				[[ -z "$C_FASTA_STAGED" ]] || fasta="$case_work/fasta/candidates/$pin/$f/${C_FASTA_PATH##*/}"
				candidate_cmd "${CANDIDATES[$i]}" "$f" "$cand_dir" "$fasta"
				print_cmd
				[[ -z "$C_REFUSAL_DIR" ]] || continue
				expected=$(expected_output "$(perl_format "$f")")
				if [[ "$C_TOOL_COUNT" == "true" ]]; then
					echo "compared as text: $expected $cand_dir/$(format_file "$f")"
					continue
				fi
				compare_cmd "$f" "$expected" "$cand_dir/$(format_file "$f")" "$OUT_DIR/$pin/$C_ID/$f.json"
				print_cmd
			done
			if [[ "$C_WARNINGS" == "true" && -z "$C_REFUSAL_DIR" ]]; then
				echo "DRY-RUN $C_ID $pin warnings"
				if [[ -n "$C_REFERENCE_DIR" ]]; then
					expected=$(stored_warnings_file)
				else
					expected="$case_work/perl/$(format_file "$C_WARNINGS_FORMAT")_warnings.txt"
				fi
				warnings_cmd "$expected" "$cand_dir/$(format_file "$C_WARNINGS_FORMAT")_warnings.txt" "$OUT_DIR/$pin/$C_ID/warnings.json"
				print_cmd
			fi
			if [[ "$C_STATS" == "true" && -z "$C_REFUSAL_DIR" ]]; then
				echo "DRY-RUN $C_ID $pin stats"
				if [[ -n "$C_REFERENCE_DIR" ]]; then
					expected=$(stored_stats_file)
				else
					expected=$(stats_report_file "$case_work/perl")
				fi
				stats_cmd "$expected" "$(stats_report_file "$cand_dir")" "$OUT_DIR/$pin/$C_ID/stats.json"
				print_cmd
			fi
		done
	done <<<"$SELECTED"
	log "dry run complete: $N_SELECTED case(s) x $N_CANDIDATES candidate(s) planned, nothing written"
	exit 0
fi

mkdir -p "$OUT_DIR" "$WORK_DIR"
OUT_DIR=$(abs_path "$OUT_DIR")
WORK_DIR=$(abs_path "$WORK_DIR")

# An interrupted run stops its cases and their engines. Every engine and comparator runs
# as a background job waited on, since bash runs a trap only once a foreground command
# has returned; under --jobs the cases are background subshells, which bash leaves running
# when the parent dies, and each engine's step timer is a child of its case's subshell and
# hands the signal on.
on_signal() {
	trap - INT TERM
	local pid
	for pid in $(jobs -p); do
		pkill -TERM -P "$pid" 2>/dev/null || true
		kill -TERM "$pid" 2>/dev/null || true
	done
	# A staged directory of a read-only fasta cell is in mode 0555 while its engine runs; each
	# gets its write bit, so the work tree an interrupted run leaves can be removed.
	chmod -R u+w "$WORK_DIR"/cases/*/fasta 2>/dev/null || true
	exit 130
}
trap on_signal INT TERM

# Decompress every .gz FASTA once per run.
while IFS= read -r src; do
	case "$src" in
	*.gz)
		target=$(decompressed_fasta "$src")
		mkdir -p "$(dirname "$target")"
		gzip -dc "$src" >"$target.tmp" || die "cannot decompress $src"
		mv -f "$target.tmp" "$target"
		# The index preflight found: <name>.fa.gz.fai, else <name>.fa.fai; both index the
		# uncompressed sequence.
		fai="$src.fai"
		[[ -f "$fai" ]] || fai="${src%.gz}.fai"
		cp -f "$fai" "$target.fai"
		;;
	esac
done <<<"$FASTA_SOURCES"

# repo_digest <image>: the first RepoDigest of a local image, or nothing.
repo_digest() {
	docker image inspect --format '{{index .RepoDigests 0}}' "$1" 2>/dev/null || true
}

# The digest of every image a live case uses, pulled when absent; recorded per report.
IMAGE_DIGESTS=""
for image in $IMAGES; do
	digest=$(repo_digest "$image")
	if [[ -z "$digest" ]]; then
		log "pulling $image"
		docker pull "$image" >/dev/null || die "docker pull $image failed"
		digest=$(repo_digest "$image")
	fi
	[[ -n "$digest" ]] || die "docker image inspect reports no digest for $image"
	log "image $image digest $digest"
	IMAGE_DIGESTS="${IMAGE_DIGESTS}${image}	${digest}
"
done

# One case: Ensembl VEP once per format, then every candidate

# run_engine <log file> <cmd...>: runs, exit status returned, stderr kept for the report and
# the step's figures left at <log file>.usage.
run_engine() {
	local errlog="$1" STEP_USAGE="$1.usage"
	shift
	run_waited "$errlog.out" "$errlog" "$@"
}
# run_perl <log file>: runs CMD, the Ensembl VEP command. Its first --nearest run on a
# cache slice writes transcript_gene_tss.txt through a temporary file named after the
# process id (AnnotationSource/Cache/Transcript.pm:288, 317), which is the same in every
# container, so two first runs would tear one file: while the slice lacks it, the
# --nearest runs of this invocation hold one lock per slice. flock execs the engine
# (--no-fork) so that the signal its step timer hands on reaches it; without flock
# the run is unlocked, which preflight allows only under --jobs 1 or where no two such
# cases share a slice.
run_perl() {
	local errlog="$1" slice_name lock
	if tss_first_run && command -v flock >/dev/null 2>&1; then
		mkdir -p "$WORK_DIR/locks"
		slice_name="$(basename "$(dirname "$C_PERL_SLICE")")_$(basename "$C_PERL_SLICE")"
		lock="$WORK_DIR/locks/$slice_name.tss.lock"
		run_engine "$errlog" flock --no-fork "$lock" "${CMD[@]}"
	else
		run_engine "$errlog" "${CMD[@]}"
	fi
}
# run_waited <stdout file> <stderr file> <cmd...>: runs the command as a job and waits for
# it, so a signal reaches the trap at once; returns the command's status. Where the caller
# names a file in STEP_USAGE the command runs under the step timer, which leaves the step's
# figures there and passes on the signals it is sent, so an interrupted run still stops the
# engine behind it.
run_waited() {
	local out="$1" err="$2" pid rc
	shift 2
	if [[ -n "${STEP_USAGE:-}" ]]; then
		"$PYTHON" -c "$STEP_USAGE_PROGRAM" "$STEP_USAGE" "$WORK_DIR" "$@" >"$out" 2>"$err" &
	else
		"$@" >"$out" 2>"$err" &
	fi
	pid=$!
	set +e
	wait "$pid"
	rc=$?
	set -e
	return "$rc"
}
error_text() { # error_text <prefix> <rc> <stderr file>
	local tail
	tail=$(tail -c 400 "$3" 2>/dev/null | tr '\n\t' '  ')
	printf '%s exited %s: %s' "$1" "$2" "${tail:-no stderr}"
}
# run_compare <name> <report> <log prefix>: runs CMD, the comparator, and prints the
# write_report argument for the result: format_report=<name>=<report>, or
# format_error=<name>=<text> when the comparator exited 2 or wrote no report.
run_compare() {
	local name="$1" report="$2" logs="$3" rc=0 STEP_USAGE="$3.compare.usage"
	run_waited "$logs.compare.out" "$logs.compare.err" "${CMD[@]}" || rc=$?
	if [[ "$rc" -eq 2 || ! -f "$report" ]]; then
		printf 'format_error=%s=%s' "$name" "$(error_text "comparator" "$rc" "$logs.compare.err")"
	else
		printf 'format_report=%s=%s' "$name" "$report"
	fi
}

# write_report <report path> <key=value ...>: composes report.json and prints the verdict;
# given discovery=<classes asked for>, a discovery pass's record, it writes the case's
# discovery.tsv beside the report, appends the case's lines of the run's table to the file
# discovery_rows=<path> names, and prints the status, then the case's summary line.
write_report() {
	"$PYTHON" - "$@" <<'PY'
import datetime
import gzip
import hashlib
import json
import sys
import zlib

report_path = sys.argv[1]
record = {"formats": {}}
reasons = []
# A discovery pass names the classes it asked for and the file its run-wide rows go to; a
# format that did not run is an error of the pass, not a difference.
discovery_request = None
discovery_rows_path = None
format_errors = []
# The findings of the checks a marked fasta cell asks for, each given as
# index_reason=<format>: fasta index: <finding>: a reason of the case beside its format's
# comparison, never a format that did not run.
index_reasons = []
# What sizing reads: each step's record (usage=<step>=<path>), the comparison program's name
# and version, and per compared format the larger of the two outputs' entry counts.
NOT_MEASURED = "not measured"
usage_files = []
sizing = {"comparison_program": NOT_MEASURED, "comparison_program_version": NOT_MEASURED}
rows_compared = {}
# Per compared format, the entries each requested mask class set aside (the warnings report
# counts its lines), read against the case's masks cell below, and per output format the
# classes its report names among those asked for and leaves out of its unused list.
applied_by_format = {}
used_by_format = {}
args = sys.argv[2:]
i = 0
while i < len(args):
    key, _, value = args[i].partition("=")
    i += 1
    if key == "format_report":
        fmt, _, path = value.partition("=")
        with open(path, encoding="utf-8") as fh:
            sub = json.load(fh)
        sub_reasons = list(sub.get("reasons", []))
        # Two outputs with nothing to compare agree vacuously; an engine that annotated
        # nothing must not pass as one that annotated everything alike.
        counts = sub.get("counts", {})
        entry_counts = [counts.get(side) for side in ("expected_entries", "actual_entries")]
        if sub.get("mode") == "compare" and all(isinstance(count, int) for count in entry_counts):
            rows_compared[fmt] = max(entry_counts)
        if sub.get("mode") == "compare" and counts.get("expected_entries", 0) == 0 and counts.get("actual_entries", 0) == 0:
            sub_reasons.append("no entry compared: both outputs carry no annotation")
        if sub.get("verdict") != "green" and not sub_reasons:
            sub_reasons.append(f"comparator verdict {sub.get('verdict')!r} without a stated reason")
        entry = {"verdict": "red" if sub_reasons else "green", "report": path, "reasons": sub_reasons}
        if "one_sided" in sub:
            entry["one_sided"] = sub["one_sided"]
        if sub.get("floors_not_applicable"):
            entry["floors_not_applicable"] = sub["floors_not_applicable"]
        if sub.get("mode") == "compare":
            masks_report = sub.get("masks", {})
            applied_by_format[fmt] = dict(masks_report.get("applied", {}))
            if "unused" in masks_report:
                used_by_format[fmt] = set(masks_report.get("requested", [])) - set(masks_report["unused"])
        elif sub.get("mode") == "warnings":
            applied_by_format[fmt] = dict(sub.get("warnings", {}).get("masks_applied", {}))
        if fmt in applied_by_format:
            entry["masks_applied"] = applied_by_format[fmt]
        record["formats"][fmt] = entry
        reasons.extend(f"{fmt}: {r}" for r in sub_reasons)
    elif key == "format_error":
        fmt, _, message = value.partition("=")
        record["formats"][fmt] = {"verdict": "red", "error": message}
        reasons.append(f"{fmt}: {message}")
        format_errors.append(f"{fmt}: {message}")
    elif key == "format_text":
        # <format>=<Ensembl VEP output>\t<candidate output>: a --count tool run writes the
        # count alone, so the two files are compared as text; an empty file is red, since
        # the tool prints a count for every input it reads.
        fmt, _, files = value.partition("=")
        texts = {}
        sub_reasons = []
        for side, label, path in zip(("expected", "actual"), ("Ensembl VEP", "candidate"), files.split("\t")):
            # A stored reference may be gzipped, as every stored output may.
            opener = gzip.open if path.endswith(".gz") else open
            try:
                with opener(path, "rt", encoding="utf-8") as fh:
                    texts[side] = fh.read().rstrip("\n")
            except (OSError, ValueError, EOFError, zlib.error) as exc:
                texts[side] = None
                sub_reasons.append(f"{label} output unreadable: {getattr(exc, 'strerror', None) or exc}")
                continue
            if not texts[side].strip():
                sub_reasons.append(f"{label} output is empty where the tool writes its count")
        if not sub_reasons and texts["expected"] != texts["actual"]:
            sub_reasons.append(f"count differs: Ensembl VEP {texts['expected']!r}, candidate {texts['actual']!r}")
        record["formats"][fmt] = {"verdict": "red" if sub_reasons else "green", "compared": "text", "text": texts, "reasons": sub_reasons}
        reasons.extend(f"{fmt}: {r}" for r in sub_reasons)
    elif key == "format_refusal":
        # <format>=<Ensembl VEP exit>\t<its text matched>\t<candidate exit>\t<matched>\t<Ensembl
        # VEP stderr tail>\t<candidate stderr tail>: both engines must exit non-zero with the
        # expected text in their stderr.
        fmt, _, fields = value.partition("=")
        perl_rc, perl_matched, cand_rc, cand_matched, perl_tail, cand_tail = fields.split("\t")
        sub_reasons = []
        entry = {"refusal": {}}
        for engine, label, rc, matched, tail in (
            ("perl", "Ensembl VEP", perl_rc, perl_matched, perl_tail),
            ("candidate", "candidate", cand_rc, cand_matched, cand_tail),
        ):
            entry["refusal"][engine] = {"exit": int(rc), "text_matched": matched == "1"}
            if rc == "0":
                sub_reasons.append(f"{label} exited 0 where a refusal was expected")
            elif matched != "1":
                sub_reasons.append(f"{label} exited {rc} without the expected refusal text: {tail or 'no stderr'}")
        entry["verdict"] = "red" if sub_reasons else "green"
        entry["reasons"] = sub_reasons
        record["formats"][fmt] = entry
        reasons.extend(f"{fmt}: {r}" for r in sub_reasons)
    elif key == "index_reason":
        reasons.append(value)
        index_reasons.append(value)
    elif key == "discovery":
        discovery_request = [] if value == "-" else value.split(",")
    elif key == "discovery_rows":
        discovery_rows_path = value
    elif key == "usage":
        step, _, path = value.partition("=")
        usage_files.append((step, path))
    elif key in ("comparison_program", "comparison_program_version"):
        sizing[key] = value or NOT_MEASURED
    elif key == "flags":
        record[key] = value.split()
    elif key == "masks":
        # The cell as written; its classes and their expected counts are checked once every
        # format is read.
        record[key] = value
    elif key == "unordered":
        record[key] = value == "1"
    elif key == "formats_requested":
        record["formats_requested"] = value.split(",")
    elif key == "reference_provenance":
        record[key] = value or None
        if value:
            try:
                with open(value, encoding="utf-8") as fh:
                    provenance = json.load(fh)
                if not isinstance(provenance, dict):
                    raise ValueError("not a JSON object")
            except (OSError, ValueError) as exc:
                reasons.append(f"stored reference provenance.json unreadable: {exc}")
                provenance = {}
            record["reference_provenance_image"] = provenance.get("vep_image")
    elif key in ("perl_image_digest", "perl_server", "reference_dir", "refusal_dir", "sources_dir", "fasta", "fasta_staged", "manifest", "candidate_dir_cache", "candidate_version", "tool"):
        record[key] = value or None
    else:
        record[key] = value
record["fasta_index_reasons"] = index_reasons
# A stored reference carries its image in its provenance; a live run records the digest
# docker reports.
if not record.get("perl_image_digest") and record.get("reference_provenance_image"):
    record["perl_image_digest"] = record["reference_provenance_image"]
# A floor a format could not count is met only where another output format of the same case
# counted it, one whose report does not set that floor aside too; a case that gives the floor
# no such format leaves it unchecked.
floors_set_aside = {
    fmt: set(entry.get("floors_not_applicable", []))
    for fmt, entry in record["formats"].items()
    if fmt in ("default", "tab", "vcf", "parquet", "json") and "report" in entry
}
for fmt, entry in sorted(record["formats"].items()):
    for name in entry.get("floors_not_applicable", []):
        counted_elsewhere = any(name not in aside for other, aside in floors_set_aside.items() if other != fmt)
        if not counted_elsewhere:
            reasons.append(f"{fmt}: floor {name} is not applicable to {fmt} and no other format of the case counted it")
# A class the cell names must be used in some compared format: it sets something aside
# there, or the format's report leaves it out of its unused list (a class a value set aside
# under another class needs, counted under that other class alone). A class with an expected
# count must apply to that many entries in every output format compared (the warnings report
# counts lines and is not held to it).
masks_cell = record.get("masks", "-")
if discovery_request is None and masks_cell != "-" and applied_by_format:
    for item in masks_cell.split(","):
        klass, _, expected_count = item.partition("=")
        counted = any(applied.get(klass) for applied in applied_by_format.values())
        if not counted and not any(klass in used for used in used_by_format.values()):
            reasons.append(f"mask {klass} is requested and applied in no compared format")
            continue
        if not expected_count:
            continue
        for fmt, applied in sorted(applied_by_format.items()):
            if fmt == "warnings":
                continue
            count = applied.get(klass, 0)
            if count != int(expected_count):
                reasons.append(f"{fmt}: mask {klass} applied to {count} entries, {expected_count} expected")
if discovery_request is None:
    record["verdict"] = "red" if reasons else "green"
    record["reasons"] = reasons
    outcome = [record["verdict"]]
else:
    # A discovery pass: per class the entries it set aside in each compared format, a class
    # an output format's report leaves out of its unused list applied there with nothing
    # counted; the classes in the order asked. No verdict is written: a format's entry keeps
    # its differences as `residual`, the record a `status`.
    named = [] if masks_cell == "-" else [item.partition("=")[0] for item in masks_cell.split(",")]
    compared = [fmt for fmt in record.get("formats_requested", []) if fmt in record["formats"]]
    compared += [fmt for fmt in record["formats"] if fmt not in compared]
    applied = {}
    for fmt in compared:
        counts = applied_by_format.get(fmt, {})
        for klass in sorted(set(counts) | used_by_format.get(fmt, set())):
            if counts.get(klass) or klass in used_by_format.get(fmt, set()):
                applied.setdefault(klass, {})[fmt] = counts.get(klass, 0)
    asked_at = {klass: index for index, klass in enumerate(discovery_request)}
    classes = sorted(applied, key=lambda klass: (asked_at.get(klass, len(asked_at)), klass))
    named_unapplied = [klass for klass in named if klass not in applied]
    applied_unnamed = [klass for klass in classes if klass not in named]
    residual = [reason for reason in reasons if reason not in format_errors]
    for entry in record["formats"].values():
        entry.pop("verdict", None)
        if "reasons" in entry:
            entry["residual"] = entry.pop("reasons")
    record["mode"] = "discover"
    record["status"] = "incomplete" if format_errors else "discovered"
    record["discovery"] = {
        "requested": discovery_request,
        "named": named,
        "applied": {klass: applied[klass] for klass in classes},
        "named_applied_nowhere": named_unapplied,
        "applied_not_named": applied_unnamed,
        "residual": {"count": len(residual), "reasons": residual},
        "errors": format_errors,
    }

    def listed(items):
        return ",".join(str(item) for item in items) or "-"

    if format_errors:
        outcome = [record["status"], "did not run: " + "; ".join(format_errors)]
    else:
        # The case's table beside its report, and its lines of the run's table.
        yes_no = {True: "yes", False: "no"}
        lines = ["format\tclass\trows\tnamed\tapplied"]
        lines += [
            f"{fmt}\t{klass}\t{rows}\t{yes_no[klass in named]}\tyes"
            for klass in classes
            for fmt, rows in applied[klass].items()
        ]
        lines += [f"-\t{klass}\t0\tyes\tno" for klass in named_unapplied]
        with open(report_path.rsplit("/", 1)[0] + "/discovery.tsv", "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        head = f"{record.get('case')}\t{record.get('candidate_pin')}"
        run_lines = [f"{head}\t{klass}\t{listed(applied[klass])}\t{listed(applied[klass].values())}" for klass in classes]
        if discovery_rows_path:
            with open(discovery_rows_path, "a", encoding="utf-8") as fh:
                fh.write("\n".join(run_lines or [f"{head}\t-\t-\t-"]) + "\n")
        outcome = [
            record["status"],
            f"applied: {listed(classes)}; named and applied nowhere: {listed(named_unapplied)}; "
            f"applied and not named: {listed(applied_unnamed)}; residual differences: {len(residual)}",
        ]
# The sizing block: a step whose record is absent or unreadable is left out, and a figure no
# step gives reads `not measured`.
steps = {}
disk_starts = []
disk_peaks = []
for step, path in usage_files:
    try:
        with open(path, encoding="utf-8") as fh:
            usage = json.load(fh)
    except (OSError, ValueError):
        continue
    if not isinstance(usage, dict):
        continue
    peak_rss = usage.get("peak_rss_kib")
    # The Ensembl VEP run is a container: the timed process is its client.
    if step.startswith("ensembl-vep.") or not isinstance(peak_rss, int):
        peak_rss = NOT_MEASURED
    wall_seconds = usage.get("wall_seconds")
    steps[step] = {
        "wall_seconds": wall_seconds if isinstance(wall_seconds, (int, float)) else NOT_MEASURED,
        "peak_rss_kib": peak_rss,
        "exit": usage.get("exit"),
    }
    if isinstance(usage.get("disk_used_start_bytes"), int):
        disk_starts.append(usage["disk_used_start_bytes"])
    if isinstance(usage.get("disk_used_peak_bytes"), int):
        disk_peaks.append(usage["disk_used_peak_bytes"])
peaks = [entry["peak_rss_kib"] for entry in steps.values() if isinstance(entry["peak_rss_kib"], int)]
try:
    stamps = [
        datetime.datetime.strptime(record[key], "%Y-%m-%dT%H:%M:%SZ") for key in ("started_utc", "finished_utc")
    ]
    wall_minutes = round((stamps[1] - stamps[0]).total_seconds() / 60, 2)
except (KeyError, TypeError, ValueError):
    wall_minutes = NOT_MEASURED
sizing.update(
    {
        "rows_compared": rows_compared or NOT_MEASURED,
        "peak_memory_kib": max(peaks) if peaks else NOT_MEASURED,
        "work_disk_peak_bytes": max(disk_peaks) if disk_peaks else NOT_MEASURED,
        "work_disk_growth_bytes": max(max(disk_peaks) - min(disk_starts), 0) if disk_peaks and disk_starts else NOT_MEASURED,
        "wall_minutes": wall_minutes,
        "steps": steps,
    }
)
record["sizing"] = sizing
with open(report_path, "w", encoding="utf-8") as fh:
    json.dump(record, fh, indent=1, sort_keys=True)
    fh.write("\n")
print("\n".join(outcome))
PY
}

run_case() {
	load_case "$1"
	case "$C_RUN_SOURCE" in
	table) log "$C_ID: compares under the run $C_RUN of its manifest, named by $CASE_RUNS: $C_RUN_REASON" ;;
	pooled) log "$C_ID: compares under no named run of its manifest, pooled by $CASE_RUNS: $C_RUN_REASON" ;;
	esac
	local case_work fasta engine_fasta perl_dir f rc started
	case_work=$(case_work_dir)
	local perl_image_digest="" perl_errors=""
	# A staged directory that an interrupted run of a read-only case left in mode 0555 cannot
	# be emptied, so every staged directory gets its write bit before the tree is removed.
	[[ ! -d "$case_work/fasta" ]] || chmod -R u+w "$case_work/fasta"
	rm -rf "$case_work"
	mkdir -p "$case_work/perl"
	: >"$case_work/status.rows"
	# What the case's report is handed beside its formats' results: the steps whose figures
	# the sizing block reads and, in a discovery pass, the classes asked for.
	local usage_args=()
	local discovery_args=()
	[[ "$DISCOVER" != "true" ]] || discovery_args=("discovery=$C_REQUEST_CLASSES" "discovery_rows=$case_work/discovery.rows")
	fasta=$(case_fasta)
	started=$(date -u +%Y-%m-%dT%H:%M:%SZ)

	# stage_fasta <directory>: gives one engine run of a case whose fasta cell is marked its
	# copy of the source, alone in the directory under the source's file name, in mode 0555 for
	# the read-only form. It is nested in run_case, as format_result is.
	stage_fasta() {
		mkdir -p "$1"
		cp "$C_FASTA_PATH" "$1/"
		[[ "$C_FASTA_STAGED" != "readonly" ]] || chmod 0555 "$1"
	}
	# fasta_index_findings <format> <Ensembl VEP's staged directory> <the candidate's> <Ensembl
	# VEP's exit status> <the candidate's>: the checks of a marked fasta cell (the header's
	# fasta column) over one output format, one finding a line as `<format>: fasta index:
	# <finding>`, which names an engine and a file name and never a path under the work tree.
	# It is nested in run_case, as format_result is.
	fasta_index_findings() {
		local f="$1" perl_staged="$2" cand_staged="$3" perl_rc="$4" cand_rc="$5"
		local name="${C_FASTA_PATH##*/}" engine dir entry index indexes="fai" head
		if [[ "$C_FASTA_STAGED" == "readonly" ]]; then
			# Each directory holds the copy alone where no engine wrote beside it.
			for engine in "Ensembl VEP:$perl_staged" "the candidate:$cand_staged"; do
				dir="${engine#*:}"
				for entry in "$dir"/* "$dir"/.*; do
					[[ -e "$entry" || -L "$entry" ]] || continue
					case "${entry##*/}" in
					. | .. | "$name") ;;
					*) printf '%s: fasta index: %s wrote %s into a directory it was given read-only\n' "$f" "${engine%%:*}" "${entry##*/}" ;;
					esac
				done
			done
			return 0
		fi
		# An engine whose run of the format failed is not asked for an index.
		[[ "$perl_rc" == "0" && "$cand_rc" == "0" ]] || return 0
		# The first 18 bytes tell a BGZF file, which is indexed by a .gzi beside its .fai, from
		# any other file, as compare_cmd tells one in a cache directory.
		head=$(head -c 18 "$C_FASTA_PATH" | od -An -v -tx1 | tr -s ' \n' ' ')
		head=" ${head# }"
		head="${head% } "
		case "$head" in
		" 1f 8b 08 "?[4567cdef]" "??" "??" "??" "??" "??" "??" 06 00 42 43 02 00 "??" "??" ") indexes="fai gzi" ;;
		esac
		for engine in "Ensembl VEP:$perl_staged" "the candidate:$cand_staged"; do
			dir="${engine#*:}"
			engine="${engine%%:*}"
			[[ -f "$dir/$name.fai" ]] || printf '%s: fasta index: %s wrote no %s.fai\n' "$f" "$engine" "$name"
			[[ "$indexes" == "fai" || -f "$dir/$name.gzi" ]] || printf '%s: fasta index: %s wrote no %s.gzi beside a BGZF file\n' "$f" "$engine" "$name"
			# Ensembl VEP writes <name>.vep.lock before it indexes (FastaSequence.pm:334) and
			# removes it once the index stands (FastaSequence.pm:356).
			[[ ! -e "$dir/$name.vep.lock" ]] || printf '%s: fasta index: %s left %s.vep.lock\n' "$f" "$engine" "$name"
			for entry in "$dir"/* "$dir"/.*; do
				case "$entry" in
				*.tmp) [[ ! -e "$entry" || "${entry##*/}" == "$name" ]] || printf '%s: fasta index: %s left %s\n' "$f" "$engine" "${entry##*/}" ;;
				esac
			done
		done
		for index in $indexes; do
			if [[ -f "$perl_staged/$name.$index" && -f "$cand_staged/$name.$index" ]] && ! cmp -s "$perl_staged/$name.$index" "$cand_staged/$name.$index"; then
				printf '%s: fasta index: %s.%s differs between Ensembl VEP and the candidate\n' "$f" "$name" "$index"
			fi
		done
	}

	# Every live Ensembl VEP run keeps its exit status per format, which the checks of a marked
	# fasta cell read; a refusal case also keeps whether the stderr carried perl.txt's text, and
	# the stderr's tail, for the candidate's formats to be judged against.
	local reference_provenance="" perl_text="" cand_text="" perl_rcs="" perl_matches="" perl_tails="" matched
	if [[ -n "$C_REFERENCE_DIR" ]]; then
		perl_dir="$C_REFERENCE_DIR"
		# A stored directory's provenance.json sits beside it; a stored filtered file's is the
		# corpus's, beside the input.
		if [[ -n "$C_REFERENCE_FILE" ]]; then
			reference_provenance="$C_CORPUS_DIR/provenance.json"
		else
			reference_provenance="$(dirname "$C_REFERENCE_DIR")/provenance.json"
		fi
		[[ -f "$reference_provenance" ]] || reference_provenance=""
		log "[$C_ID] Ensembl VEP: stored reference $perl_dir"
	else
		perl_dir="$case_work/perl"
		perl_image_digest=$(tsv_lookup "$IMAGE_DIGESTS" "$C_PERL_IMAGE")
		if [[ -n "$C_REFUSAL_DIR" ]]; then
			perl_text=$(<"$C_REFUSAL_DIR/perl.txt")
			cand_text=$(<"$C_REFUSAL_DIR/candidate.txt")
		fi
		for f in $C_PERL_FORMATS; do
			engine_fasta="$fasta"
			if [[ -n "$C_FASTA_STAGED" ]]; then
				engine_fasta="$case_work/fasta/perl/$f/${C_FASTA_PATH##*/}"
				stage_fasta "$case_work/fasta/perl/$f"
			fi
			perl_cmd "$f" "$perl_dir" "$engine_fasta"
			log "[$C_ID] Ensembl VEP $f ($C_PERL_MODE${C_PERL_SERVER:+ $C_PERL_SERVER})"
			rc=0
			run_perl "$case_work/perl.$f.stderr" || rc=$?
			# The read-only form's directory gets its write bit back as soon as its run has ended,
			# so no removal of the work tree meets a directory it cannot empty.
			[[ "$C_FASTA_STAGED" != "readonly" ]] || chmod u+w "$case_work/fasta/perl/$f"
			perl_rcs="${perl_rcs}${f}	${rc}
"
			if [[ -n "$C_REFUSAL_DIR" ]]; then
				matched=0
				if [[ "$rc" -ne 0 && "$(<"$case_work/perl.$f.stderr")" == *"$perl_text"* ]]; then
					matched=1
				fi
				perl_matches="${perl_matches}${f}	${matched}
"
				perl_tails="${perl_tails}${f}	$(tail -c 400 "$case_work/perl.$f.stderr" 2>/dev/null | tr '\n\t' '  ')
"
				log "[$C_ID] Ensembl VEP $f exited $rc; refusal text matched: $matched"
			elif [[ "$rc" -ne 0 ]]; then
				perl_errors="${perl_errors}${f}	$(error_text "Ensembl VEP" "$rc" "$case_work/perl.$f.stderr")
"
				echo "WARN: [feature_oracle] $C_ID: Ensembl VEP $f exited $rc" >&2
			elif [[ ! -e "$perl_dir/$(format_file "$f")" ]]; then
				perl_errors="${perl_errors}${f}	Ensembl VEP exited 0 and wrote no $(format_file "$f")
"
			fi
		done
	fi

	# format_result <candidate> <format> <candidate dir> <report dir> <fasta or empty>: the
	# candidate's run of one output format of the loaded case, then its comparison; prints the
	# write_report argument for the result (format_refusal, format_error, format_text or
	# run_compare's line). Reads run_case's locals (case_work, pin, perl_errors, perl_rcs and,
	# for a refusal case, cand_text, perl_matches and perl_tails) through dynamic scoping, and
	# is defined inside run_case so that the script's top-level functions stay its interface.
	# For a case whose fasta cell is marked the run reads a staged copy of its own, and the
	# findings of the format's checks go to <candidate dir>/<format>.index, one a line.
	format_result() {
		local bin="$1" f="$2" cand_dir="$3" report_dir="$4" fasta="$5"
		local rc=0 pf perl_error expected actual matched perl_rc perl_matched perl_tail cand_tail
		local staged_dir="$case_work/fasta/candidates/$pin/$f"
		if [[ -n "$C_FASTA_STAGED" ]]; then
			fasta="$staged_dir/${C_FASTA_PATH##*/}"
			stage_fasta "$staged_dir"
		fi
		candidate_cmd "$bin" "$f" "$cand_dir" "$fasta"
		run_engine "$cand_dir/$f.stderr" "${CMD[@]}" || rc=$?
		[[ "$C_FASTA_STAGED" != "readonly" ]] || chmod u+w "$staged_dir"
		# The run's exit status, for the warnings and stats verdicts, which read the default
		# run's: a run that fails after opening its output leaves the files.
		printf '%s' "$rc" >"$cand_dir/$f.rc"
		pf=$(perl_format "$f")
		if [[ -n "$C_FASTA_STAGED" ]]; then
			fasta_index_findings "$f" "$case_work/fasta/perl/$pf" "$staged_dir" "$(tsv_lookup "$perl_rcs" "$pf")" "$rc" >"$cand_dir/$f.index"
		fi
		if [[ -n "$C_REFUSAL_DIR" ]]; then
			matched=0
			if [[ "$rc" -ne 0 && "$(<"$cand_dir/$f.stderr")" == *"$cand_text"* ]]; then
				matched=1
			fi
			perl_rc=$(tsv_lookup "$perl_rcs" "$pf")
			perl_matched=$(tsv_lookup "$perl_matches" "$pf")
			perl_tail=$(tsv_lookup "$perl_tails" "$pf")
			cand_tail=$(tail -c 400 "$cand_dir/$f.stderr" 2>/dev/null | tr '\n\t' '  ')
			printf 'format_refusal=%s=%s\t%s\t%s\t%s\t%s\t%s' \
				"$f" "$perl_rc" "$perl_matched" "$rc" "$matched" "$perl_tail" "$cand_tail"
			return
		fi
		if [[ "$rc" -ne 0 ]]; then
			printf 'format_error=%s=%s' "$f" "$(error_text "candidate" "$rc" "$cand_dir/$f.stderr")"
			return
		fi
		perl_error=$(tsv_lookup "$perl_errors" "$pf")
		if [[ -n "$perl_error" ]]; then
			printf 'format_error=%s=%s' "$f" "$perl_error"
			return
		fi
		expected=$(expected_output "$pf")
		actual="$cand_dir/$(format_file "$f")"
		if [[ "$C_TOOL_COUNT" == "true" ]]; then
			printf 'format_text=%s=%s\t%s' "$f" "$expected" "$actual"
			return
		fi
		compare_cmd "$f" "$expected" "$actual" "$report_dir/$f.json"
		run_compare "$f" "$report_dir/$f.json" "$cand_dir/$f"
	}
	# An interrupted format job stops the engine or comparator it waits on.
	on_format_job_signal() {
		trap - INT TERM
		local pid
		for pid in $(jobs -p); do
			kill -TERM "$pid" 2>/dev/null || true
		done
		exit 130
	}
	local i pin bin cand_dir report_dir perl_error expected actual report_args job_tail verdict digest cand_default_rc
	local index_reason
	local warned_file cand_warned_rc
	for ((i = 0; i < N_CANDIDATES; i++)); do
		pin="${PINS[$i]}"
		bin="${CANDIDATES[$i]}"
		cand_dir="$case_work/candidates/$pin"
		report_dir="$OUT_DIR/$pin/$C_ID"
		rm -rf "$report_dir"
		mkdir -p "$cand_dir" "$report_dir"
		# A format's result lands in <cand_dir>/<format>.result. The format's start is logged
		# by the run once the job has a slot, not inside the job, so the log shows the formats
		# in flight while --format-jobs runs them; what a job itself writes goes to
		# <format>.job.out and .job.err beside its outputs and is replayed in format order once
		# every job has ended.
		for f in $C_OUTPUT_FORMATS; do
			if [[ "$FORMAT_JOBS" -eq 1 ]]; then
				log "[$C_ID] $pin $f"
				format_result "$bin" "$f" "$cand_dir" "$report_dir" "$fasta" >"$cand_dir/$f.result"
			else
				while [[ "$(running_jobs)" -ge "$FORMAT_JOBS" ]]; do sleep 1; done
				log "[$C_ID] $pin $f"
				(
					trap on_format_job_signal INT TERM
					format_result "$bin" "$f" "$cand_dir" "$report_dir" "$fasta" >"$cand_dir/$f.result"
				) >"$cand_dir/$f.job.out" 2>"$cand_dir/$f.job.err" &
			fi
		done
		if [[ "$FORMAT_JOBS" -gt 1 ]]; then
			wait
			for f in $C_OUTPUT_FORMATS; do
				cat "$cand_dir/$f.job.out"
				cat "$cand_dir/$f.job.err" >&2
			done
		fi
		# A job that died leaves its result file empty (the redirection opened it first); its
		# stderr's tail travels in the reason, since the candidate directory goes with the case.
		report_args=()
		for f in $C_OUTPUT_FORMATS; do
			if [[ -s "$cand_dir/$f.result" ]]; then
				report_args+=("$(<"$cand_dir/$f.result")")
			else
				job_tail=$(tail -c 400 "$cand_dir/$f.job.err" 2>/dev/null | tr '\n\t' '  ')
				report_args+=("format_error=$f=the format's job ended without a result: ${job_tail:-no stderr}")
			fi
			# The findings of a marked fasta cell's checks over the format, each a reason of its
			# own beside the format's result.
			if [[ -s "$cand_dir/$f.index" ]]; then
				while IFS= read -r index_reason; do
					report_args+=("index_reason=$index_reason")
				done <"$cand_dir/$f.index"
			fi
		done
		# A format job that died before its engine ended leaves no status: the run did not complete.
		cand_default_rc=0
		[[ ! -s "$cand_dir/default.rc" ]] || cand_default_rc=$(<"$cand_dir/default.rc")
		if [[ "$C_WARNINGS" == "true" ]]; then
			# The run whose warnings the case compares, and its exit status.
			warned_file=$(format_file "$C_WARNINGS_FORMAT")
			cand_warned_rc=0
			[[ ! -s "$cand_dir/$C_WARNINGS_FORMAT.rc" ]] || cand_warned_rc=$(<"$cand_dir/$C_WARNINGS_FORMAT.rc")
			perl_error=$(tsv_lookup "$perl_errors" "$C_WARNINGS_FORMAT")
			if [[ -n "$perl_error" ]]; then
				report_args+=("format_error=warnings=$perl_error")
			elif [[ "$cand_warned_rc" -ne 0 ]]; then
				report_args+=("format_error=warnings=the candidate's $C_WARNINGS_FORMAT run exited $cand_warned_rc")
			elif [[ ! -e "$cand_dir/$warned_file" ]]; then
				report_args+=("format_error=warnings=the candidate's $C_WARNINGS_FORMAT run did not complete")
			else
				if [[ -n "$C_REFERENCE_DIR" ]]; then
					expected=$(stored_warnings_file)
				else
					expected="$perl_dir/${warned_file}_warnings.txt"
					[[ -e "$expected" ]] || : >"$expected"
				fi
				actual="$cand_dir/${warned_file}_warnings.txt"
				[[ -e "$actual" ]] || : >"$actual"
				warnings_cmd "$expected" "$actual" "$report_dir/warnings.json"
				report_args+=("$(run_compare warnings "$report_dir/warnings.json" "$cand_dir/warnings")")
			fi
		fi
		# Unlike a warnings file, whose absence states no warnings, a statistics report is
		# written by every default run under --stats_text, so a missing one is a failure.
		if [[ "$C_STATS" == "true" ]]; then
			perl_error=$(tsv_lookup "$perl_errors" default)
			if [[ -n "$C_REFERENCE_DIR" ]]; then
				expected=$(stored_stats_file)
			else
				expected=$(stats_report_file "$perl_dir")
			fi
			actual=$(stats_report_file "$cand_dir")
			if [[ -n "$perl_error" ]]; then
				report_args+=("format_error=stats=$perl_error")
			elif [[ ! -e "$expected" ]]; then
				report_args+=("format_error=stats=Ensembl VEP exited 0 and wrote no default.txt_summary.txt")
			elif [[ "$cand_default_rc" -ne 0 ]]; then
				report_args+=("format_error=stats=the candidate's default run exited $cand_default_rc")
			elif [[ ! -e "$cand_dir/default.txt" ]]; then
				report_args+=("format_error=stats=the candidate's default run did not complete")
			elif [[ ! -e "$actual" ]]; then
				report_args+=("format_error=stats=the candidate exited 0 and wrote no default.txt_summary.txt")
			else
				stats_cmd "$expected" "$actual" "$report_dir/stats.json"
				report_args+=("$(run_compare stats "$report_dir/stats.json" "$cand_dir/stats")")
			fi
		fi
		# The steps whose figures the report's sizing block reads: the live Ensembl VEP runs,
		# shared by every candidate, then this candidate's runs and comparisons.
		usage_args=("comparison_program=$COMPARISON_PROGRAM" "comparison_program_version=$COMPARISON_PROGRAM_VERSION")
		if [[ -z "$C_REFERENCE_DIR" ]]; then
			for f in $C_PERL_FORMATS; do
				usage_args+=("usage=ensembl-vep.$f=$case_work/perl.$f.stderr.usage")
			done
		fi
		for f in $C_OUTPUT_FORMATS; do
			usage_args+=("usage=candidate.$f=$cand_dir/$f.stderr.usage" "usage=comparison.$f=$cand_dir/$f.compare.usage")
		done
		[[ "$C_WARNINGS" != "true" ]] || usage_args+=("usage=comparison.warnings=$cand_dir/warnings.compare.usage")
		[[ "$C_STATS" != "true" ]] || usage_args+=("usage=comparison.stats=$cand_dir/stats.compare.usage")
		# A marked fasta cell is recorded by its source and its marker: the staged copies go with
		# the case's outputs.
		verdict=$(write_report "$report_dir/report.json" \
			"case=$C_ID" "tier=$C_TIER" "input=$C_INPUT_PATH" "cache=$C_CACHE" "species=$C_SPECIES" \
			"candidate_json_cache=$C_JSON_CACHE" "candidate_dir_cache=$C_DIR_CACHE" "fasta=${fasta:-$C_FASTA_PATH}" \
			"fasta_staged=${C_FASTA_STAGED:+${C_FASTA%%:*}}" "flags=$C_FLAGS" \
			"formats_requested=$C_FORMATS" "manifest=$C_MANIFEST" "masks=$C_MASKS" \
			"floors=$C_FLOORS" "field_pairs=$C_FIELD_PAIRS" "unordered=$C_UNORDERED" \
			"perl_image=$C_PERL_IMAGE" "perl_image_digest=$perl_image_digest" "perl_mode=$C_PERL_MODE" \
			"perl_server=$C_PERL_SERVER" "reference_release=$C_REFERENCE_RELEASE" \
			"reference_dir=$C_REFERENCE_DIR" "reference_provenance=$reference_provenance" \
			"refusal_dir=$C_REFUSAL_DIR" "sources_dir=$SOURCES_DIR" "tool=${C_TOOL:+$C_TOOL_CELL}" \
			"candidate=$bin" "candidate_pin=$pin" "candidate_version=$("$bin" --version 2>/dev/null | head -1 || true)" \
			"run_utc=$RUN_UTC" "started_utc=$started" "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
			${discovery_args[@]+"${discovery_args[@]}"} "${usage_args[@]}" "${report_args[@]}")
		digest=$(sha256_of "$report_dir/report.json")
		# A discovery pass's record prints two lines: the status, then the case's summary.
		printf '%s\t%s\t%s\t%s\t%s\n' "$C_ID" "$pin" "${verdict%%$'\n'*}" "$digest" "$RUN_UTC" >>"$case_work/status.rows"
		if [[ "$DISCOVER" == "true" ]]; then
			printf '%s\t%s\n' "$pin" "${verdict#*$'\n'}" >>"$case_work/discovery.summary"
			log "[$C_ID] $pin: ${verdict%%$'\n'*}; ${verdict#*$'\n'}"
		else
			log "[$C_ID] $pin: $verdict"
		fi
		[[ "$KEEP_OUTPUTS" == "true" ]] || rm -rf "$cand_dir"
	done
	# A format job that dies while its engine runs leaves its read-only directory in mode 0555.
	[[ ! -d "$case_work/fasta" ]] || chmod -R u+w "$case_work/fasta"
	[[ "$KEEP_OUTPUTS" == "true" ]] || rm -rf "$case_work/perl" "$case_work/fasta"
}

# All cases, --jobs at a time
running_jobs() {
	local n
	n=$(jobs -rp | grep -c .) || n=0
	echo "$n"
}
CASE_NAMES=""
while IFS= read -r line; do
	[[ -n "$line" ]] || continue
	cid=$(cell "$line" case)
	CASE_NAMES="$CASE_NAMES $cid"
	if [[ "$JOBS" -eq 1 ]]; then
		run_case "$line"
	else
		while [[ "$(running_jobs)" -ge "$JOBS" ]]; do sleep 1; done
		mkdir -p "$WORK_DIR/logs"
		(run_case "$line") >"$WORK_DIR/logs/$cid.log" 2>&1 &
	fi
done <<<"$SELECTED"
if [[ "$JOBS" -gt 1 ]]; then
	wait
	for cid in $CASE_NAMES; do
		cat "$WORK_DIR/logs/$cid.log"
	done
fi

# status.tsv: one block per candidate, the cases in table order; a case that left no
# row for a candidate (the runner stopped inside it) is red. A discovery pass writes the
# same table with its own two words, discovered and incomplete, counts the rows that are
# not `discovered` in N_RED, and writes the run's discovery.tsv in the same order.
STATUS="$OUT_DIR/status.tsv"
printf 'case\tcandidate_pin\tverdict\treport_digest\trun_utc\n' >"$STATUS.tmp"
DONE_WORD="green"
UNDONE_WORD="red"
ROW_NOUN="verdict"
if [[ "$DISCOVER" == "true" ]]; then
	DONE_WORD="discovered"
	UNDONE_WORD="incomplete"
	ROW_NOUN="status"
	printf 'case\tcandidate_pin\tclass\tapplied_formats\trows\n' >"$OUT_DIR/discovery.tsv.tmp"
fi
N_RED=0
N_ROWS=0
for ((i = 0; i < N_CANDIDATES; i++)); do
	pin="${PINS[$i]}"
	for cid in $CASE_NAMES; do
		# The pin is compared as text (`p ""`), here and below: two pins that are numerals of
		# one number (1, 1.0) are equal to awk as numbers.
		row=$(awk -F'\t' -v p="$pin" '$2 == p ""' "$WORK_DIR/cases/$cid/status.rows" 2>/dev/null || true)
		if [[ -z "$row" ]]; then
			row=$(printf '%s\t%s\t%s\t-\t%s' "$cid" "$pin" "$UNDONE_WORD" "$RUN_UTC")
			echo "ERROR: [feature_oracle] $cid: no $ROW_NOUN for $pin; the case did not complete" >&2
		fi
		printf '%s\n' "$row" >>"$STATUS.tmp"
		N_ROWS=$((N_ROWS + 1))
		[[ "$row" == "$cid	$pin	$DONE_WORD	"* ]] || N_RED=$((N_RED + 1))
		[[ "$DISCOVER" != "true" ]] || awk -F'\t' -v p="$pin" '$2 == p ""' "$WORK_DIR/cases/$cid/discovery.rows" >>"$OUT_DIR/discovery.tsv.tmp" 2>/dev/null || true
	done
done
mv -f "$STATUS.tmp" "$STATUS"
[[ "$DISCOVER" != "true" ]] || mv -f "$OUT_DIR/discovery.tsv.tmp" "$OUT_DIR/discovery.tsv"

"$PYTHON" - "$OUT_DIR/run.json" "$RUN_UTC" "$TIER" "$CASES" "$REFERENCE_DIR" "$SOURCES_DIR" "$PERL_FORK" "$FORK" "$IMAGE_DIGESTS" "$N_ROWS" "$N_RED" "$DISCOVER" "$DISCOVER_CLASSES" "$COMPARISON_PROGRAM" "$COMPARISON_PROGRAM_VERSION" "${PINS[@]}" <<'PY'
import csv
import json
import os
import sys

(path, run_utc, tier, cases, reference_dir, sources_dir, perl_fork, fork, digests, rows, red, discover, classes) = sys.argv[1:14]
program, program_version = sys.argv[14:16]
record = {
    "run_utc": run_utc,
    "tier": tier,
    "case_table": cases,
    "reference_dir": reference_dir or None,
    "sources_dir": sources_dir or None,
    "perl_fork": int(perl_fork),
    "candidate_fork": int(fork) if fork else None,
    "image_digests": dict(line.split("\t", 1) for line in digests.splitlines() if line),
    "candidate_pins": sys.argv[16:],
    "status_rows": int(rows),
    "comparison_program": program,
    "comparison_program_version": program_version,
}
# Per status row, the sizing block of the case's report; a case that left no report, or a
# report without the block, reads `not measured`.
out_dir = os.path.dirname(path)
record["case_sizes"] = []
with open(os.path.join(out_dir, "status.tsv"), encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        sizing = None
        # A row without a digest left no report: a file at its path is an earlier run's.
        if row["report_digest"] != "-":
            try:
                with open(os.path.join(out_dir, row["candidate_pin"], row["case"], "report.json"), encoding="utf-8") as report:
                    sizing = json.load(report).get("sizing")
            except (OSError, ValueError, AttributeError):
                sizing = None
        record["case_sizes"].append(
            {
                "case": row["case"],
                "candidate_pin": row["candidate_pin"],
                "sizing": sizing if isinstance(sizing, dict) else "not measured",
            }
        )
# A discovery pass counts no red: its rows are discovered or incomplete.
if discover == "true":
    record["mode"] = "discover"
    record["discover_classes"] = classes.split(",") if classes else None
    record["discovered"] = int(rows) - int(red)
    record["incomplete"] = int(red)
else:
    record["mode"] = "verdict"
    record["red"] = int(red)
with open(path, "w", encoding="utf-8") as fh:
    json.dump(record, fh, indent=1, sort_keys=True)
    fh.write("\n")
PY

if [[ "$DISCOVER" == "true" ]]; then
	# One summary line per case and candidate, in the tables' order.
	for ((i = 0; i < N_CANDIDATES; i++)); do
		pin="${PINS[$i]}"
		for cid in $CASE_NAMES; do
			summary=$(awk -F'\t' -v p="$pin" '$1 == p "" {print $2}' "$WORK_DIR/cases/$cid/discovery.summary" 2>/dev/null || true)
			log "discovery $cid $pin: ${summary:-did not run: the case did not complete}"
		done
	done
	log "feature_oracle discovery complete: $((N_ROWS - N_RED)) case(s) discovered, $N_RED incomplete, no verdict given; $OUT_DIR/discovery.tsv"
else
	log "feature_oracle complete: $N_ROWS verdict(s), $N_RED red; $STATUS"
fi
[[ "$N_RED" -eq 0 ]] || exit 1
