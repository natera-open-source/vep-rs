# Changelog

All notable changes to vep-rs, a Rust reimplementation of Ensembl VEP, are recorded in this file, one section per release.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Dependencies and tooling
  - The Bioconda recipe `bioconda/vep-rs` names 0.6.0 and drops its `duckdb-cli` run requirement: 0.6.0 writes Parquet without the `duckdb` program. (#PR)

## [0.6.0] - 2026-11-19

### Added

- Caches, reference FASTA and the cache tools
  - A transcript's codon table selects any of the 18 genetic codes BioPerl's `Bio::Tools::CodonTable` ships, each with its start and stop codons; a table it lacks is refused. (#115)
  - `vep-cache-builder --gtf` or `--gff3` with `--genome-fasta` builds a cache from any GTF or GFF3, read as Ensembl VEP reads a `--gtf` or `--gff` source, no protein FASTA needed. (#115)
  - A converted RefSeq or merged cache switches `--use_transcript_ref` on with Ensembl VEP's notice, writes `BAM_EDIT`, reads each edited model's cDNA and `source_cache`; `--use_given_ref` overrides. (#115)
  - `storable_to_json.pl` reads an edited model's mapper pairs from the Storable mapper object and writes `source_type` into `info.json`; a declared set the flags do not select is refused. (#115)
  - A cache of another species lists its chromosomes in its karyotype's order, keeps its own names at run time, and gets an assembly-hint warning only when a plugin loads. (#115)
  - A cache built from the genome alone keeps the `artifact` and readthrough transcripts Ensembl VEP's `--gff` source annotates, which a build with the release's protein FASTA or core dumps leaves out, and a GTF transcript line's `ccds_id` is read as its CCDS identifier, as a GFF3's `ccdsid` is. (#115)
  - `vep-cache-builder --auto <species> <assembly> --release <N>` fetches a release's inputs, checks each against Ensembl's `CHECKSUMS`, builds the cache in Ensembl VEP's layout and places the indexed FASTA inside it. (#115)
  - On an edited RefSeq model, an allele wholly inside bases the model's alignment deleted takes the region's term, the warning naming the model once. (#115)
  - An edited RefSeq model's HGVSc spells a shifted insertion's bases as the shift rotates them, as Ensembl VEP does on a coding row under `--shift_3prime 1` alone. (#115)
  - `storable_to_json.pl` and `vep-cache-converter` write the cache's regulatory tree as `regulatory/<chr>/<start>-<end>.json` shards and set the `regulatory` marker of `info.json` when they do. (#78)
  - `scripts/data/storable_to_json.pl` writes regulatory features as activity codes and motifs naming one `regulatory/binding_matrices.json` table, and a tree whose `info.json` lacks `regulatory_schema` is refused. (#115)
  - `--transcript_filter` and `--exclude_predicted` on a cache built from a GTF or GFF3 keep the transcripts Ensembl VEP's `--gtf` and `--gff` sources keep; the keys a file lacks are documented. (#115)
  - `vep-cache-builder --regulatory-gff`, `--emars-gff` and `--motif-gff` bake those files into a cache's regulatory tree, each motif's matrix a stub without counts, and warn on a non-integer promoter bound or a repeated `ID` as `--regulatory_gff` does. (#115)
- Transcript attributes, predictions and consequences
  - `--regulatory` (and `--everything`) writes the `RegulatoryFeature` and `MotifFeature` rows of a cache converted with its regulatory tree. (#78, #61)
  - The regulatory and motif rows carry their terms, `BIOTYPE`, `STRAND`, `MOTIF_NAME`, `MOTIF_POS` and `TRANSCRIPTION_FACTORS`, and an allele overlapping no transcript keeps its intergenic companion row beside them. (#78)
  - `--regulatory` writes `HIGH_INF_POS` and `MOTIF_SCORE_CHANGE` on motif rows from the cache's binding matrices, as Ensembl VEP computes them, and orders motif rows by database id. (#115)
  - Regulatory and motif rows follow Ensembl VEP under `--per_gene`, the pick flags, plugin feature types, structural classes, breakend positions and `--shift_3prime`, `SHIFT_LENGTH` written on transcript rows alone; plugins run on the row kinds they declare. (#115)
  - A structural variant's regulatory and motif rows come in the cache's order, or file by file in the `--regulatory_gff` flag's source order and then in file order; the pick family ranks them in that written order, and a pick among tied rows keeps the first listed. (#115)
  - `--cell_type <list>` prints `CELL_TYPE`, each listed epigenome's activity on every regulatory-feature row, read against the cache's epigenome names, which may contain commas; an unknown name is refused. (#115)
  - `--regulatory_gff file=,motifs=,emars=` reads the Ensembl regulatory GFF files in place of the cache's regulatory data, as Ensembl VEP does; `--extended_promoters` widens their promoters. (#115)
  - `--regulatory_gff` refuses a record whose start is 0 or whose end is below its start, naming the record, and a skipped-type line whose start or end no number leads or whose end is below zero, where Ensembl VEP loses the records after such an end. (#115)
  - A `--regulatory_gff` file refused as it is opened, or a refused `--cell_type` name, ends the run before the output file and the statistics report exist, and `--show_cache_info` prints the `regulatory_gff` entry with the cache's source versions and refuses a missing GFF before printing, as Ensembl VEP does. (#115)
  - Two `--regulatory_gff` records of one file sharing an `ID` give a variant inside both one row, the later record's, with Ensembl VEP's duplicate-identifier warning once per identifier and file, and `--extended_promoters` ignores an `extended_start` or `extended_end` that is not a whole number with one warning naming the record and the attribute. (#115)
- Output formats, columns and headers
  - `--refseq` writes `REFSEQ_MATCH` and `REFSEQ_OFFSET`, and `--merged` writes `REFSEQ_MATCH`, `SOURCE` and `REFSEQ_OFFSET`, in every format, and the `Transcript-assembly mismatch` warning is written. (#115)
  - `--parquet_partition chrom` lays the Parquet parts out in one `chrom=<name>/` directory per chromosome as spelled, `chrom` still a column. (#115)
  - `--stats_text` writes Ensembl VEP's run statistics report, `<output>_summary.txt`, with its sections and charts; `--no_stats` skips it, `--stats_file` names it, `--stats_html` is accepted. (#115)
  - The statistics report is written by default as a self-contained HTML page, `<output>_summary.html`, with Ensembl VEP's sections, element ids, titles and tables, every chart drawn as inline SVG. (#115)
  - The statistics report counts a `--minimal` split feature writing no row that `--check_frequency` drops, lists each `--custom` source under `Annotation sources` and prints `-i stdin` as given. (#115)
  - `--overlaps` reports `OverlapBP` and `OverlapPC` on regulatory-feature and motif-feature rows of a structural variant as well as on transcript rows. (#78)
  - `Annotator::use_transcript_ref` and `Annotator::bam_edited` return the two flags a BAM-edited cache switches on, so an embedder's `FieldOptions` write `GIVEN_REF` and `USED_REF` as the `vep` binary does. (#115)
- Picking, filtering and summaries
  - `--refseq` and `--merged` annotate the RefSeq and merged transcript sets, under `--json_cache` too; `--all_refseq` keeps every model and `_source_cache` is a `--transcript_filter` key. (#115)
  - `vep filter` filters a VCF or tab output by its fields as Ensembl VEP's `filter_vep` does, `--ontology` answered from the built-in Sequence Ontology table and rewriting the first `Consequence is` unit's exact text wherever it recurs; every INFO item under the consequence field's name is read, so a `--keep_csq` output passes on its new annotation, `--only_matched` keeps each item's passing entries, and a repeated scalar or switch option keeps its last value, as Getopt::Long does. (#115)
  - `vep filter` reads a line, a command-line word, an `in` list file and a `match` pattern as bytes, as `filter_vep` does, writing a line that is not UTF-8 as the bytes it read and comparing its fields byte by byte, and `--list` sorts by byte. (#115)
  - A `vep filter` `match` pattern reads as Perl reads it: case folds fully under Unicode rules, so `ss` matches `ß`, braces, octal escapes and flags are Perl's, the `u`, `a` and `aa` flags rule their own group alone, and `\N{U+}` reads hexadecimal numbers and dotted sequences; an abbreviated field name is the pattern Ensembl VEP interpolates, `Inf` and `NaN` are numbers, and a word without `=` in a tab line's `Extra` column is undefined. (#115)
  - `filter_vep`'s one-pattern reading of a VCF line's consequence entries, an empty entry's next item and a key ending in the name, and its reading of a consequence field a short entry lacks as the INFO key of the same name are documented classes; `vep filter` reads the exact INFO key and the field as undefined. (#115)
  - `--pick`, `--per_gene`, `--summary`, `--most_severe`, `--no_intergenic` and `--coding_only` rank, pass through, fold in or drop regulatory and motif rows as Ensembl VEP does. (#78)
- Command line, configuration and errors
  - An existing statistics report is refused before the run unless `--force_overwrite`, checked before the output file; `--quiet` leaves the report complete. (#115)
  - `vep recode` translates an identifier, HGVS, SPDI, VCF, region or Ensembl-format input into its identifiers and notations from the cache, as Ensembl VEP's `variant_recoder` does, a repeated scalar or switch option keeping its last value; a `vep` option the recoder's command line never declares, which `variant_recoder` reads from a configuration file and applies, is a documented intended divergence, `recode_config_file_vep_option`. (#115)
  - `vep recode` reads `--config <file>`, `<dir>/vep.ini` and the `VEP_*` variables as `variant_recoder` does, lower sources yielding to higher, the command line's value winning and its presets kept, accepts `--debug`, reads a merged cache without `--merged`, as `variant_recoder` reads the merged set on every run, and refuses `--refseq` with `Can't use --refseq and --merged together`. (#115)
  - `vep recode --ga4gh_vrs` writes an empty `ga4gh_vrs` list, not `[null]`, for an SPDI off an `NC_` accession, and `--fields` writes `hgvs_offset` beside `hgvsc` alone, as Ensembl VEP does. (#115)
  - `--no_progress` switches off the `Processed N variants` line written to standard error every 10,000 variants and nothing else; `--quiet` still silences every line on standard error. (#115)
  - `examples/homo_sapiens_GRCh37.vcf` and `homo_sapiens_GRCh38.vcf`, 175 records each, annotate against the bundled golden caches with no download; the quick start and the release smoke tests run the GRCh37 file. (#115)
  - Hidden `--parse_threads` sizes the pool the cache shards parse on, hidden `--prewarm_ahead` sets how many chromosomes are warmed ahead of the reader, or `all`, and hidden `--reader_core` (`auto`) reserves a physical core for the reader thread on hosts with hardware threads once the input reads as inflate-bound; `--resolved_config` prints all three. (#115)
- Plugins
  - `--plugin Name,param1,param2` parameters follow one grammar across the built-in plugins: all positional or all `key=value`; a mixed list or a key given twice is refused at plugin initialization naming it, and an undeclared key too, except by CADD and SpliceVault. (#115)
  - Plugin ABI version 2: a shared library declares its fields, which get header lines and `--tab`, VCF and Parquet columns, and the row kinds it runs on, intergenic included. (#115)
  - `dbNSFP` reads the plugin's options (`consequence=`, `output_prefix=`, `pep_match=`, `transcript_match=`, `ALL`, a replacement-logic file) and the file prepared for either assembly, as Ensembl VEP's plugin does; its field names carry the `dbNSFP_` prefix under `output_prefix=1` alone. (#115)
  - `--custom` strings parse as Ensembl VEP parses them, both forms, every default; a malformed string, remote file or GFF or GTF transcript model is refused before output. (#115)
  - `--custom` VCF and BED sources annotate as Ensembl VEP's do, with per-allele INFO values, FILTER and summary statistics; `--custom_multi_allelic` and `--custom_suppress_filter` are accepted. (#115)
  - A `--custom` VCF source whose `##fileformat` major version is not 4, or names no version, is refused before output with Ensembl VEP's message. (#115)
  - `--custom` sources print Ensembl VEP's keys, header lines and column order in every format: names, fields and statistics per source, `SOURCE` added, a record-level `custom_annotations` object in JSON. (#115)
  - bigWig sources are read with Ensembl VEP's matching, statistics and span rules, and `--bigwig`, `--phyloP`, `--phastCons`, `--ucsc_assembly` and `--ucsc_data_root` expand them from a local UCSC tree. (#115)
  - pLI, AncestralAllele, Conservation, MaxEntScan, TSSDistance, UTRAnnotator and the four parameterless plugins read a surplus positional parameter as Ensembl VEP's plugins do, logging it instead of refusing the plugin, and dbscSNV reads an empty second parameter as no assembly, keeping the run's. (#115)
  - `--plugin LoF`, the name of LoFTEE's Ensembl VEP module, loads the LoFTEE built-in, so an Ensembl VEP argument string loads it unchanged. (#115)
  - dbNSFP skips a transcript's own `.` list entry under `transcript_match=1`, refuses `transcript_match=1` on a `chr`-prefixed file, refuses a replacement-logic line opening with a blank or split by non-ASCII blanks, and reads a non-UTF-8 replacement-logic file, as Ensembl VEP's plugin does. (#115)
  - `--plugin CADD,sv=<file>` scores `<DEL>`, `<INS>` and `<DUP>` records from a CADD-SV file by class and start, as Ensembl VEP's plugin does. (#115)
  - `--plugin GWAS,type=sstate,file=<path>` reads a study's harmonised summary-statistics file, compressed or not, its study, PubMed id and trait accession from the file name. (#115)
  - `--plugin PrimateAI`, `--plugin ClinPred` and `--plugin BayesDel` are built in, reading the user's own tabix score files as the Ensembl VEP plugins do. (#115)
  - `--plugin UTRAnnotator,file=<uORF file>[,max_overlap=<percent>]` annotates 5' UTR rows with gained or lost upstream start and stop codons, frameshifts, Kozak context and the existing uORF counts; the format comparator sets aside its hash-ordered effect in Ensembl VEP as `utrannotator_consequence_precedence`, beside its gained-start window class. (#115)
  - `--plugin LOEUF` and `--plugin DosageSensitivity` are built in: the gnomAD LOEUF of the transcript or the gene minimum, and the pHaplo and pTriplo dosage scores, `cover=1` included. (#115)
  - `--plugin AVI` annotates single-base substitutions with AlphaGenome Variant Impact PHRED scores from the tabix table, with `raw=1` and `links=1` for the raw score and Atlas link. (#115)
  - Two built-in plugins: `TranscriptAnnotator` appends a tabix table's columns to the transcript rows it names, and `gnomADMt` reports gnomAD mitochondrial frequencies with haplogroup and population filters (`*_filter=all` accepted), refusing at initialization a file whose `##fileformat` major version is not 4, or names no version, with Ensembl VEP's message. (#115)
  - `--plugin NMD` marks a stop-gained, frameshift or splice-site row `NMD=NMD_escaping_variant` when the variant escapes nonsense-mediated decay under Ensembl VEP's four exon-position rules; a `--minimal` split allele is judged over its own span, where Ensembl VEP judges the rejoined record's, a documented intended divergence, `nmd_minimal_rejoined_end`. (#115)
  - `--plugin AncestralAllele,<fasta>` writes `AA`, the ancestral allele over a record's span from Ensembl's ancestral-sequence FASTA, as Ensembl VEP's plugin writes it. (#115)
  - `--plugin MaxEntScan,<dir>[,SWA][,NCSS][,verbose]` scores donor and acceptor splice sites from the MaxEntScan matrices as Ensembl VEP's plugin does, to the same three decimals. (#115)
  - `--plugin SpliceVault` is built in, reading the user's own SpliceVault table as the Ensembl VEP plugin does, its seven fields nested under `splicevault` in JSON; every `key=value` parameter but `transcript_match` names a data file, every file read, and the `splicevault` object Ensembl VEP's JSON loses on a multi-allelic record's later alleles and on a repeated record, under the default match too, is a documented class, `splicevault_json_later_rows`. (#115)
  - `--plugin OpenTargets` adds the Open Targets Platform's GWAS and QTL associations from its export for Ensembl VEP, with `cols`, `study_types` and `lead_variants_only` filters, grouped per record in JSON. (#115)
  - `--plugin Conservation,<file.bw>[,MAX]` writes the GERP conservation score of a record's span from an Ensembl Compara bigWig, as Ensembl VEP's plugin writes it; an empty parameter list or a file that is no bigWig is refused at load, and its page states the text-order rule where `MAX` differs. (#115)
  - `--plugin Blosum62` writes `BLOSUM62`, the substitution score of a missense change from the plugin's embedded matrix, as Ensembl VEP's plugin writes it. (#115)
  - `--plugin Downstream` translates a frameshift's changed coding sequence into `DownstreamProtein` and `ProteinLengthChange`, as Ensembl VEP's plugin does. (#115)
  - `--plugin HGVSIntronOffset` writes the HGVSc coordinates' intron offsets as `HGVS_IntronStartOffset` and `HGVS_IntronEndOffset`; a run without `--hgvs` is refused. (#115)
  - `--plugin NearestExonJB` writes each overlapped exon's nearest junction with `max_range=` and `intronic=1`, items in transcript order where Ensembl VEP's follow hash order. (#115)
  - Three more behaviours of Ensembl VEP's NearestExonJB are documented intended divergences: one string cached per transcript and span, handed to a later allele or record (`nearest_exon_jb_shared_span_cache`); every exon without a stable id filed under one key (`nearest_exon_jb_id_less_exon_fold`); and, under `--shift_3prime 1`, the exons taken from the shifted span and the distances from the record's (`nearest_exon_jb_shifted_span`). (#115)
  - `--plugin TSSDistance` writes the distance to the transcription start site on upstream rows, on every row under `both_direction=1`, structural variants included. (#115)
  - TranscriptAnnotator, LOEUF, AVI and SpliceVault write an empty table cell as the empty string, an empty tab cell and `""` in JSON, and SpliceVault a missing cell as JSON `null`, as their Ensembl VEP modules do. (#115)
  - ClinPred, BayesDel (on the intergenic row too) and DosageSensitivity write a row short of its score cell as the empty value and an empty cell as the empty string, and SpliceVault an empty or missing events cell as the empty list, an empty tab cell and `{}` in JSON, as their Ensembl VEP modules do. (#115)
  - A fractional or negative `num_records` of a `--custom` source warns once, and a record whose chromosome only a `--custom` or `--regulatory_gff` file lists is admitted and annotated from those files, with Ensembl VEP's bare-`M` rename, where it was skipped as not found. (#115)
  - A `--custom` key resolves by its whole short name, so sources named `A` and `A_B` both print their keys, and lists the records of the whole input record before the row's allele records; a source named `SOURCE` prints its record names under the `SOURCE` field in the default, tab and VCF formats, and a `COSV` record two sources share under one short name is listed once with the fields of both, as in Ensembl VEP. (#115)
  - As Ensembl VEP does, `--custom` GFF, BED and VCF sources list a record whose start is a minus sign and digits (`-5`) from the first base, under that number, and read a start or end with an exponent or a fraction (`1e2`, `100.5`), outside a VCF `type=exact` match; a GFF or BED record whose end `tabix` fails to parse (`xyz`, `-9`) is not listed. (#115)
  - `--custom format=gff,gff_type=gencode_promoter` overlays GENCODE promoter windows as Ensembl VEP does: one record per matching window, named by its coordinates, with `type`, `feature_id` and `associated_gene` fields, a window without `ID` carrying a `null` `feature_id`; the escapes in `ID` and `gene_id` are decoded, where Ensembl VEP prints them as written, a documented divergence. (#115)
  - A `--phyloP` source is named `phyloP<N>way`, where Ensembl VEP misspells the key `phlyoP<N>way`, and a `--custom` source's `PC` field header describes the value written, where Ensembl VEP's two texts are swapped, two documented intended divergences; `--json` writes a bigWig source's statistics object on a record no interval matched (`{"count":0}`, `{}`, or `null` for `sum` alone), as Ensembl VEP does. (#115)
  - `--plugin MaxEntScan` skips the blanks around a number on a matrix line, as Ensembl VEP's plugin does, and prints an infinite score as `Inf`. (#115)
  - The comparison tools read a run's SpliceAI score files under `--spliceai-table` and set aside rows scored by a record of several comma-joined entries, which both SpliceAI corpora carry. (#115)
- Documentation, test corpora and comparison tools
  - A golden corpus manifest or comparison case names a file inside the corpus as `{corpus}/<path>` or a shared source as `{sources}/<path>`, read by both engines without an absolute path. (#115)
  - The `vep-cli` library surface is stated as a contract on the embedding page and proven to write the `vep --json` bytes on the GRCh37 dataset. (#115)
  - `generate_reference_output.sh --corpus` keeps Ensembl VEP's statistics report as `expected/<format>.summary.log` (`--no-stats` opts out), the `GRCh38-statistics` corpus carrying one per run; `--corpus --flag-set <name>` writes one flag set's expected outputs; and `--perl-env NAME=VALUE` sets a variable in the container and records it under `perl_env`, `PERL_HASH_SEED=0` fixing a plugin's hash-order output. (#115)
  - Golden corpora of the fruit fly and the mouse compare both engines under the default flags, `--everything` and `--hgvs`: Drosophila melanogaster BDGP6.54 on 2L, the mitochondrion (genetic code 5), chromosomes 4 and Y, rDNA and seven of the scaffolds, and Mus musculus GRCm39 on chromosome 19, the mitochondrion (genetic code 2) and the scaffolds, its `--hgvs` set with `--sift b`. (#115)
  - Golden corpora of yeast and Arabidopsis compare both engines under the default flags, `--everything` and `--hgvs`: Saccharomyces cerevisiae R64-1-1, the mitochondrion read under genetic code 3, and Arabidopsis thaliana TAIR10 on chromosome 4, the mitochondrion and the plastid (genetic code 11). (#115)
  - `compare_vep_stats.py` compares two `--stats_text` statistics reports section by section, and a golden-corpus run may expect a `stats_file`, which the harness compares the same way. (#115)
  - `scripts/fixtures/make_plugin_fixture.py` writes a corpus's synthetic dbNSFP and dbscSNV files, prepared and indexed per assembly as Ensembl VEP's plugins document. (#115)
  - The release page states each dataset's peak RSS per architecture beside its wall time, and the peak of a whole genome annotated as one input file when measured. (#115)
  - The installation page tables every platform vep-rs runs on, states what a run holds in memory, and states how Apptainer and Singularity pull and run the image, with `--bind` for a cache or output directory outside the host's home, current directory and `/tmp`; the workflow-managers page maps the nf-core, Ensembl and Snakemake wrappers to vep-rs. (#115)
  - The installation and embedding pages state the crates.io route: a maintainer publishes the nine crates after each release is tagged, and `cargo install vep-cli` then installs that release's `vep`. (#115)
  - The README and the species page state the genomes the concordance measurement covers, human on both assemblies and the four other golden corpora, and what any other species gets. (#115)
  - Golden corpora compare both engines on `--custom` VCF, BED, GFF and bigWig sources, on a cache built from a GTF slice against Ensembl VEP's `--gtf` route, and on the converted release 116 RefSeq and merged caches of GRCh38 and GRCh37. (#115)
  - The GRCh38 regulatory golden corpus compares both engines on its records and structural variants, an Ensembl-format twin with its bi-allelic single-base substitutions on the reverse strand among its flag sets, with every regulatory record class represented, and the GRCh38 plugin-stack corpus runs LOEUF, DosageSensitivity, NMD and gnomADc in one batch, in load order and in reverse. (#115)
  - `tools/vep-outputs-compare` builds the tuple comparator in Rust, byte-identical to `compare_vep_outputs.py` in its three modes and across a host's cores; the concordance scripts score with it by default, and `VEP_OUTPUTS_COMPARATOR=python` selects the Python. (#115)
  - `docs/plugins.md` lists all 87 Ensembl VEP plugins with data file, terms, assemblies and one of three support statuses. (#115)
  - `scripts/divergences/reevaluate.py` re-runs every intended divergence's reproduction on Ensembl VEP and vep-rs and reports each entry still present, closed or changed; the index gains a `Verified against` column. (#115)
  - Ensembl VEP's HGVSc for an allele whose span runs over a transcript's edge, the whole alternate allele read against the bases inside, is a documented divergence with a comparator class, `hgvs_clipped_reference_against_whole_alternate`, which accepts an `N` reference base in the substitution vep-rs writes over a masked reference window. (#115)
  - The tabix-backed plugins read a mitochondrial record's file under `MT`, `chrMT`, `M` and `chrM` where the names are one sequence, GRCh38 and AlphaMissense on GRCh37; a documented divergence. (#115)
  - The damaged-cache page states the refusal of a cache whose `var_type` or `variation_cols` contradicts its variation files, where Ensembl VEP writes every record without co-located variants. (#115)
  - The format comparator applies two more built-in classes, `mature_mirna_unshifted_span` and `stop_codon_unshifted_span`, for the mature miRNA and stop-codon tests Ensembl VEP reads at the unshifted span under `--shift_3prime 1`; two divergence pages describe them. (#115)
  - The contributor guide gains recipes for fixing a divergence, a performance change and porting a plugin. (#115)
  - The nf-core `ensemblvep/vep` module's command line is tested on the run-control corpus, sarek's default arguments included; the workflow-managers page states what its version probes read. (#115, #110)
  - The installation page's Bioconda section says what the repository's recipe builds, and the limitations page points at it. (#115)
  - `scripts/concordance/fastvep_plugin_tuples.py` writes the supplementary columns of fastVEP's tab output as plugin tuples through the map `fastvep_field_map.tsv`, for scoring against Ensembl VEP's plugin and `--custom` fields. (#115)
  - `run_clone_measurement.sh --fastvep-extra-args "<flags>"` appends words to every fastVEP `annotate` run and records them in `provenance.json`; words carrying `--sa-dir` take cells of one assembly. (#115)
  - The format comparator sets aside plugin fields under thirteen more documented differences, reading the GWAS catalog under `--gwas-catalog` and a named record's alleles from `--input`. (#115)
  - `scripts/concordance/score_plugin_tuples.py` scores each plugin's fields against an Ensembl VEP output as (Location, Allele, Feature, field, value) tuples, writing the counts behind a raw and an adjusted F1. (#115)
  - `scripts/parity/plugin_gate_rows.py` derives from the plugin support table what a release scores each supported plugin on and the plugins not ported, grouped by reason. (#115)
  - The release page states raw and adjusted F1 per plugin and assembly and counts the plugins not ported by reason, when the release record carries the plugin block. (#115)
  - `scripts/docs/render_release_figure.py` draws the release's comparison figure and its values file from the provenance records, for the tools `scripts/docs/comparators.tsv` lists. (#115)
  - The format comparator sets aside a `NearestExonJB` value and a `5UTR_annotation` whose items Ensembl VEP writes in hash order, under `nearest_exon_jb_item_order` and `utrannotator_annotation_pair_order`. (#115)
  - The format comparator and the corpus builder decide `downstream_utr_frame_past_coding_end` and `downstream_start_past_own_codon` by recomputing both engines' `DownstreamProtein` from the candidate cache; the plugin tuple scorer takes `--json-cache`. (#115)
  - The format comparator sets aside the HGVSIntronOffset offsets Ensembl VEP writes for an insertion directly beside a transcript, under `hgvs_intron_offset_abutting_insertion`. (#115)
  - The format comparator sets aside a known variant's `PUBMED` or `VAR_SYNONYMS` list that the cache cut at 1,024 bytes, under `colocated_list_cut_at_cache_cap`. (#115)
  - The format comparator sets aside the `HGVSc` of an insertion standing on the bases an edited RefSeq model holds, under `refseq_edited_model_seam_insertion_hgvsc`. (#115)
  - `tools/vep-plugin-pairs` extracts, sorts and pairs the plugin scorer's values on every core, and `score_plugin_tuples.py --pairing-tool` settles the rest, writing the same report and residual files byte for byte. (#115)
  - The format comparator sets aside the `HGVSp` vep-rs writes for an indel shifted from an intron into the coding sequence, under `hgvsp_shifted_indel_written_span_gate`. (#115)
  - `tools/vep-formats-remainder` reads both outputs of a default or tab format comparison on every core, withholds only record groups written alike and passes the rest on. (#115)
  - `compare_vep_formats.py` hands a comparison over to the `tools/vep-formats-remainder` build `VEP_FORMATS_REMAINDER` names, with the same report byte for byte; unset, it compares the whole outputs. (#115)
  - `scripts/concordance/plugin_class_table.tsv` declares the documented classes of plugin fields, and `score_plugin_tuples.py` reads it unless `--class-table` names another, applying every row but those reading `not met`. (#115)
  - The format comparator decides `hgvs_shift_shared_across_alleles` on a later allele's `CDS_position`, `Protein_position` and `Amino_acids` from the cache, and reads the record's own CDS pair for its Downstream and NMD values. (#115)
  - `scripts/concordance/plugin_gate_runs.toml` lists the runs a release scores its supported plugins on, and `reproduce_plugin_concordance.py` reproduces one row: both engines, the scorer and the counts in the plugin record's form. (#115)
  - `scripts/data/write_loftee_stand_ins.py` writes stand-ins for the scripts of LoFTEE's `maxEntScan` directory from the module's own call sites, so the module runs with its splice predictions off. (#115)
  - The page Reproducing the plugin concordance is a runbook for one row: the data per plugin, LoFTEE's stand-ins, Ensembl VEP's output, the driver's commands and what differs between runs. (#115)
  - The Test job runs the plugin concordance driver on one committed corpus through the pairing tool, with no container, and holds the printed rows to the committed ones. (#115)
  - `reproduce_plugin_concordance.py` gains the step `utrannotator-sweep`, which writes the seeded Ensembl VEP outputs the two combined plugin runs are scored with. (#115)
  - The format comparator and the plugin class table set aside UTRAnnotator's two fields on a structural allele's row under `utrannotator_structural_allele_text`. (#115)
  - An upstream ORF Ensembl VEP's UTRAnnotator leaves off its frameshift list is documented as an intended divergence, `utrannotator_frameshift_at_stop_boundary`, with its comparison class. (#115)
  - A lost stop codon Ensembl VEP's UTRAnnotator reports as a frameshift is documented as an intended divergence, `utrannotator_lost_stop_written_before_stop_codon`, with its comparison class. (#115)
  - `plugin_data_layout.sh` names the LOEUF and DosageSensitivity tables, `setup_plugin_data.sh` prints their sources and `run_concordance.sh` passes both plugins to the two engines. (#115)
  - `reproduce_plugin_concordance.py` asks a whole run for `--input` and a reference only where a named step reads them, so the `utrannotator-sweep` step runs alone on `--sweep-input`. (#115)
  - The embedding page gives the three sources of the assembly the binary writes and says a BAM-edited cache leaves `--use_transcript_ref` off under `--use_given_ref`. (#115)
  - The format comparator reads the assembly a corpus manifest names and, on GRCh37, sets aside under `plugin_contig_mitochondrial_spelling` the AlphaMissense fields alone. (#115)
  - The MaxEntScan section of `docs/plugins/transcript-model.md` names the matrix lines vep-rs refuses where Ensembl VEP's plugin reads on. (#115)
  - The two MaxEntScan fixture corpora each carry a record in two splice windows, and the format comparator decides `maxentscan_overlapping_windows` from the candidate cache's introns, setting nothing aside without `--json-cache`. (#115)
  - The release figure's plugin panel sums a real-data row over the genome-wide ClinVar, gnomAD and 1000 Genomes datasets of its assembly, from the plugin record's `population_data` entries. (#115)
  - The comparison tools read every SpliceAI record line at a row's allele: `spliceai_multi_entry_cell` sets aside rows scored by repeated lines or by lines that differ. (#115)
  - The plugin scorer's summary counts, under `--spliceai-table`, the rows `spliceai_multi_entry_cell` was asked about and refused, by the ground of each refusal. (#115)
  - `plugin_score_fixtures.py` writes a SpliceAI site's record lines more than once (`--repeat-lines`) or again with other scores (`--rescore-lines`). (#115)
- Dependencies and tooling
  - `x86_64-apple-darwin` joins the release archives, built natively at `x86-64-v3` for every Intel Mac macOS 26 supports. (#115)
  - The Bioconda recipe `bioconda/vep-rs` builds `vep`, `vep-cache-builder` and `vep-cache-converter` from the release source tarball at the generic CPU floors. (#109)
  - The Bioconda recipe installs its binaries from the release directory `CARGO_BUILD_TARGET` selects, lists `NOTICE` among its license files and requires `duckdb-cli` at run time. (#115, #114)
  - `scripts/bioconda/check_recipe.py` holds the Bioconda recipe to the newest published release, so it passes on a release branch, and its tests run in the Comparator Tests job. (#115)

### Changed

- HGVS, SPDI and identifiers
  - HGVSc on an indel-edited RefSeq or merged transcript model is numbered against the edited model, so an intronic offset counts from the edited exon end. (#115)
  - Frameshift, extension, deletion and duplication peptides on an organelle transcript are translated with the transcript's codon table, as its substitutions are. (#115)
  - The HGVSc of an insertion standing on an edited RefSeq model's own bases names their replacement, placed 3' on the model, where Ensembl VEP's `ins` keeps them. (#115)
  - The HGVSc of an insertion the genome shows as a duplication, standing on an edited RefSeq model's own bases, names their replacement (`refseq_edited_model_seam_insertion_hgvsc`). (#115)
  - Beside an edited RefSeq model's own bases, the HGVSc of an allele that keeps a shared reference base inserts its other bases next to that base. (#115)
  - On an edited RefSeq model the `HGVSc` of an allele sharing a base with its reference spells the alleles the clip leaves (`refseq_edited_model_unminimised_allele_hgvsc`). (#115)
  - The `HGVSc` of a substitution clipped at a base an edited RefSeq model spells differently, or beside bases the genome lacks, is a `delins` over those bases (`refseq_edited_model_unminimised_allele_hgvsc`). (#115)
  - On an edited RefSeq model the `HGVSc` of a `delins` or `del` clipped beside bases the genome lacks spans those bases (`refseq_edited_model_unminimised_allele_hgvsc`). (#115)
  - On a failed-edit RefSeq model (`BAM_EDIT=FAILED`) a row whose `HGVSc` the 3' rule moves beside a base the allele keeps prints `HGVS_OFFSET`. (#115)
  - An allele keeping one shared base beside an edited RefSeq model's own bases takes its minimised spelling's HGVSc where the 3' rule reaches or passes the exon's last base. (#115)
  - On a failed-edit RefSeq model (`BAM_EDIT=FAILED`) the `HGVSc` of an allele that keeps a shared base is the one the same bases take written as a minimised insertion. (#115)
  - `HGVSc`, `HGVSp` and `HGVS_OFFSET` of an insertion on bases an edited RefSeq model holds and the genome lacks are read there even where the 3' shift moves it (`refseq_edited_model_seam_shift_and_duplication_hgvs`). (#115)
  - The `HGVSc` of a duplication whose copied bases lie on both sides of bases an edited RefSeq model holds and the genome lacks is numbered on the model (`refseq_edited_model_seam_shift_and_duplication_hgvs`). (#115)
  - `HGVSc` of a duplication spelt as a doubled reference names the replacement of an edited RefSeq model's span where the reference crosses a genome base the model drops (`refseq_edited_model_seam_shift_and_duplication_hgvs`). (#115)
  - `HGVSc` of an insertion on bases an edited RefSeq model holds is placed within the exon, past a genome base the model drops, as its doubled reference allele is. (#115)
  - HGVSc of an insertion on an edited RefSeq model's own bases is carried 3' past the exon's last base, onto the model's bases there or into the intron (`refseq_edited_model_seam_insertion_hgvsc`). (#115)
  - An insertion left directly after an exon's last base on an edited RefSeq model's own bases is written between that base and the intron's first position (`refseq_edited_model_seam_insertion_hgvsc`). (#115)
  - `--hgvs` writes the HGVSp of an insertion or deletion shifted 3' from an intron into a complete coding sequence, Ensembl VEP only on transcripts with a frameshift intron (`hgvsp_shifted_indel_written_span_gate`). (#115)
  - `--hgvs` writes the HGVSp of a change at the stop codon from the two proteins, where Ensembl VEP describes it from the peptides of the allele's codon window (`hgvsp_stop_codon_protein_comparison`). (#115)
  - `--hgvs` writes no HGVSp where the description would name the stop alone on a coding sequence that no stop codon closes; Ensembl VEP names a stop there (`hgvsp_stop_named_on_incomplete_coding_sequence`). (#115)
  - `--hgvs` writes no HGVSp for an insertion left between an intron and a coding exon, where Ensembl VEP reads the inserted bases as the exon's own (`hgvsp_intron_end_insertion`). (#115)
  - The HGVSp of an insertion the 3' shift moves off the bases an edited RefSeq model holds at the start of its coding sequence is read on those bases. (#115)
- Caches, reference FASTA and the cache tools
  - `vep-cache-builder` settles each region's codon table from the core dumps or `--codon-table <region>=<table>`; a region named like an organelle genome that neither names refuses the build. (#115)
  - `vep-cache-builder` refuses a GFF3 or GTF line cut short of its start column, where it skipped one; a blank line and a GFF3 sequence section stay unread. (#115)
  - `vep-cache-builder` refuses a GFF3 or GTF line of four to eight columns whose start or end no number leads or whose end is negative; 0.5.0 skipped the GFF3 line. (#115)
  - A stray tab before a GFF3 or GTF record, alone or behind white space, refuses the `vep-cache-builder` build, which counts the columns the file writes. (#115)
  - `vep-cache-builder` refuses a GFF3 record inside the file's sequence section, behind its first `>` header, naming its line; it read the record. (#115)
  - Every refusal of `vep-cache-builder`'s GTF, GFF3 and protein FASTA readers opens with the file and the line; the seqid and span refusals named the line alone, the others neither. (#115)
  - `--gencode_basic` and `--gencode_primary` on a cache whose transcripts carry neither attribute warn once that every transcript is filtered; the human-specific flags on another species are documented. (#115)
  - `vep-cache-builder` reads a skipped GFF3 line's identifier without its attribute table, builds a GTF tag list only where a record keeps it, and parses a version without copying it. (#115)
- Transcript attributes, predictions and consequences
  - `TranscriptConsequence` drops every public field that copied a value of its transcript (`transcript_id`, `gene_symbol`, `biotype`, `flags`, `mane_select` and twenty more); an embedder calls the method of the same name. (#115)
  - `TranscriptConsequence` holds a shared handle on its transcript and reads its identifiers, biotype, flags, cross-references and pick criteria through accessors instead of copying them per row. (#115)
  - On an edited RefSeq model, an insertion as long as the model's own bases it stands on is written as a `missense_variant`; Ensembl VEP writes `intergenic_variant` there. (#115)
  - An insertion at a base a failed-edit RefSeq model (`BAM_EDIT=FAILED`) holds and the genome lacks is read as an insertion on the genome, a `frameshift_variant` (`refseq_failed_edit_model_insertion_term`). (#115)
  - On a failed-edit RefSeq model (`BAM_EDIT=FAILED`) an allele that keeps the bases on both sides of the model's own base takes the term and notations of an insertion. (#115)
  - A deletion inside a stop codon an intron splits is tested on the codon's own bases and reads `stop_retained_variant` where a stop codon remains (`stop_codon_overlap_genomic_span`). (#115)
  - Where an `_rna_edit` completes a model's stop codon, the stop codon overlap reads only the stop bases the genome holds, as on a model with the whole codon (`stop_codon_overlap_genomic_span`). (#115)
  - An insertion or deletion carries `stop_retained_variant` only where the edited stop codon is tested; a span without transcript coordinates, or an allele letter naming no base, keeps nothing (`stop_retained_without_stop_codon_test`). (#115)
  - An insertion before the stop codon carries `stop_retained_variant` where it leaves every residue before the codon and a stop codon at the old stop position, on either strand (`stop_codon_overlap_insertion_range`). (#115)
  - A tested insertion or deletion carries `stop_retained_variant` only where every codon before the old stop position keeps its residue; a 3' UTR stop codon brought there keeps nothing (`stop_retained_three_base_stop_codon_test`). (#115)
  - An insertion or deletion whose span holds the whole stop codon without transcript coordinates reads `stop_lost` unless the sequence the record gives shows the stop codon kept (`stop_codon_removed_without_stop_codon_test`). (#115)
  - A deletion whose span holds the whole stop codon without transcript coordinates keeps `stop_retained_variant` where it can be written outside the coding sequence (`stop_codon_removed_without_stop_codon_test`). (#115)
  - Over the stop codon, an alternate allele holding every reference base is read as its added bases: `stop_retained_variant` where a stop codon stands in place, `stop_lost` where none does. (#115)
  - A deletion whose span holds the whole stop codon is read at another place only where known reference bases repeat; an unknown reference base ends the search (`stop_codon_removed_without_stop_codon_test`). (#115)
- Output formats, columns and headers
  - `--output_format parquet` writes typed Parquet natively: ZSTD, row groups sorted by contig and position, Bloom filters on identifier columns, no external converter. (#115)
  - Every Parquet tab column has a declared type (`LIST<DOUBLE>` frequencies, `BOOLEAN` flags, `STRUCT` positions) and a cell contract that projects a row back to its `--tab` text. (#115)
  - `--output_format parquet` refuses an input out of chromosome and position order, naming the first record; `--no_check_variants_order` writes it unchecked, in input order, with one warning. (#115)
  - `contig_index` ranks the reference FASTA's sequences the cache lacks in the index's order, after the cache's chromosomes. (#115)
  - `--output_format parquet` writes a directory of parts, one per chromosome, closed at 512 MiB, with `_vep_rs_manifest.json` listing them. (#115)
  - `--force_overwrite` on a Parquet output removes only the parts and manifest vep-rs wrote and refuses a directory holding anything else. (#115)
  - `--json` numbers carry at most 15 significant digits, the precision Ensembl VEP's JSON output writes a numeric value with. (#115)
  - `--compress_output bgzip` deflates its blocks on worker threads and writes them in order; the hidden `--compression_threads` sets how many, and every count writes the same bytes. (#115)
  - `--output_format parquet` encodes each block of rows on a thread of its own while the writer takes the next rows; rows, row groups and footers are unchanged. (#115)
  - `PUBMED` and `VAR_SYNONYMS` end before an id the cache cut at 1,024 characters, which Ensembl VEP prints as an id; a last id that may be whole goes too. (#115)
  - `--total_length` prints the cDNA's length beside `cDNA_position` on a model whose `_rna_edit` attributes add or remove bases: the exons and the edits together, not the exon sum (`refseq_edited_model_tail_term`). (#115)
- Picking, filtering and summaries
  - `--check_frequency` and `--filter_common` judge a population frequency of exactly 0 and read the matched known allele's frequency, as Ensembl VEP 116.3 does. (#115)
- Command line, configuration and errors
  - `--regulatory` on a cache converted without its regulatory tree warns once and writes no `RegulatoryFeature` or `MotifFeature` row. (#78)
  - `--stats_text` alone writes the text report alone, both flags write both and `--no_stats` writes nothing; a run that wants no report file passes `--no_stats`. (#115)
  - `--plugin` naming one of the 56 Ensembl VEP plugins vep-rs does not provide, or a Perl module, ends the run with an error stating the reason, `--safe` or not. (#115)
  - `--decompression_threads` defaults to `auto` on an aarch64 host and to `0` on any other, so an inflate-bound indexed input gets its wider pool without the flag. (#115)
  - `--debug` is accepted without effect: Ensembl VEP reads it only to print the plugin directories and modules it loads, and the annotation output is unchanged. (#115)
  - `--gff` and `--gtf` are refused with an error naming the cache route: a cache built from the GFF3 or GTF with `vep-cache-builder`. (#115, #91)
  - `--bam` is refused with an error naming the cache route: Ensembl's RefSeq cache under `--refseq`, which carries the alignment edits. (#115, #90)
  - The process ends as soon as its outputs are closed, leaving the loaded cache to the operating system instead of freeing it piecemeal; `MIMALLOC_SHOW_STATS` keeps the ordinary exit. (#115)
  - A built-in plugin's remote data file (`file=https://...`) is refused by name inside `Failed to instantiate plugin`, a warning and an error under `--safe`. (#115)
  - `--fork 0` on an indexed input whose records inflate to 8,000 bytes or more leaves the inflaters their CPUs, the annotators never fewer than half the capped count. (#115)
  - The `cache` error component's meaning names a regulatory shard that cannot be read and a regulatory chromosome directory holding files but no shard. (#115)
  - Run statistics are counted into one accumulator per pool thread for the whole run, the VCF parser builds a record's variants from its chopped alleles' buffers and one shared copy of its strings, each intron's splice windows are evaluated once per allele, and the cache loader keeps each SIFT and PolyPhen matrix gzip-compressed until its first lookup. (#115)
- Plugins
  - AlphaMissense writes `am_class` and `am_pathogenicity` on the variant's missense transcript rows whose reference and alternate amino acids match its record's, and on no other row. (#78)
  - LoFTEE prints `LoF_info` numbers as Ensembl VEP's plugin does, reads an insertion's CDS end as its lower coordinate, and keeps a call without an exon number high confidence. (#115)
  - `LoFTEE` refuses `run_splice_predictions=1`, which asks for the splice-prediction layer vep-rs does not implement, instead of warning and continuing. (#115)
  - `LoFTEE` refuses a `conservation_file` naming a database, which asks for the PhyloCSF items vep-rs does not compute, instead of accepting and ignoring it. (#115)
  - On GRCh37 `LoFTEE` refuses a run that does not set `run_splice_predictions` to `0`, where Ensembl VEP's plugin runs the splice-prediction layer by default, instead of running without it. (#115)
  - gnomADc prints Ensembl VEP's `gnomAD_<key>_cov` fields, their keys chosen by the coverage file's column count, and every value as the four-decimal mean over the variant's span. (#115)
  - `force_annotate=1` lets CADD annotate structural variants from `snv=` and `indels=` files and regions answering over 100,000 rows, with Ensembl VEP's warnings; the first matching file's scores stand. (#115)
  - CADD refuses a data file whose header does not name the `--assembly` build, and its `CADD_PHRED` and `CADD_RAW` descriptions carry the data's non-commercial clause, as Ensembl VEP's do. (#115)
  - dbNSFP reads the dbNSFP version from the data file's name as Ensembl VEP's plugin does: a name carrying none fails the plugin's load, a bare version without an `a` or `c` suffix is accepted with its warning. (#115)
  - `--plugin LoFtool` writes its percentile under `LoFtool` as Perl's `%g` prints it, matches symbols case-insensitively, and reads the bundled scores file from `--dir_plugins` when given none. (#115)
  - `--plugin pLI` writes `pLI_gene_value` with two decimals, three significant figures under `sigfig`, or `pLI_transcript_value` in `transcript` mode, keyed case-insensitively, the bundled file read from `--dir_plugins`. (#115)
  - `dbscSNV` reads the download prepared for either assembly (`pos` or `hg38_pos` by name), takes the assembly as its second parameter and writes `ada_score` and `rf_score`. (#115)
  - `--plugin REVEL` writes the field `REVEL` on missense rows whose reference and alternate amino acids and, unless `no_match`, transcript its record names; `no_match` reads by position or key. (#115)
  - REVEL reads its header as Ensembl VEP's plugin does, the last leading `#` line split on whitespace, and refuses a table whose header line lacks its `#` with `could not read headers`; the documented recipe puts the `#` on. (#115)
  - `--plugin AlphaMissense` annotates missense rows per transcript under `file=`, `cols=` (`all` or `:`-joined columns) and `transcript_match=`, the fields nested under `alphamissense` in JSON. (#115)
  - GWAS reads its parameters as `key=value` alone, refusing a bare path and `db=1`, and matches the risk allele as Ensembl VEP's plugin does, without the complement. (#115)
  - GWAS writes every field of a matching record as Ensembl VEP's plugin does, empty fields included, and reports one record where several match a site. (#115)
  - `SpliceAI` prints one `SpliceAI_pred` field in the default and tab formats, nine `SpliceAI_pred_*` fields in VCF or under `split_output=1`, and a `spliceai` JSON object, as Ensembl VEP's plugin. (#115)
  - `SpliceAI` `cutoff=` flags each record `SpliceAI_cutoff` `PASS` or `FAIL` by its largest delta score instead of dropping the records below it, and refuses a value outside 0 to 1. (#115)
  - Built-in plugins that write per row annotate transcript rows alone; CADD, gnomADc and GWAS keep their record-level values on every row. (#78)
  - `--dir_plugins` is expanded as Ensembl VEP's `glob` expands it: `*`, `?`, `[...]`, `{a,b}`, `~` and `~user`; every directory match is searched, the last match first. (#115)
  - CADD reads its parameters as Ensembl VEP's plugin does: `key=value` pairs, every key but `force_annotate` naming a data file, or files by position with a trailing `1` as `force_annotate`. (#115)
  - CADD trims the bases both alleles share at either end before matching, so a substitution reducing to one base change takes that change's row. (#115)
  - GWAS reads a haplotype or interaction row SNP by SNP, each at its location with its own risk allele, where Ensembl VEP's plugin gives every location to every SNP. (#115)
  - `SpliceAI` reads each comma-joined entry of a score record's INFO cell under its own gene, where Ensembl VEP's plugin prints the whole cell as the first gene's value. (#115)
  - Plugins and `--custom` sources are prefetched for the next batch while the current one is annotated; a plugin whose prefetch reads consequence rows waits for them. (#115)
  - `--plugin Downstream` starts `DownstreamProtein` at the codon a deletion from a codon's third base changes, where Ensembl VEP's plugin starts one codon later. (#115)
  - `--plugin Downstream` reads the 3' UTR in frame after a deletion from a codon's third base to the last coding base, where Ensembl VEP reads it out of frame. (#115)
  - `--plugin Downstream` reads the record's own CDS span and allele under `--hgvs`, where Ensembl VEP's plugin reads the HGVS-shifted span (`plugin_hgvs_shifted_cds_coordinates`), so `--hgvs` changes neither of its fields. (#115)
  - `--plugin Downstream` writes an insertion's bases in place of the bases an edited RefSeq model holds and the genome lacks, the span `CDS_position` prints. (#115)
  - UTRAnnotator annotates a symbolic deletion inside one 5' UTR exon from its span and gives plain `<INS>` the uORF frameshift its stated length decides; other structural alleles keep `-`. (#115)
  - The intended divergence `utrannotator_structural_allele_text` covers the two fields vep-rs fills on a symbolic deletion and on a plain insertion's frameshift as well as those it leaves `-`. (#115)
  - `--plugin AVI` warns on a GRCh37 run that it checks no assembly and matches the file's rows by position and alleles, so the file must be on GRCh37 coordinates. (#115)
  - An unported plugin's refusal gives FATHMM the reason `Runs an external Python script against a MySQL database.` and MPC `GRCh37-only regional missense constraint scores.` (#115)
- Documentation, test corpora and comparison tools
  - NearestExonJB's `max_range=len`, which names the plugin's length key, and a LOEUF, pLI or Conservation mode word outside its set are documented intended divergences, `nearest_exon_jb_len_key_collision` and `plugin_unknown_mode_word_refused`. (#115)
  - A `cols` name holding a regular-expression metacharacter, which Ensembl VEP's AlphaMissense, OpenTargets and TranscriptAnnotator accept as a pattern, is a documented intended divergence, `plugin_cols_names_matched_exactly`. (#115)
  - A GFF3 gene, transcript or CDS `version` that is not a whole number refuses a cache build; Ensembl VEP prints its text: `gff_version_not_a_number`. (#115)
  - A GTF line repeating a key read as one value refuses a cache build; Ensembl VEP prints the values joined, an address, or loses the record: `gtf_repeated_attribute_key`. (#115)
  - Each documentation page carries its own search description; the home page's title names vep-rs fast Ensembl VEP in Rust and its schema.org record a drop-in replacement for Ensembl VEP. (#82)
  - The README states fastVEP and Ensembl VEP ClinVar times with ratios, drops vep-rs v0.4.0's re-scored rows, places Installation and Quick start before the figures and indexes the documentation site. (#77)
  - Converting an Ensembl VEP Storable cache has its own page, Cache conversion, holding the converter, its field rule, the variation files and the published conversions. (#115)
  - The Parquet page states the typed schema, its inverse to `--tab`, the ordering precondition, parts, `--parquet_partition`, the overwrite rule and reader recipes for DuckDB, Polars, Spark and Athena. (#115)
  - Concordance and golden scripts read the species from the case table, the corpus manifest and the dataset row, `homo_sapiens` by default, and run another genome's cache and FASTA. (#115)
  - `scripts/adapters/parquet_to_vep_tab.py` reads Parquet output through `pyarrow`, inverting the writer's cell contract, with `--canonical` for a two-sided comparison against `--tab`, and `compare_vep_formats.py --format parquet` compares a frequency as a number, so `0.0010` or `1.47e-05` in the tab agrees with the Parquet value. (#115)
  - The comparison oracle takes per-case stored references, unordered rows, expected refusals and a cache root in Ensembl VEP's layout, and its default and tab comparisons name the input beside a manifest, so the manifest's flag set over that input is the one applied; a JSON floor on a record-level column is not applicable. (#115)
  - `feature_oracle.sh --format-jobs N` runs and compares a case's output formats N at a time, each format's log kept beside its outputs; the report reads as a serial run's. (#115)
  - `prune_json_cache.py` resolves a corpus's records through its `--synonyms` file and keeps an empty indexed variation slice for a chromosome whose cache names no variant. (#115)
  - The concordance harness reads the user's own dbNSFP file from `DBNSFP_FILE`; both engines' plugins refuse a file name that carries no dbNSFP version. (#115)
  - The everything and pick-summary golden corpora carry their records' regulatory tree and compare the `RegulatoryFeature` and `MotifFeature` rows, the pick-summary most-severe set under `--regulatory`. (#115)
  - The `--everything` help and the regulatory pages state every row and field the flag writes; the divergence page keeps the warning, the SPDI chromosome name and the GFF index. (#115)
  - The format comparator, the golden harness, the golden corpus builder and the oracle set aside the corrected behaviours' shapes under one list of classes, each documented as an intended divergence, and `compare_vep_formats.py` places a floor where each format carries its column: a plugin's JSON object, a custom source's record-level annotations, a tab or VCF column neither header declares. (#115)
  - The plugin, input, output and shared-library pages state Ensembl VEP 116.2's behaviour beside each vep-rs rule they had stated alone; every plugin and feature page states the Ensembl VEP defects it does not reproduce as classes with their Perl lines, tiers and worked examples run on both engines, and every ported plugin has a page section linked from the support table. (#115)
  - The page of what vep-rs does not do lists only the Ensembl VEP features vep-rs lacks; what the RefSeq, custom-source, filter and recode features do is on their own pages. (#115)
  - The intended-divergence index groups every entry by category and orders each group by the measured frequency of its triggering shape, each row naming its tier and count; the counts are over 392 committed flag sets, and every count the pages state is backed by a committed record, test or data file, or the sentence keeps its claim without it, two provenance records carrying the counts. (#115)
  - The intended-divergences index lists the entries only malformed input reaches in a last section, every reproduction records its input class (`observed`, `unobserved` or `malformed`) with its basis, and `reevaluate.py --check` and `check_divergence_order.py` hold each entry's input class to its tier and to the index section its row stands in. (#115)
  - The divergence index says the 1000 Genomes high-coverage chromosome 21 dataset is a concordance dataset, not one of the cuts its rates are measured on. (#115)
  - The divergence index counts a flag set by the command line its run is given, the corpus FASTA as `--fasta` included: nine entries that need `--fasta` gain their runs. (#115)
  - The format comparator applies the structural-variant comparator's own predicates for `reference_skipped_record`, `breakend_mate_local_read` and `bnd_own_chromosome_rows_lost`, reading the input VCF and the candidate cache, and sets aside the lent `HGNC_ID` under `hgnc_id_copy`, which the VCF format carries without `--symbol`. (#115)
  - A configuration-file value keeps a literal `___SPACE___`, which Ensembl VEP turns into a blank, and a key in another case or abbreviated is dropped without Ensembl VEP's header echo. (#115)
  - The format comparator sets aside two RefSeq and merged-cache classes, each a documented intended divergence: an `HGVS_OFFSET` Ensembl VEP omits on an edited RefSeq model (`refseq_edited_model_hgvs_offset`) and the transcript rows it drops when a merged or RefSeq cache carries one stable identifier at two loci (`merged_same_stable_id_collapse`). (#115)
  - Three more RefSeq edited-model classes are set aside: `refseq_edited_model_deleted_span_term`, the default term Ensembl VEP writes for a variant inside the bases a RefSeq model's alignment deleted, and `refseq_edited_model_tail_term` and `refseq_edited_model_deleted_span_codons`, which two divergence pages describe. (#115)
  - On an edited RefSeq model, the unrotated bases Ensembl VEP spells in a shifted insertion's HGVSc (`refseq_edited_model_shifted_insertion_bases`) and, under `--hgvs`, an alternate allele equal to the model's edited base, which Ensembl VEP writes as the base changed to itself (`edited_model_no_change_hgvsc`), are documented intended divergences. (#115)
  - The format comparator applies five more classes, each named on a divergence page: the manifest classes `minimal_symbolic_alleles`, `shift_3prime_shared_predicates` and `shift_start_codon_test`, and the built-in classes `hgvs_shifted_insertion_clamped_end` and `shift_intron_term_unshifted_span`, for insertions and deletions shifted under `--shift_3prime 1`. (#115)
  - The format comparator applies two more built-in classes, `hgvs_shifted_indel_written_span_gate` and `hgvs_shifted_deletion_clamped_span`, for the HGVSc of shifted insertions and deletions at a transcript's edges; two divergence pages describe them. (#115)
  - LoFTEE's plugin page and the divergence index name the GERP bigWig zoom-level summation of Ensembl VEP's LoF plugin, where vep-rs sums the exact per-base values, and AlphaMissense's prediction handed to a sibling allele of a multi-allelic record is a documented intended divergence, `alphamissense_sibling_allele_rows`: Ensembl VEP matches every allele of the record. (#115)
  - SpliceAI's score from a row whose alleles share a base beyond an indel's anchor, which Ensembl VEP compares as written (`spliceai_trimmed_table_row_match`), and the rejoined record Ensembl VEP hands every plugin under `--minimal`, where vep-rs hands each split allele its own span (`minimal_plugin_sees_the_rejoined_record`), are documented intended divergences. (#115)
  - OpenTargets' repeated `QTL` or export without `isLead`, AncestralAllele's uncompressed-FASTA clamp and GWAS's blank-shifted summary-statistics index are documented intended divergences; `plugins_assembly_before_cache` names CADD. (#115)
  - A plugin key or `--custom` short name the tab and VCF column list already holds, which Ensembl VEP lists twice, is a documented intended divergence, `duplicate_output_column_collapsed`. (#115)
  - Ensembl VEP's bare `inf` and `nan` tokens in the JSON output, where vep-rs writes the cell's text as a string (`json_non_finite_number_token`), and its reading of the `--output_format` spelling as typed in two output tests, where vep-rs reads the value once, lower-cased (`output_format_flag_read_in_its_given_case`), are documented intended divergences. (#115)
  - Ensembl VEP's `ins` for a duplication on a reference FASTA that ends inside the transcript is a documented class, `hgvs_dup_read_as_insertion_past_the_fasta_end`, vep-rs writing the `dup` on either reference, and the `individual_trimmable_indel` class covers a `1/1` sample whose chopped alternate equals the minimised reference (`ATT/AT`): `HOMREF` on the reference row in Ensembl VEP, `HOM` in vep-rs. (#115)
  - The process id Ensembl VEP writes before a forked child's warning is a documented intended divergence, `fork_warning_pid_prefix`, set aside in the comparator's warnings mode; vep-rs writes the condition alone. (#115)
  - Two tabix-reader behaviours of Ensembl VEP are documented intended divergences: seven tabix-backed plugins' empty field on a record at position 1 or 2 of a chromosome, where its window below position 1 dies (`tabix_window_below_position_one`), and a `--custom` or plugin VCF row whose `END` lies below its position or is no number, which its tabix reader never returns (`tabix_row_with_an_end_below_its_start`). (#115)
  - A bgzipped, indexed `--custom` or `--regulatory_gff` source named without `.gz` or `.bgz` is a documented intended divergence, `tabix_source_named_without_gz`: Ensembl VEP refuses it by name. (#115)
  - A bigWig record's name under `coords` and second listing under `min` and `max` are documented intended divergences, `bigwig_coords_record_name` and `bigwig_min_max_record_listed_twice`: Ensembl VEP prints a 0-based start and walks twice. (#115)
  - The `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` fields Ensembl VEP leaves empty in a VCF output selected with `--output_format vcf` are a documented intended divergence, `output_format_vcf_fields_left_empty`. (#115)
  - Ensembl VEP's exit on a JSON output under `--terms display` or `--terms NCBI`, which vep-rs writes, is a documented intended divergence, `json_terms_label_without_rank`. (#115)
  - The hash order of Ensembl VEP's `FREQS` list under `--check_frequency`, where vep-rs keeps the record's allele order, is a documented intended divergence, `check_frequency_freqs_list_order`. (#115)
  - The `LoF` fields Ensembl VEP's GRCh37 LoFTEE plugin leaves unwritten where its GERP walk crosses a whole exon without `conservation_file` are a documented intended divergence, `loftee_gerp_walk_without_conservation_database`. (#115)
  - A cache built from a GTF or GFF3 honours `--gencode_basic` and `--gencode_primary`, which Ensembl VEP's file sources accept and ignore (`file_cache_honours_gencode_basic`), and prints an Ensembl GFF3's `ccdsid` under `--ccds` and ranks it under `--pick`, where Ensembl VEP's `--gff` source drops it (`file_cache_reads_ccdsid`), two documented intended divergences. (#115)
  - A built cache annotates the transcripts of an Ensembl GFF3's `ncRNA_gene` lines, which Ensembl VEP 116.2's `--gff` source skips as unlisted (`file_cache_reads_ncrna_gene_lines`), and a run on it writes no `SOURCE` field, which Ensembl VEP's `--gtf` and `--gff` sources write on every transcript row (`file_cache_writes_no_source_field`), two documented intended divergences. (#115)
  - A GFF3 transcript whose `Parent` list opens with an empty field is its named gene's; Ensembl VEP gives it the gene or none by hash order: `gff_parent_list_with_an_empty_field`. (#115)
  - A GFF3 translation whose first CDS line names no protein gets no HGVSp from a built cache, where Ensembl VEP writes one on the empty identifier (`gff_cds_without_protein_id_hgvsp`), and two CDS lines of one GFF3 transcript sharing a start are both read, where Ensembl VEP keeps one of them alone (`gff_cds_lines_sharing_a_start`). (#115)
  - A GTF, GFF3 or regulatory GFF line whose coordinate is malformed or zero-padded is a documented intended divergence, `gxf_malformed_coordinate_line_refused`: Ensembl VEP builds without it or on the span as written. (#115)
  - A built cache decodes a GFF3 attribute value holding an escape, where Ensembl VEP's `--gff` source prints the escape, `A%2CB` for `A,B` (`gff_attribute_escapes_read_raw`), and reads unpadded a GTF or GFF3 line's last value that is padded with blanks and closes without `;`, where Ensembl VEP keeps the blanks (`gxf_last_value_trailing_blanks`). (#115)
  - `--transcript_filter` lists read as written where Ensembl VEP misreads them: an `in` list mixing a number and other text reads by its first written entry, where Ensembl VEP reads it by a hash-order entry (`transcript_filter_in_list_numeric_by_hash_order`); a unit whose first word is `0`, a list file named with a mode character and a list naming a directory also read as written. (#115)
  - `CONTRIBUTING.md` says the 30-word cap binds the bullet a pull request adds; a release's section may merge bullets into a longer one. (#115)
  - Documentation page titles end with the project's name alone, without the version, and the site build refuses two pages that share a title or a description. (#115)
  - The comparison oracle admits `--dir` under `{corpus}` or `{sources}` in a flags cell, refusing a `vep.ini` there, and refuses a runner-set word in any dash style or case. (#115)
  - `compare_plugin_outputs.py` counts the fields of every supported plugin from a registry held to the golden corpora's references, with gnomADc, OpenTargets and dbNSFP under the names Ensembl VEP writes. (#115)
  - The comparison oracle refuses a flags word that either engine reads as a runner-set option, abbreviated or under another of its names, resolved against both engines' option lists. (#115)
  - `run_clone_measurement.sh` refuses a cell list or `--sa-dir` selection before it touches the output directory, and tests a `--fastvep-binary` named without a slash where `PATH` resolves it. (#115)
  - `fastvep_plugin_tuples.py` keeps a carriage return inside a cell as `%0D` in its row, and refuses a header naming a mapped column twice and an empty `--fields` name. (#115)
  - The format comparator's two UTRAnnotator classes accept only what the two gained-start windows explain and what Ensembl VEP's seeded runs write (`--utrannotator-sweep`, also in the plugin tuple scorer). (#115)
  - The format comparator reads the HGVSIntronOffset fields under `hgvs_shifted_deletion_clamped_span` and holds them to a `del` notation's coordinates under `hgvs_shifted_indel_written_span_gate`. (#115)
  - The format comparator decides `hgvs_shift_sequence_edge_window` without a manifest from the reference FASTA (`--reference-fasta`) and the candidate cache, and the comparison oracle passes both. (#115)
  - The format comparator decides `organelle_peptide_paths_translation_table` without a manifest from the candidate cache's codon tables (`--json-cache`). (#115)
  - The plugin tuple scorer holds two numbers equal only when they differ by at most its tolerance as a share of the larger and outright. (#115)
  - With `--json-cache`, `hgvs_clipped_reference_against_whole_alternate` takes a substitution over a transcript's edge and an intron whose two `delins` notations show the clamp against the cache's exon spans. (#115)
  - `gwas_haplotype_cross_product` also sets aside a GWAS field vep-rs alone fills from the site's own catalog row, where the two rows are annotated under different studies. (#115)
  - `splice_family_swap` also sets aside the Downstream plugin's two fields that Ensembl VEP alone writes where the two engines' Consequence sets differ by the splice-region term. (#115)
  - The two GWAS classes of the format comparator read a VCF or JSON row at a site only where its span fits its allele, as the text formats do. (#115)
  - The three sibling-allele classes of the format comparator count a sibling only where the candidate's row of that allele carries the value too. (#115)
  - The hash order of NearestExonJB's exon items and of UTRAnnotator's `5UTR_annotation` pairs, which changes between Ensembl VEP's runs, is documented as two intended divergences, `nearest_exon_jb_item_order` and `utrannotator_annotation_pair_order`. (#115)
  - The page of `utrannotator_consequence_precedence` states the class by the terms Ensembl VEP's runs under other Perl hash seeds write on the row. (#115)
  - A page states what a comparison sets aside for UTRAnnotator's gained-start window class, `utrannotator_gained_start_window`. (#115)
  - The page of `hgvs_shifted_deletion_clamped_span` says how HGVSIntronOffset's `HGVS_IntronStartOffset` differs on a row of the class. (#115)
  - The order of Ensembl VEP's `DOMAINS` list on a transcript stored in two region files, which changes with the buffer and `--fork`, is a documented intended divergence, `domains_order_region_file_copy`. (#115)
  - The page of `hgvs_shift_sequence_edge_window` shows the class on the human mitochondrial genome and counts it on the mitochondrial records of the ClinVar files. (#115)
  - The page of `organelle_peptide_paths_translation_table` counts the class on the mitochondrial records of the ClinVar dataset and states the comparator's rule for an input without a manifest. (#115)
  - Ensembl VEP's Downstream fields on a `--minimal` split insertion, computed over the rejoined record's span in the default and tab formats, are a documented intended divergence, `downstream_minimal_rejoined_span`. (#115)
  - Ensembl VEP's header-only output at `--buffer_size 0`, where vep-rs annotates every record, is a documented intended divergence, `buffer_size_zero_writes_the_header_alone`. (#115)
  - Ensembl VEP's exit on an input with more than 100 records out of position order within a buffer fill, which vep-rs annotates, is a documented intended divergence, `unsorted_input_ends_the_run`. (#115)
  - Ensembl VEP's refusal of `--check_svs` in an offline run, where vep-rs reads Ensembl's structural-variant VCF, is a documented intended divergence, `check_svs_refused_offline`. (#115)
  - Ensembl VEP's offline exit under `--hgvsg_use_accession` on a chromosome name with no `NC_` synonym of its own, where vep-rs writes every row, is a documented intended divergence, `hgvsg_use_accession_ends_the_run_offline`. (#115)
  - The pages of eight plugin entries state what the format comparator sets aside for their classes: SpliceVault, GWAS, AlphaMissense, REVEL, two of NearestExonJB, the mitochondrial contig spelling and SpliceAI. (#115)
  - Two counts are corrected: ten residues on the organelle codon-table page and five outside bases on the clipped-reference page. (#115)
  - The LoFTEE page says which three of its five classes the format comparator codes, and the annotation-store page which four of its seven further classes it sets aside. (#115)
  - The cache conversion page says what a RefSeq model the alignment's correction failed on keeps and loses: its attributes stay, the alignment's edits go, and no cDNA is cached. (#115)
  - The index of intended divergences says which rates are counted from output rows, and the page gate holds the index to 4,500 prose words. (#115)
  - The GWAS entry's page says when the format comparator sets aside a `GWAS_` field vep-rs alone fills: where both engines annotate a paired row under different studies. (#115)
  - The reference codons and amino acids of an insertion standing on bases an edited RefSeq model holds and the genome lacks are a documented intended divergence, `refseq_edited_model_inserted_span_codons`. (#115)
  - The splice-region difference's page says that Ensembl VEP's Downstream plugin writes its two fields on a frameshift row of the class and that the format comparator sets them aside. (#115)
  - The `HGVS_IntronStartOffset` and `HGVS_IntronEndOffset` Ensembl VEP's HGVSIntronOffset plugin writes for an insertion directly beside a transcript, on a row without HGVSc, are a documented intended divergence, `hgvs_intron_offset_abutting_insertion`. (#115)
  - The coding position Ensembl VEP gives a genomic stop base of a RefSeq mitochondrial model whose stop codon an `_rna_edit` completes is a documented intended divergence, `refseq_edited_model_appended_edit_pairs`. (#115)
  - An insertion of the very bases an edited RefSeq model holds and the genome lacks is `synonymous_variant`, where Ensembl VEP writes its default `intergenic_variant`; a documented divergence. (#115)
  - The `DownstreamProtein` Ensembl VEP's Downstream plugin starts one codon late for a deletion from a codon's third base is a documented intended divergence, `downstream_start_past_own_codon`. (#115)
  - The out-of-frame `DownstreamProtein` of Ensembl VEP's Downstream plugin for a deletion from a codon's third base to the last coding base is a documented intended divergence, `downstream_utr_frame_past_coding_end`. (#115)
  - The `--minimal` split allele's HGVSp page says Ensembl VEP's `Allele` column shows rotated bases on a split insertion's frameshift rows under `--hgvs`, where vep-rs prints the record's. (#115)
  - `score_plugin_tuples.py` hands the two Downstream classes an empty `DownstreamProtein` as a value, and the Downstream rule takes the pair of consequence sets `start_cooccurrence_swap` sets aside. (#115)
  - `score_plugin_tuples.py` pairs by record the values left under a key two input records share. (#115)
  - A `score_plugin_tuples.py` run whose sorted file is not ascending, or whose value has no row under its key, ends in an error line. (#115)
  - `compare_vep_outputs.py` stops with `KeyFileOrderError` and no summary on a sorted key file that does not ascend, and keeps the key files and the rows written before the stop. (#115)
  - `tools/vep-outputs-compare` ends with an error and no result file where a side's sorted keys do not ascend, and keeps the key files and the rows written before the stop. (#115)
  - The tuple comparator reads a Location as a record Ensembl VEP dropped only when every row there names a record and no Ensembl VEP row carries any such name. (#115)
  - The empty `HGVSp` Ensembl VEP writes for an insertion or deletion shifted from an intron into the coding sequence is a documented intended divergence, `hgvsp_shifted_indel_written_span_gate`. (#115)
  - The `DownstreamProtein`, `ProteinLengthChange` and `NMD` that Ensembl VEP's Downstream and NMD plugins compute at the HGVS-shifted CDS coordinates under `--hgvs` are a documented intended divergence, `plugin_hgvs_shifted_cds_coordinates`. (#115)
  - Ensembl VEP's `coding_sequence_variant` for a deletion inside a stop codon that an intron splits, where a stop codon remains, is a documented intended divergence, `stop_codon_overlap_genomic_span`. (#115)
  - The `HGVSc` Ensembl VEP writes for an allele that keeps a shared base on an edited RefSeq model is a documented intended divergence, `refseq_edited_model_unminimised_allele_hgvsc`. (#115)
  - `stop_retained_variant` written by Ensembl VEP for an insertion or deletion whose stop codon test it does not make is a documented intended divergence, `stop_retained_without_stop_codon_test`. (#115)
  - Ensembl VEP's `coding_sequence_variant` or `frameshift_variant` for a forward-strand insertion that keeps the stop codon, where vep-rs writes `stop_retained_variant`, is a documented intended divergence, `stop_codon_overlap_insertion_range`. (#115)
  - The `synonymous_variant` Ensembl VEP writes for an insertion at the twice-read position of a failed-edit RefSeq model is a documented intended divergence, `refseq_failed_edit_model_insertion_term`. (#115)
  - The `HGVSp` Ensembl VEP writes for a change at the stop codon is a documented intended divergence, `hgvsp_stop_codon_protein_comparison`. (#115)
  - The stop Ensembl VEP's `HGVSp` names on an incomplete coding sequence is a documented intended divergence, `hgvsp_stop_named_on_incomplete_coding_sequence`. (#115)
  - The `HGVSp` Ensembl VEP writes for an insertion left between an intron and a coding exon is a documented intended divergence, `hgvsp_intron_end_insertion`. (#115)
  - The `HGVSp` Ensembl VEP writes for an insertion before an intron's end, on a transcript whose short intron widens its coding test, is a documented intended divergence, `hgvsp_coding_test_widened_by_a_short_intron`. (#115)
  - Ensembl VEP's `stop_retained_variant` on an insertion or deletion whose protein changes before a stop codon left at the old stop position is a documented intended divergence, `stop_retained_three_base_stop_codon_test`. (#115)
  - Ensembl VEP's untested `stop_retained_variant` for an insertion or deletion whose span holds the whole stop codon without its coordinates is a documented intended divergence, `stop_codon_removed_without_stop_codon_test`. (#115)
  - Ensembl VEP's notations for an insertion shifted off an edited RefSeq model's added bases, or a reference doubled across them, are a documented intended divergence, `refseq_edited_model_seam_shift_and_duplication_hgvs`. (#115)
  - The entry `refseq_edited_model_tail_term` documents a second field: the unedited exon sum Ensembl VEP prints as the total beside `cDNA_position` under `--total_length` on an edited RefSeq model. (#115)
  - Ensembl VEP's `feature_elongation` on a structural insertion point just outside an edited RefSeq model's renumbered end base is a documented intended divergence, `refseq_edited_model_edge_insertion_elongation`. (#115)
  - The page of `refseq_edited_model_appended_edit_pairs` counts, under `--hgvs`, the rows of its six models where Ensembl VEP alone, both engines differently, or vep-rs alone print an `HGVSp`. (#115)
  - Ensembl VEP's default `intergenic_variant` on an insertion three bases short of an edited RefSeq model's added bases, where vep-rs writes `inframe_deletion`, is a documented intended divergence, `refseq_edited_model_shorter_insertion_term`. (#115)
  - The page of `downstream_start_past_own_codon` says which value it keeps where a deletion's repeat reaches an exon's last base, names the other reading and counts those ClinVar rows. (#115)
  - The coding reference Ensembl VEP reads across a whole intron of an edited RefSeq model is a documented intended divergence, `refseq_edited_model_intron_span_coding_reference`. (#115)
  - The page of a `--minimal` split allele's HGVSp says what Ensembl VEP's second pass writes for a split insertion whose notations do not shift. (#115)
  - The `--use_transcript_ref` page says where Ensembl VEP's reference allele differs on the cache's own genome: after the base a failed-edit RefSeq model holds and the genome lacks. (#115)
  - The page of `hgvs_shift_sequence_edge_window` names the `SHIFT_LENGTH` that differs under `--shift_3prime 1 --shift_length` on a deletion near the mitochondrial genome's end, which stays outside the class. (#115)
  - The page of `shift_intron_term_unshifted_span` names a deletion that ends past the transcript's last base, which keeps Ensembl VEP's `intron_variant`, has no HGVSc there and stays outside the class. (#115)
  - The page of `shift_start_codon_test` states the conditions of its overlap test and covers a deletion starting before the first exon, one shifted inside the intron, and a whole-transcript deletion. (#115)
  - The page of `stop_codon_unshifted_span` covers a one-base insertion after the stop codon's first base, which keeps Ensembl VEP's `stop_retained_variant` beside `3_prime_UTR_variant` and stays outside the class. (#115)
  - The protein notation of an insertion the 3' shift moves onto bases an edited RefSeq model holds and the genome lacks is a documented intended divergence, `refseq_edited_model_inserted_span_hgvsp`. (#115)
  - The entry `refseq_edited_model_unminimised_allele_hgvsc` has a page of its own, `docs/divergences/refseq-edited-model-unminimised-allele-hgvsc.md`, with its worked records and four further forms. (#115)
  - The page of `refseq_edited_model_seam_insertion_hgvsc` documents the same insertion in an untranslated region and one the 3' shift moves onto a model's added bases; nothing is set aside for either. (#115)
  - The tuple comparators set aside the five stop codon pair classes of the 116.2 registry, read from the two rows of a pair, in `compare_vep_outputs.py`, `vep-outputs-compare` and `compare_sv_concordance.py`. (#115)
  - The HGVSp of an insertion the 3' shift moves off the bases that open an edited RefSeq model's coding sequence has a page, an index entry and a reproduction. (#115)
  - The tuple comparators and the golden corpus builder take, under `stop_codon_overlap_insertion_range`, an insertion between neighbouring positions alone and, under `stop_retained_three_base_stop_codon_test`, at most one in-frame or protein-altering term beside `stop_retained_variant`. (#115)
  - `hgvs_intron_offset_abutting_insertion` also sets aside the offsets 0 and -1 Ensembl VEP writes for an insertion directly before a forward-strand transcript's first base. (#115)
  - `hgvs_clipped_reference_against_whole_alternate` also sets aside the HGVSIntronOffset plugin's two fields that Ensembl VEP writes beside the insertion or duplication notation of an allele's outside bases. (#115)
  - The plugin tuple scorer's class table declares `hgvs_clipped_reference_against_whole_alternate` for the HGVSIntronOffset plugin's two fields. (#115)
  - A splice or intron term Ensembl VEP writes from positions past the end of a record's span is a documented intended divergence, `splice_term_past_the_record_span`. (#115)
  - The tuple comparators, the format comparator and the golden corpus builder set aside `splice_term_past_the_record_span`, reading the record's alleles from the first column or the input records. (#115)
  - The plugin tuple scorer sets aside, under `plugin_value_follows_splice_term_past_the_record_span`, a Downstream, NMD or LoFTEE field one side alone writes beside that class's Consequence. (#115)
  - Seven golden corpora document their keys under `splice_term_past_the_record_span`. (#115)
  - The contributor guide and CONTRIBUTING.md name the Advisory Reachability check's two commands, and a test holds the workflow to an audit and a minimum-toolchain build of every tool root. (#115)
  - An upstream start codon Ensembl VEP's UTRAnnotator reports lost after an insertion that leaves it whole is documented as an intended divergence, `utrannotator_lost_start_after_insertion`. (#115)
  - A class that reads a record's rows finds the record of a row under an identifier several input lines carry by the line that writes the row. (#115)
  - Ensembl VEP's REVEL score from a table row of another reference amino acid, which vep-rs does not write, is a documented intended divergence, `revel_other_reference_amino_acid_rows`. (#115)
  - The plugin tuple scorer sets aside, under `plugin_value_on_reference_skipped_record`, a DosageSensitivity or UTRAnnotator field vep-rs alone writes on a record Ensembl VEP drops. (#115)
  - The comparators' `plugin_hgvs_shifted_cds_coordinates` reads, on a later allele's row, the CDS pair an earlier allele's notation left, for the Downstream and NMD plugins' fields. (#115)
  - The comparators' `hgvs_clipped_reference_against_whole_alternate` takes the HGVSIntronOffset plugin's two fields vep-rs alone writes, or both engines write, beside an insertion; the HGVSc beside vep-rs's fields stays a difference. (#115)
  - The plugin tuple scorer sets aside, under `blosum62_hgvs_shift_shared_across_alleles`, the `BLOSUM62` Ensembl VEP reads at a moved span on a later allele's row. (#115)
  - The HGVSc, `HGVS_OFFSET` and HGVSIntronOffset offsets Ensembl VEP writes for the bases an allele adds to a reference it keeps are a documented intended divergence, `hgvs_added_bases_left_at_the_span_end`. (#115)
  - `hgvs_added_bases_left_at_the_span_end` covers added bases Ensembl VEP leaves at the end of the record's span where the 3' rule places them at a transcript's first or last base. (#115)
  - The HGVSc, `HGVS_OFFSET` and HGVSIntronOffset fields vep-rs writes for inserted bases written outside a transcript and placed on it by the 3' rule are a documented intended divergence, `hgvs_added_bases_written_outside_the_transcript`. (#115)
  - The notation Ensembl VEP writes where its 3' shift stops on its bound, for an insertion whose run of placements goes on, is a documented intended divergence, `hgvs_insertion_left_where_the_shift_stops`. (#115)
  - `hgvs_clipped_reference_against_whole_alternate` also takes the rows of a deletion written with a base it keeps, running over a transcript's first or last base and writable nowhere wholly inside the transcript. (#115)
  - On the row of a deletion over a transcript's first or last base, `hgvs_clipped_reference_against_whole_alternate` sets aside the `del`, `delins` or substitution Ensembl VEP writes where vep-rs writes no HGVSc. (#115)
  - `hgvs_clipped_reference_against_whole_alternate` sets aside the HGVSIntronOffset plugin's two fields beside the `del`, `delins` or substitution Ensembl VEP writes for a deletion over a transcript's first or last base. (#115)
  - The page of `hgvs_intron_offset_abutting_insertion` counts the rows the entry keeps, those on which vep-rs writes no HGVSc. (#115)
  - The pages of `hgvs_intron_offset_abutting_insertion`, `hgvs_shifted_insertion_clamped_end` and `hgvs_shifted_indel_written_span_gate` state the HGVSc vep-rs writes on a row outside a transcript whose inserted bases the HGVS 3' rule places on it. (#115)
  - The tuple comparators read a coordinate in the digits `0` to `9` alone: a position, span end, mate position or Location in other digits is no number. (#115)
  - A test of the tuple comparator runs its Rust build beside it over fabricated output pairs and holds both to the same reports; the Test job runs it. (#115)
  - The concordance README's rows for the format comparator and the plugin scorer say what each reads from the input records, and under which classes. (#115)
  - The documented divergence `utrannotator_lost_start_after_insertion` covers a start codon that an insertion written after its first or second base only moves on. (#115)
  - The record and pages of `utrannotator_gained_start_window` show a start codon the inserted bases make, and its comparison class refuses a start codon the insertion only moves on. (#115)
  - A documented divergence, `utrannotator_gained_start_remade_in_place`, covers the start codon Ensembl VEP reports gained when an insertion written inside it leaves it in place. (#115)
  - The format comparator's `hgvs_insertion_left_where_the_shift_stops` takes an insertion's row whose Ensembl VEP duplication, a copy around an intron, counts fewer bases than the allele holds, two fewer excepted. (#115)
  - `utrannotator_structural_allele_text` also sets aside the two UTRAnnotator fields vep-rs fills on a symbolic deletion inside one 5' UTR exon and on a plain `<INS>`'s frameshift. (#115)
  - The comparators' `hgvs_clipped_reference_against_whole_alternate` also sets aside the HGVSIntronOffset end offset Ensembl VEP writes as -1 beside an insertion from one base before a forward-strand transcript's first base. (#115)
  - The documentation home page opens with the release's largest whole-genome speed-up and comparison figure; `scripts/docs/render_home_page.py` writes every number on it from the release's records. (#115)
  - The release's comparison figure runs to the content column's width on the Benchmarks page, up to its own 1080 pixels, on a white ground in every theme. (#115)
  - The page of `splice_term_past_the_record_span` says the structural-variant comparator finds a row's record under `--input-dir`. (#115)
  - `setup_plugin_data.sh` stages the open plugin data sources alone; a plugin whose data cannot be redistributed reads the user's own file, which `run_concordance.sh` takes from an environment variable. (#115)
  - `compare_vep_formats.py` finds a record's rows by the input line that writes them: a record of many rows under an identifier of several lines takes its rows times its lines. (#115)
  - The page of `stop_codon_unshifted_span` shows the difference on a deletion that `--shift_3prime 1` carries over a transcript's stop codon. (#115)
  - The page of `hgvs_shift_shared_across_alleles` states the one consequence term by which the engines differ on the second allele alone. (#115)
  - The stop codon class of the format comparator reads Ensembl VEP's residue and stop inserted before the last residue against vep-rs's stop unchanged, and its page shows the shape. (#115)
- Dependencies and tooling
  - The nine crates carry crates.io metadata naming the project (description, homepage, keywords, categories, MSRV) and pin siblings to the exact workspace version, so `cargo publish --workspace` publishes them. (#83)
  - Each of the nine crate packages carries `LICENSE` and `NOTICE`, and `cargo install vep-cli` installs `vep` alone: `bench_plugins` builds under the default-off `bench` feature. (#115)
  - Plugin libraries load through `libloading` 0.9; a plugin that fails to load or lacks a symbol is reported with the dynamic loader's own reason. (#115)
  - The cache builder downloads through `reqwest` 0.13 with native TLS without ALPN as its backend, so the Linux binary links OpenSSL and downloads speak HTTP/1.1. (#115)
  - A manual run of the Docker workflow refuses a tag older than 0.6.0 by name: an older `vep` needs the `duckdb` CLI for Parquet output, which the image lacks. (#115)
  - The release page may run to 800 lines, raised from 400, and its renderer refuses a page over 120,000 characters. (#115)
  - The format comparator, the golden harness and the comparison oracle read a corpus item's `vep_rs_consequence_sets_with_reference_sequence` for a run with a reference sequence, and `--run` reads one run's items alone. (#115)

### Removed

- Output formats, columns and headers
  - `--parquet_shape`, `--keep_intermediate` and `--parquet_row_group_size` are refused as unknown arguments; the Parquet writer has one flat row shape and no intermediate file. (#115)
- Dependencies and tooling
  - The `duckdb` CLI is dropped from the container image: Parquet output needs no external tool. (#115)

### Fixed

- Input formats and records
  - A `<DEL:*>` the minimisation of a bi-allelic record leaves writes no row, as in Ensembl VEP; `C <DEL:*>C` was annotated as a one-base deletion. (#115)
- HGVS, SPDI and identifiers
  - `--hgvs` writes HGVSc on a `transcript_ablation` row for a deletion covering the whole transcript, clamped to its first and last bases, as Ensembl VEP does. (#115)
  - A shifted allele rotates on the transcript's strand, so terms, codons, `Allele`, HGVSp and `USED_REF` read the variant's bases; Ensembl VEP errs on long shifts and reverse-strand records. (#115)
  - A 3'-shifted insertion or deletion carried into or over a transcript's edge gains HGVSc and `HGVS_OFFSET`, and a deletion straddling the first base is described over its whole shifted span. (#115)
  - The HGVSc of an allele whose span runs over a transcript's edge describes the change inside the transcript, none when there is none; Ensembl VEP reads the whole alternate allele against the bases inside. (#115)
  - The `--process_ref_homs` reference row of a reference-homozygous copy `--lookup_ref` rewrote carries no `SPDI`, as in Ensembl VEP: the rewritten string lists the reference once. (#115)
  - `--hgvsg_use_accession` names a chromosome written under a synonym by its `NC_` accession, and a scaffold or patch by its RefSeq name, as Ensembl VEP's database mode does on human. (#115)
  - `--hgvsg_use_accession` names a sequence the synonyms table gives no RefSeq name as the cache names it, and `vep recode` writes a scaffold's or patch's SPDI on its RefSeq name. (#115)
  - `vep recode` names an SPDI's sequence as the cache does on a row off any transcript, and `--hgvsg_use_accession` takes an accession's versioned form where the synonyms table lists both. (#115)
  - The help of `--hgvsg_use_accession` names the human chromosome as the one whose accession Ensembl VEP prints with a database. (#115)
  - On the mouse cache `--hgvsg_use_accession` named sequences by the previous assembly's accession version, whichever the synonyms table listed first; it takes the highest version. (#115)
  - `--ga4gh_vrs` names its object by the accession's latest form, the versioned one over one without a version, as the HGVSg does under `--hgvsg_use_accession`. (#115)
  - `HGVSp` under `--hgvs` is left out for an insertion the 3' shift leaves after the last genomic base of a model whose stop codon an `_rna_edit` completes. (#115)
  - `HGVSp` of an untested stop codon insertion or deletion keeps Ensembl VEP's notation where the change runs across exons and is a frameshift or extension inside one exon. (#115)
  - Under `--hgvs`, bases an allele adds at a transcript's last base that repeat its last bases are described as a duplication of those bases. (#115)
  - Under `--hgvs`, bases a run of placements carries across a transcript's first base are described at their most 3' placement, as a duplication or an insertion. (#115)
  - Under `--hgvs`, the HGVSIntronOffset plugin's two fields beside bases an allele adds at a transcript's end are those of the notation vep-rs writes for them. (#115)
  - Ensembl VEP writes bases an allele adds at a transcript's end as an insertion one base outside the transcript, an insertion held at the record's span, or nothing. (#115)
  - Under `--hgvs` and `--shift_3prime 0`, inserted bases a record writes directly past a transcript's last base or before its first get an HGVSc and an `HGVS_OFFSET`. (#115)
  - The HGVSc and `HGVS_OFFSET` of inserted bases a record writes directly outside a transcript are those the HGVS 3' rule gives them on the transcript. (#115)
  - A row whose inserted bases the HGVS 3' rule places on the transcript stays `upstream_gene_variant` or `downstream_gene_variant` beside its HGVSc. (#115)
  - Under `--hgvs`, an insertion whose run of placements is longer than the 3' shift walks is described at the end of its run. (#115)
  - Under `--hgvs`, an insertion whose run of placements is longer than the 3' shift walks gets no HGVSc where its run ends past the transcript's last base. (#115)
  - Under `--hgvs`, a deletion written with a base it keeps, running over a transcript's first or last base and writable nowhere wholly inside the transcript, gets no HGVSc. (#115)
  - The genomic description alone holds a deletion that gets no HGVSc at a transcript's first or last base. (#115)
  - `--hgvs` writes the stop at its residue number as unchanged where an insertion before the stop codon leaves the predicted protein unchanged, as the HGVS nomenclature gives it. (#115)
  - `--hgvs` writes the stop at its residue number as unchanged where an insertion inside the last residue's codon, not a whole number of codons, leaves the predicted protein unchanged. (#115)
- Caches, reference FASTA and the cache tools
  - `vep-cache-builder --include-variations` refuses a variation VCF whose `##fileformat` major version is not 4 with Ensembl VEP's message, as its dumper does; it was read. (#115)
  - `vep-cache-builder` refuses a GTF or GFF3 `version` that is no whole number or above 4294967295, and a protein FASTA header version above it, naming the file and line. (#115)
  - `vep-cache-builder` refuses a GTF line repeating a key read as one value, and an escape spelling bytes that are not UTF-8, naming the key; the last value was kept. (#115)
  - `vep-cache-builder` refuses a transcript whose translation starts or ends in no exon, on either strand, naming the file, the CDS line and the base; a cache was written. (#115)
  - `vep-cache-builder` refuses a transcript whose translation ends at an offset outside the exon it is counted in, naming the CDS line, the offset and exon; a cache was written. (#115)
  - `convert-plugin-data` verifies a store against the table through the tabix reader, leaves no failed store, and refuses a cell longer than a length field records; the length was truncated silently. (#115)
  - `convert-plugin-data` removes a store whose verification fails to run, as one that mismatches, and refuses a position of more rows than the index counts; the unchecked store stayed. (#115)
  - `vep-cache-builder --include-variations` reads a variation VCF header as Ensembl VEP's dumper does: only a `##fileformat` line outside major version 4 ends the build; one without the line was refused. (#115)
  - A cache built from a GTF or GFF3 takes a translation's identifier and version from its 5'-most CDS line, as Ensembl VEP does, where it took the lowest-start line's, and with a protein FASTA takes the version from the FASTA header when the CDS line has none, where it had none or 1; a header whose identifier is the CDS line's `protein_id` whole states no version. (#115)
  - A built cache reads a GTF value as Ensembl VEP does, every quote removed and the `%XX` escapes decoded, so `"A%2CB"` is `A,B`, where the outer quotes alone came off, and a GFF3 attribute value with a raw second `=` ends before it, in a built cache and under `--regulatory_gff`, where the rest was kept. (#115)
  - A GFF3 or GTF attribute whose value is `0` is absent, as it is to Ensembl VEP, so `ID=0;Name=foo` names the record `foo`; the `0` was read. (#115)
  - `--regulatory_gff`, `--custom` GFF sources and `vep-cache-builder` decode every `%XX` escape of a GFF3 attribute value, in either hex case; six upper-case escapes alone were decoded. (#115)
  - `vep-cache-builder` reads a GFF3 transcript model's type and source columns decoded, as Ensembl VEP's `--gff` source does; a `g%65ne` line was skipped as an unlisted type. (#115)
  - `vep-cache-builder` reads the Ensembl `gene:` and `transcript:` prefixes in any letter case wherever it tests them, as it strips them; `Gene:` was no prefix. (#115)
  - `vep-cache-builder` exits 1, naming the file, when no transcript is read from a GFF3 or GTF; it failed at `info.json` or exited 0 with a cache `vep` cannot open. (#115)
  - `vep-cache-builder` exits 1, naming the file and the counts, when a build with a protein FASTA or core dumps sets every transcript aside; it failed at `info.json`. (#115)
  - `vep-cache-builder` reads a GFF3, GTF or protein FASTA with byte-order marks as the file without them; a mark stayed in a record's seqid and lost a FASTA header's protein. (#115)
  - `vep-cache-builder` logs its per-chromosome lines in name order, and the contig order of `info.json` and the order of a shard's transcripts are the same on every run. (#115)
  - A GTF transcript's single `tag` value is split on commas, as Ensembl VEP splits it, where `tag "basic,Ensembl_canonical"` was one tag, and a single GTF `tag` of `"0"` or `""` lists no tag, as it lists none to Ensembl VEP, where a `"0"` listed the tag `0`. (#115)
  - A built cache splits a GFF3 `tag` list on the comma alone and matches each member whole, as Ensembl VEP does, where `tag=gencode_basic, Ensembl_canonical` was canonical, and reads a GFF3 `tag` value of `0` as no tag, where it read the tag `0`. (#115)
  - A GFF3 or GTF record whose parent list opens with an empty field is its own record, or the transcript of a gene it names. (#115)
  - `vep-cache-builder` skips a GFF3 or GTF line of an unlisted type, as Ensembl VEP does, refusing only a start or end its index build fails; a zero start refused it. (#115)
  - The cache builder's warning on a parent in no line states what it does with the records naming it; it named every child as skipped. (#115)
  - `--fasta` of `vep-cache-builder` reads a bgzipped protein FASTA whole; a file of several gzip members stopped at its first block. (#115)
  - A merged cache's RefSeq model sharing a dbID with an Ensembl model elsewhere on the chromosome loads instead of being dropped as a duplicate. (#115)
  - `vep-cache-converter` skips a `variations/` directory or file whose name is not UTF-8 with its unexpected-name warning instead of copying it under a chromosome name no record reaches. (#78)
  - `convert-plugin-data` converts a cell that is not UTF-8 as U+FFFD, the form a tabix query gives it, where the data file was refused, and reads a header row carrying no `#`, as the plugins do, where the row's start cell failed the conversion. (#115)
  - `convert-plugin-data` removes a store's index before it rewrites the data file and refuses an output named `.vpdi`; a write cut short left a store that opened. (#115)
  - A merged cache's RefSeq models print the `HGNC_ID` of a same-symbol Ensembl model, as Ensembl VEP does; the column was empty on their rows. (#115)
  - A native cache built from core dumps without the phenotype tables stores `gene_phenotype` 0 on every transcript, as the Ensembl VEP dumper does; the key was omitted. (#115)
  - `vep-cache-builder` matches a protein FASTA record by its whole `transcript:` identifier; a transcript whose identifier holds a dot of its own took no record. (#115)
  - `vep-cache-builder --gtf` beside `--gff3` gives a transcript the tags of the line naming it whole; transcripts sharing an identifier's first component shared `cds_start_NF` and the other kept tags. (#115)
  - `vep-cache-builder` sets aside a protein FASTA record whose length is neither its transcript's whole codons nor one fewer, translating from the genome; four TAIR10 mitochondrial transcripts took such records. (#115)
  - `convert-plugin-data` verifies a store inside each record's span, at its ends and neighbours, and through batch queries, at fixed starts; a store losing a spanning record passed. (#115, #104)
  - On a RefSeq model whose coding sequence begins or ends in bases an `_rna_edit` adds and the genome lacks, coding rows carry `CDS_position`, `Protein_position`, `Amino_acids` and `Codons`. (#115)
  - The cache builder reads every gzip member of a compressed core table dump; a dump written as several members ended at its first member's last row with no error. (#115)
- Transcript attributes, predictions and consequences
  - `TSL` carries a GFF3 level's digit run as text, as Ensembl VEP prints it; two levels past 64 bits rank by exact value, an intended divergence. (#115)
  - A structural deletion, duplication or copy-number record inside a frameshift intron writes Ensembl VEP's terms (`feature_truncation` or `feature_elongation`, and the region term); intron or context terms were written. (#115)
  - An insertion, inversion, tandem-repeat or breakend record inside a frameshift intron writes Ensembl VEP's region term, and a span across such an intron and its exons no `intron_variant`. (#115)
  - A structural record on only an intron's first or last two bases writes `intergenic_variant`, as Ensembl VEP does; a duplication wrote an empty term, a breakend no transcript row. (#115)
  - Given `--fasta`, on a RefSeq mitochondrial model whose stop codon an `_rna_edit` completes, a record at a stop base the genome holds reads at that base's coding position. (#115)
  - An insertion three bases short of an edited RefSeq model's own bases it stands on gets its row, an `inframe_deletion`; Ensembl VEP writes `protein_altering_variant` or `intergenic_variant` there (`refseq_edited_model_shorter_insertion_term`). (#115)
  - An `Annotator` built by hand over transcripts loaded without the reference sequence pairs a stop codon an `_rna_edit` completes as `vep --fasta` does, and `bench_plugins` takes `--fasta`. (#115)
  - `LazyTranscriptIndexes::with_reference` is public: a library caller that builds its own transcript index hands it the reference sequence and reads an edited RefSeq model's stop bases as `vep --fasta` does. (#115)
  - A longer alternate allele's bases past the record's span give no splice donor, acceptor, donor 5th base, donor region, polypyrimidine tract or intron term; Ensembl VEP writes one. (#115)
- Output formats, columns and headers
  - A record's rows follow Ensembl VEP's transcript key, the stable id with its database id, so a versioned id prefixing another keeps Ensembl VEP's order. (#115)
  - `--vcf` output writes its header block for a VCF input that ends before a `#CHROM` line or a record, as Ensembl VEP does; the file was empty. (#115)
  - The VCF output writes a non-UTF-8 record line from its stored columns, as Ensembl VEP does; trailing empty columns and a structural QUAL `0` were written as read. (#115)
  - `--json` renames a plugin's lowercased `chr` key to `seq_region_name`, a string, as Ensembl VEP's JSON writer renames every key of its table. (#115)
  - `--json` writes numbers as Ensembl VEP does: a padded numeric value (`93 `) as the number, an integer past 2^63 and an exponent form naming a whole number as the integer, and a cell of the ten bytes `0 but true` as the number `0`; strings and a double were written. (#115)
  - `--individual` or `--individual_zyg` with an empty list (`,`, `""`) declares no `IND` or `ZYG` column in any header, as Ensembl VEP declares none. (#115)
  - `--shift_genomic 1 --json` writes no `location` on a structural variant's transcript entries, as Ensembl VEP's writer leaves them; sequence variants keep theirs. (#115)
  - `--output_format parquet` writes beside its output and renames it into place on success, so a failed or refused replacement leaves the earlier output whole; `--force_overwrite` removed it first. (#115, #105)
  - `PUBMED` and `CLIN_SIG` split a known variant's stored list as Ensembl VEP does: a list ending in a comma writes no trailing separator and no empty JSON element. (#115)
  - A stored `pubmed` or `clin_sig` of `0`, of commas alone or of false elements prints no row field, and commas alone an empty JSON list, as in Ensembl VEP. (#115)
  - A known variant's `clin_sig` or `pubmed` element keeps its spaces, as in Ensembl VEP; they were trimmed. (#115)
  - A frequency pair prints the value between its first and second colon, as Ensembl VEP reads it: a stored `G:0.25:` prints `0.25`, where `0.25:` was written. (#115)
  - A frequency column ending in a comma names no empty allele, so a known variant whose allele string opens with `/` keeps its `AF` under `--no_check_alleles`. (#115)
  - `JsonOutputFormatter::format_record_value` returns the object carrying `vep_rs_version` alone for a record with no allele; a library caller's empty slice panicked. (#115)
- Picking, filtering and summaries
  - A full tie under the pick family and a `--minimal` split allele's rows follow Ensembl VEP's keys: a transcript's stable id with its database id, a motif's database id. (#115)
  - `--transcript_filter` reads `_gene_symbol_id` and `_bam_edit_status` on a RefSeq or merged cache, as Ensembl VEP's filter does; a filter naming either dropped every transcript. (#115)
  - `--check_frequency` filters on a cache that declares variation columns and holds no variation file, as Ensembl VEP does: `include` drops every record, where the filter was inactive. (#115)
- Command line, configuration and errors
  - `--show_cache_info` opens the `--custom` files and refuses one that is missing, unindexed or not bgzipped before printing, as Ensembl VEP does; it printed the cache's lines and exited 0. (#115)
  - `--show_cache_info` ends on a warnings file that cannot be written once a warning is drawn, as Ensembl VEP does; it printed the cache's lines and exited 0 with the warning lost. (#115)
  - `--config` naming a directory reads nothing, and `--config 0` is refused with the empty name, as Ensembl VEP reads both; the directory was refused, `0` opened. (#115)
  - `--everything` with `--hgvsg`, `--check_ref`, `--lookup_ref` or `--use_transcript_ref` and no FASTA, and `--var_synonyms` with a cache or filter refusal, print their `INFO:` line before the error, as Ensembl VEP does. (#115)
  - The `Processed N variants` line is written at every 10,000 variants in the default pipeline mode; it followed batch hand-outs, so a short run could print none. (#115)
  - A run whose standard error refuses writes, a full device or a pipe without a reader, completes and exits 0; the first line it could not write aborted it. (#115)
  - The help of `--af_gnomad`, `--freq_pop` and `--freq_filter` covers a cache that declares no variation columns: flag accepted, every population name accepted, filter inactive, as in Ensembl VEP. (#115)
  - `--shift_3prime` and `--shift_genomic` take any integer and read any but 0 as on, as Ensembl VEP does; a value other than `0` or `1` was an argument syntax error. (#115)
  - `--convert` with a value Perl reads as false, `0` or the empty value, proceeds, as in Ensembl VEP; it was refused as deprecated. (#115)
- Warnings and skipped variants
  - A record's `Transcript-assembly mismatch` and `Unable to assign consequence type` warnings are written in input order in both pipeline modes, where streaming wrote them as its tasks finished, and the mismatch warning is drawn from the rows the output prints, so a marked row `--pick`, `--per_gene`, `--most_severe` or `--summary` drops warns nothing, as in Ensembl VEP. (#115)
- Plugins
  - LoFtool, pLI, NearestExonJB and AVI read a `0` file parameter as no file, as Perl does; LoFtool and pLI then read the bundled file, NearestExonJB accepts `parent_only`. (#115)
  - A `key=value` plugin parameter takes the second `=`-field alone as its value, as Ensembl VEP's `split('=')` does (LoFTEE's `=` form excepted); `gff3=0=x` kept `0=x` and dropped NearestExonJB. (#115)
  - LoFTEE derives its default GERP source from `loftee_path` as Ensembl VEP's plugin does, and accepts `gerp_file`; a `LoF` word naming no GERP file silently skipped `END_TRUNC`. (#115)
  - LoFTEE reads a `key:value` word at its colon before any `=`, dropping a trailing colon, as Ensembl VEP's plugin does; a `loftee_path` directory holding `=` was an unknown key. (#115)
  - A LoFTEE `key:value` word of three or more colon fields (`key:value:x`, a URL) is dropped and logged, as Ensembl VEP's plugin drops it; it was read whole. (#115)
  - LoFTEE keeps the blanks in a `key:value` word's fields, as Ensembl VEP's plugin does; `run_splice_predictions: 0` was trimmed to `0`, false, where that plugin reads it true. (#115)
  - On GRCh37 LoFTEE ignores `use_gerp_end_trunc`, `gerp_end_trunc_cutoff` and `filter_position`, as Ensembl VEP's plugin there does, and judges `END_TRUNC` by the GERP rule at its cutoff of 180. (#115)
  - On GRCh37 LoFTEE honours `skip_lof_info`, dropping `LoF_info`, and accepts `tabix_path`, as Ensembl VEP's plugin there does; both were logged as unknown. (#115)
  - LoFTEE judges `run_splice_predictions` by its last value, as Ensembl VEP's plugin keeps one value per key; `run_splice_predictions:1,run_splice_predictions:0` was refused. (#115)
  - Built-in plugins skip a structural variant's rows unless their Ensembl VEP module declares the kind (CADD, DosageSensitivity, TSSDistance, UTRAnnotator); LoFtool, pLI and gnomADc wrote on them. (#115)
  - A `.vpd` store spans each record as the tabix reader does, INFO `END` and reference allele included, and a table read by other columns refuses it; a query matched starts. (#115)
  - A `.vpd` store returns the rows at one start in file order, as tabix does, and one of an earlier index version is refused; rows were sorted by end. (#115)
  - A plugin table without `<file>.tbi` beside it is refused naming that index, as Ensembl VEP refuses it; a `.csi` or stem-named index alone was read. (#115)
  - AlphaMissense `cols=` refuses a literal bare `am_` as a column the header lacks, as Ensembl VEP does; only a trailing empty name is dropped. (#115)
  - A `.vpd` store is read only beside its table: alone, its plugin is not loaded and a warning names the table; 0.5.0 opened a store without its table. (#115, #106)
  - REVEL and AlphaMissense match a substitution written with an unchanged leading base to its single-base row, trimming both sides' alleles as Ensembl VEP does. (#115)
  - Plugin parameter values are read as their Ensembl VEP modules read them, by number (`intronic=01`, `max_range=50.0`, `max_overlap=80%`) or Perl's truth (`trim=2`, `transcript_match=yes`); each dropped the plugin. (#115)
  - GWAS reads every gzip member of a `bgzip`-compressed data file, where only the first member loaded, and loads the v1.0 catalog download, and any header lacking a column a field's value comes from, with that field empty as Ensembl VEP writes it, where the header was refused. (#115)
  - `SpliceAI` matches a record on its position as well as its alleles, takes a gene's record only on that gene's transcripts, scores an unanchored multi-allelic indel allele, and keeps a truncated score entry with the parts it has, the missing parts empty and `null` in JSON, as Ensembl VEP does, where the entry was dropped. (#115)
  - LoFtool and pLI read a score cell as Ensembl VEP's `sprintf` numifies it, its leading number, warn on a word score, and give a gene listed twice its later line, a later `NA` included; pLI takes the last line containing `gene` or `transcript` as its header, as Ensembl VEP's plugin does, a comment line before or after the table included. (#115)
  - A plugin's flat file, a tabix header and a `--custom` VCF header are read as bytes, as Ensembl VEP reads them; a byte that is not UTF-8 refused the file. (#115)
  - A plugin failing at run time costs only the variants whose own span fails its fields, with one `Plugin '<Name>' went wrong` warning; the plugins after it still run. (#115)
  - A built-in plugin's run-time failure is one `went wrong` warning per plugin per run carrying the variants affected, as Ensembl VEP writes one; each malformed row wrote its own. (#115)
  - LoFTEE flags `NAGNAG_SITE` for `AG.AG` at any offset of the acceptor's 9-base context, where only a motif opening the context was flagged, and on a one-base acceptor deletion, reads `check_complete_cds` as any case of `true`, and reads an exon or intron range such as `6-7/9`, as Ensembl VEP does. (#115)
  - LoFTEE resolves its GERP source per key as Ensembl VEP's plugin does, and a failed GERP read fails its row alone under the plugin's `went wrong` warning; the last key set decided, and a failed read counted 0. (#115)
  - LoFTEE warns once on a chromosome its tabix GERP file lacks, as for a bigWig, and a plugin's `went wrong` warning states which rows lost its fields. (#115)
  - LoFTEE's off word `false` folds case on `use_gerp_end_trunc` and every data-file word, as Ensembl VEP's plugin reads `human_ancestor_fa`; `False` was a path and dropped the plugin. (#115)
  - LoFTEE writes `GERP_DIST` and `BP_DIST` before `PERCENTILE` in `LoF_info` on GRCh37, as the plugin Ensembl VEP runs there does; `PERCENTILE` was written first. (#115)
  - SpliceVault warns under `--refseq` or `--merged` with `transcript_match` on, and it, AlphaMissense, REVEL and pLI print that warning, and dbNSFP its transcript-set and column warnings and its per-row list errors, on standard error alone, under `--quiet` too, as Ensembl VEP does. (#115)
  - Conservation and LoFTEE open a bigWig themselves, so a file that is no bigWig draws the plugin's refusal alone and `--quiet` leaves standard error empty. (#115)
  - The tabix-backed plugins and a `.vpd` store keep a row's trailing whitespace, printing a last cell with it and a whitespace-only last cell as a value, as Ensembl VEP does, and gnomADc reads a coverage row's trailing whitespace as Ensembl VEP's `split` leaves it, through the table as through a `.vpd` store, where the table path trimmed it. (#115)
  - A plugin table with no header row whose first start cell carries a suffix (`100x`) has no header, as tabix indexed it; the row was taken as the header. (#115)
  - A VCF row short of its INFO cell spans its reference allele, as tabix returns it, in plugin tables, custom sources and the structural-variant check; it ended at its first base. (#115)
  - CADD and REVEL leave a variant unannotated when matching reaches a row with no allele cell, as Ensembl VEP does; the row was skipped. (#115)
  - A matched row's empty or missing cell is written as Ensembl VEP writes it: dbscSNV's empty or absent `ada_score` or `rf_score` cell and dbNSFP's empty cell as the empty value, CADD's and AlphaMissense's missing trailing cell as the empty value and their empty cell as the empty string, and REVEL's missing score cell as the empty value; the field or key was dropped, and REVEL took a later row's score. (#115)
  - A dbNSFP column named `TSL` or `APPRIS`, and any plugin key named as a core field, print the plugin's value over the cache's in the text and JSON formats, as Ensembl VEP does. (#115)
  - dbNSFP reads its README as Ensembl VEP does: any file named `*dbnsfp*readme.txt*` in any case, where a `.bak` was skipped; non-UTF-8 bytes read, only ASCII blanks and word characters counted, and a name ended at the first colon a blank follows; and the continuation lines of a line that opens no `name: description` entry dropped, where they were appended to the column before it. (#115)
  - A `.vpd` store's record names its chromosome as the store spells it, as a tabix record names the file's; dbNSFP's `chr` column and detection read the query's spelling. (#115)
  - Position cells are matched as Ensembl VEP's plugins match them: SpliceAI, dbscSNV and dbNSFP as Perl's `==` numifies the cell, where a cell with a fraction matched at its leading integer; AlphaMissense, CADD and REVEL as the text Ensembl VEP's key holds, stepped past a shared base as each module steps it; and the tabix-backed plugins read a start cell to its leading integer, as tabix does, where the row was dropped. (#115)
  - AlphaMissense declares a column `cols=` names twice once, one header line and one cell, as Ensembl VEP's hash over the list does, and reads its table from two bases before a variant's start, so a multi-base row anchored one base earlier is trimmed and matched. (#115)
  - REVEL matches a row's ninth cell as its transcript ids under an eight-name header too, as Ensembl VEP reads every row; the header's count decided. (#115)
  - A bare `--plugin` name resolves only to a `lib<Name>` library under the plugins directories, so a working-directory entry of that name is not a plugin and an unported plugin's refusal fires; a library named by its file name in the working directory (`--plugin libecho.so`) loads, where the loader was handed the bare name, which it searched for on its own path. (#115)
  - A `.vpd` store is read only when its index's recorded SHA256 and length match its table and data, else the table; a replaced table kept the earlier answers. (#115, #106)
  - A plugin table's row at the last position its index holds (2^29 for a `.tbi`) is returned, as tabix returns it; 0.5.0's reader stopped one position short. (#115)
  - UTRAnnotator reports no lost start codon for an insertion directly before the codon, which leaves it whole; Ensembl VEP reports it lost. (#115)
  - REVEL writes a table row's score only where the row's reference amino acid, as well as its alternate, is the output row's; Ensembl VEP tests the alternate alone. (#115)
  - `--plugin UTRAnnotator` reports an upstream start codon that an insertion written inside it only moves on as neither gained nor lost, reading its frame effects directly before the codon. (#115)
  - `--plugin UTRAnnotator` reports no gained start for an upstream start codon that an insertion written inside it leaves in place, and reads the frameshift directly after that codon. (#115)
  - UTRAnnotator leaves out of the frameshift a uORF that still ends on a stop codon at its own index beside a start codon an insertion only moves on. (#115)
  - UTRAnnotator reports a uORF's lost stop codon on every way of writing an insertion beside a moved start codon whose bases can stand inside that stop codon. (#115)
  - On GRCh37 a tabix-indexed plugin table other than AlphaMissense's, dbNSFP's or gnomADMt's is read for a mitochondrial record within `MT`/`chrMT` or `chrM`/`M` alone, two sequences; 0.5.0 read across them. (#115)
  - A plugin table, `--custom` file or cache variation file joined from two compressed streams is read past the empty block inside it; 0.5.0 refused a plugin table's next chunk. (#115)
  - `--plugin Conservation` prints the stored score of the one scored base a record covers, so a stored negative zero is written `-0.000`, as Ensembl VEP's plugin writes it. (#115)
- Documentation, test corpora and comparison tools
  - The compatibility page's prefix example is `--bio` is `--biotype`; `--bi`, which it named, is ambiguous between `--biotype` and `--bigwig` in both engines (exit 2). (#115)
  - The golden corpus tools refuse a missing input by name, `build_golden_corpus.py` a `--cache` without `transcripts/` or an unreadable reference and `prune_json_cache.py` a `--vcf` that is not a file, and `build_golden_corpus.py` resolves a VCF record whose ID is `;rs2`, `0;rs2`, empty or `0` to its Ensembl VEP rows by the name Ensembl VEP gives it. (#115)
  - The concordance harness reads every option spelling Ensembl VEP accepts (single dash, abbreviation, any letter case), so `-nearest` or `--NEAREST` mounts the cache writable and takes the first-run lock. (#115)
  - `run_concordance.sh` records the plugin argument strings in `provenance.json`, and `run_clone_measurement.sh` takes `--cells` and refuses an image whose digest is not the pinned one. (#115)
  - The structural-variant comparator reads INFO/SVTYPE by Ensembl VEP's truth rule, so `SVTYPE=0` names no allele, and discards a symbolic `<BND>`'s END as Ensembl VEP does. (#115)
  - `generate_reference_output.sh` reads a flag set's `--format <input format>` word in place of its own `--format vcf`, so a corpus run over an Ensembl-format, region or HGVS input regenerates; the word was refused. (#115)
  - `build_golden_corpus.py classify` classes every Ensembl VEP row on a key vep-rs writes no row for as `unexplained_residual`, and a script test refuses a committed corpus manifest carrying that class. (#115)
  - The admission check reads a model declaration only from the guide's two forms; a co-author trailer that names a person is not read as one. (#115)
  - `compare_plugin_outputs.py` honours `--max-mismatches`: the JSON sample keeps the first N mismatches per plugin per file pair, 0 keeps none and a negative value is refused. (#115)
  - The documentation site shows an image whose file lies under `docs/` from the site's own copy; the link preprocessor sent it to a GitHub page. (#115)
  - The golden corpus builder names a tuple under a stop codon class where the tuple comparator's rule takes every pair of its rows, all five classes alike. (#115)
  - The plugin tuple scorer hands the three sibling-allele classes the engine's rows of the same input record, on its own path and through `tools/vep-plugin-pairs`, which keeps them by Location. (#115)
  - `run_concordance.sh`'s cache build hint names `--core-dumps`, and `--codon-table MT=2` for a build without the dumps. (#115)
  - `tools/vep-outputs-compare` counts a coordinate at or past 2,147,478,648 under `bnd_mate_32bit_buffer_intergenic` whatever its length, as `compare_vep_outputs.py` does. (#115)
  - `compare_vep_outputs.py` reads a coordinate of more than 4,300 digits, where it ended with an error. (#115)
  - The plugin scorer and the format comparator find a row's input record under the first `;`-joined element of its ID column, the name both engines print. (#115)
  - The structural-variant comparator reads a row's reference bases by position from the input files under `--input-dir`, so it sets aside `splice_term_past_the_record_span` under an identifier. (#115)
  - `scripts/docs/build_release_site.sh` sends a redirect to its destination page when the release's own page lacks the anchor; the post-build check failed the site of v0.5.0 on three such redirects. (#115)
  - `run_concordance.sh`'s cache build hint names the run's assembly in its command, FASTA name and dumps directory; 0.5.0 named GRCh38 for every run. (#115)
- Dependencies and tooling
  - Converter tests take their own scratch directories, the `/dev/full` test is compiled on Linux alone, the derived-cells size bound follows the platform's `Once`, and three `chunks_exact` calls become `as_chunks`. (#94)
  - The container image carries `ps` from Debian's `procps` package: without it a Nextflow task that collects a trace, report or timeline exits 1 before `vep` starts. (#115)
  - On a pull request the Docker workflow builds the tree's binaries by the release recipe, glibc floor included, in a job with read-only permissions, and packages and smoke-tests them. (#115)
  - The format comparator counts a documented item once on an input that is not a VCF, so each format applies a class to the same entries. (#115)

### Security

- Security and robustness
  - Cleared from the dependency tree: RUSTSEC-2026-0442 (arrow-data 60.0.0, an out-of-bounds read validating run-end-encoded arrays, which vep-rs never builds), through arrow and parquet 60.1.0. (#115)
  - The advisory check allows RUSTSEC-2026-0543, which reports `chrono` as unmaintained and names no defect, under a dated entry. (#115)

## [0.5.0] - 2026-10-25

### Added

- Input formats and records
  - `--delimiter` splits an Ensembl-format or VCF input's columns (a tab, a space, runs of spaces, `\s+`, literal text or a regular expression); under `--format guess` the detected one applies. (#73)
  - `--check_ref` skips a record whose reference allele differs from the FASTA, `--lookup_ref` replaces it from the FASTA, and `CHECK_REF` prints `failed` on a kept record, as Ensembl VEP does. (#73)
  - `--chr` restricts the run to the chromosomes listed, matched against the input's own spelling, and reads the `1-5,10,MT` range grammar Ensembl VEP's help shows but its parser does not split. (#73)
  - Every input record passes Ensembl VEP's `validate_vf` checks before annotation: a bad coordinate, allele string, insertion or reference-allele length skips it with Ensembl VEP's warning text. (#73)
  - `--allow_non_variant` keeps a VCF record whose ALT is `.` in every output format: its line in the VCF output, one object in the JSON output, no row elsewhere. (#73)
  - Region input needs `--fasta`, is refused with `--check_ref`, and reads its reference allele from the FASTA, as Ensembl VEP's region parser does. (#73, #71, #72)
  - `--use_transcript_ref` reads each row's reference allele from the transcript, as Ensembl VEP does, and adds the `GIVEN_REF` and `USED_REF` columns; without `--fasta` the run exits. (#73)
  - Under `--check_ref` or `--lookup_ref` a contig the FASTA lacks skips the record, where Ensembl VEP annotates a one-base record at position 1 against a placeholder `N`. (#73)
  - Compressed standard input is decompressed and a stream cut inside a gzip member is exit 1, where Ensembl VEP reads standard input raw and annotates the prefix. (#73)
- HGVS, SPDI and identifiers
  - `--shift_3prime 0|1`, `--shift_hgvs 0|1`, `--shift_genomic 0|1` and `--shift_length` move an insertion or deletion 3' as Ensembl VEP does, with `SHIFT_LENGTH`, shifted `Location` and `Allele`, and an unshifted HGVS mode. (#73)
  - `--hgvsp_use_prediction` writes HGVSp in the predicted form, `p.(Arg97Ter)`, and `--remove_hgvsp_version` leaves the translation version off its reference name, as Ensembl VEP does. (#73)
  - `--spdi` writes each allele's SPDI notation as Ensembl VEP's `spdi_genomic` does, from the record's alleles as the parser leaves them, with no reference FASTA required. (#73)
  - `--ga4gh_vrs` writes the GA4GH VRS allele Ensembl VEP writes with a database connection, the `NC_` accession read from the chromosome synonyms table; refused with `--vcf`. (#73)
  - `--hgvsg_use_accession` names the HGVSg chromosome by the first `NC_` synonym of the input name, as Ensembl VEP does, and HGVSg follows `--shift_hgvs 0`. (#73)
  - `--format spdi` reads SPDI notation (`sequence:position:deletion:insertion`) as Ensembl VEP does, the reference allele checked against or read from `--fasta`, which the format needs. (#73)
  - An SPDI record is named by its notation and prints under the chromosome name the cache uses, resolved through the synonyms. (#73)
  - A notation Ensembl VEP's adaptor refuses is warned about with its text and skipped, the run going on; a span over `--max_sv_size` is skipped with Ensembl VEP's warning. (#73)
  - `--format id` resolves a list of variant identifiers from the cache's variation files in one pass, one record per genomic mapping, as Ensembl VEP's database lookup does. (#73)
  - An identifier the cache lacks is warned about as Ensembl VEP warns (`No variant found with ID`) and the list goes on; a mapping failing coordinate checks is skipped. (#73)
  - HGVS transcript notation (`c.`, `n.`) resolves offline from the transcript cache and `--fasta`, by stable id or gene name through Ensembl VEP's pick order, with UTR and intronic positions. (#73)
  - A run with HGVS transcript notations reads the cache's chromosome and scaffold identifiers once to find them, alternate loci and patches only for an id found on none. (#73)
  - HGVS protein notation (`p.`) resolves offline from the cache and `--fasta` as Ensembl VEP resolves it: substitutions and delins take the shortest nucleotide path, deletions their codons, frameshifts refused. (#73)
  - A protein notation finds its transcript by translation id, transcript id or its gene's transcripts, and the `--ambiguous_hgvs` refusal names protein notations among the inputs resolved to one variant. (#73)
  - `HGVSg` is each row's own allele's notation on a record whose ALT list repeats an allele or names the reference, where Ensembl VEP prints another entry's; documented as intended. (#73)
  - A record whose alleles shift by different lengths under `--shift_3prime 1` is annotated per allele, where Ensembl VEP shares one allele's predicates; an intended divergence page states it. (#73)
  - A version-less `NC_` accession and a span past the chromosome end in SPDI input are skipped, where Ensembl VEP's database resolves both; a divergence page states it. (#73)
- Caches, reference FASTA and the cache tools
  - `storable_to_json.pl` converts every chromosome directory of a Storable cache and exports the UniProt isoforms, protein-feature labels, SIFT and PolyPhen matrices, `chr_synonyms.txt` and `contig_order`. (#73)
  - `storable_to_json.pl` refuses a Sereal cache, a shard-less cache, an unlisted chromosome and any shard that fails to read or write, with one `ERROR:` line and exit 1. (#73)
  - `cache_parity_audit.pl` classifies array elements and the new exports, so a conversion audits with no unexplained path. (#73)
  - Every chromosome directory of a cache loads, patches, alternate loci, scaffolds and `LRG_` regions included; a cache with no transcript shard is refused with an error naming its serialiser. (#73)
  - `--cache` takes Ensembl VEP's optional directory value, the cache root when `--dir_cache` is absent. (#73)
  - `--full_cache_dir` (an alias for `--json_cache`), `--show_cache_info` (sorted `key<TAB>value` source versions, exit 0), `--no_fasta`, and the cache directory's first indexed `.fa` as the reference. (#73)
  - `storable_to_json.pl` and `vep-cache-converter` copy a cache's variation files (`all_vars.gz` with its index, or the `_var.gz` regions) into `variations/<chr>/`. (#73)
  - `vep-cache-builder --core-dumps <dir>` reads the Ensembl core and variation table dumps, so a native cache carries the symbol source, HGNC id, xrefs, protein features, phenotype counts and attributes. (#73)
  - Native caches carry Ensembl's exon phases, the N-padded translateable sequence, both UTRs, genomic-order sorted exons, transcript-order introns and codon table 2 on MT, as Ensembl VEP's does. (#73)
  - `--variation-dumps` reads the ClinVar allele and somatic classification columns from the variation database dumps, which no Ensembl VCF carries. (#73)
  - `--synonyms <file>` and the cache's `chr_synonyms.txt` resolve an input chromosome named by accession or `chr` prefix to the cache's name, as Ensembl VEP does. (#73)
  - `vep-cache-builder` writes `info.json` with `contig_order`, `sift`, `polyphen` and `variation_cols` set from what the build wrote, and `chr_synonyms.txt` from the core dumps. (#73)
  - `--fasta` reads a bgzipped FASTA through its `.fai` and `.gzi`, checking the index against the file at open, and reads a directory as every indexed `.fa` and `.fa.gz`. (#73)
  - `--fasta_dir` selects the FASTA named for the species and assembly, `toplevel` and bgzipped first, as Ensembl VEP does, and warns when none qualifies. (#73)
  - `vep-cache-builder --include-variations` writes Ensembl VEP's `variations/<chr>/all_vars.gz` and `.csi` from the per-chromosome, somatic, phenotype, clinical and 1000 Genomes VCFs, replacing the JSON variation shards. (#73)
  - vep-cache-builder and vep-cache-converter write `CHECKSUMS.sha256`; vep checks every shard, index and metadata file it reads against it and refuses a mismatch; a cache without one loads with a warning. (#73)
- Co-located variants and frequencies
  - `--check_existing` reads Ensembl VEP's own variation files, `variations/<chr>/all_vars.gz` through its index on first use, with columns from the cache's `variation_cols`; the JSON variation shards are gone. (#73)
  - `--old_maf` is accepted and rewrites each co-located variant's AFR, AMR, ASN and EUR values as Ensembl VEP does, so those columns print empty on a release 116 cache. (#73)
  - An unreadable variation file ends the run with `ERROR: [cache]` before its batch is written; the library's per-row API returns the failure on that batch's every row. (#73)
  - `--no_check_alleles`, `--exclude_null_alleles` and `--failed <int>` are accepted and match known variants on their span alone, drop null-allele records and admit failed ones as Ensembl VEP does. (#73)
  - `--clin_sig_allele`, `--clinvar_somatic_classification` and `--var_synonyms` are accepted, the last switching `--check_existing` on with Ensembl VEP's notice; `VAR_SYNONYMS` and `CLINVAR_SOMATIC_CLASSIFICATION` print as it does. (#73)
  - `--af_gnomad` is accepted as the older spelling of `--af_gnomade`, refused beside it, and fails on a cache without gnomAD exome columns as Ensembl VEP fails. (#73)
  - `--check_svs <file>` reports the known structural variants overlapping each variant (`SV`, JSON `colocated_structural_variants`) from Ensembl's structural-variant VCF offline, a breakend allele's mate position included. (#73)
  - `--check_frequency`, `--freq_pop`, `--freq_freq`, `--freq_gt_lt`, `--freq_filter` and `--filter_common` filter the input by co-located variant frequency as Ensembl VEP does, with `FREQS`. (#73)
- Transcript attributes, predictions and consequences
  - Converted caches load MANE Select, MANE Plus Clinical, APPRIS and TSL from transcript attributes and the GENCODE sets by attribute presence, so the attribute columns and `--pick` tiers fill. (#73)
  - `--overlaps` adds the `OverlapBP` and `OverlapPC` columns to the tab and VCF outputs, as in Ensembl VEP, on transcript rows of non-breakend structural variants. (#73)
  - `--gencode_basic` filters the transcript set as the cache loads, and `--gencode_primary` and `--flag_gencode_primary` arrive with Ensembl VEP's incompatibility rule, GRCh38 warning and `GENCODE_PRIMARY` column. (#73)
  - `--transcript_version` and `--gene_version` append the stable-id version to the Feature and Gene columns, as Ensembl VEP does. (#73)
  - `--nearest transcript|gene|symbol` adds `NEAREST`, the protein-coding transcription start sites nearest each record, found as Ensembl VEP's transcript tree finds them; every transcript at a shared site is listed. (#73)
  - `--domains` lists the protein features the variant's residues overlap on every coding row, as `label:identifier`, from converted and native caches alike. (#73)
  - `--humdiv` reads the PolyPhen column from the HumDiv matrix instead of HumVar; a cache without that matrix prints an empty column, as in Ensembl VEP. (#73)
  - `--nearest` lists every protein-coding transcript of a shared start site in stable-id order, where Ensembl VEP prints the one its tree search returns; a divergence page states it. (#73)
- Output formats, columns and headers
  - Every output format names the vep-rs version: `## vep-rs <version>` in the default and tab headers, `##vep-rs=<version>` in the VCF header, `vep_rs_version` in every JSON object. (#73)
  - `--fields` replaces the tab columns and the VCF CSQ subfields with the list given, a name nothing produces printing `-` under a `?` description, as Ensembl VEP does. (#73)
  - `--keep_csq` keeps an input VCF's existing annotation, header lines and INFO entry beside the new one; `--compress_output gzip|bgzip` writes the output through an in-process gzip or BGZF encoder. (#73)
  - `--terms` (`-t`), `--no_escape`, `--uploaded_allele` and `--ambiguity` do what Ensembl VEP's do: display or NCBI vocabularies, unescaped HGVSp and symbols, the uploaded allele string, the IUPAC code. (#73)
  - VCF output from Ensembl, HGVS, SPDI, region or ID input writes the line Ensembl VEP builds: an indel anchored one base back, reverse-strand alleles complemented, class words as `<DEL>`. (#73)
  - An empty `CSQ` entry in an input INFO (`CSQ=;DP=4`) is removed alone and `DP=4` kept, where Ensembl VEP's substitution takes the following field; a divergence page states it. (#73)
  - An empty entry under the `--vcf_info_field` name (`ANN=;CSQ=...`) is removed alone, where Ensembl VEP's substitution takes the next field; the empty-`CSQ` divergence page states it. (#73)
- Picking, filtering and summaries
  - `--transcript_filter` honours Ensembl VEP's filter grammar: every operator, parenthesised groups, `not`, `in` lists and files, repeated flags joined by `and`, and `is_child` over a Sequence Ontology table. (#73)
  - `--flag_pick`, `--flag_pick_allele` and `--flag_pick_allele_gene` mark once over a `--minimal` split record, where Ensembl VEP's first pass marks each split feature; a divergence page states it. (#73)
  - Under `--minimal --no_intergenic` vep-rs keeps the transcript row of a multi-allelic record's genic allele, where Ensembl VEP drops the whole record; an intended divergence page states it. (#73)
  - An intended divergence page states where `--transcript_filter` differs from Ensembl VEP: an abbreviated field on the first transcript, an ambiguous prefix, the regex dialect and `is_child`. (#73)
- Per-sample genotype modes
  - `--individual all|<list>` annotates a VCF record once per named sample, with `IND` and `ZYG` on every row, one VCF line and one JSON object per individual. (#73)
  - `--process_ref_homs` adds the reference allele's row to every feature of every record, numbered 0 and first, with the terms an allele equal to the reference takes. (#73)
  - `--individual_zyg all|<list>` annotates a VCF record once with every named sample's zygosity in one `ZYG` list, the reference allele alone for a single reference-homozygous sample. (#73)
- Command line, configuration and errors
  - `--grch37`, `-o` repeated (last wins with Ensembl VEP's warning), `STDOUT` in any case, and `docs/cli-reference.md` generated from the command line, the option tables and the error table. (#73)
  - `--safe` makes an unknown `--plugin` name or a missing plugins directory an error instead of a warning and skips the incompatible-pair check. (#73)
  - Read the `--config` file, `<dir>/vep.ini` and `VEP_*` variables in Ensembl VEP's priority, the command line winning. (#73)
  - Long options are read in any case, as a unique prefix and after one dash, as Getopt::Long does; an Ensembl VEP option this version lacks is refused by name. (#73)
  - `--input_data` (alias `--id`) replaces the input file; BCF and CAid input are refused by name; the `## VEP command-line:` header names each option by its resolved name. (#73)
  - `--decompression_threads` takes `auto` (opt-in, the default unchanged), reading an indexed BGZF input's record shape; hidden `--resolved_config` prints the resolved knobs and that shape before any cache opens. (#73)
- Warnings and skipped variants
  - `--warning_file` names the warnings file (default `<output_file>_warnings.txt`); every warning is written there and echoed on standard error unless `--quiet`, so `STDERR` carries it twice, as Ensembl VEP does. (#73)
  - `--skipped_variants_file` lists every record the run skipped in Ensembl VEP's layout; an existing file is refused before the run without `--force_overwrite`. (#73)
  - The `--skipped_variants_file` header names standard input `STDIN` in any spelling and under `--input_data`, and the record column is 40 characters; a divergence page states both against Ensembl VEP's. (#73)
- Plugins
  - Tabix-indexed data opens through a `.csi` as well as a `.tbi`, in htslib's order, and a record whose span reaches a queried position is returned as `tabix` returns it. (#73)
  - The plugin-facing variant JSON carries `input_chr`, the input's spelling, on a record renamed from `M` to `MT`; absent otherwise. (#73)
- Documentation, test corpora and comparison tools
  - Golden-corpus manifests carry `covers` (the parity-matrix rows a corpus exercises) and `flag_sets` (further runs under other flags, inputs or stdin), and the golden tests run every set. (#73)
  - `check_parity_matrix.py` accepts a `scripts/` pytest function as a process row's `test1`, requires its `perl_reference` to be `-`, and refuses a `covers` name that is not a row. (#73)
  - `scripts/concordance/compare_vep_formats.py` compares every column of every row against Ensembl VEP in the default, tab, VCF, JSON and Parquet formats, with `--self-parity` and `--warnings` modes, replacing `compare_output_formats.py`. (#73)
  - `tools/vep-fields-compare` builds the per-column field comparator in Rust, byte-identical to `compare_vep_fields.py`, and `run_clone_measurement.sh` runs it when `VEP_FIELDS_COMPARATOR` names it. (#73)
  - `generate_reference_output.sh` takes `--flags`, `--format default|tab|vcf|json|all`, `--stats` and `--perl-image-digest`, and `--corpus` writes a golden corpus's four expected outputs and `provenance.json`. (#73)
  - `scripts/concordance/feature_oracle.sh` runs a case table's comparisons, Ensembl VEP against candidate `vep` binaries in every output format, and writes per-case reports and a `status.tsv`. (#73)
  - `prune_json_cache.py --variations` cuts a corpus's variation slice as bgzip with a `.csi`, `--exclude-variation` drops records by name, and `info.json` keeps `variation_cols` only beside a slice. (#73)
  - `run_clone_measurement.sh --vep-extra-args "<flags>"` appends arguments to every vep-rs invocation through a wrapper and records the binary, its sha256 and every argument vector in `provenance.json`. (#73)
  - `build_golden_corpus.py select --keep-samples N` keeps `FORMAT` and the first N sample columns and samples records by the genotype classes they show, which `check` asserts are all covered. (#73)
  - Golden-corpus flag sets carry `unordered` (`build_golden_corpus.py flag-set --unordered`) for sets whose rows Ensembl VEP emits in hash order; the `GRCh38-pick-summary` corpus runs twelve pick and summary sets. (#73)
  - The `GRCh38-colocated-alt` golden corpus: span-only matches, null alleles, failed records and allele-specific ClinVar entries with Ensembl VEP's default-flag output. (#73)
  - The `GRCh38-colocated`, `GRCh37-colocated` and `GRCh38-colocated-alt` corpora gain flag sets running `--check_existing` with the co-location flags against Ensembl VEP 116.2's output. (#73)
  - The `GRCh38-colocated`, `GRCh37-colocated` and `GRCh38-colocated-alt` corpora run their flag sets with the frequency flags in every format, and `GRCh38-oldmaf` gains an `--old_maf` set. (#73)
  - The `GRCh38-input-compressed` corpus runs one input as BGZF (`.gz` and `.bgz`), plain gzip, standard input, CRLF and a Latin-1 header against one expected output. (#73)
  - The `GRCh38-checksvs` golden corpus: a slice of Ensembl's structural-variant dump with Ensembl VEP 116.2's `--check_svs` output from its database in four formats. (#73)
  - `prune_json_cache.py --drop-except <file>` keeps every `--drop` key on the transcripts the file names, so a corpus carries matrices and protein features for its own transcripts only. (#73)
  - The `GRCh38-freqfilter` and `GRCh38-freqfilter-include` corpora compare both engines under the exclude, include and `--filter_common` forms in every output format. (#73)
  - The `GRCh38-input-ensembl` corpus compares Ensembl-format input, tab- and space-delimited, and a space-delimited VCF against Ensembl VEP 116.2's output. (#73)
  - The `GRCh38-input-spdi` corpus compares SPDI input, under `--format spdi` and under format detection, against Ensembl VEP 116.2's database-mode output. (#73)
  - The `GRCh38-input-hgvs` corpus compares HGVS g., m., c., n. and p. input with Ensembl VEP 116.2's database-mode output in the default, tab, VCF and JSON formats. (#73)
  - The `GRCh38-everything` golden corpus compares `--everything` with Ensembl VEP 116.2 on 204 chromosome 21 records; the reference's regulatory rows are masked, since vep-rs annotates no regulatory feature. (#73)
  - Golden-corpus manifests carry `masks`, and `compare_vep_formats.py --masks regulatory_feature_rows` and `build_golden_corpus.py classify` apply the same class: the reference's RegulatoryFeature and MotifFeature rows are set aside. (#73)
  - The GRCh38-native-builder corpus runs its own native variation files under `--check_existing --exclude_null_alleles --af --af_1kg --clinvar_somatic_classification` against Ensembl VEP 116.2's output. (#73)
  - `prune_json_cache.py --variation-info` writes a variation slice's own `variation_cols` and `var_type` and refuses a declaration whose width differs from the slice's records. (#73)
  - The `GRCh38-vcf-check-ref` golden corpus: reference checks, skipped records, `--dont_skip`, `--chr`, `--lookup_ref` and `--allow_non_variant` against Ensembl VEP 116.2's output. (#73)
  - The `GRCh38-input-id` corpus compares identifier input, under `--format id` and under detection, against Ensembl VEP 116.2 connected to the Ensembl database. (#73)
  - The `GRCh38-input-region` corpus compares region input, under `--format region` and under format detection, against Ensembl VEP 116.2's output. (#73)
  - `compare_vep_formats.py --masks order_only` compares `DOMAINS`, `FREQS` and `ZYG` as multisets, counting every order-only agreement per column, since Ensembl VEP prints none in a meaningful order. (#73)
  - The `GRCh38-vcf-edge-cases` golden corpus: `<DEL:*>`, reference-equal and repeated ALTs, `*`, the anchor chop, gVCF alleles, breakends, the incomplete deletion, chromosome spellings, POS 0, five-column, CRLF, Latin-1 and stdin twins. (#73)
  - `build_golden_corpus.py classify` keeps a zero-row record whose reference JSON object carries `most_severe_consequence` `?`, a feature with no allele to annotate, instead of refusing it. (#73)
  - Golden-corpus manifests may mask `DOMAINS`, `FREQS` and `ZYG` order: a manifest or flag set naming `order_only` under `masks` compares their items as multisets in the golden tests. (#73)
  - Golden corpora `GRCh38-vcf-minimal` (23 records under `--minimal --uploaded_allele`) and `GRCh38-vcf-tandem-repeats` (15 `<CNV:TR>` records without a FASTA, with one, and under `--max_sv_size 10`). (#73)
  - The `GRCh38-individual` corpus runs eight sample columns of the 1000 Genomes high-coverage chromosome 21 panel plain and under five per-sample flag sets, and `GRCh38-vcf-minimal` gains a reference-row flag set. (#73)
  - The `GRCh38-individual-zyg` golden corpus runs eight sample columns of the 1000 Genomes high-coverage chromosome 21 panel under four `--individual_zyg` flag sets, the `ZYG` list compared as a set. (#73)
  - The `GRCh38-vcf-check-ref` corpus gains two `--use_transcript_ref` flag sets, plain and with `--process_ref_homs`, compared with Ensembl VEP 116.2 in tab, VCF and JSON. (#73)
  - The `GRCh38-minor-allele-columns` golden corpus: a fabricated variation slice declaring `minor_allele` and `minor_allele_freq`, compared with Ensembl VEP 116.2 under `--check_existing --af_1kg` in every format. (#73)
- Dependencies and tooling
  - Container image for linux/arm64 and linux/amd64 from generic-CPU archives; the measured binaries ship as the `X.Y.Z-neoverse-v1` and `X.Y.Z-x86-64-v3` tags. (#73, #24)
  - Container image carries `bgzip` and `tabix` from Debian's `tabix` package, so `--compress_output bgzip` VCF output can be indexed in the container. (#73, #68)

### Changed

- Input formats and records
  - An Ensembl-format allele given on the reverse strand is analysed as its forward-strand form, so its terms, codons and HGVS notations match the forward spelling. (#73, #63, #64)
  - The input format is detected from the first data line as Ensembl VEP detects it, and compression from the file's bytes (BGZF or plain gzip, whatever the name). (#73)
  - `-i STDIN` needs `--format`; a missing input, an undetectable first line, a region input without FASTA and an HGVS `r.` line are refused with Ensembl VEP's texts. (#73)
  - VCF input is read as bytes: a non-UTF-8 header line or sample name reaches the VCF output unchanged, and a data line with such bytes is annotated. (#73)
  - Empty and whitespace-only VCF lines are skipped and warned about in the words of Ensembl VEP's `skip_empty_lines`; the record numbering does not count them. (#73)
  - Ensembl-format input reads its strand as Ensembl VEP does: `-` and `-1` are reverse, any other token or an absent column forward, so four-column lines are accepted. (#73)
  - An Ensembl-format allele column without a `/` (`DUP`, `DEL`, `INS`, `INV`) is a structural variant over the line's span; an unknown word is warned about and annotated as written. (#73)
  - Ensembl-format alleles are upper-cased and a bi-allelic indel written un-minimised (`AC/A`) is trimmed to Ensembl VEP's representation before annotation. (#73)
  - Empty lines in an Ensembl-format input are skipped with the `skip_empty_lines` warning, where Ensembl VEP's parser dies on one; a skipped Ensembl-format record's warning quotes its line. (#73)
  - `--dont_skip` keeps a record that fails validation and annotates it, as Ensembl VEP does; a malformed line is skipped with its warning and the run completes with exit 0. (#73)
  - Lower-case alleles in an Ensembl-format record are upper-cased before annotation, every allele of a multi-allelic string included, under the name built from the line as written. (#73)
  - Region input is matched against Ensembl VEP's region expression: `start-end` is required, the allele may follow a `:`, and the strand token is `:1`, `:+1` or `:-1`. (#73)
  - A region allele outside `ACGTN-` (`DUP`, `DEL`, `INS:ME:ALU`) names a structural variant, and a line outside the expression is skipped with a warning. (#73)
  - A region record is named by its line as written in `Uploaded_variation` and the JSON `id`, as Ensembl VEP names it. (#73)
  - VCF alleles follow Ensembl VEP's overlap-allele rules: an ALT equal to the reference or a repeat takes no row and no `ALLELE_NUM`, `<DEL:*>` is numbered without a row. (#73)
  - A `<DEL>` with neither `END` nor `SVLEN`, or whose `END` is at or before `POS`, is warned `deletion looks incomplete` and dropped as in Ensembl VEP; `--dont_skip` keeps it. (#73)
  - A VCF record of five to seven columns is read with `.` for missing QUAL, FILTER and INFO, as Ensembl VEP reads it; fewer than five columns are refused. (#73)
  - A record at POS 0 keeps its coordinates in range: `A G` annotates at `21:0`, `A AT` is skipped with Ensembl VEP's insertion-coordinate warning, an empty REF the same. (#73)
  - `--minimal` trims every bi-allelic record from both ends and splits a multi-allelic record whose alleles trim, rejoined as Ensembl VEP rejoins it; it switches on `--allele_number`. (#73)
  - A `<CNV:TR>` record with a FASTA expands to explicit alleles as Ensembl VEP expands it; without one it stays the symbolic gain over the run with Ensembl VEP's warning. (#73)
  - An alternate allele equal to the reference, looked up or written, writes no row and takes no number; `FREQS` lists a repeated alternate allele once. (#73)
  - `--max_sv_size` measures every record as Ensembl VEP does; a negative tandem length is refused as `start > end+1`; a structural insertion past its own flank writes no row. (#73)
  - The region, `--chr` and HGVS transcript pages state an unsupported structural word's run end, a chromosome literally named `1-5`, and the four transcript readings, each with a test. (#73)
  - A VCF record's INFO column is split into its fields once per record, and every lookup on the record reads that split rather than scanning the column again. (#73)
  - The `--dont_skip` page states that a `#` line in region, Ensembl or HGVS input counts toward the line numbers, a bare one an empty line. (#73)
  - The Defect K, Defect C and structural-variant pages give `SVLEN=NaN` its default reading: dropped as incomplete, the two-base interval only under `--dont_skip`. (#73)
- HGVS, SPDI and identifiers
  - `--hgvsg` writes each row's HGVS genomic notation as Ensembl VEP's `hgvs_genomic` does; it needs a reference FASTA and exits with Ensembl VEP's error without one. (#73, #62)
  - `--hgvs` or `--hgvsg` without a FASTA stops with Ensembl VEP's message; `--everything` drops HGVS with its `INFO` line; `--shift_3prime 1` without one is refused. (#73)
  - `--format spdi` without `--fasta` is refused with `Cannot use SPDI format in offline mode without a FASTA file`, where every SPDI input was refused before. (#73)
  - HGVS genomic input (`g.`, `m.`) resolves against `--fasta` as Ensembl VEP's database does: `del`, `dup`, `inv` take genome bases, `delins` is trimmed, `ins` swaps positions. (#73)
  - An HGVS notation whose reference allele the genome contradicts, a repeat, an unknown shape or a variant `validate_vf` refuses is skipped with Ensembl VEP's warning. (#73)
  - HGVS input needs `--fasta`, is named by the notation in `Uploaded_variation` and the JSON `id`, and prints an accession or `chr` spelling under the FASTA's chromosome name. (#73)
  - `--shift_3prime` takes a value, `0` or `1`, as Ensembl VEP's option does; the bare switch of 0.4.0 is refused with the value list, exit 2. (#73)
  - The SPDI page and the `--format` help state that a version-less accession in an HGVS genomic notation is skipped, where Ensembl VEP resolves it as a synonym prefix. (#73)
  - The region page's `UPLOADED_ALLELE` section covers HGVS, SPDI and identifier records, whose second `--minimal` pass prints the once-trimmed string in Ensembl VEP (`TCA>TG`: `CA/G`). (#73)
  - The `--format` help scopes the any-case whole-name match to SPDI and HGVS input. (#73)
  - The regulatory page states that a regulatory-only record's intergenic `SPDI` names the chromosome as the cache does in Ensembl VEP and as the FASTA does here. (#73)
  - The output page and the `--spdi` help state that a record whose reference allele ends outside `ACGTN-` (`AR/CA`) gets no `SPDI`, where Ensembl VEP's `spdi_genomic` ends the run. (#73)
  - The HGVS divergence page states that a notation every core group refuses carries the first group's message alone, where Ensembl VEP joins every group's, with tests. (#73)
- Caches, reference FASTA and the cache tools
  - `--dir_cache`, `--species`, `--cache_version`, `--refseq` and `--merged` select the cache in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout; a run with no cache exits 1 before reading input. (#73)
  - `--assembly` must equal the cache directory's and `info.json`'s assembly (`hg19` and `hg38` read as `GRCh37` and `GRCh38`); the cache's assembly feeds headers and plugins. (#73)
  - Converted caches carry the Ensembl VEP variation files, about 23.5 GB more than 0.4.0's on GRCh38 and 21.8 on GRCh37; the `minor_allele`, `minor_allele_freq`, `AA` and `EA` fields are gone. (#73)
  - A cache converted by an earlier release needs the Ensembl VEP cache's `all_vars.gz` and `.csi` files copied into `variations/<chr>/`, or a new conversion. (#73)
  - `vep-cache-builder` reads the GFF3's `ccdsid`, `transcript_support_level`, `tag` codes, gene `version` and CDS `protein_id`, keeps every non-coding gene and drops `artifact` and readthrough transcripts. (#73)
  - A native cache carries MANE, TSL and APPRIS as transcript attributes alone, as the Ensembl cache does, and `is_canonical` only on the canonical transcript. (#73)
  - Native variation records take Ensembl VEP's form: the export's alleles kept at their position, an insertion's start past its end, `phenotype_or_disease` from the evidence flag. (#73)
  - `vep-cache-converter` checks `info.json` (`species`, `assembly`, its list keys) and copies `chr_synonyms.txt` beside it. (#73)
  - A plain-gzip `--fasta` is refused with the `bgzip` and `samtools faidx` remedy; the cache directory's first indexed `.fa.gz` counts as its FASTA too. (#73)
  - Native variation records take the variation database's form: `TSA=insertion` records two bases right of the VCF anchor, COSMIC and HGMD-PUBLIC variants as `COSMIC_MUTATION`/`HGMD_MUTATION`, `clin_sig` from the `CLIN_<term>` flags. (#73)
  - A native cache built by an earlier release holds JSON variation files nothing reads; build it again with `vep-cache-builder` 0.5.0 for co-located variants and the 1000 Genomes columns. (#73)
  - A native cache lists `cds_start_NF` before `cds_end_NF` on every transcript carrying both, where the Ensembl cache's order puts `cds_end_NF` first on seven GRCh38 transcripts. (#73)
  - The cache page states the `--fasta_dir` tie order, the assembly that arms the pseudoautosomal redirect, a non-UTF-8 `--synonyms` file and a two-hop synonym admitted after a miss, with tests. (#73)
  - The cache page names the FASTA record aliases vep-rs builds in (`MT`, `M`, `chrM`, `chrMT`, `chr` variants), which Ensembl VEP lacks, and the warning neither engine prints. (#73)
  - A cache holding no transcript shard has a divergence section and index row: refused here, run to the end by Ensembl VEP with every record `intergenic_variant`. (#73)
  - The cache page tabulates `--hgvs` over a chromosome the FASTA lacks: Ensembl VEP prints `N` as the reference base, vep-rs no HGVSc. (#73)
  - The cache page states that an insertion on a contig the FASTA lacks passes `--check_ref` and `--lookup_ref` in both engines. (#73)
  - The FASTA refusal table lists the `.gzi` that does not describe its file, the third refusal at open. (#73)
  - The compatibility page and the cache page's FASTA table state that `--fasta 0` or `--fasta ''` is a path here, where Ensembl VEP reads the flag as unset. (#73)
  - The cache page states the `--fasta_dir` filter reads the cache's assembly, where Ensembl VEP under `--offline --full_cache_dir` without `--assembly` or `--grch37` filters on the species alone. (#73)
- Co-located variants and frequencies
  - Under `--check_existing` a cache without variation data warns once, in the warnings file and on standard error unless `--quiet`, that `Existing_variation` and the frequency columns stay empty; Ensembl VEP annotates silently. (#73)
  - `--check_existing` matches known variants by allele after trimming from either end, strand-aware and span-exact, in Ensembl VEP's order; structural variants are not co-located. (#73)
  - A co-located variant matched by allele prints only on the rows of the alleles it matched, so a multi-allelic record's rows differ as Ensembl VEP's do. (#73)
  - `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED` print as Ensembl VEP does: allele-specific ClinVar entries by default, one flag per variant with all-zero lists dropped, duplicates kept. (#73)
  - JSON `colocated_variants` entries carry every set cache column as Ensembl VEP writes them, frequencies nested per allele; the record-level `existing_variation` array is gone. (#73)
  - `AF` prints with four decimals and the frequency columns report the record allele matched to the row's own allele, merged across co-located variants, as Ensembl VEP does. (#73)
  - `MAX_AF` prints its maximum over every frequency group, combined populations excluded, and `MAX_AF_POPS` lists every tied population in Ensembl VEP's order. (#73)
  - `--freq_pop 1KG_ALL` and `--filter_common` read the cache's `AF` column, where Ensembl VEP 116.2 finds no 1000 Genomes data and filters nothing. (#73)
  - The reach of a `--check_frequency` drop's rejoin reset under `--minimal` is stated on each side: a `--buffer_size` batch here, Ensembl VEP's chromosome buffer or fork chunk. (#73)
- Transcript attributes, predictions and consequences
  - `--distance up,0` reads the downstream distance as the upstream one, as Ensembl VEP does; a value that is not one or two integers exits with an error. (#73)
  - MANE, APPRIS, TSL, CCDS, RefSeq, UniProt, GENE_PHENO and miRNA print Ensembl VEP's values: accessions, APPRIS principal and alternative grades, digits, lists, no `-` sentinel. (#73)
  - `--sift` and `--polyphen` predict every single-residue substitution in the coding sequence, start codon changes included; a cache naming no version of the tool stops the run. (#73)
  - `--sift s` and `--polyphen s` print the bare score and set `sift_score` alone in JSON, a whole score an integer; JSON `domains` follow Ensembl VEP's `(\w+):(\w+)` re-split. (#73)
  - The consequence path reads biotype facts once per transcript, hoists the coding predicate once per allele, and keeps differing regions and the term list inline. (#73)
  - A reverse-strand record's alleles are read against each transcript's strand in every coding predicate; only its nucleotide alleles are reverse-complemented, a structural-variant word or symbolic allele staying as written. (#73)
  - `Codons` complements an IUPAC letter of an allele on the strand opposite its transcript's (`gcY`, not `gcR`), as Ensembl VEP's `reverse_comp` does. (#73)
  - `vep_core::codon::complement_base` complements the IUPAC letters as Ensembl's `reverse_comp` does, so a natively built cache's reverse-strand sequences agree with Ensembl's; `vep_effects::coding::reverse_complement_pub` is removed. (#73)
  - The `--minimal` split-allele page names the `miRNA` column: the default and tab outputs of Ensembl VEP compute a split allele's hairpin terms over the rejoined record. (#73)
- Output formats, columns and headers
  - `--output_format` is matched in any case, `ensembl` naming the default format, and an unknown value is an argument syntax error (exit 2) where 0.4.0 wrote the default format; plugin header descriptions reach the default, tab and VCF headers; `--fields` or `--compress_output` with Parquet is refused. (#73)
  - `--total_length` appends the cDNA, CDS and peptide lengths to the position columns, and `VARIANT_CLASS` is the whole record's class, so a multi-allelic record prints one. (#73)
  - A record on chromosome `M` prints `MT` only where Ensembl VEP renames it: the cache knows `MT` and neither `M` nor `chrM`. (#73)
  - The VCF output writes a VCF input's line back as the bytes it was read from, under any `--delimiter`, so a record that is not UTF-8 keeps every byte. (#73)
  - The default and tab renderers write rows as bytes: an allele's columns once per record, Extra values without intermediate strings, feature rows ordered by id rank. (#73)
- Picking, filtering and summaries
  - The APPRIS and TSL `--pick` tiers grade attribute values by Ensembl VEP's rules, so `principalM` and a zero carry no grade and `alternative` adds 10 only in lower case. (#73)
  - The pick family (`--pick`, `--pick_allele`, `--per_gene`, `--pick_allele_gene` and the `--flag_pick*` flags) acts per record as Ensembl VEP does: one row per record, allele, gene or both, intergenic rows too. (#73)
  - `--most_severe` and `--summary` print one row per record, `-` in the allele, feature and position columns, variant-level keys alone in Extra, and `intergenic_consequences` in JSON. (#73)
  - `--pick_order` is honoured with eleven criteria in the given order and Ensembl VEP's default list; `length` reads the translateable or spliced length, `ensembl` and `refseq` the transcript's source cache. (#73)
  - `--coding_only` keeps the rows whose allele affects a transcript's coding sequence and `--no_intergenic` drops the intergenic rows, as Ensembl VEP's output factory filters them, before the pick family. (#73)
  - `AnnotationResources` carries the compiled `--transcript_filter` and `FilterConfig` the `--coding_only` and `--no_intergenic` row filters, so an embedding caller sees the rows the `vep` binary writes. (#73)
- Command line, configuration and errors
  - Ensembl VEP's option rules apply as data: `--everything` is its exact set (no `ALLELE_NUM` or `HGVSg`), option sets override explicit values, incompatible pairs such as `--vcf --json` are refused. (#73)
  - Every refusal prints one `ERROR: [<component>] <message>` line and exits 1; argument syntax errors stay clap's and exit 2; an empty command line prints the help. (#73)
  - `--database`, `--host`, `--port`, `--user`, `--password`, `--registry`, `--db_version`, `--genomes`, `--is_multispecies`, `--lrg`, `--convert`, `--gp`, `--ambiguous_hgvs` and `--af_exac` are refused by name. (#73)
  - `--decompression_threads` above the host's logical CPUs is lowered to that count, so a huge value reads the input instead of exhausting threads. (#73)
  - The pipeline streams by default (the hidden `--pipeline_mode batched` restores whole-batch passes); the hidden `--pipeline_depth`, `--render_task_records`, `--writer_strategy` and `--output_buffer_size` knobs write the same bytes at every setting. (#73)
  - `--buffer_size 0` annotates every record, where Ensembl VEP writes the header alone; a negative or underscored value is an argument syntax error (exit 2). (#73)
  - The compatibility page states Ensembl VEP's reading of a repeatable `--config` key alone on its line and of `--sift 0`, `--polyphen 0` and `--nearest 0`, with tests. (#73)
  - The compatibility and cache pages state the `--synonyms` value is a file name as given, where Ensembl VEP's two-argument `open` reads `-` as standard input and `0` as unset. (#73)
- Warnings and skipped variants
  - A record Ensembl VEP skips (unknown chromosome, invalid alleles, unsupported or oversize SV, incomplete deletion) gets its `line N skipped` warning, N the record ordinal. (#73)
  - A skipped record's warning quotes its input line in every input format, as Ensembl VEP's `line N skipped (...)` does. (#73)
  - A skipped record of the Ensembl, HGVS, SPDI or region format is numbered counting the empty lines before it, as Ensembl VEP numbers it. (#73)
  - The `--warning_file` and `--skipped_variants_file` values are file names as given; the compatibility page states Ensembl VEP's `-`, whitespace, `>` and false-value readings. (#73)
  - The `--skipped_variants_file` help names Ensembl VEP's end-of-run check of an existing file beside the before-the-run refusal here. (#73)
- Plugins
  - `--dir_plugins` defaults to `<--dir>/Plugins`; a directory that does not exist is a warnings-file line when `--plugin` names a plugin. (#73)
  - Tabix-indexed VCF data reads a record's end from its INFO `END`, as `tabix -p vcf` indexes it, so a long record is found from every position it spans. (#73)
  - A `--dir_plugins` value holding a glob pattern or a `~user` prefix is read as one literal path, where Ensembl VEP expands it with `glob`. (#73)
- Documentation, test corpora and comparison tools
  - `compare_vep_formats.py --floors` counts a column under its JSON key in the JSON format and reports a record-level column as not applicable there instead of red. (#73)
  - The documentation site opens on a landing page and holds six sections of short pages, with every old page anchor redirecting to its new page. (#65)
  - The README is shortened to what vep-rs is, why, installation, quick start and a map of the site; its generated release sections are unchanged. (#65)
  - `build_golden_corpus.py flag-set` keeps a set's hand-written `masks` and `reference_warnings` when it rewrites the entry and validates both as the golden tests do. (#73)
  - The documentation gains "What vep-rs does not do", listing the Ensembl VEP features vep-rs lacks and what a run does with their flags. (#73)
  - Five divergence pages: a breakend mate beyond 5,000 bases under `--distance`, a `--check_svs` mate on another chromosome, two region notations, the `HGNC_ID` copy, a VCF without `#CHROM`. (#73)
  - The documented exclusions of `--check_ref`, `--fasta_dir`, identifier lists, skip warnings, `--spdi` under `--minimal`, the header time, `assembly_name` and `--gp` are stated on their pages, each with a test. (#73)
  - Seven more departures are stated: `--sift b,`, `-o` spellings, an unwritable `--warning_file`, a tab in a delimited field, HGVS spans past the end, repeated plugin keys, `--assembly` aliases. (#73)
- Dependencies and tooling
  - The container image and the CI test job install the DuckDB 2.0.0 CLI from its `duckdb-cli-linux-*.tar.gz` archive (1.x shipped a zip), both verified by sha256. (#73)
  - The release page groups its notes by area with nested bullets, marks breaking changes inline and may run to 400 lines; `CHANGELOG.md` nests the area groups under each kind. (#73)
  - The minimum supported Rust version is 1.89; the noodles crates move to 0.115, whose BGZF reader rejects a block that declares more than 65536 bytes. (#73)
  - `--output_format parquet` needs DuckDB 1.3.2 or newer, the first release whose Parquet `COPY` accepts `WRITE_BLOOM_FILTER`; the check compares the patch version. (#73)
- Per-sample genotype modes
  - The per-individual trimmed-record page covers every bi-allelic record `--minimal` trims past the anchor chop (`GA/GT`), where a `0/0` sample reads `HOM` in Ensembl VEP and `HOMREF` here. (#73)
  - A `#` line after the `#CHROM` line is ignored, where Ensembl VEP reads it as the header and names no sample; a new page states it. (#73)
  - The per-individual trimmed-record page names the record shape as the trigger: an indel trimmed past the anchor chop with or without `--minimal`. (#73)

### Removed

- Input formats and records
  - `vep_io::input::vcf::VcfParser` and its module; the `vep` binary's line parser is the one VCF reader. (#73)
- Caches, reference FASTA and the cache tools
  - The transcript `description` key, which Ensembl VEP's cache dumper deletes: the loader ignores it and `--transcript_filter` reads it as undefined on every cache, 0.4.0 native caches included. (#73)
- Command line, configuration and errors
  - `--transcript_index_impl` and the `coitrees` and sorted-array transcript indexes: UCSC binning is the one overlap index, and the flag is refused as an unknown argument. (#73)
  - `vep_core::DEFAULT_DISTANCE`, which nothing read; the `--distance` default of 5000 is the flag's own. (#73)
  - `vep_core::DEFAULT_BUFFER_SIZE` and `vep_core::DEFAULT_CACHE_REGION_SIZE`, which nothing read; the `--buffer_size` default and the cache's region width are declared where they are read. (#73)

### Fixed

- Input formats and records
  - A breakend mate's chromosome resolves from its written name against the cache, as in Ensembl VEP; a `chrM` mate without a synonym stays unresolved rather than reading as `MT`. (#73)
  - CAid input is refused with Ensembl VEP's own text, `Cannot use CAid format in offline mode`, under the `input` component. (#73)
  - A `*` or `<DEL:*>` allele of an Ensembl-format line writes no row and keeps its allele number; a line with no other allele is the feature with none. (#73)
  - `--allow_non_variant` reads a VCF record whose first ALT is `.` as a non-variant, dropped without the flag and carried bare with it, as Ensembl VEP does. (#73)
  - A VCF data line's leading blanks stay in its chromosome column, so the record is skipped as Ensembl VEP skips it, with its `Chromosome ... not found` warning. (#73)
  - `--delimiter` takes a regular expression and splits as Perl's `split` splits on it, capture groups included; a pattern that can match the empty string is refused. (#73)
  - A structural record whose SVLEN overflows the coordinate range is dropped after its oversize warning instead of being annotated to the chromosome's end. (#73)
  - A symbolic ALT whose letters name no class prints Ensembl VEP's type string in `Allele` and `VARIANT_CLASS`: `<CPX>` as written, `CPX` under `SVTYPE=CPX`. (#73)
  - A VCF whose `##fileformat` major version is not 4, or names no `<digits>.<digits>` version, is refused with Ensembl VEP's message before any output is written. (#73)
  - A VCF record with a `0` or empty ID and lower-case alleles, and each `--individual` copy, is named from its alleles as written, as Ensembl VEP does. (#73)
  - A record `--dont_skip` keeps after the `--chr` list or the chromosome test refused it is annotated with its alleles as written, as Ensembl VEP leaves them. (#73)
  - A VCF line's trailing empty columns are dropped, and a structural record echoes QUAL, FILTER and INFO as `.` where absent, empty or `0`, as Ensembl VEP stores them. (#73)
  - A tandem repeat or duplication whose SVLEN list ends in a comma reads it as Ensembl VEP does, so an all-negative `SVLEN=-300,` is refused as `start > end+1`. (#73)
  - A `*` or `.` beside another ALT under a true INFO/SVTYPE writes its own row of the record's class, numbered at its written place, as Ensembl VEP writes it. (#73)
  - Region input keeps the `N` placeholder as `UPLOADED_ALLELE` under `--minimal`, as Ensembl VEP prints it for every bi-allelic record. (#73)
  - `--lookup_ref` removes the `-` from every allele before the rejoin and drops a trailing empty allele, so a region allele `A-` annotates as `A`. (#73)
  - A structural record whose INFO/END is not a number is skipped with Ensembl VEP's `Invalid start '<POS+1>' or end '<text>' coordinate` line, after `deletion looks incomplete` when Perl's number for the text lies below the start. (#73)
  - A BGZF input followed by zero padding or a text trailer is read whole, as `gzip -dc` reads it, instead of failing with `invalid frame size`. (#73)
  - An Ensembl-format record `--dont_skip` keeps after the `--chr` list or the chromosome test refused it prints its alleles as written, as Ensembl VEP does. (#73)
  - An ALT beginning with `D` or `I` is chopped as an indel, Ensembl VEP's reading of the pileup spelling, so `DA DG` annotates as the SNV `A/G`. (#73)
  - An Ensembl-format record under `--minimal` is trimmed in its written case before upper-casing, so `AAAG/AGGg` keeps its `G/g` pair as Ensembl VEP keeps it. (#73)
  - Region input under `--minimal` measures `--max_sv_size` on every record's trimmed alleles, so a same-length allele with shared `N`s is kept as Ensembl VEP keeps it. (#73)
  - An HGVS notation naming no alternate bases (`DUP` under `--lookup_ref`, `ins` with no letters) annotates nothing, as Ensembl VEP's feature holds the reference alone. (#73)
  - An HGVS notation under `--minimal` is trimmed in its written case before upper-casing, so `delinsGGa` keeps its `A/a` pair as Ensembl VEP keeps it. (#73)
  - An HGVS notation prints `UPLOADED_ALLELE` as the pair it builds, before trim and upper-casing, as Ensembl VEP does when the first trim leaves equal lengths or a `-`. (#73)
  - An `--individual` copy of a lower-case record with a Perl-false ID prints its own alleles as written in the default and tab `Uploaded_variation`, as Ensembl VEP names it. (#73)
  - A deletion whose INFO/END is no coordinate warns `deletion looks incomplete` only when Perl's number for the text lies below the start, as Ensembl VEP compares it. (#73)
  - The `--individual` refusal names the header's samples as Ensembl VEP lists them, one per name in column order, a trailing tab naming no sample. (#73)
  - `-i ""` and `-i 0` are refused as `No input file given`, Ensembl VEP's text for a name Perl reads as false. (#73)
  - A BGZF stream followed by a plain gzip member is read whole, as `gzip -dc` reads every member for Ensembl VEP. (#73)
  - An empty ALT between two others keeps its empty part in the record's working string (`C/G//T`) and its `-` allele. (#73)
  - A VCF record's one ALT is anchor-chopped whatever it carries, and a `*` the minimisation leaves writes no row, as in Ensembl VEP. (#73)
  - A `*` the `--minimal` split's trim leaves in a multi-allelic record writes no row and keeps its allele number, as in Ensembl VEP. (#73)
  - A refused region line is quoted with its white space collapsed, so a kept carriage return reaches the warning as a space. (#73)
  - Region, Ensembl and HGVS input count a `#` line toward the line numbers in warnings, as Ensembl VEP does; `#` alone is an empty line in region and HGVS. (#73)
- HGVS, SPDI and identifiers
  - HGVSp and shifted codons of a moved insertion turn the allele as the record wrote it, unturned on a reverse-strand transcript when the shift exceeds its length. (#73)
  - `SHIFT_LENGTH` is written on sequence variants' transcript rows alone under `--shift_3prime 1 --shift_length`; a structural record's rows carry none, as in Ensembl VEP. (#73)
  - An HGVS protein notation whose reference contains `ENST` but names a translation (another species' id) resolves through the translation lookup, as Ensembl VEP's first lookup does. (#73)
  - Under `--shift_3prime 1` a coding row of a shifted insertion prints the inserted bases rotated along the shift as its `Allele`, as Ensembl VEP does. (#73)
  - Under `--shift_3prime 1` the HGVSp of a deletion shifted onto the stop codon reads its predicates on the shifted allele, `p.Ter687=` where an extension was printed. (#73)
  - A `--minimal` split allele shifts its HGVS notations as a bi-allelic record does; the record's mixed allele classes do not stop the shift. (#73)
  - `--spdi` expands a parenthesised tandem allele such as `(AC)3` before testing its characters, as Ensembl VEP does, so the row carries its notation. (#73)
  - SPDI input past the end of its sequence is skipped naming the span and the sequence's length, where Ensembl VEP annotates over `N` padding. (#73)
  - An HGVS `c.`, `n.` or `p.` input naming a transcript the cache cannot read is refused under the `cache` component, not `input`. (#73)
  - `HGVSg` of a trimmed allele of a `--minimal` multi-allelic record is the record-level notation Ensembl VEP prints, where it printed nothing. (#73)
  - An SPDI deletion of length 0 is skipped naming the empty span, not `Sequence name ... not valid`; Ensembl VEP keeps a `/-` feature that writes no row. (#73)
  - The `--format` help names the zero-length SPDI deletion skip and its warning beside the span-past-end skip. (#73)
- Caches, reference FASTA and the cache tools
  - A cache whose chromosome directories are symlinks loads; it was refused as holding no transcript shards. (#73)
  - `vep-cache-builder` writes a `TSA=indel` variation record with the export's alleles at its position, so a left-aligned input co-locates with a native cache as Ensembl VEP co-locates it. (#73)
  - The library `Annotator` resolves input chromosomes and breakend mates through the cache's synonyms, so a `chr21` row annotates as the `vep` binary annotates it. (#73)
  - `--synonyms` reads a line opening with whitespace as Ensembl VEP does: the reference is empty and the line's names are not linked to each other. (#73)
  - A chromosome directory of the cache holding no transcript shard annotates its records `intergenic_variant`, and a dot-prefixed directory is no chromosome, as Ensembl VEP lists them. (#73)
  - A reference read on `Y` inside a GRCh37 or GRCh38 pseudoautosomal region is served from `X`, as Ensembl VEP serves it, so `--check_ref` and HGVS see the bases. (#73)
  - A directory given as `--synonyms` or as the cache's `chr_synonyms.txt` is refused, where Ensembl VEP opens it and reads an empty table; the cache page states it. (#73)
  - A present but empty `assembly` in a cache's `info.json` is compared with `--assembly` and refused with the mismatch text. (#73)
  - A transcript shard the cache cannot read is refused under the `cache` component on every input path, not as a run failure. (#73)
  - `storable_to_json.pl` keeps a cache's `sift` or `polyphen` marker only beside its `source_<tool>` line, the one Ensembl VEP reads, so a converted cache refuses the flag where Ensembl VEP does. (#73)
- Co-located variants and frequencies
  - Under `--minimal`, every row of a split record lists the co-located variants and `SV` names of its first alternate allele's feature, as Ensembl VEP prints them. (#73)
  - A region file's empty `clin_sig_allele` or `clinical_impact` cell is a value, as Ensembl VEP reads it: `CLIN_SIG` is suppressed and `CLINVAR_SOMATIC_CLASSIFICATION` prints empty. (#73)
  - The JSON `frequencies` of a co-located variant are keyed by the record's own alleles, so a `--minimal` split allele's trimmed spelling finds none, as Ensembl VEP writes. (#73)
  - `CLIN_SIG` on a `--minimal` split row reads the ClinVar reference against the reference `REF_ALLELE` prints, not the allele's own trimmed one, as Ensembl VEP compares them. (#73)
  - JSON `colocated_variants` objects echo `minor_allele` and `minor_allele_freq` from a cache declaring them, as Ensembl VEP's JSON factory copies the record. (#73)
  - `--check_frequency` under `--minimal` judges each trimmed feature of a split record against its own co-located frequencies, as Ensembl VEP does. (#73)
  - Under `--freq_filter exclude` the `--process_ref_homs` reference row prints the record's rewritten allele string and class. (#73)
- Transcript attributes, predictions and consequences
  - `--mane_select` alone lists no MANE Plus Clinical transcript, UNIPARC and UNIPROT_ISOFORM carry the cache's accessions, and a VCF SYMBOL writes whitespace as `%20`. (#73)
  - `--distance` up to `18446744073709551615` annotates structural variants without overflowing; a debug build panicked and a release build wrapped the window past a transcript's end. (#73)
  - `OverlapPC` rounds a two-decimal tie as Ensembl VEP's `100 * (overlap / length)` rounds it: 14.37, not 14.38, for 23 of 160 bases. (#73)
  - A reference-length-mismatch record `--dont_skip` keeps is evaluated over its span, so one over a stop codon prints `stop_lost` without `frameshift_variant`. (#73)
- Output formats, columns and headers
  - The Parquet writer keeps its numeric casts out of DuckDB's CSV reader, so a malformed numeric token is NULL under DuckDB 2.0.0 as it was under 1.x. (#73)
  - An output file vep-rs cannot create, a Parquet intermediate or an empty name included, ends the run with Ensembl VEP's `ERROR: [output] Could not write to output file <file>`. (#73)
  - The tab header describes every `--fields` name in Ensembl VEP's description table; `ID`, `CELL_TYPE`, the ESP, ExAC and RefSeq keys and `BAM_EDIT` printed `?`. (#73)
  - `--tab --regulatory` without `--biotype` prints the `BIOTYPE` column as `-`, the first field of Ensembl VEP's regulatory group. (#73)
  - `--total_length` appends the CDS and peptide totals behind `-` for an insertion between a UTR and the CDS, as Ensembl VEP's coding predicate admits it. (#73)
  - `VARIANT_CLASS` classes the alleles `--check_frequency` kept, so a record whose dropped allele was its `-` prints `SNV`, as Ensembl VEP does. (#73)
  - `NEAREST` on a `--minimal` multi-allelic record is the first split allele's on every row of the default, tab and VCF formats, as in JSON and Ensembl VEP. (#73)
  - `--vcf_info_field` refuses a value that is not a VCF INFO key, and a removed INFO entry takes its closing `;` with it. (#73)
- Picking, filtering and summaries
  - `--pick` on a `--minimal` split ranks the alleles in Ensembl VEP's rejoin key order, so a repeated trimmable alternate allele keeps the row Ensembl VEP keeps. (#73)
  - `--transcript_filter` reads `uniprot_isoform`, the cache's UniProt isoform ids, matching any id of a list as Ensembl VEP does. (#73)
  - `--exclude_predicted` is accepted and appends Ensembl VEP's `not stable_id match ^X._` transcript filter after the `--transcript_filter` values; it was refused as an unexpected argument. (#73, #67)
  - `--coding_only` judges a `--minimal` split record's rows over the record's span in the default and tab outputs and prints the record's own `REF_ALLELE` where its first feature wrote none. (#73)
  - A repeated `--pick_order` keeps its last list, as Getopt::Long keeps the last value of a scalar option; the lists were concatenated. (#73)
  - `--summary` beside `--most_severe` under `--safe` lists every term, as Ensembl VEP's summary branch does; the most severe term alone was printed. (#73)
  - `--transcript_filter` admits an `is_child` value outside the Sequence Ontology and matches it as a regular expression first, as Ensembl VEP does; it was refused. (#73)
  - `--pick_order ""` is the empty order, keeping the first candidate as Ensembl VEP does; it was refused as a syntax error. (#73)
- Per-sample genotype modes
  - `--process_ref_homs` reads a reference row's codon from the transcript, as Ensembl VEP does, so a multi-base reference row whose bases mismatch the genome is `synonymous_variant`, not `missense_variant`. (#73)
  - Under `--safe --individual --minimal`, an individual's copy carrying one alternate allele of a multi-allelic record keeps the record's chopped alleles, as Ensembl VEP writes it. (#73)
  - Under `--process_ref_homs --shift_3prime 1` the reference row's positions, codons and coding terms are Ensembl VEP's; its non-coding terms are the shifted span's, an intended divergence. (#73)
  - Under `--individual`, a copy whose genotype names an ALT after a leading `.` is annotated as that allele, as Ensembl VEP annotates it. (#73)
  - A reference-homozygous `--individual` copy under `--lookup_ref --process_ref_homs` keeps its written reference as allele 1 beside the looked-up reference row, as Ensembl VEP does. (#73)
- Command line, configuration and errors
  - The `## VEP command-line:` header names `--config`, `vep.ini` and `VEP_*` flags, leaves default values out and prints an `=i` value as Getopt::Long stores it; a non-UTF-8 argument runs. (#73)
  - A trailing comma in a list flag (`--chr 21,`, `--distance 5000,`, `--pick_order canonical,`, `--individual S1,`) names no empty item, as Ensembl VEP's `split` reads it. (#73)
  - A bare `--` on the command line, and every token after it, is refused as an unexpected parameter with Ensembl VEP's text unless `--help` precedes it. (#73)
  - `--cache_version 0`, an empty or `0` `--dir_cache` or `--dir_plugins`, and the `--dir` default without `HOME` fall back to their defaults, as Ensembl VEP's `||` reads them. (#73)
  - `--help` and `--version` win over every other command-line token, as the `vep` script prints its usage before refusing extra parameters. (#73)
  - `--quiet` silences `--verbose` under `--safe`, as `status_msg` returns on `--quiet` before printing. (#73)
  - The command-line header prints `--full_cache_dir` under its own key from the command line, a `--config` file or `VEP_FULL_CACHE_DIR`. (#73)
  - The existing-output check reads the `-o` value as spelled, so a file named `STDOUT` refuses `-o STDOUT` without `--force_overwrite`. (#73)
  - `--input_data` given as empty or `0` beside `-i` reads the file and reports a skipped record as `variant skipped`, without a line number. (#73)
  - A `var_synonyms 0` from a `--config` file, `vep.ini` or `VEP_VAR_SYNONYMS` switches `check_existing` on with the notice, without synonyms, as Ensembl VEP reads a defined `0`. (#73)
  - `--help` after a flag Ensembl VEP lacks prints the usage, and after `--gff`, `--gtf` or `--bam` is that option's value, as Getopt::Long reads it. (#73)
  - A `VEP_*` variable whose value is not UTF-8 sets its flag with the replacement character for the byte, where it was dropped. (#73)
  - `-o 0` writes its warnings to `vep_warnings.txt`, the file Ensembl VEP names for an output Perl reads as false. (#73)
  - `-o -` writes the annotation to standard output, as Ensembl VEP's two-argument open does, except under `--compress_output`. (#73)
  - `--show_cache_info` prints the cache information before an unreadable `--synonyms` file is refused, as Ensembl VEP reads the file only in a run. (#73)
  - A skipped non-VCF line ending `\r\r\n` is quoted with its lone `\r` as a trailing space, as Ensembl VEP quotes it. (#73)
  - The warning for a repeated --output_file ends with a space before the newline, as Ensembl VEP's does. (#73)
  - `af_exac 1` from the `--config` file, `vep.ini` or `VEP_AF_EXAC=1` reaches the `--af_exac` refusal (exit 1) instead of being dropped silently. (#73)
  - The `--af_exac` refusal names Ensembl VEP release 107, the release that dropped the flag, in place of release 115. (#73)
- Warnings and skipped variants
  - A refused Ensembl-format line or identifier mapping is skipped with the reason alone; the warning opened with `VCF parse error:`. (#73)
  - `--warning_file STDERR` carries each warning twice unless `--quiet`, as `print $fh` then `warn` write it; an empty or `0` value is the default file. (#73)
  - Under `--individual`, the identical skip warning of a record's sample copies is written once, as Ensembl VEP's seen set keys the whole text. (#73)
  - The line a skip warning quotes is cut by bytes at the last character boundary at or before byte 46, where Ensembl VEP's `substr` splits the character. (#73)
- Plugins
  - A tabix query answered from several chunks returns every record once, and a line cut between chunks is read whole instead of losing a record or yielding a fragment. (#73)
  - A tabix index built for another data file is refused with both files named instead of stopping the process. (#73)
  - A reverse-strand Ensembl-format or region record is matched to CADD, dbNSFP, LoFTEE and every tabix-backed plugin file on its forward-strand alleles, as Ensembl VEP matches it. (#73)
  - A plugin whose constructor fails is dropped with Ensembl VEP's `Failed to instantiate plugin` warning and the run completes; `--safe` makes it the plugin component's error. (#73)
  - A `--plugin` argument's trailing comma passes no parameter, as Perl's `split` passes none; a middle empty field stays. (#73)
  - A `--plugin` library that does not open is dropped with Ensembl VEP's `Failed to compile plugin` warning, the plugin component's error under `--safe`; it ended the run. (#73)
  - Plugin-loading warnings follow the `--plugin` arguments' order, as Ensembl VEP's one loop writes them; built-in failures came first. (#73)
- Documentation, test corpora and comparison tools
  - `feature_oracle.sh` decompresses each `.gz` FASTA into its own directory, so two cases whose FASTAs share a file name run in one table. (#73)
  - The structural-variant comparator sets aside the literal rows a run with a FASTA writes for the real allele behind a leading `RN=0`, beside the symbolic rows. (#73)
  - `build_golden_corpus.py classify` numbers a flag set run over its own input by that file's records and keeps a record `--coding_only` or `--no_intergenic` left rowless on both sides. (#73)
  - The lone `.` ALT page states that vep-rs writes no warning for it, and the empty-lines page that Ensembl VEP's Ensembl-format parser dies on an empty line. (#73)
  - Defect N's page states that Ensembl VEP's anchor chop also numifies a sequence indel's non-numeric POS (`12abc` annotated at 21:13, `abc` refused as `'abd'`), where vep-rs refuses every form. (#73)
  - The divergence index gains the one-column `--fields` list whose row is `0`, at which Ensembl VEP's output loop ends and vep-rs writes every row. (#73)
  - The compatibility page's prefix example is `--bi` is `--biotype`; `--sy`, which it named, is ambiguous between `--symbol` and `--synonyms` in both engines (exit 2). (#73)
  - `--no_progress` and `--debug` are listed among the flags vep-rs refuses as unknown arguments (exit 2) on the page of what vep-rs does not do. (#73)
  - `compare_vep_formats.py` keys a JSON object by the first five columns of its `input`, so a structural variant whose QUAL `0` Ensembl VEP writes as `.` keeps its record. (#73)
  - `compare_vep_formats.py --unordered` compares a flag set's entries without their order, as the test harness does for the three sets Ensembl VEP emits in hash order. (#73)
- Security and robustness
  - The container image disables the two CA roots the eMudhra intermediates chain from (emSign ECC G3, Comodo AAA Services); cache downloads are checked against the rest. (#70)
  - Input files drive no unbounded allocation: line lengths, index sizes, cache sequence lengths and plugin-data counts are checked, and control bytes are escaped in tab and VCF output. (#73)
  - Inputs that crashed the tools or escaped the output directory are refused or read: unterminated symbolic ALTs, oversized BGZF blocks, mixed exon ranks, path-escaping chromosome names, overlong FASTA indexes. (#73)
  - Plugin catalogues, cache-builder and converter inputs, variation files and FASTA indexes refuse a line past their line cap instead of reading it whole. (#73)
  - A transcript's reference allele read from the FASTA is not sized from a mapper's cDNA ends, which a malformed cache can set any distance apart. (#73)
  - The cache loader refuses a transcript whose mapper's coding end precedes its coding start. (#73)
  - A tabix data line or converter input whose zero-based start has no one-based form is skipped or refused instead of overflowing the shift. (#73)
- Dependencies and tooling
  - The release page's caption labels a cross-generation class pair the newest available class per architecture, says when ClinVar shares the whole genome's class, and refuses a missing record field. (#73)
  - Linux release archives run on glibc 2.17 or later, the floor read from each binary's ELF version needs and printed on the release page. (#73, #23)
  - The container image's documentation label points at the installation page's container section. (#73)

### Security

- Security and robustness
  - Cleared from the dependency tree: RUSTSEC-2026-0346, RUSTSEC-2026-0347, RUSTSEC-2026-0348 and RUSTSEC-2026-0349 (rustls 0.23.46) and the RUSTSEC-2026-0253 unsoundness notice (lru 0.18.5 via mysql 28.0.4), none reachable from the `vep` binary. (#73)
  - Cleared from the dependency tree: RUSTSEC-2026-0382 (rustls-webpki 0.103.16, an optional `reqwest` TLS backend no binary compiles), not reachable from the `vep` binary. (#73)
  - Cleared from the dependency tree: RUSTSEC-2026-0386 (tokio 1.54.1, the `AsyncFd` constructors' I/O-safety notice), reached only through `bigtools` and `reqwest`, which never construct one. (#73)

## [0.4.0] - 2026-10-16

### Added

- `--protein_version` appends the translation version to `ENSP`, `--mask_header_cache_path` masks the cache directory in the output headers, and `--regulatory_gff`, `--extended_promoters` and `--custom_suppress_filter` are refused, naming what each would do. (#44)
- The plugin-visible `InputVariant` JSON carries `input_line_name`, `record_so_term` and `allele_num` (the record's input-line name, class term and printed ALLELE_NUM), and `mate_id` only on a record naming a mate. (#44)
- Comparators take `--reference-release {115.2,116.2}` (default 116.2) with a registry of records Ensembl VEP drops or annotates differently; release pages tabulate whole-genome and structural-variant F1; documentation: https://natera-open-source.github.io/vep-rs/. (#43, #44, #50, #57)

### Changed

- Ensembl VEP 116.2 is the reference release: headers and banner read `v116.2`, `--cache_version` defaults to 116, `tests/golden/116/` replaces `tests/golden/115/`, and `docs/intended-divergences.md` is rewritten against it. (#38, #44, #45)
- Stop-codon consequences follow Ensembl VEP release 116: an edit beginning inside the stop codon is `stop_lost` or `stop_retained_variant` without `frameshift_variant`; indels decide `stop_retained_variant` on the genomic span. (#44)
- Breakend records follow Ensembl VEP 116: an `intergenic_variant` row per allele when no transcript is reached, own breakend `N.` without MATEID at ALLELE_NUM 0, `<BND>` mates named through CHR2. (#44)
- A multi-allelic structural variant takes one class from its joined ALTs and writes one row per ALT, each labelled with the record's class term, as Ensembl VEP 116 does. (#44)
- Records Ensembl VEP 116 drops before annotation are annotated, the divergence the golden manifests class `reference_skipped_record`, and `--max_sv_size -1` lifts the size limit in every output format. (#44)
- JSON `cdna_end`, `cds_end` and `protein_end` are omitted when the end is undefined (`445-?` writes `cdna_start` alone), and a structural record's `allele_string` carries every ALT as written (Ensembl VEP 116). (#44)
- `Uploaded_variation` reads a VCF record's first `;`-joined ID, and a record whose ID is `.`, empty, `0` or absent is named from its input line (Ensembl VEP 116). (#44, #51)
- `--hgvs` follows Ensembl VEP 116: a deletion past a transcript's edge written over the bases inside it, a stop-codon frameshift as an extension, a stop-codon indel by its peptides. (#44)

### Fixed

- Structural-variant records are routed, typed and spanned as Ensembl VEP's parser reads them, in one module shared by the command-line parser and the `vep-io` reader; no coordinate wraps. (#44)
- `<CNV:TR>` alleles follow Ensembl VEP 116: the length reads RUC before RB and sums the RN sequences, and a repeat longer than `--max_sv_size`, or any under `-1`, stays symbolic. (#44)
- A structural deletion or copy-number loss over the stop-codon window takes `stop_lost` unless its first base in transcript order lies in the incomplete terminal codon (Ensembl VEP 116). (#58)
- `--protein` fills ENSP on every structural-variant row of a translated transcript, and the `## VEP command-line` header keeps `--max_sv_size -1` and the `key=` prefix of a masked `--custom` path. (#44)
- JSON `id` of a record whose first ID is empty or `0` is the generated name, and the header cache path without `--json_cache` follows Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout. (#51)
- The concordance harness passes `run_splice_predictions:0` to Perl LoFTEE and keys its Perl output cache on plugin arguments, and `--plugin LoFTEE` warns when asked for that unimplemented splice layer. (#47, #49)

## [0.3.2] - 2026-10-06

### Changed

- `--fork 0` (the default) uses every logical CPU up to 32 instead of every
  logical CPU: one reader feeds the worker pool, so the wall time stops falling
  past about 32 workers and rises on wider hosts (gnomAD v4.1 chromosome 2 took
  1.8 times longer at 192 workers than at 32 on a 192-vCPU host). An explicit
  `--fork N` is never capped.
- A batch is annotated, filtered, rendered and freed in one pass over the worker
  pool, each finished variant freed by the thread that annotated it and the
  task size following the pool's thread count, instead of five pool barriers per
  batch with the coordinator freeing every variant serially between batches. The
  output is byte-identical. Together with the libdeflate inflate below, measured
  against 0.3.1 under the same whole-genome protocol (the median of 20 machines
  per chromosome at 16 threads, summed over the chromosomes), the whole of gnomAD
  v4.1 falls from 974 to 724 s on ARM Graviton4 and from 1,219 to 1,169 s on x86
  and the whole of gnomAD v2.1.1 from 586 to 282 s and from 835 to 476 s, with
  CPU time down 19 to 44 percent on those inputs and peak memory at most 2.4 GB
  on any chromosome against 4.9 GB.
- bgzip blocks inflate through libdeflate (`noodles-bgzf` with its `libdeflate`
  feature), and the cache builder's MySQL client takes its Rust deflate backend
  so that no crate of the workspace links C zlib: the single-inflater floor on
  gnomAD v4.1 chromosome 21 falls from 73.8 s to 34.5 s.

## [0.3.1] - 2026-10-02

### Added

- The SNP/indel comparator applies the transcript-selection class of
  `docs/intended-divergences.md` (section 3) to the records of a whole-genome
  input that Ensembl VEP marks skipped, under the structural-variant
  comparator's two guards, and reports it as `sv_transcript_selection`; the
  scope test is defined once and shared by the two comparators.

### Fixed

Each entry names the Ensembl VEP mechanism vep-rs follows.

- A symbolic insertion whose `END` is its position or whose `SVLEN` is 0 is the
  insertion between the two bases, `start = POS + 1, end = POS`, as ensembl-io's
  `get_end` reads it: its Location names both flanks, its predicates read the
  pair as `_intron_effects` reads a sequence insertion (both flanks inside a
  region, the two insertion special cases at the intron edges, the UTR special
  cases at the coding region's edges, `within_cdna` through `map_insert`), and
  the row carries no `OverlapBP` or `OverlapPC` (an overlap length of 0). A
  symbolic record with neither field stays the single base at `POS + 1`, which
  is also ensembl-io's reading.
- UTR term on an insertion at the edge of the coding region: `_before_coding`
  and `_after_coding` special-case an insertion whose start is the coding
  region start or whose end is the coding region end, so an insertion between
  the last base of an intron and the first base of a coding exon (or between
  the coding region's last base and the intron after it) is `5_prime_UTR_variant` or
  `3_prime_UTR_variant` by strand beside its splice term, where it carried the
  splice term alone.
- A frameshift intron on the UTR side of the coding region: `within_cds`
  requires the span to overlap the coding region, and `_bvfo_preds` sets `utr`
  rather than `coding` for a span that does not, so a variant inside an intron
  of twelve bases or fewer that lies before the coding region is
  `5_prime_UTR_variant` (or `3_prime_UTR_variant` after it) rather than
  `coding_sequence_variant`.
- Structural `stop_lost` behind the `coding` pre-predicate: the term carries
  `include {coding => 1}`, and `_bvfo_preds` sets `coding` only for a span that
  overlaps the coding region and an exon, so a deletion inside the intron that
  splits a stop codon is `non_coding` and carries no `stop_lost`, although the
  codon window drawn on genomic coordinates reaches into that intron.
- No UTR term on a structural span that starts on the first base of a
  transcript whose coding region starts there (`cds_start_NF`, and likewise at
  the last base with `cds_end_NF`): `_before_coding` is
  `overlap(start, end, transcript start, coding_region_start - 1)`, an inverted
  window that admits only a span starting before the transcript.

## [0.3.0] - 2026-10-01

### Added

- A container image, `ghcr.io/natera-open-source/vep-rs`, for each published
  release from 0.2.0: the release's attested x86_64 Linux binaries and the
  `duckdb` CLI on Debian 13, smoke-tested on the GRCh37 corpus before it is
  pushed, with a build-provenance attestation on the pushed digest
  (`.github/workflows/docker.yml`).
- `docs/intended-divergences.md`: every class the adjusted concordance sets
  aside, each with the record, both engines' rows, the mechanism in Ensembl
  VEP's source, and a command that reproduces it.
- Structural-variant comparator masks for the second side of two classes and
  for three more: the `intergenic_variant` row vep-rs writes for a record
  whose VEP tuples all name transcripts on another chromosome
  (`filter_cross_chromosome_orphan_intergenic`); the mate allele of a
  breakend above `--max_sv_size` that VEP never wrote, gated on the mate
  lying within 5 kb of the transcript span
  (`filter_giant_breakend_mate_divergences`); the SNP/indel registry's
  excluding rules applied to structural pairs (`filter_registry_swap_pairs`);
  VEP's batch-dependent Transcript tuples on a `<NON_REF>` record
  (`filter_non_ref_batch_divergences`); and a mate-side row VEP evaluated at
  the local coordinate (`filter_breakend_mate_local_read_pairs`), with
  `scripts/validation/breakend_mate_context.py` deriving the mate-coordinate
  and local-coordinate readings from the JSON cache. The transcript-selection
  mask's unsupported-type scope names `<CPX>` and `<CTX>`.
- The golden classifier names `breakend_mate_local_read`, computed by the
  same derivation over the corpus's own pruned cache.
- The comparators write the adjusted residual: `discordant_open.tsv` beside
  `discordant.tsv` for the SNP and indel sets and `discordant_adjusted.tsv`
  for the structural-variant sets, each the one-sided rows every mask leaves
  in place. Per-class mode caches the per-term VEP totals of a reference
  output (`--vep-totals-cache`).

### Changed

- JSON cache format: each transcript's `variation_effect_feature_cache`
  carries `seq_edits`, the translation's `Bio::EnsEMBL::SeqEdit` list
  (`start`, `end` and `alt_seq` in protein coordinates, with `code`, `name`
  and `length_diff`), which `scripts/data/storable_to_json.pl` writes from
  the Storable cache. vep-rs applies the edits overlapping a variant's
  translation span to the reference allele's peptide, as
  `TranscriptVariationAllele::peptide` does. A cache converted without the
  key still loads: no edit is applied, and a variant at codon 1 of a
  transcript with an alternative initiation codon is called from the cached
  peptide. Convert such a cache again to get VEP's calls there.
- The deletion-shaped evaluation of a same-chromosome breakend span borrows
  the variant for each overlapping transcript instead of copying it with the
  consequences accumulated so far, which grew with the square of the
  transcript count. Output is unchanged.

### Fixed

Each entry names the Ensembl VEP mechanism vep-rs follows.

- Start codon on an alternative initiation codon: with `seq_edits` in the
  cache, an `initial_met` or `amino_acid_sub` edit over a non-ATG start
  codon gives a reference `M` against the alternate allele's literal
  residue, so the variant is `start_lost` rather than `synonymous_variant`
  and `Amino_acids` prints the edited residue.
- Splice terms across a frameshift intron: `_intron_effects` skips an intron
  of twelve bases or fewer for every differing region that overlaps it, in
  both of its loops, so a deletion running from one exon across such an
  intron into the next carries no splice term; the coding predicates alone
  decide the row. The exonic splice windows apply the same skip, so a donor
  window that crosses such an intron yields no `splice_donor_5th_base_variant`.
- `transcript_ablation` requires Perl's `deletion` pre-predicate (the
  reference span longer than the ALT sequence) beside `complete_overlap`. A
  same-length substitution or a longer ALT covering a whole transcript takes
  the tier-3 predicates instead: the biotype's context term with the UTR,
  intron and splice terms of the span.
- `splice_polypyrimidine_tract_variant` on structural alleles: a structural
  allele whose span reaches the 15 bp window inside an intron's acceptor end
  carries the term unless the span touches an exon of the transcript, as
  Ensembl's predicate does for every BaseVariationFeature. A `<CNV:TR>`
  allele keeps its endpoint test.
- A paired breakend's bracket allele gets a row for a transcript only when
  the mate coordinate lies within `MAX_DISTANCE_FROM_TRANSCRIPT` (5 kb) of
  it (`StructuralVariationOverlap::_close_to_feature`); the local POS is
  carried by the single-breakend allele.
- A single-breakend allele above `--max_sv_size` that partially covers a
  coding transcript takes the region predicates like every other
  same-chromosome span; VEP requires `complete_overlap_feature` for
  `coding_transcript_variant` at any span.
- The context term of a transcript a structural allele engulfs follows the
  biotype, as VEP's include hashes do: `coding_transcript_variant` only for
  `protein_coding`, `NMD_transcript_variant` only for
  `nonsense_mediated_decay`, `non_coding_transcript_variant` only without a
  translation, and `intergenic_variant` for a translated transcript of any
  other biotype. One rule serves the inversion, copy-number, insertion and
  breakend arms.
- Structural insertions, mobile-element insertions and duplications whose
  span covers a start-codon base with both ends in exons emit `start_lost`,
  the genomic-overlap predicate the deletion arm applies.
  `start_retained_variant`, which VEP co-emits there without reading
  sequence, is not emitted.
- The mate-side Transcript row of a paired breakend describes the mate
  breakend at the mate coordinate: `feature_truncation` with the region term
  when the mate lies inside the transcript, `upstream_gene_variant` or
  `downstream_gene_variant` with the distance when it lies within the window
  outside, and no row beyond it. VEP evaluates every predicate but
  `feature_truncation` at the local coordinate there.
- A gVCF `<NON_REF>` record is a reference-confidence block, not an
  alternate allele. It is written as one `intergenic_variant` row spanning
  POS+1 to END with the allele as written and `IMPACT=MODIFIER`, with a
  warning per record; that is the row VEP writes when no other record in the
  batch loads the region. A `<NON_REF>` carrying a supported `INFO/SVTYPE`
  keeps that type.

## [0.2.0] - 2026-09-29

### Fixed

- HGVSp is Perl VEP's `hgvs_protein`, ported whole. In-frame insertions
  between codons, duplications, stop-loss extensions (`extTer{n}` and
  `extTer?`), insertions carrying a stop codon and deletions shifted at the
  peptide level print VEP's string where they printed nothing or a different
  one. The coding gate is read from the allele as annotated and the
  frameshift test from its 3'-shifted span, as VEP reads them, and the
  transcript peptide as `Transcript::translate` writes it (`Met1` for any
  start codon of the table). The one intended difference: on a variant that
  keeps the start codon intact and is reported as `start_retained_variant`,
  VEP prints `p.Met1?` from its own `start_lost` call while vep-rs prints
  the peptide change.
- HGVSc is Perl VEP's `hgvs_transcript`, ported whole: the 3' shift within
  1,000 bases of flank, the variant typing and duplication lookup of
  `hgvs_variant_notation`, allele clipping and `_get_cDNA_position`.
  Deletions that cross an exon boundary or the start or stop codon,
  reverse-strand `delins` ranges, insertions without a readable reference
  and variants outside the transcript span now match VEP.
- `HGVS_OFFSET` is emitted beside `HGVSc` and `HGVSp` when an insertion or
  deletion was shifted, signed by the transcript strand as VEP signs it.

### Added

- A GRCh37 chromosome 21 golden corpus generated with `--hgvs`
  (`tests/golden/115/GRCh37-hgvs/`), carrying the complete chromosome 21
  reference so the golden tests compare `HGVSc`, `HGVSp` and `HGVS_OFFSET`
  against Ensembl VEP's output in every format with no external file.
- Release measurements beside the paper's: the README's "Current release"
  section states the released version's ten-suite concordance and its wall
  times (20 independent machines per cell on both architectures) beneath the
  paper's tables, which stay as published, from the release's provenance
  record under `docs/concordance-provenance/`.

### Changed

- Releases: the workflow creates a draft that a maintainer publishes, with a
  body rendered from the repository's files (`scripts/release/render_release_notes.py`);
  the archives are `vep-X.Y.Z-<target>.tar.gz`, each unpacking to its own
  directory; a `git archive` source tarball ships with its `.sha256`; one
  `SHA256SUMS` covers every archive and every file in it carries a build
  provenance attestation (`gh attestation verify`); `MD5SUMS` is not published.
  `workflow_dispatch` runs the same build as a dry run without a release.
- Dependencies: noodles 0.109 (the last release declaring rust-version 1.88;
  the bgzf readers and writer and the tabix index loader moved to their
  `io` and `fs` modules), mysql 28 (vep-cache-builder), rand 0.10
  (vep-cache-converter), criterion 0.8 (benches), rustls 0.23.45, and the
  minor and patch releases of bigtools, bstr, clap, crossbeam-channel,
  flate2, indexmap, rayon, rustc-hash, serde, serde_json, smallvec, tempfile,
  thiserror and tracing-subscriber. The declared rust-version stays 1.88.
- Workflows: every action is pinned to a commit sha with its tag beside it,
  each rust-toolchain step names its Rust version through the `toolchain:`
  input, and the CI workflow's token is read-only. Dependabot ignores
  noodles 0.110 and later (rust-version 1.89) and the pinned toolchain refs.
- Golden corpus provenance records the digest of the released 0.1.0 binary
  that classified the manifests.
- The measurement harness (`scripts/concordance/run_clone_measurement.sh`)
  flushes and deletes the discarded warmup's output before every timed run,
  so the timed run starts with the memory the warmup had; the warmup's
  pages otherwise had to be reclaimed while the timed run wrote its own,
  and that reclaim landed in the timed wall time while user and system time
  stayed unchanged.

### Security

- Cleared from the dependency tree: RUSTSEC-2026-0235 (rkyv),
  RUSTSEC-2026-0173 (proc-macro-error2), RUSTSEC-2026-0002 (lru 0.12) and
  RUSTSEC-2026-0097 (rand 0.8), all through the mysql and rand updates, and
  RUSTSEC-2026-0285 (rustls) through the rustls update. None was reachable
  from the `vep` binary; all sat under vep-cache-builder or
  vep-cache-converter.

## [0.1.0] - 2026-09-22

First public release.

### Added

- `vep`, a from-scratch Rust implementation of the Ensembl Variant Effect
  Predictor algorithm. It reads a JSON transcript cache converted from the
  Ensembl VEP cache, accepts VCF, Ensembl default, HGVS and region input, and
  annotates SNVs, indels, MNPs and structural variants (deletions,
  duplications, insertions, inversions, copy-number variants and breakends)
  with Sequence Ontology consequence terms, HGVS notation, protein predictions
  and co-located variation.
- Output in the four Perl VEP formats (VEP default, VCF, JSON, Tab), each
  field for field Ensembl VEP 115.2's for the same input and flags, for every
  row outside the keys the golden manifests list as divergent (header lines carrying a time,
  a path or Perl API component versions excepted), plus Parquet, a vep-rs
  extension: five key columns (`chrom`, `pos`, `end`, `ref`, `alt`) then
  exactly the tab columns, `-` stored as NULL, Parquet v2, ZSTD 9, sorted row
  groups, a dictionary and a Bloom filter on every column, and the run's
  identity in the footer metadata.
- Golden corpora under `tests/golden/<release>/<assembly>/`: per Ensembl
  release and assembly, exemplar records for every consequence-set combination
  observed on the benchmark datasets, Ensembl VEP's output for them in the
  default, tab, VCF and JSON formats, a pruned JSON cache, and a manifest of
  the consequence keys documented to differ. `crates/vep-cli/tests/golden.rs`
  and `format_parity.rs` compare every column of every format against them on
  each `cargo test`, and the Parquet output round-trips through
  `scripts/adapters/parquet_to_vep_tab.py` to the tab output. Generator and
  pruner under `scripts/golden/`.
- `scripts/concordance/compare_vep_fields.py`: per-column and per-Extra-key
  agreement between two default-format outputs over their shared consequence
  keys, written as `fields.json` beside every F1 report the measurement
  harness produces.
- `--max_sv_size` (default 10,000,000): as in VEP, a wider structural variant
  keeps its VCF line without consequences and is absent from the JSON output;
  the default and tab formats list every transcript it overlaps.
- Parallel annotation with `--fork N`; output is deterministic and identical at
  every `--fork` value.
- Built-in plugins for CADD, REVEL, gnomADc, AlphaMissense, dbscSNV, LoFTEE,
  LoFtool and pLI on both GRCh37 and GRCh38, GWAS and SpliceAI on GRCh38, and
  dbNSFP, built in but unsupported;
  a binary annotation-store format for plugin data; and a dynamic-library
  plugin ABI (`vep-plugin`) for plugins built outside this repository.
- `vep-cache-builder`, which builds the JSON transcript cache directly from
  Ensembl source data (GFF3, protein FASTA, variation VCF, and with `--gtf`
  the release's GTF for the `cds_start_NF` / `cds_end_NF` transcript
  attributes and the assembly and genebuild versions) without a Perl
  installation, and `vep-cache-converter`, which validates an existing JSON
  cache for runtime use and converts tabix-indexed plugin data to the binary
  annotation-store format. Every JSON cache carries its `source_versions` in
  `info.json`.
- A concordance harness that compares vep-rs output against Perl VEP
  (release 115) tuple by tuple and classifies every divergence, and
  structural-variant validation VCFs for both assemblies under
  `tests/sv_validation/`.
- Documentation under `docs/` covering the CLI, output formats, plugins and
  cache setup, the runbook for reproducing the published concordance under
  `scripts/`, and the data and figure sources for the accompanying manuscript
  under `manuscript/`.
- `scripts/check_advisory_reachability.sh`, which fails when a known
  third-party advisory is reachable from the `vep` binary rather than merely
  present somewhere in the workspace, with an expiring allowlist for accepted
  risks.
- Continuous integration on GitHub Actions: tests, Clippy, formatting, a
  build on the declared minimum supported Rust version, a module-documentation
  check, a Developer Certificate of Origin check on pull
  requests, and a release workflow that publishes Linux (x86_64, aarch64) and
  macOS (aarch64) tarballs with sha256 and md5 checksum files.

[Unreleased]: https://github.com/natera-open-source/vep-rs/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/natera-open-source/vep-rs/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/natera-open-source/vep-rs/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/natera-open-source/vep-rs/compare/v0.3.2...v0.4.0
[0.3.2]: https://github.com/natera-open-source/vep-rs/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/natera-open-source/vep-rs/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/natera-open-source/vep-rs/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/natera-open-source/vep-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/natera-open-source/vep-rs/releases/tag/v0.1.0
