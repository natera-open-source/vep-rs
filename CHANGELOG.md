# Changelog

All notable changes to this project are recorded in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- The README states fastVEP and Ensembl VEP ClinVar times with ratios, drops vep-rs v0.4.0's re-scored rows, places Installation and Quick start before the figures and indexes the documentation site. (#77)

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

[Unreleased]: https://github.com/natera-open-source/vep-rs/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/natera-open-source/vep-rs/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/natera-open-source/vep-rs/compare/v0.3.2...v0.4.0
[0.3.2]: https://github.com/natera-open-source/vep-rs/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/natera-open-source/vep-rs/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/natera-open-source/vep-rs/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/natera-open-source/vep-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/natera-open-source/vep-rs/releases/tag/v0.1.0
