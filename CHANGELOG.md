# Changelog

All notable changes to this project are recorded in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Golden-corpus manifests carry `covers` (the parity-matrix rows a corpus exercises) and `flag_sets` (further runs under other flags, inputs or stdin), and the golden tests run every set. (#N)
- `check_parity_matrix.py` accepts a `scripts/` pytest function as a process row's `test1`, requires its `perl_reference` to be `-`, and refuses a `covers` name that is not a row. (#N)
- `scripts/concordance/compare_vep_formats.py` compares every column of every row against Ensembl VEP in the default, tab, VCF, JSON and Parquet formats, with `--self-parity` and `--warnings` modes, replacing `compare_output_formats.py`. (#N)
- Tabix-indexed data opens through a `.csi` as well as a `.tbi`, in htslib's order, and a record whose span reaches a queried position is returned as `tabix` returns it. (#N)
- `tools/vep-fields-compare` builds the per-column field comparator in Rust, byte-identical to `compare_vep_fields.py`, and `run_clone_measurement.sh` runs it when `VEP_FIELDS_COMPARATOR` names it. (#N)
- `generate_reference_output.sh` takes `--flags`, `--format default|tab|vcf|json|all`, `--stats` and `--perl-image-digest`, and `--corpus` writes a golden corpus's four expected outputs and `provenance.json`. (#N)
- `--grch37`, `-o` repeated (last wins with Ensembl VEP's warning), `STDOUT` in any case, and `docs/cli-reference.md` generated from the command line, the option tables and the error table. (#N)
- `storable_to_json.pl` converts every chromosome directory of a Storable cache and exports the UniProt isoforms, protein-feature labels, SIFT and PolyPhen matrices, `chr_synonyms.txt` and `contig_order`. (#N)
- `storable_to_json.pl` refuses a Sereal cache, a shard-less cache, an unlisted chromosome and any shard that fails to read or write, with one `ERROR:` line and exit 1. (#N)
- `cache_parity_audit.pl` classifies array elements and the new exports, so a conversion audits with no unexplained path. (#N)
- `scripts/concordance/feature_oracle.sh` runs a case table's comparisons, Ensembl VEP against candidate `vep` binaries in every output format, and writes per-case reports and a `status.tsv`. (#N)
- `--transcript_filter` honours Ensembl VEP's filter grammar: every operator, parenthesised groups, `not`, `in` lists and files, repeated flags joined by `and`, and `is_child` over a Sequence Ontology table. (#N)
- Converted caches load MANE Select, MANE Plus Clinical, APPRIS and TSL from transcript attributes and the GENCODE sets by attribute presence, so the attribute columns and `--pick` tiers fill. (#N)
- Every chromosome directory of a cache loads, patches, alternate loci, scaffolds and `LRG_` regions included; a cache with no transcript shard is refused with an error naming its serialiser. (#N)
- `--cache` takes Ensembl VEP's optional directory value, the cache root when `--dir_cache` is absent. (#N)
- `--full_cache_dir` (an alias for `--json_cache`), `--show_cache_info` (sorted `key<TAB>value` source versions, exit 0), `--no_fasta`, and the cache directory's first indexed `.fa` as the reference. (#N)
- `--check_existing` reads Ensembl VEP's own variation files, `variations/<chr>/all_vars.gz` through its index on first use, with columns from the cache's `variation_cols`; the JSON variation shards are gone. (#N)
- `--old_maf` is accepted and rewrites each co-located variant's AFR, AMR, ASN and EUR values as Ensembl VEP does, so those columns print empty on a release 116 cache. (#N)
- `storable_to_json.pl` and `vep-cache-converter` copy a cache's variation files (`all_vars.gz` with its index, or the `_var.gz` regions) into `variations/<chr>/`. (#N)
- `prune_json_cache.py --variations` cuts a corpus's variation slice as bgzip with a `.csi`, `--exclude-variation` drops records by name, and `info.json` keeps `variation_cols` only beside a slice. (#N)
- An unreadable variation file ends the run with `ERROR: [cache]` before its batch is written; the library's per-row API returns the failure on that batch's every row. (#N)
- `vep-cache-builder --core-dumps <dir>` reads the Ensembl core and variation table dumps, so a native cache carries the symbol source, HGNC id, xrefs, protein features, phenotype counts and attributes. (#N)
- `--warning_file` names the warnings file (default `<output_file>_warnings.txt`, `STDERR` for standard error alone); every warning is written there once and echoed on standard error unless `--quiet`. (#N)
- `--skipped_variants_file` lists every record the run skipped in Ensembl VEP's layout; an existing file is refused before the run without `--force_overwrite`. (#N)
- `--safe` makes an unknown `--plugin` name or a missing plugins directory an error instead of a warning and skips the incompatible-pair check. (#N)
- Native caches carry Ensembl's exon phases, the N-padded translateable sequence, both UTRs, genomic-order sorted exons, transcript-order introns and codon table 2 on MT, as Ensembl VEP's does. (#N)
- `--variation-dumps` reads the ClinVar allele and somatic classification columns from the variation database dumps, which no Ensembl VCF carries. (#N)
- Every output format names the vep-rs version: `## vep-rs <version>` in the default and tab headers, `##vep-rs=<version>` in the VCF header, `vep_rs_version` in every JSON object. (#N)
- `--overlaps` adds the `OverlapBP` and `OverlapPC` columns to the tab and VCF outputs, as in Ensembl VEP, on transcript rows of non-breakend structural variants. (#N)
- `run_clone_measurement.sh --vep-extra-args "<flags>"` appends arguments to every vep-rs invocation through a wrapper and records the binary, its sha256 and every argument vector in `provenance.json`. (#N)
- `build_golden_corpus.py select --keep-samples N` keeps `FORMAT` and the first N sample columns and samples records by the genotype classes they show, which `check` asserts are all covered. (#N)
- `--fields` replaces the tab columns and the VCF CSQ subfields with the list given, a name nothing produces printing `-` under a `?` description, as Ensembl VEP does. (#N)
- `--keep_csq` keeps an input VCF's existing annotation, header lines and INFO entry beside the new one; `--compress_output gzip|bgzip` writes the output through an in-process gzip or BGZF encoder. (#N)
- `--gencode_basic` filters the transcript set as the cache loads, and `--gencode_primary` and `--flag_gencode_primary` arrive with Ensembl VEP's incompatibility rule, GRCh38 warning and `GENCODE_PRIMARY` column. (#N)
- Golden-corpus flag sets carry `unordered` (`build_golden_corpus.py flag-set --unordered`) for sets whose rows Ensembl VEP emits in hash order; the `GRCh38-pick-summary` corpus runs twelve pick and summary sets. (#N)
- `--no_check_alleles`, `--exclude_null_alleles` and `--failed <int>` are accepted and match known variants on their span alone, drop null-allele records and admit failed ones as Ensembl VEP does. (#N)
- The `GRCh38-colocated-alt` golden corpus: span-only matches, null alleles, failed records and allele-specific ClinVar entries with Ensembl VEP's default-flag output. (#N)
- Read the `--config` file, `<dir>/vep.ini` and `VEP_*` variables in Ensembl VEP's priority, the command line winning. (#N)
- Long options are read in any case, as a unique prefix and after one dash, as Getopt::Long does; an Ensembl VEP option this version lacks is refused by name. (#N)
- `--shift_3prime 0|1`, `--shift_hgvs 0|1`, `--shift_genomic 0|1` and `--shift_length` move an insertion or deletion 3' as Ensembl VEP does, with `SHIFT_LENGTH`, shifted `Location` and `Allele`, and an unshifted HGVS mode. (#N)
- `--synonyms <file>` and the cache's `chr_synonyms.txt` resolve an input chromosome named by accession or `chr` prefix to the cache's name, as Ensembl VEP does. (#N)
- `vep-cache-builder` writes `info.json` with `contig_order`, `sift`, `polyphen` and `variation_cols` set from what the build wrote, and `chr_synonyms.txt` from the core dumps. (#N)
- `--fasta` reads a bgzipped FASTA through its `.fai` and `.gzi`, checking the index against the file at open, and reads a directory as every indexed `.fa` and `.fa.gz`. (#N)
- `--fasta_dir` selects the FASTA named for the species and assembly, `toplevel` and bgzipped first, as Ensembl VEP does, and warns when none qualifies. (#N)
- `--hgvsp_use_prediction` writes HGVSp in the predicted form, `p.(Arg97Ter)`, and `--remove_hgvsp_version` leaves the translation version off its reference name, as Ensembl VEP does. (#N)
- `--clin_sig_allele`, `--clinvar_somatic_classification` and `--var_synonyms` are accepted, the last switching `--check_existing` on with Ensembl VEP's notice; `VAR_SYNONYMS` and `CLINVAR_SOMATIC_CLASSIFICATION` print as it does. (#N)
- The `GRCh38-colocated`, `GRCh37-colocated` and `GRCh38-colocated-alt` corpora gain flag sets running `--check_existing` with the co-location flags against Ensembl VEP 116.2's output. (#N)
- `--terms` (`-t`), `--no_escape`, `--uploaded_allele` and `--ambiguity` do what Ensembl VEP's do: display or NCBI vocabularies, unescaped HGVSp and symbols, the uploaded allele string, the IUPAC code. (#N)
- `--input_data` (alias `--id`) replaces the input file; BCF and CAid input are refused by name; the `## VEP command-line:` header names each option by its resolved name. (#N)
- `--transcript_version` and `--gene_version` append the stable-id version to the Feature and Gene columns, as Ensembl VEP does. (#N)
- `--spdi` writes each allele's SPDI notation as Ensembl VEP's `spdi_genomic` does, from the record's alleles as the parser leaves them, with no reference FASTA required. (#N)
- `--ga4gh_vrs` writes the GA4GH VRS allele Ensembl VEP writes with a database connection, the `NC_` accession read from the chromosome synonyms table; refused with `--vcf`. (#N)
- `--hgvsg_use_accession` names the HGVSg chromosome by the first `NC_` synonym of the input name, as Ensembl VEP does, and HGVSg follows `--shift_hgvs 0`. (#N)
- `vep-cache-builder --include-variations` writes Ensembl VEP's `variations/<chr>/all_vars.gz` and `.csi` from the per-chromosome, somatic, phenotype, clinical and 1000 Genomes VCFs, replacing the JSON variation shards. (#N)
- `--nearest transcript|gene|symbol` adds `NEAREST`, the protein-coding transcription start sites nearest each record, found as Ensembl VEP's transcript tree finds them; every transcript at a shared site is listed. (#N)
- `--af_gnomad` is accepted as the older spelling of `--af_gnomade`, refused beside it, and fails on a cache without gnomAD exome columns as Ensembl VEP fails. (#N)
- The `GRCh38-colocated`, `GRCh37-colocated` and `GRCh38-colocated-alt` corpora run their flag sets with the frequency flags in every format, and `GRCh38-oldmaf` gains an `--old_maf` set. (#N)
- The plugin-facing variant JSON carries `input_chr`, the input's spelling, on a record renamed from `M` to `MT`; absent otherwise. (#N)
- The `GRCh38-input-compressed` corpus runs one input as BGZF (`.gz` and `.bgz`), plain gzip, standard input, CRLF and a Latin-1 header against one expected output. (#N)
- `--check_svs <file>` reports the known structural variants overlapping each variant (`SV`, JSON `colocated_structural_variants`) from Ensembl's structural-variant VCF offline, a breakend allele's mate position included. (#N)
- The `GRCh38-checksvs` golden corpus: a slice of Ensembl's structural-variant dump with Ensembl VEP 116.2's `--check_svs` output from its database in four formats. (#N)
- `--domains` lists the protein features the variant's residues overlap on every coding row, as `label:identifier`, from converted and native caches alike. (#N)
- `--humdiv` reads the PolyPhen column from the HumDiv matrix instead of HumVar; a cache without that matrix prints an empty column, as in Ensembl VEP. (#N)
- `prune_json_cache.py --drop-except <file>` keeps every `--drop` key on the transcripts the file names, so a corpus carries matrices and protein features for its own transcripts only. (#N)
- `--check_frequency`, `--freq_pop`, `--freq_freq`, `--freq_gt_lt`, `--freq_filter` and `--filter_common` filter the input by co-located variant frequency as Ensembl VEP does, with `FREQS`. (#N)
- The `GRCh38-freqfilter` and `GRCh38-freqfilter-include` corpora compare both engines under the exclude, include and `--filter_common` forms in every output format. (#N)
- `--delimiter` splits an Ensembl-format or VCF input's columns (a tab, a space, a run of spaces, `\s+` or literal text); under `--format guess` the detected one applies. (#N)
- The `GRCh38-input-ensembl` corpus compares Ensembl-format input, tab- and space-delimited, and a space-delimited VCF against Ensembl VEP 116.2's output. (#N)
- `--format spdi` reads SPDI notation (`sequence:position:deletion:insertion`) as Ensembl VEP does, the reference allele checked against or read from `--fasta`, which the format needs. (#N)
- An SPDI record is named by its notation and prints under the chromosome name the cache uses, resolved through the synonyms. (#N)
- A notation Ensembl VEP's adaptor refuses is warned about with its text and skipped, the run going on; a span over `--max_sv_size` is skipped with Ensembl VEP's warning. (#N)
- The `GRCh38-input-spdi` corpus compares SPDI input, under `--format spdi` and under format detection, against Ensembl VEP 116.2's database-mode output. (#N)
- `--decompression_threads` takes `auto` (opt-in, the default unchanged), reading an indexed BGZF input's record shape; hidden `--resolved_config` prints the resolved knobs and that shape before any cache opens. (#N)
- The `GRCh38-input-hgvs` corpus compares HGVS g. and m. input with Ensembl VEP 116.2's database-mode output in the default, tab and JSON formats. (#N)
- The `GRCh38-everything` golden corpus compares `--everything` with Ensembl VEP 116.2 on 204 chromosome 21 records; the reference's regulatory rows are masked, since vep-rs annotates no regulatory feature. (#N)
- Golden-corpus manifests carry `masks`, and `compare_vep_formats.py --masks regulatory_feature_rows` and `build_golden_corpus.py classify` apply the same class: the reference's RegulatoryFeature and MotifFeature rows are set aside. (#N)
- The GRCh38-native-builder corpus runs its own native variation files under `--check_existing --exclude_null_alleles --af --af_1kg --clinvar_somatic_classification` against Ensembl VEP 116.2's output. (#N)
- `prune_json_cache.py --variation-info` writes a variation slice's own `variation_cols` and `var_type` and refuses a declaration whose width differs from the slice's records. (#N)
- `--check_ref` skips a record whose reference allele differs from the FASTA, `--lookup_ref` replaces it from the FASTA, and `CHECK_REF` prints `failed` on a kept record, as Ensembl VEP does. (#N)
- `--chr` restricts the run to the chromosomes listed, matched against the input's own spelling, with Ensembl VEP's `1-5,10,MT` range grammar. (#N)
- Every input record passes Ensembl VEP's `validate_vf` checks before annotation: a bad coordinate, allele string, insertion or reference-allele length skips it with Ensembl VEP's warning text. (#N)
- `--allow_non_variant` keeps a VCF record whose ALT is `.` in every output format: its line in the VCF output, one object in the JSON output, no row elsewhere. (#N)
- Region input needs `--fasta`, is refused with `--check_ref`, and reads its reference allele from the FASTA, as Ensembl VEP's region parser does. (#N)
- The `GRCh38-vcf-check-ref` golden corpus: reference checks, skipped records, `--dont_skip`, `--chr`, `--lookup_ref` and `--allow_non_variant` against Ensembl VEP 116.2's output. (#N)
- Container image for linux/arm64 and linux/amd64 from generic-CPU archives; the measured binaries ship as the `X.Y.Z-neoverse-v1` and `X.Y.Z-x86-64-v3` tags. (#N)
- `--format id` resolves a list of variant identifiers from the cache's variation files in one pass, one record per genomic mapping, as Ensembl VEP's database lookup does. (#N)
- An identifier the cache lacks is warned about as Ensembl VEP warns (`No variant found with ID`) and the list goes on; a mapping failing coordinate checks is skipped. (#N)
- The `GRCh38-input-id` corpus compares identifier input, under `--format id` and under detection, against Ensembl VEP 116.2 connected to the Ensembl database. (#N)

### Changed

- Ensembl VEP's option rules apply as data: `--everything` is its exact set (no `ALLELE_NUM` or `HGVSg`), option sets override explicit values, incompatible pairs such as `--vcf --json` are refused. (#N)
- Every refusal prints one `ERROR: [<component>] <message>` line and exits 1; argument syntax errors stay clap's and exit 2; an empty command line prints the help. (#N)
- `--database`, `--host`, `--port`, `--user`, `--password`, `--registry`, `--db_version`, `--genomes`, `--is_multispecies`, `--lrg`, `--convert`, `--gp`, `--ambiguous_hgvs` and `--af_exac` are refused by name. (#N)
- `--dir_cache`, `--species`, `--cache_version`, `--refseq` and `--merged` select the cache in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout; a run with no cache exits 1 before reading input. (#N)
- `--assembly` must equal the cache directory's and `info.json`'s assembly (`hg19` and `hg38` read as `GRCh37` and `GRCh38`); the cache's assembly feeds headers and plugins. (#N)
- A converted cache carries the Ensembl VEP cache's variation files (about 23.5 GB more on GRCh38, 21.8 on GRCh37); the `minor_allele`, `minor_allele_freq`, `AA` and `EA` fields are gone. (#N)
- Under `--check_existing` a cache without variation data warns once that `Existing_variation` and the frequency columns stay empty; Ensembl VEP annotates silently. (#N)
- A cache converted by an earlier release needs the Ensembl VEP cache's `all_vars.gz` and `.csi` files copied into `variations/<chr>/`, or a new conversion. (#N)
- `vep-cache-builder` reads the GFF3's `ccdsid`, `transcript_support_level`, `tag` codes, gene `version` and CDS `protein_id`, keeps every non-coding gene and drops `artifact` and readthrough transcripts. (#N)
- A native cache carries MANE, TSL and APPRIS as transcript attributes alone, as the Ensembl cache does, and `is_canonical` only on the canonical transcript. (#N)
- A record Ensembl VEP skips (unknown chromosome, invalid alleles, unsupported or oversize SV, incomplete deletion) gets its `line N skipped` warning, N the record ordinal. (#N)
- `--dir_plugins` defaults to `<--dir>/Plugins`; a directory that does not exist is a warnings-file line when `--plugin` names a plugin. (#N)
- Native variation records take Ensembl VEP's form: the export's alleles kept at their position, an insertion's start past its end, `phenotype_or_disease` from the evidence flag. (#N)
- `--distance up,0` reads the downstream distance as the upstream one, as Ensembl VEP does; a value that is not one or two integers exits with an error. (#N)
- `--hgvsg` writes each row's HGVS genomic notation as Ensembl VEP's `hgvs_genomic` does; it needs a reference FASTA and exits with Ensembl VEP's error without one. (#N)
- `--output_format` is matched in any case, `ensembl` naming the default format; plugin header descriptions reach the default, tab and VCF headers; `--fields` or `--compress_output` with Parquet is refused. (#N)
- The APPRIS and TSL `--pick` tiers grade attribute values by Ensembl VEP's rules, so `principalM` and a zero carry no grade and `alternative` adds 10 only in lower case. (#N)
- The pick family (`--pick`, `--pick_allele`, `--per_gene`, `--pick_allele_gene` and the `--flag_pick*` flags) acts per record as Ensembl VEP does: one row per record, allele, gene or both, intergenic rows too. (#N)
- `--most_severe` and `--summary` print one row per record, `-` in the allele, feature and position columns, variant-level keys alone in Extra, and `intergenic_consequences` in JSON. (#N)
- `--pick_order` is honoured with eleven criteria in the given order and Ensembl VEP's default list; `length` reads the translateable or spliced length, `ensembl` and `refseq` the transcript's source cache. (#N)
- `--check_existing` matches known variants by allele after trimming from either end, strand-aware and span-exact, in Ensembl VEP's order; structural variants are no longer co-located. (#N)
- A co-located variant matched by allele prints only on the rows of the alleles it matched, so a multi-allelic record's rows differ as Ensembl VEP's do. (#N)
- `--hgvs` or `--hgvsg` without a FASTA stops with Ensembl VEP's message; `--everything` drops HGVS with its `INFO` line; `--shift_3prime 1` without one is refused. (#N)
- An Ensembl-format allele given on the reverse strand is analysed as its forward-strand form, so its terms, codons and HGVS notations match the forward spelling. (#N)
- `vep-cache-converter` checks `info.json` (`species`, `assembly`, its list keys) and copies `chr_synonyms.txt` beside it. (#N)
- A plain-gzip `--fasta` is refused with the `bgzip` and `samtools faidx` remedy; the cache directory's first indexed `.fa.gz` counts as its FASTA too. (#N)
- `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED` print as Ensembl VEP does: allele-specific ClinVar entries by default, one flag per variant with all-zero lists dropped, duplicates kept. (#N)
- JSON `colocated_variants` entries carry every set cache column as Ensembl VEP writes them, frequencies nested per allele; the record-level `existing_variation` array is gone. (#N)
- `--coding_only` keeps the rows whose allele affects a transcript's coding sequence and `--no_intergenic` drops the intergenic rows, as Ensembl VEP's output factory filters them, before the pick family. (#N)
- `AnnotationResources` carries the compiled `--transcript_filter` and `FilterConfig` the `--coding_only` and `--no_intergenic` row filters, so an embedding caller sees the rows the `vep` binary writes. (#N)
- `--total_length` appends the cDNA, CDS and peptide lengths to the position columns, and `VARIANT_CLASS` is the whole record's class, so a multi-allelic record prints one. (#N)
- The input format is detected from the first data line as Ensembl VEP detects it, and compression from the file's bytes (BGZF or plain gzip, whatever the name). (#N)
- `-i STDIN` needs `--format`; a missing input, an undetectable first line, a region input without FASTA and an HGVS `r.` line are refused with Ensembl VEP's texts. (#N)
- MANE, APPRIS, TSL, CCDS, RefSeq, UniProt, GENE_PHENO and miRNA print Ensembl VEP's values: accessions, `P1`/`A2`/`PM` grades, digits, lists, no `-` sentinel. (#N)
- Native variation records take the variation database's form: `TSA=insertion` records two bases right of the VCF anchor, COSMIC and HGMD-PUBLIC variants as `COSMIC_MUTATION`/`HGMD_MUTATION`, `clin_sig` from the `CLIN_<term>` flags. (#N)
- A native cache built by an earlier release holds JSON variation files nothing reads; build it again with `vep-cache-builder` 0.5.0 for co-located variants and the 1000 Genomes columns. (#N)
- `AF` prints with four decimals and the frequency columns report the record allele matched to the row's own allele, merged across co-located variants, as Ensembl VEP does. (#N)
- `MAX_AF` prints its maximum over every frequency group, combined populations excluded, and `MAX_AF_POPS` lists every tied population in Ensembl VEP's order. (#N)
- VCF input is read as bytes: a non-UTF-8 header line or sample name reaches the VCF output unchanged, and a data line with such bytes is annotated. (#N)
- Empty and whitespace-only VCF lines are skipped and warned about in the words of Ensembl VEP's `skip_empty_lines`; the record numbering does not count them. (#N)
- A record on chromosome `M` prints `MT` only where Ensembl VEP renames it: the cache knows `MT` and neither `M` nor `chrM`. (#N)
- `--decompression_threads` above the host's logical CPUs is lowered to that count, so a huge value reads the input instead of exhausting threads. (#N)
- Tabix-indexed VCF data reads a record's end from its INFO `END`, as `tabix -p vcf` indexes it, so a long record is found from every position it spans. (#N)
- `--sift` and `--polyphen` predict every single-residue substitution in the coding sequence, start codon changes included; a cache naming no version of the tool stops the run. (#N)
- `--sift s` and `--polyphen s` print the bare score and set `sift_score` alone in JSON, a whole score an integer; JSON `domains` follow Ensembl VEP's `(\w+):(\w+)` re-split. (#N)
- The consequence path reads biotype facts once per transcript, hoists the coding predicate once per allele, and keeps differing regions and the term list inline. (#N)
- The pipeline streams by default (`--pipeline_mode batched` restores whole-batch passes); hidden `--pipeline_depth`, `--render_task_records`, `--writer_strategy` and `--output_buffer_size` knobs write the same bytes at every setting. (#N)
- A reverse-strand record's alleles are read against each transcript's strand in every coding predicate; only its nucleotide alleles are reverse-complemented, a structural-variant word or symbolic allele staying as written. (#N)
- `Codons` complements an IUPAC letter of an allele on the strand opposite its transcript's (`gcY`, not `gcR`), as Ensembl VEP's `reverse_comp` does. (#N)
- `vep_core::codon::complement_base` complements the IUPAC letters as Ensembl's `reverse_comp` does, so a natively built cache's reverse-strand sequences agree with Ensembl's; `vep_effects::coding::reverse_complement_pub` is removed. (#N)
- `--freq_pop 1KG_ALL` and `--filter_common` read the cache's `AF` column, where Ensembl VEP 116.2 finds no 1000 Genomes data and filters nothing. (#N)
- Ensembl-format input reads its strand as Ensembl VEP does: `-` and `-1` are reverse, any other token or an absent column forward, so four-column lines are accepted. (#N)
- An Ensembl-format allele column without a `/` (`DUP`, `DEL`, `INS`, `INV`) is a structural variant over the line's span; an unknown word is warned about and annotated as written. (#N)
- Ensembl-format alleles are upper-cased and a bi-allelic indel written un-minimised (`AC/A`) is trimmed to Ensembl VEP's representation before annotation. (#N)
- Empty lines in an Ensembl-format input are skipped with the `skip_empty_lines` warning; a skipped Ensembl-format record's warning quotes its line. (#N)
- `--format spdi` without `--fasta` is refused with `Cannot use SPDI format in offline mode without a FASTA file`, where every SPDI input was refused before. (#N)
- HGVS genomic input (`g.`, `m.`) resolves against `--fasta` as Ensembl VEP's database does: `del`, `dup`, `inv` take genome bases, `delins` is trimmed, `ins` swaps positions. (#N)
- An HGVS notation whose reference allele the genome contradicts, a repeat, an unknown shape or a variant `validate_vf` refuses is skipped with Ensembl VEP's warning. (#N)
- HGVS input needs `--fasta`, is named by the notation in `Uploaded_variation` and the JSON `id`, and prints an accession or `chr` spelling under the FASTA's chromosome name. (#N)
- `--dont_skip` keeps a record that fails validation and annotates it, as Ensembl VEP does; a malformed line is skipped with its warning and the run completes with exit 0. (#N)
- A skipped record's warning quotes its input line in every input format, as Ensembl VEP's `line N skipped (...)` does. (#N)
- Lower-case alleles in an Ensembl-format record are upper-cased before annotation, every allele of a multi-allelic string included, under the name built from the line as written. (#N)

### Removed

- `vep_io::input::vcf::VcfParser` and its module; the `vep` binary's line parser is the one VCF reader. (#N)
- The transcript `description` key, which Ensembl VEP's cache dumper deletes: the loader ignores it and `--transcript_filter` reads it as undefined on every cache, 0.4.0 native caches included. (#N)
- `--transcript_index_impl` and the `coitrees` and sorted-array transcript indexes: UCSC binning is the one overlap index, and the flag is refused as an unknown argument. (#N)

### Fixed

- A tabix query answered from several chunks returns every record once, and a line cut between chunks is read whole instead of losing a record or yielding a fragment. (#N)
- A tabix index built for another data file is refused with both files named instead of stopping the process. (#N)
- The release page's caption labels a cross-generation class pair the newest available class per architecture, says when ClinVar shares the whole genome's class, and refuses a missing record field. (#N)
- A cache whose chromosome directories are symlinks loads; it was refused as holding no transcript shards. (#N)
- `--mane_select` alone lists no MANE Plus Clinical transcript, UNIPARC and UNIPROT_ISOFORM carry the cache's accessions, and a VCF SYMBOL writes whitespace as `%20`. (#N)
- `vep-cache-builder` writes a `TSA=indel` variation record with the export's alleles at its position, so a left-aligned input co-locates with a native cache as Ensembl VEP co-locates it. (#N)
- Linux release archives run on glibc 2.17 or later, the floor read from each binary's ELF version needs and printed on the release page. (#N)
- A reverse-strand Ensembl-format or region record is matched to CADD, dbNSFP, LoFTEE and every tabix-backed plugin file on its forward-strand alleles, as Ensembl VEP matches it. (#N)
- HGVSp and shifted codons of a moved insertion turn the allele as the record wrote it, unturned on a reverse-strand transcript when the shift exceeds its length. (#N)

### Security

- Cleared from the dependency tree: RUSTSEC-2026-0346, RUSTSEC-2026-0347, RUSTSEC-2026-0348 and RUSTSEC-2026-0349 (rustls 0.23.46) and the RUSTSEC-2026-0253 unsoundness notice (lru 0.18.5 via mysql 28.0.4), none reachable from the `vep` binary. (#N)

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

[Unreleased]: https://github.com/natera-open-source/vep-rs/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/natera-open-source/vep-rs/compare/v0.3.2...v0.4.0
[0.3.2]: https://github.com/natera-open-source/vep-rs/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/natera-open-source/vep-rs/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/natera-open-source/vep-rs/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/natera-open-source/vep-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/natera-open-source/vep-rs/releases/tag/v0.1.0
