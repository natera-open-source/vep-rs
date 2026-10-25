// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CLI argument parsing for VEP using clap 4 derive API.
//!
//! All flags and aliases match the Ensembl VEP CLI interface so that users
//! can switch between implementations with minimal changes. Ensembl VEP
//! uses underscores in flag names (e.g., `--input_file`, `--dir_cache`),
//! so every underscore-containing flag has an explicit `long = "..."`.
//!
//! A flag the error table of [`crate::config`] governs (refused, or accepted
//! without effect) carries no help text here: [`Args::command_with_error_table`]
//! applies the table's, so `--help` and the generated CLI reference read one
//! source. Every flag lives under a help heading, the section of the reference
//! it is documented in.
//!
//! The command line reaches this parser through [`crate::argv::prepass`], which
//! adds the flags of the `--config` file, `vep.ini` and the `VEP_*` variables
//! and spells single-dash and other-case long options as declared here; a
//! unique prefix of a long name is accepted here as well (`infer_long_args`),
//! as Getopt::Long's `auto_abbrev` accepts it.

use clap::{ArgAction, CommandFactory, Parser};

/// Ensembl VEP's `--max_sv_size` default: the span, end minus start, above which a
/// record is oversize, and the repeat-sequence length above which a tandem repeat
/// stays symbolic.
pub const DEFAULT_MAX_SV_SIZE: i64 = vep_core::sv_record::DEFAULT_MAX_SV_SIZE as i64;

/// The help's first line is the binary's name and version, as `--version` prints.
const HELP_TEMPLATE: &str =
    "{name} {version}\n{about-with-newline}\n{usage-heading} {usage}\n\n{all-args}{after-help}";

#[derive(Parser, Debug)]
#[command(
    name = "vep",
    version = env!("CARGO_PKG_VERSION"),
    about = "Ensembl Variant Effect Predictor (Rust implementation)",
    long_about = "Predict functional effects of genomic variants.\n\n\
                  Usage: vep -i input.vcf -o output.txt --json_cache <dir>",
    help_template = HELP_TEMPLATE,
    args_override_self = true,
    infer_long_args = true
)]
// The field doc comments are the command's help text, which clap prints as
// written: a `<value>` placeholder or `[PATH]` is literal there, not rustdoc markup.
#[allow(rustdoc::invalid_html_tags, rustdoc::broken_intra_doc_links)]
pub struct Args {
    /// Configuration file, one flag per line as `<name> <value>` without the dashes (`format vcf`, `offline 1`); a flag the command line also gives keeps the command line's value
    #[arg(help_heading = "Input", long, value_name = "FILE")]
    pub config: Option<String>,

    /// Input file path, or "STDIN" (in any case) for standard input, which needs --format.
    /// A gzip or BGZF file is decompressed whatever its name; BCF is refused.
    /// A VCF whose `##fileformat` line names a major version other than 4, or a value with no `<digits>.<digits>` version, is refused with Ensembl VEP's message before any output is written; the minor version is not checked, and a file without the line or with an empty value is read
    #[arg(
        help_heading = "Input",
        short = 'i',
        long = "input_file",
        default_value = "STDIN"
    )]
    pub input_file: String,

    /// Input records as a string, read in place of --input_file; one record per line, the format detected as for a file
    #[arg(help_heading = "Input", long = "input_data", visible_alias = "id")]
    pub input_data: Option<String>,

    /// Input format. "guess" detects it from the first data line, as Ensembl VEP does: a single field is a region, SPDI, CAid, HGVS or identifier, several fields a VCF or Ensembl record.
    /// An identifier list (one variant name per line) is resolved in one pass over the cache's variation files, every name held in memory and every file read once: a name maps to each variation feature of that name in the cache (rsIDs and the other names the cache carries; structural variant names, synonyms and subsnp identifiers are not among them), a name the cache lacks is warned about and skipped, and a cache without variation files refuses the format.
    /// CAid input is looked up in the Ensembl variation database, which vep-rs does not reach, and is refused; a region, SPDI or HGVS input needs --fasta, which fills a region record's reference allele and fills or checks an SPDI record's.
    /// An SPDI span past the end of its sequence, an insertion after a position the sequence lacks or a deletion running past its last base, is skipped with `Span <start>-<end> lies past the end of sequence <name> (<length> bases)`, and an HGVS genomic (g., m.) span past the last base with `The reference FASTA has no bases at <name>:<start>-<end>`, where Ensembl VEP reads `N` for the missing bases and annotates the record at a position the chromosome lacks.
    /// An SPDI deletion of length 0 (`<name>:<position>:0:`) is skipped with `Deleted sequence length (0) leaves nothing to delete: span <start>-<end> is empty`, where Ensembl VEP keeps a `/-` feature that writes no row.
    /// An SPDI or HGVS sequence name is matched whole, in any letter case, against the cache's chromosomes, the FASTA's records and their synonyms, never as a prefix, so a version-less accession (`NC_000021`) is skipped, in SPDI with `Sequence name <name> not valid` and in HGVS with `The reference FASTA has no sequence named <name>`, where Ensembl VEP's slice adaptor resolves the name as a prefix of a synonym (`NC_000021.9`) and annotates the record.
    /// HGVS input is resolved against the FASTA and the transcript cache as Ensembl VEP resolves it through its database: a genomic notation (g., m.) reads its bases from the genome, a transcript notation (c., n.) finds its transcript by stable id, else by gene name through Ensembl VEP's pick order with its warning, places the position through the transcript's exons and keeps the transcript's strand, the reference allele is checked, a deletion, duplication or inversion takes its bases from the genome, and a notation Ensembl VEP cannot resolve is skipped with its warning; a protein notation (p.) finds its transcript by translation id, by transcript id or, every transcript of the gene tried in turn, by gene name, and resolves a substitution or delins to the one shortest nucleotide path to the new residue (several equally short paths are skipped with Ensembl VEP's message, as --ambiguous_hgvs is not accepted) and a deletion to its codons, frameshifts refused as Ensembl VEP refuses them; a run with transcript or protein notations reads the cache's assembled chromosomes and scaffolds once to find them, its alternate loci and patches only for an id found on none; an r. notation is refused
    #[arg(
        help_heading = "Input",
        long,
        default_value = "guess",
        value_parser = ["ensembl", "vcf", "hgvs", "id", "spdi", "region", "guess"]
    )]
    pub format: String,

    /// Delimiter between the columns of an Ensembl-format or VCF input, a tab unless given: \t for a tab (a JSON `input` then joins its fields with a tab, where Ensembl VEP writes the two characters), a space to split on whitespace, a space and a plus (` +`) to split on runs of spaces, `\s+`, text with no regular-expression character to split on that text, or any other text as the regular expression Ensembl VEP hands to Perl's split, read byte by byte as Perl reads it (a value its syntax refuses, one that can match the empty string, or one that could match part of a UTF-8 character is an error); a VCF record read under another delimiter is re-joined with tabs, so a field holding a tab is cut at it in the JSON `input` and the VCF line, where Ensembl VEP keeps the field whole. With --format guess the delimiter is detected from the first data line instead, as Ensembl VEP detects it, and this flag is not read
    #[arg(help_heading = "Input", long, value_name = "TEXT")]
    pub delimiter: Option<String>,

    #[arg(help_heading = "Input", long)]
    pub gp: bool,

    #[arg(help_heading = "Input", long = "ambiguous_hgvs")]
    pub ambiguous_hgvs: bool,

    #[arg(help_heading = "Input", long)]
    pub phased: bool,

    /// Annotate a VCF record once per individual, `all` for every sample column or a comma-separated list
    /// of sample names from the #CHROM header line; a single-`#` line after that line is ignored, where
    /// Ensembl VEP takes the header's last single-`#` line as the column line and so names no sample. An
    /// unknown name exits with Ensembl VEP's error naming both lists, checked before the first record,
    /// where Ensembl VEP checks it at its first sequence record and so annotates a VCF of structural
    /// records alone. Each row of the copy
    /// carries IND and ZYG (HET when the genotype holds more than one
    /// distinct allele, else HOM, with REF appended for a reference-homozygous individual); a copy carries
    /// the alleles its genotype names, so `./1` and a haploid `1` annotate the one allele as HOM, `1/2` both
    /// alleles as HET, `*` counts as no alternate allele, and `./.` writes nothing for that individual.
    /// A reference-homozygous individual's copy writes no row unless --process_ref_homs: the VCF output
    /// carries its line without consequences and the JSON output its object without them. Every output
    /// repeats the record once per individual, one VCF line and one JSON object each. A VCF without a
    /// #CHROM line, or whose header names no sample, writes nothing under `all`, where a VCF without that
    /// line stops Ensembl VEP with a Perl error. Implies
    /// --allow_non_variant; cannot be combined with --minimal or --individual_zyg
    #[arg(help_heading = "Input", long, value_name = "all|LIST")]
    pub individual: Option<String>,

    /// Annotate a VCF record once, with every named sample's zygosity in one ZYG list, `all` for every
    /// sample column or a comma-separated list of sample names from the #CHROM header line; an unknown
    /// name exits with Ensembl VEP's error naming both lists. Each row carries ZYG as `sample:HET`,
    /// `sample:HOM`, `sample:HETREF` or `sample:HOMREF` entries (HET when the genotype holds more than
    /// one distinct allele, REF appended when the sample carries no alternate allele, `*` counting as
    /// none), one per sample whose genotype is called, in the order the run names them (the header's
    /// column order under `all`) where Ensembl VEP prints its hash order; a sample whose genotype is
    /// `./.` has no entry, and a record no named sample has a genotype for carries no ZYG. When the
    /// run names one sample (a list of one name, or `all` over a single-sample VCF) and that sample
    /// carries no alternate allele, the record annotates its reference allele alone, as
    /// --process_ref_homs annotates it, with the allele string REF/REF. A structural variant carries
    /// no ZYG. Implies --allow_non_variant; cannot be combined with --individual
    #[arg(
        help_heading = "Input",
        long = "individual_zyg",
        value_name = "all|LIST"
    )]
    pub individual_zyg: Option<String>,

    /// Add the reference allele's row to every feature of every record, numbered 0 under --allele_number
    /// and listed before the alternate alleles' rows, with the consequence terms an allele equal to the
    /// reference takes (synonymous_variant in a coding sequence), no HGVSc or HGVSp, and the alternate
    /// allele's HGVSg on a bi-allelic record; a reference-homozygous individual under --individual
    /// annotates that row instead of writing none. A structural variant gains no reference row
    #[arg(help_heading = "Input", long = "process_ref_homs")]
    pub process_ref_homs: bool,

    #[arg(help_heading = "Input", long = "no_check_variants_order")]
    pub no_check_variants_order: bool,

    /// Annotate the records of these chromosomes alone, a comma-separated list matched against each
    /// record's chromosome as the input spells it (`21` does not admit `chr21`); a pair of integers
    /// joined by `-` names every chromosome between them (`1-5,10,MT`), which Ensembl VEP matches
    /// literally. A record outside the list is dropped without a warning, or kept by --dont_skip
    #[arg(help_heading = "Input", long, value_name = "LIST")]
    pub chr: Option<String>,

    /// Check each record's reference allele against the --fasta sequence over its span, on its strand;
    /// a mismatch skips the record with Ensembl VEP's warning, or keeps it under --dont_skip with
    /// CHECK_REF `failed`. Needs --fasta; cannot be combined with --lookup_ref
    #[arg(help_heading = "Input", long = "check_ref")]
    pub check_ref: bool,

    /// Replace each record's reference allele with the --fasta sequence over its span, on its strand,
    /// before annotation; the record's name keeps the allele as written. Needs --fasta; cannot be
    /// combined with --check_ref
    #[arg(help_heading = "Input", long = "lookup_ref")]
    pub lookup_ref: bool,

    /// Output file path ("STDOUT" in any case is standard output, as is "-" unless
    /// --compress_output writes a file, or Parquet a directory, of that name); the last of
    /// several wins
    #[arg(
        help_heading = "Output",
        short = 'o',
        long = "output_file",
        action = ArgAction::Append,
        default_value = "variant_effect_output.txt"
    )]
    pub output_file: Vec<String>,

    /// Output format, matched in any case; `ensembl` is another name for the default `vep` format. `parquet` is vep-rs's own: it writes a TSV intermediate and post-processes it
    /// to Parquet through a `duckdb` CLI subprocess; the `duckdb` binary must be on PATH
    #[arg(
        help_heading = "Output",
        long = "output_format",
        default_value = "vep",
        ignore_case = true,
        value_parser = ["vep", "ensembl", "vcf", "json", "tab", "parquet"]
    )]
    pub output_format: String,

    /// Parquet row shape (only effective when --output_format=parquet).
    ///   nested (default): one row per variant, CSQ subfields are LIST<T>.
    ///   flat: one row per (variant × transcript-consequence) tuple, CSQ
    ///   subfields are scalar.
    #[arg(
        help_heading = "Output",
        long = "parquet_shape",
        default_value = "nested",
        value_parser = ["nested", "flat"]
    )]
    pub parquet_shape: String,

    /// Keep the Parquet intermediate TSV after DuckDB finalize succeeds.
    /// Debug aid: normally the intermediate is deleted once the Parquet
    /// directory is written.
    #[arg(
        help_heading = "Output",
        long = "keep_intermediate",
        default_value_t = false,
        hide = true
    )]
    pub keep_intermediate: bool,

    /// Rows per Parquet row group (only effective when --output_format=parquet).
    /// Rounded up to a multiple of 2,048, DuckDB's vector size, which is also
    /// the dictionary size limit. Test aid for producing several row groups
    /// from a small input.
    #[arg(
        help_heading = "Output",
        long = "parquet_row_group_size",
        default_value_t = 122_880,
        hide = true
    )]
    pub parquet_row_group_size: u64,

    /// Shorthand for --output_format vcf; also switches on --symbol, --biotype and --numbers
    #[arg(help_heading = "Output", long)]
    pub vcf: bool,

    /// Shorthand for --output_format json
    #[arg(help_heading = "Output", long)]
    pub json: bool,

    /// Shorthand for --output_format tab
    #[arg(help_heading = "Output", long)]
    pub tab: bool,

    /// Compress the output, STDOUT included: `gzip` writes a single-member gzip stream, `bgzip` a BGZF stream (a gzip file that bgzip, tabix and samtools read); the output file keeps the name given. Not accepted with --output_format parquet, which carries its own compression
    #[arg(help_heading = "Output", long = "compress_output", value_parser = ["gzip", "bgzip"])]
    pub compress_output: Option<String>,

    /// Overwrite output file if it already exists
    #[arg(
        help_heading = "Output",
        long = "force_overwrite",
        visible_alias = "force"
    )]
    pub force_overwrite: bool,

    /// Suppress header lines in output
    #[arg(help_heading = "Output", long = "no_headers")]
    pub no_headers: bool,

    /// Print the cache directory in the output headers as [PATH]/<leaf>
    #[arg(help_heading = "Output", long = "mask_header_cache_path")]
    pub mask_header_cache_path: bool,

    /// File the warnings of the run are written to, created at the first warning;
    /// default <output_file>_warnings.txt beside the output (STDOUT_warnings.txt in the
    /// working directory for -o STDOUT, vep_warnings.txt for an output named 0 or
    /// empty). STDERR, in any case, keeps them on standard
    /// error, each written twice unless --quiet, as Ensembl VEP writes it. A warning
    /// is also echoed on standard error unless --quiet
    #[arg(help_heading = "Output", long = "warning_file", value_name = "FILE")]
    pub warning_file: Option<String>,

    /// File listing every input record the run skipped: a header line naming the
    /// input, then one line per record with its line number, the record's first 40
    /// characters and the reason. An existing file is refused before the run without
    /// --force_overwrite, where Ensembl VEP checks when the run ends; both engines write
    /// the file as the run ends
    #[arg(
        help_heading = "Output",
        long = "skipped_variants_file",
        value_name = "FILE"
    )]
    pub skipped_variants_file: Option<String>,

    /// The columns of the --tab output and the subfields of the VCF CSQ entry, comma-separated, in this order and replacing the whole default list; any field name a flag or a plugin produces is valid, and a name nothing produces prints `-` (tab) or an empty subfield (VCF) under Ensembl VEP's description of the name, or `?` where it has none. The default and JSON formats ignore it; --output_format parquet refuses it
    #[arg(help_heading = "Output", long)]
    pub fields: Option<String>,

    /// The INFO key the consequence entry of a VCF output is written under and its ##INFO
    /// header line declares. It must be a VCF INFO key (a letter or underscore, then letters,
    /// digits, underscores or dots, or 1000G); any other value, the empty string and 0
    /// included, is an argument syntax error, where Ensembl VEP accepts it, declares the given
    /// name in the header and keys the entry under it, or under CSQ for the empty string and 0
    #[arg(
        help_heading = "Output",
        long = "vcf_info_field",
        default_value = "CSQ",
        value_parser = vcf_info_key
    )]
    pub vcf_info_field: String,

    /// The vocabulary consequence terms are printed in: SO (Sequence Ontology terms), display
    /// (Ensembl display labels such as NON_SYNONYMOUS_CODING) or NCBI (labels such as missense),
    /// the SO term standing in for a consequence without a label in the chosen vocabulary.
    /// Severity ranks stay on the SO term, so the most severe term (--most_severe, JSON's
    /// most_severe_consequence) is chosen by SO rank and printed in the vocabulary; a display or
    /// NCBI JSON output is one Ensembl VEP does not produce, since it exits when it ranks a label
    /// against another term there
    #[arg(
        help_heading = "Output",
        short = 't',
        long,
        default_value = "SO",
        value_parser = ["SO", "display", "NCBI"]
    )]
    pub terms: String,

    /// Do not percent-encode the two characters Ensembl VEP escapes: HGVSp keeps its = (else %3D)
    /// and whitespace in a gene symbol is left to the VCF subfield rule, which writes it as _
    /// (else %20); implied by --json and by --output_format json
    #[arg(help_heading = "Output", long = "no_escape")]
    pub no_escape: bool,

    /// Keep an input VCF's existing annotation: its `##INFO` and `##VEP` header lines and the INFO entry under --vcf_info_field's name stay, and the new entry is appended under the same name, so the line carries both
    #[arg(help_heading = "Output", long = "keep_csq")]
    pub keep_csq: bool,

    #[arg(help_heading = "Output", long = "keep_ann")]
    pub keep_ann: bool,

    #[arg(help_heading = "Output", long)]
    pub solr: bool,

    #[arg(help_heading = "Output", long, value_name = "FORMAT")]
    pub convert: Option<String>,

    /// Annotate from the local cache, which vep-rs always does; a directory value is the cache root, as --dir_cache, when that flag is absent
    #[arg(
        help_heading = "Cache / Offline",
        long,
        num_args = 0..=1,
        default_missing_value = crate::config::CACHE_FLAG_ALONE,
        value_name = "DIR"
    )]
    pub cache: Option<String>,

    #[arg(help_heading = "Cache / Offline", long)]
    pub offline: bool,

    /// Top-level VEP data directory: the cache root when neither --dir_cache nor a --cache directory is given
    #[arg(help_heading = "Cache / Offline", long, default_value = "~/.vep")]
    pub dir: String,

    /// Cache root: the cache read is <dir_cache>/<species>[_refseq|_merged]/<cache_version>_<assembly>, as Ensembl VEP lays it out (defaults to --dir)
    #[arg(help_heading = "Cache / Offline", long = "dir_cache")]
    pub dir_cache: Option<String>,

    /// Cache version: the prefix of the cache directory under the species directory (defaults to the VEP release)
    #[arg(help_heading = "Cache / Offline", long = "cache_version")]
    pub cache_version: Option<u32>,

    /// The cache directory itself, the one holding info.json and transcripts/ (built by vep-cache-builder or converted from an Ensembl VEP cache); it is read in place of the --dir_cache layout, and only --assembly is checked against its info.json
    #[arg(
        help_heading = "Cache / Offline",
        long = "json_cache",
        visible_alias = "full_cache_dir"
    )]
    pub json_cache: Option<String>,

    /// Print the cache's source versions as sorted key<TAB>value lines and exit
    #[arg(help_heading = "Cache / Offline", long = "show_cache_info")]
    pub show_cache_info: bool,

    /// FASTA file for reference sequence lookup: plain with its .fai, or bgzipped with the .fai and .gzi that samtools faidx writes; a gzip file that is not BGZF is refused and no index is built. A directory stands for every indexed .fa and .fa.gz in it, read together. Without the flag, the first indexed .fa or .fa.gz by name in the cache directory is used
    #[arg(help_heading = "Cache / Offline", long, visible_alias = "fa")]
    pub fasta: Option<String>,

    /// Directory holding FASTA files of several species or assemblies: the .fa or .fa.gz whose name carries --species and the assembly is used, a toplevel file before any other, a bgzipped file before a plain one and the first by name of two tied on both (where Ensembl VEP keeps the directory listing's order), in place of --fasta and of the cache directory's FASTA; when none qualifies the run warns and falls back to --fasta or the cache directory's FASTA
    #[arg(help_heading = "Cache / Offline", long = "fasta_dir")]
    pub fasta_dir: Option<String>,

    /// Do not use a FASTA file found in the cache directory
    #[arg(help_heading = "Cache / Offline", long = "no_fasta")]
    pub no_fasta: bool,

    /// Chromosome synonyms file: on each line, whitespace-separated names of one chromosome, the first the reference name; an input chromosome the cache does not name is looked up through its synonyms. Replaces the cache's own chr_synonyms.txt, which is read when this flag is absent
    #[arg(help_heading = "Cache / Offline", long, value_name = "FILE")]
    pub synonyms: Option<String>,

    /// Move each insertion or deletion to its most 3' position on the transcript, over
    /// a 1,000-base window of the reference, before its consequence terms, positions,
    /// codons and amino acids are computed (1), or annotate it where the input placed it
    /// (0); needs a FASTA and stops the run without one. The HGVS notations shift under
    /// --shift_hgvs, not this flag. The shift runs along the genome only: a transcript
    /// whose cached sequence is edited away from the genome is not shifted again along
    /// its own sequence
    #[arg(
        help_heading = "Cache / Offline",
        long = "shift_3prime",
        value_parser = ["0", "1"]
    )]
    pub shift_3prime: Option<String>,

    /// Move each insertion or deletion to its most 3' position on the transcript before
    /// its HGVSc and HGVSp are built (1, the default), or describe it where the input
    /// placed it (0), which also sets --shift_3prime to 0 and cannot be combined with
    /// --shift_3prime or --shift_genomic
    #[arg(
        help_heading = "Cache / Offline",
        long = "shift_hgvs",
        value_parser = ["0", "1"]
    )]
    pub shift_hgvs: Option<String>,

    /// Under --shift_3prime 1, print the shifted position in Location and the shifted
    /// sequence in Allele on every row (1), intergenic rows moving towards higher
    /// coordinates, or keep the input position (0); without --shift_3prime 1 nothing
    /// moves
    #[arg(
        help_heading = "Cache / Offline",
        long = "shift_genomic",
        value_parser = ["0", "1"]
    )]
    pub shift_genomic: Option<String>,

    /// Add the SHIFT_LENGTH column: under --shift_3prime 1, the bases each insertion
    /// or deletion moved along the transcript, 0 on every other transcript row of a
    /// sequence variant and none on a structural record's rows
    #[arg(help_heading = "Cache / Offline", long = "shift_length")]
    pub shift_length: bool,

    #[arg(help_heading = "Cache / Offline", long)]
    pub sereal: bool,

    /// Species: names the cache directory under --dir_cache; the Latin name only, with no alias such as human
    #[arg(
        help_heading = "Species / Assembly",
        short = 's',
        long,
        default_value = "homo_sapiens"
    )]
    pub species: String,

    /// Genome assembly (e.g., GRCh38): chooses among several caches of one version and must equal the cache's own name exactly; hg19 and hg38 are read as GRCh37 and GRCh38
    #[arg(help_heading = "Species / Assembly", short = 'a', long)]
    pub assembly: Option<String>,

    /// Annotate against GRCh37: sets --assembly GRCh37, over an explicit --assembly, as Ensembl VEP's option set does
    #[arg(help_heading = "Species / Assembly", long)]
    pub grch37: bool,

    #[arg(help_heading = "Database", long)]
    pub database: bool,

    #[arg(help_heading = "Database", long, value_name = "HOST")]
    pub host: Option<String>,

    #[arg(help_heading = "Database", long, value_name = "PORT")]
    pub port: Option<String>,

    #[arg(help_heading = "Database", short = 'u', long, value_name = "USER")]
    pub user: Option<String>,

    #[arg(
        help_heading = "Database",
        long,
        visible_alias = "pass",
        value_name = "PASSWORD"
    )]
    pub password: Option<String>,

    #[arg(help_heading = "Database", long, value_name = "FILE")]
    pub registry: Option<String>,

    #[arg(help_heading = "Database", long = "db_version", value_name = "VERSION")]
    pub db_version: Option<String>,

    #[arg(help_heading = "Database", long)]
    pub genomes: bool,

    #[arg(
        help_heading = "Database",
        long = "is_multispecies",
        value_name = "0|1"
    )]
    pub is_multispecies: Option<String>,

    /// Number of variants to process per batch
    #[arg(help_heading = "Runtime", long = "buffer_size", default_value_t = 5000)]
    pub buffer_size: usize,

    /// A record whose end minus start exceeds this many bases is oversize, a sequence
    /// record included; -1 lifts the limit, as in Ensembl VEP 116. With a FASTA, a
    /// <CNV:TR> record whose repeat sequence is longer than this stays a symbolic
    /// tandem repeat rather than expanding to its alleles (under -1 every one does,
    /// and 0 reads as 5000), as in Ensembl VEP 116. Ensembl VEP 116 drops an
    /// oversize record before annotation and writes nothing for it in any format;
    /// vep-rs annotates it: rows in the default, tab and parquet outputs, its VCF
    /// line without consequences, no JSON object (docs/output-formats.md and
    /// docs/intended-divergences.md describe the divergence). A region, SPDI or
    /// HGVS transcript record whose span exceeds the limit is skipped with Ensembl
    /// VEP's warning, as every oversize record is in Ensembl VEP 116; the skip is
    /// made as the line is parsed, so it holds under --dont_skip too, where Ensembl
    /// VEP keeps the feature (--dont_skip names the records it cannot keep). An HGVS
    /// genomic notation's end minus its start must be smaller than this value, 0
    /// reading as 5000 and -1 as -1 itself, so under -1 a notation written in
    /// genome order is skipped, as Ensembl VEP 116 skips it
    #[arg(
        help_heading = "Runtime",
        long = "max_sv_size",
        default_value_t = DEFAULT_MAX_SV_SIZE,
        allow_negative_numbers = true,
        value_parser = clap::value_parser!(i64).range(-1..)
    )]
    pub max_sv_size: i64,

    /// Number of parallel annotation worker threads (Ensembl VEP forks processes; vep-rs runs threads and is parallel by default).
    /// Default 0 = auto: uses the number of available logical CPUs, capped at
    /// 32 (one reader thread feeds the workers, so extra workers wait at the
    /// batch barrier). Set to 1 to force single-threaded annotation, or any
    /// positive value to pin the worker count; user-supplied values are not
    /// capped.
    #[arg(
        help_heading = "Runtime",
        long,
        visible_alias = "threads",
        default_value_t = 0
    )]
    pub fork: usize,

    /// Number of BGZF decompression worker threads for a BGZF-compressed input
    /// (a plain gzip input is decompressed on one thread).
    /// Default 0: uses the number of available logical CPUs, capped at 10.
    /// `auto` chooses the count from the input's inflated bytes per record: a
    /// BGZF input with a .tbi or .csi index beside it whose records inflate to
    /// at least 8,000 bytes gets 32 workers on an aarch64 host, capped at the
    /// available logical CPUs; any other input, and every input on another
    /// architecture, gets the 0 rule. `auto` is opt-in, never the default. A
    /// positive value pins the worker count, e.g. `--decompression_threads=1`
    /// for single-threaded; a value above the number of available logical CPUs
    /// is lowered to that number.
    #[arg(
        help_heading = "Runtime",
        long = "decompression_threads",
        default_value = "0",
        value_parser = crate::config::decompression_threads_text
    )]
    pub decompression_threads: String,

    /// How the batches move through the pipeline. `streaming` hands each batch to
    /// the worker pool as independent tasks that annotate and render their own
    /// records (and parse them, for a VCF input) and takes the next batch without
    /// waiting for them; the writer puts the rendered output back in input order.
    /// `batched` parses, annotates and renders each batch on the pool as a whole
    /// before the next one starts.
    /// A run with plugins, or one reading known variants from the cache's variation
    /// files or from the --check_svs file, uses `batched` whatever this says: plugins
    /// see a whole batch at once, and each of those sources prefetches one batch's
    /// regions at a time
    #[arg(
        help_heading = "Runtime",
        long = "pipeline_mode",
        default_value = "streaming",
        hide = true,
        value_parser = ["batched", "streaming"]
    )]
    pub pipeline_mode: String,

    /// Batches a pipeline stage may hold ahead of the next one, at most 1024; 0 is
    /// the built-in depth
    #[arg(
        help_heading = "Runtime",
        long = "pipeline_depth",
        default_value_t = 0,
        hide = true,
        value_parser = clap::value_parser!(u64).range(0..=1024)
    )]
    pub pipeline_depth: u64,

    /// Input records per pool task in both pipeline modes (per render task on the
    /// plugin path), at most 1,000,000; 0 is the built-in sizing
    #[arg(
        help_heading = "Runtime",
        long = "render_task_records",
        default_value_t = 0,
        hide = true,
        value_parser = clap::value_parser!(u64).range(0..=1_000_000)
    )]
    pub render_task_records: u64,

    /// How rendered output reaches the file or standard output. `buffered` collects
    /// it in a buffer of --output_buffer_size bytes; `direct` writes each rendered
    /// chunk with one call; `vectored` writes the chunks that are ready together with
    /// one vectored call. Every strategy writes the same uncompressed bytes
    #[arg(
        help_heading = "Runtime",
        long = "writer_strategy",
        default_value = "buffered",
        hide = true,
        value_parser = ["buffered", "direct", "vectored"]
    )]
    pub writer_strategy: String,

    /// Bytes of the `buffered` writer strategy's buffer, at most 1073741824 (1 GiB);
    /// 0 is the built-in size. The other strategies ignore it
    #[arg(
        help_heading = "Runtime",
        long = "output_buffer_size",
        default_value_t = 0,
        hide = true,
        value_parser = clap::value_parser!(u64).range(0..=1_073_741_824)
    )]
    pub output_buffer_size: u64,

    /// Print the resolved runtime knobs and the input's shape as name<TAB>value
    /// lines and exit without annotating
    #[arg(help_heading = "Runtime", long = "resolved_config", hide = true)]
    pub resolved_config: bool,

    /// Upstream and downstream distance in bp: one integer, or the pair "up,down"
    /// (a downstream part written 0 reads as the upstream value)
    #[arg(help_heading = "Runtime", long, default_value = "5000")]
    pub distance: String,

    /// Shorthand for Ensembl VEP's --everything set: --sift b, --polyphen b, --ccds, --hgvs, --symbol,
    /// --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg,
    /// --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris,
    /// --variant_class, --gene_phenotype and --mirna, each overriding an explicit value as in Ensembl VEP.
    /// Every column of Ensembl VEP's --everything except the regulatory columns: --regulatory is
    /// switched on and annotates nothing, so no RegulatoryFeature or MotifFeature row is written
    #[arg(help_heading = "Output Control", short = 'e', long)]
    pub everything: bool,

    /// Add gene symbol to output
    #[arg(help_heading = "Output Control", long)]
    pub symbol: bool,

    /// Add transcript biotype to output
    #[arg(help_heading = "Output Control", long)]
    pub biotype: bool,

    /// Flag canonical transcripts
    #[arg(help_heading = "Output Control", long)]
    pub canonical: bool,

    /// Add the MANE, MANE_SELECT and MANE_PLUS_CLINICAL columns: the MANE sets the transcript
    /// belongs to and the RefSeq accession each set names (switches on --mane_select)
    #[arg(help_heading = "Output Control", long)]
    pub mane: bool,

    /// Add the MANE and MANE_SELECT columns alone: MANE lists MANE_Select on a MANE Select
    /// transcript and MANE_SELECT carries its RefSeq accession
    #[arg(help_heading = "Output Control", long = "mane_select")]
    pub mane_select: bool,

    /// Add the TSL column: the digit of the transcript's support level, empty for tslNA (the GRCh37
    /// gene set carries no level)
    #[arg(help_heading = "Output Control", long)]
    pub tsl: bool,

    /// Add the APPRIS column with the grade abbreviated as Ensembl VEP prints it: principal1 is P1,
    /// alternative2 is A2, principalM is PM
    #[arg(help_heading = "Output Control", long)]
    pub appris: bool,

    /// Add the CCDS column: the transcript's CCDS identifier (a native cache carries it from the
    /// GFF3; no RefSeq or merged cache is built natively)
    #[arg(help_heading = "Output Control", long)]
    pub ccds: bool,

    /// Add the ENSP column: the Ensembl protein identifier
    #[arg(help_heading = "Output Control", long)]
    pub protein: bool,

    /// Add the translation version to the Ensembl protein ID (ENSP00000286808.3); needs --protein
    #[arg(help_heading = "Output Control", long = "protein_version")]
    pub protein_version: bool,

    /// Add the SWISSPROT, TREMBL, UNIPARC and UNIPROT_ISOFORM columns: the translation's UniProt
    /// cross-references (a native cache carries them only when built with --core-dumps)
    #[arg(help_heading = "Output Control", long)]
    pub uniprot: bool,

    /// Add the RefSeq column: the RefSeq mRNA identifiers cross-referenced to the transcript (a
    /// native cache carries them only when built with --core-dumps)
    #[arg(help_heading = "Output Control", long = "xref_refseq")]
    pub xref_refseq: bool,

    /// Append the transcript's version to the Feature column (ENST00000366667.4) unless the
    /// identifier already ends in one
    #[arg(help_heading = "Output Control", long = "transcript_version")]
    pub transcript_version: bool,

    /// Append the gene's version to the Gene column (ENSG00000142168.14) unless the identifier
    /// already ends in one
    #[arg(help_heading = "Output Control", long = "gene_version")]
    pub gene_version: bool,

    /// Add HGVS coding and protein notations; needs a FASTA named with --fasta or found in the cache directory and stops the run without one
    #[arg(help_heading = "Output Control", long)]
    pub hgvs: bool,

    /// Add the HGVS genomic notation (HGVSg); needs a FASTA named with --fasta or found in the cache directory and stops the run without one
    #[arg(help_heading = "Output Control", long)]
    pub hgvsg: bool,

    /// Name the chromosome in HGVSg by its RefSeq accession: the first NC_ name among
    /// the input name's own synonyms in the chromosome synonyms table (one hop: only a
    /// line pairing the name as written with the accession renames it); an input name
    /// with no NC_ synonym keeps its own name, where Ensembl VEP stops with an error
    #[arg(help_heading = "Output Control", long = "hgvsg_use_accession")]
    pub hgvsg_use_accession: bool,

    /// Add the SPDI column: each allele in Sequence Position Deletion Insertion
    /// notation, the position 0-based, the deleted bases (their count when more than
    /// 20) and the inserted bases, as the input wrote the record; the chromosome is
    /// named as the cache names it, and on a row overlapping no transcript as the
    /// FASTA names it. Ensembl VEP needs a FASTA for those rows; this reads none. A
    /// record whose reference allele ends in a character outside `ACGTN-` (`AR/CA`) is
    /// annotated with no SPDI, where Ensembl VEP stops the run on it
    #[arg(help_heading = "Output Control", long)]
    pub spdi: bool,

    /// Add the GA4GH_VRS column: the allele's SPDI as a GA4GH Variation Representation
    /// Specification Allele object on the sequence refseq:<NC_ accession>, the accession
    /// being the first NC_ synonym of the cache's chromosome name in the chromosome
    /// synonyms table (Ensembl VEP writes the object only where its SPDI is named by an
    /// NC_ accession: with a database, or on a row overlapping no transcript whose FASTA
    /// record carries that name; a chromosome with no NC_ synonym gets none);
    /// JSON-encoded in the text outputs, nested in the JSON output. Cannot be combined
    /// with --vcf
    #[arg(help_heading = "Output Control", long = "ga4gh_vrs")]
    pub ga4gh_vrs: bool,

    /// Write HGVSp in the predicted form, the description after p. in parentheses
    /// (ENSP00000286808.3:p.(Arg97Ter)); the = of a synonymous description is still
    /// escaped as %3D in the text formats
    #[arg(help_heading = "Output Control", long = "hgvsp_use_prediction")]
    pub hgvsp_use_prediction: bool,

    /// Leave the translation version off the HGVSp reference name
    /// (ENSP00000286808:p.Arg97Ter); HGVSc keeps its version
    #[arg(help_heading = "Output Control", long = "remove_hgvsp_version")]
    pub remove_hgvsp_version: bool,

    /// Add the SIFT column on every single-residue substitution inside the coding sequence
    /// (start_lost included, a stop change never): p the prediction, s the score, b both as
    /// prediction(score). Needs a cache whose metadata names a SIFT version; a cache without
    /// one stops the run, or under --everything or --safe drops the flag with an INFO line (a
    /// native cache carries the matrices only when built with --include-predictions)
    #[arg(help_heading = "Output Control", long, value_parser = ["s", "p", "b"])]
    pub sift: Option<String>,

    /// Add the PolyPhen column on every single-residue substitution inside the coding sequence
    /// (start_lost included, a stop change never): p the prediction, s the score, b both as
    /// prediction(score). Needs a cache whose metadata names a PolyPhen version; a cache
    /// without one stops the run, or under --everything or --safe drops the flag with an INFO
    /// line (a native cache carries the matrices only when built with --include-predictions)
    #[arg(help_heading = "Output Control", long, value_parser = ["s", "p", "b"])]
    pub polyphen: Option<String>,

    /// Read the PolyPhen column from the HumDiv-trained matrix instead of the default HumVar one;
    /// a cache without the HumDiv matrix prints an empty column with no error, as in Ensembl VEP
    #[arg(help_heading = "Output Control", long)]
    pub humdiv: bool,

    #[arg(help_heading = "Output Control", long, value_name = "CONDEL")]
    pub condel: Option<String>,

    /// Add the NEAREST column: the protein-coding transcription start site(s) nearest each
    /// variant as a transcript stable id, gene stable id or gene symbol; the site above the
    /// variant and the site below are measured to the nearer end of the variant, both are
    /// reported when they tie, and every transcript starting at a reported site is listed in
    /// stable-id order where Ensembl VEP lists one
    #[arg(
        help_heading = "Output Control",
        long,
        value_name = "TYPE",
        value_parser = ["transcript", "gene", "symbol"]
    )]
    pub nearest: Option<String>,

    /// Add exon/intron numbers
    #[arg(help_heading = "Output Control", long)]
    pub numbers: bool,

    /// Add the DOMAINS column: on a coding row, the protein features whose residues the variant
    /// overlaps, each as label:identifier (Pfam:PF00080, Gene3D:1.10.510.10; a native cache
    /// carries the features only when built with --core-dumps)
    #[arg(help_heading = "Output Control", long)]
    pub domains: bool,

    /// Accepted for VEP compatibility; regulatory-feature annotation is not implemented and the flag warns
    #[arg(help_heading = "Output Control", long)]
    pub regulatory: bool,

    #[arg(
        help_heading = "Output Control",
        long = "regulatory_gff",
        value_name = "FILE"
    )]
    pub regulatory_gff: Option<String>,

    #[arg(help_heading = "Output Control", long = "extended_promoters")]
    pub extended_promoters: bool,

    /// Add variant class from SO
    #[arg(help_heading = "Output Control", long = "variant_class")]
    pub variant_class: bool,

    /// Add the OverlapBP and OverlapPC columns to the tab, VCF and Parquet outputs: the
    /// bases of a structural variant inside the transcript and the percentage of the
    /// transcript they cover (transcript rows only; regulatory features are not annotated)
    #[arg(help_heading = "Output Control", long)]
    pub overlaps: bool,

    /// Check for co-located known variants
    #[arg(help_heading = "Output Control", long = "check_existing")]
    pub check_existing: bool,

    /// Report a known variant as co-located when its start and end equal the input variant's, whatever its alleles; such a variant prints on every allele's row
    #[arg(help_heading = "Output Control", long = "no_check_alleles")]
    pub no_check_alleles: bool,

    /// Never report a known variant whose alleles the cache does not hold (COSMIC and HGMD records, dbSNP records without variation); without the flag such a variant is reported when its start and end equal the input variant's
    #[arg(help_heading = "Output Control", long = "exclude_null_alleles")]
    pub exclude_null_alleles: bool,

    /// Report known variants whose failed value is at most this: 0 excludes variants that failed Ensembl's quality checks, 1 includes them
    #[arg(
        help_heading = "Output Control",
        long,
        value_name = "INT",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    pub failed: i64,

    /// Report the ClinVar significance of the row's own allele where a co-located variant carries allele-specific entries (1, the default); 0 reports every variant's significance whatever the allele
    #[arg(
        help_heading = "Output Control",
        long = "clin_sig_allele",
        value_name = "INT",
        default_value_t = 1,
        allow_negative_numbers = true
    )]
    pub clin_sig_allele: i64,

    /// Add the ClinVar somatic classification of co-located variants (CLINVAR_SOMATIC_CLASSIFICATION), the cache's string as written
    #[arg(
        help_heading = "Output Control",
        long = "clinvar_somatic_classification"
    )]
    pub clinvar_somatic_classification: bool,

    /// Add the known synonyms of co-located variants and their sources (VAR_SYNONYMS); switches on --check_existing
    #[arg(help_heading = "Output Control", long = "var_synonyms")]
    pub var_synonyms: bool,

    /// Report the known structural variants overlapping each variant (SV) from Ensembl's indexed structural-variant VCF (homo_sapiens_structural_variations.vcf.gz with its .csi or .tbi beside it), Ensembl's dump of the database Ensembl VEP reads them from; Ensembl VEP takes the flag without a value and refuses it offline
    #[arg(
        help_heading = "Output Control",
        long = "check_svs",
        value_name = "FILE"
    )]
    pub check_svs: Option<String>,

    /// Keep only the rows whose allele lies in a transcript's coding sequence (a
    /// CDS position or a coding consequence term, a variant covering a translated
    /// transcript whole included); intergenic and non-coding rows go, and a record
    /// left with no row is written as an unannotated one
    #[arg(help_heading = "Output Control", long = "coding_only")]
    pub coding_only: bool,

    /// Write no row for a record that overlaps no transcript: no intergenic row
    /// in the default and tab formats, the VCF line without a CSQ entry, and in
    /// JSON a most_severe_consequence of "?" with no intergenic_consequences
    #[arg(help_heading = "Output Control", long = "no_intergenic")]
    pub no_intergenic: bool,

    #[arg(help_heading = "Output Control", long)]
    pub hgnc: bool,

    /// Add global allele frequency from 1000 Genomes Phase 3
    #[arg(help_heading = "Allele Frequencies", long)]
    pub af: bool,

    /// Add continental allele frequencies from 1000 Genomes
    #[arg(help_heading = "Allele Frequencies", long = "af_1kg")]
    pub af_1kg: bool,

    /// Add gnomAD exome allele frequencies
    #[arg(help_heading = "Allele Frequencies", long = "af_gnomade")]
    pub af_gnomade: bool,

    /// Add gnomAD exome allele frequencies under the older spelling of --af_gnomade; the two spellings cannot be combined, and this one alone fails on a cache without gnomAD exome data
    #[arg(help_heading = "Allele Frequencies", long = "af_gnomad")]
    pub af_gnomad: bool,

    /// Add gnomAD genome allele frequencies
    #[arg(help_heading = "Allele Frequencies", long = "af_gnomadg")]
    pub af_gnomadg: bool,

    /// Add highest allele frequency across populations
    #[arg(help_heading = "Allele Frequencies", long = "max_af")]
    pub max_af: bool,

    /// Add PubMed IDs for co-located variants
    #[arg(help_heading = "Allele Frequencies", long)]
    pub pubmed: bool,

    /// Rewrite the AFR, AMR, ASN and EUR frequency values of each co-located variant as a bare minor-allele frequency: the allele is stripped, the first value kept and a value above 0.5 inverted; the allele-keyed frequency columns then find no allele in them and print empty
    #[arg(help_heading = "Allele Frequencies", long = "old_maf")]
    pub old_maf: bool,

    /// Filter the input by the frequency of co-located known variants in the population --freq_pop names, against the threshold --freq_freq on the side --freq_gt_lt gives, dropping the alleles that pass or keeping the records whose alleles all pass as --freq_filter says; switches on --check_existing and adds the FREQS field (the frequencies the filter read, NA where an allele has none)
    #[arg(help_heading = "Allele Frequencies", long = "check_frequency")]
    pub check_frequency: bool,

    /// The population --check_frequency reads: 1KG_ALL (the 1000 Genomes AF column of the cache), a frequency column of the cache in any case (AF, AFR, gnomADe_NFE; the gnomAD genome columns are refused, as Ensembl VEP refuses them), or 1KG_ followed by such a column
    #[arg(
        help_heading = "Allele Frequencies",
        long = "freq_pop",
        value_name = "POP",
        default_value = "1KG_ALL"
    )]
    pub freq_pop: String,

    /// The frequency --check_frequency compares against
    #[arg(
        help_heading = "Allele Frequencies",
        long = "freq_freq",
        value_name = "FLOAT",
        default_value_t = 0.01,
        value_parser = real_number,
        allow_negative_numbers = true
    )]
    pub freq_freq: f64,

    /// Whether an allele passes --check_frequency with a frequency at or above (gt) or at or below (lt) --freq_freq
    #[arg(
        help_heading = "Allele Frequencies",
        long = "freq_gt_lt",
        value_name = "gt|lt",
        default_value = "gt",
        value_parser = ["gt", "lt"]
    )]
    pub freq_gt_lt: String,

    /// Whether --check_frequency drops the alleles that pass (exclude) or keeps a record only when every allele passes (include)
    #[arg(
        help_heading = "Allele Frequencies",
        long = "freq_filter",
        value_name = "exclude|include",
        default_value = "exclude",
        value_parser = ["exclude", "include"]
    )]
    pub freq_filter: String,

    /// Drop the alleles with a 1000 Genomes frequency of 0.01 or more: --check_frequency with --freq_pop 1KG_ALL --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude
    #[arg(help_heading = "Allele Frequencies", long = "filter_common")]
    pub filter_common: bool,

    /// Add the GENE_PHENO column: 1 on a transcript whose gene the cache associates with a
    /// phenotype, disease or trait, empty otherwise
    #[arg(help_heading = "Allele Frequencies", long = "gene_phenotype")]
    pub gene_phenotype: bool,

    /// Pick one consequence per variant
    #[arg(help_heading = "Filtering / Picking", long)]
    pub pick: bool,

    /// Pick one consequence per variant allele
    #[arg(help_heading = "Filtering / Picking", long = "pick_allele")]
    pub pick_allele: bool,

    /// Pick one consequence per gene
    #[arg(help_heading = "Filtering / Picking", long = "per_gene")]
    pub per_gene: bool,

    /// Pick one consequence per allele per gene
    #[arg(help_heading = "Filtering / Picking", long = "pick_allele_gene")]
    pub pick_allele_gene: bool,

    /// Report only the single most severe consequence
    #[arg(help_heading = "Filtering / Picking", long = "most_severe")]
    pub most_severe: bool,

    /// Report a summary of all consequences
    #[arg(help_heading = "Filtering / Picking", long)]
    pub summary: bool,

    /// Flag (but don't remove) non-picked consequences
    #[arg(help_heading = "Filtering / Picking", long = "flag_pick")]
    pub flag_pick: bool,

    /// Flag non-picked consequences per allele
    #[arg(help_heading = "Filtering / Picking", long = "flag_pick_allele")]
    pub flag_pick_allele: bool,

    /// Flag non-picked consequences per allele per gene
    #[arg(help_heading = "Filtering / Picking", long = "flag_pick_allele_gene")]
    pub flag_pick_allele_gene: bool,

    /// The criteria the pick family compares, in order, as a comma-separated list; a
    /// shorter list compares its criteria and then keeps the first transcript in
    /// stable-id order, and an empty value keeps it at once. Default: mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq
    #[arg(
        help_heading = "Filtering / Picking",
        long = "pick_order",
        value_parser = PickOrderParser
    )]
    pub pick_order: Option<PickOrder>,

    /// Annotate only the transcripts that pass a filter, in Ensembl VEP's filter
    /// grammar: units of `field operator value` joined by `and` and `or`, grouped
    /// by parentheses and inverted by `not`; repeated flags are joined by `and`.
    /// The fields are the cache's transcript keys, named exactly, in another
    /// case, without the leading `_`, or by a unique prefix: `stable_id`,
    /// `version`, `biotype`, `start`, `end`, `strand`, `is_canonical`, `source`,
    /// `dbID`, `cdna_coding_start`, `cdna_coding_end`, `_gene_stable_id`,
    /// `_gene_version`, `_gene_symbol`, `_gene_symbol_source`, `_gene_hgnc_id`,
    /// `_ccds`, `_refseq`, `_protein`, `_swissprot`, `_trembl`, `_uniparc`,
    /// `_uniprot_isoform`, `_gene_phenotype`; the object keys `translation`
    /// (coding transcripts), `slice`, `attributes`, `_gene`, `_trans_exon_array`,
    /// `_variation_effect_feature_cache` and `_vep_lazy_loaded`, which are defined
    /// and read meaningfully by `exists` and `nex` alone; and `description`,
    /// `coding_region_start` and `coding_region_end`, which the cache leaves
    /// undefined.
    /// `_source_cache` is absent, as it is from every cache but a merged one. The
    /// operators are `is` (`eq`, `=`), `ne` (`!=`), `>`, `<`, `>=`, `<=`, `exists`
    /// (`defined`, the default), `nex`, `match` (`matches`, `regex`, `re`), `nre`,
    /// `in` (a comma-separated list or a file with one value per line) and
    /// `is_child`. `match` and `nre` take the `regex` crate's dialect,
    /// case-insensitive and unanchored: a back-reference or an escape the crate
    /// lacks is an error naming the pattern. A record whose every transcript is
    /// filtered out is intergenic.
    #[arg(
        help_heading = "Filtering / Picking",
        long = "transcript_filter",
        value_name = "FILTER"
    )]
    pub transcript_filter: Vec<String>,

    /// Drop RefSeq's predicted transcript models, the `XM_` and `XR_` accessions,
    /// by adding the filter `not stable_id match ^X._` after the --transcript_filter
    /// values, as Ensembl VEP's option set does. An Ensembl cache carries no
    /// predicted model, so the flag changes nothing there; it takes effect on a
    /// converted RefSeq or merged cache (no RefSeq or merged cache is built
    /// natively).
    #[arg(help_heading = "Filtering / Picking", long = "exclude_predicted")]
    pub exclude_predicted: bool,

    /// Add allele number from input VCF
    #[arg(help_heading = "Filtering / Picking", long = "allele_number")]
    pub allele_number: bool,

    /// Show reference allele in output
    #[arg(help_heading = "Filtering / Picking", long = "show_ref_allele")]
    pub show_ref_allele: bool,

    /// Add the UPLOADED_ALLELE column: the record's allele string as uploaded, REF and every ALT
    /// joined by /
    #[arg(help_heading = "Filtering / Picking", long = "uploaded_allele")]
    pub uploaded_allele: bool,

    /// Add the AMBIGUITY column: the IUPAC ambiguity code of the record's alleles (R for A/G);
    /// absent on an indel or any other allele set without a code
    #[arg(help_heading = "Filtering / Picking", long)]
    pub ambiguity: bool,

    /// Trim every bi-allelic record's alleles from both ends and split a multi-allelic
    /// record whose alleles trim into one allele per row, each on its own trimmed span,
    /// as Ensembl VEP minimises; switches on --allele_number as in Ensembl VEP
    #[arg(help_heading = "Filtering / Picking", long)]
    pub minimal: bool,

    /// Append /<total length> to the cDNA, CDS and protein positions: the transcript's cDNA
    /// length, its coding sequence length and its peptide length. JSON keeps Ensembl VEP's
    /// shape: a single base prints the string "428/1200" as both cdna_start and cdna_end, and a
    /// range prints its start as both
    #[arg(help_heading = "Filtering / Picking", long = "total_length")]
    pub total_length: bool,

    /// Load plugin(s) - can be specified multiple times
    #[arg(help_heading = "Plugins", long)]
    pub plugin: Vec<String>,

    /// Custom annotation source(s) (accepted for VEP compatibility; not implemented, sources are not loaded)
    #[arg(help_heading = "Plugins", long)]
    pub custom: Vec<String>,

    #[arg(help_heading = "Plugins", long = "custom_suppress_filter")]
    pub custom_suppress_filter: bool,

    /// Directory the plugin libraries are looked up in; default <--dir>/Plugins. A
    /// directory that does not exist is warned about when --plugin names a plugin
    /// (an error under --safe). A value holding a glob pattern (`*`, `?`, `[...]`,
    /// `{a,b}`) or a `~user` prefix is read as one literal path; Ensembl VEP expands
    /// it with `glob`
    #[arg(help_heading = "Plugins", long = "dir_plugins", value_name = "DIR")]
    pub dir_plugins: Option<String>,

    /// Make a plugin that cannot be loaded, or a missing plugins directory, an error
    /// instead of a warning, and skip the check of incompatible flag pairs
    #[arg(help_heading = "Plugins", long)]
    pub safe: bool,

    /// Reads the cache under the species directory suffixed _refseq (the transcripts annotated are that cache's contents), and resolves an HGVS gene name among the RefSeq transcripts before the Ensembl ones, as Ensembl VEP's --refseq orders its core groups; under --json_cache only that order applies
    #[arg(help_heading = "Transcript Sets", long)]
    pub refseq: bool,

    /// Reads the cache under the species directory suffixed _merged (the transcripts annotated are that cache's contents); no effect under --json_cache
    #[arg(help_heading = "Transcript Sets", long)]
    pub merged: bool,

    /// Annotate against the GENCODE basic transcript set alone: a transcript without a
    /// gencode_basic attribute is dropped as the cache loads, before the transcript index
    /// is built
    #[arg(help_heading = "Transcript Sets", long = "gencode_basic")]
    pub gencode_basic: bool,

    /// Annotate against the GENCODE primary transcript set alone: a transcript without a
    /// gencode_primary attribute is dropped as the cache loads; the set exists for human
    /// GRCh38 only, so --assembly naming another assembly is warned about
    #[arg(help_heading = "Transcript Sets", long = "gencode_primary")]
    pub gencode_primary: bool,

    /// Add the GENCODE_PRIMARY column: 1 on a transcript of the GENCODE primary set, empty
    /// otherwise; the set exists for human GRCh38 only, so --assembly naming another
    /// assembly is warned about
    #[arg(help_heading = "Transcript Sets", long = "flag_gencode_primary")]
    pub flag_gencode_primary: bool,

    /// Read each row's reference allele from the transcript instead of the input, and add
    /// GIVEN_REF (the input's reference allele) and USED_REF (the row's) to the output. The
    /// row's allele is the transcript's spliced cDNA over the variant where the variant lies
    /// wholly inside the exons, `-` for an insertion, and the input's allele elsewhere; an
    /// alternate allele equal to the transcript's sequence takes the terms an allele equal to
    /// the reference takes, and the reference row of --process_ref_homs is the transcript's
    /// allele. Needs --fasta, which supplies a non-coding transcript's cDNA; --use_given_ref
    /// switches it off. Not read: the edited transcript sequences of a BAM-edited (RefSeq)
    /// cache, which Ensembl VEP reads with this flag switched on by the cache itself
    #[arg(help_heading = "Transcript Sets", long = "use_transcript_ref")]
    pub use_transcript_ref: bool,

    #[arg(help_heading = "Transcript Sets", long = "use_given_ref")]
    pub use_given_ref: bool,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written, with or without it)
    #[arg(help_heading = "Statistics", long = "no_stats")]
    pub no_stats: bool,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written)
    #[arg(help_heading = "Statistics", long = "stats_file", visible_alias = "sf")]
    pub stats_file: Option<String>,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written)
    #[arg(help_heading = "Statistics", long = "stats_text")]
    pub stats_text: bool,

    /// Verbose output (debug-level logging)
    #[arg(help_heading = "Verbosity", short = 'v', long)]
    pub verbose: bool,

    /// Suppress the status and progress lines and the standard-error copy of each
    /// warning; the warnings file is still written
    #[arg(help_heading = "Verbosity", short = 'q', long)]
    pub quiet: bool,

    /// Keep and annotate a record that fails a validation check (an unknown chromosome, an invalid
    /// allele string, a reference allele --check_ref rejects, a chromosome --chr leaves out) instead of
    /// skipping it; the skip warning is still written, except for a chromosome --chr leaves out, which
    /// is dropped without a warning and kept without one. Two records are skipped under the flag too,
    /// each with its warning, where Ensembl VEP keeps the feature: a record whose start is past its end
    /// plus one (Ensembl VEP's run then ends on the first transcript the record reaches), and a region,
    /// SPDI or HGVS transcript record whose span exceeds --max_sv_size, refused as its line is parsed.
    /// A line that cannot be parsed is skipped either way
    #[arg(help_heading = "Miscellaneous", long = "dont_skip")]
    pub dont_skip: bool,

    /// Keep a VCF record whose first ALT is `.` (a reference-only call): the VCF output carries its line with
    /// any earlier CSQ removed unless --keep_csq, the JSON output writes its object with no consequence,
    /// and the default and tab outputs write no row for it; under --individual a copy whose genotype names
    /// a later ALT of the record is annotated as that allele
    #[arg(help_heading = "Miscellaneous", long = "allow_non_variant")]
    pub allow_non_variant: bool,

    /// Add the miRNA column: the secondary-structure elements of a miRNA transcript the variant
    /// overlaps, as the SO terms miRNA_stem and miRNA_loop (a native cache carries the structure
    /// only when built with --core-dumps)
    #[arg(help_heading = "Miscellaneous", long)]
    pub mirna: bool,

    #[arg(help_heading = "Miscellaneous", long)]
    pub lrg: bool,

    #[arg(help_heading = "Miscellaneous", long, value_name = "DIR")]
    pub tmpdir: Option<String>,
}

/// A real number as Ensembl VEP's option parser reads one (Getopt::Long `=f`,
/// `Config.pm:138`): an optional sign, digits or a leading `.`, a fraction of at
/// least one digit and an optional exponent, with `_` allowed among the digits
/// and dropped, so `.5`, `1e-2`, `+0.5` and `1_0` pass and `nan`, `inf` and `1.`
/// are refused as the option's value.
fn real_number(text: &str) -> Result<f64, String> {
    const NOT_A_REAL_NUMBER: &str = "real number expected";
    let unsigned = text.strip_prefix(['+', '-']).unwrap_or(text);
    let mantissa = unsigned
        .split_once(['e', 'E'])
        .map_or(unsigned, |(mantissa, _exponent)| mantissa);
    let has_real_number_shape = mantissa.starts_with(|c: char| c.is_ascii_digit() || c == '.')
        && mantissa
            .split_once('.')
            .is_none_or(|(_, fraction)| !fraction.is_empty());
    if !has_real_number_shape {
        return Err(NOT_A_REAL_NUMBER.to_string());
    }
    text.replace('_', "")
        .parse()
        .map_err(|_| NOT_A_REAL_NUMBER.to_string())
}

/// A VCF INFO key, `[A-Za-z_][0-9A-Za-z_.]*` or `1000G`, the one key the VCF
/// specification reserves outside that grammar: any other `--vcf_info_field`
/// value would write an INFO entry no VCF reader can key. Ensembl VEP takes any
/// string and, for one Perl reads as false (the empty string, `0`), prints it in
/// the `##INFO` line (`OutputFactory/VCF.pm:268`) but keys the entry, and strips
/// an input's earlier entry, under `CSQ` (`:219`, `:312` `|| 'CSQ'`).
fn vcf_info_key(text: &str) -> Result<String, String> {
    let mut chars = text.chars();
    let is_key = chars
        .next()
        .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.');
    if is_key || text == "1000G" {
        Ok(text.to_string())
    } else {
        Err(
            "a VCF INFO key is expected: a letter or underscore, then letters, digits, \
             underscores or dots, or 1000G"
                .to_string(),
        )
    }
}

/// The `--pick_order` criteria as Ensembl VEP reads the option: the fields of
/// Perl's `split(',', ...)` (`Config.pm:761`), so a trailing comma names no
/// criterion and an empty value is the empty order, under which the pick family
/// keeps the first row (`OutputFactory.pm:808`). The list is one clap value (a
/// `Vec` field would be appended to on a repeat), so a repeat keeps the last
/// list, as the scalar `pick_order=s` does (`Config.pm:117`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PickOrder(pub Vec<String>);

impl std::ops::Deref for PickOrder {
    type Target = [String];

    fn deref(&self) -> &[String] {
        &self.0
    }
}

/// The value parser of `--pick_order`: every field must be one of Ensembl VEP's
/// `%VALID` tokens (`Config.pm:601`, checked at `:946`), each checked by clap's
/// own enumeration parser, so an unknown token or a middle empty field is the
/// invalid-value error an enumerated flag gives, its suggestion included.
#[derive(Clone)]
struct PickOrderParser;

impl clap::builder::TypedValueParser for PickOrderParser {
    type Value = PickOrder;

    fn parse_ref(
        &self,
        cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<PickOrder, clap::Error> {
        let Some(text) = value.to_str() else {
            return Err(clap::Error::new(clap::error::ErrorKind::InvalidUtf8).with_cmd(cmd));
        };
        let field_parser = clap::builder::PossibleValuesParser::new(
            crate::option_rules::valid_values("pick_order"),
        );
        let order = crate::option_rules::split_list_flag(text);
        for field in &order {
            field_parser.parse_ref(cmd, arg, std::ffi::OsStr::new(field))?;
        }
        Ok(PickOrder(order))
    }

    fn possible_values(
        &self,
    ) -> Option<Box<dyn Iterator<Item = clap::builder::PossibleValue> + '_>> {
        Some(Box::new(
            crate::option_rules::valid_values("pick_order")
                .iter()
                .map(clap::builder::PossibleValue::new),
        ))
    }
}

impl Args {
    /// The command with the error table's help text applied to every flag the
    /// table governs; the parser `main` runs and the source of `--help` and the
    /// CLI reference.
    pub fn command_with_error_table() -> clap::Command {
        Args::command().mut_args(|arg| {
            let entry = arg
                .get_long()
                .and_then(crate::config::error_table_entry)
                .filter(|e| e.declared);
            match entry {
                Some(entry) => arg.help(entry.help()),
                None => arg,
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{Config, ParseFailure};

    #[test]
    fn test_default_args() {
        let args = Args::parse_from(["vep"]);
        assert_eq!(args.input_file, "STDIN");
        assert_eq!(args.output_file, ["variant_effect_output.txt"]);
        assert_eq!(args.output_format, "vep");
        assert_eq!(args.format, "guess");
        assert_eq!(args.species, "homo_sapiens");
        assert_eq!(args.buffer_size, 5000);
        assert_eq!(args.fork, 0);
        assert!(args.cache.is_none());
        assert!(!args.offline);
        assert!(!args.everything);
    }

    /// `--cache` takes an optional value (`cache:s`, Config.pm:238): alone it is
    /// on, and a following path is its directory.
    #[test]
    fn cache_takes_an_optional_directory() {
        let args = Args::parse_from(["vep", "--cache", "--offline"]);
        assert_eq!(args.cache.as_deref(), Some("1"));
        assert!(args.offline);
        let args = Args::parse_from(["vep", "--cache", "/opt/vep", "--offline"]);
        assert_eq!(args.cache.as_deref(), Some("/opt/vep"));
        let args = Args::parse_from(["vep", "--cache=/opt/vep"]);
        assert_eq!(args.cache.as_deref(), Some("/opt/vep"));
    }

    /// `--synonyms` takes a file path (`synonyms=s`, Config.pm:250), in either
    /// spelling, and is absent by default.
    #[test]
    fn synonyms_flag_parses() {
        let args = Args::parse_from(["vep", "--synonyms", "chr_synonyms.txt"]);
        assert_eq!(args.synonyms.as_deref(), Some("chr_synonyms.txt"));
        let args = Args::parse_from(["vep", "--synonyms=/data/syn.txt"]);
        assert_eq!(args.synonyms.as_deref(), Some("/data/syn.txt"));
        assert!(Args::parse_from(["vep"]).synonyms.is_none());
        assert!(Args::try_parse_from(["vep", "--synonyms"]).is_err());
    }

    #[test]
    fn test_input_output_flags() {
        let args = Args::parse_from([
            "vep",
            "-i",
            "test.vcf",
            "-o",
            "out.txt",
            "--cache",
            "--offline",
        ]);
        assert_eq!(args.input_file, "test.vcf");
        assert_eq!(args.output_file, ["out.txt"]);
        assert_eq!(args.cache.as_deref(), Some("1"));
        assert!(args.offline);
    }

    /// A repeated single-value flag keeps its last value, as Getopt::Long does;
    /// `-o` accumulates (`output_file|o=s@`) so the duplicate can be reported.
    /// `--pick_order` is a scalar option too (`pick_order=s`, Config.pm:117, not
    /// in `@ALLOW_MULTIPLE`, :340-346): its last occurrence is the whole list.
    #[test]
    fn repeated_flags_keep_the_last_value() {
        let args = Args::parse_from(["vep", "-i", "a.vcf", "-i", "b.vcf", "-o", "x", "-o", "y"]);
        assert_eq!(args.input_file, "b.vcf");
        assert_eq!(args.output_file, ["x", "y"]);
        let args = Args::parse_from([
            "vep", "--fork", "2", "--fork", "3", "--sift", "p", "--sift", "b",
        ]);
        assert_eq!(args.fork, 3);
        assert_eq!(args.sift.as_deref(), Some("b"));
        let args = Args::parse_from(["vep", "--pick_order", "length", "--pick_order", "rank"]);
        assert_eq!(args.pick_order.as_deref(), Some(&["rank".to_string()][..]));
        let args = Args::parse_from([
            "vep",
            "--pick_order",
            "length,canonical",
            "--pick_order=rank,tsl",
        ]);
        assert_eq!(
            args.pick_order.as_deref(),
            Some(&["rank".to_string(), "tsl".to_string()][..])
        );
    }

    /// `--pick_order` is read with Perl's `split` (`Config.pm:761`): a trailing
    /// comma names no criterion, an empty value is the empty order (given, and a
    /// list of nothing), and a middle empty field or an unknown token is clap's
    /// invalid-value error, in clap's words for the piece, listing Ensembl VEP's
    /// `%VALID` tokens (`:601`).
    #[test]
    fn pick_order_splits_as_perl_does() {
        let parse = |value: &str| Args::try_parse_from(["vep", "--pick_order", value]);
        assert_eq!(
            parse("canonical,").unwrap().pick_order.as_deref(),
            Some(&["canonical".to_string()][..])
        );
        assert_eq!(parse("").unwrap().pick_order.as_deref(), Some(&[][..]));
        assert_eq!(parse(",,").unwrap().pick_order.as_deref(), Some(&[][..]));
        for (value, message) in [
            (
                "canonical,,rank",
                "a value is required for '--pick_order <PICK_ORDER>' but none was supplied",
            ),
            (
                "canonical,x",
                "invalid value 'x' for '--pick_order <PICK_ORDER>'",
            ),
        ] {
            let error = parse(value).unwrap_err();
            assert_eq!(
                error.kind(),
                clap::error::ErrorKind::InvalidValue,
                "{value}"
            );
            let text = error.to_string();
            assert!(
                text.starts_with(&format!("error: {message}")),
                "{value}: {text}"
            );
            assert!(
                text.contains(
                    "[possible values: mane_select, mane_plus_clinical, canonical, appris, \
                     tsl, biotype, ccds, rank, length, ensembl, refseq]"
                ),
                "{value}: {text}"
            );
        }
        let near_miss = parse("canonica").unwrap_err().to_string();
        assert!(
            near_miss.contains("tip: a similar value exists: 'canonical'"),
            "{near_miss}"
        );
    }

    #[test]
    fn test_force_alias() {
        let args = Args::parse_from(["vep", "--force"]);
        assert!(args.force_overwrite);
    }

    #[test]
    fn test_threads_alias() {
        let args = Args::parse_from(["vep", "--threads", "4"]);
        assert_eq!(args.fork, 4);
    }

    #[test]
    fn test_fa_alias() {
        let args = Args::parse_from(["vep", "--fa", "/path/to/ref.fa"]);
        assert_eq!(args.fasta.as_deref(), Some("/path/to/ref.fa"));
    }

    /// `-u` and `--pass` are Ensembl VEP's short forms of `--user` and `--password`
    /// (`Config.pm:90-91`), so they reach the same refusal.
    #[test]
    fn database_short_forms_parse_as_their_flags() {
        let args = Args::parse_from(["vep", "-u", "anonymous", "--pass", "secret"]);
        assert_eq!(args.user.as_deref(), Some("anonymous"));
        assert_eq!(args.password.as_deref(), Some("secret"));
    }

    #[test]
    fn test_everything_flag() {
        let args = Args::parse_from(["vep", "--everything"]);
        assert!(args.everything);
    }

    #[test]
    fn test_output_format_shortcuts() {
        let args = Args::parse_from(["vep", "--vcf"]);
        assert!(args.vcf);

        let args = Args::parse_from(["vep", "--json"]);
        assert!(args.json);

        let args = Args::parse_from(["vep", "--tab"]);
        assert!(args.tab);
    }

    /// `--plugin` accumulates (`plugin=s@`, Config.pm:271) and the plugins load in
    /// command-line order (Runner.pm:925 iterates the list as given).
    #[test]
    fn repeated_plugin_flags_load_in_order() {
        let args = Args::parse_from([
            "vep",
            "--plugin",
            "CADD,/path/to/cadd.tsv.gz",
            "--plugin",
            "dbNSFP,/path/to/dbnsfp.gz",
        ]);
        assert_eq!(
            args.plugin,
            ["CADD,/path/to/cadd.tsv.gz", "dbNSFP,/path/to/dbnsfp.gz"]
        );
    }

    /// The warning sink's flags parse: `--safe` (Config.pm:272), `--warning_file`
    /// (`:173`), `--skipped_variants_file` (`:174`).
    #[test]
    fn warning_sink_flags_parse() {
        let args = Args::parse_from([
            "vep",
            "--safe",
            "--warning_file",
            "w.txt",
            "--skipped_variants_file",
            "s.txt",
        ]);
        assert!(args.safe);
        assert_eq!(args.warning_file.as_deref(), Some("w.txt"));
        assert_eq!(args.skipped_variants_file.as_deref(), Some("s.txt"));
        let args = Args::parse_from(["vep"]);
        assert!(!args.safe);
        assert!(args.warning_file.is_none());
        assert!(args.skipped_variants_file.is_none());
    }

    #[test]
    fn test_underscore_flag_names() {
        let args = Args::parse_from([
            "vep",
            "--dir_cache",
            "/cache",
            "--dir_plugins",
            "/plugins",
            "--cache_version",
            "115",
            "--buffer_size",
            "1000",
            "--output_format",
            "json",
            "--force_overwrite",
            "--no_headers",
        ]);
        assert_eq!(args.dir_cache.as_deref(), Some("/cache"));
        assert_eq!(args.dir_plugins.as_deref(), Some("/plugins"));
        assert_eq!(args.cache_version, Some(115));
        assert_eq!(args.buffer_size, 1000);
        assert_eq!(args.output_format, "json");
        assert!(args.force_overwrite);
        assert!(args.no_headers);
    }

    #[test]
    fn test_output_control_underscore_flags() {
        let args = Args::parse_from([
            "vep",
            "--mane_select",
            "--xref_refseq",
            "--variant_class",
            "--af_1kg",
            "--af_gnomade",
            "--af_gnomadg",
            "--max_af",
            "--gene_phenotype",
            "--check_existing",
            "--coding_only",
            "--no_intergenic",
        ]);
        assert!(args.mane_select);
        assert!(args.xref_refseq);
        assert!(args.variant_class);
        assert!(args.af_1kg);
        assert!(args.af_gnomade);
        assert!(args.af_gnomadg);
        assert!(args.max_af);
        assert!(args.gene_phenotype);
        assert!(args.check_existing);
        assert!(args.coding_only);
        assert!(args.no_intergenic);
    }

    /// `--failed` is Ensembl VEP's `failed=i` (Config.pm:119, default 0 at :308): an
    /// integer, negative included, and a syntax error for anything else; the two
    /// co-location switches parse beside it.
    #[test]
    fn failed_takes_an_integer() {
        let config = Config::from_argv(["vep", "--check_existing"]).unwrap();
        assert_eq!(config.failed, 0);
        assert!(!config.no_check_alleles);
        assert!(!config.exclude_null_alleles);
        let config = Config::from_argv([
            "vep",
            "--check_existing",
            "--failed",
            "1",
            "--no_check_alleles",
            "--exclude_null_alleles",
        ])
        .unwrap();
        assert_eq!(config.failed, 1);
        assert!(config.no_check_alleles);
        assert!(config.exclude_null_alleles);
        assert_eq!(
            Config::from_argv(["vep", "--failed=-1"]).unwrap().failed,
            -1
        );
        assert_eq!(
            Config::from_argv(["vep", "--failed", "-1"]).unwrap().failed,
            -1
        );
        let err = Config::from_argv(["vep", "--failed", "yes"]).unwrap_err();
        assert!(
            matches!(err, ParseFailure::Usage(_)),
            "a non-integer is clap's own error: {err}"
        );
    }

    /// `--check_svs` names the structural-variant file to read, where Ensembl
    /// VEP's `check_svs` (Config.pm:125) is a switch reading its database; the
    /// flag beside `--offline`, which Ensembl VEP refuses (`%INCOMPATIBLE`,
    /// Config.pm:628), is accepted here, the file needing no database.
    #[test]
    fn check_svs_takes_the_sv_file() {
        assert_eq!(Config::from_argv(["vep"]).unwrap().check_svs, None);
        let config = Config::from_argv(["vep", "--check_svs", "svs.vcf.gz"]).unwrap();
        assert_eq!(config.check_svs.as_deref(), Some("svs.vcf.gz"));
        let offline = Config::from_argv(["vep", "--offline", "--check_svs=svs.vcf.gz"]).unwrap();
        assert_eq!(offline.check_svs.as_deref(), Some("svs.vcf.gz"));
        assert!(offline.offline);
        let err = Config::from_argv(["vep", "--check_svs"]).unwrap_err();
        assert!(
            matches!(err, ParseFailure::Usage(_)),
            "the bare flag is clap's own error: {err}"
        );
    }

    /// `--af_gnomad` is Ensembl VEP's older spelling of `--af_gnomade`
    /// (Config.pm:132), its own flag with the same output fields
    /// (Constants.pm:120-121) and the same option set, which switches
    /// `--check_existing` on (:449-452).
    #[test]
    fn af_gnomad_parses_as_its_own_flag() {
        let config = Config::from_argv(["vep", "--af_gnomad"]).unwrap();
        assert!(config.af_gnomad);
        assert!(!config.af_gnomade);
        assert!(config.check_existing);
        let config = Config::from_argv(["vep", "--af_gnomade"]).unwrap();
        assert!(config.af_gnomade);
        assert!(!config.af_gnomad);
        let plain = Config::from_argv(["vep"]).unwrap();
        assert!(!plain.af_gnomad && !plain.af_gnomade && !plain.check_existing);
    }

    /// The two spellings together are refused as `%INCOMPATIBLE` pairs them
    /// (Config.pm:619), with Ensembl VEP's message.
    #[test]
    fn af_gnomad_and_af_gnomade_together_is_an_error() {
        let err = Config::from_argv(["vep", "--af_gnomade", "--af_gnomad"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [cli] Can't use --af_gnomade and --af_gnomad together"
        );
    }

    /// The frequency filter's flags take Ensembl VEP's defaults (Config.pm:319-322:
    /// `freq_pop 1KG_ALL`, `freq_freq 0.01`, `freq_gt_lt gt`, `freq_filter exclude`),
    /// `--check_frequency` switches `--check_existing` on (:449-452) and
    /// `--filter_common` switches `--check_frequency` on (:442-446). A
    /// `--freq_gt_lt` or `--freq_filter` value outside the two words is clap's
    /// own error: Ensembl VEP takes any word, but its `check_pass`
    /// (BaseCacheVariation.pm:467-468) passes an allele only under `gt` or `lt`,
    /// so any other word passes nothing, and `frequency_check_buffer` (:319)
    /// reads any `--freq_filter` but `include` as `exclude`. `--freq_freq` takes
    /// the real-number shapes Getopt::Long takes (`.5`, `1e-2`, `+0.5`, `-1`,
    /// `1_0` read as 10 with the `_` dropped) and refuses `nan`, `inf`, `1.`
    /// and `_1`, which Ensembl VEP 116.2 refuses as the option's value.
    #[test]
    fn frequency_filter_flags_take_ensembl_veps_defaults() {
        let plain = Config::from_argv(["vep"]).unwrap();
        assert!(!plain.check_frequency && !plain.check_existing);
        assert_eq!(plain.freq_pop, "1KG_ALL");
        assert_eq!(plain.freq_freq, 0.01);
        assert_eq!(plain.freq_gt_lt, "gt");
        assert_eq!(plain.freq_filter, "exclude");
        let on = Config::from_argv(["vep", "--check_frequency"]).unwrap();
        assert!(on.check_frequency && on.check_existing);
        let common = Config::from_argv(["vep", "--filter_common"]).unwrap();
        assert!(common.check_frequency && common.check_existing);
        let given = Config::from_argv([
            "vep",
            "--check_frequency",
            "--freq_pop",
            "gnomADe_AFR",
            "--freq_freq",
            "0.05",
            "--freq_gt_lt",
            "lt",
            "--freq_filter",
            "include",
        ])
        .unwrap();
        assert_eq!(given.freq_pop, "gnomADe_AFR");
        assert_eq!(given.freq_freq, 0.05);
        assert_eq!(given.freq_gt_lt, "lt");
        assert_eq!(given.freq_filter, "include");
        for (text, value) in [
            (".5", 0.5),
            ("1e-2", 0.01),
            ("+0.5", 0.5),
            ("-1", -1.0),
            ("1_0", 10.0),
            ("0.0_1", 0.01),
        ] {
            let config =
                Config::from_argv(["vep", "--check_frequency", "--freq_freq", text]).unwrap();
            assert_eq!(config.freq_freq, value, "{text}");
        }
        for argv in [
            ["vep", "--check_frequency", "--freq_gt_lt", "ge"],
            ["vep", "--check_frequency", "--freq_filter", "Include"],
            ["vep", "--check_frequency", "--freq_freq", "one"],
            ["vep", "--check_frequency", "--freq_freq", "nan"],
            ["vep", "--check_frequency", "--freq_freq", "inf"],
            ["vep", "--check_frequency", "--freq_freq", "1."],
            ["vep", "--check_frequency", "--freq_freq", "5e"],
            ["vep", "--check_frequency", "--freq_freq", "_1"],
        ] {
            let err = Config::from_argv(argv).unwrap_err();
            assert!(
                matches!(err, ParseFailure::Usage(_)),
                "{}: clap's own error: {err}",
                argv.join(" ")
            );
        }
    }

    /// `--max_af` parses as a bare flag and switches `--check_existing` on
    /// (Config.pm:449-452).
    #[test]
    fn max_af_flag_parses() {
        let config = Config::from_argv(["vep", "--max_af"]).unwrap();
        assert!(config.max_af);
        assert!(config.check_existing);
        assert!(!Config::from_argv(["vep"]).unwrap().max_af);
        let err = Config::from_argv(["vep", "--max_af", "1"]).unwrap_err();
        assert!(
            matches!(err, ParseFailure::Usage(_)),
            "a value is clap's own error: {err}"
        );
    }

    #[test]
    fn protein_version_and_mask_header_cache_path_parse() {
        let config =
            Config::from_argv(["vep", "--protein_version", "--mask_header_cache_path"]).unwrap();
        assert!(config.protein_version);
        assert!(config.mask_header_cache_path);
    }

    /// `--transcript_version` and `--gene_version` parse as the version flags
    /// (Config.pm:191-192 at release/116.2), and `--mane` sets
    /// `mane_plus_clinical` as its option set does (Config.pm:580-586) while
    /// `--mane_select` alone leaves it off.
    #[test]
    fn version_flags_and_mane_option_set_parse() {
        let versioned =
            Config::from_argv(["vep", "--transcript_version", "--gene_version"]).unwrap();
        assert!(versioned.transcript_version && versioned.gene_version);
        let plain = Config::from_argv(["vep"]).unwrap();
        assert!(!plain.transcript_version && !plain.gene_version && !plain.mane_plus_clinical);
        let mane = Config::from_argv(["vep", "--mane"]).unwrap();
        assert!(mane.mane && mane.mane_select && mane.mane_plus_clinical);
        let select = Config::from_argv(["vep", "--mane_select"]).unwrap();
        assert!(select.mane_select && !select.mane && !select.mane_plus_clinical);
    }

    /// The help of the attribute flags whose column Ensembl VEP fills from a
    /// source some caches lack states the exclusion the CLI reference prints:
    /// `--tsl` the GRCh37 gene set, `--uniprot`, `--xref_refseq` and `--mirna`
    /// a native cache built without `--core-dumps`, `--ccds` the RefSeq and
    /// merged cache types no native build produces.
    #[test]
    fn attribute_flags_state_their_documented_exclusions() {
        let mut command = Args::command_with_error_table();
        command.build();
        let help = |flag: &str| {
            command
                .get_arguments()
                .find(|a| a.get_long() == Some(flag))
                .and_then(|a| a.get_long_help().or(a.get_help()))
                .map(|h| h.to_string())
                .unwrap()
        };
        assert!(help("tsl").contains("GRCh37 gene set carries no level"));
        for flag in ["uniprot", "xref_refseq", "mirna"] {
            let text = help(flag);
            assert!(
                text.contains("only when built with --core-dumps"),
                "{flag}: {text}"
            );
        }
        assert!(help("ccds").contains("no RefSeq or merged cache is built natively"));
    }

    /// `--sift` takes Ensembl VEP's `%VALID` values `s`, `p` and `b`
    /// (Config.pm:599) and nothing else: a value outside them is refused when the
    /// command line is read, the valid values listed, as every `%VALID` flag is;
    /// the option-rule check behind the parser carries Ensembl VEP's own text for
    /// the same value (Config.pm:945-957).
    #[test]
    fn sift_value_outside_spb_is_an_error() {
        for value in ["s", "p", "b"] {
            assert_eq!(
                Config::from_argv(["vep", "--sift", value])
                    .unwrap()
                    .sift
                    .as_deref(),
                Some(value)
            );
        }
        let ParseFailure::Usage(error) = Config::from_argv(["vep", "--sift", "x"]).unwrap_err()
        else {
            panic!("--sift x is refused by the parser");
        };
        assert_eq!(error.kind(), clap::error::ErrorKind::InvalidValue);
        let text = error.to_string();
        assert!(text.contains("'x'") && text.contains("s, p, b"), "{text}");
        assert!(matches!(
            Config::from_argv(["vep", "--sift", "both"]).unwrap_err(),
            ParseFailure::Usage(_)
        ));
        let mut state = crate::option_rules::FlagState::default();
        state.set("sift", "x");
        assert_eq!(
            state.check(false),
            Err("\"x\" is not a valid value for --sift\nValid values: s, p, b".to_string())
        );
    }

    /// `--polyphen` takes the same three values (Config.pm:600); `--humdiv`
    /// parses as a flag and rewrites `polyphen_analysis` to `humdiv` through its
    /// option set (Config.pm:517-522), nothing more.
    #[test]
    fn polyphen_value_outside_spb_is_an_error() {
        for value in ["s", "p", "b"] {
            assert_eq!(
                Config::from_argv(["vep", "--polyphen", value])
                    .unwrap()
                    .polyphen
                    .as_deref(),
                Some(value)
            );
        }
        let ParseFailure::Usage(error) =
            Config::from_argv(["vep", "--polyphen", "humdiv"]).unwrap_err()
        else {
            panic!("--polyphen humdiv is refused by the parser");
        };
        assert_eq!(error.kind(), clap::error::ErrorKind::InvalidValue);
        let mut state = crate::option_rules::FlagState::default();
        state.set("polyphen", "P");
        assert_eq!(
            state.check(false),
            Err("\"P\" is not a valid value for --polyphen\nValid values: s, p, b".to_string())
        );
        let humdiv = Config::from_argv(["vep", "--polyphen", "b", "--humdiv"]).unwrap();
        assert!(humdiv.humdiv);
        assert_eq!(humdiv.polyphen.as_deref(), Some("b"));
        let alone = Config::from_argv(["vep", "--humdiv"]).unwrap();
        assert!(alone.humdiv && alone.polyphen.is_none());
        assert!(
            !Config::from_argv(["vep", "--polyphen", "b"])
                .unwrap()
                .humdiv
        );
    }

    /// `--regulatory_gff` is refused: regulatory annotation is not supported.
    #[test]
    fn regulatory_gff_is_refused() {
        let err = Config::from_argv(["vep", "--regulatory_gff", "reg.gff3.gz"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [cli] --regulatory_gff is not supported: it would annotate regulatory features \
             and motifs from a GFF3 file in place of the cache, and regulatory annotation \
             is not supported"
        );
    }

    /// `--extended_promoters` is refused: regulatory annotation is not supported.
    #[test]
    fn extended_promoters_is_refused() {
        let err = Config::from_argv(["vep", "--extended_promoters"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [cli] --extended_promoters is not supported: it would widen the promoters read \
             with --regulatory_gff to their extended bounds, and regulatory annotation is \
             not supported"
        );
    }

    /// `--custom_suppress_filter` is refused: custom annotation is not supported.
    #[test]
    fn custom_suppress_filter_is_refused() {
        let err = Config::from_argv(["vep", "--custom_suppress_filter"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [cli] --custom_suppress_filter is not supported: it would keep a --custom VCF \
             source's FILTER column out of the output unless the --custom source's fields \
             argument names it, and custom annotation is not supported"
        );
    }

    /// `--transcript_filter` repeats, each value a whole filter string with its
    /// spaces, as Ensembl VEP's `transcript_filter=s@` takes it.
    #[test]
    fn transcript_filter_repeats() {
        let args = Args::parse_from([
            "vep",
            "--transcript_filter",
            "biotype is protein_coding",
            "--transcript_filter",
            "not stable_id match ENST0000040",
        ]);
        assert_eq!(
            args.transcript_filter,
            [
                "biotype is protein_coding",
                "not stable_id match ENST0000040"
            ]
        );
        assert!(Args::parse_from(["vep"]).transcript_filter.is_empty());
    }

    #[test]
    fn test_filtering_underscore_flags() {
        let args = Args::parse_from([
            "vep",
            "--pick_allele",
            "--per_gene",
            "--pick_allele_gene",
            "--most_severe",
            "--flag_pick",
            "--flag_pick_allele",
            "--flag_pick_allele_gene",
            "--pick_order",
            "canonical,biotype,rank",
            "--allele_number",
            "--show_ref_allele",
            "--total_length",
            "--dont_skip",
            "--allow_non_variant",
        ]);
        assert!(args.pick_allele);
        assert!(args.per_gene);
        assert!(args.pick_allele_gene);
        assert!(args.most_severe);
        assert!(args.flag_pick);
        assert!(args.flag_pick_allele);
        assert!(args.flag_pick_allele_gene);
        assert_eq!(
            args.pick_order.as_deref(),
            Some(
                &[
                    "canonical".to_string(),
                    "biotype".to_string(),
                    "rank".to_string()
                ][..]
            )
        );
        assert!(args.allele_number);
        assert!(args.show_ref_allele);
        assert!(args.total_length);
        assert!(args.dont_skip);
        assert!(args.allow_non_variant);
    }

    /// `shift_hgvs=i` (Config.pm:201) takes any integer in Ensembl VEP, which
    /// BaseRunner.pm:498 matches against `/(0|1)/` as a substring, so `2` runs
    /// with the default and `10` sets the shifting variables to 10. Here `0` and
    /// `1` alone are read: the bare flag and any other value are syntax errors,
    /// so neither passes as a silent default, the departure the CLI reference's
    /// compatibility page states.
    #[test]
    fn shift_hgvs_takes_0_or_1() {
        assert_eq!(
            Args::parse_from(["vep", "--shift_hgvs", "0"])
                .shift_hgvs
                .as_deref(),
            Some("0")
        );
        assert_eq!(
            Args::parse_from(["vep", "--shift_hgvs=1"])
                .shift_hgvs
                .as_deref(),
            Some("1")
        );
        assert_eq!(Args::parse_from(["vep"]).shift_hgvs, None);
        for argv in [
            &["vep", "--shift_hgvs"][..],
            &["vep", "--shift_hgvs", "2"][..],
            &["vep", "--shift_hgvs", "yes"][..],
        ] {
            let err = Args::try_parse_from(argv).unwrap_err();
            assert_eq!(err.exit_code(), 2, "{argv:?}: {err}");
        }
    }

    /// `shift_3prime=i` and `shift_genomic=i` (Config.pm:231-232) take a value as
    /// `--shift_hgvs` does; the bare `--shift_3prime` is a syntax error, and
    /// `--shift_length` is a plain switch (Config.pm:233).
    #[test]
    fn shift_3prime_takes_0_or_1() {
        let args = Args::parse_from(["vep", "--shift_3prime", "1", "--shift_genomic", "0"]);
        assert_eq!(args.shift_3prime.as_deref(), Some("1"));
        assert_eq!(args.shift_genomic.as_deref(), Some("0"));
        assert!(!args.shift_length);
        assert_eq!(
            Args::parse_from(["vep", "--shift_3prime=0"])
                .shift_3prime
                .as_deref(),
            Some("0")
        );
        assert!(Args::parse_from(["vep", "--shift_length"]).shift_length);
        for argv in [
            &["vep", "--shift_3prime"][..],
            &["vep", "--shift_3prime", "2"][..],
            &["vep", "--shift_genomic"][..],
            &["vep", "--shift_genomic", "-1"][..],
        ] {
            let err = Args::try_parse_from(argv).unwrap_err();
            assert_eq!(err.exit_code(), 2, "{argv:?}: {err}");
        }
    }

    /// `nearest=s` takes one of `transcript`, `gene` and `symbol` (Config.pm:149,
    /// 602); any other value, and the bare flag, is a syntax error.
    #[test]
    fn nearest_rejects_other_values() {
        for value in ["transcript", "gene", "symbol"] {
            assert_eq!(
                Args::parse_from(["vep", "--nearest", value])
                    .nearest
                    .as_deref(),
                Some(value)
            );
        }
        assert!(Args::parse_from(["vep"]).nearest.is_none());
        for argv in [
            &["vep", "--nearest"][..],
            &["vep", "--nearest", "Symbol"][..],
            &["vep", "--nearest", "exon"][..],
            &["vep", "--nearest=1"][..],
        ] {
            let err = Args::try_parse_from(argv).unwrap_err();
            assert_eq!(err.exit_code(), 2, "{argv:?}: {err}");
        }
    }

    /// Every flag has help text: its own doc comment, or the error table's.
    #[test]
    fn every_flag_has_help() {
        let mut command = Args::command_with_error_table();
        command.build();
        let without: Vec<&str> = command
            .get_arguments()
            .filter(|a| a.get_help().is_none() && a.get_long_help().is_none())
            .filter_map(|a| a.get_long())
            .collect();
        assert!(without.is_empty(), "flags without help: {without:?}");
    }

    /// Every flag is documented under a help heading, one of the reference's sections.
    #[test]
    fn every_flag_has_a_help_heading() {
        let mut command = Args::command();
        command.build();
        let without: Vec<&str> = command
            .get_arguments()
            .filter(|a| a.get_long().is_some_and(|l| l != "help" && l != "version"))
            .filter(|a| a.get_help_heading().is_none())
            .filter_map(|a| a.get_long())
            .collect();
        assert!(
            without.is_empty(),
            "flags without a help heading: {without:?}"
        );
    }
}
