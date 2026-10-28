// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of `tests/parity/matrix.tsv` whose exclusion a public page states:
//! one test per row, or per family of rows one page states together, each
//! asserting the page's own sentences, so a rewrite that drops the exclusion
//! fails here.

mod common;

/// The text of the page at `path` below `docs/`.
fn page_text(path: &str) -> String {
    std::fs::read_to_string(format!("{}/../../docs/{path}", env!("CARGO_MANIFEST_DIR")))
        .unwrap_or_else(|e| panic!("docs/{path}: {e}"))
}

fn assert_states(what: &str, text: &str, phrases: &[&str]) {
    for phrase in phrases {
        assert!(text.contains(phrase), "{what} lacks {phrase:?}");
    }
}

/// Asserts that the page at `path` below `docs/` states every phrase, its
/// whitespace collapsed so that a sentence wrapped over several lines is
/// matched as one.
fn assert_page_states(path: &str, phrases: &[&str]) {
    let text = page_text(path)
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    assert_states(&format!("docs/{path}"), &text, phrases);
}

/// Asserts that the flag table row of a generated CLI page that opens with
/// `| \`--flag\`` states every phrase.
fn assert_flag_row_states(path: &str, flag: &str, phrases: &[&str]) {
    let text = page_text(path);
    let prefix = format!("| `{flag}`");
    let row = text
        .lines()
        .find(|l| l.starts_with(&prefix))
        .unwrap_or_else(|| panic!("docs/{path} has no {flag} row"));
    assert_states(&format!("the {flag} row of docs/{path}"), row, phrases);
}

/// `check_ref`: a contig the FASTA lacks skips the record at every position with
/// the sub-slice and mismatch warnings, where BaseVEP.pm:424-434 gives such a
/// contig a one-base placeholder slice reading `N` (FastaSequence.pm:597-601), so
/// Parser.pm:614 passes a one-base record at position 1 written `N` and :624-631
/// refuses any other with the mismatch line alone; an insertion names no
/// reference base, so Parser.pm:600-602 (`--check_ref`) and :636 (`--lookup_ref`)
/// pass it, as vep-rs does.
#[test]
fn check_ref_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Under `--check_ref` the same contig skips the record at every position with the \
             sub-slice warning and then the mismatch warning; an insertion, naming no reference \
             base, passes under either flag in both engines (`Parser.pm` lines 600-602 and 636)",
            "Ensembl VEP annotates a one-base record there whose reference is written `N`, and \
             refuses one written otherwise with the mismatch warning alone (`Parser.pm` lines \
             614-631)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &["under the first two flags, an insertion excepted in both engines"],
    );
}

/// `id-input-format`: an identifier list is resolved from the cache's variation
/// files, where Parser/ID.pm:100 refuses the format offline, and a list line `0`
/// is a name warned `No variant found`, where ID.pm:146 reads Perl's false string
/// as no record and the buffer fill ends (Parser.pm:186-187, InputBuffer.pm:212,
/// Runner.pm:423-424).
#[test]
fn id_input_exclusion_is_documented() {
    assert_page_states(
        "divergences/id-input.md",
        &[
            "Ensembl VEP reads an identifier list (`--format id`) only with a database \
             connection",
            "vep-rs resolves the list from the cache's variation files",
            "A structural variant name (`nsv`, `esv`), a synonym, a subsnp identifier and a \
             variant without a genomic mapping all take the `No variant found` warning",
            "A list line `0` is a name here and takes the `No variant found` warning",
            "Ensembl VEP reads that line, Perl's false string, as no record (`Parser/ID.pm` line \
             146) and ends the buffer fill there",
            "a `0` after a name is dropped without a word, and one that opens a fill, the list's \
             first line among them, ends the run with the names after it unread",
        ],
    );
}

/// `default-format-header-lines`: the `## Output produced at` time is UTC, where
/// Utils.pm:467-468 reads `localtime`, and the Perl API, database and `--custom`
/// header lines have no counterpart.
#[test]
fn default_format_header_lines_exclusion_is_documented() {
    assert_page_states(
        "output/vep-and-tab.md",
        &[
            "the time, in UTC where Ensembl VEP prints the process's local time (`Utils.pm` \
             lines 467-468)",
            "Two kinds of line Ensembl VEP prints there have no counterpart",
            "The Perl API component versions (`## ensembl version 116.0d85231`, `## \
             ensembl-variation version ...`) are not written, since no Perl API runs here",
            "The `## Connected to <database>` line of a database run is not written",
            "The `## <short name> : [PATH]/<file>` lines of `--custom` sources are not written \
             either",
        ],
    );
}

/// `vcf-csq-info-header`: the `##VEP=` line carries no Perl API component
/// versions and no `--custom` `##INFO` line is written (OutputFactory/VCF.pm:234-263,
/// 273), and a header without a `#CHROM` line gets the eight-column line, where
/// VCF.pm:224 pops the last `#` line of the input whatever it is.
#[test]
fn vcf_csq_info_header_exclusion_is_documented() {
    assert_page_states(
        "output/vcf.md",
        &[
            "Ensembl VEP's `##VEP=` line also carries its Perl API component versions \
             (`ensembl=116.0d85231 ...`), which have no counterpart here",
            "Its `##INFO=<ID=<short>_<field>,...>` lines for `--custom` sources are not \
             written, since custom annotation sources are not supported",
        ],
    );
    assert_page_states(
        "divergences/vcf-without-chrom-line.md",
        &[
            "Ensembl VEP's VCF output copies the input's `#` lines into its header and takes \
             the last of them as the column line (`OutputFactory/VCF.pm` line 224)",
            "A header without a `#CHROM` line hands its last meta line to that place",
            "vep-rs takes the input's `#CHROM` line as the column line when the header has one \
             and writes the eight-column line `#CHROM POS ID REF ALT QUAL FILTER INFO` \
             otherwise; every meta line stays in the meta block",
        ],
    );
}

/// `single-dash-long-options`: a multi-letter one-dash name reads as Getopt::Long
/// reads it under Config.pm:62, and the tokens that differ (`-qv`, `-ofile.txt`,
/// an upper-case single-letter short, `+name`) are clap's readings.
#[test]
fn single_dash_long_options_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A short with its value attached (`-ofile.txt`) and a cluster of switches (`-qv`) \
             are clap's readings, where Getopt::Long refuses the token",
            "A `+` before an option name, which Getopt::Long accepts, is not an option here",
            "A single-letter short is clap's and is matched by case",
            "Ensembl VEP refuses `-h`, a prefix of eight of its options, and has no `--version`",
        ],
    );
}

/// `stdin-input`, `gzip-bgzip-compressed-input`, `gzip-bgzip-input`: gzip or
/// BGZF bytes on standard input or a pipe are decompressed, where Parser.pm:242
/// decompresses a regular file alone, and a stream cut inside a member ends the
/// run, where Ensembl VEP annotates the records before the cut and exits 0.
#[test]
fn compressed_input_exclusion_is_documented() {
    assert_page_states(
        "divergences/compressed-input.md",
        &[
            "Ensembl VEP decompresses an input only when it is a regular file whose bytes \
             Perl's `-B` test calls binary (`Parser.pm` line 242)",
            "vep-rs decompresses gzip or BGZF bytes on standard input, and on a pipe named by \
             path, as it decompresses a file, so the piped stream annotates exactly as the \
             file does",
            "`--format` is still required on standard input, as in Ensembl VEP",
            "vep-rs ends the run at the cut: one line `ERROR: [run] Failed to read input ...` on \
             standard error, exit 1, and nothing past the header written",
        ],
    );
}

/// `attributes-cds-start-end-nf`: a native cache lists `cds_start_NF` before
/// `cds_end_NF` on every transcript carrying both, where the Ensembl cache keeps
/// its database's order (Dumper/Core.pm:203-207) and OutputFactory.pm:1477 prints
/// it, `cds_end_NF` first on a few GRCh38 transcripts.
#[test]
fn cds_flag_order_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "A native cache lists `cds_start_NF` before `cds_end_NF` on every transcript \
             carrying both, whatever order the core dumps or the GTF carried",
            "The Ensembl cache keeps the order its database returned, `cds_end_NF` first on 7 \
             of the 2,946 GRCh38 transcripts that carry both at release 116",
            "Ensembl VEP prints that order",
        ],
    );
}

/// `buffer_size`: `--buffer_size 0` annotates every record where Ensembl VEP's
/// buffer admits nothing (InputBuffer.pm:191, 212), and a negative or underscored
/// value Getopt::Long's `=i` reads (Config.pm:118) is an argument syntax error.
#[test]
fn buffer_size_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "`--buffer_size 0` annotates every record, where Ensembl VEP's buffer admits a \
             record only while it holds fewer than the value",
            "so there a zero writes the header alone",
            "A value with a minus sign or an underscore (`-1`, `1_0`) is an argument syntax \
             error (exit 2)",
        ],
    );
}

/// `cache-fasta-auto-detection`: an unindexed `.fa` or `.fa.gz` in the cache
/// directory is passed over, where CacheDir.pm:235-251 takes the first such
/// file and Ensembl VEP builds the missing index.
#[test]
fn cache_fasta_auto_detection_exclusion_is_documented() {
    assert_flag_row_states(
        "cli/cache-species-runtime.md",
        "--fasta",
        &[
            "Without the flag, the first indexed .fa or .fa.gz by name in the cache directory \
             is used",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &["An unindexed file there is not read (Ensembl VEP builds a missing index itself)"],
    );
}

/// `af-1kg-columns`: a native cache writes no 1000 Genomes entry for an allele
/// at 0 in every population and derives `AF` from the counts, where
/// Dumper/Variation.pm:370-378 copies the release's `AF` and writes `X:0`.
#[test]
fn af_1kg_columns_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "The file states no panel-wide frequency, so `AF` is the allele's count over the \
             panel's 5,008 chromosomes",
            "An allele whose frequency is 0 in every population gets no entry",
            "For the latter the Ensembl cache prints `0`",
            "The release's own panel-wide frequency, which the Ensembl cache prints, differs \
             from the recovered count's for 443 of those records",
        ],
    );
}

/// `allele_string`: a native cache holds the export's right-trimmed indel
/// alleles (ensembl-variation VariationFeature.pm:2390), where the Ensembl cache
/// holds the database's string over the whole run (Database/Variation.pm:193).
#[test]
fn allele_string_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "The Ensembl cache keeps dbSNP's alleles over the whole run an indel can shift \
             along (`CCCCC/CCCC` at 5030241-5030245), which the export trims to `CC/C` at \
             5030241",
            "A native cache holds the trimmed form",
            "An input anchored on the run's first base (`CC/C` at 5030241) minimises to the \
             record's second base",
            "so the input co-locates with the native record alone",
        ],
    );
}

/// `byte_transparent_input`: a byte that is not UTF-8 is carried into the VCF
/// output and printed as U+FFFD in the default, tab and JSON outputs, where
/// Ensembl VEP writes the byte in every format (Parser.pm:235-242).
#[test]
fn byte_transparent_input_exclusion_is_documented() {
    assert_page_states(
        "divergences/non-utf8-input.md",
        &[
            "Ensembl VEP opens its input with no decoding layer",
            "comes out as that byte in every format",
            "So a VCF written from a VCF carries every input byte",
            "A record is otherwise read as UTF-8 with each invalid sequence replaced by \
             U+FFFD, sample columns included",
            "prints the replacement character in `Uploaded_variation`, in the JSON `input` \
             string",
        ],
    );
}

/// `cache`: a run naming none of `--database`, `--cache`, `--offline` and
/// `--custom` proceeds on the cache it is given, where Config.pm:993-1017 dies.
#[test]
fn cache_flag_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A run naming none of `--database`, `--cache`, `--offline` and `--custom` is not \
             refused (`Config.pm:993-1017`)",
            "vep-rs reads a local cache in every run, so those flags select nothing and the \
             run proceeds on the cache it is given",
        ],
    );
}

/// The page's sentences on a configuration-file key alone on its line: a scalar
/// or switch key alone is dropped by both engines (Config.pm:1100 `next unless
/// is_valid_param($config, $key, $split[0]);`, whose Getopt::Long call sets no
/// option for the scalar's undefined value, and :1101 `$config->{$key} ||=
/// $split[0];`, which leaves the switch unset), and a repeatable key alone is
/// dropped here (`read_config_file` in `argv.rs`) where Config.pm:1094-1097
/// pushes the empty string, which Runner.pm:935 warns `Failed to compile plugin`
/// for and AnnotationSourceAdaptor.pm:376 dies on.
const KEY_ALONE_ON_ITS_LINE: &[&str] = &[
    "A key alone on its line is dropped",
    "Ensembl VEP drops a scalar or switch key alone on its line too (`Config.pm:1101`) and \
     keeps an empty value for a repeatable key alone on its line (`Config.pm:1094-1097`): \
     `plugin` warns `Failed to compile plugin` there, and `custom` and `transcript_filter` \
     are fatal",
];

/// `config`: a lower source never overwrites a higher source's value, and a
/// value the parser refuses is an argument syntax error, where Config.pm:1101's
/// `||=` lets a false value yield and :1046-1047 drops the line; a repeatable
/// key alone on its line is dropped where Config.pm:1094-1097 keeps `''`.
#[test]
fn config_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A value from the `--config` file, `vep.ini` or the environment never overwrites \
             one a higher source gives, even a `0`, an empty value or a bare `--cache`",
            "Ensembl VEP's `||=` (`Config.pm:1101, 1158`) lets a lower source overwrite a \
             value Perl reads as false",
            "A value from a file or the environment that the parser refuses (`fork many`, an \
             empty `format`) is an argument syntax error (exit 2)",
            "Ensembl VEP's `is_valid_param` drops the line (`Config.pm:1046-1047`)",
        ],
    );
    assert_page_states("cli/compatibility.md", KEY_ALONE_ON_ITS_LINE);
}

/// `config-file`: the same merge rule read from `vep.ini` and the environment
/// (Config.pm:1101, 1158), with the page's own `shift_hgvs` example, and the
/// same reading of a repeatable key alone on its line.
#[test]
fn config_file_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "so `--shift_hgvs 0` on its command line yields to a file's `shift_hgvs 1`",
            "A flag read from a file or the environment acts exactly as if typed, refusals \
             included",
        ],
    );
    assert_page_states("cli/compatibility.md", KEY_ALONE_ON_ITS_LINE);
}

/// `dont_skip`: four refused shapes are skipped under the flag, where
/// Parser.pm:189 returns a record whatever `validate_vf` returned; and a `#`
/// line of the list-based formats counts toward the line numbers, as Ensembl
/// VEP numbers the record it reads there (Parser/Region.pm:167), a region or
/// HGVS `#` line with no other word character an empty line (Parser.pm:273-275).
#[test]
fn dont_skip_exclusion_is_documented() {
    assert_page_states(
        "divergences/dont-skip-cannot-keep.md",
        &[
            "`--dont_skip` returns a record whatever `validate_vf` returned for it \
             (`Parser.pm` line 189)",
            "Four refused shapes have no record vep-rs can hold, and vep-rs skips them under \
             the flag too",
            "a coordinate that is not a number",
            "a start past the end plus one",
            "a line with too few columns to name a variant, or one the parser refuses for its \
             shape",
            "It reads an Ensembl-format, region or HGVS line opening with `#` as a comment, \
             counted toward the line numbers that follow, as Ensembl VEP counts the record it \
             reads there, and otherwise skipped without a word",
            "A region or HGVS `#` line with no other word character is an empty line to both \
             engines (`Parser.pm` lines 273-275), counted and warned about as one",
            "a region, SPDI or HGVS transcript record whose span exceeds `--max_sv_size`",
        ],
    );
}

/// `dir_plugins`: the value is one literal path, where Runner.pm:908 passes it
/// through `glob()`.
#[test]
fn dir_plugins_exclusion_is_documented() {
    assert_flag_row_states(
        "cli/plugins.md",
        "--dir_plugins",
        &[
            "A value holding a glob pattern (`*`, `?`, `[...]`, `{a,b}`) or a `~user` prefix is \
             read as one literal path; Ensembl VEP expands it with `glob`",
        ],
    );
}

/// `flag_pick`, `flag_pick_allele`, `flag_pick_allele_gene`, `spdi`: under
/// `--minimal` the pick family runs once over the rejoined record, where Ensembl
/// VEP's rejoin filters each split feature first (OutputFactory.pm:2289, 2303
/// reaching :514 and :609-626) and the marks persist into the merged pass (:301),
/// so a split record carries one PICK per split feature; and a split allele's
/// SPDI and GA4GH VRS are its own, where :1398 builds the table from the rejoined
/// allele string and :1400 looks it up by the trimmed Allele.
#[test]
fn flag_pick_minimal_exclusion_is_documented() {
    assert_page_states(
        "divergences/minimal-split-hgvsp.md",
        &[
            "Each row's notation is then looked up by the row's `Allele`, the split allele's \
             trimmed string (lines 1347 and 1400)",
            "vep-rs builds each split allele's notation from the allele itself",
            "| `GT`, trimmed `-` | `21:38218615:TT:`, the `G` allele's | `21:38218616:T:` |",
            "So do the `SPDI` and `GA4GH_VRS` of a split allele and the `PICK` marks of the \
             flagging variants",
            "A split feature holds one alternate allele, so the flagging branches (lines \
             609-626) set `PICK` on one row of it",
            "So Ensembl VEP marks one row per split feature",
            "Under `--flag_pick` a record `--minimal` splits carries one mark per split \
             allele, where a record it does not split carries one",
            "such a record carries two marks per allele string, or two per allele string and \
             gene, in place of one",
            "vep-rs runs the pick family once over the rejoined record, grouping the rows by \
             the Allele column, so it marks one row per record under `--flag_pick` and one \
             per allele string, or per allele string and gene, under the other two",
            "| `min_gt_tg`, `--flag_pick_allele` | Ensembl VEP | `T` ENST00000291552 at 1 \
             and at 2 |",
            "| `min_gt_tg`, `--flag_pick_allele` | vep-rs | `T` ENST00000291552 at 2 |",
        ],
    );
}

/// `no_intergenic`: under `--minimal --no_intergenic` the genic allele of a
/// multi-allelic record keeps its transcript row, where the rejoin leaves the
/// rejoined feature's intergenic overlap defined (OutputFactory.pm:2323-2338)
/// and :482 drops every overlap of the record.
#[test]
fn no_intergenic_minimal_exclusion_is_documented() {
    assert_page_states(
        "divergences/minimal-no-intergenic.md",
        &[
            "`--no_intergenic` then drops every overlap of a feature that carries an intergenic \
             overlap (line 482), so the rejoined record writes nothing, the genic allele's \
             transcript row included",
            "vep-rs judges each allele alone",
            "removes the intergenic row of the allele that overlaps no transcript and keeps the \
             other allele's transcript row, its CSQ entry and its JSON `transcript_consequences`",
            "Without `--minimal` the two engines agree, since a record's alleles then share one \
             span and one transcript set",
        ],
    );
}

/// `hgvsg`: a row prints its own allele's notation, where Ensembl VEP indexes
/// the record's allele string by an allele number counted without repeated and
/// reference-equal ALTs (OutputFactory.pm:1386-1393, VariationFeatureOverlap.pm:463,
/// 481), so a `C` row of `A A,C` prints the reference allele's `dup`.
#[test]
fn hgvsg_exclusion_is_documented() {
    assert_page_states(
        "divergences/hgvsg-allele-lookup.md",
        &[
            "A row prints the notation stored under one entry of that string: the second \
             entry when the string has two, else the entry at the row's allele number \
             (`OutputFactory.pm` lines 1386-1393)",
            "The allele number is counted over a list from which a repeated allele \
             (`VariationFeatureOverlap.pm` line 463) and an allele equal to the reference \
             (line 481) are dropped first",
            "the row prints another allele's notation",
            "vep-rs computes each row's notation from the row's own allele at its trimmed \
             span, so a repeated or reference-equal entry changes no other row's value",
            "| `refeq`, `C` | 1 | `21:g.34775850dup` | `21:g.34775850A>C` |",
            "| `dupalt`, `G` | 2 | `21:g.34775850A>C` | `21:g.34775850A>G` |",
        ],
    );
}

/// `transcript_filter` and `transcript-filter-grammar`: an abbreviated field
/// reads every transcript's value, the first included, where FilterSet.pm:567-572
/// records the synonym and reads nothing on the first call; `e` names `end` on
/// a transcript without a symbol, where the undefined `external_*` keys
/// Dumper/Core.pm:172 keeps make it ambiguous (FilterSet.pm:574-580); a pattern
/// takes the `regex` crate's dialect, where Perl's engine runs `\Z` (:714-715);
/// `is_child` walks a Sequence Ontology table where `vep` has no ontology
/// adaptor and dies (:736, AnnotationSource.pm:87-90).
#[test]
fn transcript_filter_exclusion_is_documented() {
    assert_page_states(
        "divergences/transcript-filter-fields.md",
        &[
            "The prefix branch records the synonym and reads nothing, so the first \
             transcript a run evaluates against an abbreviated field has no value",
            "vep-rs resolves an abbreviation once, when the filter is parsed, over a fixed \
             list of the cache's transcript keys, and reads every transcript's value, the \
             first included",
            "A transcript without a gene symbol therefore keeps all five, undefined",
            "vep-rs's key list leaves the five out, so `e` names `end` on every transcript",
            "vep-rs compiles it with the `regex` crate, which lacks back-references, \
             look-around and the `\\Z` and `\\K` escapes",
            "Under `--transcript_filter` the operator therefore passes a transcript through \
             the regex alone",
            "vep-rs matches the regex first and otherwise walks its built-in Sequence \
             Ontology table of consequence terms",
        ],
    );
}

/// `shift_3prime`: in Ensembl VEP the alleles of one record share one allele's
/// span predicates (BaseVariationFeatureOverlapAllele.pm:363) and cached
/// coordinates (TranscriptVariationAllele.pm:502-528, BaseTranscriptVariation.pm:142),
/// so a multi-allelic indel's later allele and the reference row of
/// `--process_ref_homs` take another allele's span; vep-rs reads each allele at
/// its own shifted span. The RefSeq edited-sequence exclusion is
/// `hgvs.rs::documented_subset_names_the_refseq_edited_transcript_paths`'s.
#[test]
fn shift_3prime_shared_predicates_exclusion_is_documented() {
    assert_page_states(
        "divergences/shift-3prime-shared-predicates.md",
        &[
            "The alleles of one record on one transcript share one transcript-variation \
             object, and two of its caches are filled once for all of them",
            "The first allele evaluated stores them at its shifted span \
             (`BaseVariationFeatureOverlapAllele.pm` line 363, with the offset of line 460)",
            "so a row reads coordinates another allele of the record left there",
            "Under `--process_ref_homs` the reference allele comes first \
             (`BaseVariationFeatureOverlap.pm` line 370) and is never shifted itself \
             (`TranscriptVariationAllele.pm` line 116), so it stores the predicates at the \
             unshifted span",
            "vep-rs reads each allele's predicates and coordinates at its own shifted span, \
             so a row depends neither on the record's other alleles nor on \
             `--process_ref_homs`",
            "Its positions, codons and coding terms therefore agree with Ensembl VEP's, and \
             its non-coding terms are the shifted span's",
        ],
    );
    assert_page_states(
        "output/vep-and-tab.md",
        &[
            "Its non-coding terms are that span's as well, where Ensembl VEP keeps the \
             unshifted span's and rewrites the alternates' rows",
            "A record whose alleles shift by different lengths shares one allele's span \
             predicates in Ensembl VEP too",
        ],
    );
}

/// `validate-vf-checks`: a POS that is not a number is refused in every form,
/// where Ensembl VEP's anchor chop of a sequence indel runs `$start++` on the
/// text before `validate_vf` (Parser/VCF.pm:310, 334; Parser.pm:475): `12abc`
/// numifies to 13 and the record is annotated, `abc` increments to `abd` and
/// is refused quoting that text.
#[test]
fn validate_vf_checks_exclusion_is_documented() {
    assert_page_states(
        "divergences/defect-n.md",
        &[
            "Ensembl VEP's own sequence path refuses the same POS on a substitution (`Invalid \
             start '.' or end '0' coordinate`) and numifies it on an anchor-chopped indel",
            "A POS that is not a number gets the sequence path's own text, `Invalid start \
             '<POS>' or end '<end>' coordinate`; a negative or fractional one gets `Invalid \
             POS: '<POS>'`",
            "An indel whose REF and ALT share their first base loses it, and the sequence \
             path's own `$start++` (`Parser/VCF.pm` lines 310 and 334) numifies the text the \
             same way, so `21 . x GA G` is annotated at `21:1`",
            "The fifth record, `bad_indel`, is skipped with the incremented text quoted: \
             `Invalid start 'abd' or end '1' coordinate`",
            "The sixth, `num_prefix`, is annotated: its POS numifies to 12 and the chop moves \
             it to 13, so Ensembl VEP writes the deletion `A/-` there",
        ],
    );
    assert_page_states(
        "concordance-provenance/116.2-divergence-classes.md",
        &[
            "the sequence path's own `$start++` numifies the text the same way: `21 12abc x GA \
             G` is `21:13 - intergenic_variant`",
            "vep-rs refuses every form: `Invalid start '<POS>' or end '<end>' coordinate` for \
             a POS that is not a number, `Invalid POS` for a negative or fractional one",
        ],
    );
}

/// `skipped-variant-semantics`: a VCF indel whose POS is not a number is skipped
/// with the POS as written. Ensembl VEP's anchor chop runs `$start++` on the
/// text first (Parser/VCF.pm:328-334), so it skips `21 abc bad_indel GA G` as
/// `Invalid start 'abd' or end '1' coordinate` and, `12abc` numified to 13,
/// passes `21 12abc num_prefix GA G` through `validate_vf` (Parser.pm:475) and
/// writes `num_prefix 21:13 - - - - intergenic_variant`. `--format vcf` is
/// given because neither engine detects the format from such a first line.
#[test]
fn anchor_chopped_indel_with_non_numeric_pos_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("badpos.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\tabc\tbad_indel\tGA\tG\t.\t.\t.\n\
         21\t12abc\tnum_prefix\tGA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.path().join("badpos.out");
    let cache = format!(
        "{}/../../tests/golden/116/GRCh38-identifiers/json_cache",
        env!("CARGO_MANIFEST_DIR")
    );
    let status = common::isolated_vep(dir.path())
        .stdin(std::process::Stdio::null())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--json_cache", &cache])
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args([
            "--format",
            "vcf",
            "--no_stats",
            "--quiet",
            "--force_overwrite",
        ])
        .status()
        .unwrap();
    assert!(status.success(), "{status}");
    let written = std::fs::read_to_string(&out).unwrap();
    assert!(written.lines().all(|l| l.starts_with('#')), "{written}");
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        "WARNING: line 1 skipped (21 abc bad_indel GA G . . .): Invalid start 'abc' or end '1' \
         coordinate\n\
         WARNING: line 2 skipped (21 12abc num_prefix GA G . . .): Invalid start '12abc' or end \
         '13' coordinate\n"
    );
}

/// `lookup_ref`: a contig the FASTA lacks skips the record at every position,
/// where BaseVEP.pm:424-434 gives such a contig a placeholder slice of length 1,
/// so FastaSequence.pm:597-601 returns `N` for a one-base record at position 1
/// and Parser.pm:645 annotates it as `N/<alt>`.
#[test]
fn lookup_ref_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Under `--lookup_ref` a contig the FASTA lacks skips the record at every position \
             with the could-not-fetch warning",
            "where Ensembl VEP's one-base placeholder slice (`BaseVEP.pm` lines 424-434) gives a \
             one-base record at position 1 the reference `N`",
        ],
    );
}

/// `skipped_variants_file`, `skipped-variant-reporting`,
/// `skipped-variant-warning-lines`: the header names standard input `STDIN` in
/// any case and under `--input_data`, the record column is 40 characters, and
/// the line a skip warning quotes is cut at a character boundary, where
/// BaseRunner.pm:371-376 prints the `input_file` parameter as typed (the `*IN`
/// glob of Runner.pm:818-820 under `--input_data`), `%-40.40s` counts the bytes
/// of the undecoded line and BaseVEP.pm:850's `substr` cuts byte 46 inside a
/// character.
#[test]
fn skipped_variants_file_exclusion_is_documented() {
    assert_page_states(
        "divergences/skipped-variants-file-header.md",
        &[
            "The header prints the `input_file` parameter as typed",
            "which prints as `*Bio::EnsEMBL::VEP::Runner::IN` (`Runner.pm` lines 818-820)",
            "so `length`, `substr` and `%-40.40s` count bytes",
            "When byte 46 falls inside a character, the cut splits it and the warning quotes the \
             character's lone lead byte",
            "The quoted line is cut by bytes too, at the last character boundary at or before \
             byte 46, so a line whose byte 46 falls inside a character loses that whole \
             character and the stray lead byte is never written",
            "The header names standard input `STDIN` whatever case `-i` spelled it in, and \
             under `--input_data` it names `STDIN`, or the `--input_file` value when one \
             stands beside the string",
            "The record column is the line's first 40 characters, padded to 40 characters, so \
             a line with a multi-byte character pads to a different width",
            "On a line of ASCII characters under `-i <path>`, `-i STDIN` or no `-i`, the two \
             files are byte-identical",
        ],
    );
}

/// `nearest`, `nearest-algorithm` and `cache-transcript-gene-tss`: every
/// protein-coding transcript of a shared start site is listed in stable-id order,
/// where Ensembl VEP prints the one object each tree search returns
/// (TranscriptTree.pm:233-236, AnnotationType/Transcript.pm:707-711).
#[test]
fn nearest_exclusion_is_documented() {
    assert_flag_row_states(
        "cli/output-control.md",
        "--nearest",
        &[
            "every transcript starting at a reported site is listed in stable-id order where \
             Ensembl VEP lists one",
        ],
    );
    assert_page_states(
        "divergences/nearest-shared-site.md",
        &[
            "A site several transcripts share is several points at one coordinate, and each \
             search returns one of them",
            "Nothing in Ensembl VEP's source chooses which",
            "The object depends on the order of the table file's lines, which `readdir` gives \
             (`Cache/Transcript.pm` line 297), and on the direction of the search that found \
             the site",
            "vep-rs lists every transcript starting at a winning site, in stable-id order",
            "so the two engines agree wherever the transcripts of the site belong to one gene",
            "Where transcripts of several genes share a site, Ensembl VEP prints the gene of \
             the one object and vep-rs each gene once",
        ],
    );
}

/// `variation-coordinates`: a native cache keeps the export's right-trimmed
/// indel alleles at the line's POS (ensembl-variation VariationFeature.pm:2390,
/// `trim_right`), so its row is shorter than the Ensembl cache's database row
/// (Database/Variation.pm:170-171) and an input anchored on the run's first
/// base co-locates with the native record alone.
#[test]
fn variation_coordinates_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "A native cache holds the trimmed form, which minimises from the right end to the \
             same alleles and position as the Ensembl cache's record",
            "So an input spelled at the left end of the run (the VCF convention, `AC/A` at \
             5030240) co-locates with it exactly as Ensembl VEP co-locates it",
            "The trimmed form also reaches that base and the Ensembl cache's record does not, so \
             the input co-locates with the native record alone",
            "An input spelled at the run's right end co-locates with neither, its start lying \
             outside the one-base window the start-indexed files answer",
        ],
    );
}

/// `vep-environment-variables`: a `VEP_<FLAG>` value never overwrites a typed
/// value, and one the parser refuses is an argument syntax error, where
/// Config.pm:1158's `||=` replaces a typed value Perl reads as false and
/// `is_valid_param` (Config.pm:1046-1047) drops the variable.
#[test]
fn vep_environment_variables_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "`VEP_<FLAG>=<value>` sets `--<flag>`, the name lower-cased, when the value is not \
             empty",
            "or the environment never overwrites one a higher source gives, even a `0`, an empty \
             value or a bare `--cache`",
            "Ensembl VEP's `||=` (`Config.pm:1101, 1158`) lets a lower source overwrite a value \
             Perl reads as false",
            "A value from a file or the environment that the parser refuses (`fork many`, an \
             empty `format`) is an argument syntax error (exit 2)",
            "Ensembl VEP's `is_valid_param` drops the line (`Config.pm:1046-1047`)",
        ],
    );
}

/// `option-prefix-abbreviation`: a prefix naming one flag of this version
/// resolves to it even where Ensembl VEP's larger `VEP_PARAMS` list
/// (Config.pm:73-276) makes it ambiguous, a token `pass_through` (Config.pm:62)
/// leaves in `@ARGV` for vep:30 to refuse.
#[test]
fn option_prefix_abbreviation_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A prefix is resolved against this version's flags first",
            "A prefix that names one flag of this version resolves to it even where Ensembl \
             VEP's larger set leaves it ambiguous",
            "A prefix of several flags that Getopt::Long resolves no better (`--fo`), or a name \
             no option has, is an argument syntax error (exit 2)",
        ],
    );
}

/// `spdi-input-format`: SPDI notation is resolved offline from the cache and
/// the FASTA, where Ensembl VEP reads the format with a database alone
/// (Parser/SPDI.pm:96). A version-less accession, which the core slice adaptor
/// resolves as a synonym prefix (SliceAdaptor.pm:2871-2884), a span past the
/// chromosome end, which it pads with `N` (FastaSequence.pm:619-639), and a
/// deletion of length 0, which the variation adaptor keeps as the feature `/-`
/// (VariationFeatureAdaptor.pm:2563-2566), are skipped.
#[test]
fn spdi_input_format_exclusion_is_documented() {
    assert_page_states(
        "divergences/spdi-input.md",
        &[
            "Ensembl VEP reads SPDI notation (`--format spdi`) only with a database \
             connection",
            "vep-rs resolves the same notations from the cache and the reference FASTA",
            "Three readings differ.",
            "queries the synonym table again with the name as a prefix",
            "vep-rs matches a name against the cache's chromosomes and their synonyms \
             exactly, regardless of ASCII case, never as a prefix",
            "tests no coordinate against the sequence's length",
            "so the record is annotated at a coordinate the chromosome lacks",
            "refuses a span the FASTA does not hold",
            "The adaptor keeps the feature `/-` from the position plus one to the position \
             (line 2566), a feature with no allele to annotate.",
            "vep-rs refuses the notation and skips the record, naming the empty span",
            "Deleted sequence length (0) leaves nothing to delete",
        ],
    );
}

/// `dir_cache`, `cache`, `json_cache`, `synonyms`, `check_svs`: a leading `~`
/// in the path is expanded here, where Ensembl VEP passes the value as typed
/// (AnnotationSourceAdaptor.pm:127 hands `dir_cache` to CacheDir.pm, whose `:318`
/// finds no directory named `~`; `check_svs` takes no value there, Config.pm:125).
#[test]
fn tilde_expansion_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A leading `~` in `--config`, `--dir`, `--dir_cache`, `--cache <dir>`, \
             `--json_cache`, `--synonyms`, `--check_svs` and a file's `dir` is expanded to the \
             home directory",
            "Ensembl VEP leaves the expansion to the shell",
            "A `~/` the shell did not expand (quoted, or from a file or the environment) in \
             `--config`, `--dir`, `--dir_cache`, `--cache`, `--synonyms` or `--full_cache_dir` \
             (`--json_cache` here) is passed on as typed there and names no configuration \
             file, cache directory or synonyms file",
            "The plugins directory under `--dir` alone is globbed there (`Runner.pm:908`)",
            "`--check_svs` takes no value there (`Config.pm:125`)",
        ],
    );
}

/// `option_rule_tables`: a comma-joined or trailing-comma value of a `%VALID`
/// scalar flag is an argument syntax error here, where Config.pm:946 splits the
/// value on commas before the `%VALID` test and OutputFactory.pm:1808 reads
/// `--sift b,` by its first letter; so is a `0` or empty value, which Config.pm:904-905
/// (`return 0 unless defined($config->{$flag}); return ... : $config->{$flag};`)
/// and :945 (`grep {$self->_is_flag_active($config, $_)} keys %VALID`) read as the
/// flag unset and OutputFactory.pm:1805 (`if (my $opt = $self->{$lc_tool})`) and
/// AnnotationType/Transcript.pm:170 (`if(my $nearest_type = $self->{nearest})`)
/// print no column for.
#[test]
fn option_rule_tables_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A trailing comma or a comma-joined value of `--sift`, `--polyphen`, `--terms`, \
             `--nearest`, `--compress_output` or `--format` (`--sift b,`) is an argument syntax \
             error (exit 2)",
            "Ensembl VEP's `check_config` splits such a value on commas before its `%VALID` test \
             (`Config.pm:946`) and reads the first token by prefix (`OutputFactory.pm:1808`)",
            "A `0` or empty value of such a flag (`--sift 0`, `--sift ''`, `--polyphen 0`, \
             `--nearest 0`) is an argument syntax error too (exit 2)",
            "Ensembl VEP reads a value Perl reads as false as the flag unset \
             (`Config.pm:904-905, 945`) and annotates without the column \
             (`OutputFactory.pm:1805`, `AnnotationType/Transcript.pm:170`)",
        ],
    );
}

/// `tab-column-order-and-placeholder`: a plugin header key the registry lists
/// twice (one plugin named twice, Runner.pm:925 and :994 loading one instance per
/// `--plugin`, or two plugins sharing a key) is one tab column, where Tab.pm:195
/// appends every key it is handed.
#[test]
fn tab_column_order_exclusion_is_documented() {
    assert_page_states(
        "output/vep-and-tab.md",
        &[
            "A plugin header key listed twice, the plugin named twice on the command line or \
             two plugins sharing a key, prints one column, where Ensembl VEP prints one per \
             instance",
        ],
    );
}

/// `hgvs-genomic-input`: a `g.` or `m.` span past the sequence's last base is
/// skipped, where Ensembl VEP's slice read pads the missing bases with `N`
/// (SequenceAdaptor.pm:148-152, :222-233) and annotates an intergenic row; and a
/// version-less accession (`NC_000021:g.30282169C>A`) is skipped by the exact
/// record lookup of `genomic.rs` (`reference.record_name`, the FASTA's records
/// and their aliases), where VariationFeatureAdaptor.pm:1889
/// (`fetch_by_region('chromosome', $reference) || fetch_by_region(undef,
/// $reference)`) reaches the slice adaptor's synonym retry with the name as a
/// prefix (SliceAdaptor.pm:2870-2884, `"$escaped_seq_region_name%"`) and
/// annotates the record as `NC_000021.9`'s.
#[test]
fn hgvs_genomic_input_exclusion_is_documented() {
    assert_page_states(
        "divergences/spdi-input.md",
        &[
            "An HGVS genomic notation (`g.`, `m.`) whose span runs past the last base is read \
             the same way by Ensembl VEP",
            "pads with `N` (`SequenceAdaptor.pm` lines 148-152 and 222-233)",
            "vep-rs reads the genome before it reads the shape and skips each with `Unable to \
             parse HGVS notation '<notation>': The reference FASTA has no bases at \
             21:46709984-46709984`",
            "An HGVS genomic notation reads the same way: `fetch_by_hgvs_notation` fetches its \
             slice through the same adaptor (`VariationFeatureAdaptor.pm` line 1889), so Ensembl \
             VEP annotates `NC_000021:g.30282169C>A` as it annotates \
             `NC_000021.9:g.30282169C>A`",
            "vep-rs matches the whole name, in any letter case, against the FASTA's records and \
             their aliases, never a prefix, and skips the notation",
            "Unable to parse HGVS notation 'NC_000021:g.30282169C>A': The reference FASTA has no \
             sequence named NC_000021",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "resolves a version-less accession as a synonym prefix, in an HGVS genomic notation \
             (`NC_000021:g.30282169C>A`) as in an SPDI one",
        ],
    );
    assert_flag_row_states(
        "cli/input-output.md",
        "--format",
        &[
            "an HGVS genomic (g., m.) span past the last base with `The reference FASTA has no bases at <name>:<start>-<end>`",
            "An SPDI or HGVS sequence name is matched whole, in any letter case, against the cache's chromosomes, the FASTA's records and their synonyms, never as a prefix, so a version-less accession (`NC_000021`) is skipped, in SPDI with `Sequence name <name> not valid` and in HGVS with `The reference FASTA has no sequence named <name>`, where Ensembl VEP's slice adaptor resolves the name as a prefix of a synonym (`NC_000021.9`) and annotates the record",
        ],
    );
}

/// `output_file`: the `-o` value names the file as given, where BaseRunner.pm:226
/// opens it with Perl's two-argument `open`, which strips surrounding whitespace
/// and reads a leading `>` as the append mode.
#[test]
fn output_file_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "The `-o` value names the output file as given",
            "Ensembl VEP opens it with Perl's two-argument `open` (`BaseRunner.pm:226`), which \
             strips surrounding whitespace and reads a leading `>` as the append mode",
        ],
    );
}

/// `assembly`: `--assembly` is compared exactly and case-sensitively, as
/// CacheDir.pm:347 and :407 compare it with `ne`, with one addition Ensembl VEP
/// lacks: `hg19` and `hg38` read as `GRCh37` and `GRCh38`.
#[test]
fn assembly_alias_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "`--assembly` is compared exactly and case-sensitively, as Ensembl VEP compares it, \
             with one addition: `hg19` and `hg38` are read as `GRCh37` and `GRCh38`",
        ],
    );
}

/// `freq_filter` and `check_frequency` under `--minimal`: a split record none of
/// whose alleles writes a row stays one record under a `--check_frequency` drop,
/// where Ensembl VEP's buffer, its rejoin cancelled by the drop
/// (BaseCacheVariation.pm:365-368 `if($filtered_count) { $buffer->reset_buffer();`,
/// InputBuffer.pm:535-537), writes each of its features as its own; and the reset
/// reaches the dropped record's chromosome run within a `--buffer_size` batch of
/// input lines here (`runs_reset_by_a_drop` in `pipeline.rs`), where it reaches
/// Ensembl VEP's buffer of one chromosome's records (InputBuffer.pm:191-199
/// `if($prev_chr && $vf->{chr} ne $prev_chr) { ... return $buffer; }`) or under
/// `--fork` a chunk of it (Runner.pm:486 `int($buffer_size / (2 * $fork_number))`).
#[test]
fn minimal_frequency_filter_exclusion_is_documented() {
    assert_page_states(
        "output/known-variants.md",
        &[
            "A split record none of whose alleles writes a row (`GG *,GG`) stays one VCF line \
             and one JSON object under a reset, where Ensembl VEP writes each of its features \
             as its own",
            "The reach of the reset differs. Here it is the dropped record's run of one \
             chromosome within a batch of `--buffer_size` input lines counted from the file's \
             first data line, skipped lines included",
            "In Ensembl VEP it is the buffer, one chromosome's records up to `--buffer_size` \
             (`InputBuffer.pm` lines 191-199), or under `--fork` a chunk of the buffer holding \
             at most `--buffer_size` over twice the fork count (`Runner.pm` line 486)",
            "Where the units differ, one engine rejoins a split record the other resets",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "a split record none of whose alleles writes a row (`GG *,GG`) stays one VCF line \
             and one JSON object here",
            "Ensembl VEP's reach is its buffer of one chromosome's records up to \
             `--buffer_size`, or under `--fork` a chunk of it, and vep-rs's the dropped \
             record's chromosome run within a batch of `--buffer_size` input lines counted \
             from the file's first line, so where the units differ one engine rejoins a split \
             record the other resets",
        ],
    );
}

/// `shift_hgvs`: `--shift_hgvs` takes `0` or `1`, where Ensembl VEP's `=i`
/// (Config.pm:201) takes any integer and BaseRunner.pm:498 matches it against
/// `/(0|1)/` as a substring, so `2` runs with the default and `10` shifts
/// HGVSc and HGVSp (TranscriptVariationAllele.pm:1441 tests `!= 0`) but not
/// HGVSg (VariationFeature.pm:1986 tests `== 1`).
#[test]
fn shift_hgvs_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "`--shift_hgvs` takes `0` or `1`; any other value is an argument syntax error \
             (exit 2)",
            "Ensembl VEP's `shift_hgvs=i` (`Config.pm:201`) takes any integer and \
             `BaseRunner.pm:498` matches it against `/(0|1)/` as a substring",
            "So `2` leaves the default shift on, and `10` sets the shifting variables to 10, \
             which shifts HGVSc and HGVSp (tested against `0`) and not HGVSg (tested against \
             `1`)",
        ],
    );
}

/// `fasta_dir`: two files tied on both preferences are ordered by name, the
/// first taken (`fasta_dir_match` in vep-fasta's `lib.rs`), where
/// BaseVEP.pm:500-509 sorts on the two preferences alone (`$a !~ /toplevel/ <=> $b
/// !~ /toplevel/ || $a !~ /gz$/ <=> $b !~ /gz$/`) over `readdir`'s order and :513
/// takes `$files[0]`. The filter reads the assembly the cache resolves to
/// (`fasta_source` in `runner.rs` passes the resolved cache's), where :497 reads
/// `param('assembly')`, which `--assembly`, the `--grch37` option set
/// (Config.pm:399-403) or the layout scan of CacheDir.pm:303 (its :356) sets and
/// BaseRunner.pm:139 fills from the database `--offline` never reaches, so under
/// `--offline --full_cache_dir` with neither flag it is unset and BaseVEP.pm:508's
/// `/$assembly/i` is an empty pattern, Perl's "last successful match", and every
/// species match passes.
#[test]
fn fasta_dir_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Two files tied on both preferences are ordered by name, the first taken, where \
             Ensembl VEP's sort on the two preferences alone keeps the directory listing's \
             order (`BaseVEP.pm` lines 500-509), so on such a directory the two engines can \
             read different files",
            "Here the filter reads the same assembly as the redirect, where Ensembl VEP reads \
             its assembly parameter alone (line 497)",
            "With that parameter unset, Ensembl VEP's filter passes every file naming the \
             species (line 508), and its preferences and the listing's order choose among \
             every assembly's files",
        ],
    );
    assert_flag_row_states(
        "cli/cache-species-runtime.md",
        "--fasta_dir",
        &[
            "a toplevel file before any other, a bgzipped file before a plain one and the first \
             by name of two tied on both (where Ensembl VEP keeps the directory listing's order)",
        ],
    );
}

/// `hgvs-transcript-input`: the five readings of the transcript branch that
/// differ. An `NM_` or `XM_` reference resolves against the
/// cache's RefSeq transcripts alone (`resolve` in vep-io's `hgvs/transcript.rs`),
/// the lookup reads the run's transcript set, a gene name resolves in
/// database-id order, and a description whose position the expression leaves
/// undefined (`c.-A>G`) is refused as one it could not parse
/// (`genomic_position` there), where VariationFeatureAdaptor.pm:1829-1849 and
/// :2146-2226 read the RefSeq and core databases and
/// `_parse_hgvs_transcript_position` (:1608, called at :1858) carries the
/// undefined position on. A notation every core group refuses is skipped with
/// the first group's message alone (`first_refusal` in `variant` there, one
/// line through `unable` in `hgvs.rs`), where Parser/HGVS.pm collects every
/// group's distinct exception text (:199, `merge_arrays` over the errors list,
/// unique strings through Utils.pm:306-307) and joins the whole list under the
/// one header (:208 copies it into the message list, :214 the `join`); the
/// first error alone (`$errors[0]`, :207) drives only the `--max_sv_size`
/// collapse of :209-212. The protein branch keeps the same first refusal
/// across the groups (`first_refusal` in `protein.rs`), where the parser's group
/// loop (:171) is not notation-type specific, so the page's protein list states
/// the reading too.
#[test]
fn hgvs_transcript_input_exclusion_is_documented() {
    assert_page_states(
        "divergences/hgvs-transcript-input.md",
        &[
            "Five readings differ:",
            "An `NM_` or `XM_` reference resolves against the cache's RefSeq transcripts, present \
             under `--refseq` or `--merged`, and is otherwise refused as a reference no \
             transcript carries (`Could not get a Transcript object`)",
            "Ensembl VEP reads such a reference from its RefSeq database whatever cache the run \
             annotates with",
            "The lookup reads the run's transcript set, so under `--gencode_basic` or \
             `--gencode_primary` a transcript outside the set is not found, where Ensembl VEP's \
             resolver reads its whole database",
            "A gene name resolves among the transcripts on the assembled chromosomes and the \
             unplaced scaffolds, in database-id order, where Ensembl VEP's `is_reference` test \
             reads the database",
            "A description whose position its regular expression leaves undefined (`c.-A>G`, a \
             sign with no digits) is refused as one it could not parse, where Ensembl VEP \
             carries the undefined position on",
            "A notation every core group refuses is skipped with the first group's message \
             alone, where Ensembl VEP joins every group's distinct exception text under its \
             `Unable to parse HGVS notation` warning (`Parser/HGVS.pm` lines 199, 208 and 214)",
            "So in its database run a refused `ENST` notation also carries the `otherfeatures` \
             group's `Could not get a Transcript object` message, and a refused gene name that \
             group's own refusal",
            "Six readings of the protein branch differ:",
            "A notation every core group refuses is skipped with the first group's message alone, \
             as in the transcript branch, where Ensembl VEP joins every group's distinct exception \
             text under its warning",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &["five readings of the transcript branch and six of the protein branch differ"],
    );
}

/// `individual_zyg`, `individual` and `individual-modes`: an individual's copy of
/// a bi-allelic record trimmed past the anchor chop, an indel whose written
/// alleles differ in length (Parser.pm:860-881, `elsif
/// (length($ref_allele_string) != length($alt)) { $is_non_minimised_indel = 1;`
/// at :876-877, with or without `--minimal`) or under `--minimal` a same-length
/// record (`GA/GT`, `AG/TG`; Parser.pm:852), is compared chopped against chopped
/// (`CalledGenotype::read` in `individuals.rs`), so a `0/0` sample is
/// reference-homozygous and `--individual_zyg` prints its reference row
/// (`zygosity_record`), where Parser/VCF.pm:691-692 take the minimised `A` as
/// `$ref` and :716 and :720 compare the chop-level `GA` with it (`grep {$_ ne
/// $ref && $_ !~ /\*/} @bits`), so Ensembl VEP prints the alternate row with
/// `ZYG=S1:HOM`. `--individual` refuses `--minimal` (Config.pm:626 `individual
/// => [qw(minimal individual_zyg)]`), so the page and the index name the record
/// shape, not the flag, as the trigger.
#[test]
fn individual_minimal_trimmed_record_exclusion_is_documented() {
    assert_page_states(
        "divergences/individual-trimmable-indel.md",
        &[
            "An individual's copy of a record trimmed past the anchor chop",
            "A difference in the alleles annotated for an individual, on a record shape no \
             measured per-individual run reads; nothing is set aside",
            "The shape is a bi-allelic indel whose anchor-chopped alleles still share a base, \
             which Ensembl VEP trims again with or without `--minimal` (`Parser.pm` lines \
             860-881), or a same-length record such as `GA/GT` under `--minimal`",
            "Ensembl VEP minimises a bi-allelic indel after the anchor chop whatever the flags \
             (`Parser.pm` lines 860-881, which pass a record whose written alleles differ in \
             length to `minimise_alleles`: `ACAC/AC` becomes `AC/-` two bases on)",
            "`--individual` refuses `--minimal` unless `--safe` lifts the check (`Config.pm` \
             lines 626 and 961), so that pass alone trims the records it copies",
            "Under `--minimal` it minimises every bi-allelic record whose alleles still share a \
             base at either end (line 852: `GA/GT` becomes `A/T` a base on, `AG/TG` becomes \
             `A/T` at the written position)",
            "vep-rs compares the chopped alleles with the chopped reference, so the `0/0` \
             individual is reference-homozygous",
            "as does a same-length record without `--minimal`",
            "So `21 34716539 GA GT` with one sample `0/0` under `--individual_zyg --minimal` is \
             the `T` row with `ZYG=S1:HOM` and `ALLELE_NUM=1` in Ensembl VEP and the `A` row \
             with `ZYG=S1:HOMREF` and `ALLELE_NUM=0` here, `MINIMISED=1` on both",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP builds an individual's copy of a bi-allelic record it trims past the \
             anchor chop, an indel such as `ACAC/AC` with or without `--minimal` or a \
             same-length record such as `GA/GT` under it, from the chopped alleles against the \
             minimised reference",
            "nothing: no measured per-individual run reads such a record",
            "[An individual's copy of a record trimmed past the anchor chop]\
             (divergences/individual-trimmable-indel.md)",
        ],
    );
    assert_page_states(
        "output-formats.md",
        &[
            "an indel such as `ACAC/AC` with or without `--minimal` and a same-length record \
             such as `GA/GT` under it",
            "The reading above of a record trimmed past the anchor chop, an indel or under \
             `--minimal` a same-length record such as `GA/GT`, applies to the list too. Ensembl \
             VEP reads a `0/0` sample of it as `HOM` and keeps the alternate allele's rows, where \
             vep-rs reads `HOMREF` and annotates the reference allele",
        ],
    );
}

/// `mirna` under `--minimal`: a split allele's hairpin terms come from its own
/// span in every format (`calculate_consequences` in vep-effects'
/// `consequences.rs` runs `mirna_structure_terms` on the allele it analyses),
/// where the default and tab factories rejoin the record (OutputFactory.pm:294)
/// and then clear the cDNA coordinates (:297-298 `reset_shifted_positions`,
/// TranscriptVariationAllele.pm:515-517 `delete($tv->{cdna_start});
/// delete($tv->{cdna_end});`), so the miRNA block (OutputFactory.pm:1630 `my
/// ($cdna_start, $cdna_end) = ($tv->cdna_start, $tv->cdna_end);`) reads the
/// rejoined record's span.
#[test]
fn mirna_minimal_split_exclusion_is_documented() {
    assert_page_states(
        "divergences/minimal-split-hgvsp.md",
        &[
            "and the `miRNA` terms of a split allele on a miRNA transcript",
            "The `miRNA` column reads the same cleared coordinates (`OutputFactory.pm` line 1630 \
             takes `cdna_start` and `cdna_end`, which `clear_shifting_variables` deleted, \
             `TranscriptVariationAllele.pm` lines 515-517). So the default and tab factories \
             compute a split allele's hairpin terms over the record's restored span, and the VCF \
             and JSON factories over the allele's own",
            "For `21 8987003 multi ACT A,ACTTTA` under `--mirna --minimal`, the `TTA` allele's \
             row on the miRNA transcript prints `miRNA=miRNA_loop` in Ensembl VEP's default \
             output and `miRNA=miRNA_loop,miRNA_stem` here; the VCF and JSON of both engines \
             print `miRNA_loop&miRNA_stem`",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "under `--mirna` the default and tab outputs compute a split allele's hairpin terms \
             over the record's span",
            "vep-rs prints the allele's own positions, terms, amino acids, miRNA terms and \
             notations in every format",
        ],
    );
}

/// `region-input-format`: a region line whose structural word `%SO_TERMS` lacks
/// is warned about and skipped (`RegionParser::parse_record` in vep-io's
/// `region.rs` fails the line with `IoError::SkippedLine`) and the run goes on,
/// where Parser/Region.pm:205-206 (`my $so_term = $self->get_SO_term($allele);
/// return [] unless defined $so_term;`) builds no feature after the warning of
/// Parser.pm:725, so a line that opens a buffer fill ends Ensembl VEP's run
/// (Parser.pm:186-187, InputBuffer.pm:212, Runner.pm:424).
#[test]
fn region_unsupported_structural_word_exclusion_is_documented() {
    assert_page_states(
        "divergences/dont-skip-cannot-keep.md",
        &[
            "A region line whose structural word is not a supported type \
             (`21:14221000-14222000/FOO`) is warned about by both engines, `FOO is not a \
             supported structural variant type`",
            "Ensembl VEP then builds no feature for it (`Parser/Region.pm` lines 205-206), so \
             when the line opens a fill it ends the run the same way, with the warning and no \
             row written, where vep-rs goes on to the lines after it",
        ],
    );
}

/// `chromosome-synonyms`: a name two hops from a cache chromosome under a chained
/// or whitespace-led `--synonyms` table is skipped as not found whatever came
/// before it (`known_as` in `pipeline.rs` tests one hop and `validate.rs` writes
/// the skip), where Parser.pm:782 (`$valid->{$_} = 1 for map {keys
/// %{$synonyms->{$_} || {}}} keys %$valid;`) widens the parser's persistent valid
/// hash (:328-331) in place at every miss, so a later record on the two-hop name
/// passes :784 and BaseVEP.pm:711-716 resolve one hop to nothing:
/// `intergenic_variant`.
#[test]
fn chromosome_synonyms_two_hop_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Three cases differ from Ensembl VEP, which checks a record against the cache's \
             names widened by their synonyms, a set it widens again in place at every miss \
             (`Parser.pm` line 782), and then looks the unwidened name up one hop",
            "The third is a name two hops from a cache chromosome under a `--synonyms` table \
             whose lines chain or open with whitespace (`21 ALPHA` then `ALPHA BETA`: `BETA`), \
             admitted once an earlier record missed",
            "Each passes its check there and annotates as intergenic; in vep-rs a bare `M` \
             resolves to `MT` and every other such name is skipped as not found, whatever \
             record came before it",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "a set each miss widens again in place, and then resolves the unwidened name one hop",
            "or a two-hop synonym of a chained or whitespace-led `--synonyms` table after an \
             earlier miss gets rows from one engine that the other does not write",
        ],
    );
}

/// `fasta`: the pseudoautosomal `Y` to `X` redirect is armed from the assembly
/// the cache resolves to (`open_fasta` in `runner.rs` sets the opened FASTA's
/// assembly to the one `cache_resolve::resolve` took, `info.json` first), where
/// BaseVEP.pm:525 (`-ASSEMBLY => $self->param('assembly'),`) and
/// FastaSequence.pm:254 (`$pars ||= $PARS{$assembly} if $assembly && ...`) arm it
/// from the assembly parameter, set by `--assembly`, `--grch37`
/// (Config.pm:399-403) or the layout scan (CacheDir.pm:303-356), which runs for
/// `--dir_cache` and `--dir` and not under `--full_cache_dir`, so
/// `--offline --full_cache_dir` without either flag leaves it unset, the redirect
/// unarmed there (FastaSequence.pm:550 `$Bio::EnsEMBL::Slice::_fasta_PARs &&`).
///
/// A FASTA record is also reached through the aliases `contig_aliases` in
/// vep-fasta's `lib.rs` builds in (`chr` added or removed, and the
/// `MT`/`M`/`chrM`/`chrMT` set, matched after `to_ascii_lowercase`), where
/// BaseVEP.pm:413 maps a chromosome onto the FASTA's record names through
/// `get_source_chr_name` (:689-738: the name, a synonym, `chr` added or removed)
/// and FastaSequence.pm:380-389 retries the same three at the read, so `MT`
/// never reaches `chrM` there; the unreached name gets the one-base placeholder
/// of BaseVEP.pm:425-426 and `--check_ref` skips its records (Parser.pm:607-610,
/// :624-631).
///
/// The redirected `Y` read on a FASTA without `X` warns in neither engine:
/// FastaSequence.pm:599 calls `warning(...)` at the default level 2000 while
/// BaseRunner.pm:491 sets `Bio::EnsEMBL::Utils::Exception::verbose(1999)` for
/// the run, and `warning` returns when the verbosity is below the level.
#[test]
fn fasta_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "The redirect is armed by the assembly the cache resolves to, its `info.json` first, \
             then its directory name, then `--assembly`; Ensembl VEP arms it from its assembly \
             parameter (`BaseVEP.pm` line 525, `Utils/FastaSequence.pm` line 254)",
            "That parameter is set by `--assembly`, `--grch37` (`Config.pm` lines 399-403) or \
             the layout scan (`CacheDir.pm` lines 303-356), which runs for `--dir_cache` and \
             `--dir`, not under `--full_cache_dir`",
            "Under `--offline --full_cache_dir` (`--json_cache` here) without `--assembly` or \
             `--grch37` it stays unset, so Ensembl VEP reads a `Y` span inside the region as \
             written and vep-rs reads `X`",
            "Beyond the table, a record is reached by its name with the prefix `chr` added or \
             removed and, for a record named `MT`, `M`, `chrM` or `chrMT`, by any of the four \
             names, matched in any letter case: aliases vep-rs builds in",
            "Ensembl VEP tries the exact name, a synonym and the `chr` variants alone \
             (`BaseVEP.pm` lines 413 and 689-738, `Utils/FastaSequence.pm` lines 380-389), so a \
             cache's `MT` over a FASTA record named `chrM`, with no synonym linking the two, \
             reaches no sequence there",
            "That chromosome gets the one-base placeholder slice described above (`BaseVEP.pm` \
             lines 425-426), so under `--check_ref` a record such as `MT 3308 T G` is skipped \
             with the sub-slice and mismatch warnings, where vep-rs annotates it",
            "Neither engine warns on the read: Ensembl VEP's `FASTA does not contain sequence \
             for X` (line 599) is raised at warning level 2000, above the verbosity of 1999 a run \
             sets (`BaseRunner.pm` line 491)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Under `--full_cache_dir` with neither `--assembly` nor `--grch37`, Ensembl VEP's \
             assembly stays unset, so its JSON prints `\"assembly_name\":null`, its \
             pseudoautosomal redirect of `Y` reads to `X` is not armed",
            "vep-rs takes the assembly the cache's `info.json` names, as it does under \
             `--json_cache`, prints it, arms the redirect from it",
            "on the FASTA side vep-rs also reaches a record by the name with the prefix `chr` \
             added or removed or, for a record named `MT`, `M`, `chrM` or `chrMT`, by any of the \
             four names in any letter case, where Ensembl VEP tries the exact name, a synonym \
             and the `chr` variants alone, so a record such as `MT 3308 T G` on a cache's `MT` \
             over a FASTA record named `chrM`, with no synonym linking them, annotates here under \
             `--check_ref` and is skipped there",
            "vep-rs returns the flanking `Y` bases with the `N`s, and neither engine prints a \
             warning for the read",
        ],
    );
}

/// `synonyms`: a `--synonyms` file holding a byte that is not UTF-8 is refused
/// (`ChromosomeSynonyms::load_from_file` in vep-core's `synonyms.rs` reads it
/// with `read_to_string`, worded `Could not read synonyms file <file>: stream
/// did not contain valid UTF-8`, exit 1), where Ensembl VEP reads the file as
/// bytes; and the flag's value is a file name as given (`load_synonyms` in
/// `cache_resolve.rs` opens it as a path, so `-`, `0`, an empty value, a
/// padded path and a `>`-led value are refused when no such file exists), where
/// BaseVEP.pm:643-644 (`if ($file) { open IN, $file or throw(...)`) reads `-` as
/// standard input, trims the whitespace and opens a `>`-led value for writing
/// through Perl's two-argument `open`, and a false value (`0`, empty) reads as
/// the flag unset at CacheDir.pm:254, which fills it with the cache's
/// `chr_synonyms.txt` before that `open` runs (Runner.pm:128, then :131).
#[test]
fn synonyms_file_encoding_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "A file holding a byte that is not UTF-8 is refused the same way, the reason \
             `stream did not contain valid UTF-8`, where Ensembl VEP reads the file as bytes",
            "The `--synonyms` value is a file name as given (`-`, `0`, an empty value, \
             surrounding whitespace, a leading `>` or `|`), where Ensembl VEP's two-argument \
             `open` reads `-` as standard input, trims the whitespace and takes a false value as \
             unset",
        ],
    );
    assert_page_states(
        "cli/compatibility.md",
        &[
            "The `--synonyms` value names its file as given too, `-`, `0`, the empty value, a \
             value with surrounding whitespace and a `>`-led value included",
            "Ensembl VEP's two-argument `open` (`BaseVEP.pm:644`) reads `-` as standard input, \
             trims the whitespace, and opens a `>`-led value for writing, emptying the file it \
             names and reading an empty table",
            "A value Perl reads as false leaves the flag unset there, so the cache's own file is \
             read (`CacheDir.pm:254`)",
        ],
    );
}

/// `chr`: a chromosome literally named `1-5` is refused silently under `--chr 1-5`
/// (`ChromosomeList::new` in `validate.rs` files the token under `ranges`), where
/// Parser.pm:471 (`grep {$vf->{chr} eq $_} @{$self->{chr}}`) admits it by
/// equality.
#[test]
fn chr_literal_range_name_exclusion_is_documented() {
    assert_page_states(
        "divergences/chr-range.md",
        &[
            "A chromosome literally named `1-5`, which Ensembl VEP admits under that token, is \
             refused silently here, as a record `--chr` leaves out is in both engines",
        ],
    );
}

/// `cache-valid-chromosomes`: a cache none of whose chromosome directories holds
/// a transcript shard is refused as it opens (`enumerate_transcript_shards` in
/// `json_cache.rs`, `ERROR: [cache] no transcript shards under <dir>
/// (serialiser_type <type>)`, exit 1), where Ensembl VEP lists every non-dot
/// directory as a valid chromosome whatever it holds (CacheDir.pm:431 `sort grep
/// {!/^\./ && -d $self->dir.'/'.$_} readdir DIR`), reads a region's file only
/// `if -e $file` (BaseSerialized.pm:261-263) and so writes every record as
/// `intergenic_variant`, exit 0.
#[test]
fn cache_without_a_transcript_shard_refusal_is_documented() {
    assert_page_states(
        "divergences/damaged-cache-shard.md",
        &[
            "Ensembl VEP's valid chromosomes are the cache's directories whose names do not start \
             with a dot, whatever they hold (`CacheDir.pm` line 431), and a region whose file \
             does not exist yields no feature (`AnnotationSource/Cache/BaseSerialized.pm` line \
             263, `if -e $file`)",
            "So a cache none of whose chromosome directories holds a transcript file runs to the \
             end: every record on a directory's name is written as `intergenic_variant`, exit 0, \
             and nothing reports the cache empty",
            "vep-rs lists the shards of every chromosome directory as the cache opens, before \
             the input is read, and refuses a cache that holds none, a dot-prefixed directory's \
             files not counted",
            "ERROR: [cache] no transcript shards under <cache>/transcripts (serialiser_type \
             storable)",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &[
            "`transcripts/` must also hold a shard under some chromosome directory: a cache with \
             none is refused, where Ensembl VEP writes every record as `intergenic_variant`",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP runs a cache none of whose chromosome directories holds a transcript \
             file, reads no feature from any region, and writes every record as \
             `intergenic_variant`, exit 0; vep-rs refuses the cache before the input is read, \
             `ERROR: [cache] no transcript shards under <dir>/transcripts (serialiser_type \
             <type>)`, exit 1",
        ],
    );
}

/// `hgvs`: under `--hgvs` a contig the FASTA lacks gives no HGVSc for a
/// substitution or a deletion, whose reference bases cannot be read, where
/// Ensembl VEP reads the transcript's slice as `N` per base
/// (FastaSequence.pm:597-601, `return 'N' x (($end - $start) + 1)`), hands that
/// sequence to `hgvs_variant_notation` as the reference
/// (TranscriptVariationAllele.pm:1474-1476) and cuts the reference allele from
/// it (Sequence.pm:526), so its substitution names `N` (`c.764N>T` for the
/// genome's `c.764G>T`) and its deletion reads as over the genome (`c.632del`);
/// an insertion's `c.307_308insC` and HGVSp agree, and neither engine warns
/// (`a_contig_the_fasta_lacks_gives_no_hgvsc_for_a_substitution_or_deletion` in
/// `input_hgvs.rs` runs the vep-rs side).
#[test]
fn hgvs_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Under `--hgvs` the same contig gives these notations; neither engine warns:",
            "| Record on the contig under `--hgvs` | Ensembl VEP, whose transcript slice reads \
             `N` per base there (`Utils/FastaSequence.pm` lines 597-601) | vep-rs, which reads no \
             reference base there |",
            "| A substitution | `c.764N>T`, the `N` standing for the genome's `G` | no HGVSc |",
            "| A deletion | `c.632del`, as over the genome | no HGVSc |",
            "| An insertion | `c.307_308insC` | `c.307_308insC` |",
            "| HGVSp | as over the genome | as over the genome |",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "[A contig the FASTA lacks under `--check_ref`, `--lookup_ref` or \
             `--hgvs`](cache-setup.md#reference-fasta)",
            "under `--hgvs` it reads a transcript's slice there as `N` per base and prints an \
             HGVSc whose reference base is `N` (`c.764N>T`)",
            "under `--hgvs` prints no HGVSc for a substitution or a deletion, where an \
             insertion's notation and HGVSp agree",
        ],
    );
}

/// `warning_file`, `skipped_variants_file`: both values are file names as given,
/// `-`, surrounding whitespace and a `>`-led value included (`target_for` in
/// `warnings.rs` builds a `PathBuf` from any value but `stderr`;
/// `reserve_skipped_variants_file` likewise), where Ensembl VEP's two-argument
/// `open` for writing, BaseVEP.pm:894 `$fh->open(">".$file) or throw(...)` and
/// BaseRunner.pm:344 `$fh->open(">$file_name") or`, sends `-` to standard
/// output, trims the whitespace and appends a `>`-led value to the file it names
/// (`warning_file_and_skipped_file_values_are_file_names_as_given` in
/// `cli_run_control.rs` runs the vep-rs side). A false `--skipped_variants_file`
/// value, `0` or empty, is a file name here too, where Runner.pm:176
/// `$self->dump_skipped_variants if $self->param('skipped_variants_file');`
/// writes nothing (`skipped_variants_file_false_value_is_a_file_name`).
#[test]
fn warning_file_as_given_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "The `--warning_file` and `--skipped_variants_file` values name their files as given \
             likewise, `-`, surrounding whitespace and a `>`-led value included",
            "Ensembl VEP opens each for writing with the two-argument `open` (`BaseVEP.pm:894`, \
             `BaseRunner.pm:344`), which sends `-` to standard output, trims the whitespace and \
             appends a `>`-led value to the file it names",
            "So `--warning_file -` prints the warnings on standard output there and writes a file \
             named `-` here, `--warning_file ' sp.txt'` writes `sp.txt` there and ` sp.txt` here, \
             and `--skipped_variants_file '>skipped.txt'` appends to `skipped.txt` there and \
             writes `>skipped.txt` here",
            "A `--skipped_variants_file` value Perl reads as false, `0` or the empty value, leaves \
             the flag unset there (`Runner.pm:176`), so no file is written. vep-rs writes a file \
             named `0`, and the empty value fails the run as it ends, the output written, with \
             `Could not write to skipped variants file` (exit 1). A false `--warning_file` value \
             is the default warnings file in both engines (`BaseVEP.pm:890`)",
        ],
    );
}

/// `force_overwrite`, `skipped_variants_file`: an existing skipped-variants file
/// is refused before the run (`reserve_skipped_variants_file` in `warnings.rs`),
/// where Ensembl VEP makes the check inside `finish` after the output is written
/// (BaseRunner.pm:337 `if(-e $file_name && !$self->param('force_overwrite')) {`,
/// reached from Runner.pm:176 `$self->dump_skipped_variants if
/// $self->param('skipped_variants_file');`), and both engines write the file as
/// the run ends (BaseRunner.pm:364-380).
#[test]
fn skipped_variants_file_timing_exclusion_is_documented() {
    assert_flag_row_states(
        "cli/input-output.md",
        "--skipped_variants_file",
        &[
            "An existing file is refused before the run without `--force_overwrite`, where \
             Ensembl VEP checks when the run ends; both engines write the file as the run ends",
        ],
    );
}

/// `individual`, `individual-modes`, `individual_zyg`: the sample names come from
/// the `#CHROM` line alone (`pipeline.rs`, `line.starts_with(b"#CHROM")` under
/// `cfg.individuals.is_some()`), where ensembl-io takes the header's last
/// single-`#` line (BaseVCF4.pm:96-97 `elsif ($line =~ /^#\s*(.+)$/) {
/// $self->{'metadata'}->{'header'} = [split("\t",$1)];`), so a `#note` line after
/// `#CHROM` empties Ensembl VEP's sample list: `--individual all` writes no row
/// where vep-rs writes the 28 rows of the `hb_multi` record
/// (`individual_all_reads_the_chrom_line_when_a_hash_line_follows_it` in
/// `individuals.rs`).
#[test]
fn hash_line_after_chrom_exclusion_is_documented() {
    assert_page_states(
        "divergences/vcf-hash-line-after-chrom.md",
        &[
            "Ensembl VEP reads a VCF header through ensembl-io, which takes every single-`#` line \
             as the column line in turn, so the last one stands (`BaseVCF4.pm` lines 96-97",
            "so a header whose `#CHROM` line is followed by `#note` names none",
            "Under `--individual all` the run writes no row and exits 0; under a list it stops \
             with Ensembl VEP's unknown-sample error, which names the list against an empty header \
             (`Parser/VCF.pm` lines 612-617); under `--individual_zyg all` every row is written \
             without `ZYG`",
            "The VCF output copies the `#CHROM` line among the meta lines and prints `#note` where \
             the column line belongs (`OutputFactory/VCF.pm` line 224)",
            "vep-rs reads the sample names from the `#CHROM` line alone and writes a single-`#` \
             line after it nowhere",
            "Under `--individual all` vep-rs writes the 28 rows the header without the stray line \
             gets, under `--individual S1` S1's four",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Fifty-three differences lie outside every class",
            "Ensembl VEP takes a VCF header's last single-`#` line as its column line, so a `#` \
             line after the `#CHROM` line leaves its sample list empty: `--individual all` writes \
             no row, a list is refused, `--individual_zyg` writes no `ZYG`, and the VCF output \
             prints that line where the column line belongs; vep-rs reads the samples from the \
             `#CHROM` line and writes the stray line nowhere",
        ],
    );
    assert_flag_row_states(
        "cli/input-output.md",
        "--individual",
        &[
            "a single-`#` line after that line is ignored, where Ensembl VEP takes the header's \
             last single-`#` line as the column line and so names no sample",
        ],
    );
}

/// `incomplete-del`: `SVLEN=NaN` numifies to 0 (`perl_numeric_f64` in
/// `sv_span.rs`), the deletion's end lands before its start and
/// `deletion_looks_incomplete` drops the record with its warning unless
/// `--dont_skip` keeps it at `21:30425650-30425651`, where Ensembl VEP's end is
/// `NaN` (ensembl-io BaseVCF4.pm:304-305 `$end = $self->get_start +
/// abs($svlen) - 1;`), its `$start > $end` test (Parser/VCF.pm:555) is false,
/// and it writes the Defect K row `21:30425651-NaN deletion - - -
/// intergenic_variant`.
#[test]
fn svlen_nan_deletion_exclusion_is_documented() {
    assert_page_states(
        "divergences/defect-k.md",
        &[
            "The `NaN` deletion, zero bases long, is dropped as incomplete with its warning unless \
             `--dont_skip` keeps its two-base interval",
            "The `NaN` deletion, zero bases long, is dropped as incomplete unless `--dont_skip` \
             keeps its two-base interval",
            "The `SVLEN=NaN` record, a zero-length deletion once numified, is dropped by default \
             with the warning `deletion looks incomplete`, where the reference wrote one \
             intergenic row; `--dont_skip` keeps it at the two-base interval `21:30425650-30425651` \
             (two transcript rows, `coding_sequence_variant` on ENST00000390690)",
            "The zero-length deletion is dropped as incomplete, or kept at a two-base interval \
             under `--dont_skip`",
            "and, under `--dont_skip`, at the two-base interval for `SVLEN=NaN`, \
             `21:30425650-30425651`; by default vep-rs writes no row for that record",
        ],
    );
    assert_page_states(
        "divergences/defect-c.md",
        &[
            "A `SVLEN` of `NaN` is the exception: Ensembl VEP's end is `NaN`, which no `start > \
             end` test catches, so it writes the record's intergenic row at `21:30425651-NaN` \
             ([Defect K](defect-k.md)); vep-rs reads the length as 0 and drops the record as \
             incomplete",
        ],
    );
    assert_page_states(
        "output/structural-variants.md",
        &[
            "A `SVLEN` of `NaN` is the exception: vep-rs reads it as 0 and drops the record, where \
             Ensembl VEP annotates it with the end `NaN` ([Defect K](../divergences/defect-k.md))",
        ],
    );
}

/// `spdi`, `regulatory`: on a regulatory-only record's rows under
/// `--regulatory`, Ensembl VEP's SPDI names the chromosome as the cache does,
/// since AnnotationType/RegFeat.pm:126 `$vf->{slice} ||= $slice;` hands the
/// variation feature the regulatory feature's slice before InputBuffer.pm:518
/// builds one from the FASTA, where vep-rs (`reference_name` in
/// `vep-effects/src/spdi.rs`) has no regulatory slice and names the
/// chromosome as the FASTA does.
#[test]
fn regulatory_spdi_chromosome_name_exclusion_is_documented() {
    assert_page_states(
        "divergences/regulatory-rows.md",
        &[
            "When the FASTA names a chromosome unlike the cache, a regulatory-only record's rows \
             under `--spdi` name it as the cache does in Ensembl VEP and as the FASTA does here, \
             `21:30234967:A:G` against `chr21:30234967:A:G`",
            "The regulatory feature hands the variant its slice there (`AnnotationType/RegFeat.pm` \
             line 126), a slice vep-rs never has",
        ],
    );
}

/// `fasta`: the three refusals at open the page's sentence counts each have a
/// table row, the `.gzi` one included (`checked_gzi` in `vep-fasta`, whose
/// `FastaError::GziMismatch` prints `BGZF index {gzi} does not describe {fasta}
/// ({reason}): index the file again with samtools faidx`).
#[test]
fn fasta_open_refusals_are_each_tabulated() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Three files are refused at open: a gzip file that is not BGZF, whose `.fai` offsets \
             need the whole file inflated; a `.gzi` that does not match the file's blocks, \
             recompressed and not re-indexed; an `.fai` naming more bases than the file holds",
            "| A gzip file that is not BGZF | `ERROR: [fasta] <file> is gzip but not BGZF: \
             recompress with bgzip and index with samtools faidx` |",
            "| A `.gzi` that does not describe the file's blocks | `ERROR: [fasta] BGZF index \
             <file>.gzi does not describe <file> (<the first disagreement>): index the file again \
             with samtools faidx` |",
            "| An `.fai` naming bases past the file's end | `ERROR: [fasta] FASTA index names <n> \
             bases beyond the end of <file>` |",
        ],
    );
}

/// `spdi`: a record whose expanded reference allele ends in a character outside
/// `ACGTN-` is annotated with no `SPDI` (`spdi_genomic` in `vep-core/src/spdi.rs`
/// returns `None` for it), where Ensembl VEP 116.2's `spdi_genomic` throws
/// (VariationFeature.pm:2159 `if( $ref_allele =~ m/[^ACGTN\-]$/ig ){`, :2160
/// `throw("No supported SPDI genomic is available for Variation Feature ...")`),
/// called with no `eval` from OutputFactory.pm:1398 under `--spdi` or
/// `--ga4gh_vrs` (:1397), so its run ends with exit 255 before any row of the
/// record's buffer is written (Runner.pm:391 fills the output buffer per input
/// buffer, :208-209 prints it); Parser.pm:556's unanchored allele test lets the
/// record through to it. The output page and the `--spdi` help state the reading.
#[test]
fn spdi_reference_allele_stop_is_documented() {
    assert_page_states(
        "output/vep-and-tab.md",
        &[
            "A record whose reference allele ends in a character outside `ACGTN-` after that \
             expansion (`AR/CA`, `(AC)/-`) is annotated with no `SPDI` on any row",
            "Ensembl VEP's `spdi_genomic` throws on such a record (`VariationFeature.pm` line \
             2159), under `--ga4gh_vrs` too, and the run ends with exit 255 before any row of the \
             record's buffer is written",
        ],
    );
    assert_flag_row_states(
        "cli/output-control.md",
        "--spdi",
        &[
            "A record whose reference allele ends in a character outside `ACGTN-` (`AR/CA`) is \
             annotated with no SPDI, where Ensembl VEP stops the run on it",
        ],
    );
}

/// `fasta`: a `--fasta` value Perl reads as false (`0`, the empty value) is a
/// path here (`fasta: args.fasta` in `config.rs`, with no truth filter), refused
/// at open, where Ensembl VEP reads the flag as unset, so CacheDir.pm:235
/// (`unless($self->param('fasta') || $self->param('no_fasta'))`) auto-detects the
/// cache directory's FASTA and BaseVEP.pm:521 (`if(my $fasta_file =
/// $self->param('fasta'))`) opens none otherwise; `Config.pm:246` declares it
/// `fasta|fa=s`. The compatibility page states the reading beside the other
/// false-value flags and the cache page's FASTA table cross-references it.
#[test]
fn fasta_false_value_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A `--fasta` value Perl reads as false, `0` or the empty value, leaves that flag unset \
             there likewise, so the cache directory's FASTA is auto-detected or the run proceeds \
             without one (`CacheDir.pm:235`, `BaseVEP.pm:521`)",
            "vep-rs opens the value as a path and refuses it with `Specified FASTA file/directory \
             <value> not found` (exit 1)",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &[
            "| `--fasta 0` or `--fasta ''`, a value Ensembl VEP reads as the flag unset \
             ([compatibility page](cli/compatibility.md#where-configpm-is-not-followed)) | \
             `ERROR: [fasta] Specified FASTA file/directory <value> not found` |",
        ],
    );
}
