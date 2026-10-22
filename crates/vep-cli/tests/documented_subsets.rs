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
/// contig a one-base placeholder slice reading `N` (FastaSequence.pm:597-600), so
/// Parser.pm:614 passes a one-base record at position 1 written `N` and :624-631
/// refuses any other with the mismatch line alone.
#[test]
fn check_ref_exclusion_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Under `--check_ref` the same contig skips the record at every position with the \
             sub-slice warning and then the mismatch warning",
            "Ensembl VEP annotates a one-base record there whose reference is written `N`, and \
             refuses one written otherwise with the mismatch warning alone (`Parser.pm` lines \
             614-631)",
        ],
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

/// `config`: a lower source never overwrites a higher source's value, and a
/// value the parser refuses is an argument syntax error, where Config.pm:1101's
/// `||=` lets a false value yield and :1046-1047 drops the line.
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
}

/// `config-file`: the same merge rule read from `vep.ini` and the environment
/// (Config.pm:1101, 1158), with the page's own `shift_hgvs` example.
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
}

/// `dont_skip`: four refused shapes are skipped under the flag, where
/// Parser.pm:189 returns a record whatever `validate_vf` returned.
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
/// so FastaSequence.pm:597-600 returns `N` for a one-base record at position 1
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
/// resolves as a synonym prefix (SliceAdaptor.pm:2871-2884), and a span past
/// the chromosome end, which it pads with `N` (FastaSequence.pm:619-639), are
/// skipped.
#[test]
fn spdi_input_format_exclusion_is_documented() {
    assert_page_states(
        "divergences/spdi-input.md",
        &[
            "Ensembl VEP reads SPDI notation (`--format spdi`) only with a database \
             connection",
            "vep-rs resolves the same notations from the cache and the reference FASTA",
            "queries the synonym table again with the name as a prefix",
            "vep-rs matches a name against the cache's chromosomes and their synonyms \
             exactly, regardless of ASCII case, never as a prefix",
            "tests no coordinate against the sequence's length",
            "so the record is annotated at a coordinate the chromosome lacks",
            "refuses a span the FASTA does not hold",
        ],
    );
}
