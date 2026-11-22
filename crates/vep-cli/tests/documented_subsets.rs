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
/// Utils.pm:467-468 reads `localtime`, the Perl API and database header lines
/// have no counterpart, and the `--custom` header lines follow the plugin lines
/// (OutputFactory/VEP_output.pm:167).
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
            "The `## <short name> : [PATH]/<file>` lines of `--custom` sources follow the \
             plugin lines, one per key of each source",
        ],
    );
}

/// `vcf-csq-info-header`: the `##VEP=` line carries no Perl API component
/// versions, the `--custom` `##INFO` lines follow the plugin lines
/// (OutputFactory/VCF.pm:234-263, 273-276), and a header without a `#CHROM`
/// line gets the eight-column line, where VCF.pm:224 pops the last `#` line of
/// the input whatever it is.
#[test]
fn vcf_csq_info_header_exclusion_is_documented() {
    assert_page_states(
        "output/vcf.md",
        &[
            "Ensembl VEP's `##VEP=` line also carries its Perl API component versions \
             (`ensembl=116.0d85231 ...`), which have no counterpart here",
            "Its `##INFO=<ID=<short>_<field>,...>` lines for `--custom` sources follow the \
             plugin lines, one per key of each source",
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

/// The text of the cache builder's README, whitespace collapsed.
fn cache_builder_readme() -> String {
    std::fs::read_to_string(format!(
        "{}/../vep-cache-builder/README.md",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
    .split_whitespace()
    .collect::<Vec<_>>()
    .join(" ")
}

/// `install-pl-auto`, `api-install-options`: `vep-cache-builder --auto` is the
/// cache and FASTA steps of INSTALL.pl's `--AUTO cf` (INSTALL.pl:296-331), the
/// cache built from the release files where INSTALL.pl:1402-1418 downloads the
/// prebuilt archive and checks its `sum` against the directory's CHECKSUMS, and
/// `--CACHEDIR`, `--CACHE_VERSION`, `--CACHEURL` and `--FASTAURL`
/// (INSTALL.pl:217-220, 273-280) have counterparts where the API, plugin and
/// build options have none.
#[test]
fn install_pl_auto_subset_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "`--auto` fetches every input from the Ensembl site, checks each against its \
             directory's `CHECKSUMS` file, builds the cache as \
             `<dir_cache>/<species>/<release>_<assembly>` and places the indexed genome FASTA \
             inside it",
            "A file whose BSD `sum` differs from its `CHECKSUMS` line stops the installation, \
             naming the file",
            "an existing release directory is never written into",
            "The cache is built from the release files, never from Ensembl VEP's prebuilt cache \
             archive, whose Storable format `vep` does not read",
            "`--include-predictions` adds the SIFT and PolyPhen matrices",
            "`--include-variations` the variation data below",
        ],
    );
    assert_states(
        "crates/vep-cache-builder/README.md",
        &cache_builder_readme(),
        &[
            "`--auto <species> <assembly>` is its `--AUTO cf --SPECIES --ASSEMBLY`, the cache \
             and FASTA steps alone",
            "`--dir_cache` is `--CACHEDIR`, `--release` is `--CACHE_VERSION` and `--site` \
             stands for both `--CACHEURL` and `--FASTAURL`",
            "`--include-predictions` and `--include-variations` add what that archive carries \
             by default",
            "`--GITHUB_TOKEN` and the `VEP_*` environment variables have no counterpart",
        ],
    );
}

/// `fasta-download-and-indexing`: the installer takes `dna_index/` before
/// `dna/` under every protocol where INSTALL.pl:1576 reads `dna_index/` over
/// FTP alone, reads GRCh37 from its own release tree where INSTALL.pl:1476
/// rewrites the URL to release 75, keeps a `dna/` file that is already BGZF
/// where INSTALL.pl:1676-1678 recompresses every one, prefers the primary
/// assembly to the toplevel file as INSTALL.pl:1611-1612 does, and writes the
/// bgzip-compressed file with `.fai` and `.gzi` only, never the plain `.fa` with
/// a `Bio::DB::Fasta` index of INSTALL.pl's `NO_HTSLIB` route (:1713-1719).
#[test]
fn fasta_download_and_indexing_subset_is_documented() {
    assert_page_states(
        "cache-setup.md",
        &[
            "the `dna_index/` file with its published `.fai` and `.gzi`, else the `dna/` file, \
             recompressed to bgzip and indexed by the builder",
            "the shape of Ensembl's `dna_index/` files and of the file `vep-cache-builder \
             --auto` places in the cache directory",
            "An index the file lacks is written beside it as `samtools faidx` writes it, as \
             Ensembl VEP writes it",
        ],
    );
    assert_states(
        "crates/vep-cache-builder/README.md",
        &cache_builder_readme(),
        &[
            "The FASTA step takes `dna_index/` before `dna/` whatever the protocol, keeps a \
             `dna/` file that is already BGZF, reads GRCh37 from `pub/grch37/release-<N>/` \
             rather than release 75, and writes a bgzip-compressed file with `.fai` and `.gzi` \
             only, never the plain `.fa` with a `Bio::DB::Fasta` index `INSTALL.pl` writes \
             without htslib",
        ],
    );
    assert_page_states(
        "getting-started/quick-start.md",
        &["The FASTA `--auto` places in the cache directory is found without a flag"],
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
            "So a VCF written from a VCF carries every input byte of the columns it stores",
            "A record is otherwise read as UTF-8 with each invalid sequence replaced by \
             U+FFFD, sample columns included",
            "prints the replacement character in `Uploaded_variation`, in the JSON `input` \
             string",
            "A `--regulatory_gff` or `--custom` GFF record's attribute values and type column \
             follow the rule too.",
            "so two identifiers or two types differing only in such bytes are one to vep-rs and \
             two to Ensembl VEP",
            "Both readers decode a type column's escapes, and an attribute value's, over the \
             record's bytes, so an escape completing a raw first byte of a character forms the \
             character as Ensembl VEP's bytes do.",
            "`vep-cache-builder` decodes the escapes the same way and refuses a line holding such \
             a byte raw, naming the line.",
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
/// (Config.pm:1101, 1158), with the page's own `shift_hgvs` example, the same
/// reading of a repeatable key alone on its line, a value's literal
/// `___SPACE___` kept where Config.pm:1082 substitutes that text for a blank
/// inside double quotes and :1089 restores every occurrence, the literal one
/// included, and a prefix or case variant of a key dropped where
/// `is_valid_param` (Config.pm:1035-1052, Getopt::Long's `ignore_case` and
/// `auto_abbrev`) passes it, :1101 stores it under the written key and
/// `full_command` (:734-745) echoes that key in the `## VEP command-line:`
/// header.
#[test]
fn config_file_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "so `--shift_hgvs 0` on its command line yields to a file's `shift_hgvs 1`",
            "A flag read from a file or the environment acts exactly as if typed, refusals \
             included",
            "A literal `___SPACE___` in a value stays as written",
            "turns a literal one into a blank too (`Config.pm:1082, 1089`)",
            "the `## VEP command-line:` header echoes the key",
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
            "The crate lacks back-references, look-around, atomic, branch-reset and conditional \
             groups, recursion, code blocks, backtracking control verbs, the `\\b{wb}` bound \
             types, the `(?l)` locale flag, the `\\X` escape and a `\\N{name}` character",
            "Under `--transcript_filter` the operator therefore passes a transcript through \
             the regex alone",
            "vep-rs matches the regex first and otherwise walks its built-in Sequence \
             Ontology table of consequence terms",
        ],
    );
}

/// `transcript_filter`, `transcript-filter-grammar` and `filter-vep`: an `in`
/// list is a hash (FilterSet.pm:391-392), and `get_input` picks one key of it
/// to stand for the filter's value before the numeric reading of a
/// `text:number` input (:597-599, a hash-order key), so a list mixing a number
/// and other text reads through the number under one seed (:600-601) and whole
/// under another (:603-604); vep-rs decides by the list's first written entry
/// (`read_list` and `unit` in `filter_set.rs`), a documented intended
/// divergence that reaches `vep filter` through the shared grammar.
#[test]
fn transcript_filter_in_list_numeric_by_hash_order_divergence_is_documented() {
    assert_page_states(
        "divergences/transcript-filter-in-list-numeric-by-hash-order.md",
        &[
            "Ensembl VEP 116.2 decides whether an `in` list compares a `text:number` or \
             `text(number)` value through its number by one entry it picks from the list's hash",
            "Lines 596-599: when the filter's value is a hash, `(keys %{$value})[0]` replaces it \
             by the first key the hash yields, under the comment `Assumes that all the values are \
             of similar format`",
            "vep-rs decides the reading by the list's first written entry, the same in every run",
            "Under `--transcript_filter 'gene_hgnc_id in HGNC:2559,1132'`, vep-rs writes the 2 \
             CXADR rows. Ensembl VEP 116.2 writes the 34 BTG3 rows under `PERL_HASH_SEED=0` and \
             the 2 CXADR rows under `PERL_HASH_SEED=11`",
            "Of the seeds 0 to 20, sixteen give the 34 BTG3 rows and five (11, 13, 15, 17 and 19) \
             the 2 CXADR rows.",
            "`filter_vep` keeps 34 rows under `PERL_HASH_SEED=0` and 2 under `PERL_HASH_SEED=11`, \
             with either written order; `vep filter` keeps 2 rows for that order and 34 for \
             `HGNC_ID in 1132,HGNC:2559`",
            "Writing the list as `2559,1132` or as `HGNC:2559,HGNC:1132` makes the engines agree \
             under every seed, 36 rows each",
            "The class is documented only: no committed flag set and no committed tool run passes \
             such a list.",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "so a list mixing a number and other text, `gene_hgnc_id in HGNC:2559,1132`, keeps \
             the BTG3 rows under one hash seed and the CXADR rows under another, in \
             `--transcript_filter` and `filter_vep` alike; vep-rs decides by the list's first \
             written entry, the same in every run | 3 · 0 of 392 flag sets · 9 flag sets, 0 \
             corpora | 116.2 |",
        ],
    );
    assert_page_states(
        "cli/filtering-picking.md",
        &[
            "An `in` list of `--transcript_filter` that mixes a number and other text, \
             `gene_hgnc_id in HGNC:2559,1132`, reads a `text:number` or `text(number)` value \
             through its number when its first written entry is a number and whole otherwise",
            "Ensembl VEP decides by an entry it picks from the list's hash, so write the entries \
             in one shape, `HGNC:2559,HGNC:1132` or `2559,1132`, and both engines keep the same \
             rows.",
        ],
    );
    assert_page_states(
        "cli/filter.md",
        &[
            "compares through its number against a numeric value, and against an `in` list when \
             the list's first written entry is a number, where `filter_vep` picks a hash-order \
             entry",
            "Write an `in` list's entries in one shape, `HGNC:2559,HGNC:1132` or `2559,1132`, and \
             both tools keep the same lines.",
        ],
    );
    assert_page_states(
        "divergences/transcript-filter-fields.md",
        &[
            "A third, in the grammar `--transcript_filter` and both tools share, is the \
             hash-order reading of an `in` list mixing a number and other text",
            "Write the list's entries in one shape, `HGNC:2559,HGNC:1132` or `2559,1132`, and the \
             engines keep the same rows.",
        ],
    );
}

/// `transcript_filter` and `filter-vep`: Ensembl VEP's parser fills a unit's
/// field slot while `!$current->{field}` holds (FilterSet.pm:260), a Perl truth
/// test false for the string `0`, so a unit whose first word is `0` loses it:
/// the operator word becomes the field, and a value is read as an operator and
/// refused (:276). vep-rs reads `0` as the field (`parse_into` in
/// `filter_set.rs`), which no transcript key matches, a documented intended
/// divergence that reaches `vep filter` through the shared grammar.
#[test]
fn transcript_filter_zero_field_divergence_is_documented() {
    assert_page_states(
        "divergences/transcript-filter-zero-field.md",
        &[
            "Ensembl VEP 116.2 tests whether a filter unit has its field yet with Perl's truth \
             test, and the string `0` is false to Perl. So a unit whose first word is `0` loses \
             that word",
            "Lines 260-262: the field slot is filled while `!$current->{field}` holds; `0` is \
             false to Perl, so the next word fills the slot again.",
            "vep-rs reads the first word of a unit as its field whatever the word is, so `0` \
             names a field called `0`, which no transcript key and no output column is named or \
             begins with, undefined everywhere.",
            "Under `--transcript_filter '0 is 1'`, Ensembl VEP 116.2 refuses the run with `ERROR: \
             No such operator \"1\"` and writes nothing; vep-rs runs and writes the record's one \
             `intergenic_variant` row",
            "Under `--transcript_filter '0 ne'`, Ensembl VEP writes that intergenic row and \
             vep-rs the 41 transcript rows.",
            "Under `0 ne`, Ensembl VEP writes the 35 intergenic rows and vep-rs the 952.",
            "`filter_vep -filter '0 is 1'` and `-filter '0 in <list>'` refuse with the same \
             messages",
            "The class is documented only: no committed flag set and no committed tool run passes \
             such a unit.",
            "A `--custom` source given `short_name=0` is named after its file on both engines",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "so a unit whose first word is `0` loses it: `0 is 1` reads `is` as the field and \
             refuses `1` as `No such operator`, and `0 ne` tests a field named `ne`; vep-rs reads \
             `0` as the field, which no transcript key matches, so `0 is 1` keeps nothing and \
             `0 ne` keeps every transcript | 3 · 0 of 392 flag sets · 9 flag sets, 0 corpora | \
             116.2 |",
        ],
    );
}

/// `transcript_filter` and `filter-vep`: FilterSet.pm:397 opens an `in` list's
/// file with Perl's two-argument `open`, which reads a mode out of the name
/// after :396 tested the written name with `-e`, so `<ids.txt` opens `ids.txt`
/// and `true|` runs the command `true` and reads its output; vep-rs opens the
/// file of the written name (`read_list` in `filter_set.rs`), a documented
/// intended divergence that reaches `vep filter` through the shared grammar.
#[test]
fn transcript_filter_list_file_two_argument_open_divergence_is_documented() {
    assert_page_states(
        "divergences/transcript-filter-list-file-two-argument-open.md",
        &[
            "Ensembl VEP 116.2 opens the file an `in` list names with Perl's two-argument `open`, \
             which reads a mode out of the name",
            "Line 397: `open IN, $list or throw(\"ERROR: Could not read from file $list\\n\")`, \
             the two-argument form",
            "vep-rs opens the file of the written name and reads its lines, whatever characters \
             the name holds.",
            "Under `--transcript_filter 'stable_id in <ids.txt'`, vep-rs keeps the one row of \
             `ENST00000624153` on the record `tf_7130369`; Ensembl VEP 116.2 reads `ids.txt` \
             instead and keeps the two rows of `ENST00001020056` on the record `tf_32493707`.",
            "`stable_id in true|` keeps that row on vep-rs, where Ensembl VEP runs the command \
             `true`, reads an empty list and drops every transcript.",
            "Under `<ids.txt`, Ensembl VEP writes 2 transcript rows and 33 intergenic, vep-rs 1 \
             and 34. Under `true|` and the file-writing name, Ensembl VEP writes 0 and 35, vep-rs \
             1 and 34.",
            "The mode reading is never replicated: a name beginning with `<` or `>`, or ending \
             with `|`, is a file name and nothing else",
            "The class is documented only: no committed flag set and no committed tool run names \
             such a file.",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "so a file named `<ids.txt` is tested for and `ids.txt` is read, a name beginning \
             with `>` truncates another file, and a file whose name ends in a vertical bar is \
             tested for and the rest of its name is run as a command, its output read as the \
             list; vep-rs opens the file of the written name and never replicates the mode \
             reading | 3 · 0 of 392 flag sets · 9 flag sets, 0 corpora | 116.2 |",
        ],
    );
}

/// `transcript_filter` and `filter-vep`: FilterSet.pm:396 finds a directory
/// with `-e`, :397's `open` succeeds on it, and :398's `readline` fails with
/// `EISDIR`, read as the end of the file, so the list is empty and every
/// transcript is dropped with no message; vep-rs refuses the run with the
/// message Ensembl VEP prints for an open that fails (`read_list` in
/// `filter_set.rs`), a documented intended divergence that reaches `vep filter`
/// through the shared grammar.
#[test]
fn transcript_filter_list_directory_divergence_is_documented() {
    assert_page_states(
        "divergences/transcript-filter-list-directory.md",
        &[
            "Ensembl VEP 116.2 reads a directory named as an `in` list file as an empty list: the \
             open succeeds, every read of the handle fails, the loop ends with nothing read",
            "Lines 398-410: `while(<IN>)` reads the handle; `read(2)` on a directory fails with \
             `EISDIR`, `readline` returns nothing, the body never runs and the list stays empty.",
            "vep-rs refuses the run with `ERROR: [filter] Could not read from file <name>`, the \
             message Ensembl VEP prints when its own open fails, exit status 1 and no output. The \
             empty reading is never replicated",
            "Ensembl VEP 116.2 writes the record `tf_7130369`, within reach of 3 transcripts, as \
             one `intergenic_variant` row and prints nothing on its error stream",
            "Over the corpus's 32 records, Ensembl VEP writes 35 intergenic rows and none of the \
             952 transcript rows.",
            "The class is documented only: no committed flag set and no committed tool run names \
             a directory as a list.",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "the open succeeding and every read of the handle failing, so every transcript is \
             dropped and each record is written as intergenic with no message; vep-rs refuses the \
             run as a list it could not read | 3 · 0 of 392 flag sets · 9 flag sets, 0 corpora | \
             116.2 |",
        ],
    );
    assert_page_states(
        "divergences/transcript-filter-fields.md",
        &[
            "Three more lie in the reading of a unit and of an `in` list file. A unit whose first \
             word is `0`, which Ensembl VEP drops, is on",
            "A list file whose name carries a mode character is on",
            "and a list that names a directory on",
        ],
    );
}

/// `pick-order-criteria` and `tsl`: the `tsl` category sorts two levels with
/// Perl's `<=>` (OutputFactory.pm:789), which holds a digit run past 64 bits as
/// a double and rounds a 64-bit integer within one step of it to the same
/// double, so the two tie and the stable sort keeps the row listed first;
/// vep-rs orders the levels by exact value (`Rank` in vep-core's `pick`), a
/// documented intended divergence the gene-name pick of an HGVS input shares
/// (VariationFeatureAdaptor.pm:2207).
#[test]
fn tsl_pick_order_past_sixty_four_bits_divergence_is_documented() {
    assert_page_states(
        "divergences/tsl-pick-order-past-sixty-four-bits.md",
        &[
            "Ensembl VEP 116.2's pick order compares two transcript support levels with Perl's \
             `<=>`, which holds a level past 64 bits as a double and rounds a 64-bit integer \
             within one step of it to the same double",
            "vep-rs orders two levels by their exact values: two levels within 64 bits as \
             integers, a level past 64 bits after every integer, and two levels past 64 bits as \
             doubles",
            "Line 789: the candidates are sorted by `$a->{$cat} <=> $b->{$cat}`, a numeric \
             comparison of those strings",
            "Under `--pick --pick_order tsl`, Ensembl VEP 116.2 keeps ENST00001110244, the row \
             listed first; vep-rs keeps ENST00001113783, the lower level.",
            "The class is documented only: no comparator codes it, since no committed flag set \
             reaches it.",
            "Of the 392 committed flag sets, 14 pass a pick-family flag with `tsl` in the order \
             and none passes such a level; no committed HGVS input names a gene.",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "so two such levels tie and its stable sort keeps the row listed first; vep-rs orders \
             the levels by exact value and keeps the lower | 3 · 0 of 392 flag sets · 14 flag \
             sets, 0 corpora | 116.2 |",
        ],
    );
    assert_page_states(
        "output/vep-and-tab.md",
        &[
            "`TSL` carries the attribute's digit run as written, as Ensembl VEP prints it: `1` to \
             `5` from an Ensembl gene set, and another shape, `007` or `999`, only from a \
             third-party GFF3 build or an edited cache; two levels past 64 bits under the pick \
             order are [`tsl_pick_order_past_sixty_four_bits`]",
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
/// `skipped-variant-warning-lines`: the header names standard input `STDIN`
/// under `--input_data`, the record column is 40 characters, and the line a
/// skip warning quotes is cut at a character boundary, where
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
            "The header names the input as `-i` spelled it, `STDIN` by default, and under \
             `--input_data` it names `STDIN`, or the `--input_file` value when one stands \
             beside the string",
            "The record column is the line's first 40 characters, padded to 40 characters, so \
             a line with a multi-byte character pads to a different width",
            "On a line of ASCII characters under `-i <path>`, `-i STDIN` in any letter case or \
             no `-i`, the two files are byte-identical",
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

/// `option-prefix-abbreviation`: a prefix is resolved against this version's
/// flags first, then as Getopt::Long resolves it against `@VEP_PARAMS`
/// (Config.pm:73-276); a prefix of several options, a token `pass_through`
/// (Config.pm:62) leaves in `@ARGV` for vep:30 to refuse, is an argument
/// syntax error here.
#[test]
fn option_prefix_abbreviation_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "A prefix is resolved against this version's flags first",
            "A prefix that names one flag of this version resolves to it (`--bio` is `--biotype`)",
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
/// appends every key it is handed; a plugin key equal to an active flag's column
/// (a dbNSFP column named `TSL` under `--tsl`) is one column too, where Tab.pm:194 lists the
/// flag's columns and :195 the plugin's.
#[test]
fn tab_column_order_exclusion_is_documented() {
    assert_page_states(
        "output/vep-and-tab.md",
        &[
            "A plugin key listed twice, by a plugin named twice or two plugins sharing a key, \
             prints one column where Ensembl VEP prints one per instance \
             ([a key listed twice](../divergences/duplicate-output-column.md))",
            "a key named as a flag field (a dbNSFP column named `TSL`) takes the flag's place under its flag, \
             once where Ensembl VEP [lists it twice](../divergences/duplicate-output-column.md)",
        ],
    );
}

/// `output_format`: OutputFactory.pm:229 lower-cases a copy of the value to pick
/// the writer while :150 keeps the spelling as typed, which :1527 (a symbol's
/// `%20` under exactly `vcf`) and :1409 (the GA4GH_VRS string unless exactly
/// `json`) read case-sensitively; vep-rs reads the value once, lower-cased.
#[test]
fn output_format_given_case_divergence_is_documented() {
    assert_page_states(
        "divergences/output-format-flag-case.md",
        &[
            "Ensembl VEP 116.2 lower-cases the `--output_format` value to choose the writer and \
             keeps the spelling as typed for two later readers",
            "`output_format_flag_read_in_its_given_case`: a `SYMBOL` in a VCF output selected \
             with `--output_format` spelled other than `vcf`, whose whitespace Ensembl VEP writes \
             as `_` and vep-rs as `%20`. It also covers a `GA4GH_VRS` in a JSON output selected \
             with `--output_format` spelled other than `json`, which Ensembl VEP writes as a \
             string and vep-rs as an object",
            "Spelled `vcf`, or given as `--vcf`, both engines write `PGR5-LIKE%20A`; under \
             `--no_escape` both write `PGR5-LIKE_A`",
            "vep-rs writes the object under both spellings, equal to Ensembl VEP's `json` object \
             member for member",
        ],
    );
    assert_page_states(
        "output-formats.md",
        &[
            "The `--output_format` value is matched in any case (`VCF`, `Json`) and read once, \
             lower-cased, by every part of the run; Ensembl VEP lower-cases it to choose the \
             writer alone and reads the spelling as typed in two places, a documented intended \
             divergence, [`output_format_flag_read_in_its_given_case`]\
             (divergences/output-format-flag-case.md)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP lower-cases the `--output_format` value to choose the writer and reads \
             the spelling as typed in two later tests, so under `--output_format VCF` a gene \
             symbol's whitespace is written as `_` in place of `%20`, and under \
             `--output_format JSON` the `GA4GH_VRS` object is written as a JSON string; vep-rs \
             reads the value once, lower-cased, so every spelling writes the lower-case \
             spelling's rows",
        ],
    );
}

/// `AVI`, `OpenTargets`: OutputFactory.pm:2291-2293 restore the record's span and
/// allele string on the first split feature and :2314 points every later allele's
/// overlaps at it before :2226 runs the plugins, so AVI.pm:130-134 reads the
/// record's span (no score on a split record) and OpenTargets.pm:349 the record's
/// alleles (a match on every allele's rows); vep-rs hands each split allele its own.
#[test]
fn minimal_plugin_rejoined_record_divergence_is_documented() {
    assert_page_states(
        "divergences/minimal-plugin-sees-the-rejoined-record.md",
        &[
            "Under `--minimal`, Ensembl VEP's output pass hands every plugin the rejoined record, \
             its span and allele string restored, so a plugin that reads the record's span or \
             alleles values every split allele by the record and not by the allele's own trimmed \
             form",
            "vep-rs writes each allele's row with its score and link; Ensembl VEP 116.2 writes \
             both rows without an AVI field, since the rejoined record spans two bases",
            "Under `--minimal` Ensembl VEP 116.2 writes the association's five fields on all six \
             rows of the record; vep-rs writes them on the first allele's three rows alone",
            "gnomADc, gnomADMt, Conservation, NearestExonJB, TSSDistance and GWAS read the same \
             feature, and vep-rs hands them the allele's own span and alleles too",
            "One committed flag set runs a covered plugin under `--minimal`, gnomADc over the \
             coverage dataset",
        ],
    );
    assert_page_states(
        "plugins/annotation-store.md",
        &[
            "under `--minimal` a split allele over its own trimmed span, where Ensembl VEP reads \
             the rejoined record and scores no allele of it",
            "under `--minimal` a split allele by its own trimmed form, where Ensembl VEP reads the \
             rejoined record and writes a matched record's fields on every allele's rows",
            "Under `--minimal`, a split allele reads its own trimmed span, the one its row's \
             alleles describe; Ensembl VEP's rejoin hands the plugin the record's span",
        ],
    );
    assert_page_states(
        "plugins/conservation.md",
        &[
            "Under `--minimal`, each allele of a multi-allelic record the flag splits reads \
             its own trimmed span, the one its row's alleles describe, as it does when it is \
             the record's only alternate",
        ],
    );
    assert_page_states(
        "plugins/transcript-model.md",
        &[
            "Under `--minimal`, a split allele reads its own trimmed span for its exon set \
             and distances; Ensembl VEP's rejoin hands the plugin the record's span",
            "Under `--minimal`, a split allele measures from its own trimmed span, the one its \
             row's alleles describe; Ensembl VEP's rejoin hands the plugin the record's span",
        ],
    );
    assert_page_states(
        "plugins/flat-file.md",
        &[
            "Under `--minimal` a split allele is matched by its own trimmed alleles over its own \
             span; Ensembl VEP matches the rejoined record's whole allele list over its span",
        ],
    );
    assert_page_states(
        "divergences/nearest-exon-jb-shared-span-cache.md",
        &[
            "a `--minimal` allele's exon set and distances from its trimmed span",
            "vep-rs writes `ENSE00001282176|0|end|272` on the `G` row and nothing on the `T` row, \
             in both orders",
        ],
    );
    assert_page_states(
        "plugins/ancestral-allele.md",
        &[
            "a plugin that reads the rejoined record's span or alleles without caching is the \
             sibling class [`minimal_plugin_sees_the_rejoined_record`]\
             (../divergences/minimal-plugin-sees-the-rejoined-record.md)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "AVI scores no allele of a split record and OpenTargets writes a matched record's \
             fields on every allele's rows; vep-rs hands each split allele its own span and \
             alleles",
        ],
    );
}

/// `fork`: Runner.pm:662 captures a child's standard error into a string and :565
/// relays it through `warning_msg` behind the child's pid, so BaseVEP.pm:806
/// writes `WARNING: <pid> : WARNING: <text>`; vep-rs writes `WARNING: <text>` at
/// every worker count, and the warnings comparator sets the prefix aside.
#[test]
fn fork_warning_pid_prefix_divergence_is_documented() {
    assert_page_states(
        "divergences/fork-warning-pid-prefix.md",
        &[
            "Under `--fork`, Ensembl VEP's parent process relays a child's warnings into the \
             warnings file as one block headed by the child's process id, so the block's first \
             line reads `WARNING: <pid> : WARNING: <text>` where the same run without the flag \
             writes `WARNING: <text>`",
            "vep-rs writes every warning as the condition alone, `WARNING: <text>`, at every \
             worker count, the shape Ensembl VEP writes without `--fork`",
            "In the format comparator's warnings mode the class, named in `--masks`, removes a \
             leading `<pid> : ` from each of the reference's lines before the two files are \
             compared and counts the lines it changed under `masks_applied`",
            "The golden datasets record each reference warning in their manifests as the \
             condition alone, so the harness compares vep-rs's lines against the serial shape",
        ],
    );
    assert_page_states(
        "output/refseq-and-merged.md",
        &[
            "Ensembl VEP under `--fork` writes the line behind the child process's id, a \
             documented intended divergence, [`fork_warning_pid_prefix`]\
             (../divergences/fork-warning-pid-prefix.md)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Under `--fork`, Ensembl VEP's parent process relays a child's warnings into the \
             warnings file as one block headed by the child's process id, so the block's first \
             line reads `WARNING: <pid> : WARNING: <text>` where a run without the flag writes \
             `WARNING: <text>`; vep-rs writes the condition alone at every worker count",
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
/// shape, not the flag, as the trigger. The same comparison runs the other way
/// on a deletion whose chopped alternate the second trim leaves as the
/// reference (`ATT/AT`, chopped `TT/T`, minimised `T/-`; `ATTGG/ATG`, chopped
/// `TTGG/TG`, minimised `TG/-`): :720 keeps no allele of a `1/1` sample and
/// Ensembl VEP reads it as `HOMREF` on the reference row (:723-729,
/// OutputFactory.pm:506-507), where vep-rs reads a homozygous carrier, `HOM` on
/// the deletion's row.
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
            "against the minimised reference (`Parser/VCF.pm` lines 603-639, after \
             `post_process_vfs` at lines 362-364), which the chopped reference never equals",
            "When the second trim leaves a reference equal to the chopped alternate, that \
             alternate counts as no alternate. It does for a deletion of one copy of a doubled \
             unit such as `ATT/AT` (chopped `TT/T`, trimmed to `T/-` one base on) and for \
             `ATTGG/ATG` (chopped `TTGG/TG`, trimmed to `TG/-`). Ensembl VEP then reads the \
             homozygous carrier (`1/1`) as reference-homozygous, a copy with no row of its own \
             unless `--process_ref_homs` adds its reference row with `ZYG` `HOMREF`",
            "vep-rs compares the chopped alleles with the chopped reference, so the `0/0` \
             individual is reference-homozygous, the `1/1` individual a homozygous carrier",
            "as does a same-length record without `--minimal`",
            "For such a record Ensembl VEP prints `HOM` for a `0/0` sample where vep-rs prints \
             `HOMREF`, and, where the chopped alternate equals the minimised reference, \
             `HOMREF` for a `1/1` sample where vep-rs prints `HOM`. With that sample alone \
             named, Ensembl VEP keeps the alternate allele's rows for the `0/0` sample and \
             annotates the reference allele alone for the `1/1` sample; vep-rs does the \
             reverse",
            "So `21 34716539 GA GT` with one sample `0/0` under `--individual_zyg --minimal` is \
             the `T` row with `ZYG=S1:HOM` and `ALLELE_NUM=1` in Ensembl VEP and the `A` row \
             with `ZYG=S1:HOMREF` and `ALLELE_NUM=0` here, `MINIMISED=1` on both",
            "And `21 34716420 ATT AT` with one sample `1/1` under `--individual_zyg \
             --allele_number`, with or without `--minimal`, is the `T` row with `ZYG=S1:HOMREF` \
             and `ALLELE_NUM=0` in Ensembl VEP and the `-` row with `ZYG=S1:HOM` and \
             `ALLELE_NUM=1` here",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP builds an individual's copy of a bi-allelic record it trims past the \
             anchor chop, an indel such as `ACAC/AC` with or without `--minimal` or a \
             same-length record such as `GA/GT` under it, from the chopped alleles against the \
             minimised reference",
            "and a `1/1` genotype of a deletion whose chopped alternate equals that reference, \
             such as `ATT/AT`, reads `HOMREF`; vep-rs compares chopped with chopped, keeps the \
             minimised alleles and reads `0/0` as reference-homozygous and `1/1` as a \
             homozygous carrier",
            "[`individual_trimmable_indel`](divergences/individual-trimmable-indel.md)",
            "records with sample columns",
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
             vep-rs reads `HOMREF` and annotates the reference allele. A `1/1` sample of a \
             deletion whose chopped alternate equals the minimised reference, such as \
             `ATT/AT`, Ensembl VEP reads as `HOMREF`, annotating the reference allele, where \
             vep-rs reads `HOM` and keeps the deletion's rows",
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

/// `info-var-type-serialiser`: a cache whose `var_type` contradicts the
/// variation files a chromosome directory holds, or whose `variation_cols`
/// names columns but no `start` or does not begin with `chr` under `var_type
/// tabix`, is refused as the variation source opens (`VariationSource::open` in
/// `variation_source.rs`, exit 1), where Ensembl VEP selects its reader by the
/// key (CacheDir.pm:184-186), opens the one file name that reader knows only
/// `if -e $file` (VariationTabix.pm:156 and 270, Variation.pm:104) and zips
/// each record's fields with the column list as written
/// (BaseCacheVariation.pm:119-124), so it writes every record on such a
/// chromosome without a co-located variant, exit 0.
#[test]
fn var_type_layout_refusal_is_documented() {
    assert_page_states(
        "divergences/damaged-cache-shard.md",
        &[
            "each reader opens the one file name it knows only if the file exists",
            "refuses a cache whose `var_type` and files disagree, in either direction",
            "ERROR: [cache] <cache>/variations/21 holds region files (<start>-<end>_var.gz) \
             while <cache>/info.json says var_type tabix",
            "ERROR: [cache] <cache>/variations/21 holds all_vars.gz while <cache>/info.json \
             says var_type (unset)",
            "a `var_type tabix` cache whose `variation_cols` does not begin with `chr` is \
             refused as the source opens, exit 1",
            "ERROR: [cache] variation_cols of <cache>/info.json does not begin with chr, the \
             first column of an all_vars.gz record",
            "The list must also name `start`, the column each known variant is placed by: a list \
             naming other columns but not `start` is refused before the `chr` and file checks, \
             exit 1",
            "ERROR: [cache] variation_cols of <cache>/info.json names no start column",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Under `--check_existing`, Ensembl VEP reads a chromosome's known variants from the \
             one file name `var_type` selects and, where a directory holds only the other form, \
             writes every record there with `Existing_variation` empty, exit 0; vep-rs refuses \
             the cache before any record is annotated, `ERROR: [cache] <dir>/variations/<chr> \
             holds <form> while <dir>/info.json says var_type <value>`, exit 1",
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
            "[`contig_absent_from_fasta`](cache-setup.md#reference-fasta)",
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
            "Two hundred and forty-eight differences lie outside every class",
            "Each of the two hundred and forty-eight [differences outside the classes]",
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

/// The regulatory-rows page names the committed runs that carry the
/// `regulatory_feature_rows` mask class, and the manifests decide which they
/// are: each flag set naming the class is a `<corpus>` whose `<flag set>` the
/// page's sentence must name, and a corpus the page names must carry it.
#[test]
fn regulatory_feature_rows_carriers_are_the_manifests() {
    let mut carriers: Vec<(String, String)> = Vec::new();
    for corpus in common::corpora() {
        let names_class = |masks: &serde_json::Value| {
            masks
                .as_array()
                .is_some_and(|m| m.iter().any(|v| v == "regulatory_feature_rows"))
        };
        if names_class(&corpus.manifest["masks"]) {
            carriers.push((corpus.name.clone(), "default".to_string()));
        }
        for flag_set in corpus.manifest["flag_sets"]
            .as_array()
            .into_iter()
            .flatten()
        {
            if names_class(&flag_set["masks"]) {
                let name = flag_set["name"].as_str().unwrap_or("").to_string();
                carriers.push((corpus.name.clone(), name));
            }
        }
    }
    carriers.sort();
    assert!(
        !carriers.is_empty(),
        "no committed run carries regulatory_feature_rows"
    );
    let text = page_text("divergences/regulatory-rows.md");
    let paragraph = text
        .lines()
        .find(|l| l.contains("Of the committed corpora, "))
        .expect("docs/divergences/regulatory-rows.md names the runs carrying the class");
    let start = paragraph.find("Of the committed corpora, ").unwrap();
    let sentence = &paragraph[start..];
    let sentence = &sentence[..sentence.find(". ").map_or(sentence.len(), |end| end + 1)];
    for (corpus, flag_set) in &carriers {
        assert!(
            sentence.contains(&format!("`{corpus}`")),
            "{sentence} lacks `{corpus}`"
        );
        assert!(
            sentence.contains(&format!("`{flag_set}` flag set")),
            "{sentence} lacks `{flag_set}` flag set"
        );
    }
    let named: Vec<&str> = sentence
        .split('`')
        .skip(1)
        .step_by(2)
        .filter(|name| name.contains('-'))
        .collect();
    let carrying: std::collections::BTreeSet<&str> =
        carriers.iter().map(|(c, _)| c.as_str()).collect();
    for corpus in named {
        assert!(
            carrying.contains(corpus),
            "{sentence} names `{corpus}`, which carries no such mask"
        );
    }
    let count = format!(
        "Of the committed corpora, {} carry the class",
        number_word(carrying.len())
    );
    assert!(sentence.contains(&count), "{sentence} lacks {count:?}");
}

/// The count word the page writes for `n` corpora.
fn number_word(n: usize) -> &'static str {
    [
        "none", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
    ][n]
}

/// `fasta`: the three refusals at open the page's sentence counts each have a
/// table row, the `.gzi` one included (`checked_gzi` in `vep-fasta`, whose
/// `FastaError::GziMismatch` prints `BGZF index {gzi} does not describe {fasta}
/// ({reason}): index the file again with samtools faidx`), and so has the
/// malformed block the page names after them (`FastaError::MalformedBgzfBlock`,
/// `malformed BGZF block at offset {offset} in {path}`).
#[test]
fn fasta_open_refusals_are_each_tabulated() {
    assert_page_states(
        "cache-setup.md",
        &[
            "Three files are refused at open: a gzip file that is not BGZF, whose `.fai` offsets \
             need the whole file inflated; a `.gzi` that does not match the file's blocks, \
             recompressed and not re-indexed; an `.fai` naming more bases than the file holds",
            "So is a BGZF file holding a malformed block.",
            "| A gzip file that is not BGZF | `ERROR: [fasta] <file> is gzip but not BGZF: \
             recompress with bgzip and index with samtools faidx` |",
            "| A `.gzi` that does not describe the file's blocks | `ERROR: [fasta] BGZF index \
             <file>.gzi does not describe <file> (<the first disagreement>): index the file again \
             with samtools faidx` |",
            "| An `.fai` naming bases past the file's end | `ERROR: [fasta] FASTA index names <n> \
             bases beyond the end of <file>` |",
            "| A BGZF block declaring more than 65,536 uncompressed bytes, the most a block holds | \
             `ERROR: [fasta] malformed BGZF block at offset <n> in <file>`, `<n>` the block's byte \
             offset in the file |",
        ],
    );
}

/// `fasta`: what a run does beside a FASTA given without an index has a table
/// row, and so have two refusals of that writing: the lock file that cannot be
/// written (`write_lock` in `vep-fasta`, whose `FastaError::LockFile` prints
/// Ensembl VEP's sentence, FastaSequence.pm:334, then the cause and the two
/// ways on) and the file no index describes (`FastaError::Unindexable`,
/// `failed to index {path}: {reason}`).
#[test]
fn fasta_index_writing_is_tabulated() {
    assert_page_states(
        "cache-setup.md",
        &[
            "| No `.fai`, or no `.gzi` beside a BGZF file | Writes the missing index into \
             `<index>.<token>.tmp` and renames it into place, so no run reads part of an index; \
             the directory must be writable |",
            "| None of `<file>.fai`, `<file>.index` and `<file>.gzi` | Writes `<file>.vep.lock` \
             before the indexes and removes it once they are written or the file is refused |",
            "| A `<file>.vep.lock`, left by a run that did not finish | Removes `<file>.fai`, \
             `<file>.index`, `<file>.gzi` and the lock, then writes the missing `.fai` and \
             `.gzi`; what it cannot remove stays, and an index that stays is read |",
            "| Other runs indexing the same file | Waits for none of them: each run writes whole \
             indexes, and each reads the file |",
            "| A lock file that cannot be written, as in a directory closed to the run | `ERROR: \
             [fasta] Could not write to FASTA lock file <file>.vep.lock: <cause>; index the FASTA \
             with samtools faidx, or link it into a directory the run can write` |",
            "| A layout no `.fai` describes: a line longer than the first of its sequence, text \
             after a blank or a shorter line, text before the first header, a last header with \
             no sequence line | `ERROR: [fasta] failed to index <file>: <reason>`, the reason \
             naming the line or the sequence |",
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

/// `api-feature-types-gating`: a shared library's `vep_plugin_feature_types` is
/// its one interest check, where OutputFactory.pm:2220 also asks
/// `check_variant_feature_type`; a library exporting no list runs on transcript
/// rows alone, as BaseVepPlugin.pm:99 defaults `feature_types` to `Transcript`.
#[test]
fn feature_types_gating_exclusion_is_documented() {
    assert_page_states(
        "developers/writing-a-plugin.md",
        &[
            "The feature type is the one interest check a shared library has. Ensembl VEP also \
             asks each plugin whether it takes the variant's kind, `VariationFeature` or \
             `StructuralVariationFeature` (`check_variant_feature_type`, `OutputFactory.pm` line \
             2220). The ABI has no such declaration, so a shared library sees sequence and \
             structural variants alike",
            "A library that does not export the symbol runs on transcript rows alone, Ensembl \
             VEP's default for a plugin that does not override `feature_types` \
             (`BaseVepPlugin.pm` line 99)",
        ],
    );
}

/// `api-filter-the-line`: a `run` that returns NULL or a non-object costs the
/// row that plugin's fields and a warning, where OutputFactory.pm:2238-2240
/// reads `undef` as a request to drop the line; :2234-2236 and :2243-2245 warn
/// and keep the line on a non-hashref or a die, as the host does.
#[test]
fn plugin_line_filter_exclusion_is_documented() {
    assert_page_states(
        "developers/writing-a-plugin.md",
        &[
            "A `run` call that returns NULL or text that is not JSON costs the row that plugin's \
             fields alone and writes the warning `Plugin '<name>' went wrong: <reason>`. A JSON \
             value that is not an object writes `Plugin '<name>' did not return a hashref, \
             output ignored!` and is dropped the same way. The row is kept, the plugins after it \
             still run on it, and the run goes on, as Ensembl VEP keeps the line when a plugin \
             returns something other than a hash or dies (`OutputFactory.pm` lines 2226-2245)",
            "Ensembl VEP reads a plugin's `undef` return as a request to drop the line from the \
             output (`OutputFactory.pm` lines 2238-2240, the `BaseVepFilterPlugin` mechanism). \
             The ABI has no filter return, so a shared library cannot remove a row",
        ],
    );
}

/// `api-filter-the-line`, the overwrite half: a plugin key naming a column or
/// field replaces its value in every format, as OutputFactory.pm:2231 writes
/// every returned key into the line hash after the core fields; a Parquet
/// column keeps its named type.
#[test]
fn plugin_core_field_overwrite_is_documented() {
    assert_page_states(
        "developers/writing-a-plugin.md",
        &[
            "A key a plugin returns that names a column or field of the output (`Codons`, \
             `TSL`, `SYMBOL`) replaces that value on the row in every format, as Ensembl \
             VEP's `$line_hash->{$key} = $plugin_results->{$key}` (`OutputFactory.pm` line \
             2231) does after the core fields set theirs. A declared flag field the plugin \
             did not write keeps the row's own value under its flag and is empty without it.",
            "A Parquet column keeps the type its name gives it \
             ([Parquet output](../output/parquet.md)), so a value outside that type, a list \
             in the INT32 `TSL` column, is stored NULL",
        ],
    );
}

/// `api-config-access`: a shared library receives its parameters and nothing of
/// the run's configuration, where BaseVepPlugin.pm:155-168 hands every plugin the
/// configuration hash.
#[test]
fn plugin_config_access_exclusion_is_documented() {
    assert_page_states(
        "developers/writing-a-plugin.md",
        &[
            "A plugin receives its `--plugin` parameters and nothing else of the run. Ensembl VEP \
             hands every plugin the whole configuration hash as well (`BaseVepPlugin.pm` \
             `config`, lines 155-168), which a plugin reads to learn the assembly, the cache \
             directory or any other flag. The ABI has no such call, so a shared library that \
             needs a run setting takes it as a parameter",
        ],
    );
}

/// `api-finish-has-cache-freeze`: the ABI has no `finish` call, where
/// Runner.pm:168-173 calls each plugin's `finish` after the last variant.
#[test]
fn plugin_finish_hook_absence_is_documented() {
    assert_page_states(
        "developers/writing-a-plugin.md",
        &[
            "The ABI has no `finish` call. Ensembl VEP calls each plugin's `finish` method, when \
             it has one, after the last variant (`Runner.pm` lines 168-173). A shared library \
             that holds a file to flush does so in `vep_plugin_destroy`, which the host calls on \
             every loaded instance when the run ends",
        ],
    );
}

/// `api-per-line-hook`: a shared library receives the consequence and variant
/// documents, never the rendered line, where BaseVepPlugin.pm:329-370 passes
/// `run` the `$line_hash` of the whole output line.
#[test]
fn plugin_per_line_hook_absence_is_documented() {
    assert_page_states(
        "developers/writing-a-plugin.md",
        &[
            "Ensembl VEP's per-line hook, the `$line_hash` of the whole output line that `run` \
             receives as its second argument (`BaseVepPlugin.pm` lines 329-370), has no \
             counterpart either. A shared library reads the consequence and variant documents \
             above, which carry the annotation the line is rendered from, and writes fields, \
             never the rendered line",
        ],
    );
}

/// `convert-cache-pl`: Ensembl VEP's `convert_cache.pl` (release/116.2, lines
/// 70-87 and 375-436) rewrites its own cache's `<start>-<end>_var.gz` region
/// files as one bgzipped, tabix-indexed `all_vars.gz` per chromosome and, under
/// `--sereal`, its transcript and regulatory shards as Sereal. The cache
/// conversion page states that it has no counterpart and needs none, the
/// converter copying the variation files in either form, and that
/// `vep-cache-converter` is another tool.
#[test]
fn convert_cache_pl_is_stated_not_applicable() {
    assert_page_states(
        "cache-conversion.md",
        &[
            "Ensembl VEP's `convert_cache.pl` rewrites its own cache's region variation files as \
             one tabix-indexed `all_vars.gz` per chromosome and, with `--sereal`, its transcript \
             shards as Sereal.",
            "It has no counterpart and none is needed: the converter copies a cache's variation \
             files in either form.",
            "`vep-cache-converter` is a different tool, which copies a flat JSON transcript cache \
             into the per-chromosome runtime layout.",
        ],
    );
}

/// `library-embedding`: the contract page states the surface (the `pub use`
/// list of `lib.rs`), the git-dependency recipe, the proof and the version
/// rule, so a rewrite that drops one of them fails here; the module doc of
/// `lib.rs` names the page.
#[test]
fn embedding_contract_is_documented() {
    assert_page_states(
        "developers/embedding.md",
        &[
            "The release workflow publishes the workspace's crates to crates.io when a \
             maintainer publishes a release page",
            "A program adds `vep-cli` from crates.io at a published version, or as a git \
             dependency pinned to a release's tag",
            "[dependencies.vep-cli] git = \"https://github.com/natera-open-source/vep-rs\" \
             tag = \"vX.Y.Z\"",
            "The surface is the list of `pub use` re-exports at the top of \
             `crates/vep-cli/src/lib.rs`",
            "everything else in `vep_cli` is a module of the binary, reachable at its module path \
             but not bound by this page",
            "The library has no skip step",
            "the binary writes no JSON object for a record whose span exceeds `--max_sv_size`",
            "`crates/vep-cli/tests/library_embedding.rs` runs this path over every record of the \
             `tests/golden/116/GRCh37` dataset",
            "the two are equal byte for byte",
            "While the version is 0.y.z, a release that changes any crate's public API, this \
             surface included, raises the minor number; a patch release keeps every public item's \
             name and signature",
        ],
    );
    let module_doc_names_the_page = include_str!("../src/lib.rs")
        .lines()
        .any(|line| line.starts_with("//!") && line.contains("`docs/developers/embedding.md`"));
    assert!(
        module_doc_names_the_page,
        "the lib.rs module doc names the page"
    );
}

/// `source-cache`, `info-bam` and `bam-edit-status-rna-edit`: the conversion
/// page states what a converted RefSeq or merged cache carries and what a run on
/// it does (the `source_cache` key where the Storable has it, the `bam` line
/// that switches `--use_transcript_ref` and `BAM_EDIT` on with Ensembl VEP's
/// notice, an edited model's cDNA and folded mapper, a failed model dumped
/// unedited) and the subset: no editing at run time, `--bam` refused. The
/// converter's own `source_type` key and its refusal are stated as this
/// layout's, not Ensembl VEP's.
#[test]
fn refseq_cache_fields_and_bam_notice_subset_is_documented() {
    assert_page_states(
        "cache-conversion.md",
        &[
            "The converter writes `source_type` into `info.json` from the species directory's \
             suffix: `refseq`, `merged` or `ensembl`, a key of its own with no counterpart in \
             Ensembl VEP's `info.txt`.",
            "A run whose `--refseq` or `--merged` disagrees with the declared set is refused \
             before any input is read, where Ensembl VEP filters whatever models the cache \
             holds by the flags alone; a cache converted without the key is checked against \
             nothing.",
            "A merged cache's transcripts carry `source_cache`, `Ensembl` or `RefSeq`, exactly \
             where the Storable has it: on the region files both caches held. A region file one \
             cache alone held carries no key, and its transcripts print no `SOURCE` in either \
             engine.",
            "A run on one prints `INFO: BAM-edited cache detected, enabling --use_transcript_ref; \
             use --use_given_ref to override this`, reads each row's reference allele from the \
             transcript and adds the `GIVEN_REF`, `USED_REF` and `BAM_EDIT` columns, as Ensembl \
             VEP does. `--use_given_ref` keeps the input's allele and drops `GIVEN_REF` and \
             `USED_REF`; `BAM_EDIT` stays on an edited transcript's rows in the default and \
             JSON formats, with no header line, as Ensembl VEP writes it.",
            "Its exon coordinate mapper's pairs fold the insertions and deletions in, read from \
             the Storable mapper object rather than rebuilt from the exons.",
            "`BAM_EDIT` prints `OK` on its rows and its reference allele is cut from the edited \
             cDNA, coding or not.",
            "A model the correction failed on, status `failed`, keeps its own attributes, which \
             can include an `_rna_edit`, loses only the alignment's edits and has no cDNA cached: \
             `BAM_EDIT` prints `FAILED` and its rows read the genome in vep-rs.",
            "No editing happens at run time: `--bam`, which edits the models against an \
             alignment file as the cache is read, is refused.",
        ],
    );
    assert_page_states(
        "output-formats.md",
        &[
            "A BAM-edited cache, a converted RefSeq or merged cache whose metadata names the \
             alignment file its models were corrected against, switches the flag on by itself \
             with an `INFO:` line, as in Ensembl VEP, and adds `BAM_EDIT`",
            "An edited transcript's allele is cut from its edited cDNA, coding or not. `--bam`, \
             which edits the models at run time, is refused.",
        ],
    );
}

/// `attributes-rseq-refseq-compare` and `refseq-match-bam-edit-source-columns`:
/// the `rseq_*` attributes print as `REFSEQ_MATCH`, the `enst_refseq_compare`
/// attribute is carried and printed by neither engine (`grep {$_->code =~
/// /^rseq/}`, OutputFactory.pm:1547), `SOURCE` and `REFSEQ_OFFSET` follow their
/// rules, and an edited model's HGVSc is numbered against the edited model,
/// with Ensembl VEP's unedited numbering stated as the intended divergence.
#[test]
fn refseq_match_source_and_offset_columns_subset_is_documented() {
    assert_page_states(
        "output/refseq-and-merged.md",
        &[
            "`--refseq` adds `REFSEQ_MATCH` and `REFSEQ_OFFSET` after `RefSeq`; `--merged` adds \
             `REFSEQ_MATCH`, `SOURCE` and `REFSEQ_OFFSET`.",
            "The cache also carries the `enst_refseq_compare` attribute of a model compared with \
             an Ensembl transcript, which neither engine prints.",
            "`SOURCE` names the cache a merged cache's transcript came from, `Ensembl` or `RefSeq`; \
             a single-source cache stamps no transcript, so the column stays empty there.",
            "It prints only when it is not 0, so a model whose edits cancel out, an allele upstream \
             of every edit and an intronic allele, which has no cDNA start, carry none.",
            "The HGVSc of a model whose alignment inserted or deleted bases is numbered against \
             the edited model, the transcript the identifier names. An exonic position takes its \
             edited coordinate and an intronic distance counts from the edited exon end.",
            "Its cache dumper strips the transcript's edit-aware mapper, so it walks the exons in \
             their unedited numbering and corrects an exonic position on an `NM_` or `XM_` model \
             by the offset above. An intronic allele, which has no offset, is anchored on the \
             unedited exon end, and an Ensembl model with such an edit is never corrected.",
            "The run then writes `WARNING: Transcript-assembly mismatch in <Uploaded_variation>` \
             to standard error and the warnings file, once per record, for a row the output \
             prints. A marked row that `--pick` or another record-level flag drops draws \
             nothing, as do `--most_severe` and `--summary`, which print no transcript row, and \
             a run under `--use_given_ref`.",
        ],
    );
}

/// `refseq` and `refseq-merged-cache-selection`: the pages state what a
/// `--refseq` run annotates (the RefSeq accessions, `AnnotationType/
/// Transcript.pm:236-248`, every model under `--all_refseq`, `:237`), that
/// under `--json_cache` the flags name the set of the cache given
/// (`CacheDir.pm:517`) and the subset: the caches are Ensembl's, converted, so
/// no RefSeq model is built natively, and `--bam`, Ensembl VEP's run-time
/// editing, is refused.
#[test]
fn refseq_exclusions_are_documented() {
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &[
            "vep-rs builds no RefSeq model natively and produces no merged cache.",
            "`--bam`, which corrects RefSeq transcript models against an alignment file as a \
             cache, GFF, GTF or database source is read, is refused with an error naming the \
             route (exit 1)",
        ],
    );
    assert_page_states(
        "cache-conversion.md",
        &[
            "A run whose `--refseq` or `--merged` disagrees with the declared set is refused \
             before any input is read",
            "Under `--json_cache` the flags name the set of the cache given, so `--refseq` on \
             a converted Ensembl cache without the `source_type` key keeps only the \
             transcripts that pass its test of a RefSeq id, as Ensembl VEP does.",
            "Of a cache of `ENST` ids those are the transcripts on `MT`, where the test passes \
             an id made of capitals and digits, so a record on any other chromosome gets an \
             intergenic row.",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &[
            "The converter reads the set from the species directory's suffix, so it is given \
             the directory inside that layout",
            "`--refseq` annotates the RefSeq accessions alone, `--merged` every Ensembl model \
             plus the RefSeq accessions among the RefSeq-sourced ones, and `--all_refseq` \
             every model",
            "vep-rs builds no RefSeq model natively and produces no merged cache.",
        ],
    );
    assert_page_states(
        "getting-started/coming-from-ensembl-vep.md",
        &[
            "`--refseq` and `--merged` select the `_refseq` or `_merged` directory and its \
             transcript set, under `--json_cache` too; `--all_refseq` keeps every model.",
        ],
    );
    assert_flag_row_states(
        "cli/filtering-picking.md",
        "--refseq",
        &[
            "of its models only the RefSeq accessions are kept",
            "unless `--all_refseq` keeps every model",
        ],
    );
    assert_flag_row_states(
        "cli/filtering-picking.md",
        "--all_refseq",
        &[
            "Keep every model of a RefSeq or merged cache",
            "changes nothing on an Ensembl cache",
        ],
    );
}

/// `merged`: the pages state what a `--merged` run annotates (every Ensembl
/// model, and of the RefSeq-sourced ones the accessions alone,
/// `AnnotationType/Transcript.pm:240`), that the `source_cache` stamp is read
/// where the Storable has it, and the subset: vep-rs produces no merged cache,
/// the specification of a producer `MergeVEP.pm` is.
#[test]
fn merged_exclusions_are_documented() {
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &["produces no merged cache"],
    );
    assert_page_states(
        "cache-conversion.md",
        &[
            "A merged cache's transcripts carry `source_cache`, `Ensembl` or `RefSeq`, exactly \
             where the Storable has it: on the region files both caches held.",
        ],
    );
    assert_flag_row_states(
        "cli/filtering-picking.md",
        "--merged",
        &[
            "every Ensembl model is kept, and of the RefSeq-sourced models (source\\_cache \
             RefSeq) only the RefSeq accessions",
        ],
    );
    assert_flag_row_states(
        "cli/filtering-picking.md",
        "--transcript_filter",
        &["`_source_cache` (`Ensembl` or `RefSeq` on a merged cache's transcripts, undefined elsewhere"],
    );
}

/// `filter_vep_consequence_field_pattern`: `vep filter` takes a VCF line's
/// consequence entries from every INFO item whose key is exactly the field's
/// name and whose value holds an entry, where filter_vep:261 finds them with
/// one unanchored `m/($vcf_info_field)\=(.+?)(\;|$|\s)/g` over the line, so an
/// empty entry before another item swallows that item and a key ending in the
/// name is read as the field, and :296 rewrites the first match under
/// `--only_matched`.
#[test]
fn filter_vep_consequence_field_pattern_is_documented() {
    assert_page_states(
        "divergences/filter-vep-consequence-field-pattern.md",
        &[
            "takes the entries from every INFO item whose key is exactly the field's name and \
             whose value holds an entry",
            "under `--only_matched` rewrites each such item to its own passing entries",
        ],
    );
}

/// `filter_vep_short_entry_reads_info_key`: `vep filter` reads a consequence
/// field the entry lacks as undefined, where filter_vep:389 starts each
/// entry's hash as a copy of the line's hash with its INFO keys, :395 writes
/// the fields the entry has and :413-416 sets a header field undefined only
/// `unless exists`, so a short entry reads the INFO key of a field's name
/// (`short_entry_field_is_undefined_not_the_info_key` in `filter_tool.rs`).
#[test]
fn filter_vep_short_entry_reads_info_key_is_documented() {
    assert_page_states(
        "divergences/filter-vep-short-entry-info-key.md",
        &[
            "sets a header field undefined only when the hash lacks the key",
            "`vep filter` reads a consequence field the entry lacks as undefined, whatever the \
             line's INFO column holds",
            "a name in the `Format` list is answered from the entry alone",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "The five tool entries count the tool runs of their tool",
            "so an entry shorter than the `Format` list reads the INFO key of the same name for a \
             field it lacks, the input's `AF` for the consequence field `AF`; `vep filter` reads a \
             field the entry lacks as undefined",
        ],
    );
    assert_page_states(
        "divergences/transcript-filter-fields.md",
        &["`vep filter` reads such a field as undefined"],
    );
}

/// `variant_recoder`: the `vep recode` page and the companion-tools section
/// state the subset the cache alone can give (`VariantRecoder.pm:96-109` forces
/// `database`, `merged` and `lrg`; `:425-437` reads the synonyms from the
/// database): no LRG notation (`--lrg` refused), one variant per protein
/// notation (`--ambiguous_hgvs` refused), the synonyms from the cache's
/// variation data, no ClinGen allele identifier input, `--refseq` refused as
/// the recoder's forced `merged` (`:98`) refuses it (`Config.pm:621`), the
/// merged set read without `--merged` wherever the cache named holds it, so
/// RefSeq notations only with a merged cache, and a GA4GH VRS allele only for
/// a chromosome whose synonyms name an `NC_` accession, the list empty for any
/// other.
#[test]
fn recode_exclusions_are_documented() {
    assert_page_states(
        "cli/recode.md",
        &[
            "The synonyms and the co-located variants are the cache's variation data, where \
             `variant_recoder` reads the Ensembl database.",
            "A RefSeq notation (`NM_`, `NP_`, `XM_`, `XP_`) appears only with a merged cache.",
            "An LRG notation never appears, as `--lrg` is refused, and a protein notation resolves \
             to one variant, as `--ambiguous_hgvs` is refused.",
            "A ClinGen allele identifier (`CA...`) is not an input.",
            "`ga4gh_vrs` holds a GA4GH VRS Allele object only for a chromosome whose synonyms name \
             an `NC_` accession, and is an empty list for any other.",
            "`--mane_select` is for human alone.",
            "`variant_recoder` reads the merged transcript set on every run, and `vep recode` \
             reads it wherever the cache named holds it, with `--merged` or without.",
            "Such a cache is a `--json_cache` directory whose `info.json` declares the merged set, \
             or a JSON cache in the `_merged` species directory under the cache root.",
            "Beside the Ensembl set's cache the run would read, it is the cache of the same \
             `<cache_version>_<assembly>` directory name, so the run keeps its cache version and \
             assembly; under a root with no directory of the Ensembl set, it is the one `--merged` \
             resolves.",
            "An Ensembl VEP Storable cache is not one, and a cache of the Ensembl set alone is read \
             as it is.",
            "Where the directory of the Ensembl set is there and holds no cache the run can read, \
             the run is refused as it is without a merged cache, and `--merged` reads the merged \
             one.",
            "`--refseq` is refused with `Can't use --refseq and --merged together`, as \
             `variant_recoder` refuses it, and with `--gencode_basic` the error names that flag \
             instead.",
            "A RefSeq cache, which is read under `--refseq` alone, is not an input; a merged cache \
             carries the RefSeq models beside the Ensembl ones.",
            "An option is read from the four sources `vep` reads its own from, the first that names \
             it winning: the command line, the `--config` file, `<dir>/vep.ini`, where `<dir>` is \
             the file's `dir` or else `~/.vep`, then the environment variable `VEP_<OPTION>`.",
            "The options `variant_recoder` settles before it reads the sources are set on the \
             command line alone, as there.",
            "`fields`, `buffer_size`, `failed`, `merged`, `lrg` and `ambiguous_hgvs` keep the values \
             it presets.",
            "`vcf_string`, `var_synonyms`, `mane_select` and `ga4gh_vrs` add their output fields \
             before the sources are read, so a source's value adds none.",
            "`pretty` and `vcf_string` are names its reader does not know.",
            "Every other option of this page but `--config` may be set from a source.",
            "A `vep` option this page lacks is dropped from every source, where `variant_recoder` \
             reads it and applies it ([the class](../divergences/recode-config-file.md)).",
        ],
    );
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &[
            "`variant_recoder` is `vep recode` ([its page](../cli/recode.md)): the same options \
             over the cache alone.",
            "`--lrg` and `--ambiguous_hgvs` are refused with the database connection flags.",
            "A RefSeq notation appears only with a merged cache, a ClinGen allele identifier is \
             not an input, and `ga4gh_vrs` is written only for a chromosome whose synonyms name \
             an `NC_` accession.",
        ],
    );
    assert_flag_row_states(
        "cli/recode.md",
        "--lrg",
        &["Refused: LRG mapping needs the Ensembl core database"],
    );
    assert_flag_row_states(
        "cli/recode.md",
        "--ambiguous_hgvs",
        &["Refused: resolving a protein notation to every matching variant needs the Ensembl database"],
    );
    assert_flag_row_states(
        "cli/recode.md",
        "--refseq",
        &["Refused: the recoder reads the merged transcript set on every run"],
    );
    assert_flag_row_states(
        "cli/recode.md",
        "--merged",
        &["a cache that holds the set is read without the flag"],
    );
}

/// `custom`: a `--custom` source is read against a cache, where Config.pm:993
/// lets the source alone satisfy the run's source check; a remote file and a
/// GFF or GTF transcript model are refused (`AnnotationSourceAdaptor.pm:378`
/// reads a `no_remote` option no flag sets; `GFF.pm:154-156` reads the model in
/// place of the cache).
#[test]
fn custom_exclusions_are_documented() {
    assert_page_states(
        "custom-annotation.md",
        &[
            "A `--custom` run reads its cache as every run does; Ensembl VEP accepts a `--custom` \
             source as a run's only annotation source",
            "Ensembl VEP reads it over the network; vep-rs reads local files only, so copy the \
             file to a local path",
            "which Ensembl VEP reads in place of its cache",
        ],
    );
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &[
            "A remote file (`http://`, `ftp://` or any URL) is refused with an error naming a \
             local copy, where Ensembl VEP reads it over the network",
            "A GFF or GTF transcript model is refused with the two-step cache recipe, where \
             Ensembl VEP reads it in place of its cache",
            "A `--custom` run still reads a cache, where Ensembl VEP accepts the source as the \
             run's only annotation source",
        ],
    );
}

/// `custom-output-fields`: `--custom` adds the `SOURCE` column (Constants.pm:101),
/// whose value is the transcript's `_source_cache` (OutputFactory.pm:1583), set
/// by a GFF or GTF model read in place of the cache and by a merged cache; a
/// plain cache's transcripts leave it empty, and a cache built from a GTF or GFF
/// carries no source name and leaves it empty too. A source named `SOURCE`
/// lists its names there (`:1301`), in the one column where Ensembl VEP's tab
/// and VCF writers append a second (Tab.pm:192-196, VCF.pm:445-449).
#[test]
fn custom_output_source_column_is_documented() {
    assert_page_states(
        "custom-annotation.md",
        &[
            "`--custom` also adds the `SOURCE` column, which names the cache a merged cache's \
             transcript came from, `Ensembl` or `RefSeq`, and is empty on a transcript of a \
             single-source cache, a cache built from a GTF or GFF included",
            "A source named `SOURCE` lists its record names in that column, under the one key; \
             Ensembl VEP's tab and VCF writers list the column twice, vep-rs once \
             ([the class](divergences/duplicate-output-column.md)).",
            "The VCF writes each key as an `##INFO` line after the plugin lines",
        ],
    );
}

/// `duplicate_output_column_collapsed`, `recode_config_file_vep_option`: a plugin
/// key or a `--custom` short name equal to a flag column prints one tab column
/// and one `CSQ` field (`fields.rs`, `fields_of_groups` and `set_plugin_fields`),
/// where Tab.pm:195 and VCF.pm:449 append the key a second time after the
/// uniquified flag columns (OutputFactory.pm:2530-2534); `vep recode` drops a
/// configuration line naming an option its own command line lacks (`argv.rs`,
/// `with_config_file`), where Config.pm:1100-1101 admits every `vep` option and
/// applies one the recoder never declares.
#[test]
fn duplicate_column_and_recode_config_file_classes_are_documented() {
    assert_page_states(
        "divergences/duplicate-output-column.md",
        &[
            "Ensembl VEP's tab and VCF writers append every plugin key and every `--custom` short \
             name to the column list after the flag columns, without testing whether the list \
             already carries the name.",
            "A plugin key or a custom short name that equals a flag column takes the flag's place \
             under its flag, and a plugin key handed twice is one column: one tab column, one \
             description line and one `CSQ` field, carrying the value the row holds.",
            "`Bio/EnsEMBL/VEP/OutputFactory/Tab.pm`, lines 187-202, `fields`: the 13 default columns \
             (line 193), the flag columns (194), then every plugin header key and every custom \
             header key (195), appended with no test against the names already listed.",
            "`Bio/EnsEMBL/VEP/OutputFactory/VCF.pm`, lines 434-455, `fields`: the flag columns are \
             filtered against the VCF's own fixed columns alone (445-447), and the plugin and custom \
             keys are appended untested (449).",
            "With REVEL named twice, Ensembl VEP 116.2's tab output has two `REVEL` columns",
            "Ensembl VEP's tab header carries `## SOURCE : Source of transcript` twice, vep-rs's once.",
            "Under `--tsl` with `TSL` asked of the plugin, Ensembl VEP 116.2's tab output has two `TSL` \
             columns and its `CSQ` list two `TSL` fields",
        ],
    );
    assert_page_states(
        "divergences/recode-config-file.md",
        &[
            "`variant_recoder` reads a configuration source through `vep`'s option list, so a line \
             naming a `vep` option its own command line never declares, `symbol 1` among them, is \
             read and applied: the recoder computes what the option asks for, and `--fields` can \
             print it.",
            "`vep recode` reads every configuration source as its own command line: an option of \
             [its page](../cli/recode.md) may be set there, the first source naming it wins, and a \
             name the page lacks is dropped from every source.",
            "A file naming the recoder's own options is read alike by both tools.",
            "`variant_recoder` prints a `gene_symbol` list beside `hgvsg`: the file switched on \
             `vep`'s `--symbol`, which its command line never offers.",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "so a source named `SOURCE`, a dbNSFP column named `TSL` under `--tsl`, or a plugin named twice heads \
             two tab columns, two description lines and two `CSQ` fields, each holding the row's one \
             value; vep-rs lists the column once",
            "`vep recode` reads a source as its own command line and drops the line",
        ],
    );
}

/// `plugin_cols_names_matched_exactly`: a `cols` name is matched to the
/// file's column names exactly (`columns_to_write` in `open_targets.rs`,
/// `select_columns` in `alpha_missense.rs`, the `cols` partition in
/// `transcript_annotator.rs`) and a name the header lacks is refused with
/// the valid names listed, where OpenTargets.pm:224, AlphaMissense.pm:158 and
/// TranscriptAnnotator.pm:136 test the name with `grep(/^$col$/, ...)` over the column names
/// and keep it as typed, so a name holding a regex metacharacter passes the
/// test and heads a column whose every value is undefined.
#[test]
fn plugin_cols_names_matched_exactly_class_is_documented() {
    assert_page_states(
        "divergences/plugin-cols-names-matched-exactly.md",
        &[
            "Ensembl VEP's AlphaMissense, OpenTargets and TranscriptAnnotator plugins test a \
             `cols` name against the file's column names as an anchored regular expression, so a \
             name holding a metacharacter heads a column of undefined values; vep-rs matches names \
             exactly and refuses one the header lacks.",
            "`OpenTargets.pm`, lines 219-224: `_parse_colnames` splits `cols` on `:`, stores the \
             names as typed and pushes a name to `@invalid_cols` unless `grep(/^$col$/, @{ \
             $self->{colnames} })` finds a column.",
            "`AlphaMissense.pm`, lines 149-163: the same test on the `am_`-prefixed names, the \
             refusal listing the valid columns.",
            "`TranscriptAnnotator.pm`, lines 134-150: `grep /^$col$/` over the column names keeps a \
             name as typed and warns about the others; the plugin dies only when no name passes.",
            "`MaveDB.pm`, line 123, and `PromoterAI.pm`, line 117, validate `cols` with the same \
             test; `EVE.pm`, line 136, tests `class_number` and `GO.pm`, line 131, tests `match` \
             the same way. None of the four is ported.",
            "OpenTargets_gwasGeneI.=NA,NA,NA,NA",
            "WARNING: Failed to instantiate plugin OpenTargets: OpenTargets: columns not found in \
             the header of <file>: gwasGeneI.",
            "Under `cols=gwas.*`, a pattern matching every GWAS column, Ensembl VEP writes one \
             column, `OpenTargets_gwas.*`, with every value `NA`",
            "Unmatched ( in regex",
            "WARNING: Failed to instantiate plugin AlphaMissense: AlphaMissense: the following \
             columns were not found in file header: am_pathogenicit.; valid columns are: am_genome, \
             am_uniprot_id, am_transcript_id, am_protein_variant, am_pathogenicity, am_class",
            "WARNING: Failed to instantiate plugin TranscriptAnnotator: TranscriptAnnotator: no \
             valid columns selected. Valid column names: SIFT_score, SIFT_pred, Comment",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Two hundred and forty-eight differences lie outside every class",
            "Each of the two hundred and forty-eight [differences outside the classes]",
            "so a name holding a regular-expression metacharacter (`gwasGeneI.`) passes the test \
             and heads a column, named by the pattern's text, whose every value is undefined; \
             vep-rs matches the name to the column names exactly and refuses one the header lacks, \
             listing the valid names | 3 · 0 of 392 flag sets · 15 flag sets, 0 corpora | 116.2 |",
        ],
    );
    assert_page_states(
        "plugins/annotation-store.md",
        &[
            "A name the header lacks is refused with the valid names listed; Ensembl VEP 116.2 \
             tests each name as a pattern ([a documented \
             class](../divergences/plugin-cols-names-matched-exactly.md)).",
            "A name the header lacks is refused with the valid names listed; Ensembl VEP 116.2 \
             tests each name as a pattern, so `cols=gwasGeneI.` loads there with every value \
             `NA` ([a documented class](../divergences/plugin-cols-names-matched-exactly.md)).",
            "A name is matched to the header names exactly; Ensembl VEP 116.2 tests each name as \
             a pattern and keeps it as typed, so `SIFT_.core` heads an empty column there \
             ([a documented class](../divergences/plugin-cols-names-matched-exactly.md)).",
        ],
    );
}

/// `tabix-and-bigwig-readers`: the tabix formats are read from local files
/// through their `.tbi` or `.csi` index, no BAM is read (`AnnotationType/
/// Transcript.pm:378-402` reads one for a GFF or GTF model) and no file is
/// fetched over `http` or `ftp` (`File.pm:73-94` opens either).
#[test]
fn custom_readers_exclusions_are_documented() {
    assert_page_states(
        "custom-annotation.md",
        &[
            "Files are read locally through their tabix index (`.tbi` or `.csi`); a BAM is not \
             read, and a remote file is refused as above",
        ],
    );
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &["A BAM is not read."],
    );
}

/// `stats-charts-distribution-per-chromosome`: the per-chromosome distribution
/// charts read the chromosome lengths from the FASTA's index alone, where
/// Ensembl VEP reads them from its database when one is connected
/// (`BaseVEP.pm:550-580`); offline without a FASTA neither engine writes one.
#[test]
fn stats_distribution_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "The `Distribution of variants on chromosome` charts need the chromosome lengths, \
             which vep-rs reads from the index of the `--fasta` file alone, one chart per \
             chromosome spelled as the FASTA names it. Ensembl VEP reads them from its database \
             when one is connected and from the FASTA otherwise (`BaseVEP.pm:550-580`). Without \
             a FASTA neither engine writes a distribution chart.",
        ],
    );
    assert_flag_row_states(
        "cli/statistics-verbosity-misc.md",
        "--stats_text",
        &["their distribution along each chromosome when a FASTA gives the lengths"],
    );
}

/// `stats-html-chart-technology`: the HTML report keeps `dump_html`'s sections,
/// ids, titles, tables and sort orders (`Stats.pm:846-915`) and draws its charts
/// as inline SVG, where Ensembl VEP's page loads a remote chart library and two
/// images (`Stats.pm:1001-1003, 1138, 1142`).
#[test]
fn stats_html_technology_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "The HTML report has the same sections, ids, titles, tables and sort orders as Ensembl \
             VEP's (`Stats.pm` `dump_html`), rendered as inline SVG without Ensembl VEP's remote \
             chart loader and images. The page is one file that loads no script, no chart library \
             and no image, so it opens without a network connection.",
        ],
    );
    assert_flag_row_states(
        "cli/statistics-verbosity-misc.md",
        "--stats_html",
        &["the default when neither `--stats_text` nor `--stats_html` is given"],
    );
}

/// `quiet-and-q`: an output named `stdout` does not force `--quiet`
/// (`Config.pm:933-937`), and `--quiet` leaves the statistics report complete,
/// as Ensembl VEP keys its counters on `--no_stats` alone (`Stats.pm:230, 277`).
#[test]
fn quiet_exclusion_is_documented() {
    assert_page_states(
        "cli/compatibility.md",
        &[
            "An output file matching `stdout` in any case does not force `--quiet` or drop \
             `--verbose` (`Config.pm:933-937`). The status lines go to standard error, so the \
             annotation on standard output is clean either way, and `--quiet --verbose` stays \
             refused with `-o STDOUT` as without. `--quiet` silences those lines and nothing \
             else: the statistics report is written in full under it, as Ensembl VEP keys its \
             counters on `--no_stats` alone (`Stats.pm:230, 277`).",
        ],
    );
}

/// `force-overwrite-scope`: the three guarded files are checked report first,
/// then output, then skipped-variants file, and created only once every check
/// passes (`Runner.pm:144-148` opens the report before `:202` opens the output;
/// `BaseRunner.pm:337-340` checks the skipped file as the run ends).
#[test]
fn force_overwrite_scope_exclusion_is_documented() {
    assert_flag_row_states(
        "cli/input-output.md",
        "--skipped_variants_file",
        &[
            "The statistics report is checked first, then the output file, then this file, and \
             the three are created only once every check passes",
        ],
    );
    assert_flag_row_states(
        "cli/input-output.md",
        "--force_overwrite",
        &["Overwrite the output file, the statistics report and the skipped-variants file when they exist"],
    );
    assert_page_states(
        "cli/compatibility.md",
        &[
            "The report is checked for an existing file before the output and the \
             skipped-variants file are, as Ensembl VEP checks it (`Runner.pm:144-148`, before \
             `:202`). The three files are created only once every check passes, so a refused \
             output leaves no empty report behind, where Ensembl VEP leaves one.",
        ],
    );
}

/// `conservation-shortcuts`: the pages state the two exclusions of the bigWig
/// shortcuts. The remote root Ensembl VEP reads by default (`Config.pm:314`) is
/// refused with the local tree's layout and the `--custom` form of the same source,
/// and a variant whose span the span rule's per-record walk covers is annotated
/// once where Ensembl VEP 116.2 lists its records twice (`BigWig.pm:195-204`).
/// The `--phyloP` short name is spelt `phyloP<N>way`, not Ensembl VEP's
/// `phlyoP<N>way` (`Config.pm:555`).
#[test]
fn conservation_exclusions_are_documented() {
    assert_page_states(
        "output/custom-bigwig.md",
        &[
            "The default root is Ensembl VEP's remote `http://hgdownload.cse.ucsc.edu/goldenpath/`.",
            "exits 1 with `ERROR: [custom] --ucsc_data_root must name a local directory laid out as \
             <root>/<assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw; the same source is \
             --custom file=<local>.bw,format=bigwig,type=exact`.",
            "A local mirror of the UCSC tree, or the bigWig itself through `--custom`, is the route.",
            "Every variant is annotated once on that path; Ensembl VEP 116.2 lists each record of a \
             shorter variant twice there.",
            "vep-rs prints the one form on both paths.",
            "The `--phyloP <N>` source's short name, its output key, is `phyloP<N>way`.",
            "Ensembl VEP's option set spells it `phlyoP<N>way` (`Config.pm:555`), a misspelling \
             vep-rs does not reproduce",
        ],
    );
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &[
            "`--phyloP` and `--phastCons` read their bigWigs from a local mirror alone: the remote \
             root Ensembl VEP reads by default is refused",
        ],
    );
    assert_flag_row_states(
        "cli/plugins.md",
        "--ucsc_data_root",
        &["The default is Ensembl VEP's remote root, which is refused: a local mirror must be named"],
    );
    assert_flag_row_states(
        "cli/plugins.md",
        "--phyloP",
        &["so the column is phyloP\\<N\\>way"],
    );
}

/// The text of the repository README, whitespace collapsed.
fn readme() -> String {
    std::fs::read_to_string(format!("{}/../../README.md", env!("CARGO_MANIFEST_DIR")))
        .unwrap()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// `species_statement`: the README and the species page state the genomes the
/// concordance measurement covers (human on both assemblies and the four other
/// golden corpora), what any other species gets (a cache and the codon tables
/// it names, no measurement), and how the flags that read human-only data
/// behave on another cache.
#[test]
fn species_statement_is_documented() {
    let statements = [
        "vep-rs is verified by concordance against Ensembl VEP 116.2 on human GRCh37 and GRCh38, the \
         whole-genome datasets of the current release, and on mouse GRCm39, fruit fly BDGP6.54, \
         Arabidopsis thaliana TAIR10 and yeast R64-1-1.",
        "Each of the four has a golden corpus of its own",
        "Any other species in Ensembl release 116 or Ensembl Genomes release 63 with a GFF3 and a \
         peptide FASTA, or with an Ensembl VEP cache to convert, builds a cache",
        "annotates under the codon tables the cache names, without a concordance measurement.",
        "`--sift` and `--polyphen` refuse on a cache without their data, and `--af_gnomad` on a \
         cache that declares variation columns without the gnomAD exome column; every other flag \
         that reads data Ensembl ships for human alone runs there as [the flag table](",
        "cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone) states.",
    ];
    assert_states("README.md", &readme(), &statements);
    assert_page_states("cache-species.md", &statements);
}

/// `os_and_architecture_matrix`: the installation page tables every platform
/// vep-rs runs on with its archive or image and what it needs, and states the
/// platforms it is not built for and the absence of a Homebrew formula.
#[test]
fn os_and_architecture_table_is_documented() {
    assert_page_states(
        "getting-started/installation.md",
        &[
            "| Platform | Archive or image | Runs on |",
            "| Linux x86_64 | `x86_64-unknown-linux-gnu` (measured, `x86-64-v3`) and \
             `x86_64-unknown-linux-gnu-v2` (generic, `x86-64-v2`) | glibc 2.17 or later |",
            "| Linux aarch64 | `aarch64-unknown-linux-gnu` (measured, `neoverse-v1`) and \
             `aarch64-unknown-linux-gnu-generic` (generic, Armv8.0-A) | glibc 2.17 or later |",
            "| macOS on Apple silicon | `aarch64-apple-darwin` | macOS 11 or later |",
            "| macOS on Intel | `x86_64-apple-darwin` (`x86-64-v3`) | an Intel Mac with AVX2 (every \
             Intel Mac macOS 26 supports) |",
            "| Container image | `ghcr.io/natera-open-source/vep-rs` on linux/amd64 and linux/arm64 | \
             Docker, Podman, Singularity or Apptainer |",
            "| Windows | none | the container image under Docker, or the Linux x86_64 archive under \
             WSL2; neither is tested on Windows |",
            "| FreeBSD, Linux with musl | none | not built |",
            "There is no Homebrew formula: macOS users install from the release archive or build from \
             source.",
        ],
    );
}

/// `docker-image`: the installation page and the README state the image's
/// platforms and tags (the multi-platform tags over linux/amd64 and
/// linux/arm64 with the generic binaries, the single-platform tags with the
/// measured ones), its contents and its attestation.
#[test]
fn container_image_platforms_are_documented() {
    assert_page_states(
        "getting-started/installation.md",
        &[
            "The `X.Y.Z`, `X.Y` and `latest` tags carry the generic Linux binaries (`x86-64-v2` on \
             linux/amd64, generic aarch64 on linux/arm64)",
            "The `X.Y.Z-x86-64-v3` and `X.Y.Z-neoverse-v1` tags carry the measured binaries for hosts \
             that meet their CPU floors.",
            "Every image holds the release's own `vep`, `vep-cache-builder` and `vep-cache-converter`",
            "add `bgzip` and `tabix`",
            "Tags are `X.Y.Z`, `X.Y` (the newest patch of that line) and `latest`, each a manifest list \
             over linux/amd64 and linux/arm64, plus the single-platform `X.Y.Z-x86-64-v3` (linux/amd64) \
             and `X.Y.Z-neoverse-v1` (linux/arm64).",
            "On Apple silicon the multi-platform tags run natively as linux/arm64.",
            "Each image carries a build-provenance attestation:",
            "gh attestation verify oci://ghcr.io/natera-open-source/vep-rs:",
            "| Container image | `ghcr.io/natera-open-source/vep-rs` on linux/amd64 and linux/arm64 | \
             Docker, Podman, Singularity or Apptainer |",
        ],
    );
    assert_states(
        "README.md",
        &readme(),
        &[
            "the `X.Y.Z`, `X.Y` and `latest` tags carry the generic Linux binaries for linux/amd64 and \
             linux/arm64 and run on Apple silicon and Graviton2 without a platform flag; \
             `X.Y.Z-x86-64-v3` and `X.Y.Z-neoverse-v1` carry the measured binaries.",
        ],
    );
}

/// `singularity-apptainer`: the installation page states how Apptainer and
/// Singularity pull the image into a single file, which binary each tag gives,
/// that the container runs as the invoking user with the current directory,
/// home and `/tmp` mounted, and that any other directory is bound with
/// `--bind`.
#[test]
fn singularity_and_apptainer_pull_is_documented() {
    assert_page_states(
        "getting-started/installation.md",
        &[
            "### Singularity and Apptainer",
            "`apptainer pull docker://ghcr.io/natera-open-source/vep-rs:X.Y.Z` writes `vep-rs_X.Y.Z.sif` \
             from the image for the host's architecture: the generic binary under the multi-platform tags, \
             the measured one under `X.Y.Z-x86-64-v3` or `X.Y.Z-neoverse-v1` on a host that meets its CPU \
             floor.",
            "The container runs as the invoking user with the host's current directory, home directory and \
             `/tmp` mounted, so the image's own user and `/data` working directory do not apply.",
            "A cache or output directory anywhere else is mounted with `--bind`:",
            "apptainer exec --bind /shared/vep-cache:/cache:ro --bind \"$PWD/results:/out\"",
            "`singularity pull` and `singularity exec` take the same arguments.",
        ],
    );
}

/// `memory_footprint_statement`: the installation page states the loading model
/// (the transcript index per chromosome, resident for the run; the variation
/// files per chromosome; an identifier input's name list) and where the two
/// measured figures live, the README's per-chromosome bound and the release
/// page's single-file bound; the README carries the per-chromosome bound as its
/// Peak memory row and points at the page.
#[test]
fn memory_statement_and_loading_model_are_documented() {
    assert_page_states(
        "getting-started/installation.md",
        &[
            "The transcript index is loaded one chromosome at a time, and every chromosome a run \
             touches stays resident until the run ends.",
            "Variation files are read per chromosome as needed, and an identifier input holds every \
             name of the list while it reads the variation files once.",
            "The index is the part of the peak that grows with the chromosomes an input touches.",
            "The README's Peak memory row is the largest per-chromosome median of the release \
             measurement, one chromosome per input file at 16 threads.",
            "A whole genome in one input file holds every chromosome's index at once.",
            "The release page states the measured peak of gnomAD v4.1 GRCh38 as a single file per \
             architecture, beside the per-dataset peaks of its wall-time table",
        ],
    );
    let readme = readme();
    let row = readme
        .split_once("| **Peak memory** | ")
        .map(|(_, rest)| rest.split(" |").next().unwrap_or(""))
        .expect("the README's Why vep-rs table has a Peak memory row");
    assert!(
        row.starts_with("At most ") && row.ends_with(" GB RSS on any chromosome"),
        "the Peak memory row states the per-chromosome bound: {row}"
    );
    assert_states(
        "README.md",
        &readme,
        &["[Installation](docs/getting-started/installation.md#memory) states what a run holds in memory"],
    );
}

/// `nextflow-nf-core-drop-in`: the workflow-managers page names the command
/// each wrapper issues and the lines that change: the nf-core module's
/// container override, `ext.args` and version probes, Ensembl's pipeline's
/// `--config` file and `vep filter` for `filter_vep`, and the Snakemake
/// wrapper's standard-output, standard-error, `--safe` and `--stats_file`
/// flags.
#[test]
fn workflow_manager_mapping_is_documented() {
    assert_page_states(
        "getting-started/workflow-managers.md",
        &[
            "--cache --cache_version <version> --dir_cache <cache> --fork <cpus>",
            "Every flag of that line is honoured by vep-rs",
            "The cache the module stages must be a vep-rs cache",
            "a `--gtf` input in place of the cache is refused",
            "container = 'ghcr.io/natera-open-source/vep-rs:X.Y.Z'",
            "ext.args = '--vcf --everything'",
            "the vep-rs image answers the `tabix` line alone",
            "`vep --help` prints no `ensembl-vep` line, so the module's `ensemblvep` version reads \
             empty",
            "The image carries no `Math::CDF`, so the `perl -MMath::CDF` probe fails, and Nextflow \
             fails the task with `Unable to evaluate output` after the annotation has been written",
            "A copy of the module patched with `nf-core modules patch ensemblvep/vep` drops the \
             `versions_perlmathcdf` output",
            "vep-rs reads the same `--config` file, one flag per line as `<name> <value>` without the \
             dashes, and `vep filter` takes `filter_vep`'s options under the same spellings",
            "its `filter_vep` command becomes `vep filter`, since the image carries no executable of \
             that name",
            "--output_file STDOUT --warning_file STDERR --safe --stats_file <stats>",
            "`--output_file STDOUT` writes the annotation on standard output, `--warning_file STDERR` \
             keeps the warnings on standard error, `--safe` makes a plugin that cannot be loaded an \
             error, and `--stats_file` names the statistics report",
            "A cache built with `vep-cache-builder --auto` into the rule's cache directory is \
             therefore found as the wrapper expects",
            "The Perl modules the `bio/vep/plugins` rule downloads into `--dir_plugins` are not loaded",
            "vep-rs is the Bioconda package `vep-rs`",
        ],
    );
}

/// `plugin-install-and-data`: the data-files page states that no plugin data is
/// downloaded for the user by either engine, as INSTALL.pl installs a plugin's
/// module and prints that the plugin requires data (`INSTALL.pl:1925-1927`).
#[test]
fn plugin_data_download_exclusion_is_documented() {
    assert_page_states(
        "plugins/data-files.md",
        &[
            "It does not download from the upstream databases itself, and neither does Ensembl VEP's \
             installer: `INSTALL.pl` fetches a plugin's module and prints that the plugin requires \
             data, leaving the download to the user, so no plugin data is downloaded for the user by \
             either engine.",
        ],
    );
}

/// `bioconda-distribution`: the limitations page states that the Bioconda
/// package's version can trail the newest release, the forms every release is
/// published in, and what a pipeline that installs Ensembl VEP from Bioconda
/// changes.
#[test]
fn bioconda_status_is_documented() {
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &[
            "The Bioconda channel carries vep-rs as the package `vep-rs`, and the version it serves \
             can trail the newest release: the channel builds a release once its own copy of the \
             recipe names it ([Installation](installation.md#bioconda)).",
            "Every release is published as pre-built binaries on the \
             [Releases](https://github.com/natera-open-source/vep-rs/releases) page, as source and \
             as the container image ([Installation](installation.md)).",
            "A pipeline that installs Ensembl VEP from Bioconda does not get vep-rs by that package: \
             its environment names `vep-rs` in place of `ensembl-vep`, or it runs vep-rs from the \
             container image or the release archive ([Workflow managers](workflow-managers.md)).",
        ],
    );
}

/// `bioconda-distribution`: the installation page states the channel's package,
/// its install line and why its version can trail a release; it and the
/// limitations page state what the repository's Bioconda recipe builds, and the
/// home page names the recipe among what the installation page covers.
#[test]
fn bioconda_recipe_is_documented() {
    assert_page_states(
        "getting-started/installation.md",
        &[
            "The repository holds a Bioconda recipe, \
             [`bioconda/vep-rs`](../../bioconda/vep-rs/meta.yaml).",
            "It builds `vep`, `vep-cache-builder` and `vep-cache-converter` from the source \
             tarball of the published release its `version` names",
            "The Bioconda channel carries vep-rs as the package `vep-rs`, built for Linux and macOS \
             on x86_64 and arm64.",
            "conda install -c conda-forge -c bioconda vep-rs",
            "The channel builds a release from its own copy of the recipe, `recipes/vep-rs` in \
             [bioconda-recipes](https://github.com/bioconda/bioconda-recipes), once that copy names \
             the release, so the version it serves can trail the newest release.",
        ],
    );
    assert_page_states(
        "getting-started/what-vep-rs-does-not-do.md",
        &[
            "The repository holds a Bioconda recipe, which builds the binaries from a published \
             release's source tarball ([Installation](installation.md#bioconda)).",
        ],
    );
    assert_page_states(
        "index.md",
        &[
            "covers the CPU floor, the container tags, the Bioconda recipe and building for \
             your own machine",
        ],
    );
}

/// `gff-gtf-transcript-source`: `filter_transcript` (AnnotationType/Transcript.pm:215-251)
/// applies `--gencode_basic` and `--gencode_primary`, and only Cache/Transcript.pm:241
/// and Database/Transcript.pm:251 call it, so Ensembl VEP's `--gff` and `--gtf`
/// sources annotate every transcript under either flag; a cache built from the
/// file honours them, a documented intended divergence.
#[test]
fn file_cache_gencode_flags_divergence_is_documented() {
    assert_page_states(
        "divergences/file-cache-honours-gencode-basic.md",
        &[
            "Ensembl VEP 116.2 applies `--gencode_basic` and `--gencode_primary` in one routine, \
             `filter_transcript`, that its cache and database sources call and its `--gff` and \
             `--gtf` sources never do, so under either file source both flags are accepted and \
             change nothing: every transcript of the file is annotated.",
            "vep-rs reads a GTF or GFF3 into its cache and annotates against that cache, so the two \
             flags apply as they do on every cache: a transcript without the attribute is left out \
             of the annotation.",
            "Under `--gencode_basic`, Ensembl VEP 116.2 reading the GTF through `--gtf` writes the \
             same two rows, and vep-rs annotating the cache built from that GTF writes the \
             `ENST00000419431` row alone.",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "`--gencode_basic` and `--gencode_primary` filter the built cache as they filter every \
             cache, where Ensembl VEP's `--gtf` and `--gff` sources accept both flags and annotate \
             every transcript, a documented intended divergence, \
             [`file_cache_honours_gencode_basic`]\
             (divergences/file-cache-honours-gencode-basic.md)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP applies `--gencode_basic` and `--gencode_primary` in the one routine its \
             cache and database sources call and its `--gff` and `--gtf` sources never do, so under \
             either file source both flags are accepted and every transcript of the file is \
             annotated; a cache vep-rs builds from the same file honours the flags as every cache \
             does, and the rows on the transcripts without the tag are not written | 2 · 876 of \
             77,675 cache transcripts · 2 flag sets, 2 corpora | 116.2 |",
        ],
    );
}

/// `gff-gtf-transcript-source`: `_record_get_parent_ids` (GFF.pm:242-256) keeps the
/// empty first field of `Parent=,gene:X`, `_get_parent_child_structure`
/// (BaseGXF.pm:385-401) files the record under the gene and the prune at :430-435
/// reads the first field alone, so `_create_transcripts` (:311-345) visits the
/// record twice and :336 leaves it the gene of the last visit, by hash order;
/// the builder makes it the named gene's transcript once, a documented intended
/// divergence.
#[test]
fn gff_parent_list_empty_field_divergence_is_documented() {
    assert_page_states(
        "divergences/gff-parent-list-with-an-empty-field.md",
        &[
            "Ensembl VEP 116.2's `--gff` source links a record to every set identifier in its \
             `Parent` list but keeps the record at the top level whenever the list's first field \
             is empty. A transcript line written `Parent=,gene:X` is therefore visited twice, as \
             the gene's child and as a record of its own, and its gene identifier, symbol and \
             gene version are those of whichever visit comes last in Perl's hash order.",
            "vep-rs reads a transcript line whose `Parent` list opens with an empty field as the \
             transcript of the gene line a later field names, once, and gives it that gene's \
             identifier, symbol and version.",
            "Ensembl VEP 116.2 reading the file through `--gff` writes the same row under one hash \
             seed and a row with `-` for the gene and no `SYMBOL` under another, with a Perl \
             warning about an undefined value.",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP's `--gff` source links a transcript line whose `Parent` list opens with an \
             empty field, `Parent=,gene:X`, to the gene it names and keeps it at the top level as \
             well, so the transcript is visited twice and its `Gene`, `SYMBOL` and gene version are \
             the gene's or empty by Perl hash order; a cache vep-rs builds from the file makes the \
             transcript the named gene's, once | 3 · 0 of 964 table rows · 14 flag sets, 0 corpora \
             | 116.2 |",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "A record whose list opens with an empty field is a record of its own, or the transcript \
             of a gene a later field names. Ensembl VEP's `--gff` source gives such a transcript \
             the gene's identifiers or none by hash order, a documented intended divergence, \
             [`gff_parent_list_with_an_empty_field`]\
             (divergences/gff-parent-list-with-an-empty-field.md).",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &[
            "| A record whose parent list opens with an empty field (`Parent=,g1`) | kept as a record \
             of its own, or as the transcript of a gene a later field names \
             ([`gff_parent_list_with_an_empty_field`]\
             (divergences/gff-parent-list-with-an-empty-field.md)) |",
        ],
    );
}

/// `gff-gtf-transcript-source`: `prefetch_transcript_ids` (Database/Transcript.pm:616-624)
/// is the only writer of a transcript's `_ccds`, which OutputFactory.pm:1534-1537
/// prints under `--ccds` and :740 ranks at the `ccds` step of `--pick`, and
/// `_add_identifiers` (BaseGXF.pm:632-703) reads no `ccdsid`, so a `--gff` run prints
/// no `CCDS` and ranks every transcript as without one; the builder reads the
/// attribute, a documented intended divergence.
#[test]
fn file_cache_ccdsid_divergence_is_documented() {
    assert_page_states(
        "divergences/file-cache-reads-ccdsid.md",
        &[
            "Ensembl VEP 116.2 fills a transcript's CCDS identifier from the database's xrefs \
             alone, when it reads the database or builds its cache, so its `--gff` and `--gtf` \
             sources never set one. A transcript built from an Ensembl GFF3 line that carries \
             `ccdsid=CCDS42939.1`, or from a GTF line that carries `ccds_id \"CCDS42939\"`, \
             therefore prints no `CCDS` under `--ccds` and ranks as a transcript without a CCDS \
             at the `ccds` step of `--pick`.",
            "vep-rs reads an Ensembl GFF3 line's `ccdsid` or a GTF line's `ccds_id` into the \
             cache as the transcript's CCDS identifier, the value as the line carries it. It \
             prints the identifier under `--ccds` and ranks the transcript at the `ccds` step of \
             `--pick` as Ensembl VEP ranks the same transcript from its cache.",
            "With the two transcripts `ENST00000896630`, `ccdsid=CCDS13686.1`, and \
             `ENST00000896631`, which has none, under `--ccds --pick`, vep-rs picks \
             `ENST00000896630` and Ensembl VEP picks `ENST00000896631`.",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "So is the CCDS identifier, an Ensembl GFF3's `ccdsid` or a GTF's `ccds_id` as the line \
             carries it, which Ensembl VEP's `--gff` and `--gtf` sources never print or rank under \
             `--pick`, a documented intended divergence, \
             [`file_cache_reads_ccdsid`](divergences/file-cache-reads-ccdsid.md)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP fills a transcript's CCDS identifier from the database's xrefs alone, so its \
             `--gff` and `--gtf` sources never set one and a transcript built from an Ensembl GFF3 \
             line carrying `ccdsid` or a GTF line carrying `ccds_id` prints no `CCDS` under `--ccds` \
             and ranks as a transcript without a CCDS at the `ccds` step of `--pick`; a cache vep-rs \
             builds from the file reads the attribute, prints it and ranks the transcript as Ensembl \
             VEP's cache does | 2 · 235 of 77,675 cache transcripts · 1 flag set, 2 corpora | 116.2 |",
        ],
    );
}

/// `gff-gtf-transcript-source`: `%INCLUDE_FEATURE_TYPES` (GFF.pm:71-116 at
/// release/116.2) lacks `ncRNA_gene`, the type an Ensembl GFF3 gives every
/// non-coding gene, so `_get_records_by_coords` (BaseGXF.pm:224) drops the gene
/// line and `_get_parent_child_structure` (BaseGXF.pm:430-435) prunes every
/// transcript whose `Parent` names it; release 116.3 adds the type, and the
/// builder reads the line as a gene, a documented intended divergence.
#[test]
fn file_cache_ncrna_gene_lines_divergence_is_documented() {
    assert_page_states(
        "divergences/file-cache-reads-ncrna-gene-lines.md",
        &[
            "Ensembl VEP 116.2's `--gff` source reads a line only when its type is on a fixed \
             list, and `ncRNA_gene`, the type an Ensembl GFF3 gives every non-coding gene, is \
             not on it. The gene line is skipped, every transcript whose `Parent` names it is \
             pruned, and a record inside the gene is written as `intergenic_variant`.",
            "vep-rs reads an `ncRNA_gene` line as a gene, as it reads `gene` and `pseudogene`, \
             and the cache it builds carries the gene's transcripts.",
            "Of the 646,577 transcript lines of the GRCh38 release-116 GFF3, 199,386 belong to \
             the 42,155 genes written as `ncRNA_gene`, and the file source annotates none of \
             them. Release 116.3 adds `ncRNA_gene` to the list.",
            "vep-rs annotating the cache built from the file writes two \
             `non_coding_transcript_exon_variant` rows; Ensembl VEP 116.2 reading the file \
             through `--gff` writes one `intergenic_variant` row.",
            "The committed caches hold 313 transcripts built from a file whose biotype an Ensembl \
             GFF3 writes under an `ncRNA_gene` line. 23 are in the cache built from a GTF, 266 in \
             the cache built from a GFF3 and 24 in a mitochondrial shard built from the GFF3 that \
             a third corpus splices into a converted cache.",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "Ensembl VEP 116.2's `--gff` source skips an `ncRNA_gene` line and drops every \
             transcript whose `Parent` names one, a documented intended divergence, \
             [`file_cache_reads_ncrna_gene_lines`](divergences/file-cache-reads-ncrna-gene-lines.md)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "so the gene line is skipped, every transcript whose `Parent` names it is pruned and a \
             record inside the gene is written as `intergenic_variant`; a cache vep-rs builds from \
             the file reads the line as a gene, as release 116.3 does, and annotates its \
             transcripts | 2 · 313 of 77,675 cache transcripts · 0 flag sets, 3 corpora | 116.2 |",
        ],
    );
}

/// `gff-gtf-transcript-source`: `--gtf` and `--gff` are `--custom` sources
/// (Config.pm:531-543) whose short name, the file's base name (File.pm:154-158),
/// `_add_identifiers` stamps on every transcript (BaseGXF.pm:635) and
/// OutputFactory.pm:1583 prints as `SOURCE`, with the source's own column and
/// header lines (Constants.pm:101, OutputFactory.pm:2435-2453); a cache built
/// from the file names no source, so `without_transcript_source` in
/// `common/mod.rs` sets the field aside, a documented intended divergence.
#[test]
fn file_cache_no_source_field_divergence_is_documented() {
    assert_page_states(
        "divergences/file-cache-writes-no-source-field.md",
        &[
            "Ensembl VEP's `--gtf` and `--gff` register the file as a custom annotation source \
             and stamp every transcript read from it with the file's base name. That name is \
             printed as a `SOURCE` field on every transcript row, with a `SOURCE` column, an \
             empty column named after the file and header lines describing both, none of which \
             a run from a cache writes.",
            "A cache vep-rs builds from a GTF or GFF3 records no source name, and a run on it \
             writes `SOURCE` only from a merged cache's `Ensembl` or `RefSeq` stamp, as Ensembl \
             VEP's own cache run does.",
            "Its tab output has 29 columns, `SOURCE` the 28th, filled on every transcript row, and \
             `chr21_genes.gtf.gz` the 29th, `-` on every row; vep-rs's has the other 27. Its VCF \
             `CSQ` entry has 31 fields, `SOURCE` the 30th and `chr21_genes.gtf.gz` the 31st, \
             empty,",
            "The class is coded: the golden harness sets the field, the column and the header \
             lines aside before comparing a run on a file-built cache with Ensembl VEP reading \
             the file.",
            "The committed caches hold 1,082 transcripts built from a file (of 77,675): 1,045 in \
             the two caches built whole from a GTF or a GFF3, and 37 in a mitochondrial shard \
             built from the GFF3 that a third corpus splices into a converted cache.",
            "A corpus with such a reference names the class in its manifest's `masks`, so the \
             comparison states what it sets aside.",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "Ensembl VEP's file source prints a `SOURCE` field naming the file on every \
             transcript row, which a run on the built cache does not, a documented intended \
             divergence, \
             [`file_cache_writes_no_source_field`](divergences/file-cache-writes-no-source-field.md)",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "a cache vep-rs builds from the file names no source, so a run on it writes none of \
             them | 2 · 1,082 of 77,675 cache transcripts · 2 flag sets, 3 corpora | 116.2 |",
        ],
    );
}

/// `gff-gtf-transcript-source`: `_add_translation` (BaseGXF.pm:719-734) reads
/// `protein_id` from the first ordered CDS record, undefined when absent (:723),
/// and builds the translation without an identifier; `hgvs_protein`
/// (TranscriptVariationAllele.pm:1748, :1759) names the reference from the
/// translation's `display_id`, the empty string (Translation.pm:907-910), plus
/// the version, so Ensembl VEP writes `.1:p.Pro6Leu`; the built cache writes no
/// HGVSp for a translation without an identifier, a documented intended
/// divergence.
#[test]
fn gff_cds_without_protein_id_hgvsp_divergence_is_documented() {
    assert_page_states(
        "divergences/gff-cds-without-protein-id-hgvsp.md",
        &[
            "Ensembl VEP 116.2's `--gff` source takes a translation's protein identifier from the \
             transcript's first CDS line in transcript order, and when that line carries no \
             `protein_id` it builds the translation with none. Its HGVSp writer then names the \
             reference sequence as the empty identifier followed by the version, so the protein \
             notation is written `.1:p.Pro6Leu`, a change on a sequence the notation does not \
             name.",
            "A cache vep-rs builds from the file and the genome writes no `HGVSp` and no `ENSP` for \
             a translation whose first CDS line names no protein, as both engines write none for a \
             transcript without a translation. A build of the cache's content that reads the core \
             table dumps takes their identifier and writes both.",
            "Ensembl VEP 116.2 reading the file through `--gff` writes `HGVSp=.1:p.Pro6Leu` and \
             no `ENSP`; vep-rs writes neither.",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "A line without `protein_id` gives the translation none, and its rows no HGVSp, unless a \
             build of the cache's content reads the core dumps, which name it; Ensembl VEP writes an \
             HGVSp on the empty identifier, a documented intended divergence, \
             [`gff_cds_without_protein_id_hgvsp`]\
             (divergences/gff-cds-without-protein-id-hgvsp.md).",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "and its HGVSp writer names the reference as the empty identifier plus the version, \
             `.1:p.Pro6Leu`; a cache vep-rs builds from the file and the genome writes no HGVSp \
             for such a translation, as both engines write no ENSP | 3 · 0 of 964 table rows · \
             0 flag sets, 0 corpora | 116.2 |",
        ],
    );
}

/// `gff-gtf-transcript-source`: `_create_transcript` stores each CDS child under
/// its start (BaseGXF.pm:511), so of two CDS lines of one transcript at one start
/// one replaces the other (the later line of the file when both follow their
/// transcript line; the earlier when two same-`ID` lines are read before it,
/// :412), and `_add_translation` takes
/// the identifiers from the first kept record (:723-734) and the translation's
/// end from the last (:749-756); the builder reads every line, the coding span
/// their union and the identifiers from the line that opens the coding sequence,
/// a documented intended divergence.
#[test]
fn gff_cds_lines_sharing_a_start_divergence_is_documented() {
    assert_page_states(
        "divergences/gff-cds-lines-sharing-a-start.md",
        &[
            "Ensembl VEP 116.2's `--gff` source stores a transcript's CDS lines in a hash keyed by \
             their start, so when two CDS lines of one transcript share a start one replaces the \
             other. The later line of the file is kept when the two follow their transcript \
             line, as a sorted file reads them. Everything the dropped line alone stated is lost \
             without a warning: its protein identifier and version when the pair opens the \
             coding sequence, and its end when the pair closes it.",
            "vep-rs reads every CDS line of the transcript: the coding span is the union of the \
             lines, so a base any line calls coding is coding. The translation's identifier and \
             version come from the line that opens the coding sequence, the lowest start on the \
             forward strand and the highest on the reverse. Where two lines share that start, \
             the forward strand takes the earlier line in the file and the reverse strand the \
             later,",
            "On `pwp2_codon186`, `21 44115741 G>A`, between the two ends, Ensembl VEP writes \
             `3_prime_UTR_variant` and `c.*3G>A` with no HGVSp, vep-rs `missense_variant`, \
             `c.557G>A` and `p.Arg186His`.",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "CDS lines sharing a start. Two CDS lines of one transcript at one start are both \
             read: every base either covers is coding, and the translation is named by the \
             earlier line of the file on the forward strand and by the later on the reverse. \
             Ensembl VEP keeps one of the two lines alone and loses what the other stated, a \
             documented intended divergence, \
             [`gff_cds_lines_sharing_a_start`](divergences/gff-cds-lines-sharing-a-start.md).",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "a cache vep-rs builds from the file reads every line, the coding span their union \
             and the identifiers from the line that opens the coding sequence | 3 · 0 of 964 \
             table rows · 14 flag sets, 0 corpora | 116.2 |",
        ],
    );
}

/// `gff-gtf-transcript-source`: `_create_transcript` takes a transcript line's
/// `version`, else `transcript_version` (BaseGXF.pm:481), `_add_identifiers`
/// the gene line's `version` (:639) and `_add_translation` the first CDS
/// record's `version`, `|| 1` (:727), each as text; the `Gene`, `Feature` and
/// `ENSP` writers (OutputFactory.pm:1469, :1465, :1560-1565) and the HGVS
/// reference names (TranscriptVariationAllele.pm:1499-1500, :1759-1760) append
/// it unless the identifier already ends in `.` and digits, so Ensembl VEP
/// prints `ENST00001113783.v2`; the builder holds every version as a number
/// and refuses to build a cache from text that is not one, a documented
/// intended divergence.
#[test]
fn gff_version_not_a_number_divergence_is_documented() {
    assert_page_states(
        "divergences/gff-version-not-a-number.md",
        &[
            "Ensembl VEP 116.2's `--gff` source takes the text of a `version` attribute as the \
             version of the gene, the transcript or the translation without reading it as a \
             number.",
            "vep-rs holds every version as a number. A `version` whose text is not a whole number \
             refuses the cache build with the error `<file> line <n>: version \"<text>\" is not \
             a whole number`, on a gene, transcript or CDS line alike, so no cache is written \
             from the file. A `transcript_version` beside such a `version` is read by neither \
             engine. A whole number above 4294967295, the largest the cache holds, is refused \
             the same way, naming the limit.",
            "A `transcript_version` of `0`, or an empty one, with no `version` beside it reaches \
             one field alone. Ensembl VEP writes `HGVSc=ENST00001113783.0:c.17C>T`, or \
             `ENST00001113783.:c.17C>T` for the empty value, and `ENST00001113783` in `Feature` \
             even under `--transcript_version`; vep-rs writes the bare identifier in both.",
            "Ensembl VEP 116.2 reading the file through `--gff` writes `ENSP=ENSP00000783588.v2` \
             and `HGVSp=ENSP00000783588.v2:p.Pro6Leu`; vep-rs refuses to build a cache from the \
             file at that CDS line.",
            "On `tmprss3_codon6`, Ensembl VEP 116.2 writes `ENST00001113783.v2` in `Feature` and \
             `HGVSc=ENST00001113783.v2:c.17C>T`; vep-rs refuses to build a cache from the file at \
             the transcript line.",
            "On `tmprss3_codon6`, Ensembl VEP 116.2 writes `ENSG00000160183.v2` in `Gene`; vep-rs \
             refuses to build a cache from the file at the gene line.",
            "Without `--transcript_version` Ensembl VEP writes `ENST00001113783` in `Feature` and \
             `.v2` in the HGVSc reference name still.",
            "The cache holds each version as a 32-bit number, so digits above 4294967295 on any \
             of the three lines end the build with an error naming the line and the limit, \
             where Ensembl VEP prints the digits",
            "vep-rs refuses to build the cache: the build ends at the first such line with an \
             error naming the file, the line and the limit, and writes no row.",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "The transcript version is `version`, else `transcript_version`, and the gene \
             version is the gene line's `version`; a GTF's `gene_version` is not read. A \
             `version` that is not a whole number, or is above 4294967295, the largest the cache \
             holds, refuses the build, naming the file and line, as does a protein FASTA \
             header's version above it; Ensembl VEP prints the text, \
             [`gff_version_not_a_number`](divergences/gff-version-not-a-number.md).",
            "else none; a `version` that is not a whole number refuses the build, where Ensembl \
             VEP prints its text, \
             [`gff_version_not_a_number`](divergences/gff-version-not-a-number.md).",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "vep-rs holds a version as a number and refuses to build a cache from a file whose \
             `version` is set to text that is no whole number, naming the line and the text \
             (`version \"v2\" is not a whole number`); `01` prints as `1`, and digits above \
             4294967295, the cache's 32-bit limit, are refused the same way, naming the limit, \
             where Ensembl VEP prints them | 3 · 0 of 964 table rows · 14 flag sets, 0 \
             corpora | 116.2 |",
        ],
    );
}

/// `gff-gtf-transcript-source`: the `--gff` source reads the ninth column
/// through `GFF3::get_attributes` (GFF3.pm:190-198), which splits the raw
/// column on `;` and `=` (:194) and decodes nothing, where the same parser's
/// `get_attribute_by_name` (:174-183) and `get_type` (:161-164) decode through
/// `decode_string`, `uri_unescape` (GXF.pm:303-307); `_add_identifiers`
/// (BaseGXF.pm:632-703) takes the symbol from that hash, so a gene line's
/// `Name=A%2CB` prints `A%2CB`, and the `gencode_promoter` overlay's
/// `_create_records` (GFF.pm:198-201) takes `ID` and `gene_id` from the same
/// hash (BaseGXF.pm:289); the builder and the overlay decode every `%XX`
/// escape, in either letter case, after the tests Ensembl VEP runs on the
/// written value, a documented intended divergence. A GTF value is decoded by
/// both engines (GTF.pm:117-135), every quote removed, and the builder reads
/// it that way.
#[test]
fn gff_attribute_escapes_read_raw_divergence_is_documented() {
    assert_page_states(
        "divergences/gff-attribute-escapes-read-raw.md",
        &[
            "Ensembl VEP 116.2's `--gff` and `--custom` GFF promoter sources read every attribute \
             value of the ninth column as written, so a value carrying a GFF3 escape prints the \
             escape in place of the character: `Name=A%2CB` prints the symbol `A%2CB` and \
             `ID=PROM%2c1%3E` the `feature_id` `PROM%2c1%3E`.",
            "The escaping is the URL convention: any byte may be written as `%XX`, the \
             hexadecimal digits in either case, and a reader decodes every such escape, `%3E` \
             for `>` as much as `%2c` for `,`.",
            "vep-rs decodes every `%XX` escape in every attribute value it reads from a GFF3, in \
             either letter case and in one pass, as the parser's own `decode_string` reads one. \
             So the symbol, identifiers, biotype and support level of a cache built from the \
             file, and the `feature_id` and `associated_gene` of a promoter window a `--custom` \
             GFF source overlays, are the text the file encodes. Ensembl VEP's tests on a value, \
             whether it is set and where its list splits on `,`, run on the written value here \
             too; the decoding follows them. The pair split is the parser's, so `Name=a=b` reads \
             `a` on both engines, and only the decoding differs.",
            "`GFF.pm`, lines 198-201, `_create_records`: a `gff_type=gencode_promoter` source's \
             `feature_id` is the hash's `ID` and its `associated_gene` the hash's `gene_id`, the \
             same hash of `BaseGXF.pm`, line 289, so a promoter window's two fields print as \
             written.",
            "The window's two fields print `PROM%2c1%3E` and `ENSG%3B1` on Ensembl VEP in every \
             format, the lower-case escape and the `%3E` as written. vep-rs prints `PROM,1>` in \
             the default, tab and JSON formats and `PROM&1>` in its `CSQ` entry, since a `>` is \
             re-spelt nowhere. It prints `ENSG;1` in the tab and JSON formats. In the default \
             format it prints `ENSG%3B1`, since the `Extra` column re-spells a `;` inside a value \
             that way on both engines, and in its `CSQ` entry the same, the spelling every `CSQ` \
             value gives a `;`.",
            "On the record `tmprss3_codon6`, `21 42395401 G>A`, a `missense_variant` on that \
             transcript, vep-rs writes `SYMBOL=A,B`. Ensembl VEP 116.2 reading the file through \
             `--gff` writes `SYMBOL=A%2CB` in every format",
            "The Saccharomyces cerevisiae R64-1-1 GFF3 names three genes with a comma \
             (`ADE5%2C7`, `ARG5%2C6`, `DUR1%2C2`) and the Arabidopsis thaliana TAIR10 GFF3 115 \
             genes with a semicolon (`CYCA2%3B3` among them).",
            "In the VCF output Ensembl VEP's `CSQ` entry carries `A%2CB` and vep-rs's `A&B`, the \
             spelling both engines give a `,` inside a `CSQ` value.",
            "so does the same comma in a GTF, `gene_name \"A%2CB\"`, which both engines read as \
             `A,B`.",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "Ensembl VEP's `--gff` source and its `--custom` GFF source of promoter windows read \
             every attribute value of the ninth column as the file writes it, so a gene or \
             transcript line's `Name`, `gene_name`, identifier, biotype or support level, or a \
             promoter window's `ID` or `gene_id`, holding a GFF3 escape prints the escape in \
             place of the character, `A%2CB` for the symbol a file spells `A%2CB`; a cache \
             vep-rs builds from the file, and its promoter overlay, decode every `%XX` escape \
             and print `A,B`, the name Ensembl's own cache carries | 3 · 0 of 1,022 table rows · \
             15 flag sets, 0 corpora | 116.2 |",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "A GTF value is read as `GTF.pm` reads it: the text from the key to the next `;`, \
             every `\"` removed wherever it stands, the `%XX` escapes decoded, so `\"A%2CB\"` is \
             `A,B` and `\"%30\"` the absent `0`. A line repeating a key read as one value is \
             refused, naming the key;",
            "A GFF3 value is the text between its pair's first and second `=`, as `GFF3.pm` cuts \
             it, its `%XX` escapes decoded; Ensembl VEP's `--gff` source prints them as \
             written, a documented intended divergence, \
             [`gff_attribute_escapes_read_raw`](divergences/gff-attribute-escapes-read-raw.md).",
        ],
    );
    assert_page_states(
        "custom-annotation.md",
        &[
            "`feature_id` (its `ID` attribute) and `associated_gene` (its `gene_id` attribute, \
             else `.`), every `%XX` escape decoded.",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &[
            "| A GTF line repeating an attribute key read as one value (`gene_id`, \
             `transcript_id`, `gene_name`, `transcript_biotype`, a version, `protein_id` and the \
             rest; `tag` is a list) | refused: `<file> line <n>: GTF attribute <key> is repeated; \
             Ensembl VEP reads a repeated key as a list`;",
            "| A GTF value whose `%XX` escapes spell bytes that are not UTF-8 (`\"%C3\"`) | \
             refused: `<file> line <n>: GTF attribute <key> \"<value>\" decodes to text that is \
             not UTF-8`;",
            "| A GFF3 value holding a `%XX` escape (`Name=A%2CB`) | decoded (`A,B`); Ensembl \
             VEP's `--gff` source prints the escape, \
             [`gff_attribute_escapes_read_raw`](divergences/gff-attribute-escapes-read-raw.md). \
             A value with a raw second `=` ends before it on both engines |",
        ],
    );
    assert_page_states(
        "divergences/regulatory-rows.md",
        &[
            "That split also ends a value at a second `=`, so `ID=a=b` names the feature `a`, in \
             both engines.",
            "A line with `ID=a=b` names the feature `a` in both engines, the text between the \
             first and the second `=`.",
        ],
    );
}

/// `gff-gtf-transcript-source`: the `--gtf` source reads the ninth column
/// through `GTF::get_attributes` (GTF.pm:116-136), which lists a key's values
/// when the key is written twice after a defined value (:128-131); the source
/// reads the list where it expects one value, so a symbol, biotype or protein
/// identifier prints joined (Utils.pm:185-193 `convert_arrayref`), a version
/// prints the list's address (OutputFactory.pm:1465, :1560-1565), a support
/// level prints no TSL (:1609-1612) and a record or exon whose identifier or
/// parent is the list is lost (GTF.pm:133-141, :156-158). The cache builder
/// refuses the line, a documented intended divergence. The same page states
/// the line-end reading: both engines read a tabix-indexed line without its
/// `\n` or `\r\n`, and Ensembl VEP keeps the blanks padding an unterminated
/// last value where vep-rs reads the padded line as the unpadded one.
#[test]
fn gtf_repeated_attribute_key_divergence_is_documented() {
    assert_page_states(
        "divergences/gtf-repeated-attribute-key.md",
        &[
            "Ensembl VEP 116.2's `--gtf` source reads a key written twice on one line as a list of \
             both values and, wherever it expects one value, prints the list joined by commas, \
             prints the list's memory address, or loses the record.",
            "`vep-cache-builder` refuses the line: `ERROR: [cache-builder] Failed to parse GTF: \
             <file> line <n>: GTF attribute <key> is repeated; Ensembl VEP reads a repeated key as \
             a list`, exit 1, and no cache is written from the file.",
            "`tag` is read as the list it is, and a repeated key no reader takes keeps its last \
             value. A key written bare and then with a value, `gene_id; gene_id \"G\"`, is one \
             value to both engines.",
            "On the record `tmprss3_codon6`, `21 42395401 G>A`, a `missense_variant` on that \
             transcript, Ensembl VEP 116.2 reading the file through `--gtf` writes \
             `SYMBOL=TMPRSS3,OTHER` in every format.",
            "| `transcript_version` | transcript | `Feature` and the HGVSc reference name \
             `ENST00001113783.ARRAY(0xaaaafe585258)`, another address in the next run |",
            "| `protein_id` | CDS | `ENSP=ARRAY(0xaaaadaff3988).1` and \
             `HGVSp=ARRAY(0xaaaadaff3988).1:p.Pro6Leu`; without `--protein_version`, \
             `ENSP=ENSP00000783588,OTHER` and an HGVSp reference name that is again an address |",
            "| `transcript_support_level` | transcript | no `TSL` field, where the control gene's \
             single `\"3\"` prints `TSL=3` |",
            "| `gene_id` or `transcript_id` | gene or transcript | the TMPRSS3 record reads \
             `intergenic_variant` with no transcript row; the warnings file names the parent not \
             found |",
            "The joined values and the address are written without a warning; the lost record \
             alone draws one, naming the parent that was not found.",
            "A third record, `tmprss3_codon60`, `21 42395000 C>T` inside the transcript's third \
             exon, is run beside them for the exon-line shape of the table below.",
            "| `transcript_id` | exon | the exon is lost: the record inside it, `tmprss3_codon60`, \
             reads `downstream_gene_variant` with `DISTANCE=324` in place of `missense_variant` |",
            "| `exon_id`, `ccds_id`, `protein_version` | exon, transcript, CDS | the control \
             file's rows |",
            "Ensembl VEP 116.2 reading the file through `--gtf` keeps the blanks in the symbol and \
             writes `SYMBOL=TMPRSS3 `; vep-rs reads the line without the blanks that pad it and \
             writes `SYMBOL=TMPRSS3`.",
            "A file whose lines end in `\\r\\n` reads as one whose lines end in `\\n` on both \
             engines, in every format, so the line end itself makes no difference. This is a \
             second mechanism, with its own entry, `gxf_last_value_trailing_blanks`, and its own \
             reproduction.",
            "The class, `gxf_last_value_trailing_blanks`, covers a `SYMBOL`, `BIOTYPE` or `ENSP` \
             value, and its JSON key, of a transcript whose GTF or GFF3 line ends in a pair that \
             carries a value, closes without `;` and is followed by blanks. Ensembl VEP's value \
             carries the blanks; its text with the trailing blanks removed must equal vep-rs's. A \
             support level is not covered: Ensembl VEP reads its digits out of the padded value, \
             so both engines print `TSL=3` for a padded `\"3\"`.",
            "Ensembl VEP's tab output carries `TMPRSS3 ` in its `SYMBOL` column, its JSON output \
             `\"gene_symbol\":\"TMPRSS3 \"`, and its VCF output `TMPRSS3%20%20%20`",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "prints a symbol, biotype or protein identifier as the values joined by commas, a \
             version (and a protein identifier under its version flag) as the list's address, no \
             `TSL` for a support level, and loses the transcript or exon whose identifier or \
             parent the list was; vep-rs refuses to build a cache from the file, naming the line \
             and the key | 3 · 0 of 964 table rows · 14 flag sets, 0 corpora | 116.2 |",
            "keep the blanks that pad a line's last attribute value when its pair closes without \
             `;`, so a symbol prints as `TMPRSS3 `, `TMPRSS3%20%20%20` in a VCF; a cache vep-rs \
             builds from the file reads the padded line as the unpadded one and prints `TMPRSS3`, \
             the name the file encodes | 3 · 0 of 964 table rows · 14 flag sets, 0 corpora | \
             116.2 |",
        ],
    );
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "A line repeating a key read as one value is refused, naming the key; Ensembl VEP \
             reads the list: \
             [`gtf_repeated_attribute_key`](divergences/gtf-repeated-attribute-key.md).",
            "A line's end and padding blanks are not read; Ensembl VEP keeps the blanks of an \
             unterminated last value",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &[
            "no `TSL` for a repeated support level, and loses a record whose identifier it was, \
             [`gtf_repeated_attribute_key`](divergences/gtf-repeated-attribute-key.md) |",
        ],
    );
    assert_page_states(
        "regulatory.md",
        &[
            "A byte that is not UTF-8 in a record is read as U+FFFD by `--regulatory_gff` and \
             printed so ([input bytes that are not UTF-8](divergences/non-utf8-input.md)), a type \
             so read skipping its record with a warning, where `vep-cache-builder` refuses the \
             file, naming the line.",
        ],
    );
}

/// `gxf_malformed_coordinate_line_refused`: Ensembl VEP 116.2's `--gff`,
/// `--gtf`, `--regulatory_gff` and `--custom` GFF sources read the records the
/// `tabix` index returns (BaseGXF.pm:211, 222-237; RegFeat.pm:291-325), and
/// the release's own `tabix` (htslib 1.9, `tbx.c:95,105,235-236`) indexes a
/// file past a line whose end no number leads or reads below zero and, at
/// the file's first line, past one whose start no number leads, printing
/// `Failed to parse TBX_GENERIC` and exiting 0, indexes a zero start and an
/// end below its start that is not below zero silently, and reads a
/// zero-padded column as octal, so the model is built
/// without the line or on the span as written. The builder refuses the build
/// at the first kept line with such a coordinate and at any line whose start
/// or end no number leads or whose end is below zero, naming the file, the
/// line, the column and the value; `--regulatory_gff` refuses the record, and
/// such a start or end on a line it skips, when its chromosome is first read;
/// a `--custom` source refuses the index a line at the file's head bloats,
/// reads a start or end as Perl reads it, an `exact` VCF match aside, and
/// lists no window whose end the index build fails; a zero-padded column is
/// read in decimal; a CDS line that sets the translation's start or end in
/// no exon of its transcript, which Ensembl VEP's sources place by its start
/// alone (BaseGXF.pm:544, 740-756), refuses the build at the line, the
/// section a page of its own holds; a documented intended divergence.
#[test]
fn gxf_malformed_coordinate_line_refused_divergence_is_documented() {
    assert_page_states(
        "divergences/gxf-cds-line-outside-its-exons.md",
        &[
            "They attach each CDS line to the exon that holds the line's start coordinate, and \
             take the translation's two ends from the first and the last attached line.",
            "Its start is the 5' end of the first CDS line in transcript order and its end the \
             3' end of the last: the lowest and the highest coding base on the forward strand, \
             the highest and the lowest on the reverse. Each lies inside an exon of the \
             transcript, or the build is refused.",
            "line 12: CDS 44107360-44107434 of transcript ENST00000456705 starts the \
             translation at 44107360, in no exon of the transcript",
            "line 6: CDS 42395324-42395500 of transcript ENST00001113783 starts the translation \
             at 42395500, in no exon of the transcript",
            "line 6: CDS 42395500-42395900 of transcript ENST00001113783 starts the translation \
             at 42395900, in no exon of the transcript",
            "| one more CDS line in an intron beside the released ones | writes the released \
             rows: the line is given to no exon and sets neither end of the translation | \
             matched: builds the cache of the released lines |",
            "A CDS line that leaves its exon between the translation's two ends is not refused, \
             the last two rows of the table.",
            "refuses the GFF3: `<file> line <n>: CDS <start>-<end> of transcript <id> ends the \
             translation at base <offset> of the <length>-base exon <start>-<end>, the last in \
             transcript order that the coding sequence overlaps` |",
            "No build writes a translation offset outside the exon it is counted in, so a \
             further malformed shape either builds offsets inside their exons or is refused.",
        ],
    );
    assert_page_states(
        "divergences/gxf-malformed-coordinate-line-refused.md",
        &[
            "| `vep-cache-builder`: a CDS line that sets the translation's start or end in no \
             exon of its transcript, [section 9](gxf-cds-line-outside-its-exons.md) | refuses \
             the line |",
            "htslib 1.9, indexes a file past a line inside a sequence whose end no number leads \
             or reads below zero, warning and exiting 0.",
            "It does the same past a start no number leads at the file's first line, and fails \
             the build when that line sits on a later sequence.",
            "A start of 0 and an end below its start, not below zero, it indexes silently.",
            "A variant in a CDS whose line was lost is a `5_prime_UTR_variant`; a transcript \
             whose gene line was lost is intergenic; a GTF transcript whose first CDS line was \
             lost has every codon renumbered.",
            "refuses the build at the first kept line whose start or end is not a number, whose \
             start is 0 or whose end is below its start.",
            "A line it skips, for its type or for holding four to eight columns, is read for \
             one thing more: a start or end no number leads, or an end below zero, refuses the \
             build there with the same error.",
            "`--regulatory_gff` refuses a record whose start or end is not a whole number, whose \
             start is 0 or whose end is below its start, naming the record.",
            "which takes blanks, a sign and a `0x` or `0` prefix, and reads on past an empty or \
             blank column into the next. Text no number so leads fails the line, and so does an \
             end that reads below zero (line 142); a start that does is read as the first base \
             (line 99).",
            "At the file's first line the id was never set, and the index names 67.7 million \
             sequences, 541 MB uncompressed.",
            "A line no number leads it skips, with the same message and exit 0, from 1.10 to \
             1.21, and refuses from 1.22.",
            "| `--regulatory_gff` at run time | refuses the record | refuses the line |",
            "At run time `--regulatory_gff` refuses such a skipped line as it does a kept one, \
             where Ensembl VEP loses the records after an end of the class in its block \
             (section 6).",
            "| `--custom` GFF overlay | refuses no line: reads each column as Perl reads it \
             (section 8), and lists no window whose end the index build fails |",
            "| a `--custom` window inside the file whose end no number leads, or a minus sign \
             leads | lists no such window, and loses the window after it for a record that \
             listed one before it | lists no such window, and every window after it; at the \
             file's head the line leaves the index of section 5 |",
            "its `--custom` sources take the number Perl reads, `exact` VCF matches aside.",
            "| `1e2`, an exponent (`+1e2`, `1.5e2`, `2000e-1`) | index reader: the integer that \
             leads, 1 (2000 for `2000e-1`); Perl: the number spelled, 100 (150, 200), so a \
             promoter window runs from 100 for the records the index returns it to | builder and \
             `--regulatory_gff`: as for `100x`; `--custom`: 100, as Perl | this entry; matched on \
             `--custom` |",
            "`--custom`: 100.5, as Perl | this entry; matched on `--custom` |",
            "`--custom`: -9223372036854775808, as Perl | this entry; matched on `--custom` |",
            "takes 32-bit integers, so a window is listed from the low 32 bits of its start | \
             `--custom`: the number written, the window listed from the first base | an Ensembl \
             VEP defect, not replicated |",
            "| `2147483648`, an end from 2^31 up (`9223372036854775808`, `1e19`) | index reader: an \
             end cut to 32 bits, below the window's start (1 for `1e19`); Perl: the number \
             written, of which the compiled `overlap` of the row above takes the low 32 bits, so \
             the window is listed for no record | `--custom`: the number written, held at \
             9223372036854775807 from 2^63 up; the window is listed for the records the reader \
             returns it to | an Ensembl VEP defect, not replicated |",
            "`--custom`: 100, and no warning: Perl's `isn't numeric` line is its own diagnostic, \
             not a row | this entry |",
            "A `--custom` start of `-5` is matched (section 8).",
            "| a `--regulatory_gff` line of a type the reader skips, an `exon` line, whose end \
             is `xyz` or `-9` | loses the records after it in its block, exit 0 (section 6) | \
             refuses the file at run time, `the end of exon record <ID> of chromosome 21 is not \
             a whole number: xyz`, and at the bake |",
            "| a `--custom` GFF window, BED record or VCF row whose start is `-5` | \
             lists it from the sequence's first base to its end, under the number Perl reads | \
             matched: lists the same record under the same number |",
            "`--custom`: -5, as Perl, the window listed from the first base | this entry; \
             matched on `--custom` |",
            "A `--custom` VCF row whose `END` lies below its position is \
             [`tabix_row_with_an_end_below_its_start`](tabix-row-with-an-end-below-its-start.md)'s.",
            "htslib 1.22.1's `tabix` refuses a column no number leads and an end below its start \
             with exit 1, so a file indexed with it meets the refusal before either engine runs.",
            "Ensembl VEP 116.2 reading the file through `--gff` writes `5_prime_UTR_variant` and \
             `HGVSc=ENST00001113783.1:c.-78C>T`, where the released file gives \
             `missense_variant` and `c.17C>T`; vep-rs refuses to build a cache from the file at \
             that CDS line.",
            "`[E::hts_idx_push] Unsorted positions on sequence #1: 44115691 followed by 1` and \
             `tbx_index_build failed`, exits 1 and writes no index | refuses the build at the \
             line |",
            "writing `5_prime_UTR_variant` and `HGVSc=ENST00001113783.1:c.-77C>T` for \
             `tmprss3_codon6`, exit 0; vep-rs refuses the build at the line.",
            "Ensembl VEP 116.2 builds the transcript without those exons and dies under its CDS, \
             `Can't call method \"seq\" on an undefined value at .../Bio/EnsEMBL/Transcript.pm \
             line 828`, exit 255; vep-rs refuses the build at the line.",
            "A `five_prime_UTR` line with the same end, a type the model skips, kills Ensembl \
             VEP's run the same way, and vep-rs refuses the build at that line too.",
            "| start `abc`, the line on a later sequence | no index to read: `tabix` exits 1 | \
             refuses the build at the line |",
            "| a line of four to eight columns whose start or end no number leads, or whose \
             end is below zero | the index build of the three rows above: `tabix` reads no \
             column past the end | refuses the build at the line, which it skips whatever its \
             type |",
            "| end `xyz` or `-9` | loses the lines after it in its block (section 3) | refuses \
             the build at the line |",
            "| a line that opens with a tab, or with white space and a tab, a type no number \
             leads where a start is read | fails the index build's parse, as the start `abc` \
             does | refuses the build at the line, the record's type kept or skipped |",
            "| start `0`; end below start, not below zero; a start a number leads without being \
             one (`42395942x`, `-5`) | the released file's rows on either sequence of a \
             two-sequence file, the line skipped by its type with `Ignoring 'five_prime_UTR' \
             feature_type` | the released file's rows, the line skipped at its type |",
            "| an empty or blank start, the end a number below every later start of its sequence \
             | the released file's rows on either sequence: `tabix` reads the end's number as \
             the start, in silence | the released file's rows, the line skipped at its type |",
            "| an `unconfirmed_transcript` line whose `ID` opens with `gene:` or `transcript:`, \
             the cache built with the release's protein FASTA or core dumps | skips it by its \
             type | keeps it, as every line so named whatever its type, and holds it to the \
             whole rule of section 2: a start of 0 refuses it |",
            "The release's `tabix` fails a line of fewer than four columns as it fails a start \
             no number leads.",
            "| a `##FASTA` directive | a comment: `tabix` skips it and the recipe's `grep` drops \
             it, so through the recipe the records behind a bare one are indexed and Ensembl VEP \
             reads them; indexed as written, a query of the sequence behind it returns nothing | \
             a comment: skips it and reads the records behind it |",
            "| the first line beginning with `>` that a sequence follows: the next line that is \
             neither blank nor a comment exists and has fewer than four columns, or a start or \
             end `tabix` fails; a column ends at a tab, here and in every row, where `tabix` \
             ends one at a NUL byte too; a line longer than 16777216 bytes is a sequence's only \
             where no tab stands in it, its 16777217th byte aside: one in which a tab stands \
             elsewhere opens no section and is refused for its length, as on any line of the \
             file | fails the `>` line's parse or, where it can parse it, places it as a \
             sequence of its own; indexed as written it then fails the line behind it, and by \
             that line the build has exited 1 where the line has fewer than four columns, as a \
             sequence line has; section 4's table gives the outcome of a line of four or more \
             | opens the sequence section: reads nothing from there on, and logs the line |",
            "| a line beginning with `>` that a line of four or more columns follows, its start \
             and end ones `tabix` parses, or that nothing follows | where it can parse it, \
             places it as a sequence of its own, the first table's last row, and indexes the \
             lines behind it; fails its parse otherwise: the build exits 1 or, on the file's \
             first data line, leaves the index of section 1 | reads it as any other line: a \
             record of its type at nine columns, skipped at four to eight where `tabix` parses \
             its start and end, and refused at fewer than four, a record cut short |",
            "| a record inside the sequence section, nine or more columns whose fourth and fifth \
             are whole numbers, neither a comment nor led by `>` | indexed as written, has \
             exited 1 before the record where the line behind the section's `>` line has fewer \
             than four columns, as a sequence line has; section 4's table gives the outcome of \
             a line of four or more | refuses it, annotation the section has ended: `<file> \
             line <n>: GFF3 annotation after the sequence section, which opens at line <m>` |",
            "| any other line of fewer than four columns | fails its parse | refuses it, a record \
             cut short: `<file> line <n>: GFF3 line ends at column <k>, before its start` (`GTF` \
             in a GTF) |",
            "| a comment behind a byte-order mark | no comment to `tabix`: the recipe's `grep` \
             drops it; indexed as written it fails its parse and, on the file's first line, \
             leaves the index of section 1, a query of the first sequence returning nothing \
             where more comments follow | skips it: the mark is taken off the start of every \
             line before the line is read, in the protein FASTA as in a GFF3, a GTF and a \
             regulatory file, so a file with marks builds the cache of the same file without \
             them |",
            "A tab-delimited file whose every data line `tabix` places opens no section.",
            "A byte-order mark before a file's first record names a sequence of its own to \
             `tabix`: Ensembl VEP loses the line, which vep-rs reads with the mark off.",
            "A kept line led by spaces and no tab is lost likewise; vep-rs reads it trimmed.",
            "So is a kept record whose seqid begins with `>`, which vep-rs reads unless a sequence \
             follows.",
            "vep-rs refuses to open the source, since the index decompresses past its \
             268435456-byte bound, and writes no row.",
            "With the enhancer's end rewritten `xyz`, Ensembl VEP 116.2 loses it and the enhancer \
             after it in its block, writing `intergenic_variant` alone for `in_lost` and \
             `in_last`, exit 0; vep-rs refuses the file, naming the record.",
            "ERROR: [cache] reg_start_zero.gff3.gz: the coordinates 0-5014200 of enhancer record \
             ENSR21_BZTD of chromosome 21 are not 1-based and ordered",
            "The 1.9 `tabix` reads a column in base 0, a leading `0` as octal up to the first \
             digit above 7, so the end reads 275.",
            "vep-rs reads the column in decimal, builds the released model and writes \
             `missense_variant` and `c.17C>T`.",
            "| `+7`, a plus sign leads | index reader: 7; Perl: 7 | 7 in every reader | matched |",
            "| `0100`, zero-padded | index reader: 64, octal; Perl: 100 | 100 in every reader | \
             the octal extension |",
        ],
    );
    assert_page_states(
        "cache-setup.md",
        &[
            "on a line the model reads; a start or end no number leads (`abc`, a lone sign, an \
             empty end before a score of `.`) or an end below zero (`-9`), on any line of the \
             model's file | refused:",
            "an empty or blank column is led by the number of the column after it, as `tabix` \
             reads one;",
            "| A line that is no comment and ends before its fourth column, the start | refused: \
             `<file> line <n>: GFF3 line ends at column <k>, before its start` (`GTF` for a GTF)",
            "a blank line is skipped, a GFF3's `##FASTA` directive is a comment, and nothing is \
             read from the first line that begins with `>` and has a sequence behind it, the \
             next line that is neither blank nor a comment having fewer than four columns, or \
             a start or end the index build fails (a line longer than 16777216 bytes, only \
             where no tab stands in it, its 16777217th byte aside), to the file's end, a line \
             the log names; a `>` line before a line of four or more columns whose start and \
             end the index build parses, or with nothing behind it, is read as any other line \
             |",
            "| A record inside a GFF3's sequence section, nine or more columns whose fourth and \
             fifth are whole numbers, neither a comment nor led by `>` | refused: `<file> line \
             <n>: GFF3 annotation after the sequence section, which opens at line <m>`; the \
             section ends the annotation, so the record would go unread |",
            "| A byte-order mark at the start of a line | taken off before the line is read, in \
             a GFF3 or GTF model, a GTF read for its tags, a regulatory file and the protein \
             FASTA: a file with marks builds the cache of the same file without them |",
            "| A line that opens with a tab | read at the columns the file writes, its seqid \
             empty, as the index build reads it:",
            "white space and a tab before a record refuse it for the type alike, the record's \
             type kept or skipped |",
            "| A file from which no transcript is read | refused before the build writes \
             anything: `No transcript was read from <file>` |",
            "Ensembl VEP builds its model without a kept line its `tabix` could not place, or \
             dies, a documented intended divergence, \
             [`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md)",
            "`<file> line <n>: CDS <start>-<end> of transcript <id> starts the translation at \
             <base>, in no exon of the transcript` or `... ends the translation at <base>, in no \
             exon of the transcript`",
            "a CDS line that leaves its exon between the two ends is read as written",
            "Ensembl VEP builds without the line or on its span as written, a documented \
             intended divergence, \
             [`gxf_malformed_coordinate_line_refused`](divergences/gxf-cds-line-outside-its-exons.md)",
        ],
    );
    assert_page_states(
        "custom-annotation.md",
        &[
            "and so is an index `tabix` built past a start that is no number, \
             [`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md).",
        ],
    );
    assert_page_states(
        "regulatory.md",
        &[
            "A record whose start or end is not a whole number, whose start is 0 or whose end is \
             below its start stops the run when its chromosome is first read, after the output \
             file is opened, naming the record and the file. So does any other line whose start \
             or end no number leads or whose end is below zero. A start or end no number leads \
             at the file's first line leaves an index vep-rs refuses as malformed; such a start \
             on a later sequence leaves none. Ensembl VEP reads on past a line its \
             `tabix` could not place and, when the end was the fault, loses the later records \
             of its block, a documented intended divergence, \
             [`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md).",
        ],
    );
    assert_page_states(
        "intended-divergences.md",
        &[
            "[`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md) | \
             Ensembl VEP 116.2's `--gff`, `--gtf`, `--regulatory_gff` and `--custom` GFF sources \
             read the records the release's own `tabix` indexes",
            "and the sources place a CDS line by its start alone, whatever exon its span leaves, \
             so the model is built without the line or on the span as written",
            "at a CDS line that sets the translation's start or end in no exon of its \
             transcript, and at any line whose start or end no number leads",
            "and at any line whose start or end no number leads or whose end is below zero, \
             refuses such a `--regulatory_gff` record or line at run time",
            "decimal; a skipped-type line the index build parses is matched |",
            "| 3 · 0 of 28,494 table rows · 21 flag sets, 0 corpora | 116.2 | [A GFF3 or GTF line \
             whose coordinates are malformed](divergences/gxf-malformed-coordinate-line-refused.md) |",
        ],
    );
}

/// The index's count of the differences outside every class is its entry
/// rows less the sixteen rows of the defects and the representation
/// difference that produce the classes, so a change that adds a row moves the
/// two count sentences, and this pin, with it.
#[test]
fn the_index_counts_its_entries_outside_the_classes() {
    let rows = page_text("intended-divergences.md")
        .lines()
        .filter(|line| line.starts_with("| [`"))
        .count();
    assert_eq!(rows, 264, "entry rows of docs/intended-divergences.md");
    assert_page_states(
        "intended-divergences.md",
        &[
            "Two hundred and forty-eight differences lie outside every class",
            "Each of the two hundred and forty-eight [differences outside the classes]",
        ],
    );
}
