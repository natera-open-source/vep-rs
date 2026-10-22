//! The rows of `tests/parity/matrix.tsv` whose exclusion a public page states:
//! one test per row, each asserting the page's own sentences, so a rewrite that
//! drops the exclusion fails here.

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
