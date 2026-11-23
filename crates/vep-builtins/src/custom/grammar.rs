// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `--custom` argument grammar, read as Ensembl VEP reads it
//! (`AnnotationSourceAdaptor.pm:327-460` and `AnnotationSource/File.pm:143-197`
//! at release/116.3).
//!
//! A string in the keyed form, `file=<f>,format=<fmt>,short_name=<n>,...`, or the
//! positional form, `<file>,<short_name>,<format>,<type>,<coords>,<field>,...`
//! (the keyed form is the one holding an `=`), becomes one [`CustomSpec`]: the
//! file, the short name the output keys carry, the format, the matching rule and
//! the reporting options, each defaulted as Ensembl VEP defaults it. A malformed
//! string is refused with Ensembl VEP's own text, and a string naming what
//! vep-rs does not read (a remote file) with a text of its own. A GFF or GTF
//! transcript model parses to a source of transcripts
//! ([`CustomSpec::is_transcript_source`]). A file holding `###CHR###` expands to one specification per
//! chromosome whose file exists (`expand_chr_placeholder`).

use std::fmt;
use std::path::Path;

use indexmap::IndexMap;

/// The keys the keyed form accepts (`AnnotationSourceAdaptor.pm:332-346`).
pub(crate) const VALID_OPTIONS: [&str; 13] = [
    "file",
    "format",
    "short_name",
    "fields",
    "type",
    "overlap_cutoff",
    "reciprocal",
    "distance",
    "coords",
    "same_type",
    "num_records",
    "summary_stats",
    "gff_type",
];

/// The placeholder a file name carries to stand for each chromosome
/// (`AnnotationSourceAdaptor.pm:438`).
pub(crate) const CHR_PLACEHOLDER: &str = "###CHR###";

/// Records a source lists per variant when `num_records` is not given
/// (`AnnotationSourceAdaptor.pm:409`).
pub(crate) const DEFAULT_NUM_RECORDS: u64 = 50;

/// The `gff_type` a GFF source takes when none is given
/// (`AnnotationSourceAdaptor.pm:392`): the source is a transcript model, whose
/// transcripts are annotated beside the cache's.
pub(crate) const GFF_TYPE_TRANSCRIPT: &str = "transcript";

/// The `gff_type` of a GFF of GENCODE promoter windows, the one mode that
/// overlays its records on the input (`GFF.pm:162-164`).
pub(crate) const GFF_TYPE_GENCODE_PROMOTER: &str = "gencode_promoter";

/// The short name of a `gff_type=gencode_promoter` source without one of its own
/// (`File.pm:154-157`).
pub(crate) const GENCODE_PROMOTER_SHORT_NAME: &str = "GENCODE_PROMOTER";

/// The fields of a `gencode_promoter` source, whatever its string names
/// (`File.pm:172-174`).
pub(crate) const GENCODE_PROMOTER_FIELDS: [&str; 3] = ["type", "feature_id", "associated_gene"];

/// The mode of a GFF source (`gff_type`, `GFF.pm:155-165`).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum GffType {
    /// A transcript model: the file's transcripts are annotated beside the
    /// cache's ([`CustomSpec::is_transcript_source`]).
    Transcript,
    /// GENCODE promoter windows, overlaid on the input as plain records
    /// (`super::gff_overlay`).
    GencodePromoter,
    /// Any other value: Ensembl VEP opens the source and annotates nothing
    /// from it, as vep-rs does.
    Other(String),
}

impl GffType {
    fn parse(value: &str) -> Self {
        match value {
            GFF_TYPE_TRANSCRIPT => GffType::Transcript,
            GFF_TYPE_GENCODE_PROMOTER => GffType::GencodePromoter,
            other => GffType::Other(other.to_string()),
        }
    }
}

/// How a GFF or GTF transcript model is annotated against without another
/// cache in the run: a cache is built from the file first. A run that names a
/// cache reads the file beside it.
#[doc(hidden)]
pub const GFF_GTF_RECIPE: &str =
    "build a cache from the file and the genome, then annotate against it: \
     vep-cache-builder --species <species> --assembly <assembly> --release <N> --gtf <file> \
     --genome-fasta <genome.fa> --output-dir <dir> (--gff3 <file> for a GFF3), then \
     vep -i <in> -o <out> --offline --json_cache <dir> --fasta <genome.fa>; \
     the page docs/cache-from-gtf-or-gff3.md states the route";

/// The formats a source file is read as (`File.pm:96-102`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum CustomFormat {
    Vcf,
    Bed,
    BigWig,
    /// A GFF: a transcript model under `gff_type=transcript`, an overlay of
    /// promoter windows under `gencode_promoter` ([`GffType`]).
    Gff,
    /// A GTF, a transcript model whatever its `gff_type` (`GTF.pm` has no
    /// other mode).
    Gtf,
}

impl CustomFormat {
    /// The format's key in the `--custom` string, lower case.
    pub fn as_str(self) -> &'static str {
        match self {
            CustomFormat::Vcf => "vcf",
            CustomFormat::Bed => "bed",
            CustomFormat::BigWig => "bigwig",
            CustomFormat::Gff => "gff",
            CustomFormat::Gtf => "gtf",
        }
    }
}

/// How a record of the source and an input variant must lie to match
/// (`File.pm:104-109`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MatchType {
    /// The two spans share a base.
    Overlap,
    /// The record lies within the variant.
    Within,
    /// The record surrounds the variant.
    Surrounding,
    /// The two spans coincide (and, for a VCF, the alleles agree).
    Exact,
}

impl MatchType {
    fn parse(name: &str) -> Option<Self> {
        match name {
            "overlap" => Some(MatchType::Overlap),
            "within" => Some(MatchType::Within),
            "surrounding" => Some(MatchType::Surrounding),
            "exact" => Some(MatchType::Exact),
            _ => None,
        }
    }

    /// The type's key in the `--custom` string.
    pub fn as_str(self) -> &'static str {
        match self {
            MatchType::Overlap => "overlap",
            MatchType::Within => "within",
            MatchType::Surrounding => "surrounding",
            MatchType::Exact => "exact",
        }
    }
}

/// A statistic `summary_stats` can ask for over a variant's records
/// (`AnnotationSourceAdaptor.pm:423`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SummaryStat {
    Min,
    Max,
    Mean,
    Count,
    Sum,
}

impl SummaryStat {
    /// The statistics in the order Ensembl VEP's refusal lists them.
    pub(crate) const ALL: [SummaryStat; 5] = [
        SummaryStat::Min,
        SummaryStat::Max,
        SummaryStat::Mean,
        SummaryStat::Count,
        SummaryStat::Sum,
    ];

    fn parse(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|stat| stat.as_str() == name)
    }

    /// The statistic's key, which is also its output suffix.
    pub fn as_str(self) -> &'static str {
        match self {
            SummaryStat::Min => "min",
            SummaryStat::Max => "max",
            SummaryStat::Mean => "mean",
            SummaryStat::Count => "count",
            SummaryStat::Sum => "sum",
        }
    }
}

/// How many matching records a source lists per variant (`num_records`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RecordLimit {
    /// At most this many, then one `...` record (`File.pm:436-440`); a
    /// fractional or negative value is read as its integer part, at least 0,
    /// under the warning [`CustomSpec::num_records_warning`] carries.
    Count(u64),
    /// Every record (`num_records=all`).
    All,
}

/// One `--custom` source as Ensembl VEP configures it, every option resolved
/// to its default where the string left it out.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct CustomSpec {
    /// The file as given; a local path.
    pub file: String,
    /// The key the source's values carry in the output: as given, else the
    /// file's last path segment, or `GENCODE_PROMOTER_SHORT_NAME` for a
    /// `gencode_promoter` GFF (`File.pm:154-157`).
    pub short_name: String,
    pub format: CustomFormat,
    /// `type`, default `overlap`.
    pub match_type: MatchType,
    /// `coords`: report a record by its coordinates in place of its identifier.
    pub report_coords: bool,
    /// `overlap_cutoff`, a percentage, default 0.
    pub overlap_cutoff: f64,
    /// `distance`: the bases either end of a record may lie from the variant's;
    /// `None` when not given.
    pub distance: Option<Distance>,
    /// `same_type`: match only records of the variant's own class.
    pub same_type: bool,
    /// `reciprocal`: the overlap percentage is the variant's share covered,
    /// not the record's.
    pub reciprocal: bool,
    pub num_records: RecordLimit,
    /// The warning a fractional or negative `num_records` earns, naming the
    /// value and the reading; written once the source opens, so a string the
    /// run refuses or never opens leaves no warnings file.
    pub num_records_warning: Option<String>,
    /// `summary_stats`, `None` when not given or `none`.
    pub summary_stats: Option<Vec<SummaryStat>>,
    /// `gff_type`, default `transcript`.
    pub gff_type: GffType,
    /// The fields reported beside the identifier (`fields=A%B`, or the
    /// positional form's trailing items); `None` when none are named. A
    /// `gencode_promoter` GFF's are `GENCODE_PROMOTER_FIELDS` whatever the
    /// string names (`File.pm:172-174`).
    pub fields: Option<Vec<String>>,
}

impl CustomSpec {
    /// Whether the file stands for one file per chromosome, to be expanded
    /// over the run's chromosomes by `expand_chr_placeholder`.
    pub fn has_chr_placeholder(&self) -> bool {
        self.file.contains(CHR_PLACEHOLDER)
    }

    /// Whether the file is a transcript model, read as a source of
    /// transcripts and annotating no record of its own: a GTF, or a GFF in
    /// its default `transcript` mode (`GFF.pm:155-164`; `GTF.pm` inherits
    /// `BaseGXF.pm`'s transcript annotation alone).
    pub fn is_transcript_source(&self) -> bool {
        match self.format {
            CustomFormat::Gtf => true,
            CustomFormat::Gff => self.gff_type == GffType::Transcript,
            CustomFormat::Vcf | CustomFormat::Bed | CustomFormat::BigWig => false,
        }
    }
}

/// `distance` as Perl reads the option's text where a difference is tested
/// against it (`File.pm:547-548`): an integer for digits alone and for an
/// exponent form whose value is whole (`1e3`), each with nothing but blanks
/// around it and inside the unsigned 64-bit range, and a float for every other
/// text, a point without an exponent (`10.0`) and trailing text (`10abc`)
/// among them.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Distance {
    /// A whole number Perl holds as an integer.
    Whole(u64),
    /// Any other number Perl reads from the text (`perl_number`), a number
    /// below zero among them, which admits nothing.
    Float(f64),
}

impl Distance {
    /// The distance Perl reads from `text`.
    pub(crate) fn read(text: &str) -> Self {
        /// 2^64, one past the unsigned 64-bit range.
        const RANGE_END: f64 = 18_446_744_073_709_551_616.0;
        let number = perl_number(text);
        let text = text.trim_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
        let unsigned = text.strip_prefix(['+', '-']).unwrap_or(text);
        let magnitude = if !unsigned.is_empty() && all_digits(unsigned) {
            unsigned.parse::<u64>().ok()
        } else if spells_an_exponent(unsigned) && number.fract() == 0.0 && number.abs() < RANGE_END
        {
            Some(number.abs() as u64)
        } else {
            None
        };
        match magnitude {
            Some(magnitude) if magnitude == 0 || !text.starts_with('-') => {
                Distance::Whole(magnitude)
            }
            _ => Distance::Float(number),
        }
    }

    /// Whether two whole ends `difference` bases apart lie within the
    /// distance: an exact test against an integer, and against any other
    /// number a test in floats, where a difference past 2^53 rounds.
    pub(crate) fn admits(self, difference: u64) -> bool {
        match self {
            Distance::Whole(whole) => difference <= whole,
            Distance::Float(number) => difference as f64 <= number,
        }
    }

    /// Whether two ends `difference` apart, one of them a float, lie within
    /// the distance: Perl tests a float difference in floats, where an
    /// integer distance past 2^53 rounds.
    pub(crate) fn admits_float(self, difference: f64) -> bool {
        match self {
            Distance::Whole(whole) => difference <= whole as f64,
            Distance::Float(number) => difference <= number,
        }
    }
}

/// Whether every byte of `text` is a decimal digit.
fn all_digits(text: &str) -> bool {
    text.bytes().all(|byte| byte.is_ascii_digit())
}

/// Whether `text`, its sign taken off, spells a number with an exponent and
/// nothing else: digits with or without a point, `e` or `E`, an optional sign
/// and digits.
fn spells_an_exponent(text: &str) -> bool {
    let Some((mantissa, exponent)) = text.split_once(['e', 'E']) else {
        return false;
    };
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let exponent = exponent.strip_prefix(['+', '-']).unwrap_or(exponent);
    all_digits(whole)
        && all_digits(fraction)
        && !(whole.is_empty() && fraction.is_empty())
        && !exponent.is_empty()
        && all_digits(exponent)
}

/// A `--custom` string that is refused: Ensembl VEP's text, less its own
/// `ERROR: ` prefix, where it has one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GrammarError {
    message: String,
}

impl GrammarError {
    fn new(message: impl Into<String>) -> Self {
        GrammarError {
            message: message.into(),
        }
    }

    /// Ensembl VEP's `throw` texts for a `--custom` string end with a `LINE:`
    /// line quoting it (`AnnotationSourceAdaptor.pm:373, 376-377, 434`).
    fn with_line(message: impl fmt::Display, argument: &str) -> Self {
        GrammarError::new(format!("{message}\nLINE: --custom {argument}"))
    }
}

impl fmt::Display for GrammarError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for GrammarError {}

/// Perl's truth of a string: false when empty or `0`.
pub(crate) fn truthy(value: &str) -> bool {
    !value.is_empty() && value != "0"
}

/// Perl's number of a string: its leading numeric prefix (a sign, digits, a
/// fraction, an exponent; `inf` and `nan` by their first letters), else 0, as
/// the comparisons of `File.pm:546-566` read `overlap_cutoff` and `distance`.
pub(crate) fn perl_number(value: &str) -> f64 {
    let text = value.trim_start();
    let bytes = text.as_bytes();
    // The index after the run of digits that starts at `from`.
    let digits_end = |from: usize| {
        from + bytes[from..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count()
    };
    // The index after a `+` or `-` at `from`, or `from` when there is none.
    let after_sign = |from: usize| match bytes.get(from) {
        Some(b'+' | b'-') => from + 1,
        _ => from,
    };

    let digits_start = after_sign(0);
    let word = text[digits_start..].to_ascii_lowercase();
    if word.starts_with("inf") {
        return if text.starts_with('-') {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    if word.starts_with("nan") {
        return f64::NAN;
    }
    let mut end = digits_end(digits_start);
    if bytes.get(end) == Some(&b'.') {
        let fraction_end = digits_end(end + 1);
        // `3.` and `.5` are numbers; a lone `.` is not.
        if fraction_end > end + 1 || end > digits_start {
            end = fraction_end;
        }
    }
    if end == digits_start {
        return 0.0;
    }
    if matches!(bytes.get(end), Some(b'e' | b'E')) {
        let exponent_digits = after_sign(end + 1);
        let exponent_end = digits_end(exponent_digits);
        if exponent_end > exponent_digits {
            end = exponent_end;
        }
    }
    text[..end].parse().unwrap_or(0.0)
}

/// Perl's `split` without a limit: the pieces of `text` at `separator`, the
/// trailing empty ones dropped, and no piece at all for an empty string.
pub(crate) fn perl_split(text: &str, separator: char) -> Vec<&str> {
    let mut pieces: Vec<&str> = text.split(separator).collect();
    while pieces.last().is_some_and(|piece| piece.is_empty()) {
        pieces.pop();
    }
    pieces
}

/// The short name a source takes when none is given (`File.pm:154-157`).
fn default_short_name(file: &str, gff_type: &GffType) -> String {
    if *gff_type == GffType::GencodePromoter {
        return GENCODE_PROMOTER_SHORT_NAME.to_string();
    }
    perl_split(file, '/')
        .last()
        .map_or_else(String::new, |segment| segment.to_string())
}

/// The count a `num_records` value other than `all` sets for `file`, and the
/// warning it earns: the Perl number of the text, its integer part, at least 0.
/// Ensembl VEP keeps the text and compares the record count with it as a
/// number (`File.pm:337, 436-439`), so a fractional or negative value lists
/// records there as no whole number does; here it is read as that count and
/// the warning names the value and the reading. `inf` and `nan` earn none.
fn record_count(value: &str, file: &str) -> (u64, Option<String>) {
    let number = perl_number(value);
    let count = number.max(0.0) as u64;
    let negative = number < 0.0;
    let fractional = number.is_finite() && number.fract() != 0.0;
    let warning = (negative || fractional).then(|| {
        format!(
            "num_records={value} for {file} is read as {count}: \
             the limit is a whole number of at least 0"
        )
    });
    (count, warning)
}

/// Whether `file` names a remote resource: a URL scheme before `://`.
#[doc(hidden)]
pub fn is_remote(file: &str) -> bool {
    let Some((scheme, _)) = file.split_once("://") else {
        return false;
    };
    let mut chars = scheme.chars();
    chars
        .next()
        .is_some_and(|first| first.is_ascii_alphabetic())
        && chars.all(|c| c.is_ascii_alphanumeric() || matches!(c, '+' | '-' | '.'))
}

/// A `--custom` string split into its options: the keyed form's pairs in the
/// order given, or the positional form's five leading items under their keys
/// with the rest as `fields` (`AnnotationSourceAdaptor.pm:353-357`).
struct Options<'a> {
    by_key: IndexMap<&'a str, &'a str>,
    positional_fields: Vec<&'a str>,
}

impl<'a> Options<'a> {
    fn parse(argument: &'a str) -> Result<Self, GrammarError> {
        let mut by_key: IndexMap<&str, &str> = IndexMap::new();
        let mut positional_fields = Vec::new();
        if argument.contains('=') {
            for param in perl_split(argument, ',') {
                // Ensembl VEP's `my ($key, $val) = split('=', $param)` keeps at
                // most the two first `=`-separated pieces.
                let mut pieces = param.split('=');
                let key = pieces.next().unwrap_or_default();
                let Some(value) = pieces.next() else {
                    return Err(GrammarError::new(format!(
                        "Failed to parse parameter {param}; Please add <VALUE_OF_PARAMETER>={param}"
                    )));
                };
                by_key.insert(key, value);
            }
            let invalid: Vec<&str> = by_key
                .keys()
                .copied()
                .filter(|key| !VALID_OPTIONS.contains(key))
                .collect();
            if !invalid.is_empty() {
                return Err(GrammarError::with_line(
                    format!(
                        "The following options are not supported for custom annotations: {}",
                        invalid.join(", ")
                    ),
                    argument,
                ));
            }
        } else {
            let items = perl_split(argument, ',');
            for (key, item) in ["file", "short_name", "format", "type", "coords"]
                .into_iter()
                .zip(items.iter())
            {
                by_key.insert(key, item);
            }
            positional_fields = items.get(5..).unwrap_or_default().to_vec();
        }
        Ok(Options {
            by_key,
            positional_fields,
        })
    }

    fn get(&self, key: &str) -> Option<&'a str> {
        self.by_key.get(key).copied()
    }

    /// The fields named: the keyed form's `fields` split at `%` when it is true
    /// (`:400`), else the positional form's trailing items when there are any
    /// (`:401`).
    fn fields(&self) -> Option<Vec<String>> {
        let pieces = match self.get("fields").filter(|value| truthy(value)) {
            Some(value) => perl_split(value, '%'),
            None => self.positional_fields.clone(),
        };
        (!pieces.is_empty()).then(|| pieces.into_iter().map(str::to_string).collect())
    }
}

/// Ensembl VEP's text for a `--custom` string without a `format`
/// (`AnnotationSourceAdaptor.pm:377`).
const NO_FORMAT: &str = "No format specified for custom annotation source.";

/// Parses one `--custom` string. A string with a `###CHR###` file comes back
/// with the placeholder in place; `expand_chr_placeholder` resolves it.
pub fn parse_custom(argument: &str) -> Result<CustomSpec, GrammarError> {
    let options = Options::parse(argument)?;
    let option = |key: &str| options.get(key);

    let Some(file) = option("file") else {
        return Err(GrammarError::with_line(
            "No file was added for custom annotation source.",
            argument,
        ));
    };
    let Some(format) = option("format") else {
        return Err(GrammarError::with_line(NO_FORMAT, argument));
    };
    if is_remote(file) {
        return Err(GrammarError::new(format!(
            "remote files are not read; copy {file} to a local path"
        )));
    }

    let match_type_name = option("type")
        .filter(|value| truthy(value))
        .unwrap_or("overlap");
    let gff_type = GffType::parse(
        option("gff_type")
            .filter(|value| truthy(value))
            .unwrap_or(GFF_TYPE_TRANSCRIPT),
    );
    let (num_records, num_records_warning) = match option("num_records") {
        None => (RecordLimit::Count(DEFAULT_NUM_RECORDS), None),
        Some("all") => (RecordLimit::All, None),
        Some(value) => {
            let (count, warning) = record_count(value, file);
            (RecordLimit::Count(count), warning)
        }
    };
    let fields = if gff_type == GffType::GencodePromoter {
        Some(GENCODE_PROMOTER_FIELDS.map(String::from).to_vec())
    } else {
        options.fields()
    };

    let summary_stats = match option("summary_stats").filter(|value| truthy(value)) {
        None | Some("none") => None,
        Some(value) => {
            let names = perl_split(value, '%');
            let invalid: Vec<&str> = names
                .iter()
                .copied()
                .filter(|name| SummaryStat::parse(name).is_none())
                .collect();
            if !invalid.is_empty() {
                let valid: Vec<&str> = SummaryStat::ALL.iter().map(|s| s.as_str()).collect();
                return Err(GrammarError::with_line(
                    format!(
                        "The following summary statistics for custom annotations are not supported: {}. Available options are: {}.",
                        invalid.join(", "),
                        valid.join(", ")
                    ),
                    argument,
                ));
            }
            Some(names.into_iter().filter_map(SummaryStat::parse).collect())
        }
    };

    // `File.pm:149-197`: the source's own checks, in its order.
    if !truthy(file) {
        return Err(GrammarError::new("No file given"));
    }
    let short_name = option("short_name")
        .filter(|value| truthy(value))
        .map_or_else(|| default_short_name(file, &gff_type), str::to_string);
    let Some(match_type) = MatchType::parse(match_type_name) else {
        return Err(GrammarError::new(format!(
            "New type \"{match_type_name}\" is not valid"
        )));
    };
    if !truthy(format) {
        return Err(GrammarError::with_line(NO_FORMAT, argument));
    }
    let format = match format.to_ascii_lowercase().as_str() {
        "vcf" => CustomFormat::Vcf,
        "bed" => CustomFormat::Bed,
        "bigwig" => CustomFormat::BigWig,
        "gff" => CustomFormat::Gff,
        "gtf" => CustomFormat::Gtf,
        other => {
            return Err(GrammarError::new(format!(
                "Unknown or unsupported format {other}"
            )));
        }
    };
    // `BaseGXF.pm:109`: a GFF or GTF source never matches exactly.
    if matches!(format, CustomFormat::Gff | CustomFormat::Gtf) && match_type == MatchType::Exact {
        return Err(GrammarError::new(
            "GXF annotation sources cannot be used as \"exact\" custom annotation type",
        ));
    }

    Ok(CustomSpec {
        file: file.to_string(),
        short_name,
        format,
        match_type,
        report_coords: option("coords").is_some_and(truthy),
        overlap_cutoff: option("overlap_cutoff")
            .filter(|value| truthy(value))
            .map_or(0.0, perl_number),
        distance: option("distance").map(Distance::read),
        same_type: option("same_type").is_some_and(truthy),
        reciprocal: option("reciprocal").is_some_and(truthy),
        num_records,
        num_records_warning,
        summary_stats,
        gff_type,
        fields,
    })
}

/// The chromosomes a `###CHR###` file stands for when the run has no reference
/// sequence names (`AnnotationSourceAdaptor.pm:440`).
#[doc(hidden)]
pub fn default_chromosomes() -> Vec<String> {
    (1..=22)
        .map(|n| n.to_string())
        .chain(["X", "Y", "MT"].into_iter().map(String::from))
        .collect()
}

/// One specification per chromosome of `chromosomes` whose file, the
/// placeholder's first occurrence replaced by the chromosome, exists
/// (`AnnotationSourceAdaptor.pm:438-456`); a specification without the
/// placeholder comes back as itself. A short name the string did not give
/// follows each expanded file. No file at all is refused with Ensembl VEP's
/// text, less its `Error: ` prefix.
pub(crate) fn expand_chr_placeholder(
    spec: CustomSpec,
    chromosomes: &[String],
) -> Result<Vec<CustomSpec>, GrammarError> {
    if !spec.has_chr_placeholder() {
        return Ok(vec![spec]);
    }
    let short_name_defaulted = spec.short_name == default_short_name(&spec.file, &spec.gff_type);
    let expanded: Vec<CustomSpec> = chromosomes
        .iter()
        .map(|chr| spec.file.replacen(CHR_PLACEHOLDER, chr, 1))
        .filter(|file| Path::new(file).exists())
        .map(|file| {
            let short_name = if short_name_defaulted {
                default_short_name(&file, &spec.gff_type)
            } else {
                spec.short_name.clone()
            };
            CustomSpec {
                file,
                short_name,
                ..spec.clone()
            }
        })
        .collect();
    if expanded.is_empty() {
        return Err(GrammarError::new(format!(
            "No files with pattern {} were found",
            spec.file
        )));
    }
    Ok(expanded)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(argument: &str) -> CustomSpec {
        parse_custom(argument).unwrap_or_else(|e| panic!("{argument}: {e}"))
    }

    fn error(argument: &str) -> String {
        parse_custom(argument)
            .map(|spec| panic!("{argument} parsed: {spec:?}"))
            .unwrap_err()
            .to_string()
    }

    /// `file=a.vcf.gz,short_name=X,format=vcf,type=exact,fields=A%B` and
    /// `a.vcf.gz,X,vcf,exact,0,A,B` build the same source
    /// (`AnnotationSourceAdaptor.pm:353-357`, `:400-401`), and a string naming
    /// the file and the format alone takes every default of `:380-412`: type
    /// `overlap`, coords 0, cutoff 0, no distance, same_type 0, reciprocal 0,
    /// 50 records, `gff_type` `transcript`, no fields, no statistics, and the
    /// file's last path segment as its short name (`File.pm:157`).
    #[test]
    fn keyed_and_positional_forms_parse_to_the_same_source() {
        let keyed = spec("file=a.vcf.gz,short_name=X,format=vcf,type=exact,fields=A%B");
        let positional = spec("a.vcf.gz,X,vcf,exact,0,A,B");
        assert_eq!(keyed, positional);
        assert_eq!(keyed.short_name, "X");
        assert_eq!(keyed.match_type, MatchType::Exact);
        assert_eq!(keyed.fields, Some(vec!["A".to_string(), "B".to_string()]));

        let defaults = spec("file=data/a.vcf.gz,format=VCF");
        assert_eq!(
            defaults,
            CustomSpec {
                file: "data/a.vcf.gz".to_string(),
                short_name: "a.vcf.gz".to_string(),
                format: CustomFormat::Vcf,
                match_type: MatchType::Overlap,
                report_coords: false,
                overlap_cutoff: 0.0,
                distance: None,
                same_type: false,
                reciprocal: false,
                num_records: RecordLimit::Count(50),
                num_records_warning: None,
                summary_stats: None,
                gff_type: GffType::Transcript,
                fields: None,
            }
        );
        assert_eq!(spec("data/a.bed.gz,,bed").short_name, "a.bed.gz");

        let full = spec(
            "file=a.bed.gz,format=bed,type=within,coords=1,overlap_cutoff=50,distance=10,\
             same_type=1,reciprocal=1,num_records=all,summary_stats=min%max",
        );
        assert_eq!(full.match_type, MatchType::Within);
        assert!(full.report_coords);
        assert_eq!(full.overlap_cutoff, 50.0);
        assert_eq!(full.distance, Some(Distance::Whole(10)));
        assert!(full.same_type && full.reciprocal);
        assert_eq!(full.num_records, RecordLimit::All);
        assert_eq!(
            full.summary_stats,
            Some(vec![SummaryStat::Min, SummaryStat::Max])
        );
        assert_eq!(
            spec("file=a.bed.gz,format=bed,num_records=0").num_records,
            RecordLimit::Count(0)
        );
        // A later key overwrites an earlier one, as a Perl hash assignment does.
        assert_eq!(
            spec("file=a.bed.gz,file=b.bed.gz,format=bed").file,
            "b.bed.gz"
        );
    }

    /// A key outside the thirteen is refused with the text of
    /// `AnnotationSourceAdaptor.pm:371-373`, every such key listed and the
    /// string quoted on the `LINE:` line; a parameter without `=` with the text
    /// of `:361`.
    #[test]
    fn unknown_key_fails_with_vep_text() {
        assert_eq!(
            error("file=a.vcf.gz,format=vcf,name=X,colour=red"),
            "The following options are not supported for custom annotations: name, colour\n\
             LINE: --custom file=a.vcf.gz,format=vcf,name=X,colour=red"
        );
        assert_eq!(
            error("file=a.vcf.gz,vcf"),
            "Failed to parse parameter vcf; Please add <VALUE_OF_PARAMETER>=vcf"
        );
    }

    /// The keyed form without `file` or without `format` is refused with the
    /// texts of `AnnotationSourceAdaptor.pm:376-377`; the positional form with
    /// fewer than three items the same way, since its missing items are
    /// undefined; an empty file with `File.pm:149`'s text and an unknown
    /// format with `:181`'s, lower-cased; an invalid type with `:260`'s.
    #[test]
    fn missing_file_or_format_fails_with_vep_text() {
        assert_eq!(
            error("format=vcf,short_name=X"),
            "No file was added for custom annotation source.\nLINE: --custom format=vcf,short_name=X"
        );
        assert_eq!(
            error("file=a.vcf.gz,short_name=X"),
            "No format specified for custom annotation source.\nLINE: --custom file=a.vcf.gz,short_name=X"
        );
        assert_eq!(
            error("a.vcf.gz,X"),
            "No format specified for custom annotation source.\nLINE: --custom a.vcf.gz,X"
        );
        assert_eq!(
            error(""),
            "No file was added for custom annotation source.\nLINE: --custom "
        );
        assert_eq!(error("file=,format=vcf"), "No file given");
        assert_eq!(
            error("file=a.txt.gz,format=TSV"),
            "Unknown or unsupported format tsv"
        );
        assert_eq!(
            error("file=a.vcf.gz,format=vcf,type=Exact"),
            "New type \"Exact\" is not valid"
        );
    }

    /// A statistic outside the five is refused with the text of
    /// `AnnotationSourceAdaptor.pm:431-434`, the invalid names in the order
    /// given; `summary_stats=none` asks for none (`:415-416`).
    #[test]
    fn invalid_summary_stat_fails_with_vep_text() {
        assert_eq!(
            error("file=a.bed.gz,format=bed,summary_stats=median%max%mode"),
            "The following summary statistics for custom annotations are not supported: median, mode. \
             Available options are: min, max, mean, count, sum.\n\
             LINE: --custom file=a.bed.gz,format=bed,summary_stats=median%max%mode"
        );
        assert_eq!(
            spec("file=a.bed.gz,format=bed,summary_stats=none").summary_stats,
            None
        );
        assert_eq!(
            spec("file=a.bed.gz,format=bed,summary_stats=count").summary_stats,
            Some(vec![SummaryStat::Count])
        );
    }

    /// A `###CHR###` file becomes one source per chromosome whose file exists,
    /// in the order the chromosomes are given, each named after its own file
    /// (`AnnotationSourceAdaptor.pm:438-456`, `File.pm:157`); a given short name
    /// is every source's; no file at all is the text of `:452`.
    #[test]
    fn chr_placeholder_expands_over_existing_files_and_dies_on_none() {
        let dir = tempfile::tempdir().unwrap();
        for chr in ["1", "2"] {
            std::fs::write(dir.path().join(format!("{chr}.vcf.gz")), b"").unwrap();
        }
        let pattern = dir
            .path()
            .join("###CHR###.vcf.gz")
            .to_str()
            .unwrap()
            .to_string();
        let chromosomes: Vec<String> = ["1", "2", "3"].into_iter().map(String::from).collect();

        let spec = parse_custom(&format!("file={pattern},format=vcf")).unwrap();
        assert_eq!(spec.short_name, "###CHR###.vcf.gz");
        let expanded = expand_chr_placeholder(spec, &chromosomes).unwrap();
        let files: Vec<&str> = expanded.iter().map(|s| s.file.as_str()).collect();
        assert_eq!(
            files,
            [
                dir.path().join("1.vcf.gz").to_str().unwrap(),
                dir.path().join("2.vcf.gz").to_str().unwrap(),
            ]
        );
        let names: Vec<&str> = expanded.iter().map(|s| s.short_name.as_str()).collect();
        assert_eq!(names, ["1.vcf.gz", "2.vcf.gz"]);

        let named = parse_custom(&format!("file={pattern},format=vcf,short_name=X")).unwrap();
        let expanded = expand_chr_placeholder(named, &chromosomes).unwrap();
        assert!(expanded.iter().all(|s| s.short_name == "X"));

        let plain = parse_custom("file=a.vcf.gz,format=vcf").unwrap();
        assert_eq!(
            expand_chr_placeholder(plain.clone(), &chromosomes).unwrap(),
            vec![plain]
        );

        let none = parse_custom(&format!("file={pattern},format=vcf")).unwrap();
        assert_eq!(
            expand_chr_placeholder(none, &["4".to_string()])
                .unwrap_err()
                .to_string(),
            format!("No files with pattern {pattern} were found")
        );
        assert_eq!(default_chromosomes().len(), 25);
        assert_eq!(default_chromosomes()[22..], ["X", "Y", "MT"]);
    }

    /// A file behind a URL scheme is refused, once the file and the format are
    /// present and before the type is read. Ensembl VEP refuses an `http` or
    /// `ftp` one only under a `no_remote` option that no flag sets
    /// (`AnnotationSourceAdaptor.pm:378`), so it reads the file over the network.
    #[test]
    fn remote_url_is_refused() {
        for url in [
            "http://example.org/a.vcf.gz",
            "ftp://example.org/a.vcf.gz",
            "https://example.org/a.vcf.gz",
        ] {
            assert_eq!(
                error(&format!("file={url},format=vcf,type=bogus")),
                format!("remote files are not read; copy {url} to a local path")
            );
        }
        assert_eq!(spec("/data/a.vcf.gz,A,vcf").file, "/data/a.vcf.gz");
    }

    /// A GTF, and a GFF in its default `transcript` mode, are transcript models
    /// (`GFF.pm:159-161`, `GTF.pm` through `BaseGXF.pm`): each parses to a
    /// source of transcripts under its short name, the file's last path
    /// segment when none is given (`File.pm:154-157`). Either is refused as
    /// `exact`, as `BaseGXF.pm:109` refuses it. A GFF in `gencode_promoter`
    /// mode is an overlay whose fields are `type`, `feature_id` and
    /// `associated_gene` whatever the string names (`File.pm:172-174`), and a
    /// mode Ensembl VEP does not know is kept as given (`:168`).
    #[test]
    fn gff_and_gtf_transcript_models_parse_as_transcript_sources() {
        let gtf = spec("file=dir/a.gtf.gz,format=gtf");
        assert_eq!(gtf.format, CustomFormat::Gtf);
        assert_eq!(gtf.gff_type, GffType::Transcript);
        assert_eq!(gtf.short_name, "a.gtf.gz");
        assert_eq!(gtf.fields, None);
        assert!(gtf.is_transcript_source());
        let gff = spec("file=a.gff.gz,format=GFF,gff_type=transcript,short_name=models");
        assert_eq!(gff.format, CustomFormat::Gff);
        assert_eq!(gff.short_name, "models");
        assert!(gff.is_transcript_source());
        assert!(spec("a.gff.gz,models,gff").is_transcript_source());
        // `GTF.pm` knows no overlay mode: the file is a transcript model under
        // any `gff_type`, with the name and fields `File.pm:154-174` gives it.
        let named_promoter = spec("file=a.gtf.gz,format=gtf,gff_type=gencode_promoter");
        assert!(named_promoter.is_transcript_source());
        assert_eq!(named_promoter.short_name, "GENCODE_PROMOTER");
        for argument in [
            "file=a.gtf.gz,format=gtf,type=exact",
            "file=a.gff.gz,format=gff,type=exact",
        ] {
            assert_eq!(
                error(argument),
                "GXF annotation sources cannot be used as \"exact\" custom annotation type",
                "{argument}"
            );
        }
        for argument in [
            "file=a.vcf.gz,format=vcf",
            "file=a.bed.gz,format=bed",
            "file=a.bw,format=bigwig",
            "file=a.gff.gz,format=gff,gff_type=gencode_promoter",
            "file=a.gff.gz,format=gff,gff_type=enhancer_window",
        ] {
            assert!(!spec(argument).is_transcript_source(), "{argument}");
        }
        let promoter = spec("file=a.gff.gz,format=gff,gff_type=gencode_promoter");
        assert_eq!(promoter.format, CustomFormat::Gff);
        assert_eq!(promoter.gff_type, GffType::GencodePromoter);
        assert_eq!(promoter.short_name, "GENCODE_PROMOTER");
        let promoter_fields = Some(vec![
            "type".to_string(),
            "feature_id".to_string(),
            "associated_gene".to_string(),
        ]);
        assert_eq!(promoter.fields, promoter_fields);
        assert_eq!(
            spec("file=a.gff.gz,format=gff,gff_type=gencode_promoter,fields=PC%ID").fields,
            promoter_fields
        );
        assert_eq!(
            error("file=a.gff.gz,format=gff,gff_type=gencode_promoter,type=exact"),
            "GXF annotation sources cannot be used as \"exact\" custom annotation type"
        );
        let unknown = spec("file=a.gff.gz,format=gff,gff_type=enhancer_window,fields=ID");
        assert_eq!(
            unknown.gff_type,
            GffType::Other("enhancer_window".to_string())
        );
        assert_eq!(unknown.short_name, "a.gff.gz");
        assert_eq!(unknown.fields, Some(vec!["ID".to_string()]));
    }

    /// Perl's split drops trailing empty pieces and keeps inner ones, and
    /// Perl's truth and number of a string are the ones the options are read
    /// with: `0.0` is true, `50abc` is 50, `abc` is 0, `inf` is infinite.
    #[test]
    fn perl_split_truth_and_number_follow_perl() {
        assert_eq!(perl_split("a,b,,", ','), ["a", "b"]);
        assert_eq!(perl_split("%A%%B", '%'), ["", "A", "", "B"]);
        assert!(perl_split("", ',').is_empty());
        assert!(truthy("0.0") && truthy("00") && !truthy("0") && !truthy(""));
        assert_eq!(perl_number("50abc"), 50.0);
        assert_eq!(perl_number(" 1e2x"), 100.0);
        assert_eq!(perl_number("-2.5"), -2.5);
        assert_eq!(perl_number(".5"), 0.5);
        assert_eq!(perl_number("abc"), 0.0);
        assert_eq!(perl_number("inf"), f64::INFINITY);
        assert_eq!(perl_number("-Infinity"), f64::NEG_INFINITY);
        assert!(perl_number("nan").is_nan());
        assert_eq!(perl_number("e5"), 0.0);
        assert_eq!(perl_number("3e"), 3.0);
        let positional = spec("a.vcf.gz,X,vcf,exact,0,,B");
        assert_eq!(
            positional.fields,
            Some(vec![String::new(), "B".to_string()])
        );
        assert_eq!(
            spec("a.vcf.gz,X,vcf,exact,0,A,").fields,
            Some(vec!["A".to_string()])
        );
        assert_eq!(
            spec("a.vcf.gz,X,vcf,exact,0,0").fields,
            Some(vec!["0".to_string()])
        );
        assert_eq!(spec("a.vcf.gz,X,vcf,,").match_type, MatchType::Overlap);
        assert_eq!(
            spec("file=a.vcf.gz,format=vcf,fields=A%B%")
                .fields
                .unwrap()
                .len(),
            2
        );
        assert_eq!(spec("file=a.vcf.gz,format=vcf,fields=0").fields, None);
        assert_eq!(
            spec("file=a.vcf.gz,format=vcf,distance=").distance,
            Some(Distance::Float(0.0))
        );
        assert_eq!(
            spec("file=a.vcf.gz,format=vcf,num_records=abc").num_records,
            RecordLimit::Count(0)
        );
    }

    /// A `distance` is an integer where Perl holds one for the text: digits
    /// alone, with blanks or a `+` around them, and an exponent form whose
    /// value is whole, up to 2^64 - 1. A point without an exponent, trailing
    /// text, a fraction, a number below zero or past the range and a text
    /// with no number are floats.
    #[test]
    fn a_distance_is_whole_where_perl_holds_an_integer() {
        use Distance::{Float, Whole};
        for (text, read) in [
            ("10", Whole(10)),
            (" +10 ", Whole(10)),
            ("0009007199254740992", Whole(1 << 53)),
            ("18446744073709551615", Whole(u64::MAX)),
            ("1e3", Whole(1000)),
            ("1E+3", Whole(1000)),
            ("1.5e1", Whole(15)),
            ("10.e0", Whole(10)),
            (".1e2", Whole(10)),
            ("9.007199254740992e15", Whole(1 << 53)),
            ("9.223372036854775808e18", Whole(1 << 63)),
            ("-0", Whole(0)),
            ("10.0", Float(10.0)),
            ("10.", Float(10.0)),
            ("10abc", Float(10.0)),
            ("1e3x", Float(1000.0)),
            ("10e", Float(10.0)),
            ("15e-1", Float(1.5)),
            ("-5", Float(-5.0)),
            ("-1e3", Float(-1000.0)),
            ("18446744073709551616", Float(18_446_744_073_709_551_616.0)),
            ("1.8446744073709552e19", Float(18_446_744_073_709_551_616.0)),
            ("", Float(0.0)),
            ("abc", Float(0.0)),
            ("inf", Float(f64::INFINITY)),
        ] {
            assert_eq!(Distance::read(text), read, "{text:?}");
        }
        assert!(matches!(Distance::read("nan"), Float(number) if number.is_nan()));
    }

    /// Ensembl VEP 116.3 keeps a `num_records` text and compares the record
    /// count with it as a number (`File.pm:337`, `:436-439`), so `2.5` lists
    /// three records and `-1` none; vep-rs reads the integer part, at least 0,
    /// and the source carries one warning naming the value and the reading. A
    /// whole number, `all`, `inf`, `nan` and a text with no number carry none.
    #[test]
    fn a_fractional_or_negative_num_records_is_read_as_a_count_with_one_warning() {
        let warning = |value: &str, count: u64| {
            Some(format!(
                "num_records={value} for a.vcf.gz is read as {count}: \
                 the limit is a whole number of at least 0"
            ))
        };
        let cases = [
            ("num_records=2.5", RecordLimit::Count(2), warning("2.5", 2)),
            ("num_records=-1", RecordLimit::Count(0), warning("-1", 0)),
            ("num_records=0.5", RecordLimit::Count(0), warning("0.5", 0)),
            ("num_records=2", RecordLimit::Count(2), None),
            ("num_records=-0", RecordLimit::Count(0), None),
            ("num_records=all", RecordLimit::All, None),
            ("num_records=inf", RecordLimit::Count(u64::MAX), None),
            ("num_records=nan", RecordLimit::Count(0), None),
            ("num_records=abc", RecordLimit::Count(0), None),
        ];
        for (option, num_records, num_records_warning) in cases {
            let spec = spec(&format!("file=a.vcf.gz,format=vcf,{option}"));
            assert_eq!(spec.num_records, num_records, "{option}");
            assert_eq!(spec.num_records_warning, num_records_warning, "{option}");
        }
    }
}
