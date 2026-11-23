// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The UTRAnnotator classes: spans, the gained-start window, the seed sweep and the row classes.

use crate::argparse::path_text;
use crate::consts::{
    CDNA_END_FIELD, HGVSC_FIELDS, HGVS_OFFSET_FIELDS, STRAND_FIELDS, STRUCTURAL_INSERTION_WORD,
    UTRANNOTATOR_ANNOTATION_FIELDS, UTRANNOTATOR_EXISTING_FIELDS,
};
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::{Error, Result};
use crate::fields::{entry_kind, field_named, record_label};
use crate::input::{line_slices, read_text};
use crate::parse::text::{parse_default_under, COMPARISON_CALLERS};
use crate::pyfmt;
use crate::pyjson;
use crate::rules::downstream::field_of;
use crate::rules::record_alleles::{alternates_in_order, record_allele_lengths};
use crate::rules::{nearest_i64, whole_text};
use crate::self_parity::reverse_complement;
use crate::shared::{blank, ABSENT};
use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::LazyLock;

pub const UTRANNOTATOR_GAINED_START_TERM: &str = "5_prime_UTR_premature_start_codon_gain_variant";

/// The position of a reported start codon inside `5UTR_annotation`: `CapDistanceToStart=64`
/// in the text formats, a member of the JSON object.
pub const CAP_DISTANCE_RE_PATTERN: &str =
    r#"CapDistanceToStart"?[\s\x1c-\x1f]*[=:][\s\x1c-\x1f]*"?(-?\d+)"#;
pub static CAP_DISTANCE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(CAP_DISTANCE_RE_PATTERN).expect("the pattern compiles"));

pub const UTRANNOTATOR_CONSEQUENCE_FIELDS: &[&str] = &["5UTR_consequence", "5utr_consequence"];

pub const UTRANNOTATOR_TERMS: &[&str] = &[
    "5_prime_UTR_premature_start_codon_gain_variant",
    "5_prime_UTR_premature_start_codon_loss_variant",
    "5_prime_UTR_uORF_frameshift_variant",
    "5_prime_UTR_uORF_stop_codon_gain_variant",
    "5_prime_UTR_uORF_stop_codon_loss_variant",
];

/// The Allele of a structural allele's row is the Sequence Ontology term of the record's class,
/// a value of `%SO_TERMS` (ensembl-variation Utils/Config.pm:1417-1438); a breakend with a
/// mate writes its own string there instead. The breakpoint class's term is left out: for
/// such a record Parser/VCF.pm:499 puts the REF column before the allele string built from
/// the ALT text unless that string opens or ends with a dot, and the plugin's reference
/// allele is the part that stands first (UTRAnnotator.pm:151-152), the REF column or the
/// first ALT text. The term does not tell the two apart, so its row stays a difference.
pub const STRUCTURAL_CLASS_WORDS: &[&str] = &[
    "Alu_deletion",
    "Alu_insertion",
    "HERV_deletion",
    "HERV_insertion",
    "LINE1_deletion",
    "LINE1_insertion",
    "SVA_deletion",
    "SVA_insertion",
    "copy_number_variation",
    "deletion",
    "duplication",
    "insertion",
    "inversion",
    "mobile_element_deletion",
    "mobile_element_insertion",
    "tandem_duplication",
    "tandem_repeat",
];

/// The class words of a deletion's row among them. A deletion of a stated span fires any of
/// the plugin's analyses but the gained stop, which is a substitution's alone
/// (UTRAnnotator.pm:474), and the frameshift only where the span is no multiple of three
/// (:923). The row of a deletion inside one exon of the 5' UTR names `5_prime_UTR_variant`
/// beside no Consequence term but these two.
///
/// Membership only.
pub const STRUCTURAL_DELETION_WORDS: &[&str] = &[
    "Alu_deletion",
    "HERV_deletion",
    "LINE1_deletion",
    "SVA_deletion",
    "deletion",
    "mobile_element_deletion",
];

pub const UTRANNOTATOR_FRAMESHIFT_TERM: &str = "5_prime_UTR_uORF_frameshift_variant";

/// Membership only.
pub const UTRANNOTATOR_DELETION_TERMS: &[&str] = &[
    "5_prime_UTR_premature_start_codon_gain_variant",
    "5_prime_UTR_premature_start_codon_loss_variant",
    "5_prime_UTR_uORF_frameshift_variant",
    "5_prime_UTR_uORF_stop_codon_loss_variant",
];

pub const STRUCTURAL_DELETION_ROW_TERMS: &[&str] =
    &["NMD_transcript_variant", "feature_truncation"];

/// The class word of a symbolic inversion's row. An inversion replaces its span by as many
/// bases, so it fires any of the plugin's analyses but the frameshift, which returns nothing
/// for two alleles of one length (UTRAnnotator.pm:889).
pub const STRUCTURAL_INVERSION_WORD: &str = "inversion";

/// Membership only.
pub const UTRANNOTATOR_REPLACEMENT_TERMS: &[&str] = &[
    "5_prime_UTR_premature_start_codon_gain_variant",
    "5_prime_UTR_premature_start_codon_loss_variant",
    "5_prime_UTR_uORF_stop_codon_gain_variant",
    "5_prime_UTR_uORF_stop_codon_loss_variant",
];

/// The class word on every row of a record that holds a copy-number loss beside a gain.
pub const STRUCTURAL_COPY_NUMBER_WORD: &str = "copy_number_variation";

/// The Consequence term a row of an inversion or of a copy-number record inside one exon of
/// the 5' UTR may name beside `5_prime_UTR_variant`.
pub const STRUCTURAL_SPAN_ROW_TERMS: &[&str] = &["NMD_transcript_variant"];

/// The ALT texts of a symbolic deletion.
///
/// Membership only.
pub const DELETION_ALTERNATE_TEXTS: &[&str] = &["<CN0>", "<CN=0>", "<DEL>"];

pub const INVERSION_ALTERNATE_TEXT: &str = "<INV>";

pub const UTRANNOTATOR_REFERENCE_UTR_KEYS: &[&str] = &[
    "Evidence",
    "KozakContext",
    "KozakStrength",
    "ref_StartDistanceToCDS",
    "ref_type",
    "ref_type_length",
];

pub const UTRANNOTATOR_ALTERNATE_TYPE: &[&str] = &["alt_type", "uORF"];

pub const UTRANNOTATOR_LOST_START_TERM: &str = "5_prime_UTR_premature_start_codon_loss_variant";

pub const START_CODON: &str = "ATG";

/// An insertion or a duplication between two 5' UTR positions as an HGVSc body writes it.
pub const UTR_HGVSC_INSERTION_RE_PATTERN: &str = r"c\.(-\d+)_(-\d+)ins[A-Za-z]+";
pub static UTR_HGVSC_INSERTION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(UTR_HGVSC_INSERTION_RE_PATTERN).expect("the pattern compiles")
});

pub const UTR_HGVSC_DUPLICATION_RE_PATTERN: &str = r"c\.(-\d+)(?:_(-\d+))?dup";
pub static UTR_HGVSC_DUPLICATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(UTR_HGVSC_DUPLICATION_RE_PATTERN).expect("the pattern compiles")
});

pub const UTRANNOTATOR_LOST_STOP_TERM: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";

pub const UTRANNOTATOR_COUNT_FIELDS: &[&str] = &["Existing_uORFs"];

static CAP_DISTANCE_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(CAP_DISTANCE_RE_PATTERN));

static UTR_HGVSC_INSERTION_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(UTR_HGVSC_INSERTION_RE_PATTERN));

static UTR_HGVSC_DUPLICATION_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(UTR_HGVSC_DUPLICATION_RE_PATTERN));

/// An `HGVS_OFFSET` that is a whole number, with or without a minus sign.
static OFFSET_WHOLE_RE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(r"-?\d+"));

/// The deepest nesting of lists and objects under a group's key that `annotation_frames`
/// writes as text: one level more is the interpreter's `RecursionError`. The interpreter's
/// own bound is what is left of its C stack, so it moves by a few levels from one process
/// to the next: CPython 3.14 under a stack limit of 10240 kB first refuses between 50,195
/// and 50,205 levels, of lists and of objects alike, and this is the middle.
const WRITTEN_MAX_DEPTH: usize = 50_199;

/// `int(text)`, the nearest `i64` for a number past either end. A text the interpreter
/// does not convert ends the run as its `ValueError` does.
fn whole_number(text: &str) -> i64 {
    if let Some(number) = pyfmt::int_of(text) {
        return number;
    }
    match pyfmt::int_checked(text) {
        Ok(decimal) if decimal.starts_with('-') => i64::MIN,
        Ok(_) => i64::MAX,
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// Whether `allele` is one or more of the bases `ACGTN` and nothing else.
fn sequence(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|base| b"ACGTN".contains(&base))
}

/// The entry's `5UTR_consequence` (`ABSENT` for none).
fn term_of(e: &Entry) -> String {
    field_of(Some(e), UTRANNOTATOR_CONSEQUENCE_FIELDS)
}

/// The entry's `5UTR_annotation` (`ABSENT` for none).
fn annotation_of(e: &Entry) -> String {
    field_of(Some(e), UTRANNOTATOR_ANNOTATION_FIELDS)
}

/// Whether the entry carries the plugin's count of the reference 5' UTR's upstream ORFs.
fn counted(e: &Entry) -> bool {
    !blank(&field_of(Some(e), UTRANNOTATOR_COUNT_FIELDS))
}

/// The `CapDistanceToStart` beside `term`: the entry's `5UTR_consequence` is `term` and
/// its `5UTR_annotation` names one position, however often. `None` otherwise.
fn start_position(e: &Entry, term: &str) -> Option<i64> {
    if term_of(e) != term {
        return None;
    }
    let annotation = annotation_of(e);
    let mut positions = CAP_DISTANCE_RE
        .captures_iter(&annotation)
        .filter_map(|found| found.get(1))
        .map(|position| position.as_str());
    let first = positions.next()?;
    positions
        .all(|position| position == first)
        .then(|| whole_number(first))
}

/// The two coordinate texts of a Location, its last one the first where it names one
/// position; `None` where either is not a number.
fn location_bounds(location: &str) -> Option<(&str, &str)> {
    let coordinates = location.split_once(':').map_or(location, |(_, rest)| rest);
    let (start, _, end) = pyfmt::partition(coordinates, "-");
    (pyfmt::is_digit(start) && (end.is_empty() || pyfmt::is_digit(end)))
        .then_some((start, if end.is_empty() { start } else { end }))
}

/// The positions from `start` to `last`, both counted.
fn span_between(start: &str, last: &str) -> i128 {
    let last = whole_number(last);
    i128::from(last) - i128::from(whole_number(start)) + 1
}

/// Whether the entry fills neither `5UTR_consequence` nor `5UTR_annotation`.
fn reports_nothing(e: &Entry) -> bool {
    blank(&term_of(e)) && blank(&annotation_of(e))
}

/// Whether the entry's two fields are ones a class reads beside the other entry's: its
/// `5UTR_consequence` is one of the plugin's five terms, or both of its fields are empty.
fn plugin_fields_fit(e: &Entry) -> bool {
    let term = term_of(e);
    if blank(&term) {
        blank(&annotation_of(e))
    } else {
        UTRANNOTATOR_TERMS.contains(&term.as_str())
    }
}

/// Whether the row's cDNA positions are an insertion's two flanking bases and its allele is
/// a sequence.
fn insertion_shaped(e: &Entry) -> bool {
    cdna_span(e).is_some_and(|(first, last)| i128::from(last) - i128::from(first) == 1)
        && sequence(&e.allele)
}

/// The size of the entry's `HGVS_OFFSET`, 0 for a row that states no whole number there.
fn offset_size(e: &Entry) -> i128 {
    let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
    if OFFSET_WHOLE_RE.is_match(&offset) {
        i128::from(whole_number(&offset)).abs()
    } else {
        0
    }
}

/// Whether the reference allele of the entry's record is empty: what the parser read from
/// the record's own alleles (`Entry.insertion`), else `record_allele_lengths`, a row that
/// states no lengths read as an insertion's.
fn reference_allele_empty(e: &Entry, record_rows: Option<&[&Entry]>) -> bool {
    e.insertion.unwrap_or_else(|| {
        record_allele_lengths(e, record_rows).is_none_or(|(reference, _)| reference == 0)
    })
}

/// The row's first and last cDNA positions, the lower first: the `a-b` or the single
/// position of the text formats' `cDNA_position`, read up to a `/total` suffix, or the
/// JSON `cdna_start` and `cdna_end`, which an insertion writes the higher first.
/// `None` where either is not a number.
pub fn cdna_span(e: &Entry) -> Option<(i64, i64)> {
    fn before_total(position: &str) -> &str {
        pyfmt::partition(position, "/").0
    }
    let (start, end) = match e.fields.get("cDNA_position") {
        Some(position) => {
            let (start, _, end) = pyfmt::partition(before_total(position), "-");
            (start, if end.is_empty() { start } else { end })
        }
        None => {
            let part = |name: &str| e.fields.get(name).map_or(ABSENT, |text| before_total(text));
            (part("cdna_start"), part(CDNA_END_FIELD))
        }
    };
    if !(pyfmt::is_digit(start) && pyfmt::is_digit(end)) {
        return None;
    }
    let (start, end) = (whole_number(start), whole_number(end));
    Some((start.min(end), start.max(end)))
}

/// The `CapDistanceToStart` of the start codon the entry reports as gained: its
/// `5UTR_consequence` is the gained-start term and its `5UTR_annotation` names one
/// position. `None` otherwise.
pub fn gained_start_position(e: &Entry) -> Option<i64> {
    start_position(e, UTRANNOTATOR_GAINED_START_TERM)
}

/// The `CapDistanceToStart` of the start codon the entry reports as lost: its
/// `5UTR_consequence` is the lost-start term and its `5UTR_annotation` names one
/// position. `None` otherwise.
pub fn lost_start_position(e: &Entry) -> Option<i64> {
    start_position(e, UTRANNOTATOR_LOST_START_TERM)
}

/// The entry's sequence allele as its transcript reads it: as written where the row's
/// `STRAND` is `1`, reverse-complemented where it is `-1`. `None` for an allele that
/// is no sequence and for a row that states no strand.
pub fn allele_on_transcript_strand(e: &Entry) -> Option<String> {
    let strand = field_of(Some(e), STRAND_FIELDS);
    if !sequence(&e.allele) {
        return None;
    }
    match strand.as_str() {
        "1" => Some(e.allele.clone()),
        "-1" => Some(reverse_complement(&e.allele)),
        _ => None,
    }
}

/// Whether the two entries of a row differ as the two gained-start windows make them
/// (`utrannotator_gained_start_window`). Both windows open two bases before the allele's
/// first base in the mutated UTR and close two bases after its last: the reference's
/// (UTRAnnotator.pm:345) counts the reference allele's bases, the candidate's the alternate
/// allele's, and each engine reports the first `ATG` of its window by the 0-based index of
/// its `A`. So the row is the class's when one engine alone reports the gained start, its
/// codon begins inside that engine's allele length and at or past the other's, and the
/// reporting engine is the one whose allele is the longer; or when the reference reports it
/// at a negative index, which only a window opened before the UTR's first base yields. An
/// insertion's row is refused where the reference reports a start codon lost and the
/// candidate's gained start stands the allele's length past that codon's index: that codon
/// is the reference's own, moved on by the inserted bases. On an insertion's row the
/// candidate alone may also report the gained start three or more bases before the base
/// after the insertion, where the allele spells the codon there
/// (`start_codon_before_insertion`): the candidate names the first start codon of the
/// changed bases of the mutated sequence's first writing, which lies before the index of a
/// record written further on, and the reference's window at that index (UTRAnnotator.pm:345)
/// does not reach it.
pub fn gained_start_window_row(e: &Entry, a: &Entry, record_rows: Option<&[&Entry]>) -> bool {
    let reference_gains = term_of(e) == UTRANNOTATOR_GAINED_START_TERM;
    let candidate_gains = term_of(a) == UTRANNOTATOR_GAINED_START_TERM;
    let Some(at) = gained_start_position(if reference_gains { e } else { a }) else {
        return false;
    };
    if reference_gains && at < 0 {
        return true;
    }
    let (lengths, span) = (record_allele_lengths(e, record_rows), cdna_span(e));
    if reference_gains == candidate_gains {
        return false;
    }
    let (Some((reference_length, allele_length)), Some((first, last))) = (lengths, span) else {
        return false;
    };
    let at = i128::from(at);
    let (reference_length, allele_length) =
        (i128::from(reference_length), i128::from(allele_length));
    let (first, last) = (i128::from(first), i128::from(last));
    // The cDNA positions are the reference allele's own, an insertion's its two flanking bases.
    let own_positions = if reference_length == 0 {
        2
    } else {
        reference_length
    };
    if last - first + 1 != own_positions {
        return false;
    }
    // A codon the reference reports lost that stands the allele's length further on was moved,
    // not made, by the inserted bases.
    let lost = lost_start_position(e);
    if candidate_gains
        && reference_length == 0
        && lost.is_some_and(|lost| at == i128::from(lost) + allele_length)
    {
        return false;
    }
    let start = if reference_length == 0 {
        first
    } else {
        first - 1
    };
    let (reporting, other) = if reference_gains {
        (reference_length, allele_length)
    } else {
        (allele_length, reference_length)
    };
    if start + other <= at && at < start + reporting {
        return true;
    }
    candidate_gains && reference_length == 0 && start_codon_before_insertion(e, start - at)
}

/// Whether the allele of the insertion's row `e`, read on the transcript's strand
/// (`allele_on_transcript_strand`; a row that states no strand is refused) and on from its
/// last base to its first, spells a start codon from the base that stands `before` bases
/// before its end, at least three: offset `(L - before mod L) mod L` of an allele of `L`
/// bases. An insertion that can be written `before` bases nearer the cap repeats its own
/// bases over that stretch, so those are the three bases the mutated sequence holds
/// `before` bases ahead of the base after the insertion as the record writes it.
fn start_codon_before_insertion(e: &Entry, before: i128) -> bool {
    if before < 3 {
        return false;
    }
    let Some(read) = allele_on_transcript_strand(e) else {
        return false;
    };
    let bases = read.as_bytes();
    let length = bases.len() as i128;
    let from = ((length - before % length) % length) as usize;
    START_CODON
        .bytes()
        .enumerate()
        .all(|(at, base)| bases[(from + at) % bases.len()] == base)
}

/// A text-format `5UTR_annotation` as its groups of `key=value` pairs, the order of
/// the groups and of each group's pairs set aside (UTRAnnotator.pm:282-283 walk both as
/// hashes).
pub fn annotation_groups(value: &str) -> Vec<Vec<String>> {
    let mut groups: Vec<Vec<String>> = value
        .split(',')
        .map(|group| {
            let mut pairs: Vec<String> = group.split(':').map(str::to_string).collect();
            pairs.sort_unstable();
            pairs
        })
        .collect();
    groups.sort_unstable();
    groups
}

/// The `5UTR_annotation` values written beside one `5UTR_consequence` term, each as its
/// groups (`annotation_groups`).
pub type WrittenAnnotations = BTreeSet<Vec<Vec<String>>>;

/// What Ensembl VEP's UTRAnnotator plugin writes on each row of an input across runs
/// under different Perl hash seeds (`PERL_HASH_SEED` with `PERL_PERTURB_KEYS=0`), each
/// of the five flags of UTRAnnotator.pm:262-267 walked last under one of them: per row of
/// the default-format outputs (the Location with the `Uploaded_variation`, the Allele and
/// the Feature), every `5UTR_consequence` term with the `5UTR_annotation` values written
/// beside it. The plugin keeps the last flag that fired, so the terms of a row are the
/// analyses that fire on it.
#[derive(Debug)]
pub struct UtrannotatorSeedSweep {
    pub written: IndexMap<(String, String, String), IndexMap<String, WrittenAnnotations>>,
}

impl UtrannotatorSeedSweep {
    pub fn new() -> Self {
        UtrannotatorSeedSweep {
            written: IndexMap::new(),
        }
    }

    /// The rows of the default-format outputs `paths`, one run each. A file whose header
    /// line does not name the columns a row is keyed by is refused: read as an output, it
    /// would write no term on any row.
    pub fn from_outputs(paths: &[&Path]) -> Result<Self> {
        let mut sweep = Self::new();
        for path in paths {
            let text = read_text(path)?;
            let columns: Vec<&str> = line_slices(&text)
                .find(|line| line.starts_with('#') && !line.starts_with("##"))
                .and_then(|header| header.strip_prefix('#'))
                .unwrap_or("")
                .split('\t')
                .collect();
            if !["Location", "Allele", "Feature"]
                .iter()
                .all(|name| columns.contains(name))
            {
                return Err(Error::comparator(format!(
                    "{}: not a default-format output (columns Location, Allele, Feature)",
                    path_text(path)
                )));
            }
            // A sweep file is parsed one frame nearer the module than a compared output.
            sweep.add_output(&parse_default_under(&text, COMPARISON_CALLERS - 1)?);
        }
        Ok(sweep)
    }

    pub fn row(e: &Entry) -> (String, String, String) {
        (record_label(e), e.allele.clone(), e.feature.clone())
    }

    /// One run's entries: each row's term, where it is one of the plugin's five.
    pub fn add_output(&mut self, entries: &[Entry]) {
        for e in entries {
            let term = term_of(e);
            if UTRANNOTATOR_TERMS.contains(&term.as_str()) {
                let annotation = annotation_groups(&annotation_of(e));
                self.written
                    .entry(Self::row(e))
                    .or_default()
                    .entry(term)
                    .or_default()
                    .insert(annotation);
            }
        }
    }

    /// The one `CapDistanceToStart` the runs write beside the gained-start term on
    /// `e`'s row; `None` where no run writes that term there or they write more than one
    /// position. A VCF or JSON entry names no Location and matches no row.
    pub fn gained_start(&self, e: &Entry) -> Option<i64> {
        e.location.as_ref()?;
        let annotations = self
            .written
            .get(&Self::row(e))?
            .get(UTRANNOTATOR_GAINED_START_TERM)?;
        // A position is told from another by the text of its pair.
        let mut pairs = annotations
            .iter()
            .flatten()
            .flatten()
            .filter_map(|pair| CAP_DISTANCE_WHOLE_RE.captures(pair));
        let first = pairs.next()?;
        if pairs.any(|pair| pair[0] != first[0]) {
            return None;
        }
        Some(whole_number(&first[1]))
    }

    /// Whether the runs write, on `e`'s row, the term of the reference's entry `e`
    /// and the term of the candidate's entry `a`, and under `with_annotation` beside
    /// each the entry's own `5UTR_annotation`. A VCF or JSON entry names no Location and
    /// matches no row.
    pub fn writes(&self, e: &Entry, a: &Entry, with_annotation: bool) -> bool {
        if e.location.is_none() {
            return false;
        }
        let Some(written) = self.written.get(&Self::row(e)) else {
            return false;
        };
        [e, a].into_iter().all(|entry| {
            written.get(&term_of(entry)).is_some_and(|annotations| {
                !with_annotation || annotations.contains(&annotation_groups(&annotation_of(entry)))
            })
        })
    }

    /// Whether a run writes `term` on `e`'s row. A VCF or JSON entry names no Location and
    /// matches no row.
    pub fn writes_term(&self, e: &Entry, term: &str) -> bool {
        e.location.is_some()
            && self
                .written
                .get(&Self::row(e))
                .is_some_and(|written| written.contains_key(term))
    }
}

impl Default for UtrannotatorSeedSweep {
    fn default() -> Self {
        Self::new()
    }
}

/// Whether the two entries' `5UTR_consequence` values are both among the plugin's five
/// terms and differ (`utrannotator_consequence_precedence`).
pub fn utrannotator_precedence_row(e: &Entry, a: &Entry) -> bool {
    // The field under either of its two spellings as written, the text formats' first.
    fn term(entry: &Entry) -> &str {
        UTRANNOTATOR_CONSEQUENCE_FIELDS
            .iter()
            .find_map(|name| entry.fields.get(*name))
            .map_or(ABSENT, String::as_str)
    }
    let (ev, av) = (term(e), term(a));
    ev != av && UTRANNOTATOR_TERMS.contains(&ev) && UTRANNOTATOR_TERMS.contains(&av)
}

/// Whether the reference's entry `e` carries the effect of a structural allele's ALT
/// text read as deleted bases and the candidate's entry `a` none
/// (`utrannotator_structural_allele_text`): the row's Allele is a structural class word on
/// both sides, the reference's `5UTR_consequence` is one of the plugin's five terms, and
/// the candidate fills neither `5UTR_consequence` nor `5UTR_annotation`.
pub fn structural_allele_text_row(e: &Entry, a: &Entry) -> bool {
    if !STRUCTURAL_CLASS_WORDS.contains(&e.allele.as_str()) || a.allele != e.allele {
        return false;
    }
    if !UTRANNOTATOR_TERMS.contains(&term_of(e).as_str()) {
        return false;
    }
    blank(&term_of(a)) && blank(&annotation_of(a))
}

/// Whether the reference's entry `e` carries what the plugin writes in a run without a
/// reference sequence and the candidate's entry `a` nothing of the plugin's
/// (`utrannotator_without_reference`): the candidate's entry names none of the plugin's five
/// fields; the reference fills neither `5UTR_consequence` nor `5UTR_annotation`, or its
/// `5UTR_consequence` is the gained start, the one effect an allele's own bases can fire;
/// and the reference's three `Existing_` counts are each `0`, or all empty on an entry the
/// plugin left alone, which reports no effect. A field the reference's entry carries under
/// two letter cases is read under both.
pub fn without_reference_row(e: &Entry, a: &Entry) -> bool {
    let names = |patterns: &[&str]| a.fields.keys().any(|name| field_named(name, patterns));
    if names(UTRANNOTATOR_CONSEQUENCE_FIELDS)
        || names(UTRANNOTATOR_ANNOTATION_FIELDS)
        || names(UTRANNOTATOR_EXISTING_FIELDS)
    {
        return false;
    }
    let held = |patterns: &[&str]| {
        e.fields
            .iter()
            .filter(|(name, _)| field_named(name, patterns))
            .map(|(_, value)| value.as_str())
            .collect::<Vec<&str>>()
    };
    let terms = held(UTRANNOTATOR_CONSEQUENCE_FIELDS);
    let annotations = held(UTRANNOTATOR_ANNOTATION_FIELDS);
    let silent = terms.iter().chain(&annotations).all(|text| blank(text));
    let gained = !terms.is_empty()
        && terms
            .iter()
            .all(|term| *term == UTRANNOTATOR_GAINED_START_TERM);
    if !silent && !gained {
        return false;
    }
    let counts: Vec<Vec<&str>> = UTRANNOTATOR_EXISTING_FIELDS
        .iter()
        .map(|name| held(&[*name]))
        .collect();
    counts
        .iter()
        .all(|spellings| !spellings.is_empty() && spellings.iter().all(|count| *count == "0"))
        || (silent && counts.iter().flatten().all(|count| blank(count)))
}

/// Whether the candidate's entry `a` carries the effect of a symbolic deletion's stated
/// span on a row where that span lies inside one exon of the 5' UTR
/// (`utrannotator_structural_allele_text`). The row tells it by its own columns: its
/// Allele is a deletion's class word on both sides; its Location's span is the distance
/// between its two cDNA positions, so both ends lie in one exon; its Consequence names
/// `5_prime_UTR_variant` beside no term but `feature_truncation` and
/// `NMD_transcript_variant`, so neither end leaves the 5' UTR. The candidate's
/// `5UTR_consequence` is a term a deletion can fire, the frameshift only where the span is
/// no multiple of three (UTRAnnotator.pm:923), beside a `5UTR_annotation`; the reference's
/// `5UTR_consequence` is one of the plugin's five terms, or both of its fields are empty.
/// A VCF or JSON entry names no Location and is not the class's.
pub fn structural_deletion_span_row(e: &Entry, a: &Entry) -> bool {
    if !STRUCTURAL_DELETION_WORDS.contains(&e.allele.as_str()) {
        return false;
    }
    span_inside_one_utr_exon(e, a, STRUCTURAL_DELETION_ROW_TERMS)
        .is_some_and(|span| deletion_term_fits(a, span))
}

/// The span of a structural row that lies inside one exon of the 5' UTR, as its own columns
/// tell it: both entries carry one Allele; the Location's span is the distance between the
/// two cDNA positions, so both ends lie in one exon; the Consequence names
/// `5_prime_UTR_variant` beside no term but those of `beside`, so neither end leaves the
/// 5' UTR; and the reference's two fields fit (`plugin_fields_fit`). `None` for any other
/// row; a VCF or JSON entry names no Location and is one.
fn span_inside_one_utr_exon(e: &Entry, a: &Entry, beside: &[&str]) -> Option<i128> {
    if a.allele != e.allele {
        return None;
    }
    let location = e.location.as_deref()?;
    let positions = cdna_span(e);
    let ((start, last), (first, second)) = location_bounds(location).zip(positions)?;
    let span = span_between(start, last);
    if span < 1 || i128::from(second) - i128::from(first) + 1 != span {
        return None;
    }
    let mut others = e
        .consequence
        .split(',')
        .filter(|term| !beside.contains(term))
        .peekable();
    if others.peek().is_none() || !others.all(|term| term == "5_prime_UTR_variant") {
        return None;
    }
    plugin_fields_fit(e).then_some(span)
}

/// Whether the candidate's entry `a` fills its two fields with a term the deletion of
/// `span` bases can fire: one of `UTRANNOTATOR_DELETION_TERMS` beside a
/// `5UTR_annotation`, the frameshift only where the span is no multiple of three
/// (UTRAnnotator.pm:923).
fn deletion_term_fits(a: &Entry, span: i128) -> bool {
    let term = term_of(a);
    if !UTRANNOTATOR_DELETION_TERMS.contains(&term.as_str()) || blank(&annotation_of(a)) {
        return false;
    }
    span % 3 != 0 || term != UTRANNOTATOR_FRAMESHIFT_TERM
}

/// The ALT texts of the record of the structural row `e` in the record's own order
/// (`alternates_in_order`), with the place among them of the alternate whose row `e` is. A
/// record writes one row for each of its alternates on a feature, in the alternates' order,
/// so the place is `e`'s own among the record's rows on its feature; no place where those
/// rows are not as many as the alternates or `e` is none of them. `None` where no record
/// is handed with the row or none states its alternates.
fn row_alternate(
    e: &Entry,
    record_rows: Option<&[&Entry]>,
) -> Option<(Vec<String>, Option<usize>)> {
    let rows = record_rows?;
    let alternates = alternates_in_order(e, Some(rows))?;
    let (kind, label) = (entry_kind(e), record_label(e));
    let written: Vec<&Entry> = rows
        .iter()
        .copied()
        .filter(|row| {
            row.feature == e.feature && entry_kind(row) == kind && record_label(row) == label
        })
        .collect();
    let place = if written.len() == alternates.len() {
        written.iter().position(|row| std::ptr::eq(*row, e))
    } else {
        None
    };
    Some((alternates, place))
}

/// Whether the candidate's entry `a` carries the effect of a symbolic inversion's stated
/// span, replaced by its reverse complement, on a row where that span lies inside one exon
/// of the 5' UTR (`utrannotator_structural_allele_text`). The row's Allele is `inversion`
/// on both sides; its Location's span is the distance between its two cDNA positions and
/// its Consequence names `5_prime_UTR_variant` beside no term but
/// `NMD_transcript_variant`. The candidate's `5UTR_consequence` is a term a replacement by
/// as many bases can fire, any of the plugin's five but the frameshift
/// (UTRAnnotator.pm:889), beside a `5UTR_annotation`; the reference's `5UTR_consequence`
/// is one of the plugin's five terms, or both of its fields are empty. Where the row's
/// record is handed and states its alternates, `<INV>` is its one alternate. A VCF or JSON
/// entry names no Location and is not the class's.
pub fn structural_inversion_span_row(e: &Entry, a: &Entry, record_rows: Option<&[&Entry]>) -> bool {
    if e.allele != STRUCTURAL_INVERSION_WORD
        || span_inside_one_utr_exon(e, a, STRUCTURAL_SPAN_ROW_TERMS).is_none()
    {
        return false;
    }
    if !UTRANNOTATOR_REPLACEMENT_TERMS.contains(&term_of(a).as_str()) || blank(&annotation_of(a)) {
        return false;
    }
    row_alternate(e, record_rows)
        .is_none_or(|(alternates, _)| alternates == [INVERSION_ALTERNATE_TEXT])
}

/// Whether the candidate's entry `a` carries, on a row of a record with several
/// alternates, the effect of the deletion of the record's span, which lies inside one exon
/// of the 5' UTR (`utrannotator_structural_allele_text`). The row's Allele is
/// `copy_number_variation` on both sides; its Location's span is the distance between its
/// two cDNA positions and its Consequence names `5_prime_UTR_variant` beside no term but
/// `NMD_transcript_variant`. The candidate's `5UTR_consequence` is a term a deletion can
/// fire, the frameshift only where the span is no multiple of three (UTRAnnotator.pm:923),
/// beside a `5UTR_annotation`; the reference's `5UTR_consequence` is one of the plugin's
/// five terms, or both of its fields are empty. Where the row's record is handed and
/// states its alternates, they are two or more and the row is that of an alternate whose
/// own text is a deletion's (`DELETION_ALTERNATE_TEXTS`). A VCF or JSON entry names no
/// Location and is not the class's.
pub fn structural_copy_number_loss_row(
    e: &Entry,
    a: &Entry,
    record_rows: Option<&[&Entry]>,
) -> bool {
    if e.allele != STRUCTURAL_COPY_NUMBER_WORD {
        return false;
    }
    if !span_inside_one_utr_exon(e, a, STRUCTURAL_SPAN_ROW_TERMS)
        .is_some_and(|span| deletion_term_fits(a, span))
    {
        return false;
    }
    row_alternate(e, record_rows).is_none_or(|(alternates, place)| {
        alternates.len() > 1
            && place.is_some_and(|at| DELETION_ALTERNATE_TEXTS.contains(&alternates[at].as_str()))
    })
}

/// Whether the candidate's entry `a` carries the frameshift a plain symbolic insertion's
/// stated length decides (`utrannotator_structural_allele_text`). The row's Allele is
/// `insertion` on both sides, its Consequence names `5_prime_UTR_variant` and its
/// Location's span, the stated length, is no multiple of three (UTRAnnotator.pm:889, :923).
/// The candidate's `5UTR_consequence` is the frameshift, and each effect of its
/// `5UTR_annotation` holds exactly the six keys read from the reference UTR, or those six
/// and `alt_type=uORF`; the reference's `5UTR_consequence` is one of the plugin's five
/// terms, or both of its fields are empty. The row's columns do not show where the candidate
/// fills it, so the shape is exact by the value it admits. A VCF or JSON entry names no
/// Location and is not the class's.
pub fn structural_insertion_length_row(e: &Entry, a: &Entry) -> bool {
    if e.allele != STRUCTURAL_INSERTION_WORD || a.allele != e.allele {
        return false;
    }
    let Some((start, last)) = e.location.as_deref().and_then(location_bounds) else {
        return false;
    };
    let span = span_between(start, last);
    if span < 1
        || span % 3 == 0
        || !e
            .consequence
            .split(',')
            .any(|term| term == "5_prime_UTR_variant")
    {
        return false;
    }
    if !plugin_fields_fit(e) {
        return false;
    }
    let annotation = annotation_of(a);
    if term_of(a) != UTRANNOTATOR_FRAMESHIFT_TERM || blank(&annotation) {
        return false;
    }
    let (alternate_key, alternate_type) = (
        UTRANNOTATOR_ALTERNATE_TYPE[0],
        UTRANNOTATOR_ALTERNATE_TYPE[1],
    );
    annotation.split(',').all(|effect| {
        let pairs: Vec<(&str, &str)> = effect
            .split(':')
            .map(|pair| {
                let (key, _, value) = pyfmt::partition(pair, "=");
                (key, value)
            })
            .collect();
        let keys: BTreeSet<&str> = pairs.iter().map(|(key, _)| *key).collect();
        let read_from_reference = keys.iter().filter(|key| **key != alternate_key);
        keys.len() == pairs.len()
            && read_from_reference.clone().count() == UTRANNOTATOR_REFERENCE_UTR_KEYS.len()
            && read_from_reference
                .clone()
                .all(|key| UTRANNOTATOR_REFERENCE_UTR_KEYS.contains(key))
            && (keys.len() == UTRANNOTATOR_REFERENCE_UTR_KEYS.len()
                || pairs.contains(&(alternate_key, alternate_type)))
    })
}

/// Whether the reference's entry reports as lost a start codon that an insertion only moves
/// on and the candidate's entry, which carries the plugin's counts, reports no lost start
/// (`utrannotator_lost_start_after_insertion`). An insertion's index in the 5' UTR is that
/// of the base after it (UTRAnnotator.pm:1285-1313), which is the row's lower cDNA position
/// where the two cDNA positions are the insertion's flanking bases. A codon reported at that
/// very index begins directly after the insertion (`:769`) and stands whole whatever the
/// allele. A codon reported one index below it holds the insertion after its first base and
/// stands whole further on when the allele, read on the transcript's strand, ends with that
/// base, `A`; one reported two below holds it after its second base and stands whole when
/// the allele so read ends with `AT` (`allele_on_transcript_strand`; a row that states no
/// strand is refused there). Where the row's record states its alleles (`Entry.insertion`,
/// else `record_allele_lengths`) the reference allele is empty; a text row under an
/// identifier with no record handed is read by its positions and its allele alone.
pub fn lost_start_after_insertion_row(
    e: &Entry,
    a: &Entry,
    record_rows: Option<&[&Entry]>,
) -> bool {
    if term_of(a) == UTRANNOTATOR_LOST_START_TERM || !counted(a) {
        return false;
    }
    let lost = lost_start_position(e);
    let span = cdna_span(e);
    let (Some(lost), Some((first, last))) = (lost, span) else {
        return false;
    };
    if i128::from(last) - i128::from(first) != 1 || !sequence(&e.allele) {
        return false;
    }
    // How many of the codon's bases stand before the insertion.
    let before = i128::from(first) - i128::from(lost);
    if before != 0 {
        let read = allele_on_transcript_strand(e);
        let whole = match before {
            1 | 2 => read.is_some_and(|read| read.ends_with(&START_CODON[..before as usize])),
            _ => false,
        };
        if !whole {
            return false;
        }
    }
    reference_allele_empty(e, record_rows)
}

/// Whether the reference reports as gained, on an insertion's row, a start codon the
/// inserted bases leave at its own index, and the candidate's entry, which carries the
/// plugin's counts, reports no gained start there
/// (`utrannotator_gained_start_remade_in_place`). The row's lower cDNA position is the 5'
/// UTR index of the base after the insertion (UTRAnnotator.pm:1285-1313). A start codon
/// reported one index below it holds the insertion after its first base, and the mutated UTR
/// reads a start codon at that index when the allele, read on the transcript's strand
/// (`allele_on_transcript_strand`; a row that states no strand is refused), opens with the
/// codon's last two bases, `TG`; one reported two below holds it after its second base
/// and the allele so read opens with `G`. The index is the one the reference's own
/// `5UTR_annotation` names beside the gained-start term; where the reference writes another
/// of the plugin's terms, the one index the runs of `utrannotator_sweep` write beside the
/// gained start on the row. A gained start at a negative index, which the reference's
/// window yields by wrapping before the UTR's first base, is no codon of the UTR and is
/// refused. Where the row's record states its alleles (`Entry.insertion`, else
/// `record_allele_lengths`) the reference allele is empty.
pub fn gained_start_remade_in_place_row(
    e: &Entry,
    a: &Entry,
    record_rows: Option<&[&Entry]>,
    utrannotator_sweep: Option<&UtrannotatorSeedSweep>,
) -> bool {
    if !counted(a) {
        return false;
    }
    let span = cdna_span(e);
    let read = allele_on_transcript_strand(e);
    let (Some((first, last)), Some(read)) = (span, read) else {
        return false;
    };
    if i128::from(last) - i128::from(first) != 1 {
        return false;
    }
    let term = term_of(e);
    let at = if term == UTRANNOTATOR_GAINED_START_TERM {
        gained_start_position(e)
    } else if UTRANNOTATOR_TERMS.contains(&term.as_str()) {
        utrannotator_sweep.and_then(|sweep| sweep.gained_start(e))
    } else {
        None
    };
    // A negative index is the wrapped window's (`utrannotator_gained_start_window`), no codon's.
    let Some(at) = at.filter(|at| *at >= 0) else {
        return false;
    };
    // How many of the codon's bases stand before the insertion.
    let before = i128::from(first) - i128::from(at);
    if !matches!(before, 1 | 2) || !read.starts_with(&START_CODON[before as usize..]) {
        return false;
    }
    if term_of(a) == UTRANNOTATOR_GAINED_START_TERM
        && gained_start_position(a).is_none_or(|own| own == at)
    {
        return false;
    }
    reference_allele_empty(e, record_rows)
}

/// The groups of a `5UTR_annotation` as maps of key to value, one group per upstream
/// ORF: the text formats' `key=value` pairs joined by `:`, with `,` between groups
/// (`&` in a VCF value), or the numbered objects of the JSON format's object. `None` for
/// a value that is neither, an empty one included.
///
/// A text the interpreter refuses as no JSON is `None`; one nested past its recursion
/// limit, or holding half a surrogate pair in a string, ends the run.
pub fn annotation_frames(value: &str) -> Option<Vec<Fields>> {
    if value.starts_with('{') {
        let parsed = match pyjson::loads(value) {
            Ok(parsed) => parsed,
            Err(error) if error.is_value_error() => return None,
            Err(error) if error.is_recursion_error() => crate::uncaught!("{error}"),
            Err(_) => crate::uncaught!(
                "a JSON object holds half a surrogate pair: {}...",
                pyfmt::prefix(value, 24)
            ),
        };
        let groups: Vec<&IndexMap<String, pyjson::Value>> = parsed
            .as_dict()?
            .values()
            .map(pyjson::Value::as_dict)
            .collect::<Option<_>>()?;
        return Some(
            groups
                .into_iter()
                .map(|group| {
                    group
                        .iter()
                        .map(|(key, item)| {
                            if pyjson::depth(item) > WRITTEN_MAX_DEPTH {
                                crate::uncaught!(
                                    "Stack overflow (used 10200 kB) while getting the repr of an object"
                                );
                            }
                            (key.clone(), pyjson::py_str(item))
                        })
                        .collect()
                })
                .collect(),
        );
    }
    if blank(value) || (value.contains(',') && value.contains('&')) {
        return None;
    }
    let separator = if value.contains('&') { '&' } else { ',' };
    let groups: Vec<Vec<&str>> = value
        .split(separator)
        .map(|group| group.split(':').collect())
        .collect();
    if !groups.iter().flatten().all(|pair| pair.contains('=')) {
        return None;
    }
    Some(
        groups
            .iter()
            .map(|group| {
                group
                    .iter()
                    .map(|pair| {
                        let (key, _, item) = pyfmt::partition(pair, "=");
                        (key.to_string(), item.to_string())
                    })
                    .collect()
            })
            .collect(),
    )
}

/// The 5' UTR base directly after the insertion of the entry's row as its record writes
/// it, in bases before the coding start, read from the row's HGVSc and `HGVS_OFFSET`. An
/// `ins` between two neighbouring 5' UTR positions names that base as its second
/// coordinate; a `dup` stands for the copy written directly after the duplicated bases, so
/// the base is the one after its last coordinate. The notation stands 3' of the record's own
/// insertion point on the transcript by the size of `HGVS_OFFSET`, which is printed
/// negative on a reverse-strand transcript (OutputFactory.pm:1764), so the record's point
/// lies that many bases further from the coding start. `None` for any other notation and
/// for a coordinate outside the 5' UTR.
pub fn utr_base_after_insertion(e: &Entry) -> Option<i64> {
    utr_bases_after_insertion(e).map(|(written, _)| written)
}

/// The base `utr_base_after_insertion` reads, with the base its notation itself names: the
/// one after the insertion written at the last index that writes the same mutated sequence,
/// nearer the coding start by the size of `HGVS_OFFSET`.
fn utr_bases_after_insertion(e: &Entry) -> Option<(i64, i64)> {
    let hgvsc = field_of(Some(e), HGVSC_FIELDS);
    let notation = pyfmt::partition(&hgvsc, ":").2;
    let shift = offset_size(e);
    let named = if let Some(inserted) = UTR_HGVSC_INSERTION_WHOLE_RE.captures(notation) {
        let first = i128::from(whole_number(&inserted[1]));
        let last = i128::from(whole_number(&inserted[2]));
        if last != first + 1 {
            return None;
        }
        -last
    } else {
        let duplicated = UTR_HGVSC_DUPLICATION_WHOLE_RE.captures(notation)?;
        let last = duplicated.get(2).or_else(|| duplicated.get(1))?.as_str();
        -i128::from(whole_number(last)) - 1
    };
    Some((nearest_i64(shift + named), nearest_i64(named)))
}

/// A `5UTR_annotation`'s groups, each as its pairs in key order, with the count of the
/// groups that hold the same pairs, in the order the groups first stand.
fn counted_groups(frames: &[Fields]) -> IndexMap<Vec<(&str, &str)>, usize> {
    let mut counted = IndexMap::new();
    for frame in frames {
        let mut pairs: Vec<(&str, &str)> = frame
            .iter()
            .map(|(key, item)| (key.as_str(), item.as_str()))
            .collect();
        pairs.sort_unstable();
        *counted.entry(pairs).or_insert(0) += 1;
    }
    counted
}

/// The groups one counted `5UTR_annotation` holds more often than another.
fn groups_beyond<'a, 'b>(
    held: &'a IndexMap<Vec<(&'b str, &'b str)>, usize>,
    other: &IndexMap<Vec<(&'b str, &'b str)>, usize>,
) -> Vec<&'a Vec<(&'b str, &'b str)>> {
    held.iter()
        .filter(|(group, count)| **count > other.get(*group).copied().unwrap_or(0))
        .map(|(group, _)| group)
        .collect()
}

/// The value of `name` among a group's pairs.
fn pair_item<'a>(pairs: &[(&'a str, &'a str)], name: &str) -> Option<&'a str> {
    pairs
        .iter()
        .find(|(key, _)| *key == name)
        .map(|(_, item)| *item)
}

/// Whether the two frameshift lists of an insertion's row differ by upstream ORFs at whose
/// stop codon the insertion is written (`utrannotator_frameshift_at_stop_boundary`).
/// UTRAnnotator.pm:921 holds the record's own index against the check point one index
/// before an ORF's stop codon (`:911-917`) and reads no base of the mutated UTR. An entry
/// lists the groups of its `5UTR_annotation` where its `5UTR_consequence` is the frameshift
/// and none where it fills neither field; at least one entry reports the frameshift, and a
/// candidate's entry that fills neither field carries the plugin's counts. The stop codon's
/// first base of a uORF is its `ref_StartDistanceToCDS` less its `ref_type_length` plus
/// three bases before the coding start; a row that states no insertion point
/// (`utr_base_after_insertion`) is refused. Two shapes, never both on one row.
///
/// The candidate lists more: an insertion written directly before the stop codon's first
/// base is past the check point and the ORF is left out, though it runs on past its stop
/// codon's index and its length changes. Every group the candidate adds is a uORF whose two
/// printed lengths differ and whose stop codon's first base is the base after the insertion
/// as the record writes it, and no group both entries list is such an ORF: there the
/// reference's own test admits what it leaves out here.
///
/// The reference lists more: an insertion written at or before the check point is admitted
/// whatever the mutated UTR reads, so an ORF is listed that ends on a stop codon at its stop
/// codon's own index and keeps its length. Every group the reference adds is a uORF in the
/// mutated UTR too (`alt_type=uORF`) whose `alt_type_length` equals its
/// `ref_type_length` and whose stop codon's first base lies at or after the base after the
/// insertion as the record writes it and no more than one base beyond the base after the
/// insertion as the row's notation names it, the last index that writes the same sequence.
pub fn frameshift_at_stop_boundary_row(e: &Entry, a: &Entry) -> bool {
    let shifted = |entry: &Entry| term_of(entry) == UTRANNOTATOR_FRAMESHIFT_TERM;
    if !(shifted(e) || shifted(a)) || !(shifted(a) || counted(a)) {
        return false;
    }
    let frames = |entry: &Entry| {
        if shifted(entry) {
            annotation_frames(&annotation_of(entry))
        } else {
            reports_nothing(entry).then(Vec::new)
        }
    };
    let (expected, actual) = (frames(e), frames(a));
    let (Some(expected), Some(actual)) = (expected, actual) else {
        return false;
    };
    let (listed, written) = (counted_groups(&expected), counted_groups(&actual));
    let dropped = groups_beyond(&listed, &written);
    let added = groups_beyond(&written, &listed);
    if dropped.is_empty() == added.is_empty() {
        return false;
    }
    let Some((after, named)) = utr_bases_after_insertion(e) else {
        return false;
    };
    let ends_at_the_insertion = |pairs: &[(&str, &str)]| {
        let start = pair_item(pairs, "ref_StartDistanceToCDS").unwrap_or(ABSENT);
        let before = pair_item(pairs, "ref_type_length").unwrap_or(ABSENT);
        if pair_item(pairs, "ref_type") != Some("uORF")
            || !(pyfmt::is_digit(start) && pyfmt::is_digit(before))
        {
            return false;
        }
        pair_item(pairs, "alt_type_length") != Some(before)
            && i128::from(whole_number(start)) - i128::from(whole_number(before)) + 3
                == i128::from(after)
    };
    let keeps_its_length = |pairs: &[(&str, &str)]| {
        let start = pair_item(pairs, "ref_StartDistanceToCDS").unwrap_or(ABSENT);
        let before = pair_item(pairs, "ref_type_length").unwrap_or(ABSENT);
        if pair_item(pairs, "ref_type") != Some("uORF")
            || pair_item(pairs, "alt_type") != Some("uORF")
            || pair_item(pairs, "alt_type_length") != Some(before)
            || !(pyfmt::is_digit(start) && pyfmt::is_digit(before))
        {
            return false;
        }
        let stop = i128::from(whole_number(start)) - i128::from(whole_number(before)) + 3;
        stop <= i128::from(after) && stop >= i128::from(named) - 1
    };
    if added.is_empty() {
        return dropped.iter().all(|pairs| keeps_its_length(pairs));
    }
    added.iter().all(|pairs| ends_at_the_insertion(pairs))
        && !written
            .keys()
            .filter(|group| listed.contains_key(*group))
            .any(|pairs| ends_at_the_insertion(pairs))
}

/// Whether the two entries of an insertion's row differ by a lost stop that only one way of
/// writing the mutated sequence shows the reference
/// (`utrannotator_lost_stop_written_before_stop_codon`). UTRAnnotator.pm:628-633 read a
/// stop codon's index only for a record written one or two bases past the stop codon's
/// first, and `:639-642` then read the three bases at that index alone. The row's cDNA
/// positions are an insertion's two flanking bases, the candidate's entry carries the
/// plugin's counts, and where the row's record states its alleles (`Entry.insertion`,
/// else `record_allele_lengths`) the reference allele is empty. Three shapes.
///
/// The candidate reports the lost stop of upstream ORFs the reference lists as shifted
/// (`lost_stop_of_listed_orfs`).
///
/// The candidate reports a lost stop the reference does not read
/// (`lost_stop_written_further_on`): the record is written at or before the stop codon's
/// first base and the same sequence can be written inside the stop codon.
///
/// The reference reports a lost stop the candidate does not
/// (`stop_codon_whole_in_frame`): the record is written inside the stop codon and the same
/// sequence can be written at its first base, where the stop codon stands whole after the
/// inserted bases, in frame.
pub fn lost_stop_written_before_stop_codon_row(
    e: &Entry,
    a: &Entry,
    record_rows: Option<&[&Entry]>,
    utrannotator_sweep: Option<&UtrannotatorSeedSweep>,
) -> bool {
    let lost = |entry: &Entry| term_of(entry) == UTRANNOTATOR_LOST_STOP_TERM;
    if lost(e) == lost(a) || !counted(a) || !insertion_shaped(e) {
        return false;
    }
    let fits = if lost(e) {
        stop_codon_whole_in_frame(e, a)
    } else {
        lost_stop_of_listed_orfs(e, a) || lost_stop_written_further_on(e, a, utrannotator_sweep)
    };
    fits && reference_allele_empty(e, record_rows)
}

/// The value of `name` in an annotation's group.
fn frame_item<'a>(frame: &'a Fields, name: &str) -> Option<&'a str> {
    frame.get(name).map(String::as_str)
}

/// Whether the reference's entry `e` lists as shifted the upstream ORFs whose lost stop the
/// candidate's entry `a` reports. For the mutated sequence written before the start codon
/// that overlaps a stop codon, UTRAnnotator.pm:921 reports the ORF shifted, with its length
/// in the mutated UTR. The reference's `5UTR_consequence` is the frameshift; each of the
/// candidate's groups carries the `KozakContext` and `KozakStrength` of exactly one uORF of
/// the reference's list, a different one for each group, and the alternative stop codon
/// where that ORF's printed `alt_type_length` ends it: at `ref_StartDistanceToCDS` plus
/// the allele's length less `alt_type_length` plus three bases before the coding start
/// where the mutated UTR closes the ORF (`alt_type=uORF`), none otherwise.
fn lost_stop_of_listed_orfs(e: &Entry, a: &Entry) -> bool {
    if term_of(e) != UTRANNOTATOR_FRAMESHIFT_TERM {
        return false;
    }
    let shifted = annotation_frames(&annotation_of(e));
    let lost = annotation_frames(&annotation_of(a));
    let (Some(shifted), Some(lost)) = (shifted, lost) else {
        return false;
    };
    let item = frame_item;
    let mut taken: BTreeSet<usize> = BTreeSet::new();
    for group in &lost {
        let (Some(context), Some(strength)) =
            (item(group, "KozakContext"), item(group, "KozakStrength"))
        else {
            return false;
        };
        let mut frames = shifted.iter().enumerate().filter(|(_, frame)| {
            item(frame, "ref_type") == Some("uORF")
                && item(frame, "KozakContext") == Some(context)
                && item(frame, "KozakStrength") == Some(strength)
        });
        let (Some((at, frame)), None) = (frames.next(), frames.next()) else {
            return false;
        };
        if !taken.insert(at) {
            return false;
        }
        let start = item(frame, "ref_StartDistanceToCDS").unwrap_or(ABSENT);
        let length = item(frame, "alt_type_length").unwrap_or(ABSENT);
        let stop = (item(group, "AltStop"), item(group, "AltStopDistanceToCDS"));
        if item(frame, "alt_type") == Some("uORF") {
            if !(pyfmt::is_digit(start) && pyfmt::is_digit(length)) {
                return false;
            }
            let ends = i128::from(whole_number(start)) + pyfmt::len(&e.allele) as i128
                - i128::from(whole_number(length))
                + 3;
            if stop != (Some("True"), Some(ends.to_string().as_str())) {
                return false;
            }
        } else if stop != (Some("False"), Some("NA")) {
            return false;
        }
    }
    true
}

/// Whether the candidate's entry `a` reports a lost stop that the reference's entry `e`,
/// written at or before the stop codon's first base, does not read (UTRAnnotator.pm:631).
/// The reference fills neither field or its `5UTR_consequence` is the lost start, the
/// gained start or the frameshift, and no group of a frameshift's list carries the
/// `KozakContext` and `KozakStrength` of a group of the candidate's
/// (`lost_stop_of_listed_orfs` holds that row). The allele's length is no multiple of
/// three, so no writing leaves the ORF its frame, and the row's `HGVS_OFFSET` is a whole
/// number other than zero: the same sequence can be written further on, inside the stop
/// codon. Each group of the candidate's `5UTR_annotation` names an ORF by its
/// `KozakContext` and `KozakStrength`. Where seeded runs are given, none writes the lost
/// stop on the row.
fn lost_stop_written_further_on(
    e: &Entry,
    a: &Entry,
    utrannotator_sweep: Option<&UtrannotatorSeedSweep>,
) -> bool {
    const READ_BEFORE_STOP: [&str; 3] = [
        UTRANNOTATOR_FRAMESHIFT_TERM,
        UTRANNOTATOR_GAINED_START_TERM,
        UTRANNOTATOR_LOST_START_TERM,
    ];
    let reference = term_of(e);
    if !(READ_BEFORE_STOP.contains(&reference.as_str()) || reports_nothing(e)) {
        return false;
    }
    if pyfmt::len(&e.allele).is_multiple_of(3) || offset_size(e) == 0 {
        return false;
    }
    let kozak = |frame: &Fields| {
        frame_item(frame, "KozakContext")
            .zip(frame_item(frame, "KozakStrength"))
            .map(|(context, strength)| (context.to_string(), strength.to_string()))
    };
    let Some(lost) = annotation_frames(&annotation_of(a)) else {
        return false;
    };
    let Some(named) = lost.iter().map(kozak).collect::<Option<BTreeSet<_>>>() else {
        return false;
    };
    if named.is_empty() {
        return false;
    }
    if reference == UTRANNOTATOR_FRAMESHIFT_TERM {
        let Some(shifted) = annotation_frames(&annotation_of(e)) else {
            return false;
        };
        if shifted
            .iter()
            .filter_map(kozak)
            .any(|pair| named.contains(&pair))
        {
            return false;
        }
    }
    utrannotator_sweep.is_none_or(|sweep| !sweep.writes_term(e, UTRANNOTATOR_LOST_STOP_TERM))
}

/// Whether the reference's entry `e` reports as lost a stop codon that stands whole in
/// frame, and the candidate's entry `a` no lost stop. UTRAnnotator.pm:639-642 read the
/// three bases at the stop codon's index for a record written one or two bases past it;
/// the same sequence written at the stop codon's first base leaves the stop codon whole,
/// as many bases on as the allele is long, and an allele whose length is a multiple of
/// three leaves it in the ORF's frame. Each group of the reference's `5UTR_annotation`
/// reads `AltStop=True` with an `AltStopDistanceToCDS` (`:679`) one or two bases more
/// than the distance of the base after the insertion as the record writes it
/// (`utr_base_after_insertion`): the stop codon the ORF meets is its own. A stop codon
/// that the inserted bases put in the ORF's frame before it stands a multiple of three
/// bases further from the coding start, so the distance may be more by any number below
/// the allele's length that is no multiple of three. The candidate's `5UTR_consequence`
/// is another of the plugin's five terms, or both of its fields are empty.
fn stop_codon_whole_in_frame(e: &Entry, a: &Entry) -> bool {
    if !plugin_fields_fit(a) {
        return false;
    }
    let length = pyfmt::len(&e.allele) as i128;
    if length % 3 != 0 {
        return false;
    }
    let Some(after) = utr_base_after_insertion(e) else {
        return false;
    };
    let Some(lost) = annotation_frames(&annotation_of(e)) else {
        return false;
    };
    !lost.is_empty()
        && lost.iter().all(|group| {
            let distance = frame_item(group, "AltStopDistanceToCDS").unwrap_or(ABSENT);
            if frame_item(group, "AltStop") != Some("True") || !pyfmt::is_digit(distance) {
                return false;
            }
            let before = i128::from(whole_number(distance)) - i128::from(after);
            before >= 1 && before < length && before % 3 != 0
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Uncaught;

    /// A default-format row at `21:1000-1001` under the identifier `ins` with `fields`.
    fn row(allele: &str, fields: &[(&str, &str)]) -> Entry {
        let mut held: Fields = Fields::new();
        held.insert("Uploaded_variation".to_string(), "ins".to_string());
        for (name, value) in fields {
            held.insert(name.to_string(), value.to_string());
        }
        Entry {
            key: format!("21:1000-1001|{allele}|ENST1|Transcript"),
            location: Some("21:1000-1001".to_string()),
            allele: allele.to_string(),
            feature: "ENST1".to_string(),
            consequence: "5_prime_UTR_variant".to_string(),
            fields: held,
            ..Entry::default()
        }
    }

    /// What a run says that ends on an input no reading of the comparison expects.
    fn uncaught<T>(call: impl FnOnce() -> T) -> String {
        let payload = std::panic::catch_unwind(std::panic::AssertUnwindSafe(call))
            .err()
            .expect("the call ends the run");
        let Uncaught(what) = payload
            .downcast_ref::<Uncaught>()
            .expect("the run ends as an uncaught exception does");
        what.clone()
    }

    fn frame(pairs: &[(&str, &str)]) -> Fields {
        pairs
            .iter()
            .map(|(key, item)| (key.to_string(), item.to_string()))
            .collect()
    }

    #[test]
    fn the_patterns_hold_a_whole_text_as_they_hold_a_part() {
        assert!(CAP_DISTANCE_WHOLE_RE.is_match("CapDistanceToStart=64"));
        assert!(CAP_DISTANCE_WHOLE_RE.is_match("CapDistanceToStart\": -7"));
        assert!(!CAP_DISTANCE_WHOLE_RE.is_match("CapDistanceToStart=64:type=uORF"));
        assert!(!CAP_DISTANCE_WHOLE_RE.is_match("CapDistanceToStart=64\n"));
        assert!(UTR_HGVSC_INSERTION_WHOLE_RE.is_match("c.-55_-54insA"));
        assert!(!UTR_HGVSC_INSERTION_WHOLE_RE.is_match("c.-55_-54insA "));
        assert!(UTR_HGVSC_DUPLICATION_WHOLE_RE.is_match("c.-56_-55dup"));
        assert!(!UTR_HGVSC_DUPLICATION_WHOLE_RE.is_match("c.-55dupA"));
        assert!(OFFSET_WHOLE_RE.is_match("-2") && !OFFSET_WHOLE_RE.is_match("+2"));
    }

    #[test]
    fn the_cdna_positions_are_read_in_either_spelling_the_lower_first() {
        let span = |fields: &[(&str, &str)]| cdna_span(&row("CAT", fields));
        assert_eq!(span(&[("cDNA_position", "64-65")]), Some((64, 65)));
        assert_eq!(span(&[("cDNA_position", "7")]), Some((7, 7)));
        assert_eq!(span(&[("cDNA_position", "336-340/1200")]), Some((336, 340)));
        assert_eq!(span(&[("cDNA_position", "65-64")]), Some((64, 65)));
        assert_eq!(
            span(&[("cdna_start", "65"), ("cdna_end", "64")]),
            Some((64, 65))
        );
        for unread in ["?-65", "64-?", "-", "", "64-65-66", " 64"] {
            assert_eq!(span(&[("cDNA_position", unread)]), None, "{unread}");
        }
        assert_eq!(span(&[("cdna_start", "65")]), None);
        assert_eq!(span(&[]), None);
        // The text formats' field is the one read where a row holds both spellings.
        assert_eq!(
            span(&[
                ("cdna_start", "1"),
                ("cdna_end", "2"),
                ("cDNA_position", "-")
            ]),
            None
        );
        // A number past what the tool holds is the nearest it holds.
        assert_eq!(
            span(&[("cDNA_position", "5-99999999999999999999")]),
            Some((5, i64::MAX))
        );
    }

    #[test]
    fn a_digit_the_interpreter_does_not_convert_ends_the_run() {
        let e = row("CAT", &[("cDNA_position", "6\u{b2}-65")]);
        assert_eq!(
            uncaught(|| cdna_span(&e)),
            "invalid literal for int() with base 10: '6\u{b2}'"
        );
    }

    #[test]
    fn an_annotation_is_its_groups_whatever_their_order() {
        assert_eq!(
            annotation_groups("b=2:a=1,c=3"),
            annotation_groups("c=3,a=1:b=2")
        );
        assert_eq!(
            annotation_groups("b=2:a=1,c=3"),
            [vec!["a=1", "b=2"], vec!["c=3"]]
        );
        assert_ne!(annotation_groups("a=1,a=1"), annotation_groups("a=1"));
        assert_eq!(annotation_groups(""), [[""]]);
    }

    #[test]
    fn two_terms_of_the_five_that_differ_are_a_precedence_row() {
        let with = |name: &str, term: &str| row("T", &[(name, term)]);
        let lost = with("5UTR_consequence", UTRANNOTATOR_LOST_START_TERM);
        let stop = with("5UTR_consequence", UTRANNOTATOR_LOST_STOP_TERM);
        assert!(utrannotator_precedence_row(&lost, &stop));
        assert!(!utrannotator_precedence_row(&lost, &lost));
        assert!(!utrannotator_precedence_row(
            &lost,
            &with("5UTR_consequence", "-")
        ));
        assert!(!utrannotator_precedence_row(&lost, &row("T", &[])));
        // The JSON format's spelling is read, another case of the name is not.
        assert!(utrannotator_precedence_row(
            &lost,
            &with("5utr_consequence", UTRANNOTATOR_LOST_STOP_TERM)
        ));
        assert!(!utrannotator_precedence_row(
            &lost,
            &with("5UTR_Consequence", UTRANNOTATOR_LOST_STOP_TERM)
        ));
    }

    #[test]
    fn the_groups_of_an_annotation_are_read_in_the_three_formats_spellings() {
        let two = vec![
            frame(&[("ref_type", "uORF"), ("ref_type_length", "36")]),
            frame(&[("ref_type", "oORF")]),
        ];
        assert_eq!(
            annotation_frames("ref_type=uORF:ref_type_length=36,ref_type=oORF"),
            Some(two.clone())
        );
        assert_eq!(
            annotation_frames("ref_type=uORF:ref_type_length=36&ref_type=oORF"),
            Some(two.clone())
        );
        assert_eq!(
            annotation_frames(
                r#"{"1": {"ref_type": "uORF", "ref_type_length": 36}, "2": {"ref_type": "oORF"}}"#
            ),
            Some(two)
        );
        // A value of the JSON format's object is written as the interpreter writes it.
        assert_eq!(
            annotation_frames(
                r#"{"1": {"a": null, "b": true, "c": 1.0, "d": [1, "x"], "e": {"k": 1e400}}}"#
            ),
            Some(vec![frame(&[
                ("a", "None"),
                ("b", "True"),
                ("c", "1.0"),
                ("d", "[1, 'x']"),
                ("e", "{'k': inf}")
            ])])
        );
        // A key written twice keeps its first place and its last value.
        assert_eq!(
            annotation_frames("a=1:b=2:a=3"),
            Some(vec![frame(&[("a", "3"), ("b", "2")])])
        );
        assert_eq!(
            annotation_frames("a=1=2"),
            Some(vec![frame(&[("a", "1=2")])])
        );
        for unread in [
            "",
            "-",
            ABSENT,
            "a=1,b=2&c=3",
            "a=1:b",
            "{",
            r#"{"1": 1}"#,
            "[1]",
            r#"{"1": {}} x"#,
        ] {
            assert_eq!(annotation_frames(unread), None, "{unread}");
        }
        assert_eq!(annotation_frames("{}"), Some(Vec::new()));
    }

    #[test]
    fn an_annotation_nested_past_what_the_interpreter_writes_ends_the_run() {
        let nested = |depth: usize| {
            format!(
                r#"{{"1": {{"k": {}{}}}}}"#,
                "[".repeat(depth),
                "]".repeat(depth)
            )
        };
        let read = annotation_frames(&nested(WRITTEN_MAX_DEPTH)).expect("the object is read");
        assert_eq!(read[0]["k"].len(), 2 * WRITTEN_MAX_DEPTH);
        assert!(
            uncaught(|| annotation_frames(&nested(WRITTEN_MAX_DEPTH + 1)))
                .contains("while getting the repr of an object")
        );
        assert!(uncaught(|| annotation_frames(&nested(pyjson::MAX_DEPTH)))
            .contains("while decoding a JSON array"));
        assert!(uncaught(|| annotation_frames(r#"{"1": {"k": "\ud800"}}"#))
            .contains("half a surrogate pair"));
    }

    #[test]
    fn a_sweep_file_that_names_no_row_columns_is_refused() {
        let directory = std::env::temp_dir().join(format!(
            "vep-formats-compare-utrannotator-sweep-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).expect("a directory for the test's files");
        let planted = |name: &str, text: &str| {
            let path = directory.join(name);
            std::fs::write(&path, text).expect("the file is written");
            path
        };
        let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tExtra\n";
        let no_rows = planted("no_rows.txt", &format!("## a comment\n{header}"));
        let sweep = UtrannotatorSeedSweep::from_outputs(&[&no_rows, &no_rows])
            .expect("an output of no rows is read");
        assert!(sweep.written.is_empty());
        for (name, text) in [
            ("empty.txt", ""),
            ("comments.txt", "## a comment\n"),
            ("vcf.txt", "#CHROM\tPOS\tID\tREF\tALT\n"),
            ("two_of_three.txt", "#Location\tAllele\n"),
            // The first line that opens with one `#` is the one read.
            ("second.txt", &format!("#Gene\n{header}")),
        ] {
            let path = planted(name, text);
            let refused = UtrannotatorSeedSweep::from_outputs(&[&no_rows, &path])
                .expect_err("the file is refused");
            assert_eq!(
                refused.to_string(),
                format!(
                    "{}: not a default-format output (columns Location, Allele, Feature)",
                    path.display()
                )
            );
        }
        let absent = directory.join("absent.txt");
        let refused =
            UtrannotatorSeedSweep::from_outputs(&[&absent]).expect_err("the file is refused");
        assert_eq!(
            refused.to_string(),
            format!(
                "cannot read {}: No such file or directory",
                absent.display()
            )
        );
        std::fs::remove_dir_all(&directory).expect("the test's files are removed");
    }

    #[test]
    fn a_row_of_the_sweep_is_its_record_its_allele_and_its_feature() {
        let e = row("CCA", &[]);
        assert_eq!(
            UtrannotatorSeedSweep::row(&e),
            (
                "21:1000-1001 ins".to_string(),
                "CCA".to_string(),
                "ENST1".to_string()
            )
        );
    }

    #[test]
    fn a_start_codons_position_is_the_one_its_annotation_names() {
        let reports = |term: &str, annotation: &str| {
            row(
                "CAT",
                &[("5UTR_consequence", term), ("5UTR_annotation", annotation)],
            )
        };
        let (gain, loss) = (UTRANNOTATOR_GAINED_START_TERM, UTRANNOTATOR_LOST_START_TERM);
        let e = reports(gain, "CapDistanceToStart=64:type=uORF");
        assert_eq!(gained_start_position(&e), Some(64));
        assert_eq!(lost_start_position(&e), None);
        let e = reports(loss, r#"{"1":{"CapDistanceToStart":63,"type":"uORF"}}"#);
        assert_eq!(
            (gained_start_position(&e), lost_start_position(&e)),
            (None, Some(63))
        );
        // One position named twice is one position; two spellings of one number are two.
        for (annotation, at) in [
            ("CapDistanceToStart=-1:type=uORF", Some(-1)),
            ("CapDistanceToStart=7,CapDistanceToStart=7", Some(7)),
            ("CapDistanceToStart=7,CapDistanceToStart=07", None),
            ("CapDistanceToStart=7,CapDistanceToStart=9", None),
            ("type=uORF", None),
            ("-", None),
        ] {
            assert_eq!(
                gained_start_position(&reports(gain, annotation)),
                at,
                "{annotation}"
            );
        }
        assert_eq!(gained_start_position(&row("CAT", &[])), None);
    }

    #[test]
    fn an_allele_is_read_on_the_strand_its_row_states() {
        let on = |allele: &str, strand: &[(&str, &str)]| {
            allele_on_transcript_strand(&row(allele, strand))
        };
        assert_eq!(on("CCA", &[("STRAND", "1")]), Some("CCA".to_string()));
        assert_eq!(on("CCA", &[("STRAND", "-1")]), Some("TGG".to_string()));
        assert_eq!(on("CCA", &[("strand", "-1")]), Some("TGG".to_string()));
        assert_eq!(on("CCA", &[]), None);
        assert_eq!(on("CCA", &[("STRAND", "-")]), None);
        assert_eq!(on("-", &[("STRAND", "1")]), None);
        assert_eq!(on("<INS>", &[("STRAND", "1")]), None);
    }

    #[test]
    fn the_base_after_an_insertion_is_read_from_the_notation_and_its_offset() {
        let after = |fields: &[(&str, &str)]| utr_base_after_insertion(&row("A", fields));
        let cases: [(&str, Option<&str>, Option<i64>); 12] = [
            ("c.-55_-54insA", None, Some(54)),
            ("c.-53_-52insA", Some("2"), Some(54)),
            ("c.-53_-52insA", Some("-2"), Some(54)),
            ("c.-55dup", None, Some(54)),
            ("c.-56_-55dup", None, Some(54)),
            ("c.-54dup", Some("1"), Some(54)),
            ("c.-54_-53insA", Some("-"), Some(53)),
            ("c.-55_-53insA", None, None),
            ("c.-55_-54delinsAA", None, None),
            ("c.1_2insA", None, None),
            ("c.-55del", None, None),
            ("c.-55_-54insA ", None, None),
        ];
        for (notation, offset, base) in cases {
            let hgvsc = format!("ENST1.1:{notation}");
            let mut fields = vec![("HGVSc", hgvsc.as_str())];
            fields.extend(offset.map(|offset| ("HGVS_OFFSET", offset)));
            assert_eq!(after(&fields), base, "{notation} {offset:?}");
        }
        assert_eq!(after(&[]), None);
        assert_eq!(after(&[("hgvsc", "c.-55dup")]), None);
        assert_eq!(
            after(&[("hgvsc", "ENST1.1:c.-55dup"), ("hgvs_offset", "3")]),
            Some(57)
        );
        // The notation itself names the base after the last writing, the offset left out.
        let both = |fields: &[(&str, &str)]| utr_bases_after_insertion(&row("A", fields));
        assert_eq!(both(&[("HGVSc", "ENST1.1:c.-55_-54insA")]), Some((54, 54)));
        for offset in ["2", "-2"] {
            assert_eq!(
                both(&[("HGVSc", "ENST1.1:c.-53_-52insA"), ("HGVS_OFFSET", offset)]),
                Some((54, 52)),
                "{offset}"
            );
        }
        assert_eq!(
            both(&[("HGVSc", "ENST1.1:c.-54dup"), ("HGVS_OFFSET", "1")]),
            Some((54, 53))
        );
        assert_eq!(both(&[("HGVSc", "ENST1.1:c.-55del")]), None);
    }

    #[test]
    fn an_allele_spells_a_start_codon_from_a_base_before_its_end() {
        let spells = |allele: &str, strand: &str, before: i128| {
            start_codon_before_insertion(&row(allele, &[("STRAND", strand)]), before)
        };
        // `GTAATG` holds the codon from the base three before its end, and nine, and fifteen.
        for before in 0..=16 {
            assert_eq!(
                spells("GTAATG", "1", before),
                [3, 9, 15].contains(&before),
                "{before}"
            );
        }
        // Three bases that spell the codon from their last base: one base more than a
        // multiple of their length before the end, never fewer than three.
        for before in 0..=10 {
            assert_eq!(
                spells("TGA", "1", before),
                [4, 7, 10].contains(&before),
                "{before}"
            );
            assert_eq!(
                spells("TCA", "-1", before),
                [4, 7, 10].contains(&before),
                "{before}"
            );
            assert!(!spells("TCA", "1", before), "{before}");
        }
        // The codon itself stands a whole length before its own end.
        assert!(spells("ATG", "1", 3) && spells("CAT", "-1", 6));
        for (allele, before) in [("A", 3), ("AT", 4), ("ATGA", 3), ("-", 3), ("<INS>", 3)] {
            assert!(!spells(allele, "1", before), "{allele}");
        }
        assert!(!start_codon_before_insertion(&row("ATG", &[]), 3));
    }

    #[test]
    fn a_structural_rows_alternate_is_its_place_among_its_records_rows() {
        let structural = |word: &str, feature: &str| {
            let mut e = row(word, &[]);
            e.key = format!("21:1000-1001|{word}|{feature}|Transcript");
            e.feature = feature.to_string();
            e
        };
        let record = |alternates: &str| Entry {
            key: "21:1000-1001|record".to_string(),
            location: Some("21:1000-1001".to_string()),
            fields: [
                ("Uploaded_variation", "ins"),
                ("REF", "N"),
                ("ALT", alternates),
            ]
            .into_iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect(),
            ..Entry::default()
        };
        let word = STRUCTURAL_COPY_NUMBER_WORD;
        let (first, second) = (structural(word, "ENST1"), structural(word, "ENST1"));
        let elsewhere = structural(word, "ENST2");
        let stated = record("<CN0>,<CN2>");
        let rows = [&first, &elsewhere, &second, &stated];
        let texts = vec!["<CN0>".to_string(), "<CN2>".to_string()];
        assert_eq!(
            row_alternate(&first, Some(&rows)),
            Some((texts.clone(), Some(0)))
        );
        assert_eq!(
            row_alternate(&second, Some(&rows)),
            Some((texts.clone(), Some(1)))
        );
        // A row of another feature has no place among two rows of its own, and an equal row
        // that is none of the record's has none either.
        assert_eq!(
            row_alternate(&elsewhere, Some(&rows)),
            Some((texts.clone(), None))
        );
        assert_eq!(
            row_alternate(&structural(word, "ENST1"), Some(&rows)),
            Some((texts, None))
        );
        // No record handed, and rows among which none states the alternates.
        assert_eq!(row_alternate(&first, None), None);
        assert_eq!(row_alternate(&first, Some(&[&first, &second])), None);
    }

    #[test]
    fn the_seeded_runs_hold_each_rows_terms_and_what_is_written_beside_them() {
        let reports = |term: &str, annotation: &str| {
            row(
                "TG",
                &[("5UTR_consequence", term), ("5UTR_annotation", annotation)],
            )
        };
        let (gain, shift) = (UTRANNOTATOR_GAINED_START_TERM, UTRANNOTATOR_FRAMESHIFT_TERM);
        let gained = reports(gain, "CapDistanceToStart=145:type=uORF");
        let shifted = reports(shift, "alt_type=uORF:ref_type=uORF");
        let mut sweep = UtrannotatorSeedSweep::new();
        sweep.add_output(&[gained.clone(), reports("-", "-")]);
        sweep.add_output(std::slice::from_ref(&shifted));
        sweep.add_output(&[reports(gain, "type=uORF:CapDistanceToStart=145")]);
        assert_eq!(sweep.written.len(), 1);
        assert_eq!(sweep.gained_start(&shifted), Some(145));
        assert!(sweep.writes(&gained, &shifted, true));
        assert!(sweep.writes(
            &gained,
            &reports(shift, "ref_type=uORF:alt_type=uORF"),
            true
        ));
        assert!(!sweep.writes(&gained, &reports(shift, "ref_type=uORF"), true));
        assert!(sweep.writes(&gained, &reports(shift, "ref_type=uORF"), false));
        assert!(!sweep.writes(&gained, &reports(UTRANNOTATOR_LOST_STOP_TERM, "-"), false));
        assert!(sweep.writes_term(&gained, gain) && sweep.writes_term(&gained, shift));
        assert!(!sweep.writes_term(&gained, UTRANNOTATOR_LOST_STOP_TERM));
        // A row the runs do not hold, and an entry that names no Location.
        let mut elsewhere = gained.clone();
        elsewhere.feature = "ENST2".to_string();
        assert_eq!(sweep.gained_start(&elsewhere), None);
        assert!(!sweep.writes(&elsewhere, &shifted, false));
        assert!(!sweep.writes_term(&elsewhere, gain));
        let mut unplaced = gained.clone();
        unplaced.location = None;
        assert_eq!(sweep.gained_start(&unplaced), None);
        assert!(!sweep.writes(&unplaced, &shifted, false));
        assert!(!sweep.writes_term(&unplaced, gain));
        // Two positions beside the gained start name no one codon.
        sweep.add_output(&[reports(gain, "CapDistanceToStart=143:type=uORF")]);
        assert_eq!(sweep.gained_start(&shifted), None);
    }
}
