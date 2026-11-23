// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The UTRAnnotator plugin's classes and the seeded runs they read.
//!
//! The plugin reports one effect per row: `5UTR_consequence` names it with one of five
//! terms and `5UTR_annotation` lists the upstream ORFs it concerns as groups of
//! `key=value` pairs, the pairs joined by `:` and the groups by `,` (`&` in a VCF value).
//! Each class is one way the two engines differ in those two fields; a rule here answers
//! whether a row's two entries differ in that way.

use super::entry::{blank, entry_kind, field_named, field_of, parse_default, Entry, Row, ABSENT};
use super::record::record_label;
use super::shared::{
    allele_on_transcript_strand, alternates_in_order, cdna_span, digits_value, is_integer,
    is_sequence, list_parts, location_ends, re_digit, read_text, record_allele_lengths,
    HGVSC_FIELDS, HGVS_OFFSET_FIELDS,
};
use super::{Consult, Refusal};
use crate::pytext::{self, Json};
use flate2::read::MultiGzDecoder;
use regex::Regex;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::LazyLock;

pub const GAINED_START_TERM: &str = "5_prime_UTR_premature_start_codon_gain_variant";
pub const GAINED_STOP_TERM: &str = "5_prime_UTR_uORF_stop_codon_gain_variant";
pub const LOST_STOP_TERM: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";
pub const LOST_START_TERM: &str = "5_prime_UTR_premature_start_codon_loss_variant";
pub const FRAMESHIFT_TERM: &str = "5_prime_UTR_uORF_frameshift_variant";
/// The terms of the plugin's five analyses (UTRAnnotator.pm:426, 562, 704, 838, 1003).
pub const TERMS: [&str; 5] = [
    GAINED_START_TERM,
    GAINED_STOP_TERM,
    LOST_STOP_TERM,
    LOST_START_TERM,
    FRAMESHIFT_TERM,
];
pub const CONSEQUENCE_FIELDS: [&str; 2] = ["5UTR_consequence", "5utr_consequence"];
pub const ANNOTATION_FIELDS: [&str; 2] = ["5UTR_annotation", "5utr_annotation"];
/// The count of the reference 5' UTR's upstream ORFs, which the plugin writes on every row
/// it annotates.
const COUNT_FIELDS: [&str; 1] = ["Existing_uORFs"];
const START_CODON: &str = "ATG";
const UTR_TERM: &str = "5_prime_UTR_variant";

/// The Allele of a structural allele's row: the Sequence Ontology term of the record's
/// class (ensembl-variation Utils/Config.pm:1417-1438). The breakpoint class's term is left
/// out: its row does not tell whether the plugin read the REF column or the ALT text.
const STRUCTURAL_CLASS_WORDS: [&str; 17] = [
    "insertion",
    "mobile_element_insertion",
    "Alu_insertion",
    "HERV_insertion",
    "LINE1_insertion",
    "SVA_insertion",
    "deletion",
    "mobile_element_deletion",
    "Alu_deletion",
    "HERV_deletion",
    "LINE1_deletion",
    "SVA_deletion",
    "tandem_repeat",
    "tandem_duplication",
    "duplication",
    "copy_number_variation",
    "inversion",
];
/// A plain insertion's class word; a mobile element's insertion has a word of its own.
const STRUCTURAL_INSERTION_WORD: &str = "insertion";
/// The class word of a record whose alternates state a loss beside a gain.
const STRUCTURAL_COPY_NUMBER_WORD: &str = "copy_number_variation";
const STRUCTURAL_INVERSION_WORD: &str = "inversion";
/// The Consequence terms the row of a deletion inside one exon of the 5' UTR may name
/// beside `5_prime_UTR_variant`.
const STRUCTURAL_DELETION_ROW_TERMS: [&str; 2] = ["feature_truncation", "NMD_transcript_variant"];
/// The Consequence term the row of an inversion or of a copy-number record inside one exon
/// of the 5' UTR may name beside `5_prime_UTR_variant`.
const STRUCTURAL_SPAN_ROW_TERMS: [&str; 1] = ["NMD_transcript_variant"];
/// The ALT texts of a symbolic deletion.
const DELETION_ALTERNATE_TEXTS: [&str; 3] = ["<DEL>", "<CN0>", "<CN=0>"];
const INVERSION_ALTERNATE_TEXT: &str = "<INV>";
/// The keys of an effect that are read from the reference UTR alone.
const REFERENCE_UTR_KEYS: [&str; 6] = [
    "Evidence",
    "KozakContext",
    "KozakStrength",
    "ref_StartDistanceToCDS",
    "ref_type",
    "ref_type_length",
];
const ALTERNATE_TYPE: (&str, &str) = ("alt_type", "uORF");

/// The position of a reported start codon inside `5UTR_annotation`:
/// `CapDistanceToStart=64` in the text formats, a member of the JSON object.
const CAP_DISTANCE: &str = r#"CapDistanceToStart"?[{space}]*[=:][{space}]*"?(-?[{digit}]+)"#;

fn cap_distance(anchored: bool) -> Regex {
    let pattern = CAP_DISTANCE
        .replace("{space}", pytext::RE_SPACE)
        .replace("{digit}", re_digit());
    let pattern = if anchored {
        format!("^(?:{pattern})$")
    } else {
        pattern
    };
    Regex::new(&pattern).expect("regex")
}

static CAP_DISTANCE_RE: LazyLock<Regex> = LazyLock::new(|| cap_distance(false));
/// `CAP_DISTANCE_RE` matched against a whole pair.
static CAP_DISTANCE_PAIR_RE: LazyLock<Regex> = LazyLock::new(|| cap_distance(true));

/// An insertion between two 5' UTR positions as an HGVSc body writes it, matched whole.
static UTR_HGVSC_INSERTION_RE: LazyLock<Regex> = LazyLock::new(|| {
    let digit = re_digit();
    Regex::new(&format!(r"^c\.(-[{digit}]+)_(-[{digit}]+)ins[A-Za-z]+$")).expect("regex")
});
/// A duplication of one 5' UTR position or of a span of them, matched whole.
static UTR_HGVSC_DUPLICATION_RE: LazyLock<Regex> = LazyLock::new(|| {
    let digit = re_digit();
    Regex::new(&format!(r"^c\.(-[{digit}]+)(?:_(-[{digit}]+))?dup$")).expect("regex")
});

/// A text-format `5UTR_annotation` with the order of its groups and of each group's pairs
/// set aside (UTRAnnotator.pm:282-283 walk both as hashes): the pairs of each group sorted
/// and joined by `:`, the groups sorted and joined by `,`. Two values give one text exactly
/// where they carry the same groups of the same pairs, each as often.
pub fn annotation_groups(value: &str) -> String {
    let mut groups: Vec<String> = value
        .split(',')
        .map(|group| {
            let mut pairs: Vec<&str> = group.split(':').collect();
            pairs.sort_unstable();
            pairs.join(":")
        })
        .collect();
    groups.sort_unstable();
    groups.join(",")
}

/// The annotations written beside one term on one row, each as `annotation_groups` gives it.
#[derive(Debug)]
struct Written {
    /// The term's index in `TERMS`.
    term: u8,
    annotations: Vec<Box<str>>,
}

/// What Ensembl VEP's UTRAnnotator plugin writes on each row of an input across runs under
/// different Perl hash seeds (`PERL_HASH_SEED` with `PERL_PERTURB_KEYS=0`), each of the
/// five flags of UTRAnnotator.pm:262-267 walked last under one of them: per row of the
/// default-format outputs (the Location with the `Uploaded_variation`, the Allele and the
/// Feature), every `5UTR_consequence` term with the `5UTR_annotation` values written
/// beside it. The plugin keeps the last flag that fired, so the terms of a row are the
/// analyses that fire on it.
#[derive(Debug, Default)]
pub struct UtrannotatorSeedSweep {
    /// Per row (`row_key`), the terms written on it.
    written: HashMap<Box<str>, Vec<Written>>,
}

/// The row an entry stands on: its record's label, its Allele and its Feature. No column
/// of an output holds a tab, so the joined text names one row.
fn row_key(e: &Entry) -> String {
    format!("{}\t{}\t{}", record_label(e), e.allele, e.feature)
}

/// The lines of a text whose line ends are `\n`, `\r\n` or `\r`, without their ends; a
/// last line without an end is a line where it holds anything.
struct Lines<R: BufRead> {
    reader: R,
    line: Vec<u8>,
    /// Whether the last line ended in `\r`, so a `\n` that follows belongs to its end.
    after_return: bool,
}

impl<R: BufRead> Lines<R> {
    fn new(reader: R) -> Lines<R> {
        Lines {
            reader,
            line: Vec::new(),
            after_return: false,
        }
    }

    fn next(&mut self) -> io::Result<Option<&[u8]>> {
        self.line.clear();
        loop {
            let buffer = self.reader.fill_buf()?;
            if buffer.is_empty() {
                return Ok((!self.line.is_empty()).then_some(self.line.as_slice()));
            }
            let skipped = usize::from(self.after_return && buffer[0] == b'\n');
            self.after_return = false;
            let rest = &buffer[skipped..];
            match memchr::memchr2(b'\n', b'\r', rest) {
                Some(at) => {
                    self.line.extend_from_slice(&rest[..at]);
                    self.after_return = rest[at] == b'\r';
                    self.reader.consume(skipped + at + 1);
                    return Ok(Some(self.line.as_slice()));
                }
                None => {
                    self.line.extend_from_slice(rest);
                    let read = buffer.len();
                    self.reader.consume(read);
                }
            }
        }
    }
}

impl UtrannotatorSeedSweep {
    /// The rows of the default-format outputs `paths`, one run each, a file whose name
    /// ends in `.gz` gunzipped. Each output is read once as a stream: no more of it is
    /// held than one row and what its rows with a term state. A file that cannot be read
    /// is refused, then a file whose first header line does not name the columns a row is
    /// keyed by (read as an output, it would write no term on any row), then a row
    /// `parse_default` refuses.
    pub fn from_outputs(paths: &[PathBuf]) -> Result<UtrannotatorSeedSweep, Refusal> {
        let mut sweep = UtrannotatorSeedSweep::default();
        for path in paths {
            match sweep.add_file(path) {
                Ok(read) => read?,
                // `read_text` words the refusal of a file that cannot be read.
                Err(_) => {
                    let text = read_text(path)?;
                    sweep
                        .add_lines(Lines::new(text.as_bytes()), path)
                        .map_err(|error| {
                            Refusal(format!("cannot read {}: {error}", path.display()))
                        })??;
                }
            }
        }
        sweep.written.shrink_to_fit();
        Ok(sweep)
    }

    /// Adds the rows of the output at `path`; the outer error is the stream's.
    fn add_file(&mut self, path: &Path) -> io::Result<Result<(), Refusal>> {
        let file = File::open(path)?;
        if path.as_os_str().as_encoded_bytes().ends_with(b".gz") {
            let lines = Lines::new(BufReader::with_capacity(1 << 16, MultiGzDecoder::new(file)));
            self.add_lines(lines, path)
        } else {
            self.add_lines(Lines::new(BufReader::with_capacity(1 << 16, file)), path)
        }
    }

    /// Adds the rows of one run's output, each read under the header line before it. The
    /// first refusal of a row stands until the stream's end: a stream that cannot be read
    /// to its end, or is no UTF-8, is the outer error whatever its rows hold, and a first
    /// header line without the key columns is refused before any row.
    fn add_lines<R: BufRead>(
        &mut self,
        mut lines: Lines<R>,
        path: &Path,
    ) -> io::Result<Result<(), Refusal>> {
        let mut header: Option<String> = None;
        let mut keyed = false;
        let mut refused: Option<Refusal> = None;
        let mut text = String::new();
        while let Some(line) = lines.next()? {
            let line = std::str::from_utf8(line)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
            if let Some(names) = line.strip_prefix('#') {
                if !names.starts_with('#') {
                    if header.is_none() {
                        let names: Vec<&str> = names.split('\t').collect();
                        keyed = ["Location", "Allele", "Feature"]
                            .iter()
                            .all(|column| names.contains(column));
                    }
                    header = Some(format!("{line}\n"));
                }
                continue;
            }
            if refused.is_some() {
                continue;
            }
            text.clear();
            text.push_str(header.as_deref().unwrap_or(""));
            text.push_str(line);
            text.push('\n');
            match parse_default(&text) {
                Ok(rows) => rows.iter().for_each(|row| self.add(row)),
                Err(refusal) => refused = Some(refusal),
            }
        }
        if !keyed {
            return Ok(Err(Refusal(format!(
                "{}: not a default-format output (columns Location, Allele, Feature)",
                path.display()
            ))));
        }
        Ok(refused.map_or(Ok(()), Err))
    }

    /// One row of a run: its term, where it is one of the plugin's five, with the
    /// annotation beside it.
    fn add(&mut self, e: &Entry) {
        let term = field_of(e, &CONSEQUENCE_FIELDS);
        let Some(term) = TERMS.iter().position(|known| *known == term) else {
            return;
        };
        let annotation = annotation_groups(field_of(e, &ANNOTATION_FIELDS));
        let terms = self.written.entry(row_key(e).into_boxed_str()).or_default();
        match terms.iter_mut().find(|held| usize::from(held.term) == term) {
            Some(held) => {
                if !held.annotations.iter().any(|known| **known == *annotation) {
                    held.annotations.push(annotation.into_boxed_str());
                }
            }
            None => terms.push(Written {
                term: term as u8,
                annotations: vec![annotation.into_boxed_str()],
            }),
        }
    }

    /// The annotations the runs write beside `term` on `e`'s row. An entry without a
    /// Location matches no row.
    fn annotations(&self, e: &Entry, term: &str) -> Option<&[Box<str>]> {
        e.location.as_ref()?;
        let term = TERMS.iter().position(|known| *known == term)?;
        self.written
            .get(row_key(e).as_str())?
            .iter()
            .find(|held| usize::from(held.term) == term)
            .map(|held| held.annotations.as_slice())
    }

    /// The one `CapDistanceToStart` the runs write beside the gained-start term on `e`'s
    /// row; `None` where no run writes that term there or they write more than one
    /// position.
    pub fn gained_start(&self, e: &Entry) -> Option<i128> {
        let mut positions = self
            .annotations(e, GAINED_START_TERM)?
            .iter()
            .flat_map(|groups| groups.split(','))
            .flat_map(|group| group.split(':'))
            .filter(|pair| CAP_DISTANCE_PAIR_RE.is_match(pair));
        let pair = positions.next()?;
        if positions.any(|other| other != pair) {
            return None;
        }
        pytext::int(CAP_DISTANCE_PAIR_RE.captures(pair)?.get(1)?.as_str())
    }

    /// Whether a run writes `term` on `e`'s row.
    pub fn writes_term(&self, e: &Entry, term: &str) -> bool {
        self.annotations(e, term).is_some()
    }

    /// Whether the runs write, on `e`'s row, the term of the reference's entry `e` and the
    /// term of the engine's entry `a`, and under `with_annotation` beside each the entry's
    /// own `5UTR_annotation`.
    pub fn writes(&self, e: &Entry, a: &Entry, with_annotation: bool) -> bool {
        [e, a].iter().all(|entry| {
            let term = field_of(entry, &CONSEQUENCE_FIELDS);
            self.annotations(e, term).is_some_and(|annotations| {
                !with_annotation || {
                    let own = annotation_groups(field_of(entry, &ANNOTATION_FIELDS));
                    annotations.iter().any(|known| **known == *own)
                }
            })
        })
    }
}

/// The `CapDistanceToStart` of the start codon the entry reports under `term`: its
/// `5UTR_consequence` is the term and its `5UTR_annotation` names one position. `None`
/// otherwise, and for a position outside `i128`.
fn reported_position(e: &Entry, term: &str) -> Option<i128> {
    if field_of(e, &CONSEQUENCE_FIELDS) != term {
        return None;
    }
    let mut positions = CAP_DISTANCE_RE
        .captures_iter(field_of(e, &ANNOTATION_FIELDS))
        .filter_map(|found| found.get(1))
        .map(|position| position.as_str());
    let position = positions.next()?;
    if positions.any(|other| other != position) {
        return None;
    }
    pytext::int(position)
}

/// The `CapDistanceToStart` of the start codon the entry reports as gained.
pub fn gained_start_position(e: &Entry) -> Option<i128> {
    reported_position(e, GAINED_START_TERM)
}

/// The `CapDistanceToStart` of the start codon the entry reports as lost.
pub fn lost_start_position(e: &Entry) -> Option<i128> {
    reported_position(e, LOST_START_TERM)
}

/// Whether the two entries of a row differ as the two gained-start windows make them
/// (`utrannotator_gained_start_window`). Both windows open two bases before the allele's
/// first base in the mutated UTR and close two bases after its last: the reference's
/// (UTRAnnotator.pm:345) counts the reference allele's bases, the engine's the alternate
/// allele's, and each reports the first start codon of its window by the 0-based index of
/// its first base. So the row is the class's when one engine alone reports the gained
/// start, its codon begins inside that engine's allele length and at or past the other's,
/// and the reporting engine is the one whose allele is the longer; or when the reference
/// reports it at a negative index, which only a window opened before the UTR's first base
/// yields. An insertion's row is refused where the reference reports a start codon lost and
/// the engine's gained start stands the allele's length past that codon's index: that
/// codon is the reference's own, moved on by the inserted bases. On an insertion's row
/// where the engine alone reports the gained start, its codon may also begin three or more
/// bases before the base after the insertion: the engine opens its window at the first
/// index the mutated sequence can be written at, and the reference's four-base window
/// (UTRAnnotator.pm:339-345) at the record's own index. The row is the class's where the
/// same sequence written that many bases earlier shows a start codon there
/// (`start_codon_written_earlier`).
pub fn gained_start_window_row(e: &Entry, a: &Entry, record_rows: Option<&[Row]>) -> bool {
    let reference_gains = field_of(e, &CONSEQUENCE_FIELDS) == GAINED_START_TERM;
    let candidate_gains = field_of(a, &CONSEQUENCE_FIELDS) == GAINED_START_TERM;
    let Some(at) = gained_start_position(if reference_gains { e } else { a }) else {
        return false;
    };
    if reference_gains && at < 0 {
        return true;
    }
    let (Some((reference_length, allele_length)), Some((first, last))) =
        (record_allele_lengths(e, record_rows), cdna_span(e))
    else {
        return false;
    };
    if reference_gains == candidate_gains {
        return false;
    }
    let within = || -> Option<bool> {
        // A Location of one position more than the largest `i128` gives that largest number
        // as the reference length, and no cDNA width equals such a span.
        if reference_length == i128::MAX
            && e.location
                .as_deref()
                .is_some_and(|location| stated_span(location).is_none())
        {
            return Some(false);
        }
        // The cDNA positions are the reference allele's own, an insertion's its two
        // flanking bases.
        let written = if reference_length == 0 {
            2
        } else {
            reference_length
        };
        if last.checked_sub(first)?.checked_add(1)? != written {
            return Some(false);
        }
        // A codon the reference reports lost that stands the allele's length further on
        // was moved, not made, by the inserted bases.
        if candidate_gains && reference_length == 0 {
            if let Some(lost) = lost_start_position(e) {
                if at == lost.checked_add(allele_length)? {
                    return Some(false);
                }
            }
        }
        let start = if reference_length == 0 {
            first
        } else {
            first.checked_sub(1)?
        };
        if candidate_gains && reference_length == 0 && at < start {
            return Some(start_codon_written_earlier(e, start.checked_sub(at)?));
        }
        let (reporting, other) = if reference_gains {
            (reference_length, allele_length)
        } else {
            (allele_length, reference_length)
        };
        Some(start.checked_add(other)? <= at && at < start.checked_add(reporting)?)
    };
    within().unwrap_or(false)
}

/// Whether an insertion's allele, read on the transcript's strand and round its end,
/// spells a start codon from the base that stands `before` bases before its end: the three
/// bases the mutated sequence holds `before` bases before the insertion wherever it can be
/// written that many bases earlier, since the bases passed over repeat the allele. A codon
/// fewer than three bases before the insertion ends inside the allele as the record writes
/// it and is refused, as is a row that states no strand (`allele_on_transcript_strand`).
fn start_codon_written_earlier(e: &Entry, before: i128) -> bool {
    if before < 3 {
        return false;
    }
    let Some(read) = allele_on_transcript_strand(e) else {
        return false;
    };
    let bases = read.as_bytes();
    let length = bases.len() as i128;
    let offset = (length - before % length) % length;
    START_CODON.bytes().enumerate().all(|(step, base)| {
        let at = (offset + step as i128) % length;
        bases[at as usize] == base
    })
}

/// Whether the two entries' `5UTR_consequence` values are both among the plugin's five
/// terms and differ (`utrannotator_consequence_precedence`).
pub fn precedence_row(e: &Entry, a: &Entry) -> bool {
    let ev = e.fields.first_of(&CONSEQUENCE_FIELDS, ABSENT);
    let av = a.fields.first_of(&CONSEQUENCE_FIELDS, ABSENT);
    ev != av && TERMS.contains(&ev) && TERMS.contains(&av)
}

/// The number of positions a Location names after its first `:`; `None` where an end is no
/// number.
fn stated_span(location: &str) -> Option<i128> {
    let (start, end) = location_ends(location);
    let first = digits_value(start)?;
    let last = if end.is_empty() {
        first
    } else {
        digits_value(end)?
    };
    last.checked_sub(first)?.checked_add(1)
}

/// Whether an entry writes one of the plugin's five terms, or neither of the two fields.
fn term_or_none(e: &Entry) -> bool {
    let term = field_of(e, &CONSEQUENCE_FIELDS);
    if blank(term) {
        blank(field_of(e, &ANNOTATION_FIELDS))
    } else {
        TERMS.contains(&term)
    }
}

/// Whether the reference's entry `e` carries the effect of a structural allele's ALT text
/// read as deleted bases and the engine's entry `a` none
/// (`utrannotator_structural_allele_text`): the row's Allele is a structural class word on
/// both sides, the reference's `5UTR_consequence` is one of the plugin's five terms, and
/// the engine fills neither `5UTR_consequence` nor `5UTR_annotation`.
pub fn structural_allele_text_row(e: &Entry, a: &Entry) -> bool {
    STRUCTURAL_CLASS_WORDS.contains(&e.allele.as_str())
        && a.allele == e.allele
        && TERMS.contains(&field_of(e, &CONSEQUENCE_FIELDS))
        && blank(field_of(a, &CONSEQUENCE_FIELDS))
        && blank(field_of(a, &ANNOTATION_FIELDS))
}

/// The number of bases a structural record's row states inside one exon of the 5' UTR,
/// where the row's two entries fit an effect of that span. The row tells it by its own
/// columns: its Allele is the same on both sides; its Location's span is the distance
/// between its two cDNA positions, so both ends lie in one exon; its Consequence names
/// `5_prime_UTR_variant` beside no term but those of `beside`, so neither end leaves the
/// 5' UTR. The engine's `5UTR_consequence` is one of the plugin's five terms beside a
/// `5UTR_annotation`; the reference's `5UTR_consequence` is one of the five, or both of
/// its fields are empty. `None` otherwise, and for an entry without a Location.
fn span_inside_one_utr_exon(e: &Entry, a: &Entry, beside: &[&str]) -> Option<i128> {
    if a.allele != e.allele {
        return None;
    }
    let span = stated_span(e.location.as_deref()?)?;
    let (first, last) = cdna_span(e)?;
    if span < 1 || last.checked_sub(first)?.checked_add(1)? != span {
        return None;
    }
    let mut terms = e.consequence.split(',');
    if !terms.all(|term| term == UTR_TERM || beside.contains(&term))
        || !e.consequence.split(',').any(|term| term == UTR_TERM)
    {
        return None;
    }
    if !term_or_none(e) {
        return None;
    }
    if !TERMS.contains(&field_of(a, &CONSEQUENCE_FIELDS)) || blank(field_of(a, &ANNOTATION_FIELDS))
    {
        return None;
    }
    Some(span)
}

/// Whether `term` is one the deletion of `span` bases can fire: any of the plugin's five
/// but the gained stop, which is a substitution's alone (UTRAnnotator.pm:474), the
/// frameshift only where the span is no multiple of three (:923).
fn deletion_term(term: &str, span: i128) -> bool {
    term != GAINED_STOP_TERM && (span % 3 != 0 || term != FRAMESHIFT_TERM)
}

/// Whether the engine's entry `a` carries the effect of a symbolic deletion's stated span
/// on a row where that span lies inside one exon of the 5' UTR
/// (`utrannotator_structural_allele_text`). The row's Allele is a deletion's class word and
/// its columns state the span as `span_inside_one_utr_exon` reads it, its Consequence
/// naming `5_prime_UTR_variant` beside no term but `feature_truncation` and
/// `NMD_transcript_variant`. The engine's `5UTR_consequence` is a term a deletion of that
/// span can fire (`deletion_term`).
pub fn structural_deletion_span_row(e: &Entry, a: &Entry) -> bool {
    let deletion = e.allele == "deletion" || e.allele.ends_with("_deletion");
    deletion
        && STRUCTURAL_CLASS_WORDS.contains(&e.allele.as_str())
        && span_inside_one_utr_exon(e, a, &STRUCTURAL_DELETION_ROW_TERMS)
            .is_some_and(|span| deletion_term(field_of(a, &CONSEQUENCE_FIELDS), span))
}

/// The ALT texts of the record of the structural row `e` in the record's own order
/// (`alternates_in_order`), with the place among them of the alternate whose row `e` is. A
/// record writes one row for each of its alternates on a Feature, in the alternates'
/// order, so the place is `e`'s own among the record's rows on its Feature; no place where
/// those rows are not as many as the alternates or `e` is none of them. `None` where no
/// record is handed with the row or none states its alternates.
fn row_alternate(e: &Entry, record_rows: Option<&[Row]>) -> Option<(Vec<String>, Option<usize>)> {
    let rows = record_rows?;
    let alternates = alternates_in_order(e, rows)?;
    let (kind, label) = (entry_kind(e), record_label(e));
    let written: Vec<&Row> = rows
        .iter()
        .filter(|row| {
            row.feature == e.feature && entry_kind(row) == kind && record_label(row) == label
        })
        .collect();
    let place = if written.len() == alternates.len() {
        written
            .iter()
            .position(|row| std::ptr::eq(Rc::as_ptr(row), e))
    } else {
        None
    };
    Some((alternates, place))
}

/// Whether the engine's entry `a` carries the effect of the deletion of a copy-number
/// record's span, which the engine writes on the row of the record's loss alternate, where
/// that span lies inside one exon of the 5' UTR (`utrannotator_structural_allele_text`).
/// The reference reads the first ALT text as deleted bases on every row of the record
/// (UTRAnnotator.pm:151-153). The row's Allele is `copy_number_variation` and its columns
/// state the span as `span_inside_one_utr_exon` reads it, its Consequence naming
/// `5_prime_UTR_variant` beside no term but `NMD_transcript_variant`. The engine's
/// `5UTR_consequence` is a term a deletion of that span can fire (`deletion_term`). Where
/// the row's record is handed and states its alternates (`row_alternate`), they are two or
/// more and the row is that of an alternate whose own text is a deletion's.
pub fn structural_copy_number_loss_row(e: &Entry, a: &Entry, record_rows: Option<&[Row]>) -> bool {
    e.allele == STRUCTURAL_COPY_NUMBER_WORD
        && span_inside_one_utr_exon(e, a, &STRUCTURAL_SPAN_ROW_TERMS)
            .is_some_and(|span| deletion_term(field_of(a, &CONSEQUENCE_FIELDS), span))
        && row_alternate(e, record_rows).is_none_or(|(alternates, place)| {
            alternates.len() > 1
                && place
                    .is_some_and(|at| DELETION_ALTERNATE_TEXTS.contains(&alternates[at].as_str()))
        })
}

/// Whether the engine's entry `a` carries the effect of a symbolic inversion's span
/// replaced by its reverse complement, on a row where that span lies inside one exon of
/// the 5' UTR (`utrannotator_structural_allele_text`). The reference reads the ALT text as
/// deleted bases (UTRAnnotator.pm:151-153). The row's Allele is `inversion` and its
/// columns state the span as `span_inside_one_utr_exon` reads it, its Consequence naming
/// `5_prime_UTR_variant` beside no term but `NMD_transcript_variant`. The engine's
/// `5UTR_consequence` is a term a replacement of equal length can fire: any of the
/// plugin's five but the frameshift (UTRAnnotator.pm:889). Where the row's record is
/// handed and states its alternates (`row_alternate`), `<INV>` is its one alternate.
pub fn structural_inversion_span_row(e: &Entry, a: &Entry, record_rows: Option<&[Row]>) -> bool {
    e.allele == STRUCTURAL_INVERSION_WORD
        && span_inside_one_utr_exon(e, a, &STRUCTURAL_SPAN_ROW_TERMS).is_some()
        && field_of(a, &CONSEQUENCE_FIELDS) != FRAMESHIFT_TERM
        && row_alternate(e, record_rows)
            .is_none_or(|(alternates, _)| alternates == [INVERSION_ALTERNATE_TEXT])
}

/// Whether the engine's entry `a` carries the frameshift a plain symbolic insertion's
/// stated length decides (`utrannotator_structural_allele_text`). The row's Allele is
/// `insertion` on both sides, its Consequence names `5_prime_UTR_variant` and its
/// Location's span, the stated length, is no multiple of three (UTRAnnotator.pm:889,
/// :923). The engine's `5UTR_consequence` is the frameshift, and each effect of its
/// `5UTR_annotation` holds exactly the six keys read from the reference UTR, or those six
/// and `alt_type=uORF`; the reference's `5UTR_consequence` is one of the plugin's five
/// terms, or both of its fields are empty. An entry without a Location is not the class's.
pub fn structural_insertion_length_row(e: &Entry, a: &Entry) -> bool {
    if e.allele != STRUCTURAL_INSERTION_WORD || a.allele != e.allele {
        return false;
    }
    let Some(span) = e.location.as_deref().and_then(stated_span) else {
        return false;
    };
    if span < 1 || span % 3 == 0 || !e.consequence.split(',').any(|term| term == UTR_TERM) {
        return false;
    }
    if !term_or_none(e) {
        return false;
    }
    let annotation = field_of(a, &ANNOTATION_FIELDS);
    if field_of(a, &CONSEQUENCE_FIELDS) != FRAMESHIFT_TERM || blank(annotation) {
        return false;
    }
    annotation.split(',').all(|effect| {
        let pairs: Vec<(&str, &str)> = effect
            .split(':')
            .map(|pair| {
                let (key, _, value) = pytext::partition(pair, "=");
                (key, value)
            })
            .collect();
        let mut keys: Vec<&str> = pairs.iter().map(|(key, _)| *key).collect();
        keys.sort_unstable();
        keys.dedup();
        if keys.len() != pairs.len() {
            return false;
        }
        let typed = keys.contains(&ALTERNATE_TYPE.0);
        let read_from_reference = keys.len() - usize::from(typed) == REFERENCE_UTR_KEYS.len()
            && REFERENCE_UTR_KEYS.iter().all(|key| keys.contains(key));
        read_from_reference && (!typed || pairs.contains(&ALTERNATE_TYPE))
    })
}

/// Whether the record of an insertion's row states an empty reference allele where it
/// states its alleles: `Entry.insertion`, else `record_allele_lengths`; a row under an
/// identifier with no record handed is read by its positions and its allele alone.
fn reference_allele_is_empty_where_stated(e: &Entry, record_rows: Option<&[Row]>) -> bool {
    match e.insertion {
        Some(insertion) => insertion,
        None => record_allele_lengths(e, record_rows).is_none_or(|(reference, _)| reference == 0),
    }
}

/// Whether the row's two cDNA positions are an insertion's two flanking bases; the lower
/// of them, the 5' UTR index of the base after the insertion (UTRAnnotator.pm:1285-1313),
/// where they are.
fn insertion_index(e: &Entry) -> Option<i128> {
    let (first, last) = cdna_span(e)?;
    (last.checked_sub(first) == Some(1)).then_some(first)
}

/// Whether the reference's entry reports as lost a start codon that an insertion only
/// moves on and the engine's entry, which carries the plugin's counts, reports no lost
/// start (`utrannotator_lost_start_after_insertion`). A codon reported at the insertion's
/// index begins directly after the insertion (UTRAnnotator.pm:769) and stands whole
/// whatever the allele. A codon reported one index below it holds the insertion after its
/// first base and stands whole further on when the allele, read on the transcript's
/// strand, ends with that base; one reported two below holds it after its second base and
/// stands whole when the allele so read ends with the codon's first two bases
/// (`allele_on_transcript_strand`; a row that states no strand is refused there). Where
/// the row's record states its alleles the reference allele is empty.
pub fn lost_start_after_insertion_row(e: &Entry, a: &Entry, record_rows: Option<&[Row]>) -> bool {
    if field_of(a, &CONSEQUENCE_FIELDS) == LOST_START_TERM || blank(field_of(a, &COUNT_FIELDS)) {
        return false;
    }
    let (Some(lost), Some(index)) = (lost_start_position(e), insertion_index(e)) else {
        return false;
    };
    if !is_sequence(&e.allele) {
        return false;
    }
    // How many of the codon's bases stand before the insertion.
    match index.checked_sub(lost) {
        Some(0) => {}
        Some(before @ 1..=2) => {
            let moved_on = allele_on_transcript_strand(e)
                .is_some_and(|read| read.ends_with(&START_CODON[..before as usize]));
            if !moved_on {
                return false;
            }
        }
        _ => return false,
    }
    reference_allele_is_empty_where_stated(e, record_rows)
}

/// Whether the reference reports as gained, on an insertion's row, a start codon the
/// inserted bases leave at its own index, and the engine's entry, which carries the
/// plugin's counts, reports no gained start there
/// (`utrannotator_gained_start_remade_in_place`). A start codon reported one index below
/// the insertion's holds the insertion after its first base, and the mutated UTR reads a
/// start codon at that index when the allele, read on the transcript's strand
/// (`allele_on_transcript_strand`; a row that states no strand is refused), opens with the
/// codon's last two bases; one reported two below holds it after its second base and the
/// allele so read opens with the codon's last base. The index is the one the reference's
/// own `5UTR_annotation` names beside the gained-start term; where the reference writes
/// another of the plugin's terms, the one index the runs of `sweep` write beside the
/// gained start on the row. A gained start at a negative index, which the reference's
/// window yields by wrapping before the UTR's first base, is no codon of the UTR and is
/// refused. Where the row's record states its alleles the reference allele is empty.
pub fn gained_start_remade_in_place_row(
    e: &Entry,
    a: &Entry,
    record_rows: Option<&[Row]>,
    sweep: Option<&UtrannotatorSeedSweep>,
) -> bool {
    if blank(field_of(a, &COUNT_FIELDS)) {
        return false;
    }
    let (Some(index), Some(read)) = (insertion_index(e), allele_on_transcript_strand(e)) else {
        return false;
    };
    let term = field_of(e, &CONSEQUENCE_FIELDS);
    let at = if term == GAINED_START_TERM {
        gained_start_position(e)
    } else if TERMS.contains(&term) {
        sweep.and_then(|sweep| sweep.gained_start(e))
    } else {
        None
    };
    let Some(at) = at.filter(|at| *at >= 0) else {
        return false;
    };
    // How many of the codon's bases stand before the insertion.
    let remade = match index.checked_sub(at) {
        Some(before @ 1..=2) => read.starts_with(&START_CODON[before as usize..]),
        _ => false,
    };
    if !remade {
        return false;
    }
    if field_of(a, &CONSEQUENCE_FIELDS) == GAINED_START_TERM
        && gained_start_position(a).is_none_or(|own| own == at)
    {
        return false;
    }
    reference_allele_is_empty_where_stated(e, record_rows)
}

/// One group of a `5UTR_annotation`: its keys with their values, a key once, in the order
/// the group first writes each.
type Frame = Vec<(String, String)>;

fn value_of<'a>(frame: &'a [(String, String)], key: &str) -> Option<&'a str> {
    frame
        .iter()
        .find_map(|(held, value)| (held == key).then_some(value.as_str()))
}

/// `repr()` of a parsed JSON value.
fn python_repr(value: &Json) -> String {
    match value {
        Json::Str(text) => pytext::repr(text),
        other => python_str(other),
    }
}

/// `str()` of a parsed JSON value.
fn python_str(value: &Json) -> String {
    match value {
        Json::Null => "None".to_string(),
        Json::Bool(true) => "True".to_string(),
        Json::Bool(false) => "False".to_string(),
        Json::Int(digits) => digits.clone(),
        Json::Float(number) => pytext::float_repr(*number),
        Json::Str(text) => text.clone(),
        Json::List(items) => {
            let items: Vec<String> = items.iter().map(python_repr).collect();
            format!("[{}]", items.join(", "))
        }
        Json::Object(members) => {
            let members: Vec<String> = members
                .iter()
                .map(|(key, item)| format!("{}: {}", pytext::repr(key), python_repr(item)))
                .collect();
            format!("{{{}}}", members.join(", "))
        }
    }
}

/// The groups of a `5UTR_annotation` as maps of key to value, one group per upstream ORF:
/// the text formats' `key=value` pairs joined by `:`, with `,` between groups (`&` in a
/// VCF value), or the numbered objects of the JSON format's object, each member's value as
/// `str()` writes it. `None` for a value that is neither, an empty one included.
pub fn annotation_frames(value: &str) -> Option<Vec<Frame>> {
    if value.starts_with('{') {
        let Ok(Json::Object(groups)) = pytext::json_loads(value) else {
            return None;
        };
        return groups
            .iter()
            .map(|(_, group)| match group {
                Json::Object(members) => Some(
                    members
                        .iter()
                        .map(|(key, item)| (key.clone(), python_str(item)))
                        .collect(),
                ),
                _ => None,
            })
            .collect();
    }
    if blank(value) || (value.contains(',') && value.contains('&')) {
        return None;
    }
    let separator = if value.contains('&') { '&' } else { ',' };
    value
        .split(separator)
        .map(|group| {
            let mut frame: Frame = Vec::new();
            for pair in group.split(':') {
                let (key, sign, item) = pytext::partition(pair, "=");
                if sign.is_empty() {
                    return None;
                }
                match frame.iter_mut().find(|(held, _)| held == key) {
                    Some(held) => held.1 = item.to_string(),
                    None => frame.push((key.to_string(), item.to_string())),
                }
            }
            Some(frame)
        })
        .collect()
}

/// The size of the row's `HGVS_OFFSET`, the bases between the record's own insertion point
/// and the last index its mutated sequence can be written at; 0 where the row writes no
/// whole number there, `None` for a number outside `i128`.
fn hgvs_shift(e: &Entry) -> Option<i128> {
    let offset = field_of(e, &HGVS_OFFSET_FIELDS);
    if is_integer(offset) {
        pytext::int(offset)?.checked_abs()
    } else {
        Some(0)
    }
}

/// The 5' UTR base directly after the insertion of the entry's row as its record writes
/// it, in bases before the coding start, read from the row's HGVSc and `HGVS_OFFSET`. An
/// `ins` between two neighbouring 5' UTR positions names that base as its second
/// coordinate; a `dup` stands for the copy written directly after the duplicated bases, so
/// the base is the one after its last coordinate. The notation stands 3' of the record's
/// own insertion point on the transcript by the size of `HGVS_OFFSET`, which is printed
/// negative on a reverse-strand transcript (OutputFactory.pm:1764), so the record's point
/// lies that many bases further from the coding start. `None` for any other notation, for
/// a coordinate outside the 5' UTR and for a number outside `i128`.
pub fn utr_base_after_insertion(e: &Entry) -> Option<i128> {
    let notation = pytext::partition(field_of(e, &HGVSC_FIELDS), ":").2;
    let shift = hgvs_shift(e)?;
    if let Some(inserted) = UTR_HGVSC_INSERTION_RE.captures(notation) {
        let first = pytext::int(inserted.get(1)?.as_str())?;
        let last = pytext::int(inserted.get(2)?.as_str())?;
        if last != first.checked_add(1)? {
            return None;
        }
        return last.checked_neg()?.checked_add(shift);
    }
    let duplicated = UTR_HGVSC_DUPLICATION_RE.captures(notation)?;
    let last = pytext::int(duplicated.get(2).or(duplicated.get(1))?.as_str())?;
    last.checked_neg()?.checked_sub(1)?.checked_add(shift)
}

/// One group of an annotation with its pairs in key order.
type Group<'a> = Vec<(&'a str, &'a str)>;

/// The groups of an annotation, each with its pairs in key order, with the number of times
/// the annotation lists it.
fn counted(frames: &[Frame]) -> BTreeMap<Group<'_>, usize> {
    let mut groups = BTreeMap::new();
    for frame in frames {
        let mut pairs: Group = frame
            .iter()
            .map(|(key, value)| (key.as_str(), value.as_str()))
            .collect();
        pairs.sort_unstable();
        *groups.entry(pairs).or_insert(0) += 1;
    }
    groups
}

/// The value of `key` in a group.
fn pair_of<'a>(group: &Group<'a>, key: &str) -> Option<&'a str> {
    group
        .iter()
        .find_map(|(held, value)| (*held == key).then_some(*value))
}

/// The groups an entry lists beside the frameshift: those of its `5UTR_annotation` where
/// its `5UTR_consequence` is the frameshift, none where it fills neither field. `None` for
/// any other entry and for an annotation `annotation_frames` does not read.
fn shifted_frames(e: &Entry) -> Option<Vec<Frame>> {
    let annotation = field_of(e, &ANNOTATION_FIELDS);
    let term = field_of(e, &CONSEQUENCE_FIELDS);
    if term == FRAMESHIFT_TERM {
        annotation_frames(annotation)
    } else if blank(term) && blank(annotation) {
        Some(Vec::new())
    } else {
        None
    }
}

/// The first base of a listed uORF's stop codon, `ref_StartDistanceToCDS` less
/// `ref_type_length` plus three bases before the coding start, with whether the group
/// prints an `alt_type_length` equal to its `ref_type_length`. `None` for a group that is
/// no uORF of the reference UTR and for one whose two numbers are not whole numbers.
fn listed_stop_codon(group: &Group) -> Option<(i128, bool)> {
    let start = pair_of(group, "ref_StartDistanceToCDS")?;
    let length = pair_of(group, "ref_type_length")?;
    if pair_of(group, "ref_type") != Some("uORF") {
        return None;
    }
    let stop = digits_value(start)?
        .checked_sub(digits_value(length)?)?
        .checked_add(3)?;
    Some((stop, pair_of(group, "alt_type_length") == Some(length)))
}

/// Whether the two entries' frameshift lists differ only by upstream ORFs at whose stop
/// codon the insertion stands (`utrannotator_frameshift_at_stop_boundary`). Each entry
/// lists the groups of its `5UTR_annotation` beside the frameshift, or none where it fills
/// neither field (`shifted_frames`), and one list holds every group of the other and more.
/// UTRAnnotator.pm:921 holds the record's own index against the check point one index
/// before an ORF's stop codon (:911-917) and reads no base. A row that states no insertion
/// point is refused (`utr_base_after_insertion`, read from the reference's row).
///
/// Where the engine lists more, an insertion written directly before the stop codon's
/// first base is past the check point and the ORF is left out, though it runs on past its
/// stop codon's index and its length changes. Each group the engine adds is a uORF whose
/// two printed lengths differ and whose stop codon's first base (`listed_stop_codon`) is
/// the base directly after the insertion. No group both entries list may be such an ORF:
/// there the reference's own test admits what it leaves out here.
///
/// Where the reference lists more, the record's index is at or before the check point and
/// the ORF is listed, though the same mutated sequence written further on leaves the stop
/// codon at its own index and the ORF its length. Each group the reference adds is a uORF
/// that stays one (`alt_type=uORF`), whose two printed lengths are equal and whose stop
/// codon's first base lies at or after the base directly after the insertion and no more
/// than one base after the base that follows the last writing of the insertion, which
/// stands nearer the coding start by the size of `HGVS_OFFSET`. An engine's entry that
/// fills neither field carries the plugin's counts.
pub fn frameshift_at_stop_boundary_row(e: &Entry, a: &Entry) -> bool {
    let (Some(expected), Some(actual)) = (shifted_frames(e), shifted_frames(a)) else {
        return false;
    };
    if actual.is_empty() && blank(field_of(a, &COUNT_FIELDS)) {
        return false;
    }
    let (listed, written) = (counted(&expected), counted(&actual));
    let beyond = |more: &BTreeMap<Group, usize>, fewer: &BTreeMap<Group, usize>| -> bool {
        more.iter()
            .any(|(group, count)| *count > fewer.get(group).copied().unwrap_or(0))
    };
    let (adds, drops) = (beyond(&written, &listed), beyond(&listed, &written));
    if adds == drops {
        return false;
    }
    let Some(after) = utr_base_after_insertion(e) else {
        return false;
    };
    if adds {
        let ends_at_the_insertion =
            |group: &Group| listed_stop_codon(group) == Some((after, false));
        return written
            .iter()
            .filter(|(group, count)| **count > listed.get(*group).copied().unwrap_or(0))
            .map(|(group, _)| group)
            .all(ends_at_the_insertion)
            && !written
                .keys()
                .filter(|group| listed.contains_key(*group))
                .any(ends_at_the_insertion);
    }
    let Some(last) = hgvs_shift(e).and_then(|shift| after.checked_sub(shift)) else {
        return false;
    };
    listed
        .iter()
        .filter(|(group, count)| **count > written.get(*group).copied().unwrap_or(0))
        .all(|(group, _)| {
            pair_of(group, "alt_type") == Some("uORF")
                && listed_stop_codon(group).is_some_and(|(stop, kept)| {
                    kept && stop <= after && last.checked_sub(1).is_some_and(|from| stop >= from)
                })
        })
}

/// The `KozakContext` and `KozakStrength` of a group.
fn kozak(group: &[(String, String)]) -> (Option<&str>, Option<&str>) {
    (
        value_of(group, "KozakContext"),
        value_of(group, "KozakStrength"),
    )
}

/// Whether the reference reports the frameshift and each of the engine's lost-stop groups
/// `lost` is an upstream ORF the reference lists among the shifted ones: for the same
/// mutated sequence written before the start codon that overlaps the stop codon the ORF is
/// reported shifted (UTRAnnotator.pm:921), with its new length. Each of the engine's groups
/// carries the `KozakContext` and `KozakStrength` of exactly one uORF of the reference's
/// list, a different one for each group, and the alternative stop codon where that ORF's
/// printed `alt_type_length` ends it: at `ref_StartDistanceToCDS` plus the allele's length
/// less `alt_type_length` plus three bases before the coding start where the ORF still
/// closes in the 5' UTR (`alt_type=uORF`), none otherwise.
fn lost_stop_of_shifted_orfs(e: &Entry, lost: &[Frame]) -> bool {
    if field_of(e, &CONSEQUENCE_FIELDS) != FRAMESHIFT_TERM {
        return false;
    }
    let Some(shifted) = annotation_frames(field_of(e, &ANNOTATION_FIELDS)) else {
        return false;
    };
    let mut taken: Vec<usize> = Vec::new();
    for group in lost {
        let named = kozak(group);
        let mut frames = shifted.iter().enumerate().filter(|(_, frame)| {
            value_of(frame, "ref_type") == Some("uORF") && kozak(frame) == named
        });
        let (Some((at, frame)), None) = (frames.next(), frames.next()) else {
            return false;
        };
        if named.0.is_none() || named.1.is_none() || taken.contains(&at) {
            return false;
        }
        taken.push(at);
        let stop = (
            value_of(group, "AltStop"),
            value_of(group, "AltStopDistanceToCDS"),
        );
        if value_of(frame, "alt_type") == Some("uORF") {
            let distance =
                digits_value(value_of(frame, "ref_StartDistanceToCDS").unwrap_or(ABSENT))
                    .zip(digits_value(
                        value_of(frame, "alt_type_length").unwrap_or(ABSENT),
                    ))
                    .and_then(|(start, length)| {
                        start
                            .checked_add(e.allele.chars().count() as i128)?
                            .checked_sub(length)?
                            .checked_add(3)
                    });
            let Some(distance) = distance else {
                return false;
            };
            if stop != (Some("True"), Some(distance.to_string().as_str())) {
                return false;
            }
        } else if stop != (Some("False"), Some("NA")) {
            return false;
        }
    }
    true
}

/// Whether the engine's lost-stop group `group` is of another upstream ORF than the
/// shifted ORF `frame` of the reference's list, whose start codon carries the same
/// `KozakContext` and `KozakStrength`. An ORF meets its next stop codon in its own frame,
/// and in the mutated UTR the listed ORF's start stands the allele's length further from
/// the coding start than its `ref_StartDistanceToCDS`. Two things hold of another ORF's
/// group. Its `AltStopDistanceToCDS` lies no whole number of codons past the listed ORF's
/// start. And it lies one codon or more, in whole codons, past the three bases read in
/// place of a stop codon that the insertion can be written directly before and inside
/// (UTRAnnotator.pm:639). Such a stop codon begins at the base after the record's own point
/// (`utr_base_after_insertion`) or at one nearer the coding start by less than the size of
/// `HGVS_OFFSET`, and the bases read in its place begin the allele's length further from
/// the coding start.
fn lost_stop_of_another_orf(e: &Entry, group: &Frame, frame: &Frame) -> bool {
    let length = e.allele.chars().count() as i128;
    let stop = digits_value(value_of(group, "AltStopDistanceToCDS").unwrap_or(ABSENT));
    let start = digits_value(value_of(frame, "ref_StartDistanceToCDS").unwrap_or(ABSENT));
    let (Some(stop), Some(start), Some(after), Some(shift)) =
        (stop, start, utr_base_after_insertion(e), hgvs_shift(e))
    else {
        return false;
    };
    let past = |from: i128| from.checked_add(length)?.checked_sub(stop);
    past(start).is_some_and(|bases| bases % 3 != 0)
        && (0..shift.min(3)).any(|nearer| {
            after
                .checked_sub(nearer)
                .and_then(past)
                .is_some_and(|bases| bases >= 3 && bases % 3 == 0)
        })
}

/// Whether the reference's entry leaves the lost stop of the engine's groups `lost` unread
/// on a record written at or before a stop codon's first base: UTRAnnotator.pm:631-633
/// pass over the stop codon of such a record, though the same mutated sequence can be
/// written inside the stop codon. The reference fills neither field or reports the lost
/// start, the gained start or the frameshift, the frameshift with no group that carries
/// the `KozakContext` and `KozakStrength` of a group of the engine's unless that group is
/// another ORF's (`lost_stop_of_another_orf`; `lost_stop_of_shifted_orfs` reads the row
/// of a listed ORF); the engine lists a group, and each of
/// its groups carries both keys. The allele's length is no multiple of three, so the ORF
/// does not end on the stop codon that stands whole after the inserted bases; the row's
/// `HGVS_OFFSET` is a whole number other than zero, so the sequence can be written further
/// on; and no run of `sweep` writes the lost stop on the row.
fn lost_stop_written_further_on(
    e: &Entry,
    lost: &[Frame],
    sweep: Option<&UtrannotatorSeedSweep>,
) -> bool {
    let named = |group: &Frame| {
        let (context, strength) = kozak(group);
        context.is_some() && strength.is_some()
    };
    if lost.is_empty() || !lost.iter().all(named) {
        return false;
    }
    let term = field_of(e, &CONSEQUENCE_FIELDS);
    let annotation = field_of(e, &ANNOTATION_FIELDS);
    let unread = if term == FRAMESHIFT_TERM {
        annotation_frames(annotation).is_some_and(|shifted| {
            !shifted.iter().any(|frame| {
                lost.iter().any(|group| {
                    kozak(group) == kozak(frame) && !lost_stop_of_another_orf(e, group, frame)
                })
            })
        })
    } else if blank(term) {
        blank(annotation)
    } else {
        term == LOST_START_TERM || term == GAINED_START_TERM
    };
    let length = e.allele.chars().count() as i128;
    unread
        && length % 3 != 0
        && hgvs_shift(e).is_some_and(|shift| shift != 0)
        && sweep.is_none_or(|sweep| !sweep.writes_term(e, LOST_STOP_TERM))
}

/// Whether the reference's lost stop is of a stop codon that stands whole after the
/// inserted bases, in frame. UTRAnnotator.pm:639-642 read the three bases at the stop
/// codon's index for a record written after its first or second base (:628-633); the same
/// mutated sequence can be written at the stop codon's first base, so the codon stands the
/// allele's length further on and, that length being a multiple of three, the ORF ends on
/// it or on a stop codon the inserted bases hold before it in frame. The reference's
/// annotation lists a group, and each of its groups reads `AltStop=True` with an
/// `AltStopDistanceToCDS` that names such a codon: more than the base directly after the
/// insertion (`utr_base_after_insertion`) by one base to the allele's length less one,
/// and by no multiple of three. The engine's `5UTR_consequence` is another of the
/// plugin's five terms, or both of its fields are empty.
fn stop_codon_stands_in_frame(e: &Entry, a: &Entry) -> bool {
    let length = e.allele.chars().count() as i128;
    if field_of(e, &CONSEQUENCE_FIELDS) != LOST_STOP_TERM
        || field_of(a, &CONSEQUENCE_FIELDS) == LOST_STOP_TERM
        || !term_or_none(a)
        || length % 3 != 0
    {
        return false;
    }
    let (Some(lost), Some(after)) = (
        annotation_frames(field_of(e, &ANNOTATION_FIELDS)),
        utr_base_after_insertion(e),
    ) else {
        return false;
    };
    !lost.is_empty()
        && lost.iter().all(|group| {
            value_of(group, "AltStop") == Some("True")
                && digits_value(value_of(group, "AltStopDistanceToCDS").unwrap_or(ABSENT))
                    .and_then(|distance| distance.checked_sub(after))
                    .is_some_and(|past| (1..length).contains(&past) && past % 3 != 0)
        })
}

/// Whether the two entries of an insertion's row differ over a lost stop as the
/// reference's reading of a stop codon's index makes them
/// (`utrannotator_lost_stop_written_before_stop_codon`). UTRAnnotator.pm:628-633 read a
/// stop codon's index only for a record written inside the stop codon, and :639-642 then
/// read only the three bases at that index. The row's cDNA positions are an insertion's
/// two flanking bases, the engine's entry carries the plugin's counts, and where the
/// row's record states its alleles the reference allele is empty. The row is the class's
/// in one of three shapes:
///
/// - the engine reports the lost stop of upstream ORFs the reference lists among the
///   shifted ones (`lost_stop_of_shifted_orfs`);
/// - the engine reports the lost stop on a record written at or before the stop codon's
///   first base, where the reference reports none, the lost start, the gained start or the
///   frameshift of other ORFs (`lost_stop_written_further_on`, which reads `sweep`);
/// - the reference reports the lost stop of a stop codon that stands whole in frame after
///   the inserted bases, and the engine reports no lost stop
///   (`stop_codon_stands_in_frame`).
pub fn lost_stop_written_before_stop_codon_row(
    e: &Entry,
    a: &Entry,
    record_rows: Option<&[Row]>,
    sweep: Option<&UtrannotatorSeedSweep>,
) -> bool {
    if blank(field_of(a, &COUNT_FIELDS)) || insertion_index(e).is_none() || !is_sequence(&e.allele)
    {
        return false;
    }
    let differs = if field_of(a, &CONSEQUENCE_FIELDS) == LOST_STOP_TERM {
        annotation_frames(field_of(a, &ANNOTATION_FIELDS)).is_some_and(|lost| {
            lost_stop_of_shifted_orfs(e, &lost) || lost_stop_written_further_on(e, &lost, sweep)
        })
    } else {
        stop_codon_stands_in_frame(e, a)
    };
    differs && reference_allele_is_empty_where_stated(e, record_rows)
}

/// Whether two `5UTR_annotation` texts carry the same groups of the same `key=value` pairs
/// with the groups, or the pairs inside one, in different orders, the engine's pairs in
/// the order of their keys inside every group (`utrannotator_annotation_pair_order`):
/// UTRAnnotator.pm:282-283 walk the groups and each group's pairs as hashes, and the engine
/// writes each group's pairs in the order of their keys.
pub fn reordered_pair_groups(expected: &str, actual: &str) -> bool {
    let (expected_parts, actual_parts) = list_parts(expected, actual);
    let groups = |parts: &[&'_ str]| -> Vec<Vec<String>> {
        parts
            .iter()
            .map(|group| group.split(':').map(str::to_string).collect())
            .collect()
    };
    let (mut expected_groups, mut actual_groups) = (groups(&expected_parts), groups(&actual_parts));
    if expected_groups == actual_groups {
        return false;
    }
    if !expected_groups
        .iter()
        .chain(&actual_groups)
        .flatten()
        .all(|pair| pair.contains('='))
    {
        return false;
    }
    let key = |pair: &str| -> String { pytext::partition(pair, "=").0.to_string() };
    if actual_groups.iter().any(|group| {
        group
            .windows(2)
            .any(|pairs| key(&pairs[0]) > key(&pairs[1]))
    }) {
        return false;
    }
    for group in expected_groups.iter_mut().chain(actual_groups.iter_mut()) {
        group.sort_unstable();
    }
    expected_groups.sort_unstable();
    actual_groups.sort_unstable();
    expected_groups == actual_groups
}

/// `utrannotator_gained_start_window` on a consult.
pub fn gained_start_window(c: &Consult) -> bool {
    c.a.is_some_and(|a| gained_start_window_row(c.e, a, c.record_rows))
}

/// `utrannotator_consequence_precedence` on a consult: UTRAnnotator.pm:262-267 walk the
/// flags of the plugin's analyses as a hash and keep the last that fired, so where several
/// fire the reported effect follows the process's hash order. Two terms on a row do not
/// show that two analyses fire there, so the class accepts nothing without the seeded
/// runs: both terms are ones those runs write on the row, and each side's
/// `5UTR_annotation` is one they write beside its term. An annotation the two entries
/// share is no differing field of its own, so the consequence is accepted only where the
/// runs write it beside both terms.
pub fn consequence_precedence(c: &Consult) -> bool {
    let (Some(a), Some(sweep)) = (c.a, c.utrannotator_sweep) else {
        return false;
    };
    if !precedence_row(c.e, a) {
        return false;
    }
    let annotations = [
        field_of(c.e, &ANNOTATION_FIELDS),
        field_of(a, &ANNOTATION_FIELDS),
    ];
    let shared = !annotations.contains(&ABSENT)
        && annotation_groups(annotations[0]) == annotation_groups(annotations[1]);
    sweep.writes(c.e, a, shared || field_named(c.name, &ANNOTATION_FIELDS))
}

/// `utrannotator_structural_allele_text` on a consult, `ev` and `av` the two values with a
/// blank one as `-`: the fields the reference alone fills on a structural class word's
/// row, or the fields the engine fills under `insertion`, under `inversion`, under
/// `copy_number_variation` or under a deletion's class word.
pub fn structural_allele_text(c: &Consult, ev: &str, av: &str) -> bool {
    let Some(a) = c.a else {
        return false;
    };
    if blank(av) {
        return !blank(ev) && structural_allele_text_row(c.e, a);
    }
    match c.e.allele.as_str() {
        STRUCTURAL_INSERTION_WORD => structural_insertion_length_row(c.e, a),
        STRUCTURAL_INVERSION_WORD => structural_inversion_span_row(c.e, a, c.record_rows),
        STRUCTURAL_COPY_NUMBER_WORD => structural_copy_number_loss_row(c.e, a, c.record_rows),
        _ => structural_deletion_span_row(c.e, a),
    }
}

/// `utrannotator_lost_start_after_insertion` on a consult.
pub fn lost_start_after_insertion(c: &Consult) -> bool {
    c.a.is_some_and(|a| lost_start_after_insertion_row(c.e, a, c.record_rows))
}

/// `utrannotator_gained_start_remade_in_place` on a consult.
pub fn gained_start_remade_in_place(c: &Consult) -> bool {
    c.a.is_some_and(|a| {
        gained_start_remade_in_place_row(c.e, a, c.record_rows, c.utrannotator_sweep)
    })
}

/// `utrannotator_frameshift_at_stop_boundary` on a consult: the annotation of a row the
/// class accepts, and its consequence where one entry reports the frameshift and the
/// other nothing.
pub fn frameshift_at_stop_boundary(c: &Consult) -> bool {
    let Some(a) = c.a else {
        return false;
    };
    let reports = |entry: &Entry| !blank(field_of(entry, &CONSEQUENCE_FIELDS));
    if field_named(c.name, &CONSEQUENCE_FIELDS) && reports(c.e) == reports(a) {
        return false;
    }
    frameshift_at_stop_boundary_row(c.e, a)
}

/// `utrannotator_lost_stop_written_before_stop_codon` on a consult.
pub fn lost_stop_written_before_stop_codon(c: &Consult) -> bool {
    c.a.is_some_and(|a| {
        lost_stop_written_before_stop_codon_row(c.e, a, c.record_rows, c.utrannotator_sweep)
    })
}
