// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The written-span gates of a shifted indel and the intron offsets abutting an insertion.

use crate::consts::CDNA_POSITION_FIELDS;
use crate::consts::CDS_POSITION_FIELDS;
use crate::consts::DELETION_NOTATION_RE_PATTERN;
use crate::consts::FLAGS_FIELDS;
use crate::consts::HGVSC_FIELDS;
use crate::consts::HGVS_INTRON_OFFSET_FIELDS;
use crate::consts::HGVS_OFFSET_FIELDS;
use crate::consts::INDEL_NOTATION_RE_PATTERN;
use crate::consts::INTEGER_RE_PATTERN;
use crate::consts::INTRON_OFFSET_RE;
use crate::consts::NON_ZERO_INTEGER_RE_PATTERN;
use crate::consts::OUTSIDE_TRANSCRIPT_TERMS;
use crate::consts::PROTEIN_POSITION_FIELDS;
use crate::consts::REFSEQ_OFFSET_FIELDS;
use crate::consts::SHIFT_LENGTH_FIELDS;
use crate::consts::STRAND_FIELDS;
use crate::entry::Entry;
use crate::error::Result;
use crate::fields::field_named;
use crate::pyfmt;
use crate::rules::clipped::insertion_left_where_the_shift_stops;
use crate::rules::clipped::parted;
use crate::rules::downstream::field_of;
use crate::rules::downstream::positive_shift;
use crate::rules::downstream::shift_field_consistent;
use crate::rules::edited_model::bam_edit_among;
use crate::rules::hgvsp_stop::STOP_EXTENSION_RE_PATTERN;
use crate::rules::hgvsp_stop::STOP_UNCHANGED_RE_PATTERN;
use crate::rules::shifted_start_test::start_test_read_on_the_written_span;
use crate::rules::shifted_start_test::ShiftedRows;
use crate::rules::start_test_past_the_allele::start_codon_test_past_the_shifted_allele;
use crate::rules::start_test_past_the_allele::StartTestRows;
use crate::rules::transcript_edges::TranscriptEdges;
use crate::rules::{first_field, nearest_i64, whole_number, whole_number_text, whole_text};
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::indel_row;
use crate::shared::insertion_row;
use crate::shared::integer_of;
use crate::shared::position_bounds;
use crate::shared::ABSENT;
use std::sync::LazyLock;

pub const SHIFTED_INDEL_GATE_FIELDS: &[&str] = &[
    "HGVSc",
    "hgvsc",
    "HGVS_OFFSET",
    "hgvs_offset",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];

pub const PROTEIN_NOTATION_RE_PATTERN: &str = r"[^:\s\x1c-\x1f]+:p\.[^\s\x1c-\x1f]+";
pub static PROTEIN_NOTATION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(PROTEIN_NOTATION_RE_PATTERN).expect("the pattern compiles"));

pub const STOP_CODON_NOTATION_RE_PATTERN: &str = r"[^:\s\x1c-\x1f]+:p\.\(?(?:Ter|\*)\d";
pub static STOP_CODON_NOTATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(STOP_CODON_NOTATION_RE_PATTERN).expect("the pattern compiles")
});

pub const STOP_CODON_FLANK_INSERTION_RE_PATTERN: &str =
    r"[^:\s\x1c-\x1f]+:p\.[^_\s\x1c-\x1f]+_(?:Ter|\*)\d+ins";
pub static STOP_CODON_FLANK_INSERTION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(STOP_CODON_FLANK_INSERTION_RE_PATTERN).expect("the pattern compiles")
});

/// A duplicated, deleted or replaced range that ends at the stop codon holds the stop's residue,
/// so it changes the stop codon as a notation that starts there does.
pub const STOP_CODON_RANGE_END_RE_PATTERN: &str =
    r"[^:\s\x1c-\x1f]+:p\.[^_\s\x1c-\x1f]+_(?:Ter|\*)\d+(?:dup|del)";
pub static STOP_CODON_RANGE_END_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(STOP_CODON_RANGE_END_RE_PATTERN).expect("the pattern compiles")
});

/// The flag of a coding sequence whose end is not found: Ensembl VEP reads no insertion or
/// deletion as overlapping the stop codon on such a transcript (VariationEffect.pm:1373,
/// 1400).
pub const CDS_END_FLAG: &str = "cds_end_NF";

/// The term of an allele that lies on the transcript before its coding region.
pub const FIVE_PRIME_UTR_TERM: &str = "5_prime_UTR_variant";

pub const CODING_COORDINATE_RE_PATTERN: &str = r"[1-9]\d*";
pub static CODING_COORDINATE_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(CODING_COORDINATE_RE_PATTERN).expect("the pattern compiles")
});

pub const ACCEPTOR_INTRON_COORDINATE_RE_PATTERN: &str = r"(\*|[1-9]\d*)-([1-9]\d*)";
pub static ACCEPTOR_INTRON_COORDINATE_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(ACCEPTOR_INTRON_COORDINATE_RE_PATTERN).expect("the pattern compiles")
});

pub const DONOR_INTRON_COORDINATE_RE_PATTERN: &str = r"([1-9]\d*)\+[1-9]\d*";
pub static DONOR_INTRON_COORDINATE_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(DONOR_INTRON_COORDINATE_RE_PATTERN).expect("the pattern compiles")
});

/// `pattern` held to the start of a text, as `re.match` holds it.
fn opening(pattern: &str) -> regex::Regex {
    regex::Regex::new(&format!(r"\A(?:{pattern})")).expect("the pattern compiles")
}

pub(crate) static INDEL_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INDEL_NOTATION_RE_PATTERN));
pub(crate) static DELETION_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(DELETION_NOTATION_RE_PATTERN));
pub(crate) static NON_ZERO_INTEGER_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NON_ZERO_INTEGER_RE_PATTERN));
static INTEGER_WHOLE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(INTEGER_RE_PATTERN));
static PROTEIN_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(PROTEIN_NOTATION_RE_PATTERN));
static STOP_CODON_OPENINGS: LazyLock<[regex::Regex; 3]> = LazyLock::new(|| {
    [
        opening(STOP_CODON_NOTATION_RE_PATTERN),
        opening(STOP_CODON_FLANK_INSERTION_RE_PATTERN),
        opening(STOP_CODON_RANGE_END_RE_PATTERN),
    ]
});
static STOP_UNCHANGED_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_UNCHANGED_RE_PATTERN));
static STOP_EXTENSION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_EXTENSION_RE_PATTERN));
static CODING_COORDINATE_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(CODING_COORDINATE_RE_PATTERN));
static ACCEPTOR_INTRON_COORDINATE_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(ACCEPTOR_INTRON_COORDINATE_RE_PATTERN));
static DONOR_INTRON_COORDINATE_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(DONOR_INTRON_COORDINATE_RE_PATTERN));
/// The position and the intron offset of a coordinate without its region sign.
static COORDINATE_PARTS_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(r"(\d*)([+-]\d+)?"));
/// A coordinate before the coding start: a base of the 5' UTR, alone or with the offset of
/// an intron base counted from it.
static UTR_COORDINATE_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(r"-[1-9]\d*(?:[+-][1-9]\d*)?"));

/// The number a run of the digits `0` to `9` writes plus `addend`, as decimal text, at any
/// length.
fn plus(digits: &str, addend: usize) -> String {
    let mut carry = addend as u128;
    let mut sum: Vec<u8> = Vec::with_capacity(digits.len() + 1);
    for digit in digits.bytes().rev() {
        carry += u128::from(digit - b'0');
        sum.push(b'0' + (carry % 10) as u8);
        carry /= 10;
    }
    while carry > 0 {
        sum.push(b'0' + (carry % 10) as u8);
        carry /= 10;
    }
    sum.reverse();
    String::from_utf8(sum).expect("decimal digits")
}

/// Whether one of the entry's terms is `upstream_gene_variant` or `downstream_gene_variant`.
fn annotated_outside_transcript(e: &Entry) -> bool {
    consequence_set(&e.consequence)
        .split(',')
        .any(|term| OUTSIDE_TRANSCRIPT_TERMS.contains(&term))
}

/// Whether one of the entry's terms is `term`.
fn names_term(e: &Entry, term: &str) -> bool {
    consequence_set(&e.consequence)
        .split(',')
        .any(|named| named == term)
}

/// Whether a notation of `kind` names an insertion: a `dup` or an `ins`.
fn names_insertion(kind: &str) -> bool {
    matches!(kind, "dup" | "ins")
}

/// The kind of an HGVSc that names the entry's transcript (its stable id, with or without a
/// version) and describes an insertion or deletion: `del`, `dup` or `ins`; `None` for
/// any other notation.
pub fn indel_notation_of(e: &Entry, hgvsc: &str) -> Option<String> {
    let (name, _, body) = pyfmt::partition(hgvsc, ":");
    let stable_id = |text: &str| text.split('.').next().map(str::to_string);
    if stable_id(name) != stable_id(&e.feature) {
        return None;
    }
    let notation = INDEL_NOTATION_WHOLE.captures(body)?;
    Some(notation[3][..3].to_string())
}

/// `hgvs_shifted_indel_written_span_gate` on one differing field: the reference's HGVSc,
/// HGVS_OFFSET and HGVSIntronOffset fields are all empty, the candidate's HGVSc is an insertion
/// or deletion notation of the row's transcript beside a non-zero HGVS_OFFSET, the differing
/// field is one of those (an intron offset an integer, and beside a `del` notation the one
/// its coordinates give: `deletion_offsets_follow_notation`), and the row is an indel's
/// annotated inside the transcript whose SHIFT_LENGTH, where printed, is positive. An allele of two
/// or more
/// bases at a Location of as many, under an identifier that does not spell the record's alleles,
/// reads the same for an insertion and a substitution in the text formats (`indel_row` refuses
/// it); there the candidate's `dup` or `ins` notation names the row an insertion's. On a row
/// carrying `BAM_EDIT` `OK` on both sides the differing field may also be the REFSEQ_OFFSET,
/// the candidate's a non-zero integer: `hgvs_transcript` returns at
/// TranscriptVariationAllele.pm:1472, before its call of `get_misalignment_offset` (:1517)
/// stores the offset (:1621), so the reference's row prints none (OutputFactory.pm:1778).
pub fn shifted_indel_written_span_gate(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> Result<bool> {
    if !blank(ev) || blank(av) {
        return Ok(false);
    }
    let Some(actual) = actual else {
        return Ok(false);
    };
    if !shift_field_consistent(e)? {
        return Ok(false);
    }
    if SHIFTED_INDEL_GATE_FIELDS
        .iter()
        .any(|field| !blank(&field_of(Some(e), &[*field])))
    {
        return Ok(false);
    }
    // A row annotated outside the transcript has no notation on either engine: one the
    // candidate writes there is a difference, not this class.
    if annotated_outside_transcript(e) {
        return Ok(false);
    }
    let hgvsc = if field_named(name, HGVSC_FIELDS) {
        av.to_string()
    } else {
        field_of(Some(actual), HGVSC_FIELDS)
    };
    let offset = if field_named(name, HGVS_OFFSET_FIELDS) {
        av.to_string()
    } else {
        field_of(Some(actual), HGVS_OFFSET_FIELDS)
    };
    let intron_offset = field_named(name, HGVS_INTRON_OFFSET_FIELDS);
    if intron_offset && !INTEGER_WHOLE.is_match(av) {
        return Ok(false);
    }
    if field_named(name, REFSEQ_OFFSET_FIELDS)
        && !(NON_ZERO_INTEGER_WHOLE.is_match(av)
            && [e, actual].iter().all(|row| bam_edit_among(row, &["OK"])))
    {
        return Ok(false);
    }
    let Some(kind) = indel_notation_of(e, &hgvsc) else {
        return Ok(false);
    };
    if !NON_ZERO_INTEGER_WHOLE.is_match(&offset) {
        return Ok(false);
    }
    if intron_offset && kind == "del" && !deletion_offsets_follow_notation(actual, &hgvsc) {
        return Ok(false);
    }
    Ok(indel_row(e)? || (names_insertion(&kind) && equal_span_allele_row(e)))
}

/// The intron offset of an HGVS cDNA coordinate as the HGVSIntronOffset plugin writes it:
/// `9` for `29+9`, `-944` for `30-944`, `0` for an exonic coordinate.
pub fn coordinate_intron_offset(coordinate: &str) -> String {
    match INTRON_OFFSET_RE.find(coordinate).ok().flatten() {
        Some(offset) => whole_number_text(offset.as_str()),
        None => "0".to_string(),
    }
}

/// Whether the candidate row's two HGVSIntronOffset fields are the intron offsets of the
/// two coordinates of its `del` notation `hgvsc`, one coordinate standing for both, in
/// either order: the plugin reads them before the coordinates are ordered for printing. A
/// coordinate before the stop codon's last base (`*-5`) prints an intron distance without
/// its exon base, so a notation carrying one says nothing here.
pub fn deletion_offsets_follow_notation(a: &Entry, hgvsc: &str) -> bool {
    let Some(notation) = DELETION_NOTATION_WHOLE.captures(pyfmt::partition(hgvsc, ":").2) else {
        return false;
    };
    let first = &notation[2];
    let last = notation.get(3).map_or(first, |end| end.as_str());
    if first.starts_with("*-") || last.starts_with("*-") {
        return true;
    }
    let mut written: Vec<String> = HGVS_INTRON_OFFSET_FIELDS
        .iter()
        .map(|field| field_of(Some(a), &[*field]))
        .collect();
    written.sort();
    let mut read = [
        coordinate_intron_offset(first),
        coordinate_intron_offset(last),
    ];
    read.sort();
    written == read
}

/// Whether a protein notation names the stop codon alone and keeps every residue before it:
/// the stop written as unchanged (`Ter84=`, the `=` as the text formats escape it too) or an
/// extension named at the stop (`Ter84GlnextTer17`, the count a number or `?`), in prediction
/// parentheses or not.
fn keeps_every_residue_before_the_stop(notation: &str) -> bool {
    let (_, separator, description) = pyfmt::partition(notation, ":p.");
    if separator.is_empty() {
        return false;
    }
    let description = description.replace("%3D", "=");
    let description = description
        .strip_prefix('(')
        .and_then(|inner| inner.strip_suffix(')'))
        .unwrap_or(&description);
    STOP_UNCHANGED_WHOLE.is_match(description) || STOP_EXTENSION_WHOLE.is_match(description)
}

/// `hgvsp_shifted_indel_written_span_gate` on a differing HGVSp: the reference's is empty
/// and the candidate's a protein notation, the row is an indel's annotated inside the
/// transcript with no CDS or protein position and no positive SHIFT_LENGTH, both engines
/// write the same HGVSc and the same non-zero HGVS_OFFSET, or the reference's are an
/// insertion's left where its shift stops against a candidate that describes it further 3'
/// or not at all (`insertion_left_where_the_shift_stops`, read without the cache's exon
/// spans: a `dup` by its count of bases alone), and the reference's HGVSc is an insertion or
/// deletion notation of the row's transcript that places the shifted allele in the coding
/// sequence. The allele as the record writes it lies where the predicate Ensembl VEP tests
/// does not read `coding`: in an intron inside the coding region (`non_coding`,
/// BaseVariationFeatureOverlapAllele.pm:535-537) or on the transcript before its coding
/// start (`utr`, :523-525). The row shows that place in one of two ways.
///
/// A row with no cDNA position is an allele written in an intron
/// (`shifted_out_of_an_intron`). Where the notation names or bounds the first base of the
/// exon the allele is shifted into as coding base 2 or later (`shifted_exon_start`), the
/// intron lies inside the coding region and the row is taken from the row alone. So is a
/// `dup` that ends on a coding base and begins before the coding start, with or without
/// `edges`: the allele as written stands in the copy or before its first base, in an intron
/// of the 5' UTR or in one inside the coding region that the copy spans. A notation of
/// coding bases alone, or with an intron coordinate counted from a coding base, whose bound
/// is below 2 is taken where `edges` holds the candidate cache's exon spans and the allele
/// as the row's Location writes it lies wholly in the intron before the exon that holds the
/// shifted allele (`placed_exon_start`), whatever that exon's first coding base. Without
/// `edges` such a row stays a difference.
///
/// A row with a cDNA position is an allele written on exon bases, taken where one of the
/// row's terms is `5_prime_UTR_variant`. Its notation is read from the row alone
/// (`reaches_the_coding_sequence`), and `edges` is not read.
///
/// A notation at the stop codon is taken in two forms, the stop written as unchanged and an
/// extension named at the stop (`keeps_every_residue_before_the_stop`): any other notation
/// that starts at the stop codon, in parentheses or not, and one that names it as a flank of
/// an insertion or as the end of a duplicated, deleted or replaced range, stay differences.
/// A row either engine flags `cds_start_NF` is read as a row without a flag: the candidate
/// writes the notation of the shifted change on a coding sequence that lacks its start and
/// ends on its stop codon. A row either engine flags `cds_end_NF` is taken with any notation
/// but those two forms at the stop codon: the candidate writes the notation of a shifted
/// change that the codons of such a sequence give, and reads its stop codon as Ensembl VEP
/// does (`CDS_END_FLAG`).
/// The two-base allele at a two-base Location is read as `shifted_indel_written_span_gate`
/// reads it.
pub fn shifted_indel_written_span_protein_gate(
    e: &Entry,
    actual: Option<&Entry>,
    ev: &str,
    av: &str,
    edges: Option<&TranscriptEdges>,
) -> Result<bool> {
    if !blank(ev) || blank(av) {
        return Ok(false);
    }
    let Some(actual) = actual else {
        return Ok(false);
    };
    if positive_shift(e)? {
        return Ok(false);
    }
    if !PROTEIN_NOTATION_WHOLE.is_match(av) {
        return Ok(false);
    }
    let stop_kept = keeps_every_residue_before_the_stop(av);
    if !stop_kept && STOP_CODON_OPENINGS.iter().any(|stop| stop.is_match(av)) {
        return Ok(false);
    }
    if stop_kept
        && [e, actual]
            .iter()
            .any(|row| field_of(Some(row), FLAGS_FIELDS).contains(CDS_END_FLAG))
    {
        return Ok(false);
    }
    if annotated_outside_transcript(e) {
        return Ok(false);
    }
    let positioned = |fields: &[&str]| {
        fields
            .iter()
            .any(|field| !blank(&field_of(Some(e), &[*field])))
    };
    if positioned(CDS_POSITION_FIELDS) || positioned(PROTEIN_POSITION_FIELDS) {
        return Ok(false);
    }
    // A cDNA position without a CDS position is an allele on exon bases outside the coding
    // region, where the predicate reads `utr` (BaseVariationFeatureOverlapAllele.pm:523-525).
    // The term tells the 5' UTR, out of which a 3' shift reaches the coding sequence.
    let on_exon_bases = positioned(CDNA_POSITION_FIELDS);
    if on_exon_bases && !names_term(e, FIVE_PRIME_UTR_TERM) {
        return Ok(false);
    }
    let hgvsc = field_of(Some(e), HGVSC_FIELDS);
    let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
    let alike = hgvsc == field_of(Some(actual), HGVSC_FIELDS)
        && offset == field_of(Some(actual), HGVS_OFFSET_FIELDS);
    if blank(&hgvsc)
        || !(alike || insertion_left_where_the_shift_stops(e, Some(actual), HGVSC_FIELDS[0], None)?)
    {
        return Ok(false);
    }
    let Some(kind) = indel_notation_of(e, &hgvsc) else {
        return Ok(false);
    };
    if !NON_ZERO_INTEGER_WHOLE.is_match(&offset) {
        return Ok(false);
    }
    let shift = whole_number(&offset).saturating_abs();
    let body = pyfmt::partition(&hgvsc, ":").2;
    let shifted_into_the_coding_sequence = if on_exon_bases {
        reaches_the_coding_sequence(body)
    } else {
        shifted_out_of_an_intron(e, body, shift, edges)?
    };
    if !shifted_into_the_coding_sequence {
        return Ok(false);
    }
    Ok(indel_row(e)? || (names_insertion(&kind) && equal_span_allele_row(e)))
}

/// `hgvsp_shifted_indel_written_span_start_test` on a differing HGVSp: the two rows' terms,
/// HGVS_OFFSET and CDS position and the candidate's HGVSc give the two values the class's
/// shape (`start_test_read_on_the_written_span`), on an insertion's or deletion's row. The
/// two-base allele at a two-base Location is read as `shifted_indel_written_span_gate` reads
/// it, from the candidate's notation. A comparison that holds no candidate row is refused.
pub fn shifted_indel_written_span_start_test(
    e: &Entry,
    actual: Option<&Entry>,
    ev: &str,
    av: &str,
) -> Result<bool> {
    let Some(actual) = actual else {
        return Ok(false);
    };
    let hgvsc = field_of(Some(actual), HGVSC_FIELDS);
    let insertion_named = indel_notation_of(e, &hgvsc).is_some_and(|kind| names_insertion(&kind));
    let coding_bounds = |row: &Entry| position_bounds(row, "CDS_position", "cds_start", "cds_end");
    let rows = ShiftedRows {
        indel: indel_row(e)? || (insertion_named && equal_span_allele_row(e)),
        terms: [&e.consequence, &actual.consequence],
        offsets: [
            &field_of(Some(e), HGVS_OFFSET_FIELDS),
            &field_of(Some(actual), HGVS_OFFSET_FIELDS),
        ],
        coding_bounds: [coding_bounds(e), coding_bounds(actual)],
        candidate_hgvsc: &hgvsc,
    };
    Ok(start_test_read_on_the_written_span(&rows, ev, av))
}

/// `hgvsp_start_codon_test_past_the_shifted_allele` on a differing HGVSp: the two rows' terms,
/// SHIFT_LENGTH and CDS position give the two values the class's shape
/// (`start_codon_test_past_the_shifted_allele`), on an insertion's or deletion's row. The
/// SHIFT_LENGTH is the field `positive_shift` reads, as each row writes it. The two-base
/// allele at a two-base Location is read as `shifted_indel_written_span_gate` reads it, from
/// the candidate's notation. A comparison that holds no candidate row is refused.
pub fn start_codon_test_past_the_shifted_allele_rows(
    e: &Entry,
    actual: Option<&Entry>,
    ev: &str,
    av: &str,
) -> Result<bool> {
    let Some(actual) = actual else {
        return Ok(false);
    };
    let hgvsc = field_of(Some(actual), HGVSC_FIELDS);
    let insertion_named = indel_notation_of(e, &hgvsc).is_some_and(|kind| names_insertion(&kind));
    let coding_bounds = |row: &Entry| position_bounds(row, "CDS_position", "cds_start", "cds_end");
    let rows = StartTestRows {
        indel: indel_row(e)? || (insertion_named && equal_span_allele_row(e)),
        terms: [&e.consequence, &actual.consequence],
        shift_lengths: [
            first_field(e, SHIFT_LENGTH_FIELDS),
            first_field(actual, SHIFT_LENGTH_FIELDS),
        ],
        coding_bounds: [coding_bounds(e), coding_bounds(actual)],
    };
    Ok(start_codon_test_past_the_shifted_allele(&rows, ev, av))
}

/// Whether an HGVSc body places an insertion or deletion in the coding sequence by its
/// coordinates alone, whatever exon holds it: a `del` whose every coordinate is a coding
/// base; a `dup` that ends on a coding base, its first coordinate a coding base or one
/// before the coding start (`UTR_COORDINATE_WHOLE`); an `ins` between two coding bases.
fn reaches_the_coding_sequence(body: &str) -> bool {
    let Some(notation) = INDEL_NOTATION_WHOLE.captures(body) else {
        return false;
    };
    if !body.starts_with("c.") {
        return false;
    }
    let first = &notation[1];
    let last = notation.get(2).map(|end| end.as_str());
    let coding = |coordinate: &str| CODING_COORDINATE_WHOLE.is_match(coordinate);
    match &notation[3][..3] {
        "del" => coding(first) && last.is_none_or(coding),
        "dup" => {
            coding(last.unwrap_or(first))
                && (last.is_none() || coding(first) || UTR_COORDINATE_WHOLE.is_match(first))
        }
        _ => coding(first) && last.is_some_and(coding),
    }
}

/// Whether a row with no cDNA position is one of an allele written in an intron and shifted
/// `shift` bases into the coding sequence. The HGVSc body alone says so in two cases. It
/// names or bounds the first base of the exon the allele is shifted into as coding base 2 or
/// later (`exon_start_shifted`): the intron lies inside the coding region. Or it is a `dup`
/// that ends on a coding base and begins before the coding start: its first coordinate is
/// the first base of the copy that stands before the shifted insertion point, and the
/// allele as written stands at or before that base where the shift is as long as the
/// allele, or inside the copy where it is shorter, so its intron is one of the 5' UTR
/// (`utr`, BaseVariationFeatureOverlapAllele.pm:523-525) or one inside the coding region
/// that the copy spans (`non_coding`, :535-537). Where the bound is below 2, the exon spans
/// of `edges` say so (`placed_exon_start`), whatever the exon's first coding base: a row
/// they do not place, and any such row without `edges`, is not taken.
fn shifted_out_of_an_intron(
    e: &Entry,
    body: &str,
    shift: i128,
    edges: Option<&TranscriptEdges>,
) -> Result<bool> {
    match exon_start_shifted(body, shift) {
        Some(exon_start) if exon_start >= 2 => Ok(true),
        // No bound: the body is no notation in the coding sequence, or it is a `dup` that
        // begins before the coding start and so names no coding base of its exon.
        None => Ok(reaches_the_coding_sequence(body)),
        Some(_) => match edges {
            Some(edges) => Ok(placed_exon_start(e, body, shift, edges)?.is_some()),
            None => Ok(false),
        },
    }
}

/// The coding coordinate of the first base of the exon into which an HGVSc body places an
/// insertion or deletion shifted `shift` bases, counted along the exon spans of the
/// candidate cache: the coordinate `shifted_exon_start` bounds, read exactly, and a number
/// below 1 for an exon that opens before the coding start, 0 for one base before it. The
/// body is a `del` or an `ins` of coding bases alone, or a `dup` that ends on a coding
/// base, whatever its first coordinate. The base the notation starts at (the first base of
/// a `del`, the base after a `dup` or an `ins`) stands `shift` bases past the allele as the
/// row's Location writes it, in the transcript's direction, and the exon that holds it
/// opens as many coding bases before it as exon bases lie between the two.
///
/// `None` where the body is no such notation, the row has no Location to read (a text
/// row's own, or the one the default format writes for a VCF or JSON entry), the cache
/// files no transcript of the row's Feature there, the shifted base lies in no exon, or the
/// allele as written does not lie wholly in the intron before that exon.
pub fn placed_exon_start(
    e: &Entry,
    body: &str,
    shift: i128,
    edges: &TranscriptEdges,
) -> Result<Option<i128>> {
    let Some(notation) = INDEL_NOTATION_WHOLE.captures(body) else {
        return Ok(None);
    };
    let first = &notation[1];
    let last = notation.get(2).map(|end| end.as_str());
    let kind = &notation[3][..3];
    let coding = |coordinate: &str| CODING_COORDINATE_WHOLE.is_match(coordinate);
    // A `dup` starts at the base after its last coordinate, so its first one is not read.
    let on_coding_bases = if kind == "dup" {
        coding(last.unwrap_or(first))
    } else {
        coding(first) && last.is_none_or(coding)
    };
    if !body.starts_with("c.") || !on_coding_bases {
        return Ok(None);
    }
    let insertion = names_insertion(kind);
    // The coding base the notation starts at, as `exon_start_shifted` reads it.
    let named = match (kind, last) {
        ("del", _) => whole_number(first),
        ("dup", _) => whole_number(last.unwrap_or(first)).saturating_add(1),
        (_, Some(last)) => whole_number(last),
        (_, None) => return Ok(None),
    };
    let Some(location) = e.location.as_deref().or(e.allele_location.as_deref()) else {
        return Ok(None);
    };
    let (chrom, span) = location.rsplit_once(':').unwrap_or(("", location));
    let (low, high) = parted(span, '-');
    if !pyfmt::is_digit(low) || (!high.is_empty() && !pyfmt::is_digit(high)) {
        return Ok(None);
    }
    let low = integer_of(low, "Location position")?;
    let high = if high.is_empty() {
        low
    } else {
        integer_of(high, "Location position")?
    };
    // An insertion's Location names the two bases that flank it.
    if low > high || (insertion && i128::from(high) - i128::from(low) != 1) {
        return Ok(None);
    }
    let Some((exons, reverse)) = edges.transcript(chrom, low, high, &e.feature)? else {
        return Ok(None);
    };
    // The base the notation starts at: on the forward strand the deletion's lowest base or
    // the base above the insertion, on the reverse strand the deletion's highest base or
    // the base below the insertion, carried the shift along the transcript.
    let written = i128::from(match (reverse, insertion) {
        (false, false) => low,
        (false, true) | (true, false) => high,
        (true, true) => low,
    });
    let shifted = if reverse {
        written.saturating_sub(shift)
    } else {
        written.saturating_add(shift)
    };
    let Some(index) = exons
        .iter()
        .position(|(start, end)| i128::from(*start) <= shifted && shifted <= i128::from(*end))
    else {
        return Ok(None);
    };
    let (start, end) = (i128::from(exons[index].0), i128::from(exons[index].1));
    // The allele as written lies wholly in the intron before that exon.
    let in_the_intron_before = if reverse {
        exons.get(index + 1).is_some_and(|(next_start, _)| {
            end < i128::from(low) && i128::from(high) < i128::from(*next_start)
        })
    } else {
        index > 0 && i128::from(exons[index - 1].1) < i128::from(low) && i128::from(high) < start
    };
    if !in_the_intron_before {
        return Ok(None);
    }
    let along = if reverse {
        end - shifted
    } else {
        shifted - start
    };
    Ok(Some(named.saturating_sub(along)))
}

/// The coding coordinate of the first base of the exon into which an HGVSc body places an
/// insertion or deletion shifted `shift` bases out of the intron before it, or a lower bound
/// of that coordinate; `None` where the body does not place the allele in the coding
/// sequence. A `del` lies there when its every coordinate is a coding base (no UTR sign, no
/// intron offset); a `dup` when it ends on a coding base, its first coordinate a coding
/// base or a base of the intron before; an `ins` between two coding bases. A `dup` that
/// ends on an intron base and an `ins` beside one are refused, and with them every
/// insertion the shift leaves between the intron's last base and the exon's first base,
/// the copy of that one last base included, and an insertion before the stop codon's last
/// base, whose intron base carries no coordinate (`*-1`). A `dup`'s first coordinate in
/// the intron names the exon's first base, the base it counts back from or the one after
/// the coding base it counts forward from, and the stop codon's last base reads as the
/// largest one there; a notation of coding bases alone gives the bound
/// `start - shift + 1`, `start` the first base of a `del` or the base after a `dup`
/// or an `ins`. A `dup` whose first coordinate lies before the coding start names no
/// coding base of the exon and gives `None`.
pub fn shifted_exon_start(body: &str, shift: i64) -> Option<i64> {
    exon_start_shifted(body, i128::from(shift)).map(nearest_i64)
}

/// `shifted_exon_start` on whole numbers past `i64`.
fn exon_start_shifted(body: &str, shift: i128) -> Option<i128> {
    let notation = INDEL_NOTATION_WHOLE.captures(body)?;
    if !body.starts_with("c.") {
        return None;
    }
    let first = &notation[1];
    let last = notation.get(2).map(|end| end.as_str());
    let kind = &notation[3];
    let coding = |coordinate: &str| CODING_COORDINATE_WHOLE.is_match(coordinate);
    // The base `after` bases past `base`, carried back over the shift, and the one after it.
    let bound = |base: &str, after: i128| {
        whole_number(base)
            .saturating_add(after)
            .saturating_sub(shift)
            .saturating_add(1)
    };
    if kind == "del" {
        return (coding(first) && last.is_none_or(coding)).then(|| bound(first, 0));
    }
    if kind == "dup" {
        let end = last.unwrap_or(first);
        if !coding(end) {
            return None;
        }
        if last.is_none() || coding(first) {
            return Some(bound(end, 1));
        }
        return named_exon_start(first);
    }
    match last {
        Some(last) if coding(first) && coding(last) => Some(bound(last, 0)),
        _ => None,
    }
}

/// The first base of the exon an intron coordinate counts from: the base an acceptor-side
/// coordinate counts back from, the largest whole number of a machine word for the stop
/// codon's last base, or the base after the one a donor-side coordinate counts forward from.
fn named_exon_start(coordinate: &str) -> Option<i128> {
    if let Some(acceptor) = ACCEPTOR_INTRON_COORDINATE_WHOLE.captures(coordinate) {
        return Some(if &acceptor[1] == "*" {
            i128::from(i64::MAX)
        } else {
            whole_number(&acceptor[1])
        });
    }
    let donor = DONOR_INTRON_COORDINATE_WHOLE.captures(coordinate)?;
    Some(whole_number(&donor[1]).saturating_add(1))
}

/// Whether the entry's allele is two or more bases over a Location of as many, with no
/// record alleles to read: the text an insertion and a substitution of that length both
/// print, which `indel_row` reads as a substitution.
pub fn equal_span_allele_row(e: &Entry) -> bool {
    if e.indel.is_some() {
        return false;
    }
    let Some(location) = e.location.as_deref() else {
        return false;
    };
    if e.allele.len() < 2 || !e.allele.bytes().all(|base| b"ACGTN".contains(&base)) {
        return false;
    }
    let coordinates = location.split_once(':').map_or(location, |(_, rest)| rest);
    let (start, _, end) = pyfmt::partition(coordinates, "-");
    if !(pyfmt::is_digit(start) && pyfmt::is_digit(end)) {
        return false;
    }
    let end = whole_number_text(end);
    let start = whole_number_text(start);
    end == plus(&start, e.allele.len() - 1)
}

/// A key ordering HGVS cDNA coordinates along the transcript: the 5' UTR (`-N`), the
/// coding sequence (`N`), the 3' UTR (`*N`), each by its position, then an intron offset
/// (`+M` after the exon base, `-M` before it). `*-M`, an intron distance before the stop
/// codon's last base printed without that base's coordinate
/// (TranscriptVariationAllele.pm:2854-2859), orders after every coding coordinate and before
/// `*1`.
pub fn cdna_coordinate_order(coordinate: &str) -> Result<(i64, i64, i64)> {
    let (region, number) = if let Some(rest) = coordinate.strip_prefix('*') {
        (1, rest)
    } else if let Some(rest) = coordinate.strip_prefix('-') {
        (-1, rest)
    } else {
        (0, coordinate)
    };
    let Some(parts) = COORDINATE_PARTS_WHOLE.captures(number) else {
        return Ok((region, 0, 0));
    };
    let digits = |group: usize, absent: &'static str| {
        parts
            .get(group)
            .map(|found| found.as_str())
            .filter(|text| !text.is_empty())
            .map_or(absent.to_string(), str::to_string)
    };
    let position = integer_of(&digits(1, "0"), "HGVS coordinate")?;
    let offset = integer_of(&digits(2, "0"), "HGVS coordinate offset")?;
    Ok((
        region,
        if region < 0 { -position } else { position },
        offset,
    ))
}

/// `hgvs_intron_offset_abutting_insertion` on one differing HGVSIntronOffset field: the
/// reference fills both of the plugin's fields and the candidate neither, on an insertion's row
/// both engines annotate outside the transcript, its one term `upstream_gene_variant` or
/// `downstream_gene_variant` beside a DISTANCE of 0, with no HGVSc and no HGVS_OFFSET on
/// either side. `HGVSIntronOffset.pm:62` asks for the notation on every Transcript row, past
/// the `within_feature` test that keeps HGVSc off such a row (`OutputFactory.pm:1730`); an
/// insertion directly beside the transcript's first or last base has one slice coordinate
/// inside (`TranscriptVariationAllele.pm:2729`) and is clamped onto that base (:2732-2733),
/// so the offsets stored at :1543-1545 belong to a notation the row never prints. The pair is
/// 0 and 0; or 1 and 0 on a downstream row of a reverse-strand transcript, where :2802 wraps
/// the position past the last base to the far exon and :2822-2826 build "last base +1"; or 0
/// and -1 on an upstream row of a forward-strand transcript, where `_clip_alleles` strips
/// the first base off the inserted allele's end and moves the notation's end to the position
/// before it (:2231-2238), which :2802 and :2822-2826 read as "first base -1" (a row without
/// STRAND admits either). Any other pair, a row with a term inside the transcript, a value on
/// the candidate's side and the fields on the candidate's side alone stay differences.
pub fn intron_offset_abutting_insertion(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> bool {
    let Some(actual) = actual else {
        return false;
    };
    if !field_named(name, HGVS_INTRON_OFFSET_FIELDS) {
        return false;
    }
    if blank(ev) || !blank(av) || field_of(Some(e), &[name]) != ev {
        return false;
    }
    if HGVS_INTRON_OFFSET_FIELDS
        .iter()
        .any(|field| !blank(&field_of(Some(actual), &[*field])))
    {
        return false;
    }
    let term = consequence_set(&e.consequence);
    if !OUTSIDE_TRANSCRIPT_TERMS.contains(&term.as_str())
        || consequence_set(&actual.consequence) != term
    {
        return false;
    }
    for entry in [e, actual] {
        if field_of(Some(entry), &["DISTANCE"]) != "0" {
            return false;
        }
        if HGVSC_FIELDS
            .iter()
            .chain(HGVS_OFFSET_FIELDS)
            .any(|field| !blank(&field_of(Some(entry), &[*field])))
        {
            return false;
        }
    }
    if !insertion_row(e) {
        return false;
    }
    let start = field_of(Some(e), &[HGVS_INTRON_OFFSET_FIELDS[0]]);
    let end = field_of(Some(e), &[HGVS_INTRON_OFFSET_FIELDS[1]]);
    let pair = (start.as_str(), end.as_str());
    if pair == ("0", "0") {
        return true;
    }
    let strand = field_of(Some(e), STRAND_FIELDS);
    if term == "downstream_gene_variant" {
        return matches!(strand.as_str(), "-1" | ABSENT) && pair == ("1", "0");
    }
    matches!(strand.as_str(), "1" | ABSENT) && pair == ("0", "-1")
}
