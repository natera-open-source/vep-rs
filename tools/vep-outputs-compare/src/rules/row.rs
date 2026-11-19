// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rules: what they read of an output row (the row's own columns and,
//! where the first column spells the record's alleles, whether they make the row's allele an
//! insertion or a deletion, whether the record's reference allele is empty once the shared
//! first base is dropped, and the reference bases the allele stands in for) and their
//! predicates over the reference's and vep-rs's row of one pair, each written from its
//! definition in `compare_vep_outputs.py`.

use regex::Regex;
use std::collections::BTreeSet;
use std::sync::LazyLock;

/// One output row as a row-aware pair rule reads it, the `PairRow` of
/// `compare_vep_outputs.py`: the Location and the Allele as the row writes them, its
/// Feature_type, Consequence, cDNA_position and CDS_position, and the allele's shape from the
/// record's own alleles, `None` where the first column spells none: whether they make it an
/// insertion or a deletion (`indel`), and whether every sequence allele is an insertion, the
/// reference allele empty once the shared first base is dropped (`insertion`). `reference` is
/// the reference bases the allele stands in for at the row's Location, as the record's own
/// alleles give them by position (`replaced_reference`): from the first column where it
/// spells the alleles, else from the input records where the comparison is given them; empty
/// for a pure insertion's row, `None` where no record is read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PairRow {
    pub location: String,
    pub allele: String,
    pub feature_type: String,
    pub consequence: String,
    pub cdna_position: String,
    pub cds_position: String,
    pub indel: Option<bool>,
    pub insertion: Option<bool>,
    pub reference: Option<String>,
}

/// The first nine columns of an output row, Uploaded_variation to CDS_position, in the form
/// a row-aware pair rule reads.
pub fn pair_row(columns: &[&str; 9]) -> PairRow {
    let [name, location, allele, _, _, feature_type, consequence, cdna_position, cds_position] =
        *columns;
    let record_alleles = uploaded_alleles(name);
    let indel = record_alleles
        .as_ref()
        .and_then(|(reference, alternates)| allele_is_indel(reference, alternates, allele));
    let insertion = record_alleles
        .as_ref()
        .map(|(reference, alternates)| reference_allele_is_empty(reference, alternates));
    let reference = record_alleles.as_ref().and_then(|(reference, alternates)| {
        let position = uploaded_position(name)?;
        replaced_reference(position, reference, alternates, location, allele)
    });
    PairRow {
        location: location.to_string(),
        allele: allele.to_string(),
        feature_type: feature_type.to_string(),
        consequence: consequence.to_string(),
        cdna_position: cdna_position.to_string(),
        cds_position: cds_position.to_string(),
        indel,
        insertion,
        reference,
    }
}

/// Ensembl VEP's Uploaded_variation for a VCF record without an ID (OutputFactory.pm:883-931):
/// the chromosome, the position and the record's own alleles joined by `/`. The classes name
/// white space as Python's `\s` does.
static UPLOADED_ALLELES: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^[^\s\x1c-\x1f]+_[0-9]+_([ACGTNacgtn*-]+(?:/[^/\s\x1c-\x1f]+)+)$").expect("regex")
});

/// The reference and alternate alleles an identifier-less Uploaded_variation spells, or
/// `None` for an identifier, which names no allele.
pub fn uploaded_alleles(uploaded: &str) -> Option<(&str, Vec<&str>)> {
    let spelt = UPLOADED_ALLELES.captures(uploaded)?.get(1)?.as_str();
    let mut parts = spelt.split('/');
    let reference = parts.next()?;
    Some((reference, parts.collect()))
}

/// The same name, the group its position.
static UPLOADED_POSITION: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^[^\s\x1c-\x1f]+_([0-9]+)_[ACGTNacgtn*-]+(?:/[^/\s\x1c-\x1f]+)+$").expect("regex")
});

/// The POS an identifier-less Uploaded_variation spells before its alleles, or `None` for an
/// identifier. A position past `i128`, or in digits other than `0` to `9`, is read as no
/// number.
pub fn uploaded_position(uploaded: &str) -> Option<i128> {
    let position = UPLOADED_POSITION.captures(uploaded)?.get(1)?;
    position.as_str().parse().ok()
}

/// The first and the last position a Location names where both are written in digits: the
/// text after the first `:`, parted at its first `-`, one position standing for both. A
/// position past `i128` is read as no number.
fn location_ends(location: &str) -> Option<(i128, i128)> {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let (start, end) = span.split_once('-').unwrap_or((span, ""));
    let end = if end.is_empty() { start } else { end };
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    if !(digits(start) && digits(end)) {
        return None;
    }
    Some((start.parse().ok()?, end.parse().ok()?))
}

/// Whether `spelt` is `head`, `written` and `tail` one after the other.
fn spells(spelt: &[char], head: &[char], written: &[char], tail: &[char]) -> bool {
    spelt.len() == head.len() + written.len() + tail.len()
        && spelt.starts_with(head)
        && spelt[head.len()..].starts_with(written)
        && spelt.ends_with(tail)
}

/// The reference bases an output row's allele stands in for, read by position from a record's
/// own alleles, `replaced_reference` of `compare_vep_outputs.py`. The record whose reference
/// allele `reference` starts at `position` writes `allele` at `location` where one of its
/// alternates is the reference allele with the Location's bases replaced by the allele (`-`
/// replaces them by nothing): those bases are returned. Where an alternate is the reference
/// allele with the allele inserted between the Location's two neighbouring positions the row
/// is a pure insertion's and the empty string is returned. `None` where the record writes no
/// such row, where two alternates read it differently, and for an allele or a Location that
/// is not a sequence's. The record's alleles are read in upper case.
pub fn replaced_reference(
    position: i128,
    reference: &str,
    alternates: &[&str],
    location: &str,
    allele: &str,
) -> Option<String> {
    if allele != "-" && !spelt_in(allele, b"ACGTN") {
        return None;
    }
    let (start, end) = location_ends(location)?;
    let (first, last) = (start - position, end - position);
    let written: Vec<char> = match allele {
        "-" => Vec::new(),
        bases => bases.chars().collect(),
    };
    let reference: Vec<char> = reference.to_uppercase().chars().collect();
    let length = reference.len() as i128;
    let mut found: Vec<String> = Vec::new();
    let mut find = |bases: String| {
        if !found.contains(&bases) {
            found.push(bases);
        }
    };
    for alternate in alternates {
        let spelt: Vec<char> = alternate.to_uppercase().chars().collect();
        if 0 <= first && first <= last && last < length {
            let (first, last) = (first as usize, last as usize);
            if spells(
                &spelt,
                &reference[..first],
                &written,
                &reference[last + 1..],
            ) {
                find(reference[first..=last].iter().collect());
            }
        }
        if !written.is_empty() && last - first == 1 && 0 <= first && first < length {
            let after = first as usize + 1;
            if spells(&spelt, &reference[..after], &written, &reference[after..]) {
                find(String::new());
            }
        }
    }
    match found.len() {
        1 => found.pop(),
        _ => None,
    }
}

fn stripped_len(allele: &str) -> usize {
    allele.trim_matches('-').chars().count()
}

fn first_char(allele: &str) -> Option<char> {
    allele.chars().next()
}

fn without_first_char(allele: &str) -> &str {
    let mut chars = allele.chars();
    chars.next();
    chars.as_str()
}

/// Whether `allele` is an insertion or a deletion of the record whose reference allele is
/// `reference` and alternate alleles `alternates`, as `allele_is_indel` of
/// `compare_vep_outputs.py` reads a text row: the alternate the allele names and the reference
/// differ in length. A record with one alternate not spelt with a `*` has that alternate's
/// shape on every row; with several, the shared first base is dropped where Ensembl VEP drops
/// it (`Parser/VCF.pm`:300-323). `None` when the allele is not a sequence, when no alternate
/// names it, or when two alternates of different shapes name it.
pub fn allele_is_indel(reference: &str, alternates: &[&str], allele: &str) -> Option<bool> {
    if allele.is_empty() || !allele.bytes().all(|b| b"ACGTN-".contains(&b)) {
        return None;
    }
    let trimmed: Vec<&str> = alternates
        .iter()
        .copied()
        .filter(|a| !a.contains('*'))
        .collect();
    if let [only] = trimmed[..] {
        return Some(stripped_len(reference) != stripped_len(only));
    }
    let reference_len = reference.chars().count();
    let holds_indel = alternates
        .iter()
        .any(|a| a.chars().count() != reference_len);
    let trim = !reference.is_empty()
        && holds_indel
        && trimmed
            .iter()
            .all(|a| first_char(a) == first_char(reference));
    let kept_reference = if trim {
        without_first_char(reference)
    } else {
        reference
    };
    let mut shapes: Vec<bool> = Vec::new();
    for alternate in alternates {
        let alt = if trim && !alternate.contains('*') {
            without_first_char(alternate)
        } else {
            alternate
        };
        if (if alt.is_empty() { "-" } else { alt }) == allele {
            let shape = stripped_len(kept_reference) != stripped_len(alt);
            if !shapes.contains(&shape) {
                shapes.push(shape);
            }
        }
    }
    match shapes[..] {
        [shape] => Some(shape),
        _ => None,
    }
}

/// Whether the reference allele of a record is empty once the record's shared first base is
/// dropped where Ensembl VEP drops it (`Parser/VCF.pm`:300-323), `reference_allele_is_empty` of
/// `compare_vep_outputs.py` on a text row's alleles. Every sequence allele of such a record is
/// an insertion.
pub fn reference_allele_is_empty(reference: &str, alternates: &[&str]) -> bool {
    let reference_len = reference.chars().count();
    let holds_indel = alternates
        .iter()
        .any(|a| a.chars().count() != reference_len);
    let trim = !reference.is_empty()
        && holds_indel
        && alternates
            .iter()
            .filter(|a| !a.contains('*'))
            .all(|a| first_char(a) == first_char(reference));
    stripped_len(if trim {
        without_first_char(reference)
    } else {
        reference
    }) == 0
}

const CODING_SEQUENCE_TERM: &str = "coding_sequence_variant";
const STOP_RETAINED_TERM: &str = "stop_retained_variant";
const SPLICE_REGION_TERM: &str = "splice_region_variant";

/// The terms of a row's Consequence as the predicates of `compare_vep_outputs.py` read them:
/// parted at `&` and `,`, the empty ones dropped, an empty Consequence reading as one empty
/// term.
fn terms(consequence: &str) -> BTreeSet<&str> {
    let terms: BTreeSet<&str> = consequence
        .split(['&', ','])
        .filter(|term| !term.is_empty())
        .collect();
    if terms.is_empty() {
        BTreeSet::from([""])
    } else {
        terms
    }
}

/// The terms the reference's row holds, those only it holds and those only vep-rs's holds.
fn term_differences<'a>(
    e: &'a PairRow,
    a: &'a PairRow,
) -> (BTreeSet<&'a str>, BTreeSet<&'a str>, BTreeSet<&'a str>) {
    let (expected, actual) = (terms(&e.consequence), terms(&a.consequence));
    let gone = expected.difference(&actual).copied().collect();
    let come = actual.difference(&expected).copied().collect();
    (expected, gone, come)
}

fn is(terms: &BTreeSet<&str>, members: &[&str]) -> bool {
    terms.len() == members.len() && members.iter().all(|member| terms.contains(member))
}

/// Whether the row is a Transcript row: its Feature_type, read after a last `|` and up to a
/// `#`, in any letter case.
fn transcript_row(e: &PairRow) -> bool {
    let kind = e.feature_type.rsplit('|').next().unwrap_or("");
    let kind = kind.split('#').next().unwrap_or("");
    kind.to_lowercase() == "transcript"
}

fn spelt_in(allele: &str, letters: &[u8]) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| letters.contains(&b))
}

/// The number of bases a Location spans where its start and end are both written in digits:
/// the text after the first `:`, parted at its first `-`. A coordinate past `i128` is read as
/// no number.
fn location_span(location: &str) -> Option<i128> {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let (start, end) = span.split_once('-')?;
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    if !(digits(start) && digits(end)) {
        return None;
    }
    let (start, end): (i128, i128) = (start.parse().ok()?, end.parse().ok()?);
    Some((end - start).saturating_add(1))
}

/// Whether a row is an insertion's or a deletion's, `indel_row` of `compare_vep_outputs.py`:
/// an allele of `-` is one; else the record's own alleles decide where the first column
/// spells them; else the Location decides, an indel's spanning another number of bases than
/// its allele of `A`, `C`, `G`, `T` and `N` has.
fn indel_row(e: &PairRow) -> bool {
    if e.allele == "-" {
        return true;
    }
    if let Some(indel) = e.indel {
        return indel;
    }
    spelt_in(&e.allele, b"ACGTN")
        && location_span(&e.location).is_some_and(|span| span != e.allele.len() as i128)
}

/// Whether a position column carries its start and its end: `?` stands for an end that maps
/// to no transcript base, one number for both, `-` for neither, each read up to a `/total`
/// suffix (OutputFactory.pm:1726-1739).
fn position_bounds(position: &str) -> (bool, bool) {
    let position = position.split_once('/').map_or(position, |(head, _)| head);
    if matches!(position, "" | "-" | "<absent>") {
        return (false, false);
    }
    match position.split_once('-') {
        Some((start, end)) => (start != "?", end != "?"),
        None => (position != "?", position != "?"),
    }
}

/// Whether the row carries the coordinates the stop codon test edits at
/// (Utils/VariationEffect.pm:1507-1508): both cDNA positions and the CDS start.
fn stop_codon_test_coordinates(e: &PairRow) -> bool {
    position_bounds(&e.cdna_position) == (true, true) && position_bounds(&e.cds_position).0
}

/// Whether the row shows that `_ins_del_stop_altered_cil` (Utils/VariationEffect.pm:1488-1539)
/// returns before it edits anything: the allele holds a letter other than `A`, `C`, `G` and
/// `T` (:1498), or the span lacks one of the test's coordinates (:1508).
fn stop_codon_test_not_made(e: &PairRow) -> bool {
    (e.allele != "-" && !spelt_in(&e.allele, b"ACGT")) || !stop_codon_test_coordinates(e)
}

/// `stop_codon_overlap_genomic_span`: a Transcript row of an insertion or deletion with both
/// cDNA positions and a CDS start on both sides, `splice_region_variant` in the reference's
/// terms, and the two sets equal but for the reference's `coding_sequence_variant` against
/// vep-rs's `stop_retained_variant` (Utils/VariationEffect.pm:1389-1431).
pub fn stop_codon_overlap_genomic_span_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !indel_row(e) {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    is(&gone, &[CODING_SEQUENCE_TERM])
        && is(&come, &[STOP_RETAINED_TERM])
        && expected.contains(SPLICE_REGION_TERM)
        && !stop_codon_test_not_made(e)
        && !stop_codon_test_not_made(a)
}

const FRAMESHIFT_TERM: &str = "frameshift_variant";

/// Whether a row is an insertion's or a deletion's as `indel_row` reads one, an allele of
/// other letters than `A`, `C`, `G`, `T` and `N` included: where the record's own alleles say
/// nothing, such an allele against a Location of another length.
fn length_changing_row(e: &PairRow) -> bool {
    if indel_row(e) {
        return true;
    }
    e.indel.is_none()
        && !e.allele.is_empty()
        && e.allele.bytes().all(|b| b.is_ascii_alphabetic())
        && location_span(&e.location).is_some_and(|span| span != e.allele.len() as i128)
}

/// `stop_retained_without_stop_codon_test`: a Transcript row of an insertion or deletion,
/// the reference's terms holding `stop_retained_variant` and vep-rs's being the same without
/// it, with nothing, `coding_sequence_variant` or `frameshift_variant` in its place, and both
/// rows showing that the stop codon test is not made (Utils/VariationEffect.pm:1318,
/// 1488-1539).
pub fn stop_retained_without_stop_codon_test_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !length_changing_row(e) {
        return false;
    }
    let (_, gone, come) = term_differences(e, a);
    is(&gone, &[STOP_RETAINED_TERM])
        && (come.is_empty() || is(&come, &[CODING_SEQUENCE_TERM]) || is(&come, &[FRAMESHIFT_TERM]))
        && stop_codon_test_not_made(e)
        && stop_codon_test_not_made(a)
}

const STOP_GAINED_TERM: &str = "stop_gained";
/// The terms `frameshift` (Utils/VariationEffect.pm:1541-1567) reads away where
/// `stop_retained` does not hold.
const FRAMESHIFT_READ_TERMS: [&str; 3] = [
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
];

/// `stop_retained_three_base_stop_codon_test`: a Transcript row of an insertion or deletion
/// whose stop codon test is made on both sides, the reference's terms holding
/// `stop_retained_variant` and vep-rs's `coding_sequence_variant` in its place, or, where the
/// reference's hold it alone or beside one of `inframe_insertion`, `inframe_deletion` and
/// `protein_altering_variant`, `frameshift_variant` with or without `stop_gained` in the place
/// of those terms (Utils/VariationEffect.pm:1488-1539).
pub fn stop_retained_three_base_stop_codon_test_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !length_changing_row(e) {
        return false;
    }
    let (_, gone, come) = term_differences(e, a);
    let terms_stand = if is(&come, &[CODING_SEQUENCE_TERM]) {
        is(&gone, &[STOP_RETAINED_TERM])
    } else if is(&come, &[FRAMESHIFT_TERM]) || is(&come, &[STOP_GAINED_TERM, FRAMESHIFT_TERM]) {
        gone.contains(STOP_RETAINED_TERM)
            && gone.len() <= 2
            && gone
                .iter()
                .all(|term| *term == STOP_RETAINED_TERM || FRAMESHIFT_READ_TERMS.contains(term))
    } else {
        false
    };
    terms_stand && !stop_codon_test_not_made(e) && !stop_codon_test_not_made(a)
}

const STOP_LOST_TERM: &str = "stop_lost";

/// `stop_codon_removed_without_stop_codon_test`: a Transcript row of an insertion or
/// deletion, both rows lacking a cDNA position at one end or a CDS start, and the two sets of
/// terms equal but for the reference's `stop_retained_variant` against vep-rs's `stop_lost`
/// (Utils/VariationEffect.pm:1234-1288, 1508).
pub fn stop_codon_removed_without_stop_codon_test_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !length_changing_row(e) {
        return false;
    }
    let (_, gone, come) = term_differences(e, a);
    is(&gone, &[STOP_RETAINED_TERM])
        && is(&come, &[STOP_LOST_TERM])
        && !stop_codon_test_coordinates(e)
        && !stop_codon_test_coordinates(a)
}

const INFRAME_INSERTION_TERM: &str = "inframe_insertion";

/// Whether a row's allele is a sequence of `A`, `C`, `G`, `T` and `N` and the record's own
/// alleles mark it an indel where the first column spells them, else the Location spans the
/// two bases around an insertion's point.
fn insertion_row(e: &PairRow) -> bool {
    if !spelt_in(&e.allele, b"ACGTN") {
        return false;
    }
    match e.indel {
        Some(indel) => indel,
        None => location_span(&e.location) == Some(2),
    }
}

/// Whether a row `insertion_row` takes is an insertion between two neighbouring positions,
/// `inserted_between_neighbours` of `compare_vep_outputs.py`: its Location is two neighbouring
/// positions, and the record's reference allele is empty once the shared first base is dropped
/// where the first column spells the alleles; under an identifier its allele has another length than two instead, so two bases
/// replaced by a sequence of another length are taken with the insertions.
fn inserted_between_neighbours(e: &PairRow) -> bool {
    location_span(&e.location) == Some(2) && e.insertion.unwrap_or(e.allele.len() != 2)
}

/// `stop_codon_overlap_insertion_range`: a Transcript row of an insertion between two
/// neighbouring positions (`inserted_between_neighbours`) with both cDNA
/// positions and a CDS start on both sides, the two sets of terms equal but for the
/// reference's `coding_sequence_variant` against vep-rs's `stop_retained_variant` away from
/// `splice_region_variant`, or the reference's `frameshift_variant` with or without
/// `stop_gained` against vep-rs's `stop_retained_variant` with or without `inframe_insertion`
/// (Utils/VariationEffect.pm:1389-1431).
pub fn stop_codon_overlap_insertion_range_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !insertion_row(e) || !inserted_between_neighbours(e) {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    let terms_stand = if is(&gone, &[CODING_SEQUENCE_TERM]) {
        is(&come, &[STOP_RETAINED_TERM]) && !expected.contains(SPLICE_REGION_TERM)
    } else if is(&gone, &[FRAMESHIFT_TERM]) || is(&gone, &[STOP_GAINED_TERM, FRAMESHIFT_TERM]) {
        is(&come, &[STOP_RETAINED_TERM]) || is(&come, &[INFRAME_INSERTION_TERM, STOP_RETAINED_TERM])
    } else {
        false
    };
    terms_stand && stop_codon_test_coordinates(e) && stop_codon_test_coordinates(a)
}

/// A set of the seven terms a row of `intron_rows` holds, one bit a term in the order of
/// [`INTRON_ROW_TERMS`].
type TermSet = u8;

/// The five splice terms whose windows lie inside an intron: the two bases at its 5' end, the
/// two at its 3' end, its fifth base, its third to sixth base, and the bases 3 to 17 from its
/// 3' end (Utils/VariationEffect.pm:558-614); the intron term, whose flag `within_intron`
/// (Utils/VariationEffect.pm:629-634) returns and whose window is the intron itself, its end
/// bases included; and `splice_region_variant`.
const INTRON_ROW_TERMS: [&str; 7] = [
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "intron_variant",
    SPLICE_REGION_TERM,
];
const DONOR: TermSet = 1;
const ACCEPTOR: TermSet = 1 << 1;
const DONOR_5TH_BASE: TermSet = 1 << 2;
const DONOR_REGION: TermSet = 1 << 3;
const TRACT: TermSet = 1 << 4;
const INTRON: TermSet = 1 << 5;
const REGION: TermSet = 1 << 6;
/// The six terms `_intron_effects` sets from the differing regions and vep-rs writes only
/// where the record's span holds a base of the term's window.
const PAST_THE_SPAN_TERMS: TermSet =
    DONOR | ACCEPTOR | DONOR_5TH_BASE | DONOR_REGION | TRACT | INTRON;
/// The four beside which `splice_region` returns 0 (Utils/VariationEffect.pm:616-627).
const SPLICE_REGION_WITHHOLDING_TERMS: TermSet = DONOR | ACCEPTOR | DONOR_5TH_BASE | DONOR_REGION;
/// A gap past which no window of an intron's end reaches a span or its regions.
const FAR: i64 = 1_000_000_000;
/// The row-aware rules of the stop codon classes, whose terms a row can differ by beside the
/// six: one row can carry both differences.
const STOP_CODON_ROW_PAIRS: [fn(&PairRow, &PairRow) -> bool; 5] = [
    stop_codon_overlap_genomic_span_pair,
    stop_retained_without_stop_codon_test_pair,
    stop_retained_three_base_stop_codon_test_pair,
    stop_codon_removed_without_stop_codon_test_pair,
    stop_codon_overlap_insertion_range_pair,
];

/// The terms of the seven among a row's terms.
fn intron_row_terms(terms: &BTreeSet<&str>) -> TermSet {
    INTRON_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, term)| terms.contains(*term))
        .fold(0, |set, (bit, _)| set | 1 << bit)
}

/// The terms of a row outside the seven.
fn other_terms<'a>(terms: &BTreeSet<&'a str>) -> BTreeSet<&'a str> {
    terms
        .iter()
        .copied()
        .filter(|term| !INTRON_ROW_TERMS.contains(term))
        .collect()
}

/// The start and the end a position column writes: an empty cell or `-` is `-` twice; any
/// other cell is read up to a `/total` suffix and parted at its first `-`, one number standing
/// for both ends. `None` where a part is not digits, `?` or `-`.
fn position_parts(position: &str) -> Option<(&str, &str)> {
    if matches!(position, "" | "-") {
        return Some(("-", "-"));
    }
    let head = position.split_once('/').map_or(position, |(head, _)| head);
    let (start, end) = head.split_once('-').unwrap_or((head, head));
    let written = |part: &str| {
        matches!(part, "?" | "-") || (!part.is_empty() && part.bytes().all(|b| b.is_ascii_digit()))
    };
    (written(start) && written(end)).then_some((start, end))
}

/// The cDNA start and end and the CDS start and end a row carries, as written (`?` for an end
/// that maps to no transcript base, `-` for a position the row lacks), `transcript_positions`
/// of `compare_vep_outputs.py` on a text row; `None` where a position is not of that form.
fn transcript_positions(e: &PairRow) -> Option<[&str; 4]> {
    let (cdna_start, cdna_end) = position_parts(&e.cdna_position)?;
    let (cds_start, cds_end) = position_parts(&e.cds_position)?;
    Some([cdna_start, cdna_end, cds_start, cds_end])
}

/// The number of positions a Location names, `location_span` of `compare_vep_outputs.py`, or
/// `None` where it names no span of whole positions.
fn positions_named(location: &str) -> Option<i128> {
    let (start, end) = location_ends(location)?;
    let span = (end - start).saturating_add(1);
    (span >= 1).then_some(span)
}

/// The number of bases of the transcript's sequence from the first cDNA position of a row to
/// its second, both counted, where both are numbers. A position past `i128` is read as no
/// number.
fn sequence_width(positions: &[&str; 4]) -> Option<i128> {
    let (start, end): (i128, i128) = (positions[0].parse().ok()?, positions[1].parse().ok()?);
    Some((end - start).saturating_add(1))
}

/// Whether a row shows, without its record, a span of three or more bases that lies in the
/// transcript's sequence from end to end under an allele of more nucleotide letters than the
/// span has bases, `longer_allele_inside_one_exon` of `compare_vep_outputs.py`: its Location
/// names the span, and its two cDNA positions are numbers as far apart as the span's two
/// ends, which a span that holds an intron's base does not have. A Location of two
/// neighbouring positions is refused: an insertion between them is written with the same
/// Location, allele and positions as two replaced bases, and only the record tells the two
/// apart.
fn longer_allele_inside_one_exon(e: &PairRow) -> bool {
    let Some(span) = positions_named(&e.location) else {
        return false;
    };
    if span < 3 || !spelt_in(&e.allele, b"ACGTN") || e.allele.len() as i128 <= span {
        return false;
    }
    transcript_positions(e).is_some_and(|positions| sequence_width(&positions) == Some(span))
}

/// The regions `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460) cuts a
/// record into, as offsets from the span's first base: where both alleles hold more than one
/// base (:418) the runs of indices at which they differ, every index past the shorter allele's
/// end differing (:428-450); else one region, the span (:452-454).
fn differing_regions(reference: &[char], allele: &[char]) -> Vec<(i64, i64)> {
    let span = vec![(0, reference.len() as i64 - 1)];
    if reference.len() <= 1 || allele.len() <= 1 {
        return span;
    }
    let mut runs: Vec<(i64, i64)> = Vec::new();
    for index in 0..reference.len().max(allele.len()) {
        if index < reference.len() && index < allele.len() && reference[index] == allele[index] {
            continue;
        }
        let index = index as i64;
        match runs.last_mut() {
            Some(run) if run.1 == index - 1 => run.1 = index,
            _ => runs.push((index, index)),
        }
    }
    if runs.is_empty() {
        span
    } else {
        runs
    }
}

/// The splice and intron terms of the reference's row and of vep-rs's for a record whose span
/// is the offsets `0` to `span - 1`, on a transcript of the given strand whose one intron near
/// the span runs from offset `start` to offset `end`, `intron_rows` of
/// `compare_vep_outputs.py`.
///
/// The reference's row is `_intron_effects` (BaseTranscriptVariationAllele.pm:99-224) on
/// `regions`, each placed with no bound at the span's end (:117-118). The intron is read
/// where the span overlaps it widened by three bases, and its two ends where the span
/// overlaps the eleven bases around each (BaseTranscriptVariation.pm:1027-1029); an intron of
/// thirteen bases or fewer is a frameshift intron (:1032-1035) and a region that overlaps it
/// sets nothing (BaseTranscriptVariationAllele.pm:138-141, :171-174). A region in the
/// intron's interior sets the intron term (:143-149) and one in a window of the five splice
/// terms its flag (:152-159, :176-204). The flag of `splice_region` is assigned for each
/// region in turn, whether the region lies in one of its four windows, until a region has
/// met the intron's first two or last two bases (:209-216; Utils/VariationEffect.pm:91-94),
/// so the last region assigned decides it. The tract term is written only off the exons
/// (Constants.pm:971-983), which a transcript with a frameshift intron reads twelve bases
/// wider (BaseTranscriptVariation.pm:861; `wide_exons` says another intron of the transcript
/// is one). The donor region is withheld beside the fifth base and `splice_region_variant`
/// beside the first four terms (Utils/VariationEffect.pm:597-627).
///
/// vep-rs's row keeps a term of the six only where the span itself holds a base of a window
/// of that term, the intron being the intron term's window, reads `splice_region` from any
/// region in one of its four windows, and withholds the donor region and
/// `splice_region_variant` by its own row's terms.
fn intron_rows(
    regions: &[(i64, i64)],
    span: i64,
    start: i64,
    end: i64,
    forward: bool,
    wide_exons: bool,
) -> (TermSet, TermSet) {
    let last = span - 1;
    let in_list = last >= start - 3 && end + 3 >= 0;
    let at_an_end = (last >= start - 3 && start + 7 >= 0) || (last >= end - 7 && end + 3 >= 0);
    let frameshift = end - start <= 12;
    let read = || {
        let outside = move |&(from, to): &(i64, i64)| !(frameshift && to >= start && from <= end);
        regions.iter().copied().filter(outside)
    };
    let met =
        |first: i64, final_base: i64| read().any(|(from, to)| to >= first && from <= final_base);
    let held = |first: i64, final_base: i64| last >= first && final_base >= 0;
    let windows = if forward {
        [
            (DONOR, start, start + 1),
            (ACCEPTOR, end - 1, end),
            (DONOR_5TH_BASE, start + 4, start + 4),
            (DONOR_REGION, start + 2, start + 5),
            (TRACT, end - 16, end - 2),
        ]
    } else {
        [
            (DONOR, end - 1, end),
            (ACCEPTOR, start, start + 1),
            (DONOR_5TH_BASE, end - 4, end - 4),
            (DONOR_REGION, end - 5, end - 2),
            (TRACT, start + 2, start + 16),
        ]
    };
    let off_the_exons = !(start > 0
        || last > end
        || ((frameshift || wide_exons) && (start > -12 || last + 12 > end)));
    let (mut set_by_a_region, mut held_by_the_span) = (0, 0);
    for (term, first, final_base) in windows {
        let read_here = if term == TRACT {
            (in_list || at_an_end) && off_the_exons
        } else {
            at_an_end
        };
        if read_here && met(first, final_base) {
            set_by_a_region |= term;
        }
        if held(first, final_base) {
            held_by_the_span |= term;
        }
    }
    if in_list && start + 2 <= end - 2 && met(start + 2, end - 2) {
        set_by_a_region |= INTRON;
    }
    if held(start, end) {
        held_by_the_span |= INTRON;
    }
    let beside = [
        (start - 3, start - 1),
        (start + 2, start + 7),
        (end - 7, end - 2),
        (end + 1, end + 3),
    ];
    let (mut last_assigned, mut in_any, mut at_a_site) = (false, false, false);
    if at_an_end {
        for (from, to) in read() {
            let in_one = beside
                .iter()
                .any(|&(first, final_base)| to >= first && from <= final_base);
            in_any |= in_one;
            at_a_site |= (to >= start && from <= start + 1) || (to >= end - 1 && from <= end);
            if !at_a_site {
                last_assigned = in_one;
            }
        }
    }
    let row = |mut terms: TermSet, region: bool| {
        if terms & DONOR_5TH_BASE != 0 {
            terms &= !DONOR_REGION;
        }
        if region && terms & SPLICE_REGION_WITHHOLDING_TERMS == 0 {
            terms |= REGION;
        }
        terms
    };
    (
        row(set_by_a_region, last_assigned),
        row(set_by_a_region & held_by_the_span, in_any),
    )
}

/// Calls `visit` with the places `(forward, start, end)` of one intron beside a span of `span`
/// bases that a row's cDNA positions allow, as offsets from the span's first base, for regions
/// that run `reach` bases past the span's end, until it returns true, and returns whether it
/// did: the places of `intron_placements` of `compare_vep_outputs.py`, in its order. A
/// numbered cDNA position is a base of the transcript's sequence and any other a base of the
/// intron; the first position is the span's first base on the forward strand and its last on
/// the reverse. Two numbers as far apart as the span put the intron after the span, one to
/// seven bases on, as far as one of its ends is still read, or before it, where an intron
/// shorter than its windows has them run on past its last base; two numbers closer than the
/// span put it inside the span. A gap wider than any window, and than the regions' reach, is
/// placed once.
fn intron_placements(
    span: i64,
    reach: i64,
    positions: &[&str; 4],
    mut visit: impl FnMut(bool, i64, i64) -> bool,
) -> bool {
    let last = span - 1;
    let before: Vec<i64> = (0..18).map(|gap| -gap).chain([-FAR]).collect();
    let gaps: BTreeSet<i64> = (0..42).chain(42.max(reach - 20)..reach + 19).collect();
    let after: Vec<i64> = gaps.iter().map(|gap| last + gap).chain([FAR]).collect();
    let numbered = |part: &str| !part.is_empty() && part.bytes().all(|b| b.is_ascii_digit());
    let (numbered_first, numbered_last) = (numbered(positions[0]), numbered(positions[1]));
    for forward in [true, false] {
        let (first_in_sequence, last_in_sequence) = match forward {
            true => (numbered_first, numbered_last),
            false => (numbered_last, numbered_first),
        };
        if first_in_sequence && last_in_sequence {
            let Some(width) = sequence_width(positions) else {
                continue;
            };
            if width == i128::from(span) {
                for start in (1..=7).map(|gap| last + gap) {
                    let far_ends = after.iter().copied().filter(|end| *end >= start);
                    let near_ends = (0..42).map(|length| start + length);
                    if far_ends
                        .chain(near_ends)
                        .any(|end| visit(forward, start, end))
                    {
                        return true;
                    }
                }
                if (-16..0).any(|start| (start..0).any(|end| visit(forward, start, end))) {
                    return true;
                }
            } else if 0 < width && width < i128::from(span) {
                let length = span - width as i64;
                if (1..=last - length).any(|start| visit(forward, start, start + length - 1)) {
                    return true;
                }
            }
        } else if first_in_sequence {
            if (1..span).any(|start| after.iter().any(|end| visit(forward, start, *end))) {
                return true;
            }
        } else if last_in_sequence {
            if (0..last).any(|end| before.iter().any(|start| visit(forward, *start, end))) {
                return true;
            }
        } else if before
            .iter()
            .any(|start| after.iter().any(|end| visit(forward, *start, *end)))
        {
            return true;
        }
    }
    false
}

/// A place of the transcript's intron, `(forward, start, end, wide_exons)` as `intron_rows`
/// reads it, at which the record's regions give the reference row's splice and intron terms
/// `expected` and the span gives vep-rs's `actual`, or `None` where the cDNA positions
/// `positions` allow no such place, `record_span_placement` of `compare_vep_outputs.py`.
/// `reference` is the reference bases the allele stands in for; the places are those of
/// `intron_placements`.
fn record_span_placement(
    positions: &[&str; 4],
    reference: &[char],
    allele: &[char],
    expected: TermSet,
    actual: TermSet,
) -> Option<(bool, i64, i64, bool)> {
    let regions = differing_regions(reference, allele);
    let span = reference.len() as i64;
    let rows = (expected, actual);
    let mut place = None;
    intron_placements(
        span,
        allele.len() as i64 - span,
        positions,
        |forward, start, end| {
            let narrow = intron_rows(&regions, span, start, end, forward, false);
            let wide_exons = if narrow == rows {
                false
            } else if narrow.0 & TRACT != 0
                && expected & TRACT == 0
                && intron_rows(&regions, span, start, end, forward, true) == rows
            {
                true
            } else {
                return false;
            };
            place = Some((forward, start, end, wide_exons));
            true
        },
    );
    place
}

/// `splice_term_past_the_record_span`: a Transcript row with the same cDNA and CDS positions
/// on both sides whose two sets of terms differ by terms of the six the reference alone holds
/// and, at most, the donor region and `splice_region_variant` vep-rs alone holds, beside the
/// terms of a stop codon class where its rule takes the rest of the difference
/// (`STOP_CODON_ROW_PAIRS`).
///
/// `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460) compares the two
/// alleles index by index to the longer allele's length (:428-429), so every index past the
/// reference allele's end differs, and `_intron_effects` places each differing region on the
/// genome from the record's start with no bound at the span's end
/// (BaseTranscriptVariationAllele.pm:117-118): on a record whose alternate allele is the
/// longer one, the flags of the five splice terms and of the intron term (:143-204) are set
/// from positions the record does not replace. vep-rs writes such a term only where a base of
/// the span itself lies in a window of that term, the intron being the intron term's window.
///
/// With the record's alleles (`PairRow::reference`) the two rows are those of a longer
/// alternate allele that is no pure insertion, both alleles holding more than one base, at a
/// place of one intron the cDNA positions allow (`record_span_placement`): no input names the
/// transcript's introns, so the place is read off the reference's row, and where two places
/// give that row and differ in what the span holds, the pair is taken at the one that gives
/// vep-rs's row. Without them the span of three or more positions lies in the transcript's
/// sequence under a longer allele (`longer_allele_inside_one_exon`), vep-rs's row holds none
/// of the six, and `splice_region_variant` alone comes, where a term that withholds it went;
/// a Location of two neighbouring positions, which an insertion's row shares, and a span that
/// holds an intron's base then stay differences.
pub fn splice_term_past_the_record_span_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) {
        return false;
    }
    let Some(positions) = transcript_positions(e) else {
        return false;
    };
    if transcript_positions(a) != Some(positions) {
        return false;
    }
    let (expected_terms, actual_terms) = (terms(&e.consequence), terms(&a.consequence));
    let expected = intron_row_terms(&expected_terms);
    let actual = intron_row_terms(&actual_terms);
    let other = other_terms(&expected_terms);
    if other != other_terms(&actual_terms) {
        // The rest of the difference is a stop codon class's: its rule reads the two rows
        // with vep-rs's splice and intron terms on both, and this one the two rows with the
        // reference's other terms on both. Those differ by terms of the seven alone, the
        // only terms of either row read below.
        let held = |bit: &usize| actual & (1 << bit) != 0;
        let splice_and_intron = (0..INTRON_ROW_TERMS.len())
            .filter(held)
            .map(|bit| INTRON_ROW_TERMS[bit]);
        let beside = PairRow {
            consequence: other
                .iter()
                .copied()
                .chain(splice_and_intron)
                .collect::<Vec<_>>()
                .join(","),
            ..e.clone()
        };
        if !STOP_CODON_ROW_PAIRS.iter().any(|pair| pair(&beside, a)) {
            return false;
        }
    }
    let (gone, come) = (expected & !actual, actual & !expected);
    if gone == 0 || gone & !PAST_THE_SPAN_TERMS != 0 || come & !(REGION | DONOR_REGION) != 0 {
        return false;
    }
    let Some(reference) = e.reference.as_ref().or(a.reference.as_ref()) else {
        return longer_allele_inside_one_exon(e)
            && longer_allele_inside_one_exon(a)
            && actual & PAST_THE_SPAN_TERMS == 0
            && (come == 0 || (come == REGION && gone & SPLICE_REGION_WITHHOLDING_TERMS != 0));
    };
    let reference: Vec<char> = reference.chars().collect();
    let allele: Vec<char> = e.allele.chars().collect();
    if reference.len() < 2 || positions_named(&e.location) != Some(reference.len() as i128) {
        return false;
    }
    if !spelt_in(&e.allele, b"ACGTN") || allele.len() <= reference.len() {
        return false;
    }
    record_span_placement(&positions, &reference, &allele, expected, actual).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_generated_name_with_a_position_in_other_digits_spells_no_record() {
        assert_eq!(uploaded_alleles("21_100_A/G"), Some(("A", vec!["G"])));
        assert_eq!(uploaded_position("21_100_A/G"), Some(100));
        assert_eq!(uploaded_alleles("21_١٠٠_A/G"), None);
        assert_eq!(uploaded_position("21_١٠٠_A/G"), None);
    }

    #[test]
    fn an_identifier_spells_no_alleles_and_a_generated_name_spells_the_records() {
        assert_eq!(uploaded_alleles("1993671"), None);
        assert_eq!(uploaded_alleles("rs_12_x"), None);
        assert_eq!(uploaded_alleles("21_100 _A/T"), None);
        assert_eq!(
            uploaded_alleles("21_10542414_A/T/AG"),
            Some(("A", vec!["T", "AG"]))
        );
        assert_eq!(
            uploaded_alleles("chr_un_7_100_-/<INS>"),
            Some(("-", vec!["<INS>"]))
        );
        assert_eq!(
            uploaded_alleles("21_100_A/T_5_C/G"),
            Some(("C", vec!["G"])),
            "the last position and alleles the name can end in"
        );
    }

    #[test]
    fn the_records_own_alleles_give_the_alleles_shape() {
        assert_eq!(allele_is_indel("GACCTG", &["G"], "-"), Some(true));
        assert_eq!(allele_is_indel("A", &["T"], "T"), Some(false));
        assert_eq!(allele_is_indel("TCT", &["TCTC"], "C"), Some(true));
        assert_eq!(allele_is_indel("A", &["AT", "G"], "AT"), Some(true));
        assert_eq!(allele_is_indel("A", &["AT", "G"], "G"), Some(false));
        assert_eq!(allele_is_indel("AT", &["A", "AG"], "-"), Some(true));
        assert_eq!(allele_is_indel("AT", &["A", "G"], "-"), None);
        assert_eq!(allele_is_indel("A", &["T", "G"], "C"), None);
        assert_eq!(allele_is_indel("A", &["<DEL>"], "deletion"), None);
        assert_eq!(allele_is_indel("A", &["T"], ""), None);
        assert_eq!(allele_is_indel("AT", &["A", "*"], "-"), Some(true));
    }

    /// A row of one transcript under an identifier: the Location, the allele, the terms and
    /// the two position columns.
    fn row(location: &str, allele: &str, consequence: &str, cdna: &str, cds: &str) -> PairRow {
        pair_row(&[
            "1993671",
            location,
            allele,
            "ENSG1",
            "ENST00000635895",
            "Transcript",
            consequence,
            cdna,
            cds,
        ])
    }

    #[test]
    fn a_position_column_carries_its_start_and_its_end() {
        for (position, bounds) in [
            ("1705", (true, true)),
            ("1704-1705", (true, true)),
            ("?-7057", (false, true)),
            ("1705-?", (true, false)),
            ("-", (false, false)),
            ("", (false, false)),
            ("-/7140", (false, false)),
            ("120-7100/7140", (true, true)),
            ("?-6980/7023", (false, true)),
        ] {
            assert_eq!(position_bounds(position), bounds, "{position}");
        }
    }

    #[test]
    fn a_location_and_an_allele_tell_an_insertion_or_a_deletion() {
        assert_eq!(location_span("3:36993662-36993663"), Some(2));
        assert_eq!(location_span("3:36993663"), None);
        assert_eq!(location_span("3:5-x"), None);
        assert_eq!(location_span("36993662-36993663"), Some(2));
        assert!(indel_row(&row("3:36993663", "-", "", "", "")));
        assert!(indel_row(&row("3:36993662-36993663", "A", "", "", "")));
        assert!(!indel_row(&row("3:36993663", "A", "", "", "")));
        assert!(!indel_row(&row("3:36993662-36993663", "AC", "", "", "")));
        assert!(!indel_row(&row("3:36993662-36993663", "R", "", "", "")));
        let mut named = row("3:36993663", "A", "", "", "");
        named.indel = Some(true);
        assert!(
            indel_row(&named),
            "the record's own alleles come before the Location"
        );
    }

    #[test]
    fn the_terms_of_a_consequence_are_a_set() {
        assert!(is(&terms("b,a&b,,"), &["a", "b"]));
        assert!(is(&terms(""), &[""]));
        let (expected, gone, come) = (terms("a,b"), terms("a,b"), terms("a,c"));
        assert!(expected == gone && gone != come);
    }

    #[test]
    fn the_genomic_span_class_reads_a_tested_row_at_an_exon_end() {
        let terms = |term: &str| format!("splice_region_variant,{term},NMD_transcript_variant");
        let pair = |e: &str, a: &str, cdna: &str, cds: &str| {
            stop_codon_overlap_genomic_span_pair(
                &row("3:36993663", "-", &terms(e), cdna, cds),
                &row("3:36993663", "-", &terms(a), cdna, cds),
            )
        };
        let (coding, retained) = ("coding_sequence_variant", "stop_retained_variant");
        assert!(pair(coding, retained, "1705", "1630"));
        assert!(pair(coding, retained, "1705/1800", "1630/1653"));
        assert!(!pair(retained, coding, "1705", "1630"));
        assert!(!pair(coding, retained, "1705-?", "1630-?"));
        assert!(!pair(coding, retained, "1705", "-"));
        assert!(!pair(
            coding,
            "stop_retained_variant,stop_lost",
            "1705",
            "1630"
        ));
        assert!(!stop_codon_overlap_genomic_span_pair(
            &row("3:36993663", "-", coding, "1705", "1630"),
            &row("3:36993663", "-", retained, "1705", "1630"),
        ));
        assert!(!stop_codon_overlap_genomic_span_pair(
            &row("3:36993663", "A", &terms(coding), "1705", "1630"),
            &row("3:36993663", "A", &terms(retained), "1705", "1630"),
        ));
    }

    #[test]
    fn the_untested_class_reads_rows_whose_stop_codon_test_is_not_made() {
        let pair = stop_retained_without_stop_codon_test_pair;
        let acceptor = |term: &str| format!("splice_acceptor_variant,{term}");
        let (location, retained) = ("19:13208046-13208050", "stop_retained_variant");
        for (cdna, cds) in [
            ("?-7057", "-"),
            ("-", "-"),
            ("1705-?", "1630-?"),
            ("120-7100", "-"),
        ] {
            let reference = row(location, "-", &acceptor(retained), cdna, cds);
            for term in ["coding_sequence_variant", "frameshift_variant"] {
                assert!(pair(
                    &reference,
                    &row(location, "-", &acceptor(term), cdna, cds)
                ));
            }
            let nothing = row(location, "-", "splice_acceptor_variant", cdna, cds);
            assert!(pair(&reference, &nothing), "{cdna} {cds}");
        }
        let lettered = |term: &str| {
            row(
                "7:152648662-152648663",
                "ACRT",
                term,
                "7038-7039",
                "6961-6962",
            )
        };
        assert!(pair(&lettered(retained), &lettered("frameshift_variant")));
        let tested = |term: &str| row(location, "-", &acceptor(term), "7055-7059", "6978-6982");
        assert!(!pair(&tested(retained), &tested("coding_sequence_variant")));
        let untested = |term: &str| row(location, "-", &acceptor(term), "?-7057", "-");
        assert!(!pair(&untested(retained), &untested("stop_lost")));
        assert!(!pair(
            &untested(retained),
            &untested("coding_sequence_variant,frameshift_variant")
        ));
        assert!(!pair(
            &untested("coding_sequence_variant"),
            &untested(retained)
        ));
        assert!(!pair(
            &untested(retained),
            &tested("coding_sequence_variant")
        ));
        let substitution = |term: &str| row("19:13208046", "A", term, "-", "-");
        assert!(!pair(
            &substitution(retained),
            &substitution("coding_sequence_variant")
        ));
    }

    #[test]
    fn the_three_base_class_reads_rows_whose_stop_codon_test_is_made() {
        let pair = stop_retained_three_base_stop_codon_test_pair;
        let location = "2:190970701-190970707";
        let utr = |term: &str, cdna: &str| {
            row(
                location,
                "-",
                &format!("{term},3_prime_UTR_variant"),
                cdna,
                "2249-?",
            )
        };
        let (retained, coding) = ("stop_retained_variant", "coding_sequence_variant");
        assert!(pair(&utr(retained, "2558-2564"), &utr(coding, "2558-2564")));
        let coded = |term: &str| row(location, "-", term, "961-967", "264-270");
        assert!(pair(&coded(retained), &coded("frameshift_variant")));
        assert!(pair(
            &coded("inframe_deletion,stop_retained_variant"),
            &coded("frameshift_variant")
        ));
        assert!(pair(
            &coded("protein_altering_variant,stop_retained_variant"),
            &coded("stop_gained,frameshift_variant")
        ));
        assert!(!pair(&utr(retained, "2558-?"), &utr(coding, "2558-?")));
        assert!(!pair(&utr(retained, "2558-2564"), &utr(coding, "2558-?")));
        assert!(!pair(
            &utr(coding, "2558-2564"),
            &utr(retained, "2558-2564")
        ));
        assert!(!pair(
            &coded("inframe_deletion,stop_retained_variant"),
            &coded(coding)
        ));
        assert!(!pair(
            &coded("missense_variant,stop_retained_variant"),
            &coded("frameshift_variant")
        ));
        assert!(!pair(
            &utr(retained, "2558-2564"),
            &utr("stop_lost", "2558-2564")
        ));
        let lettered = |term: &str| {
            row(
                "7:152648662-152648663",
                "ACGTACGN",
                term,
                "1213-1214",
                "1100-1101",
            )
        };
        assert!(!pair(&lettered(retained), &lettered("frameshift_variant")));
    }

    #[test]
    fn the_removed_codon_class_reads_rows_without_the_tests_coordinates() {
        let pair = stop_codon_removed_without_stop_codon_test_pair;
        let location = "19:13208046-13208050";
        let acceptor = |term: &str, cdna: &str, cds: &str| {
            row(
                location,
                "-",
                &format!("splice_acceptor_variant,{term}"),
                cdna,
                cds,
            )
        };
        let (retained, lost) = ("stop_retained_variant", "stop_lost");
        for (cdna, cds) in [
            ("?-7057", "-"),
            ("-", "-"),
            ("1705-?", "1630-?"),
            ("120-7100", "-"),
            ("21-99", "?-69"),
        ] {
            assert!(
                pair(&acceptor(retained, cdna, cds), &acceptor(lost, cdna, cds)),
                "{cdna} {cds}"
            );
        }
        for (cdna, cds) in [("7055-7059", "6978-6982"), ("7055-7059", "6978-?")] {
            assert!(
                !pair(&acceptor(retained, cdna, cds), &acceptor(lost, cdna, cds)),
                "{cdna} {cds}"
            );
        }
        assert!(!pair(
            &acceptor(retained, "?-7057", "-"),
            &acceptor("stop_lost,frameshift_variant", "?-7057", "-")
        ));
        assert!(!pair(
            &acceptor(lost, "?-7057", "-"),
            &acceptor(retained, "?-7057", "-")
        ));
        assert!(!pair(
            &acceptor(retained, "?-7057", "-"),
            &acceptor("coding_sequence_variant", "?-7057", "-")
        ));
        let substitution = |term: &str| row("19:13208046", "A", term, "-", "-");
        assert!(!pair(&substitution(retained), &substitution(lost)));
    }

    #[test]
    fn the_insertion_range_class_reads_an_insertion_that_keeps_the_stop_codon() {
        let pair = stop_codon_overlap_insertion_range_pair;
        let inserted = |allele: &str, terms: &str, cdna: &str, cds: &str| {
            row("19:13208046-13208047", allele, terms, cdna, cds)
        };
        let whole = |terms: &str| inserted("TTGA", terms, "7054-7055", "6977-6978");
        for (reference, candidate) in [
            (
                "inframe_insertion,coding_sequence_variant",
                "inframe_insertion,stop_retained_variant",
            ),
            (
                "stop_gained,frameshift_variant",
                "inframe_insertion,stop_retained_variant",
            ),
            (
                "frameshift_variant",
                "inframe_insertion,stop_retained_variant",
            ),
            ("stop_gained,frameshift_variant", "stop_retained_variant"),
        ] {
            assert!(pair(&whole(reference), &whole(candidate)), "{reference}");
            assert!(!pair(&whole(candidate), &whole(reference)), "{reference}");
            assert!(pair(
                &inserted("TNGA", reference, "7054-7055", "6978-?"),
                &inserted("TNGA", candidate, "7054-7055", "6978-?")
            ));
            assert!(!pair(
                &inserted("TTGA", reference, "7054-?", "6977-?"),
                &inserted("TTGA", candidate, "7054-?", "6977-?")
            ));
            assert!(!pair(
                &inserted("TRGA", reference, "7054-7055", "6977-6978"),
                &inserted("TRGA", candidate, "7054-7055", "6977-6978")
            ));
            let deleted =
                |terms: &str| row("19:13208046-13208049", "-", terms, "7055-7058", "6978-6981");
            assert!(!pair(&deleted(reference), &deleted(candidate)));
        }
        let region = |term: &str| whole(&format!("inframe_insertion,splice_region_variant,{term}"));
        assert!(!pair(
            &region("coding_sequence_variant"),
            &region("stop_retained_variant")
        ));
        assert!(stop_codon_overlap_genomic_span_pair(
            &region("coding_sequence_variant"),
            &region("stop_retained_variant")
        ));
        assert!(!pair(
            &whole("frameshift_variant"),
            &whole("protein_altering_variant,stop_retained_variant")
        ));
        assert!(!pair(
            &whole("stop_gained"),
            &whole("stop_retained_variant")
        ));
    }

    #[test]
    fn a_reference_allele_is_empty_once_the_shared_first_base_is_dropped() {
        for (reference, alternates) in [
            ("-", vec!["TGA"]),
            ("A", vec!["ATGA"]),
            ("-", vec!["CG", "T"]),
            ("A", vec!["AT", "AG"]),
            ("A", vec!["AT", "*"]),
        ] {
            assert!(
                reference_allele_is_empty(reference, &alternates),
                "{reference}"
            );
        }
        for (reference, alternates) in [
            ("AC", vec!["GTT"]),
            ("AC", vec!["AGT"]),
            ("AC", vec!["A", "T"]),
            ("A", vec!["AC", "T"]),
            ("A", vec!["T"]),
            ("TCT", vec!["TCTC"]),
            ("GACCTG", vec!["G"]),
        ] {
            assert!(
                !reference_allele_is_empty(reference, &alternates),
                "{reference}"
            );
        }
    }

    #[test]
    fn the_insertion_range_class_takes_an_insertion_between_two_positions_alone() {
        let taken = |name: &str, location: &str, allele: &str| {
            let row = |terms: &str| {
                pair_row(&[
                    name,
                    location,
                    allele,
                    "ENSG1",
                    "ENST00000635895",
                    "Transcript",
                    terms,
                    "7054-7055",
                    "6977-6978",
                ])
            };
            stop_codon_overlap_insertion_range_pair(
                &row("coding_sequence_variant"),
                &row("stop_retained_variant"),
            )
        };
        let two = "19:13208046-13208047";
        for (name, location, allele, shape) in [
            (
                "19_13208046_AC/GTT",
                two,
                "GTT",
                "two bases replaced by three",
            ),
            (
                "19_13208046_AC/AGT",
                two,
                "AGT",
                "the same, the first base shared",
            ),
            (
                "19_13208046_AC/A/T",
                two,
                "A",
                "a deletion written as a sequence",
            ),
            (
                "19_13208046_A/AC/T",
                "19:13208046",
                "AC",
                "an insertion kept with its first base",
            ),
            (
                "19_13208046_-/TGA",
                "19:13208046",
                "TGA",
                "a spelt insertion at one position",
            ),
            (
                "19_13208046_AC/A/T",
                two,
                "GTT",
                "an allele no alternate names, the reference of two",
            ),
            (
                "1993671",
                two,
                "CG",
                "a two-base allele under an identifier",
            ),
            ("1993671", "19:13208046-13208048", "TGA", "three positions"),
        ] {
            assert!(!taken(name, location, allele), "{shape}");
        }
        for (name, allele, shape) in [
            (
                "19_13208046_-/TGA",
                "TGA",
                "an insertion, the reference allele empty",
            ),
            (
                "19_13208046_A/ATGA",
                "TGA",
                "an insertion spelt with its shared first base",
            ),
            (
                "19_13208046_-/CG/T",
                "CG",
                "a two-base insertion, the alleles spelt",
            ),
            ("1993671", "TGA", "an insertion under an identifier"),
            (
                "1993671",
                "G",
                "under an identifier, two bases replaced by one",
            ),
            (
                "1993671",
                "CGT",
                "under an identifier, two bases replaced by three",
            ),
        ] {
            assert!(taken(name, two, allele), "{shape}");
        }
    }

    #[test]
    fn the_three_base_class_takes_the_retained_term_alone_or_beside_one_read_away_term() {
        let pair = stop_retained_three_base_stop_codon_test_pair;
        let coded = |terms: &str| row("2:190970701-190970707", "-", terms, "961-967", "264-270");
        let retained = "stop_retained_variant";
        for candidate in ["frameshift_variant", "stop_gained,frameshift_variant"] {
            assert!(pair(&coded(retained), &coded(candidate)), "alone");
            for one in FRAMESHIFT_READ_TERMS {
                let beside = format!("{one},{retained}");
                assert!(pair(&coded(&beside), &coded(candidate)), "{one}");
            }
            for several in [
                "inframe_insertion,inframe_deletion",
                "inframe_insertion,protein_altering_variant",
                "inframe_deletion,protein_altering_variant",
                "inframe_insertion,inframe_deletion,protein_altering_variant",
            ] {
                let beside = format!("{several},{retained}");
                assert!(!pair(&coded(&beside), &coded(candidate)), "{several}");
            }
        }
    }

    #[test]
    fn a_row_carries_its_columns_and_the_shape_of_its_allele() {
        let named = pair_row(&[
            "19_13208045_GACCTG/G",
            "19:13208046-13208050",
            "-",
            "ENSG1",
            "ENST00000635895",
            "Transcript",
            "stop_lost",
            "?-7057",
            "-",
        ]);
        assert_eq!(
            named,
            PairRow {
                location: "19:13208046-13208050".to_string(),
                allele: "-".to_string(),
                feature_type: "Transcript".to_string(),
                consequence: "stop_lost".to_string(),
                cdna_position: "?-7057".to_string(),
                cds_position: "-".to_string(),
                indel: Some(true),
                insertion: Some(false),
                reference: Some("ACCTG".to_string()),
            }
        );
        let identified = pair_row(&[
            "1993671",
            "19:13208046",
            "A",
            "ENSG1",
            "ENST00000635895",
            "Transcript",
            "stop_lost",
            "-",
            "-",
        ]);
        assert_eq!(
            (identified.indel, identified.insertion, identified.reference),
            (None, None, None)
        );
    }

    #[test]
    fn a_position_column_writes_its_start_and_its_end() {
        for (position, parts) in [
            ("999-1000", Some(("999", "1000"))),
            ("999-1000/2000", Some(("999", "1000"))),
            ("1705", Some(("1705", "1705"))),
            ("1705/1800", Some(("1705", "1705"))),
            ("1000-?", Some(("1000", "?"))),
            ("?-1000", Some(("?", "1000"))),
            ("?", Some(("?", "?"))),
            ("-", Some(("-", "-"))),
            ("", Some(("-", "-"))),
            ("-/7140", None),
            ("999-", None),
            ("-1000", None),
            ("999-1000-1001", None),
            ("9x9-1000", None),
            ("<absent>", None),
        ] {
            assert_eq!(position_parts(position), parts, "{position}");
        }
        assert_eq!(
            transcript_positions(&row("7:100998-101000", "GTGTG", "", "998-1000/2000", "-")),
            Some(["998", "1000", "-", "-"])
        );
        assert_eq!(
            transcript_positions(&row("7:100998-101000", "GTGTG", "", "998-1000", "9x9")),
            None
        );
    }

    #[test]
    fn a_generated_name_spells_the_records_position_before_its_alleles() {
        assert_eq!(uploaded_position("7_100999_TA/GTGTG"), Some(100999));
        assert_eq!(uploaded_position("chr_un_7_100_-/<INS>"), Some(100));
        assert_eq!(
            uploaded_position("21_100_A/T_5_C/G"),
            Some(5),
            "the last position and alleles the name can end in"
        );
        assert_eq!(uploaded_position("1993671"), None);
        assert_eq!(
            uploaded_position(&format!("7_{}_TA/GTGTG", "9".repeat(40))),
            None,
            "a position past 128 bits"
        );
    }

    #[test]
    fn a_location_names_its_positions() {
        assert_eq!(positions_named("7:100998-101000"), Some(3));
        assert_eq!(positions_named("7:101000"), Some(1));
        assert_eq!(positions_named("100999-101000"), Some(2));
        assert_eq!(positions_named("7:101000-100998"), None);
        assert_eq!(positions_named("7:100998-x"), None);
    }

    const TWO_POSITIONS: &str = "7:100999-101000";

    fn bases(sequence: &str) -> Option<String> {
        Some(sequence.to_string())
    }

    #[test]
    fn a_replacements_allele_stands_in_for_the_bases_of_its_location() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100999, "TA", &["GTGTG"], TWO_POSITIONS, "GTGTG"),
            bases("TA")
        );
        assert_eq!(
            replaced(100998, "CTA", &["CGTGTG"], TWO_POSITIONS, "GTGTG"),
            bases("TA"),
            "behind a shared first base"
        );
        assert_eq!(
            replaced(100999, "ta", &["gtgtg"], TWO_POSITIONS, "GTGTG"),
            bases("TA"),
            "the record's alleles in lower case"
        );
        assert_eq!(
            replaced(100999, "TA", &["C", "GTGTG"], TWO_POSITIONS, "GTGTG"),
            bases("TA"),
            "beside another alternate"
        );
        assert_eq!(
            replaced(101000, "T", &["GTG"], "7:101000", "GTG"),
            bases("T"),
            "a Location of one position"
        );
    }

    #[test]
    fn a_pure_insertions_allele_stands_in_for_no_base() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100999, "T", &["TGTG"], TWO_POSITIONS, "GTG"),
            bases("")
        );
        assert_eq!(
            replaced(100998, "CT", &["CTGTG"], TWO_POSITIONS, "GTG"),
            bases(""),
            "after the record's last base"
        );
        assert_eq!(
            replaced(100999, "TAC", &["TGTGAC"], "7:100999-101001", "GTG"),
            None,
            "a Location of three positions"
        );
    }

    #[test]
    fn a_deletions_allele_stands_in_for_the_bases_it_removes() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100997, "CTAC", &["C"], "7:100998-101000", "-"),
            bases("TAC")
        );
        assert_eq!(
            replaced(100998, "CTA", &["C"], TWO_POSITIONS, "-"),
            bases("TA")
        );
        assert_eq!(
            replaced(100999, "TA", &["TA"], TWO_POSITIONS, "-"),
            None,
            "nothing inserted between two neighbouring positions"
        );
    }

    #[test]
    fn a_record_that_writes_no_such_row_gives_no_bases() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100999, "TA", &["GTGTG", "C"], TWO_POSITIONS, "GTG"),
            None,
            "no alternate writes the allele"
        );
        assert_eq!(
            replaced(100999, "TA", &["GTG", "TGTGA"], TWO_POSITIONS, "GTG"),
            None,
            "one alternate replaces the two bases and one inserts between them"
        );
        assert_eq!(
            replaced(101000, "TA", &["GTG"], TWO_POSITIONS, "GTG"),
            None,
            "a Location that starts before the record"
        );
        assert_eq!(
            replaced(100999, "TA", &["GTR"], TWO_POSITIONS, "GTR"),
            None,
            "an allele that is not a sequence"
        );
        assert_eq!(replaced(100999, "TA", &["GTG"], "7:100999-x", "GTG"), None);
    }

    fn letters(sequence: &str) -> Vec<char> {
        sequence.chars().collect()
    }

    #[test]
    fn a_record_is_cut_into_the_runs_of_indices_at_which_its_alleles_differ() {
        let regions = |reference: &str, allele: &str| {
            differing_regions(&letters(reference), &letters(allele))
        };
        assert_eq!(regions("TA", "GAGT"), [(0, 0), (2, 3)]);
        assert_eq!(regions("TA", "GTGTG"), [(0, 4)]);
        assert_eq!(regions("TA", "TA"), [(0, 1)], "no index differs");
        assert_eq!(regions("T", "GTG"), [(0, 0)], "a reference of one base");
        assert_eq!(regions("TAC", "G"), [(0, 2)], "an allele of one base");
    }

    /// Whether the class takes the reference's and vep-rs's row of one transcript, the
    /// Location, the allele and the two position columns the same on both: the reference
    /// writes the donor term and vep-rs `splice_region_variant` in its place.
    fn past_the_span_rows(location: &str, allele: &str, cdna: &str, cds: &str) -> bool {
        let row = |terms: &str| row(location, allele, terms, cdna, cds);
        splice_term_past_the_record_span_pair(
            &row("splice_donor_variant,frameshift_variant"),
            &row("frameshift_variant,splice_region_variant"),
        )
    }

    /// Whether the class takes two rows of three positions replaced by five letters in the
    /// transcript's sequence, the reference's terms and vep-rs's as given.
    fn past_the_span_terms(reference: &str, candidate: &str) -> bool {
        let row = |terms: &str| row("7:100998-101000", "GTGTG", terms, "998-1000", "898-900");
        splice_term_past_the_record_span_pair(&row(reference), &row(candidate))
    }

    #[test]
    fn the_past_the_span_class_takes_a_donor_term_gone_and_the_region_term_come() {
        assert!(past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "998-1000",
            "898-900"
        ));
        assert!(past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "998-1000/2000",
            "898-900/1800"
        ));
    }

    #[test]
    fn the_past_the_span_class_takes_terms_of_the_six_gone_and_none_come() {
        assert!(past_the_span_terms(
            "splice_donor_variant,frameshift_variant,splice_donor_region_variant",
            "frameshift_variant"
        ));
        assert!(past_the_span_terms(
            "splice_polypyrimidine_tract_variant,frameshift_variant",
            "frameshift_variant"
        ));
        assert!(past_the_span_terms(
            "splice_donor_variant,frameshift_variant,intron_variant",
            "frameshift_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_span_across_an_exons_end() {
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "999-?",
            "899-?"
        ));
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "?-1000",
            "?-900"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_span_inside_an_intron() {
        assert!(!past_the_span_rows("7:100998-101000", "GTGTG", "-", "-"));
    }

    #[test]
    fn the_past_the_span_class_refuses_an_allele_no_longer_than_the_span() {
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GTG",
            "998-1000",
            "898-900"
        ));
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GT",
            "998-1000",
            "898-900"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_location_of_one_position() {
        assert!(!past_the_span_rows("7:101000", "GTGTG", "1000", "900"));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_location_of_two_neighbouring_positions() {
        assert!(!past_the_span_rows(
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            "899-900"
        ));
    }

    #[test]
    fn the_first_column_tells_two_replaced_bases_from_an_insertion_between_them() {
        let spelt = |name: &str| {
            let row = |terms: &str| {
                pair_row(&[
                    name,
                    TWO_POSITIONS,
                    "GTGTG",
                    "ENSG1",
                    "ENST00000635895",
                    "Transcript",
                    terms,
                    "999-1000",
                    "899-900",
                ])
            };
            splice_term_past_the_record_span_pair(
                &row("splice_donor_variant,frameshift_variant"),
                &row("frameshift_variant,splice_region_variant"),
            )
        };
        assert!(spelt("7_100999_TA/GTGTG"), "two bases replaced by five");
        assert!(!spelt("7_100999_T/TGTGTG"), "an insertion of five");
    }

    #[test]
    fn the_past_the_span_class_refuses_a_deletion() {
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "-",
            "998-1000",
            "898-900"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_candidate_that_keeps_one_of_the_six() {
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant,splice_donor_region_variant",
            "frameshift_variant,splice_donor_region_variant"
        ));
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant,intron_variant",
            "frameshift_variant,intron_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_donor_region_term_that_comes() {
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant",
            "frameshift_variant,splice_donor_region_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_the_region_term_where_the_tract_term_alone_went() {
        assert!(!past_the_span_terms(
            "splice_polypyrimidine_tract_variant,frameshift_variant",
            "frameshift_variant,splice_region_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_further_differing_term() {
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant",
            "stop_gained,frameshift_variant,splice_region_variant"
        ));
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant,stop_lost",
            "frameshift_variant,splice_region_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_rows_whose_transcript_positions_differ() {
        let row =
            |terms: &str, cdna: &str, cds: &str| row("7:100998-101000", "GTGTG", terms, cdna, cds);
        let reference = row(
            "splice_donor_variant,frameshift_variant",
            "998-1000",
            "898-900",
        );
        let candidate = "frameshift_variant,splice_region_variant";
        assert!(!splice_term_past_the_record_span_pair(
            &reference,
            &row(candidate, "997-999", "898-900")
        ));
        assert!(!splice_term_past_the_record_span_pair(
            &reference,
            &row(candidate, "998-1000", "897-899")
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_the_terms_the_other_way_round() {
        assert!(!past_the_span_terms(
            "frameshift_variant,splice_region_variant",
            "splice_donor_variant,frameshift_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_takes_a_row_that_differs_by_a_stop_codon_classes_terms_as_well() {
        let row = |terms: &str| row("7:100998-101000", "GTGTG", terms, "998-1000", "-");
        let reference = row("splice_donor_variant,stop_retained_variant");
        let candidate = row("stop_lost,splice_region_variant");
        assert!(splice_term_past_the_record_span_pair(
            &reference, &candidate
        ));
        assert!(
            stop_codon_removed_without_stop_codon_test_pair(
                &row("stop_retained_variant,splice_region_variant"),
                &candidate
            ),
            "the rest of the difference is that class's"
        );
        assert!(
            !splice_term_past_the_record_span_pair(
                &reference,
                &row("stop_lost,missense_variant,splice_region_variant")
            ),
            "a further term outside both classes"
        );
    }

    /// The terms of a reference row whose donor, donor-region and intron terms come from
    /// positions past the span, and of vep-rs's row with `splice_region_variant` in their
    /// place.
    const DONOR_SIDE_TERMS: &str =
        "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant";
    const REGION_TERM_BACK: &str = "frameshift_variant,splice_region_variant";

    /// A row of one transcript whose allele stands in for the reference bases `reference`.
    fn record_row(
        reference: &str,
        location: &str,
        allele: &str,
        terms: &str,
        cdna: &str,
    ) -> PairRow {
        PairRow {
            reference: bases(reference),
            ..row(location, allele, terms, cdna, "-")
        }
    }

    /// Whether the class takes the reference's and vep-rs's row of one transcript on the
    /// record that replaces `reference` at `location` by `allele`, their terms as given.
    fn past_the_span_record(
        reference: &str,
        location: &str,
        allele: &str,
        cdna: &str,
        terms: (&str, &str),
    ) -> bool {
        let row = |terms: &str| record_row(reference, location, allele, terms, cdna);
        splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1))
    }

    #[test]
    fn with_its_record_the_class_takes_a_span_inside_an_exon_of_two_positions_and_of_three() {
        let terms = (DONOR_SIDE_TERMS, REGION_TERM_BACK);
        assert!(past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            terms
        ));
        assert!(past_the_span_record(
            "TAC",
            "7:100998-101000",
            "GTGTGT",
            "998-1000",
            terms
        ));
        assert!(
            past_the_span_record(
                "TA",
                TWO_POSITIONS,
                "GTGT",
                "901-902",
                (
                    "splice_acceptor_variant,frameshift_variant",
                    REGION_TERM_BACK
                )
            ),
            "the acceptor term, on the reverse strand"
        );
    }

    #[test]
    fn with_its_record_the_class_reads_the_reference_bases_of_either_row() {
        let reference = row(TWO_POSITIONS, "GTGTG", DONOR_SIDE_TERMS, "999-1000", "-");
        let candidate = record_row("TA", TWO_POSITIONS, "GTGTG", REGION_TERM_BACK, "999-1000");
        assert!(splice_term_past_the_record_span_pair(
            &reference, &candidate
        ));
    }

    #[test]
    fn with_its_record_the_class_takes_a_span_across_an_exons_end() {
        assert!(past_the_span_record(
            "TA",
            "7:101000-101001",
            "GTGTG",
            "1000-?",
            (
                "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                "splice_donor_variant,coding_sequence_variant,intron_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_takes_a_span_in_an_intron() {
        assert!(past_the_span_record(
            "TA",
            "7:101003-101004",
            "GTGT",
            "-",
            (
                "splice_donor_5th_base_variant,intron_variant",
                "splice_donor_region_variant,intron_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_takes_the_intron_term_gone_one_to_three_bases_before_the_intron() {
        for (location, allele, cdna, intron_start) in [
            (TWO_POSITIONS, "GTGTG", "999-1000", 2),
            ("7:100998-100999", "GTGTGTG", "998-999", 3),
            ("7:100997-100998", "GTGTGTGT", "997-998", 4),
        ] {
            assert!(
                past_the_span_record(
                    "TA",
                    location,
                    allele,
                    cdna,
                    (DONOR_SIDE_TERMS, REGION_TERM_BACK)
                ),
                "{cdna}"
            );
            let (start, end) = cdna.split_once('-').unwrap();
            let place = record_span_placement(
                &[start, end, "-", "-"],
                &letters("TA"),
                &letters(allele),
                DONOR | DONOR_REGION | INTRON,
                REGION,
            );
            assert_eq!(
                place.map(|(forward, start, ..)| (forward, start)),
                Some((true, intron_start)),
                "{cdna}"
            );
        }
    }

    #[test]
    fn with_its_record_the_class_takes_a_donor_region_term_gone_beside_a_kept_donor_term() {
        assert!(past_the_span_record(
            "TA",
            "7:101001-101002",
            "GTGT",
            "-",
            (
                "splice_donor_variant,splice_donor_region_variant,intron_variant",
                "splice_donor_variant,intron_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_a_pure_insertion() {
        assert!(!past_the_span_record(
            "",
            TWO_POSITIONS,
            "GTG",
            "999-1000",
            ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_an_allele_as_long_as_the_reference_bases() {
        assert!(!past_the_span_record(
            "TAC",
            "7:100998-101000",
            "GTG",
            "998-1000",
            ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_one_reference_base() {
        assert!(!past_the_span_record(
            "T",
            "7:101000",
            "GTG",
            "1000",
            ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
        ));
    }

    /// The record replaces two bases across an exon's end by four letters whose second is the
    /// reference's: the intron's first base lies in the span at an index where the alleles
    /// agree, so vep-rs keeps the donor term.
    #[test]
    fn with_its_record_the_class_refuses_a_term_gone_whose_window_holds_a_base_of_the_span() {
        let record = |candidate: &str| {
            past_the_span_record(
                "TA",
                "7:101000-101001",
                "GAGT",
                "1000-?",
                (
                    "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                    candidate,
                ),
            )
        };
        assert!(record(
            "splice_donor_variant,coding_sequence_variant,intron_variant"
        ));
        assert!(!record("coding_sequence_variant,intron_variant"));
    }

    #[test]
    fn with_its_record_the_class_refuses_rows_no_place_of_the_intron_gives() {
        assert!(!past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (DONOR_SIDE_TERMS, "frameshift_variant")
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_a_further_differing_term() {
        assert!(!past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (
                DONOR_SIDE_TERMS,
                "stop_gained,frameshift_variant,splice_region_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_rows_whose_transcript_positions_differ() {
        let row = |terms: &str, cdna: &str| record_row("TA", TWO_POSITIONS, "GTGTG", terms, cdna);
        assert!(!splice_term_past_the_record_span_pair(
            &row(DONOR_SIDE_TERMS, "999-1000"),
            &row(REGION_TERM_BACK, "499-500")
        ));
    }

    /// The places of an intron beside a span that two cDNA positions allow, in their order.
    fn places(span: i64, cdna: (&str, &str)) -> Vec<(bool, i64, i64)> {
        let mut places = Vec::new();
        intron_placements(
            span,
            2,
            &[cdna.0, cdna.1, "-", "-"],
            |forward, start, end| {
                places.push((forward, start, end));
                false
            },
        );
        places
    }

    #[test]
    fn two_numbers_as_far_apart_as_the_span_put_the_intron_after_the_span_or_before_it() {
        let places = places(3, ("998", "1000"));
        let starts = |strand: bool, after: bool| -> BTreeSet<i64> {
            let on = |(forward, start, _): &(bool, i64, i64)| {
                *forward == strand && (*start > 0) == after
            };
            places
                .iter()
                .filter(|place| on(place))
                .map(|place| place.1)
                .collect()
        };
        for strand in [true, false] {
            assert_eq!(
                starts(strand, true),
                (3..10).collect(),
                "one to seven bases on"
            );
            assert_eq!(starts(strand, false), (-16..0).collect(), "before the span");
        }
        assert!(places.contains(&(true, 3, 3)), "an intron of one base");
        assert!(
            places.contains(&(false, 9, 10)),
            "a two-base intron seven bases on"
        );
        assert!(
            places.contains(&(true, 3, FAR)),
            "an intron longer than any window"
        );
        assert!(
            places.contains(&(true, -2, -1)),
            "a two-base intron that ends beside the span"
        );
        let before = places.iter().filter(|place| place.1 < 0);
        assert!(before
            .clone()
            .all(|(_, start, end)| start <= end && *end < 0));
        assert_eq!(before.count(), 2 * (1..=16).sum::<usize>());
        assert_eq!(places.len(), 1406);
    }

    #[test]
    fn two_numbers_closer_than_the_span_put_an_intron_of_the_missing_length_inside_it() {
        assert_eq!(
            places(5, ("998", "1000")),
            [(true, 1, 2), (true, 2, 3), (false, 1, 2), (false, 2, 3)],
            "two bases of five, the span's ends kept"
        );
        assert_eq!(
            places(3, ("999", "1000")),
            [(true, 1, 1), (false, 1, 1)],
            "the middle base of three"
        );
        assert_eq!(
            places(2, ("1000", "1000")),
            [],
            "no base inside a span of two"
        );
        assert_eq!(
            places(3, ("997", "1000")),
            [],
            "numbers further apart than the span"
        );
        assert_eq!(
            places(3, ("1000", "998")),
            [],
            "the second number below the first"
        );
    }

    /// Eight bases inside an intron that begins two bases before the span, on the reverse
    /// strand: the span's first six bases lie in a window of `splice_region` and all of them
    /// in the tract's.
    #[test]
    fn the_last_region_assigned_decides_the_references_region_term_and_any_region_vep_rs_s() {
        let rows = |regions: &[(i64, i64)]| intron_rows(regions, 8, -2, FAR, false, false);
        let two_runs = differing_regions(&letters("TACGATCA"), &letters("GACGATCAGT"));
        assert_eq!(two_runs, [(0, 0), (8, 9)]);
        let (reference, candidate) = rows(&two_runs);
        assert_eq!(
            (reference, candidate),
            (TRACT | INTRON, TRACT | INTRON | REGION),
            "the first run in a window and the last in none"
        );
        assert_eq!(candidate & TRACT, TRACT, "the span holds the tract's bases");
        assert_eq!(
            rows(&[(0, 0), (4, 9)]),
            (TRACT | INTRON | REGION, TRACT | INTRON | REGION),
            "the last run in a window"
        );
        assert_eq!(
            rows(&[(7, 9)]),
            (TRACT | INTRON, TRACT | INTRON),
            "no run in a window"
        );
        assert_eq!(
            intron_rows(&[(0, 0), (2, 2), (6, 8)], 6, 2, 101, true, false),
            (
                DONOR | DONOR_5TH_BASE | INTRON,
                DONOR | DONOR_REGION | INTRON
            ),
            "a run on the intron's first two bases between two runs in a window"
        );
    }

    /// Eight bases of an intron replaced by eleven letters alike at four of them, five runs in
    /// all: the last run lies in a window of `splice_region` and the span holds bases of the
    /// tract's window at every place that gives the reference's row.
    #[test]
    fn with_its_record_the_class_refuses_the_tract_term_gone_from_a_span_of_five_runs() {
        let (reference, allele) = ("TACGATCA", "GATGCTGAGTG");
        assert_eq!(
            differing_regions(&letters(reference), &letters(allele)),
            [(0, 0), (2, 2), (4, 4), (6, 6), (8, 10)]
        );
        assert!(!past_the_span_record(
            reference,
            "7:101007-101014",
            allele,
            "-",
            (
                "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
                "splice_region_variant,intron_variant"
            )
        ));
    }

    const FIFTH_BASE_TERMS: &str = "splice_donor_5th_base_variant,frameshift_variant";

    /// The windows of an intron of two bases run on past its last base: its fifth base is the
    /// first base past a span of two that follows it.
    #[test]
    fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_a_two_base_intron_before_the_span(
    ) {
        assert_eq!(
            intron_rows(&[(0, 4)], 2, -2, -1, true, false),
            (DONOR_5TH_BASE, DONOR_REGION)
        );
        assert!(past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (
                FIFTH_BASE_TERMS,
                "splice_donor_region_variant,frameshift_variant"
            )
        ));
    }

    /// On the reverse strand the fifth base of an intron of two bases that begins six bases
    /// after a span of two lies three bases past the span, and the eleven bases around the
    /// intron's end hold the span's last base.
    #[test]
    fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_a_two_base_intron_six_bases_after_the_span(
    ) {
        assert_eq!(
            intron_rows(&[(0, 5)], 2, 7, 8, false, false),
            (DONOR_5TH_BASE, REGION)
        );
        assert!(past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTGT",
            "999-1000",
            (FIFTH_BASE_TERMS, REGION_TERM_BACK)
        ));
    }

    /// Six positions whose two cDNA positions are four apart hold an intron of two bases: the
    /// alleles agree on it, and its fifth base is the first base past the span.
    #[test]
    fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_an_intron_inside_the_span() {
        let record = |cdna: &str| {
            past_the_span_record(
                "TACGAT",
                "7:100999-101004",
                "GACGATGT",
                cdna,
                (
                    "splice_donor_5th_base_variant,coding_sequence_variant",
                    "splice_donor_region_variant,coding_sequence_variant",
                ),
            )
        };
        assert_eq!(
            intron_rows(&[(0, 0), (6, 7)], 6, 2, 3, true, false),
            (DONOR_5TH_BASE, DONOR_REGION)
        );
        assert!(record("999-1002"));
        assert!(!record("999-998"), "the second number below the first");
    }
}
