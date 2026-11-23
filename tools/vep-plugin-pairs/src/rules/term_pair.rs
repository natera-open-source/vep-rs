// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rules over a row whose two Consequence sets differ: the pairs of sets the consequence
//! comparison sets aside, the row rules of the stop codon classes, of
//! `splice_term_past_the_record_span`, of `start_codon_removed_without_five_prime_utr` and of
//! `start_codon_kept_without_five_prime_utr`, and the plugin fields that follow such a pair
//! of sets or a pair a class of a coding end takes.

use super::coding_end::{coding_end_class, landed_rule_takes_beside_a_cache};
use super::coding_ends::terms::LandedRule;
use super::downstream::codes::{genetic_code, nucleotide_codon_residue};
use super::entry::{
    blank, consequence_set, entry_kind, feature_kind, field_named, Entry, Row, ABSENT,
};
use super::shared::{
    at_its_location, dash, indel_row, insertion_row, integer_of, location_ends, location_span,
    CODING_SEQUENCE_TERM, FRAMESHIFT_TERM, STOP_GAINED_TERM, STOP_RETAINED_TERM,
};
use super::{Consult, Res};
use crate::pytext;
use std::collections::BTreeSet;

pub const START_COOCCURRENCE_SWAP: &str = "start_cooccurrence_swap";
pub const SPLICE_FAMILY_SWAP: &str = "splice_family_swap";

const SPLICE_REGION_TERM: &str = "splice_region_variant";
const STOP_LOST_TERM: &str = "stop_lost";
const MISSENSE_TERM: &str = "missense_variant";
const START_LOST_TERM: &str = "start_lost";
const INFRAME_INSERTION_TERM: &str = "inframe_insertion";
const INTRON_TERM: &str = "intron_variant";

pub const IMPACT_FIELDS: [&str; 2] = ["IMPACT", "impact"];
pub const HGVSP_FIELDS: [&str; 2] = ["HGVSp", "hgvsp"];
/// The plugin fields written only on a row whose terms let the plugin run. Downstream.pm:97-99
/// asks for `frameshift_variant` and no term holding `splice`; NMD.pm:65 and :83 and
/// LOFTEE's LoF.pm:166-167 for one of `stop_gained`, `frameshift_variant` and the two
/// splice site terms.
pub const DOWNSTREAM_FIELDS: [&str; 2] = ["DownstreamProtein", "ProteinLengthChange"];
pub const TRUNCATING_TERM_FIELDS: [&str; 5] = ["NMD", "LoF", "LoF_info", "LoF_filter", "LoF_flags"];
pub const TRUNCATING_TERMS: [&str; 4] = [
    STOP_GAINED_TERM,
    FRAMESHIFT_TERM,
    "splice_donor_variant",
    "splice_acceptor_variant",
];
/// The plugin fields written only on a row with `missense_variant`: REVEL's (REVEL.pm:162)
/// and AlphaMissense's (AlphaMissense.pm:242), each a full name or a prefix ending in `_`.
pub const MISSENSE_TERM_FIELDS: [&str; 3] = ["REVEL", "am_", "alphamissense"];
/// The plugin fields written only on a row with one of `PROTEIN_CHANGE_TERMS`: ClinPred's
/// (ClinPred.pm:69 and :106) and the columns dbNSFP writes under its prefix
/// (dbNSFP.pm:159) where the run leaves the plugin's `consequence` parameter at its
/// default (dbNSFP.pm:157, :290-291 and :581). A column dbNSFP writes under its own name
/// is not among them: the name alone does not tell it from another plugin's field.
pub const PROTEIN_CHANGE_TERM_FIELDS: [&str; 2] = ["ClinPred", "dbNSFP_"];
pub const PROTEIN_CHANGE_TERMS: [&str; 4] = [
    MISSENSE_TERM,
    STOP_LOST_TERM,
    STOP_GAINED_TERM,
    START_LOST_TERM,
];

/// The Sequence Ontology terms of ensembl-variation's `%OVERLAP_CONSEQUENCES`
/// (Utils/Constants.pm) in the order of their rank, the most severe first, each with its
/// impact.
pub const TERMS_BY_RANK: [(&str, &str); 41] = [
    ("transcript_ablation", "HIGH"),
    ("splice_acceptor_variant", "HIGH"),
    ("splice_donor_variant", "HIGH"),
    ("stop_gained", "HIGH"),
    ("frameshift_variant", "HIGH"),
    ("stop_lost", "HIGH"),
    ("start_lost", "HIGH"),
    ("transcript_amplification", "HIGH"),
    ("feature_elongation", "HIGH"),
    ("feature_truncation", "HIGH"),
    ("inframe_insertion", "MODERATE"),
    ("inframe_deletion", "MODERATE"),
    ("missense_variant", "MODERATE"),
    ("protein_altering_variant", "MODERATE"),
    ("splice_donor_5th_base_variant", "LOW"),
    ("splice_region_variant", "LOW"),
    ("splice_donor_region_variant", "LOW"),
    ("splice_polypyrimidine_tract_variant", "LOW"),
    ("incomplete_terminal_codon_variant", "LOW"),
    ("start_retained_variant", "LOW"),
    ("stop_retained_variant", "LOW"),
    ("synonymous_variant", "LOW"),
    ("coding_sequence_variant", "MODIFIER"),
    ("mature_miRNA_variant", "MODIFIER"),
    ("5_prime_UTR_variant", "MODIFIER"),
    ("3_prime_UTR_variant", "MODIFIER"),
    ("non_coding_transcript_exon_variant", "MODIFIER"),
    ("intron_variant", "MODIFIER"),
    ("NMD_transcript_variant", "MODIFIER"),
    ("non_coding_transcript_variant", "MODIFIER"),
    ("coding_transcript_variant", "MODIFIER"),
    ("upstream_gene_variant", "MODIFIER"),
    ("downstream_gene_variant", "MODIFIER"),
    ("TFBS_ablation", "MODERATE"),
    ("TFBS_amplification", "MODIFIER"),
    ("TF_binding_site_variant", "MODIFIER"),
    ("regulatory_region_ablation", "MODIFIER"),
    ("regulatory_region_amplification", "MODIFIER"),
    ("regulatory_region_variant", "MODIFIER"),
    ("intergenic_variant", "MODIFIER"),
    ("sequence_variant", "MODIFIER"),
];

/// The terms of a Consequence text as a set: parted at `&` and `,`, the empty ones
/// dropped, an empty Consequence reading as one empty term.
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

/// The terms of a comma-joined text as a set, as it stands.
fn joined_terms(text: &str) -> BTreeSet<&str> {
    text.split(',').collect()
}

fn is(terms: &BTreeSet<&str>, members: &[&str]) -> bool {
    terms.len() == members.len() && members.iter().all(|member| terms.contains(member))
}

/// The terms the reference's row holds, those it alone holds and those the candidate's
/// alone holds.
fn term_differences<'a>(
    e: &'a Entry,
    a: &'a Entry,
) -> (BTreeSet<&'a str>, BTreeSet<&'a str>, BTreeSet<&'a str>) {
    let (expected, actual) = (terms(&e.consequence), terms(&a.consequence));
    let gone = expected.difference(&actual).copied().collect();
    let come = actual.difference(&expected).copied().collect();
    (expected, gone, come)
}

/// Whether a tuple's allele names a structural allele: Ensembl VEP writes a sequence allele
/// as its bases or `-` and a symbolic structural allele as its Sequence Ontology class, so
/// anything outside `-` and the nucleotide alphabet is structural.
fn is_structural_allele(allele: &str) -> bool {
    allele != "-"
        && !allele
            .to_uppercase()
            .chars()
            .all(|letter| "ACGTN".contains(letter))
}

/// Ensembl VEP's emission of `start_lost` and `start_retained_variant` together on one
/// allele (Utils/VariationEffect.pm:851-881, :947-963 and :1028-1075): the reference's set
/// holds both start terms and the candidate's is the same set without the erroneous member
/// for the allele kind, `start_retained_variant` on a structural allele and `start_lost` on
/// a sequence allele.
fn start_cooccurrence_swap(expected: &str, actual: &str, allele: &str) -> bool {
    let (reference, candidate) = (joined_terms(expected), joined_terms(actual));
    if !(reference.contains("start_lost") && reference.contains("start_retained_variant")) {
        return false;
    }
    if !candidate.is_subset(&reference) {
        return false;
    }
    let dropped = if is_structural_allele(allele) {
        "start_retained_variant"
    } else {
        "start_lost"
    };
    reference.difference(&candidate).eq([&dropped])
}

/// The one shape Ensembl VEP's `_intron_effects` (BaseTranscriptVariationAllele.pm:99-224)
/// can produce by the assignment at :215: the candidate's set is the reference's with
/// `splice_region_variant` added, with `intron_variant` or alone.
fn covered_splice_region_swap(expected: &str, actual: &str) -> bool {
    let (reference, candidate) = (joined_terms(expected), joined_terms(actual));
    if !reference.is_subset(&candidate) {
        return false;
    }
    let mut added = candidate.difference(&reference);
    let region = added.clone().any(|term| *term == SPLICE_REGION_TERM);
    region && added.all(|term| [SPLICE_REGION_TERM, INTRON_TERM].contains(term))
}

/// The registry class among `classes` whose excluding rule accepts the pair of consequence
/// sets, each a sorted comma-joined text (`consequence_set`), or `None`: the rules by which
/// the consequence comparison's adjusted figure sets a pair aside, `start_cooccurrence_swap`
/// first.
pub fn registry_class(
    expected_consequence: &str,
    actual_consequence: &str,
    allele: &str,
    classes: &[&str],
) -> Option<&'static str> {
    if classes.contains(&START_COOCCURRENCE_SWAP)
        && start_cooccurrence_swap(expected_consequence, actual_consequence, allele)
    {
        return Some(START_COOCCURRENCE_SWAP);
    }
    if classes.contains(&SPLICE_FAMILY_SWAP)
        && covered_splice_region_swap(expected_consequence, actual_consequence)
    {
        return Some(SPLICE_FAMILY_SWAP);
    }
    None
}

/// Whether the entry is a transcript's row.
fn transcript_row(e: &Entry) -> bool {
    feature_kind(entry_kind(e)) == "transcript"
}

fn spelt_in(allele: &str, letters: &[u8]) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| letters.contains(&b))
}

/// Whether the entry carries the start and the end of a transcript position: the two bounds
/// of the text formats' `column` (`?` for an end that maps to no transcript base, one
/// number standing for both, `-` for neither, each read up to the `/total` suffix
/// `--total_length` appends, OutputFactory.pm:1735-1748), or the keys `start_key` and
/// `end_key`, which an unmapped end leaves out.
fn position_bounds(e: &Entry, column: &str, start_key: &str, end_key: &str) -> (bool, bool) {
    let Some(position) = e.fields.get(column) else {
        return (
            !blank(e.fields.get_or(start_key, ABSENT)),
            !blank(e.fields.get_or(end_key, ABSENT)),
        );
    };
    let position = pytext::partition(position, "/").0;
    if blank(position) {
        return (false, false);
    }
    let (start, sep, end) = pytext::partition(position, "-");
    (
        start != "?",
        (if sep.is_empty() { start } else { end }) != "?",
    )
}

/// Whether the entry carries the coordinates the stop codon test edits at
/// (Utils/VariationEffect.pm:1507-1508): both cDNA positions and the CDS start. A CDS end
/// alone, the `?-N` of a span that begins in the 5' UTR, is not one of them.
pub fn stop_codon_test_coordinates(e: &Entry) -> bool {
    position_bounds(e, "cDNA_position", "cdna_start", "cdna_end") == (true, true)
        && position_bounds(e, "CDS_position", "cds_start", "cds_end").0
}

/// Whether the entry shows that `_ins_del_stop_altered_cil` (Utils/VariationEffect.pm
/// :1488-1539) returns before it edits anything: the allele holds a letter other than `A`,
/// `C`, `G` and `T` (:1498), the span lacks a cDNA position at one end, or it has no CDS
/// start (:1508).
pub fn stop_codon_test_not_made(e: &Entry) -> bool {
    (e.allele != "-" && !spelt_in(&e.allele, b"ACGT")) || !stop_codon_test_coordinates(e)
}

/// `stop_codon_overlap_genomic_span`: a Transcript entry of an insertion or deletion with
/// both cDNA positions and a CDS start on both sides, `splice_region_variant` in the
/// reference's terms, and the two sets equal but for the reference's
/// `coding_sequence_variant` against the candidate's `stop_retained_variant`
/// (Utils/VariationEffect.pm:1389-1431).
pub fn stop_codon_overlap_genomic_span_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) || !indel_row(e)? {
        return Ok(false);
    }
    let (expected, gone, come) = term_differences(e, a);
    Ok(is(&gone, &[CODING_SEQUENCE_TERM])
        && is(&come, &[STOP_RETAINED_TERM])
        && expected.contains(SPLICE_REGION_TERM)
        && !stop_codon_test_not_made(e)
        && !stop_codon_test_not_made(a))
}

/// Whether an entry is an insertion's or a deletion's as `indel_row` reads one, an allele
/// of other nucleotide letters than `A`, `C`, `G`, `T` and `N` included: where the record's
/// own alleles say nothing, such an allele against a Location of another length.
fn length_changing_row(e: &Entry) -> Res<bool> {
    if indel_row(e)? {
        return Ok(true);
    }
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    if e.indel.is_some()
        || e.allele.is_empty()
        || !e.allele.bytes().all(|b| b.is_ascii_alphabetic())
    {
        return Ok(false);
    }
    let (start, end) = location_ends(location);
    if !(pytext::is_ascii_digits(start) && pytext::is_ascii_digits(end)) {
        return Ok(false);
    }
    let between = integer_of(end, "Location position")? - integer_of(start, "Location position")?;
    Ok(between.saturating_add(1) != e.allele.len() as i128)
}

/// `stop_retained_without_stop_codon_test`: a Transcript entry of an insertion or deletion,
/// the reference's terms holding `stop_retained_variant` and the candidate's being the same
/// without it, with nothing, `coding_sequence_variant` or `frameshift_variant` in its
/// place, and both entries showing that the stop codon test is not made
/// (Utils/VariationEffect.pm:1318, :1488-1539).
pub fn stop_retained_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) || !length_changing_row(e)? {
        return Ok(false);
    }
    let (_, gone, come) = term_differences(e, a);
    Ok(is(&gone, &[STOP_RETAINED_TERM])
        && (come.is_empty() || is(&come, &[CODING_SEQUENCE_TERM]) || is(&come, &[FRAMESHIFT_TERM]))
        && stop_codon_test_not_made(e)
        && stop_codon_test_not_made(a))
}

/// The terms `frameshift` (Utils/VariationEffect.pm:1541-1567) reads away where
/// `stop_retained` does not hold.
const FRAMESHIFT_READ_TERMS: [&str; 3] = [
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
];

/// `stop_retained_three_base_stop_codon_test`: a Transcript entry of an insertion or
/// deletion whose stop codon test is made on both sides, the reference's terms holding
/// `stop_retained_variant` and the candidate's `coding_sequence_variant` in its place, or,
/// where the reference's hold it alone or beside one of `inframe_insertion`,
/// `inframe_deletion` and `protein_altering_variant`, `frameshift_variant` with or without
/// `stop_gained` in the place of those terms (Utils/VariationEffect.pm:1488-1539).
pub fn stop_retained_three_base_stop_codon_test_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) || !length_changing_row(e)? {
        return Ok(false);
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
    Ok(terms_stand && !stop_codon_test_not_made(e) && !stop_codon_test_not_made(a))
}

/// `stop_codon_removed_without_stop_codon_test`: a Transcript entry of an insertion or
/// deletion, both entries lacking a cDNA position at one end or a CDS start, and the two
/// sets of terms equal but for the reference's `stop_retained_variant` against the
/// candidate's `stop_lost` (Utils/VariationEffect.pm:1234-1288, :1508).
pub fn stop_codon_removed_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) || !length_changing_row(e)? {
        return Ok(false);
    }
    let (_, gone, come) = term_differences(e, a);
    Ok(is(&gone, &[STOP_RETAINED_TERM])
        && is(&come, &[STOP_LOST_TERM])
        && !stop_codon_test_coordinates(e)
        && !stop_codon_test_coordinates(a))
}

/// Whether an entry `insertion_row` takes is an insertion between two neighbouring
/// positions, the allele Ensembl VEP stores with a start one above its end and the only one
/// whose start `_overlaps_stop_codon_cil` moves (Utils/VariationEffect.pm:1408-1415). Its
/// Location, where the entry has one, is two neighbouring positions, and the first of these
/// that applies decides:
///
/// - its reference bases (`Entry.reference`) are empty: a pure insertion's entry, whatever
///   the allele's length and however the record pads its alleles;
/// - the parser read the record's own alleles: its reference allele is empty once the
///   shared first base is dropped (`Entry.insertion`);
/// - the lines of the input give the entry bases that are not empty: a replacement's entry;
/// - no line gives the entry its bases: it has a Location and its allele has another length
///   than two, which refuses the two-base substitution of the same Location.
///
/// Under an identifier without the lines of the input, two bases replaced by a sequence of
/// another length are therefore taken with the insertions and an insertion of two bases is
/// refused. The lines that carry the entry's identifier give its bases by Location and
/// allele: an insertion none of them places is refused where one of them replaces the same
/// two bases by the same allele.
fn inserted_between_neighbours(e: &Entry) -> bool {
    if let Some(location) = e.location.as_deref() {
        let (start, end) = location_ends(location);
        if !(pytext::is_ascii_digits(start) && pytext::is_ascii_digits(end)) {
            return false;
        }
        match (pytext::int(start), pytext::int(end)) {
            (Some(start), Some(end)) if start.checked_add(1) == Some(end) => {}
            _ => return false,
        }
    }
    match (e.reference.borrow().as_deref(), e.insertion) {
        (Some(""), _) => true,
        (_, Some(insertion)) => insertion,
        (Some(_), None) => false,
        (None, None) => e.location.is_some() && e.allele.chars().count() != 2,
    }
}

/// `stop_codon_overlap_insertion_range`: a Transcript entry of an insertion between two
/// neighbouring positions (`inserted_between_neighbours`) with both cDNA positions and a
/// CDS start on both sides, the two sets of terms equal but for the reference's
/// `coding_sequence_variant` against the candidate's `stop_retained_variant`, or the
/// reference's `frameshift_variant` with or without `stop_gained` against the candidate's
/// `stop_retained_variant` with or without `inframe_insertion`
/// (Utils/VariationEffect.pm:1389-1431). Beside `splice_region_variant` the first pair is
/// `stop_codon_overlap_genomic_span_pair`'s where the allele is of `A`, `C`, `G` and `T`
/// alone, and this rule's where it holds another letter, the rows that one refuses: no
/// peptide is read for such an allele (TranscriptVariationAllele.pm:701), so
/// `stop_retained` takes its second arm (Utils/VariationEffect.pm:1318) and the moved start
/// fails the overlap there too.
pub fn stop_codon_overlap_insertion_range_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) || !insertion_row(e) || !inserted_between_neighbours(e) {
        return Ok(false);
    }
    let (expected, gone, come) = term_differences(e, a);
    let terms_stand = if is(&gone, &[CODING_SEQUENCE_TERM]) {
        let unknown_base = [e, a]
            .iter()
            .any(|row| row.allele != "-" && !spelt_in(&row.allele, b"ACGT"));
        is(&come, &[STOP_RETAINED_TERM]) && (!expected.contains(SPLICE_REGION_TERM) || unknown_base)
    } else if is(&gone, &[FRAMESHIFT_TERM]) || is(&gone, &[STOP_GAINED_TERM, FRAMESHIFT_TERM]) {
        is(&come, &[STOP_RETAINED_TERM]) || is(&come, &[INFRAME_INSERTION_TERM, STOP_RETAINED_TERM])
    } else {
        false
    };
    Ok(terms_stand && stop_codon_test_coordinates(e) && stop_codon_test_coordinates(a))
}

/// A set of the seven terms a row of `intron_rows` holds, one bit a term in the order of
/// `INTRON_ROW_TERMS`.
type TermSet = u8;

/// The five splice terms whose windows lie inside an intron: the two bases at its 5' end,
/// the two at its 3' end, its fifth base, its third to sixth base, and the bases 3 to 17
/// from its 3' end (Utils/VariationEffect.pm:558-614); the intron term, whose flag
/// `within_intron` (:629-634) returns and whose window is the intron itself, its end bases
/// included; and `splice_region_variant`.
const INTRON_ROW_TERMS: [&str; 7] = [
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    INTRON_TERM,
    SPLICE_REGION_TERM,
];
const DONOR: TermSet = 1;
const ACCEPTOR: TermSet = 1 << 1;
const DONOR_5TH_BASE: TermSet = 1 << 2;
const DONOR_REGION: TermSet = 1 << 3;
const TRACT: TermSet = 1 << 4;
const INTRON: TermSet = 1 << 5;
const REGION: TermSet = 1 << 6;
/// The six terms `_intron_effects` sets from the differing regions and the candidate writes
/// only where the record's span holds a base of the term's window.
const PAST_THE_SPAN_TERMS: TermSet =
    DONOR | ACCEPTOR | DONOR_5TH_BASE | DONOR_REGION | TRACT | INTRON;
/// The four beside which `splice_region` returns 0 (Utils/VariationEffect.pm:616-627).
const SPLICE_REGION_WITHHOLDING_TERMS: TermSet = DONOR | ACCEPTOR | DONOR_5TH_BASE | DONOR_REGION;
/// A gap past which no window of an intron's end reaches a span or its regions.
const FAR: i64 = 1_000_000_000;
/// The names `--terms display` prints for the seven terms, each with the terms that print it
/// (the `display_term` of each term, Utils/Constants.pm): one name for the donor and the
/// acceptor term, one for the fifth base, the donor region, the tract and
/// `splice_region_variant`, one for the intron term. No other term prints one of the three.
const DISPLAY_ROW_TERMS: [(&str, TermSet); 3] = [
    ("ESSENTIAL_SPLICE_SITE", DONOR | ACCEPTOR),
    (
        "SPLICE_SITE",
        DONOR_5TH_BASE | DONOR_REGION | TRACT | REGION,
    ),
    ("INTRONIC", INTRON),
];
/// The row rules of the stop codon classes, whose terms a row can differ by beside the six:
/// one row can carry both differences.
const STOP_CODON_ROW_PAIRS: [fn(&Entry, &Entry) -> Res<bool>; 5] = [
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

/// The display names of the seven among a row's terms, one bit a name in the order of
/// [`DISPLAY_ROW_TERMS`].
fn display_row_names(terms: &BTreeSet<&str>) -> u8 {
    DISPLAY_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, (name, _))| terms.contains(name))
        .fold(0, |set, (bit, _)| set | 1 << bit)
}

/// The display names a row of the seven prints, as [`display_row_names`] reads them.
fn displayed(row: TermSet) -> u8 {
    DISPLAY_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, (_, members))| row & members != 0)
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

/// Whether a part of a position is digits of `0` to `9`, `?` or `-`.
fn position_part(part: &str) -> bool {
    matches!(part, "?" | "-") || pytext::is_ascii_digits(part)
}

/// The cDNA start and end and the CDS start and end an entry carries, as written (`?` for
/// an end that maps to no transcript base, `-` for a position the entry lacks), or `None`
/// where a position is not of that form. A text format's column is read up to the `/total`
/// suffix `--total_length` appends, one number standing for both ends; the keys of an
/// entry without the column are read as they are.
fn transcript_positions(e: &Entry) -> Option<[&str; 4]> {
    let mut positions = [""; 4];
    let columns = [
        ("cDNA_position", "cdna_start", "cdna_end"),
        ("CDS_position", "cds_start", "cds_end"),
    ];
    for (at, (column, start_key, end_key)) in columns.into_iter().enumerate() {
        let (start, end) = match e.fields.get(column) {
            None => (
                dash(e.fields.get_or(start_key, "-")),
                dash(e.fields.get_or(end_key, "-")),
            ),
            Some("" | "-") => ("-", "-"),
            Some(position) => {
                let (start, sep, end) = pytext::partition(pytext::partition(position, "/").0, "-");
                (start, if sep.is_empty() { start } else { end })
            }
        };
        if !(position_part(start) && position_part(end)) {
            return None;
        }
        (positions[2 * at], positions[2 * at + 1]) = (start, end);
    }
    Some(positions)
}

/// The number of bases of the transcript's sequence from the first cDNA position of a row
/// to its second, both counted; a number past the largest `i128` is read as that one.
fn sequence_width(positions: &[&str; 4]) -> Res<i128> {
    let last = integer_of(positions[1], "cDNA position")?;
    Ok((last - integer_of(positions[0], "cDNA position")?).saturating_add(1))
}

/// Whether an entry shows, without its record, a span of three or more bases that lies in
/// the transcript's sequence from end to end under an allele of more nucleotide letters
/// than the span has bases: its Location names the span, and its two cDNA positions are
/// numbers as far apart as the span's two ends, which a span that holds an intron's base
/// does not have. A Location of two neighbouring positions is refused: an insertion between
/// them is written with the same Location, allele and positions as two replaced bases, and
/// only the record tells the two apart.
fn longer_allele_inside_one_exon(e: &Entry) -> Res<bool> {
    let Some(span) = location_span(e)? else {
        return Ok(false);
    };
    if span < 3 || !spelt_in(&e.allele, b"ACGTN") || e.allele.len() as i128 <= span {
        return Ok(false);
    }
    let Some(positions) = transcript_positions(e) else {
        return Ok(false);
    };
    if !(pytext::is_ascii_digits(positions[0]) && pytext::is_ascii_digits(positions[1])) {
        return Ok(false);
    }
    Ok(sequence_width(&positions)? == span)
}

/// The regions `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460) cuts a
/// record into, as offsets from the span's first base: where both alleles hold more than
/// one base (:418) the runs of indices at which they differ, every index past the shorter
/// allele's end differing (:428-450); else one region, the span (:452-454).
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

/// The splice and intron terms of the reference's row and of the candidate's for a record
/// whose span is the offsets `0` to `span - 1`, on a transcript of the given strand whose
/// one intron near the span runs from offset `start` to offset `end`.
///
/// The reference's row is `_intron_effects` (BaseTranscriptVariationAllele.pm:99-224) on
/// `regions`, each placed with no bound at the span's end (:117-118). The intron is read
/// where the span overlaps it widened by three bases, and its two ends where the span
/// overlaps the eleven bases around each (BaseTranscriptVariation.pm:1027-1029); an intron
/// of thirteen bases or fewer is a frameshift intron (:1032-1035) and a region that
/// overlaps it sets nothing (BaseTranscriptVariationAllele.pm:138-141, :171-174). A region
/// in the intron's interior sets the intron term (:143-150) and one in a window of the five
/// splice terms its flag (:152-161, :176-206). The flag of `splice_region` is assigned for
/// each region in turn, whether the region lies in one of its four windows, until a region
/// has met the intron's first two or last two bases (:209-216; Utils/VariationEffect.pm
/// :91-94), so the last region assigned decides it. The tract term is written only off the
/// exons (Constants.pm:971-983), which a transcript with a frameshift intron reads twelve
/// bases wider (BaseTranscriptVariation.pm:861; `wide_exons` says another intron of the
/// transcript is one). The donor region is withheld beside the fifth base and
/// `splice_region_variant` beside the first four terms (Utils/VariationEffect.pm:597-627).
///
/// The candidate's row keeps a term of the six only where the span itself holds a base of a
/// window of that term, the intron being the intron term's window, reads `splice_region`
/// from any region in one of its four windows, and withholds the donor region and
/// `splice_region_variant` by its own row's terms. For an allele that opens with its whole
/// reference allele (`insertion_after_the_span`) it is the row of the insertion the allele
/// writes (`insertion_after_the_span_row`).
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
    let candidate = if insertion_after_the_span(regions, span) {
        insertion_after_the_span_row(span, start, end, forward, wide_exons)
    } else {
        row(set_by_a_region & held_by_the_span, in_any)
    };
    (row(set_by_a_region, last_assigned), candidate)
}

/// Whether `regions` are those of an allele that opens with its whole reference allele of
/// `span` bases, two or more, and holds further bases: one run that begins at the offset
/// `span`, no base of the span differing (VariationFeatureOverlapAllele.pm:418, :428-450).
/// Such an allele changes no base of its span, and the candidate reads it as the insertion it
/// writes, between the span's last base and the base after it.
fn insertion_after_the_span(regions: &[(i64, i64)], span: i64) -> bool {
    span >= 2 && matches!(regions, [(first, _)] if *first == span)
}

/// Whether `allele` opens with its reference allele `reference` of one base and holds
/// further bases. The reference does not compare such a pair of alleles and reads the one
/// base as changed (VariationFeatureOverlapAllele.pm:418, :452-454), so the regions do not
/// show the shape. Such an allele changes no base of its span either, and the candidate
/// reads it as the insertion it writes, between that base and the base after it.
fn insertion_after_one_base(reference: &[char], allele: &[char]) -> bool {
    matches!(reference, [base] if allele.len() > 1 && allele[0] == *base)
}

/// The splice and intron terms of the candidate's row for an allele
/// `insertion_after_the_span` takes, on the intron `intron_rows` reads: those of the
/// insertion between the span's last base, the offset `span - 1`, and the base after it, read
/// by those two neighbours as `_intron_effects` reads an insertion
/// (BaseTranscriptVariationAllele.pm:126).
///
/// The intron is read where the two neighbours overlap it widened by three bases, and its two
/// ends where they overlap the eleven bases around each (BaseTranscriptVariation.pm:1027-1029);
/// a frameshift intron that holds both sets nothing (BaseTranscriptVariationAllele.pm:138-141,
/// :171-174). A window of the donor, acceptor or donor region term that holds both neighbours
/// sets its flag (:176-206), which the one base of the fifth base's window never does. The
/// tract's window sets its flag where it holds either neighbour (:152-161), and the term is
/// written off the exons the two neighbours overlap, which a transcript with a frameshift
/// intron reads twelve bases wider. The intron term stands where the intron's interior holds
/// both neighbours, or the second is its third base, or the first its third base from the end
/// (:143-150). `splice_region` stands where one of its four windows holds both neighbours, or
/// the second neighbour is the intron's first or third base, or the first its last base or its
/// third from the end (Utils/VariationEffect.pm:91-102), and is withheld beside the donor,
/// acceptor and donor region terms.
fn insertion_after_the_span_row(
    span: i64,
    start: i64,
    end: i64,
    forward: bool,
    wide_exons: bool,
) -> TermSet {
    let (before, after) = (span - 1, span);
    let in_list = after >= start - 3 && before <= end + 3;
    let at_an_end =
        (after >= start - 3 && before <= start + 7) || (after >= end - 7 && before <= end + 3);
    let frameshift = end - start <= 12;
    let holds_both = |first: i64, final_base: i64| before >= first && after <= final_base;
    if !(in_list || at_an_end) || (frameshift && holds_both(start, end)) {
        return 0;
    }
    let (windows, tract) = if forward {
        (
            [
                (DONOR, start, start + 1),
                (ACCEPTOR, end - 1, end),
                (DONOR_REGION, start + 2, start + 5),
            ],
            (end - 16, end - 2),
        )
    } else {
        (
            [
                (DONOR, end - 1, end),
                (ACCEPTOR, start, start + 1),
                (DONOR_REGION, end - 5, end - 2),
            ],
            (start + 2, start + 16),
        )
    };
    let mut terms = 0;
    for (term, first, final_base) in windows {
        if at_an_end && holds_both(first, final_base) {
            terms |= term;
        }
    }
    let on_an_exon = before < start
        || after > end
        || ((frameshift || wide_exons) && (before < start + 12 || after + 12 > end));
    if after >= tract.0 && before <= tract.1 && !on_an_exon {
        terms |= TRACT;
    }
    if in_list && (holds_both(start + 2, end - 2) || after == start + 2 || before == end - 2) {
        terms |= INTRON;
    }
    let beside = [
        (start + 2, start + 7),
        (end - 7, end - 2),
        (start - 3, start - 1),
        (end + 1, end + 3),
    ];
    let region = beside
        .iter()
        .any(|&(first, final_base)| holds_both(first, final_base))
        || after == start
        || before == end
        || after == start + 2
        || before == end - 2;
    if at_an_end && region && terms & SPLICE_REGION_WITHHOLDING_TERMS == 0 {
        terms |= REGION;
    }
    terms
}

/// Calls `visit` with the places `(forward, start, end)` of one intron beside a span of
/// `span` bases that a row's cDNA positions allow, as offsets from the span's first base,
/// for regions that run `reach` bases past the span's end, until it returns true, and
/// returns whether it did. A numbered cDNA position is a base of the transcript's sequence
/// and any other a base of the intron; the first position is the span's first base on the
/// forward strand and its last on the reverse. Two numbers as far apart as the span put the
/// intron after the span, one to seven bases on, as far as one of its ends is still read,
/// or before it, where an intron shorter than its windows has them run on past its last
/// base; two numbers closer than the span put it inside the span. A gap wider than any
/// window, and than the regions' reach, is placed once.
fn intron_placements(
    span: i64,
    reach: i64,
    positions: &[&str; 4],
    mut visit: impl FnMut(bool, i64, i64) -> bool,
) -> Res<bool> {
    let last = span - 1;
    let before: Vec<i64> = (0..18).map(|gap| -gap).chain([-FAR]).collect();
    let gaps: BTreeSet<i64> = (0..42).chain(42.max(reach - 20)..reach + 19).collect();
    let after: Vec<i64> = gaps.iter().map(|gap| last + gap).chain([FAR]).collect();
    let numbered_first = pytext::is_ascii_digits(positions[0]);
    let numbered_last = pytext::is_ascii_digits(positions[1]);
    let width = if numbered_first && numbered_last {
        Some(sequence_width(positions)?)
    } else {
        None
    };
    for forward in [true, false] {
        let (first_in_sequence, last_in_sequence) = match forward {
            true => (numbered_first, numbered_last),
            false => (numbered_last, numbered_first),
        };
        if let Some(width) = width {
            if width == i128::from(span) {
                for start in (1..=7).map(|gap| last + gap) {
                    let far_ends = after.iter().copied().filter(|end| *end >= start);
                    let near_ends = (0..42).map(|length| start + length);
                    if far_ends
                        .chain(near_ends)
                        .any(|end| visit(forward, start, end))
                    {
                        return Ok(true);
                    }
                }
                if (-16..0).any(|start| (start..0).any(|end| visit(forward, start, end))) {
                    return Ok(true);
                }
            } else if 0 < width && width < i128::from(span) {
                let length = span - width as i64;
                if (1..=last - length).any(|start| visit(forward, start, start + length - 1)) {
                    return Ok(true);
                }
            }
        } else if first_in_sequence {
            if (1..span).any(|start| after.iter().any(|end| visit(forward, start, *end))) {
                return Ok(true);
            }
        } else if last_in_sequence {
            if (0..last).any(|end| before.iter().any(|start| visit(forward, *start, end))) {
                return Ok(true);
            }
        } else if before
            .iter()
            .any(|start| after.iter().any(|end| visit(forward, *start, *end)))
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// A place of the transcript's intron, `(forward, start, end, wide_exons)` as `intron_rows`
/// reads it, at which the record's regions give the reference row's splice and intron terms
/// `expected` and the span, or the insertion written after it, gives the candidate's
/// `actual`, or `None` where the cDNA positions `positions` allow no such place
/// (`intron_placements`). The exons are read at their own width or, for a tract term a row
/// lacks, twelve bases wider. `reference` is the reference bases the allele stands in for.
/// For an allele that opens with its one reference base (`insertion_after_one_base`) the
/// reference's row is that of the base read as changed, and the candidate's the row of the
/// insertion after it (`insertion_after_the_span_row`).
fn record_span_placement(
    positions: &[&str; 4],
    reference: &[char],
    allele: &[char],
    expected: TermSet,
    actual: TermSet,
) -> Res<Option<(bool, i64, i64, bool)>> {
    let regions = differing_regions(reference, allele);
    let span = reference.len() as i64;
    let after_one_base = insertion_after_one_base(reference, allele);
    let rows = (expected, actual);
    let mut place = None;
    let visit = |forward, start, end| {
        let read = |wide_exons: bool| {
            placed_rows(
                &regions,
                span,
                after_one_base,
                (forward, start, end, wide_exons),
            )
        };
        let narrow = read(false);
        let wide_exons = if narrow == rows {
            false
        } else if ((narrow.0 & TRACT != 0 && expected & TRACT == 0)
            || (narrow.1 & TRACT != 0 && actual & TRACT == 0))
            && read(true) == rows
        {
            true
        } else {
            return false;
        };
        place = Some((forward, start, end, wide_exons));
        true
    };
    intron_placements(span, allele.len() as i64 - span, positions, visit)?;
    Ok(place)
}

/// The reference's row and the candidate's at one place `(forward, start, end, wide_exons)`
/// of the intron: `intron_rows`, with the candidate's row that of the insertion after the base
/// (`insertion_after_the_span_row`) for an allele that opens with its one reference base.
fn placed_rows(
    regions: &[(i64, i64)],
    span: i64,
    after_one_base: bool,
    (forward, start, end, wide_exons): (bool, i64, i64, bool),
) -> (TermSet, TermSet) {
    let (reference_row, span_row) = intron_rows(regions, span, start, end, forward, wide_exons);
    if after_one_base {
        let insertion = insertion_after_the_span_row(span, start, end, forward, wide_exons);
        (reference_row, insertion)
    } else {
        (reference_row, span_row)
    }
}

/// Whether the cDNA positions `positions` allow a place of the transcript's intron at which
/// the reference's splice and intron terms print the display names `expected` and the
/// candidate's the names `actual` (`displayed`). The exons are read twelve bases wider as
/// well where a row at their own width holds the tract term, whose name three other terms
/// print.
fn displayed_span_placement(
    positions: &[&str; 4],
    reference: &[char],
    allele: &[char],
    expected: u8,
    actual: u8,
) -> Res<bool> {
    let regions = differing_regions(reference, allele);
    let span = reference.len() as i64;
    let after_one_base = insertion_after_one_base(reference, allele);
    let shown = |(reference_row, candidate_row): (TermSet, TermSet)| {
        (displayed(reference_row), displayed(candidate_row)) == (expected, actual)
    };
    let visit = |forward, start, end| {
        let read = |wide_exons: bool| {
            placed_rows(
                &regions,
                span,
                after_one_base,
                (forward, start, end, wide_exons),
            )
        };
        let narrow = read(false);
        shown(narrow) || ((narrow.0 | narrow.1) & TRACT != 0 && shown(read(true)))
    };
    intron_placements(span, allele.len() as i64 - span, positions, visit)
}

/// The class on two rows printed under `--terms display`, or `None` where neither row holds
/// a display name of the seven (`DISPLAY_ROW_TERMS`). The rows' other terms are equal and
/// their names of the seven differ; the candidate's row alone does not hold the name of the
/// donor and the acceptor term, as it alone holds neither term under the Sequence Ontology
/// names; the record's alleles show a longer alternate allele that is no pure insertion, as
/// under those names, and a place of one intron gives both rows' names
/// (`displayed_span_placement`). Where the allele is not written after its span, the
/// reference's row holds a name the candidate's lacks, as it holds a term of the six alone
/// under those names. Without the record's alleles no such row is taken: three names stand
/// for seven terms, and the rows alone do not tell them apart.
fn display_vocabulary_rows(
    e: &Entry,
    a: &Entry,
    positions: &[&str; 4],
    expected_terms: &BTreeSet<&str>,
    actual_terms: &BTreeSet<&str>,
) -> Res<Option<bool>> {
    let expected = display_row_names(expected_terms);
    let actual = display_row_names(actual_terms);
    if expected | actual == 0 {
        return Ok(None);
    }
    let others = |terms: &BTreeSet<&str>| -> Vec<String> {
        let named = |term: &&str| DISPLAY_ROW_TERMS.iter().any(|(name, _)| name == term);
        terms
            .iter()
            .filter(|term| !named(term))
            .map(|term| term.to_string())
            .collect()
    };
    let site_name = displayed(DONOR | ACCEPTOR);
    if expected == actual
        || actual & !expected & site_name != 0
        || others(expected_terms) != others(actual_terms)
    {
        return Ok(Some(false));
    }
    let reference = e
        .reference
        .borrow()
        .clone()
        .or_else(|| a.reference.borrow().clone());
    let Some(reference) = reference else {
        return Ok(Some(false));
    };
    let reference: Vec<char> = reference.chars().collect();
    let allele: Vec<char> = e.allele.chars().collect();
    let after_one_base = insertion_after_one_base(&reference, &allele);
    if (reference.len() < 2 && !after_one_base)
        || location_span(e)? != Some(reference.len() as i128)
        || !spelt_in(&e.allele, b"ACGTN")
        || allele.len() <= reference.len()
    {
        return Ok(Some(false));
    }
    let after_the_span = after_one_base
        || insertion_after_the_span(
            &differing_regions(&reference, &allele),
            reference.len() as i64,
        );
    if !after_the_span && expected & !actual == 0 {
        return Ok(Some(false));
    }
    Ok(Some(displayed_span_placement(
        positions, &reference, &allele, expected, actual,
    )?))
}

/// `splice_term_past_the_record_span` on two entries at their Locations: a Transcript entry
/// with the same cDNA and CDS positions on both sides whose two sets of terms differ by
/// terms of the six the reference alone holds and, at most, the donor region and
/// `splice_region_variant` the candidate alone holds, beside the terms of a stop codon class
/// where its rule takes the rest of the difference (`STOP_CODON_ROW_PAIRS`).
///
/// `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460) compares the two
/// alleles index by index to the longer allele's length (:428-429), so every index past the
/// reference allele's end differs, and `_intron_effects` places each differing region on
/// the genome from the record's start with no bound at the span's end
/// (BaseTranscriptVariationAllele.pm:117-118): on a record whose alternate allele is the
/// longer one, the flags of the five splice terms and of the intron term (:143-206) are set
/// from positions the record does not replace. The candidate writes such a term only where
/// a base of the span itself lies in a window of that term.
///
/// An allele that opens with its whole reference allele of two or more bases changes no base
/// of its span: every index that differs lies past the reference allele's end (:428-429), and
/// the reference reads the further bases as a change of the bases after the span
/// (BaseTranscriptVariationAllele.pm:117-118). The candidate reads such an allele as the
/// insertion it writes, between the span's last base and the base after it, and writes that
/// insertion's terms (`insertion_after_the_span_row`): on its row the candidate alone may
/// hold the tract term and the intron term as well, and the reference need hold no term
/// alone.
///
/// An allele that opens with its reference allele of one base changes no base of its span
/// either. The reference does not compare the two alleles (:418) and reads the one base as
/// changed (:452-454), so its terms are those of that base, and the candidate's those of the
/// insertion after it (`insertion_after_one_base`): on its row the reference alone may hold
/// `splice_region_variant` as well.
///
/// With the record's alleles (`Entry.reference`) the two rows are those of a longer
/// alternate allele that is no pure insertion, both alleles holding more than one base or
/// the allele opening with its one reference base, at
/// a place of one intron the cDNA positions allow (`record_span_placement`). Without them
/// the span of three or more positions lies in the transcript's sequence under a longer
/// allele (`longer_allele_inside_one_exon`), the candidate's row holds none of the six, and
/// `splice_region_variant` alone comes, where a term that withholds it went; no row is then
/// read as that of an allele that opens with its whole reference allele, which the record
/// alone shows.
///
/// Two rows printed under `--terms display` are read by the names the seven terms print
/// there, with the record's alleles alone (`display_vocabulary_rows`).
fn splice_term_past_the_record_span_rows(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) {
        return Ok(false);
    }
    let Some(positions) = transcript_positions(e) else {
        return Ok(false);
    };
    if transcript_positions(a) != Some(positions) {
        return Ok(false);
    }
    let (expected_terms, actual_terms) = (terms(&e.consequence), terms(&a.consequence));
    if let Some(taken) = display_vocabulary_rows(e, a, &positions, &expected_terms, &actual_terms)?
    {
        return Ok(taken);
    }
    let expected = intron_row_terms(&expected_terms);
    let actual = intron_row_terms(&actual_terms);
    let other = other_terms(&expected_terms);
    if other != other_terms(&actual_terms) {
        // The rest of the difference is a stop codon class's: its rule reads the two rows
        // with the candidate's splice and intron terms on both, and this one the two rows
        // with the reference's other terms on both. Those differ by terms of the seven
        // alone, the only terms of either row read below.
        let held = |bit: &usize| actual & (1 << bit) != 0;
        let splice_and_intron = (0..INTRON_ROW_TERMS.len())
            .filter(held)
            .map(|bit| INTRON_ROW_TERMS[bit]);
        let beside = Entry {
            consequence: other
                .iter()
                .copied()
                .chain(splice_and_intron)
                .collect::<Vec<_>>()
                .join(","),
            ..e.clone()
        };
        let mut taken = false;
        for pair in STOP_CODON_ROW_PAIRS {
            if pair(&beside, a)? {
                taken = true;
                break;
            }
        }
        if !taken {
            return Ok(false);
        }
    }
    let (gone, come) = (expected & !actual, actual & !expected);
    let reference = e
        .reference
        .borrow()
        .clone()
        .or_else(|| a.reference.borrow().clone());
    let reference: Option<Vec<char>> = reference.map(|reference| reference.chars().collect());
    let allele: Vec<char> = e.allele.chars().collect();
    let after_one_base = reference
        .as_ref()
        .is_some_and(|reference| insertion_after_one_base(reference, &allele));
    let after_the_span = after_one_base
        || reference.as_ref().is_some_and(|reference| {
            insertion_after_the_span(
                &differing_regions(reference, &allele),
                reference.len() as i64,
            )
        });
    let (may_come, differ) = if after_the_span {
        (REGION | DONOR_REGION | TRACT | INTRON, gone | come != 0)
    } else {
        (REGION | DONOR_REGION, gone != 0)
    };
    let may_go = if after_one_base {
        PAST_THE_SPAN_TERMS | REGION
    } else {
        PAST_THE_SPAN_TERMS
    };
    if !differ || gone & !may_go != 0 || come & !may_come != 0 {
        return Ok(false);
    }
    let Some(reference) = reference else {
        return Ok(longer_allele_inside_one_exon(e)?
            && longer_allele_inside_one_exon(a)?
            && actual & PAST_THE_SPAN_TERMS == 0
            && (come == 0 || (come == REGION && gone & SPLICE_REGION_WITHHOLDING_TERMS != 0)));
    };
    if location_span(e)? != Some(reference.len() as i128)
        || (reference.len() < 2 && !after_one_base)
    {
        return Ok(false);
    }
    if !spelt_in(&e.allele, b"ACGTN") || allele.len() <= reference.len() {
        return Ok(false);
    }
    Ok(record_span_placement(&positions, &reference, &allele, expected, actual)?.is_some())
}

/// `splice_term_past_the_record_span` on the two entries of one key, each read at the
/// Location its record gives its allele (`at_its_location`).
pub fn splice_term_past_the_record_span_pair(e: &Row, a: &Row) -> Res<bool> {
    splice_term_past_the_record_span_rows(&at_its_location(e), &at_its_location(a))
}

/// The IMPACT a row of these terms carries, its most severe term's (OutputFactory.pm:1337
/// and :1344 read the first consequence by rank), or `None` for a term without a rank and
/// for no term.
pub fn terms_impact(terms: &BTreeSet<&str>) -> Option<&'static str> {
    let rank = |term: &str| TERMS_BY_RANK.iter().position(|(ranked, _)| *ranked == term);
    let mut most_severe: Option<usize> = None;
    for term in terms {
        let rank = rank(term)?;
        most_severe = Some(most_severe.map_or(rank, |held| held.min(rank)));
    }
    most_severe.map(|rank| TERMS_BY_RANK[rank].1)
}

/// Whether the plugin that writes the field `name` writes it on a row of these terms, or
/// `None` for a field no term decides.
fn terms_write_field(name: &str, terms: &BTreeSet<&str>) -> Option<bool> {
    if field_named(name, &DOWNSTREAM_FIELDS) {
        return Some(
            terms.contains(FRAMESHIFT_TERM) && !terms.iter().any(|term| term.contains("splice")),
        );
    }
    if field_named(name, &TRUNCATING_TERM_FIELDS) {
        return Some(TRUNCATING_TERMS.iter().any(|term| terms.contains(term)));
    }
    if field_named(name, &MISSENSE_TERM_FIELDS) {
        return Some(terms.contains(MISSENSE_TERM));
    }
    if field_named(name, &PROTEIN_CHANGE_TERM_FIELDS) {
        return Some(PROTEIN_CHANGE_TERMS.iter().any(|term| terms.contains(term)));
    }
    None
}

/// Whether an HGVSp is a frameshift's notation: `fs`, or `ext` where the frameshift begins
/// at the stop codon.
fn hgvsp_frameshift_form(value: &str) -> bool {
    let change = pytext::partition(value, ":p.").2;
    change.contains("fs") || change.contains("ext")
}

/// Whether the differing values `ev` and `av` of the field `name` on the reference entry
/// `e` and the candidate entry `a` of one key are the ones the two entries' Consequence
/// sets give: the IMPACT of each set (its most severe term's); a field of a plugin that runs
/// on a row's terms (`terms_write_field`: Downstream, NMD, LOFTEE, REVEL, AlphaMissense,
/// ClinPred and dbNSFP under its prefix) present exactly on the side whose terms its plugin
/// writes on and absent on the other;
/// or an HGVSp that is a frameshift's notation exactly on the side without
/// `stop_retained_variant`, whose predicate keeps `frameshift` from holding
/// (Utils/VariationEffect.pm:1549, which the notation's type reads at
/// TranscriptVariationAllele.pm:2053), both naming one protein. A cell that differs in any
/// other way, and every field no term decides, is not taken.
pub fn term_pair_following_cell(e: &Entry, a: &Entry, name: &str, ev: &str, av: &str) -> bool {
    let (expected_terms, actual_terms) = (terms(&e.consequence), terms(&a.consequence));
    if expected_terms == actual_terms || ev == av {
        return false;
    }
    if field_named(name, &IMPACT_FIELDS) {
        return terms_impact(&expected_terms) == Some(ev)
            && terms_impact(&actual_terms) == Some(av);
    }
    if let Some(expected_written) = terms_write_field(name, &expected_terms) {
        let actual_written = terms_write_field(name, &actual_terms);
        return Some(expected_written) != actual_written
            && blank(ev) != expected_written
            && Some(blank(av)) != actual_written;
    }
    if field_named(name, &HGVSP_FIELDS) {
        let expected_kept = expected_terms.contains(STOP_RETAINED_TERM);
        let actual_kept = actual_terms.contains(STOP_RETAINED_TERM);
        if expected_kept == actual_kept
            || blank(ev)
            || blank(av)
            || pytext::partition(ev, ":p.").0 != pytext::partition(av, ":p.").0
        {
            return false;
        }
        return hgvsp_frameshift_form(ev) != expected_kept
            && hgvsp_frameshift_form(av) != actual_kept;
    }
    false
}

/// `splice_family_swap` on a Downstream field: the reference alone writes the field, on a
/// frameshift row none of whose terms names a splice site (Downstream.pm:95-99), and the
/// two Consequence sets are a pair the registry sets aside under the class.
pub fn splice_family_swap(c: &Consult, ev: &str, av: &str) -> Res<bool> {
    let Some(a) = c.a else {
        return Ok(false);
    };
    if blank(ev) || !blank(av) {
        return Ok(false);
    }
    let expected = consequence_set(&c.e.consequence);
    let mut terms = expected.split(',');
    if !terms.clone().any(|term| term == FRAMESHIFT_TERM)
        || terms.any(|term| term.contains("splice"))
    {
        return Ok(false);
    }
    let pair = registry_class(
        &expected,
        &consequence_set(&a.consequence),
        &c.e.allele,
        &[SPLICE_FAMILY_SWAP],
    );
    Ok(pair == Some(SPLICE_FAMILY_SWAP))
}

/// `plugin_value_follows_splice_term_past_the_record_span`: a field of Downstream, NMD or
/// LOFTEE that one engine alone writes on a row `splice_term_past_the_record_span_pair`
/// takes, where `term_pair_following_cell` reads the two values from the two sets of terms.
pub fn follows_splice_term_past_the_record_span(c: &Consult, ev: &str, av: &str) -> Res<bool> {
    let Some(a) = c.a else {
        return Ok(false);
    };
    if !splice_term_past_the_record_span_pair(c.e, a)? {
        return Ok(false);
    }
    Ok(term_pair_following_cell(c.e, a, c.name, ev, av))
}

/// The codons cell and the residues cell of the text formats and of a JSON entry.
const CODONS_FIELDS: [&str; 2] = ["Codons", "codons"];
const RESIDUE_FIELDS: [&str; 2] = ["Amino_acids", "amino_acids"];
/// The table id of the standard genetic code.
const STANDARD_CODE: u8 = 1;

/// The alternate peptide of the alternate half of a `Codons` cell that holds a letter
/// outside A, C, G and T after bases that settle a stop codon: a stop codon stands among
/// its whole codons before the first of them that holds such a letter. The peptide is the
/// translation of the whole codons under the standard genetic code, with the `X` both
/// engines append for bases past the last whole codon unless the peptide is the stop alone.
/// `None` for an alternate half of A, C, G and T alone, for one whose whole codons hold a
/// letter that names no base, and for one whose known bases settle no stop codon before
/// the first codon that holds another letter.
fn settled_alternate_peptide(alternate_codons: &str) -> Option<String> {
    let code = genetic_code(STANDARD_CODE)?;
    let bases = alternate_codons.to_ascii_uppercase();
    let bases = bases.as_bytes();
    let known = |base: &u8| b"ACGT".contains(base);
    if bases.iter().all(known) {
        return None;
    }
    let residue = |codon: &[u8]| {
        let letters = [codon[0], codon[1], codon[2]].map(char::from);
        nucleotide_codon_residue(code, letters)
    };
    let mut settled = false;
    for codon in bases.chunks_exact(3) {
        if !codon.iter().all(known) {
            break;
        }
        if residue(codon) == Some('*') {
            settled = true;
            break;
        }
    }
    if !settled {
        return None;
    }
    let mut peptide: String = bases.chunks_exact(3).map(residue).collect::<Option<_>>()?;
    if !bases.len().is_multiple_of(3) && peptide != "*" {
        peptide.push('X');
    }
    Some(peptide)
}

/// Whether two sets of terms differ as those of a stop codon the inserted bases keep: the
/// reference's terms (`expected`), those the reference alone holds (`gone`) and those
/// the candidate alone holds (`come`). The candidate alone holds `inframe_insertion`,
/// which needs both peptides (Utils/VariationEffect.pm:1113-1115) and cuts the alternate
/// one after its first stop (:1126). Where the reference holds `stop_retained_variant`,
/// from the second arm of `stop_retained` (:1318), whose edit test is not made for a letter
/// that names no base (:1498), nothing else differs. Otherwise the candidate alone holds
/// `stop_retained_variant` as well, the reference alone `frameshift_variant` (:1564) or
/// `coding_sequence_variant` (:1616-1649), and the entry is an insertion between two
/// neighbouring positions with the stop codon test's coordinates on both sides, whose moved
/// start fails the reference's overlap (:1408-1415).
fn kept_stop_terms(
    e: &Entry,
    a: &Entry,
    expected: &BTreeSet<&str>,
    gone: &BTreeSet<&str>,
    come: &BTreeSet<&str>,
) -> bool {
    if expected.contains(STOP_RETAINED_TERM) {
        return gone.is_empty() && is(come, &[INFRAME_INSERTION_TERM]);
    }
    (is(gone, &[FRAMESHIFT_TERM]) || is(gone, &[CODING_SEQUENCE_TERM]))
        && is(come, &[INFRAME_INSERTION_TERM, STOP_RETAINED_TERM])
        && insertion_row(e)
        && inserted_between_neighbours(e)
        && stop_codon_test_coordinates(e)
        && stop_codon_test_coordinates(a)
}

/// `stop_codon_settled_before_unknown_bases`: a Transcript entry whose `Codons` cell is the
/// same on both sides and has an alternate half `settled_alternate_peptide` reads; no
/// residues in the reference and, in the candidate, a reference half and that peptide; and
/// terms that differ by the gained stop alone or by the kept stop alone. For a gained stop,
/// where the reference writes `frameshift_variant` the candidate adds `stop_gained` to the
/// same terms; otherwise the
/// reference alone writes `coding_sequence_variant` and the candidate alone `stop_gained`,
/// with at most one of `inframe_insertion`, `inframe_deletion` and
/// `protein_altering_variant`. No peptide is read for an allele that holds a letter outside
/// A, C, G and T (TranscriptVariationAllele.pm:701), so `stop_gained`
/// (Utils/VariationEffect.pm:1211-1232) is false and `coding_unknown` (:1616-1649) holds
/// where the length is whole codons.
///
/// The reference can write `stop_retained_variant` in the place of either term: for an
/// insertion or deletion without a peptide `stop_retained` holds where the overlap with the
/// stop codon holds, read with an insertion's start moved by the inserted length, and the
/// stop codon test returns before it reads a base (Utils/VariationEffect.pm:1318,
/// :1408-1415, :1498). On such a row, an insertion's or a deletion's
/// (`length_changing_row`) with both entries showing that the test is not made
/// (`stop_codon_test_not_made`), the reference alone writes `stop_retained_variant` and the
/// candidate alone `stop_gained`, with `frameshift_variant`, one of the three in-frame terms
/// or no other term.
///
/// For a stop codon the inserted bases keep, the two sets differ as `kept_stop_terms` reads
/// them.
pub fn stop_codon_settled_before_unknown_bases_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) {
        return Ok(false);
    }
    let codons = e.fields.first_of(&CODONS_FIELDS, ABSENT);
    if blank(codons) || codons != a.fields.first_of(&CODONS_FIELDS, ABSENT) {
        return Ok(false);
    }
    let Some(peptide) = codons
        .split_once('/')
        .and_then(|(_, alternate)| settled_alternate_peptide(alternate))
    else {
        return Ok(false);
    };
    if !blank(e.fields.first_of(&RESIDUE_FIELDS, ABSENT)) {
        return Ok(false);
    }
    match a.fields.first_of(&RESIDUE_FIELDS, ABSENT).split_once('/') {
        Some((reference, alternate)) if !reference.is_empty() && alternate == peptide => {}
        _ => return Ok(false),
    }
    let (expected, gone, come) = term_differences(e, a);
    if kept_stop_terms(e, a, &expected, &gone, &come) {
        return Ok(true);
    }
    if expected.contains(FRAMESHIFT_TERM) {
        return Ok(gone.is_empty() && is(&come, &[STOP_GAINED_TERM]));
    }
    if !come.contains(STOP_GAINED_TERM) || come.len() > 2 {
        return Ok(false);
    }
    let beside_the_stop = come.iter().find(|term| **term != STOP_GAINED_TERM);
    let inframe = beside_the_stop.is_none_or(|term| FRAMESHIFT_READ_TERMS.contains(term));
    if is(&gone, &[STOP_RETAINED_TERM]) {
        return Ok(length_changing_row(e)?
            && stop_codon_test_not_made(e)
            && stop_codon_test_not_made(a)
            && (inframe || beside_the_stop == Some(&FRAMESHIFT_TERM)));
    }
    Ok(is(&gone, &[CODING_SEQUENCE_TERM]) && inframe)
}

const START_RETAINED_TERM: &str = "start_retained_variant";
const INFRAME_DELETION_TERM: &str = "inframe_deletion";

/// The cDNA end and the CDS end of an entry whose cDNA position and CDS position both
/// begin at 1 and end at a number: the coding sequence opens on the transcript's first
/// base and the span on the first coding base. `None` for any other entry.
fn ends_of_a_span_from_the_first_base(e: &Entry) -> Option<(&str, &str)> {
    let [cdna_start, cdna_end, cds_start, cds_end] = transcript_positions(e)?;
    (cdna_start == "1"
        && cds_start == "1"
        && pytext::is_ascii_digits(cdna_end)
        && pytext::is_ascii_digits(cds_end))
    .then_some((cdna_end, cds_end))
}

/// `start_codon_removed_without_five_prime_utr`: a Transcript entry of a deletion of whole
/// codons from the first coding base of a transcript whose coding sequence opens on its
/// first base. The candidate's terms are the reference's and `start_lost`, both sets
/// holding `inframe_deletion` and neither `start_retained_variant`; the cDNA position and
/// the CDS position begin at 1 and end at a number, the same on both sides; and the
/// residues, the same on both sides, have `-` for their alternate half. For such a
/// deletion no line of `start_lost` (Utils/VariationEffect.pm:851-904) finds the start
/// codon kept and none names it lost: the length test is set aside for an in-frame
/// deletion (:870), `_inv_start_altered` (:906-945) returns where the transcript has no
/// 5' UTR (:926) and the peptide test returns on the empty alternate peptide (:875).
pub fn start_codon_removed_without_five_prime_utr_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) {
        return Ok(false);
    }
    let (expected, gone, come) = term_differences(e, a);
    if !gone.is_empty()
        || !is(&come, &[START_LOST_TERM])
        || !expected.contains(INFRAME_DELETION_TERM)
        || expected.contains(START_RETAINED_TERM)
    {
        return Ok(false);
    }
    let ends = ends_of_a_span_from_the_first_base(e);
    if ends.is_none() || ends != ends_of_a_span_from_the_first_base(a) {
        return Ok(false);
    }
    let residues = e.fields.first_of(&RESIDUE_FIELDS, ABSENT);
    if residues != a.fields.first_of(&RESIDUE_FIELDS, ABSENT) {
        return Ok(false);
    }
    Ok(matches!(
        residues.split_once('/'),
        Some((removed, "-")) if !blank(removed)
    ))
}

/// The groups of terms the reference writes on a row of
/// `frameshift_inside_a_held_stop_codon`, one of them whole in the place of the candidate's
/// `frameshift_variant`.
const HELD_STOP_UNSHIFTED_TERMS: [&[&str]; 3] = [
    &[CODING_SEQUENCE_TERM],
    &[INFRAME_DELETION_TERM],
    &[INFRAME_INSERTION_TERM, CODING_SEQUENCE_TERM],
];

/// `frameshift_inside_a_held_stop_codon`: a Transcript entry of an insertion or a deletion
/// whose length change is no multiple of three and whose first base lies inside the codon of
/// a stop that the transcript's peptide holds as a stop and reads on past. `frameshift`
/// returns 0 where the reference peptide of the allele's codons begins with a stop
/// (Utils/VariationEffect.pm:1554), so the reference writes no `frameshift_variant` and the
/// candidate writes it in the place of one group of `HELD_STOP_UNSHIFTED_TERMS`, every
/// other term standing on both sides.
///
/// The two entries show the shape by the cells they share: the same residues, whose
/// reference half opens with `*` and whose alternate half holds `X`, and the same `Codons`
/// cell, whose two halves differ in length by no multiple of three, a half written `-`
/// holding no base. A reference entry that holds `stop_retained_variant` or `stop_lost` is
/// refused.
pub fn frameshift_inside_a_held_stop_codon_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) {
        return Ok(false);
    }
    let (expected, gone, come) = term_differences(e, a);
    if expected.contains(STOP_RETAINED_TERM) || expected.contains(STOP_LOST_TERM) {
        return Ok(false);
    }
    if !is(&come, &[FRAMESHIFT_TERM])
        || !HELD_STOP_UNSHIFTED_TERMS
            .iter()
            .any(|group| is(&gone, group))
    {
        return Ok(false);
    }
    let residues = e.fields.first_of(&RESIDUE_FIELDS, ABSENT);
    if residues != a.fields.first_of(&RESIDUE_FIELDS, ABSENT) {
        return Ok(false);
    }
    match residues.split_once('/') {
        Some((reference, alternate)) if reference.starts_with('*') && alternate.contains('X') => {}
        _ => return Ok(false),
    }
    let codons = e.fields.first_of(&CODONS_FIELDS, ABSENT);
    if codons != a.fields.first_of(&CODONS_FIELDS, ABSENT) {
        return Ok(false);
    }
    let Some((reference, alternate)) = codons.split_once('/') else {
        return Ok(false);
    };
    let bases = |half: &str| match half {
        "-" => 0,
        _ => half.chars().count(),
    };
    Ok(!bases(reference)
        .abs_diff(bases(alternate))
        .is_multiple_of(3))
}

/// The terms of an insertion or deletion inside the coding sequence.
const LENGTH_TERMS: [&str; 3] = [
    FRAMESHIFT_TERM,
    INFRAME_INSERTION_TERM,
    INFRAME_DELETION_TERM,
];

/// The CDS start and end of an entry whose cDNA span is its CDS span and opens on a base
/// of the start codon: the coding sequence opens on the transcript's first base, and the
/// smaller of the two bounds is 1, 2 or 3. `None` for any other entry, a bound that is no
/// number without a leading zero included.
fn span_on_the_start_codon(e: &Entry) -> Option<(&str, &str)> {
    let [cdna_start, cdna_end, cds_start, cds_end] = transcript_positions(e)?;
    let number = |text: &str| pytext::is_ascii_digits(text) && !text.starts_with('0');
    if cdna_start != cds_start || cdna_end != cds_end || !number(cds_start) || !number(cds_end) {
        return None;
    }
    let opens = if (cds_start.len(), cds_start) <= (cds_end.len(), cds_end) {
        cds_start
    } else {
        cds_end
    };
    matches!(opens, "1" | "2" | "3").then_some((cds_start, cds_end))
}

/// `start_codon_kept_without_five_prime_utr`: a Transcript entry of an insertion or
/// deletion that overlaps the start codon of a transcript whose coding sequence opens on
/// its first base. The candidate's terms are the reference's without `start_lost` and with
/// `start_retained_variant`, which the reference does not hold, every other term the same;
/// both sets hold `frameshift_variant`, `inframe_insertion` or `inframe_deletion`; and the
/// cDNA span is the CDS span and opens on one of the first three coding bases, the same on
/// both sides. `_ins_del_start_altered` (Utils/VariationEffect.pm:1028-1075) reads a start
/// as kept inside `if ($utr)` alone (:1061-1065), so on such a transcript Ensembl VEP
/// writes no `start_retained_variant` (:947-963), and `start_lost` where the allele is a
/// frameshift (:870) or where its peptide test holds (:878-882).
pub fn start_codon_kept_without_five_prime_utr_pair(e: &Entry, a: &Entry) -> Res<bool> {
    if !transcript_row(e) {
        return Ok(false);
    }
    let (expected, gone, come) = term_differences(e, a);
    let actual = terms(&a.consequence);
    if !is(&come, &[START_RETAINED_TERM])
        || !(gone.is_empty() || is(&gone, &[START_LOST_TERM]))
        || actual.contains(START_LOST_TERM)
        || ![&expected, &actual]
            .iter()
            .all(|held| LENGTH_TERMS.iter().any(|term| held.contains(term)))
    {
        return Ok(false);
    }
    let span = span_on_the_start_codon(e);
    Ok(span.is_some() && span == span_on_the_start_codon(a))
}

/// A row rule of a stop codon class.
type RowPair = fn(&Entry, &Entry) -> Res<bool>;

/// The row rules of the stop codon classes beside which a plugin's field follows the terms.
/// Three stand with the rule's name among the classes of a codon that stands, by which a
/// comparison beside the cache holds the rule to a transcript that shows the codon. The
/// rules of `stop_codon_settled_before_unknown_bases` and of
/// `frameshift_inside_a_held_stop_codon` stand with none: they read nothing of the row's
/// transcript.
const STOP_CODON_PLUGIN_VALUE_PAIRS: [(RowPair, Option<LandedRule>); 5] = [
    (
        stop_retained_without_stop_codon_test_pair,
        Some(LandedRule::StopRetainedWithoutStopCodonTest),
    ),
    (
        stop_codon_overlap_insertion_range_pair,
        Some(LandedRule::StopCodonOverlapInsertionRange),
    ),
    (
        stop_retained_three_base_stop_codon_test_pair,
        Some(LandedRule::StopRetainedThreeBaseStopCodonTest),
    ),
    (stop_codon_settled_before_unknown_bases_pair, None),
    (frameshift_inside_a_held_stop_codon_pair, None),
];

/// `plugin_value_follows_stop_codon_term_pair`: a field of a plugin that runs on a row's
/// terms, which one engine alone writes on a row whose two Consequence sets differ by one
/// of ten causes, where `term_pair_following_cell` reads the two values from the two
/// sets of terms.
///
/// Five causes are the classes of a coding end (`coding_end_class`): each takes the pair
/// by the transcript's coding ends in the candidate cache, so none takes one where the
/// class is handed no cache or the cache holds no reading of the transcript. Five are
/// the rules of `STOP_CODON_PLUGIN_VALUE_PAIRS`. Beside the cache three of them take the
/// pair only on a transcript that shows the rule's stop codon
/// (`landed_rule_takes_beside_a_cache`), and without the cache by the two rows alone. The
/// rules of `stop_codon_settled_before_unknown_bases` and of
/// `frameshift_inside_a_held_stop_codon` take the pair by the two rows alone, beside the
/// cache and without it.
pub fn follows_stop_codon_term_pair(c: &Consult, ev: &str, av: &str) -> Res<bool> {
    let Some(a) = c.a else {
        return Ok(false);
    };
    let beside_the_cache = match c.downstream_sequences {
        Some(sequences) => Some(sequences.coding_ends(c.e)?),
        None => None,
    };
    let mut taken = beside_the_cache.is_some_and(|ends| coding_end_class(c.e, a, ends).is_some());
    for (pair, rule) in STOP_CODON_PLUGIN_VALUE_PAIRS {
        if taken {
            break;
        }
        taken = pair(c.e, a)?
            && rule.is_none_or(|rule| {
                beside_the_cache
                    .is_none_or(|ends| landed_rule_takes_beside_a_cache(rule, c.e, a, ends))
            });
    }
    Ok(taken && term_pair_following_cell(c.e, a, c.name, ev, av))
}

#[cfg(test)]
mod tests;
