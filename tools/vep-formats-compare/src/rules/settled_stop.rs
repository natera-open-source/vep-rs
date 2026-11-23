// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row class of an alternate allele that holds unknown bases after bases
//! that settle a stop codon: `stop_codon_settled_before_unknown_bases`.

use crate::consts::HGVSC_FIELDS;
use crate::consts::HGVSP_FIELDS;
use crate::consts::HGVS_OFFSET_FIELDS;
use crate::consts::NON_ZERO_INTEGER_RE_PATTERN;
use crate::consts::OUTSIDE_TRANSCRIPT_TERMS;
use crate::entry::Entry;
use crate::error::Result;
use crate::fields::field_named;
use crate::genetic_codes::{genetic_code, CodonTable, Reading};
use crate::genetic_codes::{GeneticCode, MODULE_CODES, STANDARD_CODE};
use crate::pyfmt;
use crate::rules::downstream::{ambiguous_residue, DownstreamSequences};
use crate::rules::edited_model::CODONS_FIELDS;
use crate::rules::first_field;
use crate::rules::following::term_pair_following_cell;
use crate::rules::hgvsp_stop::STOP_FRAMESHIFT_RE_PATTERN;
use crate::rules::term_pairs::INFRAME_DELETION_TERM;
use crate::rules::term_pairs::PROTEIN_ALTERING_TERM;
use crate::rules::transcript_model::TranscriptModel;
use crate::rules::whole_text;
use crate::self_parity::reverse_complement;
use crate::shared::_ten_digits;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::feature_kind;
use crate::shared::inserted_between_neighbours;
use crate::shared::insertion_row;
use crate::shared::length_changing_row;
use crate::shared::position_bounds;
use crate::shared::stop_codon_test_coordinates;
use crate::shared::stop_codon_test_not_made;
use crate::shared::CODING_SEQUENCE_TERM;
use crate::shared::FRAMESHIFT_TERM;
use crate::shared::INFRAME_INSERTION_TERM;
use crate::shared::STOP_GAINED_TERM;
use crate::shared::STOP_RETAINED_TERM;
use std::collections::BTreeSet;
use std::sync::LazyLock;

/// The residues cell of the text formats and of a JSON entry.
pub const RESIDUE_FIELDS: &[&str] = &["Amino_acids", "amino_acids"];

/// The three position columns of the text formats, each with the two keys a JSON entry
/// writes its bounds under.
const POSITION_COLUMNS: [(&str, &str, &str); 3] = [
    ("cDNA_position", "cdna_start", "cdna_end"),
    ("CDS_position", "cds_start", "cds_end"),
    ("Protein_position", "protein_start", "protein_end"),
];

static FRAMESHIFT_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_FRAMESHIFT_RE_PATTERN));
static NON_ZERO_INTEGER_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NON_ZERO_INTEGER_RE_PATTERN));

/// The terms the candidate alone writes beside `stop_gained` on a row whose
/// length is whole codons. Membership only.
pub const INFRAME_TERMS: &[&str] = &[
    INFRAME_INSERTION_TERM,
    INFRAME_DELETION_TERM,
    PROTEIN_ALTERING_TERM,
];

/// The bases a nucleotide letter stands for; empty for a letter that names
/// none.
fn bases_of(letter: u8) -> &'static [u8] {
    match letter {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' | b'X' => b"ACGT",
        _ => b"",
    }
}

/// The residue of one codon as both engines' translation reads it: the
/// code's for a codon of A, C, G and T; for a codon that holds another
/// nucleotide letter the residue every reading of it shares, `B` for
/// aspartate or asparagine, `Z` for glutamate or glutamine, `X` otherwise.
/// `None` for a letter that names no base.
fn residue(codon: &[u8], code: &GeneticCode) -> Option<u8> {
    if let Some(residue) = code.amino_acid(codon) {
        return Some(residue);
    }
    let &[first, second, third] = codon else {
        return None;
    };
    let mut read: BTreeSet<u8> = BTreeSet::new();
    for &a in bases_of(first) {
        for &b in bases_of(second) {
            for &c in bases_of(third) {
                read.insert(code.amino_acid(&[a, b, c])?);
            }
        }
    }
    let read: Vec<u8> = read.into_iter().collect();
    match read.as_slice() {
        [] => None,
        [one] => Some(*one),
        [b'D', b'N'] => Some(b'B'),
        [b'E', b'Q'] => Some(b'Z'),
        _ => Some(b'X'),
    }
}

/// The alternate peptide of the alternate half of a `Codons` cell that holds
/// a letter outside A, C, G and T after bases that settle a stop codon: a
/// stop codon stands among its whole codons before the first of them that
/// holds such a letter. The peptide is the translation of the whole codons
/// under the standard genetic code, with the `X` both engines append for bases past
/// the last whole codon unless the peptide is the stop alone. `None` for an
/// alternate half of A, C, G and T alone, for one with a letter that names no
/// base, and for one whose known bases settle no stop codon before the first
/// codon that holds another letter.
pub fn settled_alternate_peptide(alternate_codons: &str) -> Option<String> {
    let code = MODULE_CODES.get(&STANDARD_CODE)?;
    let bases = alternate_codons.to_ascii_uppercase();
    let bases = bases.as_bytes();
    let known = |base: &u8| b"ACGT".contains(base);
    if bases.iter().all(known) {
        return None;
    }
    let mut settled = false;
    for codon in bases.chunks_exact(3) {
        if !codon.iter().all(known) {
            break;
        }
        if code.amino_acid(codon) == Some(b'*') {
            settled = true;
            break;
        }
    }
    if !settled {
        return None;
    }
    let mut peptide: Vec<u8> = bases
        .chunks_exact(3)
        .map(|codon| residue(codon, code))
        .collect::<Option<_>>()?;
    if !bases.len().is_multiple_of(3) && peptide != b"*" {
        peptide.push(b'X');
    }
    String::from_utf8(peptide).ok()
}

/// Whether two sets of terms differ as those of a stop codon the inserted
/// bases keep, given the terms the reference alone holds and those the
/// candidate alone holds, each in order, and whether the reference's hold
/// `stop_retained_variant` (`retained`). The candidate alone holds
/// `inframe_insertion`, which needs both peptides
/// (Utils/VariationEffect.pm:1113-1115) and cuts the alternate one after its
/// first stop (:1126). Where the reference holds `stop_retained_variant`,
/// from the second arm of `stop_retained` (:1318), whose edit test is not
/// made for a letter that names no base (:1498), nothing else differs.
/// Otherwise the candidate alone holds `stop_retained_variant` as well, the
/// reference alone `frameshift_variant` (:1564) or `coding_sequence_variant`
/// (:1616-1649), and the entry is an insertion between two neighbouring
/// positions with the stop codon test's coordinates on both sides, whose
/// moved start fails the reference's overlap (:1408-1415).
fn kept_stop_terms(
    e: &Entry,
    a: &Entry,
    retained: bool,
    reference_only: &[&str],
    candidate_only: &[&str],
) -> bool {
    if retained {
        return reference_only.is_empty() && candidate_only == [INFRAME_INSERTION_TERM];
    }
    (reference_only == [FRAMESHIFT_TERM] || reference_only == [CODING_SEQUENCE_TERM])
        && candidate_only == [INFRAME_INSERTION_TERM, STOP_RETAINED_TERM]
        && insertion_row(e)
        && inserted_between_neighbours(e)
        && stop_codon_test_coordinates(e)
        && stop_codon_test_coordinates(a)
}

/// Whether the reference entry `e` and the candidate entry `a` of one key
/// have the shape of `stop_codon_settled_before_unknown_bases`: a Transcript
/// entry whose `Codons` cell is the same on both sides and has an alternate
/// half `settled_alternate_peptide` reads; no residues in the reference and,
/// in the candidate, a reference half and that peptide; and terms that differ
/// by the gained stop alone or by the kept stop alone. For a gained stop, where
/// the reference writes `frameshift_variant` the candidate adds `stop_gained`
/// to the same terms; otherwise the reference alone writes
/// `coding_sequence_variant` and the candidate alone `stop_gained`, with at
/// most one of `INFRAME_TERMS`.
///
/// The reference can write `stop_retained_variant` in the place of either
/// term: for an insertion or deletion without a peptide `stop_retained` holds
/// where the overlap with the stop codon holds, read with an insertion's start
/// moved by the inserted length, and the stop codon test returns before it
/// reads a base (Utils/VariationEffect.pm:1318, :1408-1415, :1498). On such a
/// row, an insertion's or a deletion's (`length_changing_row`) with both
/// entries showing that the test is not made (`stop_codon_test_not_made`), the
/// reference alone writes `stop_retained_variant` and the candidate alone
/// `stop_gained`, with `frameshift_variant`, one of `INFRAME_TERMS` or no
/// other term.
///
/// For a kept stop the two sets differ as `kept_stop_terms` reads them.
pub fn stop_codon_settled_before_unknown_bases_pair(e: &Entry, a: &Entry) -> Result<bool> {
    let last = e.key.rsplit('|').next().unwrap_or("");
    if feature_kind(last.split('#').next().unwrap_or("")) != "transcript" {
        return Ok(false);
    }
    let codons = first_field(e, CODONS_FIELDS);
    if blank(codons) || codons != first_field(a, CODONS_FIELDS) {
        return Ok(false);
    }
    let Some(peptide) = codons
        .split_once('/')
        .and_then(|(_, alternate)| settled_alternate_peptide(alternate))
    else {
        return Ok(false);
    };
    if !blank(first_field(e, RESIDUE_FIELDS)) {
        return Ok(false);
    }
    match first_field(a, RESIDUE_FIELDS).split_once('/') {
        Some((reference, alternate)) if !reference.is_empty() && alternate == peptide => {}
        _ => return Ok(false),
    }
    let (expected_set, actual_set) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    let expected: BTreeSet<&str> = expected_set.split(',').collect();
    let actual: BTreeSet<&str> = actual_set.split(',').collect();
    let reference_only: Vec<&str> = expected.difference(&actual).copied().collect();
    let candidate_only: Vec<&str> = actual.difference(&expected).copied().collect();
    let retained = expected.contains(STOP_RETAINED_TERM);
    if kept_stop_terms(e, a, retained, &reference_only, &candidate_only) {
        return Ok(true);
    }
    if expected.contains(FRAMESHIFT_TERM) {
        return Ok(reference_only.is_empty() && candidate_only == [STOP_GAINED_TERM]);
    }
    let beside_the_stop: Vec<&str> = candidate_only
        .iter()
        .copied()
        .filter(|term| *term != STOP_GAINED_TERM)
        .collect();
    if !candidate_only.contains(&STOP_GAINED_TERM) || beside_the_stop.len() > 1 {
        return Ok(false);
    }
    if reference_only == [STOP_RETAINED_TERM] {
        return Ok(length_changing_row(e)?
            && stop_codon_test_not_made(e)
            && stop_codon_test_not_made(a)
            && beside_the_stop
                .iter()
                .all(|term| *term == FRAMESHIFT_TERM || INFRAME_TERMS.contains(term)));
    }
    Ok(reference_only == [CODING_SEQUENCE_TERM]
        && beside_the_stop
            .iter()
            .all(|term| INFRAME_TERMS.contains(term)))
}

/// Whether the differing values `ev` and `av` of the field `name` on a row
/// `stop_codon_settled_before_unknown_bases_pair` takes are the ones the
/// class reads: the candidate's residues, which the pair's predicate holds to
/// the codons; a protein notation that names a stop where the reference
/// writes none; the HGVS_OFFSET written beside that notation alone; and the
/// cells that follow the two sets of terms (`term_pair_following_cell`).
pub fn stop_codon_settled_before_unknown_bases_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    ev: &str,
    av: &str,
) -> bool {
    if field_named(name, RESIDUE_FIELDS) {
        return blank(ev) && !blank(av);
    }
    if field_named(name, HGVSP_FIELDS) {
        return blank(ev) && pyfmt::partition(av, ":p.").2.contains("Ter");
    }
    if field_named(name, HGVS_OFFSET_FIELDS) {
        return blank(ev)
            && av.parse::<i64>().is_ok()
            && blank(first_field(e, HGVSP_FIELDS))
            && pyfmt::partition(first_field(a, HGVSP_FIELDS), ":p.")
                .2
                .contains("Ter");
    }
    term_pair_following_cell(e, a, name, ev, av)
}

/// What the cells of a row `carried_insertion` takes state: the inserted bases on the
/// forward strand, in upper case; the description of the candidate's protein notation; its
/// HGVS_OFFSET; and the two positions of the Location.
struct CarriedInsertion<'a> {
    allele: String,
    description: &'a str,
    offset: i64,
    flanks: (i64, i64),
}

/// Whether the key of an entry names a Transcript feature.
fn transcript_row(e: &Entry) -> bool {
    let last = e.key.rsplit('|').next().unwrap_or("");
    feature_kind(last.split('#').next().unwrap_or("")) == "transcript"
}

/// The two neighbouring positions an entry's Location names, a text row's own or the one
/// its record gives the allele of a VCF or JSON entry; `None` for any other Location.
fn flanks(e: &Entry) -> Option<(i64, i64)> {
    let location = e.location.as_deref().or(e.allele_location.as_deref())?;
    let (first, _, last) = pyfmt::partition(pyfmt::rpartition(location, ":").2, "-");
    if !(_ten_digits(first) && _ten_digits(last)) {
        return None;
    }
    let (first, last) = (pyfmt::int_of(first)?, pyfmt::int_of(last)?);
    (first.checked_add(1) == Some(last)).then_some((first, last))
}

/// The reading of the reference entry `e` and the candidate entry `a` of one key as a row
/// of an insertion with unknown bases that both engines annotate in an intron and the
/// candidate alone describes where the HGVS 3' shift ends; `None` for any other row:
///
/// - a Transcript entry on which the reference writes no HGVSc, no HGVSp and no
///   HGVS_OFFSET (TranscriptVariationAllele.pm:1397 and :1684 return for an allele with a
///   letter other than A, C, G, T and `-`; OutputFactory.pm:1764 and :1776 write the offset
///   beside a notation alone), and the candidate no HGVSc, a protein notation that names a
///   stop and an HGVS_OFFSET that is a whole number other than zero;
/// - two equal sets of terms, none of `OUTSIDE_TRANSCRIPT_TERMS` among them, and no cDNA,
///   CDS or protein position and no `Codons` cell on either side;
/// - an insertion's row (`insertion_row`) between two neighbouring positions
///   (`inserted_between_neighbours`, which refuses a record that replaces the two bases),
///   whose allele holds a letter other than A, C, G and T, each of its letters naming bases.
fn carried_insertion<'a>(e: &Entry, a: &'a Entry) -> Option<CarriedInsertion<'a>> {
    if !transcript_row(e) {
        return None;
    }
    let unwritten = |entry: &Entry, names: &[&str]| blank(first_field(entry, names));
    if !(unwritten(e, HGVSC_FIELDS)
        && unwritten(e, HGVSP_FIELDS)
        && unwritten(e, HGVS_OFFSET_FIELDS)
        && unwritten(a, HGVSC_FIELDS))
    {
        return None;
    }
    let description = pyfmt::partition(first_field(a, HGVSP_FIELDS), ":p.").2;
    let offset = first_field(a, HGVS_OFFSET_FIELDS);
    if !description.contains("Ter") || !NON_ZERO_INTEGER_WHOLE.is_match(offset) {
        return None;
    }
    let offset = offset.parse::<i64>().ok()?;
    let terms = consequence_set(&e.consequence);
    if terms != consequence_set(&a.consequence)
        || terms
            .split(',')
            .any(|term| OUTSIDE_TRANSCRIPT_TERMS.contains(&term))
    {
        return None;
    }
    let annotated_in_an_intron = |entry: &Entry| {
        unwritten(entry, CODONS_FIELDS)
            && POSITION_COLUMNS.iter().all(|(column, start, end)| {
                position_bounds(entry, column, start, end) == (false, false)
            })
    };
    if !annotated_in_an_intron(e)
        || !annotated_in_an_intron(a)
        || !insertion_row(e)
        || !inserted_between_neighbours(e)
    {
        return None;
    }
    let allele = pyfmt::upper(&e.allele);
    let known_alone = allele.bytes().all(|base| b"ACGT".contains(&base));
    let names_bases = allele.bytes().all(|base| !bases_of(base).is_empty());
    if known_alone || !names_bases {
        return None;
    }
    Some(CarriedInsertion {
        allele,
        description,
        offset,
        flanks: flanks(e)?,
    })
}

/// Where the shifted point of an insertion stands in a transcript's coding sequence, as
/// its `translateable_seq` spells it: the count of the sequence's bases before the point,
/// and the count of those, directly before the point, that are coding bases of the exon
/// that holds it, from the exon's first base or, in an exon that opens before the coding
/// sequence, from the first coding base. The `N` the sequence opens with for the start
/// exon's phase are bases of the first count and of no exon.
struct ShiftedPoint {
    bases_before: usize,
    of_its_exon: usize,
}

/// The shifted point of an insertion written between the two neighbouring positions
/// `flanks` and moved `offset` positions along the genome. `None` unless both positions
/// lie outside every exon as written and, moved, lie inside one exon on two neighbouring
/// bases of the coding sequence.
fn shifted_point(model: &TranscriptModel, flanks: (i64, i64), offset: i64) -> Option<ShiftedPoint> {
    let exon_of = |position: i64| {
        model
            .exons
            .iter()
            .position(|&(start, end)| start <= position && position <= end)
    };
    let (low, high) = flanks;
    if exon_of(low).is_some() || exon_of(high).is_some() {
        return None;
    }
    let (low, high) = (low.checked_add(offset)?, high.checked_add(offset)?);
    let exon = exon_of(low)?;
    if exon_of(high)? != exon {
        return None;
    }
    // The exon's first base in the transcript's direction is its highest on the reverse
    // strand.
    let (start, end) = model.exons[exon];
    let (before, after, opening) = if model.strand == 1 {
        (low, high, start)
    } else {
        (high, low, end)
    };
    let (before, after) = (model.cdna(before)?, model.cdna(after)?);
    if before.checked_add(1) != Some(after)
        || before < model.coding_start
        || after > model.coding_end
    {
        return None;
    }
    let first_coding = model.cdna(opening)?.max(model.coding_start);
    let bases_before =
        i128::from(before) - i128::from(model.coding_start) + 1 + i128::from(model.padding);
    let of_its_exon = i128::from(before) - i128::from(first_coding) + 1;
    Some(ShiftedPoint {
        bases_before: usize::try_from(bases_before).ok()?,
        of_its_exon: usize::try_from(of_its_exon)
            .ok()
            .filter(|count| *count > 0)?,
    })
}

/// The codon of `sequence` at `place`, counted in whole codons from 0; `None` past its
/// last whole codon.
fn codon_at(sequence: &str, place: usize) -> Option<&[u8]> {
    sequence.as_bytes().get(3 * place..3 * place + 3)
}

/// The first codon an insertion changes: its place, counted in codons of the coding
/// sequence from 0, the residue the reference reads there (`None` past the sequence's last
/// whole codon) and the residue the changed codon reads.
struct ChangedCodon {
    place: usize,
    reference: Option<u8>,
    residue: u8,
}

/// The three-letter name a protein description gives the residue `letter`, read in either
/// case, as `rules::frameshift_count` names the residues of a frameshift's description;
/// `None` for a stop, for the unknown residue `X` and for a letter that names no residue.
fn residue_name(letter: u8) -> Option<&'static str> {
    Some(match letter.to_ascii_uppercase() {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'Q' => "Gln",
        b'E' => "Glu",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        b'U' => "Sec",
        b'O' => "Pyl",
        b'B' => "Asx",
        b'Z' => "Glx",
        b'J' => "Xle",
        _ => return None,
    })
}

/// The stop codon an insertion settles: `bases`, on the transcript's strand, are written
/// into the coding sequence `coding` behind its first `point` bases, and the whole codons
/// are read from the codon that holds the point, which opens with its reference bases
/// where the point splits it, to the codon that holds the last inserted base. The result
/// is the place of the first stop codon of `code` among them, counted in codons of the
/// coding sequence from 0, where it and every codon read before it are of A, C, G and T
/// alone; with it comes the first codon read whose residue differs from the reference's,
/// with both residues, `None` where each one up to that stop codon reads the reference's
/// residue. `None` where no such stop codon stands there.
fn settled_stop_codon(
    coding: &str,
    point: usize,
    bases: &str,
    code: GeneticCode,
) -> Option<(usize, Option<ChangedCodon>)> {
    if !coding.is_ascii() || bases.is_empty() || point >= coding.len() {
        return None;
    }
    let changed = [&coding[..point], bases, &coding[point..]].concat();
    let mut first_changed = None;
    for place in point / 3..=(point + bases.len() - 1) / 3 {
        // The code holds no residue for a codon with a letter other than A, C, G and T.
        let residue = code.amino_acid(codon_at(&changed, place)?)?;
        if first_changed.is_none() {
            let reference = codon_at(coding, place).map(|held| {
                code.amino_acid(held)
                    .unwrap_or_else(|| ambiguous_residue(code, held))
            });
            if reference != Some(residue) {
                first_changed = Some(ChangedCodon {
                    place,
                    reference,
                    residue,
                });
            }
        }
        if residue == b'*' {
            return Some((place, first_changed));
        }
    }
    None
}

/// Whether the transcript's model and coding sequence give the candidate's cells of `row`:
///
/// - the offset's sign is the transcript's strand, and its size is smaller than the count
///   of known bases the allele opens with on that strand, so the shift stops at a known
///   base and passes the allele's first bases alone;
/// - the insertion as written lies between the two positions of the Location, both outside
///   every exon, and moved by the offset both lie inside one exon on bases of the coding
///   sequence (`shifted_point`);
/// - the shift passes every coding base of that exon that stands before the point, and
///   those bases of the coding sequence are, letter for letter with the case set aside, the
///   last of the allele's bases the shift passes on the transcript's strand, as a shift
///   that moves the allele out of the intron and past the exon's first bases leaves them;
/// - the coding base behind the point is another letter, the case set aside, than the first
///   of the allele's bases the shift does not pass, so the shift ends at that point;
/// - the allele on the transcript's strand, as many of its first bases moved to its end as
///   the shift passes, settles a stop codon of `code` at that point (`settled_stop_codon`);
/// - a frameshift's description (`Xaa<number>XaafsTer<count>`, a pair of prediction
///   parentheses set aside) names the reference's residue at the first codon whose residue
///   differs from the reference's, numbers that codon, names the residue the changed codon
///   reads there (`residue_name`) and counts from it to that stop codon, the first counted.
///   Any other description is read by the conditions above alone.
fn settled_where_the_shift_ends(
    row: &CarriedInsertion,
    model: &TranscriptModel,
    coding: &str,
    code: GeneticCode,
) -> Option<()> {
    if row.offset.signum() != model.strand {
        return None;
    }
    let passed = usize::try_from(row.offset.unsigned_abs()).ok()?;
    let read = if model.strand == -1 {
        reverse_complement(&row.allele)
    } else {
        row.allele.clone()
    };
    let opening = read
        .bytes()
        .take_while(|base| b"ACGT".contains(base))
        .count();
    if passed >= opening {
        return None;
    }
    let point = shifted_point(model, row.flanks, row.offset)?;
    // A shift shorter than the exon's coding bases before the point is refused here.
    let passed_before_them = passed.checked_sub(point.of_its_exon)?;
    let first_of_them = point.bases_before.checked_sub(point.of_its_exon)?;
    let of_its_exon = coding.as_bytes().get(first_of_them..point.bases_before)?;
    if !read.as_bytes()[passed_before_them..passed].eq_ignore_ascii_case(of_its_exon) {
        return None;
    }
    let behind = coding.as_bytes().get(point.bases_before)?;
    if behind.eq_ignore_ascii_case(read.as_bytes().get(passed)?) {
        return None;
    }
    let turned = format!("{}{}", &read[passed..], &read[..passed]);
    let (stop, changed) = settled_stop_codon(coding, point.bases_before, &turned, code)?;
    let described = row
        .description
        .strip_prefix('(')
        .and_then(|inner| inner.strip_suffix(')'))
        .unwrap_or(row.description);
    let Some(parts) = FRAMESHIFT_WHOLE.captures(described) else {
        return Some(());
    };
    let changed = changed?;
    let number = |text: &str| text.parse::<usize>().ok();
    let named = |residue: Option<u8>| residue.and_then(residue_name);
    (named(changed.reference) == Some(&parts[1])
        && number(&parts[2]) == Some(changed.place + 1)
        && named(Some(changed.residue)) == Some(&parts[3])
        && number(&parts[4]) == Some(stop - changed.place + 1))
    .then_some(())
}

/// Whether the differing values `ev` and `av` of the field `name` are the protein notation
/// or the HGVS_OFFSET of a row without a `Codons` cell that the class reads beside the
/// transcripts of the candidate's cache (`sequences`): an insertion with unknown bases that
/// is written in an intron, where both engines annotate it, and that the HGVS 3' shift
/// carries into the coding sequence, where its known bases settle a stop codon before its
/// first unknown base. The reference writes neither cell, the two entries are such a row by
/// their own cells (`carried_insertion`), and the row's transcript in the cache, one that
/// states its model under a genetic code a table is held for, gives the candidate's two
/// cells (`settled_where_the_shift_ends`). The two cells of one row are taken together.
pub fn carried_insertion_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    ev: &str,
    av: &str,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    if !(field_named(name, HGVSP_FIELDS) || field_named(name, HGVS_OFFSET_FIELDS))
        || !blank(ev)
        || blank(av)
    {
        return Ok(false);
    }
    let Some(row) = carried_insertion(e, a) else {
        return Ok(false);
    };
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let (Some(model), CodonTable::Held(table)) = (transcript.model.as_ref(), &transcript.table)
    else {
        return Ok(false);
    };
    let Some(code) = genetic_code(*table, Reading::VepRs) else {
        return Ok(false);
    };
    Ok(settled_where_the_shift_ends(&row, model, &transcript.coding, code).is_some())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_stop_codon_before_the_first_unknown_codon_gives_the_peptide() {
        // A sense codon, a stop codon, an unknown codon and a codon whose
        // first base is unknown.
        assert_eq!(
            settled_alternate_peptide("GAATAGNNNNCA").as_deref(),
            Some("E*XX")
        );
        // The same with one base past the last whole codon.
        assert_eq!(
            settled_alternate_peptide("GAATAGNNNN").as_deref(),
            Some("E*XX")
        );
        // A codon every reading of which is one residue keeps it.
        assert_eq!(
            settled_alternate_peptide("TAACTNGGA").as_deref(),
            Some("*LG")
        );
        assert_eq!(settled_alternate_peptide("taaN").as_deref(), Some("*"));
        // Aspartate or asparagine reads `B`, glutamate or glutamine `Z`, and a
        // codon every reading of which is a stop codon the stop.
        assert_eq!(
            settled_alternate_peptide("TAGRATSAATAR").as_deref(),
            Some("*BZ*")
        );
    }

    #[test]
    fn no_stop_codon_before_the_first_unknown_codon_gives_none() {
        assert_eq!(settled_alternate_peptide("GAAGCTNNNNCA"), None);
        // The stop codon stands after the unknown bases.
        assert_eq!(settled_alternate_peptide("GAANNNTAG"), None);
        // The stop codon itself holds the unknown base.
        assert_eq!(settled_alternate_peptide("GAATANGGA"), None);
        // A stop codon of the vertebrate mitochondrial code alone.
        assert_eq!(settled_alternate_peptide("GAAAGANNNN"), None);
        // Known bases alone are no allele of the class.
        assert_eq!(settled_alternate_peptide("GAATAGGGA"), None);
        // A letter that names no base.
        assert_eq!(settled_alternate_peptide("TAGJJJ"), None);
        assert_eq!(settled_alternate_peptide(""), None);
    }

    /// An entry that carries the protein notation `notation`; none for an
    /// empty text.
    fn noted(notation: &str) -> Entry {
        let mut entry = Entry::default();
        if !notation.is_empty() {
            entry
                .fields
                .insert("HGVSp".to_string(), notation.to_string());
        }
        entry
    }

    #[test]
    fn an_offset_is_read_beside_a_notation_that_names_a_stop() {
        let offset = |e: &Entry, a: &Entry| {
            stop_codon_settled_before_unknown_bases_cell(e, a, "HGVS_OFFSET", "-", "4")
        };
        let stop = noted("ENSP1.1:p.Ile142_Leu143insTer");
        assert!(offset(&noted(""), &stop));
        // The notation of a stop codon left unchanged names it.
        assert!(offset(&noted(""), &noted("ENSP1.1:p.Ter534%3D")));
        // A notation that names no stop is not one the class reads.
        assert!(!offset(&noted(""), &noted("ENSP1.1:p.Asn141_Ile142insGlu")));
        assert!(!offset(&noted(""), &noted("")));
        // The reference writes a notation of its own.
        assert!(!offset(&stop, &stop));
    }
}
