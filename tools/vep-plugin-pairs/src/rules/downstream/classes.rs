// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rules of the Downstream, NMD and Blosum62 classes: each recomputes what the two
//! engines write from the transcripts of the candidate cache.

use super::codes::{blosum62, codon_residue};
use super::readings::{
    cds_pair_readings, earlier_allele_shift, first_allele_shift, later_allele_prints_shifted_pair,
    moves, non_zero_offset, printed_pair_is_moved, shifted_pair_readings, CdsPair,
};
use super::transcript::DownstreamTranscript;
use super::DownstreamSequences;
use crate::pytext;
use crate::rules::coding_end::coding_end_class;
use crate::rules::entry::{blank, consequence_set, field_named, field_of, Entry, Row, ABSENT};
use crate::rules::shared::{
    digits_value, integer_of, is_integer, notation_indel_length, reverse_complement,
    FRAMESHIFT_TERM, HGVSC_FIELDS, HGVS_OFFSET_FIELDS, STRAND_FIELDS,
};
use crate::rules::term_pair::{
    registry_class, start_codon_kept_without_five_prime_utr_pair,
    stop_retained_without_stop_codon_test_pair, START_COOCCURRENCE_SWAP,
};
use crate::rules::{
    Consult, Res, DOWNSTREAM_START_PAST_OWN_CODON, DOWNSTREAM_UTR_FRAME_PAST_CODING_END,
    PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
};

pub const PROTEIN_LENGTH_CHANGE_FIELDS: [&str; 2] = ["ProteinLengthChange", "proteinlengthchange"];
pub const NMD_FIELDS: [&str; 2] = ["NMD", "nmd"];
pub const NMD_ESCAPING: &str = "NMD_escaping_variant";
/// The terms one of which NMD.pm:65 and :83 ask of a row.
pub const NMD_GATE_TERMS: [&str; 4] = [
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
];
pub const CDS_POSITION_FIELDS: [&str; 3] = ["CDS_position", "cds_start", "cds_end"];
pub const PROTEIN_POSITION_FIELDS: [&str; 3] = ["Protein_position", "protein_start", "protein_end"];
pub const AMINO_ACIDS_FIELDS: [&str; 1] = ["Amino_acids"];
pub const BLOSUM62_FIELDS: [&str; 1] = ["BLOSUM62"];

/// `downstream_reverse_strand_rotation`: a `DownstreamProtein` or `ProteinLengthChange`
/// both engines fill on a `frameshift_variant` row of a `STRAND` -1 transcript whose
/// `HGVS_OFFSET` exceeds the number of bases its HGVSc inserts, duplicates or deletes and
/// is not a multiple of it.
pub fn reverse_strand_rotation(e: &Entry, ev: &str, av: &str) -> Res<bool> {
    if blank(ev)
        || blank(av)
        || !consequence_set(&e.consequence)
            .split(',')
            .any(|term| term == FRAMESHIFT_TERM)
    {
        return Ok(false);
    }
    if !STRAND_FIELDS
        .iter()
        .any(|field| e.fields.get(field) == Some("-1"))
    {
        return Ok(false);
    }
    let offset = field_of(e, &HGVS_OFFSET_FIELDS);
    let Some(length) = notation_indel_length(field_of(e, &HGVSC_FIELDS)) else {
        return Ok(false);
    };
    if !is_integer(offset) || length <= 0 {
        return Ok(false);
    }
    // The smallest `i128` has no `i128` magnitude, so both are counted unsigned.
    let bases_shifted = integer_of(offset, "HGVS_OFFSET")?.unsigned_abs();
    let length = length.unsigned_abs();
    Ok(bases_shifted > length && bases_shifted % length != 0)
}

/// The bases a plugin reads under --hgvs for the allele as written, `allele`, on a row the
/// notation moved by `shift`: `shift_feature_seqs` (TranscriptVariationAllele.pm:1330-1345)
/// turns the allele left one base as many times as the shift's length on a forward-strand
/// transcript and as the allele's length less it on a reverse-strand one, and `feature_seq`
/// reads the result on the transcript's strand.
pub fn hgvs_rotated_allele(allele: &str, shift: i128, strand: i8) -> String {
    let letters: Vec<char> = allele.chars().collect();
    let length = letters.len() as u128;
    let shifted = shift.unsigned_abs();
    // On the reverse strand a shift of the allele's length or more leaves no turn.
    let turns = if strand == 1 {
        Some(shifted)
    } else {
        length.checked_sub(shifted)
    };
    let turned: String = match turns {
        Some(turns) if length > 0 && turns > 0 => {
            let turns = (turns % length) as usize;
            letters[turns..].iter().chain(&letters[..turns]).collect()
        }
        _ => allele.to_string(),
    };
    if strand == -1 {
        reverse_complement(&turned)
    } else {
        turned
    }
}

/// The class of the two Downstream start classes whose shape the own pair `start` to `end`
/// has with the bases `allele` in its place, or `None`: a span that starts on a codon's
/// third base, `downstream_start_past_own_codon` where the offset of Downstream.pm:113 lies
/// inside the changed sequence or at its end, `downstream_utr_frame_past_coding_end` where
/// it lies past that end and the transcript has a 3' UTR.
pub fn downstream_start_rule_class(
    transcript: &DownstreamTranscript,
    start: i128,
    end: i128,
    allele: &str,
) -> Option<&'static str> {
    let coding = transcript.coding.len() as i128;
    if start.rem_euclid(3) != 0 || start > end || end > coding {
        return None;
    }
    if start <= start - 1 + allele.chars().count() as i128 + coding - end {
        return Some(DOWNSTREAM_START_PAST_OWN_CODON);
    }
    (!transcript.utr.is_empty()).then_some(DOWNSTREAM_UTR_FRAME_PAST_CODING_END)
}

/// The one class that accepts the Downstream field `name` (a DownstreamProtein or a
/// ProteinLengthChange) whose values are `expected_value` on the reference entry `e` and
/// `actual_value` on the candidate entry `a`, with the class that must be declared beside
/// it (`None` for a row one cause changes), or `None`.
///
/// Both are values, the entries share one consequence set, or two that differ as the
/// registry's `start_cooccurrence_swap` admits, as the row rule of
/// `start_codon_kept_without_five_prime_utr` takes them
/// (`start_codon_kept_without_five_prime_utr_pair`) or as a class of a coding end takes
/// them on the transcript's record in the cache (`coding_end_class`), with
/// `frameshift_variant` and no splice term and one CDS span without a gap, and the cache
/// holds the transcript under a table of the codes. Two CDS spans are one on the row of a later allele of `record_rows` whose
/// reference entry prints the record's span moved by its HGVS_OFFSET or by an earlier
/// allele's (`later_allele_prints_shifted_pair`): the own pair is then the candidate's.
///
/// The candidate's value is its rule on the own pair (`DownstreamTranscript::
/// downstream_value`) and the reference's is Downstream.pm:104-149 on the pair the plugin
/// read: the pair and the allele the notation shifted where the reference's non-zero
/// HGVS_OFFSET reaches the plugin (`shifted_pair_readings`, `hgvs_rotated_allele`),
/// whatever earlier pair the row prints; on a later allele's row without an offset, the
/// printed pair an earlier allele's notation left (`earlier_allele_shift`), with the allele
/// as written; else the own pair.
///
/// Where Downstream.pm gives another value on the own pair the class is
/// `plugin_hgvs_shifted_cds_coordinates`: alone where the candidate's value is that
/// own-pair value, the value Ensembl VEP writes without the flag, and beside the own pair's
/// class (`downstream_start_rule_class`) where it is not, since the start rule changes the
/// value as well. Where Downstream.pm gives one value on both pairs, a DownstreamProtein on
/// a row with one ProteinLengthChange is the own pair's class.
pub fn downstream_class(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    sequences: &DownstreamSequences,
    record_rows: Option<&[Row]>,
) -> Res<Option<(&'static str, Option<&'static str>)>> {
    if [expected_value, actual_value]
        .iter()
        .any(|value| matches!(*value, "-" | ABSENT))
    {
        return Ok(None);
    }
    let length_change = field_named(name, &PROTEIN_LENGTH_CHANGE_FIELDS);
    if length_change && !(is_integer(expected_value) && is_integer(actual_value)) {
        return Ok(None);
    }
    let (terms, candidate_terms) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    // One consequence set, the pair the registry sets aside as `start_cooccurrence_swap`,
    // the pair of `start_codon_kept_without_five_prime_utr`, or a pair a class of a coding
    // end takes on the transcript.
    if terms != candidate_terms
        && registry_class(
            &terms,
            &candidate_terms,
            &e.allele,
            &[START_COOCCURRENCE_SWAP],
        )
        .is_none()
        && !start_codon_kept_without_five_prime_utr_pair(e, a)?
        && coding_end_class(e, a, sequences.coding_ends(e)?).is_none()
    {
        return Ok(None);
    }
    if !terms.split(',').any(|term| term == FRAMESHIFT_TERM) || terms.contains("splice") {
        return Ok(None);
    }
    let own_pairs = cds_pair_readings(a);
    let two_pairs = own_pairs != cds_pair_readings(e);
    let Some(&(Some(start), Some(end))) = own_pairs.first() else {
        return Ok(None);
    };
    let allele = if e.allele == "-" {
        String::new()
    } else {
        e.allele.to_uppercase()
    };
    if e.allele != "-" && (allele.is_empty() || !allele.bytes().all(|b| b"ACGTN".contains(&b))) {
        return Ok(None);
    }
    let shift = non_zero_offset(e)?.unwrap_or(0);
    // A later allele's row that prints another pair reads the pair at its own offset, or
    // with none the pair the last earlier allele's notation left.
    let earlier = if two_pairs {
        earlier_allele_shift(e, record_rows)?.unwrap_or(0)
    } else {
        0
    };
    let read_at = if shift != 0 { shift } else { earlier };
    // The cache is read for a row the notation shifted or whose span starts on a codon's
    // third base.
    if read_at == 0 && (two_pairs || length_change || start % 3 != 0 || start > end) {
        return Ok(None);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(None);
    };
    if transcript.table.is_none() {
        return Ok(None);
    }
    let model = transcript.model.as_ref();
    if read_at != 0 && model.is_none() {
        return Ok(None);
    }
    if two_pairs && !later_allele_prints_shifted_pair(e, a, model, shift, record_rows)? {
        return Ok(None);
    }
    let strand = match model {
        Some(model) => model.strand,
        None if field_of(e, &STRAND_FIELDS) == "-1" => -1,
        None => 1,
    };
    let own_allele = if strand == -1 {
        reverse_complement(&allele)
    } else {
        allele.clone()
    };
    let readings = match model {
        Some(model) if read_at != 0 => shifted_pair_readings(e, model, read_at, &own_pairs),
        _ => vec![((Some(start), Some(end)), None, 0)],
    };
    let own_length_change = field_of(e, &PROTEIN_LENGTH_CHANGE_FIELDS);
    let one_length_change = !blank(own_length_change)
        && own_length_change == field_of(a, &PROTEIN_LENGTH_CHANGE_FIELDS);
    let moves = moves(shift, Some(earlier));
    for (own, shifted, _) in readings {
        let (Some(start), Some(end)) = own else {
            continue;
        };
        // A later allele's reading is the one whose pair, moved by the row's offset or by
        // the earlier allele's, the reference prints.
        if two_pairs
            && !model.is_some_and(|model| {
                moves
                    .iter()
                    .any(|by| printed_pair_is_moved(e, a, model, *by, Some(own)))
            })
        {
            continue;
        }
        let candidate = transcript.downstream_value(start, end, &own_allele, length_change, false);
        if candidate.as_deref() != Some(actual_value) {
            continue;
        }
        let unshifted = transcript.downstream_value(start, end, &own_allele, length_change, true);
        let reference = match shifted {
            Some((first, last)) => transcript.downstream_value(
                first,
                last,
                &hgvs_rotated_allele(&allele, shift, strand),
                length_change,
                true,
            ),
            None => unshifted.clone(),
        };
        if reference.as_deref() != Some(expected_value) {
            continue;
        }
        let start_rule = if length_change {
            None
        } else {
            downstream_start_rule_class(transcript, start, end, &own_allele)
        };
        if reference != unshifted {
            if unshifted.as_deref() == Some(actual_value) {
                return Ok(Some((PLUGIN_HGVS_SHIFTED_CDS_COORDINATES, None)));
            }
            if start_rule.is_some() {
                return Ok(Some((PLUGIN_HGVS_SHIFTED_CDS_COORDINATES, start_rule)));
            }
        } else if one_length_change {
            if let Some(start_rule) = start_rule {
                return Ok(Some((start_rule, None)));
            }
        }
    }
    Ok(None)
}

/// Whether `value` is how the entry's format spells a field the row does not carry: no
/// Extra pair of a default-format row (`ABSENT`; `-` where the Extra text is compared pair
/// by pair), the `-` of a tab cell, and the empty cell of an entry without a Location and
/// without a `variant_allele`, which holds every field of its header. An empty value
/// anywhere else is a value.
pub fn field_absent(e: &Entry, value: &str) -> bool {
    matches!(value, "-" | ABSENT)
        || (value.is_empty() && e.location.is_none() && !e.fields.contains("variant_allele"))
}

/// Whether the NMD `expected_value` of the reference entry `e` and `actual_value` of the
/// candidate entry `a` have the shape of `plugin_hgvs_shifted_cds_coordinates`: one is
/// `NMD_escaping_variant` and the other no value (`field_absent`), the entries share one
/// consequence set, or two that differ as `stop_retained_without_stop_codon_test_pair`
/// admits, each with a term of `NMD_GATE_TERMS`, and one CDS span, or the two of a later
/// allele of `record_rows` (`later_allele_prints_shifted_pair`), whose own pair is the
/// candidate's and whose printed pair is the span moved by the row's offset or an earlier
/// allele's, the reference's non-zero HGVS_OFFSET reaches the plugin
/// (`shifted_pair_readings`), the reference carries the value exactly when the rule of
/// NMD.pm:98-109 holds on the shifted pair and the candidate exactly when it holds on the
/// own pair: both ends of the pair are coordinates and the CDS end is at most 101 or the
/// record's genomic end passes `variant_exon_check`
/// (`TranscriptModel::ends_where_nmd_escapes`).
pub fn nmd_reads_shifted_pair(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
    sequences: &DownstreamSequences,
    record_rows: Option<&[Row]>,
) -> Res<bool> {
    let (reference, candidate) = (expected_value == NMD_ESCAPING, actual_value == NMD_ESCAPING);
    let absent = if reference {
        field_absent(a, actual_value)
    } else {
        field_absent(e, expected_value)
    };
    if reference == candidate || !absent {
        return Ok(false);
    }
    let (terms, candidate_terms) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    // One consequence set, or the pair `stop_retained_without_stop_codon_test_pair` admits.
    if terms != candidate_terms && !stop_retained_without_stop_codon_test_pair(e, a)? {
        return Ok(false);
    }
    let gated = |side: &str| side.split(',').any(|term| NMD_GATE_TERMS.contains(&term));
    if !(gated(&terms) && gated(&candidate_terms)) {
        return Ok(false);
    }
    let (own_pairs, printed) = (cds_pair_readings(a), cds_pair_readings(e));
    let Some(offset) = non_zero_offset(e)? else {
        return Ok(false);
    };
    if own_pairs.is_empty() {
        return Ok(false);
    }
    let Some(model) = sequences
        .transcript(e)?
        .and_then(|transcript| transcript.model.as_ref())
    else {
        return Ok(false);
    };
    let two_pairs = own_pairs != printed;
    if two_pairs && !later_allele_prints_shifted_pair(e, a, Some(model), offset, record_rows)? {
        return Ok(false);
    }
    let moves = moves(offset, earlier_allele_shift(e, record_rows)?);
    for (own, shifted, genomic_end) in shifted_pair_readings(e, model, offset, &own_pairs) {
        let Some((_, shifted_end)) = shifted else {
            continue;
        };
        if two_pairs
            && !moves
                .iter()
                .any(|by| printed_pair_is_moved(e, a, model, *by, Some(own)))
        {
            continue;
        }
        let escapes = model.ends_where_nmd_escapes(genomic_end);
        let on_own = matches!(own, (Some(_), Some(end)) if end <= 101 || escapes);
        if (shifted_end <= 101 || escapes) == reference && on_own == candidate {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The three classes that read the candidate cache's transcripts, on the consult `c`: an
/// NMD by `nmd_reads_shifted_pair`, a Downstream field by `downstream_class`, which names
/// the one class that takes the value and the class that must be declared beside it. The
/// values are read as the entries carry them: an empty one is the plugin's, a `-` or an
/// absent one is no value.
pub fn transcript_reading(c: &Consult) -> Res<bool> {
    let (Some(a), Some(sequences)) = (c.a, c.downstream_sequences) else {
        return Ok(false);
    };
    if field_named(c.name, &NMD_FIELDS) {
        return nmd_reads_shifted_pair(c.e, a, c.expected, c.actual, sequences, c.record_rows);
    }
    let decided = downstream_class(
        c.e,
        a,
        c.name,
        c.expected,
        c.actual,
        sequences,
        c.record_rows,
    )?;
    Ok(decided.is_some_and(|(klass, beside)| {
        klass == c.klass && beside.is_none_or(|beside| c.declared.contains(beside))
    }))
}

/// The first and last codon an entry's protein position names, the lower first: the
/// `Protein_position` read up to a `/total` suffix, or the keys `protein_start` and
/// `protein_end`; `None` where either is no number.
pub fn protein_span(e: &Entry) -> Option<(i128, i128)> {
    let (first, last) = match e.fields.get("Protein_position") {
        Some(position) => {
            let (first, _, last) = pytext::partition(pytext::partition(position, "/").0, "-");
            (first, if last.is_empty() { first } else { last })
        }
        None => (
            pytext::partition(e.fields.get_or("protein_start", ABSENT), "/").0,
            pytext::partition(e.fields.get_or("protein_end", ABSENT), "/").0,
        ),
    };
    let (first, last) = (digits_value(first)?, digits_value(last)?);
    Some((first.min(last), first.max(last)))
}

/// The first and last codon of a CDS pair, the lower first; an insertion's pair names the
/// codons of its two flanks.
pub fn codon_span(pair: (i128, i128)) -> (i128, i128) {
    // The floor of a third of the position plus two, without the sum.
    let codon = |position: i128| position.div_euclid(3) + i128::from(position.rem_euclid(3) != 0);
    (codon(pair.0.min(pair.1)), codon(pair.0.max(pair.1)))
}

/// The coding sequence's translation over its codons `codons` (first, last) under the
/// transcript's table, a stop `*`; `None` for a codon past the sequence's end or one the
/// table does not read.
pub fn codon_residues(transcript: &DownstreamTranscript, codons: (i128, i128)) -> Option<String> {
    let code = transcript.code()?;
    let mut residues = String::new();
    let mut codon = codons.0;
    while codon <= codons.1 {
        // A codon whose bases lie outside `i128` holds no base of the sequence.
        let end = codon.checked_mul(3)?;
        let letters = transcript
            .coding
            .slice(Some(end.checked_sub(3)?), Some(end));
        residues.push(codon_residue(code, letters)?);
        codon += 1;
    }
    Some(residues)
}

/// The two residues of an `Amino_acids` value that is one letter of `A` to `Z` over one.
fn peptide_pair(amino_acids: &str) -> Option<(char, char)> {
    match amino_acids.as_bytes() {
        [first @ b'A'..=b'Z', b'/', second @ b'A'..=b'Z'] => {
            Some((char::from(*first), char::from(*second)))
        }
        _ => None,
    }
}

/// `hgvs_shift_shared_across_alleles` on the differing field `name` of a later allele's
/// row: the reference entry `e` is a second or later alternate's whose first allele's row
/// carries a non-zero HGVS_OFFSET (`first_allele_shift`), the candidate entry `a` prints
/// the record's own CDS pair and the reference that span moved by the first allele's
/// offset (`shifted_pair_readings`). A CDS position is the class's there. A protein
/// position is where each entry's is the codons of its pair (`codon_span`) and the two
/// differ. An `Amino_acids` is on such a row where the two values share their alternate
/// part and each reference part is the translation of its entry's codons
/// (`codon_residues`).
pub fn later_allele_shared_shift(
    e: &Entry,
    a: &Entry,
    name: &str,
    record_rows: Option<&[Row]>,
    sequences: &DownstreamSequences,
) -> Res<bool> {
    let shift = first_allele_shift(e, record_rows)?;
    let (own_pairs, printed) = (cds_pair_readings(a), cds_pair_readings(e));
    let Some(shift) = shift else {
        return Ok(false);
    };
    if !matches!(own_pairs.first(), Some((Some(_), Some(_)))) || own_pairs == printed {
        return Ok(false);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let Some(model) = transcript.model.as_ref() else {
        return Ok(false);
    };
    for (own, moved, _) in shifted_pair_readings(e, model, shift, &own_pairs) {
        let (Some(moved), (Some(own_start), Some(own_end))) = (moved, own) else {
            continue;
        };
        let own = (own_start, own_end);
        if moved == own || !printed.contains(&(Some(moved.0), Some(moved.1))) {
            continue;
        }
        if field_named(name, &CDS_POSITION_FIELDS) {
            return Ok(true);
        }
        let (own_codons, moved_codons) = (codon_span(own), codon_span(moved));
        if own_codons == moved_codons
            || protein_span(a) != Some(own_codons)
            || protein_span(e) != Some(moved_codons)
        {
            continue;
        }
        if field_named(name, &PROTEIN_POSITION_FIELDS) {
            return Ok(true);
        }
        if transcript.table.is_none() {
            continue;
        }
        let (reference, _, changed) = pytext::partition(field_of(e, &AMINO_ACIDS_FIELDS), "/");
        let (own_reference, _, own_changed) =
            pytext::partition(field_of(a, &AMINO_ACIDS_FIELDS), "/");
        if changed.is_empty() || changed != own_changed {
            continue;
        }
        if codon_residues(transcript, moved_codons).as_deref() == Some(reference)
            && codon_residues(transcript, own_codons).as_deref() == Some(own_reference)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `blosum62_hgvs_shift_shared_across_alleles` on a differing BLOSUM62: the reference
/// entry `e` writes a value and the candidate entry `a` none, on a later allele's row whose
/// CDS position and protein position `later_allele_shared_shift` takes (the candidate
/// prints the record's own pair and the reference that span moved by the first allele's
/// offset, each protein position the codons of its pair). The reference's `Amino_acids` is
/// one residue over one (Blosum62.pm:110), the first the translation of the moved pair's
/// one codon (`codon_residues`) and the second the alternate part of the candidate's own
/// `Amino_acids`, which fails that test; and the value is the matrix cell of the
/// reference's pair.
pub fn blosum62_later_allele_shared_shift(
    e: &Entry,
    a: &Entry,
    ev: &str,
    av: &str,
    record_rows: Option<&[Row]>,
    sequences: &DownstreamSequences,
) -> Res<bool> {
    if blank(ev) || !blank(av) || !blank(field_of(a, &BLOSUM62_FIELDS)) {
        return Ok(false);
    }
    for name in [CDS_POSITION_FIELDS[0], PROTEIN_POSITION_FIELDS[0]] {
        if !later_allele_shared_shift(e, a, name, record_rows, sequences)? {
            return Ok(false);
        }
    }
    let own = field_of(a, &AMINO_ACIDS_FIELDS);
    let Some(scored) = peptide_pair(field_of(e, &AMINO_ACIDS_FIELDS)) else {
        return Ok(false);
    };
    if peptide_pair(own).is_some() || pytext::partition(own, "/").2 != scored.1.to_string() {
        return Ok(false);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let (Some(model), Some(_), Some(shift)) = (
        transcript.model.as_ref(),
        transcript.table,
        first_allele_shift(e, record_rows)?,
    ) else {
        return Ok(false);
    };
    let printed = cds_pair_readings(e);
    let own_pairs: Vec<CdsPair> = cds_pair_readings(a);
    let scored_at_a_moved_codon = shifted_pair_readings(e, model, shift, &own_pairs)
        .into_iter()
        .filter_map(|(_, moved, _)| moved)
        .filter(|moved| printed.contains(&(Some(moved.0), Some(moved.1))))
        .any(|moved| {
            codon_residues(transcript, codon_span(moved)).as_deref()
                == Some(scored.0.to_string().as_str())
        });
    if !scored_at_a_moved_codon {
        return Ok(false);
    }
    Ok(blosum62(scored.0, scored.1).is_some_and(|cell| ev == cell.to_string()))
}

/// The Blosum62 class on the consult `c`, which accepts nothing without the candidate's row
/// and the candidate cache's transcripts.
pub fn blosum62_shift_shared_across_alleles(c: &Consult, ev: &str, av: &str) -> Res<bool> {
    let (Some(a), Some(sequences)) = (c.a, c.downstream_sequences) else {
        return Ok(false);
    };
    blosum62_later_allele_shared_shift(c.e, a, ev, av, c.record_rows, sequences)
}
