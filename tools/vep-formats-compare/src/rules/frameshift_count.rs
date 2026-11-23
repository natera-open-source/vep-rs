// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a frameshift's stop count past a stop codon inside the coding sequence
//! that the transcript's peptide reads as a residue or holds as a stop and reads on past:
//! the count of the frameshift's own description, and the count of the extension that
//! names a frameshift whose first changed residue is a stop.

use crate::consts::STRAND_FIELDS;
use crate::entry::Entry;
use crate::error::Result;
use crate::genetic_codes::{genetic_code, CodonTable, GeneticCode, Reading};
use crate::pyfmt;
use crate::rules::downstream::{
    ambiguous_residue, cds_pair_readings, field_of, DownstreamSequences, DownstreamTranscript,
};
use crate::rules::hgvsp_stop::{STOP_EXTENSION_RE_PATTERN, STOP_FRAMESHIFT_RE_PATTERN};
use crate::rules::whole_text;
use crate::self_parity::reverse_complement;
use crate::shared::{blank, consequence_set};
use std::sync::LazyLock;

static FRAMESHIFT_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_FRAMESHIFT_RE_PATTERN));
static EXTENSION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_EXTENSION_RE_PATTERN));

/// The three-letter name a protein description gives the residue `letter`, read in either
/// case; `None` for a stop, for the unknown residue `X`, which a description writes as a
/// stop, and for a letter that names no residue.
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

/// A frameshift description's parts: the first changed residue's name and number, the new
/// residue's name and the count to the new stop, `?` for an unknown one. `extension` for
/// the extension that names a frameshift whose first changed residue is a stop
/// (`Ter<N><Xaa>extTer<count>`), whose count leaves that stop out.
struct Frameshift<'a> {
    residue: &'a str,
    number: &'a str,
    new_residue: &'a str,
    count: &'a str,
    extension: bool,
}

fn frameshift(description: &str) -> Option<Frameshift<'_>> {
    if let Some(parts) = FRAMESHIFT_WHOLE.captures(description) {
        let part = |at: usize| parts.get(at).map_or("", |part| part.as_str());
        return Some(Frameshift {
            residue: part(1),
            number: part(2),
            new_residue: part(3),
            count: part(4),
            extension: false,
        });
    }
    let parts = EXTENSION_WHOLE.captures(description)?;
    let part = |at: usize| parts.get(at).map_or("", |part| part.as_str());
    Some(Frameshift {
        residue: "Ter",
        number: part(1),
        new_residue: part(2),
        count: part(3),
        extension: true,
    })
}

/// The whole codons of `sequence`, each with its residue under `code`: a codon with an
/// ambiguous letter by the amino acids its readings code, as the translation of the
/// Downstream classes reads it.
fn codons(sequence: &str, code: GeneticCode) -> Vec<(&[u8], u8)> {
    sequence
        .as_bytes()
        .chunks_exact(3)
        .map(|codon| {
            let residue = code
                .amino_acid(codon)
                .unwrap_or_else(|| ambiguous_residue(code, codon));
            (codon, residue)
        })
        .collect()
}

/// The transcript's coding sequence with `allele`, on the transcript's strand, over the
/// coding pair `start` and `end` (an insertion's start one past its end) and the 3' UTR
/// after it. `None` where the pair lies outside the coding sequence and where the change
/// keeps the reading frame. Both sequences are ASCII.
fn out_of_frame_bases(
    transcript: &DownstreamTranscript,
    (start, end): (i64, i64),
    allele: &str,
) -> Option<String> {
    let coding = transcript.coding.as_str();
    let before = usize::try_from(start.saturating_sub(1)).ok()?;
    let after = usize::try_from(end.max(start.saturating_sub(1))).ok()?;
    if before > coding.len() || after > coding.len() || after < before {
        return None;
    }
    let replaced = after - before;
    if (allele.len() + 3 - replaced % 3).is_multiple_of(3) {
        return None;
    }
    let mut changed = String::with_capacity(coding.len() + allele.len() + transcript.utr.len());
    changed.push_str(&coding[..before]);
    changed.push_str(allele);
    changed.push_str(&coding[after..]);
    changed.push_str(&transcript.utr);
    Some(changed)
}

/// Whether the transcript's sequences give the count `counted` for a frameshift whose
/// first changed residue is number `first_changed`, the residue `named.0` read as
/// `named.1`, where the allele `allele`, on the transcript's strand, stands over the coding
/// pair `start` and `end` (an insertion's start one past its end):
///
/// - the transcript's peptide is the protein of its coding sequence: at each codon it has
///   the codon's residue, a stop where it holds a stop codon as one and names the residues
///   after it, or the codon is a stop codon that the peptide reads as a residue; its first
///   residue may be a methionine whatever the first codon;
/// - a stop codon of either kind stands before the first changed residue, so a translation
///   of the alternate bases ends there;
/// - the alternate coding sequence, with the 3' UTR after it, is out of frame, reads the
///   reference's codons before number `first_changed`, and at it the peptide's residue
///   `named.0`, which is no stop, against `named.1`;
/// - its first stop from that residue on, sought through the 3' UTR and no farther, stands
///   `counted` residues on, the first changed one counted, and is another codon than every
///   stop codon the peptide reads as a residue.
fn sequence_gives_the_count(
    transcript: &DownstreamTranscript,
    code: GeneticCode,
    (start, end): (i64, i64),
    allele: &str,
    first_changed: usize,
    named: (&str, &str),
    counted: &str,
) -> bool {
    let coding = transcript.coding.as_str();
    if !coding.is_ascii() || !transcript.utr.is_ascii() {
        return false;
    }
    let Some(peptide) = transcript.peptide.as_deref().map(str::as_bytes) else {
        return false;
    };
    let reference = codons(coding, code);
    let residues = peptide.len().min(reference.len());
    let Some(first) = first_changed.checked_sub(1) else {
        return false;
    };
    if first >= residues {
        return false;
    }
    let fits = (0..residues).all(|at| {
        let residue = reference[at].1;
        residue == peptide[at]
            || (residue == b'*' && peptide[at] != b'*')
            || (at == 0 && peptide[0] == b'M')
    });
    let stands_before = reference[..first]
        .iter()
        .any(|(_, residue)| *residue == b'*');
    if !fits || !stands_before {
        return false;
    }
    let read_as_residues: Vec<&[u8]> = reference[..residues]
        .iter()
        .zip(peptide)
        .filter(|((_, residue), held)| *residue == b'*' && **held != b'*')
        .map(|((codon, _), _)| *codon)
        .collect();
    let Some(changed) = out_of_frame_bases(transcript, (start, end), allele) else {
        return false;
    };
    let alternate = codons(&changed, code);
    if alternate.len() <= first
        || alternate[..first]
            .iter()
            .zip(&reference[..first])
            .any(|(alternate, reference)| alternate.1 != reference.1)
        || residue_name(peptide[first]) != Some(named.0)
        || residue_name(alternate[first].1) != Some(named.1)
    {
        return false;
    }
    let Some(to_the_stop) = alternate[first..]
        .iter()
        .position(|(_, residue)| *residue == b'*')
    else {
        return false;
    };
    let new_stop = alternate[first + to_the_stop].0;
    (to_the_stop + 1).to_string() == counted
        && !read_as_residues
            .iter()
            .any(|codon| codon.eq_ignore_ascii_case(new_stop))
}

/// Whether the transcript's sequences give the count `counted` for the extension that
/// names a frameshift whose first changed residue is the stop of residue number `lost`,
/// read as `new_residue`, where the allele `allele`, on the transcript's strand, stands
/// over the coding pair `start` and `end` (an insertion's start one past its end):
///
/// - the transcript's peptide holds a stop and is the protein of its coding sequence: at
///   each codon it has the codon's residue, its first residue a methionine whatever the
///   first codon;
/// - residue `lost` is a stop the peptide holds, or the stop codon that closes the coding
///   sequence, one past the peptide's last residue;
/// - a stop codon stands before it, so a translation of the alternate bases ends there;
/// - the alternate coding sequence, with the 3' UTR after it, is out of frame, reads the
///   reference's codons before number `lost`, and at it `new_residue`, which is no stop;
/// - its first stop after that residue, sought through the 3' UTR and no farther, stands
///   `counted` residues on from the lost stop.
fn sequence_gives_the_extension_count(
    transcript: &DownstreamTranscript,
    code: GeneticCode,
    pair: (i64, i64),
    allele: &str,
    lost: usize,
    new_residue: &str,
    counted: &str,
) -> bool {
    let coding = transcript.coding.as_str();
    if !coding.is_ascii() || !transcript.utr.is_ascii() {
        return false;
    }
    let Some(peptide) = transcript.peptide.as_deref().map(str::as_bytes) else {
        return false;
    };
    let Some(at) = lost.checked_sub(1) else {
        return false;
    };
    let reference = codons(coding, code);
    let residues = peptide.len().min(reference.len());
    let fits = (0..residues)
        .all(|place| reference[place].1 == peptide[place] || (place == 0 && peptide[0] == b'M'));
    let held = at < residues && peptide[at] == b'*';
    let closing = at == peptide.len()
        && coding.len() == 3 * (at + 1)
        && reference
            .get(at)
            .is_some_and(|(_, residue)| *residue == b'*');
    let stands_before = reference
        .iter()
        .take(at)
        .any(|(_, residue)| *residue == b'*');
    if !peptide.contains(&b'*') || !fits || !(held || closing) || !stands_before {
        return false;
    }
    let Some(changed) = out_of_frame_bases(transcript, pair, allele) else {
        return false;
    };
    let alternate = codons(&changed, code);
    if alternate.len() <= at
        || alternate[..at]
            .iter()
            .zip(&reference[..at])
            .any(|(alternate, reference)| alternate.1 != reference.1)
        || residue_name(alternate[at].1) != Some(new_residue)
    {
        return false;
    }
    alternate[at..]
        .iter()
        .position(|(_, residue)| *residue == b'*')
        .is_some_and(|residues_on| residues_on.to_string() == counted)
}

/// `hgvsp_frameshift_count_past_a_stop_read_as_a_residue` on a differing HGVSp: both
/// values name one protein; the reference's description, a pair of prediction parentheses
/// set aside, is a frameshift with the unknown count and the candidate's the same
/// frameshift, at the same residue and with the same new residue, with a count; and the
/// row's transcript in the candidate's cache (`sequences`) gives that count for the
/// candidate row's allele at one reading of its coding position
/// (`sequence_gives_the_count`). A pair of extensions named at one stop with one new
/// residue, the reference's of unknown count and the candidate's with a count, is read
/// the same way on a reference row whose terms hold `frameshift_variant`
/// (`sequence_gives_the_extension_count`). The cache is read for a pair of either shape
/// alone.
pub fn frameshift_count_past_a_stop_read_as_a_residue(
    e: &Entry,
    a: &Entry,
    ev: &str,
    av: &str,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    if blank(ev) || blank(av) {
        return Ok(false);
    }
    let (reference_name, separator, reference) = pyfmt::partition(ev, ":p.");
    let (candidate_name, _, candidate) = pyfmt::partition(av, ":p.");
    if separator.is_empty() || reference_name != candidate_name {
        return Ok(false);
    }
    let predicted = |description: &str| description.starts_with('(') && description.ends_with(')');
    let (mut reference, mut candidate) = (reference, candidate);
    if predicted(reference) && predicted(candidate) {
        reference = &reference[1..reference.len() - 1];
        candidate = &candidate[1..candidate.len() - 1];
    }
    let (Some(unknown), Some(counted)) = (frameshift(reference), frameshift(candidate)) else {
        return Ok(false);
    };
    if unknown.count != "?"
        || counted.count == "?"
        || (
            unknown.extension,
            unknown.residue,
            unknown.number,
            unknown.new_residue,
        ) != (
            counted.extension,
            counted.residue,
            counted.number,
            counted.new_residue,
        )
    {
        return Ok(false);
    }
    let frameshift_row = consequence_set(&e.consequence)
        .split(',')
        .any(|term| term == "frameshift_variant");
    if counted.extension && (counted.new_residue == "Ter" || !frameshift_row) {
        return Ok(false);
    }
    let Ok(first_changed) = counted.number.parse::<usize>() else {
        return Ok(false);
    };
    let allele = if a.allele == "-" {
        String::new()
    } else {
        pyfmt::upper(&a.allele)
    };
    if !allele.bytes().all(|base| b"ACGT".contains(&base)) {
        return Ok(false);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let CodonTable::Held(table) = &transcript.table else {
        return Ok(false);
    };
    let Some(code) = genetic_code(*table, Reading::VepRs) else {
        return Ok(false);
    };
    let strand = match transcript.model.as_ref() {
        Some(model) => model.strand,
        None if field_of(Some(e), STRAND_FIELDS) == "-1" => -1,
        None => 1,
    };
    let allele = if strand == -1 {
        reverse_complement(&allele)
    } else {
        allele
    };
    Ok(cds_pair_readings(a).into_iter().any(|pair| {
        let (Some(start), Some(end)) = pair else {
            return false;
        };
        if counted.extension {
            return sequence_gives_the_extension_count(
                &transcript,
                code,
                (start, end),
                &allele,
                first_changed,
                counted.new_residue,
                counted.count,
            );
        }
        sequence_gives_the_count(
            &transcript,
            code,
            (start, end),
            &allele,
            first_changed,
            (counted.residue, counted.new_residue),
            counted.count,
        )
    }))
}
