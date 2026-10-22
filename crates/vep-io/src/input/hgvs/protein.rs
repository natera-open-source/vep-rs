// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein branch (`p.`) of Ensembl VEP's `fetch_by_hgvs_notation`
//! (ensembl-variation release/116, VariationFeatureAdaptor.pm:1900-1986): the
//! transcript whose translation the notation names is found, the described
//! residue change is turned into the nucleotide change that reaches it in the
//! fewest steps (`_parse_hgvs_protein_position`, `:2232-2319`), or a deletion
//! of its codons (`_parse_hgvs_protein_position_del`, `:2323-2354`), and the
//! variant is built from that change as the other branches build theirs
//! (`components::from_alleles`), on the transcript's strand.
//!
//! A description of a frameshift (letters, digits, letters and `fs`, `:1903`)
//! is refused before anything is looked up. The transcript is
//! `fetch_by_translation_stable_id` on the reference, its version removed from
//! an Ensembl identifier (`:1909`); a reference naming a transcript (`ENST`)
//! is read as that transcript (`:1922-1925`); any other is read as a gene
//! name, and every transcript of the gene is tried in database order, the
//! first to yield a variant taken (`_get_gene_transcripts` with `multiple_ok`,
//! `:1933`, `:1946-1985`). Ensembl VEP also reads an `NP_`, `XP_` or `NM_`
//! reference and a UniProt accession through its external-name tables
//! (`:1915-1916`, `:1929`); here a RefSeq protein or transcript resolves
//! where the cache carries it (a RefSeq cache, by translation or transcript
//! id) and no other external name is read, so a UniProt accession is read
//! as a gene name and refused.
//!
//! A substitution is read by the regular expression at `:2236`: an optional
//! residue, the position, an optional `delins` and the new residue (a word,
//! `*`, `=` or `?`); `=` names the written residue again (`:2237`) and a
//! `delins` without a written residue takes the reference's (`:2250`).
//! Three-letter codes are read case-sensitively (`Arg`, not `ARG`, `:2256-
//! 2257`). The codon is read from the genome through the transcript's mapper
//! (`pep2genomic`, TranscriptMapper.pm:372-388, the coding start less the
//! start phase), on the transcript's strand, and must translate to the written
//! residue (`:2262-2267`); the new residue's codons are listed (`revtranslate`)
//! and, for each, the codon positions that change are the path to it, a change
//! of the first and third base alone counting the middle base too
//! (`:2274-2292`). The shortest paths are kept (`:2295-2296`); when one is
//! left, the variant spans its changed bases, the reference allele read from
//! the genome over that span and the alternate the new bases (`:2302-2317`);
//! when several are, the notation is refused as one whose nucleotide change
//! cannot be uniquely determined (`:1956`), since `--ambiguous_hgvs`, which
//! would list them all, is not accepted. A deletion of one residue or a range
//! of them takes the codons' bases as its reference allele and the first exon
//! piece's span as its position (`:2351`), so a deletion whose codons cross an
//! exon boundary is refused at the reference check (`:2100`), as it is in
//! Ensembl VEP.
//!
//! Refusals carry Ensembl VEP's texts. A position the mapper places in no
//! exon, or a range written end first, is refused as one whose genomic
//! reference sequence cannot be found (`:2242`), where Ensembl VEP's resolver
//! dies in `get_cds_reference` on the gap; a description the regular
//! expressions leave without a position is refused with `pep2genomic`'s
//! `Must call with start and end` (TranscriptMapper.pm:376), the error Ensembl
//! VEP's resolver throws for it. When every transcript of a gene refuses a
//! notation, the distinct refusals are joined on one line in the order met,
//! where Ensembl VEP joins them in hash order (`:1983`).

use std::collections::BTreeSet;

use vep_core::codon::translate_codon_with_table;
use vep_core::coordinate::Strand;
use vep_core::transcript::{MappedPiece, Transcript};
use vep_core::variant::InputVariant;

use crate::error::IoError;
use crate::input::{ReferenceSequence, TranscriptSource};

use super::components::{self, genome_on_strand, Placement, ReadFlags};
use super::transcript::mapper_of;
use super::{cache_unreadable, unable, Notation};

/// The variant a protein notation describes, read against `reference` on the
/// transcript `transcripts` resolves its reference to, or on the first
/// transcript of the gene so named that yields one, its alleles trimmed as
/// written under `flags.minimal`; refused with Ensembl VEP's message when none
/// does.
pub(super) fn variant(
    hgvs: &str,
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    transcripts: &dyn TranscriptSource,
    flags: ReadFlags,
) -> Result<InputVariant, IoError> {
    let refuse = |message: String| unable(hgvs, Some(&message));
    if is_frameshift(&notation.description) {
        return Err(refuse(
            "Frameshifts are not supported for HGVS protein input".to_string(),
        ));
    }
    let not_found = || {
        format!(
            "Could not get a Transcript object for '{}'",
            notation.reference
        )
    };
    let mut first_refusal = None;
    for group in resolve(&notation.reference, transcripts)? {
        if group.is_empty() {
            first_refusal.get_or_insert_with(not_found);
            continue;
        }
        match on_transcripts(hgvs, notation, reference, &group, flags) {
            Ok(variant) => return Ok(variant),
            Err(refusal) => {
                first_refusal.get_or_insert(refusal);
            }
        }
    }
    Err(refuse(first_refusal.unwrap_or_else(not_found)))
}

/// The transcripts `reference` names, resolved as Ensembl VEP's protein branch
/// resolves them (`:1908-1941`) with one lookup moved: the transcript an id
/// finds, in a list of its own, else one list per core group of every
/// transcript of the gene so named, in database order (`:1933`). A reference
/// naming a transcript (`ENST`, or `NM`, whose external-name lookup in Ensembl
/// VEP's RefSeq database finds the RefSeq transcript itself, `:1915-1916`)
/// goes to the transcript lookup first (`:1922-1925`), ahead of the
/// translation lookup Ensembl VEP makes first (`:1909`), since a translation
/// rarely carries such an id and a miss there would read the cache's
/// alternate tier for nothing; a reference the transcript lookup misses goes
/// on to the translation lookup (a translation id whose species code begins
/// with `T` contains `ENST`), then an `ENST` reference is one empty list, and
/// an `NM` one goes on to the gene lookup, which follows a missed
/// external-name lookup (`:1929-1933`), so a gene symbol containing `NM` is
/// read as a gene name. Any other reference goes to the translation lookup
/// (`:1909`; `XP` and `NP` are translation ids, and an `XM` reference is a
/// gene name, as the test at `:1915` has no `XM`), then to the gene lookup. An
/// empty list is Ensembl VEP's `Could not get a Transcript object`. A cache
/// that cannot be read stops the run.
fn resolve(
    reference: &str,
    transcripts: &dyn TranscriptSource,
) -> Result<Vec<Vec<Transcript>>, IoError> {
    if reference.contains("ENST") || reference.contains("NM") {
        if let Some(transcript) = transcripts
            .by_stable_id(reference)
            .map_err(cache_unreadable)?
        {
            return Ok(vec![vec![transcript]]);
        }
    }
    if let Some(transcript) = transcripts
        .by_translation_id(reference)
        .map_err(cache_unreadable)?
    {
        return Ok(vec![vec![transcript]]);
    }
    if reference.contains("ENST") {
        return Ok(vec![Vec::new()]);
    }
    transcripts
        .gene_transcripts(reference)
        .map_err(cache_unreadable)
}

/// The variant the first of `transcripts` to yield one gives (`:1946-1985`);
/// when none does, every distinct refusal met, joined on one line.
fn on_transcripts(
    hgvs: &str,
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    transcripts: &[Transcript],
    flags: ReadFlags,
) -> Result<InputVariant, String> {
    let mut refusals: Vec<String> = Vec::new();
    let mut refuse = |message: String| {
        if !refusals.contains(&message) {
            refusals.push(message);
        }
    };
    for transcript in transcripts {
        let changes = match changes_of(notation, reference, transcript) {
            Ok(changes) => changes,
            Err(message) => {
                refuse(message);
                continue;
            }
        };
        if changes.len() > 1 {
            refuse(format!(
                "Could not uniquely determine nucleotide change from {}",
                notation.cleaned
            ));
            continue;
        }
        for change in changes {
            let placement = Placement {
                chr: transcript.chr.clone(),
                start: i128::from(change.start),
                end: i128::from(change.end),
                strand: change.strand,
            };
            match components::from_alleles(
                hgvs,
                notation,
                reference,
                placement,
                change.reference_allele,
                change.alt_allele,
                flags,
            ) {
                Ok(variant) => return Ok(variant),
                Err(message) => refuse(message),
            }
        }
    }
    Err(refusals.join("; "))
}

/// A nucleotide change a description resolves to on one transcript, as
/// `_parse_hgvs_protein_position` lists it: the alleles on the transcript's
/// strand and the genomic span they occupy.
struct Change {
    reference_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    start: i64,
    end: i64,
    strand: Strand,
}

/// The changes `notation`'s description resolves to on `transcript`: the
/// deletion of its codons for a description ending in `del` (`:1950-1952`),
/// else the shortest nucleotide paths to the written residue.
fn changes_of(
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    transcript: &Transcript,
) -> Result<Vec<Change>, String> {
    let description = notation.description.as_str();
    if description.ends_with("del") {
        deletion(description, &notation.reference, reference, transcript)
    } else {
        substitution(description, &notation.reference, reference, transcript)
    }
}

/// `_parse_hgvs_protein_position` (`:2232-2319`).
fn substitution(
    description: &str,
    reference_name: &str,
    reference: &dyn ReferenceSequence,
    transcript: &Transcript,
) -> Result<Vec<Change>, String> {
    let Some((from, position, to)) = substitution_captures(description) else {
        return Err("Must call with start and end".to_string());
    };
    let mut from = from.map(str::to_string);
    let to = if to == "=" {
        from.clone()
    } else {
        Some(to.to_string())
    };
    let codon = CodonBases::read(reference, reference_name, transcript, position, None)?;
    let table = codon_table(transcript);
    let translated = translate(&codon.bases, table);
    if from.is_none() && description.contains("delins") {
        from = Some(translated.clone());
    }
    // Perl's truth test on the captures: a residue is a word, the position
    // fails as the string `0` alone.
    let truthy = |text: &str| !text.is_empty() && text != "0";
    let (from, to) = match (&from, &to) {
        (Some(from), Some(to)) if truthy(from) && truthy(position) && truthy(to) => {
            (one_letter(from), one_letter(to))
        }
        _ => {
            return Err(format!(
                "Could not parse HGVS protein notation {reference_name}:p.{description}"
            ))
        }
    };
    if translated != from {
        return Err(format!(
            "Sequence translated from reference ({} -> {translated}) does not match input sequence ({from})",
            String::from_utf8_lossy(&codon.bases)
        ));
    }
    let paths = shortest_paths(&codon.bases, &reverse_translate(&to, table));
    if paths.is_empty() {
        return Err(format!(
            "Could not determine nucleotide change from peptide change {from} -> {to}"
        ));
    }
    let mut positions: Vec<i64> = codon.spans.iter().flat_map(|&(s, e)| s..=e).collect();
    positions.sort_unstable();
    positions.dedup();
    if codon.strand == Strand::Reverse {
        positions.reverse();
    }
    paths
        .into_iter()
        .map(|path| {
            let changed: Vec<i64> = path.iter().map(|step| positions[step.index]).collect();
            let (start, end) = (
                *changed.iter().min().unwrap(),
                *changed.iter().max().unwrap(),
            );
            let reference_allele = genome_on_strand(
                reference,
                reference_name,
                &transcript.chr,
                i128::from(start),
                i128::from(end),
                codon.strand,
            )?;
            Ok(Change {
                reference_allele,
                alt_allele: path.iter().map(|step| step.to).collect(),
                start,
                end,
                strand: codon.strand,
            })
        })
        .collect()
}

/// `_parse_hgvs_protein_position_del` (`:2323-2354`): the codons of the
/// residue, or of the range `first_second`, deleted; the span is the first
/// exon piece's (`:2351`).
fn deletion(
    description: &str,
    reference_name: &str,
    reference: &dyn ReferenceSequence,
    transcript: &Transcript,
) -> Result<Vec<Change>, String> {
    let (position, end_position) = if description.contains('_') {
        let mut parts = description.split('_');
        let first = parts.next().unwrap_or_default();
        let second = parts.next().unwrap_or_default();
        (
            deletion_first_captures(first).map(|(_, position)| position),
            deletion_captures(second).map(|(_, position, _)| position),
        )
    } else {
        (
            deletion_captures(description).map(|(_, position, _)| position),
            None,
        )
    };
    let position = position.ok_or_else(|| "Must call with start and end".to_string())?;
    let codons = CodonBases::read(
        reference,
        reference_name,
        transcript,
        position,
        end_position,
    )?;
    let (start, end) = codons.spans[0];
    Ok(vec![Change {
        reference_allele: codons.bases,
        alt_allele: b"-".to_vec(),
        start,
        end,
        strand: codons.strand,
    }])
}

/// `get_cds_reference` (`:2357-2379`): the bases of the codons from one
/// residue position to another, read on the transcript's strand from each
/// exon piece the mapper places them on, with those pieces' spans.
struct CodonBases {
    bases: Vec<u8>,
    spans: Vec<(i64, i64)>,
    strand: Strand,
}

impl CodonBases {
    fn read(
        reference: &dyn ReferenceSequence,
        reference_name: &str,
        transcript: &Transcript,
        position: &str,
        end_position: Option<&str>,
    ) -> Result<Self, String> {
        let unable_to_find = || {
            format!("Unable to find the genomic reference sequence for protein {reference_name}")
        };
        let pieces = codon_pieces(transcript, position, end_position.unwrap_or(position))
            .ok_or_else(unable_to_find)?;
        let mut bases = Vec::new();
        let mut spans = Vec::new();
        let mut strand = None;
        for piece in pieces {
            let MappedPiece::Coordinate {
                start,
                end,
                strand: piece_strand,
            } = piece
            else {
                return Err(unable_to_find());
            };
            if start > end {
                return Err(unable_to_find());
            }
            bases.extend(genome_on_strand(
                reference,
                reference_name,
                &transcript.chr,
                i128::from(start),
                i128::from(end),
                piece_strand,
            )?);
            spans.push((start, end));
            strand.get_or_insert(piece_strand);
        }
        Ok(Self {
            bases,
            spans,
            strand: strand.ok_or_else(unable_to_find)?,
        })
    }
}

/// `pep2genomic` (TranscriptMapper.pm:372-388): the genomic pieces of the
/// codons from residue `start` to residue `end`, each the digits of the
/// description read as a number; `None` for a number the arithmetic cannot
/// hold.
fn codon_pieces(transcript: &Transcript, start: &str, end: &str) -> Option<Vec<MappedPiece>> {
    let mapper = mapper_of(transcript);
    let shift = i64::from(mapper.start_phase.max(0));
    let offset = i64::try_from(mapper.cdna_coding_start).ok()? - 1 - shift;
    let cdna_start = start
        .parse::<i64>()
        .ok()?
        .checked_mul(3)?
        .checked_sub(2)?
        .checked_add(offset)?;
    let cdna_end = end
        .parse::<i64>()
        .ok()?
        .checked_mul(3)?
        .checked_add(offset)?;
    Some(mapper.exon_coord_mapper.cdna2genomic(cdna_start, cdna_end))
}

/// One changed codon base of a path (`:2282`): its index in the codon and
/// the base before and after.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Step {
    index: usize,
    from: u8,
    to: u8,
}

/// The paths from `codon` to each of `to_codons` with the fewest steps
/// (`:2272-2296`): a step per base that differs, every base when the codons
/// are the same, and the middle base added between a changed first and third
/// (`:2286-2288`); each distinct shortest path once.
fn shortest_paths(codon: &[u8], to_codons: &[[u8; 3]]) -> Vec<Vec<Step>> {
    let mut paths: Vec<Vec<Step>> = Vec::new();
    for to_codon in to_codons {
        let same_codon = codon.eq_ignore_ascii_case(to_codon);
        let mut steps: Vec<Step> = codon
            .iter()
            .zip(to_codon)
            .enumerate()
            .filter(|(_, (from, to))| same_codon || !from.eq_ignore_ascii_case(to))
            .map(|(index, (from, to))| Step {
                index,
                from: from.to_ascii_uppercase(),
                to: to.to_ascii_uppercase(),
            })
            .collect();
        if steps.len() == 2 && steps[0].index == 0 && steps[1].index == 2 {
            steps.insert(
                1,
                Step {
                    index: 1,
                    from: codon[1],
                    to: codon[1],
                },
            );
        }
        if !paths.contains(&steps) {
            paths.push(steps);
        }
    }
    let Some(shortest) = paths.iter().map(Vec::len).min() else {
        return paths;
    };
    paths.retain(|path| path.len() == shortest);
    paths
}

/// `/[A-Za-z]+[0-9]+[A-Za-z]+fs/` (`:1903`): some `fs` directly after letters,
/// digits and letters.
fn is_frameshift(description: &str) -> bool {
    let bytes = description.as_bytes();
    let run_before = |end: usize, is_in: fn(&u8) -> bool| {
        bytes[..end].iter().rev().take_while(|b| is_in(b)).count()
    };
    description.match_indices("fs").any(|(at, _)| {
        let letters = run_before(at, u8::is_ascii_alphabetic);
        let digits = run_before(at - letters, u8::is_ascii_digit);
        letters > 0 && digits > 0 && run_before(at - letters - digits, u8::is_ascii_alphabetic) > 0
    })
}

/// `^([A-Za-z]+?)?(\d+)(?:delins)?(\w+|\*|\=|\?)$` (`:2236`), as Perl's engine
/// matches it: the residue before the position is tried from one letter up
/// and then left out, the position from all its digits down, a `delins` taken
/// before being left, and the rest must be a word, `*`, `=` or `?` to the end.
/// The captures: the written residue, the position's digits, the new residue.
fn substitution_captures(description: &str) -> Option<(Option<&str>, &str, &str)> {
    let letters = description
        .bytes()
        .take_while(u8::is_ascii_alphabetic)
        .count();
    let residue_lengths = (1..=letters).map(Some).chain(std::iter::once(None));
    for residue_len in residue_lengths {
        let rest = &description[residue_len.unwrap_or(0)..];
        let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
        for taken in (1..=digits).rev() {
            let after_digits = &rest[taken..];
            for tail in [after_digits.strip_prefix("delins"), Some(after_digits)]
                .into_iter()
                .flatten()
            {
                let is_word_run = !tail.is_empty() && tail.bytes().all(|b| is_word(&b));
                if is_word_run || matches!(tail, "*" | "=" | "?") {
                    return Some((
                        residue_len.map(|len| &description[..len]),
                        &rest[..taken],
                        tail,
                    ));
                }
            }
        }
    }
    None
}

/// `^(\w+?)(\d+)(\w+?|\*)$` (`:2337`, `:2342`): the shortest word run, then
/// the longest digits that leave a word run or a `*` to the end.
fn deletion_captures(text: &str) -> Option<(&str, &str, &str)> {
    let words = text.bytes().take_while(is_word).count();
    for residue_len in 1..=words {
        let rest = &text[residue_len..];
        let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
        for taken in (1..=digits).rev() {
            let tail = &rest[taken..];
            if (!tail.is_empty() && tail.bytes().all(|b| is_word(&b))) || tail == "*" {
                return Some((&text[..residue_len], &rest[..taken], tail));
            }
        }
    }
    None
}

/// `^(\w+?)(\d+)$` (`:2334`): the shortest word run followed by digits to the
/// end.
fn deletion_first_captures(text: &str) -> Option<(&str, &str)> {
    let words = text.bytes().take_while(is_word).count();
    (1..=words).find_map(|residue_len| {
        let rest = &text[residue_len..];
        (!rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()))
            .then(|| (&text[..residue_len], rest))
    })
}

/// Perl's `\w`: an ASCII letter, digit or underscore.
fn is_word(byte: &u8) -> bool {
    byte.is_ascii_alphanumeric() || *byte == b'_'
}

/// `$Bio::SeqUtils::ONECODE{$code} || $code` (`:2256-2257`): the one-letter
/// code of a three-letter one written as BioPerl spells it (`Arg`, `Ter`,
/// `Sec`), any other text as written.
fn one_letter(code: &str) -> String {
    let one = match code {
        "Ala" => "A",
        "Arg" => "R",
        "Asn" => "N",
        "Asp" => "D",
        "Asx" => "B",
        "Cys" => "C",
        "Gln" => "Q",
        "Glu" => "E",
        "Glx" => "Z",
        "Gly" => "G",
        "His" => "H",
        "Ile" => "I",
        "Leu" => "L",
        "Lys" => "K",
        "Met" => "M",
        "Phe" => "F",
        "Pro" => "P",
        "Pyl" => "O",
        "Sec" => "U",
        "Ser" => "S",
        "Ter" => "*",
        "Thr" => "T",
        "Trp" => "W",
        "Tyr" => "Y",
        "Val" => "V",
        "Xaa" => "X",
        "Xle" => "J",
        other => other,
    };
    one.to_string()
}

/// The codon table of the transcript's sequence region (`:2245-2248`): the
/// cache's, read as the standard table when absent or zero.
fn codon_table(transcript: &Transcript) -> u8 {
    transcript.vefc.as_ref().map_or(1, |vefc| {
        if vefc.codon_table == 0 {
            1
        } else {
            vefc.codon_table
        }
    })
}

/// The nucleotides an IUPAC base stands for, as BioPerl's `iupac_iub` lists
/// them; `None` for a byte that is no base.
fn iupac_bases(base: u8) -> Option<&'static [u8]> {
    Some(match base.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'M' => b"AC",
        b'R' => b"AG",
        b'W' => b"AT",
        b'S' => b"CG",
        b'Y' => b"CT",
        b'K' => b"GT",
        b'V' => b"ACG",
        b'H' => b"ACT",
        b'D' => b"AGT",
        b'B' => b"CGT",
        b'X' | b'N' => b"GATC",
        _ => return None,
    })
}

/// `Bio::Tools::CodonTable::translate` of one codon under `table`: the amino
/// acid the codon codes; for a codon with an ambiguous base, the one amino acid
/// every reading gives, `B` for aspartate and asparagine, `Z` for glutamate and
/// glutamine, else `X`; an empty text for a codon not three bases long.
fn translate(codon: &[u8], table: u8) -> String {
    if codon.len() != 3 {
        return String::new();
    }
    let readings: Option<Vec<&[u8]>> = codon.iter().map(|&base| iupac_bases(base)).collect();
    let Some(readings) = readings else {
        return "X".to_string();
    };
    let mut amino_acids = BTreeSet::new();
    for &first in readings[0] {
        for &second in readings[1] {
            for &third in readings[2] {
                amino_acids.insert(translate_codon_with_table(&[first, second, third], table));
            }
        }
    }
    let translated = match amino_acids.len() {
        1 => amino_acids.into_iter().next().unwrap_or(b'X'),
        2 if amino_acids == BTreeSet::from(*b"DN") => b'B',
        2 if amino_acids == BTreeSet::from(*b"EQ") => b'Z',
        _ => b'X',
    };
    (translated as char).to_string()
}

/// `Bio::Tools::CodonTable::revtranslate` under `table`: the codons coding
/// `residue`, a one-letter code in upper case or a three-letter code in any
/// case, in the table's order (`T`, `C`, `A`, `G` at each position); an
/// ambiguity code lists the codons of each amino acid it stands for; any other
/// text, and a residue no codon codes, lists none.
fn reverse_translate(residue: &str, table: u8) -> Vec<[u8; 3]> {
    let one = if residue.len() == 3 && residue.is_ascii() {
        let mut spelled = residue.to_ascii_lowercase();
        spelled[..1].make_ascii_uppercase();
        match one_letter(&spelled) {
            one if one.len() == 1 => one,
            _ => return Vec::new(),
        }
    } else {
        residue.to_string()
    };
    let [one] = one.as_bytes() else {
        return Vec::new();
    };
    if !b"ACDEFGHIJKLMNOPQRSTUVWYBZX*".contains(one) {
        return Vec::new();
    }
    let amino_acids: &[u8] = match one {
        b'B' => b"DN",
        b'Z' => b"EQ",
        b'J' => b"IL",
        _ => std::slice::from_ref(one),
    };
    let mut codons = Vec::new();
    for amino_acid in amino_acids {
        for first in b"TCAG" {
            for second in b"TCAG" {
                for third in b"TCAG" {
                    let codon = [*first, *second, *third];
                    if translate_codon_with_table(&codon, table) == *amino_acid {
                        codons.push(codon);
                    }
                }
            }
        }
    }
    codons
}

#[cfg(test)]
mod tests {
    use super::super::transcript::tests::{
        base, complement, context, context_bounded, forward, non_coding, read_with, reverse, Cache,
        Genome,
    };
    use super::super::HgvsParser;
    use super::*;
    use crate::input::{InputParser, ParserContext};
    use std::io::Cursor;
    use std::sync::{Arc, Mutex};

    /// On the `ACGT` chromosome the forward transcript's codons cycle `ACG`
    /// (Thr), `TAC` (Tyr), `GTA` (Val), `CGT` (Arg) from genomic 121, and the
    /// reverse transcript's `GTA` (Val), `CGT` (Arg), `ACG` (Thr), `TAC` (Tyr)
    /// from genomic 550 downward on its strand; the forward transcript's
    /// translation is `ENSP00000000001`, the reverse one's `ENSP00000000002`.
    fn read(data: &str) -> Vec<Result<InputVariant, String>> {
        read_with(vec![forward(), reverse(), non_coding()], false, data)
    }

    /// `(start, end, strand, allele string)` of one parsed notation.
    fn placed(data: &str) -> (u64, u64, Strand, String) {
        let v = read(data).remove(0).unwrap();
        (v.start, v.end, v.strand, v.allele_string)
    }

    fn refused(data: &str) -> String {
        read(data).remove(0).unwrap_err()
    }

    fn unable(hgvs: &str, message: &str) -> String {
        format!("Unable to parse HGVS notation '{hgvs}': {message}")
    }

    /// A substitution resolves to the one nucleotide path with the fewest
    /// changes from the genome's codon to a codon of the new residue
    /// (`:2272-2317`): `Thr1Ala` changes the first base alone (`ACG` to `GCG`,
    /// the three other alanine codons two or three bases away); `Thr1Trp`
    /// changes two consecutive bases to tryptophan's one codon; `Tyr2Met`
    /// changes all three to methionine's; `Arg4Trp` changes the first and
    /// third, so the unchanged middle base is carried in both alleles
    /// (`:2286-2288`); `Tyr2=` is the synonymous change to the one tyrosine
    /// codon a base away (`:2237`); `Thr1Ter` the two-base change to `TAG`,
    /// the other stops three bases away. The record is named by the notation,
    /// keeps the transcript's strand, and resolves through the translation id
    /// with or without a version, through the transcript id written with
    /// `p.`, with one-letter codes, and with a lower-case three-letter new
    /// residue, which `revtranslate` reads in any case where the written
    /// residue is read as spelled. On the reverse transcript `Val1Ala` changes
    /// the middle base, genomic 549, and the alleles are on the transcript's
    /// strand.
    #[test]
    fn substitution_unique_codon_path() {
        let v = read("ENSP00000000001:p.Thr1Ala").remove(0).unwrap();
        assert_eq!(
            (
                v.start,
                v.end,
                v.strand,
                v.chr.as_str(),
                v.allele_string.as_str()
            ),
            (121, 121, Strand::Forward, "T", "A/G")
        );
        assert_eq!(v.id.as_deref(), Some("ENSP00000000001:p.Thr1Ala"));
        assert_eq!(v.raw_input.as_deref(), Some("ENSP00000000001:p.Thr1Ala"));
        for notation in [
            "ENSP00000000001.3:p.Thr1Ala",
            "ENST00000000001:p.Thr1Ala",
            "ENSP00000000001:p.T1A",
            "ENSP00000000001:p.Thr1ala",
        ] {
            assert_eq!(
                placed(notation),
                (121, 121, Strand::Forward, "A/G".to_string()),
                "{notation}"
            );
        }
        assert_eq!(
            placed("ENSP00000000001:p.Thr1Trp"),
            (121, 122, Strand::Forward, "AC/TG".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.Tyr2Met"),
            (124, 126, Strand::Forward, "TAC/ATG".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.Arg4Trp"),
            (130, 132, Strand::Forward, "CGT/TGG".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.Tyr2="),
            (126, 126, Strand::Forward, "C/T".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.Thr1Ter"),
            (121, 122, Strand::Forward, "AC/TA".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.T1*"),
            placed("ENSP00000000001:p.Thr1Ter")
        );
        let v = read("ENSP00000000002:p.Val1Ala").remove(0).unwrap();
        assert_eq!(
            (v.start, v.end, v.strand, v.allele_string.as_str()),
            (549, 549, Strand::Reverse, "T/C")
        );
        assert_eq!(complement(base(549)), b'T');
    }

    /// A residue change several codon paths reach in the same number of
    /// steps is refused with the resolver's message (`:1956`), `--ambiguous_hgvs`
    /// not being accepted: `Thr1=` has three threonine codons a base away,
    /// `Tyr2Ter` two stop codons a base away, `Thr1Phe` two phenylalanine
    /// codons three bases away.
    #[test]
    fn substitution_ambiguous_paths_refused_as_perl_without_flag() {
        for notation in [
            "ENSP00000000001:p.Thr1=",
            "ENSP00000000001:p.Tyr2Ter",
            "ENSP00000000001:p.Thr1Phe",
        ] {
            assert_eq!(
                refused(notation),
                unable(
                    notation,
                    &format!("Could not uniquely determine nucleotide change from {notation}")
                )
            );
        }
    }

    /// A `delins` resolves as the substitution does, through the trim of
    /// `_hgvs_from_components` (`:2079-2082`), the written residue checked when
    /// given and read from the genome when not (`:2250`). A `del` deletes the
    /// codons of the residue or of the range, written residues unchecked
    /// (`:2323-2354`), on either strand; a range whose codons cross an exon
    /// boundary takes the first exon piece's span and the bases of both, so the
    /// reference check refuses it (`:2351`, `:2100`); a range written end first
    /// maps to no bases and is refused as a sequence that cannot be found.
    #[test]
    fn delins_and_del() {
        assert_eq!(
            placed("ENSP00000000001:p.Thr1delinsTrp"),
            (121, 122, Strand::Forward, "AC/TG".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.1delinsTrp"),
            placed("ENSP00000000001:p.Thr1delinsTrp")
        );
        assert_eq!(
            placed("ENSP00000000001:p.Thr1del"),
            (121, 123, Strand::Forward, "ACG/-".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.Thr1_Tyr2del"),
            (121, 126, Strand::Forward, "ACGTAC/-".to_string())
        );
        assert_eq!(
            placed("ENSP00000000001:p.T1_Y2del"),
            placed("ENSP00000000001:p.Thr1_Tyr2del")
        );
        assert_eq!(
            placed("ENSP00000000001:p.Gly1del"),
            placed("ENSP00000000001:p.Thr1del")
        );
        assert_eq!(
            placed("ENSP00000000002:p.Val1del"),
            (548, 550, Strand::Reverse, "GTA/-".to_string())
        );
        assert_eq!(
            refused("ENSP00000000001:p.Tyr14del"),
            unable(
                "ENSP00000000001:p.Tyr14del",
                "Reference allele extracted from ENSP00000000001:160-160 (T) does not match \
                 reference allele given by HGVS notation ENSP00000000001:p.Tyr14del (TAC)"
            )
        );
        assert_eq!(
            refused("ENSP00000000001:p.Tyr2_Thr1del"),
            unable(
                "ENSP00000000001:p.Tyr2_Thr1del",
                "Unable to find the genomic reference sequence for protein ENSP00000000001"
            )
        );
    }

    /// A description of a frameshift is refused before any lookup (`:1903-1905`),
    /// whatever its reference; one written as residue, position and `fs` alone
    /// escapes that test, as it does in Ensembl VEP, and fails on `fs` being
    /// no residue.
    #[test]
    fn frameshift_refused_as_perl() {
        for notation in [
            "ENSP00000000001:p.Thr1Alafs*12",
            "ENSP00000000001:p.T1Afs",
            "ENSP00000000001:p.Thr1AlafsTer12",
            "ENSP00000000099:p.Thr1Alafs",
        ] {
            assert_eq!(
                refused(notation),
                unable(
                    notation,
                    "Frameshifts are not supported for HGVS protein input"
                ),
                "{notation}"
            );
        }
        assert_eq!(
            refused("ENSP00000000001:p.Thr1fs"),
            unable(
                "ENSP00000000001:p.Thr1fs",
                "Could not determine nucleotide change from peptide change T -> fs"
            )
        );
    }

    /// A reference that is no translation's or transcript's id is read as a
    /// gene name, and every transcript of the gene is tried in order, the
    /// first to yield a variant taken (`:1933`, `:1946-1985`): `GENE1:p.Thr1Ala`
    /// resolves on the forward transcript, `Val1Ala` on the reverse one after
    /// the forward refuses it, and a residue no transcript carries is refused
    /// with every transcript's message, the non-coding transcript's position
    /// mapping to no exon through Perl's arithmetic on its undefined coding
    /// start. The gene-identifier warning is written for a resolved gene name.
    /// A reference no transcript, translation or gene carries is refused with
    /// the resolver's message. A gene symbol containing `NM` (`DNMT3A`) is
    /// read as a gene name once the transcript lookup misses it, as the gene
    /// lookup follows the external-name one (`:1929-1933`); one containing
    /// `XM` is a gene name too, since the external-name test (`:1915`) names
    /// `NP`, `XP` and `NM` alone; one containing `ENST` is not (`:1922-1925`),
    /// unless a translation carries it, which the translation lookup finds
    /// (`:1909`): `ENSTGUP00000000001`, a translation id whose species code
    /// begins with `T`.
    #[test]
    fn gene_name_tries_every_transcript_in_order() {
        let gene_named = |symbol: &str| {
            let mut named = forward();
            named.gene_symbol = Some(symbol.into());
            vec![named]
        };
        for symbol in ["DNMT3A", "XMGENE"] {
            let notation = format!("{symbol}:p.Thr1Ala");
            assert_eq!(
                read_with(gene_named(symbol), false, &notation)
                    .remove(0)
                    .unwrap()
                    .start,
                121,
                "{notation}"
            );
        }
        let mut other_species = forward();
        other_species.stable_id = "ENSTGUT00000000001".into();
        other_species.translation.as_mut().unwrap().stable_id = "ENSTGUP00000000001".into();
        assert_eq!(
            read_with(vec![other_species], false, "ENSTGUP00000000001:p.Thr1Ala")
                .remove(0)
                .unwrap()
                .start,
            121
        );
        assert_eq!(
            read_with(gene_named("ENSTGENE"), false, "ENSTGENE:p.Thr1Ala")
                .remove(0)
                .unwrap_err(),
            unable(
                "ENSTGENE:p.Thr1Ala",
                "Could not get a Transcript object for 'ENSTGENE'"
            )
        );
        assert_eq!(
            placed("GENE1:p.Thr1Ala"),
            (121, 121, Strand::Forward, "A/G".to_string())
        );
        assert_eq!(
            placed("GENE1:p.Val1Ala"),
            (549, 549, Strand::Reverse, "T/C".to_string())
        );
        assert_eq!(
            refused("GENE1:p.Gly1Ala"),
            unable(
                "GENE1:p.Gly1Ala",
                "Sequence translated from reference (ACG -> T) does not match input sequence (G); \
                 Sequence translated from reference (GTA -> V) does not match input sequence (G); \
                 Unable to find the genomic reference sequence for protein GENE1"
            )
        );
        let mut parser = HgvsParser::new(
            Cursor::new(b"GENE1:p.Thr1Ala\n".to_vec()),
            &context(vec![forward(), reverse()]),
            false,
        );
        assert!(parser.next_variant().unwrap().is_ok());
        assert_eq!(
            parser.take_warnings(),
            ["Possible invalid use of gene or protein identifier 'GENE1' as HGVS reference; most \
              likely transcript will be selected"]
        );
        assert_eq!(
            refused("ENSP00000000099:p.Thr1Ala"),
            unable(
                "ENSP00000000099:p.Thr1Ala",
                "Could not get a Transcript object for 'ENSP00000000099'"
            )
        );
        assert_eq!(
            refused("ENST00000000099:p.Thr1Ala"),
            unable(
                "ENST00000000099:p.Thr1Ala",
                "Could not get a Transcript object for 'ENST00000000099'"
            )
        );
    }

    /// A RefSeq protein or transcript id resolves against a cache that carries
    /// the RefSeq model, the version falling back as the id lookups fall back:
    /// `NP_000001.1` through the translation id, `NM_000001.4` through the
    /// transcript id (Ensembl VEP's external-name lookup in its RefSeq
    /// database, `:1915-1916`); in a cache without it both are refused, the
    /// transcript id as one no transcript carries.
    #[test]
    fn refseq_ids_resolve_where_the_cache_carries_the_model() {
        let mut refseq = forward();
        refseq.stable_id = "NM_000001".into();
        refseq.version = Some(4);
        refseq.source = "RefSeq".into();
        if let Some(translation) = &mut refseq.translation {
            translation.stable_id = "NP_000001".into();
            translation.version = Some(1);
        }
        let with_model = |data: &str| read_with(vec![refseq.clone()], false, data).remove(0);
        assert_eq!(with_model("NP_000001.1:p.Thr1Ala").unwrap().start, 121);
        assert_eq!(with_model("NM_000001.4:p.Thr1Ala").unwrap().start, 121);
        assert_eq!(
            with_model("NM_000001.5:p.Thr1Ala").unwrap_err(),
            unable(
                "NM_000001.5:p.Thr1Ala",
                "Could not get a Transcript object for 'NM_000001.5'"
            )
        );
        assert_eq!(
            refused("NP_000001.1:p.Thr1Ala"),
            unable(
                "NP_000001.1:p.Thr1Ala",
                "Could not get a Transcript object for 'NP_000001.1'"
            )
        );
    }

    /// The codon table is the transcript's (`:2245-2248`): under the vertebrate
    /// mitochondrial table `ATA` codes methionine, so `Val3Met` on `GTA` is the
    /// one-base change to `ATA`, where the standard table's one methionine codon
    /// `ATG` is three bases away and carried whole.
    #[test]
    fn codon_table_two_changes_the_codon_paths() {
        let mut mitochondrial = forward();
        mitochondrial.vefc.as_mut().unwrap().codon_table = 2;
        let under_table_two = read_with(vec![mitochondrial], false, "ENSP00000000001:p.Val3Met")
            .remove(0)
            .unwrap();
        assert_eq!(
            (
                under_table_two.start,
                under_table_two.end,
                under_table_two.allele_string.as_str()
            ),
            (127, 127, "G/A")
        );
        assert_eq!(
            placed("ENSP00000000001:p.Val3Met"),
            (127, 129, Strand::Forward, "GTA/ATG".to_string())
        );
    }

    /// A transcript whose cache record carries no mapper places a residue
    /// through one built from its exons, with the first exon's phase and the
    /// transcript's coding start (TranscriptMapper.pm:95-124): the same bases
    /// as the cached mapper gives, on either strand.
    #[test]
    fn mapper_built_from_the_exons_when_the_cache_carries_none() {
        let mut uncached = vec![forward(), reverse()];
        for transcript in &mut uncached {
            transcript.vefc = None;
        }
        let read_uncached = |data: &str| {
            let v = read_with(uncached.clone(), false, data).remove(0).unwrap();
            (v.start, v.end, v.strand, v.allele_string)
        };
        assert_eq!(
            read_uncached("ENSP00000000001:p.Thr1Ala"),
            placed("ENSP00000000001:p.Thr1Ala")
        );
        assert_eq!(
            read_uncached("ENSP00000000002:p.Val1del"),
            placed("ENSP00000000002:p.Val1del")
        );
    }

    /// A first exon of phase 2 opens mid-codon, so `pep2genomic` reads residue
    /// `n` from cDNA `3n - 2 - 2` (TranscriptMapper.pm:372-388, the coding
    /// start less the start phase): on a transcript coding from cDNA 1 with
    /// that phase, residue 2 is cDNA 2-4, genomic 102-104 (`CGT`, Arg), where a
    /// phase of 0 would read cDNA 4-6 (`TAC`, Tyr). The cached mapper and one
    /// built from the exons place it alike.
    #[test]
    fn a_positive_start_phase_moves_the_codon_frame() {
        let mut cached = forward();
        cached.cdna_coding_start = Some(1);
        let mapper = cached.vefc.as_mut().unwrap().mapper.as_mut().unwrap();
        mapper.cdna_coding_start = 1;
        mapper.start_phase = 2;
        let mut uncached = cached.clone();
        uncached.vefc = None;
        uncached.exons[0].phase = 2;
        for transcript in [cached, uncached] {
            let v = read_with(vec![transcript], false, "ENSP00000000001:p.Arg2Gly")
                .remove(0)
                .unwrap();
            assert_eq!(
                (v.start, v.end, v.strand, v.allele_string.as_str()),
                (102, 102, Strand::Forward, "C/G")
            );
        }
    }

    /// The lookups a reference makes, in order: a reference naming a transcript
    /// (`ENST`, `NM`) goes to the transcript lookup first, since a translation
    /// rarely carries such an id and a missed translation lookup would read the
    /// cache's alternate tier for nothing, then, on a miss, to the translation
    /// lookup, after which an `ENST` reference stops and an `NM` one goes on to
    /// the gene lookup; a translation id goes to the translation lookup alone;
    /// a gene name, `XM` included, goes to the translation lookup, then the
    /// gene lookup (`:1909-1933`).
    #[test]
    fn transcript_references_skip_the_translation_lookup() {
        struct Recording {
            cache: Cache,
            lookups: Mutex<Vec<&'static str>>,
        }
        impl TranscriptSource for Recording {
            fn by_stable_id(&self, stable_id: &str) -> Result<Option<Transcript>, String> {
                self.lookups.lock().unwrap().push("stable_id");
                self.cache.by_stable_id(stable_id)
            }
            fn by_translation_id(&self, stable_id: &str) -> Result<Option<Transcript>, String> {
                self.lookups.lock().unwrap().push("translation_id");
                self.cache.by_translation_id(stable_id)
            }
            fn gene_transcripts(&self, name: &str) -> Result<Vec<Vec<Transcript>>, String> {
                self.lookups.lock().unwrap().push("gene");
                self.cache.gene_transcripts(name)
            }
        }
        let lookups_of = |hgvs: &str| -> (Vec<&'static str>, bool) {
            let source = Arc::new(Recording {
                cache: Cache(vec![forward()]),
                lookups: Mutex::new(Vec::new()),
            });
            let context = ParserContext {
                fasta: Some(Arc::new(Genome)),
                transcripts: Some(Arc::clone(&source) as Arc<dyn TranscriptSource>),
                ..ParserContext::default()
            };
            let mut parser =
                HgvsParser::new(Cursor::new(hgvs.as_bytes().to_vec()), &context, false);
            let resolved = parser.next_variant().unwrap().is_ok();
            let lookups = source.lookups.lock().unwrap().clone();
            (lookups, resolved)
        };
        assert_eq!(
            lookups_of("ENST00000000001:p.Thr1Ala"),
            (vec!["stable_id"], true)
        );
        assert_eq!(
            lookups_of("ENST00000000099:p.Thr1Ala"),
            (vec!["stable_id", "translation_id"], false)
        );
        assert_eq!(
            lookups_of("NM_000001.4:p.Thr1Ala"),
            (vec!["stable_id", "translation_id", "gene"], false)
        );
        assert_eq!(
            lookups_of("ENSP00000000001:p.Thr1Ala"),
            (vec!["translation_id"], true)
        );
        assert_eq!(
            lookups_of("GENE1:p.Thr1Ala"),
            (vec!["translation_id", "gene"], true)
        );
        assert_eq!(
            lookups_of("XM_000001.1:p.Thr1Ala"),
            (vec!["translation_id", "gene"], false)
        );
    }

    /// A cache a lookup cannot read fails the notation with
    /// [`IoError::CacheUnreadable`] carrying the cache's refusal after the
    /// parser's text, whichever lookup meets it: the transcript lookup of an
    /// `ENST` reference, the translation lookup of an `ENSP` one, the gene
    /// lookup of a gene name.
    #[test]
    fn an_unreadable_cache_fails_the_notation_with_its_refusal() {
        struct Failing(&'static str);
        impl Failing {
            fn refusal<T>(&self, lookup: &str, hit: T) -> Result<T, String> {
                if self.0 == lookup {
                    Err("shard 21/1-1000000.json.gz is cut".into())
                } else {
                    Ok(hit)
                }
            }
        }
        impl TranscriptSource for Failing {
            fn by_stable_id(&self, _: &str) -> Result<Option<Transcript>, String> {
                self.refusal("stable_id", None)
            }
            fn by_translation_id(&self, _: &str) -> Result<Option<Transcript>, String> {
                self.refusal("translation_id", None)
            }
            fn gene_transcripts(&self, _: &str) -> Result<Vec<Vec<Transcript>>, String> {
                self.refusal("gene", vec![Vec::new()])
            }
        }
        for (lookup, hgvs) in [
            ("stable_id", "ENST00000000001:p.Thr1Ala"),
            ("translation_id", "ENSP00000000001:p.Thr1Ala"),
            ("gene", "GENE1:p.Thr1Ala"),
        ] {
            let context = ParserContext {
                fasta: Some(Arc::new(Genome)),
                transcripts: Some(Arc::new(Failing(lookup)) as Arc<dyn TranscriptSource>),
                ..ParserContext::default()
            };
            let mut parser =
                HgvsParser::new(Cursor::new(hgvs.as_bytes().to_vec()), &context, false);
            let Some(Err(IoError::CacheUnreadable(message))) = parser.next_variant() else {
                panic!("{hgvs}: the notation must fail with the cache's refusal");
            };
            assert_eq!(
                message,
                "Cannot read the transcript cache for HGVS input: shard 21/1-1000000.json.gz is cut",
                "{hgvs}"
            );
        }
    }

    /// The refusals of `_parse_hgvs_protein_position`, each with Ensembl VEP's
    /// text: a residue the genome's codon does not translate to (a wrong
    /// residue, a lower-case three-letter code, which `ONECODE` does not
    /// read); a new residue no codon codes (`X`, `?`, selenocysteine); a
    /// description without a written residue and without `delins`, or with the
    /// position `0`, which fail the truth test after the codon is read; a
    /// position past the transcript, which maps to no exon; a description the
    /// regular expressions do not match (a substitution written with `>`, a
    /// deletion range whose first residue has no position), which reaches
    /// `pep2genomic` with no position.
    #[test]
    fn unresolvable_residues_and_positions_are_refused() {
        let cases = [
            (
                "ENSP00000000001:p.Ala1Thr",
                "Sequence translated from reference (ACG -> T) does not match input sequence (A)",
            ),
            (
                "ENSP00000000001:p.thr1Ala",
                "Sequence translated from reference (ACG -> T) does not match input sequence (thr)",
            ),
            (
                "ENSP00000000001:p.Thr1Xaa",
                "Could not determine nucleotide change from peptide change T -> X",
            ),
            (
                "ENSP00000000001:p.Thr1?",
                "Could not determine nucleotide change from peptide change T -> ?",
            ),
            (
                "ENSP00000000001:p.Thr1Sec",
                "Could not determine nucleotide change from peptide change T -> U",
            ),
            (
                "ENSP00000000001:p.1Ala",
                "Could not parse HGVS protein notation ENSP00000000001:p.1Ala",
            ),
            (
                "ENSP00000000001:p.Thr0Ala",
                "Could not parse HGVS protein notation ENSP00000000001:p.Thr0Ala",
            ),
            (
                "ENSP00000000001:p.Thr99Ala",
                "Unable to find the genomic reference sequence for protein ENSP00000000001",
            ),
            (
                "ENSP00000000001:p.Thr99del",
                "Unable to find the genomic reference sequence for protein ENSP00000000001",
            ),
            ("ENSP00000000001:p.Thr1>Ala", "Must call with start and end"),
            ("ENSP00000000001:p.1_Thrdel", "Must call with start and end"),
        ];
        for (notation, message) in cases {
            assert_eq!(refused(notation), unable(notation, message), "{notation}");
        }
    }

    /// `validate_vf`'s size test applies to a protein record as to any other:
    /// a deletion of thirteen codons within the first exon spans 38 bases and
    /// is skipped under `--max_sv_size 30` with Perl's text, kept under 38.
    #[test]
    fn protein_record_over_max_sv_size_is_skipped() {
        let read_bounded = |max_sv_size: Option<u64>| {
            let mut parser = HgvsParser::new(
                Cursor::new(b"ENSP00000000001:p.Thr1_Thr13del\n".to_vec()),
                &context_bounded(vec![forward()], max_sv_size),
                false,
            );
            parser.next_variant().unwrap().map_err(|e| e.to_string())
        };
        assert_eq!(
            read_bounded(Some(30)).unwrap_err(),
            "variant size (38) is bigger than --max_sv_size (30)"
        );
        assert_eq!(read_bounded(Some(38)).unwrap().end, 159);
    }

    /// The captures of the resolver's regular expressions, as Perl's engine
    /// takes them (`:2236`, `:2334`, `:2337`, `:2342`): the digits backtrack to
    /// leave a new residue, so `97` alone reads as position 9 and residue `7`
    /// and `Gln1111` as position 111 and residue `1`; a range's underscore
    /// text is a word, so `Arg97_Gly98delinsAla` reads as one residue `_Gly98delinsAla`;
    /// `delins` is taken when present; a deletion's written residue may be
    /// digits.
    #[test]
    fn captures_follow_perls_regular_expressions() {
        let sub = |d: &str| {
            substitution_captures(d)
                .map(|(from, pos, to)| (from.map(str::to_string), pos.to_string(), to.to_string()))
        };
        let owned = |from: Option<&str>, pos: &str, to: &str| {
            Some((from.map(str::to_string), pos.to_string(), to.to_string()))
        };
        assert_eq!(sub("Arg97Gly"), owned(Some("Arg"), "97", "Gly"));
        assert_eq!(sub("R97G"), owned(Some("R"), "97", "G"));
        assert_eq!(sub("97delinsAla"), owned(None, "97", "Ala"));
        assert_eq!(sub("Met97delinsAla"), owned(Some("Met"), "97", "Ala"));
        assert_eq!(sub("Arg97*"), owned(Some("Arg"), "97", "*"));
        assert_eq!(sub("Arg97="), owned(Some("Arg"), "97", "="));
        assert_eq!(sub("Arg97?"), owned(Some("Arg"), "97", "?"));
        assert_eq!(sub("97"), owned(None, "9", "7"));
        assert_eq!(
            sub("Arg97_Gly98delinsAla"),
            owned(Some("Arg"), "97", "_Gly98delinsAla")
        );
        assert_eq!(sub("Arg0Gly"), owned(Some("Arg"), "0", "Gly"));
        assert_eq!(sub("97delins"), owned(None, "97", "delins"));
        assert_eq!(sub("Arg97>Gly"), None);
        assert_eq!(sub("Arg"), None);
        assert_eq!(
            deletion_captures("Lys1110del"),
            Some(("Lys", "1110", "del"))
        );
        assert_eq!(deletion_captures("1110del"), Some(("1", "110", "del")));
        assert_eq!(deletion_captures("Gln1111"), Some(("Gln", "111", "1")));
        assert_eq!(deletion_captures("Lys1110*"), Some(("Lys", "1110", "*")));
        assert_eq!(deletion_captures("del"), None);
        assert_eq!(deletion_first_captures("Lys1110"), Some(("Lys", "1110")));
        assert_eq!(deletion_first_captures("11110"), Some(("1", "1110")));
        assert_eq!(deletion_first_captures("Lys"), None);
        for (description, frameshift) in [
            ("Arg97Glyfs*12", true),
            ("R97Gfs*12", true),
            ("Arg97Glyfs", true),
            ("Arg97GlyfsTer12", true),
            ("Arg97fs", false),
            ("Met1fs", false),
            ("97fs", false),
            ("Argfs12", false),
            ("Arg97Ter", false),
        ] {
            assert_eq!(is_frameshift(description), frameshift, "{description}");
        }
    }

    /// BioPerl's codon table as the resolver uses it: `revtranslate` lists an
    /// amino acid's codons in the table's order under the standard table and
    /// the vertebrate mitochondrial one, reads a three-letter code in any case
    /// and a one-letter code in upper case alone, and lists nothing for a
    /// residue no codon codes; `translate` reads an ambiguous base through
    /// every nucleotide it stands for.
    #[test]
    fn codon_table_reads_as_bioperl() {
        let spelled = |codons: Vec<[u8; 3]>| -> Vec<String> {
            codons
                .into_iter()
                .map(|c| String::from_utf8(c.to_vec()).unwrap())
                .collect()
        };
        assert_eq!(
            spelled(reverse_translate("R", 1)),
            ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
        );
        assert_eq!(
            spelled(reverse_translate("arg", 1)),
            spelled(reverse_translate("R", 1))
        );
        assert_eq!(spelled(reverse_translate("Ter", 1)), ["TAA", "TAG", "TGA"]);
        assert_eq!(
            spelled(reverse_translate("*", 2)),
            ["TAA", "TAG", "AGA", "AGG"]
        );
        assert_eq!(spelled(reverse_translate("Trp", 2)), ["TGA", "TGG"]);
        assert_eq!(
            spelled(reverse_translate("B", 1)),
            ["GAT", "GAC", "AAT", "AAC"]
        );
        for none in ["r", "X", "U", "fs", "?", "=", "delins", "Xaa"] {
            assert!(reverse_translate(none, 1).is_empty(), "{none}");
        }
        assert_eq!(translate(b"CGC", 1), "R");
        assert_eq!(translate(b"cgc", 1), "R");
        assert_eq!(translate(b"ATA", 2), "M");
        assert_eq!(translate(b"TGA", 2), "W");
        assert_eq!(translate(b"CGN", 1), "R");
        assert_eq!(translate(b"RAY", 1), "B");
        assert_eq!(translate(b"SAR", 1), "Z");
        assert_eq!(translate(b"NNN", 1), "X");
        assert_eq!(translate(b"CG", 1), "");
    }
}
