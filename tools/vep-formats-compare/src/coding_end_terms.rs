// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The coding terms vep-rs owes on a row at a coding end of a transcript, derived from the
//! row and from the reading of the transcript's cache record (`CodingEnds`).
//!
//! Ensembl VEP writes thirteen terms from the predicates of its coding evaluation
//! (`Utils/VariationEffect.pm`), and vep-rs evaluates the same predicates
//! (`crates/vep-effects/src/coding.rs`, `PerlCodingEval`). The engines part in the tests
//! that ask for a stop codon or a start codon: Ensembl VEP reads the last and the first
//! three coding positions as those codons on every transcript without `cds_end_NF` or
//! `cds_start_NF` (VariationEffect.pm:1373-1381, 975-986), vep-rs only where the codon
//! stands, and on a stop codon that polyadenylation completes vep-rs reads the row on the
//! completed coding sequence where Ensembl VEP reads an incomplete codon
//! (VariationEffect.pm:1588-1614). `expected_terms` evaluates the thirteen predicates
//! with those tests answered as vep-rs answers them, on the cells both rows print equal
//! or, on a completed stop codon, on the cached coding sequence followed by adenines. An
//! `Amino_acids` and a `Codons` cell may differ as the class
//! `alternate_codon_past_the_coding_end` documents, vep-rs's ending where the coding
//! sequence ends: Ensembl VEP's are read there. A row it cannot derive is "not derived",
//! with the reason.
//!
//! The file names no type of a comparison tool: a row comes in as plain values and the
//! terms go out as plain values. It asks the two shapes of that class of the module
//! `coding_end`, the name under which every tool root that compiles this file compiles
//! `rules/coding_end.rs`.

use super::{CodingEnds, CodingModel, StartEnd, StopEnd, UnshownStartSequence};
use std::collections::BTreeSet;
use std::fmt;

/// A term of the coding family: the terms Ensembl VEP's coding evaluation writes, which
/// vep-rs writes and replaces as one set (`crates/vep-effects/src/coding.rs`,
/// `is_perl_coding_term`). The order is the order vep-rs evaluates them in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CodingTerm {
    StopGained,
    FrameshiftVariant,
    StopLost,
    StartLost,
    InframeInsertion,
    InframeDeletion,
    MissenseVariant,
    ProteinAlteringVariant,
    IncompleteTerminalCodonVariant,
    StartRetainedVariant,
    StopRetainedVariant,
    SynonymousVariant,
    CodingSequenceVariant,
}

impl CodingTerm {
    /// Every term of the family.
    pub const ALL: [CodingTerm; 13] = [
        CodingTerm::StopGained,
        CodingTerm::FrameshiftVariant,
        CodingTerm::StopLost,
        CodingTerm::StartLost,
        CodingTerm::InframeInsertion,
        CodingTerm::InframeDeletion,
        CodingTerm::MissenseVariant,
        CodingTerm::ProteinAlteringVariant,
        CodingTerm::IncompleteTerminalCodonVariant,
        CodingTerm::StartRetainedVariant,
        CodingTerm::StopRetainedVariant,
        CodingTerm::SynonymousVariant,
        CodingTerm::CodingSequenceVariant,
    ];

    /// The term as a `Consequence` cell writes it.
    pub fn name(self) -> &'static str {
        match self {
            CodingTerm::StopGained => "stop_gained",
            CodingTerm::FrameshiftVariant => "frameshift_variant",
            CodingTerm::StopLost => "stop_lost",
            CodingTerm::StartLost => "start_lost",
            CodingTerm::InframeInsertion => "inframe_insertion",
            CodingTerm::InframeDeletion => "inframe_deletion",
            CodingTerm::MissenseVariant => "missense_variant",
            CodingTerm::ProteinAlteringVariant => "protein_altering_variant",
            CodingTerm::IncompleteTerminalCodonVariant => "incomplete_terminal_codon_variant",
            CodingTerm::StartRetainedVariant => "start_retained_variant",
            CodingTerm::StopRetainedVariant => "stop_retained_variant",
            CodingTerm::SynonymousVariant => "synonymous_variant",
            CodingTerm::CodingSequenceVariant => "coding_sequence_variant",
        }
    }

    /// The term a `Consequence` cell writes as `name`, `None` for a term outside the
    /// family.
    pub fn of_name(name: &str) -> Option<CodingTerm> {
        CodingTerm::ALL.into_iter().find(|term| term.name() == name)
    }
}

impl fmt::Display for CodingTerm {
    fn fmt(&self, text: &mut fmt::Formatter<'_>) -> fmt::Result {
        text.write_str(self.name())
    }
}

/// The terms of the coding family in a `Consequence` cell, its terms joined by commas.
pub fn coding_family(consequence: &str) -> BTreeSet<CodingTerm> {
    consequence
        .split(',')
        .filter_map(|term| CodingTerm::of_name(term.trim()))
        .collect()
}

/// The terms of a `Consequence` cell that lie outside the coding family.
pub fn outside_the_coding_family(consequence: &str) -> BTreeSet<&str> {
    consequence
        .split(',')
        .map(str::trim)
        .filter(|term| !term.is_empty() && CodingTerm::of_name(term).is_none())
        .collect()
}

/// The five cells of one row that the predicates read, each as the default format
/// writes it; a cell the row leaves blank is `-` or the empty text. A position cell may
/// end in `/` and a length (`--total_length`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cells<'a> {
    pub cdna_position: &'a str,
    pub cds_position: &'a str,
    pub protein_position: &'a str,
    pub amino_acids: &'a str,
    pub codons: &'a str,
}

/// A position cell as the default format writes it, of the two members a JSON entry
/// carries it as (`cds_start` and `cds_end`, and their kin): one number where both are
/// equal, the two joined by a dash, `?` for a member the entry lacks beside one it
/// carries, and `-` where it carries neither (OutputFactory/JSON.pm:332-341 undoes
/// `format_coords` so).
pub fn position_cell(start: Option<&str>, end: Option<&str>) -> String {
    match (start, end) {
        (None, None) => "-".to_string(),
        (Some(start), Some(end)) if start == end => start.to_string(),
        (start, end) => format!("{}-{}", start.unwrap_or("?"), end.unwrap_or("?")),
    }
}

/// One pair's row: what both engines' rows agree on by their key, and Ensembl VEP's
/// terms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RowFacts<'a> {
    /// Whether the row's `Feature_type` is `Transcript`.
    pub transcript_row: bool,
    /// The row's `Allele` cell: bases, `-` for none, or the class word of a structural
    /// allele.
    pub allele: &'a str,
    /// The lowest and the highest genomic base of the row's Location; an insertion's are
    /// the two bases it stands between.
    pub location: (i64, i64),
    /// Ensembl VEP's `Consequence` cell, its terms joined by commas.
    pub ensembl_vep_terms: &'a str,
}

/// What a comparison read beside the pair's key and terms.
#[derive(Clone, Copy)]
pub struct WhatWasRead<'a> {
    /// The five cells of Ensembl VEP's row and of vep-rs's, in that order; `None` where
    /// the output carries no such columns (`--fields`, `--summary`, `--most_severe`).
    pub cells: Option<(Cells<'a>, Cells<'a>)>,
    /// The count of reference bases the allele stands in for at the row's Location, 0
    /// for an insertion; `None` where the input record was not read.
    pub reference_length: Option<u64>,
    /// The reference base at a genomic position of the row's sequence, `None` past its
    /// ends; `None` where the comparison has no reference sequence. One shape of row asks
    /// it: a deletion that begins in a completed stop codon's transcript and runs past
    /// the transcript's last base.
    pub reference_bases: Option<&'a dyn Fn(i64) -> Option<u8>>,
    /// Whether the input record is given on the reverse strand (the strand column of
    /// the default input format). vep-rs reads an allele reversed and complemented where
    /// the record's strand is not the transcript's; a row of such a record is not
    /// derived where the allele's bases are read, on a completed stop codon.
    pub record_on_the_reverse_strand: bool,
}

/// Why a row's terms are not derived.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NotDerived {
    /// The row is no Transcript row.
    NotATranscriptRow,
    /// No cache record answers for the row's transcript, or the record is not read.
    NoVerdict,
    /// The output carries none of the five cell columns.
    NoCells,
    /// The input record was not read, so the allele's reference length is not known.
    NoReferenceLength,
    /// The allele is neither bases of the four alone, nor `-`, nor the class word of a
    /// structural allele: it holds a letter that names no base.
    AlleleNamesNoBases,
    /// The named cell differs between the two rows, otherwise than `equal_cells` admits.
    CellsDiffer(&'static str),
    /// The named cell is in no form `format_coords` or the peptide and codon writers
    /// give.
    CellNotRead(&'static str),
    /// Ensembl VEP writes both start terms where a start codon stands: vep-rs keeps one
    /// of them there, by another rule.
    BothStartTermsWhereAStartCodonStands,
    /// A whole stop codon stands, the peptides do not decide the stop terms, and the
    /// allele does not lie clear of the stop codon: vep-rs reads the edited sequence
    /// there, which the cells do not show.
    AtAStandingStopCodon,
    /// A start codon stands, the allele meets it, and a term the start terms depend on
    /// is derived otherwise than Ensembl VEP writes it.
    AtAStandingStartCodon,
    /// The input record is given on the reverse strand, on a completed stop codon.
    RecordOnTheReverseStrand,
    /// The record's coordinates are not kept for the transcript, and the row needs them.
    ModelNotKept,
    /// An insertion on a model that holds bases the genome lacks, whose coordinates are
    /// not kept: its place cannot be told from the cells.
    InsertionOnAnEditedModel,
    /// The record's sequence is not kept for the transcript, and the row needs it.
    SequenceNotKept,
    /// The codons the allele spans hold a residue an edit of the sequence writes, or
    /// protein position 1 where the first codon is not read as any other.
    SequenceEdit,
    /// A shape of allele on a completed stop codon that the derivation does not cover.
    ShapeNotCovered,
    /// A span that lacks a cDNA or CDS coordinate on a completed stop codon, other than a
    /// deletion that runs past the transcript's last base.
    SpanLacksACoordinate,
    /// The deletion past the transcript's last base can be written with its end inside
    /// the transcript, where the bases after it are read.
    PlacementInsideTheTranscript,
    /// The row holds no cell of the kind asked, or a shape of notation that is not
    /// derived.
    CellNotDerived,
}

impl fmt::Display for NotDerived {
    fn fmt(&self, text: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NotDerived::NotATranscriptRow => text.write_str("no Transcript row"),
            NotDerived::NoVerdict => text.write_str("no coding ends read for the transcript"),
            NotDerived::NoCells => {
                text.write_str("the output carries no position, residue or codon cells")
            }
            NotDerived::NoReferenceLength => text.write_str("the input record was not read"),
            NotDerived::AlleleNamesNoBases => {
                text.write_str("the allele is neither bases nor a structural class word")
            }
            NotDerived::CellsDiffer(cell) => write!(text, "{cell} differs between the rows"),
            NotDerived::CellNotRead(cell) => write!(text, "{cell} is in a form that is not read"),
            NotDerived::BothStartTermsWhereAStartCodonStands => {
                text.write_str("both start terms where a start codon stands")
            }
            NotDerived::AtAStandingStopCodon => {
                text.write_str("the allele is not clear of a stop codon that stands")
            }
            NotDerived::AtAStandingStartCodon => {
                text.write_str("the allele meets a start codon that stands")
            }
            NotDerived::RecordOnTheReverseStrand => {
                text.write_str("the record is given on the reverse strand")
            }
            NotDerived::ModelNotKept => text.write_str("the transcript's coordinates are not kept"),
            NotDerived::InsertionOnAnEditedModel => {
                text.write_str("an insertion on a model that holds bases the genome lacks")
            }
            NotDerived::SequenceNotKept => text.write_str("the transcript's sequence is not kept"),
            NotDerived::SequenceEdit => text.write_str("the codons hold an edited residue"),
            NotDerived::ShapeNotCovered => text.write_str("a shape that is not covered"),
            NotDerived::SpanLacksACoordinate => text.write_str("the span lacks a coordinate"),
            NotDerived::PlacementInsideTheTranscript => {
                text.write_str("the deletion can be written with its end inside the transcript")
            }
            NotDerived::CellNotDerived => text.write_str("no such cell is derived for the row"),
        }
    }
}

/// What the terms read at the stop end of a transcript.
#[derive(Clone, Copy, PartialEq, Eq)]
enum StopCodon {
    None { flagged: bool },
    Whole,
    Completed,
}

fn stop_codon(stop: StopEnd) -> Option<StopCodon> {
    match stop {
        StopEnd::Unread => None,
        StopEnd::NoneShown { flagged } => Some(StopCodon::None { flagged }),
        StopEnd::BegunByTheTail { .. } => Some(StopCodon::Completed),
        StopEnd::Whole => Some(StopCodon::Whole),
    }
}

/// A row as it is derived: on Ensembl VEP's terms alone (a structural allele), on the
/// cells, or on the completed coding sequence.
enum Derivation {
    Structural(BTreeSet<CodingTerm>),
    Cells(BTreeSet<CodingTerm>),
    Completed(OnTheCompletedSequence),
}

/// The derivation of a row, step by step:
///
/// 1. a Transcript row, with both coding ends of its transcript read;
/// 2. a structural allele is read on Ensembl VEP's terms (`structural_terms`);
/// 3. a sequence allele needs the five cells of both rows and the length of the
///    record's reference bases;
/// 4. on a transcript whose stop codon polyadenylation completes it is read on the
///    completed coding sequence (`on_the_completed_sequence`), elsewhere on the cells
///    (`on_the_cells`).
fn derive(
    row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Result<Derivation, NotDerived> {
    if !row.transcript_row {
        return Err(NotDerived::NotATranscriptRow);
    }
    let ends = ends.ok_or(NotDerived::NoVerdict)?;
    let (Some(stop), false) = (stop_codon(ends.stop), ends.start == StartEnd::Unread) else {
        return Err(NotDerived::NoVerdict);
    };
    let allele = match allele_of(row.allele)? {
        Allele::Structural => return structural_terms(row, ends, stop).map(Derivation::Structural),
        Allele::Bases(allele) => allele,
    };
    let cells = read.cells.as_ref().ok_or(NotDerived::NoCells)?;
    let reference = read.reference_length.ok_or(NotDerived::NoReferenceLength)?;
    let lengths = Lengths {
        reference: i64::try_from(reference).unwrap_or(i64::MAX),
        allele: allele.len() as i64,
    };
    if stop != StopCodon::Completed {
        return on_the_cells(row, ends, stop, cells, lengths).map(Derivation::Cells);
    }
    if read.record_on_the_reverse_strand {
        return Err(NotDerived::RecordOnTheReverseStrand);
    }
    on_the_completed_sequence(row, ends, cells, lengths, allele, read.reference_bases)
        .map(Derivation::Completed)
}

/// The coding terms vep-rs owes on the row: the terms of the coding family a correct
/// vep-rs writes on it, the other terms being Ensembl VEP's own. `ends` is the reading
/// of the cache record of the row's transcript, `None` where no record answers.
pub fn expected_terms(
    row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Result<BTreeSet<CodingTerm>, NotDerived> {
    Ok(match derive(row, ends, read)? {
        Derivation::Structural(terms) | Derivation::Cells(terms) => terms,
        Derivation::Completed(completed) => completed.terms,
    })
}

fn completed(
    row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Result<OnTheCompletedSequence, NotDerived> {
    match derive(row, ends, read)? {
        Derivation::Completed(completed) => Ok(completed),
        Derivation::Structural(_) | Derivation::Cells(_) => Err(NotDerived::CellNotDerived),
    }
}

/// The `Amino_acids` cell vep-rs owes on a row of a transcript whose stop codon
/// polyadenylation completes: the reference peptide of the completed coding sequence,
/// and the alternate after a `/` where it differs. Ensembl VEP's cell there is the
/// incomplete codon's.
pub fn amino_acids_on_a_completed_stop_codon(
    row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Result<String, NotDerived> {
    completed(row, ends, read).map(|completed| completed.amino_acids)
}

/// The `Codons` cell vep-rs owes on a row of a transcript whose stop codon
/// polyadenylation completes: the two codons of the completed coding sequence in lower
/// case, each allele's own bases in upper case.
pub fn codons_on_a_completed_stop_codon(
    row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Result<String, NotDerived> {
    completed(row, ends, read).map(|completed| completed.codons)
}

/// The protein notation vep-rs owes on a row of a transcript whose stop codon
/// polyadenylation completes, from `p.` on, without the protein's id before it. Three
/// shapes are derived, each on the completed message, whose poly(A) tail holds no
/// further stop codon:
///
/// - the stop codon is the first residue that changes, on a row with `stop_lost` and
///   no frameshift or with a frameshift and no `stop_lost`: the extension with no new
///   stop, `p.Ter<n><Res>extTer?`, `n` the stop's residue number and `Res` the residue
///   of the three bases the edited sequence holds at the old stop position;
/// - a frameshift whose first changed residue lies before the stop codon:
///   `p.<Ref><n><Alt>fsTer?`;
/// - a substitution of as many bases that changes the last residue and the stop
///   codon: `p.<Ref><n>_Ter<n+1>delins<Alt><Alt>`.
///
/// Any other row is not derived.
pub fn protein_notation_on_a_completed_stop_codon(
    row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Result<String, NotDerived> {
    completed(row, ends, read)?
        .notation
        .ok_or(NotDerived::CellNotDerived)
}

/// What `protein_notation` reads of a row on the completed coding sequence.
struct ProteinChange<'a> {
    /// The transcript's peptide.
    peptide: &'a [u8],
    whole_codons: usize,
    /// The completed coding sequence with the allele written in and adenines after it.
    alternate_sequence: Option<&'a [u8]>,
    /// The first codon the span names.
    first_codon: Option<i64>,
    substitution: bool,
    frameshift: bool,
    stop_lost: bool,
}

fn protein_notation(change: &ProteinChange<'_>, code: &super::GeneticCode) -> Option<String> {
    let alternate = change.alternate_sequence?;
    let first = usize::try_from(change.first_codon?)
        .ok()
        .filter(|first| *first >= 1)?;
    let stop_place = change.whole_codons + 1;
    if change.peptide.len() != change.whole_codons || first > stop_place {
        return None;
    }
    let mut reference = change.peptide[first - 1..].to_vec();
    reference.push(b'*');
    let read: Vec<u8> = alternate
        .get(3 * (first - 1)..)?
        .chunks_exact(3)
        .map(|codon| code.residue(codon))
        .collect();
    // A stop codon in the sequence read gives a count of residues to it, a shape that
    // is not derived.
    if read.contains(&b'*') {
        return None;
    }
    let changed = reference
        .iter()
        .zip(&read)
        .position(|(reference, read)| reference != read)?;
    let place = first + changed;
    let name = |residue: u8| three_letters(residue);
    if place == stop_place {
        return (change.stop_lost != change.frameshift)
            .then(|| name(read[changed]))
            .flatten()
            .map(|residue| format!("p.Ter{place}{residue}extTer?"));
    }
    if change.frameshift && !change.stop_lost {
        return Some(format!(
            "p.{}{place}{}fsTer?",
            name(reference[changed])?,
            name(read[changed])?
        ));
    }
    if change.substitution && change.stop_lost && place + 1 == stop_place {
        return Some(format!(
            "p.{}{place}_Ter{stop_place}delins{}{}",
            name(reference[changed])?,
            name(read[changed])?,
            name(*read.get(changed + 1)?)?
        ));
    }
    None
}

/// What a frameshift deletes where the transcript's sequence ends before its first
/// changed codon is whole, as the frameshift routine of both engines words it
/// (`crates/vep-effects/src/coding.rs`, `hgvsp_fs_peptides`, on
/// TranscriptVariationAllele.pm's `_get_fs_peptides`, :2342-2346).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CutShortFrameshift {
    /// The deletion after `p.`: `<Ref><n>del`, or `<First><n>_<Last><m>del` where more
    /// than one residue is named.
    deletion: String,
    /// The first of the deleted residues with its number, `<First><n>`, as the deletion
    /// opens.
    first_residue: String,
}

impl CutShortFrameshift {
    /// The protein notation vep-rs owes on the row, after `p.`, and `None` for a row on
    /// which it owes none: the short frameshift form of the HGVS nomenclature,
    /// `<First><n>fs`, the first of the deleted residues and nothing after it
    /// (`crates/vep-effects/src/coding.rs`, the `Del` arm of `hgvsp_format`). The
    /// alternate sequence holds no whole codon from the first changed one on to name a
    /// residue or a stop by, and no stop codon stands at the coding end for the
    /// frameshift to lose.
    pub fn vep_rs_notation(&self) -> Option<String> {
        Some(format!("{}fs", self.first_residue))
    }

    /// The protein notation Ensembl VEP writes on the row, after `p.`: the deletion
    /// followed by an extension of unknown length, because its lost-stop test holds on
    /// the row and its alternate translation meets no stop codon
    /// (TranscriptVariationAllele.pm:1916-1936 with `_stop_loss_extra_AA`).
    pub fn ensembl_vep_notation(&self) -> String {
        format!("{}extTer?", self.deletion)
    }
}

/// The allele's bases as the transcript reads them: reversed and complemented on a
/// transcript of the reverse strand.
fn as_the_transcript_reads(allele: &[u8], strand: Option<i64>) -> Vec<u8> {
    if strand != Some(-1) {
        return allele.to_vec();
    }
    allele
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// The peptide Ensembl VEP writes for the codons `codons` of an allele
/// (TranscriptVariationAllele.pm's `peptide`): the residues of the whole codons under
/// `code` and the unknown residue for bases that fill no codon after them, which a
/// peptide that is one stop does not take.
fn peptide_of_codons(codons: &[u8], code: &super::GeneticCode) -> Vec<u8> {
    let whole = codons.len() - codons.len() % 3;
    let mut peptide: Vec<u8> = codons[..whole]
        .chunks_exact(3)
        .map(|codon| code.residue(codon))
        .collect();
    if whole < codons.len() && peptide != b"*" {
        peptide.push(b'X');
    }
    peptide
}

/// The frameshift of a row at a coding end where no stop codon stands and the
/// transcript carries no `cds_end_NF`, on which Ensembl VEP writes `frameshift_variant`
/// beside `stop_lost` and vep-rs `frameshift_variant` alone, where the transcript's
/// sequence ends before the codon of the row's translation start is whole. The
/// alternate sequence is the cached coding sequence with the allele written over the
/// CDS span and the 3' UTR after it; where it holds fewer whole codons than the
/// translation start counts, the frameshift routine of both engines names the row's
/// reference peptide as deleted and reads no further
/// (TranscriptVariationAllele.pm:2342-2346).
///
/// The reference peptide is read from the cached coding sequence: the codons from the
/// translation start to the translation end, a codon the coding sequence cuts short as
/// the unknown residue. The alternate peptide of such a row is at most the unknown
/// residue of a codon cut short, so the two share no residue but that one at their
/// end, which is clipped with the end of the span (`_clip_alleles`, :2228-2242). One
/// residue is named with the translation start, more than one as the first with the
/// translation start and the last with the end. Ensembl VEP writes that deletion with an
/// extension of unknown length after it, and vep-rs the first residue with the
/// translation start as a frameshift.
///
/// Not derived: a row whose terms are not those, on either side or as derived; a
/// transcript whose stop codon stands, that carries the flag, or whose sequence is not
/// kept; a record given on the reverse strand; a span without both CDS ends and both
/// translation ends; an alternate sequence that holds the codon of the translation
/// start whole; a residue of a whole codon that the cached peptide names otherwise; a
/// reference or an alternate peptide that is not the one Ensembl VEP's `Amino_acids`
/// cell writes; a reference peptide of which no residue is left; and a residue outside
/// the twenty.
pub fn frameshift_the_coding_end_cuts_short(
    row: &RowFacts<'_>,
    vep_rs_terms: &str,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Result<CutShortFrameshift, NotDerived> {
    let refused = NotDerived::CellNotDerived;
    let Derivation::Cells(owed) = derive(row, ends, read)? else {
        return Err(refused);
    };
    let ends = ends.ok_or(NotDerived::NoVerdict)?;
    let alone = BTreeSet::from([CodingTerm::FrameshiftVariant]);
    let beside_a_lost_stop = BTreeSet::from([CodingTerm::FrameshiftVariant, CodingTerm::StopLost]);
    if ends.stop != (StopEnd::NoneShown { flagged: false })
        || owed != alone
        || coding_family(vep_rs_terms) != alone
        || coding_family(row.ensembl_vep_terms) != beside_a_lost_stop
    {
        return Err(refused);
    }
    if read.record_on_the_reverse_strand {
        return Err(NotDerived::RecordOnTheReverseStrand);
    }
    let kept = ends
        .unclosed_end
        .as_deref()
        .ok_or(NotDerived::SequenceNotKept)?;
    let (Some((theirs, _)), Some(reference_length), Ok(Allele::Bases(allele))) = (
        read.cells.as_ref(),
        read.reference_length,
        allele_of(row.allele),
    ) else {
        return Err(refused);
    };
    let code = super::CODES.get(&kept.table).ok_or(refused)?;
    let lengths = Lengths {
        reference: i64::try_from(reference_length).unwrap_or(i64::MAX),
        allele: allele.len() as i64,
    };
    let neighbours = between_neighbours(ends, row.location, lengths)?;
    let (Some(cds_start), Some(cds_end)) =
        positions(theirs.cds_position, "CDS_position", neighbours)?
    else {
        return Err(refused);
    };
    let (Some(first), Some(last)) =
        positions(theirs.protein_position, "Protein_position", neighbours)?
    else {
        return Err(refused);
    };
    if cds_start < 1 || cds_end < cds_start - 1 || first < 1 || last < first {
        return Err(refused);
    }
    let allele = as_the_transcript_reads(allele, ends.strand);
    let place = |position: i64| {
        usize::try_from(position)
            .unwrap_or(0)
            .min(kept.coding.len())
    };
    let mut with_the_allele = kept.coding[..place(cds_start - 1)].to_vec();
    with_the_allele.extend_from_slice(&allele);
    with_the_allele.extend_from_slice(&kept.coding[place(cds_end)..]);
    // The sequence the frameshift routine translates: fewer than three bases of coding
    // sequence are dropped, and the 3' UTR follows.
    let translated = match with_the_allele.len() < 3 {
        true => kept.three_prime_utr.len(),
        false => with_the_allele.len() + kept.three_prime_utr.len(),
    };
    if first <= (translated / 3) as i64 {
        return Err(refused);
    }
    // The codons of the two alleles (TranscriptVariationAllele.pm's `codon`): the
    // reference's from the first base of the translation start's codon to the last of
    // the translation end's, the alternate's as many bases less the bases the allele
    // takes away.
    let before_the_window = (i128::from(first) - 1) * 3;
    let window = place(nearest(before_the_window));
    let codons = (i128::from(last) - i128::from(first) + 1) * 3;
    let taken_away = cds_end - cds_start + 1 - allele.len() as i64;
    let window_end = place(nearest(i128::from(last) * 3));
    let reference = peptide_of_codons(&kept.coding[window..window_end], code);
    let alternate_end = place(nearest(
        before_the_window + (codons - i128::from(taken_away)).max(0),
    ));
    let alternate = peptide_of_codons(
        with_the_allele
            .get(window..alternate_end.min(with_the_allele.len()))
            .unwrap_or_default(),
        code,
    );
    // A residue an edit of the sequence writes is the cached peptide's, not its codon's.
    let edited = reference.iter().enumerate().any(|(offset, residue)| {
        *residue != b'X'
            && kept
                .peptide
                .get(place(first - 1) + offset)
                .is_some_and(|cached| cached != residue)
    });
    if edited {
        return Err(NotDerived::SequenceEdit);
    }
    if peptide_alleles(theirs.amino_acids)? != Some((reference.clone(), alternate.clone())) {
        return Err(refused);
    }
    let (mut named, mut end) = (reference, last);
    if alternate == b"X" && named.last() == Some(&b'X') {
        named.pop();
        end -= 1;
    }
    let names: Vec<&str> = named
        .iter()
        .map(|residue| three_letters(*residue).filter(|name| *name != "Ter"))
        .collect::<Option<_>>()
        .ok_or(refused)?;
    let (opening, deletion) = match names.as_slice() {
        [] => return Err(refused),
        [one] => (one, format!("{one}{first}del")),
        [opening, .., closing] => (opening, format!("{opening}{first}_{closing}{end}del")),
    };
    Ok(CutShortFrameshift {
        first_residue: format!("{opening}{first}"),
        deletion,
    })
}

/// The protein notation vep-rs owes on a row at the start of a transcript that shows no
/// start codon, from `p.` on, without the protein's id before it. Three kinds of row are
/// derived:
///
/// - a row whose derived terms hold `frameshift_variant`, on which Ensembl VEP writes a
///   start term: the frameshift of its printed span
///   (`frameshift_where_no_start_codon_is_shown`);
/// - a row whose derived terms hold no `frameshift_variant` because its `CDS_position`
///   is open at its start: a deletion, an allele without bases, that begins before the
///   coding start, whose two rows carry one `HGVS_OFFSET` above 0, whose span moved 3'
///   by that offset lies inside the coding sequence, and whose reference bases, as many
///   as the moved span has positions, are no multiple of three. Both engines make the
///   notation on the moved deletion, so the row is read as the frameshift of the moved
///   span, by the same function (`span_the_notation_reads`);
/// - an allele of as many bases as it stands in for whose two peptides, in the
///   `Amino_acids` cell both rows print equal, are one residue each at protein position
///   1: `p.<Ref>1<Alt>`, with `Ter` for a stop, or `p.<Ref>1=` where the two are equal
///   (the default format writes the `=` as `%3D`).
///
/// Ensembl VEP's is the lost start at residue 1, `p.<Ref>1?`
/// (TranscriptVariationAllele.pm:2161). Any other row is not derived. `offsets` is the
/// `HGVS_OFFSET` cell of Ensembl VEP's row and of vep-rs's, each as the row writes it,
/// blank for a row that carries none.
pub fn protein_notation_where_no_start_codon_is_shown(
    row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
    offsets: (&str, &str),
) -> Result<String, NotDerived> {
    let Derivation::Cells(terms) = derive(row, ends, read)? else {
        return Err(NotDerived::CellNotDerived);
    };
    let not_shown = ends.filter(|ends| matches!(ends.start, StartEnd::NotShown { .. }));
    let (Some((theirs, _)), Some(reference_length), Some(ends)) =
        (read.cells.as_ref(), read.reference_length, not_shown)
    else {
        return Err(NotDerived::CellNotDerived);
    };
    let open_at_its_start = matches!(
        positions(theirs.cds_position, "CDS_position", false),
        Ok((None, Some(_)))
    );
    if terms.contains(&CodingTerm::FrameshiftVariant) || open_at_its_start {
        return frameshift_where_no_start_codon_is_shown(
            row,
            ends,
            theirs,
            reference_length,
            read,
            offsets,
        );
    }
    let same_length = allele_of(row.allele).is_ok_and(
        |allele| matches!(allele, Allele::Bases(bases) if bases.len() as u64 == reference_length),
    );
    let at_the_first_residue =
        positions(theirs.protein_position, "Protein_position", false)? == (Some(1), Some(1));
    match peptide_alleles(theirs.amino_acids)? {
        Some((reference, alternate))
            if same_length
                && at_the_first_residue
                && reference.len() == 1
                && alternate.len() == 1
                && reference != b"*" =>
        {
            let named = three_letters(reference[0]).ok_or(NotDerived::CellNotDerived)?;
            if reference == alternate {
                return Ok(format!("p.{named}1="));
            }
            let alternate = three_letters(alternate[0]).ok_or(NotDerived::CellNotDerived)?;
            Ok(format!("p.{named}1{alternate}"))
        }
        _ => Err(NotDerived::CellNotDerived),
    }
}

/// The coding span on which the protein notation of an insertion or a deletion is made.
#[derive(Clone, Copy)]
struct NotationSpan {
    /// The first coding position of the span, one past `end` for an insertion between
    /// two neighbouring positions.
    start: i64,
    end: i64,
    /// Whether the span is that of a deletion moved 3' by the row's `HGVS_OFFSET`.
    moved: bool,
}

impl NotationSpan {
    /// The count of coding positions the span holds, 0 for an insertion between two
    /// neighbouring positions.
    fn length(self) -> i128 {
        i128::from(self.end) - i128::from(self.start) + 1
    }

    /// The residue whose codon holds the span's start, which is the translation start of
    /// the span (`int((cds_start + 2) / 3)`).
    fn first_residue(self) -> Option<usize> {
        usize::try_from(self.start.checked_add(2)? / 3)
            .ok()
            .filter(|residue| *residue >= 1)
    }

    /// Whether an allele of `bases` bases written over the span moves the reading frame:
    /// the positions it takes away less the bases it adds are no multiple of three
    /// (VariationEffect.pm:1564).
    fn moves_the_reading_frame(self, bases: usize) -> bool {
        (bases as i128 - self.length()).abs() % 3 != 0
    }
}

/// The coding span on which both engines make the protein notation of a row whose
/// `CDS_position` cell is `cds_position`, a span of as many coding positions as the
/// record has reference bases (`lengths.reference`):
///
/// - a cell that writes both ends is that span, an insertion's start one past its end
///   (`between_neighbours`), whatever the row's `HGVS_OFFSET`: every place at which one
///   insertion or deletion can be written inside the coding sequence leaves one
///   alternate sequence;
/// - a deletion whose cell is open at its start, because it begins before the coding
///   start, and whose two rows carry one `HGVS_OFFSET` above 0 (`offsets`) is read on
///   the deletion moved 3' by that offset: the span ends at the cell's end plus the
///   offset and holds as many positions as the record has reference bases.
///
/// `None` for any other row: a blank cell, a cell open at its end, a whole span that
/// holds another count of positions than the record has reference bases, as a deletion
/// that takes an intron with it does and an insertion on bases the model alone holds,
/// an open start beside no offset or beside two that differ, an allele that has bases
/// over an open start, and a moved span that leaves the `coding_length` positions of the
/// coding sequence.
fn span_the_notation_reads(
    cds_position: &str,
    offsets: (&str, &str),
    lengths: Lengths,
    between_neighbours: bool,
    coding_length: usize,
) -> Result<Option<NotationSpan>, NotDerived> {
    let whole_number_above_zero = |cell: &str| {
        cell.bytes()
            .all(|digit| digit.is_ascii_digit())
            .then(|| cell.parse::<i64>().ok())
            .flatten()
            .filter(|offset| *offset > 0)
    };
    Ok(
        match positions(cds_position, "CDS_position", between_neighbours)? {
            (Some(start), Some(end)) => Some(NotationSpan {
                start,
                end,
                moved: false,
            })
            .filter(|span| span.length() == i128::from(lengths.reference)),
            (None, Some(printed_end)) if lengths.allele == 0 && lengths.reference > 0 => {
                match (
                    whole_number_above_zero(offsets.0),
                    whole_number_above_zero(offsets.1),
                ) {
                    (Some(offset), Some(other)) if offset == other => printed_end
                        .checked_add(offset)
                        .zip(i64::try_from(coding_length).ok())
                        .and_then(|(end, coding_length)| {
                            let start = end.checked_sub(lengths.reference)?.checked_add(1)?;
                            (start >= 1 && end <= coding_length).then_some(NotationSpan {
                                start,
                                end,
                                moved: true,
                            })
                        }),
                    _ => None,
                }
            }
            _ => None,
        },
    )
}

/// What the frameshift walk of both engines finds on a transcript's sequence
/// (`crates/vep-effects/src/coding.rs`, `hgvsp_fs_peptides` on
/// TranscriptVariationAllele.pm's `_get_fs_peptides`): the first residue at which the
/// translation of the alternate sequence differs from the transcript's peptide.
struct FrameshiftWalk {
    /// The residue the transcript's peptide holds there.
    was: u8,
    /// The residue's number.
    place: usize,
    /// The residue the alternate translation holds there.
    becomes: u8,
    /// The number of the first stop of the alternate translation, `None` where it holds
    /// none.
    stop: Option<usize>,
}

/// The frameshift walk on the sequence `kept` of a transcript at whose coding start no
/// start codon stands, from the residue `first` on, for `allele`, the allele as the
/// transcript reads it, written over `span`. The alternate sequence is the cached coding
/// sequence with the allele in the span's place and the 3' UTR after it
/// (`_get_alternate_cds`, TranscriptVariationAllele.pm:2394-2439); it is translated in
/// whole codons under `code` and held to the transcript's peptide followed by a stop.
/// `None` where the allele does not move the reading frame; where fewer than three
/// coding bases are left, which both engines drop before the 3' UTR, so that the frame
/// read is not the transcript's; where the translation ends before a residue differs;
/// where the first difference lies at or past the reference's stop; and where a codon
/// from the first difference to the first stop holds a letter that is no base, since
/// such a codon may itself be a stop codon and the sequence then gives no count.
fn frameshift_walk(
    kept: &UnshownStartSequence,
    code: &super::GeneticCode,
    allele: &[u8],
    span: NotationSpan,
    first: usize,
) -> Option<FrameshiftWalk> {
    if !span.moves_the_reading_frame(allele.len()) {
        return None;
    }
    let clip = |position: i64| {
        usize::try_from(position)
            .unwrap_or(0)
            .min(kept.coding.len())
    };
    let mut alternate = kept.coding[..clip(span.start - 1)].to_vec();
    alternate.extend_from_slice(allele);
    alternate.extend_from_slice(&kept.coding[clip(span.end)..]);
    if alternate.len() < 3 {
        return None;
    }
    alternate.extend_from_slice(&kept.three_prime_utr);
    let translated: Vec<u8> = alternate
        .chunks_exact(3)
        .map(|codon| code.residue(codon))
        .collect();
    let mut reference = kept.peptide.clone();
    reference.push(b'*');
    let mut place = first;
    let (was, becomes) = loop {
        let at = place.checked_sub(1)?;
        let becomes = *translated.get(at)?;
        let was = *reference.get(at)?;
        if was == b'*' {
            return None;
        }
        if was != becomes {
            break (was, becomes);
        }
        place += 1;
    };
    let stop = translated.iter().position(|residue| *residue == b'*');
    let read = |codon: &[u8]| codon.iter().all(|base| super::rank(*base).is_some());
    if let Some(stop) = stop {
        let mut counted_through = alternate.chunks_exact(3).take(stop).skip(place - 1);
        if !counted_through.all(read) {
            return None;
        }
    }
    Some(FrameshiftWalk {
        was,
        place,
        becomes,
        stop: stop.map(|stop| stop + 1),
    })
}

/// The protein notation of a frameshift on a row on which Ensembl VEP writes a start
/// term and the lost start's notation, on a transcript that shows no start codon and
/// carries no `cds_start_NF`, as vep-rs writes a frameshift's
/// (`crates/vep-effects/src/coding.rs`, the frameshift arm of `hgvsp_format` on
/// `frameshift_walk`). The walk reads the span of `span_the_notation_reads` from the
/// span's first translation position on: the one `Protein_position` writes, or for a
/// deletion moved 3' the residue of the moved span's first codon. The first residue
/// that differs is named with its number and the residue it becomes; `fsTer` and the
/// count of residues from it to the first stop of the alternate translation follow, `?`
/// where the translation holds none. A residue that becomes a stop is written as that
/// stop alone. Not derived: a row without a start term of Ensembl VEP's, a transcript
/// that carries `cds_start_NF`, a record whose sequence is not kept or that is given on
/// the reverse strand, a span that `span_the_notation_reads` does not read, a whole span
/// without a translation start, a walk that names no residue, and a residue outside the
/// twenty.
fn frameshift_where_no_start_codon_is_shown(
    row: &RowFacts<'_>,
    ends: &CodingEnds,
    theirs: &Cells<'_>,
    reference_length: u64,
    read: &WhatWasRead<'_>,
    offsets: (&str, &str),
) -> Result<String, NotDerived> {
    let refused = NotDerived::CellNotDerived;
    let theirs_terms = coding_family(row.ensembl_vep_terms);
    let at_the_start = theirs_terms.contains(&CodingTerm::StartLost)
        || theirs_terms.contains(&CodingTerm::StartRetainedVariant);
    let without_the_flag = matches!(ends.start, StartEnd::NotShown { flagged: false, .. });
    let (Some(kept), true, true, false, Ok(Allele::Bases(allele))) = (
        ends.unshown_start.as_deref(),
        without_the_flag,
        at_the_start,
        read.record_on_the_reverse_strand,
        allele_of(row.allele),
    ) else {
        return Err(refused);
    };
    let code = super::CODES.get(&kept.table).ok_or(refused)?;
    let lengths = Lengths {
        reference: i64::try_from(reference_length).unwrap_or(i64::MAX),
        allele: allele.len() as i64,
    };
    let neighbours = between_neighbours(ends, row.location, lengths)?;
    let span = span_the_notation_reads(
        theirs.cds_position,
        offsets,
        lengths,
        neighbours,
        kept.coding.len(),
    )?
    .ok_or(refused)?;
    let first = match span.moved {
        true => span.first_residue(),
        false => truthy(positions(theirs.protein_position, "Protein_position", neighbours)?.0)
            .and_then(|first| usize::try_from(first).ok()),
    }
    .ok_or(refused)?;
    let allele = as_the_transcript_reads(allele, ends.strand);
    let walk = frameshift_walk(kept, code, &allele, span, first).ok_or(refused)?;
    let place = walk.place;
    let was = three_letters(walk.was).ok_or(refused)?;
    let becomes = three_letters(walk.becomes).ok_or(refused)?;
    if becomes == "Ter" {
        return Ok(format!("p.{was}{place}Ter"));
    }
    let count = walk
        .stop
        .filter(|stop| *stop >= place)
        .map_or("?".to_string(), |stop| (stop + 1 - place).to_string());
    Ok(format!("p.{was}{place}{becomes}fsTer{count}"))
}

/// The two protein notations of a frameshift named at residue 1 on a transcript at whose
/// coding start no start codon stands, where the new reading frame meets a stop codon
/// (`frameshift_counted_at_the_first_residue`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CountedAtTheFirstResidue {
    /// The notation after `p.` up to the count: `<Ref>1<Alt>fsTer`.
    frameshift: String,
    /// The number of the first stop codon of the new reading frame, which is the count
    /// of residues from residue 1 to it.
    count: usize,
}

impl CountedAtTheFirstResidue {
    /// The protein notation vep-rs owes on the row, after `p.`: the frameshift with the
    /// count of residues from residue 1 to the first stop codon of the new reading frame
    /// (`crates/vep-effects/src/coding.rs`, `hgvsp_stop_loss_extra_aa`).
    pub fn vep_rs_notation(&self) -> String {
        format!("{}{}", self.frameshift, self.count)
    }

    /// The protein notation Ensembl VEP writes on the row, after `p.`: the same
    /// frameshift with the unknown count, because its count is not made for the position
    /// before the first residue (TranscriptVariationAllele.pm:2483, handed 0 at :2003,
    /// and :2006).
    pub fn ensembl_vep_notation(&self) -> String {
        format!("{}?", self.frameshift)
    }

    /// Whether `ensembl_vep_cell` and `vep_rs_cell`, the `HGVSp` cells of the two rows,
    /// are the pair: both name one protein before `:p.`, and their descriptions, a pair
    /// of prediction parentheses set aside where both carry one, are
    /// `ensembl_vep_notation` and `vep_rs_notation`.
    pub fn takes(&self, ensembl_vep_cell: &str, vep_rs_cell: &str) -> bool {
        descriptions_of_one_protein(ensembl_vep_cell, vep_rs_cell).is_some_and(|(theirs, ours)| {
            theirs == self.ensembl_vep_notation() && ours == self.vep_rs_notation()
        })
    }
}

/// The descriptions of two `HGVSp` cells that name one protein before `:p.`, each
/// without its prediction parentheses where both carry a pair.
fn descriptions_of_one_protein<'a>(theirs: &'a str, ours: &'a str) -> Option<(&'a str, &'a str)> {
    let (their_protein, theirs) = theirs.split_once(":p.")?;
    let (our_protein, ours) = ours.split_once(":p.")?;
    if their_protein != our_protein {
        return None;
    }
    let predicted = |description: &'a str| {
        description
            .strip_prefix('(')
            .and_then(|description| description.strip_suffix(')'))
    };
    Some(match (predicted(theirs), predicted(ours)) {
        (Some(theirs), Some(ours)) => (theirs, ours),
        _ => (theirs, ours),
    })
}

/// Whether two `HGVSp` cells differ as the pair of `CountedAtTheFirstResidue` differs,
/// read without the transcript's record: both name one protein, Ensembl VEP's
/// description is a frameshift named at residue 1 with the unknown count,
/// `<Ref>1<Alt>fsTer?`, and vep-rs's is the same text with a number in the place of
/// the `?`.
pub fn count_alone_differs_at_the_first_residue(ensembl_vep_cell: &str, vep_rs_cell: &str) -> bool {
    let Some((theirs, ours)) = descriptions_of_one_protein(ensembl_vep_cell, vep_rs_cell) else {
        return false;
    };
    let Some(frameshift) = theirs.strip_suffix('?') else {
        return false;
    };
    let named_at_residue_1 = frameshift.strip_suffix("fsTer").is_some_and(|change| {
        let letters = |name: &[u8]| name.len() == 3 && name.iter().all(u8::is_ascii_alphabetic);
        let change = change.as_bytes();
        change.len() == 7 && letters(&change[..3]) && change[3] == b'1' && letters(&change[4..])
    });
    named_at_residue_1
        && ours.strip_prefix(frameshift).is_some_and(|count| {
            !count.is_empty() && count.bytes().all(|digit| digit.is_ascii_digit())
        })
}

/// The frameshift named at residue 1 on a row of a transcript at whose coding start no
/// start codon stands, flagged `cds_start_NF` or not, on which neither engine reads a
/// lost start. vep-rs counts the residues from residue 1 to the first stop codon of the
/// new reading frame there. Ensembl VEP writes the unknown count for every frameshift
/// whose first changed residue is residue 1, also where the new frame meets a stop
/// codon: `_stop_loss_extra_AA` returns at `return undef unless $ref_var_pos`
/// (TranscriptVariationAllele.pm:2483) for the position 0 the frameshift branch hands it
/// (:2003).
///
/// The derivation: the start end of the record reads `NotShown`; neither row's terms
/// hold `start_lost`, and no more is asked of the terms; both rows print one
/// `CDS_position`; the span is the one of `span_the_notation_reads`, as many positions
/// as the record has reference bases (`offsets` is the `HGVS_OFFSET` cell of Ensembl
/// VEP's row and of vep-rs's, blank for a row that carries none); the walk of
/// `frameshift_walk` opens at the residue of the span's first codon and stops at
/// residue 1; the residue there and the one it becomes are among the twenty amino
/// acids; and the alternate translation, the coding sequence and then the 3' UTR, holds
/// a stop codon, whose number is the count.
///
/// Not derived: a row that is no Transcript's; no record, or one that is not read; a
/// start codon that stands, where vep-rs keeps the unknown count; a lost start on
/// either row, where Ensembl VEP writes the lost start's notation
/// (TranscriptVariationAllele.pm:2161-2164); a structural allele or one that names no
/// bases; no cells or no reference length; a record given on the reverse strand; a
/// `CDS_position` that differs between the rows; a sequence that is not kept
/// (`CodingEnds::unshown_start`); a span that is not read, or that holds another count
/// of positions than the record has reference bases; an allele that does not move the
/// reading frame; fewer than three coding bases left; a walk that passes residue 1, as
/// on a first codon the start exon's phase pads, whose padding letter reads as the
/// unknown residue on both sides; a residue outside the twenty, a stop among them: a
/// residue 1 that becomes a stop is written as that stop alone by both engines
/// (:1997-1999); a new frame that meets no stop codon, where both engines write the
/// unknown count; and a codon before that stop codon that holds a letter that is no
/// base.
pub fn frameshift_counted_at_the_first_residue(
    row: &RowFacts<'_>,
    vep_rs_terms: &str,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
    offsets: (&str, &str),
) -> Result<CountedAtTheFirstResidue, NotDerived> {
    let refused = NotDerived::CellNotDerived;
    if !row.transcript_row {
        return Err(NotDerived::NotATranscriptRow);
    }
    let ends = ends.ok_or(NotDerived::NoVerdict)?;
    match ends.start {
        StartEnd::Unread => return Err(NotDerived::NoVerdict),
        StartEnd::Shown => return Err(refused),
        StartEnd::NotShown { .. } => {}
    }
    let loses_the_start = |terms: &str| coding_family(terms).contains(&CodingTerm::StartLost);
    if loses_the_start(row.ensembl_vep_terms) || loses_the_start(vep_rs_terms) {
        return Err(refused);
    }
    let Allele::Bases(allele) = allele_of(row.allele)? else {
        return Err(refused);
    };
    let (theirs, ours) = read.cells.as_ref().ok_or(NotDerived::NoCells)?;
    let reference_length = read.reference_length.ok_or(NotDerived::NoReferenceLength)?;
    if read.record_on_the_reverse_strand {
        return Err(NotDerived::RecordOnTheReverseStrand);
    }
    if without_length(theirs.cds_position) != without_length(ours.cds_position) {
        return Err(NotDerived::CellsDiffer("CDS_position"));
    }
    let kept = ends
        .unshown_start
        .as_deref()
        .ok_or(NotDerived::SequenceNotKept)?;
    let code = super::CODES.get(&kept.table).ok_or(refused)?;
    let lengths = Lengths {
        reference: i64::try_from(reference_length).unwrap_or(i64::MAX),
        allele: allele.len() as i64,
    };
    let neighbours = between_neighbours(ends, row.location, lengths)?;
    let span = span_the_notation_reads(
        theirs.cds_position,
        offsets,
        lengths,
        neighbours,
        kept.coding.len(),
    )?
    .ok_or(refused)?;
    let first = span.first_residue().ok_or(refused)?;
    let allele = as_the_transcript_reads(allele, ends.strand);
    let walk = frameshift_walk(kept, code, &allele, span, first).ok_or(refused)?;
    let among_the_twenty = |residue: u8| three_letters(residue).filter(|name| *name != "Ter");
    let (1, Some(was), Some(becomes), Some(count)) = (
        walk.place,
        among_the_twenty(walk.was),
        among_the_twenty(walk.becomes),
        walk.stop,
    ) else {
        return Err(refused);
    };
    Ok(CountedAtTheFirstResidue {
        frameshift: format!("{was}1{becomes}fsTer"),
        count,
    })
}

/// The description Ensembl VEP writes for a frameshift on which it reads a lost start,
/// after `p.`, where its frameshift walk stopped at another residue than the one the
/// row's translation span ends on. `owed` is the description vep-rs owes on the row
/// (`protein_notation_where_no_start_codon_is_shown`, after `p.`), which opens with the
/// residue and the number the walk stopped at. The walk moves the notation's start
/// alone (TranscriptVariationAllele.pm:2340, 2360-2361); the lost start then sets `?`
/// in the alternate's place and empties the notation's type (:2161-2164), and a
/// notation without a type whose start is not its end is written with both numbers,
/// `<Ref><n>_?<end>` (:2030-2031). The end is the end of the translation span that
/// Ensembl VEP's `Protein_position` cell writes, which is the span the notation reads
/// on a row that is not shifted 3'. `None` where the walk stopped at the span's end,
/// where Ensembl VEP writes `<Ref><n>?` (:2036), and for a row whose allele, reference
/// length or cells are not read.
pub fn lost_start_range_of_a_frameshift(
    owed: &str,
    row: &RowFacts<'_>,
    ends: &CodingEnds,
    read: &WhatWasRead<'_>,
) -> Option<String> {
    let digits = owed.bytes().skip(3).take_while(u8::is_ascii_digit).count();
    let opening = owed.get(..3 + digits).filter(|_| digits > 0)?;
    let stopped_at: i64 = opening.get(3..)?.parse().ok()?;
    let (theirs, _) = read.cells.as_ref()?;
    let Ok(Allele::Bases(allele)) = allele_of(row.allele) else {
        return None;
    };
    let lengths = Lengths {
        reference: i64::try_from(read.reference_length?).unwrap_or(i64::MAX),
        allele: allele.len() as i64,
    };
    let neighbours = between_neighbours(ends, row.location, lengths).ok()?;
    let (_, end) = positions(theirs.protein_position, "Protein_position", neighbours).ok()?;
    let end = end?;
    (end != stopped_at).then(|| format!("{opening}_?{end}"))
}

/// A class of pair at a coding end: a difference between Ensembl VEP's row and vep-rs's
/// that one reading of Ensembl VEP charges, each named as the class file names it. The
/// order is the order of asking.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CodingEndClass {
    /// Ensembl VEP writes `stop_lost` for a structural deletion over the last three
    /// genomic coding positions, whatever stands there (VariationEffect.pm:1268-1280).
    StructuralDeletionStopLost,
    /// Ensembl VEP reads a stop codon that polyadenylation completes as an incomplete
    /// terminal codon (VariationEffect.pm:1588-1614 with 1240, 1296, 1548).
    StopCodonCompleted,
    /// Ensembl VEP reads the last three coding positions as the stop codon in its
    /// overlap tests (VariationEffect.pm:1373-1381 with 1466-1480; 1400, 1408-1427 and
    /// 1508 read at 1318).
    StopCodonOverlap,
    /// Ensembl VEP reads a stop codon brought behind the last residue as the
    /// transcript's own kept (VariationEffect.pm:1344 and the first clause of 1357).
    StopCodonBehindLastResidue,
    /// Ensembl VEP reads the first three coding positions as the start codon
    /// (VariationEffect.pm:975-986).
    StartCodonOverlap,
}

impl CodingEndClass {
    /// The classes in the order they are asked.
    pub const IN_ORDER: [CodingEndClass; 5] = [
        CodingEndClass::StructuralDeletionStopLost,
        CodingEndClass::StopCodonCompleted,
        CodingEndClass::StopCodonOverlap,
        CodingEndClass::StopCodonBehindLastResidue,
        CodingEndClass::StartCodonOverlap,
    ];

    /// The class's name in a class file.
    pub fn name(self) -> &'static str {
        match self {
            CodingEndClass::StructuralDeletionStopLost => {
                "structural_deletion_stop_lost_without_stop_codon"
            }
            CodingEndClass::StopCodonCompleted => "stop_codon_completed_by_polyadenylation",
            CodingEndClass::StopCodonOverlap => "stop_codon_overlap_without_stop_codon",
            CodingEndClass::StopCodonBehindLastResidue => {
                "stop_codon_behind_last_residue_without_stop_codon"
            }
            CodingEndClass::StartCodonOverlap => "start_codon_overlap_without_start_codon",
        }
    }

    /// The class a class file names `name`.
    pub fn of_name(name: &str) -> Option<CodingEndClass> {
        CodingEndClass::IN_ORDER
            .into_iter()
            .find(|class| class.name() == name)
    }

    /// Whether the class takes the pair of Ensembl VEP's row `row` and vep-rs's row,
    /// whose `Consequence` cell is `vep_rs_terms`. Every class asks that the row is a
    /// Transcript row with both coding ends read, that the terms outside the coding
    /// family are equal on the two rows, that the two rows' terms differ, and that
    /// vep-rs's terms of the coding family are the derived ones (`expected_terms`).
    /// Beside that each asks its premise on the coding ends and its reading of Ensembl
    /// VEP's row:
    ///
    /// - `StructuralDeletionStopLost`: a structural allele; no stop codon stands, with
    ///   the flag or without, or polyadenylation completes one; Ensembl VEP's row
    ///   holds `stop_lost` and the derived set none.
    /// - `StopCodonCompleted`: a sequence allele; polyadenylation completes the stop
    ///   codon; Ensembl VEP read the incomplete codon: its row holds
    ///   `incomplete_terminal_codon_variant`, or the reference side of its
    ///   `Amino_acids` ends on the unknown residue and the larger number of its
    ///   `Protein_position` is one past the whole codons.
    /// - `StopCodonOverlap`: a sequence allele; no stop codon stands and the
    ///   transcript carries no `cds_end_NF`, or polyadenylation completes one; Ensembl
    ///   VEP's row holds a stop term by its overlap tests (`stop_lost` with its
    ///   peptides undefined or the unknown residue in the alternate one;
    ///   `stop_retained_variant` with them undefined, or the alternate one empty or
    ///   holding the unknown residue); the derived set holds neither stop term.
    /// - `StopCodonBehindLastResidue`: a sequence allele; no stop codon stands, with
    ///   the flag or without; Ensembl VEP's row holds `stop_retained_variant` by its
    ///   peptide comparison (the alternate peptide defined, not empty and without the
    ///   unknown residue) though the two peptides do not hold their first stop at one
    ///   place (VariationEffect.pm:1358); the derived set holds no
    ///   `stop_retained_variant`.
    /// - `StartCodonOverlap`: any allele; no start codon is shown and the transcript
    ///   carries no `cds_start_NF`; Ensembl VEP's row holds a start term and the
    ///   derived set none.
    pub fn takes(
        self,
        row: &RowFacts<'_>,
        vep_rs_terms: &str,
        ends: Option<&CodingEnds>,
        read: &WhatWasRead<'_>,
    ) -> bool {
        let Some(ends) = ends else {
            return false;
        };
        let theirs = coding_family(row.ensembl_vep_terms);
        let ours = coding_family(vep_rs_terms);
        if theirs == ours
            || outside_the_coding_family(row.ensembl_vep_terms)
                != outside_the_coding_family(vep_rs_terms)
        {
            return false;
        }
        let Ok(owed) = expected_terms(row, Some(ends), read) else {
            return false;
        };
        let (Some(stop), Ok(allele)) = (stop_codon(ends.stop), allele_of(row.allele)) else {
            return false;
        };
        if ours != owed {
            return false;
        }
        let structural = matches!(allele, Allele::Structural);
        let holds = |terms: &BTreeSet<CodingTerm>, term: CodingTerm| terms.contains(&term);
        // Ensembl VEP's peptides and translation positions, of its own cells.
        let peptides = read
            .cells
            .as_ref()
            .and_then(|(theirs, _)| peptide_alleles(theirs.amino_acids).ok())
            .flatten();
        let unknown = |peptide: &[u8]| peptide.contains(&b'X');
        match self {
            CodingEndClass::StructuralDeletionStopLost => {
                structural
                    && stop != StopCodon::Whole
                    && holds(&theirs, CodingTerm::StopLost)
                    && !holds(&owed, CodingTerm::StopLost)
            }
            CodingEndClass::StopCodonCompleted => {
                let past_the_whole_codons = read.cells.as_ref().is_some_and(|(theirs, _)| {
                    positions(theirs.protein_position, "Protein_position", false).is_ok_and(
                        |(first, last)| first.max(last) == Some(ends.coding_length / 3 + 1),
                    )
                });
                !structural
                    && stop == StopCodon::Completed
                    && (holds(&theirs, CodingTerm::IncompleteTerminalCodonVariant)
                        || (past_the_whole_codons
                            && peptides
                                .as_ref()
                                .is_some_and(|(reference, _)| reference.last() == Some(&b'X'))))
            }
            CodingEndClass::StopCodonOverlap => {
                let lost_by_the_overlap = holds(&theirs, CodingTerm::StopLost)
                    && peptides
                        .as_ref()
                        .is_none_or(|(_, alternate)| unknown(alternate));
                let retained_by_the_overlap = holds(&theirs, CodingTerm::StopRetainedVariant)
                    && peptides
                        .as_ref()
                        .is_none_or(|(_, alternate)| alternate.is_empty() || unknown(alternate));
                !structural
                    && matches!(
                        stop,
                        StopCodon::None { flagged: false } | StopCodon::Completed
                    )
                    && (lost_by_the_overlap || retained_by_the_overlap)
                    && !holds(&owed, CodingTerm::StopLost)
                    && !holds(&owed, CodingTerm::StopRetainedVariant)
            }
            CodingEndClass::StopCodonBehindLastResidue => {
                let stop_place =
                    |peptide: &[u8]| peptide.iter().position(|residue| *residue == b'*');
                let by_the_peptides_without_one_place =
                    peptides.as_ref().is_some_and(|(reference, alternate)| {
                        !alternate.is_empty()
                            && !unknown(alternate)
                            && !(stop_place(reference).is_some()
                                && stop_place(alternate) == stop_place(reference))
                    });
                !structural
                    && matches!(stop, StopCodon::None { .. })
                    && holds(&theirs, CodingTerm::StopRetainedVariant)
                    && by_the_peptides_without_one_place
                    && !holds(&owed, CodingTerm::StopRetainedVariant)
            }
            CodingEndClass::StartCodonOverlap => {
                matches!(ends.start, StartEnd::NotShown { flagged: false, .. })
                    && (holds(&theirs, CodingTerm::StartLost)
                        || holds(&theirs, CodingTerm::StartRetainedVariant))
                    && !holds(&owed, CodingTerm::StartLost)
                    && !holds(&owed, CodingTerm::StartRetainedVariant)
            }
        }
    }
}

impl fmt::Display for CodingEndClass {
    fn fmt(&self, text: &mut fmt::Formatter<'_>) -> fmt::Result {
        text.write_str(self.name())
    }
}

/// The class that takes the pair of Ensembl VEP's row `row` and vep-rs's row, whose
/// `Consequence` cell is `vep_rs_terms`: the first of `CodingEndClass::IN_ORDER` whose
/// `takes` holds, `None` where none takes it. A pair that carries a charged term of
/// both ends is counted under the stop end's class.
pub fn class_of_pair(
    row: &RowFacts<'_>,
    vep_rs_terms: &str,
    ends: Option<&CodingEnds>,
    read: &WhatWasRead<'_>,
) -> Option<CodingEndClass> {
    CodingEndClass::IN_ORDER
        .into_iter()
        .find(|class| class.takes(row, vep_rs_terms, ends, read))
}

/// A class of pair at a stop codon or a start codon that stands, each named as the
/// class file names it. Its rule reads the two rows' terms, the kind of allele and what
/// it reads of the allele's place; `landed_rule_may_take` says on which coding ends it
/// takes beside a cache.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LandedRule {
    StopCodonRemovedWithoutStopCodonTest,
    StopRetainedWithoutStopCodonTest,
    StopRetainedThreeBaseStopCodonTest,
    StopCodonOverlapInsertionRange,
    StopCodonOverlapGenomicSpan,
    StopCodonUnshiftedSpan,
    StartCooccurrenceSwap,
    HgvspStopNamedOnIncompleteCodingSequence,
}

impl LandedRule {
    pub const ALL: [LandedRule; 8] = [
        LandedRule::StopCodonRemovedWithoutStopCodonTest,
        LandedRule::StopRetainedWithoutStopCodonTest,
        LandedRule::StopRetainedThreeBaseStopCodonTest,
        LandedRule::StopCodonOverlapInsertionRange,
        LandedRule::StopCodonOverlapGenomicSpan,
        LandedRule::StopCodonUnshiftedSpan,
        LandedRule::StartCooccurrenceSwap,
        LandedRule::HgvspStopNamedOnIncompleteCodingSequence,
    ];

    /// The rule's name in a class file.
    pub fn name(self) -> &'static str {
        match self {
            LandedRule::StopCodonRemovedWithoutStopCodonTest => {
                "stop_codon_removed_without_stop_codon_test"
            }
            LandedRule::StopRetainedWithoutStopCodonTest => "stop_retained_without_stop_codon_test",
            LandedRule::StopRetainedThreeBaseStopCodonTest => {
                "stop_retained_three_base_stop_codon_test"
            }
            LandedRule::StopCodonOverlapInsertionRange => "stop_codon_overlap_insertion_range",
            LandedRule::StopCodonOverlapGenomicSpan => "stop_codon_overlap_genomic_span",
            LandedRule::StopCodonUnshiftedSpan => "stop_codon_unshifted_span",
            LandedRule::StartCooccurrenceSwap => "start_cooccurrence_swap",
            LandedRule::HgvspStopNamedOnIncompleteCodingSequence => {
                "hgvsp_stop_named_on_incomplete_coding_sequence"
            }
        }
    }

    /// The rule a class file names `name`.
    pub fn of_name(name: &str) -> Option<LandedRule> {
        LandedRule::ALL.into_iter().find(|rule| rule.name() == name)
    }
}

/// Whether a comparison that was given the cache may let `rule` take a pair on the row,
/// beside the rule's own reading of the two rows: where the coding ends of the row's
/// transcript show that the codon the rule speaks of stands. Where it stands every pair
/// of the rule is the rule's as it was, and what the rule reads of the allele's place is
/// its own predicate's; this asks nothing of the row or of what was read beside it. A
/// comparison given no cache does not ask.
///
/// - A rule of a stop codon (`stop_codon_removed_without_stop_codon_test`,
///   `stop_retained_without_stop_codon_test`, `stop_retained_three_base_stop_codon_test`,
///   `stop_codon_overlap_insertion_range`, `stop_codon_overlap_genomic_span`,
///   `stop_codon_unshifted_span`) takes where a whole stop codon stands or
///   polyadenylation completes one, and nothing where no stop codon stands, where the
///   record is not read, and where no record answers for the row.
///   `stop_retained_three_base_stop_codon_test` takes nothing on a completed stop
///   codon: its premise is one triplet at one old stop position read two ways, and the
///   engines read two different triplets there.
///   `stop_retained_without_stop_codon_test` takes only where the genome carries the
///   stop codon whole: nothing on a completed stop codon, and nothing on a whole one of
///   which a base has no genome position, a base an edit of the model adds. On a span
///   without its coordinates over such a codon vep-rs reads the sequence the record
///   gives and writes `stop_lost`, the pair of
///   `stop_codon_removed_without_stop_codon_test`.
/// - `start_cooccurrence_swap` takes where a start codon is shown, and nothing where
///   none is, where the record is not read, and where no record answers.
/// - `hgvsp_stop_named_on_incomplete_coding_sequence` reads as it does without the
///   coding ends, but for a completed stop codon, on which vep-rs owes a notation.
pub fn landed_rule_may_take(
    rule: LandedRule,
    _row: &RowFacts<'_>,
    ends: Option<&CodingEnds>,
    _read: &WhatWasRead<'_>,
) -> bool {
    let stop = ends.and_then(|ends| stop_codon(ends.stop));
    match rule {
        LandedRule::HgvspStopNamedOnIncompleteCodingSequence => stop != Some(StopCodon::Completed),
        LandedRule::StartCooccurrenceSwap => {
            stop.is_some() && ends.is_some_and(|ends| ends.start == StartEnd::Shown)
        }
        LandedRule::StopRetainedThreeBaseStopCodonTest => stop == Some(StopCodon::Whole),
        LandedRule::StopRetainedWithoutStopCodonTest => {
            stop == Some(StopCodon::Whole)
                && ends.is_some_and(|ends| ends.last_codon.iter().all(Option::is_some))
        }
        LandedRule::StopCodonRemovedWithoutStopCodonTest
        | LandedRule::StopCodonOverlapInsertionRange
        | LandedRule::StopCodonOverlapGenomicSpan
        | LandedRule::StopCodonUnshiftedSpan => {
            matches!(stop, Some(StopCodon::Whole | StopCodon::Completed))
        }
    }
}

/// The class words an `Allele` cell holds on a row of a structural allele: the Sequence
/// Ontology terms of the structural variant classes
/// (`Bio::EnsEMBL::Variation::Utils::Constants`).
pub const STRUCTURAL_CLASS_WORDS: [&str; 30] = [
    "insertion",
    "deletion",
    "duplication",
    "tandem_duplication",
    "inversion",
    "copy_number_variation",
    "copy_number_gain",
    "copy_number_loss",
    "mobile_element_insertion",
    "mobile_element_deletion",
    "novel_sequence_insertion",
    "Alu_insertion",
    "Alu_deletion",
    "LINE1_insertion",
    "LINE1_deletion",
    "SVA_insertion",
    "SVA_deletion",
    "HERV_insertion",
    "HERV_deletion",
    "chromosome_breakpoint",
    "interchromosomal_breakpoint",
    "intrachromosomal_breakpoint",
    "translocation",
    "interchromosomal_translocation",
    "intrachromosomal_translocation",
    "complex_structural_alteration",
    "complex_chromosomal_rearrangement",
    "tandem_repeat",
    "short_tandem_repeat_variation",
    "loss_of_heterozygosity",
];

/// A row's allele: the bases of a sequence allele, none for `-`, or a structural
/// allele, whose cell is one of `STRUCTURAL_CLASS_WORDS`. Any other cell is no allele a
/// row is derived for: one with an ambiguity letter, the unknown base, a lower-case
/// base or a star is a sequence allele that names no bases, for which both engines
/// leave their sequence tests early (VariationEffect.pm:1001, 1038, 1498).
enum Allele<'a> {
    Bases(&'a [u8]),
    Structural,
}

fn allele_of(cell: &str) -> Result<Allele<'_>, NotDerived> {
    if cell == "-" {
        return Ok(Allele::Bases(b""));
    }
    if STRUCTURAL_CLASS_WORDS.contains(&cell) {
        return Ok(Allele::Structural);
    }
    let letters = cell.as_bytes();
    if letters.is_empty()
        || !letters
            .iter()
            .all(|letter| matches!(letter, b'A' | b'C' | b'G' | b'T'))
    {
        return Err(NotDerived::AlleleNamesNoBases);
    }
    Ok(Allele::Bases(letters))
}

/// The alleles' lengths and the class Ensembl VEP reads from them: of one length, an
/// insertion (the allele is the longer) or a deletion (the reference bases are).
#[derive(Clone, Copy)]
struct Lengths {
    reference: i64,
    allele: i64,
}

impl Lengths {
    fn same(self) -> bool {
        self.reference == self.allele
    }

    fn insertion(self) -> bool {
        self.reference < self.allele
    }

    fn deletion(self) -> bool {
        self.reference > self.allele
    }

    /// An insertion that stands in for no reference base: its Location is the two bases
    /// it stands between.
    fn between_two_bases(self) -> bool {
        self.reference == 0
    }
}

fn blank(cell: &str) -> bool {
    cell.is_empty() || cell == "-"
}

/// A position cell without the `/` and length that `--total_length` appends.
fn without_length(cell: &str) -> &str {
    cell.split_once('/').map_or(cell, |(position, _)| position)
}

/// The start and the end a position cell writes (`format_coords`, Utils.pm:141-166):
/// one number is both, `?` is an undefined end, and a blank cell is two. Of two numbers
/// the cell writes the smaller first, so an insertion between two neighbouring bases,
/// whose start lies one past its end, reads its start from the second
/// (`between_neighbours`).
fn positions(
    cell: &str,
    name: &'static str,
    between_neighbours: bool,
) -> Result<(Option<i64>, Option<i64>), NotDerived> {
    let cell = without_length(cell);
    if blank(cell) {
        return Ok((None, None));
    }
    let number = |text: &str| match text {
        "?" => Ok(None),
        digits if !digits.is_empty() && digits.bytes().all(|digit| digit.is_ascii_digit()) => {
            digits
                .parse::<i64>()
                .map(Some)
                .map_err(|_| NotDerived::CellNotRead(name))
        }
        _ => Err(NotDerived::CellNotRead(name)),
    };
    match cell.split_once('-') {
        None => {
            let position = number(cell)?;
            Ok((position, position))
        }
        Some((first, second)) => match (number(first)?, number(second)?) {
            (Some(first), Some(second)) if between_neighbours => {
                Ok((Some(first.max(second)), Some(first.min(second))))
            }
            (None, None) => Err(NotDerived::CellNotRead(name)),
            ends => Ok(ends),
        },
    }
}

/// The reference and the alternate peptide an `Amino_acids` cell writes
/// (`pep_allele_string`): one text is both, `-` on a side is the empty peptide, and a
/// blank cell is Ensembl VEP's empty list of `_get_peptide_alleles`
/// (VariationEffect.pm:772-798).
type Alleles = (Vec<u8>, Vec<u8>);

fn peptide_alleles(cell: &str) -> Result<Option<Alleles>, NotDerived> {
    if blank(cell) {
        return Ok(None);
    }
    let side = |text: &str| match text {
        "-" => Ok(Vec::new()),
        residues
            if !residues.is_empty()
                && residues
                    .bytes()
                    .all(|residue| residue.is_ascii_uppercase() || residue == b'*') =>
        {
            Ok(residues.as_bytes().to_vec())
        }
        _ => Err(NotDerived::CellNotRead("Amino_acids")),
    };
    match cell.split_once('/') {
        None => Ok(Some((side(cell)?, side(cell)?))),
        Some((reference, alternate)) => Ok(Some((side(reference)?, side(alternate)?))),
    }
}

/// The reference and the alternate codon a `Codons` cell writes
/// (`display_codon_allele_string`), in upper case: the cell writes each in lower case
/// with the allele's own bases in upper case, and `-` for an empty one.
fn codon_alleles(cell: &str) -> Result<Option<Alleles>, NotDerived> {
    if blank(cell) {
        return Ok(None);
    }
    let side = |text: &str| match text {
        "-" => Ok(Vec::new()),
        bases if !bases.is_empty() && bases.bytes().all(|base| base.is_ascii_alphabetic()) => {
            Ok(bases.to_ascii_uppercase().into_bytes())
        }
        _ => Err(NotDerived::CellNotRead("Codons")),
    };
    let (reference, alternate) = cell
        .split_once('/')
        .ok_or(NotDerived::CellNotRead("Codons"))?;
    Ok(Some((side(reference)?, side(alternate)?)))
}

/// `Bio::EnsEMBL::Variation::Utils::VariationEffect::overlap`.
fn overlap(first_start: i64, first_end: i64, second_start: i64, second_end: i64) -> bool {
    first_end >= second_start && first_start <= second_end
}

/// Perl's truth of a number: defined and not zero.
fn truthy(value: Option<i64>) -> Option<i64> {
    value.filter(|value| *value != 0)
}

/// The cDNA position of a genomic base, `None` for a base no mapper pair holds.
fn cdna_of(model: &CodingModel, genomic: i64) -> Option<i128> {
    model
        .pairs
        .iter()
        .find(|pair| pair.2 <= genomic && genomic <= pair.3)
        .map(|&(from_start, _, to_start, to_end, orientation)| {
            let (from_start, genomic) = (i128::from(from_start), i128::from(genomic));
            if orientation >= 0 {
                from_start + (genomic - i128::from(to_start))
            } else {
                from_start + (i128::from(to_end) - genomic)
            }
        })
}

/// Whether the allele is an insertion between two neighbouring bases of the model: it
/// stands in for no reference base, and the two genome bases it stands between map to
/// neighbouring cDNA positions or one of them to none. An insertion between two genome
/// bases whose cDNA positions lie further apart stands on the bases the model holds
/// between them, and its cells write that span.
fn between_neighbours(
    ends: &CodingEnds,
    (lowest, highest): (i64, i64),
    lengths: Lengths,
) -> Result<bool, NotDerived> {
    if !lengths.between_two_bases() {
        return Ok(false);
    }
    let Some(model) = ends.model.as_deref() else {
        // Without the coordinates the place of an insertion on a model that holds
        // bases the genome lacks cannot be told from the cells.
        if ends.edited_model {
            return Err(NotDerived::InsertionOnAnEditedModel);
        }
        return Ok(true);
    };
    Ok(match (cdna_of(model, lowest), cdna_of(model, highest)) {
        (Some(one), Some(other)) => (one - other).abs() == 1,
        _ => true,
    })
}

/// Ensembl VEP's `within_cds` behind `coding_unknown` (VariationEffect.pm:1616-1639), as
/// vep-rs evaluates it: the span projects onto a coding position, which a `CDS_position`
/// cell that holds a number shows; or, where the cell is blank, a span of reference
/// bases holds the genomic base of a coding position between two ends that are none;
/// or the span meets an intron of the record whose ends lie at most twelve bases apart
/// and the coding region. An insertion whose cell is blank stands on no coding position.
/// A span of reference bases that holds the record's whole span gives no
/// `coding_sequence_variant` (VariationEffect.pm:1621, `complete_overlap_feature`), as
/// vep-rs has it (`crates/vep-effects/src/consequences.rs`, `complete_overlap`).
fn within_cds(
    cds_position: &str,
    ends: &CodingEnds,
    (lowest, highest): (i64, i64),
    lengths: Lengths,
) -> Result<bool, NotDerived> {
    if !lengths.between_two_bases()
        && ends
            .span
            .is_some_and(|(start, end)| lowest <= start && end <= highest)
    {
        return Ok(false);
    }
    if without_length(cds_position)
        .bytes()
        .any(|letter| letter.is_ascii_digit())
    {
        return Ok(true);
    }
    let model = ends.model.as_deref().ok_or(NotDerived::ModelNotKept)?;
    let holds_a_coding_base = || {
        model.pairs.iter().any(|pair| {
            let (first, last) = (pair.0.max(model.coding_start), pair.1.min(model.coding_end));
            first <= last && {
                let image = |cdna: i64| {
                    let along = i128::from(cdna) - i128::from(pair.0);
                    match pair.4 >= 0 {
                        true => i128::from(pair.2) + along,
                        false => i128::from(pair.3) - along,
                    }
                };
                let (one, other) = (image(first), image(last));
                i128::from(highest) >= one.min(other) && i128::from(lowest) <= one.max(other)
            }
        })
    };
    if !lengths.between_two_bases() && holds_a_coding_base() {
        return Ok(true);
    }
    // Ensembl VEP reads an insertion's span with its start one past its end.
    let (start, end) = match lengths.between_two_bases() {
        true => (highest, lowest),
        false => (lowest, highest),
    };
    if let Some((first, last)) = model.coding_region {
        if model
            .short_introns
            .iter()
            .any(|&(intron_start, intron_end)| overlap(start, end, intron_start, intron_end))
        {
            return Ok(overlap(start, end, first, last));
        }
    }
    Ok(false)
}

/// Ensembl VEP's `partial_codon` (VariationEffect.pm:1588-1614): the codon that holds
/// the span's translation start is the coding sequence's last and is shorter than three
/// bases.
fn partial_codon(translation_start: Option<i64>, coding_length: i64) -> bool {
    translation_start.is_some_and(|start| coding_length % 3 > 0 && start == coding_length / 3 + 1)
}

/// The values the predicates read of one allele: from the cells `equal_cells` holds, or
/// computed on the completed coding sequence.
struct Values<'a> {
    partial_codon: bool,
    lengths: Lengths,
    cds: (Option<i64>, Option<i64>),
    /// `_get_peptide_alleles`; `None` is Perl's empty list.
    peptides: Option<&'a Alleles>,
    /// The two codons before `_get_codon_alleles` withholds them behind a frameshift.
    codons: Option<&'a Alleles>,
    /// Whether the alternate allele's peptide is undefined, which the cells show as no
    /// peptide alleles.
    alternate_peptide_undefined: bool,
}

/// The four terms a coding end decides.
#[derive(Clone, Copy)]
struct Charged {
    stop_lost: bool,
    stop_retained: bool,
    start_lost: bool,
    start_retained: bool,
}

/// The terms of the coding family that vep-rs emits (`PerlCodingEval::terms`) from the
/// four charged terms and Ensembl VEP's other predicates, each evaluated on `values` as
/// the line cited reads it. `within_cds` is asked only where `coding_unknown` reads it.
fn terms_of(
    values: &Values<'_>,
    charged: Charged,
    within_cds: impl FnOnce() -> Result<bool, NotDerived>,
) -> Result<BTreeSet<CodingTerm>, NotDerived> {
    let Charged {
        stop_lost,
        stop_retained,
        start_lost,
        start_retained,
    } = charged;
    let partial = values.partial_codon;
    let lengths = values.lengths;
    let stop = |peptide: &[u8]| peptide.contains(&b'*');
    let unknown = |peptide: &[u8]| peptide.contains(&b'X');
    let allele_length = i128::from(lengths.allele);
    let spanned = match values.cds {
        (Some(start), Some(end)) => Some(i128::from(end) - i128::from(start) + 1),
        _ => None,
    };

    // `frameshift` (VariationEffect.pm:1541-1565): not on a partial codon (1548), not
    // beside `stop_retained_variant` (1549), both CDS ends defined (1551), a reference
    // peptide that does not open with a stop (1554), and an allele whose length less the
    // CDS span's is no multiple of three (1564).
    let frameshift = !partial
        && !stop_retained
        && !values
            .peptides
            .is_some_and(|(reference, _)| reference.first() == Some(&b'*'))
        && spanned.is_some_and(|spanned| (allele_length - spanned).abs() % 3 != 0);

    // `stop_gained` (VariationEffect.pm:1211-1232): not beside `stop_retained_variant`
    // (1221) or `stop_lost` (1222); the alternate peptide holds a stop and the reference
    // none (1228).
    let stop_gained = !stop_retained
        && !stop_lost
        && values
            .peptides
            .is_some_and(|(reference, alternate)| stop(alternate) && !stop(reference));

    // `_get_codon_alleles` (VariationEffect.pm:812-830): no codons behind a frameshift
    // (815).
    let codons = values.codons.filter(|_| !frameshift);

    // `inframe_insertion` (VariationEffect.pm:1100-1128): not beside `start_lost`
    // (1109); the alternate codon is the longer (1111); not beside
    // `start_retained_variant` where the alternate peptide ends on the reference's
    // (1119); not a lone stop on both sides (1121); the alternate peptide, cut after its
    // first stop, begins or ends with the reference's (1126-1127).
    let inframe_insertion = !start_lost
        && codons.is_some_and(|(reference, alternate)| alternate.len() > reference.len())
        && values.peptides.is_some_and(|(reference, alternate)| {
            let cut = match alternate.iter().position(|residue| *residue == b'*') {
                Some(stop) => &alternate[..=stop],
                None => &alternate[..],
            };
            !(start_retained && alternate.ends_with(reference))
                && !(reference == b"*" && alternate == b"*")
                && (cut.starts_with(reference) || cut.ends_with(reference))
        });

    // `inframe_deletion` (VariationEffect.pm:1161-1186): not on a partial codon (1168);
    // the alternate codon is the shorter (1173); the reference peptide is not the lone
    // stop (1175); the reference codon begins or ends with the alternate, or with the
    // bases they share at the start and then at the end cut off the alternate is empty
    // and what is left of the reference is whole codons (1178-1185).
    let inframe_deletion = !partial
        && codons.is_some_and(|(reference, alternate)| {
            alternate.len() < reference.len()
                && !values.peptides.is_some_and(|(peptide, _)| peptide == b"*")
                && (reference.starts_with(alternate) || reference.ends_with(alternate) || {
                    let shorter = alternate.len();
                    let shared_start = reference
                        .iter()
                        .zip(alternate)
                        .take_while(|(one, other)| one == other)
                        .count();
                    let shared_end = reference
                        .iter()
                        .rev()
                        .zip(alternate.iter().rev())
                        .take(shorter - shared_start)
                        .take_while(|(one, other)| one == other)
                        .count();
                    shared_start + shared_end == shorter && (reference.len() - shorter) % 3 == 0
                })
        });

    // `missense_variant` (VariationEffect.pm:1085-1098): not beside `start_lost`,
    // `stop_lost`, `stop_gained`, a partial codon or `stop_retained_variant`
    // (1090-1094); the peptides differ and are of one length (1097).
    let missense = !(start_lost || stop_lost || stop_gained || partial || stop_retained)
        && values.peptides.is_some_and(|(reference, alternate)| {
            reference != alternate && reference.len() == alternate.len()
        });

    // `protein_altering_variant` (VariationEffect.pm:376-393): the peptides differ in
    // length, neither opens with a stop, the alternate neither begins nor ends with the
    // reference, and the row is no `inframe_deletion`, `start_lost` or frameshift.
    let protein_altering = values.peptides.is_some_and(|(reference, alternate)| {
        reference.len() != alternate.len()
            && reference.first() != Some(&b'*')
            && alternate.first() != Some(&b'*')
            && !(alternate.starts_with(reference) || alternate.ends_with(reference))
    }) && !inframe_deletion
        && !start_lost
        && !frameshift;

    // `synonymous_variant` (VariationEffect.pm:1077-1083): a reference peptide that the
    // alternate equals, not beside `stop_retained_variant` or `start_retained_variant`,
    // and no unknown residue in either.
    let synonymous = values.peptides.is_some_and(|(reference, alternate)| {
        !reference.is_empty() && reference == alternate && !unknown(reference)
    }) && !stop_retained
        && !start_retained;

    // `coding_unknown` (VariationEffect.pm:1616-1639): within the coding sequence; the
    // alternate peptide or the peptide alleles undefined, or an unknown residue in
    // either peptide (1627-1630); and none of the frameshift, `inframe_deletion`,
    // `protein_altering_variant`, the two start terms and the two stop terms
    // (1633-1635), each as the predicate reads whatever the allele's class.
    let coding_unknown = (values.alternate_peptide_undefined
        || values
            .peptides
            .is_none_or(|(reference, alternate)| unknown(reference) || unknown(alternate)))
        && !(frameshift
            || inframe_deletion
            || protein_altering
            || start_retained
            || start_lost
            || stop_retained
            || stop_lost)
        && within_cds()?;

    // An insertion that stands on bases the model holds and the genome lacks is asked
    // the predicate of the class its CDS span gives it
    // (`stands_on_as_many_model_bases`, `stands_on_more_model_bases`).
    let on_as_many_model_bases =
        lengths.insertion() && spanned.is_some_and(|spanned| allele_length == spanned);
    let on_more_model_bases =
        lengths.insertion() && spanned.is_some_and(|spanned| allele_length < spanned);

    Ok([
        (CodingTerm::StopGained, stop_gained),
        (CodingTerm::FrameshiftVariant, frameshift && !lengths.same()),
        (CodingTerm::StopLost, stop_lost),
        (CodingTerm::StartLost, start_lost),
        (
            CodingTerm::InframeInsertion,
            inframe_insertion && lengths.insertion(),
        ),
        (
            CodingTerm::InframeDeletion,
            inframe_deletion && (lengths.deletion() || on_more_model_bases),
        ),
        (
            CodingTerm::MissenseVariant,
            missense && (lengths.same() || on_as_many_model_bases),
        ),
        (CodingTerm::ProteinAlteringVariant, protein_altering),
        (CodingTerm::IncompleteTerminalCodonVariant, partial),
        (CodingTerm::StartRetainedVariant, start_retained),
        (CodingTerm::StopRetainedVariant, stop_retained),
        (CodingTerm::SynonymousVariant, synonymous),
        (CodingTerm::CodingSequenceVariant, coding_unknown),
    ]
    .into_iter()
    .filter_map(|(term, holds)| holds.then_some(term))
    .collect())
}

/// The two start terms. Every arm of Ensembl VEP's `start_lost` and
/// `start_retained_variant` stands behind `_overlaps_start_codon`
/// (VariationEffect.pm:862 and 957 on 965-990), which vep-rs answers "no" where no
/// start codon stands: both terms are false there. Where one stands the engines read
/// the start by one rule and each term is Ensembl VEP's own; a row on which Ensembl VEP
/// writes both is not derived, since vep-rs keeps one of the two there by another rule.
fn start_terms(start: StartEnd, theirs: &BTreeSet<CodingTerm>) -> Result<(bool, bool), NotDerived> {
    match start {
        StartEnd::NotShown { .. } => Ok((false, false)),
        StartEnd::Shown | StartEnd::Unread => {
            let lost = theirs.contains(&CodingTerm::StartLost);
            let retained = theirs.contains(&CodingTerm::StartRetainedVariant);
            if lost && retained {
                return Err(NotDerived::BothStartTermsWhereAStartCodonStands);
            }
            Ok((lost, retained))
        }
    }
}

/// A structural allele's row: Ensembl VEP's own terms of the family, less `stop_lost`
/// where the span takes no base of a stop codon and less the two start terms where no
/// start codon stands. Ensembl VEP's structural arm of `stop_lost` holds a deletion's
/// span against the last three genomic coding positions and reads no flag and no base
/// (VariationEffect.pm:1268-1280); vep-rs writes the term only where a stop codon
/// stands and the span holds a genomic base of it, the tail's alone where the tail
/// begins it (`crates/vep-effects/src/sv/mod.rs`, `structural_stop_lost`). The
/// structural start terms stand behind the start codon overlap (VariationEffect.pm:
/// 886-896 behind 862). No other term of a structural row reads one of the four
/// (VariationEffect.pm:1643), and no cell is read.
fn structural_terms(
    row: &RowFacts<'_>,
    ends: &CodingEnds,
    stop: StopCodon,
) -> Result<BTreeSet<CodingTerm>, NotDerived> {
    let mut owed = coding_family(row.ensembl_vep_terms);
    // Where a start codon stands and Ensembl VEP writes both start terms, vep-rs keeps
    // one of the two by another rule.
    start_terms(ends.start, &owed)?;
    let (lowest, highest) = row.location;
    let holds_a_tail_base = ends
        .tail
        .iter()
        .flatten()
        .any(|base| lowest <= *base && *base <= highest);
    let takes_a_stop_codon_base = match stop {
        StopCodon::Whole => true,
        StopCodon::Completed => holds_a_tail_base,
        StopCodon::None { .. } => false,
    };
    if !takes_a_stop_codon_base {
        owed.remove(&CodingTerm::StopLost);
    }
    if matches!(ends.start, StartEnd::NotShown { .. }) {
        owed.remove(&CodingTerm::StartLost);
        owed.remove(&CodingTerm::StartRetainedVariant);
    }
    Ok(owed)
}

/// The cells of the two rows, each held to be equal on both: a rule reads no cell of
/// vep-rs alone. `Amino_acids` and `Codons` may also differ as the class
/// `alternate_codon_past_the_coding_end` documents, vep-rs's alternate codon ending where
/// the coding sequence ends and Ensembl VEP's reading on into the 3' UTR: the two shapes
/// are that class's own, asked as its rule asks them, with the first residue of the
/// reference's `Protein_position`. Any other difference of the two cells is refused.
fn equal_cells<'a>(
    theirs: &Cells<'a>,
    ours: &Cells<'_>,
    positions_alone: bool,
) -> Result<(), NotDerived> {
    let same = |one: &str, other: &str| one == other || (blank(one) && blank(other));
    for (name, one, other) in [
        ("cDNA_position", theirs.cdna_position, ours.cdna_position),
        ("CDS_position", theirs.cds_position, ours.cds_position),
        (
            "Protein_position",
            theirs.protein_position,
            ours.protein_position,
        ),
    ] {
        if !same(without_length(one), without_length(other)) {
            return Err(NotDerived::CellsDiffer(name));
        }
    }
    if positions_alone {
        return Ok(());
    }
    if !same(theirs.amino_acids, ours.amino_acids)
        && !crate::coding_end::residues_read_past_the_coding_end(
            theirs.codons,
            ours.codons,
            theirs.amino_acids,
            ours.amino_acids,
        )
    {
        return Err(NotDerived::CellsDiffer("Amino_acids"));
    }
    let first_residue = without_length(theirs.protein_position)
        .split('-')
        .next()
        .unwrap_or_default();
    if !same(theirs.codons, ours.codons)
        && !crate::coding_end::codons_read_past_the_coding_end(
            first_residue,
            theirs.codons,
            ours.codons,
        )
    {
        return Err(NotDerived::CellsDiffer("Codons"));
    }
    Ok(())
}

/// The terms of a sequence allele's row on a transcript whose stop codon no
/// polyadenylation completes, read on Ensembl VEP's cells, which `equal_cells` holds to
/// vep-rs's.
fn on_the_cells(
    row: &RowFacts<'_>,
    ends: &CodingEnds,
    stop_end: StopCodon,
    (theirs, ours): &(Cells<'_>, Cells<'_>),
    lengths: Lengths,
) -> Result<BTreeSet<CodingTerm>, NotDerived> {
    equal_cells(theirs, ours, false)?;
    let neighbours = between_neighbours(ends, row.location, lengths)?;
    let cds = positions(theirs.cds_position, "CDS_position", neighbours)?;
    let protein = positions(theirs.protein_position, "Protein_position", neighbours)?;
    let peptides = peptide_alleles(theirs.amino_acids)?;
    let codons = codon_alleles(theirs.codons)?;
    let partial = partial_codon(protein.0, ends.coding_length);
    let theirs_terms = coding_family(row.ensembl_vep_terms);
    let stop = |peptide: &[u8]| peptide.iter().position(|residue| *residue == b'*');
    let unknown = |peptide: &[u8]| peptide.contains(&b'X');

    // The stop end. Ensembl VEP's `stop_lost` (VariationEffect.pm:1234-1288) and
    // `stop_retained` (1290-1324) each read the peptides where the alternate peptide
    // is defined and holds no unknown residue (1259-1260; 1313-1315 on 1326-1361), and
    // the overlap of the span with the last three coding positions otherwise (1263 on
    // 1363-1385 and 1433-1484; 1318 on 1389-1431 and 1488-1539).
    let (stop_lost, stop_retained) = match stop_end {
        // No stop codon stands: every overlap test is false, and of the peptide
        // comparison of `ref_eq_alt_sequence` the two tests that see no stop codon in
        // the reference (1344, the first clause of 1357) are false too. What stands is
        // the stop a peptide shows inside the coding sequence: lost where the
        // reference peptide holds one and the alternate none (1260), kept where both
        // hold their first at one place (1358).
        StopCodon::None { .. } => {
            let lost = !partial
                && peptides.as_ref().is_some_and(|(reference, alternate)| {
                    !unknown(alternate) && stop(alternate).is_none() && stop(reference).is_some()
                });
            let retained = !partial
                && !lost
                && protein.0.is_some()
                && protein.1.is_some()
                && peptides.as_ref().is_some_and(|(reference, alternate)| {
                    !alternate.is_empty()
                        && !unknown(alternate)
                        && stop(reference).is_some()
                        && stop(alternate) == stop(reference)
                });
            (lost, retained)
        }
        // A stop codon stands: the engines read the end by one rule where the peptides
        // decide both terms, and where the allele lies clear of the stop codon no
        // overlap test holds in either; each term is then Ensembl VEP's own. Elsewhere
        // vep-rs reads the edited sequence, which the cells do not show, and the row
        // is not derived.
        StopCodon::Whole | StopCodon::Completed => {
            let by_the_peptides = peptides
                .as_ref()
                .is_some_and(|(_, alternate)| !alternate.is_empty() && !unknown(alternate));
            if !by_the_peptides && !clear_of_the_stop_codon(ends, row.location, lengths) {
                return Err(NotDerived::AtAStandingStopCodon);
            }
            (
                theirs_terms.contains(&CodingTerm::StopLost),
                theirs_terms.contains(&CodingTerm::StopRetainedVariant),
            )
        }
    };
    let (start_lost, start_retained) = start_terms(ends.start, &theirs_terms)?;
    let terms = terms_of(
        &Values {
            partial_codon: partial,
            lengths,
            cds,
            peptides: peptides.as_ref(),
            codons: codons.as_ref(),
            alternate_peptide_undefined: peptides.is_none(),
        },
        Charged {
            stop_lost,
            stop_retained,
            start_lost,
            start_retained,
        },
        || within_cds(theirs.cds_position, ends, row.location, lengths),
    )?;
    at_a_standing_start_codon(row, ends, &theirs_terms, &terms)?;
    Ok(terms)
}

/// Whether the allele lies clear of a whole stop codon: its Location holds none
/// of the stop codon's own genomic bases and none of the last three genomic positions
/// of the coding region, and an insertion stands more than its own length away from
/// each of them, the reach of Ensembl VEP's genomic overlap for an insertion
/// (VariationEffect.pm:1408-1415). Positions the reading does not hold are not clear of.
fn clear_of_the_stop_codon(
    ends: &CodingEnds,
    (lowest, highest): (i64, i64),
    lengths: Lengths,
) -> bool {
    let (Some(own), Some((first, last))) = (
        ends.last_codon.into_iter().collect::<Option<Vec<i64>>>(),
        ends.coding_region_end,
    ) else {
        return false;
    };
    let reach = match lengths.between_two_bases() {
        true => lengths.allele,
        false => 0,
    };
    own.into_iter()
        .chain(first..=last)
        .all(|base| base < lowest.saturating_sub(reach) || highest.saturating_add(reach) < base)
}

/// Where a start codon stands the two start terms are Ensembl VEP's own, which holds
/// while the terms they depend on are too: Ensembl VEP's `start_lost` asks
/// `inframe_insertion` and `inframe_deletion` (VariationEffect.pm:870), which read no
/// codon behind a frameshift (815), which is withheld behind `stop_retained_variant`
/// (1549). A row whose allele meets the first codon, and on which `stop_lost`,
/// `stop_retained_variant` or `frameshift_variant` is derived otherwise than Ensembl VEP
/// writes it, is not derived. The allele meets the first codon where its Location holds
/// a genomic base of the first three coding positions, or where the reading holds no
/// genomic base for one of them.
fn at_a_standing_start_codon(
    row: &RowFacts<'_>,
    ends: &CodingEnds,
    theirs: &BTreeSet<CodingTerm>,
    owed: &BTreeSet<CodingTerm>,
) -> Result<(), NotDerived> {
    if ends.start != StartEnd::Shown {
        return Ok(());
    }
    let (lowest, highest) = row.location;
    let meets = ends
        .first_codon
        .iter()
        .any(|base| base.is_none_or(|base| lowest <= base && base <= highest));
    let moved = [
        CodingTerm::StopLost,
        CodingTerm::StopRetainedVariant,
        CodingTerm::FrameshiftVariant,
    ]
    .iter()
    .any(|term| theirs.contains(term) != owed.contains(term));
    if meets && moved {
        return Err(NotDerived::AtAStandingStartCodon);
    }
    Ok(())
}

/// A whole number as `i64`, the nearest one where it lies past either end.
fn nearest(number: i128) -> i64 {
    i64::try_from(number).unwrap_or(if number < 0 { i64::MIN } else { i64::MAX })
}

/// Perl's `substr($text, $offset, $length)`: `None` past the end of the text, clipped
/// otherwise.
fn perl_substr(text: &[u8], offset: i64, length: i64) -> Option<&[u8]> {
    let offset = usize::try_from(offset)
        .ok()
        .filter(|offset| *offset <= text.len())?;
    let end = match usize::try_from(length) {
        Ok(length) => offset.saturating_add(length).min(text.len()),
        Err(_) => usize::try_from(text.len() as i64 + length)
            .unwrap_or(0)
            .max(offset),
    };
    Some(&text[offset..end])
}

/// Perl's `substr($text, $offset, $length) = $replacement`, clipped as Perl clips it.
fn perl_substr_assign(text: &mut Vec<u8>, offset: i64, length: i64, replacement: &[u8]) {
    let offset = usize::try_from(offset).unwrap_or(0).min(text.len());
    let end = offset
        .saturating_add(usize::try_from(length).unwrap_or(0))
        .min(text.len());
    text.splice(offset..end, replacement.iter().copied());
}

/// What the stop codon edit test finds (`PerlCodingEval::edited_stop_codon` on
/// VariationEffect.pm:1521-1529).
enum Edited {
    /// The test is not made: the allele changes no length, or its span lacks a cDNA
    /// coordinate at one end or a CDS start (VariationEffect.pm:1508).
    NotMade,
    /// The edited sequence is shorter than the coding sequence (VariationEffect.pm:1524).
    Shorter,
    /// The three bases the edited sequence holds at the old stop position, and whether
    /// every codon before them reads the reference's residue.
    Read { codon: [u8; 3], residues_kept: bool },
}

/// What the sequence a record gives shows at the old stop position, for a span whose
/// edit test cannot be made (`PerlCodingEval::removed_stop_codon`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Given {
    Lost,
    Retained,
}

/// How many positions a deletion is followed along the reference sequence in each
/// direction for the other places the same edit can be written at, as vep-rs follows it.
const EQUIVALENT_SITES_LIMIT: i64 = 10_000;

/// The names of the twenty amino acids and of the stop in a protein notation.
fn three_letters(residue: u8) -> Option<&'static str> {
    Some(match residue {
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
        b'*' => "Ter",
        _ => return None,
    })
}

/// A row read on the completed coding sequence: its terms and the cells that vep-rs
/// writes otherwise than Ensembl VEP there.
struct OnTheCompletedSequence {
    terms: BTreeSet<CodingTerm>,
    amino_acids: String,
    codons: String,
    notation: Option<String>,
}

/// The terms and cells of a sequence allele's row on a transcript whose record shows
/// that polyadenylation completes its stop codon. vep-rs reads every term and cell of
/// such a row on the transcript with the completed codon in its coding sequence and
/// the adenines of the poly(A) tail after it (`completed_by_polyadenylation`), where
/// the stop codon is whole: the cached coding sequence followed by the adenines that
/// complete its last codon, a 3' UTR of adenines as long as that sequence, and the
/// coding end that many cDNA positions later. Ensembl VEP reads the record, on which
/// the last codon is incomplete (VariationEffect.pm:1588-1614). The predicates are
/// evaluated here as `PerlCodingEval` evaluates them on the completed transcript.
fn on_the_completed_sequence(
    row: &RowFacts<'_>,
    ends: &CodingEnds,
    (theirs, ours): &(Cells<'_>, Cells<'_>),
    lengths: Lengths,
    allele: &[u8],
    reference_bases: Option<&dyn Fn(i64) -> Option<u8>>,
) -> Result<OnTheCompletedSequence, NotDerived> {
    let kept = ends
        .completed
        .as_deref()
        .ok_or(NotDerived::SequenceNotKept)?;
    let model = ends.model.as_deref().ok_or(NotDerived::ModelNotKept)?;
    let code = super::CODES
        .get(&kept.table)
        .ok_or(NotDerived::SequenceNotKept)?;
    // The completing adenines have no genomic base and no allele lies on them, so the
    // three position cells are the record's on both rows.
    equal_cells(theirs, ours, true)?;

    let cached_length = kept.coding.len();
    let (whole_codons, tail) = (cached_length / 3, cached_length % 3);
    let mut coding = kept.coding.clone();
    coding.resize(cached_length + 3 - tail, b'A');
    let length = coding.len();
    let adenines = vec![b'A'; length];
    let coding_end = i128::from(model.coding_end) + 3 - tail as i128;
    let reverse = ends.strand == Some(-1);
    // The allele as the transcript reads it. The `Allele` cell of a record given on the
    // forward strand is on the genome's forward strand, and vep-rs reads it reversed
    // and complemented on a reverse-strand transcript (`feature_seq`).
    let allele: Vec<u8> = match reverse {
        false => allele.to_vec(),
        true => allele
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect(),
    };
    let neighbours = between_neighbours(ends, row.location, lengths)?;
    let cdna = positions(theirs.cdna_position, "cDNA_position", neighbours)?;
    let cds = positions(theirs.cds_position, "CDS_position", neighbours)?;
    let protein = positions(theirs.protein_position, "Protein_position", neighbours)?;
    let (lowest, highest) = row.location;

    // One shape is left out: an insertion of a multiple of three bases between the two
    // tail bases. A deletion that begins at the second of two tail bases and runs past
    // the transcript's last base is read as the sequence given (`given`): it leaves the
    // first tail base before the adenines.
    let at_the_second_tail_base = tail == 2 && cds.0 == Some(cached_length as i64);
    if at_the_second_tail_base && lengths.between_two_bases() && lengths.allele % 3 == 0 {
        return Err(NotDerived::ShapeNotCovered);
    }

    // `TranscriptVariationAllele::codon`: the reference codon is the codons the span's
    // translation positions name, the alternate the same window of the coding sequence
    // with the allele written over the CDS span and the 3' UTR after it, longer or
    // shorter by the allele's length less the span's. A coding sequence the allele
    // leaves with fewer than three bases is read as empty. An empty window is `-`.
    let window = match (truthy(protein.0), truthy(protein.1), cds) {
        (Some(first), Some(last), (Some(start), Some(end))) => {
            let codon_start = i128::from(first) * 3 - 2;
            let codon_length = i128::from(last) * 3 - codon_start + 1;
            Some((codon_start, codon_length, start, end))
        }
        _ => None,
    };
    let alternate_sequence = window.map(|(_, _, start, end)| {
        let clip = |position: i64| usize::try_from(position).unwrap_or(0).min(length);
        let mut alternate = coding[..clip(start - 1)].to_vec();
        alternate.extend_from_slice(&allele);
        alternate.extend_from_slice(&coding[clip(end)..]);
        if alternate.len() < 3 {
            alternate.clear();
        }
        alternate.extend_from_slice(&adenines);
        alternate
    });
    let or_dash = |codon: Option<&[u8]>| match codon {
        Some(codon) if !codon.is_empty() => codon.to_vec(),
        _ => b"-".to_vec(),
    };
    let reference_codon = window.map(|(codon_start, codon_length, _, _)| {
        or_dash(perl_substr(
            &coding,
            nearest(codon_start - 1),
            nearest(codon_length),
        ))
    });
    let alternate_codon = window.zip(alternate_sequence.as_deref()).map(
        |((codon_start, codon_length, start, end), alternate)| {
            let spanned = i128::from(end) - i128::from(start) + 1;
            or_dash(perl_substr(
                alternate,
                nearest(codon_start - 1),
                nearest(codon_length + i128::from(lengths.allele) - spanned),
            ))
        },
    );

    // `TranscriptVariationAllele::peptide`: the whole codons translated, the unknown
    // residue for a partial last codon unless the peptide is the lone stop, and `-` for
    // nothing translated. The reference peptide takes the residues an edit of the
    // sequence writes, which are read from the cached peptide: a window that holds a
    // residue the cached peptide states otherwise than its codon is not derived, nor
    // one that holds protein position 1, unless the first codon is whole, shows no
    // start codon and carries no `cds_start_NF`: vep-rs reads that codon as any other.
    let translate = |codon: &[u8]| -> Vec<u8> {
        if codon == b"-" {
            return b"-".to_vec();
        }
        let whole = codon.len() / 3 * 3;
        let mut peptide: Vec<u8> = codon[..whole]
            .chunks(3)
            .map(|codon| code.residue(codon))
            .collect();
        if whole < codon.len() && peptide != b"*" {
            peptide.push(b'X');
        }
        if peptide.is_empty() {
            peptide.push(b'-');
        }
        peptide
    };
    let alternate_peptide = alternate_codon.as_deref().map(translate);
    let reference_peptide = reference_codon.as_deref().map(translate);
    if let (Some(codon), (Some(first), Some(last))) = (reference_codon.as_deref(), protein) {
        if codon != b"-" {
            let (first, last) = (first.min(last), first.max(last));
            let as_any_other = StartEnd::NotShown {
                flagged: false,
                padded: false,
            };
            if first < 1 || (first == 1 && ends.start != as_any_other) {
                return Err(NotDerived::SequenceEdit);
            }
            for position in first..=last.min(whole_codons as i64) {
                let position = position as usize;
                let codon = &kept.coding[3 * position - 3..3 * position];
                if kept
                    .peptide
                    .get(position - 1)
                    .is_some_and(|cached| *cached != code.residue(codon))
                {
                    return Err(NotDerived::SequenceEdit);
                }
            }
        }
    }
    let strip = |text: &Vec<u8>| match text.as_slice() {
        b"-" => Vec::new(),
        text => text.to_vec(),
    };
    let peptides: Option<Alleles> = reference_peptide
        .as_ref()
        .zip(alternate_peptide.as_ref())
        .map(|(reference, alternate)| (strip(reference), strip(alternate)));
    let codons: Option<Alleles> = reference_codon
        .as_ref()
        .zip(alternate_codon.as_ref())
        .map(|(reference, alternate)| (strip(reference), strip(alternate)));

    // The stop codon edit test (VariationEffect.pm:1521-1529): the allele written into
    // the coding sequence and 3' UTR at the CDS start, over the cDNA length of the span.
    let edited = match (
        lengths.insertion() || lengths.deletion(),
        truthy(cdna.0),
        truthy(cdna.1),
        truthy(cds.0),
    ) {
        (true, Some(cdna_start), Some(cdna_end), Some(cds_start)) => {
            let mut edited = coding.clone();
            edited.extend_from_slice(&adenines);
            perl_substr_assign(
                &mut edited,
                cds_start - 1,
                cdna_end - cdna_start + 1,
                &allele,
            );
            if edited.len() < length {
                Edited::Shorter
            } else {
                let stop = length - 3;
                let first = (usize::try_from(cds_start - 1).unwrap_or(0) / 3 * 3).min(stop);
                let residues_kept = (first..stop).step_by(3).all(|codon| {
                    let (reference, read) = (&coding[codon..codon + 3], &edited[codon..codon + 3]);
                    reference == read || {
                        let residue = code.residue(read);
                        residue != b'X' && residue == code.residue(reference)
                    }
                });
                Edited::Read {
                    codon: [edited[stop], edited[stop + 1], edited[stop + 2]],
                    residues_kept,
                }
            }
        }
        _ => Edited::NotMade,
    };
    // `_overlaps_stop_codon` (VariationEffect.pm:1363-1385) on the whole stop codon of
    // the completed sequence, and `_ins_del_stop_altered` behind it (1433-1484).
    let overlaps_the_stop_codon = truthy(cdna.0)
        .zip(truthy(cdna.1))
        .is_some_and(|(start, end)| {
            i128::from(end) >= coding_end - 2 && i128::from(start) <= coding_end
        });
    let altered = overlaps_the_stop_codon
        && match edited {
            Edited::NotMade => false,
            Edited::Shorter => true,
            Edited::Read { codon, .. } => code.residue(&codon) != b'*',
        };
    // The genomic bases of the stop codon are the tail's: the completing adenines have
    // none.
    let tail_bases: Vec<i64> = ends.tail.iter().flatten().copied().collect();
    // A span that lacks a coordinate of the edit test and holds every base the genome
    // has of the stop codon from its first base on is read on the sequence the record
    // gives (`PerlCodingEval::given_stop_codon`).
    let coordinates =
        truthy(cdna.0).is_some() && truthy(cdna.1).is_some() && truthy(cds.0).is_some();
    let holds_the_stop_codon = !lengths.between_two_bases()
        && tail_bases
            .iter()
            .any(|base| lowest <= *base && *base <= highest)
        && tail_bases.iter().all(|base| match reverse {
            false => *base <= highest,
            true => *base >= lowest,
        });
    let given =
        match (lengths.insertion() || lengths.deletion()) && !coordinates && holds_the_stop_codon {
            true => deletion_past_the_transcript_end(
                DeletionPastTheEnd {
                    location: row.location,
                    span: ends.span,
                    reverse,
                    model,
                    coding: &coding,
                    cached_length,
                    cds_start: truthy(cds.0),
                    shape: lengths.deletion()
                        && allele.is_empty()
                        && truthy(cdna.0).is_some()
                        && cdna.1.is_none(),
                },
                code,
                reference_bases,
            )?,
            false => None,
        };

    // `stop_lost` (VariationEffect.pm:1234-1288) and `stop_retained` (1290-1324) as
    // vep-rs evaluates them where the stop codon stands (`PerlCodingEval::stop_lost`,
    // `stop_retained`).
    let stop_lost = match &peptides {
        Some((reference, alternate)) if !alternate.contains(&b'X') => {
            !alternate.contains(&b'*') && reference.contains(&b'*')
        }
        _ => altered || given == Some(Given::Lost),
    };
    let stop_retained = !stop_lost
        && match &peptides {
            Some((reference, alternate)) if !alternate.is_empty() && !alternate.contains(&b'X') => {
                reference_equals_alternate_sequence(reference, alternate, protein, &kept.peptide)
            }
            _ => {
                // `_overlaps_stop_codon_cil` (VariationEffect.pm:1389-1431) as vep-rs
                // reads it: an insertion lies at or before the stop codon and leaves
                // every residue before it; any other span holds a genomic base of the
                // stop codon. Behind it, the edit test finds a stop codon at the old
                // stop position after the reference's residues.
                let overlaps = match lengths.between_two_bases() {
                    true => truthy(cds.0).is_some_and(|follows| {
                        1 <= follows
                            && follows <= length as i64
                            && match allele.is_empty() {
                                true => follows >= length as i64 - 2,
                                false => matches!(
                                    edited,
                                    Edited::Read {
                                        residues_kept: true,
                                        ..
                                    }
                                ),
                            }
                    }),
                    false => {
                        let mut bases = tail_bases.clone();
                        bases.sort_unstable();
                        bases
                            .chunk_by(|one, next| i128::from(*next) - i128::from(*one) == 1)
                            .any(|run| overlap(lowest, highest, run[0], run[run.len() - 1]))
                    }
                };
                let kept_by_the_edit = matches!(
                    edited,
                    Edited::Read { codon, residues_kept: true }
                        if codon.iter().all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
                            && code.residue(&codon) == b'*'
                );
                ((lengths.insertion() || lengths.deletion()) && overlaps && kept_by_the_edit)
                    || given == Some(Given::Retained)
            }
        };
    let theirs_terms = coding_family(row.ensembl_vep_terms);
    let (start_lost, start_retained) = start_terms(ends.start, &theirs_terms)?;
    let terms = terms_of(
        &Values {
            partial_codon: false,
            lengths,
            cds,
            peptides: peptides.as_ref(),
            codons: codons.as_ref(),
            alternate_peptide_undefined: alternate_peptide.is_none(),
        },
        Charged {
            stop_lost,
            stop_retained,
            start_lost,
            start_retained,
        },
        || within_cds(theirs.cds_position, ends, row.location, lengths),
    )?;
    at_a_standing_start_codon(row, ends, &theirs_terms, &terms)?;

    // The cells (`pep_allele_string`, `display_codon_allele_string`): the reference
    // peptide, and the alternate after a `/` where it differs; the two codons in lower
    // case with each allele's own bases in upper case from the span's place in its
    // codon.
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    let amino_acids = match (&reference_peptide, &alternate_peptide) {
        (Some(reference), Some(alternate)) if reference == alternate => text(reference),
        (Some(reference), Some(alternate)) => format!("{}/{}", text(reference), text(alternate)),
        _ => "-".to_string(),
    };
    let place = cdna.0.map(|start| {
        let before = i128::from(start) - i128::from(model.coding_start);
        (before + i128::from(model.padding)).rem_euclid(3)
    });
    let display = |codon: &[u8], own_bases: i64| {
        let mut display = codon.to_ascii_lowercase();
        if let (Some(place), true, true) = (place, own_bases > 0, codon != b"-") {
            let from = (place as usize).min(display.len());
            let to = from.saturating_add(own_bases as usize).min(display.len());
            display[from..to].make_ascii_uppercase();
        }
        text(&display)
    };
    let codons_cell = match (&reference_codon, &alternate_codon) {
        (Some(reference), Some(alternate)) => format!(
            "{}/{}",
            display(reference, lengths.reference),
            display(alternate, lengths.allele)
        ),
        _ => "-".to_string(),
    };

    let notation = protein_notation(
        &ProteinChange {
            peptide: &kept.peptide,
            whole_codons,
            alternate_sequence: alternate_sequence.as_deref(),
            first_codon: protein.0.zip(protein.1).map(|(one, other)| one.min(other)),
            substitution: lengths.same(),
            frameshift: terms.contains(&CodingTerm::FrameshiftVariant),
            stop_lost: terms.contains(&CodingTerm::StopLost),
        },
        code,
    );
    Ok(OnTheCompletedSequence {
        terms,
        amino_acids,
        codons: codons_cell,
        notation,
    })
}

/// `ref_eq_alt_sequence` (VariationEffect.pm:1326-1361) where a stop codon stands, on
/// the transcript's peptide: the span starts beyond the peptide and the alternate
/// peptide opens with a stop (1344); the peptide is unchanged up to its old length and
/// a stop follows it (1357); or both peptide alleles hold their first stop at one place
/// (1358).
fn reference_equals_alternate_sequence(
    reference: &[u8],
    alternate: &[u8],
    protein: (Option<i64>, Option<i64>),
    peptide: &[u8],
) -> bool {
    let (Some(start), Some(end)) = protein else {
        return false;
    };
    if reference == b"X" && alternate == b"X" {
        return false;
    }
    if start > peptide.len() as i64 && alternate.first() == Some(&b'*') {
        return true;
    }
    let mut changed = peptide.to_vec();
    perl_substr_assign(&mut changed, start - 1, end - start + 1, alternate);
    if changed.len() > peptide.len()
        && changed[..peptide.len()] == *peptide
        && changed[peptide.len()] == b'*'
    {
        return true;
    }
    let stop = |residues: &[u8]| residues.iter().position(|residue| *residue == b'*');
    stop(reference).is_some() && stop(alternate) == stop(reference)
}

/// A span on a completed stop codon that lacks a coordinate of the edit test.
struct DeletionPastTheEnd<'a> {
    location: (i64, i64),
    span: Option<(i64, i64)>,
    reverse: bool,
    model: &'a CodingModel,
    /// The completed coding sequence.
    coding: &'a [u8],
    /// The length of the cached coding sequence, without the completing adenines.
    cached_length: usize,
    cds_start: Option<i64>,
    /// Whether the span is the one shape derived: a deletion with no alternate base,
    /// a cDNA and a CDS start, and no cDNA end.
    shape: bool,
}

/// What the sequence the record gives shows at the old stop position
/// (`PerlCodingEval::given_stop_codon`), for the one shape derived: a deletion that
/// begins in the coding sequence and runs past the transcript's last base. The
/// reference bases tell how many places the same deletion can be written further along
/// the transcript and how many back; a comparison without a reference sequence reads
/// the deletion at the record's own place, as vep-rs does when it is run without one.
/// A class takes a pair only where vep-rs's row holds the terms derived
/// (`CodingEndClass::takes`), so a row on which vep-rs reads other terms at another
/// place is not set aside. Written where it takes no coding base, it keeps
/// the stop codon. Written with its end inside the transcript, the bases after it are
/// read, which is not derived. Otherwise the transcript's sequence ends where the
/// deletion begins at its furthest place, and the adenines of the poly(A) tail follow:
/// the stop codon is kept where that sequence is the reference's up to the coding end
/// or holds a stop codon at the old stop position after the reference's residues, lost
/// where it holds none there. `None` where the mapper pairs are not the transcript's
/// exons in cDNA order, on which vep-rs reads no given sequence.
fn deletion_past_the_transcript_end(
    asked: DeletionPastTheEnd<'_>,
    code: &super::GeneticCode,
    reference_bases: Option<&dyn Fn(i64) -> Option<u8>>,
) -> Result<Option<Given>, NotDerived> {
    let (lowest, highest) = asked.location;
    let wide = i128::from;
    let past_the_end = match (asked.span, asked.reverse) {
        (Some((_, end)), false) => wide(highest) - wide(end),
        (Some((start, _)), true) => wide(start) - wide(lowest),
        (None, _) => 0,
    };
    let (true, Some(cds_start), true) = (asked.shape, asked.cds_start, past_the_end >= 1) else {
        return Err(NotDerived::SpanLacksACoordinate);
    };
    // The mapper pairs as exons on an axis that runs in the transcript's direction.
    let strand = if asked.reverse { -1 } else { 1 };
    let mut pairs = asked.model.pairs.clone();
    pairs.sort_unstable();
    let (mut next_cdna, mut before) = (1, None);
    for &(from_start, from_end, to_start, to_end, orientation) in &pairs {
        let (first, last) = match asked.reverse {
            false => (wide(to_start), wide(to_end)),
            true => (-wide(to_end), -wide(to_start)),
        };
        if orientation != strand
            || wide(from_start) != next_cdna
            || wide(from_end) - wide(from_start) != wide(to_end) - wide(to_start)
            || before.is_some_and(|before| before >= first)
        {
            return Ok(None);
        }
        (next_cdna, before) = (wide(from_end) + 1, Some(last));
    }
    let no_base = |_: i64| None;
    let base = reference_bases.unwrap_or(&no_base);
    let (start, end) = match asked.reverse {
        false => (wide(lowest), wide(highest)),
        true => (-wide(highest), -wide(lowest)),
    };
    let at = |position: i128| {
        super::narrow(if asked.reverse { -position } else { position }).and_then(base)
    };
    let repeats = |from: i128, onto: i128, step: i128| {
        (0..EQUIVALENT_SITES_LIMIT)
            .take_while(|place| {
                let moved = step * wide(*place);
                matches!(
                    (at(from + moved), at(onto + moved)),
                    (Some(one), Some(other))
                        if one == other
                            && matches!(one.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
                )
            })
            .count() as i64
    };
    let (along, back) = (repeats(start, end + 1, 1), repeats(end, start - 1, -1));
    // Written `along` places further on, the deletion takes that many coding bases
    // fewer: none where it can be written past every coding base the record holds.
    let coding_removed = asked.cached_length as i64 - cds_start + 1;
    if along >= coding_removed {
        return Ok(Some(Given::Retained));
    }
    if wide(back) >= past_the_end {
        return Err(NotDerived::PlacementInsideTheTranscript);
    }
    // The sequence given: the coding bases before the deletion at its furthest place,
    // then the adenines of the poly(A) tail.
    let length = asked.coding.len();
    let kept = usize::try_from(cds_start - 1 + along)
        .unwrap_or(0)
        .min(length);
    let mut given = asked.coding[..kept].to_vec();
    given.resize(length, b'A');
    if given == asked.coding {
        return Ok(Some(Given::Retained));
    }
    if code.residue(&given[length - 3..]) != b'*' {
        return Ok(Some(Given::Lost));
    }
    // A stop codon at the old stop position keeps the stop codon where every codon
    // before it reads the reference's residue.
    let before = length - 3;
    let lead = usize::try_from(asked.model.padding)
        .unwrap_or(0)
        .min(before);
    let whole = (lead + (3 - lead % 3) % 3).min(before);
    let same = given[..whole] == asked.coding[..whole]
        && given[whole..before]
            .chunks(3)
            .zip(asked.coding[whole..before].chunks(3))
            .all(|(read, reference)| {
                read == reference || {
                    let residue = code.residue(read);
                    residue != b'X' && residue == code.residue(reference)
                }
            });
    match same {
        true => Ok(Some(Given::Retained)),
        false => Err(NotDerived::ShapeNotCovered),
    }
}
