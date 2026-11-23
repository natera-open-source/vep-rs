// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The situations a record exercises, one function a category. A function reads one
//! transcript, or one place of a sequence, and gives the change that puts a record in its
//! situation, drawing every choice from the stream it is handed; `None` says the
//! transcript holds no such place.

use vep_core::codon::{complement_base, is_start_codon, known_table, translate_codon_with_table};

use crate::cache::{Coding, Span, Strand, Transcript};
use crate::random::SplitMix64;

/// How far from a transcript's end a record counts as upstream or downstream of it, in
/// bases: Ensembl VEP's default distance.
pub const FLANK: u64 = 5000;

/// How many bases a change inside an exon of a multi-exon transcript keeps from each end
/// of that exon, which puts it outside the splice region.
pub const EXON_MARGIN: u64 = 6;

/// How many bases an intron record keeps from each end of its intron.
pub const INTRON_MARGIN: u64 = 20;

/// The shortest intron a splice site record is placed in.
pub const SPLICE_INTRON: u64 = 20;

const BASES: [u8; 4] = *b"ACGT";

/// A situation of the genome check. A record's ID is its category's name and an ordinal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Category {
    /// A codon the named genetic code reads differently from the standard code.
    ChangedCodon(u8),
    /// The first codon of a coding sequence on the named code.
    StartCodon(u8),
    /// A codon turned into a stop codon of the named code.
    StopGained(u8),
    /// A stop codon of the named code turned into a sense codon.
    StopLost(u8),
    /// A codon the named code reads differently, made by a substitution that keeps the
    /// residue.
    MadeCodonSynonymous(u8),
    /// A codon the named code reads differently, made by a substitution that changes
    /// the residue.
    MadeCodonMissense(u8),
    /// A codon the named code reads differently, made by a substitution in the closing
    /// stop codon.
    MadeCodonStop(u8),
    FirstSequence,
    LastSequence,
    Organelle,
    Synonym,
    Missense(Strand),
    Synonymous(Strand),
    SpliceDonor,
    SpliceAcceptor,
    FivePrimeUtr,
    ThreePrimeUtr,
    Intron,
    Upstream,
    Downstream,
    Intergenic,
    FrameshiftInsertion,
    FrameshiftDeletion,
    InframeDeletion,
    FirstBase,
    LastBase,
    Origin,
}

impl Category {
    /// The categories of one genetic code, in the summary's order: four under every
    /// code, then the code's made-codon categories ([`made_codon_categories`]).
    pub fn of_code(code: u8) -> Vec<Category> {
        let mut categories = vec![
            Category::ChangedCodon(code),
            Category::StartCodon(code),
            Category::StopGained(code),
            Category::StopLost(code),
        ];
        categories.extend(made_codon_categories(code));
        categories
    }

    /// The categories that name no genetic code, in the summary's order.
    pub const FIXED: [Category; 22] = [
        Category::FirstSequence,
        Category::LastSequence,
        Category::Organelle,
        Category::Synonym,
        Category::Missense(Strand::Forward),
        Category::Missense(Strand::Reverse),
        Category::Synonymous(Strand::Forward),
        Category::Synonymous(Strand::Reverse),
        Category::SpliceDonor,
        Category::SpliceAcceptor,
        Category::FivePrimeUtr,
        Category::ThreePrimeUtr,
        Category::Intron,
        Category::Upstream,
        Category::Downstream,
        Category::Intergenic,
        Category::FrameshiftInsertion,
        Category::FrameshiftDeletion,
        Category::InframeDeletion,
        Category::FirstBase,
        Category::LastBase,
        Category::Origin,
    ];

    /// The category's name: the prefix of its records' IDs and its line of the summary.
    pub fn name(&self) -> String {
        let strand = |strand: &Strand| match strand {
            Strand::Forward => "forward",
            Strand::Reverse => "reverse",
        };
        match self {
            Category::ChangedCodon(code) => format!("code{code}_changed_codon"),
            Category::StartCodon(code) => format!("code{code}_start_codon"),
            Category::StopGained(code) => format!("code{code}_stop_gained"),
            Category::StopLost(code) => format!("code{code}_stop_lost"),
            Category::MadeCodonSynonymous(code) => format!("code{code}_made_codon_synonymous"),
            Category::MadeCodonMissense(code) => format!("code{code}_made_codon_missense"),
            Category::MadeCodonStop(code) => format!("code{code}_made_codon_stop"),
            Category::FirstSequence => "first_sequence".to_owned(),
            Category::LastSequence => "last_sequence".to_owned(),
            Category::Organelle => "organelle".to_owned(),
            Category::Synonym => "synonym".to_owned(),
            Category::Missense(s) => format!("missense_{}", strand(s)),
            Category::Synonymous(s) => format!("synonymous_{}", strand(s)),
            Category::SpliceDonor => "splice_donor".to_owned(),
            Category::SpliceAcceptor => "splice_acceptor".to_owned(),
            Category::FivePrimeUtr => "five_prime_utr".to_owned(),
            Category::ThreePrimeUtr => "three_prime_utr".to_owned(),
            Category::Intron => "intron".to_owned(),
            Category::Upstream => "upstream".to_owned(),
            Category::Downstream => "downstream".to_owned(),
            Category::Intergenic => "intergenic".to_owned(),
            Category::FrameshiftInsertion => "frameshift_insertion".to_owned(),
            Category::FrameshiftDeletion => "frameshift_deletion".to_owned(),
            Category::InframeDeletion => "inframe_deletion".to_owned(),
            Category::FirstBase => "first_base".to_owned(),
            Category::LastBase => "last_base".to_owned(),
            Category::Origin => "origin".to_owned(),
        }
    }

    /// The category's place in the summary and among records of one position: the
    /// categories of a code by code, then the others in the order of [`Category::FIXED`].
    pub fn rank(&self) -> (u8, u16, u8) {
        match self {
            Category::ChangedCodon(code) => (0, u16::from(*code), 0),
            Category::StartCodon(code) => (0, u16::from(*code), 1),
            Category::StopGained(code) => (0, u16::from(*code), 2),
            Category::StopLost(code) => (0, u16::from(*code), 3),
            Category::MadeCodonSynonymous(code) => (0, u16::from(*code), 4),
            Category::MadeCodonMissense(code) => (0, u16::from(*code), 5),
            Category::MadeCodonStop(code) => (0, u16::from(*code), 6),
            other => {
                let place = Category::FIXED
                    .iter()
                    .position(|category| category == other)
                    .unwrap_or(Category::FIXED.len());
                (1, place as u16, 0)
            }
        }
    }
}

/// How a record's alternate allele follows from its reference allele.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Alt {
    /// These bases, on the forward strand.
    Bases(Vec<u8>),
    /// One of the three bases the one-base reference allele is not, by its number from 0.
    Other(u8),
    /// The reference allele's first base alone: the rest is deleted.
    First,
    /// The reference allele with these bases after it: an insertion.
    Append(Vec<u8>),
}

/// A record before its alleles are settled: where it stands, the reference bases the
/// cache gives for it and how its alternate allele is made.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Change {
    /// The position of the reference allele's first base.
    pub position: u64,
    /// The reference allele on the forward strand, base by base; `None` where the
    /// transcript the change was read from gives no base.
    pub reference: Vec<Option<u8>>,
    pub alt: Alt,
    /// The length of the sequence where the reference allele runs past its last base and
    /// on from its first.
    pub wraps_at: Option<u64>,
}

impl Change {
    fn substitution(position: u64, reference: Option<u8>, alt: Alt) -> Change {
        Change {
            position,
            reference: vec![reference],
            alt,
            wraps_at: None,
        }
    }

    /// The sequence position of each base of the reference allele.
    pub fn positions(&self) -> Vec<u64> {
        (0..self.reference.len() as u64)
            .map(|offset| match self.wraps_at {
                Some(length) if self.position > length => self.position - length + offset,
                Some(length) if offset > length - self.position => {
                    offset - (length - self.position)
                }
                _ => self.position + offset,
            })
            .collect()
    }

    /// The alternate allele that goes with `reference`, the settled reference allele.
    pub fn alt_for(&self, reference: &[u8]) -> Vec<u8> {
        match &self.alt {
            Alt::Bases(bases) => bases.clone(),
            Alt::Other(number) => vec![other_base(reference[0], *number)],
            Alt::First => reference[..1].to_vec(),
            Alt::Append(bases) => {
                let mut alt = reference.to_vec();
                alt.extend_from_slice(bases);
                alt
            }
        }
    }
}

/// The base numbered `number` among the three that `base` is not, in the order A, C, G, T.
pub fn other_base(base: u8, number: u8) -> u8 {
    let others: Vec<u8> = BASES.iter().copied().filter(|b| *b != base).collect();
    others[usize::from(number) % others.len()]
}

/// A base drawn from the stream.
pub fn drawn_base(stream: &mut SplitMix64) -> u8 {
    BASES[stream.below(4) as usize]
}

fn forward(strand: Strand, base: u8) -> u8 {
    match strand {
        Strand::Forward => base,
        Strand::Reverse => complement_base(base),
    }
}

/// A substitution at `position` whose reference base the cache does not give here and
/// whose alternate base is drawn among the other three.
fn unread_substitution(position: u64, stream: &mut SplitMix64) -> Change {
    Change::substitution(position, None, Alt::Other(stream.below(3) as u8))
}

/// The position numbered `number` among the positions of `spans`, counted along them.
fn nth_position(spans: &[Span], mut number: u64) -> Option<u64> {
    for span in spans {
        if number < span.length() {
            return Some(span.start + number);
        }
        number -= span.length();
    }
    None
}

/// A position drawn among the positions of `spans`, every position equally likely.
fn draw_position(spans: &[Span], stream: &mut SplitMix64) -> Option<u64> {
    let total: u64 = spans.iter().map(Span::length).sum();
    if total == 0 {
        return None;
    }
    nth_position(spans, stream.below(total))
}

/// The part of `span` that lies in `within`.
fn overlap(span: Span, within: Span) -> Option<Span> {
    let start = span.start.max(within.start);
    let end = span.end.min(within.end);
    (start <= end).then_some(Span { start, end })
}

/// A transcript whose codons the cache gives ([`Transcript::reads_codons`]) on a genetic
/// code the annotator has a table for, with what
/// the codon categories read from it.
struct Reading<'a> {
    transcript: &'a Transcript,
    coding: &'a Coding,
    /// Each exon's span on the spliced transcript, cut back by [`EXON_MARGIN`] at both
    /// ends in a multi-exon transcript: where a change is clear of its exon's ends.
    clear: Vec<Span>,
    /// The exons' whole spans on the spliced transcript.
    exons: Vec<Span>,
    /// How many whole codons the coding sequence holds.
    codons: u64,
}

impl<'a> Reading<'a> {
    fn new(transcript: &'a Transcript) -> Option<Reading<'a>> {
        if !transcript.reads_codons() || !known_table(transcript.code()) {
            return None;
        }
        let coding = transcript.coding.as_ref()?;
        let exons = transcript.exon_cdna_spans();
        let margin = if exons.len() > 1 { EXON_MARGIN } else { 0 };
        let clear = exons
            .iter()
            .filter(|exon| exon.length() > 2 * margin)
            .map(|exon| Span {
                start: exon.start + margin,
                end: exon.end - margin,
            })
            .collect();
        Some(Reading {
            transcript,
            coding,
            clear,
            exons,
            codons: coding.sequence.len() as u64 / 3,
        })
    }

    fn code(&self) -> u8 {
        self.transcript.code()
    }

    fn codon(&self, number: u64) -> [u8; 3] {
        let at = (number * 3) as usize;
        [
            self.coding.sequence[at],
            self.coding.sequence[at + 1],
            self.coding.sequence[at + 2],
        ]
    }

    /// The codon's span on the spliced transcript.
    fn span(&self, number: u64) -> Span {
        let start = self.coding.cdna_start + number * 3;
        Span {
            start,
            end: start + 2,
        }
    }

    fn lies_in(&self, number: u64, spans: &[Span]) -> bool {
        let codon = self.span(number);
        spans
            .iter()
            .any(|span| span.start <= codon.start && codon.end <= span.end)
    }

    /// The substitution that writes `base`, on the transcript's strand, over base
    /// `offset` of codon `number`.
    fn substitution(&self, number: u64, offset: usize, base: u8) -> Option<Change> {
        let strand = self.transcript.strand;
        let position = self
            .transcript
            .genomic_position(self.span(number).start + offset as u64)?;
        Some(Change::substitution(
            position,
            Some(forward(strand, self.codon(number)[offset])),
            Alt::Bases(vec![forward(strand, base)]),
        ))
    }

    /// A substitution in one codon from `first` to `last` that `admits` accepts, given
    /// the codon before and after: the scan starts at a drawn codon and takes the first
    /// that has an accepted substitution, drawing among that codon's.
    fn pick(
        &self,
        first: u64,
        last: u64,
        spans: &[Span],
        stream: &mut SplitMix64,
        admits: impl Fn(&[u8; 3], &[u8; 3]) -> bool,
    ) -> Option<Change> {
        if last < first || last >= self.codons {
            return None;
        }
        let count = last - first + 1;
        let begin = stream.below(count);
        for step in 0..count {
            let number = first + (begin + step) % count;
            if !self.lies_in(number, spans) {
                continue;
            }
            let before = self.codon(number);
            let mut accepted: Vec<(usize, u8)> = Vec::new();
            for offset in 0..3 {
                for base in BASES {
                    if base == before[offset] {
                        continue;
                    }
                    let mut after = before;
                    after[offset] = base;
                    if admits(&before, &after) {
                        accepted.push((offset, base));
                    }
                }
            }
            if !accepted.is_empty() {
                let (offset, base) = accepted[stream.below(accepted.len() as u64) as usize];
                return self.substitution(number, offset, base);
            }
        }
        None
    }

    /// A substitution in a codon that is neither the first nor the last and lies clear
    /// of its exon's ends.
    fn pick_inside(
        &self,
        stream: &mut SplitMix64,
        admits: impl Fn(&[u8; 3], &[u8; 3]) -> bool,
    ) -> Option<Change> {
        if self.codons < 3 {
            return None;
        }
        self.pick(1, self.codons - 2, &self.clear, stream, admits)
    }

    fn residue(&self, codon: &[u8; 3]) -> u8 {
        translate_codon_with_table(codon, self.code())
    }

    /// Where on the sequence an indel can be anchored: the runs of anchor positions from
    /// which `touched` bases, the anchor included, stay inside the coding part of one
    /// exon, clear of its ends and of the first and the last codon.
    fn anchors(&self, touched: u64) -> Vec<Span> {
        let length = self.coding.sequence.len() as u64;
        if length < 12 {
            return Vec::new();
        }
        let inner = Span {
            start: self.coding.cdna_start + 3,
            end: self.coding.cdna_start + length - 4,
        };
        let mut anchors = Vec::new();
        for clear in &self.clear {
            let Some(part) = overlap(*clear, inner) else {
                continue;
            };
            let (Some(a), Some(b)) = (
                self.transcript.genomic_position(part.start),
                self.transcript.genomic_position(part.end),
            ) else {
                continue;
            };
            let (low, high) = (a.min(b), a.max(b));
            if high - low + 1 >= touched {
                anchors.push(Span {
                    start: low,
                    end: high + 1 - touched,
                });
            }
        }
        anchors.sort();
        anchors
    }

    /// The forward-strand bases of `count` positions from `position` on.
    fn bases(&self, position: u64, count: u64) -> Option<Vec<Option<u8>>> {
        (position..position + count)
            .map(|at| self.transcript.forward_base(at).map(Some))
            .collect()
    }

    /// The number, from 0, of the coding sequence's base at sequence position `genomic`.
    fn coding_index(&self, genomic: u64) -> Option<u64> {
        self.transcript
            .cdna_position(genomic)?
            .checked_sub(self.coding.cdna_start)
    }
}

/// How a genetic code differs from the standard code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodeDifference {
    /// At least one codon translates to another residue.
    Residues,
    /// Every codon translates alike and at least one codon's start mark differs.
    StartMarks,
    /// Nothing differs.
    None,
}

fn every_codon() -> impl Iterator<Item = [u8; 3]> {
    BASES.into_iter().flat_map(|first| {
        BASES
            .into_iter()
            .flat_map(move |second| BASES.into_iter().map(move |third| [first, second, third]))
    })
}

fn reads_differently(codon: &[u8; 3], code: u8) -> bool {
    translate_codon_with_table(codon, code) != translate_codon_with_table(codon, 1)
}

fn starts_differently(codon: &[u8; 3], code: u8) -> bool {
    is_start_codon(codon, code) != is_start_codon(codon, 1)
}

/// How `code` differs from the standard code, over the 64 codons. A code the annotator
/// has no table for differs in nothing that can be read.
pub fn code_difference(code: u8) -> CodeDifference {
    if !known_table(code) {
        CodeDifference::None
    } else if every_codon().any(|codon| reads_differently(&codon, code)) {
        CodeDifference::Residues
    } else if every_codon().any(|codon| starts_differently(&codon, code)) {
        CodeDifference::StartMarks
    } else {
        CodeDifference::None
    }
}

/// A substitution inside a codon the transcript's code reads differently from the
/// standard code. Under a code that differs in start marks alone it is a substitution in
/// the first codon, where that codon's start mark differs before or after the change.
pub fn changed_codon(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    let code = reading.code();
    match code_difference(code) {
        CodeDifference::Residues => {
            reading.pick_inside(stream, |before, _| reads_differently(before, code))
        }
        CodeDifference::StartMarks => {
            reading.pick(0, 0, &reading.exons, stream, |before, after| {
                starts_differently(before, code) || starts_differently(after, code)
            })
        }
        CodeDifference::None => None,
    }
}

/// What a made-codon category asks of the codon its substitution starts from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Made {
    /// A sense codon whose residue the substitution keeps.
    Synonymous,
    /// A sense codon whose residue the substitution changes into another residue.
    Missense,
    /// A stop codon.
    Stop,
}

impl Made {
    /// Whether the substitution that turns `before` into `after` serves the category
    /// under `code`: the code translates `before` as the standard code does and `after`
    /// differently.
    fn admits(self, before: &[u8; 3], after: &[u8; 3], code: u8) -> bool {
        if reads_differently(before, code) || !reads_differently(after, code) {
            return false;
        }
        let was = translate_codon_with_table(before, code);
        let is = translate_codon_with_table(after, code);
        match self {
            Made::Synonymous => was == is && !matches!(was, b'*' | b'X'),
            Made::Missense => {
                was != is && !matches!(was, b'*' | b'X') && !matches!(is, b'*' | b'X')
            }
            Made::Stop => was == b'*',
        }
    }

    /// Whether a substitution of one base of some codon serves the category under `code`.
    fn has_a_codon(self, code: u8) -> bool {
        every_codon().any(|before| {
            (0..3).any(|offset| {
                BASES.into_iter().any(|base| {
                    let mut after = before;
                    after[offset] = base;
                    self.admits(&before, &after, code)
                })
            })
        })
    }
}

/// The made-codon categories of `code`, in the summary's order: each one that a
/// substitution of one base of some codon serves. A code that translates every codon as
/// the standard code does, or that the annotator has no table for, has none.
pub fn made_codon_categories(code: u8) -> Vec<Category> {
    if code_difference(code) != CodeDifference::Residues {
        return Vec::new();
    }
    [
        (Made::Synonymous, Category::MadeCodonSynonymous(code)),
        (Made::Missense, Category::MadeCodonMissense(code)),
        (Made::Stop, Category::MadeCodonStop(code)),
    ]
    .into_iter()
    .filter(|(made, _)| made.has_a_codon(code))
    .map(|(_, category)| category)
    .collect()
}

/// A substitution that keeps the residue of a sense codon and makes a codon the
/// transcript's code reads differently from the standard code.
pub fn made_codon_synonymous(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    let code = reading.code();
    reading.pick_inside(stream, |before, after| {
        Made::Synonymous.admits(before, after, code)
    })
}

/// A substitution that changes the residue of a sense codon and makes a sense codon the
/// transcript's code reads differently from the standard code.
pub fn made_codon_missense(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    let code = reading.code();
    reading.pick_inside(stream, |before, after| {
        Made::Missense.admits(before, after, code)
    })
}

/// A substitution in the stop codon ending the coding sequence that makes a codon the
/// transcript's code reads differently from the standard code. The codon made is a sense
/// codon or, under a code that reads a sense codon of the standard code as a stop codon,
/// a stop codon.
pub fn made_codon_stop(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    if reading.coding.sequence.len() % 3 != 0 || reading.codons < 2 {
        return None;
    }
    let code = reading.code();
    let last = reading.codons - 1;
    reading.pick(last, last, &reading.exons, stream, |before, after| {
        Made::Stop.admits(before, after, code)
    })
}

/// A substitution inside the first codon of the coding sequence.
pub fn start_codon(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    reading.pick(0, 0, &reading.exons, stream, |_, _| true)
}

/// A substitution that turns a sense codon into a stop codon of the transcript's code.
pub fn stop_gained(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    reading.pick_inside(stream, |before, after| {
        !matches!(reading.residue(before), b'*' | b'X') && reading.residue(after) == b'*'
    })
}

/// A substitution that turns the stop codon ending the coding sequence into a sense
/// codon of the transcript's code.
pub fn stop_lost(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    if reading.coding.sequence.len() % 3 != 0 || reading.codons < 2 {
        return None;
    }
    let last = reading.codons - 1;
    reading.pick(last, last, &reading.exons, stream, |before, after| {
        reading.residue(before) == b'*' && !matches!(reading.residue(after), b'*' | b'X')
    })
}

/// A substitution that changes the residue a codon translates to, no stop codon before
/// or after.
pub fn missense(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    reading.pick_inside(stream, |before, after| {
        let (was, is) = (reading.residue(before), reading.residue(after));
        was != is && !matches!(was, b'*' | b'X') && !matches!(is, b'*' | b'X')
    })
}

/// A substitution that leaves the residue a sense codon translates to as it was.
pub fn synonymous(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    reading.pick_inside(stream, |before, after| {
        let was = reading.residue(before);
        was == reading.residue(after) && !matches!(was, b'*' | b'X')
    })
}

/// A substitution at one of the two intron bases next to an exon: at the intron's 5' end
/// for a donor, at its 3' end for an acceptor.
fn splice_site(transcript: &Transcript, stream: &mut SplitMix64, donor: bool) -> Option<Change> {
    let introns: Vec<Span> = transcript
        .introns()
        .into_iter()
        .filter(|intron| intron.length() >= SPLICE_INTRON)
        .collect();
    if introns.is_empty() {
        return None;
    }
    let intron = introns[stream.below(introns.len() as u64) as usize];
    let step = stream.below(2);
    let at_low_end = donor == (transcript.strand == Strand::Forward);
    let position = if at_low_end {
        intron.start + step
    } else {
        intron.end - step
    };
    Some(unread_substitution(position, stream))
}

/// A substitution at one of the first two bases of an intron.
pub fn splice_donor(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    splice_site(transcript, stream, true)
}

/// A substitution at one of the last two bases of an intron.
pub fn splice_acceptor(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    splice_site(transcript, stream, false)
}

/// A substitution inside an intron, [`INTRON_MARGIN`] bases or more from both its ends.
pub fn intron(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let deep: Vec<Span> = transcript
        .introns()
        .into_iter()
        .filter(|intron| intron.length() > 2 * INTRON_MARGIN)
        .map(|intron| Span {
            start: intron.start + INTRON_MARGIN,
            end: intron.end - INTRON_MARGIN,
        })
        .collect();
    let position = draw_position(&deep, stream)?;
    Some(unread_substitution(position, stream))
}

/// A substitution in an untranslated region of a transcript whose codons the cache
/// gives, clear of its exon's ends: before the coding sequence or after it.
fn untranslated(transcript: &Transcript, stream: &mut SplitMix64, five: bool) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    let length = transcript.cdna_length();
    let region = if five {
        (reading.coding.cdna_start > 1).then(|| Span {
            start: 1,
            end: reading.coding.cdna_start - 1,
        })
    } else {
        (reading.coding.cdna_end < length).then(|| Span {
            start: reading.coding.cdna_end + 1,
            end: length,
        })
    }?;
    let parts: Vec<Span> = reading
        .clear
        .iter()
        .filter_map(|clear| overlap(*clear, region))
        .collect();
    let cdna = draw_position(&parts, stream)?;
    let position = transcript.genomic_position(cdna)?;
    Some(Change::substitution(
        position,
        transcript.forward_base(position),
        Alt::Other(stream.below(3) as u8),
    ))
}

/// A substitution in the 5' untranslated region.
pub fn five_prime_utr(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    untranslated(transcript, stream, true)
}

/// A substitution in the 3' untranslated region.
pub fn three_prime_utr(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    untranslated(transcript, stream, false)
}

/// An insertion of one or two bases inside the coding part of one exon.
pub fn frameshift_insertion(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    let anchor = draw_position(&reading.anchors(2), stream)?;
    let inserted = (0..stream.between(1, 2))
        .map(|_| drawn_base(stream))
        .collect();
    Some(Change {
        position: anchor,
        reference: reading.bases(anchor, 1)?,
        alt: Alt::Append(inserted),
        wraps_at: None,
    })
}

/// A deletion of one or two bases inside the coding part of one exon.
pub fn frameshift_deletion(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    let deleted = stream.between(1, 2);
    let anchor = draw_position(&reading.anchors(deleted + 1), stream)?;
    Some(Change {
        position: anchor,
        reference: reading.bases(anchor, deleted + 1)?,
        alt: Alt::First,
        wraps_at: None,
    })
}

/// A deletion of one whole codon inside the coding part of one exon.
pub fn inframe_deletion(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let reading = Reading::new(transcript)?;
    let anchors = reading.anchors(4);
    let total: u64 = anchors.iter().map(Span::length).sum();
    if total == 0 {
        return None;
    }
    // The deleted bases are a codon when the one read first, the lowest position on the
    // forward strand and the highest on the reverse, opens a codon.
    let opens_codon = |anchor: u64| {
        let read_first = match transcript.strand {
            Strand::Forward => anchor + 1,
            Strand::Reverse => anchor + 3,
        };
        reading
            .coding_index(read_first)
            .is_some_and(|index| index % 3 == 0)
    };
    let begin = stream.below(total);
    let anchor = (0..total.min(3))
        .filter_map(|step| nth_position(&anchors, (begin + step) % total))
        .find(|anchor| opens_codon(*anchor))?;
    Some(Change {
        position: anchor,
        reference: reading.bases(anchor, 4)?,
        alt: Alt::First,
        wraps_at: None,
    })
}

/// A substitution inside an exon of the transcript: at a base the cache gives where the
/// transcript has a coding model, anywhere in its exons where it has none. A transcript
/// with a coding model whose bases the cache does not give as the genome's has no such
/// place.
pub fn exon_substitution(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let Some(coding) = &transcript.coding else {
        let length = transcript.cdna_length();
        if length == 0 {
            return None;
        }
        let position = transcript.genomic_position(stream.between(1, length))?;
        return Some(unread_substitution(position, stream));
    };
    if !transcript.reads_codons() {
        return None;
    }
    let mut read = vec![Span {
        start: coding.cdna_start,
        end: coding.cdna_end,
    }];
    if coding.five_prime_utr.is_some() && coding.cdna_start > 1 {
        read.push(Span {
            start: 1,
            end: coding.cdna_start - 1,
        });
    }
    let length = transcript.cdna_length();
    if coding.three_prime_utr.is_some() && coding.cdna_end < length {
        read.push(Span {
            start: coding.cdna_end + 1,
            end: length,
        });
    }
    read.sort();
    let position = transcript.genomic_position(draw_position(&read, stream)?)?;
    Some(Change::substitution(
        position,
        Some(transcript.forward_base(position)?),
        Alt::Other(stream.below(3) as u8),
    ))
}

/// A substitution a drawn distance, [`FLANK`] bases at most, before the transcript's 5'
/// end; `None` where that falls before the sequence's first base.
pub fn upstream(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let distance = stream.between(1, FLANK);
    let position = match transcript.strand {
        Strand::Forward => transcript.start.checked_sub(distance).filter(|at| *at >= 1),
        Strand::Reverse => transcript.end.checked_add(distance),
    }?;
    Some(unread_substitution(position, stream))
}

/// A substitution a drawn distance, [`FLANK`] bases at most, past the transcript's 3'
/// end; `None` where that falls before the sequence's first base.
pub fn downstream(transcript: &Transcript, stream: &mut SplitMix64) -> Option<Change> {
    let distance = stream.between(1, FLANK);
    let position = match transcript.strand {
        Strand::Forward => transcript.end.checked_add(distance),
        Strand::Reverse => transcript.start.checked_sub(distance).filter(|at| *at >= 1),
    }?;
    Some(unread_substitution(position, stream))
}

/// A substitution at a drawn position of `gap`, a run of bases more than [`FLANK`] bases
/// from every transcript.
pub fn intergenic(gap: Span, stream: &mut SplitMix64) -> Change {
    let position = stream.between(gap.start, gap.end);
    unread_substitution(position, stream)
}

/// A substitution at `position`, for the first and the last base of a sequence.
pub fn edge_base(position: u64, stream: &mut SplitMix64) -> Change {
    unread_substitution(position, stream)
}

/// The shortest circular sequence the record across the origin is written on.
pub const ORIGIN_SEQUENCE: u64 = 4;

/// The record across the origin of a circular sequence of `length` bases: a deletion
/// written at the base before the last, whose reference allele is the last two bases and
/// then the first two, and whose alternate allele is its first base. Both engines read it
/// as a deletion of the three bases from the last one on, along a straight line.
pub fn origin(length: u64) -> Change {
    Change {
        position: length - 1,
        reference: vec![None; 4],
        alt: Alt::First,
        wraps_at: Some(length),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_category_has_its_name() {
        let mut names: Vec<String> = Category::FIXED.iter().map(Category::name).collect();
        names.extend(Category::of_code(2).iter().map(Category::name));
        let count = names.len();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), count);
        assert_eq!(count, 29);
        // The standard code, a code that differs from it in start codons alone and a
        // code without a table have no codon a substitution could make.
        assert_eq!(Category::of_code(1).len(), 4);
        assert_eq!(Category::of_code(11).len(), 4);
        assert_eq!(Category::of_code(7).len(), 4);
        assert_eq!(Category::MadeCodonStop(4).name(), "code4_made_codon_stop");
        assert!(Category::StopLost(4).rank() < Category::MadeCodonSynonymous(4).rank());
        assert!(Category::MadeCodonStop(4).rank() < Category::ChangedCodon(5).rank());
        assert_eq!(Category::StopLost(11).name(), "code11_stop_lost");
        assert!(Category::StopLost(2).rank() < Category::ChangedCodon(3).rank());
        assert!(Category::StopLost(24).rank() < Category::FirstSequence.rank());
        assert!(Category::Intron.rank() < Category::Origin.rank());
    }

    #[test]
    fn the_codes_differ_from_the_standard_code_as_their_tables_say() {
        assert_eq!(code_difference(1), CodeDifference::None);
        assert_eq!(code_difference(2), CodeDifference::Residues);
        assert_eq!(code_difference(3), CodeDifference::Residues);
        assert_eq!(code_difference(5), CodeDifference::Residues);
        assert_eq!(code_difference(11), CodeDifference::StartMarks);
        assert_eq!(code_difference(7), CodeDifference::None);
    }

    #[test]
    fn a_code_has_the_made_codon_categories_a_substitution_can_serve() {
        let names = |code: u8| -> Vec<String> {
            made_codon_categories(code)
                .iter()
                .map(Category::name)
                .collect()
        };
        assert_eq!(
            names(4),
            [
                "code4_made_codon_synonymous",
                "code4_made_codon_missense",
                "code4_made_codon_stop"
            ]
        );
        assert_eq!(names(2).len(), 3);
        // Code 12 reads CTG, a leucine codon of the standard code, as serine. No serine
        // codon and no stop codon is one substitution from it.
        assert_eq!(names(12), ["code12_made_codon_missense"]);
        // Code 23 reads TTA, a leucine codon of the standard code, as a stop codon, so
        // no substitution that makes it ends on a sense codon.
        assert_eq!(names(23), ["code23_made_codon_stop"]);
        for code in [1, 11, 24, 7, 0, 25] {
            assert!(names(code).is_empty(), "code {code}");
        }
    }

    #[test]
    fn the_other_bases_are_the_three_a_base_is_not() {
        let others: Vec<u8> = (0..3).map(|number| other_base(b'C', number)).collect();
        assert_eq!(others, b"AGT");
    }

    #[test]
    fn the_origin_record_runs_from_the_last_base_to_the_first() {
        let change = origin(500);
        assert_eq!(change.position, 499);
        assert_eq!(change.positions(), [499, 500, 1, 2]);
        assert_eq!(origin(u64::MAX).positions(), [u64::MAX - 1, u64::MAX, 1, 2]);
        assert_eq!(change.alt_for(b"GACT"), b"G");
    }

    #[test]
    fn positions_are_counted_along_spans() {
        let spans = [Span { start: 5, end: 6 }, Span { start: 10, end: 12 }];
        let all: Vec<Option<u64>> = (0..6).map(|n| nth_position(&spans, n)).collect();
        assert_eq!(all, [Some(5), Some(6), Some(10), Some(11), Some(12), None]);
    }
}
