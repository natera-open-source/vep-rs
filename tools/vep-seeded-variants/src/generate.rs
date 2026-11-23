// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The run over a cache. The sequences are read in the cache's order, one sequence's
//! transcripts in memory at a time; every transcript of a sequence the run places
//! records on is offered to each category it can serve, a category keeps its candidates
//! in a bounded reservoir, and at the end the kept candidates become the records and
//! the summary.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt;

use vep_core::codon::known_table;
use vep_core::synonyms::ChromosomeSynonyms;

use crate::cache::{Cache, CacheError, SequenceRegion, Span, Strand, Transcript, Transcripts};
use crate::fasta::Fasta;
use crate::random::{keyed_hash, keyed_stream, SplitMix64};
use crate::reservoir::{Key, Reservoir};
use crate::select::{self, Category, Change, CodeDifference, FLANK};
use crate::vcf::{Record, SummaryLine};

/// How many records a category gets in the wide mode unless the caller names a number.
pub const WIDE_RECORDS: usize = 12;

/// How many organelle or plasmid sequences get a record in the light mode.
pub const ORGANELLE_SEQUENCES: usize = 3;

/// The most sequences a cache holds for the rule "every sequence but the longest is an
/// organelle or a plasmid" to apply.
pub const FEW_SEQUENCES: usize = 12;

/// How many candidates a category keeps beyond the records it is asked for, which take
/// the place of a candidate that gives the record another one already gave.
const SPARE: usize = 4;

/// How many draws one transcript gets for one category before it is left out.
const DRAWS: usize = 8;

/// How many records each category is asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Records {
    /// One record a category: about thirty records a genome.
    Light,
    /// This many records a category, spread over sequences and biotypes.
    Wide(usize),
}

impl Records {
    fn want(&self) -> usize {
        match self {
            Records::Light => 1,
            Records::Wide(count) => *count,
        }
    }

    /// The mode as the output's header states it.
    pub fn describe(&self) -> String {
        match self {
            Records::Light => "light".to_owned(),
            Records::Wide(count) => format!("wide per_category={count}"),
        }
    }
}

/// What a run is asked for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Options {
    pub seed: u64,
    pub records: Records,
    /// The sequences to place records on; empty for every sequence of the cache. The
    /// rule for organelle and plasmid sequences reads every sequence of the cache
    /// whichever are named.
    pub sequences: Vec<String>,
}

/// What a run gives: the records in output order and the summary, one line a category.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Generated {
    pub records: Vec<Record>,
    pub summary: Vec<SummaryLine>,
}

/// Why a run stopped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Failure(pub String);

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Failure {}

impl From<CacheError> for Failure {
    fn from(error: CacheError) -> Failure {
        Failure(error.to_string())
    }
}

/// A kept candidate with its alleles settled.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Pick {
    /// The sequence by its number in the cache's reading order.
    sequence: usize,
    position: u64,
    reference: Vec<u8>,
    alt: Vec<u8>,
}

/// The picks without those that repeat an earlier one, `want` at most.
fn distinct(picks: Vec<Pick>, want: usize) -> Vec<Pick> {
    let mut kept: Vec<Pick> = Vec::new();
    for pick in picks {
        if kept.len() < want && !kept.contains(&pick) {
            kept.push(pick);
        }
    }
    kept
}

/// What the cache says of the base at a position.
enum Base {
    /// A transcript's cached sequence gives it.
    Known(u8),
    /// The position lies in a coding sequence the cache does not give as the genome's.
    Unread,
    /// No coding sequence covers the position.
    Free,
}

/// One sequence with its transcripts in memory, sorted by start.
struct Place<'a> {
    index: usize,
    name: &'a str,
    transcripts: &'a [Transcript],
    /// For each transcript, the highest end among it and those before it.
    reach: Vec<u64>,
    /// The last position known to exist: the sequence's length where the side file gives
    /// it, else the highest transcript end.
    bound: u64,
    /// The last position that takes no record, 0 where every position may: a transcript
    /// stored across the origin keeps records away from the sequence's first base up to
    /// [`FLANK`] bases past its end. There Ensembl VEP writes a row for that transcript on
    /// which no predicate holds, under the default term
    /// (`BaseVariationFeatureOverlapAllele.pm:285`).
    kept_clear: u64,
    /// The first position from which no record is taken up to the sequence's last base,
    /// `None` where every position may: the same transcript keeps records away from
    /// the [`FLANK`] bases before its first base, and from the base before those, since
    /// an insertion after that base stands among them. There Ensembl VEP writes no row
    /// for that transcript, whose stored span ends among the sequence's first bases
    /// (`AnnotationType/Transcript.pm:103-108`). It is 0, the whole sequence, where
    /// neither a mapper pair nor the sequence's length says where the transcript begins.
    kept_clear_from: Option<u64>,
    /// The name the reference file knows the sequence by.
    fasta_name: Option<String>,
}

impl<'a> Place<'a> {
    fn new(
        index: usize,
        name: &'a str,
        transcripts: &'a [Transcript],
        region: Option<&SequenceRegion>,
        fasta_name: Option<String>,
    ) -> Place<'a> {
        let mut reach = Vec::with_capacity(transcripts.len());
        let mut far = 0;
        let mut kept_clear = 0;
        let mut kept_clear_from: Option<u64> = None;
        for transcript in transcripts {
            far = far.max(transcript.end);
            reach.push(far);
            if let Some(before_origin) = &transcript.before_origin {
                kept_clear = kept_clear.max(transcript.end.saturating_add(FLANK));
                let from = before_origin
                    .first_base(region.map(|region| region.length))
                    .map_or(0, |first| first.saturating_sub(FLANK + 1));
                kept_clear_from = Some(kept_clear_from.map_or(from, |lowest| lowest.min(from)));
            }
        }
        Place {
            index,
            name,
            transcripts,
            reach,
            bound: region.map_or(far.max(1), |region| region.length),
            kept_clear,
            kept_clear_from,
            fasta_name,
        }
    }

    /// The highest transcript end, 0 on a sequence without transcripts.
    fn extent(&self) -> u64 {
        self.reach.last().copied().unwrap_or(0)
    }

    /// The transcripts whose span holds `position`.
    fn overlapping(&self, position: u64) -> impl Iterator<Item = &'a Transcript> + '_ {
        let after = self
            .transcripts
            .partition_point(|transcript| transcript.start <= position);
        (0..after)
            .rev()
            .take_while(move |at| self.reach[*at] >= position)
            .map(|at| &self.transcripts[at])
            .filter(move |transcript| transcript.end >= position)
    }

    fn in_transcript(&self, position: u64) -> bool {
        self.overlapping(position).next().is_some()
    }

    /// The base at `position` as the sequence's transcripts give it. Two transcripts
    /// that give different bases leave it unread.
    fn base(&self, position: u64) -> Base {
        let mut known: Option<u8> = None;
        let mut unread = false;
        for transcript in self.overlapping(position) {
            let Some(coding) = &transcript.coding else {
                continue;
            };
            let Some(cdna) = transcript.cdna_position(position) else {
                continue;
            };
            match transcript.forward_base(position) {
                Some(base) => match known {
                    Some(other) if other != base => return Base::Unread,
                    _ => known = Some(base),
                },
                None => unread |= (coding.cdna_start..=coding.cdna_end).contains(&cdna),
            }
        }
        match known {
            Some(base) => Base::Known(base),
            None if unread => Base::Unread,
            None => Base::Free,
        }
    }

    /// The runs of bases more than [`FLANK`] bases from every transcript, up to the
    /// last position known to exist that takes a record.
    fn gaps(&self) -> Vec<Span> {
        let last = self
            .kept_clear_from
            .map_or(self.bound, |from| self.bound.min(from.saturating_sub(1)));
        let mut gaps = Vec::new();
        let mut covered = self.kept_clear;
        for transcript in self.transcripts {
            let low = transcript.start.saturating_sub(FLANK).max(1);
            if low > covered + 1 {
                gaps.push(Span {
                    start: covered + 1,
                    end: low - 1,
                });
            }
            covered = covered.max(transcript.end.saturating_add(FLANK));
        }
        if covered < last {
            gaps.push(Span {
                start: covered + 1,
                end: last,
            });
        }
        gaps.into_iter()
            .filter(|gap| gap.start <= last)
            .map(|gap| Span {
                start: gap.start,
                end: gap.end.min(last),
            })
            .collect()
    }
}

/// What one candidate is offered as.
struct Offer<'a> {
    /// The name the candidate's priority and stream are keyed by.
    label: &'a str,
    key: Key,
    biotype: &'a str,
    /// The record must lie outside every transcript.
    outside: bool,
}

/// How a candidate came out.
enum Settled {
    Pick(Pick),
    /// The transcript holds no such place, or no draw found one the cache can vouch for.
    NoPlace,
    /// The reference file is masked there, or its bases are not the cache's.
    PassedOver,
}

/// Makes a candidate and settles its alleles. The reference bases the change leaves open
/// are taken from the sequence's other transcripts where one gives them; a place inside a
/// coding sequence that no transcript gives is drawn again. With a reference file every
/// base is read from it instead, and must be one of the four bases and equal to every
/// base the cache gives; without one an open base is drawn from the stream.
fn settle(
    place: &Place,
    fasta: &mut Option<&mut Fasta>,
    outside: bool,
    stream: &mut SplitMix64,
    make: &dyn Fn(&mut SplitMix64) -> Option<Change>,
) -> Settled {
    for _ in 0..DRAWS {
        let Some(change) = make(stream) else {
            return Settled::NoPlace;
        };
        let positions = change.positions();
        if change.position <= place.kept_clear
            || positions.iter().any(|at| {
                *at < 1
                    || *at > place.bound
                    || place.kept_clear_from.is_some_and(|from| *at >= from)
                    || (outside && place.in_transcript(*at))
            })
        {
            continue;
        }
        let mut reference = change.reference.clone();
        let mut unread = false;
        for (base, at) in reference.iter_mut().zip(&positions) {
            if base.is_none() {
                match place.base(*at) {
                    Base::Known(known) => *base = Some(known),
                    Base::Unread => unread = true,
                    Base::Free => {}
                }
            }
        }
        if unread {
            continue;
        }
        let reference: Vec<u8> = match fasta.as_deref_mut() {
            Some(fasta) => {
                let Some(name) = &place.fasta_name else {
                    return Settled::PassedOver;
                };
                let mut genome = Vec::with_capacity(positions.len());
                for at in &positions {
                    match fasta.bases(name, *at, 1).ok().as_deref() {
                        Some([base]) if b"ACGT".contains(base) => genome.push(*base),
                        _ => return Settled::PassedOver,
                    }
                }
                let unlike = reference
                    .iter()
                    .zip(&genome)
                    .any(|(cached, read)| cached.is_some_and(|cached| cached != *read));
                if unlike {
                    return Settled::PassedOver;
                }
                genome
            }
            None => reference
                .into_iter()
                .map(|base| base.unwrap_or_else(|| select::drawn_base(stream)))
                .collect(),
        };
        let alt = change.alt_for(&reference);
        return Settled::Pick(Pick {
            sequence: place.index,
            position: change.position,
            reference,
            alt,
        });
    }
    Settled::NoPlace
}

/// Offers one candidate to a reservoir: its priority is a hash of the seed and of what
/// it is, and it is made only where that priority has a place.
fn offer(
    reservoir: &mut Reservoir<Pick>,
    passed_over: &mut u64,
    seed: u64,
    place: &Place,
    fasta: &mut Option<&mut Fasta>,
    offer: Offer,
    make: &dyn Fn(&mut SplitMix64) -> Option<Change>,
) {
    let (start, end) = (offer.key.start.to_le_bytes(), offer.key.end.to_le_bytes());
    let parts = |purpose: &'static [u8]| -> [&[u8]; 6] {
        [
            purpose,
            offer.label.as_bytes(),
            place.name.as_bytes(),
            offer.key.id.as_bytes(),
            &start,
            &end,
        ]
    };
    let priority = keyed_hash(seed, &parts(b"priority"));
    let stratum = (place.index, offer.biotype.to_owned());
    if !reservoir.admits(priority, &offer.key, &stratum) {
        return;
    }
    let mut stream = keyed_stream(seed, &parts(b"change"));
    match settle(place, fasta, offer.outside, &mut stream, make) {
        Settled::Pick(pick) => reservoir.insert(priority, offer.key, stratum, pick),
        Settled::PassedOver => *passed_over += 1,
        Settled::NoPlace => {}
    }
}

type Maker = fn(&Transcript, &mut SplitMix64) -> Option<Change>;

/// The categories a transcript with a coding sequence is offered to whatever its
/// genetic code.
fn coding_categories(transcript: &Transcript) -> [(Category, Maker); 11] {
    let code = transcript.code();
    let strand = transcript.strand;
    [
        (Category::ChangedCodon(code), select::changed_codon),
        (Category::StartCodon(code), select::start_codon),
        (Category::StopGained(code), select::stop_gained),
        (Category::StopLost(code), select::stop_lost),
        (Category::Missense(strand), select::missense),
        (Category::Synonymous(strand), select::synonymous),
        (Category::FivePrimeUtr, select::five_prime_utr),
        (Category::ThreePrimeUtr, select::three_prime_utr),
        (Category::FrameshiftInsertion, select::frameshift_insertion),
        (Category::FrameshiftDeletion, select::frameshift_deletion),
        (Category::InframeDeletion, select::inframe_deletion),
    ]
}

/// The made-codon categories of `code`, which a transcript with a coding sequence on
/// that code is offered to besides.
fn made_codon_categories(code: u8) -> Vec<(Category, Maker)> {
    select::made_codon_categories(code)
        .into_iter()
        .map(|category| {
            let maker: Maker = match category {
                Category::MadeCodonSynonymous(_) => select::made_codon_synonymous,
                Category::MadeCodonMissense(_) => select::made_codon_missense,
                _ => select::made_codon_stop,
            };
            (category, maker)
        })
        .collect()
}

/// The categories a transcript of more than one exon is offered to.
const SPLICED: [(Category, Maker); 3] = [
    (Category::SpliceDonor, select::splice_donor),
    (Category::SpliceAcceptor, select::splice_acceptor),
    (Category::Intron, select::intron),
];

/// The categories of a place beside a transcript.
const BESIDE: [(Category, Maker); 2] = [
    (Category::Upstream, select::upstream),
    (Category::Downstream, select::downstream),
];

/// The label of a sequence's own sample of exon substitutions.
const SAMPLE: &str = "sequence";

/// A record before its ordinal is known.
struct Draft {
    category: Category,
    /// The sequence's name as the record spells it.
    chrom: String,
    pick: Pick,
}

struct Run<'a> {
    seed: u64,
    want: usize,
    light: bool,
    fasta: Option<&'a mut Fasta>,
    /// Every sequence of the cache, in reading order.
    names: &'a [String],
    /// The sequences this run places records on, by number.
    selected: Vec<usize>,
    reservoirs: BTreeMap<Category, Reservoir<Pick>>,
    passed_over: BTreeMap<Category, u64>,
    /// The kept samples of exon substitutions, by sequence, each with how many
    /// candidates the reference file made it pass over.
    samples: BTreeMap<usize, (Reservoir<Pick>, u64)>,
    /// The sequences whose sample is kept whatever their codes: the first, the last and
    /// those chosen for a synonym.
    kept: BTreeSet<usize>,
    /// The names a chosen sequence can be written by.
    synonyms: BTreeMap<usize, Vec<String>>,
    /// The sequences chosen for a synonym, best first.
    chosen: Vec<usize>,
    /// How many transcripts of the cache name each genetic code, on every sequence.
    code_counts: BTreeMap<u8, u64>,
    /// The codes that transcripts with a coding sequence name, on the sequences this
    /// run places records on.
    coding_codes: BTreeSet<u8>,
    /// Those of them under which a transcript's codons are read.
    read_codes: BTreeSet<u8>,
    /// The made-codon categories of each code of `coding_codes`, each with its
    /// function.
    made_codon: BTreeMap<u8, Vec<(Category, Maker)>>,
    /// For each code, the first sequences in reading order that name it, of those this
    /// run places records on.
    by_code: BTreeMap<u8, Vec<usize>>,
    /// In a cache of few sequences, the length and the highest transcript end of each
    /// sequence of the cache.
    extents: BTreeMap<usize, (Option<u64>, u64)>,
    forward_seen: bool,
    reverse_seen: bool,
    spliced_seen: bool,
    length_seen: bool,
    /// A sequence the side file calls circular is long enough for the origin record.
    circular_seen: bool,
}

impl Run<'_> {
    /// Whether the cache holds few sequences, however many of them this run places
    /// records on.
    fn few(&self) -> bool {
        self.names.len() <= FEW_SEQUENCES
    }

    fn organelle_cap(&self) -> usize {
        if self.light {
            ORGANELLE_SEQUENCES
        } else {
            self.want
        }
    }

    fn offer(
        &mut self,
        category: Category,
        place: &Place,
        key: Key,
        biotype: &str,
        outside: bool,
        make: &dyn Fn(&mut SplitMix64) -> Option<Change>,
    ) {
        let want = self.want;
        let reservoir = self
            .reservoirs
            .entry(category)
            .or_insert_with(|| Reservoir::new(want.saturating_add(SPARE)));
        let passed_over = self.passed_over.entry(category).or_default();
        let label = category.name();
        let offered = Offer {
            label: &label,
            key,
            biotype,
            outside,
        };
        offer(
            reservoir,
            passed_over,
            self.seed,
            place,
            &mut self.fasta,
            offered,
            make,
        );
    }

    /// Offers one sequence's transcripts and places to every category.
    fn sequence(
        &mut self,
        index: usize,
        transcripts: &[Transcript],
        region: Option<&SequenceRegion>,
        fasta_name: Option<String>,
    ) {
        let name = &self.names[index];
        let place = Place::new(index, name, transcripts, region, fasta_name);
        let mut sample = Reservoir::new(self.want.saturating_mul(3).saturating_add(SPARE));
        let mut sample_passed_over = 0;
        let mut codes = BTreeSet::new();
        for transcript in transcripts {
            let code = transcript.code();
            *self.code_counts.entry(code).or_default() += 1;
            codes.insert(code);
            match transcript.strand {
                Strand::Forward => self.forward_seen = true,
                Strand::Reverse => self.reverse_seen = true,
            }
            let key = Key {
                sequence: index,
                start: transcript.start,
                end: transcript.end,
                id: transcript.stable_id.clone(),
            };
            let biotype = transcript.biotype.as_str();
            if transcript.coding.is_some() {
                self.coding_codes.insert(code);
                if transcript.reads_codons() {
                    self.read_codes.insert(code);
                }
                let difference = select::code_difference(code);
                for (category, maker) in coding_categories(transcript) {
                    if category == Category::ChangedCodon(code)
                        && difference == CodeDifference::None
                    {
                        continue;
                    }
                    self.offer(category, &place, key.clone(), biotype, false, &|stream| {
                        maker(transcript, stream)
                    });
                }
                let made = self
                    .made_codon
                    .entry(code)
                    .or_insert_with(|| made_codon_categories(code))
                    .clone();
                for (category, maker) in made {
                    self.offer(category, &place, key.clone(), biotype, false, &|stream| {
                        maker(transcript, stream)
                    });
                }
            }
            if transcript.exons.len() > 1 {
                self.spliced_seen = true;
                for (category, maker) in SPLICED {
                    self.offer(category, &place, key.clone(), biotype, false, &|stream| {
                        maker(transcript, stream)
                    });
                }
            }
            for (category, maker) in BESIDE {
                self.offer(category, &place, key.clone(), biotype, true, &|stream| {
                    maker(transcript, stream)
                });
            }
            let offered = Offer {
                label: SAMPLE,
                key,
                biotype,
                outside: false,
            };
            offer(
                &mut sample,
                &mut sample_passed_over,
                self.seed,
                &place,
                &mut self.fasta,
                offered,
                &|stream| select::exon_substitution(transcript, stream),
            );
        }
        let at = |start: u64, end: u64| Key {
            sequence: index,
            start,
            end,
            id: String::new(),
        };
        for gap in place.gaps() {
            self.offer(
                Category::Intergenic,
                &place,
                at(gap.start, gap.end),
                "",
                true,
                &|stream| Some(select::intergenic(gap, stream)),
            );
        }
        self.offer(
            Category::FirstBase,
            &place,
            at(1, 1),
            "",
            false,
            &|stream| Some(select::edge_base(1, stream)),
        );
        if let Some(region) = region {
            self.length_seen = true;
            let length = region.length;
            self.offer(
                Category::LastBase,
                &place,
                at(length, length),
                "",
                false,
                &|stream| Some(select::edge_base(length, stream)),
            );
            if region.circular == Some(true) && length >= select::ORIGIN_SEQUENCE {
                self.circular_seen = true;
                self.offer(
                    Category::Origin,
                    &place,
                    at(length, length),
                    "",
                    false,
                    &|_| Some(select::origin(length)),
                );
            }
        }
        if sample.is_empty() {
            // A sequence none of whose transcripts gives a place still gets a record:
            // a substitution near its first base.
            let last = place.bound.min(100);
            let offered = Offer {
                label: SAMPLE,
                key: at(1, last),
                biotype: "",
                outside: false,
            };
            offer(
                &mut sample,
                &mut sample_passed_over,
                self.seed,
                &place,
                &mut self.fasta,
                offered,
                &|stream| Some(select::edge_base(stream.between(1, last), stream)),
            );
        }
        let mut keep = self.kept.contains(&index) || self.few();
        let cap = self.organelle_cap();
        for code in codes {
            let named = self.by_code.entry(code).or_default();
            if named.len() < cap {
                named.push(index);
                keep = true;
            }
        }
        if self.few() {
            self.extents
                .insert(index, (region.map(|region| region.length), place.extent()));
        }
        if keep {
            self.samples.insert(index, (sample, sample_passed_over));
        }
    }

    /// Reads a sequence this run places no record on, for what the rule of
    /// `Run::organelles` counts over the whole cache: the genetic codes its transcripts
    /// name and, in a cache of few sequences, its length and highest transcript end.
    fn pass_by(
        &mut self,
        index: usize,
        transcripts: Transcripts,
        region: Option<&SequenceRegion>,
    ) -> Result<(), CacheError> {
        let mut far = 0;
        for transcript in transcripts {
            let transcript = transcript?;
            *self.code_counts.entry(transcript.code()).or_default() += 1;
            far = far.max(transcript.end);
        }
        if self.few() {
            self.extents
                .insert(index, (region.map(|region| region.length), far));
        }
        Ok(())
    }

    /// The organelle and plasmid sequences this run places records on, in reading
    /// order. The rule reads the whole cache: every sequence one of whose transcripts
    /// names another genetic code than most transcripts of the cache; where no sequence
    /// of the cache does and the cache holds few sequences, every sequence but the
    /// cache's longest.
    fn organelles(&self) -> Vec<usize> {
        let majority = self
            .code_counts
            .iter()
            .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
            .map(|(code, _)| *code);
        let mut found: BTreeSet<usize> = self
            .by_code
            .iter()
            .filter(|(code, _)| Some(**code) != majority)
            .flat_map(|(_, sequences)| sequences.iter().copied())
            .collect();
        let one_code = self.code_counts.keys().all(|code| Some(*code) == majority);
        if one_code && self.few() && self.names.len() >= 2 {
            let every_length = self.extents.values().all(|(length, _)| length.is_some());
            let size = |index: &usize| match self.extents.get(index) {
                Some((Some(length), _)) if every_length => *length,
                Some((_, extent)) => *extent,
                None => 0,
            };
            let longest = (0..self.names.len()).max_by(|a, b| size(a).cmp(&size(b)).then(b.cmp(a)));
            found = self
                .selected
                .iter()
                .copied()
                .filter(|index| Some(*index) != longest)
                .collect();
        }
        found.into_iter().take(self.organelle_cap()).collect()
    }

    /// The picks of a sequence's sample for one of its three users, which start at
    /// different picks so that two categories on one sequence get different records
    /// where the sequence has enough.
    fn sample_picks(&self, index: usize, user: usize) -> Vec<Pick> {
        let Some((reservoir, _)) = self.samples.get(&index) else {
            return Vec::new();
        };
        let mut picks = distinct(reservoir.picks(), self.want.saturating_mul(3));
        if !picks.is_empty() {
            let shift = (user * (self.want % picks.len())) % picks.len();
            picks.rotate_left(shift);
        }
        picks
    }

    /// `total` picks at most from the samples of `sequences`, taken in turns.
    fn in_turns(&self, sequences: &[usize], user: usize, total: usize) -> Vec<Pick> {
        let pools: Vec<Vec<Pick>> = sequences
            .iter()
            .map(|index| self.sample_picks(*index, user))
            .collect();
        let deepest = pools.iter().map(Vec::len).max().unwrap_or(0);
        let mut picks = Vec::new();
        for round in 0..deepest {
            for pool in &pools {
                if picks.len() < total {
                    picks.extend(pool.get(round).cloned());
                }
            }
        }
        picks
    }

    fn sample_passed_over(&self, sequences: &[usize]) -> u64 {
        sequences
            .iter()
            .filter_map(|index| self.samples.get(index))
            .map(|(_, passed_over)| passed_over)
            .sum()
    }

    /// Why a category wrote no record, or fewer than it wanted.
    fn reason(
        &self,
        category: Category,
        wanted: usize,
        written: usize,
        passed_over: u64,
    ) -> String {
        let reason = match category {
            _ if wanted == 0 => match category {
                Category::ChangedCodon(_) => "the standard code",
                _ => "no organelle or plasmid sequence",
            },
            _ if written >= wanted => "-",
            _ if written > 0 => "few candidates",
            _ if passed_over > 0 => "reference masked",
            Category::ChangedCodon(code)
            | Category::StartCodon(code)
            | Category::StopGained(code)
            | Category::StopLost(code)
                if !known_table(code) =>
            {
                "no table for this code"
            }
            Category::ChangedCodon(_) => "no codon this code reads differently",
            Category::StopLost(_) => "no complete stop codon",
            Category::MadeCodonSynonymous(code)
            | Category::MadeCodonMissense(code)
            | Category::MadeCodonStop(code) => {
                if self.read_codes.contains(&code) {
                    "no substitution makes a codon this code reads differently"
                } else {
                    "no complete coding sequence"
                }
            }
            Category::StartCodon(_)
            | Category::StopGained(_)
            | Category::FrameshiftInsertion
            | Category::FrameshiftDeletion
            | Category::InframeDeletion => "no complete coding sequence",
            Category::LastSequence if self.selected.len() == 1 => "one sequence",
            Category::Synonym if self.chosen.is_empty() => "no synonym",
            Category::FirstSequence
            | Category::LastSequence
            | Category::Organelle
            | Category::Synonym
            | Category::FirstBase => "no place on the sequence",
            Category::LastBase if self.length_seen => "no place on the sequence",
            Category::Origin if self.circular_seen => "no place on the sequence",
            Category::Missense(strand) | Category::Synonymous(strand) => {
                let seen = match strand {
                    Strand::Forward => self.forward_seen,
                    Strand::Reverse => self.reverse_seen,
                };
                if seen {
                    "no complete coding sequence on this strand"
                } else {
                    "no transcript on this strand"
                }
            }
            Category::SpliceDonor | Category::SpliceAcceptor | Category::Intron => {
                if self.spliced_seen {
                    "no intron long enough"
                } else {
                    "no multi-exon transcript"
                }
            }
            Category::FivePrimeUtr | Category::ThreePrimeUtr => "no untranslated region",
            Category::Upstream | Category::Downstream => "no room beside a transcript",
            Category::Intergenic => "no gap between genes",
            Category::LastBase => "no sequence length",
            Category::Origin => {
                if self.length_seen {
                    "no circular sequence"
                } else {
                    "no sequence length"
                }
            }
        };
        reason.to_owned()
    }

    /// Turns the kept candidates into the records and the summary.
    fn finish(self) -> Generated {
        let mut drafts: Vec<Draft> = Vec::new();
        let mut wanted: BTreeMap<Category, usize> = BTreeMap::new();
        let mut passed_over = self.passed_over.clone();
        let name_of = |index: usize| self.names[index].clone();
        for (category, reservoir) in &self.reservoirs {
            for pick in distinct(reservoir.picks(), self.want) {
                drafts.push(Draft {
                    category: *category,
                    chrom: name_of(pick.sequence),
                    pick,
                });
            }
        }
        let mut place = |category: Category, picks: Vec<Pick>, over: &[usize]| {
            passed_over.insert(category, self.sample_passed_over(over));
            for pick in picks {
                drafts.push(Draft {
                    category,
                    chrom: name_of(pick.sequence),
                    pick,
                });
            }
        };
        let first = self.selected[0];
        let last = self.selected[self.selected.len() - 1];
        place(
            Category::FirstSequence,
            self.in_turns(&[first], 0, self.want),
            &[first],
        );
        if last != first {
            place(
                Category::LastSequence,
                self.in_turns(&[last], 0, self.want),
                &[last],
            );
        }
        let organelles = self.organelles();
        let organelle_records = match (self.light, organelles.is_empty()) {
            (_, true) => 0,
            (true, false) => organelles.len(),
            (false, false) => self.want,
        };
        wanted.insert(Category::Organelle, organelle_records);
        place(
            Category::Organelle,
            self.in_turns(&organelles, 1, organelle_records),
            &organelles,
        );
        passed_over.insert(Category::Synonym, self.sample_passed_over(&self.chosen));
        for (ordinal, pick) in self
            .in_turns(&self.chosen, 2, self.want)
            .into_iter()
            .enumerate()
        {
            let sequence = &self.names[pick.sequence];
            let names = &self.synonyms[&pick.sequence];
            let drawn = keyed_hash(
                self.seed,
                &[
                    b"synonym name",
                    sequence.as_bytes(),
                    &(ordinal as u64).to_le_bytes(),
                ],
            );
            drafts.push(Draft {
                category: Category::Synonym,
                chrom: names[(drawn % names.len() as u64) as usize].clone(),
                pick,
            });
        }

        drafts.sort_by(|a, b| {
            let key = |draft: &'_ Draft| {
                (
                    draft.pick.sequence,
                    draft.pick.position,
                    draft.category.rank(),
                )
            };
            key(a)
                .cmp(&key(b))
                .then_with(|| a.pick.reference.cmp(&b.pick.reference))
                .then_with(|| a.pick.alt.cmp(&b.pick.alt))
                .then_with(|| a.chrom.cmp(&b.chrom))
        });
        let mut written: BTreeMap<Category, usize> = BTreeMap::new();
        let text = |bases: &[u8]| String::from_utf8_lossy(bases).into_owned();
        let records = drafts
            .into_iter()
            .map(|draft| {
                let ordinal = written.entry(draft.category).or_default();
                *ordinal += 1;
                Record {
                    sequence: draft.chrom,
                    position: draft.pick.position,
                    id: format!("{}_{}", draft.category.name(), ordinal),
                    reference: text(&draft.pick.reference),
                    alt: text(&draft.pick.alt),
                }
            })
            .collect();

        let categories = self
            .coding_codes
            .iter()
            .flat_map(|code| Category::of_code(*code))
            .chain(Category::FIXED);
        let summary = categories
            .map(|category| {
                let standard = matches!(category, Category::ChangedCodon(code)
                    if known_table(code) && select::code_difference(code) == CodeDifference::None);
                let wanted = match wanted.get(&category) {
                    Some(wanted) => *wanted,
                    None if standard => 0,
                    None => self.want,
                };
                let written = written.get(&category).copied().unwrap_or(0);
                let passed_over = passed_over.get(&category).copied().unwrap_or(0);
                SummaryLine {
                    category: category.name(),
                    wanted,
                    written,
                    passed_over,
                    reason: self.reason(category, wanted, written, passed_over),
                }
            })
            .collect();
        Generated { records, summary }
    }
}

/// For each of `selected`, the names `chr_synonyms.txt` gives the sequence that neither
/// the cache's own name nor the `chr` rule reaches: a name the cache holds no sequence
/// by, that the synonyms resolve to the sequence and that an empty table does not.
fn synonym_names(cache: &Cache, selected: &[usize]) -> BTreeMap<usize, Vec<String>> {
    let sequences = cache.sequences();
    let held: HashSet<&str> = sequences.iter().map(String::as_str).collect();
    let holds = |name: &str| held.contains(name);
    let mut table = ChromosomeSynonyms::new();
    let mut partners: HashMap<&str, Vec<&str>> = HashMap::new();
    for (first, later) in cache.synonym_pairs() {
        table.add_synonym(first, later);
        partners.entry(first).or_default().push(later);
        partners.entry(later).or_default().push(first);
    }
    let unlinked = ChromosomeSynonyms::new();
    let mut found = BTreeMap::new();
    for index in selected {
        let sequence = sequences[*index].as_str();
        let mut names: Vec<String> = partners
            .get(sequence)
            .into_iter()
            .flatten()
            .filter(|name| {
                !holds(name)
                    && table.source_chr_name(name, holds) == sequence
                    && unlinked.source_chr_name(name, holds) != sequence
            })
            .map(|name| (*name).to_owned())
            .collect();
        names.sort();
        names.dedup();
        if !names.is_empty() {
            found.insert(*index, names);
        }
    }
    found
}

/// Generates the records of `cache` under `options`. With `fasta`, every reference
/// allele is read from that file, and a read of it that fails stops the run.
pub fn generate(
    cache: &Cache,
    fasta: Option<&mut Fasta>,
    options: &Options,
) -> Result<Generated, Failure> {
    let names = cache.sequences();
    let mut selected: Vec<usize> = if options.sequences.is_empty() {
        (0..names.len()).collect()
    } else {
        let mut selected = Vec::new();
        for wanted in &options.sequences {
            match names.iter().position(|name| name == wanted) {
                Some(index) => selected.push(index),
                None => {
                    return Err(Failure(format!(
                        "{}: the cache holds no sequence {wanted}",
                        cache.dir().display()
                    )))
                }
            }
        }
        selected
    };
    selected.sort();
    selected.dedup();
    if selected.is_empty() {
        return Err(Failure(format!(
            "{}: the cache names no sequence",
            cache.dir().display()
        )));
    }
    let want = options.records.want();
    let seed = options.seed;

    let mut synonyms = synonym_names(cache, &selected);
    let mut ranked: Vec<usize> = synonyms.keys().copied().collect();
    ranked.sort_by_key(|index| {
        (
            keyed_hash(seed, &[b"synonym", names[*index].as_bytes()]),
            *index,
        )
    });
    let mut chosen = Vec::new();
    let mut bare = Vec::new();
    for index in ranked {
        if chosen.len() >= want {
            break;
        }
        if cache.shards(&names[index])?.is_empty() {
            bare.push(index);
        } else {
            chosen.push(index);
        }
    }
    chosen.extend(bare);
    chosen.truncate(want);
    synonyms.retain(|index, _| chosen.contains(index));

    let regions: HashMap<&str, &SequenceRegion> = cache
        .sequence_regions()
        .unwrap_or_default()
        .iter()
        .map(|region| (region.name.as_str(), region))
        .collect();
    let mut table = ChromosomeSynonyms::new();
    for (first, later) in cache.synonym_pairs() {
        table.add_synonym(first, later);
    }
    let mut kept: BTreeSet<usize> = chosen.iter().copied().collect();
    kept.insert(selected[0]);
    kept.insert(selected[selected.len() - 1]);
    let mut run = Run {
        seed,
        want,
        light: options.records == Records::Light,
        fasta,
        names,
        selected: selected.clone(),
        reservoirs: BTreeMap::new(),
        passed_over: BTreeMap::new(),
        samples: BTreeMap::new(),
        kept,
        synonyms,
        chosen,
        code_counts: BTreeMap::new(),
        coding_codes: BTreeSet::new(),
        read_codes: BTreeSet::new(),
        made_codon: BTreeMap::new(),
        by_code: BTreeMap::new(),
        extents: BTreeMap::new(),
        forward_seen: false,
        reverse_seen: false,
        spliced_seen: false,
        length_seen: false,
        circular_seen: false,
    };
    for (index, name) in names.iter().enumerate() {
        let region = regions.get(name.as_str()).copied();
        if selected.binary_search(&index).is_err() {
            run.pass_by(index, cache.transcripts(name)?, region)?;
            continue;
        }
        let mut transcripts = Vec::new();
        for transcript in cache.transcripts(name)? {
            transcripts.push(transcript?);
        }
        transcripts
            .sort_by(|a, b| (a.start, &a.stable_id, a.end).cmp(&(b.start, &b.stable_id, b.end)));
        let fasta_name = run.fasta.as_deref().and_then(|fasta| {
            let known = table.source_chr_name(name, |candidate| fasta.has(candidate));
            fasta.has(&known).then_some(known)
        });
        run.sequence(index, &transcripts, region, fasta_name);
        if let Some(failure) = run.fasta.as_deref().and_then(Fasta::failure) {
            return Err(Failure(failure.to_owned()));
        }
    }
    Ok(run.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::BeforeOrigin;
    use crate::select::Alt;

    /// A transcript without exons from `start` to `end`.
    fn transcript(
        id: &str,
        start: u64,
        end: u64,
        before_origin: Option<BeforeOrigin>,
    ) -> Transcript {
        Transcript {
            stable_id: id.to_owned(),
            gene_stable_id: None,
            start,
            end,
            strand: Strand::Forward,
            biotype: String::new(),
            exons: Vec::new(),
            codon_table: None,
            coding: None,
            cds_start_nf: false,
            cds_end_nf: false,
            edited: false,
            across_origin: before_origin.is_some(),
            before_origin,
        }
    }

    fn region(length: u64) -> SequenceRegion {
        SequenceRegion {
            name: "ring".to_owned(),
            length,
            circular: Some(true),
        }
    }

    /// A change of `bases` reference bases from `position` on.
    fn change(position: u64, bases: usize, alt: Alt) -> Change {
        Change {
            position,
            reference: vec![None; bases],
            alt,
            wraps_at: None,
        }
    }

    fn takes(place: &Place, change: Change) -> bool {
        let mut stream = SplitMix64::new(1);
        let settled = settle(place, &mut None, false, &mut stream, &|_| {
            Some(change.clone())
        });
        matches!(settled, Settled::Pick(_))
    }

    #[test]
    fn a_transcript_across_the_origin_keeps_both_its_sides_clear() {
        // A sequence of 30,000 bases. One transcript is stored from 199 bases before the
        // origin to base 700, so it begins at base 30,000 - 199 = 29,801; another lies
        // from 12,000 to 12,999.
        let across = BeforeOrigin {
            bases: 200,
            pair_start: Some(29_801),
        };
        let transcripts = [
            transcript("across", 1, 700, Some(across)),
            transcript("other", 12_000, 12_999, None),
        ];
        let ring = region(30_000);
        let place = Place::new(0, "ring", &transcripts, Some(&ring), None);
        // Clear to 700 + 5,000, and from the base before 29,801 - 5,000 = 24,801.
        assert_eq!(
            (place.kept_clear, place.kept_clear_from),
            (5_700, Some(24_800))
        );
        // More than 5,000 bases from both: after 5,700 up to 12,000 - 5,001, and after
        // 12,999 + 5,000 up to the last base that takes a record.
        assert_eq!(
            place.gaps(),
            [
                Span {
                    start: 5_701,
                    end: 6_999
                },
                Span {
                    start: 18_000,
                    end: 24_799
                }
            ]
        );

        let substitution = |position| change(position, 1, Alt::Other(0));
        assert!(!takes(&place, substitution(5_700)));
        assert!(takes(&place, substitution(5_701)));
        assert!(takes(&place, substitution(24_799)));
        assert!(!takes(&place, substitution(24_800)));
        assert!(!takes(&place, substitution(29_801)));
        assert!(!takes(&place, substitution(30_000)));
        // A deletion is refused by the last base of its reference allele.
        assert!(takes(&place, change(24_797, 3, Alt::First)));
        assert!(!takes(&place, change(24_798, 3, Alt::First)));
        // An insertion after base 24,800 would stand at 24,801.
        let insertion = |position| change(position, 1, Alt::Append(b"A".to_vec()));
        assert!(takes(&place, insertion(24_799)));
        assert!(!takes(&place, insertion(24_800)));
        // The record across the origin holds the last two bases and the first two.
        assert!(!takes(&place, select::origin(30_000)));

        // Without the transcript across the origin every such place takes its record.
        let plain = Place::new(0, "ring", &transcripts[1..], Some(&ring), None);
        assert_eq!((plain.kept_clear, plain.kept_clear_from), (0, None));
        for position in [1, 5_700, 24_800, 29_801, 30_000] {
            assert!(takes(&plain, substitution(position)), "{position}");
        }
        assert!(takes(&plain, select::origin(30_000)));
        assert_eq!(
            plain.gaps(),
            [
                Span {
                    start: 1,
                    end: 6_999
                },
                Span {
                    start: 18_000,
                    end: 30_000
                }
            ]
        );
    }

    #[test]
    fn the_stretch_among_the_last_bases_starts_where_the_transcript_is_said_to_begin() {
        let from = |pair_start: Option<u64>, length: Option<u64>| {
            let across = BeforeOrigin {
                bases: 200,
                pair_start,
            };
            let transcripts = [transcript("across", 1, 700, Some(across))];
            let ring = length.map(region);
            Place::new(0, "ring", &transcripts, ring.as_ref(), None).kept_clear_from
        };
        // The mapper pair alone, the length alone (30,000 + 1 - 200 = 29,801), and both.
        assert_eq!(from(Some(29_801), None), Some(24_800));
        assert_eq!(from(None, Some(30_000)), Some(24_800));
        assert_eq!(from(Some(29_801), Some(30_000)), Some(24_800));
        // Neither says where the transcript begins: the whole sequence is kept clear.
        assert_eq!(from(None, None), Some(0));
        // A transcript that begins within 5,001 bases of base 1 leaves no base free.
        assert_eq!(from(Some(5_001), None), Some(0));
        assert_eq!(from(Some(5_002), None), Some(1));

        // Two such transcripts: the stretch starts for the one that begins lower, at the
        // base before 27,001 - 5,000 = 22,001.
        let begins = |pair_start| {
            Some(BeforeOrigin {
                bases: 200,
                pair_start: Some(pair_start),
            })
        };
        let transcripts = [
            transcript("first", 1, 700, begins(29_801)),
            transcript("second", 1, 300, begins(27_001)),
        ];
        let place = Place::new(0, "ring", &transcripts, None, None);
        assert_eq!(
            (place.kept_clear, place.kept_clear_from),
            (5_700, Some(22_000))
        );
        assert!(place.gaps().is_empty());
        // The same with the two in the other order.
        let [first, second] = transcripts;
        let swapped = [second, first];
        let place = Place::new(0, "ring", &swapped, None, None);
        assert_eq!(
            (place.kept_clear, place.kept_clear_from),
            (5_700, Some(22_000))
        );
    }

    #[test]
    fn a_sequence_without_such_a_transcript_takes_a_record_at_the_highest_position() {
        // A side file can state a length of `u64::MAX` bases: with no transcript across
        // the origin the last base takes its record and the last gap ends at it.
        let transcripts = [transcript("other", 12_000, 12_999, None)];
        let ring = region(u64::MAX);
        let place = Place::new(0, "ring", &transcripts, Some(&ring), None);
        assert_eq!(place.kept_clear_from, None);
        assert!(takes(&place, change(u64::MAX, 1, Alt::Other(0))));
        assert!(takes(&place, select::origin(u64::MAX)));
        assert_eq!(place.gaps().last().map(|gap| gap.end), Some(u64::MAX));
    }
}
