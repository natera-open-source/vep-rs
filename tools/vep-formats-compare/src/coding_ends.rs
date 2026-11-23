// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What stands at the two ends of a cache transcript's coding sequence, as vep-rs reads the
//! transcript's cache record: whether a stop codon closes the coding sequence, whether a
//! start codon opens it, and where the genome holds the bases of both ends.
//!
//! vep-rs writes a stop codon's terms and a start codon's terms only where that codon
//! stands (`crates/vep-effects/src/coding.rs`, `stop_codon_at_the_coding_end` and
//! `start_codon_stands_at_the_coding_start`); Ensembl VEP reads the last and the first three
//! coding positions as those codons on every transcript without `cds_end_NF` or
//! `cds_start_NF` (VariationEffect.pm:1373-1381, 975-986). `coding_ends` is that reading of
//! vep-rs, made step for step from the facts of one record, and
//! `tests/fixtures/coding_ends/verdicts.tsv` holds both readings to each other on the
//! transcripts of committed corpus caches.
//!
//! The file names no type of a comparison tool, so that each tool root can compile this one
//! source: the facts of a record come in as plain values (`RecordFacts`) and the reading
//! goes out as plain values (`CodingEnds`). Its module `terms` (`coding_end_terms.rs`,
//! beside this file) derives from that reading and a row the coding terms vep-rs owes, and
//! asks the root's module `coding_end` (`rules/coding_end.rs`, which reads
//! `genetic_codes.rs` as the root's `genetic_codes`) for the shapes of one class.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::LazyLock;

#[path = "coding_end_terms.rs"]
pub mod terms;

/// The stop codon at the coding end of a transcript.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopEnd {
    /// The coding sequence is a whole number of codons and its last codon is a stop codon
    /// of the transcript's genetic code.
    Whole,
    /// The coding sequence ends `bases` bases, one or two, into a codon on the transcript's
    /// last base, those bases followed by adenines are a stop codon of the code, and the
    /// record shows the model complete, by a genetic code of
    /// `CODES_THAT_SHOW_A_COMPLETED_STOP_CODON` or by a consensus coding sequence id:
    /// polyadenylation completes the stop codon.
    BegunByTheTail { bases: u8 },
    /// No stop codon: the transcript carries `cds_end_NF` (`flagged`), its last whole codon
    /// names a residue, or its coding sequence ends inside a codon that a 3' UTR follows,
    /// that adenines do not complete to a stop codon, or that no record shows completed.
    NoneShown { flagged: bool },
    /// The record is not read: it holds no coding sequence, or states a part of itself in
    /// a form vep-rs reads otherwise or refuses.
    Unread,
}

/// The start codon at the coding start of a transcript.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StartEnd {
    /// The transcript carries no `cds_start_NF`, its first codon is whole, and its
    /// translation opens with methionine or its record carries a consensus coding sequence
    /// id.
    Shown,
    /// No start codon: the transcript carries `cds_start_NF` (`flagged`), the start exon's
    /// phase pads its first codon (`padded`), or neither its translation nor its record
    /// shows one.
    NotShown { flagged: bool, padded: bool },
    /// The record is not read, as `StopEnd::Unread`.
    Unread,
}

impl fmt::Display for StopEnd {
    fn fmt(&self, text: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            StopEnd::Whole => text.write_str("whole"),
            StopEnd::BegunByTheTail { bases } => write!(
                text,
                "begun by the tail, {bases} {}, shown",
                if bases == 1 { "base" } else { "bases" }
            ),
            StopEnd::NoneShown { flagged: true } => text.write_str("none shown, flagged"),
            StopEnd::NoneShown { flagged: false } => text.write_str("none shown"),
            StopEnd::Unread => text.write_str("unread"),
        }
    }
}

impl fmt::Display for StartEnd {
    fn fmt(&self, text: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            StartEnd::Shown => text.write_str("shown"),
            StartEnd::NotShown { flagged, padded } => {
                text.write_str("not shown")?;
                if flagged {
                    text.write_str(", flagged")?;
                }
                if padded {
                    text.write_str(", padded")?;
                }
                Ok(())
            }
            StartEnd::Unread => text.write_str("unread"),
        }
    }
}

/// The genome's bases on either side of a span of a model that has no genome base, a block
/// of bases an edit of the model inserts: `before` is the base whose cDNA position is the
/// last before the span and `after` the one whose cDNA position is the first after it,
/// `None` where the model holds no genome base on that side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EditFlanks {
    pub before: Option<i64>,
    pub after: Option<i64>,
}

/// The coordinates of a transcript one of whose coding ends shows no codon and carries no
/// flag: what tells an insertion between two neighbouring bases of the model from one
/// across bases the model holds and the genome lacks, and a span that holds a coding base
/// from one that holds none.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CodingModel {
    /// The mapper's `cdna_coding_start`.
    pub coding_start: i64,
    /// The mapper's `cdna_coding_end`.
    pub coding_end: i64,
    /// The count of bases the start exon's phase pads the coding sequence with, which
    /// have no cDNA position.
    pub padding: i64,
    /// The mapper's pairs, each a run of cDNA coordinates, the genomic bases it maps to
    /// and their orientation: (`from_start`, `from_end`, `to_start`, `to_end`, `ori`).
    pub pairs: Vec<(i64, i64, i64, i64, i64)>,
    /// The genomic positions of the first and the last coding base a pair holds, in
    /// genomic order.
    pub coding_region: Option<(i64, i64)>,
    /// The record's introns whose two ends lie at most twelve bases apart, each as its
    /// first and its last genomic base. Ensembl VEP reads a span that meets such an
    /// intron and the coding region as inside the coding sequence
    /// (VariationEffect.pm:1616-1639 on BaseTranscriptVariation.pm's `_introns`).
    pub short_introns: Vec<(i64, i64)>,
}

/// The sequence of a transcript whose record shows that polyadenylation completes its
/// stop codon (`StopEnd::BegunByTheTail`), on which vep-rs reads every term
/// and cell of a sequence allele's row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompletedStopCodon {
    /// The cached coding sequence in upper case, without the adenines that complete its
    /// last codon.
    pub coding: Vec<u8>,
    /// The transcript's peptide as vep-rs reads it (`crates/vep-effects/src/coding.rs`,
    /// `transcript_peptide`, on the completed coding sequence): the cached one with its
    /// first residue read as methionine where the first codon is a start codon of the
    /// code, or without a cached one the translation of the whole codons less the closing
    /// stop codon.
    pub peptide: Vec<u8>,
    /// The id of the transcript's genetic code.
    pub table: i64,
}

/// The sequence of a transcript at whose coding start no start codon stands, with
/// `cds_start_NF` or without, on which vep-rs writes the protein notation of a frameshift
/// at its first codon to the stop codon of the new frame. Ensembl VEP writes a lost start
/// there where it reads one, on a transcript without the flag, and elsewhere the same
/// frameshift with the unknown count where residue 1 is the first that changes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnshownStartSequence {
    /// The cached coding sequence in upper case, the padding letters of the start exon's
    /// phase included.
    pub coding: Vec<u8>,
    /// The cached 3' UTR in upper case, empty for a transcript that ends at its coding
    /// end.
    pub three_prime_utr: Vec<u8>,
    /// The transcript's peptide as vep-rs reads it (`CompletedStopCodon::peptide`, on the
    /// cached coding sequence).
    pub peptide: Vec<u8>,
    /// The id of the transcript's genetic code.
    pub table: i64,
}

/// The sequence of a transcript at whose coding end no stop codon stands and that carries
/// no `cds_end_NF`, on which vep-rs writes the protein notation of a frameshift that the
/// transcript's end cuts short where Ensembl VEP writes it beside a lost stop.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnclosedEndSequence {
    /// The cached coding sequence in upper case.
    pub coding: Vec<u8>,
    /// The cached 3' UTR in upper case, empty for a transcript that ends at its coding
    /// end.
    pub three_prime_utr: Vec<u8>,
    /// The transcript's peptide as vep-rs reads it (`CompletedStopCodon::peptide`, on the
    /// cached coding sequence).
    pub peptide: Vec<u8>,
    /// The id of the transcript's genetic code.
    pub table: i64,
}

/// The two ends of one transcript's coding sequence. A coding position counts from the
/// first base of the cached coding sequence, the padding of the start exon's phase
/// included, and L is the sequence's length; a position the genome holds no base for is
/// `None`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CodingEnds {
    pub stop: StopEnd,
    pub start: StartEnd,
    /// The count of coding bases after the last whole codon: L modulo three, and 0 for a
    /// record without a coding sequence.
    pub tail_length: u8,
    /// The transcript's strand, 1 or -1.
    pub strand: Option<i64>,
    /// The genomic positions of coding positions 1, 2 and 3.
    pub first_codon: [Option<i64>; 3],
    /// Where coding positions 1, 2 and 3 all lie in one span of the model that has no
    /// genome base, the genome's bases on either side of that span.
    pub first_codon_edit: Option<EditFlanks>,
    /// The genomic positions of coding positions L-2, L-1 and L.
    pub last_codon: [Option<i64>; 3],
    /// The lowest and the highest of the last three genomic positions of the coding
    /// region, which lie beside each other whatever intron splits the last codon
    /// (VariationEffect.pm:1418-1427).
    pub coding_region_end: Option<(i64, i64)>,
    /// The genomic positions of the tail's bases, the last `tail_length` coding positions,
    /// in coding order.
    pub tail: [Option<i64>; 2],
    /// The tail's bases as the cache holds them, in upper case.
    pub tail_bases: [Option<u8>; 2],
    /// L, the length of the cached coding sequence; 0 for a record that is not read.
    pub coding_length: i64,
    /// The record's own first and last genomic base (`start`, `end`), by which a row is
    /// answered with the record whose span holds the row's Location.
    pub span: Option<(i64, i64)>,
    /// Whether two of the mapper's pairs that follow each other on the cDNA leave a cDNA
    /// position between them: the model holds bases there that the genome lacks, and an
    /// insertion between the genome bases on either side stands on those bases.
    pub edited_model: bool,
    /// The transcript's coordinates, kept where the stop end reads `NoneShown` without
    /// the flag or `BegunByTheTail`, or the start end `NotShown` without the flag.
    pub model: Option<Box<CodingModel>>,
    /// The transcript's sequence, kept where the stop end reads `BegunByTheTail` under a
    /// genetic code the module holds a table for.
    pub completed: Option<Box<CompletedStopCodon>>,
    /// The transcript's sequence, kept where the start end reads `NotShown`, with the
    /// flag or without, no polyadenylation completes the stop codon, the genetic code is
    /// one the module holds a table for, and the reader handed the record's 3' UTR or the
    /// transcript has none (`coding_ends_with_three_prime_utr`).
    pub unshown_start: Option<Box<UnshownStartSequence>>,
    /// The transcript's sequence, kept where the stop end reads `NoneShown` without the
    /// flag, the start end is read, the genetic code is one the module holds a table for,
    /// and the reader handed the record's 3' UTR or the transcript has none
    /// (`coding_ends_with_three_prime_utr`).
    pub unclosed_end: Option<Box<UnclosedEndSequence>>,
}

impl CodingEnds {
    /// The ends of a record that is not read, whose coding sequence has `tail_length` bases
    /// after its last whole codon.
    pub fn unread(tail_length: u8) -> CodingEnds {
        CodingEnds {
            stop: StopEnd::Unread,
            start: StartEnd::Unread,
            tail_length,
            strand: None,
            first_codon: [None; 3],
            first_codon_edit: None,
            last_codon: [None; 3],
            coding_region_end: None,
            tail: [None; 2],
            tail_bases: [None; 2],
            coding_length: 0,
            span: None,
            edited_model: false,
            model: None,
            completed: None,
            unshown_start: None,
            unclosed_end: None,
        }
    }
}

/// A cache transcript's `variation_effect_feature_cache.mapper` as vep-rs loads it.
#[derive(Debug, Clone, Copy)]
pub struct MapperFacts<'a> {
    /// `start_phase`, 0 where the mapper states none.
    pub start_phase: i64,
    /// `cdna_coding_start`, 0 where the mapper states none.
    pub coding_start: i64,
    /// `cdna_coding_end`, 0 where the mapper states none.
    pub coding_end: i64,
    /// `pairs`, each a run of cDNA coordinates, the genomic bases it maps to and their
    /// orientation: (`from_start`, `from_end`, `to_start`, `to_end`, `ori`).
    pub pairs: &'a [(i64, i64, i64, i64, i64)],
}

/// The facts of one cache transcript object that the reading asks, each as the record
/// states it; `None` for a member the record lacks or holds as `null`.
#[derive(Debug, Clone, Copy, Default)]
pub struct RecordFacts<'a> {
    /// The `code` of each item of `attributes`, whatever the item's value.
    pub attribute_codes: &'a [&'a str],
    /// `ccds`.
    pub consensus_id: Option<&'a str>,
    /// `strand`, 1 or -1.
    pub strand: Option<i64>,
    /// `variation_effect_feature_cache.codon_table`.
    pub codon_table: Option<i64>,
    /// `variation_effect_feature_cache.translateable_seq`.
    pub coding: Option<&'a str>,
    /// `variation_effect_feature_cache.peptide`.
    pub peptide: Option<&'a str>,
    /// `variation_effect_feature_cache.mapper`.
    pub mapper: Option<MapperFacts<'a>>,
    /// The record's `start` and `end`, its first and last genomic base.
    pub span: Option<(i64, i64)>,
    /// The record's introns as vep-rs loads them, each as its first and its last genomic
    /// base: those of `variation_effect_feature_cache.introns`, or where that holds none
    /// the stretches between the record's `exons`.
    pub introns: &'a [(i64, i64)],
}

/// The genetic codes under which a reference sequence record annotates a coding sequence
/// that ends inside its stop codon as completed by polyadenylation: the vertebrate,
/// invertebrate, echinoderm and ascidian mitochondrial codes.
pub const CODES_THAT_SHOW_A_COMPLETED_STOP_CODON: [i64; 4] = [2, 5, 9, 13];

const STANDARD_CODE: i64 = 1;

/// The nuclear code `Bio::Tools::CodonTable` holds no table for and vep-rs reads
/// (`crates/vep-core/src/codon.rs`): the standard code but for one codon of leucine, which
/// it reads as alanine, with that codon and the standard code's codon of methionine as its
/// start codons. The two codons are named by their places in a code's 64.
const CODE_THE_MODULE_LACKS: i64 = 26;
const ALANINE_OF_THAT_CODE: usize = 19;
const METHIONINE: usize = 35;

/// The 64 codons under every genetic code of `Bio::Tools::CodonTable`, under `#` header
/// lines one line a codon: the table id, the codon, its amino acid, its start mark and its
/// stop mark, tab-separated.
const MODULE_CODES_TABLE: &str =
    include_str!("../../../crates/vep-core/tests/data/genetic_codes.tsv");

/// One genetic code: the amino acid of each codon, a stop `*`, and whether the codon is a
/// start codon; a codon's place counts its first base as the most significant, the bases
/// `T`, `C`, `A` and `G` as 0 to 3.
#[derive(Clone)]
struct GeneticCode {
    amino_acids: [u8; 64],
    starts: [bool; 64],
}

/// The place of a base that is one of the four, in either case, `U` read as `T`.
fn rank(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// The bases a nucleotide letter stands for, as places, in either case; none for a byte
/// that is no nucleotide letter.
fn readings(base: u8) -> &'static [usize] {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => &[0],
        b'C' => &[1],
        b'A' => &[2],
        b'G' => &[3],
        b'Y' => &[0, 1],
        b'W' => &[0, 2],
        b'K' => &[0, 3],
        b'M' => &[1, 2],
        b'S' => &[1, 3],
        b'R' => &[2, 3],
        b'H' => &[0, 1, 2],
        b'B' => &[0, 1, 3],
        b'D' => &[0, 2, 3],
        b'V' => &[1, 2, 3],
        b'N' | b'X' => &[0, 1, 2, 3],
        _ => &[],
    }
}

impl GeneticCode {
    /// `Bio::Tools::CodonTable::_codon_is`, as vep-rs reads it: the codon is three
    /// nucleotide letters and every codon of the four bases it stands for satisfies `marked`.
    fn every_reading(&self, codon: &[u8], marked: impl Fn(usize) -> bool) -> bool {
        let &[first, second, third] = codon else {
            return false;
        };
        let (firsts, seconds, thirds) = (readings(first), readings(second), readings(third));
        if firsts.is_empty() || seconds.is_empty() || thirds.is_empty() {
            return false;
        }
        firsts.iter().all(|a| {
            seconds
                .iter()
                .all(|b| thirds.iter().all(|c| marked((a << 4) | (b << 2) | c)))
        })
    }

    fn is_start(&self, codon: &[u8]) -> bool {
        self.every_reading(codon, |place| self.starts[place])
    }

    fn is_stop(&self, codon: &[u8]) -> bool {
        self.every_reading(codon, |place| self.amino_acids[place] == b'*')
    }

    /// `Bio::Tools::CodonTable::translate` of one codon: a codon of the four bases by the
    /// code, three gaps as a gap, and a codon with an ambiguous letter by the amino acids
    /// its readings code: one, that one; aspartate and asparagine `B`; glutamate and
    /// glutamine `Z`; anything else `X`.
    fn residue(&self, codon: &[u8]) -> u8 {
        let &[first, second, third] = codon else {
            return b'X';
        };
        if let (Some(a), Some(b), Some(c)) = (rank(first), rank(second), rank(third)) {
            return self.amino_acids[(a << 4) | (b << 2) | c];
        }
        if codon.iter().all(|base| *base == b'-') {
            return b'-';
        }
        let mut distinct: Vec<u8> = Vec::new();
        for a in readings(first) {
            for b in readings(second) {
                for c in readings(third) {
                    let amino_acid = self.amino_acids[(a << 4) | (b << 2) | c];
                    if !distinct.contains(&amino_acid) {
                        distinct.push(amino_acid);
                    }
                }
            }
        }
        match distinct.as_slice() {
            [amino_acid] => *amino_acid,
            [b'D', b'N'] | [b'N', b'D'] => b'B',
            [b'E', b'Q'] | [b'Q', b'E'] => b'Z',
            _ => b'X',
        }
    }

    /// The residue `Bio::EnsEMBL::Transcript::translate` writes for a peptide's first
    /// codon: methionine for a start codon of the code, else the codon's own.
    fn first_residue(&self, codon: &[u8]) -> u8 {
        if self.is_start(codon) {
            b'M'
        } else {
            self.residue(codon)
        }
    }
}

/// The genetic codes vep-rs translates with: the module's by table id, read from
/// `MODULE_CODES_TABLE`, and `CODE_THE_MODULE_LACKS`. Id 0, which the module answers as
/// the standard code, is left out. A line that is no codon of a table, or a table short of
/// a codon, is a panic, since the text is part of the build.
static CODES: LazyLock<BTreeMap<i64, GeneticCode>> = LazyLock::new(|| {
    let mut read: BTreeMap<i64, ([u8; 64], [bool; 64])> = BTreeMap::new();
    for line in MODULE_CODES_TABLE.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let place = cells.get(1).and_then(|codon| match codon.as_bytes() {
            &[first, second, third] => {
                Some((rank(first)? << 4) | (rank(second)? << 2) | rank(third)?)
            }
            _ => None,
        });
        let (Some(id), Some(place), Some(&[amino_acid]), Some(start)) = (
            cells.first().and_then(|id| id.parse::<i64>().ok()),
            place,
            cells.get(2).map(|residue| residue.as_bytes()),
            cells.get(3),
        ) else {
            panic!("the genetic-code table holds a line that names no codon of a table: {line}")
        };
        let (amino_acids, starts) = read.entry(id).or_insert(([0; 64], [false; 64]));
        amino_acids[place] = amino_acid;
        starts[place] = *start == "1";
    }
    read.remove(&0);
    let mut codes: BTreeMap<i64, GeneticCode> = BTreeMap::new();
    for (id, (amino_acids, starts)) in read {
        assert!(
            amino_acids.iter().all(|amino_acid| *amino_acid != 0),
            "the genetic-code table holds fewer than 64 codons under id {id}"
        );
        codes.insert(
            id,
            GeneticCode {
                amino_acids,
                starts,
            },
        );
    }
    let mut lacking = codes[&STANDARD_CODE].clone();
    lacking.amino_acids[ALANINE_OF_THAT_CODE] = b'A';
    lacking.starts = [false; 64];
    lacking.starts[ALANINE_OF_THAT_CODE] = true;
    lacking.starts[METHIONINE] = true;
    codes.insert(CODE_THE_MODULE_LACKS, lacking);
    codes
});

/// The table `Bio::Tools::CodonTable` reads the held code `table` under: itself, or the
/// module's last table for the code the module lacks (CodonTable.pm:328-338, :412-460). A
/// dumper on that module writes a transcript's cached peptide under the table returned.
fn module_table(table: i64) -> i64 {
    if table == CODE_THE_MODULE_LACKS {
        CODES
            .keys()
            .copied()
            .filter(|id| *id != CODE_THE_MODULE_LACKS)
            .max()
            .unwrap_or(table)
    } else {
        table
    }
}

/// Whether the transcript's peptide opens with methionine, as vep-rs reads the peptide
/// (`crates/vep-effects/src/coding.rs`, `transcript_peptide`). A cached peptide is first
/// read again under the transcript's own code where the module wrote it under another
/// (`crates/vep-cli/src/json_cache.rs`, `peptide_under_own_code`) and then opens with
/// methionine where its first residue is one or the first codon is a start codon of the
/// code. Without a cached peptide the translation of the whole codons, less a closing stop
/// codon, opens with methionine where its first codon names it or is a start codon.
fn opens_with_methionine(
    peptide: Option<&[u8]>,
    coding: &[u8],
    table: i64,
    code: &GeneticCode,
) -> bool {
    let first_codon = coding.get(..3);
    let Some(peptide) = peptide else {
        let whole = &coding[..coding.len() - coding.len() % 3];
        let Some(last_codon) = whole.len().checked_sub(3) else {
            return false;
        };
        let residues = whole.len() / 3 - usize::from(code.is_stop(&whole[last_codon..]));
        return residues > 0
            && first_codon
                .is_some_and(|codon| code.residue(codon) == b'M' || code.is_start(codon));
    };
    let Some(&cached) = peptide.first() else {
        return false;
    };
    let dumped_under = module_table(table);
    let first = match (first_codon, CODES.get(&dumped_under)) {
        (Some(codon), Some(module_code))
            if dumped_under != table && cached == module_code.first_residue(codon) =>
        {
            code.first_residue(codon)
        }
        _ => cached,
    };
    first == b'M' || first_codon.is_some_and(|codon| code.is_start(codon))
}

/// The transcript's peptide as vep-rs reads it on the coding sequence `completed`, which
/// a stop codon closes (`crates/vep-effects/src/coding.rs`, `transcript_peptide`): the
/// cached peptide, or without one the translation of the whole codons less the closing
/// stop codon, each with its first residue read as methionine where the first codon is a
/// start codon of the code.
fn peptide_on_the_completed_sequence(
    cached: Option<&[u8]>,
    completed: &[u8],
    code: &GeneticCode,
) -> Vec<u8> {
    let whole = &completed[..completed.len() - completed.len() % 3];
    let mut peptide = match cached {
        Some(cached) => cached.to_vec(),
        None => {
            let mut translated: Vec<u8> =
                whole.chunks(3).map(|codon| code.residue(codon)).collect();
            if whole.len() >= 3 && code.is_stop(&whole[whole.len() - 3..]) {
                translated.pop();
            }
            translated
        }
    };
    if let (Some(first), Some(codon)) = (peptide.first_mut(), completed.get(..3)) {
        if *first != b'M' && code.is_start(codon) {
            *first = b'M';
        }
    }
    peptide
}

fn narrow(value: i128) -> Option<i64> {
    i64::try_from(value).ok()
}

impl MapperFacts<'_> {
    /// The genomic base of the cDNA position `cdna` inside the pair `pair`.
    fn image(pair: &(i64, i64, i64, i64, i64), cdna: i64, reverse: bool) -> Option<i64> {
        let &(from_start, _, to_start, to_end, _) = pair;
        let along = i128::from(cdna) - i128::from(from_start);
        narrow(if reverse {
            i128::from(to_end) - along
        } else {
            i128::from(to_start) + along
        })
    }

    /// The genomic base of the coding position `cds`, counted with `padding` bases before
    /// the first coding base of the cDNA; `None` for a position outside the coding bounds
    /// or the pairs.
    fn genomic(&self, padding: i64, cds: i128) -> Option<i64> {
        let cdna = narrow(cds - i128::from(padding) + i128::from(self.coding_start) - 1)?;
        if !(self.coding_start <= cdna && cdna <= self.coding_end) {
            return None;
        }
        self.pairs
            .iter()
            .find(|pair| pair.0 <= cdna && cdna <= pair.1)
            .and_then(|pair| Self::image(pair, cdna, pair.4 == -1))
    }

    /// The genomic positions of the first and the last coding base a pair holds, in
    /// genomic order (`crates/vep-cli/src/json_cache.rs`, `compute_coding_region`): a
    /// coding bound in bases the genome lacks is the nearest coding base the genome has.
    fn coding_region(&self) -> Option<(i64, i64)> {
        if self.coding_start <= 0 || self.coding_end <= 0 {
            return None;
        }
        let (first_cdna, last_cdna) = (
            self.coding_start.min(self.coding_end),
            self.coding_start.max(self.coding_end),
        );
        let coding = || {
            self.pairs
                .iter()
                .filter(|pair| pair.0 <= pair.1 && pair.1 >= first_cdna && pair.0 <= last_cdna)
        };
        let first = coding().min_by_key(|pair| pair.0)?;
        let last = coding().max_by_key(|pair| pair.1)?;
        let start = Self::image(first, first_cdna.max(first.0), first.4 < 0)?;
        let end = Self::image(last, last_cdna.min(last.1), last.4 < 0)?;
        Some((start.min(end), start.max(end)))
    }

    /// Whether two pairs that follow each other on the cDNA leave a cDNA position
    /// between them.
    fn skips_a_cdna_position(&self) -> bool {
        let mut runs: Vec<(i64, i64)> = self.pairs.iter().map(|pair| (pair.0, pair.1)).collect();
        runs.sort_unstable();
        runs.windows(2)
            .any(|runs| i128::from(runs[1].0) > i128::from(runs[0].1) + 1)
    }

    /// The genome's bases on either side of the cDNA positions `first` to `last`, where no
    /// pair holds one of them.
    fn flanks(&self, first: i64, last: i64) -> Option<EditFlanks> {
        if self
            .pairs
            .iter()
            .any(|pair| pair.0 <= last && first <= pair.1)
        {
            return None;
        }
        let before = self
            .pairs
            .iter()
            .filter(|pair| pair.1 < first)
            .max_by_key(|pair| pair.1);
        let after = self
            .pairs
            .iter()
            .filter(|pair| pair.0 > last)
            .min_by_key(|pair| pair.0);
        Some(EditFlanks {
            before: before.and_then(|pair| Self::image(pair, pair.1, pair.4 == -1)),
            after: after.and_then(|pair| Self::image(pair, pair.0, pair.4 == -1)),
        })
    }
}

/// The coding ends of the transcript whose cache record states `record`, each step as
/// vep-rs takes it on the loaded record.
///
/// The stop end (`stop_codon_at_the_coding_end`): a `cds_end_NF` attribute, whatever its
/// value, shows none; whole codons are closed where the last is a stop codon of the code;
/// a tail of one or two bases begins a stop codon where the mapper's largest cDNA
/// coordinate is a coding end above 0, the tail followed by adenines is a stop codon,
/// and the record shows it completed, under a code of
/// `CODES_THAT_SHOW_A_COMPLETED_STOP_CODON` or by a consensus id that is neither empty nor
/// a dash; a tail no record shows completed begins none. The start end (`start_codon_stands_at_the_coding_start`): no `cds_start_NF`
/// attribute, a start phase that is not above 0, and a peptide that opens with methionine
/// or such a consensus id. A `codon_table` the record lacks, or of 0, is the standard code.
///
/// Both ends are `Unread` for a record whose coding sequence is absent, empty or not
/// ASCII, and for a genetic code no table is held for, which vep-rs refuses to load; the
/// start end is `Unread` for a cached peptide that is not ASCII. A record that is not
/// read keeps its span.
///
/// Beside the two ends the reading keeps what a derivation of the row's terms asks of
/// the record (`terms::expected_terms`): the coding length and the record's span of every
/// transcript, the coordinates of one with an end that shows no codon and carries no flag
/// (`CodingModel`), and the sequence of one whose stop codon polyadenylation completes
/// (`CompletedStopCodon`).
pub fn coding_ends(record: &RecordFacts<'_>) -> CodingEnds {
    coding_ends_with_three_prime_utr(record, None)
}

/// `coding_ends` of a record whose cached 3' UTR sequence
/// (`variation_effect_feature_cache.three_prime_utr`) the reader hands over, `None` for
/// a record that caches none. With it the reading keeps the sequence of a transcript
/// at whose coding start no start codon stands, flagged `cds_start_NF` or not
/// (`UnshownStartSequence`), and of one at whose coding end no stop codon stands and that
/// carries no `cds_end_NF` (`UnclosedEndSequence`). vep-rs reads the cached 3' UTR where the
/// record holds one that is not empty; a transcript whose largest cDNA coordinate is its
/// coding end has none; any other record's 3' UTR vep-rs reads from the reference
/// sequence, and its sequence is not kept.
pub fn coding_ends_with_three_prime_utr(
    record: &RecordFacts<'_>,
    three_prime_utr: Option<&str>,
) -> CodingEnds {
    let mut ends = read_coding_ends(record);
    let table = match record.codon_table {
        None | Some(0) => STANDARD_CODE,
        Some(table) => table,
    };
    let unshown_start = matches!(
        (ends.start, ends.stop),
        (
            StartEnd::NotShown { .. },
            StopEnd::Whole | StopEnd::NoneShown { .. }
        )
    );
    let unclosed_end =
        ends.stop == (StopEnd::NoneShown { flagged: false }) && ends.start != StartEnd::Unread;
    let (true, Some(coding), Some(code), true) = (
        unshown_start || unclosed_end,
        record.coding,
        CODES.get(&table),
        module_table(table) == table,
    ) else {
        return ends;
    };
    let ends_at_the_coding_end = record.mapper.as_ref().is_some_and(|mapper| {
        mapper.coding_end > 0
            && mapper.pairs.iter().map(|pair| pair.1).max() == Some(mapper.coding_end)
    });
    let three_prime_utr = match three_prime_utr.filter(|bases| !bases.is_empty()) {
        Some(bases) if bases.is_ascii() => bases.to_ascii_uppercase().into_bytes(),
        None if ends_at_the_coding_end => Vec::new(),
        _ => return ends,
    };
    let coding = coding.to_ascii_uppercase().into_bytes();
    let peptide =
        peptide_on_the_completed_sequence(record.peptide.map(str::as_bytes), &coding, code);
    if unclosed_end {
        ends.unclosed_end = Some(Box::new(UnclosedEndSequence {
            coding: coding.clone(),
            three_prime_utr: three_prime_utr.clone(),
            peptide: peptide.clone(),
            table,
        }));
    }
    if unshown_start {
        ends.unshown_start = Some(Box::new(UnshownStartSequence {
            coding,
            three_prime_utr,
            peptide,
            table,
        }));
    }
    ends
}

fn read_coding_ends(record: &RecordFacts<'_>) -> CodingEnds {
    let unread = |tail_length: u8| CodingEnds {
        span: record.span,
        ..CodingEnds::unread(tail_length)
    };
    let Some(coding) = record.coding.filter(|coding| !coding.is_empty()) else {
        return unread(0);
    };
    let tail_length = (coding.chars().count() % 3) as u8;
    let table = match record.codon_table {
        None | Some(0) => STANDARD_CODE,
        Some(table) => table,
    };
    let (true, Some(code)) = (coding.is_ascii(), CODES.get(&table)) else {
        return unread(tail_length);
    };
    let coding = coding.as_bytes();
    let tail = usize::from(tail_length);
    let carries = |code: &str| record.attribute_codes.contains(&code);
    let consensus = record
        .consensus_id
        .is_some_and(|id| !id.is_empty() && id != "-");
    let mapper = record.mapper.as_ref();

    let stop = if carries("cds_end_NF") {
        StopEnd::NoneShown { flagged: true }
    } else if tail == 0 {
        match code.is_stop(&coding[coding.len() - 3..]) {
            true => StopEnd::Whole,
            false => StopEnd::NoneShown { flagged: false },
        }
    } else {
        let ends_at_the_coding_end = mapper.is_some_and(|mapper| {
            mapper.coding_end > 0
                && mapper.pairs.iter().map(|pair| pair.1).max() == Some(mapper.coding_end)
        });
        let mut completed = [b'A'; 3];
        completed[..tail].copy_from_slice(&coding[coding.len() - tail..]);
        let shown_complete = CODES_THAT_SHOW_A_COMPLETED_STOP_CODON.contains(&table) || consensus;
        if ends_at_the_coding_end && code.is_stop(&completed) && shown_complete {
            StopEnd::BegunByTheTail { bases: tail_length }
        } else {
            StopEnd::NoneShown { flagged: false }
        }
    };

    let padded = mapper.is_some_and(|mapper| mapper.start_phase > 0);
    let flagged = carries("cds_start_NF");
    let start = match record.peptide {
        Some(peptide) if !peptide.is_ascii() => StartEnd::Unread,
        peptide => {
            let opens = opens_with_methionine(peptide.map(str::as_bytes), coding, table, code);
            if !padded && (opens || consensus) && !flagged {
                StartEnd::Shown
            } else {
                StartEnd::NotShown { flagged, padded }
            }
        }
    };

    let padding = mapper.map_or(0, |mapper| mapper.start_phase.max(0));
    let length = coding.len() as i128;
    let place = |cds: i128| mapper.and_then(|mapper| mapper.genomic(padding, cds));
    let tail_place = |index: usize| match index < tail {
        true => place(length - tail as i128 + 1 + index as i128),
        false => None,
    };
    let tail_base = |index: usize| match index < tail {
        true => Some(coding[coding.len() - tail + index].to_ascii_uppercase()),
        false => None,
    };
    CodingEnds {
        stop,
        start,
        tail_length,
        strand: record.strand,
        first_codon: [place(1), place(2), place(3)],
        first_codon_edit: mapper
            .filter(|mapper| padding == 0 && mapper.coding_start > 0)
            .and_then(|mapper| {
                mapper.flanks(mapper.coding_start, mapper.coding_start.checked_add(2)?)
            }),
        last_codon: [place(length - 2), place(length - 1), place(length)],
        coding_region_end: mapper.and_then(MapperFacts::coding_region).and_then(
            |(lowest, highest)| match record.strand {
                Some(1) => Some((highest.checked_sub(2)?, highest)),
                Some(-1) => Some((lowest, lowest.checked_add(2)?)),
                _ => None,
            },
        ),
        tail: [tail_place(0), tail_place(1)],
        coding_length: coding.len() as i64,
        span: record.span,
        edited_model: mapper.is_some_and(MapperFacts::skips_a_cdna_position),
        model: mapper
            .filter(|_| {
                matches!(
                    stop,
                    StopEnd::NoneShown { flagged: false } | StopEnd::BegunByTheTail { .. }
                ) || matches!(start, StartEnd::NotShown { flagged: false, .. })
            })
            .map(|mapper| {
                Box::new(CodingModel {
                    coding_start: mapper.coding_start,
                    coding_end: mapper.coding_end,
                    padding,
                    pairs: mapper.pairs.to_vec(),
                    coding_region: mapper.coding_region(),
                    short_introns: record
                        .introns
                        .iter()
                        .copied()
                        .filter(|(first, last)| {
                            (i128::from(*last) - i128::from(*first)).abs() <= 12
                        })
                        .collect(),
                })
            }),
        completed: match (stop, start) {
            (StopEnd::BegunByTheTail { .. }, StartEnd::Shown | StartEnd::NotShown { .. })
                if module_table(table) == table =>
            {
                let mut completed = coding.to_ascii_uppercase();
                completed.resize(coding.len() + 3 - tail, b'A');
                let peptide = peptide_on_the_completed_sequence(
                    record.peptide.map(str::as_bytes),
                    &completed,
                    code,
                );
                Some(Box::new(CompletedStopCodon {
                    coding: coding.to_ascii_uppercase(),
                    peptide,
                    table,
                }))
            }
            _ => None,
        },
        unshown_start: None,
        unclosed_end: None,
        tail_bases: [tail_base(0), tail_base(1)],
    }
}
