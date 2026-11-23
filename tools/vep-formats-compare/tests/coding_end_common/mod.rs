// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the classes of a coding end share: a planted cache of made-up
//! transcripts and the rows of made-up records on them, written in each format.
//!
//! Every transcript lies on chromosome 7 with its first cDNA base at the genomic base
//! `FIRST`, on the forward strand unless its id ends in `_REVERSE`, where its last cDNA base
//! lies there, so one record has a row on each: six bases of 5' UTR unless the transcript
//! states another count, the coding sequence, then the 3' UTR. A row's cells are Ensembl
//! VEP's on the record, written by hand.

#![allow(dead_code)]

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use vep_formats_compare::cli;
use vep_formats_compare::consts::CONSEQUENCE_RANK;
use vep_formats_compare::pyjson::{self, Value};

pub const CHROMOSOME: &str = "7";
pub const FIRST: i64 = 1001;
pub const LEADER: i64 = 6;

/// Codons of methionine, alanine, glycine and glutamate, and a stop codon of the standard
/// and the vertebrate mitochondrial code.
pub const MET: &str = "ATG";
pub const ALA: &str = "GCT";
pub const GLY: &str = "GGA";
pub const GLU: &str = "GAA";
pub const STOP: &str = "TAA";

pub const STRUCTURAL: &str = "structural_deletion_stop_lost_without_stop_codon";
pub const COMPLETED: &str = "stop_codon_completed_by_polyadenylation";
pub const OVERLAP: &str = "stop_codon_overlap_without_stop_codon";
pub const BEHIND: &str = "stop_codon_behind_last_residue_without_stop_codon";
pub const START: &str = "start_codon_overlap_without_start_codon";
pub const FIVE: [&str; 5] = [STRUCTURAL, COMPLETED, OVERLAP, BEHIND, START];

/// One made-up transcript of the planted cache.
#[derive(Clone)]
pub struct Transcript {
    pub id: &'static str,
    pub coding: String,
    pub trailer: i64,
    pub flags: Vec<&'static str>,
    pub consensus_id: Option<&'static str>,
    pub codon_table: Option<i64>,
    pub peptide: Option<&'static str>,
    /// The 3' UTR the cache holds for the transcript, where it holds one.
    pub three_prime_utr: Option<&'static str>,
    /// The genomic base of the first cDNA base.
    pub first: i64,
    /// The cDNA position after which an intron of a hundred bases lies, 0 for none.
    pub intron_after: i64,
    /// The count of bases of 5' UTR.
    pub leader: i64,
    /// The count of bases the start exon's phase pads the coding sequence with: the
    /// first `start_phase` letters of `coding` have no cDNA position.
    pub start_phase: i64,
    /// The strand, 1 or -1; a transcript of the reverse strand has no intron.
    pub strand: i64,
}

fn transcript(id: &'static str, codons: &[&str], trailer: i64) -> Transcript {
    Transcript {
        id,
        coding: codons.concat(),
        trailer,
        flags: Vec::new(),
        consensus_id: None,
        codon_table: None,
        peptide: None,
        three_prime_utr: None,
        first: FIRST,
        intron_after: 0,
        leader: LEADER,
        start_phase: 0,
        strand: 1,
    }
}

/// Whole codons whose last names a residue, before a 3' UTR: no stop codon stands.
pub fn open() -> Transcript {
    transcript("ENST_OPEN", &[MET, ALA, GLY], 6)
}

/// Whole codons closed by a stop codon, before a 3' UTR: a stop codon and a start codon
/// stand.
pub fn closed() -> Transcript {
    transcript("ENST_CLOSED", &[MET, ALA, GLY, STOP], 6)
}

/// `open` with `cds_end_NF`.
pub fn flagged() -> Transcript {
    Transcript {
        id: "ENST_FLAGGED",
        flags: vec!["cds_end_NF"],
        ..open()
    }
}

/// A coding sequence that ends one base into a stop codon on the transcript's last base,
/// under the vertebrate mitochondrial code: polyadenylation completes the stop codon.
pub fn shown() -> Transcript {
    Transcript {
        codon_table: Some(2),
        ..transcript("ENST_SHOWN", &[MET, ALA, GLY, &STOP[..1]], 0)
    }
}

/// `shown` with `cds_end_NF`: the flag shows no stop codon.
pub fn shown_flagged() -> Transcript {
    Transcript {
        id: "ENST_SHOWN_FLAGGED",
        flags: vec!["cds_end_NF"],
        ..shown()
    }
}

/// The end of `shown` under the standard code, with no consensus id: no record shows the
/// stop codon completed, and none stands.
pub fn tail_not_shown() -> Transcript {
    Transcript {
        id: "ENST_TAIL",
        codon_table: None,
        ..shown()
    }
}

/// A coding sequence that ends two bases into a codon before a 3' UTR: no stop codon
/// stands.
pub fn open_tail() -> Transcript {
    transcript("ENST_OPEN_TAIL", &[MET, ALA, GLY, &STOP[..2]], 6)
}

/// Whole codons whose last names a residue on the transcript's last base: no stop codon
/// stands and no 3' UTR follows.
pub fn open_end() -> Transcript {
    transcript("ENST_OPEN_END", &[MET, ALA, GLY, GLU], 0)
}

/// `open_end` with `cds_end_NF`.
pub fn open_end_flagged() -> Transcript {
    Transcript {
        id: "ENST_OPEN_END_FLAGGED",
        flags: vec!["cds_end_NF"],
        ..open_end()
    }
}

/// `open_end` closed by a stop codon: a stop codon stands.
pub fn closed_end() -> Transcript {
    transcript("ENST_CLOSED_END", &[MET, ALA, GLY, GLU, STOP], 0)
}

/// The coding sequence of `open_end` before a 3' UTR of six bases that the cache holds.
pub fn open_end_before_a_utr() -> Transcript {
    Transcript {
        three_prime_utr: Some("CCCCCC"),
        ..transcript("ENST_OPEN_END_UTR", &[MET, ALA, GLY, GLU], 6)
    }
}

/// The same before a 3' UTR that the cache does not hold.
pub fn open_end_before_a_utr_not_cached() -> Transcript {
    transcript("ENST_OPEN_END_UTR_NOT_CACHED", &[MET, ALA, GLY, GLU], 6)
}

/// A coding sequence that opens one base into a codon, with `cds_start_NF` and no 5' UTR,
/// and ends two bases into a codon on the transcript's last base; adenines after those
/// two spell no stop codon.
pub fn padded_tail_end() -> Transcript {
    Transcript {
        flags: vec!["cds_start_NF"],
        leader: 0,
        start_phase: 1,
        ..transcript("ENST_PADDED_TAIL_END", &["NGA", ALA, GLY, "TC"], 0)
    }
}

/// A first codon that is no start codon, closed by a stop codon: no start codon is shown.
pub fn start_open() -> Transcript {
    transcript("ENST_START_OPEN", &[GLU, ALA, GLY, STOP], 6)
}

/// `start_open` with `cds_start_NF`.
pub fn start_flagged() -> Transcript {
    Transcript {
        id: "ENST_START_FLAGGED",
        flags: vec!["cds_start_NF"],
        ..start_open()
    }
}

/// A first codon that is no start codon before three more codons and a stop codon, with
/// the 3' UTR the cache holds: the sequence on which a frameshift at the start is read to
/// the stop of its new frame.
pub fn start_read_on() -> Transcript {
    Transcript {
        three_prime_utr: Some("CTGACC"),
        ..transcript("ENST_START_READ", &["GAA", "AAA", "GGA", "TGG", "TAA"], 6)
    }
}

/// `start_read_on` under the vertebrate mitochondrial code.
pub fn start_read_on_under_another_code() -> Transcript {
    Transcript {
        id: "ENST_START_READ_MT",
        codon_table: Some(2),
        ..start_read_on()
    }
}

/// A first codon that is no start codon before the codons of alanine, glycine and
/// tryptophan and a stop codon, with the 3' UTR the cache holds: a deletion of its third
/// coding base keeps the first residue, changes the second, and can be written at no
/// place further 3'.
pub fn start_read_past() -> Transcript {
    Transcript {
        three_prime_utr: Some("CCCCCC"),
        ..transcript("ENST_START_PAST", &[GLU, ALA, GLY, "TGG", STOP], 6)
    }
}

/// The codons of the transcripts on which a frameshift named at the first residue is
/// counted: proline, isoleucine, serine and glycine before a stop codon. One base on, the
/// frame reads glutamine and a stop codon; two bases on, it reads asparagine, lysine,
/// arginine and isoleucine and runs on into the 3' UTR.
pub const COUNTED_CODONS: [&str; 5] = ["CCA", "ATA", "AGC", GLY, STOP];

/// A 3' UTR of six bases that holds a stop codon in the frame two bases on.
pub const STOP_IN_THE_UTR: &str = concat!("GC", "TAA", "C");

/// `COUNTED_CODONS` with `cds_start_NF`, no 5' UTR and the 3' UTR `STOP_IN_THE_UTR`, which
/// the cache holds: no start codon stands.
pub fn counted() -> Transcript {
    Transcript {
        flags: vec!["cds_start_NF"],
        leader: 0,
        three_prime_utr: Some(STOP_IN_THE_UTR),
        ..transcript("ENST_COUNTED", &COUNTED_CODONS, 6)
    }
}

/// `counted` behind six bases of 5' UTR.
pub fn counted_with_a_leader() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_LEADER",
        leader: LEADER,
        ..counted()
    }
}

/// `counted` on the reverse strand.
pub fn counted_reverse() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_REVERSE",
        strand: -1,
        ..counted()
    }
}

/// The four residues of `COUNTED_CODONS` with `cds_start_NF` and `cds_end_NF`, no 5' UTR
/// and no 3' UTR: the transcript ends at its coding end.
pub fn counted_flagged_at_both_ends() -> Transcript {
    Transcript {
        flags: vec!["cds_start_NF", "cds_end_NF"],
        leader: 0,
        ..transcript("ENST_COUNTED_BOTH", &COUNTED_CODONS[..4], 0)
    }
}

/// `counted` without the flag: its first codon is no start codon, so none is shown.
pub fn counted_without_the_flag() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_UNFLAGGED",
        flags: Vec::new(),
        ..counted()
    }
}

/// `counted_without_the_flag` behind six bases of 5' UTR.
pub fn counted_without_the_flag_with_a_leader() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_UNFLAGGED_LEADER",
        leader: LEADER,
        ..counted_without_the_flag()
    }
}

/// `counted` without the flag and with a first codon of methionine: a start codon stands.
pub fn counted_start_stands() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_STANDS",
        flags: Vec::new(),
        coding: [MET, "ATA", "AGC", GLY, STOP].concat(),
        ..counted()
    }
}

/// `counted` without the flag and with a consensus coding sequence id, which shows a start
/// codon.
pub fn counted_start_shown_by_an_id() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_ID",
        flags: Vec::new(),
        consensus_id: Some("CCDS1.1"),
        ..counted()
    }
}

/// `counted` before a cached 3' UTR that holds no stop codon in any frame.
pub fn counted_no_stop_in_the_utr() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_OPEN_UTR",
        three_prime_utr: Some(concat!("CCC", "CCC")),
        ..counted()
    }
}

/// `counted` before a 3' UTR that the cache does not hold.
pub fn counted_utr_not_cached() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_UTR_NOT_CACHED",
        three_prime_utr: None,
        ..counted()
    }
}

/// The residues of `COUNTED_CODONS` with `cds_start_NF` before one base of a stop codon on
/// the transcript's last base, with a consensus coding sequence id: polyadenylation
/// completes the stop codon.
pub fn counted_completed() -> Transcript {
    Transcript {
        flags: vec!["cds_start_NF"],
        leader: 0,
        consensus_id: Some("CCDS1.1"),
        ..transcript(
            "ENST_COUNTED_COMPLETED",
            &["CCA", "ATA", "AGC", GLY, &STOP[..1]],
            0,
        )
    }
}

/// `counted` with a first codon the start exon's phase pads with `phase` letters, 1 or 2.
pub fn counted_padded(phase: i64) -> Transcript {
    let (id, first_codon) = match phase {
        1 => ("ENST_COUNTED_PHASE_1", "NCA"),
        _ => ("ENST_COUNTED_PHASE_2", "NNA"),
    };
    Transcript {
        id,
        start_phase: phase,
        coding: [first_codon, "ATA", "AGC", GLY, STOP].concat(),
        ..counted()
    }
}

/// A flagged start under the vertebrate mitochondrial code: two bases on, the frame opens
/// with a codon of methionine under that code and of isoleucine under the standard one,
/// before a stop codon of both.
pub fn counted_under_another_code() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_MT",
        codon_table: Some(2),
        coding: ["CCA", "TAT", "AGC", GLY, STOP].concat(),
        ..counted()
    }
}

/// `counted` with a first codon that is a start codon of the standard code other than the
/// usual one, and a cached peptide that opens with the codon's own residue, leucine.
pub fn counted_other_start_codon() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_OTHER_START",
        coding: ["CTG", "ATA", "AGC", GLY, STOP].concat(),
        peptide: Some("LISG"),
        ..counted()
    }
}

/// `counted` with a letter that is no base in its fourth codon, which the frame two bases
/// on reads before its stop codon.
pub fn counted_unread_codon() -> Transcript {
    Transcript {
        id: "ENST_COUNTED_UNREAD",
        coding: ["CCA", "ATA", "AGC", "GNA", STOP].concat(),
        ..counted()
    }
}

/// Two codons with both flags before a cached 3' UTR that opens with a codon of lysine
/// and a stop codon.
pub fn counted_short() -> Transcript {
    Transcript {
        flags: vec!["cds_start_NF", "cds_end_NF"],
        leader: 0,
        three_prime_utr: Some(concat!("AAG", "TAA")),
        ..transcript("ENST_COUNTED_SHORT", &["CCA", "ATA"], 6)
    }
}

/// `closed` with an intron inside its stop codon, after the codon's second base.
pub fn split_stop() -> Transcript {
    Transcript {
        id: "ENST_SPLIT",
        intron_after: LEADER + 11,
        ..closed()
    }
}

/// One stable id on two records of the chromosome: one shaped as `open` at `FIRST`, one
/// shaped as `closed` four thousand bases on.
pub fn two_records() -> [Transcript; 2] {
    [
        Transcript {
            id: "ENST_TWO",
            ..open()
        },
        Transcript {
            id: "ENST_TWO",
            first: FIRST + 4000,
            ..closed()
        },
    ]
}

impl Transcript {
    pub fn cdna_length(&self) -> i64 {
        self.leader + self.coding.len() as i64 - self.start_phase + self.trailer
    }

    /// The genomic base of a cDNA position, on past the transcript's last base.
    pub fn genomic(&self, cdna: i64) -> i64 {
        if self.strand == -1 {
            return self.first + self.cdna_length() - cdna;
        }
        let intron = if self.intron_after > 0 && cdna > self.intron_after {
            100
        } else {
            0
        };
        self.first + cdna - 1 + intron
    }

    /// The base the forward strand holds where the transcript holds `base`.
    fn on_the_forward_strand(&self, base: char) -> char {
        match (self.strand, base) {
            (-1, 'A') => 'T',
            (-1, 'C') => 'G',
            (-1, 'G') => 'C',
            (-1, 'T') => 'A',
            _ => base,
        }
    }

    /// The transcript's own base at a cDNA position, a cytosine in either UTR and past
    /// the last base.
    pub fn base(&self, cdna: i64) -> char {
        usize::try_from(cdna - self.leader - 1 + self.start_phase)
            .ok()
            .filter(|_| cdna > self.leader)
            .and_then(|place| self.coding.chars().nth(place))
            .unwrap_or('C')
    }

    /// The transcript's own bases at the cDNA positions `first` to `last`.
    pub fn bases(&self, first: i64, last: i64) -> String {
        (first..=last).map(|cdna| self.base(cdna)).collect()
    }

    fn record(&self) -> String {
        let last = self.cdna_length();
        let (start, end) = (self.genomic(1), self.genomic(last));
        let coding_end = self.leader + self.coding.len() as i64 - self.start_phase;
        let (start, end) = (start.min(end), start.max(end));
        let lowest = |from: i64, to: i64| self.genomic(from).min(self.genomic(to));
        let highest = |from: i64, to: i64| self.genomic(from).max(self.genomic(to));
        let pair = |from: i64, to: i64| {
            format!(
                r#"{{"from_start": {from}, "from_end": {to}, "to_start": {}, "to_end": {}, "ori": {}}}"#,
                lowest(from, to),
                highest(from, to),
                self.strand
            )
        };
        let exon = |from: i64, to: i64| {
            format!(
                r#"{{"start": {}, "end": {}}}"#,
                lowest(from, to),
                highest(from, to)
            )
        };
        let spans: Vec<(i64, i64)> = match self.intron_after {
            0 => vec![(1, last)],
            after => vec![(1, after), (after + 1, last)],
        };
        let pairs: Vec<String> = spans.iter().map(|(from, to)| pair(*from, *to)).collect();
        let exons: Vec<String> = spans.iter().map(|(from, to)| exon(*from, *to)).collect();
        let attributes: Vec<String> = self
            .flags
            .iter()
            .map(|flag| format!(r#"{{"code": "{flag}", "value": "1"}}"#))
            .collect();
        let consensus = self
            .consensus_id
            .map_or(String::new(), |id| format!(r#", "ccds": "{id}""#));
        let table = self.codon_table.map_or(String::new(), |table| {
            format!(r#", "codon_table": {table}"#)
        });
        let peptide = self.peptide.map_or(String::new(), |peptide| {
            format!(r#", "peptide": "{peptide}""#)
        });
        let three_prime_utr = self.three_prime_utr.map_or(String::new(), |bases| {
            format!(r#", "three_prime_utr": "{bases}""#)
        });
        format!(
            r#"{{"stable_id": "{id}", "strand": {strand}, "start": {start}, "end": {end},
                "cdna_coding_start": {coding_start}, "cdna_coding_end": {coding_end},
                "exons": [{exons}], "attributes": [{attributes}]{consensus},
                "variation_effect_feature_cache": {{"translateable_seq": "{coding}"{table}{peptide}{three_prime_utr},
                "mapper": {{"start_phase": {start_phase}, "cdna_coding_start": {coding_start},
                "cdna_coding_end": {coding_end}, "pairs": [{pairs}]}}}}}}"#,
            id = self.id,
            strand = self.strand,
            start_phase = self.start_phase,
            coding_start = self.leader + 1,
            exons = exons.join(", "),
            attributes = attributes.join(", "),
            coding = self.coding,
            pairs = pairs.join(", "),
        )
    }
}

/// Plants a cache of `transcripts` under `directory` and gives its path. A record without
/// a coding sequence, `ENST_NONCODING`, stands beside them.
pub fn planted_cache(directory: &Path, transcripts: &[Transcript]) -> PathBuf {
    planted_cache_on(directory, CHROMOSOME, transcripts)
}

/// `planted_cache` with the transcripts on the chromosome `chromosome`.
pub fn planted_cache_on(directory: &Path, chromosome: &str, transcripts: &[Transcript]) -> PathBuf {
    let cache = directory.join("cache");
    let shard = cache.join("transcripts").join(chromosome);
    std::fs::create_dir_all(&shard).expect("the shard directory is made");
    let mut records: Vec<String> = transcripts.iter().map(Transcript::record).collect();
    records.push(
        r#"{"stable_id": "ENST_NONCODING", "strand": 1, "start": 1001, "end": 1021,
            "variation_effect_feature_cache": {}}"#
            .to_string(),
    );
    std::fs::write(
        shard.join("1-1000000.json"),
        format!("[{}]", records.join(",\n")),
    )
    .expect("the shard is written");
    cache
}

/// Every made-up transcript, planted as one cache.
pub fn whole_cache(directory: &Path) -> PathBuf {
    let [first, second] = two_records();
    planted_cache(
        directory,
        &[
            open(),
            closed(),
            flagged(),
            shown(),
            shown_flagged(),
            tail_not_shown(),
            open_tail(),
            open_end(),
            open_end_flagged(),
            closed_end(),
            open_end_before_a_utr(),
            open_end_before_a_utr_not_cached(),
            padded_tail_end(),
            start_open(),
            start_flagged(),
            start_read_on(),
            start_read_on_under_another_code(),
            start_read_past(),
            split_stop(),
            first,
            second,
            counted(),
            counted_with_a_leader(),
            counted_reverse(),
            counted_flagged_at_both_ends(),
            counted_without_the_flag(),
            counted_without_the_flag_with_a_leader(),
            counted_start_stands(),
            counted_start_shown_by_an_id(),
            counted_no_stop_in_the_utr(),
            counted_utr_not_cached(),
            counted_completed(),
            counted_padded(1),
            counted_padded(2),
            counted_under_another_code(),
            counted_other_start_codon(),
            counted_unread_codon(),
            counted_short(),
        ],
    )
}

/// One engine's row of a record on a transcript: its terms, the five cells in the default
/// format's spelling, and further fields by name (`HGVSp`, a plugin's).
#[derive(Clone)]
pub struct Cells {
    pub terms: String,
    pub cells: [String; 5],
    pub more: Vec<(String, String)>,
}

pub fn cells(terms: &str, five: [&str; 5]) -> Cells {
    Cells {
        terms: terms.to_string(),
        cells: five.map(str::to_string),
        more: Vec::new(),
    }
}

pub const BLANK: [&str; 5] = ["-", "-", "-", "-", "-"];

impl Cells {
    pub fn with(mut self, name: &str, value: &str) -> Cells {
        self.more.push((name.to_string(), value.to_string()));
        self
    }

    pub fn terms(mut self, terms: &str) -> Cells {
        self.terms = terms.to_string();
        self
    }

    /// The cells with `/length` after each position that is written, as `--total_length`
    /// prints them.
    pub fn total_length(mut self, lengths: [i64; 3]) -> Cells {
        for (cell, length) in self.cells.iter_mut().zip(lengths) {
            if cell != "-" {
                *cell = format!("{cell}/{length}");
            }
        }
        self
    }

    fn impact(&self) -> &'static str {
        impact_of(&self.terms)
    }
}

/// The lowest-ranked of the terms `terms` joins by commas (OutputFactory/JSON.pm:387-388).
pub fn most_severe(terms: &str) -> &str {
    let rank = |term: &str| {
        CONSEQUENCE_RANK
            .iter()
            .find(|(ranked, _)| *ranked == term)
            .map_or(1_000_000, |(_, rank)| *rank)
    };
    terms
        .split(',')
        .min_by_key(|term| (rank(term), *term))
        .unwrap_or("?")
}

/// The IMPACT of the most severe of `terms`.
pub fn impact_of(terms: &str) -> &'static str {
    match most_severe(terms) {
        "stop_gained"
        | "frameshift_variant"
        | "stop_lost"
        | "start_lost"
        | "feature_truncation"
        | "splice_acceptor_variant"
        | "splice_donor_variant" => "HIGH",
        "inframe_insertion"
        | "inframe_deletion"
        | "missense_variant"
        | "protein_altering_variant" => "MODERATE",
        "incomplete_terminal_codon_variant"
        | "start_retained_variant"
        | "stop_retained_variant"
        | "synonymous_variant"
        | "splice_region_variant" => "LOW",
        _ => "MODIFIER",
    }
}

/// One made-up record and its rows: the genomic bases its allele stands in for (`first` to
/// `last`, an insertion's the two bases it stands between with no reference bases), the
/// allele, and on each transcript Ensembl VEP's row and vep-rs's.
#[derive(Clone)]
pub struct Record {
    pub name: String,
    pub first: i64,
    pub last: i64,
    /// The bases the allele stands in for, empty for an insertion.
    pub reference: String,
    /// The base before the allele's bases, which a VCF line spells first.
    pub anchor: char,
    pub allele: String,
    pub rows: Vec<(&'static str, Cells, Cells)>,
    /// The `most_severe_consequence` each engine's JSON object states, where it is not the
    /// most severe of the engine's own terms.
    pub most_severe: Option<(String, String)>,
}

impl Record {
    /// A record on `on` that replaces the cDNA positions `first` to `last` by `allele`
    /// (`-` for none), with Ensembl VEP's row `theirs` and vep-rs's row `ours`.
    pub fn replacing(
        name: &str,
        on: &Transcript,
        first: i64,
        last: i64,
        allele: &str,
        theirs: Cells,
        ours: Cells,
    ) -> Record {
        if on.strand == -1 {
            // The record spells the forward strand: the transcript's bases reversed and
            // complemented behind the base that follows them on the transcript. `allele`
            // is the row's cell, which holds the forward strand's bases too.
            return Record {
                name: name.to_string(),
                first: on.genomic(last),
                last: on.genomic(first),
                reference: on
                    .bases(first, last)
                    .chars()
                    .rev()
                    .map(|base| on.on_the_forward_strand(base))
                    .collect(),
                anchor: on.on_the_forward_strand(on.base(last + 1)),
                allele: allele.to_string(),
                rows: vec![(on.id, theirs, ours)],
                most_severe: None,
            };
        }
        Record {
            name: name.to_string(),
            first: on.genomic(first),
            last: on.genomic(first) + last - first,
            reference: on.bases(first, last),
            anchor: on.base(first - 1),
            allele: allele.to_string(),
            rows: vec![(on.id, theirs, ours)],
            most_severe: None,
        }
    }

    /// A record on `on` that inserts `allele` between the cDNA positions `after` and the
    /// next.
    pub fn inserting(
        name: &str,
        on: &Transcript,
        after: i64,
        allele: &str,
        theirs: Cells,
        ours: Cells,
    ) -> Record {
        let (first, anchor) = match on.strand {
            -1 => (
                on.genomic(after + 1),
                on.on_the_forward_strand(on.base(after + 1)),
            ),
            _ => (on.genomic(after), on.base(after)),
        };
        Record {
            name: name.to_string(),
            first,
            last: first + 1,
            reference: String::new(),
            anchor,
            allele: allele.to_string(),
            rows: vec![(on.id, theirs, ours)],
            most_severe: None,
        }
    }

    /// The record with one more row, on the transcript `id`.
    pub fn and(mut self, id: &'static str, theirs: Cells, ours: Cells) -> Record {
        self.rows.push((id, theirs, ours));
        self
    }

    /// The record with its one row on the transcript `id` instead.
    pub fn on(mut self, id: &'static str) -> Record {
        self.rows[0].0 = id;
        self
    }

    /// The record with vep-rs's row of its first transcript replaced.
    pub fn ours(mut self, ours: Cells) -> Record {
        self.rows[0].2 = ours;
        self
    }

    /// The record moved `by` bases along the chromosome.
    pub fn moved(mut self, by: i64) -> Record {
        self.first += by;
        self.last += by;
        self
    }

    fn insertion(&self) -> bool {
        self.reference.is_empty() && self.structural().is_none()
    }

    /// The class word of a structural allele.
    fn structural(&self) -> Option<&str> {
        self.allele
            .chars()
            .any(|letter| letter.is_ascii_lowercase())
            .then_some(self.allele.as_str())
    }

    fn location(&self) -> String {
        if self.first == self.last {
            format!("{CHROMOSOME}:{}", self.first)
        } else {
            format!("{CHROMOSOME}:{}-{}", self.first, self.last)
        }
    }

    /// The record's line of the input VCF.
    pub fn input_line(&self) -> String {
        if self.structural().is_some() {
            return format!(
                "{CHROMOSOME}\t{}\t{}\t{}\t<DEL>\t.\t.\tSVTYPE=DEL;END={}",
                self.first - 1,
                self.name,
                self.anchor,
                self.last
            );
        }
        let (position, reference, alternate) = if self.insertion() {
            (
                self.first,
                self.anchor.to_string(),
                format!("{}{}", self.anchor, self.allele),
            )
        } else if self.allele == "-" || self.allele.len() != self.reference.len() {
            let kept = if self.allele == "-" { "" } else { &self.allele };
            (
                self.first - 1,
                format!("{}{}", self.anchor, self.reference),
                format!("{}{kept}", self.anchor),
            )
        } else {
            (self.first, self.reference.clone(), self.allele.clone())
        };
        format!(
            "{CHROMOSOME}\t{position}\t{}\t{reference}\t{alternate}\t.\t.\t.",
            self.name
        )
    }

    fn side(&self, row: &(&'static str, Cells, Cells), ours: bool) -> Cells {
        if ours {
            row.2.clone()
        } else {
            row.1.clone()
        }
    }

    fn default_row(&self, id: &str, cells: &Cells) -> String {
        let mut extra = vec![
            format!("IMPACT={}", cells.impact()),
            format!("STRAND={}", strand_of(id)),
        ];
        extra.extend(
            cells
                .more
                .iter()
                .map(|(name, value)| format!("{name}={value}")),
        );
        format!(
            "{}\t{}\t{}\tENSG1\t{id}\tTranscript\t{}\t{}\t-\t{}",
            self.name,
            self.location(),
            self.allele,
            cells.terms,
            cells.cells.join("\t"),
            extra.join(";")
        )
    }

    fn tab_row(&self, id: &str, cells: &Cells, more: &[String]) -> String {
        let further: Vec<String> = more
            .iter()
            .map(|name| {
                cells
                    .more
                    .iter()
                    .find(|(held, _)| held == name)
                    .map_or("-".to_string(), |(_, value)| value.clone())
            })
            .collect();
        format!(
            "{}\t{}\t{}\tENSG1\t{id}\tTranscript\t{}\t{}\t-\t{}\t{}{}",
            self.name,
            self.location(),
            self.allele,
            cells.terms,
            cells.cells.join("\t"),
            cells.impact(),
            strand_of(id),
            further
                .iter()
                .map(|value| format!("\t{value}"))
                .collect::<String>()
        )
    }

    fn csq(&self, id: &str, cells: &Cells, more: &[String]) -> String {
        let blank = |cell: &str| {
            if cell == "-" {
                String::new()
            } else {
                cell.to_string()
            }
        };
        let mut values = vec![
            self.allele.clone(),
            cells.terms.replace(',', "&"),
            cells.impact().to_string(),
            "Transcript".to_string(),
            id.to_string(),
        ];
        values.extend(cells.cells.iter().map(|cell| blank(cell)));
        values.push(strand_of(id).to_string());
        for name in more {
            let value = cells
                .more
                .iter()
                .find(|(held, _)| held == name)
                .map_or(String::new(), |(_, value)| value.clone());
            values.push(value);
        }
        values.join("|")
    }

    fn json_consequence(&self, id: &str, cells: &Cells) -> String {
        let quoted = |text: &str| pyjson::dumps_str(text, true);
        let terms: Vec<String> = cells.terms.split(',').map(quoted).collect();
        let mut members = vec![
            format!(r#""variant_allele": {}"#, quoted(&self.allele)),
            format!(r#""consequence_terms": [{}]"#, terms.join(", ")),
            format!(r#""impact": {}"#, quoted(cells.impact())),
            format!(r#""transcript_id": {}"#, quoted(id)),
            r#""gene_id": "ENSG1""#.to_string(),
            format!(r#""strand": {}"#, strand_of(id)),
        ];
        for (cell, start, end) in [
            (&cells.cells[0], "cdna_start", "cdna_end"),
            (&cells.cells[1], "cds_start", "cds_end"),
            (&cells.cells[2], "protein_start", "protein_end"),
        ] {
            if cell == "-" {
                continue;
            }
            let (first, last) = cell.split_once('-').unwrap_or((cell, cell));
            for (name, bound) in [(start, first), (end, last)] {
                if bound != "?" {
                    members.push(format!(r#""{name}": {bound}"#));
                }
            }
        }
        for (cell, name) in [
            (&cells.cells[3], "amino_acids"),
            (&cells.cells[4], "codons"),
        ] {
            if cell != "-" {
                members.push(format!(r#""{name}": {}"#, quoted(cell)));
            }
        }
        for (name, value) in &cells.more {
            let value = value.replace("%3D", "=");
            members.push(format!(r#""{}": {}"#, name.to_lowercase(), quoted(&value)));
        }
        format!("{{{}}}", members.join(", "))
    }

    fn json_object(&self, ours: bool) -> String {
        let quoted = |text: &str| pyjson::dumps_str(text, true);
        let consequences: Vec<String> = self
            .rows
            .iter()
            .map(|row| self.json_consequence(row.0, &self.side(row, ours)))
            .collect();
        let own: Vec<String> = self
            .rows
            .iter()
            .map(|row| self.side(row, ours).terms)
            .collect();
        let own = own.join(",");
        let stated = self
            .most_severe
            .as_ref()
            .map(|(theirs, written)| if ours { written } else { theirs });
        let severe = stated.map_or(most_severe(&own), String::as_str);
        let (start, end) = if self.insertion() {
            (self.last, self.first)
        } else {
            (self.first, self.last)
        };
        let reference = if self.reference.is_empty() {
            "-"
        } else {
            &self.reference
        };
        // A structural record's object spells the record's symbolic ALT.
        let allele_string = match self.structural() {
            Some(_) => "<DEL>".to_string(),
            None => format!("{reference}/{}", self.allele),
        };
        format!(
            r#"{{"input": {}, "id": {}, "seq_region_name": "{CHROMOSOME}", "start": {start}, "end": {end}, "strand": 1, "allele_string": {}, "most_severe_consequence": {}, "transcript_consequences": [{}]}}"#,
            quoted(&self.input_line()),
            quoted(&self.name),
            quoted(&allele_string),
            quoted(severe),
            consequences.join(", ")
        )
    }
}

/// The strand a row of the transcript `id` states: -1 for an id that ends in `_REVERSE`.
fn strand_of(id: &str) -> i64 {
    if id.ends_with("_REVERSE") {
        -1
    } else {
        1
    }
}

const DEFAULT_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";
const TAB_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tSTRAND";
const CSQ_FIELDS: &str = "Allele|Consequence|IMPACT|Feature_type|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|STRAND";

/// The names of the further fields the rows of `records` carry, in the order first met.
fn further_fields(records: &[Record]) -> Vec<String> {
    let mut names: Vec<String> = Vec::new();
    for record in records {
        for (_, theirs, ours) in &record.rows {
            for (name, _) in theirs.more.iter().chain(&ours.more) {
                if !names.contains(name) {
                    names.push(name.clone());
                }
            }
        }
    }
    names
}

/// The output of one engine over `records` in the format `fmt`.
pub fn output(records: &[Record], fmt: &str, ours: bool) -> String {
    let more = further_fields(records);
    let mut lines: Vec<String> = Vec::new();
    match fmt {
        "default" => lines.push(DEFAULT_HEADER.to_string()),
        "tab" => lines.push(
            std::iter::once(TAB_HEADER.to_string())
                .chain(more.iter().cloned())
                .collect::<Vec<_>>()
                .join("\t"),
        ),
        "vcf" => {
            let fields: Vec<&str> = std::iter::once(CSQ_FIELDS)
                .chain(more.iter().map(String::as_str))
                .collect();
            lines.push(format!(
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {}\">",
                fields.join("|")
            ));
            lines.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string());
        }
        _ => {}
    }
    for record in records {
        match fmt {
            "default" => lines.extend(
                record
                    .rows
                    .iter()
                    .map(|row| record.default_row(row.0, &record.side(row, ours))),
            ),
            "tab" => lines.extend(
                record
                    .rows
                    .iter()
                    .map(|row| record.tab_row(row.0, &record.side(row, ours), &more)),
            ),
            "vcf" => {
                let entries: Vec<String> = record
                    .rows
                    .iter()
                    .map(|row| record.csq(row.0, &record.side(row, ours), &more))
                    .collect();
                let line = record.input_line();
                // An INFO of `.` gives way to the CSQ; any other gains it as one more item.
                let open = match line.strip_suffix('.') {
                    Some(open) => open.to_string(),
                    None => format!("{line};"),
                };
                lines.push(format!("{open}CSQ={}", entries.join(",")));
            }
            _ => lines.push(record.json_object(ours)),
        }
    }
    lines.join("\n") + "\n"
}

/// The input VCF of `records`.
pub fn input(records: &[Record]) -> String {
    let mut lines = vec![
        "##fileformat=VCFv4.2".to_string(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
    ];
    lines.extend(records.iter().map(Record::input_line));
    lines.join("\n") + "\n"
}

pub fn text(path: &Path) -> String {
    path.to_str().expect("a path of text").to_string()
}

/// The arguments that compare the two engines' outputs over `records` in the format
/// `fmt`, planted under `directory`, with `more` after them: `--input` is given for the
/// VCF and JSON formats, and for a text format where `with_input` asks.
pub fn planted(
    directory: &Path,
    records: &[Record],
    fmt: &str,
    with_input: bool,
    more: &[&str],
) -> Vec<OsString> {
    std::fs::create_dir_all(directory).expect("the directory is made");
    let reference = directory.join(format!("reference.{fmt}"));
    let candidate = directory.join(format!("candidate.{fmt}"));
    std::fs::write(&reference, output(records, fmt, false)).expect("the reference is written");
    std::fs::write(&candidate, output(records, fmt, true)).expect("the candidate is written");
    let mut argv = vec![
        "--format".to_string(),
        fmt.to_string(),
        "--expected".to_string(),
        text(&reference),
        "--actual".to_string(),
        text(&candidate),
    ];
    if with_input || matches!(fmt, "vcf" | "json") {
        let records_path = directory.join("input.vcf");
        std::fs::write(&records_path, input(records)).expect("the input is written");
        argv.extend(["--input".to_string(), text(&records_path)]);
    }
    argv.extend(more.iter().map(|argument| argument.to_string()));
    argv.into_iter().map(OsString::from).collect()
}

/// The report of the command line `argv`, or what the comparison refuses it with.
pub fn compared(argv: &[OsString]) -> Result<Value, String> {
    let args = cli::parse_args(&cli::build_parser(), argv).map_err(|_| "no parse".to_string())?;
    cli::run(&args, None).map_err(|refused| refused.to_string())
}

fn compact(value: &Value) -> String {
    pyjson::dumps(
        value,
        &pyjson::DumpOptions {
            indent: None,
            sort_keys: true,
            ensure_ascii: false,
        },
    )
}

/// What a report counts as applied and its reasons to be red, each as compact JSON with
/// its keys in order.
pub fn applied_and_reasons(report: &Value) -> (String, String) {
    let member = |path: &[&str]| {
        path.iter().fold(report, |value, name| {
            value
                .get(name)
                .unwrap_or_else(|| panic!("the report has no member {name}"))
        })
    };
    (
        compact(member(&["masks", "applied"])),
        compact(member(&["reasons"])),
    )
}

/// `applied_and_reasons` of the comparison of `records` in `fmt` beside the planted cache
/// of every made-up transcript, the classes `masks` named and `--input` given.
pub fn beside_the_cache(
    name: &str,
    records: &[Record],
    fmt: &str,
    masks: &str,
) -> (String, String) {
    beside(name, records, fmt, masks, true)
}

/// The same with or without the cache.
pub fn beside(
    name: &str,
    records: &[Record],
    fmt: &str,
    masks: &str,
    cache: bool,
) -> (String, String) {
    let directory = scratch(name);
    let planted_at = text(&whole_cache(&directory));
    let mut more = vec!["--masks", masks];
    if cache {
        more.extend(["--json-cache", planted_at.as_str()]);
    }
    let argv = planted(&directory, records, fmt, true, &more);
    let report = compared(&argv).unwrap_or_else(|refused| panic!("{name}: {refused}"));
    applied_and_reasons(&report)
}

/// An empty directory of this test's own.
pub fn scratch(name: &str) -> PathBuf {
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("vep-formats-compare-coding-end-{name}"));
    if directory.exists() {
        std::fs::remove_dir_all(&directory).expect("the old scratch directory is removed");
    }
    std::fs::create_dir_all(&directory).expect("the scratch directory is made");
    directory
}

/// The reasons of a report with nothing left, and the applied counts of one class.
pub const GREEN: &str = "[]";

pub fn one(class: &str) -> String {
    format!("{{\"{class}\": 1}}")
}

/// The two sides of a row that both engines print with the same cells.
pub fn pair(theirs: &str, ours: &str, five: [&str; 5]) -> (Cells, Cells) {
    (cells(theirs, five), cells(ours, five))
}

// The rows of the kinds the classes take.

/// A deletion of one base in the last whole codon of `on` (the stop codon's overlap gives
/// Ensembl VEP a kept stop beside an in-frame term; vep-rs owes a frameshift).
pub fn kept_stop_deletion(on: &Transcript) -> Record {
    let (theirs, ours) = pair(
        "inframe_deletion,stop_retained_variant",
        "frameshift_variant",
        ["14", "8", "3", "G/X", "gGa/ga"],
    );
    Record::replacing("kept", on, 14, 14, "-", theirs, ours)
}

/// The same deletion where Ensembl VEP writes a lost stop beside the frameshift.
pub fn lost_stop_deletion(on: &Transcript) -> Record {
    let (theirs, ours) = pair(
        "frameshift_variant,stop_lost",
        "frameshift_variant",
        ["14", "8", "3", "G/X", "gGa/ga"],
    );
    Record::replacing("lost", on, 14, 14, "-", theirs, ours)
}

/// A deletion from the last whole codon of `on` past the transcript's last base.
pub fn deletion_past_the_end(on: &Transcript) -> Record {
    let (theirs, ours) = pair(
        "stop_retained_variant,3_prime_UTR_variant",
        "coding_sequence_variant,3_prime_UTR_variant",
        ["14-?", "8-?", "3-?", "-", "-"],
    );
    Record::replacing("past", on, 14, 30, "-", theirs, ours)
}

/// An insertion of a stop codon's bases behind the last residue of `on`.
pub fn stop_behind_the_last_residue(on: &Transcript) -> Record {
    let allele = format!("C{}", &STOP[..2]);
    let codons = format!("gga/gg{allele}a");
    let five = ["14-15", "8-9", "3", "G/G*", codons.as_str()];
    let (theirs, ours) = pair(
        "inframe_insertion,stop_retained_variant",
        "stop_gained,inframe_insertion",
        five,
    );
    Record::inserting("behind", on, 14, &allele, theirs, ours)
}

/// A substitution of the one base of the stop codon `shown` holds: Ensembl VEP reads an
/// incomplete codon, vep-rs the completed stop codon lost.
pub fn completed_stop_lost(on: &Transcript) -> Record {
    let theirs = cells(
        "incomplete_terminal_codon_variant,coding_sequence_variant",
        ["16", "10", "4", "X", "T/C"],
    );
    let ours = cells("stop_lost", ["16", "10", "4", "*/Q", "Taa/Caa"])
        .with("HGVSp", "ENSP_SHOWN.1:p.Ter4GlnextTer?");
    Record::replacing("completed", on, 16, 16, "C", theirs, ours)
}

/// A structural deletion over the coding end of `on`.
pub fn structural_deletion(on: &Transcript) -> Record {
    let (theirs, ours) = pair(
        "stop_lost,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
        "coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
        BLANK,
    );
    Record {
        name: "structural".to_string(),
        first: on.genomic(10),
        last: on.genomic(40),
        reference: String::new(),
        anchor: on.base(9),
        allele: "deletion".to_string(),
        rows: vec![(on.id, theirs, ours)],
        most_severe: None,
    }
}

/// A substitution in the first codon of `on` that changes its residue: Ensembl VEP writes
/// a lost start, vep-rs the missense change.
pub fn first_residue_changed(on: &Transcript) -> Record {
    let five = ["9", "3", "1", "E/D", "gaA/gaT"];
    let theirs = cells("start_lost", five).with("HGVSp", "ENSP_START.1:p.Glu1?");
    let ours = cells("missense_variant", five).with("HGVSp", "ENSP_START.1:p.Glu1Asp");
    Record::replacing("first", on, 9, 9, "T", theirs, ours)
}
