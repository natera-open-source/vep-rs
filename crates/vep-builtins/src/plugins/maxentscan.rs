// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! MaxEntScan plugin: maximum-entropy splice-site scores of the reference and
//! alternate sequences around a variant, computed from the MaxEntScan matrices
//! as Ensembl VEP's `MaxEntScan.pm`
//! computes them.
//!
//! **Usage:** `--plugin MaxEntScan,/path/to/fordownload[,SWA][,NCSS][,verbose]`
//!
//! The first parameter is the directory the plugin's documentation has unpacked
//! from the MaxEntScan archive (MaxEntScan.pm:40-43). The plugin reads `me2x5`
//! and `splicemodels/splice5sequences` there (MaxEntScan.pm:104-105) for the
//! 5' (donor) model, which scores a
//! 9-mer of the last 3 exonic and first 6 intronic bases, and
//! `splicemodels/me2x3acc1` to `me2x3acc9` for the 3' (acceptor) model, which
//! scores a 23-mer of the last 20 intronic and first 3 exonic bases. The plugin
//! needs the reference FASTA (`--fasta`), as the Ensembl VEP plugin does offline.
//!
//! Three groups of fields:
//!
//! - `MaxEntScan_ref`, `MaxEntScan_alt`, `MaxEntScan_diff` (reference minus
//!   alternate), for a single-base substitution inside a canonical donor or
//!   acceptor window of one of the transcript's introns (`run_MES`,
//!   MaxEntScan.pm:256-332).
//! - With `SWA`, a sliding window over the sequence around the variant
//!   (`run_SWA`, :334-418): the highest-scoring 9-mer and 23-mer of the
//!   reference and of the alternate context, a reference comparison score and
//!   the difference, as `MES-SWA_donor_*` and `MES-SWA_acceptor_*`.
//! - With `NCSS`, the nearest canonical splice sites up- and downstream of the
//!   variant's exon or intron (`run_NCSS`, :420-511), as `MES-NCSS_*`.
//!
//! `verbose` adds the sequences, k-mer positions and contexts behind the scores.
//! Every score prints with three decimals (`:245`).
//!
//! Per-transcript plugin over [`crate::BuiltinPlugin::run_batch`]: one reference
//! window around each variant serves every transcript row of the variant, read
//! on the transcript's strand.

use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use vep_core::codon::{reverse_complement, seq_is_dna};
use vep_core::consequence::{LofteeContext, TranscriptConsequence};
use vep_core::transcript::{Exon, Intron};
use vep_core::variant::InputVariant;

use crate::params::{warn_unread_parameters, PluginParams};
use crate::tabix::{lossy_text, BoundedLines};
use crate::traits::{BuiltinPlugin, PluginError};

const NAME: &str = "MaxEntScan";

/// Base frequencies and consensus probabilities of `score5_scoreconsensus` and
/// `score3_scoreconsensus` (MaxEntScan.pm:771-792, 871-892), indexed A, C, G, T.
const BACKGROUND: [f64; 4] = [0.27, 0.23, 0.23, 0.27];
const DONOR_CONSENSUS_1: [f64; 4] = [0.004, 0.0032, 0.9896, 0.0032];
const DONOR_CONSENSUS_2: [f64; 4] = [0.0034, 0.0039, 0.0042, 0.9884];
const ACCEPTOR_CONSENSUS_1: [f64; 4] = [0.9903, 0.0032, 0.0034, 0.0030];
const ACCEPTOR_CONSENSUS_2: [f64; 4] = [0.0027, 0.0037, 0.9905, 0.0030];

/// The entries the reader requires of `me2x5` and the lines it requires of
/// `splice5sequences`: one per 7-mer over ACGT, as the plugin dies on a 7-mer
/// without a line (MaxEntScan.pm:683-687).
const DONOR_TABLE_LEN: usize = 1 << (2 * 7);

/// The nine acceptor tables: file name, k-mer length and the slice of the 21
/// non-consensus bases each scores (`score3_maxentscore`, MaxEntScan.pm:851-859).
const ACCEPTOR_TABLES: [(&str, usize, usize, usize); 9] = [
    ("me2x3acc1", 7, 0, 7),
    ("me2x3acc2", 7, 7, 14),
    ("me2x3acc3", 7, 14, 21),
    ("me2x3acc4", 7, 4, 11),
    ("me2x3acc5", 7, 11, 18),
    ("me2x3acc6", 3, 4, 7),
    ("me2x3acc7", 4, 7, 11),
    ("me2x3acc8", 3, 11, 14),
    ("me2x3acc9", 4, 14, 18),
];

/// The flank of reference sequence read once per variant: the widest window any
/// group reads around the variant is the sliding-window acceptor context of 22
/// bases each side (`:380`), and the canonical windows searched for a
/// substitution lie within 21 bases of it (`:270`) plus their own 3-base exonic
/// flank; 43 covers both.
const WINDOW_FLANK: u64 = 43;

/// A scoring function of the tables: `score5` over a 9-mer or `score3` over a 23-mer.
type Score = fn(&Tables, &[u8]) -> f64;

/// Whether every byte is one of the four bases (`/^[ACGT]+$/`).
fn is_acgt(seq: &[u8]) -> bool {
    !seq.is_empty() && seq.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

fn base_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    }
}

/// The base-4 hash of an ACGT k-mer (`score3_hashseq`, MaxEntScan.pm:803-820):
/// A 0, C 1, G 2, T 3, the first base most significant.
fn hash(seq: &[u8]) -> usize {
    seq.iter().fold(0, |acc, &b| acc * 4 + base_index(b))
}

/// `log($val) / log(2)` (MaxEntScan.pm:794-797): the natural logarithms' quotient,
/// not a base-2 logarithm, so the last binary place agrees with the plugin's.
fn log2(value: f64) -> f64 {
    value.ln() / 2f64.ln()
}

/// Ensembl's `overlap`: two closed spans sharing a base.
fn overlap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> bool {
    a_end >= b_start && a_start <= b_end
}

/// The canonical donor window of an intron, a closed genomic span read on the
/// transcript's strand: the last 3 exonic bases and the first 6 intronic.
fn donor_of_intron(intron: &Intron, strand: i8) -> (i64, i64) {
    let (start, end) = (intron.start as i64, intron.end as i64);
    if strand > 0 {
        (start - 3, start + 5)
    } else {
        (end - 5, end + 3)
    }
}

/// The canonical acceptor window of an intron, a closed genomic span read on
/// the transcript's strand: the last 20 intronic bases and the first 3 exonic.
fn acceptor_of_intron(intron: &Intron, strand: i8) -> (i64, i64) {
    let (start, end) = (intron.start as i64, intron.end as i64);
    if strand > 0 {
        (end - 19, end + 3)
    } else {
        (start - 3, start + 19)
    }
}

/// The donor window at an exon's 3' end: its last 3 bases and the 6 after them.
fn donor_of_exon(exon: &Exon, strand: i8) -> (i64, i64) {
    let (start, end) = (exon.start as i64, exon.end as i64);
    if strand > 0 {
        (end - 2, end + 6)
    } else {
        (start - 6, start + 2)
    }
}

/// The acceptor window at an exon's 5' end: the 20 bases before it and its
/// first 3.
fn acceptor_of_exon(exon: &Exon, strand: i8) -> (i64, i64) {
    let (start, end) = (exon.start as i64, exon.end as i64);
    if strand > 0 {
        (start - 20, start + 2)
    } else {
        (end - 2, end + 20)
    }
}

/// The matrices as the plugin holds them after loading, each indexed by the
/// base-4 hash of the k-mer it scores.
struct Tables {
    /// `me2x5` entries by the hash of the 7-mer whose line of
    /// `splice5sequences` names them (`score5_makesequencematrix` and
    /// `score5_makescorematrix`, MaxEntScan.pm:732-762).
    donor: Vec<f64>,
    /// `me2x3acc1` to `me2x3acc9`, each by the hash of its k-mer
    /// (`score3_makemaxentscores`, :822-843).
    acceptor: Vec<Vec<f64>>,
}

impl Tables {
    fn load(dir: &Path) -> Result<Self, PluginError> {
        let scores = read_numbers(&dir.join("me2x5"), DONOR_TABLE_LEN)?;
        let sequences_path = dir.join("splicemodels").join("splice5sequences");
        let sequences = read_lines(&sequences_path, DONOR_TABLE_LEN)?;
        // The plugin keys the scores by the sequence file's text and a repeated
        // 7-mer keeps its last line (`$matrix{$_} = $n`, MaxEntScan.pm:732-762).
        let mut donor = vec![0.0; DONOR_TABLE_LEN];
        let mut listed = vec![false; DONOR_TABLE_LEN];
        for (n, seq) in sequences.iter().enumerate() {
            if seq.len() != 7 || !is_acgt(seq.as_bytes()) {
                return Err(PluginError::Init(format!(
                    "{NAME}: {} line {} is not a 7-mer over ACGT: '{seq}'",
                    sequences_path.display(),
                    n + 1
                )));
            }
            donor[hash(seq.as_bytes())] = scores[n];
            listed[hash(seq.as_bytes())] = true;
        }
        if listed.contains(&false) {
            return Err(PluginError::Init(format!(
                "{NAME}: {} does not list every 7-mer over ACGT once",
                sequences_path.display()
            )));
        }
        let acceptor = ACCEPTOR_TABLES
            .iter()
            .map(|(file, k, _, _)| read_numbers(&dir.join("splicemodels").join(file), 1 << (2 * k)))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { donor, acceptor })
    }

    /// `score5` (MaxEntScan.pm:665-696) of a 9-mer over ACGT: the consensus factor
    /// of bases 3 and 4 times the `me2x5` entry of the other seven bases, as log2.
    fn score5(&self, seq: &[u8]) -> f64 {
        let (b3, b4) = (base_index(seq[3]), base_index(seq[4]));
        let consensus =
            DONOR_CONSENSUS_1[b3] * DONOR_CONSENSUS_2[b4] / (BACKGROUND[b3] * BACKGROUND[b4]);
        let rest = [seq[0], seq[1], seq[2], seq[5], seq[6], seq[7], seq[8]];
        log2(consensus * self.donor[hash(&rest)])
    }

    /// `score3` (MaxEntScan.pm:698-726) of a 23-mer over ACGT: the consensus factor
    /// of bases 18 and 19 times the maximum-entropy ratio of the other 21 bases,
    /// five tables over four in the plugin's order (`:860`), as log2.
    fn score3(&self, seq: &[u8]) -> f64 {
        let (b18, b19) = (base_index(seq[18]), base_index(seq[19]));
        let consensus = ACCEPTOR_CONSENSUS_1[b18] * ACCEPTOR_CONSENSUS_2[b19]
            / (BACKGROUND[b18] * BACKGROUND[b19]);
        let mut rest = [0u8; 21];
        rest[..18].copy_from_slice(&seq[..18]);
        rest[18..].copy_from_slice(&seq[20..23]);
        let mut sc = [0.0; 9];
        for (slot, ((_, _, from, to), table)) in sc
            .iter_mut()
            .zip(ACCEPTOR_TABLES.iter().zip(&self.acceptor))
        {
            *slot = table[hash(&rest[*from..*to])];
        }
        let ratio = sc[0] * sc[1] * sc[2] * sc[3] * sc[4] / (sc[5] * sc[6] * sc[7] * sc[8]);
        log2(consensus * ratio)
    }

    /// The highest-scoring k-mer of `sequence` under `score`: its 1-based
    /// position (`frame`) and score (`get_max_donor` and `get_max_acceptor`,
    /// MaxEntScan.pm:517-551); the first of equal maxima wins.
    fn best_kmer(&self, sequence: &[u8], k: usize, score: Score) -> (usize, f64) {
        let mut best: Option<(usize, f64)> = None;
        for (i, kmer) in sequence.windows(k).enumerate() {
            let s = score(self, kmer);
            if best.is_none_or(|(_, max)| s > max) {
                best = Some((i, s));
            }
        }
        let (i, max) = best.expect("a context holds at least one k-mer");
        (i + 1, max)
    }
}

/// Perl's `\s` on the byte string the module reads, and the whitespace its
/// numeric conversion skips: ASCII whitespace and the vertical tab.
fn is_perl_space(c: char) -> bool {
    c.is_ascii_whitespace() || c == '\x0B'
}

/// The lines of a matrix file, each with its first whitespace character removed
/// (`chomp; s/\s//`, MaxEntScan.pm:739-740; `\s` on the byte string the module
/// reads is ASCII whitespace and the vertical tab), checked for the line count
/// the plugin's lookups imply, a line per k-mer (MaxEntScan.pm:683-687 for the
/// donor files, the base-4 hash of MaxEntScan.pm:803-813 for the acceptor
/// tables). The file is read as the bytes the module
/// reads, a line's length capped as every flat file's is ([`BoundedLines`]); a
/// line holding a byte that is not UTF-8 is then refused by the check on its
/// content, as a number or a 7-mer, that every line meets.
fn read_lines(path: &Path, expected: usize) -> Result<Vec<String>, PluginError> {
    let file = std::fs::File::open(path).map_err(|e| {
        PluginError::Init(format!("{NAME}: failed to read {}: {e}", path.display()))
    })?;
    let mut lines = Vec::new();
    for line in BoundedLines::new(BufReader::new(file), path) {
        lines.push(lossy_text(line?).replacen(is_perl_space, "", 1));
    }
    if lines.len() != expected {
        return Err(PluginError::Init(format!(
            "{NAME}: {} has {} lines, expected {expected}",
            path.display(),
            lines.len()
        )));
    }
    Ok(lines)
}

/// The number of a kept line that holds a decimal number written with a
/// point or an exponent: the whitespace before and after it is skipped, and
/// the number is the one Perl's arithmetic gives that text where the plugin
/// multiplies the line (`$a * $d`, MaxEntScan.pm:689; the product at :860).
fn perl_number(text: &str) -> Option<f64> {
    text.trim_matches(is_perl_space).parse().ok()
}

/// A matrix file as numbers, one per line, each read as [`perl_number`] reads
/// it; a line that is no such number is refused by its file and line.
fn read_numbers(path: &Path, expected: usize) -> Result<Vec<f64>, PluginError> {
    read_lines(path, expected)?
        .iter()
        .enumerate()
        .map(|(n, line)| {
            perl_number(line).ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: {} line {} is not a number: '{line}'",
                    path.display(),
                    n + 1
                ))
            })
        })
        .collect()
}

/// The reference bases around one variant, read once and sliced per transcript
/// strand.
struct Window {
    chr_start: i64,
    bases: Vec<u8>,
}

impl Window {
    /// The forward-strand bases of the closed genomic span `start..=end`, or `None`
    /// when the span leaves the window, as `sub_Slice` returns `undef` for a span
    /// leaving the chromosome.
    fn forward(&self, start: i64, end: i64) -> Option<&[u8]> {
        let chr_end = self.chr_start + self.bases.len() as i64 - 1;
        if start < self.chr_start || end > chr_end || end < start {
            return None;
        }
        let from = (start - self.chr_start) as usize;
        Some(&self.bases[from..=from + (end - start) as usize])
    }
}

/// A score as the plugin prints it (`sprintf('%.3f', ...)`, MaxEntScan.pm:245):
/// three decimals, and Perl's `Inf` and `-Inf` for an infinity, which a table
/// entry that is one gives.
fn format_score(score: f64) -> String {
    if score.is_infinite() {
        return if score > 0.0 { "Inf" } else { "-Inf" }.to_string();
    }
    format!("{score:.3}")
}

/// The rendered fields of one group: every score as `%.3f` (`:245`) and, under
/// `verbose`, the sequences, k-mer positions and contexts behind them.
type Fields = Vec<(&'static str, String)>;

/// The fields of one side of the sliding window.
struct SwaKeys {
    reference: &'static str,
    alt: &'static str,
    ref_comp: &'static str,
    diff: &'static str,
    ref_seq: &'static str,
    ref_frame: &'static str,
    ref_context: &'static str,
    alt_seq: &'static str,
    alt_frame: &'static str,
    alt_context: &'static str,
    ref_comp_seq: &'static str,
}

/// The donor or acceptor side of the sliding window: its context flank, k-mer
/// length, scoring function and fields.
struct SwaSite {
    flank: i64,
    k: usize,
    score: Score,
    keys: SwaKeys,
}

const SWA_DONOR: SwaSite = SwaSite {
    flank: 8,
    k: 9,
    score: Tables::score5,
    keys: SwaKeys {
        reference: "MES-SWA_donor_ref",
        alt: "MES-SWA_donor_alt",
        ref_comp: "MES-SWA_donor_ref_comp",
        diff: "MES-SWA_donor_diff",
        ref_seq: "MES-SWA_donor_ref_seq",
        ref_frame: "MES-SWA_donor_ref_frame",
        ref_context: "MES-SWA_donor_ref_context",
        alt_seq: "MES-SWA_donor_alt_seq",
        alt_frame: "MES-SWA_donor_alt_frame",
        alt_context: "MES-SWA_donor_alt_context",
        ref_comp_seq: "MES-SWA_donor_ref_comp_seq",
    },
};
const SWA_ACCEPTOR: SwaSite = SwaSite {
    flank: 22,
    k: 23,
    score: Tables::score3,
    keys: SwaKeys {
        reference: "MES-SWA_acceptor_ref",
        alt: "MES-SWA_acceptor_alt",
        ref_comp: "MES-SWA_acceptor_ref_comp",
        diff: "MES-SWA_acceptor_diff",
        ref_seq: "MES-SWA_acceptor_ref_seq",
        ref_frame: "MES-SWA_acceptor_ref_frame",
        ref_context: "MES-SWA_acceptor_ref_context",
        alt_seq: "MES-SWA_acceptor_alt_seq",
        alt_frame: "MES-SWA_acceptor_alt_frame",
        alt_context: "MES-SWA_acceptor_alt_context",
        ref_comp_seq: "MES-SWA_acceptor_ref_comp_seq",
    },
};

/// The nearest canonical sites' score fields and their `_seq` fields.
const NCSS_UPSTREAM_DONOR: (&str, &str) =
    ("MES-NCSS_upstream_donor", "MES-NCSS_upstream_donor_seq");
const NCSS_UPSTREAM_ACCEPTOR: (&str, &str) = (
    "MES-NCSS_upstream_acceptor",
    "MES-NCSS_upstream_acceptor_seq",
);
const NCSS_DOWNSTREAM_DONOR: (&str, &str) =
    ("MES-NCSS_downstream_donor", "MES-NCSS_downstream_donor_seq");
const NCSS_DOWNSTREAM_ACCEPTOR: (&str, &str) = (
    "MES-NCSS_downstream_acceptor",
    "MES-NCSS_downstream_acceptor_seq",
);

/// The feature at 1-based `ordinal` in transcript order (5' to 3', the order
/// `get_all_Exons` and `get_all_Introns` return and the `EXON` and `INTRON`
/// columns count), among features listed in either genomic order; `None` past
/// the list. The list is read in place: a list in neither order is searched by
/// rank of start.
fn nth_in_transcript_order<T>(
    features: &[T],
    strand: i8,
    ordinal: usize,
    start: impl Fn(&T) -> u64,
) -> Option<&T> {
    let n = features.len();
    if ordinal == 0 || ordinal > n {
        return None;
    }
    let ascending = features.windows(2).all(|w| start(&w[0]) <= start(&w[1]));
    let descending = features.windows(2).all(|w| start(&w[0]) >= start(&w[1]));
    let from_low = if strand > 0 { ordinal - 1 } else { n - ordinal };
    if ascending {
        return features.get(from_low);
    }
    if descending {
        return features.get(n - 1 - from_low);
    }
    features
        .iter()
        .find(|f| features.iter().filter(|g| start(g) < start(f)).count() == from_low)
}

/// The first ordinal of an `EXON` or `INTRON` value (`3/5`, `2-3/5`) and the total.
fn first_ordinal(value: &str) -> Option<(usize, usize)> {
    let (ordinals, total) = value.split_once('/')?;
    let first = ordinals.split('-').next()?.trim().parse().ok()?;
    Some((first, total.trim().parse().ok()?))
}

/// The key of a row's nearest canonical sites: the transcript, whether the row
/// lies in an intron (`INTRON` set) or an exon, and the first ordinal; `None` for
/// a row in neither.
fn ncss_key(tc: &TranscriptConsequence) -> Option<(Arc<str>, bool, usize)> {
    let id = tc
        .transcript_id_handle()
        .map_or_else(|| Arc::from(""), Arc::clone);
    if let Some((number, _)) = tc.intron.as_deref().and_then(first_ordinal) {
        return Some((id, true, number));
    }
    let (number, _) = tc.exon.as_deref().and_then(first_ordinal)?;
    Some((id, false, number))
}

/// A canonical window a substitution lies in: its place in transcript order (the
/// intron's rank, the donor before the acceptor), its genomic span and its model.
struct MesWindow {
    order: (i64, u8),
    span: (i64, i64),
    score: Score,
}

/// One variant read on one transcript strand: what every group of a row reads.
struct Row<'a> {
    tables: &'a Tables,
    fasta: &'a vep_fasta::IndexedFasta,
    window: &'a Window,
    variant: &'a InputVariant,
    strand: i8,
    /// `feature_seq`: the row's alternate allele on the transcript's strand, the
    /// empty sequence for a deletion (`seq_length` 0, MaxEntScan.pm:657).
    feature_seq: &'a [u8],
    verbose: bool,
}

impl Row<'_> {
    /// Whether the row is a single-base substitution, the one shape `run_MES`
    /// scores (`:260`) and the sliding window compares in the alternate k-mer's
    /// frame (`:366`, `:404`).
    fn is_snv(&self) -> bool {
        self.variant.start == self.variant.end
            && self.feature_seq.len() == 1
            && is_acgt(self.feature_seq)
    }

    /// `get_seqs` (MaxEntScan.pm:642-663): the reference span on the transcript
    /// strand and the same span with the allele substituted at the variant's
    /// offset, `None` when the span leaves the chromosome.
    fn seqs(&self, start: i64, end: i64) -> Option<(Vec<u8>, Vec<u8>)> {
        let forward = self.window.forward(start, end)?;
        let reference = if self.strand < 0 {
            reverse_complement(forward)
        } else {
            forward.to_vec()
        };
        let (vf_start, vf_end) = (self.variant.start as i64, self.variant.end as i64);
        let offset = if self.strand > 0 {
            vf_start - start
        } else {
            end - vf_end
        };
        let offset = offset.clamp(0, reference.len() as i64) as usize;
        let after = (offset + (vf_end - vf_start + 1).max(0) as usize).min(reference.len());
        let mut alternate = Vec::with_capacity(reference.len() + self.feature_seq.len());
        alternate.extend_from_slice(&reference[..offset]);
        alternate.extend_from_slice(self.feature_seq);
        alternate.extend_from_slice(&reference[after..]);
        Some((reference, alternate))
    }

    /// The reference bases of the closed span on the transcript's strand, read
    /// from the FASTA (a canonical site may lie far outside the window).
    fn read(&self, start: i64, end: i64) -> Option<Vec<u8>> {
        if start < 1 {
            return None;
        }
        let bases = self
            .fasta
            .sequence(&self.variant.chr, start as u64, end as u64)?;
        Some(if self.strand < 0 {
            reverse_complement(&bases)
        } else {
            bases
        })
    }

    fn score(&self, fields: &mut Fields, key: &'static str, score: f64) {
        fields.push((key, format_score(score)));
    }

    /// A sequence field, kept under `verbose` alone.
    fn text(&self, fields: &mut Fields, key: &'static str, seq: &[u8]) {
        if self.verbose {
            fields.push((key, String::from_utf8_lossy(seq).into_owned()));
        }
    }

    fn frame(&self, fields: &mut Fields, key: &'static str, frame: usize) {
        if self.verbose {
            fields.push((key, frame.to_string()));
        }
    }

    /// Whether `run_MES` walks `intron` (MaxEntScan.pm:270): the list it walks
    /// is `_overlapped_introns`, memoised on its first call
    /// (`BaseTranscriptVariation.pm:802-822`). For a variant overlapping the
    /// transcript the consequence predicates make that call first, with the
    /// variant's own span (`BaseVariationFeatureOverlapAllele.pm:480-500`), so
    /// the plugin's 21-base window widens nothing and the list holds the
    /// introns whose span widened by three bases on each side
    /// (`BaseTranscriptVariation.pm:1027`) overlaps the variant; an intron
    /// with its start past its end is in no list (`:1014`). For a variant
    /// outside the transcript the plugin's own call is the first, and its
    /// window holds every intron a canonical window could reach.
    fn on_overlapped_intron_list(&self, intron: &Intron, within_feature: bool) -> bool {
        let (start, end) = (intron.start as i64, intron.end as i64);
        if start > end {
            return false;
        }
        if !within_feature {
            return true;
        }
        let (vf_start, vf_end) = (self.variant.start as i64, self.variant.end as i64);
        overlap(vf_start, vf_end, start - 3, end + 3)
    }

    /// `run_MES` (MaxEntScan.pm:256-332): the first canonical donor or acceptor
    /// window, in transcript order with the donor before the acceptor of one
    /// intron, that the substitution lies in, scored on the reference and the
    /// alternate sequence; a window the reference cannot be read for, or one
    /// holding a base outside ACGT, yields nothing. The introns walked are
    /// the ones on Ensembl VEP's overlapped-intron list
    /// ([`Self::on_overlapped_intron_list`]).
    fn mes(&self, ctx: &LofteeContext, fields: &mut Fields) {
        if !self.is_snv() {
            return;
        }
        let (vf_start, vf_end) = (self.variant.start as i64, self.variant.end as i64);
        let strand = self.strand;
        // `within_feature` (BaseVariationFeatureOverlapAllele.pm:480): the
        // variant overlaps the transcript's span, the extent of its exons.
        let transcript_start = ctx.exons.iter().map(|exon| exon.start as i64).min();
        let transcript_end = ctx.exons.iter().map(|exon| exon.end as i64).max();
        let within_feature = transcript_start
            .zip(transcript_end)
            .is_some_and(|(tr_start, tr_end)| overlap(vf_start, vf_end, tr_start, tr_end));
        let mut hit: Option<MesWindow> = None;
        for intron in ctx
            .introns
            .iter()
            .filter(|intron| self.on_overlapped_intron_list(intron, within_feature))
        {
            let start = intron.start as i64;
            let rank = if strand > 0 { start } else { -start };
            let donor = donor_of_intron(intron, strand);
            let acceptor = acceptor_of_intron(intron, strand);
            for (side, span, score) in [
                (0u8, donor, Tables::score5 as Score),
                (1u8, acceptor, Tables::score3 as Score),
            ] {
                let candidate = MesWindow {
                    order: (rank, side),
                    span,
                    score,
                };
                if overlap(vf_start, vf_end, span.0, span.1)
                    && hit.as_ref().is_none_or(|h| candidate.order < h.order)
                {
                    hit = Some(candidate);
                }
            }
        }
        let Some(MesWindow { span, score, .. }) = hit else {
            return;
        };
        let Some((reference, alternate)) = self.seqs(span.0, span.1) else {
            return;
        };
        if !is_acgt(&reference) || !is_acgt(&alternate) {
            return;
        }
        let (ref_score, alt_score) = (
            score(self.tables, &reference),
            score(self.tables, &alternate),
        );
        self.score(fields, "MaxEntScan_ref", ref_score);
        self.text(fields, "MaxEntScan_ref_seq", &reference);
        self.score(fields, "MaxEntScan_alt", alt_score);
        self.text(fields, "MaxEntScan_alt_seq", &alternate);
        self.score(fields, "MaxEntScan_diff", ref_score - alt_score);
    }

    /// `run_SWA` (MaxEntScan.pm:334-418) for one site type: the highest-scoring
    /// k-mer of the reference and of the alternate context, the reference
    /// comparison (the reference context in the alternate k-mer's frame for a
    /// substitution, else the reference's own best k-mer) and the difference.
    /// The results depend on the variant and the strand alone, so one
    /// computation serves every row of the variant on that strand.
    fn swa(&self, site: &SwaSite, fields: &mut Fields) {
        let (vf_start, vf_end) = (self.variant.start as i64, self.variant.end as i64);
        let Some((ref_context, alt_context)) =
            self.seqs(vf_start - site.flank, vf_end + site.flank)
        else {
            return;
        };
        let keys = &site.keys;
        self.text(fields, keys.ref_context, &ref_context);
        let mut ref_best = None;
        if is_acgt(&ref_context) {
            let (frame, score) = self.tables.best_kmer(&ref_context, site.k, site.score);
            self.text(
                fields,
                keys.ref_seq,
                &ref_context[frame - 1..frame - 1 + site.k],
            );
            self.frame(fields, keys.ref_frame, frame);
            self.score(fields, keys.reference, score);
            ref_best = Some(frame);
        }
        self.text(fields, keys.alt_context, &alt_context);
        if !is_acgt(&alt_context) {
            return;
        }
        let (frame, score) = self.tables.best_kmer(&alt_context, site.k, site.score);
        self.text(
            fields,
            keys.alt_seq,
            &alt_context[frame - 1..frame - 1 + site.k],
        );
        self.frame(fields, keys.alt_frame, frame);
        self.score(fields, keys.alt, score);
        let Some(ref_frame) = ref_best else {
            return;
        };
        let comparison_frame = if self.is_snv() { frame } else { ref_frame };
        let comparison = &ref_context[comparison_frame - 1..comparison_frame - 1 + site.k];
        let comparison_score = (site.score)(self.tables, comparison);
        self.text(fields, keys.ref_comp_seq, comparison);
        self.score(fields, keys.ref_comp, comparison_score);
        self.score(fields, keys.diff, comparison_score - score);
    }

    /// One nearest canonical site: its score when the sequence is all ACGT, and
    /// the sequence itself (`:438-446`).
    fn site(
        &self,
        fields: &mut Fields,
        keys: (&'static str, &'static str),
        span: (i64, i64),
        score: Score,
    ) {
        let Some(seq) = self.read(span.0, span.1) else {
            return;
        };
        if is_acgt(&seq) {
            self.score(fields, keys.0, score(self.tables, &seq));
        }
        self.text(fields, keys.1, &seq);
    }

    /// `run_NCSS` (MaxEntScan.pm:420-511): the canonical donor and acceptor
    /// sequences nearest the row's intron or exon, in transcript order, read from
    /// the reference on the transcript's strand. The results depend on the
    /// transcript and the ordinal alone ([`ncss_key`]), so one computation serves
    /// every row of a batch with that key.
    fn ncss(&self, tc: &TranscriptConsequence, ctx: &LofteeContext, fields: &mut Fields) {
        let strand = self.strand;
        if let Some((number, total)) = tc.intron.as_deref().and_then(first_ordinal) {
            let intron = |ordinal: usize| {
                nth_in_transcript_order(&ctx.introns, strand, ordinal, |i: &Intron| i.start)
            };
            let Some(this) = intron(number) else {
                return;
            };
            self.site(
                fields,
                NCSS_UPSTREAM_DONOR,
                donor_of_intron(this, strand),
                Tables::score5,
            );
            self.site(
                fields,
                NCSS_DOWNSTREAM_ACCEPTOR,
                acceptor_of_intron(this, strand),
                Tables::score3,
            );
            if number != 1 {
                if let Some(upstream) = intron(number - 1) {
                    self.site(
                        fields,
                        NCSS_UPSTREAM_ACCEPTOR,
                        acceptor_of_intron(upstream, strand),
                        Tables::score3,
                    );
                }
            }
            if number != total {
                if let Some(downstream) = intron(number + 1) {
                    self.site(
                        fields,
                        NCSS_DOWNSTREAM_DONOR,
                        donor_of_intron(downstream, strand),
                        Tables::score5,
                    );
                }
            }
        } else if let Some((number, total)) = tc.exon.as_deref().and_then(first_ordinal) {
            let exon = |ordinal: usize| {
                nth_in_transcript_order(&ctx.exons, strand, ordinal, |e: &Exon| e.start)
            };
            let Some(this) = exon(number) else {
                return;
            };
            if number != 1 {
                if let Some(upstream) = exon(number - 1) {
                    self.site(
                        fields,
                        NCSS_UPSTREAM_DONOR,
                        donor_of_exon(upstream, strand),
                        Tables::score5,
                    );
                }
                self.site(
                    fields,
                    NCSS_UPSTREAM_ACCEPTOR,
                    acceptor_of_exon(this, strand),
                    Tables::score3,
                );
            }
            if number != total {
                self.site(
                    fields,
                    NCSS_DOWNSTREAM_DONOR,
                    donor_of_exon(this, strand),
                    Tables::score5,
                );
                if let Some(downstream) = exon(number + 1) {
                    self.site(
                        fields,
                        NCSS_DOWNSTREAM_ACCEPTOR,
                        acceptor_of_exon(downstream, strand),
                        Tables::score3,
                    );
                }
            }
        }
    }
}

/// MaxEntScan maximum-entropy splice-site scores.
pub struct MaxEntScanPlugin {
    tables: Option<Tables>,
    swa: bool,
    ncss: bool,
    verbose: bool,
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
}

impl MaxEntScanPlugin {
    /// Create a new uninitialized MaxEntScan plugin.
    pub fn new() -> Self {
        Self {
            tables: None,
            swa: false,
            ncss: false,
            verbose: false,
            reference_fasta: None,
        }
    }

    /// `feature_seq` of a row: the alternate allele on the forward strand,
    /// reverse-complemented for a transcript on the reverse strand when it is DNA;
    /// the empty sequence for a deletion.
    fn feature_seq(variant: &InputVariant, strand: i8) -> Vec<u8> {
        let (_, alt) = variant.forward_strand_alleles();
        if alt.as_ref() == b"-" {
            return Vec::new();
        }
        if strand < 0 && seq_is_dna(&alt) {
            reverse_complement(&alt)
        } else {
            alt.into_owned()
        }
    }

    fn header_fields(&self) -> Vec<(&'static str, &'static str)> {
        let mut fields = vec![
            ("MaxEntScan_ref", "MaxEntScan reference sequence score"),
            ("MaxEntScan_alt", "MaxEntScan alternate sequence score"),
            ("MaxEntScan_diff", "MaxEntScan score difference"),
        ];
        if self.verbose {
            fields.extend([
                ("MaxEntScan_ref_seq", "MaxEntScan reference sequence"),
                ("MaxEntScan_alt_seq", "MaxEntScan alternate sequence"),
            ]);
        }
        if self.swa {
            fields.extend([
                ("MES-SWA_donor_ref", "Highest splice donor reference sequence score"),
                ("MES-SWA_donor_alt", "Highest splice donor alternate sequence score"),
                ("MES-SWA_donor_ref_comp", "Donor reference comparison sequence score"),
                (
                    "MES-SWA_donor_diff",
                    "Difference between the donor reference comparison and alternate sequence scores",
                ),
                ("MES-SWA_acceptor_ref", "Highest splice acceptor reference sequence score"),
                ("MES-SWA_acceptor_alt", "Highest splice acceptor alternate sequence score"),
                ("MES-SWA_acceptor_ref_comp", "Acceptor reference comparison sequence score"),
                (
                    "MES-SWA_acceptor_diff",
                    "Difference between the acceptor reference comparison and alternate sequence scores",
                ),
            ]);
            if self.verbose {
                fields.extend([
                    (
                        "MES-SWA_donor_ref_seq",
                        "Highest splice donor reference sequence",
                    ),
                    (
                        "MES-SWA_donor_ref_frame",
                        "Position of the highest splice donor reference sequence",
                    ),
                    (
                        "MES-SWA_donor_ref_context",
                        "Selected donor sequence context containing the reference allele",
                    ),
                    (
                        "MES-SWA_donor_alt_seq",
                        "Highest splice donor alternate sequence",
                    ),
                    (
                        "MES-SWA_donor_alt_frame",
                        "Position of the highest splice donor alternate sequence",
                    ),
                    (
                        "MES-SWA_donor_alt_context",
                        "Selected donor sequence context containing the alternate allele",
                    ),
                    (
                        "MES-SWA_donor_ref_comp_seq",
                        "Donor reference comparison sequence",
                    ),
                    (
                        "MES-SWA_acceptor_ref_seq",
                        "Highest splice acceptor reference sequence",
                    ),
                    (
                        "MES-SWA_acceptor_ref_frame",
                        "Position of the highest splice acceptor reference sequence",
                    ),
                    (
                        "MES-SWA_acceptor_ref_context",
                        "Selected acceptor sequence context containing the reference allele",
                    ),
                    (
                        "MES-SWA_acceptor_alt_seq",
                        "Highest splice acceptor alternate sequence",
                    ),
                    (
                        "MES-SWA_acceptor_alt_frame",
                        "Position of the highest splice acceptor alternate sequence",
                    ),
                    (
                        "MES-SWA_acceptor_alt_context",
                        "Selected acceptor sequence context containing the alternate allele",
                    ),
                    (
                        "MES-SWA_acceptor_ref_comp_seq",
                        "Acceptor reference comparison sequence",
                    ),
                ]);
            }
        }
        if self.ncss {
            fields.extend([
                (
                    "MES-NCSS_upstream_acceptor",
                    "Nearest upstream canonical splice acceptor sequence score",
                ),
                (
                    "MES-NCSS_upstream_donor",
                    "Nearest upstream canonical splice donor sequence score",
                ),
                (
                    "MES-NCSS_downstream_acceptor",
                    "Nearest downstream canonical splice acceptor sequence score",
                ),
                (
                    "MES-NCSS_downstream_donor",
                    "Nearest downstream canonical splice donor sequence score",
                ),
            ]);
            if self.verbose {
                fields.extend([
                    (
                        "MES-NCSS_upstream_acceptor_seq",
                        "Nearest upstream canonical splice acceptor sequence",
                    ),
                    (
                        "MES-NCSS_upstream_donor_seq",
                        "Nearest upstream canonical splice donor sequence",
                    ),
                    (
                        "MES-NCSS_downstream_acceptor_seq",
                        "Nearest downstream canonical splice acceptor sequence",
                    ),
                    (
                        "MES-NCSS_downstream_donor_seq",
                        "Nearest downstream canonical splice donor sequence",
                    ),
                ]);
            }
        }
        fields
    }
}

impl Default for MaxEntScanPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for MaxEntScanPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The fields of the groups switched on (`get_header_info`, MaxEntScan.pm:124-219).
    fn header_info(&self) -> Vec<(String, String)> {
        self.header_fields()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    /// The first positional parameter is the matrix directory, which must exist
    /// (MaxEntScan.pm:98-100); the bare tokens `SWA`, `NCSS` and `verbose` switch
    /// the extra groups on (`:110-115`). Any other token is not read, as the
    /// Ensembl VEP plugin looks the three up in a hash of the tokens and never
    /// at the rest; it is logged.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &[], params)?;
        let mut positional = params.positional().iter();
        let dir = positional
            .next()
            .ok_or_else(|| PluginError::Init(format!("{NAME}: directory not specified")))?;
        let dir = Path::new(dir);
        if !dir.is_dir() {
            return Err(PluginError::Init(format!(
                "{NAME}: directory not found: '{}'",
                dir.display()
            )));
        }
        for token in positional {
            match token.as_str() {
                "SWA" => self.swa = true,
                "NCSS" => self.ncss = true,
                "verbose" => self.verbose = true,
                _ => warn_unread_parameters(
                    NAME,
                    std::slice::from_ref(token),
                    "the matrix directory and the switches SWA, NCSS and verbose",
                ),
            }
        }
        self.tables = Some(Tables::load(dir)?);
        Ok(())
    }

    fn requires_reference_fasta(&self) -> bool {
        true
    }

    fn set_reference_fasta(&mut self, fasta: Option<Arc<vep_fasta::IndexedFasta>>) {
        self.reference_fasta = fasta;
    }

    /// `run` (MaxEntScan.pm:221-254) over the batch: one reference window per
    /// variant, the sliding-window fields once per variant and strand, the
    /// nearest-site fields once per transcript and ordinal, the canonical-window
    /// fields per row.
    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        let (Some(tables), Some(fasta)) = (self.tables.as_ref(), self.reference_fasta.as_ref())
        else {
            return Ok(());
        };
        let mut ncss_cache: HashMap<(Arc<str>, bool, usize), Fields> = HashMap::new();
        for variant in variants.iter_mut() {
            // The plugin declares `VariationFeature` rows alone (`variant_feature_types`,
            // BaseVepPlugin.pm:100), so Ensembl VEP never runs it on a structural variant.
            if variant.is_structural
                || !variant
                    .transcript_consequences
                    .iter()
                    .any(|tc| self.annotates_row(tc))
            {
                continue;
            }
            let chr_start = variant
                .start
                .min(variant.end)
                .saturating_sub(WINDOW_FLANK)
                .max(1);
            let chr_end = variant.start.max(variant.end) + WINDOW_FLANK;
            let Some(bases) = fasta
                .sequence(&variant.chr, chr_start, chr_end)
                .or_else(|| {
                    let length = fasta.sequence_length(&variant.chr)?;
                    fasta.sequence(&variant.chr, chr_start, chr_end.min(length))
                })
            else {
                continue;
            };
            let window = Window {
                chr_start: chr_start as i64,
                bases,
            };
            let feature_seqs = [
                Self::feature_seq(variant, 1),
                Self::feature_seq(variant, -1),
            ];
            let mut swa_cache: [Option<Fields>; 2] = [None, None];
            for tc_idx in 0..variant.transcript_consequences.len() {
                let tc = &variant.transcript_consequences[tc_idx];
                if !self.annotates_row(tc) {
                    continue;
                }
                let Some(ctx) = tc.loftee_ctx.as_deref() else {
                    continue;
                };
                let side = usize::from(ctx.strand < 0);
                let row = Row {
                    tables,
                    fasta,
                    window: &window,
                    variant,
                    strand: ctx.strand,
                    feature_seq: &feature_seqs[side],
                    verbose: self.verbose,
                };
                let mut fields = Fields::new();
                row.mes(ctx, &mut fields);
                if self.swa {
                    let swa = swa_cache[side].get_or_insert_with(|| {
                        let mut swa = Fields::new();
                        row.swa(&SWA_DONOR, &mut swa);
                        row.swa(&SWA_ACCEPTOR, &mut swa);
                        swa
                    });
                    fields.extend(swa.iter().cloned());
                }
                if self.ncss {
                    if let Some(key) = ncss_key(tc) {
                        let ncss = ncss_cache.entry(key).or_insert_with(|| {
                            let mut ncss = Fields::new();
                            row.ncss(tc, ctx, &mut ncss);
                            ncss
                        });
                        fields.extend(ncss.iter().cloned());
                    }
                }
                let plugin_data = &mut variant.transcript_consequences[tc_idx].plugin_data;
                plugin_data.reserve(fields.len());
                for (key, value) in fields {
                    plugin_data.insert(key.to_string(), value);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::IndexMap;
    use std::path::PathBuf;
    use std::process::{Command, Stdio};
    use std::sync::OnceLock;
    use vep_core::consequence::FeatureType;

    /// The seed-1 matrices of `scripts/fixtures/maxentscan_fixture_matrices.py`,
    /// written once per test process. Every expected score below is that
    /// generator's own `score5` or `score3` of the named k-mer on this seed
    /// (`--dir <dir> --score <kmer>`), and the k-mers with their scores are the
    /// ones the Ensembl VEP plugin printed over the same matrices for the
    /// GRCh38 records the sequences below are copied from. The generator
    /// draws a spelling for every number line and writes files beside the
    /// eleven tables, so each test over this directory reads those too.
    fn fixture_dir() -> PathBuf {
        static DIR: OnceLock<tempfile::TempDir> = OnceLock::new();
        DIR.get_or_init(|| write_fixture(&[]))
            .path()
            .join("synthetic_dir")
    }

    /// The same seed under `--plain-lines`: every number as its five-place
    /// decimal with nothing around it, and the eleven tables alone.
    fn plain_fixture_dir() -> PathBuf {
        static DIR: OnceLock<tempfile::TempDir> = OnceLock::new();
        DIR.get_or_init(|| write_fixture(&["--plain-lines"]))
            .path()
            .join("synthetic_dir")
    }

    fn write_fixture(options: &[&str]) -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        let script = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../scripts/fixtures/maxentscan_fixture_matrices.py");
        let status = Command::new("python3")
            .arg(script)
            .args(["--seed", "1"])
            .args(options)
            .arg("--out")
            .arg(dir.path().join("synthetic_dir"))
            .stdout(Stdio::null())
            .status()
            .expect("python3 runs the matrix generator");
        assert!(status.success(), "the matrix generator writes seed 1");
        dir
    }

    fn init_message(result: Result<(), PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    fn plugin(tokens: &[&str]) -> MaxEntScanPlugin {
        let mut params = vec![fixture_dir().display().to_string()];
        params.extend(tokens.iter().map(|t| t.to_string()));
        let mut p = MaxEntScanPlugin::new();
        p.init(&params).unwrap();
        p
    }

    fn tables() -> Tables {
        Tables::load(&fixture_dir()).unwrap()
    }

    fn score_text(score: f64) -> String {
        format!("{score:.3}")
    }

    /// A deterministic filler base for position `pos`, so k-mers of the filler
    /// are distinct and never form the splice motifs the layouts place.
    fn filler(pos: usize) -> u8 {
        b"ACGT"[(pos * 7 + pos / 3) % 4]
    }

    /// A FASTA of one record `name` of `length` bases: the filler everywhere
    /// but the `(start, bases)` spans, written with its index.
    fn fasta(
        dir: &Path,
        name: &str,
        length: usize,
        spans: &[(usize, &[u8])],
    ) -> Arc<vep_fasta::IndexedFasta> {
        let mut sequence: Vec<u8> = (1..=length).map(filler).collect();
        for (start, bases) in spans {
            sequence[start - 1..start - 1 + bases.len()].copy_from_slice(bases);
        }
        let fa = dir.join(format!("{name}.fa"));
        let mut contents = format!(">{name}\n").into_bytes();
        contents.extend_from_slice(&sequence);
        contents.push(b'\n');
        std::fs::write(&fa, contents).unwrap();
        std::fs::write(
            format!("{}.fai", fa.display()),
            format!(
                "{name}\t{length}\t{}\t{length}\t{}\n",
                name.len() + 2,
                length + 1
            ),
        )
        .unwrap();
        Arc::new(vep_fasta::IndexedFasta::from_path(&fa).unwrap())
    }

    /// A transcript context over the `(start, end)` exon spans, given in any
    /// order: the exons ranked in transcript order, and the introns as the gaps
    /// between the exons in genomic order, ranked in that order, the shape the
    /// cache loader gives a transcript without its own list.
    fn context(strand: i8, exons: &[(u64, u64)]) -> LofteeContext {
        let mut by_position = exons.to_vec();
        by_position.sort_unstable();
        let introns: Vec<Intron> = by_position
            .windows(2)
            .enumerate()
            .map(|(i, pair)| Intron {
                start: pair[0].1 + 1,
                end: pair[1].0 - 1,
                rank: i as u32 + 1,
            })
            .collect();
        let mut in_transcript_order = by_position;
        if strand < 0 {
            in_transcript_order.reverse();
        }
        let exons: Vec<Exon> = in_transcript_order
            .iter()
            .enumerate()
            .map(|(i, (start, end))| Exon {
                start: *start,
                end: *end,
                rank: i as u32 + 1,
                phase: -1,
                end_phase: -1,
            })
            .collect();
        LofteeContext {
            introns,
            exons,
            strand,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            cds_start_nf: false,
            cds_end_nf: false,
        }
    }

    /// A transcript row of `ctx` with the `EXON` and `INTRON` values given.
    fn row(ctx: &LofteeContext, exon: Option<&str>, intron: Option<&str>) -> TranscriptConsequence {
        TranscriptConsequence {
            feature_type: FeatureType::Transcript,
            strand: ctx.strand,
            exon: exon.map(str::to_string),
            intron: intron.map(str::to_string),
            loftee_ctx: Some(Box::new(ctx.clone())),
            ..TranscriptConsequence::default()
        }
    }

    fn variant(
        chr: &str,
        start: u64,
        end: u64,
        reference: &[u8],
        alt: &[u8],
        tc: TranscriptConsequence,
    ) -> InputVariant {
        let mut v = InputVariant::new(chr.into(), start, end, reference.to_vec(), alt.to_vec());
        v.transcript_consequences.push(tc);
        v
    }

    fn annotate(p: &MaxEntScanPlugin, mut v: InputVariant) -> IndexMap<String, String> {
        p.run_batch(std::slice::from_mut(&mut v)).unwrap();
        v.transcript_consequences.remove(0).plugin_data
    }

    /// The forward-strand layout: chromosome `21` of 300 bases carrying, at the
    /// coordinates the tests name, the sequence of the GRCh38 records
    /// 21:14216282 (exon 1's last base; the 45-base acceptor context at
    /// 101-145), 21:14219561 and 21:14219563 (the acceptor of intron 1 at
    /// 201-223), the donor of exon 2 at 258-266 and the acceptor of exon 3 at
    /// 270-292. Exons (60,123), (221,260), (290,300).
    fn forward_layout(dir: &Path) -> (Arc<vep_fasta::IndexedFasta>, LofteeContext) {
        let fa = fasta(
            dir,
            "21",
            300,
            &[
                (101, b"TTCTGTACGAGCTGTTCCTTCAGGTACCGTCTCTGGGAAGGGGCG"),
                (201, b"AATTTTTAAGTTTCTTATAGGCG"),
                (258, b"TTGGTAGGT"),
                (270, b"TTACTCTGTTTCTTTCAAAGGGA"),
            ],
        );
        (fa, context(1, &[(60, 123), (221, 260), (290, 300)]))
    }

    /// The reverse-strand layout: chromosome `22` of 300 bases carrying the
    /// reverse complement of the sequence the Ensembl VEP plugin reads on the
    /// transcript strand of the GRCh38 records 21:14383095 and 21:14383093 (the
    /// donor of intron 1, 45 bases at 201-245), 21:14381545 and 21:14381543
    /// (the acceptor of exon 2 at 137-159), the donor of exon 2 at 94-102 and
    /// the acceptor of exon 3 at 58-80. Exons in transcript order (223,260),
    /// (100,139), (40,60).
    fn reverse_layout(dir: &Path) -> (Arc<vep_fasta::IndexedFasta>, LofteeContext) {
        let rc = reverse_complement;
        let fa = fasta(
            dir,
            "22",
            300,
            &[
                (201, &rc(b"GGCCAGAGAGATGACGATCTTAGGTGGGTTCCCGGGGGTCCTTGT")),
                (137, &rc(b"TTAATTTTCCTTTATTTCAGGAT")),
                (94, &rc(b"AAGGTAAGT")),
                (58, &rc(b"TTTGCAAATCCTATTCACAGGTT")),
            ],
        );
        (fa, context(-1, &[(223, 260), (100, 139), (40, 60)]))
    }

    #[test]
    fn directory_is_the_first_positional_parameter_and_must_exist() {
        let mut p = MaxEntScanPlugin::new();
        assert_eq!(
            init_message(p.init(&[])),
            "MaxEntScan: directory not specified"
        );
        let mut p = MaxEntScanPlugin::new();
        assert_eq!(
            init_message(p.init(&["/no/such/fordownload_synthetic".to_string()])),
            "MaxEntScan: directory not found: '/no/such/fordownload_synthetic'"
        );
        let dir = fixture_dir().display().to_string();
        let mut p = MaxEntScanPlugin::new();
        assert!(init_message(p.init(&[format!("dir={dir}")]))
            .starts_with("MaxEntScan: unknown parameter 'dir'"));
        let p = plugin(&["NCSS", "SWA", "verbose"]);
        assert!(p.swa && p.ncss && p.verbose && p.tables.is_some());
        let p = plugin(&[]);
        assert!(!p.swa && !p.ncss && !p.verbose);
    }

    /// A token after the directory that is none of the three switches is not
    /// read (`MaxEntScan.pm:110-115` looks the three up in a hash of the tokens
    /// and never at the rest): `swa` switches nothing on, and the switches
    /// around an unread token still count.
    #[test]
    fn a_token_that_is_no_switch_is_not_read() {
        let p = plugin(&["swa"]);
        assert!(!p.swa && !p.ncss && !p.verbose && p.tables.is_some());
        let p = plugin(&["SWA", "extra", "verbose"]);
        assert!(p.swa && !p.ncss && p.verbose);
    }

    #[test]
    fn matrices_with_wrong_line_counts_are_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let copy = dir.path().join("synthetic_dir");
        std::fs::create_dir_all(copy.join("splicemodels")).unwrap();
        let mut names = vec![
            "me2x5".to_string(),
            "splicemodels/splice5sequences".to_string(),
        ];
        names.extend(
            ACCEPTOR_TABLES
                .iter()
                .map(|(file, _, _, _)| format!("splicemodels/{file}")),
        );
        for name in &names {
            std::fs::copy(plain_fixture_dir().join(name), copy.join(name)).unwrap();
        }
        let short = std::fs::read_to_string(copy.join("splicemodels/me2x3acc6")).unwrap();
        let lines: Vec<&str> = short.lines().collect();
        std::fs::write(
            copy.join("splicemodels/me2x3acc6"),
            lines[..63].join("\n") + "\n",
        )
        .unwrap();
        let message = init_message(MaxEntScanPlugin::new().init(&[copy.display().to_string()]));
        assert!(
            message.ends_with("splicemodels/me2x3acc6 has 63 lines, expected 64"),
            "{message}"
        );

        std::fs::write(copy.join("splicemodels/me2x3acc6"), &short).unwrap();
        let sequences =
            std::fs::read_to_string(copy.join("splicemodels/splice5sequences")).unwrap();
        let mut lines: Vec<&str> = sequences.lines().collect();
        lines[1] = lines[0];
        std::fs::write(
            copy.join("splicemodels/splice5sequences"),
            lines.join("\n") + "\n",
        )
        .unwrap();
        let message = init_message(MaxEntScanPlugin::new().init(&[copy.display().to_string()]));
        assert!(
            message.ends_with("splice5sequences does not list every 7-mer over ACGT once"),
            "{message}"
        );

        std::fs::write(copy.join("splicemodels/splice5sequences"), &sequences).unwrap();
        let scores = std::fs::read_to_string(copy.join("me2x5")).unwrap();
        let mut lines: Vec<&str> = scores.lines().collect();
        lines[2] = "high";
        std::fs::write(copy.join("me2x5"), lines.join("\n") + "\n").unwrap();
        let message = init_message(MaxEntScanPlugin::new().init(&[copy.display().to_string()]));
        assert!(
            message.ends_with("me2x5 line 3 is not a number: 'high'"),
            "{message}"
        );
    }

    /// A byte that is not UTF-8 (Latin-1 `é`, 0xE9) is read as a byte of its
    /// line, as `<SCOREF>` reads it (MaxEntScan.pm:738, 754, 832): the line is
    /// then refused by the check on its content that every line meets, a
    /// number or a 7-mer over ACGT, with the byte written as U+FFFD. The one
    /// whitespace character `s/\s//` (`:740`) removes is ASCII on the module's
    /// byte string, so a line opening with U+00A0 keeps it and is no number.
    #[test]
    fn a_byte_that_is_not_utf8_is_refused_by_the_content_check() {
        let dir = tempfile::tempdir().unwrap();
        let copy = dir.path().join("synthetic_dir");
        std::fs::create_dir_all(copy.join("splicemodels")).unwrap();
        for (file, _, _, _) in ACCEPTOR_TABLES {
            let name = format!("splicemodels/{file}");
            std::fs::copy(plain_fixture_dir().join(&name), copy.join(&name)).unwrap();
        }
        let sequences_name = "splicemodels/splice5sequences";
        std::fs::copy(
            plain_fixture_dir().join(sequences_name),
            copy.join(sequences_name),
        )
        .unwrap();
        let scores = std::fs::read_to_string(plain_fixture_dir().join("me2x5")).unwrap();
        let mut bytes = scores.as_bytes().to_vec();
        bytes.insert(scores.match_indices('\n').nth(2).unwrap().0, 0xE9);
        std::fs::write(copy.join("me2x5"), &bytes).unwrap();
        let message = init_message(MaxEntScanPlugin::new().init(&[copy.display().to_string()]));
        assert!(
            message.ends_with(&format!(
                "me2x5 line 3 is not a number: '{}\u{FFFD}'",
                scores.lines().nth(2).unwrap()
            )),
            "{message}"
        );
        std::fs::write(copy.join("me2x5"), &scores).unwrap();
        let sequences = std::fs::read_to_string(plain_fixture_dir().join(sequences_name)).unwrap();
        let mut bytes = sequences.as_bytes().to_vec();
        bytes.insert(3, 0xE9);
        std::fs::write(copy.join(sequences_name), &bytes).unwrap();
        let message = init_message(MaxEntScanPlugin::new().init(&[copy.display().to_string()]));
        assert!(
            message.ends_with(&format!(
                "splice5sequences line 1 is not a 7-mer over ACGT: '{}\u{FFFD}{}'",
                &sequences[..3],
                &sequences[3..7]
            )),
            "{message}"
        );
        std::fs::write(copy.join(sequences_name), &sequences).unwrap();
        let mut text = scores.clone();
        text.insert(0, '\u{A0}');
        std::fs::write(copy.join("me2x5"), &text).unwrap();
        let message = init_message(MaxEntScanPlugin::new().init(&[copy.display().to_string()]));
        assert!(
            message.ends_with(&format!(
                "me2x5 line 1 is not a number: '\u{A0}{}'",
                scores.lines().next().unwrap()
            )),
            "{message}"
        );
    }

    /// A matrix line one byte over the cap is refused as the plugin opens,
    /// naming the file and the line, before the line count is checked.
    #[test]
    fn matrix_line_past_the_cap_is_refused() {
        let cap = crate::tabix::MAX_LINE_BYTES;
        let dir = tempfile::tempdir().unwrap();
        let copy = dir.path().join("synthetic_dir");
        std::fs::create_dir_all(&copy).unwrap();
        let mut text = std::fs::read(plain_fixture_dir().join("me2x5")).unwrap();
        text.extend(std::iter::repeat_n(b'x', cap + 1));
        text.push(b'\n');
        std::fs::write(copy.join("me2x5"), text).unwrap();
        let message = init_message(MaxEntScanPlugin::new().init(&[copy.display().to_string()]));
        assert_eq!(
            message,
            format!(
                "{}: line {} is longer than the {cap}-byte limit for one line",
                copy.join("me2x5").display(),
                DONOR_TABLE_LEN + 1
            )
        );
    }

    /// A one-table file of `lines`, read as the plugin's number tables are.
    fn numbers_of(lines: &[&str]) -> Result<Vec<f64>, PluginError> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("table");
        std::fs::write(&path, lines.join("\n") + "\n").unwrap();
        read_numbers(&path, lines.len())
    }

    /// Perl's numeric conversion skips the whitespace before and after a
    /// number (`$a * $d`, MaxEntScan.pm:689, over the text `s/\s//` leaves,
    /// :756-757), so a line between blanks is its number however many blanks
    /// stand there, and a carriage return before the line feed is one of them.
    #[test]
    fn blanks_around_a_number_are_not_part_of_it() {
        let lines = [
            "0.5",
            " 0.5",
            "  0.5",
            "0.5  ",
            " 0.5 ",
            "\t0.5",
            "\t\t0.5\t",
            "    0.5 \t",
            "\x0B\x0B0.5",
            "\x0C0.5\x0C",
        ];
        assert_eq!(numbers_of(&lines).unwrap(), vec![0.5; lines.len()]);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("table");
        std::fs::write(&path, "0.5\r\n  0.5\r\n0.5 \r\n").unwrap();
        assert_eq!(read_numbers(&path, 3).unwrap(), [0.5; 3]);
    }

    /// A number with an exponent is the number of its decimal spelling, under
    /// either sign of the exponent and either case of its letter.
    #[test]
    fn an_exponent_spelling_is_the_same_number() {
        let lines = [
            "0.00125",
            "1.25e-03",
            "1.25E-03",
            "125e-5",
            "0.000125e+01",
            "+1.25e-3",
            ".125e-2",
        ];
        assert_eq!(numbers_of(&lines).unwrap(), vec![0.00125; lines.len()]);
    }

    /// A line Perl reads only under a warning (an empty line, a numeric prefix
    /// before other text, a no-break space before the number) is refused by
    /// its table and line, quoting the text `s/\s//` leaves.
    #[test]
    fn a_line_perl_reads_under_a_warning_is_refused_by_its_line() {
        for text in [
            "",
            "high",
            "0.5x",
            "0.5 0.25",
            "0.5e",
            "0x10",
            "1_000",
            "0,5",
            "infinite",
            "\u{A0}0.5",
        ] {
            let message = init_message(numbers_of(&["0.5", text]).map(|_| ()));
            let kept = text.replacen(|c: char| c.is_ascii_whitespace() || c == '\x0B', "", 1);
            assert!(
                message.ends_with(&format!("table line 2 is not a number: '{kept}'")),
                "{text:?}: {message}"
            );
        }
    }

    /// Perl converts further spellings of an infinity or a NaN without a
    /// warning, and its own `0 but true`; `f64::from_str` reads none of them,
    /// and a line holding one is refused by its table and line like any text
    /// that is no number here.
    #[test]
    fn a_further_spelling_of_an_infinity_or_a_nan_is_refused_by_its_line() {
        for text in [
            "nanq",
            "NaNS",
            "qnan",
            "sNaN",
            "nan(1)",
            "nan(0x7ff)",
            "1.#INF",
            "-1.#INF",
            "1#INF",
            "1.#IND",
            "1.#QNAN",
            "0 but true",
        ] {
            let message = init_message(numbers_of(&["0.5", text]).map(|_| ()));
            let kept = text.replacen(|c: char| c.is_ascii_whitespace() || c == '\x0B', "", 1);
            assert!(
                message.ends_with(&format!("table line 2 is not a number: '{kept}'")),
                "{text:?}: {message}"
            );
        }
    }

    /// Perl converts the words for an infinity and a NaN, and a number past
    /// the largest double, without a warning, and `sprintf('%.3f')`
    /// (MaxEntScan.pm:245) prints a score made of one as `Inf` or `NaN`: a
    /// table of such entries is read, and the rows that read them carry those
    /// words.
    #[test]
    fn an_infinity_and_a_nan_are_read_and_printed_as_perl_prints_them() {
        let read = numbers_of(&["inf", "-Infinity", "1e400", "+INF", "nan", "-NaN"]).unwrap();
        assert_eq!(
            read[..4],
            [f64::INFINITY, -f64::INFINITY, f64::INFINITY, f64::INFINITY]
        );
        assert!(read[4].is_nan() && read[5].is_nan());

        let dir = tempfile::tempdir().unwrap();
        let (fa, ctx) = forward_layout(dir.path());
        let mut names = vec!["splicemodels/splice5sequences".to_string()];
        names.extend(
            ACCEPTOR_TABLES
                .iter()
                .map(|(file, _, _, _)| format!("splicemodels/{file}")),
        );
        for (entry, expected) in [("inf", ["Inf", "Inf", "NaN"]), ("nan", ["NaN"; 3])] {
            let copy = dir.path().join(entry).join("synthetic_dir");
            std::fs::create_dir_all(copy.join("splicemodels")).unwrap();
            for name in &names {
                std::fs::copy(plain_fixture_dir().join(name), copy.join(name)).unwrap();
            }
            std::fs::write(
                copy.join("me2x5"),
                format!("{entry}\n").repeat(DONOR_TABLE_LEN),
            )
            .unwrap();
            let mut p = MaxEntScanPlugin::new();
            p.init(&[copy.display().to_string()]).unwrap();
            p.set_reference_fasta(Some(fa.clone()));
            // The last exonic base of exon 1 lies in the donor window of intron 1.
            let out = annotate(
                &p,
                variant("21", 123, 123, b"G", b"A", row(&ctx, Some("1/3"), None)),
            );
            let written = ["MaxEntScan_ref", "MaxEntScan_alt", "MaxEntScan_diff"]
                .map(|key| out.get(key).cloned().unwrap_or_default());
            assert_eq!(written, expected, "{entry}");
        }
    }

    /// Every spelling the generator draws for a line (a decimal or an
    /// exponent of either sign, the lines of `me2x5` between spaces or tabs)
    /// reads as the number the plain directory holds, and the files written
    /// beside the eleven tables, two of them named as a table in the directory
    /// the plugin does not read it from, are never opened.
    #[test]
    fn every_spelling_of_a_table_reads_as_the_plain_table() {
        let me2x5 = std::fs::read_to_string(fixture_dir().join("me2x5")).unwrap();
        let blanks = |line: &str| line.len() - line.trim_matches([' ', '\t']).len();
        assert!(me2x5.lines().all(|line| blanks(line) >= 1));
        assert!(me2x5
            .lines()
            .any(|line| blanks(line) >= 2 && line.starts_with(' ')));
        assert!(me2x5.lines().any(|line| line.ends_with("  ")));
        assert!(me2x5.lines().any(|line| line.starts_with('\t')));
        for spelling in ["e-", "e+"] {
            assert!(me2x5.lines().any(|line| line.contains(spelling)));
        }
        assert!(me2x5.lines().any(|line| !line.contains('e')));
        for other in ["splicemodels/me2x5", "splice5sequences", "notes.txt"] {
            assert!(fixture_dir().join(other).is_file(), "{other}");
            assert!(!plain_fixture_dir().join(other).exists(), "{other}");
        }

        let (spelled, plain) = (tables(), Tables::load(&plain_fixture_dir()).unwrap());
        assert_eq!(spelled.donor, plain.donor);
        assert_eq!(spelled.acceptor, plain.acceptor);
    }

    #[test]
    fn score5_scores_known_9mers_on_the_made_up_tables() {
        let t = tables();
        for (kmer, expected) in [
            ("CAGGTACCG", "0.312"),
            ("CAAGTACCG", "6.055"),
            ("CAGGCACCG", "-7.442"),
            ("TTGGTAGGT", "1.725"),
            ("TCAGGTACC", "1.066"),
            ("TAGGTGGGT", "3.450"),
            ("TAAGTGGGT", "4.541"),
            ("AAGGTAAGT", "8.192"),
        ] {
            assert_eq!(score_text(t.score5(kmer.as_bytes())), expected, "{kmer}");
        }
    }

    #[test]
    fn score3_scores_known_23mers_on_the_made_up_tables() {
        let t = tables();
        for (kmer, expected) in [
            ("AATTTTTAAGTTTCTTATAGGCG", "-1.563"),
            ("AATTTTTAAGTTTCTTATGGGCG", "-9.518"),
            ("AATTTTTAAGTTTCTTATAGACG", "-3.332"),
            ("TGTACGAGCTGTTCCTTCAGGTA", "25.228"),
            ("AAGTACCGTCTCTGGGAAGGGGC", "22.130"),
            ("AGGTACCGTCTCTGGGAAGGGGC", "21.683"),
            ("TTACTCTGTTTCTTTCAAAGGGA", "10.637"),
            ("TTAATTTTCCTTTATTTCAGGAT", "7.343"),
            ("GAGATGACGATCTTAGGTGGGTT", "17.338"),
            ("TTTGCAAATCCTATTCACAGGTT", "-7.316"),
        ] {
            assert_eq!(score_text(t.score3(kmer.as_bytes())), expected, "{kmer}");
        }
    }

    #[test]
    fn hash_is_the_plugins_base_four_hash() {
        assert_eq!(hash(b"CAGAAGT"), 4619);
        assert_eq!(hash(b"AAA"), 0);
        assert_eq!(hash(b"TTTT"), 255);
    }

    #[test]
    fn mes_windows_on_both_strands() {
        let dir = tempfile::tempdir().unwrap();
        let mut p = plugin(&["verbose"]);
        let (fa, ctx) = forward_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let mes = |out: &IndexMap<String, String>| {
            [
                "MaxEntScan_ref_seq",
                "MaxEntScan_alt_seq",
                "MaxEntScan_ref",
                "MaxEntScan_alt",
                "MaxEntScan_diff",
            ]
            .map(|k| out.get(k).cloned().unwrap_or_default())
        };
        // The last exonic base of exon 1 lies in the donor window of intron 1.
        let out = annotate(
            &p,
            variant("21", 123, 123, b"G", b"A", row(&ctx, Some("1/3"), None)),
        );
        assert_eq!(
            mes(&out),
            ["CAGGTACCG", "CAAGTACCG", "0.312", "6.055", "-5.743"]
        );
        // The second intronic base of the donor.
        let out = annotate(
            &p,
            variant("21", 125, 125, b"T", b"C", row(&ctx, None, Some("1/2"))),
        );
        assert_eq!(
            mes(&out),
            ["CAGGTACCG", "CAGGCACCG", "0.312", "-7.442", "7.754"]
        );
        // The second-to-last intronic base of the acceptor, and the first exonic base after it.
        let out = annotate(
            &p,
            variant("21", 219, 219, b"A", b"G", row(&ctx, None, Some("1/2"))),
        );
        assert_eq!(
            mes(&out),
            [
                "AATTTTTAAGTTTCTTATAGGCG",
                "AATTTTTAAGTTTCTTATGGGCG",
                "-1.563",
                "-9.518",
                "7.955"
            ]
        );
        let out = annotate(
            &p,
            variant("21", 221, 221, b"G", b"A", row(&ctx, Some("2/3"), None)),
        );
        assert_eq!(
            mes(&out),
            [
                "AATTTTTAAGTTTCTTATAGGCG",
                "AATTTTTAAGTTTCTTATAGACG",
                "-1.563",
                "-3.332",
                "1.769"
            ]
        );
        // Eight bases into the intron no window reaches; a deletion is never scored.
        let out = annotate(
            &p,
            variant("21", 131, 131, b"G", b"A", row(&ctx, None, Some("1/2"))),
        );
        assert!(!out.contains_key("MaxEntScan_ref"), "{out:?}");
        let out = annotate(
            &p,
            variant("21", 124, 125, b"GT", b"-", row(&ctx, None, Some("1/2"))),
        );
        assert!(!out.contains_key("MaxEntScan_ref"), "{out:?}");

        // The reverse strand reads the windows at the introns' other ends, on the
        // transcript's strand.
        let mut p = plugin(&["verbose"]);
        let (fa, ctx) = reverse_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let out = annotate(
            &p,
            variant("22", 223, 223, b"C", b"T", row(&ctx, Some("1/3"), None)),
        );
        assert_eq!(
            mes(&out),
            ["TAGGTGGGT", "TAAGTGGGT", "3.450", "4.541", "-1.092"]
        );
        let out = annotate(
            &p,
            variant("22", 221, 221, b"A", b"G", row(&ctx, None, Some("1/2"))),
        );
        assert_eq!(
            mes(&out),
            ["TAGGTGGGT", "TAGGCGGGT", "3.450", "-4.304", "7.754"]
        );
        let out = annotate(
            &p,
            variant("22", 141, 141, b"T", b"C", row(&ctx, None, Some("1/2"))),
        );
        assert_eq!(
            mes(&out),
            [
                "TTAATTTTCCTTTATTTCAGGAT",
                "TTAATTTTCCTTTATTTCGGGAT",
                "7.343",
                "-0.612",
                "7.955"
            ]
        );
        let out = annotate(
            &p,
            variant("22", 139, 139, b"C", b"A", row(&ctx, Some("2/3"), None)),
        );
        assert_eq!(
            mes(&out),
            [
                "TTAATTTTCCTTTATTTCAGGAT",
                "TTAATTTTCCTTTATTTCAGTAT",
                "7.343",
                "8.473",
                "-1.130"
            ]
        );
    }

    /// Ensembl VEP 116.2 walks the overlapped-intron list memoised by the
    /// consequence predicates with the variant's own span over introns widened
    /// three bases on each side (`BaseTranscriptVariation.pm:802-822`, `:1027`;
    /// `BaseVariationFeatureOverlapAllele.pm:500`; `MaxEntScan.pm:270`), so a
    /// substitution in a canonical window of an intron shorter than 17 bases
    /// that lies more than three bases outside the intron is left blank, while
    /// one within three bases of it is scored; a variant outside the transcript
    /// meets no memoised list, and the plugin's own 21-base window reaches the
    /// intron.
    #[test]
    fn short_intron_windows_score_only_within_three_bases_of_the_intron() {
        let dir = tempfile::tempdir().unwrap();
        let mes_seq_len =
            |out: &IndexMap<String, String>| out.get("MaxEntScan_ref_seq").map(String::len);
        // Forward strand: a two-base intron 201-202 and a 20-base intron 301-320.
        let mut p = plugin(&["verbose"]);
        let fa = fasta(dir.path(), "23", 400, &[]);
        p.set_reference_fasta(Some(fa));
        let ctx = context(1, &[(101, 200), (203, 300), (321, 400)]);
        // 197 lies in the acceptor 23-mer (183-205) of the two-base intron and
        // four bases before it; 206 in its donor 9-mer (198-206) and four bases
        // after it: neither is scored.
        for pos in [197u64, 206] {
            let out = annotate(
                &p,
                variant("23", pos, pos, b"A", b"C", row(&ctx, Some("1/3"), None)),
            );
            assert!(!out.contains_key("MaxEntScan_ref"), "{pos}: {out:?}");
        }
        // 198 and 205, three bases from the intron, are scored on its donor window.
        for pos in [198u64, 205] {
            let out = annotate(
                &p,
                variant("23", pos, pos, b"A", b"C", row(&ctx, Some("1/3"), None)),
            );
            assert_eq!(mes_seq_len(&out), Some(9), "{pos}: {out:?}");
        }
        // The 20-base intron's acceptor window (301-323) lies inside its widened span.
        let out = annotate(
            &p,
            variant("23", 316, 316, b"A", b"C", row(&ctx, None, Some("2/2"))),
        );
        assert_eq!(mes_seq_len(&out), Some(23), "{out:?}");

        // Reverse strand: the two-base intron's donor 9-mer is 197-205 and its
        // acceptor 23-mer 198-220.
        let ctx = context(-1, &[(101, 200), (203, 300)]);
        for pos in [197u64, 210] {
            let out = annotate(
                &p,
                variant("23", pos, pos, b"A", b"C", row(&ctx, Some("2/2"), None)),
            );
            assert!(!out.contains_key("MaxEntScan_ref"), "{pos}: {out:?}");
        }
        let out = annotate(
            &p,
            variant("23", 204, 204, b"A", b"C", row(&ctx, Some("1/2"), None)),
        );
        assert_eq!(mes_seq_len(&out), Some(9), "{out:?}");

        // Outside the transcript (a one-base first exon at 101, a two-base
        // intron 102-103): the upstream substitution at 96 lies in the acceptor
        // 23-mer (84-106) and is scored.
        let ctx = context(1, &[(101, 101), (104, 200)]);
        let out = annotate(&p, variant("23", 96, 96, b"A", b"C", row(&ctx, None, None)));
        assert_eq!(mes_seq_len(&out), Some(23), "{out:?}");
    }

    #[test]
    fn swa_picks_the_highest_kmer_and_snv_reference_frame() {
        let dir = tempfile::tempdir().unwrap();
        let mut p = plugin(&["SWA", "verbose"]);
        let (fa, ctx) = forward_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let out = annotate(
            &p,
            variant("21", 123, 123, b"G", b"A", row(&ctx, Some("1/3"), None)),
        );
        let field = |k: &str| out.get(k).cloned().unwrap_or_default();
        assert_eq!(field("MES-SWA_donor_ref_context"), "TTCCTTCAGGTACCGTC");
        assert_eq!(field("MES-SWA_donor_ref_frame"), "6");
        assert_eq!(field("MES-SWA_donor_ref_seq"), "TCAGGTACC");
        assert_eq!(field("MES-SWA_donor_ref"), "1.066");
        assert_eq!(field("MES-SWA_donor_alt_context"), "TTCCTTCAAGTACCGTC");
        assert_eq!(field("MES-SWA_donor_alt_frame"), "7");
        assert_eq!(field("MES-SWA_donor_alt_seq"), "CAAGTACCG");
        assert_eq!(field("MES-SWA_donor_alt"), "6.055");
        // A substitution compares the reference in the alternate k-mer's frame.
        assert_eq!(field("MES-SWA_donor_ref_comp_seq"), "CAGGTACCG");
        assert_eq!(field("MES-SWA_donor_ref_comp"), "0.312");
        assert_eq!(field("MES-SWA_donor_diff"), "-5.743");
        assert_eq!(
            field("MES-SWA_acceptor_ref_context"),
            "TTCTGTACGAGCTGTTCCTTCAGGTACCGTCTCTGGGAAGGGGCG"
        );
        assert_eq!(field("MES-SWA_acceptor_ref_frame"), "4");
        assert_eq!(field("MES-SWA_acceptor_ref_seq"), "TGTACGAGCTGTTCCTTCAGGTA");
        assert_eq!(field("MES-SWA_acceptor_ref"), "25.228");
        assert_eq!(field("MES-SWA_acceptor_alt_frame"), "22");
        assert_eq!(field("MES-SWA_acceptor_alt_seq"), "AAGTACCGTCTCTGGGAAGGGGC");
        assert_eq!(field("MES-SWA_acceptor_alt"), "22.130");
        assert_eq!(
            field("MES-SWA_acceptor_ref_comp_seq"),
            "AGGTACCGTCTCTGGGAAGGGGC"
        );
        assert_eq!(field("MES-SWA_acceptor_ref_comp"), "21.683");
        assert_eq!(field("MES-SWA_acceptor_diff"), "-0.447");

        // A deletion compares against the reference's own best k-mer, and its
        // alternate context is the reference context without the deleted bases.
        let out = annotate(
            &p,
            variant("21", 124, 125, b"GT", b"-", row(&ctx, None, Some("1/2"))),
        );
        let field = |k: &str| out.get(k).cloned().unwrap_or_default();
        assert_eq!(field("MES-SWA_donor_ref_context"), "TCCTTCAGGTACCGTCTC");
        assert_eq!(field("MES-SWA_donor_alt_context"), "TCCTTCAGACCGTCTC");
        assert_eq!(
            field("MES-SWA_donor_ref_comp_seq"),
            field("MES-SWA_donor_ref_seq")
        );
        assert_eq!(field("MES-SWA_donor_ref_comp"), field("MES-SWA_donor_ref"));
        let t = tables();
        let (frame, score) = t.best_kmer(b"TCCTTCAGACCGTCTC", 9, Tables::score5);
        let best = &b"TCCTTCAGACCGTCTC"[frame - 1..frame - 1 + 9];
        assert_eq!(
            field("MES-SWA_donor_alt_seq"),
            String::from_utf8_lossy(best)
        );
        assert_eq!(field("MES-SWA_donor_alt_frame"), frame.to_string());
        assert_eq!(field("MES-SWA_donor_alt"), score_text(score));
        assert_eq!(
            field("MES-SWA_donor_diff"),
            score_text(t.score5(field("MES-SWA_donor_ref_seq").as_bytes()) - score)
        );

        // The reverse strand reads both contexts on the transcript's strand.
        let mut p = plugin(&["SWA", "verbose"]);
        let (fa, ctx) = reverse_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let out = annotate(
            &p,
            variant("22", 223, 223, b"C", b"T", row(&ctx, Some("1/3"), None)),
        );
        let field = |k: &str| out.get(k).cloned().unwrap_or_default();
        assert_eq!(field("MES-SWA_donor_ref_context"), "CGATCTTAGGTGGGTTC");
        assert_eq!(field("MES-SWA_donor_alt_context"), "CGATCTTAAGTGGGTTC");
        assert_eq!(field("MES-SWA_donor_ref_frame"), "7");
        assert_eq!(field("MES-SWA_donor_alt_frame"), "7");
        assert_eq!(field("MES-SWA_donor_ref"), "3.450");
        assert_eq!(field("MES-SWA_donor_alt"), "4.541");
        assert_eq!(field("MES-SWA_donor_diff"), "-1.092");
        assert_eq!(
            field("MES-SWA_acceptor_ref_context"),
            "GGCCAGAGAGATGACGATCTTAGGTGGGTTCCCGGGGGTCCTTGT"
        );
        assert_eq!(field("MES-SWA_acceptor_ref_seq"), "GAGATGACGATCTTAGGTGGGTT");
        assert_eq!(field("MES-SWA_acceptor_ref"), "17.338");
        assert_eq!(field("MES-SWA_acceptor_alt_seq"), "GAGATGACGATCTTAAGTGGGTT");
        assert_eq!(field("MES-SWA_acceptor_alt"), "14.307");
        assert_eq!(field("MES-SWA_acceptor_diff"), "3.031");
    }

    #[test]
    fn ncss_from_intron_number_and_from_exon_number() {
        let dir = tempfile::tempdir().unwrap();
        let mut p = plugin(&["NCSS", "verbose"]);
        let (fa, ctx) = forward_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let ncss = |out: &IndexMap<String, String>| -> Vec<(String, String)> {
            out.iter()
                .filter(|(k, _)| k.starts_with("MES-NCSS_") && !k.ends_with("_seq"))
                .map(|(k, v)| (k.clone(), format!("{v} {}", out[&format!("{k}_seq")])))
                .collect()
        };
        let pairs = |items: &[(&str, &str)]| -> Vec<(String, String)> {
            items
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        };
        // Exon 1: no upstream site; the donor of exon 1 and the acceptor of exon 2.
        let out = annotate(
            &p,
            variant("21", 123, 123, b"G", b"A", row(&ctx, Some("1/3"), None)),
        );
        assert_eq!(
            ncss(&out),
            pairs(&[
                ("MES-NCSS_downstream_donor", "0.312 CAGGTACCG"),
                (
                    "MES-NCSS_downstream_acceptor",
                    "-1.563 AATTTTTAAGTTTCTTATAGGCG"
                ),
            ])
        );
        // Intron 1: its own donor and acceptor, the next intron's donor, no upstream acceptor.
        let out = annotate(
            &p,
            variant("21", 125, 125, b"T", b"C", row(&ctx, None, Some("1/2"))),
        );
        assert_eq!(
            ncss(&out),
            pairs(&[
                ("MES-NCSS_upstream_donor", "0.312 CAGGTACCG"),
                (
                    "MES-NCSS_downstream_acceptor",
                    "-1.563 AATTTTTAAGTTTCTTATAGGCG"
                ),
                ("MES-NCSS_downstream_donor", "1.725 TTGGTAGGT"),
            ])
        );
        // Exon 2: every side.
        let out = annotate(
            &p,
            variant("21", 221, 221, b"G", b"A", row(&ctx, Some("2/3"), None)),
        );
        assert_eq!(
            ncss(&out),
            pairs(&[
                ("MES-NCSS_upstream_donor", "0.312 CAGGTACCG"),
                (
                    "MES-NCSS_upstream_acceptor",
                    "-1.563 AATTTTTAAGTTTCTTATAGGCG"
                ),
                ("MES-NCSS_downstream_donor", "1.725 TTGGTAGGT"),
                (
                    "MES-NCSS_downstream_acceptor",
                    "10.637 TTACTCTGTTTCTTTCAAAGGGA"
                ),
            ])
        );
        // The last intron has no downstream donor; the last exon no downstream site.
        let out = annotate(
            &p,
            variant("21", 275, 275, b"A", b"C", row(&ctx, None, Some("2/2"))),
        );
        assert_eq!(
            ncss(&out),
            pairs(&[
                ("MES-NCSS_upstream_donor", "1.725 TTGGTAGGT"),
                (
                    "MES-NCSS_downstream_acceptor",
                    "10.637 TTACTCTGTTTCTTTCAAAGGGA"
                ),
                (
                    "MES-NCSS_upstream_acceptor",
                    "-1.563 AATTTTTAAGTTTCTTATAGGCG"
                ),
            ])
        );
        let out = annotate(
            &p,
            variant("21", 295, 295, b"A", b"C", row(&ctx, Some("3/3"), None)),
        );
        assert_eq!(
            ncss(&out),
            pairs(&[
                ("MES-NCSS_upstream_donor", "1.725 TTGGTAGGT"),
                (
                    "MES-NCSS_upstream_acceptor",
                    "10.637 TTACTCTGTTTCTTTCAAAGGGA"
                ),
            ])
        );
        // A row with neither number (an upstream variant) gets no site.
        let out = annotate(&p, variant("21", 30, 30, b"A", b"C", row(&ctx, None, None)));
        assert!(ncss(&out).is_empty());

        // The reverse strand: sites at the introns' other ends, read on the transcript's strand.
        let mut p = plugin(&["NCSS", "verbose"]);
        let (fa, ctx) = reverse_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let out = annotate(
            &p,
            variant("22", 221, 221, b"A", b"G", row(&ctx, None, Some("1/2"))),
        );
        assert_eq!(
            ncss(&out),
            pairs(&[
                ("MES-NCSS_upstream_donor", "3.450 TAGGTGGGT"),
                (
                    "MES-NCSS_downstream_acceptor",
                    "7.343 TTAATTTTCCTTTATTTCAGGAT"
                ),
                ("MES-NCSS_downstream_donor", "8.192 AAGGTAAGT"),
            ])
        );
        let out = annotate(
            &p,
            variant("22", 139, 139, b"C", b"A", row(&ctx, Some("2/3"), None)),
        );
        assert_eq!(
            ncss(&out),
            pairs(&[
                ("MES-NCSS_upstream_donor", "3.450 TAGGTGGGT"),
                (
                    "MES-NCSS_upstream_acceptor",
                    "7.343 TTAATTTTCCTTTATTTCAGGAT"
                ),
                ("MES-NCSS_downstream_donor", "8.192 AAGGTAAGT"),
                (
                    "MES-NCSS_downstream_acceptor",
                    "-7.316 TTTGCAAATCCTATTCACAGGTT"
                ),
            ])
        );
    }

    #[test]
    fn verbose_adds_sequences_and_frames() {
        let dir = tempfile::tempdir().unwrap();
        let (fa, ctx) = forward_layout(dir.path());
        let mut plain = plugin(&["SWA", "NCSS"]);
        plain.set_reference_fasta(Some(fa.clone()));
        let mut verbose = plugin(&["SWA", "NCSS", "verbose"]);
        verbose.set_reference_fasta(Some(fa));
        assert_eq!(plain.header_info().len(), 3 + 8 + 4);
        assert_eq!(verbose.header_info().len(), 3 + 2 + 8 + 14 + 4 + 4);
        assert_eq!(plugin(&[]).header_info().len(), 3);
        assert_eq!(plugin(&["verbose"]).header_info().len(), 5);

        let make = || variant("21", 123, 123, b"G", b"A", row(&ctx, Some("1/3"), None));
        let plain_out = annotate(&plain, make());
        let verbose_out = annotate(&verbose, make());
        let mut plain_keys: Vec<&String> = plain_out.keys().collect();
        plain_keys.sort();
        assert_eq!(
            plain_keys,
            [
                "MES-NCSS_downstream_acceptor",
                "MES-NCSS_downstream_donor",
                "MES-SWA_acceptor_alt",
                "MES-SWA_acceptor_diff",
                "MES-SWA_acceptor_ref",
                "MES-SWA_acceptor_ref_comp",
                "MES-SWA_donor_alt",
                "MES-SWA_donor_diff",
                "MES-SWA_donor_ref",
                "MES-SWA_donor_ref_comp",
                "MaxEntScan_alt",
                "MaxEntScan_diff",
                "MaxEntScan_ref",
            ]
        );
        for (key, value) in &plain_out {
            assert_eq!(verbose_out.get(key), Some(value), "{key}");
        }
        let mut extra: Vec<&String> = verbose_out
            .keys()
            .filter(|k| !plain_out.contains_key(*k))
            .collect();
        extra.sort();
        assert_eq!(
            extra,
            [
                "MES-NCSS_downstream_acceptor_seq",
                "MES-NCSS_downstream_donor_seq",
                "MES-SWA_acceptor_alt_context",
                "MES-SWA_acceptor_alt_frame",
                "MES-SWA_acceptor_alt_seq",
                "MES-SWA_acceptor_ref_comp_seq",
                "MES-SWA_acceptor_ref_context",
                "MES-SWA_acceptor_ref_frame",
                "MES-SWA_acceptor_ref_seq",
                "MES-SWA_donor_alt_context",
                "MES-SWA_donor_alt_frame",
                "MES-SWA_donor_alt_seq",
                "MES-SWA_donor_ref_comp_seq",
                "MES-SWA_donor_ref_context",
                "MES-SWA_donor_ref_frame",
                "MES-SWA_donor_ref_seq",
                "MaxEntScan_alt_seq",
                "MaxEntScan_ref_seq",
            ]
        );
    }

    #[test]
    fn scores_are_formatted_to_three_decimals() {
        assert_eq!(score_text(0.3124), "0.312");
        assert_eq!(score_text(-13.86049), "-13.860");
        assert_eq!(score_text(2.0), "2.000");
        assert_eq!(score_text(-0.0004), "-0.000");
        let dir = tempfile::tempdir().unwrap();
        let mut p = plugin(&["SWA", "NCSS"]);
        let (fa, ctx) = forward_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let out = annotate(
            &p,
            variant("21", 221, 221, b"G", b"A", row(&ctx, Some("2/3"), None)),
        );
        assert_eq!(out.len(), 3 + 8 + 4);
        for (key, value) in &out {
            let (whole, fraction) = value
                .split_once('.')
                .unwrap_or_else(|| panic!("{key}={value}"));
            assert!(
                whole
                    .trim_start_matches('-')
                    .bytes()
                    .all(|b| b.is_ascii_digit()),
                "{key}={value}"
            );
            assert!(
                fraction.len() == 3 && fraction.bytes().all(|b| b.is_ascii_digit()),
                "{key}={value}"
            );
        }
    }

    #[test]
    fn missing_fasta_is_an_init_error() {
        let dir = fixture_dir().display().to_string();
        let mut registry = crate::BuiltinRegistry::new();
        let outcome =
            registry.load_from_args_with_context(&[format!("MaxEntScan,{dir}")], None, None);
        let [(name, error)] = &outcome.failed[..] else {
            panic!("one failure, got {:?}", outcome.failed);
        };
        assert_eq!(name, "MaxEntScan");
        assert_eq!(
            error.to_string(),
            "plugin initialization failed: MaxEntScan: a reference FASTA is needed; pass --fasta"
        );
        assert!(registry.is_empty());

        let tmp = tempfile::tempdir().unwrap();
        let (fa, _) = forward_layout(tmp.path());
        let mut registry = crate::BuiltinRegistry::new();
        let outcome = registry.load_from_args_with_context(
            &[format!("MaxEntScan,{dir},SWA")],
            Some(fa),
            None,
        );
        assert!(outcome.failed.is_empty(), "{:?}", outcome.failed);
        assert_eq!(registry.plugin_names(), ["MaxEntScan"]);
    }

    #[test]
    fn rows_without_transcript_context_and_non_transcript_rows_are_left_alone() {
        let dir = tempfile::tempdir().unwrap();
        let mut p = plugin(&["SWA"]);
        let (fa, ctx) = forward_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let mut without_context = row(&ctx, Some("1/3"), None);
        without_context.loftee_ctx = None;
        assert!(annotate(&p, variant("21", 123, 123, b"G", b"A", without_context)).is_empty());
        let mut regulatory = row(&ctx, None, None);
        regulatory.feature_type = FeatureType::RegulatoryFeature;
        assert!(annotate(&p, variant("21", 123, 123, b"G", b"A", regulatory)).is_empty());
    }

    #[test]
    fn structural_variants_are_left_alone() {
        let dir = tempfile::tempdir().unwrap();
        let mut p = plugin(&["SWA", "NCSS", "verbose"]);
        let (fa, ctx) = forward_layout(dir.path());
        p.set_reference_fasta(Some(fa));
        let deleted = b"GCTGTTCCTTCAGGTACCGT";
        let mut symbolic = variant(
            "21",
            111,
            130,
            deleted,
            b"<DEL>",
            row(&ctx, Some("1/3"), None),
        );
        symbolic.is_structural = true;
        assert!(annotate(&p, symbolic).is_empty());
        // The same span given as a sequence deletion is scored.
        let sequence = variant("21", 111, 130, deleted, b"-", row(&ctx, Some("1/3"), None));
        assert!(!annotate(&p, sequence).is_empty());
    }
}
