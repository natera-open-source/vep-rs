// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Coding sequences cut off at an edge of their sequence.
//!
//! Each test lays a made-up gene model on a contig of its own, with the
//! transcript's first base or its last on an edge of the contig, and gives the
//! model as a cache carries it: the first exon's phase, the coding sequence
//! padded with one `N` per base its first codon lacks, and the `cds_start_NF`
//! and `cds_end_NF` flags or neither. The expected terms and columns are the
//! outcome of Ensembl VEP's predicates, cited at each assertion as
//! `Module.pm:line`, but for the terms of a start codon and of a stop codon.
//!
//! Ensembl VEP reads the first three coding positions as the start codon and
//! the last three as the stop codon on every transcript without the flag of
//! that end, whatever they hold (VariationEffect.pm:975-986, 1373-1381).
//! vep-rs reads a start codon and a stop codon only where one stands. A model
//! here has none at the end its contig's edge cuts, so with the flag of that
//! end or without it a change on those positions loses no start and no stop.

use crate::consequences::{calculate_consequences, EffectsConfig};
use crate::test_helpers::make_single_exon_transcript;
use std::sync::Arc;
use vep_core::codon::reverse_complement;
use vep_core::consequence::TranscriptConsequence;
use vep_core::coordinate::Strand;
use vep_core::transcript::{Exon, ExonCoordMapper, Intron, MapperPair, Transcript};
use vep_core::variant::InputVariant;

/// The contig every model lies on.
const CONTIG: &str = "scaffold_12";

/// The length of the intron between two exons of a model.
const INTRON: u64 = 60;

/// The 3' UTR of a model whose coding sequence a stop codon closes.
const UTR3: &str = "CTGATTCCAGGAGTCTCCACCAGGA";

/// The terms of a change inside a partial last codon that leaves both alleles
/// an unknown residue.
const IN_A_PARTIAL_LAST_CODON: [&str; 2] = [
    "incomplete_terminal_codon_variant",
    "coding_sequence_variant",
];

/// A gene model and the room its contig leaves around it.
struct Shape<'a> {
    strand: Strand,
    /// The exons in transcript order, each as the bases the transcript reads.
    exons: &'a [&'a str],
    /// The contig's bases before the transcript's first base and after its
    /// last, in the transcript's direction: `0` puts that end of the transcript
    /// on an edge of the contig.
    before: u64,
    after: u64,
    /// The cDNA positions of the first and the last coding base.
    coding: (u64, u64),
    /// The first exon's phase: the bases its first codon lacks, `0` where the
    /// coding sequence opens on a codon.
    start_phase: i8,
    flags: &'a [&'a str],
}

/// Where the exons of a shape lie on its contig.
struct Layout {
    reverse: bool,
    before: u64,
    exon_lengths: Vec<u64>,
    contig_length: u64,
}

impl Layout {
    /// The genomic position of cDNA position `cdna`.
    fn genomic(&self, cdna: u64) -> u64 {
        let mut along = self.before;
        let mut first = 1;
        for &length in &self.exon_lengths {
            if cdna < first + length {
                along += cdna - first;
                break;
            }
            along += length + INTRON;
            first += length;
        }
        if self.reverse {
            self.contig_length - along
        } else {
            along + 1
        }
    }
}

/// A shape's transcript over a reference that holds its contig alone.
struct Model {
    _dir: tempfile::TempDir,
    fasta: Arc<vep_fasta::IndexedFasta>,
    transcript: Arc<Transcript>,
    layout: Layout,
    /// The transcript's bases, exon after exon.
    cdna: Vec<u8>,
}

impl Model {
    fn new(shape: &Shape<'_>) -> Self {
        let reverse = shape.strand == Strand::Reverse;
        let mut laid = vec![b'T'; shape.before as usize];
        let mut cdna = Vec::new();
        for (i, exon) in shape.exons.iter().enumerate() {
            if i > 0 {
                laid.extend(std::iter::repeat_n(b'T', INTRON as usize));
            }
            laid.extend_from_slice(exon.as_bytes());
            cdna.extend_from_slice(exon.as_bytes());
        }
        laid.extend(std::iter::repeat_n(b'T', shape.after as usize));
        let layout = Layout {
            reverse,
            before: shape.before,
            exon_lengths: shape.exons.iter().map(|exon| exon.len() as u64).collect(),
            contig_length: laid.len() as u64,
        };
        let contig = if reverse {
            reverse_complement(&laid)
        } else {
            laid
        };

        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("ref.fa");
        let header = format!(">{CONTIG}\n");
        let mut file = header.clone().into_bytes();
        file.extend_from_slice(&contig);
        file.push(b'\n');
        std::fs::write(&path, file).expect("write FASTA");
        let length = layout.contig_length;
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            format!(
                "{CONTIG}\t{length}\t{}\t{length}\t{}\n",
                header.len(),
                length + 1
            ),
        )
        .expect("write FASTA index");
        let fasta = vep_fasta::IndexedFasta::from_path_mmap(&path).expect("load indexed FASTA");

        let transcript = transcript_of(shape, &layout, &cdna);
        Model {
            _dir: dir,
            fasta: Arc::new(fasta),
            transcript: Arc::new(transcript),
            layout,
            cdna,
        }
    }

    /// The row of `variant` on the model's transcript, with the contig as the
    /// reference sequence. Ensembl VEP's coding predicates read the sequences
    /// the transcript carries (BaseTranscriptVariation.pm:1083-1115), so the
    /// row without a reference sequence has the same terms.
    fn row(&self, variant: &InputVariant) -> TranscriptConsequence {
        let with_reference = EffectsConfig {
            reference_fasta: Some(Arc::clone(&self.fasta)),
            ..EffectsConfig::default()
        };
        let row = calculate_consequences(variant, &self.transcript, &with_reference)
            .expect("a row of the transcript");
        let without = calculate_consequences(variant, &self.transcript, &EffectsConfig::default())
            .expect("a row of the transcript");
        assert_eq!(terms_of(&row), terms_of(&without));
        row
    }

    /// `bases`, which the transcript reads, as the forward strand holds them.
    fn on_forward_strand(&self, bases: &[u8]) -> Vec<u8> {
        if self.layout.reverse {
            reverse_complement(bases)
        } else {
            bases.to_vec()
        }
    }

    /// The base at cDNA position `cdna` replaced by `to`, as the transcript
    /// reads it.
    fn substitution(&self, cdna: u64, to: u8) -> InputVariant {
        let at = self.layout.genomic(cdna);
        InputVariant::new(
            CONTIG.into(),
            at,
            at,
            self.on_forward_strand(&[self.cdna[cdna as usize - 1]]),
            self.on_forward_strand(&[to]),
        )
    }

    /// The cDNA positions `from..=to` deleted.
    fn deletion(&self, from: u64, to: u64) -> InputVariant {
        let (a, b) = (self.layout.genomic(from), self.layout.genomic(to));
        InputVariant::new(
            CONTIG.into(),
            a.min(b),
            a.max(b),
            self.on_forward_strand(&self.cdna[from as usize - 1..to as usize]),
            b"-".to_vec(),
        )
    }

    /// `bases`, as the transcript reads them, inserted after cDNA position
    /// `cdna`.
    fn insertion_after(&self, cdna: u64, bases: &[u8]) -> InputVariant {
        let at = self.layout.genomic(cdna);
        let (start, end) = if self.layout.reverse {
            (at, at - 1)
        } else {
            (at + 1, at)
        };
        InputVariant::new(
            CONTIG.into(),
            start,
            end,
            b"-".to_vec(),
            self.on_forward_strand(bases),
        )
    }
}

/// The transcript of `shape` as a cache carries it: the first exon's phase on
/// the exon and on the mapper (TranscriptMapper.pm:108-113), and the coding
/// sequence padded with one `N` per base of that phase (Transcript.pm:917-920).
fn transcript_of(shape: &Shape<'_>, layout: &Layout, cdna: &[u8]) -> Transcript {
    let (coding_start, coding_end) = shape.coding;
    let phase_at = |position: u64| -> i8 {
        if position < coding_start || position > coding_end {
            -1
        } else if position == coding_start {
            shape.start_phase
        } else {
            ((shape.start_phase.max(0) as u64 + position - coding_start) % 3) as i8
        }
    };
    let mut pairs = Vec::new();
    let mut exons = Vec::new();
    let mut first = 1;
    for (i, &length) in layout.exon_lengths.iter().enumerate() {
        let last = first + length - 1;
        let (a, b) = (layout.genomic(first), layout.genomic(last));
        let (start, end) = (a.min(b), a.max(b));
        pairs.push(MapperPair {
            from_start: first,
            from_end: last,
            to_start: start,
            to_end: end,
            ori: if layout.reverse { -1 } else { 1 },
        });
        exons.push(Exon {
            start,
            end,
            rank: i as u32 + 1,
            phase: phase_at(first),
            end_phase: if last < coding_end {
                phase_at(last + 1)
            } else {
                -1
            },
        });
        first = last + 1;
    }
    let introns: Vec<Intron> = exons
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            let (start, end) = if layout.reverse {
                (pair[1].end + 1, pair[0].start - 1)
            } else {
                (pair[0].end + 1, pair[1].start - 1)
            };
            Intron {
                start,
                end,
                rank: i as u32 + 1,
            }
        })
        .collect();
    let exon_of = |position: u64| {
        pairs
            .iter()
            .position(|p| p.from_start <= position && position <= p.from_end)
            .expect("a cDNA position on an exon")
    };

    let text = |bases: &[u8]| String::from_utf8(bases.to_vec()).expect("ASCII bases");
    let mut cds = "N".repeat(shape.start_phase.max(0) as usize);
    cds.push_str(&text(&cdna[coding_start as usize - 1..coding_end as usize]));
    let utr5 = text(&cdna[..coding_start as usize - 1]);
    let utr3 = text(&cdna[coding_end as usize..]);

    let mut tx = make_single_exon_transcript(shape.strand, &cds, &utr3);
    tx.chr = CONTIG.into();
    tx.start = exons.iter().map(|e| e.start).min().expect("an exon");
    tx.end = exons.iter().map(|e| e.end).max().expect("an exon");
    tx.exons = exons.clone();
    tx.introns = introns.clone();
    let (at_start, at_end) = (layout.genomic(coding_start), layout.genomic(coding_end));
    tx.cdna_coding_start = Some(coding_start);
    tx.cdna_coding_end = Some(coding_end);
    tx.coding_region_start = Some(at_start.min(at_end));
    tx.coding_region_end = Some(at_start.max(at_end));
    tx.translation_start = Some(at_start);
    tx.translation_end = Some(at_end);
    tx.flags = shape
        .flags
        .iter()
        .map(|flag| (*flag).to_string())
        .collect::<Vec<_>>()
        .into();
    let (start_exon, end_exon) = (exon_of(coding_start), exon_of(coding_end));
    let translation = tx.translation.as_mut().expect("a translation");
    translation.start = coding_start - pairs[start_exon].from_start + 1;
    translation.end = coding_end - pairs[end_exon].from_start + 1;
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.sorted_exons = exons;
    vefc.introns = introns;
    vefc.five_prime_utr = (!utr5.is_empty()).then_some(utr5);
    vefc.three_prime_utr = (!utr3.is_empty()).then_some(utr3);
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    mapper.start_phase = shape.start_phase;
    mapper.cdna_coding_start = coding_start;
    mapper.cdna_coding_end = coding_end;
    mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
    tx
}

/// The row's terms, sorted.
fn terms_of(row: &TranscriptConsequence) -> Vec<String> {
    let mut terms: Vec<String> = row.consequences.iter().map(|c| c.to_string()).collect();
    terms.sort();
    terms
}

/// Asserts that the row of `case` carries exactly `expected`, in any order.
fn assert_terms(row: &TranscriptConsequence, expected: &[&str], case: &str) {
    let mut expected: Vec<String> = expected.iter().map(|term| (*term).to_string()).collect();
    expected.sort();
    assert_eq!(terms_of(row), expected, "{case}");
}

/// Asserts the row's `CDS_position`, `Protein_position`, `Amino_acids` and
/// `Codons` columns, in that order.
fn assert_columns(row: &TranscriptConsequence, expected: [&str; 4], case: &str) {
    let columns = [
        row.cds_position.as_deref(),
        row.protein_position.as_deref(),
        row.amino_acids.as_deref(),
        row.codons.as_deref(),
    ];
    assert_eq!(columns, expected.map(Some), "{case}");
}

/// The two flag sets a shape is given with: none, and `flag`.
fn with_and_without(flag: &'static str) -> [Vec<&'static str>; 2] {
    [Vec::new(), vec![flag]]
}

/// A model of two exons whose first base lies on the edge of its contig and
/// whose coding sequence opens there, `start_phase` bases into a codon: the
/// rest of that codon, eleven whole codons (Ala, Arg, Asp, Asn, Cys, Glu, His,
/// Ile, Lys, Val, Leu), a stop codon and a 3' UTR. With phase `0` the coding
/// sequence opens on the Ala codon, which is no start codon.
fn begins_at_the_edge(strand: Strand, start_phase: i8, flags: &[&str]) -> Model {
    let opening = match start_phase {
        2 => "C",
        1 => "CC",
        _ => "",
    };
    let exon_1 = format!("{opening}GCACGTGACAACTGCGAG");
    let exon_2 = format!("CATATCAAGGTCCTGTAA{UTR3}");
    Model::new(&Shape {
        strand,
        exons: &[&exon_1, &exon_2],
        before: 0,
        after: 50,
        coding: (1, opening.len() as u64 + 36),
        start_phase,
        flags,
    })
}

/// A model of two exons whose last base lies on the edge of its contig and
/// ends its coding sequence, which no stop codon closes: a 12-base 5' UTR, a
/// start codon, the eleven whole codons of [`begins_at_the_edge`] and `tail`,
/// the bases of a codon the edge cuts. Its last whole codon is cDNA 46 to 48.
fn ends_at_the_edge(strand: Strand, tail: &str, flags: &[&str]) -> Model {
    let exon_2 = format!("GAGCATATCAAGGTCCTG{tail}");
    Model::new(&Shape {
        strand,
        exons: &["GAGTCTCCAGGAATGGCACGTGACAACTGC", &exon_2],
        before: 40,
        after: 0,
        coding: (13, 48 + tail.len() as u64),
        start_phase: -1,
        flags,
    })
}

/// The first base of a transcript that begins mid-codon is the first base of
/// its contig on the forward strand and the last on the reverse strand. Its
/// first whole codon is residue 2 at coding bases 4 to 6, with the flag and
/// without.
#[test]
fn a_coding_sequence_that_begins_mid_codon_at_base_one_is_read_in_its_phase() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for start_phase in [1, 2] {
            for flags in with_and_without("cds_start_NF") {
                let case = format!("{strand:?}, phase {start_phase}, flags {flags:?}");
                let model = begins_at_the_edge(strand, start_phase, &flags);
                let edge = match strand {
                    Strand::Forward => 1,
                    Strand::Reverse => model.layout.contig_length,
                };
                assert_eq!(model.layout.genomic(1), edge, "{case}");
                // The cDNA position of the partial codon's last base.
                let partial = 3 - start_phase as u64;

                // `cds_start` and `genomic2pep` count the bases the first
                // codon lacks (BaseTranscriptVariation.pm:261-264,
                // TranscriptMapper.pm:515-529); `missense_variant`
                // (VariationEffect.pm:1085-1098).
                let row = model.row(&model.substitution(partial + 2, b'T'));
                assert_terms(&row, &["missense_variant"], &case);
                assert_columns(&row, ["5", "2", "A/V", "gCa/gTa"], &case);

                // Both alleles of the partial codon read the unknown residue:
                // `coding_unknown` (VariationEffect.pm:1629-1630), no
                // `start_lost` (VariationEffect.pm:875, and :975 with the
                // flag), and no `partial_codon`, which holds for the last
                // codon alone (VariationEffect.pm:1606-1610).
                let row = model.row(&model.substitution(partial, b'T'));
                assert_terms(&row, &["coding_sequence_variant"], &case);
                let codons = if start_phase == 2 {
                    "nnC/nnT"
                } else {
                    "ncC/ncT"
                };
                assert_columns(&row, ["3", "1", "X", codons], &case);
            }
        }
    }
}

/// With `cds_start_NF`, a substitution, a deletion and an insertion on the
/// first three coding bases are no lost start, whether the coding sequence
/// opens inside a codon or on a whole codon that is no start codon.
#[test]
fn a_coding_sequence_without_a_start_codon_has_no_start_lost_with_the_flag() {
    // `_overlaps_start_codon` returns 0 for the flag (VariationEffect.pm:975),
    // so `start_lost` returns 0 at its first test (VariationEffect.pm:862).
    for strand in [Strand::Forward, Strand::Reverse] {
        for start_phase in [2, 1] {
            let case = format!("{strand:?}, phase {start_phase}");
            let model = begins_at_the_edge(strand, start_phase, &["cds_start_NF"]);
            let partial = 3 - start_phase as u64;
            let row = model.row(&model.substitution(partial, b'T'));
            assert_terms(&row, &["coding_sequence_variant"], &case);
            let row = model.row(&model.deletion(partial, partial));
            assert_terms(&row, &["frameshift_variant"], &case);
            let row = model.row(&model.deletion(partial + 1, partial + 1));
            assert_terms(&row, &["frameshift_variant"], &case);
            let row = model.row(&model.insertion_after(partial, b"A"));
            assert_terms(&row, &["frameshift_variant"], &case);
        }

        let case = format!("{strand:?}, a whole first codon");
        let model = begins_at_the_edge(strand, 0, &["cds_start_NF"]);
        let row = model.row(&model.substitution(2, b'T'));
        assert_terms(&row, &["missense_variant"], &case);
        let row = model.row(&model.deletion(2, 2));
        assert_terms(&row, &["frameshift_variant"], &case);
        let row = model.row(&model.insertion_after(1, b"A"));
        assert_terms(&row, &["frameshift_variant"], &case);
    }
}

/// Without the flag Ensembl VEP reads the first three coding bases as the
/// start codon (VariationEffect.pm:977-986); these changes on them are no lost
/// start all the same.
#[test]
fn a_coding_sequence_without_a_start_codon_has_no_start_lost_without_the_flag_where_ensembl_vep_writes_none(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        for start_phase in [2, 1] {
            let case = format!("{strand:?}, phase {start_phase}");
            let model = begins_at_the_edge(strand, start_phase, &[]);
            let partial = 3 - start_phase as u64;
            // `start_lost` returns 0 for an alternate peptide that is the
            // unknown residue (VariationEffect.pm:875).
            let row = model.row(&model.substitution(partial, b'T'));
            assert_terms(&row, &["coding_sequence_variant"], &case);
            // The first base of the first whole codon is cDNA position 2 or 3
            // and lies at residue 2, where `start_lost` asks for residue 1
            // (VariationEffect.pm:879).
            let row = model.row(&model.substitution(partial + 1, b'A'));
            assert_terms(&row, &["missense_variant"], &case);
            assert_columns(&row, ["4", "2", "A/T", "Gca/Aca"], &case);
        }

        let case = format!("{strand:?}, a whole first codon");
        let model = begins_at_the_edge(strand, 0, &[]);
        // The alternate peptide is the reference peptide
        // (VariationEffect.pm:880-881).
        let row = model.row(&model.substitution(3, b'G'));
        assert_terms(&row, &["synonymous_variant"], &case);
        assert_columns(&row, ["3", "1", "A", "gcA/gcG"], &case);
        // The deletion of the whole first codon is `inframe_deletion`
        // (VariationEffect.pm:870) and leaves no alternate residue
        // (VariationEffect.pm:875).
        let row = model.row(&model.deletion(1, 3));
        assert_terms(&row, &["inframe_deletion"], &case);
    }
}

/// Without the flag, a substitution that changes the first residue of a
/// coding sequence that opens on a whole codon that is no start codon is the
/// missense variant of
/// `a_coding_sequence_without_a_start_codon_has_no_start_lost_with_the_flag`:
/// no start codon stands there. Ensembl VEP writes `start_lost` for it.
#[test]
fn a_changed_first_residue_of_a_coding_sequence_without_a_start_codon_or_the_flag_is_a_missense_variant(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}");
        let model = begins_at_the_edge(strand, 0, &[]);
        // In Ensembl VEP the allele overlaps the first three coding bases
        // (VariationEffect.pm:975-986), its span begins at residue 1 and the
        // peptide changes: `start_lost` (VariationEffect.pm:878-882).
        let row = model.row(&model.substitution(2, b'T'));
        assert_terms(&row, &["missense_variant"], &case);
        assert_columns(&row, ["2", "1", "A/V", "gCa/gTa"], &case);
    }
}

/// The last base of a transcript whose coding sequence ends inside a codon is
/// the last base of its contig on the forward strand and the first on the
/// reverse strand. A substitution in that codon is
/// `incomplete_terminal_codon_variant`, with the flag and without.
#[test]
fn a_coding_sequence_that_ends_in_a_partial_codon_at_the_last_base_gets_incomplete_terminal_codon_variant(
) {
    // `partial_codon` holds for the last codon where it has one or two bases
    // (VariationEffect.pm:1588-1611) and reads no flag. Both alleles read the
    // unknown residue (TranscriptVariationAllele.pm:774-776), so
    // `coding_unknown` holds as well (VariationEffect.pm:1629), and the
    // missense and stop tests return 0 on a partial codon
    // (VariationEffect.pm:1093, 1240, 1296).
    for strand in [Strand::Forward, Strand::Reverse] {
        for (tail, columns) in [
            ("G", ["37", "13", "X", "G/C"]),
            ("GA", ["38", "13", "X", "gA/gC"]),
        ] {
            for flags in with_and_without("cds_end_NF") {
                let case = format!("{strand:?}, tail {}, flags {flags:?}", tail.len());
                let model = ends_at_the_edge(strand, tail, &flags);
                let last = 48 + tail.len() as u64;
                let edge = match strand {
                    Strand::Forward => model.layout.contig_length,
                    Strand::Reverse => 1,
                };
                assert_eq!(model.layout.genomic(last), edge, "{case}");
                let row = model.row(&model.substitution(last, b'C'));
                assert_terms(&row, &IN_A_PARTIAL_LAST_CODON, &case);
                assert_columns(&row, columns, &case);
            }
        }

        // One exon that is the whole contig, cut at both ends: the first
        // exon's phase is 2 and two bases follow the last whole codon. The
        // coding sequence's length counts its two `N`
        // (VariationEffect.pm:1603): 20 bases, the last two at residue 7.
        for flags in [vec![], vec!["cds_start_NF", "cds_end_NF"]] {
            let case = format!("{strand:?}, cut at both ends, flags {flags:?}");
            let model = Model::new(&Shape {
                strand,
                exons: &["CGCACGTGACAACTGCGA"],
                before: 0,
                after: 0,
                coding: (1, 18),
                start_phase: 2,
                flags: &flags,
            });
            assert_eq!(model.layout.contig_length, 18, "{case}");
            let row = model.row(&model.substitution(18, b'C'));
            assert_terms(&row, &IN_A_PARTIAL_LAST_CODON, &case);
            assert_columns(&row, ["20", "7", "X", "gA/gC"], &case);
        }
    }
}

/// With `cds_end_NF`, a substitution in the last whole codon and a deletion on
/// the last three coding bases are no lost stop, whether the coding sequence
/// ends on a whole codon that is no stop codon or inside a codon.
#[test]
fn a_coding_sequence_without_a_stop_codon_has_no_stop_lost_with_the_flag() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (tail, deletion_of_the_last_base) in [
            ("", &["frameshift_variant"][..]),
            (
                "G",
                &[
                    "incomplete_terminal_codon_variant",
                    "protein_altering_variant",
                ][..],
            ),
            ("GA", &IN_A_PARTIAL_LAST_CODON[..]),
        ] {
            let case = format!("{strand:?}, tail {}", tail.len());
            let model = ends_at_the_edge(strand, tail, &["cds_end_NF"]);
            let last = 48 + tail.len() as u64;
            // The reference peptide holds no stop to lose
            // (VariationEffect.pm:1259-1260).
            let row = model.row(&model.substitution(47, b'C'));
            assert_terms(&row, &["missense_variant"], &case);
            // The alternate peptide holds an unknown residue, and the stop
            // codon edit test that decides then (VariationEffect.pm:1263)
            // returns 0 for the flag (VariationEffect.pm:1373, 1444).
            let row = model.row(&model.deletion(48, 48));
            assert_terms(&row, &["frameshift_variant"], &case);
            // Where a partial codon ends the coding sequence, `stop_lost`
            // returns 0 on it (VariationEffect.pm:1240). The deletion of its
            // only base leaves no alternate residue and is
            // `protein_altering_variant` (VariationEffect.pm:375-393).
            let row = model.row(&model.deletion(last, last));
            assert_terms(&row, deletion_of_the_last_base, &case);
        }
    }
}

/// Without the flag Ensembl VEP reads the last three coding bases as the stop
/// codon (VariationEffect.pm:1375-1381); these changes on them are no lost
/// stop all the same.
#[test]
fn a_coding_sequence_without_a_stop_codon_has_no_stop_lost_without_the_flag_where_ensembl_vep_writes_none(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}, a whole last codon");
        let model = ends_at_the_edge(strand, "", &[]);
        // Both peptides are known and the reference peptide holds no stop to
        // lose (VariationEffect.pm:1259-1260).
        let row = model.row(&model.substitution(47, b'C'));
        assert_terms(&row, &["missense_variant"], &case);
        assert_columns(&row, ["35", "12", "L/P", "cTg/cCg"], &case);
        let row = model.row(&model.deletion(46, 48));
        assert_terms(&row, &["inframe_deletion"], &case);

        let case = format!("{strand:?}, a partial last codon");
        let model = ends_at_the_edge(strand, "GA", &[]);
        // `stop_lost` returns 0 on a partial codon (VariationEffect.pm:1240).
        let row = model.row(&model.substitution(50, b'C'));
        assert_terms(&row, &IN_A_PARTIAL_LAST_CODON, &case);
        let row = model.row(&model.deletion(50, 50));
        assert_terms(&row, &IN_A_PARTIAL_LAST_CODON, &case);
    }
}

/// A one-base insertion or deletion at a cut-off end is a frameshift alone,
/// with the flag of that end and without it: no stop codon stands at the
/// coding end and no start codon at the coding start. Without the flag
/// Ensembl VEP writes `stop_lost` or `start_lost` beside the frameshift where
/// the allele lies on the three coding bases it reads as the stop or the
/// start codon.
#[test]
fn a_frameshift_at_a_cut_off_end_is_a_frameshift_alone_with_the_flag_and_without_it() {
    // `frameshift` compares the allele's length with its coding span
    // (VariationEffect.pm:1556-1564).
    let frameshift = &["frameshift_variant"][..];
    for strand in [Strand::Forward, Strand::Reverse] {
        for flags in with_and_without("cds_end_NF") {
            // The alternate peptide holds an unknown residue, so Ensembl
            // VEP's `stop_lost` is the stop codon edit test
            // (VariationEffect.pm:1263), which is not made for the flag
            // (VariationEffect.pm:1373). Without the flag it reads a lost
            // stop where a deletion leaves the edited sequence shorter than
            // the coding sequence (VariationEffect.pm:1469) and where an
            // insertion leaves three bases that are no stop codon at the
            // coding end (VariationEffect.pm:1473-1480).
            let case = format!("{strand:?}, a whole last codon, flags {flags:?}");
            let model = ends_at_the_edge(strand, "", &flags);
            let row = model.row(&model.deletion(47, 47));
            assert_terms(&row, frameshift, &case);
            let row = model.row(&model.insertion_after(47, b"A"));
            assert_terms(&row, frameshift, &case);
            // One codon before, the allele misses Ensembl VEP's overlap
            // (VariationEffect.pm:1375-1381).
            let row = model.row(&model.deletion(44, 44));
            assert_terms(&row, frameshift, &case);

            // The last three coding bases are cDNA 48 to 50: the last base of
            // the last whole codon and the two bases of the partial codon.
            let case = format!("{strand:?}, a partial last codon, flags {flags:?}");
            let model = ends_at_the_edge(strand, "GA", &flags);
            let row = model.row(&model.deletion(48, 48));
            assert_terms(&row, frameshift, &case);
            let row = model.row(&model.deletion(47, 47));
            assert_terms(&row, frameshift, &case);
            // `frameshift` returns 0 on a partial codon
            // (VariationEffect.pm:1548).
            let row = model.row(&model.deletion(49, 49));
            assert_terms(&row, &IN_A_PARTIAL_LAST_CODON, &case);
        }

        for flags in with_and_without("cds_start_NF") {
            // No allele overlaps a start codon under the flag
            // (VariationEffect.pm:975). Without the flag Ensembl VEP reads a
            // lost start where a deletion on the first three coding bases
            // leaves the edited sequence shorter than the coding sequence
            // (VariationEffect.pm:870, 1069).
            //
            // The first exon's phase is 2: cDNA 1 is the partial codon's base
            // and cDNA 2 to 4 the first whole codon.
            let case = format!("{strand:?}, a partial first codon, flags {flags:?}");
            let model = begins_at_the_edge(strand, 2, &flags);
            let row = model.row(&model.deletion(1, 1));
            assert_terms(&row, frameshift, &case);
            let row = model.row(&model.deletion(2, 2));
            assert_terms(&row, frameshift, &case);
            let row = model.row(&model.deletion(6, 6));
            assert_terms(&row, frameshift, &case);
        }
    }
}
