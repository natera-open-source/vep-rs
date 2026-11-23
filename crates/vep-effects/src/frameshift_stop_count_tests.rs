// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count to the new stop in the protein notation of a frameshift, on a
//! transcript whose peptide goes on past a stop codon before the change: it
//! reads the codon as a residue (a selenocysteine codon) or holds it as a
//! stop and names the residues after it (a polymorphic pseudogene).
//!
//! Each test lays a made-up coding sequence of 24 codons, 23 where no stop
//! codon closes it, on the first exon of the compact test transcript and
//! holds the row's HGVSp, and its HGVSc in some tests, to values derived
//! from those bases. Codon 10 is a stop codon that the peptide reads as a
//! selenocysteine or holds as a stop. Ensembl VEP counts to the first stop
//! of the alternate sequence's translation, which is that codon, or the
//! stop codon that some tests lay at residue 5, before the change, and
//! writes the unknown count on every frameshift after it
//! (TranscriptVariationAllele.pm:2492, 2503, 2508 and 2519).

use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture, make_compact_test_transcript};
use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use vep_core::transcript::{SeqEdit, Transcript};
use vep_core::variant::InputVariant;

/// The cDNA position of the first coding base; exon 1 lies at genomic 100 to
/// 399, so coding base `b` is genomic `149 + b`.
const FIRST_CODING_CDNA: u64 = 51;

/// The residue number of the stop codon that the peptide reads as a residue.
const SELENOCYSTEINE: usize = 10;

/// Codons 21 to 23 of the coding sequence: threonine, leucine and serine. Read
/// one base on, their bases and the next spell a proline codon and a stop
/// codon that is another codon than codon 10.
const TAIL: [&str; 3] = ["ACC", "CTA", "AGC"];

/// The stop codon that closes the coding sequence.
const CLOSING_STOP: &str = "TAG";

/// The 3' UTR: read two bases on from the closing stop's last base, its third
/// to fifth bases are the same stop codon as codon 10.
const UTR: [&str; 3] = ["CA", "TGA", "TTTT"];

/// The coding sequence without its closing stop: a start codon, eight alanine
/// codons, the stop codon read as a selenocysteine, ten alanine codons and
/// the three codons of `TAIL`: 23 codons, 69 bases.
fn open_coding_sequence() -> String {
    format!(
        "ATG{}TGA{}{}",
        "GCT".repeat(8),
        "GCT".repeat(10),
        TAIL.concat()
    )
}

/// The peptide of the coding sequence: its 23 residues, the tenth a
/// selenocysteine.
fn peptide() -> String {
    format!("M{}U{}TLS", "A".repeat(8), "A".repeat(10))
}

/// The compact test transcript with `cds` laid from its first coding base,
/// `utr` after it, the peptide of [`peptide`] and the sequence edit that
/// makes residue 10 a selenocysteine.
fn selenoprotein(cds: &str, utr: Option<&str>) -> Transcript {
    let mut tx = make_compact_test_transcript();
    let last = FIRST_CODING_CDNA + cds.len() as u64 - 1;
    let genomic_last = 99 + last;
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.translateable_seq = Some(cds.to_string());
    vefc.three_prime_utr = utr.map(str::to_string);
    vefc.peptide = Some(peptide());
    vefc.seq_edits = vec![SeqEdit {
        start: SELENOCYSTEINE as u64,
        end: SELENOCYSTEINE as u64,
        alt_seq: "U".to_string(),
    }];
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    mapper.cdna_coding_end = last;
    tx.cdna_coding_end = Some(last);
    tx.coding_region_end = Some(genomic_last);
    tx.translation_end = Some(genomic_last);
    tx
}

/// The transcript whose coding sequence its stop codon closes, with the 3'
/// UTR after it.
fn closed() -> Transcript {
    selenoprotein(
        &format!("{}{CLOSING_STOP}", open_coding_sequence()),
        Some(&UTR.concat()),
    )
}

/// The row of a deletion of the coding bases `first` to `last` of `tx`, on a
/// reference that holds the transcript's coding sequence and 3' UTR.
fn deletion(tx: &Transcript, first: usize, last: usize) -> HgvsNotation {
    let vefc = tx.vefc.as_ref().expect("a coding model");
    let on_genome = format!(
        "{}{}",
        vefc.translateable_seq
            .as_deref()
            .expect("a coding sequence"),
        vefc.three_prime_utr.as_deref().unwrap_or_default()
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(150, &on_genome));
    let deleted = on_genome.as_bytes()[first - 1..last].to_vec();
    let variant = InputVariant::new(
        "1".into(),
        149 + first as u64,
        149 + last as u64,
        deleted,
        b"-".to_vec(),
    );
    generate_hgvs(&variant, tx, Some(&fasta), HgvsOptions::default())
}

/// Coding base 43 is the first base of codon 15, after the selenocysteine
/// codon. Without it the frame reads five leucine codons from the alanine
/// codons 15 to 19, a leucine codon over the end of codon 20, a proline codon
/// and, at residue 22, a stop codon that is another codon than codon 10:
/// eight residues counted from the first changed one. Ensembl VEP writes the
/// unknown count, since the first stop of its alternate translation is
/// residue 10 and 10 - 14 is not positive.
#[test]
fn a_frameshift_past_a_selenocysteine_codon_counts_to_a_stop_of_another_codon() {
    let row = deletion(&closed(), 43, 43);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.43del"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15LeufsTer8")
    );
}

/// Without coding bases 43 and 44 the frame reads five cysteine codons, a
/// tyrosine, a proline, a lysine and a leucine codon, an alanine codon over
/// the last base of the closing stop and, at residue 25, the third to fifth
/// bases of the 3' UTR: the same stop codon as the one the peptide reads as a
/// selenocysteine. Whether that codon ends the protein there is not known
/// from the model, and the count stays unknown.
#[test]
fn the_count_stays_unknown_where_the_new_stop_is_the_selenocysteine_codon() {
    let row = deletion(&closed(), 43, 44);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.43_44del"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15CysfsTer?")
    );
}

/// Coding base 13 is the first base of codon 5, before the selenocysteine
/// codon, whose bases the new frame reads across two codons. The first stop
/// of the alternate translation is the stop at residue 22, 18 residues
/// counted from the first changed one, and no stop codon read as a residue
/// stands before the change.
#[test]
fn a_frameshift_before_the_selenocysteine_codon_keeps_its_count() {
    let row = deletion(&closed(), 13, 13);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.13del"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala5LeufsTer18")
    );
}

/// A coding sequence whose end is not found holds no closing stop and no 3'
/// UTR. The stop the new frame meets at residue 22 lies inside the sequence
/// the model holds, so the count is the same eight.
#[test]
fn a_coding_sequence_whose_end_is_not_found_still_counts_to_a_stop_inside_it() {
    let mut tx = selenoprotein(&open_coding_sequence(), None);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let row = deletion(&tx, 43, 43);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15LeufsTer8")
    );
}

/// The count is read past an earlier stop codon only where the peptide is
/// the coding sequence's translation at every other codon. Here the peptide
/// names a glycine at residue 3, an alanine codon, so it is not the reference
/// protein of these bases and the count stays unknown.
#[test]
fn the_count_stays_unknown_where_the_peptide_is_not_the_translation_elsewhere() {
    let mut tx = closed();
    let vefc = tx.vefc.as_mut().expect("a coding model");
    let mut residues = peptide().into_bytes();
    residues[2] = b'G';
    vefc.peptide = Some(String::from_utf8(residues).expect("ASCII residues"));
    let row = deletion(&tx, 43, 43);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15LeufsTer?")
    );
}

/// The peptide's first residue is a methionine over a first codon that is
/// no methionine codon: a translation writes one over another start codon of
/// the table, and an `initial_met` edit writes one over any other codon. The
/// peptide is still the protein of these bases either way, and the count
/// past the selenocysteine codon is the same eight.
#[test]
fn a_first_codon_that_is_no_methionine_codon_leaves_the_count() {
    let cds = format!("{}{CLOSING_STOP}", open_coding_sequence());
    for (first_codon, edited) in [("CTG", false), ("GCT", true)] {
        let mut tx = selenoprotein(&format!("{first_codon}{}", &cds[3..]), Some(&UTR.concat()));
        if edited {
            tx.vefc
                .as_mut()
                .expect("a coding model")
                .seq_edits
                .push(SeqEdit {
                    start: 1,
                    end: 1,
                    alt_seq: "M".to_string(),
                });
        }
        let row = deletion(&tx, 43, 43);
        assert_eq!(
            row.hgvsp.as_deref(),
            Some("ENSP00000000001.1:p.Ala15LeufsTer8"),
            "{first_codon}"
        );
    }
}

/// With a second stop codon read as a selenocysteine at residue 17, coding
/// base 49 is that codon's first base. Without it the first changed residue
/// is the selenocysteine itself, read as a glutamate from the codon's last
/// two bases and the next; the new frame then reads three leucine codons and
/// a proline codon and meets, at residue 22, a stop codon that is neither
/// selenocysteine codon: six residues counted from the first changed one.
/// The first stop of Ensembl VEP's alternate translation is residue 10.
#[test]
fn a_frameshift_inside_a_later_selenocysteine_codon_counts_past_the_first() {
    let cds = format!(
        "ATG{}TGA{}TGA{}{}{CLOSING_STOP}",
        "GCT".repeat(8),
        "GCT".repeat(6),
        "GCT".repeat(3),
        TAIL.concat()
    );
    let mut tx = selenoprotein(&cds, Some(&UTR.concat()));
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.peptide = Some(format!(
        "M{}U{}U{}TLS",
        "A".repeat(8),
        "A".repeat(6),
        "A".repeat(3)
    ));
    vefc.seq_edits.push(SeqEdit {
        start: 17,
        end: 17,
        alt_seq: "U".to_string(),
    });
    let row = deletion(&tx, 49, 49);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.49del"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Sec17GlufsTer6")
    );
}

/// The peptide of a model that holds the stop codon of residue 10 as a stop
/// and reads on past it: 23 letters, the tenth a stop.
fn peptide_that_holds_the_stop() -> String {
    format!("M{}*{}TLS", "A".repeat(8), "A".repeat(10))
}

/// The compact test transcript with `cds` laid from its first coding base and
/// `utr` after it, as a polymorphic pseudogene's record gives it: the peptide
/// of [`peptide_that_holds_the_stop`] and no sequence edit.
fn pseudogene(cds: &str, utr: Option<&str>) -> Transcript {
    let mut tx = selenoprotein(cds, utr);
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.peptide = Some(peptide_that_holds_the_stop());
    vefc.seq_edits.clear();
    tx
}

/// The pseudogene whose coding sequence its stop codon closes, with the 3'
/// UTR after it.
fn closed_pseudogene() -> Transcript {
    pseudogene(
        &format!("{}{CLOSING_STOP}", open_coding_sequence()),
        Some(&UTR.concat()),
    )
}

/// The peptide holds the stop codon of residue 10 as a stop and names the 13
/// residues after it, so residue 15 is a residue of the reference protein.
/// Without coding base 43 the new frame meets its first stop at residue 22,
/// eight residues counted from the first changed one. Ensembl VEP writes the
/// unknown count, since the first stop of its alternate translation is
/// residue 10 and 10 - 14 is not positive.
#[test]
fn a_frameshift_past_a_stop_the_peptide_holds_counts_to_the_new_frames_stop() {
    let row = deletion(&closed_pseudogene(), 43, 43);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.43del"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15LeufsTer8")
    );
}

/// Without coding bases 43 and 44 the new frame runs past the closing stop
/// and meets its first stop at residue 25, the third to fifth bases of the 3'
/// UTR: eleven residues counted from the first changed one. No stop codon of
/// this model is read as a residue, so that codon is a stop of the new frame
/// whatever its bases.
#[test]
fn the_count_past_a_stop_the_peptide_holds_runs_into_the_three_prime_utr() {
    let row = deletion(&closed_pseudogene(), 43, 44);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.43_44del"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15CysfsTer11")
    );
}

/// The stop is sought through the coding sequence and the transcript's 3'
/// UTR and no farther. This model's coding sequence has no closing stop and
/// its 3' UTR, read from the reference, holds no stop codon in the frame that
/// the loss of coding bases 43 and 44 opens, so the count stays unknown.
#[test]
fn the_count_stays_unknown_where_the_new_frame_meets_no_stop_in_the_transcript() {
    let mut tx = pseudogene(&open_coding_sequence(), None);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let row = deletion(&tx, 43, 44);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15CysfsTer?")
    );
}

/// A frameshift before the stop the peptide holds keeps the count both
/// engines write: the first stop of the alternate translation is the new
/// frame's, at residue 22.
#[test]
fn a_frameshift_before_a_stop_the_peptide_holds_keeps_its_count() {
    let row = deletion(&closed_pseudogene(), 13, 13);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala5LeufsTer18")
    );
}

/// A peptide that holds a stop over a codon that is no stop codon is not the
/// translation of these bases, and the count stays unknown.
#[test]
fn the_count_stays_unknown_where_the_peptide_holds_a_stop_over_another_codon() {
    let mut tx = closed_pseudogene();
    let vefc = tx.vefc.as_mut().expect("a coding model");
    let mut residues = peptide_that_holds_the_stop().into_bytes();
    residues[2] = b'*';
    vefc.peptide = Some(String::from_utf8(residues).expect("ASCII residues"));
    let row = deletion(&tx, 43, 43);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15LeufsTer?")
    );
}

/// A model that reads the stop codon of residue 5 as a selenocysteine and
/// holds the stop codon of residue 10, another codon, as a stop. Without
/// coding base 43 the new frame's first stop, at residue 22, has the bases of
/// neither, and the count is eight. Without coding bases 43 and 44 its first
/// stop, in the 3' UTR, has the bases of the codon read as a selenocysteine:
/// whether that codon ends the protein there is not known from the model, and
/// the count stays unknown.
#[test]
fn a_new_stop_with_the_bases_of_a_codon_read_as_a_residue_stays_unknown_beside_a_held_stop() {
    let cds = format!(
        "ATG{}TGA{}TAG{}{}{CLOSING_STOP}",
        "GCT".repeat(3),
        "GCT".repeat(4),
        "GCT".repeat(10),
        TAIL.concat()
    );
    let mut tx = selenoprotein(&cds, Some(&UTR.concat()));
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.peptide = Some(format!(
        "M{}U{}*{}TLS",
        "A".repeat(3),
        "A".repeat(4),
        "A".repeat(10)
    ));
    vefc.seq_edits = vec![SeqEdit {
        start: 5,
        end: 5,
        alt_seq: "U".to_string(),
    }];
    let counted = deletion(&tx, 43, 43);
    assert_eq!(
        counted.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15LeufsTer8")
    );
    let unknown = deletion(&tx, 43, 44);
    assert_eq!(
        unknown.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala15CysfsTer?")
    );
}

/// The coding sequence of a model whose peptide holds two stops, without its
/// closing stop: a start codon, three alanine codons, a stop codon at residue
/// 5, four alanine codons, another stop codon at residue 10, ten alanine
/// codons and the three codons of `TAIL`: 23 codons, 69 bases.
fn open_coding_sequence_with_two_stops() -> String {
    format!(
        "ATG{}TGA{}TAG{}{}",
        "GCT".repeat(3),
        "GCT".repeat(4),
        "GCT".repeat(10),
        TAIL.concat()
    )
}

/// The peptide of that coding sequence as a polymorphic pseudogene's record
/// gives it: 23 letters, the fifth and the tenth a stop.
fn peptide_that_holds_two_stops() -> String {
    format!("M{}*{}*{}TLS", "A".repeat(3), "A".repeat(4), "A".repeat(10))
}

/// The compact test transcript with `cds` laid from its first coding base and
/// `utr` after it, the peptide of [`peptide_that_holds_two_stops`] and no
/// sequence edit.
fn pseudogene_with_two_stops(cds: &str, utr: Option<&str>) -> Transcript {
    let mut tx = pseudogene(cds, utr);
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.peptide = Some(peptide_that_holds_two_stops());
    tx
}

/// The pseudogene of two held stops whose coding sequence its stop codon
/// closes, with the 3' UTR after it.
fn closed_pseudogene_with_two_stops() -> Transcript {
    pseudogene_with_two_stops(
        &format!("{}{CLOSING_STOP}", open_coding_sequence_with_two_stops()),
        Some(&UTR.concat()),
    )
}

/// The row of an insertion of `inserted` between the coding bases `after` and
/// `after + 1` of `tx`, on a reference that holds the transcript's coding
/// sequence and 3' UTR.
fn insertion(tx: &Transcript, after: usize, inserted: &str) -> HgvsNotation {
    let vefc = tx.vefc.as_ref().expect("a coding model");
    let on_genome = format!(
        "{}{}",
        vefc.translateable_seq
            .as_deref()
            .expect("a coding sequence"),
        vefc.three_prime_utr.as_deref().unwrap_or_default()
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(150, &on_genome));
    let variant = InputVariant::new(
        "1".into(),
        150 + after as u64,
        149 + after as u64,
        b"-".to_vec(),
        inserted.as_bytes().to_vec(),
    );
    generate_hgvs(&variant, tx, Some(&fasta), HgvsOptions::default())
}

/// One base inserted between codon 9 and the stop codon of residue 10 leaves
/// the nine residues before it and makes the stop the first changed residue,
/// a leucine. Both engines write such a frameshift as an extension named at
/// the stop. The new frame runs past the closing stop and meets its first
/// stop at residue 26, the third to fifth bases of the 3' UTR: sixteen
/// residues from the lost stop. Ensembl VEP writes the unknown count, since
/// the first stop of its alternate translation is the stop of residue 5,
/// which the peptide holds before the change, and 5 - 9 is not positive
/// (TranscriptVariationAllele.pm:2503 and 2508, written by 2007 to 2010).
#[test]
fn a_frameshift_that_changes_a_later_stop_the_peptide_holds_counts_its_extension() {
    let row = insertion(&closed_pseudogene_with_two_stops(), 27, "C");
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.27_28insC"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10LeuextTer16")
    );
}

/// Without coding base 27, the last base of codon 9, that codon takes the
/// first base of the stop codon and stays an alanine codon, and the stop of
/// residue 10 is the first changed residue, an arginine. The new frame meets
/// its first stop at residue 22, inside the coding sequence: twelve residues
/// from the lost stop.
#[test]
fn the_extension_of_a_frameshift_counts_to_a_stop_inside_the_coding_sequence() {
    let row = deletion(&closed_pseudogene_with_two_stops(), 27, 27);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10ArgextTer12")
    );
}

/// The stop codon that closes the coding sequence, behind the stop the
/// peptide holds at residue 10. One base inserted before it makes it the
/// first changed residue, a leucine, and the new frame meets its first stop
/// two residues on, in the 3' UTR. Ensembl VEP's alternate translation ends at
/// residue 10, and 10 - 23 is not positive.
#[test]
fn the_extension_of_a_frameshift_at_the_closing_stop_counts_behind_a_held_stop() {
    let row = insertion(&closed_pseudogene(), 69, "C");
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter24LeuextTer2")
    );
}

/// A frameshift whose first changed residue is the first stop the peptide
/// holds keeps the count both engines write: no stop stands before it, so
/// the first stop of the alternate translation is the new frame's, at
/// residue 26, and 26 - 9 less the stop itself is sixteen.
#[test]
fn a_frameshift_at_the_first_stop_the_peptide_holds_keeps_its_extension_count() {
    let row = insertion(&closed_pseudogene(), 27, "C");
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10LeuextTer16")
    );
}

/// The stop is sought through the coding sequence and the transcript's 3'
/// UTR and no farther. This model's coding sequence has no closing stop and
/// no 3' UTR, and the frame that the inserted base opens meets no stop codon
/// before the transcript's last base, so the count stays unknown.
#[test]
fn the_extension_count_stays_unknown_where_the_new_frame_meets_no_stop_in_the_transcript() {
    let mut tx = pseudogene_with_two_stops(&open_coding_sequence_with_two_stops(), None);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let row = insertion(&tx, 27, "C");
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10LeuextTer?")
    );
}

/// A peptide that differs from the translation of its coding sequence at a
/// codon that is no stop codon is not the protein of those bases, and the
/// count stays unknown.
#[test]
fn the_extension_count_stays_unknown_where_the_peptide_is_not_the_translation_elsewhere() {
    let mut tx = closed_pseudogene_with_two_stops();
    let vefc = tx.vefc.as_mut().expect("a coding model");
    let mut residues = peptide_that_holds_two_stops().into_bytes();
    residues[2] = b'G';
    vefc.peptide = Some(String::from_utf8(residues).expect("ASCII residues"));
    let row = insertion(&tx, 27, "C");
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10LeuextTer?")
    );
}

/// A model that reads the stop codon of residue 5 as a selenocysteine and
/// holds the stop codon of residue 10 as a stop. Whether a stop codon of the
/// new frame ends the protein is not known from such a model, and the count
/// stays unknown, as Ensembl VEP writes it.
#[test]
fn the_extension_count_stays_unknown_behind_a_stop_codon_read_as_a_residue() {
    let mut tx = selenoprotein(
        &format!("{}{CLOSING_STOP}", open_coding_sequence_with_two_stops()),
        Some(&UTR.concat()),
    );
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.peptide = Some(format!(
        "M{}U{}*{}TLS",
        "A".repeat(3),
        "A".repeat(4),
        "A".repeat(10)
    ));
    vefc.seq_edits = vec![SeqEdit {
        start: 5,
        end: 5,
        alt_seq: "U".to_string(),
    }];
    let row = insertion(&tx, 27, "C");
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10LeuextTer?")
    );
}
