// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Targeted concordance tests against Perl VEP.
//!
//! Each test targets one divergence shape and asserts the **Perl-correct**
//! behavior. A test marked `#[ignore]` pins a divergence this crate does not
//! reproduce; run those with:
//! ```sh
//! cargo test -p vep-effects -- --ignored
//! ```
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`); `TranscriptMapper` is in ensembl core release/115.

use crate::consequences::{calculate_consequences, EffectsConfig};
use crate::test_helpers::{make_test_transcript, make_test_transcript_with_flags};
use tempfile::tempdir;
use vep_core::consequence::Consequence;
use vep_core::coordinate::Strand;
use vep_core::transcript::*;
use vep_core::variant::InputVariant;

// Transcript fixture helpers

/// Clone `make_test_transcript()` and add a 3' UTR sequence containing a stop codon downstream.
/// A downstream stop in a cached UTR is the condition Perl never sees in cache mode.
fn make_transcript_with_utr() -> Transcript {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        // A 3' UTR with a stop codon (TAA) ~30bp in.
        let utr = "GCAGCAGCAGCAGCAGCAGCAGCAGCAGCTAAGCAGCAGCAGCAGCAGCA";
        vefc.three_prime_utr = Some(utr.to_string());
    }
    tx
}

/// Clone `make_test_transcript()` and replace the CDS with one that ends in a
/// real terminal stop codon. Several CDS-end concordance cases need a biologically
/// valid stop-codon context rather than the default incomplete terminal codon.
fn make_transcript_with_terminal_stop() -> Transcript {
    let mut tx = make_test_transcript();

    let mut cds = String::with_capacity(849);
    cds.push_str("ATG"); // codon 1: Met
    cds.push_str("GCT"); // codon 2: Ala
    cds.push_str("GGA"); // codon 3: Gly
    cds.push_str("AAA"); // codon 4: Lys
    cds.push_str("TTC"); // codon 5: Phe
    cds.push_str("GAT"); // codon 6: Asp
    while cds.len() < 846 {
        cds.push_str("GCT");
    }
    cds.truncate(846);
    cds.push_str("TAA"); // terminal stop codon at CDS positions 847-849

    if let Some(ref mut vefc) = tx.vefc {
        vefc.translateable_seq = Some(cds);
        if let Some(ref mut mapper) = vefc.mapper {
            mapper.cdna_coding_end = 899;
        }
    }

    tx.cdna_coding_end = Some(899);
    tx.coding_region_end = Some(25_004_298);
    tx.translation_end = Some(25_004_298);
    if let Some(ref mut translation) = tx.translation {
        translation.end = 299;
    }

    tx
}

/// Clone `make_test_transcript()` and patch the exon-boundary codon spanning
/// CDS positions 250-252 so that a SNV at genomic 25_000_299 changes `CAG` to `TAG`.
fn make_transcript_with_boundary_stop_gain() -> Transcript {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut cds) = vefc.translateable_seq {
            let mut bytes = cds.clone().into_bytes();
            bytes[249] = b'C';
            bytes[250] = b'A';
            bytes[251] = b'G';
            *cds = String::from_utf8(bytes).expect("patched CDS should remain valid ASCII");
        }
    }
    tx
}

/// Clone `make_test_transcript()` and patch codon 3 (CDS positions 7-9) to `TCT`.
/// This creates a serine codon where a pure insertion inside it can turn the
/// first altered codon into `TAG` and yield a leading stop.
fn make_transcript_with_codon3_tct() -> Transcript {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut cds) = vefc.translateable_seq {
            let mut bytes = cds.clone().into_bytes();
            bytes[6] = b'T';
            bytes[7] = b'C';
            bytes[8] = b'T';
            *cds = String::from_utf8(bytes).expect("patched CDS should remain valid ASCII");
        }
    }
    tx
}

/// Clone `make_test_transcript()` and patch CDS positions 5-7 to `ATG`.
/// This creates a downstream in-frame ATG that should not trigger
/// start_retained_variant once start_lost has already fired for a start-codon indel.
fn make_transcript_with_downstream_inframe_atg() -> Transcript {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut cds) = vefc.translateable_seq {
            let mut bytes = cds.clone().into_bytes();
            bytes[4] = b'A';
            bytes[5] = b'T';
            bytes[6] = b'G';
            *cds = String::from_utf8(bytes).expect("patched CDS should remain valid ASCII");
        }
    }
    tx
}

/// Build a non-coding miRNA transcript with mature-miRNA subfeature at cDNA 50-100.
fn make_mirna_transcript() -> Transcript {
    let mut tx = make_non_coding_transcript();
    tx.biotype = "miRNA".into();
    tx.attributes.push(Attribute {
        code: "miRNA".to_string(),
        value: "50-100".to_string(),
    });
    tx
}

/// Clone `make_test_transcript()` and convert to a non-coding transcript.
/// Clears translation, zeros CDS coords, sets biotype to "processed_transcript".
fn make_non_coding_transcript() -> Transcript {
    let mut tx = make_test_transcript();
    tx.biotype = "processed_transcript".into();
    tx.translation = None;
    tx.protein_id = None;
    tx.cdna_coding_start = None;
    tx.cdna_coding_end = None;
    tx.coding_region_start = None;
    tx.coding_region_end = None;
    tx.translation_start = None;
    tx.translation_end = None;
    if let Some(ref mut vefc) = tx.vefc {
        vefc.translateable_seq = None;
        vefc.peptide = None;
        if let Some(ref mut mapper) = vefc.mapper {
            mapper.cdna_coding_start = 0;
            mapper.cdna_coding_end = 0;
        }
    }
    tx
}

/// Clone `make_test_transcript()` and flip to reverse strand.
/// Swaps exon coordinates so that the transcript runs 3'->5' on forward genome.
fn make_reverse_strand_transcript() -> Transcript {
    let mut tx = make_test_transcript();
    tx.strand = Strand::Reverse;
    // The mapper pairs already encode the mapping; only strand and ori flip.
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut mapper) = vefc.mapper {
            for pair in &mut mapper.exon_coord_mapper.pairs {
                pair.ori = -1;
            }
        }
    }
    tx
}

/// Shift the standard test transcript down to small chr1 coordinates so tests can
/// build compact temporary FASTA fixtures without needing multi-megabase reference files.
fn make_compact_test_transcript() -> Transcript {
    let mut tx = make_test_transcript();
    let shift = 24_999_900;

    tx.chr = "1".into();
    tx.start -= shift;
    tx.end -= shift;
    for exon in &mut tx.exons {
        exon.start -= shift;
        exon.end -= shift;
    }
    for intron in &mut tx.introns {
        intron.start -= shift;
        intron.end -= shift;
    }
    if let Some(coord) = tx.coding_region_start.as_mut() {
        *coord -= shift;
    }
    if let Some(coord) = tx.coding_region_end.as_mut() {
        *coord -= shift;
    }
    if let Some(coord) = tx.translation_start.as_mut() {
        *coord -= shift;
    }
    if let Some(coord) = tx.translation_end.as_mut() {
        *coord -= shift;
    }
    if let Some(ref mut vefc) = tx.vefc {
        for exon in &mut vefc.sorted_exons {
            exon.start -= shift;
            exon.end -= shift;
        }
        for intron in &mut vefc.introns {
            intron.start -= shift;
            intron.end -= shift;
        }
        if let Some(ref mut mapper) = vefc.mapper {
            for pair in &mut mapper.exon_coord_mapper.pairs {
                pair.to_start -= shift;
                pair.to_end -= shift;
            }
        }
    }

    tx
}

fn make_compact_reverse_strand_transcript() -> Transcript {
    let mut tx = make_compact_test_transcript();
    tx.strand = Strand::Reverse;
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut mapper) = vefc.mapper {
            for pair in &mut mapper.exon_coord_mapper.pairs {
                pair.ori = -1;
            }
        }
    }
    tx
}

fn load_test_fasta(chr: &str, sequence: &str) -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
    let dir = tempdir().expect("tempdir");
    let fa = dir.path().join("ref.fa");
    let fai = dir.path().join("ref.fa.fai");
    let fasta_contents = format!(">{chr}\n{sequence}\n");
    std::fs::write(&fa, fasta_contents).expect("write FASTA");
    let header_offset = chr.len() + 2;
    std::fs::write(
        &fai,
        format!(
            "{chr}\t{}\t{header_offset}\t{}\t{}\n",
            sequence.len(),
            sequence.len(),
            sequence.len() + 1
        ),
    )
    .expect("write FASTA index");

    let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("load indexed FASTA");
    (dir, fasta)
}

/// Frameshift + stop_gained: a cached 3' UTR must not supply a downstream stop,
/// because in cache mode Perl's UTR is undef.
///
/// Expected (Perl): `frameshift_variant` only, not `stop_gained`.
#[test]
fn concordance_frameshift_ignores_cached_utr_stop_for_stop_gained() {
    let tx = make_transcript_with_utr();
    let config = EffectsConfig::default();

    // 1bp insertion near the end of the CDS to trigger a frameshift.
    // CDS ends at cDNA 900 = genomic 25_004_299.
    // Insert at genomic 25_004_280 (near CDS end, within exon 3).
    let variant = InputVariant::new(
        "21".into(),
        25_004_281,
        25_004_280, // end < start = VEP insertion convention
        b"-".to_vec(),
        b"A".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "Should contain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopGained),
        "should NOT contain stop_gained when UTR provides downstream stop \
         (Perl doesn't see cached UTR in offline mode)"
    );
}

/// Frameshift stop_gained reads only Perl's codon window: a stop the shifted
/// frame spells far downstream in the CDS is never translated. The insertion
/// here sits at a codon boundary, so the window is the inserted base alone
/// (`X`), as `compute_codon_window_peptide_alleles` builds it.
///
/// Expected (Perl): `frameshift_variant` only, not `stop_gained`.
#[test]
fn concordance_frameshift_stop_gained_ignores_stop_outside_codon_window() {
    // Build a transcript where the shifted reading frame has a stop codon far
    // downstream (alt peptide position 9), well outside the codon window.
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut cds) = vefc.translateable_seq {
            // Embed nucleotides at CDS positions 53-55 (0-indexed) that produce
            // a TAA stop codon in the +1 shifted reading frame when a 1bp insertion
            // is placed at CDS position 30.
            let mut bytes = cds.clone().into_bytes();
            bytes[53] = b'T';
            bytes[54] = b'A';
            bytes[55] = b'A';
            *cds = String::from_utf8(bytes).expect("patched CDS should remain valid ASCII");
        }
    }
    let config = EffectsConfig::default();

    // 1bp insertion at CDS position 30 (codon 10), within exon 1.
    // CDS starts at genomic 25_000_050 (cdna_coding_start=51, exon 1 starts at 25_000_000).
    // CDS position 30 = genomic 25_000_079.
    // VEP insertion convention: start = pos+1, end = pos.
    let variant = InputVariant::new(
        "21".into(),
        25_000_080,
        25_000_079, // end < start = VEP insertion convention
        b"-".to_vec(),
        b"A".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "Should contain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopGained),
        "should NOT contain stop_gained when the stop is only found far downstream \
         in the shifted reading frame, outside the codon window"
    );
}

/// Frameshift stop_gained for a 1bp insertion after the last base of a codon:
/// Perl's codon window is the inserted base alone, which
/// `TranscriptVariationAllele::peptide` resolves as `X`, so the stop the shifted
/// frame spells in the next codon is never read.
///
/// Expected (Perl): `frameshift_variant` only, not `stop_gained`.
#[test]
fn concordance_frameshift_insertion_after_codon_does_not_read_next_codon_stop() {
    // Build a transcript where a 1bp insertion at the last position of a codon
    // creates a stop in the shifted reading frame of the next codon.
    // CDS: ...CGT AAC... (codons N, N+1). Insert T at position 2 of codon N:
    // Alt: ...CGT T AAC... → the shifted frame spells TAA in the next codon, but
    // the window Perl translates is the inserted `T` alone (`X`).
    let mut tx = make_test_transcript();
    if let Some(vefc) = tx.vefc.as_mut() {
        if let Some(seq) = vefc.translateable_seq.as_mut() {
            // Set CDS at codon 10 (positions 27-29, 0-indexed) to CGT (R)
            // and codon 11 (positions 30-32) to AAC (N).
            // A 1bp insertion of T after position 29 produces:
            //   CGT [T] AAC → CGT = R, TAA = *, C?? = X → R*X
            let mut bytes = seq.clone().into_bytes();
            bytes[27] = b'C';
            bytes[28] = b'G';
            bytes[29] = b'T';
            bytes[30] = b'A';
            bytes[31] = b'A';
            bytes[32] = b'C';
            *seq = String::from_utf8(bytes).expect("patched CDS should remain valid ASCII");
        }
    }
    let config = EffectsConfig::default();

    // CDS position 30 (1-based) = codon 10, position 2 (last base of codon).
    // CDS starts at genomic 25_000_050 in the test transcript.
    // Genomic position for CDS pos 30: 25_000_050 + 29 = 25_000_079.
    // VEP insertion convention: start = pos+1, end = pos.
    let variant = InputVariant::new(
        "21".into(),
        25_000_080,
        25_000_079, // end < start = VEP insertion convention
        b"-".to_vec(),
        b"T".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "Should contain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopGained),
        "1bp insertion at last codon position should NOT \
         contain stop_gained: the stop is in the shifted-frame NEXT codon, \
         outside the window Perl translates. Got: {:?}",
        tc.consequences
    );
}

/// Splice region for insertions.
///
/// Expected (Perl): no `splice_region_variant` for a pure insertion 3bp into exon.
/// Perl's `_get_differing_regions` returns an inverted span for insertions, and
/// `_intron_overlap()`'s raw overlap check therefore does not count this exonic
/// insertion point as overlapping the acceptor-side 3bp window.
#[test]
fn concordance_exonic_insertion_three_bases_into_exon_no_splice_region() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Insertion 3bp into exon 2 from the acceptor side.
    // Exon 2 starts at 25_002_000. Position 25_002_002 = 3bp in.
    let variant = InputVariant::new(
        "21".into(),
        25_002_003,
        25_002_002, // insertion between pos 2 and 3 of exon 2
        b"-".to_vec(),
        b"A".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "pure exonic insertion 3bp from the boundary should NOT include \
         splice_region_variant under Perl's raw overlap semantics, got: {:?}",
        tc.consequences
    );
}

/// Large insertion must not over-call splice_region: its cDNA position is within
/// 3bp of an exon boundary but its genomic anchor (narrowed via
/// `exonic_splice_differing_region`) is not, and Perl's `_intron_overlap()` uses
/// the narrowed span.
///
/// Expected (Perl): no splice_region_variant.
#[test]
fn concordance_large_insertion_far_from_exon_boundary_no_splice_region() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Large coding insertion at genomic 25_000_290 (10bp from exon 1 end at 25_000_299).
    // Exon 1 ends at 25_000_299, intron 1 starts at 25_000_300.
    // The exonic donor window is [25_000_297, 25_000_299] (3bp before intron start).
    // Position 25_000_290 is 10bp away from the exon boundary, well outside the
    // 3bp window. After narrowing, the effective span is (25_000_290, 25_000_290)
    // for a pure insertion.
    let variant = InputVariant::new(
        "21".into(),
        25_000_291, // insertion between 290 and 291
        25_000_290,
        b"-".to_vec(),
        b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(), // 30bp insertion
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "large insertion far from exon boundary should NOT get \
         splice_region_variant, got: {:?}",
        tc.consequences
    );
}

/// Complex indel overcalls splice_region when full ref span
/// reaches the exonic splice window but the actual insertion point doesn't.
///
/// When `prefix + suffix >= ref_allele.len()`, `exonic_splice_differing_region()`
/// must narrow to the insertion point `(lo + prefix, hi - suffix)` as Perl's
/// `_get_differing_regions` does; the full ref span can reach the exonic splice
/// window when the net-insertion point is far away.
///
/// Example: `ref=GCCTAAG, alt=GCCCTAAG` (7bp→8bp, net +1bp after prefix=2, suffix=4)
/// at exon positions 25_000_293-25_000_299 (ending 1bp from intron boundary). The full
/// span [25_000_293, 25_000_299] overlaps the donor exonic window [25_000_297, 25_000_299],
/// but the insertion point [25_000_295, 25_000_296] does not.
///
/// Expected (Perl): no splice_region_variant.
#[test]
fn concordance_complex_indel_consumed_ref_no_splice_region() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Exon 1 ends at 25_000_299. Intron 1 starts at 25_000_300.
    // Donor exonic window: [25_000_297, 25_000_299] (3bp before intron start).
    //
    // Complex indel spanning 7 exonic bases ending at the splice boundary:
    //   Genomic: 25_000_293-25_000_299 (ref = GCCTAAG)
    //   Alt:     GCCCTAAG (8bp, net +1bp insertion of 'C' after prefix 'GCC')
    //
    // Prefix: GCC (3 bytes match)
    // Suffix: TAAG (4 bytes match from the end)
    // prefix + suffix = 7 = ref_len → entire ref consumed
    //
    // Perl's _get_differing_regions: start+3=25_000_296, end-4=25_000_295
    //   → inverted (25_000_296, 25_000_295) = insertion point
    //   → does not overlap [25_000_297, 25_000_299]
    //
    // The full span (25_000_293, 25_000_299) would overlap.
    let variant = InputVariant::new(
        "21".into(),
        25_000_293,
        25_000_299,
        b"GCCTAAG".to_vec(),
        b"GCCCTAAG".to_vec(), // 7bp → 8bp, net +1bp
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "complex indel with consumed ref should NOT overcall \
         splice_region_variant when the insertion point is outside the exonic \
         splice window, got: {:?}",
        tc.consequences
    );
}

/// Complex indel should fire splice_region when
/// its differing region is within the exonic splice window.
///
/// The narrowing above must not suppress a legitimate call: a 2bp complex indel
/// two bases before the intron, whose XOR region falls within the 3bp window.
#[test]
fn concordance_complex_indel_consumed_ref_splice_region_when_near_boundary() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Exon 1 ends at 25_000_299. Intron 1 starts at 25_000_300.
    // Donor exonic window: [25_000_297, 25_000_299].
    //
    // Complex indel at the two exonic bases before the last one:
    //   Genomic: 25_000_297-25_000_298 (ref = AG)
    //   Alt:     ACG (3bp, net +1bp)
    //
    // Perl's _get_differing_regions XORs the alleles as long as the longer one:
    // A^A = 0, G^C != 0, ""^G != 0, so the single region is offsets 1..2, genomic
    // (25_000_298, 25_000_299): inside [25_000_297, 25_000_299] and short of the
    // donor dinucleotide at 25_000_300, so splice_region fires and nothing
    // suppresses it. One base further right the region would reach the donor,
    // which Perl calls splice_donor_variant and suppresses splice_region for.
    let variant = InputVariant::new(
        "21".into(),
        25_000_297,
        25_000_298,
        b"AG".to_vec(),
        b"ACG".to_vec(), // 2bp → 3bp, net +1bp
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "complex indel with insertion point inside the \
         exonic splice window should include splice_region_variant, got: {:?}",
        tc.consequences
    );
}

/// Polypyrimidine tract: Perl narrows the overlap with `_get_differing_regions`,
/// so a substitution's differing region can miss the polypyrimidine window that
/// its full span overlaps.
///
/// Expected (Perl): `splice_region_variant` only, not `splice_polypyrimidine_tract_variant`.
#[test]
fn concordance_wide_substitution_differing_region_misses_polypyrimidine_tract() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Wide same-length intronic substitution whose full span reaches the acceptor-side
    // polypyrimidine tract, but whose only differing base is at donor+6.
    //
    // Intron 1: 25_000_300..25_001_999
    // donor-side splice_region window: 25_000_302..25_000_307
    // donor-side donor_region window: 25_000_302..25_000_305
    // acceptor-side polypyrimidine window: 25_001_983..25_001_997
    //
    // The full span (25_000_306..25_001_985) overlaps both splice_region and PPT,
    // but after trimming shared suffixes the effective differing span is just
    // 25_000_306, which should yield splice_region only.
    let len = (25_001_985u64 - 25_000_306u64 + 1) as usize;
    let ref_allele = vec![b'A'; len];
    let mut alt_allele = ref_allele.clone();
    alt_allele[0] = b'C';

    let variant = InputVariant::new("21".into(), 25_000_306, 25_001_985, ref_allele, alt_allele);

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "Should contain splice_region_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences
            .contains(&Consequence::SplicePolypyrimidineTractVariant),
        "should NOT contain splice_polypyrimidine_tract_variant \
         (Perl's differing-region narrows the effective span)"
    );
}

/// Stop/missense at an exon-boundary codon.
///
/// Expected (Perl): `stop_gained` (not `missense_variant`).
#[test]
fn concordance_exon_boundary_codon_snv_is_stop_gained() {
    let tx = make_transcript_with_boundary_stop_gain();
    let config = EffectsConfig::default();

    // SNV at the last CDS base of exon 1 (exon boundary).
    // Exon 1 CDS: genomic 25_000_050 (cDNA 51) to 25_000_299 (cDNA 300).
    // Last coding base of exon 1 = 25_000_299. This is codon position 250,
    // within codon 84 (position 1). The codon spans the exon boundary: 1 base
    // in exon 1 + 2 bases in exon 2. Patch the reference codon to `CAG` so a
    // SNV at the boundary yields `TAG` (stop_gained).
    let variant = InputVariant::new(
        "21".into(),
        25_000_299,
        25_000_299,
        b"C".to_vec(), // boundary base of patched `CAG`
        b"T".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    // The assertion is agreement with Perl on the consequence for an
    // exon-boundary codon.
    assert!(
        tc.consequences.contains(&Consequence::StopGained),
        "exon-boundary SNV should produce stop_gained (not missense), \
         got: {:?}",
        tc.consequences
    );
}

/// Stop-retained / inframe classification: Perl translates the alt CDS rather
/// than reading position alone.
///
/// Expected (Perl): `inframe_insertion` + `stop_retained_variant`.
#[test]
fn concordance_inframe_insertion_near_stop_is_stop_retained() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();

    // 3bp inframe insertion at the penultimate codon (near CDS end).
    // CDS: cDNA 51-899 = 849bp. Last codon at cDNA 897-899 (genomic 25_004_296-25_004_298).
    // Penultimate codon: cDNA 894-896 (genomic 25_004_293-25_004_295).
    // Insert 3bp (inframe) at genomic 25_004_295.
    let variant = InputVariant::new(
        "21".into(),
        25_004_296,
        25_004_295, // insertion
        b"-".to_vec(),
        b"GCT".to_vec(), // 3bp = inframe
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::InframeInsertion),
        "should contain inframe_insertion, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StopRetainedVariant),
        "should contain stop_retained_variant \
         (Perl translates alt CDS and confirms stop is preserved)"
    );
}

/// terminal stop-codon frameshift should keep
/// `frameshift_variant` when Perl also reports `stop_lost`.
#[test]
fn concordance_frameshift_stop_lost_keeps_frameshift() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();

    // Delete one base from the terminal stop codon (TAA at genomic
    // 25_004_296-25_004_298). Perl reports frameshift_variant + stop_lost.
    let variant = InputVariant::new(
        "21".into(),
        25_004_297,
        25_004_297,
        b"A".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "terminal stop frameshift should keep frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StopLost),
        "terminal stop frameshift should also contain stop_lost, got: {:?}",
        tc.consequences
    );
}

/// UTR term for spanning deletions, per Perl's `within_cdna`.
///
/// Expected (Perl): includes `3_prime_UTR_variant`.
#[test]
fn concordance_utr_spanning_deletion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Deletion from coding region into 3' UTR.
    // CDS ends at genomic 25_004_299 (cDNA 900). 3' UTR starts at 25_004_300.
    // Delete from 25_004_295 to 25_004_305 (spans CDS end into UTR).
    let variant = InputVariant::new(
        "21".into(),
        25_004_295,
        25_004_305,
        b"ACGTACGTACG".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "spanning deletion into 3' UTR should include \
         3_prime_UTR_variant, got: {:?}",
        tc.consequences
    );
}

/// Some cache transcripts carry mapper coding bounds but omit the top-level
/// genomic coding_region_start/end fields. Perl still adds the UTR term.
#[test]
fn concordance_utr_spanning_deletion_mapper_only_bounds() {
    let mut tx = make_test_transcript();
    let config = EffectsConfig::default();
    tx.coding_region_start = None;
    tx.coding_region_end = None;
    tx.translation_start = None;
    tx.translation_end = None;

    let variant = InputVariant::new(
        "21".into(),
        25_004_295,
        25_004_305,
        b"ACGTACGTACG".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "mapper-only bounds: spanning deletion into 3' UTR should include \
         3_prime_UTR_variant even when genomic coding bounds are absent, got: {:?}",
        tc.consequences
    );
}

/// UTR term for cds_end_NF transcript (forward strand).
/// Perl adds 3_prime_UTR_variant when a deletion extends past the transcript
/// boundary on a transcript with incomplete CDS end (cds_end_NF flag).
#[test]
fn concordance_utr_spanning_deletion_cds_end_nf() {
    let mut tx = make_test_transcript_with_flags(&["cds_end_NF"]);
    let config = EffectsConfig::default();
    tx.coding_region_end = Some(tx.end);
    tx.translation_end = Some(tx.end);

    let variant = InputVariant::new(
        "21".into(),
        25_005_990,
        25_008_000,
        b"ACGTACGT".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some());
    let tc = result.unwrap();
    assert!(
        tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "cds_end_NF: deletion past transcript end should add 3_prime_UTR_variant, got: {:?}",
        tc.consequences
    );
}

/// 3_prime_UTR_variant for reverse-strand cds_end_NF transcript.
///
/// With `cds_end_NF` on the reverse strand the CDS extends to the transcript 3'
/// end = low genomic = tx.start.
/// A deletion extending past tx.start (into the "past-3-prime" region) still
/// gets `3_prime_UTR_variant` in Perl via `_before_coding` + non-normalizing
/// overlap: the empty normalized UTR range `(tx.start, crs-1=tx.start-1)`
/// still matches when the variant extends below tx.start.
#[test]
fn concordance_utr_spanning_deletion_cds_end_nf_reverse() {
    let mut tx = make_test_transcript_with_flags(&["cds_end_NF"]);
    let config = EffectsConfig::default();
    tx.strand = vep_core::coordinate::Strand::Reverse;
    tx.coding_region_start = Some(tx.start);
    tx.translation_start = Some(tx.start);

    let variant = InputVariant::new(
        "21".into(),
        24_998_000,
        25_000_010,
        b"ACGTACGT".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some());
    let tc = result.unwrap();
    assert!(
        tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "cds_end_NF reverse: deletion below transcript start should add \
         3_prime_UTR_variant, got: {:?}",
        tc.consequences
    );
}

/// 5_prime_UTR_variant for reverse-strand cds_start_NF transcript.
///
/// Mirror scenario: `cds_start_NF` on reverse strand means CDS 5' end extends
/// to transcript 5' end = high genomic = tx.end. A deletion extending past
/// tx.end should get `5_prime_UTR_variant` via the NF boundary check.
#[test]
fn concordance_utr_spanning_deletion_cds_start_nf_reverse() {
    let mut tx = make_test_transcript_with_flags(&["cds_start_NF"]);
    let config = EffectsConfig::default();
    tx.strand = vep_core::coordinate::Strand::Reverse;
    tx.coding_region_end = Some(tx.end);
    tx.translation_end = Some(tx.end);

    let variant = InputVariant::new(
        "21".into(),
        25_005_990,
        25_008_000,
        b"ACGTACGT".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some());
    let tc = result.unwrap();
    assert!(
        tc.consequences.contains(&Consequence::FivePrimeUtrVariant),
        "cds_start_NF reverse: deletion above transcript end should add \
         5_prime_UTR_variant, got: {:?}",
        tc.consequences
    );
}

/// A spanning-deletion MNV must not emit stop_lost even when one endpoint lands
/// in CDS and the other in 3' UTR.
///
/// For an MNV spanning from the CDS into the 3' UTR, Perl's
/// `_ins_del_stop_altered` returns 0 (length unchanged), and `_get_peptide_alleles`
/// returns empty (tv_tr_end undefined because the variant extends past CDS).
/// Expected: coding_sequence_variant + 3_prime_UTR_variant, no stop_lost.
#[test]
fn concordance_mnv_spanning_utr_no_stop_lost() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // CDS ends at genomic 25_004_299 (cDNA 900). Create a 2bp MNV spanning
    // the CDS end into the 3' UTR: pos 25_004_299 (last CDS base) and 25_004_300
    // (first UTR base). Ref "GG", alt "TT": same length, length-preserving.
    let variant = InputVariant::new(
        "21".into(),
        25_004_299,
        25_004_300,
        b"GG".to_vec(),
        b"TT".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some());
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::StopLost),
        "MNV spanning coding→UTR should not emit stop_lost, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "MNV spanning coding→UTR should emit 3_prime_UTR_variant, got: {:?}",
        tc.consequences
    );
}

/// A large pure deletion extending past the stop codon into the 3' UTR
/// emits stop_lost + 3_prime_UTR_variant (matching Perl's
/// `_ins_del_stop_altered` via length-shorter branch).
#[test]
fn concordance_large_deletion_stop_lost_far_from_cds_start() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();

    // Terminal stop is at CDS 847-849 (genomic 25_004_296-25_004_298 on exon 3).
    // Place a deletion starting well before the stop codon (cds_pos far from
    // cds_len) that spans through the stop codon and well into 3' UTR.
    // This still emits stop_lost; a `cds_pos >= cds_len - 50` proximity gate
    // would reject it.
    //
    // Exon 3 starts at 25_004_000 (cDNA 601, CDS 551). Delete from cDNA ~750
    // (CDS ~700) to well past CDS end: genomic 25_004_150 (CDS 701) to
    // 25_005_000 (701 bp into the 3' UTR).
    let ref_seq = vec![b'A'; 851]; // 25_004_150..25_005_000 inclusive = 851bp
    let variant = InputVariant::new("21".into(), 25_004_150, 25_005_000, ref_seq, b"-".to_vec());

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some());
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::StopLost),
        "large deletion past stop codon should emit stop_lost, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "large deletion past stop codon should emit 3_prime_UTR_variant, got: {:?}",
        tc.consequences
    );
}

/// A spanning deletion on a cds_start_NF transcript must not emit start_lost:
/// Perl's `_overlaps_start_codon` early-returns 0 when cds_start_NF is set, and
/// Perl emits `coding_sequence_variant,5_prime_UTR_variant`.
#[test]
fn concordance_spanning_deletion_cds_start_nf_no_start_lost() {
    let mut tx = make_test_transcript_with_flags(&["cds_start_NF"]);
    let config = EffectsConfig::default();
    // Force CDS 5' boundary to the transcript 5' boundary (forward strand).
    tx.coding_region_start = Some(tx.start);
    tx.translation_start = Some(tx.start);
    if let Some(ref mut mapper) = tx.vefc.as_mut().and_then(|v| v.mapper.as_mut()) {
        mapper.cdna_coding_start = 1;
    }
    tx.cdna_coding_start = Some(1);
    if let Some(ref mut translation) = tx.translation {
        translation.start = 1;
    }

    // Deletion spanning from upstream (past transcript 5' end) into CDS.
    // Variant.start < tx.start so start endpoint is Upstream, end endpoint Coding.
    let ref_seq = vec![b'A'; 20];
    let variant = InputVariant::new("21".into(), 24_999_990, 25_000_009, ref_seq, b"-".to_vec());

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some());
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::StartLost),
        "cds_start_NF spanning deletion must not emit start_lost, got: {:?}",
        tc.consequences
    );
}

/// A spanning deletion on a cds_end_NF transcript must not
/// emit stop_lost. Perl's `_overlaps_stop_codon` early-returns 0 when
/// cds_end_NF is set. And `add_utr_for_overlapping_span` should still emit
/// 3_prime_UTR_variant via the NF boundary check.
#[test]
fn concordance_spanning_deletion_cds_end_nf_no_stop_lost() {
    let mut tx = make_test_transcript_with_flags(&["cds_end_NF"]);
    let config = EffectsConfig::default();
    // Force CDS 3' boundary to the transcript 3' boundary (forward strand).
    tx.coding_region_end = Some(tx.end);
    tx.translation_end = Some(tx.end);

    // Deletion spanning from CDS body past transcript 3' end.
    let ref_seq = vec![b'A'; 20];
    let variant = InputVariant::new("21".into(), 25_005_990, 25_006_009, ref_seq, b"-".to_vec());

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some());
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::StopLost),
        "cds_end_NF spanning deletion must not emit stop_lost, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "cds_end_NF spanning deletion should still emit 3_prime_UTR_variant via \
         NF boundary check, got: {:?}",
        tc.consequences
    );
}

/// Non-coding exon vs transcript naming: Perl uses
/// `non_coding_transcript_exon_variant` for a spanning deletion crossing into an
/// exon.
///
/// Expected (Perl): `non_coding_transcript_exon_variant`.
#[test]
fn concordance_non_coding_spanning_deletion_into_exon_is_exon_variant() {
    let tx = make_non_coding_transcript();
    let config = EffectsConfig::default();

    // Spanning deletion in non-coding transcript crossing from intron into exon.
    // Intron 1: 25_000_300-25_001_999. Exon 2: 25_002_000-25_002_299.
    // Delete from intron 1 into exon 2: 25_001_990 to 25_002_010.
    let variant = InputVariant::new(
        "21".into(),
        25_001_990,
        25_002_010,
        b"ACGTACGTACGTACGTACGTAC".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant),
        "spanning deletion in non-coding tx crossing into exon should \
         produce non_coding_transcript_exon_variant, got: {:?}",
        tc.consequences
    );
}

/// Splice donor/acceptor sub-terms for spanning deletions.
///
/// Expected (Perl): both `splice_donor_variant` and `splice_acceptor_variant`.
#[test]
fn concordance_spanning_deletion_emits_donor_and_acceptor() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Large deletion spanning from exon 1 through intron 1 into intron 2.
    // This crosses: exon 1 end (donor of intron 1), intron 1, exon 2, intron 2 start (acceptor of intron 2).
    // Exon 1: ..25_000_299, Intron 1: 25_000_300..25_001_999, Exon 2: 25_002_000..25_002_299,
    // Intron 2: 25_002_300..25_003_999.
    // Delete from 25_000_290 to 25_002_310.
    let variant = InputVariant::new(
        "21".into(),
        25_000_290,
        25_002_310,
        b"A".to_vec(), // simplified ref
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::SpliceDonorVariant),
        "large spanning deletion should include splice_donor_variant, \
         got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences
            .contains(&Consequence::SpliceAcceptorVariant),
        "large spanning deletion should include splice_acceptor_variant, \
         got: {:?}",
        tc.consequences
    );
}

/// Synonymous/missense at an exon-boundary codon.
///
/// Expected (Perl): `synonymous_variant` (not `missense_variant`).
#[test]
fn concordance_exon_boundary_wobble_snv_is_synonymous() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Wobble-position SNV at an exon boundary: the degenerate third position of a
    // codon split across exons.
    // Last coding base of exon 1 = 25_000_299 (cDNA 300).
    // CDS offset: 300 - 51 + 1 = 250. Codon 84, position 1 (250 = 83*3 + 1).
    // A synonymous call needs a third-position wobble in a split codon.
    // Use position 25_000_298 = cDNA 299 = CDS offset 249 = codon 83, position 3 (wobble).
    let variant = InputVariant::new(
        "21".into(),
        25_000_298,
        25_000_298,
        b"T".to_vec(),
        b"C".to_vec(), // wobble: GCT -> GCC = Ala -> Ala (synonymous)
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::SynonymousVariant),
        "wobble-position SNV at exon boundary should be synonymous, \
         got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::MissenseVariant),
        "should not be missense_variant (CDS splice point issue)"
    );
}

/// Protein-altering vs inframe: Perl's `protein_altering_variant` catch-all.
///
/// Expected (Perl): `protein_altering_variant` (not `inframe_insertion`).
#[test]
fn concordance_complex_inframe_indel_is_protein_altering_not_inframe_insertion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Complex inframe indel: 3bp -> 6bp at a specific codon position.
    // This is a net +3bp insertion that doesn't cleanly classify as pure insertion at codon level.
    // CDS start at genomic 25_000_050. Codon 3 at CDS offset 7-9 = genomic 25_000_056-25_000_058.
    let variant = InputVariant::new(
        "21".into(),
        25_000_056,
        25_000_058,
        b"GGA".to_vec(),
        b"TTCCAA".to_vec(), // 3bp -> 6bp = net +3 inframe
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "complex inframe indel should produce protein_altering_variant, \
         got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion),
        "should NOT be inframe_insertion \
         (Perl uses protein_altering_variant for ambiguous inframe indels)"
    );
}

/// A frameshift at the start codon also emits start_lost: Perl says
/// `frameshift_variant,start_lost`, and the codon-window bounds the later
/// start_lost check relies on can fail for an insertion at CDS position 1, so a
/// direct check after FrameshiftVariant is needed.
///
/// Expected (Perl): `frameshift_variant` + `start_lost`.
#[test]
fn concordance_frameshift_at_start_codon_emits_start_lost() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // 2bp insertion at CDS position 1 (genomic 25_000_050): A → AGG (net +2bp = frameshift).
    // This destroys the start codon reading frame.
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_050,
        b"A".to_vec(),
        b"AGG".to_vec(), // 1bp → 3bp = net +2bp frameshift at start codon
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "Should contain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StartLost),
        "Frameshift at start codon should also emit start_lost, got: {:?}",
        tc.consequences
    );
}

/// A net-insertion delins over the start codon that puts three bases in front
/// of the ATG (`ATG` to `GCAATG`) keeps the start: Perl co-emits `start_lost`
/// and `start_retained_variant`, and vep-rs keeps `start_retained_variant`.
///
/// Perl's `_ins_del_start_altered` edits the 5' UTR + CDS string and, the edit
/// being longer, compares the CDS with the edited string's tail: `GCA` is
/// inserted before the ATG, so the tail is the CDS unchanged and the start is
/// not altered, which makes `start_retained_variant` 1. `start_lost` then
/// fires through `_inv_start_altered`, which reads the triplet at the old CDS
/// start, `GCA`, without looking past it to the ATG three bases on. The
/// coding sequence is intact, so `start_lost` is the erroneous member of that
/// pair and `perl_coding_terms` drops it. `_ins_del_start_altered` reads the
/// 5' UTR: the CDS starts at cDNA 51, so the transcript carries its 50-base
/// UTR here as a cached transcript does.
#[test]
fn concordance_start_codon_delins_keeping_the_atg_is_start_retained_only() {
    let mut tx = make_test_transcript();
    tx.vefc.as_mut().unwrap().five_prime_utr = Some("GCCACC".repeat(8) + "AG");
    let config = EffectsConfig::default();

    // CDS starts at genomic 25_000_050: ATG (codon 1) -> GCAATG, net +3 bp.
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_052,
        b"ATG".to_vec(),
        b"GCAATG".to_vec(),
    );

    let tc = calculate_consequences(&variant, &tx, &config).expect("Should produce consequences");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["start_retained_variant"]);
}

/// Pure insertion (ref="-") classified as protein_altering_variant
/// when alt peptide does not contain ref peptide as prefix/suffix.
///
/// Perl's `VariationEffect::protein_altering_variant` checks peptide
/// containment before calling inframe_insertion. For pure insertions that disrupt
/// the codon reading frame alignment (e.g., inserting bases that shift the amino
/// acid at the insertion point), the alt peptide won't start/end with ref peptide.
/// Perl returns protein_altering_variant.
///
/// CDS: ATG GCT GGA AAA ... Insert "TAC" inside codon 3 (GGA = Gly).
/// Insertion between 1st and 2nd base of codon 3 (CDS pos 7, genomic 25_000_056).
/// Mapper maps start=25_000_057 to cds_pos=8, idx=7, rel=1 within codon [6..9].
/// ref codon = GGA (G), alt codon = G + TAC + GA = GTACGA → GTA(Val) CGA(Arg) = [V, R].
/// alt starts_with [G]? No (V≠G). alt ends_with [G]? No (R≠G) → ProteinAltering.
#[test]
fn concordance_pure_insertion_disrupting_codon_is_protein_altering() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Pure 3bp insertion inside codon 3 (GGA=Gly, CDS pos 7-9).
    // Insert between CDS positions 7 and 8 = between genomic 25_000_056 and 25_000_057.
    // VEP insertion convention: start = pos_after, end = pos_before.
    let variant = InputVariant::new(
        "21".into(),
        25_000_057, // start (after insertion point)
        25_000_056, // end (before insertion point)
        b"-".to_vec(),
        b"TAC".to_vec(), // 3bp = inframe, disrupts codon alignment
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "pure inframe insertion with disrupted codon should be \
         protein_altering_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion),
        "should NOT be inframe_insertion when alt peptide \
         doesn't contain ref peptide"
    );
}

/// Frameshift should not be downgraded to inframe_insertion+stop_retained
///
/// Perl treats the partial boundary codon at the insertion/CDS junction as 'X',
/// so a stop there is invisible to it; the scan limit from the frameshift branch
/// is floor division (`net_ins_nt / 3`), as in
/// `frameshift_stop_gained_in_codon_window` in `coding.rs`, and a `div_ceil(3)`
/// limit would find those shifted-frame stops.
///
/// Expected (Perl): `frameshift_variant` (not `inframe_insertion + stop_retained_variant`).
#[test]
fn concordance_frameshift_not_downgraded_to_stop_retained() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // 2bp net insertion (frameshift) near end of exon 1 (within CDS).
    // CDS at genomic 25_000_050. Exon 1 coding region: 25_000_050 - 25_000_299.
    // Pick a position deep in the CDS: CDS pos ~200 = genomic 25_000_249.
    // C>CGA: ref_len=1, alt_len=3, net=+2bp → frameshift.
    // net_insertion_nt=2, ceil(2/3)=1, floor(2/3)=0.
    // With ceil: scan_limit=1+1=2 (scans into boundary codon).
    // With floor: scan_limit=1+0=1 (only checks primary codon → no false stop).
    let variant = InputVariant::new(
        "21".into(),
        25_000_249,
        25_000_249,
        b"G".to_vec(),
        b"GTA".to_vec(), // 1bp → 3bp = net +2bp frameshift
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "frameshift should remain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion),
        "frameshift should not be downgraded to inframe_insertion"
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopRetainedVariant),
        "frameshift should not have stop_retained_variant"
    );
}

/// Frameshift near stop codon should not be downgraded to
/// inframe_insertion + stop_retained_variant.
///
/// The narrow-alt-codon guard on condition 3 of `is_stop_retained_insertion`:
/// for a frameshift insertion near the terminal stop codon the wide codon window
/// (via CdsSpanBounds) may find '*' at the same position in both ref and alt
/// peptides. But if the local codon at
/// `translation_start` in the alt CDS is disrupted (not a stop), the match
/// is a coincidental downstream shifted-frame stop, so condition 3 is
/// suppressed, keeping the variant as frameshift_variant.
///
/// The test transcript's CDS is 850bp (283 complete codons + 1 extra base).
/// A 1bp insertion at the last complete codon (codon 283, CDS pos 847-849)
/// creates a +1bp frameshift that disrupts the stop boundary.
#[test]
fn concordance_frameshift_near_stop_not_downgraded_to_stop_retained() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Place a 1bp frameshift insertion near the end of exon 3.
    // CDS byte 847 = codon 283 start, genomic = 25_004_000 + (847 - 551) = 25_004_296.
    // G→GA: ref_len=1, alt_len=2, net=+1bp → frameshift.
    let variant = InputVariant::new(
        "21".into(),
        25_004_296,
        25_004_296,
        b"G".to_vec(),
        b"GA".to_vec(), // net +1bp frameshift
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    // Should be frameshift, not stop_retained.
    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "frameshift near stop should remain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopRetainedVariant),
        "frameshift near stop should not have stop_retained_variant, got: {:?}",
        tc.consequences
    );
}

/// reverse-strand anchored net insertion should normalize
/// to the transcript-sense differing region before peptide containment.
///
/// Perl trims shared flanks in transcript space before checking
/// `protein_altering_variant` vs `inframe_insertion`. On reverse strand, a
/// genomic 2bp->5bp anchored indel can collapse to a pure 3bp insertion in
/// transcript sense; the raw anchored substitution window would read as
/// `protein_altering_variant`.
///
/// Expected (Perl): `inframe_insertion`.
#[test]
fn concordance_reverse_strand_internal_insertion_normalizes_before_peptide_check() {
    let tx = make_reverse_strand_transcript();
    let config = EffectsConfig::default();

    // Reverse-strand transcript, CDS positions 4-5 ("GC" in transcript sense).
    // Raw genomic allele is suffix-anchored: GC -> GCTTT.
    // In transcript sense this is GC -> AAAGC, which trims to a pure AAA insertion
    // before the retained GCT codon, so Perl emits inframe_insertion.
    let variant = InputVariant::new(
        "21".into(),
        25_000_245,
        25_000_246,
        b"GC".to_vec(),
        b"GCTTT".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::InframeInsertion),
        "reverse-strand anchored insertion should be inframe_insertion, \
         got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "should NOT remain protein_altering_variant after \
         transcript-sense normalization"
    );
}

/// downstream in-frame ATG must not resurrect
/// start_retained_variant after a start-codon indel has already caused start_lost.
///
/// Perl's `start_retained_variant` for indels does not co-emit once
/// `start_lost` is true; scanning the entire altered CDS for a later in-frame
/// `ATG` would produce `frameshift_variant,start_lost,start_retained_variant`.
#[test]
fn concordance_start_lost_does_not_coemit_start_retained() {
    let tx = make_transcript_with_downstream_inframe_atg();
    let config = EffectsConfig::default();

    // 2bp insertion at CDS position 1 (genomic 25_000_050): A -> AGG.
    // The altered CDS still contains an in-frame ATG later on, but the start
    // codon itself is destroyed, so Perl emits start_lost without
    // start_retained_variant.
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_050,
        b"A".to_vec(),
        b"AGG".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "should remain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StartLost),
        "start codon indel should emit start_lost, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StartRetainedVariant),
        "downstream ATG should not cause start_retained_variant, got: {:?}",
        tc.consequences
    );
}

/// Upstream-extending deletion: a deletion whose reference span extends upstream of
/// the transcript across the 5'UTR/CDS boundary must not overcall `start_lost`.
///
/// With `paired_position` `Upstream` and a length change, `is_start_lost` would
/// fire and `coding::replicate_ins_del_start_altered` cannot answer for an
/// upstream extension: it short-circuits to `None` without a 5'UTR sequence or
/// when `cdna_start == 0` (Perl's mapper sentinel for "before transcript"). Perl
/// never reaches `_ins_del_start_altered` in this shape because its call site
/// requires both `cds_start` and `cdna_start` defined, and emits only
/// `coding_sequence_variant`.
///
/// Expected (Perl): no `start_lost` for upstream-extending deletions across
/// the 5'UTR/CDS boundary.
#[test]
fn concordance_upstream_extending_deletion_no_start_lost_overcall() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Forward-strand transcript starts at 25_000_000, CDS at 25_000_050.
    // Construct a 60bp deletion whose ref span extends 5bp upstream of the
    // transcript (24_999_995) into 5bp of CDS (25_000_054). The other
    // endpoint (24_999_995) is upstream of the transcript -> paired_position
    // resolves to `Upstream`, length-changing fires, `is_start_lost` evaluates
    // true via the FivePrimeUtr/Upstream gate.
    let ref_allele = vec![b'A'; 60];
    let variant = InputVariant::new(
        "21".into(),
        24_999_995,
        25_000_054,
        ref_allele,
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::StartLost),
        "Upstream-extending deletion across the \
         5'UTR/CDS boundary must NOT emit start_lost (Perl emits only \
         coding_sequence_variant), got: {:?}",
        tc.consequences
    );
}

/// Upstream-extending deletion (companion): the upstream-extension shape also
/// requires `Consequence::FivePrimeUtrVariant` alongside the
/// `coding_sequence_variant` fallthrough, Perl's
/// `5_prime_UTR_variant,coding_sequence_variant` pair. With `is_start_lost`
/// suppressed the variant takes the generic `else` branch, and
/// `add_utr_for_overlapping_span` does not re-fire for this paired-position
/// shape, so `FivePrimeUtrVariant` is added when the upstream-extension
/// predicate is true and `paired_position` is `FivePrimeUtr`.
///
/// Expected (Perl): both `5_prime_UTR_variant` and `coding_sequence_variant`
/// for upstream-extending deletions; no `start_lost`.
#[test]
fn concordance_upstream_extending_deletion_emits_5utr_term() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Forward-strand transcript starts at 25_000_000, CDS at 25_000_050.
    // Construct a 60bp deletion whose ref span extends 5bp upstream of the
    // transcript (24_999_995) into 5bp of CDS (25_000_054). One endpoint
    // (24_999_995) is upstream; the other (25_000_054) lies inside CDS, so
    // the `paired_position` for the upstream endpoint resolves to `FivePrimeUtr`
    // (the cdna step from upstream into the transcript first hits the 5'UTR).
    let ref_allele = vec![b'A'; 60];
    let variant = InputVariant::new(
        "21".into(),
        24_999_995,
        25_000_054,
        ref_allele,
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FivePrimeUtrVariant),
        "Upstream-extending deletion across the \
         5'UTR/CDS boundary must emit 5_prime_UTR_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences
            .contains(&Consequence::CodingSequenceVariant),
        "Upstream-extending deletion across the \
         5'UTR/CDS boundary must emit coding_sequence_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StartLost),
        "Upstream-extending deletion \
         across the 5'UTR/CDS boundary must NOT emit start_lost, got: {:?}",
        tc.consequences
    );
}

// HGVS notation parity

/// Intronic duplication notation: Perl detects duplications in intronic
/// insertions.
///
/// Expected (Perl): HGVSc contains `dup` not `ins`.
#[test]
fn concordance_hgvsc_intronic_insertion_uses_dup_notation() {
    use crate::hgvs::generate_hgvsc;

    let tx = make_compact_test_transcript();
    let mut sequence = vec![b'A'; 7000];
    sequence[410] = b'G';
    sequence[411] = b'A';
    let sequence = String::from_utf8(sequence).expect("ASCII reference");
    let (_dir, fasta) = load_test_fasta("1", &sequence);

    // Insertion in intron 1 that right-shifts once into a duplicated G.
    // Compact transcript intron 1: 400-2099. Insert between 410 and 411.
    // Base 411 is G, so Perl-style HGVS should emit dup after 3' normalization.
    let variant = InputVariant::new(
        "1".into(),
        411,
        410, // insertion
        b"-".to_vec(),
        b"G".to_vec(),
    );

    let hgvsc = generate_hgvsc(&variant, &tx, Some(&fasta));
    assert!(
        hgvsc.is_some(),
        "Should generate HGVSc for intronic insertion"
    );
    let notation = hgvsc.unwrap();

    assert!(
        notation.contains("dup"),
        "intronic insertion duplicating preceding bases should use \
         'dup' notation, got: {notation}"
    );
    assert!(
        !notation.contains("ins"),
        "should not use 'ins' notation for duplication"
    );
}

/// Frameshift Xaa convention: Perl's notation for unknown amino acids in
/// frameshifts.
///
/// Expected (Perl): HGVSp uses Perl's Xaa convention.
#[test]
fn concordance_hgvsp_frameshift_uses_fs_notation() {
    use crate::hgvs::generate_hgvsp;

    let tx = make_test_transcript();

    // Frameshift insertion in CDS.
    // CDS at genomic 25_000_056 = codon 3 (Gly).
    let variant = InputVariant::new(
        "21".into(),
        25_000_057,
        25_000_056, // insertion
        b"-".to_vec(),
        b"T".to_vec(), // 1bp insertion = frameshift
    );

    let hgvsp = generate_hgvsp(&variant, &tx, None);
    assert!(hgvsp.is_some(), "Should generate HGVSp for frameshift");
    let notation = hgvsp.unwrap();

    // Perl uses "fs" notation for frameshifts. Check it matches Perl convention.
    assert!(
        notation.contains("fs"),
        "frameshift HGVSp should contain 'fs', got: {notation}"
    );
}

/// Ter position in frameshift notation: the distance to the downstream stop
/// codon.
///
/// Expected (Perl): `fsTerN` where N matches Perl's calculation.
#[test]
fn concordance_hgvsp_frameshift_includes_ter_position() {
    use crate::hgvs::generate_hgvsp;

    let tx = make_test_transcript();

    // Frameshift near the start of the CDS where a downstream stop is findable.
    // Insert at codon 4 (Lys, CDS offset 10-12, genomic 25_000_059-25_000_061).
    let variant = InputVariant::new(
        "21".into(),
        25_000_060,
        25_000_059, // insertion
        b"-".to_vec(),
        b"CC".to_vec(), // 2bp insertion = frameshift
    );

    let hgvsp = generate_hgvsp(&variant, &tx, None);
    assert!(hgvsp.is_some(), "Should generate HGVSp for frameshift");
    let notation = hgvsp.unwrap();

    // Check that the Ter position is present and numeric.
    assert!(
        notation.contains("Ter"),
        "frameshift should include Ter position, got: {notation}"
    );
}

/// Insertion rotation after 3' shift
///
/// Expected (Perl): alt allele sequence is rotated after 3' normalization.
#[test]
fn concordance_hgvsc_reverse_strand_insertion_rotates_after_shift() {
    use crate::hgvs::generate_hgvsc;

    let tx = make_compact_reverse_strand_transcript();
    let mut sequence = vec![b'A'; 7000];
    sequence[199] = b'T';
    sequence[198] = b'C';
    sequence[200] = b'C';
    let sequence = String::from_utf8(sequence).expect("ASCII reference");
    let (_dir, fasta) = load_test_fasta("1", &sequence);

    // Reverse-strand insertion that shifts one base to the left on the genome.
    // Raw genomic alt AT rotates to TA after transcript-3' normalization.
    let variant = InputVariant::new(
        "1".into(),
        201,
        200, // insertion
        b"-".to_vec(),
        b"AT".to_vec(),
    );

    let hgvsc = generate_hgvsc(&variant, &tx, Some(&fasta));
    assert!(
        hgvsc.is_some(),
        "should generate HGVSc for reverse-strand insertion"
    );
    let notation = hgvsc.unwrap();
    assert!(
        notation.contains("insTA"),
        "shifted reverse-strand insertion should rotate alt to TA, got: {notation}"
    );
    assert!(
        !notation.contains("insAT"),
        "should not keep the unrotated insertion sequence after 3' shift"
    );
    assert!(
        !notation.contains("dup"),
        "this fixture should stay in insertion notation, got: {notation}"
    );
}

#[test]
fn concordance_hgvsc_forward_strand_insertion_rotates_after_shift() {
    use crate::hgvs::generate_hgvsc;

    let tx = make_compact_test_transcript();
    let mut sequence = vec![b'A'; 7000];
    sequence[200] = b'A';
    sequence[201] = b'C';
    let sequence = String::from_utf8(sequence).expect("ASCII reference");
    let (_dir, fasta) = load_test_fasta("1", &sequence);

    let variant = InputVariant::new("1".into(), 201, 200, b"-".to_vec(), b"AT".to_vec());

    let hgvsc = generate_hgvsc(&variant, &tx, Some(&fasta));
    assert!(
        hgvsc.is_some(),
        "Should generate HGVSc for forward-strand insertion"
    );
    let notation = hgvsc.unwrap();

    assert!(
        notation.contains("insTA"),
        "Forward-strand 3' shift should rotate insertion sequence to TA, got: {notation}"
    );
}

/// Protein position.
///
/// Expected (Perl): protein position matches exactly.
#[test]
fn concordance_hgvsp_protein_position_matches_codon_number() {
    use crate::hgvs::generate_hgvsp;

    let tx = make_test_transcript();

    // Inframe deletion of 3bp (1 codon) at codon 3.
    // Codon 3 = CDS offset 7-9 = genomic 25_000_056-25_000_058.
    let variant = InputVariant::new(
        "21".into(),
        25_000_056,
        25_000_058,
        b"GGA".to_vec(),
        b"-".to_vec(),
    );

    let hgvsp = generate_hgvsp(&variant, &tx, None);

    // Codon 3 = amino acid 3 = Gly, and the residue after it (Lys) differs, so
    // the deletion does not move: Perl says p.Gly3del.
    assert_eq!(hgvsp.as_deref(), Some("ENSP00000000001.1:p.Gly3del"));
}

/// Inversion notation: Perl detects inversions (CT to AG) and uses `inv`.
///
/// Expected (Perl): HGVSc contains `inv` not `delins`.
#[test]
fn concordance_hgvsc_inversion_uses_inv_notation() {
    use crate::hgvs::generate_hgvsc;

    let tx = make_test_transcript();

    // 2bp inversion in CDS: CT→AG at a known position.
    // CDS offset 5-6 = genomic 25_000_054-25_000_055.
    // The test transcript has GCT at codon 2 (pos 4-6), so positions 5-6 = CT.
    // Reverse complement of CT = AG.
    let variant = InputVariant::new(
        "21".into(),
        25_000_054,
        25_000_055,
        b"CT".to_vec(),
        b"AG".to_vec(), // reverse complement = inversion
    );

    let hgvsc = generate_hgvsc(&variant, &tx, None);
    assert!(hgvsc.is_some(), "Should generate HGVSc for inversion");
    let notation = hgvsc.unwrap();

    assert!(
        notation.contains("inv"),
        "2bp inversion (CT→AG) should use 'inv' notation, got: {notation}"
    );
    assert!(
        !notation.contains("delins"),
        "should not use 'delins' for a true inversion"
    );
}

// HGVSp parity with Perl's `hgvs_protein` (`coding::perl_hgvs_protein`). Every
// expected string is what `TranscriptVariationAllele::hgvs_protein` prints for the
// same fixture; the derivations cite its helpers. `make_test_transcript()` has no
// cached UTR strings and no FASTA, so the alternate CDS carries `N`-padded UTRs.

/// `hgvsp_stop(utr)`: the terminal-stop fixture with a 3' UTR of `utr`.
fn make_transcript_with_terminal_stop_and_utr(utr: &str) -> Transcript {
    let mut tx = make_transcript_with_terminal_stop();
    if let Some(ref mut vefc) = tx.vefc {
        vefc.three_prime_utr = Some(utr.to_string());
    }
    tx
}

fn hgvsp_on(
    tx: &Transcript,
    start: u64,
    end: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
) -> Option<String> {
    let variant = InputVariant::new(
        "21".into(),
        start,
        end,
        ref_allele.to_vec(),
        alt_allele.to_vec(),
    );
    crate::hgvs::generate_hgvsp(&variant, tx, None)
}

/// H07. In-frame insertion between codons 3 and 4 whose reference peptide window
/// is empty: `_get_hgvs_protein_type` reads `ins`, `_get_surrounding_peptides`
/// names the flanking residues.
#[test]
fn concordance_hgvsp_h07_insertion_between_codons_names_flanking_residues() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_059, 25_000_058, b"-", b"AGA").as_deref(),
        Some("ENSP00000000001.1:p.Gly3_Lys4insArg")
    );
}

/// H08. The inserted residue equals the residue before the site:
/// `_check_for_peptide_duplication` turns the insertion into a duplication.
#[test]
fn concordance_hgvsp_h08_insertion_equal_to_preceding_residue_is_dup() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_059, 25_000_058, b"-", b"GGA").as_deref(),
        Some("ENSP00000000001.1:p.Gly3dup")
    );
}

/// H09. Two inserted residues equal to residues 2 and 3: a ranged duplication.
#[test]
fn concordance_hgvsp_h09_two_residue_duplication_is_ranged() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_059, 25_000_058, b"-", b"GCTGGA").as_deref(),
        Some("ENSP00000000001.1:p.Ala2_Gly3dup")
    );
}

/// H10. The inserted residue equals the residue after the site: `_shift_3prime`
/// moves the insertion past Lys4, and the duplication check then matches it.
#[test]
fn concordance_hgvsp_h10_insertion_shifts_past_matching_residue_then_dups() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_059, 25_000_058, b"-", b"AAA").as_deref(),
        Some("ENSP00000000001.1:p.Lys4dup")
    );
}

/// H11a. An inserted stop: `*` becomes `X`, then `Xaa`, then `Ter`.
#[test]
fn concordance_hgvsp_h11a_inserted_stop_is_ter() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_059, 25_000_058, b"-", b"TAA").as_deref(),
        Some("ENSP00000000001.1:p.Gly3_Lys4insTer")
    );
}

/// H11b. A residue then a stop: both are reported.
#[test]
fn concordance_hgvsp_h11b_residue_then_stop_reports_both() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_059, 25_000_058, b"-", b"GATTAA").as_deref(),
        Some("ENSP00000000001.1:p.Gly3_Lys4insAspTer")
    );
}

/// H11c. A stop then a residue: `s/Ter\w+/Ter/` drops what follows the stop.
#[test]
fn concordance_hgvsp_h11c_stop_then_residue_reports_stop_only() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_059, 25_000_058, b"-", b"TAAGAT").as_deref(),
        Some("ENSP00000000001.1:p.Gly3_Lys4insTer")
    );
}

/// H12a. Stop-loss substitution with the next stop in the 3' UTR: the first stop
/// of the alternate translation sits at residue 286, and
/// `_stop_loss_extra_AA` counts 286 - 1 - 282 = 3.
#[test]
fn concordance_hgvsp_h12a_stop_loss_substitution_counts_to_utr_stop() {
    let tx = make_transcript_with_terminal_stop_and_utr("GCAGCATAA");
    assert_eq!(
        hgvsp_on(&tx, 25_004_296, 25_004_296, b"T", b"C").as_deref(),
        Some("ENSP00000000001.1:p.Ter283GlnextTer3")
    );
}

/// H12b. Stop-loss substitution with no stop in the alternate translation: `?`.
#[test]
fn concordance_hgvsp_h12b_stop_loss_substitution_without_downstream_stop() {
    let tx = make_transcript_with_terminal_stop_and_utr("GCAGCAGCA");
    assert_eq!(
        hgvsp_on(&tx, 25_004_296, 25_004_296, b"T", b"C").as_deref(),
        Some("ENSP00000000001.1:p.Ter283GlnextTer?")
    );
}

/// H13a. In-frame deletion of the last residue and the stop: the reference
/// peptide is `A*`, the type `del`, and the stop-loss branch of
/// `_get_hgvs_protein_format` counts 284 - 1 - 282 = 1 residue to the UTR stop.
#[test]
fn concordance_hgvsp_h13a_deletion_of_last_residue_and_stop() {
    let tx = make_transcript_with_terminal_stop_and_utr("GCAGCATAA");
    assert_eq!(
        hgvsp_on(&tx, 25_004_293, 25_004_298, b"GCTTAA", b"-").as_deref(),
        Some("ENSP00000000001.1:p.Ala282_Ter283delextTer1")
    );
}

/// H13b. In-frame deletion of the stop codon alone: 285 - 1 - 282 = 2.
#[test]
fn concordance_hgvsp_h13b_deletion_of_stop_codon() {
    let tx = make_transcript_with_terminal_stop_and_utr("GCAGCATAA");
    assert_eq!(
        hgvsp_on(&tx, 25_004_296, 25_004_298, b"TAA", b"-").as_deref(),
        Some("ENSP00000000001.1:p.Ter283delextTer2")
    );
}

/// H13c. An insertion inside the stop codon that recreates a stop: `_clip_alleles`
/// returns `=` on the leading stops, `_get_hgvs_protein_type` then reads `X`
/// against `X*` as `delins`, and the alt `TerTer` collapses to `Ter`.
#[test]
fn concordance_hgvsp_h13c_insertion_in_stop_codon_recreating_stop() {
    let tx = make_transcript_with_terminal_stop_and_utr("GCAGCATAA");
    assert_eq!(
        hgvsp_on(&tx, 25_004_297, 25_004_296, b"-", b"AGT").as_deref(),
        Some("ENSP00000000001.1:p.Ter283delinsTer")
    );
}

/// H14. In-frame deletion of one codon whose following residue differs:
/// `_shift_3prime` does not move it.
#[test]
fn concordance_hgvsp_h14_single_codon_deletion_does_not_shift() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_056, 25_000_058, b"GGA", b"-").as_deref(),
        Some("ENSP00000000001.1:p.Gly3del")
    );
}

/// H15. Two codons replaced by two others: `delins` over the range.
#[test]
fn concordance_hgvsp_h15_two_codon_replacement_is_delins() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_053, 25_000_058, b"GCTGGA", b"AAAAAA").as_deref(),
        Some("ENSP00000000001.1:p.Ala2_Gly3delinsLysLys")
    );
}

/// H16. Deleting the start codon with a 5' UTR that ends in `ATG`: Perl's
/// `start_lost` predicate holds and it prints `Met1?`; this engine keeps
/// `start_retained_variant` (the CDS is a suffix of the edited sequence) and
/// drops `start_lost`, so the short-circuit is suppressed and the deletion is
/// described as one. Intended divergence.
#[test]
fn concordance_hgvsp_h16_start_codon_deletion_with_retained_start_is_del() {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        vefc.five_prime_utr = Some(format!("{}GCCACATG", "GCCACC".repeat(7)));
    }
    assert_eq!(
        hgvsp_on(&tx, 25_000_050, 25_000_052, b"ATG", b"-").as_deref(),
        Some("ENSP00000000001.1:p.Met1del")
    );
}

/// H17. An in-frame insertion of a second Met after the start codon: start
/// retained only, `_clip_alleles` leaves an insertion of `M` at 2, and the
/// duplication check matches the residue before it.
#[test]
fn concordance_hgvsp_h17_insertion_duplicating_the_start_methionine() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_052, 25_000_051, b"-", b"GAT").as_deref(),
        Some("ENSP00000000001.1:p.Met1dup")
    );
}

/// H18. Deleting one codon of a long Ala run: `_shift_3prime` walks the deletion
/// to the run's last residue.
#[test]
fn concordance_hgvsp_h18_deletion_in_a_repeat_shifts_to_its_end() {
    let tx = make_test_transcript();
    assert_eq!(
        hgvsp_on(&tx, 25_000_068, 25_000_070, b"GCT", b"-").as_deref(),
        Some("ENSP00000000001.1:p.Ala283del")
    );
}

/// No spurious UTR term on an SNV at the CDS/UTR boundary: the per-endpoint
/// mapping already classifies a single position as coding or UTR, so
/// `add_utr_for_overlapping_span()` must not fire for it.
///
/// Expected (Perl): `missense_variant` only, no spurious `3_prime_UTR_variant`.
#[test]
fn concordance_snv_at_cds_utr_boundary_no_spurious_utr_term() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // SNV at the last coding position before the 3' UTR boundary.
    // coding_region_end = 25_004_299 (last CDS base on forward strand).
    // A SNV here is fully within CDS; it must not get 3_prime_UTR_variant.
    let variant = InputVariant::new(
        "21".into(),
        25_004_299,
        25_004_299,
        b"G".to_vec(),
        b"A".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
        "SNV at last CDS position should NOT contain 3_prime_UTR_variant, \
         got: {:?}",
        tc.consequences
    );

    // Also verify for a SNV at the first position of the 3' UTR itself.
    // coding_region_end + 1 = 25_004_300 is first 3' UTR position.
    // This should get 3_prime_UTR_variant but not any coding consequence.
    let utr_variant = InputVariant::new(
        "21".into(),
        25_004_300,
        25_004_300,
        b"G".to_vec(),
        b"A".to_vec(),
    );

    let utr_result = calculate_consequences(&utr_variant, &tx, &config);
    assert!(utr_result.is_some(), "UTR SNV should produce consequences");
    let utr_tc = utr_result.unwrap();

    assert!(
        utr_tc
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant),
        "SNV in 3' UTR should contain 3_prime_UTR_variant, got: {:?}",
        utr_tc.consequences
    );

    // Verify 5' UTR boundary: SNV at last position before CDS.
    // coding_region_start = 25_000_050. Position 25_000_049 = last 5' UTR base.
    let five_utr_variant = InputVariant::new(
        "21".into(),
        25_000_049,
        25_000_049,
        b"G".to_vec(),
        b"A".to_vec(),
    );

    let five_result = calculate_consequences(&five_utr_variant, &tx, &config);
    assert!(
        five_result.is_some(),
        "5' UTR SNV should produce consequences"
    );
    let five_tc = five_result.unwrap();

    assert!(
        five_tc
            .consequences
            .contains(&Consequence::FivePrimeUtrVariant),
        "SNV in 5' UTR should contain 5_prime_UTR_variant, got: {:?}",
        five_tc.consequences
    );

    // SNV at the first CDS position should not get 5_prime_UTR_variant.
    let cds_start_variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_050,
        b"A".to_vec(),
        b"T".to_vec(),
    );

    let cds_start_result = calculate_consequences(&cds_start_variant, &tx, &config);
    assert!(
        cds_start_result.is_some(),
        "CDS start SNV should produce consequences"
    );
    let cds_start_tc = cds_start_result.unwrap();

    assert!(
        !cds_start_tc
            .consequences
            .contains(&Consequence::FivePrimeUtrVariant),
        "SNV at first CDS position should NOT contain 5_prime_UTR_variant, \
         got: {:?}",
        cds_start_tc.consequences
    );
}

/// Complex indel at exon boundary: the reference span (lo..hi) of a 3bp ref, 10bp
/// alt at the last exonic bases is entirely exonic, but the insertion's **effect**
/// extends into the splice donor/acceptor/5th_base region, and Perl's
/// `_overlapped_introns` considers the full effect span.
///
/// Expected (Perl): includes splice_donor_variant, splice_donor_5th_base_variant,
/// and intron_variant for a complex indel at the last exonic bases extending into
/// the splice donor region.
#[test]
fn concordance_complex_indel_at_exon_end_reaches_splice_donor() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Complex indel at the last 3 bases of exon 1.
    // Exon 1 ends at 25_000_299. Intron 1 starts at 25_000_300.
    // Splice donor site: intron positions 25_000_300-25_000_301 (first 2bp).
    // Splice donor 5th base: intron position 25_000_304.
    //
    // Ref: 3bp at 25_000_297-25_000_299 (last 3 exonic bases)
    // Alt: 10bp (net +7bp extending into the splice donor region)
    //
    // The effective_hi = 25_000_299 + 7 = 25_000_306, which overlaps intron 1
    // (starts at 25_000_300). This triggers splice site evaluation.
    let variant = InputVariant::new(
        "21".into(),
        25_000_297,
        25_000_299,
        b"TCT".to_vec(),
        b"TCCCTGAAAA".to_vec(), // 10bp alt, net +7bp
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::SpliceDonorVariant),
        "complex indel at exon boundary should include splice_donor_variant, \
         got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant),
        "complex indel at exon boundary should include \
         splice_donor_5th_base_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::IntronVariant),
        "complex indel at exon boundary should include intron_variant, \
         got: {:?}",
        tc.consequences
    );
}

/// MNV stop supplement: stop detection is limited to the codons the MNV
/// directly affects (`compute_codon_window_peptide_alleles`), Perl's local scope;
/// a full-CDS translation with flank trim would find stops beyond the window.
///
/// This test places a 2bp MNV at a codon boundary (codons 6-7) that produces missense changes
/// but no stop. The MNV stop supplement should not add `stop_gained`.
#[test]
fn concordance_mnv_stop_check_reads_codon_window_only() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Standard test transcript CDS layout (forward strand):
    //   Codon 5: TTC (Phe) at CDS pos 13-15, genomic 25_000_062-25_000_064
    //   Codon 6: GAT (Asp) at CDS pos 16-18, genomic 25_000_065-25_000_067
    //   Codon 7: GCT (Ala) at CDS pos 19-21, genomic 25_000_068-25_000_070
    //
    // 2bp MNV at last base of codon 6 + first base of codon 7:
    //   genomic 25_000_067-25_000_068, ref=TG, alt=CA
    //   Codon 6: GA[T->C] = GAC (Asp, D): synonymous
    //   Codon 7: [G->A]CT = ACT (Thr, T): missense
    //   Neither codon becomes a stop.
    let variant = InputVariant::new(
        "21".into(),
        25_000_067,
        25_000_068,
        b"TG".to_vec(),
        b"CA".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences for 2bp MNV");
    let tc = result.unwrap();

    // The codon-window approach should not find stop_gained since neither
    // affected codon (GAC, ACT) is a stop codon.
    assert!(
        !tc.consequences.contains(&Consequence::StopGained),
        "2bp MNV at codon boundary (GAT/GCT -> GAC/ACT) should not \
         produce stop_gained via codon-window analysis, got: {:?}",
        tc.consequences
    );

    // The variant should produce a missense or protein_altering consequence.
    assert!(
        tc.consequences.contains(&Consequence::MissenseVariant)
            || tc
                .consequences
                .contains(&Consequence::ProteinAlteringVariant),
        "2bp MNV should produce missense_variant or protein_altering_variant, \
         got: {:?}",
        tc.consequences
    );
}

// mature_miRNA overcall: miRNA rows should not include non_coding_transcript_exon_variant

/// Perl `within_mature_miRNA` tests the raw variant span against each mature
/// segment's genomic range with the non-normalising `overlap()`, so an
/// insertion (start = end + 1) between the last mature base and the next base
/// has `start > c.end` and lies outside; `non_coding_exon_variant` then holds.
#[test]
fn concordance_insertion_after_last_mature_base_is_exon_variant() {
    let tx = make_mirna_transcript();
    let config = EffectsConfig::default();

    // Mature range cDNA 50-100 is genomic 25_000_049-25_000_099 on this forward
    // transcript; the insertion anchors at 25_000_100 with end 25_000_099.
    let variant = InputVariant::new(
        "21".into(),
        25_000_100,
        25_000_099,
        b"-".to_vec(),
        b"A".to_vec(),
    );

    let tc = calculate_consequences(&variant, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["non_coding_transcript_exon_variant"]);
}

/// The whole span decides `within_mature_miRNA`: a deletion whose endpoints both
/// lie outside the mature range but cover it is inside, and Perl's
/// `non_coding_exon_variant` returns 0 whenever `within_mature_miRNA` holds.
#[test]
fn concordance_deletion_covering_mature_mirna_is_mature_mirna_variant() {
    let tx = make_mirna_transcript();
    let config = EffectsConfig::default();

    let variant = InputVariant::new(
        "21".into(),
        25_000_040,
        25_000_110,
        b"A".repeat(71),
        b"-".to_vec(),
    );

    let tc = calculate_consequences(&variant, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["mature_miRNA_variant"]);
}

/// Reverse strand: the mature range's 3' end (cDNA 100) is its lowest genomic
/// base, so the insertion immediately below it in genomic coordinates is the
/// one whose raw `end` falls short of `c.start`.
#[test]
fn concordance_reverse_strand_insertion_below_mature_range_is_exon_variant() {
    let mut tx = make_reverse_strand_transcript();
    tx.biotype = "miRNA".into();
    tx.translation = None;
    tx.protein_id = None;
    tx.cdna_coding_start = None;
    tx.cdna_coding_end = None;
    tx.coding_region_start = None;
    tx.coding_region_end = None;
    tx.translation_start = None;
    tx.translation_end = None;
    if let Some(ref mut vefc) = tx.vefc {
        vefc.translateable_seq = None;
        vefc.peptide = None;
        if let Some(ref mut mapper) = vefc.mapper {
            mapper.cdna_coding_start = 0;
            mapper.cdna_coding_end = 0;
        }
    }
    tx.attributes.push(Attribute {
        code: "miRNA".to_string(),
        value: "50-100".to_string(),
    });
    let config = EffectsConfig::default();

    // With ori -1 on the first mapper pair, cDNA 1 is genomic 25_000_299, so the
    // mature range cDNA 50-100 is genomic 25_000_200-25_000_250.
    let outside = InputVariant::new(
        "21".into(),
        25_000_200,
        25_000_199,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&outside, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["non_coding_transcript_exon_variant"]);

    let covering = InputVariant::new(
        "21".into(),
        25_000_190,
        25_000_260,
        b"A".repeat(71),
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&covering, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["mature_miRNA_variant"]);
}

/// Perl keeps a transcript for a variant only when
/// `overlap(vf.start, vf.end, tr.start - 5000, tr.end + 5000)` holds with the
/// non-normalising `overlap()`, whose left-hand test reads `vf.end >= tr.start - 5000`.
/// An insertion anchored exactly 5000 bases before a forward-strand transcript
/// has `end = tr.start - 5001` and gets no row; one base closer is upstream.
#[test]
fn concordance_insertion_at_upstream_window_edge_forward() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    let outside = InputVariant::new(
        "21".into(),
        24_995_000,
        24_994_999,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    assert!(
        calculate_consequences(&outside, &tx, &config).is_none(),
        "an insertion anchored 5000 bases before the transcript start is outside Perl's window"
    );

    let inside = InputVariant::new(
        "21".into(),
        24_995_001,
        24_995_000,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    let tc = calculate_consequences(&inside, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["upstream_gene_variant"]);
    assert_eq!(tc.distance, Some(4999));
}

/// Downstream side of a forward-strand transcript: the right-hand test reads
/// `vf.start <= tr.end + 5000`, so the anchor itself decides.
#[test]
fn concordance_insertion_at_downstream_window_edge_forward() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    let inside = InputVariant::new(
        "21".into(),
        25_011_000,
        25_010_999,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    let tc = calculate_consequences(&inside, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["downstream_gene_variant"]);
    assert_eq!(tc.distance, Some(4999));

    let outside = InputVariant::new(
        "21".into(),
        25_011_001,
        25_011_000,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    assert!(
        calculate_consequences(&outside, &tx, &config).is_none(),
        "an insertion anchored 5001 bases after the transcript end is outside Perl's window"
    );
}

/// Reverse strand, genomic-left side (downstream in transcript terms): the
/// window edge still reads off `vf.end`.
#[test]
fn concordance_insertion_at_downstream_window_edge_reverse() {
    let tx = make_reverse_strand_transcript();
    let config = EffectsConfig::default();

    let outside = InputVariant::new(
        "21".into(),
        24_995_000,
        24_994_999,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    assert!(
        calculate_consequences(&outside, &tx, &config).is_none(),
        "an insertion anchored 5000 bases below a reverse transcript's start is outside Perl's window"
    );

    let inside = InputVariant::new(
        "21".into(),
        24_995_001,
        24_995_000,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    let tc = calculate_consequences(&inside, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["downstream_gene_variant"]);
    assert_eq!(tc.distance, Some(4999));
}

/// Reverse strand, genomic-right side (upstream in transcript terms).
#[test]
fn concordance_insertion_at_upstream_window_edge_reverse() {
    let tx = make_reverse_strand_transcript();
    let config = EffectsConfig::default();

    let inside = InputVariant::new(
        "21".into(),
        25_011_000,
        25_010_999,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    let tc = calculate_consequences(&inside, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["upstream_gene_variant"]);
    assert_eq!(tc.distance, Some(4999));

    let outside = InputVariant::new(
        "21".into(),
        25_011_001,
        25_011_000,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    assert!(
        calculate_consequences(&outside, &tx, &config).is_none(),
        "an insertion anchored 5001 bases above a reverse transcript's end is outside Perl's window"
    );
}

/// Perl's `non_coding_exon_variant` re-checks true exon overlap after the
/// frameshift-intron stretch admits the position, so an SNV inside a frameshift
/// intron of a non-coding transcript carries only `non_coding_transcript_variant`:
/// no exon term, and no intron or splice term because `_intron_effects` skips
/// the frameshift intron.
#[test]
fn concordance_snv_in_frameshift_intron_of_non_coding_transcript() {
    let mut tx = make_non_coding_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(intron) = vefc.introns.get_mut(0) {
            intron.start = 25_000_300;
            intron.end = 25_000_303;
        }
    }
    if let Some(intron) = tx.introns.get_mut(0) {
        intron.start = 25_000_300;
        intron.end = 25_000_303;
    }
    let config = EffectsConfig::default();

    let variant = InputVariant::new(
        "21".into(),
        25_000_301,
        25_000_301,
        b"A".to_vec(),
        b"T".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["non_coding_transcript_variant"]);
}

/// A deleting allele spanning the whole transcript is still `transcript_ablation`
/// alone, including a delins whose alt is merely shorter than its ref.
#[test]
fn concordance_shortening_delins_spanning_transcript_is_ablation() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    let variant = InputVariant::new(
        "21".into(),
        24_999_990,
        25_006_010,
        b"A".repeat(6021),
        b"CCC".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).unwrap();
    crate::test_helpers::assert_consequence_set_eq(&tc, &["transcript_ablation"]);
}

/// stop_retained from a boundary codon: condition 1 of `is_stop_retained_insertion`
/// scans ref AAs + inserted AAs only; a `+ 1` would include a boundary codon AA
/// that Perl's narrower codon window never sees.
///
/// This test verifies the scan limit formula by constructing an insertion where:
/// - ref_pep is 1 AA ("G") matching alt_pep[0]
/// - The inserted sequence contains a stop codon (TAA) within the insertion window
/// - The function correctly detects the stop (positive case)
///
/// The companion test below pins the scan-window width for a 3 bp insertion; no
/// fixture here places a stop codon beyond the window.
#[test]
fn concordance_stop_retained_scan_limit_covers_ref_and_inserted_amino_acids() {
    use crate::consequences::is_stop_retained_insertion;
    use crate::mapper::CdsSpanBounds;

    let tx = make_test_transcript();

    // Positive case: 6bp pure insertion at CDS position 8 (inside codon 3: GGA = G).
    // Inserted "GATAAT" -> alt codon window contains a stop within the insertion:
    //   alt_window = "G" + "GATAAT" + "GA" = "GGATAAT GA"
    //   Codon 1: GGA (G): matches ref_pep[0]
    //   Codon 2: TAA (*): stop within insertion window
    //   Codon 3: TGA (*): the boundary codon, also inside the scan
    //
    // insertion_aa_count = 6/3 = 2, scan_limit = 1 + 2 = 3 (covers all 3 codons).
    // The stop at index 1 is well within the scan limit -> should return true.
    let variant_positive = InputVariant::new(
        "21".into(),
        25_000_008,
        25_000_007,
        b"-".to_vec(),
        b"GATAAT".to_vec(), // 6bp: creates TAA stop in the inserted sequence
    );

    let bounds = CdsSpanBounds {
        cds_start: 8,
        cds_end: 7,
        translation_start: 3,
        translation_end: 3,
    };

    let result_positive = is_stop_retained_insertion(&variant_positive, &tx, &bounds, true);
    assert!(
        result_positive,
        "stop within insertion window (TAA at alt_pep index 1) \
         should be detected by stop_retained"
    );

    // Boundary case: 6bp pure insertion where the inserted AAs have no stop,
    // but the boundary codon (formed from end of insertion + downstream CDS) is a stop.
    // Inserted "GAGCTT" -> alt codon window:
    //   alt_window = "G" + "GAGCTT" + "GA" = "GGAGCTT GA"
    //   Codon 1: GGA (G): matches ref_pep[0]
    //   Codon 2: GCT (A): no stop
    //   Codon 3: TGA (*): boundary stop
    //
    // insertion_aa_count = 2, scan_limit = 1 + 2 = 3.
    // alt_pep = [G, A, *]. The '*' is at index 2, which is < scan_limit 3, so it is
    // scanned and the call returns true: the scan covers ref_pep.len() +
    // insertion_aa_count positions (0, 1, 2). A stop at index 3, one position
    // beyond the insertion contribution, would lie outside it.
    let variant_boundary = InputVariant::new(
        "21".into(),
        25_000_008,
        25_000_007,
        b"-".to_vec(),
        b"GAGCTT".to_vec(), // 6bp: boundary codon TGA is at index 2 (within scan)
    );

    let result_boundary = is_stop_retained_insertion(&variant_boundary, &tx, &bounds, true);
    assert!(
        result_boundary,
        "stop at boundary codon (TGA at alt_pep index 2) is still \
         within the scan window (scan_limit = 3, covers indices 0-2)"
    );
}

/// The scan window of a 3bp insertion is two positions.
///
/// For a 3bp insertion (insertion_aa_count = 1) the scan_limit is 1 + 1 = 2,
/// covering only alt_pep[0..2] (the ref AA + 1 inserted AA); a stop codon at
/// index 2 or later does not trigger stop_retained. In this fixture the boundary
/// codon's stop lands at index 1, inside the window, so the verdict is true.
#[test]
fn concordance_stop_retained_three_bp_insertion_scan_window_is_two_positions() {
    use crate::consequences::is_stop_retained_insertion;
    use crate::mapper::CdsSpanBounds;

    let tx = make_test_transcript();

    // 3bp insertion at CDS position 8 (codon 3: GGA = G).
    // Inserted "GCT" (Ala): no stop in the insertion itself.
    //   alt_window = "G" + "GCT" + "GA" = "GGCTGA"
    //   Codon 1: GGC. Both GGC and GGA encode Glycine, so this matches
    //   ref_pep[0]: as amino acids they are the same byte 'G'.
    //   ref_pep = [G], alt_pep = [G, ...].
    //   Codon 2: TGA (*): a stop, but from the boundary, not the insertion.
    //
    // insertion_aa_count = 1, scan_limit = 1 + 1 = 2.
    // alt_pep = [G, *]. The '*' is at index 1, which is < scan_limit 2, so
    // condition 1 fires: a 3bp insertion that spells a stop in the immediately
    // following codon has a stop inside the scan window.
    let variant = InputVariant::new(
        "21".into(),
        25_000_008,
        25_000_007,
        b"-".to_vec(),
        b"GCT".to_vec(), // 3bp: no stop in insertion, but boundary forms TGA
    );

    let bounds = CdsSpanBounds {
        cds_start: 8,
        cds_end: 7,
        translation_start: 3,
        translation_end: 3,
    };

    // Verify via compute_codon_window_peptide_alleles what the actual peptides are.
    use crate::coding::compute_codon_window_peptide_alleles;
    let peps = compute_codon_window_peptide_alleles(&variant, &tx, &bounds);
    assert!(peps.is_some(), "peptide alleles should resolve");
    let (ref_pep, alt_pep) = peps.unwrap();
    assert_eq!(
        ref_pep.len(),
        1,
        "ref_pep should be single AA for mid-codon insertion"
    );

    // Only positions 0..scan_limit are scanned; any stop within that range is
    // detected.
    let result = is_stop_retained_insertion(&variant, &tx, &bounds, true);

    // If alt_pep[0] matches ref_pep[0] and there's a stop within the scan window,
    // condition 1 should fire. If alt_pep[0] doesn't match, condition 1 is skipped.
    if ref_pep[0] == alt_pep[0] {
        let scan_limit = ref_pep.len() + 1; // 1 ref AA + 1 inserted AA (3bp/3)
        let scan_end = scan_limit.min(alt_pep.len());
        let has_stop_in_window = alt_pep[..scan_end].contains(&b'*');
        assert_eq!(
            result, has_stop_in_window,
            "stop_retained should match scan window. \
             ref_pep={:?}, alt_pep={:?}, scan_limit={}, scan_end={}",
            ref_pep, alt_pep, scan_limit, scan_end
        );
    }
}

/// start_retained_variant should only fire when ATG appears
/// at a codon-aligned position in the alt CDS (indices 0, 3, 6, ...).
///
/// Perl only checks codon boundaries (reading frame aligned to CDS position 0); a
/// `windows(3).any(|w| w == b"ATG")` scan would also match out-of-frame positions.
///
/// Test: 1bp insertion at start codon that shifts the ATG out of frame.
/// CDS: ATG GCT ... → Insert 1bp "C" before ATG → C + ATG GCT ...
/// Alt CDS = "CAT GGC TGG A..." The original ATG is now split across
/// codons (positions 1-3), not at a codon boundary (position 0).
/// With the repetitive GCT fill, no codon-aligned ATG exists anywhere
/// in the alt CDS.
/// Perl: start_lost only (no start_retained).
#[test]
fn concordance_start_retained_must_be_codon_aligned() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Frameshift insertion of 1bp "C" at CDS position 1 (start codon).
    // CDS position 1 = genomic 25_000_050 (cdna_coding_start=51, exon1 starts at 25_000_000).
    // Pure insertion: ref="-", start=25_000_051, end=25_000_050 (VEP convention).
    // Alt CDS becomes: C + ATG GCT GGA ... = "CATGGCTGGA...".
    // Codon-aligned positions: CAT(0), GGC(3), TGG(6), ...: no ATG at any boundary.
    let variant = InputVariant::new(
        "21".into(),
        25_000_051, // start (after insertion point)
        25_000_050, // end (before insertion point): insertion before CDS pos 1
        b"-".to_vec(),
        b"C".to_vec(), // 1bp = frameshift
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "1bp insertion at start codon should be frameshift, got: {:?}",
        tc.consequences
    );
    // The original ATG is at position 1 in the alt CDS (out of frame), so no
    // codon-boundary ATG exists and no start_retained; a windows(3) scan would
    // match position 1.
    assert!(
        !tc.consequences.contains(&Consequence::StartRetainedVariant),
        "start_retained_variant should NOT fire when no ATG exists at a \
         codon-aligned position in the alt CDS, got: {:?}",
        tc.consequences
    );
}

/// Positive test: start_retained should fire when an indel overlapping the start
/// codon preserves `ATG` as the first codon of the altered CDS.
/// 3bp inframe insertion at the start codon: ATG -> ATG GCA (net +3bp).
/// Alt CDS still starts with Met, so Perl emits `start_retained_variant` rather
/// than `start_lost`.
#[test]
fn concordance_start_retained_fires_for_in_frame_atg() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // 3bp inframe insertion that keeps ATG as codon 1.
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_052,
        b"ATG".to_vec(),
        b"ATGGCA".to_vec(), // ATG -> ATGGCA = net +3bp, start codon preserved
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::StartRetainedVariant),
        "start_retained should fire when ATG remains as codon 1 in the alt CDS, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StartLost),
        "preserved start codon should not co-emit start_lost, got: {:?}",
        tc.consequences
    );
}

/// inframe_deletion far from stop codon should not gain
/// stop_gained from the transcript's terminal stop leaking through full-CDS
/// translation.
///
/// The full CDS always contains the terminal stop, so trimming common flanks
/// could expose it asymmetrically; codon-window peptides scope to the affected
/// codons only.
///
/// Test: 3bp inframe deletion at codon 4 (well within the CDS, far from stop).
/// CDS codon 4 = "AAA" (Lys). Delete 3bp → removes one codon.
/// CDS position 10-12. Genomic: 25_000_050 + 9 = 25_000_059 to 25_000_061.
/// This should not produce stop_gained.
#[test]
fn concordance_inframe_deletion_far_from_stop_no_stop_gained() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // 3bp inframe deletion at codon 4 (AAA=Lys, CDS positions 10-12).
    // Genomic: CDS starts at 25_000_050, so positions 10-12 = 25_000_059-25_000_061.
    let variant = InputVariant::new(
        "21".into(),
        25_000_059,
        25_000_061,
        b"AAA".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::InframeDeletion),
        "3bp deletion should be inframe_deletion, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopGained),
        "inframe_deletion far from stop should NOT gain stop_gained \
         (terminal stop leak from full-CDS), got: {:?}",
        tc.consequences
    );
}

/// Perl's `VariationEffect::protein_altering_variant` predicate
/// returns 0 when alt_pep starts with '*'. So when a stop codon appears at
/// position 0 of the alt peptide, Perl emits stop_gained alone (without
/// protein_altering_variant).
///
/// Test: replace "GGA" (Gly, codon 3) with "TAAGCT" (TAA=stop + GCT=Ala).
/// Alt peptide = "*A". Since alt_pep starts with '*', Perl suppresses
/// protein_altering and only emits stop_gained.
#[test]
fn concordance_leading_stop_delins_is_stop_gained_without_protein_altering() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Replace 3bp "GGA" with 6bp "TAAGCT" at codon 3 (CDS positions 7-9).
    // Genomic: 25_000_050 + 6 = 25_000_056 to 25_000_058.
    // Net +3bp inframe. Creates stop at alt_pep position 0.
    let variant = InputVariant::new(
        "21".into(),
        25_000_056,
        25_000_058,
        b"GGA".to_vec(),
        b"TAAGCT".to_vec(), // GGA → TAAGCT = net +3bp, TAA stop at position 0
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    // Perl parity: alt_pep starts with '*' → protein_altering suppressed,
    // only stop_gained emitted.
    assert!(
        tc.consequences.contains(&Consequence::StopGained),
        "alt peptide starting with stop → stop_gained, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "Perl suppresses protein_altering when alt_pep starts with '*', got: {:?}",
        tc.consequences
    );
}

/// Negative test: protein_altering_variant without a stop should not emit
/// stop_gained.
#[test]
fn concordance_protein_altering_without_stop_no_stop_gained() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Replace 3bp "GGA" with 6bp "GCTGCT" at codon 3 (no stop introduced).
    // GGA(Gly) → GCT(Ala) GCT(Ala). Net +3bp inframe.
    // Alt peptide = "AA", ref peptide = "G". "AA" doesn't contain "G" → protein_altering.
    let variant = InputVariant::new(
        "21".into(),
        25_000_056,
        25_000_058,
        b"GGA".to_vec(),
        b"GCTGCT".to_vec(), // GGA → GCTGCT = net +3bp, no stop
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    // protein_altering or inframe_insertion according to peptide containment;
    // either way, no stop.
    assert!(
        !tc.consequences.contains(&Consequence::StopGained),
        "protein_altering without stop in alt peptide should NOT emit \
         stop_gained, got: {:?}",
        tc.consequences
    );
}

/// transcript-sense pure insertions that create a stop
/// in the first altered codon should emit `stop_gained` without the extra
/// generic indel term.
#[test]
fn concordance_leading_stop_pure_insertion_is_stop_gained_only() {
    let tx = make_transcript_with_codon3_tct();
    let config = EffectsConfig::default();

    // Insert 9bp after the first base of codon 3 (TCT -> TAGTTGAAATCT).
    // The altered peptide starts with a gained stop, so Perl reports
    // stop_gained rather than protein_altering_variant/inframe_insertion.
    let variant = InputVariant::new(
        "21".into(),
        25_000_057,
        25_000_056,
        b"-".to_vec(),
        b"AGTTGAAAT".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::StopGained),
        "leading-stop insertion should emit stop_gained, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "leading-stop insertion should not keep protein_altering_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion),
        "leading-stop insertion should not keep inframe_insertion, got: {:?}",
        tc.consequences
    );
}

/// When `translateable_seq` is missing (None), `classify_inframe_indel_by_peptide`
/// returns None. Perl falls back to `coding_sequence_variant`, not
/// `inframe_insertion` or `inframe_deletion`.
#[test]
fn concordance_coding_seq_variant_fallback_for_insertion_without_translateable_seq() {
    let mut tx = make_test_transcript();
    // Remove translateable_seq to force classify_inframe_indel_by_peptide to return None
    if let Some(ref mut vefc) = tx.vefc {
        vefc.translateable_seq = None;
    }
    let config = EffectsConfig::default();

    // Inframe insertion: ref="-", alt="GCA" (3bp, in-frame) at CDS position
    // within exon2 (25_000_250 is mid-exon2 for test transcript).
    let variant = InputVariant::new(
        "21".into(),
        25_000_250,
        25_000_249, // insertion: end < start
        b"-".to_vec(),
        b"GCA".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences
            .contains(&Consequence::CodingSequenceVariant),
        "insertion without translateable_seq should fall back to \
         coding_sequence_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion),
        "should not assign inframe_insertion when peptide check is unavailable, got: {:?}",
        tc.consequences
    );
}

/// Same fallback for complex indels (ref.len() != alt.len(), neither is dash)
/// when translateable_seq is missing.
#[test]
fn concordance_coding_seq_variant_fallback_for_complex_indel() {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        vefc.translateable_seq = None;
    }
    let config = EffectsConfig::default();

    // Complex indel: ref="GGA" (3bp) → alt="GCAGCT" (6bp), net +3bp inframe
    let variant = InputVariant::new(
        "21".into(),
        25_000_056,
        25_000_058,
        b"GGA".to_vec(),
        b"GCAGCT".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences
            .contains(&Consequence::CodingSequenceVariant),
        "complex indel without translateable_seq should fall back to \
         coding_sequence_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion)
            && !tc.consequences.contains(&Consequence::InframeDeletion),
        "should not assign inframe_insertion/deletion when peptide check unavailable, got: {:?}",
        tc.consequences
    );
}

/// Frameshift-only insertion at stop-codon
/// boundary must emit `inframe_insertion,stop_retained_variant` when Perl's
/// `ref_eq_alt_sequence` condition 2 fires.
///
/// Perl classifies `inframe_insertion,stop_retained_variant` for a 2bp
/// insertion at the stop-codon boundary (translation_start past ref peptide
/// length, overflow=1 < 3, condition 2 holds), never `frameshift_variant`.
///
/// Reference shape: 11:5225600 A>AGT (HBB, reverse strand), a 2bp insertion just
/// before the terminal TAA stop codon on the transcript.
#[test]
fn concordance_frameshift_insertion_at_stop_boundary_is_stop_retained() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();

    // 2bp insertion between CDS pos 846 and 847 (between codon 282 GCT and
    // the terminal TAA stop at codon 283). CDS len = 849.
    // Exon 3 starts at genomic 25_004_000 with CDS pos 551, so CDS pos 846
    // -> genomic 25_004_000 + 295 = 25_004_295 and CDS pos 847 -> 25_004_296.
    // VEP-style insertion: start = 25_004_296, end = 25_004_295 (end<start).
    let variant = InputVariant::new(
        "21".into(),
        25_004_296,
        25_004_295,
        b"-".to_vec(),
        b"AG".to_vec(), // 2bp net insertion (frameshift)
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::InframeInsertion),
        "2bp insertion at stop boundary should classify as \
         inframe_insertion via ref_eq_alt_sequence condition 2, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StopRetainedVariant),
        "stop_retained_variant must fire when ref_eq_alt_sequence \
         condition 2 holds (overflow < 3), got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::FrameshiftVariant),
        "frameshift must be suppressed when stop_retained holds \
         (Perl's frameshift returns 0 when stop_retained is true), got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopLost),
        "2bp insertion before stop codon should NOT emit stop_lost \
         when the stop is preserved, got: {:?}",
        tc.consequences
    );
}

/// Negative case: a 4bp frameshift insertion in the
/// middle of the CDS (far from the terminal stop) must remain
/// `frameshift_variant` and not be reclassified as inframe_insertion +
/// stop_retained. This guards against condition 2 firing on far-from-stop
/// insertions and against condition 1 using a too-wide alt peptide scan.
///
/// Reference variant class: `10:18534227 G>GGTAA` (4bp insertion deep in a
/// 633aa CDS) where Perl reports `frameshift_variant`.
#[test]
fn concordance_mid_cds_frameshift_insertion_is_not_stop_retained() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();

    // 4bp insertion at CDS pos ~400 (well within 849bp CDS, far from stop).
    // CDS pos 400 is in exon 2 (CDS 251-550). Genomic for CDS pos 400 =
    // 25_002_000 + (400 - 251) = 25_002_149.
    let variant = InputVariant::new(
        "21".into(),
        25_002_150,
        25_002_149,
        b"-".to_vec(),
        b"GTAA".to_vec(), // 4bp frameshift, not at stop
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "mid-CDS: 4bp insertion far from stop must remain \
         frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion),
        "mid-CDS: frameshift far from stop should not be downgraded \
         to inframe_insertion, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StopRetainedVariant),
        "mid-CDS: frameshift far from stop should not gain \
         stop_retained_variant from condition 2 leakage, got: {:?}",
        tc.consequences
    );
}

/// Same fallback for ambiguous in-frame indels. Perl's peptide() returns undef
/// unless the allele sequence is unambiguous DNA, so consequence assignment
/// falls back to coding_sequence_variant rather than peptide-based terms.
#[test]
fn concordance_coding_seq_variant_fallback_for_ambiguous_indel() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Inframe complex indel with ambiguous DNA and an embedded stop motif: Perl
    // treats the peptide alleles as unavailable and falls back to
    // coding_sequence_variant, so the unambiguous-DNA gate must fire before any
    // peptide classification or stop_gained.
    let variant = InputVariant::new(
        "21".into(),
        25_000_056,
        25_000_058,
        b"GGA".to_vec(),
        b"GGATAANNN".to_vec(), // net +6bp, but contains ambiguous bases
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences
            .contains(&Consequence::CodingSequenceVariant),
        "ambiguous inframe indel should fall back to coding_sequence_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences
            .contains(&Consequence::ProteinAlteringVariant)
            && !tc.consequences.contains(&Consequence::InframeInsertion)
            && !tc.consequences.contains(&Consequence::InframeDeletion)
            && !tc.consequences.contains(&Consequence::StopGained),
        "ambiguous indel should not emit peptide-based consequences, got: {:?}",
        tc.consequences
    );
}

/// Frameshift + stop_gained across multiple ref codons: Perl's `stop_gained`
/// scans the entire alt peptide (`$alt_pep =~ /\*/`) with no position cap, so a
/// small nt delta spanning multiple ref codons (ref=3nt, alt=8nt, tr_start !=
/// tr_end) whose window yields `LLV*X` has its `*` counted. The codon window
/// peptide (`translate_codon_window_to_peptide`) is already clipped to Perl's
/// window width, so `frameshift_stop_gained_in_codon_window` scans the whole
/// peptide; a `1 + net_insertion_nt/3` cap would miss the stop.
///
/// Shape: chr9:6556200-6556202 TAA>AACCAGGA (reverse strand): ref codon window
/// `cATCaa` (HQ), alt `cTCCTGGTTTAaa` (LLV*X).
///
/// Expected (Perl): both `frameshift_variant` and `stop_gained`.
#[test]
fn concordance_frameshift_stop_gained_multicodon_span() {
    // Set up a transcript where a 3nt ref → 8nt alt (net +5nt) spans 2 ref codons
    // and produces an alt codon window containing a stop codon at index 2.
    //
    // Default make_test_transcript() fills codons 7+ with GCT (Ala). So:
    //   CDS positions (0-based): 18=G, 19=C, 20=T, 21=G, 22=C, 23=T, ...
    //   codon 7 (1-based 19-21): GCT (A)
    //   codon 8 (1-based 22-24): GCT (A)
    //
    // Variant: genomic positions mapping to CDS pos 20-22 (1-based), which is
    // 0-based bytes 19-21 = "CTG". Alt = "CCGCCTAA" (8nt, net +5nt).
    //
    // Codon window: tr_start=7, tr_end=8. codon_cds_start=19, codon_cds_end=24.
    //   window_start_idx=18, ref_window=cds[18..24]="GCTGCT" → "AA" (ref_pep)
    //
    // Alt CDS splice: cds[0..19] + "CCGCCTAA" + cds[22..]
    //   alt_cds[18..29] = "G" + "CCGCCTAA" + "CT" = "GCCGCCTAACT" (11 chars)
    //   translate: GCC=A, GCC=A, TAA=*, partial CT → append X (since pep[0] != *)
    //   → alt_pep = "AA*X" (4 chars)
    //
    // Scanning the full alt_pep "AA*X" finds `*` at index 2, so both
    // `frameshift_variant` and `stop_gained`; a cap of 1 + floor(5/3) = 2 would
    // scan "AA" only and miss it.
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // CDS pos 20 (1-based) = genomic 25_000_050 + 19 = 25_000_069.
    // CDS pos 22 (1-based) = genomic 25_000_071.
    let variant = InputVariant::new(
        "21".into(),
        25_000_069,
        25_000_071,
        b"CTG".to_vec(),
        b"CCGCCTAA".to_vec(), // 8nt = net +5nt frameshift, stop at codon 3 of alt window
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "should contain frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StopGained),
        "should contain stop_gained: alt peptide `AA*X` has `*` at position 2, which \
         a `1 + net_ins_nt/3 = 2` cap would miss; Perl scans the entire alt_pep. \
         Got: {:?}",
        tc.consequences
    );
}

// start_lost / start_retained_variant asymmetry: Perl VEP's
// `VariationEffect::inframe_insertion` returns 0 when start_lost fires, so the
// Met-loss gate in classify_inframe_indel_by_peptide suppresses
// inframe_insertion, while `inframe_deletion` has no start_lost gate: an inframe
// deletion at the start codon that preserves codon prefix/suffix alignment
// reports inframe_deletion + start_lost. The Met-loss gate therefore sits in the
// insertion branch only.

/// A pure 3bp insertion at CDS pos 1 where the inserted codon
/// shifts ATG out of the first position; with the Met-loss gate in the
/// insertion branch, inframe_insertion is suppressed in favour of
/// ProteinAltering so that the downstream start_lost emission in
/// consequences.rs stands alone for peptide-based classification.
///
/// Insert "CCC" immediately before the original ATG. Alt CDS becomes
/// "CCCATGGCT...", so the first codon is CCC (Pro): Met lost in the
/// first codon.
#[test]
fn concordance_insertion_before_start_codon_met_loss_suppresses_inframe_insertion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Pure insertion before CDS pos 1 (genomic 25_000_050). VEP insertion
    // convention: ref="-", end < start. Inserting at the CDS start puts
    // the new codon CCC in front of ATGGCT → alt first codon is CCC.
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_049,
        b"-".to_vec(),
        b"CCC".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::InframeInsertion),
        "start-codon insertion with Met loss must not emit inframe_insertion, got: {:?}",
        tc.consequences
    );
}

/// A length-change deletion at the start codon where codons prefix/suffix
/// match must emit inframe_deletion (the codon-containment test of Perl's
/// `inframe_deletion`), not protein_altering: the Met-loss gate is in the
/// insertion branch only.
#[test]
fn concordance_length_change_deletion_at_start_codon_keeps_inframe_deletion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // 6bp deletion starting at CDS pos 1. Ref="ATGGCT" → Alt="ACT" yields a
    // net -3bp deletion where Perl's trim_sequences collapses to ref="ATG"
    // alt="" (3bp deletion, divisible by 3), which emits inframe_deletion even
    // though Met is lost in the peptide translation of the local codon
    // window: the Met-loss gate must not fire before the codon-containment
    // check in classify_inframe_indel_by_peptide.
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_055,
        b"ATGGCT".to_vec(),
        b"ACT".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    // The consequence set should contain inframe_deletion (Perl-correct) and
    // may additionally contain start_lost from the downstream gate in
    // consequences.rs. It must not contain protein_altering_variant alongside
    // start_lost because Perl's protein_altering_variant predicate is
    // suppressed when start_lost or inframe_deletion fires.
    let has_inframe_del = tc.consequences.contains(&Consequence::InframeDeletion);
    let has_protein_altering = tc
        .consequences
        .contains(&Consequence::ProteinAlteringVariant);

    assert!(
        has_inframe_del || !has_protein_altering,
        "length-change deletion must emit inframe_deletion or suppress \
         protein_altering_variant, got: {:?}",
        tc.consequences
    );
}

/// splice_donor_region_variant vs splice_donor_5th_base_variant swap
/// on complex delins where the 5th-base position matches between REF and ALT.
///
/// Perl's `_intron_effects` iterates the differing regions from
/// `_get_differing_regions` and sets `fifth_base_splice_site` only if a
/// differing region overlaps it, so a multi-base delin whose matching bases sit
/// at the 5th-base position gets donor_region, not 5th_base; the full variant
/// span (extended by net_ext) bridges the matching gap and would set 5th_base.
///
/// Shape (TP53, reverse strand, intron 7670716-7673534):
///   Variant: 17:7673526-7673531 ACTTAG to GGTGAAA
///   Perl differing regions: [7673526-7673527], [7673529], [7673531-7673532]
///   5th_base_reverse position: 7673530 (intron_end - 4), no region overlaps
///   donor_region_reverse window: 7673529-7673532, regions 2 and 3 overlap
///   Perl: intron_variant,splice_donor_region_variant
#[test]
fn concordance_splice_donor_region_vs_5th_base_forward_delins() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Forward strand test transcript, intron 1 is 25_000_300-25_001_999.
    // 5th base: intron_start + 4 = 25_000_304 (single position)
    // Donor region: intron_start + 2..intron_start + 5 = 25_000_302..25_000_305
    //
    // Delins at 25_000_302-25_000_307 (6bp ref), same-length alt (6bp).
    // REF[2] and ALT[2] both 'G' so XOR matches → position 25_000_304 (5th base).
    // REF[3] and ALT[3] both 'C' so XOR matches → position 25_000_305.
    // Other positions differ.
    //
    // Differing regions (Perl semantics):
    //   - [25_000_302, 25_000_303] (positions 0-1)
    //   - [25_000_306, 25_000_307] (positions 4-5)
    //
    // 5th base (25_000_304): no region overlaps → fifth_base_splice_site = 0
    // Donor region (25_000_302-25_000_305): Region 1 overlaps → donor_region = 1
    //
    // Perl result: splice_donor_region_variant + intron_variant
    let variant = InputVariant::new(
        "21".into(),
        25_000_302,
        25_000_307,
        b"ATGCAG".to_vec(), // REF[2]='G', REF[3]='C'
        b"CGGCTT".to_vec(), // ALT[2]='G' (match), ALT[3]='C' (match)
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences
            .contains(&Consequence::SpliceDonorRegionVariant),
        "complex delins with matching base at 5th-base position should emit \
         splice_donor_region_variant (Perl parity), got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant),
        "5th base position matches between REF and ALT so per-region overlap \
         from Perl's _get_differing_regions should NOT set fifth_base_splice_site. \
         Got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::IntronVariant),
        "delins reaching intron positions still emits intron_variant, got: {:?}",
        tc.consequences
    );
}

/// Companion to the TP53 shape above, on the forward test transcript, whose
/// intron 1 is 25_000_300-25_001_999.
/// This places a complex delins
///
/// with a single matching base at the 5th-base position,
/// to cover the case where the differing regions are non-contiguous and span
/// the 5th-base window on both sides without hitting it.
#[test]
fn concordance_splice_donor_region_vs_5th_base_noncontiguous() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Intron 1: 25_000_300-25_001_999 (forward strand).
    // 5th base: 25_000_304, donor region: 25_000_302-25_000_305.
    //
    // Delins at 25_000_300-25_000_306 (7bp ref), same-length alt.
    // REF[4] and ALT[4] both 'G' → position 25_000_304 (5th base) matches.
    // Other positions differ.
    //
    // Differing regions:
    //   - [25_000_300, 25_000_303] (positions 0-3)
    //   - [25_000_305, 25_000_306] (positions 5-6)
    //
    // 5th base (25_000_304): no region overlaps.
    // Donor region (25_000_302-25_000_305): Region 1 overlaps 25_000_302-25_000_303.
    // Region 2 overlaps 25_000_305.
    // Splice donor (25_000_300-25_000_301): Region 1 overlaps.
    //
    // Perl result: splice_donor_variant,splice_donor_region_variant,intron_variant
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_306,
        b"ACGTGAT".to_vec(), // REF[4]='G'
        b"TGCAGCA".to_vec(), // ALT[4]='G' (match)
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::SpliceDonorVariant),
        "variant covering intron start should include splice_donor_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences
            .contains(&Consequence::SpliceDonorRegionVariant),
        "differing regions flanking the matched 5th base should still hit the \
         donor region window, got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant),
        "matching base at 5th-base position must suppress \
         splice_donor_5th_base_variant (Perl parity), got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::IntronVariant),
        "intronic overlap always fires intron_variant, got: {:?}",
        tc.consequences
    );
}

// 3'UTR stop_retained vs stop_lost: replicate_ins_del_stop_altered

/// Build a transcript whose CDS ends in a TAA stop codon followed by a 3'UTR
/// beginning with "AAGAAA...". A 4bp deletion that removes the second/third
/// base of the stop plus the first two UTR bases shifts the next two UTR bases
/// (A,A) into the stop codon position, so the new codon at the original stop
/// position is `T + A + A = TAA` (still a stop). Perl emits
/// `stop_retained_variant,3_prime_UTR_variant` for this shape.
fn make_utr_stop_shift_transcript() -> Transcript {
    use crate::test_helpers::make_test_transcript;
    let mut tx = make_test_transcript();
    // Rewrite the last 3 bases of the translateable_seq (CDS) to TAA. The
    // test CDS is 850bp long; positions 848-850 form the stop codon.
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(seq) = vefc.translateable_seq.as_mut() {
            let n = seq.len();
            assert!(n >= 3);
            seq.replace_range((n - 3)..n, "TAA");
        }
    }
    tx
}

/// Build a synthetic FASTA whose layout matches `make_test_transcript()`:
/// exons at 25_000_000-25_000_299, 25_002_000-25_002_299, 25_004_000-25_006_000.
/// Intergenic filler is 'N'; only exon bases matter. The CDS ends at genomic
/// 25_004_299 (cDNA 900). The stop codon TAA sits at genomic
/// 25_004_297-25_004_299 and the first 4 UTR bases are AAGA, so a 4bp deletion
/// starting at 25_004_298 removes AAGA and shifts UTR bases into the CDS. The
/// 2 bases past the deletion at the stop position, genomic 25_004_302-25_004_303,
/// are A,A so that the new codon at cDNA 898-900 is T + A + A = TAA, still a
/// stop (Perl: stop_retained_variant).
fn build_utr_stop_shift_fasta() -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
    // Build a chr21 sequence spanning 25_000_000 through 25_006_100 with N
    // padding; overwrite the exonic bases where the test needs specific ones.
    let total_len: usize = 25_006_100;
    let mut seq = vec![b'N'; total_len];

    // Exon 1 (25_000_000-25_000_299): any bases fine
    for pos in 25_000_000..=25_000_299usize {
        seq[pos - 1] = b'A';
    }
    // Exon 2 (25_002_000-25_002_299)
    for pos in 25_002_000..=25_002_299usize {
        seq[pos - 1] = b'C';
    }
    // Exon 3 (25_004_000-25_006_000)
    for pos in 25_004_000..=25_006_000usize {
        seq[pos - 1] = b'G';
    }
    // Place TAA at stop codon position (cDNA 898-900 = genomic
    // 25_004_297-25_004_299; because CDS begins at cDNA 51, stop codon is
    // cds_pos 848-850, cDNA 898-900, genomic 25_004_297-25_004_299).
    seq[25_004_297 - 1] = b'T';
    seq[25_004_298 - 1] = b'A';
    seq[25_004_299 - 1] = b'A';
    // UTR bases at genomic 25_004_300..25_004_303 = "AAGA".
    // The 4bp deletion (TAAGA>T at genomic 25_004_298..25_004_302) removes
    // bases at 25_004_299..25_004_302 = "AGAA" (VEP convention: anchor T
    // kept, 4 bases deleted).
    //
    // The record used below is anchored at 25_004_297: REF="TAAGA" ALT="T",
    // which in VEP coordinates becomes a 4bp deletion at 25_004_298..25_004_301
    // of bases AAGA. After normalization ref="AAGA" alt="-" start=25_004_298
    // end=25_004_301.  Stop codon in cDNA is 898-900 = genomic
    // 25_004_297-25_004_299 = T,A,A.  The deletion removes
    // genomic 25_004_298,25_004_299,25_004_300,25_004_301 = A,A,?,? .
    // After edit, genomic 25_004_298 takes on the value at 25_004_302 and
    // genomic 25_004_299 takes on the value at 25_004_303. For the codon at
    // the original stop position (cDNA 898-900 = genomic 25_004_297-299) to
    // remain TAA, genomic 25_004_302 and 25_004_303 must both be A.
    // Deletion REF="AAGA" covers genomic 25_004_298..25_004_301, so the
    // bases at 25_004_298..25_004_301 must be A,A,G,A. The first two are
    // already A,A (stop codon 2nd and 3rd base from TAA placement above).
    seq[25_004_300 - 1] = b'G';
    seq[25_004_301 - 1] = b'A';
    // After deletion, genomic 25_004_298 takes on what was at 25_004_302
    // and 25_004_299 takes on what was at 25_004_303. For the codon at
    // the original stop position (genomic 25_004_297..25_004_299 =
    // cDNA 898-900) to remain `T + A + A = TAA`, 25_004_302 and 25_004_303 must
    // both be A.
    seq[25_004_302 - 1] = b'A';
    seq[25_004_303 - 1] = b'A';

    // Write FASTA + FAI.
    let dir = tempdir().expect("tempdir");
    let fa = dir.path().join("ref.fa");
    let fai = dir.path().join("ref.fa.fai");
    let header = ">21\n";
    let mut contents = Vec::with_capacity(header.len() + seq.len() + 1);
    contents.extend_from_slice(header.as_bytes());
    contents.extend_from_slice(&seq);
    contents.push(b'\n');
    std::fs::write(&fa, &contents).expect("write FASTA");
    std::fs::write(
        &fai,
        format!(
            "21\t{}\t{}\t{}\t{}\n",
            total_len,
            header.len(),
            total_len,
            total_len + 1
        ),
    )
    .expect("write FASTA index");
    let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("load indexed FASTA");
    (dir, fasta)
}

/// UTR-shifted stop codon preservation.
///
/// Transcript has stop codon `TAA` at the end of CDS. A 4bp deletion
/// starting at cds_pos 849 (middle of stop) extends into 3'UTR. The two
/// UTR bases that shift into the stop position are `A,A`, so the edited
/// codon at the original stop position is `T+A+A = TAA`, still a stop.
/// Perl emits `stop_retained_variant,3_prime_UTR_variant`.
#[test]
fn concordance_utr_shifted_stop_codon_is_stop_retained() {
    let tx = make_utr_stop_shift_transcript();
    let (_tmp, fasta) = build_utr_stop_shift_fasta();
    let config = EffectsConfig {
        reference_fasta: Some(std::sync::Arc::new(fasta)),
        ..Default::default()
    };

    // 4bp deletion at genomic 25_004_298..25_004_301 (VCF: 25_004_297 TAAGA>T).
    // Per VEP convention, the anchor T is preserved and 4 bases AAGA are
    // deleted: start=25_004_298 end=25_004_301 ref="AAGA" alt="-".
    let variant = InputVariant::new(
        "21".into(),
        25_004_298,
        25_004_301,
        b"AAGA".to_vec(),
        b"-".to_vec(),
    );

    let result =
        calculate_consequences(&variant, &tx, &config).expect("should produce consequences");

    assert!(
        result
            .consequences
            .contains(&Consequence::StopRetainedVariant),
        "UTR-shifted stop preservation must emit stop_retained_variant (got {:?})",
        result.consequences
    );
    assert!(
        !result.consequences.contains(&Consequence::StopLost),
        "must not emit stop_lost when the shifted codon is still a stop (got {:?})",
        result.consequences
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant),
        "3' UTR span must co-emit 3_prime_UTR_variant (got {:?})",
        result.consequences
    );
}

/// UTR bases that do not form a stop yield stop_lost.
///
/// Same transcript + variant, but the UTR bases that shift into the stop
/// codon position make a codon that is not a stop (`C,C`). Perl emits
/// `stop_lost,3_prime_UTR_variant`.
#[test]
fn concordance_utr_shifted_codon_that_is_not_a_stop_is_stop_lost() {
    let tx = make_utr_stop_shift_transcript();
    let (dir, _fasta) = build_utr_stop_shift_fasta();

    // Rewrite the UTR bases that shift into the stop codon position so the
    // new codon is not a stop. After the 4bp deletion at genomic
    // 25_004_298..25_004_301 (removing "AAGA"), the first two UTR bases
    // (originally at 25_004_302 and 25_004_303) shift into CDS positions
    // 849 and 850. `replicate_ins_del_stop_altered` splices 4 bytes out of
    // `combined = CDS + 3'UTR` at CDS
    // anchor 848. After the 4-byte splice the codon at cds.len()-3..cds.len()
    // (original stop at combined[847..850]) reads:
    //   combined[847] = T (still first base of the original TAA stop)
    //   combined[848] = first UTR base that shifted in = original 25_004_302
    //   combined[849] = second UTR base that shifted in = original 25_004_303
    //
    // So the negative test must rewrite 25_004_302 and 25_004_303 to bases
    // that do not spell a stop codon together with T. T+C+C = TCC (Ser) is
    // not a stop, so rewrite both to 'C'.
    let fa = dir.path().join("ref.fa");
    let mut contents = std::fs::read(&fa).expect("read FASTA");
    let header_len = b">21\n".len();
    let pos_302_idx = header_len + (25_004_302 - 1);
    let pos_303_idx = header_len + (25_004_303 - 1);
    // Sanity: both positions are 'A' per build_utr_stop_shift_fasta.
    assert_eq!(
        contents[pos_302_idx], b'A',
        "expected genomic 25_004_302 = 'A' before rewrite"
    );
    assert_eq!(
        contents[pos_303_idx], b'A',
        "expected genomic 25_004_303 = 'A' before rewrite"
    );
    contents[pos_302_idx] = b'C';
    contents[pos_303_idx] = b'C';
    std::fs::write(&fa, &contents).expect("rewrite FASTA");
    let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("reload FASTA");
    assert_eq!(
        fasta.base("21", 25_004_302),
        Some(b'C'),
        "IndexedFasta must see rewritten 'C' at 25_004_302"
    );
    assert_eq!(
        fasta.base("21", 25_004_303),
        Some(b'C'),
        "IndexedFasta must see rewritten 'C' at 25_004_303"
    );

    let config = EffectsConfig {
        reference_fasta: Some(std::sync::Arc::new(fasta)),
        ..Default::default()
    };

    let variant = InputVariant::new(
        "21".into(),
        25_004_298,
        25_004_301,
        b"AAGA".to_vec(),
        b"-".to_vec(),
    );

    let result =
        calculate_consequences(&variant, &tx, &config).expect("should produce consequences");

    assert!(
        result.consequences.contains(&Consequence::StopLost),
        "shifted codon TCC is not a stop, so stop_lost must be emitted (got {:?})",
        result.consequences
    );
    assert!(
        !result
            .consequences
            .contains(&Consequence::StopRetainedVariant),
        "must not emit stop_retained when shifted codon is not a stop (got {:?})",
        result.consequences
    );
}

/// Perl parity requires splicing the full cDNA-based length of the edit into
/// the combined CDS+3'UTR string, not a length clamped at the CDS end: Perl's
/// `VariationEffect::_ins_del_stop_altered` hard-codes
/// `($cdna_end - $cdna_start) + 1` as the splice length while anchoring at
/// `$cds_start - 1`, so the splice is longer than the remaining CDS. A 4bp
/// deletion whose last byte crosses into the 3'UTR spliced as only 3 bytes
/// (the CDS portion) makes `combined[stop_idx..stop_idx+3]` read positions
/// that differ from Perl's 4-byte splice.
#[test]
fn concordance_stop_altered_splice_length_not_clamped_at_cds_end() {
    let tx = make_utr_stop_shift_transcript();
    let (_tmp, fasta) = build_utr_stop_shift_fasta();
    let config = EffectsConfig {
        reference_fasta: Some(std::sync::Arc::new(fasta)),
        ..Default::default()
    };

    // Same variant as concordance_utr_shifted_stop_codon_is_stop_retained: a 4bp
    // deletion that ends 2 bases past the CDS stop codon. With the full cDNA
    // length (4), combined[848..850] reads the bases 2 positions further along
    // at 25_004_302/303 = A,A, so T+A+A = TAA stop_retained. A splice clamped
    // at cds end (3 bytes) would read the deleted-but-shifted UTR bases at
    // genomic 25_004_300/301, G,A, and produce `stop_lost`.
    let variant = InputVariant::new(
        "21".into(),
        25_004_298,
        25_004_301,
        b"AAGA".to_vec(),
        b"-".to_vec(),
    );

    let result =
        calculate_consequences(&variant, &tx, &config).expect("should produce consequences");
    assert!(
        result
            .consequences
            .contains(&Consequence::StopRetainedVariant),
        "full-cDNA splice must emit stop_retained_variant for this variant (got {:?})",
        result.consequences
    );
    assert!(
        !result.consequences.contains(&Consequence::StopLost),
        "a clamped splice length would emit stop_lost here (got {:?})",
        result.consequences
    );
}

// 5'UTR start_lost DNA check (FASTA-backed `_ins_del_start_altered` parity):
// Perl's `start_lost` predicate relies on `_ins_del_start_altered`, which
// checks whether the ATG triplet at `length($utr->seq)` is still 'ATG' after
// the edit (and the 5'UTR bases are unchanged), so a UTR-spanning
// length-changing deletion that leaves the ATG intact does not emit
// `start_lost`.

/// Build a synthetic FASTA that matches `make_test_transcript()`'s coding
/// model with a known 5'UTR: cDNA 1..50 lives in exon 1 at genomic
/// 25_000_000..25_000_049. The CDS ATG occupies cDNA 51..53 = genomic
/// 25_000_050..25_000_052.  UTR bases are 'T' (so a deletion entirely inside
/// the UTR that leaves the ATG in place is easy to construct).
fn build_five_prime_utr_atg_fasta() -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
    let total_len: usize = 25_006_100;
    let mut seq = vec![b'N'; total_len];

    // 5'UTR bases T (cDNA 1..50 = genomic 25_000_000..25_000_049).
    for pos in 25_000_000..=25_000_049usize {
        seq[pos - 1] = b'T';
    }
    // ATG at cDNA 51..53 = genomic 25_000_050..25_000_052.
    seq[25_000_050 - 1] = b'A';
    seq[25_000_051 - 1] = b'T';
    seq[25_000_052 - 1] = b'G';
    // Fill the rest of exon 1 (cDNA 54..300 = genomic 25_000_053..25_000_299)
    // with C so the peptide content doesn't accidentally form another start
    // codon.
    for pos in 25_000_053..=25_000_299usize {
        seq[pos - 1] = b'C';
    }
    for pos in 25_002_000..=25_002_299usize {
        seq[pos - 1] = b'C';
    }
    for pos in 25_004_000..=25_006_000usize {
        seq[pos - 1] = b'G';
    }

    let dir = tempdir().expect("tempdir");
    let fa = dir.path().join("ref.fa");
    let fai = dir.path().join("ref.fa.fai");
    let header = ">21\n";
    let mut contents = Vec::with_capacity(header.len() + seq.len() + 1);
    contents.extend_from_slice(header.as_bytes());
    contents.extend_from_slice(&seq);
    contents.push(b'\n');
    std::fs::write(&fa, &contents).expect("write FASTA");
    std::fs::write(
        &fai,
        format!(
            "21\t{}\t{}\t{}\t{}\n",
            total_len,
            header.len(),
            total_len,
            total_len + 1
        ),
    )
    .expect("write FASTA index");
    let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("load FASTA");
    (dir, fasta)
}

/// A deletion entirely inside the 5'UTR leaves the ATG at cDNA 51..53
/// untouched. After the edit, `combined = UTR[1..46] + CDS`, so the last 850
/// bytes equal the original CDS and Perl's `_ins_del_start_altered` falls
/// through to the tail-comparison check and returns 0 (start not altered):
/// no `start_lost`.
///
/// The variant spans into CDS (end coordinate 25_000_050 = cDNA 51) so that
/// `is_start_lost` fires before the FASTA-backed helper runs; otherwise the
/// `spans_non_coding` branch would not evaluate the start-codon path.
#[test]
fn concordance_five_prime_utr_deletion_preserving_atg_no_start_lost() {
    let tx = {
        use crate::test_helpers::make_test_transcript;
        make_test_transcript()
    };
    let (_tmp, fasta) = build_five_prime_utr_atg_fasta();
    let config = EffectsConfig {
        reference_fasta: Some(std::sync::Arc::new(fasta)),
        ..Default::default()
    };

    // 4bp deletion at genomic 25_000_046..25_000_049 (cDNA 47..50, entirely
    // inside the 5'UTR). REF="TTTT" ALT="-". After edit combined =
    // UTR[1..46] + CDS (46 + 850 = 896 bytes). Last 850 bytes = CDS
    // unchanged → Perl returns 0, start not altered.
    //
    // Because the variant end (25_000_049) is the last UTR base and the
    // primary endpoint of the is_start_lost gate maps through CDS anchor
    // logic, this straddles the UTR/CDS boundary enough to fire the
    // `spans_non_coding` + `is_start_lost` branch.
    let variant = InputVariant::new(
        "21".into(),
        25_000_046,
        25_000_049,
        b"TTTT".to_vec(),
        b"-".to_vec(),
    );

    let result =
        calculate_consequences(&variant, &tx, &config).expect("should produce consequences");

    assert!(
        !result.consequences.contains(&Consequence::StartLost),
        "UTR-only deletion with intact ATG must NOT emit start_lost. \
         got {:?}",
        result.consequences
    );
}

/// A deletion spanning the UTR-to-CDS boundary that consumes the 'A' of the
/// ATG start codon. After the edit, `combined[50..53]` is not 'ATG' and the
/// tail comparison fails, so Perl returns 1, start altered: `start_lost`.
#[test]
fn concordance_five_prime_utr_deletion_destroying_atg_emits_start_lost() {
    let tx = {
        use crate::test_helpers::make_test_transcript;
        make_test_transcript()
    };
    let (_tmp, fasta) = build_five_prime_utr_atg_fasta();
    let config = EffectsConfig {
        reference_fasta: Some(std::sync::Arc::new(fasta)),
        ..Default::default()
    };

    // 4bp deletion at genomic 25_000_047..25_000_050 (cDNA 48..51, straddles
    // UTR→CDS boundary). REF="TTTA" ALT="-". Removes UTR[48..50] (bases
    // T,T,T) plus CDS[1] (the 'A' of ATG). After edit combined =
    // UTR[1..47] + CDS[2..] = 47 + 849 = 896 bytes. Last 850 bytes start
    // with 'T' (from CDS[2] = 'T' of ATG), not 'A' → tail != CDS → altered.
    let variant = InputVariant::new(
        "21".into(),
        25_000_047,
        25_000_050,
        b"TTTA".to_vec(),
        b"-".to_vec(),
    );

    let result =
        calculate_consequences(&variant, &tx, &config).expect("should produce consequences");

    assert!(
        result.consequences.contains(&Consequence::StartLost),
        "deletion consuming the 'A' of ATG must emit start_lost. \
         got {:?}",
        result.consequences
    );
}

// Frameshift / inframe preserved across a splice. The predicate these
// reproduce is Perl VEP's own, in `Utils::VariationEffect`; each case below
// states the span shape it pins.
//
// Test transcript layout:
//   Exon 1: genomic 25_000_000..25_000_299 (cDNA 1..300)
//   Intron 1: genomic 25_000_300..25_001_999
//   Exon 2: genomic 25_002_000..25_002_299 (cDNA 301..600)
//   Intron 2: genomic 25_002_300..25_003_999
//   Exon 3: genomic 25_004_000..25_006_000 (cDNA 601..2601)
//   cDNA CDS: 51..900, so CDS pos = cdna_pos - 50
//
// The target shape is [Coordinate, Gap, Coordinate]: a deletion spanning an
// intron with both endpoints in CDS exons. Perl's
// `BaseTranscriptVariation::cds_start` sets `cds_start` / `cds_end` from the
// first/last Coord, both non-Gap, so the `defined` gate of
// `VariationEffect::frameshift` opens. `map_genomic_span_to_cds_bounds` clamps at
// CDS boundaries and masks the Gap-vs-Coord distinction, so the gate reads the
// segment shape from `map_genomic_span_to_cds_projection` (CdsSpanProjection).

/// A deletion spanning exon 1 through intron 1 into exon 2, with both
/// endpoints in CDS exons, where `(var_len - allele_len) % 3 != 0`: Perl
/// emits `frameshift_variant`.
///
/// Span: genomic 25_000_249..25_002_049 = cDNA 250..350 = CDS 200..300.
/// var_len = 300 - 200 + 1 = 101. Pure deletion (alt = "-") allele_len = 0.
/// abs(101 - 0) = 101. 101 % 3 = 2 → frameshift.
#[test]
fn concordance_frameshift_across_intron_both_endpoints_in_cds() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Construct a deletion with a long ref allele (1801 bp) across an intron.
    // VEP convention for deletions: ref = deleted sequence, alt = "-".
    let ref_len = (25_002_049u64 - 25_000_249u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let variant = InputVariant::new(
        "21".into(),
        25_000_249,
        25_002_049,
        ref_bytes,
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("frameshift across an intron should produce consequences");

    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "[C, G, C] shape with var_len=101, allele_len=0, diff%3=2 \
         must emit frameshift_variant. got {:?}",
        result.consequences
    );
    // Perl's `coding_sequence_variant` = `VariationEffect::coding_unknown`
    // fires only when `frameshift` / `inframe_*` /
    // `protein_altering_variant` / `start_*` / `stop_*` all return 0, so
    // when frameshift fires, coding_sequence_variant must not be emitted.
    assert!(
        !result
            .consequences
            .contains(&Consequence::CodingSequenceVariant),
        "when FrameshiftVariant is emitted, CodingSequenceVariant \
         must be suppressed (Perl's `coding_unknown` predicate excludes \
         frameshift). got {:?}",
        result.consequences
    );
}

/// A delins spanning exon 1 through intron 1 into exon 2, both endpoints in
/// CDS exons, with a net length change divisible by 3 whose alternate bases do
/// not match the codons they replace: Perl emits `protein_altering_variant`.
///
/// `VariationEffect::frameshift` is 0 (`abs(98 - 101) % 3 == 0`), but
/// `inframe_deletion` compares CODONS: the alt codon window
/// (`codon_len + allele_len - vf_nt_len` bases of the alternate CDS,
/// `TranscriptVariationAllele::codon`) is 99 `G`s, neither a prefix nor a suffix
/// of the ref window and not reducible by `trim_sequences`, so it returns 0 and
/// `protein_altering_variant` takes the allele: the peptides differ
/// in length, neither starts with `*`, and the alt peptide (all Gly) contains
/// the ref peptide (all Ala) at neither end.
///
/// Span: cDNA 250..350, CDS 200..300, var_len = 101, allele_len = 98.
#[test]
fn concordance_inframe_delins_across_intron_is_protein_altering() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    let ref_len = (25_002_049u64 - 25_000_249u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let alt_bytes = vec![b'G'; 98];
    let variant = InputVariant::new("21".into(), 25_000_249, 25_002_049, ref_bytes, alt_bytes);

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("inframe delins across an intron should produce consequences");

    assert!(
        result
            .consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "[C, G, C] delins with diff%3=0 whose alt bases match no ref codon \
         must emit protein_altering_variant. got {:?}",
        result.consequences
    );
    for absent in [
        Consequence::InframeDeletion,
        Consequence::FrameshiftVariant,
        Consequence::CodingSequenceVariant,
    ] {
        assert!(
            !result.consequences.contains(&absent),
            "{absent:?} must not accompany protein_altering_variant. got {:?}",
            result.consequences
        );
    }
}

/// A pure deletion of whole codons spanning exon 1 through intron 1 into exon
/// 2, both endpoints in CDS exons: Perl emits `inframe_deletion`.
///
/// `cds_start` / `cds_end` are 199 / 300 (both `Coordinate`s, the intron a
/// `Gap` between them), `var_len` is 102 so `frameshift` is 0, and the alt codon
/// window is empty (`codon_len + 0 - 102 == 0`), which
/// `VariationEffect::inframe_deletion` accepts as a prefix of the ref window.
#[test]
fn concordance_pure_codon_deletion_across_intron_is_inframe_deletion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // cDNA 249 (exon 1, genomic 25_000_248) .. cDNA 350 (exon 2, genomic
    // 25_002_049) = CDS 199..300, 34 whole codons.
    let ref_len = (25_002_049u64 - 25_000_248u64 + 1) as usize;
    let variant = InputVariant::new(
        "21".into(),
        25_000_248,
        25_002_049,
        vec![b'A'; ref_len],
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("pure codon deletion across an intron should produce consequences");

    assert!(
        result.consequences.contains(&Consequence::InframeDeletion),
        "whole-codon [C, G, C] deletion must emit inframe_deletion. got {:?}",
        result.consequences
    );
    for absent in [
        Consequence::ProteinAlteringVariant,
        Consequence::FrameshiftVariant,
        Consequence::CodingSequenceVariant,
    ] {
        assert!(
            !result.consequences.contains(&absent),
            "{absent:?} must not accompany inframe_deletion. got {:?}",
            result.consequences
        );
    }
}

/// The `[C, G]` shape (exon 1 CDS into intron 1) must not emit frameshift.
/// Perl's `frameshift` blocks this on its `defined cds_end` gate because the
/// last segment is a Gap, so `cds_end = undef`.
#[test]
fn concordance_exon_to_intron_deletion_no_frameshift() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Deletion from exon 1 CDS (25_000_200) into intron 1 (25_000_500).
    // Start endpoint is InCds, end endpoint is InIntron → gate blocks.
    let ref_len = (25_000_500u64 - 25_000_200u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let variant = InputVariant::new(
        "21".into(),
        25_000_200,
        25_000_500,
        ref_bytes,
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("exon→intron deletion should produce consequences");

    assert!(
        !result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "[C, G] shape must NOT emit frameshift_variant \
         (end endpoint in intron → Perl cds_end=undef → gate blocks). \
         got {:?}",
        result.consequences
    );
    assert!(
        !result.consequences.contains(&Consequence::InframeInsertion),
        "[C, G] shape must not emit inframe_insertion. got {:?}",
        result.consequences
    );
    assert!(
        !result.consequences.contains(&Consequence::InframeDeletion),
        "[C, G] shape must not emit inframe_deletion. got {:?}",
        result.consequences
    );
}

/// The `[G, C]` shape (intron 1 into exon 2 CDS) must not emit frameshift.
/// Same gate as the exon-to-intron case but from the other side: first
/// segment is Gap, so `cds_start = undef`.
#[test]
fn concordance_intron_to_exon_deletion_no_frameshift() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Deletion from intron 1 (25_001_500) into exon 2 CDS (25_002_100).
    let ref_len = (25_002_100u64 - 25_001_500u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let variant = InputVariant::new(
        "21".into(),
        25_001_500,
        25_002_100,
        ref_bytes,
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("intron→exon deletion should produce consequences");

    assert!(
        !result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "[G, C] shape must NOT emit frameshift_variant \
         (start endpoint in intron → Perl cds_start=undef → gate blocks). \
         got {:?}",
        result.consequences
    );
    assert!(
        !result.consequences.contains(&Consequence::InframeDeletion),
        "[G, C] shape must not emit inframe_deletion. got {:?}",
        result.consequences
    );
}

/// A wholly intronic span, the `[Gap]` shape, must not emit any coding
/// consequence. Both endpoints InIntron.
#[test]
fn concordance_wholly_intronic_deletion_no_frameshift() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Deletion wholly within intron 1 (25_000_500..25_001_500).
    let ref_len = (25_001_500u64 - 25_000_500u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let variant = InputVariant::new(
        "21".into(),
        25_000_500,
        25_001_500,
        ref_bytes,
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    // Result may be None for wholly intronic non-coding-altering spans, or
    // Some with only intron_variant. Either way: no frameshift.
    if let Some(result) = result {
        assert!(
            !result
                .consequences
                .contains(&Consequence::FrameshiftVariant),
            "wholly intronic must not emit frameshift_variant. got {:?}",
            result.consequences
        );
    }
}

/// The `[C, G]` shape (exon 1 CDS into intron 1) with a size-3 deletion
/// must not emit inframe_deletion. Even though the net length change is a
/// multiple of 3, the gate blocks because end_endpoint is InIntron.
#[test]
fn concordance_exon_to_intron_delins_no_inframe_deletion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Span 301 bp (cDNA 201..300 + intron 1), alt = 298 bp → diff=3, %3=0.
    let ref_len = (25_000_500u64 - 25_000_200u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let alt_bytes = vec![b'G'; ref_len - 3];
    let variant = InputVariant::new("21".into(), 25_000_200, 25_000_500, ref_bytes, alt_bytes);

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("exon→intron delins should produce consequences");

    assert!(
        !result.consequences.contains(&Consequence::InframeDeletion),
        "[C, G] shape must not emit inframe_deletion. got {:?}",
        result.consequences
    );
    assert!(
        !result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "[C, G] shape must not emit frameshift_variant. got {:?}",
        result.consequences
    );
}

/// A deletion spanning CDS into 3'UTR (exon 3 CDS to 3'UTR), the
/// `[Coordinate, Gap]` shape where the Gap is trailing UTR: the is_stop_lost
/// branch handles it and the spans_non_coding else branch must not
/// double-emit frameshift.
#[test]
fn concordance_cds_to_utr_deletion_no_spurious_frameshift() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Deletion from CDS end (exon 3, genomic 25_004_295) into 3'UTR
    // (25_004_305). 11 bp ref, alt = "-", is_stop_lost path handles.
    let variant = InputVariant::new(
        "21".into(),
        25_004_295,
        25_004_305,
        b"ACGTACGTACG".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("CDS→3'UTR deletion should produce consequences");

    // Must keep the 3_prime_UTR_variant and must not add a frameshift_variant.
    assert!(
        result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant),
        "must keep 3_prime_UTR_variant emission. got {:?}",
        result.consequences
    );
    assert!(
        !result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "CDS→3'UTR must NOT emit frameshift_variant (is_stop_lost \
         path owns this case; spans_non_coding else gate must block). \
         got {:?}",
        result.consequences
    );
}

/// An exon-spanning delins over two introns whose net length change is 3 and
/// whose 449 alternate bases match none of the replaced codons: Perl's
/// `inframe_deletion` codon-containment test fails and
/// `protein_altering_variant` is emitted (both in `VariationEffect`).
#[test]
fn concordance_inframe_delins_across_three_exons_is_protein_altering() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Span exon 1 -> exon 3 (two introns). CDS 200..651. var_len = 452.
    // Alt = 449 bp -> diff = 3, divisible by 3: not a frameshift.
    let ref_len = (25_004_100u64 - 25_000_249u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let alt_bytes = vec![b'G'; 449];
    let variant = InputVariant::new("21".into(), 25_000_249, 25_004_100, ref_bytes, alt_bytes);

    let result = calculate_consequences(&variant, &tx, &config)
        .expect("inframe delins across three exons should produce consequences");

    assert!(
        result
            .consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "[C, G, C, G, C] delins with diff%3=0 and unmatched alt codons must emit \
         protein_altering_variant. got {:?}",
        result.consequences
    );
    for absent in [
        Consequence::InframeDeletion,
        Consequence::FrameshiftVariant,
        Consequence::CodingSequenceVariant,
    ] {
        assert!(
            !result.consequences.contains(&absent),
            "{absent:?} must not accompany protein_altering_variant. got {:?}",
            result.consequences
        );
    }
}

/// Guard against `[C, G, C]` emitting frameshift when the variant overlaps
/// the stop codon: `is_stop_lost` owns that case, not the `spans_non_coding`
/// else gate.
#[test]
fn concordance_deletion_to_stop_codon_never_emits_both_stop_lost_and_frameshift() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Deletion exon 2 CDS → exon 3 CDS (crosses intron 2), hitting the stop
    // codon. Span: genomic 25_002_200..25_004_299 (CDS 451..850 → entire
    // second half of CDS ending at stop codon).
    let ref_len = (25_004_299u64 - 25_002_200u64 + 1) as usize;
    let ref_bytes = vec![b'A'; ref_len];
    let variant = InputVariant::new(
        "21".into(),
        25_002_200,
        25_004_299,
        ref_bytes,
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    // The stop_lost logic or the frameshift gate may fire, but never both
    // stop_lost and frameshift_variant.
    if let Some(result) = result {
        let has_stop_lost = result.consequences.contains(&Consequence::StopLost);
        let has_frameshift = result
            .consequences
            .contains(&Consequence::FrameshiftVariant);
        assert!(
            !(has_stop_lost && has_frameshift),
            "stop-lost guard: must not emit both stop_lost AND \
             frameshift_variant. got {:?}",
            result.consequences
        );
    }
}

// Splice subterm overcall on small variants at an exon boundary

/// A small donor-side insertion must not overcall splice subterms.
///
/// An insertion at the last exonic base or just inside an intron (small
/// variant, span <= 20bp) gets only `splice_donor_variant` from Perl; the
/// `add_splice_for_overlapping_introns_insertion` pass must not extend the
/// variant span by `alt_len` and trigger `splice_donor_5th_base_variant`,
/// `splice_donor_region_variant` or `intron_variant` at positions the
/// insertion's genomic anchor does not reach.
///
/// Shape: insertion at the last exonic base of exon 1 (pos 25_000_299),
/// inserting 11bp. Perl-correct: splice_donor_variant only (the insertion
/// sits at the donor, not at the 5th base or in the core intron).
#[test]
fn concordance_small_donor_insertion_no_subterm_overcall() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // Exon 1 ends at 25_000_299. Intron 1 = 25_000_300..25_001_999.
    // Insertion between 25_000_299 (last exonic base) and 25_000_300 (donor +1).
    // VEP convention: ref="-", start = end + 1.
    let variant = InputVariant::new(
        "21".into(),
        25_000_300, // start
        25_000_299, // end (insertion convention: end < start)
        b"-".to_vec(),
        b"CTTACTTCCCG".to_vec(), // 11bp insertion
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    // No splice_donor_5th_base_variant: the 5th base position is 25_000_304,
    // the insertion is at 25_000_299/25_000_300 which is donor +0/+1.
    assert!(
        !tc.consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant),
        "Small insertion at donor side should NOT emit \
         splice_donor_5th_base_variant; got: {:?}",
        tc.consequences
    );
    // No splice_donor_region_variant: the donor region window is +2..+5, the
    // insertion is at +0/+1.
    assert!(
        !tc.consequences
            .contains(&Consequence::SpliceDonorRegionVariant),
        "Small insertion at donor side should NOT emit \
         splice_donor_region_variant; got: {:?}",
        tc.consequences
    );
    // No intron_variant: the insertion is at the boundary, not in the core
    // intron region intron_start+2..intron_end-2.
    assert!(
        !tc.consequences.contains(&Consequence::IntronVariant),
        "Small insertion at donor side should NOT emit intron_variant; \
         got: {:?}",
        tc.consequences
    );
}

// The alt-CDS ATG scanner runs before the codon-window StartLost check

/// A 1bp deletion at the third base of ATG must emit
/// `start_retained_variant`, not `start_lost`, when the next codon starts
/// with G (the alt CDS still has ATG at codon-aligned position 0).
///
/// Perl emits `frameshift_variant,start_retained_variant` for codons
/// `atG/at` -> `M/X`: it scans the alt CDS for a codon-aligned ATG, finds one
/// when the next codon starts with G, and emits start_retained. The alt-CDS
/// ATG scanner therefore runs before the codon-window peptide check
/// (alt_starts_m=false because the truncated codon translates to X), which
/// would otherwise push StartLost and gate the scanner off.
///
/// `make_test_transcript()` builds a CDS starting `ATGGCT...` so deleting
/// the third base (the G in ATG) produces alt CDS `ATGCT...` -> ATG at
/// codon-aligned position 0 -> Perl emits start_retained_variant.
#[test]
fn concordance_inframe_deletion_at_start_retained_not_lost() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // CDS pos 3 (third base of ATG) on the test transcript maps to genomic
    // position 25_000_052 (mapper pair 1: cdna 1..300 = genomic 25_000_000..
    // 25_000_299, with cdna_coding_start = 51, so CDS 1 = cdna 51 = genomic
    // 25_000_050; CDS 3 = cdna 53 = genomic 25_000_052).
    // 1bp deletion of the G: ref="G", alt="-", start=end=25_000_052.
    let variant = InputVariant::new(
        "21".into(),
        25_000_052,
        25_000_052,
        b"G".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    // After deleting cds[2] (the G), alt CDS is "ATGCTGGAAAATTCGAT..." which
    // has ATG at codon-aligned position 0. Perl emits start_retained_variant.
    assert!(
        tc.consequences.contains(&Consequence::StartRetainedVariant),
        "1bp deletion at CDS pos 3 with codon-2-starts-with-G should \
         emit start_retained_variant; got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StartLost),
        "1bp deletion at CDS pos 3 with codon-2-starts-with-G should \
         NOT emit start_lost (alt CDS retains ATG at codon-aligned position 0); \
         got: {:?}",
        tc.consequences
    );
}

/// A coding-only inframe deletion at the start codon must co-emit
/// `inframe_deletion` and `start_lost`.
///
/// Perl emits `start_lost,inframe_deletion` for a reverse-strand deletion
/// with `cds=2-4, codons=aTGAca/aca`. On a reverse-strand transcript
/// `apply_position` is called with `variant.start` (lower genomic), which
/// maps to the higher cds position, so `cds_pos = 4` and a `*cds_pos <= 3`
/// gate alone fails; the end anchor's `cds_pos = 2` never reaches the block
/// because `populate_fields = false` short-circuits first. The post-frameshift
/// StartLost block in `consequences.rs` therefore also accepts variants whose
/// CDS bounds (from `map_genomic_span_to_cds_bounds`) have
/// `min(cds_start, cds_end) <= 3`.
///
/// This forward-strand test pins the simple coding-deletion-of-ATG path.
///
/// Expected (Perl): `inframe_deletion + start_lost`.
#[test]
fn concordance_inframe_deletion_at_start_codon_emits_start_lost() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();

    // 3bp inframe deletion of the ATG start codon: cds 1-3 -> deleted.
    // CDS starts at genomic 25_000_050 (cdna 51).
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_052,
        b"ATG".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::InframeDeletion),
        "3bp inframe deletion at start codon should emit \
         inframe_deletion, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StartLost),
        "3bp coding-only inframe deletion of ATG start codon must \
         co-emit start_lost (Perl: inframe_deletion + start_lost). \
         got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StartRetainedVariant),
        "Alt CDS does not contain a codon-aligned ATG, so \
         start_retained must NOT fire; got: {:?}",
        tc.consequences
    );
}

/// A coding-only deletion at the stop codon must co-emit
/// `frameshift_variant` and `stop_lost`.
///
/// Perl emits `frameshift_variant,stop_lost` for a reverse-strand deletion
/// with `cds=1415-1416, codons=tAA/t`. On a reverse-strand transcript whose
/// deletion overlaps the last CDS positions, the start anchor `*cds_pos` is
/// the higher CDS index of the span, `cds_len`; with `ref_len >= 2`,
/// `compute_full_peptide_alleles_impl` evaluates
/// `idx + nominal_ref_len > cds.len()` -> true and returns None, so the
/// post-frameshift StopLost block in `consequences.rs` feeds
/// `compute_peptide_alleles` the 5'-most CDS index of the span
/// (`analysis_cds_pos_5prime`) instead.
///
/// This forward-strand test pins a 1bp deletion of the last CDS base of the
/// test transcript's terminal stop codon.
///
/// Expected (Perl): `frameshift_variant + stop_lost`.
#[test]
fn concordance_deletion_at_stop_codon_emits_stop_lost() {
    // Use the test-transcript variant ending with a real TAA stop codon at
    // CDS positions 847-849 (genomic 25_004_296-25_004_298 forward strand).
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();

    // 1bp deletion of the third nucleotide of the TAA stop (cds 849, genomic
    // 25_004_298). ref = "A", alt = "-". Deletion frameshifts and removes
    // the terminal stop -> Perl: frameshift_variant + stop_lost.
    let variant = InputVariant::new(
        "21".into(),
        25_004_298,
        25_004_298,
        b"A".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "1bp deletion at stop codon must emit frameshift_variant, \
         got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StopLost),
        "1bp coding-only deletion at the last CDS base (stop \
         codon) must co-emit stop_lost (Perl: frameshift_variant + \
         stop_lost). got: {:?}",
        tc.consequences
    );
}

/// A 2bp deletion at the stop codon of a reverse-strand transcript must co-emit
/// `frameshift_variant + stop_lost`.
///
/// On a reverse-strand transcript `variant.start` (the lower genomic coordinate)
/// maps to the higher cdna/cds position, so for a stop at `cds_len-2..cds_len` the
/// start anchor lands at `cds_pos = cds_len`. Feeding that anchor to
/// `compute_peptide_alleles` computes `idx = cds_len - 1`, then
/// `idx + ref_len = cds_len + 1 > cds.len()`, the non-clamped implementation
/// returns None, and the post-frameshift StopLost block no-ops without pushing
/// `stop_lost`.
///
/// The anchor resolves it, not a clamped fallback: `analysis_cds_pos_5prime`
/// feeds the 5'-most CDS index of the span, so the peptide path returns defined
/// peptides (ref `*`, alt non-`*`). Perl fires `stop_lost` here through the
/// peptide path, not `_ins_del_stop_altered`. A clamped-peptide fallback in the
/// block instead adds false `stop_lost`, because the block carries Perl
/// invariants a clamp does not preserve.
#[test]
fn concordance_reverse_strand_deletion_at_stop_emits_stop_lost() {
    let mut tx = make_test_transcript();
    // Patch CDS to end with a real TAA stop codon at cds 847-849 (cdna 897-899).
    if let Some(ref mut vefc) = tx.vefc {
        let mut cds = String::with_capacity(849);
        cds.push_str("ATG");
        while cds.len() < 846 {
            cds.push_str("GCT");
        }
        cds.truncate(846);
        cds.push_str("TAA");
        vefc.translateable_seq = Some(cds);
        if let Some(ref mut mapper) = vefc.mapper {
            mapper.cdna_coding_end = 899;
            // Flip mapper orientation to reverse.
            for pair in &mut mapper.exon_coord_mapper.pairs {
                pair.ori = -1;
            }
        }
    }
    tx.cdna_coding_end = Some(899);
    if let Some(ref mut translation) = tx.translation {
        translation.end = 299;
    }
    tx.strand = Strand::Reverse;

    let config = EffectsConfig::default();

    // On reverse strand, cdna 1 maps to the highest genomic position in
    // exon 1. Mapper pair 3 covers cdna 601-2601 -> genomic 25_004_000-25_006_000
    // with ori=-1. So cdna 897 -> genomic 25_006_000 - (897-601) = 25_005_704;
    // cdna 898 -> 25_005_703; cdna 899 -> 25_005_702.
    //
    // Delete cdna 898-899 (the AA of the TAA stop codon, cds 848-849).
    // On the forward genomic strand the ref is the reverse-complement of cdna
    // "AA" = "TT" at genomic 25_005_702-25_005_703.
    let variant = InputVariant::new(
        "21".into(),
        25_005_702,
        25_005_703,
        b"TT".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::FrameshiftVariant),
        "Reverse strand: 2bp deletion at stop codon must emit \
         frameshift_variant, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StopLost),
        "Reverse strand: 2bp coding-only deletion at stop codon \
         on reverse-strand transcript must co-emit stop_lost: the StopLost block \
         anchors on the 5'-most CDS index of the span, because the span END is \
         the 3'-most index there and makes compute_peptide_alleles return None. \
         got: {:?}",
        tc.consequences
    );
}

/// A 3bp in-frame deletion fully overlapping the ATG start codon on a
/// reverse-strand transcript must co-emit `inframe_deletion + start_lost`.
///
/// Two coordinates have to agree.
/// The StartLost gate reads `start_5prime_cds_pos <= 3`,
/// the 5'-most CDS index of the span, because on reverse strand the start anchor
/// (`*cds_pos`) is the span's 3'-most index: for a 3bp deletion at cds 1-3 it
/// reads 6, so a `*cds_pos <= 3` gate never opens. The inner
/// `alt_has_codon_zero_atg` scanner derives `idx` from `start_5prime_cds_pos - 1`
/// for the same reason, since building the alt CDS from `*cds_pos - 1` excises the
/// wrong bases, preserves the original ATG prefix, and makes the scanner emit
/// `start_retained_variant` instead. Widening the gate alone converts the miss
/// into that wrong emission and adds false positives on the 5'UTR start-codon
/// path.
///
/// Perl emits `start_lost` here via the peptide path. The mid-CDS control below
/// must stay excluded, which shows the gate discriminates on position rather
/// than on peptide content.
#[test]
fn concordance_reverse_strand_inframe_deletion_at_atg_emits_start_lost() {
    // Build the standard test transcript and flip it to reverse strand.
    // The CDS sequence (translateable_seq) is in 5'→3' transcript direction,
    // so it always begins with "ATG" regardless of strand. Only the mapper
    // orientation and `tx.strand` need to flip.
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut mapper) = vefc.mapper {
            for pair in &mut mapper.exon_coord_mapper.pairs {
                pair.ori = -1;
            }
        }
    }
    tx.strand = Strand::Reverse;

    let config = EffectsConfig::default();

    // On reverse strand with mapper pair 1 (from_start=1, to_start=25_000_000,
    // to_end=25_000_299, ori=-1), cdna position N maps to genomic position
    // `to_end - (N - from_start)`. cds 1 = cdna 51 -> genomic 25_000_249;
    // cds 2 = cdna 52 -> 25_000_248; cds 3 = cdna 53 -> 25_000_247.
    //
    // Deleting cds 1-3 (the ATG codon): on the forward genomic strand the
    // deletion covers 25_000_247..=25_000_249 in genomic 5'→3' order. The
    // forward-genome REF bases there are the reverse-complement of cdna
    // "ATG" = "CAT" (genomic positions 247=C, 248=A, 249=T).
    let variant = InputVariant::new(
        "21".into(),
        25_000_247,
        25_000_249,
        b"CAT".to_vec(),
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        tc.consequences.contains(&Consequence::InframeDeletion),
        "Reverse strand: 3bp coding-only deletion of ATG on a \
         reverse-strand transcript must emit inframe_deletion, got: {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences.contains(&Consequence::StartLost),
        "Reverse strand: 3bp coding-only deletion of ATG on a \
         reverse-strand transcript must co-emit start_lost (Perl: \
         inframe_deletion + start_lost): the StartLost gate accepts \
         `min(cds_start, cds_end) <= 3`, the 5'-most CDS index of the span, \
         because `*cds_pos` is the END coordinate on reverse strand. got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StartRetainedVariant),
        "Reverse strand: alt CDS no longer contains a \
         codon-aligned ATG (the ATG is fully deleted), so \
         start_retained_variant must NOT fire; got: {:?}",
        tc.consequences
    );
}

/// control (reverse-strand mid-CDS): a 3bp in-frame deletion on a
/// reverse-strand transcript that does not overlap the start codon must emit
/// `inframe_deletion` only, no `start_lost`.
///
/// On `10:88836354 GTCA>G`, three transcripts get the identical codon change
/// `aTGAca/aca` (AA `MT/T`), but Perl emits `start_lost` only on the one where
/// the deletion hits CDS 2-4; the mid-CDS transcripts get `inframe_deletion`
/// alone. The StartLost gate keys on `start_5prime_cds_pos <= 3`, so its
/// 5'-most-index discriminator does not overcall on mid-CDS deletions with an
/// otherwise start-codon-shaped peptide window; without the `<= 3` bound every
/// mid-CDS inframe deletion on a reverse-strand transcript would gain
/// `start_lost`.
#[test]
fn concordance_reverse_strand_midcds_deletion_no_start_lost() {
    // Reverse-strand transcript; delete cds 100-102 (mid-CDS, far from ATG).
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        if let Some(ref mut mapper) = vefc.mapper {
            for pair in &mut mapper.exon_coord_mapper.pairs {
                pair.ori = -1;
            }
        }
    }
    tx.strand = Strand::Reverse;

    let config = EffectsConfig::default();

    // mapper pair 1: from_start=1, to_start=25_000_000, to_end=25_000_299,
    // ori=-1. cdna N -> genomic 25_000_299 - (N - 1). CDS starts at cdna 51,
    // so cds 100 = cdna 150 -> genomic 25_000_150; cds 101 -> 25_000_149;
    // cds 102 -> 25_000_148. Delete cds 100-102 (genomic 25_000_148..=150).
    let variant = InputVariant::new(
        "21".into(),
        25_000_148,
        25_000_150,
        b"GCT".to_vec(), // forward-genome bases (rc of cds window); exact bytes don't gate start_lost
        b"-".to_vec(),
    );

    let result = calculate_consequences(&variant, &tx, &config);
    assert!(result.is_some(), "Should produce consequences");
    let tc = result.unwrap();

    assert!(
        !tc.consequences.contains(&Consequence::StartLost),
        "Control: a mid-CDS (cds 100-102) reverse-strand inframe \
         deletion must NOT emit start_lost: the start codon is untouched. The \
         5'-most-index gate (start_5prime_cds_pos <= 3) must exclude it; \
         got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::StartRetainedVariant),
        "Control: mid-CDS deletion must not touch the start codon, \
         so start_retained must NOT fire; got: {:?}",
        tc.consequences
    );
}

/// `incomplete_terminal_codon_variant` on a span that leaves the CDS.
///
/// Perl's `VariationEffect::partial_codon` serves two roles: it
/// suppresses missense/inframe_deletion/stop_retained/frameshift, and it emits
/// `incomplete_terminal_codon_variant` in its own right, including on the
/// span-leaves-CDS path, which Perl reports as `CDS_position = N-?`.
///
/// The stock fixture's CDS is 850 bp and 850 % 3 == 1, so CDS position 850 is a
/// 1-base partial terminal codon; `coding_region_end` is 25_004_299 and exon 3
/// runs to 25_006_000, giving real 3'UTR for the span to run into.
#[test]
fn concordance_cds_to_three_prime_utr_span_partial_terminal_codon() {
    let tx = make_test_transcript_with_flags(&["cds_end_NF"]);
    let variant = InputVariant::new(
        "21".to_string(),
        25_004_299,
        25_004_302,
        b"GCTG".to_vec(),
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default())
        .expect("deletion spanning the partial terminal codon into 3'UTR must annotate");
    assert!(
        tc.consequences
            .contains(&Consequence::IncompleteTerminalCodonVariant),
        "a span that starts in a partial terminal codon and exits the CDS must emit \
         incomplete_terminal_codon_variant (Perl partial_codon, CDS_position N-?); \
         got: {:?}",
        tc.consequences
    );
    assert!(
        !tc.consequences.contains(&Consequence::FrameshiftVariant),
        "Perl's frameshift returns 0 when partial_codon holds; got: {:?}",
        tc.consequences
    );
}

/// Guard: a leading Gap must suppress the term.
///
/// Perl's `partial_codon` returns early (`return 0 unless defined
/// $bvfo->translation_start`), and `translation_start` is undef when the 5'-most
/// mapper segment is a Gap. Perl renders that shape as `CDS_position = ?-N`. This
/// pins the `segments.first()` guard: without it the term fires on exactly this
/// geometry.
#[test]
fn concordance_leading_gap_no_partial_terminal_codon_term() {
    let tx = make_test_transcript_with_flags(&["cds_end_NF"]);
    // Start inside intron 2 (exon 2 ends 25_002_299, exon 3 starts 25_004_000) so
    // the 5'-most segment is a Gap, and run into the partial terminal codon.
    let variant = InputVariant::new(
        "21".to_string(),
        25_003_000,
        25_004_299,
        vec![b'A'; 1300],
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default())
        .expect("intron-to-CDS deletion must annotate");
    assert!(
        !tc.consequences
            .contains(&Consequence::IncompleteTerminalCodonVariant),
        "a leading Gap leaves Perl's translation_start undef, so partial_codon \
         short-circuits and the term must NOT fire; got: {:?}",
        tc.consequences
    );
}

// Short-intron splice_region spill into the flanking exon

/// Shrink intron 1 to `len` bases, keeping exon 1 fixed and pulling exon 2 in to
/// meet it, so the transcript stays self-consistent (mapper, sorted_exons, and
/// both intron lists agree).
///
/// Perl marks any intron with `abs(end - start) <= 12` as `_frameshift`
/// (`BaseTranscriptVariation::_create_intron_trees`), which is the regime this
/// exercises.
fn make_transcript_with_short_first_intron(len: u64) -> Transcript {
    let mut tx = make_test_transcript();
    let intron_start = 25_000_300;
    let intron_end = intron_start + len - 1;
    let exon2_start = intron_end + 1;
    let exon2_end = exon2_start + 299;

    tx.exons[1].start = exon2_start;
    tx.exons[1].end = exon2_end;
    tx.introns[0].start = intron_start;
    tx.introns[0].end = intron_end;

    if let Some(vefc) = tx.vefc.as_mut() {
        vefc.sorted_exons[1].start = exon2_start;
        vefc.sorted_exons[1].end = exon2_end;
        vefc.introns[0].start = intron_start;
        vefc.introns[0].end = intron_end;
        if let Some(mapper) = vefc.mapper.as_mut() {
            let mut pairs = mapper.exon_coord_mapper.pairs.clone();
            pairs[1].to_start = exon2_start;
            pairs[1].to_end = exon2_end;
            mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
        }
    }
    tx
}

/// An exonic SNV 4-7bp past a sub-8bp intron must get `splice_region_variant`.
///
/// Perl's `VariationEffect::_intron_overlap` tests
/// `overlap($vf_start, $vf_end, $intron_start + 2, $intron_start + 7)` with no
/// clamp to the intron. For a 1bp intron that window is `[start+2, start+7]` while
/// the intron occupies only `start`, so bases 1 through 7 of the downstream exon
/// fall inside it and Perl emits the term for a variant that never touches the
/// intron at all.
///
/// The variant here sits 7bp past a 1bp intron: outside the intron, outside the
/// 3bp acceptor window (`[intron_end+1, intron_end+3]`), inside the spill. Offset
/// +7 specifically, because Perl's own precedence claims the near half of the
/// spill: `VariationEffect::splice_region` returns 0 when
/// `splice_donor_region_variant` holds, and that window is `[intron_start+2,
/// intron_start+5]`, with `splice_donor_5th_base_variant` taking `intron_start+4`.
/// So only offsets +6 and +7 surface as `splice_region_variant`.
#[test]
fn concordance_exonic_snv_near_one_bp_intron_gets_splice_region() {
    let tx = make_transcript_with_short_first_intron(1);
    // Intron 1 is exactly 25_000_300. Exon 2 starts 25_000_301; 25_000_307 is
    // intron_start + 7, the far edge of [start+2, start+7].
    let variant = InputVariant::new(
        "21".to_string(),
        25_000_307,
        25_000_307,
        b"G".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default())
        .expect("exonic SNV near a short intron must annotate");
    assert!(
        tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "Perl's unclamped [intron_start+2, intron_start+7] window spills out of a \
         1bp intron into the next exon, so an exonic SNV at intron_start + 7 must \
         emit splice_region_variant; got: {:?}",
        tc.consequences
    );
}

/// Guard: the same offset past a long intron must not get the term.
///
/// This is the discriminator that makes the test above non-vacuous. With an 8bp
/// or longer intron, `[intron_start+2, intron_start+7]` is fully contained in the
/// intron, nothing spills, and an exonic variant 5bp past the boundary is 4bp
/// outside the 3bp acceptor window. Without the `intron_len <= 7` guard, this
/// test fires.
#[test]
fn concordance_exonic_snv_near_eight_bp_intron_no_splice_region() {
    let tx = make_transcript_with_short_first_intron(8);
    // Intron 1 spans 25_000_300..=25_000_307; exon 2 starts 25_000_308.
    // 25_000_312 is 5bp into exon 2 and intron_start + 12, outside every window.
    let variant = InputVariant::new(
        "21".to_string(),
        25_000_312,
        25_000_312,
        b"G".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default())
        .expect("exonic SNV near an 8bp intron must annotate");
    assert!(
        !tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "an 8bp intron fully contains [intron_start+2, intron_start+7], so nothing \
         spills and an exonic SNV 5bp past the boundary must NOT emit \
         splice_region_variant; got: {:?}",
        tc.consequences
    );
}

/// Guard: an insertion anchored near a short intron must not get the term.
///
/// Perl reaches insertions through `_intron_overlap`'s dedicated `$insertion`
/// clause, four exact-equality tests against the intron
/// edges, not through the spill windows. An insertion's VEP anchor pair
/// (`ref = "-"`, `end = start - 1`) straddles a boundary by construction, so
/// admitting insertions to the spill produces false positives at insertion
/// anchors. This pins the `is_insertion` carve-out.
#[test]
fn concordance_insertion_near_short_intron_no_splice_region_spill() {
    let tx = make_transcript_with_short_first_intron(1);
    // Insert between 25_000_304 and 25_000_305: the anchor pair brackets
    // intron_start + 5, inside the spill window, but Perl does not apply it here.
    let variant = InputVariant::new(
        "21".to_string(),
        25_000_305,
        25_000_304,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default())
        .expect("insertion near a short intron must annotate");
    assert!(
        !tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "insertions go through Perl's exact-equality $insertion clause, not the \
         spill windows, so an insertion anchored in the spill must NOT emit \
         splice_region_variant; got: {:?}",
        tc.consequences
    );
}

/// Guard: a variant inside the short intron keeps its intronic classification.
///
/// The spill windows are for variants wholly outside the intron. A variant that
/// overlaps the intron is the intronic path's business, and for a 1bp intron that
/// single base is both `intron_start` and `intron_start + 1`, so Perl's
/// `start_splice_site` fires and `splice_region` is explicitly suppressed
/// (the `unless` clause of `BaseTranscriptVariationAllele::_intron_effects`).
/// This pins the `outside_intron` guard.
#[test]
fn concordance_variant_inside_short_intron_keeps_intronic_terms() {
    let tx = make_transcript_with_short_first_intron(1);
    let variant = InputVariant::new(
        "21".to_string(),
        25_000_300,
        25_000_300,
        b"G".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default())
        .expect("variant inside a short intron must annotate");
    assert!(
        !tc.consequences.contains(&Consequence::SpliceRegionVariant),
        "Perl suppresses splice_region when start_splice_site fires, which it does \
         for the sole base of a 1bp intron; got: {:?}",
        tc.consequences
    );
}

/// Build a forward-strand transcript whose CDS ends only 3 bp before the
/// transcript's 3' boundary, so a deletion anchored in the CDS can span the
/// stop codon, cross the short 3'UTR, and extend past the transcript end.
///
/// That geometry is what produces `TranscriptPosition::Downstream` for the
/// paired endpoint, which is the shape excluded from `is_stop_lost`.
/// `make_test_transcript()` cannot express it: its CDS ends at 25_004_299 with
/// a 1,701 bp 3'UTR, so a deletion long enough to leave the transcript is
/// nowhere near the stop codon.
fn make_transcript_with_cds_at_three_prime_edge() -> Transcript {
    let mut tx = make_test_transcript();
    // Shrink the last exon and the transcript so the CDS end (25_004_299) sits
    // 3 bp inside the 3' boundary. cDNA/CDS coords are unchanged: the CDS is
    // still cDNA 51..900 and the deleted tail was pure 3'UTR.
    let new_end = 25_004_302u64;
    tx.end = new_end;
    if let Some(last) = tx.exons.last_mut() {
        last.end = new_end;
    }
    if let Some(vefc) = tx.vefc.as_mut() {
        if let Some(last) = vefc.sorted_exons.last_mut() {
            last.end = new_end;
        }
        // The mapper is the exon table Perl's `TranscriptMapper` reads; it has to
        // shrink with the exon or the span maps as still inside the transcript.
        if let Some(mapper) = vefc.mapper.as_mut() {
            let mut pairs = mapper.exon_coord_mapper.pairs.clone();
            if let Some(last) = pairs.last_mut() {
                last.from_end -= last.to_end - new_end;
                last.to_end = new_end;
            }
            mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
        }
    }
    tx
}

/// A deletion that starts in the CDS and extends past the transcript's
/// 3' boundary must not emit `stop_lost`.
///
/// Perl reaches `stop_lost` for such an allele only through
/// `VariationEffect::_ins_del_stop_altered`, because
/// `_get_peptide_alleles` is undef for a span that leaves the CDS so
/// `stop_lost` delegates to it.
/// `_ins_del_stop_altered` opens with
/// `return 0 unless _overlaps_stop_codon(@_)`, and `_overlaps_stop_codon`
/// has `return 0 unless $cdna_start && $cdna_end`. For a span that
/// runs off the transcript, `cdna_coords`' last entry is a
/// `Bio::EnsEMBL::Mapper::Gap`, so `BaseTranscriptVariation::cdna_end` is
/// undef and that guard closes. With every stop predicate at 0,
/// `coding_unknown` fires and Perl emits
/// `3_prime_UTR_variant,coding_sequence_variant`.
///
/// The mapper snaps the out-of-transcript endpoint to the nearest exon
/// boundary, so `cdna_hi` is defined and a gate that admitted
/// `TranscriptPosition::Downstream` alongside `ThreePrimeUtr` would open. This
/// is the 3' mirror of the `extends_upstream_of_transcript` exclusion on
/// `is_start_lost`.
#[test]
fn concordance_cds_deletion_past_transcript_end_no_stop_lost() {
    let tx = make_transcript_with_cds_at_three_prime_edge();
    let config = EffectsConfig::default();

    // Delete from inside the CDS (25_004_290, well before the stop codon at
    // 25_004_297..25_004_299) out to 25_004_400, 98 bp beyond the transcript
    // end at 25_004_302. Paired endpoint therefore resolves to Downstream.
    let ref_len = (25_004_400u64 - 25_004_290u64 + 1) as usize;
    let variant = InputVariant::new(
        "21".to_string(),
        25_004_290,
        25_004_400,
        vec![b'A'; ref_len],
        b"-".to_vec(),
    );

    let tc = calculate_consequences(&variant, &tx, &config)
        .expect("CDS→past-transcript-end deletion must annotate");

    assert!(
        !tc.consequences.contains(&Consequence::StopLost),
        "A deletion extending past the transcript 3' boundary leaves \
         Perl's cdna_end undef, closing _overlaps_stop_codon and hence \
         _ins_del_stop_altered, so stop_lost must NOT be emitted; got {:?}",
        tc.consequences
    );
    assert!(
        tc.consequences
            .contains(&Consequence::CodingSequenceVariant),
        "With every stop predicate at 0, Perl's coding_unknown emits \
         coding_sequence_variant; got {:?}",
        tc.consequences
    );
}

/// Correct-side guard: the same transcript, with the deletion stopping
/// inside the 3'UTR instead of running off the end, must still emit `stop_lost`.
///
/// Here `cdna_end` is defined (the endpoint maps to a real 3'UTR cDNA
/// position), so `_overlaps_stop_codon` returns 1, `_ins_del_stop_altered`
/// evaluates the re-translated codon, and Perl emits `stop_lost`. This pins the
/// exclusion to the Downstream shape only: a `ThreePrimeUtr`-paired deletion
/// over the stop codon must be untouched.
#[test]
fn concordance_cds_deletion_into_three_prime_utr_keeps_stop_lost() {
    let tx = make_transcript_with_cds_at_three_prime_edge();
    let config = EffectsConfig::default();

    // Delete 25_004_290..25_004_301: covers the stop codon and ends inside the
    // transcript (end 25_004_302), so the paired endpoint is ThreePrimeUtr.
    let ref_len = (25_004_301u64 - 25_004_290u64 + 1) as usize;
    let variant = InputVariant::new(
        "21".to_string(),
        25_004_290,
        25_004_301,
        vec![b'A'; ref_len],
        b"-".to_vec(),
    );

    let tc =
        calculate_consequences(&variant, &tx, &config).expect("CDS→3'UTR deletion must annotate");

    assert!(
        tc.consequences.contains(&Consequence::StopLost),
        "The gate must NOT suppress the ThreePrimeUtr-paired shape: cdna_end is \
         defined there, so Perl's _overlaps_stop_codon opens and stop_lost \
         fires; got {:?}",
        tc.consequences
    );
}

/// An insertion with one intronic flank must not collapse both CDS
/// endpoints onto a single position, because that off-by-one flips the
/// boundary-insertion predicate and swaps `frameshift_variant` against
/// `inframe_insertion,stop_retained_variant`.
///
/// `map_genomic_span_to_cds_bounds` maps each flank through
/// `genomic_to_cdna_or_nearest_boundary`, which snaps an intronic flank onto the
/// nearest exon-boundary cDNA position. When the other flank is the exonic base
/// at that same boundary, both endpoints resolve to the same cDNA position, so
/// `cds_a == cds_b` and `hi = max(cds_a, cds_b)` under-reports Perl's
/// `cds_start` by exactly one. Perl's `TranscriptMapper::genomic2cds` (ensembl
/// core) never collapses: it returns the real adjacent pair, so Perl's
/// `BaseTranscriptVariation::cds_start` / `cds_end` differ by one and its
/// `translation_start`/`translation_end` come from an independent `genomic2pep`
/// mapping (`BaseTranscriptVariation::translation_start` ->
/// `translation_coords` -> `TranscriptMapper::genomic2pep`,
/// `pep = int((cds + shift + 2) / 3)`).
///
/// Because `is_boundary_insertion` (coding.rs) is
/// `cds_start > cds_end && translation_start != translation_end`, the collapsed
/// value shifts which codon window is used. Perl emits a ranged
/// `Protein_position` for an adjacent insertion pair iff `cds_end % 3 == 0`;
/// the collapsed bounds satisfy that predicate iff `cds_end % 3 == 1`, so they
/// are correct on one phase only.
///
/// This test pins the CDS/peptide bounds for the snapped geometry. The
/// `make_test_transcript` exon 1 spans 25_000_000..=25_000_299 (cDNA 1..=300)
/// with CDS starting at cDNA 51, so the last exonic base is CDS
/// 300 - 51 + 1 = 250 and the first intronic base is 25_000_300. An insertion
/// there (VEP convention `start == end + 1`) must yield cds_start 251 /
/// cds_end 250, not the collapsed 250 / 249.
///
/// Adding 1 to the snapped `cds_start` matches Perl's `cds_start` but is not by
/// itself a parity model, so this crate keeps the snapped bounds and this test is
/// `#[ignore]`d. Un-snapping the bounds has two effects, not one:
///
///   1. The boundary predicate `is_boundary_insertion = (tl_start != tl_end)`
///      flips at phase0 and phase1 only.
///   2. `splice_start = cds_start - 1` (both splice-window helpers in `coding.rs`) shifts the
///      alt-CDS reading frame at every phase, regardless of that predicate.
///
/// Because the window extents are chosen by a branch on `is_boundary_insertion`,
/// phase0/phase1 rows get a new window and a new frame, while phase2 is the
/// unique phase where extents are unchanged and the frame shift is the only
/// change. A narrow window can also find a `stop_gained` stop that the
/// correctly-widened window loses, a window-content effect the `+1` does not
/// model.
#[test]
#[ignore = "pins Perl's uncollapsed cds_start for a snapped boundary insertion, which this crate does not reproduce: the +1 alone is not the complete parity model, since splice_start = cds_start - 1 shifts the alt-CDS frame at every phase, not just where the boundary predicate flips"]
fn concordance_snapped_insertion_cds_bounds_not_collapsed() {
    let tx = make_test_transcript();

    // Insertion at the exon1/intron1 boundary: start = end + 1.
    // span_start 25_000_300 is intronic (snaps to exon 1 end, cDNA 300),
    // span_end 25_000_299 is the exonic last base of exon 1 (cDNA 300).
    let bounds = crate::mapper::map_genomic_span_to_cds_bounds(25_000_300, 25_000_299, &tx)
        .expect("snapped boundary insertion should still produce CDS bounds");

    // Perl parity: the adjacent pair is (250, 251), so under insertion semantics
    // cds_start = 251 > cds_end = 250 (vf_nt_len = 0); collapsing both flanks onto
    // CDS 250 gives (250, 249).
    assert_eq!(
        bounds.cds_start, 251,
        "snapped insertion must not collapse cds_start (a collapsed cds_start would read 250)"
    );
    assert_eq!(
        bounds.cds_end, 250,
        "snapped insertion cds_end must be cds_start - 1"
    );

    // cds_end 250 has 250 % 3 == 1, so Perl reports a single Protein_position
    // and this is not a boundary insertion: pep(251) == pep(250) == 84.
    assert_eq!(bounds.translation_start, 84, "tl_start = (251-1)/3+1 = 84");
    assert_eq!(bounds.translation_end, 84, "tl_end = (250-1)/3+1 = 84");
    assert_eq!(
        bounds.translation_start, bounds.translation_end,
        "cds_end % 3 == 1 => Perl emits a single Protein_position, so the \
         boundary-insertion predicate must be FALSE. A collapsed pair would give \
         tl_start 84 / tl_end 83, making it spuriously TRUE and swapping \
         frameshift_variant for inframe_insertion,stop_retained_variant."
    );
}

// Perl codon / peptide window (`coding::perl_coding_terms`)

/// Move the CDS start of `make_test_transcript()` from cDNA 51 to cDNA 49 so
/// the last base of exon 1 (cDNA 300) is CDS 252, a codon boundary. The
/// translateable sequence keeps its 850 bp.
fn with_cds_start_at_cdna_49(mut tx: Transcript) -> Transcript {
    if let Some(vefc) = tx.vefc.as_mut() {
        if let Some(mapper) = vefc.mapper.as_mut() {
            mapper.cdna_coding_start = 49;
            mapper.cdna_coding_end = 898;
        }
    }
    tx.cdna_coding_start = Some(49);
    tx.cdna_coding_end = Some(898);
    match tx.strand {
        Strand::Forward => {
            tx.coding_region_start = Some(25_000_048);
            tx.coding_region_end = Some(25_004_297);
        }
        Strand::Reverse => {
            tx.coding_region_start = Some(25_001_703);
            tx.coding_region_end = Some(25_005_952);
        }
    }
    tx
}

/// A reverse-strand transcript whose exon blocks descend genomically with
/// cDNA position, the layout a real reverse-strand transcript has:
///   exon 1 (cDNA 1..300)    = 25_005_701..25_006_000
///   exon 2 (cDNA 301..600)  = 25_003_701..25_004_000
///   exon 3 (cDNA 601..2601) = 25_000_000..25_002_000
/// CDS cDNA 51..900, so the coding region is 25_001_701..25_005_950.
fn make_descending_reverse_strand_transcript() -> Transcript {
    let mut tx = make_test_transcript();
    tx.strand = Strand::Reverse;
    let blocks: [(u64, u64); 3] = [
        (25_005_701, 25_006_000),
        (25_003_701, 25_004_000),
        (25_000_000, 25_002_000),
    ];
    for (exon, (s, e)) in tx.exons.iter_mut().zip(blocks) {
        exon.start = s;
        exon.end = e;
    }
    tx.introns = vec![
        Intron {
            start: 25_004_001,
            end: 25_005_700,
            rank: 1,
        },
        Intron {
            start: 25_002_001,
            end: 25_003_700,
            rank: 2,
        },
    ];
    tx.coding_region_start = Some(25_001_701);
    tx.coding_region_end = Some(25_005_950);
    tx.translation_start = Some(25_005_950);
    tx.translation_end = Some(25_001_701);
    if let Some(vefc) = tx.vefc.as_mut() {
        vefc.sorted_exons = tx.exons.clone();
        vefc.introns = tx.introns.clone();
        if let Some(mapper) = vefc.mapper.as_mut() {
            let mut pairs = mapper.exon_coord_mapper.pairs.clone();
            for (pair, (s, e)) in pairs.iter_mut().zip(blocks) {
                pair.to_start = s;
                pair.to_end = e;
                pair.ori = -1;
            }
            mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
        }
    }
    tx
}

/// `perl_span` places an insertion between the last exon base and the intron
/// after that base (`Bio::EnsEMBL::Mapper::map_insert` increments the start of
/// the exonic coordinate): CDS 251 / 250 for the exon 1 / intron 1 boundary of
/// `make_test_transcript()`, not the collapsed 250 / 249 that snapping both
/// flanks onto the boundary produces.
#[test]
fn concordance_perl_span_donor_boundary_insertion_is_after_last_exon_base() {
    let tx = make_test_transcript();
    let span = crate::coding::perl_span(&tx, 25_000_300, 25_000_299)
        .expect("coding transcript with a mapper");
    assert_eq!((span.cds_start, span.cds_end), (Some(251), Some(250)));
    assert_eq!((span.tl_start, span.tl_end), (Some(84), Some(84)));
    assert_eq!((span.cdna_start, span.cdna_end), (Some(301), Some(300)));

    // Acceptor side of intron 1: before the first base of exon 2 (cDNA 301, CDS 251).
    let span = crate::coding::perl_span(&tx, 25_002_000, 25_001_999)
        .expect("coding transcript with a mapper");
    assert_eq!((span.cds_start, span.cds_end), (Some(251), Some(250)));
    assert_eq!((span.cdna_start, span.cdna_end), (Some(301), Some(300)));
}

/// Perl's codon window for a codon-boundary insertion is the inserted bases
/// alone (`codon_len` 0 plus `allele_len`), read from the alternate CDS, and its
/// peptide carries `X` for the partial codon: `-/GTAA` gives `-/VX`.
#[test]
fn concordance_perl_codon_window_boundary_insertion_is_allele_only() {
    let tx = with_cds_start_at_cdna_49(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_299,
        b"-".to_vec(),
        b"GTAA".to_vec(),
    );
    let (ref_codon, alt_codon, ref_pep, alt_pep) =
        crate::coding::perl_codon_peptides(&variant, &tx, 25_000_300, 25_000_299, None)
            .expect("coding allele");
    assert_eq!(ref_codon.as_deref(), Some(b"-".as_slice()));
    assert_eq!(alt_codon.as_deref(), Some(b"GTAA".as_slice()));
    assert_eq!(ref_pep.as_deref(), Some(b"-".as_slice()));
    assert_eq!(alt_pep.as_deref(), Some(b"VX".as_slice()));
}

/// A 4 bp insertion after the last base of a CDS exon, at a codon boundary, is
/// `frameshift_variant` (Perl: `10:18823156-18823157 GTAA ENST00000282343`,
/// `-/VX`). The window is the inserted bases alone, so no stop can be read
/// from the exon's last codon and `stop_retained` has nothing to retain.
#[test]
fn concordance_donor_boundary_frameshift_insertion_is_not_stop_retained() {
    let tx = with_cds_start_at_cdna_49(make_test_transcript());
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_299,
        b"-".to_vec(),
        b"GTAA".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["frameshift_variant", "splice_region_variant"],
    );
}

/// The same boundary on the reverse strand: the inserted bases are read in
/// transcript sense (`GTAA` -> `TTAC`, `LX`) and the verdict is the same.
#[test]
fn concordance_donor_boundary_frameshift_insertion_reverse_strand() {
    let tx = with_cds_start_at_cdna_49(make_descending_reverse_strand_transcript());
    let config = EffectsConfig::default();
    // Exon 1's transcript-3' end is genomic 25_005_701; intron 1 begins at 25_005_700.
    let variant = InputVariant::new(
        "21".into(),
        25_005_701,
        25_005_700,
        b"-".to_vec(),
        b"GTAA".to_vec(),
    );
    let (ref_codon, alt_codon, _, alt_pep) =
        crate::coding::perl_codon_peptides(&variant, &tx, 25_005_701, 25_005_700, None)
            .expect("coding allele");
    assert_eq!(ref_codon.as_deref(), Some(b"-".as_slice()));
    assert_eq!(alt_codon.as_deref(), Some(b"TTAC".as_slice()));
    assert_eq!(alt_pep.as_deref(), Some(b"LX".as_slice()));
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["frameshift_variant", "splice_region_variant"],
    );
}

/// A 6 bp insertion at the same codon boundary is `inframe_insertion` alone:
/// the alt peptide `VG` trivially contains the empty ref peptide, and
/// `ref_eq_alt_sequence` finds no stop (Perl: `16:16173335-16173336 GTAGGA
/// ENST00000345148`).
#[test]
fn concordance_donor_boundary_inframe_insertion_is_not_stop_retained() {
    let tx = with_cds_start_at_cdna_49(make_test_transcript());
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_299,
        b"-".to_vec(),
        b"GTAGGA".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["inframe_insertion", "splice_region_variant"],
    );
}

/// A 3 bp insertion at the codon boundary is `inframe_insertion`, not
/// `protein_altering_variant`: with an empty ref peptide the containment test
/// of `inframe_insertion` passes and that of `protein_altering_variant`
/// suppresses it (Perl: `10:73585593-73585594 ACA ENST00000394934`).
#[test]
fn concordance_donor_boundary_codon_insertion_is_inframe_not_protein_altering() {
    let tx = with_cds_start_at_cdna_49(make_test_transcript());
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_299,
        b"-".to_vec(),
        b"ACA".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["inframe_insertion", "splice_region_variant"],
    );
}

/// An insertion after the last exon base that is not a codon boundary keeps
/// the codon it interrupts: `make_test_transcript()` exon 1 ends at CDS 250,
/// the first base of codon 84, so Perl's window is that codon plus the 4
/// inserted bases. With codon 84 = `GAA` and the insertion `CCTA` the window
/// reads `GCC TAA A` = `A*X`: `frameshift_variant,stop_gained` (Perl:
/// `11:6636096-6636097 ATCA ENST00000299427`). Anchoring the insertion one base
/// early reads `CCT AGA A` = `PRX` and loses the stop.
#[test]
fn concordance_donor_boundary_frameshift_insertion_reads_stop_in_perl_window() {
    let mut tx = make_test_transcript();
    if let Some(seq) = tx.vefc.as_mut().and_then(|v| v.translateable_seq.as_mut()) {
        // Codon 84 = CDS 250..252 (0-based 249..251): the exon 1 last base and the
        // first two bases of exon 2.
        seq.replace_range(249..252, "GAA");
    }
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_299,
        b"-".to_vec(),
        b"CCTA".to_vec(),
    );
    let (ref_codon, alt_codon, ref_pep, alt_pep) =
        crate::coding::perl_codon_peptides(&variant, &tx, 25_000_300, 25_000_299, None)
            .expect("coding allele");
    assert_eq!(ref_codon.as_deref(), Some(b"GAA".as_slice()));
    assert_eq!(alt_codon.as_deref(), Some(b"GCCTAAA".as_slice()));
    assert_eq!(ref_pep.as_deref(), Some(b"E".as_slice()));
    assert_eq!(alt_pep.as_deref(), Some(b"A*X".as_slice()));
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["frameshift_variant", "stop_gained", "splice_region_variant"],
    );
}

/// A 4 bp insertion between the last sense codon and the stop codon: the alt
/// peptide `VX` appended after the whole reference protein leaves the protein
/// unchanged with fewer than 3 residues past its end, so `ref_eq_alt_sequence`
/// holds, `stop_retained` is 1, `frameshift` is therefore 0 and
/// `inframe_insertion` passes on the empty ref peptide (Perl:
/// `12:7031565-7031566 GTGA ENST00000229277`).
#[test]
fn concordance_insertion_before_stop_codon_is_inframe_insertion_stop_retained() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();
    // CDS 846 | 847 (the TAA at 847..849): cDNA 896 | 897 = genomic 25_004_295 | 25_004_296.
    let variant = InputVariant::new(
        "21".into(),
        25_004_296,
        25_004_295,
        b"-".to_vec(),
        b"GTGA".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["inframe_insertion", "stop_retained_variant"],
    );
}

/// A delins replacing the first base of ATG with seven bases is `start_lost`
/// alone: `_ins_del_start_altered` is true and neither inframe predicate
/// rescues it (`SDL` contains `M` at neither end), and `start_lost` then
/// suppresses `protein_altering_variant` (Perl: `12:133226468 TCAGACC
/// ENST00000503265`).
#[test]
fn concordance_start_codon_delins_is_start_lost_without_protein_altering() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_050,
        b"A".to_vec(),
        b"TCAGACC".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["start_lost"]);
}

/// A 1 bp insertion inside the incomplete terminal codon of an 850 bp CDS:
/// `partial_codon` is 1 (so no frameshift), the ref window is empty, the alt
/// window is the single inserted base (`X`), `inframe_insertion` passes on the
/// empty ref and `coding_unknown` fires on the `X` (Perl: `10:50738771-50738772
/// A ENST00000462247`).
#[test]
fn concordance_insertion_in_incomplete_terminal_codon_keeps_inframe_insertion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();
    // Between CDS 849 and 850 (the lone base of the partial codon): cDNA 899 | 900 =
    // genomic 25_004_298 | 25_004_299.
    let variant = InputVariant::new(
        "21".into(),
        25_004_299,
        25_004_298,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "coding_sequence_variant",
            "incomplete_terminal_codon_variant",
            "inframe_insertion",
        ],
    );
}

/// An insertion between the 5' UTR and the first CDS base maps to a CDS gap
/// (`genomic2cds`: the insertion's end lies before `cdna_coding_start`), so no
/// coding predicate fires; Perl reports only the UTR term (Perl:
/// `17:2573570-2573571 A ENST00000451360`).
#[test]
fn concordance_insertion_between_utr_and_start_codon_has_no_coding_term() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_050,
        25_000_049,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["5_prime_UTR_variant"]);
}

/// An insertion between the stop codon's last base and the 3' UTR is likewise
/// a CDS gap (`start > cdna_coding_end`); nothing coding is emitted (Perl:
/// `1:243663044-243663045 TTAT ENST00000336199`).
#[test]
fn concordance_insertion_after_stop_codon_has_no_coding_term() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();
    // CDS 849 | UTR: cDNA 899 | 900 = genomic 25_004_298 | 25_004_299.
    let variant = InputVariant::new(
        "21".into(),
        25_004_299,
        25_004_298,
        b"-".to_vec(),
        b"TTAT".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["3_prime_UTR_variant"]);
}

/// A 3 bp MNV over the last sense codon and the stop codon that keeps a stop at
/// the same peptide index (`GCT TAA` -> `GCA TGA`, `A*` -> `A*`) is
/// `stop_retained_variant`, and `synonymous_variant` is suppressed by it
/// (Perl: `5:256470-256472 GTG/ATA ENST00000509564`).
#[test]
fn concordance_mnv_preserving_stop_is_stop_retained_not_synonymous() {
    let tx = make_transcript_with_terminal_stop();
    let config = EffectsConfig::default();
    // CDS 846..848 = cDNA 896..898 = genomic 25_004_295..25_004_297: `TTA` -> `ATG`.
    let variant = InputVariant::new(
        "21".into(),
        25_004_295,
        25_004_297,
        b"TTA".to_vec(),
        b"ATG".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["stop_retained_variant"]);
}

/// `perl_span` on the reverse strand: `map_insert` swaps the two flanks' roles
/// (`($c2, $c1) = @coords` for strand -1), so an insertion 3' of an exon's last
/// transcript base still lands after it (donor side, exon 1 / intron 1) and an
/// insertion 5' of an exon's first transcript base lands before it (acceptor
/// side, intron 1 / exon 2).
#[test]
fn concordance_perl_span_reverse_strand_boundary_insertions() {
    let tx = make_descending_reverse_strand_transcript();
    // Exon 1 (cDNA 1..300) is 25_005_701..25_006_000; its last transcript base is
    // genomic 25_005_701 = cDNA 300 = CDS 250.
    let donor = crate::coding::perl_span(&tx, 25_005_701, 25_005_700).expect("coding transcript");
    assert_eq!((donor.cdna_start, donor.cdna_end), (Some(301), Some(300)));
    assert_eq!((donor.cds_start, donor.cds_end), (Some(251), Some(250)));
    // Exon 2 (cDNA 301..600) is 25_003_701..25_004_000; its first transcript base
    // is genomic 25_004_000 = cDNA 301 = CDS 251.
    let acceptor =
        crate::coding::perl_span(&tx, 25_004_001, 25_004_000).expect("coding transcript");
    assert_eq!(
        (acceptor.cdna_start, acceptor.cdna_end),
        (Some(301), Some(300))
    );
    assert_eq!(
        (acceptor.cds_start, acceptor.cds_end),
        (Some(251), Some(250))
    );
    // A deletion of the whole of exon 2 maps to [Gap, Coord, Gap] on the CDS
    // side: both `cds_start` and `cds_end` are undef.
    let span = crate::coding::perl_span(&tx, 25_003_700, 25_004_001).expect("coding transcript");
    assert_eq!((span.cds_start, span.cds_end), (None, None));
    assert_eq!(span.cds_coords.len(), 3);
}

// Differing-region semantics of `BaseTranscriptVariationAllele::_intron_effects`
//
// Every window is tested per differing region: for a variant with either allele
// of 1 base or less that is the raw reference span, otherwise the XOR of the two
// alleles, as long as the longer one. The introns tested are the ones the
// memoised interval-tree fetch over the raw variant span returns
// (`BaseTranscriptVariation::_overlapped_introns`, `_overlapped_introns_boundary`,
// `_create_intron_trees`), and `splice_donor_region_variant` and `splice_region`
// (`VariationEffect`) are suppressed by the transcript-wide flags. Each test
// below pins one shape against those predicates.

/// Strip the coding model from a transcript: `processed_transcript`, no
/// translation, no coding coordinates. Exonic verdicts become
/// non_coding_transcript_exon_variant and intronic ones carry
/// non_coding_transcript_variant, which keeps the splice-term assertions free of
/// peptide arithmetic.
fn strip_coding_model(mut tx: Transcript) -> Transcript {
    tx.biotype = "processed_transcript".into();
    tx.translation = None;
    tx.cdna_coding_start = None;
    tx.cdna_coding_end = None;
    tx.coding_region_start = None;
    tx.coding_region_end = None;
    tx.translation_start = None;
    tx.translation_end = None;
    tx.protein_id = None;
    if let Some(vefc) = tx.vefc.as_mut() {
        vefc.translateable_seq = None;
        if let Some(mapper) = vefc.mapper.as_mut() {
            mapper.cdna_coding_start = 0;
            mapper.cdna_coding_end = 0;
        }
    }
    tx
}

fn consequence_of(
    variant: &InputVariant,
    tx: &Transcript,
) -> vep_core::consequence::TranscriptConsequence {
    calculate_consequences(variant, tx, &EffectsConfig::default())
        .expect("variant overlapping the transcript must annotate")
}

/// `10:14951238 A>TACTT...` (200 bp): REF of 1 base, so `_get_differing_regions`
/// returns `{0, 0}` and Perl tests the single exonic position; the inserted bases
/// never reach intron 1 (10 bp away). Extending the span by the net insertion
/// length would add intron_variant plus a splice site.
#[test]
fn concordance_one_base_ref_long_alt_tests_one_position() {
    let tx = strip_coding_model(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_000_290,
        25_000_290,
        b"C".to_vec(),
        b"A".repeat(41),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(&tc, &["non_coding_transcript_exon_variant"]);
}

/// `10:8106103 T>GCTTACTTCCC` on ENST00000346208: the single position is
/// intron_start+1, inside the donor dinucleotide and nothing else; a span
/// extended by the 10 inserted bases would reach the 5th base and the intron core.
#[test]
fn concordance_one_base_ref_long_alt_at_donor_plus_one_is_donor_only() {
    let tx = strip_coding_model(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_000_301,
        25_000_301,
        b"T".to_vec(),
        b"GCTTACTTCCC".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["splice_donor_variant", "non_coding_transcript_variant"],
    );
}

/// `11:34988343 TGG>CGGACAACCCAATGCAGTGGTAGTGT` on ENST00000227868: both alleles
/// exceed 1 base, so the XOR region runs 23 bases past the reference span and
/// covers the next intron's donor and 5th base. Perl still emits nothing intronic:
/// `_overlapped_introns` is memoised on the raw span (34988343-34988345), 17 bp
/// short of the intron tree window `[intron_start-3, intron_end+3]`, so no intron
/// is fetched for the region loop to test.
#[test]
fn concordance_delins_far_from_intron_fetches_no_intron() {
    let tx = strip_coding_model(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_000_280,
        25_000_282,
        b"AAA".to_vec(),
        b"C".repeat(26),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(&tc, &["non_coding_transcript_exon_variant"]);
}

/// `11:45975175 G>AGA` on ENST00000257821: the single position is the second
/// exonic base before the intron, inside `_intron_overlap`'s
/// `[intron_start-3, intron_start-1]` window and outside the donor dinucleotide.
/// Letting the two inserted bases reach the donor would emit splice_donor and
/// drop splice_region in its favour.
#[test]
fn concordance_one_base_ref_two_before_donor_is_splice_region_not_donor() {
    let tx = strip_coding_model(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_000_298,
        25_000_298,
        b"A".to_vec(),
        b"GGG".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "non_coding_transcript_exon_variant",
            "splice_region_variant",
        ],
    );
}

/// `17:7808415 C>AA` on ENST00000330494: the single position is intron_end-2,
/// inside the polypyrimidine tract and the `[intron_end-7, intron_end-2]`
/// splice_region window, one base short of the acceptor dinucleotide that a
/// span extended by the inserted base would reach.
#[test]
fn concordance_one_base_ref_at_acceptor_minus_two_is_ppt_and_splice_region() {
    let tx = strip_coding_model(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_001_997,
        25_001_997,
        b"C".to_vec(),
        b"AA".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "intron_variant",
            "splice_polypyrimidine_tract_variant",
            "splice_region_variant",
            "non_coding_transcript_variant",
        ],
    );
}

/// `11:57379268-57379275` (8 bp REF, 40 bp ALT) on ENST00000340687: the raw span
/// lies deep in the intron, so only the intron list holds the intron and the
/// boundary list is empty; the XOR region crosses the acceptor into the next exon
/// but the boundary loop never runs. Perl keeps intron_variant and the
/// polypyrimidine tract (the region does reach `[intron_end-16, intron_end-2]`)
/// and never emits the acceptor the extended region touches.
#[test]
fn concordance_delins_deep_in_intron_reaching_acceptor_is_ppt_only() {
    let tx = strip_coding_model(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_001_975,
        25_001_977,
        b"AAA".to_vec(),
        b"C".repeat(30),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "intron_variant",
            "splice_polypyrimidine_tract_variant",
            "non_coding_transcript_variant",
        ],
    );
}

/// `12:32994929-33003861` on ENST00000070846: a deletion spanning two introns
/// where only the second reaches the 5th base. `splice_donor_region_variant`
/// (`VariationEffect`) returns 0 whenever the transcript-wide 5th-base flag
/// is set, so the donor-region window the first intron alone satisfies emits
/// nothing; suppressing per intron would emit both.
#[test]
fn concordance_fifth_base_in_one_intron_suppresses_donor_region_in_another() {
    let tx = strip_coding_model(make_test_transcript());
    // Intron 1 positions +5..; exon 2; intron 2 positions ..+4 (the 5th base).
    let variant = InputVariant::new(
        "21".into(),
        25_000_305,
        25_002_304,
        b"A".repeat(2000),
        b"-".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "intron_variant",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "splice_donor_5th_base_variant",
            "non_coding_transcript_exon_variant",
        ],
    );
}

/// `19:11213339-11217362` on ENST00000252444: a deletion from deep in intron 1
/// (past its 5th base) through exon 2 that stops 3 bp before intron 2, inside
/// its exonic `[intron_start-3, intron_start-1]` window. `splice_region`
/// (`VariationEffect`) returns 0 once any intron set a donor or
/// acceptor flag, so the window intron 2 satisfies emits nothing.
#[test]
fn concordance_acceptor_in_one_intron_suppresses_splice_region_from_another() {
    let tx = strip_coding_model(make_test_transcript());
    let variant = InputVariant::new(
        "21".into(),
        25_001_990,
        25_002_297,
        b"A".repeat(308),
        b"-".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "intron_variant",
            "splice_acceptor_variant",
            "non_coding_transcript_exon_variant",
        ],
    );
}

/// `11:108173722-108173736` (15 bp REF, 40 bp ALT) upstream of ENST00000529588:
/// the XOR region runs 25 bp into the transcript, but `within_feature`
/// (`BaseVariationFeatureOverlapAllele::_bvfo_preds`) tests the raw span, and
/// without it no intron is fetched and no splice predicate runs, so the intron
/// the inserted bases reach is not annotated beside upstream_gene_variant.
#[test]
fn concordance_long_alt_upstream_of_transcript_stays_upstream() {
    let tx = make_test_transcript();
    let mut alt = b"A".repeat(15);
    alt.extend(b"C".repeat(385));
    let variant = InputVariant::new("21".into(), 24_999_985, 24_999_999, b"A".repeat(15), alt);
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(&tc, &["upstream_gene_variant"]);
}

/// `15:72105929 C>-` on ENST00000326995 and `11:118898437 C>A` on
/// ENST00000357590: an exonic base adjacent to a 1 bp intron. The intron is a
/// frameshift intron, but `_intron_effects` skips it (`next`) only for a region
/// that overlaps it; this region does not, so the acceptor window
/// `[intron_end-1, intron_end]`, which for a 1 bp intron includes the preceding
/// exonic base, is tested and holds. Skipping every frameshift intron would
/// report the exonic splice_region window instead.
#[test]
fn concordance_snv_before_one_bp_intron_is_acceptor() {
    let tx = make_transcript_with_short_first_intron(1);
    // cDNA 300 = CDS 250, first base of a GCT codon: G>A gives ACT (Thr).
    let variant = InputVariant::new(
        "21".into(),
        25_000_299,
        25_000_299,
        b"G".to_vec(),
        b"A".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["missense_variant", "splice_acceptor_variant"],
    );
}

/// The counterpart of the test above on the other side of the 1 bp intron
/// (`17:48227385 C>G` on ENST00000316878, reverse strand, is the acceptor form):
/// the donor window `[intron_start, intron_start+1]` includes the first base of
/// the following exon.
#[test]
fn concordance_snv_after_one_bp_intron_is_donor() {
    let tx = make_transcript_with_short_first_intron(1);
    // cDNA 301 = CDS 251, second base of a GCT codon: C>A gives GAT (Asp).
    let variant = InputVariant::new(
        "21".into(),
        25_000_301,
        25_000_301,
        b"C".to_vec(),
        b"A".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["missense_variant", "splice_donor_variant"],
    );
}

/// `11:118939941 C>CT` on ENST00000300793: an insertion between the last exonic
/// base and a 2 bp intron. The inverted region `(intron_start, intron_start-1)`
/// does not overlap the frameshift intron, so it is tested: the `$insertion`
/// clause `r_end == intron_end-2` of `_intron_effects`
/// holds because `intron_end-2` is the exonic anchor, and `_intron_overlap`'s
/// `vf_start == intron_start` holds. Skipping the intron, or guarding the core
/// window with `core_start <= core_end`, would drop both terms.
#[test]
fn concordance_insertion_before_two_bp_intron_is_intronic_and_splice_region() {
    let tx = strip_coding_model(make_transcript_with_short_first_intron(2));
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_299,
        b"-".to_vec(),
        b"T".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "intron_variant",
            "splice_region_variant",
            "non_coding_transcript_variant",
        ],
    );
}

/// `20:35807790 G>GCTTACAG...` on ENST00000441008: an insertion 6 bp before a
/// 1 bp intron. Both anchors sit inside `_intron_overlap`'s `[intron_end-7,
/// intron_end-2]` window, which for a 1 bp intron spills into the exon, and the
/// non-normalising `overlap` holds for the inverted pair. The intron is in the
/// boundary list (`[intron_end-7, intron_end+3]`) and not in the intron list
/// (`[intron_start-3, intron_end+3]`); the boundary loop runs regardless.
#[test]
fn concordance_insertion_six_bases_before_one_bp_intron_gets_splice_region() {
    let tx = strip_coding_model(make_transcript_with_short_first_intron(1));
    let variant = InputVariant::new(
        "21".into(),
        25_000_295,
        25_000_294,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "non_coding_transcript_exon_variant",
            "splice_region_variant",
        ],
    );
}

/// `22:36587845 A>ACTCT` on ENST00000332987 (reverse strand) in forward form: an
/// insertion between a 1 bp intron and the following exon. The inverted pair
/// `(intron_start+1, intron_start)` does not overlap the intron (`r_start <=
/// intron_end` fails), and both anchors sit inside the donor window
/// `[intron_start, intron_start+1]`.
#[test]
fn concordance_insertion_after_one_bp_intron_is_donor() {
    let tx = strip_coding_model(make_transcript_with_short_first_intron(1));
    let variant = InputVariant::new(
        "21".into(),
        25_000_301,
        25_000_300,
        b"-".to_vec(),
        b"CTCT".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["splice_donor_variant", "non_coding_transcript_variant"],
    );
}

/// `15:23890698 G>GA` on ENST00000314233 (reverse strand) in forward form: an
/// insertion two exonic bases past a 1 bp intron. Both anchors sit inside the
/// donor-region window `[intron_start+2, intron_start+5]`, the 5th base
/// `intron_start+4` is outside the pair, and donor_region suppresses the
/// splice_region window `[intron_start+2, intron_start+7]` the pair also fills.
#[test]
fn concordance_insertion_two_bases_after_one_bp_intron_is_donor_region() {
    let tx = strip_coding_model(make_transcript_with_short_first_intron(1));
    let variant = InputVariant::new(
        "21".into(),
        25_000_303,
        25_000_302,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &[
            "non_coding_transcript_exon_variant",
            "splice_donor_region_variant",
        ],
    );
}

/// `15:72105930 C>T` on ENST00000326995 (processed_transcript): an SNV that is
/// the 1 bp intron. The region overlaps the frameshift intron, so
/// `_intron_effects` skips it in both loops and no intronic flag is set;
/// `VariationEffect::non_coding_exon_variant` needs the raw span to
/// overlap an exon and fails, so `within_non_coding_gene` yields
/// non_coding_transcript_variant, not the exon term.
#[test]
fn concordance_snv_inside_one_bp_intron_of_non_coding_transcript_is_transcript_variant() {
    let tx = strip_coding_model(make_transcript_with_short_first_intron(1));
    let variant = InputVariant::new(
        "21".into(),
        25_000_300,
        25_000_300,
        b"G".to_vec(),
        b"A".to_vec(),
    );
    let tc = consequence_of(&variant, &tx);
    crate::test_helpers::assert_consequence_set_eq(&tc, &["non_coding_transcript_variant"]);
}

// CDS-boundary spans and the UTR terms.
//
// Perl's `VariationEffect::within_5_prime_utr` / `within_3_prime_utr`
// test `overlap(bvf_start, bvf_end, region_start, region_end)` against the
// genomic UTR region with no check that the region is non-empty. On a
// transcript whose CDS reaches its boundary the region is inverted, and the
// test holds for a span that crosses the boundary point and fails for one that
// only touches it; the `cds_*_NF` flags never enter.

/// The test transcript with its CDS extended to the transcript's 3' end, so the
/// 3' UTR region `cre + 1 .. transcript.end` is empty.
fn make_transcript_without_three_prime_utr() -> Transcript {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        let mut cds = vefc.translateable_seq.take().unwrap_or_default();
        while cds.len() < 2551 {
            cds.push_str("GCT");
        }
        cds.truncate(2551);
        vefc.translateable_seq = Some(cds);
        if let Some(ref mut mapper) = vefc.mapper {
            mapper.cdna_coding_end = 2601;
        }
    }
    tx.cdna_coding_end = Some(2601);
    tx.coding_region_end = Some(25_006_000);
    tx.translation_end = Some(25_006_000);
    tx
}

/// The test transcript with its CDS extended to the transcript's 5' end, so the
/// 5' UTR region `transcript.start .. crs - 1` is empty.
fn make_transcript_without_five_prime_utr() -> Transcript {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        let tail = vefc.translateable_seq.take().unwrap_or_default();
        let mut cds = String::with_capacity(900);
        while cds.len() < 48 {
            cds.push_str("GCT");
        }
        cds.push_str("GC");
        cds.push_str(&tail);
        vefc.translateable_seq = Some(cds);
        if let Some(ref mut mapper) = vefc.mapper {
            mapper.cdna_coding_start = 1;
        }
    }
    tx.cdna_coding_start = Some(1);
    tx.coding_region_start = Some(25_000_000);
    tx.translation_start = Some(25_000_000);
    tx
}

/// A deletion from the last CDS bases across the transcript's 3' end, on a
/// transcript with no 3' UTR, is `3_prime_UTR_variant,coding_sequence_variant`
/// (Perl: `15:25584272-25584286 - ENST00000397954`): the
/// trailing mapper gap leaves `cds_end` and `cdna_end` undefined, so no
/// frameshift or stop_lost, while `_after_coding` holds on the empty region.
#[test]
fn concordance_deletion_crossing_transcript_end_without_utr_gets_three_prime_utr() {
    let tx = make_transcript_without_three_prime_utr();
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_005_990,
        25_006_010,
        b"N".repeat(21),
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["3_prime_UTR_variant", "coding_sequence_variant"],
    );
}

/// A deletion ending exactly on the transcript's last base of a `cds_end_NF`
/// transcript touches the empty 3' UTR region without crossing it, so Perl
/// emits `frameshift_variant` alone (`11:17574982-17574991 - ENST00000428619`);
/// the NF flag adds nothing.
#[test]
fn concordance_deletion_touching_transcript_end_of_cds_end_nf_has_no_three_prime_utr() {
    let mut tx = make_transcript_without_three_prime_utr();
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_005_991,
        25_006_000,
        b"N".repeat(10),
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["frameshift_variant"]);
}

/// A deletion from far upstream (past the 5 kb flank, so the start endpoint maps
/// nowhere) into the first CDS bases of a transcript with no 5' UTR is
/// `5_prime_UTR_variant,coding_sequence_variant` (Perl:
/// `MT:5782-13922 - ENST00000361567`).
#[test]
fn concordance_deletion_crossing_transcript_start_without_utr_gets_five_prime_utr() {
    let tx = make_transcript_without_five_prime_utr();
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        24_994_000,
        25_000_010,
        b"N".repeat(6011),
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["5_prime_UTR_variant", "coding_sequence_variant"],
    );
}

/// A 2 bp deletion of the first two transcript bases of a `cds_start_NF`
/// transcript with no 5' UTR is `frameshift_variant` alone (Perl:
/// `16:27585228-27585229 - ENST00000568258`).
#[test]
fn concordance_deletion_touching_transcript_start_of_cds_start_nf_has_no_five_prime_utr() {
    let mut tx = make_transcript_without_five_prime_utr();
    tx.flags = vec!["cds_start_NF".to_string()].into();
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_000,
        25_000_001,
        b"GC".to_vec(),
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["frameshift_variant"]);
}

// SeqEdits on the reference peptide.
//
// `TranscriptVariationAllele::peptide` applies the translation's SeqEdits to the
// reference peptide, so a selenocysteine `TGA` reads `U`. The JSON cache has no
// `seq_edits`, but its `peptide` is the edited translation, and the edited
// residue is read from it.

/// The test transcript with codon 4 (`AAA`) replaced by `TGA` and the cached
/// translation carrying `residue` there, as Ensembl's `Translation->seq` does for
/// a `_selenocysteine` (`U`) or amino-acid-substitution (`X`) SeqEdit.
fn make_transcript_with_edited_codon_4(residue: u8) -> Transcript {
    let mut tx = make_test_transcript();
    if let Some(ref mut vefc) = tx.vefc {
        let mut cds = vefc.translateable_seq.take().unwrap_or_default();
        cds.replace_range(9..12, "TGA");
        let mut pep: Vec<u8> = cds
            .as_bytes()
            .chunks_exact(3)
            .map(|c| vep_core::codon::translate_codon_with_table(c, 1))
            .collect();
        pep[3] = residue;
        vefc.peptide = Some(String::from_utf8(pep).expect("ascii"));
        vefc.translateable_seq = Some(cds);
    }
    tx
}

/// `TGA` -> `TGG` at a selenocysteine codon is `missense_variant` with amino
/// acids `U/W`, not `stop_lost` (Perl: `19:48283989 G ENST00000593892`;
/// `1:26139282 G ENST00000361547`).
#[test]
fn concordance_snv_at_selenocysteine_codon_is_missense_not_stop_lost() {
    let tx = make_transcript_with_edited_codon_4(b'U');
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_061,
        25_000_061,
        b"A".to_vec(),
        b"G".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["missense_variant"]);
    assert_eq!(tc.amino_acids.as_deref(), Some("U/W"));
}

/// `TGA` -> `TAA` at a selenocysteine codon is `stop_gained`, not
/// `stop_retained_variant` (Perl: `1:26139281 A ENST00000361547`).
#[test]
fn concordance_snv_turning_selenocysteine_codon_into_stop_is_stop_gained() {
    let tx = make_transcript_with_edited_codon_4(b'U');
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_060,
        25_000_060,
        b"G".to_vec(),
        b"A".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["stop_gained"]);
    assert_eq!(tc.amino_acids.as_deref(), Some("U/*"));
}

/// An `X` SeqEdit on the reference residue makes Perl's `coding_unknown` hold,
/// and `missense_variant` (`ref_pep ne alt_pep`) alongside it (Perl:
/// `1:161305878 G ENST00000672602`).
#[test]
fn concordance_snv_at_x_edited_codon_is_coding_sequence_and_missense() {
    let tx = make_transcript_with_edited_codon_4(b'X');
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_061,
        25_000_061,
        b"A".to_vec(),
        b"G".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["coding_sequence_variant", "missense_variant"],
    );
}

/// Deleting a whole selenocysteine codon is `inframe_deletion` alone: the
/// reference peptide is `U`, so Perl's peptide-route `stop_lost` does not hold.
#[test]
fn concordance_deletion_of_selenocysteine_codon_is_inframe_deletion_without_stop_lost() {
    let tx = make_transcript_with_edited_codon_4(b'U');
    let config = EffectsConfig::default();
    let variant = InputVariant::new(
        "21".into(),
        25_000_059,
        25_000_061,
        b"TGA".to_vec(),
        b"-".to_vec(),
    );
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, &["inframe_deletion"]);
    assert_eq!(tc.amino_acids.as_deref(), Some("U/-"));
}

// Structural deletions inside one exon.
//
// The structural-variant arms of Perl's `VariationEffect::frameshift`
// and `inframe_deletion` fire for a deletion lying completely
// within one exon, by its own length mod 3, and either excludes `coding_unknown`;
// the same arm of `stop_lost` is a genomic overlap with the three stop-codon bases.

/// A 68 bp `<DEL>` inside the last CDS exon ending on the first stop-codon base
/// is `feature_truncation,frameshift_variant,stop_lost` (Perl:
/// `21:33488228-33488295 deletion ENST00000314399`).
#[test]
fn concordance_structural_deletion_within_one_exon_is_frameshift_and_stop_lost() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();
    let mut variant = InputVariant::new(
        "21".into(),
        25_004_230,
        25_004_297,
        b"N".to_vec(),
        b"-".to_vec(),
    );
    variant.variant_class = vep_core::variant::VariantClass::StructuralDeletion;
    variant.is_structural = true;
    variant.sv_end = Some(25_004_297);
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["feature_truncation", "frameshift_variant", "stop_lost"],
    );
}

/// A 624 bp `<DEL>` inside one CDS exon is `feature_truncation,inframe_deletion`
/// (Perl: `21:44574108-44574731 deletion ENST00000400374`).
#[test]
fn concordance_structural_deletion_within_one_exon_multiple_of_three_is_inframe_deletion() {
    let tx = make_test_transcript();
    let config = EffectsConfig::default();
    let mut variant = InputVariant::new(
        "21".into(),
        25_002_010,
        25_002_075,
        b"N".to_vec(),
        b"-".to_vec(),
    );
    variant.variant_class = vep_core::variant::VariantClass::StructuralDeletion;
    variant.is_structural = true;
    variant.sv_end = Some(25_002_075);
    let tc = calculate_consequences(&variant, &tx, &config).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(
        &tc,
        &["feature_truncation", "inframe_deletion"],
    );
}

// Start-codon SNVs: Perl's three `start_lost` routes and `start_retained_variant`
// on the SNV path (`coding::perl_coding_terms` for an SNV in the first codon).

/// `make_test_transcript()` with its CDS opening on `first_codons` and, when
/// `utr_len` is 0, the CDS beginning at cDNA 1 (no 5' UTR; the 850 bp CDS then
/// ends at cDNA 850, genomic 25_004_249). `table` is the codon table and
/// `peptide` the cached transcript peptide, or none.
fn make_start_codon_transcript(
    first_codons: &str,
    utr_len: u64,
    table: u8,
    peptide: Option<&str>,
) -> Transcript {
    let mut tx = make_test_transcript();
    let mut cds = String::with_capacity(850);
    cds.push_str(first_codons);
    while cds.len() < 850 {
        cds.push_str("GCT");
    }
    cds.truncate(850);
    let vefc = tx
        .vefc
        .as_mut()
        .expect("test transcript has a coding model");
    vefc.translateable_seq = Some(cds);
    vefc.codon_table = table;
    vefc.peptide = peptide.map(|p| {
        let mut pep = String::from(p);
        while pep.len() < 283 {
            pep.push('A');
        }
        pep
    });
    if utr_len == 0 {
        let mapper = vefc.mapper.as_mut().expect("mapper");
        mapper.cdna_coding_start = 1;
        mapper.cdna_coding_end = 850;
        tx.cdna_coding_start = Some(1);
        tx.cdna_coding_end = Some(850);
        tx.coding_region_start = Some(25_000_000);
        tx.coding_region_end = Some(25_004_249);
        tx.translation_start = Some(25_000_000);
        tx.translation_end = Some(25_004_249);
        if let Some(ref mut translation) = tx.translation {
            translation.start = 1;
            translation.end = 250;
        }
    } else {
        assert_eq!(utr_len, 50, "the fixture's 5' UTR is 50 bp");
    }
    if table == 2 {
        tx.chr = "MT".into();
    }
    tx
}

fn assert_start_codon_snv(
    tx: &Transcript,
    pos: u64,
    ref_base: &str,
    alt_base: &str,
    want: &[&str],
) {
    let variant = InputVariant::new(
        tx.chr.to_string(),
        pos,
        pos,
        ref_base.as_bytes().to_vec(),
        alt_base.as_bytes().to_vec(),
    );
    let tc = calculate_consequences(&variant, tx, &EffectsConfig::default()).expect("annotates");
    crate::test_helpers::assert_consequence_set_eq(&tc, want);
}

/// `_overlaps_start_codon` is 0 on a `cds_start_NF` transcript, so no start
/// predicate fires and a codon-1 `ATG>GTG` is `missense_variant` (Perl:
/// `11:7654098 G ENST00000530081`, CDS at cDNA 1; `1:216595678 C
/// ENST00000307340`, CDS behind a 387 bp UTR); a verdict read from the Met
/// alone would be `start_lost`.
#[test]
fn concordance_start_codon_snv_on_cds_start_nf_transcript_is_missense() {
    let mut no_utr = make_start_codon_transcript("ATGCCTCCA", 0, 1, None);
    no_utr.flags = vec!["cds_start_NF".to_string(), "cds_end_NF".to_string()].into();
    assert_start_codon_snv(&no_utr, 25_000_000, "A", "G", &["missense_variant"]);

    let with_utr = make_test_transcript_with_flags(&["cds_start_NF", "cds_end_NF"]);
    assert_start_codon_snv(&with_utr, 25_000_050, "A", "G", &["missense_variant"]);
}

/// A non-ATG start codon behind a 5' UTR: `_inv_start_altered` edits the UTR +
/// CDS string and reads `start_lost` because the codon at the UTR length is not
/// `ATG`, whatever the residues; `missense_variant` is then suppressed (Perl:
/// `11:22647057 C ENST00000428556`, `GTT>CTT` V/L; `11:32456890 G
/// ENST00000332351`, WT1 `CTG>CCG` L/P); a verdict read from the residues alone
/// would be `missense_variant`.
#[test]
fn concordance_non_atg_start_codon_missense_snv_with_utr_is_start_lost() {
    let tx = make_start_codon_transcript("CTGCAGGAC", 50, 1, None);
    assert_start_codon_snv(&tx, 25_000_051, "T", "C", &["start_lost"]);
}

/// The same route on a synonymous change: `start_lost` from the edited codon
/// not being `ATG`, and `synonymous_variant` alongside it because the residues
/// match and `start_retained_variant` (edited codon eq `ATG`) is false (Perl:
/// `11:32456889 A ENST00000332351`, WT1 `CTG>CTT` L; `20:30640228 T
/// ENST00000375852`, `CTG>TTG` L); a verdict read from the residues alone would
/// be `synonymous_variant`.
#[test]
fn concordance_non_atg_start_codon_synonymous_snv_with_utr_is_start_lost_and_synonymous() {
    let tx = make_start_codon_transcript("CTGCAGGAC", 50, 1, None);
    assert_start_codon_snv(
        &tx,
        25_000_052,
        "G",
        "T",
        &["start_lost", "synonymous_variant"],
    );
}

/// With no 5' UTR `_inv_start_altered` returns 0 and the peptide route decides:
/// `translation_start == 1` and the alt residue neither starts nor ends with
/// the ref residue, so a codon-1 missense on a non-ATG start is `start_lost`
/// (Perl: `16:2017276 T ENST00000321392`, `GGT>GTT` G/V; `1:237433797 A
/// ENST00000542537`, `GAT>AAT` D/N); a verdict read from the residues alone
/// would be `missense_variant`.
#[test]
fn concordance_non_atg_start_codon_missense_snv_without_utr_is_start_lost() {
    let tx = make_start_codon_transcript("GGTGGCGGG", 0, 1, None);
    assert_start_codon_snv(&tx, 25_000_001, "G", "T", &["start_lost"]);
}

/// `_snp_start_altered` compares the edited codon with the literal `ATG`, not
/// with the codon table's start codons: `ATG>ATA` on a mitochondrial transcript
/// keeps Met under table 2 yet is `synonymous_variant`, not
/// `start_retained_variant` (Perl: `MT:8529 A ENST00000361899`); a verdict read
/// from the codon table's start codons would be `start_retained_variant`.
#[test]
fn concordance_mt_start_codon_atg_to_ata_is_synonymous_not_start_retained() {
    let tx = make_start_codon_transcript("ATGAACGAA", 0, 2, Some("MNE"));
    assert_start_codon_snv(&tx, 25_000_002, "G", "A", &["synonymous_variant"]);
}

/// A mitochondrial `ATA` start (Met under table 2) mutated to `AGA` (a table-2
/// stop): `stop_gained`, and `start_lost` from the peptide route because `*`
/// contains no `M` (Perl: `MT:3308 G ENST00000361390`); a verdict read from the
/// residues alone would be `stop_gained` alone.
#[test]
fn concordance_mt_start_codon_to_stop_is_start_lost_and_stop_gained() {
    let tx = make_start_codon_transcript("ATACCCATG", 0, 2, Some("MPM"));
    assert_start_codon_snv(&tx, 25_000_001, "T", "G", &["start_lost", "stop_gained"]);
}

/// An `initial_met` SeqEdit is read back from the cached peptide: the peptide
/// starts with `M` while `GTG` translates to `V` and is not a table-1 start
/// codon, so the reference residue is `M` and `GTG>ATG` is
/// `start_retained_variant` alone, the peptide route seeing `M` in `M` (Perl:
/// `20:58853261 A ENST00000306120`); a verdict read from the codon translation
/// alone would be `missense_variant`.
#[test]
fn concordance_initial_met_edit_gtg_to_atg_is_start_retained() {
    let tx = make_start_codon_transcript("GTGTTATGG", 0, 1, Some("MLW"));
    assert_start_codon_snv(&tx, 25_000_000, "G", "A", &["start_retained_variant"]);
}

/// The same `GTG>ATG` on a transcript whose peptide starts with `V` carries no
/// edit: Perl reads `start_lost` (V is not in M) together with
/// `start_retained_variant` (the new codon is `ATG`), a co-emission on a
/// sequence variant, and vep-rs keeps `start_retained_variant` of that pair:
/// the codon after the edit is `ATG`, so the start is retained and the
/// peptide-route `start_lost` is the erroneous member. A `CTG` start is a
/// table-1 start codon `Transcript::translate` reads as `M` itself, so a
/// peptide starting with `M` there is no evidence of an edit and the ref
/// residue stays `L`.
#[test]
fn concordance_start_codon_snv_to_atg_without_edit_is_start_retained_only() {
    let no_edit = make_start_codon_transcript("GTGTTATGG", 0, 1, Some("VLW"));
    assert_start_codon_snv(&no_edit, 25_000_000, "G", "A", &["start_retained_variant"]);

    let forced_met = make_start_codon_transcript("CTGCAGGAC", 50, 1, Some("MQD"));
    assert_start_codon_snv(
        &forced_met,
        25_000_052,
        "G",
        "T",
        &["start_lost", "synonymous_variant"],
    );
}

/// `perl_coding_terms` takes an SNV only when its position maps into cDNA
/// `cdna_coding_start ..= cdna_coding_start + 2`: the three first-codon bases
/// (cDNA 51-53) return a term set, the last UTR base and the first base of
/// codon 2 return `None` and leave the SNV path's own verdict in place.
#[test]
fn concordance_perl_coding_terms_gates_snvs_on_the_first_codon() {
    let tx = make_test_transcript();
    let snv = |pos: u64, r: &str, a: &str| {
        let variant = InputVariant::new(
            "21".into(),
            pos,
            pos,
            r.as_bytes().to_vec(),
            a.as_bytes().to_vec(),
        );
        crate::coding::perl_coding_terms(&variant, &tx, pos, pos, None)
    };
    assert!(snv(25_000_049, "A", "G").is_none());
    assert_eq!(
        snv(25_000_050, "A", "G"),
        Some(vec![Consequence::StartLost])
    );
    assert_eq!(
        snv(25_000_051, "T", "C"),
        Some(vec![Consequence::StartLost])
    );
    assert_eq!(
        snv(25_000_052, "G", "A"),
        Some(vec![Consequence::StartLost])
    );
    assert!(snv(25_000_053, "G", "T").is_none());
}
