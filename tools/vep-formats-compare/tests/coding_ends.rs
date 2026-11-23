// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The coding ends read of made-up cache records (`coding_ends`): each stop end and each
//! start end, the steps that decide them, and the genomic positions of both ends on each
//! strand, with a last codon an intron splits and a first codon inside a block of bases
//! the genome lacks.
//!
//! A record is one exon of `leader` cDNA bases before its coding sequence and `trailer`
//! after it, whose first cDNA base on the forward strand, or last on the reverse, is the
//! genomic base `FIRST_BASE`, unless a test gives it mapper pairs of its own.

use vep_formats_compare::coding_ends::{
    coding_ends, CodingEnds, EditFlanks, MapperFacts, RecordFacts, StartEnd, StopEnd,
};

const FIRST_BASE: i64 = 1001;

/// A start codon of every code, a codon of alanine and a stop codon of the standard code.
const START: &str = "ATG";
const SENSE: &str = "GCT";
const STOP: &str = "TAA";

/// A mapper pair: a run of cDNA coordinates, the genomic bases it maps to and their
/// orientation.
type Pair = (i64, i64, i64, i64, i64);

#[derive(Clone)]
struct Made {
    codes: Vec<&'static str>,
    consensus_id: Option<&'static str>,
    codon_table: Option<i64>,
    coding: Option<String>,
    peptide: Option<&'static str>,
    strand: i64,
    start_phase: i64,
    leader: i64,
    trailer: i64,
    pairs: Option<Vec<Pair>>,
    mapper: bool,
}

/// A record of the forward strand whose coding sequence is `pieces` joined, with a cached
/// peptide that opens with methionine, ten cDNA bases before the coding sequence and none
/// after it.
fn made(pieces: &[&str]) -> Made {
    Made {
        codes: Vec::new(),
        consensus_id: None,
        codon_table: None,
        coding: Some(pieces.concat()),
        peptide: Some("MA"),
        strand: 1,
        start_phase: 0,
        leader: 10,
        trailer: 0,
        pairs: None,
        mapper: true,
    }
}

fn read(made: &Made) -> CodingEnds {
    let on_the_cdna = made
        .coding
        .as_deref()
        .map_or(0, |coding| coding.len() as i64)
        - made.start_phase.max(0);
    let (coding_start, coding_end) = (made.leader + 1, made.leader + on_the_cdna);
    let last = coding_end + made.trailer;
    let one_exon = [(1, last, FIRST_BASE, FIRST_BASE + last - 1, made.strand)];
    coding_ends(&RecordFacts {
        attribute_codes: &made.codes,
        consensus_id: made.consensus_id,
        strand: Some(made.strand),
        codon_table: made.codon_table,
        coding: made.coding.as_deref(),
        peptide: made.peptide,
        mapper: made.mapper.then_some(MapperFacts {
            start_phase: made.start_phase,
            coding_start,
            coding_end,
            pairs: made.pairs.as_deref().unwrap_or(&one_exon),
        }),
        span: None,
        introns: &[],
    })
}

const NONE: StopEnd = StopEnd::NoneShown { flagged: false };
const NOT_SHOWN: StartEnd = StartEnd::NotShown {
    flagged: false,
    padded: false,
};

/// The stop end of a tail of `bases` bases that adenines complete to a stop codon: begun
/// by the tail where the record shows the model complete, and none where it does not.
fn begun(bases: u8, shown: bool) -> StopEnd {
    match shown {
        true => StopEnd::BegunByTheTail { bases },
        false => NONE,
    }
}

#[test]
fn whole_codons_are_closed_by_a_stop_codon_of_the_code_alone() {
    assert_eq!(read(&made(&[START, SENSE, STOP])).stop, StopEnd::Whole);
    assert_eq!(read(&made(&[START, SENSE, SENSE])).stop, NONE);
    // A codon that closes the standard code and names a residue under the vertebrate
    // mitochondrial code, and one the other way about.
    for (last, standard, mitochondrial) in [
        ("TGA", StopEnd::Whole, NONE),
        ("AGA", NONE, StopEnd::Whole),
        // Every reading of an ambiguous codon is a stop codon, or the codon is none.
        ("TAR", StopEnd::Whole, StopEnd::Whole),
        ("TRA", StopEnd::Whole, NONE),
        ("TNA", NONE, NONE),
    ] {
        let mut record = made(&[START, SENSE, last]);
        assert_eq!(read(&record).stop, standard, "{last}");
        // A record without a `codon_table`, and one of 0, is read under the standard code.
        record.codon_table = Some(0);
        assert_eq!(read(&record).stop, standard, "{last}");
        record.codon_table = Some(1);
        assert_eq!(read(&record).stop, standard, "{last}");
        record.codon_table = Some(2);
        assert_eq!(read(&record).stop, mitochondrial, "{last}");
    }
    // A stop codon inside the coding sequence closes nothing, and the last codon closes a
    // sequence that holds one.
    assert_eq!(read(&made(&[START, STOP, SENSE])).stop, NONE);
    assert_eq!(read(&made(&[START, STOP, STOP])).stop, StopEnd::Whole);
    // A lower-case sequence reads as its upper case.
    assert_eq!(read(&made(&[START, SENSE, "taa"])).stop, StopEnd::Whole);
}

#[test]
fn a_tail_begins_a_stop_codon_where_adenines_complete_it_on_the_last_base() {
    for (tail, bases) in [("T", 1), ("TA", 2), ("TG", 2), ("t", 1)] {
        let mut record = made(&[START, SENSE, tail]);
        // A record that does not show the model complete shows no stop codon.
        assert_eq!(read(&record).stop, NONE, "{tail}");
        record.consensus_id = Some("CCDS1.1");
        assert_eq!(read(&record).stop, begun(bases, true), "{tail}");
        assert_eq!(read(&record).tail_length, bases, "{tail}");
        // A 3' UTR follows the tail: the transcript does not end inside the codon.
        record.trailer = 1;
        assert_eq!(read(&record).stop, NONE, "{tail} before a 3' UTR");
        record.trailer = 0;
        // No mapper, a coding end the mapper does not state, and a cDNA that ends before
        // the coding end.
        record.mapper = false;
        assert_eq!(read(&record).stop, NONE, "{tail} without a mapper");
        record.mapper = true;
        record.pairs = Some(vec![(1, 5, FIRST_BASE, FIRST_BASE + 4, 1)]);
        assert_eq!(read(&record).stop, NONE, "{tail} on a shorter cDNA");
    }
    // Adenines complete these tails to a codon that names a residue.
    for tail in ["G", "GC", "AG", "N", "TN"] {
        let mut record = made(&[START, SENSE, tail]);
        record.consensus_id = Some("CCDS1.1");
        assert_eq!(read(&record).stop, NONE, "{tail}");
    }
    // A coding end of 0 is one the mapper does not state.
    let mut unstated = made(&["T"]);
    unstated.consensus_id = Some("CCDS1.1");
    unstated.leader = -1;
    unstated.pairs = Some(vec![(0, 0, FIRST_BASE, FIRST_BASE, 1)]);
    assert_eq!(read(&unstated).stop, NONE);
}

#[test]
fn a_begun_stop_codon_is_shown_by_a_mitochondrial_code_or_by_a_consensus_id() {
    for (codon_table, shown) in [
        (None, false),
        (Some(1), false),
        (Some(2), true),
        (Some(4), false),
        (Some(5), true),
        (Some(9), true),
        (Some(11), false),
        (Some(13), true),
        (Some(26), false),
    ] {
        for (tail, bases) in [("T", 1), ("TA", 2)] {
            let mut record = made(&[START, SENSE, tail]);
            record.codon_table = codon_table;
            assert_eq!(
                read(&record).stop,
                begun(bases, shown),
                "{codon_table:?} {tail}"
            );
            // A consensus id shows the model complete under any code; an empty one and a
            // dash are none.
            record.consensus_id = Some("CCDS1.1");
            assert_eq!(
                read(&record).stop,
                begun(bases, true),
                "{codon_table:?} {tail}"
            );
            for absent in ["", "-"] {
                record.consensus_id = Some(absent);
                assert_eq!(
                    read(&record).stop,
                    begun(bases, shown),
                    "{codon_table:?} {tail}"
                );
            }
        }
    }
    // A tail that begins a stop codon of the vertebrate mitochondrial code alone, and one
    // of the standard code alone.
    for (tail, standard, mitochondrial) in
        [("AG", NONE, begun(2, true)), ("TG", begun(2, false), NONE)]
    {
        let mut record = made(&[START, SENSE, tail]);
        assert_eq!(read(&record).stop, standard, "{tail}");
        record.codon_table = Some(2);
        assert_eq!(read(&record).stop, mitochondrial, "{tail}");
    }
    // Under a code in which the completed codon names a residue no stop codon is begun,
    // whatever the record carries.
    let mut ciliate = made(&[START, SENSE, "T"]);
    ciliate.codon_table = Some(6);
    ciliate.consensus_id = Some("CCDS1.1");
    assert_eq!(read(&ciliate).stop, NONE);
}

#[test]
fn the_end_flag_shows_no_stop_codon_whatever_the_sequence_holds() {
    for pieces in [
        [START, SENSE, STOP],
        [START, SENSE, SENSE],
        [START, SENSE, "T"],
    ] {
        let mut record = made(&pieces);
        record.codes = vec!["TSL", "cds_end_NF"];
        record.consensus_id = Some("CCDS1.1");
        record.codon_table = Some(2);
        assert_eq!(read(&record).stop, StopEnd::NoneShown { flagged: true });
        // The flag of the other end and a code that only begins like the flag are none.
        record.codes = vec!["cds_start_NF", "cds_end_NF_like", "cds_end"];
        assert_ne!(read(&record).stop, StopEnd::NoneShown { flagged: true });
    }
}

#[test]
fn a_start_codon_is_shown_by_the_peptide_by_the_code_or_by_a_consensus_id() {
    // A cached peptide that opens with methionine, on any first codon.
    assert_eq!(read(&made(&[SENSE, SENSE, STOP])).start, StartEnd::Shown);
    // A cached peptide that opens with another residue.
    let mut other = made(&[SENSE, SENSE, STOP]);
    other.peptide = Some("AA");
    assert_eq!(read(&other).start, NOT_SHOWN);
    other.consensus_id = Some("CCDS1.1");
    assert_eq!(read(&other).start, StartEnd::Shown);
    for absent in ["", "-"] {
        other.consensus_id = Some(absent);
        assert_eq!(read(&other).start, NOT_SHOWN, "{absent:?}");
    }
    // A first codon that is a start codon of the code gives the peptide its methionine:
    // the standard code holds three, its strict form one, the plastid code every codon of
    // their shape.
    for (first, codon_table, shown) in [
        ("CTG", 1, true),
        ("TTG", 1, true),
        ("GTG", 1, false),
        ("CTG", 24, false),
        ("GTG", 11, true),
        ("NTG", 11, true),
        ("NTG", 1, false),
        ("ATA", 2, true),
        ("ATA", 1, false),
    ] {
        let mut record = made(&[first, SENSE, STOP]);
        record.peptide = Some("LA");
        record.codon_table = Some(codon_table);
        let expected = if shown { StartEnd::Shown } else { NOT_SHOWN };
        assert_eq!(read(&record).start, expected, "{first} under {codon_table}");
    }
    // An empty cached peptide opens with nothing, whatever the first codon.
    let mut empty = made(&[START, SENSE, STOP]);
    empty.peptide = Some("");
    assert_eq!(read(&empty).start, NOT_SHOWN);
}

#[test]
fn a_record_without_a_peptide_is_read_by_the_translation_of_its_whole_codons() {
    for (pieces, codon_table, shown) in [
        (vec![START, SENSE, STOP], 1, true),
        (vec![START, SENSE, SENSE], 1, true),
        (vec![START, SENSE, "T"], 1, true),
        (vec![SENSE, SENSE, STOP], 1, false),
        (vec!["CTG", SENSE, STOP], 1, true),
        (vec!["CTG", SENSE, STOP], 24, false),
        // The one residue of a start codon before a stop codon; a stop codon alone
        // translates to nothing, and so do a tail alone and a first codon cut short.
        (vec![START, STOP], 1, true),
        (vec![STOP], 1, false),
        (vec!["AT"], 1, false),
        // An ambiguous first codon whose readings name two residues.
        (vec!["ATN", SENSE, STOP], 1, false),
    ] {
        let mut record = made(&pieces);
        record.peptide = None;
        record.codon_table = Some(codon_table);
        let expected = if shown { StartEnd::Shown } else { NOT_SHOWN };
        assert_eq!(
            read(&record).start,
            expected,
            "{pieces:?} under {codon_table}"
        );
    }
}

#[test]
fn the_start_flag_and_a_padded_first_codon_show_no_start_codon() {
    let not_shown = |flagged, padded| StartEnd::NotShown { flagged, padded };
    let mut record = made(&[START, SENSE, STOP]);
    record.consensus_id = Some("CCDS1.1");
    assert_eq!(read(&record).start, StartEnd::Shown);
    record.codes = vec!["cds_start_NF"];
    assert_eq!(read(&record).start, not_shown(true, false));
    // A start phase above 0 pads the first codon; one of 0 or below pads nothing.
    for (phase, padded) in [(1, true), (2, true), (0, false), (-1, false)] {
        let unknown = "N".repeat(phase.max(0) as usize);
        let mut record = made(&[unknown.as_str(), START, SENSE, STOP]);
        record.start_phase = phase;
        record.consensus_id = Some("CCDS1.1");
        let expected = if padded {
            not_shown(false, true)
        } else {
            StartEnd::Shown
        };
        assert_eq!(read(&record).start, expected, "phase {phase}");
        record.codes = vec!["cds_start_NF"];
        assert_eq!(
            read(&record).start,
            not_shown(true, padded),
            "phase {phase}"
        );
    }
    // The flag of the other end shows a start codon still.
    let mut other = made(&[START, SENSE, STOP]);
    other.codes = vec!["cds_end_NF"];
    assert_eq!(read(&other).start, StartEnd::Shown);
}

#[test]
fn the_code_the_module_lacks_is_read_as_its_own_code() {
    // Its first start codon names leucine under the table a dumper writes the peptide
    // with, and is a start codon of the code itself; the standard code's other start
    // codon of leucine is none of its start codons.
    for (first, cached, shown) in [
        ("CTG", "LA", true),
        ("CTG", "AA", true),
        ("TTG", "LA", false),
        (START, "MA", true),
        (SENSE, "AA", false),
    ] {
        let mut record = made(&[first, SENSE, STOP]);
        record.codon_table = Some(26);
        record.peptide = Some(cached);
        let expected = if shown { StartEnd::Shown } else { NOT_SHOWN };
        assert_eq!(read(&record).start, expected, "{first} cached as {cached}");
        assert_eq!(read(&record).stop, StopEnd::Whole);
    }
}

#[test]
fn a_record_the_engine_reads_otherwise_is_unread() {
    let unread = |record: &Made| {
        let ends = read(record);
        (ends.stop, ends.start, ends.tail_length)
    };
    let both = |tail| (StopEnd::Unread, StartEnd::Unread, tail);
    let mut record = made(&[START, SENSE, "TA"]);
    // A genetic code no table is held for keeps the tail's length and reads nothing more.
    for codon_table in [7, 8, 17, 20, 25, 27, -1] {
        record.codon_table = Some(codon_table);
        assert_eq!(unread(&record), both(2), "{codon_table}");
        assert_eq!(read(&record), CodingEnds::unread(2), "{codon_table}");
    }
    record.codon_table = None;
    record.coding = None;
    assert_eq!(unread(&record), both(0));
    record.coding = Some(String::new());
    assert_eq!(unread(&record), both(0));
    record.codes = vec!["cds_end_NF", "cds_start_NF"];
    assert_eq!(unread(&record), both(0));
    record.codes = Vec::new();
    record.coding = Some([START, "\u{e9}"].concat());
    assert_eq!(unread(&record), both(1));
    // A cached peptide that is not ASCII leaves the start end unread and the stop end read.
    let mut record = made(&[START, SENSE, STOP]);
    record.peptide = Some("\u{e9}");
    assert_eq!(unread(&record), (StopEnd::Whole, StartEnd::Unread, 0));
}

#[test]
fn the_genomic_positions_of_both_ends_follow_the_strand() {
    let mut record = made(&[START, SENSE, SENSE, STOP]);
    record.trailer = 5;
    let forward = read(&record);
    assert_eq!(forward.strand, Some(1));
    assert_eq!(forward.first_codon, [Some(1011), Some(1012), Some(1013)]);
    assert_eq!(forward.last_codon, [Some(1020), Some(1021), Some(1022)]);
    assert_eq!(forward.coding_region_end, Some((1020, 1022)));
    assert_eq!(forward.first_codon_edit, None);
    assert_eq!((forward.tail, forward.tail_bases), ([None; 2], [None; 2]));
    record.strand = -1;
    let reverse = read(&record);
    assert_eq!(reverse.strand, Some(-1));
    assert_eq!(reverse.first_codon, [Some(1017), Some(1016), Some(1015)]);
    assert_eq!(reverse.last_codon, [Some(1008), Some(1007), Some(1006)]);
    assert_eq!(reverse.coding_region_end, Some((1006, 1008)));
    assert_eq!(reverse.first_codon_edit, None);
}

#[test]
fn the_tail_is_placed_with_its_bases_on_each_strand() {
    for (tail, bases) in [("ta", [Some(b'T'), Some(b'A')]), ("g", [Some(b'G'), None])] {
        let mut record = made(&[START, SENSE, tail]);
        let last = FIRST_BASE + 9 + 6 + tail.len() as i64;
        let forward = read(&record);
        assert_eq!(forward.tail_length as usize, tail.len());
        assert_eq!(forward.tail_bases, bases);
        assert_eq!(
            forward.last_codon,
            [Some(last - 2), Some(last - 1), Some(last)]
        );
        record.strand = -1;
        let reverse = read(&record);
        assert_eq!(reverse.tail_bases, bases);
        assert_eq!(
            reverse.last_codon,
            [Some(FIRST_BASE + 2), Some(FIRST_BASE + 1), Some(FIRST_BASE)]
        );
        match tail.len() {
            2 => {
                assert_eq!(forward.tail, [Some(last - 1), Some(last)]);
                assert_eq!(reverse.tail, [Some(FIRST_BASE + 1), Some(FIRST_BASE)]);
            }
            _ => {
                assert_eq!(forward.tail, [Some(last), None]);
                assert_eq!(reverse.tail, [Some(FIRST_BASE), None]);
            }
        }
        // The last three genomic positions of the coding region are those of the last
        // three coding positions where no intron splits them.
        assert_eq!(forward.coding_region_end, Some((last - 2, last)));
        assert_eq!(
            reverse.coding_region_end,
            Some((FIRST_BASE, FIRST_BASE + 2))
        );
    }
}

#[test]
fn a_padded_first_codon_has_genome_bases_for_its_cached_bases_alone() {
    let mut record = made(&["N", "TG", SENSE, STOP]);
    record.start_phase = 1;
    let forward = read(&record);
    // The first coding position is the padding; the cDNA's first coding base is the second.
    assert_eq!(forward.first_codon, [None, Some(1011), Some(1012)]);
    assert_eq!(forward.last_codon, [Some(1016), Some(1017), Some(1018)]);
    assert_eq!(forward.first_codon_edit, None);
    record.strand = -1;
    assert_eq!(read(&record).first_codon, [None, Some(1008), Some(1007)]);
}

#[test]
fn a_last_codon_an_intron_splits_is_apart_from_the_end_of_the_coding_region() {
    let mut record = made(&[START, SENSE, SENSE, STOP]);
    record.trailer = 5;
    // The first exon ends on the first base of the last codon.
    record.pairs = Some(vec![(1, 20, 1001, 1020, 1), (21, 27, 1101, 1107, 1)]);
    let forward = read(&record);
    assert_eq!(forward.stop, StopEnd::Whole);
    assert_eq!(forward.last_codon, [Some(1020), Some(1101), Some(1102)]);
    assert_eq!(forward.coding_region_end, Some((1100, 1102)));
    record.strand = -1;
    record.pairs = Some(vec![(1, 20, 1101, 1120, -1), (21, 27, 1001, 1007, -1)]);
    let reverse = read(&record);
    assert_eq!(reverse.first_codon, [Some(1110), Some(1109), Some(1108)]);
    assert_eq!(reverse.last_codon, [Some(1101), Some(1007), Some(1006)]);
    assert_eq!(reverse.coding_region_end, Some((1006, 1008)));
}

#[test]
fn a_first_codon_inside_an_edit_names_the_genome_bases_beside_the_edit() {
    let mut record = made(&[START, SENSE, SENSE, STOP]);
    record.trailer = 38;
    // Nine bases the genome lacks open the coding sequence.
    record.pairs = Some(vec![(1, 10, 1001, 1010, 1), (20, 60, 1011, 1051, 1)]);
    let forward = read(&record);
    assert_eq!(forward.first_codon, [None; 3]);
    assert_eq!(
        forward.first_codon_edit,
        Some(EditFlanks {
            before: Some(1010),
            after: Some(1011)
        })
    );
    assert_eq!(forward.last_codon, [Some(1011), Some(1012), Some(1013)]);
    record.strand = -1;
    record.pairs = Some(vec![(1, 10, 1051, 1060, -1), (20, 60, 1010, 1050, -1)]);
    assert_eq!(
        read(&record).first_codon_edit,
        Some(EditFlanks {
            before: Some(1051),
            after: Some(1050)
        })
    );
    // The model opens with the edit: no genome base stands before it.
    record.strand = 1;
    record.leader = 0;
    record.trailer = 48;
    record.pairs = Some(vec![(10, 60, 1001, 1051, 1)]);
    assert_eq!(
        read(&record).first_codon_edit,
        Some(EditFlanks {
            before: None,
            after: Some(1001)
        })
    );
    // An edit that holds a part of the first codon is no span that holds it.
    record.leader = 10;
    record.trailer = 38;
    record.pairs = Some(vec![(1, 11, 1001, 1011, 1), (14, 60, 1012, 1058, 1)]);
    let partly = read(&record);
    assert_eq!(partly.first_codon, [Some(1011), None, None]);
    assert_eq!(partly.first_codon_edit, None);
}

#[test]
fn each_end_is_written_in_the_words_of_the_committed_table() {
    for (stop, words) in [
        (StopEnd::Whole, "whole"),
        (begun(1, true), "begun by the tail, 1 base, shown"),
        (begun(2, true), "begun by the tail, 2 bases, shown"),
        (NONE, "none shown"),
        (StopEnd::NoneShown { flagged: true }, "none shown, flagged"),
        (StopEnd::Unread, "unread"),
    ] {
        assert_eq!(stop.to_string(), words);
    }
    let not_shown = |flagged, padded| StartEnd::NotShown { flagged, padded };
    for (start, words) in [
        (StartEnd::Shown, "shown"),
        (NOT_SHOWN, "not shown"),
        (not_shown(true, false), "not shown, flagged"),
        (not_shown(false, true), "not shown, padded"),
        (not_shown(true, true), "not shown, flagged, padded"),
        (StartEnd::Unread, "unread"),
    ] {
        assert_eq!(start.to_string(), words);
    }
}
