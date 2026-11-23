// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvs_shift_sequence_edge_window` on the HGVSp of a frameshift's row: a walk over the
//! clipped slice's cut that runs one base past the shift the sequence admits moves the
//! insertion or deletion Ensembl VEP translates, so its protein notation differs from the
//! candidate's with its HGVSc and its HGVS_OFFSET. The class takes the three fields where
//! each side's protein notation is the one its walk gives on the transcript's coding
//! sequence, and leaves every other protein notation a difference.

mod compare_common;
mod support;

use compare_common::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use vep_formats_compare::compare::{compare_entries, CompareOptions};
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::rules::sequence_edges::SequenceEdges;
use vep_formats_compare::rules::Side;

const EDGE_CLASS: &str = "hgvs_shift_sequence_edge_window";

/// The twelve codons of the made-up coding sequence, the last a stop.
const CODONS: [&str; 12] = [
    "ATG", "GCC", "GCC", "AAA", "CCC", "GCC", "CTA", "AGC", "GCC", "GCC", "GCC", "TAA",
];

/// The eleven residues the coding sequence translates to under the standard code.
const PEPTIDE: [&str; 11] = ["M", "A", "A", "K", "P", "A", "L", "S", "A", "A", "A"];

/// Fifty bases: three before the coding sequence, the first of them the first base of the
/// allele the tests insert, the coding sequence at bases 4 to 39, and eleven of one base.
fn past_bases() -> String {
    "GTT".to_string() + &CODONS.concat() + &"C".repeat(11)
}

/// The other strand of `bases`, read from its last base.
fn other_strand(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// A reference of three regions: `past`, the fifty bases alone; `inner`, the same bases with
/// 1,000 more after them; `back`, the other strand of the fifty bases.
fn write_reference(directory: &Path) -> PathBuf {
    let regions = [
        ("past", past_bases()),
        ("inner", past_bases() + &"C".repeat(1000)),
        ("back", other_strand(&past_bases())),
    ];
    let (mut text, mut index) = (String::new(), String::new());
    for (name, bases) in regions {
        text += &format!(">{name}\n");
        let length = bases.len();
        index += &format!(
            "{name}\t{length}\t{}\t{length}\t{}\n",
            text.len(),
            length + 1
        );
        text += &format!("{bases}\n");
    }
    planted(directory, "reference.fa.fai", &index);
    planted(directory, "reference.fa", &text)
}

/// A candidate cache of coding transcripts of one exon, each coding from its first base to
/// its last: on `past`, one with its coding sequence and peptide under the standard code,
/// one whose cache states neither and one under the vertebrate mitochondrial code; the
/// first again on `inner`; and on `back` the reverse-strand transcript of the same coding
/// sequence.
fn write_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let read = |table: i64| {
        format!(
            r#", "variation_effect_feature_cache": {{"translateable_seq": "{}", "peptide": "{}", "codon_table": {table}}}"#,
            CODONS.concat(),
            PEPTIDE.concat()
        )
    };
    let record = |name: &str, start: i64, strand: i64, read: &str| {
        let end = start + 35;
        format!(
            r#"{{"stable_id": "{name}", "version": 1, "start": {start}, "end": {end}, "strand": {strand}, "biotype": "protein_coding", "exons": [{{"start": {start}, "end": {end}}}], "cdna_coding_start": 1, "cdna_coding_end": 36{read}}}"#
        )
    };
    let shards = [
        (
            "past",
            vec![
                record("ENST00000000021", 4, 1, &read(1)),
                record("ENST00000000023", 4, 1, ""),
                record("ENST00000000024", 4, 1, &read(2)),
            ],
        ),
        ("inner", vec![record("ENST00000000022", 4, 1, &read(1))]),
        ("back", vec![record("ENST00000000025", 12, -1, &read(1))]),
    ];
    for (region, records) in shards {
        let shard = cache.join("transcripts").join(region);
        std::fs::create_dir_all(&shard).expect("the region's directory is made");
        planted(
            &shard,
            "1-1000000.json",
            &format!("[{}]", records.join(", ")),
        );
    }
    cache
}

/// One side of a frameshift's row: the HGVSc body, the protein the HGVSp names with its
/// description, and the offset, an empty text for an offset the side does not write.
#[derive(Clone, Copy)]
struct Written<'a> {
    notation: &'a str,
    protein: &'a str,
    description: &'a str,
    offset: &'a str,
}

/// What Ensembl VEP writes for the two-base insertion after the ninth coding base: its
/// walk moves one base, so the bases stand turned after the tenth, the fourth codon reads
/// lysine as before and the fifth is the first to change.
const REFERENCE: Written = Written {
    notation: "c.10_11insAG",
    protein: "ENSP00000000021.1",
    description: "p.Pro5AsnfsTer4",
    offset: "1",
};

/// What the candidate writes for it: the sequence admits no shift, so the bases stand
/// after the ninth coding base and the fourth codon is the first to change.
const CANDIDATE: Written = Written {
    notation: "c.9_10insGA",
    protein: "ENSP00000000021.1",
    description: "p.Lys4GlufsTer5",
    offset: "",
};

/// Where a row stands: its Location, its allele, its strand and its transcript.
#[derive(Clone, Copy)]
struct Where<'a> {
    location: &'a str,
    allele: &'a str,
    strand: &'a str,
    feature: &'a str,
}

/// The insertion of two bases between bases 12 and 13 of `past`, the ninth and tenth
/// coding bases of the forward-strand transcript.
const PAST: Where = Where {
    location: "past:12-13",
    allele: "GA",
    strand: "1",
    feature: "ENST00000000021",
};

/// One default-format Transcript row of a frameshift with the three fields a side writes.
fn row(written: Written, at: Where) -> String {
    let mut extra = format!(
        "IMPACT=HIGH;STRAND={};HGVSc={}.1:{};HGVSp={}:{}",
        at.strand, at.feature, written.notation, written.protein, written.description
    );
    if !written.offset.is_empty() {
        extra += &format!(";HGVS_OFFSET={}", written.offset);
    }
    default_text(&[format!(
        "past1\t{}\t{}\tENSG1\t{}\tTranscript\tframeshift_variant\t9-10\t9-10\t3-4\t-\t-\t-\t{extra}",
        at.location, at.allele, at.feature
    )])
}

/// The verdict of the two rows under the class beside the made reference and cache, and
/// what was applied.
fn applied(
    reference: Written,
    candidate: Written,
    at: Where,
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(&[EDGE_CLASS]);
    let report = compare_entries(
        &parse_default(&row(reference, at)).expect("the reference parses"),
        &parse_default(&row(candidate, at)).expect("the candidate parses"),
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

fn edges_beside(directory: &Path) -> Side {
    let edges = SequenceEdges::new(&write_reference(directory), &write_cache(directory))
        .expect("the reference and the cache read");
    Side {
        sequence_edges: Some(edges),
        ..Side::default()
    }
}

/// An insertion of two bases after the ninth base of a coding sequence, 38 bases before
/// its region's end on a forward-strand transcript, admits no shift: the base after it is
/// not the allele's first. The walk over the clipped slice's last bases reads the region's
/// first base, which is, and moves one base (TranscriptVariationAllele.pm:427-437). Ensembl
/// VEP so translates the turned bases after the tenth coding base, where the fourth residue
/// stands and the fifth changes, four residues before the new stop
/// (TranscriptVariationAllele.pm:2320-2365, :2394-2439, :1995-2015); the candidate
/// translates the bases after the ninth, where the fourth residue changes five before it.
/// The three fields of the row are the class, in a comparison and through the `field-class`
/// command, and so are those of the same insertion on the reverse-strand transcript of the
/// other strand, walked towards the region's first base.
#[test]
fn a_frameshift_protein_notation_is_the_class_where_each_walk_gives_it() {
    let directory = scratch("sequence-edges-protein");
    let side = edges_beside(&directory);
    let once = ("green".to_string(), count_map(&[(EDGE_CLASS, 1)]));
    assert_eq!(applied(REFERENCE, CANDIDATE, PAST, &side), once);
    let back = Where {
        location: "back:38-39",
        allele: "TC",
        strand: "-1",
        feature: "ENST00000000025",
    };
    let reverse = Written {
        offset: "-1",
        ..REFERENCE
    };
    assert_eq!(applied(reverse, CANDIDATE, back, &side), once);

    // The questions the golden corpus tool puts for the row: one a differing field, the
    // strand on the reference's row alone, and no shift length or position on either row.
    let cells = |written: Written, strand: &str| {
        let offset = if written.offset.is_empty() {
            "-"
        } else {
            written.offset
        };
        format!(
            "Transcript\t{strand}\t{}.1:{}\t{}:{}\t{offset}\t-\t-\t-\t-",
            PAST.feature, written.notation, written.protein, written.description
        )
    };
    let question = |name: &str, reference: &str, candidate: &str| {
        format!(
            "{}\t{}\t{}\t{name}\t{reference}\t{candidate}\t{}\t{}\n",
            PAST.location,
            PAST.allele,
            PAST.feature,
            cells(REFERENCE, PAST.strand),
            cells(CANDIDATE, "-")
        )
    };
    let protein = |written: Written| format!("{}:{}", written.protein, written.description);
    let questions = question(
        "HGVSc",
        "ENST00000000021.1:c.10_11insAG",
        "ENST00000000021.1:c.9_10insGA",
    ) + &question("HGVSp", &protein(REFERENCE), &protein(CANDIDATE))
        + &question("HGVS_OFFSET", "1", "-");
    assert!(questions.lines().all(|line| line.split('\t').count() == 24));
    let asked = planted(&directory, "questions.tsv", &questions);
    let answers = directory.join("answers.txt");
    let text = |path: &Path| path.to_str().expect("a path of text").to_string();
    let ran = support::run(
        &directory,
        &[
            "field-class".to_string(),
            format!("--class={EDGE_CLASS}"),
            "--reference-fasta".to_string(),
            text(&directory.join("reference.fa")),
            "--json-cache".to_string(),
            text(&directory.join("cache")),
            "--questions".to_string(),
            text(&asked),
            "--answers".to_string(),
            text(&answers),
        ],
    );
    assert_eq!(
        (ran.code, ran.out(), ran.err()),
        (0, String::new(), String::new())
    );
    assert_eq!(
        std::fs::read_to_string(&answers).expect("the answers are written"),
        "1\n1\n1\n"
    );
}

/// A deleted base is read where its walk leaves it, as an inserted one is. The bases at 19
/// and at 26 of `past`, the first of the sixth codon and the second of the eighth, are each
/// the region's first base and not the base after them: the sequence admits no shift, and
/// the walk over the clipped slice's last bases moves one base. Without its seventeenth
/// base the coding sequence reads alanine in the sixth codon as before and a stop in the
/// seventh, which Ensembl VEP writes alone in leucine's place
/// (TranscriptVariationAllele.pm:1997-1999); without its sixteenth it reads proline in the
/// sixth, two before the same stop. Without its twenty-fourth base it reads arginine in the
/// eighth codon, without its twenty-third threonine, and neither reading meets a stop
/// before the coding sequence ends, the cache stating no 3' UTR, so each count is `?`
/// (:2003-2006). The first deletion is the class on the reverse-strand transcript of the
/// other strand too.
#[test]
fn a_deletions_protein_notation_is_the_class_at_an_immediate_stop_and_without_a_stop() {
    let directory = scratch("sequence-edges-protein-deleted");
    let side = edges_beside(&directory);
    let once = ("green".to_string(), count_map(&[(EDGE_CLASS, 1)]));
    let deleted = |location: &'static str| Where {
        location,
        allele: "-",
        ..PAST
    };
    let stopped = Written {
        notation: "c.17del",
        description: "p.Leu7Ter",
        ..REFERENCE
    };
    let counted = Written {
        notation: "c.16del",
        description: "p.Ala6ProfsTer2",
        ..CANDIDATE
    };
    assert_eq!(applied(stopped, counted, deleted("past:19"), &side), once);
    let back = Where {
        location: "back:32",
        allele: "-",
        strand: "-1",
        feature: "ENST00000000025",
    };
    let reverse = Written {
        offset: "-1",
        ..stopped
    };
    assert_eq!(applied(reverse, counted, back, &side), once);
    let unstopped = Written {
        notation: "c.24del",
        description: "p.Ser8ArgfsTer?",
        ..REFERENCE
    };
    let uncounted = Written {
        notation: "c.23del",
        description: "p.Ser8ThrfsTer?",
        ..CANDIDATE
    };
    assert_eq!(
        applied(unstopped, uncounted, deleted("past:26"), &side),
        once
    );
}

/// The nearest rows that stay differences: a count to the new stop the changed sequence
/// does not give, on either side; the two sides' notations exchanged; another residue; two
/// protein names; the same three values 1,000 bases or more before the region's end, where
/// the cut flank is the allele's own; a transcript whose cache states no coding sequence;
/// one under another genetic code, whose changed sequence Ensembl VEP translates under the
/// standard one; and an insertion of three bases, which shifts no frame.
#[test]
fn a_protein_notation_the_walks_do_not_give_stays_a_difference() {
    let directory = scratch("sequence-edges-protein-refused");
    let side = edges_beside(&directory);
    let stays = |reference: Written, candidate: Written, at: Where| {
        assert_eq!(applied(reference, candidate, at, &side).0, "red");
    };
    let described = |written: Written<'static>, description: &'static str| Written {
        description,
        ..written
    };
    stays(REFERENCE, described(CANDIDATE, "p.Lys4GlufsTer4"), PAST);
    stays(described(REFERENCE, "p.Pro5AsnfsTer5"), CANDIDATE, PAST);
    stays(
        described(REFERENCE, CANDIDATE.description),
        described(CANDIDATE, REFERENCE.description),
        PAST,
    );
    stays(REFERENCE, described(CANDIDATE, "p.Lys4AspfsTer5"), PAST);
    stays(
        Written {
            protein: "ENSP00000000029.1",
            ..REFERENCE
        },
        CANDIDATE,
        PAST,
    );
    stays(
        REFERENCE,
        CANDIDATE,
        Where {
            location: "inner:12-13",
            feature: "ENST00000000022",
            ..PAST
        },
    );
    for feature in ["ENST00000000023", "ENST00000000024"] {
        stays(REFERENCE, CANDIDATE, Where { feature, ..PAST });
    }
    stays(
        Written {
            notation: "c.10_11insACG",
            description: "p.Lys4_Pro5insThr",
            ..REFERENCE
        },
        Written {
            notation: "c.9_10insGAC",
            description: "p.Ala3_Lys4insAsp",
            ..CANDIDATE
        },
        Where {
            allele: "GAC",
            ..PAST
        },
    );
}
