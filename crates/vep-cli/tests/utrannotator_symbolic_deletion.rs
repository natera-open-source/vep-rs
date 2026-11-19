// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--plugin UTRAnnotator` through the `vep` binary on symbolic structural
//! alleles in a 5' UTR. A symbolic deletion whose stated span lies inside one
//! 5' UTR exon is annotated as the sequence deletion of that span: its row
//! carries the five keys of the same deletion written with explicit alleles.
//! A plain symbolic insertion whose stated length is no multiple of three, at
//! a point inside a reading frame of the 5' UTR, carries the frameshift with
//! the pairs its length and the reference UTR decide. Every other symbolic
//! allele carries the three counts and `-`. Ensembl VEP reads the allele's
//! text as the reference allele (UTRAnnotator.pm:151-153) and removes as many
//! bases of the UTR as the text has characters (UTRAnnotator.pm:1189-1199),
//! whatever span the record states.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

mod common;

const GAINED_START: &str = "5_prime_UTR_premature_start_codon_gain_variant";
const LOST_STOP: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";
const LOST_START: &str = "5_prime_UTR_premature_start_codon_loss_variant";
const FRAMESHIFT: &str = "5_prime_UTR_uORF_frameshift_variant";
const NONE: &str = "-";

/// A forward-strand transcript of the corpus's cache whose 5' UTR has four
/// exons: 21:10521553-10521694, 10524580-10524688, 10527355-10527412 and
/// 10538681-10538723.
const FORWARD: &str = "ENST00000427445";
/// A reverse-strand transcript whose 5' UTR is 21:14210846-14210955.
const REVERSE: &str = "ENST00000400211";

/// The corpus whose cache carries both transcripts, with the reference and the
/// uORF evidence table the plugin reads.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-plugin-utr-annotator")
}

/// The plugin's five keys of one row.
type Keys = BTreeMap<String, String>;

/// Runs the binary on `records`, whole VCF lines, with the plugin and returns
/// the record name and the plugin's keys of every row on `feature`, in output
/// order.
fn annotate(records: &[&str], feature: &str) -> Vec<(String, Keys)> {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in records {
        text.push_str(record);
        text.push('\n');
    }
    std::fs::write(&input, text).unwrap();
    let output = tmp.path().join("out.txt");
    let corpus = corpus();
    let run = common::isolated_vep_anywhere()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .arg("--plugin")
        .arg(format!(
            "UTRAnnotator,file={}",
            corpus
                .join("data/uORF_5UTR_GRCh38_PUBLIC.chr21.txt")
                .display()
        ))
        .output()
        .expect("vep binary runs");
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|columns| columns[4] == feature)
        .map(|columns| {
            let keys = columns[13]
                .split(';')
                .filter_map(|pair| pair.split_once('='))
                .filter(|(key, _)| key.starts_with("5UTR_") || key.starts_with("Existing_"))
                .map(|(key, value)| (key.to_owned(), value.to_owned()))
                .collect();
            (columns[0].to_owned(), keys)
        })
        .collect()
}

/// The keys of the one row of the record named `name`.
fn row_of<'a>(rows: &'a [(String, Keys)], name: &str) -> &'a Keys {
    let mut found = rows.iter().filter(|(record, _)| record == name);
    let (_, keys) = found.next().unwrap_or_else(|| panic!("no row of {name}"));
    assert!(found.next().is_none(), "more than one row of {name}");
    keys
}

/// A symbolic deletion inside one 5' UTR exon, the same deletion written with
/// explicit alleles, the transcript, and the term and the annotation both rows
/// carry there.
struct Deletion {
    symbolic: &'static str,
    explicit: &'static str,
    feature: &'static str,
    term: &'static str,
    annotation: &'static str,
}

/// On each strand: a deletion inside a uORF whose length is no multiple of
/// three, one that removes a uORF's stop codon, one that opens a start codon,
/// one that brings a start codon together, and one that touches neither.
const DELETIONS: [Deletion; 10] = [
    Deletion {
        symbolic: "21\t10524588\tfwd_frameshift\tN\t<DEL>\t.\t.\tEND=10524628",
        explicit: "21\t10524588\tfwd_frameshift_explicit\tGGACCGACGACACAAGACCTCAGACTTGTGTTATTCTAGCA\tG\t.\t.\t.",
        feature: FORWARD,
        term: FRAMESHIFT,
        annotation: "Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:alt_type=uORF:alt_type_length=66:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51",
    },
    Deletion {
        symbolic: "21\t10524619\tfwd_stop_loss\tN\t<DEL>\t.\t.\tEND=10524631",
        explicit: "21\t10524619\tfwd_stop_loss_explicit\tTATTCTAGCAGCT\tT\t.\t.\t.",
        feature: FORWARD,
        term: LOST_STOP,
        annotation: "AltStop=True:AltStopDistanceToCDS=159:Evidence=False:FrameWithCDS=inFrame:KozakContext=AGAATGT:KozakStrength=Moderate",
    },
    Deletion {
        symbolic: "21\t10524582\tfwd_start_loss\tN\t<DEL>\t.\t.\tEND=10524584",
        explicit: "21\t10524582\tfwd_start_loss_explicit\tAAT\tA\t.\t.\t.",
        feature: FORWARD,
        term: LOST_START,
        annotation: "CapDistanceToStart=145:DistanceToCDS=207:DistanceToStop=48:Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:type=uORF",
    },
    Deletion {
        symbolic: "21\t10521583\tfwd_start_gain\tN\t<DEL>\t.\t.\tEND=10521590",
        explicit: "21\t10521583\tfwd_start_gain_explicit\tACAACAGC\tA\t.\t.\t.",
        feature: FORWARD,
        term: GAINED_START,
        annotation: "CapDistanceToStart=30:DistanceToCDS=315:DistanceToStop=66:KozakContext=GGCATGC:KozakStrength=Moderate:type=uORF",
    },
    Deletion {
        symbolic: "21\t10521552\tfwd_none\tN\t<DEL>\t.\t.\tEND=10521555",
        explicit: "21\t10521552\tfwd_none_explicit\tGCGT\tG\t.\t.\t.",
        feature: FORWARD,
        term: NONE,
        annotation: NONE,
    },
    Deletion {
        symbolic: "21\t14210850\trev_frameshift\tN\t<DEL>\t.\t.\tEND=14210851",
        explicit: "21\t14210850\trev_frameshift_explicit\tAA\tA\t.\t.\t.",
        feature: REVERSE,
        term: FRAMESHIFT,
        annotation: "Evidence=False:KozakContext=GTGATGA:KozakStrength=Moderate:alt_type=OutOfFrame_oORF:alt_type_length=90:ref_StartDistanceToCDS=87:ref_type=InFrame_oORF:ref_type_length=1215",
    },
    Deletion {
        symbolic: "21\t14210863\trev_stop_loss\tN\t<DEL>\t.\t.\tEND=14210875",
        explicit: "21\t14210863\trev_stop_loss_explicit\tCAAAAACAGCACT\tC\t.\t.\t.",
        feature: REVERSE,
        term: LOST_STOP,
        annotation: "AltStop=False:AltStopDistanceToCDS=NA:Evidence=False:FrameWithCDS=outOfFrame:KozakContext=TCAATGA:KozakStrength=Weak",
    },
    Deletion {
        symbolic: "21\t14210893\trev_start_loss\tN\t<DEL>\t.\t.\tEND=14210905",
        explicit: "21\t14210893\trev_start_loss_explicit\tAATTCCGTAACTC\tA\t.\t.\t.",
        feature: REVERSE,
        term: LOST_START,
        annotation: "CapDistanceToStart=48:DistanceToCDS=62:DistanceToStop=30:Evidence=False:KozakContext=TCAATGA:KozakStrength=Weak:type=uORF",
    },
    Deletion {
        symbolic: "21\t14210870\trev_start_gain\tN\t<DEL>\t.\t.\tEND=14210879",
        explicit: "21\t14210870\trev_start_gain_explicit\tAGCACTCAAT\tA\t.\t.\t.",
        feature: REVERSE,
        term: GAINED_START,
        annotation: "CapDistanceToStart=75:DistanceToCDS=26:DistanceToStop=27:KozakContext=GTGATGT:KozakStrength=Moderate:type=OutOfFrame_oORF",
    },
    Deletion {
        symbolic: "21\t14210845\trev_none\tN\t<DEL>\t.\t.\tEND=14210848",
        explicit: "21\t14210845\trev_none_explicit\tTTTG\tT\t.\t.\t.",
        feature: REVERSE,
        term: NONE,
        annotation: NONE,
    },
];

/// The symbolic record's row carries the term and the annotation of the
/// explicit record's row, and the two rows agree on all five keys.
#[test]
fn a_symbolic_deletion_inside_one_utr_exon_carries_the_row_of_its_explicit_form() {
    for feature in [FORWARD, REVERSE] {
        let cases: Vec<&Deletion> = DELETIONS
            .iter()
            .filter(|case| case.feature == feature)
            .collect();
        let records: Vec<&str> = cases
            .iter()
            .flat_map(|case| [case.symbolic, case.explicit])
            .collect();
        let rows = annotate(&records, feature);
        for case in cases {
            let name = case.symbolic.split('\t').nth(2).unwrap();
            let symbolic = row_of(&rows, name);
            let explicit = row_of(&rows, &format!("{name}_explicit"));
            assert_eq!(symbolic["5UTR_consequence"], case.term, "{name}");
            assert_eq!(symbolic["5UTR_annotation"], case.annotation, "{name}");
            assert_eq!(symbolic.len(), 5, "{name}: {symbolic:?}");
            assert_eq!(symbolic, explicit, "{name}");
        }
    }
}

/// A mobile-element deletion, a copy-number allele of no copies, a deletion
/// that states its length in place of its end and one whose confidence
/// intervals have no width state the same span as the 40-base deletion and
/// carry its row.
#[test]
fn a_deletion_in_another_text_of_the_same_span_carries_the_same_row() {
    let frameshift = &DELETIONS[0];
    let records = [
        frameshift.explicit,
        "21\t10524588\tmobile_element\tN\t<DEL:ME:ALU>\t.\t.\tEND=10524628",
        "21\t10524588\tcopy_number_loss\tN\t<CN0>\t.\t.\tEND=10524628",
        "21\t10524588\tlength_without_end\tN\t<DEL>\t.\t.\tSVLEN=-40",
        "21\t10524588\tinterval_of_no_width\tN\t<DEL>\t.\t.\tEND=10524628;CIPOS=0,0;CIEND=0,0",
    ];
    let rows = annotate(&records, FORWARD);
    let explicit = row_of(&rows, "fwd_frameshift_explicit");
    assert_eq!(explicit["5UTR_consequence"], FRAMESHIFT);
    for name in [
        "mobile_element",
        "copy_number_loss",
        "length_without_end",
        "interval_of_no_width",
    ] {
        assert_eq!(row_of(&rows, name), explicit, "{name}");
    }
}

/// A deletion whose record carries the `IMPRECISE` flag states no exact span,
/// with or without a confidence interval: its row carries the counts and `-`
/// where the same record without the flag carries the frameshift.
#[test]
fn a_deletion_marked_imprecise_with_no_interval_carries_the_counts_alone() {
    let records = [
        "21\t10524588\tprecise\tN\t<DEL>\t.\t.\tEND=10524628",
        "21\t10524588\timprecise\tN\t<DEL>\t.\t.\tEND=10524628;IMPRECISE",
    ];
    let rows = annotate(&records, FORWARD);
    assert_eq!(row_of(&rows, "precise")["5UTR_consequence"], FRAMESHIFT);
    let imprecise = row_of(&rows, "imprecise");
    assert_eq!(imprecise["5UTR_consequence"], NONE);
    assert_eq!(imprecise["5UTR_annotation"], NONE);
    assert_eq!(imprecise["Existing_uORFs"], "6");
    assert_eq!(imprecise["Existing_OutOfFrame_oORFs"], "0");
    assert_eq!(imprecise["Existing_InFrame_oORFs"], "0");
    assert_eq!(imprecise.len(), 5, "{imprecise:?}");
}

/// A plain symbolic insertion inside a reading frame of the 5' UTR, the same
/// insertion written with as many explicit bases that bring no start and no
/// stop codon, the transcript, and the annotation the symbolic row carries.
struct Insertion {
    symbolic: &'static str,
    explicit: &'static str,
    feature: &'static str,
    annotation: &'static str,
}

/// On each strand: a point inside one frame that a stop codon past it ends
/// whatever is inserted, one inside one frame that no such stop codon ends,
/// and points inside two or more frames.
const INSERTIONS: [Insertion; 8] = [
    Insertion {
        symbolic: "21\t10524588\tfwd_one_frame_decided\tN\t<INS>\t.\t.\tEND=10524588;SVLEN=4",
        explicit: "21\t10524588\tfwd_one_frame_decided_explicit\tN\tNCCCC\t.\t.\t.",
        feature: FORWARD,
        annotation: "Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51",
    },
    Insertion {
        symbolic: "21\t10527361\tfwd_one_frame_open\tN\t<INS>\t.\t.\tEND=10527361;SVLEN=2",
        explicit: "21\t10527361\tfwd_one_frame_open_explicit\tN\tNCC\t.\t.\t.",
        feature: FORWARD,
        annotation: "Evidence=False:KozakContext=GTTATGG:KozakStrength=Strong:ref_StartDistanceToCDS=97:ref_type=uORF:ref_type_length=36",
    },
    Insertion {
        symbolic: "21\t10527387\tfwd_two_frames\tN\t<INS>\t.\t.\tEND=10527387;SVLEN=1",
        explicit: "21\t10527387\tfwd_two_frames_explicit\tN\tNC\t.\t.\t.",
        feature: FORWARD,
        annotation: "Evidence=False:KozakContext=GTTATGG:KozakStrength=Strong:alt_type=uORF:ref_StartDistanceToCDS=97:ref_type=uORF:ref_type_length=36,Evidence=False:KozakContext=TGTATGC:KozakStrength=Weak:ref_StartDistanceToCDS=71:ref_type=uORF:ref_type_length=57",
    },
    Insertion {
        symbolic: "21\t10538694\tfwd_four_frames\tN\t<INS>\t.\t.\tEND=10538694;SVLEN=1",
        explicit: "21\t10538694\tfwd_four_frames_explicit\tN\tNC\t.\t.\t.",
        feature: FORWARD,
        annotation: "Evidence=False:KozakContext=TGTATGC:KozakStrength=Weak:ref_StartDistanceToCDS=71:ref_type=uORF:ref_type_length=57,Evidence=False:KozakContext=ACTATGA:KozakStrength=Moderate:ref_StartDistanceToCDS=56:ref_type=uORF:ref_type_length=42,Evidence=False:KozakContext=ACAATGG:KozakStrength=Strong:ref_StartDistanceToCDS=47:ref_type=uORF:ref_type_length=33,Evidence=False:KozakContext=CAAATGA:KozakStrength=Weak:ref_StartDistanceToCDS=32:ref_type=uORF:ref_type_length=18",
    },
    Insertion {
        symbolic: "21\t14210908\trev_one_frame_decided\tN\t<INS>\t.\t.\tEND=14210908;SVLEN=1",
        explicit: "21\t14210908\trev_one_frame_decided_explicit\tN\tNG\t.\t.\t.",
        feature: REVERSE,
        annotation: "Evidence=False:KozakContext=GTGATGA:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=87:ref_type=InFrame_oORF:ref_type_length=1215",
    },
    Insertion {
        symbolic: "21\t14210846\trev_one_frame_open\tN\t<INS>\t.\t.\tEND=14210846;SVLEN=1",
        explicit: "21\t14210846\trev_one_frame_open_explicit\tN\tNG\t.\t.\t.",
        feature: REVERSE,
        annotation: "Evidence=False:KozakContext=GTGATGA:KozakStrength=Moderate:ref_StartDistanceToCDS=87:ref_type=InFrame_oORF:ref_type_length=1215",
    },
    Insertion {
        symbolic: "21\t14210878\trev_two_frames\tN\t<INS>\t.\t.\tEND=14210878;SVLEN=2",
        explicit: "21\t14210878\trev_two_frames_explicit\tN\tNGG\t.\t.\t.",
        feature: REVERSE,
        annotation: "Evidence=False:KozakContext=GTGATGA:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=87:ref_type=InFrame_oORF:ref_type_length=1215,Evidence=False:KozakContext=TCAATGA:KozakStrength=Weak:ref_StartDistanceToCDS=62:ref_type=uORF:ref_type_length=33",
    },
    Insertion {
        symbolic: "21\t14210882\trev_two_frames_decided\tN\t<INS>\t.\t.\tEND=14210882;SVLEN=2",
        explicit: "21\t14210882\trev_two_frames_decided_explicit\tN\tNGG\t.\t.\t.",
        feature: REVERSE,
        annotation: "Evidence=False:KozakContext=GTGATGA:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=87:ref_type=InFrame_oORF:ref_type_length=1215,Evidence=False:KozakContext=TCAATGA:KozakStrength=Weak:alt_type=uORF:ref_StartDistanceToCDS=62:ref_type=uORF:ref_type_length=33",
    },
];

/// The pairs of an explicit insertion's annotation that its length and the
/// reference UTR decide: every pair but `alt_type_length` and an `alt_type`
/// other than `uORF`, which the inserted bases decide.
fn decided_by_length(annotation: &str) -> String {
    annotation
        .split(',')
        .map(|effect| {
            effect
                .split(':')
                .filter(|pair| {
                    !pair.starts_with("alt_type_length=")
                        && (!pair.starts_with("alt_type=") || *pair == "alt_type=uORF")
                })
                .collect::<Vec<_>>()
                .join(":")
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// A plain symbolic insertion whose stated length is no multiple of three, at
/// a point inside a reading frame of the 5' UTR, carries the frameshift
/// (UTRAnnotator.pm:889, :921-923). Its annotation holds the pairs of the
/// explicit insertion's row that the length and the reference UTR decide, and
/// the three counts are the explicit row's.
#[test]
fn a_symbolic_insertion_inside_a_reading_frame_carries_the_frameshift_of_its_length() {
    for feature in [FORWARD, REVERSE] {
        let cases: Vec<&Insertion> = INSERTIONS
            .iter()
            .filter(|case| case.feature == feature)
            .collect();
        let records: Vec<&str> = cases
            .iter()
            .flat_map(|case| [case.symbolic, case.explicit])
            .collect();
        let rows = annotate(&records, feature);
        for case in cases {
            let name = case.symbolic.split('\t').nth(2).unwrap();
            let symbolic = row_of(&rows, name);
            let explicit = row_of(&rows, &format!("{name}_explicit"));
            assert_eq!(symbolic["5UTR_consequence"], FRAMESHIFT, "{name}");
            assert_eq!(symbolic["5UTR_annotation"], case.annotation, "{name}");
            assert_eq!(explicit["5UTR_consequence"], FRAMESHIFT, "{name}");
            assert_eq!(
                symbolic["5UTR_annotation"],
                decided_by_length(&explicit["5UTR_annotation"]),
                "{name}"
            );
            assert!(
                !symbolic["5UTR_annotation"].contains("alt_type_length"),
                "{name}"
            );
            for count in [
                "Existing_uORFs",
                "Existing_OutOfFrame_oORFs",
                "Existing_InFrame_oORFs",
            ] {
                assert_eq!(symbolic[count], explicit[count], "{name}");
            }
            assert_eq!(symbolic.len(), 5, "{name}: {symbolic:?}");
        }
    }
}

/// A record that writes no end, or its position plus its length as its end,
/// states the same insertion as one that writes its position as its end: the
/// three rows are one.
#[test]
fn an_insertion_with_no_written_end_or_its_span_s_end_carries_the_same_row() {
    let records = [
        "21\t10524588\tposition_as_end\tN\t<INS>\t.\t.\tEND=10524588;SVLEN=4",
        "21\t10524588\tno_end\tN\t<INS>\t.\t.\tSVLEN=4",
        "21\t10524588\tspan_as_end\tN\t<INS>\t.\t.\tEND=10524592;SVLEN=4",
    ];
    let rows = annotate(&records, FORWARD);
    let position_as_end = row_of(&rows, "position_as_end");
    assert_eq!(position_as_end["5UTR_consequence"], FRAMESHIFT);
    assert_eq!(row_of(&rows, "no_end"), position_as_end);
    assert_eq!(row_of(&rows, "span_as_end"), position_as_end);
}

/// Every other symbolic allele carries the three counts of the reference 5' UTR
/// and `-`: a deletion whose record states a confidence interval around its
/// position and its end, one that runs from one UTR exon into the next, into
/// the intron or in from before the transcript, each alternate of a record with
/// two, a duplication, a tandem duplication, a copy-number gain and an
/// inversion. So does an insertion whose length is a multiple of three, whose
/// point lies in no reading frame, whose record writes another end, no length,
/// a confidence interval, the `IMPRECISE` flag or a second alternate, and one
/// whose text names a mobile element or another kind of insertion.
#[test]
fn every_other_symbolic_allele_carries_the_counts_alone() {
    let records = [
        "21\t10524588\tinterval\tN\t<DEL>\t.\t.\tEND=10524628;CIPOS=-5,5;CIEND=-5,5",
        "21\t10524588\timprecise_with_interval\tN\t<DEL>\t.\t.\tEND=10524628;IMPRECISE;CIPOS=-5,5;CIEND=-5,5",
        "21\t10524680\tinto_next_part\tN\t<DEL>\t.\t.\tEND=10527360",
        "21\t10524680\tinto_intron\tN\t<DEL>\t.\t.\tEND=10524700",
        "21\t10521540\tfrom_before_transcript\tN\t<DEL>\t.\t.\tEND=10521560",
        "21\t10524588\ttwo_alternates\tN\t<CN0>,<CN2>\t.\t.\tEND=10524628",
        "21\t10524588\tduplication\tN\t<DUP>\t.\t.\tEND=10524628;SVLEN=40",
        "21\t10524588\ttandem_duplication\tN\t<DUP:TANDEM>\t.\t.\tEND=10524628",
        "21\t10524588\tcopy_number_gain\tN\t<CN2>\t.\t.\tEND=10524628",
        "21\t10524588\tinversion\tN\t<INV>\t.\t.\tEND=10524628",
        "21\t10524588\tinsertion_of_three\tN\t<INS>\t.\t.\tEND=10524588;SVLEN=3",
        "21\t10524582\tinsertion_outside_frames\tN\t<INS>\t.\t.\tEND=10524582;SVLEN=4",
        "21\t10524588\tinsertion_other_end\tN\t<INS>\t.\t.\tEND=10524590;SVLEN=4",
        "21\t10524588\tinsertion_without_length\tN\t<INS>\t.\t.\tEND=10524588",
        "21\t10524588\tinsertion_with_interval\tN\t<INS>\t.\t.\tEND=10524588;SVLEN=4;CIPOS=-5,5;CIEND=-5,5",
        "21\t10524588\tinsertion_marked_imprecise\tN\t<INS>\t.\t.\tEND=10524588;SVLEN=4;IMPRECISE",
        "21\t10524588\tinsertion_two_alternates\tN\t<INS>,<DUP>\t.\t.\tEND=10524588;SVLEN=4",
        "21\t10524588\tmobile_insertion\tN\t<INS:ME:ALU>\t.\t.\tEND=10524588;SVLEN=4",
        "21\t10524588\tother_insertion_text\tN\t<INS:MT>\t.\t.\tEND=10524588;SVLEN=4",
    ];
    let rows = annotate(&records, FORWARD);
    let names: Vec<&str> = rows.iter().map(|(name, _)| name.as_str()).collect();
    assert_eq!(names.len(), 21, "{names:?}");
    for two in ["two_alternates", "insertion_two_alternates"] {
        assert_eq!(names.iter().filter(|name| **name == two).count(), 2);
    }
    for (name, keys) in &rows {
        assert_eq!(keys["5UTR_consequence"], NONE, "{name}");
        assert_eq!(keys["5UTR_annotation"], NONE, "{name}");
        assert_eq!(keys["Existing_uORFs"], "6", "{name}");
        assert_eq!(keys["Existing_OutOfFrame_oORFs"], "0", "{name}");
        assert_eq!(keys["Existing_InFrame_oORFs"], "0", "{name}");
        assert_eq!(keys.len(), 5, "{name}: {keys:?}");
    }
}
