// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `stop_codon_overlap_insertion_range` beside `splice_region_variant`: the reference's
//! `coding_sequence_variant` against the candidate's `stop_retained_variant` on an insertion
//! between two neighbouring positions whose allele holds a letter outside `A`, `C`, `G` and
//! `T`, the rows `stop_codon_overlap_genomic_span` refuses. The predicate is asked on planted
//! entries, and two outputs are compared in the default, tab, VCF and JSON formats.

mod compare_common;
mod term_pairs_common;

use std::collections::BTreeMap;

use compare_common::*;
use term_pairs_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::rules::term_pairs::stop_codon_overlap_genomic_span_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_overlap_insertion_range_pair;
use vep_formats_compare::rules::term_pairs::stop_retained_without_stop_codon_test_pair;

const CLASS: &str = "stop_codon_overlap_insertion_range";
const GENOMIC_SPAN: &str = "stop_codon_overlap_genomic_span";
const UNTESTED: &str = "stop_retained_without_stop_codon_test";
const FEATURE: &str = "ENST00000399326";
/// The reference's terms of the row and the candidate's.
const UNKNOWN: &str = "splice_region_variant,coding_sequence_variant";
const KEPT: &str = "splice_region_variant,stop_retained_variant";
const INPUT: &str = "21\t35201927\t.\tG\tGAAN\t.\t.\t.";

/// The cells of one Transcript row that the class reads beside its terms.
#[derive(Clone, Copy)]
struct Row<'a> {
    /// The first column: an identifier, or the record's position and alleles.
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    cdna: &'a str,
    cds: &'a str,
    feature_type: &'a str,
    /// The Extra text after the IMPACT.
    extra: &'a str,
}

/// Three bases, the last unknown, between an intron's last base and the exon base that ends
/// a stop codon: the first of them is that exon base's letter, so the stop codon is kept.
const INSERTION: Row = Row {
    name: "21_35201927_G/GAAN",
    location: "21:35201927-35201928",
    allele: "AAN",
    cdna: "3048-3049",
    cds: "3047-3048",
    feature_type: "Transcript",
    extra: "STRAND=1",
};

impl Row<'_> {
    /// The row's line of a default-format output under `terms` and their `impact`.
    fn line(&self, terms: &str, impact: &str) -> String {
        format!(
            "{}\t{}\t{}\tENSG1\t{FEATURE}\t{}\t{terms}\t{}\t{}\t1016\t-\t-\t-\tIMPACT={impact};{}",
            self.name,
            self.location,
            self.allele,
            self.feature_type,
            self.cdna,
            self.cds,
            self.extra
        )
    }

    fn entry(&self, terms: &str) -> Entry {
        default_entry(&self.line(terms, "LOW"))
    }

    /// Whether the class's predicate takes the row under the reference's and the
    /// candidate's terms.
    fn taken(&self, reference: &str, candidate: &str) -> bool {
        holds(
            stop_codon_overlap_insertion_range_pair,
            &self.entry(reference),
            &self.entry(candidate),
        )
    }

    /// The verdict of the row's two default-format outputs under `classes` and what was
    /// applied.
    fn compared(
        &self,
        reference: (&str, &str),
        candidate: (&str, &str),
        classes: &[&str],
    ) -> (String, BTreeMap<String, i64>) {
        let (report, documented) = compare_default(
            &default_text(&[self.line(reference.0, reference.1)]),
            &default_text(&[self.line(candidate.0, candidate.1)]),
            classes,
        );
        (verdict(&report).to_string(), counts(&documented.applied))
    }
}

fn once() -> BTreeMap<String, i64> {
    count_map(&[(CLASS, 1)])
}

fn nothing() -> BTreeMap<String, i64> {
    count_map(&[])
}

/// The row's two entries in the tab, VCF and JSON formats, the reference's first.
fn in_format(fmt: &str, reference: &str, candidate: &str) -> (Vec<Entry>, Vec<Entry>) {
    let inputs = [INPUT.to_string()];
    let side = |terms: &str| match fmt {
        "tab" => {
            let header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND");
            let row = INSERTION
                .line(terms, "LOW")
                .replace("IMPACT=LOW;STRAND=1", "LOW\t1");
            parse_tab(&format!("{header}\n{row}\n")).expect("the output parses")
        }
        "vcf" => {
            let csq = [
                INSERTION.allele,
                &terms.replace(',', "&"),
                "LOW",
                FEATURE,
                "Transcript",
                INSERTION.cdna,
                INSERTION.cds,
            ]
            .join("|");
            let text = format!(
                "{}{}CSQ={csq}\n",
                csq_header(SPLICE_SPAN_CSQ_FIELDS),
                open_info(INPUT)
            );
            parse_vcf(&text, &inputs).expect("the output parses")
        }
        _ => {
            let held: Vec<&str> = terms.split(',').collect();
            let text = format!(
                "{{\"input\": {}, \"seq_region_name\": \"21\", \"start\": 35201928, \"end\": 35201927, \"allele_string\": \"-/AAN\", \"most_severe_consequence\": {}, \"transcript_consequences\": [{{\"variant_allele\": \"AAN\", \"transcript_id\": {}, \"consequence_terms\": {}, \"impact\": \"LOW\", \"cdna_start\": 3048, \"cdna_end\": 3049, \"cds_start\": 3047, \"cds_end\": 3048}}]}}\n",
                quoted(INPUT),
                quoted(most_severe(&held)),
                quoted(FEATURE),
                quoted_list(&held),
            );
            parse_json(&text, &inputs).expect("the output parses")
        }
    };
    (side(reference), side(candidate))
}

/// The verdict of the row in one of the tab, VCF and JSON formats under `classes` and what
/// was applied.
fn compared_in(
    fmt: &str,
    reference: &str,
    candidate: &str,
    classes: &[&str],
) -> (String, BTreeMap<String, i64>) {
    let (expected, actual) = in_format(fmt, reference, candidate);
    let mut documented = documented_with(classes);
    let report = compare(&expected, &actual, &mut documented);
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// The reference's coding_sequence_variant against the candidate's stop_retained_variant
/// beside splice_region_variant is the class's on an insertion whose allele holds an unknown
/// base, wherever the base stands, where the first column spells the record's alleles and
/// under an identifier, with the `/total` suffix of --total_length and beside a further term
/// both sides hold. The split codon's predicate refuses each of those rows, and takes the
/// same row of known bases, which this one refuses.
#[test]
fn the_first_pair_is_read_beside_splice_region_variant_where_the_allele_holds_an_unknown_base() {
    assert!(INSERTION.taken(UNKNOWN, KEPT));
    let rows = [
        INSERTION,
        Row {
            name: "rs1",
            ..INSERTION
        },
        Row {
            allele: "ANA",
            name: "rs1",
            ..INSERTION
        },
        Row {
            allele: "NAAAAT",
            name: "21_35201927_G/GNAAAAT",
            ..INSERTION
        },
        Row {
            cdna: "3048-3049/3600",
            cds: "3047-3048/3048",
            ..INSERTION
        },
        Row {
            cds: "3047-?",
            ..INSERTION
        },
    ];
    let nmd = |terms: &str| format!("{terms},NMD_transcript_variant");
    for row in rows {
        for (reference, candidate) in [
            (UNKNOWN.to_string(), KEPT.to_string()),
            (nmd(UNKNOWN), nmd(KEPT)),
        ] {
            let label = format!("{} {} {reference}", row.name, row.allele);
            assert!(row.taken(&reference, &candidate), "{label}");
            assert!(
                !holds(
                    stop_codon_overlap_genomic_span_pair,
                    &row.entry(&reference),
                    &row.entry(&candidate)
                ),
                "{label}"
            );
        }
    }
    for (name, allele) in [("21_35201927_G/GAAT", "AAT"), ("rs1", "AAA")] {
        let known = Row {
            name,
            allele,
            ..INSERTION
        };
        assert!(!known.taken(UNKNOWN, KEPT), "{allele}");
        assert!(
            holds(
                stop_codon_overlap_genomic_span_pair,
                &known.entry(UNKNOWN),
                &known.entry(KEPT)
            ),
            "{allele}"
        );
    }
}

/// Through the comparison the row is set aside once under the class in the default, tab, VCF
/// and JSON formats; without the class, and under the split codon's class, it is a
/// difference.
#[test]
fn the_row_with_an_unknown_base_is_set_aside_once_in_every_format() {
    let (reference, candidate) = ((UNKNOWN, "LOW"), (KEPT, "LOW"));
    assert_eq!(
        INSERTION.compared(reference, candidate, &[CLASS]),
        ("green".to_string(), once())
    );
    let named = Row {
        name: "rs1",
        ..INSERTION
    };
    assert_eq!(
        named.compared(reference, candidate, &[CLASS]),
        ("green".to_string(), once())
    );
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            compared_in(fmt, UNKNOWN, KEPT, &[CLASS]),
            ("green".to_string(), once()),
            "{fmt}"
        );
        assert_eq!(
            compared_in(fmt, UNKNOWN, KEPT, &[GENOMIC_SPAN]),
            ("red".to_string(), nothing()),
            "{fmt}"
        );
    }
    assert_eq!(
        INSERTION.compared(reference, candidate, &[]),
        ("red".to_string(), nothing())
    );
    assert_eq!(
        INSERTION.compared(reference, candidate, &[GENOMIC_SPAN]),
        ("red".to_string(), nothing())
    );
}

/// The sides exchanged, the reference's stop_retained_variant against the candidate's
/// coding_sequence_variant, are no row of the class in any format: that pair is the shape
/// of a stop codon test that is not made.
#[test]
fn the_row_with_an_unknown_base_stays_a_difference_with_the_sides_exchanged() {
    assert!(INSERTION.taken(UNKNOWN, KEPT));
    assert!(!INSERTION.taken(KEPT, UNKNOWN));
    assert!(holds(
        stop_retained_without_stop_codon_test_pair,
        &INSERTION.entry(KEPT),
        &INSERTION.entry(UNKNOWN)
    ));
    assert_eq!(
        INSERTION.compared((KEPT, "LOW"), (UNKNOWN, "LOW"), &[CLASS]),
        ("red".to_string(), nothing())
    );
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            compared_in(fmt, KEPT, UNKNOWN, &[CLASS]),
            ("red".to_string(), nothing()),
            "{fmt}"
        );
    }
    assert_eq!(
        INSERTION.compared((KEPT, "LOW"), (UNKNOWN, "LOW"), &[UNTESTED]),
        ("green".to_string(), count_map(&[(UNTESTED, 1)]))
    );
}

/// A row off the shape stays a difference: an allele holding a letter that is neither a base
/// nor `N`, a deletion, a row without a cDNA end or a CDS start, a Location of one position
/// or of three, a record whose reference allele keeps a base, another feature's row, and
/// terms that are not the pair.
#[test]
fn a_row_beside_splice_region_variant_off_the_shape_stays_a_difference() {
    assert!(INSERTION.taken(UNKNOWN, KEPT));
    let off = [
        Row {
            allele: "AAR",
            name: "rs1",
            ..INSERTION
        },
        Row {
            allele: "-",
            name: "rs1",
            location: "21:35201928-35201930",
            ..INSERTION
        },
        Row {
            cdna: "3048-?",
            cds: "3047-?",
            ..INSERTION
        },
        Row {
            cds: "?-3048",
            ..INSERTION
        },
        Row {
            cdna: "-",
            cds: "-",
            ..INSERTION
        },
        Row {
            location: "21:35201928",
            name: "rs1",
            ..INSERTION
        },
        Row {
            location: "21:35201927-35201929",
            name: "rs1",
            ..INSERTION
        },
        Row {
            name: "21_35201927_GA/GAAN",
            ..INSERTION
        },
        Row {
            feature_type: "RegulatoryFeature",
            ..INSERTION
        },
    ];
    for row in off {
        let label = format!("{} {} {}", row.name, row.location, row.cdna);
        assert!(!row.taken(UNKNOWN, KEPT), "{label}");
        assert_eq!(
            row.compared((UNKNOWN, "LOW"), (KEPT, "LOW"), &[CLASS]),
            ("red".to_string(), nothing()),
            "{label}"
        );
    }
    for (reference, candidate) in [
        (UNKNOWN, "splice_region_variant,stop_lost"),
        (UNKNOWN, "stop_retained_variant"),
        (
            UNKNOWN,
            "inframe_insertion,splice_region_variant,stop_retained_variant",
        ),
        ("splice_region_variant,stop_gained", KEPT),
        ("splice_region_variant", KEPT),
        (UNKNOWN, UNKNOWN),
    ] {
        assert!(
            !INSERTION.taken(reference, candidate),
            "{reference} {candidate}"
        );
    }
}

/// The rejoined rows of an allele that `--minimal` split from a record of two alternates
/// print the record's one-base Location and no transcript position in the default format, so
/// neither pair of the class reads them: the three added bases, and the two added bases
/// against the reference's frameshift_variant.
#[test]
fn the_rejoined_rows_of_a_split_allele_stay_differences_in_the_default_format() {
    assert!(INSERTION.taken(UNKNOWN, KEPT));
    let rejoined = Row {
        location: "21:35201927",
        cdna: "-",
        cds: "-",
        extra: "STRAND=1;MINIMISED=1",
        ..INSERTION
    };
    let stop_codon_classes = [CLASS, GENOMIC_SPAN, UNTESTED];
    let three = Row {
        name: "21_35201926_AG/AGAAN/AGNAG",
        ..rejoined
    };
    assert!(!three.taken(UNKNOWN, KEPT));
    assert_eq!(
        three.compared((UNKNOWN, "LOW"), (KEPT, "LOW"), &stop_codon_classes),
        ("red".to_string(), nothing())
    );
    let two = Row {
        name: "21_35201926_AG/AGA/AGAN",
        allele: "AN",
        ..rejoined
    };
    let shifted = "frameshift_variant,splice_region_variant";
    assert!(!two.taken(shifted, KEPT));
    assert_eq!(
        two.compared((shifted, "HIGH"), (KEPT, "LOW"), &stop_codon_classes),
        ("red".to_string(), nothing())
    );
}
