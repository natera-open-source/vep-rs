// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pairs of entries the term-pair classes take and the cells that follow two sets of
//! terms, on planted rows of each format.

mod compare_common;
mod term_pairs_common;

use std::collections::BTreeSet;

use term_pairs_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::coding_ends::terms::CodingEndClass;
use vep_formats_compare::consts::CONSEQUENCE_RANK;
use vep_formats_compare::consts::TERM_PAIR_FOLLOWING_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::rules::following::MISSENSE_TERM_FIELDS;
use vep_formats_compare::rules::following::PROTEIN_CHANGE_TERM_FIELDS;
use vep_formats_compare::rules::frameshift_at_a_held_stop::frameshift_inside_a_held_stop_codon_pair;
use vep_formats_compare::rules::settled_stop::stop_codon_settled_before_unknown_bases_pair;
use vep_formats_compare::rules::term_pairs::splice_term_past_the_record_span_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_overlap_genomic_span_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_overlap_insertion_range_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_removed_without_stop_codon_test_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_test_coordinates;
use vep_formats_compare::rules::term_pairs::stop_codon_test_not_made;
use vep_formats_compare::rules::term_pairs::stop_retained_three_base_stop_codon_test_pair;
use vep_formats_compare::rules::term_pairs::stop_retained_without_stop_codon_test_pair;
use vep_formats_compare::rules::term_pairs::term_pair_following_cell;
use vep_formats_compare::rules::term_pairs::terms_impact;
use vep_formats_compare::rules::term_pairs::DOWNSTREAM_FIELDS;
use vep_formats_compare::rules::term_pairs::STOP_CODON_PLUGIN_VALUE_PAIRS;
use vep_formats_compare::rules::term_pairs::TERM_IMPACT;
use vep_formats_compare::rules::term_pairs::TERM_PAIR_CLASSES;
use vep_formats_compare::rules::term_pairs::TRUNCATING_TERM_FIELDS;

/// A row of `stop_codon_row_at` as the entry the default format gives it.
fn stop_entry(allele: &str, terms: &str, cdna: &str, cds: &str, location: &str) -> Entry {
    default_entry(&stop_codon_row_at(allele, terms, cdna, cds, location))
}

fn acceptor(term: &str) -> String {
    format!("splice_acceptor_variant,{term}")
}

const STOP_INPUT: &str = "19\t13208045\t1993671\tGACCTG\tG\t.\t.\t.";

/// The reference's stop_retained_variant against nothing, coding_sequence_variant or
/// frameshift_variant on the candidate's side is the shape on an insertion's or deletion's
/// Transcript row that shows the stop codon test is not made: a cDNA position missing at one
/// end or whole, no CDS start, or an allele holding a letter that names no base. A position
/// that is missing reads the same with the `/total` suffix of --total_length after its `-`.
#[test]
fn stop_retained_without_stop_codon_test_pair_reads_a_row_whose_test_was_not_made() {
    let pair: Pair = stop_retained_without_stop_codon_test_pair;
    let row =
        |terms: &str, cdna: &str, cds: &str| stop_entry("-", terms, cdna, cds, STOP_CODON_LOCATION);
    for (cdna, cds) in [
        ("?-7057", "-"),
        ("-", "-"),
        ("1705-?", "1630-?"),
        ("?-7057", "?-6980"),
        ("120-7100", "-"),
        ("-/7140", "-/7023"),
        ("120-7100/7140", "-/7023"),
    ] {
        let reference = row(&acceptor("stop_retained_variant"), cdna, cds);
        for term in ["coding_sequence_variant", "frameshift_variant"] {
            assert!(
                holds(pair, &reference, &row(&acceptor(term), cdna, cds)),
                "{cdna} {cds} {term}"
            );
        }
        assert!(
            holds(pair, &reference, &row("splice_acceptor_variant", cdna, cds)),
            "{cdna} {cds}"
        );
    }
    let inserted = |letters: &str, terms: &str| {
        stop_entry(
            letters,
            terms,
            "7038-7039",
            "6961-6962",
            "7:152648662-152648663",
        )
    };
    for letters in ["ACNT", "ACRT", "N"] {
        assert!(
            holds(
                pair,
                &inserted(letters, "stop_retained_variant"),
                &inserted(letters, "frameshift_variant")
            ),
            "{letters}"
        );
    }
    let tested = row(&acceptor("stop_retained_variant"), "7055-7059", "6978-6982");
    assert!(!holds(
        pair,
        &tested,
        &row(
            &acceptor("coding_sequence_variant"),
            "7055-7059",
            "6978-6982"
        )
    ));
    let tested_into_the_utr = row(
        "stop_retained_variant,3_prime_UTR_variant",
        "7055-7059",
        "6978-?",
    );
    assert!(!holds(
        pair,
        &tested_into_the_utr,
        &row(
            "coding_sequence_variant,3_prime_UTR_variant",
            "7055-7059",
            "6978-?"
        )
    ));
    assert!(!holds(
        pair,
        &inserted("ACGT", "stop_retained_variant"),
        &inserted("ACGT", "frameshift_variant")
    ));
    let reference = row(&acceptor("stop_retained_variant"), "?-7057", "-");
    for terms in [
        "stop_lost",
        "coding_sequence_variant,frameshift_variant",
        "stop_retained_variant",
        "inframe_deletion",
        "stop_retained_variant,coding_sequence_variant",
    ] {
        assert!(
            !holds(pair, &reference, &row(&acceptor(terms), "?-7057", "-")),
            "{terms}"
        );
    }
    assert!(!holds(
        pair,
        &reference,
        &row("coding_sequence_variant", "?-7057", "-")
    ));
    assert!(!holds(
        pair,
        &row(&acceptor("coding_sequence_variant"), "?-7057", "-"),
        &reference
    ));
    assert!(!holds(
        pair,
        &reference,
        &row(
            &acceptor("coding_sequence_variant"),
            "7055-7059",
            "6978-6982"
        )
    ));
    let substitution = |terms: &str| stop_entry("A", terms, "-", "-", "19:13208046");
    assert!(!holds(
        pair,
        &substitution("stop_retained_variant"),
        &substitution("coding_sequence_variant")
    ));
    let regulatory = |terms: &str| {
        default_entry(&stop_codon_row_of(
            "-",
            terms,
            "-",
            "-",
            STOP_CODON_LOCATION,
            "RegulatoryFeature",
        ))
    };
    assert!(!holds(
        pair,
        &regulatory("stop_retained_variant"),
        &regulatory("coding_sequence_variant")
    ));
}

/// The same shape on the entries of the VCF and JSON formats, which write a missing
/// position as an empty field and as a key left out.
#[test]
fn stop_retained_without_stop_codon_test_pair_reads_the_entries_of_the_vcf_and_json_formats() {
    let pair: Pair = stop_retained_without_stop_codon_test_pair;
    let inputs = [STOP_INPUT.to_string()];
    let csq_entry = |term: &str, cdna: &str, cds: &str| {
        let text = format!(
            "{}{}CSQ=-|splice_acceptor_variant&{term}|ENST00000635895|Transcript|{cdna}|{cds}\n",
            csq_header("Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position"),
            open_info(STOP_INPUT)
        );
        vcf_feature_entries(&text, &inputs).remove(0)
    };
    assert!(holds(
        pair,
        &csq_entry("stop_retained_variant", "?-7057", ""),
        &csq_entry("coding_sequence_variant", "?-7057", "")
    ));
    assert!(!holds(
        pair,
        &csq_entry("stop_retained_variant", "7055-7059", "6978-6982"),
        &csq_entry("coding_sequence_variant", "7055-7059", "6978-6982")
    ));

    let json_entry = |term: &str, positions: &[(&str, i64)]| {
        let positions: String = positions
            .iter()
            .map(|(name, position)| format!(", \"{name}\": {position}"))
            .collect();
        let text = format!(
            "{{\"input\": {}, \"id\": \"1993671\", \"allele_string\": \"ACCTG/-\", \"most_severe_consequence\": \"splice_acceptor_variant\", \"transcript_consequences\": [{{\"variant_allele\": \"-\", \"transcript_id\": \"ENST00000635895\", \"consequence_terms\": [\"splice_acceptor_variant\", {}]{positions}}}]}}\n",
            quoted(STOP_INPUT),
            quoted(term)
        );
        json_feature_entries(&text, &inputs).remove(0)
    };
    let untested: [&[(&str, i64)]; 3] = [
        &[("cdna_end", 7057)],
        &[],
        &[("cdna_start", 120), ("cdna_end", 7100)],
    ];
    for positions in untested {
        assert!(
            holds(
                pair,
                &json_entry("stop_retained_variant", positions),
                &json_entry("coding_sequence_variant", positions)
            ),
            "{positions:?}"
        );
    }
    let whole = [
        ("cdna_start", 7055),
        ("cdna_end", 7059),
        ("cds_start", 6978),
        ("cds_end", 6982),
    ];
    assert!(!holds(
        pair,
        &json_entry("stop_retained_variant", &whole),
        &json_entry("coding_sequence_variant", &whole)
    ));
}

/// The reference's coding_sequence_variant against the candidate's stop_retained_variant
/// is the shape on an insertion's or deletion's Transcript row with both cDNA positions, a
/// CDS start and splice_region_variant on both sides.
#[test]
fn stop_codon_overlap_genomic_span_pair_reads_a_tested_row_at_an_exon_end() {
    let pair: Pair = stop_codon_overlap_genomic_span_pair;
    let terms = |term: &str| format!("splice_region_variant,{term},NMD_transcript_variant");
    let row = |term: &str, cdna: &str, cds: &str, allele: &str, location: &str| {
        stop_entry(allele, &terms(term), cdna, cds, location)
    };
    let plain = |term: &str| row(term, "1705", "1630", "-", "3:36993663");
    let positioned = |term: &str, cdna: &str, cds: &str| row(term, cdna, cds, "-", "3:36993663");
    assert!(holds(
        pair,
        &plain("coding_sequence_variant"),
        &plain("stop_retained_variant")
    ));
    let inserted = |term: &str| row(term, "1704-1705", "1629-1630", "A", "3:36993662-36993663");
    assert!(holds(
        pair,
        &inserted("coding_sequence_variant"),
        &inserted("stop_retained_variant")
    ));
    assert!(!holds(
        pair,
        &plain("stop_retained_variant"),
        &plain("coding_sequence_variant")
    ));
    for (cdna, cds, taken) in [
        ("1705-?", "1630-?", false),
        ("1705", "-", false),
        ("1705/1800", "1630/1653", true),
        ("1705/1800", "-/1653", false),
    ] {
        assert_eq!(
            holds(
                pair,
                &positioned("coding_sequence_variant", cdna, cds),
                &positioned("stop_retained_variant", cdna, cds)
            ),
            taken,
            "{cdna} {cds}"
        );
    }
    assert!(!holds(
        pair,
        &plain("coding_sequence_variant"),
        &plain("stop_retained_variant,stop_lost")
    ));
    assert!(!holds(
        pair,
        &plain("coding_sequence_variant"),
        &plain("stop_lost")
    ));
    let substitution = |term: &str| row(term, "1705", "1630", "A", "3:36993663");
    assert!(!holds(
        pair,
        &substitution("coding_sequence_variant"),
        &substitution("stop_retained_variant")
    ));
    let bare = |term: &str| {
        stop_entry(
            "-",
            &format!("{term},NMD_transcript_variant"),
            "1705",
            "1630",
            "3:36993663",
        )
    };
    assert!(!holds(
        pair,
        &bare("coding_sequence_variant"),
        &bare("stop_retained_variant")
    ));
}

/// The columns of a row of the three-base stop codon test, the deletion's by default.
#[derive(Clone, Copy)]
struct Tested<'a> {
    cdna: &'a str,
    cds: &'a str,
    allele: &'a str,
    location: &'a str,
    feature_type: &'a str,
}

impl Tested<'_> {
    fn row(&self, terms: &str) -> Entry {
        default_entry(&stop_codon_row_of(
            self.allele,
            terms,
            self.cdna,
            self.cds,
            self.location,
            self.feature_type,
        ))
    }
}

/// On a row whose first column spells a reference allele and two alternates, none of them
/// the row's allele once the shared first base is dropped, the Location decides whether the
/// row changes the length, as under an identifier; an alternate that is the row's allele
/// decides before it.
#[test]
fn stop_retained_three_base_stop_codon_test_pair_reads_an_allele_no_alternate_names_by_its_location(
) {
    let pair: Pair = stop_retained_three_base_stop_codon_test_pair;
    let taken = |location: &str, allele: &str| {
        let row = |terms: &str| {
            let row = stop_codon_row_at(allele, terms, "7054-7056", "6977-6979", location);
            default_entry(&row.replacen("1993671", "19_13208044_ACGTT/ACA/AG", 1))
        };
        holds(
            pair,
            &row("stop_retained_variant"),
            &row("coding_sequence_variant"),
        )
    };
    assert!(taken("19:13208046-13208048", "A"));
    assert!(!taken("19:13208046", "A"));
    assert!(taken("19:13208045", "G"));
}

/// The reference's stop_retained_variant against coding_sequence_variant, or, alone or
/// beside one of inframe_insertion, inframe_deletion and protein_altering_variant, against
/// frameshift_variant with or without stop_gained in the place of those terms, is the shape
/// on an insertion's or deletion's Transcript row that shows the stop codon test is made:
/// both cDNA positions, a CDS start and an allele of bases.
#[test]
fn stop_retained_three_base_stop_codon_test_pair_reads_a_tested_row() {
    let pair: Pair = stop_retained_three_base_stop_codon_test_pair;
    let utr = |terms: &str| format!("{terms},3_prime_UTR_variant");
    let deletion = Tested {
        cdna: "2558-2564",
        cds: "2249-?",
        allele: "-",
        location: "2:190970701-190970707",
        feature_type: "Transcript",
    };
    let at = |cdna, cds| Tested {
        cdna,
        cds,
        ..deletion
    };
    let coding = at("961-967", "264-270");
    assert!(holds(
        pair,
        &deletion.row(&utr("stop_retained_variant")),
        &deletion.row(&utr("coding_sequence_variant"))
    ));
    assert!(holds(
        pair,
        &coding.row("stop_retained_variant"),
        &coding.row("frameshift_variant")
    ));
    assert!(holds(
        pair,
        &coding.row("inframe_deletion,stop_retained_variant"),
        &coding.row("frameshift_variant")
    ));
    let totals = at("2558-2564/2700", "2249-?/2253");
    assert!(holds(
        pair,
        &totals.row("stop_retained_variant"),
        &totals.row("coding_sequence_variant")
    ));
    let inserted = Tested {
        cdna: "1213-1214",
        cds: "1100-1101",
        allele: "ACGTACGT",
        location: "7:152648662-152648663",
        feature_type: "Transcript",
    };
    for reference in [
        "inframe_insertion,stop_retained_variant",
        "protein_altering_variant,stop_retained_variant",
        "stop_retained_variant",
    ] {
        for candidate in ["stop_gained,frameshift_variant", "frameshift_variant"] {
            assert!(
                holds(pair, &inserted.row(reference), &inserted.row(candidate)),
                "{reference} {candidate}"
            );
        }
    }
    for (cdna, cds) in [
        ("2558-?", "2249-?"),
        ("?-2564", "-"),
        ("-", "-"),
        ("2558-2564", "-"),
        ("2558-?/2700", "2249-?/2253"),
    ] {
        let untested = at(cdna, cds);
        assert!(
            !holds(
                pair,
                &untested.row("stop_retained_variant"),
                &untested.row("coding_sequence_variant")
            ),
            "{cdna} {cds}"
        );
    }
    let lettered = Tested {
        allele: "ACGTACGN",
        ..inserted
    };
    assert!(!holds(
        pair,
        &lettered.row("stop_retained_variant"),
        &lettered.row("frameshift_variant")
    ));
    let reference = deletion.row(&utr("stop_retained_variant"));
    for terms in [
        "stop_retained_variant",
        "stop_lost",
        "coding_sequence_variant,frameshift_variant",
        "stop_gained",
        "inframe_deletion",
        "stop_retained_variant,coding_sequence_variant",
        "stop_lost,frameshift_variant",
    ] {
        assert!(
            !holds(pair, &reference, &deletion.row(&utr(terms))),
            "{terms}"
        );
    }
    assert!(!holds(
        pair,
        &reference,
        &deletion.row("coding_sequence_variant")
    ));
    assert!(!holds(
        pair,
        &reference,
        &deletion.row("3_prime_UTR_variant")
    ));
    assert!(!holds(
        pair,
        &deletion.row(&utr("coding_sequence_variant")),
        &reference
    ));
    assert!(!holds(
        pair,
        &coding.row("inframe_deletion,stop_retained_variant"),
        &coding.row("coding_sequence_variant")
    ));
    assert!(!holds(
        pair,
        &coding.row("missense_variant,stop_retained_variant"),
        &coding.row("frameshift_variant")
    ));
    let substitution = Tested {
        allele: "A",
        location: "2:190970701",
        ..deletion
    };
    assert!(!holds(
        pair,
        &substitution.row("stop_retained_variant"),
        &substitution.row("coding_sequence_variant")
    ));
    let regulatory = Tested {
        feature_type: "RegulatoryFeature",
        ..deletion
    };
    assert!(!holds(
        pair,
        &regulatory.row("stop_retained_variant"),
        &regulatory.row("coding_sequence_variant")
    ));
    assert!(!holds(
        stop_retained_without_stop_codon_test_pair,
        &reference,
        &deletion.row(&utr("coding_sequence_variant"))
    ));
}

/// The reference's stop_retained_variant against the candidate's stop_lost is the shape on
/// an insertion's or deletion's Transcript row that lacks a cDNA position at one end or
/// whole, or a CDS start; the row the first stop codon class takes is refused.
#[test]
fn stop_codon_removed_without_stop_codon_test_pair_reads_a_row_without_its_coordinates() {
    let pair: Pair = stop_codon_removed_without_stop_codon_test_pair;
    let untested = Tested {
        cdna: "?-7057",
        cds: "-",
        allele: "-",
        location: STOP_CODON_LOCATION,
        feature_type: "Transcript",
    };
    let at = |cdna, cds| Tested {
        cdna,
        cds,
        ..untested
    };
    for (cdna, cds) in [
        ("?-7057", "-"),
        ("-", "-"),
        ("1705-?", "1630-?"),
        ("120-7100", "-"),
        ("-/7140", "-/7023"),
        ("?-7057/7140", "?-6980/7023"),
    ] {
        let columns = at(cdna, cds);
        assert!(
            holds(
                pair,
                &columns.row(&acceptor("stop_retained_variant")),
                &columns.row(&acceptor("stop_lost"))
            ),
            "{cdna} {cds}"
        );
    }
    let without_a_start = at("120-7100", "-");
    assert!(holds(
        pair,
        &without_a_start.row("start_lost,stop_retained_variant"),
        &without_a_start.row("stop_lost,start_lost")
    ));
    let replaced = Tested {
        allele: "A",
        ..at("1705-?", "1630-?")
    };
    assert!(holds(
        pair,
        &replaced.row("stop_retained_variant"),
        &replaced.row("stop_lost")
    ));
    for (cdna, cds) in [
        ("7055-7059", "6978-6982"),
        ("7055-7059", "6978-?"),
        ("7055-7059/7140", "6978-?/7023"),
    ] {
        let columns = at(cdna, cds);
        assert!(
            !holds(
                pair,
                &columns.row(&acceptor("stop_retained_variant")),
                &columns.row(&acceptor("stop_lost"))
            ),
            "{cdna} {cds}"
        );
    }
    let reference = untested.row(&acceptor("stop_retained_variant"));
    for terms in [
        "coding_sequence_variant",
        "stop_retained_variant",
        "stop_lost,frameshift_variant",
        "stop_lost,coding_sequence_variant",
        "frameshift_variant",
    ] {
        assert!(
            !holds(pair, &reference, &untested.row(&acceptor(terms))),
            "{terms}"
        );
    }
    assert!(!holds(
        pair,
        &reference,
        &untested.row("splice_acceptor_variant")
    ));
    assert!(!holds(pair, &reference, &untested.row("stop_lost")));
    assert!(!holds(
        pair,
        &untested.row(&acceptor("stop_lost")),
        &reference
    ));
    let substitution = Tested {
        allele: "A",
        location: "19:13208046",
        ..untested
    };
    assert!(!holds(
        pair,
        &substitution.row("stop_retained_variant"),
        &substitution.row("stop_lost")
    ));
    let regulatory = Tested {
        feature_type: "RegulatoryFeature",
        ..untested
    };
    assert!(!holds(
        pair,
        &regulatory.row("stop_retained_variant"),
        &regulatory.row("stop_lost")
    ));
    assert!(holds(
        stop_retained_without_stop_codon_test_pair,
        &reference,
        &untested.row(&acceptor("coding_sequence_variant"))
    ));
    assert!(!holds(
        stop_retained_without_stop_codon_test_pair,
        &reference,
        &untested.row(&acceptor("stop_lost"))
    ));
}

/// The reference's coding_sequence_variant against the candidate's stop_retained_variant, or
/// its frameshift_variant with or without stop_gained against stop_retained_variant with or
/// without inframe_insertion, is the shape on an insertion's Transcript row with both cDNA
/// positions and a CDS start on both sides, an allele holding N included; the split
/// codon's and the three-base predicates take their own.
#[test]
fn stop_codon_overlap_insertion_range_pair_reads_an_insertion_that_keeps_the_stop_codon() {
    let pair: Pair = stop_codon_overlap_insertion_range_pair;
    let insertion = Tested {
        cdna: "7054-7055",
        cds: "6977-6978",
        allele: "TTGA",
        location: "19:13208046-13208047",
        feature_type: "Transcript",
    };
    for (reference, candidate) in [
        (
            "inframe_insertion,coding_sequence_variant",
            "inframe_insertion,stop_retained_variant",
        ),
        (
            "coding_sequence_variant,NMD_transcript_variant",
            "stop_retained_variant,NMD_transcript_variant",
        ),
        (
            "stop_gained,frameshift_variant",
            "inframe_insertion,stop_retained_variant",
        ),
        (
            "frameshift_variant",
            "inframe_insertion,stop_retained_variant",
        ),
        (
            "frameshift_variant,NMD_transcript_variant",
            "stop_retained_variant,NMD_transcript_variant",
        ),
        ("stop_gained,frameshift_variant", "stop_retained_variant"),
    ] {
        let takes = |columns: Tested| holds(pair, &columns.row(reference), &columns.row(candidate));
        assert!(takes(insertion), "{reference} {candidate}");
        assert!(
            takes(Tested {
                allele: "TNGA",
                ..insertion
            }),
            "{reference} {candidate}"
        );
        assert!(
            takes(Tested {
                cdna: "7054-7055/7140",
                cds: "6977-6978/7023",
                ..insertion
            }),
            "{reference} {candidate}"
        );
        assert!(
            takes(Tested {
                cds: "6978-?",
                ..insertion
            }),
            "{reference} {candidate}"
        );
        assert!(
            !holds(pair, &insertion.row(candidate), &insertion.row(reference)),
            "{reference} {candidate}"
        );
        assert!(
            !takes(Tested {
                cdna: "7055-7058",
                cds: "6978-6981",
                allele: "-",
                location: "19:13208046-13208049",
                ..insertion
            }),
            "{reference} {candidate}"
        );
        for (cdna, cds) in [
            ("7054-?", "6977-?"),
            ("?-7055", "?-6978"),
            ("7054-7055", "-"),
            ("7054-7055", "?-6978"),
            ("-", "-"),
            ("7054-7055/7140", "?-6978/7023"),
        ] {
            assert!(
                !takes(Tested {
                    cdna,
                    cds,
                    ..insertion
                }),
                "{reference} {candidate} {cdna} {cds}"
            );
        }
        assert!(
            !takes(Tested {
                allele: "TRGA",
                ..insertion
            }),
            "{reference} {candidate}"
        );
        assert!(
            !takes(Tested {
                feature_type: "RegulatoryFeature",
                ..insertion
            }),
            "{reference} {candidate}"
        );
    }
    let row = |terms: &str| insertion.row(terms);
    let region = |term: &str| {
        row(&format!(
            "inframe_insertion,splice_region_variant,{term},NMD_transcript_variant"
        ))
    };
    let split = (
        region("coding_sequence_variant"),
        region("stop_retained_variant"),
    );
    assert!(!holds(pair, &split.0, &split.1));
    assert!(holds(
        stop_codon_overlap_genomic_span_pair,
        &split.0,
        &split.1
    ));
    assert!(holds(
        pair,
        &row("splice_region_variant,frameshift_variant"),
        &row("inframe_insertion,splice_region_variant,stop_retained_variant")
    ));
    let three_base = (
        row("inframe_insertion,stop_retained_variant"),
        row("stop_gained,frameshift_variant"),
    );
    assert!(!holds(pair, &three_base.0, &three_base.1));
    assert!(holds(
        stop_retained_three_base_stop_codon_test_pair,
        &three_base.0,
        &three_base.1
    ));
    assert!(!holds(
        stop_retained_three_base_stop_codon_test_pair,
        &three_base.1,
        &three_base.0
    ));
    for (reference, candidate) in [
        (
            "inframe_insertion,coding_sequence_variant",
            "inframe_insertion,stop_retained_variant,stop_lost",
        ),
        (
            "coding_sequence_variant",
            "inframe_insertion,stop_retained_variant",
        ),
        (
            "frameshift_variant",
            "protein_altering_variant,stop_retained_variant",
        ),
        (
            "missense_variant,frameshift_variant",
            "stop_retained_variant",
        ),
        ("stop_gained", "stop_retained_variant"),
    ] {
        assert!(
            !holds(pair, &row(reference), &row(candidate)),
            "{reference} {candidate}"
        );
    }
}

/// Ensembl VEP leaves its stop codon test out for a span without a cDNA position at one end
/// or a CDS start. A deletion from the 5' UTR to the stop codon's last base prints both cDNA
/// positions and the CDS end alone: its test is not made, so the removed-codon and the
/// untested predicates take the row and the tested ones refuse it. A row with the CDS start
/// and no CDS end is one whose test is made.
#[test]
fn the_stop_codon_test_is_left_out_for_want_of_the_cds_start() {
    let row = |terms: &str, cds: &str, cdna: &str| {
        stop_entry("-", terms, cdna, cds, "21:40000021-40000099")
    };
    for (cdna, cds) in [("21-99", "?-69"), ("21-99/105", "?-69/69")] {
        let entry = row("start_lost,stop_retained_variant", cds, cdna);
        assert!(!stop_codon_test_coordinates(&entry), "{cds}");
        assert!(stop_codon_test_not_made(&entry), "{cds}");
        assert!(
            holds(
                stop_codon_removed_without_stop_codon_test_pair,
                &entry,
                &row("stop_lost,start_lost", cds, cdna)
            ),
            "{cds}"
        );
        assert!(
            holds(
                stop_retained_without_stop_codon_test_pair,
                &entry,
                &row("start_lost", cds, cdna)
            ),
            "{cds}"
        );
        let tested = (
            row("stop_retained_variant", cds, cdna),
            row("coding_sequence_variant", cds, cdna),
        );
        assert!(
            !holds(
                stop_retained_three_base_stop_codon_test_pair,
                &tested.0,
                &tested.1
            ),
            "{cds}"
        );
        assert!(
            holds(
                stop_retained_without_stop_codon_test_pair,
                &tested.0,
                &tested.1
            ),
            "{cds}"
        );
        let region = (
            row("splice_region_variant,coding_sequence_variant", cds, cdna),
            row("splice_region_variant,stop_retained_variant", cds, cdna),
        );
        assert!(
            !holds(stop_codon_overlap_genomic_span_pair, &region.0, &region.1),
            "{cds}"
        );
    }
    for cds in ["60-?", "60-69", "69", "60-?/69"] {
        let entry = row("start_lost,stop_retained_variant", cds, "21-99");
        assert!(stop_codon_test_coordinates(&entry), "{cds}");
        assert!(!stop_codon_test_not_made(&entry), "{cds}");
        assert!(
            !holds(
                stop_codon_removed_without_stop_codon_test_pair,
                &entry,
                &row("stop_lost,start_lost", cds, "21-99")
            ),
            "{cds}"
        );
        assert!(
            !holds(
                stop_retained_without_stop_codon_test_pair,
                &entry,
                &row("start_lost", cds, "21-99")
            ),
            "{cds}"
        );
        assert!(
            holds(
                stop_retained_three_base_stop_codon_test_pair,
                &row("stop_retained_variant", cds, "21-99"),
                &row("coding_sequence_variant", cds, "21-99")
            ),
            "{cds}"
        );
    }
}

/// A deletion's Transcript row with these terms and this Extra text.
fn following_entry(consequence: &str, extra: &str) -> Entry {
    default_entry(&format!(
        "1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{consequence}\t7055-7059\t6978-?\t-\t-\t-\t-\t{extra}"
    ))
}

/// `term_pair_following_cell` on the field `name` of two rows of one key, a value of
/// `None` being a field the row leaves out.
fn following(
    name: &str,
    reference_terms: &str,
    reference_value: Option<&str>,
    candidate_terms: &str,
    candidate_value: Option<&str>,
) -> bool {
    let entry = |terms: &str, value: Option<&str>| {
        let extra = match value {
            None => "STRAND=-1".to_string(),
            Some(value) => format!("{name}={value};STRAND=-1"),
        };
        following_entry(terms, &extra)
    };
    let e = entry(reference_terms, reference_value);
    let a = entry(candidate_terms, candidate_value);
    term_pair_following_cell(&e, &a, name, field(&e, name), field(&a, name))
}

fn spliced(term: &str) -> String {
    format!("splice_acceptor_variant,{term},intron_variant")
}

/// Each term with a rank has an impact, and a row's IMPACT is its most severe term's.
#[test]
fn term_impact_names_every_ranked_term() {
    let ranked: BTreeSet<&str> = CONSEQUENCE_RANK.iter().map(|(term, _)| *term).collect();
    let with_impact: BTreeSet<&str> = TERM_IMPACT.iter().map(|(term, _)| *term).collect();
    assert_eq!(with_impact, ranked);
    assert_eq!(TERM_IMPACT.len(), with_impact.len());
    let impacts: BTreeSet<&str> = TERM_IMPACT.iter().map(|(_, impact)| *impact).collect();
    assert_eq!(
        impacts,
        BTreeSet::from(["HIGH", "MODERATE", "LOW", "MODIFIER"])
    );
    assert_eq!(
        terms_impact(&["3_prime_UTR_variant", "stop_retained_variant"]).as_deref(),
        Some("LOW")
    );
    assert_eq!(
        terms_impact(&[
            "NMD_transcript_variant",
            "frameshift_variant",
            "stop_gained"
        ])
        .as_deref(),
        Some("HIGH")
    );
    assert_eq!(terms_impact(&["stop_retained_variant", "not_a_term"]), None);
    assert_eq!(terms_impact(&[]), None);
}

/// The IMPACT cell follows where each side carries the impact of its own most severe term;
/// an impact that is not its side's terms' and a row whose terms agree stay differences.
#[test]
fn term_pair_following_cell_takes_the_impact_of_each_sides_terms() {
    assert!(following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        UNKNOWN,
        Some("MODIFIER")
    ));
    assert!(following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        SHIFTED,
        Some("HIGH")
    ));
    assert!(following(
        "IMPACT",
        "inframe_deletion,stop_retained_variant",
        Some("MODERATE"),
        "frameshift_variant",
        Some("HIGH")
    ));
    assert!(following(
        "IMPACT",
        UNKNOWN,
        Some("MODIFIER"),
        RETAINED,
        Some("LOW")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        UNKNOWN,
        Some("HIGH")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("MODERATE"),
        UNKNOWN,
        Some("MODIFIER")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        RETAINED,
        Some("MODIFIER")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        "not_a_term,3_prime_UTR_variant",
        Some("MODIFIER")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        None,
        UNKNOWN,
        Some("MODIFIER")
    ));
}

/// DownstreamProtein follows where it stands on the side whose terms hold frameshift_variant
/// and no splice term, and is left out on the other.
#[test]
fn term_pair_following_cell_takes_downstream_protein_where_frameshift_comes() {
    let name = "DownstreamProtein";
    assert!(following(name, RETAINED, None, SHIFTED, Some("MKV")));
    assert!(following(name, SHIFTED, Some("MKV"), RETAINED, None));
    assert!(!following(name, RETAINED, None, UNKNOWN, Some("MKV")));
    assert!(!following(name, RETAINED, Some("MKV"), SHIFTED, None));
    assert!(!following(
        name,
        RETAINED,
        Some("MKL"),
        SHIFTED,
        Some("MKV")
    ));
    assert!(!following(
        name,
        &spliced("stop_retained_variant"),
        None,
        &spliced("frameshift_variant"),
        Some("MKV")
    ));
    assert!(!following(name, SHIFTED, Some("MKL"), SHIFTED, Some("MKV")));
}

/// ProteinLengthChange follows as DownstreamProtein does, by the same plugin's terms.
#[test]
fn term_pair_following_cell_takes_protein_length_change_where_frameshift_comes() {
    let name = "ProteinLengthChange";
    assert!(following(name, RETAINED, None, SHIFTED, Some("-12")));
    assert!(following(name, SHIFTED, Some("-12"), RETAINED, None));
    assert!(!following(name, RETAINED, None, UNKNOWN, Some("-12")));
    assert!(!following(name, RETAINED, Some("4"), SHIFTED, Some("-12")));
    assert!(!following(
        name,
        &spliced("stop_retained_variant"),
        None,
        &spliced("frameshift_variant"),
        Some("-12")
    ));
}

/// NMD follows where it stands on the side whose terms hold stop_gained, frameshift_variant
/// or a splice site term and is left out on the side whose terms hold none.
#[test]
fn term_pair_following_cell_takes_nmd_where_a_truncating_term_comes() {
    let value = Some("NMD_escaping_variant");
    assert!(following("NMD", RETAINED, None, SHIFTED, value));
    assert!(following(
        "NMD",
        RETAINED,
        None,
        "stop_gained,frameshift_variant",
        value
    ));
    assert!(following("NMD", SHIFTED, value, RETAINED, None));
    let (kept, unknown) = (
        spliced("stop_retained_variant"),
        spliced("coding_sequence_variant"),
    );
    assert!(!following("NMD", &kept, value, &unknown, None));
    assert!(!following("NMD", &kept, None, &unknown, value));
    assert!(!following("NMD", RETAINED, None, UNKNOWN, value));
    assert!(!following("NMD", RETAINED, value, SHIFTED, None));
}

/// LoF follows as NMD does; two confidence values on a row both of whose term sets allow
/// LOFTEE to run stay a difference.
#[test]
fn term_pair_following_cell_takes_lof_where_a_truncating_term_comes() {
    assert!(following("LoF", RETAINED, None, SHIFTED, Some("HC")));
    assert!(following("LoF", SHIFTED, Some("LC"), RETAINED, None));
    assert!(!following(
        "LoF",
        SHIFTED,
        Some("HC"),
        "stop_gained,frameshift_variant",
        Some("LC")
    ));
    assert!(!following("LoF", RETAINED, None, UNKNOWN, Some("HC")));
    assert!(!following("LoF", RETAINED, Some("HC"), SHIFTED, None));
}

/// LoF_info follows as LoF does; a text both sides carry stays a difference.
#[test]
fn term_pair_following_cell_takes_lof_info_where_a_truncating_term_comes() {
    let name = "LoF_info";
    assert!(following(
        name,
        RETAINED,
        None,
        SHIFTED,
        Some("PERCENTILE:0.97,GERP_DIST:12.5")
    ));
    assert!(!following(
        name,
        SHIFTED,
        Some("PERCENTILE:0.97"),
        "stop_gained,frameshift_variant",
        Some("PERCENTILE:0.96")
    ));
    assert!(!following(
        name,
        RETAINED,
        None,
        UNKNOWN,
        Some("PERCENTILE:0.97")
    ));
    assert!(!following(
        name,
        RETAINED,
        Some("PERCENTILE:0.97"),
        UNKNOWN,
        None
    ));
}

/// LoF_filter follows as LoF does; a list both sides carry stays a difference.
#[test]
fn term_pair_following_cell_takes_lof_filter_where_a_truncating_term_comes() {
    let name = "LoF_filter";
    assert!(following(name, RETAINED, None, SHIFTED, Some("END_TRUNC")));
    assert!(!following(
        name,
        SHIFTED,
        Some("END_TRUNC"),
        "stop_gained,frameshift_variant",
        Some("END_TRUNC,GC_TO_GT_DONOR")
    ));
    assert!(!following(name, RETAINED, None, UNKNOWN, Some("END_TRUNC")));
}

/// LoF_flags follows as LoF does, on whichever side's terms allow LOFTEE to run.
#[test]
fn term_pair_following_cell_takes_lof_flags_where_a_truncating_term_comes() {
    let name = "LoF_flags";
    assert!(following(
        name,
        RETAINED,
        None,
        SHIFTED,
        Some("SINGLE_EXON")
    ));
    assert!(following(
        name,
        "stop_gained,frameshift_variant",
        Some("SINGLE_EXON,PHYLOCSF_WEAK"),
        RETAINED,
        None
    ));
    assert!(!following(
        name,
        RETAINED,
        Some("SINGLE_EXON"),
        SHIFTED,
        None
    ));
    assert!(!following(
        name,
        RETAINED,
        None,
        UNKNOWN,
        Some("SINGLE_EXON")
    ));
    assert!(!following(
        name,
        SHIFTED,
        Some("SINGLE_EXON"),
        "stop_gained,frameshift_variant",
        Some("SINGLE_EXON,PHYLOCSF_WEAK")
    ));
}

const PROTEIN: &str = "ENSP00000490161.1:p.";

/// HGVSp follows where the side holding stop_retained_variant writes a notation that is no
/// frameshift's and the side without it a frameshift's, `fs` or the `ext` of one beginning
/// at the stop codon, both naming one protein.
#[test]
fn term_pair_following_cell_takes_a_frameshift_notation_beside_a_lost_retained_stop() {
    let notation = |change: &str| format!("{PROTEIN}{change}");
    let hgvsp = |reference_terms: &str,
                 reference: Option<String>,
                 candidate_terms: &str,
                 candidate: Option<String>| {
        following(
            "HGVSp",
            reference_terms,
            reference.as_deref(),
            candidate_terms,
            candidate.as_deref(),
        )
    };
    assert!(hgvsp(
        RETAINED,
        Some(notation("Tyr22_Ter23delinsCys")),
        UNKNOWN,
        Some(notation("Tyr22CysfsTer4"))
    ));
    assert!(hgvsp(
        RETAINED,
        Some(notation("Ter23=")),
        UNKNOWN,
        Some(notation("Ter23TyrextTer4"))
    ));
    assert!(hgvsp(
        RETAINED,
        Some(notation("Ter23%3D")),
        SHIFTED,
        Some(notation("Ter23TyrextTer?"))
    ));
    assert!(hgvsp(
        UNKNOWN,
        Some(notation("Tyr22CysfsTer4")),
        RETAINED,
        Some(notation("Tyr22_Ter23delinsCys"))
    ));
    assert!(!hgvsp(
        RETAINED,
        Some(notation("Tyr22_Ter23delinsCys")),
        UNKNOWN,
        Some(notation("Tyr22_Ter23delinsSer"))
    ));
    assert!(!hgvsp(
        RETAINED,
        Some(notation("Tyr22CysfsTer4")),
        UNKNOWN,
        Some(notation("Tyr22CysfsTer9"))
    ));
    assert!(!hgvsp(
        RETAINED,
        Some(notation("Tyr22_Ter23delinsCys")),
        UNKNOWN,
        Some("ENSP00000490162.1:p.Tyr22CysfsTer4".to_string())
    ));
    assert!(!hgvsp(
        RETAINED,
        None,
        UNKNOWN,
        Some(notation("Tyr22CysfsTer4"))
    ));
    assert!(!hgvsp(
        UNKNOWN,
        Some(notation("Tyr22_Ter23delinsCys")),
        SHIFTED,
        Some(notation("Tyr22CysfsTer4"))
    ));
    assert!(!following(
        "HGVSc",
        RETAINED,
        Some("ENST00000635895.1:c.6978_*2del"),
        UNKNOWN,
        Some("ENST00000635895.1:c.6979_*3del")
    ));
}

/// The cells of a JSON entry follow under the keys the JSON format writes.
#[test]
fn term_pair_following_cell_reads_the_keys_of_a_json_entry() {
    let inputs = [STOP_INPUT.to_string()];
    let json_entry = |term: &str, fields: &[(&str, &str)]| {
        let fields: String = fields
            .iter()
            .map(|(name, value)| format!(", {}: {}", quoted(name), quoted(value)))
            .collect();
        let text = format!(
            "{{\"input\": {}, \"id\": \"1993671\", \"allele_string\": \"ACCTG/-\", \"most_severe_consequence\": {term}, \"transcript_consequences\": [{{\"variant_allele\": \"-\", \"transcript_id\": \"ENST00000635895\", \"consequence_terms\": [{term}, \"3_prime_UTR_variant\"]{fields}}}]}}\n",
            quoted(STOP_INPUT),
            term = quoted(term)
        );
        json_feature_entries(&text, &inputs).remove(0)
    };
    let kept = format!("{PROTEIN}Ter23=");
    let extended = format!("{PROTEIN}Ter23TyrextTer4");
    let e = json_entry(
        "stop_retained_variant",
        &[("impact", "LOW"), ("hgvsp", &kept)],
    );
    let a = json_entry(
        "frameshift_variant",
        &[
            ("impact", "HIGH"),
            ("hgvsp", &extended),
            ("nmd", "NMD_escaping_variant"),
        ],
    );
    for name in ["impact", "hgvsp", "nmd"] {
        assert!(
            term_pair_following_cell(&e, &a, name, field(&e, name), field(&a, name)),
            "{name}"
        );
    }
    assert!(!term_pair_following_cell(
        &e,
        &a,
        "variant_allele",
        "-",
        "A"
    ));
}

/// A row of each of three stop codon pairs the plugin-value class reads: its
/// cDNA_position, CDS_position, Allele and Location.
const UNTESTED_DELETION: [&str; 4] = ["1705-?", "1630-?", "-", "19:13208046-13208050"];
const TESTED_DELETION: [&str; 4] = ["961-967", "264-270", "-", "2:190970701-190970707"];
const KEPT_STOP_INSERTION: [&str; 4] = ["7054-7055", "6977-6978", "TTGA", "19:13208046-13208047"];

/// Whether the plugin-value class accepts the field `name` of two rows of one key at
/// `span`, each side its terms and its value, `None` being a field the row leaves out.
fn plugin_value_taken(
    name: &str,
    span: [&str; 4],
    reference: (&str, Option<&str>),
    candidate: (&str, Option<&str>),
) -> bool {
    let [cdna, cds, allele, location] = span;
    let entry = |(terms, value): (&str, Option<&str>)| {
        let extra = match value {
            None => "STRAND=-1".to_string(),
            Some(value) => format!("{name}={value};STRAND=-1"),
        };
        default_entry(&format!(
            "1993671\t{location}\t{allele}\tENSG1\tENST00000635895\tTranscript\t{terms}\t{cdna}\t{cds}\t-\t-\t-\t-\t{extra}"
        ))
    };
    let (e, a) = (entry(reference), entry(candidate));
    class_accepts(
        PLUGIN_VALUE_CLASS,
        &e,
        name,
        field(&e, name),
        field(&a, name),
        Some(&a),
    )
}

const PLUGIN_VALUES: [(&str, &str); 7] = [
    ("DownstreamProtein", "MKV"),
    ("ProteinLengthChange", "-12"),
    ("NMD", "NMD_escaping_variant"),
    ("LoF", "HC"),
    ("LoF_info", "PERCENTILE:0.97,GERP_DIST:12.5"),
    ("LoF_filter", "END_TRUNC"),
    ("LoF_flags", "SINGLE_EXON"),
];

/// The reference's and the candidate's entry of an insertion whose alternate allele settles
/// a stop codon before its unknown bases: the allele is the alternate half of the Codons
/// cell, the reference writes the first of `terms` and no residues, the candidate the
/// second and the residues `residues`, and each side its Extra text of `extra`, the
/// reference's first.
fn settled_stop_entries(
    allele: &str,
    terms: (&str, &str),
    residues: &str,
    extra: (&str, &str),
) -> (Entry, Entry) {
    let entry = |terms: &str, residues: &str, extra: &str| {
        default_entry(&format!(
            "v\t21:10567700-10567701\t{allele}\tENSG1\tENST00000622113\tTranscript\t{terms}\t775-776\t423-424\t141-142\t{residues}\t-/{allele}\t-\t{extra}"
        ))
    };
    (
        entry(terms.0, "-", extra.0),
        entry(terms.1, residues, extra.1),
    )
}

/// An insertion of twelve bases, whole codons, and one of ten: a sense codon and a stop
/// codon before four unknown bases.
const SETTLED_INFRAME: &str = "GAATAGNNNNCA";
const SETTLED_FRAMESHIFT: &str = "GAATAGNNNN";

/// The class of a plugin field beside a stop codon pair reads the fields
/// `term_pair_following_cell` decides from a plugin's terms, those of Downstream, of NMD
/// and LOFTEE, of REVEL and AlphaMissense and of ClinPred and dbNSFP, is the last field
/// class, and names five stop codon predicates, of which each of its five rows is one's
/// alone.
#[test]
fn plugin_value_follows_stop_codon_term_pair_names_its_fields_and_its_five_pairs() {
    let tables = classes::tables();
    let fields: Vec<&str> = DOWNSTREAM_FIELDS
        .iter()
        .chain(TRUNCATING_TERM_FIELDS)
        .chain(MISSENSE_TERM_FIELDS)
        .chain(PROTEIN_CHANGE_TERM_FIELDS)
        .copied()
        .collect();
    assert_eq!(
        fields[7..],
        ["REVEL", "am_", "alphamissense", "ClinPred", "dbNSFP_"]
    );
    assert_eq!(
        tables.field_classes_of(PLUGIN_VALUE_CLASS),
        Some(&fields[..])
    );
    assert_eq!(
        tables.field_classes().last().map(|(klass, _)| *klass),
        Some(PLUGIN_VALUE_CLASS)
    );
    let named: [Pair; 5] = [
        stop_retained_without_stop_codon_test_pair,
        stop_codon_overlap_insertion_range_pair,
        stop_retained_three_base_stop_codon_test_pair,
        stop_codon_settled_before_unknown_bases_pair,
        frameshift_inside_a_held_stop_codon_pair,
    ];
    // After them stand the rules of the five classes of a coding end, which take a pair
    // on the reading of the candidate's cache alone.
    assert_eq!(STOP_CODON_PLUGIN_VALUE_PAIRS.len(), named.len() + 5);
    let (retained, shifted) = ("stop_retained_variant", "frameshift_variant");
    let held = |span: [&str; 4], reference: &str, candidate: &str| {
        let [cdna, cds, allele, location] = span;
        let e = stop_entry(allele, reference, cdna, cds, location);
        let a = stop_entry(allele, candidate, cdna, cds, location);
        let of_the_named: Vec<usize> = (0..named.len())
            .filter(|index| holds(named[*index], &e, &a))
            .collect();
        let held: Vec<usize> = (0..STOP_CODON_PLUGIN_VALUE_PAIRS.len())
            .filter(|index| holds_without_a_cache(STOP_CODON_PLUGIN_VALUE_PAIRS[*index], &e, &a))
            .collect();
        assert_eq!(held, of_the_named, "{span:?}");
        held
    };
    assert_eq!(held(UNTESTED_DELETION, retained, shifted), [0]);
    assert_eq!(held(TESTED_DELETION, retained, shifted), [2]);
    assert_eq!(held(KEPT_STOP_INSERTION, shifted, retained), [1]);
    let (e, a) = settled_stop_entries(
        SETTLED_INFRAME,
        ("coding_sequence_variant", "stop_gained,inframe_insertion"),
        "-/E*XX",
        ("STRAND=1", "STRAND=1"),
    );
    let settled: Vec<usize> = (0..named.len())
        .filter(|index| holds_without_a_cache(STOP_CODON_PLUGIN_VALUE_PAIRS[*index], &e, &a))
        .collect();
    assert_eq!(settled, [3]);
    // An insertion inside the codon of a stop the peptide holds, on which the candidate
    // alone writes `frameshift_variant`.
    let held_stop = |terms: &str| {
        default_entry(&format!(
            "v\t1:134-135\tT\tENSG1\tENST00000000001\tTranscript\t{terms}\t34-35\t28-29\t10\t*/LX\ttga/tTga\t-\tSTRAND=1"
        ))
    };
    let (e, a) = (
        held_stop("coding_sequence_variant"),
        held_stop("frameshift_variant"),
    );
    let inside: Vec<usize> = (0..named.len())
        .filter(|index| holds_without_a_cache(STOP_CODON_PLUGIN_VALUE_PAIRS[*index], &e, &a))
        .collect();
    assert_eq!(inside, [4]);
}

/// Beside the pair of a stop codon settled before unknown bases the class accepts an NMD
/// or LOFTEE field the candidate alone writes where the reference's terms hold no
/// truncating term. A Downstream field, a field on the reference's side, a field both
/// sides write beside `frameshift_variant`, and the same field on a row whose residues are
/// not the codons' translation stay differences.
#[test]
fn plugin_value_follows_stop_codon_term_pair_takes_a_field_beside_a_settled_stop() {
    let unknown = "coding_sequence_variant";
    let gained = "stop_gained,inframe_insertion";
    let accepted = |e: &Entry, a: &Entry, name: &str| {
        class_accepts(
            PLUGIN_VALUE_CLASS,
            e,
            name,
            field(e, name),
            field(a, name),
            Some(a),
        )
    };
    for (name, value) in PLUGIN_VALUES {
        let written = format!("{name}={value};STRAND=1");
        let (e, a) = settled_stop_entries(
            SETTLED_INFRAME,
            (unknown, gained),
            "-/E*XX",
            ("STRAND=1", &written),
        );
        let follows = TRUNCATING_TERM_FIELDS.contains(&name);
        assert_eq!(accepted(&e, &a, name), follows, "{name}");
        let (e, a) = settled_stop_entries(
            SETTLED_INFRAME,
            (unknown, gained),
            "-/E*XX",
            (&written, "STRAND=1"),
        );
        assert!(!accepted(&e, &a, name), "{name} on the reference's side");
        let (e, a) = settled_stop_entries(
            SETTLED_INFRAME,
            (unknown, gained),
            "-/E*XK",
            ("STRAND=1", &written),
        );
        assert!(!accepted(&e, &a, name), "{name} beside other residues");
        let other = format!("{name}=other;STRAND=1");
        let (e, a) = settled_stop_entries(
            SETTLED_FRAMESHIFT,
            ("frameshift_variant", "stop_gained,frameshift_variant"),
            "-/E*XX",
            (&other, &written),
        );
        assert!(
            holds(stop_codon_settled_before_unknown_bases_pair, &e, &a),
            "{name}"
        );
        assert!(!accepted(&e, &a, name), "{name} on both sides");
    }
}

/// The class accepts a Downstream, NMD or LOFTEE field that one side alone writes, the
/// side whose terms allow the plugin to run, on a row whose two entries are a pair of one
/// of three stop codon predicates it names: the untested stop codon, the three-base test
/// and the insertion that keeps the stop codon.
#[test]
fn plugin_value_follows_stop_codon_term_pair_takes_a_plugin_field_beside_each_of_its_pairs() {
    let (retained, shifted) = ("stop_retained_variant", "frameshift_variant");
    let gained = format!("stop_gained,{shifted}");
    let inserted = format!("inframe_insertion,{retained}");
    for (name, value) in PLUGIN_VALUES {
        for span in [UNTESTED_DELETION, TESTED_DELETION] {
            assert!(
                plugin_value_taken(name, span, (retained, None), (shifted, Some(value))),
                "{name} {span:?}"
            );
            assert_eq!(
                plugin_value_taken(name, span, (retained, None), (&gained, Some(value))),
                span == TESTED_DELETION,
                "{name} {span:?}"
            );
        }
        for reference in [shifted, gained.as_str()] {
            for candidate in [retained, inserted.as_str()] {
                assert!(
                    plugin_value_taken(
                        name,
                        KEPT_STOP_INSERTION,
                        (reference, Some(value)),
                        (candidate, None)
                    ),
                    "{name} {reference} {candidate}"
                );
            }
        }
    }
}

/// The class reads a JSON entry's fields under the keys the JSON format writes.
#[test]
fn plugin_value_follows_stop_codon_term_pair_reads_the_keys_of_a_json_entry() {
    let inputs = [STOP_INPUT.to_string()];
    let json_entry = |term: &str, fields: &[(&str, &str)]| {
        let fields: String = fields
            .iter()
            .map(|(name, value)| format!(", {}: {}", quoted(name), quoted(value)))
            .collect();
        let text = format!(
            "{{\"input\": {}, \"id\": \"1993671\", \"allele_string\": \"ACCTG/-\", \"most_severe_consequence\": {term}, \"transcript_consequences\": [{{\"variant_allele\": \"-\", \"transcript_id\": \"ENST00000635895\", \"consequence_terms\": [{term}], \"cdna_start\": 1705, \"cds_start\": 1630{fields}}}]}}\n",
            quoted(STOP_INPUT),
            term = quoted(term)
        );
        json_feature_entries(&text, &inputs).remove(0)
    };
    let e = json_entry("stop_retained_variant", &[]);
    let a = json_entry(
        "frameshift_variant",
        &[
            ("nmd", "NMD_escaping_variant"),
            ("downstreamprotein", "MKV"),
        ],
    );
    for name in ["nmd", "downstreamprotein"] {
        assert!(
            class_accepts(
                PLUGIN_VALUE_CLASS,
                &e,
                name,
                vep_formats_compare::shared::ABSENT,
                field(&a, name),
                Some(&a)
            ),
            "{name}"
        );
    }
}

/// A field on the side whose terms keep the plugin out, a field both sides carry, a field
/// beside two term sets that both allow the plugin to run, a cell that is no plugin's
/// (IMPACT, HGVSp), a row without the candidate's entry and a following cell on a row no
/// stop codon predicate of the class takes stay differences.
#[test]
fn plugin_value_follows_stop_codon_term_pair_keeps_every_other_cell_a_difference() {
    let (retained, shifted) = ("stop_retained_variant", "frameshift_variant");
    let taken = plugin_value_taken;
    for span in [UNTESTED_DELETION, TESTED_DELETION] {
        for (name, reference, candidate) in [
            (
                "DownstreamProtein",
                (retained, Some("MKV")),
                (shifted, None),
            ),
            (
                "DownstreamProtein",
                (retained, Some("MKL")),
                (shifted, Some("MKV")),
            ),
            (
                "NMD",
                (retained, Some("NMD_escaping_variant")),
                (shifted, None),
            ),
            (
                "LoF_flags",
                (retained, Some("SINGLE_EXON")),
                (shifted, None),
            ),
            ("IMPACT", (retained, Some("LOW")), (shifted, Some("HIGH"))),
            ("LoF", (shifted, Some("HC")), (shifted, Some("LC"))),
            ("LoF", (shifted, Some("HC")), (retained, None)),
        ] {
            assert!(
                !taken(name, span, reference, candidate),
                "{name} {span:?} {reference:?} {candidate:?}"
            );
        }
    }
    let kept = format!("{PROTEIN}Ter23=");
    let extended = format!("{PROTEIN}Ter23TyrextTer4");
    assert!(!taken(
        "HGVSp",
        UNTESTED_DELETION,
        (retained, Some(&kept)),
        (shifted, Some(&extended))
    ));
    let spliced = |term: &str| format!("splice_acceptor_variant,{term}");
    for (name, value) in [
        ("NMD", "NMD_escaping_variant"),
        ("LoF", "HC"),
        ("DownstreamProtein", "MKV"),
    ] {
        assert!(
            !taken(
                name,
                UNTESTED_DELETION,
                (&spliced(retained), None),
                (&spliced(shifted), Some(value))
            ),
            "{name}"
        );
    }
    assert!(!taken(
        "LoF_flags",
        KEPT_STOP_INSERTION,
        (shifted, None),
        (retained, Some("SINGLE_EXON"))
    ));
    assert!(!taken(
        "DownstreamProtein",
        KEPT_STOP_INSERTION,
        (shifted, None),
        (retained, Some("MKV"))
    ));
    let substitution = ["1705", "1630", "A", "19:13208046"];
    assert!(!taken(
        "DownstreamProtein",
        substitution,
        (retained, None),
        (shifted, Some("MKV"))
    ));
    assert!(!taken(
        "NMD",
        UNTESTED_DELETION,
        ("inframe_deletion", None),
        (shifted, Some("NMD_escaping_variant"))
    ));
    assert!(!taken(
        "NMD",
        TESTED_DELETION,
        (retained, None),
        ("stop_lost,frameshift_variant", Some("NMD_escaping_variant"))
    ));
    let entry = default_entry(&stop_codon_row("-", retained, "1705-?", "1630-?"));
    assert!(!class_accepts(
        PLUGIN_VALUE_CLASS,
        &entry,
        "NMD",
        vep_formats_compare::shared::ABSENT,
        "NMD_escaping_variant",
        None
    ));
}

/// The class takes an insertion between two neighbouring positions and no wider row. A
/// default-format row has a Location of two neighbouring positions and reference bases that
/// are empty, however its record pads the insertion; else, where its first column spells the
/// record's alleles, a reference allele that is empty once the shared first base is dropped.
/// Under an identifier the records of the input give the bases: an insertion of two bases is
/// taken and two replaced bases are refused. Without them the allele has another length than
/// two, so two bases replaced by one or by three are taken with the insertions.
#[test]
fn stop_codon_overlap_insertion_range_pair_takes_an_insertion_between_two_positions_alone() {
    use vep_formats_compare::rules::term_pairs::read_entry_references;
    let pair: Pair = stop_codon_overlap_insertion_range_pair;
    let two = "19:13208046-13208047";
    let entries = |(name, location, allele): (&str, &str, &str)| {
        let entry = |terms: &str| {
            default_entry(
                &stop_codon_row_at(allele, terms, "7054-7055", "6977-6978", location)
                    .replacen("1993671", name, 1),
            )
        };
        (
            vec![entry("coding_sequence_variant")],
            vec![entry("stop_retained_variant")],
        )
    };
    let text = |row: (&str, &str, &str)| {
        let (expected, actual) = entries(row);
        holds(pair, &expected[0], &actual[0])
    };
    let directory = compare_common::scratch("term-pairs-insertion-range-records");
    let beside = |(lines, allele): (&str, &str)| {
        let input = compare_common::planted(
            &directory,
            "input.vcf",
            &format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{lines}"
            ),
        );
        let (mut expected, mut actual) = entries(("1993671", two, allele));
        read_entry_references(&mut expected, &mut actual, &input).expect("the input reads");
        holds(pair, &expected[0], &actual[0])
    };
    let wider = [
        ("19_13208046_AC/GTT", two, "GTT"),
        ("19_13208046_AC/AGT", two, "AGT"),
        ("19_13208046_AC/A/T", two, "A"),
        ("19_13208046_A/AC/T", "19:13208046", "AC"),
        ("19_13208046_-/TGA", "19:13208046", "TGA"),
        ("1993671", two, "CG"),
        ("19_13208046_AC/A/T", two, "GTT"),
    ];
    let taken: Vec<_> = wider.into_iter().filter(|row| text(*row)).collect();
    assert!(taken.is_empty(), "{taken:?}");
    let inserted = [
        ("19_13208046_-/TGA", two, "TGA"),
        ("19_13208046_A/ATGA", two, "TGA"),
        ("19_13208046_-/CG/T", two, "CG"),
        ("1993671", two, "TGA"),
        ("1993671", two, "G"),
        ("1993671", two, "CGT"),
        ("19_13208045_TA/TATGA", two, "TGA"),
        ("19_13208046_AC/ATGAC", two, "TGA"),
        ("19_13208046_A/ATGA/C", two, "TGA"),
    ];
    let refused: Vec<_> = inserted.into_iter().filter(|row| !text(*row)).collect();
    assert!(refused.is_empty(), "{refused:?}");

    // The last record inserts before its own first base and gives the row no bases.
    let beside_an_insertion = [
        ("19\t13208046\t1993671\tA\tACG\t.\t.\t.\n", "CG"),
        ("19\t13208045\t1993671\tTA\tTACG\t.\t.\t.\n", "CG"),
        ("19\t13208046\t1993671\tA\tATGA\t.\t.\t.\n", "TGA"),
        ("19\t13208047\t1993671\tA\tCGTA\t.\t.\t.\n", "CGT"),
    ];
    let refused: Vec<_> = beside_an_insertion
        .into_iter()
        .filter(|row| !beside(*row))
        .collect();
    assert!(refused.is_empty(), "{refused:?}");
    let replaced = [
        ("19\t13208046\t1993671\tAC\tCGT\t.\t.\t.\n", "CGT"),
        ("19\t13208046\t1993671\tAC\tG\t.\t.\t.\n", "G"),
        ("19\t13208046\t1993671\tAC\tCG\t.\t.\t.\n", "CG"),
        (
            "19\t13208046\t1993672\tAC\tCGT\t.\t.\t.\n19\t13208047\t1993671\tA\tCGTA\t.\t.\t.\n",
            "CGT",
        ),
    ];
    let taken: Vec<_> = replaced.into_iter().filter(|row| beside(*row)).collect();
    assert!(taken.is_empty(), "{taken:?}");
}

/// A VCF or JSON entry carries no Location and is read by the record's own alleles.
#[test]
fn stop_codon_overlap_insertion_range_pair_reads_the_alleles_of_a_vcf_or_json_record() {
    let pair: Pair = stop_codon_overlap_insertion_range_pair;
    let vcf = |(reference, alternate, allele): (&str, &str, &str)| {
        let input = format!("19\t13208046\t1993671\t{reference}\t{alternate}\t.\t.\t.");
        let inputs = [input.clone()];
        let entry = |terms: &str| {
            let text = format!(
                "{}{}CSQ={allele}|{terms}|ENST00000635895|Transcript|7054-7055|6977-6978\n",
                csq_header("Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position"),
                open_info(&input)
            );
            let entry = vcf_feature_entries(&text, &inputs)
                .into_iter()
                .find(|entry| entry.feature == "ENST00000635895")
                .expect("the transcript's entry");
            assert!(entry.location.is_none());
            entry
        };
        holds(
            pair,
            &entry("coding_sequence_variant"),
            &entry("stop_retained_variant"),
        )
    };
    for row in [
        ("A", "ATGA", "TGA"),
        ("A", "ACG", "CG"),
        ("A", "ATGA,AG", "TGA"),
    ] {
        assert!(vcf(row), "{row:?}");
    }
    let wider = [
        ("AC", "GTT", "GTT"),
        ("AC", "AGT", "GT"),
        ("AC", "A,T", "A"),
        ("A", "AC,T", "AC"),
        ("AC", "A,T", "GTT"),
    ];
    let taken: Vec<_> = wider.into_iter().filter(|row| vcf(*row)).collect();
    assert!(taken.is_empty(), "{taken:?}");

    let input = "19\t13208046\t1993671\tA\tATGA\t.\t.\t.";
    let inputs = [input.to_string()];
    let obj = |(allele_string, allele): (&str, &str)| {
        let entry = |terms: &str| {
            let text = format!(
                "{{\"input\": {}, \"id\": \"1993671\", \"allele_string\": {}, \"most_severe_consequence\": {terms}, \"transcript_consequences\": [{{\"variant_allele\": {}, \"transcript_id\": \"ENST00000635895\", \"consequence_terms\": [{terms}], \"cdna_start\": 7054, \"cdna_end\": 7055, \"cds_start\": 6977, \"cds_end\": 6978}}]}}\n",
                quoted(input),
                quoted(allele_string),
                quoted(allele),
                terms = quoted(terms)
            );
            json_feature_entries(&text, &inputs)
                .into_iter()
                .find(|entry| entry.feature == "ENST00000635895")
                .expect("the transcript's entry")
        };
        holds(
            pair,
            &entry("coding_sequence_variant"),
            &entry("stop_retained_variant"),
        )
    };
    assert!(obj(("-/TGA", "TGA")) && obj(("-/CG", "CG")));
    let wider = [("AC/GTT", "GTT"), ("C/GT", "GT"), ("AC/A/T", "A")];
    let taken: Vec<_> = wider.into_iter().filter(|row| obj(*row)).collect();
    assert!(taken.is_empty(), "{taken:?}");
}

/// Without its record, a row under an identifier, the shape is the reference's splice
/// donor, splice acceptor, donor fifth base, donor region, polypyrimidine tract or intron
/// term against none of the six on the candidate's side, with or without
/// splice_region_variant where a term that withholds it went, on a Transcript row whose
/// Location holds three or more positions, whose allele holds more nucleotide letters than
/// the span has bases and whose two cDNA positions, the same on both sides as the CDS
/// positions are, are numbers as far apart as the span's two ends.
#[test]
fn splice_term_past_the_record_span_pair_reads_a_longer_allele_over_a_span_in_the_transcript_sequence(
) {
    let pair: Pair = splice_term_past_the_record_span_pair;
    let row = |terms: &str, columns: &SpanColumns| {
        let entry = default_entry(&splice_span_row(terms, columns));
        assert_eq!(entry.reference, None);
        entry
    };
    let plain = SpanColumns::default();
    assert!(holds(
        pair,
        &row(DONOR_PAST_THE_SPAN, &plain),
        &row(REGION_IN_ITS_PLACE, &plain)
    ));
    assert!(holds(
        pair,
        &row(DONOR_PAST_THE_SPAN, &plain),
        &row("frameshift_variant", &plain)
    ));
    for gone in [
        "splice_acceptor_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_donor_variant,splice_donor_region_variant",
        "splice_donor_variant,splice_donor_5th_base_variant",
        "splice_acceptor_variant,splice_polypyrimidine_tract_variant",
        "splice_donor_variant,intron_variant",
        "splice_donor_variant,splice_donor_region_variant,intron_variant",
    ] {
        let reference = row(&format!("{gone},frameshift_variant"), &plain);
        assert!(
            holds(pair, &reference, &row(REGION_IN_ITS_PLACE, &plain)),
            "{gone}"
        );
        assert!(
            holds(pair, &reference, &row("frameshift_variant", &plain)),
            "{gone}"
        );
    }
    for gone in [
        "splice_polypyrimidine_tract_variant",
        "intron_variant",
        "splice_polypyrimidine_tract_variant,intron_variant",
    ] {
        let reference = row(&format!("{gone},frameshift_variant"), &plain);
        assert!(
            holds(pair, &reference, &row("frameshift_variant", &plain)),
            "{gone}"
        );
        assert!(
            !holds(pair, &reference, &row(REGION_IN_ITS_PLACE, &plain)),
            "{gone}"
        );
    }
    let taken = [
        SpanColumns {
            cdna: "998-1000/1200",
            cds: "898-900/1000",
            ..plain
        },
        SpanColumns {
            cdna: "38-40",
            cds: "-",
            ..plain
        },
        SpanColumns {
            allele: "ACNTA",
            ..plain
        },
        SpanColumns {
            allele: "ACGT",
            ..plain
        },
        SpanColumns {
            cdna: "996-1000",
            cds: "896-900",
            allele: "ACGTACGT",
            location: "7:100996-101000",
            ..plain
        },
    ];
    for columns in &taken {
        assert!(
            holds(
                pair,
                &row(DONOR_PAST_THE_SPAN, columns),
                &row(REGION_IN_ITS_PLACE, columns)
            ),
            "{columns:?}"
        );
    }
    let refused = [
        (
            "a span across an exon's end",
            SpanColumns {
                cdna: "999-?",
                cds: "899-?",
                ..plain
            },
        ),
        (
            "a span that begins before an exon",
            SpanColumns {
                cdna: "?-1000",
                cds: "?-900",
                ..plain
            },
        ),
        (
            "a span in an intron",
            SpanColumns {
                cdna: "-",
                cds: "-",
                ..plain
            },
        ),
        (
            "a span with an intron inside it",
            SpanColumns {
                cdna: "998-999",
                cds: "898-899",
                ..plain
            },
        ),
        (
            "an allele as long as the span",
            SpanColumns {
                allele: "ACG",
                ..plain
            },
        ),
        (
            "an allele shorter than the span",
            SpanColumns {
                allele: "AC",
                ..plain
            },
        ),
        (
            "a deletion",
            SpanColumns {
                allele: "-",
                ..plain
            },
        ),
        (
            "an allele of a letter that names no base",
            SpanColumns {
                allele: "ACRTA",
                ..plain
            },
        ),
        (
            "two neighbouring positions",
            SpanColumns {
                cdna: "999-1000",
                cds: "899-900",
                location: "7:100999-101000",
                ..plain
            },
        ),
        (
            "one position",
            SpanColumns {
                cdna: "1000",
                cds: "900",
                location: "7:101000",
                ..plain
            },
        ),
        (
            "a regulatory row",
            SpanColumns {
                feature_type: "RegulatoryFeature",
                ..plain
            },
        ),
    ];
    for (shape, columns) in &refused {
        assert!(
            !holds(
                pair,
                &row(DONOR_PAST_THE_SPAN, columns),
                &row(REGION_IN_ITS_PLACE, columns)
            ),
            "{shape}"
        );
    }
    for kept in ["splice_donor_5th_base_variant", "intron_variant"] {
        assert!(
            !holds(
                pair,
                &row(
                    &format!("splice_donor_variant,{kept},frameshift_variant"),
                    &plain
                ),
                &row(&format!("{kept},frameshift_variant"), &plain)
            ),
            "{kept}"
        );
    }
    let with_a_further_term = format!("{REGION_IN_ITS_PLACE},stop_gained");
    for terms in [
        with_a_further_term.as_str(),
        "splice_region_variant",
        "frameshift_variant,intron_variant",
        "splice_donor_region_variant,frameshift_variant",
        DONOR_PAST_THE_SPAN,
    ] {
        assert!(
            !holds(pair, &row(DONOR_PAST_THE_SPAN, &plain), &row(terms, &plain)),
            "{terms}"
        );
    }
    assert!(!holds(
        pair,
        &row(REGION_IN_ITS_PLACE, &plain),
        &row(DONOR_PAST_THE_SPAN, &plain)
    ));
    for columns in [
        SpanColumns {
            cdna: "997-999",
            ..plain
        },
        SpanColumns {
            cds: "899-901",
            ..plain
        },
        SpanColumns {
            cdna: "999-?",
            cds: "899-?",
            ..plain
        },
    ] {
        assert!(
            !holds(
                pair,
                &row(DONOR_PAST_THE_SPAN, &plain),
                &row(REGION_IN_ITS_PLACE, &columns)
            ),
            "{columns:?}"
        );
        assert!(
            !holds(
                pair,
                &row(DONOR_PAST_THE_SPAN, &columns),
                &row(REGION_IN_ITS_PLACE, &plain)
            ),
            "{columns:?}"
        );
    }
}

/// The records of `splice_span_records` and the two replaced bases, each with its name.
fn records_with_two_replaced() -> Vec<(&'static str, SpliceSpanRecord)> {
    let mut records = splice_span_records();
    records.push(("two replaced bases", two_replaced()));
    records
}

/// What a record's entries of `fmt` carry and whether the pair takes them, for the test of
/// each format.
fn the_record_is_read_in(fmt: &str) {
    let pair: Pair = splice_term_past_the_record_span_pair;
    for (shape, record) in records_with_two_replaced() {
        let (e, a) = record.in_format(fmt);
        let bases = &record.line[2][1..];
        assert_eq!(e.reference.as_deref(), Some(bases), "{shape} {fmt}");
        assert_eq!(a.reference.as_deref(), Some(bases), "{shape} {fmt}");
        for entry in [&e, &a] {
            let location = entry.location.as_ref().or(entry.allele_location.as_ref());
            assert_eq!(location, Some(&record.location), "{shape} {fmt}");
        }
        assert_eq!(e.location.is_none(), fmt != "default", "{shape} {fmt}");
        assert!(holds(pair, &e, &a), "{shape} {fmt}");
        assert!(!holds(pair, &a, &e), "{shape} {fmt}");
    }
    let (e, a) = inserted_between_two().in_format(fmt);
    assert_eq!(e.reference.as_deref(), Some(""), "{fmt}");
    assert_eq!(a.reference.as_deref(), Some(""), "{fmt}");
    let location = e.location.as_ref().or(e.allele_location.as_ref());
    assert_eq!(location, Some(&two_replaced().location), "{fmt}");
    assert!(!holds(pair, &e, &a), "{fmt}");
}

/// An entry carries the reference bases its allele stands in for where its record's
/// alleles are at hand: a default row from a first column that spells them. With them the
/// shape is two or more reference bases under a longer allele beside one intron that gives
/// the reference's terms from the record's differing regions and the candidate's from the
/// span, or from the insertion an allele that opens with its whole reference allele writes
/// after it. A pure insertion is refused where two replaced bases at its Location are taken,
/// and so are a candidate's row the span does not give, one replaced base, an allele not
/// longer than the reference bases and a row under an identifier.
#[test]
fn splice_term_past_the_record_span_pair_reads_the_record_where_an_entry_carries_it() {
    let pair: Pair = splice_term_past_the_record_span_pair;
    the_record_is_read_in("default");
    for (shape, record) in records_with_two_replaced() {
        let (e, a) = record.entries(Some("rec2"));
        assert_eq!(e.reference, None, "{shape}");
        assert!(!holds(pair, &e, &a), "{shape}");
    }
    let exon = splice_span_record("two positions inside an exon");
    let mut one_replaced_base = exon.clone().with_line(["7", "100999", "AC", "ACGTA"]);
    one_replaced_base.location = "7:101000".to_string();
    one_replaced_base.allele = "CGTA".to_string();
    one_replaced_base.cdna = "1000".to_string();
    one_replaced_base.cds = "900".to_string();
    let mut as_long_as_the_bases = exon.clone().with_line(["7", "100999", "CG", "TA"]);
    as_long_as_the_bases.allele = "TA".to_string();
    let refused = [
        (
            "nothing in the place of the donor term, where the record's regions put splice_region_variant",
            exon.clone().with_terms(DONOR_PAST_THE_SPAN, "frameshift_variant"),
        ),
        (
            "the intron term on the candidate's side of a span inside an exon",
            exon.clone().with_terms(
                &format!("{DONOR_PAST_THE_SPAN},intron_variant"),
                &format!("{REGION_IN_ITS_PLACE},intron_variant"),
            ),
        ),
        ("one replaced base", one_replaced_base),
        ("an allele as long as the reference bases", as_long_as_the_bases),
    ];
    for (shape, record) in refused {
        let (e, a) = record.entries(None);
        assert!(
            e.reference
                .as_deref()
                .is_some_and(|bases| !bases.is_empty()),
            "{shape}"
        );
        assert!(!holds(pair, &e, &a), "{shape}");
    }
    let (e, a) = exon.entries(None);
    assert!(holds(pair, &e, &a));
    let mut other_positions = exon.clone();
    other_positions.cdna = "998-999".to_string();
    assert!(!holds(
        pair,
        &e,
        &default_entry(&other_positions.rows(None).1)
    ));
}

/// A VCF entry carries the bases from its line's POS, REF and ALT and a JSON entry from its
/// object's span and allele string, both with the Location the default format writes. An
/// allele the line's ALT does not write at the line's span states no record: the entry
/// carries neither the bases nor a Location.
#[test]
fn splice_term_past_the_record_span_pair_reads_the_record_of_a_vcf_or_json_entry() {
    let pair: Pair = splice_term_past_the_record_span_pair;
    for fmt in ["vcf", "json"] {
        the_record_is_read_in(fmt);
        let mut unwritten = splice_span_record("two positions inside an exon");
        unwritten.allele = "GTAC".to_string();
        let (e, a) = unwritten.in_format(fmt);
        assert_eq!(e.reference, None, "{fmt}");
        assert_eq!(e.allele_location, None, "{fmt}");
        assert!(!holds(pair, &e, &a), "{fmt}");
    }
}

/// The index of the row class `klass` among the built-in row classes.
fn term_pair_class(klass: &str) -> usize {
    TERM_PAIR_CLASSES
        .iter()
        .position(|(name, _)| *name == klass)
        .unwrap_or_else(|| panic!("{klass} is no row class"))
}

/// The row classes are the class table's, in its order; each class whose cells follow its
/// terms is a field class of the cells `term_pair_following_cell` reads, the two classes
/// whose rows are derived on a completed coding sequence and the class of a stop codon
/// settled before unknown bases each of four more and the classes of a start codon removed
/// without a 5' UTR and of a start codon kept without a 5' UTR of the IMPACT and the protein
/// notation alone, the class read at a stop codon's unshifted span stands before every one
/// of them, the class of the splice term past the record's span stands after the stop codon
/// classes, and the class of a stop codon settled before unknown bases, the class of a start
/// codon removed without a 5' UTR, the class of a frameshift inside a held stop codon and
/// the class of a start codon kept without a 5' UTR are the last four.
#[test]
fn the_row_classes_stand_in_the_order_of_the_class_table() {
    let tables = classes::tables();
    let names: Vec<&str> = TERM_PAIR_CLASSES.iter().map(|(name, _)| *name).collect();
    assert_eq!(names, tables.term_pair_classes());
    let following = tables.term_pair_following_classes();
    assert_eq!(
        following,
        [
            classes::STRUCTURAL_DELETION_STOP_LOST_WITHOUT_STOP_CODON,
            classes::STOP_CODON_COMPLETED_BY_POLYADENYLATION,
            classes::STOP_CODON_OVERLAP_WITHOUT_STOP_CODON,
            classes::STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON,
            classes::START_CODON_OVERLAP_WITHOUT_START_CODON,
            classes::STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
            classes::STOP_CODON_OVERLAP_GENOMIC_SPAN,
            classes::STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
            classes::STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
            classes::STOP_CODON_OVERLAP_INSERTION_RANGE,
            classes::REFSEQ_EDITED_MODEL_EQUAL_SPAN_INSERTION_TERM,
            classes::SPLICE_TERM_PAST_THE_RECORD_SPAN,
            classes::STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES,
            classes::START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR,
            classes::FRAMESHIFT_INSIDE_A_HELD_STOP_CODON,
            classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR,
        ]
    );
    let unshifted = term_pair_class(classes::STOP_CODON_UNSHIFTED_SPAN);
    // The two classes whose rows are derived on a completed coding sequence read the
    // candidate's residues and codons as well.
    let derived_cells = [
        classes::STOP_CODON_COMPLETED_BY_POLYADENYLATION,
        classes::STOP_CODON_OVERLAP_WITHOUT_STOP_CODON,
    ];
    let with_the_cells: Vec<&str> = TERM_PAIR_FOLLOWING_FIELDS
        .iter()
        .copied()
        .chain(["Amino_acids", "amino_acids", "Codons", "codons"])
        .collect();
    // The class of a stop codon settled before unknown bases reads, after the cells that
    // follow its terms, the residues and the offset of the notation its rows gain.
    let settled_stop_fields: Vec<&str> = TERM_PAIR_FOLLOWING_FIELDS
        .iter()
        .copied()
        .chain(["Amino_acids", "amino_acids", "HGVS_OFFSET", "hgvs_offset"])
        .collect();
    for klass in following {
        assert!(unshifted < term_pair_class(klass), "{klass}");
        let fields = if *klass == classes::STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES {
            &settled_stop_fields[..]
        } else if *klass == classes::START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR {
            // The terms of its rows let no plugin run on one side alone.
            &["IMPACT", "impact", "HGVSp", "hgvsp"][..]
        } else if derived_cells.contains(klass) {
            &with_the_cells[..]
        } else if *klass == classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR {
            // Its two sets of terms differ by the start terms alone, which let no plugin
            // run on one side alone.
            &["IMPACT", "impact", "HGVSp", "hgvsp"][..]
        } else {
            TERM_PAIR_FOLLOWING_FIELDS
        };
        assert_eq!(tables.field_classes_of(klass), Some(fields), "{klass}");
    }
    // The five classes of a coding end stand in the order they are asked, before the
    // classes of a stop codon that stands.
    let coding_end: Vec<&str> = CodingEndClass::IN_ORDER
        .iter()
        .map(|class| class.name())
        .collect();
    assert_eq!(names[unshifted + 1..unshifted + 6], coding_end[..]);
    assert_eq!(following[..5], coding_end[..]);
    assert_eq!(
        names[names.len() - 5..],
        [
            SPLICE_SPAN_CLASS,
            classes::STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES,
            classes::START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR,
            classes::FRAMESHIFT_INSIDE_A_HELD_STOP_CODON,
            classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR
        ]
    );
    assert_eq!(classes::SPLICE_TERM_PAST_THE_RECORD_SPAN, SPLICE_SPAN_CLASS);
    assert!(
        term_pair_class(classes::STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST)
            < term_pair_class(SPLICE_SPAN_CLASS)
    );
    assert_eq!(
        tables.record_reference_classes(),
        [
            SPLICE_SPAN_CLASS,
            SPLICE_SPAN_PLUGIN_VALUE_CLASS,
            classes::STOP_CODON_OVERLAP_INSERTION_RANGE,
            PLUGIN_VALUE_CLASS,
            classes::STOP_CODON_COMPLETED_BY_POLYADENYLATION,
            classes::STOP_CODON_OVERLAP_WITHOUT_STOP_CODON,
            classes::STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON,
            classes::START_CODON_OVERLAP_WITHOUT_START_CODON,
            classes::REFSEQ_EDITED_MODEL_EQUAL_SPAN_INSERTION_TERM,
            classes::HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE,
        ]
    );
    assert_eq!(
        classes::PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN,
        SPLICE_SPAN_PLUGIN_VALUE_CLASS
    );
    let plugin_fields: Vec<&str> = DOWNSTREAM_FIELDS
        .iter()
        .chain(TRUNCATING_TERM_FIELDS)
        .copied()
        .collect();
    assert_eq!(
        tables.field_classes_of(SPLICE_SPAN_PLUGIN_VALUE_CLASS),
        Some(&plugin_fields[..])
    );
    let last_two: Vec<&str> = tables
        .field_classes()
        .iter()
        .rev()
        .take(2)
        .rev()
        .map(|(klass, _)| *klass)
        .collect();
    assert_eq!(
        last_two,
        [SPLICE_SPAN_PLUGIN_VALUE_CLASS, PLUGIN_VALUE_CLASS]
    );
}

/// Each stop codon class and the class of the splice term past the record's span decide a
/// row by the predicate of their own name.
#[test]
fn each_row_class_asks_the_predicate_of_its_name() {
    let named: [(&str, Pair); 6] = [
        (
            classes::STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
            stop_retained_without_stop_codon_test_pair,
        ),
        (
            classes::STOP_CODON_OVERLAP_GENOMIC_SPAN,
            stop_codon_overlap_genomic_span_pair,
        ),
        (
            classes::STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
            stop_retained_three_base_stop_codon_test_pair,
        ),
        (
            classes::STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
            stop_codon_removed_without_stop_codon_test_pair,
        ),
        (
            classes::STOP_CODON_OVERLAP_INSERTION_RANGE,
            stop_codon_overlap_insertion_range_pair,
        ),
        (
            classes::SPLICE_TERM_PAST_THE_RECORD_SPAN,
            splice_term_past_the_record_span_pair,
        ),
    ];
    // One row each class takes, in the order above, and every class is asked about every row.
    let plain = SpanColumns::default();
    let rows = [
        (
            stop_entry("-", RETAINED, "?-7057", "-", STOP_CODON_LOCATION),
            stop_entry("-", UNKNOWN, "?-7057", "-", STOP_CODON_LOCATION),
        ),
        (
            stop_entry(
                "-",
                "splice_region_variant,coding_sequence_variant,NMD_transcript_variant",
                "1705",
                "1630",
                "3:36993663",
            ),
            stop_entry(
                "-",
                "splice_region_variant,stop_retained_variant,NMD_transcript_variant",
                "1705",
                "1630",
                "3:36993663",
            ),
        ),
        (
            stop_entry(
                "-",
                RETAINED,
                "2558-2564",
                "2249-?",
                "2:190970701-190970707",
            ),
            stop_entry("-", UNKNOWN, "2558-2564", "2249-?", "2:190970701-190970707"),
        ),
        (
            stop_entry(
                "-",
                "stop_retained_variant",
                "?-7057",
                "-",
                STOP_CODON_LOCATION,
            ),
            stop_entry("-", "stop_lost", "?-7057", "-", STOP_CODON_LOCATION),
        ),
        (
            stop_entry(
                "NNNN",
                "stop_gained,frameshift_variant",
                "7054-7055",
                "6977-6978",
                "19:13208046-13208047",
            ),
            stop_entry(
                "NNNN",
                "inframe_insertion,stop_retained_variant",
                "7054-7055",
                "6977-6978",
                "19:13208046-13208047",
            ),
        ),
        (
            default_entry(&splice_span_row(DONOR_PAST_THE_SPAN, &plain)),
            default_entry(&splice_span_row(REGION_IN_ITS_PLACE, &plain)),
        ),
    ];
    for (taken, (klass, pair)) in named.into_iter().enumerate() {
        let (_, of_the_table) = TERM_PAIR_CLASSES[term_pair_class(klass)];
        for (row, (e, a)) in rows.iter().enumerate() {
            let answer = holds(pair, e, a);
            assert_eq!(
                holds_without_a_cache(of_the_table, e, a),
                answer,
                "{klass} on row {row}"
            );
            assert_eq!(answer, row == taken, "{klass} on row {row}");
        }
    }
}

/// The entries of both sides whose key's two Consequence texts differ take the reference
/// bases their allele stands in for from the records of the input: a record under a `chr`
/// prefix is the row's, two records that read a row differently give it nothing, bases an
/// entry already carries stay, and of two candidate entries under one key the last is the
/// key's. An input that cannot be read is the comparator's error where an entry asks, and is
/// not opened where none does.
#[test]
fn the_reference_bases_of_differing_rows_are_read_from_the_records_of_the_input() {
    use vep_formats_compare::rules::term_pairs::read_entry_references;
    let directory = compare_common::scratch("term-pairs-entry-references");
    let input = compare_common::planted(
        &directory,
        "input.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n7\t100998\trec2\tACG\tACGTAC\t.\t.\t.\n7\t400998\trec5\tACG\tAGTA\t.\t.\t.\n7\t400999\trec6\tC\tCGTA\t.\t.\t.\n",
    );
    let row = |key: &str, location: &str, allele: &str, terms: &str, bases: Option<&str>| Entry {
        key: key.to_string(),
        location: Some(location.to_string()),
        allele: allele.to_string(),
        feature: "ENST1".to_string(),
        consequence: terms.to_string(),
        reference: bases.map(str::to_string),
        ..Entry::default()
    };
    let exon = "7:100999-101000";
    let rows = || {
        let expected = vec![
            row("a", exon, "CGTAC", "splice_donor_variant", None),
            row("b", exon, "CGTAC", "intron_variant", None),
            row("c", "7:400999-401000", "GTA", "splice_donor_variant", None),
            row(
                "d",
                "7:10099\u{b2}-101000",
                "CGTAC",
                "splice_donor_variant",
                None,
            ),
            row(
                "e",
                "chr7:100999-101000",
                "CGTAC",
                "splice_donor_variant",
                Some("ACG"),
            ),
            row("f", exon, "CGTAC", "splice_donor_variant", None),
            row("g", exon, "CGTAC", "a,b", None),
        ];
        let actual = vec![
            row("g", "chr7:100999-101000", "CGTAC", "b,a", None),
            row("a", exon, "CGTAC", "splice_region_variant", None),
            row("b", exon, "CGTAC", "intron_variant", None),
            row("c", "7:400999-401000", "GTA", "splice_region_variant", None),
            row("d", exon, "CGTAC", "splice_region_variant", None),
            row("e", exon, "CGTAC", "splice_region_variant", None),
            row("a", exon, "CGTAC", "splice_donor_variant", None),
        ];
        (expected, actual)
    };
    let bases = |entries: &[Entry]| -> Vec<Option<String>> {
        entries
            .iter()
            .map(|entry| entry.reference.clone())
            .collect()
    };
    let read = |texts: [Option<&str>; 7]| texts.map(|text| text.map(str::to_string)).to_vec();
    let (mut expected, mut actual) = rows();
    read_entry_references(&mut expected, &mut actual, &input).expect("the input reads");
    assert_eq!(
        bases(&expected),
        read([None, None, None, None, Some("ACG"), None, Some("CG")])
    );
    assert_eq!(
        bases(&actual),
        read([Some("CG"), None, None, None, Some("CG"), Some("CG"), None])
    );
    let missing = directory.join("none.vcf");
    let (mut expected, mut actual) = rows();
    let refused = read_entry_references(&mut expected, &mut actual, &missing)
        .expect_err("the input does not open");
    assert_eq!(
        refused.to_string(),
        format!(
            "cannot read {}: No such file or directory",
            missing.display()
        )
    );
    let mut agreeing = vec![row("b", exon, "CGTAC", "intron_variant", None)];
    let mut alike = agreeing.clone();
    read_entry_references(&mut agreeing, &mut alike, &missing).expect("no entry asks");
    assert_eq!(bases(&agreeing), [None]);
}

/// The bases the records of a made-up input give the reference's row of one key: a row at
/// `location` with the allele `C`, under the first column `name`, whose two Consequence
/// texts differ. The input holds, at four places of chromosome 7, two records that read
/// such a row differently (`ins1` inserts the allele, `sub1` replaces two bases by it),
/// two records that read it alike (`same1`, `same2`), one record (`only`), and two
/// records under one identifier that read it differently (`twice`).
fn bases_of_the_row_named(test: &str, name: &str, location: &str) -> Option<String> {
    use vep_formats_compare::rules::term_pairs::read_entry_references;
    let directory = compare_common::scratch(&format!("term-pairs-named-record-{test}"));
    let records = [
        "7\t400998\tins1\tA\tAC",
        "7\t400998\tsub1\tAT\tC",
        "7\t500998\tsame1\tA\tAC",
        "7\t500998\tsame2\tA\tAC",
        "7\t600998\tonly\tA\tAC",
        "7\t700998\ttwice\tA\tAC",
        "7\t700998\ttwice\tAT\tC",
    ];
    let lines: Vec<String> = records
        .iter()
        .map(|record| format!("{record}\t.\t.\t.\n"))
        .collect();
    let input = compare_common::planted(
        &directory,
        "input.vcf",
        &format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}",
            lines.concat()
        ),
    );
    let row = |terms: &str| Entry {
        key: "k".to_string(),
        location: Some(location.to_string()),
        allele: "C".to_string(),
        feature: "ENST1".to_string(),
        consequence: terms.to_string(),
        fields: [("Uploaded_variation".to_string(), name.to_string())]
            .into_iter()
            .collect(),
        ..Entry::default()
    };
    let mut expected = vec![row("frameshift_variant")];
    let mut actual = vec![row("inframe_insertion")];
    read_entry_references(&mut expected, &mut actual, &input).expect("the input reads");
    assert_eq!(expected[0].reference, actual[0].reference, "{name}");
    expected.remove(0).reference
}

const TWO_READINGS: &str = "7:400998-400999";

/// Where two records read a row differently, the record whose identifier is the row's own
/// first column decides.
#[test]
fn a_row_two_records_read_differently_takes_the_bases_of_the_record_it_names() {
    let bases = |name: &str| bases_of_the_row_named("names-one", name, TWO_READINGS);
    assert_eq!(bases("ins1"), Some(String::new()));
    assert_eq!(bases("sub1"), Some("AT".to_string()));
}

/// A row that names neither of two records that read it differently takes no bases, and
/// neither does a row whose identifier two records carry that read it differently.
#[test]
fn a_row_that_names_none_of_the_records_that_read_it_differently_takes_no_bases() {
    for name in ["other", "7_400998_A/AC", "INS1", "", "."] {
        assert_eq!(
            bases_of_the_row_named("names-none", name, TWO_READINGS),
            None,
            "{name:?}"
        );
    }
    assert_eq!(
        bases_of_the_row_named("names-two", "twice", "7:700998-700999"),
        None
    );
}

/// A row the records read alike takes their bases, whatever its first column names.
#[test]
fn a_row_the_records_read_alike_takes_their_bases_whatever_it_names() {
    for name in ["same1", "same2", "other", "ins1"] {
        assert_eq!(
            bases_of_the_row_named("alike", name, "7:500998-500999"),
            Some(String::new()),
            "{name}"
        );
    }
}

/// A row one record reads takes its bases, whatever its first column names.
#[test]
fn a_row_one_record_reads_takes_its_bases_whatever_it_names() {
    for name in ["only", "other", "sub1"] {
        assert_eq!(
            bases_of_the_row_named("one-record", name, "7:600998-600999"),
            Some(String::new()),
            "{name}"
        );
    }
}
