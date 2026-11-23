// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `refseq_edited_model_hgvsc_offset` on a deletion across the bases an edited model lacks:
//! the reference's notation closes later than the candidate's by those bases.

mod compare_common;
mod shifted_edited_common;

use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::consts::HGVSC_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson::Value;

const HGVSC_OFFSET: &str = "refseq_edited_model_hgvsc_offset";

/// Ensembl VEP 116.3's rows and vep-rs's rows of two records under `--refseq --hgvs`,
/// default format: `rs1179352929` of gnomAD v4.1 and ClinVar variation `1338712`, each on
/// `NM_030582.4` and `NM_130444.3`, whose alignment deleted the nine genomic bases
/// `21:45504526-45504534`.
const REFERENCE_ROWS: &str =
    include_str!("data/edited_model_deleted_bases/ensembl-vep.default.txt");
const CANDIDATE_ROWS: &str = include_str!("data/edited_model_deleted_bases/vep-rs.default.txt");

const MODEL: &str = "NM_990004.1";

/// The entry's own HGVSc value.
fn hgvsc_of(entry: &Entry) -> &str {
    HGVSC_FIELDS
        .iter()
        .find_map(|field| entry.fields.get(*field))
        .map_or("-", String::as_str)
}

/// Whether the class accepts the candidate entry's HGVSc against the reference entry's.
fn accepted(reference: &Entry, candidate: &Entry) -> bool {
    let name = HGVSC_FIELDS
        .iter()
        .find(|field| reference.fields.contains_key(**field))
        .expect("the entry carries an HGVSc");
    class_accepts(
        HGVSC_OFFSET,
        reference,
        name,
        hgvsc_of(reference),
        hgvsc_of(candidate),
        Some(candidate),
    )
}

/// What a planted row states beside its notation.
struct Row {
    cdna_position: &'static str,
    given_ref_bases: usize,
    bam_edit: &'static str,
    strand: &'static str,
}

impl Default for Row {
    /// A deletion of eighteen genomic bases that holds nine model bases.
    fn default() -> Self {
        Self {
            cdna_position: "3385-3393",
            given_ref_bases: 18,
            bam_edit: "OK",
            strand: "1",
        }
    }
}

/// A default-format output of one deletion row on the made-up model with `notation`.
fn planted_text(notation: &str, row: &Row) -> String {
    let given_ref = match row.given_ref_bases {
        0 => String::new(),
        bases => format!(
            "GIVEN_REF={};USED_REF={};",
            "A".repeat(bases),
            "A".repeat(bases)
        ),
    };
    let bam_edit = match row.bam_edit {
        "" => String::new(),
        status => format!("BAM_EDIT={status};"),
    };
    let strand = match row.strand {
        "" => String::new(),
        strand => format!("STRAND={strand};"),
    };
    let extra = format!("IMPACT=MODERATE;{strand}{given_ref}{bam_edit}HGVSc={MODEL}:{notation}");
    default_text(&[format!(
        "del_1\t21:1001-1018\t-\tENSG1\t{MODEL}\tTranscript\tinframe_deletion\t{}\t3376-3384\t1126-1128\t-\t-\t-\t{extra}",
        row.cdna_position
    )])
}

/// The verdict of the planted row with `reference` against the same row with `candidate`.
fn planted_verdict(reference: &str, candidate: &str, row: &Row) -> String {
    verdict_of(
        &planted_text(reference, row),
        &planted_text(candidate, row),
        &[HGVSC_OFFSET],
    )
}

/// On the stored rows of two release records the reference's deletion opens where the
/// candidate's opens and closes nine positions later, the bases the model lacks: each of the
/// four pairs is the class's, and none is with the sides reversed.
#[test]
fn the_renumbered_notation_class_takes_the_stored_rows_of_a_deletion_across_deleted_bases() {
    let reference = parse_default(REFERENCE_ROWS).expect("the reference rows parse");
    let candidate = parse_default(CANDIDATE_ROWS).expect("the candidate rows parse");
    let notations: Vec<(&str, &str)> = reference
        .iter()
        .zip(&candidate)
        .map(|(e, a)| (hgvsc_of(e), hgvsc_of(a)))
        .collect();
    assert_eq!(
        notations,
        [
            ("NM_030582.4:c.3372_3388del", "NM_030582.4:c.3372_3379del"),
            ("NM_130444.3:c.4077_4093del", "NM_130444.3:c.4077_4084del"),
            ("NM_030582.4:c.3376_3393del", "NM_030582.4:c.3376_3384del"),
            ("NM_130444.3:c.4081_4098del", "NM_130444.3:c.4081_4089del"),
        ]
    );
    for (e, a) in reference.iter().zip(&candidate) {
        assert!(accepted(e, a), "{} {}", hgvsc_of(e), hgvsc_of(a));
        assert!(!accepted(a, e), "{} {}", hgvsc_of(a), hgvsc_of(e));
    }
}

/// Through the comparison the class sets the pair aside, counted once, and without the
/// class the pair is a difference.
#[test]
fn the_comparison_sets_aside_a_deletion_that_closes_later_by_the_deleted_bases() {
    let row = Row::default();
    let (written, model_span) = ("c.3376_3393del", "c.3376_3384del");
    assert_eq!(
        default_under(
            HGVSC_OFFSET,
            &planted_text(written, &row),
            &planted_text(model_span, &row)
        ),
        ("green".to_string(), once(HGVSC_OFFSET))
    );
    assert_eq!(
        verdict_of(
            &planted_text(written, &row),
            &planted_text(model_span, &row),
            &[]
        ),
        "red"
    );
}

/// A pair of another number of positions stays a difference unless every part of the shape
/// holds: two lone deletions over coding positions that open at one position, the
/// reference's covering as many positions as GIVEN_REF has bases and the candidate's as
/// many as the row's cDNA position, on a forward-strand row carrying `BAM_EDIT` `OK`.
#[test]
fn the_renumbered_notation_class_refuses_any_other_pair_of_different_lengths() {
    let row = Row::default();
    for (reference, candidate) in [
        // Not a deletion.
        ("c.3376_3393dup", "c.3376_3384dup"),
        ("c.3376_3393delinsA", "c.3376_3384delinsA"),
        ("c.3376_3393inv", "c.3376_3384inv"),
        // The two values open at different coordinates.
        ("c.3376_3393del", "c.3377_3385del"),
        ("c.3367_3384del", "c.3376_3384del"),
        // The reference covers another number of positions than the genomic bases.
        ("c.3376_3392del", "c.3376_3384del"),
        ("c.3376_3394del", "c.3376_3384del"),
        // The candidate covers another number of positions than the row's cDNA position.
        ("c.3376_3393del", "c.3376_3385del"),
        ("c.3376_3393del", "c.3376_3383del"),
        // The sides reversed.
        ("c.3376_3384del", "c.3376_3393del"),
        // Coordinates with an intron offset.
        ("c.3376+1_3376+18del", "c.3376+1_3376+9del"),
        // Positions before the coding start, across it or past the stop codon, and `n.`
        // numbering.
        ("c.-30_-13del", "c.-30_-22del"),
        ("c.-3_15del", "c.-3_6del"),
        ("c.*5_*22del", "c.*5_*13del"),
        ("n.3385_3402del", "n.3385_3393del"),
    ] {
        assert_eq!(
            planted_verdict(reference, candidate, &row),
            "red",
            "{reference} {candidate}"
        );
    }
    let (written, model_span) = ("c.3376_3393del", "c.3376_3384del");
    for (case, row) in [
        (
            "an edit that failed",
            Row {
                bam_edit: "FAILED",
                ..Row::default()
            },
        ),
        (
            "no BAM_EDIT",
            Row {
                bam_edit: "",
                ..Row::default()
            },
        ),
        (
            "a reverse-strand row",
            Row {
                strand: "-1",
                ..Row::default()
            },
        ),
        (
            "no STRAND",
            Row {
                strand: "",
                ..Row::default()
            },
        ),
        (
            "no GIVEN_REF",
            Row {
                given_ref_bases: 0,
                ..Row::default()
            },
        ),
        (
            "no cDNA position",
            Row {
                cdna_position: "-",
                ..Row::default()
            },
        ),
        (
            "a cDNA position open on one side",
            Row {
                cdna_position: "?-3393",
                ..Row::default()
            },
        ),
        (
            "a cDNA position of one base",
            Row {
                cdna_position: "3385",
                ..Row::default()
            },
        ),
    ] {
        assert_eq!(planted_verdict(written, model_span, &row), "red", "{case}");
    }
    // A deletion over bases the model holds and the genome lacks: the genomic bases are
    // the fewer.
    let held = Row {
        cdna_position: "3385-3402",
        given_ref_bases: 9,
        ..Row::default()
    };
    assert_eq!(planted_verdict(model_span, written, &held), "red");
    assert_eq!(planted_verdict(written, model_span, &held), "red");
}

/// The row's cDNA position states a count where both sides are numbers and the end does
/// not lie before the start, with or without the transcript's length after a slash
/// (`--total_length`). A side written with a digit that is not a decimal one (a
/// superscript), an end before the start and a span of no position state none: the pair
/// is refused, whatever the candidate's pair covers.
#[test]
fn the_cdna_position_states_a_count_as_two_numbers_in_order() {
    let written = "c.3376_3393del";
    let with_position = |cdna_position: &'static str| Row {
        cdna_position,
        ..Row::default()
    };
    assert_eq!(
        planted_verdict(written, "c.3376_3384del", &with_position("3385-3393/9000")),
        "green"
    );
    for (candidate, cdna_position) in [
        ("c.3376_3384del", "3385-339\u{b3}"),
        ("c.3376_3368del", "3393-3385"),
        ("c.3376_3375del", "3385-3384"),
    ] {
        assert_eq!(
            planted_verdict(written, candidate, &with_position(cdna_position)),
            "red",
            "{candidate} {cdna_position}"
        );
    }
}

/// The row is read under each format's own names: `STRAND`, `GIVEN_REF` and
/// `cDNA_position` in the tab and VCF formats, `strand`, `given_ref`, `cdna_start` and
/// `cdna_end` in the JSON format.
#[test]
fn the_deletion_across_deleted_bases_is_read_in_every_format() {
    let notations = ["c.3376_3393del", "c.3376_3384del", "c.3376_3385del"];
    // The class takes the second notation against the first and refuses the third.
    let holds = |entry: &dyn Fn(&str) -> Entry| {
        let [reference, model_span, wrong] = notations.map(entry);
        accepted(&reference, &model_span) && !accepted(&reference, &wrong)
    };
    let bases = "A".repeat(18);
    let header = tab_header(&[
        "IMPACT",
        "STRAND",
        "GIVEN_REF",
        "USED_REF",
        "BAM_EDIT",
        "HGVSc",
    ]);
    let tab_row = format!(
        "del_1\t21:1001-1018\t-\tENSG1\t{MODEL}\tTranscript\tinframe_deletion\t3385-3393\t3376-3384\t1126-1128\t-\t-\t-\tMODERATE\t1\t{bases}\t{bases}\tOK\t{MODEL}:"
    );
    assert!(holds(&|notation| {
        tab_entries(&format!("{header}\n{tab_row}{notation}\n")).remove(0)
    }));
    let inputs = vec![format!("21\t1000\tdel_1\tC{bases}\tC\t.\t.\t.")];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|STRAND|GIVEN_REF|USED_REF|BAM_EDIT|HGVSc",
    );
    let csq_row = format!(
        "{}CSQ=-|inframe_deletion|{MODEL}|Transcript|3385-3393|3376-3384|1|{bases}|{bases}|OK|{MODEL}:",
        record_columns(&inputs[0])
    );
    assert!(holds(&|notation| {
        transcript_entry(
            vcf_entries(&format!("{header}{csq_row}{notation}\n"), &inputs),
            MODEL,
        )
    }));
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("del_1")),
        ("allele_string", text(&format!("{bases}/-"))),
        ("most_severe_consequence", text("inframe_deletion")),
    ];
    // The entry of the JSON row with `notation`, the member `missing` left out.
    let json_entry = |notation: &str, missing: Option<&str>| {
        let consequence = object(&[
            ("variant_allele", text("-")),
            ("transcript_id", text(MODEL)),
            (
                "consequence_terms",
                Value::List(vec![text("inframe_deletion")]),
            ),
            ("strand", Value::Int(1)),
            ("given_ref", text(&bases)),
            ("used_ref", text(&bases)),
            ("bam_edit", text("OK")),
            ("cdna_start", Value::Int(3385)),
            ("cdna_end", Value::Int(3393)),
            ("hgvsc", text(&format!("{MODEL}:{notation}"))),
        ]);
        let consequence = match missing {
            Some(name) => without_member(&consequence, name),
            None => consequence,
        };
        transcript_entry(
            json_entries(&json_line(&record, consequence), &inputs),
            MODEL,
        )
    };
    assert!(holds(&|notation| json_entry(notation, None)));
    let model_span = json_entry(notations[1], None);
    for missing in ["cdna_start", "cdna_end", "given_ref", "bam_edit", "strand"] {
        assert!(
            !accepted(&json_entry(notations[0], Some(missing)), &model_span),
            "{missing}"
        );
    }
}
