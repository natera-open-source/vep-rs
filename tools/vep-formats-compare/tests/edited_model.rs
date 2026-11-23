// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The edited-model classes and the clamped span of a shifted deletion, on planted rows.

mod compare_common;
mod shifted_edited_common;

use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::consts::HGVSC_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::edited_model::cdna_pair_spans;
use vep_formats_compare::rules::edited_model::edited_model_clipped_alleles;
use vep_formats_compare::rules::edited_model::edited_model_seam_insertion_replacement;
use vep_formats_compare::rules::edited_model::shared_ends;
use vep_formats_compare::shared::ABSENT;

const HGVSC_OFFSET: &str = "refseq_edited_model_hgvsc_offset";
const TURNED_BASES: &str = "refseq_edited_model_shifted_insertion_bases";
const SEAM: &str = "refseq_edited_model_seam_insertion_hgvsc";
const UNMINIMISED: &str = "refseq_edited_model_unminimised_allele_hgvsc";
const HGVS_OFFSET: &str = "refseq_edited_model_hgvs_offset";
const CLAMPED_SPAN: &str = "hgvs_shifted_deletion_clamped_span";

/// The entry's own HGVSc value, `-` for an entry without one: the form a blank value is
/// handed to a rule in.
fn hgvsc_of(entry: &Entry) -> &str {
    HGVSC_FIELDS
        .iter()
        .find_map(|field| entry.fields.get(*field))
        .map_or("-", String::as_str)
}

fn default_entries(text: &str) -> Vec<Entry> {
    parse_default(text).expect("the text parses")
}

/// An HGVSc on a BAM_EDIT OK row whose two notations name the same reference and differ in
/// their coordinate tokens alone is the class; a row without BAM_EDIT OK, a changed allele
/// text or another reference stays a difference.
#[test]
fn the_renumbered_notation_class_accepts_other_coordinates_on_an_edited_model_row() {
    let class = [HGVSC_OFFSET];
    let extra = "IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;BAM_EDIT=OK;HGVSc=NM_030582.4:c.3546_3553+23del;HGVS_OFFSET=7";
    let reference = default_text(&[default_row(
        "rs2037126010",
        "21:45505273-45505303",
        "-",
        "NM_030582.4",
        "splice_donor_variant,intron_variant",
        extra,
    )]);
    let (written, renumbered) = ("c.3546_3553+23del", "c.3544+2_3544+32del");
    assert_eq!(
        default_under(
            HGVSC_OFFSET,
            &reference,
            &reference.replace(written, renumbered)
        ),
        ("green".to_string(), once(HGVSC_OFFSET))
    );
    assert_eq!(swapped(&reference, written, renumbered, &[]), "red");
    let unedited = reference.replace("BAM_EDIT=OK;", "");
    assert_eq!(swapped(&unedited, written, renumbered, &class), "red");
    assert_eq!(
        swapped(&reference, written, "c.3546_3553+23dup", &class),
        "red"
    );
    assert_eq!(
        swapped(
            &reference,
            "NM_030582.4:c.3546",
            "NM_030582.5:c.3546",
            &class
        ),
        "red"
    );
    let snv = reference.replace(written, "c.56+4A>G");
    assert_eq!(swapped(&snv, "c.56+4A>G", "c.57+4A>G", &class), "green");
    assert_eq!(swapped(&snv, "c.56+4A>G", "c.56+4A>T", &class), "red");
}

/// A pair of coordinates that covers another number of positions in the candidate than in
/// the reference stays a difference; the deletion across bases a model lacks, the one pair
/// of other lengths the class takes, is held by `edited_model_deleted_bases.rs`. A pair
/// that mixes a plain position with an offset, as the class's own row does, names two
/// boundaries or runs over the stop codon is not measured and is accepted on its coordinate
/// tokens alone.
#[test]
fn the_renumbered_notation_class_refuses_values_that_cover_different_numbers_of_positions() {
    let extra = "IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK;HGVSc=NM_030582.4:c.69_74dup";
    let reference = default_text(&[default_row(
        "rs1",
        "21:45505273-45505274",
        "ACGACG",
        "NM_030582.4",
        "inframe_insertion",
        extra,
    )]);
    let read = |expected: &str, actual: &str| {
        let left = reference.replace("c.69_74dup", expected);
        swapped(&left, expected, actual, &[HGVSC_OFFSET])
    };
    for (expected, actual) in [
        ("c.69_74dup", "c.72_77dup"),
        ("c.-3_2del", "c.1_5del"),
        ("c.*5_*9del", "c.*8_*12del"),
        ("c.636+18_636+19insA", "c.639+18_639+19insA"),
        ("c.3546_3553+23del", "c.3544+2_3544+32del"),
        ("c.100+5_101-3del", "c.103+5_104-3del"),
        ("c.1499_*2del", "c.1502_*2del"),
    ] {
        assert_eq!(read(expected, actual), "green", "{expected} {actual}");
    }
    for (expected, actual) in [
        ("c.69_74dup", "c.75_77dup"),
        ("c.74_75insA", "c.74_78insA"),
        ("c.-3_2del", "c.1_4del"),
        ("c.636+18_636+19insA", "c.75_77insA"),
        ("c.*5_*9del", "c.*8_*9del"),
    ] {
        assert_eq!(read(expected, actual), "red", "{expected} {actual}");
    }
}

/// The count a pair of coordinates states is read from two plain positions (1 follows -1)
/// or from two offsets off one exon boundary; any other pair states none, and a notation of
/// one coordinate holds no pair.
#[test]
fn the_positions_between_two_coordinates_are_counted_where_the_pair_states_them() {
    assert_eq!(cdna_pair_spans("c.-3_2del"), [Some(4)]);
    assert_eq!(cdna_pair_spans("c.636+18_636+19insA"), [Some(1)]);
    assert_eq!(cdna_pair_spans("c.636-9_636-5del"), [Some(4)]);
    assert_eq!(cdna_pair_spans("c.3546_3553+23del"), [None]);
    assert_eq!(cdna_pair_spans("c.100+5_101-3del"), [None]);
    assert_eq!(cdna_pair_spans("c.1499_*2del"), [None]);
    assert_eq!(cdna_pair_spans("c.56+4A>G"), []);
}

/// Under --hgvs, an HGVSc on a BAM_EDIT OK insertion row whose HGVS_OFFSET is one non-zero
/// integer on both sides, the reference's insertion carrying the row's Allele on the
/// transcript's strand and the candidate's the same coordinates with those bases turned 3' by
/// the offset, is the class, counted once per row in every format; a row without BAM_EDIT OK,
/// a deletion's or a substitution's row, another transcript, other coordinates, a duplication,
/// bases that are not the allele's, a turn that is not the offset's, a missing, zero or
/// differing offset, and the reversed sides stay differences.
#[test]
fn the_turned_bases_class_accepts_the_turned_bases_on_an_edited_model_row() {
    let class = [TURNED_BASES];
    let extra = "IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_006933.7:c.204_205insGA;HGVSp=NP_008864.4:p.Ala69ArgfsTer64;HGVS_OFFSET=1";
    let row = default_row(
        "ins_cds",
        "21:35467700-35467701",
        "GA",
        "NM_006933.7",
        "frameshift_variant",
        extra,
    );
    let reference = default_text(&[row.replacen(
        "\t-\t-\t-\t-\t-\t-\t",
        "\t707-708\t203-204\t68\tL/LX\tctg/ctGAg\t-\t",
        1,
    )]);
    let written = "c.204_205insGA";
    let turned = |notation: &str| reference.replace(written, notation);

    let (report, documented) = compare_default(&reference, &turned("c.204_205insAG"), &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(TURNED_BASES));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &masked_once(TURNED_BASES)
    );
    assert_eq!(swapped(&reference, written, "c.204_205insAG", &[]), "red");
    assert_eq!(
        swapped(&reference, written, "c.204_205insAG", &[HGVSC_OFFSET]),
        "red"
    );
    for notation in [
        "c.203_204insAG",
        "c.204_205insCA",
        "c.204_205insGAG",
        "c.204_205dup",
        "n.204_205insAG",
        "-",
    ] {
        assert_eq!(
            swapped(&reference, written, notation, &class),
            "red",
            "{notation}"
        );
    }
    assert_eq!(
        verdict_of(&turned("c.204_205insAG"), &reference, &class),
        "red"
    );
    assert_eq!(
        verdict_of(
            &reference,
            &turned("c.204_205insAG").replace("HGVSc=NM_006933.7", "HGVSc=NM_006933.8"),
            &class
        ),
        "red"
    );
    let unedited = reference.replace("BAM_EDIT=OK;", "");
    assert_eq!(swapped(&unedited, written, "c.204_205insAG", &class), "red");
    // The reference's bases are the allele's; a notation spelling other bases is not the shape.
    let other_bases = reference.replace(written, "c.204_205insGG");
    assert_eq!(
        swapped(&other_bases, "c.204_205insGG", "c.204_205insAG", &class),
        "red"
    );
    // The turn is the offset modulo the allele's length, towards the transcript's 3' end.
    let three = reference
        .replace("\tGA\t", "\tGAT\t")
        .replace(written, "c.204_205insGAT");
    for (offset, notation, expected) in [
        ("1", "c.204_205insATG", "green"),
        ("2", "c.204_205insTGA", "green"),
        ("4", "c.204_205insATG", "green"),
        ("1", "c.204_205insTGA", "red"),
        ("3", "c.204_205insATG", "red"),
    ] {
        let moved = three.replace("HGVS_OFFSET=1", &format!("HGVS_OFFSET={offset}"));
        assert_eq!(
            swapped(&moved, "c.204_205insGAT", notation, &class),
            expected,
            "{offset} {notation}"
        );
    }
    for offset in ["0", "7.5", ""] {
        let field = if offset.is_empty() {
            String::new()
        } else {
            format!(";HGVS_OFFSET={offset}")
        };
        let without = reference.replace(";HGVS_OFFSET=1", &field);
        assert_eq!(
            swapped(&without, written, "c.204_205insAG", &class),
            "red",
            "{offset}"
        );
    }
    assert_eq!(
        verdict_of(
            &reference,
            &turned("c.204_205insAG").replace("HGVS_OFFSET=1", "HGVS_OFFSET=2"),
            &class
        ),
        "red"
    );
    let deletion = reference
        .replace("\tGA\t", "\t-\t")
        .replace(written, "c.204_205del");
    assert_eq!(
        swapped(&deletion, "c.204_205del", "c.204_205delinsAG", &class),
        "red"
    );
    let substitution = reference.replace("21:35467700-35467701", "21:35467700");
    assert_eq!(
        swapped(&substitution, written, "c.204_205insAG", &class),
        "red"
    );
    // On a reverse-strand transcript the allele is read reverse-complemented; the allele as
    // written there is not the shape, nor is the reverse complement on a forward strand.
    let reverse = reference
        .replace("STRAND=1", "STRAND=-1")
        .replace("\tGA\t", "\tTC\t")
        .replace("HGVS_OFFSET=1", "HGVS_OFFSET=-1");
    assert_eq!(
        swapped(&reverse, written, "c.204_205insAG", &class),
        "green"
    );
    let as_written = reverse.replace(written, "c.204_205insTC");
    assert_eq!(
        swapped(&as_written, "c.204_205insTC", "c.204_205insCT", &class),
        "red"
    );
    let complemented = reference.replace(written, "c.204_205insTC");
    assert_eq!(
        swapped(&complemented, "c.204_205insTC", "c.204_205insCT", &class),
        "red"
    );
    // A row without STRAND admits either reading of the allele.
    let unstranded = as_written.replace("STRAND=-1;", "");
    assert_eq!(
        swapped(&unstranded, "c.204_205insTC", "c.204_205insCT", &class),
        "green"
    );
    // An allele whose turn equals its reverse complement keeps the sides apart: the
    // reference must spell the allele on the transcript's strand.
    let palindromic = reference
        .replace("\tGA\t", "\tATCG\t")
        .replace(written, "c.204_205insATCG")
        .replace("HGVS_OFFSET=1", "HGVS_OFFSET=2");
    let turned_palindromic = palindromic.replace("c.204_205insATCG", "c.204_205insCGAT");
    assert_eq!(
        verdict_of(&palindromic, &turned_palindromic, &class),
        "green"
    );
    assert_eq!(verdict_of(&turned_palindromic, &palindromic, &class), "red");
    // Intronic and 3' UTR coordinates are insertions too.
    for coordinates in ["c.204+5_204+6", "c.*12_*13", "c.-3_-2"] {
        let moved = reference.replace(written, &format!("{coordinates}insGA"));
        assert_eq!(
            swapped(
                &moved,
                &format!("{coordinates}insGA"),
                &format!("{coordinates}insAG"),
                &class
            ),
            "green",
            "{coordinates}"
        );
    }
    // The tab, VCF and JSON formats count the row once.
    let header = tab_header(&["IMPACT", "STRAND", "BAM_EDIT", "HGVSc", "HGVS_OFFSET"]);
    let tab_row = "ins_cds\t21:35467700-35467701\tGA\t6526\tNM_006933.7\tTranscript\tframeshift_variant\t707-708\t203-204\t68\tL/LX\tctg/ctGAg\t-\tHIGH\t1\tOK\tNM_006933.7:c.204_205ins";
    let tab = |tail: &str| tab_entries(&format!("{header}\n{tab_row}{tail}\n"));
    assert_eq!(
        under_class(TURNED_BASES, &tab("GA\t1"), &tab("AG\t1")),
        ("green".to_string(), once(TURNED_BASES))
    );
    let inputs = vec!["21\t35467700\tins_cds\tT\tTGA\t.\t.\t.".to_string()];
    let header =
        csq_header("Allele|Consequence|Feature|Feature_type|STRAND|BAM_EDIT|HGVSc|HGVS_OFFSET");
    let csq_row = format!(
        "{}CSQ=GA|frameshift_variant|NM_006933.7|Transcript|1|OK|NM_006933.7:c.204_205ins",
        record_columns(&inputs[0])
    );
    let vcf = |tail: &str| vcf_entries(&format!("{header}{csq_row}{tail}\n"), &inputs);
    assert_eq!(
        under_class(TURNED_BASES, &vcf("GA|1"), &vcf("AG|1")),
        ("green".to_string(), once(TURNED_BASES))
    );
    assert_eq!(
        under_class(TURNED_BASES, &vcf("GA|1"), &vcf("AG|")).0,
        "red"
    );
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("ins_cds")),
        ("allele_string", text("-/GA")),
        ("most_severe_consequence", text("frameshift_variant")),
    ];
    let object_of = |hgvsc: &str, fields: &[(&str, Value)]| {
        let consequence = object(&[
            ("variant_allele", text("GA")),
            ("transcript_id", text("NM_006933.7")),
            (
                "consequence_terms",
                Value::List(vec![text("frameshift_variant")]),
            ),
            ("strand", Value::Int(1)),
            ("bam_edit", text("OK")),
            ("hgvsc", text(&format!("NM_006933.7:{hgvsc}"))),
            ("hgvs_offset", Value::Int(1)),
        ]);
        json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        )
    };
    assert_eq!(
        under_class(
            TURNED_BASES,
            &object_of("c.204_205insGA", &[]),
            &object_of("c.204_205insAG", &[])
        ),
        ("green".to_string(), once(TURNED_BASES))
    );
    assert_eq!(
        under_class(
            TURNED_BASES,
            &object_of("c.204_205insGA", &[]),
            &object_of("c.204_205insAG", &[("bam_edit", text("FAILED"))])
        )
        .0,
        "red"
    );
}

const SEAM_MODEL: &str = "NM_990001.2";

/// What a row of the made-up model `NM_990001.2` prints beside its allele and notation: an
/// insertion between two genomic bases, by default the two around the coding bases 9 to 11
/// (three bases of the codons of coding bases 7 to 12), which the model holds and the genome
/// lacks; `name` is the sequence name of its notation.
#[derive(Clone, Copy)]
struct SeamRow<'a> {
    cds: &'a str,
    protein: &'a str,
    extra: &'a str,
    location: &'a str,
    name: &'a str,
}

impl Default for SeamRow<'_> {
    fn default() -> Self {
        SeamRow {
            cds: "9-11",
            protein: "3-4",
            extra: "IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK",
            location: "1:157-158",
            name: SEAM_MODEL,
        }
    }
}

/// The entry of one row of the made-up model: `hgvsc` is the notation after the sequence
/// name, `-` for a row without one.
fn seam_entry(allele: &str, hgvsc: &str, codons: &str, row: SeamRow) -> Entry {
    let extra = if hgvsc == "-" {
        row.extra.to_string()
    } else {
        format!("{};HGVSc={}:{hgvsc}", row.extra, row.name)
    };
    let planted = default_row(
        "seam_ins",
        row.location,
        allele,
        SEAM_MODEL,
        "missense_variant",
        &extra,
    );
    default_entry(&planted.replacen(
        "\t-\t-\t-\t-\t-\t-\t",
        &format!(
            "\t59-61\t{}\t{}\tGP/GL\t{codons}\t-\t",
            row.cds, row.protein
        ),
        1,
    ))
}

/// The rule on the two entries' own HGVSc values, a missing one handed over as `-`.
fn seam_rule(reference: &Entry, candidate: &Entry) -> bool {
    edited_model_seam_insertion_replacement(
        reference,
        Some(candidate),
        hgvsc_of(reference),
        hgvsc_of(candidate),
    )
}

/// Whether the rule accepts the notation `candidate` on a row of `allele` whose reference
/// notation is `reference`; `row` changes the two sides' rows alike.
fn seam_accepts_against(
    allele: &str,
    codons: &str,
    candidate: &str,
    reference: &str,
    row: SeamRow,
) -> bool {
    seam_rule(
        &seam_entry(allele, reference, codons, row),
        &seam_entry(allele, candidate, codons, row),
    )
}

/// `seam_accepts_against` the `ins` of the allele between coding bases 11 and 12.
fn seam_accepts(allele: &str, codons: &str, candidate: &str, row: SeamRow) -> bool {
    seam_accepts_against(
        allele,
        codons,
        candidate,
        &format!("c.11_12ins{allele}"),
        row,
    )
}

/// Against the reference's `ins` of the allele, at plain coordinates or with intron
/// offsets, the candidate's HGVSc of an insertion standing on the three bases an edited
/// model holds at coding bases 9 to 11 is accepted when, applied to the row's reference
/// codons, it gives the alternate codons: a `delins` of the span, a deletion of the bases
/// the allele leaves out, a substitution of the one base it changes, an insertion beside the
/// bases it keeps, a duplication for the bases twice, an inversion for their reverse
/// complement, and no notation for the bases themselves. Another notation of the same change
/// is accepted beside the first, among them a deletion or a duplication the 3' rule writes
/// further on inside the codons. The two sides reversed stay a difference.
#[test]
fn the_seam_class_accepts_each_notation_of_the_replaced_bases() {
    let plain = SeamRow::default();
    let forms: [(&str, &[&str]); 14] = [
        ("G", &["c.9_11delinsG", "c.8_11delinsGG"]),
        ("GGG", &["c.9_11delinsGGG"]),
        ("TA", &["c.9_11delinsTA"]),
        ("TAC", &["c.9_10delinsTA", "c.9_11delinsTAC"]),
        ("AGGC", &["c.10delinsGG", "c.10_11delinsGGC"]),
        ("A", &["c.10_11del", "c.9_11delinsA"]),
        ("AC", &["c.10del", "c.11del"]),
        ("GCC", &["c.9A>G", "c.9delinsG"]),
        ("ATC", &["c.10C>T"]),
        ("ACCG", &["c.11_12insG"]),
        ("GACC", &["c.8_9insG", "c.8dup"]),
        (
            "ACCACC",
            &["c.9_11dup", "c.10_12dup", "c.11_12insACC", "c.9_11[2]"],
        ),
        ("GGT", &["c.9_11inv", "c.9_11delinsGGT"]),
        ("ACC", &["-"]),
    ];
    for (allele, notations) in forms {
        let codons = format!("ggacca/gg{allele}a");
        for notation in notations {
            assert!(
                seam_accepts(allele, &codons, notation, plain),
                "{allele} {notation}"
            );
        }
        for coordinates in ["c.2034+2_2034+3", "c.-3_-2", "c.*12_*13"] {
            assert!(
                seam_accepts_against(
                    allele,
                    &codons,
                    notations[0],
                    &format!("{coordinates}ins{allele}"),
                    plain
                ),
                "{allele} {coordinates}"
            );
        }
        assert!(
            !seam_rule(
                &seam_entry(allele, notations[0], &codons, plain),
                &seam_entry(allele, &format!("c.11_12ins{allele}"), &codons, plain)
            ),
            "{allele}"
        );
    }
    // The 3' rule writes a deletion further on while the next base is the one deleted: one
    // base of a run of three is deleted at any of its three positions.
    let two_bases = SeamRow {
        cds: "9-10",
        ..plain
    };
    for notation in ["c.10del", "c.11del", "c.12del"] {
        assert!(
            seam_accepts("C", "ggcaaa/ggCaa", notation, two_bases),
            "{notation}"
        );
    }
    // An inversion of the bases the two sequences do not share.
    let four_bases = SeamRow {
        cds: "7-10",
        ..plain
    };
    assert!(seam_accepts(
        "GTTC",
        "gaacca/GTTCca",
        "c.8_9inv",
        four_bases
    ));
    // The reference's own Codons, which read other bases on such a row, are not the rule's.
    assert!(seam_rule(
        &seam_entry("G", "c.11_12insG", "gagtac/gaGc", plain),
        &seam_entry("G", "c.9_11delinsG", "ggacca/ggGa", plain)
    ));
}

/// On a seam row the candidate's HGVSc is accepted for what it says of the row's codons and
/// for nothing else: the reference's own notation, an `ins` between positions that are not
/// neighbours, a notation that gives another sequence, one that names a base it does not
/// find, one that reaches outside the codons, a position that is not a plain coding
/// position, a text that is none of the forms, another numbering or sequence name, any
/// notation for the model's own bases and none for any other allele stay differences.
#[test]
fn the_seam_class_keeps_a_wrong_value_on_a_seam_row_a_difference() {
    let plain = SeamRow::default();
    let codons = "ggacca/ggGa";
    assert!(seam_accepts("G", codons, "c.9_11delinsG", plain));
    for wrong in [
        "c.11_12insG",
        "c.8_12insG",
        "c.9_11delinsT",
        "c.9_10delinsG",
        "c.9_11del",
        "c.9_11dup",
        "c.9_11inv",
        "c.9A>G",
        "c.9_11delinsg",
        "c.9_11delins",
        "c.9_11",
        "n.9_11delinsG",
        "-",
    ] {
        assert!(!seam_accepts("G", codons, wrong, plain), "{wrong}");
    }
    let other_name = SeamRow {
        name: "NM_990001.3",
        ..plain
    };
    assert!(!seam_rule(
        &seam_entry("G", "c.11_12insG", codons, plain),
        &seam_entry("G", "c.9_11delinsG", codons, other_name)
    ));
    let wrong_values: [(&str, &str, &[&str]); 6] = [
        // A deletion of other bases, and one written end before start.
        (
            "A",
            "c.10_11del",
            &["c.9_10del", "c.11_12del", "c.11_10del", "c.8_12insA"],
        ),
        // A substitution naming another base than the one at its position, or two positions.
        (
            "GCC",
            "c.9A>G",
            &["c.9C>G", "c.10A>G", "c.9A>T", "c.9_11A>G", "c.9_9A>G"],
        ),
        // An `ins` of the one new base between positions that are not neighbours, at one
        // position or at another place, and a `delins` written end before start.
        (
            "ACCG",
            "c.11_12insG",
            &[
                "c.8_12insG",
                "c.9_12insG",
                "c.11insG",
                "c.8_12insACCG",
                "c.12_13insG",
                "c.12_11delinsG",
            ],
        ),
        // A duplication that needs the base at 13, past the codons' last, and a repeat count
        // of more copies than the alternate codons hold.
        (
            "ACCACC",
            "c.10_12dup",
            &[
                "c.11_13dup",
                "c.12_14dup",
                "c.9_10dup",
                "c.9_11[3]",
                "c.9_11[999999999999999999]",
            ],
        ),
        ("GGT", "c.9_11inv", &["c.9_10inv", "c.10_11inv", "c.9inv"]),
        // The model's own bases change nothing, so any notation for them is a difference.
        (
            "ACC",
            "-",
            &[
                "c.9_11=",
                "c.9_11delinsACC",
                "c.9_11dup",
                "c.9A>A",
                "c.9_11[1]",
            ],
        ),
    ];
    for (allele, described, wrong_notations) in wrong_values {
        let codons = format!("ggacca/gg{allele}a");
        assert!(seam_accepts(allele, &codons, described, plain), "{allele}");
        for wrong in wrong_notations {
            assert!(
                !seam_accepts(allele, &codons, wrong, plain),
                "{allele} {wrong}"
            );
        }
    }
    // A deletion the 3' rule would write at 13, past the codons' last base, and positions
    // with an intron offset, past the stop codon or before the start codon.
    let two_bases = SeamRow {
        cds: "9-10",
        ..plain
    };
    for wrong in ["c.13del", "c.12_13del", "c.10+1del", "c.*1del", "c.-10del"] {
        assert!(
            !seam_accepts("C", "ggcaaa/ggCaa", wrong, two_bases),
            "{wrong}"
        );
    }
    // Before coding base 1 an insertion names `c.-1`, a base outside the coding sequence.
    let first_codons = SeamRow {
        cds: "1-3",
        protein: "1-2",
        ..plain
    };
    assert!(seam_accepts(
        "GATG",
        "atggca/GATGgca",
        "c.1_3delinsGATG",
        first_codons
    ));
    for wrong in ["c.-1_1insG", "c.0_1insG"] {
        assert!(
            !seam_accepts("GATG", "atggca/GATGgca", wrong, first_codons),
            "{wrong}"
        );
    }
    // An `ins` beyond the base after the codons' last or before the base before their
    // first, a `delins` reaching before their first base and an inversion at one position.
    let one = SeamRow {
        cds: "1059",
        protein: "353",
        ..plain
    };
    for (allele, codons, wrong, reference) in [
        ("GA", "ctg/ctGA", "c.1061_1062insA", "c.1060_1061insGA"),
        ("AG", "ctg/ctAG", "c.1055_1056insA", "c.1060_1061insAG"),
        ("T", "ctg/ctT", "c.1056_1059delinsT", "c.1060_1061insT"),
        ("C", "ctg/ctC", "c.1059inv", "c.1060_1061insC"),
    ] {
        assert!(
            !seam_accepts_against(allele, codons, wrong, reference, one),
            "{wrong}"
        );
    }
}

/// The seam is read from the row, so a row that does not show one keeps its difference
/// whatever the candidate writes: an ordinary insertion on an edited model, an insertion one
/// base before or after the model's bases, a seam row whose notation is a duplication of the
/// base after it, a row without BAM_EDIT, a row without a CDS position, a protein position or
/// Codons, codons that do not hold the allele in the place of the span's bases, a reference
/// notation that is not an `ins` of the row's allele in `c.` numbering on the row's
/// transcript, a deletion's or a substitution's row, and values that are not the two rows'
/// own HGVSc.
#[test]
fn the_seam_class_keeps_a_value_on_a_row_off_a_seam_a_difference() {
    let plain = SeamRow::default();
    let codons = "ggacca/ggGa";
    let described = "c.9_11delinsG";
    assert!(seam_accepts("G", codons, described, plain));
    let failed = SeamRow {
        extra: "IMPACT=MODERATE;STRAND=1;BAM_EDIT=FAILED",
        ..plain
    };
    assert!(seam_accepts("G", codons, described, failed));
    // An ordinary insertion between coding bases 74 and 75, then one before and one after
    // the model's bases at 9 to 11.
    let ordinary = SeamRow {
        cds: "74-75",
        protein: "25",
        ..plain
    };
    for notation in ["c.74_75delinsG", "c.74_75insG", "-"] {
        assert!(
            !seam_accepts_against("G", "gat/gaGt", notation, "c.77_78insG", ordinary),
            "{notation}"
        );
    }
    let before = SeamRow {
        cds: "7-8",
        protein: "3",
        ..plain
    };
    for notation in ["c.7_8insG", "c.7_8delinsG"] {
        assert!(
            !seam_accepts_against("G", "gga/gGga", notation, "c.10_11insG", before),
            "{notation}"
        );
    }
    let after = SeamRow {
        cds: "12-13",
        protein: "4-5",
        ..plain
    };
    for notation in ["c.12_13insG", "c.12delinsAG"] {
        assert!(
            !seam_accepts_against("G", "cca/ccaG", notation, "c.15_16insG", after),
            "{notation}"
        );
    }
    // On the seam, a duplication of the base after the model's bases gives another sequence
    // than the row's alternate codons.
    assert!(!seam_accepts("A", "ggacca/ggAa", "c.12dup", plain));
    // A mitochondrial model's row carries no BAM_EDIT; the field on one side alone is not
    // the shape.
    let unedited = SeamRow {
        extra: "IMPACT=MODERATE;STRAND=1",
        ..plain
    };
    assert!(!seam_accepts("G", codons, described, unedited));
    assert!(!seam_rule(
        &seam_entry("G", "c.11_12insG", codons, plain),
        &seam_entry("G", described, codons, unedited)
    ));
    assert!(!seam_rule(
        &seam_entry("G", "c.11_12insG", codons, unedited),
        &seam_entry("G", described, codons, plain)
    ));
    // A span in an untranslated region or on a non-coding model: the row carries no CDS
    // position and no Codons to read the model's bases from.
    let unplaced = SeamRow {
        cds: "-",
        protein: "-",
        ..plain
    };
    assert!(!seam_accepts_against(
        "G",
        "-",
        "c.-44_-42delinsG",
        "c.-41_-40insG",
        unplaced
    ));
    assert!(!seam_accepts_against(
        "G",
        "-",
        "n.21_23delinsG",
        "n.24_25insG",
        unplaced
    ));
    assert!(!seam_accepts("G", "-", described, plain));
    assert!(!seam_rule(
        &seam_entry("G", "c.11_12insG", codons, plain),
        &seam_entry("G", described, "-", plain)
    ));
    for position in [
        SeamRow {
            protein: "-",
            ..plain
        },
        SeamRow {
            protein: "?-4",
            ..plain
        },
        SeamRow { cds: "-", ..plain },
        SeamRow {
            cds: "?-11",
            ..plain
        },
        SeamRow {
            cds: "9-?",
            ..plain
        },
        SeamRow {
            cds: "11-9",
            ..plain
        },
    ] {
        assert!(
            !seam_accepts("G", codons, described, position),
            "{} {}",
            position.cds,
            position.protein
        );
    }
    // The span lies outside the codons the protein position names, runs past their last base
    // or starts before their first.
    for protein in ["4", "1-4"] {
        assert!(
            !seam_accepts("G", codons, described, SeamRow { protein, ..plain }),
            "{protein}"
        );
    }
    for cds in ["9-11", "6-9"] {
        let one_codon = SeamRow {
            cds,
            protein: "3",
            ..plain
        };
        assert!(!seam_accepts("G", "gga/ggG", "c.9A>G", one_codon), "{cds}");
    }
    // Codons without reference bases, as an insertion between two codons prints them, hold
    // no model bases.
    let between_codons = SeamRow {
        cds: "10",
        protein: "4",
        ..plain
    };
    assert!(!seam_accepts("G", "-/G", "c.10delinsG", between_codons));
    // The alternate codons do not hold the allele in the place of the span's bases.
    for other in ["ggacca/ggTa", "ggacca/ggaGcca", "ggacca/gGa", "ggacca"] {
        assert!(!seam_accepts("G", other, described, plain), "{other}");
    }
    // The reference's notation is an `ins` of the row's allele, in `c.` numbering, on the
    // row's transcript.
    for reference in [
        "c.11_12insT",
        "c.11_12insGG",
        "c.11_12dup",
        "c.11_12delinsG",
        "c.11insG",
        "n.11_12insG",
        "-",
    ] {
        assert!(
            !seam_accepts_against("G", codons, described, reference, plain),
            "{reference}"
        );
    }
    let other_model = SeamRow {
        name: "NM_990002.1",
        ..plain
    };
    assert!(!seam_rule(
        &seam_entry("G", "c.11_12insG", codons, other_model),
        &seam_entry("G", described, codons, other_model)
    ));
    // A deletion's row and a substitution's row.
    assert!(!seam_accepts_against(
        "-",
        "ggacca/gga",
        "c.9_11del",
        "c.11_12del",
        plain
    ));
    let one_base = SeamRow {
        location: "1:157",
        ..plain
    };
    assert!(!seam_accepts("G", codons, described, one_base));
    // The same two notations under another key are not the rows' HGVSc, nor is the Extra
    // text that carries them, and a row the candidate does not write has no notation to read.
    let edited = "IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK";
    let under_another_key = |notation: &str| {
        let extra = format!("{edited};OTHER={SEAM_MODEL}:{notation}");
        seam_entry(
            "G",
            "-",
            codons,
            SeamRow {
                extra: &extra,
                ..plain
            },
        )
    };
    let reference = under_another_key("c.11_12insG");
    let candidate = under_another_key(described);
    assert!(!edited_model_seam_insertion_replacement(
        &reference,
        Some(&candidate),
        &reference.fields["OTHER"],
        &candidate.fields["OTHER"]
    ));
    let reference = seam_entry("G", "c.11_12insG", codons, plain);
    let candidate = seam_entry("G", described, codons, plain);
    assert!(edited_model_seam_insertion_replacement(
        &reference,
        Some(&candidate),
        &reference.fields["HGVSc"],
        &candidate.fields["HGVSc"]
    ));
    assert!(!edited_model_seam_insertion_replacement(
        &reference,
        Some(&candidate),
        &reference.fields["Extra"],
        &candidate.fields["Extra"]
    ));
    assert!(!edited_model_seam_insertion_replacement(
        &reference,
        None,
        &reference.fields["HGVSc"],
        &candidate.fields["HGVSc"]
    ));
}

/// The tab format and the VCF format's CSQ carry the span, the protein position and the
/// codons under the default format's names, the JSON format as `cds_start`, `cds_end`,
/// `protein_start` and `codons`: the rule accepts the same row in each and keeps the same
/// wrong value a difference.
#[test]
fn the_seam_class_reads_the_row_in_every_format() {
    let notations = ["c.11_12insG", "c.9_11delinsG", "c.9_11delinsT"];
    // The rule takes the second notation against the first and refuses the third.
    let holds = |entry: &dyn Fn(&str) -> Entry| {
        let [reference, derived, wrong] = notations.map(entry);
        seam_rule(&reference, &derived) && !seam_rule(&reference, &wrong)
    };
    let header = tab_header(&["IMPACT", "STRAND", "BAM_EDIT", "HGVSc"]);
    let tab_row = format!(
        "seam_ins\t1:157-158\tG\tENSG1\t{SEAM_MODEL}\tTranscript\tmissense_variant\t59-61\t9-11\t3-4\tGP/GL\tggacca/ggGa\t-\tMODERATE\t1\tOK\t{SEAM_MODEL}:"
    );
    assert!(holds(&|notation| {
        tab_entries(&format!("{header}\n{tab_row}{notation}\n")).remove(0)
    }));
    let inputs = vec!["1\t157\tseam_ins\tT\tTG\t.\t.\t.".to_string()];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|CDS_position|Protein_position|Codons|STRAND|BAM_EDIT|HGVSc",
    );
    let csq_row = format!(
        "{}CSQ=G|missense_variant|{SEAM_MODEL}|Transcript|9-11|3-4|ggacca/ggGa|1|OK|{SEAM_MODEL}:",
        record_columns(&inputs[0])
    );
    assert!(holds(&|notation| {
        transcript_entry(
            vcf_entries(&format!("{header}{csq_row}{notation}\n"), &inputs),
            SEAM_MODEL,
        )
    }));
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("seam_ins")),
        ("allele_string", text("-/G")),
        ("most_severe_consequence", text("missense_variant")),
    ];
    // The entry of the JSON row with `notation`, the member `missing` left out.
    let json_entry = |notation: &str, missing: Option<&str>| {
        let consequence = object(&[
            ("variant_allele", text("G")),
            ("transcript_id", text(SEAM_MODEL)),
            (
                "consequence_terms",
                Value::List(vec![text("missense_variant")]),
            ),
            ("strand", Value::Int(1)),
            ("bam_edit", text("OK")),
            ("cds_start", Value::Int(9)),
            ("cds_end", Value::Int(11)),
            ("protein_start", Value::Int(3)),
            ("protein_end", Value::Int(4)),
            ("codons", text("ggacca/ggGa")),
            ("hgvsc", text(&format!("{SEAM_MODEL}:{notation}"))),
        ]);
        let consequence = match missing {
            Some(name) => without_member(&consequence, name),
            None => consequence,
        };
        transcript_entry(
            json_entries(&json_line(&record, consequence), &inputs),
            SEAM_MODEL,
        )
    };
    assert!(holds(&|notation| json_entry(notation, None)));
    let reference = json_entry(notations[0], None);
    for missing in ["cds_end", "protein_start", "codons", "bam_edit"] {
        assert!(
            !seam_rule(&reference, &json_entry("c.9_11delinsG", Some(missing))),
            "{missing}"
        );
    }
}

/// Through the comparison itself the class sets aside the HGVSc of the made-up model's seam
/// row, counted once, and nothing where its condition does not hold: the same pair on a row
/// without BAM_EDIT, the candidate's notation missing where the allele is not the model's
/// bases, and a notation that gives another sequence stay differences.
#[test]
fn the_seam_class_is_a_built_in_class() {
    let edited = "IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK";
    let planted = |hgvsc: &str, extra: &str| {
        let extra = if hgvsc == "-" {
            extra.to_string()
        } else {
            format!("{extra};HGVSc={SEAM_MODEL}:{hgvsc}")
        };
        let row = default_row(
            "seam_ins",
            "1:157-158",
            "G",
            SEAM_MODEL,
            "missense_variant",
            &extra,
        );
        default_text(&[row.replacen(
            "\t-\t-\t-\t-\t-\t-\t",
            "\t59-61\t9-11\t3-4\tGP/GL\tggacca/ggGa\t-\t",
            1,
        )])
    };
    assert_eq!(classes::tables().field_classes_of(SEAM), Some(HGVSC_FIELDS));
    let (reference, candidate) = (
        planted("c.11_12insG", edited),
        planted("c.9_11delinsG", edited),
    );
    let (report, documented) = compare_default(&reference, &candidate, &[SEAM]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(SEAM));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &masked_once(SEAM)
    );
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let unedited = "IMPACT=MODERATE;STRAND=1";
    for (reference, candidate) in [
        (
            planted("c.11_12insG", unedited),
            planted("c.9_11delinsG", unedited),
        ),
        (planted("c.11_12insG", edited), planted("-", edited)),
        (
            planted("c.11_12insG", edited),
            planted("c.9_11delinsT", edited),
        ),
    ] {
        assert_eq!(
            default_under(SEAM, &reference, &candidate),
            ("red".to_string(), nothing())
        );
    }
}

const CLIPPED_MODEL: &str = "NM_990003.1";

/// What a row of the made-up edited model `NM_990003.1` prints beside its alleles and
/// notation: the strand (an empty text for a row without STRAND), the `BAM_EDIT` value (an
/// empty text for a row without one), the sequence name of its notation and the record's own
/// reference allele where the model's differs from it.
#[derive(Clone, Copy)]
struct ClippedRow<'a> {
    strand: &'a str,
    status: &'a str,
    name: &'a str,
    given_ref: Option<&'a str>,
}

impl Default for ClippedRow<'_> {
    fn default() -> Self {
        ClippedRow {
            strand: "1",
            status: "OK",
            name: CLIPPED_MODEL,
            given_ref: None,
        }
    }
}

/// The entry of one row of the made-up edited model: the allele `allele` against the
/// reference allele `used_ref`, both as the row prints them; `hgvsc` is the notation after
/// the sequence name.
fn clipped_entry(hgvsc: &str, used_ref: &str, allele: &str, row: ClippedRow) -> Entry {
    let mut held = vec!["IMPACT=MODIFIER".to_string()];
    if !row.strand.is_empty() {
        held.push(format!("STRAND={}", row.strand));
    }
    held.push(format!(
        "GIVEN_REF={};USED_REF={used_ref}",
        row.given_ref.unwrap_or(used_ref)
    ));
    if !row.status.is_empty() {
        held.push(format!("BAM_EDIT={}", row.status));
    }
    held.push(format!("HGVSc={}:{hgvsc}", row.name));
    default_entry(&default_row(
        "clipped",
        "1:160-161",
        allele,
        CLIPPED_MODEL,
        "intron_variant",
        &held.join(";"),
    ))
}

/// The rule on the two entries' own HGVSc values, a missing one handed over as `-`.
fn clipped_rule(reference: &Entry, candidate: Option<&Entry>) -> bool {
    edited_model_clipped_alleles(
        reference,
        candidate,
        hgvsc_of(reference),
        candidate.map_or("-", hgvsc_of),
        false,
    )
    .expect("the rule reads the row")
}

/// Whether the rule accepts the notation `candidate` against the notation `reference` on a
/// row of `allele` against `used_ref`; `row` changes the two sides' rows alike.
fn clipped_accepts(
    used_ref: &str,
    allele: &str,
    reference: &str,
    candidate: &str,
    row: ClippedRow,
) -> bool {
    clipped_rule(
        &clipped_entry(reference, used_ref, allele, row),
        Some(&clipped_entry(candidate, used_ref, allele, row)),
    )
}

/// The models, in the rows' order, whose reference row the rule takes against the
/// candidate's row of the same place.
fn taken_features(reference: &[Entry], candidate: &[Entry]) -> Vec<String> {
    assert_eq!(reference.len(), candidate.len());
    reference
        .iter()
        .zip(candidate)
        .filter(|(e, a)| clipped_rule(e, Some(a)))
        .map(|(e, _)| e.feature.clone())
        .collect()
}

const SECOND_ALTERNATE_MODELS: [&str; 5] = [
    "NM_001136154.1",
    "NM_001243428.1",
    "NM_001243432.2",
    "NM_001291391.1",
    "NM_004449.4",
];

/// The rows of a record whose two-base allele keeps the base it shares with the reference,
/// a second alternate sharing none: five reverse-strand models carrying `BAM_EDIT` `OK`
/// whose intron holds it, with `notation` after the coordinate, a model annotated
/// upstream and a model without `BAM_EDIT`.
fn second_alternate_text(notation: &str) -> String {
    let edited = "IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=TA;USED_REF=TA;BAM_EDIT=OK";
    let row = |model: &str, extra: &str, consequence: &str| {
        default_row_of(
            "21_39960821_TA/AA/T",
            "21:39960821-39960822",
            "AA",
            model,
            consequence,
            extra,
            "Transcript",
            "2078",
        )
    };
    let mut rows: Vec<String> = SECOND_ALTERNATE_MODELS
        .iter()
        .map(|model| {
            row(
                model,
                &format!("{edited};HGVSc={model}:c.-149-3952{notation}"),
                "intron_variant",
            )
        })
        .collect();
    rows.push(row(
        "NM_001243429.1",
        &edited.replace("STRAND=-1", "DISTANCE=3952;STRAND=-1"),
        "upstream_gene_variant",
    ));
    let unedited = edited.replace(";BAM_EDIT=OK", "");
    rows.push(row(
        "NR_111949.2",
        &format!("{unedited};HGVSc=NR_111949.2:n.127-3952A>T"),
        "intron_variant,non_coding_transcript_variant",
    ));
    default_text(&rows)
}

/// A record of two alternates keeps the base its two-base allele shares with the reference.
/// On the five reverse-strand models carrying `BAM_EDIT` `OK` whose intron holds it, Ensembl
/// VEP spells the row's USED_REF and Allele whole on the transcript's strand at the one
/// position the clip leaves, and the candidate spells the one base that changes: the rule
/// accepts each of the five rows. The record's row on the model without `BAM_EDIT`, where
/// both engines spell the clipped alleles, the edited model's row with no HGVSc and the two
/// sides reversed are not the rule's.
#[test]
fn the_clipped_alleles_rule_accepts_the_two_base_allele_of_a_record_with_a_second_alternate() {
    let reference = default_entries(&second_alternate_text("TA>TT"));
    let candidate = default_entries(&second_alternate_text("A>T"));
    assert_eq!(
        taken_features(&reference, &candidate),
        SECOND_ALTERNATE_MODELS
    );
    assert!(taken_features(&candidate, &reference).is_empty());
}

/// The comparison keeps the five rows of the record with a second alternate a difference,
/// with no class and under the class that renumbers an edited model's HGVSc.
#[test]
fn the_comparison_keeps_the_whole_alleles_of_a_second_alternate_a_difference() {
    let (whole, clipped) = (second_alternate_text("TA>TT"), second_alternate_text("A>T"));
    let (report, _) = compare_default(&whole, &clipped, &[]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "HGVSc", "differing"]), 5);
    assert_eq!(verdict_of(&whole, &clipped, &[HGVSC_OFFSET]), "red");
    assert_eq!(verdict_of(&whole, &whole, &[]), "green");
}

/// The clip takes one split of the bases a pair shares: the longest front, then the longest
/// back of what is left. Of a repeat it takes the front and leaves the rest inserted; the
/// same bases turned by one, which a base taken from each end would leave, are not the
/// clip's and stay a difference. Bases that repeat the kept base lose it at the front, and a
/// pair that shares a base at each end loses both, so no `delins` at one position is left
/// of it.
#[test]
fn the_clipped_alleles_rule_refuses_a_repeat_turned_by_one_base() {
    let plain = ClippedRow::default();
    assert!(!clipped_accepts(
        "AC",
        "ACACAC",
        "c.11_12insACACAC",
        "c.11_12insCACA",
        plain
    ));
    assert!(clipped_accepts(
        "AC",
        "ACACAC",
        "c.11_12insACACAC",
        "c.11_12insACAC",
        plain
    ));
    assert_eq!(shared_ends("AC", "ACACAC"), (2, 0));
    assert!(!clipped_accepts(
        "A",
        "ACCA",
        "c.11_12insACCA",
        "c.11_12insACC",
        plain
    ));
    assert!(clipped_accepts(
        "A",
        "ACCA",
        "c.11_12insACCA",
        "c.11_12insCCA",
        plain
    ));
    assert_eq!(shared_ends("A", "ACCA"), (1, 0));
    assert!(clipped_accepts(
        "TT",
        "TGT",
        "c.11_12insTGT",
        "c.11_12insG",
        plain
    ));
    assert_eq!(shared_ends("TT", "TGT"), (1, 1));
    for candidate in ["c.12delinsGT", "c.12delinsTG"] {
        assert!(
            !clipped_accepts("TT", "TGT", "c.12delinsTGT", candidate, plain),
            "{candidate}"
        );
    }
    // The front is taken before the back, and a base both alleles are is taken once.
    for (old, new, taken) in [
        ("TA", "TT", (1, 0)),
        ("T", "GGGT", (0, 1)),
        ("TACT", "TGGT", (1, 1)),
        ("GT", "TTT", (0, 1)),
        ("T", "T", (1, 0)),
        ("TT", "T", (1, 0)),
        ("T", "G", (0, 0)),
    ] {
        assert_eq!(shared_ends(old, new), taken, "{old} {new}");
    }
}

/// The rows of a record whose two-base allele keeps the base it shares with the reference,
/// on three models of one gene: `inserted` after each model's coordinates, under the
/// `BAM_EDIT` value `status`.
fn failed_edit_text(inserted: &str, status: &str) -> String {
    let extra = format!("IMPACT=HIGH;STRAND=1;GIVEN_REF=G;USED_REF=G;BAM_EDIT={status}");
    let rows: Vec<String> = [
        ("NM_015068.3", "c.956_957"),
        ("NM_001184961.1", "c.1058_1059"),
        ("NM_001172437.2", "c.1184_1185"),
    ]
    .iter()
    .map(|(model, coordinates)| {
        default_row(
            "7_94664512_eq_multi",
            "7:94664512",
            "GC",
            model,
            "frameshift_variant",
            &format!("{extra};HGVSc={model}:{coordinates}ins{inserted}"),
        )
    })
    .collect();
    default_text(&rows)
}

/// The rule is an applied edit's: a row carrying `BAM_EDIT` `FAILED` is refused. On three
/// models of one gene whose edit failed, Ensembl VEP spells the allele whole and the
/// candidate spells the one inserted base at the same coordinates. The same pair on a row
/// whose edit was applied is the rule's shape, and a row that carries `FAILED` on one side
/// alone is refused as well.
#[test]
fn the_clipped_alleles_rule_refuses_the_row_of_a_model_whose_edit_failed() {
    let entries =
        |inserted: &str, status: &str| default_entries(&failed_edit_text(inserted, status));
    let (whole, clipped) = (entries("GC", "FAILED"), entries("C", "FAILED"));
    assert!(taken_features(&whole, &clipped).is_empty());
    assert_eq!(
        taken_features(&entries("GC", "OK"), &entries("C", "OK")).len(),
        3
    );
    assert!(taken_features(&entries("GC", "OK"), &clipped).is_empty());
    assert!(taken_features(&whole, &entries("C", "OK")).is_empty());
}

/// The rows of the models whose edit failed stay a difference in the comparison, printed as
/// the engines wrote them.
#[test]
fn the_comparison_keeps_the_rows_of_a_model_whose_edit_failed_a_difference() {
    let (report, _) = compare_default(
        &failed_edit_text("GC", "FAILED"),
        &failed_edit_text("C", "FAILED"),
        &[],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "HGVSc", "differing"]), 3);
}

/// A candidate that drops a base the alleles do not share, or more or fewer bases than the
/// coordinates leave, a reference that does not spell the row's pair, other coordinates,
/// another operation or sequence name, a notation that prints no allele, a row without
/// `BAM_EDIT`, GIVEN_REF or USED_REF, a USED_REF the two sides do not share and a value that
/// is not the entry's own HGVSc stay differences.
#[test]
fn the_clipped_alleles_rule_keeps_another_value_a_difference() {
    let plain = ClippedRow::default();
    let reference = "c.11_12insTGGG";
    for candidate in [
        "c.11_12insGG",
        "c.11_12insTGG",
        "c.11_12insGGC",
        "c.11_12insTGGG",
        "c.12_13insGGG",
        "c.11_12delinsGGG",
        "c.11_12dup",
        "n.11_12insGGG",
        "-",
    ] {
        assert!(
            !clipped_accepts("T", "TGGG", reference, candidate, plain),
            "{candidate}"
        );
    }
    for (used_ref, allele, reference, candidate) in [
        // The coordinates leave none of an `ins` reference allele and one base of a lone
        // `delins`.
        ("TT", "TGT", "c.11_12insTGT", "c.11_12insGT"),
        ("TT", "TGT", "c.12delinsTGT", "c.12delinsG"),
        ("TACT", "TGGT", "c.12_13delinsTGGT", "c.12_13delinsGGT"),
        // A substitution names one position and one base a side.
        ("TTT", "TGT", "c.12_13TTT>TGT", "c.12_13T>G"),
        ("TTT", "TGT", "c.12TTT>TGT", "c.12TT>GT"),
        ("TTT", "TGT", "c.12TTT>TGT", "c.12T>C"),
        // The reference spells the row's own pair.
        ("T", "TGGG", "c.11_12insTGGC", "c.11_12insGGC"),
        ("TTT", "TGT", "c.12TTA>TGA", "c.12T>G"),
        ("TTT", "TGT", "c.12TT>TG", "c.12T>G"),
        // A deletion and a duplication print no allele, and an insertion's row has no
        // reference bases.
        ("TTT", "T", "c.12_13del", "c.12del"),
        ("-", "GGG", "c.11_12insGGG", "c.11_12insGG"),
    ] {
        assert!(
            !clipped_accepts(used_ref, allele, reference, candidate, plain),
            "{reference} {candidate}"
        );
    }
    // The row: no `BAM_EDIT` on either side or on one, no USED_REF, another on the other
    // side.
    let unedited = ClippedRow {
        status: "",
        ..plain
    };
    assert!(!clipped_accepts(
        "T",
        "TGGG",
        reference,
        "c.11_12insGGG",
        unedited
    ));
    let whole = clipped_entry(reference, "T", "TGGG", plain);
    let clipped = clipped_entry("c.11_12insGGG", "T", "TGGG", plain);
    assert!(clipped_rule(&whole, Some(&clipped)));
    assert!(!clipped_rule(
        &whole,
        Some(&clipped_entry("c.11_12insGGG", "T", "TGGG", unedited))
    ));
    assert!(!clipped_rule(
        &clipped_entry(reference, "T", "TGGG", unedited),
        Some(&clipped)
    ));
    assert!(!clipped_rule(
        &whole,
        Some(&clipped_entry("c.11_12insGGG", "A", "TGGG", plain))
    ));
    for missing in ["GIVEN_REF", "USED_REF"] {
        let mut without = clipped_entry(reference, "T", "TGGG", plain);
        without.fields.shift_remove(missing);
        assert!(!clipped_rule(&without, Some(&clipped)), "{missing}");
    }
    let other_name = ClippedRow {
        name: "NM_990004.1",
        ..plain
    };
    assert!(!clipped_rule(
        &whole,
        Some(&clipped_entry("c.11_12insGGG", "T", "TGGG", other_name))
    ));
    // The two notations under another key are not the rows' HGVSc, and a row the candidate
    // does not write has no notation to read.
    let rule = |a: Option<&Entry>, expected: &str, actual: &str| {
        edited_model_clipped_alleles(&whole, a, expected, actual, false)
            .expect("the rule reads the row")
    };
    assert!(rule(
        Some(&clipped),
        &whole.fields["HGVSc"],
        &clipped.fields["HGVSc"]
    ));
    assert!(!rule(
        Some(&clipped),
        &whole.fields["Extra"],
        &clipped.fields["Extra"]
    ));
    assert!(!rule(
        Some(&clipped),
        &whole.fields["HGVSc"],
        &whole.fields["HGVSc"]
    ));
    assert!(!rule(
        None,
        &whole.fields["HGVSc"],
        &clipped.fields["HGVSc"]
    ));
}

/// The tab format and the VCF format's CSQ carry the allele, the strand, GIVEN_REF, USED_REF
/// and `BAM_EDIT` under the default format's names, the JSON format as `variant_allele`,
/// `strand`, `given_ref`, `used_ref` and `bam_edit`: the rule accepts the same row in each
/// and keeps the same wrong value a difference.
#[test]
fn the_clipped_alleles_rule_reads_the_row_in_every_format() {
    let model = "NM_004449.4";
    let notations = ["c.-149-3952TA>TT", "c.-149-3952A>T", "c.-149-3952T>A"];
    // The rule takes the second notation against the first and refuses the third.
    let holds = |entry: &dyn Fn(&str) -> Entry| {
        let [reference, derived, wrong] = notations.map(entry);
        clipped_rule(&reference, Some(&derived)) && !clipped_rule(&reference, Some(&wrong))
    };
    let header = tab_header(&[
        "IMPACT",
        "STRAND",
        "GIVEN_REF",
        "USED_REF",
        "BAM_EDIT",
        "HGVSc",
    ]);
    let tab_row = format!(
        "21_39960821_TA/AA/T\t21:39960821-39960822\tAA\t2078\t{model}\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-1\tTA\tTA\tOK\t{model}:"
    );
    assert!(holds(&|notation| {
        tab_entries(&format!("{header}\n{tab_row}{notation}\n")).remove(0)
    }));
    let inputs = vec!["21\t39960821\t.\tTA\tAA,T\t.\t.\t.".to_string()];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|STRAND|GIVEN_REF|USED_REF|BAM_EDIT|HGVSc",
    );
    let csq_row = format!(
        "{}CSQ=AA|intron_variant|{model}|Transcript|-1|TA|TA|OK|{model}:",
        record_columns(&inputs[0])
    );
    assert!(holds(&|notation| {
        transcript_entry(
            vcf_entries(&format!("{header}{csq_row}{notation}\n"), &inputs),
            model,
        )
    }));
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("21_39960821_TA/AA/T")),
        ("allele_string", text("TA/AA/T")),
        ("most_severe_consequence", text("intron_variant")),
    ];
    // The entry of the JSON row with `notation`: `fields` set, the member `missing` left out.
    let json_entry = |notation: &str, fields: &[(&str, Value)], missing: Option<&str>| {
        let consequence = object(&[
            ("variant_allele", text("AA")),
            ("transcript_id", text(model)),
            (
                "consequence_terms",
                Value::List(vec![text("intron_variant")]),
            ),
            ("strand", Value::Int(-1)),
            ("bam_edit", text("OK")),
            ("given_ref", text("TA")),
            ("used_ref", text("TA")),
            ("hgvsc", text(&format!("{model}:{notation}"))),
        ]);
        let consequence = with_members(&consequence, fields);
        let consequence = match missing {
            Some(name) => without_member(&consequence, name),
            None => consequence,
        };
        transcript_entry(
            json_entries(&json_line(&record, consequence), &inputs),
            model,
        )
    };
    assert!(holds(&|notation| json_entry(notation, &[], None)));
    let reference = json_entry(notations[0], &[], None);
    for missing in ["given_ref", "used_ref", "bam_edit"] {
        assert!(
            !clipped_rule(
                &reference,
                Some(&json_entry("c.-149-3952A>T", &[], Some(missing)))
            ),
            "{missing}"
        );
    }
    assert!(!clipped_rule(
        &reference,
        Some(&json_entry(
            "c.-149-3952A>T",
            &[("bam_edit", text("FAILED"))],
            None
        ))
    ));
}

/// An HGVS_OFFSET on a coding BAM_EDIT OK row whose HGVSc names a duplication, which the
/// reference leaves empty and the candidate fills with a non-zero integer, is the class,
/// counted once per row in every format; a row without BAM_EDIT OK, a substitution's row, a
/// plain insertion's or a deletion's row (their notations move with the shift), a
/// duplication outside the coding sequence, a row without an HGVSc, a zero, a non-integer,
/// an offset the reference fills, a differing HGVSc beside the offset, and the field under
/// the HGVSc class alone stay differences.
#[test]
fn the_empty_offset_class_accepts_an_offset_the_reference_leaves_empty_on_an_edited_model_row() {
    let class = [HGVS_OFFSET];
    let extra = "IMPACT=MODERATE;STRAND=1;REFSEQ_OFFSET=-15;BAM_EDIT=OK;HGVSc=NM_198689.3:c.986_1006dup;HGVSp=NP_941962.1:p.Cys335_Gln336insArgAlaProThrSerSerCys";
    let positions = "\t1010-1011\t985-986\t329\tC/CAPTSSCR\ttgc/tGCGCCCCCACCTCCTCCTGCCgc\t-\t";
    let allele = "GCGCCCCCACCTCCTCCTGCC";
    let coding_row = |name: &str, location: &str, allele: &str, consequence: &str, fields: &str| {
        default_row(name, location, allele, "NM_198689.3", consequence, fields).replacen(
            "\t-\t-\t-\t-\t-\t-\t",
            positions,
            1,
        )
    };
    let reference = default_text(&[coding_row(
        "rs782070127",
        "21:46021521-46021522",
        allele,
        "inframe_insertion",
        extra,
    )]);
    // `text` with the offset `offset` after the Extra fields `fields`.
    let with_offset = |text: &str, fields: &str, offset: &str| {
        text.replace(fields, &format!("{fields};HGVS_OFFSET={offset}"))
    };
    let candidate = with_offset(&reference, extra, "21");
    let (report, documented) = compare_default(&reference, &candidate, &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(HGVS_OFFSET));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &masked_once(HGVS_OFFSET)
    );
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    assert_eq!(verdict_of(&reference, &candidate, &[HGVSC_OFFSET]), "red");
    let unedited_extra = extra.replace("BAM_EDIT=OK;", "");
    let unedited = reference.replace(extra, &unedited_extra);
    assert_eq!(
        verdict_of(
            &unedited,
            &with_offset(&unedited, &unedited_extra, "21"),
            &class
        ),
        "red"
    );
    for offset in ["0", "-0", "7.5", "dup"] {
        assert_eq!(
            verdict_of(&reference, &with_offset(&reference, extra, offset), &class),
            "red",
            "{offset}"
        );
    }
    // The shift is signed: a reverse-strand transcript's offset is negative on both engines.
    assert_eq!(
        verdict_of(&reference, &with_offset(&reference, extra, "-3"), &class),
        "green"
    );
    // Only a duplication keeps its notation whatever the shift; a substitution is never
    // shifted.
    for notation in [
        "c.985T>C",
        "c.985_986insGCG",
        "c.986_988del",
        "c.986_988delinsAT",
    ] {
        let other = reference.replace("c.986_1006dup", notation);
        let other_extra = extra.replace("c.986_1006dup", notation);
        assert_eq!(
            verdict_of(&other, &with_offset(&other, &other_extra, "3"), &class),
            "red",
            "{notation}"
        );
    }
    // The edited-model branch is entered only for an allele with CDS coordinates.
    let utr_extra = extra.replace("c.986_1006dup", "c.*12_*14dup");
    let utr = default_text(&[default_row(
        "rs_utr",
        "21:46021521-46021522",
        "GCG",
        "NM_198689.3",
        "3_prime_UTR_variant",
        &utr_extra,
    )]);
    assert_eq!(
        verdict_of(&utr, &with_offset(&utr, &utr_extra, "3"), &class),
        "red"
    );
    let no_hgvsc = "IMPACT=MODIFIER;STRAND=1;BAM_EDIT=OK";
    let intron = default_text(&[default_row(
        "rs_intron",
        "21:46020000-46020001",
        "-",
        "NM_198689.3",
        "intron_variant",
        no_hgvsc,
    )]);
    assert_eq!(
        verdict_of(&intron, &with_offset(&intron, no_hgvsc, "4"), &class),
        "red"
    );
    assert_eq!(
        verdict_of(
            &with_offset(&reference, extra, "3"),
            &with_offset(&reference, extra, "21"),
            &class
        ),
        "red"
    );
    assert_eq!(
        verdict_of(&with_offset(&reference, extra, "21"), &reference, &class),
        "red"
    );
    // A differing HGVSc beside the offset is a difference of its own: red under this class
    // alone, and green only when the HGVSc class accepts it too, each class credited once.
    let renumbered = candidate.replace("c.986_1006dup", "c.1001_1021dup");
    assert_eq!(verdict_of(&reference, &renumbered, &class), "red");
    let (report, documented) =
        compare_default(&reference, &renumbered, &[HGVS_OFFSET, HGVSC_OFFSET]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(HGVS_OFFSET, 1), (HGVSC_OFFSET, 1)])
    );
    // The tab, VCF and JSON formats carry the offset as a field of its own and count the
    // row once.
    let header = tab_header(&[
        "IMPACT",
        "STRAND",
        "REFSEQ_OFFSET",
        "BAM_EDIT",
        "HGVSc",
        "HGVSp",
        "HGVS_OFFSET",
    ]);
    let tab_row = format!(
        "rs782070127\t21:46021521-46021522\t{allele}\tENSG1\tNM_198689.3\tTranscript\tinframe_insertion{positions}MODERATE\t1\t-15\tOK\tNM_198689.3:c.986_1006dup\tNP_941962.1:p.Cys335_Gln336insArgAlaProThrSerSerCys\t"
    );
    let tab = |offset: &str| tab_entries(&format!("{header}\n{tab_row}{offset}\n"));
    assert_eq!(
        under_class(HGVS_OFFSET, &tab("-"), &tab("21")),
        ("green".to_string(), once(HGVS_OFFSET))
    );
    let inputs = vec![format!("21\t46021521\trs782070127\tT\tT{allele}\t.\t.\t.")];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|CDS_position|BAM_EDIT|HGVSc|HGVS_OFFSET",
    );
    let csq_row = format!(
        "{}CSQ={allele}|inframe_insertion|NM_198689.3|Transcript|985-986|OK|NM_198689.3:c.986_1006dup|",
        record_columns(&inputs[0])
    );
    let vcf = |offset: &str| vcf_entries(&format!("{header}{csq_row}{offset}\n"), &inputs);
    assert_eq!(
        under_class(HGVS_OFFSET, &vcf(""), &vcf("21")),
        ("green".to_string(), once(HGVS_OFFSET))
    );
    assert_eq!(under_class(HGVS_OFFSET, &vcf(""), &vcf("0")).0, "red");
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("rs782070127")),
        ("most_severe_consequence", text("inframe_insertion")),
    ];
    let consequence = object(&[
        ("variant_allele", text(allele)),
        ("transcript_id", text("NM_198689.3")),
        (
            "consequence_terms",
            Value::List(vec![text("inframe_insertion")]),
        ),
        ("cds_start", Value::Int(986)),
        ("cds_end", Value::Int(985)),
        ("bam_edit", text("OK")),
        ("hgvsc", text("NM_198689.3:c.986_1006dup")),
    ]);
    let object_of = |fields: &[(&str, Value)]| {
        json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        )
    };
    assert_eq!(
        under_class(
            HGVS_OFFSET,
            &object_of(&[]),
            &object_of(&[("hgvs_offset", Value::Int(21))])
        ),
        ("green".to_string(), once(HGVS_OFFSET))
    );
    assert_eq!(
        under_class(
            HGVS_OFFSET,
            &object_of(&[]),
            &object_of(&[("hgvs_offset", Value::Int(0))])
        )
        .0,
        "red"
    );
}

/// A deletion's row both engines give a non-zero HGVS_OFFSET, whose two `del` notations
/// name one transcript and numbering and end at the same coordinate while the reference's
/// starts later, is the class, counted once per row in every format, with or without a
/// SHIFT_LENGTH field and across the 5' UTR, coding, intronic and 3' UTR numberings; under
/// --shift_length an n. candidate must equal the row's cDNA position and a same-region c.
/// span must be as long as it; another end, a later candidate start, a duplication or
/// substitution on either side, another transcript, an empty notation, a missing or zero
/// offset, an insertion's row and a zero SHIFT_LENGTH stay differences.
#[test]
fn the_clamped_span_class_accepts_the_longer_span_ending_with_the_clamped_one() {
    let class = [CLAMPED_SPAN];
    let extra = "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000785034.1:n.34_43del;HGVS_OFFSET=33;SHIFT_LENGTH=33";
    let row = default_row(
        "21:13037980:TTGCCTAACG:T",
        "21:13038013-13038033",
        "-",
        "ENST00000785034",
        "non_coding_transcript_exon_variant",
        extra,
    );
    let reference =
        default_text(&[row.replacen("\t-\t-\t-\t-\t-\t-\t", "\t23-43\t-\t-\t-\t-\t-\t", 1)]);
    let written = "n.34_43del";
    let candidate = |notation: &str| reference.replace(written, notation);

    let (report, documented) = compare_default(&reference, &candidate("n.23_43del"), &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(CLAMPED_SPAN));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &masked_once(CLAMPED_SPAN)
    );
    assert_eq!(swapped(&reference, written, "n.23_43del", &[]), "red");
    // A one-base reference span against the candidate's two, and the candidate clamped to
    // n.1 on a row whose cDNA position is still open on the first base.
    for (reference_notation, candidate_notation, cdna) in [
        ("n.14del", "n.13_14del", "13-14"),
        ("n.29_31del", "n.1_31del", "?-31"),
        ("n.2_4del", "n.1_4del", "1-4"),
    ] {
        let planted = reference
            .replace(written, reference_notation)
            .replace("\t23-43\t", &format!("\t{cdna}\t"));
        assert_eq!(
            swapped(&planted, reference_notation, candidate_notation, &class),
            "green",
            "{reference_notation}"
        );
    }
    // In n. numbering the candidate's coordinates are the row's cDNA position under
    // --shift_length; a span off by one on either side stays a difference, as does a
    // position the row reads differently. In c. numbering a same-region span must be as
    // long as the row's cDNA span.
    for (notation, cdna) in [
        ("n.22_43del", "23-43"),
        ("n.23_43del", "24-43"),
        ("n.23_43del", "23-44"),
        ("n.1_31del", "2-31"),
        ("n.22_43del", "23-43/1500"),
    ] {
        let planted = reference.replace("\t23-43\t", &format!("\t{cdna}\t"));
        assert_eq!(
            swapped(&planted, written, notation, &class),
            "red",
            "{notation} {cdna}"
        );
    }
    let total_length = reference.replace("\t23-43\t", "\t23-43/1500\t");
    assert_eq!(
        swapped(&total_length, written, "n.23_43del", &class),
        "green"
    );
    let coding = reference
        .replace(written, "c.-355_-327del")
        .replace("\t23-43\t", "\t7-41\t");
    assert_eq!(
        swapped(&coding, "c.-355_-327del", "c.-361_-327del", &class),
        "green"
    );
    for notation in ["c.-362_-327del", "c.-360_-327del"] {
        assert_eq!(
            swapped(&coding, "c.-355_-327del", notation, &class),
            "red",
            "{notation}"
        );
    }
    // A c. span across the c.-1/c.1 boundary or into an intron (the distance before the stop
    // codon's last base included), and a row with an open side, say nothing.
    for (reference_notation, candidate_notation, cdna) in [
        ("c.-2_3del", "c.-6_3del", "7-41"),
        ("c.-301_-228+19del", "c.-307_-228+19del", "7-41"),
        ("c.*-3_*2del", "c.*-5_*2del", "7-41"),
        ("c.*12_*20del", "c.*8_*20del", "?-41"),
    ] {
        let planted = reference
            .replace(written, reference_notation)
            .replace("\t23-43\t", &format!("\t{cdna}\t"));
        assert_eq!(
            swapped(&planted, reference_notation, candidate_notation, &class),
            "green",
            "{reference_notation}"
        );
    }
    // The 5' UTR numbering counts down towards c.-1, an intron offset follows its exon base,
    // the 3' UTR counts up from c.*1 (no SHIFT_LENGTH on these rows, so the cDNA position is
    // the written span's and says nothing).
    for (reference_notation, candidate_notation) in [
        ("c.-355_-327del", "c.-361_-327del"),
        ("c.-301_-228+19del", "c.-307_-228+19del"),
        ("c.-2_3del", "c.-6_3del"),
        ("c.*12_*20del", "c.*8_*20del"),
        ("c.100+2_100+9del", "c.100+1_100+9del"),
    ] {
        let planted = reference
            .replace(written, reference_notation)
            .replace(";SHIFT_LENGTH=33", "");
        assert_eq!(
            swapped(&planted, reference_notation, candidate_notation, &class),
            "green",
            "{reference_notation}"
        );
    }
    // Under --hgvs alone the cDNA position is the written span's, which says nothing of the
    // candidate.
    let unlabelled = reference
        .replace(";SHIFT_LENGTH=33", "")
        .replace("\t23-43\t", "\t?-10\t");
    assert_eq!(swapped(&unlabelled, written, "n.23_43del", &class), "green");
    for notation in [
        "n.23_42del",
        "n.35_43del",
        "n.23_43dup",
        "n.43del",
        "c.23_43del",
        "-",
    ] {
        assert_eq!(
            swapped(&reference, written, notation, &class),
            "red",
            "{notation}"
        );
    }
    assert_eq!(
        verdict_of(&candidate("n.23_43del"), &reference, &class),
        "red"
    );
    assert_eq!(
        verdict_of(
            &reference,
            &candidate("n.23_43del").replace("HGVSc=ENST00000785034.1", "HGVSc=ENST00000785035.1"),
            &class
        ),
        "red"
    );
    let substitution_reference = reference.replace(written, "n.43T>A");
    assert_eq!(
        swapped(&substitution_reference, "n.43T>A", "n.23_43del", &class),
        "red"
    );
    for (variant, replacement) in [
        ("HGVS_OFFSET=33;", ""),
        ("HGVS_OFFSET=33;", "HGVS_OFFSET=0;"),
        ("SHIFT_LENGTH=33", "SHIFT_LENGTH=0"),
    ] {
        let planted = reference.replace(variant, replacement);
        assert_eq!(
            swapped(&planted, written, "n.23_43del", &class),
            "red",
            "{variant} {replacement}"
        );
    }
    let insertion = reference
        .replace("\t-\t", "\tT\t")
        .replace("21:13038013-13038033", "21:13038013-13038014");
    assert_eq!(swapped(&insertion, written, "n.23_43del", &class), "red");
    // The tab, VCF and JSON formats count the row once.
    let header = tab_header(&["IMPACT", "STRAND", "HGVSc", "HGVS_OFFSET", "SHIFT_LENGTH"]);
    let tab_row = "21:13037980:TTGCCTAACG:T\t21:13038013-13038033\t-\tENSG1\tENST00000785034\tTranscript\tnon_coding_transcript_exon_variant\t23-43\t-\t-\t-\t-\t-\tMODIFIER\t1\tENST00000785034.1:";
    let tab = |notation: &str| tab_entries(&format!("{header}\n{tab_row}{notation}\t33\t33\n"));
    assert_eq!(
        under_class(CLAMPED_SPAN, &tab("n.34_43del"), &tab("n.23_43del")),
        ("green".to_string(), once(CLAMPED_SPAN))
    );
    let inputs = vec![
        "21\t13037979\t21:13037980:TTGCCTAACG:T\tTTGCCTAACGGTTCTGCTCAGC\tT\t.\t.\t.".to_string(),
    ];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|HGVSc|HGVS_OFFSET|SHIFT_LENGTH",
    );
    let csq_row = format!(
        "{}CSQ=-|non_coding_transcript_exon_variant|ENST00000785034|Transcript|23-43|ENST00000785034.1:",
        record_columns(&inputs[0])
    );
    let vcf =
        |notation: &str| vcf_entries(&format!("{header}{csq_row}{notation}|33|33\n"), &inputs);
    assert_eq!(
        under_class(CLAMPED_SPAN, &vcf("n.34_43del"), &vcf("n.23_43del")),
        ("green".to_string(), once(CLAMPED_SPAN))
    );
    // The text `json.dumps` writes for the record's object with the notation `hgvsc`.
    let object_text = |hgvsc: &str| {
        let record = [
            ("input", text(&inputs[0])),
            ("id", text("21:13037980:TTGCCTAACG:T")),
            ("allele_string", text("TGCCTAACGGTTCTGCTCAGC/-")),
            (
                "most_severe_consequence",
                text("non_coding_transcript_exon_variant"),
            ),
        ];
        let consequence = object(&[
            ("variant_allele", text("-")),
            ("transcript_id", text("ENST00000785034")),
            (
                "consequence_terms",
                Value::List(vec![text("non_coding_transcript_exon_variant")]),
            ),
            ("cdna_start", Value::Int(23)),
            ("cdna_end", Value::Int(43)),
            ("hgvsc", text(hgvsc)),
            ("hgvs_offset", Value::Int(33)),
            ("shift_length", Value::Int(33)),
        ]);
        json_line(&record, consequence)
    };
    let read = |expected: &str, actual: &str| {
        under_class(
            CLAMPED_SPAN,
            &json_entries(expected, &inputs),
            &json_entries(actual, &inputs),
        )
    };
    let clamped = object_text("ENST00000785034.1:n.34_43del");
    assert_eq!(
        read(&clamped, &object_text("ENST00000785034.1:n.23_43del")),
        ("green".to_string(), once(CLAMPED_SPAN))
    );
    assert_eq!(
        read(&clamped, &object_text("ENST00000785034.1:n.35_43del")).0,
        "red"
    );
    // The JSON object writes the two cDNA sides apart and leaves an open one out; the end
    // still binds the candidate's, as the text `?-43` does.
    let open_start = clamped.replace("\"cdna_start\": 23, ", "");
    assert_ne!(open_start, clamped);
    assert_eq!(
        read(&open_start, &open_start.replace(written, "n.1_43del")).0,
        "green"
    );
    assert_eq!(
        read(
            &open_start.replace(written, "n.34_42del"),
            &open_start.replace(written, "n.1_42del")
        )
        .0,
        "red"
    );
    let moved_end = clamped.replace("\"cdna_end\": 43", "\"cdna_end\": 42");
    assert_ne!(moved_end, clamped);
    assert_eq!(
        read(&moved_end, &moved_end.replace(written, "n.23_43del")).0,
        "red"
    );
}

/// One default-format row of a deletion under an identifier: IMPACT and STRAND, then
/// `pairs` in the Extra text.
fn offset_row(location: &str, feature: &str, consequence: &str, pairs: &str) -> String {
    default_text(&[default_row(
        "shifted_deletion",
        location,
        "-",
        feature,
        consequence,
        &format!("IMPACT=MODIFIER;STRAND=1{pairs}"),
    )])
}

/// The Extra fields of a notation, its offset and the plugin's two fields.
fn intron_offset_pairs(hgvsc: &str, offset: i64, start: i64, end: i64) -> String {
    format!(
        ";HGVSc={hgvsc};HGVS_OFFSET={offset};HGVS_IntronStartOffset={start};HGVS_IntronEndOffset={end}"
    )
}

/// A deletion straddling the first base of a transcript and shifted three bases: Ensembl
/// VEP's notation starts one intron base later than the candidate's, so the HGVSIntronOffset
/// plugin's start offset reads 2 against 1 while the end offset agrees. The start offset is
/// the class's with the notation, the row counted once, as is 9 against the 0 of an exonic
/// start. A candidate offset that is not its notation's start coordinate's, a differing end
/// offset, the same two values beside one notation on both sides, on a row the reference
/// gives no notation and where the candidate lacks the field stay differences.
#[test]
fn the_clamped_span_class_accepts_the_start_offset_of_each_notation() {
    let row = |pairs: &str| {
        offset_row(
            "18:2000202-2000212",
            "ENST00000000035",
            "splice_donor_variant,coding_sequence_variant,intron_variant",
            pairs,
        )
    };
    let notation = |coordinates: &str| format!("ENST00000000035.1:c.{coordinates}");
    let reference = row(&intron_offset_pairs(&notation("2+2_2+11del"), 3, 2, 11));
    let candidate = row(&intron_offset_pairs(&notation("2+1_2+11del"), 3, 1, 11));
    let (report, documented) = compare_default(&reference, &candidate, &[CLAMPED_SPAN]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(CLAMPED_SPAN));
    for column in ["HGVS_IntronStartOffset", "HGVSc"] {
        assert_eq!(
            at(&report, &["columns", column, "masked"]),
            &masked_once(CLAMPED_SPAN),
            "{column}"
        );
    }
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let exonic_reference = row(&intron_offset_pairs(
        "ENST00000000035.1:n.29+9_30-944del",
        37,
        9,
        -944,
    ));
    let exonic_candidate = row(&intron_offset_pairs(
        "ENST00000000035.1:n.1_30-944del",
        37,
        0,
        -944,
    ));
    assert_eq!(
        default_under(CLAMPED_SPAN, &exonic_reference, &exonic_candidate),
        ("green".to_string(), once(CLAMPED_SPAN))
    );
    for (coordinates, start, end, applied) in [
        ("2+1_2+11del", 0, 11, once(CLAMPED_SPAN)),
        ("2+1_2+11del", 1, 10, once(CLAMPED_SPAN)),
        ("2+2_2+11del", 1, 11, nothing()),
    ] {
        let other = row(&intron_offset_pairs(&notation(coordinates), 3, start, end));
        assert_eq!(
            default_under(CLAMPED_SPAN, &reference, &other),
            ("red".to_string(), applied),
            "{coordinates} {start} {end}"
        );
    }
    let missing =
        row(";HGVSc=ENST00000000035.1:c.2+1_2+11del;HGVS_OFFSET=3;HGVS_IntronEndOffset=11");
    let (report, _) = compare_default(&reference, &missing, &[CLAMPED_SPAN]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        number(&report, &["columns", "HGVS_IntronStartOffset", "differing"]),
        1
    );
    let (e, a) = (
        default_entries(&reference).remove(0),
        default_entries(&candidate).remove(0),
    );
    let start_offset = "HGVS_IntronStartOffset";
    assert!(class_accepts(
        CLAMPED_SPAN,
        &e,
        start_offset,
        "2",
        "1",
        Some(&a)
    ));
    assert!(!class_accepts(
        CLAMPED_SPAN,
        &e,
        start_offset,
        "2",
        "1",
        None
    ));
    assert!(!class_accepts(
        CLAMPED_SPAN,
        &e,
        "HGVS_IntronEndOffset",
        "11",
        "10",
        Some(&a)
    ));
    let unwritten = default_entries(&row("")).remove(0);
    assert!(!class_accepts(
        CLAMPED_SPAN,
        &unwritten,
        start_offset,
        ABSENT,
        "1",
        Some(&a)
    ));
}

/// Through the comparison itself the class sets aside the HGVSc of a row of a model whose
/// edit the alignment applied, where Ensembl VEP spells the row's whole alleles and the
/// candidate the alleles the clip leaves, counted once. The same pair on a row of a model
/// whose edit failed, a notation the candidate leaves out and a notation of another allele
/// stay differences.
#[test]
fn the_unminimised_allele_class_is_a_built_in_class() {
    assert_eq!(
        classes::REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC,
        UNMINIMISED
    );
    assert_eq!(
        classes::tables().field_classes_of(UNMINIMISED),
        Some(HGVSC_FIELDS)
    );
    let (whole, clipped, wrong) = ("TA>TT", "A>T", "T>A");
    let planted = |notation: &str, status: &str| {
        let mut extra =
            format!("IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=TA;USED_REF=TA;BAM_EDIT={status}");
        if !notation.is_empty() {
            extra.push_str(&format!(";HGVSc={CLIPPED_MODEL}:c.-149-3952{notation}"));
        }
        default_text(&[default_row(
            "clipped",
            "21:39960821-39960822",
            "AA",
            CLIPPED_MODEL,
            "intron_variant",
            &extra,
        )])
    };
    let (report, documented) = compare_default(
        &planted(whole, "OK"),
        &planted(clipped, "OK"),
        &[UNMINIMISED],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(UNMINIMISED));
    assert_eq!(
        at(&report, &["columns", "HGVSc", "masked"]),
        &masked_once(UNMINIMISED)
    );
    assert_eq!(
        verdict_of(&planted(whole, "OK"), &planted(clipped, "OK"), &[]),
        "red"
    );
    for (reference, candidate) in [
        (planted(whole, "FAILED"), planted(clipped, "FAILED")),
        (planted(whole, "OK"), planted("", "OK")),
        (planted(whole, "OK"), planted(wrong, "OK")),
    ] {
        assert_eq!(
            default_under(UNMINIMISED, &reference, &candidate),
            ("red".to_string(), nothing())
        );
    }
}
