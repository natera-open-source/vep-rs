// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Built-in classes read through the comparison of two entry lists, on planted rows.

mod compare_common;

use compare_common::*;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::consts::HGVSC_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::fields::entry_kind;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::record_rows::RecordRows;
use vep_formats_compare::rules::edited_model::edited_model_clipped_alleles;
use vep_formats_compare::rules::edited_model::edited_model_seam_insertion_replacement;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::nearest_exon::nearest_exon_items;
use vep_formats_compare::rules::record_alleles::first_alternate;
use vep_formats_compare::rules::record_alleles::record_allele_lengths;
use vep_formats_compare::rules::record_alleles::vcf_alleles_as_spelt;
use vep_formats_compare::rules::record_alleles::vcf_first_alternate;
use vep_formats_compare::rules::utrannotator::cdna_span;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

/// The one entry of a default-format output of one row.
fn one_entry(row: &str) -> Entry {
    let mut entries = parse_default(&default_text(&[row])).expect("the row parses");
    assert_eq!(entries.len(), 1);
    entries.remove(0)
}

/// The entry's own HGVSc value, `-` for an entry without one: the form a blank value is
/// handed to a rule in.
fn hgvsc_of(entry: &Entry) -> &str {
    HGVSC_FIELDS
        .iter()
        .find_map(|f| entry.fields.get(*f))
        .map_or("-", String::as_str)
}

const SEAM_MODEL: &str = "NM_990001.2";

/// What a row of the made-up model `NM_990001.2` prints beside its allele and notation: an
/// insertion between two genomic bases, by default the two around the coding bases 9 to 11
/// (three bases of the codons of coding bases 7 to 12), which the model holds and the
/// genome lacks.
#[derive(Clone, Copy)]
struct SeamRow<'a> {
    cds: &'a str,
    protein: &'a str,
    extra: &'a str,
}

impl Default for SeamRow<'_> {
    fn default() -> Self {
        SeamRow {
            cds: "9-11",
            protein: "3-4",
            extra: "IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK",
        }
    }
}

/// The entry of one row of the made-up model: `hgvsc` is the notation after the sequence
/// name, `-` for a row without one.
fn seam_entry(allele: &str, hgvsc: &str, codons: &str, row: SeamRow) -> Entry {
    let extra = if hgvsc == "-" {
        row.extra.to_string()
    } else {
        format!("{};HGVSc={SEAM_MODEL}:{hgvsc}", row.extra)
    };
    let text = default_row(
        "seam_ins",
        "1:157-158",
        allele,
        SEAM_MODEL,
        "missense_variant",
        &extra,
    );
    let text = text.replacen(
        "\t-\t-\t-\t-\t-\t-\t",
        &format!(
            "\t59-61\t{}\t{}\tGP/GL\t{codons}\t-\t",
            row.cds, row.protein
        ),
        1,
    );
    one_entry(&text)
}

/// Whether the rule accepts the notation `candidate` on a row of `allele` whose reference
/// notation is `reference`.
fn seam_accepts(
    allele: &str,
    codons: &str,
    candidate: &str,
    reference: &str,
    row: SeamRow,
) -> bool {
    let e = seam_entry(allele, reference, codons, row);
    let a = seam_entry(allele, candidate, codons, row);
    edited_model_seam_insertion_replacement(&e, Some(&a), hgvsc_of(&e), hgvsc_of(&a))
}

/// On a reverse-strand model the reference's `ins` spells the row's Allele
/// reverse-complemented and the candidate replaces the model's bases with that reading: an
/// allele of one base three times on three of the same base at coding bases 2847 to 2849 is
/// an inversion. The allele as written there, and its reverse complement on a forward
/// strand, stay differences; a row without STRAND takes the reading the reference spells.
#[test]
fn the_seam_class_reads_the_allele_on_the_transcripts_strand() {
    let reverse = SeamRow {
        cds: "2847-2849",
        protein: "949-950",
        extra: "IMPACT=MODERATE;STRAND=-1;BAM_EDIT=OK",
    };
    let between = |bases: &str| format!("c.2847+2_2847+3ins{bases}");
    let accepts = |allele: &str, codons: &str, candidate: &str, inserted: &str, row: SeamRow| {
        seam_accepts(allele, codons, candidate, &between(inserted), row)
    };
    assert!(accepts(
        "AAA",
        "gaaaag/gaTTTg",
        "c.2847_2849inv",
        "TTT",
        reverse
    ));
    assert!(accepts(
        "C",
        "gaaaag/gaGg",
        "c.2847_2849delinsG",
        "G",
        reverse
    ));
    assert!(accepts("T", "gaaaag/gaAg", "c.2848_2849del", "A", reverse));
    assert!(accepts("CTT", "gaaaag/gaAAGg", "c.2849A>G", "AAG", reverse));
    assert!(accepts("TTT", "gaaaag/gaAAAg", "-", "AAA", reverse));
    assert!(!accepts("AAA", "gaaaag/gaAAAg", "-", "AAA", reverse));
    assert!(!accepts(
        "AAA",
        "gaaaag/gaTTTg",
        "c.2847_2849inv",
        "AAA",
        reverse
    ));
    let forward = SeamRow {
        extra: "IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK",
        ..reverse
    };
    assert!(accepts(
        "TTT",
        "gaaaag/gaTTTg",
        "c.2847_2849inv",
        "TTT",
        forward
    ));
    assert!(!accepts(
        "AAA",
        "gaaaag/gaTTTg",
        "c.2847_2849inv",
        "TTT",
        forward
    ));
    let unstranded = SeamRow {
        extra: "IMPACT=MODERATE;BAM_EDIT=OK",
        ..reverse
    };
    assert!(accepts(
        "AAA",
        "gaaaag/gaTTTg",
        "c.2847_2849inv",
        "TTT",
        unstranded
    ));
    assert!(accepts("AAA", "gaaaag/gaAAAg", "-", "AAA", unstranded));
    assert!(!accepts(
        "AAA",
        "gaaaag/gaAAAg",
        "c.2847_2849inv",
        "TTT",
        unstranded
    ));
}

/// A span of one base prints one CDS position, which no ordinary insertion prints. A span
/// of two bases prints the two positions an ordinary insertion between them prints, and the
/// alternate codons tell the two rows apart: on a span they hold the allele in the bases'
/// place, on an ordinary insertion they are longer than the reference codons by the whole
/// allele.
#[test]
fn the_seam_class_accepts_spans_of_one_and_two_bases() {
    let one = SeamRow {
        cds: "1059",
        protein: "353",
        ..SeamRow::default()
    };
    let one_base_forms = [
        ("T", "c.1059G>T"),
        ("C", "c.1059G>C"),
        ("GG", "c.1059dup"),
        ("TA", "c.1059delinsTA"),
        ("GA", "c.1059_1060insA"),
        ("AG", "c.1058_1059insA"),
        ("G", "-"),
    ];
    for (allele, notation) in one_base_forms {
        assert!(
            seam_accepts(
                allele,
                &format!("ctg/ct{allele}"),
                notation,
                &format!("c.1060_1061ins{allele}"),
                one
            ),
            "{allele}"
        );
    }
    // An allele of the base three times is an insertion of two more copies after it; a repeat
    // count is read with two positions.
    assert!(seam_accepts(
        "GGG",
        "ctg/ctGGG",
        "c.1059_1060insGG",
        "c.1060_1061insGGG",
        one
    ));
    assert!(!seam_accepts(
        "GGG",
        "ctg/ctGGG",
        "c.1059[3]",
        "c.1060_1061insGGG",
        one
    ));
    let two = SeamRow {
        cds: "784-785",
        protein: "262",
        ..SeamRow::default()
    };
    let two_base_forms = [
        ("TT", "c.784_785delinsTT"),
        ("GT", "c.784_785inv"),
        ("T", "c.784_785delinsT"),
        ("A", "c.785del"),
        ("AT", "c.785C>T"),
        ("ACAC", "c.784_785dup"),
        ("ACACAC", "c.784_785[3]"),
        ("TAC", "c.783_784insT"),
        ("AC", "-"),
    ];
    for (allele, notation) in two_base_forms {
        assert!(
            seam_accepts(
                allele,
                &format!("acg/{allele}g"),
                notation,
                &format!("c.786_787ins{allele}"),
                two
            ),
            "{allele}"
        );
    }
    // An ordinary insertion between coding bases 784 and 785 prints the same CDS position.
    for notation in ["c.784_785delinsTT", "c.784_785insTT", "-"] {
        assert!(
            !seam_accepts("TT", "acg/aTTcg", notation, "c.786_787insTT", two),
            "{notation}"
        );
    }
    // Under --total_length each position carries its sequence's length after a slash.
    let lengths = SeamRow {
        cds: "784-785/1500",
        protein: "262/499",
        ..SeamRow::default()
    };
    assert!(seam_accepts(
        "TT",
        "acg/TTg",
        "c.784_785delinsTT",
        "c.786_787insTT",
        lengths
    ));
}

const CLIPPED_MODEL: &str = "NM_990003.1";

/// What a row of the made-up edited model `NM_990003.1` prints beside its alleles and
/// notation: the strand (an empty text for a row without STRAND), the `BAM_EDIT` value (an
/// empty text for a row without one) and the record's own reference allele where the
/// model's differs from it.
#[derive(Clone, Copy)]
struct ClippedRow<'a> {
    strand: &'a str,
    status: &'a str,
    given_ref: Option<&'a str>,
}

impl Default for ClippedRow<'_> {
    fn default() -> Self {
        ClippedRow {
            strand: "1",
            status: "OK",
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
    held.push(format!("HGVSc={CLIPPED_MODEL}:{hgvsc}"));
    one_entry(&default_row(
        "clipped",
        "1:160-161",
        allele,
        CLIPPED_MODEL,
        "intron_variant",
        &held.join(";"),
    ))
}

/// Whether the rule accepts the notation `candidate` against the notation `reference` on a
/// row of `allele` against `used_ref`.
fn clipped_accepts(
    used_ref: &str,
    allele: &str,
    reference: &str,
    candidate: &str,
    renumbered: bool,
    row: ClippedRow,
) -> bool {
    let e = clipped_entry(reference, used_ref, allele, row);
    let a = clipped_entry(candidate, used_ref, allele, row);
    edited_model_clipped_alleles(&e, Some(&a), hgvsc_of(&e), hgvsc_of(&a), renumbered)
        .expect("the rule reads the row")
}

/// Against the reference's whole alleles the candidate's are accepted when they are the
/// row's pair without the bases it shares at an end, as many as the coordinates leave: the
/// inserted bases of an `ins` whose reference allele is all shared, the bases that replace
/// the rest of a `delins` span, the one base a substitution changes. The pair is read on the
/// transcript's strand, a row without STRAND admitting either reading, on a row whose
/// `BAM_EDIT` is `OK`, in either numbering and at intron or untranslated coordinates.
#[test]
fn the_clipped_alleles_class_accepts_each_notation_the_clip_leaves() {
    let plain = ClippedRow::default();
    let forms = [
        ("T", "TGGG", "c.11_12insTGGG", "c.11_12insGGG"),
        ("T", "GGGT", "c.10_11insGGGT", "c.10_11insGGG"),
        ("TT", "TGT", "c.11_12insTGT", "c.11_12insG"),
        ("TT", "TGGG", "c.12delinsTGGG", "c.12delinsGGG"),
        ("TAC", "TGG", "c.12_13delinsTGG", "c.12_13delinsGG"),
        ("TACT", "TGGT", "c.12_13delinsTGGT", "c.12_13delinsGG"),
        (
            "TAC",
            "TGG",
            "c.100+5_101-3delinsTGG",
            "c.100+5_101-3delinsGG",
        ),
        ("TTT", "TGT", "c.12TTT>TGT", "c.12T>G"),
        ("TA", "AA", "c.11TA>AA", "c.11T>A"),
    ];
    let failed = ClippedRow {
        status: "FAILED",
        ..plain
    };
    for (used_ref, allele, reference, candidate) in forms {
        let case = format!("{reference} {candidate}");
        assert!(
            clipped_accepts(used_ref, allele, reference, candidate, false, plain),
            "{case}"
        );
        assert!(
            !clipped_accepts(used_ref, allele, candidate, reference, false, plain),
            "{case}"
        );
        assert!(
            !clipped_accepts(used_ref, allele, reference, candidate, false, failed),
            "{case}"
        );
        assert!(
            clipped_accepts(
                used_ref,
                allele,
                &reference.replace("c.", "n."),
                &candidate.replace("c.", "n."),
                false,
                plain
            ),
            "{case}"
        );
    }
    for coordinates in ["c.636+18_636+19", "c.-3_-2", "c.*12_*13"] {
        assert!(
            clipped_accepts(
                "T",
                "TGGG",
                &format!("{coordinates}insTGGG"),
                &format!("{coordinates}insGGG"),
                false,
                plain
            ),
            "{coordinates}"
        );
    }
    // On a reverse-strand model the pair is reverse-complemented and then shares its first
    // base. The pair as written is not the shape there, nor is the reverse complement on a
    // forward strand; a row without STRAND admits either.
    let reverse = ClippedRow {
        strand: "-1",
        ..plain
    };
    let unstranded = ClippedRow {
        strand: "",
        ..plain
    };
    let complemented = ("c.-149-3952TA>TT", "c.-149-3952A>T");
    let as_written = ("c.-149-3952TA>AA", "c.-149-3952T>A");
    assert!(clipped_accepts(
        "TA",
        "AA",
        complemented.0,
        complemented.1,
        false,
        reverse
    ));
    assert!(!clipped_accepts(
        "TA",
        "AA",
        as_written.0,
        as_written.1,
        false,
        reverse
    ));
    assert!(!clipped_accepts(
        "TA",
        "AA",
        complemented.0,
        complemented.1,
        false,
        plain
    ));
    for (reference, candidate) in [complemented, as_written] {
        assert!(
            clipped_accepts("TA", "AA", reference, candidate, false, unstranded),
            "{reference}"
        );
    }
    assert!(clipped_accepts(
        "A",
        "CCCA",
        "c.636+18_636+19insTGGG",
        "c.636+18_636+19insGGG",
        false,
        reverse
    ));
    // Where the allele spans bases the model holds and the genome lacks, USED_REF is longer
    // than GIVEN_REF: the clip reads the genome's two bases and leaves one position, and the
    // reference allele of a substitution keeps the model's bases on both sides.
    let longer_model = |given_ref: &'static str| ClippedRow {
        given_ref: Some(given_ref),
        ..reverse
    };
    assert!(clipped_accepts(
        "GGGAAG",
        "GA",
        "c.-1002CTTCCC>TC",
        "c.-1002CTTCC>T",
        false,
        longer_model("GG")
    ));
    assert!(clipped_accepts(
        "AATTT",
        "AAA",
        "c.2840delinsTTT",
        "c.2840delinsTT",
        false,
        longer_model("AC")
    ));
    assert!(!clipped_accepts(
        "GGGAAG",
        "GA",
        "c.-1002CTTCCC>TC",
        "c.-1002CTTCC>T",
        false,
        reverse
    ));
}

/// Two values whose coordinates differ stay a difference. With `renumbered` the rule
/// takes coordinates that differ as the class that renumbers an edited model's HGVSc reads
/// them, the same tokens over as many positions on a `BAM_EDIT` `OK` row, beside the
/// clipped alleles; coordinates over another number of positions, a `FAILED` row and a
/// candidate whose alleles are not the clipped pair stay differences.
#[test]
fn the_clipped_alleles_class_reads_renumbered_coordinates_only_when_asked() {
    let plain = ClippedRow::default();
    let (whole, clipped) = ("c.636+18_636+19insTGGG", "c.639+18_639+19insGGG");
    let renumbered = |used_ref: &str, allele: &str, reference: &str, candidate: &str| {
        clipped_accepts(used_ref, allele, reference, candidate, true, plain)
    };
    assert!(!clipped_accepts("T", "TGGG", whole, clipped, false, plain));
    assert!(renumbered("T", "TGGG", whole, clipped));
    assert!(renumbered("T", "TGGG", whole, "c.636+18_636+19insGGG"));
    assert!(renumbered("TTT", "TGT", "c.56+4TTT>TGT", "c.57+4T>G"));
    assert!(renumbered(
        "TAC",
        "TGG",
        "c.12_13delinsTGG",
        "c.15_16delinsGG"
    ));
    assert!(!renumbered("T", "TGGG", whole, "c.639+18_639+20insGGG"));
    assert!(!renumbered(
        "TAC",
        "TGG",
        "c.12_13delinsTGG",
        "c.15_17delinsGG"
    ));
    let failed = ClippedRow {
        status: "FAILED",
        ..plain
    };
    assert!(!clipped_accepts("T", "TGGG", whole, clipped, true, failed));
    assert!(!renumbered("T", "TGGG", whole, "c.639+18_639+19insGG"));
    assert!(!renumbered("T", "TGGG", whole, "c.639+18_639+19insTGGG"));
}

const CLAMPED_END: &str = "hgvs_shifted_insertion_clamped_end";

/// A default-format output of one shifted insertion's row of a non-coding transcript, its
/// cDNA position 14-15, with the notation `notation`.
fn clamped_reference(notation: &str) -> String {
    let extra = format!(
        "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000617755.1:{notation};HGVS_OFFSET=14;SHIFT_LENGTH=14"
    );
    let row = default_row(
        "21:25928754:C:CT",
        "21:25928767-25928768",
        "T",
        "ENST00000617755",
        "non_coding_transcript_exon_variant",
        &extra,
    );
    default_text(&[row.replacen("\t-\t-\t-\t-\t-\t-\t", "\t14-15\t-\t-\t-\t-\t-\t", 1)])
}

/// Under --shift_3prime 1 --shift_length, an insertion's row with a positive SHIFT_LENGTH and
/// an HGVS_OFFSET whose reference HGVSc is written at one position P (a substitution, a
/// delins, or an insertion of the bases before P at P-1_P) while the candidate's is a
/// duplication ending at P-1 or an insertion between P-1 and P, of as many bases as the
/// allele, is the class, counted once per row in every format, on a row whose cDNA position
/// starts at the shift's length; an empty HGVSc on either side, a zero shift, a missing offset,
/// a cDNA position elsewhere, a deletion's row, another transcript, another numbering, a
/// duplication ending elsewhere or of another length, other or more inserted bases, an
/// unclipped reference insertion and a substitution's row stay differences.
#[test]
fn the_clamped_end_class_accepts_the_shifted_notation_against_the_clamped_one() {
    let class = [CLAMPED_END];
    let reference = clamped_reference("n.15A>T");
    let candidate = |notation: &str| reference.replace("n.15A>T", notation);
    // The verdict of a text against itself with one notation in place of another.
    let swapped = |text: &str, from: &str, to: &str| -> String {
        verdict_of(text, &text.replace(from, to), &class)
    };

    let (report, documented) = compare_default(&reference, &candidate("n.14dup"), &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(CLAMPED_END, 1)]));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &json(r#"{"hgvs_shifted_insertion_clamped_end": 1}"#)
    );
    assert_eq!(verdict_of(&reference, &candidate("n.14dup"), &[]), "red");
    assert_eq!(swapped(&reference, "n.15A>T", "n.14_15insT"), "green");
    let two_bases = reference
        .replace("n.15A>T", "n.15delinsTT")
        .replace("\tT\t", "\tTT\t");
    assert_eq!(swapped(&two_bases, "n.15delinsTT", "n.13_14dup"), "green");
    // The 5' UTR numbering counts down to c.-1 and then c.1: the base before c.-1578 is c.-1579.
    let utr = reference.replace("n.15A>T", "c.-1578C>A");
    assert_eq!(swapped(&utr, "c.-1578C>A", "c.-1579dup"), "green");
    let first_coding = reference.replace("n.15A>T", "c.1A>T");
    assert_eq!(swapped(&first_coding, "c.1A>T", "c.-1dup"), "green");
    // The reference's insertion of the bases before P at P-1_P is the clamped shape too.
    let clipped = reference
        .replace("n.15A>T", "c.-258_-257insCGAGGGA")
        .replace("\tT\t", "\tCGAGGGAG\t");
    let clipped_notation = "c.-258_-257insCGAGGGA";
    assert_eq!(
        swapped(&clipped, clipped_notation, "c.-265_-258dup"),
        "green"
    );
    assert_eq!(
        swapped(&clipped, clipped_notation, "c.-258_-257insCGAGGGAG"),
        "green"
    );
    for notation in [
        "n.13dup",
        "n.15dup",
        "n.13_14dup",
        "n.13_14insT",
        "n.15_16insT",
        "n.14A>T",
        "n.14delinsTT",
        "c.14dup",
        "-",
    ] {
        assert_eq!(
            swapped(&reference, "n.15A>T", notation),
            "red",
            "{notation}"
        );
    }
    assert_eq!(swapped(&two_bases, "n.15delinsTT", "n.14dup"), "red");
    assert_eq!(swapped(&clipped, clipped_notation, "c.-264_-258dup"), "red");
    // A duplication spanning the c.-1 to c.1 boundary covers one base fewer than its bounds say.
    let boundary = reference
        .replace("n.15A>T", "c.2delinsTTT")
        .replace("\tT\t", "\tTTT\t");
    assert_eq!(swapped(&boundary, "c.2delinsTTT", "c.-2_1dup"), "green");
    assert_eq!(swapped(&boundary, "c.2delinsTTT", "c.-3_1dup"), "red");
    // A candidate insertion carries the reference's bases: the whole of a substitution's or a
    // delins's, one more than a clipped insertion's.
    for notation in ["n.14_15insG", "n.14_15insTT"] {
        assert_eq!(
            swapped(&reference, "n.15A>T", notation),
            "red",
            "{notation}"
        );
    }
    assert_eq!(swapped(&two_bases, "n.15delinsTT", "n.14_15insTT"), "green");
    assert_eq!(swapped(&two_bases, "n.15delinsTT", "n.14_15insTTT"), "red");
    for notation in ["c.-258_-257insCGAGGGAGG", "c.-258_-257insGCGAGGGAG"] {
        assert_eq!(
            swapped(&clipped, clipped_notation, notation),
            "red",
            "{notation}"
        );
    }
    // A reference insertion carrying the whole allele is not the clipped form: an ordinary shifted
    // insertion's notation, against which a longer candidate stays a difference.
    let whole = clipped.replace(clipped_notation, "c.-258_-257insCGAGGGAG");
    for notation in ["c.-265_-258dup", "c.-258_-257insCGAGGGAGA"] {
        assert_eq!(
            swapped(&whole, "c.-258_-257insCGAGGGAG", notation),
            "red",
            "{notation}"
        );
    }
    // The shift carried the insertion from before the first base, so its cDNA position starts at
    // the shift's length; a row shifted elsewhere stays a difference, as does a position of 0.
    let elsewhere = reference.replace("\t14-15\t", "\t20-21\t");
    assert_eq!(swapped(&elsewhere, "n.15A>T", "n.14dup"), "red");
    assert_eq!(swapped(&first_coding, "c.1A>T", "c.0dup"), "red");
    let own_notation = "HGVSc=ENST00000617755.1:n.15A>T;";
    let other_transcript =
        candidate("n.14dup").replace("HGVSc=ENST00000617755.1", "HGVSc=ENST00000617756.1");
    assert_eq!(verdict_of(&reference, &other_transcript, &class), "red");
    assert_eq!(verdict_of(&candidate("n.14dup"), &reference, &class), "red");
    assert_eq!(
        verdict_of(
            &reference.replace(own_notation, ""),
            &candidate("n.14dup"),
            &class
        ),
        "red"
    );
    assert_eq!(
        verdict_of(&reference, &reference.replace(own_notation, ""), &class),
        "red"
    );
    for (variant, without) in [
        ("SHIFT_LENGTH=14", "SHIFT_LENGTH=0"),
        (";HGVS_OFFSET=14", ""),
    ] {
        let unshifted = reference.replace(variant, without);
        assert_eq!(
            swapped(&unshifted, "n.15A>T", "n.14dup"),
            "red",
            "{variant}"
        );
    }
    let deletion = reference.replace("\tT\t", "\t-\t");
    assert_eq!(swapped(&deletion, "n.15A>T", "n.14dup"), "red");
    let substitution = reference.replace("21:25928767-25928768", "21:25928768");
    assert_eq!(swapped(&substitution, "n.15A>T", "n.14dup"), "red");
    // The VCF and JSON formats read the insertion from the record's alleles and count the row once.
    let inputs = vec!["21\t25928753\t21:25928754:C:CT\tC\tCT\t.\t.\t.".to_string()];
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|cDNA_position|HGVSc|HGVS_OFFSET|SHIFT_LENGTH\">\n";
    let record = &inputs[0][..inputs[0].len() - 1];
    let vcf_entries = |tail: &str| -> Vec<Entry> {
        let text = format!(
            "{csq_header}{record}CSQ=T|non_coding_transcript_exon_variant|ENST00000617755|Transcript|14-15|ENST00000617755.1:{tail}\n"
        );
        parse_vcf(&text, &inputs).expect("the text parses")
    };
    let mut documented = documented_with(&class);
    let report = compare(
        &vcf_entries("n.15A>T|14|14"),
        &vcf_entries("n.14dup|14|14"),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(CLAMPED_END, 1)]));
    let report = compare(
        &vcf_entries("n.15A>T|14|14"),
        &vcf_entries("n.14dup||14"),
        &mut documented_with(&class),
    );
    assert_eq!(verdict(&report), "red");

    let json_entries = |hgvsc: &str, shift_length: i64| -> Vec<Entry> {
        let text = dumped(&format!(
            r#"{{"input": "21\t25928753\t21:25928754:C:CT\tC\tCT\t.\t.\t.", "id": "21:25928754:C:CT", "allele_string": "-/T", "most_severe_consequence": "non_coding_transcript_exon_variant",
                "transcript_consequences": [{{"variant_allele": "T", "transcript_id": "ENST00000617755", "consequence_terms": ["non_coding_transcript_exon_variant"],
                 "cdna_start": 14, "cdna_end": 15, "hgvsc": "ENST00000617755.1:{hgvsc}", "hgvs_offset": 14, "shift_length": {shift_length}}}]}}"#
        )) + "\n";
        parse_json(&text, &inputs).expect("the text parses")
    };
    let mut documented = documented_with(&class);
    let report = compare(
        &json_entries("n.15A>T", 14),
        &json_entries("n.14dup", 14),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(CLAMPED_END, 1)]));
    let report = compare(
        &json_entries("n.15A>T", 14),
        &json_entries("n.14dup", 0),
        &mut documented_with(&class),
    );
    assert_eq!(verdict(&report), "red");
}

/// A reference sequence masked over the row's window gives the clamped substitution an
/// unknown reference base, since the notation reads the base after the position from the
/// transcript slice's sequence: such a substitution against the candidate's duplication or
/// insertion is the class, as one with a known base is, counted once per row. An unknown
/// alternate base is not a notation of the class, and the pair stays a difference without
/// the class.
#[test]
fn the_clamped_end_class_accepts_a_masked_reference_base() {
    let class = [CLAMPED_END];
    let reference = clamped_reference("n.15N>T");
    let duplication = reference.replace("n.15N>T", "n.14dup");
    let (report, documented) = compare_default(&reference, &duplication, &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(CLAMPED_END, 1)]));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &json(r#"{"hgvs_shifted_insertion_clamped_end": 1}"#)
    );
    assert_eq!(
        verdict_of(
            &reference,
            &reference.replace("n.15N>T", "n.14_15insT"),
            &class
        ),
        "green"
    );
    assert_eq!(verdict_of(&reference, &duplication, &[]), "red");
    assert_eq!(
        verdict_of(
            &reference.replace("n.15N>T", "n.15A>N"),
            &duplication,
            &class
        ),
        "red"
    );
}

/// Under --shift_3prime 1 --shift_length, a Consequence set the reference writes as the
/// candidate's plus intron_variant, on a row with a positive SHIFT_LENGTH whose HGVSc places
/// the shifted insertion or deletion at the acceptor's two bases or in the exon (the row then
/// without a splice-site term) is the class, counted once per row in every format and
/// leaving every other field compared; a zero shift, a deletion reaching an intron base away
/// from the splice site or spanning a whole intron, a row without an HGVSc unless an
/// insertion's carrying the acceptor term, the term on the candidate's side, a set differing
/// by another term or by two terms and a swap of the two engines stay differences.
#[test]
fn the_unshifted_span_class_accepts_an_intron_term_the_reference_alone_keeps() {
    let klass = "shift_intron_term_unshifted_span";
    let class = [klass];
    let acceptor_terms = "splice_acceptor_variant,intron_variant,non_coding_transcript_variant";
    let extra =
        "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000820008.1:n.469-2dup;HGVS_OFFSET=2;SHIFT_LENGTH=2";
    let reference = default_text(&[default_row(
        "21:36460365:G:GA",
        "21:36460367-36460368",
        "A",
        "ENST00000820008",
        acceptor_terms,
        extra,
    )]);
    let candidate = reference.replace(
        acceptor_terms,
        "splice_acceptor_variant,non_coding_transcript_variant",
    );
    // The verdict of a text against itself without the reference's intron term.
    let without_intron =
        |text: &str| -> String { verdict_of(text, &text.replace(",intron_variant", ""), &class) };
    let (report, documented) = compare_default(&reference, &candidate, &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
    assert_eq!(
        at(&report, &["columns", "Consequence", "masked"]),
        &json(r#"{"shift_intron_term_unshifted_span": 1}"#)
    );
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    assert_eq!(verdict_of(&candidate, &reference, &class), "red");
    assert_eq!(
        verdict_of(
            &reference,
            &candidate.replace("HGVS_OFFSET=2", "HGVS_OFFSET=1"),
            &class
        ),
        "red"
    );
    let unshifted = reference.replace("SHIFT_LENGTH=2", "SHIFT_LENGTH=0");
    assert_eq!(without_intron(&unshifted), "red");
    assert_eq!(
        verdict_of(
            &reference.replace(";SHIFT_LENGTH=2", ""),
            &candidate.replace(";SHIFT_LENGTH=2", ""),
            &class
        ),
        "red"
    );
    for terms in [
        "splice_acceptor_variant",
        "splice_acceptor_variant,intron_variant",
        "splice_acceptor_variant,splice_region_variant,non_coding_transcript_variant",
        "intron_variant,non_coding_transcript_variant",
    ] {
        assert_eq!(
            verdict_of(
                &reference,
                &reference.replace(acceptor_terms, terms),
                &class
            ),
            "red",
            "{terms}"
        );
    }
    // The shift may carry the span past the acceptor onto the exon's first bases, where the row
    // has no acceptor term: an exonic notation is the class on a row without a splice-site term,
    // in every coordinate shape the notation takes, and so is an insertion between the
    // acceptor's last base and the exon.
    let exon_row = |location: &str,
                    cdna: &str,
                    terms: &str,
                    notation: &str,
                    allele: &str,
                    shift: &str|
     -> String {
        let extra = format!(
            "IMPACT=LOW;STRAND=1;HGVSc=ENST00000435744.5:{notation};HGVS_OFFSET={shift};SHIFT_LENGTH={shift}"
        );
        let row = default_row(
            "rs544638589",
            location,
            allele,
            "ENST00000435744",
            terms,
            &extra,
        );
        default_text(&[row.replacen(
            "\t-\t-\t-\t-\t-\t-\t",
            &format!("\t{cdna}\t-\t-\t-\t-\t-\t"),
            1,
        )])
    };
    let exon_terms = "splice_region_variant,non_coding_transcript_exon_variant,intron_variant";
    let exon = exon_row("21:13065559", "1288", exon_terms, "n.1288del", "-", "8");
    assert_eq!(without_intron(&exon), "green");
    let non_coding = "splice_region_variant,intron_variant,non_coding_transcript_variant";
    let placed_rows = [
        exon_row(
            "21:13065559-13065563",
            "1288-1292",
            exon_terms,
            "n.1288_1292del",
            "-",
            "8",
        ),
        exon_row(
            "21:45504011-45504015",
            "2637-2641",
            "splice_region_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
            "c.*1390_*1394del",
            "-",
            "4",
        ),
        exon_row(
            "21:46125198-46125199",
            "1773-1774",
            "frameshift_variant,intron_variant,NMD_transcript_variant",
            "c.1773_1774insCAGGCCCCAGGCC",
            "CAGGCCCCAGGCC",
            "9",
        ),
        exon_row(
            "21:44406596-44406597",
            "3152-3153",
            "frameshift_variant,intron_variant",
            "c.2791-7_2793dup",
            "CCTGCAGATG",
            "10",
        ),
        exon_row(
            "21:46054102-46054103",
            "-",
            non_coding,
            "n.288-3_288-1dup",
            "CAG",
            "3",
        ),
        exon_row(
            "21:13065558-13065559",
            "-",
            non_coding,
            "n.1288-1_1288insA",
            "A",
            "9",
        ),
    ];
    for placed in &placed_rows {
        assert_eq!(without_intron(placed), "green", "{placed}");
    }
    for wider in [
        "n.1288-3_1292del",
        "n.1288-1_1292+1del",
        "n.1292+2del",
        "n.1288-3dup",
    ] {
        let row = exon_row(
            "21:13065559-13065563",
            "1288-1292",
            exon_terms,
            wider,
            "-",
            "8",
        );
        assert_eq!(without_intron(&row), "red", "{wider}");
    }
    // A deletion over a whole intron is numbered over the exons alone and covers both splice
    // sites; the term it carries is the shifted span's on both engines.
    for site_terms in [
        "splice_acceptor_variant,splice_donor_variant,",
        "splice_donor_variant,",
        "splice_acceptor_variant,",
    ] {
        let whole_intron = exon_row(
            "21:13065559-13065663",
            "1288-1292",
            &format!("{site_terms}{exon_terms}"),
            "n.1288_1292del",
            "-",
            "8",
        );
        assert_eq!(without_intron(&whole_intron), "red", "{site_terms}");
    }
    let unshifted_exon = exon.replace("SHIFT_LENGTH=8", "SHIFT_LENGTH=0");
    assert_eq!(without_intron(&unshifted_exon), "red");
    assert_eq!(
        verdict_of(&exon.replace(",intron_variant", ""), &exon, &class),
        "red"
    );
    assert_eq!(
        verdict_of(
            &exon,
            &exon.replace(
                "non_coding_transcript_exon_variant,intron_variant",
                "non_coding_transcript_exon_variant,splice_polypyrimidine_tract_variant"
            ),
            &class
        ),
        "red"
    );
    // A deletion shifted onto an acceptor base shows the same shape; one reaching an intron base
    // away from the splice site carries the term from its shifted span on both engines.
    let deletion = reference
        .replace("\tA\t", "\t-\t")
        .replace("21:36460367-36460368", "21:36460367")
        .replace("n.469-2dup", "n.469-2del")
        .replace(
            "HGVS_OFFSET=2;SHIFT_LENGTH=2",
            "HGVS_OFFSET=1;SHIFT_LENGTH=1",
        );
    assert_eq!(without_intron(&deletion), "green");
    for notation in [
        "n.469-3_469-2del",
        "n.469-10_469-1del",
        "n.469-2_469-1delinsA",
        "n.100+2del",
    ] {
        let wider = deletion.replace("n.469-2del", notation);
        assert_eq!(without_intron(&wider), "red", "{notation}");
    }
    // A duplication's first coordinate runs back into the tract; a 5' UTR exon's own minus is no offset.
    for notation in ["n.469-6_469-2dup", "c.-3-1_-3insA", "c.-3-2del"] {
        let placed = deletion.replace("n.469-2del", notation);
        assert_eq!(without_intron(&placed), "green", "{notation}");
    }
    // Without an HGVSc the allele and the terms decide: an insertion's row carrying the acceptor
    // term is read at the acceptor, a deletion's and an insertion's without the term are not.
    let bare = reference.replace("HGVSc=ENST00000820008.1:n.469-2dup;HGVS_OFFSET=2;", "");
    assert_eq!(without_intron(&bare), "green");
    let bare_deletion = deletion.replace("HGVSc=ENST00000820008.1:n.469-2del;HGVS_OFFSET=1;", "");
    assert_eq!(without_intron(&bare_deletion), "red");
    let bare_exon = bare.replace(
        acceptor_terms,
        "splice_region_variant,non_coding_transcript_exon_variant,intron_variant",
    );
    assert_eq!(without_intron(&bare_exon), "red");
    let inputs = vec!["21\t36460365\t21:36460365:G:GA\tG\tGA\t.\t.\t.".to_string()];
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|HGVSc|HGVS_OFFSET|SHIFT_LENGTH\">\n";
    let record = &inputs[0][..inputs[0].len() - 1];
    let vcf_entries = |intron: &str| -> Vec<Entry> {
        let text = format!(
            "{csq_header}{record}CSQ=A|splice_acceptor_variant&{intron}non_coding_transcript_variant|ENST00000820008|Transcript|ENST00000820008.1:n.469-2dup|2|2\n"
        );
        parse_vcf(&text, &inputs).expect("the text parses")
    };
    let mut documented = documented_with(&class);
    let report = compare(
        &vcf_entries("intron_variant&"),
        &vcf_entries(""),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));

    let json_entries = |terms: &str| -> Vec<Entry> {
        let text = dumped(&format!(
            r#"{{"input": "21\t36460365\t21:36460365:G:GA\tG\tGA\t.\t.\t.", "id": "21:36460365:G:GA", "allele_string": "-/A", "most_severe_consequence": "splice_acceptor_variant",
                "transcript_consequences": [{{"variant_allele": "A", "transcript_id": "ENST00000820008", "consequence_terms": [{terms}],
                 "hgvsc": "ENST00000820008.1:n.469-2dup", "hgvs_offset": 2, "shift_length": 2}}]}}"#
        )) + "\n";
        parse_json(&text, &inputs).expect("the text parses")
    };
    let mut documented = documented_with(&class);
    let report = compare(
        &json_entries(
            r#""splice_acceptor_variant", "intron_variant", "non_coding_transcript_variant""#,
        ),
        &json_entries(r#""splice_acceptor_variant", "non_coding_transcript_variant""#),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
}

const CSQ_DESCRIPTION: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: ";

/// The CSQ header line of a VCF output whose entries carry `fields`.
fn csq_header(fields: &str) -> String {
    format!("{CSQ_DESCRIPTION}{fields}\">\n")
}

/// The input record `line` as an output line: its INFO column, the last, left to follow.
fn record_columns(line: &str) -> &str {
    &line[..line.len() - 1]
}

fn tab_entries(text: &str) -> Vec<Entry> {
    parse_tab(text).expect("the text parses")
}

/// The report's verdict and what was applied of two entry lists under one class.
fn under_class(
    klass: &str,
    expected: &[Entry],
    actual: &[Entry],
) -> (String, std::collections::BTreeMap<String, i64>) {
    let mut documented = documented_with(&[klass]);
    let report = compare(expected, actual, &mut documented);
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// The tab column, the VCF value (`+` inside an item, `&` between items) and the JSON
/// string are taken alike, and a changed item stays a difference in each.
#[test]
fn the_exon_item_order_class_reads_the_tab_vcf_and_json_spellings() {
    let klass = "nearest_exon_jb_item_order";
    let once = ("green".to_string(), count_map(&[(klass, 1)]));
    let tab_header =
        "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tNearestExonJB";
    let tab = format!(
        "{tab_header}\ndel\t21:101\t-\tENST1\tTranscript\tintron_variant\tENSE2|3|end|207,ENSE1|3|start|117\n"
    );
    let swapped = tab.replace(
        "ENSE2|3|end|207,ENSE1|3|start|117",
        "ENSE1|3|start|117,ENSE2|3|end|207",
    );
    assert_eq!(
        under_class(klass, &tab_entries(&tab), &tab_entries(&swapped)),
        once
    );
    let changed = tab_entries(&swapped.replace("|117", "|118"));
    assert_eq!(under_class(klass, &tab_entries(&tab), &changed).0, "red");
    let inputs = vec!["21\t100\tdel\tAT\tA\t.\t.\t.".to_string()];
    let vcf_entries = |text: &str| parse_vcf(text, &inputs).expect("the text parses");
    let vcf_ref = format!(
        "{}{}CSQ=-|intron_variant|ENST1|Transcript|ENSE2+3+end+207&ENSE1+3+start+117\n",
        csq_header("Allele|Consequence|Feature|Feature_type|NearestExonJB"),
        record_columns(&inputs[0])
    );
    let vcf_cand = vcf_ref.replace(
        "ENSE2+3+end+207&ENSE1+3+start+117",
        "ENSE1+3+start+117&ENSE2+3+end+207",
    );
    assert_eq!(
        under_class(klass, &vcf_entries(&vcf_ref), &vcf_entries(&vcf_cand)),
        once
    );
    let changed = vcf_entries(&vcf_cand.replace("+117", "+118"));
    assert_eq!(
        under_class(klass, &vcf_entries(&vcf_ref), &changed).0,
        "red"
    );

    let json_entries = |text: &str| parse_json(text, &inputs).expect("the text parses");
    let object = |items: &str| {
        dumped(&format!(
            r#"{{"input": "21\t100\tdel\tAT\tA\t.\t.\t.", "id": "del", "allele_string": "T/-", "most_severe_consequence": "intron_variant", "transcript_consequences": [
                {{"variant_allele": "-", "transcript_id": "ENST1", "consequence_terms": ["intron_variant"], "nearestexonjb": "{items}"}}]}}"#
        )) + "\n"
    };
    let json_ref = object("ENSE2|3|end|207,ENSE1|3|start|117");
    let json_cand = object("ENSE1|3|start|117,ENSE2|3|end|207");
    assert_eq!(
        under_class(klass, &json_entries(&json_ref), &json_entries(&json_cand)),
        once
    );
    let changed = json_entries(&json_cand.replace("|117", "|118"));
    assert_eq!(
        under_class(klass, &json_entries(&json_ref), &changed).0,
        "red"
    );
}

/// The seven pairs of a substitution in the 5' UTR of one transcript that destroys an upstream
/// start codon, and one order Ensembl VEP writes them in; over repeated runs of one
/// command it writes them in other orders.
const UTR_ANNOTATION_PAIRS: [(&str, &str); 7] = [
    ("CapDistanceToStart", "407"),
    ("DistanceToCDS", "23"),
    ("DistanceToStop", "24"),
    ("Evidence", "False"),
    ("KozakContext", "ATAATGG"),
    ("KozakStrength", "Strong"),
    ("type", "OutOfFrame_oORF"),
];
const UTR_ANNOTATION_HASH_ORDER: [&str; 7] = [
    "DistanceToStop",
    "CapDistanceToStart",
    "KozakContext",
    "KozakStrength",
    "DistanceToCDS",
    "Evidence",
    "type",
];

/// One `5UTR_annotation` group: the pairs of `keys` in that order.
fn utr_annotation(keys: &[&str]) -> String {
    let pairs: Vec<String> = keys
        .iter()
        .map(|key| {
            let value = UTR_ANNOTATION_PAIRS
                .iter()
                .find(|(name, _)| name == key)
                .expect("a pair of the annotation")
                .1;
            format!("{key}={value}")
        })
        .collect();
    pairs.join(":")
}

/// The tab column and the VCF value (`&` between groups) are taken alike; the JSON
/// format writes an object, which is compared with its keys sorted, so the class applies to
/// none and a changed value there stays a difference.
#[test]
fn the_annotation_pair_order_class_reads_the_tab_and_vcf_spellings_and_no_json_object() {
    let klass = "utrannotator_annotation_pair_order";
    let once = ("green".to_string(), count_map(&[(klass, 1)]));
    let hash_order = utr_annotation(&UTR_ANNOTATION_HASH_ORDER);
    let mut sorted_keys = UTR_ANNOTATION_HASH_ORDER;
    sorted_keys.sort_unstable();
    let key_order = utr_annotation(&sorted_keys);
    let tab_header =
        "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\t5UTR_annotation";
    let tab = format!(
        "{tab_header}\nsnv\t21:100\tG\tENST1\tTranscript\t5_prime_UTR_variant\t{hash_order}\n"
    );
    assert_eq!(
        under_class(
            klass,
            &tab_entries(&tab),
            &tab_entries(&tab.replace(&hash_order, &key_order))
        ),
        once
    );
    let inputs = vec!["21\t100\tsnv\tA\tG\t.\t.\t.".to_string()];
    let vcf_entries = |text: &str| parse_vcf(text, &inputs).expect("the text parses");
    let vcf_ref = format!(
        "{}{}CSQ=G|5_prime_UTR_variant|ENST1|Transcript|k=2:t=v&t=u:k=1\n",
        csq_header("Allele|Consequence|Feature|Feature_type|5UTR_annotation"),
        record_columns(&inputs[0])
    );
    let vcf_cand = vcf_ref.replace("k=2:t=v&t=u:k=1", "k=1:t=u&k=2:t=v");
    assert_eq!(
        under_class(klass, &vcf_entries(&vcf_ref), &vcf_entries(&vcf_cand)),
        once
    );
    let changed = vcf_entries(&vcf_cand.replace("k=1", "k=3"));
    assert_eq!(
        under_class(klass, &vcf_entries(&vcf_ref), &changed).0,
        "red"
    );

    let json_entries = |annotation: &str| {
        let text = dumped(&format!(
            r#"{{"input": "21\t100\tsnv\tA\tG\t.\t.\t.", "id": "snv", "allele_string": "A/G", "most_severe_consequence": "5_prime_UTR_variant", "transcript_consequences": [
                {{"variant_allele": "G", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"], "5utr_annotation": {{"1": {annotation}}}}}]}}"#
        )) + "\n";
        parse_json(&text, &inputs).expect("the text parses")
    };
    let json_ref = json_entries(r#"{"type": "uORF", "KozakStrength": "Weak"}"#);
    let json_cand = json_entries(r#"{"KozakStrength": "Weak", "type": "uORF"}"#);
    assert_eq!(
        under_class(klass, &json_ref, &json_cand),
        ("green".to_string(), count_map(&[]))
    );
    let changed = json_entries(r#"{"KozakStrength": "Strong", "type": "uORF"}"#);
    assert_eq!(under_class(klass, &json_ref, &changed).0, "red");
}

/// Under --merged Ensembl VEP keeps one transcript per stable identifier among the models a
/// buffer loads, so a model the cache carries at two loci loses its rows at one of them: a
/// Transcript row the candidate alone writes, SOURCE RefSeq (or none, under --refseq), on a
/// record the reference wrote a transcript row for is the class in every format, and the JSON
/// record's most_severe_consequence is left uncompared; an Ensembl-sourced row, a row on a
/// record the reference wrote nothing or an intergenic row for, and a regulatory row stay
/// differences.
#[test]
fn the_collapse_class_accepts_a_transcript_row_the_candidate_alone_writes_on_an_annotated_record() {
    let klass = "merged_same_stable_id_collapse";
    let class = [klass];
    let model = "compmerge.264.Heart.chr21";
    let exon_term = "non_coding_transcript_exon_variant";
    let kept_copy = default_row(
        "rs1036535345",
        "21:32612111",
        "C",
        model,
        exon_term,
        "IMPACT=MODIFIER;STRAND=-1;SOURCE=RefSeq",
    );
    let ensembl_model = default_row(
        "rs745798146",
        "21:39462549",
        "T",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl",
    );
    let lost_copy = default_row(
        "rs745798146",
        "21:39462549",
        "T",
        model,
        exon_term,
        "IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq",
    );
    let reference = default_text(&[&kept_copy, &ensembl_model]);
    let candidate = default_text(&[&kept_copy, &ensembl_model, &lost_copy]);
    let (report, documented) = compare_default(&reference, &candidate, &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
    assert_eq!(
        texts(at(&report, &["mask_examples", klass])),
        ["21:39462549|T|compmerge.264.Heart.chr21|Transcript: candidate-only row, 'non_coding_transcript_exon_variant'"]
    );
    let (plain, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(number(&plain, &["one_sided", "actual_only"]), 1);
    // A RefSeq cache writes no SOURCE: the row is the class with the field absent, where a row
    // carries the REFSEQ_MATCH that only `--refseq` and `--merged` write (Constants.pm:97-98).
    let without_source = |text: &str| {
        text.replace(";SOURCE=RefSeq", "")
            .replace(";SOURCE=Ensembl", "")
    };
    let matched = |text: &str| text.replace("STRAND=-1", "STRAND=-1;REFSEQ_MATCH=rseq_mrna_match");
    assert_eq!(
        verdict_of(
            &matched(&without_source(&reference)),
            &matched(&without_source(&candidate)),
            &class
        ),
        "green"
    );
    // Without either field the run is an Ensembl cache's, so the row stays a difference.
    assert_eq!(
        verdict_of(
            &without_source(&reference),
            &without_source(&candidate),
            &class
        ),
        "red"
    );
    // What a candidate that stays a difference comes to: red, its rows without a counterpart
    // and what the class counted.
    let refused = |reference: &str, candidate: &str, actual_only: i64, applied: &[(&str, i64)]| {
        let (red, documented) = compare_default(reference, candidate, &class);
        assert_eq!(verdict(&red), "red");
        assert_eq!(number(&red, &["one_sided", "actual_only"]), actual_only);
        assert_eq!(counts(&documented.applied), count_map(applied));
        red
    };
    // An Ensembl-sourced model is unique per locus, so a candidate-only row of one is a difference.
    let ensembl_row = default_text(&[
        &kept_copy,
        &ensembl_model,
        &lost_copy.replace("STRAND=1;SOURCE=RefSeq", "STRAND=1;SOURCE=Ensembl"),
    ]);
    refused(&reference, &ensembl_row, 1, &[]);
    // A record the reference wrote nothing for is another shape, and so is a candidate-only
    // regulatory row.
    refused(&default_text(&[&kept_copy]), &candidate, 2, &[]);
    let regulatory = default_row_of(
        "rs745798146",
        "21:39462549",
        "T",
        "ENSR1",
        "regulatory_region_variant",
        "IMPACT=MODIFIER;SOURCE=RefSeq",
        "RegulatoryFeature",
        "-",
    );
    refused(
        &reference,
        &default_text(&[&kept_copy, &ensembl_model, &regulatory]),
        1,
        &[],
    );
    // The collapse drops a model's rows and never adds one: another allele's row and a duplicate
    // of a row of a model the reference kept on the record stay differences.
    let other_allele = default_row(
        "rs1036535345",
        "21:32612111",
        "G",
        model,
        exon_term,
        "IMPACT=MODIFIER;STRAND=-1;SOURCE=RefSeq",
    );
    refused(
        &reference,
        &default_text(&[&kept_copy, &other_allele, &ensembl_model]),
        1,
        &[],
    );
    refused(
        &reference,
        &default_text(&[&kept_copy, &kept_copy, &ensembl_model]),
        1,
        &[],
    );
    // Two input records at one site with the same alternate allele share a text-format key, so
    // the second record's lost copy carries the duplicate counter; both are the class (their
    // Uploaded_variation tells the records apart), where a repeat of one record's row is not.
    let second_record = default_row(
        "rs745798146_second",
        "21:39462549",
        "T",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl",
    );
    let second_lost = default_row(
        "rs745798146_second",
        "21:39462549",
        "T",
        model,
        exon_term,
        "IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq",
    );
    let (report, documented) = compare_default(
        &default_text(&[&kept_copy, &ensembl_model, &second_record]),
        &default_text(&[
            &kept_copy,
            &ensembl_model,
            &lost_copy,
            &second_record,
            &second_lost,
        ]),
        &class,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 2)]));
    refused(
        &reference,
        &default_text(&[&kept_copy, &ensembl_model, &lost_copy, &lost_copy]),
        1,
        &[(klass, 1)],
    );
    // A record whose every overlapping model is a lost copy gets an intergenic row from the
    // reference, no transcript row: neither side's row is the class.
    let intergenic = default_row_of(
        "rs745798146",
        "21:39462549",
        "T",
        "-",
        "intergenic_variant",
        "IMPACT=MODIFIER",
        "-",
        "-",
    );
    let red = refused(
        &default_text(&[&kept_copy, &intergenic]),
        &default_text(&[&kept_copy, &lost_copy]),
        1,
        &[],
    );
    assert_eq!(number(&red, &["one_sided", "expected_only"]), 1);
    // The VCF and JSON formats key the row by record: one CSQ entry, one consequence object.
    let inputs = vec![
        "21\t32612111\trs1036535345\tA\tC\t.\t.\t.".to_string(),
        "21\t39462549\trs745798146\tC\tT\t.\t.\t.".to_string(),
    ];
    let kept_csq = format!("C|{exon_term}|MODIFIER|Transcript|{model}|-1|RefSeq");
    let ensembl_csq = "T|intron_variant|MODIFIER|Transcript|ENST1|1|Ensembl";
    let lost_csq = format!("T|{exon_term}|MODIFIER|Transcript|{model}|1|RefSeq");
    let vcf_entries = |second_record_csq: &str| -> Vec<Entry> {
        let header = csq_header("Allele|Consequence|IMPACT|Feature_type|Feature|STRAND|SOURCE");
        let text = format!(
            "{header}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}CSQ={kept_csq}\n{}CSQ={second_record_csq}\n",
            record_columns(&inputs[0]),
            record_columns(&inputs[1])
        );
        parse_vcf(&text, &inputs).expect("the text parses")
    };
    let both_rows = format!("{ensembl_csq},{lost_csq}");
    assert_eq!(
        under_class(klass, &vcf_entries(ensembl_csq), &vcf_entries(&both_rows)),
        ("green".to_string(), count_map(&[(klass, 1)]))
    );
    let plain = compare(
        &vcf_entries(ensembl_csq),
        &vcf_entries(&both_rows),
        &mut no_documented(),
    );
    assert_eq!(verdict(&plain), "red");

    let object = |i: usize, severe: &str, consequences: &[&str]| -> String {
        let input = inputs[i].replace('\t', "\\t");
        dumped(&format!(
            r#"{{"input": "{input}", "id": "r{i}", "most_severe_consequence": "{severe}", "transcript_consequences": [{}]}}"#,
            consequences.join(", ")
        )) + "\n"
    };
    let kept_object = format!(
        r#"{{"variant_allele": "C", "transcript_id": "{model}", "consequence_terms": ["{exon_term}"], "strand": -1, "source": "RefSeq"}}"#
    );
    let ensembl_object = r#"{"variant_allele": "T", "transcript_id": "ENST1", "consequence_terms": ["intron_variant"], "strand": 1, "source": "Ensembl"}"#;
    let lost_object = format!(
        r#"{{"variant_allele": "T", "transcript_id": "{model}", "consequence_terms": ["{exon_term}"], "strand": 1, "source": "RefSeq"}}"#
    );
    let reference_json =
        object(0, exon_term, &[&kept_object]) + &object(1, "intron_variant", &[ensembl_object]);
    // The lost copy's term outranks the kept row's, so the candidate's record summary differs too.
    let candidate_json = object(0, exon_term, &[&kept_object])
        + &object(1, exon_term, &[ensembl_object, &lost_object]);
    let json_entries = |text: &str| parse_json(text, &inputs).expect("the text parses");
    assert_eq!(
        under_class(
            klass,
            &json_entries(&reference_json),
            &json_entries(&candidate_json)
        ),
        ("green".to_string(), count_map(&[(klass, 1)]))
    );
    let bare = compare(
        &json_entries(&reference_json),
        &json_entries(&candidate_json),
        &mut no_documented(),
    );
    assert_eq!(verdict(&bare), "red");
    assert_eq!(number(&bare, &["one_sided", "actual_only"]), 1);
    assert_eq!(
        number(&bare, &["columns", "most_severe_consequence", "differing"]),
        1
    );
}

/// `chr13 60000000 A [CHR17:7675000[A` on an Ensembl cache, no `--merged`: the CSQ format carries
/// SYMBOL_SOURCE and no SOURCE or REFSEQ_MATCH, so the candidate's mate rows inside TP53 stay
/// unexpected under the collapse class as without it.
#[test]
fn the_collapse_class_leaves_candidate_only_rows_of_a_run_without_refseq_fields() {
    let klass = "merged_same_stable_id_collapse";
    let inputs = vec!["chr13\t60000000\tbnd2\tA\t[CHR17:7675000[A\t.\tPASS\t.".to_string()];
    let own_row = "[CHR17:7675000[A|intron_variant|MODIFIER|PCDH20|Transcript|ENST00000409186|HGNC";
    let mate_rows = [
        "[CHR17:7675000[A|intron_variant|MODIFIER|TP53|Transcript|ENST00000269305|HGNC",
        "[CHR17:7675000[A|intron_variant|MODIFIER|TP53|Transcript|ENST00000359597|HGNC",
    ];
    let vcf_entries = |csq: &str| -> Vec<Entry> {
        let header =
            csq_header("Allele|Consequence|IMPACT|SYMBOL|Feature_type|Feature|SYMBOL_SOURCE");
        let text = format!(
            "{header}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}CSQ={csq}\n",
            record_columns(&inputs[0])
        );
        parse_vcf(&text, &inputs).expect("the text parses")
    };
    let reference = vcf_entries(own_row);
    let candidate = vcf_entries(&format!("{own_row},{}", mate_rows.join(",")));
    for masks in [&[klass][..], &[]] {
        let mut documented = documented_with(masks);
        let report = compare(&reference, &candidate, &mut documented);
        assert_eq!(verdict(&report), "red");
        assert_eq!(number(&report, &["one_sided", "actual_only"]), 2);
        assert_eq!(counts(&documented.applied), count_map(&[]));
    }
}

const GNOMAD_VEP_CUT: &str = "G|synonymous_variant|LOW|MT-ND1|ENSG00000198888|Transcript|ENST00000361390|protein_coding|1/1||ENST00000361390.2:c.135A>G|ENSP00000354687.2:p.Leu45";
const GNOMAD_VEP_TAIL: &str = "=|135|135|45|L|ctA/ctG|1||1|SNV||HGNC|HGNC:7455|YES";

fn gnomad_vep_row(value: &str) -> String {
    let extra = if value.is_empty() {
        "IMPACT=LOW;STRAND=1".to_string()
    } else {
        format!("IMPACT=LOW;STRAND=1;gnomAD_vep={value}")
    };
    default_text(&[default_row(
        "nd1_3441_A_G",
        "MT:3441",
        "G",
        "ENST00000361390",
        "synonymous_variant",
        &extra,
    )])
}

/// A mitochondrial substitution under fields=vep: the file's vep value holds = inside its
/// protein notation, the reference's value stops before that sign and the candidate writes
/// the whole value; the percent-encoded sign and the JSON gnomad_vep key read the same.
#[test]
fn the_second_equals_class_accepts_a_value_cut_at_its_equals_sign() {
    let klass = "gnomadmt_info_value_second_equals";
    let whole = format!("{GNOMAD_VEP_CUT}{GNOMAD_VEP_TAIL}");
    let once = count_map(&[(klass, 1)]);
    let (report, documented) = compare_default(
        &gnomad_vep_row(GNOMAD_VEP_CUT),
        &gnomad_vep_row(&whole),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(
        at(&report, &["columns", "gnomAD_vep", "masked"]),
        &json(r#"{"gnomadmt_info_value_second_equals": 1}"#)
    );
    assert_eq!(
        verdict_of(
            &gnomad_vep_row(GNOMAD_VEP_CUT),
            &gnomad_vep_row(&whole),
            &[]
        ),
        "red"
    );
    let encoded = gnomad_vep_row(&whole.replace('=', "%3D"));
    assert_eq!(
        verdict_of(&gnomad_vep_row(GNOMAD_VEP_CUT), &encoded, &[klass]),
        "green"
    );
    let inputs = vec!["MT\t3441\tnd1_3441_A_G\tA\tG\t.\t.\t.".to_string()];
    let record = |value: &str| -> Vec<Entry> {
        let text = dumped(&format!(
            r#"{{"input": "MT\t3441\tnd1_3441_A_G\tA\tG\t.\t.\t.", "id": "nd1_3441_A_G", "most_severe_consequence": "synonymous_variant", "transcript_consequences": [
                {{"variant_allele": "G", "transcript_id": "ENST00000361390", "consequence_terms": ["synonymous_variant"], "gnomad_vep": "{value}"}}]}}"#
        )) + "\n";
        parse_json(&text, &inputs).expect("the text parses")
    };
    assert_eq!(
        under_class(klass, &record(GNOMAD_VEP_CUT), &record(&whole)),
        ("green".to_string(), once)
    );
}

/// Default-format rows of one record on one transcript: (allele, Consequence, the LoF
/// pairs of the Extra text) each.
fn lof_rows(name: &str, location: &str, feature: &str, rows: &[(&str, &str, &str)]) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, consequence, lof)| {
            let extra = if lof.is_empty() {
                "IMPACT=HIGH;STRAND=1".to_string()
            } else {
                format!("IMPACT=HIGH;STRAND=1;{lof}")
            };
            default_row(name, location, allele, feature, consequence, &extra)
        })
        .collect();
    default_text(&rows)
}

/// A record of two alternates on one transcript: one row is a stop gain whose base is the
/// ancestral one, which the reference leaves unfiltered beside a sibling allele and the
/// candidate filters ANC_ALLELE; where the item is the only filter the row's LoF reads HC
/// against LC, and the row is still counted once. The VCF format joins the items with &.
#[test]
fn the_ancestral_allele_class_accepts_the_filter_the_candidate_alone_applies() {
    let klass = "loftee_anc_allele_multiallelic_record";
    let once = count_map(&[(klass, 1)]);
    let masked = json(r#"{"loftee_anc_allele_multiallelic_record": 1}"#);
    let sibling = ("A", "synonymous_variant", "");
    let record = |lof: &str| {
        lof_rows(
            "21_32314600_C/T/A",
            "21:32314600",
            "ENST00000339944",
            &[("T", "stop_gained", lof), sibling],
        )
    };
    let reference = record("LoF=LC;LoF_filter=END_TRUNC");
    let candidate = record("LoF=LC;LoF_filter=END_TRUNC,ANC_ALLELE");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(at(&report, &["columns", "LoF_filter", "masked"]), &masked);
    assert!(at(&report, &["columns", "LoF"]).get("masked").is_none());
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let (report, documented) = compare_default(
        &record("LoF=HC"),
        &record("LoF=LC;LoF_filter=ANC_ALLELE"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(at(&report, &["columns", "LoF", "masked"]), &masked);
    assert_eq!(at(&report, &["columns", "LoF_filter", "masked"]), &masked);
    let inputs = vec!["21\t32314600\ttwo\tC\tT,A\t.\t.\t.".to_string()];
    let vcf = format!(
        "{}{}CSQ=T|stop_gained|ENST00000339944|Transcript|LC|END_TRUNC,A|synonymous_variant|ENST00000339944|Transcript||\n",
        csq_header("Allele|Consequence|Feature|Feature_type|LoF|LoF_filter"),
        record_columns(&inputs[0])
    );
    let vcf_entries = |text: &str| parse_vcf(text, &inputs).expect("the text parses");
    assert_eq!(
        under_class(
            klass,
            &vcf_entries(&vcf),
            &vcf_entries(&vcf.replace("END_TRUNC", "END_TRUNC&ANC_ALLELE"))
        ),
        ("green".to_string(), once)
    );
}

/// The ANC_ALLELE item stays a difference on a one-allele record, on default-format rows
/// under an identifier, on the reference's side, beside another differing item, on an
/// insertion's row and on a row outside LoFTEE's terms; a LoF that differs without the item,
/// or does not read LC beside it and HC without any, stays one too.
#[test]
fn the_ancestral_allele_class_refuses_every_other_filter_difference() {
    let klass = "loftee_anc_allele_multiallelic_record";
    let sibling = ("A", "synonymous_variant", "");
    // The verdict of the first allele's row under two LoF texts beside `siblings`, the class
    // counted exactly where the verdict is green.
    let colour = |reference_lof: &str,
                  candidate_lof: &str,
                  siblings: &[(&str, &str, &str)],
                  name: &str,
                  location: &str,
                  consequence: &str|
     -> String {
        let rows = |lof: &str| {
            let mut rows = vec![("T", consequence, lof)];
            rows.extend_from_slice(siblings);
            lof_rows(name, location, "ENST00000339944", &rows)
        };
        let (report, documented) =
            compare_default(&rows(reference_lof), &rows(candidate_lof), &[klass]);
        let applied = if verdict(&report) == "red" {
            count_map(&[])
        } else {
            count_map(&[(klass, 1)])
        };
        assert_eq!(counts(&documented.applied), applied);
        verdict(&report).to_string()
    };
    let name = "21_32314600_C/T/A";
    let location = "21:32314600";
    let own = |reference_lof: &str, candidate_lof: &str| {
        colour(
            reference_lof,
            candidate_lof,
            &[sibling],
            name,
            location,
            "stop_gained",
        )
    };
    let filtered = "LoF=LC;LoF_filter=END_TRUNC,ANC_ALLELE";
    let unfiltered = "LoF=LC;LoF_filter=END_TRUNC";
    assert_eq!(own(unfiltered, filtered), "green");
    assert_eq!(
        colour(
            unfiltered,
            filtered,
            &[],
            "21_32314600_C/T",
            location,
            "stop_gained"
        ),
        "red"
    );
    assert_eq!(
        colour(
            unfiltered,
            filtered,
            &[sibling],
            "two",
            location,
            "stop_gained"
        ),
        "red"
    );
    // The record's other allele on another feature alone is no sibling row of this one.
    let elsewhere = lof_rows(
        name,
        location,
        "ENST00000339944",
        &[("T", "stop_gained", unfiltered)],
    ) + &default_row(
        name,
        location,
        "A",
        "ENST2",
        "synonymous_variant",
        "IMPACT=LOW;STRAND=1",
    ) + "\n";
    let (report, documented) = compare_default(
        &elsewhere,
        &elsewhere.replace(unfiltered, filtered),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert!(counts(&documented.applied).is_empty());
    assert_eq!(own(filtered, unfiltered), "red");
    assert_eq!(own(unfiltered, "LoF=LC;LoF_filter=ANC_ALLELE"), "red");
    assert_eq!(
        own(
            unfiltered,
            "LoF=LC;LoF_filter=END_TRUNC,ANC_ALLELE,ANC_ALLELE"
        ),
        "red"
    );
    assert_eq!(
        colour(
            unfiltered,
            filtered,
            &[("TA", "frameshift_variant", "")],
            "21_32314601_-/T/TA",
            "21:32314600-32314601",
            "stop_gained"
        ),
        "red"
    );
    assert_eq!(
        colour(
            unfiltered,
            filtered,
            &[sibling],
            name,
            location,
            "missense_variant"
        ),
        "red"
    );
    assert_eq!(
        own("LoF=HC;LoF_filter=END_TRUNC", "LoF=LC;LoF_filter=END_TRUNC"),
        "red"
    );
    assert_eq!(own("LoF=LC", "LoF=HC;LoF_filter=ANC_ALLELE"), "red");
    assert_eq!(own("LoF=HC", "LoF=HC;LoF_filter=ANC_ALLELE"), "red");
}

const DONOR_FILTERED: &str =
    "LoF=LC;LoF_filter=GC_TO_GT_DONOR;LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:261";
const DONOR_SPARED: &str = "LoF=HC;LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:261";
const DONOR_RECORD: &str = "21_46126243_C/T/A";

/// Default-format rows of one record of two alternates on one transcript, its alleles spelt
/// by the Uploaded_variation `name`: (allele, the LoF pairs of the Extra text) each.
fn donor_rows(rows: &[(&str, &str)], consequence: &str, name: &str) -> String {
    let rows: Vec<(&str, &str, &str)> = rows
        .iter()
        .map(|(allele, lof)| (*allele, consequence, *lof))
        .collect();
    lof_rows(name, "21:46126243", "ENST00001028357", &rows)
}

/// A record of two alternates on a transcript whose intron opens with the bases the filter
/// reads: the reference filters the second alternate's row GC_TO_GT_DONOR by the first
/// alternate, the candidate writes HC and no filter; with the alleles in the other order the
/// reference spares the row and the candidate filters it. Either way the row is counted
/// once, its LoF with its LoF_filter.
#[test]
fn the_first_alternate_donor_class_accepts_the_filter_read_from_the_first_alternate() {
    let klass = "loftee_first_alt_donor_filter";
    let once = count_map(&[(klass, 1)]);
    let masked = json(r#"{"loftee_first_alt_donor_filter": 1}"#);
    let donor = "splice_donor_variant";
    let reference = donor_rows(
        &[("T", DONOR_FILTERED), ("A", DONOR_FILTERED)],
        donor,
        DONOR_RECORD,
    );
    let candidate = donor_rows(
        &[("T", DONOR_FILTERED), ("A", DONOR_SPARED)],
        donor,
        DONOR_RECORD,
    );
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(at(&report, &["columns", "LoF", "masked"]), &masked);
    assert_eq!(at(&report, &["columns", "LoF_filter", "masked"]), &masked);
    assert!(texts(at(&report, &["mask_examples", klass]))[0]
        .starts_with("21:46126243|A|ENST00001028357|Transcript: "));
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let other_order = "21_46126243_C/A/T";
    let spared_first = donor_rows(
        &[("A", DONOR_SPARED), ("T", DONOR_SPARED)],
        donor,
        other_order,
    );
    let filtered_own = donor_rows(
        &[("A", DONOR_SPARED), ("T", DONOR_FILTERED)],
        donor,
        other_order,
    );
    let (report, documented) = compare_default(&spared_first, &filtered_own, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    let inputs = vec!["21\t46126243\t497375\tC\tT,A\t.\t.\t.".to_string()];
    let record = |second: &str| -> Vec<Entry> {
        let text = dumped(&format!(
            r#"{{"input": "21\t46126243\t497375\tC\tT,A\t.\t.\t.", "id": "497375", "allele_string": "C/T/A", "most_severe_consequence": "splice_donor_variant", "transcript_consequences": [
                {{"variant_allele": "T", "transcript_id": "ENST00001028357", "consequence_terms": ["splice_donor_variant"], "lof": "LC", "lof_filter": "GC_TO_GT_DONOR"}},
                {{"variant_allele": "A", "transcript_id": "ENST00001028357", "consequence_terms": ["splice_donor_variant"], {second}}}]}}"#
        )) + "\n";
        parse_json(&text, &inputs).expect("the text parses")
    };
    assert_eq!(
        under_class(
            klass,
            &record(r#""lof": "LC", "lof_filter": "GC_TO_GT_DONOR""#),
            &record(r#""lof": "HC""#)
        ),
        ("green".to_string(), once)
    );
}

/// The GC_TO_GT_DONOR item stays a difference on a one-allele record, on default-format
/// rows under an identifier, on rows handed as one record's that state no allele order, on a
/// row that is not a splice-donor row, beside another differing item, and where the row's
/// LoF does not read LC beside the item and HC without it.
#[test]
fn the_first_alternate_donor_class_refuses_every_other_filter_difference() {
    let klass = "loftee_first_alt_donor_filter";
    let donor = "splice_donor_variant";
    let t_row = ("T", DONOR_FILTERED);
    let filtered = [t_row, ("A", DONOR_FILTERED)];
    let own = |second: &str| donor_rows(&[t_row, ("A", second)], donor, DONOR_RECORD);
    let pairs = [
        (
            donor_rows(&[("A", DONOR_FILTERED)], donor, "21_46126243_C/A"),
            donor_rows(&[("A", DONOR_SPARED)], donor, "21_46126243_C/A"),
        ),
        (
            donor_rows(&filtered, donor, "497375"),
            donor_rows(&[t_row, ("A", DONOR_SPARED)], donor, "497375"),
        ),
        (
            own(DONOR_FILTERED),
            donor_rows(
                &[("T", DONOR_SPARED), ("A", DONOR_FILTERED)],
                donor,
                DONOR_RECORD,
            ),
        ),
        (
            donor_rows(&filtered, "splice_acceptor_variant", DONOR_RECORD),
            donor_rows(
                &[t_row, ("A", DONOR_SPARED)],
                "splice_acceptor_variant",
                DONOR_RECORD,
            ),
        ),
        (
            own(DONOR_FILTERED),
            own(&DONOR_SPARED.replace("LoF=HC", "LoF=LC;LoF_filter=END_TRUNC")),
        ),
        (
            own(DONOR_FILTERED),
            own(&DONOR_SPARED.replace("LoF=HC", "LoF=OS")),
        ),
        (
            own(&DONOR_FILTERED.replace("LoF=LC", "LoF=HC")),
            own(&DONOR_SPARED.replace("LoF=HC", "LoF=LC")),
        ),
        (
            own(DONOR_FILTERED),
            own(&DONOR_FILTERED.replace("INTRON_SIZE:261", "INTRON_SIZE:262")),
        ),
    ];
    for (reference, candidate) in &pairs {
        let (report, documented) = compare_default(reference, candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert!(counts(&documented.applied).is_empty(), "{candidate}");
    }

    // Rows under an identifier state no allele order: where a caller vouches that they are
    // one record's, neither the first alternate's own row nor a later allele's is accepted.
    let entries_of = |text: &str| parse_default(text).expect("the text parses");
    let named = entries_of(&donor_rows(&filtered, donor, "497375"));
    let named_rows: Vec<&Entry> = named.iter().collect();
    let side = Side::default();
    let context = RuleContext {
        record_rows: Some(&named_rows),
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    let candidates = [
        [("T", DONOR_SPARED), ("A", DONOR_FILTERED)],
        [t_row, ("A", DONOR_SPARED)],
    ];
    for (index, lofs) in candidates.iter().enumerate() {
        let e = &named[index];
        let a = &entries_of(&donor_rows(lofs, donor, "497375"))[index];
        assert!(first_alternate(e, &named_rows).is_none());
        let accepts = |name: &str, reference: &str, candidate: &str| {
            field_class_accepts(klass, e, name, reference, candidate, Some(a), &context)
                .expect("the class reads the row")
        };
        assert!(!accepts("LoF_filter", "GC_TO_GT_DONOR", ABSENT), "{index}");
        assert!(!accepts("LoF", "LC", "HC"), "{index}");
    }

    // A one-allele insertion whose rows spell the allele two ways is no multi-allelic record:
    // the other spelling sits on another feature.
    let insertion = |lof: &str| {
        default_text(&[
            default_row(
                "21_26021981_-/AAG",
                "21:26021980-26021981",
                "GAA",
                "ENST1",
                donor,
                &format!("IMPACT=HIGH;STRAND=1;{lof}"),
            ),
            default_row(
                "21_26021981_-/AAG",
                "21:26021980-26021981",
                "AAG",
                "ENST2",
                "intron_variant",
                "IMPACT=MODIFIER;STRAND=1",
            ),
        ])
    };
    let (report, documented) = compare_default(
        &insertion(DONOR_FILTERED),
        &insertion(DONOR_SPARED),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert!(counts(&documented.applied).is_empty());
    // A VCF line whose alleles share their first base beside a deletion writes its entries
    // without that base: the first alternate's row is refused, and the deletion's row is a
    // later allele's, accepted.
    let inputs = vec!["21\t46126242\tshared_base\tGC\tGT,G\t.\t.\t.".to_string()];
    let line = |t_lof: &str, deletion_lof: &str| -> Vec<Entry> {
        let text = format!(
            "{}{}CSQ=T|splice_donor_variant|ENST1|Transcript|{t_lof},-|splice_donor_variant|ENST1|Transcript|{deletion_lof}\n",
            csq_header("Allele|Consequence|Feature|Feature_type|LoF|LoF_filter"),
            record_columns(&inputs[0])
        );
        parse_vcf(&text, &inputs).expect("the text parses")
    };
    let (kept, spared) = ("LC|GC_TO_GT_DONOR", "HC|");
    let both = line(kept, kept);
    assert_eq!(
        under_class(klass, &both, &line(spared, kept)),
        ("red".to_string(), count_map(&[]))
    );
    assert_eq!(
        under_class(klass, &both, &line(kept, spared)),
        ("green".to_string(), count_map(&[(klass, 1)]))
    );
    let alternates =
        |alleles: &[&str]| -> Vec<String> { alleles.iter().map(|a| a.to_string()).collect() };
    assert_eq!(vcf_first_alternate("GC", &alternates(&["GT", "G"])), "T");
    assert_eq!(vcf_first_alternate("GC", &alternates(&["G", "GT"])), "-");
    // A default row names such a record by the line's REF and ALT: the same two rows there.
    let named_by_the_line = |t_lof: &str, deletion_lof: &str| {
        lof_rows(
            "21_46126242_GC/GT/G",
            "21:46126243",
            "ENST1",
            &[("T", donor, t_lof), ("-", donor, deletion_lof)],
        )
    };
    let both = named_by_the_line(DONOR_FILTERED, DONOR_FILTERED);
    let (report, documented) = compare_default(
        &both,
        &named_by_the_line(DONOR_SPARED, DONOR_FILTERED),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert!(counts(&documented.applied).is_empty());
    let (report, documented) = compare_default(
        &both,
        &named_by_the_line(DONOR_FILTERED, DONOR_SPARED),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
    assert_eq!(vcf_first_alternate("C", &alternates(&["T", "A"])), "T");
    assert_eq!(vcf_first_alternate("AC", &alternates(&["AT", "AG"])), "AT");
    // A * of one base beside a two-base REF is an ALT of another length: the record is chopped.
    assert_eq!(
        vcf_first_alternate("GC", &alternates(&["GT", "GA", "*"])),
        "T"
    );
    assert_eq!(vcf_first_alternate("C", &alternates(&["T", "*"])), "T");
}

/// One default-format output of one row whose Extra text ends with the NearestExonJB value,
/// left out when `value` is empty.
fn nearest_exon_row(row: (&str, &str, &str, &str, &str, &str), value: &str) -> String {
    let (name, location, allele, feature, consequence, extra) = row;
    let extra = if value.is_empty() {
        extra.to_string()
    } else {
        format!("{extra};NearestExonJB={value}")
    };
    default_text(&[default_row(
        name,
        location,
        allele,
        feature,
        consequence,
        &extra,
    )])
}

/// A substitution nine bases before the first exon of a transcript with two-base introns:
/// the reference writes no NearestExonJB on the upstream row of a position-sorted file and
/// the candidate writes the first exon's item; the JSON format carries the two fields as
/// distance and nearestexonjb.
#[test]
fn the_stretch_priming_class_accepts_a_value_within_twelve_bases_outside_the_transcript() {
    let klass = "nearest_exon_jb_stretch_priming";
    let once = count_map(&[(klass, 1)]);
    let row = (
        "nine_before_first_exon",
        "21:6354255",
        "T",
        "ENST00000624847",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=9;STRAND=1",
    );
    let reference = nearest_exon_row(row, "");
    let candidate = nearest_exon_row(row, "ENSE00003757672|9|start|259");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(
        at(&report, &["columns", "NearestExonJB", "masked"]),
        &json(r#"{"nearest_exon_jb_stretch_priming": 1}"#)
    );
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let at_the_stretch = (
        row.0,
        row.1,
        row.2,
        row.3,
        row.4,
        "IMPACT=MODIFIER;DISTANCE=12;STRAND=1",
    );
    assert_eq!(
        verdict_of(
            &nearest_exon_row(at_the_stretch, ""),
            &nearest_exon_row(at_the_stretch, "ENSE00003757672|12|start|259"),
            &[klass]
        ),
        "green"
    );
    let inputs = vec!["21\t6354255\tnine_before_first_exon\tG\tT\t.\t.\t.".to_string()];
    let record = |value: &str| -> Vec<Entry> {
        let text = dumped(&format!(
            r#"{{"input": "21\t6354255\tnine_before_first_exon\tG\tT\t.\t.\t.", "id": "nine_before_first_exon", "most_severe_consequence": "upstream_gene_variant", "transcript_consequences": [
                {{"variant_allele": "T", "transcript_id": "ENST00000624847", "consequence_terms": ["upstream_gene_variant"], "distance": 9{value}}}]}}"#
        )) + "\n";
        parse_json(&text, &inputs).expect("the text parses")
    };
    assert_eq!(
        under_class(
            klass,
            &record(""),
            &record(r#", "nearestexonjb": "ENSE00003757672|9|start|259""#)
        ),
        ("green".to_string(), once)
    );
}

const SHIFT_PROBE: (&str, &str, &str, &str, &str, &str) = (
    "shift_probe",
    "21:37480826",
    "-",
    "ENST00000338785",
    "splice_donor_variant",
    "IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1",
);
const ONTO_EXON: (&str, &str, &str, &str, &str, &str) = (
    "onto_exon",
    "21:37480637",
    "-",
    "ENST00000338785",
    "frameshift_variant,splice_region_variant",
    "IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1",
);

/// Under --shift_3prime 1 --shift_length: a deletion that shifts off the last base of an
/// exon, so the reference writes no NearestExonJB where the candidate writes the exon; one
/// that shifts onto an exon's first base, so the reference alone writes it; two values
/// naming different exons are the class's too, in the VCF format's spelling as well.
#[test]
fn the_shifted_span_class_accepts_the_exon_set_read_at_the_shifted_span() {
    let klass = "nearest_exon_jb_shifted_span";
    let once = count_map(&[(klass, 1)]);
    let reference = nearest_exon_row(SHIFT_PROBE, "");
    let candidate = nearest_exon_row(SHIFT_PROBE, "ENSE00003584551|0|end|189");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(
        at(&report, &["columns", "NearestExonJB", "masked"]),
        &json(r#"{"nearest_exon_jb_shifted_span": 1}"#)
    );
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let (report, documented) = compare_default(
        &nearest_exon_row(ONTO_EXON, "ENSE00003584551|1|start|189"),
        &nearest_exon_row(ONTO_EXON, ""),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once);
    let other_exon = nearest_exon_row(
        SHIFT_PROBE,
        "ENSE00003584551|1|start|189,ENSE00003675684|0|end|100",
    );
    assert_eq!(
        verdict_of(
            &nearest_exon_row(SHIFT_PROBE, "ENSE00003584551|1|start|189"),
            &other_exon,
            &[klass]
        ),
        "green"
    );
    let inputs = vec!["21\t37480825\tshift_probe\tAG\tA\t.\t.\t.".to_string()];
    let vcf = format!(
        "{}{}CSQ=-|splice_donor_variant|ENST00000338785|Transcript|1|ENSE00003584551+1+start+189\n",
        csq_header("Allele|Consequence|Feature|Feature_type|SHIFT_LENGTH|NearestExonJB"),
        record_columns(&inputs[0])
    );
    let vcf_candidate = vcf.replace(
        "ENSE00003584551+1+start+189",
        "ENSE00003584551+1+start+189&ENSE00003675684+0+end+100",
    );
    let vcf_entries = |text: &str| parse_vcf(text, &inputs).expect("the text parses");
    assert_eq!(
        under_class(klass, &vcf_entries(&vcf), &vcf_entries(&vcf_candidate)),
        ("green".to_string(), once)
    );
    let item = |exon: &str, distance: &str, end: &str, length: &str| {
        (
            exon.to_string(),
            distance.to_string(),
            end.to_string(),
            length.to_string(),
        )
    };
    assert_eq!(
        nearest_exon_items("ENSE00003584551+1+start+189&ENSE00003675684+0+end+100"),
        Some(vec![
            item("ENSE00003584551", "1", "start", "189"),
            item("ENSE00003675684", "0", "end", "100")
        ])
    );
    assert_eq!(
        nearest_exon_items("|106|end|135,|2780|start|109"),
        Some(vec![
            item("", "106", "end", "135"),
            item("", "2780", "start", "109")
        ])
    );
    assert_eq!(
        nearest_exon_items("ENSE00003595154|87|start_end|175"),
        Some(vec![item("ENSE00003595154", "87", "start_end", "175")])
    );
    assert!(nearest_exon_items("ENSE00003584551|1|start").is_none());
    assert!(nearest_exon_items("-").is_none());
}

const UTR_GAIN: &str = "5_prime_UTR_premature_start_codon_gain_variant";
const UTR_LOSS: &str = "5_prime_UTR_premature_start_codon_loss_variant";
const UTR_INPUTS: [&str; 3] = [
    "21\t25734954\tutr_ins\tA\tACAT\t.\t.\t.",
    "21\t500\tmnv\tAC\tGT\t.\t.\t.",
    "21\t600\tsnv\tA\tT\t.\t.\t.",
];

/// A 5' UTR row of the default format with its cDNA_position and no UTRAnnotator value.
fn utr_row(name: &str, location: &str, allele: &str, feature: &str, cdna: &str) -> String {
    let extra = "IMPACT=MODIFIER;STRAND=1;5UTR_consequence=-;5UTR_annotation=-";
    default_row(
        name,
        location,
        allele,
        feature,
        "5_prime_UTR_variant",
        extra,
    )
    .replacen(
        "\t-\t-\t-\t-\t-\t-\t",
        &format!("\t{cdna}\t-\t-\t-\t-\t-\t"),
        1,
    )
}

/// The VCF and JSON entries carry no Location: the record's own alleles say the
/// insertion's row is one, and the same value pair on a substitution's record stays a
/// difference. The JSON format writes an insertion's cdna_start past its cdna_end and the
/// annotation as an object.
#[test]
fn the_gained_start_window_class_reads_the_record_in_the_vcf_and_json_formats() {
    let klass = "utrannotator_gained_start_window";
    let once = count_map(&[(klass, 1)]);
    let inputs = vec![
        UTR_INPUTS[0].to_string(),
        "21\t400\tsnv\tA\tG\t.\t.\t.".to_string(),
    ];
    let lost = format!("{UTR_LOSS}|CapDistanceToStart=63:type=uORF");
    let vcf_ref = format!(
        "{}{}CSQ=CAT|5_prime_UTR_variant|ENST00000000033|Transcript|64-65|{lost}\n{}CSQ=G|5_prime_UTR_variant|ENST3|Transcript|64|{lost}\n",
        csq_header(
            "Allele|Consequence|Feature|Feature_type|cDNA_position|5UTR_consequence|5UTR_annotation"
        ),
        record_columns(&inputs[0]),
        record_columns(&inputs[1])
    );
    let vcf_cand = vcf_ref.replace(
        &lost,
        &format!("{UTR_GAIN}|CapDistanceToStart=64:type=uORF"),
    );
    let vcf_entries = |text: &str| parse_vcf(text, &inputs).expect("the text parses");
    let mut documented = documented_with(&[klass]);
    let report = compare(
        &vcf_entries(&vcf_ref),
        &vcf_entries(&vcf_cand),
        &mut documented,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(
        number(&report, &["columns", "5UTR_consequence", "differing"]),
        1
    );

    let object = |i: usize,
                  allele_string: &str,
                  allele: &str,
                  feature: &str,
                  start: i64,
                  end: i64,
                  utr: &str,
                  at: i64|
     -> String {
        let input = inputs[i].replace('\t', "\\t");
        dumped(&format!(
            r#"{{"input": "{input}", "id": "r{i}", "allele_string": "{allele_string}", "most_severe_consequence": "5_prime_UTR_variant", "transcript_consequences": [
                {{"variant_allele": "{allele}", "transcript_id": "{feature}", "consequence_terms": ["5_prime_UTR_variant"], "cdna_start": {start}, "cdna_end": {end},
                 "5utr_consequence": "{utr}", "5utr_annotation": {{"1": {{"CapDistanceToStart": {at}, "type": "uORF"}}}}}}]}}"#
        )) + "\n"
    };
    let records = |utr: &str, at: i64| -> Vec<Entry> {
        let text = object(0, "-/CAT", "CAT", "ENST00000000033", 65, 64, utr, at)
            + &object(1, "A/G", "G", "ENST3", 64, 64, utr, at);
        parse_json(&text, &inputs).expect("the text parses")
    };
    let mut documented = documented_with(&[klass]);
    let report = compare(
        &records(UTR_LOSS, 63),
        &records(UTR_GAIN, 64),
        &mut documented,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(counts(&documented.applied), once);
    assert_eq!(
        number(&report, &["columns", "5utr_consequence", "differing"]),
        1
    );
    // The lengths Ensembl VEP hands the plugin: a record's shared first base is dropped, a
    // deletion's Location and a Location of one base or of three or more state the reference
    // allele, and a two-base Location under an identifier states none.
    let entries = parse_default(&default_text(&[
        utr_row("21_100_AC/ACGT", "21:101", "CGT", "ENST1", "7"),
        utr_row("del", "21:1000-1004", "-", "ENST2", "336-340"),
        utr_row(
            "utr_ins",
            "21:25734954-25734955",
            "CAT",
            "ENST00000000033",
            "64-65",
        ),
        utr_row("rs1", "21:700-702", "A", "ENST1", "7-9"),
    ]))
    .expect("the rows parse");
    let lengths: Vec<Option<(i64, i64)>> = entries
        .iter()
        .map(|e| record_allele_lengths(e, None))
        .collect();
    assert_eq!(lengths, [Some((1, 3)), Some((5, 0)), None, Some((3, 1))]);
    let input_lines: Vec<String> = UTR_INPUTS.iter().map(|line| line.to_string()).collect();
    let mut rows = RecordRows::new(&entries, Some(&input_lines));
    let record = rows
        .of(&entries[2])
        .expect("the rows are read")
        .expect("the identifier names one line");
    assert_eq!(
        record_allele_lengths(&entries[2], Some(&record)),
        Some((0, 3))
    );
    let alleles = |spelt: &[&str]| -> Vec<String> { spelt.iter().map(|a| a.to_string()).collect() };
    assert_eq!(
        vcf_alleles_as_spelt("A", &alleles(&["AT"])),
        ("-".to_string(), alleles(&["T"]))
    );
    assert_eq!(
        vcf_alleles_as_spelt("AC", &alleles(&["AT", "AG"])),
        ("AC".to_string(), alleles(&["AT", "AG"]))
    );
    let spans: Vec<Option<(i64, i64)>> = entries.iter().map(cdna_span).collect();
    assert_eq!(
        spans,
        [Some((7, 7)), Some((336, 340)), Some((64, 65)), Some((7, 9))]
    );
}

/// Default-format rows of one record of two alternates on one transcript, its alleles spelt
/// by the Uploaded_variation `name`: (allele, Amino_acids, REVEL) each.
fn revel_rows(rows: &[(&str, &str, &str)], name: &str) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, amino_acids, score)| {
            format!(
                "{name}\t21:13610388\t{allele}\tENSG1\tENST00000299443\tTranscript\tmissense_variant\t-\t-\t-\t{amino_acids}\t-\t-\tIMPACT=MODERATE;STRAND=1;REVEL={score}"
            )
        })
        .collect();
    default_text(&rows)
}

const SPLICEAI_VALUE: &str = "SpliceAI_pred=GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1";

/// Default-format rows of one record whose alleles share their first base, on one
/// transcript: (allele, the SpliceAI pair of the Extra text or nothing) each.
fn chopped_rows(name: &str, rows: &[(&str, &str)]) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, spliceai)| {
            let extra = if spliceai.is_empty() {
                "IMPACT=MODIFIER;STRAND=1".to_string()
            } else {
                format!("IMPACT=MODIFIER;STRAND=1;{spliceai}")
            };
            default_row(
                name,
                "21:26965500",
                allele,
                "ENST00000284987",
                "intron_variant",
                &extra,
            )
        })
        .collect();
    default_text(&rows)
}

/// A default row under an identifier does not say whether its Location's rows are one
/// record's; the input VCF does. With its lines a class that reads a record's alleles applies
/// to a row whose identifier one line alone carries, and still refuses one two lines share;
/// the donor class reads the record's first alternate from that line.
#[test]
fn the_allele_reading_classes_read_a_named_records_alleles_from_the_input() {
    let klass = "revel_sibling_allele_rows";
    let reference = revel_rows(
        &[("C", "M/L", "0.025"), ("T", "M/L", "0.025")],
        "two_leucines",
    );
    let candidate = revel_rows(
        &[("C", "M/L", "0.025"), ("T", "M/L", "0.027")],
        "two_leucines",
    );
    // The verdict and what was applied of two texts under one class with the input's lines.
    let compared = |lines: Option<&[&str]>, expected: &str, actual: &str, name: &str| {
        let lines: Option<Vec<String>> =
            lines.map(|lines| lines.iter().map(|line| line.to_string()).collect());
        let mut documented = documented_with(&[name]);
        let side = Side::default();
        let options = CompareOptions {
            input_lines: lines.as_deref(),
            ..CompareOptions::new(&side)
        };
        let report = compare_entries(
            &parse_default(expected).expect("the reference parses"),
            &parse_default(actual).expect("the candidate parses"),
            &mut documented,
            &options,
        )
        .expect("the entries compare");
        (verdict(&report).to_string(), counts(&documented.applied))
    };
    let green = |name: &str| ("green".to_string(), count_map(&[(name, 1)]));
    let red = ("red".to_string(), count_map(&[]));
    let one_record = ["21\t13610388\ttwo_leucines\tA\tC,T\t.\t.\t."];
    let split_pair = [
        "21\t13610388\ttwo_leucines\tA\tC\t.\t.\t.",
        "21\t13610388\ttwo_leucines\tA\tT\t.\t.\t.",
    ];
    assert_eq!(
        compared(Some(&one_record), &reference, &candidate, klass),
        green(klass)
    );
    assert_eq!(compared(None, &reference, &candidate, klass), red);
    assert_eq!(
        compared(Some(&split_pair), &reference, &candidate, klass),
        red
    );
    assert_eq!(
        compared(
            Some(&["21\t13610388\tanother\tA\tC,T\t.\t.\t."]),
            &reference,
            &candidate,
            klass
        ),
        red
    );
    let donor = "loftee_first_alt_donor_filter";
    let term = "splice_donor_variant";
    let line = ["21\t46126243\t497375\tC\tT,A\t.\t.\t."];
    let both = donor_rows(
        &[("T", DONOR_FILTERED), ("A", DONOR_FILTERED)],
        term,
        "497375",
    );
    assert_eq!(
        compared(
            Some(&line),
            &both,
            &donor_rows(
                &[("T", DONOR_FILTERED), ("A", DONOR_SPARED)],
                term,
                "497375"
            ),
            donor
        ),
        green(donor)
    );
    assert_eq!(
        compared(
            Some(&line),
            &both,
            &donor_rows(
                &[("T", DONOR_SPARED), ("A", DONOR_FILTERED)],
                term,
                "497375"
            ),
            donor
        ),
        red
    );
    let entries = parse_default(&reference).expect("the reference parses");
    let input_lines: Vec<String> = one_record.iter().map(|line| line.to_string()).collect();
    let mut rows = RecordRows::new(&entries, Some(&input_lines));
    let group = rows
        .of(&entries[1])
        .expect("the rows are read")
        .expect("the identifier names one line");
    let kinds: Vec<String> = group.iter().map(|row| entry_kind(row)).collect();
    assert_eq!(kinds, ["Transcript", "Transcript", "record"]);
    assert_eq!(group[2].fields.get("ALT").map(String::as_str), Some("C,T"));
    // The line's REF and ALT say what a row under an identifier cannot: beside a deletion of
    // the second base, an allele of two bases at the one-base Location is an insertion's and
    // an allele of one base a substitution's.
    let spliceai = "spliceai_chopped_multiallelic";
    for (alternates, allele, outcome) in
        [("G,GA", "A", green(spliceai)), ("G,GCA", "CA", red.clone())]
    {
        let chopped = format!("21\t26965499\tchopped_multiallelic\tGC\t{alternates}\t.\t.\t.");
        let unscored = chopped_rows("chopped_multiallelic", &[("-", ""), (allele, "")]);
        let scored = chopped_rows(
            "chopped_multiallelic",
            &[("-", ""), (allele, SPLICEAI_VALUE)],
        );
        assert_eq!(
            compared(Some(&[chopped.as_str()]), &unscored, &scored, spliceai),
            outcome,
            "{alternates}"
        );
    }
}
