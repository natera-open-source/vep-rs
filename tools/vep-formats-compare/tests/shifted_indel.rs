// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The written-span gates of a shifted indel and the intron offsets abutting an insertion,
//! on planted rows.

mod compare_common;
mod shifted_edited_common;

use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::consts::HGVS_INTRON_OFFSET_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::shifted_indel::coordinate_intron_offset;
use vep_formats_compare::rules::shifted_indel::intron_offset_abutting_insertion;
use vep_formats_compare::shared::ABSENT;

const GATE: &str = "hgvs_shifted_indel_written_span_gate";
const PROTEIN_GATE: &str = "hgvsp_shifted_indel_written_span_gate";
const ABUTTING: &str = "hgvs_intron_offset_abutting_insertion";

/// What `klass` says of the field `field` on a pair of default-format rows.
fn offset_grid_verdict(klass: &str, reference: &str, candidate: &str, field: &str) -> bool {
    let (e, a) = (default_entry(reference), default_entry(candidate));
    class_accepts(
        klass,
        &e,
        field,
        &field_of(Some(&e), &[field]),
        &field_of(Some(&a), &[field]),
        Some(&a),
    )
}

/// An insertion's or deletion's row Ensembl VEP leaves without HGVSc and HGVS_OFFSET, where
/// the candidate writes an indel notation of the row's transcript and a non-zero offset, is the
/// class, counted once per row in every format, with or without a SHIFT_LENGTH field and with
/// the HGVSIntronOffset plugin's two fields the notation brings, and on a two-base insertion at a
/// two-base Location, which the candidate's dup or ins notation marks an insertion's row; a row
/// the reference gives a notation or an offset, a candidate without the offset or with a zero
/// one, a substitution's notation, another transcript's notation, a substitution's row, a row
/// annotated upstream or downstream and a zero SHIFT_LENGTH stay differences.
#[test]
fn the_written_span_gate_accepts_the_notation_and_offset_the_reference_lacks() {
    let class = [GATE];
    let shifted = "IMPACT=MODIFIER;STRAND=1;SHIFT_LENGTH=17";
    let row = default_row(
        "21:16420341:C:CT",
        "21:16420357-16420358",
        "T",
        "ENST00000602620",
        "non_coding_transcript_exon_variant",
        shifted,
    );
    let reference =
        default_text(&[row.replacen("\t-\t-\t-\t-\t-\t-\t", "\t15-16\t-\t-\t-\t-\t-\t", 1)]);
    // `text` with the notation and the offset after its STRAND, each left out where empty.
    let written = |hgvsc: &str, offset: &str, text: &str| -> String {
        let mut fields = String::new();
        if !hgvsc.is_empty() {
            fields.push_str(&format!(";HGVSc={hgvsc}"));
        }
        if !offset.is_empty() {
            fields.push_str(&format!(";HGVS_OFFSET={offset}"));
        }
        text.replace("STRAND=1", &format!("STRAND=1{fields}"))
    };
    let candidate = |hgvsc: &str, offset: &str| written(hgvsc, offset, &reference);
    let duplication = "ENST00000602620.2:n.15dup";

    let (report, documented) = compare_default(&reference, &candidate(duplication, "17"), &class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(GATE));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &masked_once(GATE)
    );
    assert_eq!(
        verdict_of(&reference, &candidate(duplication, "17"), &[]),
        "red"
    );
    // A deletion carried in, one still straddling the first base and one carried over the
    // last base.
    let deletion = reference
        .replace("\tT\t", "\t-\t")
        .replace("21:16420357-16420358", "21:13038013-13038033")
        .replace("ENST00000602620", "ENST00000447861");
    for notation in [
        "ENST00000447861.6:n.5_25del",
        "ENST00000447861.6:n.1_16del",
        "ENST00000447861.6:n.348del",
        "ENST00000447861.6:c.*2333_*2540del",
    ] {
        assert_eq!(
            verdict_of(&deletion, &written(notation, "33", &deletion), &class),
            "green",
            "{notation}"
        );
    }
    // A two-base insertion at a two-base Location reads as a substitution from the row's
    // text; the candidate's duplication or insertion notation marks the row an insertion's,
    // a substitution's notation does not.
    let two_bases = reference.replace("\tT\t", "\tTT\t");
    for notation in [
        "ENST00000602620.2:c.-49_-48dup",
        "ENST00000602620.2:n.1_2insTG",
    ] {
        assert_eq!(
            verdict_of(&two_bases, &written(notation, "17", &two_bases), &class),
            "green",
            "{notation}"
        );
    }
    assert_eq!(
        verdict_of(
            &two_bases,
            &written("ENST00000602620.2:n.16A>T", "17", &two_bases),
            &class
        ),
        "red"
    );
    let one_base = reference.replace("21:16420357-16420358", "21:16420358");
    assert_eq!(
        verdict_of(&one_base, &written(duplication, "17", &one_base), &class),
        "red"
    );
    // Under --hgvs alone no SHIFT_LENGTH is printed; the offset alone says the row shifted.
    let unlabelled = deletion.replace(";SHIFT_LENGTH=17", "");
    let carried_over = written("ENST00000447861.6:n.348del", "12", &unlabelled);
    assert_eq!(verdict_of(&unlabelled, &carried_over, &class), "green");
    // The HGVSIntronOffset plugin writes its two offsets beside the notation.
    let with_plugin = carried_over.replace(
        "HGVS_OFFSET=12",
        "HGVS_OFFSET=12;HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0",
    );
    assert_eq!(
        default_under(GATE, &unlabelled, &with_plugin),
        ("green".to_string(), once(GATE))
    );
    for (hgvsc, offset) in [
        (duplication, ""),
        (duplication, "0"),
        ("", "17"),
        ("ENST00000602620.2:n.16A>T", "17"),
        ("ENST00000698005.1:n.17dup", "17"),
    ] {
        assert_eq!(
            verdict_of(&reference, &candidate(hgvsc, offset), &class),
            "red",
            "{hgvsc} {offset}"
        );
    }
    assert_eq!(
        verdict_of(
            &candidate(duplication, "17"),
            &candidate(duplication, "17").replace("n.15dup", "n.14dup"),
            &class
        ),
        "red"
    );
    assert_eq!(
        verdict_of(&candidate("", "17"), &candidate(duplication, "17"), &class),
        "red"
    );
    let unshifted = reference.replace("SHIFT_LENGTH=17", "SHIFT_LENGTH=0");
    assert_eq!(
        verdict_of(&unshifted, &written(duplication, "17", &unshifted), &class),
        "red"
    );
    let substitution = reference.replace("21:16420357-16420358", "21:16420358");
    assert_eq!(
        verdict_of(
            &substitution,
            &written(duplication, "17", &substitution),
            &class
        ),
        "red"
    );
    // A row annotated outside the transcript has no notation on either engine.
    for term in ["upstream_gene_variant", "downstream_gene_variant"] {
        let outside = reference.replace("non_coding_transcript_exon_variant", term);
        assert_eq!(
            verdict_of(&outside, &written(duplication, "17", &outside), &class),
            "red",
            "{term}"
        );
    }
    // The tab, VCF and JSON formats carry the fields apart and count the row once.
    let header = tab_header(&["IMPACT", "STRAND", "HGVSc", "HGVS_OFFSET", "SHIFT_LENGTH"]);
    let tab_row = "21:16420341:C:CT\t21:16420357-16420358\tT\tENSG1\tENST00000602620\tTranscript\tnon_coding_transcript_exon_variant\t15-16\t-\t-\t-\t-\t-\tMODIFIER\t1\t";
    let tab = |tail: &str| tab_entries(&format!("{header}\n{tab_row}{tail}\n"));
    assert_eq!(
        under_class(
            GATE,
            &tab("-\t-\t17"),
            &tab("ENST00000602620.2:n.15dup\t17\t17")
        ),
        ("green".to_string(), once(GATE))
    );
    assert_eq!(
        under_class(
            GATE,
            &tab("-\t-\t17"),
            &tab("ENST00000602620.2:n.15dup\t-\t17")
        )
        .0,
        "red"
    );
    let inputs = vec!["21\t16420340\t21:16420341:C:CT\tC\tCT\t.\t.\t.".to_string()];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|HGVSc|HGVS_OFFSET|SHIFT_LENGTH",
    );
    let csq_row = format!(
        "{}CSQ=T|non_coding_transcript_exon_variant|ENST00000602620|Transcript|15-16|",
        record_columns(&inputs[0])
    );
    let vcf = |tail: &str| vcf_entries(&format!("{header}{csq_row}{tail}\n"), &inputs);
    assert_eq!(
        under_class(GATE, &vcf("||17"), &vcf("ENST00000602620.2:n.15dup|17|17")),
        ("green".to_string(), once(GATE))
    );
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("21:16420341:C:CT")),
        ("allele_string", text("-/T")),
        (
            "most_severe_consequence",
            text("non_coding_transcript_exon_variant"),
        ),
    ];
    let consequence = object(&[
        ("variant_allele", text("T")),
        ("transcript_id", text("ENST00000602620")),
        (
            "consequence_terms",
            Value::List(vec![text("non_coding_transcript_exon_variant")]),
        ),
        ("cdna_start", Value::Int(15)),
        ("cdna_end", Value::Int(16)),
        ("shift_length", Value::Int(17)),
    ]);
    let object_of = |fields: &[(&str, Value)]| {
        json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        )
    };
    let notation = ("hgvsc", text("ENST00000602620.2:n.15dup"));
    assert_eq!(
        under_class(
            GATE,
            &object_of(&[]),
            &object_of(&[notation.clone(), ("hgvs_offset", Value::Int(17))])
        ),
        ("green".to_string(), once(GATE))
    );
    assert_eq!(
        under_class(GATE, &object_of(&[]), &object_of(&[notation])).0,
        "red"
    );
}

const ABUTTING_EXTRA: &str = "IMPACT=MODIFIER;DISTANCE=0;STRAND=-1";

/// The plugin's two fields, the start's then the end's.
fn offset_fields() -> [&'static str; 2] {
    [HGVS_INTRON_OFFSET_FIELDS[0], HGVS_INTRON_OFFSET_FIELDS[1]]
}

/// What a row of an insertion directly beside a transcript prints.
#[derive(Clone, Copy)]
struct Abutting<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    feature: &'a str,
    term: &'a str,
    extra: &'a str,
}

impl Default for Abutting<'_> {
    fn default() -> Self {
        Abutting {
            name: "ins1",
            location: "2:1000443-1000444",
            allele: "T",
            feature: "ENST00000000021",
            term: "upstream_gene_variant",
            extra: ABUTTING_EXTRA,
        }
    }
}

/// The default-format entry of an insertion directly beside a transcript, with the
/// HGVSIntronOffset plugin's two fields where `start` and `end` give them.
fn abutting_row(start: Option<&str>, end: Option<&str>, row: Abutting) -> Entry {
    let mut extra = row.extra.to_string();
    for (field, value) in [
        ("HGVS_IntronEndOffset", end),
        ("HGVS_IntronStartOffset", start),
    ] {
        if let Some(value) = value {
            extra.push_str(&format!(";{field}={value}"));
        }
    }
    default_entry(&default_row(
        row.name,
        row.location,
        row.allele,
        row.feature,
        row.term,
        &extra,
    ))
}

/// What `intron_offset_abutting_insertion` says of the start and the end offset of the
/// reference entry `e` against the candidate entry `a`.
fn abutting_verdicts(e: &Entry, a: Option<&Entry>) -> [bool; 2] {
    offset_fields().map(|field| {
        intron_offset_abutting_insertion(
            e,
            a,
            field,
            &field_of(Some(e), &[field]),
            &field_of(a, &[field]),
        )
    })
}

/// The verdicts on a reference row carrying `start` and `end` against the same row without
/// them.
fn abutting_pair(start: Option<&str>, end: Option<&str>, row: Abutting) -> [bool; 2] {
    abutting_verdicts(
        &abutting_row(start, end, row),
        Some(&abutting_row(None, None, row)),
    )
}

const NEITHER: [bool; 2] = [false, false];
const BOTH: [bool; 2] = [true, true];

/// Four made-up insertions directly beside a transcript's first or last base, one per side
/// and strand: Ensembl VEP writes the plugin's two fields on the row, 0 and 0, or 1 and 0
/// past a reverse-strand transcript's last base, and the candidate writes neither. Both
/// fields are accepted, with or without a STRAND on the row.
#[test]
fn the_abutting_insertion_class_accepts_an_insertion_beside_a_transcript() {
    for (name, location, allele, feature, term, strand, start) in [
        (
            "ins1",
            "2:1000443-1000444",
            "T",
            "ENST00000000021",
            "upstream_gene_variant",
            "-1",
            "0",
        ),
        (
            "ins2",
            "1:2000434-2000435",
            "G",
            "ENST00000000022",
            "upstream_gene_variant",
            "1",
            "0",
        ),
        (
            "ins3",
            "2:3000890-3000891",
            "TTCT",
            "ENST00000000023",
            "downstream_gene_variant",
            "-1",
            "1",
        ),
        (
            "ins4",
            "3:4000590-4000591",
            "T",
            "ENST00000000024",
            "downstream_gene_variant",
            "1",
            "0",
        ),
    ] {
        let stranded = format!("IMPACT=MODIFIER;DISTANCE=0;STRAND={strand}");
        let row = Abutting {
            name,
            location,
            allele,
            feature,
            term,
            extra: &stranded,
        };
        assert_eq!(abutting_pair(Some(start), Some("0"), row), BOTH, "{name}");
        let unstranded = Abutting {
            extra: "IMPACT=MODIFIER;DISTANCE=0",
            ..row
        };
        assert_eq!(
            abutting_pair(Some(start), Some("0"), unstranded),
            BOTH,
            "{name}"
        );
    }
}

/// The two fields in the reference alone on a row with any term inside the transcript stay
/// differences: with the abutting row's own Extra, without a DISTANCE, and beside the HGVSc
/// such a row carries.
#[test]
fn the_abutting_insertion_class_refuses_a_row_inside_a_transcript() {
    for term in [
        "intron_variant",
        "non_coding_transcript_exon_variant",
        "splice_region_variant,intron_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "frameshift_variant",
        "upstream_gene_variant,intron_variant",
        "downstream_gene_variant,3_prime_UTR_variant",
    ] {
        let row = Abutting {
            term,
            ..Abutting::default()
        };
        assert_eq!(abutting_pair(Some("0"), Some("0"), row), NEITHER, "{term}");
        let inside = "IMPACT=MODIFIER;STRAND=-1";
        let without_distance = Abutting {
            extra: inside,
            ..row
        };
        assert_eq!(
            abutting_pair(Some("0"), Some("0"), without_distance),
            NEITHER,
            "{term}"
        );
        let notation = format!("{inside};HGVSc=ENST00000000021.1:n.1_2insA");
        let with_notation = Abutting {
            extra: &notation,
            ..row
        };
        assert_eq!(
            abutting_pair(Some("0"), Some("0"), with_notation),
            NEITHER,
            "{term}"
        );
    }
}

/// A value the candidate writes in either field, a reference pair other than the two the
/// clamped notation gives, a lone reference field and a non-integer stay differences, as
/// does 1 and 0 anywhere but past a reverse-strand transcript's last base.
#[test]
fn the_abutting_insertion_class_leaves_a_planted_value_a_difference() {
    let plain = Abutting::default();
    let reference = abutting_row(Some("0"), Some("0"), plain);
    for (start, end) in [
        (Some("5"), Some("0")),
        (Some("0"), Some("0")),
        (Some("0"), None),
        (None, Some("0")),
        (Some("-"), Some("-3")),
    ] {
        assert_eq!(
            abutting_verdicts(&reference, Some(&abutting_row(start, end, plain))),
            NEITHER,
            "{start:?} {end:?}"
        );
    }
    for (start, end) in [
        (Some("2"), Some("0")),
        (Some("0"), Some("1")),
        (Some("-1"), Some("0")),
        (Some("1"), Some("1")),
        (Some("1"), Some("0")),
        (Some("x"), Some("0")),
        (Some("0"), Some("-")),
        (Some("0"), None),
        (None, Some("0")),
    ] {
        assert_eq!(
            abutting_pair(start, end, plain),
            NEITHER,
            "{start:?} {end:?}"
        );
    }
    let downstream = Abutting {
        term: "downstream_gene_variant",
        ..plain
    };
    assert_eq!(abutting_pair(Some("1"), Some("0"), downstream), BOTH);
    for (start, end) in [("2", "0"), ("1", "1"), ("0", "1"), ("-1", "0")] {
        assert_eq!(
            abutting_pair(Some(start), Some(end), downstream),
            NEITHER,
            "{start} {end}"
        );
    }
    let forward = Abutting {
        extra: "IMPACT=MODIFIER;DISTANCE=0;STRAND=1",
        ..downstream
    };
    assert_eq!(abutting_pair(Some("1"), Some("0"), forward), NEITHER);
    assert_eq!(abutting_pair(Some("0"), Some("0"), forward), BOTH);
}

/// The same two fields on the candidate's row alone, the reference writing neither, stay
/// differences, and so does a row without its candidate.
#[test]
fn the_abutting_insertion_class_refuses_the_fields_in_the_candidate_alone() {
    let plain = Abutting::default();
    assert_eq!(
        abutting_verdicts(
            &abutting_row(None, None, plain),
            Some(&abutting_row(Some("0"), Some("0"), plain))
        ),
        NEITHER
    );
    let downstream = Abutting {
        term: "downstream_gene_variant",
        ..plain
    };
    assert_eq!(
        abutting_verdicts(
            &abutting_row(None, None, downstream),
            Some(&abutting_row(Some("1"), Some("0"), downstream))
        ),
        NEITHER
    );
    assert_eq!(
        abutting_verdicts(&abutting_row(Some("0"), Some("0"), plain), None),
        NEITHER
    );
}

/// Away from the transcript's edge (a DISTANCE other than 0 or none), on a deletion, on a
/// substitution, beside an HGVSc or an HGVS_OFFSET on either side, where the candidate's
/// term differs, and for a field that is not the plugin's, the fields stay differences.
#[test]
fn the_abutting_insertion_class_refuses_a_row_that_is_not_an_abutting_insertion() {
    let plain = Abutting::default();
    let filled = abutting_row(Some("0"), Some("0"), plain);
    for extra in [
        "IMPACT=MODIFIER;DISTANCE=1;STRAND=-1",
        "IMPACT=MODIFIER;DISTANCE=10;STRAND=-1",
        "IMPACT=MODIFIER;STRAND=-1",
    ] {
        let away = Abutting { extra, ..plain };
        assert_eq!(
            abutting_pair(Some("0"), Some("0"), away),
            NEITHER,
            "{extra}"
        );
        assert_eq!(
            abutting_verdicts(&filled, Some(&abutting_row(None, None, away))),
            NEITHER,
            "{extra}"
        );
    }
    for shape in [
        Abutting {
            allele: "-",
            ..plain
        },
        Abutting {
            allele: "-",
            location: "2:1000444-1000446",
            ..plain
        },
        Abutting {
            location: "2:1000443",
            ..plain
        },
        Abutting {
            location: "2:1000443-1000445",
            ..plain
        },
        Abutting {
            name: "2_1000443_GA/TT",
            allele: "TT",
            ..plain
        },
        Abutting {
            allele: "<INS>",
            ..plain
        },
    ] {
        assert_eq!(
            abutting_pair(Some("0"), Some("0"), shape),
            NEITHER,
            "{} {} {}",
            shape.name,
            shape.location,
            shape.allele
        );
    }
    let spelt = Abutting {
        name: "2_1000444_-/TT",
        allele: "TT",
        ..plain
    };
    assert_eq!(abutting_pair(Some("0"), Some("0"), spelt), BOTH);
    for pair in ["HGVSc=ENST00000000021.1:n.1G>T", "HGVS_OFFSET=2"] {
        let extra = format!("{ABUTTING_EXTRA};{pair}");
        let with_pair = Abutting {
            extra: &extra,
            ..plain
        };
        assert_eq!(
            abutting_verdicts(
                &abutting_row(Some("0"), Some("0"), with_pair),
                Some(&abutting_row(None, None, plain))
            ),
            NEITHER,
            "{pair}"
        );
        assert_eq!(
            abutting_verdicts(&filled, Some(&abutting_row(None, None, with_pair))),
            NEITHER,
            "{pair}"
        );
    }
    for term in ["downstream_gene_variant", "intron_variant"] {
        let other = Abutting { term, ..plain };
        assert_eq!(
            abutting_verdicts(&filled, Some(&abutting_row(None, None, other))),
            NEITHER,
            "{term}"
        );
    }
    let candidate = abutting_row(None, None, plain);
    for other in [
        "HGVSc",
        "HGVS_OFFSET",
        "DISTANCE",
        "HGVS_IntronStartOffsetX",
    ] {
        assert!(
            !intron_offset_abutting_insertion(&filled, Some(&candidate), other, "0", ABSENT),
            "{other}"
        );
    }
    assert!(!intron_offset_abutting_insertion(
        &filled,
        Some(&candidate),
        "HGVS_IntronStartOffset",
        "1",
        ABSENT
    ));
}

/// The tab, VCF and JSON formats carry the row's fields apart, the JSON under lower-case
/// keys; the same row is accepted in each, and a term inside the transcript refused.
#[test]
fn the_abutting_insertion_class_reads_the_tab_vcf_and_json_formats() {
    let header = tab_header(&[
        "IMPACT",
        "DISTANCE",
        "STRAND",
        "HGVSc",
        "HGVS_OFFSET",
        "HGVS_IntronStartOffset",
        "HGVS_IntronEndOffset",
    ]);
    let tab_row = "ins1\t2:1000443-1000444\tT\tENSG1\tENST00000000021\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t0\t-1\t-\t-\t";
    let tab =
        |row: &str, offsets: &str| tab_entries(&format!("{header}\n{row}{offsets}\n")).remove(0);
    assert_eq!(
        abutting_verdicts(&tab(tab_row, "0\t0"), Some(&tab(tab_row, "-\t-"))),
        BOTH
    );
    assert_eq!(
        abutting_verdicts(&tab(tab_row, "-\t-"), Some(&tab(tab_row, "0\t0"))),
        NEITHER
    );
    let inside = tab_row.replace("upstream_gene_variant", "intron_variant");
    assert_eq!(
        abutting_verdicts(&tab(&inside, "0\t0"), Some(&tab(&inside, "-\t-"))),
        NEITHER
    );

    let inputs = vec!["2\t1000443\tins1\tG\tGT\t.\t.\t.".to_string()];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|DISTANCE|STRAND|HGVSc|HGVS_OFFSET|HGVS_IntronStartOffset|HGVS_IntronEndOffset",
    );
    let of_transcript = |entries: Vec<Entry>| transcript_entry(entries, "ENST00000000021");
    let csq = |term: &str, offsets: &str| {
        let line = format!(
            "{}CSQ=T|{term}|ENST00000000021|Transcript|0|-1|||{offsets}\n",
            record_columns(&inputs[0])
        );
        of_transcript(vcf_entries(&format!("{header}{line}"), &inputs))
    };
    assert_eq!(
        abutting_verdicts(
            &csq("upstream_gene_variant", "0|0"),
            Some(&csq("upstream_gene_variant", "|"))
        ),
        BOTH
    );
    assert_eq!(
        abutting_verdicts(
            &csq("upstream_gene_variant", "|"),
            Some(&csq("upstream_gene_variant", "0|0"))
        ),
        NEITHER
    );
    assert_eq!(
        abutting_verdicts(
            &csq("intron_variant", "0|0"),
            Some(&csq("intron_variant", "|"))
        ),
        NEITHER
    );

    let object_of = |term: &str, fields: &[(&str, Value)]| {
        let consequence = object(&[
            ("variant_allele", text("T")),
            ("transcript_id", text("ENST00000000021")),
            ("consequence_terms", Value::List(vec![text(term)])),
            ("distance", Value::Int(0)),
            ("strand", Value::Int(-1)),
        ]);
        let record = [
            ("input", text(&inputs[0])),
            ("id", text("ins1")),
            ("allele_string", text("-/T")),
            ("most_severe_consequence", text(term)),
        ];
        of_transcript(json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        ))
    };
    let offsets = [
        ("hgvs_intronstartoffset", Value::Int(0)),
        ("hgvs_intronendoffset", Value::Int(0)),
    ];
    let upstream = "upstream_gene_variant";
    assert_eq!(
        abutting_verdicts(
            &object_of(upstream, &offsets),
            Some(&object_of(upstream, &[]))
        ),
        BOTH
    );
    assert_eq!(
        abutting_verdicts(
            &object_of(upstream, &[]),
            Some(&object_of(upstream, &offsets))
        ),
        NEITHER
    );
    assert_eq!(
        abutting_verdicts(
            &object_of("intron_variant", &offsets),
            Some(&object_of("intron_variant", &[]))
        ),
        NEITHER
    );
    let other_pair = [
        ("hgvs_intronstartoffset", Value::Int(3)),
        ("hgvs_intronendoffset", Value::Int(0)),
    ];
    assert_eq!(
        abutting_verdicts(
            &object_of(upstream, &other_pair),
            Some(&object_of(upstream, &[]))
        ),
        NEITHER
    );
}

/// Through the comparison itself the class sets aside both fields of the made-up upstream
/// row, counted once per row, and nothing where its condition does not hold: the same pair
/// on a row inside the transcript, a field missing on the reference's side and a pair the
/// mechanism does not give stay differences.
#[test]
fn the_abutting_insertion_class_is_a_built_in_class() {
    let filled = ";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0";
    let row = |term: &str, offsets: &str| {
        default_text(&[default_row(
            "ins1",
            "2:1000443-1000444",
            "T",
            "ENST00000000021",
            term,
            &format!("{ABUTTING_EXTRA}{offsets}"),
        )])
    };
    let upstream = "upstream_gene_variant";
    assert_eq!(
        classes::tables().field_classes_of(ABUTTING),
        Some(HGVS_INTRON_OFFSET_FIELDS)
    );
    let (report, documented) =
        compare_default(&row(upstream, filled), &row(upstream, ""), &[ABUTTING]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(ABUTTING));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &masked_once(ABUTTING)
    );
    assert_eq!(
        verdict_of(&row(upstream, filled), &row(upstream, ""), &[]),
        "red"
    );
    for (reference, candidate) in [
        (row("intron_variant", filled), row("intron_variant", "")),
        (
            row(upstream, ";HGVS_IntronStartOffset=0"),
            row(upstream, ""),
        ),
        (
            row(upstream, ";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=5"),
            row(upstream, ""),
        ),
    ] {
        assert_eq!(
            default_under(ABUTTING, &reference, &candidate),
            ("red".to_string(), nothing())
        );
    }
}

const ABUTTING_GRID_PAIRS: [(&str, &str); 8] = [
    ("0", "0"),
    ("1", "0"),
    ("0", "-1"),
    ("-1", "0"),
    ("0", "1"),
    ("1", "-1"),
    ("-1", "-1"),
    ("2", "0"),
];
const ABUTTING_GRID_TERMS: [&str; 3] = [
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
];
const ABUTTING_GRID_NOTATIONS: [Option<(&str, &str)>; 5] = [
    None,
    Some(("reference", "HGVSc=ENST00000000041.1:n.1_2insA")),
    Some(("candidate", "HGVSc=ENST00000000041.1:n.1_2insA")),
    Some(("reference", "HGVS_OFFSET=2")),
    Some(("candidate", "HGVS_OFFSET=2")),
];

/// One cell of the made-up grid of the abutting class: the plugin's pair, the row's term,
/// its strand, its DISTANCE, an HGVSc or an HGVS_OFFSET on either side, the candidate's two
/// fields blank or filled, and an insertion's Location or a one-base one.
#[derive(Clone, Copy, Debug)]
struct AbuttingCell {
    pair: (&'static str, &'static str),
    term: &'static str,
    strand: Option<&'static str>,
    distance: &'static str,
    notation: Option<(&'static str, &'static str)>,
    filled: bool,
    insertion: bool,
}

/// Every combination of the grid's values, the last varying fastest.
fn abutting_grid() -> Vec<AbuttingCell> {
    let mut cells = Vec::new();
    for pair in ABUTTING_GRID_PAIRS {
        for term in ABUTTING_GRID_TERMS {
            for strand in [Some("1"), Some("-1"), None] {
                for distance in ["0", "5"] {
                    for notation in ABUTTING_GRID_NOTATIONS {
                        for filled in [false, true] {
                            for insertion in [true, false] {
                                cells.push(AbuttingCell {
                                    pair,
                                    term,
                                    strand,
                                    distance,
                                    notation,
                                    filled,
                                    insertion,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: the row is an insertion's outside the
/// transcript at a DISTANCE of 0, no HGVSc and no HGVS_OFFSET stands on either side, the
/// candidate writes neither field, and the pair is 0 and 0, or 1 and 0 on a downstream row
/// whose strand is not forward, or 0 and -1 on an upstream row whose strand is not reverse.
fn abutting_grid_set_aside(cell: &AbuttingCell) -> bool {
    if !cell.insertion
        || cell.term == "intron_variant"
        || cell.distance != "0"
        || cell.notation.is_some()
        || cell.filled
    {
        return false;
    }
    match cell.pair {
        ("0", "0") => true,
        ("1", "0") => cell.term == "downstream_gene_variant" && cell.strand != Some("1"),
        ("0", "-1") => cell.term == "upstream_gene_variant" && cell.strand != Some("-1"),
        _ => false,
    }
}

/// The reference's and the candidate's default-format rows of a cell, under a Location
/// and a Feature of the cell's own.
fn abutting_grid_rows(cell: &AbuttingCell, number: usize) -> (String, String) {
    let (start, end) = cell.pair;
    let position = 1_000_000 + 10 * number;
    let location = if cell.insertion {
        format!("2:{position}-{}", position + 1)
    } else {
        format!("2:{position}")
    };
    let mut extra = format!("IMPACT=MODIFIER;DISTANCE={}", cell.distance);
    if let Some(strand) = cell.strand {
        extra.push_str(&format!(";STRAND={strand}"));
    }
    let mut reference = extra.clone();
    let mut candidate = extra;
    match cell.notation {
        Some(("reference", held)) => reference.push_str(&format!(";{held}")),
        Some((_, held)) => candidate.push_str(&format!(";{held}")),
        None => {}
    }
    reference.push_str(&format!(
        ";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}"
    ));
    if cell.filled {
        candidate.push_str(";HGVS_IntronEndOffset=7;HGVS_IntronStartOffset=7");
    }
    let row = |extra: &str| {
        default_row(
            &format!("cell{number}"),
            &location,
            "T",
            &format!("ENST{number:011}"),
            cell.term,
            extra,
        )
    };
    (row(&reference), row(&candidate))
}

/// Over the made-up grid the class sets aside both of the plugin's fields exactly where
/// the rule's statement does: the pair 0 and 0; 1 and 0 past a reverse-strand transcript's
/// last base; 0 and -1 before a forward-strand transcript's first base; a row without
/// STRAND admits either. Every other cell stays a difference.
#[test]
fn the_abutting_insertion_class_takes_its_three_pairs_and_nothing_wider() {
    let wrong: Vec<AbuttingCell> = abutting_grid()
        .into_iter()
        .enumerate()
        .filter(|(number, cell)| {
            let (reference, candidate) = abutting_grid_rows(cell, *number);
            let read = offset_fields()
                .map(|field| offset_grid_verdict(ABUTTING, &reference, &candidate, field));
            read != [abutting_grid_set_aside(cell); 2]
        })
        .map(|(_, cell)| cell)
        .collect();
    assert!(wrong.is_empty(), "{wrong:?}");
}

const PROTEIN_NOTATION: &str = "ENSP00000000001.1:p.Ala84ValfsTer?";

/// `text` with the protein notation `notation` before its HGVS_OFFSET.
fn with_protein(text: &str, notation: &str) -> String {
    text.replace(";HGVS_OFFSET=", &format!(";HGVSp={notation};HGVS_OFFSET="))
}

/// An insertion's or deletion's row annotated in an intron (no cDNA, CDS or protein position)
/// whose HGVSc and non-zero HGVS_OFFSET both engines write alike, the HGVSc placing the shifted
/// allele in the coding sequence, where Ensembl VEP leaves the HGVSp empty and the candidate
/// writes a protein notation, is the class, counted once per row in the default, tab, VCF and
/// JSON formats: a deletion on a coding base, a duplication ending on one with its first
/// coordinate back in the intron, and an insertion between two coding bases; the offset of
/// either sign.
#[test]
fn the_protein_gate_sets_aside_the_notation_the_reference_lacks() {
    let reference = |hgvsc: &str, offset: &str, allele: &str, location: &str| {
        let extra = format!("IMPACT=HIGH;STRAND=1;HGVSc={hgvsc};HGVS_OFFSET={offset}");
        default_text(&[default_row(
            "intron_end_indel",
            location,
            allele,
            "ENST00000000001",
            "splice_acceptor_variant",
            &extra,
        )])
    };
    let plain = reference("ENST00000000001.1:c.251del", "1", "-", "21:25001999");
    let candidate = |text: &str| with_protein(text, PROTEIN_NOTATION);

    let (report, documented) = compare_default(&plain, &candidate(&plain), &[PROTEIN_GATE]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(PROTEIN_GATE));
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &masked_once(PROTEIN_GATE)
    );
    assert_eq!(verdict_of(&plain, &candidate(&plain), &[]), "red");
    for (hgvsc, allele, location) in [
        (
            "ENST00000000001.1:c.251_252del",
            "-",
            "21:25001998-25001999",
        ),
        ("ENST00000000001.1:c.251dup", "C", "21:25001998-25001999"),
        (
            "ENST00000000001.1:c.251-2_252dup",
            "TTTA",
            "21:25001996-25001997",
        ),
        (
            "ENST00000000001.1:c.251_252insTG",
            "TG",
            "21:25001998-25001999",
        ),
    ] {
        for offset in ["3", "-3"] {
            let shifted = reference(hgvsc, offset, allele, location);
            assert_eq!(
                verdict_of(&shifted, &candidate(&shifted), &[PROTEIN_GATE]),
                "green",
                "{hgvsc} {offset}"
            );
        }
    }
    // The tab, VCF and JSON formats carry the fields apart and count the row once.
    let header = tab_header(&["IMPACT", "STRAND", "HGVSc", "HGVSp", "HGVS_OFFSET"]);
    let tab_row = "intron_end_indel\t21:25001999\t-\tENSG1\tENST00000000001\tTranscript\tsplice_acceptor_variant\t-\t-\t-\t-\t-\t-\tHIGH\t1\tENST00000000001.1:c.251del\t";
    let tab = |tail: &str| tab_entries(&format!("{header}\n{tab_row}{tail}\n"));
    assert_eq!(
        under_class(
            PROTEIN_GATE,
            &tab("-\t1"),
            &tab(&format!("{PROTEIN_NOTATION}\t1"))
        ),
        ("green".to_string(), once(PROTEIN_GATE))
    );
    let inputs = vec!["21\t25001998\tintron_end_indel\tAC\tA\t.\t.\t.".to_string()];
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|HGVSc|HGVSp|HGVS_OFFSET",
    );
    let csq_row = format!(
        "{}CSQ=-|splice_acceptor_variant|ENST00000000001|Transcript||||ENST00000000001.1:c.251del|",
        record_columns(&inputs[0])
    );
    let vcf = |tail: &str| vcf_entries(&format!("{header}{csq_row}{tail}\n"), &inputs);
    assert_eq!(
        under_class(
            PROTEIN_GATE,
            &vcf("|1"),
            &vcf(&format!("{PROTEIN_NOTATION}|1"))
        ),
        ("green".to_string(), once(PROTEIN_GATE))
    );
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("intron_end_indel")),
        ("allele_string", text("C/-")),
        ("most_severe_consequence", text("splice_acceptor_variant")),
    ];
    let consequence = object(&[
        ("variant_allele", text("-")),
        ("transcript_id", text("ENST00000000001")),
        (
            "consequence_terms",
            Value::List(vec![text("splice_acceptor_variant")]),
        ),
        ("hgvsc", text("ENST00000000001.1:c.251del")),
        ("hgvs_offset", Value::Int(1)),
    ]);
    let object_of = |fields: &[(&str, Value)]| {
        json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        )
    };
    let written = [("hgvsp", text(PROTEIN_NOTATION))];
    assert_eq!(
        under_class(PROTEIN_GATE, &object_of(&[]), &object_of(&written)),
        ("green".to_string(), once(PROTEIN_GATE))
    );
    let report = compare(&object_of(&[]), &object_of(&written), &mut no_documented());
    assert_eq!(verdict(&report), "red");
}

/// What a row of an indel shifted out of an intron prints, an empty text for a field the
/// row leaves out.
#[derive(Clone, Copy)]
struct Shifted<'a> {
    hgvsc: &'a str,
    offset: &'a str,
    hgvsp: &'a str,
    consequence: &'a str,
    allele: &'a str,
    more: &'a str,
}

impl Default for Shifted<'_> {
    fn default() -> Self {
        Shifted {
            hgvsc: "ENST00000000001.1:c.251del",
            offset: "1",
            hgvsp: "",
            consequence: "splice_acceptor_variant",
            allele: "-",
            more: "",
        }
    }
}

fn shifted_row(row: Shifted) -> String {
    let mut extra = "IMPACT=HIGH;STRAND=1".to_string();
    for (name, value) in [
        ("HGVSc", row.hgvsc),
        ("HGVSp", row.hgvsp),
        ("HGVS_OFFSET", row.offset),
    ] {
        if !value.is_empty() {
            extra.push_str(&format!(";{name}={value}"));
        }
    }
    extra.push_str(row.more);
    default_text(&[default_row(
        "intron_end_indel",
        "21:25001999",
        row.allele,
        "ENST00000000001",
        row.consequence,
        &extra,
    )])
}

/// The nearest shapes stay differences: the same row with the sides exchanged; the allele
/// shifted into an exon that may open on the first coding base, out of an intron of the
/// 5' UTR, which the rows alone do not place: a notation naming coding base 1 as the base
/// the intron ends before, or one of coding bases alone with a shift at least as long as
/// its coordinate; an insertion the shift leaves before the exon's first base, a
/// duplication ending on the intron's last base or an insertion between that base and the
/// exon's first, of any bases; an insertion before the stop codon's last base, or a
/// candidate notation that replaces the stop by a residue with no extension. So do a
/// candidate value that is no protein notation, a row
/// without an offset or with a zero one, an HGVSc or an offset the engines write
/// differently, a shifted HGVSc that leaves the allele in the intron, in a UTR or on a
/// non-coding transcript, a substitution's notation, another transcript's notation, a row
/// with a cDNA, CDS or protein position, a row annotated outside the transcript, a row
/// whose SHIFT_LENGTH is positive and a substitution's row.
#[test]
fn the_protein_gate_keeps_the_nearest_shapes_as_differences() {
    let plain = Shifted::default();
    let written = Shifted {
        hgvsp: PROTEIN_NOTATION,
        ..plain
    };
    let read = |reference: Shifted, candidate: Shifted| {
        verdict_of(
            &shifted_row(reference),
            &shifted_row(candidate),
            &[PROTEIN_GATE],
        )
    };
    // The verdict of a row against itself with the protein notation on the candidate's side.
    let taken = |row: Shifted| {
        read(
            row,
            Shifted {
                hgvsp: PROTEIN_NOTATION,
                ..row
            },
        )
    };
    assert_eq!(taken(plain), "green");
    // The sides exchanged: the notation Ensembl VEP writes and the candidate lacks.
    assert_eq!(read(written, plain), "red");
    let other_residue = PROTEIN_NOTATION.replace("84", "85");
    assert_eq!(
        read(
            written,
            Shifted {
                hgvsp: &other_residue,
                ..plain
            }
        ),
        "red"
    );
    // The exon may open on coding base 1: the notation names it, or its coding bases lie no
    // further in than the shift is long. An exon opening on base 2 or later is the class,
    // and so is a copy that ends on a coding base and begins before the coding start.
    for (hgvsc, offset) in [
        ("c.1del", "1"),
        ("c.2_3del", "2"),
        ("c.2_3del", "-2"),
        ("c.1-1dup", "1"),
        ("c.3dup", "4"),
        ("c.1-3_1dup", "4"),
        ("c.1-3_1-1dup", "3"),
        ("c.251-3_252-1dup", "3"),
        ("c.249+3_251-1dup", "9"),
        ("c.250+3_251+4dup", "1"),
        ("c.250+1_250+2insTG", "1"),
        ("c.1-1_1insTG", "1"),
        ("c.3_4insTG", "4"),
        ("c.251del", "251"),
        ("c.251del", "-300"),
    ] {
        let shifted = format!("ENST00000000001.1:{hgvsc}");
        let row = Shifted {
            hgvsc: &shifted,
            offset,
            ..plain
        };
        assert_eq!(taken(row), "red", "{hgvsc} {offset}");
    }
    for (hgvsc, offset) in [
        ("c.2del", "1"),
        ("c.3dup", "3"),
        ("c.251del", "-250"),
        ("c.2-2_2dup", "4"),
        ("c.1+3_2dup", "9"),
        ("c.-5+3_2dup", "9"),
    ] {
        let shifted = format!("ENST00000000001.1:{hgvsc}");
        let row = Shifted {
            hgvsc: &shifted,
            offset,
            ..plain
        };
        assert_eq!(taken(row), "green", "{hgvsc} {offset}");
    }
    assert_eq!(
        read(
            plain,
            Shifted {
                hgvsp: "ENST00000000001.1:c.251del",
                ..plain
            }
        ),
        "red"
    );
    // The stop replaced by a residue with no extension is neither of the class's stop forms.
    for stop in [
        PROTEIN_NOTATION.replace("Ala84Val", "Ter84Gln"),
        PROTEIN_NOTATION.replace("p.Ala84Val", "p.*84Gln"),
    ] {
        assert_eq!(
            read(
                plain,
                Shifted {
                    hgvsp: &stop,
                    ..plain
                }
            ),
            "red",
            "{stop}"
        );
    }
    for offset in ["", "0"] {
        assert_eq!(taken(Shifted { offset, ..plain }), "red", "{offset}");
    }
    assert_eq!(
        read(
            plain,
            Shifted {
                hgvsc: "ENST00000000001.1:c.252del",
                ..written
            }
        ),
        "red"
    );
    assert_eq!(
        read(
            plain,
            Shifted {
                offset: "2",
                ..written
            }
        ),
        "red"
    );
    assert_eq!(taken(Shifted { hgvsc: "", ..plain }), "red");
    for hgvsc in [
        "ENST00000000001.1:c.251-3del",
        "ENST00000000001.1:c.251-2_251-1del",
        "ENST00000000001.1:c.250+1del",
        "ENST00000000001.1:c.-5del",
        "ENST00000000001.1:c.*5del",
        "ENST00000000001.1:n.251del",
        "ENST00000000001.1:c.251-2dup",
        "ENST00000000001.1:c.251-2_251-1insTG",
        "ENST00000000001.1:c.250+1_250+2insTG",
        "ENST00000000001.1:c.251C>T",
        "ENST00000000002.1:c.251del",
    ] {
        assert_eq!(taken(Shifted { hgvsc, ..plain }), "red", "{hgvsc}");
    }
    for column in 0..3 {
        let mut cells = ["-", "-", "-"];
        cells[column] = "84";
        let positioned = shifted_row(plain).replacen(
            "\t-\t-\t-\t-\t-\t-\t",
            &format!("\t{}\t-\t-\t-\t", cells.join("\t")),
            1,
        );
        assert_eq!(
            verdict_of(
                &positioned,
                &with_protein(&positioned, PROTEIN_NOTATION),
                &[PROTEIN_GATE]
            ),
            "red",
            "{column}"
        );
    }
    for consequence in ["upstream_gene_variant", "downstream_gene_variant"] {
        assert_eq!(
            taken(Shifted {
                consequence,
                ..plain
            }),
            "red",
            "{consequence}"
        );
    }
    assert_eq!(
        taken(Shifted {
            more: ";SHIFT_LENGTH=1",
            ..plain
        }),
        "red"
    );
    // An insertion the shift leaves before the exon's first base is a change at the
    // junction, which reads as the intron's new end as well as the exon's new start: a
    // duplication ending on the intron's last base, or an insertion between that base and
    // the exon's first, of any bases, stays a difference. So do an insertion before the stop
    // codon's last base and one after an exon's last base.
    for hgvsc in [
        "c.251-2_251-1dup",
        "c.251-6_251-1dup",
        "c.250+3_251-1dup",
        "c.4-4_4-1dup",
        "c.*-3_*-1dup",
        "c.251-1_251insTAG",
        "c.251-1_251insAG",
        "c.251-1_251insAC",
        "c.*-1_251insCAG",
        "c.*-1dup",
        "c.*-1_251insTG",
        "c.251-1dup",
        "c.251-1_251insA",
        "c.251-1_251insG",
        "c.250_250+1insTG",
    ] {
        let shifted = format!("ENST00000000001.1:{hgvsc}");
        let row = Shifted {
            hgvsc: &shifted,
            offset: "3",
            ..plain
        };
        assert_eq!(taken(row), "red", "{hgvsc}");
    }
    assert_eq!(
        taken(Shifted {
            more: ";SHIFT_LENGTH=0",
            ..plain
        }),
        "green"
    );
    assert_eq!(
        taken(Shifted {
            allele: "T",
            ..plain
        }),
        "red"
    );
}

/// The reference's and the candidate's default-format text of one row of an indel shifted
/// out of an intron: the same HGVSc and HGVS_OFFSET on both, the HGVSp on the candidate's
/// alone.
fn shifted_protein_pair(hgvsc: &str, offset: &str, hgvsp: &str) -> (String, String) {
    let extra = format!("IMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:{hgvsc}");
    let row = |protein: &str| {
        default_text(&[default_row(
            "intron_end_indel",
            "21:25001999",
            "-",
            "ENST00000000001",
            "splice_acceptor_variant",
            &format!("{extra}{protein};HGVS_OFFSET={offset}"),
        )])
    };
    (row(""), row(&format!(";HGVSp={hgvsp}")))
}

/// Under the class, the row of an insertion the shift leaves between an intron's last base
/// and the exon's first base is a difference with a protein notation on the candidate's side
/// alone and with one on the reference's side alone, at a shift of one base and of more, of
/// either sign, and the class is credited with nothing.
fn assert_an_insertion_at_an_introns_end_stays_a_difference(hgvsc: &str) {
    for offset in ["1", "-1", "2", "3", "-9"] {
        let (without, with_notation) = shifted_protein_pair(hgvsc, offset, PROTEIN_NOTATION);
        for (reference, candidate) in [(&without, &with_notation), (&with_notation, &without)] {
            assert_eq!(
                default_under(PROTEIN_GATE, reference, candidate),
                ("red".to_string(), nothing()),
                "{hgvsc} {offset}"
            );
        }
    }
}

/// The row of an indel the shift carries wholly into the coding sequence, the notation on
/// the candidate's side alone, is the class, counted once.
fn assert_the_shifted_indel_class_takes(hgvsc: &str, offset: &str) {
    let (reference, candidate) = shifted_protein_pair(hgvsc, offset, PROTEIN_NOTATION);
    assert_eq!(
        default_under(PROTEIN_GATE, &reference, &candidate),
        ("green".to_string(), once(PROTEIN_GATE)),
        "{hgvsc} {offset}"
    );
}

/// The copy of an intron's one last base, written back in the intron and carried to the
/// exon's edge, is a change at the junction: both engines leave its HGVSp empty, so a
/// notation on one side alone is a difference at any shift. The copy of the exon's first
/// base, which the shift carries past the junction, is the class.
#[test]
fn the_protein_gate_leaves_a_copy_of_the_introns_last_base_a_difference() {
    assert_an_insertion_at_an_introns_end_stays_a_difference("c.251-1dup");
    for offset in ["1", "-1", "3"] {
        assert_the_shifted_indel_class_takes("c.251dup", offset);
    }
}

/// Bases that are no copy of the intron's end, one or more, which the shift leaves between
/// the intron's last base and the exon's first base, are a change at the junction: a
/// notation on one side alone is a difference at any shift. The same bases carried between
/// two coding bases are the class.
#[test]
fn the_protein_gate_leaves_bases_inserted_at_an_introns_end_a_difference() {
    for hgvsc in [
        "c.251-1_251insA",
        "c.251-1_251insG",
        "c.251-1_251insTG",
        "c.251-1_251insTGT",
    ] {
        assert_an_insertion_at_an_introns_end_stays_a_difference(hgvsc);
    }
    for hgvsc in ["c.251_252insA", "c.251_252insTG", "c.251_252insTGT"] {
        for offset in ["1", "-1", "3"] {
            assert_the_shifted_indel_class_takes(hgvsc, offset);
        }
    }
}

/// A duplication of two or more bases that ends on the intron's last base is a change at
/// the junction: a notation on one side alone is a difference at any shift. A duplication
/// that starts in the intron and ends on a coding base is the class.
#[test]
fn the_protein_gate_leaves_a_duplication_ending_on_the_introns_last_base_a_difference() {
    for hgvsc in ["c.251-2_251-1dup", "c.251-6_251-1dup", "c.250+3_251-1dup"] {
        assert_an_insertion_at_an_introns_end_stays_a_difference(hgvsc);
    }
    for hgvsc in ["c.251-2_252dup", "c.251-1_251dup"] {
        for offset in ["1", "-1", "3"] {
            assert_the_shifted_indel_class_takes(hgvsc, offset);
        }
    }
}

/// A duplicated, deleted or replaced range that ends at the stop codon holds the stop's
/// residue, and a notation in parentheses that replaces the stop by a residue with no
/// extension starts there: each changes the stop codon and stays a difference, in either
/// spelling of the stop. The same forms on residues before the stop, a notation in
/// parentheses that starts before it and an extension named at the stop in parentheses are
/// the class.
#[test]
fn the_protein_gate_refuses_a_range_ending_at_the_stop_codon_and_a_stop_in_parentheses() {
    let protein = "ENSP00000000001.1:p.";
    let read = |change: &str| {
        let (reference, candidate) =
            shifted_protein_pair("c.251del", "1", &format!("{protein}{change}"));
        default_under(PROTEIN_GATE, &reference, &candidate)
    };
    for change in [
        "Gln83_Ter84dup",
        "Gln83_Ter84del",
        "Gln83_Ter84delinsArg",
        "Gln83_*84dup",
        "Gln83_*84del",
        "Gln83_*84delinsArg",
        "(Gln83_Ter84dup)",
        "(Ter84Gln)",
        "(*84Gln)",
    ] {
        assert_eq!(read(change), ("red".to_string(), nothing()), "{change}");
    }
    for change in [
        "Gln82_Gln83dup",
        "Gln82_Gln83del",
        "Gln82_Gln83delinsArg",
        "(Gln82_Gln83dup)",
        "(Ala84ValfsTer?)",
        "(Ter84GlnextTer?)",
    ] {
        assert_eq!(
            read(change),
            ("green".to_string(), once(PROTEIN_GATE)),
            "{change}"
        );
    }
}

/// A candidate notation that names the stop codon as a flank of an insertion stays a
/// difference in the default, tab, VCF and JSON formats, the stop codon spelled `Ter` or
/// `*`: bases inserted inside the stop codon that leave a stop codon after the new residues
/// read like bases inserted before it, and the class takes the stop codon in two forms
/// alone, the stop unchanged and an extension named at it. The same row with the insertion
/// named between two residues is the class.
#[test]
fn the_protein_gate_keeps_an_insertion_beside_the_stop_codon_as_a_difference() {
    let hgvsc = "ENST00000000001.1:c.252_253insTGT";
    let between_residues = "ENSP00000000001.1:p.Ala84_Gly85insVal";
    let beside_the_stop = [
        "ENSP00000000001.1:p.Ala84_Ter85insVal",
        "ENSP00000000001.1:p.Ala84_*85insVal",
    ];
    let extra = format!("IMPACT=HIGH;STRAND=1;HGVSc={hgvsc};HGVS_OFFSET=3");
    let reference = default_text(&[default_row(
        "intron_end_indel",
        "21:25001998-25001999",
        "TGT",
        "ENST00000000001",
        "splice_acceptor_variant",
        &extra,
    )]);
    assert_eq!(
        default_under(
            PROTEIN_GATE,
            &reference,
            &with_protein(&reference, between_residues)
        ),
        ("green".to_string(), once(PROTEIN_GATE))
    );
    for notation in beside_the_stop {
        assert_eq!(
            default_under(
                PROTEIN_GATE,
                &reference,
                &with_protein(&reference, notation)
            ),
            ("red".to_string(), nothing()),
            "{notation}"
        );
    }
    let header = tab_header(&["IMPACT", "STRAND", "HGVSc", "HGVSp", "HGVS_OFFSET"]);
    let tab_row = format!(
        "intron_end_indel\t21:25001998-25001999\tTGT\tENSG1\tENST00000000001\tTranscript\tsplice_acceptor_variant\t-\t-\t-\t-\t-\t-\tHIGH\t1\t{hgvsc}\t"
    );
    let inputs = vec!["21\t25001998\tintron_end_indel\tA\tATGT\t.\t.\t.".to_string()];
    let csq = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|HGVSc|HGVSp|HGVS_OFFSET",
    );
    let csq_row = format!(
        "{}CSQ=TGT|splice_acceptor_variant|ENST00000000001|Transcript||||{hgvsc}|",
        record_columns(&inputs[0])
    );
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("intron_end_indel")),
        ("allele_string", text("-/TGT")),
        ("most_severe_consequence", text("splice_acceptor_variant")),
    ];
    let consequence = object(&[
        ("variant_allele", text("TGT")),
        ("transcript_id", text("ENST00000000001")),
        (
            "consequence_terms",
            Value::List(vec![text("splice_acceptor_variant")]),
        ),
        ("hgvsc", text(hgvsc)),
        ("hgvs_offset", Value::Int(3)),
    ]);
    let object_of = |fields: &[(&str, Value)]| {
        json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        )
    };
    // What each format's pair of outputs reads as: the verdict and what was applied.
    let verdicts = |notation: &str| {
        [
            under_class(
                PROTEIN_GATE,
                &tab_entries(&format!("{header}\n{tab_row}-\t3\n")),
                &tab_entries(&format!("{header}\n{tab_row}{notation}\t3\n")),
            ),
            under_class(
                PROTEIN_GATE,
                &vcf_entries(&format!("{csq}{csq_row}|3\n"), &inputs),
                &vcf_entries(&format!("{csq}{csq_row}{notation}|3\n"), &inputs),
            ),
            under_class(
                PROTEIN_GATE,
                &object_of(&[]),
                &object_of(&[("hgvsp", text(notation))]),
            ),
        ]
    };
    let taken = ("green".to_string(), once(PROTEIN_GATE));
    assert_eq!(
        verdicts(between_residues),
        [taken.clone(), taken.clone(), taken]
    );
    let left = ("red".to_string(), nothing());
    for notation in beside_the_stop {
        assert_eq!(
            verdicts(notation),
            [left.clone(), left.clone(), left.clone()],
            "{notation}"
        );
    }
}

/// A row on a transcript flagged `cds_end_NF`, alone or beside `cds_start_NF`, whose notation
/// changes a residue before the stop is a row of the class, as the same row without a flag
/// or flagged `cds_start_NF` alone is: the candidate writes the notation that the codons of
/// such a sequence give, in the default, tab, VCF and JSON formats. With a flag on one
/// engine's row alone the HGVSp is the class's and the `FLAGS` field stays a difference.
#[test]
fn the_protein_gate_takes_a_row_on_a_transcript_flagged_for_either_end_of_its_coding_sequence() {
    let hgvsc = "ENST00000000001.1:c.251del";
    let notation = "ENSP00000000001.1:p.Leu84fs";
    let reference = |flags: &str| {
        let extra = format!("IMPACT=HIGH;STRAND=1{flags};HGVSc={hgvsc};HGVS_OFFSET=1");
        default_text(&[default_row(
            "intron_end_indel",
            "21:25001999",
            "-",
            "ENST00000000001",
            "splice_acceptor_variant",
            &extra,
        )])
    };
    let candidate = |text: &str| with_protein(text, notation);
    let taken = ("green".to_string(), once(PROTEIN_GATE));

    let plain = reference("");
    assert_eq!(
        default_under(PROTEIN_GATE, &plain, &candidate(&plain)),
        taken
    );
    let flag_on_one_row = ("red".to_string(), once(PROTEIN_GATE));
    for (flags, both_rows, one_row) in [
        ("cds_start_NF", &taken, &flag_on_one_row),
        ("cds_end_NF", &taken, &flag_on_one_row),
        ("cds_start_NF,cds_end_NF", &taken, &flag_on_one_row),
    ] {
        let flagged = reference(&format!(";FLAGS={flags}"));
        for (reference_text, candidate_text, expected) in [
            (&flagged, candidate(&flagged), both_rows),
            (&plain, candidate(&flagged), one_row),
            (&flagged, candidate(&plain), one_row),
        ] {
            assert_eq!(
                &default_under(PROTEIN_GATE, reference_text, &candidate_text),
                expected,
                "{flags}"
            );
        }
    }

    let header = tab_header(&["IMPACT", "STRAND", "FLAGS", "HGVSc", "HGVSp", "HGVS_OFFSET"]);
    let tab = |flags: &str, hgvsp: &str| {
        let row = format!(
            "intron_end_indel\t21:25001999\t-\tENSG1\tENST00000000001\tTranscript\tsplice_acceptor_variant\t-\t-\t-\t-\t-\t-\tHIGH\t1\t{flags}\t{hgvsc}\t{hgvsp}\t1"
        );
        tab_entries(&format!("{header}\n{row}\n"))
    };
    let inputs = vec!["21\t25001998\tintron_end_indel\tAC\tA\t.\t.\t.".to_string()];
    let csq = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|FLAGS|HGVSc|HGVSp|HGVS_OFFSET",
    );
    let vcf = |flags: &str, hgvsp: &str| {
        vcf_entries(
            &format!(
                "{csq}{}CSQ=-|splice_acceptor_variant|ENST00000000001|Transcript||||{flags}|{hgvsc}|{hgvsp}|1\n",
                record_columns(&inputs[0])
            ),
            &inputs,
        )
    };
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("intron_end_indel")),
        ("allele_string", text("C/-")),
        ("most_severe_consequence", text("splice_acceptor_variant")),
    ];
    let consequence = object(&[
        ("variant_allele", text("-")),
        ("transcript_id", text("ENST00000000001")),
        (
            "consequence_terms",
            Value::List(vec![text("splice_acceptor_variant")]),
        ),
        ("hgvsc", text(hgvsc)),
        ("hgvs_offset", Value::Int(1)),
    ]);
    let object_of = |fields: &[(&str, Value)]| {
        json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        )
    };
    let written = ("hgvsp", text(notation));

    assert_eq!(
        under_class(PROTEIN_GATE, &tab("-", "-"), &tab("-", notation)),
        taken
    );
    assert_eq!(
        under_class(PROTEIN_GATE, &vcf("", ""), &vcf("", notation)),
        taken
    );
    assert_eq!(
        under_class(
            PROTEIN_GATE,
            &object_of(&[]),
            &object_of(std::slice::from_ref(&written))
        ),
        taken
    );
    for (flags, expected) in [
        (vec!["cds_start_NF"], &taken),
        (vec!["cds_end_NF"], &taken),
        (vec!["cds_start_NF", "cds_end_NF"], &taken),
    ] {
        let joined = flags.join(",");
        assert_eq!(
            &under_class(PROTEIN_GATE, &tab(&joined, "-"), &tab(&joined, notation)),
            expected,
            "{joined}"
        );
        let escaped = flags.join("&");
        assert_eq!(
            &under_class(PROTEIN_GATE, &vcf(&escaped, ""), &vcf(&escaped, notation)),
            expected,
            "{joined}"
        );
        let listed = (
            "flags",
            Value::List(flags.iter().map(|flag| text(flag)).collect()),
        );
        assert_eq!(
            &under_class(
                PROTEIN_GATE,
                &object_of(std::slice::from_ref(&listed)),
                &object_of(&[listed, written.clone()])
            ),
            expected,
            "{joined}"
        );
    }
}

/// One row of an insertion or deletion the record writes on exon bases of the 5' UTR, where
/// coding base 1 is base 6001 of chromosome 21 and base 201 of the cDNA: the three position
/// columns, the shifted allele's HGVSc and HGVS_OFFSET, the protein change the candidate
/// writes for it and any further Extra text.
#[derive(Clone, Copy)]
struct UtrRow<'a> {
    location: &'a str,
    allele: &'a str,
    consequence: &'a str,
    positions: [&'a str; 3],
    hgvsc: &'a str,
    offset: &'a str,
    change: &'a str,
    more: &'a str,
}

/// A deletion of one base written three bases before the coding start and carried onto
/// coding base 1.
const UTR_DELETION: UtrRow = UtrRow {
    location: "21:5998",
    allele: "-",
    consequence: "5_prime_UTR_variant",
    positions: ["198", "-", "-"],
    hgvsc: "c.1del",
    offset: "3",
    change: "Met1?",
    more: "",
};

/// An insertion of one base written between the fourth and the third base before the coding
/// start and carried four bases on, a copy of coding base 1.
const UTR_DUPLICATION: UtrRow = UtrRow {
    location: "21:5997-5998",
    allele: "A",
    consequence: "5_prime_UTR_variant",
    positions: ["197-198", "-", "-"],
    hgvsc: "c.1dup",
    offset: "4",
    change: "Met1?",
    more: "",
};

/// An insertion of three bases written 108 bases before the place the shift leaves it,
/// between coding bases 3 and 4.
const UTR_INSERTION: UtrRow = UtrRow {
    location: "21:5895-5896",
    allele: "GCT",
    consequence: "5_prime_UTR_variant",
    positions: ["95-96", "-", "-"],
    hgvsc: "c.3_4insGCT",
    offset: "108",
    change: "Met1_Val2insAla",
    more: "",
};

/// The rows of an allele written on exon bases of the 5' UTR that the shift carries into the
/// coding sequence: the three above; a copy of fourteen bases and one of three, each
/// beginning two bases before the coding start on a row of the splice region, the second
/// leaving the start codon as it is; a deletion of nine coding bases; and an insertion
/// between coding bases 9 and 10.
const UTR_ROWS: [UtrRow; 7] = [
    UTR_DELETION,
    UTR_DUPLICATION,
    UTR_INSERTION,
    UtrRow {
        location: "21:5998-5999",
        allele: "CTGCTGCTGCTGCT",
        consequence: "splice_region_variant,5_prime_UTR_variant",
        positions: ["198-199", "-", "-"],
        hgvsc: "c.-2_12dup",
        offset: "14",
        change: "Gly5AlafsTer12",
        more: "",
    },
    UtrRow {
        location: "21:5998-5999",
        allele: "CTG",
        consequence: "splice_region_variant,5_prime_UTR_variant",
        positions: ["198-199", "-", "-"],
        hgvsc: "c.-2_1dup",
        offset: "3",
        change: "Met1%3D",
        more: "",
    },
    UtrRow {
        location: "21:5985-5993",
        allele: "-",
        consequence: "5_prime_UTR_variant",
        positions: ["185-193", "-", "-"],
        hgvsc: "c.4_12del",
        offset: "19",
        change: "Val2_Ser4del",
        more: "",
    },
    UtrRow {
        location: "21:5985-5986",
        allele: "CAG",
        consequence: "5_prime_UTR_variant",
        positions: ["185-186", "-", "-"],
        hgvsc: "c.9_10insCAG",
        offset: "24",
        change: "Leu3_Ser4insGln",
        more: "",
    },
];

/// The default-format text of `row`, its HGVSp the protein change `change` where one is
/// given.
fn utr_text(row: UtrRow, change: &str) -> String {
    let mut extra = format!(
        "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000000001.1:{}",
        row.hgvsc
    );
    if !change.is_empty() {
        extra.push_str(&format!(";HGVSp=ENSP00000000001.1:p.{change}"));
    }
    if !row.offset.is_empty() {
        extra.push_str(&format!(";HGVS_OFFSET={}", row.offset));
    }
    extra.push_str(row.more);
    let [cdna, cds, protein] = row.positions;
    default_text(&[format!(
        "utr_indel\t{}\t{}\tENSG1\tENST00000000001\tTranscript\t{}\t{cdna}\t{cds}\t{protein}\t-\t-\t-\t{extra}",
        row.location, row.allele, row.consequence
    )])
}

/// What the class reads of `row` with the protein change on the candidate's side alone: the
/// verdict and what was applied.
fn utr_pair_under_the_class(row: UtrRow) -> (String, std::collections::BTreeMap<String, i64>) {
    default_under(PROTEIN_GATE, &utr_text(row, ""), &utr_text(row, row.change))
}

/// An allele the record writes on exon bases of the 5' UTR reads `utr` and not `coding`
/// under the predicate Ensembl VEP tests before it writes an HGVSp
/// (BaseVariationFeatureOverlapAllele.pm:523-525, TranscriptVariationAllele.pm:1731), so
/// where the shift carries it into the coding sequence the reference's row has the shifted
/// HGVSc and HGVS_OFFSET and no HGVSp. Such a row carries a cDNA position, no CDS or protein
/// position and the term `5_prime_UTR_variant`, and is the class from the row alone, counted
/// once, with an offset of either sign: a deletion of coding bases, a duplication that ends
/// on a coding base and begins on one or before the coding start, and an insertion between
/// two coding bases. The candidate's notation may be the start codon's own, lost or written
/// as unchanged.
#[test]
fn the_protein_gate_takes_an_allele_written_on_exon_bases_of_the_five_prime_utr() {
    for row in UTR_ROWS {
        let negative = format!("-{}", row.offset);
        for offset in [row.offset, negative.as_str()] {
            assert_eq!(
                utr_pair_under_the_class(UtrRow { offset, ..row }),
                ("green".to_string(), once(PROTEIN_GATE)),
                "{} {offset}",
                row.hgvsc
            );
        }
        assert_eq!(
            verdict_of(&utr_text(row, ""), &utr_text(row, row.change), &[]),
            "red",
            "{}",
            row.hgvsc
        );
    }
}

/// The nearest shapes of such a row stay differences: the sides exchanged, and two
/// notations that differ; a shifted HGVSc that leaves the allele outside the coding
/// sequence, a duplication that ends before the coding start or begins in an intron of the
/// coding region or in the 3' UTR, a deletion or an insertion with a coordinate before the
/// coding start, and a notation of a non-coding transcript; a row with a CDS or a protein
/// position; a row with a cDNA position whose terms lack `5_prime_UTR_variant`; a row whose
/// SHIFT_LENGTH is positive; a row without an offset or with a zero one; an HGVSc or an
/// offset the engines write differently; and a substitution's row.
#[test]
fn the_protein_gate_keeps_the_nearest_shapes_of_a_five_prime_utr_row_as_differences() {
    let taken = ("green".to_string(), once(PROTEIN_GATE));
    let left = ("red".to_string(), nothing());
    let read = |reference: &str, candidate: &str| default_under(PROTEIN_GATE, reference, candidate);
    let deletion = UTR_DELETION;
    assert_eq!(utr_pair_under_the_class(deletion), taken);
    assert_eq!(
        read(
            &utr_text(deletion, deletion.change),
            &utr_text(deletion, "")
        ),
        left
    );
    assert_eq!(
        read(&utr_text(deletion, "Met1?"), &utr_text(deletion, "Met1%3D")),
        left
    );
    for (row, hgvsc) in [
        (UTR_DUPLICATION, "c.-3_-1dup"),
        (UTR_DUPLICATION, "c.-1dup"),
        (UTR_DUPLICATION, "c.4-2_6dup"),
        (UTR_DUPLICATION, "c.*2_6dup"),
        (UTR_DUPLICATION, "n.201dup"),
        (deletion, "c.-2_1del"),
        (deletion, "c.-5_-1del"),
        (deletion, "c.-1del"),
        (deletion, "c.4-2_6del"),
        (deletion, "n.201del"),
        (UTR_INSERTION, "c.-1_1insGCT"),
        (UTR_INSERTION, "c.-2_-1insGCT"),
    ] {
        assert_eq!(
            utr_pair_under_the_class(UtrRow { hgvsc, ..row }),
            left,
            "{hgvsc}"
        );
    }
    for positions in [["198", "1", "-"], ["198", "-", "1"], ["198", "1", "1"]] {
        assert_eq!(
            utr_pair_under_the_class(UtrRow {
                positions,
                ..deletion
            }),
            left,
            "{positions:?}"
        );
    }
    for consequence in [
        "splice_region_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "5_prime_UTR_premature_start_codon_gain_variant",
    ] {
        assert_eq!(
            utr_pair_under_the_class(UtrRow {
                consequence,
                ..deletion
            }),
            left,
            "{consequence}"
        );
    }
    assert_eq!(
        utr_pair_under_the_class(UtrRow {
            more: ";SHIFT_LENGTH=3",
            ..deletion
        }),
        left
    );
    assert_eq!(
        utr_pair_under_the_class(UtrRow {
            more: ";SHIFT_LENGTH=0",
            ..deletion
        }),
        taken
    );
    for offset in ["", "0"] {
        assert_eq!(
            utr_pair_under_the_class(UtrRow { offset, ..deletion }),
            left,
            "{offset}"
        );
    }
    let written = utr_text(deletion, deletion.change);
    for candidate in [
        written.replace("c.1del", "c.2del"),
        written.replace("HGVS_OFFSET=3", "HGVS_OFFSET=4"),
    ] {
        assert_eq!(read(&utr_text(deletion, ""), &candidate), left);
    }
    assert_eq!(
        utr_pair_under_the_class(UtrRow {
            allele: "T",
            ..deletion
        }),
        left
    );
}

/// A row of an allele written on exon bases of the 5' UTR is the class in the tab, VCF and
/// JSON formats too, counted once: the cDNA position in its column, CSQ field or pair of
/// keys, the terms joined as each format joins them, and the start codon written as
/// unchanged with its `=` escaped in the tab and VCF formats and plain in JSON. The same
/// row without the term `5_prime_UTR_variant` stays a difference in each.
#[test]
fn the_protein_gate_takes_a_five_prime_utr_row_in_the_tab_vcf_and_json_formats() {
    let hgvsc = "ENST00000000001.1:c.-2_1dup";
    let (escaped, plain) = ("ENSP00000000001.1:p.Met1%3D", "ENSP00000000001.1:p.Met1=");
    let header = tab_header(&["IMPACT", "STRAND", "HGVSc", "HGVSp", "HGVS_OFFSET"]);
    let tab = |terms: &[&str], hgvsp: &str| {
        tab_entries(&format!(
            "{header}\nutr_indel\t21:5998-5999\tCTG\tENSG1\tENST00000000001\tTranscript\t{}\t198-199\t-\t-\t-\t-\t-\tMODIFIER\t1\t{hgvsc}\t{hgvsp}\t3\n",
            terms.join(",")
        ))
    };
    let inputs = vec!["21\t5998\tutr_indel\tA\tACTG\t.\t.\t.".to_string()];
    let csq = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|HGVSc|HGVSp|HGVS_OFFSET",
    );
    let vcf = |terms: &[&str], hgvsp: &str| {
        vcf_entries(
            &format!(
                "{csq}{}CSQ=CTG|{}|ENST00000000001|Transcript|198-199|||{hgvsc}|{hgvsp}|3\n",
                record_columns(&inputs[0]),
                terms.join("&")
            ),
            &inputs,
        )
    };
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("utr_indel")),
        ("allele_string", text("-/CTG")),
        ("most_severe_consequence", text("splice_region_variant")),
    ];
    let json = |terms: &[&str], hgvsp: Option<&str>| {
        let mut members = vec![
            ("variant_allele", text("CTG")),
            ("transcript_id", text("ENST00000000001")),
            (
                "consequence_terms",
                Value::List(terms.iter().map(|term| text(term)).collect()),
            ),
            ("cdna_start", Value::Int(199)),
            ("cdna_end", Value::Int(198)),
            ("hgvsc", text(hgvsc)),
            ("hgvs_offset", Value::Int(3)),
        ];
        if let Some(hgvsp) = hgvsp {
            members.push(("hgvsp", text(hgvsp)));
        }
        json_entries(&json_line(&record, object(&members)), &inputs)
    };
    let taken = ("green".to_string(), once(PROTEIN_GATE));
    let left = ("red".to_string(), nothing());
    for (terms, expected) in [
        (
            &["splice_region_variant", "5_prime_UTR_variant"][..],
            &taken,
        ),
        (&["splice_region_variant"][..], &left),
    ] {
        assert_eq!(
            &under_class(PROTEIN_GATE, &tab(terms, "-"), &tab(terms, escaped)),
            expected,
            "tab {terms:?}"
        );
        assert_eq!(
            &under_class(PROTEIN_GATE, &vcf(terms, ""), &vcf(terms, escaped)),
            expected,
            "vcf {terms:?}"
        );
        assert_eq!(
            &under_class(PROTEIN_GATE, &json(terms, None), &json(terms, Some(plain))),
            expected,
            "json {terms:?}"
        );
    }
}

/// One default-format row under an identifier: IMPACT and STRAND, then `pairs` in the
/// Extra text.
fn offset_row(
    location: &str,
    feature: &str,
    consequence: &str,
    pairs: &str,
    allele: &str,
) -> String {
    default_text(&[default_row(
        "shifted_deletion",
        location,
        allele,
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

/// A deletion the shift carries over the last base of a transcript: Ensembl VEP writes no
/// HGVSc, no HGVS_OFFSET and neither field of the HGVSIntronOffset plugin; the candidate
/// writes a del notation, its offset and the two fields. The fields are the class's with the
/// notation, the row counted once, where they are the intron offsets of the notation's two
/// coordinates, in either order. The same values stay differences beside a notation whose
/// coordinates give others, on a substitution's row, on a row annotated upstream and on a
/// row without the candidate's notation; so does a field the candidate lacks, or writes with
/// another value, where the reference writes it.
#[test]
fn the_written_span_gate_holds_the_plugin_offsets_to_a_del_notation() {
    let row = |pairs: &str| {
        offset_row(
            "1:1000108-1000111",
            "ENST00000000034",
            "non_coding_transcript_exon_variant",
            pairs,
            "-",
        )
    };
    let notation = |coordinates: &str| format!("ENST00000000034.1:n.{coordinates}");
    let reference = row("");
    let candidate = row(&intron_offset_pairs(&notation("348del"), 2, 0, 0));
    let (report, documented) = compare_default(&reference, &candidate, &[GATE]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), once(GATE));
    for column in ["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"] {
        assert_eq!(
            at(&report, &["columns", column, "masked"]),
            &masked_once(GATE),
            "{column}"
        );
    }
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    for (start, end) in [(9, -944), (-944, 9)] {
        let intronic = row(&intron_offset_pairs(
            &notation("29+9_30-944del"),
            2,
            start,
            end,
        ));
        assert_eq!(
            verdict_of(&reference, &intronic, &[GATE]),
            "green",
            "{start} {end}"
        );
    }
    for (coordinates, start, end) in [
        ("348del", 5, 0),
        ("29+9_30-944del", 9, 0),
        ("29+9_30-944del", 0, -944),
    ] {
        let other = row(&intron_offset_pairs(&notation(coordinates), 2, start, end));
        let (report, documented) = compare_default(&reference, &other, &[GATE]);
        assert_eq!(verdict(&report), "red", "{coordinates} {start} {end}");
        assert_eq!(
            counts(&documented.applied),
            once(GATE),
            "{coordinates} {start} {end}"
        );
        for column in ["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"] {
            assert_eq!(number(&report, &["columns", column, "differing"]), 1);
        }
    }
    // A coordinate before the stop codon's last base prints no exon base: nothing to hold to.
    let late = row(&intron_offset_pairs(
        "ENST00000000034.1:c.*-5_*3del",
        2,
        -5,
        0,
    ));
    assert_eq!(verdict_of(&reference, &late, &[GATE]), "green");
    let written = intron_offset_pairs(&notation("348del"), 2, 0, 0);
    for (location, consequence, allele) in [
        ("1:1000108", "non_coding_transcript_exon_variant", "T"),
        ("1:1000108-1000111", "upstream_gene_variant", "-"),
    ] {
        let other =
            |pairs: &str| offset_row(location, "ENST00000000034", consequence, pairs, allele);
        assert_eq!(
            default_under(GATE, &other(""), &other(&written)),
            ("red".to_string(), nothing()),
            "{location} {consequence}"
        );
    }
    assert_eq!(
        default_under(
            GATE,
            &reference,
            &row(";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0")
        ),
        ("red".to_string(), nothing())
    );
    for planted in [
        intron_offset_pairs(&notation("348del"), 2, 5, 0),
        ";HGVSc=ENST00000000034.1:n.348del;HGVS_OFFSET=2".to_string(),
    ] {
        let (report, documented) = compare_default(&candidate, &row(&planted), &[GATE]);
        assert_eq!(verdict(&report), "red", "{planted}");
        assert_eq!(counts(&documented.applied), nothing(), "{planted}");
        assert_eq!(
            number(&report, &["columns", "HGVS_IntronStartOffset", "differing"]),
            1
        );
    }
}

/// The intron offset of a coordinate is the signed number after its exon base, 0 for a
/// coordinate that names an exon base.
#[test]
fn the_intron_offset_of_a_coordinate_is_the_number_after_its_exon_base() {
    assert_eq!(coordinate_intron_offset("29+9"), "9");
    assert_eq!(coordinate_intron_offset("30-944"), "-944");
    assert_eq!(coordinate_intron_offset("-12"), "0");
    assert_eq!(coordinate_intron_offset("*7"), "0");
}
