// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rule base on small outputs: the rows the parser builds, the rows of a record, the
//! helpers the plugins' rules share and the classes that are one test each. The expected
//! values are those `tools/vep-formats-compare` and `tools/vep-outputs-compare` give on the
//! same inputs.

use std::collections::BTreeSet;
use std::path::Path;
use std::rc::Rc;
use vep_plugin_pairs::rules::entry::{
    allele_is_indel, as_number, canonical_json_text, consequence_set, decode_extra, dedup_keys,
    entry_kind, feature_kind, field_named, field_of, field_of_any, parse_default,
    reference_allele_is_empty, replaced_reference, text_lines, uploaded_alleles, uploaded_position,
    DefaultHeader, Entry, Fields, Row, ABSENT,
};
use vep_plugin_pairs::rules::record::{
    line_writes_row, perl_split, record_label, vep_rs_pos, RecordRows, VEP_RS_COORDINATE_CEILING,
};
use vep_plugin_pairs::rules::shared::{
    allele_on_transcript_strand, at_its_location, cache_integer, cdna_span, cdna_span_closed,
    first_alternate, indel_row, insertion_row, integer_of, is_integer, is_non_zero_integer,
    list_parts, location_span, notation_indel_length, positive_shift, read_shard, read_text,
    record_allele_lengths, record_alleles, record_chromosome, record_indel_row, reverse_complement,
    sequence_location, shard_transcripts, sibling_rows, vcf_alleles_as_spelt, vcf_first_alternate,
    STRAND_FIELDS,
};
use vep_plugin_pairs::rules::{
    order_only_fields, same_items, BuiltIn, Consult, Refusal, Rules, JSON_ARRAY,
};

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
const UTR_LOSS: &str = "5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;5UTR_annotation=CapDistanceToStart=63:DistanceToCDS=110:KozakStrength=Moderate:type=uORF";
const UTR_NONE: &str = "5UTR_consequence=-;5UTR_annotation=-";
const UTR_INPUTS: [&str; 3] = [
    "21\t25734954\tutr_ins\tA\tACAT\t.\t.\t.",
    "21\t500\tmnv\tAC\tGT\t.\t.\t.",
    "21\t600\tsnv\tA\tT\t.\t.\t.",
];

/// A row of a transcript in the 5' UTR with the UTRAnnotator fields `utr`.
fn utr_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    cdna: &str,
    utr: &str,
) -> String {
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t5_prime_UTR_variant\t{cdna}\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;{utr}\n")
}

fn insertion_row_text(allele: &str, utr: &str) -> String {
    utr_row(
        "utr_ins",
        "21:25734954-25734955",
        allele,
        "ENST00000000033",
        "64-65",
        utr,
    )
}

fn row(columns: [&str; 8]) -> String {
    let [name, location, allele, feature, kind, consequence, cdna, extra] = columns;
    format!("{name}\t{location}\t{allele}\tG1\t{feature}\t{kind}\t{consequence}\t{cdna}\t-\t-\t-\t-\t-\t{extra}\n")
}

/// Eight rows of three records: one whose first column spells its alleles, two under
/// identifiers.
fn rows() -> Vec<Row> {
    let text = [
        HEADER.to_string(),
        row([
            "21_100_A/AT/G",
            "21:100-101",
            "T",
            "ENST1",
            "Transcript",
            "missense_variant",
            "-",
            "STRAND=-1;SHIFT_LENGTH=3;HGVSc=ENST1.1:c.5_6insT;DISTANCE=4",
        ]),
        row([
            "21_100_A/AT/G",
            "21:100",
            "G",
            "ENST1",
            "Transcript",
            "missense_variant",
            "-",
            "STRAND=1;SHIFT_LENGTH=0",
        ]),
        row([
            "21_100_A/AT/G",
            "21:100",
            "G",
            "ENST2",
            "Transcript",
            "missense_variant",
            "-",
            "STRAND=1",
        ]),
        row([
            "rs9",
            "21:200-202",
            "-",
            "ENST1",
            "Transcript",
            "missense_variant",
            "?-12/900",
            "STRAND=-1;SHIFT_LENGTH=2",
        ]),
        row([
            "rs9",
            "21:200-202",
            "AC",
            "ENST1",
            "Transcript",
            "missense_variant",
            "10-12",
            "IMPACT=LOW",
        ]),
        row([
            "rs9",
            "21:200",
            "deletion",
            "ENST1",
            "-",
            "missense_variant",
            "7",
            "-",
        ]),
        row([
            "rs8",
            "X:300-301",
            "ACG",
            "-",
            "-",
            "intergenic_variant",
            "-",
            "-",
        ]),
        row([
            "rs8",
            "X:300-301",
            "ACG",
            "-",
            "-",
            "intergenic_variant&upstream_gene_variant,intergenic_variant",
            "-",
            "-",
        ]),
    ]
    .concat();
    parse_default(&text).unwrap()
}

fn keys(rows: &[Row]) -> Vec<&str> {
    rows.iter().map(|e| e.key.as_str()).collect()
}

fn refusal<T: std::fmt::Debug>(result: Result<T, Refusal>) -> String {
    result.unwrap_err().0
}

#[test]
fn the_parser_reads_the_committed_reference_output() {
    let corpus = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/expected/default.txt.gz");
    let entries = parse_default(&read_text(&corpus).unwrap()).unwrap();
    assert_eq!(entries.len(), 977);
    let first = &entries[0];
    assert_eq!(first.key, "21:11000001|N.|-|-");
    assert_eq!(first.location.as_deref(), Some("21:11000001"));
    assert_eq!(first.consequence, "intergenic_variant");
    assert_eq!(
        first.fields.get("Uploaded_variation"),
        Some("bnd_intergenic_pair_a")
    );
    assert_eq!(first.fields.get("Extra"), Some("IMPACT=MODIFIER"));
    assert_eq!(first.fields.get("IMPACT"), Some("MODIFIER"));
    let ablation = entries
        .iter()
        .find(|e| {
            e.fields.get("Uploaded_variation") == Some("del_10_6_mb")
                && e.feature == "ENST00000262354"
        })
        .unwrap();
    assert_eq!(ablation.consequence, "transcript_ablation");
    assert_eq!(
        ablation.fields.get("Extra"),
        Some("IMPACT=HIGH;STRAND=-1;OverlapBP=23932;OverlapPC=100.00")
    );
    assert_eq!(ablation.fields.get("OverlapBP"), Some("23932"));
    assert_eq!(ablation.fields.get("OverlapPC"), Some("100.00"));
    assert_eq!(ablation.fields.get("STRAND"), Some("-1"));
    let distinct: BTreeSet<&str> = keys(&entries).into_iter().collect();
    assert_eq!(distinct.len(), 977);
}

#[test]
fn ga4gh_vrs_is_read_as_a_parsed_object() {
    let reference = concat!(
        r#"{"type":"Allele","state":{"sequence":"C","type":"SequenceState"},"#,
        r#""location":{"type":"SequenceLocation","sequence_id":"refseq:NC_000021.9","#,
        r#""interval":{"type":"SimpleInterval","end":10538885,"start":"10538884"}}}"#
    );
    let engine = concat!(
        r#"{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval","#,
        r#""start":"10538884","end":10538885},"sequence_id":"refseq:NC_000021.9"},"#,
        r#""state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}"#
    );
    let canonical = canonical_json_text(reference).unwrap();
    assert_eq!(canonical, canonical_json_text(engine).unwrap());
    assert_eq!(
        canonical,
        concat!(
            r#"{"location":{"interval":{"end":10538885,"start":"10538884","type":"SimpleInterval"},"#,
            r#""sequence_id":"refseq:NC_000021.9","type":"SequenceLocation"},"#,
            r#""state":{"sequence":"C","type":"SequenceState"},"type":"Allele"}"#
        )
    );
    for as_written in ["21:10538884:A:C", "-", "1e5", "[1, 2]", "{\"a\":1,}"] {
        assert_eq!(canonical_json_text(as_written).unwrap(), as_written);
    }
    let with = |value: &str| {
        let extra = format!("IMPACT=MODIFIER;GA4GH_VRS={value};STRAND=1");
        let text = [
            HEADER.to_string(),
            row([
                "rs1",
                "21:10538885",
                "C",
                "T1",
                "Transcript",
                "intron_variant",
                "-",
                &extra,
            ]),
        ]
        .concat();
        Rc::clone(&parse_default(&text).unwrap()[0])
    };
    let (reference_entry, engine_entry) = (with(reference), with(engine));
    assert_eq!(
        reference_entry.fields.get("GA4GH_VRS"),
        Some(canonical.as_str())
    );
    assert_eq!(
        engine_entry.fields.get("GA4GH_VRS"),
        Some(canonical.as_str())
    );
    assert_eq!(
        reference_entry.fields.get("Extra"),
        engine_entry.fields.get("Extra")
    );
    assert_eq!(
        reference_entry.fields.get("Extra"),
        Some(format!("IMPACT=MODIFIER;GA4GH_VRS={canonical};STRAND=1").as_str())
    );
    assert_eq!(reference_entry.fields.get("STRAND"), Some("1"));
    // A pair whose value is not an object, and a bare key, are rebuilt as written.
    let dash_entry = with("-");
    assert_eq!(
        dash_entry.fields.get("Extra"),
        Some("IMPACT=MODIFIER;GA4GH_VRS=-;STRAND=1")
    );
    assert_eq!(dash_entry.fields.get("GA4GH_VRS"), Some("-"));
    let bare = [
        HEADER.to_string(),
        row([
            "rs1",
            "21:10538885",
            "C",
            "T1",
            "Transcript",
            "intron_variant",
            "-",
            "IMPACT=MODIFIER;GA4GH_VRS;STRAND=1",
        ]),
    ]
    .concat();
    let bare_entry = &parse_default(&bare).unwrap()[0];
    assert_eq!(
        bare_entry.fields.get("Extra"),
        Some("IMPACT=MODIFIER;GA4GH_VRS;STRAND=1")
    );
    assert_eq!(bare_entry.fields.get("GA4GH_VRS"), Some(""));
    let column = format!(
        "#Uploaded_variation\tLocation\tAllele\tGA4GH_VRS\nrs1\t21:10538885\tC\t{reference}\n"
    );
    assert_eq!(
        parse_default(&column).unwrap()[0].fields.get("GA4GH_VRS"),
        Some(canonical.as_str())
    );
    let deep = format!("{}1{}", "{\"a\":".repeat(498), "}".repeat(498));
    assert_eq!(
        refusal(canonical_json_text(&deep)),
        "a JSON object is nested too deep to read: {\"a\":{\"a\":{\"a\":{\"a\":{\"a\"..."
    );
}

#[test]
fn exact_duplicates_of_a_key_are_counted() {
    let mut entries: Vec<Entry> = (0..3)
        .map(|_| Entry {
            key: "k".to_string(),
            ..Entry::default()
        })
        .collect();
    entries.push(Entry {
        key: "other".to_string(),
        ..Entry::default()
    });
    dedup_keys(&mut entries);
    let found: Vec<&str> = entries.iter().map(|e| e.key.as_str()).collect();
    assert_eq!(found, ["k", "k#2", "k#3", "other"]);
    assert_eq!(
        keys(
            &parse_default(
                "#Location\tAllele\tFeature\tFeature_type\n1:5\tA\tf\tt\r\n1:5\tA\tf\tt\n"
            )
            .unwrap()
        ),
        ["1:5|A|f|t", "1:5|A|f|t#2"]
    );
}

#[test]
fn a_row_of_another_column_count_is_refused() {
    assert_eq!(
        refusal(parse_default("#A\tB\n1\t2\t3\n")),
        "column count on line: 1\t2\t3"
    );
    assert_eq!(
        refusal(parse_default("1\t2\n")),
        "column count on line: 1\t2"
    );
    assert!(parse_default("").unwrap().is_empty());
    assert!(parse_default("## x\n#A\n").unwrap().is_empty());
    assert_eq!(text_lines("a\r\nb\n\nc"), ["a", "b", "", "c"]);
    assert_eq!(text_lines("a\n"), ["a"]);
    assert!(text_lines("").is_empty());
}

#[test]
fn a_row_read_under_its_header_is_the_entry_the_parser_reads_of_the_two_lines() {
    let fourteen = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";
    let headers = [
        fourteen,
        "#Location\tExtra",
        "#Location\tExtra\r",
        "##Location\tExtra",
        "Location\tExtra",
        "",
        "#",
        "#Location\tExtra\n#Allele\tExtra",
    ];
    let rows = [
        "rs1\t1:100\tT\tENSG1\tENST1\tTranscript\tmissense_variant\t1\t-\t-\t-\t-\t-\tIMPACT=LOW;PLANT_score=0.5",
        "1_100_A/T\t1:100\tT\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
        "1:100\tPLANT_score=1;PLANT_label=a%3Db",
        "1:100\tExtra=taken;PLANT_score=1",
        "1:100\tA=1;GA4GH_VRS={\"b\": 1, \"a\": 2};Z=3",
        "1:100\t-",
        "1:100\t-\r",
        "1:100\tA=1\tB=2",
        "1:100",
        "",
        "\r",
        "#Allele\tExtra",
        "##note",
        "1:100\tA=1\n1:200\tA=2",
        "1:100\tA=1\n",
    ];
    for header in headers {
        let read = DefaultHeader::new(header);
        for row in rows {
            let parsed = parse_default(&format!("{header}\n{row}"))
                .map(|entries| entries.into_iter().next());
            assert_eq!(
                format!("{:?}", read.row(row)),
                format!("{parsed:?}"),
                "{header:?} over {row:?}"
            );
        }
    }
    let entry = DefaultHeader::new("#Location\tExtra")
        .row("1:100\tExtra=taken;PLANT_score=1")
        .unwrap()
        .unwrap();
    assert_eq!(
        entry.fields.iter().collect::<Vec<_>>(),
        [
            ("Location", "1:100"),
            ("Extra", "Extra=taken;PLANT_score=1"),
            ("PLANT_score", "1")
        ]
    );
}

#[test]
fn a_name_written_twice_keeps_its_first_place_and_its_last_value() {
    let entries =
        parse_default("#Location\tLocation\tExtra\n1:1\t2:2\tLocation=3:3;Extra=q;A;B=;=C;A=2\n")
            .unwrap();
    assert_eq!(entries[0].key, "3:3|||");
    let fields: Vec<(&str, &str)> = entries[0].fields.iter().collect();
    assert_eq!(
        fields,
        [
            ("Location", "3:3"),
            ("Extra", "Location=3:3;Extra=q;A;B=;=C;A=2"),
            ("A", "2"),
            ("B", ""),
            ("", "C")
        ]
    );
    let mut many = Fields::new();
    for number in 0..200 {
        many.insert(&format!("name{}", number % 120), number.to_string());
    }
    many.set_default("name3", "kept out");
    many.set_default("new", "set");
    assert_eq!(many.len(), 121);
    assert_eq!(many.get("name3"), Some("123"));
    assert_eq!(many.get("name119"), Some("119"));
    assert_eq!(many.get("new"), Some("set"));
    assert!(many.contains("name0") && !many.contains("name120"));
    assert_eq!(many.iter().nth(1), Some(("name1", "121")));
    assert_eq!(many.first_of(&["absent", "name5"], ABSENT), "125");
    assert_eq!(many.get_or("absent", "-"), "-");
}

#[test]
fn the_parser_reads_a_records_alleles_from_a_first_column_that_spells_them() {
    let rows = rows();
    assert_eq!(
        keys(&rows),
        [
            "21:100-101|T|ENST1|Transcript",
            "21:100|G|ENST1|Transcript",
            "21:100|G|ENST2|Transcript",
            "21:200-202|-|ENST1|Transcript",
            "21:200-202|AC|ENST1|Transcript",
            "21:200|deletion|ENST1|-",
            "X:300-301|ACG|-|-",
            "X:300-301|ACG|-|-#2"
        ]
    );
    let shapes: Vec<(Option<bool>, Option<bool>, Option<String>)> = rows
        .iter()
        .map(|e| (e.indel, e.insertion, e.reference.borrow().clone()))
        .collect();
    assert_eq!(shapes[0], (None, Some(false), Some(String::new())));
    assert_eq!(shapes[1], (Some(false), Some(false), Some("A".to_string())));
    assert_eq!(shapes[2], (Some(false), Some(false), Some("A".to_string())));
    assert!(shapes[3..].iter().all(|shape| *shape == (None, None, None)));
    let names: Vec<&str> = rows[0].fields.iter().map(|(name, _)| name).collect();
    assert_eq!(
        names[13..],
        ["Extra", "STRAND", "SHIFT_LENGTH", "HGVSc", "DISTANCE"]
    );
    assert_eq!(
        rows[0].fields.get("Extra"),
        Some("STRAND=-1;SHIFT_LENGTH=3;HGVSc=ENST1.1:c.5_6insT;DISTANCE=4")
    );
    let labels: Vec<String> = rows.iter().map(|e| record_label(e)).collect();
    assert_eq!(
        labels,
        [
            "21:100-101 21_100_A/AT/G",
            "21:100 21_100_A/AT/G",
            "21:100 21_100_A/AT/G",
            "21:200-202 rs9",
            "21:200-202 rs9",
            "21:200 rs9",
            "X:300-301 rs8",
            "X:300-301 rs8"
        ]
    );
    let kinds: Vec<&str> = rows.iter().map(|e| entry_kind(e)).collect();
    assert_eq!(
        kinds,
        [
            "Transcript",
            "Transcript",
            "Transcript",
            "Transcript",
            "Transcript",
            "-",
            "-",
            "-"
        ]
    );
    let lists: Vec<String> = rows.iter().map(|e| feature_kind(entry_kind(e))).collect();
    assert_eq!(lists[4], "transcript");
    assert_eq!(lists[5], "intergenic");
    let sets: Vec<String> = rows
        .iter()
        .map(|e| consequence_set(&e.consequence))
        .collect();
    assert_eq!(sets[6], "intergenic_variant");
    assert_eq!(sets[7], "intergenic_variant,upstream_gene_variant");
}

#[test]
fn a_generated_record_name_spells_its_position_and_alleles() {
    let alleles = |name: &str| {
        uploaded_alleles(name).map(|(reference, alternates)| (reference, alternates.join(",")))
    };
    let spelt =
        |reference: &str, alternates: &str| Some((reference.to_string(), alternates.to_string()));
    assert_eq!(alleles("21_10542414_A/T/AG"), spelt("A", "T,AG"));
    assert_eq!(
        alleles("chr1_KI270706v1_random_100_A/AAC"),
        spelt("A", "AAC")
    );
    assert_eq!(alleles("1_100_A/T_5_C/G"), spelt("C", "G"));
    assert_eq!(alleles("21_100_A/T\n"), spelt("A", "T"));
    assert_eq!(alleles("21_100_*/T"), spelt("*", "T"));
    assert_eq!(alleles("21_100_a-/<DEL>"), spelt("a-", "<DEL>"));
    assert_eq!(alleles("21_100_A/G"), spelt("A", "G"));
    assert_eq!(uploaded_position("21_10542414_A/T/AG"), Ok(Some(10542414)));
    assert_eq!(uploaded_position("1_100_A/T_5_C/G"), Ok(Some(5)));
    assert_eq!(uploaded_position("7_100998_ACG/GTTCA"), Ok(Some(100998)));
    assert_eq!(uploaded_position("7_100997_tacg/t/-"), Ok(Some(100997)));
    assert_eq!(
        uploaded_position("chrUn_made_up_100998_A/T"),
        Ok(Some(100998))
    );
    assert_eq!(uploaded_position("21_100_A/G"), Ok(Some(100)));
    for identifier in [
        "rs123",
        "del_padded",
        "21_100_A/T\n\n",
        "21_100_A/",
        "21_100_/T",
        "2 1_100_A/T",
        "x_1_A/T/",
        "_1_A/T",
        "21_100_A/T G",
        "21_\u{661}\u{660}\u{660}_A/G",
        "1993671",
        "rs1993671",
        "7_100998",
        "made_up_7",
    ] {
        assert_eq!(uploaded_alleles(identifier), None, "{identifier:?}");
        assert_eq!(uploaded_position(identifier), Ok(None), "{identifier:?}");
    }
    assert_eq!(
        refusal(uploaded_position(&format!("21_{}_A/G", "9".repeat(40)))),
        "Uploaded_variation position '999999999999999999999999...' is not an integer"
    );
}

#[test]
fn an_alleles_shape_is_read_from_the_records_alleles_as_ensembl_vep_spells_them() {
    let indel = |reference: &str, alternates: &[&str], allele: &str| {
        allele_is_indel(reference, alternates, allele, false)
    };
    assert_eq!(indel("A", &["AT"], "T"), Some(true));
    assert_eq!(indel("A", &["G"], "G"), Some(false));
    assert_eq!(indel("AT", &["A"], "-"), Some(true));
    assert_eq!(indel("A", &["<DEL>"], "deletion"), None);
    assert_eq!(indel("G", &["[21:25881361[G"], "[21:25881361[G"), None);
    // A lone alternate keeps its shape however Ensembl VEP shortens the allele; the parts of
    // an allele string are already shortened, so they are read as they are.
    assert_eq!(indel("TTTTT", &["TTTTTT"], "T"), Some(true));
    assert_eq!(indel("TCT", &["TCTC"], "C"), Some(true));
    assert_eq!(allele_is_indel("T", &["TT", "-"], "TT", true), Some(true));
    assert_eq!(allele_is_indel("T", &["TT", "-"], "-", true), Some(true));
    // Ensembl VEP drops the shared leading base only when every alternate shares it.
    assert_eq!(indel("A", &["AT", "G"], "AT"), Some(true));
    assert_eq!(indel("A", &["AT", "G"], "G"), Some(false));
    assert_eq!(indel("A", &["AT", "G"], "T"), None);
    assert_eq!(indel("A", &["AT", "T"], "AT"), Some(true));
    assert_eq!(indel("A", &["AT", "T"], "T"), Some(false));
    assert_eq!(indel("A", &["AT", "*"], "T"), Some(true));
    // An allele spelt with a `*` is exempt from the dropped base whatever else it carries.
    assert_eq!(indel("A", &["AG", "<DEL:*>"], "G"), Some(true));
    assert_eq!(
        allele_is_indel("-", &["G", "<DEL:*>"], "G", true),
        Some(true)
    );
    assert_eq!(indel("A", &["T", "<DEL:*>"], "T"), Some(false));
    assert_eq!(indel("AC", &["ACT", "AG"], "CT"), Some(true));
    assert_eq!(indel("AC", &["ACT", "AG"], "G"), Some(false));
    assert!(reference_allele_is_empty("A", &["AT", "*"], false));
    assert!(reference_allele_is_empty("A", &["AT", "AG"], false));
    assert!(!reference_allele_is_empty("TCT", &["TCTC"], false));
    assert!(reference_allele_is_empty("-", &["T"], true));
    assert!(!reference_allele_is_empty("A", &["AT"], true));
    assert!(!reference_allele_is_empty("GACCTG", &["G"], false));
    assert!(!reference_allele_is_empty("A", &["T", "G"], false));
}

/// The shape, POS, REF, the alternates, Location, Allele, the reference bases.
type Read<'a> = (
    &'a str,
    i128,
    &'a str,
    &'a [&'a str],
    &'a str,
    &'a str,
    &'a str,
);
/// The shape, POS, REF, the alternates, Location, Allele.
type NotRead<'a> = (
    &'a str,
    i128,
    &'a str,
    &'a [&'a str],
    Option<&'a str>,
    &'a str,
);

#[test]
fn the_bases_an_allele_stands_in_for_are_read_by_position_from_the_records_alleles() {
    let read: [Read; 8] = [
        (
            "three bases replaced by five",
            100998,
            "ACG",
            &["GTTCA"],
            "7:100998-101000",
            "GTTCA",
            "ACG",
        ),
        (
            "three bases replaced by five behind a shared first base",
            100997,
            "TACG",
            &["TGTTCA"],
            "7:100998-101000",
            "GTTCA",
            "ACG",
        ),
        (
            "one base replaced by one",
            101000,
            "A",
            &["G"],
            "7:101000",
            "G",
            "A",
        ),
        (
            "the second of two alternates",
            100998,
            "ACG",
            &["TTT", "GTTCA"],
            "7:100998-101000",
            "GTTCA",
            "ACG",
        ),
        (
            "a pure insertion",
            100999,
            "A",
            &["AGTT"],
            "7:100999-101000",
            "GTT",
            "",
        ),
        (
            "a pure insertion behind two shared bases",
            100998,
            "AC",
            &["ACGTT"],
            "7:100999-101000",
            "GTT",
            "",
        ),
        (
            "a deletion",
            100997,
            "TACG",
            &["T"],
            "7:100998-101000",
            "-",
            "ACG",
        ),
        (
            "lower-case alleles",
            100998,
            "acg",
            &["gttca"],
            "7:100998-101000",
            "GTTCA",
            "ACG",
        ),
    ];
    for (shape, position, reference, alternates, location, allele, bases) in read {
        assert_eq!(
            replaced_reference(position, reference, alternates, Some(location), allele),
            Ok(Some(bases.to_string())),
            "{shape}"
        );
    }
    let not_read: [NotRead; 6] = [
        (
            "an alternate that writes no such row",
            100998,
            "ACG",
            &["GTTCAT"],
            Some("7:100998-101000"),
            "GTTCA",
        ),
        (
            "a Location past the record's bases",
            100998,
            "ACG",
            &["GTTCA"],
            Some("7:100999-101001"),
            "GTTCA",
        ),
        (
            "two alternates that read the row differently",
            100999,
            "AC",
            &["GTT", "AGTTC"],
            Some("7:100999-101000"),
            "GTT",
        ),
        (
            "an allele that is not a sequence",
            100998,
            "ACG",
            &["<DEL>"],
            Some("7:100998-101000"),
            "deletion",
        ),
        (
            "a Location that names no whole positions",
            100998,
            "ACG",
            &["GTTCA"],
            Some("7:100998-101000.5"),
            "GTTCA",
        ),
        (
            "a row without a Location",
            100998,
            "ACG",
            &["GTTCA"],
            None,
            "GTTCA",
        ),
    ];
    for (shape, position, reference, alternates, location, allele) in not_read {
        assert_eq!(
            replaced_reference(position, reference, alternates, location, allele),
            Ok(None),
            "{shape}"
        );
    }
    assert_eq!(
        refusal(replaced_reference(
            1,
            "A",
            &["G"],
            Some(&format!("7:{}", "1".repeat(40))),
            "G"
        )),
        "Location position '111111111111111111111111...' is not an integer"
    );
}

#[test]
fn a_position_is_read_as_vep_rs_reads_it() {
    for text in [
        "30425650",
        "3.042565e7",
        "30425650.0",
        "+30425650",
        "30425650.",
        " 30425650 ",
        "030425650",
    ] {
        assert_eq!(vep_rs_pos(text), Some(30425650), "{text:?}");
    }
    assert_eq!(vep_rs_pos("3e7"), Some(30000000));
    assert_eq!(VEP_RS_COORDINATE_CEILING, (1 << 62) - 1);
    assert_eq!(
        vep_rs_pos(&VEP_RS_COORDINATE_CEILING.to_string()),
        Some(VEP_RS_COORDINATE_CEILING)
    );
    assert_eq!(
        vep_rs_pos(&((1i128 << 56) + 1).to_string()),
        Some((1 << 56) + 1)
    );
    assert_eq!(
        vep_rs_pos(&format!(" {} ", (1i128 << 56) + 1)),
        Some(1 << 56)
    );
    let past = (VEP_RS_COORDINATE_CEILING + 1).to_string();
    let padded = format!(" {VEP_RS_COORDINATE_CEILING} ");
    for text in [
        ".",
        "",
        "abc",
        "-30425650",
        "30_425_650",
        "\u{ff13}\u{ff10}\u{ff14}\u{ff12}\u{ff15}\u{ff16}\u{ff15}\u{ff10}",
        "0x1D04A1A",
        "inf",
        "nan",
        "1e400",
        "30425650.5",
        "++30425650",
        past.as_str(),
        padded.as_str(),
    ] {
        assert_eq!(vep_rs_pos(text), None, "{text:?}");
    }
}

#[test]
fn a_list_is_split_as_perl_splits_it_and_a_span_is_written_as_both_engines_write_it() {
    assert_eq!(perl_split("<DEL>,", ","), ["<DEL>"]);
    assert_eq!(perl_split("G,,T,,", ","), ["G", "", "T"]);
    assert!(perl_split("", ",").is_empty());
    assert!(perl_split(",", ",").is_empty());
    assert_eq!(perl_split("N/G/<NON_REF>/", "/"), ["N", "G", "<NON_REF>"]);
    assert_eq!(sequence_location("21", 101, 100), "21:100-101");
    assert_eq!(sequence_location("21", 101, 102), "21:101-102");
    assert_eq!(sequence_location("chr21", 100, 100), "chr21:100");
}

#[test]
fn rows_under_an_identifier_several_input_lines_carry_find_their_record_by_place_and_alleles() {
    let two = [
        HEADER.to_string(),
        insertion_row_text("CAT", UTR_LOSS),
        insertion_row_text("G", UTR_NONE),
    ]
    .concat();
    let entries = parse_default(&two).unwrap();
    let own = UTR_INPUTS[0];
    let other_allele = "21\t25734954\tutr_ins\tA\tAG\t.\t.\t.";
    let handed = |records: &RecordRows, entry: &Row| -> Option<Vec<(String, Option<String>)>> {
        records.of(entry).unwrap().map(|rows| {
            rows.iter()
                .map(|e| (e.allele.clone(), e.fields.get("ALT").map(str::to_string)))
                .collect()
        })
    };
    let pair = |allele: &str, alt: Option<&str>| (allele.to_string(), alt.map(str::to_string));
    // The rows of two records that share the identifier and the Location are each handed
    // their own record alone; one line of two alternates keeps both of its rows.
    let two_records = RecordRows::new(&entries, &[own, other_allele]);
    for (entry, alternate) in entries.iter().zip(["ACAT", "AG"]) {
        assert_eq!(
            handed(&two_records, entry),
            Some(vec![pair(&entry.allele, None), pair("", Some(alternate))]),
            "{alternate}"
        );
    }
    let one_record = RecordRows::new(&entries, &["21\t25734954\tutr_ins\tA\tACAT,AG\t.\t.\t."]);
    assert_eq!(
        handed(&one_record, &entries[0]),
        Some(vec![
            pair("CAT", None),
            pair("G", None),
            pair("", Some("ACAT,AG"))
        ])
    );
    let twice = RecordRows::new(&entries, &[own, own]);
    assert_eq!(handed(&twice, &entries[0]), None);
    assert_eq!(twice.one_record(&entries[0]), Ok(false));
    // A line that joins the identifier with another stands under the first of them.
    let joined = own.replace("utr_ins", "utr_ins;rs2");
    let told_apart = RecordRows::new(&entries, &[joined.as_str(), other_allele]);
    for (entry, alternate) in entries.iter().zip(["ACAT", "AG"]) {
        assert_eq!(
            handed(&told_apart, entry),
            Some(vec![pair(&entry.allele, None), pair("", Some(alternate))]),
            "{alternate}"
        );
    }
    // The line's own entry is a new entry at every call.
    let first = one_record.of(&entries[0]).unwrap().unwrap();
    let second = one_record.of(&entries[0]).unwrap().unwrap();
    assert!(Rc::ptr_eq(&first[0], &entries[0]) && Rc::ptr_eq(&first[0], &second[0]));
    assert!(!Rc::ptr_eq(&first[2], &second[2]));
    assert_eq!(first[2].key, "21:25734954-25734955|record");
    assert_eq!(first[2].location.as_deref(), Some("21:25734954-25734955"));
}

#[test]
fn a_row_under_the_first_identifier_of_a_several_identifier_line_finds_its_record() {
    let reference = [HEADER.to_string(), insertion_row_text("CAT", UTR_LOSS)].concat();
    let entries = parse_default(&reference).unwrap();
    let line = UTR_INPUTS[0].replace("utr_ins", "utr_ins;rs2;rs3;rs4;rs5");
    let records = RecordRows::new(&entries, &[line]);
    assert_eq!(records.one_record(&entries[0]), Ok(true));
    let found = records.line_of(&entries[0]).unwrap().unwrap();
    assert_eq!(found[..5], ["21", "25734954", "utr_ins", "A", "ACAT"]);
    let handed: Vec<(String, Option<String>, Option<String>)> = records
        .of(&entries[0])
        .unwrap()
        .unwrap()
        .iter()
        .map(|e| {
            (
                e.allele.clone(),
                e.fields.get("Uploaded_variation").map(str::to_string),
                e.fields.get("ALT").map(str::to_string),
            )
        })
        .collect();
    assert_eq!(
        handed,
        [
            ("CAT".to_string(), Some("utr_ins".to_string()), None),
            (
                String::new(),
                Some("utr_ins".to_string()),
                Some("ACAT".to_string())
            )
        ]
    );
    // A line that carries the identifier as a later element, or behind a first element that
    // is a dot, empty or `0`, is no record of the row; nor is a line of fewer than five columns.
    for name in ["rs2;utr_ins", ".;utr_ins", ";utr_ins", "0;utr_ins"] {
        let records = RecordRows::new(&entries, &[UTR_INPUTS[0].replace("utr_ins", name)]);
        assert_eq!(records.one_record(&entries[0]), Ok(false), "{name}");
        assert!(records.of(&entries[0]).unwrap().is_none(), "{name}");
    }
    let short = RecordRows::new(&entries, &["21\t25734954\tutr_ins\tA"]);
    assert_eq!(short.one_record(&entries[0]), Ok(false));
    let none: [&str; 0] = [];
    assert!(RecordRows::new(&entries, &none)
        .of(&entries[0])
        .unwrap()
        .is_none());
}

#[test]
fn an_input_line_stands_under_the_first_identifier_of_its_id_column_alone() {
    // Under no later element and not under the whole column, and under no name where the
    // first element is a dot, empty or `0`.
    let names = [
        "rs1;rs2", "rs3", ".;rs4", ";rs5", "0;rs6", "0", "", ".", "rs1",
    ];
    let lines: Vec<String> = names
        .iter()
        .enumerate()
        .map(|(number, name)| format!("21\t{}\t{name}\tA\tT\t.\t.\t.", 100 + number))
        .collect();
    let row_at = |name: &str, position: usize| -> Row {
        Rc::new(Entry {
            key: format!("21:{position}|T|ENST1|Transcript"),
            location: Some(format!("21:{position}")),
            allele: "T".to_string(),
            feature: "ENST1".to_string(),
            fields: [("Uploaded_variation", name)].into_iter().collect(),
            ..Entry::default()
        })
    };
    let none: [Row; 0] = [];
    let records = RecordRows::new(&none, &lines);
    let found = |name: &str, position: usize| -> Option<Vec<String>> {
        records
            .line_of(&row_at(name, position))
            .unwrap()
            .map(|line| line[..3].to_vec())
    };
    let line = |position: &str, name: &str| {
        Some(vec![
            "21".to_string(),
            position.to_string(),
            name.to_string(),
        ])
    };
    assert_eq!(found("rs3", 5), line("101", "rs3"));
    assert_eq!(found("rs1", 100), line("100", "rs1"));
    assert_eq!(found("rs1", 108), line("108", "rs1"));
    assert_eq!(found("rs1", 104), None);
    for name in ["rs2", "rs4", "rs5", "rs6", "0", "", ".", "rs1;rs2", ".;rs4"] {
        assert_eq!(found(name, 100), None, "{name:?}");
    }
}

#[test]
fn each_row_under_several_lines_is_handed_the_rows_its_own_line_writes() {
    // Rows on four features at one Location: of three one-allele records under one
    // identifier, of three alleles under a second identifier whose three lines write the
    // first once, the second twice and the third never, and of two alleles under a third
    // that one line carries.
    let features = ["ENST1", "ENST2", "ENST3", "ENST4"];
    let carried = [
        ("rs1", "C"),
        ("rs1", "G"),
        ("rs1", "T"),
        ("rs2", "C"),
        ("rs2", "G"),
        ("rs2", "G"),
        ("rs3", "T"),
    ];
    let lines: Vec<String> = carried
        .iter()
        .map(|(name, alternate)| format!("21\t1000\t{name}\tA\t{alternate}\t.\t.\t."))
        .collect();
    let written = [
        ("rs1", "C"),
        ("rs2", "G"),
        ("rs1", "G"),
        ("rs3", "C"),
        ("rs2", "C"),
        ("rs1", "T"),
        ("rs2", "T"),
        ("rs3", "T"),
    ];
    let mut entries: Vec<Row> = Vec::new();
    for feature in features {
        for (name, allele) in written {
            entries.push(Rc::new(Entry {
                key: format!("21:1000|{allele}|{feature}|Transcript"),
                location: Some("21:1000".to_string()),
                allele: allele.to_string(),
                feature: feature.to_string(),
                fields: [("Uploaded_variation", name)].into_iter().collect(),
                ..Entry::default()
            }));
        }
    }
    let records = RecordRows::new(&entries, &lines);
    let handed = |e: &Row| -> Option<Vec<(String, String, Option<String>)>> {
        records.of(e).unwrap().map(|rows| {
            rows.iter()
                .map(|row| {
                    (
                        row.allele.clone(),
                        row.feature.clone(),
                        row.fields.get("ALT").map(str::to_string),
                    )
                })
                .collect()
        })
    };
    let expected =
        |alleles: &[&str], alternate: &str| -> Option<Vec<(String, String, Option<String>)>> {
            let mut rows: Vec<(String, String, Option<String>)> = Vec::new();
            for feature in features {
                for allele in alleles {
                    rows.push((allele.to_string(), feature.to_string(), None));
                }
            }
            rows.push((String::new(), String::new(), Some(alternate.to_string())));
            Some(rows)
        };
    for _ in 0..2 {
        for (number, e) in entries.iter().enumerate() {
            let found = handed(e);
            let wanted = match written[number % written.len()] {
                ("rs1", allele) => expected(&[allele], allele),
                ("rs2", "C") => expected(&["C"], "C"),
                ("rs2", _) => None,
                (_, _) => expected(&["C", "T"], "T"),
            };
            assert_eq!(found, wanted, "row {number}");
        }
    }
    assert_eq!(
        line_writes_row(
            &[
                "21".to_string(),
                "1000".to_string(),
                "rs1".to_string(),
                "A".to_string(),
                "C".to_string()
            ],
            &entries[0]
        ),
        Ok(true)
    );
    assert_eq!(
        line_writes_row(
            &[
                "chr21".to_string(),
                "1000".to_string(),
                "rs1".to_string(),
                "A".to_string(),
                "G,C".to_string()
            ],
            &entries[0]
        ),
        Ok(true)
    );
    assert_eq!(
        line_writes_row(
            &[
                "22".to_string(),
                "1000".to_string(),
                "rs1".to_string(),
                "A".to_string(),
                "C".to_string()
            ],
            &entries[0]
        ),
        Ok(false)
    );
    assert_eq!(
        line_writes_row(
            &[
                "21".to_string(),
                "x".to_string(),
                "rs1".to_string(),
                "A".to_string(),
                "C".to_string()
            ],
            &entries[0]
        ),
        Ok(false)
    );
}

#[test]
fn the_rows_of_a_record_are_those_under_its_label() {
    let rows = rows();
    let inputs = [
        "21\t200\trs9\tGACC\tG,GAC\t.\t.\t.",
        "X\t300\trs8;rs7\tA\tAACG\t.\t.\t.",
        "21\t100\t.\tA\tAT,G",
    ];
    let records = RecordRows::new(&rows, &inputs);
    let of = |number: usize| records.of(&rows[number]).unwrap().unwrap();
    assert_eq!(keys(&of(0)), ["21:100-101|T|ENST1|Transcript"]);
    assert_eq!(
        keys(&of(1)),
        ["21:100|G|ENST1|Transcript", "21:100|G|ENST2|Transcript"]
    );
    assert_eq!(
        keys(&of(3)),
        [
            "21:200-202|-|ENST1|Transcript",
            "21:200-202|AC|ENST1|Transcript",
            "21:200-202|record"
        ]
    );
    assert_eq!(keys(&of(5)), ["21:200|deletion|ENST1|-", "21:200|record"]);
    assert_eq!(
        keys(&of(7)),
        [
            "X:300-301|ACG|-|-",
            "X:300-301|ACG|-|-#2",
            "X:300-301|record"
        ]
    );
    let record: Vec<(String, String)> = of(6)[2]
        .fields
        .iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect();
    let named = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect()
    };
    assert_eq!(
        record,
        named(&[
            ("Uploaded_variation", "rs8"),
            ("CHROM", "X"),
            ("POS", "300"),
            ("REF", "A"),
            ("ALT", "AACG")
        ])
    );
    for row in &rows {
        assert_eq!(records.one_record(row), Ok(true));
    }
    // What the helpers read of a record's rows, without them and with them.
    let alleles = |e: &Row, record: &[Row]| {
        record_alleles(e, record)
            .into_iter()
            .collect::<Vec<_>>()
            .join(" ")
    };
    assert_eq!(alleles(&rows[3], &rows), "- AC");
    assert_eq!(alleles(&rows[5], &rows), "deletion");
    assert_eq!(alleles(&rows[3], &of(3)), "- AC");
    assert_eq!(alleles(&rows[0], &rows), "T");
    assert_eq!(
        keys(&sibling_rows(&rows[3], &rows)),
        ["21:200-202|AC|ENST1|Transcript"]
    );
    assert_eq!(
        keys(&sibling_rows(&rows[4], &rows)),
        ["21:200-202|-|ENST1|Transcript"]
    );
    assert!(sibling_rows(&rows[1], &rows).is_empty() && sibling_rows(&rows[6], &rows).is_empty());
    let chromosomes: Vec<Option<String>> = rows
        .iter()
        .map(|e| record_chromosome(e, Some(&rows)))
        .collect();
    assert_eq!(chromosomes[0].as_deref(), Some("21"));
    assert_eq!(chromosomes[7].as_deref(), Some("X"));
    let shapes: Vec<bool> = rows.iter().map(|e| indel_row(e).unwrap()).collect();
    assert_eq!(shapes, [true, false, false, true, true, false, true, true]);
    let from_record: Vec<bool> = rows
        .iter()
        .map(|e| record_indel_row(e, &rows).unwrap())
        .collect();
    assert_eq!(from_record, shapes);
    let with_lines: Vec<bool> = (0..rows.len())
        .map(|number| record_indel_row(&rows[number], &of(number)).unwrap())
        .collect();
    assert_eq!(
        with_lines,
        [true, false, false, true, true, false, true, true]
    );
    let insertions: Vec<bool> = rows.iter().map(|e| insertion_row(e)).collect();
    assert_eq!(
        insertions,
        [true, false, false, false, false, false, true, true]
    );
    let spans: Vec<Option<i128>> = rows.iter().map(|e| location_span(e).unwrap()).collect();
    assert_eq!(
        spans,
        [
            Some(2),
            Some(1),
            Some(1),
            Some(3),
            Some(3),
            Some(1),
            Some(2),
            Some(2)
        ]
    );
    let lengths: Vec<Option<(i128, i128)>> = rows
        .iter()
        .map(|e| record_allele_lengths(e, Some(&rows)))
        .collect();
    assert_eq!(
        lengths,
        [
            None,
            Some((1, 1)),
            Some((1, 1)),
            Some((3, 0)),
            Some((3, 2)),
            None,
            None,
            None
        ]
    );
    let lengths_with_lines: Vec<Option<(i128, i128)>> = (0..rows.len())
        .map(|number| record_allele_lengths(&rows[number], Some(&of(number))))
        .collect();
    assert_eq!(
        lengths_with_lines,
        [
            None,
            Some((1, 1)),
            Some((1, 1)),
            Some((3, 0)),
            Some((3, 2)),
            None,
            Some((0, 3)),
            Some((0, 3))
        ]
    );
    let firsts: Vec<Option<String>> = rows.iter().map(|e| first_alternate(e, &rows)).collect();
    assert_eq!(
        firsts[..3],
        [
            Some("AT".to_string()),
            Some("AT".to_string()),
            Some("AT".to_string())
        ]
    );
    assert!(firsts[3..].iter().all(Option::is_none));
    assert_eq!(first_alternate(&rows[3], &of(3)).as_deref(), Some("-"));
    assert_eq!(first_alternate(&rows[6], &of(6)).as_deref(), Some("ACG"));
    let on_strand: Vec<Option<String>> = rows
        .iter()
        .map(|e| allele_on_transcript_strand(e))
        .collect();
    assert_eq!(
        on_strand[..3],
        [
            Some("A".to_string()),
            Some("G".to_string()),
            Some("G".to_string())
        ]
    );
    assert!(on_strand[3..].iter().all(Option::is_none));
    let shifts: Vec<bool> = rows.iter().map(|e| positive_shift(e).unwrap()).collect();
    assert_eq!(
        shifts,
        [true, false, false, true, false, false, false, false]
    );
    let cdna: Vec<Option<(i128, i128)>> = rows.iter().map(|e| cdna_span(e)).collect();
    assert_eq!(
        cdna,
        [
            None,
            None,
            None,
            None,
            Some((10, 12)),
            Some((7, 7)),
            None,
            None
        ]
    );
    let closed: Vec<bool> = rows.iter().map(|e| cdna_span_closed(e)).collect();
    assert_eq!(
        closed,
        [false, false, false, false, true, false, false, false]
    );
    let strands: Vec<&str> = rows.iter().map(|e| field_of(e, &STRAND_FIELDS)).collect();
    assert_eq!(
        strands,
        ["-1", "1", "1", "-1", ABSENT, ABSENT, ABSENT, ABSENT]
    );
    let by_prefix: Vec<&str> = rows.iter().map(|e| field_of(e, &["shift_"])).collect();
    assert_eq!(
        by_prefix,
        ["3", "0", ABSENT, "2", ABSENT, ABSENT, ABSENT, ABSENT]
    );
    assert_eq!(field_of_any(None, &["x"]), ABSENT);
    assert_eq!(
        field_of_any(Some(&rows[0]), &["hgvsc"]),
        "ENST1.1:c.5_6insT"
    );
    assert!(Rc::ptr_eq(&at_its_location(&rows[0]), &rows[0]));
    let placed = at_its_location(&Rc::new(Entry {
        allele_location: Some("1:5".to_string()),
        ..Entry::default()
    }));
    assert_eq!(placed.location.as_deref(), Some("1:5"));
}

#[test]
fn the_lengths_ensembl_vep_hands_a_plugin_are_read_from_the_record_or_the_row() {
    let text = [
        HEADER.to_string(),
        utr_row("21_100_AC/ACGT", "21:101", "CGT", "ENST1", "7", UTR_NONE),
        utr_row("del", "21:1000-1004", "-", "ENST2", "336-340", UTR_NONE),
        insertion_row_text("CAT", UTR_NONE),
        utr_row("rs1", "21:700-702", "A", "ENST1", "7-9", UTR_NONE),
    ]
    .concat();
    let entries = parse_default(&text).unwrap();
    // A record's shared first base is dropped, a deletion's Location and a Location of one
    // base or of three or more state the reference allele, and a two-base Location under an
    // identifier states none.
    let lengths: Vec<Option<(i128, i128)>> = entries
        .iter()
        .map(|e| record_allele_lengths(e, None))
        .collect();
    assert_eq!(lengths, [Some((1, 3)), Some((5, 0)), None, Some((3, 1))]);
    let record = RecordRows::new(&entries, &UTR_INPUTS)
        .of(&entries[2])
        .unwrap()
        .unwrap();
    assert_eq!(
        record_allele_lengths(&entries[2], Some(&record)),
        Some((0, 3))
    );
    let spelt = |reference: &str, alternates: &[&str]| {
        let (reference, alternates) = vcf_alleles_as_spelt(reference, alternates);
        format!("{reference} {}", alternates.join(","))
    };
    assert_eq!(spelt("A", &["AT"]), "- T");
    assert_eq!(spelt("AC", &["AT", "AG"]), "AC AT,AG");
    assert_eq!(spelt("A", &["AT", "*"]), "- T,*");
    assert_eq!(spelt("A", &["T", "<DEL:*>", "AG"]), "A T,<DEL:*>,AG");
    assert_eq!(spelt("", &["A"]), " A");
    assert_eq!(spelt("A", &["A"]), "A A");
    assert_eq!(spelt("AT", &["A", "ATT"]), "T -,TT");
    assert_eq!(spelt("A", &[""]), "A ");
    assert_eq!(
        vcf_first_alternate("AT", &["A", "ATT"]).as_deref(),
        Some("-")
    );
    let none: [&str; 0] = [];
    assert_eq!(vcf_first_alternate("AT", &none), None);
    let spans: Vec<Option<(i128, i128)>> = entries.iter().map(|e| cdna_span(e)).collect();
    assert_eq!(
        spans,
        [Some((7, 7)), Some((336, 340)), Some((64, 65)), Some((7, 9))]
    );
    let closed: Vec<bool> = entries.iter().map(|e| cdna_span_closed(e)).collect();
    assert_eq!(closed, [false, true, true, true]);
    assert_eq!(
        entries
            .iter()
            .map(|e| (e.indel, e.insertion))
            .collect::<Vec<_>>(),
        [
            (Some(true), Some(false)),
            (None, None),
            (None, None),
            (None, None)
        ]
    );
    assert_eq!(entries[0].reference.borrow().as_deref(), Some("C"));
    let insertions: Vec<bool> = entries.iter().map(|e| insertion_row(e)).collect();
    assert_eq!(insertions, [true, false, true, false]);
    let spans: Vec<Option<i128>> = entries.iter().map(|e| location_span(e).unwrap()).collect();
    assert_eq!(spans, [Some(1), Some(5), Some(2), Some(3)]);
}

#[test]
fn notations_lists_numbers_and_field_names_are_read_by_their_helpers() {
    for (hgvsc, length) in [
        ("ENST1.1:c.5_6insT", Some(1)),
        ("ENST1.1:c.5_6insTTA", Some(3)),
        ("c.5dup", None),
        ("x:c.5del", Some(1)),
        ("x:c.5_9del", Some(5)),
        ("x:c.-3_2dup", Some(5)),
        ("x:c.*5_*9del", None),
        ("x:c.5+1_9del", None),
        ("x:n.5_9dup", Some(5)),
        ("x:c.5_9delins", None),
        ("x:c.5A>T", None),
        ("x:c.-9_-5del", Some(5)),
        ("nocolon", None),
        ("x:c.\u{665}_\u{669}del", Some(5)),
        ("x:c.5_6insTN", None),
        // The length is the smallest `i128` where the positions' difference lies below it.
        (
            "x:c.1_-170141183460469231731687303715884105728del",
            Some(i128::MIN),
        ),
        ("x:c.2_-170141183460469231731687303715884105728del", None),
    ] {
        assert_eq!(notation_indel_length(hgvsc), length, "{hgvsc}");
    }
    assert_eq!(list_parts("a,b", "b,a"), (vec!["a", "b"], vec!["b", "a"]));
    assert_eq!(list_parts("a&b", "b&a"), (vec!["a", "b"], vec!["b", "a"]));
    assert_eq!(list_parts("a,b&c", "x"), (vec![], vec![]));
    assert_eq!(list_parts("a", "b,c&d"), (vec![], vec![]));
    assert_eq!(list_parts("", ""), (vec![""], vec![""]));
    assert_eq!(list_parts("a", "b&c"), (vec!["a"], vec!["b&c"]));
    for (allele, complement) in [
        ("ACGTN", "NACGT"),
        ("acgtn", "nacgt"),
        ("A-C", "G-T"),
        ("", ""),
        ("deletion", "deletion"),
        ("ACX", "ACX"),
        ("-", "-"),
    ] {
        assert_eq!(reverse_complement(allele), complement, "{allele}");
    }
    for (value, number) in [
        ("1", Some(1.0)),
        ("1.5", Some(1.5)),
        (" 2 ", Some(2.0)),
        ("abc", None),
        ("nan", None),
        ("inf", None),
        ("-", None),
        ("", None),
        ("1e3", Some(1000.0)),
        ("1_0", Some(10.0)),
        (ABSENT, None),
    ] {
        assert_eq!(as_number(value), number, "{value:?}");
    }
    for (name, names, named) in [
        ("HGVSc", &["hgvsc"][..], true),
        ("gnomAD_AF", &["gnomad_"][..], true),
        ("gnomAD", &["gnomad_"][..], false),
        ("SpliceAI_pred", &["SpliceAI_", "spliceai"][..], true),
        ("spliceai", &["SpliceAI_", "spliceai"][..], true),
        ("5UTR_annotation", &["5utr_"][..], true),
        ("am_class", &["AM_"][..], true),
        ("\u{212a}ELVIN", &["kelvin"][..], true),
        ("\u{17f}", &["S"][..], false),
        ("\u{130}x", &["i\u{307}x"][..], true),
    ] {
        assert_eq!(field_named(name, names), named, "{name}");
    }
    for (feature_type, kind) in [
        ("Transcript", "transcript"),
        ("RegulatoryFeature", "regulatory_feature"),
        ("MotifFeature", "motif_feature"),
        ("", "intergenic"),
        ("-", "intergenic"),
        ("regulatory_feature", "regulatory_feature"),
        ("intergenic", "intergenic"),
        ("Feature", "_feature"),
    ] {
        assert_eq!(feature_kind(feature_type), kind, "{feature_type:?}");
    }
    assert_eq!(decode_extra("a%3Db%3Bc"), "a=b;c");
    assert_eq!(decode_extra("%3d"), "%3d");
    assert_eq!(decode_extra("plain"), "plain");
    assert_eq!(decode_extra("%3D%3B%3D"), "=;=");
    assert_eq!(integer_of("12", "x"), Ok(12));
    assert_eq!(
        integer_of("111111111111111111111111111111", "Location position"),
        Ok(111111111111111111111111111111)
    );
    assert_eq!(
        refusal(integer_of("\u{b2}", "SHIFT_LENGTH")),
        "SHIFT_LENGTH '\u{b2}' is not an integer"
    );
    assert_eq!(refusal(integer_of("", "p")), "p '' is not an integer");
    // A number of more digits than `int` converts is refused whatever its value.
    assert_eq!(integer_of(&"0".repeat(4300), "x"), Ok(0));
    assert_eq!(
        refusal(integer_of(&"0".repeat(4301), "Location position")),
        "Location position '000000000000000000000000...' is not an integer"
    );
    assert!(
        is_integer("-12")
            && is_integer("\u{665}")
            && !is_integer("")
            && !is_integer("-")
            && !is_integer("+1")
            && !is_integer("1.0")
    );
    assert!(
        is_non_zero_integer("-12")
            && is_non_zero_integer("1\u{665}")
            && !is_non_zero_integer("0")
            && !is_non_zero_integer("012")
            && !is_non_zero_integer("\u{661}")
    );
    let shifted = |value: &str| {
        let fields: Fields = [("shift_length", value)].into_iter().collect();
        positive_shift(&Entry {
            fields,
            ..Entry::default()
        })
    };
    assert_eq!(shifted("2"), Ok(true));
    assert_eq!(shifted("0"), Ok(false));
    assert_eq!(shifted("-1"), Ok(false));
    assert_eq!(
        refusal(shifted("\u{b2}")),
        "SHIFT_LENGTH '\u{b2}' is not an integer"
    );
}

#[test]
fn list_valued_texts_agree_as_multisets_under_the_fields_separator() {
    assert!(same_items("a,b", "b,a", ",") && !same_items("a,b", "a,a,b", ","));
    assert!(same_items("-", "-", ",") && !same_items("-", "a", ","));
    assert!(!same_items("", "-", ",") && same_items(ABSENT, ABSENT, ","));
    assert!(same_items(
        r#"{"db":"x","name":"1"},{"db":"y","name":"2"}"#,
        r#"{"db":"y","name":"2"},{"db":"x","name":"1"}"#,
        JSON_ARRAY
    ));
    assert!(!same_items("a,b", "b,a", JSON_ARRAY));
    assert!(same_items("1, 2", "2,1", JSON_ARRAY));
    assert!(same_items(
        r#""a",{"b":1,"a":2}"#,
        r#"{"a":2,"b":1},"a""#,
        JSON_ARRAY
    ));
    assert!(!same_items("1,2", "2,1,", JSON_ARRAY) && !same_items("1.0", "1", JSON_ARRAY));
    assert!(same_items("k=1:j=2,x=1", "x=1,j=2:k=1", "pairs:,"));
    assert!(!same_items("k=1:j=2,x=1", "x=1,j=2:k=3", "pairs:,"));
    assert!(same_items("a&b", "b&a", "pairs:&"));
    let listed =
        |fmt: &str| -> Vec<(String, String)> { order_only_fields(fmt).into_iter().collect() };
    let named = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect()
    };
    assert_eq!(
        listed("default"),
        named(&[
            ("5UTR_annotation", "pairs:,"),
            ("DOMAINS", ","),
            ("FREQS", ","),
            ("NearestExonJB", ","),
            ("ZYG", ",")
        ])
    );
    assert_eq!(listed("tab"), listed("default"));
    assert_eq!(
        listed("vcf"),
        named(&[
            ("5UTR_annotation", "pairs:&"),
            ("DOMAINS", "&"),
            ("FREQS", "&"),
            ("NearestExonJB", "&"),
            ("ZYG", "&")
        ])
    );
    assert_eq!(
        listed("json"),
        named(&[("domains", "json"), ("nearestexonjb", ","), ("zyg", ",")])
    );
}

/// One consult of the built-in rules over a row of no fields.
fn accepts(klass: &str, name: &str, expected: &str, actual: &str) -> Result<bool, Refusal> {
    let e: Row = Rc::new(Entry::default());
    let declared = BTreeSet::new();
    BuiltIn.accepts(&Consult {
        klass,
        e: &e,
        name,
        expected,
        actual,
        a: None,
        record_rows: None,
        candidate_rows: None,
        declared: &declared,
        assembly: None,
        source: None,
        gwas_pairs: None,
        revel_table: None,
        utrannotator_sweep: None,
        downstream_sequences: None,
        splice_windows: None,
        spliceai_table: None,
        exon_spans: None,
    })
}

#[test]
fn the_classes_of_one_test_accept_their_shape_and_nothing_else() {
    // A value the engine alone writes.
    for (klass, field) in [
        ("alphamissense_two_digit_transcript_version", "am_class"),
        (
            "alphamissense_two_digit_transcript_version",
            "AM_PATHOGENICITY",
        ),
        ("splicevault_json_later_rows", "splicevault"),
    ] {
        assert_eq!(
            accepts(klass, field, ABSENT, "likely_benign"),
            Ok(true),
            "{klass}"
        );
        assert_eq!(accepts(klass, field, "-", "0.1"), Ok(true), "{klass}");
        assert_eq!(accepts(klass, field, "", "0.1"), Ok(true), "{klass}");
        assert_eq!(accepts(klass, field, "0.2", "0.1"), Ok(false), "{klass}");
        assert_eq!(accepts(klass, field, "0.2", ABSENT), Ok(false), "{klass}");
        assert_eq!(accepts(klass, field, ABSENT, ""), Ok(false), "{klass}");
        assert_eq!(accepts(klass, "REVEL", ABSENT, "0.1"), Ok(false), "{klass}");
    }
    assert_eq!(
        accepts(
            "splicevault_json_later_rows",
            "SpliceVault_site_type",
            ABSENT,
            "donor"
        ),
        Ok(false)
    );
    // A minimum the engine writes as zero where the reference writes none.
    assert_eq!(
        accepts("loeuf_zero_minimum", "LOEUF", ABSENT, "0"),
        Ok(true)
    );
    assert_eq!(accepts("loeuf_zero_minimum", "loeuf", "-", "0"), Ok(true));
    assert_eq!(
        accepts("loeuf_zero_minimum", "LOEUF", ABSENT, "0.0"),
        Ok(false)
    );
    assert_eq!(
        accepts("loeuf_zero_minimum", "LOEUF", "0.5", "0"),
        Ok(false)
    );
    assert_eq!(
        accepts("loeuf_zero_minimum", "LOEUF", ABSENT, ABSENT),
        Ok(false)
    );
    assert_eq!(
        accepts("loeuf_zero_minimum", "LOEUF_bin", ABSENT, "0"),
        Ok(false)
    );
    // A value the reference cuts at its second `=`.
    let klass = "gnomadmt_info_value_second_equals";
    assert_eq!(
        accepts(klass, "gnomAD_hap_defining_variant", "H", "H=1"),
        Ok(true)
    );
    assert_eq!(accepts(klass, "gnomAD_x", "a%3Db", "a%3Db%3Dc"), Ok(false));
    assert_eq!(accepts(klass, "gnomAD_x", "a", "a%3Db=c"), Ok(true));
    assert_eq!(accepts(klass, "gnomAD_x", "a%3Bb", "a;b=c"), Ok(true));
    assert_eq!(accepts(klass, "gnomAD_x", "H", "H"), Ok(false));
    assert_eq!(accepts(klass, "gnomAD_x", "H", "G=1"), Ok(false));
    assert_eq!(accepts(klass, "gnomAD_x", ABSENT, "=1"), Ok(false));
    assert_eq!(accepts(klass, "gnomAD_x", "H", ABSENT), Ok(false));
    assert_eq!(accepts(klass, "gnomAD", "H", "H=1"), Ok(false));
    assert_eq!(accepts(klass, "AF", "H", "H=1"), Ok(false));
}

#[test]
fn a_class_without_a_rule_is_refused_once_it_reads_the_field() {
    assert_eq!(
        accepts("hgnc_id_copy", "am_class", "HGNC:1", "HGNC:2"),
        Ok(false)
    );
    assert_eq!(
        refusal(accepts("hgnc_id_copy", "HGNC_ID", "HGNC:1", "HGNC:2")),
        "hgnc_id_copy: no rule of this build"
    );
    assert_eq!(
        refusal(accepts("no_such_class", "REVEL", "0.5", "-")),
        "no_such_class: no rule of this build"
    );
    assert_eq!(
        accepts(
            "plugin_value_on_reference_skipped_record",
            "pHaplo",
            ABSENT,
            "0.9"
        ),
        Ok(false)
    );
    assert!(BuiltIn.knows("hgnc_id_copy") && BuiltIn.knows("loeuf_zero_minimum"));
    assert!(!BuiltIn.knows("order_only") && !BuiltIn.knows("no_such_class"));
    assert!(BuiltIn.transcript_reading("plugin_hgvs_shifted_cds_coordinates"));
    assert!(BuiltIn.sibling_allele("splicevault_sibling_allele_rows"));
    assert!(BuiltIn.record_reference("plugin_value_follows_splice_term_past_the_record_span"));
    assert!(!BuiltIn.record_reference("revel_sibling_allele_rows"));
}

struct Scratch(std::path::PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-rules-{}-{tag}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn file(&self, name: &str, content: &[u8]) -> std::path::PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, content).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn a_file_is_read_as_text_with_one_line_ending_and_refused_with_its_reason() {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    let scratch = Scratch::new("read-text");
    let plain = scratch.file("plain.txt", b"a\r\nb\rc\nd\xc3\xa9");
    assert_eq!(read_text(&plain).unwrap(), "a\nb\nc\nd\u{e9}");
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(b"x\ty\r\nz\n").unwrap();
    let packed = encoder.finish().unwrap();
    let gz = scratch.file("packed.txt.gz", &packed);
    assert_eq!(read_text(&gz).unwrap(), "x\ty\nz\n");
    assert_eq!(read_text(&scratch.file("empty.gz", b"")).unwrap(), "");
    let reason = |path: &Path| {
        let message = refusal(read_text(path));
        message
            .strip_prefix(&format!("cannot read {}: ", path.display()))
            .unwrap()
            .to_string()
    };
    assert_eq!(
        reason(&scratch.0.join("absent.txt")),
        "No such file or directory"
    );
    assert_eq!(reason(&scratch.0), "Is a directory");
    assert_eq!(
        reason(&scratch.file("plain.gz", b"ab\n")),
        "Not a gzipped file (b'ab')"
    );
    assert_eq!(
        reason(&scratch.file("odd.gz", b"\x00'")),
        "Not a gzipped file (b\"\\x00'\")"
    );
    assert_eq!(
        reason(&scratch.file("one.gz", b"\x1f")),
        "Not a gzipped file (b'\\x1f')"
    );
    assert_eq!(
        reason(&scratch.file("both.gz", b"'\"")),
        "Not a gzipped file (b'\\'\"')"
    );
    assert_eq!(
        reason(&scratch.file("cut.gz", &packed[..packed.len() - 9])),
        "Compressed file ended before the end-of-stream marker was reached"
    );
    for (bytes, message) in [
        (
            &b"ab\xff"[..],
            "'utf-8' codec can't decode byte 0xff in position 2: invalid start byte",
        ),
        (
            &b"\xc0\x80"[..],
            "'utf-8' codec can't decode byte 0xc0 in position 0: invalid start byte",
        ),
        (
            &b"a\xc3"[..],
            "'utf-8' codec can't decode byte 0xc3 in position 1: unexpected end of data",
        ),
        (
            &b"\xc3("[..],
            "'utf-8' codec can't decode byte 0xc3 in position 0: invalid continuation byte",
        ),
        (
            &b"\xe2\x82"[..],
            "'utf-8' codec can't decode bytes in position 0-1: unexpected end of data",
        ),
        (
            &b"\xe2\x82A"[..],
            "'utf-8' codec can't decode bytes in position 0-1: invalid continuation byte",
        ),
        (
            &b"\xe0\x80\x80"[..],
            "'utf-8' codec can't decode byte 0xe0 in position 0: invalid continuation byte",
        ),
        (
            &b"\xed\xa0\x80"[..],
            "'utf-8' codec can't decode byte 0xed in position 0: invalid continuation byte",
        ),
        (
            &b"\xf0\x90\x80A"[..],
            "'utf-8' codec can't decode bytes in position 0-2: invalid continuation byte",
        ),
        (
            &b"\xf0\x90"[..],
            "'utf-8' codec can't decode bytes in position 0-1: unexpected end of data",
        ),
        (
            &b"\xf4\x90\x80\x80"[..],
            "'utf-8' codec can't decode byte 0xf4 in position 0: invalid continuation byte",
        ),
        (
            &b"xx\xf0(\x8c\xbc"[..],
            "'utf-8' codec can't decode byte 0xf0 in position 2: invalid continuation byte",
        ),
    ] {
        assert_eq!(
            reason(&scratch.file("bytes.txt", bytes)),
            message,
            "{bytes:?}"
        );
    }
}

#[test]
fn a_cache_shard_is_an_array_of_transcripts_or_a_map_of_arrays() {
    let scratch = Scratch::new("shard");
    let ids = |path: &Path| -> Result<Vec<String>, Refusal> {
        let transcripts = shard_transcripts(path, read_shard(path)?)?;
        Ok(transcripts
            .iter()
            .map(|t| {
                format!(
                    "{:?} {:?}",
                    cache_integer(t.get("strand")),
                    cache_integer(t.get("start"))
                )
            })
            .collect())
    };
    let bare = scratch.file(
        "bare.json",
        br#"[{"strand": -1, "start": "100"}, null, {"strand": "1", "start": 1.5}]"#,
    );
    assert_eq!(ids(&bare).unwrap(), ["Some(-1) Some(100)", "Some(1) None"]);
    let keyed = scratch.file(
        "keyed.json",
        "{\"1-1000\": [{\"strand\": true}], \"1001-2000\": [null, {\"start\": \"\u{665}\"}]}"
            .as_bytes(),
    );
    assert_eq!(ids(&keyed).unwrap(), ["None None", "None Some(5)"]);
    let shown = |path: &Path| refusal(ids(path)).replace(&path.display().to_string(), "<shard>");
    assert_eq!(
        shown(&scratch.file("map.json", br#"{"1-1000": {"a": 1}}"#)),
        "<shard>: a range-keyed map value is not a transcript array"
    );
    assert_eq!(
        shown(&scratch.file("number.json", b"5")),
        "<shard> is neither a transcript array nor a range-keyed map"
    );
    assert_eq!(
        shown(&scratch.file("items.json", b"[{}, 5]")),
        "<shard>: a transcript is not an object"
    );
    assert!(shown(&scratch.file("broken.json", b"[{")).starts_with("<shard> is not JSON: "));
}
