// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the intron classes: the gaps and the two-window positions of the cache
//! transcripts under `tests/fixtures/cache_shards/`, the forms of an HGVSc body, the gaps
//! an allele names inside a record's span, and the three classes on rows made by hand.

use super::super::fields::{classify_field_divergence, FieldClass};
use super::*;
use crate::json::loads;
use crate::TempDir;

/// Eighteen transcripts of chromosome 1. Transcripts 1 to 14 (forward) and 21 (reverse) have
/// three exons of 100 bases and code from their 51st base to their 250th, so the second
/// exon opens on coding base 51 behind an intron of 100 bases (forward) or 200 (reverse).
/// Transcript 22 codes from its 151st base to its 250th: its second exon opens in the 5' UTR
/// and its fourth in the 3' UTR. Transcripts 23 and 24 code over their whole length behind
/// introns of 13, 2 and 14 bases, the fourth exon opening on coding base 286. The exons of
/// transcript `n` lie past `n` times 10,000.
const INTRON_CLASS_TRANSCRIPTS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/cache_shards/intron_classes.json"
));

/// Three transcripts of one exon layout: an exon of three bases, 201 to 203, between an
/// intron from 101 to 200 and one from 204 to 300, on the forward strand (ENST1) and on the
/// reverse strand (ENST2); ENST3 has one intron from 101 to 300.
const MAXENTSCAN_TRANSCRIPTS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/cache_shards/maxentscan.json"
));

const PROTEIN: &str = "ENSP00000000001.1:p.Gly17ArgfsTer5";
const TERMS: &str = "frameshift_variant,splice_region_variant";

fn object(text: &str) -> Object {
    match loads(text) {
        Ok(Value::Object(object)) => object,
        other => panic!("no JSON object: {other:?}"),
    }
}

/// The transcripts of a shard with their numbers, read off the digits of their stable ids.
fn transcripts(shard: &str) -> Vec<(Position, Object)> {
    let Ok(Value::List(records)) = loads(shard) else {
        panic!("no JSON list");
    };
    records
        .into_iter()
        .map(|record| {
            let Value::Object(transcript) = record else {
                panic!("no JSON object");
            };
            let stable_id = transcript.get("stable_id").and_then(Value::as_str).unwrap();
            let number = stable_id.trim_start_matches("ENST").parse().unwrap();
            (number, transcript)
        })
        .collect()
}

fn set<const N: usize>(positions: [Position; N]) -> BTreeSet<Position> {
    BTreeSet::from(positions)
}

/// A scratch cache that holds `shard` as its one shard of chromosome 1.
fn cache(tag: &str, shard: &str) -> (TempDir, String) {
    let scratch =
        TempDir::create(&std::env::temp_dir(), &format!("golden_introns_{tag}_")).unwrap();
    let cache = scratch.path.join("json_cache");
    let directory = cache.join("transcripts").join("1");
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(directory.join("1-1000000.json"), shard).unwrap();
    let cache = cache.to_str().unwrap().to_string();
    (scratch, cache)
}

fn key(location: &str, allele: &str, feature: &str, terms: &str) -> FieldKey {
    FieldKey {
        location: location.to_string(),
        allele: allele.to_string(),
        feature: feature.to_string(),
        feature_type: "Transcript".to_string(),
        consequences: terms.to_string(),
    }
}

fn row(fields: &[(&str, &str)]) -> FieldRow {
    fields
        .iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

/// `fields` with `changed` set over them.
fn with(fields: &[(&str, &str)], changed: &[(&str, &str)]) -> FieldRow {
    let mut changed_row = row(fields);
    changed_row.extend(row(changed));
    changed_row
}

/// The class `own` names for a question, `None` where it passes; the whole dispatch names
/// the same class.
fn named(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
    own: FieldClass,
) -> Option<&'static str> {
    let class = match own(question, models) {
        Verdict::Class(class, None) => Some(class),
        Verdict::Pass => None,
        other => panic!("{other:?}"),
    };
    let dispatched = classify_field_divergence(question, models);
    assert_eq!(dispatched.map(|(class, _)| class), class);
    class
}

/// The class named for the field `name` of two rows on a key.
fn named_for_rows(
    name: &str,
    key: &FieldKey,
    rows: (&FieldRow, &FieldRow),
    models: &FieldModels,
    own: FieldClass,
) -> Option<&'static str> {
    let question = FieldQuestion {
        name,
        key,
        vep: &rows.0[name],
        vep_rs: &rows.1[name],
        strand: "1",
        flags: &[],
        perl_row: Some(rows.0),
        rust_row: Some(rows.1),
        minimal_split_record: false,
        genomic_span: None,
    };
    named(&question, models, own)
}

#[test]
fn two_window_positions_are_the_comparators() {
    for strand in [1, -1] {
        for width in 1..=6 {
            let tx = object(&format!(
                r#"{{"stable_id": "ENST1", "strand": {strand}, "exons": [{{"start": 1, "end": 100}},
                    {{"start": 201, "end": {}}}, {{"start": {}, "end": 400}}]}}"#,
                200 + width,
                301 + width
            ));
            // The exon's bases and the bases before it that both windows hold.
            let held: BTreeSet<Position> = (198 + width..=203).collect();
            assert_eq!(
                transcript_two_window_positions(&tx),
                held,
                "{strand} {width}"
            );
            assert_eq!(!held.is_empty(), width <= 5);
        }
    }
    let short_intron = object(
        r#"{"strand": 1, "exons": [{"start": 1, "end": 100}, {"start": 103, "end": 200}, {"start": 301, "end": 400}]}"#,
    );
    assert_eq!(transcript_two_window_positions(&short_intron), set([]));
    let no_integer = object(r#"{"strand": 1, "exons": [{"start": "x", "end": 100}]}"#);
    assert_eq!(transcript_two_window_positions(&no_integer), set([]));
}

#[test]
fn intron_end_forms_and_gaps_are_the_comparators() {
    assert_eq!(INTRON_END_COPY_BODY, r"c\.([1-9]\d*)-1dup");
    assert_eq!(
        INTRON_END_DUPLICATION_BODY,
        r"c\.([1-9]\d*)([+-])([1-9]\d*)_([1-9]\d*)-1dup"
    );
    assert_eq!(
        INTRON_END_INSERTION_BODY,
        r"c\.([1-9]\d*)-1_([1-9]\d*)ins([A-Z]+)"
    );
    assert_eq!(
        INTRON_END_PAIRS.map(|pair| pair.map(char::from)),
        [['A', 'G'], ['A', 'C']]
    );
    assert_eq!(SHORT_INTRON_BASES, 13);
    for (number, transcript) in transcripts(INTRON_CLASS_TRANSCRIPTS) {
        let base = number * 10_000;
        let gaps = match number {
            1..=14 => set([base + 1199, base + 1499]),
            21 => set([base + 1099, base + 1299]),
            22 => set([base + 1399]),
            _ => set([base + 1313]),
        };
        assert_eq!(transcript_intron_end_gaps(&transcript), gaps, "{number}");
    }
    let bodies = [
        "c.51-1dup",
        "c.51-3_51-1dup",
        "c.50+7_51-1dup",
        "c.51-1_51insT",
        "c.51-1_51insTAT",
        "c.51-1_51insTAG",
        "c.51-1_51insTAC",
        "c.51-1_51-1dup",
        "c.49+7_51-1dup",
        "c.-50-1dup",
        "c.*51-1dup",
        "c.100dup",
        "n.51-1dup",
    ];
    let forms = [
        Some("copy"),
        Some("duplication"),
        Some("duplication"),
        Some("insertion"),
        Some("insertion"),
        Some("pair"),
        Some("pair"),
    ];
    for (at, body) in bodies.iter().enumerate() {
        let form = forms.get(at).copied().flatten();
        assert_eq!(intron_end_insertion_form(body), form, "{body}");
        assert_eq!(
            copy_of_the_intron_last_base(body),
            form == Some("copy"),
            "{body}"
        );
        assert_eq!(
            duplication_of_the_intron_last_bases(body),
            form == Some("duplication"),
            "{body}"
        );
        assert_eq!(
            insertion_that_cannot_end_the_intron(body),
            form == Some("insertion"),
            "{body}"
        );
        assert_eq!(
            insertion_that_ends_in_an_intron_end_pair(body),
            form == Some("pair"),
            "{body}"
        );
    }
}

#[test]
fn intron_start_form_and_gaps_are_the_comparators() {
    assert_eq!(
        INTRON_START_INSERTION_BODY,
        r"c\.([1-9]\d*)_([1-9]\d*)\+1ins([A-Z]+)"
    );
    for (number, transcript) in transcripts(INTRON_CLASS_TRANSCRIPTS) {
        let base = number * 10_000;
        let three_exons = vec![
            (base + 1000, base + 1099),
            (base + 1200, base + 1299),
            (base + 1500, base + 1599),
        ];
        let (forward, exons, low, high, gaps) = match number {
            1..=14 => (
                true,
                three_exons,
                base + 1050,
                base + 1549,
                set([base + 1099, base + 1299]),
            ),
            21 => (
                false,
                three_exons,
                base + 1050,
                base + 1549,
                set([base + 1199, base + 1499]),
            ),
            22 => (
                true,
                vec![
                    (base + 1000, base + 1099),
                    (base + 1200, base + 1299),
                    (base + 1400, base + 1499),
                    (base + 1600, base + 1699),
                ],
                base + 1250,
                base + 1449,
                set([base + 1299]),
            ),
            _ => (
                true,
                vec![
                    (base + 1000, base + 1099),
                    (base + 1113, base + 1199),
                    (base + 1202, base + 1299),
                    (base + 1314, base + 1399),
                ],
                base + 1000,
                base + 1399,
                set([base + 1099, base + 1199, base + 1299]),
            ),
        };
        let coding = CodingExons {
            forward,
            exons,
            low,
            high,
        };
        assert_eq!(
            coding_transcript_exons(&transcript),
            Some(coding),
            "{number}"
        );
        assert_eq!(transcript_intron_start_gaps(&transcript), gaps, "{number}");
    }
    let bodies = [
        "c.50_50+1insT",
        "c.50_50+1insTAT",
        "c.50_51+1insTAT",
        "c.49_50insTAT",
        "c.50dup",
        "c.-51_-51+1insTAT",
        "c.*50_*50+1insTAT",
        "c.50_50+2insTAT",
        "c.50-1_50insTAT",
        "n.50_50+1insTAT",
    ];
    for (at, body) in bodies.iter().enumerate() {
        assert_eq!(insertion_after_the_exon_last_base(body), at < 2, "{body}");
    }
}

#[test]
fn a_transcript_without_a_coding_region_or_with_no_integer_has_no_gap() {
    let coding = |members: &str| {
        coding_transcript_exons(&object(&format!(
            r#"{{"exons": [{{"start": 100, "end": 199}}, {{"start": 300, "end": 399}}]{members}}}"#
        )))
    };
    let read = |forward: bool, low: Position, high: Position| {
        Some(CodingExons {
            forward,
            exons: vec![(100, 199), (300, 399)],
            low,
            high,
        })
    };
    // The strand is forward where the transcript names none, and a coordinate is read as
    // `int()` reads a number, a text of digits and a truth value.
    let coded = r#", "cdna_coding_start": "11", "cdna_coding_end": 150.9"#;
    assert_eq!(coding(coded), read(true, 110, 349));
    assert_eq!(
        coding(&format!(r#"{coded}, "strand": -1"#)),
        read(false, 150, 389)
    );
    assert_eq!(
        coding(&format!(r#"{coded}, "strand": true"#)),
        read(true, 110, 349)
    );
    assert_eq!(
        coding(&format!(r#"{coded}, "strand": " -1 ""#)),
        read(false, 150, 389)
    );
    for members in [
        "",
        r#", "cdna_coding_start": 11"#,
        r#", "cdna_coding_start": 11, "cdna_coding_end": 201"#,
        r#", "cdna_coding_start": 0, "cdna_coding_end": 150"#,
        r#", "cdna_coding_start": null, "cdna_coding_end": 150"#,
        r#", "cdna_coding_start": "x", "cdna_coding_end": 150"#,
        r#", "cdna_coding_start": [11], "cdna_coding_end": 150"#,
        r#", "cdna_coding_start": 11, "cdna_coding_end": 150, "strand": null"#,
        r#", "cdna_coding_start": 11, "cdna_coding_end": 150, "strand": "+""#,
        r#", "cdna_coding_start": 11, "cdna_coding_end": NaN"#,
        r#", "cdna_coding_start": 11, "cdna_coding_end": 9223372036854775808"#,
    ] {
        assert_eq!(coding(members), None, "{members}");
    }
    for exons in [
        r#"{"start": 100}"#,
        r#"[{"start": 100, "end": "x"}]"#,
        r#"[[100, 199]]"#,
        r#"["100-199"]"#,
        r#""100-199""#,
        "7",
        "true",
        "[null]",
    ] {
        let tx = object(&format!(
            r#"{{"exons": {exons}, "cdna_coding_start": 11, "cdna_coding_end": 50}}"#
        ));
        assert_eq!(coding_transcript_exons(&tx), None, "{exons}");
        assert_eq!(transcript_intron_end_gaps(&tx), set([]), "{exons}");
        assert_eq!(transcript_intron_start_gaps(&tx), set([]), "{exons}");
        assert_eq!(transcript_two_window_positions(&tx), set([]), "{exons}");
    }
    for exons in ["[]", "null", "0", r#""""#, "{}", "false"] {
        let tx = object(&format!(
            r#"{{"exons": {exons}, "cdna_coding_start": 11, "cdna_coding_end": 50}}"#
        ));
        assert_eq!(coding_transcript_exons(&tx), None, "{exons}");
        assert_eq!(transcript_two_window_positions(&tx), set([]), "{exons}");
    }
}

#[test]
fn held_insertion_gaps_are_the_comparators() {
    let held = [
        ((100, 1, "GTAT", "TAT", "1"), vec![100]),
        ((100, 1, "TATG", "TAT", "1"), vec![99]),
        ((100, 1, "GATA", "TAT", "-1"), vec![100]),
        ((100, 1, "GATA", "TAT", "-"), vec![100]),
        ((100, 1, "GATA", "TAT", "1"), vec![]),
        ((100, 1, "GGGG", "GGG", "1"), vec![99, 100]),
        ((100, 2, "CGTAT", "TAT", "1"), vec![101]),
        ((100, 1, "GTAT", "T", "1"), vec![]),
        ((100, 1, "TAT", "TAT", "1"), vec![]),
        ((100, 1, "GTAT", "", "1"), vec![]),
    ];
    for ((first, length, allele, inserted, strand), gaps) in held {
        assert_eq!(
            held_insertion_gaps(first, length, allele, inserted, strand),
            gaps,
            "{allele} {inserted} {strand}"
        );
    }
    for (bases, turned) in [
        ("TAT", "ATA"),
        ("GTAT", "ATAC"),
        ("GGGG", "CCCC"),
        ("CGTAT", "ATACG"),
    ] {
        assert_eq!(reverse_complement(bases), turned);
    }
    let bodies = ["c.51-1_51insTAT", "c.50_50+1insT", "c.51-1dup", "c.50dup"];
    let inserted = |pattern: &Regex| bodies.map(|body| notation_inserted_bases(pattern, body));
    assert_eq!(inserted(&INTRON_END_INSERTION), ["TAT", "", "", ""]);
    assert_eq!(inserted(&INTRON_START_INSERTION), ["", "T", "", ""]);
}

#[test]
fn lone_insertion_gaps_are_the_comparators() {
    // The first position and the length of a span, the Allele, the bases a notation
    // inserts, the positions it duplicates and the bases the record moved.
    type Asked<'a> = (Position, Position, &'a str, &'a str, Position, Position);
    let lone: [(Asked, Vec<Position>); 15] = [
        ((100, 1, "TAT", "TAT", 0, 0), vec![99, 100]),
        ((100, 1, "ATA", "TAT", 0, 0), vec![99, 100]),
        ((100, 3, "TAT", "TAT", 0, 0), vec![99, 100, 101, 102]),
        ((100, 1, "CCG", "TAT", 0, 0), vec![]),
        ((100, 1, "TAT", "TTA", 0, 2), vec![99, 100]),
        ((100, 1, "TAT", "TTA", 0, -5), vec![99, 100]),
        ((100, 1, "TAT", "TAT", 0, 2), vec![99, 100]),
        ((100, 1, "ATA", "TAT", 0, -2), vec![99, 100]),
        ((100, 1, "ATT", "TTA", 0, 2), vec![]),
        ((100, 1, "CAG", "", 3, 0), vec![99, 100]),
        ((100, 1, "CAG", "", 1, 0), vec![]),
        ((100, 1, "CAG", "", 3, 1), vec![]),
        ((100, 0, "TAT", "TAT", 0, 0), vec![]),
        ((100, 1, "DELETION", "", 8, 0), vec![]),
        ((100, 1, "-", "", 1, 0), vec![]),
    ];
    for ((first, length, allele, inserted, duplicated, moved), gaps) in lone {
        let found: Vec<Position> =
            lone_insertion_gaps(first, length, allele, inserted, duplicated, moved).collect();
        assert_eq!(found, gaps, "{allele} {inserted} {duplicated} {moved}");
    }
    let bodies = [
        "c.51-1dup",
        "c.51-3_51-1dup",
        "c.50+7_51-1dup",
        "c.51+3_51-1dup",
        "c.52-3_51-1dup",
        "c.51-1_51insTAT",
        "c.50_50+1insTAT",
    ];
    assert_eq!(
        bodies.map(notation_duplicated_positions),
        [1, 3, 0, 0, 0, 0, 0]
    );
}

#[test]
fn classify_names_the_intron_end_insertion_class() {
    let (_scratch, cache) = cache("end", INTRON_CLASS_TRANSCRIPTS);
    let gaps = intron_end_gaps(&cache).unwrap();
    let of = |feature: &str| gaps.get(&("1".to_string(), feature.to_string())).cloned();
    assert_eq!(of("ENST00000000001"), Some(set([11199, 11499])));
    assert_eq!(of("ENST00000000021"), Some(set([211099, 211299])));
    assert_eq!(of("ENST00000000022"), Some(set([221399])));
    assert_eq!(of("ENST00000000023"), Some(set([231313])));
    let fields = [
        ("HGVSc", "ENST00000000001.1:c.51-1dup"),
        ("HGVSp", PROTEIN),
        ("HGVS_OFFSET", "-"),
    ];
    let own = key("1:11199-11200", "G", "ENST00000000001", TERMS);
    let models = |intron_ends: Option<Gaps>| FieldModels {
        introns: Models {
            intron_ends,
            ..Models::default()
        },
        ..FieldModels::default()
    };
    let held = models(Some(gaps));
    let named =
        |perl: &FieldRow, rust: &FieldRow, models: &FieldModels, name: &str, key: &FieldKey| {
            named_for_rows(name, key, (perl, rust), models, intron_end_insertion)
        };
    let (perl, rust) = (row(&fields), with(&fields, &[("HGVSp", "-")]));
    let class = Some(HGVSP_INTRON_END_INSERTION);
    assert_eq!(named(&perl, &rust, &held, "HGVSp", &own), class);
    let prefixed = key("chr1:11199-11200", "G", "ENST00000000001", TERMS);
    assert_eq!(named(&perl, &rust, &held, "HGVSp", &prefixed), class);
    assert_eq!(named(&perl, &rust, &models(None), "HGVSp", &own), None);
    assert_eq!(
        named(&perl, &rust, &models(Some(Gaps::new())), "HGVSp", &own),
        None
    );
    let moved = with(&fields, &[("HGVSp", "-"), ("HGVS_OFFSET", "1")]);
    assert_eq!(named(&perl, &moved, &held, "HGVSp", &own), None);
    let reference_moved = with(&fields, &[("HGVS_OFFSET", "1")]);
    assert_eq!(named(&reference_moved, &rust, &held, "HGVSp", &own), None);
    assert_eq!(named(&perl, &rust, &held, "HGVSc", &own), None);
    let without_offset = row(&fields[..2]);
    let unwritten = with(&fields[..2], &[("HGVSp", "-")]);
    assert_eq!(
        named(&without_offset, &unwritten, &held, "HGVSp", &own),
        None
    );
}

#[test]
fn classify_names_the_intron_start_insertion_class() {
    let (_scratch, cache) = cache("start", INTRON_CLASS_TRANSCRIPTS);
    let gaps = intron_start_gaps(&cache).unwrap();
    let of = |feature: &str| gaps.get(&("1".to_string(), feature.to_string())).cloned();
    assert_eq!(of("ENST00000000001"), Some(set([11099, 11299])));
    assert_eq!(of("ENST00000000021"), Some(set([211199, 211499])));
    assert_eq!(of("ENST00000000022"), Some(set([221299])));
    assert_eq!(of("ENST00000000023"), Some(set([231099, 231199, 231299])));
    let fields = [
        ("HGVSc", "ENST00000000001.1:c.50_50+1insTAT"),
        ("HGVSp", PROTEIN),
        ("HGVS_OFFSET", "2"),
    ];
    let at = |location: &str, feature: &str| key(location, "TAT", feature, TERMS);
    let own = at("1:11097-11098", "ENST00000000001");
    let models = |intron_starts: Option<Gaps>| FieldModels {
        introns: Models {
            intron_starts,
            ..Models::default()
        },
        ..FieldModels::default()
    };
    let held = models(Some(gaps));
    let named =
        |perl: &FieldRow, rust: &FieldRow, models: &FieldModels, name: &str, key: &FieldKey| {
            named_for_rows(name, key, (perl, rust), models, intron_start_insertion)
        };
    let (perl, rust) = (row(&fields), with(&fields, &[("HGVSp", "-")]));
    let class = Some(HGVSP_INTRON_START_INSERTION);
    assert_eq!(named(&perl, &rust, &held, "HGVSp", &own), class);
    let versioned = at("chr1:11097-11098", "ENST00000000001.4");
    assert_eq!(named(&perl, &rust, &held, "HGVSp", &versioned), class);
    let at_the_gap = at("1:11099-11100", "ENST00000000001");
    assert_eq!(named(&perl, &rust, &held, "HGVSp", &at_the_gap), class);
    let at_neither = at("1:11095-11096", "ENST00000000001");
    assert_eq!(named(&perl, &rust, &held, "HGVSp", &at_neither), None);
    assert_eq!(named(&perl, &rust, &models(None), "HGVSp", &own), None);
    assert_eq!(
        named(&perl, &rust, &models(Some(Gaps::new())), "HGVSp", &own),
        None
    );
    let no_offset = [("HGVSp", "-"), ("HGVS_OFFSET", "-")];
    assert_eq!(
        named(&perl, &with(&fields, &no_offset), &held, "HGVSp", &own),
        None
    );
    let unmoved = with(&fields, &[("HGVS_OFFSET", "-")]);
    assert_eq!(
        named(&unmoved, &with(&fields, &no_offset), &held, "HGVSp", &own),
        None
    );
    let back = with(&fields, &[("HGVS_OFFSET", "-2")]);
    let back_unwritten = with(&fields, &[("HGVSp", "-"), ("HGVS_OFFSET", "-2")]);
    assert_eq!(named(&back, &back_unwritten, &held, "HGVSp", &own), None);
    assert_eq!(named(&perl, &rust, &held, "HGVSc", &own), None);
    let without_offset = row(&fields[..2]);
    let unwritten = with(&fields[..2], &[("HGVSp", "-")]);
    assert_eq!(
        named(&without_offset, &unwritten, &held, "HGVSp", &own),
        None
    );
}

#[test]
fn an_offset_is_a_non_zero_integer_in_ascii_digits() {
    assert_eq!(non_zero_integer("2"), Some(2));
    assert_eq!(non_zero_integer("-12"), Some(-12));
    for text in ["", "-", "0", "-0", "02", "+2", "2.0", " 2", "--2", "1_0"] {
        assert_eq!(non_zero_integer(text), None, "{text}");
    }
    assert_eq!(
        non_zero_integer(&"9".repeat(38)).map(|n| n.to_string()),
        Some("9".repeat(38))
    );
    assert_eq!(non_zero_integer(&"9".repeat(39)), None);
}

#[test]
fn classify_names_maxentscan_overlapping_windows_where_the_cache_places_the_record_in_two_windows()
{
    let (_scratch, cache) = cache("windows", MAXENTSCAN_TRANSCRIPTS);
    let positions = two_window_positions(&cache).unwrap();
    let both: Gaps = ["ENST1", "ENST2"]
        .into_iter()
        .map(|feature| (("1".to_string(), feature.to_string()), set([201, 202, 203])))
        .collect();
    assert_eq!(positions, both);
    let models = |two_windows: Option<Gaps>| FieldModels {
        introns: Models {
            two_windows,
            ..Models::default()
        },
        ..FieldModels::default()
    };
    let held = models(Some(positions));
    let own = key("1:202", "G", "ENST1", "splice_region_variant");
    let takes = |name: &str, vep: &str, vep_rs: &str, key: &FieldKey, models: &FieldModels| {
        let question = FieldQuestion {
            name,
            key,
            vep,
            vep_rs,
            strand: "1",
            flags: &[],
            perl_row: None,
            rust_row: None,
            minimal_split_record: false,
            genomic_span: None,
        };
        named(&question, models, maxentscan_overlapping_windows).is_some()
    };
    assert!(takes("MaxEntScan_ref", "8.100", "6.400", &own, &held));
    let versioned = key("chr1:202", "G", "ENST1.3", "splice_region_variant");
    assert!(takes("MaxEntScan_ref", "8.100", "6.400", &versioned, &held));
    assert!(!takes(
        "MaxEntScan_ref",
        "8.100",
        "6.400",
        &own,
        &models(None)
    ));
    assert!(!takes(
        "MaxEntScan_ref",
        "8.100",
        "6.400",
        &own,
        &models(Some(Gaps::new()))
    ));
    assert!(!takes("MaxEntScan_ref", "8.100", "-", &own, &held));
    assert!(!takes("MaxEntScan_ref", "8.100", "abc", &own, &held));
    assert!(!takes("MaxEntScan_ref", "8.100", "Inf", &own, &held));
    let span = key("1:201-203", "G", "ENST1", "splice_region_variant");
    assert!(!takes("MaxEntScan_ref", "8.100", "6.400", &span, &held));
    let elsewhere = key("2:202", "G", "ENST1", "splice_region_variant");
    assert!(!takes(
        "MaxEntScan_ref",
        "8.100",
        "6.400",
        &elsewhere,
        &held
    ));
    assert!(!takes("MES-SWA_donor_ref", "8.100", "6.400", &own, &held));
    // A donor window's nine bases and an acceptor window's twenty-three.
    let (donor, acceptor) = ("AC".repeat(4) + "G", "AC".repeat(11) + "G");
    assert!(takes("MaxEntScan_ref_seq", &donor, &acceptor, &own, &held));
    assert!(!takes("MaxEntScan_ref_seq", &acceptor, &donor, &own, &held));
    let turned: String = donor.chars().rev().collect();
    assert!(!takes("MaxEntScan_alt_seq", &donor, &turned, &own, &held));
    assert!(!takes("MaxEntScan_ref_seq", &donor, "-", &own, &held));
}

#[test]
fn a_feature_is_looked_up_without_its_version() {
    for (feature, name) in [
        ("ENST1.3", "ENST1"),
        ("ENST1", "ENST1"),
        ("ENST1.", "ENST1."),
        ("NM_1.2.34", "NM_1.2"),
        (".5", ""),
        ("12", "12"),
        ("ENST1.3a", "ENST1.3a"),
    ] {
        assert_eq!(without_version(feature), name);
    }
}

#[test]
fn a_models_load_reads_the_two_windows_for_a_maxentscan_field_alone() {
    let (_scratch, cache) = cache("load", MAXENTSCAN_TRANSCRIPTS);
    let manifest = Object::new();
    let load = |fields: &[&str]| {
        let fields: Vec<String> = fields.iter().map(|name| name.to_string()).collect();
        Models::load(&FieldRun {
            corpus: "c",
            cache_dir: &cache,
            manifest: &manifest,
            fields: &fields,
            spliceai_tables: &[],
            reference_fasta: None,
        })
        .unwrap()
    };
    let plain = load(&["HGVSp"]);
    assert_eq!(plain.intron_ends, Some(Gaps::new()));
    assert_eq!(plain.intron_starts, Some(Gaps::new()));
    assert_eq!(plain.two_windows, None);
    let scored = load(&["HGVSp", "MaxEntScan_diff"]);
    assert_eq!(scored.two_windows.map(|held| held.len()), Some(2));
}
