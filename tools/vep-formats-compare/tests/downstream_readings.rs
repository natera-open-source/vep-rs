// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The readers the transcript-reading classes are built of, each held to the values it
//! must return for its arguments: the transcript model, the
//! translation, the Downstream plugin's value, the pairs an entry spells and the texts a
//! number is read from.

use std::panic;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Uncaught;
use vep_formats_compare::genetic_codes::CodonTable;
use vep_formats_compare::genetic_codes::Reading;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::cds_pair_readings;
use vep_formats_compare::rules::downstream::codon_residues;
use vep_formats_compare::rules::downstream::codon_span;
use vep_formats_compare::rules::downstream::codon_table_id;
use vep_formats_compare::rules::downstream::field_absent;
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::downstream::hgvs_rotated_allele;
use vep_formats_compare::rules::downstream::organelle_peptide_path;
use vep_formats_compare::rules::downstream::peptide_description_shape;
use vep_formats_compare::rules::downstream::positive_shift;
use vep_formats_compare::rules::downstream::protein_span;
use vep_formats_compare::rules::downstream::record_span_readings;
use vep_formats_compare::rules::downstream::shift_field_consistent;
use vep_formats_compare::rules::downstream::shifted_pair_readings;
use vep_formats_compare::rules::downstream::translate_until_stop;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::downstream::DownstreamTranscript;
use vep_formats_compare::rules::transcript_model::cache_integer;
use vep_formats_compare::rules::transcript_model::TranscriptModel;
use vep_formats_compare::shared::ABSENT;

const TOO_MANY_DIGITS: &str = "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 digits; use sys.set_int_max_str_digits() to increase the limit";

/// What a run says that ends on an input no reading of the comparison expects.
fn uncaught<T>(run: impl FnOnce() -> T) -> String {
    let Err(payload) = panic::catch_unwind(panic::AssertUnwindSafe(run)) else {
        panic!("the run did not end");
    };
    payload
        .downcast::<Uncaught>()
        .expect("the run ends on an input no reading of the comparison expects")
        .0
}

/// A row of a text format at `location` with the fields `held`.
fn row(location: &str, allele: &str, held: &[(&str, &str)]) -> Entry {
    Entry {
        key: format!("{location}|{allele}|T|Transcript"),
        location: Some(location.to_string()),
        allele: allele.to_string(),
        feature: "T".to_string(),
        consequence: "frameshift_variant".to_string(),
        fields: held
            .iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect(),
        ..Entry::default()
    }
}

/// An entry of the VCF or JSON format: a record's ordinal and no Location.
fn record_entry(kind: &str, held: &[(&str, &str)]) -> Entry {
    Entry {
        key: format!("0|A|T|{kind}"),
        record: Some(0),
        location: None,
        ..row("", "A", held)
    }
}

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

/// Two exons on the forward strand, coding bases 1 to 9 at 3001-3009 and the rest from
/// 3014; and one exon on the reverse strand, coding base c at 2043-c.
fn two_exons() -> TranscriptModel {
    TranscriptModel {
        strand: 1,
        exons: vec![(3001, 3009), (3014, 3046)],
        coding_start: 1,
        coding_end: 27,
        padding: 0,
        pairs: vec![(1, 9, 3001, 3009, 1), (10, 42, 3014, 3046, 1)],
    }
}

fn reverse_strand() -> TranscriptModel {
    TranscriptModel {
        strand: -1,
        exons: vec![(2001, 2042)],
        coding_start: 1,
        coding_end: 27,
        padding: 0,
        pairs: vec![(1, 42, 2001, 2042, -1)],
    }
}

#[test]
fn a_model_maps_a_span_to_its_cds_pair_and_back() {
    let model = two_exons();
    assert_eq!(model.cdna(3009), Some(9));
    assert_eq!(model.cdna(3010), None);
    assert_eq!(model.genomic(10), Some(3014));
    assert_eq!(model.genomic(28), None);
    assert_eq!(model.cds_pair(3009, 3014), (Some(9), Some(10)));
    // An insertion after the first exon's last base and one before the second exon's first
    // print one pair; one inside the intron has no coordinate.
    assert_eq!(model.cds_pair(3010, 3009), (Some(10), Some(9)));
    assert_eq!(model.cds_pair(3014, 3013), (Some(10), Some(9)));
    assert_eq!(model.cds_pair(3012, 3011), (None, None));
    assert_eq!(model.genomic_span(10, 9), None);
    assert_eq!(model.genomic_span(9, 10), Some((3009, 3014)));
    assert_eq!(model.genomic_span(5, 4), Some((3005, 3004)));
    let reverse = reverse_strand();
    assert_eq!(reverse.cdna(2042), Some(1));
    assert_eq!(reverse.genomic(1), Some(2042));
    assert_eq!(reverse.cds_pair(2034, 2035), (Some(8), Some(9)));
    assert_eq!(reverse.cds_pair(2035, 2034), (Some(9), Some(8)));
    assert_eq!(reverse.genomic_span(8, 9), Some((2034, 2035)));
    assert_eq!(reverse.genomic_span(9, 8), Some((2035, 2034)));
}

#[test]
fn a_record_escapes_decay_in_the_last_exon_or_the_51_bases_before_it() {
    let model = two_exons();
    assert!(model.ends_where_nmd_escapes(3014));
    assert!(model.ends_where_nmd_escapes(2958));
    assert!(!model.ends_where_nmd_escapes(2957));
    assert!(reverse_strand().ends_where_nmd_escapes(2001));
    let no_exons = TranscriptModel {
        exons: Vec::new(),
        ..two_exons()
    };
    assert_eq!(
        uncaught(|| no_exons.ends_where_nmd_escapes(5)),
        "tuple index out of range"
    );
}

#[test]
fn a_model_is_read_from_numbers_or_their_texts_and_from_nothing_else() {
    let transcript = |strand: &str, coding_start: &str| {
        vep_formats_compare::pyjson::loads(&format!(
            r#"{{"strand": {strand}, "cdna_coding_start": {coding_start}, "cdna_coding_end": "27",
                "exons": [{{"start": 3014, "end": 3046}}, {{"start": "3001", "end": 3009}}],
                "variation_effect_feature_cache": {{"mapper": {{"pairs": [
                    {{"from_start": 1, "from_end": 9, "to_start": 3001, "to_end": 3009, "ori": 1}},
                    {{"from_start": "10", "from_end": 42, "to_start": 3014, "to_end": 3046, "ori": 1}}]}}}}}}"#
        ))
        .expect("the transcript is JSON")
    };
    assert_eq!(
        TranscriptModel::from_cache(&transcript("1", "1"), 0),
        Some(two_exons())
    );
    let reversed = TranscriptModel::from_cache(&transcript("\"-1\"", "\"1\""), 2).expect("a model");
    assert_eq!(reversed.exons, [(3014, 3046), (3001, 3009)]);
    assert_eq!((reversed.strand, reversed.padding), (-1, 2));
    for (strand, coding_start) in [
        ("2", "1"),
        ("true", "1"),
        ("1.0", "1"),
        ("1", "null"),
        ("1", "\"x\""),
    ] {
        assert_eq!(
            TranscriptModel::from_cache(&transcript(strand, coding_start), 0),
            None,
            "{strand} {coding_start}"
        );
    }
    assert_eq!(cache_integer(&Value::Int(-7)), Some(-7));
    assert_eq!(cache_integer(&text("-07")), Some(-7));
    assert_eq!(cache_integer(&text("٥")), Some(5));
    assert_eq!(cache_integer(&text("+5")), None);
    assert_eq!(cache_integer(&Value::Bool(true)), None);
    assert_eq!(cache_integer(&Value::Float(5.0)), None);
    assert_eq!(
        uncaught(|| cache_integer(&text(&"1".repeat(4301)))),
        TOO_MANY_DIGITS
    );
}

/// A whole number is held in 64 bits; a number they do not hold is read as no number.
#[test]
fn a_number_past_64_bits_is_read_as_no_number() {
    assert_eq!(cache_integer(&text("9223372036854775807")), Some(i64::MAX));
    assert_eq!(cache_integer(&text("9223372036854775808")), None);
    assert_eq!(
        cache_integer(&Value::BigInt("9223372036854775808".to_string())),
        None
    );
    let long = "12345678901234567890";
    assert_eq!(
        cds_pair_readings(&row("21:5", "-", &[("CDS_position", long)])),
        []
    );
    assert_eq!(
        cds_pair_readings(&row("21:5", "-", &[("CDS_position", &format!("5-{long}"))])),
        [(Some(5), None)]
    );
    assert_eq!(
        record_span_readings(&row(&format!("21:{long}"), "-", &[])),
        []
    );
    assert_eq!(
        protein_span(&row("21:5", "-", &[("Protein_position", long)])),
        None
    );
}

#[test]
fn a_translation_follows_the_table_and_ends_at_the_first_stop() {
    for reading in [Reading::EnsemblVep, Reading::VepRs] {
        let translated =
            |sequence: &str, table: i64| translate_until_stop(sequence, table, reading);
        assert_eq!(
            translated("ATGGAYAARGARCARRAYSARNNNatgATGTARCC", 1),
            "MDKEQBZXXM"
        );
        assert_eq!(translated("atg", 1), "X");
        assert_eq!(translated("ATGATAAGATGA", 2), "MM");
        assert_eq!(translated("AT", 1), "");
        assert_eq!(uncaught(|| translated("ATG", 7)), "7");
    }
    let tables =
        |values: &[Value]| -> Vec<CodonTable> { values.iter().map(codon_table_id).collect() };
    let lacking = |spelt: &str| CodonTable::Lacking(spelt.to_string());
    assert_eq!(
        tables(&[
            text("٢"),
            text("02"),
            text(" 2"),
            text("-2"),
            Value::Float(2.0)
        ]),
        [
            CodonTable::Held(2),
            CodonTable::Held(2),
            lacking(" 2"),
            lacking("-2"),
            lacking("2.0")
        ]
    );
    assert_eq!(
        uncaught(|| codon_table_id(&text(&"1".repeat(4301)))),
        TOO_MANY_DIGITS
    );
}

#[test]
fn the_downstream_value_is_each_engines_rule_on_one_pair() {
    let transcript = DownstreamTranscript {
        coding: "ATGAAAACTCTTTAA".to_string(),
        utr: "GGGTAGCC".to_string(),
        table: CodonTable::Held(1),
        peptide_length: None,
        peptide: None,
        model: None,
    };
    let value = |start: i64, end: i64, allele: &str, length_change: bool, ensembl: bool| {
        transcript.downstream_value(start, end, allele, length_change, ensembl)
    };
    let protein = |residues: &str| Some(residues.to_string());
    // A start of 0 is one base before the first, which a slice counts from the sequence's end.
    assert_eq!(value(0, 0, "", false, false), protein(""));
    assert_eq!(value(0, 0, "", false, true), protein("MKTL"));
    assert_eq!(value(6, 6, "", false, true), protein("LFKGS"));
    assert_eq!(value(6, 6, "", false, false), protein("KLFKGS"));
    assert_eq!(value(5, 4, "GC", true, true), protein("4"));
    assert_eq!(value(20, 20, "", false, true), None);
    assert_eq!(value(20, 20, "", false, false), protein("G"));
    assert_eq!(value(16, 15, "A", false, true), protein("RVA"));
    assert_eq!(value(17, 16, "A", false, true), None);
    let no_table = DownstreamTranscript {
        coding: "ATG".to_string(),
        utr: String::new(),
        table: CodonTable::Lacking("7".to_string()),
        peptide_length: None,
        peptide: None,
        model: None,
    };
    assert_eq!(no_table.downstream_value(9, 9, "", false, true), None);
    assert_eq!(
        uncaught(|| no_table.downstream_value(1, 1, "", false, false)),
        "7"
    );
    for reading in [Reading::EnsemblVep, Reading::VepRs] {
        assert_eq!(
            [(1, 2), (5, 5), (5, 6), (0, 1), (2, 1)].map(|codons| codon_residues(
                &transcript,
                codons,
                reading
            )),
            [protein("MK"), protein("*"), None, None, protein("")]
        );
        assert_eq!(uncaught(|| codon_residues(&no_table, (1, 1), reading)), "7");
    }
    assert_eq!(
        [(5, 4), (4, 5), (7, 7), (0, 1), (-1, -4)].map(codon_span),
        [(2, 2), (2, 2), (3, 3), (0, 1), (-1, 0)]
    );
}

#[test]
fn the_allele_is_turned_by_the_shift_and_read_on_the_transcripts_strand() {
    for (allele, shift, strand, read) in [
        ("ACG", 1, 1, "CGA"),
        ("ACG", 4, 1, "CGA"),
        ("ACG", 1, -1, "GTC"),
        ("ACG", 5, -1, "CGT"),
        ("ACG", -2, -1, "TCG"),
        ("", 3, 1, ""),
        ("AC", 0, -1, "GT"),
    ] {
        assert_eq!(
            hgvs_rotated_allele(allele, shift, strand),
            read,
            "{allele} {shift} {strand}"
        );
    }
}

#[test]
fn an_entry_spells_its_pairs_and_spans_in_the_forms_the_formats_write() {
    let pairs = |held: &[(&str, &str)], allele: &str| cds_pair_readings(&row("21:5", allele, held));
    assert_eq!(
        pairs(&[("CDS_position", "5-6/90")], "A"),
        [(Some(5), Some(6)), (Some(6), Some(5))]
    );
    assert_eq!(
        pairs(&[("CDS_position", "5-6/90")], "-"),
        [(Some(5), Some(6))]
    );
    assert_eq!(pairs(&[("CDS_position", "?-6")], "A"), [(None, Some(6))]);
    assert_eq!(
        pairs(&[("cds_start", "٥"), ("cds_end", "٦")], "-"),
        [(Some(5), Some(6))]
    );
    assert_eq!(pairs(&[("CDS_position", "-")], "A"), []);
    assert_eq!(pairs(&[], "A"), []);
    assert_eq!(
        uncaught(|| pairs(&[("CDS_position", "²")], "A")),
        "invalid literal for int() with base 10: '²'"
    );
    assert_eq!(
        record_span_readings(&row("chr:1:7-8", "A", &[])),
        [(7, 8), (8, 7)]
    );
    assert_eq!(record_span_readings(&row("21:7", "A", &[])), [(7, 7)]);
    assert_eq!(record_span_readings(&record_entry("Transcript", &[])), []);
    assert_eq!(
        uncaught(|| record_span_readings(&row("21:²", "A", &[]))),
        "invalid literal for int() with base 10: '²'"
    );
    assert_eq!(
        protein_span(&row("21:5", "A", &[("Protein_position", "9-4/100")])),
        Some((4, 9))
    );
    assert_eq!(
        protein_span(&row("21:5", "A", &[("protein_start", "3")])),
        None
    );
}

#[test]
fn a_shifted_reading_is_the_own_pair_the_moved_pair_and_the_records_end() {
    let e = row("21:2034-2035", "GA", &[("CDS_position", "8-9")]);
    assert_eq!(
        shifted_pair_readings(&e, &reverse_strand(), -1, &cds_pair_readings(&e)),
        [
            ((Some(8), Some(9)), Some((9, 10)), 2035),
            ((Some(9), Some(8)), Some((10, 9)), 2034),
        ]
    );
    // An entry without a Location reads the one span each own pair maps back to: the
    // insertion's pair across the intron maps back to none.
    assert_eq!(
        shifted_pair_readings(
            &record_entry("Transcript", &[]),
            &two_exons(),
            5,
            &[(Some(10), Some(9)), (Some(9), Some(10)), (None, Some(5))]
        ),
        [((Some(9), Some(10)), Some((10, 15)), 3014)]
    );
}

#[test]
fn a_shift_length_is_positive_as_a_whole_number_of_any_length() {
    let shift = |value: &str| positive_shift(&row("21:5", "A", &[("SHIFT_LENGTH", value)]));
    assert_eq!(shift(&"9".repeat(640)), Ok(true));
    assert_eq!(
        ["0", "000", "3", "-3", "", "٣"].map(|value| shift(value).expect("a shift length")),
        [false, false, true, false, false, true]
    );
    assert_eq!(
        shift("²").expect_err("no integer").to_string(),
        "SHIFT_LENGTH '²' is not an integer"
    );
    let consistent = |held: &[(&str, &str)]| {
        shift_field_consistent(&row("21:5", "A", held)).expect("a shift length")
    };
    assert!(consistent(&[]));
    assert!(!consistent(&[("shift_length", "0")]));
    assert!(consistent(&[("Shift_Length", "0")]));
    assert!(consistent(&[("SHIFT_LENGTH", "2")]));
}

#[test]
fn a_field_is_found_by_its_name_or_prefix_and_an_absent_value_by_its_format() {
    let e = row(
        "21:5",
        "A",
        &[("STRAND", "-1"), ("strand", "1"), ("gnomADe_AF", "0.1")],
    );
    assert_eq!(field_of(Some(&e), &["strand"]), "-1");
    assert_eq!(field_of(Some(&e), &["gnomade_"]), "0.1");
    assert_eq!(field_of(Some(&e), &["gnomad_"]), ABSENT);
    assert_eq!(field_of(Some(&e), &[]), ABSENT);
    assert_eq!(field_of(None, &["STRAND"]), ABSENT);
    assert_eq!(
        ["-", ABSENT, "", "x"].map(|value| field_absent(&e, value)),
        [true, true, false, false]
    );
    // The empty cell of a CSQ entry is no value; the empty string of a JSON entry is one.
    assert!(field_absent(&record_entry("Transcript", &[]), ""));
    assert!(!field_absent(
        &record_entry("transcript", &[("variant_allele", "A")]),
        ""
    ));
    assert_eq!(
        [
            "Met3dup",
            "Met3_Gln4insMet",
            "Ala2delinsGly",
            "delinsins",
            "",
            "insdelins"
        ]
        .map(peptide_description_shape),
        ["dup", "ins", "other", "ins", "other", "ins"]
    );
}

/// An inserted allele of two bases shifted past its own length by a length that is no
/// multiple of it stays a difference on a row that states no forward strand; the shift is
/// read as a whole number of any length the interpreter converts.
#[test]
fn a_shift_of_any_length_is_weighed_against_the_inserted_alleles_length() {
    let taken = |offset: &str, strand: &str| {
        let e = row(
            "MT:104",
            "CA",
            &[("HGVS_OFFSET", offset), ("STRAND", strand)],
        );
        organelle_peptide_path(&e, &e, "P:p.Ala2GlyfsTer?", "P:p.Ala2GlyfsTer3")
    };
    let even = format!("1{}", "0".repeat(30));
    let odd = format!("1{}1", "0".repeat(29));
    assert!(taken(&even, "-1"));
    assert!(!taken(&format!("-{odd}"), "-1"));
    assert!(taken(&odd, "1"));
    assert!(taken("2", "-1") && taken("4", "-1") && !taken("3", "-1"));
    assert_eq!(uncaught(|| taken(&"9".repeat(4301), "-1")), TOO_MANY_DIGITS);
}

#[test]
fn a_cache_without_a_transcripts_directory_is_refused_and_its_readings_are_kept_by_stable_id() {
    let root = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("downstream-readings");
    let _ = std::fs::remove_dir_all(&root);
    let refused = DownstreamSequences::new(&root.join("nowhere")).expect_err("no directory");
    assert_eq!(
        refused.to_string(),
        format!(
            "--json-cache: {} has no transcripts/ directory",
            root.join("nowhere").display()
        )
    );
    let shard = |chromosome: &str, name: &str, text: &str| {
        let path = root.join("transcripts").join(chromosome).join(name);
        std::fs::create_dir_all(path.parent().expect("a directory")).expect("the directory");
        std::fs::write(path, text).expect("the shard");
    };
    let transcript = |stable_id: &str, start: i64, coding: &str| {
        format!(
            r#"{{"stable_id": "{stable_id}", "strand": 1, "cdna_coding_start": 1, "cdna_coding_end": 6,
                "exons": [{{"start": {start}, "end": {end}}}],
                "variation_effect_feature_cache": {{"translateable_seq": "{coding}", "three_prime_utr": null,
                    "codon_table": 0, "peptide": "MK",
                    "mapper": {{"pairs": [{{"from_start": 1, "from_end": 9, "to_start": {start}, "to_end": {end}, "ori": 1}}]}}}}}}"#,
            end = start + 8
        )
    };
    shard(
        "X",
        "1-5.json",
        &format!(
            r#"{{"1-5": [null, {}, {}, {}]}}"#,
            transcript("SAME", 101, "nnGaaa"),
            transcript("TWO_MODELS", 101, "ATGAAA"),
            transcript("TWO_CODINGS", 101, "ATGAAA")
        ),
    );
    shard(
        "Y",
        "1-5.json",
        &format!(
            "[{}, {}, {}]",
            transcript("SAME", 101, "NNGAAA"),
            transcript("TWO_MODELS", 501, "ATGAAA"),
            transcript("TWO_CODINGS", 501, "ATGAAG")
        ),
    );
    let sequences = DownstreamSequences::new(&root).expect("the cache has transcripts");
    let read = |location: Option<&str>, feature: &str| {
        let e = Entry {
            feature: feature.to_string(),
            location: location.map(str::to_string),
            ..row("", "A", &[])
        };
        sequences.transcript(&e).expect("the cache reads")
    };
    let same = read(Some("X:5"), "SAME").expect("the transcript is held");
    assert_eq!((same.coding.as_str(), same.utr.as_str()), ("NNGAAA", ""));
    assert_eq!(
        (&same.table, same.peptide_length),
        (&CodonTable::Held(1), Some(2))
    );
    assert_eq!(same.model.as_ref().map(|model| model.padding), Some(2));
    // A Feature carrying a version the stable id lacks is looked up without it.
    assert!(read(Some("X:5"), "SAME.7").is_some());
    assert!(read(Some("X:5"), "SAME.x").is_none());
    // An entry that names a chromosome reads its shards alone; one that names none, or
    // one the cache has no directory for, reads every shard: one reading under two models
    // is kept without a model, two readings are left out.
    assert!(read(Some("Y:5"), "TWO_MODELS").is_some_and(|transcript| transcript.model.is_some()));
    for location in [None, Some("Z:5"), Some("")] {
        assert!(read(location, "SAME").is_some_and(|transcript| transcript.model.is_some()));
        assert!(read(location, "TWO_MODELS").is_some_and(|transcript| transcript.model.is_none()));
        assert!(read(location, "TWO_CODINGS").is_none());
    }
    assert!(read(Some("X:5"), "TWO_CODINGS").is_some());
    shard("Z", "1-5.json", "[");
    let unreadable = DownstreamSequences::new(&root).expect("the cache has transcripts");
    let e = Entry {
        location: Some("Z:5".to_string()),
        ..row("", "A", &[])
    };
    assert!(unreadable
        .transcript(&e)
        .expect_err("the shard is not JSON")
        .to_string()
        .contains("1-5.json is not JSON: "));
    std::fs::remove_dir_all(&root).expect("the files are removed");
}
