// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The UTRAnnotator classes: the rules on small rows, the reader of the seeded runs on
//! small files, and the `score` command over small outputs under the classes' own rows
//! of the class table. The rules' expected values are those `tools/vep-formats-compare`
//! gives on the same inputs.

use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Mutex;
use vep_plugin_pairs::registry::BuiltInRegistry;
use vep_plugin_pairs::rules::fields::class_fields;
use vep_plugin_pairs::rules::utrannotator::{
    annotation_frames, annotation_groups, reordered_pair_groups, utr_base_after_insertion,
    UtrannotatorSeedSweep,
};
use vep_plugin_pairs::rules::{
    blank, parse_default, BuiltIn, Consult, RecordRows, Refusal, Row, Rules,
};
use vep_plugin_pairs::score::run;

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
/// The tree's class table.
const COMMITTED_CLASS_TABLE: &str =
    include_str!("../../../scripts/concordance/plugin_class_table.tsv");

const PAIR_ORDER: &str = "utrannotator_annotation_pair_order";
const WINDOW: &str = "utrannotator_gained_start_window";
const PRECEDENCE: &str = "utrannotator_consequence_precedence";
const STRUCTURAL: &str = "utrannotator_structural_allele_text";
const LOST_AFTER: &str = "utrannotator_lost_start_after_insertion";
const REMADE: &str = "utrannotator_gained_start_remade_in_place";
const BOUNDARY: &str = "utrannotator_frameshift_at_stop_boundary";
const BEFORE_STOP: &str = "utrannotator_lost_stop_written_before_stop_codon";

const GAIN: &str = "5_prime_UTR_premature_start_codon_gain_variant";
const LOSS: &str = "5_prime_UTR_premature_start_codon_loss_variant";
const STOP_LOSS: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";
const STOP_GAIN: &str = "5_prime_UTR_uORF_stop_codon_gain_variant";
const FRAMESHIFT: &str = "5_prime_UTR_uORF_frameshift_variant";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-utrannotator-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_string()
    }

    fn file(&self, name: &str, content: impl AsRef<[u8]>) -> String {
        fs::write(self.0.join(name), content).unwrap();
        self.path(name)
    }

    /// A default-format output: the column header and `rows`.
    fn output(&self, name: &str, rows: &[String]) -> String {
        self.file(name, format!("{HEADER}{}", rows.concat()))
    }

    /// An input VCF of `records`, each a CHROM, POS, ID, REF and ALT.
    fn vcf(&self, name: &str, records: &[[&str; 5]]) -> String {
        let lines: String = records
            .iter()
            .map(|record| format!("{}\t.\t.\t.\n", record.join("\t")))
            .collect();
        self.file(name, format!("{VCF_HEADER}{lines}"))
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn gzip(content: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(content).unwrap();
    encoder.finish().unwrap()
}

/// A row on a transcript of the gene `ENSG1` with the Consequence `missense_variant`.
fn row(
    location: &str,
    allele: &str,
    feature: &str,
    extra: &str,
    cdna: &str,
    uploaded: &str,
) -> String {
    format!("{uploaded}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t{cdna}\t-\t-\t-\t-\t-\t{extra}\n")
}

/// A 5' UTR row with its cDNA_position and the two UTRAnnotator fields.
fn utr_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    cdna: &str,
    effect: (&str, &str),
) -> String {
    let (consequence, annotation) = effect;
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t5_prime_UTR_variant\t{cdna}\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;5UTR_consequence={consequence};5UTR_annotation={annotation}\n")
}

fn default_text(rows: &[String]) -> String {
    format!("{HEADER}{}", rows.concat())
}

fn entries(rows: &[String]) -> Vec<Row> {
    parse_default(&default_text(rows)).unwrap()
}

/// The built-in class's answer for the field `name` of the reference's row `e`.
fn accepts(
    klass: &str,
    e: &Row,
    name: &str,
    expected: &str,
    actual: &str,
    a: Option<&Row>,
    sweep: Option<&UtrannotatorSeedSweep>,
) -> bool {
    BuiltIn
        .accepts(&Consult {
            klass,
            e,
            name,
            expected,
            actual,
            a,
            record_rows: None,
            candidate_rows: None,
            declared: &BTreeSet::new(),
            assembly: None,
            source: None,
            gwas_pairs: None,
            revel_table: None,
            utrannotator_sweep: sweep,
            downstream_sequences: None,
            splice_windows: None,
            spliceai_table: None,
            exon_spans: None,
        })
        .unwrap()
}

/// The seed sweep of default-format outputs, one text per seed.
fn seed_sweep(scratch: &Scratch, tag: &str, outputs: &[&str]) -> UtrannotatorSeedSweep {
    let paths: Vec<PathBuf> = outputs
        .iter()
        .enumerate()
        .map(|(number, text)| PathBuf::from(scratch.file(&format!("{tag}-seed{number}.txt"), text)))
        .collect();
    UtrannotatorSeedSweep::from_outputs(&paths).unwrap()
}

fn refusal(read: Result<UtrannotatorSeedSweep, Refusal>) -> String {
    match read {
        Ok(_) => panic!("the outputs were read"),
        Err(Refusal(message)) => message,
    }
}

const ANNOTATION_PAIRS: [(&str, &str); 7] = [
    ("CapDistanceToStart", "407"),
    ("DistanceToCDS", "23"),
    ("DistanceToStop", "24"),
    ("Evidence", "False"),
    ("KozakContext", "ATAATGG"),
    ("KozakStrength", "Strong"),
    ("type", "OutOfFrame_oORF"),
];
const HASH_ORDER: [&str; 7] = [
    "DistanceToStop",
    "CapDistanceToStart",
    "KozakContext",
    "KozakStrength",
    "DistanceToCDS",
    "Evidence",
    "type",
];
const KEY_ORDER: [&str; 7] = [
    "CapDistanceToStart",
    "DistanceToCDS",
    "DistanceToStop",
    "Evidence",
    "KozakContext",
    "KozakStrength",
    "type",
];

/// One `5UTR_annotation` group: the pairs of `keys` in that order.
fn utr_annotation(keys: &[&str]) -> String {
    let pairs: Vec<String> = keys
        .iter()
        .map(|key| {
            let (_, value) = ANNOTATION_PAIRS
                .iter()
                .find(|(name, _)| name == key)
                .unwrap();
            format!("{key}={value}")
        })
        .collect();
    pairs.join(":")
}

#[test]
fn annotation_groups_set_the_order_of_groups_and_pairs_aside() {
    assert_eq!(
        annotation_groups("b=2:a=1,c=3"),
        annotation_groups("c=3,a=1:b=2")
    );
    assert_ne!(
        annotation_groups("b=2:a=1,c=3"),
        annotation_groups("a=1,b=2:c=3")
    );
    assert_ne!(annotation_groups("a=1,a=1"), annotation_groups("a=1"));
}

#[test]
fn pair_order_refuses_a_value_carrying_both_list_separators() {
    assert!(!reordered_pair_groups("x=1,b=2&a=3", "a=3&x=1,b=2"));
}

#[test]
fn pair_order_takes_the_same_pairs_only_where_the_engine_writes_them_in_key_order() {
    assert!(!reordered_pair_groups("b=1:a=2", "b=1:a=2"));
    assert!(!reordered_pair_groups("a=2:b=1", "b=1:a=2"));
    // The same value twice, in key order, is no reordering.
    assert!(!reordered_pair_groups("a=2:b=1", "a=2:b=1"));
    // Key order is the keys' own: `k` before `k2`, where the whole pairs sort the other way.
    assert!(reordered_pair_groups("k2=y:k=x", "k=x:k2=y"));
    assert!(!reordered_pair_groups("k=x:k2=y", "k2=y:k=x"));
    // A text that is not key=value pairs is no group of pairs, whatever its order.
    assert!(!reordered_pair_groups("b:a", "a:b"));
    assert!(reordered_pair_groups("b=1:a=2", "a=2:b=1"));
}

#[test]
fn pair_order_takes_groups_in_another_order_and_refuses_a_pair_replaced_missing_or_duplicated() {
    let (hash_order, key_order) = (utr_annotation(&HASH_ORDER), utr_annotation(&KEY_ORDER));
    assert!(reordered_pair_groups(&hash_order, &key_order));
    // Two groups: the groups in another order too, each with its own pairs; `&` joins the
    // groups of a VCF value.
    assert!(reordered_pair_groups("k=2:t=v,t=u:k=1", "k=1:t=u,k=2:t=v"));
    assert!(reordered_pair_groups("k=2:t=v&t=u:k=1", "k=1:t=u&k=2:t=v"));
    // A pair that stands in another group.
    assert!(!reordered_pair_groups("k=2:t=v,t=u:k=1", "k=1:t=v,k=2:t=u"));
    // A pair replaced, missing on either side, or written once more.
    assert!(!reordered_pair_groups(
        &hash_order,
        &key_order.replace("DistanceToCDS=23", "DistanceToCDS=24")
    ));
    assert!(!reordered_pair_groups(
        &hash_order,
        &utr_annotation(&KEY_ORDER[..6])
    ));
    assert!(!reordered_pair_groups(
        &utr_annotation(&HASH_ORDER[1..]),
        &key_order
    ));
    assert!(!reordered_pair_groups(
        &hash_order,
        &format!("{key_order},{key_order}")
    ));
    assert!(!reordered_pair_groups(
        &format!("{hash_order},{hash_order}"),
        &key_order
    ));
    // The engine's value out of key order, and the two values in each other's place.
    let reversed: Vec<&str> = KEY_ORDER.iter().rev().copied().collect();
    assert!(!reordered_pair_groups(
        &hash_order,
        &utr_annotation(&reversed)
    ));
    assert!(!reordered_pair_groups(&key_order, &hash_order));
    // A value one engine does not write.
    assert!(!reordered_pair_groups(&hash_order, "-"));
    assert!(!reordered_pair_groups("-", &key_order));
}

#[test]
fn pair_order_reads_the_annotation_field_alone() {
    let (hash_order, key_order) = (utr_annotation(&HASH_ORDER), utr_annotation(&KEY_ORDER));
    let rows = entries(&[utr_row(
        "utr_substitution",
        "21:2000984",
        "G",
        "ENST00000000032",
        "410",
        (LOSS, &hash_order),
    )]);
    let entry = &rows[0];
    let asked = |name: &str| accepts(PAIR_ORDER, entry, name, &hash_order, &key_order, None, None);
    assert!(asked("5UTR_annotation"));
    for name in [
        "5UTR_consequence",
        "NearestExonJB",
        "DOMAINS",
        "FREQS",
        "Existing_uORFs",
        "SpliceVault_site_type",
        "5UTR_annotation_x",
        "Extra",
    ] {
        assert!(!asked(name), "{name}");
    }
}

const LOST_START: (&str, &str) = (
    LOSS,
    "CapDistanceToStart=290:KozakStrength=Moderate:type=uORF",
);
const LOST_STOP: (&str, &str) = (
    STOP_LOSS,
    "AltStop=False:KozakContext=GGGATGC:KozakStrength=Moderate",
);

/// A substitution's row under an identifier on which two of the plugin's analyses fire.
fn substitution(effect: (&str, &str)) -> String {
    utr_row("rs1", "21:35748195", "T", "ENST1", "293", effect)
}

#[test]
fn precedence_accepts_two_terms_the_seeded_runs_write_on_the_row() {
    let scratch = Scratch::new("precedence-accepts");
    let agreeing = utr_row(
        "rs2",
        "21:400",
        "C",
        "ENST2",
        "40",
        (
            GAIN,
            "CapDistanceToStart=38:DistanceToCDS=110:KozakStrength=Weak:type=uORF",
        ),
    );
    let reference = default_text(&[substitution(LOST_START), agreeing.clone()]);
    let candidate = default_text(&[substitution(LOST_STOP), agreeing.clone()]);
    // The second run writes its annotation's pairs in another order.
    let other_seed = default_text(&[
        substitution((
            STOP_LOSS,
            "KozakStrength=Moderate:AltStop=False:KozakContext=GGGATGC",
        )),
        agreeing,
    ]);
    let sweep = seed_sweep(&scratch, "both", &[&reference, &other_seed]);
    let (e, a) = (
        &parse_default(&reference).unwrap()[0],
        &parse_default(&candidate).unwrap()[0],
    );
    let asked = |name: &str, sweep: Option<&UtrannotatorSeedSweep>, a: Option<&Row>| {
        let (expected, actual) = match name {
            "5UTR_consequence" => (LOST_START.0, LOST_STOP.0),
            _ => (LOST_START.1, LOST_STOP.1),
        };
        accepts(PRECEDENCE, e, name, expected, actual, a, sweep)
    };
    for name in ["5UTR_consequence", "5UTR_annotation"] {
        assert!(asked(name, Some(&sweep), Some(a)), "{name}");
        // Nothing without the seeded runs, and nothing without the engine's row.
        assert!(!asked(name, None, Some(a)), "{name}");
        assert!(!asked(name, Some(&sweep), None), "{name}");
    }
    assert!(!asked("Existing_uORFs", Some(&sweep), Some(a)));
}

#[test]
fn precedence_refuses_a_term_the_seeded_runs_never_write_on_the_row() {
    let scratch = Scratch::new("precedence-refuses");
    let reference = default_text(&[substitution(LOST_START)]);
    let candidate = default_text(&[substitution(LOST_STOP)]);
    let elsewhere = default_text(&[utr_row(
        "rs1",
        "21:35748196",
        "T",
        "ENST1",
        "294",
        LOST_STOP,
    )]);
    let another_record = default_text(&[utr_row(
        "rs9",
        "21:35748195",
        "T",
        "ENST1",
        "293",
        LOST_STOP,
    )]);
    let e = &parse_default(&reference).unwrap()[0];
    let asked = |name: &str, a: &Row, sweep: &UtrannotatorSeedSweep| {
        let field = |entry: &Row| entry.fields.get(name).unwrap_or("-").to_string();
        accepts(
            PRECEDENCE,
            e,
            name,
            &field(e),
            &field(a),
            Some(a),
            Some(sweep),
        )
    };
    let a = &parse_default(&candidate).unwrap()[0];
    // The runs write the reference's term alone, hold another row, hold the row of another
    // record at the Location, or write the engine's term alone.
    let sweeps = [
        seed_sweep(&scratch, "one-term", &[&reference, &reference]),
        seed_sweep(&scratch, "elsewhere", &[&elsewhere]),
        seed_sweep(&scratch, "another-record", &[&reference, &another_record]),
        seed_sweep(&scratch, "engine-term", &[&candidate]),
    ];
    for (number, sweep) in sweeps.iter().enumerate() {
        for name in ["5UTR_consequence", "5UTR_annotation"] {
            assert!(!asked(name, a, sweep), "{number} {name}");
        }
    }
    let both_terms = seed_sweep(&scratch, "both", &[&reference, &candidate]);
    assert!(asked("5UTR_consequence", a, &both_terms));
    // A `-` on the engine's side, a value outside the five terms and the same term on
    // both sides stay differences beside runs that write both terms.
    for wrong in [
        substitution(("-", "-")),
        substitution(("5_prime_UTR_variant", LOST_STOP.1)),
        substitution(LOST_START),
    ] {
        let wrong = &entries(&[wrong])[0];
        assert!(!asked("5UTR_consequence", wrong, &both_terms));
    }
    // An annotation the runs do not write beside the engine's term stays a difference
    // under an accepted consequence.
    let other_annotation =
        &parse_default(&candidate.replace("AltStop=False", "AltStop=True")).unwrap()[0];
    assert!(asked("5UTR_consequence", other_annotation, &both_terms));
    assert!(!asked("5UTR_annotation", other_annotation, &both_terms));
    // The reference's own annotation beside the engine's term, as written or with its
    // pairs in another order, is one the runs write beside the reference's term alone.
    for shared in [
        LOST_START.1,
        "type=uORF:CapDistanceToStart=290:KozakStrength=Moderate",
    ] {
        let sharing = &entries(&[substitution((STOP_LOSS, shared))])[0];
        assert!(!asked("5UTR_consequence", sharing, &both_terms), "{shared}");
    }
}

#[test]
fn seeded_runs_name_the_one_index_of_a_gained_start() {
    let scratch = Scratch::new("gained-start");
    let insertion = |effect: (&str, &str)| {
        row(
            "1:4000-4001",
            "TG",
            "ENST4",
            &format!(
                "IMPACT=MODIFIER;STRAND=1;5UTR_annotation={};5UTR_consequence={};Existing_uORFs=1",
                effect.1, effect.0
            ),
            "300-301",
            "swept",
        )
    };
    let shifted = insertion((FRAMESHIFT, "alt_type=uORF:ref_type=uORF"));
    let gained = insertion((GAIN, "CapDistanceToStart=299:type=uORF"));
    let elsewhere = insertion((GAIN, "type=uORF:CapDistanceToStart=298"));
    let e = &entries(std::slice::from_ref(&shifted))[0];
    let one = seed_sweep(
        &scratch,
        "one",
        &[
            &default_text(std::slice::from_ref(&shifted)),
            &default_text(std::slice::from_ref(&gained)),
        ],
    );
    assert_eq!(one.gained_start(e), Some(299));
    // Two runs that name two indexes name none, as do runs that write no gained start.
    let two = seed_sweep(
        &scratch,
        "two",
        &[&default_text(&[gained]), &default_text(&[elsewhere])],
    );
    assert_eq!(two.gained_start(e), None);
    let none = seed_sweep(&scratch, "none", &[&default_text(&[shifted])]);
    assert_eq!(none.gained_start(e), None);
}

#[test]
fn annotation_frames_read_the_text_and_the_json_spellings() {
    let frame = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect()
    };
    assert_eq!(
        annotation_frames("a=1:b=2,c=3"),
        Some(vec![frame(&[("a", "1"), ("b", "2")]), frame(&[("c", "3")])])
    );
    assert_eq!(
        annotation_frames("a=1:b=2&c=3"),
        Some(vec![frame(&[("a", "1"), ("b", "2")]), frame(&[("c", "3")])])
    );
    // A later pair of a key replaces the earlier one's value in its place.
    assert_eq!(
        annotation_frames("a=1:b=2:a=3"),
        Some(vec![frame(&[("a", "3"), ("b", "2")])])
    );
    assert_eq!(
        annotation_frames(
            r#"{"1": {"type": "uORF", "CapDistanceToStart": 64, "Evidence": false, "more": [1, "x"], "none": null, "ratio": 0.5}}"#
        ),
        Some(vec![frame(&[
            ("type", "uORF"),
            ("CapDistanceToStart", "64"),
            ("Evidence", "False"),
            ("more", "[1, 'x']"),
            ("none", "None"),
            ("ratio", "0.5"),
        ])])
    );
    for none in [
        "-",
        "",
        "<absent>",
        "a=1,b=2&c=3",
        "a=1:b",
        "{",
        "{\"1\": 2}",
        "[1]",
    ] {
        assert_eq!(annotation_frames(none), None, "{none}");
    }
}

#[test]
fn base_after_an_insertion_is_read_from_the_coding_notation_and_its_offset() {
    let noted = |extra: &str| {
        let rows = entries(&[row(
            "1:1000-1001",
            "A",
            "ENST1",
            &format!("IMPACT=MODIFIER;{extra}"),
            "146-147",
            "boundary",
        )]);
        utr_base_after_insertion(&rows[0])
    };
    assert_eq!(noted("HGVSc=ENST1.1:c.-55_-54insA"), Some(54));
    assert_eq!(noted("HGVSc=ENST1.1:c.-53_-52insA;HGVS_OFFSET=2"), Some(54));
    assert_eq!(
        noted("HGVSc=ENST1.1:c.-53_-52insA;HGVS_OFFSET=-2"),
        Some(54)
    );
    // A duplication stands for the copy written after the duplicated bases.
    assert_eq!(noted("HGVSc=ENST1.1:c.-56dup"), Some(55));
    assert_eq!(noted("HGVSc=ENST1.1:c.-58_-56dup"), Some(55));
    for none in [
        "HGVSc=ENST1.1:c.-55_-53insA",
        "HGVSc=ENST1.1:c.55_56insA",
        "HGVSc=ENST1.1:c.-55del",
        "HGVSc=c.-55_-54insA",
        "STRAND=1",
    ] {
        assert_eq!(noted(none), None, "{none}");
    }
}

#[test]
fn seeded_runs_are_read_plain_or_gzip_under_any_line_end() {
    let scratch = Scratch::new("reader");
    let reference = default_text(&[substitution(LOST_START)]);
    let candidate = default_text(&[substitution(LOST_STOP)]);
    let (e, a) = (
        &parse_default(&reference).unwrap()[0],
        &parse_default(&candidate).unwrap()[0],
    );
    let commented =
        format!("## ENSEMBL VARIANT EFFECT PREDICTOR\n## Column descriptions:\n{candidate}");
    let cases: Vec<(&str, Vec<u8>)> = vec![
        ("plain.txt", candidate.clone().into_bytes()),
        ("commented.txt", commented.clone().into_bytes()),
        ("packed.txt.gz", gzip(commented.as_bytes())),
        ("carriage.txt", candidate.replace('\n', "\r\n").into_bytes()),
        ("returns.txt", candidate.replace('\n', "\r").into_bytes()),
        (
            "unended.txt",
            candidate.trim_end_matches('\n').as_bytes().to_vec(),
        ),
        (
            "two-members.txt.gz",
            [
                gzip(HEADER.as_bytes()),
                gzip(substitution(LOST_STOP).as_bytes()),
            ]
            .concat(),
        ),
    ];
    let first = PathBuf::from(scratch.file("reference.txt", &reference));
    for (name, content) in cases {
        let second = PathBuf::from(scratch.file(name, content));
        let sweep = UtrannotatorSeedSweep::from_outputs(&[first.clone(), second]).unwrap();
        assert!(sweep.writes(e, a, true), "{name}");
    }
    // One output alone writes one term.
    let alone = UtrannotatorSeedSweep::from_outputs(&[first]).unwrap();
    assert!(alone.writes(e, e, true) && !alone.writes(e, a, false));
    // No output writes nothing.
    let nothing = UtrannotatorSeedSweep::from_outputs(&[]).unwrap();
    assert!(!nothing.writes(e, e, false));
}

#[test]
fn seeded_runs_refuse_a_file_that_is_no_default_format_output() {
    let scratch = Scratch::new("reader-refusals");
    let reference = default_text(&[substitution(LOST_START)]);
    let first = PathBuf::from(scratch.file("reference.txt", &reference));
    let not_an_output = "not a default-format output (columns Location, Allele, Feature)";
    let short_row = "rs1\t21:1\tT";
    let cases: Vec<(&str, Vec<u8>, String)> = vec![
        (
            "input.vcf",
            format!("{VCF_HEADER}21\t35748195\trs1\tG\tT\t.\t.\t.\n").into_bytes(),
            not_an_output.to_string(),
        ),
        ("empty.txt", Vec::new(), not_an_output.to_string()),
        ("empty.txt.gz", Vec::new(), not_an_output.to_string()),
        (
            "headless.txt",
            substitution(LOST_START).into_bytes(),
            not_an_output.to_string(),
        ),
        (
            "short-row.txt",
            format!("{reference}{short_row}\n").into_bytes(),
            format!("column count on line: {short_row}"),
        ),
        (
            "blank-row.txt",
            format!("{reference}\n").into_bytes(),
            "column count on line: ".to_string(),
        ),
        (
            "row-before-header.txt",
            format!("{short_row}\n{reference}").into_bytes(),
            format!("column count on line: {short_row}"),
        ),
        // The first header line decides, whatever a later one names.
        (
            "second-header.txt",
            format!("#CHROM\tPOS\n1\t2\n{reference}").into_bytes(),
            not_an_output.to_string(),
        ),
    ];
    for (name, content, message) in cases {
        let path = scratch.file(name, content);
        let read = UtrannotatorSeedSweep::from_outputs(&[first.clone(), PathBuf::from(&path)]);
        let expected = if message.starts_with("column count") {
            message
        } else {
            format!("{path}: {message}")
        };
        assert_eq!(refusal(read), expected, "{name}");
    }
}

#[test]
fn seeded_runs_refuse_a_file_that_cannot_be_read_before_any_of_its_rows() {
    let scratch = Scratch::new("reader-unread");
    let reference = default_text(&[substitution(LOST_START)]);
    let short_row = "rs1\t21:1\tT\n";
    let absent = scratch.path("absent.txt");
    assert_eq!(
        refusal(UtrannotatorSeedSweep::from_outputs(&[PathBuf::from(
            &absent
        )])),
        format!("cannot read {absent}: No such file or directory")
    );
    // Bytes that are no UTF-8 after a row that would be refused: the file is unread.
    let mut bytes = format!("{reference}{short_row}").into_bytes();
    let at = bytes.len();
    bytes.extend_from_slice(b"rs2\t21:2\t\xff\n");
    let ill_formed = scratch.file("ill-formed.txt", &bytes);
    assert_eq!(
        refusal(UtrannotatorSeedSweep::from_outputs(&[PathBuf::from(
            &ill_formed
        )])),
        format!(
            "cannot read {ill_formed}: 'utf-8' codec can't decode byte 0xff in position {}: invalid start byte",
            at + 9
        )
    );
    let not_packed = scratch.file("plain.txt.gz", &reference);
    assert_eq!(
        refusal(UtrannotatorSeedSweep::from_outputs(&[PathBuf::from(
            &not_packed
        )])),
        format!("cannot read {not_packed}: Not a gzipped file (b'#U')")
    );
    let packed = gzip(format!("{reference}{short_row}").as_bytes());
    let cut = scratch.file("cut.txt.gz", &packed[..packed.len() - 12]);
    assert_eq!(
        refusal(UtrannotatorSeedSweep::from_outputs(&[PathBuf::from(&cut)])),
        format!(
            "cannot read {cut}: Compressed file ended before the end-of-stream marker was reached"
        )
    );
}

/// The exit status, the lines of standard output and standard error of one command line.
fn run_main(arguments: &[&str]) -> (i32, Vec<String>, String) {
    let arguments: Vec<String> = arguments.iter().map(|word| word.to_string()).collect();
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = run(&arguments, &BuiltInRegistry, &BuiltIn, &mut out, &mut error);
    let out = String::from_utf8(out).unwrap();
    (
        code,
        out.lines().map(str::to_string).collect(),
        String::from_utf8(error).unwrap(),
    )
}

/// One scored run of a reference against the engine `engine`.
struct Scored {
    report: Value,
    residuals: PathBuf,
}

impl Scored {
    fn entry(&self, source: &str, projection: &str) -> &Value {
        let found: Vec<&Value> = self.report["scores"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|item| {
                item["source"] == source
                    && item["engine"] == "engine"
                    && item["projection"] == projection
            })
            .collect();
        assert_eq!(found.len(), 1, "{source} {projection}");
        found[0]
    }

    fn masked(&self, source: &str, projection: &str) -> Value {
        self.entry(source, projection)["masked"].clone()
    }

    /// The reference's values, the engine's and the pairs, raw and adjusted.
    fn counts(&self, source: &str) -> [u64; 6] {
        let item = self.entry(source, "feature");
        [
            "perl",
            "rust",
            "intersection",
            "adj_perl",
            "adj_rust",
            "adj_intersection",
        ]
        .map(|key| item[key].as_u64().unwrap())
    }

    fn residual_lines(&self, source: &str, projection: &str) -> Vec<String> {
        let path = self
            .residuals
            .join(format!("{source}.engine.{projection}.tsv"));
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(str::to_string)
            .collect()
    }
}

/// One run of the `score` command: its class table and plugins, and the run's input and
/// the seeded runs where given.
struct Run<'a> {
    name: &'a str,
    table: &'a str,
    plugins: &'a str,
    input: Option<&'a str>,
    sweep: &'a [&'a str],
}

/// Scores `engine` against `reference` as `run` says.
fn score(scratch: &Scratch, reference: &str, engine: &str, run: &Run) -> Scored {
    let Run {
        name,
        table,
        plugins,
        input,
        sweep,
    } = *run;
    let table = scratch.file(&format!("classes-{name}.tsv"), table);
    let report = scratch.path(&format!("report-{name}.json"));
    let residuals = scratch.path(&format!("residuals-{name}"));
    let work = scratch.path(&format!("work-{name}"));
    let engine = format!("engine={engine}");
    let mut line = vec![
        "--reference",
        reference,
        "--engine",
        &engine,
        "--plugins",
        plugins,
        "--class-table",
        &table,
        "--work-dir",
        &work,
        "--report",
        &report,
        "--residuals",
        &residuals,
        "--pairing-threads",
        "2",
    ];
    if let Some(input) = input {
        line.extend(["--input", input]);
    }
    for output in sweep {
        line.extend(["--utrannotator-sweep", output]);
    }
    let (code, _, error) = run_main(&line);
    assert_eq!((code, error.as_str()), (0, ""), "{name}");
    Scored {
        report: serde_json::from_str(&fs::read_to_string(&report).unwrap()).unwrap(),
        residuals: PathBuf::from(residuals),
    }
}

fn sorted(mut lines: Vec<String>) -> Vec<String> {
    lines.sort();
    lines
}

/// The residual lines of one side's row: one per field of `fields`, each `name=value`
/// pairs joined by `;`, in the order of the names, `LOEUF` left out. `feature` is the
/// row's Feature with its tab, or empty for the allele projection.
fn left(side: &str, location: &str, allele: &str, feature: &str, fields: &[&str]) -> Vec<String> {
    let mut pairs: BTreeMap<&str, &str> = BTreeMap::new();
    for pair in fields.iter().flat_map(|text| text.split(';')) {
        let (name, value) = pair.split_once('=').unwrap();
        pairs.insert(name, value);
    }
    pairs
        .iter()
        .filter(|(name, _)| **name != "LOEUF")
        .map(|(name, value)| format!("{side}\t{location}\t{allele}\t{feature}{name}\t{value}"))
        .collect()
}

/// The committed table's row of the structural class, which declares the plugin's two
/// fields met.
fn structural_table() -> String {
    let declared: Vec<&str> = COMMITTED_CLASS_TABLE
        .lines()
        .filter(|line| line.starts_with(&format!("{STRUCTURAL}\t")))
        .collect();
    assert_eq!(declared.len(), 1);
    let columns: Vec<&str> = declared[0].split('\t').collect();
    assert_eq!(
        columns[1..4],
        ["UTRAnnotator", "5UTR_consequence,5UTR_annotation", "met"]
    );
    format!("{}\n", declared[0])
}

const SPLICED: &str =
    "5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;5UTR_annotation=alt_type=uORF:ref_type=uORF";
const SPLICED_REORDERED: &str =
    "5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;5UTR_annotation=ref_type=uORF:alt_type=uORF";
const STRUCTURAL_LOST_START: &str =
    "5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;5UTR_annotation=type=uORF";

/// A 5' UTR row under `allele` whose Extra carries `fields` after its IMPACT.
fn utr_structural_row(location: &str, allele: &str, fields: &[&str]) -> String {
    let extra = [&["IMPACT=MODIFIER"], fields].concat().join(";");
    row(location, allele, "ENST1", &extra, "-", "sv")
}

/// The two fields the reference alone fills on a row whose Allele is a structural class
/// word are set aside, on each row of a key two rows share, with the seeded runs given or
/// not. A sequence allele's row, a breakend's own string, a row the engine fills, a value
/// with no engine row left to stand against, the plugin's `Existing_uORFs` and another
/// plugin's field stay residuals.
fn structural_allele_rows_are_set_aside(seeded: bool) {
    let scratch = Scratch::new(if seeded {
        "structural-seeded"
    } else {
        "structural"
    });
    let (existing, two_rows, one_row_less) = ("Existing_uORFs=6", "21:7000-7100", "21:8000-8100");
    let reference = scratch.output(
        "reference.txt",
        &[
            utr_structural_row("21:1000-1040", "deletion", &[SPLICED, existing]),
            utr_structural_row("21:2000-2100", "duplication", &[SPLICED]),
            utr_structural_row("21:3000-3004", "-", &[SPLICED]),
            utr_structural_row("21:4000-4001", "insertion", &[]),
            utr_structural_row("21:5000-5100", "inversion", &[SPLICED, existing]),
            utr_structural_row("21:6000-6001", "Alu_insertion", &[SPLICED, "LOEUF=0.31"]),
            utr_structural_row(two_rows, "copy_number_variation", &[SPLICED]),
            utr_structural_row(two_rows, "copy_number_variation", &[SPLICED_REORDERED]),
            utr_structural_row(one_row_less, "copy_number_variation", &[SPLICED]),
            utr_structural_row(one_row_less, "copy_number_variation", &[SPLICED]),
            utr_structural_row("21:9000", "N[21:9500[", &[SPLICED]),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            utr_structural_row("21:1000-1040", "deletion", &[existing]),
            utr_structural_row("21:2000-2100", "duplication", &[STRUCTURAL_LOST_START]),
            utr_structural_row("21:3000-3004", "-", &[]),
            utr_structural_row("21:4000-4001", "insertion", &[SPLICED]),
            utr_structural_row("21:5000-5100", "inversion", &[]),
            utr_structural_row("21:6000-6001", "Alu_insertion", &[]),
            utr_structural_row(two_rows, "copy_number_variation", &[]),
            utr_structural_row(two_rows, "copy_number_variation", &[]),
            utr_structural_row(one_row_less, "copy_number_variation", &[]),
            utr_structural_row("21:9000", "N[21:9500[", &[]),
        ],
    );
    let sweep: Vec<&str> = if seeded { vec![&reference] } else { Vec::new() };
    let scored = score(
        &scratch,
        &reference,
        &engine,
        &Run {
            name: "run",
            table: &structural_table(),
            plugins: "UTRAnnotator,LOEUF",
            input: None,
            sweep: &sweep,
        },
    );
    let stays: [(&str, &str, &str, &[&str]); 7] = [
        (
            "engine",
            "21:2000-2100",
            "duplication",
            &[STRUCTURAL_LOST_START],
        ),
        ("reference", "21:2000-2100", "duplication", &[SPLICED]),
        ("reference", "21:3000-3004", "-", &[SPLICED]),
        ("engine", "21:4000-4001", "insertion", &[SPLICED]),
        (
            "reference",
            "21:5000-5100",
            "inversion",
            &["Existing_uORFs=6"],
        ),
        (
            "reference",
            one_row_less,
            "copy_number_variation",
            &[SPLICED],
        ),
        ("reference", "21:9000", "N[21:9500[", &[SPLICED]),
    ];
    let lines = |feature: &str| -> Vec<String> {
        stays
            .iter()
            .flat_map(|(side, location, allele, fields)| {
                left(side, location, allele, feature, fields)
            })
            .collect()
    };
    assert_eq!(scored.counts("UTRAnnotator"), [22, 5, 1, 10, 5, 1]);
    assert_eq!(
        scored.masked("UTRAnnotator", "feature"),
        json!({STRUCTURAL: 10})
    );
    assert_eq!(
        sorted(scored.residual_lines("UTRAnnotator", "feature")),
        sorted(lines("ENST1\t"))
    );
    assert_eq!(
        scored.masked("UTRAnnotator", "allele"),
        json!({STRUCTURAL: 8})
    );
    assert_eq!(
        sorted(scored.residual_lines("UTRAnnotator", "allele")),
        sorted(lines(""))
    );
    assert_eq!(scored.masked("LOEUF", "feature"), json!({}));
    assert_eq!(
        scored.residual_lines("LOEUF", "feature"),
        ["reference\t21:6000-6001\tAlu_insertion\tENST1\tLOEUF\t0.31"]
    );
}

#[test]
fn fields_on_a_structural_allele_row_are_set_aside_without_the_seeded_runs() {
    structural_allele_rows_are_set_aside(false);
}

#[test]
fn fields_on_a_structural_allele_row_are_set_aside_beside_the_seeded_runs() {
    structural_allele_rows_are_set_aside(true);
}

const SPAN_EFFECT: &str =
    "5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;5UTR_annotation=AltStop=True:FrameWithCDS=inFrame";
const SPAN_EFFECT_OTHER: &str =
    "5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;5UTR_annotation=AltStop=False:FrameWithCDS=inFrame";

/// A deletion's row: its Location, Allele and cDNA positions, the reference's fields, the
/// engine's fields and the Consequence terms.
type DeletionShape<'a> = (
    &'a str,
    &'a str,
    &'a str,
    Vec<&'a str>,
    Vec<&'a str>,
    &'a str,
);
/// An insertion's row: its Location and Allele, the reference's fields, the engine's
/// fields and the Consequence terms.
type InsertionShape<'a> = (&'a str, &'a str, Vec<&'a str>, Vec<&'a str>, &'a str);

/// A structural record's row under `allele` with its cDNA positions and the Consequence
/// `terms`, whose Extra carries `fields` after its IMPACT.
fn utr_span_row(location: &str, allele: &str, cdna: &str, fields: &[&str], terms: &str) -> String {
    let extra = [&["IMPACT=MODIFIER"], fields].concat().join(";");
    let text = row(location, allele, "ENST1", &extra, cdna, "sv");
    assert_eq!(text.matches("\tmissense_variant\t").count(), 1);
    text.replace("\tmissense_variant\t", &format!("\t{terms}\t"))
}

/// The two fields the engine fills on the row of a symbolic deletion inside one exon of
/// the 5' UTR are set aside, beside an empty reference, beside the reference's own term
/// and, where the term is equal, for the annotation alone. A row the engine fills under
/// another class word, across two exons, into the intron, with cDNA positions open at an
/// end, with the gained stop or with a frameshift on a span that is a multiple of three
/// stays a residual, as do the plugin's `Existing_uORFs` and another plugin's field.
fn deletion_inside_one_utr_exon_is_set_aside(seeded: bool) {
    let scratch = Scratch::new(if seeded {
        "deletion-seeded"
    } else {
        "deletion"
    });
    let stop_gain = format!("5UTR_consequence={STOP_GAIN};5UTR_annotation=ref_type=uORF");
    let frameshift = format!("5UTR_consequence={FRAMESHIFT};5UTR_annotation=ref_type=uORF");
    let within = "feature_truncation,5_prime_UTR_variant";
    let shapes: Vec<DeletionShape> = vec![
        (
            "21:1001-1040",
            "deletion",
            "101-140",
            vec![],
            vec![SPAN_EFFECT],
            within,
        ),
        (
            "21:2001-2040",
            "deletion",
            "101-140",
            vec![SPLICED],
            vec![SPAN_EFFECT],
            within,
        ),
        (
            "21:3001-3040",
            "mobile_element_deletion",
            "101-140",
            vec![SPAN_EFFECT_OTHER],
            vec![SPAN_EFFECT],
            within,
        ),
        (
            "21:4001-4040",
            "deletion",
            "101-140",
            vec!["Existing_uORFs=6"],
            vec![SPAN_EFFECT, "LOEUF=0.31"],
            within,
        ),
        (
            "21:5001-5040",
            "duplication",
            "101-140",
            vec![],
            vec![SPAN_EFFECT],
            within,
        ),
        (
            "21:6001-6040",
            "deletion",
            "101-120",
            vec![],
            vec![SPAN_EFFECT],
            within,
        ),
        (
            "21:7001-7040",
            "deletion",
            "?-140",
            vec![SPLICED],
            vec![SPAN_EFFECT],
            within,
        ),
        (
            "21:8001-8040",
            "deletion",
            "101-140",
            vec![],
            vec![SPAN_EFFECT],
            "feature_truncation,5_prime_UTR_variant,intron_variant",
        ),
        (
            "21:9001-9040",
            "deletion",
            "101-140",
            vec![],
            vec![&stop_gain],
            within,
        ),
        (
            "21:9101-9139",
            "deletion",
            "101-139",
            vec![],
            vec![&frameshift],
            within,
        ),
        (
            "21:9201-9240",
            "deletion",
            "101-140",
            vec![],
            vec![&frameshift],
            within,
        ),
        (
            "21:9301-9340",
            "deletion",
            "101-140",
            vec![SPLICED],
            vec![],
            within,
        ),
    ];
    let side = |engine: bool| -> Vec<String> {
        shapes
            .iter()
            .map(|(location, allele, cdna, reference, written, terms)| {
                let fields = if engine { written } else { reference };
                utr_span_row(location, allele, cdna, fields, terms)
            })
            .collect()
    };
    let reference = scratch.output("reference.txt", &side(false));
    let engine = scratch.output("engine.txt", &side(true));
    let sweep: Vec<&str> = if seeded { vec![&reference] } else { Vec::new() };
    let scored = score(
        &scratch,
        &reference,
        &engine,
        &Run {
            name: "run",
            table: &structural_table(),
            plugins: "UTRAnnotator,LOEUF",
            input: None,
            sweep: &sweep,
        },
    );
    let stays: [(&str, &str, &str, &[&str]); 8] = [
        (
            "reference",
            "21:4001-4040",
            "deletion",
            &["Existing_uORFs=6"],
        ),
        ("engine", "21:5001-5040", "duplication", &[SPAN_EFFECT]),
        ("engine", "21:6001-6040", "deletion", &[SPAN_EFFECT]),
        ("engine", "21:7001-7040", "deletion", &[SPAN_EFFECT]),
        ("reference", "21:7001-7040", "deletion", &[SPLICED]),
        ("engine", "21:8001-8040", "deletion", &[SPAN_EFFECT]),
        ("engine", "21:9001-9040", "deletion", &[&stop_gain]),
        ("engine", "21:9101-9139", "deletion", &[&frameshift]),
    ];
    let lines: Vec<String> = stays
        .iter()
        .flat_map(|(side, location, allele, fields)| {
            left(side, location, allele, "ENST1\t", fields)
        })
        .collect();
    // Set aside: two fields each on the rows at 1001, 2001, 4001, 9201 and 9301, the
    // annotation at 3001.
    assert_eq!(
        scored.masked("UTRAnnotator", "feature"),
        json!({STRUCTURAL: 11})
    );
    assert_eq!(
        sorted(scored.residual_lines("UTRAnnotator", "feature")),
        sorted(lines)
    );
    assert_eq!(
        scored.residual_lines("LOEUF", "feature"),
        ["engine\t21:4001-4040\tdeletion\tENST1\tLOEUF\t0.31"]
    );
}

#[test]
fn fields_the_engine_fills_on_a_deletion_inside_one_utr_exon_are_set_aside_without_the_seeded_runs()
{
    deletion_inside_one_utr_exon_is_set_aside(false);
}

#[test]
fn fields_the_engine_fills_on_a_deletion_inside_one_utr_exon_are_set_aside_beside_the_seeded_runs()
{
    deletion_inside_one_utr_exon_is_set_aside(true);
}

const UTR_FRAMESHIFT: &str = "5UTR_consequence=5_prime_UTR_uORF_frameshift_variant";
const REFERENCE_KEYS: &str = "5UTR_annotation=Evidence=False:KozakContext=GGGATGC:KozakStrength=Moderate:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51";
const REFERENCE_KEYS_AND_TYPE: &str = "5UTR_annotation=Evidence=False:KozakContext=GGGATGC:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51";
const TEXT_LENGTH: &str = "5UTR_annotation=alt_type_length=39:KozakContext=GGGATGC:KozakStrength=Moderate:ref_type_length=51:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=207";

/// On the row of a plain symbolic insertion whose stated length is no multiple of three,
/// the annotation the engine writes beside the frameshift, of the six keys read from the
/// reference UTR with or without `alt_type=uORF` in each effect, is set aside beside the
/// reference's own frameshift, and both fields beside another of the plugin's terms and
/// beside an empty reference, whatever other Consequence terms the row names. An effect
/// with `alt_type_length`, another `alt_type`, five of the six keys or an eighth key,
/// another term, the frameshift beside no annotation, a span that is a multiple of three,
/// a mobile element's word, a duplication, an inversion and a row without
/// `5_prime_UTR_variant` stay residuals.
fn insertion_of_a_stated_length_is_set_aside(seeded: bool) {
    let scratch = Scratch::new(if seeded {
        "insertion-seeded"
    } else {
        "insertion"
    });
    let text_length = vec![UTR_FRAMESHIFT, TEXT_LENGTH];
    let six = vec![UTR_FRAMESHIFT, REFERENCE_KEYS];
    let six_and_type = vec![UTR_FRAMESHIFT, REFERENCE_KEYS_AND_TYPE];
    let two_effects = format!(
        "{REFERENCE_KEYS},{}",
        REFERENCE_KEYS_AND_TYPE.split_once('=').unwrap().1
    );
    let stop_loss = "5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant";
    let other_effects: [(&str, String); 4] = [
        (
            "21:6001-6113",
            format!("{REFERENCE_KEYS_AND_TYPE}:alt_type_length=66"),
        ),
        (
            "21:6201-6313",
            REFERENCE_KEYS.replace(
                ":ref_StartDistanceToCDS",
                ":alt_type=OutOfFrame_oORF:ref_StartDistanceToCDS",
            ),
        ),
        (
            "21:6401-6513",
            REFERENCE_KEYS.replace("Evidence=False:", ""),
        ),
        (
            "21:6601-6713",
            format!("{REFERENCE_KEYS_AND_TYPE}:type=uORF"),
        ),
    ];
    let other_words = [
        ("21:9001-9113", "Alu_insertion"),
        ("21:9201-9313", "mobile_element_insertion"),
        ("21:9401-9513", "duplication"),
        ("21:9601-9713", "inversion"),
    ];
    let utr = "feature_elongation,5_prime_UTR_variant";
    let nmd = "feature_elongation,5_prime_UTR_variant,NMD_transcript_variant";
    let intron = "feature_elongation,5_prime_UTR_variant,intron_variant";
    let no_utr = "feature_elongation,intron_variant";
    let mut shapes: Vec<InsertionShape> = vec![
        (
            "21:1001-1113",
            "insertion",
            text_length.clone(),
            six.clone(),
            utr,
        ),
        (
            "21:2001-2113",
            "insertion",
            text_length.clone(),
            vec![UTR_FRAMESHIFT, &two_effects],
            utr,
        ),
        (
            "21:3001-3113",
            "insertion",
            vec![STRUCTURAL_LOST_START],
            six.clone(),
            utr,
        ),
        ("21:4001-4113", "insertion", vec![], six_and_type, utr),
        (
            "21:5001-5113",
            "insertion",
            text_length.clone(),
            six.clone(),
            nmd,
        ),
        (
            "21:5201-5313",
            "insertion",
            text_length.clone(),
            six.clone(),
            intron,
        ),
    ];
    for (location, effect) in &other_effects {
        shapes.push((
            location,
            "insertion",
            text_length.clone(),
            vec![UTR_FRAMESHIFT, effect],
            utr,
        ));
    }
    shapes.push((
        "21:7001-7113",
        "insertion",
        text_length.clone(),
        vec![stop_loss, REFERENCE_KEYS],
        utr,
    ));
    shapes.push((
        "21:7201-7313",
        "insertion",
        text_length.clone(),
        vec![UTR_FRAMESHIFT],
        utr,
    ));
    shapes.push((
        "21:8001-8114",
        "insertion",
        text_length.clone(),
        six.clone(),
        utr,
    ));
    for (location, word) in other_words {
        shapes.push((location, word, text_length.clone(), six.clone(), utr));
    }
    shapes.push((
        "21:9801-9913",
        "insertion",
        text_length.clone(),
        six.clone(),
        no_utr,
    ));
    let side = |engine: bool| -> Vec<String> {
        shapes
            .iter()
            .map(|(location, allele, reference, written, terms)| {
                let fields = if engine { written } else { reference };
                utr_span_row(location, allele, "101-102", fields, terms)
            })
            .collect()
    };
    let reference = scratch.output("reference.txt", &side(false));
    let engine = scratch.output("engine.txt", &side(true));
    let sweep: Vec<&str> = if seeded { vec![&reference] } else { Vec::new() };
    let scored = score(
        &scratch,
        &reference,
        &engine,
        &Run {
            name: "run",
            table: &structural_table(),
            plugins: "UTRAnnotator",
            input: None,
            sweep: &sweep,
        },
    );
    // Beside the reference's frameshift the consequence is shared, and the annotations
    // stay.
    let mut stays: Vec<(&str, &str, &str, Vec<&str>)> = Vec::new();
    for (location, effect) in &other_effects {
        stays.push(("reference", location, "insertion", vec![TEXT_LENGTH]));
        stays.push(("engine", location, "insertion", vec![effect]));
    }
    stays.push((
        "reference",
        "21:7001-7113",
        "insertion",
        text_length.clone(),
    ));
    stays.push((
        "engine",
        "21:7001-7113",
        "insertion",
        vec![stop_loss, REFERENCE_KEYS],
    ));
    stays.push(("reference", "21:7201-7313", "insertion", vec![TEXT_LENGTH]));
    stays.push(("reference", "21:8001-8114", "insertion", vec![TEXT_LENGTH]));
    stays.push(("engine", "21:8001-8114", "insertion", vec![REFERENCE_KEYS]));
    for (location, word) in other_words {
        stays.push(("reference", location, word, vec![TEXT_LENGTH]));
        stays.push(("engine", location, word, vec![REFERENCE_KEYS]));
    }
    stays.push(("reference", "21:9801-9913", "insertion", vec![TEXT_LENGTH]));
    stays.push(("engine", "21:9801-9913", "insertion", vec![REFERENCE_KEYS]));
    let lines: Vec<String> = stays
        .iter()
        .flat_map(|(side, location, allele, fields)| {
            left(side, location, allele, "ENST1\t", fields)
        })
        .collect();
    // Set aside: the annotation on the rows at 1001, 2001, 5001 and 5201, two fields each
    // at 3001 and 4001.
    assert_eq!(
        scored.masked("UTRAnnotator", "feature"),
        json!({STRUCTURAL: 8})
    );
    assert_eq!(
        sorted(scored.residual_lines("UTRAnnotator", "feature")),
        sorted(lines)
    );
}

#[test]
fn fields_the_engine_fills_on_an_insertion_of_a_stated_length_are_set_aside_without_the_seeded_runs(
) {
    insertion_of_a_stated_length_is_set_aside(false);
}

#[test]
fn fields_the_engine_fills_on_an_insertion_of_a_stated_length_are_set_aside_beside_the_seeded_runs()
{
    insertion_of_a_stated_length_is_set_aside(true);
}

/// The classes that read the plugin's two fields, in the committed table's order, then
/// those of `annotation_alone`, which read `5UTR_annotation` alone.
fn class_table(two_fields: &[&str], annotation_alone: &[&str]) -> String {
    let mut table = String::new();
    for klass in two_fields {
        table.push_str(&format!(
            "{klass}\tUTRAnnotator\t5UTR_consequence,5UTR_annotation\tmet\t-\n"
        ));
        if *klass == REMADE {
            for klass in annotation_alone {
                table.push_str(&format!("{klass}\tUTRAnnotator\t5UTR_annotation\tmet\t-\n"));
            }
        }
    }
    table
}

/// The Extra of an insertion's row: the two fields and the plugin's count, after the
/// IMPACT and the strand where the row states one.
fn utr_fields(strand: Option<&str>, effect: (&str, &str)) -> String {
    let (consequence, annotation) = effect;
    let strand = strand.map_or(String::new(), |strand| format!("STRAND={strand};"));
    format!("IMPACT=MODIFIER;{strand}5UTR_annotation={annotation};5UTR_consequence={consequence};Existing_uORFs=1")
}

fn at_index(index: u32) -> String {
    format!("CapDistanceToStart={index}:type=uORF")
}

const SHIFT: &str = "alt_type=uORF:ref_type=uORF";
const NO_EFFECT: (&str, &str) = ("-", "-");

/// One insertion on a transcript of its own: the record's identifier, the row's Location,
/// Allele, strand, Feature and cDNA positions, and the two fields of each side.
struct Insertion<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    strand: Option<&'a str>,
    feature: &'a str,
    cdna: &'a str,
    reference: (&'a str, String),
    engine: (&'a str, String),
}

impl Insertion<'_> {
    fn line(&self, effect: (&str, &str)) -> String {
        row(
            self.location,
            self.allele,
            self.feature,
            &utr_fields(self.strand, effect),
            self.cdna,
            self.name,
        )
    }

    fn reference_line(&self) -> String {
        self.line((self.reference.0, &self.reference.1))
    }

    fn engine_line(&self) -> String {
        self.line((self.engine.0, &self.engine.1))
    }
}

/// The insertion of `record` (its identifier, Location, Allele, Feature and cDNA
/// positions) on a row that states `strand`, with the two fields of each side.
fn insertion<'a>(
    record: [&'a str; 5],
    strand: Option<&'a str>,
    reference: (&'a str, &str),
    engine: (&'a str, &str),
) -> Insertion<'a> {
    let [name, location, allele, feature, cdna] = record;
    Insertion {
        name,
        location,
        allele,
        strand,
        feature,
        cdna,
        reference: (reference.0, reference.1.to_string()),
        engine: (engine.0, engine.1.to_string()),
    }
}

/// The two outputs of `rows` and the seeded runs of `seeds`, each its rows.
fn insertion_outputs(
    scratch: &Scratch,
    rows: &[Insertion],
    seeds: [Vec<String>; 2],
) -> (String, String, [String; 2]) {
    let reference: Vec<String> = rows.iter().map(Insertion::reference_line).collect();
    let engine: Vec<String> = rows.iter().map(Insertion::engine_line).collect();
    let [seed_a, seed_b] = seeds;
    (
        scratch.output("reference.txt", &reference),
        scratch.output("engine.txt", &engine),
        [
            scratch.output("seed-a.txt", &seed_a),
            scratch.output("seed-b.txt", &seed_b),
        ],
    )
}

#[test]
fn lost_start_after_an_insertion_is_set_aside_after_the_window_and_the_precedence_classes() {
    let scratch = Scratch::new("lost-start");
    // On five rows the reference reports lost the start codon at the row's lower cDNA
    // position, the codon directly after the insertion: the engine reports nothing
    // (`whole`), a frameshift the seeded runs write on the row too (`swept`), a frameshift
    // on a row no seeded run holds (`unswept`), and a start gained inside its allele under
    // an identifier one input line carries (`gaining`) and under one that two carry
    // (`shared`). On `parted` the reference reports the codon one base before that
    // position, on a row that states no strand, and the engine reports nothing.
    let gained = (GAIN, at_index(302));
    let shifted = (FRAMESHIFT, SHIFT.to_string());
    let rows = [
        insertion(
            ["whole", "1:1000-1001", "CCA", "ENST1", "146-147"],
            None,
            (LOSS, &at_index(146)),
            NO_EFFECT,
        ),
        insertion(
            ["swept", "1:2000-2001", "C", "ENST2", "40-41"],
            None,
            (LOSS, &at_index(40)),
            (shifted.0, &shifted.1),
        ),
        insertion(
            ["unswept", "1:3000-3001", "C", "ENST3", "50-51"],
            None,
            (LOSS, &at_index(50)),
            (shifted.0, &shifted.1),
        ),
        insertion(
            ["gaining", "1:4000-4001", "CCATGC", "ENST4", "300-301"],
            None,
            (LOSS, &at_index(300)),
            (gained.0, &gained.1),
        ),
        insertion(
            ["shared", "1:5000-5001", "CCATGC", "ENST5", "300-301"],
            None,
            (LOSS, &at_index(300)),
            (gained.0, &gained.1),
        ),
        insertion(
            ["parted", "1:6000-6001", "CCA", "ENST6", "146-147"],
            None,
            (LOSS, &at_index(145)),
            NO_EFFECT,
        ),
    ];
    let seeds = [vec![rows[1].reference_line()], vec![rows[1].engine_line()]];
    let (reference, engine, [seed_a, seed_b]) = insertion_outputs(&scratch, &rows, seeds);
    let input = scratch.vcf(
        "input.vcf",
        &[
            ["1", "1000", "whole", "G", "GCCA"],
            ["1", "2000", "swept", "G", "GC"],
            ["1", "3000", "unswept", "G", "GC"],
            ["1", "4000", "gaining", "G", "GCCATGC"],
            ["1", "5000", "shared", "G", "GCCATGC"],
            ["1", "5000", "shared", "G", "GT"],
            ["1", "6000", "parted", "G", "GCCA"],
        ],
    );
    let table = class_table(&[WINDOW, PRECEDENCE, LOST_AFTER], &[]);
    let parted = [
        "reference\t1:6000-6001\tCCA\tENST6\t5UTR_annotation\tCapDistanceToStart=145:type=uORF"
            .to_string(),
        format!("reference\t1:6000-6001\tCCA\tENST6\t5UTR_consequence\t{LOSS}"),
    ];
    // With the run's input and the seeded runs, the two rows whose gained start the window
    // class explains from the record's alleles and the row whose two terms the seeded runs
    // write are those classes'; without them all five are this class's.
    let cases = [
        (
            "given",
            Some(input.as_str()),
            vec![seed_a.as_str(), seed_b.as_str()],
            json!({WINDOW: 4, PRECEDENCE: 2, LOST_AFTER: 4}),
        ),
        (
            "without",
            None,
            vec![],
            json!({WINDOW: 0, PRECEDENCE: 0, LOST_AFTER: 10}),
        ),
    ];
    for (name, given, sweep, masked) in cases {
        let scored = score(
            &scratch,
            &reference,
            &engine,
            &Run {
                name,
                table: &table,
                plugins: "UTRAnnotator",
                input: given,
                sweep: &sweep,
            },
        );
        assert_eq!(scored.masked("UTRAnnotator", "feature"), masked, "{name}");
        assert_eq!(
            scored.counts("UTRAnnotator"),
            [18, 14, 6, 8, 6, 6],
            "{name}"
        );
        assert_eq!(
            scored.residual_lines("UTRAnnotator", "feature"),
            parted,
            "{name}"
        );
    }
}

#[test]
fn lost_start_at_a_codon_an_insertion_only_moves_on_is_set_aside() {
    let scratch = Scratch::new("moved-start");
    // On every row the reference reports lost a start codon one or two bases before the
    // row's lower cDNA position, the insertion standing inside it. On five the allele,
    // read on the transcript's strand, ends with the codon's first base or its first two:
    // the engine reports nothing (`moved`, and `moved_two` at two bases), a frameshift the
    // seeded runs write on the row too, on a reverse-strand transcript (`swept`), a start
    // gained inside its allele (`making`), and a start gained at the lost codon's index
    // plus the allele's length, the moved codon itself (`named`). On `parted` the allele
    // ends with neither, and `strandless` is `moved` on a row that states no strand.
    let forward = Some("1");
    let rows = [
        insertion(
            ["moved", "1:1000-1001", "CCA", "ENST1", "146-147"],
            forward,
            (LOSS, &at_index(145)),
            NO_EFFECT,
        ),
        insertion(
            ["swept", "1:2000-2001", "TGG", "ENST2", "146-147"],
            Some("-1"),
            (LOSS, &at_index(145)),
            (FRAMESHIFT, SHIFT),
        ),
        insertion(
            ["moved_two", "1:3000-3001", "CAT", "ENST3", "40-41"],
            forward,
            (LOSS, &at_index(38)),
            NO_EFFECT,
        ),
        insertion(
            ["making", "1:4000-4001", "CATGA", "ENST4", "300-301"],
            forward,
            (LOSS, &at_index(299)),
            (GAIN, &at_index(301)),
        ),
        insertion(
            ["named", "1:5000-5001", "CCA", "ENST5", "300-301"],
            forward,
            (LOSS, &at_index(299)),
            (GAIN, &at_index(302)),
        ),
        insertion(
            ["parted", "1:6000-6001", "CCC", "ENST6", "146-147"],
            forward,
            (LOSS, &at_index(145)),
            NO_EFFECT,
        ),
        insertion(
            ["strandless", "1:7000-7001", "CCA", "ENST7", "146-147"],
            None,
            (LOSS, &at_index(145)),
            NO_EFFECT,
        ),
    ];
    let seeds = [vec![rows[1].reference_line()], vec![rows[1].engine_line()]];
    let (reference, engine, [seed_a, seed_b]) = insertion_outputs(&scratch, &rows, seeds);
    let input = scratch.vcf(
        "input.vcf",
        &[
            ["1", "1000", "moved", "G", "GCCA"],
            ["1", "2000", "swept", "G", "GTGG"],
            ["1", "3000", "moved_two", "G", "GCAT"],
            ["1", "4000", "making", "G", "GCATGA"],
            ["1", "5000", "named", "G", "GCCA"],
            ["1", "6000", "parted", "G", "GCCC"],
            ["1", "7000", "strandless", "G", "GCCA"],
        ],
    );
    let table = class_table(&[WINDOW, PRECEDENCE, LOST_AFTER], &[]);
    let stays = [
        "reference\t1:6000-6001\tCCC\tENST6\t5UTR_annotation\tCapDistanceToStart=145:type=uORF"
            .to_string(),
        format!("reference\t1:6000-6001\tCCC\tENST6\t5UTR_consequence\t{LOSS}"),
        "reference\t1:7000-7001\tCCA\tENST7\t5UTR_annotation\tCapDistanceToStart=145:type=uORF"
            .to_string(),
        format!("reference\t1:7000-7001\tCCA\tENST7\t5UTR_consequence\t{LOSS}"),
    ];
    // With the run's input and the seeded runs, the gained start of a codon the inserted
    // bases make is the window class's and the frameshift the seeded runs write the
    // precedence class's; the gained start that names the moved codon, which the window
    // class refuses, is this class's.
    let cases = [
        (
            "given",
            Some(input.as_str()),
            vec![seed_a.as_str(), seed_b.as_str()],
            json!({WINDOW: 2, PRECEDENCE: 2, LOST_AFTER: 6}),
        ),
        (
            "without",
            None,
            vec![],
            json!({WINDOW: 0, PRECEDENCE: 0, LOST_AFTER: 10}),
        ),
    ];
    for (name, given, sweep, masked) in cases {
        let scored = score(
            &scratch,
            &reference,
            &engine,
            &Run {
                name,
                table: &table,
                plugins: "UTRAnnotator",
                input: given,
                sweep: &sweep,
            },
        );
        assert_eq!(scored.masked("UTRAnnotator", "feature"), masked, "{name}");
        assert_eq!(
            scored.counts("UTRAnnotator"),
            [21, 13, 7, 11, 7, 7],
            "{name}"
        );
        assert_eq!(
            scored.residual_lines("UTRAnnotator", "feature"),
            stays,
            "{name}"
        );
    }
}

#[test]
fn gained_start_of_a_codon_left_at_its_own_index_is_set_aside() {
    let scratch = Scratch::new("remade-start");
    // On five rows the inserted bases, read on the transcript's strand, open with the last
    // bases of a start codon one or two bases before the row's lower cDNA position and
    // leave it at its own index. The reference reports that codon gained on `kept` (the
    // engine reports nothing), on `kept_rev`, a reverse-strand transcript (the engine a
    // frameshift the seeded runs write on the row too), on `shifted` at two bases (the
    // engine a frameshift no seeded run holds) and on `making` (the engine a start gained
    // further on, inside its allele). On `swept` both report the frameshift, the engine's
    // listing one ORF more, and a seeded run reports the gained start. On `parted` the
    // allele opens with neither, and `strandless` is `kept` on a row that states no strand.
    let own = format!("{SHIFT},alt_type=inFrame_oORF:ref_type=uORF");
    let forward = Some("1");
    let rows = [
        insertion(
            ["kept", "1:1000-1001", "TGC", "ENST1", "146-147"],
            forward,
            (GAIN, &at_index(145)),
            NO_EFFECT,
        ),
        insertion(
            ["kept_rev", "1:2000-2001", "GCA", "ENST2", "146-147"],
            Some("-1"),
            (GAIN, &at_index(145)),
            (FRAMESHIFT, SHIFT),
        ),
        insertion(
            ["shifted", "1:3000-3001", "G", "ENST3", "40-41"],
            forward,
            (GAIN, &at_index(38)),
            (FRAMESHIFT, &own),
        ),
        insertion(
            ["swept", "1:4000-4001", "TG", "ENST4", "300-301"],
            forward,
            (FRAMESHIFT, SHIFT),
            (FRAMESHIFT, &own),
        ),
        insertion(
            ["making", "1:5000-5001", "TGATGC", "ENST5", "300-301"],
            forward,
            (GAIN, &at_index(299)),
            (GAIN, &at_index(302)),
        ),
        insertion(
            ["parted", "1:6000-6001", "CCC", "ENST6", "146-147"],
            forward,
            (GAIN, &at_index(145)),
            NO_EFFECT,
        ),
        insertion(
            ["strandless", "1:7000-7001", "TGC", "ENST7", "146-147"],
            None,
            (GAIN, &at_index(145)),
            NO_EFFECT,
        ),
    ];
    let seeds = [
        vec![rows[1].reference_line(), rows[3].reference_line()],
        vec![rows[1].engine_line(), rows[3].line((GAIN, &at_index(299)))],
    ];
    let (reference, engine, [seed_a, seed_b]) = insertion_outputs(&scratch, &rows, seeds);
    let input = scratch.vcf(
        "input.vcf",
        &[
            ["1", "1000", "kept", "G", "GTGC"],
            ["1", "2000", "kept_rev", "G", "GGCA"],
            ["1", "3000", "shifted", "G", "GG"],
            ["1", "4000", "swept", "G", "GTG"],
            ["1", "5000", "making", "G", "GTGATGC"],
            ["1", "6000", "parted", "G", "GCCC"],
            ["1", "7000", "strandless", "G", "GTGC"],
        ],
    );
    let table = class_table(&[WINDOW, PRECEDENCE, LOST_AFTER, REMADE], &[]);
    let stays = vec![
        "reference\t1:6000-6001\tCCC\tENST6\t5UTR_annotation\tCapDistanceToStart=145:type=uORF"
            .to_string(),
        format!("reference\t1:6000-6001\tCCC\tENST6\t5UTR_consequence\t{GAIN}"),
        "reference\t1:7000-7001\tTGC\tENST7\t5UTR_annotation\tCapDistanceToStart=145:type=uORF"
            .to_string(),
        format!("reference\t1:7000-7001\tTGC\tENST7\t5UTR_consequence\t{GAIN}"),
    ];
    let swept = vec![
        "reference\t1:4000-4001\tTG\tENST4\t5UTR_annotation\talt_type=uORF:ref_type=uORF".to_string(),
        "engine\t1:4000-4001\tTG\tENST4\t5UTR_annotation\talt_type=uORF:ref_type=uORF,alt_type=inFrame_oORF:ref_type=uORF".to_string(),
    ];
    // With the seeded runs the frameshift they write beside the gained start is the
    // precedence class's, and the annotation of the row on which a seeded run reports the
    // gained start is this class's; without them the first frameshift is this class's too
    // and that annotation stays a residual on both sides.
    let cases = [
        (
            "given",
            Some(input.as_str()),
            vec![seed_a.as_str(), seed_b.as_str()],
            json!({WINDOW: 0, PRECEDENCE: 2, LOST_AFTER: 0, REMADE: 6}),
            [21, 15, 9, 13, 9, 9],
            stays.clone(),
        ),
        (
            "without",
            None,
            vec![],
            json!({WINDOW: 0, PRECEDENCE: 0, LOST_AFTER: 0, REMADE: 7}),
            [21, 15, 9, 14, 10, 9],
            [swept, stays].concat(),
        ),
    ];
    for (name, given, sweep, masked, counted, residual) in cases {
        let scored = score(
            &scratch,
            &reference,
            &engine,
            &Run {
                name,
                table: &table,
                plugins: "UTRAnnotator",
                input: given,
                sweep: &sweep,
            },
        );
        assert_eq!(scored.masked("UTRAnnotator", "feature"), masked, "{name}");
        assert_eq!(scored.counts("UTRAnnotator"), counted, "{name}");
        assert_eq!(
            sorted(scored.residual_lines("UTRAnnotator", "feature")),
            sorted(residual),
            "{name}"
        );
    }
}

const THROUGH: &str =
    "alt_type=uORF:alt_type_length=150:ref_StartDistanceToCDS=120:ref_type=OutOfFrame_oORF:ref_type_length=NA";
const ENDS: &str =
    "alt_type=uORF:alt_type_length=45:ref_StartDistanceToCDS=87:ref_type=uORF:ref_type_length=36";

#[test]
fn frameshift_list_at_a_stop_codon_boundary_is_set_aside() {
    let scratch = Scratch::new("stop-boundary");
    // The class reads the plugin's two fields; the table of this run hands it the
    // annotation alone.
    assert_eq!(
        class_fields(BOUNDARY),
        Some(
            &[
                "5UTR_consequence",
                "5utr_consequence",
                "5UTR_annotation",
                "5utr_annotation"
            ][..]
        )
    );
    // On every row both report the frameshift and the engine lists one uORF more, whose
    // stop codon begins 54 bases before the coding start. The row's HGVSc places the base
    // directly after the insertion there on `boundary`, on `shifted_back` once its
    // `HGVS_OFFSET` of two bases is taken back, and on `swept`, where a seeded run reports
    // a start codon gained at the index the inserted base leaves it at; one base nearer
    // the coding start on `elsewhere`; and `unstated` carries no HGVSc.
    let rows = [
        (
            "boundary",
            "1:1000-1001",
            "A",
            "ENST1",
            "HGVSc=ENST1.1:c.-55_-54insA",
        ),
        (
            "shifted_back",
            "1:2000-2001",
            "A",
            "ENST2",
            "HGVSc=ENST2.1:c.-53_-52insA;HGVS_OFFSET=2",
        ),
        (
            "swept",
            "1:3000-3001",
            "G",
            "ENST3",
            "HGVSc=ENST3.1:c.-55_-54insG",
        ),
        (
            "elsewhere",
            "1:4000-4001",
            "A",
            "ENST4",
            "HGVSc=ENST4.1:c.-54_-53insA",
        ),
        ("unstated", "1:5000-5001", "A", "ENST5", ""),
    ];
    let line = |index: usize, annotation: &str, consequence: &str| -> String {
        let (name, location, allele, feature, notation) = rows[index];
        let mut fields = utr_fields(Some("1"), (consequence, annotation));
        if !notation.is_empty() {
            fields = format!("{fields};{notation}");
        }
        row(location, allele, feature, &fields, "146-147", name)
    };
    let longer = format!("{THROUGH},{ENDS}");
    let all = |annotation: &str| -> Vec<String> {
        (0..rows.len())
            .map(|index| line(index, annotation, FRAMESHIFT))
            .collect()
    };
    let reference = scratch.output("reference.txt", &all(THROUGH));
    let engine = scratch.output("engine.txt", &all(&longer));
    let records: Vec<[String; 5]> = rows
        .iter()
        .map(|(name, location, allele, _, _)| {
            let position = location
                .split_once(':')
                .unwrap()
                .1
                .split_once('-')
                .unwrap()
                .0;
            [
                "1".to_string(),
                position.to_string(),
                name.to_string(),
                "G".to_string(),
                format!("G{allele}"),
            ]
        })
        .collect();
    let records: Vec<[&str; 5]> = records
        .iter()
        .map(|record| record.each_ref().map(String::as_str))
        .collect();
    let input = scratch.vcf("input.vcf", &records);
    let seed_a = scratch.output("seed-a.txt", &[line(2, THROUGH, FRAMESHIFT)]);
    let seed_b = scratch.output(
        "seed-b.txt",
        &[line(2, "CapDistanceToStart=144:type=uORF", GAIN)],
    );
    let table = class_table(&[WINDOW, PRECEDENCE, LOST_AFTER, REMADE], &[BOUNDARY]);
    let stays: Vec<String> = [("1:4000-4001", "ENST4"), ("1:5000-5001", "ENST5")]
        .iter()
        .flat_map(|(location, feature)| {
            [
                format!("reference\t{location}\tA\t{feature}\t5UTR_annotation\t{THROUGH}"),
                format!("engine\t{location}\tA\t{feature}\t5UTR_annotation\t{THROUGH},{ENDS}"),
            ]
        })
        .collect();
    // With the seeded runs, the annotation of the row on which one of them reports the
    // start codon the inserted base leaves in place is the remade-in-place class's;
    // without them it is this class's, as the two rows whose HGVSc places the insertion at
    // the added uORF's stop codon are either way.
    let cases = [
        (
            "given",
            Some(input.as_str()),
            vec![seed_a.as_str(), seed_b.as_str()],
            json!({WINDOW: 0, PRECEDENCE: 0, LOST_AFTER: 0, REMADE: 1, BOUNDARY: 2}),
        ),
        (
            "without",
            None,
            vec![],
            json!({WINDOW: 0, PRECEDENCE: 0, LOST_AFTER: 0, REMADE: 0, BOUNDARY: 3}),
        ),
    ];
    for (name, given, sweep, masked) in cases {
        let scored = score(
            &scratch,
            &reference,
            &engine,
            &Run {
                name,
                table: &table,
                plugins: "UTRAnnotator",
                input: given,
                sweep: &sweep,
            },
        );
        assert_eq!(scored.masked("UTRAnnotator", "feature"), masked, "{name}");
        assert_eq!(
            sorted(scored.residual_lines("UTRAnnotator", "feature")),
            sorted(stays.clone()),
            "{name}"
        );
        // Of the fifteen keys each side writes, ten are equal, three are set aside on both
        // sides and two are left a side.
        assert_eq!(
            scored.counts("UTRAnnotator"),
            [15, 15, 10, 12, 12, 10],
            "{name}"
        );
    }
}

#[test]
fn lost_stop_written_before_the_stop_codon_is_set_aside() {
    let scratch = Scratch::new("before-stop");
    assert_eq!(
        class_fields(BEFORE_STOP),
        Some(
            &[
                "5UTR_consequence",
                "5utr_consequence",
                "5UTR_annotation",
                "5utr_annotation"
            ][..]
        )
    );
    // On every row, of two bases between cDNA positions 146 and 147, the reference reports
    // the frameshift and lists a uORF that starts 87 bases before the coding start and is
    // 45 bases long once shifted, and the engine reports the lost stop. The engine's group
    // names that uORF's next stop codon, 47 bases before the coding start, on `before`, on
    // `swept`, where the seeded runs write both terms on the row, and on `replaced`, whose
    // record replaces one base; it names a stop codon one base off on `elsewhere`.
    let shifted = format!("{THROUGH},KozakContext=GCCATGG:KozakStrength=Strong:{ENDS}");
    let lost = |distance: u32| {
        format!("AltStop=True:AltStopDistanceToCDS={distance}:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong")
    };
    let rows = [
        ("before", "1:1000-1001", "ENST1", 47, ["G", "GAT"]),
        ("swept", "1:2000-2001", "ENST2", 47, ["G", "GAT"]),
        ("elsewhere", "1:3000-3001", "ENST3", 48, ["G", "GAT"]),
        ("replaced", "1:4000-4001", "ENST4", 47, ["GC", "GAT"]),
    ];
    let line = |index: usize, engine: bool| -> String {
        let (name, location, feature, distance, _) = rows[index];
        let written = lost(distance);
        let effect = if engine {
            (STOP_LOSS, written.as_str())
        } else {
            (FRAMESHIFT, shifted.as_str())
        };
        row(
            location,
            "AT",
            feature,
            &utr_fields(Some("1"), effect),
            "146-147",
            name,
        )
    };
    let all = |engine: bool| -> Vec<String> {
        (0..rows.len()).map(|index| line(index, engine)).collect()
    };
    let reference = scratch.output("reference.txt", &all(false));
    let engine = scratch.output("engine.txt", &all(true));
    let records: Vec<[&str; 5]> = rows
        .iter()
        .map(|(name, location, _, _, alleles)| {
            let position = location
                .split_once(':')
                .unwrap()
                .1
                .split_once('-')
                .unwrap()
                .0;
            ["1", position, *name, alleles[0], alleles[1]]
        })
        .collect();
    let input = scratch.vcf("input.vcf", &records);
    let seed_a = scratch.output("seed-a.txt", &[line(1, false)]);
    let seed_b = scratch.output("seed-b.txt", &[line(1, true)]);
    let table = class_table(
        &[WINDOW, PRECEDENCE, LOST_AFTER, REMADE, BEFORE_STOP],
        &[BOUNDARY],
    );
    let stays = |location: &str, feature: &str, distance: u32| -> Vec<String> {
        vec![
            format!("reference\t{location}\tAT\t{feature}\t5UTR_annotation\t{shifted}"),
            format!("reference\t{location}\tAT\t{feature}\t5UTR_consequence\t{FRAMESHIFT}"),
            format!(
                "engine\t{location}\tAT\t{feature}\t5UTR_annotation\t{}",
                lost(distance)
            ),
            format!("engine\t{location}\tAT\t{feature}\t5UTR_consequence\t{STOP_LOSS}"),
        ]
    };
    let elsewhere = stays("1:3000-3001", "ENST3", 48);
    let replaced = stays("1:4000-4001", "ENST4", 47);
    // With the run's input and the seeded runs, the row on which those runs write both
    // terms is the precedence class's, the row whose record states a reference allele that
    // is not empty stays a residual, and the row whose group names the listed uORF's next
    // stop codon is this class's with both fields. Without them the three are this
    // class's.
    let cases = [
        (
            "given",
            Some(input.as_str()),
            vec![seed_a.as_str(), seed_b.as_str()],
            json!({WINDOW: 0, PRECEDENCE: 2, LOST_AFTER: 0, REMADE: 0, BOUNDARY: 0, BEFORE_STOP: 2}),
            [12, 12, 4, 8, 8, 4],
            [elsewhere.clone(), replaced].concat(),
        ),
        (
            "without",
            None,
            vec![],
            json!({WINDOW: 0, PRECEDENCE: 0, LOST_AFTER: 0, REMADE: 0, BOUNDARY: 0, BEFORE_STOP: 6}),
            [12, 12, 4, 6, 6, 4],
            elsewhere,
        ),
    ];
    for (name, given, sweep, masked, counted, residual) in cases {
        let scored = score(
            &scratch,
            &reference,
            &engine,
            &Run {
                name,
                table: &table,
                plugins: "UTRAnnotator",
                input: given,
                sweep: &sweep,
            },
        );
        assert_eq!(scored.masked("UTRAnnotator", "feature"), masked, "{name}");
        assert_eq!(scored.counts("UTRAnnotator"), counted, "{name}");
        assert_eq!(
            sorted(scored.residual_lines("UTRAnnotator", "feature")),
            sorted(residual),
            "{name}"
        );
    }
}

const KEPT: &str =
    "alt_type=uORF:alt_type_length=36:ref_StartDistanceToCDS=87:ref_type=uORF:ref_type_length=36";
const LOST_GROUP: &str =
    "AltStop=True:AltStopDistanceToCDS=47:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong";
const TWO_FIELDS: [&str; 2] = ["5UTR_consequence", "5UTR_annotation"];

/// One question about a row: the class, the reference's entry and the engine's, and what
/// the run hands the class.
struct Asked<'a> {
    klass: &'a str,
    e: &'a Row,
    a: &'a Row,
    sweep: Option<&'a UtrannotatorSeedSweep>,
    record_rows: Option<&'a [Row]>,
}

impl Asked<'_> {
    /// The class's answer for the field `name`, each side's value as its entry writes it.
    fn field(&self, name: &str) -> bool {
        let written = |entry: &Row| entry.fields.get(name).unwrap_or("-").to_string();
        BuiltIn
            .accepts(&Consult {
                klass: self.klass,
                e: self.e,
                name,
                expected: &written(self.e),
                actual: &written(self.a),
                a: Some(self.a),
                record_rows: self.record_rows,
                candidate_rows: None,
                declared: &BTreeSet::new(),
                assembly: None,
                source: None,
                gwas_pairs: None,
                revel_table: None,
                utrannotator_sweep: self.sweep,
                downstream_sequences: None,
                splice_windows: None,
                spliceai_table: None,
                exon_spans: None,
            })
            .unwrap()
    }

    /// The answers for `5UTR_consequence` and for `5UTR_annotation`.
    fn fields(&self) -> [bool; 2] {
        TWO_FIELDS.map(|name| self.field(name))
    }
}

fn asked<'a>(klass: &'a str, e: &'a Row, a: &'a Row) -> Asked<'a> {
    Asked {
        klass,
        e,
        a,
        sweep: None,
        record_rows: None,
    }
}

/// The row of an insertion between cDNA positions 146 and 147 of a forward-strand
/// transcript: its Allele, the two fields, the plugin's count and `more` after them.
fn inserted_line(allele: &str, effect: (&str, &str), more: &str) -> String {
    let mut fields = utr_fields(Some("1"), effect);
    if !more.is_empty() {
        fields = format!("{fields};{more}");
    }
    row("1:1000-1001", allele, "ENST1", &fields, "146-147", "ins")
}

fn inserted(allele: &str, effect: (&str, &str), more: &str) -> Row {
    entries(&[inserted_line(allele, effect, more)]).remove(0)
}

/// `inserted` without the plugin's count.
fn uncounted(allele: &str, effect: (&str, &str)) -> Row {
    let line = inserted_line(allele, effect, "");
    assert!(line.contains(";Existing_uORFs=1"));
    entries(&[line.replace(";Existing_uORFs=1", "")]).remove(0)
}

fn coding_notation(body: &str) -> String {
    format!("HGVSc=ENST1.1:c.{body}")
}

#[test]
fn frameshift_of_an_orf_that_keeps_its_length_is_set_aside() {
    // The reference lists a uORF whose stop codon begins 54 bases before the coding start
    // and whose two lengths are equal; the engine leaves it out, reporting nothing or the
    // frameshift of the other ORF alone.
    let at_stop = coding_notation("-55_-54insAC");
    let listed = |groups: &str, more: &str| inserted("AC", (FRAMESHIFT, groups), more);
    let nothing = inserted("AC", NO_EFFECT, "");
    let fewer = inserted("AC", (FRAMESHIFT, THROUGH), "");
    let two = format!("{THROUGH},{KEPT}");
    assert_eq!(
        asked(BOUNDARY, &listed(KEPT, &at_stop), &nothing).fields(),
        [true, true]
    );
    // Where both report the frameshift the annotation alone is the class's.
    assert_eq!(
        asked(BOUNDARY, &listed(&two, &at_stop), &fewer).fields(),
        [false, true]
    );
    // The stop codon begins at or after the base after the insertion as the record writes
    // it and no more than one base after the one that follows its last writing.
    for (body, accepted) in [
        ("-55_-54insAC", true),
        ("-56_-55insAC", true),
        ("-57_-56insAC", false),
        ("-54_-53insAC", false),
        ("-53_-52insAC;HGVS_OFFSET=3", true),
        ("-53_-52insAC;HGVS_OFFSET=-3", true),
        ("-50_-49insAC;HGVS_OFFSET=2", false),
        ("-58_-57insAC;HGVS_OFFSET=1", false),
        ("-56dup", true),
        ("-57dup", false),
    ] {
        let reference = listed(KEPT, &coding_notation(body));
        assert_eq!(
            asked(BOUNDARY, &reference, &nothing).fields(),
            [accepted; 2],
            "{body}"
        );
    }
    // A group the engine drops whose two lengths differ, that closes outside the 5' UTR
    // once shifted or that is no uORF of the reference UTR.
    for dropped in [
        ENDS.to_string(),
        KEPT.replace("alt_type=uORF", "alt_type=OutOfFrame_oORF"),
        KEPT.replace("ref_type=uORF", "ref_type=InFrame_oORF"),
        KEPT.replace("alt_type_length=36:", ""),
        KEPT.replace("ref_type_length=36", "ref_type_length=NA"),
    ] {
        let reference = listed(&dropped, &at_stop);
        assert_eq!(
            asked(BOUNDARY, &reference, &nothing).fields(),
            [false; 2],
            "{dropped}"
        );
    }
    // One such group beside one that is not, a list on each side that holds a group the
    // other lacks, an engine's entry without the plugin's counts or with another term, a
    // row that states no insertion point, and two entries that report nothing.
    let refused = [
        (listed(&format!("{KEPT},{ENDS}"), &at_stop), nothing.clone()),
        (
            listed(&two, &at_stop),
            inserted("AC", (FRAMESHIFT, ENDS), ""),
        ),
        (listed(KEPT, &at_stop), uncounted("AC", NO_EFFECT)),
        (
            listed(KEPT, &at_stop),
            inserted("AC", (LOSS, &at_index(145)), ""),
        ),
        (listed(KEPT, ""), nothing.clone()),
        (inserted("AC", NO_EFFECT, &at_stop), nothing.clone()),
    ];
    for (number, (reference, engine)) in refused.iter().enumerate() {
        assert_eq!(
            asked(BOUNDARY, reference, engine).fields(),
            [false; 2],
            "{number}"
        );
    }
}

#[test]
fn frameshift_at_a_stop_codon_boundary_beside_an_empty_reference_is_set_aside() {
    // The reference reports nothing and the engine the frameshift of a uORF whose stop
    // codon begins 54 bases before the coding start, directly after the insertion, and
    // whose two lengths differ.
    let empty = |body: &str| inserted("AC", NO_EFFECT, &coding_notation(body));
    let shifted = |groups: &str| inserted("AC", (FRAMESHIFT, groups), "");
    let at_stop = empty("-55_-54insAC");
    assert_eq!(
        asked(BOUNDARY, &at_stop, &shifted(ENDS)).fields(),
        [true, true]
    );
    assert_eq!(
        asked(
            BOUNDARY,
            &empty("-53_-52insAC;HGVS_OFFSET=2"),
            &shifted(ENDS)
        )
        .fields(),
        [true, true]
    );
    // A stop codon one base off the insertion either way, a group whose two lengths are
    // equal, a group that is no uORF beside one that fits, and a reference that reports
    // another term or fills the annotation alone.
    let refused = [
        (empty("-56_-55insAC"), shifted(ENDS)),
        (empty("-54_-53insAC"), shifted(ENDS)),
        (empty("-55_-54insAC"), shifted(KEPT)),
        (empty("-55_-54insAC"), shifted(&format!("{ENDS},{THROUGH}"))),
        (
            inserted(
                "AC",
                (LOSS, &at_index(145)),
                &coding_notation("-55_-54insAC"),
            ),
            shifted(ENDS),
        ),
        (
            inserted("AC", ("-", THROUGH), &coding_notation("-55_-54insAC")),
            shifted(ENDS),
        ),
    ];
    for (number, (reference, engine)) in refused.iter().enumerate() {
        assert_eq!(
            asked(BOUNDARY, reference, engine).fields(),
            [false; 2],
            "{number}"
        );
    }
}

#[test]
fn lost_stop_of_a_record_written_at_or_before_the_stop_codon_is_set_aside() {
    let scratch = Scratch::new("unread-stop");
    // Two bases the record writes one base before the last index that gives the same
    // sequence. The engine reports a lost stop; the reference nothing, the lost start, the
    // gained start or the frameshift of another ORF.
    let further = coding_notation("-54_-53insAC;HGVS_OFFSET=1");
    let engine = inserted("AC", (STOP_LOSS, LOST_GROUP), "");
    let other_orf = format!("KozakContext=AAAATGC:KozakStrength=Weak:{ENDS}");
    let same_orf = format!("KozakContext=GCCATGG:KozakStrength=Strong:{ENDS}");
    let (lost_start, gained_start) = (at_index(145), at_index(140));
    let references: [(&str, &str); 4] = [
        NO_EFFECT,
        (LOSS, &lost_start),
        (GAIN, &gained_start),
        (FRAMESHIFT, &other_orf),
    ];
    for effect in references {
        let reference = inserted("AC", effect, &further);
        assert_eq!(
            asked(BEFORE_STOP, &reference, &engine).fields(),
            [true, true],
            "{}",
            effect.0
        );
    }
    let empty = inserted("AC", NO_EFFECT, &further);
    // A frameshift list that holds the engine's ORF is read as the lost stop of a listed
    // ORF: its group names the stop codon the listed length ends on, 47 bases before the
    // coding start, and no other.
    let listing = inserted("AC", (FRAMESHIFT, &same_orf), &further);
    assert_eq!(asked(BEFORE_STOP, &listing, &engine).fields(), [true, true]);
    let elsewhere = inserted(
        "AC",
        (
            STOP_LOSS,
            &LOST_GROUP.replace("AltStopDistanceToCDS=47", "AltStopDistanceToCDS=48"),
        ),
        "",
    );
    assert_eq!(
        asked(BEFORE_STOP, &listing, &elsewhere).fields(),
        [false; 2]
    );
    assert_eq!(
        asked(BEFORE_STOP, &empty, &elsewhere).fields(),
        [true, true]
    );
    // An allele whose length is a multiple of three, a row without an offset or with one
    // of zero, a reference that reports the gained stop or its own lost stop or fills the
    // annotation alone, and an engine's entry without the plugin's counts, whose group
    // lacks a Kozak key or whose annotation lists no group.
    let in_frame = [
        inserted("ACG", NO_EFFECT, &further),
        inserted("ACGTCA", NO_EFFECT, &further),
    ];
    for reference in &in_frame {
        let engine = inserted(&reference.allele, (STOP_LOSS, LOST_GROUP), "");
        assert_eq!(
            asked(BEFORE_STOP, reference, &engine).fields(),
            [false; 2],
            "{}",
            reference.allele
        );
    }
    let refused = [
        (
            inserted("AC", NO_EFFECT, &coding_notation("-54_-53insAC")),
            engine.clone(),
        ),
        (
            inserted(
                "AC",
                NO_EFFECT,
                &coding_notation("-54_-53insAC;HGVS_OFFSET=0"),
            ),
            engine.clone(),
        ),
        (
            inserted("AC", (STOP_GAIN, "ref_type=uORF"), &further),
            engine.clone(),
        ),
        (
            inserted("AC", (STOP_LOSS, LOST_GROUP), &further),
            elsewhere.clone(),
        ),
        (inserted("AC", ("-", THROUGH), &further), engine.clone()),
        (empty.clone(), uncounted("AC", (STOP_LOSS, LOST_GROUP))),
        (
            empty.clone(),
            inserted(
                "AC",
                (STOP_LOSS, &LOST_GROUP.replace(":KozakStrength=Strong", "")),
                "",
            ),
        ),
        (empty.clone(), inserted("AC", (STOP_LOSS, "{}"), "")),
    ];
    for (number, (reference, engine)) in refused.iter().enumerate() {
        assert_eq!(
            asked(BEFORE_STOP, reference, engine).fields(),
            [false; 2],
            "{number}"
        );
    }
    // Seeded runs that write no lost stop on the row leave it the class's; one run that
    // writes it there takes the row from the class.
    let reference_line = inserted_line("AC", NO_EFFECT, &further);
    let engine_line = inserted_line("AC", (STOP_LOSS, LOST_GROUP), "");
    let quiet = seed_sweep(
        &scratch,
        "quiet",
        &[&default_text(std::slice::from_ref(&reference_line))],
    );
    let writing = seed_sweep(
        &scratch,
        "writing",
        &[
            &default_text(&[reference_line]),
            &default_text(&[engine_line]),
        ],
    );
    let with = |sweep| Asked {
        sweep: Some(sweep),
        ..asked(BEFORE_STOP, &empty, &engine)
    };
    assert_eq!(with(&quiet).fields(), [true, true]);
    assert_eq!(with(&writing).fields(), [false; 2]);
}

#[test]
fn lost_stop_of_a_stop_codon_that_stands_whole_in_frame_is_set_aside() {
    // Six bases whose record writes the base after them 54 bases before the coding start.
    // The reference reports a lost stop and names the next stop codon; the engine reports
    // no lost stop.
    let standing = |distance: &str| {
        format!("AltStop=True:AltStopDistanceToCDS={distance}:FrameWithCDS=inFrame:KozakContext=GCCATGG:KozakStrength=Strong")
    };
    let lost = |allele: &str, groups: &str, body: &str| {
        inserted(allele, (STOP_LOSS, groups), &coding_notation(body))
    };
    let nothing = inserted("ACGTCA", NO_EFFECT, "");
    // The ORF's own stop codon, moved on by the six bases, lies one or two bases before
    // the base after the insertion, and a stop codon the inserted bases hold in frame
    // before it three bases further from the coding start.
    for (distance, accepted) in [
        ("53", false),
        ("54", false),
        ("55", true),
        ("56", true),
        ("57", false),
        ("58", true),
        ("59", true),
        ("60", false),
        ("61", false),
        ("NA", false),
    ] {
        let reference = lost("ACGTCA", &standing(distance), "-55_-54insACGTCA");
        assert_eq!(
            asked(BEFORE_STOP, &reference, &nothing).fields(),
            [accepted; 2],
            "{distance}"
        );
    }
    let reference = lost("ACGTCA", &standing("56"), "-55_-54insACGTCA");
    assert_eq!(
        asked(
            BEFORE_STOP,
            &lost("ACGTCA", &standing("56"), "-54_-53insACGTCA;HGVS_OFFSET=1"),
            &nothing
        )
        .fields(),
        [true, true]
    );
    // The engine may report another of the plugin's terms.
    let gaining = inserted("ACGTCA", (GAIN, &at_index(147)), "");
    assert_eq!(
        asked(BEFORE_STOP, &reference, &gaining).fields(),
        [true, true]
    );
    // An allele whose length is no multiple of three, with the stop codon anywhere.
    for allele in ["ACGTC", "ACGT", "AC"] {
        let body = format!("-55_-54ins{allele}");
        let engine = inserted(allele, NO_EFFECT, "");
        for distance in ["55", "56"] {
            assert_eq!(
                asked(
                    BEFORE_STOP,
                    &lost(allele, &standing(distance), &body),
                    &engine
                )
                .fields(),
                [false; 2],
                "{allele} {distance}"
            );
        }
    }
    // One group that fits beside one that does not, a group without an alternative stop
    // codon, an annotation that lists no group, a row that states no insertion point, and
    // an engine's entry that reports the lost stop, a term outside the plugin's five, an
    // annotation without a term or no counts.
    let two = format!("{},{}", standing("56"), standing("57"));
    let no_stop = standing("56").replace("AltStop=True", "AltStop=False");
    let refused = [
        (lost("ACGTCA", &two, "-55_-54insACGTCA"), nothing.clone()),
        (
            lost("ACGTCA", &no_stop, "-55_-54insACGTCA"),
            nothing.clone(),
        ),
        (lost("ACGTCA", "{}", "-55_-54insACGTCA"), nothing.clone()),
        (
            inserted("ACGTCA", (STOP_LOSS, &standing("56")), ""),
            nothing.clone(),
        ),
        (
            reference.clone(),
            inserted("ACGTCA", (STOP_LOSS, &standing("59")), ""),
        ),
        (
            reference.clone(),
            inserted("ACGTCA", ("5_prime_UTR_variant", "type=uORF"), ""),
        ),
        (
            reference.clone(),
            inserted("ACGTCA", ("-", "type=uORF"), ""),
        ),
        (reference.clone(), uncounted("ACGTCA", NO_EFFECT)),
    ];
    for (number, (reference, engine)) in refused.iter().enumerate() {
        assert_eq!(
            asked(BEFORE_STOP, reference, engine).fields(),
            [false; 2],
            "{number}"
        );
    }
}

#[test]
fn gained_start_named_before_the_insertion_is_set_aside_where_the_allele_spells_it() {
    // An insertion whose record spells its alleles, the base after it at index 300 of the
    // 5' UTR. The engine alone reports a start codon gained at `at`.
    let line = |allele: &str, strand: Option<&str>, effect: (&str, &str)| {
        row(
            "1:4000-4001",
            allele,
            "ENST4",
            &utr_fields(strand, effect),
            "300-301",
            &format!("1_4001_-/{allele}"),
        )
    };
    let pair = |allele: &str, strand: Option<&str>, at: u32, reference: (&str, &str)| {
        let rows = entries(&[
            line(allele, strand, reference),
            line(allele, strand, (GAIN, &at_index(at))),
        ]);
        asked(WINDOW, &rows[0], &rows[1]).fields()
    };
    let forward = Some("1");
    // Three or more bases before the insertion, the codon the allele spells round its end
    // from as many bases before its end; on a reverse-strand transcript the allele is read
    // reverse-complemented.
    for (allele, strand, at, accepted) in [
        ("CATG", forward, 297, true),
        ("TGCA", forward, 295, true),
        ("CATG", forward, 293, true),
        ("ATG", forward, 291, true),
        ("CATC", Some("-1"), 297, true),
        ("CATC", forward, 297, false),
        ("CATG", Some("-1"), 297, true),
        ("CATG", forward, 296, false),
        ("CATG", forward, 295, false),
        ("TGCA", forward, 297, false),
        // Fewer than three bases before the insertion, though the allele spells the codon.
        ("GCAT", forward, 298, false),
        ("TGAC", forward, 299, false),
        ("CATG", None, 297, false),
    ] {
        assert_eq!(
            pair(allele, strand, at, NO_EFFECT),
            [accepted; 2],
            "{allele} {strand:?} {at}"
        );
    }
    // Inside the allele's length the codon is the published window's.
    assert_eq!(pair("CATG", forward, 301, NO_EFFECT), [true, true]);
    assert_eq!(pair("CATG", forward, 304, NO_EFFECT), [false; 2]);
    // The reference may report any other of the plugin's terms; a reference that reports
    // a gained start of its own is not the class's.
    let lost_start = at_index(299);
    let references: [(&str, &str); 3] = [
        (LOSS, &lost_start),
        (FRAMESHIFT, SHIFT),
        (STOP_LOSS, LOST_GROUP),
    ];
    for reference in references {
        assert_eq!(
            pair("CATG", forward, 297, reference),
            [true, true],
            "{}",
            reference.0
        );
    }
    assert_eq!(
        pair("CATG", forward, 297, (GAIN, &at_index(301))),
        [false; 2]
    );
    // The reference alone reporting the codon before the insertion is not the class's.
    let swapped = entries(&[
        line("CATG", forward, (GAIN, &at_index(297))),
        line("CATG", forward, NO_EFFECT),
    ]);
    assert_eq!(asked(WINDOW, &swapped[0], &swapped[1]).fields(), [false; 2]);
}

/// A structural record's row under the identifier `name`.
fn named_span_row(
    name: &str,
    location: &str,
    allele: &str,
    cdna: &str,
    fields: &[&str],
    terms: &str,
) -> String {
    let text = utr_span_row(location, allele, cdna, fields, terms);
    assert!(text.starts_with("sv\t"));
    text.replacen("sv", name, 1)
}

/// The rows of one record as a run hands them with the record's input line: the rows and
/// an entry of the line's CHROM, POS, REF and ALT.
fn handed(rows: &[Row], line: &str) -> Vec<Row> {
    RecordRows::new(rows, &[line])
        .of(&rows[0])
        .unwrap()
        .unwrap()
}

#[test]
fn effect_of_an_inversion_inside_one_utr_exon_is_set_aside() {
    let stop_gain = format!("5UTR_consequence={STOP_GAIN};5UTR_annotation=ref_type=uORF");
    let start_gain = format!("5UTR_consequence={GAIN};5UTR_annotation=type=uORF");
    let frameshift = format!("5UTR_consequence={FRAMESHIFT};5UTR_annotation=ref_type=uORF");
    let term_alone = format!("5UTR_consequence={STOP_GAIN}");
    let utr = "5_prime_UTR_variant";
    let nmd = "5_prime_UTR_variant,NMD_transcript_variant";
    let pair = |shape: [&str; 4], reference: &[&str], engine: &[&str]| {
        let [location, allele, cdna, terms] = shape;
        let rows = entries(&[
            utr_span_row(location, allele, cdna, reference, terms),
            utr_span_row(location, allele, cdna, engine, terms),
        ]);
        asked(STRUCTURAL, &rows[0], &rows[1]).fields()
    };
    let three = ["21:1001-1003", "inversion", "101-103", utr];
    // The engine reports a term a replacement of equal length can fire, beside the
    // reference's reading of the ALT text as deleted bases or beside nothing.
    for engine in [
        stop_gain.as_str(),
        start_gain.as_str(),
        STRUCTURAL_LOST_START,
        SPAN_EFFECT,
    ] {
        assert_eq!(pair(three, &[SPLICED], &[engine]), [true, true], "{engine}");
        assert_eq!(pair(three, &[], &[engine]), [true, true], "{engine}");
    }
    assert_eq!(
        pair(
            ["21:1001-1003", "inversion", "101-103", nmd],
            &[SPLICED],
            &[&stop_gain]
        ),
        [true, true]
    );
    // The frameshift, whatever the span; a span that is not the distance between the cDNA
    // positions or has an end outside an exon; a row that leaves the 5' UTR or names a
    // deletion's term; a term without an annotation; and another class word.
    let refused: [([&str; 4], &str); 9] = [
        (["21:1001-1004", "inversion", "101-104", utr], &frameshift),
        (three, &frameshift),
        (["21:1001-1040", "inversion", "101-103", utr], &stop_gain),
        (["21:1001-1003", "inversion", "?-103", utr], &stop_gain),
        (
            [
                "21:1001-1003",
                "inversion",
                "101-103",
                "5_prime_UTR_variant,intron_variant",
            ],
            &stop_gain,
        ),
        (
            [
                "21:1001-1003",
                "inversion",
                "101-103",
                "feature_truncation,5_prime_UTR_variant",
            ],
            &stop_gain,
        ),
        (
            ["21:1001-1003", "inversion", "101-103", "intron_variant"],
            &stop_gain,
        ),
        (three, &term_alone),
        (["21:1001-1003", "duplication", "101-103", utr], &stop_gain),
    ];
    for (number, (shape, engine)) in refused.iter().enumerate() {
        assert_eq!(pair(*shape, &[SPLICED], &[engine]), [false; 2], "{number}");
    }
    // Where the run hands the record with the row, `<INV>` is its one alternate.
    let rows = entries(&[
        named_span_row(
            "inv",
            "21:1001-1003",
            "inversion",
            "101-103",
            &[SPLICED],
            utr,
        ),
        named_span_row(
            "inv",
            "21:1001-1003",
            "inversion",
            "101-103",
            &[&stop_gain],
            utr,
        ),
    ]);
    for (alternates, accepted) in [("<INV>", true), ("<DEL>,<INV>", false), ("<DUP>", false)] {
        let record = handed(&rows[..1], &format!("21\t1000\tinv\tN\t{alternates}"));
        let with_record = Asked {
            record_rows: Some(&record),
            ..asked(STRUCTURAL, &rows[0], &rows[1])
        };
        assert_eq!(with_record.fields(), [accepted; 2], "{alternates}");
    }
}

#[test]
fn effect_of_a_copy_number_loss_inside_one_utr_exon_is_set_aside() {
    let frameshift = format!("{UTR_FRAMESHIFT};{REFERENCE_KEYS}");
    let stop_gain = format!("5UTR_consequence={STOP_GAIN};5UTR_annotation=ref_type=uORF");
    let utr = "5_prime_UTR_variant";
    let nmd = "5_prime_UTR_variant,NMD_transcript_variant";
    let word = "copy_number_variation";
    let pair = |shape: [&str; 3], reference: &[&str], engine: &[&str]| {
        let [location, cdna, terms] = shape;
        let rows = entries(&[
            utr_span_row(location, word, cdna, reference, terms),
            utr_span_row(location, word, cdna, engine, terms),
        ]);
        asked(STRUCTURAL, &rows[0], &rows[1]).fields()
    };
    let four = ["21:2001-2004", "101-104", utr];
    // The engine reports the effect of the span's deletion beside the reference's reading
    // of the first ALT text, with the reference's own term or another, or beside nothing.
    assert_eq!(pair(four, &[SPLICED], &[&frameshift]), [true, true]);
    for engine in [frameshift.as_str(), STRUCTURAL_LOST_START, SPAN_EFFECT] {
        assert_eq!(pair(four, &[], &[engine]), [true, true], "{engine}");
        assert_eq!(
            pair(four, &[STRUCTURAL_LOST_START], &[engine]),
            [true, true],
            "{engine}"
        );
    }
    assert_eq!(
        pair(["21:2001-2004", "101-104", nmd], &[SPLICED], &[SPAN_EFFECT]),
        [true, true]
    );
    // A span that crosses exons, so that the cDNA positions lie nearer than its ends; the
    // frameshift on a span that is a multiple of three; the gained stop; a row that leaves
    // the 5' UTR or names a deletion's term; and an end outside an exon.
    let refused: [([&str; 3], &str); 6] = [
        (["21:2001-2040", "101-120", utr], SPAN_EFFECT),
        (["21:2001-2003", "101-103", utr], &frameshift),
        (four, &stop_gain),
        (
            [
                "21:2001-2004",
                "101-104",
                "5_prime_UTR_variant,intron_variant",
            ],
            SPAN_EFFECT,
        ),
        (
            [
                "21:2001-2004",
                "101-104",
                "feature_truncation,5_prime_UTR_variant",
            ],
            SPAN_EFFECT,
        ),
        (["21:2001-2004", "101-?", utr], SPAN_EFFECT),
    ];
    for (number, (shape, engine)) in refused.iter().enumerate() {
        assert_eq!(pair(*shape, &[SPLICED], &[engine]), [false; 2], "{number}");
    }
    // Where the run hands the record with its rows, one for each alternate in the
    // alternates' order, the row is that of an alternate whose own text is a deletion's
    // and the record states another beside it.
    let reference = entries(&[
        named_span_row("cnv", "21:2001-2004", word, "101-104", &[SPLICED], utr),
        named_span_row("cnv", "21:2001-2004", word, "101-104", &[SPLICED], utr),
    ]);
    let engine = &entries(&[named_span_row(
        "cnv",
        "21:2001-2004",
        word,
        "101-104",
        &[SPAN_EFFECT],
        utr,
    )])[0];
    for (alternates, rows, accepted) in [
        ("<CN0>,<CN2>", 2, [true, false]),
        ("<CN2>,<CN0>", 2, [false, true]),
        ("<DUP>,<DEL>", 2, [false, true]),
        ("<CN=0>,<CN=2>", 2, [true, false]),
        ("<CN2>,<CN3>", 2, [false, false]),
        ("<CN0>,<CN2>,<CN3>", 2, [false, false]),
        ("<CN0>", 1, [false, false]),
        ("<CNV>", 1, [false, false]),
    ] {
        let record = handed(
            &reference[..rows],
            &format!("21\t2000\tcnv\tN\t{alternates}"),
        );
        for (place, accepted) in accepted.iter().enumerate().take(rows) {
            let with_record = Asked {
                record_rows: Some(&record),
                ..asked(STRUCTURAL, &reference[place], engine)
            };
            assert_eq!(with_record.fields(), [*accepted; 2], "{alternates} {place}");
        }
    }
}

/// An insertion on a transcript of its own: its record, Allele, cDNA positions and HGVSc
/// body with the offset after it, then the reference's two fields and the engine's.
type InsertedShape<'a> = (
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    (&'a str, &'a str),
    (&'a str, &'a str),
);

#[test]
fn further_shapes_of_an_insertion_row_are_set_aside_under_the_classes_in_table_order() {
    let scratch = Scratch::new("further-shapes");
    let standing = |distance: u32| {
        format!("AltStop=True:AltStopDistanceToCDS={distance}:FrameWithCDS=inFrame:KozakContext=GCCATGG:KozakStrength=Strong")
    };
    let two = format!("{THROUGH},{KEPT}");
    let (in_frame, off_frame, earlier) = (standing(56), standing(57), at_index(297));
    let at_stop = "c.-55_-54ins";
    let before_stop = "c.-54_-53ins";
    // Each row is an insertion on a transcript of its own, the base after it 54 bases
    // before the coding start. The first six differ in a shape of the three classes;
    // `shifted` drops a group whose two lengths differ, `swept` reports a lost stop that a
    // seeded run writes on the row, `off` names a stop codon three bases before the base
    // after the insertion, and on `unspelt` the allele spells no start codon three bases
    // before its end.
    let shapes: [InsertedShape; 10] = [
        (
            "kept",
            "AC",
            "146-147",
            at_stop,
            "",
            (FRAMESHIFT, KEPT),
            NO_EFFECT,
        ),
        (
            "fewer",
            "AC",
            "146-147",
            at_stop,
            "",
            (FRAMESHIFT, &two),
            (FRAMESHIFT, THROUGH),
        ),
        (
            "added",
            "AC",
            "146-147",
            at_stop,
            "",
            NO_EFFECT,
            (FRAMESHIFT, ENDS),
        ),
        (
            "unread",
            "AC",
            "146-147",
            before_stop,
            ";HGVS_OFFSET=1",
            NO_EFFECT,
            (STOP_LOSS, LOST_GROUP),
        ),
        (
            "standing",
            "ACGTCA",
            "146-147",
            at_stop,
            "",
            (STOP_LOSS, &in_frame),
            NO_EFFECT,
        ),
        (
            "earlier",
            "CATG",
            "300-301",
            at_stop,
            "",
            NO_EFFECT,
            (GAIN, &earlier),
        ),
        (
            "shifted",
            "AC",
            "146-147",
            at_stop,
            "",
            (FRAMESHIFT, ENDS),
            NO_EFFECT,
        ),
        (
            "swept",
            "AC",
            "146-147",
            before_stop,
            ";HGVS_OFFSET=1",
            NO_EFFECT,
            (STOP_LOSS, LOST_GROUP),
        ),
        (
            "off",
            "ACGTCA",
            "146-147",
            at_stop,
            "",
            (STOP_LOSS, &off_frame),
            NO_EFFECT,
        ),
        (
            "unspelt",
            "CATC",
            "300-301",
            at_stop,
            "",
            NO_EFFECT,
            (GAIN, &earlier),
        ),
    ];
    let location = |number: usize| format!("1:{0}000-{0}001", number + 1);
    let feature = |number: usize| format!("ENST{}", number + 1);
    let line = |number: usize, engine: bool| -> String {
        let (name, allele, cdna, body, offset, reference, written) = shapes[number];
        let effect = if engine { written } else { reference };
        let fields = format!(
            "{};HGVSc={}.1:{body}{allele}{offset}",
            utr_fields(Some("1"), effect),
            feature(number)
        );
        row(
            &location(number),
            allele,
            &feature(number),
            &fields,
            cdna,
            name,
        )
    };
    let all = |engine: bool| -> Vec<String> {
        (0..shapes.len())
            .map(|number| line(number, engine))
            .collect()
    };
    let reference = scratch.output("reference.txt", &all(false));
    let engine = scratch.output("engine.txt", &all(true));
    let records: Vec<[String; 5]> = shapes
        .iter()
        .enumerate()
        .map(|(number, shape)| {
            [
                "1".to_string(),
                format!("{}000", number + 1),
                shape.0.to_string(),
                "G".to_string(),
                format!("G{}", shape.1),
            ]
        })
        .collect();
    let records: Vec<[&str; 5]> = records
        .iter()
        .map(|record| record.each_ref().map(String::as_str))
        .collect();
    let input = scratch.vcf("input.vcf", &records);
    let seed_a = scratch.output("seed-a.txt", &all(false));
    let seed_b = scratch.output("seed-b.txt", &[line(7, true)]);
    let table = class_table(
        &[
            WINDOW,
            PRECEDENCE,
            LOST_AFTER,
            REMADE,
            BOUNDARY,
            BEFORE_STOP,
        ],
        &[],
    );
    let left = |number: usize, engine: bool| -> Vec<String> {
        let (_, allele, _, _, _, reference, written) = shapes[number];
        let (side, effect) = if engine {
            ("engine", written)
        } else {
            ("reference", reference)
        };
        let key = format!(
            "{side}\t{}\t{allele}\t{}",
            location(number),
            feature(number)
        );
        vec![
            format!("{key}\t5UTR_annotation\t{}", effect.1),
            format!("{key}\t5UTR_consequence\t{}", effect.0),
        ]
    };
    // With the run's input and the seeded runs, the gained start is the window class's
    // from the record's alleles and the lost stop a seeded run writes stays a residual.
    // Without them the gained start stays, and the lost stop is the last class's.
    let cases = [
        (
            "given",
            Some(input.as_str()),
            vec![seed_a.as_str(), seed_b.as_str()],
            json!({WINDOW: 2, PRECEDENCE: 0, LOST_AFTER: 0, REMADE: 0, BOUNDARY: 5, BEFORE_STOP: 4}),
            [left(6, false), left(7, true), left(8, false), left(9, true)].concat(),
        ),
        (
            "without",
            None,
            vec![],
            json!({WINDOW: 0, PRECEDENCE: 0, LOST_AFTER: 0, REMADE: 0, BOUNDARY: 5, BEFORE_STOP: 6}),
            [left(5, true), left(6, false), left(8, false), left(9, true)].concat(),
        ),
    ];
    for (name, given, sweep, masked, residual) in cases {
        let scored = score(
            &scratch,
            &reference,
            &engine,
            &Run {
                name,
                table: &table,
                plugins: "UTRAnnotator",
                input: given,
                sweep: &sweep,
            },
        );
        assert_eq!(scored.masked("UTRAnnotator", "feature"), masked, "{name}");
        // Of the twenty values the reference writes and the engine's twenty-two, eleven
        // are equal; five of the reference's and seven of the engine's are set aside.
        assert_eq!(
            scored.counts("UTRAnnotator"),
            [20, 22, 11, 15, 15, 11],
            "{name}"
        );
        assert_eq!(
            sorted(scored.residual_lines("UTRAnnotator", "feature")),
            sorted(residual),
            "{name}"
        );
    }
}

#[test]
fn effects_on_inversion_and_copy_number_rows_are_set_aside_by_the_score_command() {
    let scratch = Scratch::new("structural-spans");
    let frameshift = format!("{UTR_FRAMESHIFT};{REFERENCE_KEYS}");
    let stop_gain = format!("5UTR_consequence={STOP_GAIN};5UTR_annotation=ref_type=uORF");
    let utr = "5_prime_UTR_variant";
    let word = "copy_number_variation";
    // A copy-number record of two alternates writes two rows under one key: the reference
    // its reading of the first ALT text on both, the engine the deletion's effect on the
    // loss alternate's row alone. `loss` lies inside one exon; the span of `crossing`
    // is longer than the distance between its cDNA positions. An inversion's row carries
    // the engine's gained stop on `turned` and a frameshift on `shifted`.
    let reference = scratch.output(
        "reference.txt",
        &[
            utr_span_row("21:2001-2004", word, "101-104", &[SPLICED], utr),
            utr_span_row("21:2001-2004", word, "101-104", &[SPLICED], utr),
            utr_span_row("21:3001-3040", word, "101-120", &[SPLICED], utr),
            utr_span_row("21:3001-3040", word, "101-120", &[SPLICED], utr),
            utr_span_row("21:4001-4003", "inversion", "101-103", &[SPLICED], utr),
            utr_span_row("21:5001-5004", "inversion", "101-104", &[SPLICED], utr),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            utr_span_row("21:2001-2004", word, "101-104", &[&frameshift], utr),
            utr_span_row("21:2001-2004", word, "101-104", &[], utr),
            utr_span_row("21:3001-3040", word, "101-120", &[SPAN_EFFECT], utr),
            utr_span_row("21:3001-3040", word, "101-120", &[], utr),
            utr_span_row("21:4001-4003", "inversion", "101-103", &[&stop_gain], utr),
            utr_span_row("21:5001-5004", "inversion", "101-104", &[&frameshift], utr),
        ],
    );
    let scored = score(
        &scratch,
        &reference,
        &engine,
        &Run {
            name: "run",
            table: &structural_table(),
            plugins: "UTRAnnotator",
            input: None,
            sweep: &[],
        },
    );
    let stays: [(&str, &str, &str, &[&str]); 4] = [
        ("reference", "21:3001-3040", word, &[SPLICED]),
        ("engine", "21:3001-3040", word, &[SPAN_EFFECT]),
        (
            "reference",
            "21:5001-5004",
            "inversion",
            &["5UTR_annotation=alt_type=uORF:ref_type=uORF"],
        ),
        ("engine", "21:5001-5004", "inversion", &[REFERENCE_KEYS]),
    ];
    let lines: Vec<String> = stays
        .iter()
        .flat_map(|(side, location, allele, fields)| {
            left(side, location, allele, "ENST1\t", fields)
        })
        .collect();
    assert_eq!(
        sorted(scored.residual_lines("UTRAnnotator", "feature")),
        sorted(lines)
    );
    // Set aside, a field once under a key: both fields of `loss`, where the reference's
    // second frameshift stands against the engine's empty row; both fields of `crossing`,
    // for the reference's second row alone, against the engine's empty row; and both
    // fields of `turned`. That is seven of the reference's values and three of the
    // engine's.
    assert_eq!(
        scored.masked("UTRAnnotator", "feature"),
        json!({STRUCTURAL: 6})
    );
    assert_eq!(scored.counts("UTRAnnotator"), [12, 8, 2, 5, 5, 2]);
}

/// Two planted classes, with what each consult was handed: the class and the address of
/// the seeded runs it was handed, if any. `planted_swept_terms` accepts
/// a pair of terms the seeded runs both write on the row; `planted_engine_zero` accepts an
/// engine's `0` beside no value of the reference.
#[derive(Default)]
struct Planted {
    handed: Mutex<Vec<(String, Option<usize>)>>,
}

impl Rules for Planted {
    fn knows(&self, klass: &str) -> bool {
        matches!(klass, "planted_swept_terms" | "planted_engine_zero")
    }

    fn accepts(&self, c: &Consult) -> Result<bool, Refusal> {
        let sweep = c.utrannotator_sweep;
        self.handed.lock().unwrap().push((
            c.klass.to_string(),
            sweep.map(|sweep| std::ptr::from_ref(sweep) as usize),
        ));
        Ok(match c.klass {
            "planted_swept_terms" => sweep
                .zip(c.a)
                .is_some_and(|(sweep, a)| sweep.writes(c.e, a, false)),
            _ => c.klass == "planted_engine_zero" && blank(c.expected) && c.actual == "0",
        })
    }

    fn transcript_reading(&self, _: &str) -> bool {
        false
    }

    fn sibling_allele(&self, _: &str) -> bool {
        false
    }

    fn record_reference(&self, _: &str) -> bool {
        false
    }
}

const SWEEP_TABLE: &str = "planted_swept_terms\tUTRAnnotator\t5UTR_consequence\tmet\t-\nplanted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n";

fn sha256(content: &[u8]) -> String {
    Sha256::digest(content)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// A reference and an engine that write different terms at two Locations, and a directory
/// of seeded runs that write both terms at the first Location alone, beside files that are
/// no output. Returns the reference, the engine and the directory.
fn sweep_inputs(scratch: &Scratch) -> (String, String, String) {
    let rows = |terms: [&str; 2]| -> Vec<String> {
        [100, 200]
            .iter()
            .zip(terms)
            .map(|(at, term)| {
                row(
                    &format!("1:{at}"),
                    "T",
                    "ENST1",
                    &format!("5UTR_consequence={term}"),
                    "-",
                    "v",
                )
            })
            .collect()
    };
    let reference = scratch.output(
        "reference.txt",
        &[
            rows([LOSS, LOSS]),
            vec![row("1:300", "T", "ENST1", "IMPACT=LOW", "-", "v")],
        ]
        .concat(),
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            rows([GAIN, GAIN]),
            vec![row("1:300", "T", "ENST1", "LOEUF=0", "-", "v")],
        ]
        .concat(),
    );
    fs::create_dir_all(scratch.0.join("sweep/tools")).unwrap();
    scratch.output("sweep/seed-a.txt", &rows([LOSS, LOSS]));
    let commented = format!(
        "## ENSEMBL VARIANT EFFECT PREDICTOR\n## Column descriptions:\n{HEADER}{}",
        rows([GAIN, STOP_LOSS]).concat()
    );
    scratch.file("sweep/seed-b.txt.gz", gzip(commented.as_bytes()));
    scratch.output("sweep/tools/seed-c.txt", &rows([GAIN, GAIN]));
    scratch.file("sweep/seed-a.command.txt", "vep --plugin UTRAnnotator\n");
    scratch.file("sweep/input.vcf", VCF_HEADER);
    scratch.file("sweep/seed-a.log", "");
    let bytes: Vec<u8> = (0..=255u8).cycle().take(256 * 400).collect();
    scratch.file("sweep/seeds.bin", bytes);
    (reference, engine, scratch.path("sweep"))
}

/// The exit status, the lines of standard output and standard error of one command line
/// under the planted classes.
fn run_planted(planted: &Planted, arguments: &[&str]) -> (i32, Vec<String>, String) {
    let arguments: Vec<String> = arguments.iter().map(|word| word.to_string()).collect();
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = run(&arguments, &BuiltInRegistry, planted, &mut out, &mut error);
    let out = String::from_utf8(out).unwrap();
    (
        code,
        out.lines().map(str::to_string).collect(),
        String::from_utf8(error).unwrap(),
    )
}

/// The seeded runs, given as a directory or as each output, are read and handed to the
/// classes of UTRAnnotator alone, the same runs at every consult; without them the class
/// that reads them sets nothing aside.
fn seeded_runs_are_handed_to_the_classes_of_utrannotator_alone(as_directory: bool) {
    let scratch = Scratch::new(if as_directory {
        "handed-directory"
    } else {
        "handed-outputs"
    });
    let (reference, engine, runs) = sweep_inputs(&scratch);
    let outputs = [
        format!("{runs}/seed-a.txt"),
        format!("{runs}/seed-b.txt.gz"),
    ];
    let table = scratch.file("classes.tsv", SWEEP_TABLE);
    let (first, second) = (format!("engine={engine}"), format!("again={engine}"));
    let scored = |name: &str, sweep: &[&str]| -> (Value, PathBuf, Planted) {
        let planted = Planted::default();
        let report = scratch.path(&format!("report-{name}.json"));
        let residuals = scratch.path(&format!("residuals-{name}"));
        let work = scratch.path(&format!("work-{name}"));
        let mut line = vec![
            "--reference",
            &reference,
            "--engine",
            &first,
            "--engine",
            &second,
            "--plugins",
            "UTRAnnotator,LOEUF",
            "--class-table",
            &table,
            "--work-dir",
            &work,
            "--report",
            &report,
            "--residuals",
            &residuals,
            "--pairing-threads",
            "2",
        ];
        for path in sweep {
            line.extend(["--utrannotator-sweep", path]);
        }
        let (code, _, error) = run_planted(&planted, &line);
        assert_eq!((code, error.as_str()), (0, ""), "{name}");
        let report = serde_json::from_str(&fs::read_to_string(&report).unwrap()).unwrap();
        (report, PathBuf::from(residuals), planted)
    };
    let item = |report: &Value, source: &str, engine: &str| -> Value {
        let found: Vec<&Value> = report["scores"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|item| {
                item["source"] == source
                    && item["engine"] == engine
                    && item["projection"] == "feature"
            })
            .collect();
        assert_eq!(found.len(), 1, "{source} {engine}");
        found[0].clone()
    };
    let counts = |item: &Value| -> [u64; 6] {
        [
            "perl",
            "rust",
            "intersection",
            "adj_perl",
            "adj_rust",
            "adj_intersection",
        ]
        .map(|key| item[key].as_u64().unwrap())
    };
    let (without, _, planted) = scored("without", &[]);
    assert!(without["utrannotator_sweep"].is_null());
    assert!(planted
        .handed
        .lock()
        .unwrap()
        .iter()
        .all(|(_, sweep)| sweep.is_none()));
    let utr = item(&without, "UTRAnnotator", "engine");
    assert_eq!(counts(&utr), [2, 2, 0, 2, 2, 0]);
    assert_eq!(utr["masked"], json!({"planted_swept_terms": 0}));
    let sweep: Vec<&str> = if as_directory {
        vec![&runs]
    } else {
        outputs.iter().map(String::as_str).collect()
    };
    let (report, residuals, planted) = scored("with", &sweep);
    let listed: Vec<Value> = outputs
        .iter()
        .map(|path| json!({"path": path, "sha256": sha256(&fs::read(path).unwrap())}))
        .collect();
    assert_eq!(report["utrannotator_sweep"], Value::Array(listed));
    for ident in ["engine", "again"] {
        let utr = item(&report, "UTRAnnotator", ident);
        assert_eq!(counts(&utr), [2, 2, 0, 1, 1, 0], "{ident}");
        assert_eq!(utr["masked"], json!({"planted_swept_terms": 1}), "{ident}");
        let left = fs::read_to_string(residuals.join(format!("UTRAnnotator.{ident}.feature.tsv")))
            .unwrap();
        assert_eq!(
            left.lines().collect::<Vec<_>>(),
            [
                format!("engine\t1:200\tT\tENST1\t5UTR_consequence\t{GAIN}"),
                format!("reference\t1:200\tT\tENST1\t5UTR_consequence\t{LOSS}"),
            ],
            "{ident}"
        );
        let zero = item(&report, "LOEUF", ident);
        assert_eq!(counts(&zero), [0, 1, 0, 0, 0, 0], "{ident}");
        assert_eq!(zero["masked"], json!({"planted_engine_zero": 1}), "{ident}");
    }
    let handed = planted.handed.lock().unwrap();
    let swept: Vec<Option<usize>> = handed
        .iter()
        .filter(|(klass, _)| klass == "planted_swept_terms")
        .map(|(_, sweep)| *sweep)
        .collect();
    assert_eq!(swept.len(), 4);
    assert!(swept[0].is_some() && swept.iter().all(|sweep| *sweep == swept[0]));
    let zero: Vec<Option<usize>> = handed
        .iter()
        .filter(|(klass, _)| klass == "planted_engine_zero")
        .map(|(_, sweep)| *sweep)
        .collect();
    assert_eq!(zero, [None, None]);
}

#[test]
fn seeded_runs_in_a_directory_are_handed_to_the_classes_of_utrannotator_alone() {
    seeded_runs_are_handed_to_the_classes_of_utrannotator_alone(true);
}

#[test]
fn seeded_runs_named_one_by_one_are_handed_to_the_classes_of_utrannotator_alone() {
    seeded_runs_are_handed_to_the_classes_of_utrannotator_alone(false);
}

#[test]
fn summary_says_that_classes_of_utrannotator_set_nothing_aside_without_the_seeded_runs() {
    let scratch = Scratch::new("summary");
    let (reference, engine, runs) = sweep_inputs(&scratch);
    let table = scratch.file("classes.tsv", SWEEP_TABLE);
    let engine = format!("engine={engine}");
    let notice =
        "no --utrannotator-sweep: a class of UTRAnnotator sets aside no row it reads the seeded runs for";
    let run = |name: &str, more: &[&str]| -> (i32, Vec<String>, String) {
        let report = scratch.path("report.json");
        let residuals = scratch.path(name);
        let work = scratch.path(&format!("work-{name}"));
        let mut line = vec![
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--plugins",
            "UTRAnnotator,LOEUF",
            "--class-table",
            &table,
            "--work-dir",
            &work,
            "--sort-memory",
            "16M",
            "--report",
            &report,
            "--residuals",
            &residuals,
            "--pairing-threads",
            "2",
        ];
        line.extend(more);
        run_planted(&Planted::default(), &line)
    };
    let (code, lines, error) = run("without", &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert_eq!(lines.last().map(String::as_str), Some(notice));
    let (seed_a, seed_b) = (
        format!("{runs}/seed-a.txt"),
        format!("{runs}/seed-b.txt.gz"),
    );
    let (code, lines, error) = run(
        "with",
        &[
            "--utrannotator-sweep",
            &seed_a,
            "--utrannotator-sweep",
            &seed_b,
        ],
    );
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!lines.iter().any(|line| line == notice));
    assert!(lines.iter().any(|line| line
        == "UTRAnnotator\tengine\tfeature\tclass planted_swept_terms set aside 1 [met, -]"));
    let (code, lines, error) = run("unloaded", &["--plugins", "LOEUF"]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!lines.iter().any(|line| line == notice));
    fs::create_dir(scratch.0.join("sweep/tools/none")).unwrap();
    scratch.file("sweep/tools/none/notes.txt", "no output here\n");
    let (missing, bare, empty, refused) = (
        scratch.path("nowhere.txt"),
        format!("{runs}/tools/.."),
        format!("{runs}/tools/none"),
        format!("{runs}/input.vcf"),
    );
    let cases = [
        ("missing", &missing, Some(format!("not a file: {missing}"))),
        ("bare", &bare, None),
        (
            "empty",
            &empty,
            Some(format!(
                "{empty}: no default-format output in the directory"
            )),
        ),
        (
            "refused",
            &refused,
            Some(format!(
                "{refused}: not a default-format output (columns Location, Allele, Feature)"
            )),
        ),
    ];
    for (name, path, message) in cases {
        let (code, lines, error) = run(name, &["--utrannotator-sweep", path]);
        let Some(message) = message else {
            assert_eq!((code, error.as_str()), (0, ""), "{name}");
            assert!(!lines.iter().any(|line| line == notice), "{name}");
            continue;
        };
        assert_eq!((code, lines.len()), (2, 0), "{name}");
        assert_eq!(
            error,
            format!("ERROR: [score_plugin_tuples] {message}\n"),
            "{name}"
        );
        assert!(
            !scratch.0.join(name).exists() && !scratch.0.join(format!("work-{name}")).exists(),
            "{name}"
        );
    }
}
