// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The four report files of hand-built comparisons, byte for byte: every count the report
//! derives from the masks over the union and per file, the release block, the markdown
//! lines and the two tables of one-sided tuples. One test runs the binary on a record of
//! several identifiers and reads the per-type table of the report it writes.

use std::process::Command;

use serde_json::{json, Value};
use vep_sv_compare::masks::{MaskResult, NamedCounts};
use vep_sv_compare::metrics::{Direction, DiscordantRecord, Metrics};
use vep_sv_compare::report::{json_report, markdown_report, write_reports, Reports};
use vep_sv_compare::rules::row_pair_rules;
use vep_sv_compare::strings::{Strings, Tuple, TupleSet};
use vep_sv_compare::TempDir;

/// The five texts of the fabricated tuple numbered `i`.
fn columns(i: u32) -> [String; 5] {
    [
        format!("21:{}-{}", 1000 + i, 2000 + i),
        "deletion".to_string(),
        format!("ENST{i:011}"),
        "Transcript".to_string(),
        "feature_truncation,intron_variant".to_string(),
    ]
}

/// The table of the texts of the tuples numbered below 1600.
fn table() -> Strings {
    Strings::from_texts((0..1600).flat_map(columns))
}

fn tuple(strings: &Strings, i: u32) -> Tuple {
    let texts = columns(i);
    strings
        .tuple(texts.each_ref().map(String::as_str))
        .expect("a tuple of the table")
}

fn layer<'a>(masks: &'a mut MaskResult, name: &str) -> &'a mut TupleSet {
    match name {
        "tsel_rust" => &mut masks.tsel_rust,
        "tsel_perl" => &mut masks.tsel_perl,
        "cnvtr_rust" => &mut masks.cnvtr_rust,
        "cnvtr_perl" => &mut masks.cnvtr_perl,
        "xchr_perl" => &mut masks.xchr_perl,
        "xchr_orphan_rust" => &mut masks.xchr_orphan_rust,
        "mate_rust" => &mut masks.mate_rust,
        "swap_rust" => &mut masks.swap_rust,
        "swap_perl" => &mut masks.swap_perl,
        "non_ref_rust" => &mut masks.non_ref_rust,
        "non_ref_perl" => &mut masks.non_ref_perl,
        "mate_read_rust" => &mut masks.mate_read_rust,
        "mate_read_perl" => &mut masks.mate_read_perl,
        "skipped_rust" => &mut masks.skipped_rust,
        "bnd_own_rust" => &mut masks.bnd_own_rust,
        "bnd_synonym_rust" => &mut masks.bnd_synonym_rust,
        "gvcf_rust" => &mut masks.gvcf_rust,
        "dot_pair_rust" => &mut masks.dot_pair_rust,
        "dot_pair_perl" => &mut masks.dot_pair_perl,
        "missing_alt_perl" => &mut masks.missing_alt_perl,
        "non_integer_end_perl" => &mut masks.non_integer_end_perl,
        "non_integer_end_rust" => &mut masks.non_integer_end_rust,
        "cnvtr_empty_rust" => &mut masks.cnvtr_empty_rust,
        "star_alone_perl" => &mut masks.star_alone_perl,
        "numified_pos_perl" => &mut masks.numified_pos_perl,
        "rn_zero_rust" => &mut masks.rn_zero_rust,
        "scored_engine_skip_residue" => &mut masks.scored_engine_skip_residue,
        other => panic!("no layer named {other}"),
    }
}

/// A mask result of `release` whose named layers hold the tuples numbered from the first
/// bound up to the second, the second left out.
fn masks(strings: &Strings, release: &str, layers: &[(&str, u32, u32)]) -> MaskResult {
    let mut result = MaskResult {
        release: release.to_string(),
        ..MaskResult::default()
    };
    for (name, from, to) in layers {
        *layer(&mut result, name) = (*from..*to).map(|i| tuple(strings, i)).collect();
    }
    result
}

const ROW_BUCKETS: [&str; 16] = [
    "structural_deletion_stop_lost_without_stop_codon",
    "stop_codon_completed_by_polyadenylation",
    "stop_codon_overlap_without_stop_codon",
    "stop_codon_behind_last_residue_without_stop_codon",
    "start_codon_overlap_without_start_codon",
    "stop_codon_overlap_genomic_span",
    "stop_retained_without_stop_codon_test",
    "stop_retained_three_base_stop_codon_test",
    "stop_codon_removed_without_stop_codon_test",
    "stop_codon_overlap_insertion_range",
    "stop_codon_settled_before_unknown_bases",
    "frameshift_inside_a_held_stop_codon",
    "splice_term_past_the_record_span",
    "pachysolen_nuclear_code_read_as_standard",
    "start_codon_removed_without_five_prime_utr",
    "start_codon_kept_without_five_prime_utr",
];

fn buckets(start: u64, covered: u64, rows: &[u64]) -> NamedCounts {
    let mut counts = vec![
        ("start_cooccurrence_swap".to_string(), start),
        ("covered_splice_region_swap".to_string(), covered),
    ];
    counts.extend(
        ROW_BUCKETS
            .iter()
            .zip(rows)
            .map(|(bucket, n)| (bucket.to_string(), *n)),
    );
    counts
}

fn kinds(cross: u64, same: u64) -> NamedCounts {
    vec![
        ("cross_chromosome".to_string(), cross),
        ("same_chromosome".to_string(), same),
    ]
}

fn row(counts: [u64; 8], ratios: [f64; 3]) -> Metrics {
    let [total_input, perl_variants, rust_variants, perl_tuples, rust_tuples, intersection, only_perl, only_rust] =
        counts;
    let [precision, recall, f1] = ratios;
    Metrics {
        total_input,
        perl_variants,
        rust_variants,
        perl_tuples,
        rust_tuples,
        intersection,
        only_perl,
        only_rust,
        precision,
        recall,
        f1,
        ..Metrics::default()
    }
}

/// A pooled row: the counts and ratios of [`row`], the three counts over the union of all
/// files and the two counts of tuples more than one file carries.
macro_rules! pooled {
    ($counts:expr, $ratios:expr, [$perl:expr, $rust:expr, $both:expr, $twice_perl:expr, $twice_rust:expr] $(,)?) => {
        Metrics {
            perl_tuples_union: Some($perl),
            rust_tuples_union: Some($rust),
            intersection_union: Some($both),
            cross_file_duplicate_perl: Some($twice_perl),
            cross_file_duplicate_rust: Some($twice_rust),
            ..row($counts, $ratios)
        }
    };
}

fn one_sided(
    strings: &Strings,
    source_file: &str,
    variant_type: &str,
    direction: Direction,
    i: u32,
) -> DiscordantRecord {
    DiscordantRecord {
        source_file: source_file.to_string(),
        variant_type: variant_type.to_string(),
        direction,
        tuple: tuple(strings, i),
    }
}

/// One compared file: its name, its raw and adjusted rows, its masks and the size of its
/// tandem-repeat class.
struct ComparedFile {
    name: &'static str,
    raw: Metrics,
    adjusted: Metrics,
    masks: MaskResult,
    cnvtr_total: (u64, u64),
}

/// One comparison as the report writer is handed it.
struct Case {
    release: &'static str,
    assembly: &'static str,
    scored_engine: &'static str,
    cache_supplied: bool,
    files: Vec<ComparedFile>,
    per_type: Vec<(String, Metrics)>,
    overall: Metrics,
    adjusted_overall: Metrics,
    masks: MaskResult,
    cnvtr_total: (u64, u64),
    discordants: Vec<DiscordantRecord>,
    adjusted_discordants: Vec<DiscordantRecord>,
}

/// The four files of a case, in the order the writer writes them.
struct Expected {
    json: &'static str,
    markdown: &'static str,
    discordant: &'static str,
    discordant_adjusted: &'static str,
}

fn assert_reports(strings: &Strings, case: Case, expected: Expected) {
    let named = |of: fn(&ComparedFile) -> Metrics| -> Vec<(String, Metrics)> {
        case.files
            .iter()
            .map(|file| (file.name.to_string(), of(file)))
            .collect()
    };
    let file_metrics = named(|file| file.raw.clone());
    let adjusted_file_metrics = named(|file| file.adjusted.clone());
    let masks_by_file: Vec<(String, MaskResult)> = case
        .files
        .iter()
        .map(|file| (file.name.to_string(), file.masks.clone()))
        .collect();
    let cnvtr_total_by_file: Vec<(String, (u64, u64))> = case
        .files
        .iter()
        .map(|file| (file.name.to_string(), file.cnvtr_total))
        .collect();
    let reports = Reports {
        strings,
        assembly: case.assembly,
        reference_release: case.release,
        scored_engine: case.scored_engine,
        cache_supplied: case.cache_supplied,
        file_metrics: &file_metrics,
        per_type: &case.per_type,
        overall: &case.overall,
        adjusted_overall: &case.adjusted_overall,
        adjusted_file_metrics: &adjusted_file_metrics,
        masks: &case.masks,
        masks_by_file: &masks_by_file,
        cnvtr_total: case.cnvtr_total,
        cnvtr_total_by_file: &cnvtr_total_by_file,
        discordants: &case.discordants,
        adjusted_discordants: &case.adjusted_discordants,
    };
    assert_eq!(json_report(&reports).unwrap(), expected.json);
    assert_eq!(markdown_report(&reports), expected.markdown);

    let scratch = TempDir::new("vep_sv_compare_reports_").unwrap();
    write_reports(&scratch.path, &reports).unwrap();
    let written = |name: &str| std::fs::read_to_string(scratch.path.join(name)).unwrap();
    assert_eq!(written("concordance_report.json"), expected.json);
    assert_eq!(written("concordance_report.md"), expected.markdown);
    assert_eq!(written("discordant.tsv"), expected.discordant);
    assert_eq!(
        written("discordant_adjusted.tsv"),
        expected.discordant_adjusted
    );
}

fn per_type(rows: Vec<(&str, Metrics)>) -> Vec<(String, Metrics)> {
    rows.into_iter()
        .map(|(name, m)| (name.to_string(), m))
        .collect()
}

/// A 116.2 comparison with a cache, every layer charged, three files compared in an order
/// that is not the sorted one, one of them named outside ASCII.
#[test]
fn every_layer_of_a_116_2_report_is_stated_over_the_union_and_per_file() {
    let registry: Vec<&str> = row_pair_rules().iter().map(|rule| rule.bucket).collect();
    assert_eq!(
        registry, ROW_BUCKETS,
        "the row-aware rules of the registry, in its order"
    );
    let strings = table();
    let s = &strings;
    let second = MaskResult {
        swap_by_bucket: buckets(1, 0, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        mate_read_by_kind: kinds(1, 0),
        skipped_records: 2,
        bare_reference_rows: 5,
        ..masks(
            s,
            "116.2",
            &[
                ("cnvtr_rust", 3, 5),
                ("cnvtr_perl", 20, 21),
                ("tsel_rust", 1, 2),
                ("xchr_perl", 11, 12),
                ("swap_rust", 1320, 1321),
                ("swap_perl", 30, 31),
                ("skipped_rust", 1400, 1403),
                ("scored_engine_skip_residue", 1520, 1521),
                ("mate_read_rust", 1340, 1341),
                ("mate_read_perl", 50, 51),
                ("non_ref_perl", 40, 42),
                ("xchr_orphan_rust", 1300, 1301),
                ("mate_rust", 1310, 1311),
                ("bnd_own_rust", 1420, 1421),
            ],
        )
    };
    let first = MaskResult {
        swap_by_bucket: buckets(0, 0, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        mate_read_by_kind: kinds(0, 0),
        ..masks(
            s,
            "116.2",
            &[
                ("cnvtr_rust", 5, 1240),
                ("tsel_rust", 2, 4),
                ("tsel_perl", 10, 12),
                ("bnd_synonym_rust", 1430, 1432),
            ],
        )
    };
    let union = MaskResult {
        xchr_unresolvable_scored: 7,
        swap_by_bucket: buckets(2, 0, &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        mate_read_by_kind: kinds(1, 1),
        skipped_records: 1_234_567,
        fasta_intergenic_records: 1001,
        bare_reference_rows: 2002,
        mate_overflow_records: 3003,
        mate_overflow_perl_tuples: 6,
        gvcf_records: 4004,
        ..masks(
            s,
            "116.2",
            &[
                ("tsel_rust", 1, 4),
                ("tsel_perl", 10, 12),
                ("cnvtr_rust", 3, 1240),
                ("cnvtr_perl", 20, 25),
                ("xchr_perl", 11, 14),
                ("xchr_orphan_rust", 1300, 1302),
                ("mate_rust", 1310, 1311),
                ("swap_rust", 1320, 1323),
                ("swap_perl", 30, 33),
                ("non_ref_rust", 1330, 1331),
                ("non_ref_perl", 40, 44),
                ("mate_read_rust", 1340, 1342),
                ("mate_read_perl", 50, 52),
                ("skipped_rust", 1400, 1410),
                ("bnd_own_rust", 1420, 1424),
                ("bnd_synonym_rust", 1430, 1435),
                ("gvcf_rust", 1440, 1446),
                ("dot_pair_rust", 1450, 1457),
                ("dot_pair_perl", 60, 68),
                ("missing_alt_perl", 70, 79),
                ("non_integer_end_perl", 80, 90),
                ("non_integer_end_rust", 1460, 1471),
                ("cnvtr_empty_rust", 1480, 1492),
                ("star_alone_perl", 100, 113),
                ("numified_pos_perl", 120, 134),
                ("rn_zero_rust", 1500, 1515),
                ("scored_engine_skip_residue", 1520, 1536),
            ],
        )
    };
    let second_row = row(
        [10, 8, 9, 100, 90, 80, 20, 10],
        [0.8888888888888888, 0.8, 0.8421052631578948],
    );
    let nothing = row([0; 8], [0.0; 3]);
    let case = Case {
        release: "116.2",
        assembly: "grch38",
        scored_engine: "fastvep",
        cache_supplied: true,
        files: vec![
            ComparedFile {
                name: "b_second",
                raw: second_row.clone(),
                adjusted: second_row,
                masks: second,
                cnvtr_total: (10, 3),
            },
            ComparedFile {
                name: "a_first",
                raw: row(
                    [1234, 1000, 999, 5000, 5100, 4900, 100, 200],
                    [0.9607843137254902, 0.98, 0.9702970297029703],
                ),
                adjusted: row(
                    [1234, 1000, 999, 4990, 4000, 3999, 991, 1],
                    [0.99975, 0.99995, 0.98999999],
                ),
                masks: first,
                cnvtr_total: (0, 5),
            },
            ComparedFile {
                name: "\u{e9}_third",
                raw: nothing.clone(),
                adjusted: nothing,
                masks: masks(s, "116.2", &[]),
                cnvtr_total: (0, 0),
            },
        ],
        per_type: per_type(vec![
            ("BND", row([4, 4, 4, 40, 40, 20, 20, 20], [0.5, 0.5, 0.5])),
            (
                "DEL",
                row([1000, 900, 900, 4000, 4000, 4000, 0, 0], [1.0, 1.0, 1.0]),
            ),
            ("DUP", row([3, 0, 0, 0, 0, 0, 0, 0], [0.0, 0.0, 0.0])),
            (
                "INS",
                row([30, 30, 30, 300, 300, 297, 3, 3], [0.99, 0.99, 0.99]),
            ),
            (
                "INV",
                row(
                    [7, 7, 7, 70, 70, 69, 1, 1],
                    [0.9857142857142858, 0.9857142857142858, 0.98999999],
                ),
            ),
            (
                "Unknown",
                row([1, 1, 1, 2, 3, 1, 1, 2], [0.00005, 0.000001, 0.0000005]),
            ),
        ]),
        overall: pooled!(
            [1244, 1008, 1008, 5100, 5190, 4980, 120, 210],
            [0.9595375722543352, 0.9764705882352941, 0.9679300291545189],
            [5090, 5185, 4975, 10, 5],
        ),
        adjusted_overall: pooled!(
            [1244, 1008, 1008, 5028, 3872, 3860, 1168, 12],
            [0.1234565, 0.123456789, 0.6666666666666666],
            [5020, 3870, 3859, 8, 2],
        ),
        masks: union,
        cnvtr_total: (2500, 30),
        discordants: vec![
            one_sided(s, "b_second", "DEL", Direction::OnlyPerl, 10),
            one_sided(s, "b_second", "Unknown", Direction::OnlyRust, 1),
            one_sided(s, "a_first", "BND", Direction::OnlyRust, 2),
        ],
        adjusted_discordants: vec![one_sided(s, "b_second", "DEL", Direction::OnlyPerl, 10)],
    };
    assert_reports(
        s,
        case,
        Expected {
            json: EVERY_LAYER_CHARGED_JSON,
            markdown: EVERY_LAYER_CHARGED_MD,
            discordant: EVERY_LAYER_CHARGED_TSV,
            discordant_adjusted: EVERY_LAYER_CHARGED_ADJUSTED_TSV,
        },
    );
}

/// A 115.2 comparison with exclusions and no cache: no release block, no assembly label
/// and no class total in the tandem-repeat line.
#[test]
fn a_115_2_report_keeps_its_shape() {
    let strings = table();
    let s = &strings;
    let no_pairs = |layers: &[(&str, u32, u32)]| MaskResult {
        swap_by_bucket: buckets(0, 0, &[]),
        mate_read_by_kind: kinds(0, 0),
        ..masks(s, "115.2", layers)
    };
    let pooled_row = [7, 7, 6, 59, 60, 56, 3, 4];
    let pooled_ratios = [0.9333333333333333, 0.9491525423728814, 0.9411764705882353];
    let case = Case {
        release: "115.2",
        assembly: "",
        scored_engine: "vep-rs",
        cache_supplied: false,
        files: vec![
            ComparedFile {
                name: "one",
                raw: row(
                    [5, 5, 5, 50, 52, 48, 2, 4],
                    [0.9230769230769231, 0.96, 0.9411764705882353],
                ),
                adjusted: row(
                    [5, 5, 5, 49, 50, 48, 1, 2],
                    [0.96, 0.9795918367346939, 0.9696969696969697],
                ),
                masks: no_pairs(&[
                    ("cnvtr_rust", 3, 5),
                    ("cnvtr_perl", 20, 22),
                    ("tsel_rust", 1, 3),
                    ("tsel_perl", 10, 11),
                ]),
                cnvtr_total: (0, 0),
            },
            ComparedFile {
                name: "two",
                raw: row(
                    [2, 2, 1, 9, 8, 8, 1, 0],
                    [1.0, 0.8888888888888888, 0.9411764705882353],
                ),
                adjusted: row([2, 2, 1, 9, 7, 7, 2, 0], [1.0, 0.7777777777777778, 0.875]),
                masks: no_pairs(&[("cnvtr_rust", 3, 4)]),
                cnvtr_total: (0, 0),
            },
        ],
        per_type: per_type(vec![("DEL", row(pooled_row, pooled_ratios))]),
        overall: pooled!(pooled_row, pooled_ratios, [59, 60, 56, 0, 0]),
        adjusted_overall: pooled!(
            [7, 7, 6, 58, 57, 55, 3, 2],
            [0.9649122807017544, 0.9482758620689655, 0.9565217391304348],
            [58, 57, 55, 0, 0],
        ),
        masks: no_pairs(&[
            ("tsel_rust", 1, 3),
            ("tsel_perl", 10, 11),
            ("cnvtr_rust", 3, 5),
            ("cnvtr_perl", 20, 22),
        ]),
        cnvtr_total: (0, 0),
        discordants: vec![
            one_sided(s, "one", "DEL", Direction::OnlyPerl, 10),
            one_sided(s, "one", "DEL", Direction::OnlyRust, 7),
        ],
        adjusted_discordants: vec![one_sided(s, "one", "DEL", Direction::OnlyRust, 7)],
    };
    assert_reports(
        s,
        case,
        Expected {
            json: PUBLISHED_SHAPE_JSON,
            markdown: PUBLISHED_SHAPE_MD,
            discordant: PUBLISHED_SHAPE_TSV,
            discordant_adjusted: PUBLISHED_SHAPE_ADJUSTED_TSV,
        },
    );
}

/// Files compared and nothing set aside: every per-file map is written, empty, and the
/// markdown has no adjusted line.
#[test]
fn a_report_with_nothing_set_aside_writes_empty_per_file_maps() {
    let strings = table();
    let s = &strings;
    let no_pairs = || MaskResult {
        swap_by_bucket: buckets(0, 0, &[]),
        mate_read_by_kind: kinds(0, 0),
        ..masks(s, "115.2", &[])
    };
    let counts = [3, 3, 3, 12, 12, 12, 0, 0];
    let case = Case {
        release: "115.2",
        assembly: "grch37",
        scored_engine: "vep-rs",
        cache_supplied: true,
        files: vec![ComparedFile {
            name: "x",
            raw: row(counts, [1.0; 3]),
            adjusted: row(counts, [1.0; 3]),
            masks: no_pairs(),
            cnvtr_total: (0, 0),
        }],
        per_type: per_type(vec![("SNV", row(counts, [1.0; 3]))]),
        overall: pooled!(counts, [1.0; 3], [12, 12, 12, 0, 0]),
        adjusted_overall: pooled!(counts, [1.0; 3], [12, 12, 12, 0, 0]),
        masks: no_pairs(),
        cnvtr_total: (0, 0),
        discordants: Vec::new(),
        adjusted_discordants: Vec::new(),
    };
    assert_reports(
        s,
        case,
        Expected {
            json: NOTHING_SET_ASIDE_JSON,
            markdown: NOTHING_SET_ASIDE_MD,
            discordant: NOTHING_SET_ASIDE_TSV,
            discordant_adjusted: NOTHING_SET_ASIDE_ADJUSTED_TSV,
        },
    );
}

/// A report over no file leaves five of the per-file maps out and writes three, empty.
#[test]
fn a_report_over_no_file_leaves_five_per_file_maps_out() {
    let strings = table();
    let s = &strings;
    let case = Case {
        release: "115.2",
        assembly: "grch37",
        scored_engine: "vep-rs",
        cache_supplied: false,
        files: Vec::new(),
        per_type: Vec::new(),
        overall: row([0; 8], [0.0; 3]),
        adjusted_overall: row([0; 8], [0.0; 3]),
        masks: masks(s, "115.2", &[]),
        cnvtr_total: (0, 0),
        discordants: Vec::new(),
        adjusted_discordants: Vec::new(),
    };
    assert_reports(
        s,
        case,
        Expected {
            json: NO_FILE_JSON,
            markdown: NO_FILE_MD,
            discordant: NO_FILE_TSV,
            discordant_adjusted: NO_FILE_ADJUSTED_TSV,
        },
    );
}

/// The binary on one record whose ID column holds two identifiers. Both outputs name its
/// rows by the first, at the Location of its span: the record's row of the per-type table
/// holds their tuples, and each one-sided tuple reads the record's type.
#[test]
fn the_binary_types_the_rows_of_a_record_with_several_identifiers() {
    let scratch = TempDir::new("vep_sv_compare_ids_").unwrap();
    let path = |name: &str| scratch.path.join(name);
    for directory in ["in", "perl", "rust", "out"] {
        std::fs::create_dir(path(directory)).unwrap();
    }
    std::fs::write(
        path("in/ids.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t1000\tsv_one;sv_two\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000\n",
    )
    .unwrap();
    let output = |features: [&str; 2]| {
        let mut text = String::from(
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n",
        );
        for feature in features {
            text.push_str(&format!(
                "sv_one\t21:1001-2000\tdeletion\tG\t{feature}\tTranscript\tfeature_truncation\n"
            ));
        }
        text
    };
    std::fs::write(path("perl/ids.txt"), output(["ENST_A", "ENST_B"])).unwrap();
    std::fs::write(path("rust/ids.txt"), output(["ENST_A", "ENST_C"])).unwrap();
    let run = Command::new(env!("CARGO_BIN_EXE_vep-sv-compare"))
        .args(["--reference-release", "115.2", "--assembly", "grch37"])
        .arg("--input-dir")
        .arg(path("in"))
        .arg("--perl-dir")
        .arg(path("perl"))
        .arg("--rust-dir")
        .arg(path("rust"))
        .arg("--output-dir")
        .arg(path("out"))
        .output()
        .unwrap();
    assert_eq!(
        run.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let read = |name: &str| std::fs::read_to_string(path("out").join(name)).unwrap();
    let report: Value = serde_json::from_str(&read("concordance_report.json")).unwrap();
    assert_eq!(
        report["per_variant_type"],
        json!({
            "Symbolic_DEL": {
                "total_input": 1,
                "perl_variants": 0,
                "rust_variants": 0,
                "perl_tuples": 2,
                "rust_tuples": 2,
                "intersection": 1,
                "only_perl": 1,
                "only_rust": 1,
                "precision": 0.5,
                "recall": 0.5,
                "f1": 0.5
            }
        })
    );
    assert!(read("concordance_report.md")
        .contains("\n| Symbolic_DEL | 1 | 2 | 2 | 1 | 1 | 0.5000 ** |\n"));
    let one_sided = read("discordant.tsv");
    let types: Vec<&str> = one_sided
        .lines()
        .skip(1)
        .map(|line| line.split('\t').nth(1).unwrap())
        .collect();
    assert_eq!(types, ["Symbolic_DEL", "Symbolic_DEL"]);
}

const EVERY_LAYER_CHARGED_JSON: &str = r##"{
  "overall": {
    "total_input": 1244,
    "perl_variants": 1008,
    "rust_variants": 1008,
    "perl_tuples": 5100,
    "rust_tuples": 5190,
    "intersection": 4980,
    "only_perl": 120,
    "only_rust": 210,
    "precision": 0.959538,
    "recall": 0.976471,
    "f1": 0.96793,
    "perl_tuples_union": 5090,
    "rust_tuples_union": 5185,
    "intersection_union": 4975,
    "cross_file_duplicate_perl": 10,
    "cross_file_duplicate_rust": 5
  },
  "per_file": {
    "b_second": {
      "total_input": 10,
      "perl_variants": 8,
      "rust_variants": 9,
      "perl_tuples": 100,
      "rust_tuples": 90,
      "intersection": 80,
      "only_perl": 20,
      "only_rust": 10,
      "precision": 0.888889,
      "recall": 0.8,
      "f1": 0.842105
    },
    "a_first": {
      "total_input": 1234,
      "perl_variants": 1000,
      "rust_variants": 999,
      "perl_tuples": 5000,
      "rust_tuples": 5100,
      "intersection": 4900,
      "only_perl": 100,
      "only_rust": 200,
      "precision": 0.960784,
      "recall": 0.98,
      "f1": 0.970297
    },
    "\u00e9_third": {
      "total_input": 0,
      "perl_variants": 0,
      "rust_variants": 0,
      "perl_tuples": 0,
      "rust_tuples": 0,
      "intersection": 0,
      "only_perl": 0,
      "only_rust": 0,
      "precision": 0.0,
      "recall": 0.0,
      "f1": 0.0
    }
  },
  "per_variant_type": {
    "BND": {
      "total_input": 4,
      "perl_variants": 4,
      "rust_variants": 4,
      "perl_tuples": 40,
      "rust_tuples": 40,
      "intersection": 20,
      "only_perl": 20,
      "only_rust": 20,
      "precision": 0.5,
      "recall": 0.5,
      "f1": 0.5
    },
    "DEL": {
      "total_input": 1000,
      "perl_variants": 900,
      "rust_variants": 900,
      "perl_tuples": 4000,
      "rust_tuples": 4000,
      "intersection": 4000,
      "only_perl": 0,
      "only_rust": 0,
      "precision": 1.0,
      "recall": 1.0,
      "f1": 1.0
    },
    "DUP": {
      "total_input": 3,
      "perl_variants": 0,
      "rust_variants": 0,
      "perl_tuples": 0,
      "rust_tuples": 0,
      "intersection": 0,
      "only_perl": 0,
      "only_rust": 0,
      "precision": 0.0,
      "recall": 0.0,
      "f1": 0.0
    },
    "INS": {
      "total_input": 30,
      "perl_variants": 30,
      "rust_variants": 30,
      "perl_tuples": 300,
      "rust_tuples": 300,
      "intersection": 297,
      "only_perl": 3,
      "only_rust": 3,
      "precision": 0.99,
      "recall": 0.99,
      "f1": 0.99
    },
    "INV": {
      "total_input": 7,
      "perl_variants": 7,
      "rust_variants": 7,
      "perl_tuples": 70,
      "rust_tuples": 70,
      "intersection": 69,
      "only_perl": 1,
      "only_rust": 1,
      "precision": 0.985714,
      "recall": 0.985714,
      "f1": 0.99
    },
    "Unknown": {
      "total_input": 1,
      "perl_variants": 1,
      "rust_variants": 1,
      "perl_tuples": 2,
      "rust_tuples": 3,
      "intersection": 1,
      "only_perl": 1,
      "only_rust": 2,
      "precision": 5e-05,
      "recall": 1e-06,
      "f1": 0.0
    }
  },
  "adjusted": {
    "total_input": 1244,
    "perl_variants": 1008,
    "rust_variants": 1008,
    "perl_tuples": 5028,
    "rust_tuples": 3872,
    "intersection": 3860,
    "only_perl": 1168,
    "only_rust": 12,
    "precision": 0.123456,
    "recall": 0.123457,
    "f1": 0.666667,
    "perl_tuples_union": 5020,
    "rust_tuples_union": 3870,
    "intersection_union": 3859,
    "cross_file_duplicate_perl": 8,
    "cross_file_duplicate_rust": 2,
    "excluded_divergences": 1318,
    "excluded_rust_extra": 1318,
    "excluded_perl_extra": 72,
    "excluded_cnvtr_rust": 1237,
    "excluded_cnvtr_perl": 5,
    "excluded_transcript_selection_rust": 3,
    "excluded_transcript_selection_perl": 2,
    "transcript_selection_cnvtr_overlap_rust": 1,
    "cnvtr_swap_pairs_total_rust": 2500,
    "cnvtr_swap_pairs_total_perl": 30,
    "excluded_cross_chromosome_perl": 3,
    "excluded_cross_chromosome_rust": 2,
    "excluded_transcript_selection_mate_rust": 1,
    "excluded_registry_swap_rust": 3,
    "excluded_registry_swap_perl": 3,
    "excluded_registry_swap_by_bucket": {
      "start_cooccurrence_swap": 2,
      "covered_splice_region_swap": 0,
      "structural_deletion_stop_lost_without_stop_codon": 0,
      "stop_codon_completed_by_polyadenylation": 0,
      "stop_codon_overlap_without_stop_codon": 0,
      "stop_codon_behind_last_residue_without_stop_codon": 0,
      "start_codon_overlap_without_start_codon": 0,
      "stop_codon_overlap_genomic_span": 1,
      "stop_retained_without_stop_codon_test": 0,
      "stop_retained_three_base_stop_codon_test": 0,
      "stop_codon_removed_without_stop_codon_test": 0,
      "stop_codon_overlap_insertion_range": 0,
      "stop_codon_settled_before_unknown_bases": 0,
      "frameshift_inside_a_held_stop_codon": 0,
      "splice_term_past_the_record_span": 0,
      "pachysolen_nuclear_code_read_as_standard": 0,
      "start_codon_removed_without_five_prime_utr": 0,
      "start_codon_kept_without_five_prime_utr": 0
    },
    "excluded_non_ref_batch_rust": 1,
    "excluded_non_ref_batch_perl": 4,
    "excluded_breakend_mate_local_read_rust": 2,
    "excluded_breakend_mate_local_read_perl": 2,
    "excluded_breakend_mate_local_read_by_kind": {
      "cross_chromosome": 1,
      "same_chromosome": 1
    },
    "transcript_selection_cross_chromosome_overlap_perl": 1,
    "cross_chromosome_check_ran": true,
    "scored_engine": "fastvep",
    "scored_engine_transcripts_absent_from_cache": 7,
    "excluded_cross_chromosome_by_file": {
      "b_second": 1
    },
    "excluded_cnvtr_by_file": {
      "a_first": 1235,
      "b_second": 2
    },
    "excluded_transcript_selection_by_file": {
      "a_first": {
        "rust": 2,
        "perl": 2
      },
      "b_second": {
        "rust": 1,
        "perl": 0
      }
    },
    "cnvtr_swap_pairs_total_by_file": {
      "a_first": {
        "rust": 0,
        "perl": 5
      },
      "b_second": {
        "rust": 10,
        "perl": 3
      }
    },
    "excluded_cross_chromosome_rust_by_file": {
      "b_second": 1
    },
    "excluded_transcript_selection_mate_by_file": {
      "b_second": 1
    },
    "excluded_registry_swap_by_file": {
      "b_second": {
        "rust": 1,
        "perl": 1,
        "start_cooccurrence_swap": 1,
        "covered_splice_region_swap": 0,
        "structural_deletion_stop_lost_without_stop_codon": 0,
        "stop_codon_completed_by_polyadenylation": 0,
        "stop_codon_overlap_without_stop_codon": 0,
        "stop_codon_behind_last_residue_without_stop_codon": 0,
        "start_codon_overlap_without_start_codon": 0,
        "stop_codon_overlap_genomic_span": 0,
        "stop_retained_without_stop_codon_test": 0,
        "stop_retained_three_base_stop_codon_test": 0,
        "stop_codon_removed_without_stop_codon_test": 0,
        "stop_codon_overlap_insertion_range": 0,
        "stop_codon_settled_before_unknown_bases": 0,
        "frameshift_inside_a_held_stop_codon": 0,
        "splice_term_past_the_record_span": 0,
        "pachysolen_nuclear_code_read_as_standard": 0,
        "start_codon_removed_without_five_prime_utr": 0,
        "start_codon_kept_without_five_prime_utr": 0
      }
    },
    "excluded_non_ref_batch_by_file": {
      "b_second": {
        "rust": 0,
        "perl": 2
      }
    },
    "excluded_breakend_mate_local_read_by_file": {
      "b_second": {
        "rust": 1,
        "perl": 1,
        "cross_chromosome": 1,
        "same_chromosome": 0
      }
    },
    "reference_release": "116.2",
    "masks_applied": [
      "filter_cnv_tr_expansion_divergences",
      "filter_registry_swap_pairs",
      "filter_breakend_mate_local_read_pairs",
      "filter_reference_skipped_records",
      "filter_bnd_own_rows_lost",
      "filter_bnd_synonym_mate_rows",
      "filter_bnd_bare_reference_allele_rows",
      "filter_gvcf_variant_sites",
      "filter_dot_breakend_pairs",
      "filter_missing_alt_annotated",
      "filter_non_integer_end_rows",
      "filter_cnv_tr_empty_expansion",
      "filter_star_alone_svtype_annotated",
      "filter_non_numeric_pos_numified",
      "filter_cnv_tr_rn_zero_allele_dropped"
    ],
    "masks_not_applied": {
      "filter_intended_divergences": "116 drops a vep_skip record before the input buffer; the reference writes no row to compare transcript sets against",
      "filter_giant_breakend_mate_divergences": "116 drops a breakend pair over --max_sv_size whole; its mate allele has no reference row",
      "filter_non_ref_batch_divergences": "116 drops a <NON_REF> record before the input buffer",
      "filter_cross_chromosome_divergences": "116.1's chromosome-aware overlap never pairs a record with a transcript on another chromosome",
      "filter_cross_chromosome_orphan_intergenic": "with no off-chromosome rows to displace it, the reference writes the intergenic row itself"
    },
    "excluded_reference_skipped_rust": 10,
    "reference_skipped_records": 1234567,
    "excluded_reference_skipped_by_file": {
      "b_second": {
        "rust": 3,
        "records": 2
      }
    },
    "excluded_bnd_own_rows_lost_rust": 4,
    "excluded_bnd_own_rows_lost_by_file": {
      "b_second": 1
    },
    "excluded_bnd_synonym_mate_rows_rust": 5,
    "excluded_bnd_synonym_mate_rows_by_file": {
      "a_first": 2
    },
    "fasta_named_slice_intergenic_records": 1001,
    "bnd_bare_reference_allele_rows": 2002,
    "bnd_bare_reference_allele_rows_by_file": {
      "b_second": 5
    },
    "bnd_mate_32bit_buffer_records": 3003,
    "bnd_mate_32bit_buffer_reference_tuples": 6,
    "excluded_gvcf_variant_site_rust": 6,
    "gvcf_variant_site_records": 4004,
    "excluded_dot_breakend_pair_rust": 7,
    "excluded_dot_breakend_pair_perl": 8,
    "excluded_missing_alt_annotated_perl": 9,
    "excluded_non_integer_end_perl": 10,
    "excluded_non_integer_end_rust": 11,
    "excluded_cnv_tr_empty_expansion_rust": 12,
    "excluded_star_alone_svtype_perl": 13,
    "excluded_non_numeric_pos_perl": 14,
    "excluded_cnv_tr_rn_zero_rust": 15,
    "scored_engine_skip_route_residue_rows": 16,
    "scored_engine_skip_route_residue_by_file": {
      "b_second": 1
    },
    "breakend_masks_check_ran": true,
    "divergence_classes": {
      "reference_skipped_record": {
        "definition": "A structural record the reference dropped before its input buffer, marked vep_skip for a span over --max_sv_size (its SVLEN, else its END), an ALT type with no Sequence Ontology term, an unsupported multi-allelic list or a deletion with neither END nor SVLEN, so it has no reference row in any format; vep-rs's rows for the record are set aside.",
        "excludes": true,
        "excluded": 10,
        "taxonomy_class": "reference_skipped_record",
        "perl_citation": "Parser.pm:491-501, 821-827; Parser/VCF.pm:477-481, 555-558, 575"
      },
      "bnd_own_chromosome_rows_lost": {
        "definition": "A bracket breakend whose mate is on another chromosome loses the rows of its own end when a mate-chromosome transcript, loaded for an earlier record of the same buffer, reaches it first and names its slice: the reference keeps the mate rows and writes nothing under the own-end allele, an order-dependent result.",
        "excludes": true,
        "excluded": 4,
        "taxonomy_class": "bnd_own_chromosome_rows_lost",
        "perl_citation": "InputBuffer.pm:191-207, 223-234; AnnotationSource.pm:235-258; AnnotationType/Transcript.pm:122-143"
      },
      "bnd_synonym_mate_rows_lost": {
        "definition": "A bracket breakend whose mate's chromosome is spelled with a prefix the mate-allele test does not resolve, Chr21 for 21, loses its mate rows in the reference and keeps them in vep-rs; a mate named by an accession on another chromosome loses them in the reference alone too, vep-rs resolving the accession through the synonym table, and one on the record's own chromosome keeps them in both engines.",
        "excludes": true,
        "excluded": 5,
        "taxonomy_class": "bnd_synonym_mate_rows_lost",
        "perl_citation": "StructuralVariationOverlap.pm:86-103, 144-160; Utils/VariationEffect.pm:111-124; InputBuffer.pm:348-357; AnnotationType/Transcript.pm:134"
      },
      "fasta_named_slice_intergenic": {
        "definition": "Under --check_ref or --lookup_ref with a FASTA whose sequence names differ from the cache's, the variant's slice is named after the FASTA and every transcript fails the slice test, so the reference writes intergenic_variant where vep-rs writes transcript rows; a property of the run, counted and never set aside.",
        "excludes": false,
        "excluded": 0,
        "taxonomy_class": "fasta_named_slice_intergenic",
        "perl_citation": "Parser.pm:749; BaseVEP.pm get_slice; AnnotationType/Transcript.pm:141"
      },
      "bnd_bare_reference_allele_row": {
        "definition": "A breakpoint record that names no mate (a bare <BND>, or one with SVLEN, POS2 or an unreadable END2 and no bracket) is held as N/<BND> and, with no breakend parsed, one allele is built per slash-separated part, the reference base N included, so the reference writes every row twice under chromosome_breakpoint, ALLELE_NUM 1 for the reference allele and 2 for the ALT; a reference allele has no variant effect, so vep-rs writes the ALT row once and the reference-allele rows are set aside: each duplicates the ALT row's tuple, so the unique-tuple key stream drops it before the raw and the adjusted figure, and both comparators count the rows so set aside per record.",
        "excludes": true,
        "excluded": 2002,
        "taxonomy_class": "bnd_bare_reference_allele_row",
        "perl_citation": "Parser/VCF.pm:499; StructuralVariationOverlap.pm:64-77; OutputFactory.pm BaseStructuralVariationOverlapAllele_to_output_hash"
      },
      "bnd_mate_32bit_buffer_intergenic": {
        "definition": "A coordinate of the reference's input buffer that, plus the --distance flank, reaches 2^31 (a breakend's mate position, a span's end or the record's own position) turns every record of the buffer intergenic, that record and the records loaded with it (the first such position is 2147478648 under the default flank of 5000 and 2147483647 under --distance 1); the coordinates enter the buffer's bounds with the flank added when the cache regions are chosen. vep-rs annotates the buffer's records against their transcripts and a mate at its coordinate. Counted on the reference's rows at the record's Location and never set aside, since the records it reaches beyond that one depend on the buffer boundaries.",
        "excludes": false,
        "excluded": 0,
        "taxonomy_class": "bnd_mate_32bit_buffer_intergenic",
        "perl_citation": "AnnotationSource.pm:190, 251-257, 306-309; InputBuffer.pm:386-390"
      },
      "gvcf_variant_site_dropped": {
        "definition": "A gVCF variant site lists <NON_REF> beside its real allele (`A G,<NON_REF>`); the <NON_REF> token routes the line to the structural path, where `G/<NON_REF>` has no Sequence Ontology term, and the reference drops the site whole, the real allele with it. vep-rs annotates the real allele and reads <NON_REF> as the absence of a variant. vep-rs's rows at the site are set aside: recognised here by the record's own name (`CHROM_POS_REF/ALT/<NON_REF>`) at a Location the reference wrote nothing for, and in the structural-variant comparator from the input record.",
        "excludes": true,
        "excluded": 6,
        "taxonomy_class": "gvcf_variant_site_dropped",
        "perl_citation": "Parser/VCF.pm:237-245, 469-481; Parser.pm:663-728"
      },
      "dot_breakend_without_svtype_pair": {
        "definition": "A dot-notation breakend without INFO/SVTYPE (`N.`, `.N`) takes the sequence rules for its coordinates: the start is POS + 1 when the ALT's length differs from the REF's and POS when it does not, the end is POS + len(REF) - 1, so the reference writes the pair `POS-(POS+1)` with the positions of an insertion for a one-base REF, the span `(POS+1)-(POS+len(REF)-1)` for a longer REF beside an ALT of another length and `POS-(POS+len(REF)-1)` beside an ALT of its own length; with `SVTYPE=BND` it writes the single base at POS + 1, as vep-rs writes it for every REF. The reference's span row and vep-rs's point row at the span's start or the base after it on the same feature are set aside together.",
        "excludes": true,
        "excluded": 15,
        "taxonomy_class": "dot_breakend_without_svtype_pair",
        "perl_citation": "BaseVCF4.pm:255-276, 297-317"
      },
      "missing_alt_annotated": {
        "definition": "A record whose ALT column is the missing value `.` is routed to the structural path by INFO/SVTYPE and annotated as an allele: a `chromosome_breakpoint` for `SVTYPE=BND`, a deletion over END for `SVTYPE=DEL`; without SVTYPE the reference writes nothing for it, as vep-rs writes nothing in both cases. The reference's rows are set aside by the structural-variant comparator, which reads the input; the SNP/indel comparator cannot tell them from a record with an allele.",
        "excludes": true,
        "excluded": 9,
        "taxonomy_class": "missing_alt_annotated",
        "perl_citation": "BaseVCF4.pm:486-490; Parser/VCF.pm:237-245; Parser.pm:710-711"
      },
      "non_integer_end_annotated": {
        "definition": "An INFO/END, SVLEN or POS that is not an integer (`END=30425949.5`, `SVLEN=300.7`, `SVLEN=NaN`, POS `30425949.0` on the sequence path) is taken as written, passes `looks_like_number`, and the record is annotated with the text in its Location (`21:30425651-30425950.7`, `21:30425651-NaN`, `21:30425949.0`) and a fractional DISTANCE; vep-rs refuses an END that is not a coordinate and reads a SVLEN or POS at its numified value (`21:30425651-30425950`, the pair `21:30425650-30425651`, `21:30425949`). The reference's rows, recognised by the Location, are set aside in both comparators, and vep-rs's rows of the same record by the structural-variant comparator, which reads the input.",
        "excludes": true,
        "excluded": 21,
        "taxonomy_class": "non_integer_end_annotated",
        "perl_citation": "BaseVCF4.pm:263, 302-309; Parser.pm:475"
      },
      "cnv_tr_empty_expansion_dropped": {
        "definition": "With a FASTA, a <CNV:TR> whose expansion yields no alternate allele, or one equal to the reference run (`RN=.`, a RUC Perl reads as true that gives zero or fewer repeats, a RUC times unit length equal to SVLEN), is expanded to an allele string with no alternate and the reference writes nothing for it, without a warning; vep-rs with a FASTA writes nothing for it either, and without one writes the symbolic rows over the run, as the reference does without one. vep-rs's rows are set aside by the structural-variant comparator when the input record has the shape and the reference wrote nothing at or beside the run; the SNP/indel comparator does not judge it.",
        "excludes": true,
        "excluded": 12,
        "taxonomy_class": "cnv_tr_empty_expansion_dropped",
        "perl_citation": "Parser/VCF.pm:382-431, 516-544; Parser.pm:553-592"
      },
      "star_alone_svtype_annotated": {
        "definition": "A record whose every ALT is `*` under a true INFO/SVTYPE is routed to the structural path and annotated from INFO alone, a deletion over END for `SVTYPE=DEL` or a `chromosome_breakpoint` at POS + 1 for `SVTYPE=BND`, although `*` names no allele of the record; vep-rs writes no row for it. The reference's rows are set aside by the structural-variant comparator, which reads the input; the SNP/indel comparator cannot tell them from a record with an allele.",
        "excludes": true,
        "excluded": 13,
        "taxonomy_class": "star_alone_svtype_annotated",
        "perl_citation": "Parser/VCF.pm:237-245, 469-475; StructuralVariationOverlap.pm:64-77"
      },
      "non_numeric_pos_numified": {
        "definition": "A POS that is not a non-negative number (`.`, empty, `0x1D04A1D`, fullwidth digits, `30_425_949`, `-5`) is numified by the structural path's increment, `.` to 0 and the record annotated from base 1, `30_425_949` from 31, `-5` from -4, where vep-rs refuses the line; the reference's rows at the record's start are set aside by the structural-variant comparator from the input, and by the SNP/indel comparator when the Location's first coordinate is negative, the one spelling no written POS reaches.",
        "excludes": true,
        "excluded": 14,
        "taxonomy_class": "non_numeric_pos_numified",
        "perl_citation": "BaseVCF4.pm:255-264; Parser.pm:475-488"
      },
      "cnv_tr_rn_zero_allele_dropped": {
        "definition": "With a FASTA, a two-ALT <CNV:TR> of which one allele has `RN=0` loses an allele: the empty expansion is an empty part of the allele string, which is not read as multi-allelic, so the first allele alone is kept when the empty part leads and the trailing empty part is dropped by the split; vep-rs writes both alleles when the empty part leads, and the structural-variant comparator sets aside its rows of the dropped allele: the symbolic rows told by their direction term when the real allele is a gain, or unchanged beside a leading empty part (vep-rs reads a symbolic <CNV:TR> allele as a gain unless its record states it shorter than its run or, with no direction term, as long as it), and the literal rows a run with a FASTA writes for the real allele behind a leading empty part.",
        "excludes": true,
        "excluded": 15,
        "taxonomy_class": "cnv_tr_rn_zero_allele_dropped",
        "perl_citation": "Parser/VCF.pm:411, 428; Parser.pm:881, 978"
      }
    }
  },
  "per_file_adjusted": {
    "b_second": {
      "total_input": 10,
      "perl_variants": 8,
      "rust_variants": 9,
      "perl_tuples": 100,
      "rust_tuples": 90,
      "intersection": 80,
      "only_perl": 20,
      "only_rust": 10,
      "precision": 0.888889,
      "recall": 0.8,
      "f1": 0.842105
    },
    "a_first": {
      "total_input": 1234,
      "perl_variants": 1000,
      "rust_variants": 999,
      "perl_tuples": 4990,
      "rust_tuples": 4000,
      "intersection": 3999,
      "only_perl": 991,
      "only_rust": 1,
      "precision": 0.99975,
      "recall": 0.99995,
      "f1": 0.99
    },
    "\u00e9_third": {
      "total_input": 0,
      "perl_variants": 0,
      "rust_variants": 0,
      "perl_tuples": 0,
      "rust_tuples": 0,
      "intersection": 0,
      "only_perl": 0,
      "only_rust": 0,
      "precision": 0.0,
      "recall": 0.0,
      "f1": 0.0
    }
  }
}
"##;

const EVERY_LAYER_CHARGED_MD: &str = r##"# SV Concordance Report

## Overall

- **Reference**: Ensembl VEP 116.2; masks applied: `filter_cnv_tr_expansion_divergences`, `filter_registry_swap_pairs`, `filter_breakend_mate_local_read_pairs`, `filter_reference_skipped_records`, `filter_bnd_own_rows_lost`, `filter_bnd_synonym_mate_rows`, `filter_bnd_bare_reference_allele_rows`, `filter_gvcf_variant_sites`, `filter_dot_breakend_pairs`, `filter_missing_alt_annotated`, `filter_non_integer_end_rows`, `filter_cnv_tr_empty_expansion`, `filter_star_alone_svtype_annotated`, `filter_non_numeric_pos_numified`, `filter_cnv_tr_rn_zero_allele_dropped`
- **116.2 classes**: `reference_skipped_record` 10 vep-rs rows on 1,234,567 records set aside; `bnd_own_chromosome_rows_lost` 4 set aside; `bnd_synonym_mate_rows_lost` 5 set aside; `fasta_named_slice_intergenic` 1,001 records counted, none set aside; `bnd_bare_reference_allele_row` 2,002 reference-allele rows set aside; `bnd_mate_32bit_buffer_intergenic` 3,003 records counted, none set aside; `gvcf_variant_site_dropped` 6 vep-rs rows on 4,004 sites set aside; `dot_breakend_without_svtype_pair` 8 reference and 7 vep-rs rows set aside; `missing_alt_annotated` 9 reference rows set aside; `non_integer_end_annotated` 10 reference and 11 vep-rs rows set aside; `cnv_tr_empty_expansion_dropped` 12 vep-rs rows set aside; `star_alone_svtype_annotated` 13 reference rows set aside; `non_numeric_pos_numified` 14 reference rows set aside; `cnv_tr_rn_zero_allele_dropped` 15 vep-rs rows set aside; `scored_engine_skip_route_residue` 16 scored-engine rows at Locations the reference wrote nothing for counted, none set aside
- **F1 (adjusted)**: 0.6667  *(excluded 1,318 Rust-extra + 72 Perl-extra transcript divergences)*
- **F1 (raw)**: 0.9679
- **Total input variants**: 1244
- **Perl tuples**: 5100
- **Rust tuples**: 5190
- **Intersection**: 4980
- **Only Perl**: 120
- **Only Rust**: 210

SV adjusted F1 additionally excludes the `<CNV:TR>` literal-expansion class (1,237 of 2,500 pairs GRCH38, a_first, b_second only). Raw F1 unaffected.

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|
| a_first | 1234 | 5000 | 5100 | 4900 | 0.9703 | 0.9900 |
| b_second | 10 | 100 | 90 | 80 | 0.8421 |  |
| é_third | 0 | 0 | 0 | 0 | 0.0000 |  |

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|
| BND | 4 | 40 | 40 | 20 | 20 | 0.5000 ** |
| DEL | 1000 | 4000 | 4000 | 0 | 0 | 1.0000 |
| DUP | 3 | 0 | 0 | 0 | 0 | 0.0000 |
| INS | 30 | 300 | 300 | 3 | 3 | 0.9900 |
| INV | 7 | 70 | 70 | 1 | 1 | 0.9900 ** |
| Unknown | 1 | 2 | 3 | 1 | 2 | 0.0000 ** |

"##;

const EVERY_LAYER_CHARGED_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
b_second	DEL	only_perl	21:1010-2010	deletion	ENST00000000010	Transcript	feature_truncation,intron_variant
b_second	Unknown	only_rust	21:1001-2001	deletion	ENST00000000001	Transcript	feature_truncation,intron_variant
a_first	BND	only_rust	21:1002-2002	deletion	ENST00000000002	Transcript	feature_truncation,intron_variant
"##;

const EVERY_LAYER_CHARGED_ADJUSTED_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
b_second	DEL	only_perl	21:1010-2010	deletion	ENST00000000010	Transcript	feature_truncation,intron_variant
"##;

const PUBLISHED_SHAPE_JSON: &str = r##"{
  "overall": {
    "total_input": 7,
    "perl_variants": 7,
    "rust_variants": 6,
    "perl_tuples": 59,
    "rust_tuples": 60,
    "intersection": 56,
    "only_perl": 3,
    "only_rust": 4,
    "precision": 0.933333,
    "recall": 0.949153,
    "f1": 0.941176,
    "perl_tuples_union": 59,
    "rust_tuples_union": 60,
    "intersection_union": 56,
    "cross_file_duplicate_perl": 0,
    "cross_file_duplicate_rust": 0
  },
  "per_file": {
    "one": {
      "total_input": 5,
      "perl_variants": 5,
      "rust_variants": 5,
      "perl_tuples": 50,
      "rust_tuples": 52,
      "intersection": 48,
      "only_perl": 2,
      "only_rust": 4,
      "precision": 0.923077,
      "recall": 0.96,
      "f1": 0.941176
    },
    "two": {
      "total_input": 2,
      "perl_variants": 2,
      "rust_variants": 1,
      "perl_tuples": 9,
      "rust_tuples": 8,
      "intersection": 8,
      "only_perl": 1,
      "only_rust": 0,
      "precision": 1.0,
      "recall": 0.888889,
      "f1": 0.941176
    }
  },
  "per_variant_type": {
    "DEL": {
      "total_input": 7,
      "perl_variants": 7,
      "rust_variants": 6,
      "perl_tuples": 59,
      "rust_tuples": 60,
      "intersection": 56,
      "only_perl": 3,
      "only_rust": 4,
      "precision": 0.933333,
      "recall": 0.949153,
      "f1": 0.941176
    }
  },
  "adjusted": {
    "total_input": 7,
    "perl_variants": 7,
    "rust_variants": 6,
    "perl_tuples": 58,
    "rust_tuples": 57,
    "intersection": 55,
    "only_perl": 3,
    "only_rust": 2,
    "precision": 0.964912,
    "recall": 0.948276,
    "f1": 0.956522,
    "perl_tuples_union": 58,
    "rust_tuples_union": 57,
    "intersection_union": 55,
    "cross_file_duplicate_perl": 0,
    "cross_file_duplicate_rust": 0,
    "excluded_divergences": 4,
    "excluded_rust_extra": 4,
    "excluded_perl_extra": 3,
    "excluded_cnvtr_rust": 2,
    "excluded_cnvtr_perl": 2,
    "excluded_transcript_selection_rust": 2,
    "excluded_transcript_selection_perl": 1,
    "transcript_selection_cnvtr_overlap_rust": 0,
    "cnvtr_swap_pairs_total_rust": 0,
    "cnvtr_swap_pairs_total_perl": 0,
    "excluded_cross_chromosome_perl": 0,
    "excluded_cross_chromosome_rust": 0,
    "excluded_transcript_selection_mate_rust": 0,
    "excluded_registry_swap_rust": 0,
    "excluded_registry_swap_perl": 0,
    "excluded_registry_swap_by_bucket": {
      "start_cooccurrence_swap": 0,
      "covered_splice_region_swap": 0
    },
    "excluded_non_ref_batch_rust": 0,
    "excluded_non_ref_batch_perl": 0,
    "excluded_breakend_mate_local_read_rust": 0,
    "excluded_breakend_mate_local_read_perl": 0,
    "excluded_breakend_mate_local_read_by_kind": {
      "cross_chromosome": 0,
      "same_chromosome": 0
    },
    "transcript_selection_cross_chromosome_overlap_perl": 0,
    "cross_chromosome_check_ran": false,
    "scored_engine": "vep-rs",
    "scored_engine_transcripts_absent_from_cache": 0,
    "excluded_cross_chromosome_by_file": {},
    "excluded_cnvtr_by_file": {
      "one": 2,
      "two": 1
    },
    "excluded_transcript_selection_by_file": {
      "one": {
        "rust": 2,
        "perl": 1
      }
    },
    "cnvtr_swap_pairs_total_by_file": {},
    "excluded_cross_chromosome_rust_by_file": {},
    "excluded_transcript_selection_mate_by_file": {},
    "excluded_registry_swap_by_file": {},
    "excluded_non_ref_batch_by_file": {},
    "excluded_breakend_mate_local_read_by_file": {}
  },
  "per_file_adjusted": {
    "one": {
      "total_input": 5,
      "perl_variants": 5,
      "rust_variants": 5,
      "perl_tuples": 49,
      "rust_tuples": 50,
      "intersection": 48,
      "only_perl": 1,
      "only_rust": 2,
      "precision": 0.96,
      "recall": 0.979592,
      "f1": 0.969697
    },
    "two": {
      "total_input": 2,
      "perl_variants": 2,
      "rust_variants": 1,
      "perl_tuples": 9,
      "rust_tuples": 7,
      "intersection": 7,
      "only_perl": 2,
      "only_rust": 0,
      "precision": 1.0,
      "recall": 0.777778,
      "f1": 0.875
    }
  }
}
"##;

const PUBLISHED_SHAPE_MD: &str = r##"# SV Concordance Report

## Overall

- **F1 (adjusted)**: 0.9565  *(excluded 4 Rust-extra + 3 Perl-extra transcript divergences)*
- **F1 (raw)**: 0.9412
- **Total input variants**: 7
- **Perl tuples**: 59
- **Rust tuples**: 60
- **Intersection**: 56
- **Only Perl**: 3
- **Only Rust**: 4

SV adjusted F1 additionally excludes the `<CNV:TR>` literal-expansion class (2 pairs, one, two only). Raw F1 unaffected.

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|
| one | 5 | 50 | 52 | 48 | 0.9412 | 0.9697 |
| two | 2 | 9 | 8 | 8 | 0.9412 | 0.8750 |

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|
| DEL | 7 | 59 | 60 | 3 | 4 | 0.9412 ** |

"##;

const PUBLISHED_SHAPE_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
one	DEL	only_perl	21:1010-2010	deletion	ENST00000000010	Transcript	feature_truncation,intron_variant
one	DEL	only_rust	21:1007-2007	deletion	ENST00000000007	Transcript	feature_truncation,intron_variant
"##;

const PUBLISHED_SHAPE_ADJUSTED_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
one	DEL	only_rust	21:1007-2007	deletion	ENST00000000007	Transcript	feature_truncation,intron_variant
"##;

const NOTHING_SET_ASIDE_JSON: &str = r##"{
  "overall": {
    "total_input": 3,
    "perl_variants": 3,
    "rust_variants": 3,
    "perl_tuples": 12,
    "rust_tuples": 12,
    "intersection": 12,
    "only_perl": 0,
    "only_rust": 0,
    "precision": 1.0,
    "recall": 1.0,
    "f1": 1.0,
    "perl_tuples_union": 12,
    "rust_tuples_union": 12,
    "intersection_union": 12,
    "cross_file_duplicate_perl": 0,
    "cross_file_duplicate_rust": 0
  },
  "per_file": {
    "x": {
      "total_input": 3,
      "perl_variants": 3,
      "rust_variants": 3,
      "perl_tuples": 12,
      "rust_tuples": 12,
      "intersection": 12,
      "only_perl": 0,
      "only_rust": 0,
      "precision": 1.0,
      "recall": 1.0,
      "f1": 1.0
    }
  },
  "per_variant_type": {
    "SNV": {
      "total_input": 3,
      "perl_variants": 3,
      "rust_variants": 3,
      "perl_tuples": 12,
      "rust_tuples": 12,
      "intersection": 12,
      "only_perl": 0,
      "only_rust": 0,
      "precision": 1.0,
      "recall": 1.0,
      "f1": 1.0
    }
  },
  "adjusted": {
    "total_input": 3,
    "perl_variants": 3,
    "rust_variants": 3,
    "perl_tuples": 12,
    "rust_tuples": 12,
    "intersection": 12,
    "only_perl": 0,
    "only_rust": 0,
    "precision": 1.0,
    "recall": 1.0,
    "f1": 1.0,
    "perl_tuples_union": 12,
    "rust_tuples_union": 12,
    "intersection_union": 12,
    "cross_file_duplicate_perl": 0,
    "cross_file_duplicate_rust": 0,
    "excluded_divergences": 0,
    "excluded_rust_extra": 0,
    "excluded_perl_extra": 0,
    "excluded_cnvtr_rust": 0,
    "excluded_cnvtr_perl": 0,
    "excluded_transcript_selection_rust": 0,
    "excluded_transcript_selection_perl": 0,
    "transcript_selection_cnvtr_overlap_rust": 0,
    "cnvtr_swap_pairs_total_rust": 0,
    "cnvtr_swap_pairs_total_perl": 0,
    "excluded_cross_chromosome_perl": 0,
    "excluded_cross_chromosome_rust": 0,
    "excluded_transcript_selection_mate_rust": 0,
    "excluded_registry_swap_rust": 0,
    "excluded_registry_swap_perl": 0,
    "excluded_registry_swap_by_bucket": {
      "start_cooccurrence_swap": 0,
      "covered_splice_region_swap": 0
    },
    "excluded_non_ref_batch_rust": 0,
    "excluded_non_ref_batch_perl": 0,
    "excluded_breakend_mate_local_read_rust": 0,
    "excluded_breakend_mate_local_read_perl": 0,
    "excluded_breakend_mate_local_read_by_kind": {
      "cross_chromosome": 0,
      "same_chromosome": 0
    },
    "transcript_selection_cross_chromosome_overlap_perl": 0,
    "cross_chromosome_check_ran": true,
    "scored_engine": "vep-rs",
    "scored_engine_transcripts_absent_from_cache": 0,
    "excluded_cross_chromosome_by_file": {},
    "excluded_cnvtr_by_file": {},
    "excluded_transcript_selection_by_file": {},
    "cnvtr_swap_pairs_total_by_file": {},
    "excluded_cross_chromosome_rust_by_file": {},
    "excluded_transcript_selection_mate_by_file": {},
    "excluded_registry_swap_by_file": {},
    "excluded_non_ref_batch_by_file": {},
    "excluded_breakend_mate_local_read_by_file": {}
  },
  "per_file_adjusted": {
    "x": {
      "total_input": 3,
      "perl_variants": 3,
      "rust_variants": 3,
      "perl_tuples": 12,
      "rust_tuples": 12,
      "intersection": 12,
      "only_perl": 0,
      "only_rust": 0,
      "precision": 1.0,
      "recall": 1.0,
      "f1": 1.0
    }
  }
}
"##;

const NOTHING_SET_ASIDE_MD: &str = r##"# SV Concordance Report

## Overall

- **F1 (raw)**: 1.0000
- **Total input variants**: 3
- **Perl tuples**: 12
- **Rust tuples**: 12
- **Intersection**: 12
- **Only Perl**: 0
- **Only Rust**: 0

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|
| x | 3 | 12 | 12 | 12 | 1.0000 |  |

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|
| SNV | 3 | 12 | 12 | 0 | 0 | 1.0000 |

"##;

const NOTHING_SET_ASIDE_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
"##;

const NOTHING_SET_ASIDE_ADJUSTED_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
"##;

const NO_FILE_JSON: &str = r##"{
  "overall": {
    "total_input": 0,
    "perl_variants": 0,
    "rust_variants": 0,
    "perl_tuples": 0,
    "rust_tuples": 0,
    "intersection": 0,
    "only_perl": 0,
    "only_rust": 0,
    "precision": 0.0,
    "recall": 0.0,
    "f1": 0.0
  },
  "per_file": {},
  "per_variant_type": {},
  "adjusted": {
    "total_input": 0,
    "perl_variants": 0,
    "rust_variants": 0,
    "perl_tuples": 0,
    "rust_tuples": 0,
    "intersection": 0,
    "only_perl": 0,
    "only_rust": 0,
    "precision": 0.0,
    "recall": 0.0,
    "f1": 0.0,
    "excluded_divergences": 0,
    "excluded_rust_extra": 0,
    "excluded_perl_extra": 0,
    "excluded_cnvtr_rust": 0,
    "excluded_cnvtr_perl": 0,
    "excluded_transcript_selection_rust": 0,
    "excluded_transcript_selection_perl": 0,
    "transcript_selection_cnvtr_overlap_rust": 0,
    "cnvtr_swap_pairs_total_rust": 0,
    "cnvtr_swap_pairs_total_perl": 0,
    "excluded_cross_chromosome_perl": 0,
    "excluded_cross_chromosome_rust": 0,
    "excluded_transcript_selection_mate_rust": 0,
    "excluded_registry_swap_rust": 0,
    "excluded_registry_swap_perl": 0,
    "excluded_registry_swap_by_bucket": {},
    "excluded_non_ref_batch_rust": 0,
    "excluded_non_ref_batch_perl": 0,
    "excluded_breakend_mate_local_read_rust": 0,
    "excluded_breakend_mate_local_read_perl": 0,
    "excluded_breakend_mate_local_read_by_kind": {},
    "transcript_selection_cross_chromosome_overlap_perl": 0,
    "cross_chromosome_check_ran": false,
    "scored_engine": "vep-rs",
    "scored_engine_transcripts_absent_from_cache": 0,
    "excluded_cnvtr_by_file": {},
    "excluded_transcript_selection_by_file": {},
    "cnvtr_swap_pairs_total_by_file": {}
  },
  "per_file_adjusted": {}
}
"##;

const NO_FILE_MD: &str = r##"# SV Concordance Report

## Overall

- **F1 (raw)**: 0.0000
- **Total input variants**: 0
- **Perl tuples**: 0
- **Rust tuples**: 0
- **Intersection**: 0
- **Only Perl**: 0
- **Only Rust**: 0

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|

"##;

const NO_FILE_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
"##;

const NO_FILE_ADJUSTED_TSV: &str = r##"source_file	variant_type	direction	location	allele	feature	feature_type	consequence
"##;
