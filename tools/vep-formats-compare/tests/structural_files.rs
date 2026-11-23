// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The files the structural-variant classes read: the records of the input VCF and the
//! shards of the candidate JSON cache, on planted bytes. Every expected value is stated for
//! those bytes.

mod compare_common;

use std::io::Write;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use std::path::PathBuf;

use compare_common::scratch;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Uncaught;
use vep_formats_compare::rules::structural::iter_cache_transcripts;
use vep_formats_compare::rules::structural::load_cache_authority;
use vep_formats_compare::rules::structural::parse_vcf;
use vep_formats_compare::rules::structural::StableId;
use vep_formats_compare::rules::structural::StructuralTuple;
use vep_formats_compare::rules::structural::StructuralVariants;

fn gzipped(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(bytes).expect("the bytes deflate");
    encoder.finish().expect("the member ends")
}

fn plant(root: &Path, name: &str, bytes: &[u8]) -> PathBuf {
    let path = root.join(name);
    std::fs::create_dir_all(path.parent().expect("a directory")).expect("the directory is made");
    std::fs::write(&path, bytes).expect("the file is written");
    path
}

/// What a call says that ends on an input no reading of the comparison expects.
fn uncaught(call: impl FnOnce()) -> String {
    static HOOK: std::sync::Once = std::sync::Once::new();
    HOOK.call_once(|| {
        let default_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if info.payload().downcast_ref::<Uncaught>().is_none() {
                default_hook(info);
            }
        }));
    });
    let payload = catch_unwind(AssertUnwindSafe(call)).expect_err("the call ends the run");
    payload
        .downcast_ref::<Uncaught>()
        .expect("the payload of such an ending")
        .0
        .clone()
}

const INPUT: &str = "##x\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     21\t100\tid1\tA\tG\t.\t.\tSVTYPE=DEL\r\n\
                     chr21\t 30425949\t.\tN\t<DEL>\r\
                     21\t3e7\tx\tACG\tA\t.\t.\n\
                     21\t5\n\
                     #late\n\
                     M\t.\tm\tN\tN.\t.\t.\tEND=7\tGT";

/// The records of `INPUT`: chromosome, position, ID, REF, ALT, type, Location, INFO and
/// the POS column as written.
const RECORDS: [[&str; 9]; 4] = [
    [
        "21",
        "100",
        "id1",
        "A",
        "G",
        "SNV",
        "21:100",
        "SVTYPE=DEL",
        "100",
    ],
    [
        "chr21",
        "30425949",
        ".",
        "N",
        "<DEL>",
        "Symbolic_DEL",
        "21:30425949",
        "",
        " 30425949",
    ],
    [
        "21",
        "30000000",
        "x",
        "ACG",
        "A",
        "Small_DEL",
        "21:30000000",
        "",
        "3e7",
    ],
    ["M", "0", "m", "N", "N.", "BND", "MT:0", "END=7", "."],
];

fn read(path: &Path) -> Vec<[String; 9]> {
    parse_vcf(path)
        .expect("the input opens")
        .into_iter()
        .map(|v| {
            [
                v.chrom,
                v.pos.to_string(),
                v.vid,
                v.r#ref,
                v.alt,
                v.variant_type,
                v.location,
                v.info,
                v.pos_text.expect("the POS column"),
            ]
        })
        .collect()
}

fn records() -> Vec<[String; 9]> {
    RECORDS
        .iter()
        .map(|record| record.map(str::to_string))
        .collect()
}

/// A line ends at a line feed, a carriage return or the two together; a line of fewer than
/// five columns and a `#` line anywhere are no record; POS is numified and INFO is the
/// eighth column where there is one.
#[test]
fn an_input_is_read_as_its_records() {
    let root = scratch("structural-files-input");
    assert_eq!(read(&plant(&root, "in.vcf", INPUT.as_bytes())), records());
}

/// A name whose suffix is `.gz` is read as `zcat` writes it, one member after another; a
/// name that is `.gz` alone or ends in `.GZ` is read as it stands; a gzip name that is
/// absent or is not gzip has no record, and any other name that cannot be opened is an
/// error with the system's words.
#[test]
fn an_input_is_inflated_by_its_name() {
    let root = scratch("structural-files-gzip");
    let raw = INPUT.as_bytes();
    let mut members = gzipped(&raw[..40]);
    members.extend(gzipped(&raw[40..]));
    members.extend(gzipped(b""));
    assert_eq!(read(&plant(&root, "in.vcf.gz", &members)), records());
    assert_eq!(read(&plant(&root, "plain/.gz", raw)), records());
    assert_eq!(read(&plant(&root, "x.vcf.GZ", raw)), records());
    assert!(read(&root.join("absent.vcf.gz")).is_empty());
    assert!(read(&plant(&root, "text.vcf.gz", raw)).is_empty());
    let words = |path: &Path| {
        let error = parse_vcf(path).expect_err("the input does not open");
        vep_formats_compare::report::strerror(&error)
    };
    assert_eq!(words(&root.join("absent.vcf")), "No such file or directory");
    assert_eq!(words(&root), "Is a directory");
    assert_eq!(
        uncaught(|| {
            let _ = parse_vcf(&plant(&root, "bad.vcf", b"21\t100\t.\tA\tG\n\xff\n"));
        }),
        "'utf-8' codec can't decode byte 0xff in position 13: invalid start byte"
    );
}

fn ids(ids: &[&str]) -> Vec<StableId> {
    ids.iter().map(|id| StableId::text(id)).collect()
}

/// The shards of a chromosome are read in the order of their names, plain and gzipped
/// together, so a transcript filed twice takes its span from the last; a transcript
/// without a start and an end has no span, one without a strand no model; a `null`, a
/// record without a `stable_id` and a range-keyed map's other values are passed over.
#[test]
fn a_cache_is_read_shard_by_shard_in_the_order_of_their_names() {
    let root = scratch("structural-files-cache");
    let shard = |name: &str, text: &str| {
        plant(&root, &format!("transcripts/21/{name}"), text.as_bytes());
    };
    shard(
        "s.json",
        r#"[{"stable_id": "T", "start": 30, "end": 40, "strand": 1}]"#,
    );
    shard(
        "t.json",
        r#"[{"stable_id": "T", "start": 50, "end": 60, "strand": 1}, null, {"start": 1, "end": 2}, {"stable_id": ""}, {"stable_id": "U", "start": 5}]"#,
    );
    plant(
        &root,
        "transcripts/21/a.json.gz",
        &gzipped(br#"[{"stable_id": "T", "start": 70, "end": 80, "strand": 1}]"#),
    );
    plant(
        &root,
        "transcripts/21/s.json.gz",
        &gzipped(
            br#"{"a": [{"stable_id": "T", "start": 1, "end": 2, "strand": 1}], "b": null, "c": "text", "d": {"k": 1}, "e": [null, {"stable_id": "V", "start": "7", "end": 9.9, "strand": true}]}"#,
        ),
    );
    plant(&root, "transcripts/a file", b"[]");
    plant(
        &root,
        "transcripts/X/1.json",
        br#"[{"stable_id": "T", "start": 1, "end": 2}]"#,
    );
    let (by_chr, spans, models) = load_cache_authority(&root);
    let chromosomes: Vec<(&str, Vec<StableId>)> = by_chr
        .iter()
        .map(|(chrom, held)| (chrom.as_str(), held.iter().cloned().collect()))
        .collect();
    assert_eq!(
        chromosomes,
        [("21", ids(&["T", "U", "V"])), ("X", ids(&["T"]))]
    );
    let span = |chrom: &str, id: &str| spans.get(&(chrom.to_string(), StableId::text(id))).cloned();
    assert_eq!(span("21", "T"), Some(("21".to_string(), 50, 60)));
    assert_eq!(span("21", "U"), None);
    assert_eq!(span("21", "V"), Some(("21".to_string(), 7, 9)));
    assert_eq!(span("X", "T"), Some(("X".to_string(), 1, 2)));
    assert_eq!(spans.len(), 3);
    let model = |chrom: &str, id: &str| models.get(&(chrom.to_string(), StableId::text(id)));
    assert_eq!(
        model("21", "T").map(|model| (model.start, model.end)),
        Some((50, 60))
    );
    assert_eq!(model("21", "V").map(|model| model.strand), Some(1));
    assert!(model("X", "T").is_none());
    let yielded: Vec<String> = iter_cache_transcripts(&root)
        .into_iter()
        .map(|(chrom, tx)| {
            format!(
                "{chrom} {} {}",
                tx.get("stable_id").and_then(|id| id.as_str()).unwrap_or(""),
                tx.get("start")
                    .map_or(String::new(), |start| format!("{start:?}"))
            )
        })
        .collect();
    assert_eq!(
        yielded,
        [
            "21 T Int(70)",
            "21 T Int(30)",
            "21 T Int(1)",
            "21 V Str(\"7\")",
            "21 T Int(50)",
            "21 U Int(5)",
            "X T Int(1)"
        ]
    );
    assert!(iter_cache_transcripts(&root.join("absent")).is_empty());
}

/// A shard that cannot be read as gzip, as UTF-8 or as JSON is passed over; a deflate
/// stream that is refused, a `stable_id` that is a list and an infinite coordinate end the
/// run, each with the text stated for it.
#[test]
fn a_shard_that_cannot_be_read_is_passed_over_and_a_refused_one_ends_the_run() {
    let root = scratch("structural-files-shards");
    let good = br#"[{"stable_id": "T", "start": 5, "end": 9, "strand": 1}]"#;
    let cache = |name: &str, shard: &str, bytes: &[u8]| -> PathBuf {
        let cache = root.join(name);
        plant(&cache, &format!("transcripts/21/{shard}"), bytes);
        cache
    };
    let spans = |cache: &Path| load_cache_authority(cache).1.len();
    let member = gzipped(good);
    assert_eq!(spans(&cache("good", "s.json.gz", &member)), 1);
    let mut padded = member.clone();
    padded.extend([0, 0, 0]);
    assert_eq!(spans(&cache("padded", "s.json.gz", &padded)), 1);
    let mut trailed = member.clone();
    trailed.extend(b"trailing");
    assert_eq!(spans(&cache("trailed", "s.json.gz", &trailed)), 0);
    let mut wrong_sum = member.clone();
    let at = wrong_sum.len() - 8;
    wrong_sum[at] ^= 1;
    assert_eq!(spans(&cache("sum", "s.json.gz", &wrong_sum)), 0);
    assert_eq!(
        spans(&cache("cut", "s.json.gz", &member[..member.len() - 6])),
        0
    );
    assert_eq!(spans(&cache("plain", "s.json.gz", good)), 0);
    assert_eq!(spans(&cache("bytes", "s.json", b"\xff[]")), 0);
    assert_eq!(spans(&cache("mark", "s.json", b"\xef\xbb\xbf[]")), 0);
    assert_eq!(spans(&cache("half", "s.json", &good[..20])), 0);
    assert_eq!(spans(&cache("number", "s.json", b"5")), 0);
    std::fs::create_dir_all(root.join("listed/transcripts/21/d.json")).expect("a directory");
    assert_eq!(spans(&root.join("listed")), 0);

    let mut refused = member[..10].to_vec();
    refused.extend([0xff; 8]);
    refused.extend([0; 8]);
    let refused = cache("refused", "s.json.gz", &refused);
    assert_eq!(
        uncaught(|| {
            load_cache_authority(&refused);
        }),
        "Error -3 while decompressing data: invalid block type"
    );
    let listed = cache(
        "list-id",
        "s.json",
        br#"[{"stable_id": [1], "start": 1, "end": 2}]"#,
    );
    assert_eq!(
        uncaught(|| {
            load_cache_authority(&listed);
        }),
        "cannot use 'list' as a set element (unhashable type: 'list')"
    );
    let infinite = cache(
        "infinite",
        "s.json",
        br#"[{"stable_id": "a", "start": "x", "end": -Infinity, "strand": 1}]"#,
    );
    assert_eq!(
        uncaught(|| {
            load_cache_authority(&infinite);
        }),
        "cannot convert float infinity to integer"
    );
    let walked = cache("walked", "s.json", br#"{"x": 2.5}"#);
    assert_eq!(
        uncaught(|| {
            load_cache_authority(&walked);
        }),
        "'float' object is not iterable"
    );
}

/// A `stable_id` that is a number is a key by its value: `1`, `1.0` and `true` are one
/// transcript, every NaN another, and each counts among the cache's transcripts.
#[test]
fn a_stable_id_that_is_a_number_is_a_key_by_its_value() {
    let root = scratch("structural-files-ids");
    plant(
        &root,
        "transcripts/21/s.json",
        br#"[{"stable_id": NaN, "start": 1, "end": 2}, {"stable_id": NaN, "start": 3, "end": 4},
             {"stable_id": 1, "start": 5, "end": 6}, {"stable_id": 1.0, "start": 7, "end": 8},
             {"stable_id": true, "start": 9, "end": 10}, {"stable_id": 2.5, "start": 1, "end": 1},
             {"stable_id": 1e400, "start": 1, "end": 1}, {"stable_id": Infinity, "start": 1, "end": 1},
             {"stable_id": 123456789012345678901234567890, "start": 1, "end": 1},
             {"stable_id": 1.2345678901234568e29, "start": 1, "end": 1},
             {"stable_id": -7, "start": 1, "end": 1}, {"stable_id": -7.0, "start": 2, "end": 2}]"#,
    );
    let (by_chr, spans, models) = load_cache_authority(&root);
    assert_eq!(by_chr["21"].len(), 7);
    assert_eq!(spans.len(), 7);
    assert!(models.is_empty());
    let span = |id: StableId| {
        spans
            .get(&("21".to_string(), id))
            .map(|span| (span.1, span.2))
    };
    assert_eq!(span(StableId::Nan), Some((3, 4)));
    assert_eq!(span(StableId::Whole("1".to_string())), Some((9, 10)));
    assert_eq!(span(StableId::Whole("-7".to_string())), Some((2, 2)));
    assert_eq!(
        span(StableId::Fraction(f64::INFINITY.to_bits())),
        Some((1, 1))
    );
    assert_eq!(span(StableId::text("1")), None);
    let table = StructuralVariants::new(&["bnd_own_chromosome_rows_lost"], None, Some(&root))
        .expect("the side table reads");
    assert_eq!(table.spans.len(), 7);
}

/// An entry of the VCF or JSON format takes its Location from its input record, whatever
/// allele it names where the record's rows share one Location; an entry of a record past
/// the input has none, a record counted from the end is read as an index below zero is,
/// and a record entry has no tuple.
#[test]
fn an_entry_without_a_location_takes_the_one_of_its_record() {
    let root = scratch("structural-files-entries");
    let input = plant(
        &root,
        "in.vcf",
        b"21\t100\tdel1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=300\n21\t200\tsnv\tA\tG,T\n21\t300\tins\tA\tAT,G\n21\t400\tref\tA\t<NON_REF>\n",
    );
    let table = StructuralVariants::new(&["reference_skipped_record"], Some(&input), None)
        .expect("the side table reads");
    let entry = |key: &str, record: Option<i64>, location: Option<&str>, allele: &str| Entry {
        key: key.to_string(),
        record,
        location: location.map(str::to_string),
        allele: allele.to_string(),
        feature: "ENST1".to_string(),
        consequence: "b&a,a".to_string(),
        ..Entry::default()
    };
    let location =
        |record: Option<i64>, allele: &str| table.location(&entry("k", record, None, allele));
    assert_eq!(location(Some(0), "deletion").as_deref(), Some("21:101-300"));
    assert_eq!(location(Some(0), "other").as_deref(), Some("21:101-300"));
    assert_eq!(location(Some(1), "G").as_deref(), Some("21:200"));
    assert_eq!(location(Some(2), "T").as_deref(), Some("21:300"));
    assert_eq!(location(Some(2), "AT").as_deref(), Some("21:300"));
    assert_eq!(location(Some(3), "G").as_deref(), Some("21:401"));
    assert_eq!(location(Some(4), "G"), None);
    assert_eq!(location(None, "G"), None);
    assert_eq!(location(Some(-1), "T").as_deref(), Some("21:401"));
    assert_eq!(location(Some(-2), "T").as_deref(), Some("21:300"));
    assert_eq!(location(Some(-4), "x").as_deref(), Some("21:101-300"));
    assert_eq!(
        uncaught(|| {
            location(Some(-5), "x");
        }),
        "list index out of range"
    );
    assert_eq!(
        table
            .location(&entry("k", Some(7), Some("chrM:0100-+0200"), "G"))
            .as_deref(),
        Some("MT:100-200")
    );
    let tuples = table.tuples(&[
        entry("rec0|record", Some(0), None, ""),
        entry(
            "rec0|deletion|ENST1|Transcript#2",
            Some(0),
            None,
            "deletion",
        ),
        entry("rec3|G|ENST1|Transcript", Some(3), None, "G"),
        entry("rec9|G|ENST1|Transcript", Some(9), None, "G"),
        entry(
            "21:5|N[chr22:7[|ENST1|-",
            None,
            Some("chr21:5"),
            "N[chr22:7[",
        ),
        entry("k|x", Some(2), None, ""),
    ]);
    let rows: Vec<(&str, &StructuralTuple)> = tuples
        .iter()
        .map(|(key, tuple)| (key.as_str(), tuple))
        .collect();
    let tuple = |location: &str, allele: &str, kind: &str| {
        (
            location.to_string(),
            allele.to_string(),
            "ENST1".to_string(),
            kind.to_string(),
            "a,b".to_string(),
        )
    };
    assert_eq!(
        rows,
        [
            (
                "rec0|deletion|ENST1|Transcript#2",
                &tuple("21:101-300", "deletion", "Transcript")
            ),
            (
                "rec3|G|ENST1|Transcript",
                &tuple("21:401", "G", "Transcript")
            ),
            ("21:5|N[chr22:7[|ENST1|-", &tuple("21:5", "N[22:7[", "-")),
            ("k|x", &tuple("21:300", "-", "x")),
        ]
    );
}
