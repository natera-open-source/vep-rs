// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the records of a run: the lines of a flag set's own input against the
//! manifest's records, the region lines that name their records, the records the
//! reference's JSON carries without a row, and the records the reference skipped.

use super::*;
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;

fn object(text: &str) -> Object {
    match loads(text) {
        Ok(Value::Object(object)) => object,
        other => panic!("no JSON object: {other:?}"),
    }
}

fn scratch(tag: &str) -> TempDir {
    TempDir::create(&std::env::temp_dir(), &format!("golden_records_{tag}_")).unwrap()
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

/// A manifest of three records: the second warned, the third a twin of the first whose
/// warning is no text.
const MANIFEST: &str = concat!(
    r#"{"records": [{"dataset": "s01", "vcf": "21\t100\trs1\tA\tG\t.\tPASS\t."}, "#,
    r#"{"dataset": "s01", "vcf": "21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900", "#,
    r#""reference_warning": "CPX is not supported"}, "#,
    r#"{"dataset": "s01", "vcf": "21\t100\trs1\tA\tG\t.\tPASS\tAC=1", "reference_warning": 5}]}"#
);

/// The lines and manifest places of the records of the flag set whose input is `name`,
/// or the line of the refusal with the corpus written `{C}`.
fn records_of(
    corpus: &TempDir,
    name: &str,
    content: &[u8],
) -> Result<Vec<(String, Option<usize>)>, String> {
    std::fs::write(corpus.path.join(name), content).unwrap();
    let directory = corpus.path.to_str().unwrap();
    let target = object(&format!(r#"{{"input": "{name}"}}"#));
    match run_records(directory, &object(MANIFEST), &target) {
        Ok(records) => Ok(records
            .into_iter()
            .map(|record| (record.vcf, record.manifest_index))
            .collect()),
        Err(refusal) => Err(refusal.line.replace(directory, "{C}")),
    }
}

fn lines(records: &[(&str, Option<usize>)]) -> Result<Vec<(String, Option<usize>)>, String> {
    Ok(records
        .iter()
        .map(|(line, place)| (line.to_string(), *place))
        .collect())
}

fn refusal(message: &str) -> Result<Vec<(String, Option<usize>)>, String> {
    Err(format!("ERROR: [build_golden_corpus] {message}"))
}

const FIRST: &str = "21\t100\trs1\tA\tG\t.\tPASS\t.";
const WARNED: &str = "21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900";
const TWIN: &str = "21\t100\trs1\tA\tG\t.\tPASS\tAC=1";

#[test]
fn the_first_five_columns_end_at_the_fifth_tab() {
    assert_eq!(first_five_columns("a\tb\tc\td\te\tf\tg"), "a\tb\tc\td\te");
    assert_eq!(first_five_columns("a\tb"), "a\tb");
    assert_eq!(first_five_columns(""), "");
    assert_eq!(first_five_columns("a\tb\tc\td\te"), "a\tb\tc\td\te");
    assert_eq!(first_five_columns("a\tb\tc\td\te\t"), "a\tb\tc\td\te");
}

#[test]
fn a_run_over_the_corpus_input_takes_the_manifests_records() {
    let manifest = object(MANIFEST);
    let expected = [
        RunRecord {
            vcf: FIRST.to_string(),
            manifest_index: Some(0),
            reference_warning: None,
        },
        RunRecord {
            vcf: WARNED.to_string(),
            manifest_index: Some(1),
            reference_warning: Some("CPX is not supported".to_string()),
        },
        RunRecord {
            vcf: TWIN.to_string(),
            manifest_index: Some(2),
            reference_warning: None,
        },
    ];
    for target in ["{}", r#"{"input": ""}"#, r#"{"input": "variants.vcf"}"#] {
        let records = run_records("absent", &manifest, &object(target)).unwrap();
        assert_eq!(records, expected, "{target}");
    }
    // What the manifest lacks is refused in the tool's own words.
    for (manifest, reason) in [
        ("{}", "the manifest holds no list of records"),
        (
            r#"{"records": {}}"#,
            "the manifest holds no list of records",
        ),
        (
            r#"{"records": [{"vcf": 5}]}"#,
            "manifest record 0 holds no vcf line",
        ),
        (
            r#"{"records": [{"vcf": "a"}, []]}"#,
            "manifest record 1 holds no vcf line",
        ),
    ] {
        let refused = run_records("absent", &object(manifest), &object("{}")).unwrap_err();
        assert_eq!(refused.line, format!("ERROR: [vep-golden-corpus] {reason}"));
    }
}

#[test]
fn the_lines_of_a_flag_sets_own_input_are_twins_records_of_their_own_or_the_manifests() {
    let corpus = scratch("input");
    // A VCF: a line that shares its first five columns with manifest records is the last
    // of them, another record is its own, in the file's order.
    assert_eq!(
        records_of(
            &corpus,
            "twin.vcf",
            b"#CHROM\tPOS\n21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n21\t7000\tnew\tC\tT\n21\t500\tsv1\tN\t<CPX>\n"
        ),
        lines(&[(TWIN, Some(2)), ("21\t7000\tnew\tC\tT", None), (WARNED, Some(1))])
    );
    // Every line end is a line end, and the last line needs none.
    assert_eq!(
        records_of(
            &corpus,
            "crlf.vcf",
            b"21\t100\trs1\tA\tG\r\n21\t500\tsv1\tN\t<CPX>\r21\t9\tx\tA\tC"
        ),
        lines(&[(TWIN, Some(2)), (WARNED, Some(1)), ("21\t9\tx\tA\tC", None)])
    );
    // A file without a VCF record is the manifest's records in another encoding.
    assert_eq!(
        records_of(
            &corpus,
            "space.vcf",
            b"#h\n21 100 rs1 A G\n21 500 sv1 N <CPX>\n21 1 x A C\n"
        ),
        lines(&[(FIRST, Some(0)), (WARNED, Some(1)), (TWIN, Some(2))])
    );
    assert_eq!(
        records_of(&corpus, "short.txt", b"21\t100\t100\tA/G\t1\trs1\n"),
        refusal(
            "{C}/short.txt holds no VCF record and its 1 records are not the manifest's 3 re-encoded"
        )
    );
    // In a VCF every line is a twin or a record: the refusal shows the line's first 40
    // characters.
    assert_eq!(
        records_of(
            &corpus,
            "mixed.vcf",
            b"21\t7000\tsnv_a\tC\tT\n21\t.\tsnv_x\tG\tA\twith a long tail that goes past forty characters\n"
        ),
        refusal(
            "{C}/mixed.vcf record 1 ('21\\t.\\tsnv_x\\tG\\tA\\twith a long tail that goe') is no VCF record in a VCF input"
        )
    );
    assert_eq!(
        records_of(&corpus, "blank.vcf", b"21\t7000\tsnv_a\tC\tT\n\n"),
        refusal("{C}/blank.vcf record 1 ('') is no VCF record in a VCF input")
    );
    // A POS is decimal digits of any script, and no other digit.
    assert_eq!(
        records_of(
            &corpus,
            "digits.vcf",
            "21\t\u{667}\tsnv_a\tC\tT\n".as_bytes()
        ),
        lines(&[("21\t\u{667}\tsnv_a\tC\tT", None)])
    );
    assert_eq!(
        records_of(
            &corpus,
            "power.vcf",
            "21\t\u{b2}\tsnv_a\tC\tT\n21\t1\tb\tA\tC\n21\t2\tb\tA\tC\n".as_bytes()
        ),
        refusal("{C}/power.vcf record 0 ('21\\t\u{b2}\\tsnv_a\\tC\\tT') is no VCF record in a VCF input")
    );
    // A REF is bases or a point.
    assert_eq!(
        records_of(
            &corpus,
            "ref.vcf",
            b"21\t5\tsnv_a\tAC/G\tT\n21\t5\tsnv_b\t.\tT\n21\t5\tsnv_c\t\tT\n"
        ),
        refusal(
            "{C}/ref.vcf record 0 ('21\\t5\\tsnv_a\\tAC/G\\tT') is no VCF record in a VCF input"
        )
    );
    assert_eq!(
        records_of(&corpus, "absent.vcf.gz", b"not gzip"),
        refusal("{C}/absent.vcf.gz is not a readable text file (Not a gzipped file (b'no'))")
    );
}

#[test]
fn a_region_line_is_one_column_of_a_place_and_an_allele() {
    // `chr:start-end`, a strand where one is written, a `/` or a `:` and the allele, bases
    // or a word of three characters or more, in either case (Utils/VEP.pm:391); the
    // position is the text before the allele's separator.
    for (line, position) in [
        ("21:100-100/G", "21:100-100"),
        ("21:100-100:G", "21:100-100"),
        ("21:100-100:1/G", "21:100-100:1"),
        ("21:100-100:+1/g", "21:100-100:+1"),
        ("21:100-100:-1:G", "21:100-100:-1"),
        ("chr21:101-100/ACGT", "chr21:101-100"),
        ("21:150-150/-", "21:150-150"),
        ("21:301-900/DEL", "21:301-900"),
        ("21:301-900/del", "21:301-900"),
        ("21:301-900:1/INS:ME:ALU", "21:301-900:1"),
        ("21:100-100:1DEL", "21:100-100"),
    ] {
        assert_eq!(region_position(line), Some(position), "{line}");
    }
    for line in [
        "",
        "rs1",
        "21:100/G",
        "21:100-100",
        "21:100-100/",
        "21:100-100:1",
        "21:100-100/XY",
        "21:100-100/G extra",
        "21:100-100/G\t.",
        "21\t100\t100\tA/G\t1\trs1",
        "21 100 rs1 A G",
        "21:g.100A>G",
        "ENST1:c.100A>G",
        // Ensembl VEP matches bytes, so a letter outside ASCII is none of the pattern's,
        // the two that fold to `s` and `k` included.
        "21:100-100/\u{17f}\u{17f}\u{17f}",
        "21:100-100/de\u{212a}",
    ] {
        assert_eq!(region_position(line), None, "{line}");
    }
}

/// The region lines of the flag set whose input is `name`, or the line of the refusal
/// with the corpus written `{C}`.
fn region_lines_of(
    corpus: &TempDir,
    name: &str,
    content: &[u8],
) -> Result<Vec<Option<String>>, String> {
    std::fs::write(corpus.path.join(name), content).unwrap();
    let directory = corpus.path.to_str().unwrap();
    let target = object(&format!(r#"{{"input": "{name}"}}"#));
    region_lines(directory, &target).map_err(|refusal| refusal.line.replace(directory, "{C}"))
}

#[test]
fn the_region_lines_of_a_run_stand_at_the_ordinals_of_their_records() {
    let corpus = scratch("region");
    let named = |lines: &[Option<&str>]| -> Result<Vec<Option<String>>, String> {
        Ok(lines.iter().map(|line| line.map(str::to_string)).collect())
    };
    // Every line that is no header has an ordinal; a region line stands at its own, as
    // written, and a line of another encoding leaves its ordinal to the VCF line.
    assert_eq!(
        region_lines_of(
            &corpus,
            "variants.region.txt",
            b"#h\n21:100-100:1/G\n21:g.150del\nchr21:301-900/del\n"
        ),
        named(&[Some("21:100-100:1/G"), None, Some("chr21:301-900/del")])
    );
    // A file that holds a VCF record names its records by their VCF lines, as a run over
    // `variants.vcf` does.
    assert_eq!(
        region_lines_of(&corpus, "own.vcf", b"21:100-100/G\n21\t100\trs1\tA\tG\n"),
        named(&[])
    );
    for target in ["{}", r#"{"input": ""}"#, r#"{"input": "variants.vcf"}"#] {
        assert_eq!(
            region_lines("absent", &object(target)).map_err(|refusal| refusal.line),
            named(&[]),
            "{target}"
        );
    }
    assert_eq!(
        region_lines_of(&corpus, "absent.txt.gz", b"not gzip"),
        Err(
            "ERROR: [build_golden_corpus] {C}/absent.txt.gz is not a readable text file (Not a gzipped file (b'no'))"
                .to_string()
        )
    );
}

#[test]
fn a_region_line_alone_names_its_record() {
    let record = |vcf: &str| RunRecord {
        vcf: vcf.to_string(),
        manifest_index: None,
        reference_warning: None,
    };
    let records = [
        record(FIRST),
        record("21\t149\t.\tGC\tG"),
        record("21\t300\tsv1\tN\t<DEL>"),
    ];
    let index = |region_lines: &[Option<&str>]| {
        let lines = region_lines
            .iter()
            .map(|line| line.map(str::to_string))
            .collect();
        let names = RunNames::of(&records, lines).unwrap();
        let mut index: Vec<(String, Vec<usize>)> = names.by_name.into_iter().collect();
        index.sort();
        index
    };
    let named = |names: &[(&str, &[usize])]| -> Vec<(String, Vec<usize>)> {
        names
            .iter()
            .map(|(name, ordinals)| (name.to_string(), ordinals.to_vec()))
            .collect()
    };
    // Without a region line every record carries the names of its VCF line.
    assert_eq!(
        index(&[]),
        named(&[
            ("21_149_GC/G", &[1]),
            ("21_150_GC/G", &[1]),
            ("rs1", &[0]),
            ("sv1", &[2])
        ])
    );
    // A region line is its record's one name; the record beside it keeps its own.
    assert_eq!(
        index(&[Some("21:100-100:1/G"), None, Some("21:301-900/DEL")]),
        named(&[
            ("21:100-100:1/G", &[0]),
            ("21:301-900/DEL", &[2]),
            ("21_149_GC/G", &[1]),
            ("21_150_GC/G", &[1])
        ])
    );
    // Two records of one line share its name, in order.
    assert_eq!(
        index(&[Some("21:100-100/G"), None, Some("21:100-100/G")]),
        named(&[
            ("21:100-100/G", &[0, 2]),
            ("21_149_GC/G", &[1]),
            ("21_150_GC/G", &[1])
        ])
    );
    let names = RunNames::of(&records, vec![Some("21:100-100/G".to_string())]).unwrap();
    assert_eq!(names.get("21:100-100/G"), Some(&vec![0]));
    assert_eq!(names.get("rs1"), None);
    assert_eq!(names.region_line(0), Some("21:100-100/G"));
    assert_eq!((names.region_line(1), names.region_line(9)), (None, None));
    // A VCF line that carries no name is refused, unless a region line names its record.
    let short = [record("21\t100")];
    let refused = RunNames::of(&short, Vec::new()).unwrap_err();
    assert_eq!(
        refused.line,
        "ERROR: [vep-golden-corpus] the record line '21\\t100' has fewer than five columns"
    );
    assert!(RunNames::of(&short, vec![Some("21:100-100/G".to_string())]).is_ok());
}

#[test]
fn a_problem_names_a_region_record_by_its_position() {
    let record = |vcf: &str, warning: Option<&str>| RunRecord {
        vcf: vcf.to_string(),
        manifest_index: None,
        reference_warning: warning.map(str::to_string),
    };
    let records = [
        record(FIRST, Some("odd")),
        record("21\t149\tdel1\tGC\tG", None),
        record("21\t200\tnv\tTG\t.", None),
        record("21\t600\tlast\tA\tT", None),
    ];
    let lines = vec![
        Some("21:100-100:1/G".to_string()),
        Some("21:150-150/-".to_string()),
        Some("21:200-201:-1/DEL".to_string()),
        None,
    ];
    let key = (
        "21:100".to_string(),
        "G".to_string(),
        "ENST1".to_string(),
        "Transcript".to_string(),
    );
    let names = BTreeSet::from(["21:100-100:1/G".to_string()]);
    let perl_names = TupleSets::from([(key, names)]);
    let by_name = RunNames::of(&records, lines).unwrap();
    let (skipped, problems) =
        reference_skipped_records(&records, &perl_names, &by_name, &BTreeSet::new());
    // The twin of the third region line is a non-variant VCF record, which keeps no region
    // record: a region line always carries an allele.
    assert_eq!(skipped, BTreeSet::from([1, 2, 3]));
    assert_eq!(
        problems,
        [
            "record 0 (21:100-100:1) has VEP rows and reference_warning 'odd'",
            "record 1 (21:150-150) has no VEP row and no reference_warning",
            "record 2 (21:200-201:-1) has no VEP row and no reference_warning",
            "record 3 (21 600 last A T) has no VEP row and no reference_warning",
        ]
    );
}

/// Four records of a run: the second and third share their first five columns.
fn run_of_four() -> Vec<RunRecord> {
    [
        FIRST,
        "21\t200\teq\tC\tC\t.\tPASS\t.",
        "21\t200\teq\tC\tC\t.\tPASS\t.\tGT\t0/1",
        "21\t300\tx\tC\t<DEL:*>",
    ]
    .iter()
    .map(|line| RunRecord {
        vcf: line.to_string(),
        manifest_index: None,
        reference_warning: None,
    })
    .collect()
}

/// The ordinals the JSON output of `bytes` carries, or the line of the refusal with the
/// directory written `{C}`.
fn carried(directory: &TempDir, bytes: &[u8]) -> Result<Vec<usize>, String> {
    let path = directory.path.join("json.jsonl.gz");
    std::fs::write(&path, bytes).unwrap();
    match carried_records(&run_of_four(), path.to_str().unwrap()) {
        Ok(ordinals) => Ok(ordinals.into_iter().collect()),
        Err(refusal) => Err(refusal
            .line
            .replace(directory.path.to_str().unwrap(), "{C}")),
    }
}

fn not_json(reason: &str) -> Result<Vec<usize>, String> {
    Err(format!(
        "ERROR: [build_golden_corpus] {{C}}/json.jsonl.gz is not a readable JSON output ({reason})"
    ))
}

const CARRIED: &str =
    r#"{"input": "21\t200\teq\tC\tC\t.\tPASS\t.\tGT\t0/1", "most_severe_consequence": "?"}"#;

#[test]
fn the_json_output_names_the_records_carried_without_a_row() {
    let directory = scratch("carried");
    let text = [
        CARRIED,
        "",
        r#"{"input": "21\t100\trs1\tA\tG", "most_severe_consequence": "missense_variant"}"#,
        r#"{"most_severe_consequence": "?", "input": 7}"#,
        "  ",
        "",
    ]
    .join("\n");
    assert_eq!(carried(&directory, &gzip(text.as_bytes())), Ok(vec![1, 2]));
    let last = r#"{"input": "21\t300\tx\tC\t<DEL:*>\tmore", "most_severe_consequence": "?"}"#;
    assert_eq!(carried(&directory, &gzip(last.as_bytes())), Ok(vec![3]));
    assert_eq!(carried(&directory, b""), Ok(vec![]));
    let absent = directory.path.join("nothing.gz");
    let ordinals = carried_records(&run_of_four(), absent.to_str().unwrap()).unwrap();
    assert!(ordinals.is_empty());
}

#[test]
fn a_json_output_that_cannot_be_read_is_refused_with_the_reason() {
    let directory = scratch("refused");
    let property = "Expecting property name enclosed in double quotes: line 1 column 2 (char 1)";
    assert_eq!(
        carried(&directory, &gzip(b"{not json\n")),
        not_json(property)
    );
    assert_eq!(carried(&directory, &gzip(b"{not json")), not_json(property));
    let lines = |rest: &str| gzip(format!("{CARRIED}\n{rest}").as_bytes());
    assert_eq!(
        carried(&directory, &lines("[]\n")),
        not_json("line 2 is not a JSON object")
    );
    // The reason is the one for the line with its line end, which the last line may lack.
    assert_eq!(
        carried(&directory, &lines("{\"a\":")),
        not_json("Expecting value: line 1 column 6 (char 5)")
    );
    assert_eq!(
        carried(&directory, &lines("{\"a\":\n")),
        not_json("Expecting value: line 2 column 1 (char 6)")
    );
    assert_eq!(
        carried(
            &directory,
            &gzip(format!("{CARRIED}\r\n{{\"a\": 1\r\n").as_bytes())
        ),
        not_json("Expecting ',' delimiter: line 2 column 1 (char 8)")
    );
    assert_eq!(
        carried(&directory, &gzip(b"\"abc\n")),
        not_json("Invalid control character at: line 1 column 5 (char 4)")
    );
    assert_eq!(
        carried(&directory, &gzip(b"\"abc")),
        not_json("Unterminated string starting at: line 1 column 1 (char 0)")
    );
    assert_eq!(
        carried(&directory, &gzip(b"{\"a\": 1} x\n")),
        not_json("Extra data: line 1 column 10 (char 9)")
    );
    assert_eq!(
        carried(&directory, &gzip(b"\xef\xbb\xbf{\"a\": 1}\n")),
        not_json("Unexpected UTF-8 BOM (decode using utf-8-sig): line 1 column 1 (char 0)")
    );
    assert_eq!(
        carried(&directory, b"not gzip"),
        not_json("Not a gzipped file (b'no')")
    );
    let whole = gzip(format!("{CARRIED}\n").as_bytes());
    assert_eq!(
        carried(&directory, &whole[..whole.len() - 6]),
        not_json("Compressed file ended before the end-of-stream marker was reached")
    );
}

#[test]
fn a_record_without_a_row_is_skipped_kept_or_a_problem() {
    let record = |vcf: &str, warning: Option<&str>| RunRecord {
        vcf: vcf.to_string(),
        manifest_index: None,
        reference_warning: warning.map(str::to_string),
    };
    let records = [
        record(FIRST, Some("it's \"odd\"")),
        record("21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX", Some("CPX")),
        record("21\t200\tnv\tTG\t.", None),
        record("21\t300\tc\tC\tC\t.\tPASS", None),
        record("21 400 w", None),
        record("21\t600\tlast\tA\tT\t.\t.\t.", None),
    ];
    let key = (
        "21:100".to_string(),
        "G".to_string(),
        "ENST1".to_string(),
        "Transcript".to_string(),
    );
    let names: BTreeSet<String> = ["rs1", "other"]
        .iter()
        .map(|name| name.to_string())
        .collect();
    let perl_names = TupleSets::from([(key, names)]);
    let by_name = RunNames {
        by_name: [("rs1", 0), ("sv1", 1), ("nv", 2), ("c", 3), ("w", 4)]
            .iter()
            .map(|(name, ordinal)| (name.to_string(), vec![*ordinal]))
            .collect(),
        region_lines: Vec::new(),
    };
    let warned = "record 0 (21 100 rs1 A G) has VEP rows and reference_warning 'it\\'s \"odd\"'";
    let rowless = |ordinal: usize, label: &str| {
        format!("record {ordinal} ({label}) has no VEP row and no reference_warning")
    };
    // A non-variant record and a carried one are kept; a record that stays a problem
    // stays among the skipped.
    let (skipped, problems) =
        reference_skipped_records(&records, &perl_names, &by_name, &BTreeSet::from([3]));
    assert_eq!(skipped, BTreeSet::from([1, 4, 5]));
    assert_eq!(
        problems,
        [
            warned.to_string(),
            rowless(4, "21 400 w"),
            rowless(5, "21 600 last A T")
        ]
    );
    let (skipped, problems) =
        reference_skipped_records(&records, &perl_names, &by_name, &BTreeSet::new());
    assert_eq!(skipped, BTreeSet::from([1, 3, 4, 5]));
    assert_eq!(
        problems,
        [
            warned.to_string(),
            rowless(3, "21 300 c C C"),
            rowless(4, "21 400 w"),
            rowless(5, "21 600 last A T")
        ]
    );
}
