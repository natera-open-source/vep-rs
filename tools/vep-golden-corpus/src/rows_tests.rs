// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of `rows`: the Extra column, the lines of an input VCF, the rows and the
//! records of a default-format output, the alleles of an input and the record-class pass.

use super::*;
use crate::classes::parse_family_columns;
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write as _;

const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-rows-").unwrap()
}

/// Writes `text` to `<directory>/<name>`, through gzip for a name ending in `.gz` or
/// `.bgz`, and returns the path.
fn file(directory: &TempDir, name: &str, text: &str) -> String {
    let path = directory.path.join(name);
    if name.ends_with(".gz") || name.ends_with(".bgz") {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(text.as_bytes()).unwrap();
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();
    } else {
        std::fs::write(&path, text).unwrap();
    }
    path.to_str().unwrap().to_string()
}

/// A Transcript row of fourteen columns.
fn row(name: &str, location: &str, allele: &str, feature: &str, consequence: &str) -> String {
    let extra = "IMPACT=MODIFIER;STRAND=1";
    feature_row(
        [name, location, allele],
        ["ENSG1", feature, "Transcript"],
        consequence,
        extra,
    )
}

/// A row of fourteen columns: Uploaded_variation, Location and Allele, then Gene, Feature
/// and Feature_type, the Consequence, six columns of `-` and the Extra column.
fn feature_row(record: [&str; 3], feature: [&str; 3], consequence: &str, extra: &str) -> String {
    let mut columns = record.to_vec();
    columns.extend(feature);
    columns.push(consequence);
    columns.extend(["-"; 6]);
    columns.push(extra);
    columns.join("\t") + "\n"
}

fn strings(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| item.to_string()).collect()
}

/// A chromosome, a name and the substitution ALTs under them.
type Held<'a> = (&'a str, &'a str, &'a [(i64, &'a str)]);

fn substitutions(entries: &[Held<'_>]) -> Substitutions {
    let mut held = Substitutions::default();
    for (chromosome, name, alts) in entries {
        held.extend(
            chromosome,
            name,
            alts.iter().map(|(pos, alt)| (*pos, alt.to_string())),
        );
    }
    held
}

fn rows_of(reader: Rows<'_>) -> Vec<Row> {
    reader.map(Result::unwrap).collect()
}

fn classes_of(reader: Rows<'_>) -> Vec<String> {
    rows_of(reader).iter().map(Row::klass).collect()
}

#[test]
fn the_extra_column_is_read_as_a_dict() {
    assert!(parse_extra("").is_empty() && parse_extra("-").is_empty());
    let extra = parse_extra("A=1;B;A=2;;=x;K=a=b");
    assert_eq!(
        extra.iter().collect::<Vec<_>>(),
        [("A", "2"), ("B", ""), ("", "x"), ("K", "a=b")]
    );
    assert_eq!(extra.len(), 4);
    assert_eq!(extra.get("A"), Some("2"));
    assert_eq!(extra.get("B"), Some(""));
    assert_eq!(extra.get("C"), None);
    assert_eq!(parse_extra("-;A=1").get("-"), Some(""));
    let keys = RowKeys {
        strand: Finder::new("STRAND"),
        flags: Finder::new("FLAGS"),
    };
    for raw in [
        "",
        "-",
        "STRAND=1",
        "XSTRAND=1",
        "A=STRAND=1",
        "STRAND",
        "STRAND;STRAND=2",
        "STRAND=1;STRAND",
        "A=1;STRAND=-1;FLAGS=cds_start_NF",
        "STRANDED=1",
        "STRAND=1=2",
        ";STRAND=1",
        "A=x;STRAND=",
        "FLAGS=a;FLAGS=b",
        "STRAND=1;A=STRAND",
        "STRAND=1;A=;STRAND",
        "FLAGS=cds_end_NF,cds_start_NF;STRAND=-1;IMPACT=LOW",
        "A=FLAGS;FLAGSTRAND=1",
    ] {
        let extra = parse_extra(raw);
        assert_eq!(extra_value(raw, &keys.strand), extra.get("STRAND"), "{raw}");
        assert_eq!(extra_value(raw, &keys.flags), extra.get("FLAGS"), "{raw}");
    }
}

#[test]
fn a_consequence_column_is_the_key_when_it_is_written_as_one() {
    for raw in [
        "missense_variant",
        "intron_variant,splice_region_variant",
        "",
        "a,b,c",
    ] {
        assert!(matches!(consequence_set(raw), Cow::Borrowed(_)), "{raw}");
        assert_eq!(consequence_set(raw), normalize_consequence_set(raw));
    }
    for raw in [
        "splice_region_variant,intron_variant",
        " missense_variant",
        "missense_variant\n",
        "a,,b",
        "a,a",
        ",",
        "a,",
        "b,a,B",
        "a, b",
    ] {
        assert!(matches!(consequence_set(raw), Cow::Owned(_)), "{raw}");
        assert_eq!(consequence_set(raw), normalize_consequence_set(raw));
    }
    let cases = [
        (
            "splice_region_variant,intron_variant",
            "intron_variant,splice_region_variant",
        ),
        (" missense_variant,,missense_variant\n", "missense_variant"),
        ("", ""),
        (",", ""),
        ("b,a,B", "B,a,b"),
    ];
    for (raw, expected) in cases {
        assert_eq!(normalize_consequence_set(raw), expected, "{raw:?}");
    }
}

#[test]
fn a_dataset_argument_names_an_id_a_reference_and_an_input() {
    assert_eq!(
        parse_dataset_arg("pool=ref.txt=in.vcf"),
        Ok((
            "pool".to_string(),
            "ref.txt".to_string(),
            "in.vcf".to_string()
        ))
    );
    assert_eq!(
        parse_dataset_arg("pool=a//b/=./c.vcf.gz"),
        Ok((
            "pool".to_string(),
            "a/b".to_string(),
            "c.vcf.gz".to_string()
        ))
    );
    assert_eq!(
        parse_dataset_arg("=="),
        Ok((String::new(), ".".to_string(), ".".to_string()))
    );
    for (bad, repr) in [
        ("pool=ref.txt", "'pool=ref.txt'"),
        ("a=b=c=d", "'a=b=c=d'"),
        ("", "''"),
        ("it's", "\"it's\""),
    ] {
        assert_eq!(
            parse_dataset_arg(bad),
            Err(failure(format!(
                "--dataset wants sid=reference=input, got {repr}"
            ))),
            "{bad}"
        );
    }
}

#[test]
fn vcf_lines_are_the_header_and_the_records_of_eight_fields() {
    let directory = scratch();
    let text = format!(
        "{HEADER}21\t100\trs1\tA\tG\t.\tPASS\t.\n\
         21\t200\trs2\tC\tT\t.\tPASS\n\
         \n\
         #note\n\
         21\t300\trs3\tG\tA\t.\tPASS\tAC=1\tGT:DP\t0/1:9\t1/1:8\r\n\
         21\t400\trs4\tT\tC\t.\tPASS\t.\t"
    );
    let header = strings(&[
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    ]);
    for name in ["in.vcf", "in.vcf.gz", "in.vcf.bgz"] {
        let path = file(&directory, name, &text);
        let all: Vec<(bool, Vec<String>)> = vcf_lines(&path, None).map(Result::unwrap).collect();
        assert_eq!(
            all,
            [
                (true, header.clone()),
                (
                    false,
                    strings(&["21", "100", "rs1", "A", "G", ".", "PASS", "."])
                ),
                (true, strings(&["#note"])),
                (
                    false,
                    strings(&[
                        "21", "300", "rs3", "G", "A", ".", "PASS", "AC=1", "GT:DP", "0/1:9",
                        "1/1:8"
                    ])
                ),
                (
                    false,
                    strings(&["21", "400", "rs4", "T", "C", ".", "PASS", ".", ""])
                ),
            ],
            "{name}"
        );
        let cut: Vec<(bool, Vec<String>)> = vcf_lines(&path, Some(8)).map(Result::unwrap).collect();
        assert_eq!(cut[0], (true, header.clone()));
        assert_eq!(cut[1].1.len(), 8);
        assert_eq!(cut[3].1[8], "GT:DP\t0/1:9\t1/1:8");
        assert_eq!(cut[4].1, all[4].1);
        let wider: Vec<(bool, Vec<String>)> =
            vcf_lines(&path, Some(10)).map(Result::unwrap).collect();
        assert_eq!(wider, all);
        assert_eq!(
            vcf_lines(&path, Some(9))
                .map(Result::unwrap)
                .nth(3)
                .unwrap()
                .1[9],
            "0/1:9\t1/1:8"
        );
        // Split into fewer than eight fields, no line is a record.
        let headers: Vec<(bool, Vec<String>)> =
            vcf_lines(&path, Some(6)).map(Result::unwrap).collect();
        assert_eq!(
            headers,
            [(true, header.clone()), (true, strings(&["#note"]))]
        );
        let mut reader = vcf_lines(&path, Some(8));
        let mut borrowed = Vec::new();
        while let Some(line) = reader.next_columns() {
            let (is_header, columns) = line.unwrap();
            borrowed.push((is_header, strings(&columns)));
        }
        assert_eq!(borrowed, cut);
        assert!(reader.next_columns().is_none() && reader.next().is_none());
    }
}

#[test]
fn an_unreadable_vcf_is_refused_with_its_path() {
    let directory = scratch();
    let absent = directory.path.join("absent.vcf");
    let absent = absent.to_str().unwrap();
    let refusal = failure(format!(
        "{absent} is not a readable text file ([Errno 2] No such file or directory: '{absent}')"
    ));
    let mut reader = vcf_lines(absent, None);
    assert_eq!(reader.next(), Some(Err(refusal.clone())));
    assert_eq!(reader.next(), None);
    assert_eq!(input_alleles(absent), Err(refusal.clone()));
    assert_eq!(substitution_alts(absent), Err(refusal.clone()));
    let mut rows = iter_rows::<&str>(absent, "pool", &[], None);
    assert_eq!(rows.next(), Some(Err(refusal.clone())));
    assert_eq!(rows.next(), None);
    let mut records = iter_records(absent, "pool", None);
    assert_eq!(records.next(), Some(Err(refusal)));
    assert_eq!(records.next(), None);

    let plain = file(
        &directory,
        "plain.vcf",
        &format!("{HEADER}21\t1\t.\tA\tG\t.\t.\t.\n"),
    );
    let named_gzip = directory.path.join("plain.vcf.gz");
    std::fs::copy(&plain, &named_gzip).unwrap();
    let named_gzip = named_gzip.to_str().unwrap();
    assert_eq!(
        vcf_lines(named_gzip, Some(8)).collect::<Vec<_>>(),
        [Err(failure(format!(
            "{named_gzip} is not a readable text file (Not a gzipped file (b'##'))"
        )))]
    );

    let records: String = (0..20_000)
        .map(|number| format!("21\t{number}\trs{number}\tA\tG\t.\tPASS\t.\n"))
        .collect();
    let whole = file(&directory, "whole.vcf.gz", &format!("{HEADER}{records}"));
    let bytes = std::fs::read(&whole).unwrap();
    let cut = directory.path.join("cut.vcf.gz");
    std::fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
    let cut = cut.to_str().unwrap();
    let lines: Vec<_> = vcf_lines(cut, Some(8)).collect();
    let (last, before) = lines.split_last().unwrap();
    assert_eq!(
        *last,
        Err(failure(format!(
            "{cut} is not a readable text file (Compressed file ended before the \
             end-of-stream marker was reached)"
        )))
    );
    assert!(before.len() > 100 && before.iter().all(Result::is_ok));
    assert!(substitution_alts(cut).is_err() && input_alleles(cut).is_err());
}

#[test]
fn kept_samples_are_the_first_names_after_format() {
    let header = strings(&[
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2", "S3",
    ]);
    let kept = |keep: usize, samples: &[&str]| {
        kept_samples(&header, keep, "in.vcf", &strings(samples), EARLIER_INPUT)
    };
    assert_eq!(kept(2, &[]), Ok(strings(&["S1", "S2"])));
    assert_eq!(kept(2, &["S1", "S2"]), Ok(strings(&["S1", "S2"])));
    assert_eq!(kept(0, &[]), Ok(Vec::new()));
    assert_eq!(kept(3, &[]), Ok(strings(&["S1", "S2", "S3"])));
    assert_eq!(
        kept(4, &[]),
        Err(failure("in.vcf has 3 sample columns, 4 are kept"))
    );
    assert_eq!(
        kept(2, &["S1", "X2"]),
        Err(failure(
            "in.vcf keeps samples S1,S2 where an earlier input keeps S1,X2: every input must \
             keep the same samples"
        ))
    );
    assert_eq!(
        kept_samples(
            &header,
            1,
            "in.vcf",
            &strings(&["S1", "S2"]),
            "the manifest"
        ),
        Err(failure(
            "in.vcf keeps samples S1 where the manifest keeps S1,S2: every input must keep the \
             same samples"
        ))
    );
    for short in [
        &header[..8],
        &[
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "S1",
        ]
        .map(str::to_string)[..],
    ] {
        assert_eq!(
            kept_samples(short, 0, "in.vcf", &[], EARLIER_INPUT),
            Err(failure("in.vcf has no FORMAT column, so no sample columns"))
        );
    }
}

#[test]
fn rows_are_read_with_the_columns_a_line_has() {
    let directory = scratch();
    let text = "## ENSEMBL VARIANT EFFECT PREDICTOR\n\
                #Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n\
                rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tsplice_region_variant,intron_variant\t\
                12\t?-7\t3\tA/G\tgCa/gGa\trs1\tIMPACT=LOW;STRAND=-1;FLAGS=cds_start_NF\textra\tmore\n\
                rs1\t21:100\tG\tENSG1\tENST2\tTranscript\tmissense_variant\t5-6\n\
                \n\
                rs2\t21:200-300\tdeletion\t-\t-\t-\tintergenic_variant\n\
                rs3\t21:400\tG\tENSG1\tENST1\tTranscript\n\
                rs4\tX:5\tA\t-\tENSR1\tRegulatoryFeature\t regulatory_region_variant \t-\t-\t-\t-\t-\t-\tSTRAND=1;STRAND=;IMPACT=MODIFIER";
    for name in ["ref.txt", "ref.txt.gz"] {
        let path = file(&directory, name, text);
        let rows = rows_of(iter_rows::<&str>(&path, "pool", &[], None));
        let first = Row {
            dataset: "pool".to_string(),
            uploaded_variation: "rs1".to_string(),
            location: "21:100".to_string(),
            allele: "G".to_string(),
            feature: "ENST1".to_string(),
            feature_type: "Transcript".to_string(),
            consequence_set: "intron_variant,splice_region_variant".to_string(),
            strand: "-1".to_string(),
            flags: "cds_start_NF".to_string(),
            positions: ["12", "?-7", "3"].map(str::to_string),
            family: Vec::new(),
            substitution: false,
        };
        assert_eq!(rows.len(), 4, "{name}");
        assert_eq!(rows[0], first);
        assert_eq!(rows[0].record_key(), ("rs1", "21:100"));
        assert_eq!((rows[0].preference(), rows[0].span()), (9, 0));
        // A row of eight columns: the others read as `-`.
        assert_eq!(rows[1].positions, ["5-6", "-", "-"].map(str::to_string));
        assert_eq!((rows[1].strand.as_str(), rows[1].flags.as_str()), ("", ""));
        assert_eq!(
            (rows[1].preference(), rows[1].klass()),
            (2, "snv".to_string())
        );
        // A row of seven columns.
        assert_eq!(rows[2].consequence_set, "intergenic_variant");
        assert_eq!(rows[2].positions, ["-", "-", "-"].map(str::to_string));
        assert_eq!(
            (rows[2].span(), rows[2].klass(), rows[2].preference()),
            (100, "symbolic_deletion".to_string(), 0)
        );
        // The last line needs no line end; a field given twice holds its last value.
        assert_eq!(rows[3].uploaded_variation, "rs4");
        assert_eq!(rows[3].consequence_set, "regulatory_region_variant");
        assert_eq!(rows[3].strand, "");
        assert_eq!(rows[3].family_stratum(), "snv");

        let specs = parse_family_columns(Some("IMPACT=,FLAGS,STRAND=")).unwrap();
        let mut reader = iter_rows(&path, "pool", &specs, None);
        let mut families = Vec::new();
        let mut strata = Vec::new();
        while let Some(view) = reader.next_view() {
            let view = view.unwrap();
            assert_eq!(view.to_row().family_stratum(), view.family_stratum());
            assert_eq!(view.to_row().preference(), view.preference());
            assert_eq!(view.to_row().span(), view.span());
            assert_eq!(view.to_row().klass(), view.klass());
            assert_eq!(view.to_row().record_key(), view.record_key());
            families.push(view.family.clone());
            strata.push(view.family_stratum());
        }
        assert_eq!(
            families,
            [
                ["IMPACT=LOW", "FLAGS+", "STRAND=-1"],
                ["IMPACT=-", "FLAGS-", "STRAND=-"],
                ["IMPACT=-", "FLAGS-", "STRAND=-"],
                ["IMPACT=MODIFIER", "FLAGS-", "STRAND=-"],
            ]
        );
        assert_eq!(strata[0], "IMPACT=LOW;FLAGS+;STRAND=-1;snv");
        assert_eq!(strata[2], "IMPACT=-;FLAGS-;STRAND=-;symbolic_deletion");
        assert!(reader.next_view().is_none() && reader.next().is_none());
    }
}

#[test]
fn a_row_measures_its_span_and_its_preference() {
    let row = |location: &str, flags: &str, positions: [&str; 3]| Row {
        dataset: String::new(),
        uploaded_variation: String::new(),
        location: location.to_string(),
        allele: "G".to_string(),
        feature: String::new(),
        feature_type: String::new(),
        consequence_set: String::new(),
        strand: String::new(),
        flags: flags.to_string(),
        positions: positions.map(str::to_string),
        family: strings(&["A+", "B=x"]),
        substitution: false,
    };
    let span = |location: &str| row(location, "", ["-"; 3]).span();
    assert_eq!(span("21:100"), 0);
    assert_eq!(span("21:100-105"), 5);
    assert_eq!(span("21:105-100"), 5);
    assert_eq!(span("100-105"), 5);
    assert_eq!(span("HLA:1:100-105"), 0);
    assert_eq!(span("21:a-b"), 0);
    assert_eq!(span("21:-5"), 0);
    assert_eq!(span("21:1-2-3"), 0);
    assert_eq!(span("21: 1_0 - +25 "), 15);
    assert_eq!(span("21:-9223372036854775808-9223372036854775807"), 0);
    let preference = |flags: &str, positions: [&str; 3]| row("21:1", flags, positions).preference();
    assert_eq!(preference("", ["-", "-", "-"]), 0);
    assert_eq!(preference("x", ["-", "-", "-"]), 4);
    assert_eq!(preference("", ["1", "?", "-"]), 3);
    assert_eq!(preference("", ["1-2", "5", "-"]), 2);
    assert_eq!(preference("x", ["?-2", "5", "-"]), 9);
    assert_eq!(row("21:1", "", ["-"; 3]).family_stratum(), "A+;B=x;snv");
}

#[test]
fn a_two_base_substitution_is_classed_from_its_input_record() {
    let directory = scratch();
    let reference = file(
        &directory,
        "ref.txt",
        &[
            row(
                "mnv1",
                "21:100-101",
                "TG",
                "ENST1",
                "coding_sequence_variant",
            ),
            row(
                "ins1",
                "21:200-201",
                "TG",
                "ENST1",
                "coding_sequence_variant",
            ),
            row(
                "21_300_AC/TG",
                "21:300-301",
                "TG",
                "ENST1",
                "coding_sequence_variant",
            ),
        ]
        .concat(),
    );
    let input = file(
        &directory,
        "in.vcf",
        &format!(
            "{HEADER}21\t100\tmnv1\tAC\tTG\t.\tPASS\t.\n\
             21\t200\tins1\tA\tATG\t.\tPASS\t.\n\
             21\t300\t.\tAC\tTG\t.\tPASS\t.\n"
        ),
    );
    let alts = substitution_alts(&input).unwrap();
    assert_eq!(
        alts,
        substitutions(&[
            ("21", "mnv1", &[(100, "TG")]),
            ("21", "21_300_AC/TG", &[(300, "TG")])
        ])
    );
    assert_eq!(alts.len(), 2);
    assert_eq!(alts.get("21", "mnv1"), Some(&[(100, "TG".to_string())][..]));
    assert_eq!(alts.get("21", "ins1"), None);
    assert_eq!(alts.get("22", "mnv1"), None);
    assert_eq!(
        classes_of(iter_rows::<&str>(&reference, "pool", &[], Some(&alts))),
        ["mnv_or_complex", "insertion", "mnv_or_complex"]
    );
    assert_eq!(
        classes_of(iter_rows::<&str>(&reference, "pool", &[], None)),
        ["insertion"; 3]
    );
    let none = Substitutions::default();
    assert_eq!(
        classes_of(iter_rows::<&str>(&reference, "pool", &[], Some(&none))),
        ["insertion"; 3]
    );
    let held = [(100, "TG".to_string())];
    assert!(holds_allele("100-101", "tg", &held));
    assert!(!holds_allele("100-101", "CA", &held));
    // An allele of another length, and an allele shorter than its span.
    assert!(!holds_allele("100-101", "TGA", &held));
    assert!(!holds_allele("101-102", "G", &held));
    assert!(!holds_allele("x-y", "TG", &held));
    assert!(holds_allele("101", "G", &held));
    assert!(!holds_allele("99-100", "AT", &held));
    assert!(!holds_allele("100-101", "TG", &[]));
    assert!(holds_allele(" 100 - 101 ", "TG", &held));

    // Bases the alleles share are dropped: the anchor base of a multi-allelic record with
    // an indel, either end under --minimal. A record that only shares the row's name, here
    // an insertion, is not read.
    let shared_reference = file(
        &directory,
        "shared.txt",
        &[
            ("21:100-101", "TG"),
            ("21:200-201", "TG"),
            ("21:301-302", "CG"),
            ("21:401-402", "GG"),
        ]
        .map(|(location, allele)| row("rs1", location, allele, "ENST1", "coding_sequence_variant"))
        .concat(),
    );
    let shared_input = file(
        &directory,
        "shared.vcf",
        &format!(
            "{HEADER}21\t100\trs1\tAC\tTG\t.\tPASS\t.\n\
             21\t200\trs1\tA\tATG\t.\tPASS\t.\n\
             21\t300\trs1\tATG\tAG,ACG,<DEL>,*\t.\tPASS\t.\n\
             21\t400\trs1\tAAAG\tAGGG\t.\tPASS\t.\n"
        ),
    );
    let shared = substitution_alts(&shared_input).unwrap();
    assert_eq!(
        shared,
        substitutions(&[("21", "rs1", &[(100, "TG"), (300, "ACG"), (400, "AGGG")])])
    );
    assert_eq!(
        classes_of(iter_rows::<&str>(
            &shared_reference,
            "pool",
            &[],
            Some(&shared)
        )),
        [
            "mnv_or_complex",
            "insertion",
            "mnv_or_complex",
            "mnv_or_complex"
        ]
    );
    let alleles = input_alleles(&shared_input).unwrap();
    let pair = |ref_allele: &str, alt: &str| (ref_allele.to_string(), alt.to_string());
    assert_eq!(
        alleles.get("21", "rs1").unwrap(),
        [
            pair("AC", "TG"),
            pair("A", "ATG"),
            pair("ATG", "AG"),
            pair("ATG", "ACG"),
            pair("ATG", "<DEL>"),
            pair("ATG", "*"),
            pair("AAAG", "AGGG"),
        ]
    );
    assert_eq!(alleles.len(), 1);
    let alleles = input_alleles(&input).unwrap();
    assert_eq!(alleles.len(), 3);
    assert_eq!(
        alleles.get("21", "21_300_AC/TG").unwrap(),
        [pair("AC", "TG")]
    );
}

#[test]
fn an_input_record_without_a_position_is_refused_where_its_names_are_read() {
    let directory = scratch();
    let input = file(
        &directory,
        "in.vcf",
        &format!(
            "{HEADER}21\t100\trs1\tA\tG\t.\tPASS\t.\n21\tx\trs2\tAC\tTG\t.\tPASS\t.\n\
             21\ty\trs3\tA\tG\t.\tPASS\t.\n"
        ),
    );
    let refusal = |pos: &str, id: &str| {
        Failure::own(format!(
            "the record 21 {pos} {id} has a POS that is no integer of 64 bits"
        ))
    };
    assert_eq!(input_alleles(&input), Err(refusal("x", "rs2")));
    assert_eq!(substitution_alts(&input), Err(refusal("x", "rs2")));
    // A record without a substitution is not read for its position.
    let input = file(
        &directory,
        "single.vcf",
        &format!("{HEADER}21\ty\trs3\tA\tG\t.\tPASS\t.\n21\tz\trs4\tAC\tT\t.\tPASS\t.\n"),
    );
    assert_eq!(substitution_alts(&input), Ok(Substitutions::default()));
    assert_eq!(input_alleles(&input), Err(refusal("y", "rs3")));
}

#[test]
fn records_are_the_consecutive_rows_of_one_name_and_location() {
    let directory = scratch();
    let reference = file(
        &directory,
        "ref.txt",
        &[
            row("rs1", "21:100", "G", "ENST1", "intron_variant"),
            row("rs1", "21:100", "G", "ENST2", "missense_variant"),
            row("rs1", "21:200", "G", "ENST1", "intron_variant"),
            row("rs2", "21:200", "G", "ENST1", "intron_variant"),
            row("rs1", "21:100", "T", "ENST1", "intron_variant"),
        ]
        .concat(),
    );
    let groups: Vec<Vec<(String, String)>> = iter_records(&reference, "pool", None)
        .map(|rows| {
            rows.unwrap()
                .iter()
                .map(|row| (row.uploaded_variation.clone(), row.feature.clone()))
                .collect()
        })
        .collect();
    let named = |name: &str, feature: &str| (name.to_string(), feature.to_string());
    assert_eq!(
        groups,
        [
            vec![named("rs1", "ENST1"), named("rs1", "ENST2")],
            vec![named("rs1", "ENST1")],
            vec![named("rs2", "ENST1")],
            vec![named("rs1", "ENST1")],
        ]
    );
    let empty = file(&directory, "empty.txt", "#Uploaded_variation\tLocation\n");
    assert_eq!(iter_records(&empty, "pool", None).count(), 0);
}

/// The step of a sample that never fills: the item goes behind the ones the bucket holds.
fn append(bucket: &mut Vec<RecordRows>, item: RecordRows, offered: usize, capacity: usize) {
    assert!(offered == bucket.len() + 1 && offered <= capacity);
    bucket.push(item);
}

#[test]
fn the_record_class_pass_classes_each_record_by_its_regulatory_rows() {
    let directory = scratch();
    let (rf, tf) = ("regulatory_region_variant", "TF_binding_site_variant");
    let regulatory = |name: &str, location: &str, allele: &str, feature: &str| {
        feature_row(
            [name, location, allele],
            ["-", feature, "RegulatoryFeature"],
            rf,
            "IMPACT=MODIFIER",
        )
    };
    let motif = |name: &str, location: &str, allele: &str, feature: &str, strand: &str| {
        let extra = format!("IMPACT=MODIFIER;STRAND={strand};MOTIF_NAME=ENSPFM1");
        feature_row(
            [name, location, allele],
            ["-", feature, "MotifFeature"],
            tf,
            &extra,
        )
    };
    let intergenic = |name: &str, location: &str, allele: &str| {
        feature_row(
            [name, location, allele],
            ["-", "-", "-"],
            "intergenic_variant",
            "IMPACT=MODIFIER",
        )
    };
    let reference = file(
        &directory,
        "ref.txt",
        &[
            row("rs_rf", "21:1000", "G", "ENST1", "intron_variant"),
            regulatory("rs_rf", "21:1000", "G", "ENSR1"),
            row("rs_rf2", "21:1100", "G", "ENST1", "intron_variant"),
            regulatory("rs_rf2", "21:1100", "G", "ENSR1"),
            row("rs_mf", "21:2000", "G", "ENST1", "intron_variant"),
            motif("rs_mf", "21:2000", "G", "ENSM1", "1"),
            row("rs_both", "21:3000", "G", "ENST1", "intron_variant"),
            regulatory("rs_both", "21:3000", "G", "ENSR2"),
            motif("rs_both", "21:3000", "G", "ENSM2", "-1"),
            regulatory("rs_regonly", "21:4000", "G", "ENSR3"),
            intergenic("rs_regonly", "21:4000", "G"),
            row(
                "21_5000_AT/GC",
                "21:5000-5001",
                "GC",
                "ENST1",
                "intron_variant",
            ),
            motif("21_5000_AT/GC", "21:5000-5001", "GC", "ENSM3", "1"),
            row("21_6000_AT/A", "21:6001", "-", "ENST1", "intron_variant"),
            regulatory("21_6000_AT/A", "21:6001", "-", "ENSR4"),
            motif("21_6000_AT/A", "21:6001", "-", "ENSM4", "1"),
            row("rs_plain", "21:7000", "G", "ENST1", "intron_variant"),
            motif("21_9000_A/G/T", "21:9000", "G", "ENSM5", "1"),
            intergenic("21_9000_A/G/T", "21:9000", "G"),
            motif("21_9000_A/G/T", "21:9000", "T", "ENSM5", "1"),
            intergenic("21_9000_A/G/T", "21:9000", "T"),
        ]
        .concat(),
    );
    let input = file(
        &directory,
        "in.vcf",
        &format!(
            "{HEADER}21\t1000\trs_rf\tA\tG\t.\tPASS\t.\n\
             21\t1100\trs_rf2\tA\tG\t.\tPASS\t.\n\
             21\t2000\trs_mf\tA\tG\t.\tPASS\t.\n\
             21\t3000\trs_both\tA\tG\t.\tPASS\t.\n\
             21\t4000\trs_regonly\tA\tG\t.\tPASS\t.\n\
             21\t5000\t.\tAT\tGC\t.\tPASS\t.\n\
             21\t6000\t.\tAT\tA\t.\tPASS\t.\n\
             21\t7000\trs_plain\tA\tG\t.\tPASS\t.\n\
             21\t9000\t.\tA\tG,T\t.\tPASS\t.\n"
        ),
    );
    let datasets = [("reg".to_string(), reference, input)];
    let pass = record_class_pass(&datasets, "regulatory_class", &mut append).unwrap();
    let classes = |name: &str, location: &str| {
        let key = ("reg".to_string(), name.to_string(), location.to_string());
        pass.classes_of[&key].iter().copied().collect::<Vec<_>>()
    };
    assert_eq!(classes("rs_rf", "21:1000"), ["rf_only"]);
    assert_eq!(classes("rs_rf2", "21:1100"), ["rf_only"]);
    assert_eq!(classes("rs_mf", "21:2000"), ["mf_only", "motif_snv"]);
    assert_eq!(
        classes("rs_both", "21:3000"),
        ["motif_reverse", "motif_snv", "rf_and_mf"]
    );
    assert_eq!(
        classes("rs_regonly", "21:4000"),
        ["regulatory_only", "rf_only"]
    );
    assert_eq!(
        classes("21_5000_AT/GC", "21:5000-5001"),
        ["mf_only", "motif_mnv"]
    );
    assert_eq!(
        classes("21_6000_AT/A", "21:6001"),
        ["motif_indel", "rf_and_mf"]
    );
    assert_eq!(
        classes("21_9000_A/G/T", "21:9000"),
        ["mf_only", "motif_snv", "regulatory_only"]
    );
    assert_eq!(pass.classes_of.len(), 8);
    assert_eq!(
        pass.seen
            .iter()
            .map(|(class, count)| (*class, *count))
            .collect::<Vec<_>>(),
        [
            ("mf_only", 3),
            ("motif_indel", 1),
            ("motif_mnv", 1),
            ("motif_reverse", 1),
            ("motif_snv", 3),
            ("regulatory_only", 2),
            ("rf_and_mf", 2),
            ("rf_only", 3),
        ]
    );
    let sampled: Vec<Vec<&str>> = pass.samples["rf_only"]
        .iter()
        .map(|rows| rows.iter().map(|row| row.feature.as_str()).collect())
        .collect();
    assert_eq!(
        sampled,
        [
            vec!["ENST1", "ENSR1"],
            vec!["ENST1", "ENSR1"],
            vec!["ENSR3", "-"]
        ]
    );
    assert_eq!(pass.samples["mf_only"][2].len(), 4);
    assert_eq!(pass.samples.len(), 8);

    assert_eq!(
        record_class_pass(&datasets, "splice_class", &mut append),
        Err(Failure::own("'splice_class' is no record-class key"))
    );
    let absent = directory.path.join("absent.vcf");
    let absent = absent.to_str().unwrap();
    let no_input = [(
        "reg".to_string(),
        "absent.txt".to_string(),
        absent.to_string(),
    )];
    assert_eq!(
        record_class_pass(&no_input, "regulatory_class", &mut append),
        Err(failure(format!(
            "{absent} is not a readable text file ([Errno 2] No such file or directory: '{absent}')"
        )))
    );
}

#[test]
fn odd_inputs_are_read_as_their_lines_say() {
    let directory = scratch();
    // A span that runs backwards by one holds the empty allele of a record at or before it.
    let held = [(100, "TG".to_string())];
    assert!(holds_allele("101-100", "", &held));
    assert!(!holds_allele("99-98", "", &held));
    assert!(!holds_allele("", "", &held));

    // After a byte order mark the header line of a reference is a row.
    let marked = file(
        &directory,
        "bom.txt",
        "\u{feff}#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n\
         rs1\t21:1\tG\t-\t-\t-\tintergenic_variant\n",
    );
    let names: Vec<String> = rows_of(iter_rows::<&str>(&marked, "pool", &[], None))
        .into_iter()
        .map(|row| row.uploaded_variation)
        .collect();
    assert_eq!(names, ["\u{feff}#Uploaded_variation", "rs1"]);

    // After a byte order mark a `##` line is no header and no record; a substitution is any
    // ALT of ASCII letters as long as its REF.
    let input = file(
        &directory,
        "bom.vcf",
        &format!("\u{feff}{HEADER}21\t5\trs1\tNN\tac,A1\t.\t.\t.\n21\t7\trs1\tAC\tNN\t.\t.\t.\n"),
    );
    assert_eq!(vcf_lines(&input, None).count(), 3);
    assert_eq!(
        substitution_alts(&input).unwrap(),
        substitutions(&[("21", "rs1", &[(5, "ac"), (7, "NN")])])
    );

    // After a byte order mark the `#CHROM` line is a record, whose POS is no number.
    let input = file(
        &directory,
        "bom2.vcf",
        "\u{feff}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t5\trs1\tA\tG\t.\t.\t.\n",
    );
    let lines: Vec<(bool, Vec<String>)> = vcf_lines(&input, None).map(Result::unwrap).collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(
        (lines[0].0, lines[0].1[0].as_str()),
        (false, "\u{feff}#CHROM")
    );
    let refusal =
        Failure::own("the record \u{feff}#CHROM POS ID has a POS that is no integer of 64 bits");
    assert_eq!(input_alleles(&input), Err(refusal.clone()));
    assert_eq!(substitution_alts(&input), Err(refusal));

    // A carriage return inside a row ends it there; what follows is a line of its own.
    let cut = file(
        &directory,
        "cr.txt",
        "rs1\t21:1\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MOD\rIFIER;STRAND=1\n",
    );
    let rows = rows_of(iter_rows::<&str>(&cut, "pool", &[], None));
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].strand, "");

    // A record whose rows come in two places is counted in both and classed by the later.
    let regulatory = |name: &str, location: &str, feature: &str, kind: &str, extra: &str| {
        feature_row(
            [name, location, "G"],
            ["-", feature, kind],
            "regulatory_region_variant",
            extra,
        )
    };
    let reference = file(
        &directory,
        "twice.txt",
        &[
            regulatory(
                "rs1",
                "21:5",
                "ENSR1",
                "RegulatoryFeature",
                "IMPACT=MODIFIER",
            ),
            regulatory(
                "rs2",
                "21:9",
                "ENSR1",
                "RegulatoryFeature",
                "IMPACT=MODIFIER",
            ),
            regulatory("rs1", "21:5", "ENSM1", "MotifFeature", "STRAND=-1"),
        ]
        .concat(),
    );
    let input = file(
        &directory,
        "twice.vcf",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t5\trs1\tA\tG\t.\t.\t.\n\
         21\t9\trs2\tA\tG\t.\t.\t.\n",
    );
    let datasets = [("pool".to_string(), reference, input)];
    let pass = record_class_pass(&datasets, "regulatory_class", &mut append).unwrap();
    let key = ("pool".to_string(), "rs1".to_string(), "21:5".to_string());
    assert_eq!(
        pass.classes_of[&key].iter().copied().collect::<Vec<_>>(),
        ["mf_only", "motif_reverse", "motif_snv", "regulatory_only"]
    );
    assert_eq!(pass.classes_of.len(), 2);
    assert_eq!(
        pass.seen
            .iter()
            .map(|(class, count)| (*class, *count))
            .collect::<Vec<_>>(),
        [
            ("mf_only", 1),
            ("motif_reverse", 1),
            ("motif_snv", 1),
            ("regulatory_only", 3),
            ("rf_only", 2),
        ]
    );
    let features: Vec<&str> = pass.samples["regulatory_only"]
        .iter()
        .map(|rows| rows[0].feature.as_str())
        .collect();
    assert_eq!(features, ["ENSR1", "ENSR1", "ENSM1"]);
}
