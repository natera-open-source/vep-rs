// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a SpliceAI score record of several entries, on planted rows and made-up
//! score files read through their tabix index.

mod compare_common;

use compare_common::at;
use compare_common::count_map;
use compare_common::counts;
use compare_common::default_row;
use compare_common::default_text;
use compare_common::documented_with;
use compare_common::scratch;
use compare_common::verdict;
use flate2::write::DeflateEncoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use indexmap::IndexMap;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Result;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::spliceai::spliceai_entries_at;
use vep_formats_compare::rules::spliceai::spliceai_multi_entry_cell_ground;
use vep_formats_compare::rules::spliceai::spliceai_refusal_lines;
use vep_formats_compare::rules::spliceai::SpliceAiTable;
use vep_formats_compare::rules::spliceai::SPLICEAI_GROUNDS;
use vep_formats_compare::rules::spliceai::SPLICEAI_MEMBERS;
use vep_formats_compare::rules::Side;

const SPLICEAI_MULTI: &str = "spliceai_multi_entry_cell";
/// The rows of the deletion 21 5000 GCA G stand here with the allele written `-`.
const SPLICEAI_SITE: &str = "21:5001-5002";
/// Three transcripts and their genes; no score entry names the third's.
const SPLICEAI_ROW_GENES: [(&str, &str); 3] =
    [("ENST1", "GENEA"), ("ENST2", "GENEB"), ("ENST3", "GENEX")];
const SCORE_HEADER: &str = "##fileformat=VCFv4.2\n##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAIv1.3 synthetic fixture, fabricated values. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

/// Made-up scores of three symbols and of any other.
fn scores(symbol: &str) -> &'static str {
    match symbol {
        "GENEA" => "0.15|0.08|0.06|0.17|-6|-41|-17|-9",
        "GENEB" => "0.02|0.12|0.10|0.05|-25|-23|37|17",
        _ => "0.31|0.00|0.44|0.09|3|-2|11|-30",
    }
}

/// One BGZF block of `data`: a gzip member whose `BC` extra field holds its size.
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut deflater = DeflateEncoder::new(Vec::new(), Compression::new(6));
    deflater.write_all(data).unwrap();
    let body = deflater.finish().unwrap();
    let mut crc = flate2::Crc::new();
    crc.update(data);
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00".to_vec();
    block.extend_from_slice(&u16::try_from(body.len() + 25).unwrap().to_le_bytes());
    block.extend_from_slice(&body);
    block.extend_from_slice(&crc.sum().to_le_bytes());
    block.extend_from_slice(&u32::try_from(data.len()).unwrap().to_le_bytes());
    block
}

/// A made-up table of `lines` under `header`, each line a BGZF block of its own, with a
/// tabix index on the second column that files every row of a chromosome in its first bin.
fn write_table(directory: &Path, lines: &[String], header: &str, name: &str) -> PathBuf {
    let path = directory.join(name);
    let mut blocks = vec![bgzf_block(format!("{header}\n").as_bytes())];
    let mut spans: IndexMap<String, (u64, u64)> = IndexMap::new();
    let mut at = blocks[0].len() as u64;
    for line in lines {
        let block = bgzf_block(format!("{line}\n").as_bytes());
        let end = at + block.len() as u64;
        let chromosome = line.split('\t').next().unwrap().to_string();
        spans.entry(chromosome).or_insert((at, at)).1 = end;
        at = end;
        blocks.push(block);
    }
    blocks.push(bgzf_block(b""));
    std::fs::write(&path, blocks.concat()).unwrap();
    let names: Vec<u8> = spans
        .keys()
        .flat_map(|chromosome| chromosome.bytes().chain(std::iter::once(0)))
        .collect();
    let mut index = b"TBI\x01".to_vec();
    let header_numbers = [
        spans.len() as i32,
        0,
        1,
        2,
        2,
        i32::from(b'#'),
        0,
        names.len() as i32,
    ];
    for number in header_numbers {
        index.extend_from_slice(&number.to_le_bytes());
    }
    index.extend_from_slice(&names);
    for (first, end) in spans.values() {
        index.extend_from_slice(&1i32.to_le_bytes());
        index.extend_from_slice(&0u32.to_le_bytes());
        index.extend_from_slice(&1i32.to_le_bytes());
        index.extend_from_slice(&(first << 16).to_le_bytes());
        index.extend_from_slice(&(end << 16).to_le_bytes());
        index.extend_from_slice(&0i32.to_le_bytes());
    }
    let mut zipped = GzEncoder::new(Vec::new(), Compression::default());
    zipped.write_all(&index).unwrap();
    std::fs::write(tabix_index_of(&path), zipped.finish().unwrap()).unwrap();
    path
}

fn tabix_index_of(path: &Path) -> PathBuf {
    PathBuf::from(format!("{}.tbi", path.display()))
}

fn opened(paths: &[PathBuf], cache: Option<&Path>) -> Result<SpliceAiTable> {
    let paths: Vec<&Path> = paths.iter().map(PathBuf::as_path).collect();
    SpliceAiTable::new(&paths, cache)
}

/// A made-up candidate cache that names the gene of each transcript of
/// `SPLICEAI_ROW_GENES`.
fn spliceai_cache(directory: &Path) -> PathBuf {
    let shards = directory.join("transcripts").join("21");
    std::fs::create_dir_all(&shards).unwrap();
    let transcripts: Vec<String> = SPLICEAI_ROW_GENES
        .iter()
        .map(|(feature, gene)| {
            format!("{{\"stable_id\": \"{feature}\", \"gene_symbol\": \"{gene}\"}}")
        })
        .collect();
    std::fs::write(
        shards.join("1-1000000.json"),
        format!("[{}]", transcripts.join(", ")),
    )
    .unwrap();
    directory.to_path_buf()
}

/// How the entries of the symbols at one allele are filed in records.
#[derive(Clone, Copy)]
enum Join {
    /// Every entry in one record.
    All,
    /// That many entries a record.
    Each(usize),
    /// A record an entry.
    Nothing,
}

/// The paths of a made-up indel score file: the records of the symbols' entries at
/// 21 5000 GCA G, joined as `join` says. `info` replaces the INFO column of an unjoined
/// first record.
fn indel_file(directory: &Path, symbols: &[&str], join: Join, info: Option<&str>) -> Vec<PathBuf> {
    let entries: Vec<String> = symbols
        .iter()
        .map(|symbol| format!("G|{symbol}|{}", scores(symbol)))
        .collect();
    let size = match join {
        Join::All => entries.len(),
        Join::Each(size) => size,
        Join::Nothing => 1,
    };
    let mut lines: Vec<String> = entries
        .chunks(size)
        .map(|joined| format!("21\t5000\t.\tGCA\tG\t.\t.\tSpliceAI={}", joined.join(",")))
        .collect();
    if let (Join::Nothing, Some(info)) = (join, info) {
        lines[0] = format!("21\t5000\t.\tGCA\tG\t.\t.\t{info}");
    }
    std::fs::create_dir_all(directory).unwrap();
    vec![write_table(
        directory,
        &lines,
        SCORE_HEADER,
        "indel.synthetic.vcf.gz",
    )]
}

/// A comparison's side tables over a made-up indel score file (`indel_file`), GENEA's and
/// GENEB's entries in one record unless a test asks otherwise.
fn spliceai_table(directory: &Path) -> Side {
    spliceai_table_of(directory, &["GENEA", "GENEB"], Join::All, None, None)
}

fn spliceai_table_of(
    directory: &Path,
    symbols: &[&str],
    join: Join,
    info: Option<&str>,
    cache: Option<&Path>,
) -> Side {
    side_of(
        opened(&indel_file(directory, symbols, join, info), cache).expect("the score file opens"),
    )
}

fn side_of(table: SpliceAiTable) -> Side {
    Side {
        spliceai_table: Some(table),
        ..Side::default()
    }
}

fn table_of(side: &Side) -> &SpliceAiTable {
    side.spliceai_table
        .as_ref()
        .expect("the side holds the score files")
}

/// What the rows of `spliceai_rows` state beside their values.
#[derive(Clone, Copy)]
struct Rows<'a> {
    symbols: bool,
    name: &'a str,
    allele: &'a str,
    site: &'a str,
}

const DELETION: Rows = Rows {
    symbols: true,
    name: "del_two_genes",
    allele: "-",
    site: SPLICEAI_SITE,
};
const UNNAMED: Rows = Rows {
    symbols: false,
    ..DELETION
};
/// The rows of the substitution 21 7000 A C stand at 21:7000 with the allele C.
const BASE: Rows = Rows {
    symbols: true,
    name: "snv_two_genes",
    allele: "C",
    site: "21:7000",
};

/// Default-format rows on ENST1, a GENEA transcript, ENST2, a GENEB transcript, and ENST3,
/// a transcript of a gene no entry names, as many as values are given: each value the row's
/// SpliceAI pairs of the Extra text, or nothing.
fn rows_of(values: &[&str], rows: Rows) -> String {
    let lines: Vec<String> = SPLICEAI_ROW_GENES
        .iter()
        .zip(values)
        .map(|((feature, gene), value)| {
            let mut extra = "IMPACT=MODIFIER".to_string();
            if rows.symbols {
                extra += &format!(";SYMBOL={gene}");
            }
            if !value.is_empty() {
                extra += &format!(";{value}");
            }
            default_row(
                rows.name,
                rows.site,
                rows.allele,
                feature,
                "intron_variant",
                &extra,
            )
        })
        .collect();
    default_text(&lines)
}

fn spliceai_rows(values: &[&str]) -> String {
    rows_of(values, DELETION)
}

fn spliceai_base_rows(values: &[&str]) -> String {
    rows_of(values, BASE)
}

/// What SpliceAI.pm:326-332 make of a record of these symbols' entries: the first symbol
/// and everything after it.
fn joined_value(symbols: &[&str]) -> String {
    let joined = symbols
        .iter()
        .map(|symbol| format!("G|{symbol}|{}", scores(symbol)))
        .collect::<Vec<_>>()
        .join(",");
    joined.split_once('|').unwrap().1.to_string()
}

fn own_value(symbol: &str) -> String {
    format!("{symbol}|{}", scores(symbol))
}

fn pred(value: &str) -> String {
    format!("SpliceAI_pred={value}")
}

/// The nine `SpliceAI_pred_` pairs SpliceAI.pm:276-287 cut from a value, and a flag.
fn split_pairs(value: &str, cutoff: &str) -> String {
    let mut pairs: Vec<String> = SPLICEAI_MEMBERS
        .iter()
        .zip(value.split('|'))
        .map(|(member, part)| format!("SpliceAI_pred_{member}={part}"))
        .collect();
    if !cutoff.is_empty() {
        pairs.push(format!("SpliceAI_cutoff={cutoff}"));
    }
    pairs.join(";")
}

/// The verdict of two entry lists under the class and what was applied.
fn compare_parsed(
    reference: Vec<Entry>,
    candidate: Vec<Entry>,
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(&[SPLICEAI_MULTI]);
    let report = compare_entries(
        &reference,
        &candidate,
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

fn compare_spliceai(
    reference: &str,
    candidate: &str,
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    compare_parsed(
        parse_default(reference).expect("the reference parses"),
        parse_default(candidate).expect("the candidate parses"),
        side,
    )
}

/// A verdict and the rows counted under the class, none for no class applied.
fn outcome(verdict: &str, taken: i64) -> (String, BTreeMap<String, i64>) {
    let applied = if taken == 0 {
        BTreeMap::new()
    } else {
        count_map(&[(SPLICEAI_MULTI, taken)])
    };
    (verdict.to_string(), applied)
}

/// A made-up score record at 21 `position` of one entry per value: the alternate allele
/// before a symbol and its scores (`own_value`).
fn spliceai_record(position: i64, reference: &str, alternate: &str, values: &[&str]) -> String {
    let entries: Vec<String> = values
        .iter()
        .map(|value| format!("{alternate}|{value}"))
        .collect();
    format!(
        "21\t{position}\t.\t{reference}\t{alternate}\t.\t.\tSpliceAI={}",
        entries.join(",")
    )
}

/// A record's value behind the first entry's allele.
fn record_value(record: &str) -> String {
    let info = record.split_once("SpliceAI=").unwrap().1;
    info.split_once('|').unwrap().1.to_string()
}

/// Made-up score files: the file of single-base records and, where given, the file of
/// insertion and deletion records.
fn spliceai_score_files(directory: &Path, snv: &[String], indel: Option<&[String]>) -> Side {
    std::fs::create_dir_all(directory).unwrap();
    let mut paths = vec![write_table(
        directory,
        snv,
        SCORE_HEADER,
        "snv.synthetic.vcf.gz",
    )];
    if let Some(indel) = indel {
        paths.push(write_table(
            directory,
            indel,
            SCORE_HEADER,
            "indel.synthetic.vcf.gz",
        ));
    }
    side_of(opened(&paths, None).expect("the score files open"))
}

#[test]
fn both_rows_of_a_record_of_several_entries_are_taken() {
    let tmp = scratch("spliceai-both-rows");
    let table = spliceai_table(&tmp.join("two"));
    let (whole, first, second) = (
        joined_value(&["GENEA", "GENEB"]),
        own_value("GENEA"),
        own_value("GENEB"),
    );
    let reference = spliceai_rows(&[&pred(&whole), ""]);
    let candidate = spliceai_rows(&[&pred(&first), &pred(&second)]);
    assert_eq!(
        compare_spliceai(&reference, &candidate, &table),
        outcome("green", 2)
    );
    // Rows that state no SYMBOL take their gene from the candidate cache.
    let cache = spliceai_cache(&tmp.join("cache"));
    let cached = spliceai_table_of(
        &tmp.join("cached"),
        &["GENEA", "GENEB"],
        Join::All,
        None,
        Some(&cache),
    );
    let unnamed = (
        rows_of(&[&pred(&whole), ""], UNNAMED),
        rows_of(&[&pred(&first), &pred(&second)], UNNAMED),
    );
    assert_eq!(
        compare_spliceai(&unnamed.0, &unnamed.1, &cached),
        outcome("green", 2)
    );
    let flagged = compare_spliceai(
        &spliceai_rows(&[&format!("{};SpliceAI_cutoff=FAIL", pred(&whole)), ""]),
        &spliceai_rows(&[
            &format!("{};SpliceAI_cutoff=FAIL", pred(&first)),
            &format!("{};SpliceAI_cutoff=PASS", pred(&second)),
        ]),
        &table,
    );
    assert_eq!(flagged, outcome("green", 2));
    let split = compare_spliceai(
        &spliceai_rows(&[&split_pairs(&whole, "FAIL"), ""]),
        &spliceai_rows(&[&split_pairs(&first, "FAIL"), &split_pairs(&second, "FAIL")]),
        &table,
    );
    assert_eq!(split, outcome("green", 2));
    let all_three = ["GENEA", "GENEB", "GENEC"];
    let three = spliceai_table_of(&tmp.join("three"), &all_three, Join::All, None, None);
    assert_eq!(
        compare_spliceai(
            &spliceai_rows(&[&pred(&joined_value(&all_three)), ""]),
            &candidate,
            &three
        ),
        outcome("green", 2)
    );
    // Two records at the one allele, GENEA's and GENEB's entries then GENEC's alone: each
    // gene one entry names is read.
    let pairs = spliceai_table_of(&tmp.join("pairs"), &all_three, Join::Each(2), None, None);
    assert_eq!(
        compare_spliceai(&reference, &candidate, &pairs),
        outcome("green", 2)
    );

    let inputs = vec!["21\t5000\tdel_two_genes\tGCA\tG\t.\t.\t.".to_string()];
    let names: Vec<String> = SPLICEAI_MEMBERS
        .iter()
        .map(|member| format!("SpliceAI_pred_{member}"))
        .collect();
    let header = format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|{}\">\n",
        names.join("|")
    );
    let vcf = |values: &[&str]| -> Vec<Entry> {
        let entries: Vec<String> = SPLICEAI_ROW_GENES
            .iter()
            .zip(values)
            .map(|((feature, _), value)| {
                let cells = if value.is_empty() {
                    "|".repeat(8)
                } else {
                    value.replace(',', "&")
                };
                format!("-|intron_variant|{feature}|Transcript|{cells}")
            })
            .collect();
        let line = &inputs[0][..inputs[0].len() - 1];
        let text = format!("{header}{line}CSQ={}\n", entries.join(","));
        parse_vcf(&text, &inputs).expect("the line parses")
    };
    let cut = whole.split('|').take(9).collect::<Vec<_>>().join("|");
    assert_eq!(
        compare_parsed(vcf(&[&cut, ""]), vcf(&[&first, &second]), &cached),
        outcome("green", 2)
    );
    assert_eq!(
        compare_parsed(vcf(&[&cut, ""]), vcf(&[&first, &second]), &table),
        outcome("red", 0)
    );

    let record = |values: &[&str]| -> Vec<Entry> {
        let text = |value: &str| Value::Str(value.to_string());
        let entries: Vec<Value> = SPLICEAI_ROW_GENES
            .iter()
            .zip(values)
            .map(|((feature, _), value)| {
                let mut entry: IndexMap<String, Value> = IndexMap::from([
                    ("variant_allele".to_string(), text("-")),
                    ("transcript_id".to_string(), text(feature)),
                    (
                        "consequence_terms".to_string(),
                        Value::List(vec![text("intron_variant")]),
                    ),
                ]);
                if !value.is_empty() {
                    let members: IndexMap<String, Value> = SPLICEAI_MEMBERS
                        .iter()
                        .zip(value.split('|'))
                        .map(|(member, part)| {
                            let written = if *member == "SYMBOL" || part.contains(',') {
                                text(part)
                            } else {
                                pyjson::loads(part).expect("a number")
                            };
                            (member.to_string(), written)
                        })
                        .collect();
                    entry.insert("spliceai".to_string(), Value::Dict(members));
                }
                Value::Dict(entry)
            })
            .collect();
        let object = Value::Dict(IndexMap::from([
            ("input".to_string(), text(&inputs[0])),
            ("id".to_string(), text("del_two_genes")),
            ("seq_region_name".to_string(), text("21")),
            ("start".to_string(), Value::Int(5001)),
            ("end".to_string(), Value::Int(5002)),
            ("allele_string".to_string(), text("CA/-")),
            (
                "most_severe_consequence".to_string(),
                text("intron_variant"),
            ),
            ("transcript_consequences".to_string(), Value::List(entries)),
        ]));
        let line = pyjson::dumps(&object, &pyjson::DumpOptions::default()) + "\n";
        parse_json(&line, &inputs).expect("the record parses")
    };
    assert_eq!(
        compare_parsed(record(&[&cut, ""]), record(&[&first, &second]), &cached),
        outcome("green", 2)
    );
}

#[test]
fn nothing_is_taken_without_the_score_files() {
    let tmp = scratch("spliceai-no-files");
    let (whole, first, second) = (
        joined_value(&["GENEA", "GENEB"]),
        own_value("GENEA"),
        own_value("GENEB"),
    );
    let reference = spliceai_rows(&[&pred(&whole), ""]);
    let candidate = spliceai_rows(&[&pred(&first), &pred(&second)]);
    assert_eq!(
        compare_spliceai(&reference, &candidate, &Side::default()),
        outcome("red", 0)
    );
    let table = spliceai_table(&tmp.join("two"));
    assert_eq!(
        compare_spliceai(&reference, &candidate, &table),
        outcome("green", 2)
    );
    // A file that holds no record at the row's variant.
    let elsewhere = (
        reference.replace(SPLICEAI_SITE, "21:6001-6002"),
        candidate.replace(SPLICEAI_SITE, "21:6001-6002"),
    );
    assert_eq!(
        compare_spliceai(&elsewhere.0, &elsewhere.1, &table),
        outcome("red", 0)
    );
    // Neither a SYMBOL nor the candidate cache names the row's gene.
    let unnamed = (
        rows_of(&[&pred(&whole), ""], UNNAMED),
        rows_of(&[&pred(&first), &pred(&second)], UNNAMED),
    );
    assert_eq!(
        compare_spliceai(&unnamed.0, &unnamed.1, &table),
        outcome("red", 0)
    );
}

#[test]
fn a_row_off_the_shape_stays_a_difference() {
    let tmp = scratch("spliceai-off-the-shape");
    let (whole, first, second) = (
        joined_value(&["GENEA", "GENEB"]),
        own_value("GENEA"),
        own_value("GENEB"),
    );
    let joined = spliceai_table(&tmp.join("two"));
    let reference = spliceai_rows(&[&pred(&whole), ""]);
    let candidate = spliceai_rows(&[&pred(&first), &pred(&second)]);
    // A comma-joined value against records of one entry each.
    let single = spliceai_table_of(
        &tmp.join("single"),
        &["GENEA", "GENEB"],
        Join::Nothing,
        None,
        None,
    );
    assert_eq!(
        compare_spliceai(&reference, &candidate, &single),
        outcome("red", 0)
    );
    // On the first entry's gene's row the reference must write the record's whole value.
    assert_eq!(
        compare_spliceai(&spliceai_rows(&["", ""]), &candidate, &joined),
        outcome("red", 1)
    );
    let longer = format!("{},G|GENEC|{}", pred(&whole), scores("GENEC"));
    assert_eq!(
        compare_spliceai(&spliceai_rows(&[&longer, ""]), &candidate, &joined),
        outcome("red", 1)
    );
    // The candidate's value must be the entry of the row's own gene.
    for other in [
        pred(&first),
        pred(&second.replace("0.02", "0.03")),
        pred(&own_value("GENEC")),
    ] {
        assert_eq!(
            compare_spliceai(
                &reference,
                &spliceai_rows(&[&pred(&first), &other]),
                &joined
            ),
            outcome("red", 1),
            "{other}"
        );
    }
    // A row whose stated gene is not the value's.
    let stated = candidate.replace("SYMBOL=GENEB", "SYMBOL=GENEC");
    assert_eq!(
        compare_spliceai(
            &reference.replace("SYMBOL=GENEB", "SYMBOL=GENEC"),
            &stated,
            &joined
        ),
        outcome("red", 1)
    );
    // The cache names each row's gene: the second gene's row writes the third entry, and
    // the row of a gene no entry names writes the second.
    let all_three = ["GENEA", "GENEB", "GENEC"];
    let cache = spliceai_cache(&tmp.join("cache"));
    let three = spliceai_table_of(
        &tmp.join("three"),
        &all_three,
        Join::All,
        None,
        Some(&cache),
    );
    let whole_three = pred(&joined_value(&all_three));
    assert_eq!(
        compare_spliceai(
            &rows_of(&[&whole_three, ""], UNNAMED),
            &rows_of(&[&pred(&first), &pred(&own_value("GENEC"))], UNNAMED),
            &three
        ),
        outcome("red", 1)
    );
    assert_eq!(
        compare_spliceai(
            &rows_of(&[&whole_three, "", ""], UNNAMED),
            &rows_of(&[&pred(&first), &pred(&second), &pred(&second)], UNNAMED),
            &three
        ),
        outcome("red", 2)
    );
    // The two flags of the first gene's row agree, and the reference writes none on the
    // second's.
    assert_eq!(
        compare_spliceai(
            &spliceai_rows(&[&format!("{};SpliceAI_cutoff=FAIL", pred(&whole)), ""]),
            &spliceai_rows(&[
                &format!("{};SpliceAI_cutoff=PASS", pred(&first)),
                &pred(&second)
            ]),
            &joined
        ),
        outcome("red", 1)
    );
    assert_eq!(
        compare_spliceai(
            &spliceai_rows(&[&pred(&whole), "SpliceAI_cutoff=FAIL"]),
            &spliceai_rows(&[
                &pred(&first),
                &format!("{};SpliceAI_cutoff=PASS", pred(&second))
            ]),
            &joined
        ),
        outcome("red", 1)
    );
    // Another allele at the site.
    let other_allele = Rows {
        allele: "T",
        ..DELETION
    };
    assert_eq!(
        compare_spliceai(
            &rows_of(&[&pred(&whole), ""], other_allele),
            &rows_of(&[&pred(&first), &pred(&second)], other_allele),
            &joined
        ),
        outcome("red", 0)
    );
    // A record whose INFO column holds a second key.
    let info = format!(
        "SpliceAI=G|GENEA|{},G|GENEB|{};DB",
        scores("GENEA"),
        scores("GENEB")
    );
    let second_key = spliceai_table_of(
        &tmp.join("second-key"),
        &["GENEA"],
        Join::Nothing,
        Some(&info),
        None,
    );
    assert_eq!(
        compare_spliceai(&reference, &candidate, &second_key),
        outcome("red", 0)
    );
}

#[test]
fn the_rows_of_a_single_base_record_of_several_entries_are_taken() {
    let tmp = scratch("spliceai-single-base");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let candidate = spliceai_base_rows(&[&pred(&first), &pred(&second)]);
    let two = spliceai_record(7000, "A", "C", &[&first, &second]);
    let reference = spliceai_base_rows(&[&pred(&record_value(&two)), ""]);
    let files = spliceai_score_files(&tmp.join("two"), std::slice::from_ref(&two), None);
    assert_eq!(
        compare_spliceai(&reference, &candidate, &files),
        outcome("green", 2)
    );
    // A record of nine entries with GENEA's first and GENEB's fifth.
    let others: Vec<String> = "DEFGHIJ"
        .chars()
        .map(|letter| format!("GENE{letter}|{}", scores("GENEC")))
        .collect();
    let mut entries: Vec<&str> = vec![&first];
    entries.extend(others[..3].iter().map(String::as_str));
    entries.push(&second);
    entries.extend(others[3..].iter().map(String::as_str));
    let nine = spliceai_record(7000, "A", "C", &entries);
    assert_eq!(record_value(&nine).matches(',').count(), 8);
    let whole_of_nine = spliceai_base_rows(&[&pred(&record_value(&nine)), ""]);
    let files = spliceai_score_files(&tmp.join("nine"), &[nine], None);
    assert_eq!(
        compare_spliceai(&whole_of_nine, &candidate, &files),
        outcome("green", 2)
    );
    // The records one base before the row's Location and the insertions and deletions
    // anchored there stand beside the record and leave both rows taken.
    let before = spliceai_record(6999, "T", "C", &[&first, &second]);
    let indels: Vec<String> = [
        (6999, "T", "TC"),
        (6999, "TA", "T"),
        (7000, "A", "AC"),
        (7000, "AG", "A"),
    ]
    .into_iter()
    .map(|(position, reference, alternate)| {
        spliceai_record(position, reference, alternate, &[&first, &second])
    })
    .collect();
    let beside = spliceai_score_files(&tmp.join("beside"), &[before, two], Some(&indels));
    assert_eq!(
        table_of(&beside)
            .records("21", 6999, 7000)
            .expect("the files read")
            .len(),
        6
    );
    assert_eq!(
        compare_spliceai(&reference, &candidate, &beside),
        outcome("green", 2)
    );
}

#[test]
fn a_single_base_row_stays_a_difference_without_the_score_files_and_off_the_shape() {
    let tmp = scratch("spliceai-single-base-refused");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let two = spliceai_record(7000, "A", "C", &[&first, &second]);
    let reference = spliceai_base_rows(&[&pred(&record_value(&two)), ""]);
    let candidate = spliceai_base_rows(&[&pred(&first), &pred(&second)]);
    assert_eq!(
        compare_spliceai(&reference, &candidate, &Side::default()),
        outcome("red", 0)
    );
    let single = spliceai_score_files(
        &tmp.join("single"),
        &[
            spliceai_record(7000, "A", "C", &[&first]),
            spliceai_record(7000, "A", "C", &[&second]),
        ],
        None,
    );
    assert_eq!(
        compare_spliceai(&reference, &candidate, &single),
        outcome("red", 0)
    );
    let joined = spliceai_score_files(&tmp.join("two"), &[two], None);
    assert_eq!(
        compare_spliceai(
            &reference,
            &spliceai_base_rows(&[&pred(&first), &pred(&first)]),
            &joined
        ),
        outcome("red", 1)
    );
    let unheld = format!("{},C|{}", pred(&first), own_value("GENEC"));
    assert_eq!(
        compare_spliceai(&spliceai_base_rows(&[&unheld, ""]), &candidate, &joined),
        outcome("red", 1)
    );
}

/// Made-up entries of GENEA with other scores than `scores` gives it: the first fails a
/// cutoff of 0.5 and the second reaches it.
const GENEA_LOW: &str = "GENEA|0.11|0.21|0.31|0.41|-1|-2|-3|-4";
const GENEA_HIGH: &str = "GENEA|0.12|0.92|0.32|0.42|9|8|7|6";

/// The one-field value of a row whose gene the `entries` of one record name: the first
/// behind its allele and each later one with the allele `C`.
fn entries_value(entries: &[&str]) -> String {
    entries.join(",C|")
}

/// The nine `SpliceAI_pred_` pairs of a row whose gene the `entries` of one record name,
/// each member its values over the entries joined by commas, and a flag.
fn entries_pairs(entries: &[&str], cutoff: &str) -> String {
    let parts: Vec<Vec<&str>> = entries
        .iter()
        .map(|entry| entry.split('|').collect())
        .collect();
    let mut pairs: Vec<String> = SPLICEAI_MEMBERS
        .iter()
        .enumerate()
        .map(|(index, member)| {
            let values: Vec<&str> = parts.iter().map(|entry| entry[index]).collect();
            format!("SpliceAI_pred_{member}={}", values.join(","))
        })
        .collect();
    if !cutoff.is_empty() {
        pairs.push(format!("SpliceAI_cutoff={cutoff}"));
    }
    pairs.join(";")
}

/// The six kinds of record in which several entries name one gene: the candidate's row
/// carries every entry of the record that names its gene, in the record's order, and the
/// reference's the record's whole value where the gene is the first entry's and nothing
/// where it is a later entry's, in the one-field and in the nine-field form.
#[test]
fn the_rows_of_a_gene_several_entries_of_one_record_name_are_taken() {
    let tmp = scratch("spliceai-one-gene-entries");
    let (first, other) = (own_value("GENEA"), own_value("GENEB"));
    let kinds: [(&str, Vec<&str>); 6] = [
        ("equal", vec![&first, &first]),
        ("differ", vec![GENEA_LOW, GENEA_HIGH]),
        ("four", vec![GENEA_LOW, &first, GENEA_HIGH, &first]),
        ("twice-then-other", vec![&first, GENEA_LOW, &other]),
        ("around-other", vec![GENEA_LOW, &other, &first]),
        ("other-then-twice", vec![&other, GENEA_LOW, GENEA_HIGH]),
    ];
    for (kind, entries) in kinds {
        let record = spliceai_record(7000, "A", "C", &entries);
        let files = spliceai_score_files(&tmp.join(kind), std::slice::from_ref(&record), None);
        let of = |gene: &str| -> Vec<&str> {
            entries
                .iter()
                .copied()
                .filter(|entry| entry.starts_with(gene))
                .collect()
        };
        let (own, others) = (of("GENEA"), of("GENEB"));
        let leads = entries[0].starts_with("GENEA");
        let whole = record_value(&record);

        // The one field: a record that names GENEA alone gives both engines one value.
        let candidate_others = if others.is_empty() {
            String::new()
        } else {
            pred(&entries_value(&others))
        };
        let candidate = spliceai_base_rows(&[&pred(&entries_value(&own)), &candidate_others]);
        let reference = if leads {
            spliceai_base_rows(&[&pred(&whole), ""])
        } else {
            spliceai_base_rows(&["", &pred(&whole)])
        };
        let differing = if others.is_empty() { 0 } else { 2 };
        assert_eq!(
            compare_spliceai(&reference, &candidate, &files),
            outcome("green", differing),
            "{kind}: one field"
        );

        // The nine fields: every member of a row of several entries differs but one whose
        // first entry writes the others' values too, and so do the two flags where a later
        // entry of the gene alone reaches the cutoff.
        let largest = |entry: &str| -> f64 {
            entry.split('|').collect::<Vec<_>>()[1..5]
                .iter()
                .map(|score| score.parse::<f64>().unwrap())
                .fold(0.0, f64::max)
        };
        let verdict_of = |reaches: bool| if reaches { "PASS" } else { "FAIL" };
        let ours_flag = verdict_of(own.iter().any(|entry| largest(entry) >= 0.5));
        let theirs_flag = verdict_of(largest(entries[0]) >= 0.5);
        let other_flag = verdict_of(others.iter().any(|entry| largest(entry) >= 0.5));
        let candidate_others = if others.is_empty() {
            String::new()
        } else {
            entries_pairs(&others, other_flag)
        };
        let candidate = spliceai_base_rows(&[&entries_pairs(&own, ours_flag), &candidate_others]);
        let cut = split_pairs(&whole, theirs_flag);
        let reference = if leads {
            spliceai_base_rows(&[&cut, ""])
        } else {
            spliceai_base_rows(&["", &cut])
        };
        let (verdict, applied) = compare_spliceai(&reference, &candidate, &files);
        assert_eq!(verdict, "green", "{kind}: nine fields");
        assert!(
            applied.get(SPLICEAI_MULTI).copied().unwrap_or(0) >= 1,
            "{kind}: {applied:?}"
        );
    }
}

/// A row of a gene several entries of one record name stays a difference off the shape:
/// a candidate that writes the later entry alone, the first alone or the entries in
/// another order, and a pair of flags no one cutoff gives.
#[test]
fn a_row_of_a_gene_several_entries_name_stays_a_difference_off_the_shape() {
    let tmp = scratch("spliceai-one-gene-entries-refused");
    let other = own_value("GENEB");
    let record = spliceai_record(7000, "A", "C", &[GENEA_LOW, &other, GENEA_HIGH]);
    let files = spliceai_score_files(&tmp.join("around"), std::slice::from_ref(&record), None);
    let whole = record_value(&record);
    let reference = spliceai_base_rows(&[&pred(&whole), ""]);
    let others = pred(&other);
    for (what, value) in [
        ("the later entry alone", GENEA_HIGH.to_string()),
        ("the first entry alone", GENEA_LOW.to_string()),
        (
            "the entries in the other order",
            entries_value(&[GENEA_HIGH, GENEA_LOW]),
        ),
    ] {
        let candidate = spliceai_base_rows(&[&pred(&value), &others]);
        let (verdict, applied) = compare_spliceai(&reference, &candidate, &files);
        assert_eq!(verdict, "red", "{what}");
        assert_eq!(applied, count_map(&[(SPLICEAI_MULTI, 1)]), "{what}");
    }

    // Flags: the first entry fails 0.5 and a later entry of the gene reaches it, so FAIL
    // beside PASS stands and PASS beside FAIL does not; over entries none of which passes
    // the first's largest score, no cutoff parts the two flags.
    let own = [GENEA_LOW, GENEA_HIGH];
    let cut = |flag: &str| spliceai_base_rows(&[&split_pairs(&whole, flag), ""]);
    let joined = |flag: &str| {
        spliceai_base_rows(&[
            &entries_pairs(&own, flag),
            &entries_pairs(&[&other], "FAIL"),
        ])
    };
    assert_eq!(
        compare_spliceai(&cut("FAIL"), &joined("PASS"), &files).0,
        "green"
    );
    assert_eq!(
        compare_spliceai(&cut("PASS"), &joined("FAIL"), &files).0,
        "red"
    );
    let level = spliceai_record(7000, "A", "C", &[GENEA_HIGH, GENEA_LOW]);
    let level_files = spliceai_score_files(&tmp.join("level"), std::slice::from_ref(&level), None);
    let level_whole = record_value(&level);
    let level_own = [GENEA_HIGH, GENEA_LOW];
    let level_cut = |flag: &str| spliceai_base_rows(&[&split_pairs(&level_whole, flag)]);
    let level_joined = |flag: &str| spliceai_base_rows(&[&entries_pairs(&level_own, flag)]);
    assert_eq!(
        compare_spliceai(&level_cut("PASS"), &level_joined("PASS"), &level_files).0,
        "green"
    );
    assert_eq!(
        compare_spliceai(&level_cut("FAIL"), &level_joined("PASS"), &level_files).0,
        "red"
    );
}

#[test]
fn every_record_line_at_the_allele_is_read() {
    let tmp = scratch("spliceai-every-line");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let changed = format!("GENEA|{}", scores("GENEC"));
    let other = format!("GENEB|{}", scores("GENEA"));
    let third = format!("GENEA|{}", scores("GENEB"));
    let last = format!("GENEB|{}", scores("GENEC"));
    let both: Vec<&str> = vec![&first, &second];
    let reversed: Vec<&str> = vec![&second, &first];
    let whole = |values: &[&str]| pred(&record_value(&spliceai_record(7000, "A", "C", values)));
    let own = (pred(&first), pred(&second));
    let none = String::new();
    // A name, the record lines, the reference's values, the candidate's and the outcome.
    type Shape<'a> = (
        &'a str,
        Vec<Vec<&'a str>>,
        Vec<String>,
        Vec<String>,
        (&'a str, i64),
    );
    let shapes: Vec<Shape> = vec![
        (
            "one-line",
            vec![both.clone()],
            vec![whole(&both), none.clone()],
            vec![own.0.clone(), own.1.clone()],
            ("green", 2),
        ),
        (
            "twice",
            vec![both.clone(), both.clone()],
            vec![whole(&both), none.clone()],
            vec![own.0.clone(), own.1.clone()],
            ("green", 2),
        ),
        (
            "three-times",
            vec![both.clone(), both.clone(), both.clone()],
            vec![whole(&both), none.clone()],
            vec![own.0.clone(), own.1.clone()],
            ("green", 2),
        ),
        (
            "other-order",
            vec![both.clone(), reversed.clone()],
            vec![whole(&both), whole(&reversed)],
            vec![own.0.clone(), own.1.clone()],
            ("green", 2),
        ),
        (
            "one-entry-last",
            vec![both.clone(), vec![&first]],
            vec![whole(&[&first]), none.clone()],
            vec![own.0.clone(), own.1.clone()],
            ("green", 1),
        ),
        (
            "one-entry-first",
            vec![vec![&first], both.clone()],
            vec![whole(&both), none.clone()],
            vec![own.0.clone(), own.1.clone()],
            ("green", 2),
        ),
        (
            "one-gene-twice",
            vec![vec![&first, &first]],
            vec![whole(&[&first, &first]), none.clone()],
            vec![own.0.clone(), none.clone()],
            ("red", 0),
        ),
        (
            "gene-unnamed",
            vec![both.clone(), both.clone()],
            vec![whole(&both), none.clone(), none.clone()],
            vec![own.0.clone(), own.1.clone(), own.1.clone()],
            ("red", 2),
        ),
        (
            "other-scores",
            vec![both.clone(), vec![&changed, &other]],
            vec![whole(&[&changed, &other]), none.clone()],
            vec![pred(&changed), pred(&other)],
            ("green", 2),
        ),
        (
            "other-scores-three-times",
            vec![both.clone(), vec![&changed, &other], vec![&third, &last]],
            vec![whole(&[&third, &last]), none.clone()],
            vec![pred(&third), pred(&last)],
            ("green", 2),
        ),
        (
            "other-scores-other-order",
            vec![both.clone(), vec![&other, &changed]],
            vec![whole(&both), whole(&[&other, &changed])],
            vec![pred(&changed), pred(&other)],
            ("green", 2),
        ),
        (
            "earlier-line",
            vec![both.clone(), vec![&changed, &other]],
            vec![whole(&[&changed, &other]), none.clone()],
            vec![own.0.clone(), pred(&other)],
            ("red", 1),
        ),
        (
            "later-entry-of-a-line",
            vec![vec![&first, &changed]],
            vec![whole(&[&first, &changed]), none.clone()],
            vec![pred(&changed), none.clone()],
            ("red", 0),
        ),
        (
            "earlier-entry-of-a-line",
            vec![vec![&first, &changed]],
            vec![whole(&[&first, &changed]), none.clone()],
            vec![own.0.clone(), none.clone()],
            ("red", 0),
        ),
    ];
    for (name, lines, reference, candidate, (verdict, taken)) in shapes {
        let records: Vec<String> = lines
            .iter()
            .map(|line| spliceai_record(7000, "A", "C", line))
            .collect();
        let table = spliceai_score_files(&tmp.join(name), &records, None);
        let values = |written: &[String]| -> String {
            spliceai_base_rows(&written.iter().map(String::as_str).collect::<Vec<_>>())
        };
        let rows = (values(&reference), values(&candidate));
        assert_eq!(
            compare_spliceai(&rows.0, &rows.1, &table),
            outcome(verdict, taken),
            "{name}"
        );
        assert_eq!(
            compare_spliceai(&rows.0, &rows.1, &Side::default()),
            outcome("red", 0),
            "{name}"
        );
    }
    // The flags: equal on the row of the gene that leads a line, the candidate's alone on
    // the other row.
    let line = spliceai_record(7000, "A", "C", &both);
    let twice = spliceai_score_files(&tmp.join("flagged"), &[line.clone(), line], None);
    let flagged = spliceai_base_rows(&[&format!("{};SpliceAI_cutoff=FAIL", whole(&both)), ""]);
    assert_eq!(
        compare_spliceai(
            &flagged,
            &spliceai_base_rows(&[
                &format!("{};SpliceAI_cutoff=FAIL", own.0),
                &format!("{};SpliceAI_cutoff=PASS", own.1)
            ]),
            &twice
        ),
        outcome("green", 2)
    );
    assert_eq!(
        compare_spliceai(
            &flagged,
            &spliceai_base_rows(&[
                &format!("{};SpliceAI_cutoff=PASS", own.0),
                &format!("{};SpliceAI_cutoff=PASS", own.1)
            ]),
            &twice
        ),
        outcome("red", 1)
    );
}

#[test]
fn the_rows_refused_are_counted_by_ground() {
    let tmp = scratch("spliceai-grounds");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let changed = format!("GENEA|{}", scores("GENEC"));
    let whole = format!("{},C|{second}", pred(&first));
    let own = pred(&first);
    let record = |position: i64, entries: &[&str]| spliceai_record(position, "A", "C", entries);
    let both: [&str; 2] = [&first, &second];
    // A Location, the gene a row states, the record lines of the site, the reference's
    // value and the candidate's.
    type Case = (&'static str, &'static str, Vec<String>, String, String);
    let cases: Vec<Case> = vec![
        (
            "21:7410",
            "GENEA",
            vec![record(7410, &both)],
            whole.clone(),
            String::new(),
        ),
        ("21", "GENEA", vec![], whole.clone(), own.clone()),
        (
            "21:7430",
            "GENEA",
            vec![record(7430, &both) + ";DB"],
            whole.clone(),
            own.clone(),
        ),
        (
            "21:7440",
            "GENEA",
            vec![record(7440, &both), spliceai_record(7440, "G", "C", &both)],
            whole.clone(),
            own.clone(),
        ),
        ("21:7450", "GENEA", vec![], whole.clone(), own.clone()),
        ("21:7455", "GENEA", vec![], whole.clone(), own.clone()),
        (
            "21:7460",
            "",
            vec![record(7460, &both)],
            whole.clone(),
            own.clone(),
        ),
        (
            "21:7470",
            "GENEX",
            vec![record(7470, &both)],
            String::new(),
            pred(&second),
        ),
        (
            "21:7490",
            "GENEA",
            vec![record(7490, &[&first])],
            own.clone(),
            pred(&changed),
        ),
        (
            "21:7500",
            "GENEA",
            vec![record(7500, &both), record(7500, &[&changed, &second])],
            format!("{},C|{second}", pred(&changed)),
            own.clone(),
        ),
        (
            "21:7510",
            "GENEA",
            vec![record(7510, &[&second, &first])],
            format!("{};SpliceAI_cutoff=FAIL", pred(&second)),
            format!("{own};SpliceAI_cutoff=PASS"),
        ),
        (
            "21:7520",
            "GENEA",
            vec![record(7520, &both)],
            pred(&second),
            own.clone(),
        ),
        (
            "21:7530",
            "GENEA",
            vec![record(7530, &both)],
            format!("{whole};SpliceAI_cutoff=FAIL"),
            format!("{own};SpliceAI_cutoff=PASS"),
        ),
        (
            "21:7550",
            "GENEA",
            vec![record(7550, &both)],
            whole.clone(),
            own.clone(),
        ),
        ("21:7550", "GENEB", vec![], String::new(), pred(&second)),
    ];
    let rows = |candidate: bool, chosen: &[Case]| -> Vec<Entry> {
        let lines: Vec<String> = chosen
            .iter()
            .map(|(location, gene, _, reference, written)| {
                let written = if candidate { written } else { reference };
                let mut extra = "IMPACT=MODIFIER".to_string();
                if !gene.is_empty() {
                    extra += &format!(";SYMBOL={gene}");
                }
                if !written.is_empty() {
                    extra += &format!(";{written}");
                }
                default_row(
                    &format!("snv_{location}_{gene}"),
                    location,
                    "C",
                    &format!("ENST_{gene}"),
                    "intron_variant",
                    &extra,
                )
            })
            .collect();
        parse_default(&default_text(&lines)).expect("the rows parse")
    };
    let lines: Vec<String> = cases
        .iter()
        .flat_map(|(_, _, lines, _, _)| lines.iter().cloned())
        .collect();
    let table = spliceai_score_files(&tmp.join("scores"), &lines, None);
    let refused = [
        ("no_candidate_value", 1),
        ("no_span", 1),
        ("other_info", 1),
        ("several_references", 1),
        ("no_record", 2),
        ("no_gene", 1),
        ("gene_not_named", 1),
        ("lone_entry", 1),
        ("candidate_value", 1),
        ("reference_not_silent", 1),
        ("reference_value", 1),
        ("reference_flag", 1),
    ];
    assert!(refused
        .iter()
        .map(|(ground, _)| ground)
        .eq(SPLICEAI_GROUNDS.iter().map(|(ground, _)| ground)));
    let compared = |expected: &[Entry], actual: &[Entry], side: &Side| -> Value {
        let mut documented = documented_with(&[SPLICEAI_MULTI]);
        compare_entries(
            expected,
            actual,
            &mut documented,
            &CompareOptions::new(side),
        )
        .expect("the entries compare")
    };
    let report = compared(&rows(false, &cases), &rows(true, &cases), &table);
    let grounds: Vec<String> = refused
        .iter()
        .map(|(ground, rows)| format!("\"{ground}\": {rows}"))
        .collect();
    let dumped = |value: &Value| pyjson::dumps(value, &pyjson::DumpOptions::default());
    assert_eq!(
        dumped(at(&report, &["spliceai_refusals"])),
        format!(
            "{{\"asked\": 15, \"refused\": 13, \"grounds\": {{{}}}}}",
            grounds.join(", ")
        )
    );
    let lines = spliceai_refusal_lines(at(&report, &["spliceai_refusals"]), "");
    assert_eq!(
        lines[0],
        format!("{SPLICEAI_MULTI} under --spliceai-table: rows asked about 15, rows refused 13")
    );
    assert_eq!(
        lines[5],
        format!("{SPLICEAI_MULTI}: rows refused 2: no record line of the score files stands at the row's variant")
    );
    assert_eq!(lines.len(), 13);
    for word in ["GENE", "ENST", "74", "75", "0.1", "snv_"] {
        assert!(!lines.iter().any(|line| line.contains(word)), "{word}");
    }
    // A comparison counts its own rows, and reports no counts where it asks the class
    // about no row or runs without the score files.
    let last_two = &cases[cases.len() - 2..];
    let taken = compared(&rows(false, last_two), &rows(true, last_two), &table);
    assert_eq!(verdict(&taken), "green");
    assert_eq!(
        dumped(at(&taken, &["spliceai_refusals"])),
        "{\"asked\": 2, \"refused\": 0, \"grounds\": {}}"
    );
    let same = compared(&rows(false, &cases), &rows(false, &cases), &table);
    assert!(same.get("spliceai_refusals").is_none());
    let without = compared(&rows(false, &cases), &rows(true, &cases), &Side::default());
    assert!(without.get("spliceai_refusals").is_none());
}

#[test]
fn the_table_reads_the_records_at_a_span_and_refuses_what_is_no_score_file() {
    let tmp = scratch("spliceai-table");
    let paths = indel_file(&tmp.join("two"), &["GENEA", "GENEB"], Join::All, None);
    let table = opened(&paths, None).expect("the score file opens");
    let record = (
        5000,
        "GCA".to_string(),
        "G".to_string(),
        format!(
            "SpliceAI=G|GENEA|{},G|GENEB|{}",
            scores("GENEA"),
            scores("GENEB")
        ),
    );
    let records = |chromosome: &str, low: i64, high: i64| {
        table
            .records(chromosome, low, high)
            .expect("the file reads")
    };
    assert_eq!(records("21", 4999, 5001), std::slice::from_ref(&record));
    assert_eq!(records("chr21", 5000, 5000), [record]);
    assert!(records("21", 5001, 5002).is_empty());
    assert!(records("22", 5000, 5000).is_empty());
    let mut both = table.paths.clone();
    both.extend(indel_file(
        &tmp.join("single"),
        &["GENEA", "GENEB"],
        Join::Nothing,
        None,
    ));
    let both = opened(&both, None).expect("the score files open");
    let entries: Vec<usize> = both
        .records("21", 5000, 5000)
        .expect("the files read")
        .iter()
        .map(|(_, _, _, info)| info.matches(',').count())
        .collect();
    assert_eq!(entries, [1, 0, 0]);
    // A refusal under the option's own name for a file that is no VCF or has no index.
    let seven = write_table(
        &tmp,
        &["21\t5000\tA\tT\tV\tL\t0.300".to_string()],
        "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
        "seven.tsv.gz",
    );
    let refused = opened(&[seven], None).unwrap_err().to_string();
    assert!(refused.starts_with("--spliceai-table: "), "{refused}");
    assert!(
        refused.ends_with("seven.tsv.gz is not a VCF (a header line of eight columns)"),
        "{refused}"
    );
    std::fs::remove_file(tabix_index_of(&table.paths[0])).unwrap();
    let refused = opened(&table.paths, None).unwrap_err().to_string();
    assert!(
        refused.starts_with("--spliceai-table: cannot read the tabix index "),
        "{refused}"
    );
    assert!(refused.contains("indel.synthetic.vcf.gz.tbi"), "{refused}");
}

#[test]
fn a_row_past_the_widest_whole_number_stands_at_no_record() {
    let tmp = scratch("spliceai-wide-position");
    let table = spliceai_table(&tmp.join("two"));
    let row = |location: &str| -> Entry {
        let text = default_text(&[default_row(
            "del_two_genes",
            location,
            "-",
            "ENST1",
            "intron_variant",
            "IMPACT=MODIFIER",
        )]);
        parse_default(&text).expect("the row parses").remove(0)
    };
    let read = |location: &str| {
        let (records, ground) =
            spliceai_entries_at(&row(location), table_of(&table)).expect("the file reads");
        (records.len(), ground)
    };
    assert_eq!(read(SPLICEAI_SITE), (1, None));
    assert_eq!(read("21:5002-5001"), (1, None));
    assert_eq!(read("21:5001"), (0, Some("no_record".to_string())));
    assert_eq!(
        read("21:99999999999999999999"),
        (0, Some("no_record".to_string()))
    );
    assert_eq!(
        read("21:5001-99999999999999999999"),
        (0, Some("no_record".to_string()))
    );
    assert_eq!(read("21:5001-"), (0, Some("no_span".to_string())));
}

#[test]
fn a_cached_gene_symbol_that_is_no_text_names_no_entry() {
    let tmp = scratch("spliceai-symbol-no-text");
    let cache = tmp.join("cache");
    let shards = cache.join("transcripts").join("21");
    std::fs::create_dir_all(&shards).unwrap();
    std::fs::write(
        shards.join("1-1000000.json"),
        "[{\"stable_id\": \"ENST1\", \"gene_symbol\": 5}]",
    )
    .unwrap();
    let cached = spliceai_table_of(
        &tmp.join("two"),
        &["GENEA", "GENEB"],
        Join::All,
        None,
        Some(&cache),
    );
    let (whole, first) = (joined_value(&["GENEA", "GENEB"]), own_value("GENEA"));
    let reference = parse_default(&rows_of(&[&pred(&whole)], UNNAMED)).expect("the row parses");
    let candidate = parse_default(&rows_of(&[&pred(&first)], UNNAMED)).expect("the row parses");
    let ground = |feature: &str| {
        let renamed = |entry: &Entry| Entry {
            feature: feature.to_string(),
            ..entry.clone()
        };
        spliceai_multi_entry_cell_ground(
            &renamed(&reference[0]),
            Some(&renamed(&candidate[0])),
            "SpliceAI_pred",
            &pred(&whole),
            &pred(&first),
            table_of(&cached),
        )
        .expect("the files read")
    };
    assert_eq!(ground("ENST1"), Some("gene_not_named".to_string()));
    assert_eq!(ground("ENST2"), Some("no_gene".to_string()));
}
