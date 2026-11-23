// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the passes over small outputs under a planted resolver and reads what they hand on:
//! the counts of every input and engine and the files of the remainder. One test per
//! extraction rule, the planted missing, extra and changed values, the rows a class reads,
//! and the same bytes on any number of threads.

use flate2::write::GzEncoder;
use flate2::Compression;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use vep_plugin_pairs::pair::{Counts, PROJECTIONS};
use vep_plugin_pairs::resolve::{Answer, Registry, Resolver, Table};
use vep_plugin_pairs::rules::{Refusal, Res};
use vep_plugin_pairs::score::pass::{self, InputFacts, Options, Passed};

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path =
            std::env::temp_dir().join(format!("vep-plugin-pairs-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn file(&self, name: &str, content: impl AsRef<[u8]>) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, content).unwrap();
        path
    }

    /// A default-format output: the column header and `rows`.
    fn output(&self, name: &str, rows: &[String]) -> PathBuf {
        self.file(name, format!("{HEADER}{}", rows.concat()))
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn row(location: &str, allele: &str, feature: &str, extra: &str) -> String {
    format!("v\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\t{extra}\n")
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// A planted registry: `PLANT_*` is the plugin Planted's, `GWAS_*` the plugin GWAS's,
/// `IMPACT` Ensembl VEP's own, `Refused` a field the registry refuses, and any other key no
/// source's. A Location it is asked to settle is the Location as written, without its
/// surrounding spaces. Every question is kept, as `field <hex>` or `location <hex>`.
#[derive(Default)]
struct Planted {
    asked: Mutex<Vec<String>>,
}

impl Resolver for Planted {
    fn field(&self, name: &[u8]) -> Res<Answer> {
        let asked = format!("field {}", hex(name));
        self.asked.lock().unwrap().push(asked);
        match name {
            name if name.starts_with(b"PLANT_") => Ok(Answer::Source("Planted".to_string())),
            name if name.starts_with(b"GWAS_") => Ok(Answer::Source("GWAS".to_string())),
            b"IMPACT" => Ok(Answer::Core),
            b"Refused" => Err(Refusal(format!(
                "the registry refuses the field {}",
                hex(name)
            ))),
            _ => Ok(Answer::Unclaimed),
        }
    }

    fn location(&self, written: &[u8]) -> Res<Vec<u8>> {
        let asked = format!("location {}", hex(written));
        self.asked.lock().unwrap().push(asked);
        Ok(String::from_utf8_lossy(written).trim().as_bytes().to_vec())
    }
}

/// A JSON string of a source's name.
fn quoted(text: &str) -> String {
    let mut out = String::from("\"");
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            control if control < ' ' => write!(out, "\\u{:04x}", control as u32).unwrap(),
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

fn object(entries: Vec<(String, String)>) -> String {
    let mut entries = entries;
    entries.sort();
    let body: Vec<String> = entries
        .into_iter()
        .map(|(key, value)| format!("{key}: {value}"))
        .collect();
    format!("{{{}}}", body.join(", "))
}

fn stated_input(facts: &InputFacts) -> String {
    let unclaimed = facts
        .extracted
        .unclaimed
        .iter()
        .map(|(name, times)| (quoted(&hex(name)), times.to_string()))
        .collect();
    format!(
        "{{\"kind\": {}, \"sha256\": {}, \"lines\": {}, \"outside\": {}, \"header\": {}, \"unclaimed\": {}}}",
        quoted(facts.kind.name()),
        quoted(&facts.sha256),
        facts.extracted.lines,
        facts.extracted.outside,
        quoted(&hex(&facts.extracted.header)),
        object(unclaimed),
    )
}

fn stated_projection(table: &Table, counts: &Counts) -> String {
    let mut values = Vec::new();
    for (field, name) in table.names.iter().enumerate() {
        let (reference, engine) = (counts.values[0][field], counts.values[1][field]);
        if reference > 0 || engine > 0 {
            let stated = format!("[{reference}, {engine}, {}]", counts.paired[field]);
            values.push((quoted(&hex(name)), stated));
        }
    }
    let shared = table
        .sources
        .iter()
        .zip(&counts.shared)
        .filter(|(_, &keys)| keys > 0)
        .map(|(source, keys)| (quoted(source), keys.to_string()))
        .collect();
    format!(
        "{{\"values\": {}, \"shared_keys\": {}, \"remainder\": [{}, {}]}}",
        object(values),
        object(shared),
        counts.left[0],
        counts.left[1],
    )
}

/// What the passes counted, as one text the tests read: per input its kind, the sha256 of
/// the file as stored, its rows or tuple lines, the tuples outside every source, the column
/// header and the unclaimed keys with their counts; per engine and projection each field's
/// values on the reference's side and on the engine's with the pairs equal as text, the
/// shared keys per source and the remainder's lines per side, then the rows written per
/// side. Field names and the header are hexadecimal, since an input can carry bytes no
/// text holds; every object lists its keys in ascending order.
fn stated(passed: &Passed) -> String {
    let mut text = String::from("{\n  \"inputs\": [\n");
    let inputs: Vec<String> = passed
        .inputs
        .iter()
        .map(|facts| format!("    {}", stated_input(facts)))
        .collect();
    text.push_str(&inputs.join(",\n"));
    text.push_str("\n  ],\n  \"engines\": [\n");
    let engines: Vec<String> = passed
        .engines
        .iter()
        .map(|facts| {
            let mut entries: Vec<String> = PROJECTIONS
                .iter()
                .zip(&facts.counts)
                .map(|(name, counts)| {
                    format!(
                        "{}: {}",
                        quoted(name),
                        stated_projection(&facts.table, counts)
                    )
                })
                .collect();
            entries.push(format!("\"rows\": [{}, {}]", facts.rows[0], facts.rows[1]));
            format!("    {{{}}}", entries.join(", "))
        })
        .collect();
    text.push_str(&engines.join(",\n"));
    text.push_str("\n  ]\n}\n");
    text
}

/// One run of the passes: what they counted or the message they refused with, the
/// directory of the remainder's files and the questions the resolver was asked.
struct Done {
    outcome: Result<Passed, String>,
    out: PathBuf,
    asked: Vec<String>,
}

impl Done {
    fn read(&self, name: &str) -> String {
        String::from_utf8(fs::read(self.out.join(name)).unwrap()).unwrap()
    }

    fn pairs(&self) -> String {
        match &self.outcome {
            Ok(passed) => stated(passed),
            Err(message) => panic!("the passes refused: {message}"),
        }
    }

    /// The message the passes refused with, the scratch directory's name left out.
    fn refused(&self, scratch: &Scratch) -> String {
        let Err(message) = &self.outcome else {
            panic!("the passes refused nothing");
        };
        assert!(!self.out.join("shards").exists());
        message.replace(&format!("{}/", scratch.0.display()), "")
    }

    /// Every file of the remainder by name, and the counts as `counts`.
    fn files(&self) -> BTreeMap<String, Vec<u8>> {
        let mut files: BTreeMap<String, Vec<u8>> = fs::read_dir(&self.out)
            .unwrap()
            .map(|entry| entry.unwrap())
            .map(|entry| {
                (
                    entry.file_name().into_string().unwrap(),
                    fs::read(entry.path()).unwrap(),
                )
            })
            .collect();
        files.insert("counts".to_string(), self.pairs().into_bytes());
        files
    }
}

fn pair(
    scratch: &Scratch,
    tag: &str,
    reference: &Path,
    engines: &[&Path],
    options: &Options,
) -> Done {
    let out = scratch.0.join(format!("out-{tag}"));
    let planted = Planted::default();
    let registry = Registry::new(&planted);
    let engines: Vec<PathBuf> = engines.iter().map(|engine| engine.to_path_buf()).collect();
    let outcome = pass::run(reference, &engines, &out, options, &registry, &|_| {});
    let asked = planted.asked.lock().unwrap().clone();
    Done {
        outcome: outcome.map_err(|Refusal(message)| message),
        out,
        asked,
    }
}

fn refusal(scratch: &Scratch, tag: &str, reference: &Path, engine: &Path) -> String {
    pair(scratch, tag, reference, &[engine], &Options::new(false)).refused(scratch)
}

fn sha256(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

fn gzip(members: &[&[u8]]) -> Vec<u8> {
    let mut packed = Vec::new();
    for member in members {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(member).unwrap();
        packed.extend(encoder.finish().unwrap());
    }
    packed
}

fn reference_rows() -> Vec<String> {
    vec![
        row(
            "1:100",
            "T",
            "ENST1",
            "IMPACT=MODERATE;PLANT_score=0.5;PLANT_label=benign;Stray=1",
        ),
        row(
            "1:100",
            "T",
            "ENST2",
            "IMPACT=MODERATE;PLANT_score=0.5;Bare;PLANT_empty=",
        ),
        row("1:100", "T", "ENST3", "IMPACT=MODIFIER"),
        row(
            "1:200-201",
            "-",
            "ENST1",
            "PLANT_score=0.25;Stray=2;Other=x",
        ),
        row("10:5", "G", "-", "PLANT_label=-;GWAS_pmid=17"),
    ]
}

#[test]
fn populated_fields_of_the_runs_sources_are_counted_and_every_other_key_is_set_apart() {
    let scratch = Scratch::new("populated");
    let reference = scratch.output("reference.txt", &reference_rows());
    let done = pair(
        &scratch,
        "self",
        &reference,
        &[&reference],
        &Options::new(false),
    );
    let digest = sha256(&fs::read(&reference).unwrap());
    let input = format!(
        "{{\"kind\": \"default\", \"sha256\": \"{digest}\", \"lines\": 5, \"outside\": 0, \"header\": \"{}\", \"unclaimed\": {{\"{}\": 1, \"{}\": 2}}}}",
        hex(HEADER.trim_end().as_bytes()),
        hex(b"Other"),
        hex(b"Stray"),
    );
    let feature = format!(
        "{{\"values\": {{\"{}\": [1, 1, 1], \"{}\": [1, 1, 1], \"{}\": [3, 3, 3]}}, \"shared_keys\": {{}}, \"remainder\": [0, 0]}}",
        hex(b"GWAS_pmid"),
        hex(b"PLANT_label"),
        hex(b"PLANT_score"),
    );
    let allele = feature.replace("[3, 3, 3]", "[2, 2, 2]");
    assert_eq!(
        done.pairs(),
        format!(
            "{{\n  \"inputs\": [\n    {input},\n    {input}\n  ],\n  \"engines\": [\n    {{\"feature\": {feature}, \"allele\": {allele}, \"rows\": [0, 0]}}\n  ]\n}}\n"
        )
    );
    let mut asked = done.asked.clone();
    asked.sort();
    let mut names: Vec<String> = [
        "IMPACT",
        "PLANT_score",
        "PLANT_label",
        "Stray",
        "Other",
        "GWAS_pmid",
    ]
    .iter()
    .map(|name| format!("field {}", hex(name.as_bytes())))
    .collect();
    names.sort();
    assert_eq!(
        asked, names,
        "a name is asked once, and only where a value is populated"
    );
    let files = done.files();
    let names: Vec<&str> = files.keys().map(String::as_str).collect();
    assert_eq!(
        names,
        [
            "counts",
            "input1.allele.engine",
            "input1.allele.reference",
            "input1.feature.engine",
            "input1.feature.reference"
        ]
    );
    assert!(files
        .iter()
        .all(|(name, bytes)| name == "counts" || bytes.is_empty()));
}

#[test]
fn planted_missing_extra_and_changed_values_come_out_in_the_remainder() {
    let scratch = Scratch::new("planted");
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=kept"),
            row(
                "1:200",
                "T",
                "ENST1",
                "PLANT_score=0.25;PLANT_label=missing",
            ),
            row("1:300", "T", "ENST1", "PLANT_score=1;PLANT_label=before"),
            row("1:300", "T", "ENST2", "PLANT_score=1;PLANT_label=before"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:300", "T", "ENST2", "PLANT_score=1;PLANT_label=before"),
            row("1:300", "T", "ENST1", "PLANT_score=1;PLANT_label=after"),
            row("1:200", "T", "ENST1", "PLANT_score=0.25"),
            row(
                "1:100",
                "T",
                "ENST1",
                "PLANT_score=0.5;PLANT_label=kept;GWAS_pmid=extra",
            ),
        ],
    );
    let done = pair(
        &scratch,
        "planted",
        &reference,
        &[&engine],
        &Options::new(false),
    );
    assert_eq!(
        done.read("input1.feature.reference"),
        "1:200\tT\tENST1\tPLANT_label\tmissing\n1:300\tT\tENST1\tPLANT_label\tbefore\n"
    );
    assert_eq!(
        done.read("input1.feature.engine"),
        "1:100\tT\tENST1\tGWAS_pmid\textra\n1:300\tT\tENST1\tPLANT_label\tafter\n"
    );
    // Without the Feature the changed value's old text still pairs: ENST2 writes it.
    assert_eq!(
        done.read("input1.allele.reference"),
        "1:200\tT\tPLANT_label\tmissing\n"
    );
    assert_eq!(
        done.read("input1.allele.engine"),
        "1:100\tT\tGWAS_pmid\textra\n1:300\tT\tPLANT_label\tafter\n"
    );
    let pairs = done.pairs();
    let label = hex(b"PLANT_label");
    assert!(
        pairs.contains(&format!("\"{label}\": [4, 3, 2]")),
        "{pairs}"
    );
    assert!(
        pairs.contains(&format!("\"{label}\": [3, 3, 2]")),
        "{pairs}"
    );
    assert!(
        pairs.contains(&format!("\"{}\": [0, 1, 0]", hex(b"GWAS_pmid"))),
        "{pairs}"
    );
    assert!(
        pairs.contains("\"remainder\": [2, 2]") && pairs.contains("\"remainder\": [1, 2]"),
        "{pairs}"
    );
}

#[test]
fn rows_that_share_a_key_are_never_collapsed_and_the_allele_projection_keeps_a_value_once() {
    let scratch = Scratch::new("shared");
    let rows = reference_rows();
    let twice = scratch.output("twice.txt", &[rows.clone(), rows.clone()].concat());
    let once = scratch.output("once.txt", &rows);
    let done = pair(
        &scratch,
        "shared",
        &twice,
        &[&once, &twice],
        &Options::new(false),
    );
    let pairs = done.pairs();
    let score = hex(b"PLANT_score");
    let engines: Vec<&str> = pairs.split("\"feature\": ").skip(1).collect();
    // Against the output written once: every value of the reference twice, one copy paired.
    assert!(engines[0].starts_with(&format!(
        "{{\"values\": {{\"{}\": [2, 1, 1], \"{}\": [2, 1, 1], \"{score}\": [6, 3, 3]}}, \"shared_keys\": {{\"GWAS\": 1, \"Planted\": 3}}, \"remainder\": [5, 0]}}",
        hex(b"GWAS_pmid"),
        hex(b"PLANT_label"),
    )));
    assert!(engines[0].contains(&format!(
        "\"allele\": {{\"values\": {{\"{}\": [1, 1, 1], \"{}\": [1, 1, 1], \"{score}\": [2, 2, 2]}}, \"shared_keys\": {{}}, \"remainder\": [0, 0]}}",
        hex(b"GWAS_pmid"),
        hex(b"PLANT_label"),
    )));
    assert!(engines[1].contains(&format!("\"{score}\": [6, 6, 6]")));
    assert!(engines[1]
        .contains("\"shared_keys\": {\"GWAS\": 1, \"Planted\": 3}, \"remainder\": [0, 0]"));
    assert_eq!(
        done.read("input1.feature.reference"),
        "10:5\tG\t-\tGWAS_pmid\t17\n1:100\tT\tENST1\tPLANT_label\tbenign\n1:100\tT\tENST1\tPLANT_score\t0.5\n1:100\tT\tENST2\tPLANT_score\t0.5\n1:200-201\t-\tENST1\tPLANT_score\t0.25\n"
    );
    assert_eq!(done.read("input2.feature.reference"), "");
}

#[test]
fn percent_encoded_value_pairs_with_its_decoded_form_in_a_tuple_file() {
    let scratch = Scratch::new("percent");
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "PLANT_label=a%3Bb%3Dc;PLANT_score=%3d",
        )],
    );
    let engine = scratch.file(
        "engine.tsv",
        "1:100\tT\tENST1\tPLANT_label\ta;b=c\n1:100\tT\tENST1\tPLANT_score\t%3d\n1:100\tT\tENST1\tPLANT_none\t-\n1:100\tT\tENST1\tIMPACT\tLOW\n1:100\tT\tENST1\tStray\t1\n",
    );
    let done = pair(
        &scratch,
        "percent",
        &reference,
        &[&engine],
        &Options::new(false),
    );
    let pairs = done.pairs();
    assert!(
        pairs.contains(&format!("\"{}\": [1, 1, 1]", hex(b"PLANT_label"))),
        "{pairs}"
    );
    assert!(
        pairs.contains(&format!("\"{}\": [1, 1, 1]", hex(b"PLANT_score"))),
        "{pairs}"
    );
    // A tuple file's lines are all counted, and those of no source of the run set apart.
    assert!(
        pairs.contains(&format!(
            "\"kind\": \"tuples\", \"sha256\": \"{}\", \"lines\": 5, \"outside\": 2,",
            sha256(&fs::read(&engine).unwrap())
        )),
        "{pairs}"
    );
    assert!(
        pairs.contains(&format!("\"unclaimed\": {{\"{}\": 1}}", hex(b"Stray"))),
        "{pairs}"
    );
    assert!(!done
        .asked
        .contains(&format!("field {}", hex(b"PLANT_none"))));
}

#[test]
fn location_coordinates_compare_as_numbers_and_key_columns_lose_their_surrounding_whitespace() {
    let scratch = Scratch::new("location");
    let reference = scratch.output(
        "reference.txt",
        &[
            row("21:100-0200", " T", "ENST1\x0b", "PLANT_score=1"),
            row("21:\u{b2}", "T", "ENST1", "PLANT_score=1"),
            row(" 21:+7 ", "T", "ENST1", "PLANT_score=1"),
            row("21:8.5", "T\x1c", "ENST1", "PLANT_score=1"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("21:100-200", "T", "ENST1", "PLANT_score=1"),
            row("21:\u{b2}", "T", "ENST1", "PLANT_score=1"),
            row("21:7", "T", "ENST1", "PLANT_score=1"),
            row("21:8.5", "T", "ENST1", "PLANT_score=1"),
        ],
    );
    let done = pair(
        &scratch,
        "location",
        &reference,
        &[&engine],
        &Options::new(false),
    );
    assert!(done
        .pairs()
        .contains(&format!("\"{}\": [4, 4, 3]", hex(b"PLANT_score"))));
    // The file separator is no whitespace of Python's bytes, so that Allele stays its own.
    assert_eq!(
        done.read("input1.feature.reference"),
        "21:8.5\tT\x1c\tENST1\tPLANT_score\t1\n"
    );
    assert_eq!(
        done.read("input1.feature.engine"),
        "21:8.5\tT\tENST1\tPLANT_score\t1\n"
    );
    let odd = format!("location {}", hex("21:\u{b2}".as_bytes()));
    assert!(
        done.asked.contains(&odd),
        "a Location outside ASCII is the resolver's to settle"
    );
    assert!(done
        .asked
        .iter()
        .filter(|asked| asked.starts_with("location "))
        .all(|asked| *asked == odd));
}

#[test]
fn headerless_output_and_tuple_file_are_told_apart_by_the_first_data_line() {
    let scratch = Scratch::new("kinds");
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let bare = scratch.file("bare.txt", format!("\n\r\n{}", rows.concat()));
    let done = pair(&scratch, "bare", &reference, &[&bare], &Options::new(false));
    let pairs = done.pairs();
    assert_eq!(pairs.matches("\"kind\": \"default\"").count(), 2, "{pairs}");
    assert!(
        pairs.contains(&format!("\"{}\": [3, 3, 3]", hex(b"PLANT_score"))),
        "{pairs}"
    );
    // A tuple file has no header, so a line of it that opens with a `#` is a tuple.
    let listed = scratch.file(
        "listed.tsv",
        "\n1:100\tT\tENST1\tPLANT_score\t0.5\n\n#1:100\tT\tENST1\tPLANT_score\t0.5\n",
    );
    let done = pair(
        &scratch,
        "listed",
        &reference,
        &[&listed],
        &Options::new(false),
    );
    assert!(done.pairs().contains("\"kind\": \"tuples\""));
    assert_eq!(
        done.read("input1.feature.engine"),
        "#1:100\tT\tENST1\tPLANT_score\t0.5\n"
    );
    let headed = scratch.file(
        "headed.tsv",
        format!("{HEADER}1:100\tT\tENST1\tPLANT_score\t0.5\n"),
    );
    assert_eq!(
        refusal(&scratch, "headed", &reference, &headed),
        "headed.tsv: line 2 has 5 columns; a default-format row has 14, a tuple file line 5"
    );
}

#[test]
fn last_column_header_is_kept_and_comments_and_line_ends_are_no_rows() {
    let scratch = Scratch::new("header");
    let rows = reference_rows();
    let other = HEADER.replace("Extra", "Further");
    let reference = scratch.file(
        "reference.txt",
        format!(
            "## a comment\twith a tab\r\n{HEADER}{}\r\r\n{other}{}## last\n",
            rows[..2].concat(),
            rows[2..].concat().replace('\n', "\r\n")
        ),
    );
    let done = pair(
        &scratch,
        "header",
        &reference,
        &[&reference],
        &Options {
            block_bytes: 50,
            ..Options::new(false)
        },
    );
    let pairs = done.pairs();
    assert!(
        pairs.contains(&format!(
            "\"lines\": 5, \"outside\": 0, \"header\": \"{}\"",
            hex(other.trim_end().as_bytes())
        )),
        "{pairs}"
    );
    assert!(
        pairs.contains(&format!("\"{}\": [3, 3, 3]", hex(b"PLANT_score"))),
        "{pairs}"
    );
}

#[test]
fn malformed_inputs_are_refused_in_the_scorers_words() {
    let scratch = Scratch::new("malformed");
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let short = scratch.file(
        "short.tsv",
        "1:100\tT\tENST1\tPLANT_score\t0.5\n1:100\tT\tENST1\tPLANT_score\n",
    );
    assert_eq!(
        refusal(&scratch, "short", &reference, &short),
        "short.tsv: line 2 has 4 columns; a default-format row has 14, a tuple file line 5"
    );
    let listed = scratch.file("listed.tsv", "1:100\tT\tENST1\tPLANT_score\t1\n");
    assert_eq!(
        refusal(&scratch, "listed", &listed, &reference),
        "listed.tsv: line 1 has 5 columns; a default-format row has 14"
    );
    let empty = scratch.output("empty.txt", &[]);
    assert_eq!(
        refusal(&scratch, "empty", &empty, &reference),
        "empty.txt: the reference has no data row"
    );
    let odd = scratch.file(
        "odd.txt",
        format!("{HEADER}{}#Uploaded_variation\tLocation\n", rows[0]),
    );
    assert_eq!(
        refusal(&scratch, "odd", &reference, &odd),
        "odd.txt: line 3, the column header, names 2 columns, not 14"
    );
    // The first fault of the file is the one named, whichever thread meets one first.
    let mut many: Vec<String> = (0..400)
        .map(|at| row(&format!("1:{at}"), "T", "ENST1", "PLANT_score=1"))
        .collect();
    many[150] = "1:150\tT\n".to_string();
    many[300] = "1:300\tT\tENST1\n".to_string();
    let faulty = scratch.output("faulty.txt", &many);
    let done = pair(
        &scratch,
        "faulty",
        &reference,
        &[&faulty],
        &Options {
            block_bytes: 200,
            threads: 8,
            ..Options::new(false)
        },
    );
    assert_eq!(
        done.refused(&scratch),
        "faulty.txt: line 152 has 2 columns; a default-format row has 14, a tuple file line 5"
    );
    // One distinct key too many comes before a later fault, though no thread meets every key.
    let mut spread: Vec<String> = (0..400)
        .map(|at| {
            let extra = match at % 40 {
                0 => format!("Key{at}=1;PLANT_score=1"),
                _ => "PLANT_score=1".to_string(),
            };
            row(&format!("1:{at}"), "T", "ENST1", &extra)
        })
        .collect();
    spread.push("1:400\tT\n".to_string());
    let spread = scratch.output("spread.txt", &spread);
    let done = pair(
        &scratch,
        "spread",
        &reference,
        &[&spread],
        &Options {
            block_bytes: 200,
            threads: 8,
            distinct_keys: 5,
            ..Options::new(false)
        },
    );
    assert_eq!(
        done.refused(&scratch),
        "spread.txt: more than 5 distinct keys that no source claims"
    );
    let packed = gzip(&[&fs::read(&reference).unwrap()]);
    let cut = scratch.file("cut.txt.gz", &packed[..packed.len() - 12]);
    assert!(refusal(&scratch, "cut", &cut, &reference)
        .starts_with("cut.txt.gz: cannot read the gzip stream: "));
    let wide = scratch.output(
        "wide.txt",
        &[row("1:100", "T", "ENST1", "One=1;Two=2;Three=3;One=4")],
    );
    let done = pair(
        &scratch,
        "wide",
        &reference,
        &[&wide],
        &Options {
            distinct_keys: 2,
            ..Options::new(false)
        },
    );
    assert_eq!(
        done.refused(&scratch),
        "wide.txt: more than 2 distinct keys that no source claims"
    );
    let done = pair(
        &scratch,
        "narrow",
        &reference,
        &[&wide],
        &Options {
            distinct_keys: 3,
            ..Options::new(false)
        },
    );
    assert!(done.pairs().contains(&format!("\"{}\": 2", hex(b"One"))));
    let refused = scratch.output("refused.txt", &[row("1:100", "T", "ENST1", "Refused=1")]);
    assert_eq!(
        refusal(&scratch, "refused", &reference, &refused),
        format!("the registry refuses the field {}", hex(b"Refused"))
    );
}

#[test]
fn gzip_input_is_read_whole_and_digested_as_stored() {
    let scratch = Scratch::new("gzip");
    let rows = reference_rows();
    let plain = scratch.output("reference.txt", &rows);
    let bytes = fs::read(&plain).unwrap();
    let mut stored = gzip(&[&bytes[..200], &bytes[200..]]);
    stored.extend([0u8; 37]);
    let packed = scratch.file("reference.txt.gz", &stored);
    let done = pair(
        &scratch,
        "gzip",
        &packed,
        &[&plain, &packed],
        &Options {
            block_bytes: 64,
            ..Options::new(false)
        },
    );
    let pairs = done.pairs();
    assert_eq!(
        pairs
            .matches(&format!("\"sha256\": \"{}\"", sha256(&stored)))
            .count(),
        2,
        "{pairs}"
    );
    assert_eq!(
        pairs
            .matches(&format!("\"sha256\": \"{}\"", sha256(&bytes)))
            .count(),
        1,
        "{pairs}"
    );
    assert_eq!(
        pairs
            .matches(&format!("\"{}\": [3, 3, 3]", hex(b"PLANT_score")))
            .count(),
        2,
        "{pairs}"
    );
    assert_eq!(pairs.matches("\"lines\": 5").count(), 3, "{pairs}");
    assert!(
        !done.out.join("shards").exists(),
        "the unpacked copy leaves with the scratch"
    );
}

#[test]
fn rows_of_the_remainder_are_written_for_a_class_to_read() {
    let scratch = Scratch::new("rows");
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=0.5"),
            row("1:200", "T", "ENST1", "PLANT_score=0.25;PLANT_label=before"),
            row("1:200", "T", "ENST2", "IMPACT=LOW"),
            row("1:0200", "C", "ENST1", "-"),
            row("1:300", "T", "ENST1", "PLANT_score=1"),
        ],
    );
    let changed = row("1:200", "T", "ENST1", "PLANT_score=0.25;PLANT_label=after");
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=0.5"),
            changed.clone(),
            row("1:200", "T", "ENST2", "IMPACT=LOW"),
            changed.clone(),
            row("1:300", "T", "ENST1", "PLANT_score=1"),
        ],
    );
    let listed = scratch.file(
        "engine.tsv",
        "1:200\tT\tENST1\tPLANT_score\t0.25\n1:200\tT\tENST1\tPLANT_label\tafter\n1:200\tT\tENST1\tPLANT_score\t0.25\n1:300\tT\tENST1\tPLANT_score\t1\n1:100\tT\tENST1\tPLANT_score\t0.5\n",
    );
    let done = pair(
        &scratch,
        "rows",
        &reference,
        &[&engine, &listed],
        &Options::new(true),
    );
    let kept = |ordinal: usize, key: &str, row: &str| format!("{key}\t{ordinal:012}\t{}", row);
    let reference_rows = [
        kept(3, "1:200\tC\tENST1", &row("1:0200", "C", "ENST1", "-")),
        kept(
            1,
            "1:200\tT\tENST1",
            &row("1:200", "T", "ENST1", "PLANT_score=0.25;PLANT_label=before"),
        ),
        kept(
            2,
            "1:200\tT\tENST2",
            &row("1:200", "T", "ENST2", "IMPACT=LOW"),
        ),
    ]
    .concat();
    // Each side's rows of the whole Location.
    assert_eq!(done.read("input1.rows.reference"), reference_rows);
    assert_eq!(
        done.read("input1.rows.engine"),
        [
            kept(1, "1:200\tT\tENST1", &changed),
            kept(3, "1:200\tT\tENST1", &changed),
            kept(
                2,
                "1:200\tT\tENST2",
                &row("1:200", "T", "ENST2", "IMPACT=LOW")
            ),
        ]
        .concat()
    );
    assert_eq!(done.read("input2.rows.reference"), reference_rows);
    // A tuple file states no rows: its tuples under the key stand for them, the paired ones too.
    assert_eq!(
        done.read("input2.rows.engine"),
        "1:200\tT\tENST1\tPLANT_label\tafter\n1:200\tT\tENST1\tPLANT_score\t0.25\n1:200\tT\tENST1\tPLANT_score\t0.25\n"
    );
    assert_eq!(done.pairs().matches("\"rows\": [3, 3]").count(), 2);
    let without = pair(
        &scratch,
        "no-rows",
        &reference,
        &[&engine],
        &Options::new(false),
    );
    assert!(without.pairs().contains("\"rows\": [0, 0]"));
    assert!(
        !without.out.join("input1.rows.reference").exists()
            && !without.out.join("input1.rows.engine").exists()
    );
}

#[test]
fn engine_rows_of_a_remainder_location_are_written_whole_for_a_class_that_reads_a_records_other_alleles(
) {
    let scratch = Scratch::new("record-rows");
    let own = row("1:500", "T", "ENST1", "PLANT_score=0.75");
    let rows = |second: &str| {
        [
            row("1:500", "C", "ENST1", "PLANT_score=0.5"),
            second.to_string(),
            row("1:500", "T", "ENST2", "IMPACT=LOW"),
            row("1:700", "C", "ENST1", "PLANT_score=0.5"),
        ]
    };
    let reference = scratch.output(
        "reference.txt",
        &rows(&row("1:500", "T", "ENST1", "PLANT_score=0.5")),
    );
    let engine = scratch.output("engine.txt", &rows(&own));
    let done = pair(
        &scratch,
        "record-rows",
        &reference,
        &[&engine],
        &Options::new(true),
    );
    assert_eq!(
        done.read("input1.feature.engine"),
        "1:500\tT\tENST1\tPLANT_score\t0.75\n"
    );
    let kept = |ordinal: usize, key: &str, row: &str| format!("{key}\t{ordinal:012}\t{}", row);
    // The other allele's row, whose value paired, and the row with no claimed field stand
    // beside the row of the remainder; a Location with no value left over has no row.
    assert_eq!(
        done.read("input1.rows.engine"),
        [
            kept(
                0,
                "1:500\tC\tENST1",
                &row("1:500", "C", "ENST1", "PLANT_score=0.5")
            ),
            kept(1, "1:500\tT\tENST1", &own),
            kept(
                2,
                "1:500\tT\tENST2",
                &row("1:500", "T", "ENST2", "IMPACT=LOW")
            ),
        ]
        .concat()
    );
    assert!(done.pairs().contains("\"rows\": [3, 3]"));
}

/// A small generator, so that the outputs are the same on every run.
struct Draws(u64);

impl Draws {
    fn below(&mut self, bound: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) as usize) % bound
    }

    fn of<'a>(&mut self, choices: &[&'a str]) -> &'a str {
        choices[self.below(choices.len())]
    }
}

fn drawn_rows(draws: &mut Draws, count: usize) -> Vec<String> {
    let locations = [
        "1:5",
        "1:50",
        "1:5-06",
        "10:5",
        "2:5",
        "MT:5",
        "1:5\x01",
        "X:1",
        "21:\u{b2}",
        " 3:07",
    ];
    let alleles = ["A", "-", "AT", "T ", "a"];
    let features = ["ENST1", "ENST10", "ENST2", "-", "ENST1.2"];
    let values = ["a", "b", "0.1", "0.10", "7", "a%3Bb", "a\x01", "", "-"];
    let names = [
        "PLANT_score",
        "PLANT_label",
        "PLANT_",
        "PLANT_score\x01",
        "GWAS_pmid",
        "IMPACT",
        "Stray",
    ];
    let mut rows = Vec::new();
    for _ in 0..count {
        let mut pairs = Vec::new();
        for name in names {
            if draws.below(10) < 6 {
                pairs.push(format!("{name}={}", draws.of(&values)));
            }
        }
        let extra = if pairs.is_empty() {
            "-".to_string()
        } else {
            pairs.join(";")
        };
        rows.push(row(
            draws.of(&locations),
            draws.of(&alleles),
            draws.of(&features),
            &extra,
        ));
    }
    rows
}

#[test]
fn two_runs_write_the_same_bytes_whatever_the_threads_the_shards_and_the_blocks() {
    let scratch = Scratch::new("threads");
    let mut draws = Draws(116);
    let reference = scratch.output("reference.txt", &drawn_rows(&mut draws, 700));
    let engine = scratch.output("engine.txt", &drawn_rows(&mut draws, 700));
    let tuples: String = drawn_rows(&mut draws, 300)
        .iter()
        .map(|row| {
            let columns: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
            format!(
                "{}\t{}\t{}\tPLANT_score\t{}\n",
                columns[1],
                columns[2],
                columns[4],
                columns[13].len()
            )
        })
        .collect();
    let listed = scratch.file("engine.tsv", tuples);
    let engines = [engine.as_path(), listed.as_path()];
    let first = pair(
        &scratch,
        "one",
        &reference,
        &engines,
        &Options {
            threads: 1,
            shards: Some(1),
            ..Options::new(true)
        },
    )
    .files();
    assert!(
        first["input1.feature.reference"].len() > 1000 && first["input1.rows.engine"].len() > 1000
    );
    assert!(first["input2.rows.engine"].len() > 100 && first["input2.allele.engine"].len() > 100);
    for (tag, threads, shards, block_bytes) in [
        ("again", 1, 1, 8388608),
        ("many", 7, 13, 64),
        ("few", 3, 4, 1000),
        ("wide", 16, 64, 1),
    ] {
        let options = Options {
            threads,
            shards: Some(shards),
            block_bytes,
            ..Options::new(true)
        };
        let other = pair(&scratch, tag, &reference, &engines, &options).files();
        let names: Vec<&String> = first.keys().collect();
        assert_eq!(other.keys().collect::<Vec<_>>(), names, "{tag}");
        for name in names {
            assert!(first[name] == other[name], "{tag}: {name} differs");
        }
    }
    // Every file is in the order `sort` gives its lines under C collation.
    for (name, bytes) in &first {
        if name != "counts" {
            let lines: Vec<&[u8]> = bytes
                .split(|&byte| byte == b'\n')
                .filter(|line| !line.is_empty())
                .collect();
            assert!(lines.windows(2).all(|pair| pair[0] <= pair[1]), "{name}");
        }
    }
}
