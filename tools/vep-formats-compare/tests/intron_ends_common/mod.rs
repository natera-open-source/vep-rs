// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the intron classes share: a cache of four transcripts, one insertion's
//! row in the four formats and its verdicts.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::Path;
use std::path::PathBuf;

use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::consts::HGVSP_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::intron_ends::hgvsp_intron_end_insertion;
use vep_formats_compare::rules::intron_ends::hgvsp_intron_start_insertion;
use vep_formats_compare::rules::intron_ends::IntronEnds;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::blank;

use crate::compare_common::counts;
use crate::compare_common::default_row;
use crate::compare_common::default_text;
use crate::compare_common::documented_with;
use crate::compare_common::json;
use crate::compare_common::verdict;
use crate::compare_common::DEFAULT_HEADER;

pub const INTRON_END: &str = "hgvsp_intron_end_insertion";
pub const INTRON_START: &str = "hgvsp_intron_start_insertion";
pub const INTRON_END_PROTEIN: &str = "ENSP00000000001.1:p.Gly17ArgfsTer5";
pub const FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// A cache transcript `ENST0000000000<number>` of `exons` in the transcript's order, coding
/// from the first to the second cDNA position of `coding`.
pub fn intron_end_transcript(
    number: i64,
    strand: i64,
    coding: (i64, i64),
    exons: &[(i64, i64)],
) -> Value {
    let exons: Vec<String> = exons
        .iter()
        .map(|(start, end)| format!(r#"{{"start": {start}, "end": {end}}}"#))
        .collect();
    json(&format!(
        r#"{{"stable_id": "ENST0000000000{number}", "strand": {strand}, "cdna_coding_start": {}, "cdna_coding_end": {}, "exons": [{}]}}"#,
        coding.0,
        coding.1,
        exons.join(", ")
    ))
}

/// Four transcripts of chromosome 1. The first (forward) and the second (reverse) have three
/// exons of 100 bases and code from their 51st base to their 250th: the second exon opens on
/// coding base 51 behind an intron of 100 bases that ends at 11199 (forward), or of 200 bases
/// whose last base is 25300 (reverse). The third codes from its 151st base to its 250th, so
/// its second exon opens in the 5' UTR, its third on coding base 51 and its fourth in the 3'
/// UTR. The fourth codes over its whole length behind introns of 13, 2 and 14 bases, its
/// fourth exon opening on coding base 286.
pub fn intron_end_transcripts() -> [Value; 4] {
    [
        intron_end_transcript(
            1,
            1,
            (51, 250),
            &[(11000, 11099), (11200, 11299), (11500, 11599)],
        ),
        intron_end_transcript(
            2,
            -1,
            (51, 250),
            &[(25500, 25599), (25200, 25299), (25000, 25099)],
        ),
        intron_end_transcript(
            3,
            1,
            (151, 250),
            &[
                (31000, 31099),
                (31200, 31299),
                (31400, 31499),
                (31600, 31699),
            ],
        ),
        intron_end_transcript(
            4,
            1,
            (1, 371),
            &[
                (41000, 41099),
                (41113, 41199),
                (41202, 41299),
                (41314, 41399),
            ],
        ),
    ]
}

pub fn dumped(value: &Value) -> String {
    pyjson::dumps(value, &pyjson::DumpOptions::default())
}

/// A cache of one shard, `transcripts/1/1-1000000.json`, holding `shard`.
pub fn one_shard_cache(root: &Path, shard: &Value) -> PathBuf {
    let directory = root.join("transcripts").join("1");
    std::fs::create_dir_all(&directory).expect("the cache's directory is made");
    std::fs::write(directory.join("1-1000000.json"), dumped(shard)).expect("the shard is written");
    root.to_path_buf()
}

/// A candidate JSON cache of the four transcripts on chromosome 1.
pub fn intron_end_cache(root: &Path) -> PathBuf {
    one_shard_cache(root, &Value::List(intron_end_transcripts().to_vec()))
}

/// One insertion's row on a transcript: `allele` inserted after `position` of chromosome
/// 1, the input line's REF `reference_base` and its ALT that base followed by the allele.
/// Both rows carry the HGVSc `notation` of the transcript and no HGVS_OFFSET, the reference's
/// the HGVSp `INTRON_END_PROTEIN` and the candidate's none.
#[derive(Clone)]
pub struct Row {
    feature: &'static str,
    position: i64,
    reference_base: &'static str,
    allele: &'static str,
    notation: &'static str,
    strand: i64,
    location: Option<&'static str>,
    record_reference: Option<&'static str>,
    record_alternates: Option<&'static str>,
    reference: Vec<(&'static str, &'static str)>,
    candidate: Vec<(&'static str, &'static str)>,
    input_line: Option<&'static str>,
}

pub fn row(
    feature: &'static str,
    position: i64,
    reference_base: &'static str,
    allele: &'static str,
    notation: &'static str,
) -> Row {
    Row {
        feature,
        position,
        reference_base,
        allele,
        notation,
        strand: 1,
        location: None,
        record_reference: None,
        record_alternates: None,
        reference: Vec::new(),
        candidate: Vec::new(),
        input_line: None,
    }
}

/// One format's reference text, candidate text and the input lines the format reads.
pub type Texts = Vec<(&'static str, String, String, Option<Vec<String>>)>;

impl Row {
    pub fn strand(mut self, strand: i64) -> Row {
        self.strand = strand;
        self
    }

    /// The text formats' Location.
    pub fn location(mut self, location: &'static str) -> Row {
        self.location = Some(location);
        self
    }

    /// The input line's REF.
    pub fn record_reference(mut self, reference: &'static str) -> Row {
        self.record_reference = Some(reference);
        self
    }

    /// The input line's ALT, for a record that holds the insertion with a shared base beside
    /// another allele.
    pub fn record_alternates(mut self, alternates: &'static str) -> Row {
        self.record_alternates = Some(alternates);
        self
    }

    /// A field of the reference's row.
    pub fn reference(mut self, name: &'static str, value: &'static str) -> Row {
        self.reference.push((name, value));
        self
    }

    /// A field of the candidate's row.
    pub fn candidate(mut self, name: &'static str, value: &'static str) -> Row {
        self.candidate.push((name, value));
        self
    }

    /// A field of both rows.
    pub fn both(self, name: &'static str, value: &'static str) -> Row {
        self.reference(name, value).candidate(name, value)
    }

    /// The input line where it is no VCF line: the VCF format writes the VCF line all the
    /// same, and the JSON format's `input` is the line as given.
    pub fn input_line(mut self, line: &'static str) -> Row {
        self.input_line = Some(line);
        self
    }

    /// The row in the four formats, as the reference and as the candidate write it, with the
    /// input lines the VCF and JSON formats read.
    pub fn texts(&self) -> Texts {
        let hgvsc = format!("{}.1:{}", self.feature, self.notation);
        let side = |protein: &str, replaced: &[(&'static str, &'static str)]| {
            let mut fields = [
                ("HGVSc", hgvsc.clone()),
                ("HGVSp", protein.to_string()),
                ("HGVS_OFFSET", String::new()),
            ];
            for (name, value) in replaced {
                let field = fields
                    .iter_mut()
                    .find(|(held, _)| held == name)
                    .expect("one of the row's three fields");
                field.1 = value.to_string();
            }
            fields
        };
        let sides = [
            side(INTRON_END_PROTEIN, &self.reference),
            side("", &self.candidate),
        ];
        let location = self.location.map_or_else(
            || format!("1:{}-{}", self.position, self.position + 1),
            str::to_string,
        );
        let record_reference = self.record_reference.unwrap_or(self.reference_base);
        let alternates = self.record_alternates.map_or_else(
            || format!("{record_reference}{}", self.allele),
            str::to_string,
        );
        let line = format!(
            "1\t{}\trec1\t{record_reference}\t{alternates}\t.\t.\t.",
            self.position
        );
        let allele_string = match self.record_alternates {
            Some(alternates) => format!("{record_reference}/{}", alternates.replace(',', "/")),
            None => format!("-/{}", self.allele),
        };
        let terms = "frameshift_variant,splice_region_variant";
        let tab_header =
            DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSc\tHGVSp\tHGVS_OFFSET");
        let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|STRAND|HGVSc|HGVSp|HGVS_OFFSET\">\n";
        let strand = self.strand.to_string();
        let values = |fields: &[(&str, String); 3]| -> Vec<String> {
            fields.iter().map(|(_, value)| value.clone()).collect()
        };

        let default = |fields: &[(&str, String); 3]| {
            let mut extra = format!("IMPACT=HIGH;STRAND={strand}");
            for (name, value) in fields.iter().filter(|(_, value)| !value.is_empty()) {
                extra.push_str(&format!(";{name}={value}"));
            }
            default_text(&[default_row(
                "rec1",
                &location,
                self.allele,
                self.feature,
                terms,
                &extra,
            )])
        };
        let tab = |fields: &[(&str, String); 3]| {
            let mut columns = vec!["HIGH".to_string(), strand.clone()];
            columns.extend(values(fields).into_iter().map(|value| {
                if value.is_empty() {
                    "-".to_string()
                } else {
                    value
                }
            }));
            format!(
                "{tab_header}\n{}\n",
                default_row(
                    "rec1",
                    &location,
                    self.allele,
                    self.feature,
                    terms,
                    &columns.join("\t")
                )
            )
        };
        let vcf = |fields: &[(&str, String); 3]| {
            let mut csq = vec![
                self.allele.to_string(),
                terms.replace(',', "&"),
                self.feature.to_string(),
                "Transcript".to_string(),
                strand.clone(),
            ];
            csq.extend(values(fields));
            format!(
                "{csq_header}{}CSQ={}\n",
                &line[..line.len() - 1],
                csq.join("|")
            )
        };
        let given = self.input_line.map_or_else(|| line.clone(), str::to_string);
        let text = |value: &str| pyjson::dumps_str(value, true);
        let object = |fields: &[(&str, String); 3]| {
            let mut entry = format!(
                r#"{{"variant_allele": {}, "transcript_id": {}, "consequence_terms": ["frameshift_variant", "splice_region_variant"], "strand": {strand}"#,
                text(self.allele),
                text(self.feature)
            );
            for (name, value) in fields.iter().filter(|(_, value)| !value.is_empty()) {
                entry.push_str(&format!(
                    ", {}: {}",
                    text(&name.to_lowercase()),
                    text(value)
                ));
            }
            entry.push('}');
            format!(
                r#"{{"input": {}, "id": "rec1", "allele_string": {}, "most_severe_consequence": "frameshift_variant", "transcript_consequences": [{entry}]}}"#,
                text(&given),
                text(&allele_string)
            ) + "\n"
        };
        vec![
            ("default", default(&sides[0]), default(&sides[1]), None),
            ("tab", tab(&sides[0]), tab(&sides[1]), None),
            (
                "vcf",
                vcf(&sides[0]),
                vcf(&sides[1]),
                Some(vec![given.clone()]),
            ),
            (
                "json",
                object(&sides[0]),
                object(&sides[1]),
                Some(vec![given]),
            ),
        ]
    }

    /// The reference's and the candidate's text of one format, and its input lines.
    pub fn text(&self, fmt: &str) -> (String, String, Option<Vec<String>>) {
        let (_, reference, candidate, inputs) = self
            .texts()
            .into_iter()
            .find(|(held, ..)| *held == fmt)
            .expect("one of the four formats");
        (reference, candidate, inputs)
    }
}

/// A format's verdict and what was applied.
pub type Verdict = (String, BTreeMap<String, i64>);
pub type FormatVerdicts = BTreeMap<&'static str, Verdict>;

/// How a row's verdicts are read: through the comparison or through the rule alone.
pub type Verdicts = fn(Option<&Path>, &Texts, &[&str]) -> FormatVerdicts;

pub fn green(klass: &str) -> Verdict {
    (
        "green".to_string(),
        BTreeMap::from([(klass.to_string(), 1)]),
    )
}

pub fn red() -> Verdict {
    ("red".to_string(), BTreeMap::new())
}

/// The row set aside under `klass` in every format, counted once.
pub fn set_aside(klass: &str) -> FormatVerdicts {
    FORMATS.iter().map(|fmt| (*fmt, green(klass))).collect()
}

/// The row a difference in every format.
pub fn difference() -> FormatVerdicts {
    FORMATS.iter().map(|fmt| (*fmt, red())).collect()
}

/// The verdicts of the two text formats, which a row's Location places, and of the VCF and
/// JSON formats, which the input line places.
pub fn placed(text: Verdict, entry: Verdict) -> FormatVerdicts {
    BTreeMap::from([
        ("default", text.clone()),
        ("tab", text),
        ("vcf", entry.clone()),
        ("json", entry),
    ])
}

/// The verdicts of the two text formats alone.
pub fn text_rows(verdicts: &FormatVerdicts) -> [&Verdict; 2] {
    [&verdicts["default"], &verdicts["tab"]]
}

pub fn parsed(fmt: &str, text: &str, inputs: Option<&[String]>) -> Vec<Entry> {
    let lines = inputs.unwrap_or_default();
    match fmt {
        "default" => parse_default(text),
        "tab" => parse_tab(text),
        "vcf" => parse_vcf(text, lines),
        _ => parse_json(text, lines),
    }
    .expect("the text parses")
}

fn intron_ends(cache: Option<&Path>, inputs: Option<&[String]>) -> Option<IntronEnds> {
    cache.map(|cache| IntronEnds::new(cache, inputs).expect("the cache has its transcripts"))
}

/// Each format's verdict and applied counts for one row under `classes`, the class reading
/// the intron ends of `cache` (none without a cache).
pub fn compared(cache: Option<&Path>, texts: &Texts, classes: &[&str]) -> FormatVerdicts {
    texts
        .iter()
        .map(|(fmt, reference, candidate, inputs)| {
            let inputs = inputs.as_deref();
            let mut documented = documented_with(classes);
            let side = Side {
                intron_ends: intron_ends(cache, inputs),
                ..Side::default()
            };
            let options = CompareOptions {
                input_lines: inputs,
                ..CompareOptions::new(&side)
            };
            let report = compare_entries(
                &parsed(fmt, reference, inputs),
                &parsed(fmt, candidate, inputs),
                &mut documented,
                &options,
            )
            .expect("the entries compare");
            (
                *fmt,
                (verdict(&report).to_string(), counts(&documented.applied)),
            )
        })
        .collect()
}

/// Each format's verdict for one row as the two rules alone give it: the row is set aside
/// under the first of `classes` whose rule accepts the row's differing HGVSp.
pub fn by_the_rule(cache: Option<&Path>, texts: &Texts, classes: &[&str]) -> FormatVerdicts {
    texts
        .iter()
        .map(|(fmt, reference, candidate, inputs)| {
            let inputs = inputs.as_deref();
            let ends = intron_ends(cache, inputs);
            let of_the_transcript = |text: &str| -> Entry {
                parsed(fmt, text, inputs)
                    .into_iter()
                    .find(|entry| entry.feature.starts_with("ENST"))
                    .expect("the row's entry")
            };
            let (e, a) = (of_the_transcript(reference), of_the_transcript(candidate));
            let shown = |entry: &Entry| {
                let value = field_of(Some(entry), HGVSP_FIELDS);
                if blank(&value) {
                    "-".to_string()
                } else {
                    value
                }
            };
            let (ev, av) = (shown(&e), shown(&a));
            assert_ne!(ev, av, "the two rows differ in their HGVSp");
            let taken = classes.iter().find(|klass| {
                match **klass {
                    INTRON_END => hgvsp_intron_end_insertion(&e, Some(&a), &ev, &av, ends.as_ref()),
                    INTRON_START => {
                        hgvsp_intron_start_insertion(&e, Some(&a), &ev, &av, ends.as_ref())
                    }
                    _ => Ok(false),
                }
                .expect("the cache reads")
            });
            (*fmt, taken.map_or_else(red, |klass| green(klass)))
        })
        .collect()
}

/// What the help says of `option`: the lines under its name, as one line.
pub fn option_help(help: &str, option: &str) -> String {
    let mut lines = help
        .lines()
        .skip_while(|line| line.trim_start().split(' ').next() != Some(option))
        .skip(1);
    let mut under = Vec::new();
    for line in lines.by_ref() {
        if !line.starts_with("    ") {
            break;
        }
        under.push(line.trim());
    }
    under.join(" ")
}

/// An input VCF of `records`, tab-separated data lines.
pub fn input_vcf(path: &Path, records: &[String]) -> PathBuf {
    let mut text =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string();
    for record in records {
        text.push_str(record);
        text.push('\n');
    }
    std::fs::write(path, text).expect("the input is written");
    path.to_path_buf()
}
