// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of `classify` share: a scratch corpus with the two outputs of a run, the
//! rows of a default-format output, and the JSON texts of a manifest's records and of the
//! entries `classify` writes.

#![allow(dead_code)]

use crate::common::{run, scratch, Run};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs;
use std::io::Write;
use std::path::Path;
use vep_golden_corpus::json::{dumps, dumps_object, Style, Value};
use vep_golden_corpus::manifest::manifest_text;
use vep_golden_corpus::TempDir;

/// A scratch directory with a corpus `c` beside the outputs `vep.txt` and `rs.txt`.
pub struct Case {
    pub scratch: TempDir,
    pub corpus: String,
}

impl Case {
    /// A case whose corpus directory is `c`.
    pub fn new(tag: &str) -> Case {
        Case::named(tag, "c")
    }

    /// A case whose corpus directory is `name`.
    pub fn named(tag: &str, name: &str) -> Case {
        let scratch = scratch(tag);
        let corpus = scratch.path.join(name);
        fs::create_dir(&corpus).unwrap();
        Case {
            corpus: corpus.to_str().unwrap().to_string(),
            scratch,
        }
    }

    /// The path of `name` beside the corpus.
    pub fn path(&self, name: &str) -> String {
        self.scratch.path.join(name).to_str().unwrap().to_string()
    }

    /// Writes `name` beside the corpus.
    pub fn write(&self, name: &str, text: &str) {
        fs::write(self.path(name), text).unwrap();
    }

    /// Writes `name` inside the corpus, with the directories above it.
    pub fn write_in_corpus(&self, name: &str, bytes: &[u8]) {
        let path = Path::new(&self.corpus).join(name);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, bytes).unwrap();
    }

    /// Writes the two outputs of the run.
    pub fn outputs(&self, vep: &str, vep_rs: &str) {
        self.write("vep.txt", vep);
        self.write("rs.txt", vep_rs);
    }

    /// Writes the manifest as the text given.
    pub fn manifest(&self, text: &str) {
        fs::write(Path::new(&self.corpus).join("manifest.json"), text).unwrap();
    }

    /// The manifest's bytes.
    pub fn text(&self) -> String {
        fs::read_to_string(Path::new(&self.corpus).join("manifest.json")).unwrap()
    }

    /// The manifest as `json.dumps` writes it on one line, with every key in its place; the
    /// file itself is held to the layout a manifest is written in.
    pub fn after(&self) -> String {
        let manifest = crate::common::manifest(&self.corpus);
        assert_eq!(
            self.text(),
            manifest_text(&manifest),
            "a manifest is written with one space an indent level"
        );
        dumps_object(&manifest, &Style::PLAIN)
    }

    /// Runs `classify` on the corpus and the two outputs with the further arguments.
    pub fn classify(&self, more: &[&str]) -> Run {
        self.classify_against("rs.txt", more)
    }

    /// Runs `classify` with the file `vep_rs` beside the corpus as vep-rs's output.
    pub fn classify_against(&self, vep_rs: &str, more: &[&str]) -> Run {
        let (vep, vep_rs) = (self.path("vep.txt"), self.path(vep_rs));
        let mut args = vec![
            "classify",
            "--corpus",
            &self.corpus,
            "--vep-default",
            &vep,
            "--vep-rs-output",
            &vep_rs,
        ];
        args.extend(more);
        run(&args)
    }
}

/// The bytes of a gzip file that holds `text`.
pub fn gzip(text: &str) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(text.as_bytes()).unwrap();
    encoder.finish().unwrap()
}

/// A row of a default-format output on a transcript of the gene `ENSG1`.
pub fn row(name: &str, location: &str, allele: &str, feature: &str, terms: &str) -> String {
    row_with(
        name,
        location,
        allele,
        feature,
        terms,
        "IMPACT=MODIFIER;STRAND=1",
    )
}

/// `row` with the Extra column given.
pub fn row_with(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    terms: &str,
    extra: &str,
) -> String {
    let columns = [
        name,
        location,
        allele,
        "ENSG1",
        feature,
        "Transcript",
        terms,
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
        extra,
    ];
    columns.join("\t") + "\n"
}

/// A row of a default-format output on a feature of the type given, with no gene.
pub fn feature_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: (&str, &str),
    terms: &str,
    extra: &str,
) -> String {
    let columns = [
        name, location, allele, "-", feature.0, feature.1, terms, "-", "-", "-", "-", "-", "-",
        extra,
    ];
    columns.join("\t") + "\n"
}

/// A transcript as the feature of a key.
pub fn transcript(name: &str) -> (&str, &str) {
    (name, "Transcript")
}

/// An entry with a further member after its last.
pub fn with_member(entry: &str, member: &str) -> String {
    format!("{}, {member}}}", entry.strip_suffix('}').unwrap())
}

/// A JSON string as `json.dumps` writes it.
pub fn js(text: &str) -> String {
    dumps(&Value::from(text), &Style::PLAIN)
}

/// A JSON object of the members given.
pub fn obj(members: &[&str]) -> String {
    format!("{{{}}}", members.join(", "))
}

/// A JSON list of the items given.
pub fn list<S: AsRef<str>>(items: &[S]) -> String {
    let items: Vec<&str> = items.iter().map(AsRef::as_ref).collect();
    format!("[{}]", items.join(", "))
}

/// A JSON list of the texts given.
pub fn texts(items: &[&str]) -> String {
    let items: Vec<String> = items.iter().map(|item| js(item)).collect();
    list(&items)
}

/// A manifest record of the line `vcf`.
pub fn record(vcf: &str) -> String {
    format!(
        r#"{{"dataset": "s01", "vcf": {}, "combinations": []}}"#,
        js(vcf)
    )
}

/// A manifest record of the line `vcf` with the further members `more`.
pub fn record_with(vcf: &str, more: &str) -> String {
    format!(
        r#"{{"dataset": "s01", "vcf": {}, "combinations": [], {more}}}"#,
        js(vcf)
    )
}

/// The member `records` of a manifest.
pub fn records(items: &[String]) -> String {
    format!(r#""records": {}"#, list(items))
}

/// The members every entry `classify` writes opens with: the key and its records.
pub fn key(location: &str, allele: &str, feature: (&str, &str), ordinals: &[usize]) -> String {
    let ordinals: Vec<String> = ordinals.iter().map(usize::to_string).collect();
    format!(
        r#""location": {}, "allele": {}, "feature": {}, "feature_type": {}, "record_indices": {}"#,
        js(location),
        js(allele),
        js(feature.0),
        js(feature.1),
        list(&ordinals)
    )
}

/// An entry of `divergences`.
pub fn divergence(key: &str, vep_set: &str, vep_rs_sets: &[&str], class: &str) -> String {
    format!(
        r#"{{{key}, "vep_consequence_set": {}, "vep_rs_consequence_sets": {}, "expected_divergence": {}}}"#,
        js(vep_set),
        texts(vep_rs_sets),
        js(class)
    )
}

/// An entry of `vep_rs_only_tuples`.
pub fn vep_rs_only(key: &str, vep_rs_sets: &[&str], class: &str) -> String {
    format!(
        r#"{{{key}, "vep_rs_consequence_sets": {}, "expected_divergence": {}}}"#,
        texts(vep_rs_sets),
        js(class)
    )
}

/// An entry of `field_divergences`: the consequence set, the field, both values, the class.
pub fn field_divergence(key: &str, set: &str, field: (&str, &str, &str), class: &str) -> String {
    format!(
        r#"{{{key}, "consequence_set": {}, "field": {}, "vep": {}, "vep_rs": {}, "expected_divergence": {}}}"#,
        js(set),
        js(field.0),
        js(field.1),
        js(field.2),
        js(class)
    )
}

/// The four members `classify` writes for the keys of a run, in its order.
pub fn written(divergences: &[String], vep_rs_only: &[String], summaries: (&str, &str)) -> String {
    format!(
        r#""divergences": {}, "vep_rs_only_tuples": {}, "divergence_summary": {}, "vep_rs_only_summary": {}"#,
        list(divergences),
        list(vep_rs_only),
        summaries.0,
        summaries.1
    )
}

/// The three members `classify` writes after those for the fields of `--field`.
pub fn fields_written(compared: &[&str], entries: &[String], summary: &str) -> String {
    format!(
        r#""fields_compared": {}, "field_divergences": {}, "field_divergence_summary": {summary}"#,
        texts(compared),
        list(entries)
    )
}

/// The members of a run with no difference on any key.
pub fn nothing_written() -> String {
    written(&[], &[], ("{}", "{}"))
}

/// The line of standard output of a run without `--field`.
pub fn counted(rows: (usize, &str), tuples: (usize, &str), skipped: &str) -> String {
    format!(
        "divergent VEP rows {} ({}); vep-rs-only tuples {} ({}); reference-skipped records {skipped}\n",
        rows.0, rows.1, tuples.0, tuples.1
    )
}
