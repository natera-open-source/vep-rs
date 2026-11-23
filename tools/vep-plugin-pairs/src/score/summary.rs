// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The summary a scored run prints on standard output: each entry's counts, one line per
//! class of it, the keys no source claims, the sources not populated, the context a run did
//! not name for an applied class, and the SpliceAI refusal counts.
//!
//! A character the stream's encoding lacks is written as its backslash escape. The encoding
//! is the locale's as CPython reads it: UTF-8 under the `C` and `POSIX` locales and where
//! the environment names a locale the system lacks, else the locale's codeset.

use super::reads::{
    DOWNSTREAM_PLUGIN, GWAS_PLUGIN, MAXENTSCAN_PLUGIN, REVEL_PLUGIN, SPLICEAI_PLUGIN,
    UTRANNOTATOR_PLUGIN,
};
use super::text::Encoding;
use serde_json::Value;
use std::ffi::{CStr, CString};
use std::io::{self, Write};

/// The class that reads the SpliceAI score files.
const SPLICEAI_MULTI_ENTRY_CELL: &str = "spliceai_multi_entry_cell";

fn text(value: &Value) -> &str {
    value.as_str().unwrap_or_default()
}

fn items(value: &Value) -> impl Iterator<Item = &Value> {
    value.as_array().into_iter().flatten()
}

fn entries(value: &Value) -> impl Iterator<Item = (&String, &Value)> {
    value.as_object().into_iter().flatten()
}

/// The summary's lines of `report`. `populations` gives a class the keys it can reach, a
/// count above which is marked; `refusal_lines` are the lines of the SpliceAI refusal
/// counts, printed last.
pub fn lines(
    report: &Value,
    populations: &[(String, u64)],
    refusal_lines: &[String],
) -> Vec<String> {
    let mut lines = Vec::new();
    let further = |class: &str| -> String {
        let columns: Vec<&str> = items(&report["classes"])
            .find(|item| text(&item["class"]) == class)
            .map(|item| items(&item["further_columns"]).map(text).collect())
            .unwrap_or_default();
        if columns.is_empty() {
            String::new()
        } else {
            format!(" [{}]", columns.join(", "))
        }
    };
    for entry in items(&report["scores"]) {
        let row = format!(
            "{}\t{}\t{}",
            text(&entry["source"]),
            text(&entry["engine"]),
            text(&entry["projection"])
        );
        let counts: Vec<String> = entries(entry)
            .filter(|(_, value)| value.is_number())
            .map(|(key, value)| format!("{key}={value}"))
            .collect();
        lines.push(format!("{row}\t{}", counts.join(" ")));
        for (class, count) in entries(&entry["masked"]) {
            let count = count.as_u64().unwrap_or(0);
            let population = populations
                .iter()
                .find(|(name, _)| name == class)
                .map(|(_, keys)| *keys);
            let note = match population {
                Some(keys) if count > keys => {
                    format!(" (above its reachable population of {keys})")
                }
                _ if count == 0 => " (none)".to_string(),
                _ => String::new(),
            };
            lines.push(format!(
                "{row}\tclass {class} set aside {count}{note}{}",
                further(class)
            ));
        }
        for (class, count) in entries(&entry["not_applied"]) {
            lines.push(format!(
                "{row}\tclass {class} is not applied and would set aside {count}: raw counts only{}",
                further(class)
            ));
        }
    }
    let engines = items(&report["engines"])
        .map(|engine| (format!("the engine {}", text(&engine["id"])), engine));
    for (whose, described) in
        std::iter::once(("the reference".to_string(), &report["reference"])).chain(engines)
    {
        let keys: Vec<String> = entries(&described["unclaimed_keys"])
            .map(|(name, times)| format!("{name} ({times})"))
            .collect();
        if !keys.is_empty() {
            lines.push(format!(
                "keys of {whose} that no source claims: {}",
                keys.join(", ")
            ));
        }
    }
    let not_populated: Vec<&str> = items(&report["not_populated"]).map(text).collect();
    if !not_populated.is_empty() {
        lines.push(format!(
            "not populated in the reference: {}",
            not_populated.join(", ")
        ));
    }
    let loaded = |plugin: &str| items(&report["plugins"]).any(|name| text(name) == plugin);
    for (plugin, key, option, clause) in [
        (
            GWAS_PLUGIN,
            "gwas_catalog",
            "--gwas-catalog",
            "that reads the catalog sets nothing aside",
        ),
        (
            REVEL_PLUGIN,
            "revel_table",
            "--revel-table",
            "that reads the table sets nothing aside",
        ),
        (
            UTRANNOTATOR_PLUGIN,
            "utrannotator_sweep",
            "--utrannotator-sweep",
            "sets aside no row it reads the seeded runs for",
        ),
        (
            DOWNSTREAM_PLUGIN,
            "json_cache",
            "--json-cache",
            "that reads the candidate cache sets nothing aside",
        ),
        (
            MAXENTSCAN_PLUGIN,
            "json_cache",
            "--json-cache",
            "that reads the candidate cache sets nothing aside",
        ),
    ] {
        let applied = items(&report["classes"]).any(|item| {
            item["applied"] == Value::Bool(true)
                && item["plugins"]
                    .as_object()
                    .is_some_and(|plugins| plugins.contains_key(plugin))
        });
        if report[key].is_null() && applied && loaded(plugin) {
            lines.push(format!("no {option}: a class of {plugin} {clause}"));
        }
    }
    let reads_score_files = items(&report["classes"]).any(|item| {
        item["applied"] == Value::Bool(true) && text(&item["class"]) == SPLICEAI_MULTI_ENTRY_CELL
    });
    if reads_score_files && loaded(SPLICEAI_PLUGIN) {
        if report["spliceai_table"].is_null() {
            lines.push(format!(
                "no --spliceai-table: {SPLICEAI_MULTI_ENTRY_CELL} sets nothing aside"
            ));
        } else if report["json_cache"].is_null() {
            lines.push(format!(
                "no --json-cache: {SPLICEAI_MULTI_ENTRY_CELL} sets aside no row that writes no SYMBOL"
            ));
        }
    }
    lines.extend(refusal_lines.iter().cloned());
    lines
}

/// The encoding of the process's standard streams.
pub fn stream_encoding() -> Encoding {
    let named = ["LC_ALL", "LC_CTYPE", "LANG"]
        .iter()
        .filter_map(|name| std::env::var(name).ok())
        .find(|value| !value.is_empty())
        .unwrap_or_default();
    if matches!(named.as_str(), "" | "C" | "POSIX") {
        return Encoding::Utf8;
    }
    let Ok(name) = CString::new(named) else {
        return Encoding::Utf8;
    };
    // SAFETY: `name` is a NUL-terminated string; the locale made is read and freed here, and
    // the codeset's text is copied before the locale is freed.
    let codeset = unsafe {
        let locale = libc::newlocale(libc::LC_CTYPE_MASK, name.as_ptr(), std::ptr::null_mut());
        if locale.is_null() {
            return Encoding::Utf8;
        }
        let codeset = CStr::from_ptr(libc::nl_langinfo_l(libc::CODESET, locale))
            .to_string_lossy()
            .into_owned();
        libc::freelocale(locale);
        codeset
    };
    match codeset
        .to_ascii_uppercase()
        .replace(['-', '_'], "")
        .as_str()
    {
        "UTF8" => Encoding::Utf8,
        "ISO88591" | "LATIN1" => Encoding::Latin1,
        _ => Encoding::Ascii,
    }
}

/// A stream that takes text and writes it in an encoding, a character the encoding lacks
/// as its backslash escape. Every write is whole text.
pub struct Encoded<W: Write> {
    stream: W,
    encoding: Encoding,
}

impl<W: Write> Encoded<W> {
    pub fn new(stream: W, encoding: Encoding) -> Encoded<W> {
        Encoded { stream, encoding }
    }
}

impl<W: Write> Write for Encoded<W> {
    fn write(&mut self, text: &[u8]) -> io::Result<usize> {
        match self.encoding {
            Encoding::Utf8 => self.stream.write_all(text)?,
            encoding => self
                .stream
                .write_all(&encoding.encode(&String::from_utf8_lossy(text)))?,
        }
        Ok(text.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.stream.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn report() -> Value {
        json!({
            "plugins": ["Planted", "GWAS", "MaxEntScan", "LOEUF", "SpliceAI"],
            "gwas_catalog": null,
            "revel_table": null,
            "utrannotator_sweep": null,
            "json_cache": null,
            "spliceai_table": null,
            "classes": [
                {"class": "planted_two_numbers", "plugins": {"MaxEntScan": ["MaxEntScan_*"]}, "applied": true, "further_columns": ["met", "-"]},
                {"class": "planted_reference_value", "plugins": {"GWAS": ["GWAS_*"]}, "applied": true, "further_columns": []},
                {"class": "planted_engine_zero", "plugins": {"LOEUF": ["LOEUF"]}, "applied": false, "further_columns": ["met", "-"]},
            ],
            "reference": {"path": "r", "unclaimed_keys": {"Custom": 1, "Custom_SIG": 1}},
            "engines": [{"id": "engine", "unclaimed_keys": {"Kéy": 2}}, {"id": "other", "unclaimed_keys": {}}],
            "scores": [
                {"source": "GWAS", "engine": "engine", "projection": "feature", "perl": 3, "rust": 3, "intersection": 3,
                 "adj_perl": 3, "adj_rust": 3, "adj_intersection": 3, "masked": {"planted_reference_value": 0}, "not_applied": {},
                 "text_only_differences": 0, "residual_rows": 0, "shared_keys": 0},
                {"source": "LOEUF", "engine": "engine", "projection": "feature", "perl": 1, "rust": 2, "intersection": 1,
                 "masked": {}, "not_applied": {"planted_engine_zero": 1},
                 "text_only_differences": 0, "residual_rows": 1, "shared_keys": 0},
                {"source": "MaxEntScan", "engine": "engine", "projection": "feature", "perl": 2, "rust": 2, "intersection": 1,
                 "adj_perl": 1, "adj_rust": 1, "adj_intersection": 1, "masked": {"planted_two_numbers": 1}, "not_applied": {},
                 "text_only_differences": 0, "residual_rows": 0, "shared_keys": 0},
            ],
            "not_populated": ["Planted", "SpliceAI"],
        })
    }

    #[test]
    fn summary_prints_each_class_count_and_marks_none_and_above_population() {
        let populations = [
            ("planted_two_numbers".to_string(), 0),
            ("planted_reference_value".to_string(), 5),
        ];
        let lines = lines(&report(), &populations, &["a refusal line".to_string()]);
        assert_eq!(
            lines,
            [
                "GWAS\tengine\tfeature\tperl=3 rust=3 intersection=3 adj_perl=3 adj_rust=3 adj_intersection=3 text_only_differences=0 residual_rows=0 shared_keys=0",
                "GWAS\tengine\tfeature\tclass planted_reference_value set aside 0 (none)",
                "LOEUF\tengine\tfeature\tperl=1 rust=2 intersection=1 text_only_differences=0 residual_rows=1 shared_keys=0",
                "LOEUF\tengine\tfeature\tclass planted_engine_zero is not applied and would set aside 1: raw counts only [met, -]",
                "MaxEntScan\tengine\tfeature\tperl=2 rust=2 intersection=1 adj_perl=1 adj_rust=1 adj_intersection=1 text_only_differences=0 residual_rows=0 shared_keys=0",
                "MaxEntScan\tengine\tfeature\tclass planted_two_numbers set aside 1 (above its reachable population of 0) [met, -]",
                "keys of the reference that no source claims: Custom (1), Custom_SIG (1)",
                "keys of the engine engine that no source claims: Kéy (2)",
                "not populated in the reference: Planted, SpliceAI",
                "no --gwas-catalog: a class of GWAS that reads the catalog sets nothing aside",
                "no --json-cache: a class of MaxEntScan that reads the candidate cache sets nothing aside",
                "a refusal line",
            ]
        );
    }

    #[test]
    fn summary_names_the_score_files_or_the_cache_a_run_did_not_give() {
        let mut report = report();
        report["classes"] = json!([{"class": "spliceai_multi_entry_cell", "plugins": {"SpliceAI": ["SpliceAI_pred"]}, "applied": true, "further_columns": []}]);
        report["scores"] = json!([]);
        report["not_populated"] = json!([]);
        report["reference"]["unclaimed_keys"] = json!({});
        report["engines"] = json!([]);
        assert_eq!(
            lines(&report, &[], &[]),
            ["no --spliceai-table: spliceai_multi_entry_cell sets nothing aside"]
        );
        report["spliceai_table"] = json!([{"path": "scores.synthetic.vcf.gz"}]);
        assert_eq!(
            lines(&report, &[], &[]),
            ["no --json-cache: spliceai_multi_entry_cell sets aside no row that writes no SYMBOL"]
        );
        report["json_cache"] = json!({"path": "cache"});
        assert!(lines(&report, &[], &[]).is_empty());
        report["spliceai_table"] = Value::Null;
        report["plugins"] = json!(["GWAS"]);
        assert!(lines(&report, &[], &[]).is_empty());
    }

    #[test]
    fn summary_escapes_what_the_stream_cannot_encode() {
        let line = "keys of the reference that no source claims: Kéy (1)\n";
        let mut written = Vec::new();
        Encoded::new(&mut written, Encoding::Ascii)
            .write_all(line.as_bytes())
            .unwrap();
        assert_eq!(
            written,
            b"keys of the reference that no source claims: K\\xe9y (1)\n"
        );
        let mut written = Vec::new();
        Encoded::new(&mut written, Encoding::Utf8)
            .write_all(line.as_bytes())
            .unwrap();
        assert_eq!(written, line.as_bytes());
    }
}
