// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the committed corpora document under the class
//! `refseq_edited_model_hgvsc_offset`, read from the corpus manifests, the
//! stored Ensembl VEP outputs and the corpus caches.
//!
//! The class is the `HGVSc` of a transcript model whose `_rna_edit` inserts or
//! deletes bases: Ensembl VEP numbers it from the unedited exons
//! (TranscriptVariationAllele.pm:2755, through a mapper rebuilt without the
//! edit, Dumper/Core.pm:164 and :169) and vep-rs against the edited model. A
//! RefSeq or merged cache prints `BAM_EDIT=OK` on such a model's rows. An
//! Ensembl-source cache prints no `BAM_EDIT`, so a row of its models is
//! documented from the cache record alone. The golden harness holds vep-rs's
//! output to the `vep_rs` value of every row read here.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{Read, Write};
use std::path::Path;

mod common;

const CLASS: &str = "refseq_edited_model_hgvsc_offset";

/// One `field_divergences` item a corpus manifest documents under [`CLASS`].
struct DocumentedRow {
    corpus: String,
    run: String,
    location: String,
    feature: String,
    field: String,
    vep: String,
    vep_rs: String,
    /// The `BAM_EDIT` value of the stored Ensembl VEP row, where it has one.
    bam_edit: Option<String>,
    /// The transcript's record in the corpus cache.
    record: serde_json::Value,
}

fn text(item: &serde_json::Value, key: &str) -> String {
    item[key].as_str().unwrap_or_default().to_string()
}

/// The manifest object whose lists document the run `name`: the manifest for
/// its own run, a flag set's entry where that carries `divergences` of its
/// own, else the manifest.
fn documenting_object<'a>(manifest: &'a serde_json::Value, name: &str) -> &'a serde_json::Value {
    manifest["flag_sets"]
        .as_array()
        .into_iter()
        .flatten()
        .find(|entry| entry["name"].as_str() == Some(name) && entry.get("divergences").is_some())
        .unwrap_or(manifest)
}

/// The `Extra` fields of the row of `stored`, a default-format output, at the
/// item's Location, Allele and Feature whose `HGVSc` is the item's `vep`
/// value, and `None` where no row at that key prints an `HGVSc`: the golden
/// harness holds the item to nothing in such a run. A key whose rows print an
/// `HGVSc` is refused unless exactly one of them carries the value.
fn stored_extra(
    stored: &str,
    label: &str,
    item: &serde_json::Value,
) -> Option<BTreeMap<String, String>> {
    let wanted = format!("HGVSc={}", text(item, "vep"));
    let printed: Vec<&str> = stored
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<&str>>())
        .filter(|columns| {
            columns.len() == 14
                && columns[1] == text(item, "location")
                && columns[2] == text(item, "allele")
                && columns[4] == text(item, "feature")
        })
        .map(|columns| columns[13])
        .filter(|extra| extra.split(';').any(|field| field.starts_with("HGVSc=")))
        .collect();
    if printed.is_empty() {
        return None;
    }
    let rows: Vec<&str> = printed
        .into_iter()
        .filter(|extra| extra.split(';').any(|field| field == wanted))
        .collect();
    assert_eq!(rows.len(), 1, "{label}: the stored rows carrying {wanted}");
    Some(
        rows[0]
            .split(';')
            .filter_map(|field| field.split_once('='))
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect(),
    )
}

/// The transcript records of one cache file, gzipped or plain: a JSON array of
/// records, or an object whose values are arrays of records, flattened.
fn cache_file_records(path: &Path) -> Vec<serde_json::Value> {
    let mut text = String::new();
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    if path.extension().is_some_and(|e| e == "gz") {
        flate2::read::MultiGzDecoder::new(file)
            .read_to_string(&mut text)
            .unwrap();
    } else {
        std::io::BufReader::new(file)
            .read_to_string(&mut text)
            .unwrap();
    }
    let parsed: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    match parsed {
        serde_json::Value::Array(records) => records,
        serde_json::Value::Object(ranges) => ranges
            .into_iter()
            .filter_map(|(_, records)| match records {
                serde_json::Value::Array(records) => Some(records),
                _ => None,
            })
            .flatten()
            .filter(|record| !record.is_null())
            .collect(),
        other => panic!("{}: {other}", path.display()),
    }
}

/// The record of the transcript `feature` in the corpus cache.
fn cache_record(corpus: &common::Corpus, feature: &str) -> serde_json::Value {
    let transcripts = corpus.cache_dir().join("transcripts");
    for sequence in fs::read_dir(&transcripts)
        .unwrap_or_else(|e| panic!("{}: {e}", transcripts.display()))
        .flatten()
    {
        for file in fs::read_dir(sequence.path()).unwrap().flatten() {
            let name = file.file_name().to_string_lossy().into_owned();
            if !(name.ends_with(".json") || name.ends_with(".json.gz")) {
                continue;
            }
            if let Some(record) = cache_file_records(&file.path())
                .into_iter()
                .find(|record| record["stable_id"].as_str() == Some(feature))
            {
                return record;
            }
        }
    }
    panic!(
        "{}/{}: no cache record of {feature}",
        corpus.release, corpus.name
    );
}

/// Every row the committed corpora document under [`CLASS`], run by run. A run
/// that does not compare the default format, or whose stored default output
/// prints no `HGVSc` at an item's key (`stored_extra`), has no row of the item.
fn documented_rows() -> Vec<DocumentedRow> {
    let mut rows = Vec::new();
    for corpus in common::corpora() {
        for run in common::runs(&corpus) {
            if !run.compares("default") {
                continue;
            }
            let lists = documenting_object(&corpus.manifest, &run.name);
            for item in lists["field_divergences"].as_array().into_iter().flatten() {
                if item["expected_divergence"].as_str() != Some(CLASS) {
                    continue;
                }
                let stored = run.expected("default.txt");
                let Some(mut extra) = stored_extra(&stored, &run.label(), item) else {
                    continue;
                };
                let feature = text(item, "feature");
                rows.push(DocumentedRow {
                    corpus: corpus.name.clone(),
                    run: run.name.clone(),
                    location: text(item, "location"),
                    field: text(item, "field"),
                    vep: text(item, "vep"),
                    vep_rs: text(item, "vep_rs"),
                    bam_edit: extra.remove("BAM_EDIT"),
                    record: cache_record(&corpus, &feature),
                    feature,
                });
            }
        }
    }
    rows
}

/// The distances a model's `_rna_edit` attributes put between its edited and
/// its unedited cDNA numbering: the running sums of the edits' length changes
/// in start order, as absolute values, 0 left out. A value is read as
/// `get_misalignment_offset` reads it (TranscriptVariationAllele.pm:1605-1618):
/// three fields substitute where their bases span exactly `end - start + 1`
/// and insert their bases otherwise, and two fields, or more than three,
/// delete `end - start + 1` bases.
fn edit_distances(record: &serde_json::Value) -> BTreeSet<i64> {
    let mut edits: Vec<(i64, i64)> = Vec::new();
    for attribute in record["attributes"].as_array().into_iter().flatten() {
        if !attribute["code"]
            .as_str()
            .is_some_and(|code| code.starts_with("_rna_edit"))
        {
            continue;
        }
        let value = attribute["value"].as_str().unwrap_or_default();
        let fields: Vec<&str> = value.trim_end_matches(' ').split(' ').collect();
        let (Some(start), Some(end)) = (
            fields.first().and_then(|f| f.parse::<i64>().ok()),
            fields.get(1).and_then(|f| f.parse::<i64>().ok()),
        ) else {
            continue;
        };
        match fields[..] {
            [_, _, bases] if end - start + 1 == bases.len() as i64 => {}
            [_, _, bases] => edits.push((start, bases.len() as i64)),
            _ => edits.push((start, -1 - (end - start))),
        }
    }
    edits.sort_unstable();
    let mut running = 0;
    let mut distances = BTreeSet::new();
    for (_, length_change) in edits {
        running += length_change;
        if running != 0 {
            distances.insert(running.abs());
        }
    }
    distances
}

/// The length of the cDNA coordinate token at the start of `body` (`-71`,
/// `250+601`, `*12-3`), 0 where none starts there.
fn coordinate_length(body: &[u8]) -> usize {
    let digits = |from: usize| {
        body[from..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count()
    };
    let mut at = usize::from(body.first() == Some(&b'*'));
    at += usize::from(body.get(at) == Some(&b'-'));
    let exon = digits(at);
    if exon == 0 {
        return 0;
    }
    at += exon;
    if matches!(body.get(at), Some(b'+' | b'-')) && digits(at + 1) > 0 {
        at += 1 + digits(at + 1);
    }
    at
}

/// An HGVSc body with every cDNA coordinate token replaced by `#`, and the
/// signed integers of those tokens: `250+601` gives 250 and 601, `310-763`
/// gives 310 and -763.
fn coordinates_masked(body: &str) -> (String, Vec<i64>) {
    let bytes = body.as_bytes();
    let (mut masked, mut integers) = (String::new(), Vec::new());
    let mut at = 0;
    while at < bytes.len() {
        let length = coordinate_length(&bytes[at..]);
        if length == 0 {
            masked.push(bytes[at] as char);
            at += 1;
            continue;
        }
        let token = body[at..at + length].trim_start_matches('*');
        let split = token[1..]
            .find(['+', '-'])
            .map_or(token.len(), |offset| offset + 1);
        integers.push(token[..split].parse().unwrap());
        if split < token.len() {
            integers.push(token[split..].trim_start_matches('+').parse().unwrap());
        }
        masked.push('#');
        at += length;
    }
    (masked, integers)
}

/// A page below `docs/`, its whitespace collapsed so that a sentence wrapped
/// over several lines reads as one.
fn page(path: &str) -> String {
    fs::read_to_string(format!("{}/../../docs/{path}", env!("CARGO_MANIFEST_DIR")))
        .unwrap_or_else(|e| panic!("docs/{path}: {e}"))
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// Every Markdown page below `dir`, by its path.
fn pages_below(dir: &Path, out: &mut Vec<std::path::PathBuf>) {
    for entry in fs::read_dir(dir).unwrap().flatten() {
        let path = entry.path();
        if path.is_dir() {
            pages_below(&path, out);
        } else if path.extension().is_some_and(|e| e == "md") {
            out.push(path);
        }
    }
}

#[test]
fn coordinate_tokens_are_read_as_the_class_reads_them() {
    assert_eq!(
        coordinates_masked("c.3546_3553+23del"),
        ("c.#_#del".to_string(), vec![3546, 3553, 23])
    );
    assert_eq!(
        coordinates_masked("c.310-763G>T"),
        ("c.#G>T".to_string(), vec![310, -763])
    );
    assert_eq!(
        coordinates_masked("c.*12-3_-71del"),
        ("c.#_#del".to_string(), vec![12, -3, -71])
    );
}

/// Three fields substitute or insert, and two fields or more than three
/// delete (TranscriptVariationAllele.pm:1606-1607). The length changes are
/// summed in start order, and a sum of 0 is no distance.
#[test]
fn edit_values_are_read_as_the_class_reads_them() {
    let distances = |values: &[&str]| {
        let attributes: Vec<serde_json::Value> = values
            .iter()
            .map(|value| serde_json::json!({ "code": "_rna_edit", "value": value }))
            .collect();
        edit_distances(&serde_json::json!({ "attributes": attributes }))
    };
    assert_eq!(distances(&["3387 3395 "]), BTreeSet::from([9]));
    assert_eq!(distances(&["647 646 G"]), BTreeSet::from([1]));
    assert_eq!(distances(&["2023 2023 T"]), BTreeSet::new());
    assert_eq!(distances(&["100 104", "50 49 AC"]), BTreeSet::from([2, 3]));
    assert_eq!(distances(&["5 4 AC", "9 10"]), BTreeSet::from([2]));
    assert_eq!(distances(&["10 12 A C"]), BTreeSet::from([3]));
}

/// A stored row documents an item where its key prints an `HGVSc`: a key
/// without one is no row of the run, and a key whose rows print one is refused
/// unless exactly one carries the item's `vep` value.
#[test]
fn a_stored_row_is_read_where_its_key_prints_an_hgvsc() {
    let item = serde_json::json!({
        "location": "1:100", "allele": "-", "feature": "T1", "vep": "T1:c.5-3del",
    });
    let row = |location: &str, extra: &str| {
        let columns = [
            "id",
            location,
            "-",
            "G1",
            "T1",
            "Transcript",
            "intron_variant",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            extra,
        ];
        columns.join("\t") + "\n"
    };
    let read = |stored: String| stored_extra(&stored, "made-up", &item);

    let documented = row("1:100", "IMPACT=MODIFIER;BAM_EDIT=OK;HGVSc=T1:c.5-3del")
        + &row("1:200", "IMPACT=MODIFIER;HGVSc=T1:c.9del");
    let extra = read(format!("## made-up output\n{documented}")).unwrap();
    assert_eq!(extra["BAM_EDIT"], "OK");
    assert_eq!(extra["HGVSc"], "T1:c.5-3del");

    let without = row("1:100", "IMPACT=MODIFIER") + &row("1:200", "HGVSc=T1:c.5-3del");
    assert_eq!(read(without), None);

    let another = row("1:100", "IMPACT=MODIFIER;HGVSc=T1:c.6-3del");
    let twice = row("1:100", "HGVSc=T1:c.5-3del") + &row("1:100", "HGVSc=T1:c.5-3del");
    for refused in [another, twice] {
        assert!(std::panic::catch_unwind(|| read(refused)).is_err());
    }
}

/// A shard files its records as a JSON array or, range by range, as the array
/// values of an object, and both forms read as one list of records.
#[test]
fn a_cache_file_is_read_in_both_shard_forms() {
    let dir = tempfile::tempdir().unwrap();
    let ids = |path: &Path| -> Vec<String> {
        cache_file_records(path)
            .iter()
            .map(|record| text(record, "stable_id"))
            .collect()
    };
    let array = dir.path().join("1-1000.json");
    fs::write(&array, r#"[{"stable_id": "T1"}, {"stable_id": "T2"}]"#).unwrap();
    assert_eq!(ids(&array), ["T1", "T2"]);

    let object = dir.path().join("1001-3000.json.gz");
    let ranges =
        r#"{"1001-2000": [{"stable_id": "T3"}, null], "2001-3000": [{"stable_id": "T4"}]}"#;
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&object).unwrap(),
        flate2::Compression::default(),
    );
    encoder.write_all(ranges.as_bytes()).unwrap();
    encoder.finish().unwrap();
    assert_eq!(ids(&object), ["T3", "T4"]);
}

/// Every documented row is an `HGVSc` whose two notations name one reference
/// and differ in their coordinates alone, by a distance the edits of the
/// transcript's cache record put between the two numberings: the corpus cache
/// files the model with an `_rna_edit` that inserts or deletes bases.
#[test]
fn every_documented_row_is_a_renumbered_hgvsc_on_a_model_whose_edit_changes_its_length() {
    let rows = documented_rows();
    assert!(!rows.is_empty(), "no corpus documents a row of {CLASS}");
    for row in &rows {
        let label = format!("{} [{}] {}", row.corpus, row.run, row.feature);
        assert_eq!(row.field, "HGVSc", "{label}");
        let (vep_name, vep_body) = row.vep.split_once(':').expect(&label);
        let (vep_rs_name, vep_rs_body) = row.vep_rs.split_once(':').expect(&label);
        assert_eq!(vep_name, vep_rs_name, "{label}: the reference named");
        let (vep_masked, vep_integers) = coordinates_masked(vep_body);
        let (vep_rs_masked, vep_rs_integers) = coordinates_masked(vep_rs_body);
        assert_ne!(vep_body, vep_rs_body, "{label}");
        assert_eq!(
            vep_masked, vep_rs_masked,
            "{label}: all but the coordinates"
        );
        let distances = edit_distances(&row.record);
        assert!(
            !distances.is_empty(),
            "{label}: the cache record carries no edit that inserts or deletes bases"
        );
        assert!(
            vep_integers.iter().any(|a| vep_rs_integers
                .iter()
                .any(|b| distances.contains(&(a - b).abs()))),
            "{label}: {} and {} are not one of {distances:?} apart",
            row.vep,
            row.vep_rs
        );
    }
}

/// The `BAM_EDIT` column follows the cache record's alignment status, and the
/// class does not need it: a documented row of a model the cache files
/// without a status prints no `BAM_EDIT`, as the yellow fever mosquito
/// corpus's row does in both of its flag sets under `--hgvs`.
#[test]
fn a_documented_row_of_a_model_without_an_alignment_status_prints_no_bam_edit() {
    let rows = documented_rows();
    for row in &rows {
        let label = format!("{} [{}] {}", row.corpus, row.run, row.feature);
        match row.record["bam_edit_status"].as_str() {
            Some("ok") => assert_eq!(row.bam_edit.as_deref(), Some("OK"), "{label}"),
            None => assert_eq!(row.bam_edit, None, "{label}"),
            Some(other) => panic!("{label}: a documented row on a model whose status is {other}"),
        }
    }
    let without: BTreeSet<(&str, &str, &str)> = rows
        .iter()
        .filter(|row| row.bam_edit.is_none())
        .map(|row| (row.corpus.as_str(), row.run.as_str(), row.feature.as_str()))
        .collect();
    for run in ["everything", "hgvs-sift"] {
        assert!(
            without.contains(&("AaegL5-aedes", run, "AAEL025507-RA")),
            "{without:?}"
        );
    }
}

/// The yellow fever mosquito corpus's row, derived from its cache record. The
/// model is on the reverse strand, and its `_rna_edit` `302 302` deletes one
/// base of the first exon. The record lies in the intron, 763 bases before the
/// last exon's first base and 2,283 after the first exon's last. The cached
/// mapper, which folds the edit in, puts that first base at cDNA 373, coding
/// base 309: vep-rs's `c.309-763G>T`. The exons' genomic lengths put it at
/// cDNA 374, coding base 310: Ensembl VEP's `c.310-763G>T`.
#[test]
fn the_mosquito_row_is_numbered_from_the_exon_start_of_the_edited_model() {
    let rows = documented_rows();
    let row = rows
        .iter()
        .find(|row| row.corpus == "AaegL5-aedes" && row.run == "hgvs-sift")
        .expect("the hgvs-sift row of AaegL5-aedes");
    assert_eq!(row.location, "NIGP01001870:3784");
    assert_eq!(row.vep, "AAEL025507-RA:c.310-763G>T");
    assert_eq!(row.vep_rs, "AAEL025507-RA:c.309-763G>T");

    let record = &row.record;
    let number = |value: &serde_json::Value| value.as_i64().unwrap();
    assert_eq!(number(&record["strand"]), -1);
    assert_eq!(edit_distances(record), BTreeSet::from([1]));
    let mut exons: Vec<(i64, i64)> = record["exons"]
        .as_array()
        .unwrap()
        .iter()
        .map(|exon| (number(&exon["start"]), number(&exon["end"])))
        .collect();
    exons.sort_unstable_by_key(|&(start, _)| std::cmp::Reverse(start));
    assert_eq!(exons, [(6067, 6439), (3018, 3021)]);
    let (first, last) = (exons[0], exons[1]);
    let position: i64 = row.location.split_once(':').unwrap().1.parse().unwrap();
    let (after_first, before_last) = (first.0 - position, position - last.1);
    assert_eq!((after_first, before_last), (2283, 763));

    let coding_start = number(&record["cdna_coding_start"]);
    let pairs = record["variation_effect_feature_cache"]["mapper"]["pairs"]
        .as_array()
        .unwrap();
    let edited = pairs
        .iter()
        .find(|pair| number(&pair["to_start"]) <= last.1 && last.1 <= number(&pair["to_end"]))
        .map(|pair| number(&pair["from_start"]) + number(&pair["to_end"]) - last.1)
        .unwrap();
    let unedited = first.1 - first.0 + 2;
    assert_eq!((edited, unedited), (373, 374));
    let notation = |cdna: i64| {
        format!(
            "{}:c.{}-{before_last}G>T",
            row.feature,
            cdna - coding_start + 1
        )
    };
    assert_eq!(notation(edited), row.vep_rs);
    assert_eq!(notation(unedited), row.vep);
    assert_eq!(
        number(&record["cdna_coding_end"]),
        edited + last.1 - last.0,
        "the model's coding sequence ends on the last exon's last base"
    );
}

/// The class's page says how a documented row is named and how a row that
/// prints no `BAM_EDIT` is set aside, and the rows it gives as its examples
/// are documented as it says. No page says that one corpus documents the
/// class or that every documented row carries the column.
#[test]
fn the_page_says_how_a_row_is_named_and_how_one_without_bam_edit_is_set_aside() {
    let text = page("divergences/refseq-edited-model-hgvsc.md");
    for phrase in [
        "A model of an Ensembl-source cache whose `_rna_edit` inserts or deletes bases prints \
         no `BAM_EDIT`, so the comparator's own rule leaves its row a difference and a corpus \
         manifest alone sets it aside.",
        "The corpus tool names the class from the corpus cache, not from the `BAM_EDIT` column.",
        "A documented row of a RefSeq or merged cache carries `BAM_EDIT=OK`, as the two rows \
         above do in the GRCh38 RefSeq corpus's `hgvs` flag set.",
        "A documented row of an Ensembl-source cache carries no `BAM_EDIT`, as the row of",
        "On every documented row only the `HGVSc` changes.",
    ] {
        assert!(text.contains(phrase), "the page lacks {phrase:?}");
    }
    let rows = documented_rows();
    let with_the_column: BTreeSet<&str> = rows
        .iter()
        .filter(|row| row.corpus == "GRCh38-refseq" && row.run == "hgvs")
        .filter(|row| row.bam_edit.as_deref() == Some("OK"))
        .map(|row| row.feature.as_str())
        .collect();
    assert_eq!(
        with_the_column,
        BTreeSet::from(["NM_030582.4", "NM_130444.3"]),
        "the two rows of the page's first record"
    );

    let mut pages = Vec::new();
    pages_below(
        Path::new(&format!("{}/../../docs", env!("CARGO_MANIFEST_DIR"))),
        &mut pages,
    );
    for path in pages {
        let text = fs::read_to_string(&path)
            .unwrap()
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ");
        for refused in [
            "the one corpus run under `--refseq --hgvs`",
            "Every changed row carries `BAM_EDIT=OK`",
        ] {
            assert!(
                !text.contains(refused),
                "{} states {refused:?}",
                path.display()
            );
        }
    }
}
