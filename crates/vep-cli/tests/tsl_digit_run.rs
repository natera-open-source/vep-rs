// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `TSL` attribute value outside the one-digit vocabulary of an Ensembl gene
//! set, as a GFF3 given to the cache builder can carry it: the level is the
//! digit run after `tsl` as text in the text formats and the number it spells
//! in the JSON format (`OutputFactory.pm:1618-1622` with `numberify`), and the
//! pick order's `tsl` category ranks the run as a number and ties it as text
//! (`:750-755`, `:789-796`).

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Stdio;

use vep_core::cache_manifest::CacheManifest;

mod common;

/// The native-builder corpus, whose cache carries `TSL` attributes alone (no
/// top-level `tsl` key).
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-native-builder")
}

fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// A copy of the corpus cache under `dir/cache` whose listed transcripts
/// carry the one `TSL` attribute `tsl<level>` (`NA` for `tslNA`), or none for
/// `None`, with the checksum manifest rewritten for the edited shards.
fn cache_with_levels(dir: &Path, levels: &[(&str, Option<&str>)]) -> PathBuf {
    let cache = dir.join("cache");
    copy_tree(&corpus().join("json_cache"), &cache);
    for shard in fs::read_dir(cache.join("transcripts/21")).unwrap() {
        let shard = shard.unwrap().path();
        let mut text = String::new();
        flate2::read::GzDecoder::new(fs::File::open(&shard).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        let mut records: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
        let mut edited = false;
        for record in &mut records {
            let Some((_, level)) = levels.iter().find(|(id, _)| record["stable_id"] == *id) else {
                continue;
            };
            let attributes = record["attributes"].as_array_mut().unwrap();
            attributes.retain(|a| a["code"] != "TSL");
            if let Some(level) = level {
                let value = format!("tsl{level}");
                attributes.push(serde_json::json!({"code": "TSL", "value": value}));
            }
            edited = true;
        }
        if !edited {
            continue;
        }
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(&shard).unwrap(),
            flate2::Compression::default(),
        );
        encoder
            .write_all(serde_json::to_string(&records).unwrap().as_bytes())
            .unwrap();
        encoder.finish().unwrap();
    }
    CacheManifest::write(&cache).unwrap();
    cache
}

/// Runs `vep` on the one record `21 <pos> . T C` against `cache` with `flags`
/// and returns the output file's text.
fn run(dir: &Path, cache: &Path, pos: u64, flags: &[&str]) -> String {
    let input = dir.join("one.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t{pos}\t.\tT\tC\t.\t.\t.\n"),
    )
    .unwrap();
    let out = dir.join("out.txt");
    let _ = fs::remove_file(&out);
    let output = common::isolated_vep(dir)
        .stdin(Stdio::null())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args([
            "--offline",
            "--json_cache",
            cache.to_str().unwrap(),
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(flags)
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out).unwrap()
}

/// The `TSL` value the default format prints for `feature`, `None` for no
/// `TSL=` in its `Extra` column.
fn default_tsl(text: &str, feature: &str) -> Option<String> {
    let row = data_row(text, feature);
    let extra = row.split('\t').nth(13).unwrap();
    extra
        .split(';')
        .find_map(|kv| kv.strip_prefix("TSL=").map(str::to_string))
}

/// The `TSL` column of the tab format for `feature`, `None` for `-`.
fn tab_tsl(text: &str, feature: &str) -> Option<String> {
    let columns: Vec<&str> = text
        .lines()
        .find(|l| l.starts_with('#') && !l.starts_with("##"))
        .unwrap()
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let column = columns.iter().position(|c| *c == "TSL").unwrap();
    let cell = data_row(text, feature).split('\t').nth(column).unwrap();
    (cell != "-").then(|| cell.to_string())
}

/// The `TSL` field of the VCF format's `CSQ` entry for `feature`, `None` when
/// empty.
fn vcf_tsl(text: &str, feature: &str) -> Option<String> {
    let format = text
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .unwrap();
    let fields: Vec<&str> = format
        .split("Format: ")
        .nth(1)
        .unwrap()
        .trim_end_matches("\">")
        .split('|')
        .collect();
    let tsl = fields.iter().position(|f| *f == "TSL").unwrap();
    let feature_field = fields.iter().position(|f| *f == "Feature").unwrap();
    let record = text.lines().find(|l| !l.starts_with('#')).unwrap();
    let csq = record
        .split('\t')
        .nth(7)
        .unwrap()
        .split(';')
        .find_map(|kv| kv.strip_prefix("CSQ="))
        .unwrap();
    let entry = csq
        .split(',')
        .find(|e| e.split('|').nth(feature_field) == Some(feature))
        .unwrap_or_else(|| panic!("a CSQ entry for {feature} in {csq}"));
    let cell = entry.split('|').nth(tsl).unwrap();
    (!cell.is_empty()).then(|| cell.to_string())
}

/// The `tsl` key of the JSON consequence for `feature`.
fn json_tsl(text: &str, feature: &str) -> Option<serde_json::Value> {
    let record: serde_json::Value = serde_json::from_str(text.lines().next().unwrap()).unwrap();
    let consequence = record["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["transcript_id"] == feature)
        .unwrap_or_else(|| panic!("a consequence for {feature} in {record}"));
    consequence.get("tsl").cloned()
}

/// The data line of the default or tab format whose `Feature` is `feature`.
fn data_row<'a>(text: &'a str, feature: &str) -> &'a str {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .find(|l| l.split('\t').nth(4) == Some(feature))
        .unwrap_or_else(|| panic!("a row for {feature} in {text}"))
}

/// The `Feature` ids of the data rows left, in order.
fn features(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(4).unwrap().to_string())
        .collect()
}

const TMPRSS15: &str = "ENST00000284885";
const TMPRSS15_RECORD: u64 = 18398250;
const CHODL_RECORD: u64 = 18256600;

/// `$hash->{TSL} = $1 if $1` keeps the matched digit run as a string, so the
/// default, tab and VCF formats print `007`, `00` (a true string), `255`,
/// `256` and `999` as the attribute spelt them, and the JSON format writes
/// the number each spells (7, 0, 255, 256, 999); `tsl0`, `tslNA` and a
/// transcript without the attribute print nothing in every format.
#[test]
fn a_level_outside_the_one_digit_vocabulary_prints_as_its_digit_run() {
    let dir = tempfile::tempdir().unwrap();
    for (level, number) in [
        ("1", 1u64),
        ("007", 7),
        ("00", 0),
        ("255", 255),
        ("256", 256),
        ("999", 999),
    ] {
        let case = dir.path().join(level);
        fs::create_dir_all(&case).unwrap();
        let cache = cache_with_levels(&case, &[(TMPRSS15, Some(level))]);
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl"]);
        assert_eq!(
            default_tsl(&text, TMPRSS15).as_deref(),
            Some(level),
            "{text}"
        );
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl", "--tab"]);
        assert_eq!(tab_tsl(&text, TMPRSS15).as_deref(), Some(level), "{text}");
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl", "--vcf"]);
        assert_eq!(vcf_tsl(&text, TMPRSS15).as_deref(), Some(level), "{text}");
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl", "--json"]);
        assert_eq!(
            json_tsl(&text, TMPRSS15),
            Some(serde_json::json!(number)),
            "{text}"
        );
    }
    for (name, level) in [("0", Some("0")), ("NA", Some("NA")), ("none", None)] {
        let case = dir.path().join(name);
        fs::create_dir_all(&case).unwrap();
        let cache = cache_with_levels(&case, &[(TMPRSS15, level)]);
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl"]);
        assert_eq!(default_tsl(&text, TMPRSS15), None, "{text}");
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl", "--tab"]);
        assert_eq!(tab_tsl(&text, TMPRSS15), None, "{text}");
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl", "--vcf"]);
        assert_eq!(vcf_tsl(&text, TMPRSS15), None, "{text}");
        let text = run(&case, &cache, TMPRSS15_RECORD, &["--tsl", "--json"]);
        assert_eq!(json_tsl(&text, TMPRSS15), None, "{text}");
    }
}

/// The nine CHODL transcripts over 21:18256600, in Ensembl VEP's row order
/// (the stable ids, the cache carrying no database ids), with translateable
/// lengths 822, 711, 699, 699, 711, 711, 765, 753 and 822.
const CHODL: [&str; 9] = [
    "ENST00000299295",
    "ENST00000338326",
    "ENST00000400127",
    "ENST00000400128",
    "ENST00000400131",
    "ENST00000400135",
    "ENST00000543733",
    "ENST00001123871",
    "ENST00001131842",
];

/// Every CHODL transcript without a `TSL` attribute but the ones `set` names.
fn chodl_levels<'a>(set: &[(&'a str, &'a str)]) -> Vec<(&'a str, Option<&'a str>)> {
    CHODL
        .iter()
        .map(|id| {
            (
                *id,
                set.iter().find(|(s, _)| s == id).map(|(_, level)| *level),
            )
        })
        .collect()
}

/// Runs `--pick --pick_order <order>` over the CHODL record with the levels
/// `set` names and returns the features left.
fn chodl_pick(dir: &Path, name: &str, set: &[(&str, &str)], order: &str) -> Vec<String> {
    let case = dir.join(name);
    fs::create_dir_all(&case).unwrap();
    let cache = cache_with_levels(&case, &chodl_levels(set));
    let text = run(
        &case,
        &cache,
        CHODL_RECORD,
        &["--pick", "--pick_order", order],
    );
    features(&text)
}

/// Under `--pick --pick_order tsl` the `tsl` category sorts the levels as
/// numbers (`<=>`): `300` loses to `2`, and to a transcript without a level,
/// whose default is 100, so the first of the eight at 100 in row order is the
/// pick; `00` is the number 0 and beats `1`.
#[test]
fn the_pick_order_ranks_a_level_as_a_number() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        chodl_pick(
            dir.path(),
            "300-vs-2",
            &[(CHODL[0], "300"), (CHODL[6], "2")],
            "tsl"
        ),
        [CHODL[6]]
    );
    assert_eq!(
        chodl_pick(dir.path(), "300-vs-none", &[(CHODL[0], "300")], "tsl"),
        [CHODL[1]]
    );
    assert_eq!(
        chodl_pick(
            dir.path(),
            "00-vs-1",
            &[(CHODL[0], "1"), (CHODL[6], "00")],
            "tsl"
        ),
        [CHODL[6]]
    );
}

/// A tie in the `tsl` category is the same text (`eq`): `007` sorts with the
/// `7`s but ends the tied prefix. With it second among them in row order the
/// first `7`, the shortest of the three, is the pick before `length` is read;
/// with it last, the two `7`s before it go on to `length`, which the longer
/// of the two wins over the longest of the three.
#[test]
fn the_pick_order_ties_a_level_as_text() {
    let dir = tempfile::tempdir().unwrap();
    // 338326 (711) `7`, 543733 (765) `007`, 1123871 (753) `7`.
    assert_eq!(
        chodl_pick(
            dir.path(),
            "007-second",
            &[(CHODL[1], "7"), (CHODL[6], "007"), (CHODL[7], "7")],
            "tsl,length",
        ),
        [CHODL[1]]
    );
    // 338326 (711) `7`, 543733 (765) `7`, 1131842 (822) `007`.
    assert_eq!(
        chodl_pick(
            dir.path(),
            "007-last",
            &[(CHODL[1], "7"), (CHODL[6], "7"), (CHODL[8], "007")],
            "tsl,length",
        ),
        [CHODL[6]]
    );
}
