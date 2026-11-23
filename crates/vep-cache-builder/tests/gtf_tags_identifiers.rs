// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A GTF given beside a GFF3 for its tags gives each transcript the tags of
//! the `transcript` line that names it whole: transcripts whose identifiers
//! share a first component (`T1.1`, `T1.2`) each take their own line's, and
//! a line that carries a version, as `transcript_version` or inside its
//! `transcript_id`, tags the model of that version alone. A line states its
//! identifier and its version apart, so a version written both ways names
//! the line once, and a model at version 1 (`T1`) and a sibling written
//! `T1.1` each take their own line.

use std::collections::BTreeMap;
use std::process::Command;

use serde_json::Value;

/// A GFF3 gene `G1` with one coding transcript per identifier of `ids`, each
/// transcript line carrying `version` when one is named.
fn model(ids: &[&str], version: Option<u32>) -> String {
    let version_attribute = version.map_or(String::new(), |v| format!(";version={v}"));
    let mut text = String::from(
        "##gff-version 3\n\
         21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n",
    );
    for id in ids {
        text += &format!(
            "21\tx\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:{id};Parent=gene:G1;\
             biotype=protein_coding;transcript_id={id}{version_attribute}\n\
             21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:{id}\n\
             21\tx\tCDS\t1100\t1108\t.\t+\t0\tParent=transcript:{id};protein_id=P-{id}\n"
        );
    }
    text
}

/// A GTF `transcript` line of gene `G1` for `id`, with `transcript_version`
/// when one is named and one `tag` per entry of `tags`.
fn gtf_line(id: &str, version: Option<u32>, tags: &[&str]) -> String {
    let version = version.map_or(String::new(), |v| format!(" transcript_version \"{v}\";"));
    let tags: String = tags.iter().map(|t| format!(" tag \"{t}\";")).collect();
    format!(
        "21\tx\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"{id}\";\
         {version} transcript_biotype \"protein_coding\";{tags}\n"
    )
}

/// Builds `model` with `gtf` read for its tags and returns the attribute
/// codes of every transcript of the cache, by stable identifier.
fn attribute_codes(model: &str, gtf: &str) -> BTreeMap<String, Vec<String>> {
    let dir = tempfile::tempdir().unwrap();
    let model_path = dir.path().join("model.gff3");
    std::fs::write(&model_path, model).unwrap();
    let gtf_path = dir.path().join("tags.gtf");
    std::fs::write(&gtf_path, gtf).unwrap();
    let peptides = dir.path().join("peptides.fa");
    std::fs::write(&peptides, ">ENSP1 pep transcript:none\nMK\n").unwrap();
    let out = dir.path().join("cache");
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh38", "--gff3"])
        .arg(&model_path)
        .arg("--gtf")
        .arg(&gtf_path)
        .arg("--fasta")
        .arg(&peptides)
        .arg("--output-dir")
        .arg(&out)
        .output()
        .expect("vep-cache-builder runs");
    let log = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{log}");

    let mut codes = BTreeMap::new();
    for shard in std::fs::read_dir(out.join("transcripts/21")).unwrap() {
        let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
        let records: Vec<Value> = serde_json::from_str(&text).unwrap();
        for record in records {
            let of_record = record["attributes"]
                .as_array()
                .into_iter()
                .flatten()
                .map(|a| a["code"].as_str().unwrap().to_string())
                .collect();
            codes.insert(record["stable_id"].as_str().unwrap().to_string(), of_record);
        }
    }
    codes
}

/// The one transcript `id` of a build with the attribute codes `codes`.
fn one(id: &str, codes: &[&str]) -> BTreeMap<String, Vec<String>> {
    BTreeMap::from([(
        id.to_string(),
        codes.iter().map(|c| c.to_string()).collect(),
    )])
}

/// Transcripts whose identifiers share their first component each carry the
/// kept tags of their own GTF line: the `cds_start_NF` of `T1.1` is not
/// `T1.2`'s, and three siblings with three different tags keep one each.
#[test]
fn sibling_transcripts_keep_the_tags_of_their_own_gtf_lines() {
    let two = attribute_codes(
        &model(&["T1.1", "T1.2"], None),
        &[
            gtf_line("T1.1", None, &["cds_start_NF"]),
            gtf_line("T1.2", None, &[]),
        ]
        .concat(),
    );
    assert_eq!(two["T1.1"], ["cds_start_NF"]);
    assert!(two["T1.2"].is_empty(), "{two:?}");

    let three = attribute_codes(
        &model(&["T1.1", "T1.2", "T1.3"], None),
        &[
            gtf_line("T1.1", None, &["cds_start_NF"]),
            gtf_line("T1.2", None, &["cds_end_NF"]),
            gtf_line("T1.3", None, &["basic"]),
        ]
        .concat(),
    );
    assert_eq!(three["T1.1"], ["cds_start_NF"]);
    assert_eq!(three["T1.2"], ["cds_end_NF"]);
    assert_eq!(three["T1.3"], ["gencode_basic"]);
}

/// The transcript `T1` at version 3 takes the tags of the GTF line of that
/// version, written as `transcript_version` or inside the `transcript_id`
/// (`T1.3`), and of a line that carries no version; a line of version 4,
/// written either way, is another version's and tags nothing. A model that
/// states no version takes no line that carries one, and a model that writes
/// the version inside its identifier takes the line that states it apart.
#[test]
fn a_gtf_line_tags_the_model_of_its_version() {
    let versioned = model(&["T1"], Some(3));
    let tagged = one("T1", &["cds_start_NF"]);
    let untagged = one("T1", &[]);
    let tag = ["cds_start_NF"];

    assert_eq!(
        attribute_codes(&versioned, &gtf_line("T1", Some(3), &tag)),
        tagged
    );
    assert_eq!(
        attribute_codes(&versioned, &gtf_line("T1.3", None, &tag)),
        tagged
    );
    assert_eq!(
        attribute_codes(&versioned, &gtf_line("T1", None, &tag)),
        tagged
    );
    assert_eq!(
        attribute_codes(&versioned, &gtf_line("T1", Some(4), &tag)),
        untagged
    );
    assert_eq!(
        attribute_codes(&versioned, &gtf_line("T1.4", None, &tag)),
        untagged
    );

    assert_eq!(
        attribute_codes(&model(&["T1"], None), &gtf_line("T1", Some(3), &tag)),
        untagged
    );
    assert_eq!(
        attribute_codes(&model(&["T1.3"], None), &gtf_line("T1", Some(3), &tag)),
        one("T1.3", &["cds_start_NF"])
    );
}

/// A GTF line that writes its version twice, inside its `transcript_id` and
/// as `transcript_version` (`T1.3` at version 3), is the line of `T1` at
/// version 3: it tags that model, the model written `T1.3` at that version or
/// at none, and no model of another version.
#[test]
fn a_version_written_twice_names_the_line_once() {
    let tag = ["cds_start_NF"];
    let twice = gtf_line("T1.3", Some(3), &tag);

    assert_eq!(
        attribute_codes(&model(&["T1"], Some(3)), &twice),
        one("T1", &tag)
    );
    assert_eq!(
        attribute_codes(&model(&["T1.3"], Some(3)), &twice),
        one("T1.3", &tag)
    );
    assert_eq!(
        attribute_codes(&model(&["T1.3"], None), &twice),
        one("T1.3", &tag)
    );

    assert_eq!(
        attribute_codes(&model(&["T1"], Some(4)), &twice),
        one("T1", &[])
    );
    assert_eq!(
        attribute_codes(&model(&["T1"], None), &twice),
        one("T1", &[])
    );
    assert_eq!(
        attribute_codes(&model(&["T1.3"], Some(4)), &twice),
        one("T1.3", &[])
    );
}

/// The transcript `T1` at version 1 and its sibling written `T1.1` each carry
/// the kept tags of their own GTF line, in either order of the lines: the
/// line that states `T1` at version 1 beside the line written `T1.1`, with a
/// version of its own or without; the line written `T1` beside the line
/// written `T1.1` in a file that writes no version; and the line written
/// `T1.1` beside the line written `T1.1.1` in a file that writes each version
/// inside the identifier. Where the file holds the line written `T1.1` alone,
/// it is the sibling's, whose identifier it is as written, and `T1` carries
/// no tag.
#[test]
fn a_versioned_model_and_a_dotted_sibling_keep_their_own_tags() {
    let siblings = model(&["T1", "T1.1"], Some(1));
    let first = ["cds_start_NF"];
    let second = ["cds_end_NF"];
    let own = BTreeMap::from([
        ("T1".to_string(), vec!["cds_start_NF".to_string()]),
        ("T1.1".to_string(), vec!["cds_end_NF".to_string()]),
    ]);

    for lines in [
        [
            gtf_line("T1", Some(1), &first),
            gtf_line("T1.1", None, &second),
        ],
        [
            gtf_line("T1.1", None, &second),
            gtf_line("T1", Some(1), &first),
        ],
        [
            gtf_line("T1", Some(1), &first),
            gtf_line("T1.1", Some(1), &second),
        ],
        [
            gtf_line("T1", None, &first),
            gtf_line("T1.1", None, &second),
        ],
        [
            gtf_line("T1.1", None, &first),
            gtf_line("T1.1.1", None, &second),
        ],
    ] {
        assert_eq!(
            attribute_codes(&siblings, &lines.concat()),
            own,
            "{lines:?}"
        );
    }

    let alone = attribute_codes(&siblings, &gtf_line("T1.1", None, &second));
    assert!(alone["T1"].is_empty(), "{alone:?}");
    assert_eq!(alone["T1.1"], ["cds_end_NF"]);
}
