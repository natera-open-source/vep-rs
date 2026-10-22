// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--hgvsg` through the `vep` binary: the refusal without a reference FASTA,
//! and the HGVSg strings Ensembl VEP 116.2 writes on the GRCh37 golden corpus's
//! reference, for records of the corpus and for multi-allelic records built on
//! it, in every output format.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

mod common;

fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// `--hgvsg` with no FASTA, named or auto-detected, exits 1 with Ensembl VEP's
/// message (Runner.pm 759) under the `fasta` component, once the input is open
/// (Runner.pm 136 before 138) and before any record is read or the output is
/// created.
#[test]
fn hgvsg_without_fasta_exits_with_vep_message() {
    let dir = corpus("GRCh38-release116");
    let out = tempfile::tempdir().unwrap();
    let out_file = out.path().join("out.txt");
    let output = common::isolated_vep_anywhere()
        .args(["--hgvsg", "--json_cache"])
        .arg(dir.join("json_cache"))
        .arg("-i")
        .arg(dir.join("variants.vcf"))
        .args(["--quiet", "-o", out_file.to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)\n"
    );
    assert!(!out_file.exists());
}

/// Records of `tests/golden/116/GRCh37-hgvs/variants.vcf` (CHROM, POS, ID, REF,
/// ALT) with the HGVSg Ensembl VEP 116.2 writes for each under `--hgvsg --fasta`
/// on the corpus's chromosome 21 reference: a substitution; deletions where the
/// record sits (`190797`, `2878128`) and moved 3' along a repeat (`1299494`,
/// `2505239`); an insertion where the record sits (`3653535`), one of 68 bases
/// moved two bases 3' and rotated (`166800`), two that become duplications
/// after the move (`1117421`, `1209862`) and two whose move runs past the
/// release/115 bound of the flank less the allele's length (`2074993`, 21 bases
/// moved 92; `571387`, 420 bases moved 421); and two deletion-insertions
/// (`1325111`, `3906334`).
const RECORDS: &[(&str, &str)] = &[
    ("21\t34775850\t1066015\tA\tC", "21:g.34775850A>C"),
    ("21\t35743145\t190797\tCCT\tC", "21:g.35743146_35743147del"),
    ("21\t34727851\t1299494\tTATGACC\tT", "21:g.34727853_34727858del"),
    ("21\t34619194\t2878128\tCA\tC", "21:g.34619195del"),
    ("21\t34929621\t2505239\tAG\tA", "21:g.34929623del"),
    ("21\t36259299\t3653535\tG\tGT", "21:g.36259299_36259300insT"),
    (
        "21\t44483171\t166800\tG\tGGGTTCCACCTCGTAGGTTGTCTGCTCCGTCTGGTTCAGCTCCTCCGGCTCTGCGAGGATGGACCCTTC",
        "21:g.44483173_44483174insTTCCACCTCGTAGGTTGTCTGCTCCGTCTGGTTCAGCTCCTCCGGCTCTGCGAGGATGGACCCTTCGG",
    ),
    ("21\t46925297\t1117421\tC\tCCCCCCTGGG", "21:g.46925315_46925323dup"),
    ("21\t35190736\t1209862\tT\tTA", "21:g.35190737dup"),
    (
        "21\t34927396\t2074993\tA\tAAGCCGCCGCAGCCGCACCCCC",
        "21:g.34927468_34927488dup",
    ),
    (
        "21\t44485538\t571387\tG\tGCCAGGCCACCCCCACGTGTGACTCCGGGGAGTCGAACCTGGCATTGGTGGGCGTCCTCACAATCTCAGCCCCCAGTGCCCGCAGCACGTCCACCTGCAGGAGGGAAAGCGGTGGCCTGCACCTTCCGCCTGGCCCAGGCACCCTCATCCCCTGCCCTATGACCCCGCCCCTGGCCACGCCCACCCACCTTCTCGGAGCTCATCTTCTCTGGCATCACGATGATGCAGCGATAGCCCCTCACTGCCGCAGCCAGGGCCAGCCCGATCCCTGAGGGCACACAGAGGGTGAGAGGGGCCCAGTGACCCCCCAAGCCCTGCCCCGCCCCTGCCTGGGACACAGGGGCACACCCCGATGCCGGTTCCCTTAGGGCCCAGGGAAGAGGGTTCTGTGGGATTCCAAAATTGCCCAACATGGCTGCTT",
        "21:g.44485540_44485959dup",
    ),
    (
        "21\t34924514\t1325111\tTTAG\tAATTTTATTCA",
        "21:g.34924514_34924517delinsAATTTTATTCA",
    ),
    (
        "21\t34927468\t3906334\tC\tGAGCCGCACCCCCAGCCGCCGG",
        "21:g.34927468delinsGAGCCGCACCCCCAGCCGCCGG",
    ),
];

/// Multi-allelic records on the corpus reference (none sits in the corpus itself)
/// with the HGVSg Ensembl VEP 116.2 writes for each allele: `34775849 CAT C,CGT`
/// (alleles `AT/-/GT` at 34775850-34775851 once the shared `C` is chopped) and
/// `34775850 AT AG,GT` (no indel, so nothing is chopped). An allele that trims to
/// a substitution is written as a deletion-insertion of the bases the record
/// writes, at the trimmed position.
const MULTI_ALLELIC_RECORDS: &[(&str, &[(&str, &str)])] = &[
    (
        "21\t34775849\tma1\tCAT\tC,CGT",
        &[
            ("-", "21:g.34775850_34775851del"),
            ("GT", "21:g.34775850delinsGT"),
        ],
    ),
    (
        "21\t34775850\tma2\tAT\tAG,GT",
        &[
            ("AG", "21:g.34775851delinsAG"),
            ("GT", "21:g.34775850delinsGT"),
        ],
    ),
];

/// The corpus reference, bgzipped beside its `.fai` and `.gzi`, as `--fasta`
/// reads it.
fn corpus_fasta() -> PathBuf {
    corpus("GRCh37-hgvs").join("reference.fa.gz")
}

/// Every record line of the two tables.
fn all_records() -> impl Iterator<Item = &'static str> {
    RECORDS
        .iter()
        .map(|(record, _)| *record)
        .chain(MULTI_ALLELIC_RECORDS.iter().map(|(record, _)| *record))
}

/// The HGVSg Ensembl VEP writes for the row of record `id` whose Allele column is
/// `allele`: one value per bi-allelic record, one per allele of a multi-allelic one.
fn expected_hgvsg(id: &str, allele: &str) -> &'static str {
    if let Some((_, expected)) = RECORDS.iter().find(|(record, _)| record_id(record) == id) {
        return expected;
    }
    MULTI_ALLELIC_RECORDS
        .iter()
        .find(|(record, _)| record_id(record) == id)
        .and_then(|(_, alleles)| alleles.iter().find(|(a, _)| *a == allele))
        .map(|(_, expected)| *expected)
        .unwrap_or_else(|| panic!("no expectation for {id} allele {allele}"))
}

/// Runs the binary on every record of the two tables in `format` with `extra`
/// flags and returns the output's lines, the `##` header lines included.
fn annotate(dir: &Path, fasta: &Path, format: &str, extra: &[&str]) -> Vec<String> {
    let input = dir.join("in.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in all_records() {
        text.push_str(record);
        text.push_str("\t.\t.\t.\n");
    }
    std::fs::write(&input, text).unwrap();
    let output = dir.join(format!("out.{format}"));
    let result = common::isolated_vep_anywhere()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args([
            "--json_cache",
            corpus("GRCh37-hgvs").join("json_cache").to_str().unwrap(),
        ])
        .args(["--offline", "--assembly", "GRCh37", "--hgvsg"])
        .args(["--fasta", fasta.to_str().unwrap()])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(extra)
        .args(if format == "default" {
            vec![]
        } else {
            vec![format!("--{format}")]
        })
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "--{format}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// The ID column of a record line.
fn record_id(record: &str) -> &str {
    record.split('\t').nth(2).unwrap()
}

/// Every row of every record carries Ensembl VEP's HGVSg for its allele: the
/// `HGVSg=` key of the default format's Extra column, the `HGVSg` column of the
/// tab format, the `HGVSg` field of every VCF CSQ entry, and the `hgvsg` key of
/// every consequence object of the JSON format.
#[test]
fn hgvsg_matches_vep_on_grch37_corpus_records_in_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = corpus_fasta();

    for format in ["default", "tab"] {
        let rows = annotate(dir.path(), &fasta, format, &[]);
        let tab_column = (format == "tab").then(|| {
            let header = rows
                .iter()
                .find(|r| r.starts_with('#') && !r.starts_with("##"))
                .unwrap();
            header.split('\t').position(|c| c == "HGVSg").unwrap()
        });
        for record in all_records() {
            let id = record_id(record);
            let row_prefix = format!("{id}\t");
            let mut seen = Vec::new();
            for row in rows.iter().filter(|r| r.starts_with(&row_prefix)) {
                let columns: Vec<&str> = row.split('\t').collect();
                let value = match tab_column {
                    Some(index) => columns[index],
                    None => columns[13]
                        .split(';')
                        .find_map(|kv| kv.strip_prefix("HGVSg="))
                        .unwrap_or_default(),
                };
                assert_eq!(
                    value,
                    expected_hgvsg(id, columns[2]),
                    "--{format} {id}: {row}"
                );
                seen.push(columns[2]);
            }
            assert_every_allele_seen(id, &seen, format);
        }
    }

    let vcf_rows = annotate(dir.path(), &fasta, "vcf", &[]);
    let csq_format = vcf_rows
        .iter()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .map(|fields| fields.trim_end_matches("\">"))
        .unwrap_or_default();
    let csq_field = |name: &str| {
        csq_format
            .split('|')
            .position(|f| f == name)
            .unwrap_or_else(|| panic!("CSQ format lacks {name}: {csq_format}"))
    };
    let (allele_field, hgvsg_field) = (csq_field("Allele"), csq_field("HGVSg"));
    for record in all_records() {
        let id = record_id(record);
        let line = vcf_rows
            .iter()
            .find(|r| !r.starts_with('#') && r.split('\t').nth(2) == Some(id))
            .unwrap_or_else(|| panic!("--vcf: no line for {id}"));
        let csq = line
            .split('\t')
            .nth(7)
            .unwrap()
            .split(';')
            .find_map(|kv| kv.strip_prefix("CSQ="))
            .unwrap_or_else(|| panic!("--vcf {id}: no CSQ"));
        let mut seen = Vec::new();
        for entry in csq.split(',') {
            let fields: Vec<&str> = entry.split('|').collect();
            assert_eq!(
                fields[hgvsg_field],
                expected_hgvsg(id, fields[allele_field]),
                "--vcf {id}: {entry}"
            );
            seen.push(fields[allele_field]);
        }
        assert_every_allele_seen(id, &seen, "vcf");
    }

    let json_rows = annotate(dir.path(), &fasta, "json", &[]);
    for record in all_records() {
        let id = record_id(record);
        let object: serde_json::Value = json_rows
            .iter()
            .map(|r| serde_json::from_str(r).unwrap())
            .find(|o: &serde_json::Value| o["id"] == *id)
            .unwrap_or_else(|| panic!("--json: no object for {id}"));
        let mut seen = Vec::new();
        for key in ["transcript_consequences", "intergenic_consequences"] {
            for c in object[key].as_array().into_iter().flatten() {
                let allele = c["variant_allele"].as_str().unwrap_or_default();
                assert_eq!(c["hgvsg"], expected_hgvsg(id, allele), "--json {id}: {c}");
                seen.push(allele);
            }
        }
        assert_every_allele_seen(id, &seen, "json");
    }
}

/// Record `id` had a row for every allele the tables expect: each listed allele
/// of a multi-allelic record, at least one row of a bi-allelic one (a dropped
/// allele would leave its expectation unconsulted).
fn assert_every_allele_seen(id: &str, seen: &[&str], format: &str) {
    match MULTI_ALLELIC_RECORDS
        .iter()
        .find(|(record, _)| record_id(record) == id)
    {
        Some((_, alleles)) => {
            for &(allele, _) in alleles.iter() {
                assert!(
                    seen.contains(&allele),
                    "--{format} {id}: no row for allele {allele}"
                );
            }
        }
        None => assert!(!seen.is_empty(), "--{format}: no row for {id}"),
    }
}

/// Under `--minimal` each allele of a multi-allelic record is annotated on its
/// trimmed span, but its HGVSg is the record's: the rejoin restores the record's
/// allele string and span before the rows are built (OutputFactory.pm:2291-2293)
/// and the notation is looked up under the allele as the record writes it, by
/// its allele number (:1386-1390). Ensembl VEP 116.2 prints the strings of the
/// table above with or without the flag: `21:g.34775850_34775851del` and
/// `21:g.34775850delinsGT` for `ma1`, whose rows print the trimmed alleles `-`
/// and `G`; `21:g.34775851delinsAG` and `21:g.34775850delinsGT` for `ma2`,
/// whose two trimmed alleles both print `G`. Every row carries its allele's
/// notation in the default, VCF and JSON formats, keyed by `ALLELE_NUM`.
#[test]
fn hgvsg_under_minimal_is_the_record_notation_of_each_written_allele() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = corpus_fasta();
    let flags = ["--minimal", "--allele_number"];
    let expected_by_number = |id: &str, number: usize| -> &'static str {
        MULTI_ALLELIC_RECORDS
            .iter()
            .find(|(record, _)| record_id(record) == id)
            .map(|(_, alleles)| alleles[number - 1].1)
            .unwrap_or_else(|| panic!("no multi-allelic record {id}"))
    };
    let multi_allelic_ids: Vec<&str> = MULTI_ALLELIC_RECORDS
        .iter()
        .map(|(record, _)| record_id(record))
        .collect();

    let rows = annotate(dir.path(), &fasta, "default", &flags);
    for id in &multi_allelic_ids {
        let row_prefix = format!("{id}\t");
        let mut numbers_seen = BTreeSet::new();
        for row in rows.iter().filter(|r| r.starts_with(&row_prefix)) {
            let columns: Vec<&str> = row.split('\t').collect();
            let extra = |key: &str| {
                columns[13]
                    .split(';')
                    .find_map(|kv| kv.strip_prefix(key))
                    .unwrap_or_else(|| panic!("--default {id}: no {key} on {row}"))
            };
            let number: usize = extra("ALLELE_NUM=").parse().unwrap();
            assert_eq!(extra("HGVSg="), expected_by_number(id, number), "{row}");
            assert_eq!(extra("MINIMISED="), "1", "{row}");
            numbers_seen.insert(number);
        }
        assert_eq!(
            numbers_seen,
            BTreeSet::from([1, 2]),
            "--default {id}: a row per allele"
        );
    }

    let vcf_rows = annotate(dir.path(), &fasta, "vcf", &flags);
    let csq_format = vcf_rows
        .iter()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .map(|fields| fields.trim_end_matches("\">"))
        .unwrap_or_default();
    let csq_field = |name: &str| csq_format.split('|').position(|f| f == name).unwrap();
    let (number_field, hgvsg_field) = (csq_field("ALLELE_NUM"), csq_field("HGVSg"));
    for id in &multi_allelic_ids {
        let line = vcf_rows
            .iter()
            .find(|r| !r.starts_with('#') && r.split('\t').nth(2) == Some(id))
            .unwrap_or_else(|| panic!("--vcf: no line for {id}"));
        let csq = line
            .split('\t')
            .nth(7)
            .unwrap()
            .split(';')
            .find_map(|kv| kv.strip_prefix("CSQ="))
            .unwrap();
        let mut numbers_seen = BTreeSet::new();
        for entry in csq.split(',') {
            let fields: Vec<&str> = entry.split('|').collect();
            let number: usize = fields[number_field].parse().unwrap();
            assert_eq!(
                fields[hgvsg_field],
                expected_by_number(id, number),
                "{entry}"
            );
            numbers_seen.insert(number);
        }
        assert_eq!(
            numbers_seen,
            BTreeSet::from([1, 2]),
            "--vcf {id}: an entry per allele"
        );
    }

    let json_rows = annotate(dir.path(), &fasta, "json", &flags);
    for id in &multi_allelic_ids {
        let object: serde_json::Value = json_rows
            .iter()
            .map(|r| serde_json::from_str(r).unwrap())
            .find(|o: &serde_json::Value| o["id"] == *id)
            .unwrap_or_else(|| panic!("--json: no object for {id}"));
        let mut numbers_seen = BTreeSet::new();
        for key in ["transcript_consequences", "intergenic_consequences"] {
            for c in object[key].as_array().into_iter().flatten() {
                let number = c["allele_num"].as_u64().unwrap() as usize;
                assert_eq!(
                    c["hgvsg"],
                    expected_by_number(id, number),
                    "--json {id}: {c}"
                );
                numbers_seen.insert(number);
            }
        }
        assert_eq!(
            numbers_seen,
            BTreeSet::from([1, 2]),
            "--json {id}: a consequence per allele"
        );
    }
}

/// The 3' walk of `get_3prime_seq_offset` (Utils/Sequence.pm 802) over `post`,
/// stopping after `bound` comparisons: the whole flank in release/116, the flank
/// less the allele's length in release/115.
fn three_prime_offset(allele: &[u8], post: &[u8], bound: usize) -> usize {
    let mut seq = allele.to_vec();
    let mut offset = 0;
    for n in 0..=bound {
        match post.get(n) {
            Some(&next) if seq.first() == Some(&next) => {
                offset += 1;
                seq.rotate_left(1);
            }
            _ => break,
        }
    }
    offset
}

/// The corpus carries insertion records whose 3' walk reaches past the
/// release/115 bound: `2074993` (a 21-base insertion after 34927396, walk 80 under
/// the old bound and 92 over the whole flank) and `571387` (a 420-base insertion
/// after 44485538, 51 and 421), so the `hgvsg` flag set's `covers` entry for the
/// whole-flank rule is true of the input. Every other indel of the corpus walks
/// the same distance under both bounds.
#[test]
fn grch37_hgvs_corpus_carries_a_full_flank_shift_record() {
    let fasta = vep_fasta::IndexedFasta::from_path(corpus_fasta()).unwrap();
    let mut differing = Vec::new();
    for line in std::fs::read_to_string(corpus("GRCh37-hgvs").join("variants.vcf"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let columns: Vec<&str> = line.split('\t').collect();
        let (pos, id, reference, alt) = (
            columns[1].parse::<u64>().unwrap(),
            columns[2],
            columns[3].as_bytes(),
            columns[4].as_bytes(),
        );
        if alt.contains(&b',') || reference.len() == alt.len() {
            continue;
        }
        // Ensembl VEP's minimisation (`trim_sequences`): the shared bases leave
        // the front, then the back; a pure insertion or deletion remains, and the
        // walk starts at the first base after it.
        let shared = reference
            .iter()
            .zip(alt)
            .take_while(|(r, a)| r == a)
            .count();
        let (mut deleted, mut inserted) = (&reference[shared..], &alt[shared..]);
        while let (Some(r), Some(a)) = (deleted.last(), inserted.last()) {
            if r != a {
                break;
            }
            deleted = &deleted[..deleted.len() - 1];
            inserted = &inserted[..inserted.len() - 1];
        }
        let (allele, after) = match (deleted.is_empty(), inserted.is_empty()) {
            (false, true) => (deleted, pos + shared as u64 + deleted.len() as u64),
            (true, false) => (inserted, pos + shared as u64),
            _ => continue,
        };
        let flank = if allele.len() > 50 {
            allele.len() + 50
        } else {
            100
        };
        let post = fasta
            .sequence("21", after, after + flank as u64 - 1)
            .unwrap();
        let whole = three_prime_offset(allele, &post, post.len());
        let old = three_prime_offset(allele, &post, post.len() - allele.len());
        if whole != old {
            differing.push((id.to_string(), old, whole));
        }
    }
    assert_eq!(
        differing,
        [
            ("2074993".to_string(), 80, 92),
            ("571387".to_string(), 51, 421)
        ]
    );
}
