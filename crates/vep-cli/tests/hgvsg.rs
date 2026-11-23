// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--hgvsg` through the `vep` binary: the refusal without a reference FASTA,
//! and the HGVSg strings Ensembl VEP 116.3 writes on the GRCh37 golden corpus's
//! reference, for records of the corpus and for multi-allelic records built on
//! it, in every output format; the name `--hgvsg_use_accession` prints for
//! three spellings of one chromosome and for an unplaced scaffold; and the
//! accession the `--ga4gh_vrs` object names beside it.

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
/// ALT) with the HGVSg Ensembl VEP 116.3 writes for each under `--hgvsg --fasta`
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
/// with the HGVSg Ensembl VEP 116.3 writes for each allele: `34775849 CAT C,CGT`
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
/// allele string and span before the rows are built (OutputFactory.pm:2300-2302)
/// and the notation is looked up under the allele as the record writes it, by
/// its allele number (:1386-1390). Ensembl VEP 116.3 prints the strings of the
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

/// The rows of the record `intron_snv_1` of the `GRCh38-identifiers` golden
/// corpus under `--hgvsg --hgvsg_use_accession`, `{name}` standing for the
/// chromosome name the input writes. Ensembl VEP 116.3 writes them offline
/// for `21` and `chr21`, which the synonyms table pairs with `NC_000021.9`,
/// and with a database connection (`--database --db_version 116`) for `Chr21`,
/// a spelling the cache resolves and the table does not hold.
const ACCESSION_ROWS: &[&str] = &[
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00000427445\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00000612746\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00000612957\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;FLAGS=cds_end_NF;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00000618007\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00000622113\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00001091160\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00001091161\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00001091162\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
    "intron_snv_1\t{name}:10538885\tC\tENSG00000274391\tENST00001131865\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C",
];

/// `--hgvsg_use_accession` names the accession of the chromosome a name
/// resolves to: `Chr21` prints the rows of Ensembl VEP 116.3's run with a
/// database connection, field for field, where its offline run ends on the
/// record (OutputFactory.pm 1368-1372 reads the synonyms of the name as written
/// alone), and `21` and `chr21` print the rows of its offline run.
#[test]
fn hgvsg_use_accession_names_the_resolved_chromosome() {
    let dir = corpus("GRCh38-identifiers");
    let out = tempfile::tempdir().unwrap();
    for name in ["Chr21", "chr21", "21"] {
        let input = out.path().join(format!("{name}.vcf"));
        std::fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 {name}\t10538885\tintron_snv_1\tA\tC\t.\t.\t.\n"
            ),
        )
        .unwrap();
        let output = out.path().join(format!("{name}.txt"));
        let result = common::isolated_vep_anywhere()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", output.to_str().unwrap()])
            .args(["--json_cache", dir.join("json_cache").to_str().unwrap()])
            .args(["--offline", "--assembly", "GRCh38"])
            .args(["--hgvsg", "--hgvsg_use_accession"])
            .args(["--fasta", dir.join("reference.fa.gz").to_str().unwrap()])
            .args(["--force_overwrite", "--no_stats", "--quiet"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let rows: Vec<String> = std::fs::read_to_string(&output)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_string)
            .collect();
        let expected: Vec<String> = ACCESSION_ROWS
            .iter()
            .map(|row| row.replace("{name}", name))
            .collect();
        assert_eq!(rows, expected, "{name}");
    }
}

/// The row of the record at `KI270728.1:17420` of the GRCh38-alt-contigs corpus,
/// under an identifier of its own, with `--hgvsg --hgvsg_use_accession`: the
/// one row Ensembl VEP 116.3 writes for it with a database connection under
/// each of the scaffold's three names. `{name}` stands for the name as the
/// input writes it.
const SCAFFOLD_ACCESSION_ROW: &str = "scaffold_snv\t{name}:17420\tC\tENSG00000274175\tENST00000616638\tTranscript\tmissense_variant\t187\t157\t53\tS/P\tTcg/Ccg\t-\tIMPACT=MODERATE;STRAND=1;HGVSg=NT_187383.1:g.17420T>C";

/// A sequence with no `NC_` synonym takes its RefSeq name: the unplaced
/// scaffold `KI270728.1` prints `NT_187383.1` under its Ensembl name, its UCSC
/// name and the RefSeq name itself, the slice's `RefSeq_genomic` synonym
/// Ensembl VEP 116.3 prints with a database connection (ensembl-variation
/// VariationFeature.pm 1912-1916); offline its lookup reads `NC_` names alone
/// (OutputFactory.pm 1369-1370). The synonyms file holds the three lines
/// Ensembl's table has for the scaffold, and the FASTA is a stand-in with the
/// record's reference base at its position, all a substitution's notation
/// reads.
#[test]
fn hgvsg_use_accession_names_a_scaffold_by_its_refseq_name() {
    let dir = corpus("GRCh38-alt-contigs");
    let out = tempfile::tempdir().unwrap();
    let synonyms = out.path().join("chr_synonyms.txt");
    std::fs::write(
        &synonyms,
        "NT_187383.1\tchr16_KI270728v1_random\n\
         KI270728.1\tchr16_KI270728v1_random\n\
         KI270728.1\tNT_187383.1\n",
    )
    .unwrap();
    let position = 17420;
    let sequence = format!("{}T{}", "N".repeat(position - 1), "N".repeat(60));
    let fasta = out.path().join("scaffold.fa");
    std::fs::write(&fasta, format!(">KI270728.1\n{sequence}\n")).unwrap();
    std::fs::write(
        out.path().join("scaffold.fa.fai"),
        format!(
            "KI270728.1\t{0}\t{1}\t{0}\t{2}\n",
            sequence.len(),
            ">KI270728.1\n".len(),
            sequence.len() + 1
        ),
    )
    .unwrap();
    for name in ["KI270728.1", "chr16_KI270728v1_random", "NT_187383.1"] {
        let input = out.path().join(format!("{name}.vcf"));
        std::fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 {name}\t{position}\tscaffold_snv\tT\tC\t.\t.\t.\n"
            ),
        )
        .unwrap();
        let output = out.path().join(format!("{name}.txt"));
        let result = common::isolated_vep_anywhere()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", output.to_str().unwrap()])
            .args(["--json_cache", dir.join("json_cache").to_str().unwrap()])
            .args(["--synonyms", synonyms.to_str().unwrap()])
            .args(["--offline", "--assembly", "GRCh38"])
            .args(["--hgvsg", "--hgvsg_use_accession"])
            .args(["--fasta", fasta.to_str().unwrap()])
            .args(["--force_overwrite", "--no_stats", "--quiet"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let rows: Vec<String> = std::fs::read_to_string(&output)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_string)
            .collect();
        assert_eq!(
            rows,
            [SCAFFOLD_ACCESSION_ROW.replace("{name}", name)],
            "{name}"
        );
    }
}

/// A sequence the synonyms table gives no RefSeq name is printed under
/// `--hgvsg_use_accession` by the cache's name for it, whatever the input
/// writes: the cache's name itself, a synonym the table pairs with it, or a
/// spelling the `chr` rule resolves. Ensembl VEP 116.3 prints the slice's own
/// name for such a sequence with a database connection (ensembl-variation
/// VariationFeature.pm 1916); offline it has no output for the record
/// (OutputFactory.pm 1368-1372). The table is built here: one line pairing
/// chromosome 21 with a name of no RefSeq form.
#[test]
fn hgvsg_use_accession_names_a_sequence_without_a_refseq_name_as_the_cache_does() {
    let out = tempfile::tempdir().unwrap();
    let synonyms = out.path().join("chr_synonyms.txt");
    std::fs::write(&synonyms, "21\tsequence_21\n").unwrap();
    let (record, expected) = RECORDS[0];
    let fields = record.strip_prefix("21\t").unwrap();
    for name in ["21", "sequence_21", "chr21"] {
        let input = out.path().join(format!("{name}.vcf"));
        std::fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 {name}\t{fields}\t.\t.\t.\n"
            ),
        )
        .unwrap();
        let output = out.path().join(format!("{name}.txt"));
        let result = common::isolated_vep_anywhere()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", output.to_str().unwrap()])
            .args([
                "--json_cache",
                corpus("GRCh37-hgvs").join("json_cache").to_str().unwrap(),
            ])
            .args(["--synonyms", synonyms.to_str().unwrap()])
            .args(["--offline", "--assembly", "GRCh37"])
            .args(["--hgvsg", "--hgvsg_use_accession"])
            .args(["--fasta", corpus_fasta().to_str().unwrap()])
            .args(["--force_overwrite", "--no_stats", "--quiet"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let written: Vec<String> = std::fs::read_to_string(&output)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .filter_map(|line| {
                line.split('\t')
                    .nth(13)?
                    .split(';')
                    .find_map(|pair| pair.strip_prefix("HGVSg="))
                    .map(str::to_string)
            })
            .collect();
        assert!(!written.is_empty(), "{name}: no row carries HGVSg");
        for value in &written {
            assert_eq!(value, expected, "{name}");
        }
    }
}

/// An `NC_` accession is taken wherever the table lists it among a name's
/// synonyms: after a RefSeq name of another prefix it is the one printed, as
/// Ensembl VEP's offline lookup keeps the `NC_` names alone (OutputFactory.pm
/// 1370). The table is built here: chromosome 21 paired first with a name of
/// another RefSeq prefix, then with an `NC_` name.
#[test]
fn hgvsg_use_accession_takes_an_nc_accession_listed_after_another_refseq_name() {
    let out = tempfile::tempdir().unwrap();
    let synonyms = out.path().join("chr_synonyms.txt");
    std::fs::write(&synonyms, "21\tNT_000000021.1\n21\tNC_000021.8\n").unwrap();
    let (record, notation) = RECORDS[0];
    let expected = format!("NC_000021.8{}", notation.strip_prefix("21").unwrap());
    let input = out.path().join("in.vcf");
    std::fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             {record}\t.\t.\t.\n"
        ),
    )
    .unwrap();
    let output = out.path().join("out.txt");
    let result = common::isolated_vep_anywhere()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args([
            "--json_cache",
            corpus("GRCh37-hgvs").join("json_cache").to_str().unwrap(),
        ])
        .args(["--synonyms", synonyms.to_str().unwrap()])
        .args(["--offline", "--assembly", "GRCh37"])
        .args(["--hgvsg", "--hgvsg_use_accession"])
        .args(["--fasta", corpus_fasta().to_str().unwrap()])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let written: Vec<String> = std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            line.split('\t')
                .nth(13)?
                .split(';')
                .find_map(|pair| pair.strip_prefix("HGVSg="))
                .map(str::to_string)
        })
        .collect();
    assert!(!written.is_empty(), "no row carries HGVSg");
    for value in &written {
        assert_eq!(value, &expected);
    }
}

/// Where the table lists one accession with and without its version, the
/// versioned form is the one printed, whichever the table lists first and
/// whichever of the three names the input writes: the sequence's name, the
/// accession without its version, the accession with it. The tables are built
/// here, one for an accession of another RefSeq prefix and one for an `NC_`
/// accession, each listing the form without a version first.
#[test]
fn hgvsg_use_accession_takes_the_versioned_form_of_an_accession_listed_twice() {
    let out = tempfile::tempdir().unwrap();
    let (record, notation) = RECORDS[0];
    let fields = record.strip_prefix("21\t").unwrap();
    for (bare, versioned) in [
        ("NT_000000021", "NT_000000021.4"),
        ("NC_000021", "NC_000021.8"),
    ] {
        let synonyms = out.path().join(format!("{bare}.synonyms.txt"));
        std::fs::write(
            &synonyms,
            format!("21\t{bare}\n21\t{versioned}\n{bare}\t{versioned}\n"),
        )
        .unwrap();
        let expected = format!("{versioned}{}", notation.strip_prefix("21").unwrap());
        for name in ["21", bare, versioned] {
            let input = out.path().join(format!("{bare}.{name}.vcf"));
            std::fs::write(
                &input,
                format!(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     {name}\t{fields}\t.\t.\t.\n"
                ),
            )
            .unwrap();
            let output = out.path().join(format!("{bare}.{name}.txt"));
            let result = common::isolated_vep_anywhere()
                .args(["-i", input.to_str().unwrap()])
                .args(["-o", output.to_str().unwrap()])
                .args([
                    "--json_cache",
                    corpus("GRCh37-hgvs").join("json_cache").to_str().unwrap(),
                ])
                .args(["--synonyms", synonyms.to_str().unwrap()])
                .args(["--offline", "--assembly", "GRCh37"])
                .args(["--hgvsg", "--hgvsg_use_accession"])
                .args(["--fasta", corpus_fasta().to_str().unwrap()])
                .args(["--force_overwrite", "--no_stats", "--quiet"])
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{bare} {name}: {}",
                String::from_utf8_lossy(&result.stderr)
            );
            let written: Vec<String> = std::fs::read_to_string(&output)
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('#'))
                .filter_map(|line| {
                    line.split('\t')
                        .nth(13)?
                        .split(';')
                        .find_map(|pair| pair.strip_prefix("HGVSg="))
                        .map(str::to_string)
                })
                .collect();
            assert!(!written.is_empty(), "{bare} {name}: no row carries HGVSg");
            for value in &written {
                assert_eq!(value, &expected, "{bare} {name}");
            }
        }
    }
}

/// Where the table lists two versions of one accession for a sequence, the
/// higher is the one printed, whichever the table lists first and whichever of
/// the three names the input writes. A cache's table links the previous
/// assembly's version of a sequence to it, and Ensembl VEP 116.3 with a
/// database connection prints the version of the assembly the run reads, the
/// higher (the slice's own `RefSeq_genomic` synonym, ensembl-variation
/// VariationFeature.pm 1912-1916). The tables are built here, one for an
/// accession of another RefSeq prefix and one for an `NC_` accession, each
/// listing the lower version first.
#[test]
fn hgvsg_use_accession_takes_the_highest_version_of_an_accession_listed_twice() {
    let out = tempfile::tempdir().unwrap();
    let (record, notation) = RECORDS[0];
    let fields = record.strip_prefix("21\t").unwrap();
    for (lower, higher) in [
        ("NT_000000021.2", "NT_000000021.3"),
        ("NC_000021.7", "NC_000021.8"),
    ] {
        let synonyms = out.path().join(format!("{lower}.synonyms.txt"));
        std::fs::write(
            &synonyms,
            format!("21\t{lower}\n21\t{higher}\n{lower}\t{higher}\n"),
        )
        .unwrap();
        let expected = format!("{higher}{}", notation.strip_prefix("21").unwrap());
        for name in ["21", lower, higher] {
            let input = out.path().join(format!("{lower}.{name}.vcf"));
            std::fs::write(
                &input,
                format!(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     {name}\t{fields}\t.\t.\t.\n"
                ),
            )
            .unwrap();
            let output = out.path().join(format!("{lower}.{name}.txt"));
            let result = common::isolated_vep_anywhere()
                .args(["-i", input.to_str().unwrap()])
                .args(["-o", output.to_str().unwrap()])
                .args([
                    "--json_cache",
                    corpus("GRCh37-hgvs").join("json_cache").to_str().unwrap(),
                ])
                .args(["--synonyms", synonyms.to_str().unwrap()])
                .args(["--offline", "--assembly", "GRCh37"])
                .args(["--hgvsg", "--hgvsg_use_accession"])
                .args(["--fasta", corpus_fasta().to_str().unwrap()])
                .args(["--force_overwrite", "--no_stats", "--quiet"])
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{lower} {name}: {}",
                String::from_utf8_lossy(&result.stderr)
            );
            let written: Vec<String> = std::fs::read_to_string(&output)
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('#'))
                .filter_map(|line| {
                    line.split('\t')
                        .nth(13)?
                        .split(';')
                        .find_map(|pair| pair.strip_prefix("HGVSg="))
                        .map(str::to_string)
                })
                .collect();
            assert!(!written.is_empty(), "{lower} {name}: no row carries HGVSg");
            for value in &written {
                assert_eq!(value, &expected, "{lower} {name}");
            }
        }
    }
}

/// The `--ga4gh_vrs` object names the sequence by the accession the row's
/// HGVSg carries under `--hgvsg_use_accession`: where the table lists the
/// chromosome's `NC_` accession in more than one form, both take the latest,
/// the versioned form over the one without a version and the higher of two
/// versions. The tables are built here, each listing the earlier form first.
#[test]
fn ga4gh_vrs_object_names_the_accession_the_hgvsg_names() {
    let out = tempfile::tempdir().unwrap();
    let (record, _) = RECORDS[0];
    for (earlier, later) in [("NC_000021", "NC_000021.8"), ("NC_000021.7", "NC_000021.8")] {
        let synonyms = out.path().join(format!("{earlier}.synonyms.txt"));
        std::fs::write(
            &synonyms,
            format!("21\t{earlier}\n21\t{later}\n{earlier}\t{later}\n"),
        )
        .unwrap();
        let input = out.path().join(format!("{earlier}.vcf"));
        std::fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 {record}\t.\t.\t.\n"
            ),
        )
        .unwrap();
        let output = out.path().join(format!("{earlier}.txt"));
        let result = common::isolated_vep_anywhere()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", output.to_str().unwrap()])
            .args([
                "--json_cache",
                corpus("GRCh37-hgvs").join("json_cache").to_str().unwrap(),
            ])
            .args(["--synonyms", synonyms.to_str().unwrap()])
            .args(["--offline", "--assembly", "GRCh37"])
            .args(["--hgvsg", "--hgvsg_use_accession", "--ga4gh_vrs"])
            .args(["--fasta", corpus_fasta().to_str().unwrap()])
            .args(["--force_overwrite", "--no_stats", "--quiet"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{earlier}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let text = std::fs::read_to_string(&output).unwrap();
        let extras: Vec<&str> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .filter_map(|line| line.split('\t').nth(13))
            .collect();
        assert!(!extras.is_empty(), "{earlier}: no row");
        for extra in extras {
            let named = |key: &str, end: &[char]| {
                let (_, rest) = extra.split_once(key)?;
                Some(rest.split(end).next().unwrap_or(rest))
            };
            assert_eq!(named("HGVSg=", &[':']), Some(later), "{earlier}: {extra}");
            assert_eq!(
                named("refseq:", &['"', '\\']),
                Some(later),
                "{earlier}: {extra}"
            );
        }
    }
}
