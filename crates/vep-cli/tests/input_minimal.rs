// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary runs under `--minimal`: the flag switches on
//! `--allele_number`, every bi-allelic record is trimmed and marked, a
//! multi-allelic record whose alleles trim is split and rejoined as Ensembl VEP
//! 116.2 rejoins it (the record's Location, the trimmed Allele, the written
//! position as ALLELE_NUM, the rows in the rejoin's key order) and the three
//! position columns follow the output format. Every expected value is Ensembl
//! VEP 116.2's output on the record under the same flags, except the
//! `--coding_only` rows of a split record, whose documented difference the last
//! test pins.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The identifiers corpus: a chromosome 21 cache slice over coding exons of TPTE.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A record whose four ALTs of three shapes split under `--minimal`: chopped to
/// `AAG/GAG/AGG/-` at 21:10542415-10542417, its alleles trim to `A/G` at 10542415
/// (ALLELE_NUM 1, `missense_variant`), `A/G` at 10542416 (2, `synonymous_variant`)
/// and `AAG/-` (3, `inframe_deletion`) on ENST00000427445.
const MULTI_THREE: &str = "21\t10542414\tmulti_three\tAAAG\tAGAG,AAGG,A\t.\t.\t.\n";

/// A bi-allelic substitution trimmed on both sides to `AA/GG` at 10542415-10542416.
const MNP: &str = "21\t10542414\tmnp_two_sided\tAAAG\tAGGG\t.\t.\t.\n";

fn vep(dir: &Path, records: &str, format: &[&str], extra: &[&str]) -> Output {
    let input = dir.join("in.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join("out");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .args([
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(format)
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{extra:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

/// The data lines of the run's output.
fn data_lines(dir: &Path) -> Vec<String> {
    std::fs::read_to_string(dir.join("out"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The rows on one transcript as `(Location, Allele, cDNA_position, CDS_position,
/// Protein_position, Consequence, Extra)`.
fn rows_on(lines: &[String], transcript: &str) -> Vec<Vec<String>> {
    lines
        .iter()
        .filter(|line| line.contains(transcript))
        .map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            [c[1], c[2], c[7], c[8], c[9], c[6], c[13]]
                .iter()
                .map(|s| s.to_string())
                .collect()
        })
        .collect()
}

/// `--minimal` switches on `--allele_number` (the option set at Config.pm:524-529):
/// every row prints `ALLELE_NUM`, and a bi-allelic record `MINIMISED=1` whether a
/// base was trimmed (`AAAG/AGGG` to `AA/GG`) or not (an SNV).
#[test]
fn minimal_implies_allele_number() {
    let dir = tempfile::tempdir().unwrap();
    let snv = "21\t10542415\tsnv\tA\tG\t.\t.\t.\n";
    vep(
        dir.path(),
        &format!("{snv}{MNP}"),
        &[],
        &["--minimal", "--uploaded_allele"],
    );
    let lines = data_lines(dir.path());
    assert!(lines.len() > 2, "{lines:?}");
    for line in &lines {
        let extra = line.split('\t').nth(13).unwrap();
        assert!(extra.contains("ALLELE_NUM=1"), "{line}");
        assert!(extra.contains("MINIMISED=1"), "{line}");
    }
    let mnp = rows_on(&lines, "ENST00000427445");
    let mnp: Vec<&Vec<String>> = mnp
        .iter()
        .filter(|r| r[6].contains("UPLOADED_ALLELE=AAAG/AGGG"))
        .collect();
    assert_eq!(mnp.len(), 1, "{lines:?}");
    assert_eq!(
        &mnp[0][..6],
        [
            "21:10542415-10542416",
            "GG",
            "438-439",
            "86-87",
            "29",
            "missense_variant",
        ]
    );
    // Without the flag the record keeps its alleles and no ALLELE_NUM is printed.
    vep(dir.path(), MNP, &[], &["--uploaded_allele"]);
    let plain = rows_on(&data_lines(dir.path()), "ENST00000427445");
    assert_eq!(plain[0][0], "21:10542414-10542417");
    assert_eq!(plain[0][1], "AGGG");
    assert!(!plain[0][6].contains("ALLELE_NUM"), "{}", plain[0][6]);
}

/// The rejoin of a split record (OutputFactory.pm:2272-2354): every row prints the
/// record's Location `21:10542415-10542417`, the trimmed Allele and the written
/// position as ALLELE_NUM, the rows sorted by the keys the merged alleles are stored
/// under (`A/G_<id>` and `AAG/-_<id>` before the first allele's bare `<id>`), each
/// `MINIMISED=1` with the chopped string as `UPLOADED_ALLELE`, the first trimmed
/// pair's class, `SNV`, as `VARIANT_CLASS`, and under `--show_ref_allele` the first
/// allele's trimmed reference, `A`, as `REF_ALLELE` (OutputFactory.pm:1357 runs on
/// the first allele before the rejoin restores the record's string and the cached
/// value serves every merged row). The default format prints the record
/// span's positions (`438-440`, `86-88`, `29-30`) on every row; the VCF format
/// prints each allele's own (`439`, `87`, `29` for the synonymous `G`), as the two
/// factories reset the cached positions after and before the rejoin
/// (OutputFactory.pm:294-299, OutputFactory/VCF.pm:305-308).
#[test]
fn a_split_record_rejoins_as_ensembl_vep_rejoins_it() {
    let dir = tempfile::tempdir().unwrap();
    let flags = ["--minimal", "--uploaded_allele", "--variant_class"];
    vep(dir.path(), MULTI_THREE, &[], &flags);
    let rows = rows_on(&data_lines(dir.path()), "ENST00000427445");
    let extra = "UPLOADED_ALLELE=AAG/GAG/AGG/-;IMPACT=";
    let shape: Vec<(String, String, String, String, String, String)> = rows
        .iter()
        .map(|r| {
            let num = r[6]
                .split(';')
                .find_map(|kv| kv.strip_prefix("ALLELE_NUM="))
                .unwrap()
                .to_string();
            (
                r[0].clone(),
                r[1].clone(),
                num,
                [r[2].as_str(), r[3].as_str(), r[4].as_str()].join("|"),
                r[5].clone(),
                r[6].clone(),
            )
        })
        .collect();
    let location = "21:10542415-10542417".to_string();
    let positions = "438-440|86-88|29-30".to_string();
    assert_eq!(
        shape
            .iter()
            .map(|(l, a, n, p, c, _)| (l.as_str(), a.as_str(), n.as_str(), p.as_str(), c.as_str()))
            .collect::<Vec<_>>(),
        [
            (
                location.as_str(),
                "G",
                "2",
                positions.as_str(),
                "synonymous_variant"
            ),
            (
                location.as_str(),
                "-",
                "3",
                positions.as_str(),
                "inframe_deletion"
            ),
            (
                location.as_str(),
                "G",
                "1",
                positions.as_str(),
                "missense_variant"
            ),
        ]
    );
    for (_, _, _, _, _, e) in &shape {
        assert!(e.contains(extra), "{e}");
        assert!(e.contains("VARIANT_CLASS=SNV;MINIMISED=1"), "{e}");
    }

    vep(
        dir.path(),
        MULTI_THREE,
        &[],
        &["--minimal", "--show_ref_allele"],
    );
    let with_reference = rows_on(&data_lines(dir.path()), "ENST00000427445");
    assert_eq!(with_reference.len(), 3);
    for row in &with_reference {
        assert!(row[6].contains("REF_ALLELE=A;"), "{}", row[6]);
    }

    vep(dir.path(), MULTI_THREE, &["--vcf"], &flags);
    let csq: Vec<String> = data_lines(dir.path())
        .pop()
        .unwrap()
        .split('\t')
        .nth(7)
        .unwrap()
        .trim_start_matches("CSQ=")
        .split(',')
        .filter(|entry| entry.contains("ENST00000427445"))
        .map(|entry| {
            let f: Vec<&str> = entry.split('|').collect();
            // Allele, cDNA_position, CDS_position, Protein_position, ALLELE_NUM.
            format!("{}|{}|{}|{}|{}", f[0], f[12], f[13], f[14], f[18])
        })
        .collect();
    assert_eq!(
        csq,
        ["G|439|87|29|2", "-|438-440|86-88|29-30|3", "G|438|86|29|1"]
    );
}

/// Three records at the acceptor site of exon 6 of ENST00000427445 whose
/// alleles fall on both sides of the CDS boundary: an SNV at 10542393 and a
/// deletion of `TACG` over 10542394-10542397, or two SNVs. Under `--minimal
/// --coding_only` vep-rs judges each split row on its allele's own span, so the
/// deletion's seven rows survive and the SNVs' none, and `REF_ALLELE` is the
/// first allele's trimmed reference. Ensembl VEP 116.2's default output writes
/// fourteen rows per record over the record's restored span (`OutputFactory.pm`
/// lines 297-323, 490-502) with the record's `GTACG` as `REF_ALLELE` on
/// `straddle_a` and `straddle_b`; the difference is documented at
/// `docs/divergences/minimal-coding-only-record-span.md`, and this test pins the
/// rows vep-rs writes while it stands.
#[test]
fn coding_only_judges_a_split_records_rows_on_their_own_spans() {
    let dir = tempfile::tempdir().unwrap();
    let records = "21\t10542393\tstraddle_a\tGTACG\tCTACG,G\t.\t.\t.\n\
                   21\t10542393\tstraddle_b\tGTACG\tCTACG,TTACG\t.\t.\t.\n\
                   21\t10542393\tstraddle_c\tGTACG\tG,CTACG\t.\t.\t.\n";
    vep(
        dir.path(),
        records,
        &[],
        &["--minimal", "--coding_only", "--show_ref_allele"],
    );
    let lines = data_lines(dir.path());
    let count = |name: &str| lines.iter().filter(|l| l.starts_with(name)).count();
    assert_eq!(
        (
            count("straddle_a\t"),
            count("straddle_b\t"),
            count("straddle_c\t")
        ),
        (7, 0, 7)
    );
    let on_transcript: Vec<(String, String, String, String, String)> =
        rows_on(&lines, "ENST00000427445")
            .iter()
            .map(|r| {
                let extra = |key: &str| {
                    r[6].split(';')
                        .find_map(|kv| kv.strip_prefix(key))
                        .unwrap()
                        .to_string()
                };
                (
                    r[0].clone(),
                    r[1].clone(),
                    r[5].clone(),
                    extra("ALLELE_NUM="),
                    extra("REF_ALLELE="),
                )
            })
            .collect();
    let location = "21:10542393-10542397";
    let consequence = "splice_acceptor_variant,coding_sequence_variant";
    assert_eq!(
        on_transcript
            .iter()
            .map(|(l, a, c, n, r)| (l.as_str(), a.as_str(), c.as_str(), n.as_str(), r.as_str()))
            .collect::<Vec<_>>(),
        [
            (location, "-", consequence, "2", "G"),
            (location, "-", consequence, "1", "TACG"),
        ]
    );
}
