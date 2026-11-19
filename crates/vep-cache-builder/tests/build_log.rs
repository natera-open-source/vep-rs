// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The build's log on standard error: one line per chromosome, in name order.

use std::process::Command;

/// A build logs one `chr <name>: <n> transcripts in <m> regions` line per
/// chromosome it writes, in the byte order of the names, whatever order the
/// model lists them in: eight chromosomes are logged `1`, `10`, `11`, `2`,
/// `20`, `3`, `X`, `Y`.
#[test]
fn chromosome_lines_are_logged_in_name_order() {
    let dir = tempfile::tempdir().unwrap();
    let mut model = String::from("##gff-version 3\n");
    for name in ["10", "2", "X", "1", "20", "3", "Y", "11"] {
        model += &format!(
            "{name}\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G{name};biotype=lncRNA\n\
             {name}\tx\tlnc_RNA\t1000\t2000\t.\t+\t.\tID=transcript:T{name};Parent=gene:G{name};biotype=lncRNA\n\
             {name}\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T{name}\n"
        );
    }
    let gff3 = dir.path().join("model.gff3");
    std::fs::write(&gff3, model).unwrap();
    let peptides = dir.path().join("peptides.fa");
    std::fs::write(&peptides, ">ENSP1 pep transcript:T1 gene:G1\nMK\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh38", "--gff3"])
        .arg(&gff3)
        .arg("--fasta")
        .arg(&peptides)
        .arg("--output-dir")
        .arg(dir.path().join("cache"))
        .output()
        .expect("vep-cache-builder runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let logged: Vec<&str> = stderr
        .lines()
        .filter(|line| line.ends_with(" transcripts in 1 regions"))
        .filter_map(|line| line.rsplit_once("chr ")?.1.split(':').next())
        .collect();
    assert_eq!(
        logged,
        ["1", "10", "11", "2", "20", "3", "X", "Y"],
        "{stderr}"
    );
}
