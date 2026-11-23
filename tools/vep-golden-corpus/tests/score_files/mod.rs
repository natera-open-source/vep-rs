// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Score files of the SpliceAI plugin as the tests write them: a bgzipped VCF of chromosome
//! 21 whose header lines are one BGZF block and whose every record line is a block of its
//! own, with the tabix index beside it.

#![allow(dead_code)]

use flate2::write::{DeflateEncoder, GzEncoder};
use flate2::{Compression, Crc};
use std::io::Write;
use std::path::Path;

/// The header lines of a score file.
pub const HEADER: [&str; 3] = [
    "##fileformat=VCFv4.2",
    "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAIv1.3 synthetic fixture, \
     fabricated values. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
];

/// Made-up scores of two symbols, each a symbol and its eight values.
pub const GENEA: &str = "GENEA|0.15|0.08|0.06|0.17|-6|-41|-17|-9";
pub const GENEB: &str = "GENEB|0.00|1.00|0.10|0.05|-25|-23|37|17";

/// A deletion of two bases, whose rows stand at 21:5001-5002 with the allele written `-`,
/// and a substitution of one base, whose rows stand at 21:7000 with the alternate base:
/// each a POS, a REF and an ALT.
pub const DELETION: (&str, &str, &str) = ("5000", "GCA", "G");
pub const SUBSTITUTION: (&str, &str, &str) = ("7000", "A", "C");

/// One BGZF block that holds `data`.
pub fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut encoder = DeflateEncoder::new(Vec::new(), Compression::new(6));
    encoder.write_all(data).unwrap();
    let body = encoder.finish().unwrap();
    let mut checksum = Crc::new();
    checksum.update(data);
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00".to_vec();
    block.extend(u16::try_from(body.len() + 25).unwrap().to_le_bytes());
    block.extend(body);
    block.extend(checksum.sum().to_le_bytes());
    block.extend(u32::try_from(data.len()).unwrap().to_le_bytes());
    block
}

/// One record line at a variant of chromosome 21: the entries of the line, each a symbol
/// and its eight values, behind the alternate allele.
pub fn record_line(variant: (&str, &str, &str), entries: &[&str]) -> String {
    let (position, reference, alternate) = variant;
    let entries: Vec<String> = entries
        .iter()
        .map(|entry| format!("{alternate}|{entry}"))
        .collect();
    format!(
        "21\t{position}\t.\t{reference}\t{alternate}\t.\t.\tSpliceAI={}",
        entries.join(",")
    )
}

/// The bytes of a tabix index of a VCF whose one sequence is `21` and whose one bin holds
/// one chunk, from the virtual offset `first` to `last`.
pub fn tabix_index(first: u64, last: u64) -> Vec<u8> {
    let mut index = b"TBI\x01".to_vec();
    for value in [1i32, 0, 1, 2, 2, i32::from(b'#'), 0, 3] {
        index.extend(value.to_le_bytes());
    }
    index.extend(b"21\0");
    index.extend(1i32.to_le_bytes());
    index.extend(0u32.to_le_bytes());
    index.extend(1i32.to_le_bytes());
    index.extend(first.to_le_bytes());
    index.extend(last.to_le_bytes());
    index.extend(0i32.to_le_bytes());
    index
}

/// Writes `bytes` gzipped at `path`.
pub fn write_gzip(path: &str, bytes: &[u8]) {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    std::fs::write(path, encoder.finish().unwrap()).unwrap();
}

/// Writes the score file `name` under `directory` with the record lines given, and its
/// index; the file's path.
pub fn write_score_file(directory: &Path, name: &str, lines: &[String]) -> String {
    std::fs::create_dir_all(directory).unwrap();
    let header = bgzf_block((HEADER.join("\n") + "\n").as_bytes());
    let body: Vec<u8> = lines
        .iter()
        .flat_map(|line| bgzf_block(format!("{line}\n").as_bytes()))
        .collect();
    let path = directory.join(name).to_str().unwrap().to_string();
    let mut file = header.clone();
    file.extend(&body);
    file.extend(bgzf_block(b""));
    std::fs::write(&path, file).unwrap();
    let index = tabix_index(
        (header.len() as u64) << 16,
        ((header.len() + body.len()) as u64) << 16,
    );
    write_gzip(&format!("{path}.tbi"), &index);
    path
}

/// The score file of a variant under `directory`: one record line of the two symbols'
/// entries (`joined`), or one record line a symbol. A substitution's file is `snv.synthetic.vcf.gz`
/// and any other `indel.synthetic.vcf.gz`.
pub fn score_file(directory: &Path, variant: (&str, &str, &str), joined: bool) -> String {
    let lines = if joined {
        vec![record_line(variant, &[GENEA, GENEB])]
    } else {
        vec![
            record_line(variant, &[GENEA]),
            record_line(variant, &[GENEB]),
        ]
    };
    write_score_file(directory, file_name(variant), &lines)
}

/// The name of the score file that holds a variant.
pub fn file_name(variant: (&str, &str, &str)) -> &'static str {
    if variant.1.len() == 1 && variant.2.len() == 1 {
        "snv.synthetic.vcf.gz"
    } else {
        "indel.synthetic.vcf.gz"
    }
}
