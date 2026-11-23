// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Files of several pieces: every thread count writes the same bytes, a plain file, a
//! gzip file of several members and a block-gzipped file of several hundred blocks read
//! the same, and the counts, digests and line sums are those of the generated text.

mod common;

use common::{line_sum, Scratch, BIN};
use flate2::write::GzEncoder;
use flate2::{Compress, Compression, FlushCompress};
use sha2::{Digest, Sha256};
use std::io::Write;
use std::path::Path;
use std::process::Command;

const THREADS: [Option<&str>; 4] = [Some("1"), Some("2"), Some("5"), None];

fn run(args: &[&str], threads: Option<&str>) -> Vec<u8> {
    let mut command = Command::new(BIN);
    command.args(args);
    match threads {
        Some(count) => command.env("VEP_FILE_CENSUS_THREADS", count),
        None => command.env_remove("VEP_FILE_CENSUS_THREADS"),
    };
    let output = command.output().unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

/// About seven mebibytes of VCF text on prefixed contigs: header lines at the start and
/// among the records, lines of every length up to a few kilobytes.
fn vcf_text() -> (Vec<u8>, u64, u64) {
    let mut text = b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec();
    let (mut records, mut with_key) = (0, 0);
    let mut at = 0u64;
    while text.len() < 7 << 20 {
        at += 1;
        if at.is_multiple_of(5000) {
            text.extend_from_slice(b"##a header line among the records\n");
            continue;
        }
        let contig = ["chr1", "chr2", "chrX"][(at / 1000 % 3) as usize];
        let filler = "x".repeat((at * 37 % 3000) as usize);
        let key = if at.is_multiple_of(3) {
            ";CSQ="
        } else {
            ";NOTE="
        };
        with_key += u64::from(at.is_multiple_of(3));
        records += 1;
        writeln!(
            text,
            "{contig}\t{at}\t.\tA\tG\t.\tPASS\tAF=0.5{key}{filler}"
        )
        .unwrap();
    }
    (text, records, with_key)
}

fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// `text` as a BGZF file: blocks of 20,000 bytes, each cut where it falls.
fn block_gzip(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for chunk in text.chunks(20_000).chain(std::iter::once(&[][..])) {
        let mut deflated = vec![0u8; chunk.len() + 1024];
        let mut compress = Compress::new(Compression::fast(), false);
        compress
            .compress(chunk, &mut deflated, FlushCompress::Finish)
            .unwrap();
        deflated.truncate(compress.total_out() as usize);
        let length = 18 + deflated.len() + 8;
        out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0]);
        out.extend_from_slice(&((length - 1) as u16).to_le_bytes());
        out.extend_from_slice(&deflated);
        out.extend_from_slice(&crc32fast::hash(chunk).to_le_bytes());
        out.extend_from_slice(&(chunk.len() as u32).to_le_bytes());
    }
    out
}

/// `text` as a gzip file of three members.
fn gzip_members(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for member in text.chunks(text.len() / 3 + 1) {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(member).unwrap();
        out.extend_from_slice(&encoder.finish().unwrap());
    }
    out
}

fn path_str(path: &Path) -> &str {
    path.to_str().unwrap()
}

#[test]
fn counts_digests_and_row_passes_are_the_same_on_every_thread_count_and_in_every_packing() {
    let scratch = Scratch::new("pieces");
    let (text, records, with_key) = vcf_text();
    let (plain, blocks, members) = (
        scratch.path("in.vcf"),
        scratch.path("in.blocks.vcf.gz"),
        scratch.path("in.members.vcf.gz"),
    );
    std::fs::write(&plain, &text).unwrap();
    std::fs::write(&blocks, block_gzip(&text)).unwrap();
    std::fs::write(&members, gzip_members(&text)).unwrap();

    let body: Vec<u8> = text
        .split_inclusive(|byte| *byte == b'\n')
        .filter(|line| !line.starts_with(b"##"))
        .flatten()
        .copied()
        .collect();
    let stripped: Vec<u8> = text
        .split_inclusive(|byte| *byte == b'\n')
        .flat_map(|line| {
            if line.starts_with(b"##contig=<ID=chr") {
                [&line[..13], &line[16..]].concat()
            } else {
                line.strip_prefix(b"chr").unwrap_or(line).to_vec()
            }
        })
        .collect();
    let counts = format!("{records}\n{with_key}\n");
    let tally: u64 = records;

    for threads in THREADS {
        for (file, gunzip) in [(&plain, false), (&blocks, true), (&members, true)] {
            let file = path_str(file);
            let packed: &[&str] = if gunzip { &["--gunzip"] } else { &[] };
            let with = |head: &[&str], tail: &[&str]| -> Vec<u8> {
                let args: Vec<&str> = head
                    .iter()
                    .chain(packed)
                    .chain(tail)
                    .copied()
                    .chain([file])
                    .collect();
                run(&args, threads)
            };
            assert_eq!(
                String::from_utf8(with(
                    &["count"],
                    &["--not-starting", "#", "--info-field", "CSQ"]
                ))
                .unwrap(),
                counts,
                "{file} on {threads:?} thread(s)"
            );
            assert_eq!(
                String::from_utf8(with(
                    &["digest", "--algorithm", "sha256"],
                    &["--without-lines-starting", "##"]
                ))
                .unwrap(),
                format!("{}  -\n", sha256_hex(&body)),
                "{file} on {threads:?} thread(s)"
            );
            if gunzip {
                assert_eq!(
                    String::from_utf8(with(&["digest", "--algorithm", "sha256"], &[])).unwrap(),
                    format!("{}  -\n", sha256_hex(&text)),
                    "{file} on {threads:?} thread(s)"
                );
            }
            assert!(
                run(&["strip-chr", file], threads) == stripped,
                "{file} on {threads:?} thread(s): strip-chr"
            );
            let contigs = String::from_utf8(with(&["contigs"], &[])).unwrap();
            let total: u64 = contigs
                .lines()
                .map(|line| line.split('\t').nth(1).unwrap().parse::<u64>().unwrap())
                .sum();
            assert_eq!(total, tally, "{file} on {threads:?} thread(s): contigs");
            assert_eq!(
                contigs
                    .lines()
                    .map(|line| line.split('\t').next().unwrap())
                    .collect::<Vec<_>>(),
                ["chr1", "chr2", "chrX"]
            );
        }
        assert_eq!(
            String::from_utf8(run(
                &["digest", "--algorithm", "sha256", path_str(&plain)],
                threads
            ))
            .unwrap(),
            format!("{}  {}\n", sha256_hex(&text), path_str(&plain))
        );
    }
}

#[test]
fn the_line_sum_is_the_same_on_every_thread_count_in_every_packing_and_in_any_order() {
    let scratch = Scratch::new("line-sum");
    let (text, records, _) = vcf_text();
    // The same lines from the last to the first, each ended by `\r\n` but every seventh,
    // which ends in a lone `\r`.
    let mut reversed = Vec::with_capacity(text.len() + text.len() / 64);
    for (at, line) in text.split(|byte| *byte == b'\n').rev().enumerate() {
        if line.is_empty() {
            continue;
        }
        reversed.extend_from_slice(line);
        reversed.extend_from_slice(if at % 7 == 0 { b"\r" } else { b"\r\n" });
    }
    let (plain, blocks, members, turned) = (
        scratch.path("in.vcf"),
        scratch.path("in.blocks.vcf.gz"),
        scratch.path("in.members.vcf.gz"),
        scratch.path("reversed.vcf"),
    );
    std::fs::write(&plain, &text).unwrap();
    std::fs::write(&blocks, block_gzip(&text)).unwrap();
    std::fs::write(&members, gzip_members(&text)).unwrap();
    std::fs::write(&turned, &reversed).unwrap();

    let lines = || {
        text.split(|byte| *byte == b'\n')
            .filter(|line| !line.is_empty())
    };
    let every = line_sum(lines());
    let rows = line_sum(lines().filter(|line| !line.starts_with(b"#")));
    assert!(rows.starts_with(&format!("{records}\t")));
    assert_ne!(every, rows);

    for threads in THREADS {
        for (file, gunzip) in [
            (&plain, false),
            (&blocks, true),
            (&members, true),
            (&turned, false),
        ] {
            let file = path_str(file);
            let packed: &[&str] = if gunzip { &["--gunzip"] } else { &[] };
            let with = |options: &[&str]| -> String {
                let args: Vec<&str> = ["line-sum"]
                    .iter()
                    .chain(packed)
                    .chain(options)
                    .copied()
                    .chain([file])
                    .collect();
                String::from_utf8(run(&args, threads)).unwrap()
            };
            assert_eq!(with(&[]), every, "{file} on {threads:?} thread(s)");
            assert_eq!(
                with(&["--without-lines-starting", "#"]),
                rows,
                "{file} on {threads:?} thread(s)"
            );
        }
        // One line a file, in the order of the operands.
        assert_eq!(
            String::from_utf8(run(
                &["line-sum", path_str(&turned), path_str(&plain)],
                threads
            ))
            .unwrap(),
            format!("{every}{every}")
        );
    }
}

#[test]
fn a_digest_in_chunks_and_a_body_after_the_header_do_not_depend_on_the_thread_count() {
    let scratch = Scratch::new("chunks");
    let (text, _, _) = vcf_text();
    let plain = scratch.path("in.vcf");
    std::fs::write(&plain, &text).unwrap();
    let header = text
        .split_inclusive(|byte| *byte == b'\n')
        .take_while(|line| line.starts_with(b"##"))
        .map(<[u8]>::len)
        .sum::<usize>();
    let folded: String = text[header..]
        .chunks(1 << 20)
        .map(|chunk| sha256_hex(chunk) + "\n")
        .collect();
    for threads in THREADS {
        assert_eq!(
            String::from_utf8(run(
                &[
                    "digest",
                    "--algorithm",
                    "sha256",
                    "--after-leading-lines-starting",
                    "##",
                    "--chunk-bytes",
                    "1048576",
                    path_str(&plain)
                ],
                threads
            ))
            .unwrap(),
            format!("{}  -\n", sha256_hex(folded.as_bytes()))
        );
        assert_eq!(
            String::from_utf8(run(
                &[
                    "digest",
                    "--algorithm",
                    "sha256",
                    "--after-leading-lines-starting",
                    "##",
                    path_str(&plain)
                ],
                threads
            ))
            .unwrap(),
            format!("{}  -\n", sha256_hex(&text[header..]))
        );
    }
}

#[test]
fn the_metrics_of_an_output_of_several_pieces_are_the_same_on_every_thread_count() {
    let scratch = Scratch::new("metrics");
    let mut text = b"## header\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n".to_vec();
    let mut rows = 0u64;
    while text.len() < 6 << 20 {
        rows += 1;
        let term = [
            "missense_variant",
            "intron_variant",
            "intergenic_variant,upstream_gene_variant",
        ][(rows % 3) as usize];
        let extra = if rows.is_multiple_of(4) {
            "IMPACT=LOW;CLIN_SIG=benign;AF=0.1"
        } else {
            "IMPACT=LOW"
        };
        let end = if rows.is_multiple_of(1000) {
            "\r\n"
        } else {
            "\n"
        };
        write!(text, "v{rows}\t1:{rows}\tA\tG1\tT1\tTranscript\t{term}\t-\t-\t-\t-\t-\trs{rows}\t{extra}{end}").unwrap();
    }
    let output = scratch.path("output.txt");
    std::fs::write(&output, &text).unwrap();
    let mut reports = Vec::new();
    for (at, threads) in THREADS.into_iter().enumerate() {
        let report = scratch.path(&format!("report-{at}.json"));
        let summary = String::from_utf8(run(
            &["metrics", path_str(&output), path_str(&report), "9"],
            threads,
        ))
        .unwrap();
        assert!(
            summary.starts_with(&format!(
                "  Total annotation lines: {rows}\n  Unique consequence types: 4\n"
            )),
            "{summary}"
        );
        reports.push(std::fs::read_to_string(&report).unwrap());
    }
    assert!(reports[0].contains(&format!("\"clin_sig_count\": {},", rows / 4)));
    assert!(reports.iter().all(|report| *report == reports[0]));
}
