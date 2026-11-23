// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary over inputs made here: files larger than the pieces an input is read in, in
//! each shape (plain, one gzip member, BGZF blocks), under several thread bounds; names
//! outside UTF-8; an input directory that is the output directory; and a directory named
//! through a link.

mod common;

use common::{assert_clean_exit, names_in, prepare, stdout, Scratch};
use flate2::{Compress, Compression, FlushCompress, Status};
use sha2::{Digest, Sha256};
use std::ffi::OsStr;
use std::fs;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

fn deflate(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; data.len() + data.len() / 8 + 64];
    let mut compress = Compress::new(Compression::fast(), false);
    let status = compress
        .compress(data, &mut out, FlushCompress::Finish)
        .unwrap();
    assert_eq!(status, Status::StreamEnd);
    out.truncate(compress.total_out() as usize);
    out
}

fn trailer(data: &[u8]) -> Vec<u8> {
    [
        crc32fast::hash(data).to_le_bytes(),
        (data.len() as u32).to_le_bytes(),
    ]
    .concat()
}

/// One gzip member that holds `data`.
fn member(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255];
    out.extend_from_slice(&deflate(data));
    out.extend_from_slice(&trailer(data));
    out
}

/// `data` as BGZF: blocks of `block_bytes` and the empty block that ends the file.
fn bgzf(data: &[u8], block_bytes: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for part in data.chunks(block_bytes).chain([&[][..]]) {
        let deflated = deflate(part);
        let length = 18 + deflated.len() + 8;
        out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0]);
        out.extend_from_slice(&((length - 1) as u16).to_le_bytes());
        out.extend_from_slice(&deflated);
        out.extend_from_slice(&trailer(part));
    }
    out
}

/// A made-up input of `records` variant lines after three header lines, and what the
/// default run writes for it: every seventh record is on a contig that is dropped, every
/// eleventh ends in CR LF, and the others are on `chr21`, which is written as `21`.
struct Generated {
    input: Vec<u8>,
    prepared: Vec<u8>,
    variants: usize,
    dropped: usize,
}

fn generated(records: usize) -> Generated {
    let header = b"##fileformat=VCFv4.2\n##source=made-up records\n#CHROM\tPOS\tID\tREF\tALT\n";
    let mut input = header.to_vec();
    let mut prepared = header.to_vec();
    let (mut variants, mut dropped) = (0, 0);
    for record in 0..records {
        let rest = format!(
            "\t{}\t.\tA\tC\tAC={};AN=152312;NOTE={}",
            5_000_000 + 3 * record,
            record % 97,
            "x".repeat(record % 211)
        );
        let end: &[u8] = if record % 11 == 0 { b"\r\n" } else { b"\n" };
        if record % 7 == 3 {
            input.extend_from_slice(b"chr21_GL383579v2_alt");
            dropped += 1;
        } else {
            input.extend_from_slice(b"chr21");
            prepared.extend_from_slice(b"21");
            prepared.extend_from_slice(rest.as_bytes());
            prepared.push(b'\n');
            variants += 1;
        }
        input.extend_from_slice(rest.as_bytes());
        input.extend_from_slice(end);
    }
    Generated {
        input,
        prepared,
        variants,
        dropped,
    }
}

fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

fn manifest(output_dir: &Path) -> String {
    fs::read_to_string(output_dir.join("manifest.json")).unwrap()
}

#[test]
fn inputs_larger_than_a_piece_are_prepared_alike_in_every_shape_and_under_every_thread_bound() {
    let scratch = Scratch::new("large");
    let made = generated(60_000);
    assert!(made.input.len() > 8 << 20, "{} bytes", made.input.len());
    let input = scratch.at("in");
    fs::create_dir(&input).unwrap();
    fs::write(input.join("blocks.vcf.gz"), bgzf(&made.input, 65_280)).unwrap();
    fs::write(input.join("small_blocks.vcf.gz"), bgzf(&made.input, 4_099)).unwrap();
    fs::write(input.join("member.vcf.gz"), member(&made.input)).unwrap();
    fs::write(input.join("plain.vcf"), &made.input).unwrap();
    let digest = sha256_hex(&made.prepared);
    for threads in [None, Some("1"), Some("2"), Some("7")] {
        let _ = fs::remove_dir_all(scratch.out());
        let output = prepare(&input, &scratch.out(), &[], threads);
        assert_clean_exit(&output, "generated inputs");
        let counts = format!(
            "(mode=keep, variants={}, headers=3, canonical_dropped={} chr_stripped=true)\n",
            made.variants, made.dropped
        );
        assert_eq!(
            stdout(&output),
            format!(
                "prepared plain.vcf -> plain.vcf {counts}prepared blocks.vcf.gz -> blocks.vcf \
                 {counts}prepared member.vcf.gz -> member.vcf {counts}prepared small_blocks.vcf.gz \
                 -> small_blocks.vcf {counts}wrote manifest: {}/manifest.json\n",
                scratch.out().display()
            )
        );
        for name in ["plain.vcf", "blocks.vcf", "member.vcf", "small_blocks.vcf"] {
            assert!(
                fs::read(scratch.out().join(name)).unwrap() == made.prepared,
                "{name} under {threads:?} threads"
            );
        }
        let manifest = manifest(&scratch.out());
        assert_eq!(
            manifest
                .matches(&format!("\"output_sha256\": \"{digest}\""))
                .count(),
            4
        );
        assert_eq!(
            manifest
                .matches(&format!("\"variants_written\": {},", made.variants))
                .count(),
            4
        );
    }
}

#[test]
fn a_large_input_that_nothing_changes_is_written_as_it_is_read() {
    let scratch = Scratch::new("unchanged");
    let made = generated(60_000);
    // Without the CR LF line ends every line stays as read once the filter is off.
    let text: Vec<u8> = made.input.iter().copied().filter(|b| *b != b'\r').collect();
    let input = scratch.at("in");
    fs::create_dir(&input).unwrap();
    fs::write(input.join("blocks.vcf.gz"), bgzf(&text, 65_280)).unwrap();
    fs::write(input.join("plain.vcf"), &text).unwrap();
    for threads in [None, Some("1")] {
        let _ = fs::remove_dir_all(scratch.out());
        let output = prepare(&input, &scratch.out(), &["--no-canonical-contigs"], threads);
        assert_clean_exit(&output, "unchanged inputs");
        assert!(fs::read(scratch.out().join("blocks.vcf")).unwrap() == text);
        assert!(fs::read(scratch.out().join("plain.vcf")).unwrap() == text);
        let digest = format!("\"output_sha256\": \"{}\"", sha256_hex(&text));
        assert_eq!(manifest(&scratch.out()).matches(&digest).count(), 2);
    }
}

#[test]
fn a_cap_keeps_the_first_variant_lines_of_a_large_input() {
    let scratch = Scratch::new("capped");
    let made = generated(60_000);
    let input = scratch.at("in");
    fs::create_dir(&input).unwrap();
    fs::write(input.join("blocks.vcf.gz"), bgzf(&made.input, 65_280)).unwrap();
    fs::write(input.join("plain.vcf"), &made.input).unwrap();
    let cap = 30_000;
    let kept_end = made
        .prepared
        .iter()
        .enumerate()
        .filter(|(_, byte)| **byte == b'\n')
        .nth(3 + cap - 1)
        .map(|(at, _)| at + 1)
        .unwrap();
    for threads in [None, Some("1")] {
        let _ = fs::remove_dir_all(scratch.out());
        let output = prepare(
            &input,
            &scratch.out(),
            &["--max-variants", "30000"],
            threads,
        );
        assert_clean_exit(&output, "capped inputs");
        for name in ["blocks.vcf", "plain.vcf"] {
            assert!(
                fs::read(scratch.out().join(name)).unwrap() == made.prepared[..kept_end],
                "{name}"
            );
        }
        // The records dropped before the cap was reached are counted, the others not.
        let dropped = (0..)
            .scan(0, |kept, record| {
                (*kept < cap).then(|| {
                    let drop = record % 7 == 3;
                    *kept += usize::from(!drop);
                    drop
                })
            })
            .filter(|drop| *drop)
            .count();
        assert!(stdout(&output).contains(&format!(
            "(mode=keep, variants=30000, headers=3, canonical_dropped={dropped} chr_stripped=true)"
        )));
    }
}

#[test]
fn a_name_outside_utf8_is_written_in_its_bytes_and_ordered_by_code_point() {
    let scratch = Scratch::new("names");
    let input = scratch.at("in");
    fs::create_dir(&input).unwrap();
    let text = b"#h\nchr1\t1\n";
    for name in [&b"\xffin.vcf"[..], "\u{e000}.vcf".as_bytes(), b"a.vcf"] {
        fs::write(input.join(OsStr::from_bytes(name)), text).unwrap();
    }
    let output = prepare(&input, &scratch.out(), &[], None);
    assert_clean_exit(&output, "names outside UTF-8");
    let counts = " (mode=keep, variants=1, headers=1, canonical_dropped=0 chr_stripped=true)\n";
    let expected = [
        b"prepared a.vcf -> a.vcf".as_slice(),
        counts.as_bytes(),
        b"prepared \xffin.vcf -> \xffin.vcf",
        counts.as_bytes(),
        "prepared \u{e000}.vcf -> \u{e000}.vcf".as_bytes(),
        counts.as_bytes(),
        format!(
            "wrote manifest: {}/manifest.json\n",
            scratch.out().display()
        )
        .as_bytes(),
    ]
    .concat();
    assert!(output.stdout == expected, "{}", stdout(&output));
    let manifest = manifest(&scratch.out());
    let escaped = concat!(r"/\udcffin.vcf", "\"");
    assert_eq!(manifest.matches(escaped).count(), 2, "{manifest}");
    let private_use = concat!(r"/", r"\", "ue000.vcf\"");
    assert_eq!(manifest.matches(private_use).count(), 2, "{manifest}");
    assert!(scratch
        .out()
        .join(OsStr::from_bytes(b"\xffin.vcf"))
        .is_file());
}

#[test]
fn an_input_that_is_its_own_output_is_cut_to_nothing_as_it_is_opened() {
    let scratch = Scratch::new("same-dir");
    let dir = scratch.at("both");
    fs::create_dir(&dir).unwrap();
    fs::write(dir.join("x.vcf"), b"#h\nchr1\t1\n").unwrap();
    fs::write(dir.join("y.vcf"), b"#h\nchr1\t1\nchr2\t2\n").unwrap();
    fs::write(dir.join("y.vcf.gz"), member(b"#h\n#i\n3\t3\n")).unwrap();
    for threads in [None, Some("1")] {
        let output = prepare(&dir, &dir, &[], threads);
        assert_clean_exit(&output, "one directory");
        // The naming of `y.vcf` is read before it is cut; on the second run it is what the
        // first run wrote from `y.vcf.gz`, a plain name, and the lines are the same.
        assert_eq!(
            stdout(&output),
            format!(
                "prepared x.vcf -> x.vcf (mode=keep, variants=0, headers=0, canonical_dropped=0)\n\
                 prepared y.vcf -> y.vcf (mode=keep, variants=0, headers=0, canonical_dropped=0)\n\
                 prepared y.vcf.gz -> y.vcf (mode=keep, variants=1, headers=2, canonical_dropped=0)\n\
                 wrote manifest: {}/manifest.json\n",
                dir.display()
            )
        );
        assert_eq!(fs::read(dir.join("x.vcf")).unwrap(), b"");
        assert_eq!(fs::read(dir.join("y.vcf")).unwrap(), b"#h\n#i\n3\t3\n");
        assert_eq!(
            names_in(&dir),
            ["manifest.json", "x.vcf", "y.vcf", "y.vcf.gz"]
        );
        let empty = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
        assert_eq!(manifest(&dir).matches(empty).count(), 2);
    }
}

#[test]
fn the_manifest_names_the_directories_behind_links_and_dots() {
    let scratch = Scratch::new("links");
    let real = scratch.at("real_in");
    fs::create_dir(&real).unwrap();
    fs::write(real.join("a.vcf"), b"#h\n1\t1\n").unwrap();
    std::os::unix::fs::symlink("real_in", scratch.at("linked_in")).unwrap();
    std::os::unix::fs::symlink(scratch.at("real_out"), scratch.at("linked_out")).unwrap();
    let input = scratch.at("linked_in/../linked_in/.");
    let output_dir = scratch.at("linked_out/below/../nested//");
    let output = prepare(
        &input,
        &output_dir,
        &["--cache-dir", "tests/../tests/fixtures/cache_chr/"],
        None,
    );
    assert_clean_exit(&output, "linked directories");
    let nested = scratch.at("real_out/nested");
    assert_eq!(
        stdout(&output),
        format!(
            "prepared a.vcf -> a.vcf (mode=add, variants=1, headers=1, canonical_dropped=0 chr_stripped=true)\n\
             wrote manifest: {}/manifest.json\n",
            nested.display()
        )
    );
    let manifest = manifest(&nested);
    for line in [
        format!("  \"input_dir\": \"{}\",\n", real.display()),
        format!("  \"output_dir\": \"{}\",\n", nested.display()),
        format!(
            "  \"cache_dir\": \"{}\",\n",
            common::fixture("cache_chr").display()
        ),
        format!("      \"input_file\": \"{}/a.vcf\",\n", real.display()),
        format!("      \"output_file\": \"{}/a.vcf\",\n", nested.display()),
        "      \"detected_cache_chrom_style\": \"chr\",\n".to_string(),
    ] {
        assert!(manifest.contains(&line), "{line} in {manifest}");
    }
    assert!(!scratch.at("real_out/below").exists());
}
