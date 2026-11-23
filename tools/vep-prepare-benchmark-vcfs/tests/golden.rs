// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the fixture directories and compares every prepared file, the
//! manifest and stdout byte for byte with the goldens under `tests/fixtures/expected/`,
//! which hold the tool's contract.

mod common;

use common::{
    assert_clean_exit, filled, fixture, fixtures, json_escaped, manifest_without_its_time,
    names_in, prepare, Scratch,
};
use std::fs;

/// A golden's name, the fixture directory that is its input and the options of its run;
/// `{F}` stands for the fixtures directory.
type Case = (&'static str, &'static str, &'static [&'static str]);

/// `contigs/`: every naming mode with the filter and without, both cache layouts, named
/// sets, an assembly that names no set, and caps on the variant lines.
const CONTIGS: &[Case] = &[
    ("contigs_default", "contigs", &[]),
    ("contigs_no_filter", "contigs", &["--no-canonical-contigs"]),
    ("contigs_keep", "contigs", &["--chrom-mode", "keep"]),
    ("contigs_strip", "contigs", &["--chrom-mode", "strip"]),
    ("contigs_add", "contigs", &["--chrom-mode", "add"]),
    (
        "contigs_strip_no_filter",
        "contigs",
        &["--chrom-mode=strip", "--no-canonical-contigs"],
    ),
    (
        "contigs_add_no_filter",
        "contigs",
        &["--chrom-mode=add", "--no-canonical-contigs"],
    ),
    (
        "contigs_plain_cache",
        "contigs",
        &["--cache-dir", "{F}/cache_plain", "--assembly", "GRCh38"],
    ),
    (
        "contigs_chr_cache",
        "contigs",
        &[
            "--cache-dir",
            "{F}/cache_chr",
            "--assembly",
            "GRCh37",
            "--canonical-contigs",
        ],
    ),
    (
        "contigs_chr_cache_no_filter",
        "contigs",
        &["--cache-dir", "{F}/cache_chr", "--no-canonical-contigs"],
    ),
    (
        "contigs_plain_cache_no_filter",
        "contigs",
        &["--cache-dir", "{F}/cache_plain", "--no-canonical-contigs"],
    ),
    (
        "contigs_cache_without_contigs",
        "contigs",
        &["--cache-dir", "{F}/shapes"],
    ),
    (
        "contigs_named_list",
        "contigs",
        &[
            "--canonical-contigs",
            "2L, 2R,X,,mitochondrion_genome",
            "--assembly",
            "BDGP6.54",
        ],
    ),
    (
        "contigs_named_file",
        "contigs",
        &["--canonical-contigs", "{F}/fly_contigs.txt"],
    ),
    (
        "contigs_named_prefixed",
        "contigs",
        &[
            "--canonical-contigs=chr1,chrX,2,chr",
            "--chrom-mode",
            "keep",
        ],
    ),
    (
        "contigs_named_prefixed_add",
        "contigs",
        &[
            "--canonical-contigs=chr1,chrX,2,chr22",
            "--chrom-mode",
            "add",
        ],
    ),
    (
        "contigs_other_assembly",
        "contigs",
        &["--assembly", "BDGP6.54"],
    ),
    (
        "contigs_empty_assembly",
        "contigs",
        &["--assembly", "", "--canonical-contigs"],
    ),
    (
        "contigs_later_flag_wins",
        "contigs",
        &[
            "--no-canonical-contigs",
            "--canonical-contigs",
            "--max-variants",
            "3",
            "--max-variants=5",
        ],
    ),
    ("contigs_cap_0", "contigs", &["--max-variants", "0"]),
    ("contigs_cap_1", "contigs", &["--max-variants", "1"]),
    (
        "contigs_cap_4",
        "contigs",
        &["--max-variants", "4", "--no-canonical-contigs"],
    ),
    ("contigs_cap_negative", "contigs", &["--max-variants", "-7"]),
    (
        "contigs_cap_past_the_end",
        "contigs",
        &["--max-variants", " +1_000_000 "],
    ),
    (
        "contigs_cap_huge",
        "contigs",
        &["--max-variants", "123456789012345678901234567890"],
    ),
];

/// `shapes/`: line ends, bytes outside UTF-8, a byte order mark, empty inputs, two inputs
/// written under one name, and names that are and are not inputs.
const SHAPES: &[Case] = &[
    ("shapes_default", "shapes", &[]),
    ("shapes_no_filter", "shapes", &["--no-canonical-contigs"]),
    ("shapes_add", "shapes", &["--chrom-mode", "add"]),
    (
        "shapes_add_no_filter",
        "shapes",
        &["--chrom-mode", "add", "--no-canonical-contigs"],
    ),
    (
        "shapes_strip_no_filter",
        "shapes",
        &["--chrom-mode", "strip", "--no-canonical-contigs"],
    ),
    ("shapes_cap_2", "shapes", &["--max-variants", "2"]),
    (
        "shapes_cap_2_no_filter",
        "shapes",
        &["--max-variants", "2", "--no-canonical-contigs"],
    ),
    (
        "shapes_cap_0_no_filter",
        "shapes",
        &["--max-variants", "0", "--no-canonical-contigs"],
    ),
];

/// `gzip/`: BGZF blocks in more runs than one, members with every header field and zero
/// padding, blocks then a member, and a block the member reader reads to another end.
const GZIP: &[Case] = &[
    ("gzip_default", "gzip", &[]),
    ("gzip_no_filter", "gzip", &["--no-canonical-contigs"]),
    (
        "gzip_add",
        "gzip",
        &["--chrom-mode", "add", "--no-canonical-contigs"],
    ),
    ("gzip_cap_5", "gzip", &["--max-variants", "5"]),
    (
        "gzip_cap_200",
        "gzip",
        &["--max-variants", "200", "--no-canonical-contigs"],
    ),
    (
        "gzip_plain_cache",
        "gzip",
        &["--cache-dir", "{F}/cache_plain", "--no-canonical-contigs"],
    ),
];

/// Runs one case under a thread bound and holds stdout, the manifest and every prepared
/// file to the case's golden.
fn holds(scratch: &Scratch, (name, directory, options): Case, threads: Option<&str>) {
    let what = format!("{name} under {threads:?} threads");
    let fixtures_dir = fixtures();
    let fixtures_text = fixtures_dir.to_str().unwrap();
    let options: Vec<String> = options
        .iter()
        .map(|option| option.replace("{F}", fixtures_text))
        .collect();
    let options: Vec<&str> = options.iter().map(String::as_str).collect();
    let input = fixture(directory);
    let out = scratch.out();
    let _ = fs::remove_dir_all(&out);
    let output = prepare(&input, &out, &options, threads);
    assert_clean_exit(&output, &what);

    let golden = fixture("expected").join(name);
    let golden_text = |file: &str| fs::read(golden.join(file)).unwrap();
    let out_text = out.to_str().unwrap();
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&filled(
            &golden_text("stdout.txt"),
            &[("OUTPUT_DIR", out_text)]
        )),
        "{what}: stdout"
    );
    let expected_names: Vec<String> = names_in(&golden)
        .into_iter()
        .filter(|file| file != "stdout.txt")
        .collect();
    assert_eq!(names_in(&out), expected_names, "{what}: the files written");
    for file in &expected_names {
        if file == "manifest.json" {
            let expected = filled(
                &golden_text(file),
                &[
                    ("OUTPUT_DIR", &json_escaped(&out)),
                    ("INPUT_DIR", &json_escaped(&input)),
                    ("FIXTURES", &json_escaped(&fixtures_dir)),
                ],
            );
            assert_eq!(
                String::from_utf8_lossy(&manifest_without_its_time(&out)),
                String::from_utf8_lossy(&expected),
                "{what}: manifest.json"
            );
        } else {
            assert!(
                fs::read(out.join(file)).unwrap() == golden_text(file),
                "{what}: {file} differs from its golden"
            );
        }
    }
}

#[test]
fn contig_naming_and_filter_cases_match_their_goldens() {
    let scratch = Scratch::new("contigs");
    for case in CONTIGS {
        holds(&scratch, *case, None);
    }
}

#[test]
fn line_shape_cases_match_their_goldens() {
    let scratch = Scratch::new("shapes");
    for case in SHAPES {
        holds(&scratch, *case, None);
    }
}

#[test]
fn gzip_cases_match_their_goldens() {
    let scratch = Scratch::new("gzip");
    for case in GZIP {
        holds(&scratch, *case, None);
    }
}

#[test]
fn every_thread_bound_writes_the_same_bytes() {
    let scratch = Scratch::new("threads");
    for threads in ["1", "2", "5", " 64 "] {
        for case in SHAPES.iter().chain(GZIP).chain(&CONTIGS[..4]) {
            holds(&scratch, *case, Some(threads));
        }
    }
}

#[test]
fn every_golden_is_a_case() {
    let mut listed: Vec<&str> = CONTIGS
        .iter()
        .chain(SHAPES)
        .chain(GZIP)
        .map(|(name, _, _)| *name)
        .collect();
    listed.sort_unstable();
    assert_eq!(names_in(&fixture("expected")), listed);
}
