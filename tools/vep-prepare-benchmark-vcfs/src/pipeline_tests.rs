// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two ways an input is prepared, on one thread and across threads, each held to the
//! goldens under `tests/fixtures/expected/`: the binary takes the fixtures, which are
//! small, on one thread, so the threaded way is run on them here.

use super::*;
use crate::contigs::ContigSet;
use crate::lines::Mode;
use crate::paths::{is_gzip_name, list_inputs, prepared_name};
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::AtomicBool;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// The contig set of a golden's run.
enum Set {
    Off,
    Human,
    Named(&'static str),
}

/// A golden, the fixture directory that is its input, the contig set and the cap of its
/// run.
const CASES: [(&str, &str, Set, Option<u64>); 22] = [
    ("gzip_default", "gzip", Set::Human, None),
    ("gzip_no_filter", "gzip", Set::Off, None),
    ("gzip_add", "gzip", Set::Off, None),
    ("gzip_cap_5", "gzip", Set::Human, Some(5)),
    ("gzip_cap_200", "gzip", Set::Off, Some(200)),
    ("gzip_plain_cache", "gzip", Set::Off, None),
    ("shapes_default", "shapes", Set::Human, None),
    ("shapes_no_filter", "shapes", Set::Off, None),
    ("shapes_add", "shapes", Set::Human, None),
    ("shapes_add_no_filter", "shapes", Set::Off, None),
    ("shapes_strip_no_filter", "shapes", Set::Off, None),
    ("shapes_cap_2", "shapes", Set::Human, Some(2)),
    ("shapes_cap_2_no_filter", "shapes", Set::Off, Some(2)),
    ("shapes_cap_0_no_filter", "shapes", Set::Off, Some(0)),
    ("contigs_default", "contigs", Set::Human, None),
    ("contigs_plain_cache", "contigs", Set::Human, None),
    ("contigs_chr_cache_no_filter", "contigs", Set::Off, None),
    (
        "contigs_named_prefixed",
        "contigs",
        Set::Named("chr1,chrX,2,chr"),
        None,
    ),
    (
        "contigs_named_prefixed_add",
        "contigs",
        Set::Named("chr1,chrX,2,chr22"),
        None,
    ),
    ("contigs_cap_0", "contigs", Set::Human, Some(0)),
    ("contigs_cap_1", "contigs", Set::Human, Some(1)),
    ("contigs_cap_4", "contigs", Set::Off, Some(4)),
];

/// What a golden manifest records for one input.
struct Entry {
    style: String,
    mode: Mode,
    counts: Counts,
    sha256: String,
}

fn entry(manifest: &str, input: &str) -> Entry {
    let start = manifest
        .find(&format!("\"input_file\": \"{{{{INPUT_DIR}}}}/{input}\","))
        .unwrap_or_else(|| panic!("the manifest lists {input}"));
    let block = &manifest[start..];
    let block = &block[..block.find("\n    }").expect("the entry ends")];
    let field = |key: &str| -> &str {
        let lead = format!("\"{key}\": ");
        let rest = &block[block.find(&lead).expect("the entry holds the key") + lead.len()..];
        rest[..rest.find([',', '\n']).unwrap_or(rest.len())].trim_matches('"')
    };
    let number = |key: &str| -> u64 { field(key).parse().expect("a count") };
    Entry {
        style: field("detected_input_chrom_style").to_string(),
        mode: match field("normalization_mode") {
            "keep" => Mode::Keep,
            "strip" => Mode::Strip,
            "add" => Mode::Add,
            other => panic!("mode {other}"),
        },
        counts: Counts {
            headers: number("headers_written"),
            variants: number("variants_written"),
            dropped: number("variants_dropped_non_canonical"),
            chr_stripped: field("chr_stripped") == "true",
        },
        sha256: field("output_sha256").to_string(),
    }
}

#[test]
fn both_ways_of_preparing_an_input_write_what_the_goldens_hold() {
    let scratch = std::env::temp_dir().join(format!(
        "vep-prepare-benchmark-vcfs-pipeline-{}",
        std::process::id()
    ));
    fs::create_dir_all(&scratch).unwrap();
    let pool = Pool::new(3);
    let mut prepared_files = 0;
    for (case, directory, set, cap) in CASES {
        let golden = fixtures().join("expected").join(case);
        let manifest = fs::read_to_string(golden.join("manifest.json")).unwrap();
        let contigs = match set {
            Set::Off => None,
            Set::Human => Some(Arc::new(ContigSet::human())),
            Set::Named(value) => Some(Arc::new(ContigSet::named(value.as_bytes()).unwrap())),
        };
        let input_dir = fixtures().join(directory);
        let names = list_inputs(input_dir.to_str().unwrap().as_bytes());
        assert!(!names.is_empty(), "{directory}");
        for name in &names {
            let name = std::str::from_utf8(name).unwrap();
            let expected = entry(&manifest, name);
            let src = input_dir.join(name);
            let gzip = is_gzip_name(name.as_bytes());
            assert_eq!(
                first_variant_style(&src, gzip).unwrap().name(),
                expected.style,
                "{case}: {name}"
            );
            let rules = Arc::new(Rules {
                mode: expected.mode,
                contigs: contigs.clone(),
            });
            let preparation = Preparation {
                rules: &rules,
                cap,
                pool: &pool,
                in_flight: 2,
            };
            for threaded in [false, true] {
                let what = format!("{case}: {name}, across threads: {threaded}");
                let dst = scratch.join("prepared.vcf");
                let (source, target) = (File::open(&src).unwrap(), File::create(&dst).unwrap());
                let prepared = if threaded {
                    preparation.across_threads(source, target, &src, gzip)
                } else {
                    preparation.on_one_thread(source, target, gzip)
                }
                .unwrap_or_else(|e| panic!("{what}: {e}"));
                assert_eq!(prepared.counts, expected.counts, "{what}");
                assert_eq!(prepared.sha256, expected.sha256, "{what}");
                let written = fs::read(&dst).unwrap();
                assert_eq!(hex(&Sha256::digest(&written)), expected.sha256, "{what}");
                // Where no other input is written under the name, the golden file is it.
                let prepared_as = std::str::from_utf8(prepared_name(name.as_bytes())).unwrap();
                let shared = names
                    .iter()
                    .filter(|other| prepared_name(other) == prepared_as.as_bytes())
                    .count()
                    > 1;
                if !shared {
                    assert!(
                        written == fs::read(golden.join(prepared_as)).unwrap(),
                        "{what}"
                    );
                }
                prepared_files += 1;
            }
        }
    }
    assert_eq!(prepared_files, 2 * (6 * 5 + 8 * 12 + 8 * 2));
    fs::remove_dir_all(&scratch).unwrap();
}

#[test]
fn a_failure_in_the_middle_of_an_input_is_returned_on_either_way() {
    let scratch = std::env::temp_dir().join(format!(
        "vep-prepare-benchmark-vcfs-pipeline-failure-{}",
        std::process::id()
    ));
    fs::create_dir_all(&scratch).unwrap();
    let good = fs::read(fixtures().join("gzip/blocks.vcf.gz")).unwrap();
    let src = scratch.join("cut_off.vcf.gz");
    let dst = scratch.join("cut_off.vcf");
    let pool = Pool::new(2);
    let rules = Arc::new(Rules {
        mode: Mode::Keep,
        contigs: None,
    });
    let preparation = Preparation {
        rules: &rules,
        cap: None,
        pool: &pool,
        in_flight: 2,
    };
    // The end of the file is missing, then the checksum of its first block is another.
    let mut wrong_byte = good.clone();
    wrong_byte[bgzf_block_len(&good).expect("a BGZF fixture") - 8] ^= 0x55;
    for bytes in [&good[..good.len() - 40], &wrong_byte[..]] {
        fs::write(&src, bytes).unwrap();
        for threaded in [false, true] {
            let (source, target) = (File::open(&src).unwrap(), File::create(&dst).unwrap());
            let prepared = if threaded {
                preparation.across_threads(source, target, &src, true)
            } else {
                preparation.on_one_thread(source, target, true)
            };
            assert!(prepared.is_err(), "across threads: {threaded}");
        }
    }
    fs::remove_dir_all(&scratch).unwrap();
}

#[test]
fn worker_threads_start_with_the_jobs_and_stay_within_the_bound() {
    let pool = Pool::new(3);
    assert_eq!(pool.threads.lock().unwrap().len(), 0);
    let ran = Arc::new(AtomicUsize::new(0));
    let inflated = Arc::new(AtomicBool::new(false));
    let (done, results) = channel();
    for _ in 0..40 {
        let (ran, inflated, done) = (Arc::clone(&ran), Arc::clone(&inflated), done.clone());
        pool.submit(Box::new(move |worker| {
            inflated.fetch_or(worker.inflater.is_some(), Ordering::Relaxed);
            ran.fetch_add(1, Ordering::Relaxed);
            let _ = done.send(());
        }));
    }
    for _ in 0..40 {
        results.recv().unwrap();
    }
    assert_eq!(ran.load(Ordering::Relaxed), 40);
    assert!(
        !inflated.load(Ordering::Relaxed),
        "no job asked for an inflater"
    );
    let started = pool.threads.lock().unwrap().len();
    assert!((1..=3).contains(&started), "{started} threads");
}
