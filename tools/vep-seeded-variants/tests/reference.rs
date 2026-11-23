// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Where a record's reference allele comes from: the cache's own bases inside a coding
//! sequence, and with `--fasta` the reference file's bases at every record.

mod common;

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Read;

use common::{
    as_read, coding_index, complete, pilot, pilot_cache, records_of, run, run_to_bytes, Genome,
    Scratch, PILOTS,
};
use flate2::read::MultiGzDecoder;
use vep_seeded_variants::fasta::Fasta;
use vep_seeded_variants::generate::Records;

/// Every sequence of a gzip FASTA file, read whole and without its index.
fn whole_fasta(path: &std::path::Path) -> BTreeMap<String, Vec<u8>> {
    let mut text = String::new();
    MultiGzDecoder::new(File::open(path).expect("the reference file"))
        .read_to_string(&mut text)
        .expect("the reference file's text");
    let mut sequences: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    let mut name = String::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            name = header.split_whitespace().next().unwrap_or("").to_owned();
        } else {
            sequences
                .entry(name.clone())
                .or_default()
                .extend(line.trim().bytes().map(|base| base.to_ascii_uppercase()));
        }
    }
    sequences
}

#[test]
fn the_reference_allele_of_an_exonic_record_is_the_caches_and_with_a_fasta_the_genomes() {
    // Without a reference file: a reference base inside a complete coding sequence is
    // the base the cache gives there.
    let mut inside = 0;
    for (corpus, _) in PILOTS {
        let dir = pilot_cache(corpus);
        let genome = Genome::load(&dir);
        for record in &records_of(&dir, 1, Records::Wide(4)).records {
            let Some(sequence) = genome.resolve(&record.sequence) else {
                panic!("{corpus}: {} names no sequence", record.id);
            };
            for (offset, base) in record.reference.bytes().enumerate() {
                let at = record.position + offset as u64;
                for transcript in genome.on(&sequence) {
                    if !complete(transcript) {
                        continue;
                    }
                    let Some(index) = coding_index(transcript, at) else {
                        continue;
                    };
                    let cached = transcript.coding.as_ref().expect("a coding model").sequence
                        [index as usize];
                    assert!(
                        cached == as_read(transcript, base),
                        "{corpus}: the reference allele of {} is not the cache's base",
                        record.id
                    );
                    inside += 1;
                }
            }
        }
    }
    assert!(
        inside > 200,
        "only {inside} bases were inside a coding sequence"
    );

    // With the yeast slice's reference file: every reference allele is the file's.
    let corpus = "R64-1-1-yeast";
    let reference = pilot(corpus).join("reference.fa.gz");
    let genome = Genome::load(&pilot_cache(corpus));
    let scratch = Scratch::new("reference");
    let fasta = reference.to_str().expect("a path");
    let (records, summary) = run_to_bytes(
        &pilot_cache(corpus),
        &scratch,
        &["--seed", "1", "--fasta", fasta],
    );
    let sequences = whole_fasta(&reference);
    let records = String::from_utf8(records).expect("text");
    assert!(records.contains("##reference_alleles=fasta\n"));
    let mut checked = 0;
    for line in records.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let sequence = genome.resolve(fields[0]).expect("a sequence of the cache");
        let position: usize = fields[1].parse().expect("a position");
        let bases = &sequences[&sequence];
        let expected = &bases[position - 1..position - 1 + fields[3].len()];
        assert!(
            fields[3].as_bytes() == expected,
            "the reference allele of {} is not the reference file's",
            fields[2]
        );
        assert!(expected.iter().all(|base| b"ACGT".contains(base)));
        checked += 1;
    }
    assert!(checked >= 20, "only {checked} records were written");

    // The file is masked outside the corpus's neighbourhoods: a masked candidate is
    // passed over, and the summary counts them.
    let summary = String::from_utf8(summary).expect("text");
    let passed_over: u64 = summary
        .lines()
        .skip(1)
        .map(|line| line.split('\t').nth(3).expect("a count"))
        .map(|count| count.parse::<u64>().expect("a number"))
        .sum();
    assert!(passed_over > 0, "no candidate was passed over");
    assert!(summary.contains("\treference masked\n"));
}

#[test]
fn the_indexed_reader_gives_the_bases_a_whole_read_gives() {
    let reference = pilot("R64-1-1-yeast").join("reference.fa.gz");
    let sequences = whole_fasta(&reference);
    let mut fasta = Fasta::open(&reference).expect("an indexed file");
    for (name, bases) in &sequences {
        assert_eq!(fasta.length(name), Some(bases.len() as u64), "{name}");
        let length = bases.len() as u64;
        for (position, count) in [
            (1, 1),
            (1, 130),
            (60, 2),
            (61, 61),
            (length / 2, 500),
            (length, 1),
            (length - 99, 100),
        ] {
            let read = fasta.bases(name, position, count).expect("bases");
            let start = (position - 1) as usize;
            assert!(
                read == bases[start..start + count as usize],
                "{name}: {count} bases from {position}"
            );
        }
        assert!(fasta.bases(name, length, 2).is_err());
        assert!(fasta.bases(name, 0, 1).is_err());
    }
    assert!(!fasta.has("no such sequence"));
    assert_eq!(fasta.failure(), None);

    // A plain copy cut off half way: a span the index does not hold is no failed read, a
    // read past the cut is one, and a run over the copy stops on it.
    let mut text = Vec::new();
    MultiGzDecoder::new(File::open(&reference).expect("the reference file"))
        .read_to_end(&mut text)
        .expect("the reference file's text");
    let scratch = Scratch::new("cut-reference");
    let cut = scratch.path().join("cut.fa");
    fs::write(&cut, &text[..text.len() / 2]).expect("a written file");
    fs::copy(
        pilot("R64-1-1-yeast").join("reference.fa.fai"),
        scratch.path().join("cut.fa.fai"),
    )
    .expect("a copied index");
    let mut fasta = Fasta::open(&cut).expect("an indexed file");
    for (name, bases) in &sequences {
        assert!(fasta.bases(name, bases.len() as u64 + 1, 1).is_err());
    }
    assert_eq!(fasta.failure(), None);
    let whole = sequences
        .iter()
        .filter(|(name, bases)| fasta.bases(name, bases.len() as u64, 1).is_ok())
        .count();
    assert!(
        0 < whole && whole < sequences.len(),
        "{whole} sequences end before the cut"
    );
    assert!(fasta
        .failure()
        .is_some_and(|failure| failure.contains("cut.fa")));

    let output = scratch.path().join("records.vcf");
    let ended = run(&[
        "--cache",
        pilot_cache("R64-1-1-yeast").to_str().expect("a path"),
        "--seed",
        "1",
        "--output",
        output.to_str().expect("a path"),
        "--fasta",
        cut.to_str().expect("a path"),
    ]);
    assert_eq!(ended.code, Some(1), "{}", ended.stderr);
    assert!(ended.stderr.starts_with("ERROR: [vep-seeded-variants] "));
    assert!(ended.stderr.contains("cut.fa"));
    assert_eq!(ended.stderr.lines().count(), 1);
    assert!(!output.exists());
}
