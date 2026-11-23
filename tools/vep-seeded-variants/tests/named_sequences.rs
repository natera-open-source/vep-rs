// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--sequence` names where the records stand, and the rule for organelle and plasmid
//! sequences reads the whole cache whichever sequences are named: a named sequence gets
//! an `organelle` record only where the rule takes it for one among all the cache's
//! sequences.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{bacterium, category_of, copy_tree, run_to_bytes, summary_line, Genome, Scratch};
use vep_seeded_variants::cache::Cache;
use vep_seeded_variants::generate::{generate, Generated, Options, Records};

/// The genetic code every transcript of the made-up bacterium names.
const BACTERIAL: u8 = 11;

/// Twelve scaffolds: with the chromosome and the plasmid, a cache of fourteen sequences.
const SCAFFOLDS: [&str; 12] = [
    "scaffold_01",
    "scaffold_02",
    "scaffold_03",
    "scaffold_04",
    "scaffold_05",
    "scaffold_06",
    "scaffold_07",
    "scaffold_08",
    "scaffold_09",
    "scaffold_10",
    "scaffold_11",
    "scaffold_12",
];

/// A made-up cache in `scratch`: the bacterium's chromosome and plasmid, then one
/// sequence for each of `more`, which holds the plasmid's two genes under the genetic
/// code named beside it. With `lengths` the side file gives every sequence's length,
/// 9,000 bases for the chromosome and 2,400 for each other sequence; without, the cache
/// holds no side file, and the highest transcript end is 9,000 on the chromosome and
/// 1,499 on each other sequence.
fn made_up(scratch: &Scratch, more: &[(&str, u8)], lengths: bool) -> PathBuf {
    let copy = scratch.path().join("cache");
    copy_tree(&bacterium(), &copy, false);
    let plasmid = fs::read_to_string(copy.join("transcripts/pMADE/1-1000000.json"))
        .expect("the plasmid's shard");
    let named = format!(r#""codon_table": {BACTERIAL}"#);
    assert_eq!(
        plasmid.matches(&named).count(),
        2,
        "the plasmid's two genes"
    );
    let mut order = String::from(r#""pMADE""#);
    let mut side = String::from("name\tlength\tcircular\nChromosome\t9000\t1\npMADE\t2400\t1\n");
    for (name, code) in more {
        let dir = copy.join("transcripts").join(name);
        fs::create_dir(&dir).expect("a directory");
        let shard = plasmid.replace(&named, &format!(r#""codon_table": {code}"#));
        fs::write(dir.join("1-1000000.json"), shard).expect("a written file");
        order.push_str(&format!(r#", "{name}""#));
        side.push_str(&format!("{name}\t2400\t0\n"));
    }
    let info = fs::read_to_string(copy.join("info.json")).expect("info.json");
    assert_eq!(info.matches(r#""pMADE""#).count(), 1);
    fs::write(copy.join("info.json"), info.replace(r#""pMADE""#, &order)).expect("a written file");
    let side_file = copy.join("sequence_regions.tsv");
    if lengths {
        fs::write(&side_file, side).expect("a written file");
    } else {
        fs::remove_file(&side_file).expect("a removed file");
    }
    copy
}

/// The fourteen-sequence cache, every transcript under the bacterial code but those of
/// `other`, which name the code beside it.
fn fourteen(scratch: &Scratch, other: Option<(&str, u8)>) -> PathBuf {
    let more: Vec<(&str, u8)> = SCAFFOLDS
        .iter()
        .map(|name| match other {
            Some((odd, code)) if odd == *name => (*name, code),
            _ => (*name, BACTERIAL),
        })
        .collect();
    made_up(scratch, &more, true)
}

/// A light run for seed 1 with records on `names`, on every sequence where none is named.
fn light(cache: &Path, names: &[&str]) -> Generated {
    let cache = Cache::open(cache).expect("a cache");
    let options = Options {
        seed: 1,
        records: Records::Light,
        sequences: names.iter().map(|name| (*name).to_owned()).collect(),
    };
    generate(&cache, None, &options).expect("a run")
}

/// The cache's own names of the sequences that hold an `organelle` record, in output
/// order.
fn organelle_sequences(cache: &Path, generated: &Generated) -> Vec<String> {
    let genome = Genome::load(cache);
    generated
        .records
        .iter()
        .filter(|record| category_of(&record.id) == "organelle")
        .map(|record| genome.resolve(&record.sequence).expect("a sequence"))
        .collect()
}

/// Holds one run to the rule. Its `organelle` records stand on `expected`, which is
/// also what the rule read a second way over the whole cache gives among `names`, three
/// sequences at most; every record stands on a named sequence; and the summary counts
/// the `organelle` records or says why there is none.
fn holds_the_rule(cache: &Path, names: &[&str], expected: &[&str]) {
    let genome = Genome::load(cache);
    let by_rule: Vec<String> = genome
        .organelles()
        .into_iter()
        .filter(|sequence| names.is_empty() || names.contains(&sequence.as_str()))
        .take(3)
        .collect();
    assert_eq!(by_rule, expected, "the rule, read a second way: {names:?}");
    let generated = light(cache, names);
    assert_eq!(
        organelle_sequences(cache, &generated),
        expected,
        "{names:?}"
    );
    for record in &generated.records {
        let sequence = genome.resolve(&record.sequence).expect("a sequence");
        assert!(
            names.is_empty() || names.contains(&sequence.as_str()),
            "{names:?}: {} stands on {sequence}",
            record.id
        );
    }
    let line = if expected.is_empty() {
        (0, 0, "no organelle or plasmid sequence")
    } else {
        (expected.len(), expected.len(), "-")
    };
    assert_eq!(summary_line(&generated, "organelle"), line, "{names:?}");
}

#[test]
fn two_named_sequences_of_a_cache_of_many_get_no_organelle_record() {
    // Fourteen sequences under one genetic code: more than the twelve up to which every
    // sequence but the longest is taken for a plasmid, so the cache has no organelle or
    // plasmid sequence, whether two of its sequences are named or none is.
    let scratch = Scratch::new("many");
    let cache = fourteen(&scratch, None);
    let (records, summary) = run_to_bytes(
        &cache,
        &scratch,
        &[
            "--seed",
            "1",
            "--sequence",
            "Chromosome",
            "--sequence",
            "pMADE",
        ],
    );
    let expected = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/expected");
    let pinned = [
        (expected.join("made-up-many-sequences.vcf"), records),
        (expected.join("made-up-many-sequences.summary.tsv"), summary),
    ];
    let update = std::env::var_os("VEP_SEEDED_VARIANTS_UPDATE").is_some();
    for (file, written) in pinned {
        if update {
            fs::write(&file, &written).expect("a written file");
            continue;
        }
        let expected = fs::read(&file).unwrap_or_default();
        assert!(
            written == expected,
            "the run on two named sequences of fourteen does not write {}:\n{}",
            file.display(),
            String::from_utf8_lossy(&written)
        );
    }
    holds_the_rule(&cache, &["Chromosome", "pMADE"], &[]);
    holds_the_rule(&cache, &[], &[]);
}

#[test]
fn a_named_plasmid_of_a_cache_of_few_sequences_gets_its_organelle_record() {
    // Two sequences under one code: the plasmid, the shorter one, is the cache's
    // plasmid whether it is named alone, with the chromosome or not at all.
    let cache = bacterium();
    holds_the_rule(&cache, &[], &["pMADE"]);
    holds_the_rule(&cache, &["pMADE"], &["pMADE"]);
    holds_the_rule(&cache, &["Chromosome", "pMADE"], &["pMADE"]);
    holds_the_rule(&cache, &["Chromosome"], &[]);
}

#[test]
fn the_longest_sequence_is_the_longest_of_the_cache() {
    // Chromosome, plasmid and a second plasmid under one code. The chromosome is the
    // longest, by the side file's lengths and by the highest transcript end, so each
    // plasmid is taken for one, the chromosome named or not.
    for lengths in [true, false] {
        let scratch = Scratch::new("three");
        let cache = made_up(&scratch, &[("pSECOND", BACTERIAL)], lengths);
        holds_the_rule(&cache, &[], &["pMADE", "pSECOND"]);
        holds_the_rule(&cache, &["pMADE", "pSECOND"], &["pMADE", "pSECOND"]);
        holds_the_rule(&cache, &["pSECOND"], &["pSECOND"]);
        holds_the_rule(&cache, &["Chromosome", "pSECOND"], &["pSECOND"]);
        holds_the_rule(&cache, &["Chromosome"], &[]);
    }
}

#[test]
fn another_code_is_another_than_most_transcripts_of_the_cache() {
    // One scaffold of the fourteen names code 4. Named with one other scaffold, each
    // holds two transcripts, and most transcripts of the cache name the bacterial code:
    // the scaffold under code 4 is the organelle.
    let scratch = Scratch::new("other-code");
    let cache = fourteen(&scratch, Some(("scaffold_04", 4)));
    holds_the_rule(&cache, &[], &["scaffold_04"]);
    holds_the_rule(&cache, &["scaffold_04", "scaffold_05"], &["scaffold_04"]);
    holds_the_rule(&cache, &["scaffold_05", "scaffold_06"], &[]);
    holds_the_rule(&cache, &["scaffold_04"], &["scaffold_04"]);

    // Three sequences, the third under code 4. A sequence of the cache names another
    // code, so no sequence is taken for a plasmid by its length, named or not.
    let scratch = Scratch::new("other-code-few");
    let cache = made_up(&scratch, &[("mito", 4)], true);
    holds_the_rule(&cache, &[], &["mito"]);
    holds_the_rule(&cache, &["Chromosome", "pMADE"], &[]);
    holds_the_rule(&cache, &["pMADE", "mito"], &["mito"]);
}

#[test]
fn three_named_organelles_at_most_get_a_light_record() {
    // Six sequences under one code: the five plasmids are the cache's plasmids, and a
    // light run writes a record on the first three of those it places records on.
    let scratch = Scratch::new("cap");
    let plasmids = [
        ("pSECOND", BACTERIAL),
        ("pTHIRD", BACTERIAL),
        ("pFOURTH", BACTERIAL),
        ("pFIFTH", BACTERIAL),
    ];
    let cache = made_up(&scratch, &plasmids, true);
    holds_the_rule(&cache, &[], &["pMADE", "pSECOND", "pTHIRD"]);
    holds_the_rule(
        &cache,
        &["pTHIRD", "pFOURTH", "pFIFTH"],
        &["pTHIRD", "pFOURTH", "pFIFTH"],
    );
    holds_the_rule(
        &cache,
        &["pSECOND", "pTHIRD", "pFOURTH", "pFIFTH"],
        &["pSECOND", "pTHIRD", "pFOURTH"],
    );
}

#[test]
fn naming_every_sequence_gives_the_records_of_naming_none() {
    let scratch = Scratch::new("every");
    let many = fourteen(&scratch, Some(("scaffold_04", 4)));
    for dir in [bacterium(), many] {
        let cache = Cache::open(&dir).expect("a cache");
        let mut every: Vec<String> = cache.sequences().to_vec();
        every.reverse();
        for records in [Records::Light, Records::Wide(12)] {
            let run = |sequences: Vec<String>| {
                let options = Options {
                    seed: 1,
                    records,
                    sequences,
                };
                generate(&cache, None, &options).expect("a run")
            };
            assert!(
                run(every.clone()) == run(Vec::new()),
                "{}: {records:?}",
                dir.display()
            );
        }
    }
}
