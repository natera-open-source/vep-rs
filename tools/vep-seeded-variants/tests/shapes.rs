// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Caches that lack a situation write fewer records and say why: a bacterial shape, a
//! sequence without transcripts, a cache without synonyms, a transcript stored across an
//! origin; and the wide mode writes each category several times over.

mod common;

use std::fs;
use std::path::PathBuf;

use common::{bacterium, category_of, copy_tree, pilot_cache, records_of, summary_line, Scratch};
use vep_seeded_variants::cache::Cache;
use vep_seeded_variants::generate::{generate, Generated, Options, Records};
use vep_seeded_variants::vcf::Record;

fn records_in<'a>(generated: &'a Generated, category: &str) -> Vec<&'a str> {
    generated
        .records
        .iter()
        .filter(|record| category_of(&record.id) == category)
        .map(|record| record.id.as_str())
        .collect()
}

#[test]
fn a_bacterial_shape_writes_fewer_records_and_says_why() {
    let generated = records_of(&bacterium(), 1, Records::Light);
    for category in ["splice_donor", "splice_acceptor", "intron"] {
        assert!(records_in(&generated, category).is_empty(), "{category}");
        assert_eq!(
            summary_line(&generated, category),
            (1, 0, "no multi-exon transcript")
        );
    }
    for category in ["five_prime_utr", "three_prime_utr"] {
        assert_eq!(
            summary_line(&generated, category),
            (1, 0, "no untranslated region")
        );
    }
    assert_eq!(
        summary_line(&generated, "intergenic"),
        (1, 0, "no gap between genes")
    );
    for category in ["first_base", "last_base", "origin", "organelle", "synonym"] {
        assert_eq!(records_in(&generated, category).len(), 1, "{category}");
        assert_eq!(summary_line(&generated, category), (1, 1, "-"));
    }
    for kind in ["changed_codon", "start_codon", "stop_gained", "stop_lost"] {
        let category = format!("code11_{kind}");
        assert_eq!(records_in(&generated, &category).len(), 1, "{category}");
    }
    let cache = Cache::open(&bacterium()).expect("a cache");
    let origin = generated
        .records
        .iter()
        .find(|record| record.id == "origin_1")
        .expect("the origin record");
    let length = cache
        .sequence_region(&origin.sequence)
        .expect("the sequence's row")
        .length;
    assert_eq!(origin.position, length - 1);
    assert_eq!((origin.reference.len(), origin.alt.len()), (4, 1));
    let last = generated
        .records
        .iter()
        .find(|record| record.id == "last_base_1")
        .expect("the last-base record");
    let length = cache
        .sequence_region(&last.sequence)
        .expect("the sequence's row")
        .length;
    assert_eq!(last.position, length);

    // The same cache without its side file: no length is known.
    let scratch = Scratch::new("no-side-file");
    let copy = scratch.path().join("cache");
    copy_tree(&bacterium(), &copy, false);
    fs::remove_file(copy.join("sequence_regions.tsv")).expect("a removed file");
    let without = records_of(&copy, 1, Records::Light);
    for category in ["last_base", "origin"] {
        assert!(records_in(&without, category).is_empty(), "{category}");
        assert_eq!(
            summary_line(&without, category),
            (1, 0, "no sequence length")
        );
    }
    assert_eq!(records_in(&without, "first_base").len(), 1);
    assert!(without.records.len() + 2 == generated.records.len());
}

#[test]
fn a_linear_sequence_gets_no_origin_record() {
    let scratch = Scratch::new("linear");
    let copy = scratch.path().join("cache");
    copy_tree(&bacterium(), &copy, false);
    fs::write(
        copy.join("sequence_regions.tsv"),
        "name\tlength\tcircular\nChromosome\t9000\t0\npMADE\t2400\t-\n",
    )
    .expect("a written file");
    let generated = records_of(&copy, 1, Records::Light);
    assert_eq!(
        summary_line(&generated, "origin"),
        (1, 0, "no circular sequence")
    );
    assert_eq!(summary_line(&generated, "last_base"), (1, 1, "-"));
}

#[test]
fn a_sequence_without_transcripts_and_a_cache_without_synonyms_are_handled() {
    let scratch = Scratch::new("bare");
    let copy = scratch.path().join("cache");
    copy_tree(&bacterium(), &copy, false);
    fs::remove_file(copy.join("chr_synonyms.txt")).expect("a removed file");
    fs::remove_file(copy.join("sequence_regions.tsv")).expect("a removed file");
    let info = fs::read_to_string(copy.join("info.json")).expect("info.json");
    let info = info.replace(r#""pMADE""#, r#""pMADE", "empty", "absent""#);
    assert!(info.contains("absent"));
    fs::write(copy.join("info.json"), info).expect("a written file");
    fs::create_dir(copy.join("transcripts/empty")).expect("a directory");

    // A name of `contig_order` without a directory is no sequence of the cache.
    let cache = Cache::open(&copy).expect("a cache");
    assert_eq!(cache.sequences(), ["Chromosome", "pMADE", "empty"]);
    assert!(cache
        .transcripts("empty")
        .expect("a sequence")
        .next()
        .is_none());

    let generated = records_of(&copy, 1, Records::Light);
    assert!(records_in(&generated, "synonym").is_empty());
    assert_eq!(summary_line(&generated, "synonym"), (1, 0, "no synonym"));
    assert!(generated
        .records
        .iter()
        .all(|record| record.sequence != "absent"));
    // The last sequence holds no transcript and still gets its record, at its first base.
    let last: Vec<_> = generated
        .records
        .iter()
        .filter(|record| category_of(&record.id) == "last_sequence")
        .collect();
    assert_eq!(last.len(), 1);
    assert_eq!((last[0].sequence.as_str(), last[0].position), ("empty", 1));
    // Three sequences that name one code: all but the longest are taken for plasmids.
    let organelles: Vec<&str> = generated
        .records
        .iter()
        .filter(|record| category_of(&record.id) == "organelle")
        .map(|record| record.sequence.as_str())
        .collect();
    assert_eq!(organelles, ["pMADE", "empty"]);

    // One sequence alone: no last sequence and no organelle.
    let options = Options {
        seed: 1,
        records: Records::Light,
        sequences: vec!["Chromosome".to_owned()],
    };
    let alone = generate(&cache, None, &options).expect("a run");
    assert!(alone
        .records
        .iter()
        .all(|record| record.sequence == "Chromosome"));
    assert_eq!(
        summary_line(&alone, "last_sequence"),
        (1, 0, "one sequence")
    );
    assert_eq!(
        summary_line(&alone, "organelle"),
        (0, 0, "no organelle or plasmid sequence")
    );
}

#[test]
fn a_transcript_stored_across_the_origin_keeps_records_off_its_low_side() {
    let scratch = Scratch::new("across");
    let copy = scratch.path().join("cache");
    copy_tree(&bacterium(), &copy, false);
    let shard = copy.join("transcripts/pMADE/1-1000000.json");
    let text = fs::read_to_string(&shard).expect("the plasmid's shard");
    // The plasmid's first gene, stored from 101 to 700, becomes one that starts 199
    // bases before the origin.
    let across = text.replacen(r#""start": 101, "strand""#, r#""start": -199, "strand""#, 1);
    assert_ne!(
        across, text,
        "the shard's first gene is not stored as expected"
    );
    fs::write(&shard, across).expect("a written file");
    for seed in 1..=10 {
        let generated = records_of(&copy, seed, Records::Light);
        for record in &generated.records {
            let on_plasmid = record.sequence == "pMADE" || record.sequence == "MADE0002.1";
            assert!(
                !on_plasmid || record.position > 700 + 5000,
                "seed {seed}: {} stands at {}:{}",
                record.id,
                record.sequence,
                record.position
            );
        }
        assert!(generated.records.len() >= 10, "seed {seed}");
        assert_eq!(
            summary_line(&generated, "last_sequence"),
            (1, 0, "no place on the sequence")
        );
    }
}

/// A copy of the made-up bacterium whose plasmid is 30,000 bases long and holds a gene
/// stored across its origin as a dumped cache stores one. The plasmid's first gene, from
/// 101 to 700, becomes one that starts 199 bases before the origin: the transcript and
/// its exon start at -199 and, with `pair_moved`, its mapper pair starts at the gene's
/// own first base, 30,000 - 199 = 29,801. The second gene, from 900 to 1,499, is moved
/// to 24,500..25,099. With `side_file` the side file gives the plasmid's length; without,
/// the copy holds no side file.
fn plasmid_gene_across_the_origin(scratch: &Scratch, pair_moved: bool, side_file: bool) -> PathBuf {
    let copy = scratch.path().join("cache");
    copy_tree(&bacterium(), &copy, false);
    let shard = copy.join("transcripts/pMADE/1-1000000.json");
    let mut text = fs::read_to_string(&shard).expect("the plasmid's shard");
    let mut replace = |from: &str, to: &str, times: usize| {
        assert_eq!(text.matches(from).count(), times, "{from}");
        text = text.replace(from, to);
    };
    replace(r#""start": 101"#, r#""start": -199"#, 2);
    if pair_moved {
        replace(r#""to_start": 101"#, r#""to_start": 29801"#, 1);
    }
    replace(r#""start": 900"#, r#""start": 24500"#, 2);
    replace(r#""end": 1499"#, r#""end": 25099"#, 2);
    replace(r#""to_start": 900"#, r#""to_start": 24500"#, 1);
    replace(r#""to_end": 1499"#, r#""to_end": 25099"#, 1);
    fs::write(&shard, text).expect("a written file");
    let side = copy.join("sequence_regions.tsv");
    if side_file {
        fs::write(
            &side,
            "name\tlength\tcircular\nChromosome\t9000\t1\npMADE\t30000\t1\n",
        )
        .expect("a written file");
    } else {
        fs::remove_file(&side).expect("a removed file");
    }
    copy
}

/// The records on the made-up plasmid, under either of its names.
fn plasmid_records(generated: &Generated) -> Vec<&Record> {
    generated
        .records
        .iter()
        .filter(|record| record.sequence == "pMADE" || record.sequence == "MADE0002.1")
        .collect()
}

#[test]
fn a_transcript_stored_across_the_origin_keeps_records_off_its_high_side() {
    // Among the first bases the plasmid takes no record up to 5,000 bases past the
    // gene's end, 700. Among the last bases the gene begins at 29,801; a record from
    // 24,801 on stands within 5,000 bases of it, and an insertion after base 24,800
    // stands at 24,801. So no reference allele holds a base outside 5,701..24,799.
    const FIRST_FREE: u64 = 5_701;
    const LAST_FREE: u64 = 24_799;
    // The gene's first base is read from its mapper pair, from the side file's length,
    // or from both.
    for (pair_moved, side_file) in [(true, true), (true, false), (false, true)] {
        let scratch = Scratch::new("across-high");
        let cache = plasmid_gene_across_the_origin(&scratch, pair_moved, side_file);
        let (mut in_second_gene, mut before_it) = (0, 0);
        for seed in 1..=5 {
            let generated = records_of(&cache, seed, Records::Wide(40));
            for record in plasmid_records(&generated) {
                let last = record.position + record.reference.len() as u64 - 1;
                assert!(
                    FIRST_FREE <= record.position && last <= LAST_FREE,
                    "pair {pair_moved}, side file {side_file}, seed {seed}: {} holds {}:{}-{last}",
                    record.id,
                    record.sequence,
                    record.position
                );
                if record.position >= 24_500 {
                    in_second_gene += 1;
                } else {
                    before_it += 1;
                }
            }
        }
        // The free stretch still takes records, the part of the second gene below the
        // stretch kept clear included.
        assert!(
            in_second_gene > 0 && before_it > 0,
            "pair {pair_moved}, side file {side_file}: {in_second_gene} and {before_it} records"
        );
    }
}

#[test]
fn a_transcript_across_the_origin_whose_first_base_is_unknown_keeps_its_sequence_clear() {
    // No mapper pair starts above the gene's end and no side file gives a length:
    // nothing says where among the last bases the gene begins.
    let scratch = Scratch::new("across-unknown");
    let cache = plasmid_gene_across_the_origin(&scratch, false, false);
    for seed in 1..=5 {
        let generated = records_of(&cache, seed, Records::Wide(40));
        let on_plasmid: Vec<&str> = plasmid_records(&generated)
            .iter()
            .map(|record| record.id.as_str())
            .collect();
        assert!(on_plasmid.is_empty(), "seed {seed}: {on_plasmid:?}");
        assert!(generated.records.len() >= 10, "seed {seed}");
    }
}

#[test]
fn the_ends_of_a_sequence_kept_clear_at_both_ends_take_no_record_and_say_so() {
    let scratch = Scratch::new("across-ends");
    let copy = plasmid_gene_across_the_origin(&scratch, true, true);
    let cache = Cache::open(&copy).expect("a cache");
    let options = Options {
        seed: 1,
        records: Records::Light,
        sequences: vec!["pMADE".to_owned()],
    };
    let generated = generate(&cache, None, &options).expect("a run");
    // The side file gives the plasmid's length and its circular flag, so neither is
    // the reason.
    for category in ["first_base", "last_base", "origin"] {
        assert!(records_in(&generated, category).is_empty(), "{category}");
        assert_eq!(
            summary_line(&generated, category),
            (1, 0, "no place on the sequence"),
            "{category}"
        );
    }
}

#[test]
fn wide_mode_writes_each_category_twelve_times_where_it_can() {
    let cache = pilot_cache("GRCm39-mouse");
    let wide = records_of(&cache, 1, Records::Wide(12));
    for category in [
        "code1_start_codon",
        "code1_stop_gained",
        "code1_stop_lost",
        "first_sequence",
        "missense_forward",
        "missense_reverse",
        "synonymous_forward",
        "synonymous_reverse",
        "splice_donor",
        "splice_acceptor",
        "five_prime_utr",
        "three_prime_utr",
        "intron",
        "upstream",
        "downstream",
        "frameshift_insertion",
        "frameshift_deletion",
        "inframe_deletion",
        "first_base",
        "organelle",
        "synonym",
    ] {
        assert_eq!(records_in(&wide, category).len(), 12, "{category}");
        assert_eq!(summary_line(&wide, category), (12, 12, "-"), "{category}");
    }
    for line in &wide.summary {
        assert!(line.written <= line.wanted, "{}", line.category);
        assert_eq!(
            line.reason == "-",
            line.written == line.wanted && line.wanted > 0
        );
        assert_eq!(records_in(&wide, &line.category).len(), line.written);
    }
    assert!(wide.records.len() > 250, "{} records", wide.records.len());
    // The records of a category are spread over sequences where several serve it.
    let mut sequences: Vec<&str> = wide
        .records
        .iter()
        .filter(|record| category_of(&record.id) == "first_base")
        .map(|record| record.sequence.as_str())
        .collect();
    sequences.dedup();
    assert_eq!(sequences.len(), 12);

    let five = records_of(&cache, 1, Records::Wide(5));
    assert_eq!(records_in(&five, "missense_forward").len(), 5);
    let light = records_of(&cache, 1, Records::Light);
    assert!(light.records.len() < 40);
    assert!(light
        .summary
        .iter()
        .all(|line| line.written <= 1 && line.wanted <= 1));

    // A count above every candidate writes them all, up to the largest count there is.
    let many = records_of(&bacterium(), 1, Records::Wide(100_000));
    let most = records_of(&bacterium(), 1, Records::Wide(usize::MAX));
    assert!(many.records.len() > 100, "{} records", many.records.len());
    assert_eq!(most.records.len(), many.records.len());
}
