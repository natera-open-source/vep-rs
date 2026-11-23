// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! An Ensembl VEP cache files a gene's transcripts in every region file the gene reaches,
//! so a region file repeats transcripts an earlier one holds. The reader gives each
//! transcript once a sequence, whichever region files hold it, and a run writes for such
//! a cache what it writes for the same cache with every transcript filed once. A record
//! under the same stable id at another place of the sequence is another transcript.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{bacterium, category_of, records_of, run_to_bytes, Genome, Scratch};
use serde_json::{json, Value};
use vep_seeded_variants::cache::Cache;
use vep_seeded_variants::generate::Records;

/// The genetic code the transcripts of the made-up cache's short sequence name.
const BACTERIAL: u8 = 11;

/// How the long sequence's two region files hold its gene.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Filed {
    /// Each transcript in the region file its span lies in.
    Once,
    /// Both transcripts in both region files, as a cache files a gene that reaches both.
    InBothRegions,
}

/// The two transcripts of the made-up plasmid: bases 101 to 700 on the forward strand
/// and bases 900 to 1,499 on the reverse strand.
fn plasmid_transcripts() -> (Value, Value) {
    let shard = fs::read_to_string(bacterium().join("transcripts/pMADE/1-1000000.json"))
        .expect("the plasmid's shard");
    let shard: Value = serde_json::from_str(&shard).expect("a shard");
    let records = shard.as_array().expect("an array shard");
    assert_eq!(records.len(), 2, "the plasmid's two transcripts");
    assert_eq!(
        (&records[0]["stable_id"], &records[0]["start"]),
        (&json!("MADE_p001_t"), &json!(101))
    );
    assert_eq!(
        (&records[1]["stable_id"], &records[1]["end"]),
        (&json!("MADE_p002_t"), &json!(1499))
    );
    (records[0].clone(), records[1].clone())
}

fn add(number: &mut Value, by: u64) {
    *number = json!(number.as_u64().expect("a position") + by);
}

/// The record `by` bases further along its sequence.
fn moved(mut record: Value, by: u64) -> Value {
    add(&mut record["start"], by);
    add(&mut record["end"], by);
    for exon in record["exons"].as_array_mut().expect("exons") {
        add(&mut exon["start"], by);
        add(&mut exon["end"], by);
    }
    let pairs = &mut record["variation_effect_feature_cache"]["mapper"]["pairs"];
    for pair in pairs.as_array_mut().expect("mapper pairs") {
        add(&mut pair["to_start"], by);
        add(&mut pair["to_end"], by);
    }
    record
}

/// The record under the stable id, gene and genetic code given.
fn renamed(mut record: Value, id: &str, gene: &str, code: u8) -> Value {
    record["stable_id"] = json!(id);
    record["gene_stable_id"] = json!(gene);
    record["variation_effect_feature_cache"]["codon_table"] = json!(code);
    record
}

/// A made-up cache of two sequences in `scratch`.
///
/// `short`, 2,400 bases, holds three transcripts under the bacterial code in one region
/// file: `MADE_p001_t`, `MADE_p002_t` and `MADE_p003_t`.
///
/// `long`, 1,002,400 bases, holds one gene of two transcripts under `long_code`:
/// `MADE_p001_t` at bases 101 to 700, in the first region, and `MADE_p002_t` at bases
/// 1,000,900 to 1,001,499, in the second. `filed` says which region files hold them.
fn made_up(scratch: &Scratch, filed: Filed, long_code: u8) -> PathBuf {
    let (first, second) = plasmid_transcripts();
    let dir = scratch.path().join(format!("cache-{filed:?}-{long_code}"));
    let write = |sequence: &str, region: &str, records: Vec<Value>| {
        let sequence = dir.join("transcripts").join(sequence);
        fs::create_dir_all(&sequence).expect("a directory");
        fs::write(sequence.join(region), Value::Array(records).to_string())
            .expect("a written file");
    };

    let third = renamed(
        moved(first.clone(), 1_500),
        "MADE_p003_t",
        "MADE_p003",
        BACTERIAL,
    );
    write(
        "short",
        "1-1000000.json",
        vec![first.clone(), second.clone(), third],
    );

    let low = renamed(first, "MADE_p001_t", "MADE_long", long_code);
    let high = renamed(
        moved(second, 1_000_000),
        "MADE_p002_t",
        "MADE_long",
        long_code,
    );
    let (in_first, in_second) = match filed {
        Filed::Once => (vec![low], vec![high]),
        Filed::InBothRegions => (vec![low.clone(), high.clone()], vec![low, high]),
    };
    write("long", "1-1000000.json", in_first);
    write("long", "1000001-2000000.json", in_second);

    let info = fs::read_to_string(bacterium().join("info.json")).expect("info.json");
    let mut info: Value = serde_json::from_str(&info).expect("info.json");
    info["contig_order"] = json!(["short", "long"]);
    fs::write(dir.join("info.json"), info.to_string()).expect("a written file");
    fs::write(
        dir.join("sequence_regions.tsv"),
        "name\tlength\tcircular\nshort\t2400\t0\nlong\t1002400\t0\n",
    )
    .expect("a written file");
    dir
}

/// The stable id and the start of each transcript the reader gives for `sequence`.
fn given(cache: &Path, sequence: &str) -> Vec<(String, u64)> {
    Cache::open(cache)
        .expect("a cache")
        .transcripts(sequence)
        .expect("a sequence")
        .map(|transcript| transcript.expect("a transcript"))
        .map(|transcript| (transcript.stable_id, transcript.start))
        .collect()
}

fn ids(given: &[(&str, u64)]) -> Vec<(String, u64)> {
    given
        .iter()
        .map(|(id, start)| ((*id).to_owned(), *start))
        .collect()
}

#[test]
fn a_transcript_two_region_files_hold_is_given_once() {
    let scratch = Scratch::new("repeat-reader");
    let both = [("MADE_p001_t", 101), ("MADE_p002_t", 1_000_900)];
    for filed in [Filed::Once, Filed::InBothRegions] {
        let cache = made_up(&scratch, filed, 4);
        // Where the second region file repeats a transcript that ends in the first region
        // and one that lies in the second, neither is given again.
        assert_eq!(given(&cache, "long"), ids(&both), "{filed:?}");
        // The same stable id on another sequence is another transcript.
        assert_eq!(
            given(&cache, "short"),
            ids(&[
                ("MADE_p001_t", 101),
                ("MADE_p002_t", 900),
                ("MADE_p003_t", 1_601)
            ]),
            "{filed:?}"
        );
    }
}

#[test]
fn a_transcript_three_region_files_hold_is_given_once() {
    // A third region file holds both transcripts again.
    let scratch = Scratch::new("repeat-three");
    let cache = made_up(&scratch, Filed::InBothRegions, 4);
    let long = cache.join("transcripts/long");
    fs::copy(
        long.join("1000001-2000000.json"),
        long.join("2000001-3000000.json"),
    )
    .expect("a copied file");
    assert_eq!(
        given(&cache, "long"),
        ids(&[("MADE_p001_t", 101), ("MADE_p002_t", 1_000_900)])
    );
}

#[test]
fn a_stable_id_at_another_place_of_a_sequence_is_another_transcript() {
    // The second region file holds one more record under the stable id of the
    // transcript at bases 101 to 700, at bases 1,001,701 to 1,002,300.
    let scratch = Scratch::new("repeat-second-place");
    let cache = made_up(&scratch, Filed::Once, 4);
    let (first, _) = plasmid_transcripts();
    let elsewhere = moved(first, 1_001_600);
    let shard = cache.join("transcripts/long/1000001-2000000.json");
    let text = fs::read_to_string(&shard).expect("a shard");
    let mut records: Value = serde_json::from_str(&text).expect("a shard");
    records.as_array_mut().expect("an array").push(elsewhere);
    fs::write(&shard, records.to_string()).expect("a written file");
    assert_eq!(
        given(&cache, "long"),
        ids(&[
            ("MADE_p001_t", 101),
            ("MADE_p002_t", 1_000_900),
            ("MADE_p001_t", 1_001_701)
        ])
    );

    // A record of a category that stands outside every transcript stands outside the
    // transcript at the second place too.
    let mut outside = 0;
    for seed in 1..=20 {
        for records in [Records::Light, Records::Wide(12)] {
            for record in &records_of(&cache, seed, records).records {
                let beside = ["upstream", "downstream", "intergenic"];
                if record.sequence != "long" || !beside.contains(&category_of(&record.id)) {
                    continue;
                }
                assert!(
                    !(1_001_701..=1_002_300).contains(&record.position),
                    "seed {seed}, {records:?}: {} stands inside the transcript at {}",
                    record.id,
                    record.position
                );
                outside += 1;
            }
        }
    }
    assert!(outside >= 40, "only {outside} records beside a transcript");

    // A record that shares the start alone or the end alone with a transcript given under
    // its stable id is another transcript: the first region file holds one more at bases
    // 41 to 700 and one more at bases 101 to 760.
    let (first, _) = plasmid_transcripts();
    let (mut starts_before, mut ends_after) = (first.clone(), first);
    starts_before["start"] = json!(41);
    ends_after["end"] = json!(760);
    let shard = cache.join("transcripts/long/1-1000000.json");
    let text = fs::read_to_string(&shard).expect("a shard");
    let mut records: Value = serde_json::from_str(&text).expect("a shard");
    records
        .as_array_mut()
        .expect("an array")
        .extend([starts_before, ends_after]);
    fs::write(&shard, records.to_string()).expect("a written file");
    let places: Vec<(u64, u64)> = Cache::open(&cache)
        .expect("a cache")
        .transcripts("long")
        .expect("a sequence")
        .map(|transcript| transcript.expect("a transcript"))
        .filter(|transcript| transcript.stable_id == "MADE_p001_t")
        .map(|transcript| (transcript.start, transcript.end))
        .collect();
    assert_eq!(
        places,
        [(41, 700), (101, 700), (101, 760), (1_001_701, 1_002_300)]
    );
}

/// The cache's own names of the sequences that hold an `organelle` record.
fn organelle_sequences(records: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(records)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let mut fields = line.split('\t');
            let sequence = fields.next()?;
            let id = fields.nth(1)?;
            (category_of(id) == "organelle").then(|| sequence.to_owned())
        })
        .collect()
}

#[test]
fn a_repeated_transcript_counts_once_among_the_transcripts_of_the_cache() {
    // Three transcripts name the bacterial code and two name code 4: most transcripts
    // name the bacterial code, so the long sequence is the organelle. A transcript given
    // twice makes it three against three, and the short sequence the organelle.
    let scratch = Scratch::new("repeat-pinned");
    let cache = made_up(&scratch, Filed::InBothRegions, 4);
    assert_eq!(Genome::load(&cache).organelles(), ["long"]);
    let (records, summary) = run_to_bytes(&cache, &scratch, &["--seed", "1"]);
    assert_eq!(organelle_sequences(&records), ["long"]);

    let expected = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/expected");
    let pinned = [
        (expected.join("made-up-repeated-transcript.vcf"), records),
        (
            expected.join("made-up-repeated-transcript.summary.tsv"),
            summary,
        ),
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
            "the run on a cache that files a transcript twice does not write {}:\n{}",
            file.display(),
            String::from_utf8_lossy(&written)
        );
    }
}

#[test]
fn a_cache_that_repeats_a_transcript_writes_what_the_cache_without_the_repeat_writes() {
    // Under code 4 the long sequence names another code than most transcripts; under
    // the bacterial code the cache names one code and the shorter sequence is taken for
    // a plasmid. Either way the repeat moves no record and no summary line.
    let scratch = Scratch::new("repeat-equal");
    for long_code in [4, BACTERIAL] {
        let once = made_up(&scratch, Filed::Once, long_code);
        let repeated = made_up(&scratch, Filed::InBothRegions, long_code);
        for mode in [
            &["--seed", "1"][..],
            &["--seed", "1", "--records", "wide"][..],
            &["--seed", "2", "--sequence", "long"][..],
            &["--seed", "2", "--sequence", "short"][..],
        ] {
            assert!(
                run_to_bytes(&once, &scratch, mode) == run_to_bytes(&repeated, &scratch, mode),
                "code {long_code}, {mode:?}: the repeat changes the output"
            );
        }
    }
}

/// A reference file that masks every base of the made-up cache's two sequences, with its
/// index: each sequence on one line.
fn masked_reference(scratch: &Scratch) -> PathBuf {
    let file = scratch.path().join("masked.fa");
    let (short, long) = (2_400, 1_002_400);
    let text = format!(
        ">short\n{}\n>long\n{}\n",
        "N".repeat(short),
        "N".repeat(long)
    );
    fs::write(&file, text).expect("a written file");
    let short_at = ">short\n".len();
    let long_at = short_at + short + 1 + ">long\n".len();
    let index = format!(
        "short\t{short}\t{short_at}\t{short}\t{}\nlong\t{long}\t{long_at}\t{long}\t{}\n",
        short + 1,
        long + 1
    );
    fs::write(scratch.path().join("masked.fa.fai"), index).expect("a written file");
    file
}

/// How many candidates of `category` the summary says the reference file passed over.
fn passed_over(summary: &[u8], category: &str) -> u64 {
    String::from_utf8_lossy(summary)
        .lines()
        .find_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (fields[0] == category).then(|| fields[3].parse().expect("a count"))
        })
        .unwrap_or_else(|| panic!("no summary line for {category}"))
}

#[test]
fn a_reference_file_passes_over_a_repeated_transcript_once() {
    // Every base is masked, so every candidate is passed over and the summary counts
    // them: five transcripts with a coding sequence, three on the short sequence and two
    // on the long one, however many region files hold the long sequence's two.
    let scratch = Scratch::new("repeat-masked");
    let reference = masked_reference(&scratch);
    let options = [
        "--seed",
        "1",
        "--fasta",
        reference.to_str().expect("a path"),
    ];
    let once = run_to_bytes(&made_up(&scratch, Filed::Once, 4), &scratch, &options);
    let repeated = run_to_bytes(
        &made_up(&scratch, Filed::InBothRegions, 4),
        &scratch,
        &options,
    );
    assert_eq!(passed_over(&repeated.1, "inframe_deletion"), 5);
    assert_eq!(passed_over(&repeated.1, "code4_stop_gained"), 2);
    assert!(
        once == repeated,
        "the repeat changes what a run with a reference file writes"
    );
}

/// Rewrites one region file of the long sequence through `change`.
fn rewrite(cache: &Path, region: &str, change: impl FnOnce(&mut Vec<Value>)) {
    let shard = cache.join("transcripts/long").join(region);
    let text = fs::read_to_string(&shard).expect("a shard");
    let mut records: Value = serde_json::from_str(&text).expect("a shard");
    change(records.as_array_mut().expect("an array"));
    fs::write(&shard, records.to_string()).expect("a written file");
}

/// Gives the long sequence's two transcripts the `dbID`s 11 and 12.
fn numbered(records: &mut [Value]) {
    for record in records.iter_mut() {
        let id = if record["stable_id"] == "MADE_p001_t" {
            11
        } else {
            12
        };
        record["dbID"] = json!(id);
    }
}

/// A record without a coding sequence under the stable id, start and end of the long
/// sequence's first transcript, with the biotype `cDNA` and the `dbID` 13.
fn aligned_copy() -> Value {
    let (first, _) = plasmid_transcripts();
    let mut record = renamed(first, "MADE_p001_t", "MADE_aligned", 4);
    record["biotype"] = json!("cDNA");
    record["dbID"] = json!(13);
    let keys = record.as_object_mut().expect("a record");
    for key in [
        "cdna_coding_start",
        "cdna_coding_end",
        "translation",
        "protein",
    ] {
        keys.remove(key);
    }
    let cached = record["variation_effect_feature_cache"]
        .as_object_mut()
        .expect("the feature cache");
    for key in ["translateable_seq", "peptide"] {
        cached.remove(key);
    }
    let mapper = cached["mapper"].as_object_mut().expect("the mapper");
    for key in ["cdna_coding_start", "cdna_coding_end"] {
        mapper.remove(key);
    }
    record
}

/// The stable id, the biotype and the start of each transcript the reader gives for the
/// long sequence.
fn given_with_biotype(cache: &Path) -> Vec<(String, String, u64)> {
    Cache::open(cache)
        .expect("a cache")
        .transcripts("long")
        .expect("a sequence")
        .map(|transcript| transcript.expect("a transcript"))
        .map(|transcript| (transcript.stable_id, transcript.biotype, transcript.start))
        .collect()
}

fn with_biotype(given: &[(&str, &str, u64)]) -> Vec<(String, String, u64)> {
    given
        .iter()
        .map(|(id, biotype, start)| ((*id).to_owned(), (*biotype).to_owned(), *start))
        .collect()
}

#[test]
fn a_record_under_another_database_id_is_another_transcript() {
    // A cache built from RefSeq models can hold an aligned cDNA under the stable id,
    // start and end of a model, under a `dbID` of its own. Both region files hold the
    // gene's two transcripts, and the second holds such a record beside them.
    let scratch = Scratch::new("repeat-database-id");
    let cache = made_up(&scratch, Filed::InBothRegions, 4);
    rewrite(&cache, "1-1000000.json", |records| numbered(records));
    rewrite(&cache, "1000001-2000000.json", |records| {
        numbered(records);
        records.push(aligned_copy());
    });
    assert_eq!(
        given_with_biotype(&cache),
        with_biotype(&[
            ("MADE_p001_t", "protein_coding", 101),
            ("MADE_p002_t", "protein_coding", 1_000_900),
            ("MADE_p001_t", "cDNA", 101),
        ])
    );
    let coding: Vec<bool> = Cache::open(&cache)
        .expect("a cache")
        .transcripts("long")
        .expect("a sequence")
        .map(|transcript| transcript.expect("a transcript").coding.is_some())
        .collect();
    assert_eq!(coding, [true, true, false]);

    // The first region file holds the record too, ahead of the model: both are given
    // there, the lower `dbID` first, and neither again from the second region file.
    let scratch = Scratch::new("repeat-database-id-one-file");
    let cache = made_up(&scratch, Filed::InBothRegions, 4);
    rewrite(&cache, "1-1000000.json", |records| {
        numbered(records);
        records.insert(0, aligned_copy());
    });
    rewrite(&cache, "1000001-2000000.json", |records| {
        numbered(records);
        records.push(aligned_copy());
    });
    assert_eq!(
        given_with_biotype(&cache),
        with_biotype(&[
            ("MADE_p001_t", "protein_coding", 101),
            ("MADE_p001_t", "cDNA", 101),
            ("MADE_p002_t", "protein_coding", 1_000_900),
        ])
    );
}
