// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests share: the pilot cache slices and the made-up caches, a runner of the
//! built binary, scratch directories, and a reading of a cache's geometry written here,
//! apart from the generator's own, so that a record's place is checked a second way.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use vep_core::synonyms::ChromosomeSynonyms;
use vep_seeded_variants::cache::{Cache, Strand, Transcript};
use vep_seeded_variants::generate::{generate, Generated, Options, Records};

/// The four pilot corpora, each with the genetic code its organelle names.
pub const PILOTS: [(&str, u8); 4] = [
    ("GRCm39-mouse", 2),
    ("BDGP6.54-fly", 5),
    ("TAIR10-arabidopsis", 11),
    ("R64-1-1-yeast", 3),
];

/// The committed corpus directory of a pilot.
pub fn pilot(corpus: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(corpus)
}

/// The committed cache slice of a pilot.
pub fn pilot_cache(corpus: &str) -> PathBuf {
    pilot(corpus).join("json_cache")
}

/// The made-up bacterial cache: a circular chromosome and a circular plasmid.
pub fn bacterium() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/bacterium")
}

/// The made-up cache of one sequence under genetic code 4: one coding sequence that
/// holds no codon the code reads unlike the standard code.
pub fn mitochondrion() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mitochondrion")
}

/// A directory of this test's own, removed when the value is dropped.
pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new(label: &str) -> Scratch {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "vep-seeded-variants-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).expect("a scratch directory");
        Scratch(dir)
    }

    pub fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A copy in `scratch` of the made-up cache under genetic code 4 whose coding sequence
/// names `code` instead.
pub fn mitochondrion_under(scratch: &Scratch, code: u8) -> PathBuf {
    let copy = scratch.path().join(format!("under-code-{code}"));
    copy_tree(&mitochondrion(), &copy, false);
    let shard = copy.join("transcripts/MT_MADE/1-1000000.json");
    let text = fs::read_to_string(&shard).expect("the shard");
    let named = r#""codon_table": 4"#;
    assert_eq!(text.matches(named).count(), 1, "the one coding sequence");
    let renamed = text.replace(named, &format!(r#""codon_table": {code}"#));
    fs::write(&shard, renamed).expect("a written file");
    copy
}

/// Copies a directory tree. With `reversed`, the entries of every directory are created
/// from the last name to the first.
pub fn copy_tree(from: &Path, to: &Path, reversed: bool) {
    fs::create_dir_all(to).expect("a directory");
    let mut entries: Vec<PathBuf> = fs::read_dir(from)
        .expect("a readable directory")
        .map(|entry| entry.expect("a directory entry").path())
        .collect();
    entries.sort();
    if reversed {
        entries.reverse();
    }
    for entry in entries {
        let target = to.join(entry.file_name().expect("a file name"));
        if entry.is_dir() {
            copy_tree(&entry, &target, reversed);
        } else {
            fs::copy(&entry, &target).expect("a copied file");
        }
    }
}

/// How a run of the built binary ended.
pub struct Ended {
    pub code: Option<i32>,
    pub stdout: String,
    pub stderr: String,
}

/// Runs the built binary.
pub fn run(args: &[&str]) -> Ended {
    let output = Command::new(env!("CARGO_BIN_EXE_vep-seeded-variants"))
        .args(args)
        .output()
        .expect("the built binary runs");
    Ended {
        code: output.status.code(),
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

/// Runs the binary on a cache and gives the bytes of its records and of its summary.
pub fn run_to_bytes(cache: &Path, scratch: &Scratch, more: &[&str]) -> (Vec<u8>, Vec<u8>) {
    let output = scratch.path().join("records.vcf");
    let summary = scratch.path().join("summary.tsv");
    let mut args = vec![
        "--cache",
        cache.to_str().expect("a path"),
        "--output",
        output.to_str().expect("a path"),
        "--summary",
        summary.to_str().expect("a path"),
    ];
    args.extend_from_slice(more);
    let ended = run(&args);
    assert_eq!(ended.code, Some(0), "{}", ended.stderr);
    (
        fs::read(&output).expect("the records"),
        fs::read(&summary).expect("the summary"),
    )
}

/// The records and the summary of a cache through the library.
pub fn records_of(cache: &Path, seed: u64, records: Records) -> Generated {
    let cache = Cache::open(cache).expect("a cache");
    let options = Options {
        seed,
        records,
        sequences: Vec::new(),
    };
    generate(&cache, None, &options).expect("a run")
}

/// The category a record's ID names: the ID without its ordinal.
pub fn category_of(id: &str) -> &str {
    id.rsplit_once('_').map_or(id, |(category, _)| category)
}

/// The summary line of a category.
pub fn summary_line<'a>(generated: &'a Generated, category: &str) -> (usize, usize, &'a str) {
    let line = generated
        .summary
        .iter()
        .find(|line| line.category == category)
        .unwrap_or_else(|| panic!("no summary line for {category}"));
    (line.wanted, line.written, line.reason.as_str())
}

/// A cache with every transcript in memory, for the checks.
pub struct Genome {
    pub cache: Cache,
    pub transcripts: BTreeMap<String, Vec<Transcript>>,
}

impl Genome {
    pub fn load(dir: &Path) -> Genome {
        let cache = Cache::open(dir).expect("a cache");
        let mut transcripts = BTreeMap::new();
        for sequence in cache.sequences() {
            let list: Vec<Transcript> = cache
                .transcripts(sequence)
                .expect("a sequence")
                .map(|transcript| transcript.expect("a transcript"))
                .collect();
            transcripts.insert(sequence.clone(), list);
        }
        Genome { cache, transcripts }
    }

    pub fn holds(&self, name: &str) -> bool {
        self.cache.sequences().iter().any(|held| held == name)
    }

    /// The cache's synonyms, as Ensembl VEP links them.
    pub fn synonyms(&self) -> ChromosomeSynonyms {
        let mut table = ChromosomeSynonyms::new();
        for (first, later) in self.cache.synonym_pairs() {
            table.add_synonym(first, later);
        }
        table
    }

    /// The cache's own name of the sequence a record names.
    pub fn resolve(&self, name: &str) -> Option<String> {
        let resolved = self
            .synonyms()
            .source_chr_name(name, |held| self.holds(held));
        self.holds(&resolved).then_some(resolved)
    }

    pub fn on(&self, sequence: &str) -> &[Transcript] {
        self.transcripts.get(sequence).map_or(&[], Vec::as_slice)
    }

    pub fn length(&self, sequence: &str) -> Option<u64> {
        self.cache
            .sequence_region(sequence)
            .map(|region| region.length)
    }

    /// The organelle and plasmid sequences by the rule the README states: a sequence
    /// one of whose transcripts names another code than most transcripts; where there
    /// is none and the cache holds 2 to 12 sequences, every sequence but the longest.
    pub fn organelles(&self) -> Vec<String> {
        let mut counts: BTreeMap<u8, usize> = BTreeMap::new();
        for transcript in self.transcripts.values().flatten() {
            *counts.entry(transcript.code()).or_default() += 1;
        }
        let mut most: Option<(u8, usize)> = None;
        for (code, count) in &counts {
            if most.is_none_or(|(_, top)| *count > top) {
                most = Some((*code, *count));
            }
        }
        let sequences = self.cache.sequences();
        let by_code: Vec<String> = sequences
            .iter()
            .filter(|sequence| {
                self.on(sequence)
                    .iter()
                    .any(|transcript| Some(transcript.code()) != most.map(|(code, _)| code))
            })
            .cloned()
            .collect();
        if !by_code.is_empty() || !(2..=12).contains(&sequences.len()) {
            return by_code;
        }
        let every_length = sequences.iter().all(|name| self.length(name).is_some());
        let size = |name: &String| {
            if every_length {
                self.length(name).unwrap_or(0)
            } else {
                self.on(name).iter().map(|t| t.end).max().unwrap_or(0)
            }
        };
        let mut longest = &sequences[0];
        for name in sequences {
            if size(name) > size(longest) {
                longest = name;
            }
        }
        sequences
            .iter()
            .filter(|name| *name != longest)
            .cloned()
            .collect()
    }
}

/// The transcript's exons from the lowest to the highest.
pub fn exons_low_to_high(transcript: &Transcript) -> Vec<(u64, u64)> {
    let mut exons: Vec<(u64, u64)> = transcript
        .exons
        .iter()
        .map(|exon| (exon.start, exon.end))
        .collect();
    exons.sort();
    exons
}

/// The introns from the lowest to the highest: the bases between neighbouring exons.
pub fn introns(transcript: &Transcript) -> Vec<(u64, u64)> {
    exons_low_to_high(transcript)
        .windows(2)
        .filter(|pair| pair[0].1 + 1 < pair[1].0)
        .map(|pair| (pair[0].1 + 1, pair[1].0 - 1))
        .collect()
}

/// The exon that holds a position.
pub fn exon_of(transcript: &Transcript, at: u64) -> Option<(u64, u64)> {
    exons_low_to_high(transcript)
        .into_iter()
        .find(|(start, end)| *start <= at && at <= *end)
}

/// The place of a sequence position on the spliced transcript, counted from 1 at the
/// transcript's 5' end.
pub fn spliced(transcript: &Transcript, at: u64) -> Option<u64> {
    let mut exons = exons_low_to_high(transcript);
    if transcript.strand == Strand::Reverse {
        exons.reverse();
    }
    let mut before = 0;
    for (start, end) in exons {
        if start <= at && at <= end {
            let into = match transcript.strand {
                Strand::Forward => at - start,
                Strand::Reverse => end - at,
            };
            return Some(before + into + 1);
        }
        before += end - start + 1;
    }
    None
}

/// The number, from 0, of the coding sequence's base at a sequence position.
pub fn coding_index(transcript: &Transcript, at: u64) -> Option<u64> {
    let coding = transcript.coding.as_ref()?;
    let place = spliced(transcript, at)?;
    (coding.cdna_start <= place && place <= coding.cdna_end).then(|| place - coding.cdna_start)
}

/// Whether the transcript carries a complete coding sequence the cache gives base for
/// base: no incomplete-coding-sequence flag, no edit, no start phase above zero, and as
/// many cached bases as the coding span is long.
pub fn complete(transcript: &Transcript) -> bool {
    transcript.coding.as_ref().is_some_and(|coding| {
        !transcript.cds_start_nf
            && !transcript.cds_end_nf
            && !transcript.edited
            && coding.start_phase <= 0
            && coding.sequence.len() as u64 == coding.cdna_end - coding.cdna_start + 1
    })
}

/// The base that pairs with `base`.
pub fn pair(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// A forward-strand base as the transcript reads it.
pub fn as_read(transcript: &Transcript, base: u8) -> u8 {
    match transcript.strand {
        Strand::Forward => base,
        Strand::Reverse => pair(base),
    }
}

/// The codon that holds a sequence position: its number from 0, the position's place in
/// it and its three bases as the transcript reads them.
pub fn codon_at(transcript: &Transcript, at: u64) -> Option<(u64, usize, [u8; 3])> {
    let index = coding_index(transcript, at)?;
    let sequence = &transcript.coding.as_ref()?.sequence;
    let first = (index - index % 3) as usize;
    let codon = sequence.get(first..first + 3)?;
    Some((
        index / 3,
        (index % 3) as usize,
        [codon[0], codon[1], codon[2]],
    ))
}

/// The codon at a position before and after a one-base substitution there, with the
/// codon's number; `None` unless the reference base is the cached one.
pub fn substituted(
    transcript: &Transcript,
    at: u64,
    reference: &[u8],
    alt: &[u8],
) -> Option<(u64, [u8; 3], [u8; 3])> {
    let ([reference], [alt]) = (reference, alt) else {
        return None;
    };
    let (number, offset, before) = codon_at(transcript, at)?;
    if before[offset] != as_read(transcript, *reference) {
        return None;
    }
    let mut after = before;
    after[offset] = as_read(transcript, *alt);
    Some((number, before, after))
}
