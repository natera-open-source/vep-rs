// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Sequence names of a cache beyond human, through the `vep` binary: a record
//! on a chromosome named by a word, on a plasmid and on an organelle is placed
//! under the cache's own name, under a synonym of the cache's
//! `chr_synonyms.txt` or of a `--synonyms` file and under the name behind
//! `chr`, and keeps the input's spelling; a name the cache lacks is skipped
//! with Ensembl VEP's warning. Each test builds its own cache in the layout of
//! `fixtures/species/json_cache` (its README describes it): three made-up
//! sequences, each holding one copy of that fixture's coding transcript under
//! a stable id of its own. Each test cites the Ensembl VEP source
//! (release/116.3) for the behaviour it proves.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

#[cfg(feature = "parquet")]
use vep_core::coordinate::is_cache_path_component;

mod common;

/// The cache's sequences and the transcript each holds: a chromosome named by
/// a word that begins with `Chr`, a plasmid named with a dot and an organelle.
const SEQUENCES: [(&str, &str); 3] = [
    ("Chromosome", "SYNT0000101"),
    ("pSYN1.2", "SYNT0000102"),
    ("MtDNA", "SYNT0000103"),
];

/// The base every record changes, the middle base of the transcript's second
/// codon: a missense change.
const POSITION: u64 = 1014;

/// Writes the cache under `dir`: one shard a sequence holding the species
/// fixture's coding transcript renamed, an `info.json` listing `contig_order`,
/// `synonyms` as the cache's `chr_synonyms.txt` when given, and the checksum
/// manifest over all of it.
fn write_cache(dir: &Path, contig_order: &[&str], synonyms: Option<&str>) -> PathBuf {
    let cache = dir.join("json_cache");
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/species/json_cache/transcripts/1/1-1000000.json");
    let shard: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(fixture).unwrap()).unwrap();
    for (name, transcript) in SEQUENCES {
        let mut model = shard[0].clone();
        model["chr"] = name.into();
        model["stable_id"] = transcript.into();
        let chromosome = cache.join("transcripts").join(name);
        fs::create_dir_all(&chromosome).unwrap();
        fs::write(
            chromosome.join("1-1000000.json"),
            serde_json::json!([model]).to_string(),
        )
        .unwrap();
    }
    let info = serde_json::json!({
        "species": "synthetic_species",
        "assembly": "SYN1",
        "cache_version": 116,
        "source_type": "ensembl",
        "source_versions": {"assembly": "SYN1", "genebuild": "synthetic_1"},
        "contig_order": contig_order,
    });
    fs::write(cache.join("info.json"), info.to_string()).unwrap();
    if let Some(text) = synonyms {
        fs::write(cache.join("chr_synonyms.txt"), text).unwrap();
    }
    vep_core::cache_manifest::CacheManifest::write(&cache).unwrap();
    cache
}

/// Runs the binary quietly against `cache` over one record a name of
/// `written`, each `at_<name>` at [`POSITION`], writing `<dir>/<label>.<ext>`
/// and `<dir>/<label>_warnings.txt`, and asserts success.
fn run(dir: &Path, cache: &Path, label: &str, ext: &str, written: &[&str], extra: &[&str]) {
    let records: String = written
        .iter()
        .map(|name| format!("{name}\t{POSITION}\tat_{name}\tC\tT\t.\t.\t.\n"))
        .collect();
    let input = dir.join(format!("{label}.vcf"));
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"),
    )
    .unwrap();
    let output: Output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(dir.join(format!("{label}.{ext}")))
        .arg("--warning_file")
        .arg(dir.join(format!("{label}_warnings.txt")))
        .arg("--json_cache")
        .arg(cache)
        .args([
            "--species",
            "synthetic_species",
            "--assembly",
            "SYN1",
            "--offline",
            "--no_stats",
            "--force_overwrite",
            "--quiet",
        ])
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{label}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// [`run`] in the default format: each row's record name, Location, Feature
/// and Consequence, and the lines of the warnings file.
fn annotate(
    dir: &Path,
    cache: &Path,
    label: &str,
    written: &[&str],
    extra: &[&str],
) -> (Vec<[String; 4]>, Vec<String>) {
    run(dir, cache, label, "out", written, extra);
    let rows = fs::read_to_string(dir.join(format!("{label}.out")))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            [columns[0], columns[1], columns[4], columns[6]].map(str::to_string)
        })
        .collect();
    let warnings = fs::read_to_string(dir.join(format!("{label}_warnings.txt")))
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    (rows, warnings)
}

/// The row of a record written `written` and placed in `transcript`: the
/// Location carries the input's spelling of the name.
fn placed(written: &str, transcript: &str) -> [String; 4] {
    [
        format!("at_{written}"),
        format!("{written}:{POSITION}"),
        transcript.to_string(),
        "missense_variant".to_string(),
    ]
}

/// The row `--dont_skip` keeps of a record on a name the cache lacks.
fn intergenic(written: &str) -> [String; 4] {
    [
        format!("at_{written}"),
        format!("{written}:{POSITION}"),
        "-".to_string(),
        "intergenic_variant".to_string(),
    ]
}

/// The warning of `Parser.pm:542-546` for the record written `written` on
/// input line `line`.
fn not_found(line: usize, written: &str) -> String {
    format!(
        "WARNING: line {line} skipped ({written} {POSITION} at_{written} C T . . .): Chromosome \
         {written} not found in annotation sources or synonyms; chromosome {written} does not \
         overlap any features"
    )
}

/// A record written with the cache's own name of each sequence is found among
/// the cache's chromosome directories as written (`CacheDir.pm:431`,
/// `Parser.pm:777`, `BaseVEP.pm:706`): one row in that sequence's transcript,
/// its Location the name as written, and no warning. A name that begins with
/// `Chr` loses none of its letters on the way.
#[test]
fn a_record_is_annotated_on_a_sequence_named_as_the_cache_names_it() {
    let dir = tempfile::tempdir().unwrap();
    let cache = write_cache(dir.path(), &[], None);
    let names = SEQUENCES.map(|(name, _)| name);
    let (rows, warnings) = annotate(dir.path(), &cache, "own_names", &names, &[]);
    assert_eq!(
        rows,
        SEQUENCES.map(|(name, transcript)| placed(name, transcript))
    );
    assert_eq!(warnings, [""; 0]);
}

/// A record named by a synonym of the cache's `chr_synonyms.txt`
/// (`CacheDir.pm:254`) is placed on the sequence the line links it to
/// (`BaseVEP.pm:711-716`) and keeps its spelling in the Location; a bare `M`
/// is the one name rewritten, to `MT`, once `MT` is a synonym of a cache
/// sequence and neither `M` nor `chrM` is (`Parser.pm:794-796`). A
/// `--synonyms` file replaces the cache's table, since `CacheDir.pm:254` reads
/// the cache's file only when the flag is unset: its names are placed and the
/// cache's own synonyms are skipped as not found, while a sequence's own name
/// needs no table.
#[test]
fn a_record_named_by_a_synonym_is_annotated_and_keeps_its_spelling() {
    let dir = tempfile::tempdir().unwrap();
    let cache = write_cache(
        dir.path(),
        &[],
        Some("Chromosome\tNZ_SYN000000.1\npSYN1.2\tNZ_SYN000002.1\nMtDNA\tMT\n"),
    );
    let [chromosome, plasmid, organelle] = SEQUENCES.map(|(_, transcript)| transcript);
    let written = ["NZ_SYN000000.1", "NZ_SYN000002.1", "MT", "M", "pSYN1.2"];
    let (rows, warnings) = annotate(dir.path(), &cache, "cache_table", &written, &[]);
    let mut renamed = placed("MT", organelle);
    renamed[0] = "at_M".to_string();
    assert_eq!(
        rows,
        [
            placed("NZ_SYN000000.1", chromosome),
            placed("NZ_SYN000002.1", plasmid),
            placed("MT", organelle),
            renamed,
            placed("pSYN1.2", plasmid),
        ]
    );
    assert_eq!(warnings, [""; 0]);

    let table = dir.path().join("synonyms.txt");
    fs::write(&table, "Chromosome\tSYN_ALIAS\n").unwrap();
    let written = ["SYN_ALIAS", "NZ_SYN000002.1", "MT", "Chromosome"];
    let (rows, warnings) = annotate(
        dir.path(),
        &cache,
        "flag_table",
        &written,
        &["--synonyms", table.to_str().unwrap()],
    );
    assert_eq!(
        rows,
        [
            placed("SYN_ALIAS", chromosome),
            placed("Chromosome", chromosome)
        ]
    );
    assert_eq!(
        warnings,
        [not_found(2, "NZ_SYN000002.1"), not_found(3, "MT")]
    );
}

/// A record on a name the cache lacks, and no synonym or `chr` spelling
/// reaches, fails every test of `_have_chr` (`Parser.pm:777`, `:784`, `:791`
/// and `:794`, so `:801` returns 0) and is skipped with the warning of
/// `:542-546`: the three mitochondrial spellings on a cache whose organelle
/// has another name, a sequence's name without its own first three letters,
/// in another case or cut short. `chr` before a whole name is stripped
/// (`:789-791`, `BaseVEP.pm:722-726`) and the record keeps the prefix in its
/// Location. Under `--dont_skip` each skipped record is kept (`Parser.pm:189`)
/// and written as `intergenic_variant` under the name as written, with the
/// same warning.
#[test]
fn a_record_on_a_sequence_the_cache_lacks_is_reported_as_ensembl_vep_reports_it() {
    let dir = tempfile::tempdir().unwrap();
    let cache = write_cache(dir.path(), &[], None);
    let [chromosome, plasmid, organelle] = SEQUENCES.map(|(_, transcript)| transcript);
    let lacked = ["MT", "M", "chrM", "omosome", "chromosome", "pSYN1", "mtdna"];
    let mut written = vec!["Chromosome"];
    written.extend(lacked);
    written.extend(["chrChromosome", "CHRpSYN1.2", "chrMtDNA"]);
    let warned: Vec<String> = lacked
        .iter()
        .enumerate()
        .map(|(i, name)| not_found(i + 2, name))
        .collect();

    let (rows, warnings) = annotate(dir.path(), &cache, "skipped", &written, &[]);
    assert_eq!(
        rows,
        [
            placed("Chromosome", chromosome),
            placed("chrChromosome", chromosome),
            placed("CHRpSYN1.2", plasmid),
            placed("chrMtDNA", organelle),
        ]
    );
    assert_eq!(warnings, warned);

    let (rows, warnings) = annotate(dir.path(), &cache, "kept", &written, &["--dont_skip"]);
    let mut expected = vec![placed("Chromosome", chromosome)];
    expected.extend(lacked.map(intergenic));
    expected.extend([
        placed("chrChromosome", chromosome),
        placed("CHRpSYN1.2", plasmid),
        placed("chrMtDNA", organelle),
    ]);
    assert_eq!(rows, expected);
    assert_eq!(warnings, warned);
}

/// The Parquet footer's `contig_order` of a run on the cache.
#[cfg(feature = "parquet")]
fn footer_contig_order(out: &Path) -> String {
    use parquet::file::reader::FileReader;
    let parts = vep_io::output::parquet::part_files(out).unwrap();
    let reader = parquet::file::serialized_reader::SerializedFileReader::new(
        fs::File::open(&parts[0]).unwrap(),
    )
    .unwrap();
    reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .unwrap()
        .iter()
        .find(|kv| kv.key == "contig_order")
        .and_then(|kv| kv.value.clone())
        .unwrap()
}

/// A cache whose `info.json` lists the three names in `contig_order` loads,
/// and a Parquet run ranks its records in that order, not by name, and lists
/// the order in its footer. Each name is one ordinary path component, the
/// rule a chromosome directory is held to; a name with a separator or a `..`
/// is none, names no directory of the cache (`CacheDir.pm:431` lists the
/// directories whose name does not start with a dot) and is skipped as not
/// found.
#[cfg(feature = "parquet")]
#[test]
fn the_sequences_are_listed_in_the_caches_contig_order() {
    let dir = tempfile::tempdir().unwrap();
    let order = ["Chromosome", "pSYN1.2", "MtDNA"];
    let cache = write_cache(dir.path(), &order, None);
    let flags = ["--output_format", "parquet"];
    run(dir.path(), &cache, "ordered", "parquet", &order, &flags);
    let out = dir.path().join("ordered.parquet");
    assert_eq!(footer_contig_order(&out), "Chromosome,pSYN1.2,MtDNA");
    assert_eq!(
        vep_io::output::parquet::sort_keys(&out).unwrap(),
        [(0, 1014), (1, 1014), (2, 1014)]
    );

    for name in order {
        assert!(is_cache_path_component(name), "{name}");
    }
    let outside = ["..", "Chromosome/..", "../Chromosome", "pSYN1.2\\MtDNA"];
    for name in outside {
        assert!(!is_cache_path_component(name), "{name}");
    }
    let mut written = vec!["Chromosome"];
    written.extend(outside);
    let (rows, warnings) = annotate(dir.path(), &cache, "outside", &written, &[]);
    assert_eq!(rows, [placed("Chromosome", SEQUENCES[0].1)]);
    let warned: Vec<String> = outside
        .iter()
        .enumerate()
        .map(|(i, name)| not_found(i + 2, name))
        .collect();
    assert_eq!(warnings, warned);
}
