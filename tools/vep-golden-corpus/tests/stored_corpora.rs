// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The committed corpora under `tests/golden/`, one case a corpus: the expected outputs of
//! every run hash to the digests the corpus's `provenance.json` records, the corpus is named
//! after its assembly and stands under its release, its cache states its assembly and its
//! sources, every coding transcript of the cache carries its peptide, and every transcript
//! of it with a UniProt identifier key carries `uniprot_isoform`. Two corpora are
//! held to what their runs are for: the shift-length set of `GRCh38-regulatory` and the
//! reference slice of `GRCh38-everything-fasta`. A table of made-up values a provenance
//! record lists is a file of its corpus, and a plugin corpus that holds a file named
//! `synthetic` ends in `-fixture`. `GRCh38-identifiers` is held to its selection: `select`
//! over its committed reference and input writes the manifest's values again. Six corpora
//! whose selection the tree cannot run again state no counts of the rows it chose. A cache's
//! `info.json` states the regulatory data, the variation columns and the sequences its slice
//! holds, but for one corpus whose record says what the file keeps of the whole cache's.

mod common;

use common::{
    committed_manifests, corpus_label, manifest, names, object, repository, run, scratch, text_of,
    under,
};
use flate2::read::{DeflateDecoder, MultiGzDecoder};
use sha2::{Digest, Sha256};
use std::fs;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use vep_golden_corpus::classify::cache::corpus_cache_dir;
use vep_golden_corpus::classify::intron_span::REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE;
use vep_golden_corpus::json::{loads, Object, Value};
use vep_golden_corpus::manifest::{is_corpus_path, is_name_list, read_manifest, target_run};
use vep_golden_corpus::pytext::py_str;
use vep_golden_corpus::rows::parse_extra;

/// The sha256 of `data` in hex digits.
fn sha256(data: &[u8]) -> String {
    Sha256::digest(data)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

fn bytes_of(path: &Path) -> Vec<u8> {
    fs::read(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The bytes a gzip file holds, every member in turn.
fn inflated(path: &Path) -> Vec<u8> {
    let mut plain = Vec::new();
    MultiGzDecoder::new(&bytes_of(path)[..])
        .read_to_end(&mut plain)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    plain
}

/// The bytes of the expected output `name` of a run: those of `name.gz` inflated where the
/// corpus stores it so, else the file's own.
fn expected_bytes(directory: &Path, name: &str) -> Vec<u8> {
    let packed = directory.join(format!("{name}.gz"));
    if packed.is_file() {
        inflated(&packed)
    } else {
        bytes_of(&directory.join(name))
    }
}

fn json_object(path: &Path) -> Object {
    object(&String::from_utf8(bytes_of(path)).expect("UTF-8"))
}

/// The value under `key`, which the object must hold.
fn held<'a>(object: &'a Object, key: &str, label: &str) -> &'a Value {
    object
        .get(key)
        .unwrap_or_else(|| panic!("{label}: no {key}"))
}

fn text<'a>(object: &'a Object, key: &str, label: &str) -> &'a str {
    held(object, key, label)
        .as_str()
        .unwrap_or_else(|| panic!("{label}: {key} is no text"))
}

/// The cache directory of a corpus as its manifest's layout places it.
fn cache_of(corpus: &Path, manifest: &Object) -> PathBuf {
    PathBuf::from(corpus_cache_dir(corpus.to_str().unwrap(), manifest))
}

/// The flag sets of a manifest, each with its name.
fn flag_sets(manifest: &Object) -> Vec<(&str, &Object)> {
    let sets = manifest.get("flag_sets").and_then(Value::as_list);
    sets.unwrap_or_default()
        .iter()
        .map(|set| {
            let set = set.as_object().expect("a flag set");
            (
                set.get("name").and_then(Value::as_str).expect("a name"),
                set,
            )
        })
        .collect()
}

/// Holds one corpus to its provenance record. The expected outputs are stored gzipped: their
/// bytes as Ensembl VEP wrote them hash to what the record holds, for the manifest's own run
/// and for every flag set it names, each under its own record in the `flag_sets` of
/// `provenance.json`. A tool corpus (a manifest with `tool_cases`) records a companion tool's
/// outputs instead: its `sha256` map is keyed by each file's path below `expected/`, every
/// case's record repeats its digest, and it carries a cache only when its manifest names one.
fn held_to_its_provenance(manifest_path: &Path) {
    let label = corpus_label(manifest_path);
    let corpus = manifest_path.parent().unwrap();
    let name = corpus.file_name().unwrap().to_str().unwrap();
    let release = corpus.parent().unwrap().file_name().unwrap();
    let release = release.to_str().unwrap();
    let provenance = json_object(&corpus.join("provenance.json"));
    let manifest = json_object(manifest_path);

    let digests_of = |record: &'_ Object, owner: &str| -> Vec<(String, String)> {
        let digests = held(record, "sha256", owner).as_object();
        let digests = digests.unwrap_or_else(|| panic!("{owner}: sha256 is no map"));
        let digest = |value: &Value| value.as_str().expect("a digest").to_string();
        let pairs = digests.iter();
        pairs
            .map(|(file, value)| (file.to_string(), digest(value)))
            .collect()
    };
    let mut runs = vec![("expected".to_string(), digests_of(&provenance, &label))];
    for (set_name, set) in flag_sets(&manifest) {
        let record = provenance
            .get("flag_sets")
            .and_then(|sets| sets.get(set_name))
            .and_then(Value::as_object)
            .unwrap_or_else(|| panic!("{name}: flag set {set_name} has no provenance record"));
        let directory = set.get("expected_dir").and_then(Value::as_str);
        let directory = directory.map_or_else(|| format!("expected_{set_name}"), str::to_string);
        runs.push((
            directory,
            digests_of(record, &format!("{label} [{set_name}]")),
        ));
    }
    for (directory, digests) in &runs {
        for (file, digest) in digests {
            assert_eq!(
                &sha256(&expected_bytes(&corpus.join(directory), file)),
                digest,
                "{name}/{directory}/{file} differs from its recorded digest"
            );
        }
    }

    // The directory is the assembly, with or without a suffix after a hyphen; an assembly
    // may itself hold hyphens, so the name is held to the manifest's.
    let assembly = text(&manifest, "assembly", &label);
    assert!(
        name == assembly || name.starts_with(&format!("{assembly}-")),
        "{name}: the directory does not start with the manifest's assembly {assembly}"
    );
    // The set directory names the Ensembl VEP release the reference was written with; the
    // manifest's release is the cache version both engines were given, which the record
    // carries as `cache_version`. A tool corpus runs no annotation and names the set's
    // release.
    let tool_cases = manifest.get("tool_cases");
    if tool_cases.is_some() {
        assert_eq!(text(&manifest, "release", &label), release, "{label}");
    } else {
        let reported = text(&provenance, "vep_version_reported", &label);
        let version = reported.rsplit(':').next().unwrap().trim();
        assert_eq!(
            version.split('.').next().unwrap(),
            release,
            "{name}: the reference was written with another Ensembl VEP release than the set's"
        );
        assert_eq!(
            text(&manifest, "release", &label),
            py_str(held(&provenance, "cache_version", &label)),
            "{name}: the manifest's release is not the record's cache_version"
        );
    }
    if let Some(cases) = tool_cases {
        let digests = held(&provenance, "sha256", &label);
        for case in cases.as_list().expect("a list of tool cases") {
            let case = case.as_object().expect("a tool case");
            let case_name = text(case, "name", &label);
            let expected = text(case, "expected", &label);
            let record = provenance
                .get("tool_cases")
                .and_then(|records| records.get(case_name))
                .and_then(Value::as_object)
                .unwrap_or_else(|| panic!("{name}: tool case {case_name} has no record"));
            let case_label = format!("{name}: tool case {case_name}");
            assert_eq!(
                text(record, "expected", &case_label),
                expected,
                "{case_label}"
            );
            let digest = sha256(&bytes_of(&corpus.join(expected)));
            let below = expected.strip_prefix("expected/").unwrap_or(expected);
            assert_eq!(
                digest,
                text(record, "sha256", &case_label),
                "{name}/{expected} differs from its recorded digest"
            );
            assert_eq!(
                Some(digest.as_str()),
                digests.get(below).and_then(Value::as_str),
                "{name}/{expected} differs from its recorded digest"
            );
        }
        if !manifest.contains_key("cache") {
            assert!(
                !corpus.join("json_cache").exists(),
                "{name}: a tool corpus without a cache in its manifest holds none"
            );
            return;
        }
    }
    let info = json_object(&cache_of(corpus, &manifest).join("info.json"));
    assert_eq!(text(&info, "assembly", &label), assembly, "{label}");
    assert!(
        info.get("source_versions").is_some_and(Value::truthy),
        "{label}: info.json must carry the cache source versions"
    );
    if let Some(fasta) = manifest.get("fasta") {
        let fasta = fasta.as_str().expect("a file name");
        let index = format!("{}.fai", fasta.strip_suffix(".gz").unwrap_or(fasta));
        assert!(corpus.join(fasta).is_file(), "{label}: no {fasta}");
        assert!(corpus.join(&index).is_file(), "{label}: no {index}");
    }
}

#[test]
fn committed_expected_files_match_their_provenance_digests() {
    let manifests = committed_manifests();
    assert!(manifests.len() >= 111, "{} corpora", manifests.len());
    for manifest in &manifests {
        held_to_its_provenance(manifest);
    }
}

/// The transcript records of a shard, gzipped or plain.
fn records_of(shard: &Path) -> Vec<Value> {
    let gzipped = shard.extension().is_some_and(|suffix| suffix == "gz");
    let bytes = if gzipped {
        inflated(shard)
    } else {
        bytes_of(shard)
    };
    let records = loads(&String::from_utf8(bytes).expect("UTF-8"));
    let records = records.unwrap_or_else(|error| panic!("{}: {error}", shard.display()));
    records
        .as_list()
        .expect("a list of transcript records")
        .to_vec()
}

/// The count of the transcripts of a shard that hold a coding sequence and no peptide.
fn without_a_peptide(shard: &Path) -> usize {
    let held = |effects: &Value, key: &str| effects.get(key).is_some_and(Value::truthy);
    records_of(shard)
        .iter()
        .filter_map(|record| record.get("variation_effect_feature_cache"))
        .filter(|effects| held(effects, "translateable_seq") && !held(effects, "peptide"))
        .count()
}

/// The count of the transcripts of a shard that carry a `swissprot`, `trembl` or `uniparc`
/// key and no `uniprot_isoform` key.
fn without_the_isoform_key(shard: &Path) -> usize {
    records_of(shard)
        .iter()
        .filter(|record| {
            ["swissprot", "trembl", "uniparc"]
                .iter()
                .any(|key| record.get(key).is_some())
        })
        .filter(|record| record.get("uniprot_isoform").is_none())
        .count()
}

/// Every committed cache as its corpus's label and the sum of `count` over its transcript
/// shards.
fn counted_over_committed_caches(count: fn(&Path) -> usize) -> Vec<(String, usize)> {
    let mut caches = Vec::new();
    for manifest_path in committed_manifests() {
        let corpus = manifest_path.parent().unwrap();
        let transcripts = cache_of(corpus, &json_object(&manifest_path)).join("transcripts");
        if !transcripts.is_dir() {
            continue;
        }
        let mut counted = 0;
        for chromosome in fs::read_dir(&transcripts).unwrap() {
            let chromosome = chromosome.unwrap().path();
            if !chromosome.is_dir() {
                continue;
            }
            for shard in fs::read_dir(&chromosome).unwrap() {
                counted += count(&shard.unwrap().path());
            }
        }
        caches.push((corpus_label(&manifest_path), counted));
    }
    caches
}

#[test]
fn committed_caches_carry_a_peptide_for_every_coding_transcript() {
    // Every committed cache carries the peptide beside each coding sequence, so the golden
    // tests read the sequence edits the measurement caches carry.
    let caches = counted_over_committed_caches(without_a_peptide);
    for (corpus, missing) in &caches {
        assert_eq!(
            *missing, 0,
            "{corpus}: coding transcripts without a peptide"
        );
    }
    assert!(caches.len() >= 110, "{} caches", caches.len());
}

#[test]
fn committed_caches_carry_the_isoform_key_beside_the_other_uniprot_keys() {
    // Ensembl VEP gives a translation its four UniProt identifier lists in one block
    // (`AnnotationSource/Database/Transcript.pm` lines 659 to 685), each `-` where the
    // translation has none. The cache converter copies each of the four a cache record
    // holds and the cache builder writes the four together, so in a cache that Ensembl VEP
    // release 102 or later dumped, a transcript with one of the first three and no
    // `uniprot_isoform` was converted without that key. The block writes
    // `_uniprot_isoform` from release 102 (lines 668 to 672 at release 102.0) and three
    // lists at release 101.0, so a sound conversion of a cache an older release dumped
    // carries the first three keys alone. The walk takes every committed cache: a corpus
    // cut from such a cache fails this test unless the walk leaves it out by its label.
    // An Ensembl Genomes cache is numbered by its own release (`Pipeline/DumpVEP/Dumper.pm`
    // line 66), not by the Ensembl VEP release that dumped it.
    let caches = counted_over_committed_caches(without_the_isoform_key);
    let without: Vec<String> = caches
        .iter()
        .filter(|(_, missing)| *missing > 0)
        .map(|(corpus, missing)| format!("{corpus} ({missing})"))
        .collect();
    assert!(
        without.is_empty(),
        "transcripts with a UniProt key and no uniprot_isoform: {}",
        without.join(", ")
    );
    assert!(caches.len() >= 110, "{} caches", caches.len());
}

/// Whether a directory of `tree` holds a file `named` accepts: the chromosome directories
/// of a cache's `regulatory` or `variations` tree.
fn holds_a_file(tree: &Path, named: fn(&str) -> bool) -> bool {
    let Ok(chromosomes) = fs::read_dir(tree) else {
        return false;
    };
    chromosomes
        .filter_map(|chromosome| fs::read_dir(chromosome.unwrap().path()).ok())
        .flatten()
        .any(|file| named(file.unwrap().file_name().to_str().unwrap()))
}

/// What the `info.json` of the cache at `cache` states that its directory does not bear
/// out, each as a phrase that opens with its key: regulatory data (`regulatory` true, a
/// `regulatory_schema` or a `cell_types` list) over no regulatory shard, `variation_cols`
/// or `var_type` over no variation file, and a `contig_order` that is absent or is not the
/// cache's transcript directories.
fn stated_beyond_the_slice(cache: &Path) -> Vec<String> {
    let info = json_object(&cache.join("info.json"));
    let mut stated = Vec::new();
    let shard: fn(&str) -> bool = |name| name.ends_with(".json") || name.ends_with(".json.gz");
    if !holds_a_file(&cache.join("regulatory"), shard) {
        if info.get("regulatory").is_some_and(Value::truthy) {
            stated.push("regulatory".to_string());
        }
        for key in ["regulatory_schema", "cell_types"] {
            if info.contains_key(key) {
                stated.push(key.to_string());
            }
        }
    }
    let variation_file: fn(&str) -> bool =
        |name| name == "all_vars.gz" || name.ends_with("_var.gz");
    if !holds_a_file(&cache.join("variations"), variation_file) {
        for key in ["variation_cols", "var_type"] {
            if info.contains_key(key) {
                stated.push(key.to_string());
            }
        }
    }
    let mut directories: Vec<String> = fs::read_dir(cache.join("transcripts"))
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_dir())
        .map(|path| path.file_name().unwrap().to_str().unwrap().to_string())
        .collect();
    directories.sort();
    match info.get("contig_order").and_then(Value::as_list) {
        None => stated.push("contig_order absent".to_string()),
        Some(order) => {
            let mut listed: Vec<String> = order
                .iter()
                .map(|name| name.as_str().expect("a sequence name").to_string())
                .collect();
            listed.sort();
            if listed != directories {
                stated.push("contig_order other than its transcript directories".to_string());
            }
        }
    }
    stated
}

/// The corpora whose cache keeps the metadata of the whole cache it is cut from, by label.
/// The record of each says under `cache_metadata` what that file keeps.
const CACHES_THAT_KEEP_THE_WHOLE_CACHES_METADATA: [&str; 1] = ["116/GRCh38-run-control"];

#[test]
fn committed_caches_state_what_their_slices_hold() {
    // `vep` takes `regulatory` for the cache's word that it holds regulatory data (a
    // `--regulatory_gff` run warns on that word that the file replaces the data) and builds
    // its frequency filter from `variation_cols` whether or not a variation file stands.
    // The cache converter states the regulatory keys with the tree it writes, and the
    // pruner states them, the variation keys and `contig_order` by what the slice holds
    // (`write_info` of `scripts/golden/prune_json_cache.py`).
    let mut caches = 0;
    let mut beyond = Vec::new();
    let mut kept = Vec::new();
    for manifest_path in committed_manifests() {
        let corpus = manifest_path.parent().unwrap();
        let cache = cache_of(corpus, &json_object(&manifest_path));
        if !cache.join("transcripts").is_dir() {
            continue;
        }
        caches += 1;
        let label = corpus_label(&manifest_path);
        let stated = stated_beyond_the_slice(&cache);
        if CACHES_THAT_KEEP_THE_WHOLE_CACHES_METADATA.contains(&label.as_str()) {
            kept.push((label, corpus.join("provenance.json"), stated));
        } else if !stated.is_empty() {
            beyond.push(format!("{label} ({})", stated.join(", ")));
        }
    }
    assert!(
        beyond.is_empty(),
        "caches whose info.json states more than the slice holds: {}",
        beyond.join("; ")
    );
    assert!(caches >= 110, "{caches} caches");
    assert_eq!(kept.len(), CACHES_THAT_KEEP_THE_WHOLE_CACHES_METADATA.len());
    for (label, record, stated) in &kept {
        assert!(
            !stated.is_empty(),
            "{label}: the cache states what its slice holds, so it leaves the list"
        );
        let record = json_object(record);
        for phrase in stated {
            let key = phrase.split(' ').next().unwrap();
            assert!(
                text(&record, "cache_metadata", label).contains(key),
                "{label}: cache_metadata does not name {key}"
            );
        }
    }
}

/// A committed corpus of release 116 by its directory name.
fn corpus_of_release_116(name: &str) -> PathBuf {
    repository().join("tests/golden/116").join(name)
}

#[test]
fn identifiers_manifest_holds_what_select_writes_over_its_committed_files() {
    // `GRCh38-identifiers` was selected from its own input, so `select` over the committed
    // reference and input, with the options the manifest records, writes the committed
    // input again and every value the committed manifest holds under a key of the
    // selection. The rows of the corpus's two-base substitutions count under
    // `mnv_or_complex`, apart from the rows of its insertions. The committed manifest holds
    // no `species`, which `select` writes for a manifest without one.
    let corpus = corpus_of_release_116("GRCh38-identifiers");
    let committed = json_object(&corpus.join("manifest.json"));
    let scratch = scratch("stored-identifiers-selection");
    let reference = scratch.path.join("default.txt");
    let rows = expected_bytes(&corpus.join("expected"), "default.txt");
    fs::write(&reference, rows).unwrap();
    let input = corpus.join("variants.vcf");
    let dataset = format!("identifiers={}={}", reference.display(), input.display());
    let out = under(&scratch.path, "out");
    fs::create_dir(&out).unwrap();
    let selected = run(&[
        "select",
        "--release",
        "116",
        "--assembly",
        "GRCh38",
        "--dataset",
        &dataset,
        "--out",
        &out,
        "--k",
        "100",
        "--seed",
        "1",
        "--covers",
        "hgvs",
        "--covers",
        "hgvs-offset-field",
        "--flag=--hgvs",
        "--fasta-name",
        "reference.fa.gz",
    ]);
    assert_eq!(selected.code, 0, "{}", selected.stderr);
    assert!(
        bytes_of(&Path::new(&out).join("variants.vcf")) == bytes_of(&input),
        "select writes another variants.vcf than the committed one"
    );
    let written = manifest(&out);
    let mut compared = Vec::new();
    for (key, value) in written.iter() {
        if key == "species" {
            continue;
        }
        assert_eq!(
            committed.get(key).map(text_of),
            Some(text_of(value)),
            "the committed {key} is not what select writes"
        );
        compared.push(key);
    }
    for key in ["datasets", "combinations", "variant_classes", "records"] {
        assert!(compared.contains(&key), "select wrote no {key}");
    }
    assert_eq!(
        text_of(held(&committed, "variant_classes", "GRCh38-identifiers")),
        r#"{"deletion": 133, "insertion": 310, "mnv_or_complex": 23, "snv": 92}"#
    );
}

/// Corpora of release 116 whose manifests state neither `variant_classes` nor
/// `family_strata`. Both maps count the rows a selection chose, which a manifest does not
/// list, so they are made again only by `select` over the selection's own reference and
/// input. `select` over the committed files of `GRCh38-input-hgvs` and
/// `GRCh38-plugin-utr-annotator` chooses other exemplars than their manifests record, and
/// the other four were selected from references the tree does not hold.
const WITHOUT_CHOSEN_ROW_COUNTS: [&str; 6] = [
    "GRCh37-merged",
    "GRCh37-refseq",
    "GRCh38-custom",
    "GRCh38-input-hgvs",
    "GRCh38-plugin-utr-annotator",
    "GRCh38-protein-features",
];

#[test]
fn a_manifest_select_cannot_write_again_states_no_counts_of_its_chosen_rows() {
    // A count of chosen rows that the tree cannot make again is not kept: each listed
    // manifest holds its records and combinations and neither map. A corpus that `select`
    // writes again states the maps `select` writes for it and leaves the list.
    for name in WITHOUT_CHOSEN_ROW_COUNTS {
        let manifest = json_object(&corpus_of_release_116(name).join("manifest.json"));
        for key in ["records", "combinations"] {
            held(&manifest, key, name);
        }
        for key in ["variant_classes", "family_strata"] {
            assert!(
                manifest.get(key).is_none(),
                "{name} states {key}, counts of chosen rows the tree cannot make again"
            );
        }
    }
}

#[test]
fn regulatory_shift_length_set_carries_the_field_on_transcript_rows_alone() {
    // The `regulatory-shift-length` set of the committed GRCh38-regulatory corpus runs
    // `--shift_3prime 1 --shift_length` beside `--regulatory`. Ensembl VEP sets the key in
    // `TranscriptVariationAllele_to_output_hash` alone (OutputFactory.pm:1704 at release
    // 116.3), so the stored default-format output carries `SHIFT_LENGTH` on a Transcript
    // row and on no RegulatoryFeature or MotifFeature row, and holds rows of both of those
    // kinds.
    let corpus = corpus_of_release_116("GRCh38-regulatory");
    let manifest = read_manifest(corpus.to_str().unwrap()).unwrap();
    let flag_set = target_run(&manifest, Some("regulatory-shift-length")).unwrap();
    assert_eq!(
        names(flag_set, "flags"),
        ["--regulatory", "--shift_3prime", "1", "--shift_length"]
    );
    let directory = corpus.join(text(flag_set, "expected_dir", "the set"));
    let output = String::from_utf8(inflated(&directory.join("default.txt.gz"))).expect("UTF-8");
    let kinds = ["Transcript", "RegulatoryFeature", "MotifFeature"];
    let mut carrying: [Vec<bool>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for line in output.split_terminator('\n') {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if let Some(kind) = kinds.iter().position(|kind| *kind == fields[5]) {
            carrying[kind].push(parse_extra(fields[13]).get("SHIFT_LENGTH").is_some());
        }
    }
    assert!(carrying[0].contains(&true));
    for kind in 1..kinds.len() {
        assert!(
            !carrying[kind].is_empty() && !carrying[kind].contains(&true),
            "{}",
            kinds[kind]
        );
    }
}

#[test]
fn refseq_intron_span_corpus_documents_the_edited_model_row_in_both_runs() {
    // The committed GRCh38-refseq-intron-span corpus holds one record, a deletion with a
    // whole intron of an edited RefSeq model between its two coding ends, as a VCF line
    // and, in the flag set `region`, as a region line that holds no base. The manifest
    // names the class in `masks` and documents one key under it; the set carries no
    // entries of its own, so a reader takes the manifest's for it. In both stored
    // default-format outputs Ensembl VEP's row of that key is the one row of an edited
    // model: it carries `BAM_EDIT=OK`, a `USED_REF` longer than its `CDS_position` range,
    // the three terms of a whole intron and the two coding terms the class names.
    let class = REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE;
    let corpus = corpus_of_release_116("GRCh38-refseq-intron-span");
    let manifest = read_manifest(corpus.to_str().unwrap()).unwrap();
    assert_eq!(names(&manifest, "flags"), ["--refseq", "--hgvs"]);
    assert_eq!(names(&manifest, "masks"), [class]);
    let entries = held(&manifest, "divergences", "the manifest")
        .as_list()
        .expect("a list");
    assert_eq!(entries.len(), 1);
    let entry = entries[0].as_object().expect("an entry");
    assert_eq!(text(entry, "expected_divergence", "the entry"), class);
    let feature = text(entry, "feature", "the entry");
    let location = text(entry, "location", "the entry");
    let region = target_run(&manifest, Some("region")).unwrap();
    assert_eq!(
        names(region, "flags"),
        ["--refseq", "--hgvs", "--format", "region"]
    );
    assert!(!region.contains_key("divergences") && !region.contains_key("masks"));
    let input = corpus.join(text(region, "input", "the set"));
    let line = String::from_utf8(bytes_of(&input)).expect("UTF-8");
    assert_eq!(line, format!("{location}:1/-\n"));
    for directory in ["expected", text(region, "expected_dir", "the set")] {
        let stored = inflated(&corpus.join(directory).join("default.txt.gz"));
        let output = String::from_utf8(stored).expect("UTF-8");
        let rows: Vec<Vec<&str>> = output
            .split_terminator('\n')
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect())
            .collect();
        let edited: Vec<&Vec<&str>> = rows
            .iter()
            .filter(|row| parse_extra(row[13]).get("BAM_EDIT").is_some())
            .collect();
        assert_eq!(edited.len(), 1, "{directory}");
        let row = edited[0];
        assert_eq!((row[1], row[2], row[4]), (location, "-", feature));
        let extra = parse_extra(row[13]);
        assert_eq!(extra.get("BAM_EDIT"), Some("OK"), "{directory}");
        let (start, end) = row[8].split_once('-').expect("a CDS range");
        let coding_bases = end.parse::<usize>().unwrap() - start.parse::<usize>().unwrap() + 1;
        let used_ref = extra.get("USED_REF").expect("USED_REF");
        assert!(used_ref.len() > coding_bases, "{directory}");
        let terms: Vec<&str> = row[6].split(',').collect();
        for term in [
            "splice_acceptor_variant",
            "splice_donor_variant",
            "intron_variant",
            "stop_lost",
            "protein_altering_variant",
        ] {
            assert!(terms.contains(&term), "{directory}: {term}");
        }
        assert!(!terms.contains(&"inframe_deletion"), "{directory}");
    }
}

/// The bases of 1-based inclusive `[start, end]` of a bgzipped FASTA, read through its two
/// indexes: the `.fai` gives the offset of a base in the decompressed text and the `.gzi`
/// the compressed block holding that offset, so the read opens the file at that block and
/// decompresses none before it.
fn indexed_fasta_bases(fasta: &Path, contig: &str, start: u64, end: u64) -> String {
    let beside = |suffix: &str| PathBuf::from(format!("{}{suffix}", fasta.display()));
    let listing = String::from_utf8(bytes_of(&beside(".fai"))).expect("UTF-8");
    let entry: Vec<&str> = listing
        .lines()
        .map(|line| line.split('\t').collect::<Vec<&str>>())
        .find(|fields| fields[0] == contig)
        .unwrap_or_else(|| panic!("{contig} is not a sequence of {}", fasta.display()));
    let number = |place: usize| entry[place].parse::<u64>().expect("a count");
    let (offset, line_bases, line_bytes) = (number(2), number(3), number(4));
    let place = |base: u64| offset + (base - 1) / line_bases * line_bytes + (base - 1) % line_bases;
    let (first, last) = (place(start), place(end));
    let index = bytes_of(&beside(".gzi"));
    let word = |at: usize| u64::from_le_bytes(index[at..at + 8].try_into().unwrap());
    // The index lists every block after the first as (compressed offset, decompressed offset).
    let listed = usize::try_from(word(0)).unwrap();
    let blocks = (0..listed).map(|block| (word(8 + 16 * block), word(16 + 16 * block)));
    let (compressed, decompressed) = std::iter::once((0, 0))
        .chain(blocks)
        .filter(|(_, decompressed)| *decompressed <= first)
        .max_by_key(|(_, decompressed)| *decompressed)
        .expect("the first block");
    let mut text = Vec::new();
    let mut handle = fs::File::open(fasta).unwrap();
    handle.seek(SeekFrom::Start(compressed)).unwrap();
    while decompressed + text.len() as u64 <= last {
        let mut header = [0; 18];
        handle.read_exact(&mut header).unwrap();
        assert_eq!(
            header[..4],
            [0x1f, 0x8b, 0x08, 0x04],
            "a bgzf member with the extra field"
        );
        let block_size = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        let mut member = vec![0; block_size - 18];
        handle.read_exact(&mut member).unwrap();
        let deflated = &member[..member.len() - 8];
        DeflateDecoder::new(deflated)
            .read_to_end(&mut text)
            .unwrap();
    }
    let slice = &text[usize::try_from(first - decompressed).unwrap()
        ..=usize::try_from(last - decompressed).unwrap()];
    let bases: Vec<u8> = slice
        .iter()
        .copied()
        .filter(|byte| *byte != b'\n')
        .collect();
    assert!(bases.is_ascii());
    String::from_utf8(bases).unwrap()
}

#[test]
fn everything_fasta_slice_holds_a_base_at_every_record() {
    // The GRCh38-everything-fasta corpus runs the records of GRCh38-everything, its
    // `variants.vcf` byte for byte, with a reference slice of its own, so that `HGVSc` and
    // `HGVSp` are among its compared fields. The slice, read through its indexes at each
    // record's position, holds that record's reference allele, one of the four bases at
    // every position of it: no record's notation is written from an unknown base. The first
    // line of the contig lies outside every region the slice keeps, and the same read
    // returns the unknown letter there, so a record in a masked region would fail the
    // comparison with its allele.
    let corpus = corpus_of_release_116("GRCh38-everything-fasta");
    let manifest = json_object(&corpus.join("manifest.json"));
    assert_eq!(names(&manifest, "flags"), ["--everything"]);
    let compared = names(&manifest, "fields_compared");
    assert!(["HGVSc", "HGVSp"]
        .iter()
        .all(|field| compared.iter().any(|held| held == field)));
    let vcf = bytes_of(&corpus.join("variants.vcf"));
    assert_eq!(
        vcf,
        bytes_of(&corpus_of_release_116("GRCh38-everything").join("variants.vcf"))
    );
    let fasta = corpus.join(text(&manifest, "fasta", "the manifest"));
    let vcf = String::from_utf8(vcf).expect("UTF-8");
    let records: Vec<Vec<&str>> = vcf
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(records.len(), 204);
    let mut held = 0;
    for record in &records {
        let (chrom, pos, reference) = (record[0], record[1], record[3]);
        let pos: u64 = pos.parse().expect("a position");
        let bases = indexed_fasta_bases(&fasta, chrom, pos, pos + reference.len() as u64 - 1);
        assert_eq!(
            bases,
            reference.to_uppercase(),
            "{chrom}:{pos}: the slice does not hold the record's reference allele"
        );
        held += usize::from(bases.bytes().all(|base| b"ACGT".contains(&base)));
    }
    assert_eq!(held, records.len());
    assert_eq!(
        indexed_fasta_bases(&fasta, records[0][0], 1, 60),
        "N".repeat(60)
    );
}

/// Every file below `directory`, at any depth, in path order.
fn files_below(directory: &Path) -> Vec<PathBuf> {
    let mut found = Vec::new();
    let mut entries: Vec<PathBuf> = fs::read_dir(directory)
        .unwrap_or_else(|error| panic!("{}: {error}", directory.display()))
        .map(|entry| entry.unwrap().path())
        .collect();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            found.extend(files_below(&path));
        } else if path.is_file() {
            found.push(path);
        }
    }
    found
}

fn name_of(path: &Path) -> &str {
    path.file_name().unwrap().to_str().unwrap()
}

#[test]
fn a_provenance_synthetic_table_is_a_file_of_its_corpus() {
    // A corpus's provenance record lists under `synthetic_tables` the tables of made-up
    // values whose names do not say `synthetic`: distinct paths inside the corpus, each a
    // file of that corpus. `GRCh38-regulatory` lists the table its pLI run reads.
    let mut listed: Vec<(String, Vec<String>)> = Vec::new();
    for path in files_below(&repository().join("tests/golden")) {
        if name_of(&path) != "provenance.json" {
            continue;
        }
        let record = json_object(&path);
        let Some(tables) = record.get("synthetic_tables") else {
            continue;
        };
        let corpus = path.parent().unwrap();
        let name = name_of(corpus);
        let paths = names(&record, "synthetic_tables");
        let distinct: std::collections::BTreeSet<&String> = paths.iter().collect();
        assert!(
            is_name_list(tables) && !paths.is_empty() && distinct.len() == paths.len(),
            "{name}: synthetic_tables is not a list of distinct paths"
        );
        for (table, text) in tables.as_list().unwrap().iter().zip(&paths) {
            assert!(
                is_corpus_path(table),
                "{name}: synthetic table {text} is not a path inside the corpus"
            );
            assert!(
                corpus.join(text).is_file(),
                "{name}: synthetic table {text} is not a file of the corpus"
            );
        }
        listed.push((name.to_string(), paths));
    }
    let regulatory = listed.iter().find(|(name, _)| name == "GRCh38-regulatory");
    assert_eq!(
        regulatory.map(|(_, tables)| tables.as_slice()),
        Some(&["plugins/pLI_values.txt".to_string()][..])
    );
}

#[test]
fn a_plugin_corpus_with_a_synthetic_file_ends_in_fixture() {
    // A plugin corpus that holds a file named `synthetic`, a made-up file in the place of
    // data the tree does not redistribute, says so in its directory name, which ends in
    // `-fixture`.
    let mut holders = Vec::new();
    for release in fs::read_dir(repository().join("tests/golden")).unwrap() {
        let release = release.unwrap().path();
        if !release.is_dir() {
            continue;
        }
        for corpus in fs::read_dir(&release).unwrap() {
            let corpus = corpus.unwrap().path();
            let holds_one = || {
                files_below(&corpus)
                    .iter()
                    .any(|file| name_of(file).contains("synthetic"))
            };
            if corpus.is_dir() && name_of(&corpus).contains("-plugin-") && holds_one() {
                holders.push(name_of(&corpus).to_string());
            }
        }
    }
    holders.sort();
    assert!(
        !holders.is_empty(),
        "no plugin corpus holds a synthetic file"
    );
    let unmarked: Vec<&String> = holders
        .iter()
        .filter(|name| !name.ends_with("-fixture"))
        .collect();
    assert!(unmarked.is_empty(), "{unmarked:?}");
}
