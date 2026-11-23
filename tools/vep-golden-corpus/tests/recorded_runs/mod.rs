// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The made-up runs kept as directories under `tests/fixtures/`: a corpus `c` as it stands
//! before the first call on it, the outputs the calls name, and the manifest call `n`
//! writes, `expected.<n>.json`. A test copies a directory to a scratch one, adds what the
//! directory does not keep (a cache shard, score files, a reference sequence) and makes
//! the calls in their order.

#![allow(dead_code)]

use crate::classify_support::counted;
use crate::common::{clean, fixture, formats_comparator, run, scratch};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use vep_golden_corpus::json::{loads, Value};
use vep_golden_corpus::TempDir;

pub const UNEXPLAINED: &str = "unexplained_residual";

/// A scratch copy of a fixture directory.
pub struct Recorded {
    scratch: TempDir,
    source: PathBuf,
}

/// One call of `classify --field` on a corpus `c`.
pub struct Call<'a> {
    /// The directory of the corpus and of the outputs inside the fixture, empty for the
    /// fixture's own.
    pub directory: &'a str,
    /// Ensembl VEP's output and vep-rs's, inside `directory`.
    pub outputs: (&'a str, &'a str),
    pub fields: &'a [&'a str],
    /// The score files of `--spliceai-table`, inside the scratch copy.
    pub tables: &'a [&'a str],
    /// The count of field divergences and their summary as standard output writes them.
    pub summary: (usize, &'a str),
    /// The manifest the call writes, inside `directory`.
    pub expected: &'a str,
}

/// A call on the fixture's own corpus and its `vep.txt` and `rs.txt`.
pub fn call<'a>(fields: &'a [&'a str], summary: (usize, &'a str), expected: &'a str) -> Call<'a> {
    Call {
        directory: "",
        outputs: ("vep.txt", "rs.txt"),
        fields,
        tables: &[],
        summary,
        expected,
    }
}

/// One entry of `field_divergences`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    pub location: String,
    pub feature: String,
    pub ordinals: Vec<i64>,
    pub field: String,
    pub vep: String,
    pub vep_rs: String,
    pub class: String,
}

fn copy_tree(source: &Path, target: &Path) {
    fs::create_dir_all(target).unwrap();
    for entry in fs::read_dir(source).unwrap() {
        let entry = entry.unwrap();
        let to = target.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &to);
        } else {
            fs::copy(entry.path(), to).unwrap();
        }
    }
}

impl Recorded {
    /// A scratch copy of the fixture directory `name`.
    pub fn new(name: &str) -> Recorded {
        let source = fixture(name);
        let scratch = scratch(&name.replace('/', "-"));
        copy_tree(&source, &scratch.path);
        Recorded { scratch, source }
    }

    /// The scratch copy's directory.
    pub fn root(&self) -> &Path {
        &self.scratch.path
    }

    /// The path of `relative` inside the scratch copy.
    pub fn path(&self, relative: &str) -> String {
        self.scratch
            .path
            .join(relative)
            .to_str()
            .unwrap()
            .to_string()
    }

    /// Writes `relative` inside the scratch copy, with the directories above it.
    pub fn write(&self, relative: &str, bytes: &[u8]) {
        let path = self.scratch.path.join(relative);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, bytes).unwrap();
    }

    /// Writes the fixture file `shard` as the one shard of `chromosome` in the cache of
    /// the corpus of `directory`.
    pub fn cache_shard(&self, directory: &str, chromosome: &str, shard: &str) {
        let bytes = fs::read(fixture("cache_shards").join(shard)).unwrap();
        let cache = Path::new(directory).join("c/json_cache/transcripts");
        let path = cache.join(chromosome).join("1-1000000.json");
        self.write(path.to_str().unwrap(), &bytes);
    }

    /// The arguments of a call.
    pub fn arguments(&self, call: &Call<'_>) -> Vec<String> {
        let inside = |name: &str| self.path(Path::new(call.directory).join(name).to_str().unwrap());
        let mut arguments = vec![
            "classify".to_string(),
            "--corpus".to_string(),
            inside("c"),
            "--vep-default".to_string(),
            inside(call.outputs.0),
            "--vep-rs-output".to_string(),
            inside(call.outputs.1),
        ];
        for field in call.fields {
            arguments.extend(["--field".to_string(), field.to_string()]);
        }
        for table in call.tables {
            arguments.extend(["--spliceai-table".to_string(), self.path(table)]);
        }
        // A call that names score files asks the format comparator's build.
        if !call.tables.is_empty() {
            arguments.extend(["--formats-comparator".to_string(), formats_comparator()]);
        }
        arguments
    }

    /// The arguments of a call beside the reference sequence `reference`, a file inside the
    /// scratch copy: the call's own, the reference and the format comparator's build, which
    /// such a call asks.
    pub fn arguments_beside(&self, call: &Call<'_>, reference: &str) -> Vec<String> {
        let mut arguments = self.arguments(call);
        arguments.extend(["--reference-fasta".to_string(), self.path(reference)]);
        if call.tables.is_empty() {
            arguments.extend(["--formats-comparator".to_string(), formats_comparator()]);
        }
        arguments
    }

    /// Runs `classify` as a call says, with no check of how it ends.
    pub fn classify(&self, call: &Call<'_>) -> crate::common::Run {
        run(&self.arguments(call))
    }

    /// Runs `classify` as a call says beside the reference sequence `reference`, with no
    /// check of how it ends.
    pub fn classify_beside(&self, call: &Call<'_>, reference: &str) -> crate::common::Run {
        run(&self.arguments_beside(call, reference))
    }

    /// The manifest of the corpus of `directory`, as text.
    pub fn manifest(&self, directory: &str) -> String {
        let manifest = Path::new(directory).join("c/manifest.json");
        fs::read_to_string(self.scratch.path.join(manifest)).unwrap()
    }

    /// Makes a call and holds it to its bytes: the run ends with exit code 0 on the line
    /// of no key divergence and the call's field divergences, and the manifest is the
    /// fixture's expected one. The entries of the expected manifest's `field_divergences`.
    pub fn holds(&self, call: &Call<'_>) -> Vec<Entry> {
        self.held(call, self.classify(call))
    }

    /// Makes a call beside the reference sequence `reference` and holds it to its bytes,
    /// as `holds` does a call without one.
    pub fn holds_beside(&self, call: &Call<'_>, reference: &str) -> Vec<Entry> {
        self.held(call, self.classify_beside(call, reference))
    }

    /// Holds the run `ran` of a call to the call's line and to its expected manifest.
    fn held(&self, call: &Call<'_>, ran: crate::common::Run) -> Vec<Entry> {
        let expected = self.source.join(call.directory).join(call.expected);
        let expected = fs::read_to_string(expected).unwrap();
        let label = format!("{} {}", self.source.display(), call.expected);
        let keys = counted((0, "{}"), (0, "{}"), "[]");
        let line = format!(
            "{}; field divergences {} ({})\n",
            keys.strip_suffix('\n').unwrap(),
            call.summary.0,
            call.summary.1
        );
        assert_eq!(ran, clean(0, &line), "{label}");
        assert_eq!(self.manifest(call.directory), expected, "{label}");
        let entries = entries(&expected);
        assert_eq!(entries.len(), call.summary.0, "{label}");
        entries
    }
}

/// The entries of the `field_divergences` of a manifest.
pub fn entries(manifest: &str) -> Vec<Entry> {
    let Ok(manifest) = loads(manifest) else {
        panic!("no JSON manifest");
    };
    let text =
        |entry: &Value, key: &str| entry.get(key).and_then(Value::as_str).unwrap().to_string();
    let listed = manifest.get("field_divergences").and_then(Value::as_list);
    listed
        .expect("a list of field divergences")
        .iter()
        .map(|entry| Entry {
            location: text(entry, "location"),
            feature: text(entry, "feature"),
            ordinals: entry
                .get("record_indices")
                .and_then(Value::as_list)
                .unwrap()
                .iter()
                .map(|ordinal| ordinal.as_i64().unwrap())
                .collect(),
            field: text(entry, "field"),
            vep: text(entry, "vep"),
            vep_rs: text(entry, "vep_rs"),
            class: text(entry, "expected_divergence"),
        })
        .collect()
}

/// The class of each entry under what `name` reads of it and its field, no two entries
/// under one name and field.
pub fn classes<'a>(
    entries: &'a [Entry],
    name: impl Fn(&'a Entry) -> &'a str,
) -> BTreeMap<(&'a str, &'a str), &'a str> {
    let classes: BTreeMap<(&str, &str), &str> = entries
        .iter()
        .map(|entry| ((name(entry), entry.field.as_str()), entry.class.as_str()))
        .collect();
    assert_eq!(classes.len(), entries.len(), "one entry a name and field");
    classes
}

/// The class of each entry under the name of its one record, of `names` in the records'
/// order, and its field.
pub fn classes_by_record<'a>(
    entries: &'a [Entry],
    names: &[&'a str],
) -> BTreeMap<(&'a str, &'a str), &'a str> {
    classes(entries, |entry| {
        assert_eq!(entry.ordinals.len(), 1, "one record an entry");
        names[usize::try_from(entry.ordinals[0]).unwrap()]
    })
}

/// The count of the entries of each class.
pub fn counts(entries: &[Entry]) -> BTreeMap<&str, usize> {
    let mut counts = BTreeMap::new();
    for entry in entries {
        *counts.entry(entry.class.as_str()).or_default() += 1;
    }
    counts
}
