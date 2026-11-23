// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The classification: the field groups two cache builds differ in on each transcript
//! both dumps hold, its table, its summary and its line on standard output.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::Write;
use std::os::unix::ffi::OsStringExt;
use std::path::PathBuf;

use crate::difflib::changed_lines;
use crate::dumps::{key_of, load, peptide_of};
use crate::jsoncache::{read, scan, shards, Event, Peptide};
use crate::py::{os_text, path_join, push_json_str, repr_str_list};

/// The name in brackets that opens every message of the classification.
const TAG: &str = "[classify_storable_transcript_differences]";

const PEPTIDE_KEYS: [&str; 2] = ["peptide", "translateable_seq"];
const PREDICTION_KEYS: [&str; 9] = [
    "sift",
    "polyphen_humdiv",
    "polyphen_humvar",
    "protein_function_predictions",
    "matrix",
    "matrix_compressed",
    "peptide_length",
    "sub_analysis",
    "translation_md5",
];
const FEATURE_KEYS: [&str; 17] = [
    "hseqname",
    "_display_label",
    "analysis",
    "protein_features",
    "hstart",
    "hend",
    "score",
    "percent_id",
    "p_value",
    "logic_name",
    "db",
    "db_version",
    "interpro_ac",
    "idesc",
    "ilabel",
    "start",
    "end",
];

/// What the command line asks for. Paths are normalised as `pathlib` writes them.
#[derive(Debug, Clone, Default)]
pub struct Request {
    pub version_a: String,
    pub dump_a: Vec<u8>,
    pub version_b: String,
    pub dump_b: Vec<u8>,
    pub out_dir: Vec<u8>,
    pub json_cache: Option<Vec<u8>>,
    pub json_cache_label: Option<String>,
    pub threads: usize,
}

/// How a run ends when it fails: the exit code and a line for standard error.
pub struct Failure {
    pub code: i32,
    pub line: String,
}

impl Failure {
    /// A failure CPython reports with a traceback, in one line of the tool's own.
    fn refused(detail: impl std::fmt::Display) -> Failure {
        Failure {
            code: 1,
            line: format!("ERROR: [vep-storable-transcript-classes] {detail}"),
        }
    }
}

fn os_path(bytes: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(bytes.to_vec()))
}

fn shown(path: &[u8]) -> String {
    String::from_utf8_lossy(path).into_owned()
}

/// A count per name, the names in the order of their first counting, as a
/// `collections.Counter` keeps them.
#[derive(Default)]
struct Counter {
    entries: Vec<(String, usize)>,
}

impl Counter {
    fn add(&mut self, name: &str) {
        match self.entries.iter_mut().find(|(known, _)| known == name) {
            Some(entry) => entry.1 += 1,
            None => self.entries.push((name.to_string(), 1)),
        }
    }

    /// `json.dumps(dict(counter))`.
    fn dumps(&self) -> String {
        let mut out = String::from("{");
        for (index, (name, count)) in self.entries.iter().enumerate() {
            if index > 0 {
                out.push_str(", ");
            }
            push_json_str(&mut out, name);
            write!(out, ": {count}").expect("a string takes text");
        }
        out.push('}');
        out
    }

    /// The same object as `json.dumps(..., indent=2)` writes it one level deep.
    fn dumps_indented(&self) -> String {
        if self.entries.is_empty() {
            return "{}".to_string();
        }
        let mut out = String::from("{");
        for (index, (name, count)) in self.entries.iter().enumerate() {
            out.push_str(if index > 0 { ",\n    " } else { "\n    " });
            push_json_str(&mut out, name);
            write!(out, ": {count}").expect("a string takes text");
        }
        out.push_str("\n  }");
        out
    }
}

type Blocks = HashMap<String, Vec<String>>;

fn load_dump(path: &[u8]) -> Result<Blocks, Failure> {
    let bytes = std::fs::read(os_path(path))
        .map_err(|error| Failure::refused(format!("{}: {}", shown(path), os_text(&error))))?;
    // `errors="replace"`: a byte sequence that is not UTF-8 reads as U+FFFD.
    load(&String::from_utf8_lossy(&bytes))
        .map_err(|_| Failure::refused(format!("{}: a block's first line has no end", shown(path))))
}

/// The classes of one transcript and the keys of its class `other`.
struct Classes {
    classes: BTreeSet<&'static str>,
    other_keys: BTreeSet<String>,
}

fn classify(a: &[String], b: &[String]) -> Classes {
    let a_lines: Vec<&str> = a.iter().map(String::as_str).collect();
    let b_lines: Vec<&str> = b.iter().map(String::as_str).collect();
    let (removed, added) = changed_lines(&a_lines, &b_lines);
    // A removed line that opens with `--` and an added one that opens with `++` read
    // as the diff's own file lines and are passed over.
    let lines = removed
        .into_iter()
        .map(|at| a_lines[at])
        .filter(|line| !line.starts_with("--"))
        .chain(
            added
                .into_iter()
                .map(|at| b_lines[at])
                .filter(|line| !line.starts_with("++")),
        );
    let mut classes = BTreeSet::new();
    let mut other_keys = BTreeSet::new();
    let mut any = false;
    for line in lines {
        any = true;
        let Some(key) = key_of(line) else {
            continue;
        };
        if PEPTIDE_KEYS.contains(&key) {
            classes.insert("peptide");
        } else if PREDICTION_KEYS.contains(&key) {
            classes.insert("predictions");
        } else if FEATURE_KEYS.contains(&key) {
            classes.insert("protein_features");
        } else {
            classes.insert("other");
            other_keys.insert(key.to_string());
        }
    }
    if !any {
        classes.insert("identical_after_normalisation");
    }
    Classes {
        classes,
        other_keys,
    }
}

/// What one shard says about the wanted ids, or why it cannot be read.
type Scanned = Result<Vec<Event>, String>;

/// The stored translation of every wanted id the cache holds: the first transcript of
/// the id in the order of the shards.
fn json_cache_peptides(
    cache: &[u8],
    wanted: &HashSet<String>,
    threads: usize,
) -> Result<HashMap<String, Peptide>, Failure> {
    let paths = shards(cache);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let scanned: Vec<std::sync::Mutex<Option<Scanned>>> =
        paths.iter().map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1).min(paths.len().max(1)) {
            scope.spawn(|| loop {
                let index = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                let Some(path) = paths.get(index) else {
                    return;
                };
                let events = read(path).and_then(|bytes| scan(&bytes, wanted));
                *scanned[index].lock().expect("a shard's slot") = Some(events);
            });
        }
    });
    let mut found: HashMap<String, Peptide> = HashMap::new();
    for (path, events) in paths.iter().zip(scanned) {
        let events = events
            .into_inner()
            .expect("a shard's slot")
            .expect("every shard is scanned")
            .map_err(|detail| Failure::refused(format!("{}: {detail}", shown(path))))?;
        for event in events {
            match event {
                Event::Found { id, peptide } => {
                    found.entry(id).or_insert(peptide);
                }
                Event::Unreadable { id } if !found.contains_key(&id) => {
                    return Err(Failure::refused(format!(
                        "{}: the variation_effect_feature_cache of {id} is no object",
                        shown(path)
                    )));
                }
                Event::Unreadable { .. } => {}
            }
        }
    }
    Ok(found)
}

/// `Path.name` of a normalised path: its last component, empty for `.` and `/`.
fn path_name(path: &[u8]) -> String {
    let name = path.rsplit(|byte| *byte == b'/').next().unwrap_or(b"");
    if name == b"." {
        String::new()
    } else {
        shown(name)
    }
}

/// Runs the classification. `stdout` takes the summary line.
pub fn run(request: &Request, stdout: &mut dyn Write) -> Result<i32, Failure> {
    // The two dumps are read side by side; a failure of the first is the one reported.
    let (a, b) = std::thread::scope(|scope| {
        let b = scope.spawn(|| load_dump(&request.dump_b));
        let a = load_dump(&request.dump_a);
        (a, b.join().expect("the reader of the second dump"))
    });
    let (a, b) = (a?, b?);
    let mut ids: Vec<&str> = a
        .keys()
        .filter(|id| b.contains_key(*id))
        .map(String::as_str)
        .collect();
    ids.sort_unstable();
    let mut missing: Vec<&str> = a
        .keys()
        .filter(|id| !b.contains_key(*id))
        .chain(b.keys().filter(|id| !a.contains_key(*id)))
        .map(String::as_str)
        .collect();
    missing.sort_unstable();
    if !missing.is_empty() {
        return Err(Failure {
            code: 1,
            line: format!(
                "ERROR: {TAG} {} ids are in one dump only: {}",
                missing.len(),
                repr_str_list(&missing[..missing.len().min(5)])
            ),
        });
    }
    std::fs::create_dir_all(os_path(&request.out_dir)).map_err(|error| {
        Failure::refused(format!("{}: {}", shown(&request.out_dir), os_text(&error)))
    })?;
    let json_peptides = match &request.json_cache {
        Some(cache) => {
            let wanted: HashSet<String> = ids.iter().map(|id| (*id).to_string()).collect();
            Some(json_cache_peptides(cache, &wanted, request.threads)?)
        }
        None => None,
    };

    // The transcripts are classed on the worker threads, each a run of the sorted ids.
    let threads = request.threads.max(1).min(ids.len().max(1));
    let run_len = ids.len().div_ceil(threads).max(1);
    let classed: Vec<Classes> = std::thread::scope(|scope| {
        let workers: Vec<_> = ids
            .chunks(run_len)
            .map(|run| {
                let (a, b) = (&a, &b);
                scope.spawn(move || {
                    run.iter()
                        .map(|id| classify(&a[*id], &b[*id]))
                        .collect::<Vec<Classes>>()
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|worker| worker.join().expect("a classing thread"))
            .collect()
    });

    let (va, vb) = (request.version_a.as_str(), request.version_b.as_str());
    let mut per_class = Counter::default();
    let mut combinations = Counter::default();
    let mut json_matches = Counter::default();
    let mut table =
        String::from("transcript_id\tclasses\tother_keys\tjson_cache_peptide_matches\n");
    for (id, classes) in ids.iter().zip(&classed) {
        // A class is counted where its first transcript stands; the classes of one
        // transcript in the order of their names.
        for class in &classes.classes {
            per_class.add(class);
        }
        let joined = classes
            .classes
            .iter()
            .copied()
            .collect::<Vec<_>>()
            .join("+");
        combinations.add(&joined);
        let mut matches = "";
        if let Some(peptides) = &json_peptides {
            let (pa, pb) = (peptide_of(&a[*id]), peptide_of(&b[*id]));
            let is = |peptide: &Peptide, dump: Option<&str>| match (peptide, dump) {
                (Peptide::Text(text), Some(dump)) => text == dump,
                _ => false,
            };
            matches = match peptides.get(*id) {
                None | Some(Peptide::Absent) => "absent",
                Some(pj) if pa == pb => {
                    if is(pj, pa) {
                        "both"
                    } else {
                        "neither"
                    }
                }
                Some(pj) if is(pj, pa) => va,
                Some(pj) if is(pj, pb) => vb,
                Some(_) => "neither",
            };
            json_matches.add(matches);
        }
        let other_keys = classes
            .other_keys
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(",");
        writeln!(table, "{id}\t{joined}\t{other_keys}\t{matches}").expect("a string takes text");
    }
    let write_file = |name: String, text: &str| -> Result<(), Failure> {
        let path = path_join(&request.out_dir, name.as_bytes());
        std::fs::write(os_path(&path), text)
            .map_err(|error| Failure::refused(format!("{}: {}", shown(&path), os_text(&error))))
    };
    write_file(format!("differing_{va}_{vb}_classes.tsv"), &table)?;

    let mut summary = String::from("{\n  \"version_a\": ");
    push_json_str(&mut summary, va);
    summary.push_str(",\n  \"version_b\": ");
    push_json_str(&mut summary, vb);
    write!(
        summary,
        ",\n  \"transcripts\": {},\n  \"per_class\": {},\n  \"class_combinations\": {}",
        ids.len(),
        per_class.dumps_indented(),
        combinations.dumps_indented()
    )
    .expect("a string takes text");
    if let Some(cache) = &request.json_cache {
        let label = match request.json_cache_label.as_deref() {
            Some(label) if !label.is_empty() => label.to_string(),
            _ => path_name(cache),
        };
        summary.push_str(",\n  \"json_cache\": ");
        push_json_str(&mut summary, &label);
        write!(
            summary,
            ",\n  \"json_cache_peptide_matches\": {}",
            json_matches.dumps_indented()
        )
        .expect("a string takes text");
    }
    summary.push_str("\n}\n");
    write_file(format!("differing_{va}_{vb}_classes.json"), &summary)?;

    let mut counts: Vec<&(String, usize)> = per_class.entries.iter().collect();
    counts.sort();
    let counts = counts
        .iter()
        .map(|(name, count)| format!("{name} {count}"))
        .collect::<Vec<_>>()
        .join(", ");
    let mut line = format!(
        "{TAG} {} transcripts: {counts}; combinations {}",
        ids.len(),
        combinations.dumps()
    );
    if request.json_cache.is_some() {
        write!(
            line,
            "; JSON cache peptide matches {}",
            json_matches.dumps()
        )
        .expect("a string takes text");
    }
    line.push('\n');
    stdout
        .write_all(line.as_bytes())
        .and_then(|()| stdout.flush())
        .map_err(|error| Failure::refused(format!("standard output: {}", os_text(&error))))?;
    Ok(0)
}
