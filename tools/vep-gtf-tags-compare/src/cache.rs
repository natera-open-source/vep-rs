// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cache side: the shards of a JSON cache in the order they are read, the transcripts
//! of one shard, and the code set of the first copy of every stable id.

use crate::gtf::{CodeSets, CODES};
use crate::input;
use crate::json::{Kind, Scanner};
use crate::pyfmt::py_name_order;
use memchr::memchr;
use std::cell::RefCell;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs;
use std::io::{self, Read};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};

/// Set on a transcript whose `attributes` value Python cannot loop over: a number other
/// than zero, or `true`. Such a transcript ends the run where its copy is the first of its
/// stable id, and counts for nothing where an earlier copy was read.
pub const NOT_ITERABLE: u8 = 0x80;

/// The transcripts of a shard in reading order: their unversioned stable ids, none empty,
/// end to end in one buffer, and a code set for each. A shard is three allocations, where an
/// id of its own for every transcript would have every thread grow its heap id by id.
#[derive(Default)]
pub struct Transcripts {
    ids: Vec<u8>,
    ends: Vec<usize>,
    codes: Vec<u8>,
}

impl Transcripts {
    fn push(&mut self, id: &[u8], codes: u8) {
        self.ids.extend_from_slice(id);
        self.ends.push(self.ids.len());
        self.codes.push(codes);
    }

    fn append(&mut self, other: &Transcripts) {
        let base = self.ids.len();
        self.ids.extend_from_slice(&other.ids);
        self.ends.extend(other.ends.iter().map(|end| base + end));
        self.codes.extend_from_slice(&other.codes);
    }

    /// Each transcript's id and code set, in reading order.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], u8)> {
        let mut start = 0;
        self.ends
            .iter()
            .zip(&self.codes)
            .map(move |(&end, &codes)| {
                let id = &self.ids[start..end];
                start = end;
                (id, codes)
            })
    }
}

/// The shard files of a cache in reading order. The walk starts at `transcripts/` under the
/// cache directory, or at the cache directory where that is no directory, and goes top
/// down: the `.json` and `.json.gz` files of a directory sorted by name, then its
/// subdirectories in the order the directory lists them. A symbolic link to a directory is
/// not followed, and a directory that cannot be listed is passed over.
pub fn list_shards(cache_dir: &Path) -> Vec<PathBuf> {
    let transcripts = cache_dir.join("transcripts");
    let root = if transcripts.is_dir() {
        transcripts
    } else {
        cache_dir.to_path_buf()
    };
    let mut shards = Vec::new();
    let mut pending = vec![root];
    while let Some(dir) = pending.pop() {
        let Some((mut files, subdirs)) = list_dir(&dir) else {
            continue;
        };
        files.sort_by(|a, b| py_name_order(a.as_bytes(), b.as_bytes()));
        for name in files {
            let bytes = name.as_bytes();
            if bytes.ends_with(b".json") || bytes.ends_with(b".json.gz") {
                shards.push(dir.join(&name));
            }
        }
        pending.extend(subdirs.into_iter().rev().map(|name| dir.join(name)));
    }
    shards
}

/// The entries of one directory that are no directory, and its subdirectories that are no
/// symbolic link, each in listing order; `None` when the listing fails.
fn list_dir(dir: &Path) -> Option<(Vec<OsString>, Vec<OsString>)> {
    let mut files = Vec::new();
    let mut subdirs = Vec::new();
    for entry in fs::read_dir(dir).ok()? {
        let entry = entry.ok()?;
        match entry.file_type() {
            Ok(kind) if kind.is_dir() => subdirs.push(entry.file_name()),
            // A link to a directory is neither read as a file nor walked.
            Ok(kind) if kind.is_symlink() && entry.path().is_dir() => {}
            _ => files.push(entry.file_name()),
        }
    }
    Some((files, subdirs))
}

/// Why a list of transcripts cannot be read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Fault {
    /// An element Python counts true that is no object.
    NotATranscript,
    /// A `stable_id` that is an array or an object.
    ContainerId,
}

impl Fault {
    fn error(self) -> io::Error {
        io::Error::new(
            io::ErrorKind::InvalidData,
            match self {
                Fault::NotATranscript => "a transcript that is no object",
                Fault::ContainerId => "a stable_id that is an array or an object",
            },
        )
    }
}

#[derive(Default)]
struct List {
    transcripts: Transcripts,
    fault: Option<Fault>,
}

enum Id {
    Text(Vec<u8>),
    Container,
}

struct ShardReader<R> {
    scanner: Scanner<R>,
    name: Vec<u8>,
    code: Vec<u8>,
}

thread_local! {
    /// The window this thread's scanners read through: one buffer a thread, where one a
    /// shard would have every thread map and unmap memory for every file.
    static WINDOW: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

/// Reads one shard file whole and returns its transcripts in reading order.
pub fn read_shard(path: &Path) -> io::Result<Transcripts> {
    let reader = input::open(path)?;
    let mut window = WINDOW.take();
    let transcripts = scan_shard(reader, &mut window);
    WINDOW.set(window);
    transcripts
}

/// The transcripts of one shard document. A top-level array holds them; a top-level object
/// holds them in each member that is an array, in member order, a repeated member name
/// keeping its first place and its last value; any other document holds none. Elements
/// Python counts false are passed over.
fn scan_shard(reader: impl Read, window: &mut Vec<u8>) -> io::Result<Transcripts> {
    let mut shard = ShardReader {
        scanner: Scanner::new(reader, std::mem::take(window)),
        name: Vec::new(),
        code: Vec::new(),
    };
    let lists = shard.document();
    *window = shard.scanner.into_window();
    let mut transcripts = Transcripts::default();
    for list in lists?.into_iter().flatten() {
        if let Some(fault) = list.fault {
            return Err(fault.error());
        }
        if transcripts.ends.is_empty() {
            transcripts = list.transcripts;
        } else {
            transcripts.append(&list.transcripts);
        }
    }
    Ok(transcripts)
}

impl<R: Read> ShardReader<R> {
    /// The lists of the document, in reading order.
    fn document(&mut self) -> io::Result<Vec<Option<List>>> {
        self.scanner.begin()?;
        let lists = match self.scanner.value_start()? {
            b'[' => vec![Some(self.list()?)],
            b'{' => self.member_lists()?,
            _ => {
                self.scanner.skip_value()?;
                Vec::new()
            }
        };
        self.scanner.finish()?;
        Ok(lists)
    }

    /// The value of every member of the top-level object, in the order Python's `dict`
    /// keeps them: a list where the value is an array, `None` otherwise.
    fn member_lists(&mut self) -> io::Result<Vec<Option<List>>> {
        let mut lists = Vec::new();
        let mut places: HashMap<Vec<u8>, usize> = HashMap::new();
        if self.scanner.enter(b'}')? {
            loop {
                // Read into a name of its own: the transcripts of the value reuse
                // `self.name`.
                let mut name = Vec::new();
                self.scanner.key(&mut name)?;
                let list = if self.scanner.value_start()? == b'[' {
                    Some(self.list()?)
                } else {
                    self.scanner.skip_value()?;
                    None
                };
                match places.entry(name) {
                    Entry::Occupied(place) => lists[*place.get()] = list,
                    Entry::Vacant(place) => {
                        place.insert(lists.len());
                        lists.push(list);
                    }
                }
                if !self.scanner.another(b'}')? {
                    break;
                }
            }
        }
        Ok(lists)
    }

    /// The array the scanner stands at, read as a list of transcripts. After a fault the
    /// rest of the array is still held to the grammar.
    fn list(&mut self) -> io::Result<List> {
        let mut list = List::default();
        if !self.scanner.enter(b']')? {
            return Ok(list);
        }
        loop {
            if self.scanner.value_start()? == b'{' {
                // An empty object is false to Python and passed over.
                if let Some((id, codes)) = self.transcript()? {
                    match id {
                        _ if list.fault.is_some() => {}
                        Id::Container => list.fault = Some(Fault::ContainerId),
                        Id::Text(id) if id.is_empty() => {}
                        Id::Text(id) => list.transcripts.push(&id, codes),
                    }
                }
            } else if self.scanner.skip_value()?.truthy && list.fault.is_none() {
                list.fault = Some(Fault::NotATranscript);
            }
            if !self.scanner.another(b']')? {
                return Ok(list);
            }
        }
    }

    /// The object the scanner stands at, read as a transcript: its unversioned id and its
    /// code set, `None` for an empty object. Of a repeated member name the last value
    /// counts, and a transcript without a `stable_id` has the empty id.
    fn transcript(&mut self) -> io::Result<Option<(Id, u8)>> {
        if !self.scanner.enter(b'}')? {
            return Ok(None);
        }
        let mut id = Id::Text(Vec::new());
        let mut codes = 0;
        loop {
            self.name.clear();
            self.scanner.key(&mut self.name)?;
            if self.name == b"stable_id" {
                id = self.stable_id()?;
            } else if self.name == b"attributes" {
                codes = self.attributes()?;
            } else {
                self.scanner.skip_value()?;
            }
            if !self.scanner.another(b'}')? {
                return Ok(Some((id, codes)));
            }
        }
    }

    /// `str(value).split(".")[0]` of a `stable_id` value: a string as it stands, any other
    /// scalar as Python writes it (`None`, `True`, `12`, `1.5`), cut at the first `.`.
    fn stable_id(&mut self) -> io::Result<Id> {
        let mut text = Vec::new();
        match self.scanner.value_start()? {
            b'"' => {
                self.scanner.string(Some(&mut text))?;
            }
            b'[' | b'{' => {
                self.scanner.skip_value()?;
                return Ok(Id::Container);
            }
            _ => {
                self.scanner.scalar(Some(&mut text))?;
            }
        }
        if let Some(dot) = memchr(b'.', &text) {
            text.truncate(dot);
        }
        Ok(Id::Text(text))
    }

    /// The code set of an `attributes` value: of an array, the kept codes of its elements
    /// that are objects; of a string, an object or a value Python counts false, none;
    /// `NOT_ITERABLE` for any other number and for `true`.
    fn attributes(&mut self) -> io::Result<u8> {
        if self.scanner.value_start()? != b'[' {
            let value = self.scanner.skip_value()?;
            let iterable = !value.truthy || !matches!(value.kind, Kind::Number | Kind::Bool);
            return Ok(if iterable { 0 } else { NOT_ITERABLE });
        }
        let mut codes = 0;
        if self.scanner.enter(b']')? {
            loop {
                if self.scanner.value_start()? == b'{' {
                    codes |= self.attribute()?;
                } else {
                    self.scanner.skip_value()?;
                }
                if !self.scanner.another(b']')? {
                    break;
                }
            }
        }
        Ok(codes)
    }

    /// The bit of the kept code an attribute object names in its `code` member, zero when
    /// the last `code` member is no string or names no kept code.
    fn attribute(&mut self) -> io::Result<u8> {
        let mut bit = 0;
        if self.scanner.enter(b'}')? {
            loop {
                self.name.clear();
                self.scanner.key(&mut self.name)?;
                if self.name != b"code" {
                    self.scanner.skip_value()?;
                } else if self.scanner.value_start()? == b'"' {
                    self.code.clear();
                    self.scanner.string(Some(&mut self.code))?;
                    bit = CODES
                        .iter()
                        .position(|code| code.as_bytes() == self.code)
                        .map_or(0, |i| 1 << i);
                } else {
                    self.scanner.skip_value()?;
                    bit = 0;
                }
                if !self.scanner.another(b'}')? {
                    break;
                }
            }
        }
        Ok(bit)
    }
}

/// The code set of every stable id over the shards in reading order. A stable id stored in
/// several shards counts once, with the code set of its first copy.
pub fn merge(shards: &[Transcripts], paths: &[PathBuf]) -> Result<CodeSets, String> {
    let mut sets = CodeSets::new();
    for (transcripts, path) in shards.iter().zip(paths) {
        for (id, codes) in transcripts.iter() {
            if sets.contains_key(id) {
                continue;
            }
            if codes & NOT_ITERABLE != 0 {
                return Err(format!(
                    "{}: the attributes of transcript {} are a number or true",
                    path.display(),
                    String::from_utf8_lossy(id)
                ));
            }
            sets.insert(id.into(), codes);
        }
    }
    Ok(sets)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shard(text: &str) -> io::Result<Vec<(String, u8)>> {
        Ok(scan_shard(text.as_bytes(), &mut Vec::new())?
            .iter()
            .map(|(id, codes)| (String::from_utf8_lossy(id).into_owned(), codes))
            .collect())
    }

    fn ids(text: &str) -> Vec<String> {
        shard(text).unwrap().into_iter().map(|(id, _)| id).collect()
    }

    #[test]
    fn a_list_shard_yields_ids_and_kept_codes_in_order() {
        let text = r#"[
            {"attributes": [{"code": "gencode_basic", "name": "n", "value": "v"},
                            {"code": "TSL", "value": "tsl1"},
                            {"code": "gencode_primary"}],
             "exons": [{"stable_id": "ENSE1", "start": 1}], "stable_id": "ENST1", "version": 2},
            {"stable_id": "ENST2.7", "attributes": [{"code": "cds_start_NF"}, {"code": "cds_end_NF"}]},
            {"stable_id": "ENST3"}
        ]"#;
        assert_eq!(
            shard(text).unwrap(),
            [
                ("ENST1".to_string(), 1 | 2),
                ("ENST2".to_string(), 4 | 8),
                ("ENST3".to_string(), 0),
            ]
        );
    }

    #[test]
    fn an_object_shard_yields_the_transcripts_of_its_array_members() {
        let text = r#"{"meta": {"list": [1, 2]}, "a": [{"stable_id": "FIRST"}], "n": 3,
                       "b": [{"stable_id": "OTHER"}, {"stable_id": "THIRD"}], "s": "text"}"#;
        assert_eq!(ids(text), ["FIRST", "OTHER", "THIRD"]);
        // A repeated name keeps its first place and takes its last value, whatever the
        // value it replaces held.
        let text = r#"{"a": [7, {"stable_id": "dropped"}], "b": [{"stable_id": "OTHER"}],
                       "a": [{"stable_id": "LAST"}]}"#;
        assert_eq!(ids(text), ["LAST", "OTHER"]);
        let text = r#"{"a": [{"stable_id": "FIRST"}], "b": [{"stable_id": "OTHER"}], "a": 5}"#;
        assert_eq!(ids(text), ["OTHER"]);
        for text in ["null", "5", r#""text""#, "{}", "[]", r#"{"a": {"b": [9]}}"#] {
            assert!(ids(text).is_empty(), "{text}");
        }
    }

    #[test]
    fn elements_python_counts_false_are_passed_over_and_other_scalars_refused() {
        let text = r#"[null, 0, -0, 0.0, 1e-400, "", [], {}, false, {"stable_id": "T1"}]"#;
        assert_eq!(ids(text), ["T1"]);
        for element in ["5", "true", r#""x""#, "[0]", "NaN", "1e-300"] {
            let text = format!(r#"[{{"stable_id": "T1"}}, {element}, {{"stable_id": "T2"}}]"#);
            assert!(shard(&text).is_err(), "{text}");
            let text = format!(r#"{{"a": [{element}], "a": [{{"stable_id": "T1"}}]}}"#);
            assert_eq!(ids(&text), ["T1"], "{text}");
        }
    }

    #[test]
    fn stable_ids_are_what_str_and_split_make_of_the_value() {
        let cases = [
            (r#""ENST1.4""#, Some("ENST1")),
            (r#""V.2.3""#, Some("V")),
            (r#""""#, None),
            (r#"".5""#, None),
            ("null", Some("None")),
            ("true", Some("True")),
            ("false", Some("False")),
            ("12", Some("12")),
            ("-0", Some("0")),
            ("1.5", Some("1")),
            ("0.5", Some("0")),
            ("100.0", Some("100")),
            ("1e5", Some("100000")),
            ("1e22", Some("1e+22")),
            ("1.0e-5", Some("1e-05")),
            ("NaN", Some("nan")),
            ("-Infinity", Some("-inf")),
        ];
        for (value, expected) in cases {
            let text = format!(r#"[{{"x": 1, "stable_id": {value}}}]"#);
            assert_eq!(
                ids(&text),
                expected.map(String::from).into_iter().collect::<Vec<_>>(),
                "{value}"
            );
        }
        // No stable_id is the empty id; of two the last counts.
        assert!(ids(r#"[{"x": 1}]"#).is_empty());
        assert_eq!(ids(r#"[{"stable_id": "A", "stable_id": "B"}]"#), ["B"]);
        assert!(ids(r#"[{"stable_id": "A", "stable_id": ""}]"#).is_empty());
    }

    #[test]
    fn a_stable_id_that_is_a_container_is_refused() {
        assert!(shard(r#"[{"stable_id": ["A"]}]"#).is_err());
        assert!(shard(r#"[{"stable_id": {}}]"#).is_err());
        assert_eq!(ids(r#"[{"stable_id": [], "stable_id": "A"}]"#), ["A"]);
    }

    #[test]
    fn attributes_count_only_objects_of_an_array_with_a_kept_code_string() {
        let codes = |attributes: &str| {
            let text = format!(r#"[{{"stable_id": "T", "attributes": {attributes}}}]"#);
            shard(&text).unwrap()[0].1
        };
        assert_eq!(codes("null"), 0);
        assert_eq!(codes("[]"), 0);
        assert_eq!(codes("{}"), 0);
        assert_eq!(codes("false"), 0);
        assert_eq!(codes("0"), 0);
        assert_eq!(codes("0.0"), 0);
        assert_eq!(codes(r#""gencode_basic""#), 0);
        assert_eq!(codes(r#"{"code": "gencode_basic"}"#), 0);
        assert_eq!(codes("7"), NOT_ITERABLE);
        assert_eq!(codes("true"), NOT_ITERABLE);
        assert_eq!(codes("NaN"), NOT_ITERABLE);
        assert_eq!(codes("-1.5"), NOT_ITERABLE);
        assert_eq!(
            codes(
                r#"["gencode_basic", ["gencode_basic"], null, 5, {"code": 5}, {"code": null},
                    {"name": "x"}, {"code": "basic"}, {"code": "GENCODE_BASIC"},
                    {"code": ["cds_start_NF"]}, {"code": "gencode_primary", "code": "other"},
                    {"code": "other", "code": "cds_end_NF"}]"#
            ),
            8
        );
        assert_eq!(
            codes(
                r#"[{"code": "gencode_basic"}, {"code": "cds_start_NF"}, {"code": "gencode_basic"}]"#
            ),
            1 | 4
        );
        let text =
            r#"[{"stable_id": "T", "attributes": [{"code": "cds_end_NF"}], "attributes": null}]"#;
        assert_eq!(shard(text).unwrap()[0].1, 0);
    }

    #[test]
    fn the_first_copy_of_a_stable_id_counts() {
        let shard = |transcripts: &[(&str, u8)]| {
            let mut shard = Transcripts::default();
            for (id, codes) in transcripts {
                shard.push(id.as_bytes(), *codes);
            }
            shard
        };
        let paths = [PathBuf::from("a.json"), PathBuf::from("b.json")];
        let shards = [
            shard(&[("T1", 1), ("T2", 2), ("T1", 8)]),
            shard(&[("T2", NOT_ITERABLE), ("T3", 4)]),
        ];
        let sets = merge(&shards, &paths).unwrap();
        assert_eq!(sets.len(), 3);
        assert_eq!(sets[&b"T1"[..]], 1);
        assert_eq!(sets[&b"T2"[..]], 2);
        assert_eq!(sets[&b"T3"[..]], 4);
        let shards = [shard(&[("T1", 1)]), shard(&[("T2", NOT_ITERABLE)])];
        let error = merge(&shards, &paths).unwrap_err();
        assert!(error.starts_with("b.json: "), "{error}");
    }

    #[test]
    fn shards_are_listed_top_down_with_each_directory_s_files_sorted() {
        let root =
            std::env::temp_dir().join(format!("vep-gtf-tags-compare-walk-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        for dir in ["flat/sub.json", "nested/transcripts/21", "nested/other"] {
            fs::create_dir_all(root.join(dir)).unwrap();
        }
        for file in [
            "flat/b.json",
            "flat/a.json.gz",
            "flat/notes.txt",
            "flat/a.json.bak",
            "flat/.json",
            "flat/sub.json/10-11.json",
            "flat/sub.json/1-2.json",
            "nested/info.json",
            "nested/other/x.json",
            "nested/transcripts/21/2.json",
            "nested/transcripts/21/1.json.gz",
        ] {
            fs::write(root.join(file), b"[]").unwrap();
        }
        std::os::unix::fs::symlink(root.join("nested/other"), root.join("flat/linked")).unwrap();
        std::os::unix::fs::symlink(root.join("flat/b.json"), root.join("flat/c.json")).unwrap();
        let listed = |cache: &str| -> Vec<String> {
            let cache = root.join(cache);
            list_shards(&cache)
                .iter()
                .map(|p| p.strip_prefix(&cache).unwrap().display().to_string())
                .collect()
        };
        assert_eq!(
            listed("flat"),
            [
                ".json",
                "a.json.gz",
                "b.json",
                "c.json",
                "sub.json/1-2.json",
                "sub.json/10-11.json"
            ]
        );
        assert_eq!(
            listed("nested"),
            ["transcripts/21/1.json.gz", "transcripts/21/2.json"]
        );
        assert!(listed("missing").is_empty());
        fs::remove_dir_all(&root).unwrap();
    }
}
