// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript shards of a candidate JSON cache.
//!
//! A cache files its transcripts under `transcripts/<chromosome>/<shard>.json` or
//! `<shard>.json.gz`. Listing reads no shard: a shard is read, and a malformed one ends
//! the run, when a look-up first needs it.

use crate::argparse::{path_text, repr};
use crate::error::{Error, Result, Uncaught};
use crate::pyjson::{self, Value};
use crate::report::strerror;
use std::ffi::OsString;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::panic;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Bounds the workers that read shards; unset, they are the cores of the host.
pub const THREADS_VARIABLE: &str = "VEP_FORMATS_COMPARE_THREADS";

/// The workers a reader of many shards runs on.
pub fn threads() -> usize {
    std::env::var(THREADS_VARIABLE)
        .ok()
        .and_then(|count| count.parse::<usize>().ok())
        .filter(|count| *count >= 1)
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |cores| cores.get()))
}

/// `base / name` of `pathlib`: an absolute `name` stands alone, an empty part and a `.`
/// part are dropped, `..` is kept, and exactly two leading `/` stay two.
pub fn joined(base: &Path, name: &str) -> PathBuf {
    let base = base.as_os_str().as_bytes();
    let mut raw = Vec::with_capacity(base.len() + 1 + name.len());
    if !name.starts_with('/') {
        raw.extend_from_slice(base);
        if !raw.is_empty() && !raw.ends_with(b"/") {
            raw.push(b'/');
        }
    }
    raw.extend_from_slice(name.as_bytes());
    let root: &[u8] = match raw.iter().take_while(|byte| **byte == b'/').count() {
        0 => b"",
        2 => b"//",
        _ => b"/",
    };
    let mut text = root.to_vec();
    for part in raw
        .split(|byte| *byte == b'/')
        .filter(|part| !part.is_empty() && *part != b".")
    {
        if text.len() > root.len() {
            text.push(b'/');
        }
        text.extend_from_slice(part);
    }
    if text.is_empty() {
        text.push(b'.');
    }
    PathBuf::from(OsString::from_vec(text))
}

/// The names in `directory` in byte order.
fn names(directory: &Path) -> std::io::Result<Vec<OsString>> {
    let mut names = Vec::new();
    for entry in std::fs::read_dir(directory)? {
        names.push(entry?.file_name());
    }
    names.sort();
    Ok(names)
}

fn suffixed<'a>(
    directory: &'a Path,
    names: &'a [OsString],
    suffix: &'a str,
) -> impl Iterator<Item = PathBuf> + 'a {
    names
        .iter()
        .filter(move |name| name.as_bytes().ends_with(suffix.as_bytes()))
        .map(move |name| directory.join(name))
}

/// The shards of `root`'s directory `chromosome`, or without one of every entry of `root`
/// in the order of the entries' names: a directory's `*.json` in the order of their names,
/// then its `*.json.gz`. A name that opens with a dot is a name like any other, a
/// directory of the suffix is listed like a file, and a directory that is absent or cannot
/// be listed has no shards. Ends the run where `root` itself cannot be listed.
pub fn shards(root: &Path, chromosome: Option<&str>) -> Result<Vec<PathBuf>> {
    let directories = match chromosome {
        Some(name) => vec![joined(root, name)],
        None => match names(root) {
            Ok(entries) => entries.iter().map(|entry| root.join(entry)).collect(),
            Err(error) => crate::uncaught!(
                "[Errno {}] {}: {}",
                error.raw_os_error().unwrap_or(0),
                strerror(&error),
                repr(&path_text(root))
            ),
        },
    };
    let mut found = Vec::new();
    for directory in &directories {
        let listed = names(directory).unwrap_or_default();
        found.extend(suffixed(directory, &listed, ".json"));
        found.extend(suffixed(directory, &listed, ".json.gz"));
    }
    Ok(found)
}

/// Every `*/*.json` of `root` in the order of the directories' names and then the shards'
/// names, then every `*/*.json.gz` in that order: the order `NearestSites.from_cache`
/// reads a cache in. A `root` that cannot be listed has no shards.
pub fn shards_by_suffix(root: &Path) -> Vec<PathBuf> {
    let directories: Vec<(PathBuf, Vec<OsString>)> = names(root)
        .unwrap_or_default()
        .iter()
        .map(|entry| root.join(entry))
        .filter(|directory| directory.is_dir())
        .map(|directory| {
            let listed = names(&directory).unwrap_or_default();
            (directory, listed)
        })
        .collect();
    let mut found = Vec::new();
    for suffix in [".json", ".json.gz"] {
        for (directory, listed) in &directories {
            found.extend(suffixed(directory, listed, suffix));
        }
    }
    found
}

/// `NearestSites.shard_transcripts` of a parsed shard the caller gives up: the transcript
/// objects of a bare array, or of a range-keyed map of arrays, a `null` in either dropped.
pub fn transcripts_of(shard: &Path, parsed: Value) -> Result<Vec<Value>> {
    let items = match parsed.into_dict() {
        Ok(ranges) => {
            let mut items = Vec::new();
            if !ranges.values().all(|wrapped| wrapped.as_list().is_some()) {
                return Err(Error::comparator(format!(
                    "{}: a range-keyed map value is not a transcript array",
                    path_text(shard)
                )));
            }
            for (_, wrapped) in ranges {
                items.extend(wrapped.into_list().unwrap_or_default());
            }
            items
        }
        Err(parsed) => match parsed.into_list() {
            Ok(items) => items,
            Err(_) => {
                return Err(Error::comparator(format!(
                    "{} is neither a transcript array nor a range-keyed map",
                    path_text(shard)
                )))
            }
        },
    };
    let transcripts: Vec<Value> = items
        .into_iter()
        .filter(|transcript| !matches!(transcript, Value::Null))
        .collect();
    if !transcripts
        .iter()
        .all(|transcript| transcript.as_dict().is_some())
    {
        return Err(Error::comparator(format!(
            "{}: a transcript is not an object",
            path_text(shard)
        )));
    }
    Ok(transcripts)
}

/// The transcript objects of one shard: `json.loads(read_text(shard))` handed to
/// `NearestSites.shard_transcripts`. A shard that is not JSON ends the run with
/// `<shard> is not JSON: <the reader's words>`.
pub fn shard_transcripts(shard: &Path) -> Result<Vec<Value>> {
    let text = crate::input::read_text(shard)?;
    match pyjson::loads(&text) {
        Ok(parsed) => transcripts_of(shard, parsed),
        Err(error) if error.is_value_error() => Err(Error::comparator(format!(
            "{} is not JSON: {error}",
            path_text(shard)
        ))),
        Err(error) => crate::uncaught!("{}: {error}", path_text(shard)),
    }
}

/// How the reading of one shard ended on a worker.
enum Ended<T> {
    Read(T),
    Refused(Error),
    Uncaught(Uncaught),
}

/// `read(shard, shard_transcripts(shard))` of every shard, on `threads` workers, in the
/// shards' order. The run ends as it ends where the shards are read one after another: on
/// the first shard in that order whose reading fails, whatever a later shard holds. No
/// more than `threads` shards are held parsed at once.
pub fn map_shards<T, F>(shards: &[PathBuf], threads: usize, read: F) -> Result<Vec<T>>
where
    T: Send,
    F: Fn(&Path, Vec<Value>) -> Result<T> + Sync,
{
    let workers = threads.clamp(1, shards.len().max(1));
    if workers == 1 {
        return shards
            .iter()
            .map(|shard| read(shard, shard_transcripts(shard)?))
            .collect();
    }
    let next = AtomicUsize::new(0);
    let first_failure = AtomicUsize::new(usize::MAX);
    let one = |shard: &Path| -> Ended<T> {
        let ended = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            read(shard, shard_transcripts(shard)?)
        }));
        match ended {
            Ok(Ok(value)) => Ended::Read(value),
            Ok(Err(error)) => Ended::Refused(error),
            Err(payload) => match payload.downcast::<Uncaught>() {
                Ok(uncaught) => Ended::Uncaught(*uncaught),
                Err(payload) => panic::resume_unwind(payload),
            },
        }
    };
    let mut ended: Vec<Option<Ended<T>>> = Vec::new();
    ended.resize_with(shards.len(), || None);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut done = Vec::new();
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= shards.len() || index > first_failure.load(Ordering::Relaxed) {
                            break;
                        }
                        let result = one(&shards[index]);
                        if !matches!(result, Ended::Read(_)) {
                            first_failure.fetch_min(index, Ordering::Relaxed);
                        }
                        done.push((index, result));
                    }
                    done
                })
            })
            .collect();
        for handle in handles {
            match handle.join() {
                Ok(done) => {
                    for (index, result) in done {
                        ended[index] = Some(result);
                    }
                }
                Err(payload) => panic::resume_unwind(payload),
            }
        }
    });
    let mut values = Vec::with_capacity(shards.len());
    for result in ended {
        match result {
            Some(Ended::Read(value)) => values.push(value),
            Some(Ended::Refused(error)) => return Err(error),
            Some(Ended::Uncaught(Uncaught(what))) => crate::uncaught!("{what}"),
            // A worker passes over a shard only after an earlier one failed, and the
            // walk returns at that earlier one.
            None => unreachable!("a shard before the first failure was read"),
        }
    }
    Ok(values)
}

/// `shard_transcripts` of every shard, read on `threads` workers, in the shards' order,
/// or the failure of the first shard in that order that has one. Every shard is held
/// parsed at once, where `map_shards` holds what its caller keeps of each; a caller that
/// can refuse a transcript does so inside `map_shards`, which keeps the order of the
/// failures.
pub fn shards_transcripts(shards: &[PathBuf], threads: usize) -> Result<Vec<Vec<Value>>> {
    map_shards(shards, threads, |_, transcripts| Ok(transcripts))
}

#[cfg(test)]
pub(crate) mod testing {
    use std::path::{Path, PathBuf};

    /// An empty directory of this process under the system's temporary one.
    pub(crate) fn scratch(name: &str) -> PathBuf {
        let directory =
            std::env::temp_dir().join(format!("vep-formats-compare-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&directory);
        std::fs::create_dir_all(&directory).unwrap();
        directory
    }

    pub(crate) fn plant(root: &Path, name: &str, text: &str) {
        let path = root.join(name);
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(path, text).unwrap();
    }

    /// A cache of one shard, `transcripts/1/1-1000000.json`, holding `shard`.
    pub(crate) fn one_shard_cache(root: &Path, shard: &str) -> PathBuf {
        plant(root, "transcripts/1/1-1000000.json", shard);
        root.to_path_buf()
    }
}

#[cfg(test)]
mod tests {
    use super::testing::{plant, scratch};
    use super::*;

    #[test]
    fn a_joined_path_reads_as_the_path_library_writes_it() {
        let joined = |base: &str, name: &str| joined(Path::new(base), name).display().to_string();
        assert_eq!(joined("cache/transcripts", "21"), "cache/transcripts/21");
        assert_eq!(joined("cache/transcripts", ""), "cache/transcripts");
        assert_eq!(joined("cache/transcripts", "."), "cache/transcripts");
        assert_eq!(joined("cache/transcripts", ".."), "cache/transcripts/..");
        assert_eq!(
            joined("cache/transcripts", "a//b/./c/"),
            "cache/transcripts/a/b/c"
        );
        assert_eq!(joined("cache/transcripts", "/etc"), "/etc");
        assert_eq!(joined("cache/transcripts", "//etc"), "//etc");
        assert_eq!(joined("cache/transcripts", "///etc"), "/etc");
        assert_eq!(joined(".", "transcripts"), "transcripts");
        assert_eq!(joined(".", ""), ".");
        assert_eq!(joined("/", "x"), "/x");
    }

    #[test]
    fn shards_are_listed_by_name_the_plain_ones_of_a_directory_first() {
        let root = scratch("shards");
        for name in [
            "a/x.json",
            "a/.json",
            "a/.h.json",
            "a/X.json",
            "a/x.JSON",
            "a/x.json.gz",
            "a/y.json.gz.json",
            "a/x.jsonx",
            "a-b/1.json",
            "a.b/1.json",
            ".hid/1.json",
            "B/1.json.gz",
            "file.json",
        ] {
            plant(&root, name, "[]");
        }
        std::fs::create_dir(root.join("a/dir.json")).unwrap();
        let relative = |found: Vec<PathBuf>| -> Vec<String> {
            found
                .iter()
                .map(|path| path.strip_prefix(&root).unwrap().display().to_string())
                .collect()
        };
        assert_eq!(
            relative(shards(&root, Some("a")).unwrap()),
            [
                "a/.h.json",
                "a/.json",
                "a/X.json",
                "a/dir.json",
                "a/x.json",
                "a/y.json.gz.json",
                "a/x.json.gz"
            ]
        );
        assert!(shards(&root, Some("missing")).unwrap().is_empty());
        assert!(shards(&root, Some("file.json")).unwrap().is_empty());
        assert_eq!(relative(shards(&root, Some("")).unwrap()), ["file.json"]);
        let every = [
            ".hid/1.json",
            "B/1.json.gz",
            "a/.h.json",
            "a/.json",
            "a/X.json",
            "a/dir.json",
            "a/x.json",
            "a/y.json.gz.json",
            "a/x.json.gz",
            "a-b/1.json",
            "a.b/1.json",
        ];
        assert_eq!(relative(shards(&root, None).unwrap()), every);
        let mut plain_first: Vec<&str> = every
            .iter()
            .copied()
            .filter(|name| !name.ends_with(".gz"))
            .collect();
        plain_first.extend(["B/1.json.gz", "a/x.json.gz"]);
        assert_eq!(relative(shards_by_suffix(&root)), plain_first);
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_parsed_shard_is_an_array_or_a_map_of_arrays_of_objects() {
        let shard = Path::new("c/transcripts/1/1-5.json");
        let of = |text: &str| transcripts_of(shard, pyjson::loads(text).unwrap());
        let refused = |text: &str| of(text).unwrap_err().to_string();
        assert_eq!(of(r#"[null, {"a": 1}, null, {"b": 2}]"#).unwrap().len(), 2);
        let wrapped =
            of(r#"{"1-5": [null, {"a": 1}], "6-9": [{"b": 2}, null, {"c": 3}]}"#).unwrap();
        assert_eq!(
            wrapped
                .iter()
                .map(|transcript| transcript
                    .as_dict()
                    .unwrap()
                    .keys()
                    .next()
                    .unwrap()
                    .as_str())
                .collect::<Vec<_>>(),
            ["a", "b", "c"]
        );
        assert!(of("[]").unwrap().is_empty() && of("{}").unwrap().is_empty());
        assert_eq!(
            refused(r#"{"1-5": [{"a": 1}], "6-9": {"b": 2}}"#),
            "c/transcripts/1/1-5.json: a range-keyed map value is not a transcript array"
        );
        assert_eq!(
            refused("5"),
            "c/transcripts/1/1-5.json is neither a transcript array nor a range-keyed map"
        );
        assert_eq!(
            refused(r#"[{"a": 1}, 5]"#),
            "c/transcripts/1/1-5.json: a transcript is not an object"
        );
        assert_eq!(
            refused(r#"{"1-5": [[1]]}"#),
            "c/transcripts/1/1-5.json: a transcript is not an object"
        );
    }

    #[test]
    fn many_shards_end_on_the_first_failure_in_their_order() {
        let root = scratch("map-shards");
        for at in 0..12 {
            plant(
                &root,
                &format!("1/{at:02}.json"),
                &format!(r#"[{{"n": {at}}}]"#),
            );
        }
        let count = |shard: &Path, transcripts: Vec<Value>| -> Result<String> {
            Ok(format!(
                "{} {}",
                shard.file_name().unwrap().to_string_lossy(),
                transcripts.len()
            ))
        };
        let all = shards(&root, Some("1")).unwrap();
        for threads in [1, 3, 32] {
            let read = map_shards(&all, threads, count).unwrap();
            assert_eq!(read.len(), 12);
            assert_eq!(read[0], "00.json 1");
            assert_eq!(read[11], "11.json 1");
        }
        plant(&root, "1/09.json", "[");
        plant(&root, "1/04.json", "5");
        plant(&root, "1/07.json", "{");
        for threads in [1, 3, 32] {
            let refused = map_shards(&all, threads, count).unwrap_err().to_string();
            assert!(
                refused.ends_with("1/04.json is neither a transcript array nor a range-keyed map"),
                "{refused}"
            );
        }
        std::fs::remove_dir_all(root).unwrap();
    }
}
