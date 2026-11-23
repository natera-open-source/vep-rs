// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Paths as Python resolves and lists them: `Path(arg).resolve()`, the transcript shards
//! `Path.rglob` finds under `transcripts/`, the variation files `Path.glob` finds under
//! `variations/`, and the two sides paired on their names in the order `sorted()` gives.

use crate::pyfmt::py_name_order;
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs;
use std::io;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};

fn as_path(bytes: &[u8]) -> &Path {
    Path::new(OsStr::from_bytes(bytes))
}

/// `Path(arg).resolve()` without `strict`, as `os.path.realpath` walks a path: every
/// symbolic link followed, `.` and empty parts dropped, `..` dropping the part before it,
/// a part that does not exist kept as written, a link that leads back to itself kept as
/// the link. Fails only where the working directory cannot be read.
pub fn py_resolve(arg: &OsStr) -> io::Result<PathBuf> {
    let filename = arg.as_bytes();
    let parts = |text: &[u8]| -> Vec<Option<Vec<u8>>> {
        text.split(|b| *b == b'/')
            .rev()
            .map(|part| Some(part.to_vec()))
            .collect()
    };
    // The unread parts, the next one last. `None` marks a link's target as walked: the
    // link's own path lies under it.
    let mut rest = parts(filename);
    let mut unread = rest.len();
    let mut path: Vec<u8> = if filename.starts_with(b"/") {
        b"/".to_vec()
    } else {
        std::env::current_dir()?.into_os_string().into_vec()
    };
    // A link's path to what it resolved to, `None` while its target is being walked.
    let mut seen: HashMap<Vec<u8>, Option<Vec<u8>>> = HashMap::new();
    while unread > 0 {
        let Some(name) = rest.pop().expect("an unread part") else {
            let link = rest.pop().flatten().expect("the link under its marker");
            seen.insert(link, Some(path.clone()));
            continue;
        };
        unread -= 1;
        if name.is_empty() || name == b"." {
            continue;
        }
        if name == b".." {
            let cut = path
                .iter()
                .rposition(|b| *b == b'/')
                .expect("an absolute path");
            path.truncate(cut.max(1));
            continue;
        }
        let mut longer = path.clone();
        if path != b"/" {
            longer.push(b'/');
        }
        longer.extend_from_slice(&name);
        let is_link =
            fs::symlink_metadata(as_path(&longer)).is_ok_and(|meta| meta.file_type().is_symlink());
        if !is_link {
            path = longer;
            continue;
        }
        if let Some(resolved) = seen.get(&longer) {
            path = resolved.clone().unwrap_or(longer);
            continue;
        }
        let Ok(target) = fs::read_link(as_path(&longer)) else {
            path = longer;
            continue;
        };
        let target = target.into_os_string().into_vec();
        if target.starts_with(b"/") {
            path = b"/".to_vec();
        }
        seen.insert(longer.clone(), None);
        rest.push(Some(longer));
        rest.push(None);
        let target_parts = parts(&target);
        unread += target_parts.len();
        rest.extend(target_parts);
    }
    Ok(PathBuf::from(OsStr::from_bytes(&path)))
}

/// Collects the entries under `dir` whose name ends in `.json` or `.json.gz`, at any
/// depth, as `rglob("*.json")` and `rglob("*.json.gz")` list them: an entry of any type
/// counts, a directory is entered unless it is a symbolic link, and a directory that
/// cannot be listed holds nothing. `relative` is the path of `dir` under the listing's
/// root, with a `/` at its end unless empty.
fn collect_shards(
    dir: &Path,
    relative: &mut Vec<u8>,
    found: &mut HashMap<Vec<u8>, (PathBuf, bool)>,
) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let name = entry.file_name();
        let name = name.as_bytes();
        let gzipped = name.ends_with(b".json.gz");
        if gzipped || name.ends_with(b".json") {
            let mut key = relative.clone();
            key.extend_from_slice(name);
            if gzipped {
                key.truncate(key.len() - b".gz".len());
            }
            // A name both a plain and a gzipped file hold is read from the gzipped one.
            match found.get(&key) {
                Some((_, true)) => {}
                _ => {
                    found.insert(key, (entry.path(), gzipped));
                }
            }
        }
        if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
            let length = relative.len();
            relative.extend_from_slice(name);
            relative.push(b'/');
            collect_shards(&entry.path(), relative, found);
            relative.truncate(length);
        }
    }
}

/// The transcript shards of a cache: the path of each under `transcripts/`, without a
/// `.gz` suffix, to the file read for it. Empty where `transcripts` is no directory.
fn shards(cache: &Path) -> HashMap<Vec<u8>, PathBuf> {
    let base = cache.join("transcripts");
    let mut found = HashMap::new();
    if fs::metadata(&base).is_ok_and(|meta| meta.is_dir()) {
        collect_shards(&base, &mut Vec::new(), &mut found);
    }
    found
        .into_iter()
        .map(|(key, (path, _))| (key, path))
        .collect()
}

/// The variation files of a cache: each chromosome directory under `variations/` that
/// holds an entry `all_vars.gz`, to that entry. Empty where `variations` is no directory.
fn variation_files(cache: &Path) -> HashMap<Vec<u8>, PathBuf> {
    let base = cache.join("variations");
    let mut found = HashMap::new();
    let Ok(entries) = fs::read_dir(&base) else {
        return found;
    };
    for entry in entries.flatten() {
        let directory = entry.path();
        if !fs::metadata(&directory).is_ok_and(|meta| meta.is_dir()) {
            continue;
        }
        let file = directory.join("all_vars.gz");
        if fs::symlink_metadata(&file).is_ok() {
            found.insert(entry.file_name().into_vec(), file);
        }
    }
    found
}

/// One name both sides are asked for: a shard's path under `transcripts/` or a
/// chromosome, with the file of each side that has one.
pub struct Pair {
    pub name: Vec<u8>,
    pub baseline: Option<PathBuf>,
    pub builder: Option<PathBuf>,
}

fn pairs(
    mut baseline: HashMap<Vec<u8>, PathBuf>,
    mut builder: HashMap<Vec<u8>, PathBuf>,
) -> Vec<Pair> {
    let mut names: Vec<Vec<u8>> = baseline.keys().chain(builder.keys()).cloned().collect();
    names.sort_unstable_by(|a, b| py_name_order(a, b));
    names.dedup();
    names
        .into_iter()
        .map(|name| Pair {
            baseline: baseline.remove(&name),
            builder: builder.remove(&name),
            name,
        })
        .collect()
}

/// The transcript shards of the two caches paired on their names, in sorted order.
pub fn shard_pairs(baseline: &Path, builder: &Path) -> Vec<Pair> {
    pairs(shards(baseline), shards(builder))
}

/// The variation files of the two caches paired on their chromosomes, in sorted order.
pub fn chromosome_pairs(baseline: &Path, builder: &Path) -> Vec<Pair> {
    pairs(variation_files(baseline), variation_files(builder))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::symlink;

    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let path = std::env::temp_dir()
                .join(format!("vep-cache-diff-paths-{}-{tag}", std::process::id()));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(&path).unwrap();
            Scratch(fs::canonicalize(path).unwrap())
        }

        fn file(&self, relative: &str) -> PathBuf {
            let path = self.0.join(relative);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(&path, b"[]").unwrap();
            path
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn resolved(arg: impl AsRef<OsStr>) -> PathBuf {
        py_resolve(arg.as_ref()).unwrap()
    }

    #[test]
    fn resolve_follows_links_and_keeps_parts_that_do_not_exist() {
        let scratch = Scratch::new("resolve");
        let root = &scratch.0;
        scratch.file("real/inner/file");
        symlink("real", root.join("link")).unwrap();
        symlink(root.join("real/inner"), root.join("absolute")).unwrap();
        symlink("loop", root.join("loop")).unwrap();
        symlink("missing/target", root.join("dangling")).unwrap();
        symlink("../real/./inner//", root.join("real/back")).unwrap();
        let at = |tail: &str| root.join(tail);
        assert_eq!(resolved(at("link/inner/file")), at("real/inner/file"));
        assert_eq!(resolved(at("absolute/file")), at("real/inner/file"));
        assert_eq!(resolved(at("real/back/file")), at("real/inner/file"));
        assert_eq!(resolved(at("link/../real/./inner//")), at("real/inner"));
        assert_eq!(resolved(at("link/inner/..")), at("real"));
        assert_eq!(resolved(at("nothing/here/../there")), at("nothing/there"));
        assert_eq!(resolved(at("dangling/x")), at("missing/target/x"));
        assert_eq!(resolved(at("loop/x")), at("loop/x"));
        assert_eq!(resolved("/"), Path::new("/"));
        assert_eq!(resolved("//"), Path::new("/"));
        assert_eq!(resolved("/.."), Path::new("/"));
        assert_eq!(resolved("/../x/../y/"), Path::new("/y"));
        let here = fs::canonicalize(".").unwrap();
        assert_eq!(resolved(""), here);
        assert_eq!(resolved("."), here);
        assert_eq!(resolved("a//b/./c/.."), here.join("a/b"));
    }

    fn names(pairs: &[Pair]) -> Vec<(String, Option<String>, Option<String>)> {
        let file = |path: &Option<PathBuf>| {
            path.as_ref()
                .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        };
        pairs
            .iter()
            .map(|pair| {
                (
                    String::from_utf8_lossy(&pair.name).into_owned(),
                    file(&pair.baseline),
                    file(&pair.builder),
                )
            })
            .collect()
    }

    #[test]
    fn shards_are_listed_as_rglob_lists_them_and_paired_on_their_names() {
        let scratch = Scratch::new("shards");
        let root = &scratch.0;
        for name in [
            "a/transcripts/21/1-2.json",
            "a/transcripts/21/1-2.json.gz",
            "a/transcripts/21/10-11.json",
            "a/transcripts/21/.json",
            "a/transcripts/.hidden/h.json",
            "a/transcripts/top.json",
            "a/transcripts/deep/er/d.json.gz",
            "a/transcripts/dir.json/in.json",
            "a/transcripts/real/r.json",
            "a/transcripts/21/x.JSON",
            "a/transcripts/21/y.json.GZ",
            "a/transcripts/21/z.jsonx",
            "a/transcripts/21/w.json.gz.bak",
            "b/transcripts/21/1-2.json",
            "b/transcripts/21/3-4.json.gz",
        ] {
            scratch.file(name);
        }
        symlink(
            root.join("a/transcripts/real"),
            root.join("a/transcripts/linked"),
        )
        .unwrap();
        symlink("nowhere", root.join("a/transcripts/21/broken.json")).unwrap();
        fs::create_dir(root.join("c")).unwrap();
        let both = shard_pairs(&root.join("a"), &root.join("b"));
        let s = |text: &str| Some(text.to_string());
        assert_eq!(
            names(&both),
            [
                (".hidden/h.json".to_string(), s("h.json"), None),
                ("21/.json".to_string(), s(".json"), None),
                ("21/1-2.json".to_string(), s("1-2.json.gz"), s("1-2.json")),
                ("21/10-11.json".to_string(), s("10-11.json"), None),
                ("21/3-4.json".to_string(), None, s("3-4.json.gz")),
                ("21/broken.json".to_string(), s("broken.json"), None),
                ("deep/er/d.json".to_string(), s("d.json.gz"), None),
                ("dir.json".to_string(), s("dir.json"), None),
                ("dir.json/in.json".to_string(), s("in.json"), None),
                ("real/r.json".to_string(), s("r.json"), None),
                ("top.json".to_string(), s("top.json"), None),
            ]
        );
        assert!(shard_pairs(&root.join("c"), &root.join("c")).is_empty());
    }

    #[test]
    fn variation_files_are_listed_as_glob_lists_them() {
        let scratch = Scratch::new("variations");
        let root = &scratch.0;
        for name in [
            "a/variations/1/all_vars.gz",
            "a/variations/1/all_vars.gz.csi",
            "a/variations/.3/all_vars.gz",
            "a/variations/2/other.gz",
            "a/variations/all_vars.gz",
            "b/variations/10/all_vars.gz",
            "b/variations/1/all_vars.gz",
        ] {
            scratch.file(name);
        }
        symlink(root.join("a/variations/1"), root.join("a/variations/4")).unwrap();
        fs::create_dir(root.join("a/variations/5")).unwrap();
        symlink("nowhere", root.join("a/variations/5/all_vars.gz")).unwrap();
        fs::create_dir_all(root.join("a/variations/6/all_vars.gz")).unwrap();
        symlink("nowhere", root.join("a/variations/7")).unwrap();
        let s = |text: &str| Some(text.to_string());
        let both = chromosome_pairs(&root.join("a"), &root.join("b"));
        assert_eq!(
            names(&both),
            [
                (".3".to_string(), s("all_vars.gz"), None),
                ("1".to_string(), s("all_vars.gz"), s("all_vars.gz")),
                ("10".to_string(), None, s("all_vars.gz")),
                ("4".to_string(), s("all_vars.gz"), None),
                ("5".to_string(), s("all_vars.gz"), None),
                ("6".to_string(), s("all_vars.gz"), None),
            ]
        );
        assert!(
            chromosome_pairs(&root.join("b/variations"), &root.join("b/variations")).is_empty()
        );
    }
}
