// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The directories of a run as `pathlib` names them, the inputs a directory holds and the
//! contig naming of a cache directory.

use crate::lines::has_chr_prefix;
use crate::pyfmt::{py_path_str, py_str_cmp};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs;
use std::io;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::Path;

pub fn as_path(bytes: &[u8]) -> &Path {
    Path::new(OsStr::from_bytes(bytes))
}

/// `dir / name` for a directory `pathlib` has resolved.
pub fn join(dir: &[u8], name: &[u8]) -> Vec<u8> {
    if dir == b"/" {
        [b"/", name].concat()
    } else {
        [dir, b"/", name].concat()
    }
}

enum Part {
    Name(Vec<u8>),
    /// Below it on the stack stood the parts of a symbolic link's target; the path
    /// reached when it is popped is what the link resolves to.
    Resolved(Vec<u8>),
}

fn reversed_names(path: &[u8]) -> impl Iterator<Item = Part> + '_ {
    path.rsplit(|byte| *byte == b'/')
        .map(|name| Part::Name(name.to_vec()))
}

/// `str(pathlib.Path(arg).resolve())`: the path made absolute against the working
/// directory, `.` and `..` components taken out, and every symbolic link on the way
/// followed. A component that does not exist, or cannot be read, is kept as it is named,
/// and so is a link that leads back to itself.
pub fn resolve(arg: &[u8]) -> io::Result<Vec<u8>> {
    let filename = py_path_str(arg);
    let mut path: Vec<u8> = if filename.starts_with(b"/") {
        b"/".to_vec()
    } else {
        std::env::current_dir()?.into_os_string().into_vec()
    };
    let mut rest: Vec<Part> = reversed_names(&filename).collect();
    // A link's path, and what it resolves to once that is known.
    let mut seen: HashMap<Vec<u8>, Option<Vec<u8>>> = HashMap::new();
    while let Some(part) = rest.pop() {
        let name = match part {
            Part::Name(name) => name,
            Part::Resolved(link) => {
                seen.insert(link, Some(path.clone()));
                continue;
            }
        };
        if name.is_empty() || name == b"." {
            continue;
        }
        if name == b".." {
            let parent = path.iter().rposition(|byte| *byte == b'/').unwrap_or(0);
            path.truncate(parent.max(1));
            continue;
        }
        let newpath = join(&path, &name);
        let is_link = fs::symlink_metadata(as_path(&newpath))
            .map(|meta| meta.file_type().is_symlink())
            .unwrap_or(false);
        if !is_link {
            path = newpath;
            continue;
        }
        match seen.get(&newpath) {
            Some(Some(resolved)) => {
                path = resolved.clone();
                continue;
            }
            Some(None) => {
                // The link is on its own way: a loop.
                path = newpath;
                continue;
            }
            None => {}
        }
        let Ok(target) = fs::read_link(as_path(&newpath)) else {
            path = newpath;
            continue;
        };
        let target = target.into_os_string().into_vec();
        if target.starts_with(b"/") {
            path = b"/".to_vec();
        }
        seen.insert(newpath.clone(), None);
        rest.push(Part::Resolved(newpath));
        rest.extend(reversed_names(&target));
    }
    Ok(path)
}

/// The names of the inputs in `dir`: every entry whose name ends in `.vcf`, then every
/// entry whose name ends in `.vcf.gz`, each group in code point order. A directory that
/// cannot be listed holds none.
pub fn list_inputs(dir: &[u8]) -> Vec<Vec<u8>> {
    let names: Vec<Vec<u8>> = fs::read_dir(as_path(dir))
        .and_then(|entries| {
            entries
                .map(|entry| entry.map(|e| e.file_name().into_vec()))
                .collect()
        })
        .unwrap_or_default();
    let mut inputs = Vec::new();
    for suffix in [&b".vcf"[..], b".vcf.gz"] {
        let mut group: Vec<Vec<u8>> = names
            .iter()
            .filter(|name| name.ends_with(suffix))
            .cloned()
            .collect();
        group.sort_by(|a, b| py_str_cmp(a, b));
        inputs.extend(group);
    }
    inputs
}

/// Whether an input is read through gzip: by its name alone.
pub fn is_gzip_name(name: &[u8]) -> bool {
    name.ends_with(b".gz")
}

/// The name an input is written under: its own, less a `.gz`.
pub fn prepared_name(name: &[u8]) -> &[u8] {
    name.strip_suffix(b".gz").unwrap_or(name)
}

/// The contig naming of a first variant line or of a cache directory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Style {
    Chr,
    Plain,
    Unknown,
}

impl Style {
    pub fn of(contig: &[u8]) -> Style {
        if has_chr_prefix(contig) {
            Style::Chr
        } else {
            Style::Plain
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Style::Chr => "chr",
            Style::Plain => "plain",
            Style::Unknown => "unknown",
        }
    }
}

/// The naming of a cache: that of the first directory, in code point order, under its
/// `transcripts` directory, else under its `variations` directory; unknown where neither
/// holds a directory.
pub fn cache_style(cache_dir: Option<&[u8]>) -> io::Result<Style> {
    let Some(cache_dir) = cache_dir else {
        return Ok(Style::Unknown);
    };
    for sub in [&b"transcripts"[..], b"variations"] {
        let root = join(cache_dir, sub);
        let root = as_path(&root);
        if !root.is_dir() {
            continue;
        }
        let mut first: Option<Vec<u8>> = None;
        for entry in fs::read_dir(root)? {
            let entry = entry?;
            if !entry.path().is_dir() {
                continue;
            }
            let name = entry.file_name().into_vec();
            if first
                .as_ref()
                .is_none_or(|least| py_str_cmp(&name, least).is_lt())
            {
                first = Some(name);
            }
        }
        if let Some(name) = first {
            return Ok(Style::of(&name));
        }
    }
    Ok(Style::Unknown)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::symlink;
    use std::path::PathBuf;

    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let path = std::env::temp_dir().canonicalize().unwrap().join(format!(
                "vep-prepare-benchmark-vcfs-paths-{}-{tag}",
                std::process::id()
            ));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(&path).unwrap();
            Scratch(path)
        }

        fn at(&self, rel: &str) -> Vec<u8> {
            self.0.join(rel).into_os_string().into_vec()
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn text(bytes: Vec<u8>) -> String {
        String::from_utf8(bytes).unwrap()
    }

    #[test]
    fn a_path_resolves_through_its_links_and_keeps_what_does_not_exist() {
        let scratch = Scratch::new("resolve");
        let root = scratch.0.to_str().unwrap().to_string();
        fs::create_dir_all(scratch.0.join("real/inner")).unwrap();
        symlink("real", scratch.0.join("link")).unwrap();
        symlink(scratch.0.join("real/inner"), scratch.0.join("absolute")).unwrap();
        symlink("../real/./inner/..", scratch.0.join("real/up")).unwrap();
        symlink("loop_b", scratch.0.join("loop_a")).unwrap();
        symlink("loop_a", scratch.0.join("loop_b")).unwrap();
        symlink("missing/target", scratch.0.join("dangling")).unwrap();
        let cases = [
            ("real", "real"),
            ("real/", "real"),
            ("./real//inner/.", "real/inner"),
            ("link", "real"),
            ("link/inner", "real/inner"),
            ("link/inner/..", "real"),
            ("absolute", "real/inner"),
            ("absolute/../x", "real/x"),
            ("real/up", "real"),
            ("link/up/inner", "real/inner"),
            ("absent", "absent"),
            ("absent/../real", "real"),
            ("absent/deeper/still", "absent/deeper/still"),
            ("real/inner/../../link/absent", "real/absent"),
            ("loop_a", "loop_a"),
            ("loop_a/x", "loop_a/x"),
            ("dangling", "missing/target"),
            ("dangling/x", "missing/target/x"),
        ];
        for (arg, expected) in cases {
            assert_eq!(
                text(resolve(&scratch.at(arg)).unwrap()),
                format!("{root}/{expected}"),
                "{arg}"
            );
        }
        assert_eq!(text(resolve(b"/").unwrap()), "/");
        assert_eq!(text(resolve(b"//").unwrap()), "/");
        assert_eq!(text(resolve(b"/..").unwrap()), "/");
        assert_eq!(text(resolve(b"/../../x/..").unwrap()), "/");
        let cwd = std::env::current_dir().unwrap();
        let cwd = cwd.to_str().unwrap();
        assert_eq!(text(resolve(b"").unwrap()), cwd);
        assert_eq!(text(resolve(b".").unwrap()), cwd);
        assert_eq!(
            text(resolve(b"absent-name/./x").unwrap()),
            format!("{cwd}/absent-name/x")
        );
    }

    #[test]
    fn inputs_are_the_vcf_names_then_the_gzipped_ones_each_in_code_point_order() {
        let scratch = Scratch::new("list");
        for name in [
            "b.vcf",
            "a.vcf.gz",
            "B.vcf",
            "a.vcf",
            ".vcf",
            "c.vcf.gz",
            "notes.txt",
            "x.VCF",
            "y.vcf.GZ",
            "z.vcf.gz.tbi",
            "vcf",
            "d.bcf",
        ] {
            fs::write(scratch.0.join(name), b"").unwrap();
        }
        fs::create_dir(scratch.0.join("dir.vcf")).unwrap();
        symlink("nowhere", scratch.0.join("broken.vcf.gz")).unwrap();
        fs::write(scratch.0.join(OsStr::from_bytes(b"\xff.vcf")), b"").unwrap();
        fs::write(scratch.0.join("\u{e000}.vcf"), b"").unwrap();
        let listed = list_inputs(scratch.0.as_os_str().as_bytes());
        let expected: Vec<Vec<u8>> = [
            &b".vcf"[..],
            b"B.vcf",
            b"a.vcf",
            b"b.vcf",
            b"dir.vcf",
            b"\xff.vcf",
            "\u{e000}.vcf".as_bytes(),
            b"a.vcf.gz",
            b"broken.vcf.gz",
            b"c.vcf.gz",
        ]
        .iter()
        .map(|name| name.to_vec())
        .collect();
        assert_eq!(listed, expected);
        assert!(list_inputs(&scratch.at("absent")).is_empty());
        assert!(list_inputs(&scratch.at("a.vcf")).is_empty());
    }

    #[test]
    fn an_input_is_written_under_its_name_less_the_gzip_suffix() {
        assert_eq!(prepared_name(b"a.vcf.gz"), b"a.vcf");
        assert_eq!(prepared_name(b"a.vcf"), b"a.vcf");
        assert_eq!(prepared_name(b".vcf.gz"), b".vcf");
        assert!(is_gzip_name(b"a.vcf.gz"));
        assert!(!is_gzip_name(b"a.vcf"));
        assert_eq!(join(b"/", b"a.vcf"), b"/a.vcf");
        assert_eq!(join(b"/in", b"a.vcf"), b"/in/a.vcf");
    }

    #[test]
    fn a_cache_is_named_as_its_first_contig_directory_is() {
        let scratch = Scratch::new("cache");
        let style = |rel: &str| cache_style(Some(&scratch.at(rel))).unwrap();
        assert_eq!(cache_style(None).unwrap(), Style::Unknown);
        assert_eq!(style("absent"), Style::Unknown);

        fs::create_dir_all(scratch.0.join("plain/transcripts/1")).unwrap();
        fs::create_dir_all(scratch.0.join("plain/transcripts/chr2")).unwrap();
        assert_eq!(style("plain"), Style::Plain);

        fs::create_dir_all(scratch.0.join("chr/transcripts/chr1")).unwrap();
        fs::create_dir_all(scratch.0.join("chr/transcripts/chrX")).unwrap();
        fs::write(scratch.0.join("chr/transcripts/1"), b"a file, not a contig").unwrap();
        assert_eq!(style("chr"), Style::Chr);

        // Capitals sort first, and the prefix is read in any letter case.
        fs::create_dir_all(scratch.0.join("upper/transcripts/CHR1")).unwrap();
        fs::create_dir_all(scratch.0.join("upper/transcripts/1")).unwrap();
        assert_eq!(style("upper"), Style::Plain);
        fs::remove_dir(scratch.0.join("upper/transcripts/1")).unwrap();
        assert_eq!(style("upper"), Style::Chr);

        // An empty `transcripts` hands on to `variations`; a link to a directory counts.
        fs::create_dir_all(scratch.0.join("second/transcripts")).unwrap();
        fs::create_dir_all(scratch.0.join("second/variations")).unwrap();
        fs::write(scratch.0.join("second/transcripts/chr1"), b"").unwrap();
        assert_eq!(style("second"), Style::Unknown);
        symlink(
            scratch.0.join("chr/transcripts/chr1"),
            scratch.0.join("second/variations/chrY"),
        )
        .unwrap();
        assert_eq!(style("second"), Style::Chr);

        fs::write(scratch.0.join("file"), b"").unwrap();
        assert_eq!(style("file"), Style::Unknown);
    }
}
