// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Paths and failed system calls as CPython's `pathlib`, `os` and `tempfile` read and word
//! them: the text of a path argument, its parent and name, the tests that never fail, the
//! directories made on the way to one, the scratch directory, and `str(OSError)`.

use super::text::repr;
use std::ffi::CString;
use std::fs;
use std::hash::{BuildHasher, Hasher};
use std::io;
use std::os::unix::fs::DirBuilderExt;
use std::path::{Path, PathBuf};

/// A failed system call as CPython's `OSError` holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OsFailure {
    pub errno: Option<i32>,
    /// `strerror` of the error number, or the error's own words where it has no number.
    pub reason: String,
    pub filename: Option<String>,
}

impl OsFailure {
    pub fn of(error: &io::Error, filename: Option<&str>) -> OsFailure {
        let told = error.to_string();
        let errno = error.raw_os_error();
        let reason = match errno {
            Some(code) => told
                .strip_suffix(&format!(" (os error {code})"))
                .unwrap_or(&told)
                .to_string(),
            None => told,
        };
        OsFailure {
            errno,
            reason,
            filename: filename.map(str::to_string),
        }
    }

    pub fn at(error: &io::Error, path: &Path) -> OsFailure {
        OsFailure::of(error, Some(&path.to_string_lossy()))
    }

    /// `str(error)`: `[Errno N] <reason>: '<file name>'`, the file name as `repr` writes it
    /// and left out where the call named none.
    pub fn text(&self) -> String {
        match (self.errno, &self.filename) {
            (Some(code), Some(name)) => format!("[Errno {code}] {}: {}", self.reason, repr(name)),
            (Some(code), None) => format!("[Errno {code}] {}", self.reason),
            (None, _) => self.reason.clone(),
        }
    }
}

/// `str(pathlib.Path(argument))`: repeated separators and `.` components dropped, a trailing
/// separator dropped, exactly two leading separators kept, `..` left alone, `.` for an empty
/// result.
pub fn path_text(argument: &str) -> String {
    let leading = argument.bytes().take_while(|&byte| byte == b'/').count();
    let root = match leading {
        0 => "",
        2 => "//",
        _ => "/",
    };
    let parts: Vec<&str> = argument
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    let joined = format!("{root}{}", parts.join("/"));
    if joined.is_empty() {
        ".".to_string()
    } else {
        joined
    }
}

/// Where the last component of a path written by `path_text` starts.
fn tail(path: &str) -> usize {
    path.rfind('/').map_or(0, |at| at + 1)
}

/// `str(Path(path).parent)` of a path written by `path_text`.
pub fn parent(path: &str) -> String {
    let at = tail(path);
    if at == 0 {
        return ".".to_string();
    }
    let kept = path[..at].trim_end_matches('/');
    if kept.is_empty() {
        // The path is its root and a name; two leading separators are a root of their own.
        return if path.starts_with("//") && !path.starts_with("///") {
            "//".to_string()
        } else {
            "/".to_string()
        };
    }
    kept.to_string()
}

/// `Path(path).name` of a path written by `path_text`: empty for `.` and for a root.
pub fn name(path: &str) -> &str {
    if path == "." {
        return "";
    }
    &path[tail(path)..]
}

pub fn is_file(path: &str) -> bool {
    fs::metadata(path).is_ok_and(|found| found.is_file())
}

pub fn is_dir(path: &str) -> bool {
    fs::metadata(path).is_ok_and(|found| found.is_dir())
}

pub fn exists(path: &str) -> bool {
    fs::metadata(path).is_ok()
}

/// `os.path.lexists`: whether the path names anything, a link that leads nowhere included.
pub fn lexists(path: &str) -> bool {
    fs::symlink_metadata(path).is_ok()
}

pub fn is_symlink(path: &str) -> bool {
    fs::symlink_metadata(path).is_ok_and(|found| found.file_type().is_symlink())
}

/// `os.access(path, os.W_OK)`.
pub fn writable(path: &str) -> bool {
    let Ok(name) = CString::new(path.as_bytes()) else {
        return false;
    };
    // SAFETY: `name` is a NUL-terminated string that outlives the call.
    unsafe { libc::access(name.as_ptr(), libc::W_OK) == 0 }
}

/// `Path(path).mkdir(parents=True, exist_ok=True)`: the directory and every missing one on
/// the way to it; the failure names the directory whose creation the system refused.
pub fn make_directories(path: &str) -> Result<(), OsFailure> {
    fn make(path: &str, parents: bool) -> Result<(), OsFailure> {
        match fs::create_dir(path) {
            Ok(()) => Ok(()),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                let above = parent(path);
                if !parents || above == path {
                    return Err(OsFailure::of(&error, Some(path)));
                }
                make(&above, true)?;
                make(path, false)
            }
            Err(_) if is_dir(path) => Ok(()),
            Err(error) => Err(OsFailure::of(&error, Some(path))),
        }
    }
    make(path, true)
}

/// The characters of a name `tempfile` makes.
const NAME_CHARACTERS: &[u8] = b"abcdefghijklmnopqrstuvwxyz0123456789_";
/// The names `tempfile` tries before it gives up.
const MOST_NAMES: usize = 10_000;

fn random_name() -> String {
    let mut bits = std::collections::hash_map::RandomState::new()
        .build_hasher()
        .finish();
    let count = NAME_CHARACTERS.len() as u64;
    (0..8)
        .map(|_| {
            let character = NAME_CHARACTERS[(bits % count) as usize];
            bits /= count;
            character as char
        })
        .collect()
}

/// `tempfile.gettempdir()`: the first of `TMPDIR`, `TEMP`, `TMP`, `/tmp`, `/var/tmp`,
/// `/usr/tmp` and the working directory in which a file can be made and written.
pub fn temporary_directory() -> Result<String, OsFailure> {
    let mut candidates: Vec<String> = ["TMPDIR", "TEMP", "TMP"]
        .iter()
        .filter_map(|name| std::env::var(name).ok())
        .filter(|named| !named.is_empty())
        .collect();
    candidates.extend(["/tmp", "/var/tmp", "/usr/tmp"].map(str::to_string));
    candidates.push(
        std::env::current_dir().map_or_else(|_| ".".to_string(), |at| at.display().to_string()),
    );
    for candidate in &candidates {
        let directory = match std::path::absolute(candidate) {
            Ok(absolute) => PathBuf::from(path_text(&absolute.to_string_lossy())),
            Err(_) => PathBuf::from(candidate),
        };
        for _ in 0..100 {
            let trial = directory.join(random_name());
            let made = fs::OpenOptions::new()
                .write(true)
                .read(true)
                .create_new(true)
                .open(&trial);
            match made {
                Ok(mut file) => {
                    let written = io::Write::write_all(&mut file, b"blat");
                    drop(file);
                    let _ = fs::remove_file(&trial);
                    if written.is_ok() {
                        return Ok(directory.to_string_lossy().into_owned());
                    }
                    break;
                }
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(_) => break,
            }
        }
    }
    let listed: Vec<String> = candidates.iter().map(|name| repr(name)).collect();
    Err(OsFailure {
        errno: Some(libc::ENOENT),
        reason: format!(
            "No usable temporary directory found in [{}]",
            listed.join(", ")
        ),
        filename: None,
    })
}

/// `tempfile.mkdtemp(prefix=prefix, dir=directory)`: a new directory with mode 0700
/// under `directory`, the system's temporary directory where none is given, named `prefix`
/// and eight random characters.
pub fn make_scratch(directory: Option<&str>, prefix: &str) -> Result<PathBuf, OsFailure> {
    let directory = match directory {
        Some(given) => given.to_string(),
        None => temporary_directory()?,
    };
    for _ in 0..MOST_NAMES {
        let path = Path::new(&directory).join(format!("{prefix}{}", random_name()));
        match fs::DirBuilder::new().mode(0o700).create(&path) {
            Ok(()) => return Ok(path),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(OsFailure::at(&error, &path)),
        }
    }
    Err(OsFailure {
        errno: Some(libc::EEXIST),
        reason: "No usable temporary directory name found".to_string(),
        filename: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_argument_is_written_as_pathlib_writes_it() {
        for (argument, text, above, last) in [
            ("a//b", "a/b", "a", "b"),
            ("./a", "a", ".", "a"),
            ("a/", "a", ".", "a"),
            ("", ".", ".", ""),
            (".", ".", ".", ""),
            ("//a", "//a", "//", "a"),
            ("///a", "/a", "/", "a"),
            ("a/./b", "a/b", "a", "b"),
            ("a/../b", "a/../b", "a/..", "b"),
            ("/", "/", "/", ""),
            ("//", "//", "//", ""),
            ("a/b/.", "a/b", "a", "b"),
            ("./", ".", ".", ""),
            ("././a//", "a", ".", "a"),
            ("/a/b", "/a/b", "/a", "b"),
            ("/a", "/a", "/", "a"),
        ] {
            let written = path_text(argument);
            assert_eq!(written, text, "{argument}");
            assert_eq!(parent(&written), above, "{argument}");
            assert_eq!(name(&written), last, "{argument}");
        }
    }

    #[test]
    fn failed_call_reads_as_the_interpreter_words_it() {
        let missing = io::Error::from_raw_os_error(2);
        assert_eq!(
            OsFailure::of(&missing, Some("/no/such")).text(),
            "[Errno 2] No such file or directory: '/no/such'"
        );
        assert_eq!(
            OsFailure::of(&io::Error::from_raw_os_error(13), Some("it's")).text(),
            "[Errno 13] Permission denied: \"it's\""
        );
        assert_eq!(
            OsFailure::of(&io::Error::from_raw_os_error(28), None).text(),
            "[Errno 28] No space left on device"
        );
        assert_eq!(
            OsFailure::of(&io::Error::other("no number"), Some("x")).text(),
            "no number"
        );
    }

    #[test]
    fn directories_are_made_on_the_way_and_the_refused_one_is_named() {
        let root =
            std::env::temp_dir().join(format!("vep-plugin-pairs-paths-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        let deep = root.join("a/b/c");
        let deep = deep.to_str().unwrap();
        make_directories(deep).unwrap();
        make_directories(deep).unwrap();
        assert!(is_dir(deep));
        let file = root.join("a/file");
        fs::write(&file, b"").unwrap();
        let file = file.to_str().unwrap();
        assert_eq!(
            make_directories(file).unwrap_err().text(),
            format!("[Errno 17] File exists: '{file}'")
        );
        let below = format!("{file}/x/y");
        assert_eq!(
            make_directories(&below).unwrap_err().text(),
            format!("[Errno 20] Not a directory: '{below}'")
        );
        let scratch = make_scratch(Some(deep), "score_plugin_tuples_").unwrap();
        let made = scratch.file_name().unwrap().to_str().unwrap().to_string();
        assert_eq!(made.len(), "score_plugin_tuples_".len() + 8);
        assert!(writable(deep) && !writable(&below));
        let nowhere = format!("{deep}/missing");
        let refused = make_scratch(Some(&nowhere), "score_plugin_tuples_").unwrap_err();
        assert!(refused.text().starts_with(&format!(
            "[Errno 2] No such file or directory: '{nowhere}/score_plugin_tuples_"
        )));
        fs::remove_dir_all(&root).unwrap();
    }
}
