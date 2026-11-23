// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Paths as `pathlib.PurePosixPath` writes and joins them, `os.path.realpath` resolves them
//! and `OSError` names them. A path is kept as the text Python's `str(path)` gives, so
//! every message and every argument an engine is given carries that text.

use crate::pyfmt::py_repr_str;
use std::io;
use std::path::Path;

/// `str(pathlib.PurePosixPath(text))`: repeated separators and `.` components dropped, a
/// trailing separator dropped, exactly two leading slashes kept, `..` left alone, `.` for
/// an empty result.
pub fn normalise(text: &str) -> String {
    let leading = text.bytes().take_while(|b| *b == b'/').count();
    let root = match leading {
        0 => "",
        2 => "//",
        _ => "/",
    };
    let parts: Vec<&str> = text
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

fn is_absolute(text: &str) -> bool {
    text.starts_with('/')
}

/// `str(Path(base) / other)`: an absolute right side stands alone.
pub fn join(base: &str, other: &str) -> String {
    if is_absolute(other) {
        return normalise(other);
    }
    if other.is_empty() {
        return normalise(base);
    }
    if base == "." {
        return normalise(other);
    }
    if base.ends_with('/') {
        return normalise(&format!("{base}{other}"));
    }
    normalise(&format!("{base}/{other}"))
}

/// The root of a normalised path: `/`, `//` or nothing.
fn root_of(path: &str) -> &str {
    if path.starts_with("//") && !path.starts_with("///") {
        "//"
    } else if path.starts_with('/') {
        "/"
    } else {
        ""
    }
}

/// `Path(path).name`: the last component, empty for a root and for `.`.
pub fn name(path: &str) -> &str {
    let rest = &path[root_of(path).len()..];
    if rest == "." {
        return "";
    }
    rest.rsplit('/').next().unwrap_or("")
}

/// `str(Path(path).parent)`.
pub fn parent(path: &str) -> String {
    let root = root_of(path);
    let rest = &path[root.len()..];
    match rest.rfind('/') {
        Some(cut) => format!("{root}{}", &rest[..cut]),
        None if root.is_empty() => ".".to_string(),
        None => root.to_string(),
    }
}

/// `str(Path(path).with_name(new))`; the error is the `ValueError` text.
pub fn with_name(path: &str, new: &str) -> Result<String, String> {
    if new.is_empty() || new.contains('/') || new == "." {
        return Err(format!("Invalid name {}", py_repr_str(new)));
    }
    if name(path).is_empty() {
        return Err(format!(
            "PosixPath({}) has an empty name",
            py_repr_str(path)
        ));
    }
    let directory = parent(path);
    Ok(if directory == "." && !path.starts_with("./") {
        new.to_string()
    } else {
        join(&directory, new)
    })
}

/// `Path(path).suffix` of a file name: the last dot and what follows it, unless only dots
/// stand before it.
pub fn suffix(file_name: &str) -> &str {
    match file_name.rfind('.') {
        Some(dot) if !file_name[..dot].trim_matches('.').is_empty() => &file_name[dot..],
        _ => "",
    }
}

/// `Path(path).stem` of a file name.
pub fn stem(file_name: &str) -> &str {
    &file_name[..file_name.len() - suffix(file_name).len()]
}

/// `os.path.realpath(path)`: absolute, every symbolic link resolved, a component that does
/// not exist kept as written and `..` after it taken lexically.
pub fn realpath(path: &str) -> String {
    let mut resolved = if path.starts_with('/') {
        "/".to_string()
    } else {
        std::env::current_dir()
            .map(|dir| dir.to_string_lossy().into_owned())
            .unwrap_or_else(|_| "/".to_string())
    };
    // The unresolved components, last first; a link's own path follows a None marker so
    // that its resolved target is recorded when the marker is reached.
    let mut rest: Vec<Option<String>> =
        path.split('/').rev().map(|p| Some(p.to_string())).collect();
    let mut seen: Vec<(String, Option<String>)> = Vec::new();
    while let Some(part) = rest.pop() {
        let Some(part) = part else {
            let link = rest.pop().flatten().expect("a link path under its marker");
            if let Some(entry) = seen.iter_mut().find(|(path, _)| *path == link) {
                entry.1 = Some(resolved.clone());
            }
            continue;
        };
        if part.is_empty() || part == "." {
            continue;
        }
        if part == ".." {
            let cut = resolved.rfind('/').expect("an absolute path");
            resolved.truncate(cut.max(1));
            continue;
        }
        let candidate = if resolved == "/" {
            format!("/{part}")
        } else {
            format!("{resolved}/{part}")
        };
        let is_link = std::fs::symlink_metadata(&candidate)
            .map(|meta| meta.file_type().is_symlink())
            .unwrap_or(false);
        if !is_link {
            resolved = candidate;
            continue;
        }
        if let Some((_, target)) = seen.iter().find(|(path, _)| *path == candidate) {
            // A link seen before: its target when resolved, itself inside a loop.
            resolved = target.clone().unwrap_or(candidate);
            continue;
        }
        let Ok(target) = std::fs::read_link(&candidate) else {
            resolved = candidate;
            continue;
        };
        let target = target.to_string_lossy().into_owned();
        if target.starts_with('/') {
            resolved = "/".to_string();
        }
        seen.push((candidate.clone(), None));
        rest.push(Some(candidate));
        rest.push(None);
        rest.extend(target.split('/').rev().map(|p| Some(p.to_string())));
    }
    resolved
}

/// `str(Path(path).absolute())`: a relative path joined to the working directory, with no
/// link resolved and `..` left alone; an absolute one as written.
pub fn absolute(path: &str) -> String {
    if is_absolute(path) {
        return path.to_string();
    }
    let directory = std::env::current_dir()
        .map(|dir| dir.to_string_lossy().into_owned())
        .unwrap_or_else(|_| "/".to_string());
    join(&directory, path)
}

/// `Path(path).resolve()` as text.
pub fn resolve(path: &str) -> String {
    realpath(path)
}

pub fn is_dir(path: &str) -> bool {
    Path::new(path).is_dir()
}

pub fn is_file(path: &str) -> bool {
    Path::new(path).is_file()
}

pub fn exists(path: &str) -> bool {
    Path::new(path).exists()
}

/// `os.strerror(errno)` of an I/O error, or its own text when the system gave no number.
fn strerror(error: &io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(cut) => text[..cut].to_string(),
        None => text,
    }
}

/// `str(OSError)` for a failure on one path: `[Errno N] <strerror>: '<path>'`.
pub fn os_error(error: &io::Error, path: &str) -> String {
    match error.raw_os_error() {
        Some(code) => format!("[Errno {code}] {}: {}", strerror(error), py_repr_str(path)),
        None => strerror(error),
    }
}

#[cfg(test)]
#[path = "pypath_tests.rs"]
mod tests;
