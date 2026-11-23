// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The repository the tool belongs to: where the reproduction records and the index page
//! are found by default, what a record's relative path is relative to, and what `{repo}`
//! and `REPO` name.

use crate::pypath;

const RECORDS: &str = "scripts/divergences/reproductions";

/// The nearest directory above the executable that holds the reproduction records, else
/// the checkout the tool was built in, with every symbolic link resolved.
pub fn repo_root() -> String {
    let executable = std::env::current_exe()
        .ok()
        .and_then(|path| path.canonicalize().ok());
    if let Some(executable) = executable {
        for ancestor in executable.ancestors().skip(1) {
            if ancestor.join(RECORDS).is_dir() {
                return ancestor.to_string_lossy().into_owned();
            }
        }
    }
    pypath::resolve(&pypath::parent(&pypath::parent(env!("CARGO_MANIFEST_DIR"))))
}
