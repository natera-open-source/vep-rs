// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The divergence re-evaluation tool: the reproduction records, the index page, the two
//! engines' runs and the comparison of what each record names.

pub mod check;
pub mod cli;
pub mod dump;
pub mod engine;
pub mod evaluate;
pub mod extract;
mod gzip;
pub mod item;
pub mod page;
pub mod pyfmt;
pub mod pyjson;
pub mod pypath;
pub mod pyre;
pub mod record;
pub mod repo;
pub mod schema;
pub mod tempdir;
#[cfg(test)]
mod testing;
pub mod value;
