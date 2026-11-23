// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's parts as a library, so that the tests reach them.
//!
//! The command line as `argparse` reads it and text as CPython reads and writes it are
//! the modules of `tools/vep-storable-digest-join`, compiled here from their own files.
//! Over them: a dump as the classification reads it (`dumps`), the lines two blocks
//! differ in (`difflib`), the stored translations of a JSON cache (`jsoncache`) and the
//! classification (`classify`).

#[path = "../../vep-storable-digest-join/src/argparse.rs"]
pub mod argparse;
#[path = "../../vep-storable-digest-join/src/py.rs"]
pub mod py;

pub mod classify;
pub mod difflib;
pub mod dumps;
pub mod jsoncache;
