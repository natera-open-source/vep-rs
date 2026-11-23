// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's parts as a library, so that the tests reach them: the command line as
//! `argparse` reads it (`argparse`), text as CPython reads and writes it (`py`), a table
//! as `csv` reads and writes it (`csv`) and the join (`join`).

pub mod argparse;
pub mod csv;
pub mod join;
pub mod py;
