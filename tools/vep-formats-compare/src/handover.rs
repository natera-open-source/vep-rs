// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs a comparison and ends a run that met an input no reading of the comparison expects.
//!
//! A comparison writes nothing, to a file or to a standard stream, before it has its whole
//! report, so a run that ends this way has left no output of its own.

use crate::argparse::standard_error;
use crate::error::Uncaught;
use crate::COMPONENT;
use std::ffi::OsString;
use std::panic;

/// The exit code of `main` over `argv`: its own, or 1 after one `ERROR:` line on the
/// standard error where it ended on an `Uncaught` payload. Any other panic goes on.
pub fn run(argv: &[OsString], main: impl FnOnce(&[OsString]) -> i32) -> i32 {
    let default_hook = panic::take_hook();
    panic::set_hook(Box::new(move |info| {
        if info.payload().downcast_ref::<Uncaught>().is_none() {
            default_hook(info);
        }
    }));
    match panic::catch_unwind(panic::AssertUnwindSafe(|| main(argv))) {
        Ok(code) => code,
        Err(payload) => match payload.downcast_ref::<Uncaught>() {
            Some(Uncaught(what)) => {
                eprintln!(
                    "{}",
                    standard_error(&format!("ERROR: [{COMPONENT}] {what}"))
                );
                1
            }
            None => panic::resume_unwind(payload),
        },
    }
}
