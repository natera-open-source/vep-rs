// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The bound on the threads of a run.

use std::ffi::OsStr;
use std::num::NonZero;

/// The environment variable that bounds the threads at work; unset or empty, the host's
/// available cores are the bound. An environment variable and not an option, since the
/// command line names the table and the bound changes no byte written.
pub const THREADS_ENV: &str = "VEP_DISCORDANT_CLASSES_THREADS";

/// The bound of this run.
pub fn thread_count() -> Result<usize, String> {
    bound(std::env::var_os(THREADS_ENV).as_deref())
}

fn bound(value: Option<&OsStr>) -> Result<usize, String> {
    let available = || std::thread::available_parallelism().map_or(1, NonZero::get);
    match value {
        None => Ok(available()),
        Some(value) if value.is_empty() => Ok(available()),
        Some(value) => value
            .to_str()
            .and_then(|text| text.trim().parse::<usize>().ok())
            .filter(|threads| *threads >= 1)
            .ok_or_else(|| {
                format!(
                    "{THREADS_ENV}={}: expected a whole number of threads, 1 or more",
                    value.to_string_lossy()
                )
            }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_unset_or_empty_bound_is_the_host_s_cores() {
        let cores = std::thread::available_parallelism().map_or(1, NonZero::get);
        assert_eq!(bound(None), Ok(cores));
        assert_eq!(bound(Some(OsStr::new(""))), Ok(cores));
    }

    #[test]
    fn a_whole_number_of_one_or_more_is_the_bound() {
        assert_eq!(bound(Some(OsStr::new("1"))), Ok(1));
        assert_eq!(bound(Some(OsStr::new(" 12 "))), Ok(12));
    }

    #[test]
    fn any_other_value_is_refused_by_name() {
        for value in ["0", "-1", "two", "1.5", " "] {
            assert_eq!(
                bound(Some(OsStr::new(value))),
                Err(format!(
                    "VEP_DISCORDANT_CLASSES_THREADS={value}: expected a whole number of threads, 1 or more"
                ))
            );
        }
    }
}
