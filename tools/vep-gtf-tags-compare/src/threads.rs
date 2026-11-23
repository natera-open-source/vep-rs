// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The worker-thread bound and the indexed fan-out the inputs are read on.

use std::num::NonZero;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// The environment variable that bounds the worker threads; unset or empty, the host's
/// available cores are the bound. An environment variable and not an option, so that the
/// command line holds the four options alone.
pub const THREADS_ENV: &str = "VEP_GTF_TAGS_COMPARE_THREADS";

/// The worker-thread bound of this run.
pub fn thread_count() -> Result<usize, String> {
    let available = || std::thread::available_parallelism().map_or(1, NonZero::get);
    match std::env::var_os(THREADS_ENV) {
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

/// Runs `task` on every index below `count`, on at most `threads` threads, and returns the
/// results in index order; after a failure no further index is started and the failure of
/// the lowest index is returned. One thread runs every index on the calling thread.
pub fn run_indexed<T: Send, E: Send>(
    threads: usize,
    count: usize,
    task: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let work = || {
        let mut done = Vec::new();
        while !failed.load(Ordering::Relaxed) {
            let index = next.fetch_add(1, Ordering::Relaxed);
            if index >= count {
                break;
            }
            let result = task(index);
            if result.is_err() {
                failed.store(true, Ordering::Relaxed);
            }
            done.push((index, result));
        }
        done
    };
    let workers = threads.min(count).max(1);
    let mut done: Vec<(usize, Result<T, E>)> = if workers == 1 {
        work()
    } else {
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers).map(|_| scope.spawn(work)).collect();
            handles
                .into_iter()
                .flat_map(|handle| {
                    handle
                        .join()
                        .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
                })
                .collect()
        })
    };
    done.sort_by_key(|(index, _)| *index);
    done.into_iter().map(|(_, result)| result).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn results_come_back_in_index_order_for_every_thread_count() {
        for threads in [1, 2, 7, 64] {
            let squares = run_indexed(threads, 100, |i| Ok::<_, String>(i * i)).unwrap();
            assert_eq!(squares, (0..100).map(|i| i * i).collect::<Vec<_>>());
        }
        assert!(run_indexed(4, 0, Ok::<usize, String>).unwrap().is_empty());
    }

    #[test]
    fn a_failure_is_returned_and_stops_the_run() {
        for threads in [1, 3] {
            let error = run_indexed(threads, 50, |i| {
                if i == 20 {
                    Err("index 20".to_string())
                } else {
                    Ok(i)
                }
            })
            .unwrap_err();
            assert_eq!(error, "index 20");
        }
    }
}
