// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The workers of one reading of the two outputs: a queue of jobs any worker takes, the
//! count of blocks a side may have in work at once, and the results put back in order.
//!
//! A panic never leaves a thread of a reading waiting: a job's panic travels to its
//! collector in the place of its result, a queue closes when its reading ends in any way,
//! and a reader passes freely once its collector has ended.

use std::collections::{BTreeMap, VecDeque};
use std::panic::{catch_unwind, resume_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, Sender};
use std::sync::{Condvar, Mutex};
use std::thread::ScopedJoinHandle;

/// What a job ended with: its result, or the panic it ended in.
pub(crate) type Outcome<T> = std::thread::Result<T>;

pub(crate) type Job<'a> = Box<dyn FnOnce() + Send + 'a>;

struct Queue<'a> {
    jobs: VecDeque<Job<'a>>,
    /// Whether a worker returns once no job is waiting.
    closed: bool,
}

/// The jobs waiting for a worker.
pub(crate) struct Jobs<'a> {
    waiting: Mutex<Queue<'a>>,
    changed: Condvar,
}

impl<'a> Jobs<'a> {
    pub(crate) fn new() -> Jobs<'a> {
        Jobs {
            waiting: Mutex::new(Queue {
                jobs: VecDeque::new(),
                closed: false,
            }),
            changed: Condvar::new(),
        }
    }

    /// Queues `job`; once the queue is closed, when no worker is left to take it, runs it.
    pub(crate) fn push(&self, job: Job<'a>) {
        let mut waiting = self.waiting.lock().expect("job queue");
        if waiting.closed {
            drop(waiting);
            job();
            return;
        }
        waiting.jobs.push_back(job);
        drop(waiting);
        self.changed.notify_one();
    }

    /// Queues `work` as the job `index` of a reading; what it ends with goes to `results`.
    pub(crate) fn queue<T: Send + 'a>(
        &self,
        index: u64,
        results: Sender<(u64, Outcome<T>)>,
        work: impl FnOnce() -> T + Send + 'a,
    ) {
        self.push(Box::new(move || {
            let outcome = catch_unwind(AssertUnwindSafe(work));
            let _ = results.send((index, outcome));
        }));
    }

    /// After this the workers return once no job is waiting.
    pub(crate) fn close(&self) {
        self.waiting.lock().expect("job queue").closed = true;
        self.changed.notify_all();
    }

    /// Runs jobs until the queue is closed and empty.
    pub(crate) fn work(&self) {
        loop {
            let job = {
                let mut waiting = self.waiting.lock().expect("job queue");
                loop {
                    if let Some(job) = waiting.jobs.pop_front() {
                        break job;
                    }
                    if waiting.closed {
                        return;
                    }
                    waiting = self.changed.wait(waiting).expect("job queue");
                }
            };
            job();
        }
    }
}

/// Closes a queue when dropped, so that the workers of a reading return however the
/// reading ends.
pub(crate) struct Closing<'q, 'a>(pub(crate) &'q Jobs<'a>);

impl Drop for Closing<'_, '_> {
    fn drop(&mut self) {
        self.0.close();
    }
}

struct Free {
    count: usize,
    /// Whether every taker passes.
    open: bool,
}

/// The blocks of one side that may be in work at once: read, scanned or waiting for their
/// turn.
pub(crate) struct Permits {
    free: Mutex<Free>,
    changed: Condvar,
}

impl Permits {
    pub(crate) fn new(count: usize) -> Permits {
        Permits {
            free: Mutex::new(Free { count, open: false }),
            changed: Condvar::new(),
        }
    }

    pub(crate) fn take(&self) {
        let mut free = self.free.lock().expect("block permits");
        while free.count == 0 && !free.open {
            free = self.changed.wait(free).expect("block permits");
        }
        if !free.open {
            free.count -= 1;
        }
    }

    pub(crate) fn give(&self) {
        self.free.lock().expect("block permits").count += 1;
        self.changed.notify_one();
    }

    /// After this every taker passes, whatever the count.
    fn open(&self) {
        self.free.lock().expect("block permits").open = true;
        self.changed.notify_all();
    }
}

/// A side's collector at work. When it ends, nothing gives a permit back, so its reader
/// passes freely; when it ends in a panic, `stop` ends the readers.
pub(crate) struct Collecting<'a> {
    pub(crate) permits: &'a Permits,
    pub(crate) stop: &'a AtomicBool,
}

impl Drop for Collecting<'_> {
    fn drop(&mut self) {
        self.permits.open();
        if std::thread::panicking() {
            self.stop.store(true, Ordering::Relaxed);
        }
    }
}

/// What a scoped thread returned; a panic it ended in goes on in this thread.
pub(crate) fn joined<T>(thread: ScopedJoinHandle<'_, T>) -> T {
    thread
        .join()
        .unwrap_or_else(|payload| resume_unwind(payload))
}

/// Hands the numbered results of `results` to `absorb` in the order of their numbers,
/// which start at 0 and have no gap.
pub(crate) fn in_order<T>(results: Receiver<(u64, T)>, mut absorb: impl FnMut(T)) {
    let mut early: BTreeMap<u64, T> = BTreeMap::new();
    let mut next = 0;
    for (index, result) in results {
        early.insert(index, result);
        while let Some(result) = early.remove(&next) {
            absorb(result);
            next += 1;
        }
    }
    assert!(early.is_empty(), "a job ended without its result");
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;
    use std::sync::mpsc::channel;

    #[test]
    fn every_job_runs_once_and_results_come_back_in_order() {
        let ran = AtomicUsize::new(0);
        let jobs = Jobs::new();
        let (results, collected) = channel();
        let mut order = Vec::new();
        std::thread::scope(|scope| {
            for _ in 0..4 {
                scope.spawn(|| jobs.work());
            }
            for index in 0..500u64 {
                let results = results.clone();
                let ran = &ran;
                jobs.push(Box::new(move || {
                    ran.fetch_add(1, Ordering::Relaxed);
                    let _ = results.send((index, index * 3));
                }));
            }
            drop(results);
            in_order(collected, |value| order.push(value));
            jobs.close();
        });
        assert_eq!(ran.load(Ordering::Relaxed), 500);
        assert!(order.iter().copied().eq((0..500).map(|index| index * 3)));
    }

    #[test]
    fn a_job_s_panic_reaches_the_collector_and_a_closed_queue_frees_its_workers() {
        let jobs = Jobs::new();
        let (results, collected) = channel();
        let mut outcomes = Vec::new();
        std::thread::scope(|scope| {
            let _closing = Closing(&jobs);
            for _ in 0..2 {
                scope.spawn(|| jobs.work());
            }
            jobs.queue(0, results.clone(), || 7);
            jobs.queue(1, results.clone(), || panic!("a job's own panic"));
            jobs.queue(2, results, || 9);
            in_order(collected, |outcome: Outcome<i32>| {
                outcomes.push(outcome.ok())
            });
        });
        assert_eq!(outcomes, [Some(7), None, Some(9)]);
        // A job queued after the workers have gone runs where it is queued.
        let (results, collected) = channel();
        jobs.queue(0, results, || 11);
        assert_eq!(collected.recv().unwrap().1.ok(), Some(11));
    }

    #[test]
    fn a_collector_that_ends_lets_its_reader_pass() {
        let permits = Permits::new(1);
        let stop = AtomicBool::new(false);
        permits.take();
        drop(Collecting {
            permits: &permits,
            stop: &stop,
        });
        permits.take();
        permits.take();
        assert!(!stop.load(Ordering::Relaxed));
    }

    #[test]
    fn a_collector_that_panics_stops_the_readers_and_lets_its_own_pass() {
        let permits = Permits::new(0);
        let stop = AtomicBool::new(false);
        std::thread::scope(|scope| {
            let collector = scope.spawn(|| {
                let _collecting = Collecting {
                    permits: &permits,
                    stop: &stop,
                };
                panic!("a collector's own panic");
            });
            assert!(collector.join().is_err());
        });
        assert!(stop.load(Ordering::Relaxed));
        permits.take();
    }

    #[test]
    fn a_permit_taken_is_free_again_once_given() {
        let permits = Permits::new(2);
        permits.take();
        permits.take();
        std::thread::scope(|scope| {
            let waiting = scope.spawn(|| permits.take());
            permits.give();
            waiting.join().unwrap();
        });
        assert_eq!(permits.free.lock().unwrap().count, 0);
    }
}
