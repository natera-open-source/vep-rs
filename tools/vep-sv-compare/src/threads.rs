// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The bound on worker threads and the fan-outs the comparison runs on.
//!
//! A fan-out returns its results in the order of its indices and reports the failure of
//! the lowest index that fails, so what it returns is the same at every thread count.
//! Fan-outs nested in one another (the files of a run, the masks of a file, the records a
//! mask walks) share the threads of one crew: see [`as_crew`].

use std::cell::RefCell;
use std::collections::{BTreeMap, VecDeque};
use std::convert::Infallible;
use std::num::NonZero;
use std::ops::ControlFlow;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, PoisonError};

use crate::failure::Failure;

/// The environment variable that bounds the worker threads; unset or empty, the host's
/// available cores are the bound. An environment variable and not an option, so that the
/// command line stays the nine options of the comparison.
pub const THREADS_ENV: &str = "VEP_SV_COMPARE_THREADS";

/// The worker-thread bound of this run.
pub fn thread_count() -> Result<usize, Failure> {
    let available = || std::thread::available_parallelism().map_or(1, NonZero::get);
    match std::env::var_os(THREADS_ENV) {
        None => Ok(available()),
        Some(value) if value.is_empty() => Ok(available()),
        Some(value) => threads_from(&value.to_string_lossy()),
    }
}

/// The bound a non-empty value of [`THREADS_ENV`] names: a whole number, 1 or more.
fn threads_from(value: &str) -> Result<usize, Failure> {
    value
        .trim()
        .parse::<usize>()
        .ok()
        .filter(|threads| *threads >= 1)
        .ok_or_else(|| {
            Failure::value(format!(
                "{THREADS_ENV}={value}: expected a whole number of threads, 1 or more"
            ))
        })
}

/// What a scoped thread returned, its panic carried on to the caller.
pub fn joined<T>(handle: std::thread::ScopedJoinHandle<'_, T>) -> T {
    handle
        .join()
        .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
}

/// The threads of a crew that no fan-out is using.
struct Crew {
    idle: Mutex<usize>,
    freed: Condvar,
}

impl Crew {
    fn idle(&self) -> std::sync::MutexGuard<'_, usize> {
        self.idle.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Takes one idle thread, if there is one.
    fn borrow(&self) -> bool {
        let mut idle = self.idle();
        let free = *idle > 0;
        if free {
            *idle -= 1;
        }
        free
    }

    fn give(&self) {
        *self.idle() += 1;
        self.freed.notify_one();
    }

    /// Takes one idle thread, waiting for one.
    fn take(&self) {
        let mut idle = self.idle();
        while *idle == 0 {
            idle = self
                .freed
                .wait(idle)
                .unwrap_or_else(PoisonError::into_inner);
        }
        *idle -= 1;
    }
}

thread_local! {
    /// The crew this thread works in.
    static CREW: RefCell<Option<Arc<Crew>>> = const { RefCell::new(None) };
}

/// Runs `body` as the first of a crew of `threads` threads. The fan-outs `body` starts
/// through [`crew_indexed_prefix`], [`crew_indexed`] and [`crew_map`], and the fan-outs
/// their tasks start in turn, share those threads: a fan-out runs on the thread that
/// starts it and borrows an idle thread of the crew for each further index it has yet to
/// start; a borrowed thread is idle again when it finds no index left, and a thread that
/// waits for the ones it borrowed gives up its own place meanwhile and takes a place again
/// before it goes on. At most `threads` threads are at work at a time, however deep the
/// nesting.
pub fn as_crew<T>(threads: usize, body: impl FnOnce() -> T) -> T {
    struct Restore(Option<Arc<Crew>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            CREW.set(self.0.take());
        }
    }
    let crew = Arc::new(Crew {
        idle: Mutex::new(threads.max(1) - 1),
        freed: Condvar::new(),
    });
    let _restore = Restore(CREW.replace(Some(crew)));
    body()
}

/// Runs `task` on the indices below `count`, on the crew of the calling thread (on that
/// thread alone when it works in none), up to the lowest index that fails: the results of
/// the indices before it, in index order, and its failure. After a failure no further
/// index is started; an index is started only once every index below it has been, so the
/// lowest failing index is the same at every thread count.
pub fn crew_indexed_prefix<T: Send, E: Send>(
    count: usize,
    task: impl Fn(usize) -> Result<T, E> + Sync,
) -> (Vec<T>, Option<E>) {
    let crew = CREW.with_borrow(Clone::clone);
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    // The first index that has not started, run to its end; `None` once there is none to
    // start.
    let step = || {
        if failed.load(Ordering::Relaxed) {
            return None;
        }
        let index = next.fetch_add(1, Ordering::Relaxed);
        if index >= count {
            return None;
        }
        let result = task(index);
        if result.is_err() {
            failed.store(true, Ordering::Relaxed);
        }
        Some((index, result))
    };
    let mut done: Vec<(usize, Result<T, E>)> = Vec::with_capacity(count);
    match &crew {
        None => done.extend(std::iter::from_fn(step)),
        Some(crew) => std::thread::scope(|scope| {
            let mut helpers = Vec::new();
            loop {
                // One thread for each index left beyond the one this thread takes next.
                let unstarted = count.saturating_sub(next.load(Ordering::Relaxed));
                for _ in 1..unstarted {
                    if failed.load(Ordering::Relaxed) || !crew.borrow() {
                        break;
                    }
                    helpers.push(scope.spawn(|| {
                        CREW.set(Some(Arc::clone(crew)));
                        let done: Vec<_> = std::iter::from_fn(step).collect();
                        crew.give();
                        done
                    }));
                }
                match step() {
                    Some(result) => done.push(result),
                    None => break,
                }
            }
            if !helpers.is_empty() {
                crew.give();
                for helper in helpers {
                    done.extend(joined(helper));
                }
                crew.take();
            }
        }),
    }
    done.sort_by_key(|(index, _)| *index);
    let mut results = Vec::with_capacity(done.len());
    for (_, result) in done {
        match result {
            Ok(value) => results.push(value),
            Err(failure) => return (results, Some(failure)),
        }
    }
    (results, None)
}

/// Runs `task` on every index below `count`, on the crew of the calling thread, and
/// returns the results in index order, or the failure of the lowest index that fails.
pub fn crew_indexed<T: Send, E: Send>(
    count: usize,
    task: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    match crew_indexed_prefix(count, task) {
        (results, None) => Ok(results),
        (_, Some(failure)) => Err(failure),
    }
}

/// The items one index of [`crew_map`] takes.
const ITEMS_A_RUN: usize = 1024;

/// `task` on every item, runs of items across the crew of the calling thread, the results
/// in the order of the items.
pub fn crew_map<I: Sync, T: Send>(items: &[I], task: impl Fn(&I) -> T + Sync) -> Vec<T> {
    if items.len() <= ITEMS_A_RUN {
        return items.iter().map(task).collect();
    }
    let runs: Vec<&[I]> = items.chunks(ITEMS_A_RUN).collect();
    let mapped = crew_indexed(runs.len(), |run| {
        Ok::<Vec<T>, Infallible>(runs[run].iter().map(&task).collect())
    });
    let Ok(mapped) = mapped;
    mapped.into_iter().flatten().collect()
}

/// Runs `task` on the indices below `count`, on at most `threads` threads, up to the lowest
/// index that fails: the results of the indices before it, in index order, and its failure.
/// After a failure no further index is started; an index is started only once every index
/// below it has been, so the lowest failing index is the same at every thread count.
pub fn run_indexed_prefix<T: Send, E: Send>(
    threads: usize,
    count: usize,
    task: impl Fn(usize) -> Result<T, E> + Sync,
) -> (Vec<T>, Option<E>) {
    as_crew(threads, || crew_indexed_prefix(count, task))
}

/// Runs `task` on every index below `count`, on at most `threads` threads, and returns the
/// results in index order, or the failure of the lowest index that fails.
pub fn run_indexed<T: Send, E: Send>(
    threads: usize,
    count: usize,
    task: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    match run_indexed_prefix(threads, count, task) {
        (results, None) => Ok(results),
        (_, Some(failure)) => Err(failure),
    }
}

/// What [`run_in_order`] holds between its threads.
struct InOrder<T> {
    /// The next index to start.
    next: usize,
    /// The indices taken so far.
    taken: usize,
    /// The results that wait to be taken.
    ready: BTreeMap<usize, T>,
    /// Whether the run ends before its last index: `take` broke off, or a thread ended on
    /// a panic of `task` and its index has no result.
    ended: bool,
}

/// Runs `task` on the indices below `count` on at most `threads` threads and hands each
/// result to `take` on the calling thread, in index order, as soon as every result before
/// it has been taken; the run ends at the first result on which `take` breaks off. No
/// index is started more than `ahead` indices past the last one taken, which bounds the
/// results held. A thread is started each time the calling thread has to wait for a
/// result, so a run whose `take` is slower than its tasks keeps to the threads that feed
/// it. With one thread, every index is run and taken in turn on the calling thread.
pub fn run_in_order<T: Send>(
    threads: usize,
    count: usize,
    ahead: usize,
    task: impl Fn(usize) -> T + Sync,
    mut take: impl FnMut(usize, T) -> ControlFlow<()>,
) {
    let workers = threads.min(count);
    if workers <= 1 {
        for index in 0..count {
            if take(index, task(index)).is_break() {
                break;
            }
        }
        return;
    }
    let ahead = ahead.max(1);
    let shared = Mutex::new(InOrder {
        next: 0,
        taken: 0,
        ready: BTreeMap::new(),
        ended: false,
    });
    let (room, arrived) = (Condvar::new(), Condvar::new());
    let held = || shared.lock().unwrap_or_else(PoisonError::into_inner);
    // A thread that ends on a panic of `task` ends the run, which would otherwise wait
    // for the result of its index.
    struct Leaving<'a, T>(&'a Mutex<InOrder<T>>, &'a Condvar, &'a Condvar);
    impl<T> Drop for Leaving<'_, T> {
        fn drop(&mut self) {
            if std::thread::panicking() {
                self.0.lock().unwrap_or_else(PoisonError::into_inner).ended = true;
                self.1.notify_all();
                self.2.notify_all();
            }
        }
    }
    let work = || {
        let _leaving = Leaving(&shared, &room, &arrived);
        loop {
            let index = {
                let mut state = held();
                while !state.ended && state.next < count && state.next >= state.taken + ahead {
                    state = room.wait(state).unwrap_or_else(PoisonError::into_inner);
                }
                if state.ended || state.next >= count {
                    return;
                }
                state.next += 1;
                if state.next == count {
                    // No index is left to start: the threads that wait for room end.
                    room.notify_all();
                }
                state.next - 1
            };
            let result = task(index);
            held().ready.insert(index, result);
            arrived.notify_one();
        }
    };
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for index in 0..count {
            let result = {
                let mut state = held();
                let mut grown = false;
                loop {
                    if let Some(result) = state.ready.remove(&index) {
                        state.taken = index + 1;
                        break Some(result);
                    }
                    if state.ended {
                        break None;
                    }
                    if !grown && handles.len() < workers {
                        grown = true;
                        drop(state);
                        handles.push(scope.spawn(work));
                        state = held();
                    } else {
                        state = arrived.wait(state).unwrap_or_else(PoisonError::into_inner);
                    }
                }
            };
            // One index more may start.
            room.notify_one();
            // Without a result, a thread ended on a panic, which its join carries on.
            if result.is_none_or(|result| take(index, result).is_break()) {
                held().ended = true;
                room.notify_all();
                break;
            }
        }
        handles.into_iter().for_each(joined);
    });
}

/// A file read once from end to end as runs of whole lines, each run parsed on whichever
/// thread is free.
pub trait BlockSource: Sync {
    /// What the read itself establishes.
    type Read: Send;
    /// What one run of lines parses to.
    type Parsed: Send;

    /// Reads the file and hands every run of whole lines, in order, to `emit`.
    fn read(&self, emit: &mut dyn FnMut(Vec<u8>)) -> Self::Read;

    /// Parses the `sequence`-th run of the file.
    fn parse(&self, sequence: u64, block: Vec<u8>) -> Self::Parsed;
}

/// The parsed runs of one source, each with its place in the file.
type Runs<P> = Mutex<Vec<(u64, P)>>;

/// What [`read_sources`] holds between its threads.
struct Handover {
    /// Runs that are read and wait to be parsed: `(source, sequence, run)`.
    queue: VecDeque<(usize, u64, Vec<u8>)>,
    /// Sources that have not been read to their end.
    unread: usize,
    /// Threads waiting for a run.
    idle: usize,
    /// Threads started.
    started: usize,
}

/// The sources of one [`read_sources`] and what its threads share.
struct Reading<'a, S: BlockSource> {
    sources: &'a [S],
    threads: usize,
    handover: Mutex<Handover>,
    wake: Condvar,
    /// The next source to read.
    next: AtomicUsize,
    parsed: Vec<Runs<S::Parsed>>,
    read: Vec<Mutex<Option<S::Read>>>,
}

impl<S: BlockSource> Reading<'_, S> {
    fn held(&self) -> std::sync::MutexGuard<'_, Handover> {
        self.handover.lock().unwrap_or_else(PoisonError::into_inner)
    }

    fn parse(&self, source: usize, sequence: u64, block: Vec<u8>) {
        let result = self.sources[source].parse(sequence, block);
        self.parsed[source]
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .push((sequence, result));
    }

    /// One thread of the reading, to its end. In a `scope`, a run that no thread waits for
    /// is handed to a thread started for it, while fewer threads than the bound have
    /// started.
    fn work<'scope>(&'scope self, scope: Option<&'scope std::thread::Scope<'scope, '_>>) {
        loop {
            let queued = self.held().queue.pop_front();
            if let Some((source, sequence, block)) = queued {
                self.parse(source, sequence, block);
                continue;
            }
            let source = self.next.fetch_add(1, Ordering::Relaxed);
            if source < self.sources.len() {
                let mut sequence = 0;
                let outcome = self.sources[source].read(&mut |block| {
                    let mut handover = self.held();
                    if handover.idle > handover.queue.len() {
                        handover.queue.push_back((source, sequence, block));
                        drop(handover);
                        self.wake.notify_one();
                    } else if let Some(scope) = scope.filter(|_| handover.started < self.threads) {
                        handover.started += 1;
                        handover.queue.push_back((source, sequence, block));
                        drop(handover);
                        scope.spawn(move || self.work(Some(scope)));
                    } else {
                        drop(handover);
                        self.parse(source, sequence, block);
                    }
                    sequence += 1;
                });
                *self.read[source]
                    .lock()
                    .unwrap_or_else(PoisonError::into_inner) = Some(outcome);
                let mut handover = self.held();
                handover.unread -= 1;
                if handover.unread == 0 {
                    drop(handover);
                    self.wake.notify_all();
                }
                continue;
            }
            let mut handover = self.held();
            loop {
                if let Some((source, sequence, block)) = handover.queue.pop_front() {
                    drop(handover);
                    self.parse(source, sequence, block);
                    break;
                }
                if handover.unread == 0 {
                    return;
                }
                handover.idle += 1;
                handover = self
                    .wake
                    .wait(handover)
                    .unwrap_or_else(PoisonError::into_inner);
                handover.idle -= 1;
            }
        }
    }
}

/// Reads every source on at most `threads` threads and returns, per source, what its read
/// established and its parsed runs in the order of the file. One thread starts for each
/// source, up to the bound, and one more for each run read while no thread waits for one,
/// up to the bound. A thread with no source left to read parses the runs the reading
/// threads hand over; at the bound, a reading thread parses its own run when no thread
/// waits for one.
pub fn read_sources<S: BlockSource>(
    threads: usize,
    sources: &[S],
) -> Vec<(S::Read, Vec<S::Parsed>)> {
    let reading = Reading {
        sources,
        threads: threads.max(1),
        handover: Mutex::new(Handover {
            queue: VecDeque::new(),
            unread: sources.len(),
            idle: 0,
            started: 0,
        }),
        wake: Condvar::new(),
        next: AtomicUsize::new(0),
        parsed: sources.iter().map(|_| Mutex::new(Vec::new())).collect(),
        read: sources.iter().map(|_| Mutex::new(None)).collect(),
    };
    if reading.threads == 1 {
        reading.work(None);
    } else {
        std::thread::scope(|scope| {
            let first = reading.threads.min(sources.len()).max(1);
            reading.held().started = first;
            let handles: Vec<_> = (0..first)
                .map(|_| scope.spawn(|| reading.work(Some(scope))))
                .collect();
            handles.into_iter().for_each(joined);
        });
    }
    let Reading { read, parsed, .. } = reading;
    read.into_iter()
        .zip(parsed)
        .map(|(read, parsed)| {
            let read = read
                .into_inner()
                .unwrap_or_else(PoisonError::into_inner)
                .expect("every source is read");
            let mut parsed = parsed.into_inner().unwrap_or_else(PoisonError::into_inner);
            parsed.sort_by_key(|(sequence, _)| *sequence);
            (read, parsed.into_iter().map(|(_, run)| run).collect())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A source of `runs` runs, the `n`-th holding the byte `n`.
    struct Counting {
        runs: u8,
    }

    impl BlockSource for Counting {
        type Read = u8;
        type Parsed = (u64, u8);

        fn read(&self, emit: &mut dyn FnMut(Vec<u8>)) -> u8 {
            for n in 0..self.runs {
                emit(vec![n]);
            }
            self.runs
        }

        fn parse(&self, sequence: u64, block: Vec<u8>) -> (u64, u8) {
            (sequence, block[0])
        }
    }

    #[test]
    fn every_source_comes_back_whole_and_in_order_for_every_thread_count() {
        let sources: Vec<Counting> = [0, 1, 200, 3, 17, 0, 64]
            .into_iter()
            .map(|runs| Counting { runs })
            .collect();
        for threads in [1, 2, 5, 32] {
            let read = read_sources(threads, &sources);
            assert_eq!(read.len(), sources.len());
            for (source, (runs, parsed)) in sources.iter().zip(read) {
                assert_eq!(runs, source.runs);
                let expected: Vec<(u64, u8)> =
                    (0..source.runs).map(|n| (u64::from(n), n)).collect();
                assert_eq!(parsed, expected, "threads {threads}");
            }
        }
        assert!(read_sources::<Counting>(4, &[]).is_empty());
    }

    #[test]
    fn results_come_back_in_index_order_for_every_thread_count() {
        for threads in [1, 2, 7, 64] {
            let squares = run_indexed(threads, 100, |i| Ok::<_, Failure>(i * i)).unwrap();
            assert_eq!(squares, (0..100).map(|i| i * i).collect::<Vec<_>>());
        }
        assert!(run_indexed(4, 0, Ok::<usize, Failure>).unwrap().is_empty());
    }

    #[test]
    fn the_lowest_failure_is_returned_with_the_results_before_it() {
        for threads in [1, 3, 16] {
            let task = |i: usize| {
                if i == 20 || i == 35 {
                    Err(format!("index {i}"))
                } else {
                    Ok(i)
                }
            };
            let (before, failure) = run_indexed_prefix(threads, 50, task);
            assert_eq!(before, (0..20).collect::<Vec<_>>());
            assert_eq!(failure.as_deref(), Some("index 20"));
            assert_eq!(run_indexed(threads, 50, task).unwrap_err(), "index 20");
        }
    }

    #[test]
    fn nested_fan_outs_return_the_same_at_every_crew_size() {
        let expected: Vec<Vec<usize>> = (0..40)
            .map(|outer| (0..2700).map(|inner| outer * 10_000 + inner).collect())
            .collect();
        for threads in [1, 2, 3, 8, 64] {
            let nested = as_crew(threads, || {
                crew_indexed(40, |outer| {
                    let items: Vec<usize> = (0..2700).collect();
                    Ok::<_, Failure>(crew_map(&items, |inner| outer * 10_000 + inner))
                })
            });
            assert_eq!(nested.unwrap(), expected, "threads {threads}");
        }
        let items: Vec<usize> = (0..3000).collect();
        let alone = crew_map(&items, |item| item + 1);
        assert_eq!(alone, (1..=3000).collect::<Vec<_>>());
    }

    #[test]
    fn a_crew_keeps_no_more_threads_at_work_than_it_has() {
        for threads in [1, 2, 4, 9] {
            let (at_work, most) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let busy = || {
                let now = at_work.fetch_add(1, Ordering::SeqCst) + 1;
                most.fetch_max(now, Ordering::SeqCst);
                std::thread::sleep(std::time::Duration::from_millis(2));
                at_work.fetch_sub(1, Ordering::SeqCst);
            };
            as_crew(threads, || {
                crew_indexed(6, |_| {
                    busy();
                    crew_indexed(6, |_| {
                        busy();
                        Ok::<_, Failure>(())
                    })
                })
            })
            .unwrap();
            let most = most.load(Ordering::SeqCst);
            assert!(most <= threads, "{most} at work in a crew of {threads}");
        }
    }

    #[test]
    fn the_lowest_failure_of_a_crew_is_returned_with_the_results_before_it() {
        for threads in [1, 3, 16] {
            let (before, failure) = as_crew(threads, || {
                crew_indexed_prefix(50, |i| {
                    if i == 20 || i == 35 {
                        Err(format!("index {i}"))
                    } else {
                        Ok(i)
                    }
                })
            });
            assert_eq!(before, (0..20).collect::<Vec<_>>());
            assert_eq!(failure.as_deref(), Some("index 20"));
        }
    }

    #[test]
    fn results_are_taken_in_index_order_and_none_is_started_too_far_ahead() {
        for (threads, ahead) in [(1, 1), (2, 1), (4, 3), (16, 5), (8, 1000)] {
            let taking = AtomicUsize::new(0);
            let mut seen = Vec::new();
            run_in_order(
                threads,
                200,
                ahead,
                |index| {
                    // `taking` is at most two behind the count of results taken.
                    let taking = taking.load(Ordering::SeqCst);
                    assert!(index <= taking + ahead + 1, "{index} started at {taking}");
                    index * 3
                },
                |index, result| {
                    taking.store(index, Ordering::SeqCst);
                    seen.push((index, result));
                    ControlFlow::Continue(())
                },
            );
            let expected: Vec<(usize, usize)> = (0..200).map(|index| (index, index * 3)).collect();
            assert_eq!(seen, expected, "threads {threads}");
        }
        run_in_order(4, 0, 2, |index| index, |_, _| panic!("nothing to take"));
    }

    #[test]
    fn a_take_that_breaks_off_ends_the_run_at_that_result() {
        for threads in [1, 2, 6, 32] {
            let mut seen = Vec::new();
            run_in_order(
                threads,
                500,
                4,
                |index| index,
                |index, result| {
                    seen.push(result);
                    match index {
                        40 => ControlFlow::Break(()),
                        _ => ControlFlow::Continue(()),
                    }
                },
            );
            assert_eq!(seen, (0..=40).collect::<Vec<_>>(), "threads {threads}");
        }
    }

    #[test]
    fn a_task_that_panics_ends_the_run_in_order() {
        let ended = std::panic::catch_unwind(|| {
            run_in_order(
                4,
                50,
                2,
                |index| assert_ne!(index, 7, "the task of index 7"),
                |_, ()| ControlFlow::Continue(()),
            );
        });
        assert!(ended.is_err());
    }

    #[test]
    fn the_bound_is_a_whole_number_of_one_or_more() {
        assert_eq!(threads_from("1").unwrap(), 1);
        assert_eq!(threads_from(" 12 ").unwrap(), 12);
        for refused in ["0", "-1", "1.5", "many", " "] {
            let failure = threads_from(refused).unwrap_err();
            assert_eq!(
                failure.message,
                format!(
                    "VEP_SV_COMPARE_THREADS={refused}: expected a whole number of threads, 1 or more"
                )
            );
        }
    }
}
