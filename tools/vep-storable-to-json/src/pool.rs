// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The host's cores across shards, with results taken in order.
//!
//! Workers take tasks in index order, one shard a worker at a time, and a single
//! consumer takes the outcomes by index, so what the run writes and prints follows the
//! shards' order whatever order they finish in. A worker does not start a task more
//! than a window ahead of the consumer, which bounds the outcomes held in memory.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};

struct State<R> {
    done: BTreeMap<usize, R>,
    /// The index the consumer takes next.
    wanted: usize,
}

/// A run of `tasks` through `work` on `threads` workers.
pub struct Ordered<'a, T, R> {
    tasks: &'a [T],
    next: AtomicUsize,
    stop: AtomicBool,
    window: usize,
    state: Mutex<State<R>>,
    ready: Condvar,
    room: Condvar,
}

impl<'a, T: Sync, R: Send> Ordered<'a, T, R> {
    pub fn new(tasks: &'a [T], threads: usize) -> Self {
        Self {
            tasks,
            next: AtomicUsize::new(0),
            stop: AtomicBool::new(false),
            window: threads.max(1) * 2,
            state: Mutex::new(State {
                done: BTreeMap::new(),
                wanted: 0,
            }),
            ready: Condvar::new(),
            room: Condvar::new(),
        }
    }

    /// One worker's loop; returns when the tasks are done or the consumer stopped.
    pub fn work(&self, work: impl Fn(&T) -> R) {
        loop {
            let index = self.next.fetch_add(1, Ordering::SeqCst);
            if index >= self.tasks.len() {
                return;
            }
            {
                let mut state = self.state.lock().expect("pool state");
                while index >= state.wanted + self.window && !self.stop.load(Ordering::SeqCst) {
                    state = self.room.wait(state).expect("pool state");
                }
            }
            if self.stop.load(Ordering::SeqCst) {
                return;
            }
            let outcome = work(&self.tasks[index]);
            let mut state = self.state.lock().expect("pool state");
            state.done.insert(index, outcome);
            self.ready.notify_all();
        }
    }

    /// The outcome of task `index`; the consumer takes every index once, in order.
    pub fn take(&self, index: usize) -> R {
        let mut state = self.state.lock().expect("pool state");
        state.wanted = index;
        self.room.notify_all();
        loop {
            if let Some(outcome) = state.done.remove(&index) {
                state.wanted = index + 1;
                self.room.notify_all();
                return outcome;
            }
            state = self.ready.wait(state).expect("pool state");
        }
    }

    /// Tells the workers to start nothing more.
    pub fn stop(&self) {
        self.stop.store(true, Ordering::SeqCst);
        let _state = self.state.lock().expect("pool state");
        self.room.notify_all();
    }
}

/// Runs `consume` on this thread while `threads` workers run `work` over `tasks`;
/// `consume` takes outcomes from the pool it is given. The workers are stopped when
/// `consume` returns.
pub fn run_ordered<T: Sync, R: Send, O>(
    tasks: &[T],
    threads: usize,
    work: impl Fn(&T) -> R + Sync,
    consume: impl FnOnce(&Ordered<'_, T, R>) -> O,
) -> O {
    let pool = Ordered::new(tasks, threads);
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1).min(tasks.len().max(1)) {
            scope.spawn(|| pool.work(&work));
        }
        let output = consume(&pool);
        pool.stop();
        output
    })
}
