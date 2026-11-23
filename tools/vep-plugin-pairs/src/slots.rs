// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The threads of a run that work at one time. A thread takes a slot before it works and
//! gives it back before it waits for another thread, so that a run works on as many
//! threads as it was given however many it starts. A slot goes to the thread that has
//! waited longest.

use std::collections::VecDeque;
use std::sync::{Condvar, Mutex};

struct State {
    free: usize,
    /// The takers that wait, the one that asked first in front.
    waiting: VecDeque<u64>,
    /// The number the next taker waits under.
    next: u64,
}

/// The slots of one run.
pub struct Slots {
    state: Mutex<State>,
    changed: Condvar,
}

impl Slots {
    /// `count` slots, one where `count` is zero.
    pub fn new(count: usize) -> Slots {
        Slots {
            state: Mutex::new(State {
                free: count.max(1),
                waiting: VecDeque::new(),
                next: 0,
            }),
            changed: Condvar::new(),
        }
    }

    /// Takes a slot, after every thread that asked before this one, and holds it until the
    /// returned value is dropped.
    pub fn take(&self) -> Slot<'_> {
        let mut state = self.state.lock().unwrap();
        let turn = state.next;
        state.next += 1;
        state.waiting.push_back(turn);
        while state.free == 0 || state.waiting.front() != Some(&turn) {
            state = self.changed.wait(state).unwrap();
        }
        state.waiting.pop_front();
        state.free -= 1;
        let more = state.free > 0 && !state.waiting.is_empty();
        drop(state);
        if more {
            self.changed.notify_all();
        }
        Slot { slots: self }
    }

    fn give_back(&self) {
        self.state.lock().unwrap().free += 1;
        self.changed.notify_all();
    }
}

/// One slot, held by the thread that took it.
pub struct Slot<'a> {
    slots: &'a Slots,
}

impl<'a> Slot<'a> {
    /// Gives the slot back until the returned value is dropped, when one is taken again:
    /// for a thread that does nothing meanwhile but wait for other threads.
    pub fn lent(&mut self) -> Lent<'_, 'a> {
        self.slots.give_back();
        Lent { slot: self }
    }
}

/// A slot its thread gave back for a wait.
pub struct Lent<'s, 'a> {
    slot: &'s mut Slot<'a>,
}

impl Drop for Lent<'_, '_> {
    fn drop(&mut self) {
        // The slot taken is the one `slot` stands for again, and gives back when it is dropped.
        std::mem::forget(self.slot.slots.take());
    }
}

impl Drop for Slot<'_> {
    fn drop(&mut self) {
        self.slots.give_back();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{mpsc, Barrier};
    use std::time::Duration;

    #[test]
    fn no_more_threads_hold_a_slot_than_the_run_has() {
        for count in [1, 2, 3] {
            let slots = Slots::new(count);
            let (holding, most) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let started = Barrier::new(8);
            std::thread::scope(|scope| {
                for _ in 0..8 {
                    scope.spawn(|| {
                        started.wait();
                        for _ in 0..50 {
                            let _slot = slots.take();
                            let held = holding.fetch_add(1, Ordering::SeqCst) + 1;
                            most.fetch_max(held, Ordering::SeqCst);
                            std::thread::yield_now();
                            holding.fetch_sub(1, Ordering::SeqCst);
                        }
                    });
                }
            });
            assert!((1..=count).contains(&most.load(Ordering::SeqCst)));
            assert_eq!(slots.state.lock().unwrap().free, count);
        }
        let single = Slots::new(0);
        drop(single.take().lent());
        assert_eq!(single.state.lock().unwrap().free, 1);
    }

    #[test]
    fn slot_goes_to_the_thread_that_waited_longest() {
        let slots = Slots::new(1);
        let held = slots.take();
        let (served, order) = mpsc::channel();
        std::thread::scope(|scope| {
            for taker in 0..4 {
                let served = served.clone();
                let slots = &slots;
                scope.spawn(move || {
                    let _slot = slots.take();
                    served.send(taker).unwrap();
                });
                // The taker waits in line before the next one is started.
                while slots.state.lock().unwrap().waiting.len() <= taker {
                    std::thread::sleep(Duration::from_millis(1));
                }
            }
            drop(held);
        });
        drop(served);
        assert_eq!(order.iter().collect::<Vec<_>>(), [0, 1, 2, 3]);
    }

    #[test]
    fn lent_slot_is_worked_on_by_another_thread_and_taken_again() {
        let slots = Slots::new(1);
        let mut own = slots.take();
        let worked = AtomicUsize::new(0);
        let lent = own.lent();
        std::thread::scope(|scope| {
            scope.spawn(|| {
                let _slot = slots.take();
                worked.fetch_add(1, Ordering::SeqCst);
            });
        });
        drop(lent);
        assert_eq!(worked.load(Ordering::SeqCst), 1);
        assert_eq!(slots.state.lock().unwrap().free, 0);
        drop(own);
        assert_eq!(slots.state.lock().unwrap().free, 1);
    }
}
