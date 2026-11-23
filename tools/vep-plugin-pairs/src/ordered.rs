// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Work made on several threads and taken in the order of its items.

use std::collections::BTreeMap;
use std::sync::{Condvar, Mutex};

/// The items a thread may stand ahead of the taking, made or in the making.
const AHEAD: usize = 4;

struct Made<T> {
    /// The next item a thread makes.
    next: usize,
    /// The items taken.
    taken: usize,
    /// The items made and not taken; `None` for one whose thread ended before it was made.
    items: BTreeMap<usize, Option<T>>,
    /// Whether the taking has ended before the last item.
    ended: bool,
}

struct Board<T> {
    made: Mutex<Made<T>>,
    moved: Condvar,
}

/// One item a thread has taken up, posted when the value is dropped.
struct Posting<'b, T> {
    board: &'b Board<T>,
    at: usize,
    item: Option<T>,
}

impl<T> Drop for Posting<'_, T> {
    fn drop(&mut self) {
        if let Ok(mut made) = self.board.made.lock() {
            made.items.insert(self.at, self.item.take());
        }
        self.board.moved.notify_all();
    }
}

/// Ends the taking when it is dropped.
struct Ending<'b, T>(&'b Board<T>);

impl<T> Drop for Ending<'_, T> {
    fn drop(&mut self) {
        if let Ok(mut made) = self.0.made.lock() {
            made.ended = true;
        }
        self.0.moved.notify_all();
    }
}

/// Makes the items `0..count` and hands each to `take` on the calling thread, in order.
/// The items are made on `threads` threads, each taking up the next item no thread has;
/// with one thread, or one item, they are made on the calling thread. A thread takes up an
/// item only while fewer than four items a thread stand made and not taken, so the items
/// held at a time are bounded. `take` returns whether the work goes on: after a `false` no
/// further item is taken up or handed over.
pub fn in_order<T: Send>(
    count: usize,
    threads: usize,
    make: &(dyn Fn(usize) -> T + Sync),
    take: &mut dyn FnMut(T) -> bool,
) {
    let threads = threads.min(count);
    if threads <= 1 {
        for at in 0..count {
            if !take(make(at)) {
                return;
            }
        }
        return;
    }
    let board = Board {
        made: Mutex::new(Made {
            next: 0,
            taken: 0,
            items: BTreeMap::new(),
            ended: false,
        }),
        moved: Condvar::new(),
    };
    std::thread::scope(|scope| {
        // The scope waits for its threads on every way out, a panic of `take` included, and
        // a thread that waits for an item to be taken ends only once the taking has.
        let _ending = Ending(&board);
        for _ in 0..threads {
            scope.spawn(|| loop {
                let at = {
                    let mut made = board.made.lock().expect("lock");
                    while !made.ended
                        && made.next < count
                        && made.next >= made.taken + AHEAD * threads
                    {
                        made = board.moved.wait(made).expect("lock");
                    }
                    if made.ended || made.next >= count {
                        return;
                    }
                    made.next += 1;
                    made.next - 1
                };
                let mut posting = Posting {
                    board: &board,
                    at,
                    item: None,
                };
                posting.item = Some(make(at));
            });
        }
        for at in 0..count {
            let item = {
                let mut made = board.made.lock().expect("lock");
                loop {
                    if let Some(item) = made.items.remove(&at) {
                        made.taken = at + 1;
                        break item;
                    }
                    made = board.moved.wait(made).expect("lock");
                }
            };
            board.moved.notify_all();
            let goes_on = match item {
                Some(item) => take(item),
                None => false,
            };
            if !goes_on {
                board.made.lock().expect("lock").ended = true;
                board.moved.notify_all();
                break;
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::mpsc;
    use std::time::Duration;

    #[test]
    fn panic_of_the_taking_ends_the_threads_and_reaches_the_caller() {
        let (ended, waited) = mpsc::channel();
        std::thread::spawn(move || {
            let taken = std::panic::catch_unwind(|| {
                in_order(1000, 2, &|at| at, &mut |item| {
                    assert!(item < 3, "the taking panics at its fourth item");
                    true
                });
            });
            let _ = ended.send(taken.is_err());
        });
        assert_eq!(waited.recv_timeout(Duration::from_secs(60)), Ok(true));
    }

    #[test]
    fn items_are_taken_in_their_order_on_any_thread_count() {
        for threads in [0, 1, 2, 5, 16] {
            for count in [0, 1, 3, 100] {
                let mut taken = Vec::new();
                in_order(count, threads, &|at| at * at, &mut |item| {
                    taken.push(item);
                    true
                });
                let expected: Vec<usize> = (0..count).map(|at| at * at).collect();
                assert_eq!(taken, expected, "{threads} threads, {count} items");
            }
        }
    }

    #[test]
    fn nothing_is_taken_after_the_taking_ends() {
        for threads in [1, 2, 5, 16] {
            let made = AtomicUsize::new(0);
            let mut taken = Vec::new();
            in_order(
                1000,
                threads,
                &|at| {
                    made.fetch_add(1, Ordering::SeqCst);
                    at
                },
                &mut |item| {
                    taken.push(item);
                    item < 7
                },
            );
            assert_eq!(taken, [0, 1, 2, 3, 4, 5, 6, 7], "{threads} threads");
            // No thread stands further ahead of the taking than its bound.
            assert!(
                made.load(Ordering::SeqCst) <= 8 + AHEAD * threads + threads,
                "{threads} threads"
            );
        }
    }
}
