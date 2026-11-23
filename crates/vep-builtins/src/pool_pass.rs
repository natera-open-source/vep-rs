// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One pass of a batch spread over the current rayon pool, with the warnings
//! its tasks raise handed to the caller in task order.
//!
//! A run holds the warnings a batch's passes raise and reports them with the
//! batch (`vep_core::warning_hold`). A hold belongs to one thread, and a pool
//! thread that waits on its own tasks runs other tasks meanwhile, another
//! batch's among them. So each task of a pass holds its own warnings, and the
//! thread that started the tasks takes them once all have ended, in the order
//! of the items and not the order the threads reached them. A pass that reads
//! ahead for a batch keeps them with what it read ([`collect_in_order`]); a
//! pass over the batch itself raises them where the batch's hold takes them
//! ([`map_in_order`]). A pass that works record by record names each record
//! on the warnings its work raises ([`about_record`]).

use rayon::prelude::*;
use vep_core::warning_hold::{self, Held};

/// `work` over every item, in parallel on the current rayon pool: the results
/// in item order, and the warnings the tasks raised, in item order, whatever
/// thread ran which task.
pub(crate) fn collect_in_order<I, R>(
    items: I,
    work: impl Fn(I::Item) -> R + Sync + Send,
) -> (Vec<R>, Vec<Held>)
where
    I: IntoParallelIterator,
    R: Send,
{
    let done: Vec<(R, Vec<Held>)> = items
        .into_par_iter()
        .map(|item| warning_hold::holding(|| work(item)))
        .collect();
    let mut results = Vec::with_capacity(done.len());
    let mut raised = Vec::new();
    for (result, held) in done {
        results.push(result);
        raised.extend(held);
    }
    (results, raised)
}

/// [`collect_in_order`], with the warnings the tasks raised joining the hold
/// of the calling thread once every task has ended ([`raise_again`]).
pub(crate) fn map_in_order<I, R>(items: I, work: impl Fn(I::Item) -> R + Sync + Send) -> Vec<R>
where
    I: IntoParallelIterator,
    R: Send,
{
    let (results, raised) = collect_in_order(items, work);
    raise_again(raised);
    results
}

/// Runs `work`, the part of a pass that is about one input record, and hands
/// on the warnings it raises ([`raise_again`]), those that name no record
/// under `input_record`.
pub(crate) fn about_record<T>(input_record: u64, work: impl FnOnce() -> T) -> T {
    let (result, mut raised) = warning_hold::holding(work);
    for held in &mut raised {
        if held.input_record == 0 {
            held.input_record = input_record;
        }
    }
    raise_again(raised);
    result
}

/// Hands `warnings` to the hold the calling thread is inside, in order; a
/// thread inside none raises each again as a `warn!` event, under the record
/// it names.
pub(crate) fn raise_again(warnings: Vec<Held>) {
    for Held { input_record, text } in warning_hold::hand_on(warnings) {
        if input_record == 0 {
            tracing::warn!("{text}");
        } else {
            tracing::warn!(input_record, "{text}");
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::time::Duration;

    use super::*;
    use crate::test_fixtures::Warnings;

    fn pool(threads: usize) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
    }

    /// The task of item `n`: the earlier the item the longer the task runs, so
    /// the tasks end in the reverse of item order on a pool of several
    /// threads. It raises one warning naming record `n + 1` and one naming
    /// none, as the layer that writes a run's warnings offers them.
    fn task(n: u64) -> u64 {
        std::thread::sleep(Duration::from_millis(5 * (8 - n)));
        warning_hold::offer(n + 1, format!("record of task {n}")).unwrap();
        warning_hold::offer(0, format!("task {n}")).unwrap();
        n * 2
    }

    fn expected() -> Vec<Held> {
        (0..8u64)
            .flat_map(|n| {
                [
                    Held {
                        input_record: n + 1,
                        text: format!("record of task {n}"),
                    },
                    Held {
                        input_record: 0,
                        text: format!("task {n}"),
                    },
                ]
            })
            .collect()
    }

    /// The results come back in item order, and the warnings the tasks raised
    /// join the caller's hold in item order after those the caller raised
    /// itself, whatever thread ran which task and in whatever order they
    /// ended.
    #[test]
    fn the_tasks_warnings_join_the_callers_hold_in_item_order() {
        for threads in [1, 2, 8] {
            let (results, raised) = pool(threads).install(|| {
                warning_hold::holding(|| {
                    warning_hold::offer(0, "before the pass".to_string()).unwrap();
                    map_in_order(0..8u64, task)
                })
            });
            assert_eq!(results, [0, 2, 4, 6, 8, 10, 12, 14], "{threads} threads");
            let mut all = vec![Held {
                input_record: 0,
                text: "before the pass".to_string(),
            }];
            all.extend(expected());
            assert_eq!(raised, all, "{threads} threads");
        }
    }

    /// A task that runs on a thread waiting inside another task's hold keeps
    /// its warnings from that hold: a pass inside a pass hands each inner
    /// task's warnings to the outer task that started it.
    #[test]
    fn a_task_leaves_no_warning_in_another_tasks_hold() {
        let (per_outer, raised) = pool(4).install(|| {
            warning_hold::holding(|| {
                map_in_order(0..4u64, |outer| {
                    let ((), own) = warning_hold::holding(|| {
                        map_in_order(0..4u64, |inner| {
                            std::thread::sleep(Duration::from_millis(2));
                            warning_hold::offer(outer + 1, format!("{outer}.{inner}")).unwrap();
                        });
                    });
                    own
                })
            })
        });
        assert!(raised.is_empty());
        for (outer, own) in per_outer.into_iter().enumerate() {
            let texts: Vec<String> = own.iter().map(|held| held.text.clone()).collect();
            let expected: Vec<String> = (0..4).map(|inner| format!("{outer}.{inner}")).collect();
            assert_eq!(texts, expected, "outer task {outer}");
            assert!(own.iter().all(|held| held.input_record == outer as u64 + 1));
        }
    }

    /// A pass that keeps its tasks' warnings gets them in item order and
    /// leaves none in the hold of the thread that started it.
    #[test]
    fn a_pass_that_keeps_its_warnings_leaves_none_in_the_callers_hold() {
        let ((results, kept), raised) =
            pool(8).install(|| warning_hold::holding(|| collect_in_order(0..8u64, task)));
        assert_eq!(results, [0, 2, 4, 6, 8, 10, 12, 14]);
        assert_eq!(kept, expected());
        assert!(raised.is_empty());
    }

    /// Work about one record hands on what it raises in the order raised: a
    /// warning that names no record under the work's record, one that names a
    /// record under its own.
    #[test]
    fn work_about_a_record_names_it_on_the_warnings_that_name_none() {
        let (answer, raised) = warning_hold::holding(|| {
            about_record(7, || {
                warning_hold::offer(0, "names none".to_string()).unwrap();
                warning_hold::offer(3, "names its own".to_string()).unwrap();
                42
            })
        });
        assert_eq!(answer, 42);
        assert_eq!(
            raised,
            [
                Held {
                    input_record: 7,
                    text: "names none".to_string(),
                },
                Held {
                    input_record: 3,
                    text: "names its own".to_string(),
                },
            ]
        );
    }

    /// A caller inside no hold raises the tasks' warnings again as events, in
    /// item order.
    #[test]
    fn a_caller_inside_no_hold_raises_the_warnings_again() {
        let warnings = Arc::new(Warnings::default());
        let single = pool(1);
        single.install(|| {
            tracing::subscriber::with_default(Arc::clone(&warnings), || {
                map_in_order(0..8u64, task);
            })
        });
        let texts: Vec<String> = expected().into_iter().map(|held| held.text).collect();
        assert_eq!(*warnings.0.lock().unwrap(), texts);
    }
}
