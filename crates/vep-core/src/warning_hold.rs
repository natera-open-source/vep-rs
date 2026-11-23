// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The hold a thread places on the warnings it raises, so that a run reports
//! them where the work stands among its warnings, in an order no thread's
//! timing decides.
//!
//! Work that wants its warnings back runs inside [`holding`]; the layer that
//! writes a run's warnings gives each one to the hold of the thread raising it
//! ([`offer`]) and writes it only when that thread is inside none. A warning
//! names the input record it is about through the event field
//! [`RECORD_FIELD`], which travels with it as [`Held::input_record`], so the
//! code that reports a batch's held warnings places each with its record.
//!
//! A hold belongs to one thread, and a pool thread that waits on tasks of its
//! own runs other tasks meanwhile, so a pass that spreads over a pool has each
//! task hold its own warnings and hands them to the hold of the thread that
//! started the tasks, in task order ([`hand_on`]).

use std::cell::RefCell;

/// The name of the event field through which a `warn!` names the input record
/// it is about, [`InputVariant::input_record`]: `warn!(input_record = n, ..)`.
/// The layer that writes a run's warnings reads the field as
/// [`Held::input_record`] and leaves it out of the warning's text.
///
/// [`InputVariant::input_record`]: crate::variant::InputVariant::input_record
pub const RECORD_FIELD: &str = "input_record";

/// One warning a thread raised inside a hold.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Held {
    /// The input record the warning is about, 0 when the code raising it named
    /// none.
    pub input_record: u64,
    /// The warning's text.
    pub text: String,
}

thread_local! {
    /// The warnings this thread has raised inside a [`holding`] call that has
    /// not ended; `None` outside one.
    static HELD: RefCell<Option<Vec<Held>>> = const { RefCell::new(None) };
}

/// Puts back the warnings held outside a [`holding`] call when the call ends,
/// by returning or by unwinding.
struct OuterHold(Option<Vec<Held>>);

impl Drop for OuterHold {
    fn drop(&mut self) {
        let outer = self.0.take();
        HELD.with(|held| *held.borrow_mut() = outer);
    }
}

/// Runs `work` with every warning this thread raises kept from the run's
/// warnings, and returns those warnings, in the order raised, beside what
/// `work` returned. The caller reports them where the work stands among the
/// run's warnings. A warning another thread raises is not held; a call inside
/// another holds its own warnings and leaves the outer call's as they were;
/// the warnings of work that unwinds are dropped with it.
pub fn holding<T>(work: impl FnOnce() -> T) -> (T, Vec<Held>) {
    let outer = OuterHold(HELD.with(|held| held.borrow_mut().replace(Vec::new())));
    let result = work();
    let raised = HELD
        .with(|held| held.borrow_mut().take())
        .unwrap_or_default();
    drop(outer);
    (result, raised)
}

/// Gives the warning `text`, about `input_record` (0 for none), to the hold
/// this thread is inside; outside a hold the text comes back for the caller to
/// write.
pub fn offer(input_record: u64, text: String) -> Result<(), String> {
    HELD.with(|held| match held.borrow_mut().as_mut() {
        Some(raised) => {
            raised.push(Held { input_record, text });
            Ok(())
        }
        None => Err(text),
    })
}

/// Gives `warnings`, in order, to the hold this thread is inside and returns
/// none; outside a hold they all come back for the caller to raise again.
pub fn hand_on(warnings: Vec<Held>) -> Vec<Held> {
    HELD.with(|held| match held.borrow_mut().as_mut() {
        Some(raised) => {
            raised.extend(warnings);
            Vec::new()
        }
        None => warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn held(input_record: u64, text: &str) -> Held {
        Held {
            input_record,
            text: text.to_string(),
        }
    }

    /// A warning offered inside a hold comes back from the call, in the order
    /// offered and with the record it names; one offered outside a hold comes
    /// back to the caller.
    #[test]
    fn a_hold_returns_the_warnings_offered_inside_it() {
        assert_eq!(offer(0, "before".to_string()), Err("before".to_string()));
        let (answer, raised) = holding(|| {
            assert_eq!(offer(0, "first".to_string()), Ok(()));
            assert_eq!(offer(7, "second".to_string()), Ok(()));
            42
        });
        assert_eq!(answer, 42);
        assert_eq!(raised, [held(0, "first"), held(7, "second")]);
        assert_eq!(offer(3, "after".to_string()), Err("after".to_string()));
    }

    /// A hold inside another keeps its own warnings and leaves the outer
    /// hold's as they were, also when the inner work unwinds.
    #[test]
    fn a_hold_inside_another_keeps_its_own() {
        let ((), outer) = holding(|| {
            offer(1, "outer first".to_string()).unwrap();
            let ((), inner) = holding(|| offer(2, "inner".to_string()).unwrap());
            assert_eq!(inner, [held(2, "inner")]);
            let unwound = std::panic::catch_unwind(|| {
                holding(|| {
                    offer(3, "dropped".to_string()).unwrap();
                    panic!("the inner work unwinds");
                })
            });
            assert!(unwound.is_err());
            offer(1, "outer second".to_string()).unwrap();
        });
        assert_eq!(outer, [held(1, "outer first"), held(1, "outer second")]);
        assert_eq!(offer(0, "after".to_string()), Err("after".to_string()));
    }

    /// A warning another thread offers is not held by this thread's hold.
    #[test]
    fn a_hold_takes_no_warning_of_another_thread() {
        let ((), raised) = holding(|| {
            offer(0, "this thread".to_string()).unwrap();
            let other = std::thread::spawn(|| offer(0, "another thread".to_string()))
                .join()
                .unwrap();
            assert_eq!(other, Err("another thread".to_string()));
        });
        assert_eq!(raised, [held(0, "this thread")]);
    }

    /// Warnings handed on inside a hold join it in the order given, after
    /// those the thread has raised; outside a hold they all come back.
    #[test]
    fn handed_on_warnings_join_the_hold_in_order() {
        let batch = vec![held(4, "of a task"), held(0, "of the next task")];
        assert_eq!(hand_on(batch.clone()), batch);
        let ((), raised) = holding(|| {
            offer(9, "raised here".to_string()).unwrap();
            assert!(hand_on(batch.clone()).is_empty());
        });
        assert_eq!(
            raised,
            [
                held(9, "raised here"),
                held(4, "of a task"),
                held(0, "of the next task")
            ]
        );
    }
}
