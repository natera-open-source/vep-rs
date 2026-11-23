// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! An interrupt or a termination signal taken once while a run scores: either sets a flag
//! the run reads between its steps, on which it removes what it wrote and ends, and a
//! repeat of either changes nothing while the run stops. A signal the process ignores
//! stays ignored, and the handlers found are put back when the run leaves. A system call
//! the run waits in when the signal arrives is started again, except under `Unrestarted`.

use std::sync::atomic::{AtomicBool, Ordering};

static REQUESTED: AtomicBool = AtomicBool::new(false);

/// Whether a stop signal arrived since the run began to take them.
pub fn requested() -> bool {
    REQUESTED.load(Ordering::SeqCst)
}

extern "C" fn note(_signal: libc::c_int) {
    REQUESTED.store(true, Ordering::SeqCst);
}

/// The handlers a run found, put back when this is dropped.
pub struct Taken {
    previous: Vec<(libc::c_int, libc::sigaction)>,
}

impl Taken {
    /// Takes the interrupt signal as a stop of the run until the returned value is dropped.
    /// The flag starts cleared.
    pub fn interrupt() -> Taken {
        REQUESTED.store(false, Ordering::SeqCst);
        let mut taken = Taken {
            previous: Vec::new(),
        };
        taken.take(libc::SIGINT);
        taken
    }

    /// Takes the termination signal as well.
    pub fn termination(&mut self) {
        self.take(libc::SIGTERM);
    }

    fn take(&mut self, signal: libc::c_int) {
        // SAFETY: `found` is a zeroed `sigaction` the first call fills; the handler stores
        // to an atomic alone, which a signal handler may do.
        unsafe {
            let mut found: libc::sigaction = std::mem::zeroed();
            if libc::sigaction(signal, std::ptr::null(), &mut found) != 0
                || found.sa_sigaction == libc::SIG_IGN
            {
                return;
            }
            let mut taken: libc::sigaction = std::mem::zeroed();
            taken.sa_sigaction = note as extern "C" fn(libc::c_int) as libc::sighandler_t;
            taken.sa_flags = libc::SA_RESTART;
            libc::sigemptyset(&mut taken.sa_mask);
            if libc::sigaction(signal, &taken, std::ptr::null_mut()) == 0 {
                self.previous.push((signal, found));
            }
        }
    }
}

impl Drop for Taken {
    fn drop(&mut self) {
        for (signal, found) in &self.previous {
            // SAFETY: `found` is the `sigaction` the system handed out for this signal.
            unsafe {
                libc::sigaction(*signal, found, std::ptr::null_mut());
            }
        }
    }
}

/// While this lives, a system call the run waits in fails with `EINTR` when a stop signal
/// the run takes arrives, where the system would otherwise start the call again: the caller
/// reads `requested` and ends its wait. The handlers are as they were once it is dropped,
/// which is before the run's `Taken` is.
pub struct Unrestarted {
    restarted: Vec<(libc::c_int, libc::sigaction)>,
}

impl Unrestarted {
    pub fn new() -> Unrestarted {
        let mut restarted = Vec::new();
        for signal in [libc::SIGINT, libc::SIGTERM] {
            // SAFETY: `found` is a zeroed `sigaction` the first call fills; the second sets
            // the handler found, the run's own, with one flag less.
            unsafe {
                let mut found: libc::sigaction = std::mem::zeroed();
                if libc::sigaction(signal, std::ptr::null(), &mut found) != 0
                    || found.sa_sigaction
                        != note as extern "C" fn(libc::c_int) as libc::sighandler_t
                {
                    continue;
                }
                let mut unrestarted = found;
                unrestarted.sa_flags &= !libc::SA_RESTART;
                if libc::sigaction(signal, &unrestarted, std::ptr::null_mut()) == 0 {
                    restarted.push((signal, found));
                }
            }
        }
        Unrestarted { restarted }
    }
}

impl Default for Unrestarted {
    fn default() -> Unrestarted {
        Unrestarted::new()
    }
}

impl Drop for Unrestarted {
    fn drop(&mut self) {
        for (signal, found) in &self.restarted {
            // SAFETY: `found` is the `sigaction` the system handed out for this signal.
            unsafe {
                libc::sigaction(*signal, found, std::ptr::null_mut());
            }
        }
    }
}
