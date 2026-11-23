// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The CPUs a thread may run on, as the kernel's scheduler affinity holds
//! them: read for the run's own set, and set when the reader thread is given a
//! core of its own ([`crate::config::reserve_reader_core`]): the reader onto
//! its CPU, every other thread of the process onto the rest. Linux carries the
//! calls; on any other system the set reads as unknown and nothing moves.

use crate::config::ReaderCoreReservation;

/// What [`apply`] did: whether the reader was pinned, the other threads it
/// moved and the ones it could not (gone by the time they were named, or
/// refused by the kernel).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct Applied {
    pub reader_pinned: bool,
    pub moved: usize,
    pub failed: usize,
}

/// The logical CPUs the calling thread may run on, `None` where the system
/// does not report them.
pub fn allowed_cpus() -> Option<Vec<usize>> {
    platform::allowed_cpus()
}

/// Puts every other thread of the process on the reservation's other CPUs,
/// then the calling thread on its reader CPU; when no other thread could be
/// moved the caller stays where it is, since a CPU of its own is a core of its
/// own only once the rest of the run has left it. A thread created afterwards
/// inherits its creator's set, so a thread the reader creates runs on the
/// reader's CPU and one any other thread creates on the rest.
pub(crate) fn apply(reservation: &ReaderCoreReservation) -> Applied {
    platform::apply(reservation)
}

#[cfg(target_os = "linux")]
mod platform {
    use std::io;
    use std::mem;

    use super::Applied;
    use crate::config::ReaderCoreReservation;

    /// The CPUs a `cpu_set_t` can name.
    const CPU_SET_CAPACITY: usize = libc::CPU_SETSIZE as usize;

    fn cpu_set(cpus: &[usize]) -> Option<libc::cpu_set_t> {
        // SAFETY: a zeroed cpu_set_t is the empty set, and CPU_SET writes a bit
        // of it for a CPU under its capacity.
        let mut set: libc::cpu_set_t = unsafe { mem::zeroed() };
        let mut any = false;
        for &cpu in cpus {
            if cpu < CPU_SET_CAPACITY {
                unsafe { libc::CPU_SET(cpu, &mut set) };
                any = true;
            }
        }
        any.then_some(set)
    }

    pub(super) fn allowed_cpus() -> Option<Vec<usize>> {
        // SAFETY: a zeroed cpu_set_t is a valid set for the kernel to fill, and
        // the size passed is the set's own.
        let mut set: libc::cpu_set_t = unsafe { mem::zeroed() };
        let status =
            unsafe { libc::sched_getaffinity(0, mem::size_of::<libc::cpu_set_t>(), &mut set) };
        if status != 0 {
            return None;
        }
        // SAFETY: CPU_ISSET reads a bit of a set the kernel filled, for a CPU
        // under the set's capacity.
        let cpus: Vec<usize> = (0..CPU_SET_CAPACITY)
            .filter(|&cpu| unsafe { libc::CPU_ISSET(cpu, &set) })
            .collect();
        (!cpus.is_empty()).then_some(cpus)
    }

    /// Restricts the thread `tid` of this process to `cpus`.
    pub(super) fn set_thread_cpus(tid: libc::pid_t, cpus: &[usize]) -> io::Result<()> {
        let set = cpu_set(cpus).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "no CPU the affinity set can name",
            )
        })?;
        // SAFETY: the set is initialised and the size passed is its own.
        let status =
            unsafe { libc::sched_setaffinity(tid, mem::size_of::<libc::cpu_set_t>(), &set) };
        if status == 0 {
            Ok(())
        } else {
            Err(io::Error::last_os_error())
        }
    }

    pub(super) fn current_thread_id() -> libc::pid_t {
        // SAFETY: gettid takes no argument and cannot fail.
        unsafe { libc::gettid() }
    }

    /// Every thread of this process, by kernel thread id.
    pub(super) fn thread_ids() -> Vec<libc::pid_t> {
        std::fs::read_dir("/proc/self/task")
            .map(|entries| {
                entries
                    .flatten()
                    .filter_map(|entry| entry.file_name().to_str()?.parse().ok())
                    .collect()
            })
            .unwrap_or_default()
    }

    pub(super) fn apply(reservation: &ReaderCoreReservation) -> Applied {
        let mut applied = Applied::default();
        let own_tid = current_thread_id();
        // The other threads first: a thread that is being moved onto the rest
        // never shares the reader's CPU with the reader once it is pinned.
        for tid in thread_ids() {
            if tid == own_tid {
                continue;
            }
            match set_thread_cpus(tid, &reservation.others) {
                Ok(()) => applied.moved += 1,
                Err(_) => applied.failed += 1,
            }
        }
        if applied.moved > 0 {
            applied.reader_pinned = set_thread_cpus(own_tid, &[reservation.reader_cpu]).is_ok();
        }
        applied
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::sync::mpsc;

        /// The `Cpus_allowed_list` line of a thread's status file, as a CPU list.
        fn status_cpus(tid: libc::pid_t) -> Vec<usize> {
            let status = std::fs::read_to_string(format!("/proc/self/task/{tid}/status")).unwrap();
            let line = status
                .lines()
                .find_map(|line| line.strip_prefix("Cpus_allowed_list:"))
                .expect("a Cpus_allowed_list line");
            line.trim()
                .split(',')
                .flat_map(|token| match token.split_once('-') {
                    Some((a, b)) => (a.parse().unwrap()..=b.parse().unwrap()).collect(),
                    None => vec![token.parse().unwrap()],
                })
                .collect()
        }

        /// Another thread of the process, moved onto one of the allowed CPUs
        /// and back, reads each set from its status file, and an empty set is
        /// refused; the calling thread reads its own id and its allowed set,
        /// which its status file agrees with.
        #[test]
        fn a_thread_is_moved_onto_a_cpu_and_back() {
            let allowed = allowed_cpus().expect("Linux reports the affinity set");
            assert!(!allowed.is_empty());
            let own_tid = current_thread_id();
            assert!(own_tid > 0);
            assert_eq!(status_cpus(own_tid), allowed);
            assert!(thread_ids().contains(&own_tid));

            let (tid_tx, tid_rx) = mpsc::channel();
            let (done_tx, done_rx) = mpsc::channel::<()>();
            let parked = std::thread::spawn(move || {
                tid_tx.send(current_thread_id()).unwrap();
                let _ = done_rx.recv();
            });
            let other = tid_rx.recv().unwrap();
            assert_ne!(other, own_tid);
            let one = *allowed.last().unwrap();
            set_thread_cpus(other, &[one]).unwrap();
            assert_eq!(status_cpus(other), vec![one]);
            set_thread_cpus(other, &allowed).unwrap();
            assert_eq!(status_cpus(other), allowed);
            assert!(set_thread_cpus(other, &[]).is_err());
            done_tx.send(()).unwrap();
            parked.join().unwrap();
        }
    }
}

#[cfg(not(target_os = "linux"))]
mod platform {
    use super::Applied;
    use crate::config::ReaderCoreReservation;

    pub(super) fn allowed_cpus() -> Option<Vec<usize>> {
        None
    }

    pub(super) fn apply(_reservation: &ReaderCoreReservation) -> Applied {
        Applied::default()
    }
}
