// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` command under an interrupt and a termination signal: one error line, exit
//! status 130, no scratch and no residual file left, a repeat of either changing nothing,
//! and the handlers found put back.
//!
//! The signals are the process's, so the tests of this file run one at a time.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use vep_plugin_pairs::registry::{BuiltInRegistry, Registry};
use vep_plugin_pairs::rules::Refusal;
use vep_plugin_pairs::score::api::{BuiltIn, Consult, Rules};
use vep_plugin_pairs::score::run;

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
const INTERRUPTED: &str =
    "ERROR: [score_plugin_tuples] interrupted before the report was written\n";
static ONE_AT_A_TIME: Mutex<()> = Mutex::new(());
/// The times the handler a test set for itself ran.
static CALLERS_HANDLER_RAN: AtomicUsize = AtomicUsize::new(0);

extern "C" fn callers(_signal: libc::c_int) {
    CALLERS_HANDLER_RAN.fetch_add(1, Ordering::SeqCst);
}

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-signals-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_string()
    }

    fn output(&self, name: &str, rows: &[String]) -> String {
        fs::write(self.0.join(name), format!("{HEADER}{}", rows.concat())).unwrap();
        self.path(name)
    }

    fn names(&self, name: &str) -> Vec<String> {
        fs::read_dir(self.0.join(name))
            .map(|entries| {
                entries
                    .map(|entry| entry.unwrap().file_name().into_string().unwrap())
                    .collect()
            })
            .unwrap_or_default()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn row(location: &str, extra: &str) -> String {
    format!(
        "v\t{location}\tT\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\t{extra}\n"
    )
}

/// The tool's registry, which sends the process one of `signals` each time it is asked for
/// the plugin of a field, until none is left.
struct Signalling {
    signals: Vec<libc::c_int>,
    asked: AtomicUsize,
}

impl Signalling {
    fn sending(signals: &[libc::c_int]) -> Signalling {
        Signalling {
            signals: signals.to_vec(),
            asked: AtomicUsize::new(0),
        }
    }
}

impl Registry for Signalling {
    fn all_plugin_names(&self) -> &[String] {
        BuiltInRegistry.all_plugin_names()
    }

    fn plugin_alias(&self, name: &str) -> Option<&str> {
        BuiltInRegistry.plugin_alias(name)
    }

    fn plugin_of_field(
        &self,
        field: &str,
        plugins: &[String],
        flags: Option<&[String]>,
    ) -> Result<Option<String>, String> {
        if let Some(&signal) = self.signals.get(self.asked.fetch_add(1, Ordering::SeqCst)) {
            // SAFETY: `raise` sends this thread a signal whose handler the run has set.
            unsafe {
                libc::raise(signal);
            }
        }
        BuiltInRegistry.plugin_of_field(field, plugins, flags)
    }
}

/// The handler the process has for `signal`.
fn handler(signal: libc::c_int) -> libc::sighandler_t {
    // SAFETY: `found` is a zeroed `sigaction` the call fills.
    unsafe {
        let mut found: libc::sigaction = std::mem::zeroed();
        libc::sigaction(signal, std::ptr::null(), &mut found);
        found.sa_sigaction
    }
}

fn set_handler(signal: libc::c_int, handler: libc::sighandler_t) {
    // SAFETY: `handler` is `SIG_DFL`, `SIG_IGN` or a function that stores to an atomic.
    unsafe {
        libc::signal(signal, handler);
    }
}

/// Runs a command line over a reference and two engines with `registry` and returns the
/// exit status and the two streams' text.
fn scored(scratch: &Scratch, registry: &dyn Registry) -> (i32, String, String) {
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "LOEUF=0.31;IMPACT=LOW"),
            row("1:200", "LOEUF=0.5"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "LOEUF=0.31;IMPACT=LOW"),
            row("1:200", "LOEUF=0.7"),
        ],
    );
    fs::write(
        scratch.0.join("classes.tsv"),
        "#class\tplugin\tfield names\n",
    )
    .unwrap();
    let arguments: Vec<String> = [
        "--reference",
        &reference,
        "--engine",
        &format!("first={engine}"),
        "--engine",
        &format!("second={engine}"),
        "--plugins",
        "LOEUF",
        "--class-table",
        &scratch.path("classes.tsv"),
        "--work-dir",
        &scratch.path("work"),
        "--report",
        &scratch.path("report.json"),
        "--residuals",
        &scratch.path("residuals"),
    ]
    .iter()
    .map(|word| word.to_string())
    .collect();
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = run(&arguments, registry, &BuiltIn, &mut out, &mut error);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(error).unwrap(),
    )
}

#[test]
fn termination_ends_in_one_error_line_and_leaves_nothing_behind() {
    let _alone = ONE_AT_A_TIME
        .lock()
        .unwrap_or_else(|held| held.into_inner());
    let scratch = Scratch::new("termination");
    let before = handler(libc::SIGTERM);
    set_handler(
        libc::SIGTERM,
        callers as extern "C" fn(libc::c_int) as libc::sighandler_t,
    );
    CALLERS_HANDLER_RAN.store(0, Ordering::SeqCst);
    let done = scored(&scratch, &Signalling::sending(&[libc::SIGTERM]));
    let restored = handler(libc::SIGTERM);
    set_handler(libc::SIGTERM, before);
    assert_eq!(done, (130, String::new(), INTERRUPTED.to_string()));
    assert!(!scratch.0.join("report.json").exists());
    assert!(scratch.names("work").is_empty());
    assert!(scratch.names("residuals").is_empty());
    assert_eq!(
        restored,
        callers as extern "C" fn(libc::c_int) as libc::sighandler_t
    );
    assert_eq!(CALLERS_HANDLER_RAN.load(Ordering::SeqCst), 0);
}

#[test]
fn stop_signal_ends_the_run_and_a_repeat_changes_nothing() {
    let _alone = ONE_AT_A_TIME
        .lock()
        .unwrap_or_else(|held| held.into_inner());
    for (what, signals) in [
        ("a termination", vec![libc::SIGTERM]),
        ("an interrupt", vec![libc::SIGINT]),
        ("a second termination", vec![libc::SIGTERM, libc::SIGTERM]),
        (
            "a termination after an interrupt",
            vec![libc::SIGINT, libc::SIGTERM],
        ),
        (
            "an interrupt after a termination",
            vec![libc::SIGTERM, libc::SIGINT],
        ),
        ("a second interrupt", vec![libc::SIGINT, libc::SIGINT]),
    ] {
        let scratch = Scratch::new("repeat");
        let before = (handler(libc::SIGINT), handler(libc::SIGTERM));
        let done = scored(&scratch, &Signalling::sending(&signals));
        assert_eq!(
            done,
            (130, String::new(), INTERRUPTED.to_string()),
            "{what}"
        );
        assert!(scratch.names("work").is_empty(), "{what}");
        assert!(scratch.names("residuals").is_empty(), "{what}");
        assert!(!scratch.0.join("report.json").exists(), "{what}");
        assert_eq!(
            (handler(libc::SIGINT), handler(libc::SIGTERM)),
            before,
            "{what}"
        );
    }
}

#[test]
fn signal_the_process_ignores_stays_ignored() {
    let _alone = ONE_AT_A_TIME
        .lock()
        .unwrap_or_else(|held| held.into_inner());
    let scratch = Scratch::new("ignored");
    let before = handler(libc::SIGINT);
    set_handler(libc::SIGINT, libc::SIG_IGN);
    let done = scored(&scratch, &Signalling::sending(&[libc::SIGINT]));
    let left = handler(libc::SIGINT);
    set_handler(libc::SIGINT, before);
    assert_eq!((done.0, done.2.as_str()), (0, ""));
    assert_eq!(left, libc::SIG_IGN);
    assert!(scratch.0.join("report.json").exists());
}

/// Class rules of one class, `planted`, that accepts nothing and sends the process a
/// termination signal the `at`-th time it is consulted.
struct Terminating {
    at: usize,
    consulted: AtomicUsize,
}

impl Rules for Terminating {
    fn knows(&self, klass: &str) -> bool {
        klass == "planted"
    }

    fn accepts(&self, _consult: &Consult) -> Result<bool, Refusal> {
        if self.consulted.fetch_add(1, Ordering::SeqCst) + 1 == self.at {
            // SAFETY: `raise` sends this thread a signal whose handler the run has set.
            unsafe {
                libc::raise(libc::SIGTERM);
            }
        }
        Ok(false)
    }

    fn transcript_reading(&self, _klass: &str) -> bool {
        false
    }

    fn sibling_allele(&self, _klass: &str) -> bool {
        false
    }

    fn record_reference(&self, _klass: &str) -> bool {
        false
    }
}

#[test]
fn termination_while_the_second_engine_is_settled_leaves_no_residual_file_of_the_first() {
    let _alone = ONE_AT_A_TIME
        .lock()
        .unwrap_or_else(|held| held.into_inner());
    let scratch = Scratch::new("second-engine");
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "LOEUF=0.31"), row("1:200", "LOEUF=0.5")],
    );
    let engine = scratch.output(
        "engine.txt",
        &[row("1:100", "LOEUF=0.31"), row("1:200", "LOEUF=0.7")],
    );
    fs::write(scratch.0.join("classes.tsv"), "planted\tLOEUF\tLOEUF\n").unwrap();
    let line = |more: &[&str]| -> Vec<String> {
        let mut words = vec![
            "--reference".to_string(),
            reference.clone(),
            "--engine".to_string(),
            format!("first={engine}"),
            "--engine".to_string(),
            format!("second={engine}"),
            "--plugins".to_string(),
            "LOEUF".to_string(),
            "--class-table".to_string(),
            scratch.path("classes.tsv"),
            "--work-dir".to_string(),
            scratch.path("work"),
            "--report".to_string(),
            scratch.path("report.json"),
        ];
        words.extend(more.iter().map(|word| word.to_string()));
        words
    };
    // The class is consulted as often for each engine: once for the one changed value.
    let counting = Terminating {
        at: usize::MAX,
        consulted: AtomicUsize::new(0),
    };
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = run(
        &line(&["--residuals", &scratch.path("counted")]),
        &BuiltInRegistry,
        &counting,
        &mut out,
        &mut error,
    );
    assert_eq!((code, String::from_utf8_lossy(&error).as_ref()), (0, ""));
    let consulted = counting.consulted.load(Ordering::SeqCst);
    assert!(consulted >= 2 && consulted.is_multiple_of(2), "{consulted}");
    assert_eq!(scratch.names("counted").len(), 4);
    fs::remove_file(scratch.0.join("report.json")).unwrap();
    let terminating = Terminating {
        at: consulted / 2 + 1,
        consulted: AtomicUsize::new(0),
    };
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = run(
        &line(&["--residuals", &scratch.path("residuals")]),
        &BuiltInRegistry,
        &terminating,
        &mut out,
        &mut error,
    );
    assert_eq!(
        (code, out.len(), String::from_utf8_lossy(&error).as_ref()),
        (130, 0, INTERRUPTED)
    );
    assert!(terminating.consulted.load(Ordering::SeqCst) > consulted / 2);
    assert!(!scratch.0.join("report.json").exists());
    assert!(scratch.names("work").is_empty());
    assert!(scratch.names("residuals").is_empty());
}
