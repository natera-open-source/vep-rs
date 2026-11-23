// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The warning sink and the skip registry: Ensembl VEP's `warning_msg`,
//! `skipped_variant_msg` and `warning_fh` (BaseVEP.pm:802-901) and its
//! `--skipped_variants_file` (BaseRunner.pm:331-380).
//!
//! Every warning of a run goes through one [`Sink`]: the `warn!` events of every
//! crate, routed here by [`WarningLayer`], and the skipped-variant events the
//! pipeline raises. An event raised while a record is annotated is held for
//! the record (`holding_events`) and reaches the sink with the record's
//! other warnings, as does one a plugin or a `--custom` source raises while a
//! whole batch is annotated (`vep_core::warning_hold`). A warning is written
//! as `WARNING: <text>` to the warnings
//! file (`--warning_file`, default `<output_file>_warnings.txt`, the literal
//! `STDERR` in any case for standard error) and echoed on standard error unless
//! `--quiet`, so a standard-error target carries it twice; a text seen before, a
//! trailing `line N` set aside, is written once. The warnings file is opened at
//! the first warning, so a run without one leaves no file. A process that never
//! configures the sink, a program embedding the library, gets each warning once
//! on standard error: the echo is a configured run's. A run that carries
//! its warnings in its own output (`vep recode`) targets the sink's memory and
//! reads them back when it ends. A built-in plugin's failures are held through
//! the run and written at its end as one `Plugin '<name>' went wrong` line per
//! plugin carrying the number of variants it failed on, as `warning_msg`'s seen
//! set leaves one line for a plugin dying the same way on every line
//! (BaseVEP.pm:810-813).
//!
//! A skipped variant is reported as `line N skipped (<line>): <description>`,
//! the line collapsed to single spaces and cut to 50 bytes, `variant skipped`
//! under `--input_data`, and recorded for the skipped-variants file, written at
//! the end of the run as a header line and one `Line N\t<line>\t<description>`
//! row per skip.

use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{LazyLock, Mutex, MutexGuard};

use tracing::field::{Field, Visit};
use tracing::{Event, Level, Subscriber};
use tracing_subscriber::layer::{Context, Layer};
use vep_core::warning_hold::RECORD_FIELD;

use crate::config::CliError;

/// The line shown in a skipped-variant warning is cut past this many bytes:
/// `$maxChar = 50` (BaseVEP.pm:848) measured by `length` on a line Ensembl VEP reads
/// with no encoding layer, so it counts bytes.
const QUOTED_LINE_BYTES: usize = 50;

/// The width of the line column of the skipped-variants file (`%-40.40s`).
const SKIPPED_LINE_WIDTH: usize = 40;

/// Where warnings are written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Target {
    /// Standard error alone: `--warning_file STDERR`, which carries a warning
    /// twice unless `--quiet`, or no run configured yet, which carries it once.
    Stderr,
    /// A file, opened at the first warning.
    File(PathBuf),
    /// Held in the sink for the run to read back ([`Sink::take_collected`]):
    /// `vep recode` carries its warnings in its output.
    Memory,
}

/// One skipped variant, as the skipped-variants file lists it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Skipped {
    /// The input line as quoted in the warning, when the event had one.
    pub line: Option<String>,
    pub line_number: Option<u64>,
    pub description: String,
}

struct State {
    target: Target,
    /// Open once the first warning is written to a file target.
    file: Option<File>,
    /// Why the warnings file could not be opened, kept for the run to fail on.
    failure: Option<String>,
    /// The run tests a file target before its first record, as a file that
    /// `--warning_file` names is tested; the file under the default name is
    /// first met by the first warning.
    file_tested: bool,
    /// A warning is written to standard error as well as to the target: a
    /// configured run without `--quiet`, never a sink no run configured.
    echo: bool,
    /// `--input_data` is given: a skip names no line number.
    input_data: bool,
    seen: HashSet<String>,
    /// The skips to list in the skipped-variants file, kept only while one is
    /// reserved: a run without `--skipped_variants_file` holds no skip.
    skipped: Vec<Skipped>,
    /// The `--skipped_variants_file` path, checked before the run and written at
    /// its end.
    skipped_file: Option<PathBuf>,
    /// The texts written under [`Target::Memory`], in order, each with the
    /// input line number of the record it skipped, when it skipped one.
    collected: Vec<(Option<u64>, String)>,
    /// Per failed built-in plugin, in the order the failures were first
    /// reported: the plugin's name and load-order instance, the text of its
    /// first failure and the number of variants it failed on over the run,
    /// written at the run's end.
    plugin_failures: Vec<PluginFailureCount>,
}

/// One built-in plugin instance's failures over a run
/// ([`Sink::plugin_failure`]).
struct PluginFailureCount {
    plugin: String,
    instance: usize,
    reason: String,
    failed_variants: usize,
}

impl Default for State {
    /// Standard error with no echo, so a process that never configures the
    /// sink reads each warning once; nothing seen or skipped.
    fn default() -> Self {
        State {
            target: Target::Stderr,
            file: None,
            failure: None,
            file_tested: false,
            echo: false,
            input_data: false,
            seen: HashSet::new(),
            skipped: Vec::new(),
            skipped_file: None,
            collected: Vec::new(),
            plugin_failures: Vec::new(),
        }
    }
}

/// The warning sink of one run.
pub struct Sink {
    state: Mutex<State>,
}

impl Default for Sink {
    fn default() -> Self {
        Sink::new()
    }
}

impl Sink {
    /// A sink no run has configured: each warning goes to standard error once.
    pub fn new() -> Self {
        Sink {
            state: Mutex::new(State::default()),
        }
    }

    fn lock(&self) -> MutexGuard<'_, State> {
        self.state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    /// Puts the sink back as a process starts it: standard error with no echo,
    /// nothing seen, skipped, collected or held, no file open.
    #[cfg(test)]
    fn forget_every_run(&self) {
        *self.lock() = State::default();
    }

    /// Points the sink at a run's target before the run starts; the state of an
    /// earlier run (its seen set and skips) is dropped. Unless `quiet`, the
    /// configured sink echoes every warning on standard error.
    pub fn configure(&self, target: Target, quiet: bool, input_data: bool) {
        *self.lock() = State {
            target,
            echo: !quiet,
            input_data,
            ..State::default()
        };
    }

    /// `warning_msg` (BaseVEP.pm:802-822): `WARNING: ` is prefixed unless the text
    /// already opens with `warn` in any case, a newline is appended unless present,
    /// a text seen before with the same dedup key is dropped, and the text goes to
    /// the target and to standard error unless quiet, so a standard-error target
    /// carries it twice, as `print $fh` (:819) then `warn` (:821) write it. A
    /// sink no run configured writes it to standard error once.
    pub fn warn(&self, text: &str) {
        self.warn_about_line(text, None);
    }

    /// [`Sink::warn`] with the input line number of the record the warning is
    /// about, kept beside the text under [`Target::Memory`].
    fn warn_about_line(&self, text: &str, line_number: Option<u64>) {
        Self::write_to(&mut self.lock(), text, line_number);
    }

    /// Writes one warning text, finished ([`finished_text`]) and not seen before
    /// (BaseVEP.pm:810-813), to the target, then to standard error in a
    /// configured run that is not quiet, a standard-error target included; a
    /// text seen before is dropped, and a sink no run configured echoes
    /// none. A file that cannot be created or written is kept as the run's failure
    /// ([`Sink::failure`], the `ERROR: [output]` line the run ends with); the
    /// standard-error echo follows `--quiet` all the same. In the `vep`
    /// binary a text standard error refuses, as the target or as the echo,
    /// is dropped ([`crate::stderr`]).
    fn write_to(state: &mut State, text: &str, line_number: Option<u64>) {
        let text = finished_text(text);
        if !state.seen.insert(dedup_key(&text)) {
            return;
        }
        match state.target.clone() {
            Target::Stderr => crate::stderr::write(format_args!("{text}")),
            Target::Memory => state.collected.push((line_number, text.clone())),
            Target::File(path) => {
                let could_not_write =
                    || format!("Could not write to warnings file {}", path.display());
                if state.file.is_none() && state.failure.is_none() {
                    match File::create(&path) {
                        Ok(file) => state.file = Some(file),
                        Err(_) => state.failure = Some(could_not_write()),
                    }
                }
                if let Some(file) = state.file.as_mut() {
                    if file.write_all(text.as_bytes()).is_err() {
                        state.failure = Some(could_not_write());
                    }
                }
            }
        }
        if state.echo {
            crate::stderr::write(format_args!("{text}"));
        }
    }

    /// `skipped_variant_msg` (BaseVEP.pm:840-863): warns about a record the run
    /// drops and records it for the skipped-variants file when one is reserved.
    /// `line` is the input line, `line_number` Ensembl VEP's count of data records
    /// read. The warning passes through the seen set like any other (`line N
    /// skipped` opens the text, so the trailing-number rule never strips the
    /// number): two records are two lines, and two copies of one record, as
    /// `--individual` makes them, are one.
    pub fn skipped_variant(&self, description: &str, line_number: Option<u64>, line: Option<&str>) {
        let line = line.map(quoted_line);
        let described = match &line {
            Some(line) => format!("({line}): {description}"),
            None => description.to_string(),
        };
        let numbered = line_number.filter(|_| !self.lock().input_data);
        let text = match numbered {
            Some(n) => format!("WARNING: line {n} skipped {described}\n"),
            None => format!("WARNING: variant skipped {described}\n"),
        };
        self.warn_about_line(&text, line_number);
        let mut state = self.lock();
        if state.skipped_file.is_some() {
            state.skipped.push(Skipped {
                line,
                line_number,
                description: description.to_string(),
            });
        }
    }

    /// The file the warnings go to, created at the first warning; `None`
    /// under any other target.
    #[cfg(feature = "parquet")]
    pub(crate) fn file(&self) -> Option<PathBuf> {
        match &self.lock().target {
            Target::File(path) => Some(path.clone()),
            Target::Stderr | Target::Memory => None,
        }
    }

    /// Has the run test the warnings file before it reads its first record
    /// (`file_to_test`), as it tests the other files it names: asked, after
    /// [`Sink::configure`], which withdraws the request, for a file that
    /// `--warning_file` names. The file under the default name is left to
    /// the first warning, which opens either file in Ensembl VEP
    /// (BaseVEP.pm:817 calls `warning_fh`, 878-901), so a run that draws no
    /// warning never meets it.
    pub fn test_file_before_first_record(&self) {
        self.lock().file_tested = true;
    }

    /// The warnings file the run tests before its first record: the file of
    /// a sink asked to have it tested ([`Sink::test_file_before_first_record`])
    /// while no warning has opened it. `None` under any other target, for a
    /// file no caller asked the test for, once a warning opened the file, and
    /// once one could not.
    pub(crate) fn file_to_test(&self) -> Option<PathBuf> {
        let state = self.lock();
        match &state.target {
            Target::File(path)
                if state.file_tested && state.file.is_none() && state.failure.is_none() =>
            {
                Some(path.clone())
            }
            _ => None,
        }
    }

    /// The skipped-variants file reserved for this run, written at its end.
    pub(crate) fn skipped_variants_file(&self) -> Option<PathBuf> {
        self.lock().skipped_file.clone()
    }

    /// The warnings-file failure of this run, if the file could not be written:
    /// the message of the `ERROR: [output]` line the run exits with.
    pub fn failure(&self) -> Option<CliError> {
        self.lock()
            .failure
            .clone()
            .map(|message| CliError::new("output", message))
    }

    /// Reserves the `--skipped_variants_file` before the run: refused with Ensembl
    /// VEP's text when the file exists and `--force_overwrite` is not given
    /// (BaseRunner.pm:337-340, a check Ensembl VEP makes at the end of the run).
    /// The file itself is created at [`Sink::finish`], as Ensembl VEP creates it
    /// (`:364-380`), so a run that fails leaves none behind.
    pub fn reserve_skipped_variants_file(
        &self,
        path: &str,
        force_overwrite: bool,
    ) -> Result<(), CliError> {
        let path = PathBuf::from(path);
        if path.exists() && !force_overwrite {
            return Err(CliError::new(
                "output",
                format!(
                    "Skipped variants file {} already exists. Specify a different output file using --skipped_variants_file or overwrite existing file with --force_overwrite",
                    path.display()
                ),
            ));
        }
        self.lock().skipped_file = Some(path);
        Ok(())
    }

    /// The skips recorded so far, in input order.
    pub fn skipped(&self) -> Vec<Skipped> {
        self.lock().skipped.clone()
    }

    /// The texts written under [`Target::Memory`] since the last call, each a
    /// finished warning line with the input line number of the record it
    /// skipped, when it skipped one; the sink keeps none of them.
    pub fn take_collected(&self) -> Vec<(Option<u64>, String)> {
        std::mem::take(&mut self.lock().collected)
    }

    /// A built-in plugin's failure on a batch: `instance` its place in load
    /// order (a built-in loaded twice is two plugins of one name, each with its
    /// own line, as Ensembl VEP loads one instance per `--plugin` argument),
    /// `reason` its error text and `failed_variants` the batch's variants it
    /// annotated nothing on (OutputFactory.pm:2252-2254: `Plugin '<Class>' went
    /// wrong: $@`, one `warning_msg` per failing line, which the seen set leaves
    /// as one line for a plugin dying the same way on every line). The
    /// failures of an instance are held and written at the run's end
    /// ([`Sink::finish`]) as one line, `Plugin '<name>' went wrong: <first
    /// reason>; <N> variants affected`, summing the batches.
    pub fn plugin_failure(
        &self,
        plugin: &str,
        instance: usize,
        reason: &str,
        failed_variants: usize,
    ) {
        let mut state = self.lock();
        match state
            .plugin_failures
            .iter_mut()
            .find(|failure| failure.plugin == plugin && failure.instance == instance)
        {
            Some(failure) => failure.failed_variants += failed_variants,
            None => state.plugin_failures.push(PluginFailureCount {
                plugin: plugin.to_string(),
                instance,
                reason: reason.to_string(),
                failed_variants,
            }),
        }
    }

    /// Writes the plugin failures held so far ([`Sink::plugin_failure`]), one
    /// line per plugin, through the seen set like any other warning, taking and
    /// writing them under one lock so two callers never split a failure between
    /// them; [`Sink::finish`] runs it, a run ending in an error runs it before it
    /// exits so the warnings file keeps the lines, and a program's batch runs it
    /// as the batch ends.
    pub fn write_plugin_failures(&self) {
        let mut state = self.lock();
        let failures = std::mem::take(&mut state.plugin_failures);
        for PluginFailureCount {
            plugin,
            reason,
            failed_variants,
            ..
        } in failures
        {
            let variants = if failed_variants == 1 {
                "variant"
            } else {
                "variants"
            };
            let text = format!(
                "Plugin '{plugin}' went wrong: {reason}; {failed_variants} {variants} affected"
            );
            Self::write_to(&mut state, &text, None);
        }
    }

    /// Ends the run: writes the plugin failures held through it
    /// ([`Sink::plugin_failure`]) and the skipped-variants file
    /// (BaseRunner.pm:364-380, the header line then one row per skip, the
    /// header alone when nothing was skipped) and closes the warnings file. The
    /// warnings-file failure, if any, is returned after the skipped-variants
    /// file is written.
    pub fn finish(&self, input_file: &str) -> Result<(), CliError> {
        self.write_plugin_failures();
        let mut state = self.lock();
        if let Some(path) = state.skipped_file.take() {
            let mut text = format!("[VEP skipped the following variants from {input_file}]\n");
            for skip in &state.skipped {
                text.push_str(&skipped_row(skip));
            }
            std::fs::write(&path, text).map_err(|_| skipped_file_unwritable(&path))?;
        }
        state.file = None;
        state.seen.clear();
        state.skipped.clear();
        match state.failure.take() {
            Some(message) => Err(CliError::new("output", message)),
            None => Ok(()),
        }
    }
}

/// The `output` component's error for a warnings file that cannot be created
/// or written, with Ensembl VEP's text (BaseVEP.pm:894): the refusal of a run
/// that tests the file before a warning opens it.
pub(crate) fn warnings_file_unwritable(path: &Path) -> CliError {
    CliError::new(
        "output",
        format!("Could not write to warnings file {}", path.display()),
    )
}

/// The `output` component's error for a skipped-variants file that cannot be
/// created or written, with Ensembl VEP's text.
pub(crate) fn skipped_file_unwritable(path: &Path) -> CliError {
    CliError::new(
        "output",
        format!(
            "Could not write to skipped variants file {}",
            path.display()
        ),
    )
}

/// A warning's finished text: `WARNING: ` prefixed unless the text already
/// opens with `warn` in any case, a newline appended unless present
/// (BaseVEP.pm:806-807).
fn finished_text(text: &str) -> String {
    let mut text = text.to_string();
    if !text
        .get(..4)
        .is_some_and(|head| head.eq_ignore_ascii_case("warn"))
    {
        text.insert_str(0, "WARNING: ");
    }
    if !text.ends_with('\n') {
        text.push('\n');
    }
    text
}

/// The text a warning is deduplicated on: the text without its newline and
/// without a trailing `line N` (`s/line \d+\n?$//`).
fn dedup_key(text: &str) -> String {
    let text = text.trim_end_matches('\n');
    if let Some(at) = text.rfind("line ") {
        let digits = &text[at + "line ".len()..];
        if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
            return text[..at].to_string();
        }
    }
    text.to_string()
}

/// The shortest prefix of a raw input line that [`quoted_line`] renders as it
/// renders the whole line: the line through its first `QUOTED_LINE_BYTES + 1`
/// bytes once whitespace runs are collapsed, cut at a character boundary. A
/// reader that keeps part of a line keeps this much for the skip warning.
pub(crate) fn quotable_prefix(line: &[u8]) -> &[u8] {
    let mut bytes = 0usize;
    let mut in_space = false;
    for (at, &byte) in line.iter().enumerate() {
        // A continuation byte stays with the character its lead byte opened.
        if byte & 0xC0 == 0x80 {
            bytes += 1;
            continue;
        }
        if bytes > QUOTED_LINE_BYTES {
            return &line[..at];
        }
        if byte.is_ascii_whitespace() || byte == 0x0B {
            if !in_space {
                bytes += 1;
            }
            in_space = true;
        } else {
            bytes += 1;
            in_space = false;
        }
    }
    line
}

/// The input line as a skipped-variant warning quotes it: every run of
/// whitespace (Perl's `\s`, the vertical tab included) one space, and a line
/// over [`QUOTED_LINE_BYTES`] bytes cut at the last character boundary at or
/// before byte `QUOTED_LINE_BYTES - 4`, then `...` (BaseVEP.pm:848-850 counts
/// bytes and may split a character there; a whole character is kept here).
pub(crate) fn quoted_line(line: &str) -> String {
    let mut collapsed = String::with_capacity(line.len());
    let mut in_space = false;
    for c in line.chars() {
        if c.is_ascii_whitespace() || c == '\u{0B}' {
            if !in_space {
                collapsed.push(' ');
            }
            in_space = true;
        } else {
            collapsed.push(c);
            in_space = false;
        }
    }
    if collapsed.len() > QUOTED_LINE_BYTES {
        let mut end = QUOTED_LINE_BYTES - 4;
        while !collapsed.is_char_boundary(end) {
            end -= 1;
        }
        collapsed.truncate(end);
        collapsed.push_str("...");
    }
    collapsed
}

/// One row of the skipped-variants file: `Line %-5s\t%-40.40s\t%s`.
fn skipped_row(skip: &Skipped) -> String {
    let number = skip.line_number.map(|n| n.to_string()).unwrap_or_default();
    let line: String = skip
        .line
        .as_deref()
        .unwrap_or("")
        .chars()
        .take(SKIPPED_LINE_WIDTH)
        .collect();
    format!(
        "Line {number:<5}\t{line:<width$}\t{}\n",
        skip.description,
        width = SKIPPED_LINE_WIDTH
    )
}

static GLOBAL: LazyLock<Sink> = LazyLock::new(Sink::new);

/// The process's sink: configured by the binary for its run, written by every
/// warning event and skip, and until configured standard error alone, one
/// copy of each warning.
pub fn global() -> &'static Sink {
    // A unit test that reaches the sink with no guard can write into the
    // warnings file of the test that holds one.
    #[cfg(test)]
    assert!(
        under_sink_test_guard(),
        "a unit test that reaches the process-wide sink holds `sink_test_guard`"
    );
    &GLOBAL
}

/// The lock [`sink_test_guard`] hands out.
#[cfg(test)]
static SINK_TESTS: Mutex<()> = Mutex::new(());

/// The thread of the test that took [`sink_test_guard`] last.
#[cfg(test)]
static SINK_TEST_THREAD: Mutex<Option<std::thread::ThreadId>> = Mutex::new(None);

/// The guard a test holds for its run when it, or the code it runs, reaches
/// the process-wide sink: to write a warning or a skip, to write the plugin
/// failures the sink holds, to read what the sink holds or to reserve its
/// skipped-variants file. Such tests run in turn, and [`global`] panics in a
/// test that reaches the sink without the guard. The holder finds the sink as
/// a process starts it, whatever the holder before it left there, so no test
/// puts the sink back. A test that only reads what one record reports gives
/// it a [`Sink`] of its own instead, which no other test can reach.
#[cfg(test)]
pub(crate) fn sink_test_guard() -> MutexGuard<'static, ()> {
    let guard = SINK_TESTS
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    *SINK_TEST_THREAD
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(std::thread::current().id());
    GLOBAL.forget_every_run();
    guard
}

/// Whether the caller reaches the process-wide sink under [`sink_test_guard`]:
/// a test holds the guard, and a caller on a test's own thread, which the
/// test harness names `<module>::tests::<test>`, is that test. A thread a test
/// starts carries no such name and is taken for the holder's.
#[cfg(test)]
fn under_sink_test_guard() -> bool {
    if !matches!(
        SINK_TESTS.try_lock(),
        Err(std::sync::TryLockError::WouldBlock)
    ) {
        return false;
    }
    let caller = std::thread::current();
    let on_a_test_thread = caller.name().is_some_and(|name| name.contains("::tests::"));
    let holder = *SINK_TEST_THREAD
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    !on_a_test_thread || holder == Some(caller.id())
}

/// Runs `work` with every `WARN` event this thread raises kept from the sink,
/// and returns those events, in the order raised, beside what `work` returned
/// (`vep_core::warning_hold::holding`, the record each event names set
/// aside). The caller reports them where the work stands among the run's
/// warnings: a record's with the record's other warnings, so a warning the
/// annotation of a record raises is listed in input order and never for a
/// batch the run did not reach. An event another thread raises is not held; a
/// call inside another holds its own events and leaves the outer call's as
/// they were; the events of work that unwinds are dropped with it.
pub(crate) fn holding_events<T>(work: impl FnOnce() -> T) -> (T, Vec<String>) {
    let (result, held) = vep_core::warning_hold::holding(work);
    (result, held.into_iter().map(|held| held.text).collect())
}

/// The `tracing` layer that hands every `WARN` event to the global sink as one
/// line, the message, then every other field as ` name=value`, or to the hold
/// the thread raising it is inside (`vep_core::warning_hold`). The field
/// [`RECORD_FIELD`] names the input record the event is about: it travels
/// with a held event and is no part of the line.
pub struct WarningLayer;

#[derive(Default)]
struct EventText {
    message: String,
    fields: Vec<(String, String)>,
    input_record: u64,
}

impl EventText {
    fn record(&mut self, field: &Field, value: String) {
        if field.name() == "message" {
            self.message = value;
        } else if field.name() == RECORD_FIELD {
            self.input_record = value.parse().unwrap_or(0);
        } else {
            self.fields.push((field.name().to_string(), value));
        }
    }

    fn line(self) -> String {
        let mut line = self.message;
        for (name, value) in self.fields {
            line.push_str(&format!(" {name}={value}"));
        }
        line
    }
}

impl Visit for EventText {
    fn record_debug(&mut self, field: &Field, value: &dyn std::fmt::Debug) {
        self.record(field, format!("{value:?}"));
    }

    fn record_str(&mut self, field: &Field, value: &str) {
        self.record(field, value.to_string());
    }
}

/// One event as the input record it names (0 for none) and its warning line:
/// its message, then every other field as ` name=value`.
fn event_line(event: &Event<'_>) -> (u64, String) {
    let mut text = EventText::default();
    event.record(&mut text);
    (text.input_record, text.line())
}

impl<S: Subscriber> Layer<S> for WarningLayer {
    fn on_event(&self, event: &Event<'_>, _ctx: Context<'_, S>) {
        if *event.metadata().level() != Level::WARN {
            return;
        }
        let (input_record, line) = event_line(event);
        if let Err(line) = vep_core::warning_hold::offer(input_record, line) {
            global().warn(&line);
        }
    }
}

/// The path of the warnings file beside an output file: `<output_file>_warnings.txt`
/// (BaseVEP.pm:890), `STDOUT_warnings.txt` for `-o STDOUT` as the user spelled it.
pub(crate) fn default_warning_file(output_file: &str) -> String {
    format!("{output_file}_warnings.txt")
}

/// The target a `--warning_file` value names: standard error for the literal
/// `stderr` in any case (BaseVEP.pm:885), else the file.
pub fn target_for(warning_file: &str) -> Target {
    if warning_file.eq_ignore_ascii_case("stderr") {
        Target::Stderr
    } else {
        Target::File(PathBuf::from(warning_file))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn file_sink(dir: &Path) -> (Sink, PathBuf) {
        let path = dir.join("out.txt_warnings.txt");
        let sink = Sink::new();
        sink.configure(Target::File(path.clone()), true, false);
        (sink, path)
    }

    /// `WARNING: ` is prefixed unless the text already opens with `warn` in any
    /// case (`$text =~ /^warn/i`), and a newline is added once.
    #[test]
    fn warning_prefix_follows_warning_msg() {
        let dir = tempfile::tempdir().unwrap();
        let (sink, path) = file_sink(dir.path());
        sink.warn("plain text");
        sink.warn("WARNING: already prefixed\n");
        sink.warn("Warned once more");
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "WARNING: plain text\nWARNING: already prefixed\nWarned once more\n"
        );
    }

    /// A built-in plugin's failures are held until the run ends, then written as
    /// one line per plugin instance: the first batch's text with the variants
    /// failed on summed over the batches, `variant` for one (`warning_msg`'s
    /// seen set leaves one line for a plugin dying the same way on every line,
    /// BaseVEP.pm:810-813); a second instance of the same name keeps its own
    /// line.
    #[test]
    fn plugin_failures_are_one_counted_line_per_plugin_at_the_end_of_the_run() {
        let sink = Sink::new();
        sink.configure(Target::Memory, true, false);
        sink.plugin_failure(
            "CADD",
            0,
            "no alt cell in a row of a_synthetic.tsv.gz; the first at 21:100",
            2,
        );
        sink.plugin_failure("REVEL", 1, "bad block; the first at 21:300", 1);
        sink.plugin_failure(
            "CADD",
            0,
            "no alt cell in a row of a_synthetic.tsv.gz; the first at 21:900",
            3,
        );
        sink.plugin_failure(
            "CADD",
            2,
            "no alt cell in a row of b_synthetic.tsv.gz; the first at 21:100",
            1,
        );
        assert!(sink.take_collected().is_empty(), "held until the run ends");
        sink.finish("in.vcf").unwrap();
        let texts: Vec<String> = sink
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect();
        assert_eq!(
            texts,
            [
                "WARNING: Plugin 'CADD' went wrong: no alt cell in a row of a_synthetic.tsv.gz; the first at \
                 21:100; 5 variants affected\n",
                "WARNING: Plugin 'REVEL' went wrong: bad block; the first at 21:300; 1 variant \
                 affected\n",
                "WARNING: Plugin 'CADD' went wrong: no alt cell in a row of b_synthetic.tsv.gz; the first at \
                 21:100; 1 variant affected\n",
            ]
        );
    }

    /// A repeated text is written once; two texts differing only in a trailing
    /// `line N` are one text (`$to_check =~ s/line \d+\n?$//`).
    #[test]
    fn warnings_deduplicate_ignoring_a_trailing_line_number() {
        let dir = tempfile::tempdir().unwrap();
        let (sink, path) = file_sink(dir.path());
        sink.warn("Skipped 1 empty line from line 3");
        sink.warn("Skipped 1 empty line from line 9");
        sink.warn("Skipped 1 empty line from line 3");
        sink.warn("another matter");
        sink.warn("another matter");
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "WARNING: Skipped 1 empty line from line 3\nWARNING: another matter\n"
        );
        assert_eq!(dedup_key("WARNING: x on line 12\n"), "WARNING: x on ");
        assert_eq!(
            dedup_key("WARNING: line 12 skipped (a): b\n"),
            "WARNING: line 12 skipped (a): b"
        );
        assert_eq!(dedup_key("WARNING: a line\n"), "WARNING: a line");
    }

    /// No warnings file exists until the first warning is written.
    #[test]
    fn warnings_file_is_opened_at_the_first_warning() {
        let dir = tempfile::tempdir().unwrap();
        let (sink, path) = file_sink(dir.path());
        assert!(!path.exists());
        sink.warn("now");
        assert!(path.is_file());
    }

    /// The run is handed a warnings file to test before its first record only
    /// where the sink was asked to have it tested, as a file `--warning_file`
    /// names is: the file under the default name is left to the first warning
    /// (BaseVEP.pm:817), and a file a warning opened, or could not open, is
    /// not handed over again. Configuring the sink for another run withdraws
    /// the request.
    #[test]
    fn the_warnings_file_is_handed_to_the_test_only_where_it_was_asked_for() {
        let dir = tempfile::tempdir().unwrap();
        let (sink, path) = file_sink(dir.path());
        assert_eq!(sink.file_to_test(), None);
        sink.test_file_before_first_record();
        assert_eq!(sink.file_to_test(), Some(path.clone()));
        sink.warn("now");
        assert!(path.is_file());
        assert_eq!(sink.file_to_test(), None);

        let unwritable = dir.path().join("absent").join("w.txt");
        sink.configure(Target::File(unwritable.clone()), true, false);
        assert_eq!(sink.file_to_test(), None);
        sink.test_file_before_first_record();
        assert_eq!(sink.file_to_test(), Some(unwritable));
        sink.warn("now");
        assert!(sink.failure().is_some());
        assert_eq!(sink.file_to_test(), None);

        sink.configure(Target::Stderr, true, false);
        sink.test_file_before_first_record();
        assert_eq!(sink.file_to_test(), None);
    }

    /// `skipped_variant_msg`: the line collapsed to single spaces and cut to 46
    /// bytes plus `...` past 50, `line N skipped` with a number, `variant
    /// skipped` without one or under `--input_data`.
    #[test]
    fn skipped_variant_lines_follow_skipped_variant_msg() {
        let dir = tempfile::tempdir().unwrap();
        let (sink, path) = file_sink(dir.path());
        let long = "21\t30165800\tmulti_ins_del\tN\t<INS>,<DEL>\t.\t.\tEND=30166000";
        sink.skipped_variant(
            "INS/DEL is not a supported structural variant type",
            Some(13),
            Some(long),
        );
        sink.skipped_variant(
            "Chromosome Z not found",
            Some(40),
            Some("Z\t1000\tunknown\tA\tG\t.\t.\t."),
        );
        sink.skipped_variant("no line at all", None, None);
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "WARNING: line 13 skipped (21 30165800 multi_ins_del N <INS>,<DEL> . . EN...): INS/DEL is not a supported structural variant type\n\
             WARNING: line 40 skipped (Z 1000 unknown A G . . .): Chromosome Z not found\n\
             WARNING: variant skipped no line at all\n"
        );
        let input_data = Sink::new();
        input_data.configure(Target::File(dir.path().join("id_warnings.txt")), true, true);
        input_data.skipped_variant(
            "Invalid allele string X/Y or possible parsing error",
            Some(1),
            Some("21\t1\t.\tX\tY"),
        );
        assert_eq!(
            std::fs::read_to_string(dir.path().join("id_warnings.txt")).unwrap(),
            "WARNING: variant skipped (21 1 . X Y): Invalid allele string X/Y or possible parsing error\n"
        );
        assert_eq!(quoted_line("a  b\t\tc"), "a b c");
        assert_eq!(quoted_line(&"x".repeat(50)), "x".repeat(50));
        assert_eq!(
            quoted_line(&"x".repeat(51)),
            format!("{}...", "x".repeat(46))
        );
    }

    /// The skipped-variants file: the header naming the input, then
    /// `Line %-5s\t%-40.40s\t%s` per skip, the header alone with no skip; an
    /// existing file is refused without `--force_overwrite`.
    #[test]
    fn skipped_variants_file_follows_dump_skipped_variants() {
        let dir = tempfile::tempdir().unwrap();
        let (sink, _) = file_sink(dir.path());
        let skipped = dir.path().join("skipped.txt");
        sink.reserve_skipped_variants_file(skipped.to_str().unwrap(), false)
            .unwrap();
        sink.skipped_variant(
            "deletion looks incomplete",
            Some(16),
            Some("21\t30165850\tdel_without_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL"),
        );
        sink.skipped_variant("short", Some(123456), Some("a b"));
        sink.finish("in.vcf").unwrap();
        assert_eq!(
            std::fs::read_to_string(&skipped).unwrap(),
            format!(
                "[VEP skipped the following variants from in.vcf]\n\
                 Line 16   \t21 30165850 del_without_end N <DEL> . . \tdeletion looks incomplete\n\
                 Line 123456\ta b{}\tshort\n",
                " ".repeat(37)
            )
        );
        let refused = sink
            .reserve_skipped_variants_file(skipped.to_str().unwrap(), false)
            .unwrap_err();
        assert_eq!(
            refused.to_string(),
            format!(
                "ERROR: [output] Skipped variants file {} already exists. Specify a different output file using --skipped_variants_file or overwrite existing file with --force_overwrite",
                skipped.display()
            )
        );
        sink.reserve_skipped_variants_file(skipped.to_str().unwrap(), true)
            .unwrap();
        sink.finish("in.vcf").unwrap();
        assert_eq!(
            std::fs::read_to_string(&skipped).unwrap(),
            "[VEP skipped the following variants from in.vcf]\n"
        );
    }

    /// The skipped-variants file is created at the end of the run, as Ensembl
    /// VEP creates it (BaseRunner.pm:364-380): a reserved path that cannot be
    /// created fails `finish` with Ensembl VEP's text (`:344-345`) under the output
    /// component, and a run that never finishes leaves no file.
    #[test]
    fn uncreatable_skipped_variants_file_fails_at_finish() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("no-such-dir").join("skipped.txt");
        let sink = Sink::new();
        sink.reserve_skipped_variants_file(path.to_str().unwrap(), false)
            .unwrap();
        assert!(!path.exists(), "nothing is created before the run ends");
        let refused = sink.finish("in.vcf").unwrap_err();
        assert_eq!(
            refused.to_string(),
            format!(
                "ERROR: [output] Could not write to skipped variants file {}",
                path.display()
            )
        );
        assert!(!path.exists());
        let kept = dir.path().join("skipped.txt");
        let sink = Sink::new();
        sink.reserve_skipped_variants_file(kept.to_str().unwrap(), false)
            .unwrap();
        sink.skipped_variant("dropped", Some(3), Some("a b"));
        assert!(
            !kept.exists(),
            "a run that fails before finish leaves no file"
        );
    }

    /// A skip is kept for the skipped-variants file only when one is reserved,
    /// and every skip line passes through the seen set, as `_seen_warnings` holds
    /// one entry per skip (BaseVEP.pm:810-813): two records skipped for one
    /// reason are two lines, their numbers apart; the identical line a second
    /// `--individual` copy of a record raises is written once; and a `variant
    /// skipped` text repeats once.
    #[test]
    fn skips_are_kept_only_for_a_reserved_file_and_enter_the_seen_set() {
        let dir = tempfile::tempdir().unwrap();
        let (sink, path) = file_sink(dir.path());
        sink.skipped_variant("same reason", Some(1), Some("a"));
        sink.skipped_variant("same reason", Some(2), Some("b"));
        sink.skipped_variant("same reason", Some(2), Some("b"));
        sink.skipped_variant("same reason", None, None);
        sink.skipped_variant("same reason", None, None);
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "WARNING: line 1 skipped (a): same reason\nWARNING: line 2 skipped (b): same reason\nWARNING: variant skipped same reason\n"
        );
        assert!(
            sink.skipped().is_empty(),
            "no skipped-variants file is reserved"
        );
        assert_eq!(
            sink.lock().seen.len(),
            3,
            "one entry per distinct skip text"
        );
    }

    /// A file target holds each distinct text once whatever `--quiet` says: the
    /// echo `print $fh` (BaseVEP.pm:819) then `warn` (:821) give standard error is
    /// not the file's.
    #[test]
    fn file_target_is_written_once_whatever_quiet_says() {
        let dir = tempfile::tempdir().unwrap();
        for quiet in [true, false] {
            let path = dir.path().join(format!("quiet-{quiet}.txt"));
            let sink = Sink::new();
            sink.configure(Target::File(path.clone()), quiet, false);
            sink.warn("once");
            sink.warn("once");
            assert_eq!(
                std::fs::read_to_string(&path).unwrap(),
                "WARNING: once\n",
                "quiet={quiet}"
            );
        }
    }

    /// A sink no run configured echoes no warning beside its standard-error
    /// copy; a configured one echoes each unless `--quiet`, whatever its target.
    #[test]
    fn only_a_configured_sink_that_is_not_quiet_echoes() {
        let sink = Sink::new();
        assert_eq!(sink.lock().target, Target::Stderr);
        assert!(!sink.lock().echo);
        for target in [Target::Stderr, Target::Memory] {
            sink.configure(target.clone(), false, false);
            assert!(sink.lock().echo, "{target:?}");
            sink.configure(target.clone(), true, false);
            assert!(!sink.lock().echo, "{target:?} under --quiet");
        }
    }

    /// `quoted_line` counts bytes, as `length` and `substr` do on Ensembl VEP's
    /// undecoded line (BaseVEP.pm:848-850): a 60-byte line of `é`s (two bytes each) is
    /// cut at byte 46, 23 characters, plus `...`; a line whose byte 46 falls inside a
    /// character is cut at the boundary before it, where Perl's `substr` keeps the
    /// character's lone lead byte (`docs/divergences/skipped-variants-file-header.md`);
    /// 50 bytes pass whole.
    #[test]
    fn quoted_line_cuts_by_bytes_at_a_character_boundary() {
        let accented = "é".repeat(30);
        assert_eq!(accented.len(), 60);
        assert_eq!(
            quoted_line(&accented),
            format!("{}...", "é".repeat(23)),
            "46 bytes then the ellipsis"
        );
        let straddling = format!("{}é{}", "x".repeat(45), "y".repeat(10));
        assert_eq!(
            quoted_line(&straddling),
            format!("{}...", "x".repeat(45)),
            "byte 46 is the second byte of é, so the cut moves to byte 45"
        );
        // A VCF record whose ID is `xy` and twenty `é`s collapses to 59 bytes, byte 46
        // the second byte of the nineteenth `é`: eighteen are kept, where Perl's
        // `substr` in Ensembl VEP writes the nineteenth's lead byte before the
        // ellipsis.
        let split_id = format!("21\t100\txy{}\tA\tG\t.\t.\t.", "é".repeat(20));
        assert_eq!(split_id.replace('\t', " ").len(), 59);
        assert_eq!(
            quoted_line(&split_id),
            format!("21 100 xy{}...", "é".repeat(18)),
            "45 bytes then the ellipsis, the split character dropped whole"
        );
        let fifty_bytes = format!("{}é", "x".repeat(48));
        assert_eq!(quoted_line(&fifty_bytes), fifty_bytes);
    }

    /// `quotable_prefix` keeps exactly what `quoted_line` renders: a line it cuts
    /// quotes as the whole line does, a line of at most 51 collapsed bytes is kept
    /// whole, and a multi-byte character is never split.
    #[test]
    fn quotable_prefix_quotes_as_the_whole_line_does() {
        let line = "21\t100\t.\tA\tG\t.\tPASS\tDP=10;AF=0.5\tGT\t0/1".to_string();
        let whitespace_run = format!("{} {}", "x".repeat(49), " ".repeat(30));
        let runs = format!("a{}b{}c", " \t ".repeat(20), "x".repeat(60));
        for (line, prefix_len) in [
            (line.clone(), line.len()),
            ("x".repeat(50), 50),
            ("x".repeat(51), 51),
            ("x".repeat(52), 51),
            // The tab that becomes the 51st byte, its run cut after it.
            (format!("{}\t\t\t{}", "x".repeat(50), "y".repeat(10)), 51),
            (whitespace_run.clone(), whitespace_run.len()),
            // 26 two-byte characters reach 52 bytes.
            ("é".repeat(60), 52),
            (format!("{}é{}", "x".repeat(50), "y".repeat(10)), 52),
            // `a`, one space for the run, `b`, then 48 of the `x`s.
            (runs, 110),
        ] {
            let prefix = quotable_prefix(line.as_bytes());
            assert_eq!(prefix.len(), prefix_len, "{line:?}");
            let prefix = std::str::from_utf8(prefix).expect("cut at a character boundary");
            assert_eq!(quoted_line(prefix), quoted_line(&line), "{line:?}");
        }
    }

    /// A `warn!` event becomes one line: its message, then every other field as
    /// ` name=value`, a string field unquoted, in the event's field order.
    #[test]
    fn warning_events_render_their_fields_after_the_message() {
        struct Capture(std::sync::Arc<Mutex<Vec<String>>>);
        impl<S: Subscriber> Layer<S> for Capture {
            fn on_event(&self, event: &Event<'_>, _ctx: Context<'_, S>) {
                self.0.lock().unwrap().push(event_line(event).1);
            }
        }
        use tracing_subscriber::layer::SubscriberExt;
        let lines = std::sync::Arc::new(Mutex::new(Vec::new()));
        let subscriber = tracing_subscriber::registry().with(Capture(lines.clone()));
        tracing::subscriber::with_default(subscriber, || {
            tracing::warn!(chr = "21", pos = 100u64, "could not read {}", "x");
            tracing::warn!("plain");
        });
        assert_eq!(
            lines.lock().unwrap().as_slice(),
            ["could not read x chr=21 pos=100", "plain"]
        );
    }

    /// Runs `work` under a subscriber that is [`WarningLayer`] alone, with the
    /// process-wide sink holding its warnings in memory, and returns what the
    /// sink took.
    fn sink_lines_of(work: impl FnOnce()) -> Vec<String> {
        use tracing_subscriber::layer::SubscriberExt;
        global().configure(Target::Memory, true, false);
        tracing::subscriber::with_default(tracing_subscriber::registry().with(WarningLayer), work);
        global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect()
    }

    /// The events a thread raises inside [`holding_events`] come back from the
    /// call in the order raised, fields rendered, and none reaches the sink;
    /// the events before and after the call reach it as they are raised.
    #[test]
    fn held_events_are_returned_in_order_and_kept_from_the_sink() {
        let _guard = sink_test_guard();
        let mut held = None;
        let sink = sink_lines_of(|| {
            tracing::warn!("before");
            held = Some(holding_events(|| {
                tracing::warn!("first");
                tracing::warn!(model = "NM_1.1", "second");
                tracing::info!("not a warning");
                7
            }));
            tracing::warn!("after");
        });
        assert_eq!(
            held.unwrap(),
            (
                7,
                vec!["first".to_string(), "second model=NM_1.1".to_string()]
            )
        );
        assert_eq!(sink, ["WARNING: before\n", "WARNING: after\n"]);
    }

    /// The field that names an event's input record is no part of its line,
    /// held or written: a held event carries the record beside its text, in
    /// the order raised, and an event raised outside a hold reaches the sink
    /// as its message and its other fields.
    #[test]
    fn the_record_field_travels_with_a_held_event_and_is_never_written() {
        use vep_core::warning_hold::{holding, Held};
        let _guard = sink_test_guard();
        let mut held = Vec::new();
        let sink = sink_lines_of(|| {
            ((), held) = holding(|| {
                tracing::warn!(input_record = 12u64, "of record twelve");
                tracing::warn!("of no record");
                tracing::warn!(input_record = 3u64, file = "x.gff", "of record three");
            });
            tracing::warn!(input_record = 5u64, file = "x.gff", "outside a hold");
        });
        let held_as = |input_record: u64, text: &str| Held {
            input_record,
            text: text.to_string(),
        };
        assert_eq!(
            held,
            [
                held_as(12, "of record twelve"),
                held_as(0, "of no record"),
                held_as(3, "of record three file=x.gff"),
            ]
        );
        assert_eq!(sink, ["WARNING: outside a hold file=x.gff\n"]);
    }

    /// A call inside another holds its own events and leaves the outer call's
    /// as they were, whether it returns or unwinds; the events of a call that
    /// unwinds are dropped, and once the outermost call has ended, by either
    /// way, an event reaches the sink again.
    #[test]
    fn a_hold_inside_another_keeps_its_own_events() {
        let _guard = sink_test_guard();
        let mut outer = Vec::new();
        let mut inner = Vec::new();
        let sink = sink_lines_of(|| {
            ((), outer) = holding_events(|| {
                tracing::warn!("outer first");
                ((), inner) = holding_events(|| tracing::warn!("inner"));
                let unwound = std::panic::catch_unwind(|| {
                    holding_events(|| {
                        tracing::warn!("dropped");
                        panic!("the held work unwinds");
                    })
                });
                assert!(unwound.is_err());
                tracing::warn!("outer second");
            });
            let unwound = std::panic::catch_unwind(|| {
                holding_events(|| {
                    tracing::warn!("dropped too");
                    panic!("the held work unwinds");
                })
            });
            assert!(unwound.is_err());
            tracing::warn!("after");
        });
        assert_eq!(inner, ["inner"]);
        assert_eq!(outer, ["outer first", "outer second"]);
        assert_eq!(sink, ["WARNING: after\n"]);
    }

    /// An event another thread raises while this one is inside
    /// [`holding_events`] is not held: it reaches the sink.
    #[test]
    fn another_threads_event_is_not_held() {
        use tracing_subscriber::layer::SubscriberExt;
        let _guard = sink_test_guard();
        let mut held = Vec::new();
        let sink = sink_lines_of(|| {
            ((), held) = holding_events(|| {
                tracing::warn!("this thread");
                std::thread::scope(|scope| {
                    scope.spawn(|| {
                        tracing::subscriber::with_default(
                            tracing_subscriber::registry().with(WarningLayer),
                            || tracing::warn!("another thread"),
                        );
                    });
                });
            });
        });
        assert_eq!(held, ["this thread"]);
        assert_eq!(sink, ["WARNING: another thread\n"]);
    }

    /// A warnings file that cannot be created is the run's failure, surfaced by
    /// `failure` and by `finish`, with Ensembl VEP's text under the output component.
    #[test]
    fn unwritable_warnings_file_fails_the_run() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("no-such-dir").join("w.txt");
        let sink = Sink::new();
        sink.configure(Target::File(path.clone()), true, false);
        assert!(sink.failure().is_none());
        sink.warn("first");
        let expected = format!(
            "ERROR: [output] Could not write to warnings file {}",
            path.display()
        );
        assert_eq!(sink.failure().unwrap().to_string(), expected);
        assert_eq!(sink.finish("in.vcf").unwrap_err().to_string(), expected);
    }

    /// The default warnings file sits beside the output under the name the user
    /// gave, and `stderr` in any case names standard error.
    #[test]
    fn warning_file_defaults_and_stderr_literal() {
        assert_eq!(default_warning_file("out.txt"), "out.txt_warnings.txt");
        assert_eq!(default_warning_file("stdout"), "stdout_warnings.txt");
        assert_eq!(target_for("STDERR"), Target::Stderr);
        assert_eq!(target_for("stderr"), Target::Stderr);
        assert_eq!(target_for("StdErr"), Target::Stderr);
        assert_eq!(target_for("w.txt"), Target::File(PathBuf::from("w.txt")));
    }

    /// A unit test that holds no guard cannot reach the process-wide sink,
    /// whether or not another test holds one; the test that holds the guard
    /// reaches it, and so does a thread that test starts.
    #[test]
    fn a_test_without_the_guard_cannot_reach_the_process_wide_sink() {
        let reaches_the_sink = || {
            global();
        };
        assert!(std::panic::catch_unwind(reaches_the_sink).is_err());

        let _in_turn = sink_test_guard();
        reaches_the_sink();
        std::thread::spawn(reaches_the_sink).join().unwrap();
        let another_test = std::thread::Builder::new()
            .name("module::tests::another_test".to_string())
            .spawn(reaches_the_sink)
            .unwrap();
        assert!(another_test.join().is_err());
    }

    /// The holder of the guard finds the process-wide sink as a process starts
    /// it, whatever the holder before it left there: a target in memory or in
    /// a file, an echo, a warnings file open or one that could not be opened,
    /// a reserved skipped-variants file, and warnings, skips and plugin
    /// failures held.
    #[test]
    fn the_guard_hands_over_the_process_wide_sink_as_a_process_starts_it() {
        let dir = tempfile::tempdir().unwrap();
        let skipped_file = dir.path().join("skipped.txt");
        let leave_a_run_behind = |target: Target| {
            let _earlier = sink_test_guard();
            let sink = global();
            sink.configure(target, false, true);
            sink.reserve_skipped_variants_file(skipped_file.to_str().unwrap(), false)
                .unwrap();
            sink.warn("left behind");
            sink.skipped_variant("left behind", Some(3), Some("21 100"));
            sink.plugin_failure("Left", 0, "behind", 1);
        };
        for target in [
            Target::Memory,
            Target::File(dir.path().join("warnings.txt")),
            Target::File(dir.path().join("absent").join("warnings.txt")),
        ] {
            leave_a_run_behind(target.clone());
            let _guard = sink_test_guard();
            let state = global().lock();
            assert_eq!(state.target, Target::Stderr, "after {target:?}");
            assert!(!state.echo && !state.input_data, "after {target:?}");
            assert!(
                state.file.is_none() && state.failure.is_none(),
                "after {target:?}"
            );
            assert!(
                state.seen.is_empty() && state.collected.is_empty(),
                "after {target:?}"
            );
            assert!(
                state.skipped.is_empty() && state.skipped_file.is_none(),
                "after {target:?}"
            );
            assert!(state.plugin_failures.is_empty(), "after {target:?}");
        }
    }
}
