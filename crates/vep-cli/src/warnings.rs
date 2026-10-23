// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The warning sink and the skip registry: Ensembl VEP's `warning_msg`,
//! `skipped_variant_msg` and `warning_fh` (BaseVEP.pm:802-901) and its
//! `--skipped_variants_file` (BaseRunner.pm:331-380).
//!
//! Every warning of a run goes through one [`Sink`]: the `warn!` events of every
//! crate, routed here by [`WarningLayer`], and the skipped-variant events the
//! pipeline raises. A warning is written as `WARNING: <text>` to the warnings
//! file (`--warning_file`, default `<output_file>_warnings.txt`, the literal
//! `STDERR` in any case for standard error) and echoed on standard error unless
//! `--quiet`, so a standard-error target carries it twice; a text seen before, a
//! trailing `line N` set aside, is written once. The warnings file is opened at
//! the first warning, so a run without one leaves no file.
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

use crate::config::CliError;

/// The line shown in a skipped-variant warning is cut past this many bytes:
/// `$maxChar = 50` (BaseVEP.pm:848) measured by `length` on a line Perl reads
/// with no encoding layer, so it counts bytes.
const QUOTED_LINE_BYTES: usize = 50;

/// The width of the line column of the skipped-variants file (`%-40.40s`).
const SKIPPED_LINE_WIDTH: usize = 40;

/// Where warnings are written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Target {
    /// Standard error alone: `--warning_file STDERR`, or no run configured yet.
    Stderr,
    /// A file, opened at the first warning.
    File(PathBuf),
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
    quiet: bool,
    /// `--input_data` is given: a skip names no line number.
    input_data: bool,
    seen: HashSet<String>,
    /// The skips to list in the skipped-variants file, kept only while one is
    /// reserved: a run without `--skipped_variants_file` holds no skip.
    skipped: Vec<Skipped>,
    /// The `--skipped_variants_file` path, checked before the run and written at
    /// its end.
    skipped_file: Option<PathBuf>,
}

impl Default for State {
    /// Standard error, not quiet, nothing seen or skipped.
    fn default() -> Self {
        State {
            target: Target::Stderr,
            file: None,
            failure: None,
            quiet: false,
            input_data: false,
            seen: HashSet::new(),
            skipped: Vec::new(),
            skipped_file: None,
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
    /// A sink writing to standard error, not quiet.
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

    /// Points the sink at a run's target before the run starts; the state of an
    /// earlier run (its seen set and skips) is dropped.
    pub fn configure(&self, target: Target, quiet: bool, input_data: bool) {
        *self.lock() = State {
            target,
            quiet,
            input_data,
            ..State::default()
        };
    }

    /// `warning_msg` (BaseVEP.pm:802-822): `WARNING: ` is prefixed unless the text
    /// already opens with `warn` in any case, a newline is appended unless present,
    /// a text seen before with the same dedup key is dropped, and the text goes to
    /// the target and to standard error unless quiet, so a standard-error target
    /// carries it twice, as `print $fh` (:819) then `warn` (:821) write it.
    pub fn warn(&self, text: &str) {
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
        if !self.lock().seen.insert(dedup_key(&text)) {
            return;
        }
        self.write(text);
    }

    /// Writes one finished warning text to the target, then to standard error
    /// unless quiet, a standard-error target included. A file that cannot be
    /// created or written is kept as the run's failure ([`Sink::failure`], the
    /// `ERROR: [output]` line the run ends with); the standard-error echo follows
    /// `--quiet` all the same.
    fn write(&self, text: String) {
        let mut state = self.lock();
        match state.target.clone() {
            Target::Stderr => eprint!("{text}"),
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
        if !state.quiet {
            eprint!("{text}");
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
        self.warn(&text);
        let mut state = self.lock();
        if state.skipped_file.is_some() {
            state.skipped.push(Skipped {
                line,
                line_number,
                description: description.to_string(),
            });
        }
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

    /// Ends the run: writes the skipped-variants file (BaseRunner.pm:364-380, the
    /// header line then one row per skip, the header alone when nothing was
    /// skipped) and closes the warnings file. The warnings-file failure, if any,
    /// is returned after the skipped-variants file is written.
    pub fn finish(&self, input_file: &str) -> Result<(), CliError> {
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

/// The `output` component's error for a skipped-variants file that cannot be
/// created or written, with Ensembl VEP's text.
fn skipped_file_unwritable(path: &Path) -> CliError {
    CliError::new(
        "output",
        format!(
            "Could not write to skipped variants file {}",
            path.display()
        ),
    )
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
/// warning event and skip, standard error alone until configured.
pub fn global() -> &'static Sink {
    &GLOBAL
}

/// The guard a test that reads the process-wide sink, or reserves its
/// skipped-variants file, holds for its run, so such tests run in turn. A test
/// that only reads what one record reports gives it a [`Sink`] of its own
/// instead, which no other test can reach.
#[cfg(test)]
pub(crate) fn sink_test_guard() -> MutexGuard<'static, ()> {
    static SINK_TESTS: Mutex<()> = Mutex::new(());
    SINK_TESTS
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// The `tracing` layer that hands every `WARN` event to the global sink as one
/// line: the message, then every other field as ` name=value`.
pub struct WarningLayer;

#[derive(Default)]
struct EventText {
    message: String,
    fields: Vec<(String, String)>,
}

impl EventText {
    fn record(&mut self, field: &Field, value: String) {
        if field.name() == "message" {
            self.message = value;
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

/// One event as a warning line: its message, then every other field as
/// ` name=value`.
fn event_line(event: &Event<'_>) -> String {
    let mut text = EventText::default();
    event.record(&mut text);
    text.line()
}

impl<S: Subscriber> Layer<S> for WarningLayer {
    fn on_event(&self, event: &Event<'_>, _ctx: Context<'_, S>) {
        if *event.metadata().level() != Level::WARN {
            return;
        }
        global().warn(&event_line(event));
    }
}

/// The path of the warnings file beside an output file: `<output_file>_warnings.txt`
/// (BaseVEP.pm:890), `STDOUT_warnings.txt` for `-o STDOUT` as the user spelled it.
pub fn default_warning_file(output_file: &str) -> String {
    format!("{output_file}_warnings.txt")
}

/// The target a `--warning_file` value names: standard error for the literal
/// `stderr` in any case (BaseVEP.pm:884), else the file.
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

    /// `quoted_line` counts bytes, as `length` and `substr` do on Perl's undecoded
    /// line (BaseVEP.pm:848-850): a 60-byte line of `é`s (two bytes each) is cut at
    /// byte 46, 23 characters, plus `...`; a line whose byte 46 falls inside a
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
        // `substr` writes the nineteenth's lead byte before the ellipsis.
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
                self.0.lock().unwrap().push(event_line(event));
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
}
