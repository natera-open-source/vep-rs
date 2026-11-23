// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Counts, digests and row passes over whole files, across the host's cores.
//!
//! Each subcommand writes on standard output the bytes of the shell pipeline it stands
//! for: a count as `grep -c` prints it, digest lines as `sha256sum`, `md5sum` and
//! `xxh128sum` print them, the lines of an `awk` pass. `line-sum` stands for no pipeline:
//! it prints the count of a file's lines and the sum of their SHA256 digests, which is
//! the same wherever the lines stand. `README.md` states the contract of each.

mod contigs;
mod count;
mod digest;
mod fetched;
mod gzip;
mod line_sum;
mod metrics;
mod names;
mod pipeline;
mod strip;

use count::Rule;
use digest::{Algorithm, Body, Plan, Status};
use line_sum::LineSum;
use pipeline::{lines, thread_bound, Ends, Reading};
use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{self, Read, Write};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::PathBuf;
use std::process::ExitCode;

const USAGE: &str = "\
usage: vep-file-census count [--gunzip] RULE [RULE ...] FILE
         RULE: --starting TEXT | --not-starting TEXT | --containing TEXT
               | --info-field NAME | --extra-key NAME
       vep-file-census digest --algorithm {sha256,md5,xxh128} [--algorithm ...]
         [--names {lines,read,nul,words}] [--directory DIR] [--gunzip]
         [--without-lines-starting TEXT | --after-leading-lines-starting TEXT]
         [--chunk-bytes N] [FILE ...]
       vep-file-census strip-chr FILE
       vep-file-census metrics OUTPUT REPORT SECONDS
       vep-file-census contigs [--gunzip] FILE
       vep-file-census slice --records N --tally FILE [--label TEXT] [--table FILE]
         [--gunzip] FILE
       vep-file-census fetched --list TSV --root DIR [--jobs N] [--report TSV]
       vep-file-census line-sum [--gunzip] [--without-lines-starting TEXT] FILE [FILE ...]
";

/// How a run ends without its output.
enum Failure {
    /// A command line the tool cannot take: exit code 2.
    Usage(String),
    /// What failed, and the exit code.
    Run(String, u8),
}

fn usage(detail: impl Into<String>) -> Failure {
    Failure::Usage(detail.into())
}

fn complain(detail: &str) {
    eprintln!("ERROR: [vep-file-census] {detail}");
}

fn shown(name: &OsStr) -> String {
    name.to_string_lossy().into_owned()
}

/// The words of a command line after the subcommand.
struct Words {
    rest: std::vec::IntoIter<OsString>,
    options_ended: bool,
}

enum Word {
    /// An option's name and the value written after `=`, if any.
    Option(String, Option<OsString>),
    Operand(OsString),
}

impl Words {
    fn next(&mut self) -> Result<Option<Word>, Failure> {
        let Some(word) = self.rest.next() else {
            return Ok(None);
        };
        let bytes = word.as_bytes();
        if self.options_ended || !bytes.starts_with(b"-") || bytes == b"-" {
            return Ok(Some(Word::Operand(word)));
        }
        if bytes == b"--" {
            self.options_ended = true;
            return self.next();
        }
        if !bytes.starts_with(b"--") {
            return Err(usage(format!("unknown option {}", shown(&word))));
        }
        let (name, value) = match memchr::memchr(b'=', bytes) {
            Some(at) => (
                &bytes[..at],
                Some(OsString::from_vec(bytes[at + 1..].to_vec())),
            ),
            None => (bytes, None),
        };
        match std::str::from_utf8(name) {
            Ok(name) => Ok(Some(Word::Option(name.to_string(), value))),
            Err(_) => Err(usage(format!("unknown option {}", shown(&word)))),
        }
    }

    /// The value of the option `name`: what stood after `=`, or the next word.
    fn value(&mut self, name: &str, inline: Option<OsString>) -> Result<OsString, Failure> {
        inline
            .or_else(|| self.rest.next())
            .ok_or_else(|| usage(format!("{name} takes a value")))
    }

    fn flag(name: &str, inline: Option<OsString>) -> Result<(), Failure> {
        match inline {
            None => Ok(()),
            Some(_) => Err(usage(format!("{name} takes no value"))),
        }
    }
}

fn text_of(name: &str, value: OsString) -> Result<Vec<u8>, Failure> {
    let bytes = value.into_vec();
    if bytes.contains(&b'\n') {
        return Err(usage(format!("{name} takes a text of one line")));
    }
    Ok(bytes)
}

fn number_of(name: &str, value: &OsStr) -> Result<u64, Failure> {
    value
        .to_str()
        .and_then(|text| text.parse::<u64>().ok())
        .ok_or_else(|| {
            usage(format!(
                "{name} takes a whole number, not '{}'",
                shown(value)
            ))
        })
}

fn one_file(operands: Vec<OsString>) -> Result<OsString, Failure> {
    let mut operands = operands.into_iter();
    match (operands.next(), operands.next()) {
        (Some(file), None) => Ok(file),
        _ => Err(usage("one FILE is read")),
    }
}

fn threads() -> Result<usize, Failure> {
    thread_bound().map_err(Failure::Usage)
}

fn written(result: io::Result<()>) -> Result<(), Failure> {
    result.map_err(|error| Failure::Run(format!("standard output: {error}"), 1))
}

fn count(mut words: Words) -> Result<u8, Failure> {
    let mut gunzip = false;
    let mut rules = Vec::new();
    let mut operands = Vec::new();
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => operands.push(operand),
            Word::Option(name, inline) => {
                let rule: fn(Vec<u8>) -> Rule = match name.as_str() {
                    "--gunzip" => {
                        Words::flag(&name, inline)?;
                        gunzip = true;
                        continue;
                    }
                    "--starting" => Rule::Starting,
                    "--not-starting" => Rule::NotStarting,
                    "--containing" => Rule::Containing,
                    "--info-field" => Rule::InfoField,
                    "--extra-key" => Rule::ExtraKey,
                    _ => return Err(usage(format!("count: unknown option {name}"))),
                };
                let value = words.value(&name, inline)?;
                rules.push(rule(text_of(&name, value)?));
            }
        }
    }
    if rules.is_empty() {
        return Err(usage("count takes a rule"));
    }
    let path = one_file(operands)?;
    let threads = threads()?;
    let mut totals = vec![0u64; rules.len()];
    let reading = if gunzip {
        Reading::Gunzip
    } else {
        Reading::Plain { from: 0 }
    };
    let mut out = io::stdout().lock();
    let file = match File::open(&path) {
        Ok(file) => file,
        Err(error) => {
            complain(&format!("{}: {error}", shown(&path)));
            if !gunzip {
                // grep prints no count of a file it cannot open.
                return Ok(2);
            }
            // gzip writes nothing of it, and grep counts that.
            for _ in &rules {
                written(writeln!(out, "0"))?;
            }
            return Ok(1);
        }
    };
    let read = lines(
        &file,
        reading,
        threads,
        Ends::Newline,
        |data| count::counts(&rules, data),
        |part| {
            for (total, more) in totals.iter_mut().zip(part) {
                *total += more;
            }
            Ok(())
        },
    );
    if let Err(error) = &read {
        complain(&format!("{}: {error}", shown(&path)));
    }
    for total in &totals {
        written(writeln!(out, "{total}"))?;
    }
    let none = rules
        .iter()
        .zip(&totals)
        .any(|(rule, total)| rule.is_grep() && *total == 0);
    Ok(match read {
        // grep leaves with 2 on a file it cannot read to its end; a pipeline whose gzip
        // fails leaves with gzip's 1.
        Err(_) if !gunzip => 2,
        Err(_) => 1,
        Ok(_) if none => 1,
        Ok(finish) if finish.garbage => 2,
        Ok(_) => 0,
    })
}

fn digest(mut words: Words) -> Result<u8, Failure> {
    let mut algorithms = Vec::new();
    let mut form = None;
    let mut gunzip = false;
    let mut body = Body::Whole;
    let mut chunk_bytes = None;
    let mut directory = None;
    let mut operands = Vec::new();
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => operands.push(operand.into_vec()),
            Word::Option(name, inline) => match name.as_str() {
                "--gunzip" => {
                    Words::flag(&name, inline)?;
                    gunzip = true;
                }
                "--directory" => directory = Some(words.value(&name, inline)?),
                "--algorithm" => {
                    let value = words.value(&name, inline)?;
                    let algorithm = value.to_str().and_then(Algorithm::parse).ok_or_else(|| {
                        usage(format!(
                            "--algorithm takes sha256, md5 or xxh128, not '{}'",
                            shown(&value)
                        ))
                    })?;
                    algorithms.push(algorithm);
                }
                "--names" => {
                    let value = words.value(&name, inline)?;
                    form = Some(value.to_str().and_then(names::Form::parse).ok_or_else(|| {
                        usage(format!(
                            "--names takes lines, read, nul or words, not '{}'",
                            shown(&value)
                        ))
                    })?);
                }
                "--without-lines-starting" | "--after-leading-lines-starting" => {
                    if body != Body::Whole {
                        return Err(usage("one body is digested"));
                    }
                    let text = text_of(&name, words.value(&name, inline)?)?;
                    body = if name == "--without-lines-starting" {
                        Body::WithoutLinesStarting(text)
                    } else {
                        Body::AfterLeadingLinesStarting(text)
                    };
                }
                "--chunk-bytes" => {
                    let value = words.value(&name, inline)?;
                    let bytes = number_of(&name, &value)?;
                    if bytes == 0 {
                        return Err(usage("--chunk-bytes takes 1 or more"));
                    }
                    chunk_bytes = Some(bytes);
                }
                _ => return Err(usage(format!("digest: unknown option {name}"))),
            },
        }
    }
    if algorithms.is_empty() {
        return Err(usage("digest takes an --algorithm"));
    }
    let from_offset = matches!(body, Body::AfterLeadingLinesStarting(_)) || chunk_bytes.is_some();
    if from_offset && (gunzip || matches!(body, Body::WithoutLinesStarting(_))) {
        return Err(usage(
            "--after-leading-lines-starting and --chunk-bytes read a file's own bytes",
        ));
    }
    let names = match form {
        Some(_) if !operands.is_empty() => {
            return Err(usage(
                "names come from standard input or as operands, not both",
            ))
        }
        Some(form) => {
            let mut text = Vec::new();
            io::stdin()
                .lock()
                .read_to_end(&mut text)
                .map_err(|error| Failure::Run(format!("standard input: {error}"), 1))?;
            names::parse(&text, form).map_err(|detail| Failure::Run(detail, 1))?
        }
        None if operands.is_empty() => return Err(usage("digest takes a FILE or --names")),
        None => operands,
    };
    let plan = Plan {
        algorithms,
        gunzip,
        body,
        chunk_bytes,
    };
    let threads = threads()?;
    if let Some(directory) = directory {
        // The names are those of a listing made inside the directory.
        std::env::set_current_dir(&directory)
            .map_err(|error| Failure::Run(format!("{}: {error}", shown(&directory)), 1))?;
    }
    let mut out = io::BufWriter::new(io::stdout().lock());
    let statuses = digest::run(
        &plan,
        &names,
        form.is_some_and(names::Form::through_xargs),
        threads,
        &mut out,
        &mut |detail: &str| complain(detail),
    )
    .and_then(|statuses| out.flush().map(|()| statuses));
    let statuses =
        statuses.map_err(|error| Failure::Run(format!("standard output: {error}"), 1))?;
    let code = |status: &Status| match status {
        Status::Digested => 0,
        Status::Warned => 2,
        Status::NoLine | Status::Failed => 1,
    };
    Ok(match form {
        // The digest program's own status over its operands.
        None => statuses
            .iter()
            .map(code)
            .fold(0, |worst, code| match (worst, code) {
                (1, _) | (_, 1) => 1,
                (worst, code) => worst.max(code),
            }),
        // xargs leaves with 123 where a start of the program failed.
        Some(form) if form.through_xargs() => {
            if statuses.iter().all(|status| *status == Status::Digested) {
                0
            } else {
                123
            }
        }
        // A loop's status is that of its last pass.
        Some(_) => statuses.last().map_or(0, code),
    })
}

fn strip_chr(mut words: Words) -> Result<u8, Failure> {
    let mut operands = Vec::new();
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => operands.push(operand),
            Word::Option(name, _) => {
                return Err(usage(format!("strip-chr: unknown option {name}")))
            }
        }
    }
    let path = one_file(operands)?;
    let threads = threads()?;
    let mut out = io::stdout().lock();
    let read = File::open(&path).and_then(|file| {
        lines(
            &file,
            Reading::Auto,
            threads,
            Ends::Newline,
            strip::strip,
            |made: Vec<u8>| out.write_all(&made),
        )
    });
    match read {
        Ok(_) => Ok(0),
        Err(error) => Err(Failure::Run(format!("{}: {error}", shown(&path)), 1)),
    }
}

fn metrics(mut words: Words) -> Result<u8, Failure> {
    let mut operands = Vec::new();
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => operands.push(operand),
            Word::Option(name, _) => return Err(usage(format!("metrics: unknown option {name}"))),
        }
    }
    let [output, report, seconds] = <[OsString; 3]>::try_from(operands)
        .map_err(|_| usage("metrics takes OUTPUT, REPORT and SECONDS"))?;
    let threads = threads()?;
    let elapsed = seconds.to_str().and_then(metrics::py_int).ok_or_else(|| {
        Failure::Run(
            format!("SECONDS is no whole number: '{}'", shown(&seconds)),
            1,
        )
    })?;
    let counts = File::open(&output)
        .and_then(|file| metrics::counts(&file, threads))
        .map_err(|error| Failure::Run(format!("{}: {error}", shown(&output)), 1))?;
    std::fs::write(&report, counts.report(&elapsed))
        .map_err(|error| Failure::Run(format!("{}: {error}", shown(&report)), 1))?;
    written(io::stdout().lock().write_all(counts.summary().as_bytes()))?;
    Ok(0)
}

/// What `contigs` and `slice` take beside their own options: `--gunzip` and one FILE.
fn reading_of(gunzip: bool) -> Reading {
    if gunzip {
        Reading::Gunzip
    } else {
        Reading::Plain { from: 0 }
    }
}

/// The exit code of a pass that writes what it read before a failure, as a pipeline whose
/// reader failed leaves with the reader's status.
fn pass_code(read: &io::Result<pipeline::Finish>, path: &OsStr) -> u8 {
    match read {
        Ok(finish) if finish.garbage => 2,
        Ok(_) => 0,
        Err(error) => {
            complain(&format!("{}: {error}", shown(path)));
            1
        }
    }
}

fn contigs(mut words: Words) -> Result<u8, Failure> {
    let mut gunzip = false;
    let mut operands = Vec::new();
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => operands.push(operand),
            Word::Option(name, inline) if name == "--gunzip" => {
                Words::flag(&name, inline)?;
                gunzip = true;
            }
            Word::Option(name, _) => return Err(usage(format!("contigs: unknown option {name}"))),
        }
    }
    let path = one_file(operands)?;
    let threads = threads()?;
    let mut tally = contigs::Tally::default();
    let read = File::open(&path).and_then(|file| {
        lines(
            &file,
            reading_of(gunzip),
            threads,
            Ends::Newline,
            contigs::Tally::of,
            |part| {
                tally.merge(part);
                Ok(())
            },
        )
    });
    let mut out = io::BufWriter::new(io::stdout().lock());
    for (contig, records) in &tally.contigs {
        written(
            out.write_all(contig)
                .and_then(|()| writeln!(out, "\t{records}")),
        )?;
    }
    written(out.flush())?;
    Ok(pass_code(&read, &path))
}

fn slice(mut words: Words) -> Result<u8, Failure> {
    let mut gunzip = false;
    let mut records = None;
    let mut tally = None;
    let mut label = Vec::new();
    let mut table = None;
    let mut operands = Vec::new();
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => operands.push(operand),
            Word::Option(name, inline) => match name.as_str() {
                "--gunzip" => {
                    Words::flag(&name, inline)?;
                    gunzip = true;
                }
                "--records" => {
                    let value = words.value(&name, inline)?;
                    records = Some(number_of(&name, &value)?);
                }
                "--tally" => tally = Some(words.value(&name, inline)?),
                "--label" => label = words.value(&name, inline)?.into_vec(),
                "--table" => table = Some(words.value(&name, inline)?),
                _ => return Err(usage(format!("slice: unknown option {name}"))),
            },
        }
    }
    let (Some(records), Some(tally)) = (records, tally) else {
        return Err(usage("slice takes --records and --tally"));
    };
    let path = one_file(operands)?;
    let threads = threads()?;
    let tally_text = std::fs::read(&tally)
        .map_err(|error| Failure::Run(format!("{}: {error}", shown(&tally)), 1))?;
    let mut slice = contigs::Slice::new(records, &tally_text);
    let mut out = io::stdout().lock();
    let mut kept = Vec::new();
    let read = File::open(&path).and_then(|file| {
        lines(
            &file,
            reading_of(gunzip),
            threads,
            Ends::Newline,
            <[u8]>::to_vec,
            |run: Vec<u8>| {
                kept.clear();
                slice.take(&run, &mut kept);
                out.write_all(&kept)
            },
        )
    });
    // awk opens the table at the first line it appends, so no line leaves no file.
    let rows = slice.table(&label);
    if let Some(table) = table.filter(|_| !rows.is_empty()) {
        std::fs::OpenOptions::new()
            .append(true)
            .create(true)
            .open(&table)
            .and_then(|mut file| file.write_all(&rows))
            .map_err(|error| Failure::Run(format!("{}: {error}", shown(&table)), 1))?;
    }
    Ok(pass_code(&read, &path))
}

fn fetched(mut words: Words) -> Result<u8, Failure> {
    let mut list = None;
    let mut root = None;
    let mut jobs = None;
    let mut report = None;
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => {
                return Err(usage(format!(
                    "fetched: unexpected word {}",
                    shown(&operand)
                )))
            }
            Word::Option(name, inline) => match name.as_str() {
                "--list" => list = Some(PathBuf::from(words.value(&name, inline)?)),
                "--root" => root = Some(PathBuf::from(words.value(&name, inline)?)),
                "--report" => report = Some(PathBuf::from(words.value(&name, inline)?)),
                "--jobs" => {
                    let value = words.value(&name, inline)?;
                    let count = value
                        .to_str()
                        .and_then(|text| text.parse::<i64>().ok())
                        .ok_or_else(|| {
                            usage(format!(
                                "--jobs takes a whole number, not '{}'",
                                shown(&value)
                            ))
                        })?;
                    jobs = Some(count);
                }
                _ => return Err(usage(format!("fetched: unknown option {name}"))),
            },
        }
    }
    let (Some(list), Some(root)) = (list, root) else {
        return Err(usage("fetched takes --list and --root"));
    };
    let jobs = match jobs {
        Some(count) if count < 1 => {
            return Err(usage(format!("--jobs takes 1 or more, not {count}")))
        }
        Some(count) => usize::try_from(count).unwrap_or(usize::MAX),
        None => threads()?.min(32),
    };
    if !root.is_dir() {
        return Err(usage(format!("the root {} is no folder", root.display())));
    }
    let rows = fetched::listed(&list).map_err(Failure::Usage)?;
    let verdicts = fetched::verdicts(&root, &rows, jobs)
        .map_err(|error| Failure::Run(error.to_string(), 2))?;
    if let Some(report) = report {
        let text: String = rows
            .iter()
            .zip(&verdicts)
            .map(|(row, word)| format!("{}\t{}\t{word}\n", row.key, row.size_text))
            .collect();
        std::fs::write(&report, text)
            .map_err(|error| Failure::Run(format!("{}: {error}", report.display()), 2))?;
    }
    let mut out = io::BufWriter::new(io::stdout().lock());
    for (row, word) in rows.iter().zip(&verdicts) {
        if word == fetched::DIFFERS || word == fetched::BYTES_ALONE {
            written(writeln!(out, "{}\t{word}", row.key))?;
        }
    }
    written(writeln!(out, "{}", fetched::held_line(&verdicts)))?;
    written(out.flush())?;
    Ok(if verdicts.iter().any(|word| word == fetched::DIFFERS) {
        3
    } else {
        0
    })
}

fn line_sum(mut words: Words) -> Result<u8, Failure> {
    let mut gunzip = false;
    let mut without = None;
    let mut operands = Vec::new();
    while let Some(word) = words.next()? {
        match word {
            Word::Operand(operand) => operands.push(operand),
            Word::Option(name, inline) => match name.as_str() {
                "--gunzip" => {
                    Words::flag(&name, inline)?;
                    gunzip = true;
                }
                "--without-lines-starting" => {
                    without = Some(text_of(&name, words.value(&name, inline)?)?);
                }
                _ => return Err(usage(format!("line-sum: unknown option {name}"))),
            },
        }
    }
    if operands.is_empty() {
        return Err(usage("line-sum takes a FILE"));
    }
    let threads = threads()?;
    let mut out = io::stdout().lock();
    let mut code = 0;
    for path in &operands {
        let mut total = LineSum::default();
        let read = File::open(path).and_then(|file| {
            lines(
                &file,
                reading_of(gunzip),
                threads,
                Ends::NewlineOrReturn,
                |data| LineSum::of(data, without.as_deref()),
                |part| {
                    total.merge(part);
                    Ok(())
                },
            )
        });
        // A sum of the lines read before a failure would pass for the file's.
        let finish = read.map_err(|error| Failure::Run(format!("{}: {error}", shown(path)), 1))?;
        written(writeln!(out, "{}\t{}", total.lines, total.hex()))?;
        if finish.garbage {
            code = 2;
        }
    }
    Ok(code)
}

fn run(args: Vec<OsString>) -> Result<u8, Failure> {
    let mut rest = args.into_iter();
    let Some(subcommand) = rest.next() else {
        return Err(usage("a subcommand is named"));
    };
    let words = Words {
        rest,
        options_ended: false,
    };
    match subcommand.to_str() {
        Some("-h" | "--help") => {
            print!("{USAGE}");
            Ok(0)
        }
        Some("count") => count(words),
        Some("digest") => digest(words),
        Some("strip-chr") => strip_chr(words),
        Some("metrics") => metrics(words),
        Some("contigs") => contigs(words),
        Some("slice") => slice(words),
        Some("fetched") => fetched(words),
        Some("line-sum") => line_sum(words),
        _ => Err(usage(format!("unknown subcommand {}", shown(&subcommand)))),
    }
}

fn main() -> ExitCode {
    match run(std::env::args_os().skip(1).collect()) {
        Ok(code) => ExitCode::from(code),
        Err(Failure::Usage(detail)) => {
            complain(&detail);
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
        Err(Failure::Run(detail, code)) => {
            complain(&detail);
            ExitCode::from(code)
        }
    }
}
