// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep filter`: filters an annotated output line by line, as Ensembl VEP's
//! `filter_vep` does (`filter_vep` at release/116.2, whose line numbers the
//! citations below name), over the filter grammar of [`crate::filter_set`].
//!
//! The input is a VCF or tab output of `vep` or of Ensembl VEP, read from
//! `--input_file` or standard input; its header lines are written through and
//! name the fields a filter may read ([`crate::filter_input`]). A VCF line
//! passes when any entry of its consequence field passes the filters, and a
//! tab line when its row does (282-291); `--only_matched` keeps the passing
//! entries alone in the field (294-297); `--soft_filter` writes every line and
//! marks the FILTER column `filter_vep_pass` or `filter_vep_fail`, adding the
//! two `##FILTER` lines ahead of the column line (201-206, 305-317). `--count`
//! prints the number of passing lines alone (343), `--list` the sorted field
//! names and stops (219-223). `--start`, `--limit` and `--test` bound the
//! passing lines written and the lines read (195, 299-301, 322).
//!
//! The command line is read as Getopt::Long reads `filter_vep`'s (72-103): an
//! option after one dash or two, in any case, abbreviated to any prefix that
//! names one option; a repeated option keeps its last value, and every
//! `--filter` given is kept (102). Every word is read as the bytes the shell
//! passed ([`filter_input::decode_word`]), so a value, a field name or a file
//! name outside UTF-8 compares byte for byte with the input, as Perl's
//! `@ARGV` byte strings do. A field name resolves against each line's
//! own keys through the synonym cache
//! ([`crate::filter_set::FilterSet::deferred`]), limited after the first
//! failed search except on a tab line with an `Extra` column (244, 273). Every
//! error is one line `ERROR: [filter] <message>` on standard error and exit 1,
//! the message Ensembl VEP's where it has one; a command-line syntax error is
//! clap's, exit 2.
//!
//! Five readings of `filter_vep` are corrected rather than reproduced: a
//! field that resolves as the prefix of one key reads its value on the first
//! line as well (`FilterSet.pm:570-571` records the synonym and reads
//! nothing); an INFO flag without `=` is a defined, empty field
//! (`filter_vep:253-257` leaves it undefined); an INFO or `Extra` value is
//! everything after its key's first `=` (`filter_vep:253, 402` keep the text
//! before the second `=` alone); `--only_matched` rewrites each INFO item
//! under the consequence field's name to its own passing entries
//! ([`crate::filter_input::ParsedLine::only_matched`]); and a VCF line's INFO
//! column is its eighth word with leading blanks skipped (`filter_vep:252`
//! keeps a leading blank's empty first field and reads the FILTER column).
//! `detect_format` at `filter_vep:425-429` classes a first data line starting
//! with a blank as a tab output; `vep filter` reads the file as VCF.
//!
//! `--ontology` rewrites `Consequence is <term>` to `Consequence is_child
//! <term>` (113-120) and answers `is_child` from the built-in Sequence
//! Ontology table ([`vep_core::so_ontology`]) instead of an Ensembl database;
//! the flags that name a database connection are refused.

use std::borrow::Cow;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, ErrorKind, Write};
use std::sync::LazyLock;

use clap::error::ErrorKind as ClapErrorKind;
use clap::{ArgAction, CommandFactory, Parser};
use regex::Regex;
use vep_core::so_ontology;

use crate::argv;
use crate::config::CliError;
use crate::filter_input::{self, Format, Headers, ParsedLine};
use crate::filter_set::FilterSet;

/// The default `--limit`, Ensembl VEP's `1e12` (63).
const NO_LIMIT: i64 = 1_000_000_000_000;

const COMPONENT: &str = "filter";

/// A `con... is` or `con... is_child` unit, in any case, which `--ontology`
/// rewrites to `Consequence is_child` (116); its letters and blank are ASCII,
/// as the pattern's are over a byte string.
static CONSEQUENCE_IS: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?-u)(?i)con[a-z]+\sis(_child)?").unwrap());

/// `--ontology`'s rewrite of one filter string (116-117): the first unit
/// [`CONSEQUENCE_IS`] finds, and every later repeat of that unit's exact
/// text, become `Consequence is_child`, as `s/$&/Consequence is_child/g`
/// replaces the matched text wherever it recurs and no other spelling.
fn ontology_rewrite(filter: &str) -> String {
    match CONSEQUENCE_IS.find(filter) {
        Some(found) => filter.replace(found.as_str(), "Consequence is_child"),
        None => filter.to_string(),
    }
}

fn error(message: impl Into<String>) -> CliError {
    CliError::new(COMPONENT, message)
}

/// The command line of `vep filter`, `filter_vep`'s options (72-103). A
/// repeated option overrides itself, as Getopt::Long assigns each occurrence
/// of a scalar option over the last (`Getopt/Long.pm:677`); `--filter`
/// appends, as its array destination does (102).
#[derive(Debug, Parser)]
#[command(
    name = "vep filter",
    about = "Filter an annotated VCF or tab output by its fields, as Ensembl VEP's filter_vep does",
    long_about = "Filter an annotated VCF or tab output by its fields, as Ensembl VEP's filter_vep does.\n\nA filter is `<field> <operator> <value>`; several are joined by `and` and `or`, grouped by parentheses and inverted by `not`. The operators are `is` (`eq`, `=`), `ne` (`!=`), `gt` (`>`), `lt` (`<`), `gte` (`>=`), `lte` (`<=`), `exists` (`ex`, `defined`), `nex`, `match` (`matches`, `regex`, `re`), `nre`, `in` and `is_child`. A field alone tests that it is defined; a value starting with `#` names another field. A VCF line passes when any entry of its consequence field passes; a tab line passes when its row does.",
    disable_version_flag = true,
    arg_required_else_help = false,
    args_override_self = true
)]
pub struct FilterArgs {
    /// The annotated output to filter, a VCF or tab file, read from standard input when absent; a `.gz` file, or any gzip stream, is decompressed
    #[arg(
        long = "input_file",
        short = 'i',
        value_name = "FILE",
        help_heading = "Input and output"
    )]
    pub input_file: Option<String>,

    /// The input format, `vcf` or `tab`, detected from the first data line when absent
    #[arg(
        long = "format",
        value_name = "vcf|tab",
        help_heading = "Input and output"
    )]
    pub format: Option<String>,

    /// Read the input as gzip-compressed whatever its name
    #[arg(long = "gz", help_heading = "Input and output")]
    pub gz: bool,

    /// Where the filtered lines go; `stdout` (in any case) is standard output
    #[arg(
        long = "output_file",
        short = 'o',
        value_name = "FILE",
        default_value = "stdout",
        help_heading = "Input and output"
    )]
    pub output_file: String,

    /// Overwrite an existing output file, which is otherwise an error
    #[arg(long = "force_overwrite", help_heading = "Input and output")]
    pub force_overwrite: bool,

    /// The INFO key the consequence annotations of a VCF input sit under
    #[arg(
        long = "vcf_info_field",
        value_name = "KEY",
        default_value = "CSQ",
        help_heading = "Input and output"
    )]
    pub vcf_info_field: String,

    /// A filter; repeated flags are joined by `and`
    #[arg(
        long = "filter",
        short = 'f',
        value_name = "FILTER",
        action = ArgAction::Append,
        allow_hyphen_values = true,
        help_heading = "Filters"
    )]
    pub filter: Vec<String>,

    /// Print the field names the input's header block names, sorted, and stop
    #[arg(long = "list", short = 'l', help_heading = "Filters")]
    pub list: bool,

    /// Print the number of passing lines alone
    #[arg(long = "count", short = 'c', help_heading = "Filters")]
    pub count: bool,

    /// In a VCF input, keep only the passing entries of the consequence field
    #[arg(long = "only_matched", help_heading = "Filters")]
    pub only_matched: bool,

    /// In a VCF input, write every line and mark its FILTER column `filter_vep_pass` or `filter_vep_fail`
    #[arg(long = "soft_filter", help_heading = "Filters")]
    pub soft_filter: bool,

    /// Do not warn about VCF lines without the consequence field
    #[arg(long = "no_missing_consequence_warning", help_heading = "Filters")]
    pub no_missing_consequence_warning: bool,

    /// Skip the first N-1 passing lines
    #[arg(
        long = "start",
        short = 's',
        value_name = "N",
        default_value_t = 1,
        allow_hyphen_values = true,
        help_heading = "Run control"
    )]
    pub start: i64,

    /// Write at most N passing lines
    #[arg(
        long = "limit",
        value_name = "N",
        allow_hyphen_values = true,
        help_heading = "Run control"
    )]
    pub limit: Option<i64>,

    /// Read at most N data lines
    #[arg(
        long = "test",
        value_name = "N",
        allow_hyphen_values = true,
        help_heading = "Run control"
    )]
    pub test: Option<i64>,

    /// Match `Consequence is <term>` against the term and every Sequence Ontology term below it, from the built-in table
    #[arg(long = "ontology", short = 'y', help_heading = "Ontology")]
    pub ontology: bool,

    /// Refused: Ensembl VEP's `filter_vep` connects to an Ensembl database for --ontology; the built-in table answers it here
    #[arg(long = "host", value_name = "HOST", help_heading = "Ontology")]
    pub host: Option<String>,

    /// Refused: the database user of --host
    #[arg(long = "user", value_name = "USER", help_heading = "Ontology")]
    pub user: Option<String>,

    /// Refused: the database password of --host
    #[arg(
        long = "pass",
        visible_alias = "password",
        value_name = "PASS",
        help_heading = "Ontology"
    )]
    pub pass: Option<String>,

    /// Refused: the database port of --host
    #[arg(long = "port", value_name = "PORT", help_heading = "Ontology")]
    pub port: Option<String>,

    /// Refused: the Ensembl database release --host connects to
    #[arg(long = "version", value_name = "N", help_heading = "Ontology")]
    pub version: Option<String>,

    /// Refused: an Ensembl registry file naming the database
    #[arg(long = "registry", value_name = "FILE", help_heading = "Ontology")]
    pub registry: Option<String>,
}

/// Runs `vep filter` over `argv`, the tokens after the `filter` word, and
/// returns the exit status: 0 on success, 1 after an `ERROR: [filter]` line,
/// 2 on a command-line syntax error.
pub fn run(argv: Vec<OsString>) -> i32 {
    let mut command = FilterArgs::command();
    command.build();
    if argv.is_empty() {
        // filter_vep prints its usage and exits 0 with no argument (106-109).
        return if command.print_long_help().is_ok() {
            0
        } else {
            1
        };
    }
    let argv: Vec<OsString> = argv
        .iter()
        .map(|word| OsString::from(filter_input::decode_word(word)))
        .collect();
    let spelled = argv::spell_as_declared(&argv, &command);
    let args = match FilterArgs::try_parse_from(
        std::iter::once(OsString::from("vep filter")).chain(spelled),
    ) {
        Ok(args) => args,
        Err(error) if error.kind() == ClapErrorKind::DisplayHelp => {
            let _ = error.print();
            return 0;
        }
        Err(error) => error.exit(),
    };
    match Filter::new(args).and_then(|filter| filter.run()) {
        Ok(()) => 0,
        Err(error) => {
            // The message may quote a key of the input, written as the bytes it was read from.
            let _ = write_line(&mut io::stderr().lock(), &error.to_string());
            1
        }
    }
}

/// A parsed command line, its filters compiled and its input and output named.
struct Filter {
    args: FilterArgs,
    filter_set: FilterSet,
}

impl Filter {
    /// Checks the command line as `configure` checks it (106-140): the
    /// database flags are refused, a filter or `--list` is required, and the
    /// filters are compiled, every `is_child` parent a term of the built-in
    /// table.
    fn new(args: FilterArgs) -> Result<Filter, CliError> {
        // The flags that name a database connection, refused by name (91-97).
        let database_flags = [
            ("host", args.host.is_some()),
            ("user", args.user.is_some()),
            ("pass", args.pass.is_some()),
            ("port", args.port.is_some()),
            ("version", args.version.is_some()),
            ("registry", args.registry.is_some()),
        ];
        if let Some((flag, _)) = database_flags.iter().find(|(_, given)| *given) {
            return Err(error(format!(
                "--{flag} is not supported: it names the Ensembl database connection filter_vep \
                 uses for --ontology; vep filter answers --ontology from its built-in table of \
                 Sequence Ontology consequence terms and connects to no database"
            )));
        }
        if args.filter.is_empty() && !args.list {
            return Err(error("No valid filters given"));
        }
        let filters: Vec<String> = args
            .filter
            .iter()
            .map(|filter| {
                if args.ontology {
                    ontology_rewrite(filter)
                } else {
                    filter.clone()
                }
            })
            .collect();
        let filter_set = FilterSet::deferred(&filters).map_err(|e| error(e.to_string()))?;
        if let Some(term) = filter_set
            .is_child_values()
            .into_iter()
            .find(|term| !so_ontology::is_term(term))
        {
            return Err(error(format!(
                "no SO term named {term} in the consequence ontology"
            )));
        }
        Ok(Filter { args, filter_set })
    }

    /// Opens the input (149-163), then the output (166-176), and filters.
    fn run(&self) -> Result<(), CliError> {
        let input: Box<dyn BufRead> = match &self.args.input_file {
            Some(path) => open_input(path, self.args.gz)?,
            None => Box::new(BufReader::new(io::stdin())),
        };
        let out_path = &self.args.output_file;
        let mut output: Box<dyn Write> = if out_path.to_ascii_lowercase().contains("stdout") {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            let out_file = filter_input::path_of(out_path);
            if out_file.exists() && !self.args.force_overwrite {
                return Err(error(format!(
                    "Output file {out_path} already exists. Specify a different output file with \
                     --output_file or overwrite existing file with --force_overwrite"
                )));
            }
            let file = File::create(&out_file)
                .map_err(|_| error(format!("Could not write to output file {out_path}")))?;
            Box::new(BufWriter::new(file))
        };
        let outcome = self
            .filter_lines(input, &mut output)
            .and_then(|()| output.flush().map_err(RunError::Io));
        match outcome {
            Ok(()) => Ok(()),
            // A reader that closed the pipe has what it asked for.
            Err(RunError::Io(e)) if e.kind() == ErrorKind::BrokenPipe => Ok(()),
            Err(RunError::Io(e)) => Err(error(e.to_string())),
            Err(RunError::Refused(e)) => Err(e),
        }
    }

    /// The loop of `main` (178-348) over the input's lines.
    fn filter_lines(
        &self,
        mut input: Box<dyn BufRead>,
        out: &mut dyn Write,
    ) -> Result<(), RunError> {
        let args = &self.args;
        let start = if args.soft_filter { 0 } else { args.start };
        // The count of the last passing line written (322); an extreme
        // `--limit` saturates instead of wrapping into an early stop.
        let last = args
            .limit
            .unwrap_or(NO_LIMIT)
            .saturating_add(start)
            .saturating_sub(1);
        let mut raw_headers: Vec<String> = Vec::new();
        let mut headers: Option<Headers> = None;
        // Settled on the first data line, not before the loop: an unknown
        // `--format` is refused only there, after the header lines are written
        // and `--list` is answered.
        let mut format: Option<Format> = None;
        let mut data_lines: i64 = 0;
        let mut count: i64 = 0;
        let mut lines_without_field: u64 = 0;
        let mut bytes = Vec::new();
        loop {
            bytes.clear();
            if input.read_until(b'\n', &mut bytes)? == 0 {
                break;
            }
            if bytes.last() == Some(&b'\n') {
                bytes.pop();
            }
            let decoded = filter_input::decode_line(&bytes);
            let line: &str = &decoded;
            if line.starts_with('#') {
                raw_headers.push(line.to_string());
                continue;
            }
            data_lines += 1;
            if args.test.is_some_and(|test| test != 0 && data_lines > test) {
                break;
            }
            if headers.is_none() {
                if raw_headers.is_empty() {
                    return Err(RunError::Refused(error("No headers found in input file")));
                }
                if args.soft_filter {
                    let column_line = raw_headers.pop().unwrap_or_default();
                    raw_headers.push(
                        "##FILTER=<ID=filter_vep_pass,Description=\"Variant passes filter_vep\">"
                            .to_string(),
                    );
                    raw_headers.push(
                        "##FILTER=<ID=filter_vep_fail,Description=\"Variant fails filter_vep\">"
                            .to_string(),
                    );
                    raw_headers.push(column_line);
                }
                if !args.count && !args.list {
                    write_lines(out, &raw_headers)?;
                }
                let parsed = Headers::parse(std::mem::take(&mut raw_headers), &args.vcf_info_field);
                if args.list {
                    print_list(&parsed)?;
                    return Ok(());
                }
                headers = Some(parsed);
            }
            let header_block = headers.as_ref().expect("parsed above");
            let line_format = match format {
                Some(line_format) => line_format,
                None => {
                    let first = match args.format.as_deref() {
                        None => filter_input::detect_format(line),
                        Some("vcf") => Format::Vcf,
                        Some("tab") => Format::Tab,
                        Some(other) => {
                            return Err(RunError::Refused(error(format!(
                                "Unable to parse data in format {other}"
                            ))));
                        }
                    };
                    format = Some(first);
                    first
                }
            };
            if line_format != Format::Vcf {
                if args.only_matched {
                    return Err(RunError::Refused(error(
                        "--only_matched is compatible only with VCF files",
                    )));
                }
                if args.soft_filter {
                    return Err(RunError::Refused(error(
                        "--soft_filter is compatible only with VCF files",
                    )));
                }
            }
            // The failed search is limited from the first line on (244, 273),
            // except on a tab line whose `Extra` keys differ line by line.
            self.filter_set.limit_synonym_search(
                line_format == Format::Vcf
                    || !header_block.fields().iter().any(|field| field == "Extra"),
            );
            let parsed: ParsedLine<'_> = match line_format {
                Format::Tab => filter_input::parse_tab_line(line, header_block),
                Format::Vcf => {
                    let parsed =
                        filter_input::parse_vcf_line(line, header_block, &args.vcf_info_field);
                    if parsed.lacks_consequence_field() {
                        lines_without_field += 1;
                    }
                    parsed
                }
            };

            let mut verdicts: Vec<bool> = Vec::new();
            for row in &parsed.rows() {
                let pass = self
                    .filter_set
                    .try_evaluate(row)
                    .map_err(|e| RunError::Refused(error(e.to_string())))?;
                verdicts.push(pass);
            }
            let line_pass = verdicts.contains(&true);

            let mut written: Cow<'_, str> = Cow::Borrowed(line);
            if args.only_matched {
                if let Some(rewritten) = parsed.only_matched(line, &verdicts) {
                    written = Cow::Owned(rewritten);
                }
            }
            if line_pass {
                count += 1;
            }
            if count < start {
                continue;
            }
            if args.soft_filter {
                write_line(out, &soft_filtered(&written, line_pass))?;
            } else if line_pass && !args.count {
                write_line(out, &written)?;
            }
            if count >= last {
                break;
            }
        }

        // An input without a data line writes its header lines (326-331), or
        // lists its fields (333-341).
        if data_lines == 0 {
            if !args.count && !args.list {
                write_lines(out, &raw_headers)?;
            }
            if args.list {
                print_list(&Headers::parse(raw_headers, &args.vcf_info_field))?;
                return Ok(());
            }
        }
        if args.count {
            writeln!(out, "{count}")?;
        }
        if lines_without_field > 0 && !args.no_missing_consequence_warning {
            out.flush()?;
            eprintln!(
                "WARNING: the VEP annotations INFO field {} was not found in {} line(s) of the input file",
                args.vcf_info_field, lines_without_field
            );
        }
        Ok(())
    }
}

/// A failure of the run: an I/O error, or a refusal with its message.
enum RunError {
    Io(io::Error),
    Refused(CliError),
}

impl From<io::Error> for RunError {
    fn from(e: io::Error) -> Self {
        RunError::Io(e)
    }
}

/// One line of input text and a newline, written as the bytes the text was
/// read from ([`filter_input::encode_line`]).
fn write_line(out: &mut dyn Write, text: &str) -> io::Result<()> {
    out.write_all(&filter_input::encode_line(text))?;
    out.write_all(b"\n")
}

/// The header lines, each on its own line; no line at all is one empty line,
/// as `join` and the newline after it write for an empty list (210-211, 329-330).
fn write_lines(out: &mut dyn Write, lines: &[String]) -> io::Result<()> {
    if lines.is_empty() {
        return out.write_all(b"\n");
    }
    for line in lines {
        write_line(out, line)?;
    }
    Ok(())
}

/// `--list`'s report on standard output (219-223), whatever `--output_file`
/// names: the field names in the order of their bytes.
fn print_list(headers: &Headers) -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = stdout.lock();
    writeln!(out, "Available fields:\n")?;
    for field in headers.listed() {
        write_line(&mut out, field)?;
    }
    out.flush()
}

/// The line under `--soft_filter` (305-317): the FILTER column replaced by the
/// verdict when it is `.`, else the verdict appended to it after `;`.
fn soft_filtered(line: &str, pass: bool) -> String {
    let verdict = if pass {
        "filter_vep_pass"
    } else {
        "filter_vep_fail"
    };
    // Perl's `split` without a limit drops the trailing empty columns.
    let mut columns: Vec<&str> = line.split('\t').collect();
    while columns.last() == Some(&"") {
        columns.pop();
    }
    let column = |i: usize| columns.get(i).copied().unwrap_or("");
    let filter = match column(6) {
        "." => verdict.to_string(),
        other => format!("{other};{verdict}"),
    };
    let mut out: Vec<&str> = (0..6).map(column).collect();
    out.push(&filter);
    out.extend(columns.iter().skip(7).copied());
    out.join("\t")
}

/// The input file (149-159): a `.gz` name, `--gz` or a gzip stream opens a
/// decompressing reader (BGZF included, a multi-member gzip stream).
fn open_input(path: &str, gz: bool) -> Result<Box<dyn BufRead>, CliError> {
    let file_path = filter_input::path_of(path);
    if !file_path.exists() {
        return Err(error(format!("Could not find input file {path}")));
    }
    let mut file = BufReader::new(
        File::open(&file_path)
            .map_err(|_| error(format!("Could not read from input file {path}")))?,
    );
    let magic = file.fill_buf().map_err(|e| error(e.to_string()))?;
    let gzip = gz || path.ends_with(".gz") || magic.starts_with(&[0x1f, 0x8b]);
    Ok(if gzip {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(file)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spelled(argv: &[&str]) -> Vec<String> {
        let mut command = FilterArgs::command();
        command.build();
        let argv: Vec<OsString> = argv.iter().map(OsString::from).collect();
        argv::spell_as_declared(&argv, &command)
            .into_iter()
            .map(|t| t.to_string_lossy().into_owned())
            .collect()
    }

    fn parsed(argv: &[&str]) -> FilterArgs {
        let spelled = spelled(argv).into_iter().map(OsString::from);
        FilterArgs::try_parse_from(std::iter::once(OsString::from("vep filter")).chain(spelled))
            .unwrap()
    }

    /// A repeated option is read as Getopt::Long reads `filter_vep`'s
    /// (72-103): a scalar or switch keeps its last occurrence
    /// (`Getopt/Long.pm:677` assigns each one), and `--filter` keeps every
    /// occurrence (102).
    #[test]
    fn repeated_options_keep_the_last_value_and_filters_accumulate() {
        let args = parsed(&[
            "-s",
            "1",
            "-s",
            "2",
            "--count",
            "--count",
            "--limit",
            "2",
            "--limit",
            "1",
            "-f",
            "IMPACT is HIGH",
            "-f",
            "SYMBOL is TPTE",
        ]);
        assert_eq!(args.start, 2);
        assert!(args.count);
        assert_eq!(args.limit, Some(1));
        assert_eq!(args.filter, ["IMPACT is HIGH", "SYMBOL is TPTE"]);
        let args = parsed(&["--start", "3", "--start=4", "-o", "a.vcf", "-o", "b.vcf"]);
        assert_eq!(args.start, 4);
        assert_eq!(args.output_file, "b.vcf");
    }

    /// Getopt::Long's spellings reach the parser's names: one dash, any case,
    /// a unique prefix, an attached value; the token after a value option is
    /// its value whatever it looks like; an unknown name is left as written.
    #[test]
    fn options_are_spelled_as_getopt_long_reads_them() {
        assert_eq!(
            spelled(&[
                "-i",
                "in.vcf",
                "-filter",
                "SYMBOL is -",
                "--COUNT",
                "--only",
                "-f=x"
            ]),
            [
                "--input_file",
                "in.vcf",
                "--filter",
                "SYMBOL is -",
                "--count",
                "--only_matched",
                "--filter=x"
            ]
        );
        assert_eq!(
            spelled(&[
                "--soft",
                "--no_missing",
                "-y",
                "-s",
                "2",
                "--lim",
                "3",
                "--password",
                "x"
            ]),
            [
                "--soft_filter",
                "--no_missing_consequence_warning",
                "--ontology",
                "--start",
                "2",
                "--limit",
                "3",
                "--pass",
                "x"
            ]
        );
        assert_eq!(spelled(&["--nope", "--", "-i"]), ["--nope", "--", "-i"]);
    }

    /// `--soft_filter` replaces a `.` FILTER and appends to any other (309-312).
    #[test]
    fn soft_filter_marks_the_filter_column() {
        assert_eq!(
            soft_filtered("21\t10\tid\tG\tA\t50\t.\tDP=1", true),
            "21\t10\tid\tG\tA\t50\tfilter_vep_pass\tDP=1"
        );
        assert_eq!(
            soft_filtered("21\t10\tid\tG\tA\t50\tPASS\tDP=1\tGT\t0/1", false),
            "21\t10\tid\tG\tA\t50\tPASS;filter_vep_fail\tDP=1\tGT\t0/1"
        );
    }

    /// `--ontology` rewrites `<con...> is` in any case to `Consequence
    /// is_child`, an `is_child` already written included (115-118). The
    /// substitution is `s/$&/.../g` over the first match's text, so a second
    /// unit spelled in another case stands, as do a unit with two blanks and
    /// one whose blank is a UTF-8 no-break space; a tab is a blank.
    #[test]
    fn ontology_rewrites_consequence_is() {
        assert_eq!(
            ontology_rewrite("Consequence is coding_sequence_variant"),
            "Consequence is_child coding_sequence_variant"
        );
        assert_eq!(
            ontology_rewrite("consequence is_child intron_variant and IMPACT is HIGH"),
            "Consequence is_child intron_variant and IMPACT is HIGH"
        );
        assert_eq!(ontology_rewrite("SYMBOL is BRCA2"), "SYMBOL is BRCA2");
        assert_eq!(
            ontology_rewrite("Consequence is X or consequence is Y or Consequence is Z"),
            "Consequence is_child X or consequence is Y or Consequence is_child Z"
        );
        assert_eq!(
            ontology_rewrite("consequence is_child A and CONSEQUENCE IS B"),
            "Consequence is_child A and CONSEQUENCE IS B"
        );
        assert_eq!(ontology_rewrite("Consequence  is X"), "Consequence  is X");
        assert_eq!(
            ontology_rewrite("Consequence\u{a0}is X"),
            "Consequence\u{a0}is X"
        );
        assert_eq!(
            ontology_rewrite("Consequence\tis X"),
            "Consequence is_child X"
        );
    }
}
