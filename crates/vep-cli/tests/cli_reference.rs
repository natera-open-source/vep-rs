// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `docs/cli-reference.md` is rendered from the command-line definition, the
//! option-rule tables and the error table, so the page, `--help` and the parser
//! cannot disagree. The test asserts the committed page is the rendering;
//! `VEP_UPDATE_CLI_REFERENCE=1` rewrites it instead.

use std::any::TypeId;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use clap::{Arg, ArgAction, Command};
use vep_cli::args::Args;
use vep_cli::argv::{
    unresolved_perl_spellings, COMMAND_LINE_ONLY, DEFAULT_DIR, ENV_PREFIX, INI_FILE,
};
use vep_cli::config::{Treatment, COMPONENTS, ENVIRONMENT, ERROR_TABLE};
use vep_cli::option_rules::{
    ALLOW_MULTIPLE, DEPRECATED, INCOMPATIBLE, INCOMPATIBLE_LIFTED, LIST_FLAGS, OPTION_SETS,
    OWN_RULES, REQUIRES, UNSORTABLE, VALID,
};

fn reference_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/cli-reference.md")
}

/// A Markdown table cell: pipes escaped, newlines and runs of blanks collapsed.
fn cell(text: &str) -> String {
    text.split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .replace('|', "\\|")
}

fn code(text: &str) -> String {
    format!("`{text}`")
}

/// A long flag name in code: `--name`.
fn flag_in_code(name: &str) -> String {
    code(&format!("--{name}"))
}

fn flags_in_code(flags: &[&str]) -> String {
    flags
        .iter()
        .map(|f| flag_in_code(f))
        .collect::<Vec<_>>()
        .join(", ")
}

/// `--long`, then `-s`, then the visible aliases, each in code.
fn flag_names(arg: &Arg) -> String {
    let mut names = Vec::new();
    if let Some(long) = arg.get_long() {
        names.push(flag_in_code(long));
    }
    if let Some(short) = arg.get_short() {
        names.push(code(&format!("-{short}")));
    }
    for alias in arg.get_visible_aliases().into_iter().flatten() {
        names.push(flag_in_code(alias));
    }
    for alias in arg.get_visible_short_aliases().into_iter().flatten() {
        names.push(code(&format!("-{alias}")));
    }
    let mut out = names.join(", ");
    if arg.is_hide_set() {
        out.push_str(" (hidden)");
    }
    out
}

/// The value type as the reference names it: `bool` for a switch, the possible
/// values when the flag is an enumeration, the Rust type otherwise, with
/// `repeatable` for a flag that accumulates and `comma-separated` for a list.
fn flag_type(arg: &Arg) -> String {
    if matches!(arg.get_action(), ArgAction::SetTrue) {
        return "bool".to_string();
    }
    let possible = arg.get_possible_values();
    let base = if !possible.is_empty() {
        possible
            .iter()
            .map(|p| code(p.get_name()))
            .collect::<Vec<_>>()
            .join(", ")
    } else {
        let id = arg.get_value_parser().type_id();
        let name = if id == TypeId::of::<String>() {
            "String"
        } else if id == TypeId::of::<usize>() {
            "usize"
        } else if id == TypeId::of::<u32>() {
            "u32"
        } else if id == TypeId::of::<u64>() {
            "u64"
        } else if id == TypeId::of::<i64>() {
            "i64"
        } else {
            panic!(
                "--{}: a value type the reference does not name",
                arg.get_id()
            )
        };
        name.to_string()
    };
    let qualifier = if arg.get_value_delimiter().is_some() {
        Some("comma-separated")
    } else if matches!(arg.get_action(), ArgAction::Append) {
        Some("repeatable")
    } else {
        None
    };
    match qualifier {
        Some(qualifier) => format!("{base} ({qualifier})"),
        None => base,
    }
}

fn flag_default(arg: &Arg) -> String {
    if matches!(arg.get_action(), ArgAction::SetTrue) {
        return "false".to_string();
    }
    let defaults: Vec<String> = arg
        .get_default_values()
        .iter()
        .map(|v| code(&v.to_string_lossy()))
        .collect();
    if defaults.is_empty() {
        "--".to_string()
    } else {
        defaults.join(", ")
    }
}

fn flag_help(arg: &Arg) -> String {
    arg.get_long_help()
        .or(arg.get_help())
        .map(|h| h.to_string())
        .unwrap_or_default()
}

fn table(out: &mut String, header: &[&str], rows: &[Vec<String>]) {
    writeln!(out, "| {} |", header.join(" | ")).unwrap();
    writeln!(
        out,
        "|{}",
        header.iter().map(|_| " --- |").collect::<String>()
    )
    .unwrap();
    for row in rows {
        writeln!(out, "| {} |", row.join(" | ")).unwrap();
    }
    out.push('\n');
}

/// The page.
fn render(command: &Command) -> String {
    let mut out = String::new();
    out.push_str("# CLI Reference\n\n");
    out.push_str(
        "The `vep` binary accepts the Ensembl VEP command line: every flag is Ensembl VEP's, does \
         what Ensembl VEP's does, or is refused with an error naming what it would do, and the \
         few flags vep-rs adds are marked as its own. A flag the error table below governs \
         (refused, or accepted without effect) carries that table's text as its help.\n\n",
    );
    out.push_str(
        "This page is generated from the command-line definition and the option tables by the \
         `cli_reference` test of `crates/vep-cli`; edit the sources, then run \
         `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference` to rewrite \
         it.\n\n",
    );
    out.push_str("```bash\n# Show all flags\nvep --help\n\n# Short form\nvep -i input.vcf -o output.txt --offline --json_cache /path/to/json_cache\n```\n\n");
    out.push_str(
        "Worked examples are in [examples.md](examples.md); the output formats and their parity \
         with Ensembl VEP in [output-formats.md](output-formats.md).\n\n",
    );

    out.push_str("## Flags\n\n");
    out.push_str(
        "Grouped by section as `--help` groups them. A repeated single-value flag keeps its last \
         value, as Getopt::Long does; `-o` accumulates and the last wins with a warning, as the \
         `vep` script does (vep:32-33). A long option is also read in any case, as a unique \
         prefix and after one dash, and from the `--config` file, `vep.ini` and the `VEP_*` \
         variables (the configuration sources below).\n\n",
    );
    let mut sections: Vec<(String, Vec<Vec<String>>)> = Vec::new();
    for arg in command.get_arguments() {
        let Some(long) = arg.get_long() else {
            continue;
        };
        if long == "help" || long == "version" {
            continue;
        }
        let heading = arg
            .get_help_heading()
            .unwrap_or_else(|| panic!("--{long} has no help heading"))
            .to_string();
        let row = vec![
            flag_names(arg),
            flag_type(arg),
            flag_default(arg),
            cell(&flag_help(arg)),
        ];
        match sections.iter_mut().find(|(h, _)| *h == heading) {
            Some((_, rows)) => rows.push(row),
            None => sections.push((heading, vec![row])),
        }
    }
    for (heading, rows) in &sections {
        writeln!(out, "### {heading}\n").unwrap();
        table(&mut out, &["Flag", "Type", "Default", "Description"], rows);
    }
    out.push_str(
        "`--help` prints this list and exits 0, as does a command line with no argument at all \
         (vep:29); `--version` prints `vep <version>` and exits 0.\n\n",
    );

    out.push_str("## Option rules\n\n");
    out.push_str(
        "Ensembl VEP's option rules are applied as data, from `Config.pm:327-645` at \
         release/116.2, in this order after parsing: a flag the error table refuses is refused \
         first; then the option sets switch on the flags they imply, in the order below, a set's \
         value overwriting an explicit one (`apply_option_sets`, `Config.pm:816-884`); then the \
         checks of `check_config` (`Config.pm:944-985`): valid values, incompatible pairs, \
         requirements, deprecations, each refused with Ensembl VEP's own message. A rule \
         naming a flag vep-rs does not declare cannot fire and is listed so the table reads as \
         Ensembl VEP's.\n\n",
    );
    out.push_str("### Option sets (`Config.pm:351-593`)\n\n");
    let rows: Vec<Vec<String>> = OPTION_SETS
        .iter()
        .map(|set| {
            vec![
                flags_in_code(set.flags),
                set.set
                    .iter()
                    .map(|(k, v)| {
                        if *v == "1" {
                            flag_in_code(k)
                        } else {
                            code(&format!("--{k} {v}"))
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(", "),
            ]
        })
        .collect();
    table(&mut out, &["Flag", "Switches on"], &rows);
    out.push_str(
        "`--everything` therefore reaches `--hgvsc` and `--hgvsp` through `--hgvs`, \
         `--mane_select` and `--mane_plus_clinical` through `--mane`, and `--check_existing` \
         through the frequency flags; it does not switch on `--allele_number` or `--hgvsg`. \
         `--regulatory` is switched on and annotates nothing, so every column of Ensembl VEP's \
         `--everything` is produced except the regulatory columns, and no RegulatoryFeature or \
         MotifFeature row is written.\n\n",
    );
    out.push_str("### Incompatible flags (`Config.pm:615-632`)\n\n");
    out.push_str(
        "`Can't use --<flag> and --<other> together`, naming the first pair in this table's \
         order (Ensembl VEP iterates a hash there, so which pair it names varies from run to \
         run); Ensembl VEP skips the check under `--safe`.\n\n",
    );
    let rows: Vec<Vec<String>> = INCOMPATIBLE
        .iter()
        .map(|(flag, others)| vec![flag_in_code(flag), flags_in_code(others)])
        .collect();
    table(&mut out, &["Flag", "Cannot be combined with"], &rows);
    for pair in INCOMPATIBLE_LIFTED {
        writeln!(
            out,
            "{} beside {} is accepted here: {}.\n",
            flag_in_code(pair.flag),
            flag_in_code(pair.other),
            cell(pair.reason)
        )
        .unwrap();
    }
    out.push_str("### Required flags (`Config.pm:607-612`)\n\n");
    out.push_str("`You must set --<required> to use --<flag>`.\n\n");
    let rows: Vec<Vec<String>> = REQUIRES
        .iter()
        .map(|(flag, required)| vec![flag_in_code(flag), flags_in_code(required)])
        .collect();
    table(&mut out, &["Flag", "Needs"], &rows);
    out.push_str("### Valid values (`Config.pm:596-604`)\n\n");
    out.push_str(
        "A declared flag carries its list as a value enumeration, so a value outside it is an \
         argument syntax error (exit 2) listing the possible values.\n\n",
    );
    let rows: Vec<Vec<String>> = VALID
        .iter()
        .map(|(flag, values)| {
            vec![
                flag_in_code(flag),
                values
                    .iter()
                    .map(|v| code(v))
                    .collect::<Vec<_>>()
                    .join(", "),
            ]
        })
        .collect();
    table(&mut out, &["Flag", "Values"], &rows);
    out.push_str("### Deprecated flags (`Config.pm:636-638`)\n\n");
    for (flag, replacement) in DEPRECATED {
        match replacement {
            Some(new) => writeln!(
                out,
                "- `--{flag}`: `--{flag} has been deprecated - please use --{new} instead`"
            )
            .unwrap(),
            None => writeln!(out, "- `--{flag}`: `--{flag} has been deprecated`").unwrap(),
        }
    }
    out.push('\n');
    out.push_str("### Lists (`Config.pm:327-346`)\n\n");
    writeln!(
        out,
        "- Comma-separated values: {}.",
        flags_in_code(LIST_FLAGS)
    )
    .unwrap();
    writeln!(
        out,
        "- Repeatable, each occurrence appending: {}.",
        flags_in_code(ALLOW_MULTIPLE)
    )
    .unwrap();
    writeln!(
        out,
        "- Input formats Ensembl VEP exempts from its sorted-input check (`Config.pm:640-645`): \
         {}. vep-rs performs no such check in any format: every record is annotated in input \
         order, and an input whose chromosomes interleave keeps more transcript indexes in \
         memory.\n",
        UNSORTABLE
            .iter()
            .map(|f| code(f))
            .collect::<Vec<_>>()
            .join(", ")
    )
    .unwrap();
    out.push_str("### Rules of vep-rs's own\n\n");
    for rule in OWN_RULES {
        writeln!(out, "- `ERROR: [{}] {}`", rule.component, rule.message).unwrap();
    }
    out.push('\n');
    out.push_str("### Where `Config.pm` is not followed\n\n");
    out.push_str(
        "- An output file matching `stdout` in any case does not force `--quiet` or drop \
         `--verbose` (`Config.pm:933-937`): the status lines go to standard error, so the \
         annotation on standard output is clean either way, and `--quiet --verbose` stays \
         refused with `-o STDOUT` as without.\n\
         - A run naming none of `--database`, `--cache`, `--offline` and `--custom` is not \
         refused (`Config.pm:993-1017`): vep-rs reads a local cache in every run, so those \
         flags select nothing and the run proceeds on the cache it is given.\n\
         - A `%name%` in an option-set value is replaced even when the flag's value is `0` or \
         empty (`Config.pm:867` leaves the reference in place then); no declared flag carries \
         such a value.\n\
         - A VCF output selected with `--output_format vcf` fills `SYMBOL`, `BIOTYPE`, `EXON` and \
         `INTRON` as `--vcf` does; Ensembl VEP's `vcf` option set applies to the `--vcf` flag \
         alone and leaves the CSQ fields empty under `--output_format vcf`.\n\
         - A value from the `--config` file, `vep.ini` or the environment never overwrites one a \
         higher source gives, even a `0`, an empty value or a bare `--cache`: Ensembl VEP's \
         `||=` (`Config.pm:1101, 1158`) lets a lower source overwrite a value Perl reads as \
         false, so `--shift_hgvs 0` on its command line yields to a file's `shift_hgvs 1`.\n\
         - A value from a file or the environment that the parser refuses (`fork many`, an empty \
         `format`) is an argument syntax error (exit 2); Ensembl VEP's `is_valid_param` drops the \
         line (`Config.pm:1046-1047`).\n\
         - A `--config` path that is a directory, or cannot be read for any reason, is the `config` \
         error; Ensembl VEP's `open` succeeds on a directory and reads nothing from it. A \
         `vep.ini` that is a directory is read as nothing here too. A byte outside UTF-8 in a \
         value becomes the replacement character, where Ensembl VEP passes the byte.\n\
         - A flag read from a file or the environment acts exactly as if typed, refusals included; \
         the status lines Ensembl VEP prints under `--verbose` while reading the sources (`Read \
         configuration from <file>`, `Set '<flag>=<value>' from environment variable <VAR>`, \
         `Config.pm:1107, 1145-1164`) are not printed.\n\
         - A leading `~` in `--config`, `--dir` and a file's `dir` is expanded to the home \
         directory, which Ensembl VEP leaves to the shell.\n\
         - `-V` is the version, where Getopt::Long's `ignore_case` reads it as `-v`, and an \
         upper-case short names no flag here; a short with its value attached (`-ofile.txt`) and a \
         cluster of switches (`-qv`) are clap's readings, where Getopt::Long refuses the token; a \
         `+` before an option name, which Getopt::Long accepts, is not an option here.\n\n",
    );

    out.push_str("## Flags refused or accepted without effect\n\n");
    out.push_str(
        "The error table: every flag vep-rs refuses, and the flags accepted without effect \
         whose only documentation is this table, with what Ensembl VEP does with each (cited to \
         its source at release/116.2). A refused flag exits 1 with the message below under the \
         `cli` component; an accepted one parses and changes nothing. A flag whose description \
         above reads `Accepted for VEP compatibility` parses as well, its description saying \
         what, if anything, it changes. The help text of each entry is this table's.\n\n",
    );
    let rows: Vec<Vec<String>> = ERROR_TABLE
        .iter()
        .map(|entry| {
            let (treatment, text) = match entry.treatment {
                Treatment::Refused { component, message } => (
                    format!("refused (`{component}`)"),
                    code(&format!("ERROR: [{component}] {message}")),
                ),
                Treatment::AcceptedWithoutEffect { reason } => {
                    ("accepted without effect".to_string(), cell(reason))
                }
            };
            let flag = if entry.declared {
                flag_in_code(entry.flag)
            } else {
                format!("`--{}` (not declared; caught before parsing)", entry.flag)
            };
            vec![flag, treatment, text, cell(entry.ensembl_vep)]
        })
        .collect();
    table(
        &mut out,
        &["Flag", "Treatment", "Message or reason", "Ensembl VEP"],
        &rows,
    );

    out.push_str("## Configuration sources\n\n");
    writeln!(
        out,
        "A flag is read from four sources, the first that names it winning, as Ensembl VEP reads \
         them (`Config.pm:693-708`): the command line; the `--config` file; `<dir>/{INI_FILE}`, \
         where `<dir>` is `--dir` (the command line's, else the `--config` file's) and \
         `{DEFAULT_DIR}` by default; and the environment variable `{ENV_PREFIX}<FLAG>` \
         (`{ENV_PREFIX}DIR_PLUGINS` sets `--dir_plugins`). A repeatable flag ({}) collects every \
         source's values, the command line's first. No source behind the command line can set \
         {}.\n",
        flags_in_code(ALLOW_MULTIPLE),
        flags_in_code(COMMAND_LINE_ONLY)
    )
    .unwrap();
    writeln!(
        out,
        "### The `--config` file and `{INI_FILE}` (`Config.pm:1068-1110`)\n"
    )
    .unwrap();
    writeln!(
        out,
        "One flag per line: the flag's long name without its dashes (one leading dash is allowed), \
         a blank, then its value. A switch is on with the value `1` and off with `0`; `cache` is \
         on with `1`, a directory or an empty quoted value and off with `0`; a key alone on its \
         line is dropped; a flag taking a value takes the first word after its name, or, for a \
         repeatable flag, the rest of the line; blanks inside double quotes are part of a value \
         and the quotes are removed. A line starting with `#`, or blank, is skipped. A key is the \
         long name exactly: a prefix, another case or a name that is no option of Ensembl VEP's \
         or of this version names nothing and the line is dropped, as Ensembl VEP drops it, and \
         so does an alias, unless the alias is an option of Ensembl VEP's own (`full_cache_dir` \
         names `--json_cache`); a key of an Ensembl VEP option this version lacks (`stats_html 1`, \
         `clin_sig_allele 0`) is refused as the flag on the command line is (exit 2), unless it \
         is a switch of Ensembl VEP's set to `0`, which is off there and adds nothing here. Of \
         several lines naming one flag, the first whose value is not `0` or empty counts. A \
         `--config` file that cannot be opened is `ERROR: [config] Could not open config \
         file \"<file>\"`; `{INI_FILE}` is read only when it is a file.\n\n\
         ```\n# flags of a run\ncache 1\noffline 1\ndir_cache /data/vep\nspecies homo_sapiens\n\
         assembly GRCh38\nplugin LoFtool,/data/LoFtool_scores.txt\n```\n"
    )
    .unwrap();
    out.push_str("### Environment variables (`Config.pm:1127-1168`)\n\n");
    writeln!(
        out,
        "`{ENV_PREFIX}<FLAG>=<value>` sets `--<flag>`, the name lower-cased, when the value is \
         not empty; the name is read as a configuration file's key is: a name that is no option \
         of Ensembl VEP's or of this version is ignored, and one of an Ensembl VEP option this \
         version lacks is refused as the flag is (exit 2). The variables are read in name \
         order.\n"
    )
    .unwrap();
    out.push_str(
        "Environment variables vep-rs reads outside this mechanism, which Ensembl VEP has no \
         counterpart for:\n\n",
    );
    for (name, help) in ENVIRONMENT {
        writeln!(out, "- `{name}`: {help}.").unwrap();
    }
    out.push('\n');
    out.push_str("### Option spellings (`Config.pm:62`)\n\n");
    out.push_str(
        "Getopt::Long's defaults apply to the command line: a long option is read in any case \
         (`--FORMAT vcf`), abbreviated to any prefix that names one option (`--form vcf`), and \
         after one dash as well as two (`-cache`, `-form vcf`); a value follows its option or is \
         attached with `=`, and the token after an option that requires a value is that value \
         whatever it looks like (`--input_file -cache` names the file `-cache`). A prefix is \
         resolved against this version's flags first and then, where a flag vep-rs adds makes it \
         ambiguous, as Getopt::Long resolves it against Ensembl VEP's `@VEP_PARAMS` \
         (`Config.pm:73-276`): `--js` is `--json` though `--json_cache` also starts with it, `--pa` \
         `--password`, `--ve` `--verbose`. An Ensembl VEP option this version lacks is refused by \
         its full name in every spelling (`-stats_html` and `--stats_h` are refused as `--stats_html`, exit \
         2), never read as a cluster of shorts. A prefix that names one flag of this version \
         resolves to it even where Ensembl VEP's larger set leaves it ambiguous (`--sy` is \
         `--symbol`); a prefix of several flags that Getopt::Long resolves no better (`--fo`), or \
         a name no option has, is an argument syntax error (exit 2). A single-letter short is \
         clap's and is matched by case: `-v` is `--verbose`, `-V` the version, `-qv` the cluster \
         `-q -v`, and `-ofile.txt` the value `file.txt` of `-o`, a spelling Getopt::Long refuses. \
         Of the spellings Getopt::Long resolves to an option this version declares, these resolve \
         differently here:\n\n",
    );
    let unresolved = unresolved_perl_spellings(command);
    if unresolved.is_empty() {
        out.push_str("- none: every such spelling names the same option.\n\n");
    } else {
        let rows: Vec<Vec<String>> = unresolved
            .iter()
            .map(|u| {
                vec![
                    flag_in_code(&u.spelling),
                    flag_in_code(&u.perl),
                    if u.candidates.is_empty() {
                        "no flag".to_string()
                    } else {
                        flags_in_code(&u.candidates.iter().map(String::as_str).collect::<Vec<_>>())
                    },
                ]
            })
            .collect();
        table(
            &mut out,
            &[
                "Spelling",
                "Ensembl VEP reads",
                "vep-rs flags it could name",
            ],
            &rows,
        );
    }

    out.push_str("## Exit status and errors\n\n");
    out.push_str(
        "- `0`: the run completed; also `--help`, `--version` and an empty command line.\n\
         - `1`: a refusal or a failed run. One line `ERROR: [<component>] <message>` is written \
         to standard error, where `<message>` is Ensembl VEP's text where it has one, without \
         its own `ERROR: ` prefix. The components:\n",
    );
    for (component, meaning) in COMPONENTS {
        writeln!(out, "  - `{component}`: {meaning}").unwrap();
    }
    out.push_str(
        "- `2`: an argument syntax error, in clap's words with its `error:` prefix: an unknown \
         flag, a prefix of several flags, a missing value, a value outside an enumeration.\n\n",
    );
    out.push_str(
        "`--help` and `--version` win over a flag the error table refuses before parsing \
         wherever they stand on the command line; over a flag the parser does not know they win \
         only when they come first, since clap reads left to right (Ensembl VEP's Getopt::Long \
         runs with `pass_through`, Config.pm:62, so any unknown flag reaches the usage there, \
         vep:29).\n\n",
    );
    out.push_str(
        "A panic aborts the process with the panic's own status. Ensembl VEP's own exit codes \
         (255 on `die`, 2 on a failed open) are not copied.\n",
    );
    out
}

#[test]
fn cli_reference_is_current() {
    let mut command = Args::command_with_error_table();
    command.build();
    let rendered = render(&command);
    let path = reference_path();
    if std::env::var_os("VEP_UPDATE_CLI_REFERENCE").is_some() {
        fs::write(&path, &rendered).unwrap();
        return;
    }
    let committed = fs::read_to_string(&path).unwrap_or_default();
    if committed != rendered {
        let first = committed
            .lines()
            .zip(rendered.lines())
            .position(|(a, b)| a != b)
            .unwrap_or(committed.lines().count().min(rendered.lines().count()));
        panic!(
            "docs/cli-reference.md is not the rendering of the command line and the tables: first \
             difference at line {}\n  committed: {:?}\n  rendered:  {:?}\nrun VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference",
            first + 1,
            committed.lines().nth(first).unwrap_or(""),
            rendered.lines().nth(first).unwrap_or("")
        );
    }
}

/// Every table-backed flag's help is the table's text, verbatim, so `--help`
/// and the page cannot drift from the refusal.
#[test]
fn help_of_table_backed_flags_is_the_tables() {
    let mut command = Args::command_with_error_table();
    command.build();
    let help: BTreeMap<&str, String> = command
        .get_arguments()
        .filter_map(|a| a.get_long().map(|l| (l, flag_help(a))))
        .collect();
    for entry in ERROR_TABLE.iter().filter(|e| e.declared) {
        assert_eq!(
            help.get(entry.flag),
            Some(&entry.help()),
            "--{}",
            entry.flag
        );
    }
}
