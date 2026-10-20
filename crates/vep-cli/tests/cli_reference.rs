// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The CLI reference pages (`docs/cli-reference.md` and the group pages under
//! `docs/cli/`) are rendered from the command-line definition, the option-rule
//! tables and the error table, so the pages, `--help` and the parser cannot
//! disagree. The test asserts every committed page is its rendering;
//! `VEP_UPDATE_CLI_REFERENCE=1` rewrites them instead.

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

const LANDING: &str = "docs/cli-reference.md";

const COMPATIBILITY: &str = "docs/cli/compatibility.md";

const REGENERATE: &str = "VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference";

/// One group page: its file under `docs/`, its title, the landing table's
/// summary of it, its hand-written introduction, and the `--help` headings it
/// carries, each as a section of the page.
struct Group {
    file: &'static str,
    title: &'static str,
    covers: &'static str,
    intro: &'static str,
    sections: &'static [(&'static str, &'static str)],
}

const GROUPS: &[Group] = &[
    Group {
        file: "docs/cli/input-output.md",
        title: "Input and output",
        covers: "The input file, its format and how it is read, the output file, its format and \
                 compression, the header lines, the warnings and skipped-variants files.",
        intro: "These flags name the input file and how it is read, and the output file, its \
                format, its compression and its header lines. The formats themselves are described \
                under [Output formats](../output-formats.md).\n\n",
        sections: &[("Input", "Input"), ("Output", "Output")],
    },
    Group {
        file: "docs/cli/cache-species-runtime.md",
        title: "Cache, species and runtime",
        covers: "The cache directory and its layout, the reference FASTA, species and assembly, \
                 the batch size, thread counts and distances.",
        intro: "These flags name the annotation cache and the reference FASTA, the species and \
                assembly, and the batch size, thread counts and distances of a run. The cache is \
                found as Ensembl VEP finds it, under `--dir_cache`, `--species` and \
                `--cache_version`, or named directly with `--json_cache`. Building or converting \
                a cache is described under [Cache setup](../cache-setup.md).\n\n",
        sections: &[
            ("Cache / Offline", "Cache and offline"),
            ("Species / Assembly", "Species and assembly"),
            ("Runtime", "Runtime"),
        ],
    },
    Group {
        file: "docs/cli/output-control.md",
        title: "Output control",
        covers: "The field flags, `--everything` among them, that choose which annotation \
                 columns appear.",
        intro: "These flags choose which annotation fields appear in the output. `--everything` \
                switches on exactly Ensembl VEP's `--everything` set, listed under \
                [option sets](../cli-reference.md#option-sets-configpm351-593).\n\n",
        sections: &[("Output Control", "Flags")],
    },
    Group {
        file: "docs/cli/allele-frequencies.md",
        title: "Allele frequencies and co-located variants",
        covers: "Co-located known variants, population frequencies, the frequency filter, PubMed \
                 identifiers and gene-phenotype associations.",
        intro: "These flags report the population allele frequencies, PubMed identifiers and \
                gene-phenotype associations of co-located known variants, and filter the input by \
                frequency. `--check_existing`, which reports the co-located variants themselves, \
                is listed under [Output control](output-control.md); the frequency flags \
                `--everything` switches on are listed under \
                [option sets](../cli-reference.md#option-sets-configpm351-593).\n\n",
        sections: &[("Allele Frequencies", "Flags")],
    },
    Group {
        file: "docs/cli/filtering-picking.md",
        title: "Filtering, picking and transcript sets",
        covers: "The pick and flag-pick families, the allele columns, the transcript filter and \
                 the transcript-set flags.",
        intro: "These flags reduce the rows reported per record, allele or gene, flag the ones a \
                pick would drop, shape the allele columns and filter transcripts by their \
                attributes. The transcript-set flags select the cache read and the GENCODE set \
                annotated.\n\n",
        sections: &[
            ("Filtering / Picking", "Filtering and picking"),
            ("Transcript Sets", "Transcript sets"),
        ],
    },
    Group {
        file: "docs/cli/plugins.md",
        title: "Plugins",
        covers: "Loading built-in and dynamic-library plugins.",
        intro: "These flags load the built-in plugins and dynamic-library plugins, and name the \
                directory the latter are searched in. The built-in plugins, their parameters and \
                their data files are described under [Plugins](../plugins.md); the C ABI a \
                dynamic-library plugin implements is described under \
                [Writing a plugin](../developers/writing-a-plugin.md).\n\n",
        sections: &[("Plugins", "Flags")],
    },
    Group {
        file: "docs/cli/statistics-verbosity-misc.md",
        title: "Statistics, verbosity and miscellaneous",
        covers: "The statistics flags, logging, and the flags that change how a record is \
                 handled.",
        intro: "These flags cover the statistics file Ensembl VEP writes and vep-rs does not, the \
                logging level, and the flags that change how a record is handled.\n\n",
        sections: &[
            ("Statistics", "Statistics"),
            ("Verbosity", "Verbosity"),
            ("Miscellaneous", "Miscellaneous"),
        ],
    },
    Group {
        file: COMPATIBILITY,
        title: "Compatibility with Ensembl VEP",
        covers: "The database flags, the flags refused or accepted without effect, where \
                 `Config.pm` is not followed, the configuration sources, the option spellings \
                 and the exit status.",
        intro: "This page gathers the Ensembl VEP behaviour vep-rs refuses or accepts without \
                effect, and how the command line is read. It covers the database flags, the \
                error table, where `Config.pm` is not followed, the configuration sources and \
                their priority, the option spellings, and the exit status and error \
                contract.\n\n",
        sections: &[("Database", "Database flags")],
    },
];

/// The declared flags whose help reads `Accepted for VEP compatibility`, with
/// what Ensembl VEP does with each. The rendering checks the list against the
/// help texts both ways, so a flag that gains an implementation leaves this
/// table and a new compatibility flag enters it.
const COMPATIBILITY_FLAGS: &[(&str, &str)] = &[
    (
        "regulatory",
        "Reports overlaps with regulatory regions and transcription factor binding sites",
    ),
    (
        "custom",
        "Adds annotation from a tabix-indexed or bigWig file",
    ),
    ("no_stats", "Writes no statistics file"),
    ("stats_file", "Names the statistics file"),
    (
        "stats_text",
        "Writes the statistics as plain text instead of HTML",
    ),
];

const COMPATIBILITY_HELP: &str = "accepted for vep compatibility";

/// The repository path of a page, from this crate's directory.
fn repo_path(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(file)
}

/// A page's link target from the landing page: its path under `docs/`.
fn site_link(file: &str) -> &str {
    file.strip_prefix("docs/")
        .unwrap_or_else(|| panic!("{file} is not a page under docs/"))
}

/// A Markdown table cell: newlines and runs of blanks collapsed, a pipe escaped
/// everywhere, and outside a code span a bare `--flag` put in code (the site's
/// smart punctuation would turn its dashes into an en dash) and `<`, `>` and
/// `_` escaped, so that `<leaf>` is not read as an HTML tag and
/// `[_refseq|_merged]` does not open an emphasis.
fn cell(text: &str) -> String {
    let joined = text.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut out = String::with_capacity(joined.len() + 16);
    for (i, segment) in joined.split('`').enumerate() {
        if i > 0 {
            out.push('`');
        }
        let in_code_span = i % 2 == 1;
        if in_code_span {
            out.push_str(&segment.replace('|', "\\|"));
        } else {
            out.push_str(&prose_cell(segment));
        }
    }
    out
}

/// Prose outside a code span: each bare `--flag` (with an attached `=value`)
/// in code, the rest escaped.
fn prose_cell(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(pos) = rest.find("--") {
        let (before, from) = rest.split_at(pos);
        out.push_str(&escape_prose(before));
        let at_word_start = before.chars().last().is_none_or(|c| !c.is_alphanumeric());
        let name_len = flag_token_len(&from[2..]);
        if at_word_start && name_len > 0 {
            let mut end = 2 + name_len;
            if from[end..].starts_with('=') {
                end += 1 + flag_token_len(&from[end + 1..]);
            }
            out.push_str(&code(&from[..end]));
            rest = &from[end..];
        } else {
            out.push_str("--");
            rest = &from[2..];
        }
    }
    out.push_str(&escape_prose(rest));
    out
}

/// The length in bytes of the flag-shaped word that opens `text`: ASCII
/// letters, digits and underscores, with a hyphen counted only between two of
/// them, so `--core-dumps` is one word and a trailing dash is not.
fn flag_token_len(text: &str) -> usize {
    let bytes = text.as_bytes();
    let is_word_byte = |byte: u8| byte.is_ascii_alphanumeric() || byte == b'_';
    let mut len = 0;
    while len < bytes.len() {
        let byte = bytes[len];
        let inner_hyphen =
            byte == b'-' && len > 0 && bytes.get(len + 1).is_some_and(|next| is_word_byte(*next));
        if !(is_word_byte(byte) || inner_hyphen) {
            break;
        }
        len += 1;
    }
    len
}

fn escape_prose(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        if matches!(c, '|' | '<' | '>' | '_') {
            out.push('\\');
        }
        out.push(c);
    }
    out
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
        } else if id == TypeId::of::<f64>() {
            "f64"
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

/// `false` for a switch, `none` for no default, a number as written, an empty
/// default as `""`, any other default in code.
fn flag_default(arg: &Arg) -> String {
    if matches!(arg.get_action(), ArgAction::SetTrue) {
        return "false".to_string();
    }
    let defaults: Vec<String> = arg
        .get_default_values()
        .iter()
        .map(|v| {
            let value = v.to_string_lossy();
            if value.is_empty() {
                code("\"\"")
            } else if value.parse::<f64>().is_ok() {
                value.into_owned()
            } else {
                code(&value)
            }
        })
        .collect();
    if defaults.is_empty() {
        "none".to_string()
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

/// One `--help` heading and the table row of each flag under it.
type FlagSection = (String, Vec<Vec<String>>);

/// The flag rows of every `--help` heading, in the command's order.
fn flag_sections(command: &Command) -> Vec<FlagSection> {
    let mut sections: Vec<FlagSection> = Vec::new();
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
    sections
}

/// The flag table of one `--help` heading. A heading `GROUPS` names that no
/// flag carries is a stale `GROUPS`.
fn flag_table(out: &mut String, heading: &str, sections: &[FlagSection]) {
    let rows = sections
        .iter()
        .find(|(h, _)| h == heading)
        .map(|(_, rows)| rows)
        .unwrap_or_else(|| {
            panic!(
                "the --help heading {heading:?} named in GROUPS has no flags: a heading renamed \
                 in args.rs is renamed in GROUPS too"
            )
        });
    table(out, &["Flag", "Type", "Default", "Description"], rows);
}

/// A group page: its title, its introduction and one flag table per `--help`
/// heading it carries.
fn render_group(group: &Group, sections: &[FlagSection]) -> String {
    let mut out = String::new();
    writeln!(out, "# {}\n", group.title).unwrap();
    out.push_str(group.intro);
    for (heading, section) in group.sections {
        writeln!(out, "## {section}\n").unwrap();
        flag_table(&mut out, heading, sections);
    }
    out
}

/// The landing page: what the command line is, the group pages, how a flag is
/// read, and Ensembl VEP's option rules applied as data.
fn render_landing() -> String {
    let mut out = String::new();
    out.push_str("# CLI reference\n\n");
    out.push_str(
        "The `vep` binary accepts the Ensembl VEP command line. Every flag is Ensembl VEP's and \
         does what Ensembl VEP's does, or is refused with an error naming what it would do; the \
         few flags vep-rs adds are marked as its own. A flag the error table governs (refused, or \
         accepted without effect) carries that table's text as its help.\n\n",
    );
    writeln!(
        out,
        "These pages are generated from the command-line definition and the option tables by \
         the `cli_reference` test of `crates/vep-cli`. Edit the sources, then run `{REGENERATE}` \
         to rewrite them.\n"
    )
    .unwrap();
    out.push_str("```bash\n# Show all flags\nvep --help\n\n# Short form\nvep -i input.vcf -o output.txt --offline --json_cache /path/to/json_cache\n```\n\n");
    out.push_str(
        "Worked examples are on [Examples](examples.md); the output formats and their parity \
         with Ensembl VEP on [Output formats](output-formats.md); the divergences from Ensembl \
         VEP that a flag's description names on [Intended divergences](intended-divergences.md).\n\n",
    );

    out.push_str("## Option groups\n\n");
    out.push_str(
        "Every flag has a row on one of these pages, grouped as `--help` groups them: its names, \
         its value type, its default and its description. A default of `false` is a switch that \
         is off; `none` marks a flag with no declared default, and where a value applies in its \
         absence the description names it.\n\n",
    );
    let rows: Vec<Vec<String>> = GROUPS
        .iter()
        .map(|g| {
            vec![
                format!("[{}]({})", g.title, site_link(g.file)),
                cell(g.covers),
            ]
        })
        .collect();
    table(&mut out, &["Group", "Covers"], &rows);
    writeln!(
        out,
        "A repeated single-value flag keeps its last value, as Getopt::Long does. `-o` \
         accumulates and the last wins with a warning, as the `vep` script does (vep:32-33). A \
         long option is also read in any case, as a unique prefix and after one dash, and from \
         the `--config` file, `vep.ini` and the `VEP_*` variables ([configuration \
         sources]({}#configuration-sources)).\n",
        site_link(COMPATIBILITY)
    )
    .unwrap();
    out.push_str(
        "`--help` prints the flag list and exits 0, as does a command line with no argument at \
         all (vep:29); `--version` prints `vep <version>` and exits 0.\n\n",
    );

    out.push_str("## Option rules\n\n");
    out.push_str(
        "Ensembl VEP's option rules are applied as data, from `Config.pm:327-645` at \
         release/116.2, in this order after parsing. A flag the error table refuses is refused \
         first. Then the option sets switch on the flags they imply, in the order below, a set's \
         value overwriting an explicit one (`apply_option_sets`, `Config.pm:816-884`). Then the \
         checks of `check_config` (`Config.pm:944-985`) run: valid values, incompatible pairs, \
         requirements, deprecations, each refused with Ensembl VEP's own message. A rule naming \
         a flag vep-rs does not declare cannot fire and is listed so the table reads as Ensembl \
         VEP's.\n\n",
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
         `--regulatory` is switched on and annotates nothing: its columns appear in the tab and \
         VCF headers and stay empty, no `Extra` key is written for them, and no RegulatoryFeature \
         or MotifFeature row is written.\n\n",
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
    if !INCOMPATIBLE_LIFTED.is_empty() {
        out.push_str("Pairs Ensembl VEP refuses that are accepted here:\n\n");
        let rows: Vec<Vec<String>> = INCOMPATIBLE_LIFTED
            .iter()
            .map(|pair| {
                vec![
                    flag_in_code(pair.flag),
                    flag_in_code(pair.other),
                    cell(pair.reason),
                ]
            })
            .collect();
        table(
            &mut out,
            &["Flag", "Beside", "Why the pair is accepted"],
            &rows,
        );
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
    out.push_str(
        "Checked after Ensembl VEP's rules on the same state, for behaviour Ensembl VEP has no \
         flag for; each is refused with the line below.\n\n",
    );
    let rows: Vec<Vec<String>> = OWN_RULES
        .iter()
        .map(|rule| {
            vec![
                code(rule.component),
                code(&format!("ERROR: [{}] {}", rule.component, rule.message)),
            ]
        })
        .collect();
    table(&mut out, &["Component", "Error"], &rows);
    out
}

/// The compatibility page: the database flags, the error table, where
/// `Config.pm` is not followed, the configuration sources, the option
/// spellings and the exit status.
fn render_compatibility(group: &Group, command: &Command, sections: &[FlagSection]) -> String {
    let mut out = String::new();
    writeln!(out, "# {}\n", group.title).unwrap();
    out.push_str(group.intro);
    let [(heading, section)] = group.sections else {
        panic!("the compatibility page carries the database flags alone");
    };
    writeln!(out, "## {section}\n").unwrap();
    out.push_str(
        "vep-rs reads a local cache and never contacts a database. Each database flag is \
         declared so that the command line names it, and refused with the error its description \
         shows.\n\n",
    );
    flag_table(&mut out, heading, sections);

    out.push_str("## Consequence terms never assigned\n\n");
    out.push_str(
        "Six regulatory-feature terms need a regulatory-feature cache that is not loaded: \
         `TFBS_ablation`, `TFBS_amplification`, `TF_binding_site_variant`, \
         `regulatory_region_ablation`, `regulatory_region_amplification` and \
         `regulatory_region_variant`. `--regulatory`, alone or through `--everything`, is \
         accepted and warns that regulatory-feature annotation is not implemented.\n\n",
    );
    out.push_str(
        "The root term `sequence_variant` has no predicate, and Ensembl VEP's consequence engine \
         never assigns it either. The other 34 of Ensembl VEP's 41 Sequence Ontology terms are \
         assigned.\n\n",
    );

    out.push_str("## Input formats\n\n");
    out.push_str(
        "The input formats `--format` reads, how `guess` detects one, and the formats refused \
         are stated in its row under [Input](input-output.md#input).\n\n",
    );

    out.push_str("## Refused flags\n\n");
    out.push_str(
        "Every flag vep-rs refuses, with the line it prints and what Ensembl VEP does with the \
         flag (cited to its source at release/116.2). A refused flag exits 1 with the message \
         below under the `cli` component, and the help text of a declared one is this table's \
         message.\n\n",
    );
    let rows: Vec<Vec<String>> = ERROR_TABLE
        .iter()
        .filter_map(|entry| match entry.treatment {
            Treatment::Refused { component, message } => Some(vec![
                error_table_flag(entry.flag, entry.declared),
                code(&format!("ERROR: [{component}] {message}")),
                cell(entry.ensembl_vep),
            ]),
            Treatment::AcceptedWithoutEffect { .. } => None,
        })
        .collect();
    table(&mut out, &["Flag", "Error", "Ensembl VEP"], &rows);

    out.push_str("## Accepted without effect\n\n");
    out.push_str(
        "Every flag below parses and changes nothing. The first table is the error table's: \
         each flag's help is its reason, and the Ensembl VEP column cites what Ensembl VEP does \
         with it at release/116.2.\n\n",
    );
    let rows: Vec<Vec<String>> = ERROR_TABLE
        .iter()
        .filter_map(|entry| match entry.treatment {
            Treatment::AcceptedWithoutEffect { reason } => Some(vec![
                error_table_flag(entry.flag, entry.declared),
                cell(reason),
                cell(entry.ensembl_vep),
            ]),
            Treatment::Refused { .. } => None,
        })
        .collect();
    table(&mut out, &["Flag", "Reason", "Ensembl VEP"], &rows);
    out.push_str(
        "The flags below are declared with the description their group page shows, which says \
         what, if anything, each changes; the second column is what Ensembl VEP does with the \
         flag.\n\n",
    );
    let rows: Vec<Vec<String>> = compatibility_flags(command)
        .iter()
        .map(|(flag, ensembl_vep)| vec![flag_in_code(flag), cell(ensembl_vep)])
        .collect();
    table(&mut out, &["Flag", "Ensembl VEP"], &rows);

    out.push_str("## Where `Config.pm` is not followed\n\n");
    out.push_str(
        "- An output file matching `stdout` in any case does not force `--quiet` or drop \
         `--verbose` (`Config.pm:933-937`). The status lines go to standard error, so the \
         annotation on standard output is clean either way, and `--quiet --verbose` stays \
         refused with `-o STDOUT` as without.\n\
         - A run naming none of `--database`, `--cache`, `--offline` and `--custom` is not \
         refused (`Config.pm:993-1017`). vep-rs reads a local cache in every run, so those flags \
         select nothing and the run proceeds on the cache it is given.\n\
         - A `%name%` in an option-set value is replaced even when the flag's value is `0` or \
         empty (`Config.pm:867` leaves the reference in place then). No declared flag carries \
         such a value.\n\
         - A VCF output selected with `--output_format vcf` fills `SYMBOL`, `BIOTYPE`, `EXON` \
         and `INTRON` as `--vcf` does. Ensembl VEP's `vcf` option set applies to the `--vcf` \
         flag alone and leaves the CSQ fields empty under `--output_format vcf`.\n\
         - A value from the `--config` file, `vep.ini` or the environment never overwrites one \
         a higher source gives, even a `0`, an empty value or a bare `--cache`. Ensembl VEP's \
         `||=` (`Config.pm:1101, 1158`) lets a lower source overwrite a value Perl reads as \
         false, so `--shift_hgvs 0` on its command line yields to a file's `shift_hgvs 1`.\n\
         - A value from a file or the environment that the parser refuses (`fork many`, an \
         empty `format`) is an argument syntax error (exit 2). Ensembl VEP's `is_valid_param` \
         drops the line (`Config.pm:1046-1047`).\n\
         - A `--config` path that is a directory, or cannot be read for any reason, is the \
         `config` error; Ensembl VEP's `open` succeeds on a directory and reads nothing from \
         it. A `vep.ini` that is a directory is read as nothing here too. A byte outside UTF-8 \
         in a value becomes the replacement character, where Ensembl VEP passes the byte.\n\
         - A flag read from a file or the environment acts exactly as if typed, refusals \
         included. The status lines Ensembl VEP prints under `--verbose` while reading the \
         sources (`Read configuration from <file>`, `Set '<flag>=<value>' from environment \
         variable <VAR>`, `Config.pm:1107, 1145-1164`) are not printed.\n\
         - A leading `~` in `--config`, `--dir` and a file's `dir` is expanded to the home \
         directory, which Ensembl VEP leaves to the shell.\n\
         - `-V` is the version, where Getopt::Long's `ignore_case` reads it as `-v`, and an \
         upper-case short names no flag here. A short with its value attached (`-ofile.txt`) \
         and a cluster of switches (`-qv`) are clap's readings, where Getopt::Long refuses the \
         token. A `+` before an option name, which Getopt::Long accepts, is not an option \
         here.\n\n",
    );

    out.push_str("## Configuration sources\n\n");
    out.push_str(
        "A flag is read from four sources, the first that names it winning, as Ensembl VEP reads \
         them (`Config.pm:693-708`):\n\n",
    );
    writeln!(
        out,
        "1. the command line;\n\
         2. the `--config` file;\n\
         3. `<dir>/{INI_FILE}`, where `<dir>` is `--dir` (the command line's, else the \
         `--config` file's) and `{DEFAULT_DIR}` by default;\n\
         4. the environment variable `{ENV_PREFIX}<FLAG>` (`{ENV_PREFIX}DIR_PLUGINS` sets \
         `--dir_plugins`).\n"
    )
    .unwrap();
    writeln!(
        out,
        "A repeatable flag ({}) collects every source's values, the command line's first. No \
         source behind the command line can set {}.\n",
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
        "- One flag per line: the flag's long name without its dashes (one leading dash is \
         allowed), a blank, then its value. A line starting with `#`, or blank, is skipped.\n\
         - A switch is on with the value `1` and off with `0`. `cache` is on with `1`, a \
         directory or an empty quoted value and off with `0`. A key alone on its line is \
         dropped.\n\
         - A flag taking a value takes the first word after its name, or, for a repeatable \
         flag, the rest of the line. Blanks inside double quotes are part of a value and the \
         quotes are removed.\n\
         - A key is the long name exactly. A prefix, another case or a name that is no option \
         of Ensembl VEP's or of this version names nothing, and the line is dropped as Ensembl \
         VEP drops it. So is an alias, unless the alias is an option of Ensembl VEP's own \
         (`full_cache_dir` names `--json_cache`).\n\
         - A key of an Ensembl VEP option this version lacks (`stats_html 1`) is refused as the \
         flag on the command line is (exit 2). The exception is a switch of Ensembl VEP's set \
         to `0`: off there, it adds nothing here.\n\
         - Of several lines naming one flag, the first whose value is not `0` or empty counts.\n\
         - A `--config` file that cannot be opened is `ERROR: [config] Could not open config \
         file \"<file>\"`; `{INI_FILE}` is read only when it is a file.\n\n\
         ```\n# flags of a run\ncache 1\noffline 1\ndir_cache /data/vep\nspecies homo_sapiens\n\
         assembly GRCh38\nplugin LoFtool,/data/LoFtool_scores.txt\n```\n"
    )
    .unwrap();
    out.push_str("### Environment variables (`Config.pm:1127-1168`)\n\n");
    writeln!(
        out,
        "`{ENV_PREFIX}<FLAG>=<value>` sets `--<flag>`, the name lower-cased, when the value is \
         not empty. The name is read as a configuration file's key is: a name that is no option \
         of Ensembl VEP's or of this version is ignored, and one of an Ensembl VEP option this \
         version lacks is refused as the flag is (exit 2). The variables are read in name \
         order.\n"
    )
    .unwrap();
    out.push_str(
        "Environment variables vep-rs reads outside this mechanism, which Ensembl VEP has no \
         counterpart for:\n\n",
    );
    let rows: Vec<Vec<String>> = ENVIRONMENT
        .iter()
        .map(|(name, help)| vec![code(name), cell(help)])
        .collect();
    table(&mut out, &["Variable", "Meaning"], &rows);

    out.push_str("## Option spellings (`Config.pm:62`)\n\n");
    out.push_str(
        "- Getopt::Long's defaults apply to the command line: a long option is read in any case \
         (`--FORMAT vcf`), abbreviated to any prefix that names one option (`--form vcf`), and \
         after one dash as well as two (`-cache`, `-form vcf`).\n\
         - A value follows its option or is attached with `=`, and the token after an option \
         that requires a value is that value whatever it looks like (`--input_file -cache` \
         names the file `-cache`).\n\
         - A prefix is resolved against this version's flags first. Where a flag vep-rs adds \
         makes it ambiguous, it is resolved as Getopt::Long resolves it against Ensembl VEP's \
         `VEP_PARAMS` list (`Config.pm:73-276`): `--js` is `--json` though `--json_cache` also \
         starts with it, `--pa` `--password`, `--ve` `--verbose`.\n\
         - An Ensembl VEP option this version lacks is refused by its full name in every \
         spelling (`-stats_html` and `--stats_h` are refused as `--stats_html`, exit 2), never \
         read as a cluster of shorts.\n\
         - A prefix that names one flag of this version resolves to it even where Ensembl VEP's \
         larger set leaves it ambiguous (`--sy` is `--symbol`). A prefix of several flags that \
         Getopt::Long resolves no better (`--fo`), or a name no option has, is an argument \
         syntax error (exit 2).\n\
         - A single-letter short is clap's and is matched by case: `-v` is `--verbose`, `-V` \
         the version, `-qv` the cluster `-q -v`, and `-ofile.txt` the value `file.txt` of `-o`, \
         a spelling Getopt::Long refuses.\n\n",
    );
    out.push_str(
        "Of the spellings Getopt::Long resolves to an option this version declares, these \
         resolve differently here:\n\n",
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
         its own `ERROR: ` prefix. The components are listed below.\n\
         - `2`: an argument syntax error, in clap's words with its `error:` prefix: an unknown \
         flag, a prefix of several flags, a missing value, a value outside an enumeration.\n\n",
    );
    let rows: Vec<Vec<String>> = COMPONENTS
        .iter()
        .map(|(component, meaning)| vec![code(component), cell(meaning)])
        .collect();
    table(&mut out, &["Component", "Meaning"], &rows);
    out.push_str(
        "`--help` and `--version` win over a flag the error table refuses before parsing, \
         wherever they stand on the command line. Over a flag the parser does not know they win \
         only when they come first, since clap reads left to right. Ensembl VEP's Getopt::Long \
         runs with `pass_through` (Config.pm:62), so there `--help` wins over an unknown flag \
         wherever it stands (vep:29) and an unknown flag alone is refused (vep:30).\n\n",
    );
    out.push_str(
        "A panic aborts the process with the panic's own status. Ensembl VEP's own exit codes \
         (255 on `die`, 2 on a failed open) are not copied.\n",
    );
    out
}

fn error_table_flag(flag: &str, declared: bool) -> String {
    if declared {
        flag_in_code(flag)
    } else {
        format!("`--{flag}` (not declared; caught before parsing)")
    }
}

/// `COMPATIBILITY_FLAGS` checked against the help texts: every declared flag
/// whose help reads `Accepted for VEP compatibility` is listed, and every
/// listed flag still reads so.
fn compatibility_flags(command: &Command) -> &'static [(&'static str, &'static str)] {
    let from_help: Vec<&str> = command
        .get_arguments()
        .filter(|a| {
            flag_help(a)
                .to_ascii_lowercase()
                .contains(COMPATIBILITY_HELP)
        })
        .filter_map(|a| a.get_long())
        .collect();
    for flag in &from_help {
        assert!(
            COMPATIBILITY_FLAGS.iter().any(|(f, _)| f == flag),
            "--{flag} is accepted for compatibility: add it to COMPATIBILITY_FLAGS with what \
             Ensembl VEP does with it"
        );
    }
    for (flag, _) in COMPATIBILITY_FLAGS {
        assert!(
            from_help.contains(flag),
            "--{flag} no longer reads as accepted for compatibility: remove it from \
             COMPATIBILITY_FLAGS"
        );
    }
    COMPATIBILITY_FLAGS
}

/// Every page as `(file under docs/, text)`: the landing page and the group
/// pages. A `--help` heading no group carries is a stale `GROUPS`.
fn render(command: &Command) -> Vec<(&'static str, String)> {
    let sections = flag_sections(command);
    let mut pages = vec![(LANDING, render_landing())];
    for group in GROUPS {
        let text = if group.file == COMPATIBILITY {
            render_compatibility(group, command, &sections)
        } else {
            render_group(group, &sections)
        };
        pages.push((group.file, text));
    }
    for (heading, _) in &sections {
        assert!(
            GROUPS
                .iter()
                .any(|g| g.sections.iter().any(|(h, _)| h == heading)),
            "the --help heading {heading:?} belongs to no group page: add it to GROUPS"
        );
    }
    pages
}

#[test]
fn cli_reference_is_current() {
    let mut command = Args::command_with_error_table();
    command.build();
    let update = std::env::var_os("VEP_UPDATE_CLI_REFERENCE").is_some();
    for (file, rendered) in render(&command) {
        let path = repo_path(file);
        if update {
            fs::write(&path, &rendered).unwrap();
            continue;
        }
        let committed = fs::read_to_string(&path).unwrap_or_default();
        if committed != rendered {
            let first = committed
                .lines()
                .zip(rendered.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(committed.lines().count().min(rendered.lines().count()));
            panic!(
                "{file} is not the rendering of the command line and the tables: first \
                 difference at line {}\n  committed: {:?}\n  rendered:  {:?}\nrun {REGENERATE}",
                first + 1,
                committed.lines().nth(first).unwrap_or(""),
                rendered.lines().nth(first).unwrap_or("")
            );
        }
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
