// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP's option rules as data: the lists, option sets, valid values,
//! requirements, incompatibilities and deprecations of `Config.pm:327-645` at
//! release/116.2, and the flag state they act on.
//!
//! The tables name every flag Ensembl VEP declares, declared by vep-rs or not, so
//! they read as the Perl source does. [`FlagState`] is the Perl configuration
//! hash: one value per flag given on the command line, option sets applied in
//! order (`Config.pm:816-884`), then the checks of `check_config` (`:944-985`)
//! with Ensembl VEP's messages.

use std::collections::{BTreeMap, BTreeSet};

use clap::{ArgAction, ArgMatches, Command};
use vep_core::sv_routing::perl_split;

/// Flags whose single value is a comma-separated list (`Config.pm:327-336`).
pub const LIST_FLAGS: &[&str] = &[
    "individual",
    "individual_zyg",
    "cell_type",
    "pick_order",
    "fields",
    "chr",
    "distance",
    "dont_export",
];

/// Flags that may be given more than once and accumulate (`Config.pm:340-346`).
pub const ALLOW_MULTIPLE: &[&str] = &[
    "plugin",
    "custom",
    "phyloP",
    "phastCons",
    "transcript_filter",
];

/// One option set: when any flag in `flags` is active, every `(key, value)` of
/// `set` is applied, a `%name%` in a value standing for the current value of
/// `name` (`Config.pm:348-350`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OptionSet {
    pub flags: &'static [&'static str],
    pub set: &'static [(&'static str, &'static str)],
}

/// The option sets of `Config.pm:351-593`, in their order of application.
pub const OPTION_SETS: &[OptionSet] = &[
    OptionSet {
        flags: &["everything"],
        set: &[
            ("sift", "b"),
            ("polyphen", "b"),
            ("ccds", "1"),
            ("hgvs", "1"),
            ("symbol", "1"),
            ("numbers", "1"),
            ("domains", "1"),
            ("regulatory", "1"),
            ("canonical", "1"),
            ("protein", "1"),
            ("biotype", "1"),
            ("af", "1"),
            ("af_1kg", "1"),
            ("af_gnomade", "1"),
            ("af_gnomadg", "1"),
            ("max_af", "1"),
            ("pubmed", "1"),
            ("uniprot", "1"),
            ("mane", "1"),
            ("tsl", "1"),
            ("appris", "1"),
            ("variant_class", "1"),
            ("gene_phenotype", "1"),
            ("mirna", "1"),
        ],
    },
    OptionSet {
        flags: &["hgvs"],
        set: &[("hgvsc", "1"), ("hgvsp", "1")],
    },
    OptionSet {
        flags: &["genomes"],
        set: &[("host", "mysql-eg-publicsql.ebi.ac.uk"), ("port", "4157")],
    },
    OptionSet {
        flags: &["grch37"],
        set: &[("port", "3337"), ("assembly", "GRCh37")],
    },
    OptionSet {
        flags: &["refseq"],
        set: &[("core_type", "otherfeatures")],
    },
    OptionSet {
        flags: &["individual"],
        set: &[("allow_non_variant", "1")],
    },
    OptionSet {
        flags: &["individual_zyg"],
        set: &[("allow_non_variant", "1")],
    },
    OptionSet {
        flags: &["cell_type"],
        set: &[("regulatory", "1")],
    },
    OptionSet {
        flags: &["regulatory_gff"],
        set: &[("regulatory", "1")],
    },
    OptionSet {
        flags: &["filter_common"],
        set: &[("check_frequency", "1")],
    },
    OptionSet {
        flags: &[
            "check_frequency",
            "af",
            "af_1kg",
            "af_gnomad",
            "af_gnomade",
            "af_gnomadg",
            "max_af",
            "pubmed",
        ],
        set: &[("check_existing", "1")],
    },
    OptionSet {
        flags: &["show_cache_info"],
        set: &[("offline", "1"), ("format", "ensembl")],
    },
    OptionSet {
        flags: &["offline"],
        set: &[("cache", "1"), ("database", "0")],
    },
    OptionSet {
        flags: &["cache"],
        set: &[("no_slice_cache", "1")],
    },
    OptionSet {
        flags: &["convert"],
        set: &[("no_stats", "1")],
    },
    OptionSet {
        flags: &["json"],
        set: &[("rest", "1"), ("output_format", "json")],
    },
    OptionSet {
        flags: &["rest"],
        set: &[("no_escape", "1")],
    },
    OptionSet {
        flags: &["vcf"],
        set: &[
            ("output_format", "vcf"),
            ("symbol", "1"),
            ("biotype", "1"),
            ("numbers", "1"),
        ],
    },
    OptionSet {
        flags: &["tab"],
        set: &[("output_format", "tab")],
    },
    OptionSet {
        flags: &["humdiv"],
        set: &[("polyphen_analysis", "humdiv")],
    },
    OptionSet {
        flags: &["minimal"],
        set: &[("allele_number", "1")],
    },
    OptionSet {
        flags: &["gff"],
        set: &[("custom", "file=%gff%,format=gff")],
    },
    OptionSet {
        flags: &["gtf"],
        set: &[("custom", "file=%gtf%,format=gtf")],
    },
    OptionSet {
        flags: &["bigwig"],
        set: &[("custom", "%bigwig%,,bigwig,exact")],
    },
    OptionSet {
        flags: &["phyloP"],
        set: &[(
            "custom",
            "%ucsc_data_root%%ucsc_assembly%/phyloP%phyloP%way/%ucsc_assembly%.phyloP%phyloP%way.bw,phlyoP%phyloP%way,bigwig,exact",
        )],
    },
    OptionSet {
        flags: &["phastCons"],
        set: &[(
            "custom",
            "%ucsc_data_root%%ucsc_assembly%/phastCons%phastCons%way/%ucsc_assembly%.phastCons%phastCons%way.bw,phastCons%phastCons%way,bigwig,exact",
        )],
    },
    OptionSet {
        flags: &["exclude_predicted"],
        set: &[("transcript_filter", "not stable_id match ^X._")],
    },
    OptionSet {
        flags: &["bam"],
        set: &[("use_transcript_ref", "1"), ("bam_edited", "1")],
    },
    OptionSet {
        flags: &["mane"],
        set: &[("mane_select", "1"), ("mane_plus_clinical", "1")],
    },
    OptionSet {
        flags: &["use_given_ref"],
        set: &[("use_transcript_ref", "0")],
    },
];

/// The values a flag accepts (`Config.pm:596-604`); a declared flag carries the
/// same list as its clap value enumeration.
pub const VALID: &[(&str, &[&str])] = &[
    (
        "format",
        &["ensembl", "vcf", "hgvs", "id", "spdi", "region", "guess"],
    ),
    ("terms", &["SO", "display", "NCBI"]),
    ("sift", &["s", "p", "b"]),
    ("polyphen", &["s", "p", "b"]),
    (
        "pick_order",
        &[
            "mane_select",
            "mane_plus_clinical",
            "canonical",
            "appris",
            "tsl",
            "biotype",
            "ccds",
            "rank",
            "length",
            "ensembl",
            "refseq",
        ],
    ),
    ("nearest", &["transcript", "gene", "symbol"]),
    ("compress_output", &["gzip", "bgzip"]),
];

/// The values [`VALID`] lists for `flag`; empty for a flag outside the table.
pub fn valid_values(flag: &str) -> &'static [&'static str] {
    VALID
        .iter()
        .find(|(name, _)| *name == flag)
        .map_or(&[], |(_, values)| values)
}

/// Flags that need another flag set (`Config.pm:607-612`).
pub const REQUIRES: &[(&str, &[&str])] = &[
    ("original", &["filters"]),
    ("phyloP", &["ucsc_assembly"]),
    ("phastCons", &["ucsc_assembly"]),
    ("custom_multi_allelic", &["custom"]),
];

/// Flags that refuse each other (`Config.pm:615-632`). The table's order decides
/// the pair an error names when several apply; Ensembl VEP iterates a hash there,
/// so its choice varies from run to run.
pub const INCOMPATIBLE: &[(&str, &[&str])] = &[
    (
        "most_severe",
        &[
            "biotype",
            "no_intergenic",
            "protein",
            "symbol",
            "sift",
            "polyphen",
            "coding_only",
            "ccds",
            "mane",
            "mane_select",
            "canonical",
            "xref_refseq",
            "numbers",
            "domains",
            "tsl",
            "appris",
            "uniprot",
            "pick",
            "flag_pick",
            "pick_allele",
            "flag_pick_allele",
            "vcf",
            "summary",
        ],
    ),
    (
        "summary",
        &[
            "biotype",
            "no_intergenic",
            "protein",
            "symbol",
            "sift",
            "polyphen",
            "coding_only",
            "ccds",
            "mane",
            "mane_select",
            "canonical",
            "xref_refseq",
            "numbers",
            "domains",
            "tsl",
            "appris",
            "uniprot",
            "pick",
            "flag_pick",
            "pick_allele",
            "flag_pick_allele",
            "vcf",
            "most_severe",
        ],
    ),
    (
        "database",
        &[
            "af_1kg",
            "af_gnomad",
            "af_gnomade",
            "af_gnomadg",
            "af_exac",
            "max_af",
            "pubmed",
            "var_synonyms",
            "offline",
            "cache",
        ],
    ),
    ("af_gnomade", &["af_gnomad"]),
    ("quiet", &["verbose"]),
    ("refseq", &["gencode_basic", "merged"]),
    ("gencode_primary", &["gencode_basic", "refseq"]),
    ("json", &["vcf", "tab"]),
    ("vcf", &["json", "tab"]),
    ("tab", &["vcf", "json"]),
    ("individual", &["minimal", "individual_zyg"]),
    ("check_ref", &["lookup_ref"]),
    ("check_svs", &["offline"]),
    ("ga4gh_vrs", &["vcf"]),
    ("shift_hgvs", &["shift_3prime", "shift_genomic"]),
    ("lrg", &["offline"]),
];

/// An [`INCOMPATIBLE`] pair that is not refused here, with the reason: the
/// pair guards a resource vep-rs does without.
pub struct LiftedPair {
    pub flag: &'static str,
    pub other: &'static str,
    pub reason: &'static str,
}

/// The pairs of [`INCOMPATIBLE`] that hold in Ensembl VEP and not here.
pub const INCOMPATIBLE_LIFTED: &[LiftedPair] = &[LiftedPair {
    flag: "check_svs",
    other: "offline",
    reason: "Ensembl VEP reads the overlapping structural variants from its database, which \
             `--offline` forgoes (`AnnotationSourceAdaptor.pm:203-206`); `--check_svs <file>` \
             reads them from Ensembl's structural-variant VCF and needs no database, so the \
             pair is accepted and the file is read offline",
}];

/// Flags Ensembl VEP refuses as deprecated, with their replacement when one
/// exists (`Config.pm:636-638`).
pub const DEPRECATED: &[(&str, Option<&str>)] = &[("convert", None)];

/// Input formats exempt from the sorted-input check (`Config.pm:640-645`); the
/// check itself is not performed, so every format is annotated in input order.
pub const UNSORTABLE: &[&str] = &["id", "hgvs", "spdi", "region"];

/// A rule of vep-rs's own, checked after Ensembl VEP's on the same state.
pub struct OwnRule {
    pub component: &'static str,
    pub message: &'static str,
    pub applies: fn(&FlagState) -> bool,
}

/// Rules vep-rs adds for behaviour Ensembl VEP has no flag for.
pub const OWN_RULES: &[OwnRule] = &[
    OwnRule {
        component: "output",
        message: "--output_format parquet needs a file path: Parquet cannot be written to STDOUT",
        applies: |state| {
            state.scalar("output_format") == Some("parquet")
                && state
                    .last("output_file")
                    .is_some_and(|f| f.eq_ignore_ascii_case("STDOUT"))
        },
    },
    OwnRule {
        component: "output",
        message: "--fields applies to --tab and --vcf output; the Parquet column set is fixed",
        applies: |state| {
            state.scalar("output_format") == Some("parquet") && state.is_active("fields")
        },
    },
    OwnRule {
        component: "output",
        message: "--compress_output applies to the default, tab, VCF and JSON outputs; Parquet output is compressed by its own codec",
        applies: |state| {
            state.scalar("output_format") == Some("parquet") && state.is_active("compress_output")
        },
    },
];

/// A flag's value as the Perl configuration hash holds it: a scalar, or the
/// list a repeatable or comma-separated flag becomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FlagValue {
    Scalar(String),
    List(Vec<String>),
}

/// The configuration hash: every flag given, by its long name, with Ensembl
/// VEP's list conventions applied (`Config.pm:757-762`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FlagState {
    values: BTreeMap<String, FlagValue>,
}

/// Perl's truth of a scalar: everything but the empty string and `0`.
fn truthy(value: &str) -> bool {
    !(value.is_empty() || value == "0")
}

/// A comma-separated flag value as Perl's `split(',', ...)` reads it, the reading
/// of every `LIST_FLAGS` value (`Config.pm:760`) and of a scalar the value check
/// of `check_config` meets (`:946`): a trailing empty field is no field and a
/// middle one is kept, so `21,` is one name and `21,,22` three.
pub fn split_list_flag(value: &str) -> Vec<String> {
    perl_split(value, &[',']).map(str::to_string).collect()
}

impl FlagState {
    /// The state of a parsed command line: every argument with a long name that
    /// was given or defaulted, keyed by that name; a repeatable flag is a list,
    /// as is each `LIST_FLAGS` value split on commas; the `ALLOW_MULTIPLE`
    /// flags are present as lists even when not given.
    pub fn from_matches(command: &Command, matches: &ArgMatches) -> Self {
        let mut state = FlagState::default();
        for arg in command.get_arguments() {
            let Some(long) = arg.get_long() else {
                continue;
            };
            let id = arg.get_id().as_str();
            match arg.get_action() {
                ArgAction::SetTrue => {
                    if matches.get_flag(id) {
                        state.set(long, "1");
                    }
                }
                ArgAction::SetFalse
                | ArgAction::Count
                | ArgAction::Help
                | ArgAction::HelpShort
                | ArgAction::HelpLong
                | ArgAction::Version => {}
                action => {
                    let Ok(Some(raw)) = matches.try_get_raw(id) else {
                        continue;
                    };
                    let values: Vec<String> =
                        raw.map(|v| v.to_string_lossy().into_owned()).collect();
                    if matches!(action, ArgAction::Append) {
                        state
                            .values
                            .insert(long.to_string(), FlagValue::List(values));
                    } else if let Some(last) = values.last() {
                        state.set(long, last);
                    }
                }
            }
        }
        for flag in ALLOW_MULTIPLE {
            state
                .values
                .entry(flag.to_string())
                .or_insert_with(|| FlagValue::List(Vec::new()));
        }
        for flag in LIST_FLAGS {
            if let Some(FlagValue::Scalar(s)) = state.values.get(*flag) {
                let list = split_list_flag(s);
                state.values.insert(flag.to_string(), FlagValue::List(list));
            }
        }
        state
    }

    /// Sets a scalar value.
    pub fn set(&mut self, flag: &str, value: &str) {
        self.values
            .insert(flag.to_string(), FlagValue::Scalar(value.to_string()));
    }

    /// Whether the state holds a value for the flag, given or defaulted,
    /// whatever that value is.
    pub fn is_given(&self, flag: &str) -> bool {
        self.values.contains_key(flag)
    }

    /// Whether the flag is set to something true: a defined scalar other than
    /// the empty string or `0`, or a non-empty list (`Config.pm:902-906`).
    pub fn is_active(&self, flag: &str) -> bool {
        match self.values.get(flag) {
            None => false,
            Some(FlagValue::Scalar(s)) => truthy(s),
            Some(FlagValue::List(l)) => !l.is_empty(),
        }
    }

    /// The scalar value of a flag, when it has one.
    pub fn scalar(&self, flag: &str) -> Option<&str> {
        match self.values.get(flag)? {
            FlagValue::Scalar(s) => Some(s.as_str()),
            FlagValue::List(_) => None,
        }
    }

    /// The flag's value, or the last element of its list.
    pub fn last(&self, flag: &str) -> Option<&str> {
        match self.values.get(flag)? {
            FlagValue::Scalar(s) => Some(s.as_str()),
            FlagValue::List(l) => l.last().map(String::as_str),
        }
    }

    /// The list value of a flag, a scalar read as a one-element list.
    pub fn list(&self, flag: &str) -> Vec<String> {
        match self.values.get(flag) {
            None => Vec::new(),
            Some(FlagValue::Scalar(s)) => vec![s.clone()],
            Some(FlagValue::List(l)) => l.clone(),
        }
    }

    /// The value a `%name%` reference stands for: the `instance`th element of a
    /// list (the first when the list is shorter), or the scalar; `None` when the
    /// flag is unset (`Config.pm:846-864`).
    fn substitute(&self, name: &str, instance: usize) -> Option<String> {
        match self.values.get(name)? {
            FlagValue::Scalar(s) => Some(s.clone()),
            FlagValue::List(l) => l.get(instance).or_else(|| l.first()).cloned(),
        }
    }

    /// Replaces every `%name%` in a set value by the flag's value; an unset
    /// flag leaves its reference in place (`Config.pm:867`).
    fn expand(&self, value: &str, instance: usize) -> String {
        let mut out = String::with_capacity(value.len());
        let mut rest = value;
        while let Some(open) = rest.find('%') {
            out.push_str(&rest[..open]);
            let after = &rest[open + 1..];
            match after.find(|c: char| c == '%' || c.is_whitespace()) {
                Some(close) if after.as_bytes()[close] == b'%' && close > 0 => {
                    let name = &after[..close];
                    match self.substitute(name, instance) {
                        Some(v) => out.push_str(&v),
                        None => {
                            out.push('%');
                            out.push_str(name);
                            out.push('%');
                        }
                    }
                    rest = &after[close + 1..];
                }
                _ => {
                    out.push('%');
                    rest = after;
                }
            }
        }
        out.push_str(rest);
        out
    }

    /// Applies [`OPTION_SETS`] in order: for each active flag of a set, once per
    /// instance of a list-valued flag, every key of the set is assigned its
    /// expanded value, appended when the key already holds a list
    /// (`Config.pm:816-884`).
    pub fn apply_option_sets(&mut self) {
        for set in OPTION_SETS {
            for flag in set.flags {
                if !self.is_active(flag) {
                    continue;
                }
                let instances = self.list(flag).len().max(1);
                for instance in 0..instances {
                    for (key, value) in set.set {
                        let expanded = self.expand(value, instance);
                        match self.values.get_mut(*key) {
                            Some(FlagValue::List(list)) => list.push(expanded),
                            _ => self.set(key, &expanded),
                        }
                    }
                }
            }
        }
    }

    /// The checks of `check_config` (`Config.pm:944-985`) on the state after the
    /// option sets: a value outside its flag's [`VALID`] list, an
    /// [`INCOMPATIBLE`] pair (skipped under `--safe`, `:961`, and the
    /// [`INCOMPATIBLE_LIFTED`] pairs), a [`REQUIRES`] rule unmet, a
    /// [`DEPRECATED`] flag. The first failure is returned as Ensembl VEP's
    /// message without its `ERROR: ` prefix.
    pub fn check(&self, safe: bool) -> Result<(), String> {
        for (flag, allowed) in VALID {
            if !self.is_active(flag) {
                continue;
            }
            let values: Vec<String> = match self.values.get(*flag) {
                Some(FlagValue::List(l)) => l.clone(),
                Some(FlagValue::Scalar(s)) => split_list_flag(s),
                None => Vec::new(),
            };
            for value in values {
                if !allowed.contains(&value.as_str()) {
                    return Err(format!(
                        "\"{value}\" is not a valid value for --{flag}\nValid values: {}",
                        allowed.join(", ")
                    ));
                }
            }
        }
        if !safe {
            for (flag, invalid) in INCOMPATIBLE {
                if !self.is_active(flag) {
                    continue;
                }
                let lifted = |other: &str| {
                    INCOMPATIBLE_LIFTED
                        .iter()
                        .any(|pair| pair.flag == *flag && pair.other == other)
                };
                if let Some(other) = invalid.iter().find(|f| self.is_active(f) && !lifted(f)) {
                    return Err(format!("Can't use --{flag} and --{other} together"));
                }
            }
        }
        for (flag, required) in REQUIRES {
            if !self.is_active(flag) {
                continue;
            }
            if let Some(missing) = required.iter().find(|f| !self.is_active(f)) {
                return Err(format!("You must set --{missing} to use --{flag}"));
            }
        }
        for (flag, replacement) in DEPRECATED {
            if !self.is_active(flag) {
                continue;
            }
            return Err(match replacement {
                Some(new) => format!("--{flag} has been deprecated - please use --{new} instead"),
                None => format!("--{flag} has been deprecated"),
            });
        }
        Ok(())
    }

    /// The flags an option set turns on when `flag` alone is given, transitively
    /// and sorted by name: `everything` reaches `hgvsc` through `hgvs` and
    /// `mane_select` through `mane`.
    pub fn implied_by(flag: &str) -> Vec<&'static str> {
        let mut state = FlagState::default();
        state.set(flag, "1");
        state.apply_option_sets();
        let implied: BTreeSet<&'static str> = OPTION_SETS
            .iter()
            .flat_map(|set| set.set.iter())
            .map(|(key, _)| *key)
            .filter(|key| *key != flag && state.is_active(key))
            .collect();
        implied.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::{CommandFactory, Parser};

    /// `@LIST_FLAGS` and `@ALLOW_MULTIPLE`, Config.pm:327-346 at release/116.2.
    const CONFIG_PM_LISTS: &str = r#"our @LIST_FLAGS = qw(
  individual
  individual_zyg
  cell_type
  pick_order
  fields
  chr
  distance
  dont_export
);

# these flags can be specified more than once on the command line
# and so get turned into listrefs
our @ALLOW_MULTIPLE = qw(
  plugin
  custom
  phyloP
  phastCons
  transcript_filter
);"#;

    /// `%VALID`, `%REQUIRES`, `%INCOMPATIBLE` and `%DEPRECATED`, Config.pm:596-638.
    const CONFIG_PM_CHECKS: &str = r#"our %VALID = (
  format          => [qw(ensembl vcf hgvs id spdi region guess)],
  terms           => [qw(SO display NCBI)],
  sift            => [qw(s p b)],
  polyphen        => [qw(s p b)],
  pick_order      => [qw(mane_select mane_plus_clinical canonical appris tsl biotype ccds rank length ensembl refseq)],
  nearest         => [qw(transcript gene symbol)],
  compress_output => [qw(gzip bgzip)],
);

# these flags require others to be set, no sensible defaults can be set by us
our %REQUIRES = (
  original  => [qw(filters)],
  phyloP    => [qw(ucsc_assembly)],
  phastCons => [qw(ucsc_assembly)],
  custom_multi_allelic => [qw(custom)]
);

# incompatible options
our %INCOMPATIBLE = (
  most_severe => [qw(biotype no_intergenic protein symbol sift polyphen coding_only ccds mane mane_select canonical xref_refseq numbers domains tsl appris uniprot pick flag_pick pick_allele flag_pick_allele vcf summary)],
  summary     => [qw(biotype no_intergenic protein symbol sift polyphen coding_only ccds mane mane_select canonical xref_refseq numbers domains tsl appris uniprot pick flag_pick pick_allele flag_pick_allele vcf most_severe)],
  database    => [qw(af_1kg af_gnomad af_gnomade af_gnomadg af_exac max_af pubmed var_synonyms offline cache)],
  af_gnomade    => [qw(af_gnomad)],
  quiet       => [qw(verbose)],
  refseq      => [qw(gencode_basic merged)],
  gencode_primary => [qw(gencode_basic refseq)],
  json        => [qw(vcf tab)],
  vcf         => [qw(json tab)],
  tab         => [qw(vcf json)],
  individual  => [qw(minimal individual_zyg)],
  check_ref   => [qw(lookup_ref)],
  check_svs   => [qw(offline)],
  ga4gh_vrs   => [qw(vcf)],
  shift_hgvs  => [qw(shift_3prime shift_genomic)],
  lrg         => [qw(offline)],
);

# deprecated/replaced flags
# e.g. 'gmaf' => 'af'
our %DEPRECATED = (
  'convert' => undef,
);"#;

    /// `@OPTION_SETS`, Config.pm:351-593.
    const CONFIG_PM_OPTION_SETS: &str = r#"our @OPTION_SETS = (
  {
    flags => ['everything'],
    set   => {
      sift           => 'b',
      polyphen       => 'b',
      ccds           => 1,
      hgvs           => 1,
      symbol         => 1,
      numbers        => 1,
      domains        => 1,
      regulatory     => 1,
      canonical      => 1,
      protein        => 1,
      biotype        => 1,
      af             => 1,
      af_1kg         => 1,
      af_gnomade     => 1,
      af_gnomadg     => 1,
      max_af         => 1,
      pubmed         => 1,
      uniprot        => 1,
      mane           => 1,
      tsl            => 1,
      appris         => 1,
      variant_class  => 1,
      gene_phenotype => 1,
      mirna          => 1,
    },
  },

  {
    flags => ['hgvs'],
    set   => {
      hgvsc => 1,
      hgvsp => 1,
    },
  },

  {
    flags => ['genomes'],
    set   => {
      host => 'mysql-eg-publicsql.ebi.ac.uk',
      port => 4157,
    },
  },

  {
    flags => ['grch37'],
    set   => {
      port => 3337,
      assembly => 'GRCh37',
    },
  },

  {
    flags => ['refseq'],
    set   => {
      core_type => 'otherfeatures',
    },
  },

  {
    flags => ['individual'],
    set   => {
      allow_non_variant => 1,
    },
  },

  {
    flags => ['individual_zyg'],
    set   => {
      allow_non_variant => 1,
    },
  },

  {
    flags => ['cell_type'],
    set   => {
      regulatory => 1,
    },
  },

  {
    flags => ['regulatory_gff'],
    set   => {
      regulatory => 1,
    },
  },

  {
    flags => ['filter_common'],
    set   => {
      check_frequency => 1,
    },
  },

  {
    flags => [qw(check_frequency af af_1kg af_gnomad af_gnomade af_gnomadg max_af pubmed)],
    set   => {
      check_existing => 1,
    },
  },

  {
    flags => [qw(show_cache_info)],
    set   => {
      offline => 1,
      format => 'ensembl',
    }
  },

  {
    flags => [qw(offline)],
    set   => {
      cache => 1,
      database => 0,
    }
  },

  {
    flags => ['cache'],
    set   => {
      no_slice_cache => 1,
    }
  },

  {
    flags => ['convert'],
    set   => {
      no_stats => 1,
    }
  },

  {
    flags => ['json'],
    set   => {
      rest => 1,
      output_format => 'json',
    }
  },

  {
    flags => ['rest'],
    set   => {
      no_escape => 1,
    }
  },

  {
    flags => ['vcf'],
    set   => {
      output_format => 'vcf',
      symbol => 1,
      biotype => 1,
      numbers => 1
    }
  },

  {
    flags => ['tab'],
    set   => {
      output_format => 'tab',
    }
  },

  {
    flags => ['humdiv'],
    set   => {
      polyphen_analysis => 'humdiv',
    }
  },

  {
    flags => ['minimal'],
    set   => {
      allele_number => 1,
    }
  },

  {
    flags => ['gff'],
    set   => {
      custom => 'file=%gff%,format=gff'
    }
  },

  {
    flags => ['gtf'],
    set   => {
      custom => 'file=%gtf%,format=gtf'
    }
  },

  {
    flags => ['bigwig'],
    set   => {
      custom => '%bigwig%,,bigwig,exact'
    }
  },

  {
    flags => ['phyloP'],
    set   => {
      custom => '%ucsc_data_root%%ucsc_assembly%/phyloP%phyloP%way/%ucsc_assembly%.phyloP%phyloP%way.bw,phlyoP%phyloP%way,bigwig,exact'
    }
  },

  {
    flags => ['phastCons'],
    set   => {
      custom => '%ucsc_data_root%%ucsc_assembly%/phastCons%phastCons%way/%ucsc_assembly%.phastCons%phastCons%way.bw,phastCons%phastCons%way,bigwig,exact'
    }
  },

  {
    flags => ['exclude_predicted'],
    set   => {
      transcript_filter => 'not stable_id match ^X._'
    }
  },

  {
    flags => ['bam'],
    set   => {
      use_transcript_ref => 1,
      bam_edited => 1,
    }
  },
  {
    flags => ['mane'],
    set   => {
      mane_select => 1,
      mane_plus_clinical => 1,
    }
  },
  {
    flags => ['use_given_ref'],
    set   => {
      use_transcript_ref => 0,
    }
  },
);"#;

    /// The words of a Perl `qw(...)` list, or the one quoted word of `['x']`.
    fn words(list: &str) -> Vec<String> {
        let inner = list
            .trim()
            .trim_start_matches('[')
            .trim_end_matches(']')
            .trim();
        let inner = inner
            .strip_prefix("qw(")
            .and_then(|s| s.strip_suffix(')'))
            .unwrap_or(inner);
        inner
            .split(|c: char| c.is_whitespace() || c == ',')
            .map(|w| w.trim_matches('\''))
            .filter(|w| !w.is_empty())
            .map(str::to_string)
            .collect()
    }

    /// The `name => [list]` entries of a Perl hash declaration named `our %NAME`.
    fn hash_of_lists(source: &str, name: &str) -> Vec<(String, Vec<String>)> {
        let start = source
            .find(&format!("our %{name} = ("))
            .unwrap_or_else(|| panic!("no %{name}"));
        let body = &source[start..];
        let end = body.find("\n);").unwrap();
        body[..end]
            .lines()
            .skip(1)
            .filter(|l| l.contains("=>"))
            .map(|l| {
                let (k, v) = l.split_once("=>").unwrap();
                (
                    k.trim().trim_matches('\'').to_string(),
                    words(v.trim().trim_end_matches(',')),
                )
            })
            .collect()
    }

    /// The `our @NAME = qw(...)` words.
    fn array(source: &str, name: &str) -> Vec<String> {
        let start = source.find(&format!("our @{name} = qw(")).unwrap();
        let body = &source[start..];
        let end = body.find(");").unwrap();
        words(&body[body.find("qw(").unwrap()..=end])
    }

    /// One option set as the Perl list declares it: its flags and its (key, value) pairs.
    type PerlOptionSet = (Vec<String>, Vec<(String, String)>);

    /// The option sets parsed from the Perl list.
    fn perl_option_sets() -> Vec<PerlOptionSet> {
        let mut sets = Vec::new();
        for block in CONFIG_PM_OPTION_SETS.split("flags =>").skip(1) {
            let (flags, rest) = block.split_once("set   =>").unwrap();
            let flags = words(flags.trim().trim_end_matches(','));
            let body = rest.split_once('{').unwrap().1;
            let body = &body[..body.find('}').unwrap()];
            let set = body
                .lines()
                .filter(|l| l.contains("=>"))
                .map(|l| {
                    let (k, v) = l.split_once("=>").unwrap();
                    (
                        k.trim().to_string(),
                        v.trim()
                            .trim_end_matches(',')
                            .trim_matches('\'')
                            .to_string(),
                    )
                })
                .collect();
            sets.push((flags, set));
        }
        sets
    }

    fn owned(table: &[(&str, &[&str])]) -> Vec<(String, Vec<String>)> {
        table
            .iter()
            .map(|(k, v)| (k.to_string(), v.iter().map(|s| s.to_string()).collect()))
            .collect()
    }

    /// The flag state of a command line parsed by the `Args` command.
    fn parse(argv: &[&str]) -> FlagState {
        let mut command = Args::command();
        command.build();
        let matches = command.clone().try_get_matches_from(argv).unwrap();
        FlagState::from_matches(&command, &matches)
    }

    /// `LIST_FLAGS` and `ALLOW_MULTIPLE` are the Perl lists, word for word.
    #[test]
    fn list_flags_split_on_comma() {
        assert_eq!(array(CONFIG_PM_LISTS, "LIST_FLAGS"), LIST_FLAGS);
        assert_eq!(array(CONFIG_PM_LISTS, "ALLOW_MULTIPLE"), ALLOW_MULTIPLE);
        let state = parse(&[
            "vep",
            "--pick_order",
            "canonical,rank",
            "--distance",
            "1000,500",
        ]);
        assert_eq!(state.list("pick_order"), ["canonical", "rank"]);
        assert_eq!(state.list("distance"), ["1000", "500"]);
        assert!(state.is_active("pick_order"));
    }

    /// `Config.pm:760` turns a `@LIST_FLAGS` value into a list with Perl's
    /// `split`, which drops trailing empty fields and keeps a middle one: `--chr
    /// 21,` names one chromosome, `--chr 21,,22` three with the middle one empty,
    /// and `--individual S1,` or `--individual_zyg S1,` one sample. The value
    /// check of `check_config` splits a scalar the same way (`:946`), so a
    /// trailing comma passes it and an empty middle field is refused by name.
    #[test]
    fn list_flags_split_as_perl_does() {
        assert_eq!(parse(&["vep", "--chr", "21,"]).list("chr"), ["21"]);
        assert_eq!(
            parse(&["vep", "--chr", "21,,22"]).list("chr"),
            ["21", "", "22"]
        );
        assert!(parse(&["vep", "--chr", ","]).list("chr").is_empty());
        assert_eq!(
            parse(&["vep", "--individual", "S1,"]).list("individual"),
            ["S1"]
        );
        assert_eq!(
            parse(&["vep", "--individual_zyg", "S1,"]).list("individual_zyg"),
            ["S1"]
        );
        let mut state = FlagState::default();
        state.set("terms", "SO,");
        assert!(state.check(false).is_ok());
        state.set("terms", "SO,,NCBI");
        assert_eq!(
            state.check(false).unwrap_err(),
            "\"\" is not a valid value for --terms\nValid values: SO, display, NCBI"
        );
    }

    /// `%INCOMPATIBLE` is the Perl table, pair for pair and in its order.
    #[test]
    fn incompatible_pairs_match_config_pm() {
        assert_eq!(
            hash_of_lists(CONFIG_PM_CHECKS, "INCOMPATIBLE"),
            owned(INCOMPATIBLE)
        );
        let err = parse(&["vep", "--quiet", "--verbose"])
            .check(false)
            .unwrap_err();
        assert_eq!(err, "Can't use --quiet and --verbose together");
        assert_eq!(
            parse(&["vep", "--json", "--vcf"]).check(false).unwrap_err(),
            "Can't use --json and --vcf together"
        );
        // Under --everything the set turns on --biotype, the first entry of the
        // most_severe list, so Ensembl VEP names that pair.
        let mut state = parse(&["vep", "--most_severe", "--everything"]);
        state.apply_option_sets();
        assert_eq!(
            state.check(false).unwrap_err(),
            "Can't use --most_severe and --biotype together"
        );
        assert!(
            state.check(true).is_ok(),
            "--safe skips the incompatibility check"
        );
    }

    /// The lifted pairs are pairs of the Perl table, and the one lifted pair
    /// passes the check while its neighbour in the table, `lrg` beside
    /// `offline`, is still refused.
    #[test]
    fn lifted_pairs_are_perl_pairs_the_check_passes() {
        for pair in INCOMPATIBLE_LIFTED {
            assert!(
                INCOMPATIBLE
                    .iter()
                    .any(|(flag, others)| *flag == pair.flag && others.contains(&pair.other)),
                "{} beside {} is not a pair of %INCOMPATIBLE",
                pair.flag,
                pair.other
            );
        }
        let mut state = parse(&["vep", "--check_svs", "svs.vcf.gz", "--offline"]);
        state.apply_option_sets();
        assert!(state.check(false).is_ok());
        let mut state = parse(&["vep", "--lrg", "--offline"]);
        state.apply_option_sets();
        assert_eq!(
            state.check(false).unwrap_err(),
            "Can't use --lrg and --offline together"
        );
    }

    /// `%REQUIRES`, `%VALID` and `%DEPRECATED` are the Perl tables.
    #[test]
    fn requires_match_config_pm() {
        assert_eq!(hash_of_lists(CONFIG_PM_CHECKS, "REQUIRES"), owned(REQUIRES));
        assert_eq!(hash_of_lists(CONFIG_PM_CHECKS, "VALID"), owned(VALID));
        let deprecated = hash_of_lists(CONFIG_PM_CHECKS, "DEPRECATED");
        assert_eq!(deprecated.len(), DEPRECATED.len());
        for ((flag, replacement), (perl_flag, perl_value)) in DEPRECATED.iter().zip(&deprecated) {
            assert_eq!(flag, perl_flag);
            let perl_replacement = (perl_value != &["undef"]).then(|| perl_value[0].clone());
            assert_eq!(replacement.map(str::to_string), perl_replacement);
        }
        let mut state = FlagState::default();
        state.set("phyloP", "100");
        assert_eq!(
            state.check(false).unwrap_err(),
            "You must set --ucsc_assembly to use --phyloP"
        );
        let mut state = FlagState::default();
        state.set("convert", "vcf");
        assert_eq!(
            state.check(false).unwrap_err(),
            "--convert has been deprecated"
        );
        let mut state = FlagState::default();
        state.set("terms", "HGVS");
        assert_eq!(
            state.check(false).unwrap_err(),
            "\"HGVS\" is not a valid value for --terms\nValid values: SO, display, NCBI"
        );
    }

    /// `OPTION_SETS` is the Perl list, set for set and in its order, and applies
    /// as `apply_option_sets` does: `--grch37` sets the assembly over an explicit
    /// `--assembly`, `--hgvs` reaches `hgvsc`, a `%name%` takes the flag's value
    /// per instance and stays when the flag is unset.
    #[test]
    fn option_sets_match_config_pm() {
        let perl = perl_option_sets();
        assert_eq!(perl.len(), OPTION_SETS.len());
        for (set, (flags, pairs)) in OPTION_SETS.iter().zip(&perl) {
            assert_eq!(set.flags, flags.as_slice());
            let ours: Vec<(String, String)> = set
                .set
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect();
            assert_eq!(ours, *pairs, "set for {flags:?}");
        }
        let mut state = parse(&["vep", "--grch37", "--assembly", "GRCh38"]);
        state.apply_option_sets();
        assert_eq!(state.scalar("assembly"), Some("GRCh37"));
        assert_eq!(state.scalar("port"), Some("3337"));
        let mut state = parse(&["vep", "--everything"]);
        state.apply_option_sets();
        for flag in [
            "hgvsc",
            "hgvsp",
            "mane_select",
            "mane_plus_clinical",
            "check_existing",
        ] {
            assert!(state.is_active(flag), "{flag}");
        }
        assert_eq!(state.scalar("sift"), Some("b"));
        let mut state = FlagState::default();
        state.values.insert(
            "gff".to_string(),
            FlagValue::List(vec!["a.gff".to_string(), "b.gff".to_string()]),
        );
        state
            .values
            .insert("custom".to_string(), FlagValue::List(Vec::new()));
        state.apply_option_sets();
        assert_eq!(
            state.list("custom"),
            ["file=a.gff,format=gff", "file=b.gff,format=gff"]
        );
        let mut state = FlagState::default();
        state.set("phyloP", "100");
        state
            .values
            .insert("custom".to_string(), FlagValue::List(Vec::new()));
        state.apply_option_sets();
        assert_eq!(
            state.list("custom"),
            ["%ucsc_data_root%%ucsc_assembly%/phyloP100way/%ucsc_assembly%.phyloP100way.bw,phlyoP100way,bigwig,exact"]
        );
    }

    /// `--plugin` and `--custom` are `ALLOW_MULTIPLE`: each occurrence appends.
    #[test]
    fn repeated_plugin_and_custom_accumulate() {
        let state = parse(&[
            "vep",
            "--plugin",
            "CADD,a",
            "--plugin",
            "REVEL,b",
            "--custom",
            "x.bed,X,bed",
        ]);
        assert_eq!(state.list("plugin"), ["CADD,a", "REVEL,b"]);
        assert_eq!(state.list("custom"), ["x.bed,X,bed"]);
        assert!(parse(&["vep"]).values.get("plugin") == Some(&FlagValue::List(Vec::new())));
    }

    /// Every `%VALID` flag vep-rs declares carries the Perl list as its clap
    /// value enumeration, so a value outside it is refused at parse time, and
    /// `--output_format`, vep-rs's own, is an enumeration as well: Ensembl VEP's
    /// `%FORMAT_MAP` names (`ensembl` and `vep` the same format) plus `parquet`.
    #[test]
    fn output_format_validated() {
        let mut command = Args::command();
        command.build();
        let declared: Vec<&clap::Arg> = command.get_arguments().collect();
        let mut checked = 0;
        for (flag, allowed) in VALID {
            let Some(arg) = declared.iter().find(|a| a.get_long() == Some(flag)) else {
                continue;
            };
            let possible: Vec<String> = arg
                .get_possible_values()
                .iter()
                .map(|p| p.get_name().to_string())
                .collect();
            assert_eq!(possible, *allowed, "--{flag}");
            checked += 1;
        }
        assert!(checked >= 5, "{checked} VALID flags declared");
        let output_format = declared
            .iter()
            .find(|a| a.get_long() == Some("output_format"))
            .unwrap();
        let possible: Vec<String> = output_format
            .get_possible_values()
            .iter()
            .map(|p| p.get_name().to_string())
            .collect();
        assert_eq!(
            possible,
            ["vep", "ensembl", "vcf", "json", "tab", "parquet"]
        );
        let err = Args::try_parse_from(["vep", "--output_format", "xml"]).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::InvalidValue);
        let err = Args::try_parse_from(["vep", "--format", "caid"]).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::InvalidValue);
    }

    /// `--grch37` is the option set of Config.pm:398-404: the assembly becomes
    /// `GRCh37`; the port it also sets belongs to the database connection.
    #[test]
    fn grch37_sets_assembly() {
        let mut state = parse(&["vep", "--grch37"]);
        assert_eq!(state.scalar("assembly"), None);
        state.apply_option_sets();
        assert_eq!(state.scalar("assembly"), Some("GRCh37"));
        assert_eq!(FlagState::implied_by("grch37"), ["assembly", "port"]);
    }
}
