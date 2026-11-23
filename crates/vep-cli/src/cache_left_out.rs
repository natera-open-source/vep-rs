// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines a run prints for what its cache does not carry and its command
//! asks for. The cache's own record states each part: the names of
//! `left_out` in `info.json`, which the builder and the converter write
//! (`vep_core::cache_info::left_out`), and the absence of `variation_cols`.
//! A cache whose record holds no `left_out` key draws none of the lines read
//! from it.

use vep_core::cache_info::{left_out, CacheInfo};

use crate::config::Config;
use crate::variation_source::VariationSource;

/// The flags `--everything` sets among those this module names
/// (`option_rules::OPTION_SETS`).
const SET_BY_EVERYTHING: [&str; 14] = [
    "symbol",
    "domains",
    "uniprot",
    "mane",
    "appris",
    "gene_phenotype",
    "mirna",
    "af",
    "af_1kg",
    "af_gnomade",
    "af_gnomadg",
    "max_af",
    "pubmed",
    "regulatory",
];

/// The flags that ask for a column a cache built without the core table
/// dumps holds no value for, each with its columns of that kind
/// (`FLAG_FIELD_GROUPS` of `vep-io`).
fn core_table_columns(config: &Config) -> [(&'static str, bool, &'static [&'static str]); 9] {
    [
        ("symbol", config.symbol, &["SYMBOL_SOURCE", "HGNC_ID"]),
        ("mane_select", config.mane_select, &["MANE_SELECT"]),
        ("mane", config.mane, &["MANE_SELECT", "MANE_PLUS_CLINICAL"]),
        ("appris", config.appris, &["APPRIS"]),
        ("mirna", config.mirna, &["miRNA"]),
        ("xref_refseq", config.xref_refseq, &["RefSeq"]),
        (
            "uniprot",
            config.uniprot,
            &["SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM"],
        ),
        ("domains", config.domains, &["DOMAINS"]),
        ("gene_phenotype", config.gene_phenotype, &["GENE_PHENO"]),
    ]
}

/// Whether a token of the command line is the flag `name`, with one dash or
/// two and with or without a value after `=`. [`Config::command_line`] lists
/// the flags of the configuration file and of the environment beside the
/// command line's own.
fn named(config: &Config, name: &str) -> bool {
    config.command_line.iter().any(|token| {
        let flag = token.trim_start_matches('-');
        token.starts_with('-') && flag.split('=').next() == Some(name)
    })
}

/// The flag of the command that switches `flag` on, as a line names it: the
/// flag itself where the command names it, else the flag that sets it
/// (`--everything`; `--mane`, which sets `mane_select`; `--vcf`, whose layout
/// holds the symbol fields; `--filter_common`), else the flag itself.
fn asking(config: &Config, flag: &str) -> String {
    let set_by = if named(config, flag) {
        flag
    } else if flag == "mane_select" && config.mane {
        return asking(config, "mane");
    } else if config.everything && SET_BY_EVERYTHING.contains(&flag) {
        "everything"
    } else if flag == "symbol" && config.output_format == "vcf" {
        "vcf"
    } else if flag == "check_frequency" && named(config, "filter_common") {
        "filter_common"
    } else {
        flag
    };
    format!("--{set_by}")
}

/// The flags of `flags` that are on, each as [`asking`] names it, once.
fn asked_by(config: &Config, flags: &[(&str, bool)]) -> Vec<String> {
    let mut asked: Vec<String> = Vec::new();
    for (flag, on) in flags {
        let flag = asking(config, flag);
        if *on && !asked.contains(&flag) {
            asked.push(flag);
        }
    }
    asked
}

/// `items` as a line lists them: `a`, `a and b`, `a, b and c`.
fn listed<S: AsRef<str>>(items: &[S]) -> String {
    let items: Vec<&str> = items.iter().map(AsRef::as_ref).collect();
    match items.split_last() {
        None => String::new(),
        Some((only, [])) => only.to_string(),
        Some((last, rest)) => format!("{} and {last}", rest.join(", ")),
    }
}

/// Whether `--check_frequency` reads the cache's `AF` column: the population
/// `1KG_ALL`, the default and the one `--filter_common` reads
/// (`frequency_filter`).
fn filter_reads_af(config: &Config) -> bool {
    config.check_frequency && config.freq_pop.eq_ignore_ascii_case("1KG_ALL")
}

/// Whether the output writes the column `name`: a tab or a VCF output under
/// `--fields` writes the columns its list names, and every other output each
/// column its flags switch on (`Runner::run_stages`, `writes_field`). A line
/// names a column only where the output writes it.
fn written(config: &Config, name: &str) -> bool {
    !matches!(config.output_format.as_str(), "tab" | "vcf")
        || config
            .fields
            .as_ref()
            .is_none_or(|fields| fields.iter().any(|field| field == name))
}

/// Whether the output writes a column of the gnomAD set `prefix` names
/// (`gnomADe`, `gnomADg`): [`written`] for the set's columns, named
/// `<prefix>_AF` and `<prefix>_<population>_AF`.
fn writes_gnomad(config: &Config, prefix: &str) -> bool {
    !matches!(config.output_format.as_str(), "tab" | "vcf")
        || config.fields.as_ref().is_none_or(|fields| {
            fields.iter().any(|field| {
                field
                    .strip_prefix(prefix)
                    .is_some_and(|rest| rest.starts_with('_'))
            })
        })
}

/// The flags whose option set adds a `--custom` source
/// (`option_rules::OPTION_SETS`).
const ADD_A_CUSTOM_SOURCE: [&str; 5] = ["gff", "gtf", "bigwig", "phyloP", "phastCons"];

/// The route to a cell the builder left out, as its closing lines name it.
const CONVERSION_CARRIES: &str =
    "a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// The caches that hold variation data, as the lines of a cache without it
/// name them.
const VARIATION_DATA_HELD_BY: &str =
    "a cache from vep-cache-builder --auto or --include-variations, or a conversion of Ensembl \
     VEP's cache, holds it (docs/cache-variation-data.md)";

/// The 1000 Genomes columns of a cache's variation files, each beside the
/// output column that prints it.
const THOUSAND_GENOMES_COLUMNS: [(&str, &str); 6] = [
    ("AF", "AF"),
    ("AFR", "AFR_AF"),
    ("AMR", "AMR_AF"),
    ("EAS", "EAS_AF"),
    ("EUR", "EUR_AF"),
    ("SAS", "SAS_AF"),
];

/// The frequency columns a cache's record declares, by the three groups a
/// release-116 cache carries.
struct FrequencyColumns<'a> {
    columns: &'a [String],
}

impl FrequencyColumns<'_> {
    fn holds(&self, column: &str) -> bool {
        self.columns.iter().any(|declared| declared == column)
    }

    /// Whether a population of the 1000 Genomes continental columns is
    /// declared; `AF`, the panel's own, is no population `MAX_AF` reads.
    fn thousand_genomes_populations(&self) -> bool {
        THOUSAND_GENOMES_COLUMNS[1..]
            .iter()
            .any(|(column, _)| self.holds(column))
    }

    /// Whether a column of the gnomAD set `prefix` names (`gnomADe`,
    /// `gnomADg`) is declared, the combined column or a population's.
    fn gnomad(&self, prefix: &str) -> bool {
        self.columns.iter().any(|declared| {
            declared
                .strip_prefix(prefix)
                .is_some_and(|rest| rest.is_empty() || rest.starts_with('_'))
        })
    }
}

/// The lines of a run on a cache that holds variation data for the cells
/// its command asks for and the cache does not carry: one line naming every
/// column that stays empty, and one for a `MAX_AF` that ranges over fewer
/// populations than Ensembl VEP's own cache of the genome gives it. The
/// frequency columns are read from the record's `variation_cols`; the
/// PubMed identifiers, the variant synonyms and the motif cells from its
/// `left_out` names, which a build writes for a cell it holds no value for.
/// The route is named where `left_out` names a cell of the line, and a
/// column is named only where the output writes it ([`written`]). A run under
/// `--regulatory_gff` reads those files in place of the cache's regulatory
/// tree (`RegulatorySource::open_for`), so the motif cells draw no line there.
fn cell_notices(config: &Config, info: &CacheInfo) -> Vec<String> {
    let mut notices = Vec::new();
    let mut columns: Vec<&str> = Vec::new();
    let mut flags: Vec<(&str, bool)> = Vec::new();
    let mut routed = false;
    let mut filter = "";
    let declared = FrequencyColumns {
        columns: &info.variation_cols,
    };
    let variation_data = info.variation_columns().is_some() && config.check_existing;
    if variation_data {
        for (flag, on, of_flag) in [
            ("af", config.af, &THOUSAND_GENOMES_COLUMNS[..1]),
            ("af_1kg", config.af_1kg, &THOUSAND_GENOMES_COLUMNS[..]),
        ] {
            let mut missing = false;
            for (column, printed) in of_flag.iter().filter(|_| on) {
                if !declared.holds(column) && written(config, printed) {
                    missing = true;
                    if !columns.contains(printed) {
                        columns.push(printed);
                    }
                }
            }
            if missing {
                flags.push((flag, true));
                routed |= info.leaves_out(left_out::THOUSAND_GENOMES_FREQUENCIES);
            }
        }
        for (flag, on, prefix, printed) in [
            (
                "af_gnomade",
                config.af_gnomade || config.af_gnomad,
                "gnomADe",
                "the gnomAD exome columns",
            ),
            (
                "af_gnomadg",
                config.af_gnomadg,
                "gnomADg",
                "the gnomAD genome columns",
            ),
        ] {
            if on && !declared.gnomad(prefix) && writes_gnomad(config, prefix) {
                columns.push(printed);
                let exomes_by_the_older_flag = flag == "af_gnomade" && !config.af_gnomade;
                flags.push((
                    if exomes_by_the_older_flag {
                        "af_gnomad"
                    } else {
                        flag
                    },
                    true,
                ));
                routed |= info.leaves_out(left_out::GNOMAD_FREQUENCIES);
            }
        }
        let populations = [
            declared.thousand_genomes_populations(),
            declared.gnomad("gnomADe"),
            declared.gnomad("gnomADg"),
        ];
        let maximum: Vec<&str> = ["MAX_AF", "MAX_AF_POPS"]
            .into_iter()
            .filter(|column| config.max_af && written(config, column))
            .collect();
        if !maximum.is_empty() && !populations.contains(&true) {
            columns.extend(&maximum);
            flags.push(("max_af", true));
            routed |= info.leaves_out(left_out::THOUSAND_GENOMES_FREQUENCIES)
                || info.leaves_out(left_out::GNOMAD_FREQUENCIES);
        }
        if filter_reads_af(config) && !declared.holds("AF") {
            flags.push(("check_frequency", true));
            routed |= info.leaves_out(left_out::THOUSAND_GENOMES_FREQUENCIES);
            filter = if config.freq_filter == "include" {
                "the frequency filter finds no AF and drops every record"
            } else {
                "the frequency filter finds no AF and keeps every record"
            };
        }
        if !maximum.is_empty() && populations.contains(&true) && populations.contains(&false) {
            let names = ["1000 Genomes", "gnomAD exome", "gnomAD genome"];
            let of = |held: bool| -> Vec<&str> {
                names
                    .iter()
                    .zip(populations)
                    .filter(|(_, declared)| *declared == held)
                    .map(|(name, _)| *name)
                    .collect()
            };
            let route = if info.leaves_out(left_out::THOUSAND_GENOMES_FREQUENCIES)
                || info.leaves_out(left_out::GNOMAD_FREQUENCIES)
            {
                format!("; {CONVERSION_CARRIES}")
            } else {
                String::new()
            };
            notices.push(format!(
                "INFO: {} ({}) {} over the {} populations alone on this cache, whose variation \
                 data holds no {} frequency{route}",
                listed(&maximum),
                asking(config, "max_af"),
                if maximum.len() == 1 {
                    "ranges"
                } else {
                    "range"
                },
                listed(&of(true)),
                of(false).join(" or "),
            ));
        }
        for (name, flag, on, of_name) in [
            (left_out::PUBMED, "pubmed", config.pubmed, "PUBMED"),
            (
                left_out::VAR_SYNONYMS,
                "var_synonyms",
                config.var_synonyms,
                "VAR_SYNONYMS",
            ),
        ] {
            if on && info.leaves_out(name) && written(config, of_name) {
                columns.push(of_name);
                flags.push((flag, true));
                routed = true;
            }
        }
    }
    if info.regulatory
        && config.regulatory
        && config.regulatory_gff.is_none()
        && info.leaves_out(left_out::MOTIF_SCORES)
    {
        let motif_cells: Vec<&str> = ["HIGH_INF_POS", "MOTIF_SCORE_CHANGE"]
            .into_iter()
            .filter(|column| written(config, column))
            .collect();
        if !motif_cells.is_empty() {
            columns.extend(motif_cells);
            flags.push(("regulatory", true));
            routed = true;
        }
    }
    if columns.is_empty() && filter.is_empty() {
        return notices;
    }
    let asked = asked_by(config, &flags);
    let route = if routed {
        format!("; {CONVERSION_CARRIES}")
    } else {
        String::new()
    };
    let line = if columns.is_empty() {
        format!(
            "INFO: {} {} AF, which this cache does not carry: {filter}{route}",
            listed(&asked),
            if asked.len() == 1 { "reads" } else { "read" },
        )
    } else {
        let stay = if columns.len() == 1 && !columns[0].starts_with("the ") {
            "the column stays"
        } else {
            "the columns stay"
        };
        let filter = if filter.is_empty() {
            String::new()
        } else {
            format!(" and {filter}")
        };
        format!(
            "INFO: {} {} for {}, which this cache does not carry: {stay} empty{filter}{route}",
            listed(&asked),
            if asked.len() == 1 { "asks" } else { "ask" },
            listed(&columns),
        )
    };
    notices.insert(0, line);
    notices
}

/// The lines of a run for what the cache's `left_out` list names and the
/// command asks for, each printed once as the cache is opened. `custom` is
/// whether the run names a `--custom` source, which adds the `SOURCE`
/// column.
///
/// - `core_tables`: the columns the command's flags ask for among those a
///   build without the core table dumps holds no value for.
/// - `source`: the `SOURCE` column, where a `--custom` source, `--merged` or
///   a `--fields` list of a tab or VCF output writes it; a source is named by
///   the flag of the command that adds it (`ADD_A_CUSTOM_SOURCE`), and by
///   `--custom` where the command names none of them. A JSON output holds
///   no such column: an entry carries `source` only where its transcript
///   names one, so a JSON run, the annotation run of `vep recode` among
///   them, draws no line.
/// - `panel_frequency`: the `AF` column, where `--af`, `--af_1kg` or a
///   frequency filter on `1KG_ALL` reads it from the cache's variation files.
/// - the cells of a cache with variation data or a regulatory tree
///   (`cell_notices`): the frequency columns the record's `variation_cols`
///   do not declare, and the cells `left_out` names.
/// - `seq_edits`: every run, since the reference peptide of each coding row
///   is read without the edits.
///
/// A column is named only where the output writes it ([`written`]): a tab or
/// a VCF output under `--fields` draws a line for the columns its list names.
pub(crate) fn notices(config: &Config, info: &CacheInfo, custom: bool) -> Vec<String> {
    let mut notices = Vec::new();
    if info.leaves_out(left_out::CORE_TABLES) {
        let groups = core_table_columns(config);
        let mut columns: Vec<&str> = Vec::new();
        let mut flags: Vec<(&str, bool)> = Vec::new();
        for (flag, on, of_flag) in &groups {
            let mut asks = false;
            for column in of_flag
                .iter()
                .filter(|column| *on && written(config, column))
            {
                asks = true;
                if !columns.contains(column) {
                    columns.push(column);
                }
            }
            flags.push((flag, asks));
        }
        if !columns.is_empty() {
            let asked = asked_by(config, &flags);
            let stay = if columns.len() == 1 {
                "the column stays"
            } else {
                "the columns stay"
            };
            let arms = if columns.contains(&"miRNA") {
                " and mature_miRNA_variant is never assigned"
            } else {
                ""
            };
            notices.push(format!(
                "INFO: {} {} for {}, which this cache does not carry: it was built without the \
                 core table dumps (vep-cache-builder --core-dumps), so {stay} empty{arms}",
                listed(&asked),
                if asked.len() == 1 { "asks" } else { "ask" },
                listed(&columns),
            ));
        }
    }
    if info.leaves_out(left_out::SOURCE) && config.output_format != "json" {
        let in_fields = matches!(config.output_format.as_str(), "tab" | "vcf")
            && config
                .fields
                .as_ref()
                .is_some_and(|fields| fields.iter().any(|field| field == "SOURCE"));
        // A source or `--merged` adds the column to an output that lists
        // none, and a field list writes it only where it names it.
        let custom = custom && written(config, "SOURCE");
        let merged = config.merged && written(config, "SOURCE");
        let added_by = ADD_A_CUSTOM_SOURCE.map(|flag| (flag, custom && named(config, flag)));
        let by_its_own_flag =
            custom && (named(config, "custom") || added_by.iter().all(|(_, on)| !on));
        let mut flags = vec![("custom", by_its_own_flag)];
        flags.extend(added_by);
        flags.extend([("merged", merged), ("fields", in_fields)]);
        let asked = asked_by(config, &flags);
        if !asked.is_empty() {
            notices.push(format!(
                "INFO: {} {} the SOURCE column, which is empty on every transcript of this \
                 cache: a cache built from a GTF or GFF3 read as the transcript model names no \
                 source",
                listed(&asked),
                if asked.len() == 1 { "writes" } else { "write" },
            ));
        }
    }
    if info.leaves_out(left_out::PANEL_FREQUENCY) && config.check_existing {
        let asked = asked_by(
            config,
            &[
                ("af", config.af && written(config, "AF")),
                ("af_1kg", config.af_1kg && written(config, "AF")),
                ("check_frequency", filter_reads_af(config)),
            ],
        );
        if !asked.is_empty() {
            notices.push(format!(
                "INFO: {} {} the AF of this cache, which is each allele's count over the 1000 \
                 Genomes panel, recovered from the continental frequencies, and can differ from \
                 the AF of Ensembl VEP's own cache \
                 (docs/cache-variation-data.md#the-1000-genomes-frequencies)",
                listed(&asked),
                if asked.len() == 1 { "reads" } else { "read" },
            ));
        }
    }
    notices.extend(cell_notices(config, info));
    if info.leaves_out(left_out::SEQ_EDITS) {
        notices.push(
            "INFO: this cache was converted from transcript files that hold no sequence edits \
             (seq_edits): a variant at codon 1 of a transcript with an alternative initiation \
             codon is called from the cached peptide; a conversion of Ensembl VEP's cache with \
             vep-storable-to-json carries the edits (docs/cache-conversion.md)"
                .to_string(),
        );
    }
    notices
}

/// The warning of a run that asks for known variants from a cache without
/// variation data: [`VariationSource::missing_warning`], then the flags of
/// the command that ask for them, what a frequency filter then does, and
/// the caches that hold the data. No frequency filter is built on a cache
/// whose record declares no variation columns, so `--check_frequency` and
/// `--filter_common` keep every record there. Where the record declares
/// them and the cache holds no variation file, the filter is built and reads
/// no known variant: `--freq_filter include` drops every record, and
/// `exclude` keeps every record.
pub(crate) fn no_variation_data_warning(
    config: &Config,
    info: &CacheInfo,
    cache_dir: &str,
) -> String {
    let asked = asked_by(
        config,
        &[
            ("af", config.af),
            ("af_1kg", config.af_1kg),
            ("af_gnomad", config.af_gnomad),
            ("af_gnomade", config.af_gnomade),
            ("af_gnomadg", config.af_gnomadg),
            ("max_af", config.max_af),
            ("pubmed", config.pubmed),
            ("var_synonyms", config.var_synonyms),
            ("check_frequency", config.check_frequency),
        ],
    );
    let asked = if asked.is_empty() {
        vec![asking(config, "check_existing")]
    } else {
        asked
    };
    let filter = if !config.check_frequency {
        ""
    } else if info.variation_columns().is_some() && config.freq_filter == "include" {
        ", and the frequency filter drops every record"
    } else {
        ", and the frequency filter keeps every record"
    };
    format!(
        "{} under {}{filter}; {VARIATION_DATA_HELD_BY}",
        VariationSource::missing_warning(cache_dir),
        listed(&asked),
    )
}

/// The line of a `vep recode` run whose fields read known variants over a
/// cache without variation data. The tool reads the cache's variation files
/// itself and its annotation run never asks for known variants, so that run
/// raises no [`no_variation_data_warning`]. `fields` are the output keys the
/// tool reads from those files, which no object of such a run holds; the
/// caches that hold the data are named as the warning names them.
pub(crate) fn no_variation_data_recode_notice(fields: &[&str], cache_dir: &str) -> String {
    format!(
        "INFO: no variation data under {cache_dir}: {}, read from a cache's known variants, {} \
         unwritten; {VARIATION_DATA_HELD_BY}",
        listed(fields),
        if fields.len() == 1 { "stays" } else { "stay" },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(args: &[&str]) -> Config {
        let argv = ["vep", "--offline", "-i", "in.vcf", "-o", "out.txt"]
            .iter()
            .chain(args)
            .map(|arg| arg.to_string());
        Config::from_argv(argv).expect("the command line parses")
    }

    fn leaving_out(names: &[&str]) -> CacheInfo {
        let mut info = CacheInfo::default();
        info.left_out = names.iter().map(|name| name.to_string()).collect();
        info
    }

    /// A record without the key draws no line, whatever the command asks
    /// for, and a record that names a part draws none from a command that
    /// asks for none of it.
    #[test]
    fn a_record_without_the_key_and_a_command_that_asks_for_nothing_draw_no_line() {
        let everything = config(&["--everything", "--xref_refseq"]);
        assert!(notices(&everything, &CacheInfo::default(), true).is_empty());
        let plain = config(&[]);
        let stated = leaving_out(&["core_tables", "source", "panel_frequency"]);
        assert!(notices(&plain, &stated, false).is_empty());
    }

    /// `--everything` asks for twelve of the thirteen columns of a build
    /// without the core table dumps, `--xref_refseq` for the thirteenth; the
    /// line names the flags the command gives and the columns in the order
    /// of the flags' groups.
    #[test]
    fn the_core_tables_line_names_the_flags_given_and_their_columns() {
        let info = leaving_out(&["core_tables"]);
        assert_eq!(
            notices(&config(&["--everything"]), &info, false),
            [
                "INFO: --everything asks for SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, \
              MANE_PLUS_CLINICAL, APPRIS, miRNA, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, \
              DOMAINS and GENE_PHENO, which this cache does not carry: it was built without the \
              core table dumps (vep-cache-builder --core-dumps), so the columns stay empty and \
              mature_miRNA_variant is never assigned"
            ]
        );
        assert_eq!(
            notices(&config(&["--xref_refseq"]), &info, false),
            [
                "INFO: --xref_refseq asks for RefSeq, which this cache does not carry: it was \
              built without the core table dumps (vep-cache-builder --core-dumps), so the column \
              stays empty"
            ]
        );
        assert_eq!(
            notices(
                &config(&["--uniprot", "--symbol", "--mane_select"]),
                &info,
                false
            ),
            [
                "INFO: --symbol, --mane_select and --uniprot ask for SYMBOL_SOURCE, HGNC_ID, \
              MANE_SELECT, SWISSPROT, TREMBL, UNIPARC and UNIPROT_ISOFORM, which this cache does \
              not carry: it was built without the core table dumps (vep-cache-builder \
              --core-dumps), so the columns stay empty"
            ]
        );
        assert_eq!(
            notices(&config(&["--vcf"]), &info, false),
            [
                "INFO: --vcf asks for SYMBOL_SOURCE and HGNC_ID, which this cache does not carry: \
              it was built without the core table dumps (vep-cache-builder --core-dumps), so the \
              columns stay empty"
            ]
        );
    }

    /// The `SOURCE` line is drawn by a `--custom` source, by `--merged` and
    /// by a `--fields` list of a tab output that names the column, and by no
    /// other command; a JSON output, which holds no such column, draws none.
    /// A source that `--bigwig` or `--gtf` adds is named by that flag.
    #[test]
    fn the_source_line_names_what_writes_the_column() {
        let info = leaving_out(&["source"]);
        let line = |flags: &str, verb: &str| {
            format!(
                "INFO: {flags} {verb} the SOURCE column, which is empty on every transcript of \
                 this cache: a cache built from a GTF or GFF3 read as the transcript model names \
                 no source"
            )
        };
        assert!(notices(&config(&["--symbol"]), &info, false).is_empty());
        assert_eq!(
            notices(&config(&[]), &info, true),
            [line("--custom", "writes")]
        );
        assert_eq!(
            notices(
                &config(&["--tab", "--fields", "Feature,SOURCE", "--merged"]),
                &info,
                false
            ),
            [line("--merged and --fields", "write")]
        );
        assert_eq!(
            notices(&config(&["--bigwig", "scores.bw"]), &info, true),
            [line("--bigwig", "writes")]
        );
        assert_eq!(
            notices(
                &config(&["--custom", "x.bed,X,bed", "--gtf", "genes.gtf.gz"]),
                &info,
                true
            ),
            [line("--custom and --gtf", "write")]
        );
        assert!(notices(&config(&["--fields", "Feature,SOURCE"]), &info, false).is_empty());
        assert!(notices(&config(&["--json", "--merged"]), &info, true).is_empty());
    }

    /// The recovered `AF` line is drawn by `--af`, by `--af_1kg`, by
    /// `--everything` and by a frequency filter on `1KG_ALL`, which
    /// `--filter_common` is; `--max_af` and a filter on another population
    /// read other columns and draw none.
    #[test]
    fn the_recovered_frequency_line_names_the_flags_that_read_af() {
        let info = leaving_out(&["panel_frequency"]);
        let line = |flags: &str, verb: &str| {
            format!(
                "INFO: {flags} {verb} the AF of this cache, which is each allele's count over \
                 the 1000 Genomes panel, recovered from the continental frequencies, and can \
                 differ from the AF of Ensembl VEP's own cache \
                 (docs/cache-variation-data.md#the-1000-genomes-frequencies)"
            )
        };
        assert_eq!(
            notices(&config(&["--af"]), &info, false),
            [line("--af", "reads")]
        );
        assert_eq!(
            notices(&config(&["--everything"]), &info, false),
            [line("--everything", "reads")]
        );
        assert_eq!(
            notices(&config(&["--af_1kg", "--filter_common"]), &info, false),
            [line("--af_1kg and --filter_common", "read")]
        );
        assert_eq!(
            notices(&config(&["--check_frequency"]), &info, false),
            [line("--check_frequency", "reads")]
        );
        assert!(notices(&config(&["--max_af"]), &info, false).is_empty());
        assert!(notices(
            &config(&["--check_frequency", "--freq_pop", "1KG_AFR"]),
            &info,
            false
        )
        .is_empty());
    }

    /// The sequence-edits line is drawn by every run on a cache whose record
    /// names them, and the lines of several names follow the order core
    /// tables, source, frequency, edits.
    #[test]
    fn the_sequence_edits_line_is_drawn_by_every_run() {
        let edits = "INFO: this cache was converted from transcript files that hold no sequence \
                     edits (seq_edits): a variant at codon 1 of a transcript with an alternative \
                     initiation codon is called from the cached peptide; a conversion of Ensembl \
                     VEP's cache with vep-storable-to-json carries the edits \
                     (docs/cache-conversion.md)";
        assert_eq!(
            notices(&config(&[]), &leaving_out(&["seq_edits"]), false),
            [edits]
        );
        let all = leaving_out(&["seq_edits", "panel_frequency", "source", "core_tables"]);
        let lines = notices(&config(&["--everything"]), &all, true);
        assert_eq!(lines.len(), 4, "{lines:?}");
        assert!(lines[0].starts_with("INFO: --everything asks for SYMBOL_SOURCE"));
        assert!(lines[1].starts_with("INFO: --custom writes the SOURCE column"));
        assert!(lines[2].starts_with("INFO: --everything reads the AF of this cache"));
        assert_eq!(lines[3], edits);
    }

    /// The columns of a cache's variation files: the fifteen every cache
    /// declares, then `frequencies`.
    fn with_variation(frequencies: &[&str], names: &[&str]) -> CacheInfo {
        let mut info = leaving_out(names);
        info.variation_cols = vep_core::cache_info::FIXED_VARIATION_COLUMNS
            .iter()
            .chain(frequencies)
            .map(|column| column.to_string())
            .collect();
        info
    }

    const PANEL: [&str; 6] = ["AF", "AFR", "AMR", "EAS", "EUR", "SAS"];
    const ROUTE: &str =
        "; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

    /// A cache the builder made for `homo_sapiens` GRCh38 declares the 1000
    /// Genomes columns alone and names the cells it holds no value for. An
    /// `--everything` run draws one line for every column that stays empty,
    /// the gnomAD columns read from `variation_cols` and the others from the
    /// names, with the route, and one line for `MAX_AF`; a cache whose
    /// record declares every column the command asks for and names nothing
    /// draws neither.
    #[test]
    fn a_built_cache_draws_one_line_for_its_empty_cells_and_one_for_max_af() {
        let mut built = with_variation(
            &PANEL,
            &[
                "gnomad_frequencies",
                "pubmed",
                "var_synonyms",
                "failed_variants",
                "cell_type_activity",
                "motif_scores",
            ],
        );
        built.regulatory = true;
        assert_eq!(
            notices(&config(&["--everything"]), &built, false),
            [
                format!(
                    "INFO: --everything asks for the gnomAD exome columns, the gnomAD genome \
                     columns, PUBMED, HIGH_INF_POS and MOTIF_SCORE_CHANGE, which this cache \
                     does not carry: the columns stay empty{ROUTE}"
                ),
                format!(
                    "INFO: MAX_AF and MAX_AF_POPS (--everything) range over the 1000 Genomes \
                     populations alone on this cache, whose variation data holds no gnomAD \
                     exome or gnomAD genome frequency{ROUTE}"
                ),
            ]
        );
        assert_eq!(
            notices(&config(&["--var_synonyms", "--pubmed"]), &built, false),
            [format!(
                "INFO: --pubmed and --var_synonyms ask for PUBMED and VAR_SYNONYMS, which this \
                 cache does not carry: the columns stay empty{ROUTE}"
            )]
        );
        assert_eq!(
            notices(&config(&["--regulatory"]), &built, false),
            [format!(
                "INFO: --regulatory asks for HIGH_INF_POS and MOTIF_SCORE_CHANGE, which this \
                 cache does not carry: the columns stay empty{ROUTE}"
            )]
        );
        assert!(notices(&config(&["--regulatory_gff", "reg.gff3.gz"]), &built, false).is_empty());
        assert!(notices(&config(&["--af", "--af_1kg", "--symbol"]), &built, false).is_empty());

        let mut gnomad: Vec<&str> = PANEL.to_vec();
        gnomad.extend(["gnomADe", "gnomADe_AFR", "gnomADg", "gnomADg_AFR"]);
        let mut converted = with_variation(&gnomad, &[]);
        converted.regulatory = true;
        assert!(notices(
            &config(&["--everything", "--var_synonyms"]),
            &converted,
            false
        )
        .is_empty());
    }

    /// A cache with known variants and no frequency column is a cache that
    /// holds no frequencies: each frequency flag draws the line with its own
    /// columns, `MAX_AF` and `MAX_AF_POPS` among them, and a frequency filter
    /// on `1KG_ALL` says what it then does to the records. The route is
    /// named only where the record names the frequencies as left out.
    #[test]
    fn a_cache_with_known_variants_and_no_frequency_column_draws_the_frequency_line() {
        let bare = with_variation(&[], &[]);
        assert_eq!(
            notices(&config(&["--af"]), &bare, false),
            ["INFO: --af asks for AF, which this cache does not carry: the column stays empty"]
        );
        assert_eq!(
            notices(
                &config(&["--af_1kg", "--max_af", "--af_gnomadg"]),
                &bare,
                false
            ),
            [
                "INFO: --af_1kg, --af_gnomadg and --max_af ask for AF, AFR_AF, AMR_AF, EAS_AF, \
              EUR_AF, SAS_AF, the gnomAD genome columns, MAX_AF and MAX_AF_POPS, which this \
              cache does not carry: the columns stay empty"
            ]
        );
        assert_eq!(
            notices(&config(&["--filter_common"]), &bare, false),
            [
                "INFO: --filter_common reads AF, which this cache does not carry: the frequency \
              filter finds no AF and keeps every record"
            ]
        );
        assert_eq!(
            notices(
                &config(&["--af", "--check_frequency", "--freq_filter", "include"]),
                &bare,
                false
            ),
            [
                "INFO: --af and --check_frequency ask for AF, which this cache does not carry: \
              the column stays empty and the frequency filter finds no AF and drops every \
              record"
            ]
        );
        let built = with_variation(&[], &["thousand_genomes_frequencies", "gnomad_frequencies"]);
        assert_eq!(
            notices(&config(&["--af", "--af_gnomade"]), &built, false),
            [format!(
                "INFO: --af and --af_gnomade ask for AF and the gnomAD exome columns, which \
                 this cache does not carry: the columns stay empty{ROUTE}"
            )]
        );
        // No variation columns at all: the run's warning says so, not this line.
        assert!(notices(&config(&["--everything"]), &CacheInfo::default(), false).is_empty());
    }

    /// `MAX_AF` on a cache that declares the populations of two of the three
    /// groups names the two it ranges over and the one it lacks, with no
    /// route where the record names nothing as left out.
    #[test]
    fn the_max_af_line_names_the_groups_the_record_declares() {
        let two = with_variation(&["AF", "AFR", "gnomADe", "gnomADe_AFR"], &[]);
        assert_eq!(
            notices(&config(&["--max_af"]), &two, false),
            [
                "INFO: MAX_AF and MAX_AF_POPS (--max_af) range over the 1000 Genomes and gnomAD \
              exome populations alone on this cache, whose variation data holds no gnomAD \
              genome frequency"
            ]
        );
    }

    /// A tab or a VCF output under `--fields` writes the columns its list
    /// names, so each line names those alone and the flags that ask for
    /// them, and a list that names none of a line's columns draws no line;
    /// a frequency filter reads `AF` whatever the list names.
    #[test]
    fn a_field_list_narrows_each_line_to_the_columns_it_writes() {
        let mut built = with_variation(
            &PANEL,
            &[
                "core_tables",
                "source",
                "panel_frequency",
                "gnomad_frequencies",
                "pubmed",
                "motif_scores",
            ],
        );
        built.regulatory = true;
        let tab = |fields: &str, flags: &[&str]| {
            let args = [&["--tab", "--fields", fields][..], flags].concat();
            notices(&config(&args), &built, true)
        };
        assert!(tab("Uploaded_variation,Consequence", &["--everything"]).is_empty());
        assert_eq!(
            tab(
                "Uploaded_variation,SWISSPROT,MAX_AF,gnomADe_AF,HIGH_INF_POS",
                &["--everything"]
            ),
            [
                "INFO: --everything asks for SWISSPROT, which this cache does not carry: it \
                 was built without the core table dumps (vep-cache-builder --core-dumps), so \
                 the column stays empty"
                    .to_string(),
                format!(
                    "INFO: --everything asks for the gnomAD exome columns and HIGH_INF_POS, \
                     which this cache does not carry: the columns stay empty{ROUTE}"
                ),
                format!(
                    "INFO: MAX_AF (--everything) ranges over the 1000 Genomes populations alone \
                     on this cache, whose variation data holds no gnomAD exome or gnomAD genome \
                     frequency{ROUTE}"
                ),
            ]
        );
        assert_eq!(
            tab("Uploaded_variation,AF", &["--af", "--symbol", "--pubmed"]),
            [
                "INFO: --af reads the AF of this cache, which is each allele's count over the \
              1000 Genomes panel, recovered from the continental frequencies, and can differ \
              from the AF of Ensembl VEP's own cache \
              (docs/cache-variation-data.md#the-1000-genomes-frequencies)"
            ]
        );
        assert_eq!(
            tab("Uploaded_variation", &["--filter_common"]),
            [
                "INFO: --filter_common reads the AF of this cache, which is each allele's count \
              over the 1000 Genomes panel, recovered from the continental frequencies, and can \
              differ from the AF of Ensembl VEP's own cache \
              (docs/cache-variation-data.md#the-1000-genomes-frequencies)"
            ]
        );
    }

    /// A VCF output under `--fields` is held to its list as a tab output
    /// is. The default format and JSON read no such list, so a command that
    /// gives one there draws each line whole.
    #[test]
    fn a_vcf_field_list_narrows_the_lines_and_another_format_reads_no_list() {
        let built = with_variation(&PANEL, &["core_tables", "gnomad_frequencies"]);
        let listing = |format: &[&str]| {
            let list = [
                "--fields",
                "Consequence,HGNC_ID,MAX_AF_POPS",
                "--everything",
            ];
            notices(&config(&[format, &list[..]].concat()), &built, false)
        };
        assert_eq!(
            listing(&["--vcf"]),
            [
                "INFO: --everything asks for HGNC_ID, which this cache does not carry: it was \
                 built without the core table dumps (vep-cache-builder --core-dumps), so the \
                 column stays empty"
                    .to_string(),
                format!(
                    "INFO: MAX_AF_POPS (--everything) ranges over the 1000 Genomes populations \
                     alone on this cache, whose variation data holds no gnomAD exome or gnomAD \
                     genome frequency{ROUTE}"
                ),
            ]
        );
        let whole = notices(&config(&["--everything"]), &built, false);
        assert!(
            whole[0].contains("SYMBOL_SOURCE, HGNC_ID, MANE_SELECT"),
            "{whole:?}"
        );
        for format in [&[][..], &["--json"]] {
            assert_eq!(listing(format), whole, "{format:?}");
        }
    }

    /// A frequency filter reads the cache's `AF` whatever a field list
    /// names: on a cache with known variants and no frequency column, a tab
    /// output whose list leaves `AF` out draws the filter's clause alone, and
    /// one whose list names `AF` draws the clause after the column.
    #[test]
    fn a_field_list_keeps_the_clause_of_a_frequency_filter() {
        let bare = with_variation(&[], &[]);
        let tab = |fields: &str| {
            let args = ["--tab", "--fields", fields, "--af", "--filter_common"];
            notices(&config(&args), &bare, false)
        };
        assert_eq!(
            tab("Uploaded_variation"),
            [
                "INFO: --filter_common reads AF, which this cache does not carry: the frequency \
              filter finds no AF and keeps every record"
            ]
        );
        assert_eq!(
            tab("Uploaded_variation,AF"),
            [
                "INFO: --af and --filter_common ask for AF, which this cache does not carry: the \
              column stays empty and the frequency filter finds no AF and keeps every record"
            ]
        );
    }

    /// The warning of a cache without variation data names the flags of the
    /// command that ask for known variants, says that a frequency filter
    /// keeps every record where one is asked for, and names the caches that
    /// hold the data. Where the record declares variation columns, a filter
    /// that includes by frequency drops every record, and the warning says
    /// that in place of the other.
    #[test]
    fn the_variation_warning_names_the_flags_the_filter_and_the_route() {
        let route = "; a cache from vep-cache-builder --auto or --include-variations, or a \
                     conversion of Ensembl VEP's cache, holds it (docs/cache-variation-data.md)";
        let opening = "WARNING: [cache] no variation data under /c: Existing_variation and the \
                       frequency columns stay empty under ";
        let undeclared = CacheInfo::default();
        assert_eq!(
            no_variation_data_warning(&config(&["--check_existing"]), &undeclared, "/c"),
            format!("{opening}--check_existing{route}")
        );
        assert_eq!(
            no_variation_data_warning(&config(&["--everything"]), &undeclared, "/c"),
            format!("{opening}--everything{route}")
        );
        assert_eq!(
            no_variation_data_warning(
                &config(&["--af", "--pubmed", "--var_synonyms"]),
                &undeclared,
                "/c"
            ),
            format!("{opening}--af, --pubmed and --var_synonyms{route}")
        );
        assert_eq!(
            no_variation_data_warning(
                &config(&["--everything", "--filter_common"]),
                &undeclared,
                "/c"
            ),
            format!(
                "{opening}--everything and --filter_common, and the frequency filter keeps every \
                 record{route}"
            )
        );
        assert_eq!(
            no_variation_data_warning(&config(&["--check_frequency"]), &undeclared, "/c"),
            format!(
                "{opening}--check_frequency, and the frequency filter keeps every record{route}"
            )
        );
        let include = config(&["--check_frequency", "--freq_filter", "include"]);
        assert_eq!(
            no_variation_data_warning(&include, &undeclared, "/c"),
            format!(
                "{opening}--check_frequency, and the frequency filter keeps every record{route}"
            )
        );
        let declared = with_variation(&[], &[]);
        assert_eq!(
            no_variation_data_warning(&include, &declared, "/c"),
            format!(
                "{opening}--check_frequency, and the frequency filter drops every record{route}"
            )
        );
        assert_eq!(
            no_variation_data_warning(&config(&["--filter_common"]), &declared, "/c"),
            format!("{opening}--filter_common, and the frequency filter keeps every record{route}")
        );
    }

    /// The line of `vep recode` names the fields that read known variants, one
    /// or both, says they stay unwritten and names the caches that hold the
    /// data in the warning's words.
    #[test]
    fn the_recode_line_names_the_fields_left_unwritten_and_the_route() {
        let route = "; a cache from vep-cache-builder --auto or --include-variations, or a \
                     conversion of Ensembl VEP's cache, holds it (docs/cache-variation-data.md)";
        assert_eq!(
            no_variation_data_recode_notice(&["id"], "/c"),
            format!(
                "INFO: no variation data under /c: id, read from a cache's known variants, \
                 stays unwritten{route}"
            )
        );
        assert_eq!(
            no_variation_data_recode_notice(&["var_synonyms"], "/c"),
            format!(
                "INFO: no variation data under /c: var_synonyms, read from a cache's known \
                 variants, stays unwritten{route}"
            )
        );
        assert_eq!(
            no_variation_data_recode_notice(&["id", "var_synonyms"], "/c"),
            format!(
                "INFO: no variation data under /c: id and var_synonyms, read from a cache's \
                 known variants, stay unwritten{route}"
            )
        );
    }
}
