// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The sources of one scored run: the plugins it loaded and its custom sources.

use super::error::{Scored, ScorerError};
use crate::registry::Registry;
use crate::rules::fields::is_core_column;
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;

/// The field names remembered with their source before the memory is emptied.
pub const RESOLVED_LIMIT: usize = 100_000;

/// Bytes as text, a byte sequence that is not UTF-8 read as the replacement character.
pub fn text(value: &[u8]) -> String {
    String::from_utf8_lossy(value).into_owned()
}

/// A chromosome's name without a `chr` prefix of any case, the form in which an input
/// record's chromosome is held to a row's.
pub fn without_chr(chromosome: &str) -> &str {
    match chromosome.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &chromosome[3..],
        _ => chromosome,
    }
}

/// What a field name belongs to.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SourceOf {
    /// A column of Ensembl VEP's own.
    Core,
    /// A key no source of the run claims.
    Unclaimed,
    /// A plugin of the run or one of its custom sources.
    Named(String),
}

/// The spellings of the flag that names a run's assembly.
pub const ASSEMBLY_FLAGS: [&str; 2] = ["--assembly", "-a"];
/// The flag that sets a run's assembly to `GRCH37`.
pub const GRCH37_FLAG: &str = "--grch37";
/// The assembly `GRCH37_FLAG` sets.
pub const GRCH37: &str = "GRCh37";

/// The sources of one scored run. A field belongs to the plugin the registry returns among
/// the run's plugins, else to the first custom source whose prefix it starts with: one
/// source at most, and never a plugin the run did not load.
pub struct Sources<'r> {
    registry: &'r dyn Registry,
    registry_names: Vec<String>,
    plugins: Vec<String>,
    custom: Vec<(String, String)>,
    names: Vec<String>,
    flags: Option<Vec<String>>,
    core: Option<HashSet<String>>,
    resolved: Mutex<HashMap<Vec<u8>, SourceOf>>,
}

impl<'r> Sources<'r> {
    /// The sources of a run that loaded `plugins`, named as the registry names them or by
    /// an alias, and that groups the fields starting with a prefix under each `(name,
    /// prefix)` of `custom`. `run_words`, the words of the run's command line, are handed
    /// to the registry as its flags. A plugin the registry lacks, a custom source named as a
    /// registry plugin and a run without any source are refusals.
    pub fn new(
        plugins: &[String],
        custom: &[(String, String)],
        run_words: Option<&[String]>,
        registry: &'r dyn Registry,
    ) -> Scored<Sources<'r>> {
        let registry_names = registry.all_plugin_names().to_vec();
        let named: Vec<String> = plugins
            .iter()
            .map(|name| registry.plugin_alias(name).unwrap_or(name).to_string())
            .collect();
        let unknown: Vec<&str> = named
            .iter()
            .filter(|name| !registry_names.contains(name))
            .map(String::as_str)
            .collect();
        if !unknown.is_empty() {
            return Err(ScorerError(format!(
                "--plugins names {}, no plugin of the registry",
                unknown.join(", ")
            )));
        }
        let mut sources: Vec<(String, String)> = Vec::new();
        for (name, prefix) in custom {
            if registry_names.contains(name) {
                return Err(ScorerError(format!(
                    "--custom-source {name} names a plugin of the registry"
                )));
            }
            match sources.iter_mut().find(|(held, _)| held == name) {
                Some(held) => held.1 = prefix.clone(),
                None => sources.push((name.clone(), prefix.clone())),
            }
        }
        if named.is_empty() && sources.is_empty() {
            return Err(ScorerError(
                "no source to score: name the run's plugins (--plugins) or a --custom-source"
                    .to_string(),
            ));
        }
        let mut loaded: Vec<String> = Vec::new();
        for name in named {
            if !loaded.contains(&name) {
                loaded.push(name);
            }
        }
        let names = loaded
            .iter()
            .chain(sources.iter().map(|(name, _)| name))
            .cloned()
            .collect();
        Ok(Sources {
            registry,
            registry_names,
            plugins: loaded,
            custom: sources,
            names,
            flags: run_words.map(<[String]>::to_vec),
            core: None,
            resolved: Mutex::new(HashMap::new()),
        })
    }

    /// The same sources with `columns`, lower-case names, as the columns of Ensembl VEP's
    /// own in place of the built-in list.
    pub fn with_core_columns(mut self, columns: &[&str]) -> Sources<'r> {
        self.core = Some(columns.iter().map(|name| name.to_string()).collect());
        self.resolved = Mutex::new(HashMap::new());
        self
    }

    /// The plugins the run loaded by their registry names, each once, in the run's order.
    pub fn plugins(&self) -> &[String] {
        &self.plugins
    }

    /// The run's custom sources, each with its prefix.
    pub fn custom(&self) -> &[(String, String)] {
        &self.custom
    }

    /// Every source of the run: its plugins, then its custom sources.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Every plugin name of the registry.
    pub fn registry_names(&self) -> &[String] {
        &self.registry_names
    }

    /// The assembly the run's command line names: GRCh37 under `--grch37`, which sets the
    /// assembly over an `--assembly` word, else the value of the last `--assembly` or `-a`
    /// word, given as the next word or after `=`. `None` for a run that gave no command
    /// line and for one whose words name no assembly.
    pub fn assembly(&self) -> Option<&str> {
        let words = self.flags.as_deref()?;
        if words.iter().any(|word| word == GRCH37_FLAG) {
            return Some(GRCH37);
        }
        let mut named = None;
        for (at, word) in words.iter().enumerate() {
            for flag in ASSEMBLY_FLAGS {
                if word == flag {
                    named = words.get(at + 1).map(String::as_str).or(named);
                } else if let Some(value) = word.strip_prefix(flag) {
                    named = value.strip_prefix('=').or(named);
                }
            }
        }
        named
    }

    /// The field's source. A field the registry refuses for the run's command line is a
    /// refusal.
    pub fn of(&self, field: &[u8]) -> Scored<SourceOf> {
        if let Some(found) = self.resolved.lock().unwrap().get(field) {
            return Ok(found.clone());
        }
        let name = text(field);
        let mut claimed = if self.plugins.is_empty() {
            None
        } else {
            self.registry
                .plugin_of_field(&name, &self.plugins, self.flags.as_deref())
                .map_err(ScorerError)?
        };
        if claimed.is_none() {
            claimed = self
                .custom
                .iter()
                .find(|(_, prefix)| name.starts_with(prefix.as_str()))
                .map(|(custom, _)| custom.clone());
        }
        let source = match claimed {
            Some(source) => SourceOf::Named(source),
            None if self.core_column(&name) => SourceOf::Core,
            None => SourceOf::Unclaimed,
        };
        let mut resolved = self.resolved.lock().unwrap();
        if resolved.len() >= RESOLVED_LIMIT {
            resolved.clear();
        }
        resolved.insert(field.to_vec(), source.clone());
        Ok(source)
    }

    /// The source a tuple's field belongs to. The passes hand over tuples of claimed
    /// fields alone, so any other field is a refusal.
    pub fn named(&self, field: &[u8]) -> Scored<String> {
        match self.of(field)? {
            SourceOf::Named(source) => Ok(source),
            SourceOf::Core | SourceOf::Unclaimed => Err(ScorerError(format!(
                "the field {} of a tuple belongs to no source of the run",
                text(field)
            ))),
        }
    }

    /// Whether the name is a column of Ensembl VEP's own in any case.
    fn core_column(&self, name: &str) -> bool {
        match &self.core {
            Some(columns) => columns.contains(&name.to_lowercase()),
            None => is_core_column(name),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One question put to the registry: the field, the run's plugins and its flags.
    type Asked = (String, Vec<String>, Option<Vec<String>>);

    /// A registry of planted rows: a field belongs to the one plugin of the run that lists it.
    struct Planted {
        names: Vec<String>,
        rows: Vec<(&'static str, Vec<&'static str>)>,
        aliases: Vec<(&'static str, &'static str)>,
        asked: Mutex<Vec<Asked>>,
    }

    impl Planted {
        fn new(rows: Vec<(&'static str, Vec<&'static str>)>) -> Planted {
            Planted {
                names: rows.iter().map(|(name, _)| name.to_string()).collect(),
                rows,
                aliases: Vec::new(),
                asked: Mutex::new(Vec::new()),
            }
        }
    }

    impl Registry for Planted {
        fn all_plugin_names(&self) -> &[String] {
            &self.names
        }

        fn plugin_alias(&self, name: &str) -> Option<&str> {
            self.aliases
                .iter()
                .find(|(alias, _)| *alias == name)
                .map(|(_, plugin)| *plugin)
        }

        fn plugin_of_field(
            &self,
            field: &str,
            plugins: &[String],
            flags: Option<&[String]>,
        ) -> Result<Option<String>, String> {
            self.asked.lock().unwrap().push((
                field.to_string(),
                plugins.to_vec(),
                flags.map(<[String]>::to_vec),
            ));
            if field == "REFUSED" {
                return Err("the run's flags make REFUSED a key of Ensembl VEP's own".to_string());
            }
            let writers: Vec<&str> = self
                .rows
                .iter()
                .filter(|(name, fields)| {
                    fields.contains(&field) && plugins.iter().any(|plugin| plugin == name)
                })
                .map(|(name, _)| *name)
                .collect();
            Ok(match writers.as_slice() {
                [one] => Some(one.to_string()),
                _ => None,
            })
        }
    }

    fn strings(words: &[&str]) -> Vec<String> {
        words.iter().map(|word| word.to_string()).collect()
    }

    fn pairs(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(name, prefix)| (name.to_string(), prefix.to_string()))
            .collect()
    }

    fn planted() -> Planted {
        Planted::new(vec![
            ("Planted", vec!["PLANT_score", "PLANT_label"]),
            ("LOEUF", vec!["LOEUF"]),
            ("dbNSFP", vec!["SIFT_score", "SIFT_pred"]),
            ("TranscriptAnnotator", vec!["SIFT_score", "Comment"]),
        ])
    }

    fn named(source: &str) -> SourceOf {
        SourceOf::Named(source.to_string())
    }

    #[test]
    fn field_two_plugins_write_belongs_to_the_runs_own_plugin() {
        let registry = planted();
        let one = Sources::new(&strings(&["TranscriptAnnotator"]), &[], None, &registry).unwrap();
        assert_eq!(one.of(b"SIFT_score").unwrap(), named("TranscriptAnnotator"));
        assert_eq!(one.of(b"Comment").unwrap(), named("TranscriptAnnotator"));
        let other = Sources::new(&strings(&["dbNSFP"]), &[], None, &registry).unwrap();
        assert_eq!(other.of(b"SIFT_score").unwrap(), named("dbNSFP"));
        assert_eq!(other.of(b"Comment").unwrap(), SourceOf::Unclaimed);
        let both = Sources::new(
            &strings(&["TranscriptAnnotator", "dbNSFP"]),
            &[],
            None,
            &registry,
        )
        .unwrap();
        assert_eq!(both.of(b"SIFT_score").unwrap(), SourceOf::Unclaimed);
        assert_eq!(both.of(b"Comment").unwrap(), named("TranscriptAnnotator"));
    }

    #[test]
    fn field_is_asked_of_the_registry_among_the_runs_plugins_alone() {
        let registry = planted();
        let sources = Sources::new(
            &strings(&["Planted", "LOEUF", "Planted"]),
            &[],
            None,
            &registry,
        )
        .unwrap();
        assert_eq!(sources.plugins(), strings(&["Planted", "LOEUF"]));
        assert_eq!(sources.names(), strings(&["Planted", "LOEUF"]));
        assert_eq!(sources.of(b"PLANT_score").unwrap(), named("Planted"));
        assert_eq!(sources.of(b"PLANT_score").unwrap(), named("Planted"));
        assert_eq!(sources.of(b"depth_x").unwrap(), SourceOf::Unclaimed);
        assert_eq!(
            *registry.asked.lock().unwrap(),
            vec![
                (
                    "PLANT_score".to_string(),
                    strings(&["Planted", "LOEUF"]),
                    None
                ),
                ("depth_x".to_string(), strings(&["Planted", "LOEUF"]), None),
            ]
        );
    }

    #[test]
    fn run_words_reach_the_registry_as_its_flags_and_its_refusal_is_the_runs() {
        let registry = planted();
        let words = strings(&["vep", "--af"]);
        let sources = Sources::new(&strings(&["Planted"]), &[], Some(&words), &registry).unwrap();
        assert_eq!(sources.of(b"PLANT_label").unwrap(), named("Planted"));
        assert_eq!(registry.asked.lock().unwrap()[0].2, Some(words));
        assert_eq!(
            sources.of(b"REFUSED"),
            Err(ScorerError(
                "the run's flags make REFUSED a key of Ensembl VEP's own".to_string()
            ))
        );
    }

    #[test]
    fn run_names_its_sources() {
        let registry = planted();
        let refused = |plugins: &[&str], custom: &[(&str, &str)]| {
            Sources::new(&strings(plugins), &pairs(custom), None, &registry)
                .err()
                .map(|error| error.0)
        };
        assert_eq!(
            refused(&["Planted", "NoSuchPlugin", "Neither"], &[]).as_deref(),
            Some("--plugins names NoSuchPlugin, Neither, no plugin of the registry")
        );
        assert_eq!(
            refused(&[], &[]).as_deref(),
            Some("no source to score: name the run's plugins (--plugins) or a --custom-source")
        );
        assert_eq!(
            refused(&["Planted"], &[("LOEUF", "LOEUF")]).as_deref(),
            Some("--custom-source LOEUF names a plugin of the registry")
        );
        let custom = Sources::new(&[], &pairs(&[("Custom", "Custom")]), None, &registry).unwrap();
        assert!(custom.plugins().is_empty());
        assert_eq!(custom.names(), strings(&["Custom"]));
        assert_eq!(custom.of(b"Custom_a").unwrap(), named("Custom"));
        assert_eq!(custom.of(b"PLANT_score").unwrap(), SourceOf::Unclaimed);
    }

    #[test]
    fn plugin_named_by_its_alias_is_the_registrys_plugin() {
        let mut registry = planted();
        registry.aliases.push(("Seeded", "Planted"));
        let sources = Sources::new(&strings(&["Seeded", "Planted"]), &[], None, &registry).unwrap();
        assert_eq!(sources.plugins(), strings(&["Planted"]));
    }

    #[test]
    fn custom_source_claims_only_what_no_plugin_of_the_run_claims() {
        let registry = planted();
        let custom = pairs(&[("Wide", "PLANT_"), ("Wider", "PLANT")]);
        let sources = Sources::new(&strings(&["Planted"]), &custom, None, &registry).unwrap();
        assert_eq!(sources.names(), strings(&["Planted", "Wide", "Wider"]));
        assert_eq!(sources.of(b"PLANT_score").unwrap(), named("Planted"));
        assert_eq!(sources.of(b"PLANT_other").unwrap(), named("Wide"));
        assert_eq!(sources.of(b"PLANTED").unwrap(), named("Wider"));
        assert_eq!(sources.of(b"Mystery").unwrap(), SourceOf::Unclaimed);
    }

    #[test]
    fn column_of_ensembl_veps_own_is_no_unclaimed_key_in_any_case() {
        let registry = planted();
        let sources = Sources::new(&strings(&["Planted"]), &[], None, &registry).unwrap();
        assert_eq!(sources.of(b"IMPACT").unwrap(), SourceOf::Core);
        assert_eq!(sources.of(b"impact").unwrap(), SourceOf::Core);
        assert_eq!(sources.of(b"gnomADg_AF").unwrap(), SourceOf::Core);
        assert_eq!(sources.of(b"Mystery").unwrap(), SourceOf::Unclaimed);
        assert_eq!(
            sources.named(b"IMPACT"),
            Err(ScorerError(
                "the field IMPACT of a tuple belongs to no source of the run".to_string()
            ))
        );
        let planted = sources.with_core_columns(&["mystery"]);
        assert_eq!(planted.of(b"Mystery").unwrap(), SourceOf::Core);
        assert_eq!(planted.of(b"IMPACT").unwrap(), SourceOf::Unclaimed);
    }

    #[test]
    fn field_that_is_no_text_is_asked_by_its_replacement_text() {
        let registry = planted();
        let sources = Sources::new(&strings(&["Planted"]), &[], None, &registry).unwrap();
        assert_eq!(sources.of(b"PLANT\xff").unwrap(), SourceOf::Unclaimed);
        assert_eq!(registry.asked.lock().unwrap()[0].0, "PLANT\u{fffd}");
    }

    #[test]
    fn chromosome_loses_a_chr_prefix_of_any_case() {
        assert_eq!(without_chr("chr21"), "21");
        assert_eq!(without_chr("CHRX"), "X");
        assert_eq!(without_chr("ChR"), "");
        assert_eq!(without_chr("21"), "21");
        assert_eq!(without_chr("ch"), "ch");
        assert_eq!(without_chr("çhr1"), "çhr1");
    }
}
