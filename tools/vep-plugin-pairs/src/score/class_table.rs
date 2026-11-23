// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The documented classes a run declares in its class table, and the ones it applies.

use super::error::{Scored, ScorerError};
use super::sources::{Sources, RESOLVED_LIMIT};
use crate::pytext::{repr, splitlines, strip};
use crate::rules::{Rules, ORDER_ONLY};
use sha2::{Digest, Sha256};
use std::collections::{BTreeSet, HashMap};
use std::fmt::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};

/// What the fourth column of a class table's row reads where a command line applies the
/// class by its name alone.
pub const NOT_MET: &str = "not met";

/// One class of the table: its name and, per plugin in the row's order, the field patterns
/// it reads there.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Declared {
    pub name: String,
    pub plugins: Vec<(String, Vec<String>)>,
    /// The row's columns after the third, carried into the report unread but for the first.
    pub further: Vec<String>,
}

impl Declared {
    fn patterns(&self, source: &str) -> Option<&[String]> {
        self.plugins
            .iter()
            .find(|(plugin, _)| plugin == source)
            .map(|(_, patterns)| patterns.as_slice())
    }
}

/// The class table a run read: its path as given and the digest of its bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TableFacts {
    pub path: String,
    pub sha256: String,
}

/// The classes that read one field of one source, in the order they are consulted.
pub type Reading = Arc<[String]>;

/// The documented classes a run declares, in the order they are consulted, and the ones the
/// run applies.
pub struct Classes {
    declared: Vec<Declared>,
    table: Option<TableFacts>,
    applied: BTreeSet<String>,
    populations: Vec<(String, u64)>,
    reading: Mutex<HashMap<(String, String), Reading>>,
}

/// Whether the whole of `field` has the shape of a class table's field pattern, `*` its
/// only wildcard: any run of characters without a line feed.
pub fn field_shape_matches(pattern: &str, field: &str) -> bool {
    let mut parts = pattern.split('*');
    let first = parts.next().unwrap_or("");
    if !field.starts_with(first) {
        return false;
    }
    let mut reached: BTreeSet<usize> = BTreeSet::from([first.len()]);
    let mut wildcard = false;
    for part in parts {
        wildcard = true;
        let mut next = BTreeSet::new();
        if let Some(&earliest) = reached.first() {
            // A run from any place reached runs on from the earliest one unless a line
            // feed lies between them.
            let mut open = BTreeSet::new();
            for (offset, character) in field[earliest..]
                .char_indices()
                .map(|(offset, character)| (earliest + offset, Some(character)))
                .chain([(field.len(), None)])
            {
                if reached.contains(&offset) {
                    open.insert(offset);
                }
                if !open.is_empty() && field[offset..].starts_with(part) {
                    next.insert(offset + part.len());
                }
                if character == Some('\n') {
                    open.clear();
                }
            }
        }
        reached = next;
    }
    if wildcard {
        reached.contains(&field.len())
    } else {
        first.len() == field.len()
    }
}

impl Classes {
    /// The classes of the table `table`, whose bytes are `table_bytes` or, without them,
    /// the file's; no class without a table. `applied` names the classes the run applies:
    /// every declared class where it is `None`, less with `met_alone` those whose fourth
    /// column reads `not met`. `populations` gives a class the keys it can reach. A row that
    /// is malformed, names a class `rules` does not know, a plugin that is no source or no
    /// field of it, and an applied or populated class the table does not declare are
    /// refusals.
    pub fn new(
        table: Option<&Path>,
        table_bytes: Option<&[u8]>,
        applied: Option<&[String]>,
        populations: &[(String, u64)],
        sources: &Sources,
        met_alone: bool,
        rules: &dyn Rules,
    ) -> Scored<Classes> {
        let mut declared: Vec<Declared> = Vec::new();
        let mut facts = None;
        let mut read = Vec::new();
        let mut shown = String::new();
        if let Some(path) = table {
            shown = path.to_string_lossy().into_owned();
            read = match table_bytes {
                Some(bytes) => bytes.to_vec(),
                None => std::fs::read(path).map_err(|error| ScorerError::os(&error, path))?,
            };
            let digest = Sha256::digest(&read);
            let mut sha256 = String::with_capacity(64);
            for byte in digest.iter() {
                write!(sha256, "{byte:02x}").unwrap();
            }
            facts = Some(TableFacts {
                path: shown.clone(),
                sha256,
            });
        }
        let lines = match std::str::from_utf8(&read) {
            Ok(text) => splitlines(text),
            Err(_) => {
                return Err(ScorerError(format!(
                    "{shown}: the class table is not UTF-8 text"
                )))
            }
        };
        let known = |plugin: &str| {
            sources.registry_names().iter().any(|name| name == plugin)
                || sources.custom().iter().any(|(name, _)| name == plugin)
        };
        for (number, line) in lines.iter().enumerate() {
            let number = number + 1;
            if strip(line).is_empty() || line.starts_with('#') {
                continue;
            }
            let columns: Vec<&str> = line.split('\t').map(strip).collect();
            if columns.len() < 3 {
                return Err(ScorerError(format!(
                    "{shown}: line {number} lacks a class, its plugins and their field names"
                )));
            }
            let name = columns[0];
            let plugins: Vec<&str> = columns[1].split(';').collect();
            let groups: Vec<&str> = columns[2].split(';').collect();
            if name != ORDER_ONLY && !rules.knows(name) {
                return Err(ScorerError(format!(
                    "{shown}: line {number} names {name}, no class of compare_vep_formats.py"
                )));
            }
            if declared.iter().any(|class| class.name == name) {
                return Err(ScorerError(format!(
                    "{shown}: line {number} declares {name} a second time"
                )));
            }
            if groups.len() != plugins.len() {
                return Err(ScorerError(format!(
                    "{shown}: line {number} gives {} plugins and {} field lists",
                    plugins.len(),
                    groups.len()
                )));
            }
            let mut class = Declared {
                name: name.to_string(),
                plugins: Vec::new(),
                further: columns[3..]
                    .iter()
                    .map(|column| column.to_string())
                    .collect(),
            };
            for (plugin, group) in plugins.iter().zip(&groups) {
                let plugin = strip(plugin);
                let patterns: Vec<String> = group
                    .split(',')
                    .map(strip)
                    .filter(|pattern| !pattern.is_empty())
                    .map(str::to_string)
                    .collect();
                if !known(plugin) || patterns.is_empty() {
                    return Err(ScorerError(format!(
                        "{shown}: line {number} names {}, no source, or no field of it",
                        repr(plugin)
                    )));
                }
                if class.patterns(plugin).is_some() {
                    return Err(ScorerError(format!(
                        "{shown}: line {number} names {} twice",
                        repr(plugin)
                    )));
                }
                class.plugins.push((plugin.to_string(), patterns));
            }
            declared.push(class);
        }
        let named = applied
            .unwrap_or(&[])
            .iter()
            .map(|name| ("--classes", name));
        let populated = populations
            .iter()
            .map(|(name, _)| ("--class-population", name));
        for (option, name) in named.chain(populated) {
            if !declared.iter().any(|class| &class.name == name) {
                return Err(ScorerError(format!(
                    "{option} names {name}, a class the class table does not declare"
                )));
            }
        }
        let applied = match applied {
            Some(names) => names.iter().cloned().collect(),
            None => declared
                .iter()
                .filter(|class| {
                    !(met_alone && class.further.first().map(String::as_str) == Some(NOT_MET))
                })
                .map(|class| class.name.clone())
                .collect(),
        };
        let mut held: Vec<(String, u64)> = Vec::new();
        for (name, population) in populations {
            match held.iter_mut().find(|(other, _)| other == name) {
                Some(entry) => entry.1 = *population,
                None => held.push((name.clone(), *population)),
            }
        }
        Ok(Classes {
            declared,
            table: facts,
            applied,
            populations: held,
            reading: Mutex::new(HashMap::new()),
        })
    }

    /// The declared classes in table order.
    pub fn declared(&self) -> &[Declared] {
        &self.declared
    }

    /// The columns after the third of the class's row; none for a class not declared.
    pub fn further(&self, name: &str) -> &[String] {
        self.declared
            .iter()
            .find(|class| class.name == name)
            .map_or(&[], |class| class.further.as_slice())
    }

    /// The table the run read; `None` for a run without one.
    pub fn table(&self) -> Option<&TableFacts> {
        self.table.as_ref()
    }

    /// The classes the run applies.
    pub fn applied(&self) -> &BTreeSet<String> {
        &self.applied
    }

    /// The class's place in the table.
    pub fn order(&self, name: &str) -> Option<usize> {
        self.declared.iter().position(|class| class.name == name)
    }

    /// The keys each class can reach, as the run gave them, a class once.
    pub fn populations(&self) -> &[(String, u64)] {
        &self.populations
    }

    /// The applied classes of `source`, in table order.
    pub fn of(&self, source: &str) -> Vec<String> {
        self.declared
            .iter()
            .filter(|class| class.patterns(source).is_some() && self.applied.contains(&class.name))
            .map(|class| class.name.clone())
            .collect()
    }

    /// The declared classes the table gives `plugin`.
    pub fn of_plugin(&self, plugin: &str) -> BTreeSet<String> {
        self.declared
            .iter()
            .filter(|class| class.patterns(plugin).is_some())
            .map(|class| class.name.clone())
            .collect()
    }

    /// The classes of `source` that read `field`: the applied ones in table order, then the
    /// others.
    pub fn reading(&self, source: &str, field: &str) -> Reading {
        let key = (source.to_string(), field.to_string());
        let mut reading = self.reading.lock().unwrap();
        if let Some(found) = reading.get(&key) {
            return Arc::clone(found);
        }
        let names: Vec<&String> = self
            .declared
            .iter()
            .filter(|class| {
                class.patterns(source).is_some_and(|patterns| {
                    patterns
                        .iter()
                        .any(|pattern| field_shape_matches(pattern, field))
                })
            })
            .map(|class| &class.name)
            .collect();
        let (applied, others): (Vec<&String>, Vec<&String>) = names
            .into_iter()
            .partition(|name| self.applied.contains(*name));
        let found: Reading = applied.into_iter().chain(others).cloned().collect();
        if reading.len() >= RESOLVED_LIMIT {
            reading.clear();
        }
        reading.insert(key, Arc::clone(&found));
        found
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::registry::Registry;
    use crate::rules::Consult;
    use crate::rules::Refusal;

    struct PlantedRegistry(Vec<String>);

    impl PlantedRegistry {
        fn new() -> &'static PlantedRegistry {
            let names = ["LOEUF", "GWAS", "MaxEntScan", "NearestExonJB", "REVEL"];
            Box::leak(Box::new(PlantedRegistry(
                names.iter().map(|name| name.to_string()).collect(),
            )))
        }
    }

    impl Registry for PlantedRegistry {
        fn all_plugin_names(&self) -> &[String] {
            &self.0
        }

        fn plugin_alias(&self, _name: &str) -> Option<&str> {
            None
        }

        fn plugin_of_field(
            &self,
            _field: &str,
            _plugins: &[String],
            _flags: Option<&[String]>,
        ) -> Result<Option<String>, String> {
            Ok(None)
        }
    }

    struct PlantedRules;

    impl Rules for PlantedRules {
        fn knows(&self, klass: &str) -> bool {
            klass.starts_with("planted_")
        }

        fn accepts(&self, _consult: &Consult) -> Result<bool, Refusal> {
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

    const TABLE: &str = "#class\tplugin\tfield names\tmet or not met\tcommit\n\
        order_only\tNearestExonJB\tNearestExonJB\tmet\t-\n\
        planted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n\
        planted_two_numbers\tMaxEntScan\tMaxEntScan_*\tmet\t-\n\
        planted_reference_value\tGWAS\tGWAS_*\tmet\t-\n\
        \n\
        planted_sibling_value\tREVEL\tREVEL\tnot met\t-\n";

    fn strings(words: &[&str]) -> Vec<String> {
        words.iter().map(|word| word.to_string()).collect()
    }

    fn sources() -> Sources<'static> {
        let custom = [("Custom".to_string(), "Custom".to_string())];
        Sources::new(&strings(&["LOEUF"]), &custom, None, PlantedRegistry::new()).unwrap()
    }

    fn classes(table: &str, applied: Option<&[&str]>, met_alone: bool) -> Result<Classes, String> {
        let applied = applied.map(strings);
        Classes::new(
            Some(Path::new("classes.tsv")),
            Some(table.as_bytes()),
            applied.as_deref(),
            &[],
            &sources(),
            met_alone,
            &PlantedRules,
        )
        .map_err(|error| error.0)
    }

    #[test]
    fn malformed_class_table_is_refused() {
        for (table, message) in [
            (
                "planted_engine_zero\tLOEUF\n",
                "classes.tsv: line 1 lacks a class, its plugins and their field names",
            ),
            (
                "#header\nno_such_class\tLOEUF\tLOEUF\n",
                "classes.tsv: line 2 names no_such_class, no class of compare_vep_formats.py",
            ),
            (
                "planted_engine_zero\tLOEUF\tLOEUF\nplanted_engine_zero\tGWAS\tGWAS_*\n",
                "classes.tsv: line 2 declares planted_engine_zero a second time",
            ),
            (
                "planted_engine_zero\tLOEUF;GWAS\tLOEUF\n",
                "classes.tsv: line 1 gives 2 plugins and 1 field lists",
            ),
            (
                "planted_engine_zero\tNoSuchPlugin\tLOEUF\n",
                "classes.tsv: line 1 names 'NoSuchPlugin', no source, or no field of it",
            ),
            (
                "planted_engine_zero\tLOEUF\t ,\n",
                "classes.tsv: line 1 names 'LOEUF', no source, or no field of it",
            ),
            (
                "planted_engine_zero\tLOEUF; LOEUF\tLOEUF;LOEUF_b\n",
                "classes.tsv: line 1 names 'LOEUF' twice",
            ),
        ] {
            assert_eq!(classes(table, None, false).err().as_deref(), Some(message));
        }
    }

    #[test]
    fn class_table_that_is_not_text_is_refused() {
        let refused = Classes::new(
            Some(Path::new("some/classes.tsv")),
            Some(b"planted_engine_zero\tLOEUF\tLOEUF\xff\n"),
            None,
            &[],
            &sources(),
            false,
            &PlantedRules,
        );
        assert_eq!(
            refused.err(),
            Some(ScorerError(
                "some/classes.tsv: the class table is not UTF-8 text".to_string()
            ))
        );
    }

    #[test]
    fn class_table_is_read_from_its_file_and_digested() {
        let directory = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-class-table-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("classes.tsv");
        std::fs::write(&path, "# no class\n").unwrap();
        let read = Classes::new(
            Some(&path),
            None,
            None,
            &[],
            &sources(),
            false,
            &PlantedRules,
        )
        .unwrap();
        std::fs::remove_dir_all(&directory).unwrap();
        assert!(read.declared().is_empty());
        assert_eq!(
            read.table(),
            Some(&TableFacts {
                path: path.to_string_lossy().into_owned(),
                sha256: "d4bc6f7b0a3f6328e2772c73fc5a33f23631afb480336d3d4f92151f6a85ed84"
                    .to_string()
            })
        );
        let missing = Classes::new(
            Some(&path),
            None,
            None,
            &[],
            &sources(),
            false,
            &PlantedRules,
        );
        assert_eq!(
            missing.err().map(|error| error.0),
            Some(format!(
                "[Errno 2] No such file or directory: '{}'",
                path.display()
            ))
        );
    }

    #[test]
    fn run_without_a_table_declares_no_class() {
        let none = Classes::new(None, None, None, &[], &sources(), false, &PlantedRules).unwrap();
        assert!(none.declared().is_empty() && none.applied().is_empty());
        assert_eq!(none.table(), None);
        assert!(none.reading("LOEUF", "LOEUF").is_empty());
    }

    #[test]
    fn classes_are_declared_in_table_order_with_their_further_columns() {
        let read = classes(TABLE, None, false).unwrap();
        let names: Vec<&str> = read
            .declared()
            .iter()
            .map(|class| class.name.as_str())
            .collect();
        assert_eq!(
            names,
            [
                "order_only",
                "planted_engine_zero",
                "planted_two_numbers",
                "planted_reference_value",
                "planted_sibling_value"
            ]
        );
        assert_eq!(
            read.further("planted_sibling_value"),
            strings(&["not met", "-"])
        );
        assert_eq!(read.further("order_only"), strings(&["met", "-"]));
        assert_eq!(read.order("planted_two_numbers"), Some(2));
        assert_eq!(read.order("no_such_class"), None);
        assert_eq!(read.applied().len(), 5);
        assert_eq!(
            read.declared()[2].plugins,
            vec![("MaxEntScan".to_string(), strings(&["MaxEntScan_*"]))]
        );
    }

    #[test]
    fn run_applies_the_classes_it_names_or_every_class_met() {
        let named = classes(TABLE, Some(&["planted_engine_zero"]), false).unwrap();
        assert_eq!(
            named.applied().iter().collect::<Vec<_>>(),
            ["planted_engine_zero"]
        );
        let none = classes(TABLE, Some(&[]), false).unwrap();
        assert!(none.applied().is_empty());
        let met = classes(TABLE, None, true).unwrap();
        assert!(!met.applied().contains("planted_sibling_value"));
        assert_eq!(met.applied().len(), 4);
        let by_name = classes(TABLE, Some(&["planted_sibling_value"]), true).unwrap();
        assert_eq!(
            by_name.applied().iter().collect::<Vec<_>>(),
            ["planted_sibling_value"]
        );
    }

    #[test]
    fn applied_or_populated_class_the_table_lacks_is_refused() {
        assert_eq!(
            classes(
                TABLE,
                Some(&["planted_engine_zero", "planted_other"]),
                false
            )
            .err()
            .as_deref(),
            Some("--classes names planted_other, a class the class table does not declare")
        );
        let populated = Classes::new(
            Some(Path::new("classes.tsv")),
            Some(TABLE.as_bytes()),
            None,
            &[
                ("planted_engine_zero".to_string(), 3),
                ("planted_other".to_string(), 1),
            ],
            &sources(),
            false,
            &PlantedRules,
        );
        assert_eq!(
            populated.err().map(|error| error.0).as_deref(),
            Some(
                "--class-population names planted_other, a class the class table does not declare"
            )
        );
    }

    #[test]
    fn class_of_two_plugins_reads_each_plugins_own_fields() {
        let table = "planted_engine_zero\tLOEUF;GWAS\tLOEUF;GWAS_pmid\n";
        let read = classes(table, None, false).unwrap();
        assert_eq!(
            read.declared()[0].plugins,
            vec![
                ("LOEUF".to_string(), strings(&["LOEUF"])),
                ("GWAS".to_string(), strings(&["GWAS_pmid"]))
            ]
        );
        assert!(read.further("planted_engine_zero").is_empty());
        assert_eq!(
            &*read.reading("LOEUF", "LOEUF"),
            strings(&["planted_engine_zero"])
        );
        assert_eq!(
            &*read.reading("GWAS", "GWAS_pmid"),
            strings(&["planted_engine_zero"])
        );
        assert!(read.reading("GWAS", "LOEUF").is_empty());
        assert_eq!(read.of("GWAS"), strings(&["planted_engine_zero"]));
        assert!(read.of("REVEL").is_empty());
        assert_eq!(read.of_plugin("LOEUF").len(), 1);
    }

    #[test]
    fn custom_source_can_carry_a_class() {
        let read = classes("planted_engine_zero\tCustom\tCustom_*\n", None, false).unwrap();
        assert_eq!(
            &*read.reading("Custom", "Custom_a"),
            strings(&["planted_engine_zero"])
        );
    }

    #[test]
    fn applied_classes_are_read_before_the_others_in_table_order() {
        let table = "planted_a\tLOEUF\tLOEUF\nplanted_b\tLOEUF\tLOEUF*\nplanted_c\tLOEUF\tLOEUF\n";
        let read = classes(table, Some(&["planted_c", "planted_b"]), false).unwrap();
        assert_eq!(
            &*read.reading("LOEUF", "LOEUF"),
            strings(&["planted_b", "planted_c", "planted_a"])
        );
        assert_eq!(&*read.reading("LOEUF", "LOEUF_b"), strings(&["planted_b"]));
        assert_eq!(read.of("LOEUF"), strings(&["planted_b", "planted_c"]));
    }

    #[test]
    fn star_is_the_class_tables_only_wildcard() {
        for (pattern, field, matches) in [
            ("LOEU?", "LOEU?", true),
            ("LOEU?", "LOEUF", false),
            ("[A-Z]OEUF", "LOEUF", false),
            ("[A-Z]OEUF", "[A-Z]OEUF", true),
            ("LO*F_b", "LOEUF_b", true),
            ("LO*F_b", "LOF_b", true),
            ("LO*F_b", "LOEUF_bc", false),
            ("L.EUF", "LOEUF", false),
            ("L.EUF", "L.EUF", true),
            ("MaxEntScan_*", "MaxEntScan_diff", true),
            ("MaxEntScan_*", "MaxEntScan_", true),
            ("MaxEntScan_*", "MaxEntScan", false),
            ("*", "", true),
            ("*", "anything", true),
            ("*", "two\nlines", false),
            ("a*b*c", "abc", true),
            ("a*b*c", "a1b2b3c", true),
            ("a*b*c", "a1b2b3", false),
            ("a*a", "a", false),
            ("a*a", "aa", true),
            ("a**a", "aa", true),
            ("a*b", "a\nb", false),
            ("a*b*c", "ab\nb1c", false),
            ("a*b*c", "a\nab1c", false),
            ("*é", "caf\u{e9}", true),
            ("", "", true),
            ("", "a", false),
        ] {
            assert_eq!(
                field_shape_matches(pattern, field),
                matches,
                "{pattern:?} on {field:?}"
            );
        }
    }
    #[test]
    fn committed_class_table_is_taken_whole_under_the_tools_own_classes_and_registry() {
        use crate::registry::BuiltInRegistry;
        use crate::rules::BuiltIn;
        let table = include_str!("../../../../scripts/concordance/plugin_class_table.tsv");
        let rows: Vec<Vec<&str>> = table
            .lines()
            .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
            .map(|line| line.split('\t').collect())
            .collect();
        let registry = BuiltInRegistry;
        let every = registry.all_plugin_names().to_vec();
        let sources = Sources::new(&every, &[], None, &registry).unwrap();
        let read = Classes::new(
            Some(Path::new("plugin_class_table.tsv")),
            Some(table.as_bytes()),
            None,
            &[],
            &sources,
            true,
            &BuiltIn,
        )
        .unwrap();
        let names: Vec<&str> = read
            .declared()
            .iter()
            .map(|class| class.name.as_str())
            .collect();
        assert_eq!(names, rows.iter().map(|cells| cells[0]).collect::<Vec<_>>());
        assert_eq!(names.len(), 47);
        for (class, cells) in read.declared().iter().zip(&rows) {
            let plugins: Vec<(String, Vec<String>)> = cells[1]
                .split(';')
                .zip(cells[2].split(';'))
                .map(|(plugin, group)| {
                    (
                        plugin.to_string(),
                        strings(&group.split(',').collect::<Vec<_>>()),
                    )
                })
                .collect();
            assert_eq!(class.plugins, plugins, "{}", class.name);
            assert_eq!(class.further, strings(&cells[3..]), "{}", class.name);
            assert_eq!(
                read.applied().contains(&class.name),
                cells[3] != NOT_MET,
                "{}",
                class.name
            );
        }
    }
}
