// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The plugin registry: which plugin writes a field name.
//!
//! The tables are the rows of `data/plugin_fields.tsv`, which describes its row kinds: the
//! field sets with their entries, the sets that are not plugins, the names Ensembl VEP's
//! command line loads a plugin under, and the plugin columns that carry the name of one of
//! Ensembl VEP's own Extra keys.

use regex::Regex;
use std::collections::HashSet;
use std::sync::OnceLock;

const TABLE: &str = include_str!("../data/plugin_fields.tsv");

/// A field set of the registry: a plugin's unless the table's `not_a_plugin` rows name it.
pub struct FieldSet {
    pub name: String,
    pub entries: Vec<String>,
}

/// The lookups a run's sources are resolved through. `BuiltInRegistry` answers from the
/// tool's own table; a test plants another.
pub trait Registry: Sync {
    /// The registry's plugins, in the table's order.
    fn all_plugin_names(&self) -> &[String];
    /// The registry's name of a plugin that Ensembl VEP's command line loads under `name`,
    /// where the two differ.
    fn plugin_alias(&self, name: &str) -> Option<&str>;
    /// The one plugin of `plugins` that writes `field`: see [`plugin_of_field`].
    fn plugin_of_field(
        &self,
        field: &str,
        plugins: &[String],
        flags: Option<&[String]>,
    ) -> Result<Option<String>, String>;
}

/// The registry of `data/plugin_fields.tsv`.
pub struct BuiltInRegistry;

impl Registry for BuiltInRegistry {
    fn all_plugin_names(&self) -> &[String] {
        all_plugin_names()
    }

    fn plugin_alias(&self, name: &str) -> Option<&str> {
        plugin_alias(name)
    }

    fn plugin_of_field(
        &self,
        field: &str,
        plugins: &[String],
        flags: Option<&[String]>,
    ) -> Result<Option<String>, String> {
        plugin_of_field(field, plugins, flags)
    }
}

/// Tells whether a field name matches one of a list of registry entries: an entry without a
/// placeholder by equality, an entry with one as a pattern over the whole name, in which
/// `<N>` stands for one or more decimal digits and any other `<word>` for one or more
/// characters other than a line feed.
pub struct FieldMatcher {
    names: HashSet<String>,
    pattern: Option<Regex>,
}

impl FieldMatcher {
    pub fn new<S: AsRef<str>>(entries: &[S]) -> FieldMatcher {
        let mut names = HashSet::new();
        let mut patterns = Vec::new();
        for entry in entries {
            let entry = entry.as_ref();
            match entry_pattern(entry) {
                Some(pattern) => patterns.push(pattern),
                None => {
                    names.insert(entry.to_string());
                }
            }
        }
        let pattern = (!patterns.is_empty()).then(|| {
            Regex::new(&format!("^(?:{})$", patterns.join("|")))
                .expect("escaped literals and two fixed classes always compile")
        });
        FieldMatcher { names, pattern }
    }

    pub fn claims(&self, field: &str) -> bool {
        self.names.contains(field)
            || self
                .pattern
                .as_ref()
                .is_some_and(|pattern| pattern.is_match(field))
    }
}

/// The pattern of an entry that holds a placeholder, `None` for a plain name. A placeholder
/// is `<`, one or more word characters as CPython's `\w` reads text (letters, numbers and
/// the underscore) and `>`.
fn entry_pattern(entry: &str) -> Option<String> {
    static PLACEHOLDER: OnceLock<Regex> = OnceLock::new();
    let placeholder = PLACEHOLDER.get_or_init(|| {
        Regex::new(r"<([\p{L}\p{N}_]+)>").expect("the placeholder pattern compiles")
    });
    let mut pattern = String::new();
    let mut literal_from = 0;
    for found in placeholder.captures_iter(entry) {
        let (whole, name) = (found.get(0)?, found.get(1)?);
        pattern.push_str(&regex::escape(&entry[literal_from..whole.start()]));
        // CPython's `\d` over text is a decimal digit of any script; its `.` is any
        // character but a line feed.
        pattern.push_str(if name.as_str() == "N" {
            r"\p{Nd}+"
        } else {
            r"[^\n]+"
        });
        literal_from = whole.end();
    }
    if literal_from == 0 {
        return None;
    }
    pattern.push_str(&regex::escape(&entry[literal_from..]));
    Some(pattern)
}

/// A plugin's column that carries the name of one of Ensembl VEP's own Extra keys, with the
/// flags that turn that key on.
pub struct CoreNamedColumn {
    pub plugin: String,
    pub column: String,
    pub flags: Vec<String>,
}

struct Tables {
    sets: Vec<FieldSet>,
    plugin_names: Vec<String>,
    aliases: Vec<(String, String)>,
    core_named: Vec<CoreNamedColumn>,
    /// Per plugin of `plugin_names`: its entries and its core-named columns. The comparison
    /// reads the entries alone, where a core-named column would count a core key as the
    /// plugin's.
    matchers: Vec<FieldMatcher>,
}

fn tables() -> &'static Tables {
    static TABLES: OnceLock<Tables> = OnceLock::new();
    TABLES.get_or_init(|| read_table(TABLE))
}

fn read_table(text: &str) -> Tables {
    let mut sets: Vec<FieldSet> = Vec::new();
    let mut not_plugins: Vec<&str> = Vec::new();
    let mut aliases: Vec<(String, String)> = Vec::new();
    let mut core_named: Vec<CoreNamedColumn> = Vec::new();
    for (number, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        match cells[..] {
            ["field", set, entry] => match sets.iter_mut().find(|known| known.name == set) {
                Some(known) => known.entries.push(entry.to_string()),
                None => sets.push(FieldSet {
                    name: set.to_string(),
                    entries: vec![entry.to_string()],
                }),
            },
            ["not_a_plugin", set] => not_plugins.push(set),
            ["alias", alias, plugin] => aliases.push((alias.to_string(), plugin.to_string())),
            ["core_named", plugin, column, flag] => {
                match core_named
                    .iter_mut()
                    .find(|known| known.plugin == plugin && known.column == column)
                {
                    Some(known) => known.flags.push(flag.to_string()),
                    None => core_named.push(CoreNamedColumn {
                        plugin: plugin.to_string(),
                        column: column.to_string(),
                        flags: vec![flag.to_string()],
                    }),
                }
            }
            _ => panic!(
                "data/plugin_fields.tsv line {}: not a row of the registry",
                number + 1
            ),
        }
    }
    let plugins = || {
        sets.iter()
            .filter(|set| !not_plugins.contains(&set.name.as_str()))
    };
    let plugin_names: Vec<String> = plugins().map(|set| set.name.clone()).collect();
    let matchers = plugins()
        .map(|set| {
            let columns = core_named
                .iter()
                .filter(|known| known.plugin == set.name)
                .map(|known| known.column.as_str());
            let entries: Vec<&str> = set
                .entries
                .iter()
                .map(String::as_str)
                .chain(columns)
                .collect();
            FieldMatcher::new(&entries)
        })
        .collect();
    Tables {
        sets,
        plugin_names,
        aliases,
        core_named,
        matchers,
    }
}

/// Every field set of the registry with its entries, in the table's order.
pub fn field_sets() -> &'static [FieldSet] {
    &tables().sets
}

/// The plugin columns that carry the name of one of Ensembl VEP's own Extra keys, in the
/// table's order. No field set lists them.
pub fn core_named_columns() -> &'static [CoreNamedColumn] {
    &tables().core_named
}

/// The registry's plugins, in the table's order: its field sets without those that are not
/// plugins.
pub fn all_plugin_names() -> &'static [String] {
    &tables().plugin_names
}

/// The registry's name of a plugin that Ensembl VEP's command line loads under `name`,
/// where the two differ.
pub fn plugin_alias(name: &str) -> Option<&'static str> {
    tables()
        .aliases
        .iter()
        .find(|(alias, _)| alias == name)
        .map(|(_, plugin)| plugin.as_str())
}

/// The one plugin of a run that writes the Extra key `field`, or `None`.
///
/// `plugins` names the plugins the run loaded, as the registry names them or as Ensembl
/// VEP's command line does; the answer is the registry's name. A name two plugins write
/// belongs to the one of them the run loaded: with both loaded the answer is `None`, as it
/// is for a key none of the run's plugins writes. A key of a field set that is not a plugin
/// belongs to no plugin.
///
/// A plugin's column that carries the name of one of Ensembl VEP's own keys answers its
/// plugin when the run loaded it. `flags` holds the words of the run's command line: a run
/// that loaded the plugin and set a flag that turns the core key of that name on is refused
/// with `Err`, since the name is then both the column's and the key's. Without `flags` no
/// run is refused.
pub fn plugin_of_field(
    field: &str,
    plugins: &[String],
    flags: Option<&[String]>,
) -> Result<Option<String>, String> {
    let tables = tables();
    let loaded = |name: &String| {
        plugins.contains(name)
            || tables
                .aliases
                .iter()
                .any(|(alias, plugin)| plugin == name && plugins.contains(alias))
    };
    let mut writers = tables
        .plugin_names
        .iter()
        .zip(&tables.matchers)
        .filter(|(name, matcher)| matcher.claims(field) && loaded(name))
        .map(|(name, _)| name);
    let (Some(writer), None) = (writers.next(), writers.next()) else {
        return Ok(None);
    };
    if let Some(flags) = flags {
        let column = tables
            .core_named
            .iter()
            .find(|known| &known.plugin == writer && known.column == field);
        for flag in column.map_or(&[][..], |known| &known.flags[..]) {
            if flags.contains(flag) {
                return Err(format!(
                    "the run loads {writer} and sets {flag}: {field} is both the plugin's column and Ensembl VEP's own key"
                ));
            }
        }
    }
    Ok(Some(writer.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(names: &[&str]) -> Vec<String> {
        names.iter().map(|name| name.to_string()).collect()
    }

    fn of(field: &str, plugins: &[&str]) -> Option<String> {
        plugin_of_field(field, &words(plugins), None).expect("no run is refused without flags")
    }

    fn of_run(field: &str, plugins: &[&str], flags: &[&str]) -> Result<Option<String>, String> {
        plugin_of_field(field, &words(plugins), Some(&words(flags)))
    }

    fn every() -> Vec<&'static str> {
        all_plugin_names().iter().map(String::as_str).collect()
    }

    fn plugin(name: &str) -> Option<String> {
        Some(name.to_string())
    }

    #[test]
    fn the_table_gives_the_field_sets_and_the_plugins_in_its_order() {
        let sets: Vec<&str> = field_sets().iter().map(|set| set.name.as_str()).collect();
        assert_eq!(sets.len(), 32);
        assert_eq!(sets[..4], ["HGVS", "CADD", "REVEL", "SpliceAI"]);
        assert_eq!(sets[31], "UTRAnnotator");
        assert_eq!(every(), sets[1..]);
        assert_eq!(field_sets()[0].entries, ["HGVSc", "HGVSp"]);
        let entries: usize = field_sets().iter().map(|set| set.entries.len()).sum();
        assert_eq!(entries, 641);
        let dbnsfp = field_sets()
            .iter()
            .find(|set| set.name == "dbNSFP")
            .unwrap();
        assert_eq!(dbnsfp.entries.len(), 458);
        assert_eq!(dbnsfp.entries[0], "dbNSFP_<column>");
        assert_eq!(dbnsfp.entries[1], "1000Gp3_AC");
        assert_eq!(dbnsfp.entries[456], "rs_dbSNP");
        assert_eq!(dbnsfp.entries[457], "NoSuchColumn");
        assert!(!dbnsfp
            .entries
            .iter()
            .any(|entry| entry == "APPRIS" || entry == "TSL"));
    }

    #[test]
    fn plugin_of_field_answers_with_the_plugin_of_the_run_that_writes_a_field() {
        let every = every();
        assert_eq!(of("gnomAD_30x_cov", &["gnomADc"]), plugin("gnomADc"));
        assert_eq!(of("gnomAD_30x_cov", &["REVEL", "gnomADMt"]), None);
        assert_eq!(of("gnomAD_median_approx_cov", &every), plugin("gnomADc"));
        assert_eq!(of("gnomAD_x_cov", &every), None);
        assert_eq!(of("gnomAD_meanx_cov", &every), None);
        assert_eq!(of("gnomAD_30x_cov_het", &every), None);
        assert_eq!(of("am_gnomAD_30x_cov", &every), None);
        assert_eq!(of("gnomADc_mean", &every), None);
        assert_eq!(of("gnomAD_AN", &every), plugin("gnomADMt"));
        assert_eq!(of("gnomAD_exomes_AF", &every), plugin("dbNSFP"));
        assert_eq!(of("dbNSFP_SIFT_score", &every), plugin("dbNSFP"));
        assert_eq!(of("dbNSFP_", &every), None);
        assert_eq!(
            of("transcript_annotator_Comment", &every),
            plugin("TranscriptAnnotator")
        );
        assert_eq!(of("IMPACT", &every), None);
        let with_the_field_set: Vec<&str> = every.iter().copied().chain(["HGVS"]).collect();
        assert_eq!(of("HGVSc", &with_the_field_set), None);
        assert_eq!(of("LoF_filter", &["LoF"]), plugin("LoFTEE"));
        assert_eq!(of("LoF_filter", &["LoFTEE"]), plugin("LoFTEE"));
        assert_eq!(
            of("LoF_filter", &["REVEL", "LoF", "LoFTEE", "LoF"]),
            plugin("LoFTEE")
        );
        assert_eq!(of("REVEL", &["REVEL", "REVEL"]), plugin("REVEL"));
        assert_eq!(of("REVEL", &["CADD", "REVEL"]), plugin("REVEL"));
        assert_eq!(of("REVEL", &["REVEL", "CADD"]), plugin("REVEL"));
    }

    #[test]
    fn a_name_two_plugins_write_is_settled_by_the_run_and_never_guessed() {
        for name in ["SIFT_pred", "SIFT_score"] {
            assert_eq!(of(name, &["dbNSFP"]), plugin("dbNSFP"));
            assert_eq!(
                of(name, &["REVEL", "TranscriptAnnotator"]),
                plugin("TranscriptAnnotator")
            );
            assert_eq!(of(name, &["dbNSFP", "TranscriptAnnotator"]), None);
            assert_eq!(of(name, &["TranscriptAnnotator", "dbNSFP"]), None);
            assert_eq!(of(name, &every()), None);
            assert_eq!(of(name, &["REVEL"]), None);
            assert_eq!(of(name, &[]), None);
        }
    }

    #[test]
    fn a_column_one_of_two_loaded_writers_alone_names_is_that_writers() {
        assert_eq!(
            of("Comment", &["TranscriptAnnotator"]),
            plugin("TranscriptAnnotator")
        );
        assert_eq!(
            of("Comment", &["TranscriptAnnotator", "dbNSFP"]),
            plugin("TranscriptAnnotator")
        );
        assert_eq!(of("Comment", &["dbNSFP"]), None);
    }

    #[test]
    fn field_matcher_reads_several_pattern_entries() {
        let claims = FieldMatcher::new(&["exact", "a_<N>", "b_<name>_c"]);
        for field in ["exact", "a_12", "b_x_c"] {
            assert!(claims.claims(field), "{field}");
        }
        for field in ["exactly", "a_", "a_1x", "a_12b_x_c", "b__c", ""] {
            assert!(!claims.claims(field), "{field}");
        }
    }

    #[test]
    fn field_matcher_reads_placeholders_and_literal_text_by_cpython_pattern_rules() {
        // A digit run is decimal digits of any script, and no other numeric character.
        let digits = FieldMatcher::new(&["d_<N>"]);
        assert!(digits.claims("d_\u{0663}\u{0664}"));
        assert!(!digits.claims("d_\u{00b2}"));
        // A name is any characters but a line feed.
        let name = FieldMatcher::new(&["n_<column>"]);
        assert!(name.claims("n_a b\tc\r"));
        assert!(!name.claims("n_a\nb"));
        // Angle brackets around anything but word characters are the entry's own text.
        let literal = FieldMatcher::new(&["a<>", "b<x y>", "c<x", "d.<N>"]);
        for field in ["a<>", "b<x y>", "c<x", "d.7"] {
            assert!(literal.claims(field), "{field}");
        }
        for field in ["a<1>", "b<x>", "dx7", "d.<N>"] {
            assert!(!literal.claims(field), "{field}");
        }
        // Two placeholders in one entry, and one that is the whole entry.
        let several = FieldMatcher::new(&["<N>x<N>", "<word>"]);
        assert!(several.claims("3x4") && several.claims("anything"));
        assert!(!FieldMatcher::new(&["<N>x<N>"]).claims("3x"));
        assert!(!FieldMatcher::new(&["<N>"]).claims(""));
        // A placeholder's name is letters, numbers and underscores of any script.
        let named = FieldMatcher::new(&["p_<\u{e9}t\u{e9}_\u{00b2}>", "q_<a-b>"]);
        assert!(named.claims("p_x") && named.claims("q_<a-b>"));
        assert!(!named.claims("q_x"));
    }

    #[test]
    fn a_run_that_loads_the_plugin_and_turns_the_core_key_on_is_refused() {
        assert_eq!(of("APPRIS", &["dbNSFP"]), plugin("dbNSFP"));
        assert_eq!(
            of_run("TSL", &["REVEL", "dbNSFP"], &["--offline", "--appris"]),
            Ok(plugin("dbNSFP"))
        );
        assert_eq!(
            of_run("APPRIS", &["dbNSFP"], &["--tsl"]),
            Ok(plugin("dbNSFP"))
        );
        assert_eq!(of_run("APPRIS", &["REVEL"], &["--everything"]), Ok(None));
        assert_eq!(
            of_run("SIFT_score", &["dbNSFP"], &["--everything"]),
            Ok(plugin("dbNSFP"))
        );
        for (key, flag) in [
            ("APPRIS", "--appris"),
            ("APPRIS", "--everything"),
            ("TSL", "--tsl"),
            ("TSL", "--everything"),
        ] {
            assert_eq!(
                of_run(key, &["dbNSFP"], &["--offline", flag]),
                Err(format!(
                    "the run loads dbNSFP and sets {flag}: {key} is both the plugin's column and Ensembl VEP's own key"
                ))
            );
        }
    }

    #[test]
    fn the_refusal_names_the_first_listed_flag_the_run_sets_and_reads_whole_words() {
        assert_eq!(
            of_run("APPRIS", &["dbNSFP"], &["--everything", "--appris"]),
            Err("the run loads dbNSFP and sets --appris: APPRIS is both the plugin's column and Ensembl VEP's own key".to_string())
        );
        assert_eq!(
            of_run(
                "TSL",
                &["dbNSFP"],
                &["vep", "--tsl=1", "--plugin", "dbNSFP,synthetic.gz,--tsl"]
            ),
            Ok(plugin("dbNSFP"))
        );
        assert_eq!(of_run("TSL", &["dbNSFP"], &[]), Ok(plugin("dbNSFP")));
        assert_eq!(
            of_run("TSL", &["LoF", "dbNSFP"], &["--everything"]),
            Err("the run loads dbNSFP and sets --everything: TSL is both the plugin's column and Ensembl VEP's own key".to_string())
        );
    }

    #[test]
    fn a_core_named_column_is_its_plugins_only_in_a_run_that_loads_the_plugin() {
        for key in ["APPRIS", "TSL"] {
            assert_eq!(of(key, &["dbNSFP", "LOEUF"]), plugin("dbNSFP"));
            assert_eq!(of(key, &["LOEUF"]), None);
            assert_eq!(of(key, &[]), None);
            let others: Vec<&str> = every()
                .into_iter()
                .filter(|name| *name != "dbNSFP")
                .collect();
            assert_eq!(of(key, &others), None);
            assert_eq!(of_run(key, &others, &["--everything"]), Ok(None));
        }
    }

    #[test]
    fn a_plugin_is_named_by_its_registry_name_or_an_alias_the_registry_lists() {
        assert_eq!(plugin_alias("LoF"), Some("LoFTEE"));
        assert_eq!(plugin_alias("LoFTEE"), None);
        assert_eq!(plugin_alias("lof"), None);
        assert_eq!(plugin_alias(""), None);
        assert!(!all_plugin_names().iter().any(|name| name == "LoF"));
        assert_eq!(of("LoF", &["LoF"]), plugin("LoFTEE"));
        assert_eq!(of("LoF", &["LoF", "LoFTEE"]), plugin("LoFTEE"));
        assert_eq!(of("LoF_info", &["REVEL"]), None);
    }

    #[test]
    fn the_answer_is_never_a_plugin_the_run_did_not_load() {
        for set in field_sets() {
            for entry in set.entries.iter().filter(|entry| !entry.contains('<')) {
                for loaded in [&["REVEL"][..], &["LOEUF", "CADD"][..]] {
                    let answer = of(entry, loaded);
                    assert!(
                        answer.as_deref().is_none_or(|name| loaded.contains(&name)),
                        "{entry} under {loaded:?} answered {answer:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_built_in_registry_answers_as_the_functions_do() {
        let registry: &dyn Registry = &BuiltInRegistry;
        assert_eq!(registry.all_plugin_names(), all_plugin_names());
        assert_eq!(registry.plugin_alias("LoF"), Some("LoFTEE"));
        assert_eq!(registry.plugin_alias("REVEL"), None);
        assert_eq!(
            registry.plugin_of_field("gnomAD_30x_cov", &words(&["gnomADc"]), None),
            Ok(plugin("gnomADc"))
        );
        assert!(registry
            .plugin_of_field("TSL", &words(&["dbNSFP"]), Some(&words(&["--tsl"])))
            .is_err());
    }

    #[test]
    fn a_planted_registry_stands_in_through_the_trait() {
        struct Planted(Vec<String>);
        impl Registry for Planted {
            fn all_plugin_names(&self) -> &[String] {
                &self.0
            }
            fn plugin_alias(&self, name: &str) -> Option<&str> {
                (name == "Loaded").then_some("Planted")
            }
            fn plugin_of_field(
                &self,
                field: &str,
                plugins: &[String],
                flags: Option<&[String]>,
            ) -> Result<Option<String>, String> {
                if flags.is_some_and(|flags| flags.iter().any(|flag| flag == "--planted")) {
                    return Err("refused".to_string());
                }
                Ok(
                    (field.starts_with("PLANT_") && plugins.contains(&self.0[0]))
                        .then(|| self.0[0].clone()),
                )
            }
        }
        let registry: &dyn Registry = &Planted(words(&["Planted"]));
        assert_eq!(registry.all_plugin_names(), ["Planted"]);
        assert_eq!(registry.plugin_alias("Loaded"), Some("Planted"));
        assert_eq!(
            registry.plugin_of_field("PLANT_score", &words(&["Planted"]), None),
            Ok(plugin("Planted"))
        );
        assert_eq!(
            registry.plugin_of_field(
                "PLANT_score",
                &words(&["Planted"]),
                Some(&words(&["--planted"]))
            ),
            Err("refused".to_string())
        );
    }

    /// One recorded answer: a line of JSON with the keys `key`, `plugins`, `flags` and
    /// `answer`, the answer `null`, `{"plugin": name}` or `{"error": message}`.
    struct Recorded {
        key: String,
        plugins: Vec<String>,
        flags: Option<Vec<String>>,
        answer: Result<Option<String>, String>,
    }

    /// A reader of the JSON the recorded lines are written in: objects, arrays, strings and
    /// `null`.
    struct Json<'a> {
        text: &'a [u8],
        at: usize,
    }

    enum Value {
        Null,
        Text(String),
        List(Vec<Value>),
        Object(Vec<(String, Value)>),
    }

    impl Json<'_> {
        fn peek(&mut self) -> u8 {
            while self.text[self.at].is_ascii_whitespace() {
                self.at += 1;
            }
            self.text[self.at]
        }

        fn take(&mut self, byte: u8) {
            assert_eq!(self.peek(), byte, "at byte {}", self.at);
            self.at += 1;
        }

        fn unit(&mut self) -> u32 {
            let digits = std::str::from_utf8(&self.text[self.at..self.at + 4]).unwrap();
            self.at += 4;
            u32::from_str_radix(digits, 16).unwrap()
        }

        fn string(&mut self) -> String {
            self.take(b'"');
            let mut bytes = Vec::new();
            loop {
                let byte = self.text[self.at];
                self.at += 1;
                match byte {
                    b'"' => break,
                    b'\\' => {
                        let escape = self.text[self.at];
                        self.at += 1;
                        let character = match escape {
                            b'n' => '\n',
                            b't' => '\t',
                            b'r' => '\r',
                            b'b' => '\u{8}',
                            b'f' => '\u{c}',
                            b'u' => {
                                let mut code = self.unit();
                                if (0xd800..0xdc00).contains(&code) {
                                    assert_eq!(&self.text[self.at..self.at + 2], b"\\u");
                                    self.at += 2;
                                    code =
                                        0x10000 + ((code - 0xd800) << 10) + (self.unit() - 0xdc00);
                                }
                                char::from_u32(code).unwrap()
                            }
                            other => other as char,
                        };
                        bytes.extend_from_slice(character.encode_utf8(&mut [0; 4]).as_bytes());
                    }
                    other => bytes.push(other),
                }
            }
            String::from_utf8(bytes).unwrap()
        }

        fn value(&mut self) -> Value {
            match self.peek() {
                b'n' => {
                    assert_eq!(&self.text[self.at..self.at + 4], b"null");
                    self.at += 4;
                    Value::Null
                }
                b'"' => Value::Text(self.string()),
                b'[' => {
                    self.take(b'[');
                    let mut items = Vec::new();
                    while self.peek() != b']' {
                        items.push(self.value());
                        if self.peek() == b',' {
                            self.take(b',');
                        }
                    }
                    self.take(b']');
                    Value::List(items)
                }
                _ => {
                    self.take(b'{');
                    let mut members = Vec::new();
                    while self.peek() != b'}' {
                        let name = self.string();
                        self.take(b':');
                        members.push((name, self.value()));
                        if self.peek() == b',' {
                            self.take(b',');
                        }
                    }
                    self.take(b'}');
                    Value::Object(members)
                }
            }
        }
    }

    fn texts(value: Value) -> Vec<String> {
        let Value::List(items) = value else {
            panic!("a list of strings");
        };
        items
            .into_iter()
            .map(|item| match item {
                Value::Text(text) => text,
                _ => panic!("a list of strings"),
            })
            .collect()
    }

    fn recorded(line: &str) -> Recorded {
        let mut json = Json {
            text: line.as_bytes(),
            at: 0,
        };
        let Value::Object(members) = json.value() else {
            panic!("a recorded answer is an object");
        };
        let (mut key, mut plugins, mut flags, mut answer) = (None, None, None, None);
        for (name, value) in members {
            match (name.as_str(), value) {
                ("key", Value::Text(text)) => key = Some(text),
                ("plugins", value) => plugins = Some(texts(value)),
                ("flags", Value::Null) => flags = Some(None),
                ("flags", value) => flags = Some(Some(texts(value))),
                ("answer", Value::Null) => answer = Some(Ok(None)),
                ("answer", Value::Object(mut members)) => {
                    let (kind, Value::Text(text)) = members.remove(0) else {
                        panic!("an answer names a plugin or an error");
                    };
                    answer = Some(match kind.as_str() {
                        "plugin" => Ok(Some(text)),
                        "error" => Err(text),
                        _ => panic!("an answer names a plugin or an error"),
                    });
                }
                _ => {}
            }
        }
        Recorded {
            key: key.unwrap(),
            plugins: plugins.unwrap(),
            flags: flags.unwrap(),
            answer: answer.unwrap(),
        }
    }

    #[test]
    fn the_reader_of_recorded_answers_reads_each_shape() {
        let read = recorded(
            r#"{"command": "x", "key": "a\u00e9\ud83d\ude00\"\\b\t", "plugins": ["P", "Q"], "flags": ["--f"], "answer": {"plugin": "P"}}"#,
        );
        assert_eq!(read.key, "a\u{e9}\u{1f600}\"\\b\t");
        assert_eq!(read.plugins, ["P", "Q"]);
        assert_eq!(read.flags, Some(words(&["--f"])));
        assert_eq!(read.answer, Ok(plugin("P")));
        let read = recorded(r#"{"key": "k", "plugins": [], "flags": null, "answer": null}"#);
        assert_eq!(
            (read.plugins.len(), read.flags, read.answer),
            (0, None, Ok(None))
        );
        let read =
            recorded(r#"{"key": "k", "plugins": ["P"], "flags": [], "answer": {"error": "no"}}"#);
        assert_eq!(read.flags, Some(Vec::new()));
        assert_eq!(read.answer, Err("no".to_string()));
    }

    /// Replays recorded answers: every line of the file `VEP_PLUGIN_PAIRS_REGISTRY_ANSWERS`
    /// names must be answered as recorded. Without the variable nothing is replayed.
    #[test]
    fn recorded_answers_are_given_again() {
        let Some(path) = std::env::var_os("VEP_PLUGIN_PAIRS_REGISTRY_ANSWERS") else {
            return;
        };
        let text = std::fs::read_to_string(&path).expect("the recorded answers are readable");
        let mut replayed = 0;
        let mut differing = Vec::new();
        for (number, line) in text.lines().enumerate() {
            if line.is_empty() {
                continue;
            }
            let read = recorded(line);
            let answer = plugin_of_field(&read.key, &read.plugins, read.flags.as_deref());
            if answer != read.answer {
                differing.push(format!(
                    "line {}: {:?} under {:?}: recorded {:?}, answered {:?}",
                    number + 1,
                    read.key,
                    read.plugins,
                    read.answer,
                    answer
                ));
            }
            replayed += 1;
        }
        println!(
            "replayed {replayed} recorded answers, {} differ",
            differing.len()
        );
        assert!(replayed > 0, "the file of recorded answers is empty");
        assert!(differing.is_empty(), "{}", differing.join("\n"));
    }
}
