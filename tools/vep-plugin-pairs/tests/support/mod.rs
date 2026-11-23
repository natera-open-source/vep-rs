// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the `score` program share: a scratch directory with small outputs,
//! tuple files and input VCFs in it, a planted plugin registry, planted classes that keep
//! every consult they were put, and the runs of `score` and of its command line over them.

#![allow(dead_code)]

use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use vep_plugin_pairs::pytext;
use vep_plugin_pairs::registry::Registry;
use vep_plugin_pairs::rules::{blank, BuiltIn, Consult, Entry, Refusal, Row, Rules};
use vep_plugin_pairs::score::{self, ClassTable, Options, Stop};

/// The column header of a default-format output, with its line end.
pub const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
/// The header lines of an input VCF.
pub const VCF_HEADER: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
/// What every error line of the program opens with.
pub const ERROR: &str = "ERROR: [score_plugin_tuples] ";
/// The count keys of a report entry: the raw three, then the adjusted three.
pub const COUNT_KEYS: [&str; 6] = [
    "perl",
    "rust",
    "intersection",
    "adj_perl",
    "adj_rust",
    "adj_intersection",
];
/// The raw count keys of a report entry.
pub const RAW_KEYS: [&str; 3] = ["perl", "rust", "intersection"];
/// The planted registry's rows: each plugin with the fields it writes.
pub const PLANTED: [(&str, &[&str]); 8] = [
    ("Planted", &["PLANT_score", "PLANT_label"]),
    ("PlantedDynamic", &["PlantedDynamic_a", "PlantedDynamic_b"]),
    ("LOEUF", &["LOEUF"]),
    ("MaxEntScan", &["MaxEntScan_diff"]),
    ("GWAS", &["GWAS_pmid"]),
    ("NearestExonJB", &["NearestExonJB"]),
    ("UTRAnnotator", &["5UTR_consequence"]),
    ("REVEL", &["REVEL"]),
];
/// The planted classes, each with the predicate `PlantedRules::planted` gives it.
pub const PLANTED_CLASSES: [&str; 5] = [
    "planted_engine_zero",
    "planted_two_numbers",
    "planted_reference_value",
    "planted_row_cdna_one",
    "planted_sibling_value",
];
/// The class table a run of `Run::new` reads.
pub const CLASS_TABLE: &str = "#class\tplugin\tfield names\tmet or not met\tcommit\n\
order_only\tNearestExonJB\tNearestExonJB\tmet\t-\n\
planted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n\
planted_two_numbers\tMaxEntScan\tMaxEntScan_*\tmet\t-\n\
planted_reference_value\tGWAS\tGWAS_*\tmet\t-\n\
planted_row_cdna_one\tUTRAnnotator\t5UTR_consequence\tmet\t-\n\
\n\
planted_sibling_value\tREVEL\tREVEL\tnot met\t-\n";
/// The registry rows that are field sets of Ensembl VEP's own and no plugin.
const NON_PLUGIN_FIELD_SETS: [&str; 1] = ["HGVS"];

static SCRATCHES: AtomicUsize = AtomicUsize::new(0);

/// A directory of a test's own, removed when the value is dropped.
pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-tests-{}-{}-{tag}",
            std::process::id(),
            SCRATCHES.fetch_add(1, Ordering::SeqCst)
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    /// The path of `name` in the directory.
    pub fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_string()
    }

    /// Writes `content` to `name`, below directories made as needed, and returns the path.
    pub fn file(&self, name: &str, content: impl AsRef<[u8]>) -> String {
        let path = self.0.join(name);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(&path, content).unwrap();
        self.path(name)
    }

    /// A default-format output: the column header and `rows`.
    pub fn output(&self, name: &str, rows: &[String]) -> String {
        self.file(name, format!("{HEADER}{}", rows.concat()))
    }

    /// A tuple file: one line per Location, Allele, Feature, field and value.
    pub fn tuples(&self, name: &str, lines: &[[&str; 5]]) -> String {
        let text: String = lines
            .iter()
            .map(|line| format!("{}\n", line.join("\t")))
            .collect();
        self.file(name, text)
    }

    /// An input VCF: the header lines and one record per chromosome, position, identifier,
    /// reference allele and alternate alleles.
    pub fn vcf(&self, name: &str, records: &[[&str; 5]]) -> String {
        let text: String = records
            .iter()
            .map(|record| format!("{}\t.\t.\t.\n", record.join("\t")))
            .collect();
        self.file(name, format!("{VCF_HEADER}{text}"))
    }

    pub fn exists(&self, name: &str) -> bool {
        self.0.join(name).exists()
    }

    /// The names in the directory `name`, sorted; none where it is no directory.
    pub fn names(&self, name: &str) -> Vec<String> {
        names(self.0.join(name))
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The names in a directory, sorted; none where it is no directory.
pub fn names(directory: impl AsRef<Path>) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(directory)
        .map(|entries| {
            entries
                .map(|entry| entry.unwrap().file_name().into_string().unwrap())
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    names
}

pub fn gzip(content: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(content).unwrap();
    encoder.finish().unwrap()
}

pub fn sha256(content: &[u8]) -> String {
    Sha256::digest(content)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

pub fn strings(words: &[&str]) -> Vec<String> {
    words.iter().map(|word| word.to_string()).collect()
}

/// A default-format row of the record `v` that writes no cDNA position.
pub fn row(location: &str, allele: &str, feature: &str, extra: &str) -> String {
    row_with(location, allele, feature, extra, "-", "v")
}

/// A default-format row with its cDNA position and its `Uploaded_variation`.
pub fn row_with(
    location: &str,
    allele: &str,
    feature: &str,
    extra: &str,
    cdna: &str,
    uploaded: &str,
) -> String {
    format!(
        "{uploaded}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t{cdna}\t-\t-\t-\t-\t-\t{extra}\n"
    )
}

/// Five rows over three Locations that write a field of five planted plugins and of the
/// custom source `Custom`.
pub fn reference_rows() -> Vec<String> {
    vec![
        row(
            "1:100",
            "T",
            "ENST1",
            "IMPACT=MODERATE;PLANT_score=0.5;PLANT_label=benign;LOEUF=0.31",
        ),
        row(
            "1:100",
            "T",
            "ENST2",
            "IMPACT=MODERATE;PLANT_score=0.5;MaxEntScan_diff=1.5",
        ),
        row("1:100", "T", "ENST3", "IMPACT=MODIFIER"),
        row(
            "1:200-201",
            "-",
            "ENST1",
            "PLANT_score=0.25;PlantedDynamic_a=x;PlantedDynamic_b=y;Custom=id1;Custom_SIG=p",
        ),
        row("10:5", "G", "-", "PLANT_label=-;GWAS_pmid=17"),
    ]
}

/// One question put to a planted registry: the field, the run's plugins and its flags.
pub type Asked = (String, Vec<String>, Option<Vec<String>>);

type Resolver =
    Box<dyn Fn(&str, &[String], Option<&[String]>) -> Result<Option<String>, String> + Send + Sync>;

/// A registry of planted rows. A field belongs to the one plugin of the run whose row lists
/// it, to none where two do, unless a resolver is planted in that lookup's place.
pub struct PlantedRegistry {
    names: Vec<String>,
    rows: Vec<(String, Vec<String>)>,
    aliases: Vec<(String, String)>,
    resolver: Option<Resolver>,
    asked: Mutex<Vec<Asked>>,
}

impl PlantedRegistry {
    /// The registry of `rows`, each a plugin and the fields it writes. A row named as a
    /// field set of Ensembl VEP's own is no plugin of the registry.
    pub fn new(rows: &[(&str, &[&str])]) -> PlantedRegistry {
        PlantedRegistry {
            names: rows
                .iter()
                .map(|(name, _)| name.to_string())
                .filter(|name| !NON_PLUGIN_FIELD_SETS.contains(&name.as_str()))
                .collect(),
            rows: rows
                .iter()
                .map(|(name, fields)| (name.to_string(), strings(fields)))
                .collect(),
            aliases: Vec::new(),
            resolver: None,
            asked: Mutex::new(Vec::new()),
        }
    }

    /// The registry of `PLANTED`.
    pub fn planted() -> PlantedRegistry {
        PlantedRegistry::new(&PLANTED)
    }

    /// The same registry with `alias` as a name a command line loads `plugin` under.
    pub fn alias(mut self, alias: &str, plugin: &str) -> PlantedRegistry {
        self.aliases.push((alias.to_string(), plugin.to_string()));
        self
    }

    /// The same registry with `resolver` in place of the lookup over the rows: it is handed
    /// the field, the run's plugins and the run's flags.
    pub fn resolving(
        mut self,
        resolver: impl Fn(&str, &[String], Option<&[String]>) -> Result<Option<String>, String>
            + Send
            + Sync
            + 'static,
    ) -> PlantedRegistry {
        self.resolver = Some(Box::new(resolver));
        self
    }

    /// Every plugin of the registry, in the rows' order.
    pub fn plugins(&self) -> Vec<String> {
        self.names.clone()
    }

    /// Every question the registry was put, in order.
    pub fn asked(&self) -> Vec<Asked> {
        self.asked.lock().unwrap().clone()
    }
}

impl Registry for PlantedRegistry {
    fn all_plugin_names(&self) -> &[String] {
        &self.names
    }

    fn plugin_alias(&self, name: &str) -> Option<&str> {
        self.aliases
            .iter()
            .find(|(alias, _)| alias == name)
            .map(|(_, plugin)| plugin.as_str())
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
        if let Some(resolver) = &self.resolver {
            return resolver(field, plugins, flags);
        }
        let writers: Vec<&String> = self
            .rows
            .iter()
            .filter(|(name, fields)| {
                fields.iter().any(|held| held == field) && plugins.contains(name)
            })
            .map(|(name, _)| name)
            .collect();
        Ok(match writers.as_slice() {
            [one] => Some(one.to_string()),
            _ => None,
        })
    }
}

/// One consult a planted class was put, its rows as copies.
#[derive(Debug, Clone)]
pub struct Call {
    pub klass: String,
    /// The reference's row.
    pub e: Entry,
    /// The field.
    pub name: String,
    pub expected: String,
    pub actual: String,
    /// The engine's row.
    pub a: Option<Entry>,
    /// The reference's rows of the row's input record.
    pub record_rows: Option<Vec<Entry>>,
    /// The engine's rows of the same input record.
    pub candidate_rows: Option<Vec<Entry>>,
    /// The applied classes the consult was handed; empty where it was handed none.
    pub declared: BTreeSet<String>,
    /// The reads the consult was handed, by their names in `Consult`, in its order.
    pub handed: Vec<&'static str>,
}

impl Call {
    fn of(consult: &Consult) -> Call {
        let copies = |rows: &[Row]| rows.iter().map(|row| Entry::clone(row)).collect();
        let handed = [
            ("gwas_pairs", consult.gwas_pairs.is_some()),
            ("revel_table", consult.revel_table.is_some()),
            ("utrannotator_sweep", consult.utrannotator_sweep.is_some()),
            (
                "downstream_sequences",
                consult.downstream_sequences.is_some(),
            ),
            ("splice_windows", consult.splice_windows.is_some()),
            ("spliceai_table", consult.spliceai_table.is_some()),
            ("exon_spans", consult.exon_spans.is_some()),
        ];
        Call {
            klass: consult.klass.to_string(),
            e: Entry::clone(consult.e),
            name: consult.name.to_string(),
            expected: consult.expected.to_string(),
            actual: consult.actual.to_string(),
            a: consult.a.map(|row| Entry::clone(row)),
            record_rows: consult.record_rows.map(copies),
            candidate_rows: consult.candidate_rows.map(copies),
            declared: consult.declared.clone(),
            handed: handed
                .iter()
                .filter(|(_, given)| *given)
                .map(|(name, _)| *name)
                .collect(),
        }
    }
}

type Accepts = Arc<dyn Fn(&Consult) -> Result<bool, Refusal> + Send + Sync>;

/// Whether the text reads as a number where CPython's `float` takes it.
pub fn is_number(value: &str) -> bool {
    pytext::float(value).is_some()
}

/// The value, or `None` where it is blank in any format's spelling.
pub fn written(value: &str) -> Option<&str> {
    (!blank(value)).then_some(value)
}

/// Planted classes: each a name, the fields it is listed with and a predicate over a
/// consult. Every consult is kept before its class's predicate is asked. A predicate is
/// asked whatever the fields its class is listed with, and a class is handed the rows and
/// reads the tool's own class of the same name is handed.
pub struct PlantedRules {
    classes: Vec<(String, Vec<String>, Accepts)>,
    calls: Mutex<Vec<Call>>,
}

impl PlantedRules {
    /// Rules that know no class.
    pub fn none() -> PlantedRules {
        PlantedRules {
            classes: Vec::new(),
            calls: Mutex::new(Vec::new()),
        }
    }

    /// The classes of `PLANTED_CLASSES`, listed with no field. `planted_engine_zero`
    /// accepts an engine's `0` beside no value of the reference; `planted_two_numbers` two
    /// values that are numbers; `planted_reference_value` any value of the reference;
    /// `planted_row_cdna_one` a reference row at cDNA position 1; `planted_sibling_value` a
    /// reference value that a row of the same record and Feature under another allele
    /// writes.
    pub fn planted() -> PlantedRules {
        PlantedRules::none()
            .class("planted_engine_zero", &[], |c| {
                written(c.expected).is_none() && written(c.actual) == Some("0")
            })
            .class("planted_two_numbers", &[], |c| {
                matches!(
                    (written(c.expected), written(c.actual)),
                    (Some(expected), Some(actual)) if is_number(expected) && is_number(actual)
                )
            })
            .class("planted_reference_value", &[], |c| {
                written(c.expected).is_some()
            })
            .class("planted_row_cdna_one", &[], |c| {
                c.e.fields.get("cDNA_position") == Some("1")
            })
            .class("planted_sibling_value", &[], |c| {
                written(c.expected)
                    .zip(c.record_rows)
                    .is_some_and(|(expected, rows)| {
                        rows.iter().any(|other| {
                            other.allele != c.e.allele
                                && other.feature == c.e.feature
                                && other.fields.get(c.name) == Some(expected)
                        })
                    })
            })
    }

    /// The same rules with the class `name`, listed with `fields`, under `accepts`; a
    /// class of that name already held is replaced in its place.
    pub fn class(
        self,
        name: &str,
        fields: &[&str],
        accepts: impl Fn(&Consult) -> bool + Send + Sync + 'static,
    ) -> PlantedRules {
        self.refusing_class(name, fields, move |consult| Ok(accepts(consult)))
    }

    /// As `class`, under a predicate that may refuse the run.
    pub fn refusing_class(
        mut self,
        name: &str,
        fields: &[&str],
        accepts: impl Fn(&Consult) -> Result<bool, Refusal> + Send + Sync + 'static,
    ) -> PlantedRules {
        let planted = (
            name.to_string(),
            strings(fields),
            Arc::new(accepts) as Accepts,
        );
        match self.classes.iter_mut().find(|(held, _, _)| held == name) {
            Some(held) => *held = planted,
            None => self.classes.push(planted),
        }
        self
    }

    /// The same rules with every class of `names`, listed with no field, under the one
    /// predicate `accepts`, which tells them apart by the consult's class.
    pub fn classes(
        mut self,
        names: &[&str],
        accepts: impl Fn(&Consult) -> bool + Send + Sync + 'static,
    ) -> PlantedRules {
        let shared: Accepts = Arc::new(move |consult| Ok(accepts(consult)));
        for name in names {
            let planted = (name.to_string(), Vec::new(), Arc::clone(&shared));
            match self.classes.iter_mut().find(|(held, _, _)| held == name) {
                Some(held) => *held = planted,
                None => self.classes.push(planted),
            }
        }
        self
    }

    /// The fields the class is listed with; `None` for a class these rules do not know.
    pub fn fields(&self, klass: &str) -> Option<Vec<String>> {
        self.classes
            .iter()
            .find(|(name, _, _)| name == klass)
            .map(|(_, fields, _)| fields.clone())
    }

    /// Every consult put to these rules, in order.
    pub fn calls(&self) -> Vec<Call> {
        self.calls.lock().unwrap().clone()
    }
}

impl Rules for PlantedRules {
    fn knows(&self, klass: &str) -> bool {
        self.classes.iter().any(|(name, _, _)| name == klass)
    }

    fn accepts(&self, consult: &Consult) -> Result<bool, Refusal> {
        self.calls.lock().unwrap().push(Call::of(consult));
        match self
            .classes
            .iter()
            .find(|(name, _, _)| name == consult.klass)
        {
            Some((_, _, accepts)) => accepts(consult),
            None => Err(Refusal(format!("{}: no planted class", consult.klass))),
        }
    }

    fn transcript_reading(&self, klass: &str) -> bool {
        BuiltIn.transcript_reading(klass)
    }

    fn sibling_allele(&self, klass: &str) -> bool {
        BuiltIn.sibling_allele(klass)
    }

    fn record_reference(&self, klass: &str) -> bool {
        BuiltIn.record_reference(klass)
    }
}

/// What a run of `score` takes beside its inputs. `Run::new` reads `CLASS_TABLE`, applies
/// every class of it, loads every plugin of the registry and writes no residual file.
#[derive(Debug, Clone)]
pub struct Run {
    /// The custom sources, each a name and its field prefix.
    pub custom: Vec<(String, String)>,
    /// The directory of the residual files.
    pub residuals: Option<String>,
    /// The class table's text; no table where `None`.
    pub table: Option<String>,
    /// A class table read as the file it is, in place of `table`.
    pub table_file: Option<String>,
    pub applied: Option<Vec<String>>,
    pub populations: Vec<(String, u64)>,
    /// The plugins the scored run loaded; every plugin of the registry where `None`.
    pub plugins: Option<Vec<String>>,
    pub catalog: Option<String>,
    pub run_input: Option<String>,
    pub sweep: Vec<String>,
    pub cache: Option<String>,
    /// The threads of the passes; two where `None`.
    pub threads: Option<i64>,
    pub revel: Option<String>,
    pub revel_no_match: bool,
    pub spliceai: Vec<String>,
    pub run_command: Option<Vec<String>>,
    pub engine_facts: Vec<(String, Vec<(String, String)>)>,
    pub met_alone: bool,
    /// The bytes of a sorted file one range of the settling step holds; the program's own
    /// figure where `None`.
    pub range_bytes: Option<u64>,
}

impl Default for Run {
    fn default() -> Run {
        Run {
            custom: Vec::new(),
            residuals: None,
            table: Some(CLASS_TABLE.to_string()),
            table_file: None,
            applied: None,
            populations: Vec::new(),
            plugins: None,
            catalog: None,
            run_input: None,
            sweep: Vec::new(),
            cache: None,
            threads: None,
            revel: None,
            revel_no_match: false,
            spliceai: Vec::new(),
            run_command: None,
            engine_facts: Vec::new(),
            met_alone: false,
            range_bytes: None,
        }
    }
}

impl Run {
    pub fn new() -> Run {
        Run::default()
    }

    pub fn custom(mut self, sources: &[(&str, &str)]) -> Run {
        self.custom = sources
            .iter()
            .map(|(name, prefix)| (name.to_string(), prefix.to_string()))
            .collect();
        self
    }

    pub fn residuals(mut self, directory: &str) -> Run {
        self.residuals = Some(directory.to_string());
        self
    }

    pub fn table(mut self, text: &str) -> Run {
        self.table = Some(text.to_string());
        self
    }

    pub fn no_table(mut self) -> Run {
        self.table = None;
        self
    }

    pub fn table_file(mut self, path: &str) -> Run {
        self.table_file = Some(path.to_string());
        self
    }

    pub fn applied(mut self, classes: &[&str]) -> Run {
        self.applied = Some(strings(classes));
        self
    }

    pub fn populations(mut self, populations: &[(&str, u64)]) -> Run {
        self.populations = populations
            .iter()
            .map(|(class, keys)| (class.to_string(), *keys))
            .collect();
        self
    }

    pub fn plugins(mut self, plugins: &[&str]) -> Run {
        self.plugins = Some(strings(plugins));
        self
    }

    pub fn catalog(mut self, path: &str) -> Run {
        self.catalog = Some(path.to_string());
        self
    }

    pub fn run_input(mut self, path: &str) -> Run {
        self.run_input = Some(path.to_string());
        self
    }

    pub fn sweep(mut self, paths: &[&str]) -> Run {
        self.sweep = strings(paths);
        self
    }

    pub fn cache(mut self, path: &str) -> Run {
        self.cache = Some(path.to_string());
        self
    }

    pub fn threads(mut self, threads: i64) -> Run {
        self.threads = Some(threads);
        self
    }

    pub fn revel(mut self, path: &str, no_match: bool) -> Run {
        self.revel = Some(path.to_string());
        self.revel_no_match = no_match;
        self
    }

    pub fn spliceai(mut self, paths: &[&str]) -> Run {
        self.spliceai = strings(paths);
        self
    }

    pub fn run_command(mut self, words: &[&str]) -> Run {
        self.run_command = Some(strings(words));
        self
    }

    /// The same run with `facts`, each a key and its value, on the report entry of the
    /// engine `ident`.
    pub fn engine_facts(mut self, ident: &str, facts: &[(&str, &str)]) -> Run {
        let facts = facts
            .iter()
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect();
        self.engine_facts.push((ident.to_string(), facts));
        self
    }

    /// The same run, which applies no class whose fourth column reads `not met` unless
    /// the run names its classes.
    pub fn met_alone(mut self) -> Run {
        self.met_alone = true;
        self
    }
}

/// Scores `reference` against each `(id, path)` of `engines` as `run` says, under
/// `registry` and `rules`. A class table given as text is written to `classes.tsv` of the
/// scratch, and the program's own scratch is made in its `work`.
pub fn try_score(
    scratch: &Scratch,
    reference: &str,
    engines: &[(&str, &str)],
    run: &Run,
    registry: &dyn Registry,
    rules: &dyn Rules,
) -> Result<Value, Stop> {
    let engines: Vec<(String, String)> = engines
        .iter()
        .map(|(ident, path)| (ident.to_string(), path.to_string()))
        .collect();
    let plugins = match &run.plugins {
        Some(plugins) => plugins.clone(),
        None => registry.all_plugin_names().to_vec(),
    };
    let options = Options {
        residuals_dir: run.residuals.clone(),
        class_table: match (&run.table_file, &run.table) {
            (Some(path), _) => Some(ClassTable::File(path.clone())),
            (None, Some(text)) => Some(ClassTable::File(scratch.file("classes.tsv", text))),
            (None, None) => None,
        },
        applied: run.applied.clone(),
        populations: run.populations.clone(),
        engine_facts: run.engine_facts.clone(),
        run_command: run.run_command.clone(),
        gwas_catalog: run.catalog.clone(),
        revel_table: run.revel.clone(),
        revel_no_match: run.revel_no_match,
        run_input: run.run_input.clone(),
        utrannotator_sweep: run.sweep.clone(),
        json_cache: run.cache.clone(),
        pairing_threads: Some(run.threads.unwrap_or(2)),
        met_alone: run.met_alone,
        spliceai_table: run.spliceai.clone(),
        settle_range_bytes: run.range_bytes,
    };
    score::score(
        reference,
        &engines,
        &plugins,
        &run.custom,
        Some(&scratch.path("work")),
        &options,
        registry,
        rules,
    )
}

/// The message of the error line a run that was stopped ends in.
pub fn message(stop: &Stop) -> String {
    match stop {
        Stop::Refused(message) => message.clone(),
        Stop::Os(failure) => failure.text(),
        Stop::Interrupted => "interrupted before the report was written".to_string(),
    }
}

/// The message a run was refused with; panics on a run that went through.
pub fn refusal(scored: Result<Value, Stop>) -> String {
    match scored {
        Ok(_) => panic!("the run was not refused"),
        Err(stop) => message(&stop),
    }
}

/// A planted registry and planted classes, as one test holds them.
pub struct Planted {
    pub registry: PlantedRegistry,
    pub rules: PlantedRules,
}

impl Planted {
    /// The registry of `PLANTED` and the classes of `PLANTED_CLASSES`.
    pub fn new() -> Planted {
        Planted {
            registry: PlantedRegistry::planted(),
            rules: PlantedRules::planted(),
        }
    }

    /// The registry of `rows` and the classes of `PLANTED_CLASSES`.
    pub fn registry(rows: &[(&str, &[&str])]) -> Planted {
        Planted {
            registry: PlantedRegistry::new(rows),
            rules: PlantedRules::planted(),
        }
    }

    pub fn try_score(
        &self,
        scratch: &Scratch,
        reference: &str,
        engines: &[(&str, &str)],
        run: &Run,
    ) -> Result<Value, Stop> {
        try_score(
            scratch,
            reference,
            engines,
            run,
            &self.registry,
            &self.rules,
        )
    }

    /// The report of a run that goes through.
    pub fn score(
        &self,
        scratch: &Scratch,
        reference: &str,
        engines: &[(&str, &str)],
        run: &Run,
    ) -> Value {
        self.try_score(scratch, reference, engines, run)
            .unwrap_or_else(|stop| panic!("the run was stopped: {}", message(&stop)))
    }

    /// Every consult put to the planted classes, in order.
    pub fn calls(&self) -> Vec<Call> {
        self.rules.calls()
    }

    /// The exit status, the lines of standard output and standard error of one command line.
    pub fn run_main(&self, arguments: &[&str]) -> (i32, Vec<String>, String) {
        run_main(arguments, &self.registry, &self.rules)
    }
}

impl Default for Planted {
    fn default() -> Planted {
        Planted::new()
    }
}

/// The exit status, the lines of standard output and standard error of one command line of
/// `score` run in process under `registry` and `rules`.
pub fn run_main(
    arguments: &[&str],
    registry: &dyn Registry,
    rules: &dyn Rules,
) -> (i32, Vec<String>, String) {
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = score::run(&strings(arguments), registry, rules, &mut out, &mut error);
    let out = String::from_utf8(out).unwrap();
    (
        code,
        pytext::splitlines(&out)
            .into_iter()
            .map(str::to_string)
            .collect(),
        String::from_utf8(error).unwrap(),
    )
}

/// Runs the built program's `score` over `arguments`, with `environment` added to the
/// caller's.
pub fn command(arguments: &[&str], environment: &[(&str, &str)]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"))
        .arg("score")
        .args(arguments)
        .envs(environment.iter().copied())
        .output()
        .unwrap()
}

/// The lines of the residual file of one source, engine and projection.
pub fn residual_lines(
    directory: &str,
    source: &str,
    engine: &str,
    projection: &str,
) -> Vec<String> {
    let path = Path::new(directory).join(format!("{source}.{engine}.{projection}.tsv"));
    let text = fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.to_string_lossy()));
    pytext::splitlines(&text)
        .into_iter()
        .map(str::to_string)
        .collect()
}

/// The residual lines of a source's feature projection against the engine `engine`.
pub fn feature_residuals(directory: &str, source: &str) -> Vec<String> {
    residual_lines(directory, source, "engine", "feature")
}

/// The report's one entry of a source, an engine and a projection.
pub fn entry<'a>(report: &'a Value, source: &str, engine: &str, projection: &str) -> &'a Value {
    let found: Vec<&Value> = report["scores"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|item| {
            item["source"] == source && item["engine"] == engine && item["projection"] == projection
        })
        .collect();
    assert_eq!(
        found.len(),
        1,
        "{source} {engine} {projection}: {}",
        report["scores"]
    );
    found[0]
}

/// The report's entry of a source's feature projection against the engine `engine`.
pub fn feature<'a>(report: &'a Value, source: &str) -> &'a Value {
    entry(report, source, "engine", "feature")
}

/// The report's entry of a source's allele projection against the engine `engine`.
pub fn allele<'a>(report: &'a Value, source: &str) -> &'a Value {
    entry(report, source, "engine", "allele")
}

/// An entry's raw and adjusted counts; panics on an entry that lacks one.
pub fn counts(item: &Value) -> [i64; 6] {
    COUNT_KEYS.map(|key| {
        item[key]
            .as_i64()
            .unwrap_or_else(|| panic!("no count {key}: {item}"))
    })
}

/// An entry's raw counts.
pub fn raw_counts(item: &Value) -> [i64; 3] {
    RAW_KEYS.map(|key| {
        item[key]
            .as_i64()
            .unwrap_or_else(|| panic!("no count {key}: {item}"))
    })
}

/// The report and the residual files, by name, of one run scored on one thread and on
/// several, once both wrote the same bytes of the report and of every file. Each run
/// writes its residual files into a directory of its own, removed after the read, and
/// leaves the program's scratch empty.
pub fn on_any_thread_count(
    scratch: &Scratch,
    reference: &str,
    engines: &[(&str, &str)],
    run: &Run,
    registry: &dyn Registry,
    rules: &dyn Rules,
) -> (Value, BTreeMap<String, Vec<u8>>) {
    let mut written = Vec::new();
    for (name, threads) in [("one", 1), ("several", 5)] {
        let residuals = scratch.path(&format!("residuals-{name}"));
        let run = run.clone().residuals(&residuals).threads(threads);
        let report = try_score(scratch, reference, engines, &run, registry, rules)
            .unwrap_or_else(|stop| panic!("the run on {name} was stopped: {}", message(&stop)));
        let files: BTreeMap<String, Vec<u8>> = names(&residuals)
            .into_iter()
            .map(|file| {
                let bytes = fs::read(Path::new(&residuals).join(&file)).unwrap();
                (file, bytes)
            })
            .collect();
        fs::remove_dir_all(&residuals).unwrap();
        assert!(scratch.0.join("work").is_dir(), "{name}");
        assert!(scratch.names("work").is_empty(), "{name}");
        written.push((score::report::dumps(&report), report, files));
    }
    let (several, one) = (written.pop().unwrap(), written.pop().unwrap());
    assert_eq!(one.0, several.0);
    assert_eq!(one.2, several.2);
    (one.1, one.2)
}
