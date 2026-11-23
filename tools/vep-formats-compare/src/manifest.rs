// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The corpus manifest and the divergences it documents.

use crate::classes;
use crate::classes::SECOND_EXPECTED_DIVERGENCE;
use crate::consts::OWN_RUN;
use crate::consts::REFERENCE_PRESENT;
use crate::consts::ROWS_ALIKE;
use crate::consts::WITH_REFERENCE_SETS;
use crate::counter::Counter;
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::Error;
use crate::error::Result;
use crate::pyfmt;
use crate::pyjson;
use crate::pyjson::DumpOptions;
use crate::pyjson::Value;
use crate::shared::consequence_set;
use crate::shared::dash;
use crate::shared::ABSENT;
use indexmap::Equivalent;
use indexmap::IndexMap;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;

/// The input a corpus manifest's own run annotates; a flag set names any other.
pub const CORPUS_INPUT: &str = "variants.vcf";

/// Undoes the percent-encoding of `=` and `;` in Extra and CSQ values.
pub fn decode_extra(s: &str) -> String {
    s.replace("%3D", "=").replace("%3B", ";")
}

/// A field value on the default format's spelling: `-` for an empty or absent
/// value, percent-encoding undone, and the `&` the VCF format writes for every comma
/// of a value (`OutputFactory/VCF.pm:401`), between the items of a list and inside a
/// name alike, read as the default format's `,`.
pub fn documented_form(s: &str) -> String {
    if s == ABSENT {
        "-".to_string()
    } else {
        dash(&decode_extra(s).replace('&', ","))
    }
}

/// One member of a tuple or a list, for `repr_tuple` and `repr_list`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Member<'a> {
    Text(&'a str),
    Number(i64),
}

fn repr_members(members: &[Member<'_>]) -> String {
    let texts: Vec<String> = members
        .iter()
        .map(|member| match member {
            Member::Text(text) => pyfmt::repr_str(text),
            Member::Number(number) => number.to_string(),
        })
        .collect();
    texts.join(", ")
}

/// `repr()` of a tuple of texts and whole numbers: `('1:100', 'A', 'ENST1', '')`, a tuple of
/// one member with its comma.
pub fn repr_tuple(members: &[Member<'_>]) -> String {
    match members {
        [_] => format!("({},)", repr_members(members)),
        _ => format!("({})", repr_members(members)),
    }
}

/// `repr()` of a list of texts and whole numbers.
pub fn repr_list(members: &[Member<'_>]) -> String {
    format!("[{}]", repr_members(members))
}

/// `repr()` of a list of texts: `['a', 'b']`.
pub fn repr_text_list<S: AsRef<str>>(texts: impl IntoIterator<Item = S>) -> String {
    let members: Vec<String> = texts
        .into_iter()
        .map(|text| pyfmt::repr_str(text.as_ref()))
        .collect();
    format!("[{}]", members.join(", "))
}

/// What a loop over a manifest value visits: the items of a list, the characters of a text,
/// the names of an object's members, and nothing for a value that reads as false.
enum Visited<'a> {
    Nothing,
    Items(&'a [Value]),
    Characters(&'a str),
    Names(&'a IndexMap<String, Value>),
}

/// A number or `true` where a list is looped over ends the run: `'int' object is not iterable`.
fn visited(value: Option<&Value>) -> Visited<'_> {
    match value {
        Some(value) if pyjson::truthy(value) => match value {
            Value::List(items) => Visited::Items(items),
            Value::Str(text) => Visited::Characters(text),
            Value::Dict(members) => Visited::Names(members),
            Value::Bool(_) => crate::uncaught!("'bool' object is not iterable"),
            Value::Float(_) => crate::uncaught!("'float' object is not iterable"),
            Value::Int(_) | Value::BigInt(_) => crate::uncaught!("'int' object is not iterable"),
            Value::Null => Visited::Nothing,
        },
        _ => Visited::Nothing,
    }
}

/// The record ordinals of a manifest item: the whole numbers of its `record_indices` that
/// are not negative. An ordinal beyond `i64` numbers no record and is left out.
pub fn _record_indices(item: &Value) -> Vec<i64> {
    match visited(item.get("record_indices")) {
        Visited::Items(items) => items
            .iter()
            .filter_map(|ordinal| match ordinal {
                Value::Int(ordinal) if *ordinal >= 0 => Some(*ordinal),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    }
}

fn text_of<'a>(item: &'a Value, key: &str) -> &'a str {
    item.get(key).and_then(Value::as_str).unwrap_or("")
}

/// The text a manifest item holds under `key`, empty where it holds none or another value.
pub fn _s(item: &Value, key: &str) -> String {
    text_of(item, key).to_string()
}

/// The consequence sets a `divergences` item states for a run that has a reference
/// sequence (`WITH_REFERENCE_SETS`), or `None` for an item without the key. The key is a
/// non-empty list of distinct sets in sorted order, each spelt as `consequence_set` spells
/// it and other than the item's `vep_rs_consequence_sets` as a whole; a list that holds
/// Ensembl VEP's set (`vep_consequence_set`) holds nothing else, and states that the two
/// engines agree on the key in such a run. A key of any other shape is refused, whatever the
/// comparison's setting.
pub fn with_reference_sets(item: &Value) -> Result<Option<Vec<String>>> {
    let Some(sets) = item.get(WITH_REFERENCE_SETS) else {
        return Ok(None);
    };
    let reason = match sets.as_list().filter(|items| !items.is_empty()) {
        None => "is not a non-empty list",
        Some(items) => {
            let texts: Option<Vec<&str>> = items
                .iter()
                .map(|set| {
                    set.as_str()
                        .filter(|set| !set.is_empty() && consequence_set(set) == *set)
                })
                .collect();
            let Some(texts) = texts else {
                return Err(malformed_reference_sets(
                    item,
                    "holds an item that is not a sorted, comma-joined set of terms",
                ));
            };
            let repeats = item
                .get("vep_rs_consequence_sets")
                .and_then(Value::as_list)
                .is_some_and(|others| {
                    others.len() == texts.len()
                        && others
                            .iter()
                            .zip(&texts)
                            .all(|(other, text)| other.as_str() == Some(text))
                });
            if !texts.windows(2).all(|pair| pair[0] < pair[1]) {
                "is not a sorted list of distinct sets"
            } else if repeats {
                "repeats vep_rs_consequence_sets and states no second outcome"
            } else if texts.len() > 1
                && texts.contains(&consequence_set(text_of(item, "vep_consequence_set")).as_str())
            {
                "names Ensembl VEP's set beside another"
            } else {
                return Ok(Some(texts.into_iter().map(str::to_string).collect()));
            }
        }
    };
    Err(malformed_reference_sets(item, reason))
}

fn malformed_reference_sets(item: &Value, reason: &str) -> Error {
    let key = ["location", "allele", "feature"].map(|name| text_of(item, name));
    Error::comparator(format!(
        "the manifest's item on {}: {WITH_REFERENCE_SETS} {reason}",
        key.join(" ")
    ))
}

fn flag_sets(manifest: &Value) -> &[Value] {
    manifest
        .get("flag_sets")
        .and_then(Value::as_list)
        .unwrap_or(&[])
}

/// Whether a flag set annotates the input named `input_name`; one without an `input`
/// annotates `CORPUS_INPUT`.
fn annotates(flag_set: &Value, input_name: &str) -> bool {
    match flag_set.get("input") {
        Some(input) => input.as_str() == Some(input_name),
        None => input_name == CORPUS_INPUT,
    }
}

/// The names of the manifest's runs over the input named `input_name`, its own run
/// (`OWN_RUN`) first where that input is `CORPUS_INPUT`; with no input named, of every run.
/// A name that is no text is given as `pyjson::py_str` writes it, a flag set without one
/// as `None`.
fn runs_over(manifest: &Value, input_name: Option<&str>) -> Vec<String> {
    let mut names = Vec::new();
    if input_name.is_none_or(|input_name| input_name == CORPUS_INPUT) {
        names.push(OWN_RUN.to_string());
    }
    names.extend(
        flag_sets(manifest)
            .iter()
            .filter(|flag_set| {
                flag_set.as_dict().is_some()
                    && input_name.is_none_or(|input_name| annotates(flag_set, input_name))
            })
            .map(|flag_set| pyjson::py_str(flag_set.get("name").unwrap_or(&Value::Null))),
    );
    names
}

/// The classes of `masks` that an item of one of the three lists of an object of `runs`
/// carries as its class, sorted. A list of another shape is not looked at.
fn documented_classes<'a>(runs: &[&Value], masks: &'a BTreeSet<String>) -> Vec<&'a str> {
    let documents = |class: &str| {
        runs.iter().any(|run| {
            ["divergences", "vep_rs_only_tuples", "field_divergences"]
                .iter()
                .filter_map(|section| run.get(section).and_then(Value::as_list))
                .flatten()
                .any(|item| {
                    item.as_dict().is_some() && text_of(item, "expected_divergence") == class
                })
        })
    };
    masks
        .iter()
        .map(String::as_str)
        .filter(|class| documents(class))
        .collect()
}

/// The object whose three lists document the run `run` of the manifest, and the run's
/// input as the manifest holds it: the manifest itself for `OWN_RUN`; for a flag set its entry
/// where the entry carries `divergences` of its own, else the manifest, as the golden harness
/// reads them. A name the manifest gives no run, and a compared input (`input_name`) other than
/// the run's, are refused: a pair is of one run.
pub fn manifest_run<'a>(
    manifest: &'a Value,
    run: &str,
    input_name: Option<&str>,
) -> Result<(&'a Value, Value)> {
    let corpus_input = || Value::Str(CORPUS_INPUT.to_string());
    let (lists, run_input) = if run == OWN_RUN {
        (manifest, corpus_input())
    } else {
        let entry = flag_sets(manifest).iter().find(|flag_set| {
            flag_set.as_dict().is_some()
                && flag_set.get("name").and_then(Value::as_str) == Some(run)
        });
        let Some(entry) = entry else {
            return Err(Error::comparator(format!(
                "--run {run}: the manifest has no flag set of that name ({OWN_RUN} names its own run)"
            )));
        };
        let lists = if entry.get("divergences").is_some() {
            entry
        } else {
            manifest
        };
        (
            lists,
            entry.get("input").cloned().unwrap_or_else(corpus_input),
        )
    };
    if let Some(input_name) = input_name {
        if run_input.as_str() != Some(input_name) {
            return Err(another_input(run, &run_input, input_name));
        }
    }
    Ok((lists, run_input))
}

fn another_input(run: &str, run_input: &Value, input_name: &str) -> Error {
    Error::comparator(format!(
        "--run {run} annotates {}, and --input names {input_name}",
        pyjson::py_str(run_input)
    ))
}

/// The name Ensembl VEP's command line loads the dbNSFP plugin under.
pub const DBNSFP_PLUGIN: &str = "dbNSFP";

/// The spellings of the flag that loads a plugin.
pub const PLUGIN_FLAGS: [&str; 2] = ["--plugin", "-plugin"];

/// The word that stands for every column of the table (dbNSFP.pm:233-235).
pub const DBNSFP_EVERY_COLUMN: &str = "ALL";

/// The opening of the word that sets the plugin's `output_prefix` (dbNSFP.pm:598-599).
pub const DBNSFP_OUTPUT_PREFIX_WORD: &str = "output_prefix=";

/// The `flags` of the runs a pair may be of: the run `run` names (`OWN_RUN` the
/// manifest's own, else the flag set of that name), and without one the manifest's own
/// where the compared input is its own or is not named, then those of each flag set over
/// the compared input.
pub fn run_flags<'a>(
    manifest: &'a Value,
    run: Option<&str>,
    input_name: Option<&str>,
) -> Vec<&'a Value> {
    let flags = |held: &'a Value| held.get("flags");
    let sets = || {
        flag_sets(manifest)
            .iter()
            .filter(|flag_set| flag_set.as_dict().is_some())
    };
    match run {
        Some(OWN_RUN) => flags(manifest).into_iter().collect(),
        Some(run) => sets()
            .find(|flag_set| flag_set.get("name").and_then(Value::as_str) == Some(run))
            .and_then(flags)
            .into_iter()
            .collect(),
        None => {
            let own = input_name.is_none_or(|input_name| input_name == CORPUS_INPUT);
            let own = flags(manifest).filter(|_| own);
            let over = sets()
                .filter(|flag_set| {
                    input_name.is_none_or(|input_name| annotates(flag_set, input_name))
                })
                .filter_map(flags);
            own.into_iter().chain(over).collect()
        }
    }
}

/// The columns the `--plugin dbNSFP,...` words of `flags` hand the plugin, each once, as
/// the words spell them. dbNSFP.pm:192-193 leave out the words that hold `=` and take the
/// first of the others as the table, and :232-241 read the rest as columns. A run that
/// names `ALL` writes every column of the table's header (:233-235), which no word states,
/// so its word gives no column. Nor does the word of a run whose last `output_prefix` is
/// neither empty nor `0`: such a run writes every column under the prefix `dbNSFP_`
/// (:257-261), and no field under the column's bare name. The word after the table names a
/// file of replacement logic where such a file exists (:220-223), which the words do not
/// tell from a column: it is kept as a column, a name no output field carries.
pub fn dbnsfp_columns(flags: &[&Value]) -> Vec<String> {
    let mut columns: Vec<String> = Vec::new();
    for flags in flags {
        let words: Vec<&str> = match flags.as_list() {
            Some(words) => words
                .iter()
                .map(|word| word.as_str().unwrap_or(""))
                .collect(),
            None => continue,
        };
        for (at, word) in words.iter().enumerate() {
            let handed = PLUGIN_FLAGS.iter().find_map(|flag| {
                if word == flag {
                    words.get(at + 1).copied()
                } else {
                    word.strip_prefix(flag)?.strip_prefix('=')
                }
            });
            let Some(handed) = handed else {
                continue;
            };
            let mut parts = handed.split(',');
            if parts.next() != Some(DBNSFP_PLUGIN) {
                continue;
            }
            // The value is the text between the first `=` and the next (:589).
            let prefixed = parts
                .clone()
                .filter_map(|part| part.strip_prefix(DBNSFP_OUTPUT_PREFIX_WORD))
                .next_back()
                .map(|value| value.split('=').next().unwrap_or(""))
                .is_some_and(|value| !value.is_empty() && value != "0");
            let named: Vec<&str> = parts
                .filter(|part| !part.is_empty() && !part.contains('='))
                .skip(1)
                .collect();
            if prefixed || named.contains(&DBNSFP_EVERY_COLUMN) {
                continue;
            }
            for column in named {
                if !columns.iter().any(|held| held == column) {
                    columns.push(column.to_string());
                }
            }
        }
    }
    columns
}

/// A key of a table of `Documented`, borrowed from an entry: the Location or the record
/// ordinal, the Allele and the Feature. It hashes as the table's own key does.
#[derive(Hash)]
struct Key<'a, P>(P, &'a str, &'a str);

impl Equivalent<(String, String, String)> for Key<'_, &str> {
    fn equivalent(&self, key: &(String, String, String)) -> bool {
        self.0 == key.0 && self.1 == key.1 && self.2 == key.2
    }
}

impl Equivalent<(i64, String, String)> for Key<'_, i64> {
    fn equivalent(&self, key: &(i64, String, String)) -> bool {
        self.0 == key.0 && self.1 == key.1 && self.2 == key.2
    }
}

/// A key of a field table of `Documented`, borrowed from an entry: `Key` and the
/// consequence set.
#[derive(Hash)]
struct FieldKey<'a, P>(P, &'a str, &'a str, &'a str);

impl Equivalent<(String, String, String, String)> for FieldKey<'_, &str> {
    fn equivalent(&self, key: &(String, String, String, String)) -> bool {
        self.0 == key.0 && self.1 == key.1 && self.2 == key.2 && self.3 == key.3
    }
}

impl Equivalent<(i64, String, String, String)> for FieldKey<'_, i64> {
    fn equivalent(&self, key: &(i64, String, String, String)) -> bool {
        self.0 == key.0 && self.1 == key.1 && self.2 == key.2 && self.3 == key.3
    }
}

/// A key of a field table as `repr()` writes the tuple.
pub trait FieldTableKey: Hash + Eq {
    fn repr(&self) -> String;
}

impl FieldTableKey for (String, String, String, String) {
    fn repr(&self) -> String {
        repr_tuple(&[
            Member::Text(&self.0),
            Member::Text(&self.1),
            Member::Text(&self.2),
            Member::Text(&self.3),
        ])
    }
}

impl FieldTableKey for (i64, String, String, String) {
    fn repr(&self) -> String {
        repr_tuple(&[
            Member::Number(self.0),
            Member::Text(&self.1),
            Member::Text(&self.2),
            Member::Text(&self.3),
        ])
    }
}

/// Where a manifest item's value is looked up: at its Location or at one of its record
/// ordinals.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Place {
    Location(String),
    Record(i64),
}

/// A key of `Documented::_documented_keys`: the place, the Allele, the Feature and the
/// field in lower case.
pub type DocumentedKey = (Place, String, String, String);

/// A `field_divergences` item that states `ROWS_ALIKE`: the item's key as a stale line
/// names it (its Location, Allele and Feature between `|`), its field, and the count of
/// rows alike it documents, 0 for a value that is no whole number from 0 up.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SharedItem {
    pub key: String,
    pub field: String,
    pub rows_alike: i64,
}

/// The fields of a key that items stating `ROWS_ALIKE` document, each with its item's place
/// in `Documented::shared_items`.
pub type SharedFields = IndexMap<String, usize>;

/// The keys a corpus manifest documents as divergent, restricted to the classes applied:
/// consequence sets per key, keys only the candidate emits, field values per key, and the
/// records the reference skipped. Keys are addressable by Location or by record ordinal.
/// `masks` is the set of classes requested, the built-in ones included.
/// `reference_sequence` is the comparison's setting where the caller states it, and
/// `run` the manifest's run whose lists alone were read where the caller names one.
/// `assembly` is the assembly the manifest names for its runs, `None` without one.
/// `dbnsfp_columns` is the columns the runs read hand the dbNSFP plugin (`run_flags`,
/// `dbnsfp_columns`), none without a manifest.
#[derive(Debug)]
pub struct Documented {
    pub masks: BTreeSet<String>,
    pub run: Option<String>,
    pub assembly: Option<String>,
    pub dbnsfp_columns: Vec<String>,
    pub reference_sequence: Option<String>,
    /// The keys an item states as equal on both engines in a run with a reference sequence,
    /// filled under `REFERENCE_PRESENT` alone: such a key is compared like any other and
    /// no class accepts a difference of its Consequence.
    pub equal_by_location: IndexMap<(String, String, String), bool>,
    pub equal_by_record: IndexMap<(i64, String, String), bool>,
    /// Per class, the manifest items counted in `applied` as they were read.
    pub read: Counter,
    /// Whether the lists read number the records of another input than the compared one: an
    /// input no run of the manifest names reads the manifest's own lists, whose ordinals
    /// count the records of `CORPUS_INPUT`. No key is then looked up by ordinal.
    pub foreign_ordinals: bool,
    pub by_location: IndexMap<(String, String, String), Vec<String>>,
    pub by_record: IndexMap<(i64, String, String), Vec<String>>,
    pub extra_by_location: BTreeSet<(String, String, String)>,
    pub extra_by_record: BTreeSet<(i64, String, String)>,
    /// A field divergence's key carries the row's consequence set when the manifest states
    /// one (empty otherwise), as the golden harness keys it: two records sharing a Location
    /// and Allele (an insertion and a substitution of the same bases) each keep their own rows.
    pub fields_by_location: IndexMap<(String, String, String, String), Fields>,
    pub fields_by_record: IndexMap<(i64, String, String, String), Fields>,
    /// The `field_divergences` items read that state `ROWS_ALIKE`, in the order they were
    /// read: several records write a row on such an item's key, and on some of those rows
    /// both engines print one value of the field.
    pub shared_items: Vec<SharedItem>,
    /// The documented fields of those items, under the same two keys as the field tables.
    pub shared_fields_by_location: IndexMap<(String, String, String, String), SharedFields>,
    pub shared_fields_by_record: IndexMap<(i64, String, String, String), SharedFields>,
    pub skipped_records: BTreeSet<i64>,
    pub available: Counter,
    pub applied: Counter,
    /// The classes a value set aside under another class needs named beside it (the second
    /// class of a value two causes change): used, though nothing is counted under them.
    pub needed: BTreeSet<String>,
    /// The classes the command line named that the comparison's reference release corrects
    /// (`classes::corrected_under`), each with the first Ensembl VEP release that corrects
    /// its defect. None of them is in `masks`.
    pub corrected: BTreeMap<String, String>,
}

impl Documented {
    pub fn new() -> Self {
        Documented {
            masks: BTreeSet::new(),
            run: None,
            assembly: None,
            dbnsfp_columns: Vec::new(),
            reference_sequence: None,
            equal_by_location: IndexMap::new(),
            equal_by_record: IndexMap::new(),
            read: Counter::new(),
            foreign_ordinals: false,
            by_location: IndexMap::new(),
            by_record: IndexMap::new(),
            extra_by_location: BTreeSet::new(),
            extra_by_record: BTreeSet::new(),
            fields_by_location: IndexMap::new(),
            fields_by_record: IndexMap::new(),
            shared_items: Vec::new(),
            shared_fields_by_location: IndexMap::new(),
            shared_fields_by_record: IndexMap::new(),
            skipped_records: BTreeSet::new(),
            available: Counter::new(),
            applied: Counter::new(),
            needed: BTreeSet::new(),
            corrected: BTreeMap::new(),
        }
    }

    /// `input_name` is the compared input's file name when the caller knows it; the runs
    /// documented over another input are then left out, since their record ordinals count
    /// that input's records. `run` names the one run of the manifest the pair belongs to
    /// (`manifest_run`): its lists alone are read, so an item of another run over the same
    /// input neither is expected on this pair's rows nor sets one aside. Without it the
    /// lists of every run over the input are read as one, which is refused where more than
    /// one run annotates the input and their lists document a class of `masks`: an item
    /// of one run would be expected on, or set aside from, the rows of another. Where they
    /// document none, every run's lists apply nothing and the pair needs no run named.
    /// `reference_sequence` states whether the
    /// candidate's run had a reference sequence: under `REFERENCE_PRESENT` an item that
    /// carries `WITH_REFERENCE_SETS` documents those sets, and one whose sets there are
    /// Ensembl VEP's own documents nothing, is not counted as applied and leaves its key to
    /// `equal_with_reference`; under `REFERENCE_ABSENT`, and where the caller states no
    /// setting, every item documents its `vep_rs_consequence_sets`.
    pub fn from_manifest(
        manifest: Option<&Value>,
        masks: &BTreeSet<String>,
        input_name: Option<&str>,
        run: Option<&str>,
        reference_sequence: Option<&str>,
    ) -> Result<Self> {
        let mut d = Documented::new();
        d.masks = masks.clone();
        d.reference_sequence = reference_sequence.map(str::to_string);
        let Some(manifest) = manifest.filter(|manifest| pyjson::truthy(manifest)) else {
            return match run {
                Some(run) => Err(Error::comparator(format!(
                    "--run {run} names a run of a manifest, and none is given"
                ))),
                None => Ok(d),
            };
        };
        d.assembly = manifest
            .get("assembly")
            .and_then(Value::as_str)
            .map(str::to_string);
        d.dbnsfp_columns = dbnsfp_columns(&run_flags(manifest, run, input_name));
        // Whether the compared input is named. A named run names its own, which the
        // manifest may hold as any value: `null` names none, and a NaN, equal to no value,
        // is refused as an input other than the run's.
        let mut input_named = input_name.is_some();
        let runs = if let Some(run) = run {
            let (lists, run_input) = manifest_run(manifest, run, input_name)?;
            input_named = match &run_input {
                Value::Null => false,
                Value::Float(number) if number.is_nan() => {
                    return Err(another_input(run, &run_input, &pyjson::py_str(&run_input)));
                }
                _ => true,
            };
            d.run = Some(run.to_string());
            vec![lists]
        } else {
            if let Some(input_name) = input_name.filter(|name| *name != CORPUS_INPUT) {
                d.foreign_ordinals = !flag_sets(manifest).iter().any(|flag_set| {
                    flag_set.as_dict().is_some() && annotates(flag_set, input_name)
                });
            }
            let runs = d.runs_read(manifest, input_name)?;
            let over = runs_over(manifest, input_name);
            let documented = documented_classes(&runs, masks);
            if over.len() > 1 && !documented.is_empty() {
                return Err(Error::comparator(format!(
                    "--run is not given, and the manifest lists {} runs{} ({}) that document \
                     {}, which --masks names: name the pair's run",
                    over.len(),
                    input_name.map_or_else(String::new, |name| format!(" over {name}")),
                    over.join(", "),
                    documented.join(", ")
                )));
            }
            runs
        };

        let present = reference_sequence == Some(REFERENCE_PRESENT);
        d.each_applied_item(
            &runs,
            "divergences",
            masks,
            input_named,
            &mut |d, _, item| {
                let sets: Vec<String> = match with_reference_sets(item)? {
                    Some(sets) if present => sets,
                    _ => match visited(item.get("vep_rs_consequence_sets")) {
                        Visited::Items(items) => items
                            .iter()
                            .filter_map(Value::as_str)
                            .map(str::to_string)
                            .collect(),
                        Visited::Characters(text) => text.chars().map(String::from).collect(),
                        Visited::Names(members) => members.keys().cloned().collect(),
                        Visited::Nothing => Vec::new(),
                    },
                };
                let (location, allele, feature) = item_key(item);
                d.by_location
                    .entry((location, allele.clone(), feature.clone()))
                    .or_default()
                    .extend(sets.iter().cloned());
                for ordinal in d.ordinals(item) {
                    d.by_record
                        .entry((ordinal, allele.clone(), feature.clone()))
                        .or_default()
                        .extend(sets.iter().cloned());
                }
                Ok(())
            },
        )?;
        d.each_applied_item(
            &runs,
            "vep_rs_only_tuples",
            masks,
            input_named,
            &mut |d, klass, item| {
                let (location, allele, feature) = item_key(item);
                d.extra_by_location
                    .insert((location, allele.clone(), feature.clone()));
                for ordinal in d.ordinals(item) {
                    d.extra_by_record
                        .insert((ordinal, allele.clone(), feature.clone()));
                    if classes::tables().in_record_dropping_classes(klass) {
                        d.skipped_records.insert(ordinal);
                    }
                }
                Ok(())
            },
        )?;
        d.each_applied_item(
            &runs,
            "field_divergences",
            masks,
            input_named,
            &mut |d, _, item| {
                let (fname, value) = (text_of(item, "field"), text_of(item, "vep_rs"));
                let (location, allele, feature) = item_key(item);
                let terms = consequence_set(text_of(item, "consequence_set"));
                let shared = item.get(ROWS_ALIKE).map(|documented| {
                    let rows_alike = match documented {
                        Value::Int(count) if *count >= 0 => *count,
                        _ => 0,
                    };
                    d.shared_items.push(SharedItem {
                        key: [location.as_str(), &allele, &feature].join("|"),
                        field: fname.to_string(),
                        rows_alike,
                    });
                    d.shared_items.len() - 1
                });
                if !location.is_empty() {
                    let key = (location, allele.clone(), feature.clone(), terms.clone());
                    if let Some(place) = shared {
                        d.shared_fields_by_location
                            .entry(key.clone())
                            .or_default()
                            .insert(fname.to_string(), place);
                    }
                    Documented::_document_field(
                        &mut d.fields_by_location,
                        key,
                        fname,
                        value,
                        true,
                    )?;
                }
                for ordinal in d.ordinals(item) {
                    let key = (ordinal, allele.clone(), feature.clone(), terms.clone());
                    if let Some(place) = shared {
                        d.shared_fields_by_record
                            .entry(key.clone())
                            .or_default()
                            .insert(fname.to_string(), place);
                    }
                    // Without the input named, the compared input may be another than the
                    // one the ordinals number; the text formats look nothing up by ordinal.
                    Documented::_document_field(
                        &mut d.fields_by_record,
                        key,
                        fname,
                        value,
                        input_named,
                    )?;
                }
                Ok(())
            },
        )?;
        Ok(d)
    }

    fn ordinals(&self, item: &Value) -> Vec<i64> {
        if self.foreign_ordinals {
            Vec::new()
        } else {
            _record_indices(item)
        }
    }

    /// Records one documented field value; with `refuse`, a key whose field the lists read
    /// document differently is refused, since a pair has one documented value of a field.
    pub fn _document_field<K: FieldTableKey>(
        table: &mut IndexMap<K, Fields>,
        key: K,
        fname: &str,
        value: &str,
        refuse: bool,
    ) -> Result<()> {
        if refuse {
            let recorded = table.get(&key).and_then(|fields| fields.get(fname));
            if let Some(recorded) = recorded.filter(|recorded| *recorded != value) {
                return Err(Error::comparator(format!(
                    "the manifest documents field {fname} on {} as both {} and {}",
                    key.repr(),
                    pyfmt::repr_str(recorded),
                    pyfmt::repr_str(value)
                )));
            }
        }
        table
            .entry(key)
            .or_default()
            .insert(fname.to_string(), value.to_string());
        Ok(())
    }

    /// The keys a manifest entry's value is looked up by: its Location and, with
    /// `by_ordinal`, each of its record ordinals, each with the Allele, the Feature and the
    /// field, lower-cased as the JSON format spells it and empty where the entry names none.
    pub fn _documented_keys(item: &Value, by_ordinal: bool) -> Vec<DocumentedKey> {
        let rest = |place: Place| {
            (
                place,
                _s(item, "allele"),
                _s(item, "feature"),
                pyfmt::lower(text_of(item, "field")),
            )
        };
        let mut keys = Vec::new();
        let location = text_of(item, "location");
        if !location.is_empty() {
            keys.push(rest(Place::Location(location.to_string())));
        }
        if by_ordinal {
            keys.extend(
                _record_indices(item)
                    .into_iter()
                    .map(|ordinal| rest(Place::Record(ordinal))),
            );
        }
        keys
    }

    /// The objects whose lists are read: the named run's (`manifest_run`), else the
    /// manifest's own run and each flag set over the compared input.
    fn runs_read<'a>(
        &self,
        manifest: &'a Value,
        input_name: Option<&str>,
    ) -> Result<Vec<&'a Value>> {
        if let Some(run) = &self.run {
            return Ok(vec![manifest_run(manifest, run, input_name)?.0]);
        }
        let selected: Vec<&Value> = flag_sets(manifest)
            .iter()
            .filter(|flag_set| {
                flag_set.as_dict().is_some()
                    && input_name.is_none_or(|input_name| annotates(flag_set, input_name))
            })
            .collect();
        let inherits = selected.is_empty()
            || selected
                .iter()
                .any(|flag_set| flag_set.get("divergences").is_none());
        let mut runs = Vec::new();
        if input_name.is_none_or(|input_name| input_name == CORPUS_INPUT) || inherits {
            runs.push(manifest);
        }
        runs.extend(selected);
        Ok(runs)
    }

    /// The (class, item) pairs of one section of the manifest's own run and of each of its
    /// `flag_sets` entries whose class `masks` names, an item documented by several runs
    /// read once; every item's class is counted in `available`, each yielded one's in
    /// `applied`. With `input_name` known, the manifest's own run counts as a run over
    /// `CORPUS_INPUT` and a flag set over another input is skipped; a selected flag set
    /// carrying no `divergences` of its own, and an input no flag set names, inherit the
    /// manifest's lists, as the golden harness reads them. An item naming a second class
    /// (`SECOND_EXPECTED_DIVERGENCE`), a value two causes change, is yielded only where
    /// `masks` names both. Where it does not, no item of the runs read that a lookup
    /// would find for the same value is yielded either: one sharing a key of
    /// `_documented_keys` with it (an ordinal names one record only where the compared
    /// input is named) under its consequence set, or with no set on either item. An item
    /// documenting the value under one class cannot be told from the one both causes change.
    /// A yielded item's second class is recorded in `needed`. Where `self.run` names a
    /// run, its lists alone are read (`manifest_run`). Under `REFERENCE_PRESENT` a
    /// `divergences` item whose `WITH_REFERENCE_SETS` is Ensembl VEP's own set is not
    /// yielded and not counted in `applied`: its key goes to `equal_by_location` and
    /// `equal_by_record`.
    pub fn _applied_items<'a>(
        &mut self,
        manifest: &'a Value,
        section: &str,
        masks: &BTreeSet<String>,
        input_name: Option<&str>,
    ) -> Result<Vec<(String, &'a Value)>> {
        let runs = self.runs_read(manifest, input_name)?;
        let mut items = Vec::new();
        self.each_applied_item(
            &runs,
            section,
            masks,
            input_name.is_some(),
            &mut |_, klass, item| {
                items.push((klass.to_string(), item));
                Ok(())
            },
        )?;
        Ok(items)
    }

    /// Calls `body` with each (class, item) pair `_applied_items` yields, one pair at a
    /// time: an item is counted, then handed to `body`, then the next is looked at.
    fn each_applied_item<'a>(
        &mut self,
        runs: &[&'a Value],
        section: &str,
        masks: &BTreeSet<String>,
        by_ordinal: bool,
        body: &mut dyn FnMut(&mut Documented, &str, &'a Value) -> Result<()>,
    ) -> Result<()> {
        let sorted_keys = DumpOptions {
            sort_keys: true,
            ..DumpOptions::default()
        };
        let mut seen: HashSet<String> = HashSet::new();
        let mut named: Vec<(&'a str, &'a Value)> = Vec::new();
        let mut withheld: HashMap<DocumentedKey, BTreeSet<String>> = HashMap::new();
        for run in runs {
            let Visited::Items(items) = visited(run.get(section)) else {
                continue;
            };
            for item in items {
                if item.as_dict().is_none() || !seen.insert(pyjson::dumps(item, &sorted_keys)) {
                    continue;
                }
                let klass = text_of(item, "expected_divergence");
                let second = text_of(item, SECOND_EXPECTED_DIVERGENCE);
                self.available.add(klass, 1);
                if !second.is_empty() && !(masks.contains(klass) && masks.contains(second)) {
                    for key in Documented::_documented_keys(item, by_ordinal) {
                        withheld
                            .entry(key)
                            .or_default()
                            .insert(consequence_set(text_of(item, "consequence_set")));
                    }
                } else if masks.contains(klass) {
                    named.push((klass, item));
                }
            }
        }
        let present = self.reference_sequence.as_deref() == Some(REFERENCE_PRESENT);
        for (klass, item) in named {
            let terms = consequence_set(text_of(item, "consequence_set"));
            let held: BTreeSet<&str> = Documented::_documented_keys(item, by_ordinal)
                .iter()
                .filter_map(|key| withheld.get(key))
                .flatten()
                .map(String::as_str)
                .collect();
            if !held.is_empty()
                && (terms.is_empty() || held.contains("") || held.contains(terms.as_str()))
            {
                continue;
            }
            if section == "divergences" && present {
                let ensembl = consequence_set(text_of(item, "vep_consequence_set"));
                if with_reference_sets(item)?.is_some_and(|sets| sets == [ensembl]) {
                    let (location, allele, feature) = item_key(item);
                    self.equal_by_location
                        .insert((location, allele.clone(), feature.clone()), true);
                    for ordinal in self.ordinals(item) {
                        self.equal_by_record
                            .insert((ordinal, allele.clone(), feature.clone()), true);
                    }
                    continue;
                }
            }
            self.applied.add(klass, 1);
            self.read.add(klass, 1);
            let second = text_of(item, SECOND_EXPECTED_DIVERGENCE);
            if !second.is_empty() {
                self.needed.insert(second.to_string());
            }
            body(self, klass, item)?;
        }
        Ok(())
    }

    /// Takes the manifest items counted as they were read out of `applied`, for a
    /// comparison none of whose entries names a Location or an input record the lists number
    /// (the VCF and JSON formats of an input that is not a VCF, or that no run of the
    /// manifest names): a manifest item is looked up by one or the other, so none meets a
    /// row there, and a row is counted where a class's own predicate takes it, once.
    pub fn leave_unmet_items_uncounted(&mut self) {
        for (klass, count) in self.read.iter() {
            self.applied.add(klass, -count);
            if self.applied.get(klass) <= 0 {
                self.applied.remove(klass);
            }
        }
        self.read.clear();
    }

    /// The value documented on the entry's key: by Location when the entry has one, then
    /// by record ordinal when it has one.
    pub fn _lookup<'a, T>(
        by_location: &'a IndexMap<(String, String, String), T>,
        by_record: &'a IndexMap<(i64, String, String), T>,
        e: &Entry,
    ) -> Option<&'a T> {
        if by_location.is_empty() && by_record.is_empty() {
            return None;
        }
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        let at_location = e
            .location
            .as_deref()
            .and_then(|location| by_location.get(&Key(location, &allele, &feature)));
        at_location.or_else(|| {
            e.record
                .and_then(|record| by_record.get(&Key(record, &allele, &feature)))
        })
    }

    /// The field values documented on the entry's key under its consequence set, else
    /// under no set (an item the manifest documents without a `consequence_set`).
    pub fn _divergent_fields(&self, e: &Entry) -> Option<&Fields> {
        if self.fields_by_location.is_empty() && self.fields_by_record.is_empty() {
            return None;
        }
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        let own = consequence_set(&e.consequence);
        for terms in [own.as_str(), ""] {
            if let Some(location) = e.location.as_deref() {
                let key = FieldKey(location, &allele, &feature, terms);
                if let Some(fields) = self.fields_by_location.get(&key) {
                    return Some(fields);
                }
            }
            if let Some(record) = e.record {
                let key = FieldKey(record, &allele, &feature, terms);
                if let Some(fields) = self.fields_by_record.get(&key) {
                    return Some(fields);
                }
            }
        }
        None
    }

    /// The documented candidate value of `fname` on the entry's key; the manifest
    /// spells the default-format name, the JSON format its lower-case form.
    pub fn divergent_field(&self, e: &Entry, fname: &str) -> Option<String> {
        let fields = self._divergent_fields(e)?;
        if let Some(value) = fields.get(fname) {
            return Some(value.clone());
        }
        let lowered = pyfmt::lower(fname);
        fields
            .iter()
            .find(|(name, _)| pyfmt::lower(name) == lowered)
            .map(|(_, value)| value.clone())
    }

    pub fn has_divergent_fields(&self, e: &Entry) -> bool {
        self._divergent_fields(e).is_some()
    }

    /// The place in `shared_items` of the item that documents `fname` on the entry's key
    /// and states `ROWS_ALIKE`, the key found as `_divergent_fields` finds it and the field
    /// named in any letter case; `None` for a field no such item documents there. On such a
    /// key a row whose two values agree is the row of another record than the documented
    /// one.
    pub fn shared_item(&self, e: &Entry, fname: &str) -> Option<usize> {
        if self.shared_items.is_empty() {
            return None;
        }
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        let own = consequence_set(&e.consequence);
        let lowered = pyfmt::lower(fname);
        let named = |fields: Option<&SharedFields>| {
            let mut fields = fields?.iter();
            let (_, place) = fields.find(|(name, _)| pyfmt::lower(name) == lowered)?;
            Some(*place)
        };
        for terms in [own.as_str(), ""] {
            if let Some(location) = e.location.as_deref() {
                let key = FieldKey(location, &allele, &feature, terms);
                if self.fields_by_location.contains_key(&key) {
                    return named(self.shared_fields_by_location.get(&key));
                }
            }
            if let Some(record) = e.record {
                let key = FieldKey(record, &allele, &feature, terms);
                if self.fields_by_record.contains_key(&key) {
                    return named(self.shared_fields_by_record.get(&key));
                }
            }
        }
        None
    }

    /// One stale line for each item of `shared_items` whose key a comparison reached and
    /// whose count of rows alike is not the number the item documents. `rows` holds, in the
    /// items' order, the rows of each item's key that carry its field and those of them on
    /// which both engines print one value: the rows of the key in a format that names a
    /// Location, the entries of the item's records in one that names records. An item whose
    /// key no row reached with the field is held to no count: the outputs compared are
    /// another run's than the item's, or lack the rows, which the comparison reports as
    /// entries of one side.
    pub fn stale_row_counts(&self, rows: &[(i64, i64)]) -> Vec<String> {
        let counted = self.shared_items.iter().zip(rows);
        counted
            .filter(|(item, (carrying, alike))| *carrying > 0 && item.rows_alike != *alike)
            .map(|(item, (_, alike))| {
                format!(
                    "{}: field {} documents {} rows alike, {} are",
                    item.key, item.field, item.rows_alike, alike
                )
            })
            .collect()
    }

    /// The documented candidate consequence sets of the entry's key, if documented.
    pub fn divergent(&self, e: &Entry) -> Option<Vec<String>> {
        Documented::_lookup(&self.by_location, &self.by_record, e).cloned()
    }

    /// Whether an item states the entry's key as equal on both engines in a run with a
    /// reference sequence, the comparison's setting: the entry is compared like any other
    /// and no class accepts a difference of its Consequence.
    pub fn equal_with_reference(&self, e: &Entry) -> bool {
        Documented::_lookup(&self.equal_by_location, &self.equal_by_record, e).is_some()
    }

    /// Whether the entry's key is documented as one only the candidate emits.
    pub fn extra(&self, e: &Entry) -> bool {
        if self.extra_by_location.is_empty() && self.extra_by_record.is_empty() {
            return false;
        }
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        if let Some(location) = &e.location {
            let key = (location.clone(), allele.clone(), feature.clone());
            if self.extra_by_location.contains(&key) {
                return true;
            }
        }
        e.record
            .is_some_and(|record| self.extra_by_record.contains(&(record, allele, feature)))
    }

    /// Whether a documented key belongs to the record, so its summary values may differ.
    pub fn record_divergent(&self, record: i64) -> bool {
        self.by_record.keys().any(|key| key.0 == record)
            || self.extra_by_record.iter().any(|key| key.0 == record)
    }

    pub fn reference_skipped(&self, record: i64) -> bool {
        self.skipped_records.contains(&record)
    }
}

/// The Location, the Allele and the Feature a manifest item names, each empty where it
/// names none.
fn item_key(item: &Value) -> (String, String, String) {
    (
        _s(item, "location"),
        _s(item, "allele"),
        _s(item, "feature"),
    )
}

impl Default for Documented {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Uncaught;

    /// A manifest of its own run over `variants.vcf` and three flag sets: one over the same
    /// input that repeats an item, one over another input, and one over a third input
    /// without lists of its own.
    const MANIFEST: &str = r#"{
 "assembly": "GRCh38",
 "divergences": [
  {"expected_divergence": "k1", "location": "1:100", "allele": "A", "feature": "ENST1", "record_indices": [0, 2], "vep_consequence_set": "intron_variant", "vep_rs_consequence_sets": ["missense_variant", "stop_gained"]},
  {"expected_divergence": "k2", "location": "1:200", "allele": "T", "feature": "ENST2", "record_indices": [1], "vep_consequence_set": "intron_variant", "vep_rs_consequence_sets": ["stop_gained"]},
  {"expected_divergence": "k1", "location": "1:300", "allele": "G", "feature": "ENST3", "record_indices": [3], "vep_consequence_set": "intron_variant", "vep_rs_consequence_sets": ["stop_gained"], "vep_rs_consequence_sets_with_reference_sequence": ["intron_variant"]}
 ],
 "vep_rs_only_tuples": [
  {"expected_divergence": "reference_skipped_record", "location": "1:400", "allele": "C", "feature": "ENST4", "record_indices": [4]},
  {"expected_divergence": "k1", "location": "1:500", "allele": "-", "feature": "-", "record_indices": [5]}
 ],
 "field_divergences": [
  {"expected_divergence": "k1", "location": "1:100", "allele": "A", "feature": "ENST1", "record_indices": [0], "consequence_set": "missense_variant", "field": "HGVSc", "vep_rs": "c.1A>G"},
  {"expected_divergence": "k1", "location": "1:100", "allele": "A", "feature": "ENST1", "record_indices": [0], "field": "IMPACT", "vep_rs": "HIGH"}
 ],
 "flag_sets": [
  {"name": "same", "divergences": [
    {"expected_divergence": "k1", "location": "1:100", "allele": "A", "feature": "ENST1", "record_indices": [0, 2], "vep_consequence_set": "intron_variant", "vep_rs_consequence_sets": ["missense_variant", "stop_gained"]},
    {"expected_divergence": "k1", "location": "1:600", "allele": "A", "feature": "ENST6", "record_indices": [6], "vep_consequence_set": "intron_variant", "vep_rs_consequence_sets": ["stop_lost"]}]},
  {"name": "other", "input": "other.vcf", "divergences": [
    {"expected_divergence": "k1", "location": "1:700", "allele": "A", "feature": "ENST7", "record_indices": [0], "vep_consequence_set": "intron_variant", "vep_rs_consequence_sets": ["start_lost"]}]},
  {"name": "inherits", "input": "third.vcf"}
 ]
}"#;

    /// Five `field_divergences` items on one key: the first names a second class.
    const TWO_CAUSES: &str = r#"{"field_divergences": [
  {"expected_divergence": "k1", "second_expected_divergence": "k2", "location": "1:100", "allele": "A", "feature": "ENST1", "record_indices": [0], "consequence_set": "missense_variant", "field": "HGVSc", "vep_rs": "both"},
  {"expected_divergence": "k1", "location": "1:100", "allele": "A", "feature": "ENST1", "consequence_set": "missense_variant", "field": "hgvsc", "vep_rs": "one"},
  {"expected_divergence": "k1", "location": "1:100", "allele": "A", "feature": "ENST1", "consequence_set": "stop_gained", "field": "HGVSc", "vep_rs": "another set"},
  {"expected_divergence": "k1", "location": "1:100", "allele": "A", "feature": "ENST1", "field": "HGVSc", "vep_rs": "no set"},
  {"expected_divergence": "k1", "location": "1:900", "allele": "A", "feature": "ENST1", "record_indices": [0], "consequence_set": "missense_variant", "field": "HGVSc", "vep_rs": "by ordinal"}]}"#;

    fn json(text: &str) -> Value {
        pyjson::loads(text).expect("the test's JSON reads")
    }

    fn masks(names: &[&str]) -> BTreeSet<String> {
        names.iter().map(|name| name.to_string()).collect()
    }

    fn documented(
        manifest: &str,
        classes: &[&str],
        input_name: Option<&str>,
        run: Option<&str>,
        reference_sequence: Option<&str>,
    ) -> Result<Documented> {
        Documented::from_manifest(
            Some(&json(manifest)),
            &masks(classes),
            input_name,
            run,
            reference_sequence,
        )
    }

    fn refusal<T: std::fmt::Debug>(outcome: Result<T>) -> String {
        outcome.expect_err("the call is refused").to_string()
    }

    fn counts(counter: &Counter) -> Vec<(&str, i64)> {
        counter.iter().collect()
    }

    fn key(place: &str, allele: &str, feature: &str) -> (String, String, String) {
        (place.to_string(), allele.to_string(), feature.to_string())
    }

    fn entry(
        location: Option<&str>,
        record: Option<i64>,
        allele: &str,
        feature: &str,
        consequence: &str,
    ) -> Entry {
        Entry {
            location: location.map(str::to_string),
            record,
            allele: allele.to_string(),
            feature: feature.to_string(),
            consequence: consequence.to_string(),
            ..Entry::default()
        }
    }

    /// The text the run ends with where `call` ends it.
    fn uncaught(call: impl FnOnce() + std::panic::UnwindSafe) -> String {
        let payload = std::panic::catch_unwind(call).expect_err("the call ends the run");
        let Uncaught(what) = payload
            .downcast_ref::<Uncaught>()
            .expect("the run ends as an uncaught exception");
        what.clone()
    }

    const ALL: &[&str] = &["k1", "k2", "reference_skipped_record"];

    #[test]
    fn percent_encoding_is_undone_for_the_equals_sign_and_the_semicolon() {
        assert_eq!(decode_extra("a%3Db&c%3Bd"), "a=b&c;d");
        assert_eq!(decode_extra("%3%3DB"), "%3=B");
    }

    #[test]
    fn a_documented_value_reads_on_the_default_formats_spelling() {
        assert_eq!(documented_form("<absent>"), "-");
        assert_eq!(documented_form(""), "-");
        assert_eq!(documented_form("a%3Db&c%3Bd"), "a=b,c;d");
    }

    #[test]
    fn a_tuple_and_a_list_print_as_the_interpreter_prints_them() {
        use Member::{Number, Text};
        assert_eq!(
            repr_tuple(&[Text("1:100"), Text("A"), Text("ENST1"), Text("")]),
            "('1:100', 'A', 'ENST1', '')"
        );
        assert_eq!(repr_tuple(&[Number(3)]), "(3,)");
        assert_eq!(repr_tuple(&[]), "()");
        assert_eq!(
            repr_tuple(&[Number(1), Text("it's"), Text("q\""), Text("bo'th\"")]),
            r#"(1, "it's", 'q"', 'bo\'th"')"#
        );
        assert_eq!(repr_list(&[Text("a"), Number(2)]), "['a', 2]");
        assert_eq!(
            repr_text_list(&["a".to_string(), "b".to_string()]),
            "['a', 'b']"
        );
        assert_eq!(repr_text_list(Vec::<String>::new()), "[]");
    }

    #[test]
    fn record_ordinals_are_the_whole_numbers_that_are_not_negative() {
        let ordinals = |item: &str| _record_indices(&json(item));
        assert_eq!(
            ordinals(r#"{"record_indices": [0, 3, -1, true, 2.0, "4", null, 7]}"#),
            [0, 3, 7]
        );
        for nothing in [
            "{}",
            r#"{"record_indices": "12"}"#,
            r#"{"record_indices": {"1": 2}}"#,
            r#"{"record_indices": 0}"#,
            r#"{"record_indices": null}"#,
        ] {
            assert_eq!(ordinals(nothing), [0i64; 0], "{nothing}");
        }
        assert_eq!(
            ordinals(r#"{"record_indices": [1, 9223372036854775808]}"#),
            [1]
        );
    }

    #[test]
    fn a_number_where_a_list_is_looped_over_ends_the_run() {
        for (item, what) in [
            (r#"{"record_indices": 5}"#, "'int' object is not iterable"),
            (
                r#"{"record_indices": 2.5}"#,
                "'float' object is not iterable",
            ),
            (
                r#"{"record_indices": true}"#,
                "'bool' object is not iterable",
            ),
        ] {
            let item = json(item);
            assert_eq!(uncaught(|| drop(_record_indices(&item))), what);
        }
    }

    #[test]
    fn a_reference_sequence_key_of_another_shape_is_refused() {
        let sets = |sets: &str, more: &str| {
            with_reference_sets(&json(&format!(
                r#"{{"location": "1:100", "allele": "A", "feature": "ENST1", "{WITH_REFERENCE_SETS}": {sets}{more}}}"#
            )))
        };
        let on = format!("the manifest's item on 1:100 A ENST1: {WITH_REFERENCE_SETS} ");
        for (value, more, reason) in [
            ("5", "", "is not a non-empty list"),
            ("[]", "", "is not a non-empty list"),
            (
                r#"["b,a"]"#,
                "",
                "holds an item that is not a sorted, comma-joined set of terms",
            ),
            (
                r#"[""]"#,
                "",
                "holds an item that is not a sorted, comma-joined set of terms",
            ),
            (r#"["x", "a"]"#, "", "is not a sorted list of distinct sets"),
            (r#"["a", "a"]"#, "", "is not a sorted list of distinct sets"),
            (
                r#"["a", "x"]"#,
                r#", "vep_rs_consequence_sets": ["a", "x"]"#,
                "repeats vep_rs_consequence_sets and states no second outcome",
            ),
            (
                r#"["a", "x"]"#,
                r#", "vep_consequence_set": "x""#,
                "names Ensembl VEP's set beside another",
            ),
        ] {
            assert_eq!(
                refusal(sets(value, more)),
                format!("{on}{reason}"),
                "{value}"
            );
        }
        assert_eq!(
            sets(r#"["x"]"#, r#", "vep_consequence_set": "x""#),
            Ok(Some(vec!["x".to_string()]))
        );
        assert_eq!(
            sets(r#"["a", "x"]"#, ""),
            Ok(Some(vec!["a".to_string(), "x".to_string()]))
        );
        assert_eq!(
            with_reference_sets(&json(r#"{"location": "1:100"}"#)),
            Ok(None)
        );
    }

    #[test]
    fn a_named_run_reads_its_own_lists_and_names_its_input() {
        let manifest = json(MANIFEST);
        let run = |run: &str, input_name: Option<&str>| {
            manifest_run(&manifest, run, input_name)
                .map(|(lists, input)| (std::ptr::eq(lists, &manifest), input))
        };
        let text = |text: &str| Value::Str(text.to_string());
        assert_eq!(run("default", None), Ok((true, text("variants.vcf"))));
        assert_eq!(
            run("same", Some("variants.vcf")),
            Ok((false, text("variants.vcf")))
        );
        assert_eq!(run("other", None), Ok((false, text("other.vcf"))));
        assert_eq!(run("inherits", None), Ok((true, text("third.vcf"))));
        assert_eq!(
            refusal(run("inherits", Some("variants.vcf"))),
            "--run inherits annotates third.vcf, and --input names variants.vcf"
        );
        assert_eq!(
            refusal(run("nope", None)),
            "--run nope: the manifest has no flag set of that name (default names its own run)"
        );
        let odd =
            json(r#"{"flag_sets": [{"name": "a", "input": [1, "x"]}, {"name": "b", "input": 5}]}"#);
        assert_eq!(
            refusal(manifest_run(&odd, "a", Some("in.vcf"))),
            "--run a annotates [1, 'x'], and --input names in.vcf"
        );
        assert_eq!(
            manifest_run(&odd, "b", None).map(|(_, input)| input),
            Ok(Value::Int(5))
        );
    }

    #[test]
    fn no_manifest_documents_nothing_and_a_named_run_needs_one() {
        let d = Documented::from_manifest(None, &masks(&["k1"]), None, None, Some("present"))
            .expect("no manifest is no error");
        assert_eq!(d.masks, masks(&["k1"]));
        assert_eq!(d.reference_sequence.as_deref(), Some("present"));
        assert!(d.by_location.is_empty() && d.applied.is_empty() && d.available.is_empty());
        let e = entry(Some("1:100"), Some(0), "A", "ENST1", "missense_variant");
        assert_eq!(d.divergent(&e), None);
        assert_eq!(d.divergent_field(&e, "HGVSc"), None);
        assert!(!d.extra(&e) && !d.equal_with_reference(&e) && !d.has_divergent_fields(&e));
        assert!(!d.record_divergent(0) && !d.reference_skipped(0));
        for manifest in [None, Some(json("{}"))] {
            let refused =
                Documented::from_manifest(manifest.as_ref(), &masks(&[]), None, Some("a"), None);
            assert_eq!(
                refusal(refused),
                "--run a names a run of a manifest, and none is given"
            );
        }
    }

    #[test]
    fn several_runs_that_document_a_named_class_are_not_read_as_one_list() {
        // No input named: every run of the manifest is read.
        assert_eq!(
            refusal(documented(MANIFEST, ALL, None, None, None)),
            "--run is not given, and the manifest lists 4 runs (default, same, other, inherits) \
             that document k1, k2, reference_skipped_record, which --masks names: name the \
             pair's run"
        );
        // The manifest's own run and the flag set `same` annotate `variants.vcf`. The class
        // `k2` is the own run's alone: read as one list it would be set aside on a pair of
        // `same`.
        assert_eq!(
            refusal(documented(
                MANIFEST,
                &["k2", "k9"],
                Some("variants.vcf"),
                None,
                None
            )),
            "--run is not given, and the manifest lists 2 runs over variants.vcf (default, same) \
             that document k2, which --masks names: name the pair's run"
        );
        // A name that is no text is listed as `py_str` writes it, a flag set without one as
        // `None`.
        assert_eq!(
            refusal(documented(
                r#"{"divergences": [{"expected_divergence": "k1"}], "flag_sets": [{"name": 7}, {}]}"#,
                &["k1"],
                Some("variants.vcf"),
                None,
                None
            )),
            "--run is not given, and the manifest lists 3 runs over variants.vcf (default, 7, None) \
             that document k1, which --masks names: name the pair's run"
        );
    }

    #[test]
    fn several_runs_are_read_unnamed_where_no_named_class_is_documented() {
        let d = documented(MANIFEST, &["k9"], None, None, None).expect("the manifest reads");
        assert_eq!(d.assembly.as_deref(), Some("GRCh38"));
        assert_eq!(d.run, None);
        // The item the flag set `same` repeats is read once.
        assert_eq!(
            counts(&d.available),
            [("k1", 7), ("k2", 1), ("reference_skipped_record", 1)]
        );
        assert!(d.applied.is_empty() && d.read.is_empty() && d.needed.is_empty());
        assert!(d.by_location.is_empty() && d.by_record.is_empty());
        assert!(d.extra_by_location.is_empty() && d.fields_by_location.is_empty());
        assert!(d.skipped_records.is_empty() && !d.foreign_ordinals);
        // No class named applies nothing: the lists of the runs over `variants.vcf` are read.
        let none = documented(MANIFEST, &[], Some("variants.vcf"), None, None).expect("reads");
        assert_eq!(
            counts(&none.available),
            [("k1", 6), ("k2", 1), ("reference_skipped_record", 1)]
        );
        assert!(none.applied.is_empty());
    }

    #[test]
    fn a_named_input_leaves_out_the_runs_over_another() {
        // The item at 1:700 is the flag set's over `other.vcf`.
        let own = documented(MANIFEST, &["k9"], Some("variants.vcf"), None, None).expect("reads");
        assert_eq!(
            counts(&own.available),
            [("k1", 6), ("k2", 1), ("reference_skipped_record", 1)]
        );

        let other = documented(MANIFEST, &["k1"], Some("other.vcf"), None, None).expect("reads");
        assert_eq!(counts(&other.available), [("k1", 1)]);
        assert_eq!(
            other.by_location.keys().cloned().collect::<Vec<_>>(),
            [key("1:700", "A", "ENST7")]
        );

        // A flag set without lists of its own inherits the manifest's.
        let third = documented(MANIFEST, &["k1"], Some("third.vcf"), None, None).expect("reads");
        assert_eq!(
            counts(&third.available),
            [("k1", 5), ("k2", 1), ("reference_skipped_record", 1)]
        );
        assert_eq!(counts(&third.applied), [("k1", 5)]);
        assert!(!third.foreign_ordinals);
        assert_eq!(third.by_record.len(), 3);

        // An input no run names reads the manifest's lists, and nothing by ordinal.
        let fourth = documented(MANIFEST, &["k1"], Some("fourth.vcf"), None, None).expect("reads");
        assert!(fourth.foreign_ordinals);
        assert_eq!(counts(&fourth.applied), [("k1", 5)]);
        assert_eq!(fourth.by_location.len(), 2);
        assert!(fourth.by_record.is_empty() && fourth.extra_by_record.is_empty());
        assert!(fourth.fields_by_record.is_empty());
        assert_eq!(fourth.fields_by_location.len(), 2);
    }

    #[test]
    fn a_named_run_documents_its_own_items_alone() {
        let same = documented(MANIFEST, &["k1"], None, Some("same"), None).expect("reads");
        assert_eq!(same.run.as_deref(), Some("same"));
        assert_eq!(counts(&same.available), [("k1", 2)]);
        assert_eq!(
            same.by_record.keys().map(|key| key.0).collect::<Vec<_>>(),
            [0, 2, 6]
        );
        assert!(same.extra_by_location.is_empty() && same.fields_by_location.is_empty());
        assert_eq!(
            refusal(documented(
                MANIFEST,
                &["k1"],
                Some("variants.vcf"),
                Some("inherits"),
                None
            )),
            "--run inherits annotates third.vcf, and --input names variants.vcf"
        );
        assert_eq!(
            refusal(documented(
                r#"{"flag_sets": [{"name": "a", "input": NaN}]}"#,
                &[],
                None,
                Some("a"),
                None
            )),
            "--run a annotates nan, and --input names nan"
        );
    }

    #[test]
    fn an_item_equal_with_a_reference_sequence_is_not_applied() {
        let d = documented(MANIFEST, ALL, None, Some("default"), Some("present")).expect("reads");
        assert_eq!(
            d.equal_by_location.keys().cloned().collect::<Vec<_>>(),
            [key("1:300", "G", "ENST3")]
        );
        assert_eq!(
            d.equal_by_record.keys().cloned().collect::<Vec<_>>(),
            [(3, "G".to_string(), "ENST3".to_string())]
        );
        assert_eq!(
            counts(&d.available),
            [("k1", 5), ("k2", 1), ("reference_skipped_record", 1)]
        );
        assert_eq!(
            counts(&d.applied),
            [("k1", 4), ("k2", 1), ("reference_skipped_record", 1)]
        );
        assert!(!d.by_location.contains_key(&key("1:300", "G", "ENST3")));
        let e = entry(Some("1:300"), None, "G", "ENST3", "intron_variant");
        assert!(d.equal_with_reference(&e));
        assert_eq!(d.divergent(&e), None);
        // Without the setting the item documents its sets like any other.
        let absent =
            documented(MANIFEST, ALL, None, Some("default"), Some("absent")).expect("reads");
        assert!(absent.equal_by_location.is_empty());
        assert_eq!(absent.divergent(&e), Some(vec!["stop_gained".to_string()]));
    }

    #[test]
    fn a_value_two_causes_change_needs_both_classes_named() {
        let field_keys = |d: &Documented| -> Vec<(String, String)> {
            d.fields_by_location
                .iter()
                .flat_map(|(key, fields)| {
                    fields
                        .values()
                        .map(|value| (format!("{} {}", key.0, key.3), value.clone()))
                })
                .collect()
        };
        let pair = |key: &str, value: &str| (key.to_string(), value.to_string());

        let both = documented(TWO_CAUSES, &["k1", "k2"], None, None, None).expect("reads");
        assert_eq!(counts(&both.applied), [("k1", 5)]);
        assert_eq!(both.needed, masks(&["k2"]));
        assert_eq!(
            field_keys(&both),
            [
                pair("1:100 missense_variant", "both"),
                pair("1:100 missense_variant", "one"),
                pair("1:100 stop_gained", "another set"),
                pair("1:100 ", "no set"),
                pair("1:900 missense_variant", "by ordinal"),
            ]
        );

        // With one class named the item is withheld, and with it the items a lookup would
        // find for the same value: the one under its set and the one under no set.
        let one = documented(TWO_CAUSES, &["k1"], None, None, None).expect("reads");
        assert_eq!(counts(&one.available), [("k1", 5)]);
        assert_eq!(counts(&one.applied), [("k1", 2)]);
        assert!(one.needed.is_empty());
        assert_eq!(
            field_keys(&one),
            [
                pair("1:100 stop_gained", "another set"),
                pair("1:900 missense_variant", "by ordinal"),
            ]
        );

        // With the input named an ordinal names one record, so the item sharing it goes too.
        let named =
            documented(TWO_CAUSES, &["k1"], Some("variants.vcf"), None, None).expect("reads");
        assert_eq!(counts(&named.applied), [("k1", 1)]);
        assert_eq!(
            field_keys(&named),
            [pair("1:100 stop_gained", "another set")]
        );
        assert!(named.fields_by_record.is_empty());
    }

    #[test]
    fn a_text_or_an_object_where_the_sets_are_listed_is_looped_over() {
        let d = documented(
            r#"{"divergences": [
                {"expected_divergence": "k1", "location": "2:1", "allele": "A", "feature": "F", "vep_rs_consequence_sets": "ab"},
                {"expected_divergence": "k1", "location": "2:2", "allele": "A", "feature": "F", "vep_rs_consequence_sets": {"k": 1, "j": 2}},
                {"expected_divergence": "k1", "location": "2:3", "allele": "A", "feature": "F", "vep_rs_consequence_sets": ["x", 1, null]},
                {"expected_divergence": "k1", "location": "2:4", "allele": "A", "feature": "F"}]}"#,
            &["k1"],
            None,
            None,
            None,
        )
        .expect("reads");
        let sets = |location: &str| d.divergent(&entry(Some(location), None, "A", "F", ""));
        let texts = |texts: &[&str]| {
            Some(
                texts
                    .iter()
                    .map(|text| text.to_string())
                    .collect::<Vec<_>>(),
            )
        };
        // The characters of a text and the member names of an object are read as sets.
        assert_eq!(sets("2:1"), texts(&["a", "b"]));
        assert_eq!(sets("2:2"), texts(&["k", "j"]));
        assert_eq!(sets("2:3"), texts(&["x"]));
        // An item that states no sets documents its key with none.
        assert_eq!(sets("2:4"), texts(&[]));
        assert_eq!(sets("2:5"), None);
        let number = json(
            r#"{"divergences": [{"expected_divergence": "k1", "vep_rs_consequence_sets": 3}]}"#,
        );
        let classes = masks(&["k1"]);
        assert_eq!(
            uncaught(|| drop(Documented::from_manifest(
                Some(&number),
                &classes,
                None,
                None,
                None
            ))),
            "'int' object is not iterable"
        );
    }

    #[test]
    fn two_values_for_one_field_are_refused() {
        let item = |location: &str, value: &str| {
            format!(
                r#"{{"expected_divergence": "k1", "location": "{location}", "allele": "A", "feature": "ENST1", "record_indices": [3], "field": "HGVSc", "vep_rs": "{value}"}}"#
            )
        };
        let clash = |second_location: &str| {
            format!(
                r#"{{"field_divergences": [{}, {}]}}"#,
                item("1:100", "it's"),
                item(second_location, "other")
            )
        };
        assert_eq!(
            refusal(documented(&clash("1:100"), &["k1"], None, None, None)),
            r#"the manifest documents field HGVSc on ('1:100', 'A', 'ENST1', '') as both "it's" and 'other'"#
        );
        assert_eq!(
            refusal(documented(
                &clash(""),
                &["k1"],
                Some("variants.vcf"),
                None,
                None
            )),
            r#"the manifest documents field HGVSc on (3, 'A', 'ENST1', '') as both "it's" and 'other'"#
        );
        // Without the input named an ordinal names a record of each input: the later value stands.
        let d = documented(&clash(""), &["k1"], None, None, None).expect("reads");
        assert_eq!(
            d.fields_by_record[&(3, "A".to_string(), "ENST1".to_string(), String::new())]["HGVSc"],
            "other"
        );
    }

    #[test]
    fn a_key_is_looked_up_by_location_then_by_record_ordinal() {
        let mut d =
            documented(MANIFEST, ALL, Some("variants.vcf"), Some("default"), None).expect("reads");
        let sets = Some(vec![
            "missense_variant".to_string(),
            "stop_gained".to_string(),
        ]);

        let at_location = entry(Some("1:100"), None, "A", "ENST1", "missense_variant");
        assert_eq!(d.divergent(&at_location), sets);
        assert!(d.has_divergent_fields(&at_location));
        assert_eq!(
            d.divergent_field(&at_location, "HGVSc").as_deref(),
            Some("c.1A>G")
        );
        assert_eq!(
            d.divergent_field(&at_location, "hgvsc").as_deref(),
            Some("c.1A>G")
        );
        assert_eq!(d.divergent_field(&at_location, "impact"), None);

        let by_ordinal = entry(None, Some(2), "A", "ENST1", "");
        assert_eq!(d.divergent(&by_ordinal), sets);
        assert!(!d.has_divergent_fields(&by_ordinal));
        assert!(d.record_divergent(2) && !d.reference_skipped(2));

        // A Location the manifest lacks falls back to the ordinal.
        let elsewhere = entry(Some("9:9"), Some(1), "T", "ENST2", "");
        assert_eq!(
            d.divergent(&elsewhere),
            Some(vec!["stop_gained".to_string()])
        );

        // An empty Allele and Feature are looked up as the default format's dash.
        let bare = entry(Some("1:500"), None, "", "", "");
        assert!(d.extra(&bare));
        assert_eq!(d.divergent(&bare), None);
        let skipped = entry(None, Some(4), "C", "ENST4", "");
        assert!(d.extra(&skipped));
        assert!(d.record_divergent(4) && d.reference_skipped(4));

        // Under another consequence set the item without a set is the one found.
        let another_set = entry(Some("1:100"), Some(0), "A", "ENST1", "stop_gained");
        assert_eq!(d.divergent_field(&another_set, "HGVSc"), None);
        assert_eq!(
            d.divergent_field(&another_set, "impact").as_deref(),
            Some("HIGH")
        );

        let unknown = entry(None, Some(7), "A", "ENST1", "");
        assert_eq!(d.divergent(&unknown), None);
        assert!(!d.extra(&unknown) && !d.record_divergent(7) && !d.reference_skipped(7));

        // The items counted as they were read leave `applied`; a row counted since stays.
        d.applied.add("k1", 2);
        d.leave_unmet_items_uncounted();
        assert_eq!(counts(&d.applied), [("k1", 2)]);
        assert!(d.read.is_empty());
    }
}
