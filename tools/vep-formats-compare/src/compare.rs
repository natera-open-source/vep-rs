// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison of two entry lists.

use crate::classes;
use crate::classes::BND_OWN_CHROMOSOME_ROWS_LOST;
use crate::classes::BREAKEND_MATE_LOCAL_READ;
use crate::classes::COLOCATED_LIST_CUT_AT_CACHE_CAP;
use crate::classes::FREQUENCY_FILTER_ZERO_VALUE;
use crate::classes::HGVS_SHIFT_SEQUENCE_EDGE_WINDOW;
use crate::classes::INDIVIDUAL_ALLELE_ORDER;
use crate::classes::MERGED_SAME_STABLE_ID_COLLAPSE;
use crate::classes::MOTIF_SCORE_REVERSE_STRAND_MNV;
use crate::classes::ORDER_ONLY;
use crate::classes::PLUGIN_HGVS_SHIFTED_CDS_COORDINATES;
use crate::classes::REFERENCE_SKIPPED_RECORD;
use crate::classes::REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS;
use crate::classes::START_COOCCURRENCE_SWAP;
use crate::classes::STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES;
use crate::classes::USE_TRANSCRIPT_REF_CACHE_CDNA;
use crate::coding_ends::terms::{CodingEndClass, LandedRule};
use crate::consts::EXAMPLES_PER_GROUP;
use crate::consts::FREQS_FIELDS;
use crate::consts::NMD_FIELDS;
use crate::consts::STRAND_FIELDS;
use crate::counter::Counter;
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::Error;
use crate::error::Result;
use crate::fields::extra_accepted;
use crate::fields::extra_class_accepted;
use crate::fields::extra_order_only_accepted;
use crate::fields::extra_pairs_accepted;
use crate::fields::field_named;
use crate::fields::populated_counts;
use crate::fields::record_label;
use crate::floors::StrListOrStr;
use crate::lists::populated_under;
use crate::lists::same_items;
use crate::manifest::documented_form;
use crate::manifest::repr_text_list;
use crate::manifest::Documented;
use crate::nearest::extra_nearest_accepted;
use crate::pyfmt;
use crate::pyjson::truthy;
use crate::pyjson::Value;
use crate::record_rows::RecordRows;
use crate::remainder::KeptBack;
use crate::rules::coding_end_pairs::{
    most_severe_terms, stands_on_the_record, AtACodingEnd, RowPair,
};
use crate::rules::cut_lists::cut_list_example;
use crate::rules::cut_lists::cut_list_reading;
use crate::rules::downstream::downstream_class_at;
use crate::rules::field_class_accepts;
use crate::rules::registry::registry_class;
use crate::rules::structural;
use crate::rules::structural::StructuralTuple;
use crate::rules::term_pairs::deleted_span_inframe_term_pair;
use crate::rules::term_pairs::freqs_all_zero;
use crate::rules::term_pairs::TERM_PAIR_CLASSES;
use crate::rules::transcript_reference::failed_edit_reference_cell;
use crate::rules::RuleContext;
use crate::rules::Side;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::feature_kind;
use crate::shared::ABSENT;
use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::LazyLock;

pub const MOTIF_SCORE_FIELDS: &[&str] = &["MOTIF_SCORE_CHANGE", "motif_score_change"];

/// A motif score change is the difference of two relative similarity scores, each in [0, 1],
/// printed to three decimals; the JSON number drops trailing zeros.
pub const MOTIF_SCORE_TEXT_PATTERN: &str = r"-?[01](?:\.[0-9]{1,3})?";
pub static MOTIF_SCORE_TEXT: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(MOTIF_SCORE_TEXT_PATTERN).expect("the pattern compiles"));

/// `MOTIF_SCORE_TEXT` held to a whole text.
static MOTIF_SCORE_WHOLE_TEXT: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(&format!(r"\A(?:{MOTIF_SCORE_TEXT_PATTERN})\z"))
        .expect("the pattern compiles")
});

/// A number with a fractional part, as a manifest spells a documented value.
static FRACTIONAL_NUMBER: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"\A-?\d+\.\d+\z").expect("the pattern compiles"));

/// The record-level token of an individual's copy: IND on a default, tab or CSQ entry
/// (OutputFactory.pm:971), ind on a JSON record object.
pub const INDIVIDUAL_FIELDS: &[&str] = &["IND", "ind"];

/// The fields carrying the copy's allele string, and the allele's index in it.
pub const ALLELE_STRING_FIELDS: &[&str] = &["allele_string", "UPLOADED_ALLELE", "uploaded_allele"];

pub const ALLELE_NUMBER_FIELDS: &[&str] = &["ALLELE_NUM", "allele_num"];

/// The fields carrying the reference allele of a --show_ref_allele row.
pub const REFERENCE_ALLELE_FIELDS: &[&str] = &["REF_ALLELE", "ref_allele"];

pub const CONSEQUENCE_FIELDS: &[&str] = &["Consequence", "consequence_terms"];

pub const SOURCE_FIELDS: &[&str] = &["SOURCE", "source"];

pub const REFSEQ_SOURCE: &str = "RefSeq";

/// The notes of a floor a format cannot count: the JSON format's row-less columns, a custom
/// source's columns in the JSON format (JSON.pm:249-258 at release/116.3 put a record's custom
/// annotations, their statistics first, on the record object alone), and a column neither
/// header of a tab or VCF output declares.
pub const FLOOR_NOT_APPLICABLE: &str = "a column the JSON format has no row count for (a record-level column is carried once per record or not at all, Feature_type names the consequence list, Extra is the text formats' container), so the floor is checked in the text formats of the same case";

pub const NEAREST_FIELDS: &[&str] = &["NEAREST", "nearest"];

/// The largest structural variant the reference annotates where a run names no other
/// bound: the bound a record's skip route is read under.
pub(crate) const MAX_SV_SIZE: i64 = 10_000_000;

/// The stale documented divergences a report lists.
const STALE_LISTED: usize = 20;

/// What `compare_entries` reads beside the two entry lists and the documented divergences.
///
/// `floors` is `None`, `field_pairs` is `None`, `check_order` is true, `examples` is
/// `EXAMPLES_PER_GROUP`, `order_only` is `None`, `floor_fields` is `None`, `input_lines` is `None`
/// and `kept` is `None` in a comparison that names none of them (`CompareOptions::new`).
#[derive(Debug)]
pub struct CompareOptions<'a> {
    pub floors: Option<&'a IndexMap<String, i64>>,
    pub field_pairs: Option<&'a IndexMap<String, BTreeSet<(String, String)>>>,
    pub check_order: bool,
    pub examples: i64,
    pub order_only: Option<&'a Fields>,
    pub floor_fields: Option<&'a IndexMap<String, Option<StrListOrStr>>>,
    pub input_lines: Option<&'a [String]>,
    pub kept: Option<&'a KeptBack>,
    pub side: &'a Side,
}

impl<'a> CompareOptions<'a> {
    /// The options of a comparison that names its side tables and nothing else: every other
    /// argument at the default of `compare_entries`.
    pub fn new(side: &'a Side) -> Self {
        CompareOptions {
            floors: None,
            field_pairs: None,
            check_order: true,
            examples: EXAMPLES_PER_GROUP,
            order_only: None,
            floor_fields: None,
            input_lines: None,
            kept: None,
            side,
        }
    }
}

/// The record an entry's key names: the text before the key's first `|`.
pub(crate) fn key_record(key: &str) -> &str {
    key.split('|').next().unwrap_or(key)
}

/// The last segment of an entry's key without its duplicate counter: the Feature_type
/// (default, tab, VCF), the JSON list's kind, or `record`.
pub(crate) fn key_kind(key: &str) -> &str {
    let last = key.rsplit('|').next().unwrap_or(key);
    last.split('#').next().unwrap_or(last)
}

/// Whether the key is a record's own entry's; a repeated record entry's key ends in its
/// counter and is none.
pub(crate) fn record_key(key: &str) -> bool {
    key.ends_with("|record")
}

/// `len(items[:n])` of a list of `len` items: a negative `n` leaves that many off the end.
fn prefix_len(len: usize, n: i64) -> usize {
    let size = usize::try_from(n.unsigned_abs()).unwrap_or(usize::MAX);
    if n >= 0 {
        len.min(size)
    } else {
        len.saturating_sub(size)
    }
}

/// `len(group) < examples`.
fn takes_example(group: &[String], examples: i64) -> bool {
    i64::try_from(group.len()).unwrap_or(i64::MAX) < examples
}

/// Whether the candidate writes the reference's own consequence set on a row of a documented
/// key whose documented `sets` hold that set beside another. Several records share such a key
/// (every `<CNV:TR>` record of one position writes one Location and Allele), the documented
/// divergence is the row of another of them, and this row is compared like any other. A key
/// documented with the reference's own set alone is set aside whole: its rows agree in their
/// terms and differ in a field.
fn agrees_on_a_shared_key(e: &Entry, a: Option<&Entry>, sets: &[String]) -> bool {
    let own = consequence_set(&e.consequence);
    let stated = |held: bool| sets.iter().any(|set| (consequence_set(set) == own) == held);
    a.is_some_and(|a| consequence_set(&a.consequence) == own) && stated(true) && stated(false)
}

fn int(value: i64) -> Value {
    Value::Int(value)
}

fn count(value: usize) -> Value {
    Value::Int(i64::try_from(value).unwrap_or(i64::MAX))
}

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

fn texts(values: &[String]) -> Value {
    Value::List(values.iter().map(|value| text(value)).collect())
}

fn dict<const N: usize>(members: [(&str, Value); N]) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect(),
    )
}

/// `f"{key}: expected {ev!r}, got {av!r}"`.
fn pair_example(key: &str, ev: &str, av: &str) -> String {
    format!(
        "{key}: expected {}, got {}",
        pyfmt::repr_str(ev),
        pyfmt::repr_str(av)
    )
}

/// Whether the entry itself carries one of `INDIVIDUAL_FIELDS` with a value.
fn carries_individual(e: &Entry) -> bool {
    INDIVIDUAL_FIELDS.iter().any(|f| {
        e.fields
            .get(*f)
            .is_some_and(|value| !value.is_empty() && value != "-")
    })
}

/// The copy an entry belongs to: its record label, name and individual.
fn copy_of(e: &Entry) -> Vec<String> {
    let field = |name: &str| e.fields.get(name).cloned().unwrap_or_default();
    let mut copy = vec![key_record(&e.key).to_string(), field("Uploaded_variation")];
    copy.extend(INDIVIDUAL_FIELDS.iter().map(|f| field(f)));
    copy
}

/// The copy and feature an entry belongs to: `copy_of` with its feature and its kind.
fn copy_feature(e: &Entry) -> Vec<String> {
    let mut copy = copy_of(e);
    copy.push(e.feature.clone());
    copy.push(key_kind(&e.key).to_string());
    copy
}

/// The entry's allele number, under either field name, or `None`.
fn allele_number(e: &Entry) -> Option<&str> {
    ALLELE_NUMBER_FIELDS
        .iter()
        .find_map(|name| e.fields.get(*name))
        .map(String::as_str)
}

/// Whether the entry is identifiably the copy's reference-allele row: numbered 0, or
/// its allele is the reference allele it prints (`-` for an insertion's); a row
/// carrying neither field is not.
fn reference_row(e: &Entry) -> bool {
    if let Some(number) = allele_number(e) {
        return number == "0";
    }
    REFERENCE_ALLELE_FIELDS
        .iter()
        .find_map(|name| e.fields.get(*name))
        .is_some_and(|reference| {
            let reference = if reference.is_empty() { "-" } else { reference };
            reference == e.allele
        })
}

/// The same reference allele followed by the same alternate alleles in any order.
fn same_alleles(expected_value: &str, actual_value: &str) -> bool {
    if blank(expected_value) || blank(actual_value) {
        return expected_value == actual_value;
    }
    let (reference, _, alternates) = pyfmt::partition(expected_value, "/");
    let (actual_reference, _, actual_alternates) = pyfmt::partition(actual_value, "/");
    reference == actual_reference && same_items(alternates, actual_alternates, "/")
}

/// Whether a text is a motif score change as either engine prints one.
fn well_formed(score: &str) -> bool {
    MOTIF_SCORE_WHOLE_TEXT.is_match(score)
        && pyfmt::float_of(score).is_some_and(|value| value.abs() <= 1.0)
}

/// Per copy, the allele numbers each of its alleles carries across the copy's rows, the
/// alleles in the order the rows first name them.
type AlleleNumbers = HashMap<Vec<String>, IndexMap<String, BTreeSet<String>>>;

/// What the readings of an individual's copy share (`INDIVIDUAL_ALLELE_ORDER`).
struct Individuals<'a> {
    allele_order: bool,
    individual_records: HashSet<&'a str>,
    expected_numbers: AlleleNumbers,
    actual_numbers: AlleleNumbers,
    order_only: &'a Fields,
}

impl Individuals<'_> {
    /// An entry of an individual's copy of a record: it carries IND itself (a default, tab
    /// or CSQ entry) or its record object carries ind (a JSON consequence entry).
    fn is_individual(&self, e: &Entry) -> bool {
        carries_individual(e) || self.individual_records.contains(key_record(&e.key))
    }

    /// Per copy, the allele numbers each of its alleles carries across the copy's rows.
    fn allele_numbers(&self, entries: &[Entry]) -> AlleleNumbers {
        let mut numbers = AlleleNumbers::new();
        for e in entries {
            if !self.is_individual(e) {
                continue;
            }
            if let Some(number) = allele_number(e) {
                numbers
                    .entry(copy_of(e))
                    .or_default()
                    .entry(e.allele.clone())
                    .or_default()
                    .insert(number.to_string());
            }
        }
        numbers
    }

    /// Whether the row's allele numbers, `expected_value` and `actual_value`, differ by
    /// one permutation of the alternates' numbers applied to every row of the copy: each
    /// allele carries one number on each side, the row's own among them, the two sides number
    /// the same alleles with the same numbers, and the allele numbered 0, the reference, is
    /// the same allele on both (Parser/VCF.pm:639 keeps the reference first).
    fn renumbered_alternates(
        &self,
        e: &Entry,
        a: &Entry,
        expected_value: &str,
        actual_value: &str,
    ) -> bool {
        let Some(expected_map) = self.expected_numbers.get(&copy_of(e)) else {
            return false;
        };
        let Some(actual_map) = self.actual_numbers.get(&copy_of(a)) else {
            return false;
        };
        if expected_map.is_empty()
            || expected_map.len() != actual_map.len()
            || !expected_map
                .keys()
                .all(|allele| actual_map.contains_key(allele))
        {
            return false;
        }
        if expected_map
            .values()
            .chain(actual_map.values())
            .any(|numbers| numbers.len() != 1)
        {
            return false;
        }
        let number_of = |map: &IndexMap<String, BTreeSet<String>>| -> IndexMap<String, String> {
            map.iter()
                .filter_map(|(allele, numbers)| Some((allele.clone(), numbers.first()?.clone())))
                .collect()
        };
        let expected_number_of = number_of(expected_map);
        let actual_number_of = number_of(actual_map);
        if expected_number_of.get(&e.allele).map(String::as_str) != Some(expected_value)
            || actual_number_of.get(&a.allele).map(String::as_str) != Some(actual_value)
        {
            return false;
        }
        let sorted_numbers = |numbers: &IndexMap<String, String>| -> Vec<String> {
            let mut values: Vec<String> = numbers.values().cloned().collect();
            values.sort();
            values
        };
        if sorted_numbers(&expected_number_of) != sorted_numbers(&actual_number_of) {
            return false;
        }
        let reference = |numbers: &IndexMap<String, String>| -> Vec<String> {
            numbers
                .iter()
                .filter(|(_, number)| *number == "0")
                .map(|(allele, _)| allele.clone())
                .collect()
        };
        reference(&expected_number_of) == reference(&actual_number_of)
    }

    /// Whether a field of an individual's entry differs in the order of the copy's
    /// alleles alone; in the Extra text every differing pair must.
    fn allele_order_accepts(&self, e: &Entry, a: &Entry, name: &str, ev: &str, av: &str) -> bool {
        if !(self.allele_order && self.is_individual(e)) {
            return false;
        }
        let mut reordered_only = |key: &str, expected_value: &str, actual_value: &str| -> bool {
            if ALLELE_STRING_FIELDS.contains(&key) {
                return same_alleles(expected_value, actual_value);
            }
            if ALLELE_NUMBER_FIELDS.contains(&key) {
                return self.renumbered_alternates(e, a, expected_value, actual_value);
            }
            self.order_only
                .get(key)
                .is_some_and(|separator| same_items(expected_value, actual_value, separator))
        };
        if name == "Extra" {
            return extra_accepted(ev, av, &mut reordered_only);
        }
        (ALLELE_STRING_FIELDS.contains(&name) || ALLELE_NUMBER_FIELDS.contains(&name))
            && reordered_only(name, ev, av)
    }

    /// The key the order comparison reads: an individual's consequence entry by its copy
    /// and feature (`copy_feature`) and whether it is the reference row, without its
    /// allele and duplicate counter, since the copy's alternate rows come in the order of
    /// its alleles while the reference row keeps its place.
    fn order_key(&self, e: &Entry) -> String {
        if !self.allele_order || record_key(&e.key) || !self.is_individual(e) {
            return e.key.clone();
        }
        let mut parts = copy_feature(e);
        let row = if reference_row(e) {
            "reference"
        } else {
            "alternate"
        };
        parts.push(row.to_string());
        parts.join("|")
    }
}

/// Whether the differing field is a MOTIF_SCORE_CHANGE the candidate writes on a
/// reverse-strand MotifFeature entry whose allele is two or more bases
/// (`MOTIF_SCORE_REVERSE_STRAND_MNV`); a score the candidate leaves out where the
/// reference writes one is a difference, and in the Extra text that pair alone may
/// differ.
fn reverse_motif_score_accepts(
    reverse_motif: bool,
    e: &Entry,
    name: &str,
    ev: &str,
    av: &str,
) -> bool {
    // The Allele column alone does not tell an insertion's bases from a substitution's, so
    // an insertion's row reads as a multi-base allele here.
    let bases = e.allele.len() >= 2 && e.allele.bytes().all(|b| b"ACGT".contains(&b));
    if !(reverse_motif && bases) {
        return false;
    }
    if feature_kind(key_kind(&e.key)) != "motif_feature" {
        return false;
    }
    if !STRAND_FIELDS
        .iter()
        .any(|f| e.fields.get(*f).is_some_and(|strand| strand == "-1"))
    {
        return false;
    }
    if name == "Extra" {
        let mut score_alone =
            |key: &str, _: &str, actual: &str| key == "MOTIF_SCORE_CHANGE" && well_formed(actual);
        if extra_accepted(ev, av, &mut score_alone) {
            return true;
        }
        // Ensembl VEP writes no score where the allele runs past the motif's end; the
        // candidate writes one.
        let scored = |pair: &str| pair.starts_with("MOTIF_SCORE_CHANGE=");
        let actual_pairs: Vec<&str> = av.split(';').collect();
        let without_score: Vec<&str> = actual_pairs
            .iter()
            .copied()
            .filter(|pair| !scored(pair))
            .collect();
        let scores: Vec<&str> = actual_pairs
            .iter()
            .copied()
            .filter(|pair| scored(pair))
            .map(|pair| pyfmt::partition(pair, "=").2)
            .collect();
        let expected_scored = ev.contains("MOTIF_SCORE_CHANGE=");
        return scores.len() == 1
            && well_formed(scores[0])
            && !expected_scored
            && without_score.join(";") == ev;
    }
    MOTIF_SCORE_FIELDS.contains(&name) && well_formed(av)
}

/// What the built-in field classes read of one comparison.
struct FieldClasses<'a> {
    field_classes: Vec<&'static str>,
    record_index: Option<RecordRows<'a>>,
    candidate_index: Option<RecordRows<'a>>,
    /// A VCF or JSON consequence entry's span is on its record's entry (`indel_span`).
    record_entries: HashMap<i64, &'a Entry>,
    term_pair_rows: HashMap<&'a str, &'static str>,
    side: &'a Side,
    assembly: Option<String>,
    dbnsfp_columns: Vec<String>,
    order_only: &'a Fields,
}

impl FieldClasses<'_> {
    /// The built-in field classes that accept the differing field, a field two classes
    /// accept being the first's in `FIELD_CLASSES` order: in the `Extra` text, one value
    /// carrying every pair, each differing pair must be accepted by one of them, agree
    /// under `order_only` or, a `NEAREST` pair, agree by its start sites
    /// (`nearest_sites`), and each class used is credited. The class an accepted
    /// Downstream value needs named beside its own (`downstream_class`) is added to
    /// `needed`.
    fn field_class_accepting(
        &mut self,
        e: &Entry,
        a: &Entry,
        name: &str,
        ev: &str,
        av: &str,
        needed: &mut BTreeSet<String>,
    ) -> Result<Vec<&'static str>> {
        if self.field_classes.is_empty() {
            return Ok(Vec::new());
        }
        let record_rows = match self.record_index.as_mut() {
            Some(index) => index.of(e)?,
            None => None,
        };
        let candidate_rows = match self.candidate_index.as_mut() {
            Some(index) => index.of(a)?,
            None => None,
        };
        let tables = classes::tables();
        // A row a class of `TERM_PAIR_FOLLOWING_CLASSES` took reads its cells under that
        // class before every other, and such a class reads the cells of no other row but
        // `STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES`, which reads after the row's own class
        // the cells of a row an earlier such class took where its rule takes the row too.
        // That class also reads, among the field classes and beside the side tables, a row
        // whose two sets of terms are equal, which its own predicate does not take:
        // `unpaired` names the class on such a row.
        let following = self
            .term_pair_rows
            .get(e.key.as_str())
            .copied()
            .filter(|klass| tables.in_term_pair_following_classes(klass));
        let settled_beside = following.is_some_and(|klass| {
            klass != STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES
                && self
                    .field_classes
                    .contains(&STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES)
        });
        let equal_terms = self
            .field_classes
            .contains(&STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES)
            && consequence_set(&e.consequence) == consequence_set(&a.consequence);
        let unpaired = equal_terms.then_some(STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES);
        let alone = RuleContext {
            record_rows: None,
            candidate_rows: None,
            declared: &[],
            assembly: None,
            dbnsfp_columns: &[],
            side: self.side,
        };
        let context = RuleContext {
            record_rows: record_rows.as_deref(),
            candidate_rows: candidate_rows.as_deref(),
            declared: &self.field_classes,
            assembly: self.assembly.as_deref(),
            dbnsfp_columns: &self.dbnsfp_columns,
            side: self.side,
        };
        let mut found: BTreeSet<String> = BTreeSet::new();
        let mut first_accepting =
            |key: &str, expected_value: &str, actual_value: &str| -> Result<Option<&'static str>> {
                if let Some(following) = following {
                    let accepted = field_class_accepts(
                        following,
                        e,
                        key,
                        expected_value,
                        actual_value,
                        Some(a),
                        &alone,
                    )?;
                    if accepted {
                        return Ok(Some(following));
                    }
                    if settled_beside
                        && field_class_accepts(
                            STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES,
                            e,
                            key,
                            expected_value,
                            actual_value,
                            Some(a),
                            &alone,
                        )?
                    {
                        return Ok(Some(STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES));
                    }
                }
                for &klass in &self.field_classes {
                    if tables.in_term_pair_following_classes(klass) && Some(klass) != unpaired {
                        continue;
                    }
                    if klass == HGVS_SHIFT_SEQUENCE_EDGE_WINDOW {
                        let record = e
                            .record
                            .and_then(|ordinal| self.record_entries.get(&ordinal).copied());
                        let accepted = self
                            .side
                            .sequence_edges
                            .as_ref()
                            .is_some_and(|edges| edges.accepts(e, Some(a), record, key));
                        if accepted {
                            return Ok(Some(klass));
                        }
                        continue;
                    }
                    let accepted = field_class_accepts(
                        klass,
                        e,
                        key,
                        expected_value,
                        actual_value,
                        Some(a),
                        &context,
                    )?;
                    if !accepted {
                        continue;
                    }
                    if klass == PLUGIN_HGVS_SHIFTED_CDS_COORDINATES && !field_named(key, NMD_FIELDS)
                    {
                        // The class accepts no such field without the candidate's transcripts.
                        if let Some(sequences) = self.side.downstream_sequences.as_ref() {
                            let decided = self.side.asked_at(e, |at| {
                                downstream_class_at(
                                    e,
                                    a,
                                    key,
                                    expected_value,
                                    actual_value,
                                    sequences,
                                    record_rows.as_deref(),
                                    at,
                                )
                            })?;
                            if let Some((_, Some(second))) = decided {
                                found.insert(second);
                            }
                        }
                    }
                    return Ok(Some(klass));
                }
                Ok(None)
            };
        let accepted: Vec<&'static str> = if name == "Extra" {
            let mut used: Vec<&'static str> = Vec::new();
            let mut failed: Option<Error> = None;
            let mut accept_pair = |key: &str, expected_value: &str, actual_value: &str| -> bool {
                if failed.is_some() {
                    return false;
                }
                match first_accepting(key, expected_value, actual_value) {
                    Ok(Some(klass)) => {
                        if !used.contains(&klass) {
                            used.push(klass);
                        }
                        true
                    }
                    Ok(None) => {
                        let same_start_sites = key == "NEAREST"
                            && self
                                .side
                                .nearest_sites
                                .as_ref()
                                .is_some_and(|sites| sites.agree(expected_value, actual_value));
                        same_start_sites
                            || self.order_only.get(key).is_some_and(|separator| {
                                same_items(expected_value, actual_value, separator)
                            })
                    }
                    Err(error) => {
                        failed = Some(error);
                        false
                    }
                }
            };
            let every_pair = extra_class_accepted(ev, av, &mut accept_pair);
            if let Some(error) = failed {
                return Err(error);
            }
            if every_pair {
                used
            } else {
                Vec::new()
            }
        } else {
            first_accepting(name, ev, av)?.into_iter().collect()
        };
        if !accepted.is_empty() {
            needed.extend(found);
        }
        Ok(accepted)
    }
}

/// The rows the candidate alone writes that a structural-variant class sets aside, with the
/// records carrying one.
#[derive(Default)]
struct CandidateOnly<'a> {
    rows: HashSet<&'a str>,
    records: HashSet<i64>,
}

/// The entries of one side grouped by consequence text in its byte order, each with its
/// count and its first keys.
fn grouped(entries: &[&Entry], examples: i64) -> Value {
    let mut by_consequence: IndexMap<&str, Vec<&Entry>> = IndexMap::new();
    for e in entries {
        by_consequence
            .entry(e.consequence.as_str())
            .or_default()
            .push(e);
    }
    by_consequence.sort_keys();
    Value::Dict(
        by_consequence
            .into_iter()
            .map(|(consequence, items)| {
                let shown = prefix_len(items.len(), examples);
                let keys = items[..shown].iter().map(|e| text(&e.key)).collect();
                (
                    consequence.to_string(),
                    dict([
                        ("count", count(items.len())),
                        ("examples", Value::List(keys)),
                    ]),
                )
            })
            .collect(),
    )
}

/// Compares expected against actual entries and returns the format's report: the
/// one-sided key counts, every field with its differing count and examples grouped by the
/// expected consequence set, the stale documented divergences, the order of the shared
/// keys, the floors, and the `verdict`. `order_only` names the fields the class compares
/// as lists, with their list separators (`order_only_fields`), each agreement counted in
/// `documented.applied` under the class; `floor_fields` maps a floor's column to the
/// fields `populated_under` counts it by in a format whose keys differ from the column
/// names, to `None` for a column the JSON format has no row count for, or to the note of
/// a floor the format cannot count for another reason; a column absent from the map, and
/// every column of a format without one, is counted by its own name over every entry.
/// `codon_tables` is the candidate cache's codon tables, which
/// `organelle_peptide_paths_translation_table` reads.
/// `downstream_sequences` is the candidate cache's transcripts, which the classes of
/// `TRANSCRIPT_READING_CLASSES` read.
/// `transcript_edges` is the candidate cache's exon spans as
/// `hgvs_clipped_reference_against_whole_alternate`,
/// `hgvs_insertion_left_where_the_shift_stops` and
/// `hgvsp_shifted_indel_written_span_gate` count them (`TranscriptEdges`).
/// `INDIVIDUAL_ALLELE_ORDER` in
/// `documented.masks` compares an individual's alleles as a set: the alternate rows of one
/// copy may come in any order within a feature, its allele strings agree on the reference
/// allele and the set of alternate alleles, and its `ALLELE_NUM` values differ by one
/// permutation of the alternates' numbers applied to every row of the copy; each agreement
/// is counted under the class, the reordered entries under `counts`.
/// `sequence_edges` is what `HGVS_SHIFT_SEQUENCE_EDGE_WINDOW` reads as a field class, the
/// reference sequence and the candidate cache's transcripts; without it the class applies
/// from a manifest alone.
/// `MOTIF_SCORE_REVERSE_STRAND_MNV` in `documented.masks` accepts a differing
/// MOTIF_SCORE_CHANGE on a reverse-strand motif entry whose allele is two or more bases,
/// counted under the class. `MERGED_SAME_STABLE_ID_COLLAPSE` in `documented.masks`
/// accepts a Transcript entry the candidate alone writes, of a model the reference wrote no row
/// for on a record it wrote a Transcript row for, its SOURCE RefSeq or absent, once per record
/// and allele, counted under the class, and leaves such a record's `most_severe_consequence`
/// uncompared. Each class of `TERM_PAIR_CLASSES` in `documented.masks` accepts the
/// Consequence of a pair its predicate admits, the row counted once and its other fields
/// compared, and leaves such a record's `most_severe_consequence` uncompared; a predicate is
/// asked with the coding ends `side.partial_codon_models` holds of the row's transcript. On
/// a record that carries a row a class of a coding end took, and no row another mechanism
/// sets aside, the `most_severe_consequence` is each side's most severe term over its own
/// entries of the record, and a difference where it is not. On a row of a
/// class of `TERM_PAIR_FOLLOWING_CLASSES` the cells `term_pair_following_cell` accepts, or
/// for a class of a coding end `coding_end_following_cell`, are
/// counted under that class as well. Where the rule of
/// `STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES` takes a row an earlier such class holds, the
/// residues, the protein notation and the HGVS_OFFSET it reads are counted under it beside
/// that class. On a row whose two sets of terms are equal that class reads, beside
/// `downstream_sequences`, the protein notation and the HGVS_OFFSET of an insertion written
/// in an intron and carried into the coding sequence
/// (`rules::settled_stop::carried_insertion_cell`). `REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS` in
/// `documented.masks` leaves uncompared the `most_severe_consequence` of a record carrying
/// a row whose Consequence is the class's term pair, and `USE_TRANSCRIPT_REF_CACHE_CDNA`
/// that of a record carrying a row whose Consequence it takes. The classes
/// of `structural_variants.classes` (`STRUCTURAL_VARIANT_CLASSES`) are decided by the
/// structural-variant comparator's predicates over the entries' tuples: a
/// `breakend_mate_local_read` pair is skipped whole as a registry pair is, a
/// `reference_skipped_record` row and a `bnd_own_chromosome_rows_lost` row are
/// candidate-only rows counted under their class, a skipped record's record entry is counted
/// as a manifest's skipped record's is, and the `most_severe_consequence` of a record
/// carrying a pair or a lost own-end row is left uncompared. `intron_ends` is the candidate
/// cache's intron ends, without which `hgvsp_intron_end_insertion` and
/// `hgvsp_intron_start_insertion` accept nothing. `splice_windows` is the candidate
/// cache's introns, without which `maxentscan_overlapping_windows` accepts nothing. A key
/// the manifest states as equal on
/// both engines in a run with a reference sequence (`Documented.equal_with_reference`) is
/// left to no built-in row class, and a difference of its Consequence is accepted by no
/// class, pair or rule. A row of a documented key whose documented sets hold the reference's
/// own set beside another, on which the candidate writes that set, is compared like any other
/// (`agrees_on_a_shared_key`). Where no entry of either side names a Location or an input record
/// the lists number, the manifest items are taken out of the applied counts
/// (`Documented.leave_unmet_items_uncounted`).
/// `kept` is what a front end kept back (`hand_over`), the entries then being the
/// remainder's: its entries are added to the entry counts, the compared pairs and the
/// populated columns, the place where the order of the shared keys differs counts the
/// kept-back entries before it, a JSON floor counts the kept-back entries under its fields
/// (`KeptBack::counted_under`), and the kept-back entries are among those asked whether
/// one names a Location or an input record.
pub fn compare_entries<'a>(
    expected: &'a [Entry],
    actual: &'a [Entry],
    documented: &mut Documented,
    options: &CompareOptions,
) -> Result<Value> {
    let no_floors = IndexMap::new();
    let floors = options.floors.unwrap_or(&no_floors);
    let no_field_pairs = IndexMap::new();
    let field_pairs = options.field_pairs.unwrap_or(&no_field_pairs);
    let no_order_only = Fields::new();
    let order_only = options.order_only.unwrap_or(&no_order_only);
    let no_floor_fields = IndexMap::new();
    let floor_fields = options.floor_fields.unwrap_or(&no_floor_fields);
    let examples = options.examples;
    let kept = options.kept;
    let side = options.side;
    let nearest_sites = side.nearest_sites.as_ref();
    let spliceai_table = side.spliceai_table.as_ref();
    if let Some(table) = spliceai_table {
        table.refusals.clear();
    }
    let actual_by_key: HashMap<&str, &Entry> = actual.iter().map(|e| (e.key.as_str(), e)).collect();
    let expected_keys: HashSet<&str> = expected.iter().map(|e| e.key.as_str()).collect();
    let unmet = match kept {
        // A kept-back entry of a VCF or JSON comparison names an input record and no
        // Location.
        Some(kept) if kept.entries != 0 && kept.line_runs.is_some() => documented.foreign_ordinals,
        _ => {
            (!expected.is_empty() || !actual.is_empty())
                && expected.iter().chain(actual).all(|e| {
                    e.location.is_none() && (e.record.is_none() || documented.foreign_ordinals)
                })
        }
    };
    if unmet {
        documented.leave_unmet_items_uncounted();
    }

    // Whether the manifest speaks for the entry's key, so that no built-in row class
    // reads the row: a documented key, or one stated as equal with a reference sequence.
    let manifest_reads = |documented: &Documented, e: &Entry| -> bool {
        documented.divergent(e).is_some() || documented.equal_with_reference(e)
    };

    let mut missing: Vec<&Entry> = Vec::new();
    let mut unexpected: Vec<&Entry> = Vec::new();
    let mut mismatches: IndexMap<String, IndexMap<String, Vec<String>>> = IndexMap::new();
    let mut mismatch_totals = Counter::new();
    let mut pairs_accepted = Counter::new();
    let mut nearest_accepted = Counter::new();
    let mut order_only_accepted = Counter::new();
    let mut order_only_examples: IndexMap<String, Vec<String>> = IndexMap::new();
    let allele_order = documented.masks.contains(INDIVIDUAL_ALLELE_ORDER);
    let mut allele_order_accepted = Counter::new();
    let mut allele_order_examples: IndexMap<String, Vec<String>> = IndexMap::new();
    let reverse_motif = documented.masks.contains(MOTIF_SCORE_REVERSE_STRAND_MNV);
    let mut reverse_motif_accepted = Counter::new();
    let mut reverse_motif_examples: IndexMap<String, Vec<String>> = IndexMap::new();
    // The built-in classes of the corrected behaviours: the registry pairs are found before
    // the entries are compared, since a record entry's summary may precede its pair; a field
    // class is credited once per entry it accepts a field of, whatever the format splits the
    // row into, so its applied count is the same in every format (a class that reads a record's
    // alleles leaves a default or tab row under an identifier, `RecordRows.one_record`).
    let tables = classes::tables();
    let registry: Vec<&str> = tables
        .registry_classes()
        .iter()
        .copied()
        .filter(|k| documented.masks.contains(*k))
        .collect();
    let field_classes: Vec<&'static str> = tables
        .field_classes()
        .iter()
        .map(|(k, _)| *k)
        .filter(|k| documented.masks.contains(*k))
        .collect();
    let record_index = field_classes
        .iter()
        .any(|k| tables.in_record_reading_classes(k))
        .then(|| RecordRows::new(expected, options.input_lines));
    let candidate_index = field_classes
        .iter()
        .any(|k| tables.in_sibling_allele_classes(k))
        .then(|| RecordRows::new(actual, options.input_lines));
    let mut record_entries: HashMap<i64, &Entry> = HashMap::new();
    if side.sequence_edges.is_some() {
        for e in expected {
            if let Some(record) = e.record {
                if record_key(pyfmt::partition(&e.key, "#").0) {
                    record_entries.insert(record, e);
                }
            }
        }
    }
    let mut registry_pairs: HashMap<&str, String> = HashMap::new();
    let mut registry_records: HashSet<i64> = HashSet::new();
    if !registry.is_empty() {
        for e in expected {
            let Some(a) = actual_by_key.get(e.key.as_str()) else {
                continue;
            };
            if record_key(&e.key) || manifest_reads(documented, e) {
                continue;
            }
            let klass = registry_class(
                &consequence_set(&e.consequence),
                &consequence_set(&a.consequence),
                &e.allele,
                &registry,
            );
            // Beside the candidate's cache the pair of both start terms is the class's only
            // on a transcript that shows a start codon.
            let rule = LandedRule::StartCooccurrenceSwap;
            let klass = match klass {
                Some(klass)
                    if klass == START_COOCCURRENCE_SWAP
                        && !side.asked_at(e, |at| Ok(stands_on_the_record(rule, at)))? =>
                {
                    None
                }
                klass => klass,
            };
            if let Some(klass) = klass {
                registry_pairs.insert(&e.key, klass);
                if let Some(record) = e.record {
                    registry_records.insert(record);
                }
            }
        }
    }
    let mut class_examples: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut class_entries: IndexMap<String, HashSet<&str>> = IndexMap::new();
    let mut class_columns: IndexMap<String, Counter> = IndexMap::new();
    // The built-in row class `merged_same_stable_id_collapse`: a Transcript entry the candidate
    // alone writes, its SOURCE RefSeq or absent, on a record the reference wrote a Transcript
    // row for, of a model the reference wrote no row for on that record (a record it skipped
    // whole, or wrote an intergenic row for because every model it overlaps is a lost copy, is
    // another shape; so are another allele's row and a duplicate of a row of a model it kept,
    // since the collapse drops a model's rows and never adds one). Found before the entries
    // are compared: only the pass over the actual entries reaches a candidate-only row, and the
    // summary of a record carrying one is left uncompared where the record entries are compared,
    // before that pass.
    let mut collapsed_rows: HashSet<&str> = HashSet::new();
    let mut collapsed_records: HashSet<i64> = HashSet::new();
    if documented.masks.contains(MERGED_SAME_STABLE_ID_COLLAPSE) {
        let transcript_entry = |e: &Entry| feature_kind(key_kind(&e.key)) == "transcript";
        let kept_models: HashSet<(String, &str)> = expected
            .iter()
            .filter(|e| transcript_entry(e))
            .map(|e| (record_label(e), e.feature.as_str()))
            .collect();
        let transcript_labels: HashSet<&str> = kept_models
            .iter()
            .map(|(label, _)| label.as_str())
            .collect();
        // A lost copy's row once per record and allele: a repeat of one is a duplicate, whatever
        // the format's key spells (the text key omits the Uploaded_variation the label carries).
        let mut accepted_rows: HashSet<(String, &str, &str)> = HashSet::new();
        for a in actual {
            if expected_keys.contains(a.key.as_str()) || documented.extra(a) || !transcript_entry(a)
            {
                continue;
            }
            let label = record_label(a);
            if !transcript_labels.contains(label.as_str())
                || kept_models.contains(&(label.clone(), a.feature.as_str()))
            {
                continue;
            }
            let source = SOURCE_FIELDS
                .iter()
                .find_map(|f| a.fields.get(*f))
                .map_or(ABSENT, String::as_str);
            let row = (label, a.allele.as_str(), a.feature.as_str());
            if (blank(source) || source == REFSEQ_SOURCE) && !accepted_rows.contains(&row) {
                accepted_rows.insert(row);
                collapsed_rows.insert(&a.key);
                if let Some(record) = a.record {
                    collapsed_records.insert(record);
                }
            }
        }
    }
    // The built-in row classes of `TERM_PAIR_CLASSES`: the pairs are found before the entries
    // are compared, as the registry pairs are, since a JSON record's summary precedes its
    // consequence entries and its `most_severe_consequence` follows the term.
    let term_pair_classes: Vec<(&'static str, RowPair)> = TERM_PAIR_CLASSES
        .iter()
        .copied()
        .filter(|(k, _)| documented.masks.contains(*k))
        .collect();
    let mut term_pair_rows: HashMap<&str, &'static str> = HashMap::new();
    let mut term_pair_records: HashSet<i64> = HashSet::new();
    // The records that carry a row a class of a coding end took: their
    // `most_severe_consequence` is held to each side's own terms.
    let mut coding_end_records: HashSet<i64> = HashSet::new();
    if !term_pair_classes.is_empty() {
        for e in expected {
            let Some(a) = actual_by_key.get(e.key.as_str()) else {
                continue;
            };
            if record_key(&e.key) || manifest_reads(documented, e) {
                continue;
            }
            let first_taking = |at: &AtACodingEnd<'_>| -> Result<Option<&'static str>> {
                for (k, predicate) in &term_pair_classes {
                    if predicate(e, a, at)? {
                        return Ok(Some(*k));
                    }
                }
                Ok(None)
            };
            // No class takes two rows of one Consequence text, so the cache is read for a
            // row whose two texts differ alone.
            let klass = if e.consequence == a.consequence {
                first_taking(&AtACodingEnd::NO_CACHE)?
            } else {
                side.asked_at(e, first_taking)?
            };
            if let Some(klass) = klass {
                term_pair_rows.insert(&e.key, klass);
                if let Some(record) = e.record {
                    term_pair_records.insert(record);
                    if CodingEndClass::of_name(klass).is_some() {
                        coding_end_records.insert(record);
                    }
                }
            }
        }
    }
    // A row whose Consequence `refseq_edited_model_deleted_span_codons` takes moves its
    // record's `most_severe_consequence` with the term, as a row of a term-pair class does.
    if documented
        .masks
        .contains(REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS)
    {
        for e in expected {
            let Some(a) = actual_by_key.get(e.key.as_str()) else {
                continue;
            };
            if record_key(&e.key) || manifest_reads(documented, e) {
                continue;
            }
            if deleted_span_inframe_term_pair(e, a, &e.consequence, &a.consequence) {
                if let Some(record) = e.record {
                    term_pair_records.insert(record);
                }
            }
        }
    }
    // So does a row whose Consequence `use_transcript_ref_cache_cdna` takes: the term each
    // side's residues give is the record's most severe one wherever no other row gives a
    // term above it.
    if documented.masks.contains(USE_TRANSCRIPT_REF_CACHE_CDNA) {
        for e in expected {
            let Some(a) = actual_by_key.get(e.key.as_str()) else {
                continue;
            };
            if record_key(&e.key) || manifest_reads(documented, e) {
                continue;
            }
            if failed_edit_reference_cell(
                e,
                a,
                CONSEQUENCE_FIELDS[0],
                &e.consequence,
                &a.consequence,
            ) {
                if let Some(record) = e.record {
                    term_pair_records.insert(record);
                }
            }
        }
    }
    let most_severe = (!coding_end_records.is_empty()).then(|| {
        (
            most_severe_terms(expected, &coding_end_records),
            most_severe_terms(actual, &coding_end_records),
        )
    });
    // The built-in structural-variant classes: the structural-variant comparator's predicates
    // over the entries read as its tuples, found before the entries are compared as the registry
    // pairs and the collapsed rows are. A mate-read pair joins the registry pairs, being skipped
    // whole in the same way; a skipped record's rows and a lost own-end row are candidate-only
    // rows, set aside where the collapsed rows are. A key a manifest documents is the manifest's.
    let mut skipped = CandidateOnly::default();
    let mut own_rows = CandidateOnly::default();
    let mut mate_rows = CandidateOnly::default();
    if let Some(variants) = side
        .structural_variants
        .as_ref()
        .filter(|variants| !variants.classes.is_empty())
    {
        let expected_tuples = variants.tuples(expected);
        let actual_tuples = variants.tuples(actual);
        let perl_tuples: BTreeSet<StructuralTuple> = expected_tuples.values().cloned().collect();
        let rust_tuples: BTreeSet<StructuralTuple> = actual_tuples.values().cloned().collect();
        let candidate_only_rows =
            |documented: &Documented,
             excluded: &BTreeSet<StructuralTuple>,
             found: &mut CandidateOnly<'a>| {
                for a in actual {
                    if expected_keys.contains(a.key.as_str())
                        || documented.extra(a)
                        || !actual_tuples
                            .get(a.key.as_str())
                            .is_some_and(|tuple| excluded.contains(tuple))
                    {
                        continue;
                    }
                    found.rows.insert(&a.key);
                    if let Some(record) = a.record {
                        found.records.insert(record);
                    }
                }
            };
        if variants.classes.contains(REFERENCE_SKIPPED_RECORD) {
            let (excluded, _records) = structural::filter_reference_skipped_records(
                &perl_tuples,
                &rust_tuples,
                &variants.input_variants,
                MAX_SV_SIZE,
            );
            candidate_only_rows(documented, &excluded, &mut skipped);
        }
        if variants.classes.contains(BREAKEND_MATE_LOCAL_READ) {
            let (excluded_rust, excluded_perl, _kinds) =
                structural::filter_breakend_mate_local_read_pairs(
                    &perl_tuples,
                    &rust_tuples,
                    &variants.models,
                );
            for e in expected {
                let Some(a) = actual_by_key.get(e.key.as_str()) else {
                    continue;
                };
                if registry_pairs.contains_key(e.key.as_str()) || manifest_reads(documented, e) {
                    continue;
                }
                let reference_excluded = expected_tuples
                    .get(e.key.as_str())
                    .is_some_and(|tuple| excluded_perl.contains(tuple));
                let candidate_excluded = actual_tuples
                    .get(a.key.as_str())
                    .is_some_and(|tuple| excluded_rust.contains(tuple));
                if reference_excluded && candidate_excluded {
                    registry_pairs.insert(&e.key, BREAKEND_MATE_LOCAL_READ.to_string());
                    if let Some(record) = e.record {
                        registry_records.insert(record);
                    }
                }
            }
        }
        if variants.classes.contains(BND_OWN_CHROMOSOME_ROWS_LOST) {
            let excluded = structural::filter_bnd_own_rows_lost(
                &perl_tuples,
                &rust_tuples,
                &variants.by_chr,
                &variants.spans,
            );
            candidate_only_rows(documented, &excluded, &mut own_rows);
            let excluded = structural::filter_bnd_synonym_mate_rows(
                &perl_tuples,
                &rust_tuples,
                &variants.by_chr,
                &variants.spans,
            );
            candidate_only_rows(documented, &excluded, &mut mate_rows);
        }
    }
    // A record the reference wrote no row for whose every candidate row is a lost row of
    // `bnd_own_chromosome_rows_lost`: the reference's VCF line carries no CSQ at all.
    let lost = |key: &str| own_rows.rows.contains(key) || mate_rows.rows.contains(key);
    // The records with a reference row, or a candidate row the class keeps, in one pass.
    let written: HashSet<i64> = expected
        .iter()
        .chain(actual.iter().filter(|a| !lost(&a.key)))
        .filter(|e| !record_key(&e.key))
        .filter_map(|e| e.record)
        .collect();
    let rows_lost_whole: HashSet<i64> = mate_rows
        .records
        .iter()
        .copied()
        .filter(|record| !written.contains(record))
        .collect();
    let mut stale: Vec<String> = Vec::new();
    // For each manifest item that states `rows_alike`: the rows of its key that carry its
    // field, and those of them on which both engines print one value.
    let mut shared_rows: Vec<(i64, i64)> = vec![(0, 0); documented.shared_items.len()];
    let mut compared: i64 = 0;
    let mut skipped_documented: i64 = 0;
    let mut annotated_without_reference: i64 = 0;
    let mut registry_set_aside: i64 = 0;

    let mut individuals = Individuals {
        allele_order,
        individual_records: expected
            .iter()
            .filter(|e| record_key(&e.key) && carries_individual(e))
            .map(|e| key_record(&e.key))
            .collect(),
        expected_numbers: AlleleNumbers::new(),
        actual_numbers: AlleleNumbers::new(),
        order_only,
    };
    if allele_order {
        individuals.expected_numbers = individuals.allele_numbers(expected);
        individuals.actual_numbers = individuals.allele_numbers(actual);
    }
    let individuals = individuals;

    let mut classes_of_fields = FieldClasses {
        field_classes,
        record_index,
        candidate_index,
        record_entries,
        term_pair_rows,
        side,
        assembly: documented.assembly.clone(),
        dbnsfp_columns: documented.dbnsfp_columns.clone(),
        order_only,
    };

    // Counts the differing field `name` of the entry, an example kept under the
    // reference's consequence set.
    let mut count_difference = |e: &Entry, name: &str, ev: &str, av: &str| {
        mismatch_totals.add(name, 1);
        let group = mismatches
            .entry(name.to_string())
            .or_default()
            .entry(e.consequence.clone())
            .or_default();
        if takes_example(group, examples) {
            group.push(pair_example(&e.key, ev, av));
        }
    };

    for e in expected {
        let a = actual_by_key.get(e.key.as_str()).copied();
        let documented_sets = documented
            .divergent(e)
            .filter(|sets| !agrees_on_a_shared_key(e, a, sets));
        if let Some(sets) = documented_sets {
            skipped_documented += 1;
            match a {
                None => {
                    if !sets.is_empty() {
                        stale.push(format!(
                            "{}: documented candidate sets {}, but the candidate emits no row",
                            e.key,
                            repr_text_list(&sets)
                        ));
                    }
                }
                Some(a) => {
                    let got = consequence_set(&a.consequence);
                    // A VCF line and a JSON object print no consequence set of their own, so a
                    // record-level row the named run documents for the record (the one row of a
                    // `--summary` run) has nothing to observe on its record entry.
                    let record_entry = documented.run.is_some()
                        && record_key(&e.key)
                        && e.consequence.is_empty()
                        && got.is_empty();
                    if !record_entry && !sets.iter().any(|s| consequence_set(s) == got) {
                        stale.push(format!(
                            "{}: documented candidate sets {}, the candidate now emits {} (expected {})",
                            e.key,
                            repr_text_list(&sets),
                            pyfmt::repr_str(&a.consequence),
                            pyfmt::repr_str(&e.consequence)
                        ));
                    }
                }
            }
            continue;
        }
        let Some(a) = a else {
            missing.push(e);
            continue;
        };
        if let Some(registry_klass) = registry_pairs.get(e.key.as_str()) {
            registry_set_aside += 1;
            documented.applied.add(registry_klass, 1);
            let shown = class_examples.entry(registry_klass.clone()).or_default();
            if takes_example(shown, examples) {
                shown.push(pair_example(&e.key, &e.consequence, &a.consequence));
            }
            continue;
        }
        compared += 1;
        let summary_record = e.record.filter(|_| record_key(&e.key));
        let other_rows_change_the_summary = summary_record.is_some_and(|record| {
            documented.record_divergent(record)
                || registry_records.contains(&record)
                || collapsed_records.contains(&record)
                || own_rows.records.contains(&record)
                || mate_rows.records.contains(&record)
        });
        // A record that carries a row a class of a coding end took, and no row another
        // mechanism sets aside, has its `most_severe_consequence` derived on each side.
        let derived_summary = summary_record
            .filter(|record| !other_rows_change_the_summary && coding_end_records.contains(record));
        let record_summary_differs = other_rows_change_the_summary
            || (derived_summary.is_none()
                && summary_record.is_some_and(|record| term_pair_records.contains(&record)));
        // The default format's Extra is the raw text of the fields compared one by one
        // below, so it differs wherever one of them is documented to.
        let extra_differs = documented.has_divergent_fields(e);
        let stated_equal = documented.equal_with_reference(e);
        let mut names: Vec<&str> = e
            .fields
            .keys()
            .chain(a.fields.keys())
            .map(String::as_str)
            .collect();
        names.sort_unstable();
        names.dedup();
        for name in names {
            if record_summary_differs && name == "most_severe_consequence" {
                continue;
            }
            if name == "CSQ_present" && summary_record.is_some_and(|r| rows_lost_whole.contains(&r))
            {
                continue;
            }
            if extra_differs && name == "Extra" {
                continue;
            }
            let ev = e.fields.get(name).map_or(ABSENT, String::as_str);
            let av = a.fields.get(name).map_or(ABSENT, String::as_str);
            // A key without documented fields documents no value of this one.
            let documented_value = if extra_differs {
                documented.divergent_field(e, name)
            } else {
                None
            };
            if let Some(documented_value) = documented_value {
                // On a key whose item states `rows_alike`, a row on which both engines print
                // one value is another record's and agrees like any other row; the item's
                // count of such rows is held once every entry is compared.
                if let Some(item) = documented.shared_item(e, name) {
                    let (carrying, alike) = &mut shared_rows[item];
                    *carrying += 1;
                    if ev == av {
                        *alike += 1;
                        continue;
                    }
                }
                // The JSON format writes a numeric plugin value as a number, so the
                // documented default-format text is also matched with the trailing zeros
                // of its fractional part dropped (`1311092.0000` is written `1311092`).
                let printed = documented_form(av);
                let recorded = documented_form(&documented_value);
                let shortened = if FRACTIONAL_NUMBER.is_match(&recorded) {
                    recorded.trim_end_matches('0').trim_end_matches('.')
                } else {
                    recorded.as_str()
                };
                if printed != recorded && printed != shortened {
                    stale.push(format!(
                        "{}: field {} documented as {}, the candidate now prints {} (expected {})",
                        e.key,
                        name,
                        pyfmt::repr_str(&documented_value),
                        pyfmt::repr_str(av),
                        pyfmt::repr_str(ev)
                    ));
                }
                continue;
            }
            if ev == av {
                continue;
            }
            if let (Some(record), "most_severe_consequence") = (derived_summary, name) {
                let own = |terms: &HashMap<i64, String>, value: &str| {
                    terms.get(&record).is_some_and(|term| term == value)
                };
                let derived = most_severe
                    .as_ref()
                    .is_some_and(|(theirs, ours)| own(theirs, ev) && own(ours, av));
                if !derived {
                    count_difference(e, name, ev, av);
                }
                continue;
            }
            if stated_equal && field_named(name, CONSEQUENCE_FIELDS) {
                count_difference(e, name, ev, av);
                continue;
            }
            // `colocated_list_cut_at_cache_cap`: the reference's value is read as the candidate
            // reads its cut lists and compared on as that reading, so an `Extra` text with
            // such a pair beside a pair another class accepts meets both.
            let reading = if documented.masks.contains(COLOCATED_LIST_CUT_AT_CACHE_CAP) {
                cut_list_reading(e, name, ev, av)
            } else {
                None
            };
            if reading.is_some() {
                let klass = COLOCATED_LIST_CUT_AT_CACHE_CAP;
                class_columns
                    .entry(klass.to_string())
                    .or_default()
                    .add(name, 1);
                class_entries
                    .entry(klass.to_string())
                    .or_default()
                    .insert(&e.key);
                let shown = class_examples.entry(klass.to_string()).or_default();
                if takes_example(shown, examples) {
                    shown.push(cut_list_example(e, name, ev, av));
                }
            }
            let ev = reading.as_deref().unwrap_or(ev);
            if reading.is_some() && ev == av {
                continue;
            }
            if field_pairs
                .get(name)
                .is_some_and(|pairs| pairs.contains(&(ev.to_string(), av.to_string())))
                || (name == "Extra" && extra_pairs_accepted(ev, av, field_pairs))
            {
                pairs_accepted.add(name, 1);
                continue;
            }
            if let Some(sites) = nearest_sites {
                if (NEAREST_FIELDS.contains(&name) && sites.agree(ev, av))
                    || (name == "Extra" && extra_nearest_accepted(ev, av, sites))
                {
                    nearest_accepted.add(name, 1);
                    continue;
                }
            }
            if order_only
                .get(name)
                .is_some_and(|separator| same_items(ev, av, separator))
                || (name == "Extra"
                    && !order_only.is_empty()
                    && extra_order_only_accepted(ev, av, order_only))
            {
                order_only_accepted.add(name, 1);
                documented.applied.add(ORDER_ONLY, 1);
                let shown = order_only_examples.entry(name.to_string()).or_default();
                if takes_example(shown, examples) {
                    shown.push(pair_example(&e.key, ev, av));
                }
                continue;
            }
            if individuals.allele_order_accepts(e, a, name, ev, av) {
                allele_order_accepted.add(name, 1);
                documented.applied.add(INDIVIDUAL_ALLELE_ORDER, 1);
                let shown = allele_order_examples.entry(name.to_string()).or_default();
                if takes_example(shown, examples) {
                    shown.push(pair_example(&e.key, ev, av));
                }
                continue;
            }
            if reverse_motif_score_accepts(reverse_motif, e, name, ev, av) {
                reverse_motif_accepted.add(name, 1);
                documented.applied.add(MOTIF_SCORE_REVERSE_STRAND_MNV, 1);
                let shown = reverse_motif_examples.entry(name.to_string()).or_default();
                if takes_example(shown, examples) {
                    shown.push(pair_example(&e.key, ev, av));
                }
                continue;
            }
            let term_pair_class = classes_of_fields
                .term_pair_rows
                .get(e.key.as_str())
                .copied()
                .filter(|_| field_named(name, CONSEQUENCE_FIELDS));
            let accepting: Vec<&str> = match term_pair_class {
                Some(klass) => vec![klass],
                None => classes_of_fields.field_class_accepting(
                    e,
                    a,
                    name,
                    ev,
                    av,
                    &mut documented.needed,
                )?,
            };
            if !accepting.is_empty() {
                for klass in accepting {
                    class_columns
                        .entry(klass.to_string())
                        .or_default()
                        .add(name, 1);
                    class_entries
                        .entry(klass.to_string())
                        .or_default()
                        .insert(&e.key);
                    let shown = class_examples.entry(klass.to_string()).or_default();
                    if takes_example(shown, examples) {
                        shown.push(format!(
                            "{}: {} expected {}, got {}",
                            e.key,
                            name,
                            pyfmt::repr_str(ev),
                            pyfmt::repr_str(av)
                        ));
                    }
                }
                continue;
            }
            count_difference(e, name, ev, av);
        }
    }
    // The built-in row shape of `frequency_filter_zero_value`: a record the candidate alone
    // writes, every consequence row of which carries a FREQS of zeros, is set aside with its
    // record entry, as a manifest's documented tuples of the class are.
    let mut zero_frequency_records: HashSet<String> = HashSet::new();
    if documented.masks.contains(FREQUENCY_FILTER_ZERO_VALUE) {
        let reference_labels: HashSet<String> = expected.iter().map(record_label).collect();
        let mut candidate_only: IndexMap<String, Vec<&Entry>> = IndexMap::new();
        for a in actual {
            if !expected_keys.contains(a.key.as_str())
                && !documented.extra(a)
                && !record_key(&a.key)
            {
                candidate_only.entry(record_label(a)).or_default().push(a);
            }
        }
        // A record the reference wrote any row for is kept by both engines, so an extra
        // candidate row on it is a difference whatever its FREQS says.
        zero_frequency_records = candidate_only
            .into_iter()
            .filter(|(label, rows)| {
                !reference_labels.contains(label)
                    && !rows.is_empty()
                    && rows.iter().all(|a| freqs_all_zero(a))
            })
            .map(|(label, _)| label)
            .collect();
    }
    for a in actual {
        if expected_keys.contains(a.key.as_str()) || documented.extra(a) {
            continue;
        }
        // Only the record-level line or object is without a counterpart; a consequence
        // entry on a skipped record is a candidate-only key, documented, built in or unexpected.
        if record_key(&a.key)
            && a.record.is_some_and(|record| {
                documented.reference_skipped(record) || skipped.records.contains(&record)
            })
        {
            annotated_without_reference += 1;
            continue;
        }
        // No label is one of an empty set.
        if !zero_frequency_records.is_empty() && zero_frequency_records.contains(&record_label(a)) {
            if record_key(&a.key) {
                annotated_without_reference += 1;
            } else {
                class_entries
                    .entry(FREQUENCY_FILTER_ZERO_VALUE.to_string())
                    .or_default()
                    .insert(&a.key);
                let shown = class_examples
                    .entry(FREQUENCY_FILTER_ZERO_VALUE.to_string())
                    .or_default();
                if takes_example(shown, examples) {
                    let freqs = FREQS_FIELDS
                        .iter()
                        .find_map(|f| a.fields.get(*f))
                        .map(String::as_str);
                    shown.push(format!(
                        "{}: candidate-only row, FREQS {}",
                        a.key,
                        pyfmt::repr_opt_str(freqs)
                    ));
                }
            }
            continue;
        }
        let set_aside = if collapsed_rows.contains(a.key.as_str()) {
            Some((MERGED_SAME_STABLE_ID_COLLAPSE, "candidate-only row"))
        } else if skipped.rows.contains(a.key.as_str()) {
            Some((
                REFERENCE_SKIPPED_RECORD,
                "candidate-only row on a record the reference skipped",
            ))
        } else if own_rows.rows.contains(a.key.as_str()) {
            Some((BND_OWN_CHROMOSOME_ROWS_LOST, "candidate-only own-end row"))
        } else if mate_rows.rows.contains(a.key.as_str()) {
            Some((BND_OWN_CHROMOSOME_ROWS_LOST, "candidate-only mate row"))
        } else {
            None
        };
        if let Some((klass, what)) = set_aside {
            class_entries
                .entry(klass.to_string())
                .or_default()
                .insert(&a.key);
            let shown = class_examples.entry(klass.to_string()).or_default();
            if takes_example(shown, examples) {
                shown.push(format!(
                    "{}: {}, {}",
                    a.key,
                    what,
                    pyfmt::repr_str(&a.consequence)
                ));
            }
            continue;
        }
        unexpected.push(a);
    }

    let shared_expected: Vec<&Entry> = expected
        .iter()
        .filter(|e| actual_by_key.contains_key(e.key.as_str()) && documented.divergent(e).is_none())
        .collect();
    let shared_actual: Vec<&Entry> = actual
        .iter()
        .filter(|a| expected_keys.contains(a.key.as_str()) && documented.divergent(a).is_none())
        .collect();
    let mut order = Value::Null;
    let mut order_differs_at: Option<i64> = None;
    let mut reordered: i64 = 0;
    if options.check_order {
        let expected_order: Vec<String> = shared_expected
            .iter()
            .map(|e| individuals.order_key(e))
            .collect();
        let actual_order: Vec<String> = shared_actual
            .iter()
            .map(|a| individuals.order_key(a))
            .collect();
        if expected_order != actual_order {
            // Two lists of which one opens the other differ at no shared place, and the
            // place named is then their first.
            let first = expected_order
                .iter()
                .zip(&actual_order)
                .position(|(x, y)| x != y)
                .unwrap_or(0);
            // The front end keeps back no row between two rows of one handed-over group, on
            // either side, and keeps its rows in one order on both, so a shared key of the
            // remainder has the same kept-back rows before it on both sides: the remainder's
            // first difference is the outputs' first, moved by those rows.
            let before = match (kept, shared_expected.get(first)) {
                (Some(kept), Some(e)) => kept.before(&e.key),
                _ => 0,
            };
            let key_at =
                |entries: &[&Entry]| entries.get(first).map_or(Value::Null, |e| text(&e.key));
            let place = i64::try_from(first).unwrap_or(i64::MAX) + before;
            order_differs_at = Some(place);
            order = dict([
                ("first_divergence", int(place)),
                ("expected", key_at(&shared_expected)),
                ("actual", key_at(&shared_actual)),
            ]);
        } else {
            let moved = shared_expected
                .iter()
                .zip(&shared_actual)
                .filter(|(e, a)| e.key != a.key)
                .count();
            reordered = i64::try_from(moved).unwrap_or(i64::MAX);
            if reordered != 0 {
                documented.applied.add(INDIVIDUAL_ALLELE_ORDER, reordered);
            }
        }
    }

    for (klass, keys) in &class_entries {
        documented
            .applied
            .add(klass, i64::try_from(keys.len()).unwrap_or(i64::MAX));
    }

    let mut expected_populated = populated_counts(expected);
    let mut actual_populated = populated_counts(actual);
    let mut kept_entries: i64 = 0;
    if let Some(kept) = kept {
        kept_entries = kept.entries;
        expected_populated.update(&kept.populated);
        actual_populated.update(&kept.populated);
    }
    let mut columns: IndexMap<String, Value> = IndexMap::new();
    // Per floor's column, whether the floor is met and the two counts it was read against.
    let mut floors_read: HashMap<&str, (Option<bool>, i64, i64)> = HashMap::new();
    let mut floors_not_applicable: Vec<String> = Vec::new();
    let mut column_names: BTreeSet<&str> = BTreeSet::new();
    column_names.extend(expected_populated.iter().map(|(name, _)| name));
    column_names.extend(actual_populated.iter().map(|(name, _)| name));
    column_names.extend(mismatch_totals.iter().map(|(name, _)| name));
    column_names.extend(floors.keys().map(String::as_str));
    for name in column_names {
        let mut expected_count = expected_populated.get(name);
        let mut actual_count = actual_populated.get(name);
        let floor = floors.get(name).copied();
        let mut floor_met: Option<bool> = None;
        let mut floor_detail: Option<(&str, Value)> = None;
        if let Some(floor) = floor {
            match floor_fields.get(name) {
                None => floor_met = Some(expected_count >= floor && actual_count >= floor),
                Some(None) => {
                    floors_not_applicable.push(name.to_string());
                    floor_detail = Some(("floor_note", text(FLOOR_NOT_APPLICABLE)));
                }
                Some(Some(StrListOrStr::Str(note))) => {
                    floors_not_applicable.push(name.to_string());
                    floor_detail = Some(("floor_note", text(note)));
                }
                Some(Some(StrListOrStr::StrList(counted))) => {
                    let counted_fields: Vec<&str> = counted.iter().map(String::as_str).collect();
                    let under = kept.map_or(0, |kept| kept.counted_under(name, counted));
                    expected_count = populated_under(expected, &counted_fields) + under;
                    actual_count = populated_under(actual, &counted_fields) + under;
                    floor_detail = Some(("floor_fields", texts(counted)));
                    floor_met = Some(expected_count >= floor && actual_count >= floor);
                }
            }
            floors_read.insert(name, (floor_met, expected_count, actual_count));
        }
        let mut col: IndexMap<String, Value> = IndexMap::new();
        let mut put = |key: &str, value: Value| {
            col.insert(key.to_string(), value);
        };
        put("expected_populated", int(expected_count));
        put("actual_populated", int(actual_count));
        put("differing", int(mismatch_totals.get(name)));
        put("pairs_accepted", int(pairs_accepted.get(name)));
        put("nearest_accepted", int(nearest_accepted.get(name)));
        put("order_only", int(order_only_accepted.get(name)));
        put(
            "individual_allele_order",
            int(allele_order_accepted.get(name)),
        );
        put(
            "motif_score_reverse_strand_mnv",
            int(reverse_motif_accepted.get(name)),
        );
        put("floor", floor.map_or(Value::Null, int));
        put("floor_met", floor_met.map_or(Value::Null, Value::Bool));
        if let Some((key, value)) = floor_detail {
            put(key, value);
        }
        if let Some(groups) = mismatches.get(name) {
            put(
                "examples",
                Value::Dict(
                    groups
                        .iter()
                        .map(|(consequence, shown)| (consequence.clone(), texts(shown)))
                        .collect(),
                ),
            );
        }
        if let Some(shown) = order_only_examples.get(name) {
            put("order_only_examples", texts(shown));
        }
        if let Some(shown) = allele_order_examples.get(name) {
            put("individual_allele_order_examples", texts(shown));
        }
        if let Some(shown) = reverse_motif_examples.get(name) {
            put("motif_score_reverse_strand_mnv_examples", texts(shown));
        }
        let mut masked: Vec<(&str, i64)> = class_columns
            .iter()
            .map(|(klass, counts)| (klass.as_str(), counts.get(name)))
            .filter(|(_, counted)| *counted != 0)
            .collect();
        masked.sort_unstable();
        if !masked.is_empty() {
            put(
                "masked",
                Value::Dict(
                    masked
                        .into_iter()
                        .map(|(klass, counted)| (klass.to_string(), int(counted)))
                        .collect(),
                ),
            );
        }
        columns.insert(name.to_string(), Value::Dict(col));
    }

    // An item that states `rows_alike` documents that count of rows on its key.
    stale.extend(documented.stale_row_counts(&shared_rows));
    let mut reasons: Vec<String> = Vec::new();
    if !missing.is_empty() {
        reasons.push(format!(
            "missing in the candidate: {} entries",
            missing.len()
        ));
    }
    if !unexpected.is_empty() {
        reasons.push(format!(
            "unexpected in the candidate: {} entries",
            unexpected.len()
        ));
    }
    let mut differing: Vec<(&str, i64)> = mismatch_totals.iter().collect();
    differing.sort_unstable();
    for (name, total) in differing {
        reasons.push(format!("field {name}: {total} differing entries"));
    }
    if !stale.is_empty() {
        reasons.push(format!(
            "documented divergences not observed: {}",
            stale.len()
        ));
    }
    if let Some(place) = order_differs_at {
        reasons.push(format!("order differs at shared position {place}"));
    }
    let mut floor_names: Vec<&str> = floors.keys().map(String::as_str).collect();
    floor_names.sort_unstable();
    for name in floor_names {
        if let Some((Some(false), expected_count, actual_count)) = floors_read.get(name).copied() {
            reasons.push(format!(
                "floor {}={}: expected populated {}, actual populated {}",
                name,
                floors.get(name).copied().unwrap_or_default(),
                expected_count,
                actual_count
            ));
        }
    }

    // The rows the class of a record of several entries was asked about under the score files.
    let refusals = spliceai_table.map(|table| table.refusals.counts());
    let mut mask_examples: Vec<(&String, &Vec<String>)> = class_examples.iter().collect();
    mask_examples.sort_unstable_by(|a, b| a.0.cmp(b.0));
    stale.truncate(STALE_LISTED);
    let mut report: IndexMap<String, Value> = IndexMap::new();
    let mut put = |key: &str, value: Value| {
        report.insert(key.to_string(), value);
    };
    put(
        "verdict",
        text(if reasons.is_empty() { "green" } else { "red" }),
    );
    put("reasons", texts(&reasons));
    put(
        "counts",
        dict([
            (
                "expected_entries",
                int(i64::try_from(expected.len()).unwrap_or(i64::MAX) + kept_entries),
            ),
            (
                "actual_entries",
                int(i64::try_from(actual.len()).unwrap_or(i64::MAX) + kept_entries),
            ),
            ("compared", int(compared + kept_entries)),
            ("documented_divergences_skipped", int(skipped_documented)),
            (
                "record_entries_on_reference_skipped_records",
                int(annotated_without_reference),
            ),
            ("individual_entries_reordered", int(reordered)),
            ("registry_pairs_set_aside", int(registry_set_aside)),
        ]),
    );
    put(
        "one_sided",
        dict([
            ("expected_only", count(missing.len())),
            ("actual_only", count(unexpected.len())),
        ]),
    );
    put("missing", grouped(&missing, examples));
    put("unexpected", grouped(&unexpected, examples));
    put("columns", Value::Dict(columns));
    put(
        "mask_examples",
        Value::Dict(
            mask_examples
                .into_iter()
                .map(|(klass, shown)| (klass.clone(), texts(shown)))
                .collect(),
        ),
    );
    put("floors_not_applicable", texts(&floors_not_applicable));
    put("stale_divergences", texts(&stale));
    put("order", order);
    if let Some(refusals) = refusals {
        if truthy(&refusals) && refusals.get("asked").is_some_and(truthy) {
            put("spliceai_refusals", refusals);
        }
    }
    Ok(Value::Dict(report))
}
