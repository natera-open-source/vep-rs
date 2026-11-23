// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One output row as the class rules read it: its fields in the order the row writes them,
//! its key, and what an identifier-less `Uploaded_variation` says of the record's alleles.

use super::shared::{dash, integer_of};
use super::Refusal;
use crate::pytext::{self, Json, JsonError};
use regex::Regex;
use std::cell::RefCell;
use std::collections::{BTreeSet, HashMap};
use std::rc::Rc;
use std::sync::LazyLock;

/// The field value of an entry that lacks the field.
pub const ABSENT: &str = "<absent>";

/// The fields whose text is a JSON object, compared in canonical form.
pub const JSON_OBJECT_FIELDS: [&str; 1] = ["GA4GH_VRS"];

/// The fields of an entry, a Python `dict` of text to text: the names in the order of their
/// first insertion, a second insertion of a name replacing its value in place.
#[derive(Clone, Default)]
pub struct Fields {
    /// Every name and value held, back to back.
    text: String,
    items: Vec<Item>,
    /// Open addressing over `items` for more names than `SCANNED`: an item's index plus
    /// one, zero for a free slot.
    slots: Vec<u32>,
}

/// One field: where its name and its value stand in the text of its `Fields`.
#[derive(Clone, Copy)]
struct Item {
    name: usize,
    name_end: usize,
    value: usize,
    value_end: usize,
}

/// The most names looked up by reading them all.
const SCANNED: usize = 8;

fn name_hash(name: &str) -> usize {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in name.bytes() {
        hash = (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
    }
    (hash ^ (hash >> 29)) as usize
}

impl Fields {
    pub fn new() -> Fields {
        Fields::default()
    }

    /// Fields with room for `names` names and `bytes` bytes of names and values.
    pub fn with_capacity(names: usize, bytes: usize) -> Fields {
        Fields {
            text: String::with_capacity(bytes),
            items: Vec::with_capacity(names),
            slots: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    fn name(&self, at: usize) -> &str {
        &self.text[self.items[at].name..self.items[at].name_end]
    }

    fn value(&self, at: usize) -> &str {
        &self.text[self.items[at].value..self.items[at].value_end]
    }

    fn position(&self, name: &str) -> Option<usize> {
        if self.slots.is_empty() {
            return (0..self.items.len()).find(|&at| self.name(at) == name);
        }
        let mask = self.slots.len() - 1;
        let mut slot = name_hash(name) & mask;
        loop {
            match self.slots[slot] {
                0 => return None,
                held if self.name(held as usize - 1) == name => return Some(held as usize - 1),
                _ => slot = (slot + 1) & mask,
            }
        }
    }

    fn index(&mut self) {
        let size = (self.items.len() * 4).next_power_of_two();
        self.slots = vec![0; size];
        for at in 0..self.items.len() {
            self.place(at);
        }
    }

    fn place(&mut self, at: usize) {
        let mask = self.slots.len() - 1;
        let mut slot = name_hash(self.name(at)) & mask;
        while self.slots[slot] != 0 {
            slot = (slot + 1) & mask;
        }
        self.slots[slot] = at as u32 + 1;
    }

    fn push(&mut self, name: &str, value: &str) {
        let start = self.text.len();
        self.text.push_str(name);
        self.text.push_str(value);
        let name_end = start + name.len();
        self.items.push(Item {
            name: start,
            name_end,
            value: name_end,
            value_end: self.text.len(),
        });
        if self.items.len() <= SCANNED {
            return;
        }
        if self.items.len() * 2 > self.slots.len() {
            self.index();
        } else {
            self.place(self.items.len() - 1);
        }
    }

    /// `fields[name] = value`.
    pub fn insert(&mut self, name: &str, value: impl AsRef<str>) {
        let value = value.as_ref();
        match self.position(name) {
            Some(at) => {
                let start = self.text.len();
                self.text.push_str(value);
                (self.items[at].value, self.items[at].value_end) = (start, self.text.len());
            }
            None => self.push(name, value),
        }
    }

    /// `fields.setdefault(name, value)`: the value is kept where the name is held.
    pub fn set_default(&mut self, name: &str, value: impl AsRef<str>) {
        if self.position(name).is_none() {
            self.push(name, value.as_ref());
        }
    }

    /// `fields.get(name)`.
    pub fn get(&self, name: &str) -> Option<&str> {
        self.position(name).map(|at| self.value(at))
    }

    /// `fields.get(name, default)`.
    pub fn get_or<'a>(&'a self, name: &str, default: &'a str) -> &'a str {
        self.get(name).unwrap_or(default)
    }

    /// `name in fields`.
    pub fn contains(&self, name: &str) -> bool {
        self.position(name).is_some()
    }

    /// `fields.items()`: the names and values in the order of insertion.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &str)> {
        (0..self.items.len()).map(|at| (self.name(at), self.value(at)))
    }

    /// The value of the first of `names` the entry holds (`next((fields[f] for f in names
    /// if f in fields), default)`).
    pub fn first_of<'a>(&'a self, names: &[&str], default: &'a str) -> &'a str {
        names
            .iter()
            .find_map(|name| self.get(name))
            .unwrap_or(default)
    }
}

/// Equal where both hold the same names with the same values in the same order.
impl PartialEq for Fields {
    fn eq(&self, other: &Fields) -> bool {
        self.len() == other.len() && self.iter().eq(other.iter())
    }
}

/// The names and values as a map in the order of insertion.
impl std::fmt::Debug for Fields {
    fn fmt(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.debug_map().entries(self.iter()).finish()
    }
}

impl<K: AsRef<str>, V: AsRef<str>> FromIterator<(K, V)> for Fields {
    fn from_iter<I: IntoIterator<Item = (K, V)>>(pairs: I) -> Fields {
        let mut fields = Fields::new();
        for (name, value) in pairs {
            fields.insert(name.as_ref(), value);
        }
        fields
    }
}

/// One row of the default or tab format, or an entry a caller builds in a row's place (the
/// record entry of an input line, a tuple file's values under one key).
///
/// `key` is the entry's identity within its output: the Location, the allele, the feature
/// and the feature type joined by `|`, plus a `#n` counter on exact duplicates.
#[derive(Debug, Clone, Default)]
pub struct Entry {
    pub key: String,
    pub record: Option<i64>,
    pub location: Option<String>,
    pub allele: String,
    pub feature: String,
    pub consequence: String,
    pub fields: Fields,
    /// Whether the allele is an insertion's or a deletion's, read from the record's own
    /// alleles (`allele_is_indel` on an identifier-less `Uploaded_variation`); `None` where
    /// they say nothing.
    pub indel: Option<bool>,
    /// Whether the same alleles make every sequence allele an insertion
    /// (`reference_allele_is_empty`); `None` where the entry's record spells no alleles.
    pub insertion: Option<bool>,
    /// The reference bases the allele stands in for at the entry's Location
    /// (`replaced_reference`), read from an identifier-less `Uploaded_variation` by the
    /// parser and filled by a caller that reads the record's input line; empty for a pure
    /// insertion's entry, `None` where no record was read.
    pub reference: RefCell<Option<String>>,
    pub line: Option<i64>,
    pub allele_location: Option<String>,
}

/// A shared entry. Two rows are the same entry where `Rc::ptr_eq` says so.
pub type Row = Rc<Entry>;

/// The lines of `text` without their terminators, no trailing empty line.
pub fn text_lines(text: &str) -> Vec<&str> {
    let mut lines: Vec<&str> = text.split('\n').collect();
    if lines.last() == Some(&"") {
        lines.pop();
    }
    for line in &mut lines {
        *line = line.strip_suffix('\r').unwrap_or(line);
    }
    lines
}

/// Suffixes the second and later occurrences of a key with `#n`.
pub fn dedup_keys(entries: &mut [Entry]) {
    let mut seen: HashMap<String, u64> = HashMap::new();
    for e in entries {
        let n = seen.entry(e.key.clone()).or_insert(0);
        *n += 1;
        if *n > 1 {
            e.key = format!("{}#{n}", e.key);
        }
    }
}

/// A JSON object's text with its members sorted (`pytext::json_compact`), or `text` as it
/// is when it is not a JSON object. Ensembl VEP encodes the `GA4GH_VRS` object with its
/// keys in hash order, which differs from row to row, so the parsed object is compared. An
/// object nested deeper than `pytext::JSON_NESTING` is refused.
pub fn canonical_json_text(text: &str) -> Result<String, Refusal> {
    if !text
        .trim_start_matches([' ', '\t', '\n', '\r'])
        .starts_with('{')
    {
        return Ok(text.to_string());
    }
    match pytext::json_loads(text) {
        Ok(value @ Json::Object(_)) => Ok(pytext::json_compact(&value)),
        Ok(_) | Err(JsonError::Invalid) => Ok(text.to_string()),
        Err(JsonError::TooDeep) => {
            let shown: String = text.chars().take(24).collect();
            Err(Refusal(format!(
                "a JSON object is nested too deep to read: {shown}..."
            )))
        }
    }
}

/// The default format: the columns of the `#` header line, the `Extra` column's
/// `key=value` pairs becoming fields of their own while `Extra` keeps its raw text, a
/// `JSON_OBJECT_FIELDS` value in canonical form in both. A row with another number of
/// columns than the header is refused.
pub fn parse_default(text: &str) -> Result<Vec<Row>, Refusal> {
    let mut cols: Vec<&str> = Vec::new();
    let mut entries: Vec<Entry> = Vec::new();
    for line in text_lines(text) {
        if let Some(header) = line.strip_prefix('#') {
            if !header.starts_with('#') {
                cols = header.split('\t').collect();
            }
            continue;
        }
        entries.push(default_entry(&cols, line)?);
    }
    dedup_keys(&mut entries);
    Ok(entries.into_iter().map(Rc::new).collect())
}

/// The column header of a default-format output, under which its rows are read one at a
/// time.
pub struct DefaultHeader {
    line: String,
    /// The column names, where the line is one column header line.
    columns: Option<Vec<String>>,
}

impl DefaultHeader {
    pub fn new(line: &str) -> DefaultHeader {
        let header = line.strip_suffix('\r').unwrap_or(line);
        let columns = header
            .strip_prefix('#')
            .filter(|names| !names.starts_with('#') && !line.contains('\n'))
            .map(|names| names.split('\t').map(str::to_string).collect());
        DefaultHeader {
            line: line.to_string(),
            columns,
        }
    }

    /// The entry of `row`, one line of the output: the first entry `parse_default` reads
    /// of the header line and the row, `None` where it reads none.
    pub fn row(&self, row: &str) -> Result<Option<Row>, Refusal> {
        let plain = !row.is_empty() && !row.contains('\n');
        let line = row.strip_suffix('\r').unwrap_or(row);
        match &self.columns {
            Some(columns) if plain && !line.starts_with('#') => {
                Ok(Some(Rc::new(default_entry(columns, line)?)))
            }
            _ => Ok(parse_default(&format!("{}\n{row}", self.line))?
                .into_iter()
                .next()),
        }
    }
}

/// One row of the default format under the column names `cols`.
fn default_entry(cols: &[impl AsRef<str>], line: &str) -> Result<Entry, Refusal> {
    let vals: Vec<&str> = line.split('\t').collect();
    if vals.len() != cols.len() {
        return Err(Refusal(format!("column count on line: {line}")));
    }
    // A row's fields are its columns and the pairs of its Extra column, whose text
    // the names and values are cut from.
    let pairs = vals
        .last()
        .map_or(0, |extra| extra.matches(';').count() + 1);
    let named: usize = cols.iter().map(|name| name.as_ref().len()).sum();
    let mut fields = Fields::with_capacity(cols.len() + pairs, named + 2 * line.len());
    for (name, value) in cols.iter().zip(&vals) {
        fields.insert(name.as_ref(), *value);
    }
    for name in JSON_OBJECT_FIELDS {
        if let Some(value) = fields.get(name) {
            let canonical = canonical_json_text(value)?;
            fields.insert(name, canonical);
        }
    }
    if let Some(extra) = fields.get("Extra").filter(|extra| *extra != "-") {
        let extra = extra.to_string();
        let mut rewritten: Option<Vec<String>> = None;
        for (number, pair) in extra.split(';').enumerate() {
            let (name, sep, value) = pytext::partition(pair, "=");
            if !sep.is_empty() && JSON_OBJECT_FIELDS.contains(&name) {
                let canonical = canonical_json_text(value)?;
                let pairs = rewritten.get_or_insert_with(|| {
                    extra.split(';').take(number).map(str::to_string).collect()
                });
                pairs.push(format!("{name}={canonical}"));
                fields.insert(name, canonical);
            } else {
                if let Some(pairs) = rewritten.as_mut() {
                    pairs.push(pair.to_string());
                }
                fields.insert(name, value);
            }
        }
        // A pair named `Extra` has taken the column's place among the fields.
        let kept = rewritten.map_or(extra, |pairs| pairs.join(";"));
        if fields.get("Extra") != Some(kept.as_str()) {
            fields.insert("Extra", kept);
        }
    }
    let location = fields.get_or("Location", "").to_string();
    let allele = fields.get_or("Allele", "").to_string();
    let feature = fields.get_or("Feature", "").to_string();
    let (mut indel, mut insertion, mut reference) = (None, None, None);
    let uploaded = fields.get_or("Uploaded_variation", "");
    if let Some((position, alleles)) = uploaded_parts(uploaded) {
        let mut alleles = alleles.split('/');
        let reference_allele = alleles.next().unwrap_or("");
        let alternate_alleles: Vec<&str> = alleles.collect();
        indel = allele_is_indel(reference_allele, &alternate_alleles, &allele, false);
        insertion = Some(reference_allele_is_empty(
            reference_allele,
            &alternate_alleles,
            false,
        ));
        let position = integer_of(position, "Uploaded_variation position")?;
        reference = replaced_reference(
            position,
            reference_allele,
            &alternate_alleles,
            Some(&location),
            &allele,
        )?;
    }
    let key = [
        location.as_str(),
        allele.as_str(),
        feature.as_str(),
        fields.get_or("Feature_type", ""),
    ]
    .join("|");
    Ok(Entry {
        key,
        record: None,
        location: Some(location),
        allele,
        feature,
        consequence: fields.get_or("Consequence", "").to_string(),
        fields,
        indel,
        insertion,
        reference: RefCell::new(reference),
        line: None,
        allele_location: None,
    })
}

/// Whether a field value is empty in any format's spelling: `-`, the empty string or absent
/// from the entry.
pub fn blank(value: &str) -> bool {
    matches!(value, "" | "-" | ABSENT)
}

/// Undoes the percent-encoding of `=` and `;` in Extra values.
pub fn decode_extra(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    s.replace("%3D", "=").replace("%3B", ";")
}

/// Whether a field name is one of `names`, each a full name or a prefix ending in `_`, the
/// case of both set aside.
pub fn field_named(name: &str, names: &[&str]) -> bool {
    if name.is_ascii() {
        return names.iter().any(|pattern| {
            if pattern.ends_with('_') {
                name.len() >= pattern.len()
                    && name.as_bytes()[..pattern.len()].eq_ignore_ascii_case(pattern.as_bytes())
            } else {
                name.eq_ignore_ascii_case(pattern)
            }
        });
    }
    let lowered = name.to_lowercase();
    names.iter().any(|pattern| {
        let pattern_lowered = pattern.to_lowercase();
        if pattern.ends_with('_') {
            lowered.starts_with(&pattern_lowered)
        } else {
            lowered == pattern_lowered
        }
    })
}

/// The value of the entry's first field named by one of `names` (`field_named`), `ABSENT`
/// for none.
pub fn field_of<'a>(e: &'a Entry, names: &[&str]) -> &'a str {
    e.fields
        .iter()
        .find_map(|(name, value)| field_named(name, names).then_some(value))
        .unwrap_or(ABSENT)
}

/// `field_of` for an entry that may be missing, `ABSENT` where it is.
pub fn field_of_any<'a>(e: Option<&'a Entry>, names: &[&str]) -> &'a str {
    e.map_or(ABSENT, |e| field_of(e, names))
}

/// Consequence terms as a sorted, de-duplicated, comma-joined set.
pub fn consequence_set(terms: &str) -> String {
    let terms: BTreeSet<&str> = terms
        .split(['&', ','])
        .filter(|term| !term.is_empty())
        .collect();
    terms.into_iter().collect::<Vec<_>>().join(",")
}

/// The last segment of an entry's key without its duplicate counter: the Feature_type, or
/// `record` for a record entry.
pub fn entry_kind(e: &Entry) -> &str {
    let last = e.key.rsplit('|').next().unwrap_or("");
    last.split('#').next().unwrap_or("")
}

/// The consequence list a feature type's entries go to in Ensembl VEP's JSON output
/// (`Transcript` to `transcript`, `RegulatoryFeature` to `regulatory_feature`, an undefined
/// type to `intergenic`); a list name is returned unchanged.
pub fn feature_kind(feature_type: &str) -> String {
    let kind = if matches!(feature_type, "" | "-") {
        "intergenic".to_string()
    } else {
        feature_type.to_lowercase()
    };
    if kind == "intergenic" || kind.contains("_feature") {
        return kind;
    }
    kind.replace("feature", "_feature")
}

/// The value as a finite number, or `None` when it is not one.
pub fn as_number(value: &str) -> Option<f64> {
    pytext::float(value).filter(|number| number.is_finite())
}

/// Ensembl VEP's Uploaded_variation for a VCF record without an ID (OutputFactory.pm:883-931):
/// the chromosome, the position and the record's own alleles joined by `/`. The groups are
/// the position and the alleles; the end is read before a final line feed, as Python's `$`
/// reads it.
static UPLOADED: LazyLock<Regex> = LazyLock::new(|| {
    let space = pytext::RE_SPACE;
    Regex::new(&format!(
        r"^[^{space}]+_([0-9]+)_([ACGTNacgtn*-]+(?:/[^/{space}]+)+)\n?$"
    ))
    .expect("regex")
});

/// The position and the `/`-joined alleles an identifier-less `Uploaded_variation` spells.
fn uploaded_parts(uploaded: &str) -> Option<(&str, &str)> {
    if !uploaded.contains('/') {
        return None;
    }
    let found = UPLOADED.captures(uploaded)?;
    Some((found.get(1)?.as_str(), found.get(2)?.as_str()))
}

/// The reference and alternate alleles an identifier-less `Uploaded_variation` spells
/// (`21_10542414_A/T/AG`), or `None` for an identifier, which names no allele.
pub fn uploaded_alleles(uploaded: &str) -> Option<(String, Vec<String>)> {
    let (_, alleles) = uploaded_parts(uploaded)?;
    let mut alleles = alleles.split('/').map(str::to_string);
    let reference = alleles.next()?;
    Some((reference, alleles.collect()))
}

/// The POS an identifier-less `Uploaded_variation` spells before its alleles, or `None` for
/// an identifier. A position `integer_of` does not read is refused.
pub fn uploaded_position(uploaded: &str) -> Result<Option<i128>, Refusal> {
    match uploaded_parts(uploaded) {
        Some((position, _)) => integer_of(position, "Uploaded_variation position").map(Some),
        None => Ok(None),
    }
}

fn length(text: &str) -> usize {
    if text.is_ascii() {
        text.len()
    } else {
        text.chars().count()
    }
}

fn stripped_length(allele: &str) -> usize {
    length(allele.trim_matches('-'))
}

fn without_first(allele: &str) -> &str {
    let mut chars = allele.chars();
    chars.next();
    chars.as_str()
}

/// Whether the record's shared first base is dropped as `Parser/VCF.pm:300-323` drop it: the
/// record holds an ALT of another length than its REF and every alternate not spelt with a
/// `*` starts with the reference's first base. The parts of an allele string already so
/// spelt (`chopped`) drop none.
pub(super) fn drops_first_base<S: AsRef<str>>(
    reference: &str,
    alternates: &[S],
    chopped: bool,
) -> bool {
    let reference_length = length(reference);
    !chopped
        && !reference.is_empty()
        && alternates
            .iter()
            .any(|a| length(a.as_ref()) != reference_length)
        && alternates
            .iter()
            .map(AsRef::as_ref)
            .filter(|a| !a.contains('*'))
            .all(|a| a.chars().next() == reference.chars().next())
}

/// Whether `allele` is an insertion or a deletion of the record whose reference allele is
/// `reference` and alternate alleles `alternates`: the alternate the allele names and the
/// reference differ in length. A record with one alternate not spelt with a `*` has that
/// alternate's shape on every row, since Ensembl VEP shortens a lone allele past the shared
/// base; with several, the shared first base is dropped where Ensembl VEP drops it. `None`
/// when the allele is not a sequence, when no alternate names it, or when two alternates of
/// different shapes name it.
pub fn allele_is_indel<S: AsRef<str>>(
    reference: &str,
    alternates: &[S],
    allele: &str,
    chopped: bool,
) -> Option<bool> {
    if allele.is_empty() || !allele.bytes().all(|b| b"ACGTN-".contains(&b)) {
        return None;
    }
    let mut plain = alternates
        .iter()
        .map(AsRef::as_ref)
        .filter(|a| !a.contains('*'));
    if let (Some(only), None) = (plain.next(), plain.next()) {
        return Some(stripped_length(reference) != stripped_length(only));
    }
    let trim = drops_first_base(reference, alternates, chopped);
    let kept = if trim {
        without_first(reference)
    } else {
        reference
    };
    let mut shapes: (bool, bool) = (false, false);
    for alternate in alternates.iter().map(AsRef::as_ref) {
        let alt = if trim && !alternate.contains('*') {
            without_first(alternate)
        } else {
            alternate
        };
        if dash(alt) == allele {
            if stripped_length(kept) != stripped_length(alt) {
                shapes.0 = true;
            } else {
                shapes.1 = true;
            }
        }
    }
    match shapes {
        (true, false) => Some(true),
        (false, true) => Some(false),
        _ => None,
    }
}

/// Whether the reference allele of a record is empty once the record's shared first base is
/// dropped, the alleles read as `allele_is_indel` reads them. Every sequence allele of such
/// a record is an insertion.
pub fn reference_allele_is_empty<S: AsRef<str>>(
    reference: &str,
    alternates: &[S],
    chopped: bool,
) -> bool {
    let kept = if drops_first_base(reference, alternates, chopped) {
        without_first(reference)
    } else {
        reference
    };
    kept.trim_matches('-').is_empty()
}

/// Whether `spelt` is `head`, `written` and `tail` one after the other.
fn spells(spelt: &[char], head: &[char], written: &[char], tail: &[char]) -> bool {
    spelt.len() == head.len() + written.len() + tail.len()
        && spelt.starts_with(head)
        && spelt[head.len()..].starts_with(written)
        && spelt.ends_with(tail)
}

/// The reference bases an output row's allele stands in for, read by position from a
/// record's own alleles. The record whose reference allele `reference` starts at `position`
/// writes `allele` at `location` where one of its alternates is the reference allele with
/// the Location's bases replaced by the allele (`-` replaces them by nothing): those bases
/// are returned. Where an alternate is the reference allele with the allele inserted between
/// the Location's two neighbouring positions the row is a pure insertion's and the empty
/// string is returned. `None` where the record writes no such row, where two alternates read
/// it differently, and for an allele or a Location that is not a sequence's. The record's
/// alleles are read in upper case. A Location position `integer_of` does not read is
/// refused.
pub fn replaced_reference<S: AsRef<str>>(
    position: i128,
    reference: &str,
    alternates: &[S],
    location: Option<&str>,
    allele: &str,
) -> Result<Option<String>, Refusal> {
    let Some(location) = location else {
        return Ok(None);
    };
    let sequence = !allele.is_empty() && allele.bytes().all(|b| b"ACGTN".contains(&b));
    if !sequence && allele != "-" {
        return Ok(None);
    }
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let (start, _, end) = pytext::partition(span, "-");
    let end = if end.is_empty() { start } else { end };
    if !(pytext::is_ascii_digits(start) && pytext::is_ascii_digits(end)) {
        return Ok(None);
    }
    let first = integer_of(start, "Location position")? - position;
    let last = integer_of(end, "Location position")? - position;
    let written: Vec<char> = if sequence {
        allele.chars().collect()
    } else {
        Vec::new()
    };
    let reference: Vec<char> = reference.to_uppercase().chars().collect();
    let length = reference.len() as i128;
    let mut found: Vec<String> = Vec::new();
    let mut find = |bases: String| {
        if !found.contains(&bases) {
            found.push(bases);
        }
    };
    for alternate in alternates {
        let spelt: Vec<char> = alternate.as_ref().to_uppercase().chars().collect();
        if 0 <= first && first <= last && last < length {
            let (first, last) = (first as usize, last as usize);
            if spells(
                &spelt,
                &reference[..first],
                &written,
                &reference[last + 1..],
            ) {
                find(reference[first..=last].iter().collect());
            }
        }
        if !written.is_empty() && 0 <= first && first < length && last == first + 1 {
            let after = first as usize + 1;
            if spells(&spelt, &reference[..after], &written, &reference[after..]) {
                find(String::new());
            }
        }
    }
    Ok(if found.len() == 1 { found.pop() } else { None })
}
