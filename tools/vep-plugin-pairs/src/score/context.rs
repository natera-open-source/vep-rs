// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows behind an unpaired pair of values, as the entries the class rules read, and the
//! classes that accept the pair.

use super::api::{Extracted, Keyword, NamedLines, PluginReads};
use super::class_table::{Classes, Reading};
use super::error::{Scored, ScorerError};
use super::merge::{Cursor, Lines, Span};
use super::sources::{text, without_chr, SourceOf, Sources, RESOLVED_LIMIT};
use crate::input::Kind;
use crate::rules::entry::{replaced_reference, DefaultHeader};
use crate::rules::record::{perl_split, vep_rs_pos};
use crate::rules::skipped_record::dropped_by_the_reference;
use crate::rules::{
    blank, decode_extra, order_only_fields, same_items, uploaded_alleles, Consult, Entry, Fields,
    Reads, RecordRows, Row, Rules, ABSENT, ORDER_ONLY, PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD,
};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::rc::Rc;

/// One unpaired pair of a field under a key: a value of each side, `None` for a side with
/// none to pair, and for each side the count its row is read by.
pub type Pair<'v> = (Option<&'v [u8]>, Option<&'v [u8]>, usize, usize);

/// A value with the number of its copies before it.
type Item<'v> = (&'v [u8], usize);

/// What the rows of an input record are asked by: a row's Feature, its record's name and
/// its allele.
type Asked = (Vec<u8>, Vec<u8>, Vec<u8>);

/// The tab-separated column `at` of a line.
fn column(line: &[u8], at: usize) -> Scored<&[u8]> {
    line.split(|&byte| byte == b'\t')
        .nth(at)
        .ok_or_else(|| short_row(line))
}

/// What follows the first `tabs` tabs of a line.
fn after(line: &[u8], tabs: usize) -> Scored<&[u8]> {
    memchr::memchr_iter(b'\t', line)
        .nth(tabs - 1)
        .map(|tab| &line[tab + 1..])
        .ok_or_else(|| short_row(line))
}

fn short_row(line: &[u8]) -> ScorerError {
    let key: Vec<&[u8]> = line.split(|&byte| byte == b'\t').take(3).collect();
    ScorerError(format!(
        "a line of the sorted rows under {} lacks a column",
        text(&key.join(&b' '))
    ))
}

/// The Location, Allele and Feature of a key, a tab after each.
fn key_columns(key: &[u8]) -> Scored<[String; 3]> {
    let body = key.strip_suffix(b"\t").unwrap_or(key);
    let columns: Vec<String> = body.split(|&byte| byte == b'\t').map(text).collect();
    <[String; 3]>::try_from(columns).map_err(|_| {
        ScorerError(format!(
            "the key {} is not a Location, an Allele and a Feature",
            text(body)
        ))
    })
}

/// A default-format input's rows at one Location: the row lines under each key and on each
/// Feature, their entries, parsed as they are asked for, and the rows of each input
/// record, found once per row's Feature, identifier and allele.
struct LocationRows {
    header: DefaultHeader,
    lines: Lines,
    /// The keys in ascending order, each as the row its rows start at, the length of the
    /// key in that row's line and its rows, which follow one another.
    keys: Vec<(usize, usize, Rc<[usize]>)>,
    /// The Features in the order their first row comes, each with its rows.
    by_feature: Vec<(Vec<u8>, Vec<usize>)>,
    records: Option<RecordRows>,
    of_record: HashMap<Asked, Option<Rc<Vec<Row>>>>,
    parsed: Vec<Option<Row>>,
}

impl LocationRows {
    fn new(header: &str) -> LocationRows {
        LocationRows {
            header: DefaultHeader::new(header),
            lines: Lines::default(),
            keys: Vec::new(),
            by_feature: Vec::new(),
            records: None,
            of_record: HashMap::new(),
            parsed: Vec::new(),
        }
    }

    /// Replaces the rows held with the lines of `cursor` that start with `location`, the
    /// input's rows at one Location.
    fn read(&mut self, cursor: &mut Cursor, location: &[u8]) -> Scored<()> {
        self.lines.clear();
        cursor.starting_into(location, &mut self.lines)?;
        self.index()
    }

    /// Files the lines held under their keys and Features. The lines come sorted, so the
    /// rows under one key follow one another and the keys ascend.
    fn index(&mut self) -> Scored<()> {
        self.keys.clear();
        self.by_feature.clear();
        self.records = None;
        self.of_record.clear();
        self.parsed.clear();
        self.parsed.resize(self.lines.len(), None);
        let mut features: HashMap<&[u8], usize> = HashMap::with_capacity(self.lines.len());
        // The first row and the key's length of the key in hand.
        let mut open: Option<(usize, usize)> = None;
        for (at, line) in self.lines.iter().enumerate() {
            let rest = after(line, 3)?;
            let key = &line[..line.len() - rest.len()];
            match open {
                Some((first, length)) if &self.lines.get(first)[..length] == key => {}
                _ => {
                    if let Some((first, length)) = open {
                        self.keys.push((first, length, (first..at).collect()));
                    }
                    open = Some((at, key.len()));
                }
            }
            let feature = column(line, 2)?;
            let place = *features.entry(feature).or_insert_with(|| {
                self.by_feature.push((feature.to_vec(), Vec::new()));
                self.by_feature.len() - 1
            });
            self.by_feature[place].1.push(at);
        }
        if let Some((first, length)) = open {
            let rows = (first..self.lines.len()).collect();
            self.keys.push((first, length, rows));
        }
        Ok(())
    }

    fn under(&self, key: &[u8]) -> Rc<[usize]> {
        let held = self
            .keys
            .binary_search_by(|&(first, length, _)| self.lines.get(first)[..length].cmp(key));
        match held {
            Ok(at) => Rc::clone(&self.keys[at].2),
            Err(_) => Rc::from(Vec::new()),
        }
    }

    fn on_feature(&self, feature: &[u8]) -> &[usize] {
        self.by_feature
            .iter()
            .find(|(held, _)| held == feature)
            .map_or(&[], |(_, rows)| rows.as_slice())
    }

    /// Every row held, Feature by Feature in the order each Feature's first row comes.
    fn held(&self) -> Vec<usize> {
        self.by_feature
            .iter()
            .flat_map(|(_, rows)| rows.iter().copied())
            .collect()
    }

    /// The row as its input wrote it.
    fn written(&self, at: usize) -> Scored<String> {
        Ok(text(after(self.lines.get(at), 4)?))
    }

    /// The `Uploaded_variation` of the row, the name of its input record.
    fn record(&self, at: usize) -> Scored<&[u8]> {
        column(self.lines.get(at), 4)
    }

    fn entry(&mut self, at: usize) -> Scored<Row> {
        if let Some(entry) = &self.parsed[at] {
            return Ok(Rc::clone(entry));
        }
        let row = String::from_utf8_lossy(after(self.lines.get(at), 4)?);
        let entry = self
            .header
            .row(&row)?
            .ok_or_else(|| ScorerError(format!("no row is read from the line: {row}")))?;
        self.parsed[at] = Some(Rc::clone(&entry));
        Ok(entry)
    }

    fn entries(&mut self, rows: &[usize]) -> Scored<Vec<Row>> {
        rows.iter().map(|&at| self.entry(at)).collect()
    }

    /// The rows of the input record of the row `at` among the rows held. Under an
    /// `Uploaded_variation` that spells the record's alleles, those on the row's Feature.
    /// Under an identifier, what `RecordRows` returns for the row, given the Location's rows
    /// and the input's lines that carry their identifiers. `None` where neither names one
    /// record, and where an allele repeats among the rows on the Feature, which rows of more
    /// than one record do.
    fn record_rows(
        &mut self,
        at: usize,
        input: Option<&dyn NamedLines>,
    ) -> Scored<Option<Rc<Vec<Row>>>> {
        let line = self.lines.get(at);
        // The record's name is a column of the row with another after it.
        after(line, 5)?;
        let (allele, feature, uploaded) = (
            column(line, 1)?.to_vec(),
            column(line, 2)?.to_vec(),
            column(line, 4)?.to_vec(),
        );
        let asked = (feature, uploaded, allele);
        if let Some(found) = self.of_record.get(&asked) {
            return Ok(found.clone());
        }
        let found = self.rows_of_record(at, &asked.0, &asked.1, input)?;
        self.of_record.insert(asked, found.clone());
        Ok(found)
    }

    fn rows_of_record(
        &mut self,
        at: usize,
        feature: &[u8],
        uploaded: &[u8],
        input: Option<&dyn NamedLines>,
    ) -> Scored<Option<Rc<Vec<Row>>>> {
        let mut group = Vec::new();
        let mut alleles = BTreeSet::new();
        for &other in self.on_feature(feature) {
            if self.record(other)? == uploaded {
                group.push(other);
                alleles.insert(column(self.lines.get(other), 1)?);
            }
        }
        if alleles.len() != group.len() {
            return Ok(None);
        }
        if uploaded_alleles(&text(uploaded)).is_some() {
            return Ok(Some(Rc::new(self.entries(&group)?)));
        }
        let Some(input) = input else {
            return Ok(None);
        };
        if self.records.is_none() {
            let held = self.held();
            let mut names = BTreeSet::new();
            for &other in &held {
                names.insert(self.record(other)?);
            }
            let carried: Vec<String> = names
                .into_iter()
                .flat_map(|name| input.named(name))
                .collect();
            let entries = self.entries(&held)?;
            self.records = Some(RecordRows::new(&entries, &carried));
        }
        let entry = self.entry(at)?;
        let found = match &self.records {
            Some(records) => records.of(&entry)?,
            None => None,
        };
        Ok(found.map(Rc::new))
    }
}

/// The rows of one unpaired pair: the reference's row, the engine's, the reference's rows of
/// the same record, and where each row stands among its side's rows at the Location, the
/// engine's only for a default-format engine.
struct PairRows {
    reference: Row,
    engine: Row,
    record_rows: Option<Rc<Vec<Row>>>,
    engine_at: Option<usize>,
    reference_at: usize,
}

/// The rows behind an unpaired pair of values and the classes that accept the pair. Keys
/// come in ascending byte order.
pub struct RowContext<'a> {
    reference: &'a Extracted,
    engine: &'a Extracted,
    /// The spans of the two sides' sorted rows that hold the rows of the keys asked for.
    ref_span: Span,
    eng_span: Span,
    sources: &'a Sources<'a>,
    classes: &'a Classes,
    rules: &'a dyn Rules,
    reads: &'a Reads,
    input: Option<&'a dyn NamedLines>,
    handed: HashMap<String, BTreeSet<Keyword>>,
    ref_cursor: Option<Cursor>,
    eng_cursor: Option<Cursor>,
    location: Option<Vec<u8>>,
    ref_rows: LocationRows,
    eng_rows: LocationRows,
    key: Option<Vec<u8>>,
    /// A default-format engine's rows under the key in hand.
    eng_at: Rc<[usize]>,
    eng_entries: Vec<Row>,
    /// A tuple file's tuple lines under the key in hand.
    eng_tuples: Vec<Vec<u8>>,
    order_only: BTreeMap<String, String>,
    dropped: Option<(Vec<u8>, bool)>,
    no_class: BTreeSet<String>,
    fields: HashMap<Vec<u8>, Rc<(String, Reading)>>,
}

impl<'a> RowContext<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference: &'a Extracted,
        engine: &'a Extracted,
        sources: &'a Sources<'a>,
        classes: &'a Classes,
        rules: &'a dyn Rules,
        reads: &'a Reads,
        plugin_reads: &PluginReads,
        input: Option<&'a dyn NamedLines>,
    ) -> RowContext<'a> {
        let mut handed: HashMap<String, BTreeSet<Keyword>> = HashMap::new();
        for (plugin, keywords) in &plugin_reads.0 {
            for class in classes.of_plugin(plugin) {
                handed.entry(class).or_default().extend(keywords);
            }
        }
        RowContext {
            reference,
            engine,
            ref_span: Span::WHOLE,
            eng_span: Span::WHOLE,
            sources,
            classes,
            rules,
            reads,
            input,
            handed,
            ref_cursor: None,
            eng_cursor: None,
            location: None,
            ref_rows: LocationRows::new(&reference.header),
            eng_rows: LocationRows::new(&engine.header),
            key: None,
            eng_at: Rc::from(Vec::new()),
            eng_entries: Vec::new(),
            eng_tuples: Vec::new(),
            order_only: order_only_fields("default"),
            dropped: None,
            no_class: BTreeSet::new(),
            fields: HashMap::new(),
        }
    }

    /// The same context for keys whose rows stand in `reference` of the reference's sorted
    /// rows and in `engine` of the engine's.
    pub fn within(mut self, reference: Span, engine: Span) -> RowContext<'a> {
        (self.ref_span, self.eng_span) = (reference, engine);
        self
    }

    fn default_engine(&self) -> bool {
        self.engine.kind == Kind::Default
    }

    /// The field's name as text and the classes of its source that read it, the applied
    /// ones first, looked up once per field name.
    fn field(&mut self, name: &[u8]) -> Scored<Rc<(String, Reading)>> {
        if let Some(held) = self.fields.get(name) {
            return Ok(Rc::clone(held));
        }
        let field = text(name);
        let reading = match self.sources.of(name)? {
            SourceOf::Named(source) => self.classes.reading(&source, &field),
            SourceOf::Core | SourceOf::Unclaimed => Reading::from(Vec::new()),
        };
        if self.fields.len() >= RESOLVED_LIMIT {
            self.fields.clear();
        }
        let held = Rc::new((field, reading));
        self.fields.insert(name.to_vec(), Rc::clone(&held));
        Ok(held)
    }

    /// Reads the reference's rows and a default-format engine's rows of the key's Location,
    /// once per Location, and a tuple file's tuples under the key.
    fn load(&mut self, key: &[u8]) -> Scored<()> {
        if self.ref_cursor.is_none() || self.eng_cursor.is_none() {
            self.ref_cursor = Some(Cursor::part(&self.reference.rows, self.ref_span)?);
            self.eng_cursor = Some(Cursor::part(&self.engine.rows, self.eng_span)?);
        }
        let (Some(ref_cursor), Some(eng_cursor)) = (&mut self.ref_cursor, &mut self.eng_cursor)
        else {
            return Ok(());
        };
        let default = self.engine.kind == Kind::Default;
        let location = match memchr::memchr(b'\t', key) {
            Some(tab) => &key[..=tab],
            None => key,
        };
        if self.location.as_deref() != Some(location) {
            self.location = Some(location.to_vec());
            self.ref_rows.read(ref_cursor, location)?;
            if default {
                self.eng_rows.read(eng_cursor, location)?;
            }
        }
        self.key = Some(key.to_vec());
        if default {
            self.eng_at = self.eng_rows.under(key);
            let rows = Rc::clone(&self.eng_at);
            self.eng_entries = self.eng_rows.entries(&rows)?;
        } else {
            self.eng_tuples = eng_cursor.starting(key)?;
        }
        Ok(())
    }

    fn engine_lines(&self) -> usize {
        if self.default_engine() {
            self.eng_at.len()
        } else {
            self.eng_tuples.len()
        }
    }

    /// The reference's rows under `key`, the rows of the key loaded. A side that writes a
    /// value of `field` under the key takes it from a row there, so a side with a value and
    /// no row under the key is a refusal: the reader of its sorted rows misses the key's
    /// rows, and its values would be left without a class.
    fn under(
        &mut self,
        key: &[u8],
        field: &str,
        ref_wrote: bool,
        eng_wrote: bool,
    ) -> Scored<Rc<[usize]>> {
        if self.key.as_deref() != Some(key) {
            self.load(key)?;
        }
        let rows = self.ref_rows.under(key);
        for (wrote, found, side) in [
            (ref_wrote, rows.len(), self.reference),
            (eng_wrote, self.engine_lines(), self.engine),
        ] {
            if wrote && found == 0 {
                let [location, allele, feature] = key_columns(key)?;
                return Err(ScorerError(format!(
                    "{}: no row under Location {location}, Allele {allele}, Feature {feature} was \
                     read from the sorted rows, where the input writes a value of {field}",
                    side.path.display()
                )));
            }
        }
        Ok(rows)
    }

    /// The index of the row that writes `field` as `value`, or for an absent value writes
    /// none, after `nth` rows that do the same; `None` without such a row.
    fn writing(entries: &[Row], field: &str, value: &str, nth: usize) -> Option<usize> {
        let mut nth = nth;
        for (index, entry) in entries.iter().enumerate() {
            let written = entry.fields.get(field).unwrap_or(ABSENT);
            let writes = if value == ABSENT {
                blank(written)
            } else {
                decode_extra(written) == value
            };
            if writes {
                if nth == 0 {
                    return Some(index);
                }
                nth -= 1;
            }
        }
        None
    }

    /// The tuples a tuple-file engine wrote under the key, as an entry; `None` where it
    /// wrote the field and `value` is absent, so that no row of it wrote none.
    fn tuple_entry(&self, field: &str, value: &str) -> Scored<Option<Row>> {
        let key = self.key.as_deref().unwrap_or_default();
        let [location, allele, feature] = key_columns(key)?;
        let mut fields = Fields::new();
        fields.insert("Location", location.as_str());
        fields.insert("Allele", allele.as_str());
        fields.insert("Feature", feature.as_str());
        for line in &self.eng_tuples {
            let tuple = line.get(key.len()..).unwrap_or_default();
            let (name, written) = match memchr::memchr(b'\t', tuple) {
                Some(tab) => (&tuple[..tab], &tuple[tab + 1..]),
                None => (tuple, &[][..]),
            };
            fields.set_default(&text(name), text(written));
        }
        if value != ABSENT {
            fields.insert(field, value);
        } else if fields.contains(field) {
            return Ok(None);
        }
        Ok(Some(Rc::new(Entry {
            key: format!("{location}|{allele}|{feature}|"),
            record: None,
            location: Some(location),
            allele,
            feature,
            consequence: String::new(),
            fields,
            ..Entry::default()
        })))
    }

    /// The rows of the pair; `None` where a side has no row under `key` that wrote its
    /// value of the pair, or wrote none where its value is absent, after the `ref_nth` and
    /// `eng_nth` rows that did and that belong to other pairs.
    fn rows(
        &mut self,
        key: &[u8],
        field: &str,
        expected: &str,
        actual: &str,
        ref_nth: usize,
        eng_nth: usize,
    ) -> Scored<Option<PairRows>> {
        let rows = self.under(key, field, expected != ABSENT, actual != ABSENT)?;
        let entries = self.ref_rows.entries(&rows)?;
        let Some(index) = Self::writing(&entries, field, expected, ref_nth) else {
            return Ok(None);
        };
        let mut engine_at = None;
        let engine = if self.default_engine() {
            let at = Self::writing(&self.eng_entries, field, actual, eng_nth);
            engine_at = at.map(|at| self.eng_at[at]);
            at.map(|at| Rc::clone(&self.eng_entries[at]))
        } else {
            self.tuple_entry(field, actual)?
        };
        let Some(engine) = engine else {
            return Ok(None);
        };
        Ok(Some(PairRows {
            reference: Rc::clone(&entries[index]),
            engine,
            record_rows: self.ref_rows.record_rows(rows[index], self.input)?,
            engine_at,
            reference_at: rows[index],
        }))
    }

    /// Gives a row under an identifier, as its entry, the reference bases its allele stands
    /// in for: read by position from the lines of the run's input that carry the identifier
    /// on the row's chromosome, a `chr` prefix set aside. A row whose `Uploaded_variation`
    /// spells the alleles carries its own from the parser. A row no such line states, or
    /// two of them read differently, keeps none.
    fn read_record(&self, entry: &Row, uploaded: Option<&[u8]>) -> Scored<()> {
        let (Some(uploaded), Some(input), Some(location)) =
            (uploaded, self.input, entry.location.as_deref())
        else {
            return Ok(());
        };
        if entry.reference.borrow().is_some() || uploaded_alleles(&text(uploaded)).is_some() {
            return Ok(());
        }
        let chromosome = without_chr(location.split(':').next().unwrap_or(location));
        let mut read: Vec<String> = Vec::new();
        for carrier in input.named(uploaded) {
            let columns: Vec<&str> = carrier.split('\t').collect();
            let &[carried_on, position, _, reference, alternates] = columns.as_slice() else {
                return Err(ScorerError(format!(
                    "the input line of {} is not five columns",
                    text(uploaded)
                )));
            };
            let Some(first) = vep_rs_pos(position) else {
                continue;
            };
            if without_chr(carried_on) != chromosome {
                continue;
            }
            let alternates = perl_split(alternates, ",");
            let found =
                replaced_reference(first, reference, &alternates, Some(location), &entry.allele)?;
            if let Some(found) = found {
                if !read.contains(&found) {
                    read.push(found);
                }
            }
        }
        if let [one] = read.as_slice() {
            *entry.reference.borrow_mut() = Some(one.clone());
        }
        Ok(())
    }

    /// Whether `key`'s Location holds the rows of a record the reference dropped, asked
    /// once per Location. A tuple file states no rows and its Locations are none; the record
    /// lines handed are those of the run's input that carry the first column of an engine
    /// row at the Location.
    fn dropped_record(&mut self, key: &[u8], field: &str) -> Scored<bool> {
        if !self.default_engine() {
            return Ok(false);
        }
        self.under(key, field, false, true)?;
        let location = match memchr::memchr(b'\t', key) {
            Some(tab) => &key[..=tab],
            None => key,
        };
        if let Some((held, dropped)) = &self.dropped {
            if held == location {
                return Ok(*dropped);
            }
        }
        let mut sides = Vec::new();
        for rows in [&self.ref_rows, &self.eng_rows] {
            let written: Scored<Vec<String>> =
                rows.held().into_iter().map(|at| rows.written(at)).collect();
            sides.push(written?);
        }
        let mut names = BTreeSet::new();
        for at in self.eng_rows.held() {
            names.insert(self.eng_rows.record(at)?);
        }
        let records: Vec<String> = match self.input {
            Some(input) => names
                .into_iter()
                .flat_map(|name| input.named(name))
                .collect(),
            None => Vec::new(),
        };
        let dropped =
            dropped_by_the_reference(&sides[0], &sides[1], &self.reference.path, &records)?;
        self.dropped = Some((location.to_vec(), dropped));
        Ok(dropped)
    }

    /// `value` as a class that reads a transcript is handed it: the empty text for a side
    /// that left no tuple of the field because its row writes the field with an empty value,
    /// which those classes read as a value; any other value as it is.
    fn as_written<'v>(entry: &Entry, field: &str, value: &'v str) -> &'v str {
        if value == ABSENT && entry.fields.get(field) == Some("") {
            ""
        } else {
            value
        }
    }

    /// Each value with the number of its copies before it; the copies are adjacent.
    fn numbered<'v>(values: &[&'v [u8]]) -> Vec<Item<'v>> {
        let mut numbered: Vec<Item<'v>> = Vec::with_capacity(values.len());
        for &value in values {
            let nth = match numbered.last() {
                Some(&(last, nth)) if last == value => nth + 1,
                _ => 0,
            };
            numbered.push((value, nth));
        }
        numbered
    }

    /// The values the two sides have left of a field under `key` as the pairs `judge`
    /// takes. Where no class reads the field, the engine is a tuple file or no rows share
    /// the key, the values pair in their order. Under a key rows share they pair by input
    /// record first, the `Uploaded_variation` of the row each value reads: a value with the
    /// other side's value of the same record, else with that record's row of the other side
    /// that wrote no value of the field; the values no record places pair in their order, a
    /// side with none to pair reading its rows that wrote no value and are not taken.
    pub fn pairs<'v>(
        &mut self,
        key: &[u8],
        name: &[u8],
        ref_left: &[&'v [u8]],
        eng_left: &[&'v [u8]],
    ) -> Scored<Vec<Pair<'v>>> {
        let (ref_items, eng_items) = (Self::numbered(ref_left), Self::numbered(eng_left));
        if self.default_engine() {
            let held = self.field(name)?;
            let (field, reading) = &*held;
            if !reading.is_empty() {
                let rows = self.under(key, field, !ref_items.is_empty(), !eng_items.is_empty())?;
                if rows.len() > 1 || self.eng_at.len() > 1 {
                    return self.by_record(field, &ref_items, &eng_items, &rows);
                }
            }
        }
        let longest = ref_items.len().max(eng_items.len());
        let mut in_order = Vec::with_capacity(longest);
        for at in 0..longest {
            let (ref_value, ref_nth) = match ref_items.get(at) {
                Some(&(value, nth)) => (Some(value), nth),
                None => (None, at - ref_items.len()),
            };
            let (eng_value, eng_nth) = match eng_items.get(at) {
                Some(&(value, nth)) => (Some(value), nth),
                None => (None, at - eng_items.len()),
            };
            in_order.push((ref_value, eng_value, ref_nth, eng_nth));
        }
        Ok(in_order)
    }

    /// Of one side's rows under a key, as entries and as the records the rows name: the
    /// record of the row each value of `items` reads, `None` without such a row, and the
    /// records of the rows that wrote no value of the field, in the rows' order.
    fn placed<'r>(
        field: &str,
        items: &[Item],
        entries: &[Row],
        records: &'r [Vec<u8>],
    ) -> (Vec<Option<&'r [u8]>>, Vec<&'r [u8]>) {
        let reads = items
            .iter()
            .map(|&(value, nth)| {
                Self::writing(entries, field, &text(value), nth).map(|at| records[at].as_slice())
            })
            .collect();
        let none = entries
            .iter()
            .zip(records)
            .filter(|(entry, _)| blank(entry.fields.get(field).unwrap_or(ABSENT)))
            .map(|(_, record)| record.as_slice())
            .collect();
        (reads, none)
    }

    /// Stands each value of `free` against its record's first row of the other side that
    /// wrote no value and is not taken, `none` being those rows' records: the values still
    /// free, and per value stood the row's count among those rows, in the values' order.
    fn stood(
        free: &[usize],
        records: &[Option<&[u8]>],
        none: &[&[u8]],
    ) -> (Vec<usize>, Vec<(usize, usize)>) {
        let mut still = Vec::new();
        let mut stood: Vec<(usize, usize)> = Vec::new();
        for &at in free {
            let row = records[at].and_then(|record| {
                (0..none.len()).find(|&row| {
                    none[row] == record && !stood.iter().any(|&(_, taken)| taken == row)
                })
            });
            match row {
                Some(row) => stood.push((at, row)),
                None => still.push(at),
            }
        }
        (still, stood)
    }

    /// The pairs of `pairs` under a key rows share; `ref_at` are the reference's rows under
    /// the key.
    fn by_record<'v>(
        &mut self,
        field: &str,
        ref_items: &[Item<'v>],
        eng_items: &[Item<'v>],
        ref_at: &[usize],
    ) -> Scored<Vec<Pair<'v>>> {
        let ref_entries = self.ref_rows.entries(ref_at)?;
        let records = |rows: &LocationRows, at: &[usize]| -> Scored<Vec<Vec<u8>>> {
            at.iter()
                .map(|&row| rows.record(row).map(<[u8]>::to_vec))
                .collect()
        };
        let ref_names = records(&self.ref_rows, ref_at)?;
        let eng_names = records(&self.eng_rows, &self.eng_at)?;
        let (ref_records, ref_none) = Self::placed(field, ref_items, &ref_entries, &ref_names);
        let (eng_records, eng_none) = Self::placed(field, eng_items, &self.eng_entries, &eng_names);
        let mut found: Vec<Pair<'v>> = Vec::new();
        let mut ref_free = Vec::new();
        let mut eng_free: Vec<usize> = (0..eng_items.len()).collect();
        for (at, &(value, nth)) in ref_items.iter().enumerate() {
            let same = ref_records[at].and_then(|record| {
                eng_free
                    .iter()
                    .position(|&other| eng_records[other] == Some(record))
            });
            match same {
                Some(place) => {
                    let (other, other_nth) = eng_items[eng_free.remove(place)];
                    found.push((Some(value), Some(other), nth, other_nth));
                }
                None => ref_free.push(at),
            }
        }
        let (ref_free, ref_stood) = Self::stood(&ref_free, &ref_records, &eng_none);
        let (eng_free, eng_stood) = Self::stood(&eng_free, &eng_records, &ref_none);
        for &(at, row) in &ref_stood {
            found.push((Some(ref_items[at].0), None, ref_items[at].1, row));
        }
        for &(at, row) in &eng_stood {
            found.push((None, Some(eng_items[at].0), row, eng_items[at].1));
        }
        // A side with none to pair reads its next row that wrote no value and is not taken;
        // past the last of them the count names no row, and the value is put to no class.
        let free_rows = |none: usize, taken: &[(usize, usize)], more: usize| -> Vec<usize> {
            (0..none)
                .filter(|row| !taken.iter().any(|(_, stood)| stood == row))
                .chain(none..none + more)
                .collect()
        };
        let ref_rows = free_rows(ref_none.len(), &eng_stood, eng_free.len());
        let eng_rows = free_rows(eng_none.len(), &ref_stood, ref_free.len());
        for at in 0..ref_free.len().max(eng_free.len()) {
            let (ref_value, ref_nth) = match ref_free.get(at) {
                Some(&free) => (Some(ref_items[free].0), ref_items[free].1),
                None => (None, ref_rows[at - ref_free.len()]),
            };
            let (eng_value, eng_nth) = match eng_free.get(at) {
                Some(&free) => (Some(eng_items[free].0), eng_items[free].1),
                None => (None, eng_rows[at - eng_free.len()]),
            };
            found.push((ref_value, eng_value, ref_nth, eng_nth));
        }
        Ok(found)
    }

    /// The class that accepts the field's two unpaired values under `key` (a side that
    /// wrote none is `None`): the first class of the field's source reading the field to
    /// accept, the applied ones consulted first. `ref_nth` and `eng_nth` count the other
    /// unpaired pairs of the key and field that carry the same value on that side, so that
    /// each pair reads rows of its own.
    pub fn judge(
        &mut self,
        key: &[u8],
        name: &[u8],
        ref_value: Option<&[u8]>,
        eng_value: Option<&[u8]>,
        ref_nth: usize,
        eng_nth: usize,
    ) -> Scored<Option<String>> {
        let held = self.field(name)?;
        let (field, reading) = &*held;
        if reading.is_empty() {
            return Ok(None);
        }
        let expected = ref_value.map_or_else(|| ABSENT.to_string(), text);
        let actual = eng_value.map_or_else(|| ABSENT.to_string(), text);
        let mut rows: Option<Option<PairRows>> = None;
        for class in reading.iter() {
            let accepted = if class == ORDER_ONLY {
                match self.order_only.get(field) {
                    Some(separator)
                        if !separator.is_empty() && ref_value.is_some() && eng_value.is_some() =>
                    {
                        same_items(&expected, &actual, separator)
                    }
                    _ => false,
                }
            } else if class == PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD {
                ref_value.is_none() && eng_value.is_some() && self.dropped_record(key, field)?
            } else {
                if rows.is_none() {
                    rows = Some(self.rows(key, field, &expected, &actual, ref_nth, eng_nth)?);
                }
                match rows.as_ref().and_then(Option::as_ref) {
                    Some(pair) => self.consult(class, field, &expected, &actual, pair)?,
                    None => false,
                }
            };
            if accepted {
                return Ok(Some(class.clone()));
            }
        }
        Ok(None)
    }

    /// Puts the pair to one class with the rows and the reads the class is handed.
    fn consult(
        &mut self,
        class: &str,
        field: &str,
        expected: &str,
        actual: &str,
        pair: &PairRows,
    ) -> Scored<bool> {
        let (mut expected, mut actual) = (expected, actual);
        if self.rules.transcript_reading(class) {
            expected = Self::as_written(&pair.reference, field, expected);
            actual = Self::as_written(&pair.engine, field, actual);
        }
        let candidate_rows = match pair.engine_at {
            Some(at) if self.rules.sibling_allele(class) => {
                self.eng_rows.record_rows(at, self.input)?
            }
            _ => None,
        };
        if self.rules.record_reference(class) {
            let named = self.ref_rows.record(pair.reference_at)?;
            self.read_record(&pair.reference, Some(named))?;
            let named = match pair.engine_at {
                Some(at) => Some(self.eng_rows.record(at)?),
                None => None,
            };
            self.read_record(&pair.engine, named)?;
        }
        let handed = self.handed.get(class);
        let hands = |keyword: Keyword| handed.is_some_and(|keywords| keywords.contains(&keyword));
        let reads = self.reads;
        let source = match self.sources.of(field.as_bytes())? {
            SourceOf::Named(source) => Some(source),
            SourceOf::Core | SourceOf::Unclaimed => None,
        };
        let consult = Consult {
            klass: class,
            e: &pair.reference,
            name: field,
            expected,
            actual,
            a: Some(&pair.engine),
            record_rows: pair.record_rows.as_ref().map(|rows| rows.as_slice()),
            candidate_rows: candidate_rows.as_ref().map(|rows| rows.as_slice()),
            declared: if hands(Keyword::Declared) {
                self.classes.applied()
            } else {
                &self.no_class
            },
            assembly: self.sources.assembly(),
            source: source.as_deref(),
            gwas_pairs: reads
                .gwas_pairs
                .as_ref()
                .filter(|_| hands(Keyword::GwasPairs)),
            revel_table: reads
                .revel_table
                .as_ref()
                .filter(|_| hands(Keyword::RevelTable)),
            utrannotator_sweep: reads
                .utrannotator_sweep
                .as_ref()
                .filter(|_| hands(Keyword::UtrannotatorSweep)),
            downstream_sequences: reads
                .downstream_sequences
                .as_deref()
                .filter(|_| hands(Keyword::DownstreamSequences)),
            splice_windows: reads
                .splice_windows
                .as_ref()
                .filter(|_| hands(Keyword::SpliceWindows)),
            spliceai_table: reads
                .spliceai_table
                .as_ref()
                .filter(|_| hands(Keyword::SpliceaiTable)),
            exon_spans: reads
                .exon_spans
                .as_ref()
                .filter(|_| hands(Keyword::ExonSpans)),
        };
        Ok(self.rules.accepts(&consult)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::parse_default;

    fn row(uploaded: &str, extra: &str) -> Row {
        let header =
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\
                      \tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\
                      \tExisting_variation\tExtra";
        let line =
            format!("{uploaded}\t1:100\tT\t-\tENST1\tTranscript\t-\t-\t-\t-\t-\t-\t-\t{extra}");
        parse_default(&format!("{header}\n{line}"))
            .unwrap()
            .remove(0)
    }

    fn items<'v>(values: &[&'v str]) -> Vec<Item<'v>> {
        let values: Vec<&[u8]> = values.iter().map(|value| value.as_bytes()).collect();
        RowContext::numbered(&values)
    }

    #[test]
    fn copies_of_a_value_are_numbered_in_their_order() {
        assert_eq!(
            items(&["a", "a", "b", "a"]),
            [
                (b"a".as_slice(), 0),
                (b"a".as_slice(), 1),
                (b"b".as_slice(), 0),
                (b"a".as_slice(), 0)
            ]
        );
        assert!(items(&[]).is_empty());
    }

    #[test]
    fn nth_row_that_writes_a_value_or_none_is_found_among_the_rows() {
        let rows = [
            row("rs1", "PLANT_score=1;IMPACT=LOW"),
            row("rs2", "IMPACT=LOW"),
            row("rs3", "PLANT_score=1"),
            row("rs4", "PLANT_score=-"),
            row("rs5", "PLANT_score=a%3Db"),
            row("rs6", "PLANT_score="),
        ];
        let writing =
            |value: &str, nth: usize| RowContext::writing(&rows, "PLANT_score", value, nth);
        assert_eq!(writing("1", 0), Some(0));
        assert_eq!(writing("1", 1), Some(2));
        assert_eq!(writing("1", 2), None);
        assert_eq!(writing("a=b", 0), Some(4));
        assert_eq!(writing("a%3Db", 0), None);
        assert_eq!(writing(ABSENT, 0), Some(1));
        assert_eq!(writing(ABSENT, 1), Some(3));
        assert_eq!(writing(ABSENT, 2), Some(5));
        assert_eq!(writing(ABSENT, 3), None);
        assert_eq!(writing("2", 0), None);
    }

    #[test]
    fn values_are_placed_on_the_records_of_the_rows_that_wrote_them() {
        let rows = [
            row("rs1", "PLANT_score=1"),
            row("rs2", "IMPACT=LOW"),
            row("rs3", "PLANT_score=1"),
            row("rs2", "PLANT_score=-"),
        ];
        let records: Vec<Vec<u8>> = ["rs1", "rs2", "rs3", "rs2"]
            .iter()
            .map(|name| name.as_bytes().to_vec())
            .collect();
        let (reads, none) = RowContext::placed(
            "PLANT_score",
            &items(&["1", "1", "1", "7"]),
            &rows,
            &records,
        );
        assert_eq!(
            reads,
            [Some(b"rs1".as_slice()), Some(b"rs3".as_slice()), None, None]
        );
        assert_eq!(none, [b"rs2".as_slice(), b"rs2".as_slice()]);
    }

    #[test]
    fn value_stands_against_its_records_first_free_row_that_wrote_none() {
        let records: Vec<Option<&[u8]>> =
            vec![Some(b"rs1"), Some(b"rs2"), None, Some(b"rs2"), Some(b"rs2")];
        let none: Vec<&[u8]> = vec![b"rs2", b"rs9", b"rs2"];
        let (still, stood) = RowContext::stood(&[0, 1, 2, 3, 4], &records, &none);
        assert_eq!(still, [0, 2, 4]);
        assert_eq!(stood, [(1, 0), (3, 2)]);
        let (still, stood) = RowContext::stood(&[4, 3], &records, &none);
        assert!(still.is_empty());
        assert_eq!(stood, [(4, 0), (3, 2)]);
    }

    #[test]
    fn empty_value_of_a_row_is_handed_to_a_transcript_reading_class_as_a_value() {
        let empty = row("rs1", "PLANT_score=");
        let lacking = row("rs1", "IMPACT=LOW");
        assert_eq!(RowContext::as_written(&empty, "PLANT_score", ABSENT), "");
        assert_eq!(
            RowContext::as_written(&lacking, "PLANT_score", ABSENT),
            ABSENT
        );
        assert_eq!(RowContext::as_written(&empty, "PLANT_score", "1"), "1");
    }

    #[test]
    fn row_line_is_read_by_its_columns() {
        let ordinal = format!("{:012}", 7);
        let line = format!("1:100\tT\tENST1\t{ordinal}\trs1\t1:100\tT");
        let line = line.as_bytes();
        assert_eq!(column(line, 1).unwrap(), b"T");
        assert_eq!(column(line, 4).unwrap(), b"rs1");
        assert_eq!(after(line, 4).unwrap(), b"rs1\t1:100\tT");
        assert_eq!(
            after(line, 3).unwrap(),
            format!("{ordinal}\trs1\t1:100\tT").as_bytes()
        );
        assert!(column(b"1:100\tT", 2).is_err());
        assert_eq!(
            after(b"1:100\tT\tENST1", 3).err(),
            Some(ScorerError(
                "a line of the sorted rows under 1:100 T ENST1 lacks a column".to_string()
            ))
        );
        assert_eq!(
            key_columns(b"1:100\tT\tENST1\t").unwrap(),
            ["1:100".to_string(), "T".to_string(), "ENST1".to_string()]
        );
        assert!(key_columns(b"1:100\tT\t").is_err());
    }

    #[test]
    fn rows_at_a_location_are_held_by_key_and_by_feature_in_their_order() {
        let mut rows = LocationRows::new("#Uploaded_variation\tLocation\tAllele\tExtra");
        for line in [
            "1:100\tA\tENST2\t000000000003\trs1\t1:100\tA\tPLANT_score=3",
            "1:100\tT\tENST1\t000000000001\trs1\t1:100\tT\tPLANT_score=1",
            "1:100\tT\tENST1\t000000000004\trs2\t1:100\tT\tPLANT_score=4",
            "1:100\tT\tENST2\t000000000002\trs1\t1:100\tT\tPLANT_score=2",
        ] {
            rows.lines.push(line.as_bytes());
        }
        rows.index().unwrap();
        assert_eq!(&*rows.under(b"1:100\tT\tENST1\t"), [1, 2]);
        assert!(rows.under(b"1:100\tT\tENST3\t").is_empty());
        assert_eq!(rows.on_feature(b"ENST2"), [0, 3]);
        assert_eq!(rows.held(), [0, 3, 1, 2]);
        assert_eq!(rows.record(2).unwrap(), b"rs2");
        assert_eq!(rows.written(1).unwrap(), "rs1\t1:100\tT\tPLANT_score=1");
        let entry = rows.entry(2).unwrap();
        assert_eq!(entry.fields.get("PLANT_score"), Some("4"));
        assert!(Rc::ptr_eq(&entry, &rows.entry(2).unwrap()));
        // Two rows of one record on a Feature with one allele each are the record's rows
        // only where the record's name spells its alleles or the input names its line.
        assert!(rows.record_rows(1, None).unwrap().is_none());
        rows.lines.clear();
        rows.index().unwrap();
        assert!(rows.held().is_empty() && rows.under(b"1:100\tT\tENST1\t").is_empty());
    }
}
