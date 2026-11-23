// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The SpliceAI table, its refusals and the multi-entry cell.

use crate::argparse::path_text;
use crate::cache;
use crate::classes::SPLICEAI_MULTI_ENTRY_CELL;
use crate::entry::Entry;
use crate::error::Error;
use crate::error::Result;
use crate::fields::as_number;
use crate::fields::field_named;
use crate::fields::python_type;
use crate::lists::json_object_member;
use crate::pyfmt;
use crate::pyjson;
use crate::pyjson::Value;
use crate::rules::downstream::field_of;
use crate::rules::record_alleles::vcf_alleles_as_spelt;
use crate::rules::revel::RevelTable;
use crate::shared::blank;
use indexmap::IndexMap;
use std::cell::RefCell;
use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;

pub const SPLICEAI_KEY: &str = "SpliceAI=";

pub const SPLICEAI_VALUE_FIELD: &str = "SpliceAI_pred";

pub const SPLICEAI_CUTOFF_FIELD: &str = "SpliceAI_cutoff";

pub const SPLICEAI_MEMBERS: &[&str] = &[
    "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
];

pub const SPLICEAI_CUTOFF_FLAGS: &[&str] = &["PASS", "FAIL"];

/// The grounds on which `spliceai_multi_entry_cell` refuses a row it is asked about under the
/// score files, in the order `spliceai_multi_entry_cell_ground` meets them: a key and what
/// the predicate read.
pub const SPLICEAI_GROUNDS: &[(&str, &str)] = &[
    (
        "no_candidate_value",
        "the candidate's row writes no SpliceAI value",
    ),
    ("no_span", "the row's Location names no span"),
    (
        "other_info",
        "a record line at the row's variant carries another INFO part beside the key or an entry of another length",
    ),
    (
        "several_references",
        "the record lines at the row's variant spell more than one REF",
    ),
    (
        "no_record",
        "no record line of the score files stands at the row's variant",
    ),
    (
        "no_gene",
        "neither a SYMBOL nor the candidate cache names the row's gene",
    ),
    (
        "gene_not_named",
        "no entry of a record line at the row's variant names the row's gene",
    ),
    (
        "lone_entry",
        "every record line that holds the gene's entry holds it alone",
    ),
    (
        "candidate_value",
        "the candidate's value is not the gene's entries of the last record line that names it",
    ),
    (
        "reference_not_silent",
        "the row's gene leads no record line and the reference's row writes a SpliceAI value or a cutoff flag, or the field is not the candidate's alone",
    ),
    (
        "reference_value",
        "the row's gene leads a record line and the reference's value is not the whole value of the last line it leads",
    ),
    (
        "reference_flag",
        "the row's gene leads a record line and no cutoff gives that line's first entry the reference's flag and the gene's entries the candidate's",
    ),
];

/// The Location a predicate over a default row's Location reads of the entry:
/// `at_its_location(e).location`, the entry's own, else the one its record gives its allele.
fn location_of(e: &Entry) -> Option<&str> {
    e.location.as_deref().or(e.allele_location.as_deref())
}

/// A transcript's `gene_symbol` as a candidate JSON cache holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneSymbol {
    Text(String),
    /// A value that is no text: it equals the symbol of no score entry.
    NoText,
}

/// The gene symbol of each transcript of a candidate JSON cache (`--json-cache`), by
/// stable id: the symbol SpliceAI.pm:314 looks a transcript's score record up under. A
/// chromosome's shards are read on the first lookup of a row whose Location names its
/// directory, every shard on the first lookup of any other entry.
#[derive(Debug)]
pub struct TranscriptSymbols {
    pub root: PathBuf,
    pub _symbols: RefCell<IndexMap<Option<String>, IndexMap<String, GeneSymbol>>>,
}

/// The stable id and gene symbol of a transcript that states both. A stable id that is no
/// text is the key of no Feature and is left out; one that is an array or an object is no
/// key at all.
fn named_gene(transcript: &Value) -> Option<(String, GeneSymbol)> {
    let stable_id = transcript
        .get("stable_id")
        .filter(|id| pyjson::truthy(id))?;
    let symbol = transcript
        .get("gene_symbol")
        .filter(|symbol| pyjson::truthy(symbol))?;
    let symbol = match symbol {
        Value::Str(text) => GeneSymbol::Text(text.clone()),
        _ => GeneSymbol::NoText,
    };
    match stable_id {
        Value::Str(id) => Some((id.clone(), symbol)),
        Value::List(_) => {
            crate::uncaught!("cannot use 'list' as a dict key (unhashable type: 'list')")
        }
        Value::Dict(_) => {
            crate::uncaught!("cannot use 'dict' as a dict key (unhashable type: 'dict')")
        }
        _ => None,
    }
}

impl TranscriptSymbols {
    pub fn new(cache: &Path) -> Result<Self> {
        let root = cache::joined(cache, "transcripts");
        if !root.is_dir() {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        Ok(TranscriptSymbols {
            root,
            _symbols: RefCell::new(IndexMap::new()),
        })
    }

    pub fn _read(&self, chromosome: Option<&str>) -> Result<IndexMap<String, GeneSymbol>> {
        let shards = cache::shards(&self.root, chromosome)?;
        let read = cache::map_shards(&shards, cache::threads(), |_, transcripts| {
            Ok(transcripts
                .iter()
                .filter_map(named_gene)
                .collect::<Vec<_>>())
        })?;
        let mut symbols = IndexMap::new();
        for (stable_id, symbol) in read.into_iter().flatten() {
            symbols.insert(stable_id, symbol);
        }
        Ok(symbols)
    }

    /// The gene symbol of the entry's transcript, or `None` for a transcript the cache
    /// does not hold or gives no symbol; a Feature carrying a version the cache's stable id
    /// lacks (`--transcript_version`) is looked up without it.
    pub fn symbol(&self, e: &Entry) -> Result<Option<GeneSymbol>> {
        let mut chromosome = location_of(e)
            .filter(|location| !location.is_empty())
            .map(|location| pyfmt::partition(location, ":").0.to_string());
        if chromosome
            .as_deref()
            .is_some_and(|name| !cache::joined(&self.root, name).is_dir())
        {
            chromosome = None;
        }
        if !self._symbols.borrow().contains_key(&chromosome) {
            let read = self._read(chromosome.as_deref())?;
            self._symbols.borrow_mut().insert(chromosome.clone(), read);
        }
        let symbols = self._symbols.borrow();
        let symbols = &symbols[&chromosome];
        let unversioned = pyfmt::partition(&e.feature, ".").0;
        Ok(symbols
            .get(&e.feature)
            .or_else(|| symbols.get(unversioned))
            .cloned())
    }
}

/// The score files a scored run's SpliceAI plugin read (`--spliceai-table`, once per
/// file: the `snv=` file and the `indel=` file), each a VCF with its tabix index beside
/// it (SpliceAI.pm:114-122), read through that index at the positions a row names and never
/// whole. `RevelTable`
/// reads the rows of a tabix-indexed table at a span of positions, and a VCF's eight columns
/// read as its eight-column layout: the position, REF and ALT cells stand in their own
/// places and the INFO column in the score's. `json_cache` is the run's candidate JSON
/// cache, whose transcripts name the gene of a row that states no `SYMBOL`
/// (`TranscriptSymbols`). `refusals` counts the rows `spliceai_multi_entry_cell` is
/// asked about under the files and those it refuses (`SpliceAiRefusals`).
#[derive(Debug)]
pub struct SpliceAiTable {
    pub paths: Vec<PathBuf>,
    pub symbols: Option<TranscriptSymbols>,
    pub refusals: SpliceAiRefusals,
    pub _files: Vec<RevelTable>,
}

impl SpliceAiTable {
    pub fn new(paths: &[&Path], json_cache: Option<&Path>) -> Result<Self> {
        let symbols = json_cache.map(TranscriptSymbols::new).transpose()?;
        let mut files = Vec::with_capacity(paths.len());
        for path in paths {
            let table = RevelTable::new(path, false).map_err(|exc| Self::_refusal(&exc))?;
            if table._columns.0 != 3 {
                return Err(Error::comparator(format!(
                    "--spliceai-table: {} is not a VCF (a header line of eight columns)",
                    path_text(path)
                )));
            }
            files.push(table);
        }
        Ok(SpliceAiTable {
            paths: paths.iter().map(|path| path.to_path_buf()).collect(),
            symbols,
            refusals: SpliceAiRefusals::new(),
            _files: files,
        })
    }

    /// A refusal of the table reader under this option's name.
    pub fn _refusal(exc: &Error) -> Error {
        let text = exc.to_string();
        let reason = pyfmt::partition(&text, ": ").2.replace(
            "is not a REVEL table (a header of 7, 8 or 9 columns)",
            "is not a VCF (a header line of eight columns)",
        );
        Error::comparator(format!("--spliceai-table: {reason}"))
    }

    /// The POS, REF, ALT and INFO cells of every record at the positions `low` to
    /// `high` of `chromosome`, the files in the order given and each in file order; a
    /// line that lacks one of the four cells is left out.
    pub fn records(
        &self,
        chromosome: &str,
        low: i64,
        high: i64,
    ) -> Result<Vec<(i64, String, String, String)>> {
        let mut found = Vec::new();
        for table in &self._files {
            let rows = table
                .rows(chromosome, low, high)
                .map_err(|exc| Self::_refusal(&exc))?;
            for row in rows {
                if !pyfmt::is_digit(&row.position) {
                    continue;
                }
                let filled = |cell: Option<String>| cell.filter(|cell| !cell.is_empty());
                let (Some(reference), Some(alternate), Some(score)) = (
                    filled(row.reference),
                    filled(row.alternate),
                    filled(row.score),
                ) else {
                    continue;
                };
                // A digit that is no decimal digit (a superscript) passes `isdigit` and
                // `int` refuses it.
                let position = match pyfmt::int_checked(&row.position) {
                    Ok(decimal) => decimal.parse::<i64>().unwrap_or(i64::MAX),
                    Err(error) => crate::uncaught!("{error}"),
                };
                found.push((position, reference, alternate, score));
            }
        }
        Ok(found)
    }
}

/// A row as `SpliceAiRefusals` counts it: its entry's key, Location, allele and feature.
pub type Row = (String, Option<String>, String, String);

/// The rows `spliceai_multi_entry_cell` is asked about under the score files and, of
/// them, the rows it refuses, each under the first ground met for it (`SPLICEAI_GROUNDS`).
/// A row is its entry's key, Location, allele and feature, and counts once however many of
/// its fields the class is asked about.
#[derive(Debug)]
pub struct SpliceAiRefusals {
    pub _asked: RefCell<BTreeSet<Row>>,
    pub _refused: RefCell<IndexMap<Row, String>>,
}

impl SpliceAiRefusals {
    pub fn new() -> Self {
        SpliceAiRefusals {
            _asked: RefCell::new(BTreeSet::new()),
            _refused: RefCell::new(IndexMap::new()),
        }
    }

    /// Counts the row of `e` as asked about and, under a `ground`, as refused.
    pub fn note(&self, e: &Entry, ground: Option<&str>) {
        let row = (
            e.key.clone(),
            e.location.clone(),
            e.allele.clone(),
            e.feature.clone(),
        );
        if let Some(ground) = ground {
            self._refused
                .borrow_mut()
                .entry(row.clone())
                .or_insert_with(|| ground.to_string());
        }
        self._asked.borrow_mut().insert(row);
    }

    /// The rows `asked` about, the rows `refused` and, under `grounds`, the refused
    /// rows of each ground that refused any, in the order of `SPLICEAI_GROUNDS`.
    pub fn counts(&self) -> Value {
        let refused = self._refused.borrow();
        let count = |rows: usize| Value::Int(i64::try_from(rows).unwrap_or(i64::MAX));
        let grounds: IndexMap<String, Value> = SPLICEAI_GROUNDS
            .iter()
            .map(|(ground, _)| {
                let rows = refused.values().filter(|met| met == ground).count();
                (*ground, rows)
            })
            .filter(|(_, rows)| *rows > 0)
            .map(|(ground, rows)| (ground.to_string(), count(rows)))
            .collect();
        Value::Dict(IndexMap::from([
            ("asked".to_string(), count(self._asked.borrow().len())),
            ("refused".to_string(), count(refused.len())),
            ("grounds".to_string(), Value::Dict(grounds)),
        ]))
    }

    pub fn clear(&self) {
        self._asked.borrow_mut().clear();
        self._refused.borrow_mut().clear();
    }
}

impl Default for SpliceAiRefusals {
    fn default() -> Self {
        Self::new()
    }
}

/// `counts[key]` of the counts a summary prints; counts that are no object, or hold no
/// such key, end the run as the interpreter's subscript does.
fn count_of<'a>(counts: &'a Value, key: &str) -> &'a Value {
    match counts {
        Value::Dict(members) => match members.get(key) {
            Some(value) => value,
            None => crate::uncaught!("{}", pyfmt::repr_str(key)),
        },
        Value::Str(_) => crate::uncaught!("string indices must be integers, not 'str'"),
        Value::List(_) => crate::uncaught!("list indices must be integers or slices, not str"),
        other => crate::uncaught!("'{}' object is not subscriptable", python_type(other)),
    }
}

/// The lines a summary prints of the counts of a `SpliceAiRefusals`: the rows the class
/// was asked about and refused, then the refused rows of each ground. `whose` names the
/// rows' side of a summary of several.
pub fn spliceai_refusal_lines(counts: &Value, whose: &str) -> Vec<String> {
    let klass = SPLICEAI_MULTI_ENTRY_CELL;
    let mut lines = vec![format!(
        "{klass} under --spliceai-table{whose}: rows asked about {}, rows refused {}",
        pyjson::py_str(count_of(counts, "asked")),
        pyjson::py_str(count_of(counts, "refused"))
    )];
    let grounds = count_of(counts, "grounds");
    let Some(grounds) = grounds.as_dict() else {
        crate::uncaught!("'{}' object has no attribute 'items'", python_type(grounds));
    };
    for (ground, rows) in grounds {
        let Some((_, read)) = SPLICEAI_GROUNDS.iter().find(|(key, _)| key == ground) else {
            crate::uncaught!("{}", pyfmt::repr_str(ground));
        };
        lines.push(format!(
            "{klass}{whose}: rows refused {}: {read}",
            pyjson::py_str(rows)
        ));
    }
    lines
}

/// A row's SpliceAI value: its `SpliceAI_pred` text where its format writes the one
/// field, else `None`; its nine members (`SPLICEAI_MEMBERS`); its cutoff flag, empty
/// for none; and whether the members are JSON values, whose numbers compare as numbers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpliceAiValue {
    pub whole: Option<String>,
    pub members: Vec<String>,
    pub cutoff: String,
    pub numbers: bool,
}

impl SpliceAiValue {
    /// Whether the row writes `value`, a `|`-joined value as the default format
    /// prints it: as its whole text where the row carries the one field, else as its nine
    /// members, the parts of `value` past the ninth being those SpliceAI.pm:276-287 drop.
    pub fn writes(&self, value: &str) -> bool {
        if let Some(whole) = &self.whole {
            return whole == value;
        }
        let parts: Vec<&str> = value.split('|').take(SPLICEAI_MEMBERS.len()).collect();
        parts.len() == self.members.len()
            && self.members.iter().zip(parts).all(|(member, part)| {
                member == part
                    || (self.numbers
                        && as_number(member).is_some_and(|number| as_number(part) == Some(number)))
            })
    }

    /// Whether the row writes `entries`, the entries of one score record that name one
    /// gene, in the record's order, each its ten members: as its whole text, the first
    /// entry behind its allele and each later one whole, joined by commas, where the row
    /// carries the one field; else as its nine members, each the entries' values of that
    /// member joined by commas. One entry is what `writes` reads of its value.
    pub fn writes_entries(&self, entries: &[&Vec<String>]) -> bool {
        if let Some(whole) = &self.whole {
            let spelt: Vec<String> = entries
                .iter()
                .enumerate()
                .map(|(place, entry)| entry[usize::from(place == 0)..].join("|"))
                .collect();
            return *whole == spelt.join(",");
        }
        self.members.len() == SPLICEAI_MEMBERS.len()
            && self.members.iter().enumerate().all(|(index, member)| {
                let part = entries
                    .iter()
                    .map(|entry| entry[1 + index].as_str())
                    .collect::<Vec<_>>()
                    .join(",");
                *member == part
                    || (self.numbers
                        && as_number(member).is_some_and(|number| as_number(&part) == Some(number)))
            })
    }
}

/// The largest of an entry's four delta scores, the members behind its allele and symbol
/// (SpliceAI.pm:336-338); a score that is no number counts as zero.
fn largest_delta_score(entry: &[String]) -> f64 {
    entry[2..6]
        .iter()
        .filter_map(|score| score.parse::<f64>().ok())
        .fold(0.0, f64::max)
}

/// Whether one cutoff gives two rows their differing flags: the reference's `theirs` by
/// the largest delta score of `first`, the first entry of the record line its value comes
/// from (SpliceAI.pm:336-338), and the candidate's `ours` by the largest over `own`, the
/// gene's entries of the candidate's line, a row passing where any of them reaches the
/// cutoff. A cutoff the one score fails and the other reaches lies between the two.
fn one_cutoff_gives_both_flags(
    theirs: &str,
    first: &[String],
    ours: &str,
    own: &[&Vec<String>],
) -> bool {
    let reference = largest_delta_score(first);
    let candidate = own
        .iter()
        .map(|entry| largest_delta_score(entry))
        .fold(0.0, f64::max);
    match (theirs, ours) {
        ("FAIL", "PASS") => reference < candidate,
        ("PASS", "FAIL") => candidate < reference,
        _ => false,
    }
}

/// The SpliceAI value of a row, read from its `SpliceAI_pred` field, else from its nine
/// `SpliceAI_pred_` fields (a CSQ field's `&` read as the comma it stands for), else from
/// the members of its JSON `spliceai` object; `None` for a row that writes none.
pub fn spliceai_row_value(e: Option<&Entry>) -> Option<SpliceAiValue> {
    let e = e?;
    let written = |name: &str| -> String {
        let value = field_of(Some(e), &[name]);
        if blank(&value) {
            String::new()
        } else {
            value
        }
    };
    let filled = |members: &[String]| members.iter().any(|member| !member.is_empty());
    let whole = written(SPLICEAI_VALUE_FIELD);
    if !whole.is_empty() {
        let members = whole
            .split('|')
            .take(SPLICEAI_MEMBERS.len())
            .map(str::to_string)
            .collect();
        return Some(SpliceAiValue {
            whole: Some(whole),
            members,
            cutoff: written(SPLICEAI_CUTOFF_FIELD),
            numbers: false,
        });
    }
    let split: Vec<String> = SPLICEAI_MEMBERS
        .iter()
        .map(|member| written(&format!("{SPLICEAI_VALUE_FIELD}_{member}")).replace('&', ","))
        .collect();
    if filled(&split) {
        return Some(SpliceAiValue {
            whole: None,
            members: split,
            cutoff: written(SPLICEAI_CUTOFF_FIELD),
            numbers: false,
        });
    }
    let members: Vec<String> = SPLICEAI_MEMBERS
        .iter()
        .map(|member| json_object_member(e, "spliceai", member).unwrap_or_default())
        .collect();
    if filled(&members) {
        return Some(SpliceAiValue {
            whole: None,
            members,
            cutoff: json_object_member(e, "spliceai", SPLICEAI_CUTOFF_FIELD).unwrap_or_default(),
            numbers: true,
        });
    }
    None
}

/// The entries of one score record, each its `|`-joined members.
pub type RecordEntries = Vec<Vec<String>>;

/// The parts of a Location the pattern `(.+):(\d+)(?:-(\d+))?` matches whole: the name
/// before the last `:`, which holds no line feed, and one or two runs of decimal digits.
fn location_span(location: &str) -> Option<(&str, &str, Option<&str>)> {
    let (name, coordinates) = location.rsplit_once(':')?;
    if name.is_empty() || name.contains('\n') {
        return None;
    }
    let (first, last) = match coordinates.split_once('-') {
        Some((first, last)) => (first, Some(last)),
        None => (coordinates, None),
    };
    (pyfmt::is_decimal(first) && last.is_none_or(pyfmt::is_decimal)).then_some((name, first, last))
}

/// `int(digits)` of a run of decimal digits: `None` for a number past `i64`. The
/// interpreter converts no run of more digits than its limit.
fn coordinate(digits: &str) -> Option<i64> {
    match pyfmt::int_checked(digits) {
        Ok(decimal) => decimal.parse().ok(),
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// The entries of the score records at a row's variant, one list per record in file order,
/// each entry its ten `|`-joined members: the allele, the symbol, four scores and four
/// positions (SpliceAI.pm:276-287, 328-332). A record stands at the row's variant when its
/// REF and ALT, spelt as Ensembl VEP
/// spells a VCF line's (`vcf_alleles_as_spelt`), give the row's Location and Allele; a
/// record whose ALT lists several alleles matches no allele (SpliceAI.pm:272). Beside them
/// stands the ground on which they are not read (a key of `SPLICEAI_GROUNDS`), `None`
/// where they are: the Location names no span, one of the records at the variant carries
/// another INFO column than the one key or an entry of another length, the records there
/// spell more than one REF, or no record stands there.
///
/// A Location whose position is past `i64` stands at no record: `SpliceAiTable::records`
/// gives every record's position as an `i64`.
pub fn spliceai_entries_at(
    e: &Entry,
    table: &SpliceAiTable,
) -> Result<(Vec<RecordEntries>, Option<String>)> {
    let refused = |ground: &str| Ok((Vec::new(), Some(ground.to_string())));
    let Some((name, first, last)) = location_span(location_of(e).unwrap_or("")) else {
        return refused("no_span");
    };
    let first = coordinate(first);
    let last = match last {
        Some(last) => coordinate(last),
        None => first,
    };
    let (Some(first), Some(last)) = (first, last) else {
        return refused("no_record");
    };
    let (low, high) = (first.min(last), first.max(last));
    let mut records: Vec<RecordEntries> = Vec::new();
    let mut references: BTreeSet<String> = BTreeSet::new();
    for (position, reference, alternate, info) in table.records(name, low - 1, high)? {
        if alternate.contains(',') {
            continue;
        }
        let (spelt, alleles) = vcf_alleles_as_spelt(&reference, &[alternate]);
        let start = i128::from(position) + i128::from(spelt != reference);
        let end = start + pyfmt::len(pyfmt::strip_chars(&spelt, "-")) as i128 - 1;
        if (start.min(end), start.max(end)) != (i128::from(low), i128::from(high))
            || alleles.first() != Some(&e.allele)
        {
            continue;
        }
        let entries: RecordEntries = pyfmt::removeprefix(&info, SPLICEAI_KEY)
            .split(',')
            .map(|entry| entry.split('|').map(str::to_string).collect())
            .collect();
        if !info.starts_with(SPLICEAI_KEY)
            || info.contains(';')
            || entries
                .iter()
                .any(|entry| entry.len() != 1 + SPLICEAI_MEMBERS.len())
        {
            return refused("other_info");
        }
        references.insert(reference);
        records.push(entries);
    }
    if references.len() > 1 {
        return refused("several_references");
    }
    let ground = records.is_empty().then(|| "no_record".to_string());
    Ok((records, ground))
}

/// The ground on which `spliceai_multi_entry_cell` refuses the differing SpliceAI field
/// `name` of a row under the score files (a key of `SPLICEAI_GROUNDS`), `None` where
/// it accepts it. Every record line of `table` at the row's variant is read, in file order
/// (`spliceai_entries_at`). It accepts where a line that holds
/// an entry naming the row's gene holds several entries; the candidate's row writes every
/// entry of the last line that names the gene which names it, in that line's order
/// (`SpliceAiValue::writes_entries`), whether or not the lines give the gene entries that
/// differ in text; and the reference's row writes
/// the whole value, behind its first allele, of the last line whose first entry names the
/// gene, SpliceAI.pm:305 keeping one line a first symbol,
/// or no SpliceAI value and no flag where no line's first entry names it. Where both rows
/// carry a value their cutoff flags are equal, or differ as one cutoff makes them
/// (`one_cutoff_gives_both_flags`): the reference's flag is the verdict of its line's
/// first entry and the candidate's that of any of the gene's entries. The row's gene is
/// the `SYMBOL` either row writes, else its transcript's symbol in the table's candidate
/// cache. Nothing is accepted on a row whose gene neither names or on a row whose gene no
/// entry names.
pub fn spliceai_multi_entry_cell_ground(
    e: &Entry,
    a: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
    table: &SpliceAiTable,
) -> Result<Option<String>> {
    let refused = |ground: &str| Ok(Some(ground.to_string()));
    let Some(ours) = spliceai_row_value(a) else {
        return refused("no_candidate_value");
    };
    let (records, ground) = spliceai_entries_at(e, table)?;
    if ground.is_some() {
        return Ok(ground);
    }
    let symbol_fields = ["SYMBOL", "gene_symbol"];
    let stated = [
        field_of(Some(e), &symbol_fields),
        field_of(a, &symbol_fields),
    ]
    .into_iter()
    .find(|value| !blank(value));
    let symbol = match (stated, &table.symbols) {
        (Some(stated), _) => Some(GeneSymbol::Text(stated)),
        (None, None) => None,
        (None, Some(symbols)) => symbols.symbol(e)?,
    };
    let symbol = match symbol {
        None => return refused("no_gene"),
        Some(GeneSymbol::NoText) => return refused("gene_not_named"),
        Some(GeneSymbol::Text(symbol)) => symbol,
    };
    let holding: Vec<&RecordEntries> = records
        .iter()
        .filter(|record| record.iter().any(|entry| entry[1] == symbol))
        .collect();
    let Some(named_last) = holding.last() else {
        return refused("gene_not_named");
    };
    let own: Vec<&Vec<String>> = named_last
        .iter()
        .filter(|entry| entry[1] == symbol)
        .collect();
    if holding.iter().all(|record| record.len() < 2) {
        return refused("lone_entry");
    }
    if !ours.writes_entries(&own) {
        return refused("candidate_value");
    }
    let theirs = spliceai_row_value(Some(e));
    let flag = field_named(name, &[SPLICEAI_CUTOFF_FIELD]);
    let led = holding.iter().rfind(|record| record[0][1] == symbol);
    let Some(led) = led else {
        let silent = theirs.is_none() && blank(&field_of(Some(e), &[SPLICEAI_CUTOFF_FIELD]));
        let taken = silent
            && blank(ev)
            && if flag {
                SPLICEAI_CUTOFF_FLAGS.contains(&av)
            } else {
                !blank(av)
            };
        return if taken {
            Ok(None)
        } else {
            refused("reference_not_silent")
        };
    };
    let joined = led
        .iter()
        .map(|entry| entry.join("|"))
        .collect::<Vec<_>>()
        .join(",");
    let whole = pyfmt::partition(&joined, "|").2;
    let Some(theirs) = theirs.filter(|theirs| theirs.writes(whole)) else {
        return refused("reference_value");
    };
    let flags_stand = if theirs.cutoff == ours.cutoff {
        !flag
    } else {
        one_cutoff_gives_both_flags(&theirs.cutoff, &led[0], &ours.cutoff, &own)
    };
    if !flags_stand {
        return refused("reference_flag");
    }
    Ok(None)
}

/// Whether the differing SpliceAI field `name` of a row is the difference a score record
/// of several entries makes (`spliceai_multi_entry_cell`): nothing without the table, and
/// under it the field `spliceai_multi_entry_cell_ground` refuses on no ground. The row
/// is counted among the table's `refusals`, under the ground where one is met.
pub fn spliceai_multi_entry_cell_row(
    e: &Entry,
    a: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
    table: Option<&SpliceAiTable>,
) -> Result<bool> {
    let Some(table) = table else {
        return Ok(false);
    };
    let ground = spliceai_multi_entry_cell_ground(e, a, name, ev, av, table)?;
    table.refusals.note(e, ground.as_deref());
    Ok(ground.is_none())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::testing::{plant, scratch};
    use crate::error::Uncaught;

    /// What a run says that ends on an input no reading of the comparison expects.
    fn uncaught<T>(run: impl FnOnce() -> T + std::panic::UnwindSafe) -> String {
        let payload = std::panic::catch_unwind(run).err().expect("the run ends");
        let Uncaught(what) = payload.downcast_ref::<Uncaught>().expect("an exception");
        what.clone()
    }

    fn row(location: Option<&str>, feature: &str, fields: &[(&str, &str)]) -> Entry {
        Entry {
            key: format!("{}|C|{feature}|Transcript", location.unwrap_or("1")),
            location: location.map(str::to_string),
            allele: "C".to_string(),
            feature: feature.to_string(),
            fields: fields
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
            ..Entry::default()
        }
    }

    const VALUE: &str = "GENEA|0.15|0.08|0.06|0.17|-6|-41|-17|-9";

    #[test]
    fn a_location_names_a_span_behind_its_last_colon() {
        assert_eq!(
            location_span("21:5001-5002"),
            Some(("21", "5001", Some("5002")))
        );
        assert_eq!(location_span("21:7000"), Some(("21", "7000", None)));
        assert_eq!(location_span("a:b:7-8"), Some(("a:b", "7", Some("8"))));
        assert_eq!(
            location_span("21:\u{664}\u{662}"),
            Some(("21", "\u{664}\u{662}", None))
        );
        for location in [
            "",
            "21",
            ":7",
            "21:",
            "21:7-",
            "21:-7",
            "21:7-8-9",
            "21:7:x",
            "2\n1:7",
            "21:\u{b2}",
        ] {
            assert_eq!(location_span(location), None, "{location:?}");
        }
        assert_eq!(coordinate("\u{664}\u{662}"), Some(42));
        assert_eq!(coordinate("99999999999999999999"), None);
        assert!(uncaught(|| coordinate(&"9".repeat(4301)))
            .starts_with("Exceeds the limit (4300 digits)"));
    }

    #[test]
    fn a_rows_value_is_its_one_field_its_nine_fields_or_its_json_object() {
        let whole = spliceai_row_value(Some(&row(
            Some("21:7"),
            "ENST1",
            &[
                ("spliceai_pred", &format!("{VALUE},C|GENEB|1")),
                ("SpliceAI_cutoff", "PASS"),
            ],
        )))
        .expect("a value");
        assert_eq!(whole.whole.as_deref(), Some(&*format!("{VALUE},C|GENEB|1")));
        assert_eq!(whole.members.len(), 9);
        assert_eq!((whole.cutoff.as_str(), whole.numbers), ("PASS", false));
        assert!(whole.writes(&format!("{VALUE},C|GENEB|1")) && !whole.writes(VALUE));
        let split = spliceai_row_value(Some(&row(
            Some("21:7"),
            "ENST1",
            &[
                ("SpliceAI_pred_SYMBOL", "GENEA"),
                ("SpliceAI_pred_DS_AG", "0.15"),
                ("SpliceAI_pred_DP_DL", "-9&C"),
            ],
        )))
        .expect("a value");
        assert_eq!(split.whole, None);
        assert_eq!(
            split.members,
            ["GENEA", "0.15", "", "", "", "", "", "", "-9,C"]
        );
        assert!(split.writes("GENEA|0.15|||||||-9,C|GENEB|1"));
        assert!(!split.writes("GENEA|0.150|||||||-9,C"));
        assert!(!split.writes("GENEA|0.15"));
        let object = spliceai_row_value(Some(&row(
            None,
            "ENST1",
            &[(
                "spliceai",
                "{\"DP_AG\":-6,\"DS_AG\":0.15,\"SYMBOL\":\"GENEA\",\"spliceai_cutoff\":\"FAIL\"}",
            )],
        )))
        .expect("a value");
        assert_eq!(
            object.members,
            ["GENEA", "0.15", "", "", "", "-6", "", "", ""]
        );
        assert_eq!((object.cutoff.as_str(), object.numbers), ("FAIL", true));
        // The members of a JSON object compare as numbers where both are numbers.
        assert!(object.writes("GENEA|0.150||||-6.0|||"));
        assert!(!object.writes("GENEA|0.15||||-6||x|"));
        for silent in [
            vec![],
            vec![("SpliceAI_pred", "-")],
            vec![("SpliceAI_cutoff", "PASS")],
        ] {
            assert_eq!(
                spliceai_row_value(Some(&row(Some("21:7"), "ENST1", &silent))),
                None
            );
        }
        assert_eq!(spliceai_row_value(None), None);
    }

    #[test]
    fn a_row_counts_once_under_the_first_ground_met_for_it() {
        let refusals = SpliceAiRefusals::new();
        let (first, second) = (
            row(Some("21:7"), "ENST1", &[]),
            row(Some("21:8"), "ENST1", &[]),
        );
        refusals.note(&first, Some("no_record"));
        refusals.note(&first, Some("no_gene"));
        refusals.note(&first, None);
        refusals.note(&second, None);
        refusals.note(&row(Some("21:9"), "ENST1", &[]), Some("no_candidate_value"));
        refusals.note(&row(Some("21:9"), "ENST2", &[]), Some("no_record"));
        let counts = refusals.counts();
        assert_eq!(
            pyjson::dumps(&counts, &pyjson::DumpOptions::default()),
            "{\"asked\": 4, \"refused\": 3, \"grounds\": {\"no_candidate_value\": 1, \"no_record\": 2}}"
        );
        assert_eq!(
            spliceai_refusal_lines(&counts, " of the reference"),
            [
                "spliceai_multi_entry_cell under --spliceai-table of the reference: rows asked about 4, rows refused 3",
                "spliceai_multi_entry_cell of the reference: rows refused 1: the candidate's row writes no SpliceAI value",
                "spliceai_multi_entry_cell of the reference: rows refused 2: no record line of the score files stands at the row's variant",
            ]
        );
        refusals.clear();
        assert_eq!(
            pyjson::dumps_compact(&refusals.counts()),
            "{\"asked\":0,\"grounds\":{},\"refused\":0}"
        );
        // Counts of another shape end the run in the interpreter's own words.
        let counts = |text: &str| pyjson::loads(text).unwrap();
        assert_eq!(
            uncaught(|| spliceai_refusal_lines(&counts("{\"asked\": 1}"), "")),
            "'refused'"
        );
        assert_eq!(
            uncaught(|| spliceai_refusal_lines(
                &counts("{\"asked\": 1, \"refused\": 1, \"grounds\": {\"other\": 1}}"),
                ""
            )),
            "'other'"
        );
        assert_eq!(
            uncaught(|| spliceai_refusal_lines(
                &counts("{\"asked\": 1, \"refused\": 1, \"grounds\": [1]}"),
                ""
            )),
            "'list' object has no attribute 'items'"
        );
        assert_eq!(
            uncaught(|| spliceai_refusal_lines(&counts("5"), "")),
            "'int' object is not subscriptable"
        );
    }

    #[test]
    fn a_refusal_of_the_table_reader_reads_under_the_score_files_option() {
        let refusal = |text: &str| SpliceAiTable::_refusal(&Error::comparator(text)).to_string();
        assert_eq!(
            refusal(
                "--revel-table: scores.synthetic.tsv.gz is not a REVEL table (a header of 7, 8 or 9 columns)"
            ),
            "--spliceai-table: scores.synthetic.tsv.gz is not a VCF (a header line of eight columns)"
        );
        assert_eq!(
            refusal("--revel-table: cannot read the tabix index a.tbi: No such file or directory"),
            "--spliceai-table: cannot read the tabix index a.tbi: No such file or directory"
        );
        assert_eq!(refusal("no option"), "--spliceai-table: ");
    }

    #[test]
    fn a_cache_names_the_gene_of_a_rows_transcript() {
        let root = scratch("spliceai-symbols");
        let refused = TranscriptSymbols::new(&root).unwrap_err().to_string();
        assert_eq!(
            refused,
            format!(
                "--json-cache: {} has no transcripts/ directory",
                root.display()
            )
        );
        plant(
            &root,
            "transcripts/21/1-1000000.json",
            r#"[{"stable_id": "ENST1", "gene_symbol": "GENEA"}, null,
                {"stable_id": "ENST2.3", "gene_symbol": "GENEB"},
                {"stable_id": "ENST4", "gene_symbol": 5},
                {"stable_id": "ENST5", "gene_symbol": ""},
                {"stable_id": 6, "gene_symbol": "GENEC"},
                {"gene_symbol": "GENEC"}]"#,
        );
        plant(
            &root,
            "transcripts/21/2.json",
            r#"{"1-5": [{"stable_id": "ENST1", "gene_symbol": "GENED"}]}"#,
        );
        plant(
            &root,
            "transcripts/22/1.json",
            r#"[{"stable_id": "ENST9", "gene_symbol": "GENEZ"}]"#,
        );
        let symbols = TranscriptSymbols::new(&root).expect("a cache");
        let symbol = |location: Option<&str>, feature: &str| {
            symbols
                .symbol(&row(location, feature, &[]))
                .expect("the shards read")
        };
        let text = |symbol: &str| Some(GeneSymbol::Text(symbol.to_string()));
        // A later shard of the directory names the transcript last.
        assert_eq!(symbol(Some("21:7"), "ENST1"), text("GENED"));
        assert_eq!(symbol(Some("21:7"), "ENST1.7"), text("GENED"));
        assert_eq!(symbol(Some("21:7"), "ENST2.3"), text("GENEB"));
        assert_eq!(symbol(Some("21:7"), "ENST2"), None);
        assert_eq!(symbol(Some("21:7"), "ENST4"), Some(GeneSymbol::NoText));
        assert_eq!(symbol(Some("21:7"), "ENST5"), None);
        assert_eq!(symbol(Some("21:7"), "6"), None);
        assert_eq!(symbol(Some("21:7"), "ENST9"), None);
        // A row whose Location names no directory of the cache reads every shard.
        assert_eq!(symbol(Some("X:7"), "ENST9"), text("GENEZ"));
        assert_eq!(symbol(None, "ENST9"), text("GENEZ"));
        let record_row = Entry {
            allele_location: Some("22:7".to_string()),
            ..row(None, "ENST9", &[])
        };
        assert_eq!(symbols.symbol(&record_row).unwrap(), text("GENEZ"));
        assert_eq!(symbols._symbols.borrow().len(), 3);

        plant(
            &root,
            "transcripts/7/1.json",
            r#"[{"stable_id": ["ENST1"], "gene_symbol": ""}, {"stable_id": ["ENST1"], "gene_symbol": "GENEA"}]"#,
        );
        assert_eq!(
            uncaught(std::panic::AssertUnwindSafe(|| symbols.symbol(&row(
                Some("7:7"),
                "ENST1",
                &[]
            )))),
            "cannot use 'list' as a dict key (unhashable type: 'list')"
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
