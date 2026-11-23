// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The SpliceAI plugin's classes and the score files they read.
//!
//! `spliceai_chopped_multiallelic`: for a record whose allele string holds `-`,
//! SpliceAI.pm:259-266 rebuild every allele with the anchor base before matching, so a
//! substitution beside a deletion is compared as a two-base pair, which equals no
//! single-base row (:272), where the engine compares each allele at its own position. The
//! shape is a SpliceAI field the engine alone fills on a substitution allele's row of a
//! record whose rows carry an allele written `-`.
//!
//! `spliceai_multi_entry_cell`: SpliceAI.pm:326-332 remove `SpliceAI=` from a score
//! record's INFO column and split the rest on `|` into the allele, the first symbol and
//! everything after it, and :305 file the record under that symbol alone. A record whose
//! key holds several comma-joined entries so gives the rows of the first entry's symbol the
//! whole remainder as their value (the nine fields of the other formats cut from it,
//! :276-287) and the rows of a later entry's symbol nothing (:314-316), where the engine
//! writes a row every entry of the record that names the row's gene, in the record's
//! order. Where a record line stands more than once at an
//! allele, :305 keep one line a first symbol, the last, so the reference writes a row the
//! whole value of the last line its gene leads, or nothing where it leads none, and the
//! engine the gene's entries of the last line that names the gene. A cutoff flag is the
//! reference's verdict on the first entry of its line (:336-338) and the engine's on any of
//! the gene's entries. The shape is read from the
//! score files themselves (`SpliceAiTable`) and nothing is accepted without them. A row's
//! gene is its `SYMBOL` where a row writes one, else its transcript's symbol in the
//! engine's JSON cache (`TranscriptSymbols`); a row whose gene neither names is not
//! accepted.

use super::cache::TranscriptShards;
use super::entry::{as_number, blank, field_named, field_of, Entry, Row};
use super::shared::{at_its_location, record_alleles, record_indel_row, vcf_alleles_as_spelt};
use super::tabix::{IndexedTable, TableError};
use super::{Consult, Refusal, Res, SPLICEAI_MULTI_ENTRY_CELL};
use crate::pytext::{self, Json};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// The prefix `entries_at` requires of a record's INFO cell and strips before its entries.
const KEY: &str = "SpliceAI=";
const VALUE_FIELD: &str = "SpliceAI_pred";
const CUTOFF_FIELD: &str = "SpliceAI_cutoff";
/// The members of a row's SpliceAI value in the order SpliceAI.pm:276-287 read them from a
/// score entry after its allele: the symbol, four scores and four positions.
pub const MEMBERS: [&str; 9] = [
    "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
];
const CUTOFF_FLAGS: [&str; 2] = ["PASS", "FAIL"];
const SYMBOL_FIELDS: [&str; 2] = ["SYMBOL", "gene_symbol"];
/// The JSON format's object of the plugin's fields.
const JSON_OBJECT: &str = "spliceai";

/// The grounds on which `spliceai_multi_entry_cell` refuses a row it is asked about under
/// the score files, in the order they are tested: a key and what the rule read.
pub const GROUNDS: [(&str, &str); 12] = [
    ("no_candidate_value", "the candidate's row writes no SpliceAI value"),
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

/// The rows a class was asked about and refused, and the refused rows under each ground in
/// the order the grounds are tested.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RefusalCounts {
    pub asked: u64,
    pub refused: u64,
    pub grounds: Vec<(String, u64)>,
}

/// What a cache gives as a transcript's gene symbol: a text, or a value that is none and
/// so equals no score entry's symbol.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Symbol {
    Named(String),
    NoText,
}

/// Whether Python reads a JSON value as true.
fn truthy(value: &Json) -> bool {
    match value {
        Json::Null => false,
        Json::Bool(held) => *held,
        Json::Int(digits) => digits.bytes().any(|b| matches!(b, b'1'..=b'9')),
        Json::Float(number) => *number != 0.0,
        Json::Str(text) => !text.is_empty(),
        Json::List(items) => !items.is_empty(),
        Json::Object(members) => !members.is_empty(),
    }
}

type Symbols = Arc<HashMap<String, Symbol>>;

/// The gene symbol of each transcript of an engine's JSON cache, by stable id: the symbol
/// SpliceAI.pm:314 looks a transcript's score record up under. A chromosome's shards are
/// read on the first lookup of a row whose Location names its directory, every shard on the
/// first lookup of any other entry.
pub struct TranscriptSymbols {
    shards: TranscriptShards,
    read: Mutex<HashMap<Option<String>, Symbols>>,
}

impl TranscriptSymbols {
    /// The symbols of the cache at `json_cache`, refused where it holds no `transcripts/`
    /// directory.
    pub fn new(json_cache: &Path) -> Res<TranscriptSymbols> {
        Ok(TranscriptSymbols {
            shards: TranscriptShards::new(json_cache)?,
            read: Mutex::new(HashMap::new()),
        })
    }

    /// The symbols of one chromosome directory's transcripts, or of every directory's: a
    /// transcript counts where its `stable_id` and its `gene_symbol` are both filled, and a
    /// later transcript of a stable id replaces an earlier one's symbol.
    fn of(&self, chromosome: Option<&str>) -> Res<Symbols> {
        let mut read = self.read.lock().expect("lock");
        let key = chromosome.map(str::to_string);
        if let Some(symbols) = read.get(&key) {
            return Ok(Arc::clone(symbols));
        }
        let mut symbols = HashMap::new();
        let mut take = |transcript: Json| {
            if let (Some(Json::Str(id)), Some(symbol)) =
                (transcript.get("stable_id"), transcript.get("gene_symbol"))
            {
                if !id.is_empty() && truthy(symbol) {
                    let symbol = match symbol {
                        Json::Str(text) => Symbol::Named(text.clone()),
                        _ => Symbol::NoText,
                    };
                    symbols.insert(id.clone(), symbol);
                }
            }
            Ok(())
        };
        match chromosome {
            Some(chromosome) => self.shards.of_chromosome(chromosome, &mut take)?,
            None => self.shards.of_every_chromosome(&mut take)?,
        }
        let symbols = Arc::new(symbols);
        read.insert(key, Arc::clone(&symbols));
        Ok(symbols)
    }

    /// The gene symbol of the entry's transcript, or `None` for a transcript the cache does
    /// not hold or gives no symbol; a Feature carrying a version the cache's stable id
    /// lacks is looked up without it.
    fn symbol(&self, e: &Row) -> Res<Option<Symbol>> {
        let located = at_its_location(e);
        let chromosome = located
            .location
            .as_deref()
            .filter(|location| !location.is_empty())
            .map(|location| location.split(':').next().unwrap_or(""))
            .filter(|chromosome| self.shards.root().join(chromosome).is_dir());
        let symbols = self.of(chromosome)?;
        let unversioned = e.feature.split('.').next().unwrap_or("");
        Ok(symbols
            .get(&e.feature)
            .or_else(|| symbols.get(unversioned))
            .cloned())
    }
}

/// The cells `SpliceAiTable::records` keeps of a record line it reads: POS, REF, ALT and INFO.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScoreRecord {
    pub position: i128,
    pub reference: String,
    pub alternate: String,
    pub info: String,
}

/// A row as the refusals count it: its entry's key, Location, allele and feature.
type CountedRow = (String, Option<String>, String, String);

/// The rows `spliceai_multi_entry_cell` is asked about under the score files and, of them,
/// the rows it refuses, each under the first ground met for it. A row counts once however
/// many of its fields the class is asked about.
#[derive(Default)]
struct Refusals {
    asked: HashSet<CountedRow>,
    refused: HashMap<CountedRow, &'static str>,
}

/// The score files a run's SpliceAI plugin read (the `snv=` file and the `indel=` file),
/// each a VCF with its tabix index beside it (SpliceAI.pm:114-122), read through that index
/// at the positions a row names and never whole. The table also holds the engine's JSON
/// cache, whose
/// transcripts name the gene of a row that states no `SYMBOL`, and counts the rows
/// `spliceai_multi_entry_cell` is asked about under the files and those it refuses.
pub struct SpliceAiTable {
    files: Vec<IndexedTable>,
    symbols: Option<TranscriptSymbols>,
    refusals: Mutex<Refusals>,
}

fn refusal(error: TableError) -> Refusal {
    Refusal(format!("--spliceai-table: {}", error.0))
}

impl SpliceAiTable {
    /// The score files at `paths`, each with its tabix index beside it, and the
    /// cache at `json_cache`. A cache without a `transcripts/` directory is refused, then
    /// the first file that has no readable index, is no BGZF file or whose header line has
    /// another count of columns than a VCF's.
    pub fn new(paths: &[PathBuf], json_cache: Option<&Path>) -> Result<SpliceAiTable, Refusal> {
        let symbols = json_cache.map(TranscriptSymbols::new).transpose()?;
        let mut files = Vec::with_capacity(paths.len());
        for path in paths {
            let table = IndexedTable::open(path).map_err(refusal)?;
            if !matches!(table.header().len(), 8 | 9) {
                return Err(Refusal(format!(
                    "--spliceai-table: {} is not a VCF (a header line of eight columns)",
                    path.display()
                )));
            }
            files.push(table);
        }
        Ok(SpliceAiTable {
            files,
            symbols,
            refusals: Mutex::new(Refusals::default()),
        })
    }

    /// Reads a chromosome's shards of the cache on `threads` threads from the next
    /// chromosome read on.
    pub fn read_on(&self, threads: usize) {
        if let Some(symbols) = &self.symbols {
            symbols.shards.read_on(threads);
        }
    }

    /// The records at the positions `low` to `high` of `chromosome`, the files in the order
    /// given and each in file order; a line that lacks one of the four cells, or whose POS
    /// is no number, is left out.
    pub fn records(&self, chromosome: &str, low: i128, high: i128) -> Res<Vec<ScoreRecord>> {
        let mut found = Vec::new();
        for table in &self.files {
            for row in table.rows(chromosome, low, high).map_err(refusal)?.iter() {
                let filled = |column: usize| row.cell(column).filter(|cell| !cell.is_empty());
                let (Some(reference), Some(alternate), Some(info)) =
                    (filled(3), filled(4), filled(7))
                else {
                    continue;
                };
                if !pytext::isdigit(&row.position) {
                    continue;
                }
                let Some(position) = pytext::int(&row.position) else {
                    continue;
                };
                found.push(ScoreRecord {
                    position,
                    reference: reference.to_string(),
                    alternate: alternate.to_string(),
                    info: info.to_string(),
                });
            }
        }
        Ok(found)
    }

    /// Forgets the rows counted so far.
    pub fn clear_refusals(&self) {
        let mut refusals = self.refusals.lock().expect("lock");
        refusals.asked.clear();
        refusals.refused.clear();
    }

    /// The rows counted since the counts were cleared: those asked about, those refused
    /// and the refused rows of each ground that refused any, in the order of `GROUNDS`.
    pub fn refusal_counts(&self) -> RefusalCounts {
        let refusals = self.refusals.lock().expect("lock");
        let rows = |ground: &str| {
            refusals
                .refused
                .values()
                .filter(|met| **met == ground)
                .count()
        };
        RefusalCounts {
            asked: refusals.asked.len() as u64,
            refused: refusals.refused.len() as u64,
            grounds: GROUNDS
                .iter()
                .map(|(ground, _)| (ground.to_string(), rows(ground) as u64))
                .filter(|(_, rows)| *rows > 0)
                .collect(),
        }
    }

    /// Counts the row of `e` as asked about and, under a `ground`, as refused.
    fn note(&self, e: &Entry, ground: Option<&'static str>) {
        let row = (
            e.key.clone(),
            e.location.clone(),
            e.allele.clone(),
            e.feature.clone(),
        );
        let mut refusals = self.refusals.lock().expect("lock");
        if let Some(ground) = ground {
            refusals.refused.entry(row.clone()).or_insert(ground);
        }
        refusals.asked.insert(row);
    }
}

/// The lines that state `counts` for the reader of a run: the rows the class was asked
/// about and refused, then the refused rows of each ground. `whose` names the rows' side
/// of a summary of several.
pub fn refusal_lines(counts: &RefusalCounts, whose: &str) -> Vec<String> {
    let mut lines = vec![format!(
        "{SPLICEAI_MULTI_ENTRY_CELL} under --spliceai-table{whose}: rows asked about {}, rows refused {}",
        counts.asked, counts.refused
    )];
    for (ground, rows) in &counts.grounds {
        let read = GROUNDS
            .iter()
            .find(|(key, _)| key == ground)
            .map_or(ground.as_str(), |(_, read)| read);
        lines.push(format!(
            "{SPLICEAI_MULTI_ENTRY_CELL}{whose}: rows refused {rows}: {read}"
        ));
    }
    lines
}

/// A row's SpliceAI value: its `SpliceAI_pred` text where its format writes the one field,
/// else `None`; its nine members (`MEMBERS`); its cutoff flag, empty for none; and whether
/// the members are JSON values, whose numbers compare as numbers.
struct RowValue {
    whole: Option<String>,
    members: Vec<String>,
    cutoff: String,
    numbers: bool,
}

impl RowValue {
    /// Whether the row writes `value`, a `|`-joined value as the default format prints it:
    /// as its whole text where the row carries the one field, else as its nine members, the
    /// parts of `value` past the ninth being those SpliceAI.pm:276-287 drop.
    fn writes(&self, value: &str) -> bool {
        if let Some(whole) = &self.whole {
            return whole == value;
        }
        let parts: Vec<&str> = value.split('|').take(MEMBERS.len()).collect();
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
    fn writes_entries(&self, entries: &[&Vec<String>]) -> bool {
        if let Some(whole) = &self.whole {
            let spelt: Vec<String> = entries
                .iter()
                .enumerate()
                .map(|(place, entry)| entry[usize::from(place == 0)..].join("|"))
                .collect();
            return *whole == spelt.join(",");
        }
        self.members.len() == MEMBERS.len()
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
/// from (SpliceAI.pm:336-338), and the engine's `ours` by the largest over `own`, the
/// gene's entries of the engine's line, a row passing where any of them reaches the
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

/// A JSON value as one comparable string: a string itself, an array its items joined with
/// `,`, anything else its compact JSON text.
fn json_scalar(value: &Json) -> String {
    match value {
        Json::Str(text) => text.clone(),
        Json::List(items) => items.iter().map(json_scalar).collect::<Vec<_>>().join(","),
        _ => pytext::json_compact(value),
    }
}

/// The objects of a field whose text is one JSON object or several joined with `,`; none
/// for any other text.
fn field_objects(text: &str) -> Vec<Vec<(String, Json)>> {
    if !text.starts_with('{') || !text.ends_with('}') {
        return Vec::new();
    }
    match pytext::json_loads(&format!("[{text}]")) {
        Ok(Json::List(items)) => items
            .into_iter()
            .filter_map(|item| match item {
                Json::Object(members) => Some(members),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// The value of the member `key`, in any case, of the first of `objects` that carries it.
fn object_member(objects: &[Vec<(String, Json)>], key: &str) -> Option<String> {
    let wanted = key.to_lowercase();
    objects
        .iter()
        .flatten()
        .find(|(name, _)| name.to_lowercase() == wanted)
        .map(|(_, value)| json_scalar(value))
}

/// The SpliceAI value of a row, read from its `SpliceAI_pred` field, else from its nine
/// `SpliceAI_pred_` fields (a CSQ field's `&` read as the comma it stands for), else from
/// the members of its JSON `spliceai` object; `None` for a row that writes none.
fn row_value(e: &Entry) -> Option<RowValue> {
    let written = |name: &str| Some(field_of(e, &[name])).filter(|value| !blank(value));
    let cutoff = || written(CUTOFF_FIELD).unwrap_or("").to_string();
    if let Some(whole) = written(VALUE_FIELD) {
        return Some(RowValue {
            whole: Some(whole.to_string()),
            members: Vec::new(),
            cutoff: cutoff(),
            numbers: false,
        });
    }
    let split: Vec<String> = MEMBERS
        .iter()
        .map(|member| {
            written(&format!("{VALUE_FIELD}_{member}"))
                .unwrap_or("")
                .replace('&', ",")
        })
        .collect();
    if split.iter().any(|member| !member.is_empty()) {
        return Some(RowValue {
            whole: None,
            members: split,
            cutoff: cutoff(),
            numbers: false,
        });
    }
    let objects = field_objects(e.fields.get_or(JSON_OBJECT, ""));
    let members: Vec<String> = MEMBERS
        .iter()
        .map(|member| object_member(&objects, member).unwrap_or_default())
        .collect();
    members
        .iter()
        .any(|member| !member.is_empty())
        .then(|| RowValue {
            whole: None,
            members,
            cutoff: object_member(&objects, CUTOFF_FIELD).unwrap_or_default(),
            numbers: true,
        })
}

/// The chromosome and the one position or the two of a Location written `chromosome:n` or
/// `chromosome:n-m`, the digits of any script; `None` for any other text.
fn location_span(location: &str) -> Option<(&str, &str, &str)> {
    let (chromosome, span) = location.rsplit_once(':')?;
    let (first, last) = span.split_once('-').unwrap_or((span, span));
    (!chromosome.is_empty()
        && !chromosome.contains('\n')
        && pytext::isdecimal(first)
        && pytext::isdecimal(last))
    .then_some((chromosome, first, last))
}

/// The entries of one score record, each its ten `|`-joined members: the allele, the
/// symbol, four scores and four positions (SpliceAI.pm:276-287, 328-332).
type Entries = Vec<Vec<String>>;

/// The entries of the score records at a row's variant, one list per record in file order.
/// A record stands at the row's variant when its REF and ALT, spelt as Ensembl VEP spells a
/// VCF line's (`vcf_alleles_as_spelt`), give the row's Location and allele; a record whose
/// ALT lists several alleles matches no allele (SpliceAI.pm:272). Beside them stands the
/// ground on which they are not read, `None` where they are: the Location names no span,
/// one of the records at the variant carries another INFO column than the one key or an
/// entry of another length, the records there spell more than one REF, or no record stands
/// there. A position outside `i128` is one no record stands at.
fn entries_at(e: &Row, table: &SpliceAiTable) -> Res<(Vec<Entries>, Option<&'static str>)> {
    let located = at_its_location(e);
    let Some((chromosome, first, last)) = location_span(located.location.as_deref().unwrap_or(""))
    else {
        return Ok((Vec::new(), Some("no_span")));
    };
    let (Some(first), Some(last)) = (pytext::int(first), pytext::int(last)) else {
        return Ok((Vec::new(), Some("no_record")));
    };
    let (low, high) = (first.min(last), first.max(last));
    let mut records: Vec<Entries> = Vec::new();
    let mut references: BTreeSet<String> = BTreeSet::new();
    for record in table.records(chromosome, low - 1, high)? {
        if record.alternate.contains(',') {
            continue;
        }
        let (spelt, alleles) = vcf_alleles_as_spelt(&record.reference, &[&record.alternate]);
        let bases = spelt.trim_matches('-').chars().count() as i128;
        // A record whose start or end lies past the largest `i128` stands at no Location.
        let start = record
            .position
            .checked_add(i128::from(spelt != record.reference));
        let end = start.and_then(|start| start.checked_add(bases - 1));
        let (Some(start), Some(end)) = (start, end) else {
            continue;
        };
        if (start.min(end), start.max(end)) != (low, high) || alleles[0] != e.allele {
            continue;
        }
        let entries: Entries = record
            .info
            .strip_prefix(KEY)
            .unwrap_or(&record.info)
            .split(',')
            .map(|entry| entry.split('|').map(str::to_string).collect())
            .collect();
        if !record.info.starts_with(KEY)
            || record.info.contains(';')
            || entries.iter().any(|entry| entry.len() != 1 + MEMBERS.len())
        {
            return Ok((Vec::new(), Some("other_info")));
        }
        references.insert(record.reference);
        records.push(entries);
    }
    if references.len() > 1 {
        return Ok((Vec::new(), Some("several_references")));
    }
    let ground = records.is_empty().then_some("no_record");
    Ok((records, ground))
}

/// The ground on which `spliceai_multi_entry_cell` refuses the differing SpliceAI field
/// `name` of a row under the score files (a key of `GROUNDS`), `None` where it accepts it.
/// Every record line of `table` at the row's variant is read, in file order. It accepts where
/// a line holding an entry that names the row's gene holds several; the engine's row writes
/// every entry of the last line that names the gene which names it, in that line's order
/// (`RowValue::writes_entries`); and the
/// reference's row writes the whole value, behind its first allele, of the last line whose
/// first entry names the gene, or no SpliceAI value and
/// no flag where no line's first entry names it. Where both rows carry a value their cutoff
/// flags are equal, or differ as one cutoff makes them (`one_cutoff_gives_both_flags`): the
/// reference's flag is the verdict of its line's first entry and the engine's that of any of
/// the gene's entries. The row's gene is the `SYMBOL` either row
/// writes, else its transcript's symbol in the table's cache. `ev` and `av` are the two
/// values with a blank one as `-`.
fn ground(
    e: &Row,
    a: Option<&Row>,
    name: &str,
    ev: &str,
    av: &str,
    table: &SpliceAiTable,
) -> Res<Option<&'static str>> {
    let Some((a, ours)) = a.and_then(|a| Some((a, row_value(a)?))) else {
        return Ok(Some("no_candidate_value"));
    };
    let (records, met) = entries_at(e, table)?;
    if met.is_some() {
        return Ok(met);
    }
    let stated = [field_of(e, &SYMBOL_FIELDS), field_of(a, &SYMBOL_FIELDS)]
        .into_iter()
        .find(|value| !blank(value));
    let symbol = match (stated, &table.symbols) {
        (Some(stated), _) => Some(Symbol::Named(stated.to_string())),
        (None, Some(symbols)) => symbols.symbol(e)?,
        (None, None) => None,
    };
    let Some(symbol) = symbol else {
        return Ok(Some("no_gene"));
    };
    let names = |entry: &Vec<String>| matches!(&symbol, Symbol::Named(gene) if *gene == entry[1]);
    let holding: Vec<&Entries> = records
        .iter()
        .filter(|record| record.iter().any(names))
        .collect();
    let Some(named_last) = holding.last() else {
        return Ok(Some("gene_not_named"));
    };
    let own: Vec<&Vec<String>> = named_last.iter().filter(|entry| names(entry)).collect();
    if holding.iter().all(|record| record.len() < 2) {
        return Ok(Some("lone_entry"));
    }
    if !ours.writes_entries(&own) {
        return Ok(Some("candidate_value"));
    }
    let theirs = row_value(e);
    let flag = field_named(name, &[CUTOFF_FIELD]);
    let Some(led) = holding.iter().rfind(|record| names(&record[0])) else {
        let silent = theirs.is_none() && blank(field_of(e, &[CUTOFF_FIELD]));
        let filled = if flag {
            CUTOFF_FLAGS.contains(&av)
        } else {
            !blank(av)
        };
        let taken = silent && blank(ev) && filled;
        return Ok((!taken).then_some("reference_not_silent"));
    };
    let joined = led
        .iter()
        .map(|entry| entry.join("|"))
        .collect::<Vec<_>>()
        .join(",");
    let whole = pytext::partition(&joined, "|").2;
    let Some(theirs) = theirs.filter(|theirs| theirs.writes(whole)) else {
        return Ok(Some("reference_value"));
    };
    let flags_stand = if theirs.cutoff == ours.cutoff {
        !flag
    } else {
        one_cutoff_gives_both_flags(&theirs.cutoff, &led[0], &ours.cutoff, &own)
    };
    if !flags_stand {
        return Ok(Some("reference_flag"));
    }
    Ok(None)
}

/// `spliceai_chopped_multiallelic`: `ev` and `av` are the two values with a blank one as
/// `-`. Nothing is accepted without the rows of the reference's record.
pub fn chopped_multiallelic(c: &Consult, ev: &str, av: &str) -> Res<bool> {
    let Some(record_rows) = c.record_rows else {
        return Ok(false);
    };
    if !blank(ev) || blank(av) {
        return Ok(false);
    }
    let substitution = !c.e.allele.is_empty() && c.e.allele.bytes().all(|b| b"ACGTN".contains(&b));
    if !substitution || record_indel_row(c.e, record_rows)? {
        return Ok(false);
    }
    Ok(record_alleles(c.e, record_rows).contains("-"))
}

/// `spliceai_multi_entry_cell`: nothing without the score files, and under them the field
/// `ground` refuses on no ground. The row is counted among the table's refusals, under
/// the ground where one is met.
pub fn multi_entry_cell(c: &Consult, ev: &str, av: &str) -> Res<bool> {
    let Some(table) = c.spliceai_table else {
        return Ok(false);
    };
    let met = ground(c.e, c.a, c.name, ev, av, table)?;
    table.note(c.e, met);
    Ok(met.is_none())
}
