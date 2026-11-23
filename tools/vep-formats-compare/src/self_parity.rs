// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The self-parity mode.

use crate::consts::{
    COMPLEMENT, CONSEQUENCE_RANK, JSON_CONSEQUENCE_LISTS, JSON_OBJECT_FIELDS, JSON_RENAME,
    POSITION_RE, RECORD_LEVEL_COLUMNS, RECORD_LEVEL_COLUMN_RE,
};
use crate::entry::Entry;
use crate::error::Result;
use crate::input::{dedup_keys, line_slices, numberified_json_text};
use crate::manifest::{decode_extra, repr_text_list};
use crate::parse::json::{
    _json_float, json_compact, json_objects, json_record_label, json_scalar, parse_json_under,
    scalar_frames, SELF_PARITY_CALLERS,
};
use crate::parse::text::{parse_default_under, parse_tab_under};
use crate::parse::vcf::{csq_format_fields, header_columns, parse_vcf, vcf_line_key, InputIndex};
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::shared::{dash, feature_kind, ABSENT};
use indexmap::{Equivalent, IndexMap};
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::LazyLock;

/// What `perl_numeric_or_text` makes of a word.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntOrStr {
    Int(i64),
    /// A whole number beyond `i64`, as its decimal text.
    BigInt(String),
    Str(String),
}

impl IntOrStr {
    fn into_value(self) -> Value {
        match self {
            IntOrStr::Int(whole) => Value::Int(whole),
            IntOrStr::BigInt(decimal) => Value::BigInt(decimal),
            IntOrStr::Str(text) => Value::Str(text),
        }
    }
}

/// Config.pm:501-508 at release/116.3: --vcf implies --symbol, --biotype and --numbers, so a
/// VCF run carries these columns where a default run without those flags does not.
///
/// Membership only.
pub const VCF_IMPLIED_COLUMNS: &[&str] = &[
    "BIOTYPE",
    "EXON",
    "HGNC_ID",
    "INTRON",
    "SYMBOL",
    "SYMBOL_SOURCE",
];

pub const POSITION_PART_RE_PATTERN: &str = r"^\d+(/\d)?\d*(?:\n?\z)";
pub static POSITION_PART_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(POSITION_PART_RE_PATTERN).expect("the pattern compiles"));

/// The letters of the range are the 52 of ASCII and the four more that the interpreter's
/// case-insensitive `[a-z]` takes: U+0130, U+0131, U+017F and U+212A.
pub const TOOL_SCORE_RE_PATTERN: &str = r"(?i)([a-zİı_]+)?\(?([\d.]+)?\)?";
pub static TOOL_SCORE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(TOOL_SCORE_RE_PATTERN).expect("the pattern compiles"));

pub const DOMAIN_RE_PATTERN: &str = r"([\p{L}\p{N}_]+):([\p{L}\p{N}_]+)";
pub static DOMAIN_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(DOMAIN_RE_PATTERN).expect("the pattern compiles"));

pub const PERL_NUMBER_RE_PATTERN: &str =
    r"^[\s\x1c-\x1f]*[+-]?(\d+\.?\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?)";
pub static PERL_NUMBER_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(PERL_NUMBER_RE_PATTERN).expect("the pattern compiles"));

/// The call depth the self-parity mode reads its outputs and compares its entries at. A
/// JSON value is written as text by calls that nest once or twice a level of its nesting,
/// and the run ends in the words of `RecursionError` where those calls would pass the
/// interpreter's recursion limit (`pyjson::RECURSION_LIMIT`) less the depth the writing
/// runs at.
const METHOD_FRAMES: usize = SELF_PARITY_CALLERS;

/// What `RecursionError` says.
const RECURSION: &str = "maximum recursion depth exceeded";

static NULL: Value = Value::Null;
static NO_OBJECT: LazyLock<Value> = LazyLock::new(|| Value::Dict(IndexMap::new()));

/// The comparison of one default-format entry with its counterpart of another format: the
/// (column, detail) of every difference; the columns it looks at join the set.
pub type CompareFields<'a> =
    dyn Fn(&Entry, &Entry, &mut BTreeSet<String>) -> Vec<(String, String)> + 'a;

/// `json_scalar(value)` called by a function that is `callers` frames deep; `None` where
/// the depth of the call is not known and no value is refused.
fn scalar_under(value: &Value, callers: Option<usize>) -> String {
    if callers.is_some_and(|callers| callers + scalar_frames(value) > pyjson::RECURSION_LIMIT) {
        crate::uncaught!("{RECURSION}");
    }
    json_scalar(value)
}

/// `0 + $s`: the leading number of a string, 0 when it has none.
///
/// The number is read without its sign, so a negative number comes back positive.
pub fn perl_numify(s: &str) -> String {
    let Some(text) = PERL_NUMBER_RE
        .captures(s)
        .and_then(|found| found.get(1))
        .map(|number| number.as_str())
    else {
        return "0".to_string();
    };
    let Some(value) = pyfmt::float_of(text) else {
        return "0".to_string();
    };
    if value.is_infinite() {
        return if value > 0.0 { "Inf" } else { "-Inf" }.to_string();
    }
    if value.fract() == 0.0 && !text.contains(['e', 'E']) {
        if let Ok(whole) = pyfmt::int_of_f64(value) {
            return whole;
        }
    }
    _json_float(value)
}

/// A `\w+` capture as `numberify` leaves it: ASCII digits alone become an integer (Perl's
/// `+ 0` on a string `looks_like_number` accepts); anything else stays text.
///
/// A run of more digits than the interpreter converts ends the run.
pub fn perl_numeric_or_text(word: &str) -> IntOrStr {
    if word.is_empty() || !word.bytes().all(|byte| byte.is_ascii_digit()) {
        return IntOrStr::Str(word.to_string());
    }
    match pyfmt::int_checked(word) {
        Ok(decimal) => match decimal.parse() {
            Ok(whole) => IntOrStr::Int(whole),
            Err(_) => IntOrStr::BigInt(decimal),
        },
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// `float(value)`, `None` for a value that gives `TypeError` or `ValueError`; a whole number
/// beyond the largest float ends the run.
fn float_or_caught(value: &Value) -> Option<f64> {
    match pyjson::py_float(value) {
        Ok(float) => Some(float),
        Err(pyjson::FloatError::Overflow(what)) => crate::uncaught!("{what}"),
        Err(_) => None,
    }
}

/// `Utils.pm:141-166`: `start`, `start-end`, or the pair low to high.
///
/// A JSON `null` is an absent bound. Each bound is printed as the interpreter prints the
/// JSON value, whatever its type.
pub fn format_coords(start: &Value, end: &Value) -> String {
    let shown = pyjson::py_str;
    match (start, end) {
        (Value::Null, Value::Null) => "-".to_string(),
        (start, Value::Null) => format!("{}-?", shown(start)),
        (Value::Null, end) => format!("?-{}", shown(end)),
        (start, end) => {
            let bounds = float_or_caught(start).and_then(|s| Some((s, float_or_caught(end)?)));
            match bounds {
                Some((s, e)) if s > e => format!("{}-{}", shown(end), shown(start)),
                Some((s, e)) if s == e => shown(start),
                _ => format!("{}-{}", shown(start), shown(end)),
            }
        }
    }
}

/// `Uploaded_variation` for a record of VCF input as `OutputFactory.pm:883-930`
/// derives it from the JSON object's `id`, `seq_region_name`, `start`,
/// `allele_string` and `input`: the ID when it is one, else `chr_pos_ref/alt`.
pub fn uploaded_variation(obj: &Value) -> String {
    uploaded_variation_under(obj, None)
}

/// `uploaded_variation(obj)` called by a function that is `callers` frames deep.
fn uploaded_variation_under(obj: &Value, callers: Option<usize>) -> String {
    let scalar = |value: &Value| scalar_under(value, callers.map(|callers| callers + 1));
    let chrom = obj.get("seq_region_name").map_or_else(String::new, scalar);
    let start = obj.get("start").unwrap_or(&NULL);
    let allele_string = obj.get("allele_string").map_or_else(String::new, scalar);
    let whole = float_or_caught(start).and_then(|float| match pyfmt::int_of_f64(float) {
        Ok(whole) => Some(whole),
        Err(pyfmt::FloatToIntError::Nan) => None,
        Err(error) => crate::uncaught!("{error}"),
    });
    let as_given = || format!("{chrom}_{}_{allele_string}", pyjson::py_str(start));
    let var_name = match whole {
        Some(whole) => format!("{chrom}_{whole}_{allele_string}"),
        None => as_given(),
    };
    let vn = obj
        .get("id")
        .filter(|id| !matches!(id, Value::Null))
        .map(scalar);
    if let Some(vn) = &vn {
        if !vn.is_empty() && vn != "." && *vn != var_name {
            return vn.clone();
        }
    }
    if let Some(line) = obj.get("input").and_then(Value::as_str) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() >= 5 {
            let alt = cols[4].replace(',', "/");
            let alleles = if alt.is_empty() {
                cols[3].to_string()
            } else {
                format!("{}/{alt}", cols[3])
            };
            return format!("{}_{}_{alleles}", cols[0], perl_numify(cols[1]));
        }
    }
    match vn {
        Some(vn) if vn != "." && vn != var_name => vn,
        _ => as_given(),
    }
}

/// `Location` as `OutputFactory.pm:931` prints it from the JSON object's
/// `seq_region_name`, `start` and `end`.
pub fn record_location(obj: &Value) -> String {
    record_location_under(obj, None)
}

/// `record_location(obj)` called by a function that is `callers` frames of the
/// interpreter's call stack deep: a `seq_region_name` nested deeper than the frames left
/// above the function ends the run in the words of `RecursionError`. `None` refuses no
/// value.
pub fn record_location_under(obj: &Value, callers: Option<usize>) -> String {
    let chrom = obj.get("seq_region_name").map_or_else(String::new, |name| {
        scalar_under(name, callers.map(|callers| callers + 1))
    });
    let coords = format_coords(
        obj.get("start").unwrap_or(&NULL),
        obj.get("end").unwrap_or(&NULL),
    );
    format!("{chrom}:{coords}")
}

/// `VCF.pm:401-404` on a value already joined for the default format.
pub fn vcf_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut in_space = false;
    for c in s.chars() {
        let space = pyfmt::char_is_space(c);
        match c {
            _ if space && in_space => {}
            _ if space => out.push('_'),
            ',' | '|' => out.push('&'),
            ';' => out.push_str("%3B"),
            c => out.push(c),
        }
        in_space = space;
    }
    out
}

/// `vcf_escape(text.replace("%3B", ";"))`.
fn as_csq(text: &str) -> Cow<'_, str> {
    let escaped = |c: char| matches!(c, ',' | ';' | '|') || pyfmt::char_is_space(c);
    if !text.contains("%3B") && !text.chars().any(escaped) {
        return Cow::Borrowed(text);
    }
    Cow::Owned(vcf_escape(&text.replace("%3B", ";")))
}

/// `text.replace("%3B", ";")`.
fn semicolons(text: &str) -> Cow<'_, str> {
    if text.contains("%3B") {
        Cow::Owned(text.replace("%3B", ";"))
    } else {
        Cow::Borrowed(text)
    }
}

/// The reverse complement of a sequence allele; any other allele unchanged.
pub fn reverse_complement(allele: &str) -> String {
    let sequence = !allele.is_empty()
        && allele
            .bytes()
            .all(|base| base == b'-' || complement(base).is_some());
    if !sequence {
        return allele.to_string();
    }
    allele
        .bytes()
        .rev()
        .map(|base| char::from(complement(base).unwrap_or(base)))
        .collect()
}

fn complement(base: u8) -> Option<u8> {
    COMPLEMENT
        .iter()
        .find(|(from, _)| *from == i64::from(base))
        .and_then(|(_, to)| u8::try_from(*to).ok())
}

/// Equality of two scalar texts, numerically when both are numbers.
pub fn values_equal(a: &str, b: &str) -> bool {
    if a == b {
        return true;
    }
    match (pyfmt::float_of(a), pyfmt::float_of(b)) {
        (Some(a), Some(b)) => a == b,
        _ => false,
    }
}

pub fn lists_equal(a: &[String], b: &[String]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| values_equal(x, y))
}

/// `lists_equal` of two lists given item by item.
fn items_equal<'a>(
    mut a: impl Iterator<Item = &'a str>,
    mut b: impl Iterator<Item = &'a str>,
) -> bool {
    loop {
        match (a.next(), b.next()) {
            (None, None) => return true,
            (Some(x), Some(y)) if values_equal(x, y) => {}
            _ => return false,
        }
    }
}

/// `decode_extra(value)`, the text itself where it holds nothing to decode.
fn decoded(value: &str) -> Cow<'_, str> {
    if value.contains('%') {
        Cow::Owned(decode_extra(value))
    } else {
        Cow::Borrowed(value)
    }
}

/// A default-format value as the list `convert_arrayref` joined, percent-encoding
/// undone; `None` for `-`.
pub fn default_list(value: &str) -> Option<Vec<String>> {
    if value == "-" || value == ABSENT {
        return None;
    }
    Some(decoded(value).split(',').map(str::to_string).collect())
}

/// A JSON value as a list of scalar texts, a string split on `,` as the default
/// format's list would be; `None` when absent.
///
/// A JSON `null` is an absent value.
pub fn json_list(value: &Value) -> Option<Vec<String>> {
    json_list_under(value, None)
}

/// `json_list(value)` called by a function that is `callers` frames deep.
fn json_list_under(value: &Value, callers: Option<usize>) -> Option<Vec<String>> {
    let scalar = |value: &Value| scalar_under(value, callers.map(|callers| callers + 1));
    match value {
        Value::Null => None,
        Value::List(items) => Some(items.iter().map(scalar).collect()),
        other => Some(scalar(other).split(',').map(str::to_string).collect()),
    }
}

/// `JSON.pm:387-388`: the lowest-ranked term, `?` when there is none.
pub fn most_severe(terms: &BTreeSet<String>) -> String {
    let rank = |term: &str| {
        CONSEQUENCE_RANK
            .iter()
            .find(|(ranked, _)| *ranked == term)
            .map_or(1_000_000, |(_, rank)| *rank)
    };
    terms
        .iter()
        .filter(|term| *term != "?")
        .min_by_key(|term| (rank(term), term.as_str()))
        .map_or_else(|| "?".to_string(), String::clone)
}

/// A name and a Location as the key of `SelfParity::objects_by_name`, borrowed.
#[derive(Hash)]
struct Named<'a>(&'a str, &'a str);

impl Equivalent<(String, String)> for Named<'_> {
    fn equivalent(&self, key: &(String, String)) -> bool {
        self.0 == key.0 && self.1 == key.1
    }
}

fn text(value: impl Into<String>) -> Value {
    Value::Str(value.into())
}

fn count(items: usize) -> Value {
    Value::Int(i64::try_from(items).unwrap_or(i64::MAX))
}

fn texts<S: AsRef<str>>(items: impl IntoIterator<Item = S>) -> Value {
    Value::List(items.into_iter().map(|item| text(item.as_ref())).collect())
}

fn object<const N: usize>(members: [(&str, Value); N]) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, member)| (name.to_string(), member))
            .collect(),
    )
}

/// `items[:n]` of a list.
fn head<T>(items: &[T], n: i64) -> &[T] {
    let length = i64::try_from(items.len()).unwrap_or(i64::MAX);
    let end = if n < 0 {
        (length + n).max(0)
    } else {
        n.min(length)
    };
    &items[..usize::try_from(end).unwrap_or(0)]
}

/// The items of `text` between commas, sorted and joined again.
fn sorted_items(text: &str) -> String {
    let mut items: Vec<&str> = text.split(',').collect();
    items.sort_unstable();
    items.join(",")
}

/// The keys of `entries` that `kept` takes, sorted.
fn sorted_keys<'a>(
    entries: &IndexMap<&'a str, &Entry>,
    kept: impl Fn(&str) -> bool,
) -> Vec<&'a str> {
    let mut keys: Vec<&str> = entries.keys().copied().filter(|key| kept(key)).collect();
    keys.sort_unstable();
    keys
}

fn noted(compared: &mut BTreeSet<String>, column: &str) {
    if !compared.contains(column) {
        compared.insert(column.to_string());
    }
}

/// Whether a wanted text and the text a JSON entry holds differ: one of them absent, or
/// two texts that are not the same value.
fn scalar_differs(want: Option<&str>, got: Option<&str>) -> bool {
    match (want, got) {
        (Some(want), Some(got)) => !values_equal(want, got),
        (None, None) => false,
        _ => true,
    }
}

/// One engine's outputs on one input, every consequence entry keyed by the record's
/// `Uploaded_variation` and `Location`, its allele, its feature and its feature kind,
/// and compared column by column against the default format.
#[derive(Debug)]
pub struct SelfParity {
    pub examples: i64,
    /// The JSON objects in file order; `records` and `objects_by_name` name an object by
    /// its place here.
    pub objects: Vec<Value>,
    /// The object of each record label.
    pub records: IndexMap<String, usize>,
    /// The objects of each name and Location, in file order.
    pub objects_by_name: IndexMap<(String, String), Vec<usize>>,
    /// For each object, the place of its name and Location in `objects_by_name`.
    named: Vec<usize>,
    /// For each object, the most frames the interpreter has on its call stack above
    /// `uploaded_variation` while it writes the object's name.
    name_frames: Vec<usize>,
    /// The name, Location and feature kind of every key of the default format, as the last
    /// format that holds the key gives them; the parts of no other key are read.
    pub parts: HashMap<String, (String, String, String)>,
    /// Whether the default format's keys are all in `parts`.
    pivot_keyed: bool,
    pub default: Vec<Entry>,
    pub tab: Vec<Entry>,
    pub parquet: Option<Vec<Entry>>,
    pub vcf_records: Vec<Entry>,
    pub vcf: Vec<Entry>,
    pub json: Vec<Entry>,
    pub tab_columns: Vec<String>,
    pub parquet_columns: Vec<String>,
    pub csq_columns: Vec<String>,
    pub default_columns: BTreeSet<String>,
}

impl SelfParity {
    /// `parquet_text` is `None` for a run without a Parquet output; `inputs` is `None` where
    /// the input lines are the JSON objects' own `input` texts; `examples` is the most
    /// examples a list of the report holds.
    pub fn new(
        default_text: &str,
        tab_text: &str,
        vcf_text: &str,
        json_text: &str,
        parquet_text: Option<&str>,
        inputs: Option<&[String]>,
        examples: i64,
    ) -> Result<Self> {
        let objects = json_objects(json_text)?;
        let echoed: Vec<String>;
        let inputs = match inputs {
            Some(inputs) => inputs,
            None => {
                echoed = objects
                    .iter()
                    .filter_map(|(_, obj)| obj.get("input").and_then(Value::as_str))
                    .map(str::to_string)
                    .collect();
                &echoed
            }
        };
        let mut index = InputIndex::new(inputs, &mut |line: &str| vcf_line_key(line));
        let mut parity = SelfParity {
            examples,
            objects: Vec::with_capacity(objects.len()),
            records: IndexMap::new(),
            objects_by_name: IndexMap::new(),
            named: Vec::with_capacity(objects.len()),
            name_frames: Vec::with_capacity(objects.len()),
            parts: HashMap::new(),
            pivot_keyed: false,
            default: Vec::new(),
            tab: Vec::new(),
            parquet: None,
            vcf_records: Vec::new(),
            vcf: Vec::new(),
            json: Vec::new(),
            tab_columns: Vec::new(),
            parquet_columns: Vec::new(),
            csq_columns: Vec::new(),
            default_columns: BTreeSet::new(),
        };
        for (n, obj) in objects {
            let (_, label) = json_record_label(&mut index, n, &obj);
            let at = parity.objects.len();
            parity.records.insert(label, at);
            let name = uploaded_variation_under(&obj, Some(METHOD_FRAMES));
            let location = record_location_under(&obj, Some(METHOD_FRAMES));
            let entry = parity.objects_by_name.entry((name, location));
            parity.named.push(entry.index());
            entry.or_default().push(at);
            let name_frames = ["seq_region_name", "allele_string", "id"]
                .iter()
                .filter_map(|member| obj.get(member))
                .map(scalar_frames)
                .max()
                .unwrap_or(0);
            parity.name_frames.push(name_frames);
            parity.objects.push(obj);
        }
        let default = parse_default_under(default_text, METHOD_FRAMES)?;
        parity.default = parity._rekey_rows(default);
        parity.pivot_keyed = true;
        let tab = parse_tab_under(tab_text, METHOD_FRAMES)?;
        parity.tab = parity._rekey_rows(tab);
        if let Some(parquet_text) = parquet_text {
            let parquet = parse_tab_under(parquet_text, METHOD_FRAMES)?;
            parity.parquet = Some(parity._rekey_rows(parquet));
        }
        let (vcf_records, vcf_entries): (Vec<Entry>, Vec<Entry>) = parse_vcf(vcf_text, inputs)?
            .into_iter()
            .partition(|e| e.key.ends_with("|record"));
        parity.vcf_records = vcf_records;
        parity.vcf = parity._rekey_by_record(vcf_entries);
        let json = parse_json_under(json_text, inputs, METHOD_FRAMES)?;
        parity.json = parity._rekey_by_record(json);
        parity.tab_columns = header_columns(tab_text);
        parity.parquet_columns = parquet_text.map(header_columns).unwrap_or_default();
        parity.csq_columns = Self::_csq_columns(vcf_text)?;
        for e in &parity.default {
            for column in e.fields.keys() {
                noted(&mut parity.default_columns, column);
            }
        }
        Ok(parity)
    }

    pub fn _csq_columns(text: &str) -> Result<Vec<String>> {
        match line_slices(text).find(|line| line.starts_with("##INFO=<ID=CSQ")) {
            Some(line) => csq_format_fields(line),
            None => Ok(Vec::new()),
        }
    }

    /// Entries keyed by (name, Location, allele, feature, kind), duplicates numbered,
    /// each key's parts kept in `self.parts` (a name may itself carry `|`): those of every
    /// key until the default format is keyed, and of the keys it holds from then on.
    pub fn _rekeyed(&mut self, keyed: Vec<(Entry, String, String, String)>) -> Vec<Entry> {
        let mut out = Vec::with_capacity(keyed.len());
        let mut parts = Vec::with_capacity(keyed.len());
        for (e, name, location, kind) in keyed {
            out.push(Entry {
                key: format!("{name}|{location}|{}|{}|{kind}", e.allele, dash(&e.feature)),
                record: e.record,
                location: e.location,
                allele: e.allele,
                feature: e.feature,
                consequence: e.consequence,
                fields: e.fields,
                ..Entry::default()
            });
            parts.push((name, location, kind));
        }
        dedup_keys(&mut out);
        for (entry, parts) in out.iter().zip(parts) {
            match self.parts.get_mut(&entry.key) {
                Some(held) => *held = parts,
                None if !self.pivot_keyed => {
                    self.parts.insert(entry.key.clone(), parts);
                }
                None => {}
            }
        }
        out
    }

    /// Default-format rows keyed by their own `Uploaded_variation` and `Location`.
    pub fn _rekey_rows(&mut self, rows: Vec<Entry>) -> Vec<Entry> {
        let keyed = rows
            .into_iter()
            .map(|row| {
                let field = |name: &str| row.fields.get(name).map_or("", String::as_str);
                let name = field("Uploaded_variation").to_string();
                let location = field("Location").to_string();
                let kind = feature_kind(field("Feature_type"));
                (row, name, location, kind)
            })
            .collect();
        self._rekeyed(keyed)
    }

    /// VCF and JSON consequence entries keyed by the name and Location of the JSON
    /// object of their record; an entry whose record has no object keeps its label.
    pub fn _rekey_by_record(&mut self, entries: Vec<Entry>) -> Vec<Entry> {
        let mut keyed = Vec::with_capacity(entries.len());
        for e in entries {
            if e.key.ends_with("|record") {
                continue;
            }
            let parts: Vec<&str> = e.key.splitn(4, '|').collect();
            let &[label, _, _, kind] = parts.as_slice() else {
                crate::uncaught!(
                    "not enough values to unpack (expected 4, got {})",
                    parts.len()
                );
            };
            let kind = feature_kind(kind.split('#').next().unwrap_or(""));
            let (name, location) = match self.records.get(label) {
                None => (label.to_string(), String::new()),
                Some(&at) => {
                    if METHOD_FRAMES + 2 + self.name_frames[at] > pyjson::RECURSION_LIMIT {
                        crate::uncaught!("{RECURSION}");
                    }
                    self.name_of(at)
                }
            };
            keyed.push((e, name, location, kind));
        }
        self._rekeyed(keyed)
    }

    /// The name and Location of the object at `at`.
    fn name_of(&self, at: usize) -> (String, String) {
        self.objects_by_name
            .get_index(self.named[at])
            .map(|(named, _)| named.clone())
            .unwrap_or_default()
    }

    pub fn _by_key(entries: &[Entry]) -> IndexMap<&str, &Entry> {
        entries.iter().map(|e| (e.key.as_str(), e)).collect()
    }

    pub fn _pair(
        &self,
        name: &str,
        pivot: &IndexMap<&str, &Entry>,
        other: &IndexMap<&str, &Entry>,
        compare_fields: &CompareFields<'_>,
    ) -> Value {
        let pivot_only = sorted_keys(pivot, |key| !other.contains_key(key));
        let other_only = sorted_keys(other, |key| !pivot.contains_key(key));
        let shared = sorted_keys(pivot, |key| other.contains_key(key));
        let mut differing: BTreeMap<String, (i64, Vec<String>)> = BTreeMap::new();
        let mut compared_columns: BTreeSet<String> = BTreeSet::new();
        for key in &shared {
            for (column, detail) in compare_fields(pivot[key], other[key], &mut compared_columns) {
                let (total, examples) = differing.entry(column).or_default();
                *total += 1;
                if i64::try_from(examples.len()).unwrap_or(i64::MAX) < self.examples {
                    examples.push(format!("{key}: {detail}"));
                }
            }
        }
        let mut reasons = Vec::new();
        if !pivot_only.is_empty() {
            reasons.push(format!(
                "{name}: {} entries of the default format are absent",
                pivot_only.len()
            ));
        }
        if !other_only.is_empty() {
            reasons.push(format!(
                "{name}: {} entries have no default-format row",
                other_only.len()
            ));
        }
        for (column, (total, _)) in &differing {
            reasons.push(format!(
                "{name}: column {column} differs on {total} entries"
            ));
        }
        let other_name = format!("{name}_only");
        let other_examples = format!("{name}_only_examples");
        let differing = differing
            .into_iter()
            .map(|(column, (total, examples))| {
                let difference =
                    object([("count", Value::Int(total)), ("examples", texts(examples))]);
                (column, difference)
            })
            .collect();
        object([
            ("shared_entries", count(shared.len())),
            (
                "one_sided",
                object([
                    ("default_only", count(pivot_only.len())),
                    (other_name.as_str(), count(other_only.len())),
                ]),
            ),
            (
                "default_only_examples",
                texts(head(&pivot_only, self.examples)),
            ),
            (
                other_examples.as_str(),
                texts(head(&other_only, self.examples)),
            ),
            ("columns_compared", texts(&compared_columns)),
            ("differing", Value::Dict(differing)),
            ("reasons", texts(reasons)),
        ])
    }

    pub fn _compare_tab_like(&self, columns: &[String]) -> Box<CompareFields<'_>> {
        let declared: Vec<String> = columns
            .iter()
            .filter(|column| *column != "Extra")
            .cloned()
            .collect();
        Box::new(move |p, o, compared| {
            let mut differences = Vec::new();
            for column in &declared {
                noted(compared, column);
                let raw = p.fields.get(column).map_or(ABSENT, String::as_str);
                let pv = semicolons(raw);
                let ov = o.fields.get(column).map_or(ABSENT, String::as_str);
                if pv == ABSENT && ov == "-" {
                    continue;
                }
                if pv != ov {
                    let detail = format!(
                        "default {}, other {}",
                        pyfmt::repr_str(raw),
                        pyfmt::repr_str(ov)
                    );
                    differences.push((column.clone(), detail));
                }
            }
            differences
        })
    }

    pub fn _compare_vcf(&self) -> Box<CompareFields<'_>> {
        let columns: Vec<&String> = self
            .csq_columns
            .iter()
            .filter(|column| {
                !VCF_IMPLIED_COLUMNS.contains(&column.as_str())
                    || self.default_columns.contains(*column)
            })
            .collect();
        Box::new(move |p, o, compared| {
            let reverse = self
                ._object_for(p)
                .get("strand")
                .is_some_and(|strand| pyjson::py_eq(strand, &Value::Int(-1)));
            let mut differences = Vec::new();
            for column in &columns {
                noted(compared, column);
                let pv = p.fields.get(*column).map_or(ABSENT, String::as_str);
                let ov = o.fields.get(*column).map_or("", String::as_str);
                let turned;
                let expected = if *column == "Allele" {
                    let allele = if pv == ABSENT { "-" } else { pv };
                    if reverse {
                        turned = reverse_complement(allele);
                        as_csq(&turned)
                    } else {
                        as_csq(allele)
                    }
                } else if pv == ABSENT || pv == "-" {
                    Cow::Borrowed("")
                } else {
                    as_csq(pv)
                };
                if expected != ov {
                    let detail = format!(
                        "default {} (as CSQ {}), vcf {}",
                        pyfmt::repr_str(pv),
                        pyfmt::repr_str(&expected),
                        pyfmt::repr_str(ov)
                    );
                    differences.push(((*column).clone(), detail));
                }
            }
            differences
        })
    }

    pub fn _compare_json(&self) -> Box<CompareFields<'_>> {
        Box::new(move |p, o, compared| self.json_differences(p, o, compared))
    }

    /// The differences of a default-format entry `p` and its JSON entry `o`.
    fn json_differences(
        &self,
        p: &Entry,
        o: &Entry,
        compared: &mut BTreeSet<String>,
    ) -> Vec<(String, String)> {
        let kind = &self.parts_of(p).2;
        let mut remaining = vec![true; o.fields.len()];
        let mut discard = |jkey: &str| {
            if let Some(at) = o.fields.get_index_of(jkey) {
                remaining[at] = false;
            }
        };
        let got_of = |jkey: &str| o.fields.get(jkey).map(String::as_str);
        let mut differences = Vec::new();
        for (column, raw) in &p.fields {
            if column == "Extra"
                || column == "Feature_type"
                || RECORD_LEVEL_COLUMNS.contains(&column.as_str())
                || RECORD_LEVEL_COLUMN_RE.is_match(column)
            {
                continue;
            }
            noted(compared, column);
            let mut differs = |jkey: &str, got: Option<&str>| {
                let detail = format!(
                    "default {}, json {jkey}={}",
                    pyfmt::repr_str(raw),
                    pyfmt::repr_opt_str(got)
                );
                differences.push((column.clone(), detail));
            };
            let key = pyfmt::lower(column);
            let filled = raw != "-" && !raw.is_empty();
            // `re.match` holds the pattern to the start of the name.
            let coordinate = POSITION_RE
                .captures(column)
                .filter(|found| found.get(0).is_some_and(|whole| whole.start() == 0))
                .and_then(|found| found.get(1));
            if let Some(coordinate) = coordinate {
                let coordinate = pyfmt::lower(coordinate.as_str());
                let start_key = format!("{coordinate}_start");
                let end_key = format!("{coordinate}_end");
                discard(&start_key);
                discard(&end_key);
                let (mut want_start, mut want_end) = (None, None);
                if filled {
                    let mut parts = raw.split('-');
                    let start_text = parts.next().unwrap_or("");
                    let end_text = parts.next().filter(|part| !part.is_empty());
                    want_start = POSITION_PART_RE.is_match(start_text).then_some(start_text);
                    let end = match end_text {
                        Some(end) if pyfmt::is_digit(end) || end.contains('?') => end,
                        _ => start_text,
                    };
                    want_end =
                        (end_text.is_none() || POSITION_PART_RE.is_match(end)).then_some(end);
                }
                for (jkey, want) in [(&start_key, want_start), (&end_key, want_end)] {
                    let got = got_of(jkey);
                    if scalar_differs(want, got) {
                        differs(jkey, got);
                    }
                }
                continue;
            }
            if key == "sift" || key == "polyphen" {
                let pred_key = format!("{key}_prediction");
                let score_key = format!("{key}_score");
                discard(&pred_key);
                discard(&score_key);
                let found = filled.then(|| TOOL_SCORE_RE.captures(raw)).flatten();
                let group = |n: usize| {
                    found
                        .as_ref()
                        .and_then(|found| found.get(n))
                        .map(|group| group.as_str())
                };
                let wanted = [
                    (
                        &pred_key,
                        group(1).filter(|prediction| !prediction.is_empty()),
                    ),
                    (&score_key, group(2)),
                ];
                for (jkey, want) in wanted {
                    let got = got_of(jkey);
                    if scalar_differs(want, got) {
                        differs(jkey, got);
                    }
                }
                continue;
            }
            if key == "domains" {
                discard("domains");
                let want = filled.then(|| domains_model(raw));
                let got = got_of("domains");
                if want.as_deref() != got {
                    differs("domains", got);
                }
                continue;
            }
            if JSON_OBJECT_FIELDS.contains(&column.as_str()) {
                discard(&key);
                let want = filled.then(|| numberified_json_text(raw));
                let got = got_of(&key);
                if want.as_deref() != got {
                    differs(&key, got);
                }
                continue;
            }
            let kind_id;
            let jkey = if key == "feature" {
                kind_id = format!("{kind}_id");
                kind_id.as_str()
            } else {
                JSON_RENAME
                    .iter()
                    .find(|(name, _)| *name == key)
                    .map_or(key.as_str(), |(_, renamed)| renamed)
            };
            discard(jkey);
            let listed;
            let want = if raw == "-" || raw == ABSENT {
                matches!(key.as_str(), "allele" | "used_ref" | "given_ref").then_some("-")
            } else {
                listed = decoded(raw);
                Some(listed.as_ref())
            };
            let got = got_of(jkey);
            let equal = match (want, got) {
                (None, None) => true,
                (Some(want), Some(got)) => items_equal(
                    want.split(',')
                        .map(|item| if item == "YES" { "1" } else { item }),
                    got.split(','),
                ),
                _ => false,
            };
            if !equal {
                differs(jkey, got);
            }
        }
        let mut left: Vec<(&String, &String)> = o
            .fields
            .iter()
            .zip(&remaining)
            .filter(|(_, remains)| **remains)
            .map(|(field, _)| field)
            .collect();
        left.sort_unstable();
        for (jkey, value) in left {
            noted(compared, jkey);
            let detail = format!(
                "json carries {jkey}={} and the default format no column for it",
                pyfmt::repr_str(value)
            );
            differences.push((jkey.clone(), detail));
        }
        differences
    }

    /// The name, Location and feature kind of an entry's key.
    fn parts_of(&self, e: &Entry) -> &(String, String, String) {
        match self.parts.get(&e.key) {
            Some(parts) => parts,
            None => crate::uncaught!("{}", pyfmt::repr_str(&e.key)),
        }
    }

    /// The JSON object whose name and Location an entry's key carries, else empty.
    pub fn _object_for(&self, e: &Entry) -> &Value {
        let (name, location, _) = self.parts_of(e);
        self.objects_by_name
            .get(&Named(name, location))
            .and_then(|objects| objects.first())
            .map_or(&NO_OBJECT, |at| &self.objects[*at])
    }

    /// Record-level values the JSON format carries once: `most_severe_consequence`
    /// against the rank order of the object's own terms (`JSON.pm:387-388`),
    /// `colocated_variants` ids against the rows' `Existing_variation`, the `zyg`
    /// value against the rows' `ZYG` (the list of `--individual_zyg` as a multiset, since
    /// its order carries no meaning; the one zygosity of `--individual` as itself), and the
    /// VCF `ID` column against the record's name.
    pub fn _record_checks(&self) -> Value {
        let mut problems: Vec<String> = Vec::new();
        let mut existing_by_name: HashMap<(&str, &str), BTreeSet<&str>> = HashMap::new();
        let mut zyg_by_name: HashMap<(&str, &str), BTreeSet<String>> = HashMap::new();
        for row in &self.default {
            let (name, location, _) = self.parts_of(row);
            let named = (name.as_str(), location.as_str());
            let field = |name: &str| row.fields.get(name).map_or("-", String::as_str);
            existing_by_name
                .entry(named)
                .or_default()
                .insert(field("Existing_variation"));
            let zyg = field("ZYG");
            if zyg != "-" {
                zyg_by_name
                    .entry(named)
                    .or_default()
                    .insert(sorted_items(zyg));
            }
        }
        let mut by_name: Vec<(&(String, String), &Vec<usize>)> =
            self.objects_by_name.iter().collect();
        by_name.sort_unstable_by_key(|(named, _)| *named);
        for ((name, location), objects) in by_name {
            let named = (name.as_str(), location.as_str());
            let existing = existing_by_name.get(&named);
            let want_ids: BTreeSet<String> = existing
                .into_iter()
                .flatten()
                .map(|value| match default_list(value) {
                    Some(mut ids) => {
                        ids.sort_unstable();
                        ids.join(",")
                    }
                    None => String::new(),
                })
                .collect();
            let mut got_ids: BTreeSet<String> = BTreeSet::new();
            let mut got_zyg: BTreeSet<String> = BTreeSet::new();
            for obj in objects.iter().map(|at| &self.objects[*at]) {
                let mut terms: BTreeSet<String> = BTreeSet::new();
                for (list_name, _) in JSON_CONSEQUENCE_LISTS {
                    for item in iterated(obj.get(list_name)) {
                        if item.as_dict().is_none() {
                            continue;
                        }
                        let listed = json_list_under(
                            item.get("consequence_terms").unwrap_or(&NULL),
                            Some(METHOD_FRAMES + 1),
                        );
                        terms.extend(listed.unwrap_or_default());
                    }
                }
                // An object with no consequence has no row to hold its zyg to: the
                // reference-homozygous copy of --individual (OutputFactory.pm:511 takes the
                // alternate alleles alone without --process_ref_homs) and a record carried
                // with no allele to annotate.
                if let Some(zyg) = obj.get("zyg").filter(|_| !terms.is_empty()) {
                    got_zyg.insert(sorted_items(&json_scalar(zyg)));
                }
                let want = most_severe(&terms);
                let got = obj
                    .get("most_severe_consequence")
                    .filter(|got| !matches!(got, Value::Null));
                if got.is_none() && terms.is_empty() && obj.get("id").is_none() {
                    continue;
                }
                if got.and_then(Value::as_str) != Some(want.as_str()) {
                    problems.push(format!(
                        "{name}|{location}: most_severe_consequence {}, the object's terms give {}",
                        got.map_or_else(|| "None".to_string(), pyjson::py_repr),
                        pyfmt::repr_str(&want)
                    ));
                }
                let mut ids: Vec<String> = obj
                    .get("colocated_variants")
                    .and_then(Value::as_list)
                    .unwrap_or_default()
                    .iter()
                    .filter(|colocated| colocated.as_dict().is_some())
                    .map(|colocated| colocated.get("id").map_or_else(String::new, json_scalar))
                    .collect();
                ids.sort_unstable();
                got_ids.insert(ids.join(","));
            }
            if existing.is_none() {
                continue;
            }
            if want_ids != got_ids {
                problems.push(format!(
                    "{name}|{location}: Existing_variation {} against colocated_variants ids {}",
                    repr_text_list(&want_ids),
                    repr_text_list(&got_ids)
                ));
            }
            let no_zyg = BTreeSet::new();
            let want_zyg = zyg_by_name.get(&named).unwrap_or(&no_zyg);
            if *want_zyg != got_zyg {
                problems.push(format!(
                    "{name}|{location}: ZYG {} against the objects' zyg {}",
                    repr_text_list(want_zyg),
                    repr_text_list(&got_zyg)
                ));
            }
        }
        for e in &self.vcf_records {
            let label = e.key.split('|').next().unwrap_or("");
            let vcf_id = e.fields.get("ID").map_or(".", String::as_str);
            let Some(&at) = self.records.get(label) else {
                continue;
            };
            let obj = &self.objects[at];
            if vcf_id == "." || obj.get("id").is_none() {
                continue;
            }
            let name = uploaded_variation_under(obj, Some(METHOD_FRAMES + 1));
            if name != vcf_id {
                problems.push(format!(
                    "{label}: VCF ID {} against the name {}",
                    pyfmt::repr_str(vcf_id),
                    pyfmt::repr_str(&name)
                ));
            }
        }
        object([("problems", texts(problems))])
    }

    pub fn run(&self) -> Value {
        let pivot = Self::_by_key(&self.default);
        let mut pairs: IndexMap<String, Value> = IndexMap::new();
        let mut pair = |name: &str, entries: &[Entry], compare: Box<CompareFields<'_>>| {
            let report = self._pair(name, &pivot, &Self::_by_key(entries), &compare);
            pairs.insert(name.to_string(), report);
        };
        pair("tab", &self.tab, self._compare_tab_like(&self.tab_columns));
        pair("vcf", &self.vcf, self._compare_vcf());
        pair("json", &self.json, self._compare_json());
        if let Some(parquet) = &self.parquet {
            pair(
                "parquet",
                parquet,
                self._compare_tab_like(&self.parquet_columns),
            );
        }
        let records = self._record_checks();
        let listed = |value: Option<&Value>| -> Vec<Value> {
            value
                .and_then(Value::as_list)
                .map(<[Value]>::to_vec)
                .unwrap_or_default()
        };
        let mut reasons: Vec<Value> = pairs
            .values()
            .flat_map(|pair| listed(pair.get("reasons")))
            .collect();
        reasons.extend(listed(records.get("problems")));
        let record_level = |column: &str| RECORD_LEVEL_COLUMNS.contains(&column);
        let keyed_or_record_checked = [
            "Uploaded_variation",
            "Location",
            "Existing_variation",
            "ZYG",
        ];
        let default_columns = || self.default_columns.iter().map(String::as_str);
        let mut vcf_implied: Vec<&str> = self
            .csq_columns
            .iter()
            .map(String::as_str)
            .filter(|column| {
                VCF_IMPLIED_COLUMNS.contains(column) && !self.default_columns.contains(*column)
            })
            .collect();
        vcf_implied.sort_unstable();
        let declared = |columns: &[String], column: &str| columns.iter().any(|c| c == column);
        let not_compared =
            object([
                (
                    "default_record_level",
                    texts(default_columns().filter(|column| {
                        (record_level(column) && !keyed_or_record_checked.contains(column))
                            || RECORD_LEVEL_COLUMN_RE.is_match(column)
                    })),
                ),
                ("vcf_implied", texts(vcf_implied)),
                (
                    "tab_undeclared",
                    texts(default_columns().filter(|column| {
                        *column != "Extra" && !declared(&self.tab_columns, column)
                    })),
                ),
                (
                    "vcf_undeclared",
                    texts(default_columns().filter(|column| {
                        *column != "Extra"
                            && !declared(&self.csq_columns, column)
                            && !record_level(column)
                    })),
                ),
            ]);
        let verdict = if reasons.is_empty() { "green" } else { "red" };
        object([
            ("mode", text("self-parity")),
            ("verdict", text(verdict)),
            ("reasons", Value::List(reasons)),
            ("records", count(self.records.len())),
            ("default_entries", count(pivot.len())),
            ("pairs", Value::Dict(pairs)),
            ("record_checks", records),
            ("not_compared", not_compared),
        ])
    }
}

/// The objects the default format's `DOMAINS` text stands for, as the JSON format's
/// `domains` entry joins them.
///
/// OutputFactory/JSON.pm:363-371 re-splits each entry with `(\w+):(\w+)` and pushes
/// `{db => $1, name => $2}` if `$1 && $2`. A failed match leaves `$1` and `$2` at the
/// list's last successful one, so such an entry repeats the previous object and is dropped
/// only before the first match. A capture of `0` is false to Perl, so its entry and every
/// unmatched one after it push nothing. `numberify` then runs over the record
/// (JSON.pm:288), so a name of digits alone is a number. `json_scalar` joins the list's
/// compact objects with `,`, which each object also contains, so the comparison is on the
/// joined text.
fn domains_model(raw: &str) -> String {
    let listed = decoded(raw);
    let mut want: Vec<String> = Vec::new();
    let mut last: Option<(&str, &str)> = None;
    for item in listed.split(',') {
        if let Some(found) = DOMAIN_RE.captures(item) {
            last = found
                .get(1)
                .zip(found.get(2))
                .map(|(db, name)| (db.as_str(), name.as_str()));
        }
        if let Some((db, name)) = last.filter(|(db, name)| *db != "0" && *name != "0") {
            let domain = object([
                ("db", perl_numeric_or_text(db).into_value()),
                ("name", perl_numeric_or_text(name).into_value()),
            ]);
            want.push(json_compact(&domain));
        }
    }
    want.join(",")
}

/// The items the interpreter's `for item in value or []` walks that can be an object: those
/// of a list. A text and an object are walked and hold no object; a number or `true` ends
/// the run.
fn iterated(value: Option<&Value>) -> &[Value] {
    match value {
        Some(Value::List(items)) => items,
        Some(value) if pyjson::truthy(value) => match value {
            Value::Bool(_) => crate::uncaught!("'bool' object is not iterable"),
            Value::Float(_) => crate::uncaught!("'float' object is not iterable"),
            Value::Int(_) | Value::BigInt(_) => crate::uncaught!("'int' object is not iterable"),
            _ => &[],
        },
        _ => &[],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn every_char() -> impl Iterator<Item = char> {
        (0..=0x10FFFF_u32).filter_map(char::from_u32)
    }

    fn word(c: char) -> bool {
        c == '_'
            || pyfmt::char_is_alpha(c)
            || pyfmt::char_is_decimal(c)
            || pyfmt::char_is_digit(c)
            || pyfmt::char_is_numeric(c)
    }

    #[test]
    fn the_domain_pattern_takes_the_word_characters_of_the_interpreter() {
        for c in every_char() {
            let found = DOMAIN_RE
                .captures(&format!("{c}:a"))
                .is_some_and(|found| found[1] == *c.to_string());
            assert_eq!(found, word(c), "U+{:04X} before the colon", c as u32);
            let found = DOMAIN_RE
                .captures(&format!("a:{c}"))
                .is_some_and(|found| found[2] == *c.to_string());
            assert_eq!(found, word(c), "U+{:04X} after the colon", c as u32);
        }
    }

    #[test]
    fn the_number_pattern_skips_the_white_space_of_the_interpreter() {
        for c in every_char() {
            if c.is_ascii_digit() || matches!(c, '+' | '-' | '.') || pyfmt::char_is_decimal(c) {
                continue;
            }
            let found = PERL_NUMBER_RE
                .captures(&format!("{c}7"))
                .is_some_and(|found| &found[1] == "7");
            assert_eq!(found, pyfmt::char_is_space(c), "U+{:04X}", c as u32);
        }
    }

    #[test]
    fn the_patterns_take_the_decimal_digits_of_the_interpreter() {
        for c in every_char() {
            let decimal = pyfmt::char_is_decimal(c);
            let digit = c.to_string();
            assert_eq!(
                POSITION_PART_RE.is_match(&digit),
                decimal,
                "U+{:04X}",
                c as u32
            );
            let number = PERL_NUMBER_RE
                .captures(&digit)
                .is_some_and(|found| found[1] == digit);
            assert_eq!(number, decimal, "U+{:04X} as a number", c as u32);
            let score = TOOL_SCORE_RE
                .captures(&digit)
                .and_then(|found| found.get(2))
                .is_some_and(|score| score.as_str() == digit);
            assert_eq!(score, decimal || c == '.', "U+{:04X} as a score", c as u32);
        }
    }

    #[test]
    fn the_tool_pattern_takes_the_letters_the_interpreter_folds_into_its_range() {
        for c in every_char() {
            let letter = c.is_ascii_alphabetic()
                || matches!(c, '_' | '\u{130}' | '\u{131}' | '\u{17f}' | '\u{212a}');
            let found = TOOL_SCORE_RE
                .captures(&c.to_string())
                .and_then(|found| found.get(1))
                .is_some();
            assert_eq!(found, letter, "U+{:04X}", c as u32);
        }
    }

    #[test]
    fn a_slice_of_examples_ends_as_the_interpreter_cuts_it() {
        let items = ["a", "b", "c"];
        assert_eq!(head(&items, 2), ["a", "b"]);
        assert_eq!(head(&items, 9), items);
        assert_eq!(head(&items, 0), [""; 0]);
        assert_eq!(head(&items, -1), ["a", "b"]);
        assert_eq!(head(&items, -5), [""; 0]);
    }

    #[test]
    fn a_reverse_complement_turns_sequence_alleles_alone() {
        assert_eq!(reverse_complement("AAC"), "GTT");
        assert_eq!(reverse_complement("acgtn-N"), "N-nacgt");
        assert_eq!(reverse_complement("-"), "-");
        assert_eq!(reverse_complement(""), "");
        assert_eq!(reverse_complement("<DEL>"), "<DEL>");
        assert_eq!(reverse_complement("ACGU"), "ACGU");
    }

    #[test]
    fn an_escaped_value_is_the_text_of_a_csq_field() {
        assert_eq!(vcf_escape("a,b;c d|e"), "a&b%3Bc_d&e");
        assert_eq!(vcf_escape("a \t\u{1c}\u{a0} b"), "a_b");
        assert_eq!(vcf_escape("a\u{200b}b"), "a\u{200b}b");
        assert_eq!(as_csq("a%3Bb c"), "a%3Bb_c");
        assert_eq!(as_csq("plain"), "plain");
    }

    #[test]
    fn two_texts_are_equal_as_text_or_as_numbers() {
        assert!(values_equal("100.00", "100"));
        assert!(!values_equal("1/1", "1"));
        assert!(values_equal("nan", "nan"));
        assert!(!values_equal("NaN", "nan"));
        assert!(values_equal("1_0", "10"));
        assert!(values_equal(" 1", "1.0"));
        assert!(values_equal("inf", "Infinity"));
        assert!(!values_equal("0x10", "16"));
        assert!(!values_equal("\u{1c}1", "1"));
        let list = |items: &[&str]| {
            items
                .iter()
                .map(|item| item.to_string())
                .collect::<Vec<_>>()
        };
        assert!(lists_equal(&list(&["1", "a"]), &list(&["1.0", "a"])));
        assert!(!lists_equal(&list(&["1"]), &list(&["1", "1"])));
        assert!(lists_equal(&[], &[]));
    }

    #[test]
    fn a_word_of_ascii_digits_is_a_number() {
        assert_eq!(perl_numeric_or_text("007"), IntOrStr::Int(7));
        assert_eq!(
            perl_numeric_or_text("PF1"),
            IntOrStr::Str("PF1".to_string())
        );
        assert_eq!(
            perl_numeric_or_text("\u{b2}"),
            IntOrStr::Str("\u{b2}".to_string())
        );
        assert_eq!(perl_numeric_or_text(""), IntOrStr::Str(String::new()));
        assert_eq!(
            perl_numeric_or_text("0099999999999999999999"),
            IntOrStr::BigInt("99999999999999999999".to_string())
        );
    }

    #[test]
    fn the_most_severe_term_is_the_lowest_ranked() {
        let terms = |items: &[&str]| items.iter().map(|item| item.to_string()).collect();
        assert_eq!(
            most_severe(&terms(&["intron_variant", "splice_region_variant", "?"])),
            "splice_region_variant"
        );
        assert_eq!(most_severe(&terms(&["?"])), "?");
        assert_eq!(most_severe(&terms(&[])), "?");
        assert_eq!(most_severe(&terms(&["zz", "aa"])), "aa");
        assert_eq!(
            most_severe(&terms(&["aa", "sequence_variant"])),
            "sequence_variant"
        );
    }

    #[test]
    fn a_default_list_is_the_decoded_items() {
        assert_eq!(default_list("-"), None);
        assert_eq!(default_list(ABSENT), None);
        assert_eq!(default_list(""), Some(vec![String::new()]));
        assert_eq!(
            default_list("a%3Db,c%3Bd"),
            Some(vec!["a=b".to_string(), "c;d".to_string()])
        );
    }
}
